
                 module diag_clouds_W_mod

use time_manager_mod,       only: time_type
use diag_cloud_mod,         only: diag_cloud_avg, diag_cloud_driver
use utilities_mod,          only: open_file, file_exist,   &
                                  check_nml_error, error_mesg,   &
                                  print_version_number, FATAL, NOTE, &
				  WARNING, get_my_pe, close_file
use rad_step_setup_mod,     only: temp, press, pflux, jabs, iabs,  &
                                  ISRAD, IERAD, JSRAD, JERAD, & 
                                  KSRAD, KERAD, land
use rad_utilities_mod,      only: Environment, environment_type
use astronomy_package_mod,  only: get_astronomy_for_clouds,  &
			          get_astronomy_for_clouds_init

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!	           diag cloud radiative properties module
!            currently a wrapper until SKYHI goes away and this
!            module can be consolidated with diag_cloud_mod
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: diag_clouds_W.F90,v 1.2 2001/08/30 15:14:22 fms Exp $'
   character(len=128)  :: tag     =  '$Name: galway $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          diag_clouds_init, diag_clouds_calc

!---------------------------------------------------------------------
!-------- namelist  ---------

integer   :: dummy = 0



namelist /diag_clouds_W_nml /     &
			       dummy


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

integer   :: nsolwg

real, dimension(:), allocatable :: latsv

!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





subroutine diag_clouds_init ( ix, iy, kx, th, ierr)

integer,       intent(in)       :: ix, iy, kx
integer,       intent(out)      :: ierr
real, dimension(:), intent(in)  :: th

      integer            :: unit, ierr, io, ierr2
      integer            :: i, j

!---------------------------------------------------------------------
!-----  read namelist  ------
  
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=diag_clouds_W_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'diag_clouds_W_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'diag_clouds_W', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=diag_clouds_W_nml)
      endif
      call close_file (unit)


!--------------------------------------------------------------------
!  retrieve module variables that come from other modules
!--------------------------------------------------------------------

      call get_astronomy_for_clouds_init (nsolwg)

!--------------------------------------------------------------------
!  define and save latitude for later use
!--------------------------------------------------------------------
      allocate (latsv(size(th,1) ) )
      latsv(:) = th(:)


end subroutine diag_clouds_init



!#################################################################

subroutine diag_clouds_calc ( camtsw, cmxolw, crndlw, ncldsw, &
			       nmxolw, nrndlw, cirabsw, &
			       cvisrfsw, cirrfsw, emmxolw, emrndlw, &
			       is, js, Time_next)

integer,                      intent(in)    :: is, js
type(time_type),              intent(in)    :: Time_next
integer, dimension(:,:),      intent(inout) :: ncldsw, nmxolw, nrndlw
real,    dimension(:,:,:),    intent(inout) :: camtsw, cmxolw, crndlw
real,    dimension(:,:,:,:),  intent(inout) :: cirabsw, cvisrfsw, &
					       cirrfsw, emmxolw, emrndlw

!--------------------------------------------------------------------
!   these arrays define the basic radiative properties of the clouds.
!
!     cmxolw  =  amounts of maximally overlapped longwave clouds in
!                layers from KSRAD to KERAD.
!     crndlw  =  amounts of randomly overlapped longwave clouds in
!                layers from KSRAD to KERAD.
!     nmxolw  =  number of maximally overlapped longwave clouds
!                in each grid column.
!     nrndlw  =  number of maximally overlapped longwave clouds
!                in each grid column.
!     emmxolw =  longwave cloud emissivity for maximally overlapped
!                clouds through layers from KSRAD to KERAD.
!     emrndlw =  longwave cloud emissivity for randomly overlapped
!                clouds through layers from KSRAD to KERAD.
!     camtsw  =  shortwave cloud amounts. the sum of the maximally
!                overlapped and randomly overlapped longwave
!                cloud amounts.
!     ncldsw  =  number of shortwave clouds in each grid column.
!     cirabsw =  absorptivity of clouds in the infrared frequency band.
!                may be zenith angle dependent.
!     cirrfsw =  reflectivity of clouds in the infrared frequency band.
!                may be zenith angle dependent.
!    cvisrfsw =  reflectivity of clouds in the visible frequency band.
!                may be zenith angle dependent.
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!   local variables
!-------------------------------------------------------------------
      real, dimension(:,:,:), allocatable    :: cosangsolar
      integer, dimension(:,:,:), allocatable :: kbtm, ktop
      real, dimension(:,:,:), allocatable    :: temp, qmix, rhum, &
						omega, lgscldelq, &
						cnvcntq
      real, dimension(:,:), allocatable      :: convprc,lat
      real, dimension(:,:,:), allocatable    :: ql, qi, cf , &
                                                cldamt, cuvab, cirab, &
						cuvrf, cirrf, emcld
      integer                                :: j, i, k
      integer                                :: ierr, kc

!--------------------------------------------------------------------
      if (Environment%running_fms) then

!---------------------------------------------------------------------
!     obtain the appropriate zenith angles that are to be used here.
!---------------------------------------------------------------------
        allocate (cosangsolar (ISRAD:IERAD, JSRAD:JERAD, nsolwg) )
        call get_astronomy_for_clouds (cosangsolar   )

!---------------------------------------------------------------------
!     obtain the cloud amounts and areas.
!---------------------------------------------------------------------
        allocate (temp (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
        allocate (qmix (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
        allocate (rhum (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
        allocate (omega (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
        allocate (lgscldelq (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
        allocate (cnvcntq   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
        allocate (convprc   (ISRAD:IERAD, JSRAD:JERAD) )
        allocate (lat       (ISRAD:IERAD, JSRAD:JERAD) )

!       call diag_cloud_avg (iabs(ISRAD), jabs(JSRAD), temp, qmix, &
        call diag_cloud_avg (is, js, temp, qmix, &
			     rhum, omega, lgscldelq, cnvcntq, convprc, &
			     ierr)

!---------------------------------------------------------------------
!     allocate and initialize the cloud radiative property arrays.
!---------------------------------------------------------------------
        if (ierr == 0) then
          allocate (ktop   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (kbtm   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (cldamt (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (cuvrf  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (cirrf  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (cuvab  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (cirab  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (emcld  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )

	  emcld = 1.
	  ktop=1
	  kbtm=0
	  cldamt=0.
	  cuvrf=0.
	  cirrf=0.
	  cuvab=0.
	  cirab=0.

          do j=jsrad,jerad
	    lat(:,j) = latsv(jabs(j))
	  end do

!---------------------------------------------------------------------
!     obtain the cloud radiative properties.
!---------------------------------------------------------------------
          call diag_cloud_driver (is, js, &
		temp, qmix, rhum, omega, lgscldelq, cnvcntq, convprc, &
		0.1*press(:,:,KSRAD:KERAD),&
	  	0.1*pflux, 0.1*pflux(:,:,KERAD+1),  &
			      cosangsolar(:,:,1), lat, Time_next, &
			      ncldsw, ktop, kbtm, cldamt, cuvrf,  &
                              cirrf, cuvab, cirab, emcld)

!---------------------------------------------------------------------
!    map the cloud-space arrays to physical space arrays
!-------------------------------------------------------------------
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
              do kc=1, ncldsw(i,j)  
   
                do k=ktop(i,j,kc), kbtm(i,j,kc)
	          camtsw(i,j,k) = cldamt(i,j,kc) 
	          cirabsw(i,j,k,:) = cirab(i,j,kc)
	          cirrfsw(i,j,k,:) = cirrf(i,j,kc)
	          cvisrfsw(i,j,k,:) = cuvrf(i,j,kc)
                end do
                do k=ktop(i,j,kc), kbtm(i,j,kc)
	          if (ktop(i,j,kc) == kbtm(i,j,kc)) then
	            crndlw(i,j,k) = cldamt(i,j,kc)
	     	    cmxolw(i,j,k) = 0.0             
		    emrndlw(i,j,k,:) = emcld(i,j,kc) 
                  else
		    cmxolw(i,j,k) = cldamt(i,j,kc)
		    crndlw(i,j,k) = 0.0
		    emmxolw(i,j,k,:) = emcld(i,j,kc) 
                  endif
	        end do
	        if (ktop(i,j,kc) == kbtm(i,j,kc)) then
	          nrndlw(i,j) = nrndlw(i,j) + 1
                else
	          nmxolw(i,j) = nmxolw(i,j) + 1
                endif
              end do
	    end do
          end do

	  deallocate (ktop  )
	  deallocate (kbtm  )
	  deallocate (cldamt)
	  deallocate (cuvrf )
	  deallocate (cirrf )
	  deallocate (cuvab )
	  deallocate (cirab )
	  deallocate (emcld )

!-----------------------------------------------------------------
! (if ierr /= 0, then default clouds will be used. this is ok on the
! first timestep, when this condition happens, but should likely be 
! changed to indicate a real error if it occurs after the first step.
!--------------------------------------------------------------------
	else 
        endif  ! (ierr)

        deallocate (temp   )
        deallocate (qmix   )
        deallocate (rhum   )
        deallocate (omega  )
        deallocate (lgscldelq)
        deallocate (cnvcntq)
        deallocate (convprc)
        deallocate (cosangsolar  )
	deallocate (lat)

      else

        call error_mesg ('diag_clouds_W',   &
          'diag_cloud only available in FMS', FATAL)

      endif ! (running_fms)




end subroutine diag_clouds_calc


!####################################################################


	       end module diag_clouds_W_mod



