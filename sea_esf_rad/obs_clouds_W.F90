
                 module obs_clouds_W_mod

use cloud_zonal_mod,     only: getcld
use cloud_obs_mod,       only: cloud_obs, cloud_obs_init
use utilities_mod,       only: open_file, file_exist,   &
                               check_nml_error, error_mesg,   &
                               print_version_number, FATAL, NOTE, &
			       WARNING, get_my_pe, close_file
use rad_step_setup_mod,  only: jabs, iabs,  &
                               ISRAD, IERAD, JSRAD, JERAD, & 
                               KSRAD, KERAD,  Rad_time_sv,  &
			       lat_sv, pflux
use rad_utilities_mod,   only: Environment, environment_type
use constants_new_mod,   only: pstd_mks

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!	   observed cloud radiative properties module
!       currently this module is a wrapper until SKYHI goes away
!       and this module can be consolidated with obs_clouds_mod
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: obs_clouds_W.F90,v 1.2 2001/08/30 15:11:56 fms Exp $'
   character(len=128)  :: tag     =  '$Name: eugene $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
	  obs_clouds_init, obs_clouds_calc

!---------------------------------------------------------------------
!-------- namelist  ---------


integer   :: dummy = 0




namelist /obs_clouds_W_nml /     &
			       dummy


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

integer, parameter  :: NOFCLDS_SP=3  ! total number of clouds per column
integer, parameter  :: NOFMXOLW=1    ! number of max overlap clouds  
integer, parameter  :: NOFRNDLW=2    ! number of random overlap clouds
 
!--------------------------------------------------------------------
!   specified low, middle, high cloud properties (index = 1, 2, 3)
!   crfvis   :  visible band reflectivity
!   crfir    :  near-ir band reflectivity
!   cldem    :  infrared emissivity
!   cabir    :  near-ir band absorptivity
!--------------------------------------------------------------------
 
real, dimension(NOFCLDS_SP)    :: crfvis_m, crfir_m,    &
				  crfvis_fms, crfir_fms,   & 
				  crfvis, crfir,  &
                                  cabir, cldem

data crfvis_m    / 0.69E+00, 0.48E+00, 0.21E+00 / 
data crfir_m     / 0.69E+00, 0.48E+00, 0.21E+00 / 
data crfvis_fms  / 0.59E+00, 0.45E+00, 0.21E+00 /
data crfir_fms   / 0.59E+00, 0.45E+00, 0.21E+00 /
data cldem       / 1.00E+00, 1.00E+00, 1.00E+00 /  
data cabir       / 0.35E-01, 0.02E+00, 0.05E-01 / 

!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





subroutine obs_clouds_init (lonb, latb)

!---------------------------------------------------------------------
real, dimension(:), intent(in)   :: lonb, latb
!---------------------------------------------------------------------


      integer             :: unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------
  
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=obs_clouds_W_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'obs_clouds_W_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'obs_clouds_W', version_number)
      if (get_my_pe() == 0) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=obs_clouds_W_nml)
      endif
      call close_file (unit)

!--------------------------------------------------------------------
!   define cloud reflectivities dependent on the set of peripherals
!   selected
!--------------------------------------------------------------------
      if (Environment%using_fms_periphs) then
        crfvis(:) = crfvis_fms(:)
        crfir (:) = crfir_fms (:)
      else if (Environment%using_sky_periphs) then
        crfvis(:) = crfvis_m(:)
        crfir (:) = crfir_m (:)
      endif

      call cloud_obs_init (lonb, latb)

end subroutine obs_clouds_init



!####################################################################

subroutine obs_clouds_calc (camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
			    nrndlw, cirabsw, cvisrfsw, cirrfsw, &
			    emmxolw, emrndlw, is, ie, js, je)

!---------------------------------------------------------------------
integer,                     intent(in)     :: is, ie, js, je
integer, dimension(:,:),     intent(inout)  :: ncldsw, nmxolw, nrndlw
real,    dimension(:,:,:),   intent(inout)  :: camtsw, cmxolw, crndlw
real,    dimension(:,:,:,:), intent(inout)  :: cirabsw, cvisrfsw, &
					       cirrfsw, emmxolw, emrndlw
 
!---------------------------------------------------------------------
!   these arrays define the basic radiative properties of the clouds.

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
!--------------------------------------------------------------------

      real, dimension(:,:,:), allocatable    :: camtsw3
      integer, dimension(:,:,:), allocatable :: ktopsw3, kbtmsw3
      real, dimension(:,:,:), allocatable    :: phaf
      integer                                :: k, j, i


!-----------------------------------------------------------------------
!     Move property arrays for large-scale ice clouds from physical 
!     space into cloud space.  assign same emissivity properties to all 
!     freq bands. (will be replaced later). shortwave cloud amount is 
!     assumed to be the sum of the two longwave cloud amounts.
!-----------------------------------------------------------------------

       nmxolw(:,:) = NOFMXOLW
       nrndlw(:,:) = NOFRNDLW
       ncldsw(:,:) = NOFCLDS_SP

       allocate (phaf  (ISRAD:IERAD, JSRAD:JERAD, KERAD+1) )
       allocate (camtsw3 (ISRAD:IERAD, JSRAD:JERAD, NOFCLDS_SP) )
       allocate (ktopsw3 (ISRAD:IERAD, JSRAD:JERAD, NOFCLDS_SP) )
       allocate (kbtmsw3 (ISRAD:IERAD, JSRAD:JERAD, NOFCLDS_SP) )

       do k=KSRAD,KERAD+1
         phaf(:,:,k) = pflux(:,:,k)*pstd_mks/pflux(:,:,KERAD+1)
       end do

       call getcld (Rad_time_sv, lat_sv, phaf, ktopsw3, kbtmsw3, &
	            camtsw3 )

       call cloud_obs(is, js, Rad_time_sv, camtsw3)

!---------------------------------------------------------------------
!    map the cloud-space arrays obtained above to model space arrays. 
!    define the remaining cloud properties from the data provided 
!    within this module.
!-------------------------------------------------------------------
       do j=JSRAD,JERAD
         do i=ISRAD,IERAD

!--------------------------------------------------------------------
!    high clouds
!--------------------------------------------------------------------
           do k=ktopsw3(i,j,1), kbtmsw3(i,j,1)-1
	     camtsw(i,j,k) = camtsw3(i,j,1)
             cirabsw(i,j,k,:)  = cabir(3)
  	     cirrfsw(i,j,k,:)  = crfir(3)
	     cvisrfsw(i,j,k,:) = crfvis(3)
           enddo
           do k=ktopsw3(i,j,1), kbtmsw3(i,j,1)-1
	     crndlw(i,j,k) = camtsw3(i,j,1)
             emrndlw(i,j,k,:)  = cldem(3)
           enddo

!--------------------------------------------------------------------
!    middle clouds
!--------------------------------------------------------------------
           do k=ktopsw3(i,j,2), kbtmsw3(i,j,2)-1
	     camtsw(i,j,k) = camtsw3(i,j,2)
  	     cirabsw(i,j,k,:)  = cabir(2)
  	     cirrfsw(i,j,k,:)  = crfir(2)
	     cvisrfsw(i,j,k,:) = crfvis(2)
           enddo
           do k=ktopsw3(i,j,2), kbtmsw3(i,j,2)-1
	     crndlw(i,j,k) = camtsw3(i,j,2)
             emrndlw(i,j,k,:)  = cldem(2)
           enddo

!--------------------------------------------------------------------
!    low clouds
!--------------------------------------------------------------------
           do k=ktopsw3(i,j,3), kbtmsw3(i,j,3)-1
	     camtsw(i,j,k) = camtsw3(i,j,3)
  	     cirabsw(i,j,k,:)  = cabir(1)
  	     cirrfsw(i,j,k,:)  = crfir(1)
	     cvisrfsw(i,j,k,:) = crfvis(1)
           enddo
           do k=ktopsw3(i,j,3), kbtmsw3(i,j,3)-1
	     cmxolw(i,j,k) = camtsw3(i,j,3)
  	      emmxolw(i,j,k,:)  = cldem(1)                 
           end do

	 end do
       end do

       deallocate (phaf)
       deallocate (camtsw3)
       deallocate (ktopsw3)
       deallocate (kbtmsw3)


end subroutine obs_clouds_calc





!####################################################################

	       end module obs_clouds_W_mod



