
                 module diag_clouds_W_mod

use time_manager_mod,       only: time_type
use diag_cloud_mod,         only: diag_cloud_avg, diag_cloud_driver
use utilities_mod,          only: open_file, file_exist,   &
                                  check_nml_error, error_mesg,   &
                                  print_version_number, FATAL, NOTE, &
				  WARNING, get_my_pe, close_file
!use rad_step_setup_mod,     only: temp, press, pflux, jabs, iabs,  &
!use rad_step_setup_mod,     only: temp, press, pflux,              &
!use rad_step_setup_mod,     only: temp, press,             &
!use rad_step_setup_mod,     only: temp, press
!                                 ISRAD, IERAD, JSRAD, JERAD, & 
!                                 KSRAD, KERAD, land, cosz
!                                               land, cosz
!                                                      cosz
use rad_utilities_mod,      only: Environment, environment_type, &
                                  longwave_control_type, Lw_control, &
				  cld_diagnostics_type,  &
                                  shortwave_control_type, Sw_control,&
                                  cloudrad_control_type, Cldrad_control
!use astronomy_package_mod,  only: get_astronomy_for_clouds,  &
!use astronomy_package_mod,  only:                            &
!			          get_astronomy_for_clouds_init
use microphys_rad_mod,      only: microphys_rad_driver, &
                                  lwemiss_calc, &
                                  microphys_rad_init

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
   character(len=128)  :: version =  '$Id: diag_clouds_W.F90,v 1.3 2002/07/16 22:34:44 fms Exp $'
   character(len=128)  :: tag     =  '$Name: havana $'



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

!integer   :: nsolwg
!character(len=10)     :: swform
!character(len=16)     :: swform
logical  :: do_lhsw, do_esfsw
logical               :: do_lwcldemiss

!real, dimension(:), allocatable :: latsv

!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





!subroutine diag_clouds_init ( ix, iy, kx, th, ierr)
!ubroutine diag_clouds_init ( ix, iy, kx,     ierr)
subroutine diag_clouds_init ( ix, iy,         ierr)

!integer,       intent(in)       :: ix, iy, kx
integer,       intent(in)       :: ix, iy
integer,       intent(out)      :: ierr
!real, dimension(:), intent(in)  :: th

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

!      swform = Sw_control%sw_form
      do_lhsw = Sw_control%do_lhsw
      do_esfsw = Sw_control%do_esfsw
      do_lwcldemiss = Lw_control%do_lwcldemiss

      if (Environment%running_gcm ) then
!       if (trim(swform) == 'esfsw99' .or. do_lwcldemiss) then
        if (do_esfsw                  .or. do_lwcldemiss) then
          call microphys_rad_init
        endif
      endif

!--------------------------------------------------------------------
!  retrieve module variables that come from other modules
!--------------------------------------------------------------------

!     call get_astronomy_for_clouds_init (nsolwg)
!     nsolwg = 1

!--------------------------------------------------------------------
!  define and save latitude for later use
!--------------------------------------------------------------------
!      allocate (latsv(size(th,1) ) )
!      latsv(:) = th(:)


end subroutine diag_clouds_init



!#################################################################

subroutine diag_clouds_calc ( is, ie, js, je, Cld_diagnostics, lat,  pflux, deltaz, &
                      cosz,  press, temp_in,     camtsw, cmxolw,  &
                            crndlw, ncldsw, &
			       nmxolw, nrndlw,    &
!		       cirabsw, &
!		       cvisrfsw, cirrfsw,   &
			       emmxolw, emrndlw, &
!			       is, js, Time_next,  &
			               Time_next,  &
			       cirabsw, cirrfsw, cvisrfsw,   &
			       cldext, cldsct, cldasymm)

integer,                      intent(in)    :: is, ie, js, je
real, dimension(:,:),         intent(in)    :: lat, cosz
real, dimension(:,:,:),         intent(in)    :: pflux, deltaz, press, &
                                                 temp_in
type(time_type),              intent(in)    :: Time_next
type(cld_diagnostics_type),   intent(inout)    :: Cld_diagnostics
integer, dimension(:,:),      intent(inout) :: ncldsw, nmxolw, nrndlw
real,    dimension(:,:,:),    intent(inout) :: camtsw, cmxolw, crndlw
real,    dimension(:,:,:,:),  intent(inout) :: emmxolw, emrndlw
real,    dimension(:,:,:,:),  intent(inout), optional :: cirabsw,   &
                                              cvisrfsw,  cirrfsw,   &
                                              cldsct, cldext, cldasymm

!--------------------------------------------------------------------
!   these arrays define the basic radiative properties of the clouds.
!
!     cmxolw  =  amounts of maximally overlapped longwave clouds in
!                layers from 1     to kx      
!     crndlw  =  amounts of randomly overlapped longwave clouds in
!                layers from 1     to kx   .
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
!     real, dimension(:,:,:), allocatable    :: cosangsolar
      real, dimension(:,:  ), allocatable    :: cosangsolar
      integer, dimension(:,:,:), allocatable :: kbtm, ktop
      real, dimension(:,:,:), allocatable    :: temp, qmix, rhum, &
						omega, lgscldelq, &
						cnvcntq
!      real, dimension(:,:), allocatable      :: convprc,lat
      real, dimension(:,:), allocatable      :: convprc
      real, dimension(:,:,:), allocatable    :: ql, qi, cf , &
                                                cldamt, cuvab, cirab, &
                                                cuvrf, cirrf, emcld, &
                                                conc_drop, conc_ice, &
                                                size_drop, size_ice
      real, dimension(:,:,:,:), allocatable    :: abscoeff, cldemiss
      integer                                :: j, i, k
      integer                                :: ierr, kc
      integer         :: ix, jx, kx 

      ix = size (camtsw, 1)
      jx = size(camtsw,2)
      kx = size(camtsw,3)


      if (Environment%running_gcm) then
!--------------------------------------------------------------------
!     if (Environment%running_fms) then

!---------------------------------------------------------------------
!     obtain the appropriate zenith angles that are to be used here.
!---------------------------------------------------------------------
!       allocate (cosangsolar (ISRAD:IERAD, JSRAD:JERAD, nsolwg) )
!       allocate (cosangsolar (ISRAD:IERAD, JSRAD:JERAD,      1) )
!       allocate (cosangsolar (ISRAD:IERAD, JSRAD:JERAD        ) )
        allocate (cosangsolar (ix, jx                          ) )
!       call get_astronomy_for_clouds (cosangsolar   )
!       cosangsolar(:,:,:) = cosz(:,:,:)
        cosangsolar(:,:  ) = cosz(:,:  )

!---------------------------------------------------------------------
!     obtain the cloud amounts and areas.
!---------------------------------------------------------------------
!       allocate (temp (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!       allocate (qmix (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!       allocate (rhum (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!       allocate (omega (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!       allocate (lgscldelq (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!       allocate (cnvcntq   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!       allocate (convprc   (ISRAD:IERAD, JSRAD:JERAD) )
!       allocate (lat       (ISRAD:IERAD, JSRAD:JERAD) )

!!!! USE input temp if diag_cloud_avg were not being called
!!! it appears that this will not happen currently -- diag_cloud_avg
!!! will return a value unless ierr == 1, in which case it is not
!! needed here.

        allocate (temp (ix, jx, kx                           ) )
         temp(:,:,:) = temp_in(:,:,1:kx)



        allocate (qmix (ix, jx, kx                           ) )
        allocate (rhum (ix, jx, kx                           ) )
        allocate (omega (ix, jx, kx                           ) )
        allocate (lgscldelq (ix, jx, kx                           ) )
        allocate (cnvcntq   (ix, jx, kx                           ) )
        allocate (convprc   (ix, jx                  ) )
!        allocate (lat       (ix, jx                  ) )

!       call diag_cloud_avg (iabs(ISRAD), jabs(JSRAD), temp, qmix, &
        call diag_cloud_avg (is, js, temp, qmix, &
			     rhum, omega, lgscldelq, cnvcntq, convprc, &
			     ierr)

!---------------------------------------------------------------------
!     allocate and initialize the cloud radiative property arrays.
!---------------------------------------------------------------------
        if (ierr == 0) then
!         allocate (ktop   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!  allocate (kbtm   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!  allocate (cldamt (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!  allocate (emcld  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
          allocate (ktop   (ix, jx, kx                           ) )
	  allocate (kbtm   (ix, jx, kx                           ) )
	  allocate (cldamt (ix, jx, kx                           ) )
	  allocate (emcld  (ix, jx, kx                           ) )
!        if (trim(swform) == 'lhsw' ) then
         if (do_lhsw                ) then
!  allocate (cuvrf  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!  allocate (cirrf  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!  allocate (cuvab  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!  allocate (cirab  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (cuvrf  (ix, jx, kx                           ) )
	  allocate (cirrf  (ix, jx, kx                           ) )
	  allocate (cuvab  (ix, jx, kx                            ) )
	  allocate (cirab  (ix, jx, kx                           ) )
         endif
!       if (trim(swform) == 'esfsw99' .or. do_lwcldemiss) then
        if (do_esfsw                  .or. do_lwcldemiss) then
!         allocate (conc_drop  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
!         allocate (conc_ice  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!         allocate (size_drop (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
!         allocate (size_ice  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
   allocate (abscoeff ( ix,jx,kx, SIZE(emmxolw,4)) )
   allocate (cldemiss ( ix,jx,kx, SIZE(emmxolw,4)) )
          allocate (conc_drop  (ix, jx, kx                           ))
          allocate (conc_ice  (ix, jx, kx                           ) )
          allocate (size_drop (ix, jx, kx                           ))
          allocate (size_ice  (ix, jx, kx                           ) )
        endif

	  emcld = 1.
	  ktop=1
	  kbtm=0
	  cldamt=0.
!  if (trim(swform) == 'lhsw' ) then
	  if (do_lhsw                ) then
            cuvrf=0.
            cirrf=0.
            cuvab=0.
            cirab=0.
          endif
!       if (trim(swform) == 'esfsw99' .or. do_lwcldemiss) then
        if (do_esfsw                  .or. do_lwcldemiss) then
	    conc_ice(:,:,:) = 0.
	    conc_drop(:,:,:) = 0.
	    size_ice(:,:,:) = 60.
	    size_drop(:,:,:) = 20.
         end if

!         do j=jsrad,jerad
!         do j=1,jx           
!    lat(:,j) = latsv(jabs(j))
!  end do

!---------------------------------------------------------------------
!     obtain the cloud radiative properties.
!---------------------------------------------------------------------
!       if (trim(swform) == 'lhsw' ) then
        if (do_lhsw                ) then
          if ( .not. do_lwcldemiss) then
          call diag_cloud_driver (is, js, &
		temp, qmix, rhum, omega, lgscldelq, cnvcntq, convprc, &
!	0.1*press(:,:,KSRAD:KERAD),&
                0.1*press(:,:,1    :kx   ),&
!  	0.1*pflux, 0.1*pflux(:,:,KERAD+1),  &
	  	0.1*pflux, 0.1*pflux(:,:,kx   +1),  &
!		      cosangsolar(:,:,1), lat, Time_next, &
			      cosangsolar       , lat, Time_next, &
			      ncldsw, ktop, kbtm, cldamt, r_uv=cuvrf,  &
                              r_nir=cirrf, ab_uv=cuvab,   &
			      ab_nir=cirab, em_lw=emcld)
          else          
              call diag_cloud_driver (is, js, &
                 temp, qmix, rhum, omega, lgscldelq, cnvcntq, convprc, &
!                0.1*press(:,:,KSRAD:KERAD),&
                 0.1*press(:,:,    1:kx   ),&
!                0.1*pflux, 0.1*pflux(:,:,KERAD+1),  &
                 0.1*pflux, 0.1*pflux(:,:,kx   +1),  &
!                cosangsolar(:,:,1), lat, Time_next, &
                 cosangsolar       , lat, Time_next, &
                 ncldsw, ktop, kbtm, cldamt, r_uv=cuvrf,  &
                 r_nir=cirrf, ab_uv=cuvab, ab_nir=cirab, em_lw=emcld, &
		 conc_drop=conc_drop, conc_ice=conc_ice, &
		 size_drop=size_drop, size_ice=size_ice)
           endif
!---------------------------------------------------------------------
!    map the cloud-space arrays to physical space arrays
!-------------------------------------------------------------------
!         do j=JSRAD,JERAD
!           do i=ISRAD,IERAD
          do j=1,jx           
            do i=1,ix           
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


!    this caters to the case when the shortwave is lh but the
!    longwave is sea_esf with the do_lwcldemiss option on
          if (do_lwcldemiss) then
            call microphys_rad_driver (is, ie, js, je, deltaz, &
!              press, temp_in, camtsw,  &
!                          emmxolw, emrndlw, Cld_diagnostics, &
	              press, temp_in,          &
	                                            Cld_diagnostics, &
                            abscoeff=abscoeff, &
                                  conc_drop_in=conc_drop,     &
                                   conc_ice_in=conc_ice,      &
                                      size_drop_in=size_drop,    &
                                    size_ice_in=size_ice)

           call lwemiss_calc(  deltaz,                            &
                           abscoeff,                         &
                         cldemiss)

!   as of 3 august 2001, we are assuming randomly overlapped clouds
!  only. the cloud and emissivity settings follow:
       emmxolw = cldemiss
        emrndlw = cldemiss

          endif

 
!        else if (trim(swform) == 'esfsw99') then
         else if (do_esfsw                 ) then

              call diag_cloud_driver (is, js, &
                 temp, qmix, rhum, omega, lgscldelq, cnvcntq, convprc, &
!                0.1*press(:,:,KSRAD:KERAD),&
                 0.1*press(:,:,    1:kx   ),&
!                0.1*pflux, 0.1*pflux(:,:,KERAD+1),  &
                 0.1*pflux, 0.1*pflux(:,:,kx   +1),  &
!                cosangsolar(:,:,1), lat, Time_next, &
                 cosangsolar       , lat, Time_next, &
                 ncldsw, ktop, kbtm, cldamt,   &
!	 cuvrf,  &
!                cirrf, cuvab, cirab, emcld)
                 conc_drop=conc_drop, conc_ice=conc_ice, &
		 size_drop=size_drop, size_ice=size_ice)

!---------------------------------------------------------------------
!   as of 14 dec 2000, all cloud properties and cloud fractions are
!   assumed to be for randomly overlapped clouds
!---------------------------------------------------------------------
!         do j=JSRAD,JERAD
!           do i=ISRAD,IERAD
          do j=1,jx         
            do i=1,ix           
              do kc=1, ncldsw(i,j)  
   
                do k=ktop(i,j,kc), kbtm(i,j,kc)
	          camtsw(i,j,k) = cldamt(i,j,kc) 
		  crndlw(i,j,k) = camtsw(i,j,k)
		  cmxolw(i,j,k) = 0.0
		  emrndlw(i,j,k,1) = emcld(i,j,kc)
		  emmxolw(i,j,k,:) = 1.0
                 end do
                 end do
                 end do
                 end do

!         crndlw = cldamt
!         cmxolw = 0.0E+00
!        camtsw = cmxolw + crndlw
!         emmxolw = 1.0
!         emrndlw(:,:,:,1) = emcld
!
!--------------------------------------------------------------------
!  if microphysics is being used for either sw or lw calculation, call
!  microphys_rad_driver to obtain cloud radiative properties.
!--------------------------------------------------------------------
       call microphys_rad_driver (is, ie, js, je, deltaz, &
!                              press, temp_in,   camtsw, emmxolw,  &
!                                emrndlw, Cld_diagnostics, &
                               press, temp_in,                     &
                                          Cld_diagnostics, &
                                  cldext=cldext, cldsct=cldsct,    &
                            abscoeff=abscoeff, &
!                                 cldasymm=cldasymm)             
                                  cldasymm=cldasymm,            &
                                  conc_drop_in=conc_drop,     &
                                  conc_ice_in=conc_ice,      &
                                  size_drop_in=size_drop,    &
                                  size_ice_in=size_ice)


!       if (do_lwcldemiss) then
           call lwemiss_calc(  deltaz,                            &
                           abscoeff,                         &
                         cldemiss)

!   as of 3 august 2001, we are assuming randomly overlapped clouds
!  only. the cloud and emissivity settings follow:
       emmxolw = cldemiss
        emrndlw = cldemiss

!     endif


     endif
	  deallocate (ktop  )
	  deallocate (kbtm  )
	  deallocate (cldamt)
	  deallocate (emcld )
!         if (trim(swform) == 'lhsw' ) then
          if (do_lhsw                ) then
	    deallocate (cuvrf )
            deallocate (cirrf )
            deallocate (cuvab )
            deallocate (cirab )
          endif
!        if (trim(swform) == 'esfsw99' .or. do_lwcldemiss ) then
         if (do_esfsw                  .or. do_lwcldemiss ) then
         deallocate (conc_drop )
         deallocate (conc_ice  )
         deallocate (size_drop )
         deallocate (size_ice  )
         deallocate (abscoeff  )
         deallocate (cldemiss  )
        endif

!-----------------------------------------------------------------
! (if ierr /= 0, then default clouds will be used. this is ok on the
! first timestep, when this condition happens, but should likely be 
! changed to indicate a real error if it occurs after the first step.
!--------------------------------------------------------------------
	else 
! needed to supply values for radiation_diag_mod
        if (do_esfsw                  .or. do_lwcldemiss) then
          allocate (conc_drop  (ix, jx, kx                           ))
          allocate (conc_ice  (ix, jx, kx                           ) )
          allocate (size_drop (ix, jx, kx                           ))
          allocate (size_ice  (ix, jx, kx                           ) )
	    conc_ice(:,:,:) = 0.
	    conc_drop(:,:,:) = 0.
	    size_ice(:,:,:) = 60.
	    size_drop(:,:,:) = 20.
       call microphys_rad_driver (is, ie, js, je, deltaz, &
!                              press, temp_in,   camtsw, emmxolw,  &
!                                emrndlw, Cld_diagnostics, &
                               press, temp_in,                     &
                                          Cld_diagnostics, &
                                  cldext=cldext, cldsct=cldsct,    &
!                                abscoeff=abscoeff, &
!                                 cldasymm=cldasymm)             
                                  cldasymm=cldasymm,            &
                                  conc_drop_in=conc_drop,     &
                                  conc_ice_in=conc_ice,      &
                                  size_drop_in=size_drop,    &
                                  size_ice_in=size_ice)
         deallocate (conc_drop )
         deallocate (conc_ice  )
         deallocate (size_drop )
         deallocate (size_ice  )
         emrndlw = 1.0
         emmxolw = 1.0
        endif
        endif  ! (ierr)

        deallocate (temp   )
        deallocate (qmix   )
        deallocate (rhum   )
        deallocate (omega  )
        deallocate (lgscldelq)
        deallocate (cnvcntq)
        deallocate (convprc)
        deallocate (cosangsolar  )
!	deallocate (lat)

!     else

!       call error_mesg ('diag_clouds_W',   &
!         'diag_cloud only available in FMS', FATAL)

!     endif ! (running_fms)

     endif  ! (running_gcm)



end subroutine diag_clouds_calc


!####################################################################


	       end module diag_clouds_W_mod



