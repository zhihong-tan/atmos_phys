
                 module strat_clouds_W_mod

use time_manager_mod,       only: time_type
use strat_cloud_mod,        only: strat_cloud_avg
use cloud_rad_mod,          only: cloud_summary
use utilities_mod,          only: open_file, file_exist,   &
                                  check_nml_error, error_mesg,   &
                                  print_version_number, FATAL, NOTE, &
				  WARNING, get_my_pe, close_file
use constants_new_mod, only :     rh2oair
!use rad_step_setup_mod,     only: temp, rh2o, press, pflux, jabs,   &
!use rad_step_setup_mod,     only: temp, rh2o, press, pflux,         &
!use rad_step_setup_mod,     only: temp, rh2o, press,                &
!use rad_step_setup_mod,     only: temp, rh2o, press, &
!use rad_step_setup_mod,     only: temp, rh2o, press
!			  iabs, ISRAD, IERAD, JSRAD, JERAD, & 
!                                 KSRAD, KERAD, land, cosz, &
!				  iabs,                            & 
!                                                land, cosz, &
!                                                      cosz, &
!				  cloud_ice, cloud_water
use rad_utilities_mod,      only: Environment, environment_type, &
                                  longwave_control_type, Lw_control, &
				  shortwave_control_type, Sw_control,&
				  cld_diagnostics_type, &
                                  cloudrad_control_type, Cldrad_control
!use astronomy_package_mod,  only: get_astronomy_for_clouds,  &
!use astronomy_package_mod,  only:                            &
!		          get_astronomy_for_clouds_init
use microphys_rad_mod,      only: microphys_rad_driver, &
                                  lwemiss_calc
!!! DOES THIS NEED TO BE CALLED FROM HERE ???
!                                  microphys_rad_init

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!	          strat cloud radiative properties module
!            currently a wrapper until SKYHI goes away and this
!            module can be consolidated with strat_cloud_mod
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: strat_clouds_W.F90,v 1.3 2002/07/16 22:37:10 fms Exp $'
   character(len=128)  :: tag     =  '$Name: havana $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          strat_clouds_init, strat_clouds_calc

!---------------------------------------------------------------------
!-------- namelist  ---------

integer   :: dummy = 0



namelist /strat_clouds_W_nml /     &
			       dummy


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!integer               :: nsolwg
!character(len=10)     :: swform
!character(len=16)     :: swform
logical        :: do_lhsw, do_esfsw
logical               :: do_lwcldemiss


!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





subroutine strat_clouds_init 

      integer            :: unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------
  
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=strat_clouds_W_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'strat_clouds_W_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'strat_clouds_W', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=strat_clouds_W_nml)
      endif
      call close_file (unit)

!      swform = Sw_control%sw_form
      do_lhsw = Sw_control%do_lhsw
      do_esfsw = Sw_control%do_esfsw
      do_lwcldemiss = Lw_control%do_lwcldemiss

!     if (Environment%running_gcm ) then
!       if (do_esfsw                  .or. do_lwcldemiss) then
!         call microphys_rad_init
!       endif
!     endif

!--------------------------------------------------------------------
!  retrieve module variables that come from other modules
!--------------------------------------------------------------------

!     call get_astronomy_for_clouds_init (nsolwg)
!     nsolwg = 1


end subroutine strat_clouds_init



!#################################################################

subroutine strat_clouds_calc (is, ie, js, je, Cld_diagnostics, pflux,  deltaz,  &
                             land, cosz, cloud_water, cloud_ice,  &
			     press, temp, rh2o,   &
			     cldamt_out, &
			     camtsw, cmxolw,   &
                              crndlw, ncldsw, &
                              nmxolw, nrndlw,                    &
                              emmxolw, emrndlw,         &
                              Time_next,                &
                              cirabsw, cirrfsw, cvisrfsw, cldext, &  
                              cldsct, cldasymm, abscoeff)

integer, intent(in) :: is, ie, js, je
real,    dimension(:,:,:),   intent(in) :: pflux, deltaz, cloud_water, &
                                           cloud_ice, press, temp, rh2o
real,    dimension(:,:),   intent(in) :: land, cosz
real,    dimension(:,:,:),   intent(inout), optional :: camtsw, cmxolw, crndlw
real,    dimension(:,:,:),   intent(out), optional :: cldamt_out            
real,    dimension(:,:,:,:), intent(inout), optional ::   &
                                               emmxolw, emrndlw
type(cld_diagnostics_type), intent(inout) :: Cld_diagnostics
type(time_type),                 intent(in), optional ::  Time_next
real,    dimension(:,:,:,:), intent(inout), optional  ::    &
                                             cirabsw, cvisrfsw, &
                                             cirrfsw
real,    dimension(:,:,:,:), intent(out), optional  ::    &
                                             cldext, cldsct, cldasymm,&
                                             abscoeff
integer, dimension(:,:),     intent(inout), optional :: ncldsw, nmxolw, nrndlw


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
!     abscoeff = longwave absorption coefficient in (NEMISS_BDS)
!                frequency bands. may be used to obtain emissivity.
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
      real, dimension(:,:,:), allocatable    :: ql, qi, cf , &
                                                cldamt, cuvab, cirab, &
						cuvrf, cirrf, emcld,&
						conc_drop, conc_ice, &
                                                size_drop, size_ice, &
                                                sp_hum
     real, dimension(:,:,:,:), allocatable ::  abscoeff_lcl,      &
                                               cldemiss_lcl

      integer                                :: j, i, k
      integer                                :: ierr, kc
      integer  :: israd, ierad, jsrad, jerad, ksrad, kerad

      israd = 1
!     ierad = size (camtsw, 1)
      ierad = size (press , 1)
      jsrad = 1
      jerad = size (press , 2)
      ksrad = 1
      kerad = size (press , 3) - 1

     if (Environment%running_gcm) then
!--------------------------------------------------------------------
!     if (Environment%running_fms) then

!---------------------------------------------------------------------
!     obtain the appropriate zenith angles that are to be used here.
!---------------------------------------------------------------------
!       allocate (cosangsolar (ISRAD:IERAD, JSRAD:JERAD, nsolwg) )
!       allocate (cosangsolar (ISRAD:IERAD, JSRAD:JERAD,      1) )
        allocate (cosangsolar (ISRAD:IERAD, JSRAD:JERAD        ) )
!       call get_astronomy_for_clouds (cosangsolar   )
!       cosangsolar(:,:,:) = cosz(:,:,:)
        cosangsolar(:,:) = cosz(:,:)

!---------------------------------------------------------------------
!     obtain the cloud amounts and areas.
!---------------------------------------------------------------------
         allocate (ql (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
         allocate (qi (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
         allocate (cf (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!        call strat_cloud_avg (iabs(ISRAD), jabs(JSRAD), ql, qi,   &
         call strat_cloud_avg (is, js,                   ql, qi,   &
                               cf, ierr)
  
!---------------------------------------------------------------------
!     allocate and initialize the cloud radiative property arrays.
!---------------------------------------------------------------------
         if (ierr == 0) then
           allocate (ktop   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
           allocate (kbtm   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
           allocate (cldamt (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
           allocate (emcld  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
           allocate (sp_hum (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!       if (trim(swform) == 'lhsw' ) then
        if (do_lhsw                ) then
	  allocate (cuvrf  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (cirrf  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (cuvab  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	  allocate (cirab  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	endif
!       if (trim(swform) == 'esfsw99' .or. do_lwcldemiss) then
        if (do_esfsw                  .or. do_lwcldemiss) then
           allocate (conc_drop  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
           allocate (conc_ice  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
           allocate (size_drop (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
           allocate (size_ice  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
	 endif
 
          emcld = 1.
          ktop=1
          kbtm=0
          cldamt=0.
!       if (trim(swform) == 'lhsw' ) then
        if (do_lhsw                ) then
	  cuvrf=0.
	  cirrf=0.
	  cuvab=0.
	  cirab=0.
        endif
!---------------------------------------------------------------------
!     define specific humidity to pass to start_cloud_mod.
!---------------------------------------------------------------------
! rh2o (  :  ,  :  ,:) = q(:,:,:)/(1.0-q(:, :, :)/rh2oair)
       sp_hum(:,:,:) = rh2o(:,:,:)/(1.0+rh2o(:,:,:)/rh2oair)
  
!---------------------------------------------------------------------
!     obtain the cloud radiative properties.
!---------------------------------------------------------------------
!      if (trim(swform) == 'lhsw' ) then
       if (do_lhsw                ) then
	 if ( .not. do_lwcldemiss) then
 
!        call cloud_summary (iabs(ISRAD), jabs(JSRAD), &
         call cloud_summary (is, js,                   &
!                            land, ql, qi, cf, rh2o,   &
                             land, ql, qi, cf, sp_hum, &
                             0.1*press(:,:,KSRAD:KERAD),&
                             0.1*pflux, temp(:,:,KSRAD:KERAD),  &
!                            cosangsolar(:,:,1),&
                             cosangsolar       ,&
                             temp(:,:,KERAD+1),        &
                             ncldsw, ktop, kbtm, cldamt,         &
                             Time=Time_next,                 &
                             r_uv=cuvrf, r_nir=cirrf, ab_uv=cuvab, &
                             ab_nir=cirab, em_lw=emcld)
           else

!        call cloud_summary (iabs(ISRAD), jabs(JSRAD), &
         call cloud_summary (is, js,                   &
!                            land, ql, qi, cf, rh2o,   &
                             land, ql, qi, cf, sp_hum, &
                             0.1*press(:,:,KSRAD:KERAD),&
                             0.1*pflux, temp(:,:,KSRAD:KERAD),  &
!                            cosangsolar(:,:,1),&
                             cosangsolar       ,&
                             temp(:,:,KERAD+1),        &
                             ncldsw, ktop, kbtm, cldamt,         &
                             Time=Time_next,                 &
                             r_uv=cuvrf, r_nir=cirrf, ab_uv=cuvab, &
                             ab_nir=cirab, em_lw=emcld,  &
                              conc_drop=conc_drop, conc_ice=conc_ice, &
                              size_drop=size_drop, size_ice=size_ice)

	   endif

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

!    this caters to the case when the shortwave is lh but the
!    longwave is sea_esf with the do_lwcldemiss option on

!    in this case the routine returns the lw emissivities (emmxolw,  
!  emrndlw) not the abs coeff (abscoeff). abscoeff_lcl is just the
!  locally defined abs coeff, used for convenience in obtaining
!  emmxolw and emrndlw.

          if (do_lwcldemiss) then

!           call microphys_rad_driver (is, ie, js, je, deltaz, &
!              press, temp,     camtsw,  &
!                                    emmxolw, emrndlw, Cld_diagnostics, &
!                                    conc_drop_in=conc_drop,     &
!                                   conc_ice_in=conc_ice,      &
!                                    size_drop_in=size_drop,    &
!                                    size_ice_in=size_ice)

allocate (abscoeff_lcl ( SIZE(emmxolw,1), SIZE(emmxolw,2),SIZE(emmxolw,3), &
                SIZE(emmxolw,4)) )
 allocate (cldemiss_lcl ( SIZE(emmxolw,1), SIZE(emmxolw,2),SIZE(emmxolw,3), &
                  SIZE(emmxolw,4)) )
         call microphys_rad_driver (                           &
                  is,ie,js,je, deltaz,press,temp, Cld_diagnostics, &
                     conc_drop_in=conc_drop, conc_ice_in=conc_ice, &
                     size_drop_in=size_drop, size_ice_in=size_ice, &
                    abscoeff=abscoeff_lcl)
       call lwemiss_calc(   deltaz,                           &
                     abscoeff_lcl,                         &
                         cldemiss_lcl)

          emmxolw = cldemiss_lcl
     emrndlw = cldemiss_lcl
  if (ALLOCATED (abscoeff_lcl)) then
   deallocate (abscoeff_lcl)
    deallocate (cldemiss_lcl)
  endif

          endif



!       else if (trim(swform) == 'esfsw99') then
        else if (do_esfsw                 ) then


          call cloud_summary (is, js,                   &
                              land, ql, qi, cf, sp_hum, &
                              0.1*press(:,:,KSRAD:KERAD),&
                              0.1*pflux, temp(:,:,KSRAD:KERAD),  &
                              cosangsolar       ,&
                              temp(:,:,KERAD+1),        &
!                             ncldsw, ktop, kbtm, cldamt,         &
                              ncldsw, ktop, kbtm,  &
			      cldamt=cldamt_out,         &
                              Time=Time_next,                &
                              conc_drop=conc_drop, conc_ice=conc_ice, &
                              size_drop=size_drop, size_ice=size_ice)
  
!  ! DROP
!
!   as of 14 dec 2000, all cloud properties and cloud fractions are
!   assumed to be for randomly overlapped clouds
!
!        crndlw = cldamt
!        cmxolw = 0.0E+00
!        camtsw = cmxolw + crndlw
!        emmxolw = 1.0
!        emrndlw(:,:,:,1) = emcld
!!  END OF DROP
!
!--------------------------------------------------------------------
!  if microphysics is being used for either sw or lw calculation, call
!  microphys_rad_driver to obtain cloud radiative properties.
!--------------------------------------------------------------------
!         call microphys_rad_driver (is, ie, js, je, deltaz, &
!      press, temp,            camtsw, emmxolw, &
!                                     emrndlw, Cld_diagnostics, &
!                                    cldext=cldext, cldsct=cldsct,    &
!                                    cldasymm=cldasymm,            &
!                                    conc_drop_in=conc_drop,     &
!                                   conc_ice_in=conc_ice,      &
!                                    size_drop_in=size_drop,    &
!                                    size_ice_in=size_ice)

         call microphys_rad_driver (                           &
               is,ie,js,je,deltaz,press,temp, Cld_diagnostics, &
                  conc_drop_in=conc_drop, conc_ice_in=conc_ice, &
                         size_drop_in=size_drop, size_ice_in=size_ice, &
            cldext=cldext, cldsct=cldsct, cldasymm=cldasymm, &
                               abscoeff=abscoeff)

        endif

!! DROP
!       deallocate (ktop)
!       deallocate (kbtm)
!       deallocate (cldamt)
!       deallocate (emcld )
!       deallocate (sp_hum)
!       if (do_lhsw                ) then
!  deallocate (cuvrf )
!  deallocate (cirrf )
!  deallocate (cuvab )
!  deallocate (cirab )
!endif
!if (do_esfsw                  .or. do_lwcldemiss ) then
!       deallocate (conc_drop )
!       deallocate (conc_ice  )
!       deallocate (size_drop )
!       deallocate (size_ice  )
!endif
! the following needed for radiation_diag_mod on first call to strat_clouds
      else ! (ierr)
!  write error mesg after first call
        if (do_esfsw                  .or. do_lwcldemiss) then
           allocate (conc_drop  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
           allocate (conc_ice  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
           allocate (size_drop (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
           allocate (size_ice  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
            conc_drop = 0.
            conc_ice = 0.
	    size_drop = 2.02
	    size_ice = 2.02
            cldamt_out=0.
          call microphys_rad_driver (is, ie, js, je, deltaz, &
!      press, temp,            camtsw, emmxolw, &
!                                     emrndlw, Cld_diagnostics, &
	      press, temp,                             &
                                               Cld_diagnostics, &
                                     cldext=cldext, cldsct=cldsct,    &
                                     cldasymm=cldasymm,            &
				     abscoeff=abscoeff, &
                                     conc_drop_in=conc_drop,     &
                                    conc_ice_in=conc_ice,      &
                                     size_drop_in=size_drop,    &
                                     size_ice_in=size_ice)
       deallocate (conc_drop )
         deallocate (conc_ice  )
         deallocate (size_drop )
        deallocate (size_ice  )
	 endif
      endif  ! (ierr)
  
!       deallocate (ql )
!       deallocate (qi )
!       deallocate (cf )
  

!  deallocate (cosangsolar  )
!! END OF DROP

  !   else  ! (not running fms)

  !     call error_mesg ('strat_clouds_W',   &
  !       'strat_cloud only available in FMS', FATAL)

!     endif ! (running_fms)

  else if (Environment%running_standalone) then ! (running_standalone)
 
     if (Cldrad_control%do_pred_cld_microphys) then

!---------------------------------------------------------------------
!     obtain the cloud amounts and areas.
!---------------------------------------------------------------------
       allocate (ql (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
      allocate (qi (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
       allocate (cf (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )

       ql = cloud_water
       qi = cloud_ice
       cf = camtsw

!---------------------------------------------------------------------
!     allocate and initialize the cloud radiative property arrays.
!---------------------------------------------------------------------
       allocate (ktop   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
       allocate (kbtm   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!      allocate (cldamt (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
!      allocate (cldamt_out (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
       allocate (sp_hum (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
       allocate (conc_drop  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )  
       allocate (conc_ice  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
       allocate (size_drop  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) ) 
       allocate (size_ice  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )

       ktop=1
       kbtm=0
!      cldamt=0.
       cldamt_out=0.

!---------------------------------------------------------------------
!     define specific humidity to pass to start_cloud_mod.
!---------------------------------------------------------------------
! rh2o (  :  ,  :  ,:) = q(:,:,:)/(1.0-q(:, :, :)/rh2oair)
       sp_hum(:,:,:) = rh2o(:,:,:)/(1.0+rh2o(:,:,:)/rh2oair)



!---------------------------------------------------------------------
!     obtain the cloud radiative properties.
!---------------------------------------------------------------------

!   call cloud_summary (iabs(ISRAD), jabs(JSRAD), &
    call cloud_summary (is, js,                   &
!                       land, ql, qi, cf, rh2o,   &
                        land, ql, qi, cf, sp_hum, &
                        0.1*press(:,:,KSRAD:KERAD),&
                        0.1*pflux, temp(:,:,KSRAD:KERAD),  &
!                       cosangsolar(:,:,1),&
                        cosangsolar       ,&
                        temp(:,:,KERAD+1),        &
!                       ncldsw, ktop, kbtm, cldamt,         &
                        ncldsw, ktop, kbtm, cldamt=cldamt_out,   &
                      Time=Time_next,                &
                        conc_drop=conc_drop, conc_ice=conc_ice, &
                        size_drop=size_drop, size_ice=size_ice)

!
!--------------------------------------------------------------------
!  if microphysics is being used for either sw or lw calculation, call
!  microphys_rad_driver to obtain cloud radiative properties.
!--------------------------------------------------------------------
     if (do_esfsw                  ) then

       call microphys_rad_driver (is, ie, js, je, deltaz,  &
!          press, temp,            camtsw, emmxolw, &
!                                 emrndlw, Cld_diagnostics, &
           press, temp,                             &
                                           Cld_diagnostics, &
                                  cldext=cldext, cldsct=cldsct,    &
                                  cldasymm=cldasymm,            &
				  abscoeff=abscoeff, &
                                conc_drop_in=conc_drop,     &
                                  conc_ice_in=conc_ice,      &
                               size_drop_in=size_drop,    &
                                size_ice_in=size_ice)
    else if (do_lwcldemiss) then
        call microphys_rad_driver (is, ie, js, je, deltaz,  &
!          press, temp,     camtsw, emmxolw, &
!                                   emrndlw, Cld_diagnostics)
                           press, temp,                      &
                                             Cld_diagnostics, &
                                abscoeff=abscoeff, &
                                conc_drop_in=conc_drop,     &
                                  conc_ice_in=conc_ice,      &
                               size_drop_in=size_drop,    &
                                size_ice_in=size_ice)
    endif
  
    deallocate (ktop)
    deallocate (kbtm)
!         deallocate (cldamt)
!         deallocate (cldamt_out)
          deallocate (sp_hum)
       deallocate (conc_drop )
         deallocate (conc_ice  )
         deallocate (size_drop )
        deallocate (size_ice  )

        deallocate (ql )
        deallocate (qi )
        deallocate (cf )


  else !  (standalone, not with pred_cld_microphys)

       call error_mesg ('strat_clouds_W',   &
    ' standalone strat_cloud only available with pred microphys', FATAL)
 
  
  endif  ! (standalone, microphys)

 endif  ! (gcm or standalone)




end subroutine strat_clouds_calc


!####################################################################


	       end module strat_clouds_W_mod



