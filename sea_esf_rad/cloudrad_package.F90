
                 module cloudrad_package_mod

use time_manager_mod,        only: time_type
use diag_manager_mod,        only: register_diag_field, send_data
use rad_output_file_mod,     only: hold_cloud, put_cloudpackage_type
use donner_ice_mod,          only: inquire_donner_ice,   &
                                   sw_albedo_zen_angle
use standalone_clouds_mod,   only: standalone_clouds_init, &
			           standalone_clouds_driver
use no_clouds_mod,           only: no_clouds_init, no_clouds_calc
use diag_clouds_W_mod,       only: diag_clouds_init, diag_clouds_calc
use obs_clouds_W_mod,        only: obs_clouds_init, obs_clouds_calc
use zonal_clouds_W_mod,      only: zonal_clouds_init, zonal_clouds_calc
use strat_clouds_W_mod,      only: strat_clouds_init, strat_clouds_calc
use rh_based_clouds_mod,     only: rh_clouds_init, rh_clouds_calc
use mgrp_prscr_clds_mod,     only: mgrp_prscr_init, mgrp_prscr_calc
use utilities_mod,           only: open_file, file_exist,   &
                                   check_nml_error, error_mesg,   &
                                   print_version_number, FATAL, NOTE, &
				   WARNING, close_file, get_my_pe, &
				   read_data, write_data
use rad_step_setup_mod,      only: jabs, iabs, pflux, IMINP, IMAXP,   &
				   JMINP, JMAXP, ISRAD, IERAD, JSRAD,  &
				   JERAD, KSRAD, KERAD
use rad_utilities_mod,       only: longwave_control_type, Lw_control, &
                                   Environment, environment_type, &
                                   shortwave_control_type, Sw_control,&
				   radiation_control_type, Rad_control,&
				   cloudrad_control_type, Cldrad_control
use longwave_setup_mod,      only: longwave_parameter_type,  &
				   Lw_parameters
use radiation_diag_mod,      only: radiag_from_clouds_lh,   &
                                   radiag_from_clouds_esf
use std_pressures_mod,       only: get_std_pressures
use astronomy_package_mod,   only: get_astronomy_for_clouds,  &
			           get_astronomy_for_clouds_init
use esfsw_parameters_mod,    only: nbands
use constants_new_mod,       only: pstd
!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!	           cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!     character(len=5), parameter  ::  version_number = 'v0.09'
      character(len=128)  :: version =  '$Id: cloudrad_package.F90,v 1.2 2001/07/05 17:28:32 fms Exp $'
      character(len=128)  :: tag     =  '$Name: galway $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
	  cloudrad_package_init, cldmarch,  clouddrvr, &
	  cloudrad_package_alloc, cloudrad_package_dealloc, &
	  get_clouds_for_lhsw, get_clouds_for_esfsw, &
          get_ncldsw, get_clouds_for_lwrad,  &
	  get_swcldprops_from_cloudrad


private          &
	  default_clouds, cloudrad_netcdf, diag_field_init, &
	  compute_isccp_clds


!---------------------------------------------------------------------
!-------- namelist  ---------


character(len=12)            :: microphys_form = '            '
character(len=12)            :: cloud_type_form = '            '

!    logical variables derived from namelist input
logical                      :: do_pred_cld_microphys = .false.
logical                      :: do_presc_cld_microphys   = .false.
logical                      :: do_no_cld_microphys   = .false.

logical                      :: do_rh_clouds=.false.
logical                      :: do_strat_clouds=.false.
logical                      :: do_zonal_clouds=.false.
logical                      :: do_mgroup_prescribed = .false.
logical                      :: do_obs_clouds=.false.
logical                      :: do_no_clouds=.false.
logical                      :: do_diag_clouds=.false.
logical                      :: do_specified_clouds=.false.

integer                      :: n_prsc_clds=3




namelist /cloudrad_package_nml /     &
                               microphys_form, &
			       cloud_type_form, &
                               n_prsc_clds


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------


!-------------------------------------------------------------------
!    cloud property variables. these will be exported to other
!    radiation modules using appropriate output subroutines
!
!-------------------------------------------------------------------
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
!     cldext  =  the parameterization band values of the cloud      
!                extinction coefficient in kilometer**(-1)          
!     cldsct  =  the parameterization band values of the cloud      
!                scattering coefficient in kilometer**(-1)          
!     cldasymm=  the parameterization band values of the asymmetry  
!                factor                                             
!   the following arrays may be needed for shortwave cloud properties
!   in legacy radiation codes :
!
!     cirabsw =  absorptivity of clouds in the infrared frequency band.
!                may be zenith angle dependent.
!     cirrfsw =  reflectivity of clouds in the infrared frequency band.
!                may be zenith angle dependent.
!    cvisrfsw =  reflectivity of clouds in the visible frequency band.
!                may be zenith angle dependent.
!--------------------------------------------------------------------

real, dimension(:,:,:,:), allocatable ::  cldext, cldasymm, &
                                          cldsct,  emmxolw, emrndlw, &
                                          cirabsw, cirrfsw, cvisrfsw
real,    dimension(:,:,:),allocatable ::  camtsw, cmxolw, crndlw
integer, dimension(:,:),  allocatable ::  ncldsw, nmxolw, nrndlw

!------------------------------------------------------------------
!    NLWCLDB is the actual number of frequency bands for which lw
!    emissitivies are defined. 
!------------------------------------------------------------------
integer            :: NLWCLDB 
 
integer            :: kmin, kmax
integer            :: nsolwg
character(len=10)  :: swform

!-------------------- diagnostics fields ------------------------------

integer :: id_tot_cld_amt, id_cld_amt,   id_em_cld_lw, id_em_cld_10u, & 
           id_high_cld_amt, id_mid_cld_amt, id_low_cld_amt,  &
           id_ext_cld_uv,   id_sct_cld_uv,  id_asymm_cld_uv, &
           id_ext_cld_vis,  id_sct_cld_vis, id_asymm_cld_vis, &
           id_ext_cld_nir,  id_sct_cld_nir, id_asymm_cld_nir, &
           id_alb_uv_cld, id_alb_nir_cld, id_abs_uv_cld, id_abs_nir_cld


character(len=8), parameter :: mod_name = 'cloudrad'

real :: missing_value = -999.



!----------------------------------------------------------------------
!----------------------------------------------------------------------



             contains 





subroutine cloudrad_package_init (kmin_in, kmax_in, &
				  th, axes, Time, qlyr, lonb, latb)

!------------------------------------------------------------------
integer,               intent(in)             :: kmin_in, kmax_in
real,    dimension(:), intent(in)             :: th   
real,    dimension(:), intent(in), optional   :: qlyr, lonb, latb
integer, dimension(4), intent(in), optional   :: axes
type(time_type),       intent(in), optional   :: Time

!------------------------------------------------------------------
! optional arguments in FMS : axes, Time, qlyr, lonb, latb
! optional arguments in skyhi: Time, qlyr
! optional arguments in standalone : none
!--------------------------------------------------------------------

!-----------------------------------------------------------------
!    local variables
!-----------------------------------------------------------------

    integer                           :: unit, ierr, io
    integer                           :: ix, iy
    real, dimension(:), allocatable   :: qlevel, pd
    real                              :: psbar

!---------------------------------------------------------------------
!-----  read namelist  ------
  
    if (file_exist('input.nml')) then
      unit =  open_file ('input.nml', action='read')
      ierr=1; do while (ierr /= 0)
      read (unit, nml=cloudrad_package_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'cloudrad_package_nml')
      enddo
10    call close_file (unit)
    endif

    unit = open_file ('logfile.out', action='append')
!   call print_version_number (unit, 'cloudrad_package', version_number)
    if (get_my_pe() == 0) then 
      write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
      write (unit,nml=cloudrad_package_nml)
    endif
    call close_file (unit)

!-------------------------------------------------------------------
!  check for consistency between the specified radiation options and
!  the specified microphysics 
!----------------------------------------------------------------------

    swform = Sw_control%sw_form

    if (trim(swform) == 'esfsw99' .and.   &
        trim(microphys_form) == 'none') then 
      call error_mesg( 'cloudrad_package_init',  &
      ' must specify microphysics when using esfsw99 shortwave.', FATAL)
    endif

    if (trim(swform) == 'lhsw' .and.     &
	      .not. Lw_control%do_lwcldemiss .and.    &
	trim(microphys_form) /= 'none') then 
      call error_mesg( 'cloudrad_package_init',  &
       ' not using microphysics -- set microphys_form to none.', FATAL)
    endif

    if (Lw_control%do_lwcldemiss .and.   &
		 trim(microphys_form) == 'none') then 
      call error_mesg( 'cloudrad_package_init',  &
       ' lwcldemiss should not be used without microphysics.', FATAL)
    endif

!-------------------------------------------------------------------
!  define type of model cloud formulation being used
!----------------------------------------------------------------------

    if (Environment%running_gcm) then
!-------------------------------------------------------------------
!  cloud fractions, heights are predicted by the model based on rel hum
!-------------------------------------------------------------------

      if (trim(cloud_type_form)  == 'rh')   then
        do_rh_clouds = .true.

!-------------------------------------------------------------------
!  cloud fractions, heights are predicted by the model based on klein  
!  scheme
!-------------------------------------------------------------------

      else if (trim(cloud_type_form) == 'strat')  then
        if (Environment%running_fms) then
          do_strat_clouds = .true.
        else if (Environment%running_skyhi) then
          call error_mesg( 'cloudrad_package_init',  &
                ' strat clouds not available in SKYHI.', FATAL)
        endif

!-------------------------------------------------------------------
!  cloud fractions, heights are prescribed in the model using fms form
!-------------------------------------------------------------------

      else if (trim(cloud_type_form) == 'zonal')  then
        if (Environment%running_fms) then
          do_zonal_clouds = .true.
        else if (Environment%running_skyhi) then
          call error_mesg( 'cloudrad_package_init',  &
                      ' zonal clouds not available in SKYHI.', FATAL)
        endif

!-------------------------------------------------------------------
!  cloud fractions, heights are based on observations used by FMS
!-------------------------------------------------------------------

      else if (trim(cloud_type_form) == 'obs')  then
        if (Environment%running_fms) then
          do_obs_clouds = .true.
        else if (Environment%running_skyhi) then
          call error_mesg( 'cloudrad_package_init',  &
                         ' obs clouds not available in SKYHI.', FATAL)
        endif

!-------------------------------------------------------------------
!  cloud fractions, heights are prescribed in the model using skyhi form
!-------------------------------------------------------------------

      else if (trim(cloud_type_form)  == 'prescribed')  then
        do_mgroup_prescribed = .true.

!-------------------------------------------------------------------
!  model is run with Gordon diagnostic clouds
!-------------------------------------------------------------------

      else if (trim(cloud_type_form)  == 'diag')  then
        do_diag_clouds = .true.

!-------------------------------------------------------------------
!  model is run without clouds
!-------------------------------------------------------------------

      else if (trim(cloud_type_form) == 'none')  then
        do_no_clouds = .true.

!-------------------------------------------------------------------
!  error condition
!-------------------------------------------------------------------

      else
        call error_mesg( 'cloudrad_package_init',  &
                       ' invalid cloud_type_form specified.', FATAL)
      endif

!-------------------------------------------------------------------
!  define type of cloud microphysics being used. verify that the option
!  is consistent with the type of cloud formulation being requested.
!----------------------------------------------------------------------
      if (trim(microphys_form) == 'predicted') then 

!----------------------------------------------------------------------
!  cloud species concentrations, sizes, etc come from a microphysics
!  model
!----------------------------------------------------------------------

        if (do_rh_clouds) then
          call error_mesg( 'cloudrad_package_init',  &
        ' predicted microphys not yet available with rh clouds.', FATAL)
        else if (do_strat_clouds) then
	  do_pred_cld_microphys = .true.
        else if (do_mgroup_prescribed) then
          call error_mesg( 'cloudrad_package_init',  &
  ' predicted microphys not available with mgroup prescribed clouds', &
    FATAL)
        else if (do_zonal_clouds) then
          call error_mesg( 'cloudrad_package_init',  &
            ' predicted microphys not available with zonal clouds',  &
	   FATAL)
        else if (do_obs_clouds) then
          call error_mesg( 'cloudrad_package_init',  &
           ' predicted microphys not available with observed clouds', &
	   FATAL)
	else if (do_diag_clouds) then
          call error_mesg( 'cloudrad_package_init',  &
        ' predicted microphys not available with gordon diag clouds', &
	   FATAL)
        else if (do_no_clouds) then
        endif

      else if (trim(microphys_form) == 'prescribed') then
!----------------------------------------------------------------------
!  cloud species concentrations, sizes, etc are prescribed (for
!  low, middle, high clouds)
!----------------------------------------------------------------------

        if (do_rh_clouds) then
          do_presc_cld_microphys = .true.
        else if (do_strat_clouds) then
          call error_mesg( 'cloudrad_package_init',  &
        ' use predicted microphys with the strat cloud module',  &
	  FATAL)
        else if (do_mgroup_prescribed) then
          do_presc_cld_microphys = .true.
        else if (do_zonal_clouds) then
          call error_mesg( 'cloudrad_package_init',  &
           ' prescribed microphys not available with zonal clouds', &
	   				                    FATAL)
        else if (do_obs_clouds) then
          call error_mesg( 'cloudrad_package_init',  &
          ' prescribed microphys not available with observed clouds', &
							    FATAL)
        else if (do_diag_clouds) then
          call error_mesg( 'cloudrad_package_init',  &
       ' prescribed microphys not available with gordon diag clouds', &
							    FATAL)
        else if (do_no_clouds) then
        endif

      else if (trim(microphys_form) == 'none') then
!----------------------------------------------------------------------
!  cloud reflectivities, absorptivities, emissivities are prescribed
!  for low, middle, high clouds (legacy code)
!----------------------------------------------------------------------

        if (do_rh_clouds) then
          do_no_cld_microphys = .true.
        else if (do_strat_clouds) then
          do_no_cld_microphys = .true.
        else if (do_mgroup_prescribed) then
          do_no_cld_microphys = .true.
        else if (do_zonal_clouds) then
          do_no_cld_microphys = .true.
        else if (do_obs_clouds) then
          do_no_cld_microphys = .true.
        else if (do_diag_clouds) then
          do_no_cld_microphys = .true.
        else if (do_no_clouds) then
        endif

      else
!----------------------------------------------------------------------
!  error condition
!----------------------------------------------------------------------
        call error_mesg( 'cloudrad_package_init',  &
                    ' microphys_form is not an acceptable value.', &
		                                                FATAL)
      endif
    endif

    if (Environment%running_standalone) then
!-------------------------------------------------------------------
!  define type of standalone cloud formulation being used
!----------------------------------------------------------------------

       if (trim(cloud_type_form) == 'strat')  then
!-------------------------------------------------------------------
!  cloud properties are predicted by the model based on klein
!  scheme
!-------------------------------------------------------------------

         do_strat_clouds = .true.


       else if (trim(cloud_type_form)  == 'specified')  then
!-------------------------------------------------------------------
!  model is run with specified clouds and cloud properties
!-------------------------------------------------------------------
        do_specified_clouds = .true.

       else
         call error_mesg( 'cloudrad_package_init',  &
               'invalid standalone cloud_type_form specified.', FATAL)
       endif
!-------------------------------------------------------------------
!  define type of cloud microphysics being used. verify that the option
!  is consistent with the type of cloud formulation being requested.
!----------------------------------------------------------------------
       if (trim(microphys_form) == 'predicted') then
  
!----------------------------------------------------------------------
!  cloud species concentrations, sizes, etc come from a microphysics
!  model
!----------------------------------------------------------------------

         if (do_strat_clouds) then
           do_pred_cld_microphys = .true.
         else if (do_specified_clouds) then
           call error_mesg( 'cloudrad_package_init',  &
              'if microphys_form is predicted, must activate &
	          &strat clouds', FATAL)
         endif
 
       else if (trim(microphys_form) == 'prescribed') then
!----------------------------------------------------------------------
!  cloud species concentrations, sizes, etc are prescribed (for
!  low, middle, high clouds)
!----------------------------------------------------------------------

         if (do_strat_clouds) then
           call error_mesg( 'cloudrad_package_init',  &
            ' prescribed microphys not available with strat clouds',  &
             FATAL)
         endif

       else
!----------------------------------------------------------------------
!  error condition
!----------------------------------------------------------------------
         call error_mesg( 'cloudrad_package_init',  &
                 ' microphys_form is not an acceptable value.', &
                    FATAL)
       endif
    endif

!---------------------------------------------------------------------
!  define microphysics type control variable components
!---------------------------------------------------------------------

    Cldrad_control%do_pred_cld_microphys  = do_pred_cld_microphys
    Cldrad_control%do_presc_cld_microphys = do_presc_cld_microphys    
    Cldrad_control%do_no_cld_microphys    = do_no_cld_microphys

!--------------------------------------------------------------------
!  define model level top and bottom indices. retrieve module variables 
!  that come from other modules
!--------------------------------------------------------------------

    kmin = kmin_in
    kmax = kmax_in
    call get_astronomy_for_clouds_init (nsolwg)

!---------------------------------------------------------------------
!    save the sigma levels that have been input for later use within 
!    this module.
!---------------------------------------------------------------------
    allocate (qlevel (KSRAD:KERAD) )
    if (present(qlyr)) then
      qlevel(:) = qlyr(KSRAD:KERAD)
    else
      allocate ( pd(KMIN:KMAX+1) )
      call get_std_pressures (pd_out=pd)
      psbar = 1.0E-03*pstd  ! convert from cgs to mb
      qlevel(KSRAD:KERAD) = pd(KSRAD:KERAD)/psbar
      deallocate (pd)
    endif

!---------------------------------------------------------------------
!    define the number of cloud emissivity bands for use in this module.
!---------------------------------------------------------------------

    NLWCLDB = Lw_parameters%NLWCLDB

!-------------------------------------------------------------------
!  initialize the specific cloud scheme selected for this run
!-------------------------------------------------------------------

    if (Environment%running_gcm) then
      if (do_mgroup_prescribed) then
        call mgrp_prscr_init
      endif

      if (do_rh_clouds) then
        call rh_clouds_init (th, qlevel)
      endif

      if (do_no_clouds) then
        call no_clouds_init 
      endif

      if (Environment%running_fms) then
        if (do_zonal_clouds) then
          call zonal_clouds_init
        endif

        if (do_obs_clouds) then
          call obs_clouds_init (lonb, latb)
        endif

        if (do_strat_clouds) then
          call strat_clouds_init
        endif

	if (do_diag_clouds) then
          ix = size(lonb) - 1
	  iy = size(latb) - 1
	  call diag_clouds_init (ix, iy, kmax, th, ierr)
	endif
      endif
    else if (Environment%running_standalone) then
      if (do_strat_clouds) then
        call strat_clouds_init
      endif

      call standalone_clouds_init (                                 &
                        th, do_strat_clouds, do_specified_clouds)
    endif

!-------------------------------------------------------------------
!  initialize diagnostics of this module
!-------------------------------------------------------------------
    if (Environment%running_gcm .and. Environment%running_fms  &
				 .and. .not. do_no_clouds) then
      call diag_field_init (Time, axes)
    endif

!-------------------------------------------------------------------
!  pass cloud package type to rad_output_file_mod
!-------------------------------------------------------------------
    call put_cloudpackage_type (cloud_type_form)

!--------------------------------------------------------------------



end subroutine cloudrad_package_init



!#################################################################

subroutine cldmarch

!---------------------------------------------------------------------
!     cldmarch calls routines to handle temporal variation of cloud
!     variables (currently null).
!---------------------------------------------------------------------


end subroutine cldmarch



!####################################################################

subroutine clouddrvr (is, ie, js, je, Time_next, kbot, mask)

integer,                 intent(in)           ::  is, ie, js, je
type(time_type),         intent(in), optional ::  Time_next
integer, dimension(:,:), intent(in), optional ::  kbot
real, dimension(:,:,:),  intent(in), optional ::  mask
!-------------------------------------------------------------------
 
!--------------------------------------------------------------------
!     initialize the cloud property arrays
!-------------------------------------------------------------------
      call default_clouds

!--------------------------------------------------------------------
!     call desired cloud package to obtain needed cloud radiative 
!     property fields.
!-------------------------------------------------------------------
      if (Environment%running_gcm) then
        if (do_rh_clouds) then
          if (trim(swform) /= 'esfsw99') then
            call rh_clouds_calc (camtsw, cmxolw, crndlw, ncldsw,   &
		                 nmxolw, nrndlw, emmxolw, emrndlw, &
	   		         cirabsw=cirabsw, cvisrfsw=cvisrfsw, &
			         cirrfsw=cirrfsw) 
          else
            call rh_clouds_calc (camtsw, cmxolw, crndlw, ncldsw,   &
	  		         nmxolw,  nrndlw, emmxolw, emrndlw, &
			         cldext=cldext, cldsct=cldsct,       &
			         cldasymm=cldasymm)
          endif
        endif
        if (do_no_clouds) then
	  call no_clouds_calc     &
	                      (ncldsw, nrndlw, camtsw, crndlw, emrndlw)
        endif
        if (do_mgroup_prescribed) then
          if (trim(swform) /= 'esfsw99') then
	    call mgrp_prscr_calc    &
	                    (camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
			     nrndlw,                             &
			     emmxolw, emrndlw, cirabsw=cirabsw, &
			     cvisrfsw=cvisrfsw, cirrfsw=cirrfsw)
          else
	    call mgrp_prscr_calc    &
	                    (camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
			     nrndlw,                             &
			     emmxolw, emrndlw, cldext=cldext, &
			     cldsct=cldsct, cldasymm=cldasymm)
          endif
        endif
        if (Environment%running_fms) then
          if (do_strat_clouds) then
	  if (trim(swform) /= 'esfsw99') then
	    call strat_clouds_calc  &
	                    (camtsw, cmxolw, crndlw, ncldsw, nmxolw,  &
			     nrndlw, emmxolw, emrndlw,cirabsw=cirabsw,&
			     cvisrfsw=cvisrfsw, cirrfsw=cirrfsw, &
			     Time_next=Time_next)
	    else
              call strat_clouds_calc  &
                            (camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
                             nrndlw, emmxolw, emrndlw,           &
                          Time_next=Time_next,                  &
                           cldext=cldext, &
                          cldsct=cldsct, cldasymm=cldasymm)
             endif
          else if (do_zonal_clouds) then
	    call zonal_clouds_calc   &
	                    (camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
			     nrndlw, cirabsw, cvisrfsw, cirrfsw, &
			     emmxolw, emrndlw)
          else if (do_obs_clouds) then
	    call obs_clouds_calc  &
	                    (camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
			     nrndlw, cirabsw, cvisrfsw, cirrfsw, &
			     emmxolw, emrndlw, is, ie, js, je)
	  else if (do_diag_clouds) then
	    call diag_clouds_calc  &
	                    (camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
			     nrndlw, cirabsw, cvisrfsw, cirrfsw, &
			     emmxolw, emrndlw, is, js, Time_next)
          endif
        endif
      else if (Environment%running_standalone) then
        if (trim(swform) /= 'esfsw99') then
          call standalone_clouds_driver (   &
                        do_strat_clouds, &
	   		      camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
		  	       nrndlw, emmxolw, emrndlw,   &
 		       Time_next=Time_next,  &
			       cirabsw=cirabsw, cvisrfsw=cvisrfsw, &
			       cirrfsw=cirrfsw)
        else
          call standalone_clouds_driver  (   &
                        do_strat_clouds,  &
			       camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
		  	       nrndlw, emmxolw, emrndlw,   &
 		       Time_next=Time_next,  &
			       cldext=cldext, cldsct=cldsct,       &
			       cldasymm=cldasymm)
        endif
      endif

!-------------------------------------------------------------------
!    send desired variables to archive package module to hold until
!    output.
!---------------------------------------------------------------------

      call hold_cloud (cmxolw, crndlw)

!-------------------------------------------------------------------
!     generate netcdf file output fields
!-------------------------------------------------------------------
      if (Environment%running_gcm .and. Environment%running_fms .and. &
                            .not. do_no_clouds) then
        call cloudrad_netcdf (is, js, Time_next)
      endif  

!---------------------------------------------------------------------


end subroutine clouddrvr     



!###################################################################

subroutine cloudrad_package_alloc

    allocate (   camtsw    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD ))
    allocate (   cmxolw    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD ))
    allocate (   crndlw    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD ))

    allocate(emmxolw (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NLWCLDB))
    allocate(emrndlw (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NLWCLDB))

    allocate (   ncldsw    (ISRAD:IERAD, JSRAD:JERAD              ))
    allocate (   nmxolw    (ISRAD:IERAD, JSRAD:JERAD              ))
    allocate(    nrndlw    (ISRAD:IERAD, JSRAD:JERAD              ))

    if (trim(swform) /= 'esfsw99' ) then

      allocate(cirabsw (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, nsolwg ))
      allocate(cirrfsw (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, nsolwg))
      allocate(cvisrfsw(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, nsolwg ))

    else

      allocate (cldext(IMINP:IMAXP,JMINP:JMAXP,KMIN:KMAX,nbands) )
      allocate (cldsct(IMINP:IMAXP,JMINP:JMAXP,KMIN:KMAX,nbands) )
      allocate (cldasymm(IMINP:IMAXP,JMINP:JMAXP,KMIN:KMAX,nbands) )

    endif

end subroutine cloudrad_package_alloc



!####################################################################

subroutine cloudrad_package_dealloc
       
    deallocate (camtsw    )
    deallocate (cmxolw    )
    deallocate (crndlw    )
    deallocate (emmxolw   )
    deallocate (emrndlw   )
    deallocate (ncldsw    )
    deallocate (nmxolw    )
    deallocate (nrndlw    )

    if (trim(swform) /= 'esfsw99' ) then
      deallocate (cirabsw   )
      deallocate (cirrfsw   )
      deallocate (cvisrfsw  )
    else
      deallocate (cldasymm) 
      deallocate (cldsct  ) 
      deallocate (cldext  ) 
    endif



end subroutine cloudrad_package_dealloc



!####################################################################

subroutine get_clouds_for_lhsw (cirabswkc, cirrfswkc,  &
				cvisrfswkc, ktopswkc, kbtmswkc,  &
				camtswkc)
      
real, dimension(:,:,:),   intent(out) :: camtswkc
real, dimension(:,:,:,:), intent(out) :: cirabswkc, cirrfswkc,  &
					 cvisrfswkc
integer, dimension(:,:,:),intent(out) :: ktopswkc, kbtmswkc
!-------------------------------------------------------------------

      logical :: do_donner_ice
      integer :: i,j,k, index_kbot, index_ktop, ngp

!---------------------------------------------------------------------
!     obtain arrays for shortwave cloud properties in "cloud space".
!     "bottom up" counting is used in conformity with the lhsw
!     radiative algorithm. (ie, index = 1 = lowest cloud (if any), etc.)
!---------------------------------------------------------------------
      
      camtswkc(:,:,:) = 0.0E+00
      cirabswkc(:,:,:,:) = 0.0E+00
      cirrfswkc(:,:,:,:) = 0.0E+00
      cvisrfswkc(:,:,:,:) = 0.0E+00
      ktopswkc(:,:,:) = KSRAD
      kbtmswkc(:,:,:) = KSRAD

      do j=JSRAD,JERAD
        do i=ISRAD,IERAD
	  index_kbot = 0
	  index_ktop = 0
!---------------------------------------------------------------------
!     in the KERAD'th layer, the presence of cloud (camtsw) implies a
!     shortwave cloud bottom at level (KERAD+1) but nothing about
!     cloud tops.
!---------------------------------------------------------------------
          if (camtsw(i,j,KERAD) .GT. 0.0E+00) then
            index_kbot = index_kbot + 1
	    kbtmswkc(i,j,index_kbot) = KERAD+1
	    camtswkc(i,j,index_kbot) = camtsw(i,j,KERAD)
	    do ngp=1,nsolwg  
	      cirabswkc(i,j,index_kbot,ngp) = cirabsw(i,j,KERAD,ngp)
              cirrfswkc(i,j,index_kbot,ngp) = cirrfsw(i,j,KERAD,ngp)
	      cvisrfswkc(i,j,index_kbot,ngp) = cvisrfsw(i,j,KERAD,ngp)
	    end do
          endif
!---------------------------------------------------------------------
!   in other layers, cloud bottoms and tops are determined according
!   to changes in shortwave cloud (and special case).
!---------------------------------------------------------------------
          do k=KERAD-1,KSRAD,-1
!---------------------------------------------------------------------
!     cloud bottoms.
!--------------------------------------------------------------------
	    if (camtsw(i,j,k) .NE. camtsw(i,j,k+1) .AND.    &
      		camtsw(i,j,k) .GT. 0.0E+00              .OR.  &
!---------------------------------------------------------------------
!     special case where shortwave cloud amounts for adjacent
!     layers are equal, but a randomly overlapped cloud exists
!     in at least one of these layers
!---------------------------------------------------------------------
                camtsw(i,j,k) .EQ. camtsw(i,j,k+1) .AND.   &
               (crndlw(i,j,k) .NE. 0.0E+00 .OR.      &
                crndlw(i,j,k+1) .NE. 0.0E+00)                ) then
              index_kbot = index_kbot + 1
	      kbtmswkc(i,j,index_kbot) = k+1
	      camtswkc(i,j,index_kbot) = camtsw(i,j,k)
	      do ngp=1,nsolwg  
	      cirabswkc(i,j,index_kbot,ngp) = cirabsw(i,j,k,ngp)
	      cirrfswkc(i,j,index_kbot,ngp) = cirrfsw(i,j,k,ngp)
	      cvisrfswkc(i,j,index_kbot,ngp) = cvisrfsw(i,j,k,ngp)
	      end do
            endif
!---------------------------------------------------------------------
!     cloud tops.
!---------------------------------------------------------------------
            if (camtsw(i,j,k) .NE. camtsw(i,j,k+1) .AND.      &
                camtsw(i,j,k+1) .GT. 0.0E+00            .OR.  &
!---------------------------------------------------------------------
!     special case where shortwave cloud amounts for adjacent
!     layers are equal, but a randomly overlapped cloud exists
!     in at least one of these layers
!---------------------------------------------------------------------
                camtsw(i,j,k) .EQ. camtsw(i,j,k+1) .AND.    &
               (crndlw(i,j,k) .NE. 0.0E+00 .OR.        &
                crndlw(i,j,k+1) .NE. 0.0E+00)                ) then
	      index_ktop = index_ktop + 1
	      ktopswkc(i,j,index_ktop) = k+1
            endif
          enddo
        enddo
      enddo

!-------------------------------------------------------------------
!     if donner iceclouds are active, send some info to that module.
!---------------------------------------------------------------------

      if (Environment%running_gcm) then
        if (Environment%using_sky_periphs) then
          call inquire_donner_ice (do_donner_ice)
          if (do_donner_ice) then
            call sw_albedo_zen_angle (ncldsw, cvisrfswkc, cirrfswkc,  &
                                     nsolwg)
          endif
        endif
      endif 
!-------------------------------------------------------------------
!     send needed lhsw cloud information to radiation diagnostics mod.
!---------------------------------------------------------------------

      do j=JSRAD,JERAD
        if (Rad_control%do_raddg(jabs (j))) then
          call radiag_from_clouds_lh(j, camtswkc, cirabswkc, cirrfswkc,&
                                     cmxolw, crndlw, cvisrfswkc,   &
				     emmxolw, emrndlw, ncldsw, nmxolw, &
				     nrndlw, ktopswkc, kbtmswkc) 
        endif  
      end do
!---------------------------------------------------------------------



end subroutine get_clouds_for_lhsw



!####################################################################

subroutine get_clouds_for_esfsw (camtsw_out)
      
!--------------------------------------------------------------------
real, dimension(:,:,:), intent(out)   :: camtsw_out
!-------------------------------------------------------------------

      integer           ::  j
!-------------------------------------------------------------------

      camtsw_out(:,:,:) = camtsw(:,:,:)

      do j=JSRAD,JERAD
        if (Rad_control%do_raddg(jabs (j))) then
          call radiag_from_clouds_esf(j, camtsw_out, cmxolw, crndlw, &
                                      emmxolw, emrndlw, ncldsw, nmxolw,&
				      nrndlw)
        endif  
      end do


end subroutine get_clouds_for_esfsw


!#####################################################################

subroutine get_ncldsw (ncldsw_out)

integer, dimension(:,:), intent(out) :: ncldsw_out


      ncldsw_out(:,:) = ncldsw(:,:)


end subroutine get_ncldsw


!###################################################################

subroutine get_clouds_for_lwrad (cmxolw_out, crndlw_out, emmxolw_out, &
				 emrndlw_out, nmxolw_out, nrndlw_out)

real, dimension(:,:,:),   intent(out) :: cmxolw_out, crndlw_out
real, dimension(:,:,:,:), intent(out) :: emmxolw_out, emrndlw_out
integer, dimension(:,:),  intent(out) :: nmxolw_out, nrndlw_out


      cmxolw_out(:,:,:)    = cmxolw(:,:,:)
      crndlw_out(:,:,:)    = crndlw(:,:,:)
      emmxolw_out(:,:,:,:) = emmxolw(:,:,:,:)
      emrndlw_out(:,:,:,:) = emrndlw(:,:,:,:)
      nmxolw_out(:,:)      = nmxolw(:,:)
      nrndlw_out(:,:)      = nrndlw(:,:)


end subroutine get_clouds_for_lwrad



!#################################################################### 

subroutine get_swcldprops_from_cloudrad (n, cldext_out, cldsct_out, &
					 cldasymm_out)
 
!---------------------------------------------------------------------
integer,                intent(in)    ::  n
real, dimension(:,:,:), intent(out)   ::  cldext_out, cldsct_out,   &
					  cldasymm_out
!---------------------------------------------------------------------

      integer                     :: i,j,k
!---------------------------------------------------------------------

      do k=KSRAD,KERAD
        do j=JSRAD,JERAD
          do i=ISRAD, IERAD
            cldext_out   (i,j,k) =   cldext   (i,j,k,n)
            cldsct_out   (i,j,k) =   cldsct   (i,j,k,n)
            cldasymm_out (i,j,k) =   cldasymm (i,j,k,n)
          enddo
        enddo
      enddo

  
end subroutine get_swcldprops_from_cloudrad 

!####################################################################

subroutine default_clouds

      
      cmxolw(:,:,:) = 0.0E+00
      crndlw(:,:,:) = 0.0E+00
      camtsw(:,:,:) = 0.0E+00
      emmxolw(:,:,:,:) = 1.0E+00
      emrndlw(:,:,:,:) = 1.0E+00
      nmxolw (:,:) = 0
      nrndlw (:,:) = 0
      ncldsw (:,:) = 0

      if (allocated (cirrfsw) ) then
        cirrfsw (:,:,:,:) = 0.0E+00
        cvisrfsw(:,:,:,:) = 0.0E+00
        cirabsw (:,:,:,:) = 0.0E+00
      else if (allocated (cldext)  ) then
        cldsct  (:,:,:,:) = 0.0E+00
        cldext  (:,:,:,:) = 0.0E+00
        cldasymm(:,:,:,:) = 0.0E+00
      endif


end subroutine default_clouds



!####################################################################

subroutine cloudrad_netcdf(is, js, Time_diag, mask)

integer,                 intent(in)           ::  is, js
real, dimension(:,:,:),  intent(in), optional ::  mask
type(time_type),         intent(in)           ::  Time_diag
!-------------------------------------------------------------------

 !------------------------------------------------------------------
 !  local variables
 !------------------------------------------------------------------
      real, dimension(:,:,:), allocatable    :: cloud, hml_ca
      real,dimension(:,:), allocatable       :: tca
      integer                                :: k, j, i
      logical                                :: used
!---------------------------------------------------------------------

      allocate (cloud(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )

!------- TOTAL CLOUD AMOUNT  -----------------
      if ( id_tot_cld_amt > 0 ) then
          allocate (tca (ISRAD:IERAD, JSRAD:JERAD) )
          tca = 1.0
	  do k=KSRAD,KERAD
	    tca(:,:) = tca(:,:)*(1.0-camtsw(:,:,k))
          end do
	  tca = 1.-tca
         tca = tca*100.
         used = send_data ( id_tot_cld_amt, tca, Time_diag, is, js)
         deallocate (tca)
      endif

!---- high,mid,low cloud diagnostics ----
      if ( id_high_cld_amt > 0 .or. id_mid_cld_amt > 0 .or. &
            id_low_cld_amt > 0 ) then
        allocate (hml_ca(ISRAD:IERAD, JSRAD:JERAD,3))
        call compute_isccp_clds (pflux, camtsw, hml_ca)
        if ( id_high_cld_amt > 0 ) used = send_data &
            ( id_high_cld_amt, hml_ca(:,:,1), Time_diag, is, js )
           if ( id_mid_cld_amt > 0 ) used = send_data &
               ( id_mid_cld_amt, hml_ca(:,:,2), Time_diag, is, js )
           if ( id_low_cld_amt > 0 ) used = send_data &
                ( id_low_cld_amt, hml_ca(:,:,3), Time_diag, is, js )
          deallocate (hml_ca)
      endif


!------- cloud amount (only do once ?) -------------------------
      if ( id_cld_amt > 0 ) then
    used = send_data ( id_cld_amt, camtsw, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- cloud emissivity ---------------------------------------

          if ( id_em_cld_10u > 0  ) then
          if (Lw_control%do_lwcldemiss) then
!    output cloud emissivity is  weighted average of the random and max
!    overlap emissivities, over the 990-1070 cm-1 band (band 5 of 7)
            cloud(:,:,:) = (crndlw(:,:,:)*emrndlw(:,:,:,5) +    &
                           cmxolw(:,:,:)*emmxolw(:,:,:,5))/   &
                            (crndlw(:,:,:) + cmxolw(:,:,:) + 1.0E-10)
      used = send_data ( id_em_cld_10u, cloud, Time_diag, is, js, 1, &
                rmask=mask )
     endif
    endif

  if (id_em_cld_lw > 0   ) then
     if (Lw_control%do_lwcldemiss) then
     else
!    output cloud emissivity is weighted average of the random and max
!    overlap emissivities, over 1 band (0-2200 cm-1)
       cloud(:,:,:) = (crndlw(:,:,:)*emrndlw(:,:,:,1) +    &
                      cmxolw(:,:,:)*emmxolw(:,:,:,1))/   &
      (crndlw(:,:,:) + cmxolw(:,:,:) + 1.0E-10)
   used = send_data ( id_em_cld_lw, cloud, Time_diag, is, js, 1, &
               rmask=mask )
  endif
endif

!     if ( id_em_cld > 0  ) then
!             cloud(:,:,:) = (crndlw(:,:,:)*emrndlw(:,:,:,1) +    &
!		     cmxolw(:,:,:)*emmxolw(:,:,:,1))/   &
!		     (crndlw(:,:,:) + cmxolw(:,:,:) + 1.0E-10)
!        used = send_data ( id_em_cld, cloud, Time_diag, is, js, 1, &
!                           rmask=mask )
!     endif

      if (trim(swform) == 'lhsw' ) then  ! diagnostics for lacis-hansen

!------- ultra-violet reflected by cloud -----------------------------
        if ( id_alb_uv_cld > 0  ) then
              cloud(:,:,:) = cvisrfsw(:,:,:,1)
          used = send_data ( id_alb_uv_cld, cloud, Time_diag, is, js,  &
			    1, rmask=mask )
        endif

!------- infra-red reflected by cloud -----------------------------
        if ( id_alb_nir_cld > 0  ) then
              cloud(:,:,:) = cirrfsw(:,:,:,1)
          used = send_data ( id_alb_nir_cld, cloud, Time_diag, is, js, &
			    1, rmask=mask )
        endif

!------- ultra-violet absorbed by cloud (not implemented)------------
!       if ( id_abs_uv_cld > 0 ) then
!         cloud=0.0
!         do j=1,jd; do i=1,id
!           do n=1,nclds(i,j)
!             cloud(i,j,ktop(i,j,n+1):kbtm(i,j,n+1)) = cuvab(i,j,n+1)
!           enddo
!         enddo; enddo
!         used = send_data ( id_abs_uv_cld, cloud, Time_diag, is, js, 1, &
!                           rmask=mask )
!       endif

!------- infra-red absorbed by cloud -----------------------------
        if ( id_abs_nir_cld > 0  ) then
              cloud(:,:,:) = cirabsw(:,:,:,1)
          used = send_data ( id_abs_nir_cld, cloud, Time_diag, is, js, &
			   1, rmask=mask )
        endif
  else if (trim(swform) == 'esfsw99') then ! diagnostics for esf sw code

     if ( id_ext_cld_uv > 0  ) then
       cloud(:,:,:) = cldext(:,:,:,22)
     used = send_data ( id_ext_cld_uv, cloud, Time_diag, is, js, &
                          1, rmask=mask )
      endif
  
   if ( id_sct_cld_uv > 0  ) then
       cloud(:,:,:) = cldsct(:,:,:,22)
       used = send_data ( id_sct_cld_uv, cloud, Time_diag, is, js, &
                          1, rmask=mask )
    endif

    if ( id_asymm_cld_uv > 0  ) then
         cloud(:,:,:) = 100.0*cldasymm(:,:,:,22)
         used = send_data ( id_asymm_cld_uv, cloud, Time_diag, is, js, &
                          1, rmask=mask )
     endif

    if ( id_ext_cld_vis > 0  ) then
      cloud(:,:,:) = cldext(:,:,:,12)
      used = send_data ( id_ext_cld_vis, cloud, Time_diag, is, js, &
                             1, rmask=mask )
    endif

    if ( id_sct_cld_vis > 0  ) then
      cloud(:,:,:) = cldsct(:,:,:,12)
      used = send_data ( id_sct_cld_vis, cloud, Time_diag, is, js, &
                           1, rmask=mask )
    endif

     if ( id_asymm_cld_vis > 0  ) then
        cloud(:,:,:) = 100.0*cldasymm(:,:,:,12)
        used = send_data ( id_asymm_cld_vis, cloud, Time_diag, is, js, &
                            1, rmask=mask )
     endif

     if ( id_ext_cld_nir > 0  ) then
      cloud(:,:,:) = cldext(:,:,:,8)
      used = send_data ( id_ext_cld_nir, cloud, Time_diag, is, js, &
                           1, rmask=mask )
     endif

      if ( id_sct_cld_nir > 0  ) then
        cloud(:,:,:) = cldsct(:,:,:,8)
       used = send_data ( id_sct_cld_nir, cloud, Time_diag, is, js, &
                           1, rmask=mask )
    endif

    if ( id_asymm_cld_nir > 0  ) then
      cloud(:,:,:) = 100.0*cldasymm(:,:,:,8)
      used = send_data ( id_asymm_cld_nir, cloud, Time_diag, is, js, &
                          1, rmask=mask )
    endif

    endif ! end lhsw/esfsw diagnostic option loop

      deallocate (cloud)
!------------------------------------------------------------------



end subroutine cloudrad_netcdf


!####################################################################

subroutine diag_field_init ( Time, axes )

!---------------------------------------------------------------------
type(time_type), intent(in) :: Time
integer        , intent(in) :: axes(4)
!---------------------------------------------------------------------


!------------ initialize diagnostic fields in this module --------------


    id_tot_cld_amt = &
    register_diag_field ( mod_name, 'tot_cld_amt', axes(1:2), Time, &
                         'total cloud amount', 'percent'            )

    id_high_cld_amt = &
     register_diag_field ( mod_name, 'high_cld_amt', axes(1:2), Time, &
                        'high cloud amount', 'percent'            )

     id_mid_cld_amt = &
    register_diag_field ( mod_name, 'mid_cld_amt', axes(1:2), Time, &
                          'mid cloud amount', 'percent'            )
  
      id_low_cld_amt = &
      register_diag_field ( mod_name, 'low_cld_amt', axes(1:2), Time, &
                      'low cloud amount', 'percent'            )

    id_cld_amt = &
    register_diag_field ( mod_name, 'cld_amt', axes(1:3), Time, &
                         'cloud amount', 'percent',             &
                         missing_value=missing_value            )

    id_em_cld_lw = &
      register_diag_field ( mod_name, 'em_cld_lw', axes(1:3), Time, &
                         'lw cloud emissivity', 'percent',        &
                           missing_value=missing_value          )
 
    id_em_cld_10u = &
    register_diag_field ( mod_name, 'em_cld_10u', axes(1:3), Time, &
                          'cloud emissivity 10 um band', 'percent',    &
                       missing_value=missing_value          )


!   id_em_cld = &
!   register_diag_field ( mod_name, 'em_cld', axes(1:3), Time, &
!                        'cloud emissivity', 'percent',        &
!                         missing_value=missing_value          )

    if (trim(swform) == 'lhsw' ) then
      id_alb_uv_cld = &
      register_diag_field ( mod_name, 'alb_uv_cld', axes(1:3), Time, &
                         'UV reflected by cloud', 'percent',       &
                          missing_value=missing_value              )

      id_alb_nir_cld = &
      register_diag_field ( mod_name, 'alb_nir_cld', axes(1:3), Time, &
                         'IR reflected by cloud', 'percent',        &
                          missing_value=missing_value               )

!   --- do not output this field ---
!     id_abs_uv_cld = &
!     register_diag_field ( mod_name, 'abs_uv_cld', axes(1:3), Time, &
!                        'UV absorbed by cloud', 'percent',        &
!                         missing_value=missing_value              )

      id_abs_nir_cld = &
      register_diag_field ( mod_name, 'abs_nir_cld', axes(1:3), Time, &
                         'IR absorbed by cloud', 'percent',         &
                          missing_value=missing_value               )
    endif

     id_ext_cld_uv = &
      register_diag_field ( mod_name, 'ext_cld_uv', axes(1:3), Time, &
                          '.27um cloud extinction coeff', 'km-1',  &
                          missing_value=missing_value          )

     id_sct_cld_uv = &
      register_diag_field ( mod_name, 'sct_cld_uv', axes(1:3), Time, &
                           '.27um cloud scattering coeff', 'km-1', &
                           missing_value=missing_value          )

     id_asymm_cld_uv = &
    register_diag_field ( mod_name, 'asymm_cld_uv', axes(1:3), Time, &
                       '.27um cloud asymmetry parameter', 'percent', & 
                         missing_value=missing_value          )

    id_ext_cld_vis = &
    register_diag_field ( mod_name, 'ext_cld_vis', axes(1:3), Time, &
                      '.55um cloud extinction coeff', 'km-1', &
                        missing_value=missing_value          )

    id_sct_cld_vis = &
     register_diag_field ( mod_name, 'sct_cld_vis', axes(1:3), Time, &
                        '.55um cloud scattering coeff', 'km-1', &
                      missing_value=missing_value          )

     id_asymm_cld_vis = &
  register_diag_field ( mod_name, 'asymm_cld_vis', axes(1:3), Time, &
                    '.55um cloud asymmetry parameter', 'percent', &  
                       missing_value=missing_value          )

     id_ext_cld_nir = &
    register_diag_field ( mod_name, 'ext_cld_nir', axes(1:3), Time, &
               '1.4um cloud extinction coeff', 'km-1', &
			       missing_value=missing_value          )

     id_sct_cld_nir = &
     register_diag_field ( mod_name, 'sct_cld_nir', axes(1:3), Time, &
                         '1.4um cloud scattering coeff', 'km-1', &
                       missing_value=missing_value          )
 
     id_asymm_cld_nir = &
    register_diag_field ( mod_name, 'asymm_cld_nir', axes(1:3), Time, &
                     '1.4um cloud asymmetry parameter', 'percent', & 
                       missing_value=missing_value          )



!-----------------------------------------------------------------------

   end subroutine diag_field_init



!#######################################################################  1357  
subroutine compute_isccp_clds ( pflux, camtsw, hml_ca)

real,  dimension(:,:,:),   intent(in)  :: pflux, camtsw
real,  dimension(:,:,:),   intent(out) :: hml_ca

!
!   define arrays giving the fractional cloudiness for clouds with
!   tops within the ISCCP definitions of high (10-440 hPa), middle
!   (440-680 hPa) and low (680-1000 hPa).
 
!    note that at this point pflux is in cgs units. change this later.

 
! ---------------------------------------------------------------------

   real,  parameter :: mid_btm = 6.8e5,high_btm = 4.4e5
  
! local array

integer :: i, j, k

 
!---- compute high, middle and low cloud amounts assuming that -----
!       independent clouds overlap randomly

    hml_ca = 1.0
 
   do j=JSRAD,JERAD
    do i=ISRAD,IERAD
      do k = KSRAD, KERAD
        if (pflux(i,j,k)  <=  high_btm) then
          hml_ca(i,j,1) = hml_ca(i,j,1) * (1. - camtsw(i,j,k))
        else if ( (pflux(i,j,k) >  high_btm) .and.  &
           (pflux(i,j,k) <=  mid_btm) ) then
         hml_ca(i,j,2) = hml_ca(i,j,2) * (1. - camtsw(i,j,k))
       else  if ( pflux(i,j,k) > mid_btm ) then
	 hml_ca(i,j,3) = hml_ca(i,j,3) * (1. - camtsw(i,j,k))
       endif
    enddo
  enddo
  enddo

    hml_ca = 1. - hml_ca
    hml_ca = 100. * hml_ca
  
end subroutine compute_isccp_clds
!#######################################################################

	       end module cloudrad_package_mod

