
                 module cloudrad_package_mod

use time_manager_mod,        only: time_type
use diag_manager_mod,        only: register_diag_field, send_data
!use rad_output_file_mod,     only: hold_cloud, put_cloudpackage_type
!use rad_output_file_mod,     only:             put_cloudpackage_type
!use donner_ice_mod,          only: inquire_donner_ice,   &
!                                   sw_albedo_zen_angle
use standalone_clouds_mod,   only: standalone_clouds_init, &
			           standalone_clouds_driver
use no_clouds_mod,           only: no_clouds_init, no_clouds_calc
use diag_clouds_W_mod,       only: diag_clouds_init, diag_clouds_calc
use obs_clouds_W_mod,        only: obs_clouds_init, obs_clouds_calc
use zonal_clouds_W_mod,      only: zonal_clouds_init, zonal_clouds_calc
use strat_clouds_W_mod,      only: strat_clouds_init, strat_clouds_calc
use donner_deep_clouds_W_mod, only: donner_deep_clouds_init,  &
                                    donner_deep_clouds_calc
use rh_based_clouds_mod,     only: rh_clouds_init, rh_clouds_calc
use mgrp_prscr_clds_mod,     only: mgrp_prscr_init, mgrp_prscr_calc
use microphys_rad_mod,       only: lwemiss_calc, comb_cldprops_calc, &
                                   microphys_rad_init
use utilities_mod,           only: open_file, file_exist,   &
                                   check_nml_error, error_mesg,   &
                                   print_version_number, FATAL, NOTE, &
				   WARNING, close_file, get_my_pe, &
				   read_data, write_data
!use rad_step_setup_mod,      only: jabs, iabs, pflux, IMINP, IMAXP,   &
!use rad_step_setup_mod,      only: jabs, iabs, pflux,                 &
!use rad_step_setup_mod,      only:             pflux,                 &
!use rad_step_setup_mod,      only:                   &
!			   JMINP, JMAXP, ISRAD, IERAD, JSRAD,  &
!			   JERAD, KSRAD, KERAD,   &
!				   get_std_pressures
use rad_utilities_mod,       only: longwave_control_type, Lw_control, &
                                   Environment, environment_type, &
                                   shortwave_control_type, Sw_control,&
				   cldrad_properties_type, &
				   cld_diagnostics_type, &
				   cld_space_properties_type, &
				   astronomy_type, &
				   atmos_input_type, &
				   radiation_control_type, Rad_control,&
				 cloudrad_control_type, Cldrad_control,&
                                   longwave_parameter_type,  &
				   Lw_parameters
!use longwave_setup_mod,      only: longwave_parameter_type,  &
!				   Lw_parameters
!use radiation_diag_mod,      only: radiag_from_clouds_lh
!use radiation_diag_mod,      only: radiag_from_clouds_lh,   &
!                                  radiag_from_clouds_esf
!use std_pressures_mod,       only: get_std_pressures
!se astronomy_package_mod,   only: get_astronomy_for_clouds,  &
!use astronomy_package_mod,   only:                            &
!			           get_astronomy_for_clouds_init
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
      character(len=128)  :: version =  '$Id: cloudrad_package.F90,v 1.3 2002/07/16 22:34:30 fms Exp $'
      character(len=128)  :: tag     =  '$Name: havana $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
	  cloudrad_package_init, cldmarch,  clouddrvr, &
!  cloudrad_package_alloc, cloudrad_package_dealloc, &
             cloudrad_package_end, &
	                          cloudrad_package_dealloc, &
				  cloudrad_package_alloc_rad, &
!  get_clouds_for_lhsw, get_clouds_for_esfsw, &
	  convert_to_cloud_space                    
!         get_ncldsw, get_clouds_for_lwrad,  &
!         get_ncldsw
!  get_swcldprops_from_cloudrad


private          &
	  default_clouds, cloudrad_netcdf, diag_field_init, &
	  cloudrad_package_alloc,  &
	  compute_isccp_clds


!---------------------------------------------------------------------
!-------- namelist  ---------


character(len=16)            :: microphys_form = '            '
character(len=16)            :: cloud_type_form = '            '
logical :: do_lwcldemiss = .false.

!    logical variables derived from namelist input
logical                      :: do_pred_cld_microphys = .false.
logical                      :: do_presc_cld_microphys   = .false.
logical                      :: do_no_cld_microphys   = .false.

logical                      :: do_rh_clouds=.false.
logical                      :: do_strat_clouds=.false.
logical                      :: do_donner_deep_clouds=.false.
logical                      :: do_zonal_clouds=.false.
logical                      :: do_mgroup_prescribed = .false.
logical                      :: do_obs_clouds=.false.
logical                      :: do_no_clouds=.false.
logical                      :: do_diag_clouds=.false.
logical                      :: do_specified_clouds=.false.

integer                      :: n_prsc_clds=3




namelist /cloudrad_package_nml /     &
                               microphys_form, &
			       do_lwcldemiss, &
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
!                layers from 1     to kx   .
!     crndlw  =  amounts of randomly overlapped longwave clouds in
!                layers from 1     to kx   .
!     nmxolw  =  number of maximally overlapped longwave clouds
!                in each grid column.
!     nrndlw  =  number of maximally overlapped longwave clouds
!                in each grid column.
!     emmxolw =  longwave cloud emissivity for maximally overlapped
!                clouds through layers from 1     to kx   .
!     emrndlw =  longwave cloud emissivity for randomly overlapped
!                clouds through layers from 1     to kx   .
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
!     abscoeff = combined absorption coefficient (Km-1) for  
!                clouds in each of the longwave frequency bands.
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

!real, dimension(:,:,:,:), allocatable ::  cldext, cldasymm, &
!                                          cldsct,  emmxolw, emrndlw, &
!                                          cirabsw, cirrfsw, cvisrfsw
!real,    dimension(:,:,:),allocatable ::  camtsw, cmxolw, crndlw
!integer, dimension(:,:),  allocatable ::  ncldsw, nmxolw, nrndlw

!    cloud property arrays for particular cloud types
!    these may be combined into combined cloud property arrays, which
!    are described in the preceding comments
!
!    cld_lsc    = fractional area of grid box in each vertical layer
!                 occupied by large-scale clouds. obtained using
!                 strat_cloud parameterization (S. Klein).
!    cldext_lsc = extinction coefficient (Km-1) for large-scale 
!                 clouds in each of the shortwave frequency bands.
!    cldsct_lsc = scattering coefficient (Km-1) for large-scale 
!                 clouds in each of the shortwave frequency bands.
!  cldasymm_lsc = asymmetry factor  for large-scale
!                 clouds in each of the shortwave frequency bands.
!  abscoeff_lsc = absorption coefficient (Km-1) for large-scale
!                 clouds in each of the longwave frequency bands.
 
!    cld_cell    = fractional area of grid box in each vertical layer
!                  occupied by cell (convective) clouds. obtained using
!                  donner_deep parameterization (L. Donner).
!    cldext_cell = extinction coefficient (Km-1) for cell 
!                  clouds in each of the shortwave frequency bands.
!    cldsct_cell = scattering coefficient (Km-1) for cell 
!                  clouds in each of the shortwave frequency bands.
!  cldasymm_cell = asymmetry factor  for cell
!                  clouds in each of the shortwave frequency bands.
!  abscoeff_cell = absorption coefficient (Km-1) for cell
!                  clouds in each of the longwave frequency bands.
 
!    cld_meso    = fractional area of grid box in each vertical layer
!                  occupied by meso (anvil) clouds. obtained using
!                  donner_deep parameterization (L. Donner).
!    cldext_meso = extinction coefficient (Km-1) for meso 
!                  clouds in each of the shortwave frequency bands.
!    cldsct_meso = scattering coefficient (Km-1) for meso 
!                  clouds in each of the shortwave frequency bands.
!  cldasymm_meso = asymmetry factor  for meso
!                  clouds in each of the shortwave frequency bands.
!  abscoeff_meso = absorption coefficient (Km-1) for meso
!                  clouds in each of the longwave frequency bands.

!!!RSH later these arrays may be made local to clouddrvr and assed to
!      _netcdf as required.

real, dimension(:,:,:), allocatable :: cld_lsc, cld_cell, cld_meso
real, dimension(:,:,:,:), allocatable :: cldext_lsc, cldsct_lsc,  &
                                         cldasymm_lsc, abscoeff_lsc
real, dimension(:,:,:,:), allocatable :: cldext_cell, cldsct_cell,  &
                                          cldasymm_cell, abscoeff_cell
real, dimension(:,:,:,:), allocatable :: cldext_meso, cldsct_meso,  &
                                          cldasymm_meso, abscoeff_meso
!------------------------------------------------------------------
!    NLWCLDB is the actual number of frequency bands for which lw
!    emissitivies are defined. 
!------------------------------------------------------------------
integer            :: NLWCLDB 
 
integer            :: ix, jx, kx
!integer            :: kmin, kmax
integer            :: nsolwg
!character(len=10)  :: swform
!character(len=16)  :: swform
logical :: do_lhsw, do_esfsw

!-------------------- diagnostics fields ------------------------------

integer :: id_tot_cld_amt, id_cld_amt,   id_em_cld_lw, id_em_cld_10u, & 
           id_abs_lsc_cld_10u, id_abs_lsc_cld_lw,  &
             id_abs_cell_cld_10u, id_abs_cell_cld_lw,  &
            id_abs_meso_cld_10u, id_abs_meso_cld_lw,  &
           id_abs_cld_10u, id_abs_cld_lw,  &
           id_high_cld_amt, id_mid_cld_amt, id_low_cld_amt,  &
       id_lsc_cld_amt, id_cell_cld_amt, id_meso_cld_amt,  &
     id_lsc_cld_ext_uv, id_lsc_cld_ext_vis, id_lsc_cld_ext_nir, &
     id_lsc_cld_sct_uv, id_lsc_cld_sct_vis, id_lsc_cld_sct_nir, &
  id_lsc_cld_asymm_uv, id_lsc_cld_asymm_vis, id_lsc_cld_asymm_nir, &
    id_cell_cld_ext_uv, id_cell_cld_ext_vis, id_cell_cld_ext_nir, &
     id_cell_cld_sct_uv, id_cell_cld_sct_vis, id_cell_cld_sct_nir, &
   id_cell_cld_asymm_uv, id_cell_cld_asymm_vis, id_cell_cld_asymm_nir, &
     id_meso_cld_ext_uv, id_meso_cld_ext_vis, id_meso_cld_ext_nir, &
    id_meso_cld_sct_uv, id_meso_cld_sct_vis, id_meso_cld_sct_nir, &
   id_meso_cld_asymm_uv, id_meso_cld_asymm_vis, id_meso_cld_asymm_nir, &
           id_ext_cld_uv,   id_sct_cld_uv,  id_asymm_cld_uv, &
           id_ext_cld_vis,  id_sct_cld_vis, id_asymm_cld_vis, &
           id_ext_cld_nir,  id_sct_cld_nir, id_asymm_cld_nir, &
           id_alb_uv_cld, id_alb_nir_cld, id_abs_uv_cld, id_abs_nir_cld


character(len=8), parameter :: mod_name = 'cloudrad'

real :: missing_value = -999.


logical :: cloudrad_package_initialized= .false.

!----------------------------------------------------------------------
!----------------------------------------------------------------------



             contains 





!subroutine cloudrad_package_init (kmin_in, kmax_in, kx_in, &
subroutine cloudrad_package_init (                         &
!			  th, pd, axes, Time, qlyr, lonb, latb)
!			      pd, axes, Time, qlyr, lonb, latb)
!			      pd, axes, Time,       lonb, latb)
!	      pd, lonb, latb, axes, Time)
		      pref, lonb, latb, axes, Time)

!------------------------------------------------------------------
!integer,               intent(in)             :: kmin_in, kmax_in, kx_in
!eal,    dimension(:), intent(in)             :: th, pd   
!real,    dimension(:), intent(in)             ::     pd, lonb, latb
real,    dimension(:), intent(in)             ::     lonb, latb
real,    dimension(:,:), intent(in)             ::     pref
!real,    dimension(:), intent(in), optional   :: qlyr, lonb, latb
!real,    dimension(:), intent(in)             ::        lonb, latb
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
!   real, dimension(:), allocatable   :: qlevel, pd
    real, dimension(:), allocatable   :: qlevel
!   real, dimension (size(latb,1)-1 ) :: th
    real                              :: psbar
    integer             :: kx_in

    if (cloudrad_package_initialized) return

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

  Lw_control%do_lwcldemiss = do_lwcldemiss
 
!--------------------------------------------------------------------
!     NLWCLDB =  number of infrared bands with varying cloud
!                properties (emissivity). at present either unity
!                or a model-dependent amount (ifdef lwcldemiss)
!--------------------------------------------------------------------
     if (Lw_control%do_lwcldemiss) then
       NLWCLDB = 7
     else
       NLWCLDB = 1
     endif
 
!     Lw_parameters%NLWCLDB = NLWCLDB
      Cldrad_control%nlwcldb = nlwcldb

!-------------------------------------------------------------------
!  check for consistency between the specified radiation options and
!  the specified microphysics 
!----------------------------------------------------------------------

!   swform = Sw_control%sw_form
    do_lhsw = Sw_control%do_lhsw
    do_esfsw = Sw_control%do_esfsw

!   if (trim(swform) == 'esfsw99' .and.   &
    if (do_esfsw                  .and.   &
        trim(microphys_form) == 'none') then 
      call error_mesg( 'cloudrad_package_init',  &
      ' must specify microphysics when using esfsw99 shortwave.', FATAL)
    endif

    if (Environment%running_gcm) then
!   if (trim(swform) == 'lhsw' .and.     &
    if (do_lhsw                .and.     &
	      .not. Lw_control%do_lwcldemiss .and.    &
	trim(microphys_form) /= 'none') then 
      call error_mesg( 'cloudrad_package_init',  &
       ' not using microphysics -- set microphys_form to none.', FATAL)
    endif
    else 
!   if (trim(swform) == 'lhsw' .and.     &
    if (do_lhsw                .and.     &
	            Lw_control%do_lwcldemiss .and.    &
	trim(cloud_type_form) /= 'strat') then 
      call error_mesg( 'cloudrad_package_init',  &
       ' cannot currently run combination of lhsw and lwcldemiss &
                   &in standalone.',  FATAL)
    endif
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
!       if (Environment%running_fms) then
          do_strat_clouds = .true.
!       else if (Environment%running_skyhi) then
!         call error_mesg( 'cloudrad_package_init',  &
!               ' strat clouds not available in SKYHI.', FATAL)
!       endif


!-------------------------------------------------------------------
!  cloud fractions, heights are predicted by the model based on donner  
!  deep cloud (cell cloud, anvil cloud) scheme
!-------------------------------------------------------------------
 
      else if (trim(cloud_type_form) == 'deep')  then
!        if (Environment%running_fms) then
           do_donner_deep_clouds = .true.
!        else if (Environment%running_skyhi) then
!          call error_mesg( 'cloudrad_package_init',  &
!                ' donner_deep clouds not available in SKYHI.', FATAL)
!        endif

!-------------------------------------------------------------------
!  cloud fractions, heights are predicted by the model based on donner  
!  deep cloud (cell cloud, anvil cloud) and klein large-scale cloud
!  scheme.
!-------------------------------------------------------------------

      else if (trim(cloud_type_form) == 'stratdeep')  then
!       if (Environment%running_fms) then
           do_strat_clouds = .true.
           do_donner_deep_clouds = .true.
!       else if (Environment%running_skyhi) then
!          call error_mesg( 'cloudrad_package_init',  &
!                ' strat and donner_deep clouds not available in SKYHI.', FATAL)
!       endif


!-------------------------------------------------------------------
!  cloud fractions, heights are prescribed in the model using fms form
!-------------------------------------------------------------------

      else if (trim(cloud_type_form) == 'zonal')  then
!       if (Environment%running_fms) then
          do_zonal_clouds = .true.
!       else if (Environment%running_skyhi) then
!         call error_mesg( 'cloudrad_package_init',  &
!                     ' zonal clouds not available in SKYHI.', FATAL)
!       endif

!-------------------------------------------------------------------
!  cloud fractions, heights are based on observations used by FMS
!-------------------------------------------------------------------

      else if (trim(cloud_type_form) == 'obs')  then
!       if (Environment%running_fms) then
          do_obs_clouds = .true.
!       else if (Environment%running_skyhi) then
!         call error_mesg( 'cloudrad_package_init',  &
!                        ' obs clouds not available in SKYHI.', FATAL)
!       endif

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
!       else if (do_strat_clouds) then
        else if (do_strat_clouds .or. do_donner_deep_clouds) then
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
!       ' predicted microphys not available with gordon diag clouds', &
!   FATAL)
    ' predicted microphys under development with gordon diag clouds', &
         NOTE)
          do_pred_cld_microphys = .true.
        else if (do_no_clouds) then
        endif

      else if (trim(microphys_form) == 'prescribed') then
!----------------------------------------------------------------------
!  cloud species concentrations, sizes, etc are prescribed (for
!  low, middle, high clouds)
!----------------------------------------------------------------------

        if (do_rh_clouds) then
          do_presc_cld_microphys = .true.
!       else if (do_strat_clouds) then
         else if (do_strat_clouds .or. do_donner_deep_clouds) then
          call error_mesg( 'cloudrad_package_init',  &
        ' use predicted microphys with the strat cloud  or donner&
         & deep cloud module',  FATAL)
!       ' use predicted microphys with the strat cloud module',  &
!  FATAL)
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
!      ' prescribed microphys not available with gordon diag clouds', &
!						    FATAL)
       ' prescribed microphys under develop with gordon diag clouds', &
                                                    NOTE)
          do_presc_cld_microphys = .true.
        else if (do_no_clouds) then
        endif

      else if (trim(microphys_form) == 'none') then
!----------------------------------------------------------------------
!  cloud reflectivities, absorptivities, emissivities are prescribed
!  for low, middle, high clouds (legacy code)
!----------------------------------------------------------------------

        if (do_rh_clouds) then
          do_no_cld_microphys = .true.
!       else if (do_strat_clouds) then
        else if (do_strat_clouds .or. do_donner_deep_clouds) then
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

  else if (trim(cloud_type_form) == 'deep')  then
!----------------------------------------------------------------------
!  cloud fractions, heights are predicted by the model based on donner  
!  deep cloud (cell cloud, anvil cloud) scheme
!-------------------------------------------------------------------

           do_donner_deep_clouds = .true.
 
    else if (trim(cloud_type_form) == 'stratdeep')  then
!-------------------------------------------------------------------
!  cloud fractions, heights are predicted by the model based on donner  
!   deep cloud (cell cloud, anvil cloud) and klein large-scale cloud
!  scheme.
!-------------------------------------------------------------------
 
         do_strat_clouds = .true.
         do_donner_deep_clouds = .true.



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

!        if (do_strat_clouds) then
         if (do_strat_clouds .or. do_donner_deep_clouds) then
           do_pred_cld_microphys = .true.
         else if (do_specified_clouds) then
           call error_mesg( 'cloudrad_package_init',  &
!             'if microphys_form is predicted, must activate &
!          &strat clouds', FATAL)
              'if microphys_form is predicted, must activate &
	          &strat clouds or donner_deep clouds', FATAL)
         endif
 
       else if (trim(microphys_form) == 'prescribed') then
!----------------------------------------------------------------------
!  cloud species concentrations, sizes, etc are prescribed (for
!  low, middle, high clouds)
!----------------------------------------------------------------------

!        if (do_strat_clouds) then
         if (do_strat_clouds .or. do_donner_deep_clouds) then
           call error_mesg( 'cloudrad_package_init',  &
            ' prescribed microphys not available with strat clouds &
             &or donner_deep_clouds', FATAL)
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

!    kmin = kmin_in
!    kmax = kmax_in
!    kmax = size(pd,1) - 1
!   call get_astronomy_for_clouds_init (nsolwg)
    nsolwg = 1

!---------------------------------------------------------------------
!    save the sigma levels that have been input for later use within 
!    this module.
!---------------------------------------------------------------------
!   allocate (qlevel (KSRAD:KERAD) )
!   kx_in = size(pd,1) - 1
    kx_in = size(pref,1) - 1
    allocate (qlevel (1:kx_in    ) )
!4  if (present(qlyr)) then
!     qlevel(:) = qlyr(KSRAD:KERAD)
!4    qlevel(:) = qlyr(1:kx_in    )
!4  else
!     allocate ( pd(KMIN:KMAX+1) )
!     call get_std_pressures (pd_out=pd)
      psbar = 1.0E-03*pstd  ! convert from cgs to mb
!     qlevel(KSRAD:KERAD) = pd(KSRAD:KERAD)/psbar
!     qlevel(1:kx_in    ) = pd(1:kx_in    )/psbar
!     qlevel(1:kx_in    ) = (pd(1:kx_in    )*1.0E-02)/psbar
      qlevel(1:kx_in    ) = (pref(1:kx_in,1  )*1.0E-02)/psbar
!     deallocate (pd)
!4  endif

!---------------------------------------------------------------------
!    define the number of cloud emissivity bands for use in this module.
!---------------------------------------------------------------------

!   NLWCLDB = Lw_parameters%NLWCLDB

!-------------------------------------------------------------------
!  initialize the specific cloud scheme selected for this run
!-------------------------------------------------------------------

    if (Environment%running_gcm) then
      if (do_mgroup_prescribed) then
!       call mgrp_prscr_init (kx_in, pd, latb)
        call mgrp_prscr_init (kx_in, pref, latb)
      endif

      if (do_rh_clouds) then
!       call rh_clouds_init (th, qlevel, latb)
        call rh_clouds_init (    qlevel, latb)
      endif

      if (do_no_clouds) then
        call no_clouds_init 
      endif

!     if (Environment%running_fms) then
        if (do_zonal_clouds) then
          call zonal_clouds_init
        endif

        if (do_obs_clouds) then
          call obs_clouds_init (lonb, latb)
        endif

        if (do_strat_clouds) then
          call strat_clouds_init
        endif

        if (do_donner_deep_clouds) then
          call donner_deep_clouds_init
        endif

        if (do_strat_clouds .OR. do_donner_deep_clouds) then
!         if (trim(swform) == 'esfsw99' .or.               &
          if (do_esfsw                  .or.               &
                          Lw_control%do_lwcldemiss) then
              call microphys_rad_init
           endif
        endif
 

	if (do_diag_clouds) then
          ix = size(lonb) - 1
	  iy = size(latb) - 1
!  call diag_clouds_init (ix, iy, kmax, th, ierr)
!  call diag_clouds_init (ix, iy, kmax,     ierr)
	  call diag_clouds_init (ix, iy,           ierr)
	endif
!     endif
    else if (Environment%running_standalone) then
      if (do_strat_clouds) then
        call strat_clouds_init
      endif

      if (do_donner_deep_clouds) then
        call donner_deep_clouds_init
       endif


!      call standalone_clouds_init (              th,                &
!           do_strat_clouds, do_donner_deep_clouds, do_specified_clouds)

      call standalone_clouds_init ( kx_in, latb, do_strat_clouds,  &
                        do_donner_deep_clouds, do_specified_clouds)
    endif

!-------------------------------------------------------------------
!  initialize diagnostics of this module
!-------------------------------------------------------------------
!   if (Environment%running_gcm .and. Environment%running_fms  &
    if (Environment%running_gcm                                &
				 .and. .not. do_no_clouds) then
      call diag_field_init (Time, axes)
    endif

!-------------------------------------------------------------------
!  pass cloud package type to rad_output_file_mod
!-------------------------------------------------------------------
!   call put_cloudpackage_type (cloud_type_form)

!--------------------------------------------------------------------


         Cldrad_control%do_rh_clouds = do_rh_clouds
         Cldrad_control%do_strat_clouds = do_strat_clouds
         Cldrad_control%do_donner_deep_clouds = do_donner_deep_clouds
         Cldrad_control%do_zonal_clouds = do_zonal_clouds
         Cldrad_control%do_mgroup_prescribed = do_mgroup_prescribed
         Cldrad_control%do_obs_clouds = do_obs_clouds
         Cldrad_control%do_no_clouds =  do_no_clouds
         Cldrad_control%do_diag_clouds = do_diag_clouds
         Cldrad_control%do_specified_clouds = do_specified_clouds



    cloudrad_package_initialized= .true.


end subroutine cloudrad_package_init



!#################################################################

subroutine cldmarch

!---------------------------------------------------------------------
!     cldmarch calls routines to handle temporal variation of cloud
!     variables (currently null).
!---------------------------------------------------------------------


end subroutine cldmarch



!####################################################################

subroutine clouddrvr (is, ie, js, je, lat,  Rad_time,   &
!                     Atmos_input, land, Astro, &
                      Atmos_input,       Astro, &
!                     pflux_in, deltaz, &
!                     land, cosz, cloud_water, cloud_ice,  press_in, &
!	              temp, rh2o,  &
                      Cldrad_props, Cld_diagnostics, Time_next,  &
                      kbot, mask)

integer,                 intent(in)           ::  is, ie, js, je
!real, dimension(:,:),    intent(in)           ::  lat, land, cosz
!real, dimension(:,:),    intent(in)           ::  lat, land       
real, dimension(:,:),    intent(in)           ::  lat
type(atmos_input_type), intent(in)             :: Atmos_input
type(astronomy_type), intent(in)              :: Astro
!real, dimension(:,:,:),    intent(in)         ::  pflux_in, deltaz, &
!                                             cloud_water, cloud_ice, &
!                                                  press_in, temp, rh2o
type(time_type),         intent(in)           ::  Rad_time
type(cldrad_properties_type), intent(inout)   :: Cldrad_props
type(cld_diagnostics_type), intent(inout)   :: Cld_diagnostics
type(time_type),         intent(in), optional ::  Time_next
integer, dimension(:,:), intent(in), optional ::  kbot
real, dimension(:,:,:),  intent(in), optional ::  mask
!-------------------------------------------------------------------
 
!    real, dimension(size(press_in,1), size(press_in,2), &
!                     size(press_in,3)) :: press, pflux

     real, dimension(size(Atmos_input%press,1),   &
                          size(Atmos_input%press,2), &
                      size(Atmos_input%press,3)) :: press, pflux, temp

     real, dimension(size(Atmos_input%press,1),   &
                          size(Atmos_input%press,2), &
                      size(Atmos_input%press,3) -1 ) :: deltaz, &
                  cloud_water, cloud_ice, rh2o

     real, dimension(size(Atmos_input%press,1),   &
                          size(Atmos_input%press,2) ) :: land

      integer   :: ix, jx, kx
      integer :: unit,idim,jdim,kdim,ndim
      logical :: al


      ix = size (Atmos_input%press,1)
      jx = size (Atmos_input%press,2)
      kx = size (Atmos_input%press,3) - 1

      deltaz = Atmos_input%deltaz
      cloud_water = Atmos_input%cloud_water
      cloud_ice = Atmos_input%cloud_ice
      temp      = Atmos_input%temp         
      rh2o      = Atmos_input%rh2o       
      land    = Atmos_input%land
!-----------------------------------------------------------------
!  convert press and pflux to cgs.
!     if (Environment%running_gcm) then
          press(:,:,:) = 10.0*Atmos_input%press (:,:,:)
          pflux(:,:,:) = 10.0*Atmos_input%pflux (:,:,:)
!     else
!         press(:,:,:) =  10.0*(    Atmos_input%press (:,:,:))
!         pflux(:,:,:) =  10.0*(    Atmos_input%pflux (:,:,:))
!     endif

!    call cloudrad_package_alloc (ie-is+1, je-js+1, kmax-kmin+1,  &
     call cloudrad_package_alloc (ix, jx, kx,                     &
                                  Cldrad_props)
!--------------------------------------------------------------------
!     initialize the cloud property arrays
!-------------------------------------------------------------------
      call default_clouds (Cldrad_props)

!--------------------------------------------------------------------
!     call desired cloud package to obtain needed cloud radiative 
!     property fields.
!-------------------------------------------------------------------
      if (Environment%running_gcm) then
        if (do_rh_clouds) then
!         if (trim(swform) /= 'esfsw99') then
          if (.not. do_esfsw           ) then
!           call rh_clouds_calc (is, ie, js, je, deltaz, cosz,  &
            call rh_clouds_calc (is, ie, js, je, Cld_diagnostics, deltaz, Astro%cosz,  &
	                        press, temp,    &
				Cldrad_props%camtsw,  &
				Cldrad_props%cmxolw,  &
	                         Cldrad_props%crndlw,   &
				     Cldrad_props%ncldsw ,   &
		                      Cldrad_props%nmxolw ,  &
				      Cldrad_props%nrndlw ,  &
				  Cldrad_props%emmxolw,  &
				   Cldrad_props%emrndlw, &
	   		         cirabsw= Cldrad_props%cirabsw, &
				 cvisrfsw= Cldrad_props%cvisrfsw, &
			         cirrfsw= Cldrad_props%cirrfsw) 
          else
!           call rh_clouds_calc (is, ie, js, je, deltaz, cosz, &
            call rh_clouds_calc (is, ie, js, je, Cld_diagnostics, deltaz, Astro%cosz, &
	                         press, temp,  Cldrad_props%camtsw,   &
				  Cldrad_props%cmxolw,  &
	                         Cldrad_props%crndlw,   &
				     Cldrad_props%ncldsw ,   &
	  		              Cldrad_props%nmxolw ,   &
				      Cldrad_props%nrndlw ,  &
				  Cldrad_props%emmxolw,   &
				  Cldrad_props%emrndlw, &
			         cldext= Cldrad_props%cldext,  &
				 cldsct= Cldrad_props%cldsct,       &
			         cldasymm= Cldrad_props%cldasymm)
          endif
        endif
        if (do_no_clouds) then
	  call no_clouds_calc     &
	                      ( Cldrad_props%ncldsw,  Cldrad_props%nrndlw,  &
			       Cldrad_props%camtsw,  Cldrad_props%crndlw,  &
			       Cldrad_props%emrndlw)
        endif
        if (do_mgroup_prescribed) then
!         if (trim(swform) /= 'esfsw99') then
          if (.not. do_esfsw           ) then
	    call mgrp_prscr_calc ( is, ie,  &
                     js, je, Cld_diagnostics, deltaz, press, temp,   &
		      Cldrad_props%camtsw,  Cldrad_props%cmxolw,   &
		      Cldrad_props%crndlw,  Cldrad_props%ncldsw, &
		      Cldrad_props%nmxolw, &
			      Cldrad_props%nrndlw,              &
			      Cldrad_props%emmxolw,  &
			      Cldrad_props%emrndlw,  &
			     cirabsw= Cldrad_props%cirabsw, &
			     cvisrfsw= Cldrad_props%cvisrfsw,  &
			     cirrfsw= Cldrad_props%cirrfsw)
          else
	    call mgrp_prscr_calc  ( is, ie,    &
                     js, je, Cld_diagnostics, deltaz, press, temp,  &
		      Cldrad_props%camtsw,  Cldrad_props%cmxolw,    &
		      Cldrad_props%crndlw,  Cldrad_props%ncldsw,  &
		      Cldrad_props%nmxolw, &
			      Cldrad_props%nrndlw,            &
			      Cldrad_props%emmxolw, &
			      Cldrad_props%emrndlw,  &
			     cldext= Cldrad_props%cldext, &
			     cldsct= Cldrad_props%cldsct,  &
			     cldasymm= Cldrad_props%cldasymm)
          endif
        endif
!       if (Environment%running_fms) then
          if (do_strat_clouds) then
!           if (trim(swform) /= 'esfsw99') then
            if (.not. do_esfsw           ) then
	    call strat_clouds_calc  &
                           (is, ie, js, je,  Cld_diagnostics, pflux, deltaz, land,  &
!		   cosz, cloud_water, cloud_ice, press, temp, &
			   Astro%cosz, cloud_water, cloud_ice, press, temp, &
			   rh2o,   &
			    camtsw=Cldrad_props%camtsw,  cmxolw=Cldrad_props%cmxolw,   &
			    crndlw=Cldrad_props%crndlw, &
                              ncldsw=Cldrad_props%ncldsw,  &
			      nmxolw=Cldrad_props%nmxolw,  &
			      nrndlw=Cldrad_props%nrndlw,  &
			      emmxolw=Cldrad_props%emmxolw,  &
			      emrndlw=Cldrad_props%emrndlw,  &
			     cirabsw= Cldrad_props%cirabsw,&
			     cvisrfsw= Cldrad_props%cvisrfsw,  &
			     cirrfsw= Cldrad_props%cirrfsw, &
			     Time_next=Time_next)
             else
              call strat_clouds_calc  &
                            (is, ie, js, je, Cld_diagnostics, pflux,  deltaz, land, &
!                          cosz, cloud_water, cloud_ice, press, temp, &
                           Astro%cosz, cloud_water, cloud_ice, press, temp, &
			            rh2o, &
                    cldamt_out=cld_lsc, &
!		      Cldrad_props%camtsw,  Cldrad_props%cmxolw,  &
!		      Cldrad_props%crndlw,  &
                               ncldsw=Cldrad_props%ncldsw,  &
!		       Cldrad_props%nmxolw, &
!                             Cldrad_props%nrndlw,  &
!		      Cldrad_props%emmxolw,  &
!		      Cldrad_props%emrndlw,           &
                          Time_next=Time_next,                  &
!                          cldext= Cldrad_props%cldext, &
!                         cldsct= Cldrad_props%cldsct,  &
!		  cldasymm= Cldrad_props%cldasymm)
                           cldext= cldext_lsc, &
                          cldsct= cldsct_lsc,  &
                           cldasymm= cldasymm_lsc, &
                            abscoeff=abscoeff_lsc)
!----------------------------------------------------------------------
!   as of 14 dec 2000, all cloud properties and cloud fractions are
!   assumed to be for randomly overlapped clouds
!----------------------------------------------------------------------
                Cldrad_props%nrndlw = Cldrad_props%ncldsw
                Cldrad_props%nmxolw = 0.0                   
             endif
           endif

  if (do_donner_deep_clouds) then
!        if (trim(swform) /= 'esfsw99') then
         if (.not. do_esfsw           ) then
         else
           al = allocated(cld_cell)
           idim=size(cld_cell,1)
           jdim=size(cld_cell,2)
           kdim=size(cld_cell,3)
          call donner_deep_clouds_calc             &
                 (is,ie,js,je,deltaz,press,temp, Cld_diagnostics, &
                cld_cell=cld_cell,   &
              cldext_cell=cldext_cell,    &
                  cldsct_cell=cldsct_cell,     &
                 cldasymm_cell=cldasymm_cell,    &
                  abscoeff_cell=abscoeff_cell,   &
                 cld_meso=cld_meso,   &
              cldext_meso=cldext_meso,    &
             cldsct_meso=cldsct_meso,    &
             cldasymm_meso=cldasymm_meso,     &
               abscoeff_meso=abscoeff_meso)
         endif
     endif


!         else if (do_zonal_clouds) then
               if (do_zonal_clouds) then
	    call zonal_clouds_calc   &
	                    (Rad_time, lat, pflux, &
			       Cldrad_props%camtsw,  Cldrad_props%cmxolw,  &
			       Cldrad_props%crndlw,  Cldrad_props%ncldsw,  &
			       Cldrad_props%nmxolw, &
			      Cldrad_props%nrndlw,  &
			      Cldrad_props%cirabsw,  &
			      Cldrad_props%cvisrfsw,  &
			      Cldrad_props%cirrfsw, &
			      Cldrad_props%emmxolw,  &
			      Cldrad_props%emrndlw)
          else if (do_obs_clouds) then
	    call obs_clouds_calc  &
	                    (Rad_time, lat, pflux,  &
			      Cldrad_props%camtsw,  Cldrad_props%cmxolw,  &
			      Cldrad_props%crndlw,  Cldrad_props%ncldsw,  &
			      Cldrad_props%nmxolw, &
			      Cldrad_props%nrndlw,  Cldrad_props%cirabsw,  &
			      Cldrad_props%cvisrfsw,  Cldrad_props%cirrfsw, &
			      Cldrad_props%emmxolw,  &
			      Cldrad_props%emrndlw,  &
			     is, ie, js, je)
	  else if (do_diag_clouds) then
!           if (trim(swform) /= 'esfsw99') then
            if (.not. do_esfsw           ) then
	    call diag_clouds_calc ( is, ie, js, je, Cld_diagnostics, &
!                     lat, pflux, deltaz, cosz, press, temp,  &
                      lat, pflux, deltaz, Astro%cosz, press, temp,  &
		       Cldrad_props%camtsw,  Cldrad_props%cmxolw,   &
		       Cldrad_props%crndlw,  Cldrad_props%ncldsw,   &
		       Cldrad_props%nmxolw, &
			      Cldrad_props%nrndlw,                             &
			      Cldrad_props%emmxolw,   &
			      Cldrad_props%emrndlw,            &
			     Time_next=Time_next,  &
			     cvisrfsw= Cldrad_props%cvisrfsw, &
                             cirrfsw= Cldrad_props%cirrfsw,  &
			     cirabsw= Cldrad_props%cirabsw   )
             else
	    call diag_clouds_calc  (is, ie, js, je, Cld_diagnostics, &
!                    lat, pflux, deltaz, cosz, press, temp,  &
                     lat, pflux, deltaz, Astro%cosz, press, temp,  &
		      Cldrad_props%camtsw,  Cldrad_props%cmxolw,   &
		      Cldrad_props%crndlw,  Cldrad_props%ncldsw,   &
		      Cldrad_props%nmxolw, &
			      Cldrad_props%nrndlw,                             &
			      Cldrad_props%emmxolw,  &
			      Cldrad_props%emrndlw,           &
			     Time_next=Time_next,  &
			     cldext= Cldrad_props%cldext, &
			     cldsct= Cldrad_props%cldsct,   &
			     cldasymm= Cldrad_props%cldasymm)
!----------------------------------------------------------------------
!   as of 14 dec 2000, all cloud properties and cloud fractions are
!   assumed to be for randomly overlapped clouds
!----------------------------------------------------------------------
                Cldrad_props%nrndlw = Cldrad_props%ncldsw
                Cldrad_props%nmxolw = 0.0                   

             endif
          endif
!       endif
      else if (Environment%running_standalone) then
!       if (trim(swform) /= 'esfsw99') then
        if (.not. do_esfsw           ) then
          call standalone_clouds_driver (   is, ie, js, je, Cld_diagnostics, pflux, &
!                       deltaz, land, cosz, cloud_water,  &
                        deltaz, land, Astro%cosz, cloud_water,  &
			cloud_ice, press, temp, rh2o,do_strat_clouds, &
                          do_donner_deep_clouds, &
	   		       Cldrad_props%camtsw,  Cldrad_props%cmxolw,    &
			       Cldrad_props%crndlw,  Cldrad_props%ncldsw,   &
			       Cldrad_props%nmxolw, &
		  	        Cldrad_props%nrndlw,  &
			        Cldrad_props%emmxolw,   &
			        Cldrad_props%emrndlw,   &
 		       Time_next=Time_next,  &
			       cirabsw= Cldrad_props%cirabsw,  &
			       cvisrfsw= Cldrad_props%cvisrfsw, &
			       cirrfsw= Cldrad_props%cirrfsw)
        else
!!! this code only works when cloudtype = strat and microphysform 
!!! == predicted . code is missing to retrieve needed values of
!!!  cloud radiative properties for cloudtype = specified and micro-
!!!!  phys_form = prescribed . in that case cld_lsc is undefined, eg.
          call standalone_clouds_driver  (is, ie, js, je,  Cld_diagnostics, pflux, &
!                       deltaz, land, cosz, cloud_water, &
                        deltaz, land, Astro%cosz, cloud_water, &
			cloud_ice, press, temp, rh2o,do_strat_clouds,  &
                         do_donner_deep_clouds, &
!			        Cldrad_props%camtsw,  Cldrad_props%cmxolw,   &
                        camtsw=cld_lsc,       cmxolw=Cldrad_props%cmxolw,   &
			        crndlw=Cldrad_props%crndlw,  &
				ncldsw=Cldrad_props%ncldsw,   &
			        nmxolw=Cldrad_props%nmxolw, &
		  	        nrndlw=Cldrad_props%nrndlw,  & 
			        emmxolw=Cldrad_props%emmxolw,   &
			        emrndlw=Cldrad_props%emrndlw,   &
 		       Time_next=Time_next,  &
!		       cldext= Cldrad_props%cldext,  &
!		       cldsct= Cldrad_props%cldsct,       &
!		       cldasymm= Cldrad_props%cldasymm)
			       cldext= cldext_lsc,  &
			       cldsct= cldsct_lsc,       &
			       cldasymm= cldasymm_lsc, &
			       abscoeff=abscoeff_lsc)
        endif
      endif

!    combine (if necessary) cloud properties from multiple cloud
!     types into 1 set of cloud property arrays
!
!     compute lw emissivity from lw abs coeff
!
!    if (trim(swform) == 'esfsw99') then
     if (do_esfsw                 ) then
       if ( do_strat_clouds .and. do_donner_deep_clouds) then
!     unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'microphys_rad', version_number)
!     if (get_my_pe() == 0)  then
!       write (unit,*) ' cloudrad_package'
!       write (unit,*) ' cld_cell'
!write (unit,*) al
!write (unit,*) idim,jdim,kdim
!write (unit,*) cld_cell
!       write (unit,*) ' cld_meso'
!write (unit,*) cld_meso
!     endif
!     call close_file (unit)
        call comb_cldprops_calc(                               &
!          cldext, cldsct, cldasymm, abscoeff,   &
           Cldrad_props%cldext, Cldrad_props%cldsct,   &
	   Cldrad_props%cldasymm, Cldrad_props%abscoeff,   &
           cld_lsc=cld_lsc,  &
           cldext_lsc=cldext_lsc, cldsct_lsc=cldsct_lsc,  &
          cldasymm_lsc=cldasymm_lsc, abscoeff_lsc=abscoeff_lsc,  &
           cld_cell=cld_cell,  &
           cldext_cell=cldext_cell, cldsct_cell=cldsct_cell,  &
           cldasymm_cell=cldasymm_cell, abscoeff_cell=abscoeff_cell,  &
           cld_meso=cld_meso,  &
         cldext_meso=cldext_meso, cldsct_meso=cldsct_meso,  &
          cldasymm_meso=cldasymm_meso, abscoeff_meso=abscoeff_meso)
     else if (do_strat_clouds) then
       call comb_cldprops_calc(                               &
!          cldext, cldsct, cldasymm, abscoeff,   &
           Cldrad_props%cldext, Cldrad_props%cldsct,   &
	   Cldrad_props%cldasymm, Cldrad_props%abscoeff,   &
            cld_lsc=cld_lsc,  &
            cldext_lsc=cldext_lsc, cldsct_lsc=cldsct_lsc,  &
            cldasymm_lsc=cldasymm_lsc, abscoeff_lsc=abscoeff_lsc)
      else if (do_donner_deep_clouds) then
!     unit = open_file ('fort.145', action='append',threading='multi')
!       write (unit,*) ' cloudrad_package before comp'
!write (unit,*) 'jabs,iabs'
!write (unit,*) jabs,iabs
!       write (unit,*) ' cld_cell'
!write (unit,*) al
!write (unit,*) idim,jdim,kdim
!write (unit,*) cld_cell
!       write (unit,*) ' cld_meso'
!write (unit,*) cld_meso
!write (unit,*) 'cldasymm_cell'
!write (unit,*) cldasymm_cell
!write (unit,*) 'cldasymm_meso'
!write (unit,*) cldasymm_meso
!write (unit,*) 'cldsct_cell'
!write (unit,*) cldsct_cell
!write (unit,*) 'cldsct_meso'
!write (unit,*) cldsct_meso
!write (unit,*) 'cldext_cell'
!write (unit,*) cldext_cell
!write (unit,*) 'cldext_meso'
!write (unit,*) cldext_meso
!     call close_file (unit)
        call comb_cldprops_calc(                               &
!          cldext, cldsct, cldasymm, abscoeff,   &
           Cldrad_props%cldext, Cldrad_props%cldsct,   &
	   Cldrad_props%cldasymm, Cldrad_props%abscoeff,   &
            cld_cell=cld_cell,  &
           cldext_cell=cldext_cell, cldsct_cell=cldsct_cell,  &
            cldasymm_cell=cldasymm_cell, abscoeff_cell=abscoeff_cell,  &
            cld_meso=cld_meso,  &
           cldext_meso=cldext_meso, cldsct_meso=cldsct_meso,  &
           cldasymm_meso=cldasymm_meso, abscoeff_meso=abscoeff_meso)
       endif

       if ( do_strat_clouds .or. do_donner_deep_clouds) then
       call lwemiss_calc ( deltaz,                              &
!                          abscoeff,                           &
!                          cldemiss)
                           Cldrad_props%abscoeff,               &
                           Cldrad_props%cldemiss)
!     unit = open_file ('fort.145', action='append',threading='multi')
!       write (unit,*) ' cloudrad_package after comp'
!idim = size(abscoeff,1)
!jdim = size(abscoeff,2)
!kdim = size(abscoeff,3)
!ndim = size(abscoeff,4)
!write (unit,*) 'jabs= ', jabs
!write (unit,*) 'iabs= ', iabs
!write (unit,*) 'idim= ', idim, 'jdim = ', jdim, 'kdim = ', kdim
!write (unit,*) ' ndim = ', ndim
!       write (unit,*) ' cldext'
!       write (unit,*)  cldext
!       write (unit,*) ' cldsct'
!       write (unit,*)  cldsct
!       write (unit,*) ' cldasymm'
!       write (unit,*)  cldasymm
!       write (unit,*) ' abscoeff'
!       write (unit,*)  abscoeff
!       write (unit,*) ' cldemiss'
!       write (unit,*)  cldemiss
!       write (unit,*) ' abscoeff_cell'
!       write (unit,*)  abscoeff_cell
!       write (unit,*) ' abscoeff_meso'
!       write (unit,*)  abscoeff_meso
!      call close_file (unit)
!
!   as of 3 august 2001, we are assuming randomly overlapped clouds
!  only. the cloud and emissivity settings follow:

!       emmxolw = cldemiss
!       emrndlw = cldemiss
!       cmxolw = 0.0
        Cldrad_props%emmxolw = Cldrad_props%cldemiss
        Cldrad_props%emrndlw = Cldrad_props%cldemiss
        Cldrad_props%cmxolw = 0.0
      if (allocated(cld_lsc) .and. allocated(cld_cell)) then
!       crndlw = cld_lsc + cld_cell + cld_meso
        Cldrad_props%crndlw = cld_lsc + cld_cell + cld_meso
      else if (allocated(cld_lsc)) then
!        crndlw = cld_lsc
         Cldrad_props%crndlw = cld_lsc
     else if (allocated(cld_cell)) then
!       crndlw = cld_cell + cld_meso
        Cldrad_props%crndlw = cld_cell + cld_meso
        endif
!  don't allow cloud amounts to exceed one. the normalization
!  doesn't affect cloud properties
!       crndlw = MIN (crndlw, 1.00)
        Cldrad_props%crndlw = MIN (CLdrad_props%crndlw, 1.00)
 
!       camtsw = cmxolw + crndlw
        Cldrad_props%camtsw = Cldrad_props%cmxolw + Cldrad_props%crndlw
    endif   ! do-strat or do-donner-deep
    endif  ! do_esfsw

! DEFINE Cldrad_props% components here !!!

!     Cldrad_props%camtsw = camtsw
!     Cldrad_props%cmxolw = cmxolw
!     Cldrad_props%crndlw = crndlw
!     Cldrad_props%ncldsw = ncldsw
!     Cldrad_props%nmxolw = nmxolw
!     Cldrad_props%nrndlw = nrndlw
!     Cldrad_props%emmxolw = emmxolw
!     Cldrad_props%emrndlw = emrndlw
!     if (allocated (cldext)) then
!       Cldrad_props%cldext = cldext  
!       Cldrad_props%cldsct = cldsct  
!       Cldrad_props%cldasymm = cldasymm
!     endif        
!     if (allocated (cirabsw)) then
!       Cldrad_props%cirrfsw = cirrfsw
!       Cldrad_props%cirabsw = cirabsw
!       Cldrad_props%cvisrfsw = cvisrfsw
!     endif        


!-------------------------------------------------------------------
!    send desired variables to archive package module to hold until
!    output.
!---------------------------------------------------------------------

!     call hold_cloud ( Cldrad_props%cmxolw,  Cldrad_props%crndlw)

!-------------------------------------------------------------------
!     generate netcdf file output fields
!-------------------------------------------------------------------
!     if (Environment%running_gcm .and. Environment%running_fms .and. &
      if (Environment%running_gcm .and.                               &
                            .not. do_no_clouds) then
!        call cloudrad_netcdf (is, js, Time_next, pflux)
        call cloudrad_netcdf (is, js, Time_next, pflux, Cldrad_props)
      endif  

!!! ADD THIS FOR NOW UNTIL THESE ARRAYS MADE LOCAL
      call cloudrad_package_dealloc (Cldrad_props)
!---------------------------------------------------------------------


end subroutine clouddrvr     



!###################################################################

subroutine cloudrad_package_alloc (ix_in, jx_in, kx_in, Cldrad_props)

integer, intent(in) :: ix_in, jx_in, kx_in
type(cldrad_properties_type), intent(inout) :: Cldrad_props

    integer :: unit,idim,jdim,kdim

     ix =ix_in
     jx = jx_in
     kx = kx_in

!    Cldrad_props%NLWCLDB = NLWCLDB

    allocate (   Cldrad_props%camtsw    (ix, jx, kx               ))
    allocate (   Cldrad_props%cmxolw    (ix, jx, kx               ))
    allocate (   Cldrad_props%crndlw    (ix, jx, kx              ))
    allocate(Cldrad_props%emmxolw (ix, jx, kx,           NLWCLDB))
    allocate(Cldrad_props%emrndlw (ix, jx, kx,           NLWCLDB))

    allocate (   Cldrad_props%ncldsw    (ix, jx                    ))
    allocate (   Cldrad_props%nmxolw    (ix, jx              ))
    allocate(    Cldrad_props%nrndlw    (ix, jx                 ))

!   if (trim(swform) /= 'esfsw99' ) then
    if (.not. do_esfsw            ) then

      allocate(Cldrad_props%cirabsw (ix, jx, kx,           nsolwg ))
      allocate(Cldrad_props%cirrfsw (ix, jx, kx,          nsolwg))
      allocate(Cldrad_props%cvisrfsw(ix, jx, kx,            nsolwg ))
      allocate (Cldrad_props%cldext(ix, jx, kx,           0) )
      allocate (Cldrad_props%cldsct(ix, jx, kx,           0) )
      allocate (Cldrad_props%cldasymm(ix, jx, kx,         0) )

    else

      allocate (Cldrad_props%cldext(ix, jx, kx,           nbands) )
      allocate (Cldrad_props%cldsct(ix, jx, kx,           nbands) )
      allocate (Cldrad_props%cldasymm(ix, jx, kx,           nbands) )
      allocate (Cldrad_props%abscoeff(ix, jx, kx,       nlwcldb   ) )
      allocate (Cldrad_props%cldemiss(ix, jx, kx,       nlwcldb   ) )

    endif

!    print *, 'before alloc 1', get_my_pe()
if (do_strat_clouds) then
!    print *, 'before alloc 2', get_my_pe()
 allocate (cld_lsc(ix,jx,kx                            ))
 allocate (cldext_lsc(ix,jx,kx,                            nbands) )
 allocate (cldsct_lsc(ix,jx,kx,                            nbands) )
 allocate (cldasymm_lsc(ix,jx,kx,                            nbands) )
 allocate (abscoeff_lsc (ix,jx,kx,                              NLWCLDB))
  endif
    if (do_donner_deep_clouds) then
 allocate (cld_cell(ix,jx,kx               ) )
 allocate (cldext_cell(ix,jx,kx,                          nbands) )
 allocate (cldsct_cell(ix,jx,kx,                            nbands) )
 allocate (cldasymm_cell(ix,jx,kx,                            nbands) )
 allocate (abscoeff_cell (ix,jx,kx,                           NLWCLDB))
 allocate (cld_meso(ix,jx,kx ) )
 allocate (cldext_meso(ix,jx,kx,                            nbands) )
 allocate (cldsct_meso(ix,jx,kx,                            nbands) )
 allocate (cldasymm_meso(ix,jx,kx,                            nbands) )
 allocate (abscoeff_meso (ix,jx,kx,                       NLWCLDB))
    endif


end subroutine cloudrad_package_alloc



!####################################################################
subroutine cloudrad_package_alloc_rad (ix_in, jx_in, kx_in, Cldrad_props)

integer, intent(in) :: ix_in, jx_in, kx_in
type(cldrad_properties_type), intent(inout) :: Cldrad_props

    integer :: unit,idim,jdim,kdim
     ix =ix_in
     jx = jx_in
     kx = kx_in

    allocate (   Cldrad_props%camtsw    (ix, jx, kx               ))
    allocate (   Cldrad_props%cmxolw    (ix, jx, kx               ))
    allocate (   Cldrad_props%crndlw    (ix, jx, kx              ))
    allocate(Cldrad_props%emmxolw (ix, jx, kx,           NLWCLDB))
    allocate(Cldrad_props%emrndlw (ix, jx, kx,           NLWCLDB))

    allocate (   Cldrad_props%ncldsw    (ix, jx                    ))
    allocate (   Cldrad_props%nmxolw    (ix, jx              ))
    allocate(    Cldrad_props%nrndlw    (ix, jx                 ))

!   if (trim(swform) /= 'esfsw99' ) then
    if (.not. do_esfsw            ) then

      allocate(Cldrad_props%cirabsw (ix, jx, kx,           nsolwg ))
      allocate(Cldrad_props%cirrfsw (ix, jx, kx,          nsolwg))
      allocate(Cldrad_props%cvisrfsw(ix, jx, kx,            nsolwg ))
      allocate (Cldrad_props%cldext(ix, jx, kx,           0) )
      allocate (Cldrad_props%cldsct(ix, jx, kx,           0) )
      allocate (Cldrad_props%cldasymm(ix, jx, kx,         0) )

    else

      allocate (Cldrad_props%cldext(ix, jx, kx,           nbands) )
      allocate (Cldrad_props%cldsct(ix, jx, kx,           nbands) )
      allocate (Cldrad_props%cldasymm(ix, jx, kx,           nbands) )
      allocate (Cldrad_props%abscoeff(ix, jx, kx,       nlwcldb   ) )
      allocate (Cldrad_props%cldemiss(ix, jx, kx,       nlwcldb   ) )

    endif

if (do_strat_clouds) then
 allocate (cld_lsc(ix,jx,kx                            ))
 allocate (cldext_lsc(ix,jx,kx,                            nbands) )
 allocate (cldsct_lsc(ix,jx,kx,                            nbands) )
 allocate (cldasymm_lsc(ix,jx,kx,                            nbands) )
 allocate (abscoeff_lsc (ix,jx,kx,                              NLWCLDB))
  endif
    if (do_donner_deep_clouds) then
 allocate (cld_cell(ix,jx,kx               ) )
 allocate (cldext_cell(ix,jx,kx,                          nbands) )
 allocate (cldsct_cell(ix,jx,kx,                            nbands) )
 allocate (cldasymm_cell(ix,jx,kx,                            nbands) )
 allocate (abscoeff_cell (ix,jx,kx,                           NLWCLDB))
 allocate (cld_meso(ix,jx,kx ) )
 allocate (cldext_meso(ix,jx,kx,                            nbands) )
 allocate (cldsct_meso(ix,jx,kx,                            nbands) )
 allocate (cldasymm_meso(ix,jx,kx,                            nbands) )
 allocate (abscoeff_meso (ix,jx,kx,                       NLWCLDB))
    endif
end subroutine cloudrad_package_alloc_rad



!####################################################################

subroutine cloudrad_package_end

    return


end subroutine cloudrad_package_end



!####################################################################

subroutine cloudrad_package_dealloc (Cldrad_props)
       
type(cldrad_properties_type), intent(in) :: Cldrad_props

!   deallocate (camtsw    )
!   deallocate (cmxolw    )
!   deallocate (crndlw    )
!   deallocate (emmxolw   )
!   deallocate (emrndlw   )
!   deallocate (ncldsw    )
!   deallocate (nmxolw    )
!   deallocate (nrndlw    )
!   deallocate (Cldrad_props%camtsw    )
!   deallocate (Cldrad_props%cmxolw    )
!   deallocate (Cldrad_props%crndlw    )
!   deallocate (Cldrad_props%emmxolw   )
!   deallocate (Cldrad_props%emrndlw   )
!   deallocate (Cldrad_props%ncldsw    )
!   deallocate (Cldrad_props%nmxolw    )
!   deallocate (Cldrad_props%nrndlw    )

!   if (trim(swform) /= 'esfsw99' ) then
!   if (.not. do_esfsw            ) then
!     deallocate (cirabsw   )
!     deallocate (cirrfsw   )
!     deallocate (cvisrfsw  )
!     deallocate (Cldrad_props%cirabsw   )
!     deallocate (Cldrad_props%cirrfsw   )
!     deallocate (Cldrad_props%cvisrfsw  )
!     deallocate (Cldrad_props%cldasymm) 
!     deallocate (Cldrad_props%cldsct  ) 
!     deallocate (Cldrad_props%cldext  ) 
!   else
!     deallocate (cldasymm) 
!     deallocate (cldsct  ) 
!     deallocate (cldext  ) 
!     deallocate (Cldrad_props%cldasymm) 
!     deallocate (Cldrad_props%cldsct  ) 
!     deallocate (Cldrad_props%cldext  ) 
!     deallocate (Cldrad_props%abscoeff) 
!     deallocate (Cldrad_props%cldemiss) 
!   endif


!    print *, 'before dealloc 1', get_my_pe()
 if (do_strat_clouds) then
!    print *, 'before dealloc 2', get_my_pe()
      deallocate (cld_lsc)
      deallocate (cldext_lsc)
      deallocate (cldsct_lsc)
      deallocate (cldasymm_lsc)
     deallocate (abscoeff_lsc)
  endif
  if (do_donner_deep_clouds) then
      deallocate (cld_cell)
       deallocate (cldext_cell)
        deallocate (cldsct_cell)
    deallocate (cldasymm_cell)
    deallocate (abscoeff_cell)
       deallocate (cld_meso)
     deallocate (cldext_meso)
         deallocate (cldsct_meso)
        deallocate (cldasymm_meso)
        deallocate (abscoeff_meso)
  endif


end subroutine cloudrad_package_dealloc



!####################################################################

!subroutine get_clouds_for_lhsw (is, ie, js, je,  Cldrad_props, &
subroutine convert_to_cloud_space (is, ie, js, je,  Cldrad_props, &
                               cirabswkc, cirrfswkc,  &
				cvisrfswkc, ktopswkc, kbtmswkc,  &
				camtswkc, Cldspace_rad)
      
integer, intent(in) :: is, ie, js, je
type(cldrad_properties_type), intent(in) :: Cldrad_props
real, dimension(:,:,:),   intent(out) :: camtswkc
real, dimension(:,:,:  ), intent(out) :: cirabswkc, cirrfswkc,  &
					 cvisrfswkc
integer, dimension(:,:,:),intent(out) :: ktopswkc, kbtmswkc
type(cld_space_properties_type), intent(inout) :: Cldspace_rad
!-------------------------------------------------------------------

      logical :: do_donner_ice
      integer :: i,j,k, index_kbot, index_ktop, ngp
      integer :: kcldsmx

!---------------------------------------------------------------------
!     obtain arrays for shortwave cloud properties in "cloud space".
!     "bottom up" counting is used in conformity with the lhsw
!     radiative algorithm. (ie, index = 1 = lowest cloud (if any), etc.)
!---------------------------------------------------------------------
      
      kcldsmx = MAXVAL(Cldrad_props%ncldsw)
      if (kcldsmx /= 0.0) then
      camtswkc(:,:,:) = 0.0E+00
      cirabswkc(:,:,:  ) = 0.0E+00
      cirrfswkc(:,:,:  ) = 0.0E+00
      cvisrfswkc(:,:,:) = 0.0E+00
!     ktopswkc(:,:,:) = KSRAD
!     kbtmswkc(:,:,:) = KSRAD
      ktopswkc(:,:,:) = 1        
      kbtmswkc(:,:,:) = 1        

!     do j=JSRAD,JERAD
!       do i=ISRAD,IERAD
      do j=1,jx         
        do i=1,ix         
	  index_kbot = 0
	  index_ktop = 0
!---------------------------------------------------------------------
!     in the kx   'th layer, the presence of cloud (camtsw) implies a
!     shortwave cloud bottom at level (kx   +1) but nothing about
!     cloud tops.
!---------------------------------------------------------------------
!         if (camtsw(i,j,KERAD) .GT. 0.0E+00) then
!         if (camtsw(i,j,kx   ) .GT. 0.0E+00) then
          if (Cldrad_props%camtsw(i,j,kx   ) .GT. 0.0E+00) then
            index_kbot = index_kbot + 1
!    kbtmswkc(i,j,index_kbot) = KERAD+1
!    camtswkc(i,j,index_kbot) = camtsw(i,j,KERAD)
	    kbtmswkc(i,j,index_kbot) = kx   +1
!    camtswkc(i,j,index_kbot) = camtsw(i,j,kx   )
	    camtswkc(i,j,index_kbot) = Cldrad_props%camtsw(i,j,kx   )
	    do ngp=1,nsolwg  
!      cirabswkc(i,j,index_kbot,ngp) = cirabsw(i,j,KERAD,ngp)
!             cirrfswkc(i,j,index_kbot,ngp) = cirrfsw(i,j,KERAD,ngp)
!      cvisrfswkc(i,j,index_kbot,ngp) = cvisrfsw(i,j,KERAD,ngp)
!!     cirabswkc(i,j,index_kbot,ngp) = cirabsw(i,j,kx   ,ngp)
!!            cirrfswkc(i,j,index_kbot,ngp) = cirrfsw(i,j,kx   ,ngp)
!!     cvisrfswkc(i,j,index_kbot,ngp) = cvisrfsw(i,j,kx   ,ngp)
	      cirabswkc(i,j,index_kbot    ) = Cldrad_props%cirabsw(i,j,kx   ,ngp)
              cirrfswkc(i,j,index_kbot    ) = Cldrad_props%cirrfsw(i,j,kx   ,ngp)
	      cvisrfswkc(i,j,index_kbot    ) = Cldrad_props%cvisrfsw(i,j,kx   ,ngp)
	    end do
          endif
!---------------------------------------------------------------------
!   in other layers, cloud bottoms and tops are determined according
!   to changes in shortwave cloud (and special case).
!---------------------------------------------------------------------
!         do k=KERAD-1,KSRAD,-1
          do k=kx-1,1    ,-1
!---------------------------------------------------------------------
!     cloud bottoms.
!--------------------------------------------------------------------
!    if (camtsw(i,j,k) .NE. camtsw(i,j,k+1) .AND.    &
!     		camtsw(i,j,k) .GT. 0.0E+00              .OR.  &
	    if (Cldrad_props%camtsw(i,j,k) .NE. Cldrad_props%camtsw(i,j,k+1) .AND.    &
      		Cldrad_props%camtsw(i,j,k) .GT. 0.0E+00              .OR.  &
!---------------------------------------------------------------------
!     special case where shortwave cloud amounts for adjacent
!     layers are equal, but a randomly overlapped cloud exists
!     in at least one of these layers
!---------------------------------------------------------------------
!               camtsw(i,j,k) .EQ. camtsw(i,j,k+1) .AND.   &
!              (crndlw(i,j,k) .NE. 0.0E+00 .OR.      &
!               crndlw(i,j,k+1) .NE. 0.0E+00)                ) then
                Cldrad_props%camtsw(i,j,k) .EQ. Cldrad_props%camtsw(i,j,k+1) .AND.   &
               (Cldrad_props%crndlw(i,j,k) .NE. 0.0E+00 .OR.      &
                Cldrad_props%crndlw(i,j,k+1) .NE. 0.0E+00)                ) then
              index_kbot = index_kbot + 1
	      kbtmswkc(i,j,index_kbot) = k+1
!      camtswkc(i,j,index_kbot) = camtsw(i,j,k)
	      camtswkc(i,j,index_kbot) = Cldrad_props%camtsw(i,j,k)
	      do ngp=1,nsolwg  
!      cirabswkc(i,j,index_kbot,ngp) = cirabsw(i,j,k,ngp)
!      cirrfswkc(i,j,index_kbot,ngp) = cirrfsw(i,j,k,ngp)
!      cvisrfswkc(i,j,index_kbot,ngp) = cvisrfsw(i,j,k,ngp)
	      cirabswkc(i,j,index_kbot    ) = Cldrad_props%cirabsw(i,j,k,ngp)
	      cirrfswkc(i,j,index_kbot    ) = Cldrad_props%cirrfsw(i,j,k,ngp)
	      cvisrfswkc(i,j,index_kbot    ) = Cldrad_props%cvisrfsw(i,j,k,ngp)
	      end do
            endif
!---------------------------------------------------------------------
!     cloud tops.
!---------------------------------------------------------------------
!           if (camtsw(i,j,k) .NE. camtsw(i,j,k+1) .AND.      &
!               camtsw(i,j,k+1) .GT. 0.0E+00            .OR.  &
            if (Cldrad_props%camtsw(i,j,k) .NE. Cldrad_props%camtsw(i,j,k+1) .AND.      &
                Cldrad_props%camtsw(i,j,k+1) .GT. 0.0E+00            .OR.  &
!---------------------------------------------------------------------
!     special case where shortwave cloud amounts for adjacent
!     layers are equal, but a randomly overlapped cloud exists
!     in at least one of these layers
!---------------------------------------------------------------------
!               camtsw(i,j,k) .EQ. camtsw(i,j,k+1) .AND.    &
!              (crndlw(i,j,k) .NE. 0.0E+00 .OR.        &
!               crndlw(i,j,k+1) .NE. 0.0E+00)                ) then
                Cldrad_props%camtsw(i,j,k) .EQ. Cldrad_props%camtsw(i,j,k+1) .AND.    &
               (Cldrad_props%crndlw(i,j,k) .NE. 0.0E+00 .OR.        &
                Cldrad_props%crndlw(i,j,k+1) .NE. 0.0E+00)                ) then
	      index_ktop = index_ktop + 1
	      ktopswkc(i,j,index_ktop) = k+1
            endif
          enddo
        enddo
      enddo

!-------------------------------------------------------------------
!     if donner iceclouds are active, send some info to that module.
!---------------------------------------------------------------------

!!! 2/05/02
!!!! REMOVE THIS CALL NOW RATHER THAN CHANGE RANK OF THE swkc ARRAYS
!!  IT WILL ULTIMATELY BE REMOVED ANYWAY 
!     if (Environment%running_gcm) then
!       if (Environment%using_sky_periphs) then
!         call inquire_donner_ice (do_donner_ice)
!         if (do_donner_ice) then
!           call sw_albedo_zen_angle (Cldrad_props%ncldsw,   &
!                                   cvisrfswkc, cirrfswkc,  &
!                                    nsolwg)
!         endif
!       endif
!     endif 

!-------------------------------------------------------------------
!     send needed lhsw cloud information to radiation diagnostics mod.
!---------------------------------------------------------------------

!     do j=js,je        
!       if (Rad_control%do_raddg(j)) then
!         call radiag_from_clouds_lh(   camtswkc, cirabswkc, cirrfswkc,&
!             Cldrad_props%cmxolw, Cldrad_props%crndlw, cvisrfswkc,   &
! Cldrad_props%emmxolw, Cldrad_props%emrndlw,     Cldrad_props%ncldsw ,  &
!        Cldrad_props%nmxolw , &
!	         Cldrad_props%nrndlw , ktopswkc, kbtmswkc, &
!			     j, j-js+1, is, ie) 
!       endif  
!     end do
!---------------------------------------------------------------------
      allocate ( Cldspace_rad%camtswkc(ie-is+1, je-js+1,  &
                                 size(camtswkc,3) ) )
      allocate ( Cldspace_rad%cirabswkc(ie-is+1, je-js+1,  &
                                 size(camtswkc,3)   ) )
      allocate ( Cldspace_rad%cirrfswkc(ie-is+1, je-js+1,  &
                                 size(camtswkc,3)   ) )
      allocate ( Cldspace_rad%cvisrfswkc(ie-is+1, je-js+1,  &
                                 size(camtswkc,3)   ) )
      allocate ( Cldspace_rad%ktopswkc(ie-is+1, je-js+1,  &
                                 size(camtswkc,3) ) )
      allocate ( Cldspace_rad%kbtmswkc(ie-is+1, je-js+1,  &
                                 size(camtswkc,3)+1 ) )
      Cldspace_rad%camtswkc = camtswkc
      Cldspace_rad%cirabswkc = cirabswkc
      Cldspace_rad%cirrfswkc = cirrfswkc
      Cldspace_rad%cvisrfswkc = cvisrfswkc
      Cldspace_rad%ktopswkc = ktopswkc
      Cldspace_rad%kbtmswkc = kbtmswkc
      else
      endif



!end subroutine get_clouds_for_lhsw
end subroutine convert_to_cloud_space



!####################################################################

!subroutine get_clouds_for_esfsw (is, ie, js, je, camtsw_out)
      
!--------------------------------------------------------------------
!integer, intent(in) :: is, ie, js, je
!real, dimension(:,:,:), intent(out)   :: camtsw_out
!-------------------------------------------------------------------

!      integer           ::  j
!-------------------------------------------------------------------

!      camtsw_out(:,:,:) = camtsw(:,:,:)

!     do j=js,je           
!       if (Rad_control%do_raddg(j)) then
!         call radiag_from_clouds_esf(   camtsw_out, cmxolw, crndlw, &
!                                     emmxolw, emrndlw, ncldsw, nmxolw,&
!			      nrndlw, j, j-js+1, is, ie)
!       endif  
!     end do


!end subroutine get_clouds_for_esfsw


!#####################################################################

!subroutine get_ncldsw (ncldsw_out)

!integer, dimension(:,:), intent(out) :: ncldsw_out


!     ncldsw_out(:,:) = ncldsw(:,:)


!end subroutine get_ncldsw


!###################################################################

!subroutine get_clouds_for_lwrad (cmxolw_out, crndlw_out, emmxolw_out, &
!				 emrndlw_out, nmxolw_out, nrndlw_out)

!real, dimension(:,:,:),   intent(out) :: cmxolw_out, crndlw_out
!real, dimension(:,:,:,:), intent(out) :: emmxolw_out, emrndlw_out
!integer, dimension(:,:),  intent(out) :: nmxolw_out, nrndlw_out


!     cmxolw_out(:,:,:)    = cmxolw(:,:,:)
!     crndlw_out(:,:,:)    = crndlw(:,:,:)
!     emmxolw_out(:,:,:,:) = emmxolw(:,:,:,:)
!     emrndlw_out(:,:,:,:) = emrndlw(:,:,:,:)
!     nmxolw_out(:,:)      = nmxolw(:,:)
!     nrndlw_out(:,:)      = nrndlw(:,:)


!end subroutine get_clouds_for_lwrad



!#################################################################### 

!subroutine get_swcldprops_from_cloudrad (n, cldext_out, cldsct_out, &
!					 cldasymm_out)
 
!---------------------------------------------------------------------
!integer,                intent(in)    ::  n
!real, dimension(:,:,:), intent(out)   ::  cldext_out, cldsct_out,   &
!					  cldasymm_out
!---------------------------------------------------------------------

!     integer                     :: i,j,k
!---------------------------------------------------------------------

!     do k=KSRAD,KERAD
!       do j=JSRAD,JERAD
!         do i=ISRAD, IERAD
!           cldext_out   (i,j,k) =   cldext   (i,j,k,n)
!           cldsct_out   (i,j,k) =   cldsct   (i,j,k,n)
!           cldasymm_out (i,j,k) =   cldasymm (i,j,k,n)
!         enddo
!       enddo
!     enddo

!           cldext_out   (:,:,:) =   cldext   (:,:,:,n)
!           cldsct_out   (:,:,:) =   cldsct   (:,:,:,n)
!           cldasymm_out (:,:,:) =   cldasymm (:,:,:,n)
  
!end subroutine get_swcldprops_from_cloudrad 

!####################################################################

subroutine default_clouds (Cldrad_props)

type(cldrad_properties_type), intent(inout) :: Cldrad_props
      
      Cldrad_props%cmxolw(:,:,:) = 0.0E+00
      Cldrad_props%crndlw(:,:,:) = 0.0E+00
       Cldrad_props%camtsw(:,:,:) = 0.0E+00
       Cldrad_props%emmxolw(:,:,:,:) = 1.0E+00
       Cldrad_props%emrndlw(:,:,:,:) = 1.0E+00
       Cldrad_props%nmxolw (:,:) = 0
       Cldrad_props%nrndlw (:,:) = 0
       Cldrad_props%ncldsw (:,:) = 0

!     if (allocated ( Cldrad_props%cirrfsw) ) then
      if (associated ( Cldrad_props%cirrfsw) ) then
         Cldrad_props%cirrfsw (:,:,:,:) = 0.0E+00
         Cldrad_props%cvisrfsw(:,:,:,:) = 0.0E+00
         Cldrad_props%cirabsw (:,:,:,:) = 0.0E+00
         Cldrad_props%cldsct  (:,:,:,:) = 0.0E+00
         Cldrad_props%cldext  (:,:,:,:) = 0.0E+00
         Cldrad_props%cldasymm(:,:,:,:) = 0.0E+00
!     else if (allocated ( Cldrad_props%cldext)  ) then
      else if (associated ( Cldrad_props%cldext)  ) then
         Cldrad_props%cldsct  (:,:,:,:) = 0.0E+00
         Cldrad_props%cldext  (:,:,:,:) = 0.0E+00
         Cldrad_props%cldasymm(:,:,:,:) = 0.0E+00
      endif


end subroutine default_clouds



!####################################################################

subroutine cloudrad_netcdf(is, js, Time_diag, pflux, Cldrad_props,mask)

integer,                 intent(in)           ::  is, js
real, dimension(:,:,:),  intent(in)           ::  pflux
real, dimension(:,:,:),  intent(in), optional ::  mask
type(cldrad_properties_type), intent(in) :: Cldrad_props
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

!     allocate (cloud(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
      allocate (cloud(ix, jx, kx                           ) )

!------- TOTAL CLOUD AMOUNT  -----------------
      if ( id_tot_cld_amt > 0 ) then
!         allocate (tca (ISRAD:IERAD, JSRAD:JERAD) )
          allocate (tca (ix, jx                  ) )
          tca = 1.0
!  do k=KSRAD,KERAD
	  do k=1,kx        
	    tca(:,:) = tca(:,:)*(1.0-Cldrad_props%camtsw(:,:,k))
          end do
	  tca = 1.-tca
         tca = tca*100.
         used = send_data ( id_tot_cld_amt, tca, Time_diag, is, js)
         deallocate (tca)
      endif

!---- high,mid,low cloud diagnostics ----
      if ( id_high_cld_amt > 0 .or. id_mid_cld_amt > 0 .or. &
            id_low_cld_amt > 0 ) then
!       allocate (hml_ca(ISRAD:IERAD, JSRAD:JERAD,3))
        allocate (hml_ca(ix, jx,                  3))
        call compute_isccp_clds (pflux, Cldrad_props%camtsw, hml_ca)
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
    used = send_data ( id_cld_amt, Cldrad_props%camtsw, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- lsc cloud amount (only do once ?) -------------------------
 
     if ( id_lsc_cld_amt > 0 ) then
       if (do_strat_clouds) then
     used = send_data ( id_lsc_cld_amt, cld_lsc, Time_diag, is, js, 1, &
                        rmask=mask )
       endif
     endif

     if ( id_lsc_cld_ext_uv > 0 ) then
       if (do_strat_clouds) then
       cloud(:,:,:) = cldext_lsc(:,:,:,22)
  used = send_data ( id_lsc_cld_ext_uv, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
      endif
     endif
 
     if ( id_lsc_cld_ext_vis > 0 ) then
       if (do_strat_clouds) then
       cloud(:,:,:) = cldext_lsc(:,:,:,12)
   used = send_data ( id_lsc_cld_ext_vis, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
       endif
     endif

     if ( id_lsc_cld_ext_nir > 0 ) then
      if (do_strat_clouds) then
     cloud(:,:,:) = cldext_lsc(:,:,:,8)
   used = send_data ( id_lsc_cld_ext_nir, cloud, Time_diag, is, js, 1, &
                         rmask=mask )
      endif
    endif

    if ( id_lsc_cld_sct_uv > 0 ) then
      if (do_strat_clouds) then
      cloud(:,:,:) = cldsct_lsc(:,:,:,22)
   used = send_data ( id_lsc_cld_sct_uv, cloud, Time_diag, is, js, 1, &
                             rmask=mask )
     endif
   endif

    if ( id_lsc_cld_sct_vis > 0 ) then
       if (do_strat_clouds) then
         cloud(:,:,:) = cldsct_lsc(:,:,:,12)
   used = send_data ( id_lsc_cld_sct_vis, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
       endif
    endif

   if ( id_lsc_cld_sct_nir > 0 ) then
      if (do_strat_clouds) then
       cloud(:,:,:) = cldsct_lsc(:,:,:,8)
  used = send_data ( id_lsc_cld_sct_nir, cloud, Time_diag, is, js, 1, &
                           rmask=mask )
      endif
     endif
 
    if ( id_lsc_cld_asymm_uv > 0 ) then
      if (do_strat_clouds) then
       cloud(:,:,:) = cldasymm_lsc(:,:,:,22)
used = send_data ( id_lsc_cld_asymm_uv, cloud, Time_diag, is, js, 1, &
                      rmask=mask )
      endif
    endif

    if ( id_lsc_cld_asymm_vis > 0 ) then
      if (do_strat_clouds) then
       cloud(:,:,:) = cldasymm_lsc(:,:,:,12)
 used = send_data ( id_lsc_cld_asymm_vis, cloud, Time_diag, is, js, 1, &
                     rmask=mask )
        endif
      endif

    if ( id_lsc_cld_asymm_nir > 0 ) then
      if (do_strat_clouds) then
     cloud(:,:,:) = cldasymm_lsc(:,:,:,8)
 used = send_data ( id_lsc_cld_asymm_nir, cloud, Time_diag, is, js, 1, &
                        rmask=mask )
       endif
     endif

!------- cell cloud amount (only do once ?) -------------------------
     if ( id_cell_cld_amt > 0 ) then
    if (do_donner_deep_clouds) then
  used = send_data ( id_cell_cld_amt, cld_cell, Time_diag, is, js, 1,&
                        rmask=mask )
    endif
    endif

    if ( id_cell_cld_ext_uv > 0 ) then
      if (do_donner_deep_clouds) then
      cloud(:,:,:) = cldext_cell(:,:,:,22)
  used = send_data ( id_cell_cld_ext_uv, cloud, Time_diag, is, js, 1, &
                        rmask=mask )
     endif
  endif

    if ( id_cell_cld_ext_vis > 0 ) then
       if (do_donner_deep_clouds) then
       cloud(:,:,:) = cldext_cell(:,:,:,12)
  used = send_data ( id_cell_cld_ext_vis, cloud, Time_diag, is, js, 1, &
                             rmask=mask )
      endif
    endif
 
   if ( id_cell_cld_ext_nir > 0 ) then
     if (do_donner_deep_clouds) then
       cloud(:,:,:) = cldext_cell(:,:,:,8)
  used = send_data ( id_cell_cld_ext_nir, cloud, Time_diag, is, js, 1, &
                          rmask=mask )
      endif
    endif

   if ( id_cell_cld_sct_uv > 0 ) then
     if (do_donner_deep_clouds) then
      cloud(:,:,:) = cldsct_cell(:,:,:,22)
  used = send_data ( id_cell_cld_sct_uv, cloud, Time_diag, is, js, 1, &
                             rmask=mask )
      endif
  endif
 
    if ( id_cell_cld_sct_vis > 0 ) then
      if (do_donner_deep_clouds) then
       cloud(:,:,:) = cldsct_cell(:,:,:,12)
  used = send_data ( id_cell_cld_sct_vis, cloud, Time_diag, is, js, 1, &
                       rmask=mask )
      endif
    endif

    if ( id_cell_cld_sct_nir > 0 ) then
      if (do_donner_deep_clouds) then
       cloud(:,:,:) = cldsct_cell(:,:,:,8)
 used = send_data ( id_cell_cld_sct_nir, cloud, Time_diag, is, js, 1, &
                           rmask=mask )
      endif
    endif

    if ( id_cell_cld_asymm_uv > 0 ) then
      if (do_donner_deep_clouds) then
       cloud(:,:,:) = cldasymm_cell(:,:,:,22)
 used = send_data ( id_cell_cld_asymm_uv, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
      endif
    endif

    if ( id_cell_cld_asymm_vis > 0 ) then
       if (do_donner_deep_clouds) then
      cloud(:,:,:) = cldasymm_cell(:,:,:,12)
used = send_data ( id_cell_cld_asymm_vis, cloud, Time_diag, is, js, 1, &
                      rmask=mask )
    endif
    endif

     if ( id_cell_cld_asymm_nir > 0 ) then
     if (do_donner_deep_clouds) then
      cloud(:,:,:) = cldasymm_cell(:,:,:,8)
used = send_data ( id_cell_cld_asymm_nir, cloud, Time_diag, is, js, 1, &
                             rmask=mask )
  endif
      endif
 
!------- meso cloud amount (only do once ?) -------------------------
     if ( id_meso_cld_amt > 0 ) then
      if (do_donner_deep_clouds) then
   used = send_data ( id_meso_cld_amt, cld_meso, Time_diag, is, js, 1,&
                            rmask=mask )
    endif
    endif

     if ( id_meso_cld_ext_uv > 0 ) then
       if (do_donner_deep_clouds) then
       cloud(:,:,:) = cldext_meso(:,:,:,22)
   used = send_data ( id_meso_cld_ext_uv, cloud, Time_diag, is, js, 1, &
                     rmask=mask )
      endif
      endif

      if ( id_meso_cld_ext_vis > 0 ) then
      if (do_donner_deep_clouds) then
      cloud(:,:,:) = cldext_meso(:,:,:,12)
 used = send_data ( id_meso_cld_ext_vis, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
       endif
       endif

       if ( id_meso_cld_ext_nir > 0 ) then
       if (do_donner_deep_clouds) then
         cloud(:,:,:) = cldext_meso(:,:,:,8)
 used = send_data ( id_meso_cld_ext_nir, cloud, Time_diag, is, js, 1, &
                          rmask=mask )
      endif
     endif
 
       if ( id_meso_cld_sct_uv > 0 ) then
       if (do_donner_deep_clouds) then
       cloud(:,:,:) = cldsct_meso(:,:,:,22)
   used = send_data ( id_meso_cld_sct_uv, cloud, Time_diag, is, js, 1, &
                             rmask=mask )
       endif
       endif
 
      if ( id_meso_cld_sct_vis > 0 ) then
       if (do_donner_deep_clouds) then
       cloud(:,:,:) = cldsct_meso(:,:,:,12)
  used = send_data ( id_meso_cld_sct_vis, cloud, Time_diag, is, js, 1, &
                             rmask=mask )
       endif
        endif

      if ( id_meso_cld_sct_nir > 0 ) then
      if (do_donner_deep_clouds) then
        cloud(:,:,:) = cldsct_meso(:,:,:,8)
  used = send_data ( id_meso_cld_sct_nir, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
       endif
       endif

       if ( id_meso_cld_asymm_uv > 0 ) then
       if (do_donner_deep_clouds) then
      cloud(:,:,:) = cldasymm_meso(:,:,:,22)
 used = send_data ( id_meso_cld_asymm_uv, cloud, Time_diag, is, js, 1, &
                             rmask=mask )
       endif
       endif

       if ( id_meso_cld_asymm_vis > 0 ) then
       if (do_donner_deep_clouds) then
       cloud(:,:,:) = cldasymm_meso(:,:,:,12)
used = send_data ( id_meso_cld_asymm_vis, cloud, Time_diag, is, js, 1, &
                            rmask=mask )
      endif
      endif

       if ( id_meso_cld_asymm_nir > 0 ) then
       if (do_donner_deep_clouds) then
       cloud(:,:,:) = cldasymm_meso(:,:,:,8)
used = send_data ( id_meso_cld_asymm_nir, cloud, Time_diag, is, js, 1, &
                         rmask=mask )
       endif
     endif




!------- cloud emissivity ---------------------------------------

          if ( id_em_cld_10u > 0  ) then
          if (Lw_control%do_lwcldemiss) then
!    output cloud emissivity is  weighted average of the random and max
!    overlap emissivities, over the 990-1070 cm-1 band (band 5 of 7)
            cloud(:,:,:) = (Cldrad_props%crndlw(:,:,:)*Cldrad_props%emrndlw(:,:,:,5) +    &
                           Cldrad_props%cmxolw(:,:,:)*Cldrad_props%emmxolw(:,:,:,5))/   &
                            (Cldrad_props%crndlw(:,:,:) + Cldrad_props%cmxolw(:,:,:) + 1.0E-10)
      used = send_data ( id_em_cld_10u, cloud, Time_diag, is, js, 1, &
                rmask=mask )
     endif
    endif

  if (id_em_cld_lw > 0   ) then
     if (Lw_control%do_lwcldemiss) then
     else
!    output cloud emissivity is weighted average of the random and max
!    overlap emissivities, over 1 band (0-2200 cm-1)
       cloud(:,:,:) = (Cldrad_props%crndlw(:,:,:)*Cldrad_props%emrndlw(:,:,:,1) +    &
                      Cldrad_props%cmxolw(:,:,:)*Cldrad_props%emmxolw(:,:,:,1))/   &
      (Cldrad_props%crndlw(:,:,:) + Cldrad_props%cmxolw(:,:,:) + 1.0E-10)
   used = send_data ( id_em_cld_lw, cloud, Time_diag, is, js, 1, &
               rmask=mask )
  endif
endif

   if ( id_abs_lsc_cld_10u > 0  ) then
     if (Lw_control%do_lwcldemiss) then
      used = send_data ( id_abs_lsc_cld_10u, abscoeff_lsc(:,:,:,5), Time_diag, is, js, 1, &
                rmask=mask )
     endif
   endif

  if ( id_abs_lsc_cld_lw > 0  ) then
    if (Lw_control%do_lwcldemiss) then
   else
    used = send_data ( id_abs_lsc_cld_lw, abscoeff_lsc(:,:,:,1), Time_diag, is, js, 1, &
                 rmask=mask )
   endif
  endif

  if ( id_abs_cell_cld_10u > 0  ) then
    if (Lw_control%do_lwcldemiss) then
     used = send_data ( id_abs_cell_cld_10u, abscoeff_cell(:,:,:,5), Time_diag, is, js, 1, &
            rmask=mask )
  endif
   endif

   if ( id_abs_cell_cld_lw > 0  ) then
    if (Lw_control%do_lwcldemiss) then
    else
     used = send_data ( id_abs_cell_cld_lw, abscoeff_cell(:,:,:,1), Time_diag, is, js, 1, &
               rmask=mask )
    endif
   endif

   if ( id_abs_meso_cld_10u > 0  ) then
     if (Lw_control%do_lwcldemiss) then
  used = send_data ( id_abs_meso_cld_10u, abscoeff_meso(:,:,:,5), Time_diag, is, js, 1, &
                 rmask=mask )
     endif
   endif
 
  if ( id_abs_meso_cld_lw > 0  ) then
    if (Lw_control%do_lwcldemiss) then
    else
     used = send_data ( id_abs_meso_cld_lw, abscoeff_meso(:,:,:,1), Time_diag, is, js, 1, &
           rmask=mask )
    endif
   endif

  if ( id_abs_cld_10u > 0  ) then
    if (Lw_control%do_lwcldemiss) then
 used = send_data ( id_abs_cld_10u, Cldrad_props%abscoeff(:,:,:,5), Time_diag, is, js, 1, &
             rmask=mask )
     endif
  endif

     if ( id_abs_cld_lw > 0  ) then
     if (Lw_control%do_lwcldemiss) then
     else
      used = send_data ( id_abs_cld_lw, Cldrad_props%abscoeff(:,:,:,1), Time_diag, is, js, 1, &
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

!     if (trim(swform) == 'lhsw' ) then  ! diagnostics for lacis-hansen
      if (do_lhsw                ) then  ! diagnostics for lacis-hansen

!------- ultra-violet reflected by cloud -----------------------------
        if ( id_alb_uv_cld > 0  ) then
              cloud(:,:,:) = Cldrad_props%cvisrfsw(:,:,:,1)
          used = send_data ( id_alb_uv_cld, cloud, Time_diag, is, js,  &
			    1, rmask=mask )
        endif

!------- infra-red reflected by cloud -----------------------------
        if ( id_alb_nir_cld > 0  ) then
              cloud(:,:,:) =  Cldrad_props%cirrfsw(:,:,:,1)
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
              cloud(:,:,:) =  Cldrad_props%cirabsw(:,:,:,1)
          used = send_data ( id_abs_nir_cld, cloud, Time_diag, is, js, &
			   1, rmask=mask )
        endif
!  else if (trim(swform) == 'esfsw99') then ! diagnostics for esf sw code
  else if (do_esfsw                  ) then ! diagnostics for esf sw code

     if ( id_ext_cld_uv > 0  ) then
       cloud(:,:,:) =  Cldrad_props%cldext(:,:,:,22)
     used = send_data ( id_ext_cld_uv, cloud, Time_diag, is, js, &
                          1, rmask=mask )
      endif
  
   if ( id_sct_cld_uv > 0  ) then
       cloud(:,:,:) =  Cldrad_props%cldsct(:,:,:,22)
       used = send_data ( id_sct_cld_uv, cloud, Time_diag, is, js, &
                          1, rmask=mask )
    endif

    if ( id_asymm_cld_uv > 0  ) then
         cloud(:,:,:) = 100.0* Cldrad_props%cldasymm(:,:,:,22)
         used = send_data ( id_asymm_cld_uv, cloud, Time_diag, is, js, &
                          1, rmask=mask )
     endif

    if ( id_ext_cld_vis > 0  ) then
      cloud(:,:,:) =  Cldrad_props%cldext(:,:,:,12)
      used = send_data ( id_ext_cld_vis, cloud, Time_diag, is, js, &
                             1, rmask=mask )
    endif

    if ( id_sct_cld_vis > 0  ) then
      cloud(:,:,:) =  Cldrad_props%cldsct(:,:,:,12)
      used = send_data ( id_sct_cld_vis, cloud, Time_diag, is, js, &
                           1, rmask=mask )
    endif

     if ( id_asymm_cld_vis > 0  ) then
        cloud(:,:,:) = 100.0* Cldrad_props%cldasymm(:,:,:,12)
        used = send_data ( id_asymm_cld_vis, cloud, Time_diag, is, js, &
                            1, rmask=mask )
     endif

     if ( id_ext_cld_nir > 0  ) then
      cloud(:,:,:) =  Cldrad_props%cldext(:,:,:,8)
      used = send_data ( id_ext_cld_nir, cloud, Time_diag, is, js, &
                           1, rmask=mask )
     endif

      if ( id_sct_cld_nir > 0  ) then
        cloud(:,:,:) =  Cldrad_props%cldsct(:,:,:,8)
       used = send_data ( id_sct_cld_nir, cloud, Time_diag, is, js, &
                           1, rmask=mask )
    endif

    if ( id_asymm_cld_nir > 0  ) then
      cloud(:,:,:) = 100.0* Cldrad_props%cldasymm(:,:,:,8)
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

     id_abs_cld_lw = &
     register_diag_field ( mod_name, 'abs_lw', axes(1:3), Time, &
                      'cloud abs coeff lw', 'percent',        &
                           missing_value=missing_value          )

     id_abs_cld_10u = &
   register_diag_field ( mod_name, 'abs_10u', axes(1:3), Time, &
                      'cloud abs coeff 10um band', 'percent',    &
                        missing_value=missing_value          )



!   id_em_cld = &
!   register_diag_field ( mod_name, 'em_cld', axes(1:3), Time, &
!                        'cloud emissivity', 'percent',        &
!                         missing_value=missing_value          )

!    if (trim(swform) == 'lhsw' ) then
    if (do_lhsw                ) then
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

      id_lsc_cld_amt = &
     register_diag_field ( mod_name, 'lsc_cld_amt', axes(1:3), Time, &
                         'lsc cloud amount', 'percent',             &
                      missing_value=missing_value            )


     id_lsc_cld_ext_uv = &
    register_diag_field ( mod_name, 'lsc_cld_ext_uv', axes(1:3), Time, &
                        '.27um lsc cloud ext coeff', 'km-1',          &
                      missing_value=missing_value            )
 

    id_lsc_cld_ext_vis = &
   register_diag_field ( mod_name, 'lsc_cld_ext_vis', axes(1:3), Time,&
                        '.55um lsc cloud ext coeff', 'km-1',          &
                        missing_value=missing_value            )
 
 
   id_lsc_cld_ext_nir = &
   register_diag_field ( mod_name, 'lsc_cld_ext_nir', axes(1:3), Time,&
                        '1.4um lsc cloud ext coeff', 'km-1',          &
                          missing_value=missing_value            )

    id_lsc_cld_sct_uv = &
    register_diag_field ( mod_name, 'lsc_cld_sct_uv', axes(1:3), Time, &
                        '.27um lsc cloud sct coeff', 'km-1',          &
                         missing_value=missing_value            )


     id_lsc_cld_sct_vis = &
   register_diag_field ( mod_name, 'lsc_cld_sct_vis', axes(1:3), Time,&
                        '.55um lsc cloud sct coeff', 'km-1',          &
                       missing_value=missing_value            )


   id_lsc_cld_sct_nir = &
   register_diag_field ( mod_name, 'lsc_cld_sct_nir', axes(1:3), Time,&
                        '1.4um lsc cloud sct coeff', 'km-1',          &
                          missing_value=missing_value            )
 
 
     id_lsc_cld_asymm_uv = &
 register_diag_field ( mod_name, 'lsc_cld_asymm_uv', axes(1:3), Time, &
                        '.27um lsc cloud asymm coeff', 'percent',     &

                       missing_value=missing_value            )


     id_lsc_cld_asymm_vis = &
 register_diag_field ( mod_name, 'lsc_cld_asymm_vis', axes(1:3), Time,&
                       '.55um lsc cloud asymm coeff', 'percent',     &
                          missing_value=missing_value            )


     id_lsc_cld_asymm_nir = &
    register_diag_field ( mod_name, 'lsc_cld_sct_nir', axes(1:3), Time,&
                         '1.4um lsc cloud sct coeff', 'percent',       &
                         missing_value=missing_value            )


    id_abs_lsc_cld_lw = &
       register_diag_field ( mod_name, 'lsc_abs_lw', axes(1:3), Time, &
                          'lsc cloud abs coeff lw', 'percent',        &
                           missing_value=missing_value          )

       id_abs_lsc_cld_10u = &
     register_diag_field ( mod_name, 'lsc_abs_10u', axes(1:3), Time, &
               'lsc cloud abs coeff 10um band', 'percent',    &
                     missing_value=missing_value          )


  id_cell_cld_amt = &
     register_diag_field ( mod_name, 'cell_cld_amt', axes(1:3), Time, &
                          'cell cloud amount', 'percent',             &
                       missing_value=missing_value            )


    id_cell_cld_ext_uv = &
   register_diag_field ( mod_name, 'cell_cld_ext_uv', axes(1:3), Time, &
                      '.27um cell cloud ext coeff', 'km-1',        &
                        missing_value=missing_value            )


    id_cell_cld_ext_vis = &
 register_diag_field ( mod_name, 'cell_cld_ext_vis', axes(1:3), Time,&
                      '.55um cell cloud ext coeff', 'km-1',        &
              missing_value=missing_value            )


     id_cell_cld_ext_nir = &
   register_diag_field ( mod_name, 'cell_cld_ext_nir', axes(1:3), Time,&
                         '1.4um cell cloud ext coeff', 'km-1',        &
                          missing_value=missing_value            )


     id_cell_cld_sct_uv = &
   register_diag_field ( mod_name, 'cell_cld_sct_uv', axes(1:3), Time, &
                         '.27um cell cloud sct coeff', 'km-1',        &
                          missing_value=missing_value            )


 id_cell_cld_sct_vis = &
  register_diag_field ( mod_name, 'cell_cld_sct_vis', axes(1:3), Time,&
                         '.55um cell cloud sct coeff', 'km-1',        &
                      missing_value=missing_value            )

    id_cell_cld_sct_nir = &
  register_diag_field ( mod_name, 'cell_cld_sct_nir', axes(1:3), Time,&
                         '1.4um cell cloud sct coeff', 'km-1',        &
                          missing_value=missing_value            )


   id_cell_cld_asymm_uv = &
 register_diag_field ( mod_name, 'cell_cld_asymm_uv', axes(1:3), Time, &
                         '.27um cell cloud asymm coeff', 'percent',   &
                      missing_value=missing_value            )

 
   id_cell_cld_asymm_vis = &
 register_diag_field ( mod_name, 'cell_cld_asymm_vis', axes(1:3), Time,&
                     '.55um cell cloud asymm coeff', 'percent',   &
                       missing_value=missing_value            )


 id_cell_cld_asymm_nir = &
  register_diag_field ( mod_name, 'cell_cld_sct_nir', axes(1:3), Time,&
                        '1.4um cell cloud sct coeff', 'percent',     &
                       missing_value=missing_value            )


  id_abs_cell_cld_lw = &
    register_diag_field ( mod_name, 'cell_abs_lw', axes(1:3), Time, &
                         'cell cloud abs coeff lw', 'percent',        &
                         missing_value=missing_value          )

    id_abs_cell_cld_10u = &
   register_diag_field ( mod_name, 'cell_abs_10u', axes(1:3), Time, &
                        'cell cloud abs coeff 10um band', 'percent',  &
                   missing_value=missing_value          )


    id_meso_cld_amt = &
    register_diag_field ( mod_name, 'meso_cld_amt', axes(1:3), Time, &
                          'meso cloud amount', 'percent',             &
                      missing_value=missing_value            )


     id_meso_cld_ext_uv = &
  register_diag_field ( mod_name, 'meso_cld_ext_uv', axes(1:3), Time, &
                        '.27um meso cloud ext coeff', 'km-1',        &
                       missing_value=missing_value            )


   id_meso_cld_ext_vis = &
  register_diag_field ( mod_name, 'meso_cld_ext_vis', axes(1:3), Time,&
                       '.55um meso cloud ext coeff', 'km-1',        &
                       missing_value=missing_value            )


    id_meso_cld_ext_nir = &
  register_diag_field ( mod_name, 'meso_cld_ext_nir', axes(1:3), Time,&
                     '1.4um meso cloud ext coeff', 'km-1',        &
             missing_value=missing_value            )


     id_meso_cld_sct_uv = &
  register_diag_field ( mod_name, 'meso_cld_sct_uv', axes(1:3), Time, &
                         '.27um meso cloud sct coeff', 'km-1',        &
                          missing_value=missing_value            )


   id_meso_cld_sct_vis = &
  register_diag_field ( mod_name, 'meso_cld_sct_vis', axes(1:3), Time,&
                       '.55um meso cloud sct coeff', 'km-1',        &
                        missing_value=missing_value            )


   id_meso_cld_sct_nir = &
   register_diag_field ( mod_name, 'meso_cld_sct_nir', axes(1:3), Time,&
                   '1.4um meso cloud sct coeff', 'km-1',        &
                        missing_value=missing_value            )


    id_meso_cld_asymm_uv = &
 register_diag_field ( mod_name, 'meso_cld_asymm_uv', axes(1:3), Time, &
                         '.27um meso cloud asymm coeff', 'percent',   &
                       missing_value=missing_value            )


    id_meso_cld_asymm_vis = &
 register_diag_field ( mod_name, 'meso_cld_asymm_vis', axes(1:3), Time,&
                     '.55um meso cloud asymm coeff', 'percent',   &
                      missing_value=missing_value            )


    id_meso_cld_asymm_nir = &
  register_diag_field ( mod_name, 'meso_cld_sct_nir', axes(1:3), Time,&
                   '1.4um meso cloud sct coeff', 'percent',     &
                      missing_value=missing_value            )


      id_abs_meso_cld_lw = &
  register_diag_field ( mod_name, 'meso_abs_lw', axes(1:3), Time, &
                     'meso cloud abs coeff lw', 'percent',        &
                       missing_value=missing_value          )

     id_abs_meso_cld_10u = &
    register_diag_field ( mod_name, 'meso_abs_10u', axes(1:3), Time, &
                   'meso cloud abs coeff 10um band', 'percent',    &
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
 
!  do j=JSRAD,JERAD
!   do i=ISRAD,IERAD
!     do k = KSRAD, KERAD
   do j=1,jx         
    do i=1,ix         
      do k = 1,kx             
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

