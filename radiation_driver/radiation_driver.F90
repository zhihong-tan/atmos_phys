                module radiation_driver_mod


use utilities_mod,         only: get_my_pe, FATAL,  WARNING, NOTE,   &
                                 close_file, read_data, write_data,   &
                                 open_file,  print_version_number,  &
                                 get_num_pes, &
!                                get_root_pe, &
                                 check_nml_error, file_exist,  &
                                 utilities_init,  error_mesg
use diag_manager_mod,      only: register_diag_field, send_data, &
                                 diag_manager_init
use time_manager_mod,      only: time_type, set_date, set_time,  &
                                 get_time,    operator(+),       &
!                                time_manager_init, &
                                 operator(-), operator(/=), get_date,&
                                 operator(<), operator(>=), operator(>)
use rad_utilities_mod,     only: radiation_control_type, Rad_control, &
                                 radiative_gases_type, &
                                 cldrad_properties_type, &
                                 astronomy_type, &
                                 cld_diagnostics_type, &
                                 atmos_input_type, &
                                 aerosol_properties_type,   &
                                 Aerosol_props, aerosol_type, &
                                 sw_output_type, lw_output_type, &
                                 rad_output_type, &
                                 longwave_control_type, Lw_control, &
                                 shortwave_control_type, Sw_control, &
                                 fsrad_output_type, &
                                 environment_type, Environment, &
                                 cloudrad_control_type, Cldrad_control,&
                                 rad_utilities_init
use strat_cloud_mod,       only: add_strat_tend, strat_init
use edt_mod,               only: edt_init, edt_tend
use entrain_mod,           only: entrain_init, entrain_tend
use diag_integral_mod,     only: diag_integral_field_init, &
                                 sum_diag_integral_field
use sat_vapor_pres_mod,    only: sat_vapor_pres_init, lookup_es
use constants_mod,         only: RDGAS, RVGAS, STEFAN, GRAV, seconds_per_day
!use mpp_mod,               only: mpp_clock_id, mpp_clock_begin, &
use fms_mod,               only: mpp_clock_id, mpp_clock_begin, &
                                 mpp_clock_end, CLOCK_MODULE,    &
                                 MPP_CLOCK_SYNC

! original_fms_rad radiation package:

use original_fms_rad_mod,  only: original_fms_rad_init,  &
                                 original_fms_rad_end, &
                                 original_fms_rad

! the following mods are components of the sea_esf_rad radiation 
! package:

use sea_esf_rad_mod,       only: sea_esf_rad_init, sea_esf_rad_end, &
                                 sea_esf_rad
use rad_output_file_mod,   only: write_rad_output_file,    &
                                 rad_output_file_init, &
                                 rad_output_file_end
use radiative_gases_mod,   only: define_radiative_gases,  &
                                 radiative_gases_init,  &
                                 radiative_gases_end
use cloudrad_package_mod,  only: clouddrvr, cloudrad_package_alloc_rad,&
                                 cloudrad_package_end, &
                                 cloudrad_package_init
use astronomy_mod,         only: astronomy_init, &
                                 astronomy_end, &
                                 annual_mean_solar, daily_mean_solar, &
                                 diurnal_solar
use aerosol_mod,           only: aerosol_init, aerosol_driver,  &
                                 aerosol_end



implicit none 
private 


!----------------------------------------------------------------------
!    radiation_driver_mod is the interface between physics_driver_mod
!    and either the original_fms_rad or sea_esf_rad radiation package.
!    it provides radiative heating rates, boundary radiative fluxes,
!    and other radiation package output fields (cosine of zenith angle)
!    needed elsewhere in the modeling system. 
!
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!------------ version number for this module --------------------------

character(len=128) :: version = &
'$Id: radiation_driver.F90,v 1.7 2003/04/09 20:58:19 fms Exp $'
character(len=128) :: tag = '$Name: inchon $'


!---------------------------------------------------------------------
!------ interfaces -----

public    radiation_driver_init, radiation_driver,   &
          radiation_driver_end,  get_dtrad 

private  & 
          allocate_average, initialize_average, read_restart_file,  &
          initialize_diagnostic_integrals, diag_field_init, &
          define_rad_times, obtain_astronomy_variables, &
          define_atmos_input_fields, calculate_auxiliary_variables,  &
          time_average_input_data, compute_average, return_average,  &
          radiation_calc, model_to_radiation_grid,   &
          radiation_to_model_grid, deallocate_calc_arrays,   &
          update_rad_fields, radiation_netcdf, deallocate_arrays


!-----------------------------------------------------------------------
!------- namelist ---------

integer ::  rad_time_step = 14400     !  radiative time step in seconds
integer ::  offset = 1                !  offset for radiative time step 
                                      !  (in secs). note if offset=1 
                                      !  radiation will be done on the 
                                      !  first time step after t=0, 
                                      !  rad_time_step, 2*rad_time_step,
                                      !  ... .  for now do not change 
                                      !  this value.
logical ::  do_average= .false.       !  are atmospheric and astronomic
                                      !  input fields time-averaged ?
logical ::  do_average_gases = .false.!  are radiative gas inputs 
                                      !  time-averaged ?
logical ::  do_average_clouds= .false.!  are cloud property input fields
                                      !  time-averaged ?
logical ::  do_clear_sky_pass= .false.!  are the clear-sky radiation
                                      !  diagnostics to be calculated ?
character(len=24) ::    &
            zenith_spec = '      '    !  string defining how zenith 
                                      !  angle is computed. acceptable
                                      !  values: 'daily_mean', 'annual_
                                      !  mean', 'diurnally_varying'
character(len=16) ::   &
                rad_package='sea_esf' !  string defining the radiation
                                      !  package being used. acceptable
                                      !  values : 'sea_esf', 
                                      !  'original_fms'     
logical ::    &
         calc_hemi_integrals = .false.!  are hemispheric integrals 
                                      !  desired ? 
logical ::     &
        all_step_diagnostics = .false.!  are lw and sw radiative bdy
                                      !  fluxes and atmospheric heating 
                                      !  rates to be output on physics 
                                      !  steps ?
logical ::     &
         renormalize_sw_fluxes=.false.!  should sw fluxes and the zenith
                                      !  angle be renormalized on each 
                                      !  timestep because of the 
                                      !  movement of earth wrt the sun ?
integer, dimension(6) ::    &
    rad_date = (/ 0, 0, 0, 0, 0, 0 /) !  fixed date for which radiation
                                      !  is to be valid (applies to
                                      !  solar info, ozone, clouds)
                                      !  [yr, mo, day, hr, min, sec]
logical  ::  &
         all_level_radiation = .true. !  is radiation to be calculated 
                                      !  at all model levels ?
integer ::    &
          topmost_radiation_level=-99 !  if all_level_radiation is 
                                      !  false., this is the lowest
                                      !  model index at which radiation
                                      !  is calculated
logical ::    &
          drop_upper_levels = .false. !  if all_level_radiation is false
                                      !  and drop_upper_levels is true,
                                      !  radiation will be calculated
                                      !  at all model levels from
                                      !  topmost_radiation_level to the
                                      !  surface
logical ::  &
         do_aerosol = .false.         !  are aerosols included in rad-
                                      !  iation calculation ?
logical ::  &
         all_column_radiation = .true.!  is radiation to be calculated
                                      !  in all model columns ?
logical :: rsd=.false.                !  (repeat same day) - call 
                                      !  radiation for the specified 
                                      !  rad_date (yr,mo,day), but run 
                                      !  through the diurnal cycle (hr,
                                      !  min,sec)

logical :: use_mixing_ratio = .false. !  assumes q is mixing ratio
                                      !  rather than specific humidity
real    :: solar_constant = 1365.0    !  annual mean solar flux at top 
                                      !  of atmosphere [ W/(m**2) ]



namelist /radiation_driver_nml/ rad_time_step, offset, do_average, &
                                do_average_gases, do_average_clouds, &
                                do_clear_sky_pass, zenith_spec,  &
                                rad_package, calc_hemi_integrals, &
                                all_step_diagnostics, &
                                renormalize_sw_fluxes, &
                                rad_date, all_level_radiation, &
                                topmost_radiation_level,    &
                                drop_upper_levels,    &
                                do_aerosol,  &
                                all_column_radiation, &
                                rsd, use_mixing_ratio, solar_constant

!---------------------------------------------------------------------
!---- public data ----



!---------------------------------------------------------------------
!---- private data ----


!---------------------------------------------------------------------
!    logical  flags.
!---------------------------------------------------------------------
logical ::  radiation_driver_initialized =   &
                                 .false.    ! module initialized?
logical ::  do_rad                          ! radiation step ?         
logical ::  use_rad_date                    ! specify time of radiation
                                            ! independent of model time?
logical ::  do_sea_esf_rad                  ! using sea_esf_rad package?
logical ::  do_clear_sky_pass               ! doing clear-sky radiative
                                            ! flux calculation ?
logical ::  do_diurnal                      ! using diurnally-varying
                                            ! solar zenith angle ?
logical ::  do_annual                       ! using annual mean solar
                                            ! zenith angle ?
logical ::  do_daily_mean                   ! using daily mean solar
                                            ! zenith angle ?

!---------------------------------------------------------------------
!    list of restart versions of radiation_driver.res readable by this 
!    module, dependent on which radiation package is activated. restart 
!    version 1 of radiation_driver.res is not readable by this module; 
!    however, restart version 1 of sea_esf_rad.res IS readable within 
!    this code.
!---------------------------------------------------------------------
integer, dimension(3) :: restart_versions     = (/ 2, 3, 4 /)
integer, dimension(3) :: restart_versions_sea = (/ 2, 3, 4 /)


!-----------------------------------------------------------------------
!    these arrays must be preserved across timesteps:
!
!    Rad_output is a rad_output_type variable with the following 
!    components:
!          tdt_rad        radiative (sw + lw) heating rate
!          flux_sw_surf   net (down-up) sw flux at surface
!          flux_lw_surf   downward lw flux at surface
!          coszen_angle   cosine of the zenith angle (used for the 
!                         last radiation calculation)
!          tdt_rad_clr    net radiative heating rate in the absence of
!                         cloud
!          tdtsw          shortwave heating rate
!          tdtsw_clr      shortwave heating rate in he absence of cloud
!          tdtlw          longwave heating rate

!    solar_save is used when renormalize_sw_fluxes is active, to save
!    the solar factor (fracday*cosz/r**2) from the previous radiation
!    step so that the radiative forcing terms may be adjusted on each
!    timestep to reflect the current solar forcing.
!
!    sw_heating_clr, tot_heating_clr_save, sw_heating_save, 
!    tot_heating_save, flux_sw_surf_save are the radiative forcing terms
!    on radiation steps which also must be saved when renormalization 
!    is activated. 
!
!    the ***sw_save arrays are currently saved so that their values may
!    be adjusted during sw renormalization for diagnostic purposes.
!                               
!    the **lw_save arrays are currently saved so that they may be output
!    in the diagnostics file on every physics step, if desired, so that
!    when renormalize_sw_fluxes is active, total radiative terms may be
!    easily generated.
!-----------------------------------------------------------------------

type(rad_output_type)               ::  Rad_output
real, allocatable, dimension(:,:)   ::  solar_save, flux_sw_surf_save
real, allocatable, dimension(:,:,:) ::  sw_heating_save,    &
                                        tot_heating_save,  &
                                        sw_heating_clr_save, &
                                        tot_heating_clr_save, &
                                        dfsw_save, ufsw_save, fsw_save,&
                                        hsw_save, dfswcf_save,   &
                                        ufswcf_save, fswcf_save, &
                                        hswcf_save

real, allocatable, dimension(:,:,:) ::  tdtlw_save, tdtlw_clr_save
real, allocatable, dimension(:,:)   ::  olr_save, lwups_save, &
                                        lwdns_save, olr_clr_save, &
                                        lwups_clr_save, lwdns_clr_save

!-----------------------------------------------------------------------
!    time-step-related constants
 
integer    :: rad_alarm      !  time interval until the next radiation 
                             !  calculation (seconds)
integer    :: num_pts        !  counter for current number of grid 
                             !  boxes processed (when num_pts=0 or 
                             !  num_pts=total_pts certain things happen)
integer    :: total_pts      !  number of grid boxes to be processed 
                             !  every time step (note: all grid boxes 
                             !  must be processed every time step)


!-----------------------------------------------------------------------
!   accumulation arrays for time averaged input data 

real,    allocatable, dimension(:,:,:)   :: psum, tsum, qsum
real,    allocatable, dimension(:,:)     :: asum, csum, ssum,   &
                                            fsum, rrsum
integer, allocatable, dimension(:,:)     :: nsum
real,    allocatable, dimension(:,:,:)   :: gas_component_sum
real,    allocatable, dimension(:,:,:,:) :: cld_component_sum



!-----------------------------------------------------------------------
!    diagnostics variables

character(len=16)            :: mod_name = 'radiation'
integer                      :: id_alb_sfc, id_cosz, id_fracday 

integer, dimension(2)        :: id_tdt_sw,   id_tdt_lw,  &
                                id_swdn_toa, id_swup_toa, id_olr, &
                                id_swup_sfc, id_swdn_sfc,         &
                                id_lwup_sfc, id_lwdn_sfc

real                         :: missing_value = -999.
character(len=8)             :: std_digits   = 'f8.3'
character(len=8)             :: extra_digits = 'f16.11'

!-----------------------------------------------------------------------
!    timing clocks       

integer                      :: misc_clock, clouds_clock, aerosol_clock

!--------------------------------------------------------------------
! miscellaneous variables and indices

integer        ::  id         !  number of grid points in x direction 
                              !  (on processor)
integer        ::  jd         !  number of grid points in y direction 
                              !  (on processor)
integer        ::  kmin = 1   !  starting k index for model grid
integer        ::  kmax       !  number of vertical model layers
integer        ::  ks         !  model grid coordinate of top level
                              !  at which radiation is calculated 
                              !  (topmost_radiation_level)
integer        ::  ke         !  model grid coordinate of bottommost
                              !  level at which radiation is calculated

integer        ::  ksrad=1    !  always set to 1
integer        ::  kerad      !  number of layers in radiation grid

real           ::  rh2o_lower_limit_orig=3.0E-06
                              !  smallest value of h2o mixing ratio 
                              !  allowed with original_fms_rad package
!real           ::  rh2o_lower_limit_seaesf=3.0E-06
real           ::  rh2o_lower_limit_seaesf=2.0E-07
                              !  smallest value of h2o mixing ratio 
                              !  allowed with sea_esf_rad package
real           ::  rh2o_lower_limit
                              !  smallest value of h2o mixing ratio 
                              !  allowed in the current experiment
real           ::  temp_lower_limit=100.0  ! [ K ]
                              !  smallest value of temperature      
                              !  allowed in the current experiment
real           ::  temp_upper_limit=370.00  ! [ K ]
                              !  largest value of temperature 
                              !  allowed in the current experiment

real           ::  surf_flx_init=50.0  ! [w / m^2 ]
                              !  value to which surface lw and sw fluxes
                              !  are set in the absence of a .res file
                              !  containing them

real           ::  coszen_angle_init=0.50
                              !  value to which cosine of zenith angle  
                              !  is set in the absence of a .res file
                              !  containing it

real           ::  log_p_at_top=2.0
                              !  assumed value of ln of ratio of pres-
                              !  sure at flux level 2 to that at model
                              !  top (needed for deltaz calculation,
                              !  is infinite for model top at p = 0.0,
                              !  this value is used to give a reasonable
                              !  deltaz)
real,parameter ::  d608 = (RVGAS-RDGAS)/RDGAS
                              !  virtual temperature factor  
real,parameter ::  d622 = RDGAS/RVGAS
                              ! ratio of gas constants - dry air to 
                              ! water vapor
real,parameter ::  d378 = 1.0 - d622  
                              ! 1 - gas constant ratio

!---------------------------------------------------------------------
!---------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!######################################################################

subroutine radiation_driver_init (lonb, latb, pref, axes, Time)

!---------------------------------------------------------------------
!   radiation_driver_init is the constructor for radiation_driver_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real, dimension(:),      intent(in)  :: lonb, latb
real, dimension(:,:),    intent(in)  :: pref
integer, dimension(4),   intent(in)  :: axes
type(time_type),         intent(in)  :: Time
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       lonb      array of model longitudes on cell boundaries [radians]
!
!       latb      array of model latitudes at cell boundaries [radians]
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [pascals]
!
!       axes      diagnostic variable axes
!
!       Time      current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables

      integer         ::   unit, io, ierr
      logical         ::   end
!---------------------------------------------------------------------

      
!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (radiation_driver_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call utilities_init
      call rad_utilities_init
!! routine not yet compliant
!      call strat_init
      if (Environment%running_gcm .or. &
          Environment%running_sa_model) then
        call diag_manager_init
      endif
!! routine not yet existent
!      call time_manager_init
      call sat_vapor_pres_init
!! routine not public         
!     call constants_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit = open_file (file='input.nml', action='read')
        ierr=1; do while (ierr /= 0)
          read  (unit, nml=radiation_driver_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'radiation_driver_nml')
        enddo
  10    call close_file (unit)
      endif

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      unit = open_file ('logfile.out', action='append')
      if ( get_my_pe() == 0 ) then
!     if ( get_my_pe() == get_root_pe() ) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit, nml=radiation_driver_nml)
      endif
      call close_file (unit)

!---------------------------------------------------------------------
!    set logical variable defining the radiation scheme desired from the
!    namelist-input character string. set lower limit to water vapor 
!    mixing ratio that the radiation code will see, to assure keeping 
!    within radiation lookup tables. exit if value is invalid.
!---------------------------------------------------------------------
      if (rad_package == 'original_fms') then
        do_sea_esf_rad = .false.
        rh2o_lower_limit = rh2o_lower_limit_orig
      else if (rad_package == 'sea_esf') then
        do_sea_esf_rad = .true.
        rh2o_lower_limit = rh2o_lower_limit_seaesf
      else
        call error_mesg ('radiation_driver_mod',  &
           'string provided for rad_package is not valid', FATAL)
      endif

!--------------------------------------------------------------------
!    set logical variables defining how the solar zenith angle is to
!    be  defined from the namelist-input character string.  exit if the
!    character string is invalid.
!--------------------------------------------------------------------
      if (zenith_spec == 'diurnally_varying') then
        do_diurnal = .true.
        do_annual = .false.
        do_daily_mean = .false.
      else if (zenith_spec == 'daily_mean') then
        do_diurnal = .false.
        do_annual = .false.
        do_daily_mean = .true.
      else if (zenith_spec == 'annual_mean') then
        do_diurnal = .false.
        do_annual = .true.
        do_daily_mean = .false.
      else
        call error_mesg ('radiation_driver_mod', &    
            'string provided for zenith_spec is invalid', FATAL)
      endif

!--------------------------------------------------------------------
!    verify that if it is desired to time-average the radiative gas 
!    input fields that the astronomy and atmospheric input fields are 
!    also being averaged.
!--------------------------------------------------------------------
      if (do_average_gases .and. .not. do_average) then
        call error_mesg( 'radiation_driver_mod', &
         ' cannot activate do_average_gases without having do_average&
	    & also true', FATAL)
      endif

!--------------------------------------------------------------------
!    verify that if it is desired to time-average the cloud radiative
!    properties input fields that the astronomy and atmospheric input 
!    fields are also being averaged.
!--------------------------------------------------------------------
      if (do_average_clouds .and. .not. do_average) then
        call error_mesg( 'radiation_driver_mod', &
         ' cannot activate do_average_clouds without having do_average&
	    & also true', FATAL)
      endif

!---------------------------------------------------------------------
!    verify that radiation has been requested at all model levels and in
!    all model columns when the original fms radiation is activated.    
!    verify that renormalize_sw_fluxes has not been requested along
!    with the original fms radiation package. verify that all_step_diag-
!    nostics has not been requested with the original fms radiation
!    package.
!---------------------------------------------------------------------
      if (.not. do_sea_esf_rad) then
        if (.not. all_level_radiation .or. &
            .not. all_column_radiation) then
          call error_mesg ( 'radiation_driver_mod', &
          ' must specify all_level_radiation and all_column_radiation&
            & as true when using original fms radiation', FATAL)
        endif
        if (renormalize_sw_fluxes) then
          call error_mesg ( 'radiation_driver_mod', &
           ' cannot renormalize shortwave fluxes with original_fms &
	     &radiation package.', FATAL)
        endif
        if (all_step_diagnostics) then
          call error_mesg ( 'radiation_driver_mod', &
            ' cannot request all_step_diagnostics with original_fms &
              &radiation package.', FATAL)
        endif
      endif

!---------------------------------------------------------------------
!    can only renormalize shortwave fluxes when diurnally_varying
!    radiation is used and when using_fms_periphs.
!---------------------------------------------------------------------
     if (renormalize_sw_fluxes .and. .not. do_diurnal) then
       call error_mesg ('radiation_driver_mod',  &
        ' can only renormalize sw fluxes when using diurnally-varying&
	  & solar radiation', FATAL)
     endif
     if (renormalize_sw_fluxes .and.     &
                            .not. Environment%using_fms_periphs) then
       call error_mesg ('radiation_driver_mod',  &
        ' cannot renormalize sw fluxes when using skyhi peripherals', &
           FATAL)
     endif

!---------------------------------------------------------------------
!    verify that a valid radiation time step has been specified.
!---------------------------------------------------------------------
      if (rad_time_step <= 0) then
        call error_mesg ('radiation_driver_mod', &
            ' radiation timestep must be set to a positive integer', &
              FATAL)
      endif

!---------------------------------------------------------------------
!    verify that averaging of radiative input data is not active with
!    renormalization of shortwave fluxes.
!---------------------------------------------------------------------
      if (renormalize_sw_fluxes .and. do_average) then
        call error_mesg ('radiation_driver_mod',  &
           ' cannot activate input data averaging and renormalizing &
	    &shortwave fluxes at same time', FATAL)
      endif

!---------------------------------------------------------------------
!    verify that a valid offset for the radiation time step has been 
!    specified.
!---------------------------------------------------------------------
      if (offset /= 1) then
        call error_mesg ( 'radiation_driver_mod',  &
          ' offset should be set to 1 -- use of other values not yet &
	   &validated.', FATAL)
      endif

!---------------------------------------------------------------------
!    define the model dimensions on the local processor.
!---------------------------------------------------------------------
      id=size(lonb,1) - 1; jd=size(latb,1) - 1
      kmax  = size(pref,1) - 1 

!---------------------------------------------------------------------
!    check for consistency if drop_upper_levels is activated.
!----------------------------------------------------------------------
      if (drop_upper_levels) then
        if (all_level_radiation) then
          call error_mesg ( 'radiation_driver_mod',  &
            ' drop_upper_levels and all_level_radiation are &
	       &incompatible', FATAL)
        endif
      endif

!---------------------------------------------------------------------
!    define the starting and ending vertical indices of the radiation
!    grid. if all_level_radiation is .true., then radiation is done
!    at all model levels. ks, ke are model-based coordinates, while
!    ksrad and kerad are radiation-grid based coordinates (ksrad always
!    is equal to 1). define the pressure arrays on the radiation grid,
!    so they may be passed into the sea_esf_rad package.
!---------------------------------------------------------------------
      if (all_level_radiation) then
        ks = 1
        ke = kmax
        kerad = kmax
        topmost_radiation_level = 1
        if (drop_upper_levels) then
          call error_mesg ( 'radiation_driver_mod',  &
            ' drop_upper_levels and all_level_radiation are &
	       &incompatible', FATAL)
        endif
      else
        if (topmost_radiation_level <= 0) then
          call error_mesg ('radiation_driver_mod', &
              ' when all_level_radiation is .false., topmost_radiation&
	      &_level must be specified as a positive integer.', FATAL)
        endif
        if (drop_upper_levels) then
          ks = topmost_radiation_level
          ke = kmax
          kerad = ke - ks + 1
          call error_mesg ( ' radiation_driver_mod', &
            ' code currently not validated for all_level_radiation = &
               &false.', FATAL)
        else
          call error_mesg ( ' radiation_driver_mod', &
            ' currently only drop_upper_levels is available as option &
	      &when all_level_radiation = false.', FATAL)
        endif
      endif
       
!---------------------------------------------------------------------
!    exit if all_column_radiation is not .true. -- this option is not
!    yet certified.
!---------------------------------------------------------------------
      if (.not. all_column_radiation) then
        call error_mesg ('radiation_driver_mod',  &
          ' code currently not validated for all_column_radiation = &
           &false.', FATAL)
      endif

!----------------------------------------------------------------------
!    be sure both reference pressure profiles have been provided.
!----------------------------------------------------------------------
      if (size(pref,2) /= 2)    &
        call error_mesg ('radiation_driver_mod', &
         'must provide two reference pressure profiles (pref).', FATAL)

!---------------------------------------------------------------------
!    the renormalize_sw_fluxes option is not allowed when running the 
!    standalone code.
!---------------------------------------------------------------------
       if ((Environment%running_standalone .or. &
            Environment%running_sa_model)  .and.    &
           renormalize_sw_fluxes) then
        call error_mesg ( 'radiation_driver_mod', &
               'renormalizing sw fluxes when running standalone&
	       & radiation code is not available', FATAL)
      endif

!---------------------------------------------------------------------
!    the do_average option is not allowed when running the 
!    standalone code.
!---------------------------------------------------------------------
      if ((Environment%running_standalone .or.  &
           Environment%running_sa_model) .and. do_average) then
        call error_mesg ( 'radiation_driver_mod', &
               'averaging of  input fields when running standalone&
	       & radiation code is not available', FATAL)
      endif

!---------------------------------------------------------------------
!    define default values for the radiation time step and when the 
!    radiation package should first be called. the time to next call
!    radiation may be modified later if the radiation timestep has 
!    changed from that used previously and supplied in the .res file.
!---------------------------------------------------------------------
      if (Environment%running_gcm) then
        if (offset > 0) then
          rad_alarm = offset
        else
          rad_alarm = rad_time_step
        endif
      endif

!---------------------------------------------------------------------
!    allocate space for variables which must be saved when sw fluxes
!    are renormalized or diagnostics are desired to be output on every
!    physics step.
!---------------------------------------------------------------------
      if (renormalize_sw_fluxes .or. all_step_diagnostics) then
        allocate (solar_save             (id,jd))
        allocate (flux_sw_surf_save      (id,jd))
        allocate (sw_heating_save        (id,jd,kmax))
        allocate (tot_heating_save       (id,jd,kmax))
        allocate (dfsw_save              (id,jd,kmax+1))
        allocate (ufsw_save              (id,jd,kmax+1))
        allocate ( fsw_save              (id,jd,kmax+1))
        allocate ( hsw_save              (id,jd,kmax))
        if (do_clear_sky_pass) then 
          allocate (sw_heating_clr_save  (id,jd,kmax))
          allocate (tot_heating_clr_save (id,jd,kmax))
          allocate (dfswcf_save          (id,jd,kmax+1))
          allocate (ufswcf_save          (id,jd,kmax+1))
          allocate ( fswcf_save          (id,jd,kmax+1))
          allocate ( hswcf_save          (id,jd,kmax))
        endif
      endif

!---------------------------------------------------------------------
!    allocate space for variables which must be saved when lw fluxes
!    are to be output on every physics step.
!---------------------------------------------------------------------
      if (all_step_diagnostics) then
        allocate (olr_save             (id,jd))
        allocate (lwups_save           (id,jd))
        allocate (lwdns_save           (id,jd))
        allocate (tdtlw_save           (id,jd,kmax))
        if (do_clear_sky_pass) then
          allocate (olr_clr_save       (id,jd))
          allocate (lwups_clr_save     (id,jd))
          allocate (lwdns_clr_save     (id,jd))
          allocate (tdtlw_clr_save     (id,jd,kmax))
        endif
      endif

      if (Environment%running_gcm .or. &
          Environment%running_sa_model) then
!---------------------------------------------------------------------
!    allocate space for module variables to contain values which must
!    be saved between timesteps (these are used on every timestep,
!    but only calculated on radiation steps).
!---------------------------------------------------------------------
        allocate (Rad_output%tdt_rad     (id,jd,kmax))
        allocate (Rad_output%tdt_rad_clr (id,jd,kmax))
        allocate (Rad_output%tdtsw       (id,jd,kmax))
        allocate (Rad_output%tdtsw_clr   (id,jd,kmax))
        allocate (Rad_output%tdtlw       (id,jd,kmax))
        allocate (Rad_output%flux_sw_surf(id,jd))
        allocate (Rad_output%flux_lw_surf(id,jd))
        allocate (Rad_output%coszen_angle(id,jd))
        Rad_output%tdtsw     = 0.0              
        Rad_output%tdtsw_clr = 0.0             
        Rad_output%tdtlw     = 0.0             

!-----------------------------------------------------------------------
!    if input fields are to be time averaged, allocate space for and 
!    initialize the accumulation arrays. 
!-----------------------------------------------------------------------
        if (do_average) then
          call allocate_average
          call initialize_average
        endif

!-----------------------------------------------------------------------
!    if two radiation restart files exist, exit.
!-----------------------------------------------------------------------
        if ( file_exist('INPUT/sea_esf_rad.res')  .and.     &
             file_exist('INPUT/radiation_driver.res') ) then 
          call error_mesg ('radiation_driver_mod',  &
            ' both sea_esf_rad.res and radiation_driver.res files are&
	     & present in INPUT directory. which one to use ?', FATAL)
        endif

!-----------------------------------------------------------------------
!    if a valid restart file exists, call read_restart_file to read it.
!-----------------------------------------------------------------------
        if ( (do_sea_esf_rad .and.   &
             (file_exist('INPUT/sea_esf_rad.res')  .or. &
              file_exist('INPUT/radiation_driver.res') )  ) .or. &
             (.not. do_sea_esf_rad .and.   &
               file_exist('INPUT/radiation_driver.res') )  ) then
          call read_restart_file 
        else
!----------------------------------------------------------------------
!    if no restart file is present, initialize the needed fields until
!    the radiation package may be called. initial surface flux is set 
!    to 100 wm-2, and is only used for initial guess of sea ice temp.
!-----------------------------------------------------------------------
          Rad_output%tdt_rad       = 0.0
          Rad_output%tdtlw         = 0.0
          Rad_output%flux_sw_surf  = surf_flx_init
          Rad_output%flux_lw_surf  = surf_flx_init
          Rad_output%coszen_angle  = coszen_angle_init
          if (get_my_pe() == 0) then
            call error_mesg ('radiation_driver_mod', &
              'no acceptable radiation restart file present; therefore&
	           & will initialize input fields', NOTE)
	  endif
        endif
      endif   ! (running_gcm)

!--------------------------------------------------------------------
!    do package-specific initialization.
!--------------------------------------------------------------------
      if (do_sea_esf_rad) then 

!---------------------------------------------------------------------
!    define a control variable indicating whether the clear-sky forcing
!    should be calculated.
!---------------------------------------------------------------------
        Rad_control%do_totcld_forcing = do_clear_sky_pass

!---------------------------------------------------------------------
!    do the initialization specific to the sea_esf_rad radiation
!    package.
!---------------------------------------------------------------------
        call radiative_gases_init (pref, latb, lonb)

	Rad_control%do_aerosol = do_aerosol
        if (Rad_control%do_aerosol) then
          call aerosol_init (lonb, latb, kmax)
        endif

        if (ks /= 1 .or. kerad /= ke) then
          call error_mesg ('radiation_driver_mod', &
            ' verify that pref being sent to sea_esf_rad is on the &
	    &radiation grid and that all missing code has been &
	    &supplied.', FATAL)
        else
          call sea_esf_rad_init (lonb, latb, pref(ks:ke+1,:))
        endif
        if (Environment%running_gcm .or.    &
            Environment%running_sa_model) then
          call cloudrad_package_init (pref(ks:ke+1,:), lonb, latb,  &
                                      axes=axes, Time=Time)
          call rad_output_file_init  (axes=axes, Time=Time)
        else if (Environment%running_standalone) then
          call cloudrad_package_init (pref(ks:ke+1,:), lonb, latb)
          call rad_output_file_init
        endif

!---------------------------------------------------------------------
!    do the initialization specific to the original fms radiation
!    package.
!---------------------------------------------------------------------
      else
        call radiative_gases_init (pref, latb, lonb)
        call original_fms_rad_init (lonb, latb, pref, axes, Time, kmax)
      endif

!--------------------------------------------------------------------
!    initialize the astronomy_package.
!--------------------------------------------------------------------
      if (do_annual) then
        call astronomy_init (latb, lonb)
      else
        call astronomy_init
      endif

      if (Environment%running_gcm .or.  &
          Environment%running_sa_model) then
!---------------------------------------------------------------------
!    initialize the total number of points and points processed on this
!    step. 
!---------------------------------------------------------------------
        total_pts=id *jd 
        num_pts=0

!-----------------------------------------------------------------------
!    check if optional radiative date should be used.
!-----------------------------------------------------------------------
        if (rad_date(1) > 1900 .and.                        &
            rad_date(2) >   0  .and. rad_date(2) < 13 .and. &
            rad_date(3) >   0  .and. rad_date(3) < 32 ) then
          use_rad_date = .true.
        else
          use_rad_date = .false.
        endif

!----------------------------------------------------------------------
!    define characteristics of desired diagnostic integrals. 
!----------------------------------------------------------------------
        call initialize_diagnostic_integrals

!----------------------------------------------------------------------
!    register the desired netcdf output variables with the 
!    diagnostics_manager.
!----------------------------------------------------------------------
        call diag_field_init (Time, axes)
      endif

!---------------------------------------------------------------------
!    set flag to indicate that module has been successfully initialized.
!---------------------------------------------------------------------
      radiation_driver_initialized = .true.

!--------------------------------------------------------------------
!    initialize clocks to time portions of the code called from 
!    radiation_driver.
!--------------------------------------------------------------------
      misc_clock =    &
            mpp_clock_id ('   Physics_down: Radiation: misc', &
                grain = CLOCK_MODULE, flags = MPP_CLOCK_SYNC)
      clouds_clock =   &
            mpp_clock_id ('   Physics_down: Radiation: clds', &
               grain = CLOCK_MODULE, flags = MPP_CLOCK_SYNC)
      aerosol_clock =  &
            mpp_clock_id ('   Physics_down: Radiation: arsl', &
               grain = CLOCK_MODULE, flags = MPP_CLOCK_SYNC)

!--------------------------------------------------------------------


end subroutine radiation_driver_init



!#####################################################################

subroutine radiation_driver (is, ie, js, je, Time, Time_next,  &
                             lat, lon, pfull, phalf, t, q, ts,  &
                             land, albedo, tdt, flux_sw, flux_lw,  &
                             coszen, mask, kbot, cloud_ice_input,  &
                             cloud_water_input, cloud_area_input, &
                             cloudtemp, cloudvapor)

!---------------------------------------------------------------------
!    radiation_driver adds the radiative heating rate to the temperature
!    tendency and obtains the radiative boundary fluxes and cosine of 
!    the solar zenith angle to be used in the other component models.
!---------------------------------------------------------------------
 

!--------------------------------------------------------------------
integer,                 intent(in)          :: is, ie, js, je
type(time_type),         intent(in)          :: Time, Time_next
real, dimension(:,:),    intent(in)          :: lat, lon, ts, land,  &
                                                albedo
real, dimension(:,:,:),  intent(in)          :: pfull, phalf, t, q
real, dimension(:,:,:),  intent(inout)       :: tdt
real, dimension(:,:),    intent(out)         :: flux_sw, flux_lw, coszen
real, dimension(:,:,:),  intent(in),optional :: mask, &
                                                cloud_water_input, &
                                                cloud_ice_input, &
                                                cloud_area_input, &
                                                cloudtemp, cloudvapor
integer, dimension(:,:), intent(in),optional :: kbot
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Time         current model time [ time_type (days, seconds) ] 
!      Time_next    time on next timestep, used as stamp for diagnostic 
!                   output  [ time_type  (days, seconds) ]  
!      lat          latitude of model points  [ radians ]
!      lon          longitude of model points [ radians ]
!      ts           surface temperature  [ deg K ]
!      land         fraction of surface which covered by land 
!                   [ dimensionless ]
!      albedo       surface albedo  [ dimensionless ]
!      pfull        pressure at full levels [ kg / (m s^2) ]
!      phalf        pressure at half levels [ kg / (m s^2) ]
!      t            temperature at full levels [ deg K]
!      q            specific humidity of water vapor at full levels
!                   [ dimensionless ]
!
!   intent(inout) variables:
!
!      tdt        temperature tendency [  deg K / sec ]
!
!   intent (out) variables:
!
!      flux_sw    net shortwave surface flux (down-up) [ w / m^2 ]
!      flux_lw    downward longwave surface flux [ w / m^2 ]
!      coszen     cosine of the zenith angle which will be used for
!                 the next ocean_albedo calculation [ dimensionless ]
!
!   intent(in), optional variables:
!
!      mask               present when running eta vertical coordinate,
!                         mask to remove points below ground
!      kbot               present when running eta vertical coordinate,
!                         index of lowest model level above ground 
!      cloud_water_input  required for sa_gcm mode, may be present in 
!                         columns mode, not currently present in full
!                         gcm mode: SPECIFIC HUMIDITY of liquid water
!                         DAN, NOTE: in columns mode, you previously
!                         have used mixing ratio, so modified input
!                         files may be needed.
!      cloud_ice_input    required for sa_gcm mode, may be present in 
!                         columns mode, not currently present in full
!                         gcm mode: SPECIFIC HUMIDITY of ice water
!                         DAN, NOTE: in columns mode, you previously
!                         have used mixing ratio, so modified input
!                         files may be needed.
!      cloud_area_input   required in sa_gcm mode, absent in columns
!                         and full gcm mode: fractional cloud coverage
!      cloudtemp          required in sa_gcm mode, absent otherwise:
!                         temperature field to be used by cloud param-
!                         eterization 
!      cloudvapor         required in sa_gcm mode, absent otherwise:
!                         water vapor field to be used by cloud param-
!                         eterization 
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      real, dimension (ie-is+1, je-js+1) :: flux_ratio
      integer                            :: dt

      type(time_type)                    :: Dt_zen, Rad_time
      type(radiative_gases_type)         :: Rad_gases_mdl
      type(cldrad_properties_type)       :: Cldrad_props_mdl
      type(aerosol_type)                 :: Aerosol_mdl
      type(astronomy_type)               :: Astro_mdl, Astro2
      type(lw_output_type)               :: Lw_output
      type(atmos_input_type)             :: Atmos_input_mdl
      type(fsrad_output_type)            :: Fsrad_output
      type(sw_output_type)               :: Sw_output
      type(cld_diagnostics_type)         :: Cld_diagnostics 
 
!-------------------------------------------------------------------
!   local variables:
!
!      dt                physics time step (frequency of calling 
!                        radiation_driver)  [ seconds ]
!      flux_ratio        value  used to renormalize sw fluxes and 
!                        heating rates to account for earth-sun motion
!                        during the radiation timestep
!      Dt_zen            time interval over which  cosine of zenith 
!                        angle is to be averaged when diurnally-varying
!                        solar radiation is activated 
!                        [ time_type (days, seconds)]
!      Rad_time          time at which the climatologically-determined,
!                        time-varying input fields to radiation should 
!                        apply 
!                        [ time_type (days, seconds)]
!      Rad_gases_mdl     radiative gas input fields on the model grid
!                        [radiative_gases_type]
!      Cldrad_props_mdl  cloud radiative properties on model grid,
!                        [cldrad_properties_type]
!      Astro_mdl         astronomical properties on model grid, usually
!                        valid over radiation timestep
!                        [astronomy_type]
!      Astro2            astronomical properties on model grid, valid 
!                        over current physics timestep
!                        [astronomy_type]
!      Lw_output         sea-esf longwave output fields on model grid,
!                        [lw_output_type]
!      Atmos_input_mdl   atmospheric input fields on model grid,
!                        [atmos_input_type]
!      Fsrad_output      original fms radiation output fields on model
!                        grid, [fsrad_output_type]
!      Sw_output         sea-esf shortwave output fields on model grid,
!                        [sw_output_type]
!      Cld_diagnostics   microphysical diagnostic fields collected 
!                        for radiation_diag_mod analysis
!                        [cld_diagnostics_type]
!
!----------------------------------------------------------------------


!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. radiation_driver_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!-------------------------------------------------------------------
!    verify that optional arguments that are required when running in
!    sa-model mode are present.
!-------------------------------------------------------------------
      if (Environment%running_sa_model) then
        if (present(cloudtemp)   .and.   &
            present(cloudvapor)   .and.   &
            present(cloud_water_input) .and.   &
            present(cloud_ice_input) .and. &
            present(cloud_area_input)) then
        else
            call error_mesg ('radiation_driver_mod', &
                    'must pass ql, qi, cf, cloudtemp and cloudvapor to &
                    &radiation_driver when running sa model', FATAL)
        endif
      endif

!------------------------------------------------------------------
!    call define_rad_times to define time-related variables needed by 
!    the radiation code. determine whether or not new radiative fluxes
!    and heating rates are to be calculated on this timestep (do_rad). 
!    define the time to be used to determine values for climatologically
!    supplied time-varying fields needed in the radiation calculation 
!    (Rad_time). define the averaging interval to be used in calculating
!    the cosine of the zenith angle when diurnally-varying solar rad-
!    iation is activated (Dt_zen). define the physics timestep (dt).
!--------------------------------------------------------------------
      call mpp_clock_begin (misc_clock)
      call define_rad_times (Time, Time_next, do_rad, Rad_time,   &
                             Dt_zen, dt)

!---------------------------------------------------------------------
!    if this is a radiation step, or if the astronomical inputs to
!    radiation (solar, cosz, fracday, rrsun) need to be recalculated 
!    because of time averaging or renormalization, call 
!    obtain_astronomy_variables to do so.
!---------------------------------------------------------------------
      if ( do_rad .or. do_average .or. renormalize_sw_fluxes ) then 
        Astro_mdl%do_diurnal = do_diurnal
        Astro_mdl%do_annual  = do_annual 
        Astro_mdl%do_daily_mean  = do_daily_mean
        Astro_mdl%solar_constant = solar_constant
        call obtain_astronomy_variables (is, ie, js, je,  dt, lat,   &
                                         lon, Rad_time, Dt_zen,   &
                                         Astro_mdl, Astro2, Rad_output)  
      endif

!--------------------------------------------------------------------
!    if time_averaging is active, or if this is a step on which to cal-
!    culate radiative tendencies, pass the model pressure (pfull, 
!    phalf), temperature (t, ts), specific humidity (q), surface albedo
!    (albedo), land surface amount (land), and optionally (currently in 
!    the standalone code, not the gcm) cloud water (cloud_water_input) 
!    and cloud ice (cloud_ice_input) to define_atmos_input_fields, which
!    will put these fields into the form desired by the radiation code 
!    and store them as components of the derived-type structure 
!    Atmos_input_mdl.
!---------------------------------------------------------------------
      if (do_average .or. do_rad) then
        if (present(cloud_water_input) .and.  &
            present(cloud_ice_input) ) then
          if (present (cloudtemp)) then     ! sa_gcm mode
            call define_atmos_input_fields     &
                                (is, ie, js, je, pfull, phalf, t, q,  &
                                 ts, albedo, land, Atmos_input_mdl,    &
                                 cloud_water_input=cloud_water_input,  &
                                 cloud_ice_input=cloud_ice_input, &
                                 cloudtemp=cloudtemp,    &
                                 cloudvapor=cloudvapor, kbot=kbot)  
          else   ! columns_mode, predicted cloud physics
            call define_atmos_input_fields     &
                                (is, ie, js, je, pfull, phalf, t, q,  &
                                 ts, albedo, land, Atmos_input_mdl,    &
                                 cloud_water_input=cloud_water_input,  &
                                 cloud_ice_input=cloud_ice_input, &
                                 kbot=kbot)  
          endif
        else if (present(cloud_water_input) .or.    &
                 present(cloud_ice_input) ) then
            call error_mesg ( 'radiation_driver_mod', &
              'must pass in both cloud_water and cloud_ice if one is&
	          & passed', FATAL)
        else                  !  full gcm mode or columns mode w/o
                              !  predicted cloud physics
          call define_atmos_input_fields     &
                              (is, ie, js, je, pfull, phalf, t, q,  &
                               ts, albedo, land, Atmos_input_mdl,   &
                               kbot=kbot)  
        endif
      endif

!---------------------------------------------------------------------
!    call define_radiative_gases to obtain the values of the radiatively
!    active trace gases on radiation steps or on all steps if these 
!    input fields are to be time-averaged. 
!---------------------------------------------------------------------
      if (do_rad .or. do_average_gases) then
        call define_radiative_gases (is, ie, js, je, Rad_time, lat, &
                                     Atmos_input_mdl, Time_next, &
                                     Rad_gases_mdl)
      endif
      call mpp_clock_end (misc_clock)

!--------------------------------------------------------------------
!  call aerosol_driver to define the (i,j,k) aerosol field to be
!  used for the radiation calculation.   
!--------------------------------------------------------------------
      call mpp_clock_begin (aerosol_clock)
      if (do_rad) then
        if (Rad_control%do_aerosol) then 
          call aerosol_driver (Aerosol_mdl, Rad_time,   &
                               Atmos_input_mdl%pflux, is, js)  
        endif
      endif
      call mpp_clock_end (aerosol_clock)

!--------------------------------------------------------------------
!    when using the sea-esf radiation, call clouddrvr to obtain the 
!    cloud-radiative properties needed for the radiation calculation. 
!    (these properties are obtained within radiation_calc when executing
!    the original fms radiation code). this call is made on radiation
!    steps and on all steps when these fields are to be time-averaged.
!--------------------------------------------------------------------
      call mpp_clock_begin (clouds_clock)
      if (do_rad .or. do_average_clouds) then
        if (do_sea_esf_rad) then
          if (present(kbot) ) then
            if (Environment%running_sa_model) then
              call clouddrvr (is, ie, js, je, lat, Rad_time,   &
                              Atmos_input_mdl, Astro_mdl,   &
                              Cldrad_props_mdl, Cld_diagnostics, &
                              Time_next=Time_next, kbot=kbot,  &
                              mask=mask,&
                              cloud_water_input=cloud_water_input,  &
                              cloud_ice_input=cloud_ice_input, &
                              cloud_area_input=cloud_area_input)
            else    ! (running_sa_model)
              call clouddrvr (is, ie, js, je, lat, Rad_time,   &
                              Atmos_input_mdl, Astro_mdl,   &
                              Cldrad_props_mdl, Cld_diagnostics, &
                              Time_next=Time_next, kbot=kbot, mask=mask)
            endif ! (running_sa_model)
          else    !  (present(kbot) 
            if (Environment%running_sa_model) then
              call clouddrvr (is, ie, js, je, lat, Rad_time,   &
                              Atmos_input_mdl, Astro_mdl,  &
                              Cldrad_props_mdl, Cld_diagnostics, &
                              Time_next=Time_next, &
                              cloud_water_input=cloud_water_input,  &
                              cloud_ice_input=cloud_ice_input, &
                              cloud_area_input=cloud_area_input)
            else   ! (running_sa_model)
              call clouddrvr (is, ie, js, je, lat, Rad_time,   &
                              Atmos_input_mdl, Astro_mdl,  &
                              Cldrad_props_mdl, Cld_diagnostics, &
                              Time_next=Time_next)
            endif  ! (running_sa_model)
          endif
        endif
      endif
      call mpp_clock_end (clouds_clock)

!--------------------------------------------------------------------
!    if the option to time-average the input fields used to calculate 
!    the radiative fluxes is activated (do_average = true), pass the 
!    instantaneous input fields that have been defined (Astro,   
!    Atmos_input, and optionally Rad_gases and Cldrad_props) to 
!    time_average_input_data. If this is not a radiation step, the 
!    values will be accumulated; if it is, the accumulated values will 
!    be averaged and those fields returned for input to radiation_calc.
!---------------------------------------------------------------------
      call mpp_clock_begin (misc_clock)
      if (Environment%running_gcm .and. do_average ) then
        if (do_average_gases) then
          if (do_average_clouds) then
            call time_average_input_data (is, ie, js, je,   &
                                          Atmos_input_mdl, &
                                          Astro_mdl,  &
                                          Rad_gases=Rad_gases_mdl, &
                                          Cldrad_props=Cldrad_props_mdl)
          else
            call time_average_input_data (is, ie, js, je,   &
                                          Atmos_input_mdl, &
                                          Astro_mdl,  &
                                          Rad_gases=Rad_gases_mdl)
          endif
        else
          if (do_average_clouds) then
            call time_average_input_data (is, ie, js, je,   &
                                          Atmos_input_mdl, &
                                          Astro_mdl,  &
                                          Cldrad_props=Cldrad_props_mdl)
          else
            call time_average_input_data (is, ie, js, je,   &
                                          Atmos_input_mdl, &
                                          Astro_mdl)
          endif
        endif
      endif
      call mpp_clock_end (misc_clock)

!---------------------------------------------------------------------
!    on radiation timesteps, call radiation_calc to determine new radia-
!    tive fluxes and heating rates.
!---------------------------------------------------------------------
      if (do_rad) then
        call radiation_calc (is, ie, js, je, phalf, Rad_time,  &
                             Time_next, lat, lon, Atmos_input_mdl,   &
                             Rad_gases_mdl, Aerosol_mdl,   &
                             Cldrad_props_mdl,   &
                             Cld_diagnostics, Astro_mdl, Rad_output, &
                             Lw_output, Sw_output, Fsrad_output,  &
                             mask=mask, kbot=kbot)       

      endif

!-------------------------------------------------------------------
!    on all timesteps, call update_rad_fields to update the temperature 
!    tendency and define the fluxes needed by other component models.
!    if the shortwave fluxes are to be renormalized because of the 
!    change in zenith angle since the last radiation timestep, that also
!    is done in this subroutine. 
!-------------------------------------------------------------------
      call mpp_clock_begin (misc_clock)
      if (Environment%running_gcm .or.  &
          Environment%running_sa_model) then
        call update_rad_fields (is, ie, js, je, Time_next, Astro2, &
                                Sw_output, Astro_mdl, Rad_output, tdt, &
                                coszen, flux_sw, flux_lw, flux_ratio)
!-------------------------------------------------------------------
!    call radiation_netcdf to produce radiation diagnostics, both 
!    fields and integrals.
!-------------------------------------------------------------------
        if (do_sea_esf_rad) then
          call radiation_netcdf (is, ie, js, je, Time_next, lat, &
                                 Atmos_input_mdl%tsfc, &
                                 albedo, flux_ratio, &
                                 Astro_mdl, Rad_output,  &
                                 Lw_output=Lw_output,&
                                 Sw_output=Sw_output)
        else
          call radiation_netcdf (is, ie, js, je, Time_next, lat, &
                                 Atmos_input_mdl%tsfc,  &
                                 albedo, flux_ratio, &
                                 Astro_mdl, Rad_output,  &
                                 Fsrad_output=Fsrad_output,&
                                 mask=mask)
        endif ! do_sea_esf_rad

!---------------------------------------------------------------------
!    call write_rad_output_file to produce a netcdf output file of 
!    radiation-package-relevant variables. note that this is called
!    only on radiation steps, so that the effects of sw renormalization
!    will not be seen in the variables of the data file written by
!    write_rad_output_file.
!---------------------------------------------------------------------
        if (do_rad) then
          if (do_sea_esf_rad) then
            if (Rad_control%do_aerosol) then
              call write_rad_output_file (is, ie, js, je,  &
                                          Atmos_input_mdl, Rad_output, &
                                          Sw_output, Lw_output,   &
                                          Rad_gases_mdl,   &
                                          Cldrad_props_mdl,  &
                                          Cld_diagnostics, Time_next, &
                                          Aerosol_mdl%aerosol)
            else
              call write_rad_output_file (is, ie, js, je,  &
                                          Atmos_input_mdl, Rad_output, &
                                          Sw_output, Lw_output,   &
                                          Rad_gases_mdl,   &
                                          Cldrad_props_mdl,  &
                                          Cld_diagnostics, Time_next)
            endif ! (do_aerosol)
          endif ! do_sea_esf_rad
        endif  ! do_rad

!---------------------------------------------------------------------
!    call deallocate_arrays to deallocate the array space associated 
!    with stack-resident derived-type variables.
!---------------------------------------------------------------------
        call deallocate_arrays (Rad_gases_mdl, Cldrad_props_mdl, &
                                Astro_mdl, Astro2, Lw_output, &
                                Atmos_input_mdl, Fsrad_output, &
                                Sw_output, Cld_diagnostics)
        if (do_rad) then
          if (Rad_control%do_aerosol) then
            deallocate (Aerosol_mdl%aerosol)
            deallocate (Aerosol_mdl%aerosol_optical_names)
            deallocate (Aerosol_mdl%optical_index        )
            deallocate (Aerosol_mdl%sulfate_index        )
          endif
         endif

!---------------------------------------------------------------------
!    complete radiation step when running within a gcm. update the 
!    points-processed counter. if all points in the processor's sub-
!    domain have been processed, set the counter to zero, set the 
!    radiation alarm to go off rad_time_step seconds from now, and
!    set do_rad to false, so that radiation will not be calculated until
!    the alarm goes off.
!--------------------------------------------------------------------
        num_pts = num_pts + size(t,1) * size(t,2)
        if (num_pts == total_pts) then
          num_pts = 0
          if (do_rad) then
            rad_alarm = rad_alarm + rad_time_step
            do_rad = .false.
          endif
        endif
      endif  ! (running_gcm)
      call mpp_clock_end (misc_clock)

!-----------------------------------------------------------------------


end subroutine radiation_driver



!#####################################################################


subroutine radiation_driver_end

!----------------------------------------------------------------------
!    radiation_driver_end is the destructor for radiation_driver_mod.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

      integer :: unit

!----------------------------------------------------------------------

!----------------------------------------------------------------------
!    when running in gcm, write a restart file. this is not done in the
!    standalone case.
!---------------------------------------------------------------------
      if (Environment%running_gcm) then
        unit = open_file (file='RESTART/radiation_driver.res', &
                          form='native', action='write')

!---------------------------------------------------------------------
!    only the root pe will write control information -- the last value 
!    in the list of restart versions and the alarm information.
!---------------------------------------------------------------------
        if (get_my_pe() == 0) then
!       if (get_my_pe() == get_root_pe() ) then
          write (unit) restart_versions(size(restart_versions))
          write (unit) rad_alarm, rad_time_step
        endif

!---------------------------------------------------------------------
!    write out the restart data.
!---------------------------------------------------------------------
        call write_data (unit, Rad_output%tdt_rad)
        call write_data (unit, Rad_output%tdtlw)
        call write_data (unit, Rad_output%flux_sw_surf)
        call write_data (unit, Rad_output%flux_lw_surf)
        call write_data (unit, Rad_output%coszen_angle)

!---------------------------------------------------------------------
!    write out the optional time average restart data. note that 
!    do_average and renormalize_sw_fluxes may not both be true.
!---------------------------------------------------------------------
        if (get_my_pe() == 0) then
!       if (get_my_pe() == get_root_pe() ) then
          write (unit) do_average, renormalize_sw_fluxes,   &
                       do_clear_sky_pass
        endif
        if (do_average) then
          call write_data (unit, nsum)
          call write_data (unit, psum)
          call write_data (unit, tsum)
          call write_data (unit, fsum)
          call write_data (unit, qsum)
          call write_data (unit, asum)
          call write_data (unit, csum)
          call write_data (unit, ssum)
          call write_data (unit, rrsum)

!---------------------------------------------------------------------
!    write out the control information defining whether the optional 
!    time average restart data is to be present in the file.
!---------------------------------------------------------------------
          if (get_my_pe() == 0) then
!         if (get_my_pe() == get_root_pe() ) then
            write (unit) do_average_gases, do_average_clouds
          endif

!---------------------------------------------------------------------
!    write out the radiative gas time averaged data, if that option
!    has been activated.
!---------------------------------------------------------------------
          if (do_average_gases) then
            call error_mesg ('radiation_driver_mod', &
            ' do_average_gases not yet implemented', FATAL)
!           call write_data (unit,   ), etc.
          endif

!---------------------------------------------------------------------
!    write out the cloud radiative properties time averaged data, 
!    if that option has been activated.
!---------------------------------------------------------------------
          if (do_average_clouds) then
            call error_mesg ('radiation_driver_mod', &
               ' do_average_clouds not yet implemented', FATAL)
!           call write_data (unit,   ), etc.
          endif

!---------------------------------------------------------------------
!    write out the optional shortwave renormalization data. 
!---------------------------------------------------------------------
        else if (renormalize_sw_fluxes) then   ! (do_average)
          call write_data (unit, solar_save)
          call write_data (unit, flux_sw_surf_save)
          call write_data (unit, sw_heating_save)
          call write_data (unit, tot_heating_save)
          call write_data (unit, dfsw_save) 
          call write_data (unit, ufsw_save) 
          call write_data (unit, fsw_save)  
          call write_data (unit, hsw_save)
          if (do_clear_sky_pass) then
            call write_data (unit, sw_heating_clr_save)
            call write_data (unit, tot_heating_clr_save)
            call write_data (unit, dfswcf_save) 
            call write_data (unit, ufswcf_save) 
            call write_data (unit, fswcf_save)  
            call write_data (unit, hswcf_save)
          endif
        endif    ! (do_average)

!---------------------------------------------------------------------
!    close the radiation_driver.res file
!---------------------------------------------------------------------
        call close_file (unit)
      endif   ! (running_gcm)

!---------------------------------------------------------------------
!    wrap up modules initialized by this module.
!---------------------------------------------------------------------
      call radiative_gases_end
      if (Rad_control%do_aerosol) then
        call aerosol_end
      endif
      call astronomy_end

!---------------------------------------------------------------------
!    wrap up modules specific to the radiation package in use.
!---------------------------------------------------------------------
      if (do_sea_esf_rad) then
        call cloudrad_package_end
        call rad_output_file_end
        call sea_esf_rad_end
      else
        call original_fms_rad_end
      endif

!---------------------------------------------------------------------
!    release space for allocatable data.
!---------------------------------------------------------------------
      if (Environment%running_gcm  .or. &
          Environment%running_sa_model) then
        deallocate (Rad_output%tdt_rad, Rad_output%tdtlw,  &
                    Rad_output%flux_sw_surf,  &
                    Rad_output%flux_lw_surf, Rad_output%coszen_angle, &
                    Rad_output%tdt_rad_clr, Rad_output%tdtsw, &
                    Rad_output%tdtsw_clr)

!---------------------------------------------------------------------
!    release space for renormalization arrays, if that option is active.
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes)  then
          deallocate (solar_save, flux_sw_surf_save, sw_heating_save, &
                      dfsw_save, ufsw_save, fsw_save, hsw_save, &
                      tot_heating_save)
          if (do_clear_sky_pass) then
            deallocate (sw_heating_clr_save, tot_heating_clr_save,  &
                        dfswcf_save, ufswcf_save, fswcf_save,   &
                        hswcf_save)
          endif
        endif

!---------------------------------------------------------------------
!    release space for all_step_diagnostics arrays, if that option is 
!    active.
!---------------------------------------------------------------------
        if (all_step_diagnostics)  then
          deallocate (olr_save, lwups_save, lwdns_save, tdtlw_save)
          if (do_clear_sky_pass) then
            deallocate (olr_clr_save, lwups_clr_save, lwdns_clr_save, &
                        tdtlw_clr_save)
          endif
         endif

!---------------------------------------------------------------------
!    release space for accumulation arrays, if time averaging is active.
!---------------------------------------------------------------------
        if (do_average) then
          deallocate (nsum, psum, tsum, qsum, asum, csum, ssum, fsum, &
                      rrsum)
          if (do_average_gases) then
!----------------------------------------------------------------------
!    deallocate space for any Rad_gases components which were to be 
!    time-averaged.
!----------------------------------------------------------------------
            call error_mesg ('radiation_driver_mod', &
              ' do_average_gases not yet implemented', FATAL)
          endif

!----------------------------------------------------------------------
!    deallocate space for any Cldrad_props components which were to be 
!    time-averaged.
!----------------------------------------------------------------------
          if (do_average_clouds) then
            call error_mesg ('radiation_driver_mod', &
               ' do_average_clouds not yet implemented', FATAL)
          endif
        endif     ! (do_average)
      endif    ! (running_gcm)

!----------------------------------------------------------------------
!    set initialization status flag.
!----------------------------------------------------------------------
      radiation_driver_initialized = .false.



end subroutine radiation_driver_end



!#######################################################################

subroutine get_dtrad (irad)

!--------------------------------------------------------------------
!    get_dtrad retrieves the radiation timestep (seconds) for export
!    to another module (specifically, radtest.f90, standalone driver 
!    module).
!--------------------------------------------------------------------

integer, intent(out)  :: irad

!--------------------------------------------------------------------
!
!  intent(out) variables:
!
!       irad    radiation timestep [ sec ]
!
!---------------------------------------------------------------------

      irad = rad_time_step

!---------------------------------------------------------------------


end subroutine get_dtrad 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################

subroutine allocate_average 

!--------------------------------------------------------------------
!    allocate_average allocates space for the arrays used to accumulate
!    the input fields being time-averaged, when that option is active.
!-------------------------------------------------------------------- 


!---------------------------------------------------------------------
!    allocate space for the atmospheric and astronomy input field
!    accumulation arrays.
!---------------------------------------------------------------------
      allocate (psum (id ,jd ,kmax+1), tsum (id ,jd ,kmax+1),    &
                qsum (id ,jd ,kmax), asum (id ,jd ), csum (id ,jd ), &
                fsum (id ,jd ), ssum(id, jd),   nsum (id ,jd ), &
                rrsum(id,jd) )

      if (do_average_gases) then
!----------------------------------------------------------------------
!    allocate space for any Rad_gases components which are to be 
!    time-averaged.
!----------------------------------------------------------------------
        call error_mesg ('radiation_driver_mod', &
         ' do_average_gases not yet implemented', FATAL)
      endif

!----------------------------------------------------------------------
!    allocate space for any Cldrad_props components which are to be 
!    time-averaged.
!----------------------------------------------------------------------
      if (do_average_clouds) then
        call error_mesg ('radiation_driver_mod', &
         ' do_average_clouds not yet implemented', FATAL)
      endif

!---------------------------------------------------------------------

end subroutine allocate_average 


!#####################################################################

subroutine initialize_average 

!---------------------------------------------------------------------
!    initialize_average initializes the accumulation arrays for time-
!    averaged input fields when that option is active.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    initialize the astronomical and atmospheric input data accumulation
!    arrays.
!---------------------------------------------------------------------
      nsum =  0 
      psum =  0. 
      tsum =  0. 
      qsum =  0. 
      asum =  0.  
      csum =  0. 
      ssum =  0.
      fsum =  0. 
      rrsum = 0.

!---------------------------------------------------------------------
!    include code here to initialize any Rad_gases components
!    which are to be time-averaged.
!---------------------------------------------------------------------
      if (do_average_gases) then
        call error_mesg ('radiation_driver_mod', &
         ' do_average_gases not yet implemented', FATAL)
      endif

!---------------------------------------------------------------------
!    include code here to initialize Cldrad_props components which are 
!    to be time-averaged.
!---------------------------------------------------------------------
      if (do_average_clouds) then
        call error_mesg ('radiation_driver_mod', &
         ' do_average_clouds not yet implemented', FATAL)
      endif

!----------------------------------------------------------------------


end subroutine initialize_average



!#####################################################################

subroutine read_restart_file 

!-------------------------------------------------------------------
!    read_restart_file reads a restart file containing radiation
!    restart information. it may be either a radiation_driver.res, or
!    an older sea_esf_rad.res file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables

      integer               :: unit
      logical               :: end
      logical               :: avg_present, renorm_present,  &
                               cldfree_present
      character(len=4)      :: chvers
      logical               :: avg_gases, avg_clouds
      integer, dimension(5) :: null
      integer               :: vers
      integer               :: new_rad_time, old_time_step

      real, dimension(size(gas_component_sum,1),      &
                      size(gas_component_sum,2),      &
                      size(gas_component_sum,3))  ::  &
		                                dum_gas_component_sum
      real, dimension(size(cld_component_sum,1),      &
                      size(cld_component_sum,2),      &
                      size(cld_component_sum,3))  ::  &
		                                dum_cld_component_sum

!---------------------------------------------------------------------
!    if one is using the sea_esf_rad package and there is a 
!    sea_esf_rad.res restart file present in the input directory, it 
!    must be version 1 in order to be readable by the current module.
!    otherwise, exit.
!---------------------------------------------------------------------
      if (do_sea_esf_rad .and. file_exist('INPUT/sea_esf_rad.res')) then
        unit = open_file (file='INPUT/sea_esf_rad.res',  &
                          form='native', action='read')
        read (unit) vers
        if ( vers /= 1 ) then
          write (chvers,'(i4)') vers
          call error_mesg ('radiation_driver_mod', &
             'restart version '//chvers//' cannot be read &
              &by this module version', FATAL)
        endif

!---------------------------------------------------------------------
!    if a radiation_driver.res file is present, then it must be one
!    of the versions listed as readable for the radiation package
!    being employed. for the original_fms package, these versions are
!    found in the array restart_versions, while for the sea_esf_rad
!    package, they are found in restart_versions_sea. if the version
!    is not acceptable, exit.
!---------------------------------------------------------------------
      else if ( file_exist('INPUT/radiation_driver.res')  )  then   
        unit = open_file (file='INPUT/radiation_driver.res',  &
                          form='native', action='read')
        read (unit) vers
        if (do_sea_esf_rad) then
          if ( .not. any(vers == restart_versions_sea) ) then
            write (chvers,'(i4)') vers
            call error_mesg ('radiation_driver_mod', &
                    'restart version '//chvers//' cannot be read &
                     &by this module version when using sea_esf_rad &
		     &radiation_package', FATAL)
          endif
        else
          if ( .not. any(vers == restart_versions) ) then
            write (chvers,'(i4)') vers
            call error_mesg ('radiation_driver_mod', &
                    'restart version '//chvers//' cannot be read &
                     &by this module version when using the original_&
		     &fms radiation package', FATAL)
          endif
        endif
      endif
!-----------------------------------------------------------------------
!    read alarm information.  if reading an sea_esf_rad.res file 
!    (version 1), recover the time step previously used, and set the
!    radiation alarm to be 1 second from now, assuring radiation 
!    recalculation on the first model step of this run. for later 
!    restarts, read the previous radiation timestep and the rad_alarm
!    that was present when the restart was written.
!-----------------------------------------------------------------------
      if (vers == 1) then
        read (unit) null 
        old_time_step = seconds_per_day*null(4) + null(3)
        rad_alarm = 1        
      else
        read (unit) rad_alarm, old_time_step    
      endif

!---------------------------------------------------------------------
!    read the radiation restart data. it consists of radiative temper-
!    ature tendencies, sw surface fluxes, lw surface fluxes and the
!    value of the cosine of the zenith angle to be used for the next
!    ocean albedo calcuation, in restart versions after version 1. for
!    restart version 1, set the cosine of the zenith angle to the value
!    used on initialization for use in diagnostics. since in this case 
!    rad_alarm has been set so that radiation is called on the next 
!    step, the proper zenith angles will be calculated and then used to
!    define the albedos.
!---------------------------------------------------------------------
      call read_data (unit, Rad_output%tdt_rad )
      if (vers >= 4) then
        call read_data (unit, Rad_output%tdtlw )
      else          
        Rad_output%tdtlw = 0.0
      endif
      call read_data (unit, Rad_output%flux_sw_surf )
      call read_data (unit, Rad_output%flux_lw_surf )
      if (vers /= 1) then
        call read_data (unit, Rad_output%coszen_angle)
      else
        Rad_output%coszen_angle = coszen_angle_init
      endif

      if (vers <= 2) then
!---------------------------------------------------------------------
!    if averaged input data is desired, determine if accumulation arrays
!    are present in the restart file. 
!---------------------------------------------------------------------
        if (do_average) then
          call read_data (unit, nsum , end)

!---------------------------------------------------------------------
!    if accumulation arrays are present, read them.
!---------------------------------------------------------------------
          if ( .not. end) then      
!---------------------------------------------------------------------
!    read in the accumulated sums. different sums are present dependent
!    on restart version.
!---------------------------------------------------------------------
            if (vers == 1) then
              call read_data (unit, psum )
              call read_data (unit, tsum(:,:,1:kmax) )
              call read_data (unit, tsum(:,:,kmax+1) )
              call read_data (unit, fsum )
              call read_data (unit, qsum )
              call read_data (unit, asum )
              call read_data (unit, csum )
              call read_data (unit, rrsum)

!----------------------------------------------------------------------
!    if the restart was written by a job which was using the original
!    fms radiation package (rrsum == 0.) and the sea_esf_rad package is
!    now being used, reinitialize the sums since all needed variables 
!    are not available.
!----------------------------------------------------------------------
              if (rrsum(1,1) == 0 .and. do_sea_esf_rad) then
                if (get_my_pe() == 0) then
                call error_mesg ( 'radiation_driver_mod', &  
                  'cannot switch from original to seaesf radiation &
                  & while averaging radiation quantities --  different &
                  & variables are averaged. sums being reinitialized. &
	          &', NOTE )
                endif
                call initialize_average

!----------------------------------------------------------------------
!    likewise, if the restart was written by a job which was using the 
!    sea_esf_rad package (rrsum /= 0.0) and the original_fms radiation 
!    package is now being used, reinitialize the sums since different 
!    variables are needed by the two packages.
!----------------------------------------------------------------------
              else if (rrsum(1,1) /= 0 .and.   &
                       .not. (do_sea_esf_rad) ) then
                if (get_my_pe() == 0) then
                call error_mesg ( 'radiation_driver_mod', &
                   'cannot switch from seaesf to original radiation &
                   & while averaging radiation quantities --  different&
	           & variables are averaged. sums being reinitialized. &
	           &', NOTE )
                call initialize_average
		endif
              endif

!----------------------------------------------------------------------
!    version 2 is a radiation_driver.res file, written by the original
!    fms radiation package. sums from it may be used only if the
!    original_fms radiation package is being run here.
!----------------------------------------------------------------------
            else if (vers == 2 ) then
              if (.not. do_sea_esf_rad) then
                call read_data (unit, psum )
                call read_data (unit, tsum )
                call read_data (unit, qsum )
                call read_data (unit, asum )
                call read_data (unit, csum )
                call read_data (unit, ssum)
              else
!---------------------------------------------------------------------
!    a version 2 radiation_driver.res file cannot be used to supply 
!    sums when running the sea_esf_rad package, since several needed
!    variables will be missing. in this case reinitialize all of the 
!    sums.
!----------------------------------------------------------------------
                if (get_my_pe() == 0) then
                call error_mesg ( 'radiation_driver_mod', &
                   ' cannot use sums from a version 2  &
		   &radiation_driver.res file with sea_esf_rad package.&
	  	   & sums reinitialized.', NOTE)
                call initialize_average
                endif       
              endif
            endif ! (vers = 1)
          else    ! (.not. end)
            if (get_my_pe() == 0) then
              call error_mesg ( 'radiation_driver_mod', &
              ' restart file has no accumulation sums -- sums are &
	        &being initialized.', NOTE)
            endif
          endif   ! (.not. end)
        endif    ! (do_average)

!----------------------------------------------------------------------
!    version 3 is the current radiation_driver.res file, written by both
!    the original_fms radiation package and the sea_esf_rad package. 
!----------------------------------------------------------------------
      else if (vers >= 3) then
!---------------------------------------------------------------------
!    determine if accumulation arrays are present in the restart file.
!    if input fields are to be time-averaged, read the values from the
!    files.  note that avg_present and renorm_present cannot both be 
!    true.
!---------------------------------------------------------------------
        read (unit) avg_present, renorm_present, cldfree_present
        if (do_average) then
          if (avg_present) then
            call read_data (unit, nsum )
            call read_data (unit, psum )
            call read_data (unit, tsum )
            call read_data (unit, fsum )
            call read_data (unit, qsum )
            call read_data (unit, asum )
            call read_data (unit, csum )
            call read_data (unit, ssum)
            call read_data (unit, rrsum)

!---------------------------------------------------------------------
!    read logicals to indicate if the optional accumulation arrays
!    are present. if avg_gases is .true., then radiative gas accumul-
!    ation arrays are present. if avg_clouds is .true., then cloud
!    radiative property accumulation arrays are present. read these
!    arrays if averaging of these fields is desired in the current run.
!---------------------------------------------------------------------
            read (unit) avg_gases, avg_clouds
            if (avg_gases) then
              if (do_average_gases) then 
                call error_mesg ('radiation_driver_mod', &
                      ' do_average_gases not yet implemented', FATAL)
!               call read_data (unit, gas_component_sum) 
!                etc. (as many as needed)
              else
                call error_mesg ('radiation_driver_mod', &
                      ' do_average_gases not yet implemented', FATAL)
!               call read_data (unit, dum_gas_component_sum)
!                etc. (as many as needed)
              endif
            endif

            if (avg_clouds) then
              if (do_average_clouds) then
                call error_mesg ('radiation_driver_mod', &
                      ' do_average_clouds not yet implemented', FATAL)
!               call read_data (unit, cld_component_sum)
!                etc. (as many as needed)
              else
                call error_mesg ('radiation_driver_mod', &
                    ' do_average_clouds not yet implemented', FATAL)
!               call read_data (unit, dum_cld_component_sum)
!                etc. (as many as needed)
              endif
            endif

!---------------------------------------------------------------------
!    if sums are not present, send message and reinitialize them.
!---------------------------------------------------------------------
          else   ! (avg_present)
            if (get_my_pe() == 0) then
            call error_mesg ( 'radiation_driver_mod', &
            ' restart file has no accumulation sums -- sums are &
	      &being initialized.', NOTE)
            endif
          endif ! (avg_present)
  
!---------------------------------------------------------------------
!    if renormalize_sw_fluxes is true and the data is present in the
!    restart file, read it.
!---------------------------------------------------------------------
        else if (renormalize_sw_fluxes) then  ! (do_average)
          if (renorm_present) then     
            call read_data (unit, solar_save)
            call read_data (unit, flux_sw_surf_save)
            call read_data (unit, sw_heating_save)
            call read_data (unit, tot_heating_save)
            call read_data (unit, dfsw_save) 
            call read_data (unit, ufsw_save)  
            call read_data (unit, fsw_save)  
            call read_data (unit, hsw_save)

!---------------------------------------------------------------------
!    if cldfree data is desired and the data is present in the
!    restart file, read it.
!---------------------------------------------------------------------
            if (do_clear_sky_pass) then
              if (cldfree_present) then
                call read_data (unit, sw_heating_clr_save)
                call read_data (unit, tot_heating_clr_save)
                call read_data (unit, dfswcf_save) 
                call read_data (unit, ufswcf_save) 
                call read_data (unit, fswcf_save)  
                call read_data (unit, hswcf_save)

!--------------------------------------------------------------------
!    if cldfree data is desired and the data is not present in the
!    restart file, force a radiation call on next model step.  
!---------------------------------------------------------------------
              else
	        rad_alarm = 1
	      endif  ! (cldfree_present)
            endif
!---------------------------------------------------------------------
!    if renormalize_sw_fluxes is true and the data is not present in the
!    restart file, force a radiation call on next model step.  
!---------------------------------------------------------------------
          else  ! (renorm_present)
            rad_alarm = 1
          endif ! (renorm_present)
        endif ! (do_average) 
      endif    ! (vers <= 2)

!--------------------------------------------------------------------
!    close the unit used to read the .res file.
!--------------------------------------------------------------------
      call close_file (unit)

!----------------------------------------------------------------------
!    set rad_alarm to 1 to force radiation call to generate lw diagnos-
!    tcs if all_step_diagnostics is active.
!----------------------------------------------------------------------
      if (rad_alarm /= 1 .and. all_step_diagnostics) then
        rad_alarm = 1
      endif

!----------------------------------------------------------------------
!    adjust radiation alarm if radiation step has changed from restart 
!    file value, if it has not already been set to the first step.
!----------------------------------------------------------------------
      if (rad_alarm == 1) then
        if (get_my_pe() == 0)  then
!       if (get_my_pe() == get_root_pe() ) then
          call error_mesg ('radiation_driver_mod',          &
               'radiation will be called on first step of run', NOTE)
        endif
      else
        if (rad_time_step /= old_time_step ) then
          new_rad_time = rad_alarm - old_time_step + rad_time_step
          if ( new_rad_time > 0 ) then
            if (get_my_pe() == 0)  then  
!           if (get_my_pe() == get_root_pe() ) then                  
              call error_mesg ('radiation_driver_mod',          &
                  'radiation time step has changed, therefore &
                  &next radiation time also changed', NOTE)
            endif
            rad_alarm = new_rad_time
          endif
        endif  
      endif   ! (rad_alarm == 1)

!--------------------------------------------------------------------


end subroutine read_restart_file



!#####################################################################

subroutine initialize_diagnostic_integrals

!---------------------------------------------------------------------
!    initialize_diagnostic_integrals registers the desired integrals 
!    with diag_integral_mod.
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    initialize standard global quantities for integral package. 
!----------------------------------------------------------------------
      call diag_integral_field_init ('olr',    std_digits)
      call diag_integral_field_init ('abs_sw', std_digits)

!----------------------------------------------------------------------
!    if hemispheric integrals and global integrals with extended signif-
!    icance are desired, inform diag_integrals_mod.
!----------------------------------------------------------------------
      if (calc_hemi_integrals) then
        call diag_integral_field_init ('sntop_tot_sh', extra_digits)
        call diag_integral_field_init ('lwtop_tot_sh', extra_digits)
        call diag_integral_field_init ('sngrd_tot_sh', extra_digits)
        call diag_integral_field_init ('lwgrd_tot_sh', extra_digits)
        call diag_integral_field_init ('sntop_tot_nh', extra_digits)
        call diag_integral_field_init ('lwtop_tot_nh', extra_digits)
        call diag_integral_field_init ('sngrd_tot_nh', extra_digits)
        call diag_integral_field_init ('lwgrd_tot_nh', extra_digits)
        call diag_integral_field_init ('sntop_tot_gl', extra_digits)
        call diag_integral_field_init ('lwtop_tot_gl', extra_digits)
        call diag_integral_field_init ('sngrd_tot_gl', extra_digits)
        call diag_integral_field_init ('lwgrd_tot_gl', extra_digits)

!---------------------------------------------------------------------
!    if clear-sky integrals are desired, include them.
!---------------------------------------------------------------------
        if (do_clear_sky_pass) then
          call diag_integral_field_init ('sntop_clr_sh', extra_digits)
          call diag_integral_field_init ('lwtop_clr_sh', extra_digits)
          call diag_integral_field_init ('sngrd_clr_sh', extra_digits)
          call diag_integral_field_init ('lwgrd_clr_sh', extra_digits)
          call diag_integral_field_init ('sntop_clr_nh', extra_digits)
          call diag_integral_field_init ('lwtop_clr_nh', extra_digits)
          call diag_integral_field_init ('sngrd_clr_nh', extra_digits)
          call diag_integral_field_init ('lwgrd_clr_nh', extra_digits)
          call diag_integral_field_init ('sntop_clr_gl', extra_digits)
          call diag_integral_field_init ('lwtop_clr_gl', extra_digits)
          call diag_integral_field_init ('sngrd_clr_gl', extra_digits)
          call diag_integral_field_init ('lwgrd_clr_gl', extra_digits)
        endif
      endif

!--------------------------------------------------------------------


end subroutine initialize_diagnostic_integrals



!#######################################################################

subroutine diag_field_init ( Time, axes )

!---------------------------------------------------------------------
!    diag_field_init registers the desired diagnostic fields with the
!    diagnostics manager.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time
integer        , intent(in) :: axes(4)

!--------------------------------------------------------------------
!  intent(in) variables
!
!      Time        current time
!      axes        data axes for use with diagnostic fields
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables


      character(len=8)  ::   clr
      character(len=16) ::   clr2
      integer           ::   i, n

!---------------------------------------------------------------------
!    determine how many passes are needed through the name generation 
!    loop. 
!---------------------------------------------------------------------
      if (do_clear_sky_pass) then
        n= 2
      else
        n= 1
      endif

!---------------------------------------------------------------------
!    generate names for standard and clear sky diagnostic fields. if 
!    clear sky values being generated, generate the clear sky names
!    on pass 1, followed by the standard names.
!---------------------------------------------------------------------
      do i = 1, n
        if ( i == n) then
          clr  = "    "
          clr2 = "          "
        else
          clr  = "_clr"
          clr2 = "clear sky "
        endif

        id_tdt_sw(i) = register_diag_field (mod_name,   &
                'tdt_sw'//trim(clr), axes(1:3), Time, & 
                trim(clr2)//'temperature tendency for SW radiation', &
                'deg_K/sec', missing_value=missing_value) 

        id_tdt_lw(i) = register_diag_field (mod_name,    &
                'tdt_lw'//trim(clr), axes(1:3), Time, &
                trim(clr2)//'temperature tendency for LW radiation', &
                'deg_K/sec', missing_value=missing_value)

        id_swdn_toa(i) = register_diag_field (mod_name,   &
                'swdn_toa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux down at TOA', &
                'watts/m2', missing_value=missing_value)

        id_swup_toa(i) = register_diag_field (mod_name,    &
                'swup_toa'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux up at TOA', &
                'watts/m2', missing_value=missing_value)

        id_olr(i) = register_diag_field (mod_name,   &
                'olr'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'outgoing longwave radiation', &
                'watts/m2', missing_value=missing_value)

        id_swup_sfc(i) = register_diag_field (mod_name,    &
                'swup_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux up at surface', &
                'watts/m2', missing_value=missing_value)

        id_swdn_sfc(i) = register_diag_field (mod_name,     &
                'swdn_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'SW flux down at surface', &
                'watts/m2', missing_value=missing_value)

        id_lwup_sfc(i) = register_diag_field (mod_name,   &
                'lwup_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'LW flux up at surface', &
                'watts/m2', missing_value=missing_value)

        id_lwdn_sfc(i) = register_diag_field (mod_name,    &
                'lwdn_sfc'//trim(clr), axes(1:2), Time, &
                trim(clr2)//'LW flux down at surface', &
                'watts/m2', missing_value=missing_value)

      end do

!----------------------------------------------------------------------
!    register fields that are not clear-sky depedent.
!----------------------------------------------------------------------
      id_alb_sfc = register_diag_field (mod_name,    &
                'alb_sfc', axes(1:2), Time, &
!               'surface albedo', 'percent')
! BUGFIX
                'surface albedo', 'percent', &
                  missing_value=missing_value) 


      id_cosz = register_diag_field (mod_name,    &
                'cosz',axes(1:2),  Time,    &
                'cosine of zenith angle',    &
                'none', missing_value=missing_value)

      id_fracday = register_diag_field (mod_name,   &
                'fracday',axes(1:2), Time,   &
                'daylight fraction of radiation timestep',   &
                'percent', missing_value=missing_value)

!-----------------------------------------------------------------------


end subroutine diag_field_init



!#####################################################################

subroutine define_rad_times (Time, Time_next, do_rad, Rad_time,   &
                             Dt_zen, dt)

!--------------------------------------------------------------------
!    subroutine define_rad_times determines whether radiation is to be 
!    calculated on the current timestep, and defines various 
!    time-related variables.
!-------------------------------------------------------------------- 

!---------------------------------------------------------------------
type(time_type), intent(in)   :: Time, Time_next
logical,         intent(out)  :: do_rad
type(time_type), intent(out)  :: Rad_time, Dt_zen
integer,         intent(out)  :: dt
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     Time         current model time  
!                  [ time_type, days and seconds]
!     Time_next    model time on the next atmospheric timestep
!                  [ time_type, days and seconds]
!     
!   intent(out) variables:
!
!     Rad_time     time at which the climatologically-determined, time-
!                  varying input fields to radiation should apply    
!                  [ time_type, days and seconds]
!     Dt_zen       time interval over which  cosine of zenith angle
!                  is to be averaged when diurnally-varying solar
!                  radiation is activated
!                  [ time_type, days and seconds ]
!     do_rad       logical indicating whether radiation is to be cal-
!                  culated on the current timestep
!     dt           physics time step (frequency of calling 
!                  radiation_driver)
!                  [ seconds ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer        :: day, sec
      integer        :: dum, tod(3)

!---------------------------------------------------------------------
!   local variables:
!
!      day            day component of atmospheric timestep
!                     [ days ]
!      sec            seconds component of atmospheric timestep
!                     [ seconds ]
!      dum            dummy variable
!      tod            hours, minutes and seconds components of current
!                     time
!                     [ hours, minutes, seconds ]
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    compute the atmospheric timestep by differencing the time on the
!    next step with the current time. if dt .le. 0, call error exit.
!--------------------------------------------------------------------
      call get_time (Time_next-Time, sec, day)
      dt = day*seconds_per_day + sec
      if (dt <= 0) call error_mesg ('radiation_driver_mod', &
                                    'Time_next <= Time', FATAL)


!-------------------------------------------------------------------
!    for the standalone case, new radiation outputs are calculated on 
!    every step, using climatological variable values at the time spec-
!    ified by the input argument Time, with diurnally varying solar
!    radiation using the cosine of the zenith angle averaged over the
!    namelist-specified radiation time step (Dt_zen). 
!-------------------------------------------------------------------
      if (Environment%running_standalone  .or.  &
          Environment%running_sa_model) then
        do_rad = .true.
        Rad_time = Time
        Dt_zen = set_time(rad_time_step, 0)

!--------------------------------------------------------------------
!    if running a gcm aplication, if this is the first call by this
!    processor on this time step to radiation_driver (i.e. num_pts = 0),
!    determine if this is a radiation time step by decrementing the time
!    to alarm by the current model timestep.  if the alarm "goes off", 
!    i.e., is .le. 0, set do_rad to true, indicating this is a radiation
!    step. otherwise set it to .false. . 
!--------------------------------------------------------------------
      else
        if (num_pts == 0)  then
          rad_alarm = rad_alarm -  dt
        endif
        if (rad_alarm <= 0) then
          do_rad = .true.
        else
          do_rad = .false.
        endif

!-------------------------------------------------------------------
!    define the time to be used in defining the time-varying input 
!    fields for the radiation calculation (Rad_time). 
!-------------------------------------------------------------------
        if (rsd) then

!--------------------------------------------------------------------
!    if this is a repeat-same-day (rsd) experiment, define Rad_time
!    as the specified year-month-day (rad_date(1:3)), and the 
!    hr-min-sec of the current time (Time).
!---------------------------------------------------------------------
          if (.not. use_rad_date)   &
            call error_mesg ('radiation_driver_mod', &  
              'if (rsd), must set rad_date(1:3) to valid date', FATAL)
            call get_date (Time, dum, dum, dum, tod(1), tod(2), tod(3))
            Rad_time = set_date (rad_date(1), rad_date(2),& 
                                 rad_date(3), tod(1), tod(2), &
                                 tod(3))

!---------------------------------------------------------------------
!    if the specified date option is active, define Rad_time to be that
!    date and time.
!----------------------------------------------------------------------
        else if (use_rad_date) then
          Rad_time = set_date (rad_date(1), rad_date(2), rad_date(3),  &
                               rad_date(4), rad_date(5), rad_date(6))

!---------------------------------------------------------------------
!    if neither of these special cases is active, define Rad_time as
!    the current time (Time).
!---------------------------------------------------------------------
        else
          Rad_time = Time
        endif  ! (rsd)

!--------------------------------------------------------------------
!    define the averaging period over which diurnally-varying radiation
!    is to apply. if input field averaging is active, it is to apply 
!    over a model timestep, since the astronomy is calculated every step
!    and will be averaged when the radiation step is reached. otherwise,
!    on radiation steps, the period is the radiation timestep.  if 
!    renormalization is active and it is not a radiation step, calculate
!    values relevant for the model time step, so that the values calcul-
!    ated on the radiation step may be renormalized..
!---------------------------------------------------------------------
        if (do_diurnal) then
          if ( (do_average)  .or.    &
               (renormalize_sw_fluxes .and. .not. do_rad)) then
            Dt_zen = set_time (dt,0)
          else if (do_rad) then
            Dt_zen = set_time (rad_time_step,0)
          endif
        endif
      endif  ! (running_standalone)

!--------------------------------------------------------------------



end subroutine define_rad_times




!######################################################################

subroutine obtain_astronomy_variables (is, ie, js, je, dt, lat, lon, &
                                       Rad_time, Dt_zen, Astro, Astro2,&
                                       Rad_output )

!---------------------------------------------------------------------
!    obtain_astronomy_variables retrieves astronomical variables 
!    at the desired locations and time, valid over the requested time 
!    intervals.
!---------------------------------------------------------------------
integer,                     intent(in)    ::  is, ie, js, je, dt
real, dimension(:,:),        intent(in)    ::  lat, lon
type(time_type),             intent(in)    ::  Rad_time, Dt_zen
type(astronomy_type),        intent(inout)   ::  Astro2
type(astronomy_type),        intent(inout) ::  Astro
type(rad_output_type),       intent(inout)   ::  Rad_output

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      dt           atmospheric time step (frequency of calling 
!                   radiation_driver)
!                   [ seconds ]
!      Rad_time     time at which the climatologically-determined, 
!                   time-varying input fields to radiation should apply
!                   [ time_type, days and seconds]
!      Dt_zen       time interval over which cosine of zenith angle
!                   is to be averaged when diurnally-varying solar
!                   radiation is activated
!                   [ time_type, days and seconds ]
!      lat          latitude of model points  
!                   [ radians ]
!      lon          longitude of model points 
!                   [ radians ]
!
!   intent(out) variables:
!
!      Astro         astronomy_type structure; contains the following
!                    components defined in this suboutine, valid over
!                    either the radiation or model timestep:
!         solar         shortwave flux factor: cosine of zenith angle *
!                       daylight fraction / (earth-sun distance squared)
!                       [ non-dimensional ]
!         cosz          cosine of zenith angle --  mean value over
!                       appropriate averaging interval
!                       [ non-dimensional ]
!         fracday       fraction of timestep during which the sun is 
!                       shining
!                       [ non-dimensional ]
!         rrsun         inverse of square of earth-sun distance, 
!                       relative to the mean earth-sun distance
!                       [ non-dimensional ]
!
!      Rad_output    rad_output_type structure, contains the 
!                    following components defined in this subroutine
!         
!         coszen_angle  cosine of zenith angle appropriate for use
!                       in defining ocean albedo -- for diurnally 
!                       varying solar input, this is the value 
!                       obtained one radiation timestep from now, 
!                       for other solar schemes, it is the zenith 
!                       angle as calculated using the current Rad_time.
!
!
!   intent(inout) variables:
!
!      Astro2        same as Astro, but contains values valid for 
!                    current physics timestep (when 
!                    renormalize_sw_fluxes is true)
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables

      type(astronomy_type)   :: Astro1
      type(time_type)        :: Dt_zen2

!--------------------------------------------------------------------
!  local variables
!
!     Astro1        astronomy_type variable with values of its elements
!                   applicable at the appropriate time for use in
!                   calculating the ocean albedo 
!            solar1        value of solar at next timestep
!                          [ non-dimensional ]
!            cosz1         value of cosz on next timestep
!                          [ non-dimensional ]
!            fracday1      value of fracday on next timestep
!                          [ non-dimensional ]
!            rrsun1        value of rrsun on next timestep
!                          [ non-dimensional ]
!
!     Dt_zen2       time-type variable containing the components of the
!                   physics time step, needed when renormalize_sw_fluxes
!                   is true
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    define the astronomical inputs to radiation  (cosine of zenith 
!    angle, daylight fraction, earth-sun distance and solar flux factor)
!    as elements of an astronomy_type variable. if using diurnally-
!    varying sw radiation, the values are applicable over the time 
!    interval from Rad_time to Rad_time + Dt_zen.
!---------------------------------------------------------------------
      allocate ( Astro%fracday(size(lat,1), size(lat,2) ) )
      allocate ( Astro%cosz   (size(lat,1), size(lat,2) ) )
      allocate ( Astro%solar  (size(lat,1), size(lat,2) ) )
      if (Astro%do_diurnal) then
        call diurnal_solar (lat, lon, Rad_time, Astro%cosz,   &
                            Astro%fracday, Astro%rrsun, dt_time=Dt_zen)
        Astro%fracday = MIN (Astro%fracday, 1.00)
        Astro%solar (:,:) = Astro%cosz(:,:  )*Astro%fracday(:,:)*  &
                            Astro%rrsun
      else if (Astro%do_annual) then
        call annual_mean_solar (js, je, lat, Astro%cosz, Astro%solar, &
                                Astro%fracday, Astro%rrsun)
      else if (Astro%do_daily_mean) then
        call daily_mean_solar (lat, Rad_time, Astro%cosz,  &
                               Astro%fracday, Astro%rrsun)
        Astro%solar = Astro%cosz*Astro%rrsun*Astro%fracday
      else
        call error_mesg('radiation_driver_mod', &
             ' no valid zenith angle specification', FATAL)
      endif

!---------------------------------------------------------------------
!    if using diurnally-varying sw radiation, allocate the components 
!    of an astronomy_type variable (Astro1) that will contain applicable
!    values one time step from now (either radiation or physics).
!---------------------------------------------------------------------
      if (Environment%running_gcm .or.   &
          Environment%running_sa_model) then
        if (do_diurnal) then
          allocate ( Astro1%fracday(size(lat,1), size(lat,2) ) )
          allocate ( Astro1%cosz   (size(lat,1), size(lat,2) ) )
          allocate ( Astro1%solar  (size(lat,1), size(lat,2) ) )

!---------------------------------------------------------------------
!    if renormalizing sw fluxes and using diurnally-varying sw radia-
!    tion, call diurnal_solar on radiation timesteps to obtain 
!    the astronomical inputs to radiation, applicable over the current
!    physics timestep. (Astro has the values over the radiation time-
!    step.) save these values in Astro2.
!---------------------------------------------------------------------
          if (renormalize_sw_fluxes ) then
            Dt_zen2 = set_time(dt,0)
            if (do_rad) then
              allocate ( Astro2%fracday(size(lat,1), size(lat,2) ) )
              allocate ( Astro2%cosz   (size(lat,1), size(lat,2) ) )
              allocate ( Astro2%solar  (size(lat,1), size(lat,2) ) )
              Astro2%do_diurnal    = Astro%do_diurnal
              Astro2%do_daily_mean = Astro%do_daily_mean
              Astro2%do_annual     = Astro%do_annual
              call diurnal_solar (lat, lon, Rad_time, Astro2%cosz,   &
                                  Astro2%fracday, Astro2%rrsun,   &
                                  dt_time=Dt_zen2)
              Astro2%fracday = MIN (Astro2%fracday, 1.00)
              Astro2%solar (:,:) = Astro2%cosz(:,:  )*   &
                                   Astro2%fracday(:,:)*Astro2%rrsun
            endif  ! (do_rad)

!---------------------------------------------------------------------
!    define cosine of zenith angle valid one physics time step from now
!    which will be used to define ocean albedo.
!----------------------------------------------------------------------
!RSHALB3
         if (do_rad) then
            Astro1%do_diurnal     = Astro%do_diurnal
            Astro1%do_daily_mean  = Astro%do_daily_mean
            Astro1%do_annual      = Astro%do_annual
!RSHALB1    call diurnal_solar (lat, lon, Rad_time+Dt_zen2,    &
            call diurnal_solar (lat, lon, Rad_time+Dt_zen,    &
                                Astro1%cosz, Astro1%fracday,   &
!RSHALB1                        Astro1%rrsun, dt_time=Dt_zen2)
                                Astro1%rrsun, dt_time=Dt_zen)
            Astro1%fracday = MIN (Astro1%fracday, 1.00)
            Astro1%solar (:,:) = Astro1%cosz(:,:)*Astro1%fracday(:,:)* &
                                 Astro1%rrsun
            Rad_output%coszen_angle(is:ie,js:je) =      &
!RSHALB2                           Astro1%cosz(:,:)*Astro1%fracday(:,:)
                                   Astro1%cosz(:,:)

!RSHALB3
         endif

!---------------------------------------------------------------------
!    if using diurnally_varying radiation without renormalizing fluxes,
!    call diurnal_solar again to obtain the cosine of zenith angle 
!    on the next timestep (either model timestep or radiation timestep) 
!    to be used in calculating the ocean albedo. 
!---------------------------------------------------------------------
          else   ! (renormalize)
            Astro1%do_diurnal    = Astro%do_diurnal
            Astro1%do_daily_mean = Astro%do_daily_mean
            Astro1%do_annual     = Astro%do_annual
            call diurnal_solar (lat, lon, Rad_time+Dt_zen, Astro1%cosz,&
                                Astro1%fracday, Astro1%rrsun,    &
                                dt_time=Dt_zen)
            Astro1%fracday = MIN (Astro1%fracday, 1.00)
            Astro1%solar (:,:) = Astro1%cosz(:,:)*Astro1%fracday(:,:)* &
                                 Astro1%rrsun
!RSHALB2    Rad_output%coszen_angle(is:ie,js:je) = Astro1%cosz(:,:)*  &
            Rad_output%coszen_angle(is:ie,js:je) = Astro1%cosz(:,:)
!RSHALB2                                             Astro1%fracday(:,:)
          endif  ! (renormalize_sw_fluxes)
          deallocate (Astro1%cosz)
          deallocate (Astro1%solar)
          deallocate (Astro1%fracday)

!---------------------------------------------------------------------
!    if this is not diurnally-varying radiation (i.e., it is either
!    daily-mean or annual-mean radiation), save the cosine of zenith 
!    angle on the current step (which is also valid for the next step)
!    to be used to calculate ocean albedo.
!---------------------------------------------------------------------
        else   ! (do_diurnal)
          Rad_output%coszen_angle(is:ie,js:je) = Astro%cosz (:,:)
        endif  ! (do_diurnal)
      endif  ! (running_gcm)

!-------------------------------------------------------------------



end subroutine obtain_astronomy_variables 



!####################################################################

subroutine define_atmos_input_fields (is, ie, js, je, pfull, phalf, &
                                      t, q, ts, albedo, land,  &
                                      Atmos_input, cloud_water_input, &
                                      cloud_ice_input, cloudtemp,  &
                                      cloudvapor, kbot)     

!---------------------------------------------------------------------
!    define_atmos_input_fields converts the model-supplied fields 
!    (pfull, phalf, t, q, ts, albedo, land, optionally cloud_water_input
!    and cloud_ice_input) to the form needed by the radiation modules, 
!    and on radiation timesteps returns radiation-ready fields of pres-
!    sure (press, psfc), temperature (temp, tsfc), water vapor mixing 
!    ratio (rh2o), surface albedo (asfc), land fraction (land), option-
!    ally cloud_water (cloud_water) and cloud_ice (cloud_ice) in the 
!    derived type structure Atmos_input.
!---------------------------------------------------------------------
     
integer,                 intent(in)              :: is, ie, js, je
real, dimension(:,:,:),  intent(in)              :: pfull, phalf, t, q
real, dimension(:,:),    intent(in)              :: ts, albedo, land 
type(atmos_input_type),  intent(inout)             :: Atmos_input
integer, dimension(:,:), intent(in), optional    :: kbot
real, dimension(:,:,:),  intent(in), optional    :: cloud_water_input, &
                                                    cloud_ice_input, &
                                                    cloudtemp,    &
                                                    cloudvapor

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      pfull        pressure at full levels [ kg / (m s^2) ]
!      phalf        pressure at half levels [ kg / (m s^2) ]
!      t            temperature at full levels [ deg K]
!      q            specific humidity of water vapor at full levels
!                   [ dimensionless ]
!      ts           surface temperature  [ deg K ]
!      albedo       surface albedo  [ dimensionless ]
!      land         fraction of grid box which is land [ dimensionless ]
!
!   intent(out) variables:
!
!      Atmos_input   atmos_input type structure, contains the 
!                    following components defined in this subroutine
!         psfc          surface pressure 
!                       [ (kg /( m s^2) ] 
!         tsfc          surface temperature
!                       [ deg K ]
!         asfc          surface albedo
!                       [ non-dimensional ]
!         land          fraction of grid box covered by land
!                       [ non-dimensional ]
!         temp          temperature at model levels (1:nlev), surface
!                       temperature is stored at value nlev+1; if eta
!                       coordinates, surface value stored in below 
!                       ground points
!                       [ deg K ]
!         press         pressure at model levels (1:nlev), surface 
!                       pressure is stored at index value nlev+1
!                       [ (kg /( m s^2) ] 
!         rh2o          mixing ratio of water vapor at model full levels
!                       [ non-dimensional ]
!         cloud_water   cloud water mixing ratio (or specific humidity 
!                       ????)
!                       [ non-dimensional ]
!         cloud_ice     cloud ice mixing ratio (or specific humidity 
!                       ????)
!                       [ non-dimensional ]
!         deltaz        model vertical grid separation
!                       [meters]
!         pflux         average of pressure at adjacent model levels
!                       [ (kg /( m s^2) ] 
!         tflux         average of temperature at adjacent model levels
!                       [ deg K ]
!         rel_hum       relative humidity
!                       [ dimensionless ]
!         cloudtemp     temperature to be seen by clouds (used in 
!                       sa_gcm feedback studies) 
!                       [ degrees K ]
!         cloudvapor    water vapor to be seen by clouds (used in 
!                       sa_gcm feedback studies) 
!                       [ nondimensional ]
!         clouddeltaz   deltaz to be used in defining cloud paths (used
!                       in sa_gcm feedback studies)
!                       [ meters ]
!
!   intent(in), optional variables:
!
!      kbot               present when running eta vertical coordinate,
!                         index of lowest model level above ground (???)
!      cloud_water_input  cloud water mixing ratio (or specific humidity
!                         ????)
!                         [ non-dimensional ]
!      cloud_ice_input    cloud ice mixing ratio (or specific humidity 
!                         ????)
!                         [ non-dimensional ]
!      cloudtemp          temperature to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!                         [ degrees K ]
!      cloudvapor         water vapor to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!                         [ nondimensional ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables
 
      integer :: i, j, k, kb

!---------------------------------------------------------------------
!  local variables
!
!     i, j, k      do loop indices
!     kb           vertical index of lowest atmospheric level (when
!                  using eta coordinates)
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    allocate space for the components of the derived type variable
!    Atmos_input.
!---------------------------------------------------------------------
      allocate ( Atmos_input%press(size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%phalf(size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%temp (size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%rh2o (size(t,1), size(t,2), size(t,3)  ) )
      allocate ( Atmos_input%rel_hum(size(t,1), size(t,2),    &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%cloud_ice(size(t,1), size(t,2),   &
                                                        size(t,3)  ) )
      allocate ( Atmos_input%cloud_water(size(t,1), size(t,2),   &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%cloudtemp(size(t,1), size(t,2),   &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%cloudvapor(size(t,1), size(t,2),   &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%clouddeltaz(size(t,1), size(t,2),   &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%deltaz(size(t,1), size(t,2), size(t,3) ) )
      allocate ( Atmos_input%pflux (size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%tflux (size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%asfc (size(t,1), size(t,2)             ) )
      allocate ( Atmos_input%psfc (size(t,1), size(t,2)             ) )
      allocate ( Atmos_input%tsfc (size(t,1), size(t,2)             ) )
      allocate ( Atmos_input%land (size(t,1), size(t,2)             ) )

!---------------------------------------------------------------------
!    be sure cloud_water and cloud_ice have been supplied if they are
!    required by the namelist settings.
!---------------------------------------------------------------------
      if (Environment%running_standalone .and. &
          .not. Environment%running_sa_model) then
        if (Cldrad_control%do_pred_cld_microphys) then
          if (.not. present(cloud_water_input) .or.  &
              .not. present(cloud_ice_input) ) then
            call error_mesg ('radiation_driver_mod', &
            ' cloud_ice and cloud_water values must be supplied when&
             & using pred_cld_microphys', FATAL)
          endif
        endif
      endif

!---------------------------------------------------------------------
!    define the cloud_water and cloud_ice components of Atmos_input.
!---------------------------------------------------------------------
      if (present (cloud_ice_input) .and. &
          present (cloud_water_input) ) then
        Atmos_input%cloud_ice(:,:,:)   = cloud_ice_input(:,:,:)
        Atmos_input%cloud_water(:,:,:) = cloud_water_input(:,:,:)
      else
        Atmos_input%cloud_ice(:,:,:)   = 0.0
        Atmos_input%cloud_water(:,:,:) = 0.0
      endif

!---------------------------------------------------------------------
!    define the cloudtemp component of Atmos_input.
!---------------------------------------------------------------------
      if (present (cloudtemp) ) then
        Atmos_input%cloudtemp(:,:,:)   = cloudtemp(:,:,:)
      else
        Atmos_input%cloudtemp(:,:,:)   = t(:,:,kmin:kmax)
      endif

!---------------------------------------------------------------------
!    define the cloudvapor component of Atmos_input.
!---------------------------------------------------------------------
      if (present (cloudvapor) ) then
        Atmos_input%cloudvapor(:,:,:)   = cloudvapor(:,:,:)
      else
        Atmos_input%cloudvapor(:,:,:)   = q(:,:,kmin:kmax)
      endif

!---------------------------------------------------------------------
!    define values of surface pressure, temperature and albedo.
!--------------------------------------------------------------------
      if (present(kbot)) then
        do j=1,je-js+1
          do i=1,ie-is+1
            kb = kbot(i,j)
            Atmos_input%psfc(i,j) = phalf(i,j,kb+1)
          end do
        end do
      else
        Atmos_input%psfc(:,:) = phalf(:,:,kmax+1)
      endif

      Atmos_input%tsfc(:,:) = ts(:,:)
      Atmos_input%asfc(:,:) = albedo(:,:)

!------------------------------------------------------------------
!    define the atmospheric pressure and temperature arrays.
!------------------------------------------------------------------
      do k=kmin,kmax 
        Atmos_input%press(:,:,k) = pfull(:,:,k)
        Atmos_input%phalf(:,:,k) = phalf(:,:,k)
        Atmos_input%temp (:,:,k) = t(:,:,k)
      end do
      Atmos_input%press(:,:,kmax+1) = phalf(:,:,kmax+1)
      Atmos_input%phalf(:,:,kmax+1) = phalf(:,:,kmax+1)
      Atmos_input%temp (:,:,kmax+1) = ts  (:,:)

!------------------------------------------------------------------
!    if in eta coordinates, fill in underground temperatures with 
!    surface value.
!------------------------------------------------------------------
      if (present(kbot)) then
        do j=1,je-js+1
          do i=1,ie-is+1
            kb = kbot(i,j)
            if (kb < kmax) then
              do k=kb+1,kmax
                Atmos_input%temp(i,j,k) = Atmos_input%temp(i,j,kmax+1)
              end do
            endif
          end do
        end do
      endif

!------------------------------------------------------------------
!    when running the gcm, convert the input water vapor specific 
!    humidity field to mixing ratio. it is assumed that water vapor 
!    mixing ratio is the input in the standalone case.
!------------------------------------------------------------------
      if (Environment%running_gcm .or. &
          Environment%running_sa_model) then
        if(use_mixing_ratio) then
          Atmos_input%rh2o (:,:,:) = q(:,:,:)
        else
          Atmos_input%rh2o (:,:,:) = q(:,:,:)/(1.0 - q(:,:,:))
          Atmos_input%cloudvapor(:,:,:) =    &
                                 Atmos_input%cloudvapor(:,:,:)/  &
                            (1.0 - Atmos_input%cloudvapor(:,:,:))
        endif
      else if (Environment%running_standalone) then
        Atmos_input%rh2o (:,:,:) = q(:,:,:)
      endif ! (running gcm)
 
!------------------------------------------------------------------
!    define the fractional land area of each grid box.
!------------------------------------------------------------------
      Atmos_input%land(:,:) = land(:,:)
     
!------------------------------------------------------------------
!    be sure that the magnitude of the water vapor mixing ratio field 
!    to be input to the radiation code is no smaller than the value of 
!    rh2o_lower_limit, which is 2.0E-07 when running the sea_esf
!    radiation code and 3.0e-06 when running the original radiation
!    code. Likewise, the temperature that the radiation code sees is
!    constrained to lie between 100K and 370K. these are the limits of
!    the tables referenced within the radiation package.
!-----------------------------------------------------------------------
      if (do_rad .or. do_average) then
        Atmos_input%rh2o(:,:,ks:ke) =    &
	            MAX(Atmos_input%rh2o(:,:,ks:ke), rh2o_lower_limit)
        Atmos_input%cloudvapor(:,:,ks:ke) =    &
	            MAX(Atmos_input%cloudvapor(:,:,ks:ke), rh2o_lower_limit)
        Atmos_input%temp(:,:,ks:ke) =     &
                    MAX(Atmos_input%temp(:,:,ks:ke), temp_lower_limit)
        Atmos_input%temp(:,:,ks:ke) =     &
                    MIN(Atmos_input%temp(:,:,ks:ke), temp_upper_limit)
        Atmos_input%cloudtemp(:,:,ks:ke) =     &
                    MAX(Atmos_input%cloudtemp(:,:,ks:ke), temp_lower_limit)
        Atmos_input%cloudtemp(:,:,ks:ke) =     &
                    MIN(Atmos_input%cloudtemp(:,:,ks:ke), temp_upper_limit)
      endif

!--------------------------------------------------------------------
!    call calculate_aulixiary_variables to compute pressure and 
!    temperature arrays at flux levels and an array of model deltaz.
!--------------------------------------------------------------------
      if (do_rad) then
        call calculate_auxiliary_variables (Atmos_input)
      endif

!----------------------------------------------------------------------


end subroutine define_atmos_input_fields 


!####################################################################


subroutine calculate_auxiliary_variables (Atmos_input)

!----------------------------------------------------------------------
!    calculate_auxiliary_variables defines values of model delta z and
!    relative humidity, and the values of pressure and temperature at
!    the grid box vertical interfaces.
!---------------------------------------------------------------------

type(atmos_input_type), intent(inout)  :: Atmos_input

!--------------------------------------------------------------------
!   intent(inout) variables
!
!      Atmos_input  atmos_input_type variable, its press and temp
!                   components are input, and its deltaz, rel_hum, 
!                   pflux and tflux components are calculated here 
!                   and output.
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables

      real, dimension (size(Atmos_input%temp, 1), &
                       size(Atmos_input%temp, 2), &
                       size(Atmos_input%temp, 3) - 1) :: &
                                                     esat, psat, qv, tv
      integer   ::  k


!--------------------------------------------------------------------
!    define flux level pressures (pflux) as midway between data level
!    (layer-mean) pressures. specify temperatures at flux levels
!    (tflux).
!--------------------------------------------------------------------
      do k=ks+1,ke
        Atmos_input%pflux(:,:,k) = 0.5E+00*  &
                (Atmos_input%press(:,:,k-1) + Atmos_input%press(:,:,k))
        Atmos_input%tflux(:,:,k) = 0.5E+00*  &
                (Atmos_input%temp (:,:,k-1) + Atmos_input%temp (:,:,k))
      end do
      Atmos_input%pflux(:,:,ks  ) = 0.0E+00
      Atmos_input%pflux(:,:,ke+1) = Atmos_input%press(:,:,ke+1)
      Atmos_input%tflux(:,:,ks  ) = Atmos_input%temp (:,:,ks  )
      Atmos_input%tflux(:,:,ke+1) = Atmos_input%temp (:,:,ke+1)

!-------------------------------------------------------------------
!    define deltaz in meters
!-------------------------------------------------------------------
      tv(:,:,:) = Atmos_input%temp(:,:,ks:ke)*    &
                  (1.0 + d608*Atmos_input%rh2o(:,:,:))
      Atmos_input%deltaz(:,:,ks) = log_p_at_top*RDGAS*tv(:,:,ks)/GRAV
      do k =ks+1,ke   
        Atmos_input%deltaz(:,:,k) = alog(Atmos_input%pflux(:,:,k+1)/  &
                                         Atmos_input%pflux(:,:,k))*   &
                                         RDGAS*tv(:,:,k)/GRAV
      end do

!-------------------------------------------------------------------
!    define deltaz in meters to be used in cloud feedback analysis
!-------------------------------------------------------------------
      tv(:,:,:) = Atmos_input%cloudtemp(:,:,ks:ke)*    &
                  (1.0 + d608*Atmos_input%cloudvapor(:,:,:))
      Atmos_input%clouddeltaz(:,:,ks) = log_p_at_top*RDGAS*  &
                                        tv(:,:,ks)/GRAV
      do k =ks+1,ke   
        Atmos_input%clouddeltaz(:,:,k) =    &
	                            alog(Atmos_input%pflux(:,:,k+1)/  &
                                         Atmos_input%pflux(:,:,k))*   &
                                         RDGAS*tv(:,:,k)/GRAV
      end do

!------------------------------------------------------------------
!    define the relative humidity.
!------------------------------------------------------------------
      call lookup_es (Atmos_input%temp(:,:,1:kmax), esat)
      do k=1,kmax
        qv(:,:,k) = Atmos_input%rh2o(:,:,k) /    &
                                       (1.0 + Atmos_input%rh2o(:,:,k))
        psat(:,:,k) = Atmos_input%press(:,:,k) - d378*esat(:,:,k)
        psat(:,:,k) = MAX (psat(:,:,k), esat(:,:,k))
        Atmos_input%rel_hum(:,:,k) = qv(:,:,k) / (d622*esat(:,:,k) / &
                                                  psat(:,:,k))
        Atmos_input%rel_hum(:,:,k) =    &
                                  MIN (Atmos_input%rel_hum(:,:,k), 1.0)
      end do

!----------------------------------------------------------------------


end subroutine calculate_auxiliary_variables



!#####################################################################

subroutine time_average_input_data (is, ie, js, je, Atmos_input,  &
                                    Astro, Rad_gases, Cldrad_props) 

!---------------------------------------------------------------------
!    when time-averaged input data is desired, time_average_input_data 
!    calls compute_average to accumulate the sum of the fields on each 
!    timestep until a radiation step is reached. on that step, 
!    return_average is called to produce a time-averaged value which 
!    will be returned for use in the radiation calculation.
!    the components of Atmos_input and Astro are always averaged when
!    averaging is activated; previously values of Rad_gases and 
!    Cldrad_props were not, but now can be, if desired, by including
!    them as optional arguments in the call to this subroutine, and
!    supplying the needed code as indicated at various places in this
!    module.
!    
!---------------------------------------------------------------------
     
integer,                      intent(in)              :: is, ie, js, je
type(atmos_input_type),       intent(inout)           :: Atmos_input
type(astronomy_type),         intent(inout)           :: Astro
type(cldrad_properties_type), intent(inout), optional :: Cldrad_props
type(radiative_gases_type),   intent(inout), optional :: Rad_gases

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(inout) variables:
!   
!      Atmos_input  atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!      Astro        astronomy_type variable containing the astronomical
!                   input fields needed by the radiation package
!
!   intent(inout), optional variables:
!
!      Cldrad_props cldrad_properties_type variable containing the 
!                   cloud radiative property input fields needed by the 
!                   radiation package
!      Rad_gases    radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!   call compute_average to add the fields input on this step to the 
!   accumulating sum. 
!-----------------------------------------------------------------------
      if (present(Rad_gases)) then
        if (present(Cldrad_props)) then
          call compute_average  (is, ie, js, je, Atmos_input, Astro, &
                                 Rad_gases=Rad_gases, &
                                 Cldrad_props=Cldrad_props) 
        else
          call compute_average  (is, ie, js, je, Atmos_input, Astro, &
                                 Rad_gases=Rad_gases) 
        endif
      else
        if (present(Cldrad_props)) then
          call compute_average  (is, ie, js, je, Atmos_input, Astro, &
                                 Cldrad_props=Cldrad_props) 
        else
          call compute_average  (is, ie, js, je, Atmos_input, Astro)
        endif
      endif

!-------------------------------------------------------------------
!    if this is a radiation timestep, call return_average to obtain
!    the time-averaged input fields.
!-------------------------------------------------------------------
      if (do_rad) then
        if (present(Rad_gases)) then
          if (present(Cldrad_props)) then
            call return_average  (is, ie, js, je, Atmos_input, Astro, &
                                  Rad_gases=Rad_gases,   &
                                  Cldrad_props=Cldrad_props) 
          else
            call return_average  (is, ie, js, je, Atmos_input, Astro, &
                                  Rad_gases=Rad_gases) 
          endif
        else
          if (present(Cldrad_props)) then
            call return_average  (is, ie, js, je, Atmos_input, Astro, &
                                  Cldrad_props=Cldrad_props) 
          else
            call return_average  (is, ie, js, je, Atmos_input, Astro)
          endif
        endif

!--------------------------------------------------------------------
!    call calculate_aulixiary_variables to compute time-averaged deltaz,
!    relative humidity, and pressure and temperature arrays at flux 
!    levels, using the time_averaged values of press, temp and water
!    vapor mixing ratio just obtained.
!--------------------------------------------------------------------
        call calculate_auxiliary_variables (Atmos_input)
      endif

!----------------------------------------------------------------------


end subroutine time_average_input_data 


!#####################################################################

subroutine compute_average (is, ie, js, je, Atmos_input, Astro,  &
                            Rad_gases, Cldrad_props)

!---------------------------------------------------------------------
!    compute_average adds the current value of the input fields to the
!    accumulating sum which will be averaged on the next radiation 
!    time step.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
integer,                      intent(in)           :: is, ie, js, je
type(atmos_input_type),       intent(in)           :: Atmos_input
type(astronomy_type),         intent(in)           :: Astro
type(radiative_gases_type),   intent(in), optional :: Rad_gases
type(cldrad_properties_type), intent(in), optional :: Cldrad_props
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Atmos_input  atmos_input_type variable containing atmos-
!                   pheric input data for the radiation package 
!                   on the model grid
!      Astro        astronomy_type variable containing astronomical
!                   input data for the radiation package on the 
!                   model grid
!
!   intent(in), optional variables:
!
!      Rad_gases    radiative_gases_type variable containing rad-
!                   iative gas input data for the radiation package
!                   on the model grid
!      Cldrad_props cldrad_properties_type variable containing 
!                   cloud radiative property input data for the 
!                   radiation package on the model grid
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    add the current values of the Atmos_input fields to the accumul-
!    ation sum. increment the counter of the number of values making 
!    up the sum.
!---------------------------------------------------------------------
      nsum(is:ie,js:je)   = nsum(is:ie,js:je)   + 1
      psum(is:ie,js:je,:) = psum(is:ie,js:je,:) +   &
                                           Atmos_input%press(:,:,:)
      tsum(is:ie,js:je,:) = tsum(is:ie,js:je,:) +   &
                                           Atmos_input%temp (:,:,:)
      qsum(is:ie,js:je,:) = qsum(is:ie,js:je,:) +   &
                                           Atmos_input%rh2o (:,:,:)
      asum(is:ie,js:je)   = asum(is:ie,js:je)   +   &
                                           Atmos_input%asfc(:,:)

!---------------------------------------------------------------------
!    add the current values of the Astro fields to the accumulation sum.
!---------------------------------------------------------------------
      csum(is:ie,js:je)   = csum (is:ie,js:je)   + Astro%cosz  (:,:)
      ssum(is:ie,js:je)   = ssum (is:ie,js:je)   + Astro%solar (:,:)
      fsum (is:ie,js:je)  = fsum (is:ie,js:je)   + Astro%fracday(:,:)
      rrsum(is:ie,js:je)  = rrsum(is:ie, js:je)  + Astro%rrsun   

!---------------------------------------------------------------------
!    add the current values of the Rad_gases fields to the accumulation 
!    sum. Code not yet provided.
!---------------------------------------------------------------------
      if (present(Rad_gases) ) then
        call error_mesg ('radiation_driver', &
         ' do_average_gases not yet implemented', FATAL)
      endif

!---------------------------------------------------------------------
!    add the current values of the Cldrad_props fields to the accum-
!    ulation sum. Code not yet provided.
!---------------------------------------------------------------------
      if (present(Cldrad_props)) then
        call error_mesg ('radiation_driver', &
         ' do_average_clouds not yet implemented', FATAL)
      endif

!-----------------------------------------------------------------------


end subroutine compute_average



!#######################################################################

subroutine return_average (is, ie, js, je, Atmos_input, Astro,    &
                           Rad_gases, Cldrad_props)

!----------------------------------------------------------------------
!    return_average computes the time-averaged input fields on radiation
!    steps and returns those values. it also reinintializes the
!    accumulation sums.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
integer,                      intent(in)            :: is, ie, js, je
type(atmos_input_type),       intent(inout)           :: Atmos_input
type(astronomy_type),         intent(inout)           :: Astro
type(radiative_gases_type),   intent(inout), optional :: Rad_gases  
type(cldrad_properties_type), intent(inout), optional :: Cldrad_props 
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!   intent(out) variables:
!
!      Atmos_input  atmos_input_type variable containing atmos-
!                   pheric input data for the radiation package 
!                   on the model grid
!      Astro        astronomy_type variable containing astronomical
!                   input data for the radiation package on the 
!                   model grid
!
!   intent(out), optional variables:
!
!      Rad_gases    radiative_gases_type variable containing rad-
!                   iative gas input data for the radiation package
!                   on the model grid
!      Cldrad_props cldrad_properties_type variable containing 
!                   cloud radiative property input data for the 
!                   radiation package on the model grid
!
!--------------------------------------------------------------------

!-----------------------------------------------------------------------
!   local variables
      real, dimension (size(Atmos_input%press,1),   &
                       size(Atmos_input%press,2)) :: dfsum
      integer   :: n, k

!---------------------------------------------------------------------
!    check for zero divide.
!---------------------------------------------------------------------
      n = count (nsum(is:ie,js:je) <= 0)
      if ( n > 0 ) then
          call error_mesg ('radiation_driver_mod',  &
           'no data present to compute average in return_average.', &
            FATAL)
      endif

!---------------------------------------------------------------------
!    compute averages of Atmos_input components.
!---------------------------------------------------------------------
      dfsum(:,:) = 1.0 / float(nsum(is:ie,js:je))
      do k=1,size(Atmos_input%press,3)
         Atmos_input%press(:,:,k) = psum(is:ie,js:je,k) * dfsum(:,:)
         ATmos_input%temp (:,:,k) = tsum(is:ie,js:je,k) * dfsum(:,:)
      enddo

      do k=1,size(Atmos_input%rh2o,3)
         Atmos_input%rh2o (:,:,k) = qsum(is:ie,js:je,k) * dfsum(:,:)
      enddo

      Atmos_input%tsfc(:,:) = Atmos_input%temp (:,:,kmax+1)
      Atmos_input%psfc(:,:) = Atmos_input%press(:,:,kmax+1)
      Atmos_input%asfc(:,:) = asum(is:ie,js:je) *dfsum(:,:)

!---------------------------------------------------------------------
!    compute averages of Astro components.
!---------------------------------------------------------------------
      Astro%cosz  (:,:) = csum(is:ie,js:je) *dfsum(:,:)
      Astro%fracday (:,:) = fsum(is:ie,js:je) *dfsum(:,:)
      where (Astro%fracday > 1.00) Astro%fracday =  1.00
      Astro%rrsun         = rrsum(is,js) *dfsum(1,1)
      Astro%solar (:,:) = ssum(is:ie,js:je) *dfsum(:,:)

!---------------------------------------------------------------------
!    include code here to average Rad_gases components which are to 
!    be time-averaged. 
!---------------------------------------------------------------------
      if (present (Rad_gases)) then
        call error_mesg ('radiation_driver_mod', &
         ' do_average_gases not yet implemented', FATAL)
      endif

!---------------------------------------------------------------------
!    include code here to average Cldrad_props components which are to 
!    be time-averaged.
!---------------------------------------------------------------------
      if (present (Cldrad_props)) then
        call error_mesg ('radiation_driver_mod', &
         ' do_average_clouds not yet implemented', FATAL)
      endif

!---------------------------------------------------------------------
!    reinitialize sums so that new accumulations may be made until
!    the next radiation step.
!---------------------------------------------------------------------
      nsum(is:ie,js:je)   = 0
      psum(is:ie,js:je,:) = 0.0
      tsum(is:ie,js:je,:) = 0.0
      qsum(is:ie,js:je,:) = 0.0
      asum(is:ie,js:je)   = 0.0
      csum(is:ie,js:je)   = 0.0
      ssum(is:ie,js:je)   = 0.0
      fsum(is:ie,js:je)   = 0.0
      rrsum(is:ie,js:je)  = 0.0

!---------------------------------------------------------------------
!    include code here to reinitialize sums of  Rad_gases components 
!    which are to be time-averaged.
!---------------------------------------------------------------------
      if (present (Rad_gases)) then
        call error_mesg ('radiation_driver_mod', &
         ' do_average_gases not yet implemented', FATAL)
      endif

!---------------------------------------------------------------------
!    include code here to reinitialize sums of  Cldrad_props components 
!    which are to be time-averaged.
!---------------------------------------------------------------------
      if (present (Cldrad_props)) then
        call error_mesg ('radiation_driver_mod', &
         ' do_average_clouds not yet implemented', FATAL)
      endif

!-----------------------------------------------------------------------



end subroutine return_average



!#######################################################################

subroutine radiation_calc (is, ie, js, je, phalf, Rad_time, Time_diag,  &
                           lat_mdl, lon_mdl, Atm_inp_mdl,  &
                           Rad_gases_mdl, Aerosol_mdl,  &
                           Cldrad_props_mdl,  &
                           Cld_diagnostics, Astro_mdl, &
                           Rad_output, Lw_output, Sw_output,   &
                           Fsrad_output, mask, kbot)

!--------------------------------------------------------------------
!    radiation_calc is called on radiation timesteps and calculates
!    the long- and short-wave radiative fluxes and heating rates, and
!    obtains any radiation output fields needed in other portions of
!    the model.
!-----------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                      intent(in)             :: is, ie, js, je
real,                         intent(in)             :: phalf(:,:,:)
type(time_type),              intent(in)             :: Rad_time,   &
                                                        Time_diag
real, dimension(:,:),         intent(in)             :: lat_mdl, lon_mdl
type(atmos_input_type),       intent(in)             :: Atm_inp_mdl
type(radiative_gases_type),   intent(in)             :: Rad_gases_mdl
type(aerosol_type),           intent(in)             :: Aerosol_mdl
type(cldrad_properties_type), intent(in)             :: Cldrad_props_mdl
type(cld_diagnostics_type),   intent(in)             :: Cld_diagnostics
type(astronomy_type),         intent(in)             :: Astro_mdl
type(rad_output_type),        intent(inout)            :: Rad_output
type(lw_output_type),         intent(inout)            :: Lw_output
type(sw_output_type),         intent(inout)            :: Sw_output
type(fsrad_output_type),      intent(inout)            :: Fsrad_output
real, dimension(:,:,:),       intent(in),   optional :: mask
integer, dimension(:,:),      intent(in),   optional :: kbot

!-----------------------------------------------------------------------
!    intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      Rad_time          time at which the radiative fluxes are to apply
!                        [ time_type (days, seconds) ] 
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ]  
!      lat_mdl           latitude of model points on model grid 
!                        [ radians ]
!      lon_mdl           longitude of model points on model grid 
!                        [ radians ]
!      Atm_inp_mdl       atmos_input_type variable containing atmos-
!                        pheric input data for the radiation package 
!                        on the model grid
!      Rad_gases_mdl     radiative_gases_type variable containing rad-
!                        iative gas input data for the radiation package
!                        on the model grid
!      Cldrad_props_mdl  cldrad_properties_type variable containing 
!                        cloud radiative property input data for the 
!                        radiation package on the model grid
!      Cld_diagnostics   cld_diagnostics_type variable containing cloud
!                        microphysical diagnostic quantities
!      Astro_mdl         astronomy_type variable containing astronomical
!                        input data for the radiation package on the 
!                        model grid
!
!    intent(out) variables:
!
!      Lw_output         lw_output_type variable containing longwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid, when that 
!                        package is active
!      Sw_output         sw_output_type variable containing shortwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid, when that 
!                        package is active
!      Fsrad_output      fsrad_output_type variable containing radiation
!                        output data from the original_fms_rad radiation
!                        package on the model grid, when that package 
!                        is active
!      Rad_output        rad_output_type variable containing radiation
!                        output data needed by other modules
!
!    intent(in), optional variables:
!
!      mask              present when running eta vertical coordinate,
!                        mask to define values at points below ground   
!      kbot              present when running eta vertical coordinate,
!                        index of lowest model level above ground 
!
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
! local variables

      type(atmos_input_type)       :: Atmos_input
      type(astronomy_type)         :: Astro
      type(radiative_gases_type)   :: Rad_gases
      type(aerosol_type)           :: Aerosol
      type(cldrad_properties_type) :: Cldrad_props
      type(lw_output_type)         :: Lw_output_rd
      type(sw_output_type)         :: Sw_output_rd
      type(fsrad_output_type)      :: Fsrad_output_rd
      integer                      :: ierad, jerad

!---------------------------------------------------------------------
!    define the (i,j,k) points over which radiation is to be calculated.
!    this code is included for possible future use, allowing radiation
!    to be calculated on a subset of model columns and / or levels. in 
!    this code release, all_column_radiation and all_level_radiation 
!    must both be .true..  note that kerad was defined in the _init 
!    subroutine.
!---------------------------------------------------------------------
      if (all_column_radiation .and. all_level_radiation) then

!--------------------------------------------------------------------
!    call routines to perform radiation calculations, either using the
!    sea_esf_rad or original_fms_rad radiation package.
!---------------------------------------------------------------------
        if (do_sea_esf_rad) then
          call sea_esf_rad (is, ie, js, je, Atm_inp_mdl, Astro_mdl,&
                            Rad_gases_mdl,  Aerosol_mdl,  &
                            Cldrad_props_mdl, Cld_diagnostics, &
                            Lw_output, Sw_output)
        else  
          call original_fms_rad (is, ie, js, je, phalf, lat_mdl,   &
                                 lon_mdl, do_clear_sky_pass, Rad_time, &
                                 Time_diag, Atm_inp_mdl, Astro_mdl, &
                                 Rad_gases_mdl, Cldrad_props_mdl,  &
                                 Fsrad_output, mask=mask, kbot=kbot) 
        endif
  
!---------------------------------------------------------------------
!    if radiation were not to be called in all columns and at all 
!    levels, determine the number of (i,j) columns in which radiation 
!    is to be calculated (ierad, jerad). 
!    add code here to define values of ierad and jerad
!---------------------------------------------------------------------
      else  
!       ierad =   ????
!       jerad =  ????

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    when this option is coded, replace this error_mesg code with
!    the code which follows it.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call error_mesg ('radiation_driver_mod', &
               ' ability to calculate radiation on subset of columns&
              & and/or levels not yet implemented',  FATAL)

!---------------------------------------------------------------------
!    allocate the variables needed to store the input fields on the rad-
!    iation grid.
!---------------------------------------------------------------------
!!! ARRAYS WILL NEED TO BE NULLIFIED OR INITIALIZED!
        allocate ( Atmos_input%press      (ierad, jerad, kerad+1) )
        allocate ( Atmos_input%temp       (ierad, jerad, kerad+1) )
        allocate ( Atmos_input%rh2o       (ierad, jerad, kerad  ) )
        allocate ( Atmos_input%rel_hum    (ierad, jerad, kerad  ) )
        allocate ( Atmos_input%cloud_ice  (ierad, jerad, kerad  ) )
        allocate ( Atmos_input%cloud_water(ierad, jerad, kerad  ) )    
        allocate ( Atmos_input%cloudtemp  (ierad, jerad, kerad  ) )    
        allocate ( Atmos_input%cloudvapor (ierad, jerad, kerad  ) )    
        allocate ( Atmos_input%clouddeltaz(ierad, jerad, kerad  ) )    
        allocate ( Atmos_input%deltaz     (ierad, jerad, kerad  ) )
        allocate ( Atmos_input%pflux      (ierad, jerad, kerad+1) )
        allocate ( Atmos_input%tflux      (ierad, jerad, kerad+1) )
        allocate ( Atmos_input%asfc       (ierad, jerad         ) )
        allocate ( Atmos_input%psfc       (ierad, jerad         ) )
        allocate ( Atmos_input%tsfc       (ierad, jerad         ) )
        allocate ( Atmos_input%land       (ierad, jerad         ) )

        allocate ( Astro%fracday(ierad, jerad ) )
        allocate ( Astro%cosz   (ierad, jerad ) )
        allocate ( Astro%solar  (ierad, jerad ) )
     
        allocate ( Rad_gases%qo3(ierad, jerad, kerad) )

        call cloudrad_package_alloc_rad (ierad, jerad, kerad,   &
                                         Cldrad_props)

!--------------------------------------------------------------------
!    call model_to_radiation_grid to define input fields on the 
!    radiation grid. 
!--------------------------------------------------------------------
!!  ???? IS CLD_DIAGNOSTICS NEEDED IN THIS CALL ?????????????
        call model_to_radiation_grid (is, ie, js, je, Atm_inp_mdl,   &
                                      Astro_mdl, Rad_gases_mdl,   &
                                      Cldrad_props_mdl,  &
				      Cld_diagnostics, Atmos_input, &
				      Astro, Rad_gases, Cldrad_props)  

!---------------------------------------------------------------------
!    allocate the variables needed to store the output fields on the 
!    radiation grid.
!---------------------------------------------------------------------
!!!!!!!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!! NIOTE :: Cld_diagnostics not fixed HERE !!!!!!!!!!!!!!!
        if (do_sea_esf_rad) then
!!! ARRAYS WILL NEED TO BE NULLIFIED OR INITIALIZED!
          allocate (Lw_output_rd%heatra (ierad, jerad, kerad) )
          allocate (Lw_output_rd%flxnet (ierad, jerad, kerad+1) )
          allocate (Sw_output_rd%dfsw   (ierad, jerad, kerad+1) )
          allocate (Sw_output_rd%ufsw   (ierad, jerad, kerad+1) )
          allocate (Sw_output_rd%fsw    (ierad, jerad, kerad+1) )
          allocate (Sw_output_rd%hsw    (ierad, jerad, kerad  ) )
          if (do_clear_sky_pass) then
            allocate (Lw_output_rd%heatracf (ierad, jerad, kerad) )
            allocate (Lw_output_rd%flxnetcf (ierad, jerad, kerad+1) )
            allocate (Sw_output_rd%dfswcf   (ierad, jerad, kerad+1) )
            allocate (Sw_output_rd%ufswcf   (ierad, jerad, kerad+1) )
            allocate (Sw_output_rd%fswcf    (ierad, jerad, kerad+1) )
            allocate (Sw_output_rd%hswcf    (ierad, jerad, kerad  ) )
          endif
        else
          allocate (Fsrad_output_rd%tdtsw    ( ierad, jerad, kerad) )
          allocate (Fsrad_output_rd%tdtlw    ( ierad, jerad, kerad) )
          allocate (Fsrad_output_rd%swdns    ( ierad, jerad       ) )
          allocate (Fsrad_output_rd%swups    ( ierad, jerad       ) )
          allocate (Fsrad_output_rd%lwdns    ( ierad, jerad       ) )
          allocate (Fsrad_output_rd%lwups    ( ierad, jerad       ) )
          allocate (Fsrad_output_rd%swin     ( ierad, jerad       ) )
          allocate (Fsrad_output_rd%swout    ( ierad, jerad       ) )
          allocate (Fsrad_output_rd%olr      ( ierad, jerad       ) )
          if (do_clear_sky_pass) then
            allocate (Fsrad_output_rd%tdtsw_clr( ierad, jerad, kerad) )
            allocate (Fsrad_output_rd%tdtlw_clr( ierad, jerad, kerad) )
            allocate (Fsrad_output_rd%swdns_clr( ierad, jerad       ) )
            allocate (Fsrad_output_rd%swups_clr( ierad, jerad       ) )
            allocate (Fsrad_output_rd%lwdns_clr( ierad, jerad       ) )
            allocate (Fsrad_output_rd%lwups_clr( ierad, jerad       ) )
            allocate (Fsrad_output_rd%swin_clr ( ierad, jerad       ) )
            allocate (Fsrad_output_rd%swout_clr( ierad, jerad       ) )
            allocate (Fsrad_output_rd%olr_clr  ( ierad, jerad       ) )
          endif
        endif   ! (do_sea_esf_rad)

!--------------------------------------------------------------------
!    now call routines to perform radiation calculations, using the 
!    just-defined variables on the radiation grid. Note that output 
!    is returned in arrays on the radiation grid.
!---------------------------------------------------------------------
        if (do_sea_esf_rad) then
          call sea_esf_rad (is, ie, js, je, Atmos_input, Astro,  &
                            Rad_gases, Aerosol, Cldrad_props,   &
                            Cld_diagnostics, &
                            Lw_output_rd, Sw_output_rd) 
        else
          call original_fms_rad (is, ie, js, je, phalf, lat_mdl, lon_mdl,   &
                                 do_clear_sky_pass, Rad_time, &
                                 Time_diag, Atmos_input, Astro,  &
                                 Rad_gases, Cldrad_props,  &
                                 Fsrad_output_rd, mask=mask, kbot=kbot) 
        endif

!-------------------------------------------------------------------
!    call radiation_to_model_grid to convert the output from sea_esf_rad
!    or original_fms_rad to the model grid.  these fields are then 
!    passed back to radiation_driver.
!-------------------------------------------------------------------
        if (do_sea_esf_rad) then
          call radiation_to_model_grid      &
                                    (is, ie, js, je,        &
                                     Lw_output_rd=Lw_output_rd,   &
                                     Sw_output_rd=Sw_output_rd,   &
                                     Lw_output=Lw_output, &
                                     Sw_output=Sw_output)
        else
          call radiation_to_model_grid       &
                                    (is, ie, js, je,        &
                                     Fsrad_output_rd=Fsrad_output_rd, &
                                     Fsrad_output=Fsrad_output)
        endif
      endif  ! (all_column .and. all_level)

!---------------------------------------------------------------------
!    define the components of Rad_output to be passed back to 
!    radiation_driver --  total and shortwave radiative heating rates 
!    for standard and clear-sky case (if desired), and surface long- 
!    and short-wave fluxes.  mask out any below ground values if 
!    necessary.
!---------------------------------------------------------------------
      if (Environment%running_gcm .or.  &
          Environment%running_sa_model) then
        if (do_sea_esf_rad) then
          Rad_output%tdtsw(is:ie,js:je,:) = Sw_output%hsw(:,:,:)/  &
                                            seconds_per_day
          if (present(mask)) then
            Rad_output%tdtlw(is:ie,js:je,:) =   &
                           (Lw_output%heatra(:,:,:)/seconds_per_day)*  &
                            mask(:,:,:)

            Rad_output%tdt_rad (is:ie,js:je,:) =   &
                           (Rad_output%tdtsw(is:ie,js:je,:) + &
                            Lw_output%heatra(:,:,:)/seconds_per_day)*  &
                            mask(:,:,:)
          else
            Rad_output%tdtlw(is:ie,js:je,:) =   &
                             Lw_output%heatra(:,:,:)/seconds_per_day
            Rad_output%tdt_rad (is:ie,js:je,:) =  &
                                  (Rad_output%tdtsw(is:ie,js:je,:) +   &
                                   Lw_output%heatra(:,:,:)/  &
                                   seconds_per_day)
          endif
          if (do_clear_sky_pass) then
            Rad_output%tdtsw_clr(is:ie,js:je,:) =   &
                                          Sw_output%hswcf(:,:,:)/  &
                                          seconds_per_day
            if (present(mask)) then
              Rad_output%tdt_rad_clr(is:ie,js:je,:) =    &
                         (Rad_output%tdtsw_clr(is:ie,js:je,:) +  &
                          Lw_output%heatracf(:,:,:)/seconds_per_day)* &
                          mask(:,:,:)
            else
              Rad_output%tdt_rad_clr(is:ie,js:je,:) =    &
                           (Rad_output%tdtsw_clr(is:ie,js:je,:) +    &
                            Lw_output%heatracf(:,:,:)/  &
                            seconds_per_day)
            endif
          endif

          Rad_output%flux_sw_surf(is:ie,js:je) =   &
                                         Sw_output%dfsw(:,:,kmax+1) - &
                                         Sw_output%ufsw(:,:,kmax+1)
          Rad_output%flux_lw_surf(is:ie,js:je) =    &
                       STEFAN*Atm_inp_mdl%temp(:,:,KMAX+1)**4 -   &
                       Lw_output%flxnet(:,:,kmax+1)
        else
          Rad_output%tdtsw(is:ie,js:je,:) = Fsrad_output%tdtsw(:,:,:)
          if (present(mask)) then
            Rad_output%tdt_rad (is:ie,js:je,:) =   &
                               (Rad_output%tdtsw(is:ie,js:je,:) + &
                                Fsrad_output%tdtlw (:,:,:))*mask(:,:,:)
          else
            Rad_output%tdt_rad (is:ie,js:je,:) =   &
                               (Rad_output%tdtsw(is:ie,js:je,:) +   &
                                Fsrad_output%tdtlw (:,:,:))
          endif
          if (do_clear_sky_pass) then
            Rad_output%tdtsw_clr(is:ie,js:je,:) =    &
                                          Fsrad_output%tdtsw_clr(:,:,:)
            if (present(mask)) then
              Rad_output%tdt_rad_clr(is:ie,js:je,:) =    &
                         (Rad_output%tdtsw_clr(is:ie,js:je,:) +  &
                          Fsrad_output%tdtlw_clr(:,:,:))*mask(:,:,:)
            else
              Rad_output%tdt_rad_clr(is:ie,js:je,:) =    &
                         (Rad_output%tdtsw_clr(is:ie,js:je,:) +    &
                          Fsrad_output%tdtlw_clr(:,:,:))
            endif
          endif
          Rad_output%flux_sw_surf(is:ie,js:je) =    &
                                             Fsrad_output%swdns(:,:) - &
                                             Fsrad_output%swups(:,:)
          Rad_output%flux_lw_surf(is:ie,js:je) = Fsrad_output%lwdns(:,:)
        endif ! (do_sea_esf_rad)
      endif ! (running_gcm)

!---------------------------------------------------------------------
!    if the radiation grid does not coincide with the model grid,
!    call deallocate_calc_arrays to deallocate the array space used
!    by the derived type variables spanning radiation-space.
!---------------------------------------------------------------------
      if ( .not. all_column_radiation .or.    &
           .not. all_level_radiation) then
        call deallocate_calc_arrays (Atmos_input, Astro, Rad_gases, &
                                     Cldrad_props, Lw_output_rd, &
                                     Sw_output_rd, Fsrad_output_rd)
      endif

!---------------------------------------------------------------------




end subroutine radiation_calc



!#####################################################################

subroutine model_to_radiation_grid (is, ie, js, je, Atmos_input,&
                                    Astro, Rad_gases,  Cldrad_props, &
                                    Cld_diagnostics, Atmos_input_rd, &
                                    Astro_rd,    &
                                    Rad_gases_rd,  Cldrad_props_rd)

!----------------------------------------------------------------------
!   model_to_radiation_grid maps the input fields to the radiation
!   package from the model grid to the radiation grid, if they differ.
!----------------------------------------------------------------------

integer,                      intent(in)  ::  is, ie, js, je
type(atmos_input_type),       intent(in)  ::  Atmos_input
type(astronomy_type),         intent(in)  ::  Astro         
type(radiative_gases_type),   intent(in)  ::  Rad_gases       
type(cldrad_properties_type), intent(in)  ::  Cldrad_props    
type(cld_diagnostics_type),   intent(in)  ::  Cld_diagnostics  
type(atmos_input_type),       intent(inout) ::  Atmos_input_rd
type(astronomy_type),         intent(inout) ::  Astro_rd         
type(radiative_gases_type),   intent(inout) ::  Rad_gases_rd       
type(cldrad_properties_type), intent(inout) ::  Cldrad_props_rd    

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      Atmos_input       atmos_input_type variable containing atmos-
!                        pheric input data for the radiation package 
!                        on the model grid
!      Astro             astronomy_type variable containing astronomical
!                        input data for the radiation package on the 
!                        model grid
!      Rad_gases         radiative_gases_type variable containing rad-
!                        iative gas input data for the radiation package
!                        on the model grid
!      Cldrad_props      cldrad_properties_type variable containing 
!                        cloud radiative property input data for the 
!                        radiation package on the model grid
!
!   intent(out) variables:
!
!      Atmos_input_rd    atmos_input_type variable containing atmos-
!                        pheric input data for the radiation package 
!                        on the radiation grid
!      Astro_rd          astronomy_type variable containing astronomical
!                        input data for the radiation package on the 
!                        radiation grid
!      Rad_gases_rd      radiative_gases_type variable containing rad-
!                        iative gas input data for the radiation package
!                        on the radiation grid
!      Cldrad_props_rd   cldrad_properties_type variable containing 
!                        cloud radiative property input data for the 
!                        radiation package on the radiation grid
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    this routine is responsible for mapping the model data onto the
!    radiation grid, when the two grids differ. The current code 
!    only provides a mapping in the vertical which could be used
!    to eliminate radiation calculation at points above 
!    topmost_radiation_level. if other mappings are developed, they
!    could be implemented from this subroutine, similarly to the way
!    drop_upper_levels is.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    map the atmos_input_type variables to the radiation grid.
!---------------------------------------------------------------------
      if (drop_upper_levels) then
        Atmos_input_rd%press(:,:,ksrad:kerad+1) =    &
                                        Atmos_input%press(:,:,ks:ke+1)
        Atmos_input_rd%temp (:,:,ksrad:kerad+1) =    &
                                        Atmos_input%temp (:,:,ks:ke+1)
        Atmos_input_rd%rh2o (:,:,ksrad:kerad  ) =    &
                                        Atmos_input%rh2o (:,:,ks:ke  )
        Atmos_input_rd%rel_hum(:,:,ksrad:kerad  ) =    &
                                        Atmos_input%rel_hum(:,:,ks:ke  )
        Atmos_input_rd%pflux(:,:,ksrad:kerad+1) =     &
                                        Atmos_input%pflux (:,:,ks:ke+1)
        if (do_sea_esf_rad) then
          Atmos_input_rd%cloud_water(:,:,ksrad:kerad) =    &
                                     Atmos_input%cloud_water(:,:,ks:ke)
          Atmos_input_rd%cloud_ice(:,:,ksrad:kerad) =      &
                                     Atmos_input%cloud_ice(:,:,ks:ke)
          Atmos_input_rd%cloudtemp(:,:,ksrad:kerad) =      &
                                     Atmos_input%cloudtemp(:,:,ks:ke)
          Atmos_input_rd%tflux(:,:,ksrad:kerad+1) =    &
                                     Atmos_input%tflux (:,:,ks:ke+1)
          Atmos_input_rd%deltaz(:,:,ksrad:kerad) =    &
                                     Atmos_input%deltaz (:,:,ks:ke)
          Atmos_input_rd%cloudvapor (:,:,ksrad:kerad) =    &
                                     Atmos_input%cloudvapor  (:,:,ks:ke)
          Atmos_input_rd%clouddeltaz(:,:,ksrad:kerad) =    &
                                     Atmos_input%clouddeltaz (:,:,ks:ke)
        endif
        Atmos_input_rd%land = Atmos_input%land
        Atmos_input_rd%asfc = Atmos_input%asfc
        Atmos_input_rd%psfc = Atmos_input%psfc
        Atmos_input_rd%tsfc = Atmos_input%tsfc

!---------------------------------------------------------------------
!    map the radiative_gases_type variables to the radiation grid.
!---------------------------------------------------------------------
        Rad_gases_rd%qo3(:,:,ksrad:kerad) = Rad_gases%qo3(:,:,ks:ke)
        Rad_gases_rd%rrvch4               = Rad_gases%rrvch4
        Rad_gases_rd%rrvn2o               = Rad_gases%rrvn2o
        Rad_gases_rd%rrvco2               = Rad_gases%rrvco2
        Rad_gases_rd%rrvf11               = Rad_gases%rrvf11
        Rad_gases_rd%rrvf12               = Rad_gases%rrvf12
        Rad_gases_rd%rrvf113              = Rad_gases%rrvf113
        Rad_gases_rd%rrvf22               = Rad_gases%rrvf22
!       Rad_gases_rd%do_co2               = Rad_gases%do_co2
!       Rad_gases_rd%do_ch4_n2o           = Rad_gases%do_ch4_n2o
!       Rad_gases_rd%do_cfc               = Rad_gases%do_cfc
        Rad_gases_rd%time_varying_co2     = Rad_gases%time_varying_co2
        Rad_gases_rd%time_varying_ch4     = Rad_gases%time_varying_ch4
        Rad_gases_rd%time_varying_n2o     = Rad_gases%time_varying_n2o
        Rad_gases_rd%time_varying_f11     = Rad_gases%time_varying_f11
        Rad_gases_rd%time_varying_f12     = Rad_gases%time_varying_f12
        Rad_gases_rd%time_varying_f113     = Rad_gases%time_varying_f113
        Rad_gases_rd%time_varying_f22     = Rad_gases%time_varying_f22

!---------------------------------------------------------------------
!    map the astronomy_type variables to the radiation grid.
!---------------------------------------------------------------------
        Astro_rd%cosz          = Astro%cosz
        Astro_rd%fracday       = Astro%fracday
        Astro_rd%solar         = Astro%solar
        Astro_rd%rrsun         = Astro%rrsun
        Astro_rd%do_diurnal    = Astro%do_diurnal
        Astro_rd%do_annual     = Astro%do_annual
        Astro_rd%do_daily_mean = Astro%do_daily_mean
        Astro_rd%solar_constant= Astro%solar_constant

!---------------------------------------------------------------------
!    map the cldrad_properties_type variables to the radiation grid. 
!---------------------------------------------------------------------
        call error_mesg ( 'radiation_driver_mod', &
          'must break out components of Cldrad_props and map to &
	       &radiation grid.', FATAL)
        Cldrad_props_rd = Cldrad_props

        call error_mesg ('radiation_driver_mod', &
         'option to not calculate radiation at all model levels not &
	    &yet implemented', FATAL)

!---------------------------------------------------------------------
!    if another mapping from model to radiation grid is desired, the
!    code for it must be supplied here.
!---------------------------------------------------------------------
      else   ! (drop_upper_levels)
        call error_mesg ('radiation_driver_mod', &
          'must supply code to map data from model grid to radiation &
           &grid.', FATAL)
      endif  !  (drop_upper_levels)

!---------------------------------------------------------------------

end subroutine model_to_radiation_grid


!####################################################################

subroutine radiation_to_model_grid (is, ie, js, je,          &
                                    Lw_output_rd, Sw_output_rd, & 
                                    Fsrad_output_rd, Lw_output,   &
                                    Sw_output, Fsrad_output)

!--------------------------------------------------------------------
!    radiation_to_model_grid maps the radiation output fields from the
!    radiation grid to the model grid, when the grids differ.
!--------------------------------------------------------------------

!---------------------------------------------------------------------- 
integer,                 intent(in)            :: is, ie, js, je
type(lw_output_type),    intent(in),  optional :: Lw_output_rd
type(sw_output_type),    intent(in),  optional :: Sw_output_rd
type(fsrad_output_type), intent(in),  optional :: Fsrad_output_rd
type(lw_output_type),    intent(inout), optional :: Lw_output
type(sw_output_type),    intent(inout), optional :: Sw_output
type(fsrad_output_type), intent(inout), optional :: Fsrad_output

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!
!   intent(in), optional variables:
!
!      Lw_output_rd      lw_output_type variable containing output from 
!                        the longwave radiation code of the
!                        sea_esf_rad package, on the radiation grid
!      Sw_output_rd      sw_output_type variable containing output from 
!                        the shortwave radiation code of the
!                        sea_esf_rad package, on the radiation grid
!      Fsrad_output_rd   fsrad_output_type variable containing 
!                        output from the original_fms_rad radiation
!                        package, on the radiation grid
!
!   intent(out), optional variables:
!
!      Lw_output         lw_output_type variable containing output from 
!                        the longwave radiation code of the
!                        sea_esf_rad package, on the model grid
!      Sw_output         sw_output_type variable containing output from 
!                        the shortwave radiation code of the
!                        sea_esf_rad package, on the model grid
!      Fsrad_output      fsrad_output_type variable containing 
!                        output from the original_fms_rad radiation
!                        package, on the model grid
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    code needed here to map the components of the derived types from
!    the radiation grid to the model grid.
!--------------------------------------------------------------------
      call error_mesg ('radiation_driver_mod', &
       'must supply code to map data from radiation grid to model &
           &grid.', FATAL)
      if (drop_upper_levels) then
        if (present(Sw_output_rd)) then
          Sw_output = Sw_output_rd     
          Lw_output = Lw_output_rd     
        else if (present(Fsrad_output_rd)) then
          Fsrad_output = Fsrad_output_rd
        else
          call error_mesg ('radiation_driver_mod', &
            'no output data supplied to radiation_to_model_grid', FATAL)
        endif
      else     ! (drop_upper_levels)
        call error_mesg ('radiation_driver_mod', &
           ' currently only drop_upper_levels is available as option',&
            FATAL)
      endif    !  (drop_upper_levels)

!----------------------------------------------------------------------

end subroutine radiation_to_model_grid



!#######################################################################

subroutine deallocate_calc_arrays (Atmos_input, Astro, Rad_gases, &
                                   Cldrad_props, Lw_output_rd, &
                                   Sw_output_rd, Fsrad_output_rd)

!---------------------------------------------------------------------
!    deallocate_calc_arrays deallocates the array space of stack-
!    resident derived-type variables in subroutine radiation_calc.
!----------------------------------------------------------------------

type(atmos_input_type),       intent(in)   :: Atmos_input
type(astronomy_type),         intent(in)   :: Astro
type(radiative_gases_type),   intent(in)   :: Rad_gases
type(cldrad_properties_type), intent(in)   :: Cldrad_props
type(lw_output_type),         intent(in)   :: Lw_output_rd
type(sw_output_type),         intent(in)   :: Sw_output_rd
type(fsrad_output_type),      intent(in)   :: Fsrad_output_rd


!-------------------------------------------------------------------
!    deallocate components of Rad_gases.
!-------------------------------------------------------------------
      deallocate (Rad_gases%qo3)

!-------------------------------------------------------------------
!    deallocate components of Astro.
!-------------------------------------------------------------------
      deallocate (Astro%solar)
      deallocate (Astro%cosz )
      deallocate (Astro%fracday)

!-------------------------------------------------------------------
!    deallocate components of Lw_output.
!-------------------------------------------------------------------
      if (do_sea_esf_rad) then
        deallocate (Lw_output_rd%heatra    )
        deallocate (Lw_output_rd%flxnet    )
        if (Rad_control%do_totcld_forcing) then
          deallocate (Lw_output_rd%heatracf  )
          deallocate (Lw_output_rd%flxnetcf  )
        endif
      endif

!-------------------------------------------------------------------
!    deallocate components of Atmos_input.
!-------------------------------------------------------------------
      deallocate (Atmos_input%press    )
      deallocate (Atmos_input%temp     )
      deallocate (Atmos_input%rh2o     )
      deallocate (Atmos_input%rel_hum  )
      deallocate (Atmos_input%pflux    )
      deallocate (Atmos_input%tflux    )
      deallocate (Atmos_input%deltaz   )
      deallocate (Atmos_input%psfc     )
      deallocate (Atmos_input%tsfc     )
      deallocate (Atmos_input%asfc     )
      deallocate (Atmos_input%land     )
      deallocate (Atmos_input%cloud_water)
      deallocate (Atmos_input%cloud_ice)
      deallocate (Atmos_input%cloudtemp)
      deallocate (Atmos_input%cloudvapor )
      deallocate (Atmos_input%clouddeltaz)

!-------------------------------------------------------------------
!    deallocate components of Sw_output.
!-------------------------------------------------------------------
      if (do_sea_esf_rad) then
        deallocate (Sw_output_rd%dfsw     )
        deallocate (Sw_output_rd%ufsw     )
        deallocate (Sw_output_rd%fsw     )
        deallocate (Sw_output_rd%hsw     )
        if (Rad_control%do_totcld_forcing) then
          deallocate (Sw_output_rd%dfswcf   )
          deallocate (Sw_output_rd%ufswcf   )
          deallocate (Sw_output_rd%fswcf   )
          deallocate (Sw_output_rd%hswcf   )
        endif
      endif

!-------------------------------------------------------------------
!    deallocate components of Cldrad_props.
!-------------------------------------------------------------------
      deallocate (Cldrad_props%camtsw    )
      deallocate (Cldrad_props%cmxolw    )
      deallocate (Cldrad_props%crndlw    )
      deallocate (Cldrad_props%ncldsw    )
      deallocate (Cldrad_props%nmxolw     )
      deallocate (Cldrad_props%nrndlw    )

      deallocate (Cldrad_props%emmxolw   )
      deallocate (Cldrad_props%emrndlw   )

      if ( associated(Cldrad_props%cldext) ) then
        deallocate (Cldrad_props%cldext    )
        deallocate (Cldrad_props%cldasymm  )
        deallocate (Cldrad_props%cldsct    )
      endif
      if ( associated(Cldrad_props%cldemiss) ) then
        deallocate (Cldrad_props%abscoeff  )
        deallocate (Cldrad_props%cldemiss  )
      endif
      if ( associated(Cldrad_props%cvisrfsw) ) then
        deallocate (Cldrad_props%cvisrfsw  )
        deallocate (Cldrad_props%cirabsw  )
        deallocate (Cldrad_props%cirrfsw  )
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Fsrad_output_rd.
!--------------------------------------------------------------------
      if (.not. do_sea_esf_rad) then
        deallocate (Fsrad_output_rd%tdtsw   )
        deallocate (Fsrad_output_rd%tdtlw   )
        deallocate (Fsrad_output_rd%swdns   )
        deallocate (Fsrad_output_rd%swups   )
        deallocate (Fsrad_output_rd%lwdns   )
        deallocate (Fsrad_output_rd%lwups   )
        deallocate (Fsrad_output_rd%swin    )
        deallocate (Fsrad_output_rd%swout   )
        deallocate (Fsrad_output_rd%olr     )
        if (do_clear_sky_pass) then
          deallocate (Fsrad_output_rd%tdtsw_clr )
          deallocate (Fsrad_output_rd%tdtlw_clr )
          deallocate (Fsrad_output_rd%swdns_clr )
          deallocate (Fsrad_output_rd%swups_clr )
          deallocate (Fsrad_output_rd%lwdns_clr )
          deallocate (Fsrad_output_rd%lwups_clr )
          deallocate (Fsrad_output_rd%swin_clr  )
          deallocate (Fsrad_output_rd%swout_clr )
          deallocate (Fsrad_output_rd%olr_clr   )
        endif 
      endif 

!---------------------------------------------------------------------


end subroutine deallocate_calc_arrays 




!######################################################################

subroutine update_rad_fields (is, ie, js, je, Time_diag, Astro2,   &
                              Sw_output, Astro, Rad_output, tdt,  &
			      coszen, flux_sw, flux_lw, flux_ratio)

!---------------------------------------------------------------------
!    update_rad_fields defines the current radiative heating rate, 
!    surface long and short wave fluxes and cosine of zenith angle
!    to be returned to physics_driver, including renormalization 
!    effects when that option is activated.
!--------------------------------------------------------------------

integer,                 intent(in)    ::  is, ie, js, je
type(time_type),         intent(in)    ::  Time_diag
type(astronomy_type),    intent(in)    ::  Astro2
type(sw_output_type),    intent(in)    ::  Sw_output
type(astronomy_type),    intent(inout) ::  Astro
type(rad_output_type),   intent(inout) ::  Rad_output
real,  dimension(:,:,:), intent(inout) ::  tdt
real,  dimension(:,:),   intent(out)   ::  flux_sw, flux_lw, coszen, &
                                           flux_ratio

!-------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data 
!                   in the physics_window being integrated
!      Time_diag    time on next timestep, used as stamp for diag-
!                   nostic output  [ time_type  (days, seconds) ]  
!      Astro2       astronomical properties on model grid, valid over 
!                   physics timestep, used when renormalizing sw fluxes
!                   [astronomy_type]
!      Sw_output    shortwave output variables on model grid,
!                   [sw_output_type]     
!
!  intent(inout) variables:
!
!      Astro        astronomical properties on model grid, usually
!                   valid over radiation timestep on entry, on exit are 
!                   valid over model timestep when renormalizing
!                   [astronomy_type]
!      Rad_output   radiation output variables on model grid, valid
!                   on entry over either physics or radiation timestep, 
!                   on exit are valid over physics step when renormal-
!                   izing sw fluxes
!                   [rad_output_type]     
!      tdt          time tendency of temperature on current step;
!                   value on exit has had radiative heating rate
!                   added to the value on entrance
!
!  intent(out) variables:
!
!      flux_sw      net (down-up) sw flux at surface
!      flux_lw      downward lw flux at surface
!      coszen       cosine of the zenith angle 
!      flux_ratio   factor to multiply the radiation step values of 
!                   sw fluxes and heating rates by in order to get
!                   current physics timestep values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables
!
      real, dimension ( size(tdt,1), size(tdt,2), &
                        size(tdt,3) )           ::  tdtlw, tdtlw_clr    
      integer   ::   i, j, k

      if (renormalize_sw_fluxes) then

!----------------------------------------------------------------------
!    if sw fluxes are to be renormalized, save the heating rates, fluxes
!    and solar factor calculated on radiation steps.
!---------------------------------------------------------------------
        if (do_rad) then
          solar_save(is:ie,js:je)  = Astro%solar(:,:)
          dfsw_save(is:ie,js:je,:) = Sw_output%dfsw(:, :,:)
          ufsw_save(is:ie,js:je,:) = Sw_output%ufsw(:, :,:)
          fsw_save(is:ie,js:je,:)  = Sw_output%fsw(:, :,:)
          hsw_save(is:ie,js:je,:)  = Sw_output%hsw(:, :,:)
          flux_sw_surf_save(is:ie,js:je) =    &
                                    Rad_output%flux_sw_surf(is:ie,js:je)
          sw_heating_save(is:ie,js:je,:) =    &
                              Rad_output%tdtsw(is:ie,js:je,:)
          tot_heating_save(is:ie,js:je,:) =    &
                              Rad_output%tdt_rad(is:ie,js:je,:)
          if (do_clear_sky_pass) then
            sw_heating_clr_save(is:ie,js:je,:) =    &
                              Rad_output%tdtsw_clr(is:ie,js:je,:)
            tot_heating_clr_save(is:ie,js:je,:) =    &
                              Rad_output%tdt_rad_clr(is:ie,js:je,:)
            dfswcf_save(is:ie,js:je,:) = Sw_output%dfswcf(:, :,:)
            ufswcf_save(is:ie,js:je,:) = Sw_output%ufswcf(:, :,:)
            fswcf_save(is:ie,js:je,:)  = Sw_output%fswcf(:, :,:)
            hswcf_save(is:ie,js:je,:)  = Sw_output%hswcf(:, :,:)
          endif 

!---------------------------------------------------------------------
!    define the ratio of the solar factor valid over this physics step
!    to that valid over the current radiation timestep.
!---------------------------------------------------------------------
          do j=1,je-js+1
            do i=1,ie-is+1
              if (solar_save(i+is-1,j+js-1) /= 0.0) then
                flux_ratio(i, j) = Astro2%solar(i,j)/   &
                                            solar_save(i+is-1,j+js-1)
              else
                flux_ratio(i,j) = 0.0
              endif

!---------------------------------------------------------------------
!    move the physics-step values(Astro2) to Astro, which will be used 
!    to calculate diagnostics. the radiation_step values (Astro) are no
!    longer needed.
!---------------------------------------------------------------------
              Astro%cosz(i,j) = Astro2%cosz(i,j)
              Astro%fracday(i,j) = Astro2%fracday(i,j)
            end do
          end do
!----------------------------------------------------------------------
!    on non-radiation steps define the ratio of the current solar factor
!    valid for this physics step to that valid for the last radiation 
!    step. 
!----------------------------------------------------------------------
        else 
          do j=1,je-js+1
            do i=1,ie-is+1
              if (solar_save(i+is-1,j+js-1) /= 0.0) then
                flux_ratio(i, j) = Astro%solar(i,j)/   &
                                            solar_save(i+is-1,j+js-1)
              else
                flux_ratio(i,j) = 0.0
              endif
            end do
          end do
        endif  ! (do_rad)

!---------------------------------------------------------------------
!    redefine the total and shortwave heating rates, along with surface
!    sw fluxes, as a result of the difference in solar factor (the 
!    relative earth-sun motion) between the current physics and current
!    radiation timesteps.
!---------------------------------------------------------------------
        tdtlw(:,:,:) = tot_heating_save(is:ie,js:je,:) -    &
                       sw_heating_save(is:ie,js:je,:)
        do k=1, size(tdt,3)
          Rad_output%tdtsw(is:ie,js:je,k) =    &
                       sw_heating_save(is:ie,js:je,k)*flux_ratio(:,:)
        end do
        Rad_output%tdt_rad(is:ie,js:je,:) = tdtlw(:,:,:) +    &
                                       Rad_output%tdtsw(is:ie,js:je,:)
        Rad_output%flux_sw_surf(is:ie,js:je) = flux_ratio(:,:)*    &
                                          flux_sw_surf_save(is:ie,js:je)
        if (do_clear_sky_pass) then
          tdtlw_clr(:,:,:) = tot_heating_clr_save  (is:ie,js:je,:) -&
                             sw_heating_clr_save (is:ie,js:je,:)
          do k=1, size(tdt,3)
            Rad_output%tdtsw_clr(is:ie,js:je,k) =   &
                           sw_heating_clr_save (is:ie,js:je,k)*   &
                           flux_ratio(:,:)
          end do
          Rad_output%tdt_rad_clr(is:ie,js:je,:) = tdtlw_clr(:,:,:) +   &
                               Rad_output%tdtsw_clr(is:ie,js:je,:)
        endif
      else if (all_step_diagnostics) then
!----------------------------------------------------------------------
!    if sw fluxes are to be output on every physics step, save the 
!    heating rates and fluxes calculated on radiation steps.
!---------------------------------------------------------------------
        if (do_rad) then
!         solar_save(is:ie,js:je)  = Astro%solar(:,:)
          dfsw_save(is:ie,js:je,:) = Sw_output%dfsw(:, :,:)
          ufsw_save(is:ie,js:je,:) = Sw_output%ufsw(:, :,:)
          fsw_save(is:ie,js:je,:)  = Sw_output%fsw(:, :,:)
          hsw_save(is:ie,js:je,:)  = Sw_output%hsw(:, :,:)
          flux_sw_surf_save(is:ie,js:je) =    &
                              Rad_output%flux_sw_surf(is:ie,js:je)
          sw_heating_save(is:ie,js:je,:) =    &
                             Rad_output%tdtsw(is:ie,js:je,:)
          tot_heating_save(is:ie,js:je,:) =    &
                               Rad_output%tdt_rad(is:ie,js:je,:)
          if (do_clear_sky_pass) then
            sw_heating_clr_save(is:ie,js:je,:) =    &
                        Rad_output%tdtsw_clr(is:ie,js:je,:)
            tot_heating_clr_save(is:ie,js:je,:) =    &
                          Rad_output%tdt_rad_clr(is:ie,js:je,:)
            dfswcf_save(is:ie,js:je,:) = Sw_output%dfswcf(:, :,:)
            ufswcf_save(is:ie,js:je,:) = Sw_output%ufswcf(:, :,:)
            fswcf_save(is:ie,js:je,:)  = Sw_output%fswcf(:, :,:)
            hswcf_save(is:ie,js:je,:)  = Sw_output%hswcf(:, :,:)
          endif
        endif
      else
        flux_ratio(:,:) = 1.0
      endif  ! (renormalize_sw_fluxes)

!---------------------------------------------------------------------
!    define the value of coszen to be returned to be used in the ocean
!    albedo calculation. 
!---------------------------------------------------------------------
      if (Environment%running_gcm .or.  &
          Environment%running_sa_model) then
        coszen(:,:) = Rad_output%coszen_angle(is:ie, js:je)
      endif

!-------------------------------------------------------------------
!    update radiative tendency and fluxes.
!-------------------------------------------------------------------
      tdt  (:,:,:) = tdt(:,:,:) + Rad_output%tdt_rad (is:ie,js:je,:)
      flux_sw(:,:) =         Rad_output%flux_sw_surf (is:ie,js:je)
      flux_lw(:,:) =         Rad_output%flux_lw_surf (is:ie,js:je)

!-------------------------------------------------------------------
!    save radiative tendency for use in strat cloud scheme. Note that 
!    if strat cloud scheme is not operating then nothing is done by 
!    this routine.
!-------------------------------------------------------------------
      if (Environment%running_gcm) then
        call add_strat_tend (is, ie, js, je,    &
                               Rad_output%tdt_rad (is:ie,js:je,:))

!-------------------------------------------------------------------
!    save longwave tendency for use in edt turbulence scheme. Note that 
!    if edt scheme is not operating then nothing is done by this 
!    routine.               
!-------------------------------------------------------------------
        call edt_tend (is, ie, js, je, Rad_output%tdtlw (is:ie,js:je,:))
!-------------------------------------------------------------------
!    save longwave tendency for use in entrainment turbulence scheme. 
!    Note that if entrainment scheme is not operating then nothing 
!    is done by this routine.               
!-------------------------------------------------------------------
        call entrain_tend (is, ie, js, je, Rad_output%tdtlw (is:ie,js:je,:))
!--------------------------------------------------------------------

      endif

!--------------------------------------------------------------------



end subroutine update_rad_fields 



!####################################################################

subroutine radiation_netcdf (is, ie, js, je, Time_diag, lat,   &
                             ts, asfc, flux_ratio, Astro, Rad_output, &
			     Lw_output, Sw_output, Fsrad_output, mask) 

!--------------------------------------------------------------------
!    radiation_netcdf produces netcdf output and global and hemispheric
!    integrals of radiation package variables.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                 intent(in)             :: is, ie, js, je
type(time_type),         intent(in)             :: Time_diag
real,dimension(:,:),     intent(in)             :: lat, ts
real,dimension(:,:),     intent(in)             :: asfc, flux_ratio
type(astronomy_type),    intent(in)             :: Astro
type(rad_output_type),   intent(in)             :: Rad_output
type(lw_output_type),    intent(in), optional   :: Lw_output
type(fsrad_output_type), intent(in), optional   :: Fsrad_output
type(sw_output_type),    intent(in), optional   :: Sw_output
real,dimension(:,:,:),   intent(in), optional   :: mask
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data 
!                   in the physics_window being integrated
!      Time_diag    time on next timestep, used as stamp for diagnostic 
!                   output  [ time_type  (days, seconds) ]  
!      lat          latitude of model points  [ radians ]
!      ts           surface temperature  [ deg K ]
!      asfc         surface albedo  [ dimensionless ]
!      flux_ratio   renormalization factor for sw fluxes and heating 
!                   rates [ dimensionless ]
!      Astro        cosine of zenith angle [ dimensionless ]
!      Rad_output   rad_output_type variable containing radiation 
!                   output fields
!
!
!    intent(in) optional variables:
!
!      Lw_output    lw_output_type variable containing output from 
!                   the longwave radiation code of the
!                   sea_esf_rad package, on the model grid
!      Sw_output    sw_output_type variable containing output from 
!                   the shortwave radiation code of the
!                   sea_esf_rad package, on the model grid
!      Fsrad_output fsrad_output_type variable containing 
!                   output from the original_fms_rad radiation
!                   package, on the model grid
!      mask         present when running eta vertical coordinate,
!                   mask to remove points below ground
!        
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables

      real, dimension (ie-is+1,je-js+1) ::           & 
                                                swin, swout, olr, &
                                                swups, swdns, lwups, &
                                                lwdns, swin_clr,   &
                                                swout_clr, olr_clr, &
                                                swups_clr, swdns_clr,&
                                                lwups_clr, lwdns_clr

      real, dimension (ie-is+1,je-js+1, kmin:kmax) ::    &
                                                tdtlw, tdtlw_clr,&
                                                hsw, hswcf

      real, dimension (ie-is+1,je-js+1, kmin:kmax+1) ::   &
                                                dfsw, ufsw,  &
                                                dfswcf, ufswcf,&
                                                fsw, fswcf

      integer           :: j, k
      integer           :: ipass
      logical           :: used
      integer           :: iind, jind

!---------------------------------------------------------------------
!    if sw flux renormalization is active, modify the fluxes calculated
!    on the last radiation step by the normalization factor based on
!    the difference in solar factor between the current model step and
!    the current radiation step.
!----------------------------------------------------------------------
      if (renormalize_sw_fluxes) then
        do k=1, kmax         
          hsw(:,:,k) = hsw_save(is:ie,js:je,k)*flux_ratio(:,:)
        end do
        do k=1, kmax+1             
          dfsw(:,:,k) = dfsw_save(is:ie,js:je,k)*flux_ratio(:,:)
          ufsw(:,:,k) = ufsw_save(is:ie,js:je,k)*flux_ratio(:,:)
          fsw(:,:,k) = fsw_save(is:ie,js:je,k)*flux_ratio(:,:)
        end do
        if (do_clear_sky_pass) then
          do k=1, kmax            
            hswcf(:,:,k) = hswcf_save(is:ie,js:je,k)*flux_ratio(:,:)
          end do
          do k=1, kmax+1            
            dfswcf(:,:,k) = dfswcf_save(is:ie,js:je,k)*flux_ratio(:,:)
            ufswcf(:,:,k) = ufswcf_save(is:ie,js:je,k)*flux_ratio(:,:)
            fswcf(:,:,k) = fswcf_save(is:ie,js:je,k)*flux_ratio(:,:)
          end do
        endif

!----------------------------------------------------------------------
!    if renormalization is not active and this is a radiation step
!    (i.e., diagnostics desired), define the variables to be output as
!    the values present in Sw_output.
!---------------------------------------------------------------------
      else if (do_rad .and. do_sea_esf_rad) then
        do k=1, kmax            
          hsw(:,:,k) = Sw_output%hsw(:,:,k)
        end do
        do k=1, kmax+1             
          dfsw(:,:,k) = Sw_output%dfsw(:,:,k)
          ufsw(:,:,k) = Sw_output%ufsw(:,:,k)
          fsw(:,:,k) = Sw_output%fsw(:,:,k)
        end do
        if (do_clear_sky_pass) then
          do k=1, kmax             
            hswcf(:,:,k) = Sw_output%hswcf(:,:,k)
          end do
          do k=1, kmax+1            
            dfswcf(:,:,k) = Sw_output%dfswcf(:,:,k)
            ufswcf(:,:,k) = Sw_output%ufswcf(:,:,k)
            fswcf(:,:,k) = Sw_output%fswcf(:,:,k)
          end do
        endif

!----------------------------------------------------------------------
!    if renormalization is not active and this is not a radiation step
!    but all_step_diagnostics is activated (i.e., diagnostics desired),
!    define the variables to be output as the values previously saved
!    in the xxx_save variables.
!---------------------------------------------------------------------
      else if (do_sea_esf_rad .and. all_step_diagnostics) then
        do k=1, kmax
          hsw(:,:,k) = hsw_save(is:ie,js:je,k)
        end do
        do k=1, kmax+1
          dfsw(:,:,k) = dfsw_save(is:ie,js:je,k)
          ufsw(:,:,k) = ufsw_save(is:ie,js:je,k)
          fsw(:,:,k) = fsw_save(is:ie,js:je,k)
        end do
        if (do_clear_sky_pass) then
          do k=1, kmax
            hswcf(:,:,k) = hswcf_save(is:ie,js:je,k)
          end do
          do k=1, kmax+1
            dfswcf(:,:,k) = dfswcf_save(is:ie,js:je,k)
            ufswcf(:,:,k) = ufswcf_save(is:ie,js:je,k)
            fswcf(:,:,k) = fswcf_save(is:ie,js:je,k)
          end do
        endif
      endif

!---------------------------------------------------------------------
!    define the sw diagnostic arrays.
!---------------------------------------------------------------------
      if (renormalize_sw_fluxes .or. do_rad .or.    &
                                            all_step_diagnostics) then
        if (do_sea_esf_rad) then
          swin (:,:) = dfsw(:,:,1)
          swout(:,:) = ufsw(:,:,1)
          swups(:,:) = ufsw(:,:,kmax+1)
          swdns(:,:) = dfsw(:,:,kmax+1)
          if (do_clear_sky_pass) then
            swin_clr (:,:) = dfswcf(:,:,1)
            swout_clr(:,:) = ufswcf(:,:,1)
            swups_clr(:,:) = ufswcf(:,:,kmax+1)
            swdns_clr(:,:) = dfswcf(:,:,kmax+1)
          endif
        else   ! original fms rad
          swin (:,:) = Fsrad_output%swin(:,:)               
          swout(:,:) = Fsrad_output%swout(:,:)         
          swups(:,:) = Fsrad_output%swups(:,:)
          swdns(:,:) = Fsrad_output%swdns(:,:)
          if (do_clear_sky_pass) then
            swin_clr (:,:) = Fsrad_output%swin_clr(:,:)               
            swout_clr(:,:) = Fsrad_output%swout_clr(:,:)         
            swups_clr(:,:) = Fsrad_output%swups_clr(:,:)
            swdns_clr(:,:) = Fsrad_output%swdns_clr(:,:)
          endif
        endif  ! do_sea_esf_rad

!---------------------------------------------------------------------
!   send standard sw diagnostics to diag_manager.
!---------------------------------------------------------------------
        if (do_clear_sky_pass) then
          ipass = 2
        else
          ipass = 1
        endif

!------- sw tendency -----------
        if (id_tdt_sw(ipass) > 0 ) then
          used = send_data (id_tdt_sw(ipass), Rad_output%tdtsw(is:ie,js:je,:),   &
                            Time_diag, is, js, 1, rmask=mask )
        endif


!------- incoming sw flux toa -------
        if (id_swdn_toa(ipass) > 0 ) then
          used = send_data (id_swdn_toa(ipass), swin,   &
                            Time_diag, is, js )
        endif

!------- outgoing sw flux toa -------
        if (id_swup_toa(ipass) > 0 ) then
          used = send_data (id_swup_toa(ipass), swout,    &
                            Time_diag, is, js )
        endif


!------- upward sw flux surface -------
        if (id_swup_sfc(ipass) > 0 ) then
          used = send_data (id_swup_sfc(ipass), swups,    &
                            Time_diag, is, js )
        endif

!------- downward sw flux surface -------
        if (id_swdn_sfc(ipass) > 0 ) then
          used = send_data (id_swdn_sfc(ipass), swdns,   &
                            Time_diag, is, js )
        endif

!----------------------------------------------------------------------
!    now pass clear-sky diagnostics, if they have been calculated.
!----------------------------------------------------------------------
        if (do_clear_sky_pass) then
          ipass = 1

!------- sw tendency -----------
          if (id_tdt_sw(ipass) > 0 ) then
            used = send_data (id_tdt_sw(ipass), Rad_output%tdtsw_clr(is:ie,js:je,:),  &
                              Time_diag, is, js, 1, rmask=mask )
          endif

!------- incoming sw flux toa -------
          if (id_swdn_toa(ipass) > 0 ) then
            used = send_data (id_swdn_toa(ipass), swin_clr,    &
                              Time_diag, is, js )
          endif

!------- outgoing sw flux toa -------
          if (id_swup_toa(ipass) > 0 ) then
            used = send_data (id_swup_toa(ipass), swout_clr,  &
                              Time_diag, is, js )
          endif

!------- upward sw flux surface -------
          if (id_swup_sfc(ipass) > 0 ) then
            used = send_data (id_swup_sfc(ipass), swups_clr,   &
                              Time_diag, is, js )
          endif

!------- downward sw flux surface -------
          if (id_swdn_sfc(ipass) > 0 ) then
            used = send_data (id_swdn_sfc(ipass), swdns_clr,    &
                              Time_diag, is, js )
          endif
        endif  ! (do_clear_sky_pass)

!-----------------------------------------------------------------------
!    send cloud-forcing-independent diagnostics to diagnostics manager.
!-----------------------------------------------------------------------

!------- surface albedo  -------------------------
        if ( id_alb_sfc > 0 ) then
          used = send_data ( id_alb_sfc, 100.*asfc, Time_diag, is, js )
!!BUGFIX:
!         used = send_data ( id_alb_sfc, 100.*asfc, Time_diag, is, js )
        endif

!------- cosine of zenith angle ----------------
        if ( id_cosz > 0 ) then
          used = send_data ( id_cosz, Astro%cosz, Time_diag, is, js )
        endif

!------- daylight fraction  --------------
        if ( id_fracday > 0 ) then
          used = send_data (id_fracday, Astro%fracday, Time_diag,   &
                            is, js )
        endif
      endif   ! (renormalize_sw_fluxes .or. do_rad .or. all_step_diagnostics)

!---------------------------------------------------------------------
!    define the longwave diagnostic arrays for the sea-esf radiation 
!    package.  convert to mks units.
!---------------------------------------------------------------------
      if (do_sea_esf_rad) then
        if (do_rad) then
          olr  (:,:)   = Lw_output%flxnet(:,:,1)
          lwups(:,:)   =   STEFAN*ts(:,:  )**4
          lwdns(:,:)   = lwups(:,:) - Lw_output%flxnet(:,:,kmax+1)
          tdtlw(:,:,:) = Lw_output%heatra(:,:,:)/ seconds_per_day

          if (do_clear_sky_pass) then
            olr_clr  (:,:)   = Lw_output%flxnetcf(:,:,1)
            lwups_clr(:,:)   =              STEFAN*ts(:,:  )**4
            lwdns_clr(:,:)   = lwups_clr(:,:) -    & 
                               Lw_output%flxnetcf(:,:,kmax+1)
            tdtlw_clr(:,:,:) = Lw_output%heatracf(:,:,:)/seconds_per_day
          endif

!---------------------------------------------------------------------
!    if diagnostics are desired on all physics steps, save the arrays 
!    for later use.
!---------------------------------------------------------------------
          if (all_step_diagnostics) then
            olr_save  (is:ie,js:je)   = olr(:,:)
            lwups_save(is:ie,js:je)   = lwups(:,:)
            lwdns_save(is:ie,js:je)   = lwdns(:,:)
            tdtlw_save(is:ie,js:je,:) = tdtlw(:,:,:)

            if (do_clear_sky_pass) then
              olr_clr_save  (is:ie,js:je)   = olr_clr(:,:)
              lwups_clr_save(is:ie,js:je)   = lwups_clr(:,:)
              lwdns_clr_save(is:ie,js:je)   = lwdns_clr(:,:)
              tdtlw_clr_save(is:ie,js:je,:) = tdtlw_clr(:,:,:)
            endif
           endif

!---------------------------------------------------------------------
!    if this is not a radiation step, but diagnostics are desired,
!    define the fields from the xxx_save variables.
!---------------------------------------------------------------------
         else if (all_step_diagnostics) then  ! (do_rad)
           olr(:,:)     = olr_save  (is:ie,js:je)
           lwups(:,:)   = lwups_save(is:ie,js:je)
           lwdns(:,:)   = lwdns_save(is:ie,js:je)
           tdtlw(:,:,:) = tdtlw_save(is:ie,js:je,:)

           if (do_clear_sky_pass) then
             olr_clr(:,:)     = olr_clr_save  (is:ie,js:je)
             lwups_clr(:,:)   = lwups_clr_save(is:ie,js:je)
             lwdns_clr(:,:)   = lwdns_clr_save(is:ie,js:je)
             tdtlw_clr(:,:,:) = tdtlw_clr_save(is:ie,js:je,:)
           endif
         endif

!---------------------------------------------------------------------
!    on radiation steps, define the longwave diagnostic arrays for the
!    original_fms_rad package.        
!---------------------------------------------------------------------
      else   ! original fms rad
        if (do_rad) then
          olr  (:,:)   = Fsrad_output%olr(:,:)
          lwups(:,:)   = Fsrad_output%lwups(:,:)
          lwdns(:,:)   = Fsrad_output%lwdns(:,:)
          tdtlw(:,:,:) = Fsrad_output%tdtlw(:,:,:)

          if (do_clear_sky_pass) then
            olr_clr  (:,:)   = Fsrad_output%olr_clr(:,:)
            lwups_clr(:,:)   = Fsrad_output%lwups_clr(:,:)
            lwdns_clr(:,:)   = Fsrad_output%lwdns_clr(:,:)
            tdtlw_clr(:,:,:) = Fsrad_output%tdtlw_clr(:,:,:)
          endif
        endif
      endif  ! do_sea_esf_rad

      if (do_rad .or. all_step_diagnostics) then
!---------------------------------------------------------------------
!   send standard lw diagnostics to diag_manager.
!---------------------------------------------------------------------
        if (do_clear_sky_pass) then
          ipass = 2
        else
          ipass = 1
        endif

!------- lw tendency -----------
        if (id_tdt_lw(ipass) > 0 ) then
          used = send_data (id_tdt_lw(ipass), tdtlw,    &
                            Time_diag, is, js, 1, rmask=mask )
        endif

!------- outgoing lw flux toa (olr) -------
        if (id_olr(ipass) > 0 ) then
          used = send_data (id_olr(ipass), olr,    &
                            Time_diag, is, js )
        endif

!------- upward lw flux surface -------
        if ( id_lwup_sfc(ipass) > 0 ) then
          used = send_data (id_lwup_sfc(ipass), lwups,    &
                            Time_diag, is, js )
        endif

!------- downward lw flux surface -------
        if (id_lwdn_sfc(ipass) > 0 ) then
          used = send_data (id_lwdn_sfc(ipass), lwdns,    &
                            Time_diag, is, js )
        endif

!----------------------------------------------------------------------
!    now pass clear-sky diagnostics, if they have been calculated.
!----------------------------------------------------------------------
        if (do_clear_sky_pass) then
          ipass = 1

!------- lw tendency -----------
          if (id_tdt_lw(ipass) > 0 ) then
            used = send_data (id_tdt_lw(ipass), tdtlw_clr,    &
                              Time_diag, is, js, 1, rmask=mask )
          endif

!------- outgoing lw flux toa (olr) -------
          if (id_olr(ipass) > 0 ) then
            used = send_data (id_olr(ipass), olr_clr,   &
                              Time_diag, is, js )
          endif

!------- upward lw flux surface -------
          if (id_lwup_sfc(ipass) > 0 ) then
            used = send_data (id_lwup_sfc(ipass), lwups_clr,   &
                              Time_diag, is, js )
          endif

!------- downward lw flux surface -------
          if (id_lwdn_sfc(ipass) > 0 ) then
            used = send_data (id_lwdn_sfc(ipass), lwdns_clr,   &
                              Time_diag, is, js )
          endif
        endif  ! (do_clear_sky_pass)
      endif  ! (do_rad .or. all_step_diagnostics)

!--------------------------------------------------------------------
!    now define various diagnostic integrals.
!--------------------------------------------------------------------
      if (do_rad) then

!--------------------------------------------------------------------
!    accumulate global integral quantities 
!--------------------------------------------------------------------
        call sum_diag_integral_field ('olr',    olr,        is, js)
        call sum_diag_integral_field ('abs_sw', swin-swout, is, js)

!--------------------------------------------------------------------
!    accumulate hemispheric integral quantities, if desired. 
!--------------------------------------------------------------------
        if (calc_hemi_integrals) then
          do j=js,je        
            jind = j - js + 1
            iind = 1  ! are assuming all i points are at same latitude

!---------------------------------------------------------------------
!    calculate southern hemisphere integrals.
!---------------------------------------------------------------------
            if (lat(iind,jind) <= 0.0) then
              call sum_diag_integral_field ('sntop_tot_sh ',   &
                                            swin-swout, is, ie, j, j)
              call sum_diag_integral_field ('lwtop_tot_sh ', olr,     &
                                            is, ie, j, j)
              call sum_diag_integral_field ('sngrd_tot_sh ',   &
                                            swdns-swups, is, ie, j, j)
              call sum_diag_integral_field ('lwgrd_tot_sh ',   &
                                          Lw_output%flxnet(:,:,kmax+1),&
                                          is, ie,  j, j)
              if (do_clear_sky_pass) then
                call sum_diag_integral_field ('sntop_clr_sh ',   &
                                              swin_clr-swout_clr, &
                                              is, ie, j, j)
                call sum_diag_integral_field ('lwtop_clr_sh ', olr_clr,&
                                              is, ie, j, j)
                call sum_diag_integral_field ('sngrd_clr_sh ',   &
                                              swdns_clr-swups_clr, &
                                              is, ie, j, j)
                call sum_diag_integral_field ('lwgrd_clr_sh ',    &
                                       Lw_output%flxnetcf(:,:,kmax+1),&
                                       is, ie, j, j)
              endif

!---------------------------------------------------------------------
!    calculate northern hemisphere integrals.
!---------------------------------------------------------------------
            else
              call sum_diag_integral_field ('sntop_tot_nh ',    &
                                            swin-swout, is, ie, j, j)
              call sum_diag_integral_field ('lwtop_tot_nh ', olr,     &
                                            is, ie, j, j)
              call sum_diag_integral_field ('sngrd_tot_nh ',   &
                                            swdns-swups, is, ie, j, j)
              call sum_diag_integral_field ('lwgrd_tot_nh ',   &
                                          Lw_output%flxnet(:,:,kmax+1),&
                                          is, ie, j, j)        
              if (do_clear_sky_pass) then
                call sum_diag_integral_field ('sntop_clr_nh ',   &
                                              swin_clr-swout_clr,  &
                                              is, ie, j, j)
                call sum_diag_integral_field ('lwtop_clr_nh ', olr_clr,&
                                              is, ie, j, j)
                call sum_diag_integral_field ('sngrd_clr_nh ',   &
                                              swdns_clr-swups_clr, &
                                              is, ie, j, j)
                call sum_diag_integral_field ('lwgrd_clr_nh ',   &
                                       Lw_output%flxnetcf(:,:,kmax+1),&
                                       is, ie, j, j)
              endif
            endif
          end do

!--------------------------------------------------------------------
!    accumulate global integral quantities 
!--------------------------------------------------------------------
          call sum_diag_integral_field ('sntop_tot_gl ', swin-swout,  &
                                        is, js)
          call sum_diag_integral_field ('lwtop_tot_gl ', olr, is, js)
          call sum_diag_integral_field ('sngrd_tot_gl ', swdns-swups, &
                                        is, js)
          call sum_diag_integral_field ('lwgrd_tot_gl ',  &
                                  Lw_output%flxnet(:,:,kmax+1), is, js)
          if (do_clear_sky_pass) then
            call sum_diag_integral_field ('sntop_clr_gl ',   &
                                          swin_clr-swout_clr, is, js)
            call sum_diag_integral_field ('lwtop_clr_gl ', olr_clr,   &
                                          is, js)
            call sum_diag_integral_field ('sngrd_clr_gl ',   &
                                          swdns_clr-swups_clr, is, js)
            call sum_diag_integral_field ('lwgrd_clr_gl ',   &
                                        Lw_output%flxnetcf(:,:,kmax+1),&
                                        is, js)
          endif
        endif   ! (calc_hemi_integrals)
      endif  ! (do_rad)

!---------------------------------------------------------------------



end subroutine radiation_netcdf



!###################################################################

subroutine deallocate_arrays ( Rad_gases_mdl, Cldrad_props_mdl, &
                               Astro_mdl, Astro2, Lw_output, &
			       Atmos_input_mdl, Fsrad_output, &
			       Sw_output, Cld_diagnostics)

!---------------------------------------------------------------------
!    deallocate_arrays deallocates the array space of local 
!    derived-type variables.
!---------------------------------------------------------------------

type(radiative_gases_type)  , intent(in)   :: Rad_gases_mdl
type(cldrad_properties_type), intent(in)   :: Cldrad_props_mdl
type(astronomy_type)        , intent(in)   :: Astro_mdl, Astro2
type(lw_output_type)        , intent(in)   :: Lw_output
type(atmos_input_type)      , intent(in)   :: Atmos_input_mdl
type(fsrad_output_type)     , intent(in)   :: Fsrad_output
type(sw_output_type)        , intent(in)   :: Sw_output
type(cld_diagnostics_type)  , intent(in)   :: Cld_diagnostics 


!--------------------------------------------------------------------
!    deallocate the variables in Rad_gases_mdl.
!--------------------------------------------------------------------
      if (do_rad .or. do_average_gases) then
        deallocate (Rad_gases_mdl%qo3)
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Astro_mdl and Astro2.
!--------------------------------------------------------------------
      if ( do_rad .or. do_average .or. renormalize_sw_fluxes ) then 
        deallocate (Astro_mdl%solar)
        deallocate (Astro_mdl%cosz )
        deallocate (Astro_mdl%fracday)
        if (Environment%running_gcm) then
          if ( do_rad .and. renormalize_sw_fluxes   &
	              .and. do_diurnal ) then 
            deallocate (Astro2%solar)
            deallocate (Astro2%cosz )
            deallocate (Astro2%fracday)
          endif
        endif
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Lw_output.
!--------------------------------------------------------------------
      if (do_sea_esf_rad) then
        if (do_rad) then
          deallocate (Lw_output%heatra    )
          deallocate (Lw_output%flxnet    )
          if (Rad_control%do_totcld_forcing) then
            deallocate (Lw_output%heatracf  )
            deallocate (Lw_output%flxnetcf  )
          endif

!--------------------------------------------------------------------
!    deallocate the variables in Sw_output.
!--------------------------------------------------------------------
          deallocate (Sw_output%dfsw     )
          deallocate (Sw_output%ufsw     )
          deallocate (Sw_output%fsw     )
          deallocate (Sw_output%hsw     )
          if (Rad_control%do_totcld_forcing) then
            deallocate (Sw_output%dfswcf   )
            deallocate (Sw_output%ufswcf   )
            deallocate (Sw_output%fswcf   )
            deallocate (Sw_output%hswcf   )
          endif
        endif
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Atmos_input_mdl.
!--------------------------------------------------------------------
      if ( do_rad .or. do_average ) then 
        deallocate (Atmos_input_mdl%press    )
        deallocate (Atmos_input_mdl%temp     )
        deallocate (Atmos_input_mdl%rh2o     )
        deallocate (Atmos_input_mdl%rel_hum  )
        deallocate (Atmos_input_mdl%pflux    )
        deallocate (Atmos_input_mdl%tflux    )
        deallocate (Atmos_input_mdl%deltaz   )
        deallocate (Atmos_input_mdl%psfc     )
        deallocate (Atmos_input_mdl%tsfc     )
        deallocate (Atmos_input_mdl%asfc     )
        deallocate (Atmos_input_mdl%land     )
        deallocate (Atmos_input_mdl%cloud_water)
        deallocate (Atmos_input_mdl%cloud_ice)
        deallocate (Atmos_input_mdl%cloudtemp)
        deallocate (Atmos_input_mdl%cloudvapor )
        deallocate (Atmos_input_mdl%clouddeltaz)
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Cldrad_props_mdl. different variables
!    exist dependent on the sw parameterization being used.
!--------------------------------------------------------------------
      if ((do_rad .or. do_average_clouds) .and. do_sea_esf_rad) then
        deallocate (Cldrad_props_mdl%camtsw    )
        deallocate (Cldrad_props_mdl%cmxolw    )
        deallocate (Cldrad_props_mdl%crndlw    )
        deallocate (Cldrad_props_mdl%ncldsw    )
        deallocate (Cldrad_props_mdl%nmxolw     )
        deallocate (Cldrad_props_mdl%nrndlw    )

        deallocate (Cldrad_props_mdl%emmxolw   )
        deallocate (Cldrad_props_mdl%emrndlw   )

        if ( associated(Cldrad_props_mdl%cldext) ) then
          deallocate (Cldrad_props_mdl%cldext    )
          deallocate (Cldrad_props_mdl%cldasymm  )
          deallocate (Cldrad_props_mdl%cldsct    )
        endif
        if ( associated(Cldrad_props_mdl%cldemiss) ) then
          deallocate (Cldrad_props_mdl%abscoeff  )
          deallocate (Cldrad_props_mdl%cldemiss  )
        endif
        if ( associated(Cldrad_props_mdl%cvisrfsw) ) then
          deallocate (Cldrad_props_mdl%cvisrfsw  )
          deallocate (Cldrad_props_mdl%cirabsw  )
          deallocate (Cldrad_props_mdl%cirrfsw  )
        endif

!--------------------------------------------------------------------
!    deallocate the variables in Cld_diagnostics.
!--------------------------------------------------------------------
        if (associated (Cld_diagnostics%lwpath) ) then
          deallocate (Cld_diagnostics%lwpath    )
          deallocate (Cld_diagnostics%iwpath    )
          deallocate (Cld_diagnostics%size_drop )
          deallocate (Cld_diagnostics%size_ice  )
        endif
        if (associated (Cld_diagnostics%cld_isccp_hi)) then
          deallocate (Cld_diagnostics%cld_isccp_hi)
          deallocate (Cld_diagnostics%cld_isccp_mid)
          deallocate (Cld_diagnostics%cld_isccp_low)
          deallocate (Cld_diagnostics%tot_clds  )
        endif
      endif

!--------------------------------------------------------------------
!    deallocate the variables in Fsrad_output.
!--------------------------------------------------------------------

      if (.not. do_sea_esf_rad .and. do_rad) then
         deallocate (Fsrad_output%tdtsw   )
         deallocate (Fsrad_output%tdtlw   )
         deallocate (Fsrad_output%swdns   )
         deallocate (Fsrad_output%swups   )
         deallocate (Fsrad_output%lwdns   )
         deallocate (Fsrad_output%lwups   )
         deallocate (Fsrad_output%swin    )
         deallocate (Fsrad_output%swout   )
         deallocate (Fsrad_output%olr     )
        if (do_clear_sky_pass) then
           deallocate (Fsrad_output%tdtsw_clr )
           deallocate (Fsrad_output%tdtlw_clr )
           deallocate (Fsrad_output%swdns_clr )
           deallocate (Fsrad_output%swups_clr )
           deallocate (Fsrad_output%lwdns_clr )
           deallocate (Fsrad_output%lwups_clr )
           deallocate (Fsrad_output%swin_clr  )
           deallocate (Fsrad_output%swout_clr )
           deallocate (Fsrad_output%olr_clr   )
        endif 
      endif 

!---------------------------------------------------------------------



end subroutine deallocate_arrays 



!#####################################################################




                 end module radiation_driver_mod




