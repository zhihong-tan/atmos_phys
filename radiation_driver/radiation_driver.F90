                module radiation_driver_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="">
!  
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    radiation_driver_mod is the interface between physics_driver_mod
!    and a specific radiation parameterization, currently either the
!    original_fms_rad or sea_esf_rad radiation package. it provides 
!    radiative heating rates, boundary radiative fluxes, and any other 
!    radiation package output fields to other component models of the
!    modeling system.
! </OVERVIEW>
! <DESCRIPTION>
!    radiation_driver_mod is the interface between physics_driver_mod
!    and a specific radiation parameterization, currently either the
!    original_fms_rad or sea_esf_rad radiation package. it provides 
!    radiative heating rates, boundary radiative fluxes, and any other 
!    radiation package output fields to other component models of the
!    modeling system.
! </DESCRIPTION>

!   shared modules:

use fms_mod,               only: fms_init, mpp_clock_id, &
                                 mpp_clock_begin, mpp_clock_end, &
                                 CLOCK_MODULE, MPP_CLOCK_SYNC, &
                                 mpp_pe, mpp_root_pe, &
                                 open_namelist_file, stdlog, &
                                 file_exist, FATAL, WARNING, NOTE, &
                                 close_file, read_data, write_data, &
                                 write_version_number, check_nml_error,&
                                 error_mesg, open_restart_file
use diag_manager_mod,      only: register_diag_field, send_data, &
                                 diag_manager_init
use time_manager_mod,      only: time_type, set_date, set_time,  &
                                 get_time,    operator(+),       &
                                 time_manager_init, &
                                 operator(-), operator(/=), get_date,&
                                 operator(<), operator(>=), operator(>)
use sat_vapor_pres_mod,    only: sat_vapor_pres_init, lookup_es
use constants_mod,         only: constants_init, RDGAS, RVGAS,   &
                                 STEFAN, GRAV, SECONDS_PER_DAY

! shared radiation package modules:

use rad_utilities_mod,     only: radiation_control_type, Rad_control, &
                                 radiative_gases_type, &
                                 check_derived_types, &
                                 cldrad_properties_type, &
                                 astronomy_type, surface_type, &
                                 cld_specification_type, &
                                 atmos_input_type, rad_utilities_init,&
                                 aerosol_properties_type, aerosol_type,&
                                 sw_output_type, lw_output_type, &
                                 rad_output_type, microphysics_type, &
                                 shortwave_control_type, Sw_control, &
                                 fsrad_output_type, &
                                 environment_type, Environment, &
                                 cloudrad_control_type, Cldrad_control,&
                                 rad_utilities_end

!  physics support modules:

use diag_integral_mod,     only: diag_integral_init, &
                                 diag_integral_field_init, &
                                 sum_diag_integral_field
use astronomy_mod,         only: astronomy_init, annual_mean_solar, &
                                 daily_mean_solar, diurnal_solar, &
                                 astronomy_end

!  component modules:

use original_fms_rad_mod,  only: original_fms_rad_init,  &
                                 original_fms_rad, &
                                 original_fms_rad_end
use sea_esf_rad_mod,       only: sea_esf_rad_init, sea_esf_rad, &
                                 sea_esf_rad_end
use rad_output_file_mod,   only: rad_output_file_init, &
                                 write_rad_output_file,    &
                                 rad_output_file_end
use cloudrad_package_mod,  only: cloudrad_package_init, &
                                 cloud_radiative_properties, &
                                 cldrad_props_dealloc, &
                                 cloudrad_package_end
use aerosolrad_package_mod, only: aerosolrad_package_init,    &
                                  aerosol_radiative_properties, &
                                  aerosolrad_package_end

!--------------------------------------------------------------------

implicit none 
private 

!----------------------------------------------------------------------
!    radiation_driver_mod is the interface between physics_driver_mod
!    and a specific radiation parameterization, currently either the
!    original_fms_rad or sea_esf_rad radiation package. it provides 
!    radiative heating rates, boundary radiative fluxes, and any other 
!    radiation package output fields to other component models of the
!    modeling system.
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!------------ version number for this module --------------------------

character(len=128) :: version = '$Id: radiation_driver.F90,v 10.0 2003/10/24 22:00:37 fms Exp $'
character(len=128) :: tagname = '$Name: jakarta $'


!---------------------------------------------------------------------
!------ interfaces -----

public    radiation_driver_init, radiation_driver,   &
          define_rad_times, define_atmos_input_fields,  &
          define_surface, surface_dealloc, atmos_input_dealloc, &
          radiation_driver_end

private  & 

! called from radiation_driver_init:
          read_restart_file, initialize_diagnostic_integrals,   &
          diag_field_init, &

! called from radiation_driver:
          obtain_astronomy_variables, radiation_calc,    &
          update_rad_fields, produce_radiation_diagnostics,  &
          deallocate_arrays, &

! called from define_atmos_input_fields:
          calculate_auxiliary_variables


!-----------------------------------------------------------------------
!------- namelist ---------

integer ::  rad_time_step = 0         !  radiative time step in seconds
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


namelist /radiation_driver_nml/ rad_time_step, do_clear_sky_pass, &
                                zenith_spec, rad_package,    &
                                calc_hemi_integrals,     &
                                all_step_diagnostics, &
                                renormalize_sw_fluxes, &
                                rad_date, all_level_radiation, &
                                topmost_radiation_level,   &
                                drop_upper_levels,  &
                                all_column_radiation, rsd,    &
                                use_mixing_ratio, solar_constant

!---------------------------------------------------------------------
!---- public data ----


!---------------------------------------------------------------------
!---- private data ----


!---------------------------------------------------------------------
!    logical  flags.

logical ::  module_is_initialized = .false. ! module initialized?
logical ::  do_rad                          ! is this a radiation step ?
logical ::  use_rad_date                    ! specify time of radiation
                                            ! independent of model time?
logical ::  do_sea_esf_rad                  ! using sea_esf_rad package?

!---------------------------------------------------------------------
!    list of restart versions of radiation_driver.res readable by this 
!    module, dependent on which radiation package is activated. restart 
!    version 1 of radiation_driver.res is not readable by this module; 
!    however, restart version 1 of sea_esf_rad.res IS readable within 
!    this code.
!     version 1:  sea_esf_rad.res file version used initially in 
!                 AM2 model series  
!     version 2:  added cosine of zenith angle as an output to
!                 radiation_driver.res  (6/27/00)
!     version 3:  added restart variables needed when sw renormalization
!                 is active. (3/21/02)
!     version 4:  added longwave heating rate as separate output 
!                 variable, since it is needed as input to edt_mod
!                 and entrain_mod. (7/17/02)
!     version 5:  removed variables associated with the former 
!                 do_average namelist option (7/23/03)
!---------------------------------------------------------------------
integer, dimension(4) :: restart_versions     = (/ 2, 3, 4, 5 /)
integer, dimension(4) :: restart_versions_sea = (/ 2, 3, 4, 5 /)


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

type(rad_output_type),save          ::  Rad_output
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
integer    :: num_pts=0      !  counter for current number of grid 
                             !  columns processed (when num_pts=0 or 
                             !  num_pts=total_pts certain things happen)
integer    :: total_pts      !  number of grid columns to be processed 
                             !  every time step (note: all grid columns
                             !  must be processed every time step)
type(time_type) :: Rad_time  !  time at which the climatologically-
                             !  determined, time-varying input fields to
                             !  radiation should apply 
                             !  [ time_type (days, seconds)]
integer    :: dt             !  physics time step (frequency of calling 
                             !  radiation_driver)  [ seconds ]


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

integer                      :: misc_clock, clouds_clock, calc_clock

!--------------------------------------------------------------------
! miscellaneous variables and indices

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
real,parameter ::  D608 = (RVGAS-RDGAS)/RDGAS
                              !  virtual temperature factor  
real,parameter ::  D622 = RDGAS/RVGAS
                              ! ratio of gas constants - dry air to 
                              ! water vapor
real,parameter ::  D378 = 1.0 - D622  
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
! <SUBROUTINE NAME="radiation_driver_init">
!  <OVERVIEW>
!   radiation_driver_init is the constructor for radiation_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   radiation_driver_init is the constructor for radiation_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_driver_init (lonb, latb, pref, axes, Time, &
!                                  aerosol_names)
!
!  </TEMPLATE>
!  <IN NAME="lonb" TYPE="real">
!    lonb      Longitude in radians for all (i.e., the global size)
!              grid box boundaries, the size of lonb should be one more
!              than the global number of longitude points along the x-axis.
!                 [real, dimension(:)]
!  </IN>
!  <IN NAME="latb" TYPE="real">
!    latb      Latitude in radians for all (i.e., the global size)
!              grid box boundaries, the size of latb should be one more
!              than the global number of latitude points along the y-axis.
!                 [real, dimension(:)]
!  </IN>
!  <IN NAME="pref" TYPE="real">
!    pref      Two reference profiles of pressure at full model levels
!              plus the surface (nlev+1). The first profile assumes a surface
!              pressure of 101325 pa, and the second profile assumes 
!              81060 pa.  [real, dimension(nlev+1,2)]
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!    axes      The axis indices that are returned by previous calls to
!              diag_axis_init. The values of this array correspond to the
!              x, y, full (p)level, and half (p)level axes. These are the
!              axes that diagnostic fields are output on.
!                 [integer, dimension(4)]
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!    Time      The current time.  [time_type]
!  </IN>
!  <IN NAME="aerosol_names" TYPE="character">
!   Aerosol names
!  </IN>
! </SUBROUTINE>
!
subroutine radiation_driver_init (lonb, latb, pref, axes, Time, &
                                  aerosol_names)

!---------------------------------------------------------------------
!   radiation_driver_init is the constructor for radiation_driver_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real, dimension(:),              intent(in)  :: lonb, latb
real, dimension(:,:),            intent(in)  :: pref
integer, dimension(4),           intent(in)  :: axes
type(time_type),                 intent(in)  :: Time
character(len=*), dimension(:), intent(in)   :: aerosol_names
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       lonb           array of model longitudes on cell boundaries 
!                      [ radians ]
!       latb           array of model latitudes at cell boundaries 
!                      [ radians ]
!       pref           array containing two reference pressure profiles 
!                      for use in defining transmission functions
!                      [ pascals ]
!       axes           diagnostic variable axes
!       Time           current time [time_type(days, seconds)]
!       aerosol_names  names associated with the activated aerosol
!                      species
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables

      integer           ::   unit, io, ierr
      integer           ::   id, jd, kmax 

!---------------------------------------------------------------------
!   local variables
! 
!        unit    io unit number for namelist file
!        io      error status returned from io operation
!        ierr    error code
!        id      number of grid points in x direction (on processor)
!        jd      number of grid points in y direction (on processor)
!        kmax    number of model layers
!                
!---------------------------------------------------------------------

      
!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call rad_utilities_init
      call diag_manager_init
      call time_manager_init
      call sat_vapor_pres_init
      call constants_init
      call diag_integral_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=radiation_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'radiation_driver_nml')
        enddo
10      call close_file (unit)
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
           write (stdlog(), nml=radiation_driver_nml)

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
        Sw_control%do_diurnal = .true.
        Sw_control%do_annual = .false.
        Sw_control%do_daily_mean = .false.
      else if (zenith_spec == 'daily_mean') then
        Sw_control%do_diurnal = .false.
        Sw_control%do_annual = .false.
        Sw_control%do_daily_mean = .true.
      else if (zenith_spec == 'annual_mean') then
        Sw_control%do_diurnal = .false.
        Sw_control%do_annual = .true.
        Sw_control%do_daily_mean = .false.
      else
        call error_mesg ('radiation_driver_mod', &    
            'string provided for zenith_spec is invalid', FATAL)
      endif

      Sw_control%solar_constant = solar_constant

!---------------------------------------------------------------------
!    set flags indicating that the Sw_control variables have been 
!    defined.
!---------------------------------------------------------------------
      Sw_control%do_diurnal_iz = .true.
      Sw_control%do_annual_iz = .true.
      Sw_control%do_daily_mean_iz = .true.

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
        ' must specify all_level_radiation and all_column_radiation'//&
            ' as true when using original fms radiation', FATAL)
        endif
        if (renormalize_sw_fluxes) then
          call error_mesg ( 'radiation_driver_mod', &
           ' cannot renormalize shortwave fluxes with original_fms '//&
                 'radiation package.', FATAL)
        endif
        if (all_step_diagnostics) then
          call error_mesg ( 'radiation_driver_mod', &
            ' cannot request all_step_diagnostics with original_fms '//&
              'radiation package.', FATAL)
        endif
      endif

!---------------------------------------------------------------------
!    can only renormalize shortwave fluxes when diurnally_varying
!    radiation is used.
!---------------------------------------------------------------------
     if (renormalize_sw_fluxes .and. .not. Sw_control%do_diurnal) then
       call error_mesg ('radiation_driver_mod',  &
       ' can only renormalize sw fluxes when using diurnally-varying'//&
                       ' solar radiation', FATAL)
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
!    define the dimensions of the local processors portion of the grid.
!---------------------------------------------------------------------
      id    = size(lonb,1) - 1 
      jd    = size(latb,1) - 1
      kmax  = size(pref,1) - 1 

!---------------------------------------------------------------------
!    check for consistency if drop_upper_levels is activated.
!----------------------------------------------------------------------
      if (drop_upper_levels .and. all_level_radiation) then
          call error_mesg ( 'radiation_driver_mod',  &
            ' drop_upper_levels and all_level_radiation are '//&
                                         'incompatible', FATAL)
      endif

!---------------------------------------------------------------------
!    define the starting and ending vertical indices of the radiation
!    grid. if all_level_radiation is .true., then radiation is done
!    at all model levels. ks, ke are model-based coordinates, while
!    ksrad and kerad are radiation-grid based coordinates (ksrad always
!    is equal to 1). 
!---------------------------------------------------------------------
      if (all_level_radiation) then
        ks = 1
        ke = kmax
        kerad = kmax
        topmost_radiation_level = 1
      else
        if (topmost_radiation_level <= 0) then
          call error_mesg ('radiation_driver_mod', &
          ' when all_level_radiation is .false., topmost_radiation'//&
              '_level must be specified as a positive integer.', FATAL)
        endif
        if (drop_upper_levels) then
          ks = topmost_radiation_level
          ke = kmax
          kerad = ke - ks + 1
          call error_mesg ( ' radiation_driver_mod', &
            ' code has not been validated for all_level_radiation = '//&
               'false. DO NOT USE!', FATAL)
        else
          call error_mesg ( ' radiation_driver_mod', &
         ' currently only drop_upper_levels is available as option '//&
                           'when all_level_radiation = false.', FATAL)
        endif
      endif
       
!---------------------------------------------------------------------
!    exit if all_column_radiation is not .true. -- this option is not
!    yet certified.
!---------------------------------------------------------------------
      if (.not. all_column_radiation) then
        call error_mesg ('radiation_driver_mod',  &
          ' code currently not validated for all_column_radiation = '//&
                                  'false. DO NOT USE!', FATAL)
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
      if (Environment%running_standalone .and.  &
          renormalize_sw_fluxes) then
        call error_mesg ( 'radiation_driver_mod', &
               'renormalizing sw fluxes when running standalone'//&
                           ' radiation code is not available', FATAL)
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
!    if two radiation restart files exist, exit.
!-----------------------------------------------------------------------
        if ( file_exist('INPUT/sea_esf_rad.res')  .and.     &
             file_exist('INPUT/radiation_driver.res') ) then 
          call error_mesg ('radiation_driver_mod',  &
         ' both sea_esf_rad.res and radiation_driver.res files are'//&
               ' present in INPUT directory. which one to use ?', FATAL)
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
!----------------------------------------------------------------------
!    if no restart file is present, initialize the needed fields until
!    the radiation package may be called. initial surface flux is set 
!    to 100 wm-2, and is only used for initial guess of sea ice temp.
!    set rad_alarm to be 1 second from now, ie., on the first step of 
!    the job.
!-----------------------------------------------------------------------
        else
          rad_alarm                = 1
          Rad_output%tdt_rad       = 0.0
          Rad_output%tdt_rad_clr   = 0.0
          Rad_output%tdtlw         = 0.0
          Rad_output%flux_sw_surf  = surf_flx_init
          Rad_output%flux_lw_surf  = surf_flx_init
          Rad_output%coszen_angle  = coszen_angle_init
          if (mpp_pe() == mpp_root_pe() ) then
            call error_mesg ('radiation_driver_mod', &
           'no acceptable radiation restart file present; therefore'//&
           ' will initialize input fields', NOTE)
          endif
        endif
      endif   ! (running_gcm)

!--------------------------------------------------------------------
!    do the initialization specific to the sea_esf_rad radiation
!    package.
!--------------------------------------------------------------------
      if (do_sea_esf_rad) then 

!---------------------------------------------------------------------
!    define control variables indicating whether the clear-sky forcing
!    should be calculated. set a flag to indicate that the variable
!    has been defined.
!---------------------------------------------------------------------
        Rad_control%do_totcld_forcing = do_clear_sky_pass
        Rad_control%do_totcld_forcing_iz = .true.

!---------------------------------------------------------------------
!    initialize the modules that are accessed from radiation_driver_mod.
!---------------------------------------------------------------------
        call sea_esf_rad_init        (lonb, latb, pref(ks:ke+1,:))
        call cloudrad_package_init   (pref(ks:ke+1,:), lonb, latb,  &
                                      axes, Time)
        call aerosolrad_package_init (aerosol_names)
        call rad_output_file_init    (axes, Time, aerosol_names)

!---------------------------------------------------------------------
!    do the initialization specific to the original fms radiation
!    package. 
!---------------------------------------------------------------------
      else
        call original_fms_rad_init (lonb, latb, pref, axes, Time, kmax)
      endif

!--------------------------------------------------------------------
!    initialize the astronomy_package.
!--------------------------------------------------------------------
      if (Sw_control%do_annual) then
        call astronomy_init (latb, lonb)
      else
        call astronomy_init
      endif

      if (Environment%running_gcm .or.  &
          Environment%running_sa_model) then
!---------------------------------------------------------------------
!    initialize the total number of columns in the processor's domain.
!---------------------------------------------------------------------
        total_pts = id*jd 

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
      calc_clock =    &
            mpp_clock_id ('   Physics_down: Radiation: calc', &
                grain = CLOCK_MODULE, flags = MPP_CLOCK_SYNC)

!---------------------------------------------------------------------
!    call check_derived_types to verify that all logical elements of
!    public derived-type variables stored in rad_utilities_mod but
!    initialized elsewhere have been initialized.
!---------------------------------------------------------------------
      if (do_sea_esf_rad) then
        call check_derived_types
      endif

!---------------------------------------------------------------------
!    set flag to indicate that module has been successfully initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------


end subroutine radiation_driver_init



!#####################################################################
! <SUBROUTINE NAME="radiation_driver">
!  <OVERVIEW>
!    radiation_driver adds the radiative heating rate to the temperature
!    tendency and obtains the radiative boundary fluxes and cosine of 
!    the solar zenith angle to be used in the other component models.
!  </OVERVIEW>
!  <DESCRIPTION>
!    radiation_driver adds the radiative heating rate to the temperature
!    tendency and obtains the radiative boundary fluxes and cosine of 
!    the solar zenith angle to be used in the other component models.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_driver (is, ie, js, je, Time, Time_next,  &
!                             lat, lon, Surface, Atmos_input, &
!                             Aerosol, Cld_spec, Rad_gases, &
!                             Lsc_microphys, Meso_microphys,    &
!                             Cell_microphys, Radiation, mask, kbot)
!
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 
!  <IN NAME="Time" TYPE="time_type">
!   current model time 
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   The time used for diagnostic output
!  </IN>
!  <IN NAME="lon" TYPE="real">
!    lon        mean longitude (in radians) of all grid boxes processed by
!               this call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <INOUT NAME="Surface" TYPE="surface_type">
!   Surface input data to radiation package
!  </INOUT>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data to radiation package
!  </INOUT>
!  <INOUT NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatological input data to radiation package
!  </INOUT>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud microphysical and physical parameters to radiation package, 
!                     contains var-
!                     iables defining the cloud distribution, passed 
!                     through to lower level routines
!  </INOUT>
!  <INOUT NAME="Rad_gases" TYPE="radiative_gases_type">
!   Radiative gases properties to radiation package, , contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale
!                      clouds
!  </INOUT>
!  <INOUT NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                      clouds associated with donner convection
!  </INOUT>
!  <INOUT NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale
!                      clouds assciated with donner convection
!  </INOUT>
!  <INOUT NAME="Radiation" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </INOUT>
!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
subroutine radiation_driver (is, ie, js, je, Time, Time_next,  &
                             lat, lon, Surface, Atmos_input, &
                             Aerosol, Cld_spec, Rad_gases, &
                             Lsc_microphys, Meso_microphys,    &
                             Cell_microphys, Radiation, mask, kbot)

!---------------------------------------------------------------------
!    radiation_driver adds the radiative heating rate to the temperature
!    tendency and obtains the radiative boundary fluxes and cosine of 
!    the solar zenith angle to be used in the other component models.
!---------------------------------------------------------------------
 
!--------------------------------------------------------------------
integer,                      intent(in)           :: is, ie, js, je
type(time_type),              intent(in)           :: Time, Time_next
real, dimension(:,:),         intent(in)           :: lat, lon
type(surface_type),           intent(inout)        :: Surface
type(atmos_input_type),       intent(inout)        :: Atmos_input
type(aerosol_type),           intent(inout)        :: Aerosol  
type(cld_specification_type), intent(inout)        :: Cld_spec
type(radiative_gases_type),   intent(inout)        :: Rad_gases
type(microphysics_type),      intent(inout)        :: Lsc_microphys,&
                                                      Meso_microphys,&
                                                      Cell_microphys
type(rad_output_type),     intent(inout), optional :: Radiation
real, dimension(:,:,:),    intent(in),    optional :: mask
integer, dimension(:,:),   intent(in),    optional :: kbot
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Time           current model time [ time_type (days, seconds) ] 
!      Time_next      time on next timestep, used as stamp for diagnos-
!                     tic output  [ time_type  (days, seconds) ]  
!      lat            latitude of model points  [ radians ]
!      lon            longitude of model points [ radians ]
!
!   intent(inout) variables:
!
!      Surface        surface_type structure, contains variables 
!                     defining the surface characteristics, including
!                     the following component referenced in this 
!                     routine:
!
!         asfc          surface albedo  [ dimensionless ]
!
!      Atmos_input    atmos_input_type structure, contains variables
!                     defining atmospheric state, including the follow-
!                     ing component referenced in this routine
!
!         tsfc          surface temperature [ deg K ]
!
!      Aerosol        aerosol_type structure, contains variables
!                     defining aerosol fields, passed through to
!                     lower level routines
!      Cld_spec       cld_specification_type structure, contains var-
!                     iables defining the cloud distribution, passed 
!                     through to lower level routines
!      Rad_gases      radiative_gases_type structure, contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!      Lsc_microphys  microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     large-scale clouds, passed through to lower
!                     level routines
!      Meso_microphys microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     meso-scale clouds, passed through to lower
!                     level routines
!      Cell_microphys microphysics_type structure, contains variables
!                     describing the microphysical properties of the
!                     convective cell-scale clouds, passed through to 
!                     lower level routines
!
!   intent(inout), optional variables:
!
!      Radiation      rad_output_type structure, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.  present when running gcm,
!                     not present when running sa_gcm or standalone
!                     columns mode. variables defined here are:
!
!        tdt_rad         radiative (sw + lw) heating rate
!                        [ deg K / sec ]
!        flux_sw_surf    net (down-up) sw surface flux 
!                        [ watts / m^^2 ]
!        flux_lw_surf    downward lw surface flux 
!                        [ watts / m^^2 ]
!        coszen_angle    cosine of the zenith angle which will be used 
!                        for the next ocean_albedo calculation 
!                        [ dimensionless ]
!        tdtlw           longwave heating rate
!                        [ deg K / sec ]
!
!   intent(in), optional variables:
!
!        mask            present when running eta vertical coordinate,
!                        mask to remove points below ground
!        kbot            present when running eta vertical coordinate,
!                        index of lowest model level above ground 
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      type(cldrad_properties_type)       :: Cldrad_props
      type(astronomy_type)               :: Astro, Astro2
      type(lw_output_type)               :: Lw_output
      type(sw_output_type)               :: Sw_output
      type(fsrad_output_type)            :: Fsrad_output
      type(aerosol_properties_type)      :: Aerosol_props

      real, dimension (ie-is+1, je-js+1) :: flux_ratio
 
!-------------------------------------------------------------------
!   local variables:
!
!      Cldrad_props      cloud radiative properties on model grid,
!                        [cldrad_properties_type]
!      Astro             astronomical properties on model grid, usually
!                        valid over radiation timestep
!                        [astronomy_type]
!      Astro2            astronomical properties on model grid, valid 
!                        over current physics timestep
!                        [astronomy_type]
!      Lw_output         sea longwave output fields on model grid,
!                        [lw_output_type]
!      Sw_output         esf shortwave output fields on model grid,
!                        [sw_output_type]
!      Fsrad_output      original fms radiation output fields on model
!                        grid, [fsrad_output_type]
!      flux_ratio        value  used to renormalize sw fluxes and 
!                        heating rates to account for earth-sun motion
!                        during the radiation timestep
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)
     
!---------------------------------------------------------------------
!    if this is a radiation step, or if the astronomical inputs to
!    radiation (solar, cosz, fracday, rrsun) need to be obtained 
!    because of time averaging or renormalization, call 
!    obtain_astronomy_variables to do so.
!---------------------------------------------------------------------
      call mpp_clock_begin (misc_clock)
      if (do_rad .or. renormalize_sw_fluxes) then 
        call obtain_astronomy_variables (is, ie, js, je, lat, lon,  &
                                         Astro, Astro2)  
      endif

      if (do_rad) then
        if (Rad_control%do_aerosol) then
          call aerosol_radiative_properties (is, ie, js, je, &
                                             Aerosol, Aerosol_props)
        endif
      endif
      call mpp_clock_end (misc_clock)

!--------------------------------------------------------------------
!    when using the sea-esf radiation, call cloud_radiative_properties
!    to obtain the cloud-radiative properties needed for the radiation 
!    calculation. (these properties are obtained within radiation_calc
!    when executing the original fms radiation code). if these fields 
!    are to be time-averaged, this call is made on all steps; otherwise
!    just on radiation steps.
!--------------------------------------------------------------------
      call mpp_clock_begin (clouds_clock)
      if (do_rad) then
        if (do_sea_esf_rad) then
          if (present(kbot) ) then
            call cloud_radiative_properties (     &
                             is, ie, js, je, Time_next, Astro,  & 
                             Atmos_input, Cld_spec, Lsc_microphys,  &
                             Meso_microphys, Cell_microphys,    &
                             Cldrad_props, kbot=kbot, mask=mask)
          else    

            call cloud_radiative_properties (      &
                             is, ie, js, je, Time_next, Astro,  & 
                             Atmos_input, Cld_spec, Lsc_microphys,   &
                             Meso_microphys, Cell_microphys,    &
                             Cldrad_props)
          endif
        endif
      endif
      call mpp_clock_end (clouds_clock)

!---------------------------------------------------------------------
!    on radiation timesteps, call radiation_calc to determine new radia-
!    tive fluxes and heating rates.
!---------------------------------------------------------------------
      call mpp_clock_begin (calc_clock)
      if (do_rad) then
        call radiation_calc (is, ie, js, je, Rad_time, Time_next, lat, &
                             lon, Atmos_input, Surface, Rad_gases,  &
                             Aerosol_props, Aerosol, Cldrad_props, &
                             Cld_spec, Astro, Rad_output, Lw_output, &
                             Sw_output, Fsrad_output, mask=mask,   &
                             kbot=kbot)       
      endif

      call mpp_clock_end (calc_clock)
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
                                Sw_output, Astro, Rad_output,    &
                                flux_ratio)

!-------------------------------------------------------------------
!    call produce_radiation_diagnostics to produce radiation 
!    diagnostics, both fields and integrals.
!-------------------------------------------------------------------
        if (do_sea_esf_rad) then
          call produce_radiation_diagnostics        &
                            (is, ie, js, je, Time_next, lat, &
                             Atmos_input%tsfc, Surface%asfc,  &
                             flux_ratio,  Astro, Rad_output,  &
                             Lw_output=Lw_output,&
                             Sw_output=Sw_output)
        else
          call produce_radiation_diagnostics        &
                            (is, ie, js, je, Time_next, lat, &
                             Atmos_input%tsfc, Surface%asfc,  &
                             flux_ratio,  Astro, Rad_output,  &
                             Fsrad_output=Fsrad_output, mask=mask)
        endif 

!---------------------------------------------------------------------
!    call write_rad_output_file to produce a netcdf output file of 
!    radiation-package-relevant variables. note that this is called
!    only on radiation steps, so that the effects of sw renormalization
!    will not be seen in the variables of the data file written by
!    write_rad_output_file.
!---------------------------------------------------------------------
        if (do_rad .and. do_sea_esf_rad) then
          if (Rad_control%do_aerosol) then
            call write_rad_output_file (is, ie, js, je,  &
                                        Atmos_input, Surface,   &
                                        Rad_output, Sw_output,  &
                                        Lw_output, Rad_gases,   & 
                                        Cldrad_props, Cld_spec, & 
                                        Time_next, Aerosol%aerosol)
          else
            call write_rad_output_file (is, ie, js, je,  &
                                        Atmos_input,Surface, &
                                        Rad_output, Sw_output,   &
                                        Lw_output, Rad_gases,   &
                                        Cldrad_props, Cld_spec, &
                                        Time_next)
          endif
        endif ! (do_rad and do_sea_esf_rad)
      endif  ! (running_gcm)

!---------------------------------------------------------------------
!    call deallocate_arrays to deallocate the array space associated 
!    with stack-resident derived-type variables.
!---------------------------------------------------------------------
        call deallocate_arrays (Cldrad_props, Astro, Astro2,    &
                                Lw_output, Fsrad_output, Sw_output) 

!---------------------------------------------------------------------
!    complete radiation step when running within a gcm. update the 
!    points-processed counter. if all points in the processor's sub-
!    domain have been processed, set the counter to zero, set the 
!    radiation alarm to go off rad_time_step seconds from now, and
!    set do_rad to false, so that radiation will not be calculated until
!    the alarm goes off.
!--------------------------------------------------------------------
      if (Environment%running_gcm .or.  &
          Environment%running_sa_model) then
        num_pts = num_pts + size(lat,1) * size(lat,2)
        if (num_pts == total_pts) then
          num_pts = 0
          if (do_rad) then
            rad_alarm = rad_alarm + rad_time_step
            do_rad = .false.
          endif
        endif
      endif  ! (running_gcm)

!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!    define the elements of the rad_output_type variable which will
!    return the needed radiation package output to the calling routine.
!    Radiation is currently present when running within a gcm, but
!    not present for other applications.
!--------------------------------------------------------------------
      if (present (Radiation)) then 
        Radiation%coszen_angle(:,:) =      &
                                  Rad_output%coszen_angle(is:ie,js:je)
        Radiation%tdt_rad(:,:,:) =   &
                                  Rad_output%tdt_rad(is:ie,js:je,:)
        Radiation%flux_sw_surf(:,:) =    &
                                  Rad_output%flux_sw_surf(is:ie,js:je)
        Radiation%flux_lw_surf(:,:)    =   &
                                  Rad_output%flux_lw_surf(is:ie,js:je)
        Radiation%tdtlw(:,:,:)         =     &
                                  Rad_output%tdtlw(is:ie,js:je,:)   
      endif
      call mpp_clock_end (misc_clock)

!---------------------------------------------------------------------




end subroutine radiation_driver



!#####################################################################
! <SUBROUTINE NAME="define_rad_times">
!  <OVERVIEW>
!    subroutine define_rad_times determines whether radiation is to be 
!    calculated on the current timestep, and defines logical variables 
!    which determine whether various input fields to radiation_driver 
!    need to be retrieved on the current step.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine define_rad_times determines whether radiation is to be 
!    calculated on the current timestep, and defines logical variables 
!    which determine whether various input fields to radiation_driver 
!    need to be retrieved on the current step.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  define_rad_times (Time, Time_next, dt_in, Rad_time_out,    &
!                             need_aerosols, need_clouds, need_gases,  &
!                             need_basic)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current model time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   The time used for diagnostic output
!  </IN>
!  <IN NAME="dt_in" TYPE="real">
!   physics time step (frequency of calling radiaiton_driver)
!  </IN>
!  <INOUT NAME="Rad_time_out" TYPE="time_type">
!   time at which the climatologically-determined,
!                     time-varying input fields to radiation should 
!                     apply    
!  </INOUT>
!  <OUT NAME="need_aerosols" TYPE="logical">
!   aersosol input data is needed on this step ?
!  </OUT>
!  <OUT NAME="need_clouds" TYPE="logical">
!   cloud input data is needed on this step ?
!  </OUT>
!  <OUT NAME="need_gases" TYPE="logical">
!   radiative gas input data is needed on this step ?
!  </OUT>
!  <OUT NAME="need_basic" TYPE="logical">
!   atmospheric input fields are needed on this step ?
!  </OUT>
! </SUBROUTINE>
!
subroutine define_rad_times (Time, Time_next, dt_in, Rad_time_out,    &
                             need_aerosols, need_clouds, need_gases,  &
                             need_basic)

!--------------------------------------------------------------------
!    subroutine define_rad_times determines whether radiation is to be 
!    calculated on the current timestep, and defines logical variables 
!    which determine whether various input fields to radiation_driver 
!    need to be retrieved on the current step.
!-------------------------------------------------------------------- 

!---------------------------------------------------------------------
type(time_type), intent(in)     ::  Time, Time_next
real,            intent(in)     ::  dt_in
type(time_type), intent(inout)  ::  Rad_time_out
logical,         intent(out)    ::  need_aerosols, need_clouds,   &
                                    need_gases, need_basic
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     Time            current model time  
!                     [ time_type, days and seconds]
!     Time_next       model time on the next atmospheric timestep
!                     [ time_type, days and seconds]
!     dt_in           physics time step (frequency of calling 
!                     radiation_driver)
!                     [ seconds ]
!     
!   intent(inout) variables:
!
!     Rad_time_out    time at which the climatologically-determined, 
!                     time-varying input fields to radiation should 
!                     apply    
!                     [ time_type, days and seconds]
!
!   intent(out) variables:
!
!     need_aerosols   aersosol input data is needed on this step ?
!     need_clouds     cloud input data is needed on this step ?
!     need_gases      radiative gas input data is needed on this step ?
!     need_basic      atmospheric input fields are needed on this step ?
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

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    store the atmospheric timestep into a module variable for later
!    use.
!--------------------------------------------------------------------
      dt = dt_in

!--------------------------------------------------------------------
!    verify that the radiation timestep is an even multiple of the 
!    physics timestep.
!---------------------------------------------------------------------
      if (MOD(rad_time_step, dt) /= 0) then
        call error_mesg ('radiation_driver_mod',  &
       ' radiation timestep is not integral multiple of physics step', &
                                                           FATAL)
      endif

!-------------------------------------------------------------------
!    for the standalone case, new radiation outputs are calculated on 
!    every step, using climatological variable values at the time spec-
!    ified by the input argument Time. 
!-------------------------------------------------------------------
      if (Environment%running_standalone) then
        do_rad = .true.
        Rad_time = Time

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
      endif  ! (running_standalone)

!--------------------------------------------------------------------
!    set a logical variable indicating whether radiative gas input data
!    is needed on this step.
!--------------------------------------------------------------------
      if (do_rad) then
        need_gases = .true.
      else
        need_gases = .false.
      endif

!--------------------------------------------------------------------
!    set a logical variable indicating whether aerosol input data
!    is needed on this step.
!--------------------------------------------------------------------
      if (do_rad  .and. Rad_control%do_aerosol) then
        need_aerosols = .true.
      else
        need_aerosols = .false.
      endif

!--------------------------------------------------------------------
!    set a logical variable indicating whether cloud input data
!    is needed on this step.
!--------------------------------------------------------------------
      if (do_sea_esf_rad .and. do_rad) then
        need_clouds = .true.
      else
        need_clouds = .false.
      endif
      
!--------------------------------------------------------------------
!    set a logical variable indicating whether atmospheric input data
!    is needed on this step.
!--------------------------------------------------------------------
      if (need_clouds .or. need_aerosols .or. need_gases) then
        need_basic = .true.
      else
        need_basic = .false.
      endif
   
!---------------------------------------------------------------------
!    place the time at which radiation is to be applied into an output
!    variable.
!---------------------------------------------------------------------
      Rad_time_out = Rad_time

!---------------------------------------------------------------------



end subroutine define_rad_times


!######################################################################
! <SUBROUTINE NAME="define_atmos_input_fields">
!  <OVERVIEW>
!    define_atmos_input_fields converts the atmospheric input fields 
!    (pfull, phalf, t, q, ts) to the form needed by the radiation 
!    modules, and when needed returns radiation-ready fields of pressure
!    (press, psfc), temperature (temp, tsfc), water vapor mixing ratio 
!    (rh2o) and several auxiliary variables in the derived type 
!    structure Atmos_input. the optional input variables are present
!    when running radiative feedback studies (sa_model), and are needed
!    to allow variation of temperature and vapor fields while holding 
!    the aerosol and cloud amounts fixed.
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_atmos_input_fields converts the atmospheric input fields 
!    (pfull, phalf, t, q, ts) to the form needed by the radiation 
!    modules, and when needed returns radiation-ready fields of pressure
!    (press, psfc), temperature (temp, tsfc), water vapor mixing ratio 
!    (rh2o) and several auxiliary variables in the derived type 
!    structure Atmos_input. the optional input variables are present
!    when running radiative feedback studies (sa_model), and are needed
!    to allow variation of temperature and vapor fields while holding 
!    the aerosol and cloud amounts fixed.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  define_atmos_input_fields (is, ie, js, je, pfull, phalf, &
!                                      t, q, ts, Atmos_input, &
!                                      cloudtemp, cloudvapor, &
!                                      aerosoltemp, aerosolvapor, &
!                                      aerosolpress, kbot)  
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!    starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!   pressure at half levels
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at full levels
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity of water vapor at full levels
!  </IN>
!  <IN NAME="ts" TYPE="real">
!   surface temperature
!  </IN>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input type structure, contains the 
!                    following components defined in this subroutine
!  </INOUT>
!  <IN NAME="cloudtemp" TYPE="real">
!    temperature to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!  </IN>
!  <IN NAME="cloudvapor" TYPE="real">
!   water vapor to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!  </IN>
!  <IN NAME="aerosoltemp" TYPE="real">
!   required in sa_gcm mode, absent otherwise:
!                         temperature field to be used by aerosol param-
!                         eterization 
!  </IN>
!  <IN NAME="aerosolvapor" TYPE="real">
!   required in sa_gcm mode, absent otherwise:
!                         water vapor field to be used by aerosol param-
!                         eterization 
!  </IN>
!  <IN NAME="aerosolpress" TYPE="real">
!   required in sa_gcm mode, absent otherwise:
!                         pressure field to be used by aerosol param-
!                         eterization
!  </IN>
!  <IN NAME="kbot" TYPE="integer">
!   present when running eta vertical coordinate,
!                         index of lowest model level above ground
!  </IN>
! </SUBROUTINE>
!
subroutine define_atmos_input_fields (is, ie, js, je, pfull, phalf, &
                                      t, q, ts, Atmos_input, &
                                      cloudtemp, cloudvapor, &
                                      aerosoltemp, aerosolvapor, &
                                      aerosolpress, kbot)     

!---------------------------------------------------------------------
!    define_atmos_input_fields converts the atmospheric input fields 
!    (pfull, phalf, t, q, ts) to the form needed by the radiation 
!    modules, and when needed returns radiation-ready fields of pressure
!    (press, psfc), temperature (temp, tsfc), water vapor mixing ratio 
!    (rh2o) and several auxiliary variables in the derived type 
!    structure Atmos_input. the optional input variables are present
!    when running radiative feedback studies (sa_model), and are needed
!    to allow variation of temperature and vapor fields while holding 
!    the aerosol and cloud amounts fixed.
!---------------------------------------------------------------------
     
integer,                 intent(in)              :: is, ie, js, je
real, dimension(:,:,:),  intent(in)              :: pfull, phalf, t, q
real, dimension(:,:),    intent(in)              :: ts
type(atmos_input_type),  intent(inout)           :: Atmos_input
integer, dimension(:,:), intent(in), optional    :: kbot
real, dimension(:,:,:),  intent(in), optional    :: cloudtemp,    &
                                                    cloudvapor, &
                                                    aerosoltemp, &
                                                    aerosolvapor, &
                                                    aerosolpress

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
!
!   intent(out) variables:
!
!      Atmos_input   atmos_input type structure, contains the 
!                    following components defined in this subroutine
!         psfc          surface pressure 
!                       [ (kg /( m s^2) ] 
!         tsfc          surface temperature
!                       [ deg K ]
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
!         aerosoltemp   temperature to be seen by aerosols (used in 
!                       sa_gcm feedback studies) 
!                       [ degrees K ]
!         aerosolvapor  water vapor to be seen by aerosols (used in 
!                       sa_gcm feedback studies) 
!                       [ nondimensional ]
!         aerosolpress  pressure field to be seen by aerosols (used in 
!                       sa_gcm feedback studies) 
!                       [ Pa ]
!         aerosolrelhum relative humidity seen by aerosol package,
!                       used in sa_gcm feedback studies
!                       [ dimensionless ]
!
!   intent(in), optional variables:
!
!      kbot               present when running eta vertical coordinate,
!                         index of lowest model level above ground (???)
!      cloudtemp          temperature to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!                         [ degrees K ]
!      cloudvapor         water vapor to be seen by clouds (used in 
!                         sa_gcm feedback studies) 
!                         [ nondimensional ]
!      aerosoltemp        required in sa_gcm mode, absent otherwise:
!                         temperature field to be used by aerosol param-
!                         eterization 
!      aerosolvapor       required in sa_gcm mode, absent otherwise:
!                         water vapor field to be used by aerosol param-
!                         eterization 
!      aerosolpress       required in sa_gcm mode, absent otherwise:
!                         pressure field to be used by aerosol param-
!                         eterization
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables
 
      integer :: i, j, k, kb
      integer :: kmax

!---------------------------------------------------------------------
!  local variables
!
!     i, j, k      do loop indices
!     kb           vertical index of lowest atmospheric level (when
!                  using eta coordinates)
!     kmax         number of model layers
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!-------------------------------------------------------------------
!    verify that optional arguments that are required when running in
!    sa-model mode are present.
!-------------------------------------------------------------------
      if (Environment%running_sa_model) then
        if (present(cloudtemp)   .and.   &
            present(cloudvapor)  .and.   &
            present(aerosoltemp)  .and.   &
            present(aerosolvapor)  .and.   &
            present(aerosolpress) ) then        
        else
            call error_mesg ('radiation_driver_mod', &
                    'must pass optional arguments  to'//&
                    'radiation_driver when running sa model', FATAL)
        endif
      endif

!----------------------------------------------------------------------
!    define the number of model layers.
!----------------------------------------------------------------------
      kmax = size(t,3)

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
      allocate ( Atmos_input%cloudtemp(size(t,1), size(t,2),   &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%cloudvapor(size(t,1), size(t,2),   &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%clouddeltaz(size(t,1), size(t,2),   &
                                                         size(t,3)  ) )
      allocate ( Atmos_input%aerosoltemp(size(t,1), size(t,2),   &
                                                    size(t,3)  ) )
      allocate ( Atmos_input%aerosolpress(size(t,1), size(t,2),    &
                                                     size(t,3)+1) )
      allocate ( Atmos_input%aerosolvapor(size(t,1), size(t,2),   &
                                                     size(t,3)  ) )
      allocate ( Atmos_input%aerosolrelhum(size(t,1), size(t,2),   &
                                                      size(t,3)  ) )
      allocate ( Atmos_input%deltaz(size(t,1), size(t,2), size(t,3) ) )
      allocate ( Atmos_input%pflux (size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%tflux (size(t,1), size(t,2), size(t,3)+1) )
      allocate ( Atmos_input%psfc (size(t,1), size(t,2)             ) )
      allocate ( Atmos_input%tsfc (size(t,1), size(t,2)             ) )

!---------------------------------------------------------------------
!    define the cloudtemp component of Atmos_input. 
!---------------------------------------------------------------------
      if (present (cloudtemp) ) then
        Atmos_input%cloudtemp(:,:,:)   = cloudtemp(:,:,:)
      else
        Atmos_input%cloudtemp(:,:,:)   = t(:,:,:)
      endif

!---------------------------------------------------------------------
!    define the cloudvapor component of Atmos_input.
!---------------------------------------------------------------------
      if (present (cloudvapor) ) then
        Atmos_input%cloudvapor(:,:,:)   = cloudvapor(:,:,:)
      else
        Atmos_input%cloudvapor(:,:,:)   = q(:,:,:)
      endif

!---------------------------------------------------------------------
!    define the aerosoltemp component of Atmos_input.
!---------------------------------------------------------------------
      if (present (aerosoltemp) ) then
        Atmos_input%aerosoltemp(:,:,:)   = aerosoltemp(:,:,:)
      else
        Atmos_input%aerosoltemp(:,:,:)   = t(:,:,:)
      endif
 
!---------------------------------------------------------------------
!    define the aerosolvapor component of Atmos_input.
!---------------------------------------------------------------------
      if (present (aerosolvapor) ) then
        Atmos_input%aerosolvapor(:,:,:)   = aerosolvapor(:,:,:)
      else
        Atmos_input%aerosolvapor(:,:,:)   = q(:,:,:)
      endif

!---------------------------------------------------------------------
!    define values of surface pressure and temperature.
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

!------------------------------------------------------------------
!    define the atmospheric pressure and temperature arrays.
!------------------------------------------------------------------
      do k=1,kmax 
        Atmos_input%press(:,:,k) = pfull(:,:,k)
        Atmos_input%phalf(:,:,k) = phalf(:,:,k)
        Atmos_input%temp (:,:,k) = t(:,:,k)
      end do
      Atmos_input%press(:,:,kmax+1) = phalf(:,:,kmax+1)
      Atmos_input%phalf(:,:,kmax+1) = phalf(:,:,kmax+1)
      Atmos_input%temp (:,:,kmax+1) = ts  (:,:)

!---------------------------------------------------------------------
!    define the aerosolpress component of Atmos_input.
!---------------------------------------------------------------------
      if (present (aerosolpress) ) then
        do k=1,kmax
          Atmos_input%aerosolpress(:,:,k)   = aerosolpress(:,:,k)
        end do
      else
        do k=1,kmax
          Atmos_input%aerosolpress(:,:,k)   = pfull(:,:,k)
        end do
      endif
      Atmos_input%aerosolpress(:,:,kmax+1)   = phalf(:,:,kmax+1)
 
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
          Atmos_input%aerosolvapor(:,:,:) =    &
                                 Atmos_input%aerosolvapor(:,:,:)/  &
                            (1.0 - Atmos_input%aerosolvapor(:,:,:))
        endif
      else if (Environment%running_standalone) then
        Atmos_input%rh2o (:,:,:) = q(:,:,:)
      endif ! (running gcm)
 
!------------------------------------------------------------------
!    be sure that the magnitude of the water vapor mixing ratio field 
!    to be input to the radiation code is no smaller than the value of 
!    rh2o_lower_limit, which is 2.0E-07 when running the sea_esf
!    radiation code and 3.0e-06 when running the original radiation
!    code. Likewise, the temperature that the radiation code sees is
!    constrained to lie between 100K and 370K. these are the limits of
!    the tables referenced within the radiation package.
!-----------------------------------------------------------------------
      if (do_rad) then
        Atmos_input%rh2o(:,:,ks:ke) =    &
            MAX(Atmos_input%rh2o(:,:,ks:ke), rh2o_lower_limit)
        Atmos_input%cloudvapor(:,:,ks:ke) =    &
            MAX(Atmos_input%cloudvapor(:,:,ks:ke), rh2o_lower_limit)
        Atmos_input%aerosolvapor(:,:,ks:ke) =    &
                    MAX(Atmos_input%aerosolvapor(:,:,ks:ke), rh2o_lower_limit)
        Atmos_input%temp(:,:,ks:ke) =     &
                    MAX(Atmos_input%temp(:,:,ks:ke), temp_lower_limit)
        Atmos_input%temp(:,:,ks:ke) =     &
                    MIN(Atmos_input%temp(:,:,ks:ke), temp_upper_limit)
        Atmos_input%cloudtemp(:,:,ks:ke) =     &
                    MAX(Atmos_input%cloudtemp(:,:,ks:ke), temp_lower_limit)
        Atmos_input%cloudtemp(:,:,ks:ke) =     &
                    MIN(Atmos_input%cloudtemp(:,:,ks:ke), temp_upper_limit)
        Atmos_input%aerosoltemp(:,:,ks:ke) =     &
                    MAX(Atmos_input%aerosoltemp(:,:,ks:ke), temp_lower_limit)
     Atmos_input%aerosoltemp(:,:,ks:ke) =     &
                    MIN(Atmos_input%aerosoltemp(:,:,ks:ke), temp_upper_limit)
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




!#####################################################################
! <SUBROUTINE NAME="define_surface">
!  <OVERVIEW>
!    define_surface stores the input values of land fraction and 
!    surface albedo in a surface_type structure Surface. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    define_surface stores the input values of land fraction and 
!    surface albedo in a surface_type structure Surface. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call define_surface (is, ie, js, je, albedo, land, Surface)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!    starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="albedo" TYPE="real">
!   surface albedo
!  </IN>
!  <IN NAME="land" TYPE="real">
!   fraction of grid box which is land 
!  </IN>
!  <INOUT NAME="Surface" TYPE="surface_type">
!   surface_type structure to be valued
!  </INOUT>
! </SUBROUTINE>
!
subroutine define_surface (is, ie, js, je, albedo, land, Surface)

!---------------------------------------------------------------------
!    define_surface stores the input values of land fraction and 
!    surface albedo in a surface_type structure Surface.  
!---------------------------------------------------------------------
     
integer,                 intent(in)              :: is, ie, js, je
real, dimension(:,:),    intent(in)              :: albedo, land 
type(surface_type),      intent(inout)           :: Surface     

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      albedo       surface albedo  [ dimensionless ]
!      land         fraction of grid box which is land [ dimensionless ]
!
!   intent(out) variables:
!
!      Surface       surface_type structure, contains the 
!                    following components defined in this subroutine
!         asfc          surface albedo
!                       [ non-dimensional ]
!         land          fraction of grid box covered by land
!                       [ non-dimensional ]
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that the module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    allocate space for the components of the derived type variable
!    Surface.     
!---------------------------------------------------------------------
      allocate (Surface%asfc (size(albedo,1), size(albedo,2)) )
      allocate (Surface%land (size(albedo,1), size(albedo,2)) )

 
!------------------------------------------------------------------
!    define the fractional land area of each grid box and the surface
!    albedo from the input argument values.
!------------------------------------------------------------------
      Surface%land(:,:) = land(:,:)
      Surface%asfc(:,:) = albedo(:,:)
     

!----------------------------------------------------------------------


end subroutine define_surface    



!#####################################################################
! <SUBROUTINE NAME="surface_dealloc">
!  <OVERVIEW>
!    surface_dealloc deallocates the array components of the
!    surface_type structure Surface.
!  </OVERVIEW>
!  <DESCRIPTION>
!    surface_dealloc deallocates the array components of the
!    surface_type structure Surface.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call surface_dealloc (Surface)
!  </TEMPLATE>
!  <INOUT NAME="Surface" TYPE="surface_type">
!   surface_type structure to be deallocated
!  </INOUT>
! </SUBROUTINE>
!
subroutine surface_dealloc (Surface)

!----------------------------------------------------------------------
!    surface_dealloc deallocates the array components of the
!    surface_type structure Surface.
!----------------------------------------------------------------------

type(surface_type), intent(inout) :: Surface

!--------------------------------------------------------------------
!   intent(inout) variable:
!
!      Surface        surface_type structure, contains variables 
!                     defining the surface albedo and land fraction
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!-------------------------------------------------------------------
!    deallocate components of surface_type structure.
!-------------------------------------------------------------------
      deallocate (Surface%asfc)
      deallocate (Surface%land)

!--------------------------------------------------------------------


end subroutine surface_dealloc 



!#####################################################################
! <SUBROUTINE NAME="atmos_input_dealloc">
!  <OVERVIEW>
!    atmos_input_dealloc deallocates the array components of the
!    atmos_input_type structure Atmos_input.
!  </OVERVIEW>
!  <DESCRIPTION>
!    atmos_input_dealloc deallocates the array components of the
!    atmos_input_type structure Atmos_input.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call atmos_input_dealloc (Atmos_input)
!  </TEMPLATE>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!      atmos_input_type structure, contains variables 
!                     defining the atmospheric pressure, temperature
!                     and moisture distribution.
!  </INOUT>
! </SUBROUTINE>
!
subroutine atmos_input_dealloc (Atmos_input)

!----------------------------------------------------------------------
!    atmos_input_dealloc deallocates the array components of the
!    atmos_input_type structure Atmos_input.
!----------------------------------------------------------------------

type(atmos_input_type), intent(inout) :: Atmos_input

!--------------------------------------------------------------------
!   intent(inout) variable:
!
!      Atmos_input    atmos_input_type structure, contains variables 
!                     defining the atmospheric pressure, temperature
!                     and moisture distribution.
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!---------------------------------------------------------------------
!    deallocate components of atmos_input_type structure.
!---------------------------------------------------------------------
      deallocate (Atmos_input%press      )
      deallocate (Atmos_input%phalf      )
      deallocate (Atmos_input%temp       )
      deallocate (Atmos_input%rh2o       )
      deallocate (Atmos_input%rel_hum    )
      deallocate (Atmos_input%pflux      )
      deallocate (Atmos_input%tflux      )
      deallocate (Atmos_input%deltaz     )
      deallocate (Atmos_input%psfc       )
      deallocate (Atmos_input%tsfc       )
      deallocate (Atmos_input%cloudtemp  )
      deallocate (Atmos_input%cloudvapor )
      deallocate (Atmos_input%clouddeltaz)
      deallocate (Atmos_input%aerosoltemp)
      deallocate (Atmos_input%aerosolvapor )
      deallocate (Atmos_input%aerosolpress )
      deallocate (Atmos_input%aerosolrelhum )

!--------------------------------------------------------------------


end subroutine atmos_input_dealloc 




!#####################################################################
! <SUBROUTINE NAME="radiation_driver_end">
!  <OVERVIEW>
!   radiation_driver_end is the destructor for radiation_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   radiation_driver_end is the destructor for radiation_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_driver_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine radiation_driver_end

!----------------------------------------------------------------------
!    radiation_driver_end is the destructor for radiation_driver_mod.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables

      integer :: unit

!-------------------------------------------------------------------
!    verify that this module has been initialized. if not, exit.
!-------------------------------------------------------------------
      if (.not. module_is_initialized)   &
          call error_mesg ('radiation_driver_mod',  &
               'module has not been initialized', FATAL)

!----------------------------------------------------------------------
!    when running in gcm, write a restart file. this is not done in the
!    standalone case.
!---------------------------------------------------------------------
      if (Environment%running_gcm) then
        unit = open_restart_file   &
                         ('RESTART/radiation_driver.res', 'write')

!---------------------------------------------------------------------
!    only the root pe will write control information -- the last value 
!    in the list of restart versions and the alarm information.
!---------------------------------------------------------------------
        if (mpp_pe() == mpp_root_pe() ) then
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
        if (mpp_pe() == mpp_root_pe() ) then
          write (unit) renormalize_sw_fluxes, do_clear_sky_pass
        endif

!---------------------------------------------------------------------
!    write out the optional shortwave renormalization data. 
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes) then   
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
        endif    ! (renormalize)

!---------------------------------------------------------------------
!    close the radiation_driver.res file
!---------------------------------------------------------------------
        call close_file (unit)
      endif   ! (running_gcm)

!---------------------------------------------------------------------
!    wrap up modules initialized by this module.
!---------------------------------------------------------------------
      call astronomy_end

!---------------------------------------------------------------------
!    wrap up modules specific to the radiation package in use.
!---------------------------------------------------------------------
      if (do_sea_esf_rad) then
        call cloudrad_package_end
        call aerosolrad_package_end
        call rad_output_file_end
        call sea_esf_rad_end
      else
        call original_fms_rad_end
      endif


!---------------------------------------------------------------------
!    release space for renormalization arrays, if that option is active.
!---------------------------------------------------------------------
      if (renormalize_sw_fluxes .or. all_step_diagnostics)  then
        deallocate (solar_save, flux_sw_surf_save, sw_heating_save, &
                    tot_heating_save, dfsw_save, ufsw_save,   &
                    fsw_save, hsw_save)
        if (do_clear_sky_pass) then
          deallocate (sw_heating_clr_save, tot_heating_clr_save,  &
                      dfswcf_save, ufswcf_save, fswcf_save,   &
                      hswcf_save)
        endif
      endif

!---------------------------------------------------------------------
!    release space needed when all_step_diagnostics is active.
!---------------------------------------------------------------------
      if (all_step_diagnostics)  then
        deallocate (olr_save, lwups_save, lwdns_save, tdtlw_save)
        if (do_clear_sky_pass) then
          deallocate (olr_clr_save, lwups_clr_save, lwdns_clr_save, &
                      tdtlw_clr_save)
        endif
      endif

!---------------------------------------------------------------------
!    release space used for module variables that hold data between
!    timesteps.
!---------------------------------------------------------------------
      if (Environment%running_gcm  .or. &
          Environment%running_sa_model) then
        deallocate (Rad_output%tdt_rad, Rad_output%tdt_rad_clr,  & 
                    Rad_output%tdtsw, Rad_output%tdtsw_clr, &
                    Rad_output%tdtlw, Rad_output%flux_sw_surf,  &
                    Rad_output%flux_lw_surf, Rad_output%coszen_angle)
      endif    ! (running_gcm)

!---------------------------------------------------------------------
!    call rad_utilities_end to uninitialize that module.
!---------------------------------------------------------------------
        call rad_utilities_end

!----------------------------------------------------------------------
!    set initialization status flag.
!----------------------------------------------------------------------
      module_is_initialized = .false.



end subroutine radiation_driver_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!#####################################################################
! <SUBROUTINE NAME="read_restart_file">
!  <OVERVIEW>
!    read_restart_file reads a restart file containing radiation
!    restart information. it may be either a radiation_driver.res, or
!    an older sea_esf_rad.res file.
!  </OVERVIEW>
!  <DESCRIPTION>
!    read_restart_file reads a restart file containing radiation
!    restart information. it may be either a radiation_driver.res, or
!    an older sea_esf_rad.res file.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_restart_file
!  </TEMPLATE>
! </SUBROUTINE>
!
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
      integer               :: kmax
      integer               :: new_rad_time, old_time_step

!--------------------------------------------------------------------
!   local variables:
!
!       unit              i/o unit number connected to .res file
!       end               logical variable indicating, if true, that
!                         end of file has been reached on the current
!                         read operation
!       avg_present       if true, time-average data is present in the
!                         restart file
!       renorm_present    if true, sw renormalization data is present 
!                         in the restart file
!       cldfree_present   if true, and if renorm_present is true, then
!                         the clear-sky sw renormalization data is
!                         present in the restart file
!       chvers            character form of restart_version (i4)
!       avg_gases         if true, then time-average data for radiative
!                         gases is present in restart file
!       avg_clouds        if true, then time-average data for clouds is
!                         present in restart file
!       null              dummy array used as location to read older
!                         restart version data into           
!       vers              version number of the restart file being read
!       kmax              number of model layers
!       new_rad_time      time remaining until next radiation calcul-
!                         ation; replaces the rad_alarm value read from
!                         restart file when the radiation timestep
!                         changes upon restart  
!       old_time_step     radiation timestep that was used in job 
!                         which wrote the restart file
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if one is using the sea_esf_rad package and there is a 
!    sea_esf_rad.res restart file present in the input directory, it 
!    must be version 1 in order to be readable by the current module.
!    otherwise, exit.
!---------------------------------------------------------------------
      if (do_sea_esf_rad .and. file_exist('INPUT/sea_esf_rad.res')) then
        unit = open_restart_file ('INPUT/sea_esf_rad.res', 'read')
        read (unit) vers
        if ( vers /= 1 ) then
          write (chvers,'(i4)') vers
          call error_mesg ('radiation_driver_mod', &
             'restart version '//chvers//' cannot be read '//&
              'by this module version', FATAL)
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
        unit = open_restart_file ('INPUT/radiation_driver.res', 'read')
        read (unit) vers
        if (do_sea_esf_rad) then
          if ( .not. any(vers == restart_versions_sea) ) then
            write (chvers,'(i4)') vers
            call error_mesg ('radiation_driver_mod', &
                    'restart version '//chvers//' cannot be read '//&
                 'by this module version when using sea_esf_rad '//&
                                            'radiation_package', FATAL)
          endif
        else
          if ( .not. any(vers == restart_versions) ) then
            write (chvers,'(i4)') vers
            call error_mesg ('radiation_driver_mod', &
                    'restart version '//chvers//' cannot be read '//&
               'by this module version when using the original_'//&
                                         'fms radiation package', FATAL)
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
        old_time_step = SECONDS_PER_DAY*null(4) + null(3)
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

!----------------------------------------------------------------------
!    versions 3 and 4 include variables needed when sw renormalization 
!    is active, and logical variables indicating which additional fields
!    are present.
!----------------------------------------------------------------------
      if (vers == 3 .or. vers == 4) then

!---------------------------------------------------------------------
!    determine if accumulation arrays are present in the restart file.
!    if input fields are to be time-averaged, read the values from the
!    files.  note that avg_present and renorm_present cannot both be 
!    true.
!---------------------------------------------------------------------
        read (unit) avg_present, renorm_present, cldfree_present
  
!---------------------------------------------------------------------
!    if renormalize_sw_fluxes is true and the data is present in the
!    restart file, read it.
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes) then  
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
        endif ! (renormalize)
      else if (vers == 5) then

!---------------------------------------------------------------------
!    determine if accumulation arrays are present in the restart file.
!    if input fields are to be time-averaged, read the values from the
!    files.  note that avg_present and renorm_present cannot both be 
!    true.
!---------------------------------------------------------------------
        read (unit) renorm_present, cldfree_present
  
!---------------------------------------------------------------------
!    if renormalize_sw_fluxes is true and the data is present in the
!    restart file, read it.
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes) then  
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
        endif ! (renormalize)
      endif    ! (vers == 3 or 4)

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
        if (mpp_pe() == mpp_root_pe() ) then
          call error_mesg ('radiation_driver_mod',          &
               'radiation will be called on first step of run', NOTE)
        endif
      else
        if (rad_time_step /= old_time_step ) then
          new_rad_time = rad_alarm - old_time_step + rad_time_step
          if ( new_rad_time > 0 ) then
            if (mpp_pe() == mpp_root_pe() ) then
              call error_mesg ('radiation_driver_mod',          &
                  'radiation time step has changed, therefore '//&
                  'next radiation time also changed', NOTE)
            endif
            rad_alarm = new_rad_time
          endif
        endif  
      endif   ! (rad_alarm == 1)

!--------------------------------------------------------------------


end subroutine read_restart_file



!#####################################################################
! <SUBROUTINE NAME="initialize_diagnostic_integrals">
!  <OVERVIEW>
!    initialize_diagnostic_integrals registers the desired integrals 
!    with diag_integral_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    initialize_diagnostic_integrals registers the desired integrals 
!    with diag_integral_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call initialize_diagnostic_integrals
!  </TEMPLATE>
! </SUBROUTINE>
!
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
! <SUBROUTINE NAME="diag_field_init">
!  <OVERVIEW>
!    diag_field_init registers the desired diagnostic fields with the
!    diagnostics manager.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diag_field_init registers the desired diagnostic fields with the
!    diagnostics manager.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_field_init ( Time, axes )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   Current time
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   diagnostic variable axes for netcdf files
!  </IN>
! </SUBROUTINE>
!
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

!--------------------------------------------------------------------
!  local variables:
!
!       clr          character string used in netcdf variable short name
!       clr2         character string used in netcdf variable long name
!       n            number of passes through name generation loop
!       i            do-loop index
!
!--------------------------------------------------------------------

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



!######################################################################
! <SUBROUTINE NAME="obtain_astronomy_variables">
!  <OVERVIEW>
!    obtain_astronomy_variables retrieves astronomical variables, valid 
!    at the requested time and over the requested time intervals.
!  </OVERVIEW>
!  <DESCRIPTION>
!    obtain_astronomy_variables retrieves astronomical variables, valid 
!    at the requested time and over the requested time intervals.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call obtain_astronomy_variables (is, ie, js, je, lat, lon,     &
!                                       Astro, Astro2)  
!
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 
!  <INOUT NAME="Astro" TYPE="astronomy_type">
!     astronomy_type structure; It will
!                    be used to determine the insolation at toa seen
!                    by the shortwave radiation code
!  </INOUT>
!  <INOUT NAME="Astro2" TYPE="astronomy_type">
!     astronomy_type structure, defined when renormal-
!                    ization is active. the same components are defined
!                    as for Astro, but they are valid over the current
!                    physics timestep.
!  </INOUT>
!  <IN NAME="lon" TYPE="real">
!    lon        mean longitude (in radians) of all grid boxes processed by
!               this call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
! </SUBROUTINE>
!
subroutine obtain_astronomy_variables (is, ie, js, je, lat, lon,     &
                                       Astro, Astro2)  

!---------------------------------------------------------------------
!    obtain_astronomy_variables retrieves astronomical variables, valid 
!    at the requested time and over the requested time intervals.
!---------------------------------------------------------------------
integer,                     intent(in)    ::  is, ie, js, je
real, dimension(:,:),        intent(in)    ::  lat, lon
type(astronomy_type),        intent(inout) ::  Astro, Astro2

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      lat          latitude of model points  
!                   [ radians ]
!      lon          longitude of model points 
!                   [ radians ]
!
!   intent(inout) variables:
!
!      Astro         astronomy_type structure; contains the following
!                    components defined in this subroutine that will
!                    be used to determine the insolation at toa seen
!                    by the shortwave radiation code 
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
!                       relative to the mean square of earth-sun 
!                       distance
!                       [ non-dimensional ]
!
!      Astro2        astronomy_type structure, defined when renormal-
!                    ization is active. the same components are defined
!                    as for Astro, but they are valid over the current
!                    physics timestep.
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      type(time_type)                   :: Dt_zen, Dt_zen2
      real, dimension(ie-is+1, je-js+1) ::                            &
                                           cosz_r, solar_r, fracday_r, &
                                           cosz_p, solar_p, fracday_p, &
                                           cosz_a, solar_a, fracday_a
      real                              :: rrsun_r, rrsun_p, rrsun_a
      

!--------------------------------------------------------------------
!  local variables:
!
!     Dt_zen        time-type variable containing the components of the
!                   radiation time step, needed unless do_average is
!                   true or this is not a radiation step and renormal-
!                   ize_sw_fluxes is true
!     Dt_zen2       time-type variable containing the components of the
!                   physics time step, needed when renormalize_sw_fluxes
!                   or do_average is true
!     cosz_r        cosine of zenith angle --  mean value over
!                   radiation time step            
!                   [ non-dimensional ]
!     solar_r       shortwave flux factor relevant over radiation time
!                   step: cosine of zenith angle * daylight fraction / 
!                   (earth-sun distance squared)
!                   [ non-dimensional ]
!     fracday_r     fraction of timestep during which the sun is 
!                   shining over radiation time step
!                   [ non-dimensional ]
!     cosz_p        cosine of zenith angle --  mean value over
!                   physics time step            
!                   [ non-dimensional ]
!     solar_p       shortwave flux factor relevant over physics time
!                   step: cosine of zenith angle * daylight fraction / 
!                   (earth-sun distance squared)
!                   [ non-dimensional ]
!     fracday_p     fraction of timestep during which the sun is 
!                   shining over physics time step
!                   [ non-dimensional ]
!     cosz_a        cosine of zenith angle --  mean value over
!                   next radiation time step            
!                   [ non-dimensional ]
!     solar_a       shortwave flux factor relevant over next radiation 
!                   time step: cosine of zenith angle * daylight 
!                   fraction / (earth-sun distance squared)
!                   [ non-dimensional ]
!     fracday_a     fraction of timestep during which the sun is 
!                   shining over next radiation time step
!                   [ non-dimensional ]
!     rrsun_r       inverse of square of earth-sun distance, 
!                   relative to the mean square of earth-sun 
!                   distance, valid over radiation time step
!                   [ non-dimensional ]
!     rrsun_p       inverse of square of earth-sun distance, 
!                   relative to the mean square of earth-sun 
!                   distance, valid over physics time step
!                   [ non-dimensional ]
!     rrsun_a       inverse of square of earth-sun distance, 
!                   relative to the mean square of earth-sun 
!                   distance, valid over next radiation time step
!                   [ non-dimensional ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    allocate the components of the astronomy_type structure which will
!    return the astronomical inputs to radiation (cosine of zenith 
!    angle, daylight fraction, solar flux factor and earth-sun distance)
!    that are to be used on the current step.
!---------------------------------------------------------------------
      allocate ( Astro%cosz   (size(lat,1), size(lat,2) ) )
      allocate ( Astro%fracday(size(lat,1), size(lat,2) ) )
      allocate ( Astro%solar  (size(lat,1), size(lat,2) ) )

!---------------------------------------------------------------------
!    case 1: diurnally-varying shortwave radiation.
!---------------------------------------------------------------------
      if (Sw_control%do_diurnal) then

!-------------------------------------------------------------------
!    convert the radiation timestep and the model physics timestep
!    to time_type variables.
!-------------------------------------------------------------------
        Dt_zen  = set_time (rad_time_step, 0)
        Dt_zen2 = set_time (dt, 0)
        
!---------------------------------------------------------------------
!    calculate the astronomical factors averaged over the radiation time
!    step between Rad_time and Rad_time + Dt_zen. these values are 
!    needed on radiation steps. output is stored in Astro_rad.
!---------------------------------------------------------------------
        if (do_rad) then
          call diurnal_solar (lat, lon, Rad_time, cosz_r, fracday_r, &
                              rrsun_r, dt_time=Dt_zen)
          fracday_r = MIN (fracday_r, 1.00)
          solar_r = cosz_r*fracday_r*rrsun_r
        endif

!---------------------------------------------------------------------
!    calculate the astronomical factors averaged over the physics time
!    step between Rad_time and Rad_time + Dt_zen2. these values are
!    needed if either renormalization or time-averaging is active. store
!    the astronomical outputs in Astro_phys.
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes) then
          call diurnal_solar (lat, lon, Rad_time, cosz_p, fracday_p, &
                              rrsun_p, dt_time=Dt_zen2)
          fracday_p = MIN (fracday_p, 1.00)
          solar_p = cosz_p*fracday_p*rrsun_p
        endif

!--------------------------------------------------------------------
!    define the astronomy_type variable(s) to be returned and used in 
!    the radiation calculation. Astro contains the values to be used
!    in the radiation calculation, Astro2 contains values relevant 
!    over the current physics timestep and is used for renormalization.
!    when renormalization is active, the physics step set is always 
!    needed, and in addition on radiation steps, the radiation step
!    values are needed. 
!---------------------------------------------------------------------
        if (renormalize_sw_fluxes) then
          if (.not. do_rad) then
            Astro%cosz    = cosz_p
            Astro%fracday = fracday_p
            Astro%solar   = solar_p
            Astro%rrsun   = rrsun_p
          else 
            Astro%cosz    = cosz_r
            Astro%fracday = fracday_r
            Astro%solar   = solar_r
            Astro%rrsun   = rrsun_r
            allocate ( Astro2%fracday(size(lat,1), size(lat,2) ) )
            allocate ( Astro2%cosz   (size(lat,1), size(lat,2) ) )
            allocate ( Astro2%solar  (size(lat,1), size(lat,2) ) )
            Astro2%cosz    = cosz_p
            Astro2%fracday = fracday_p
            Astro2%solar   = solar_p
            Astro2%rrsun   = rrsun_p
          endif

!---------------------------------------------------------------------
!    if renormalization is active, then only the values applicable over
!    radiation steps are needed. 
!---------------------------------------------------------------------
        else                 
          Astro%cosz    = cosz_r
          Astro%fracday = fracday_r
          Astro%solar   = solar_r
          Astro%rrsun   = rrsun_r
        endif

!---------------------------------------------------------------------
!    when in the gcm and on a radiation calculation step, define cosine
!    of zenith angle valid over the next radiation step. this is needed 
!    so that the ocean albedo (function of zenith angle) may be properly
!    defined and provided as input to the radiation package on the next
!    timestep.
!----------------------------------------------------------------------
        if (Environment%running_gcm) then       
          if (do_rad) then
            call diurnal_solar (lat, lon, Rad_time+Dt_zen, cosz_a,   &
                                fracday_a, rrsun_a, dt_time=Dt_zen)
            Rad_output%coszen_angle(is:ie,js:je) = cosz_a(:,:)
          endif  ! (do_rad)
        endif ! (running_gcm)

!---------------------------------------------------------------------
!    case 2: annual-mean shortwave radiation.
!---------------------------------------------------------------------
      else if (Sw_control%do_annual) then
        call annual_mean_solar (js, je, lat, Astro%cosz, Astro%solar, &
                                Astro%fracday, Astro%rrsun)

!---------------------------------------------------------------------
!    save the cosine of zenith angle on the current step to be used to 
!    calculate ocean albedo for use on the next radiation timestep.
!---------------------------------------------------------------------
        if (Environment%running_gcm) then       
          Rad_output%coszen_angle(is:ie,js:je) = Astro%cosz(:,:)
        endif

!---------------------------------------------------------------------
!    case 3: daily-mean shortwave radiation.
!---------------------------------------------------------------------
      else if (Sw_control%do_daily_mean) then
        call daily_mean_solar (lat, Rad_time, Astro%cosz,  &
                               Astro%fracday, Astro%rrsun)
        Astro%solar = Astro%cosz*Astro%rrsun*Astro%fracday

!---------------------------------------------------------------------
!    save the cosine of zenith angle on the current step to be used to 
!    calculate ocean albedo for use on the next radiation timestep.
!---------------------------------------------------------------------
        if (Environment%running_gcm) then       
          Rad_output%coszen_angle(is:ie,js:je) = Astro%cosz(:,:)
        endif

!----------------------------------------------------------------------
!    if none of the above options are active, write an error message and
!    stop execution.
!----------------------------------------------------------------------
      else
        call error_mesg('radiation_driver_mod', &
             ' no valid zenith angle specification', FATAL)
      endif

!-------------------------------------------------------------------


end subroutine obtain_astronomy_variables 



!####################################################################
! <SUBROUTINE NAME="radiation_calc">
!  <OVERVIEW>
!    radiation_calc is called on radiation timesteps and calculates
!    the long- and short-wave radiative fluxes and heating rates, and
!    obtains the radiation output fields needed in other portions of
!    the model.
!  </OVERVIEW>
!  <DESCRIPTION>
!    radiation_calc is called on radiation timesteps and calculates
!    the long- and short-wave radiative fluxes and heating rates, and
!    obtains the radiation output fields needed in other portions of
!    the model.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call radiation_calc (is, ie, js, je, Rad_time, Time_diag,  &
!                           lat, lon, Atmos_input, Surface, Rad_gases, &
!                           Aerosol_props, Aerosol, Cldrad_props,   &
!                           Cld_spec, Astro, Rad_output, Lw_output,   &
!                           Sw_output, Fsrad_output, mask, kbot)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 
!  <IN NAME="Rad_time" TYPE="time_type">
!      Rad_time          time at which the radiative fluxes are to apply
!                        [ time_type (days, seconds) ] 
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ] 
!  </IN>
!  <IN NAME="lon" TYPE="real">
!    lon        mean longitude (in radians) of all grid boxes processed by
!               this call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="Surface" TYPE="surface_type">
!   Surface input data to radiation package
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!   Atmospheric input data to radiation package
!  </IN>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatological input data to radiation package
!  </IN>
!  <INOUT NAME="Aerosol_props" TYPE="aerosol_properties_type">
!   Aerosol radiative properties
!  </INOUT>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   Cloud radiative properties
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   Cloud microphysical and physical parameters to radiation package, 
!                     contains var-
!                     iables defining the cloud distribution, passed 
!                     through to lower level routines
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!   astronomical input data for the radiation package
!  </IN>
!  <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!   Radiative gases properties to radiation package, , contains var-
!                     iables defining the radiatively active gases, 
!                     passed through to lower level routines
!  </IN>
!  <INOUT NAME="Rad_output" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </INOUT>
!  <INOUT NAME="Lw_output" TYPE="lw_output_type">
!      longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!  </INOUT>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </INOUT>
!  <INOUT NAME="Fsrad_output" TYPE="Fsrad_output_type">
!   radiation output data from the original_fms_rad
!                        radiation package, when that package 
!                        is active
!  </INOUT>
!  <IN NAME="kbot" TYPE="integer">
!   OPTIONAL: present when running eta vertical coordinate,
!                        index of lowest model level above ground
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
subroutine radiation_calc (is, ie, js, je, Rad_time, Time_diag,  &
                           lat, lon, Atmos_input, Surface, Rad_gases, &
                           Aerosol_props, Aerosol, Cldrad_props,   &
                           Cld_spec, Astro, Rad_output, Lw_output,   &
                           Sw_output, Fsrad_output, mask, kbot)

!--------------------------------------------------------------------
!    radiation_calc is called on radiation timesteps and calculates
!    the long- and short-wave radiative fluxes and heating rates, and
!    obtains the radiation output fields needed in other portions of
!    the model.
!-----------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                      intent(in)             :: is, ie, js, je
type(time_type),              intent(in)             :: Rad_time,   &
                                                        Time_diag
real, dimension(:,:),         intent(in)             :: lat, lon
type(atmos_input_type),       intent(in)             :: Atmos_input 
type(surface_type),           intent(in)             :: Surface
type(radiative_gases_type),   intent(in)             :: Rad_gases
type(aerosol_type),           intent(in)             :: Aerosol
type(aerosol_properties_type),intent(inout)          :: Aerosol_props
type(cldrad_properties_type), intent(in)             :: Cldrad_props
type(cld_specification_type), intent(in)             :: Cld_spec
type(astronomy_type),         intent(in)             :: Astro
type(rad_output_type),        intent(inout)          :: Rad_output
type(lw_output_type),         intent(inout)          :: Lw_output
type(sw_output_type),         intent(inout)          :: Sw_output
type(fsrad_output_type),      intent(inout)          :: Fsrad_output
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
!      lat               latitude of model points on model grid 
!                        [ radians ]
!      lon               longitude of model points on model grid 
!                        [ radians ]
!      Atmos_input       atmospheric input data for the radiation 
!                        package 
!                        [ atmos_input_type ]
!      Surface           surface input data to the radiation package
!                        [ surface_type ]
!      Rad_gases         radiative gas input data for the radiation 
!                        package
!                        [ radiative_gases_type ]
!      Aerosol_props     aerosol radiative property input data for the 
!                        radiation package
!                        [ aerosol_properties_type ]
!      Aerosol           aerosol input data to the radiation package
!                        [ aerosol_type ]
!      Cldrad_props      cloud radiative property input data for the 
!                        radiation package 
!                        [ cldrad_properties_type ]
!      Cld_spec          cloud specification input data for the 
!                        radiation package
!                        [ cld_specification_type ]
!      Astro             astronomical input data for the radiation 
!                        package 
!                        [ astronomy_type ]
!
!    intent(out) variables:
!
!      Rad_output        radiation output data needed by other modules
!                        [ rad_output_type ]
!      Lw_output         longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!                        [ lw_output_type ]
!          The following are the components of Lw_output:
!                 flxnet    net longwave flux at model flux levels 
!                           (including the ground and the top of the 
!                           atmosphere).
!                 heatra    longwave heating rates in model layers.
!                 flxnetcf  net longwave flux at model flux levels 
!                           (including the ground and the top of the 
!                           atmosphere) computed for cloud-free case.
!                 heatra    longwave heating rates in model layers 
!                           computed for cloud-free case.
!      Sw_output         shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!                        [ sw_output_type ]
!      Fsrad_output      radiation output data from the original_fms_rad
!                        radiation package, when that package 
!                        is active
!                        [ fsrad_output_type ]
!
!    intent(in), optional variables:
!
!      mask              present when running eta vertical coordinate,
!                        mask to define values at points below ground   
!      kbot              present when running eta vertical coordinate,
!                        index of lowest model level above ground 
!
!----------------------------------------------------------------------

      integer :: kmax

!---------------------------------------------------------------------
!    all_column_radiation and all_level_radiation are included as 
!    future controls which may be utiliized to execute the radiation
!    code on a grid other than the model grid. in the current release
!    however, both must be .true.. 
!---------------------------------------------------------------------
      if (all_column_radiation .and. all_level_radiation) then

!--------------------------------------------------------------------
!    call routines to perform radiation calculations, either using the
!    sea_esf_rad or original_fms_rad radiation package.
!---------------------------------------------------------------------
        if (do_sea_esf_rad) then
          call sea_esf_rad (is, ie, js, je, Atmos_input, Surface,  &
                            Astro, Rad_gases, Aerosol,  Aerosol_props,&
                            Cldrad_props, Cld_spec, Lw_output,   &
                            Sw_output)
        else  
          call original_fms_rad (is, ie, js, je, Atmos_input%phalf,  &
                                 lat, lon, do_clear_sky_pass,   &
                                 Rad_time, Time_diag, Atmos_input,  &
                                 Surface, Astro, Rad_gases,   &
                                 Cldrad_props, Cld_spec,    &
                                 Fsrad_output, mask=mask, kbot=kbot) 
        endif
      else  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    when this option is coded, replace this error_mesg code with
!    code which will map the input fields from the model grid to
!    the desired radiation grid. A preliminary version of code to per-
!    form this task (at least some of it) is found with the inchon
!    tagged version of this module. it is removed here, since it has
!    not been tested or validated and is considered undesirable in
!    a code being prepared for public release. no immediate need for
!    it is seen at this time, but it will be added back when such need
!    arises.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        call error_mesg ('radiation_driver_mod', &
               ' ability to calculate radiation on subset of columns'//&
              ' and/or levels not yet implemented',  FATAL)

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
                                            SECONDS_PER_DAY
          if (present(mask)) then
            Rad_output%tdtlw(is:ie,js:je,:) =   &
                           (Lw_output%heatra(:,:,:)/SECONDS_PER_DAY)*  &
                            mask(:,:,:)

            Rad_output%tdt_rad (is:ie,js:je,:) =   &
                           (Rad_output%tdtsw(is:ie,js:je,:) + &
                            Lw_output%heatra(:,:,:)/SECONDS_PER_DAY)*  &
                            mask(:,:,:)
          else
            Rad_output%tdtlw(is:ie,js:je,:) =   &
                             Lw_output%heatra(:,:,:)/SECONDS_PER_DAY
            Rad_output%tdt_rad (is:ie,js:je,:) =  &
                                  (Rad_output%tdtsw(is:ie,js:je,:) +   &
                                   Lw_output%heatra(:,:,:)/  &
                                   SECONDS_PER_DAY)
          endif
          if (do_clear_sky_pass) then
            Rad_output%tdtsw_clr(is:ie,js:je,:) =   &
                                          Sw_output%hswcf(:,:,:)/  &
                                          SECONDS_PER_DAY
            if (present(mask)) then
              Rad_output%tdt_rad_clr(is:ie,js:je,:) =    &
                         (Rad_output%tdtsw_clr(is:ie,js:je,:) +  &
                          Lw_output%heatracf(:,:,:)/SECONDS_PER_DAY)* &
                          mask(:,:,:)
            else
              Rad_output%tdt_rad_clr(is:ie,js:je,:) =    &
                           (Rad_output%tdtsw_clr(is:ie,js:je,:) +    &
                            Lw_output%heatracf(:,:,:)/  &
                            SECONDS_PER_DAY)
            endif
          endif

          kmax = size (Rad_output%tdtsw,3)
          Rad_output%flux_sw_surf(is:ie,js:je) =   &
                                         Sw_output%dfsw(:,:,kmax+1) - &
                                         Sw_output%ufsw(:,:,kmax+1)
          Rad_output%flux_lw_surf(is:ie,js:je) =    &
                       STEFAN*Atmos_input%temp(:,:,kmax+1)**4 -   &
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




end subroutine radiation_calc




!######################################################################
! <SUBROUTINE NAME="update_rad_fields">
!  <OVERVIEW>
!    update_rad_fields defines the current radiative heating rate, 
!    surface long and short wave fluxes and cosine of zenith angle
!    to be returned to physics_driver, including renormalization 
!    effects when that option is activated.
!  </OVERVIEW>
!  <DESCRIPTION>
!    update_rad_fields defines the current radiative heating rate, 
!    surface long and short wave fluxes and cosine of zenith angle
!    to be returned to physics_driver, including renormalization 
!    effects when that option is activated.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call update_rad_fields (is, ie, js, je, Time_diag, Astro2,   &
!                              Sw_output, Astro, Rad_output, flux_ratio)
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
!      Time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ] 
!  </IN>
!  <INOUT NAME="Astro" TYPE="astronomy_type">
!   astronomical properties on model grid, usually
!                   valid over radiation timestep on entry, on exit are 
!                   valid over model timestep when renormalizing
!  </INOUT>
!  <IN NAME="Astro2" TYPE="astronomy_type">
!   astronomical properties on model grid, valid over 
!                   physics timestep, used when renormalizing sw fluxes
!  </IN>
!  <INOUT NAME="Rad_output" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </INOUT>
!  <IN NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </IN>
!  <OUT NAME="flux_ratio" TYPE="real">
!   factor to multiply the radiation step values of 
!                   sw fluxes and heating rates by in order to get
!                   current physics timestep values
!  </OUT>
! </SUBROUTINE>
!
subroutine update_rad_fields (is, ie, js, je, Time_diag, Astro2,   &
                              Sw_output, Astro, Rad_output, flux_ratio)

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
real,  dimension(:,:),   intent(out)   ::  flux_ratio

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
!
!  intent(out) variables:
!
!      flux_ratio   factor to multiply the radiation step values of 
!                   sw fluxes and heating rates by in order to get
!                   current physics timestep values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
!
      real, dimension (size(Rad_output%tdt_rad,1),                    &
                       size(Rad_output%tdt_rad,2),                    &
                       size(Rad_output%tdt_rad,3))  ::  tdtlw, tdtlw_clr
      integer   ::   i, j, k

!---------------------------------------------------------------------
!  local variables:
!
!     tdtlw              longwave heating rate
!                        [ deg K sec(-1) ]
!     tdtlw_clr          longwave heating rate under clear sky 
!                        conditions
!                        [ deg K sec(-1) ]
!     i,j,k              do-loop indices
!
!---------------------------------------------------------------------

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
        do k=1, size(Rad_output%tdt_rad,3)
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
          do k=1, size(Rad_output%tdt_rad,3)
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

!--------------------------------------------------------------------



end subroutine update_rad_fields 



!####################################################################
! <SUBROUTINE NAME="produce_radiation_diagnostics">
!  <OVERVIEW>
!    produce_radiation_diagnostics produces netcdf output and global 
!    and hemispheric integrals of radiation package variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    produce_radiation_diagnostics produces netcdf output and global 
!    and hemispheric integrals of radiation package variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call produce_radiation_diagnostics          &
!                            (is, ie, js, je, Time_diag, lat, ts, asfc, &
!                             flux_ratio, Astro, Rad_output, Lw_output, &
!                             Sw_output, Fsrad_output, mask) 
!  </TEMPLATE>
!  <IN NAME="is, ie, js, je" TYPE="integer">
!   starting/ending i,j indices in global storage arrays
!  </IN> 

!  <IN NAME="Time_diag" TYPE="time_type">
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ] 
!  </IN>

!  <IN NAME="lat" TYPE="real">
!    lat        mean latitude (in radians) of all grid boxes processed by this
!               call to radiation_driver   [real, dimension(:,:)]
!  </IN>
!  <IN NAME="ts" TYPE="real">
!   Surface skin temperature
!  </IN>
!  <IN NAME="asfc" TYPE="real">
!   surface albedo
!  </IN>
!  <IN NAME="flux_ratio" TYPE="real">
!   renormalization factor for sw fluxes and heating 
!                   rates
!  </IN>
!  <IN NAME="Astro" TYPE="astronomy_type">
!   astronomical input data for the radiation package
!  </IN>
!  <IN NAME="Rad_output" TYPE="rad_output_type">
!   Radiation output from radiation package, contains variables
!                     which are output from radiation_driver to the 
!                     calling routine, and then used elsewhere within
!                     the component models.
!  </IN>
!  <IN NAME="Lw_output" TYPE="lw_output_type">
!      longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!  </IN>
!  <IN NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
subroutine produce_radiation_diagnostics          &
                            (is, ie, js, je, Time_diag, lat, ts, asfc, &
                             flux_ratio, Astro, Rad_output, Lw_output, &
                             Sw_output, Fsrad_output, mask) 

!--------------------------------------------------------------------
!    produce_radiation_diagnostics produces netcdf output and global 
!    and hemispheric integrals of radiation package variables.
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
!      Astro        astronomical  variables input to the radiation
!                   package [ dimensionless ]
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

      real, dimension (ie-is+1,je-js+1, size(Rad_output%tdtsw,3)) ::  &
                                                tdtlw, tdtlw_clr,&
                                                hsw, hswcf

      real, dimension (ie-is+1,je-js+1, size(Rad_output%tdtsw,3)+1) :: &
                                                dfsw, ufsw,  &
                                                dfswcf, ufswcf,&
                                                fsw, fswcf

      integer           :: j, k
      integer           :: ipass
      logical           :: used
      integer           :: iind, jind
      integer           :: kmax

!---------------------------------------------------------------------
!    if sw flux renormalization is active, modify the fluxes calculated
!    on the last radiation step by the normalization factor based on
!    the difference in solar factor between the current model step and
!    the current radiation step.
!----------------------------------------------------------------------
      kmax = size (Rad_output%tdtsw,3)
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
          used = send_data (id_tdt_sw(ipass),    &
                            Rad_output%tdtsw(is:ie,js:je,:),   &
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
            used = send_data (id_tdt_sw(ipass),   &
                              Rad_output%tdtsw_clr(is:ie,js:je,:),  &
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
      endif   ! (renormalize_sw_fluxes .or. do_rad .or.   
              !  all_step_diagnostics)

!---------------------------------------------------------------------
!    define the longwave diagnostic arrays for the sea-esf radiation 
!    package.  convert to mks units.
!---------------------------------------------------------------------
      if (do_sea_esf_rad) then
        if (do_rad) then
          olr  (:,:)   = Lw_output%flxnet(:,:,1)
          lwups(:,:)   =   STEFAN*ts(:,:  )**4
          lwdns(:,:)   = lwups(:,:) - Lw_output%flxnet(:,:,kmax+1)
          tdtlw(:,:,:) = Lw_output%heatra(:,:,:)/ SECONDS_PER_DAY

          if (do_clear_sky_pass) then
            olr_clr  (:,:)   = Lw_output%flxnetcf(:,:,1)
            lwups_clr(:,:)   =              STEFAN*ts(:,:  )**4
            lwdns_clr(:,:)   = lwups_clr(:,:) -    & 
                               Lw_output%flxnetcf(:,:,kmax+1)
            tdtlw_clr(:,:,:) = Lw_output%heatracf(:,:,:)/SECONDS_PER_DAY
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



end subroutine produce_radiation_diagnostics



!###################################################################
! <SUBROUTINE NAME="deallocate_arrays">
!  <OVERVIEW>
!    deallocate_arrays deallocates the array space of local 
!    derived-type variables.
!  </OVERVIEW>
!  <DESCRIPTION>
!    deallocate_arrays deallocates the array space of local 
!    derived-type variables.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call deallocate_arrays (Cldrad_props, Astro, Astro2, Lw_output, &
!                              Fsrad_output, Sw_output)
!  </TEMPLATE>
!  <INOUT NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   Cloud radiative properties
!  </INOUT>
!  <INOUT NAME="Astro" TYPE="astronomy_type">
!   astronomical data for the radiation package
!  </INOUT>
!  <INOUT NAME="Astro2" TYPE="astronomy_type">
!   astronomical data for the radiation package
!  </INOUT>
!  <INOUT NAME="Fsrad_output" TYPE="rad_output_type">
!   radiation output data from the 
!                        original_fms_rad radiation package, when that 
!                        package is active
!  </INOUT>
!  <INOUT NAME="Lw_output" TYPE="lw_output_type">
!      longwave radiation output data from the 
!                        sea_esf_rad radiation package, when that 
!                        package is active
!  </INOUT>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output data from the 
!                        sea_esf_rad radiation package  when that 
!                        package is active
!  </INOUT>
! </SUBROUTINE>
!
subroutine deallocate_arrays (Cldrad_props, Astro, Astro2, Lw_output, &
                              Fsrad_output, Sw_output)

!---------------------------------------------------------------------
!    deallocate_arrays deallocates the array space of local 
!    derived-type variables.
!---------------------------------------------------------------------

type(cldrad_properties_type), intent(inout)   :: Cldrad_props
type(astronomy_type)        , intent(inout)   :: Astro, Astro2
type(lw_output_type)        , intent(inout)   :: Lw_output
type(fsrad_output_type)     , intent(inout)   :: Fsrad_output
type(sw_output_type)        , intent(inout)   :: Sw_output

!--------------------------------------------------------------------
!    deallocate the variables in Astro and Astro2.
!--------------------------------------------------------------------
      if ( do_rad .or. renormalize_sw_fluxes ) then 
        deallocate (Astro%solar)
        deallocate (Astro%cosz )
        deallocate (Astro%fracday)
        if (Environment%running_gcm) then
          if ( do_rad .and. renormalize_sw_fluxes   &
              .and. Sw_control%do_diurnal ) then 
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
!    call cldrad_props_dealloc to deallocate the variables in 
!    Cldrad_props. 
!--------------------------------------------------------------------
      if (do_rad .and. do_sea_esf_rad) then
        call cldrad_props_dealloc (Cldrad_props)
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
! <SUBROUTINE NAME="calculate_auxiliary_variables">
!  <OVERVIEW>
!    calculate_auxiliary_variables defines values of model delta z and
!    relative humidity, and the values of pressure and temperature at
!    the grid box vertical interfaces.
!  </OVERVIEW>
!  <DESCRIPTION>
!    calculate_auxiliary_variables defines values of model delta z and
!    relative humidity, and the values of pressure and temperature at
!    the grid box vertical interfaces.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call calculate_auxiliary_variables (Atmos_input)
!  </TEMPLATE>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!   atmos_input_type variable, its press and temp
!                   components are input, and its deltaz, rel_hum, 
!                   pflux, tflux and aerosolrelhum components are 
!                   calculated here and output.
!  </INOUT>
! </SUBROUTINE>
!
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
!                   pflux, tflux and aerosolrelhum components are 
!                   calculated here and output.
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables

      real, dimension (size(Atmos_input%temp, 1), &
                       size(Atmos_input%temp, 2), &
                       size(Atmos_input%temp, 3) - 1) :: &
                                                     esat, psat, qv, tv
      integer   ::  k
      integer   ::  kmax


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
!    define deltaz in meters.
!-------------------------------------------------------------------
      tv(:,:,:) = Atmos_input%temp(:,:,ks:ke)*    &
                  (1.0 + D608*Atmos_input%rh2o(:,:,:))
      Atmos_input%deltaz(:,:,ks) = log_p_at_top*RDGAS*tv(:,:,ks)/GRAV
      do k =ks+1,ke   
        Atmos_input%deltaz(:,:,k) = alog(Atmos_input%pflux(:,:,k+1)/  &
                                         Atmos_input%pflux(:,:,k))*   &
                                         RDGAS*tv(:,:,k)/GRAV
      end do

!-------------------------------------------------------------------
!    define deltaz in meters to be used in cloud feedback analysis.
!-------------------------------------------------------------------
      tv(:,:,:) = Atmos_input%cloudtemp(:,:,ks:ke)*    &
                  (1.0 + D608*Atmos_input%cloudvapor(:,:,:))
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
      kmax = size(Atmos_input%temp,3) - 1
      call lookup_es (Atmos_input%temp(:,:,1:kmax), esat)
      do k=1,kmax
        qv(:,:,k) = Atmos_input%rh2o(:,:,k) /    &
                                       (1.0 + Atmos_input%rh2o(:,:,k))
        psat(:,:,k) = Atmos_input%press(:,:,k) - D378*esat(:,:,k)
        psat(:,:,k) = MAX (psat(:,:,k), esat(:,:,k))
        Atmos_input%rel_hum(:,:,k) = qv(:,:,k) / (D622*esat(:,:,k) / &
                                                  psat(:,:,k))
        Atmos_input%rel_hum(:,:,k) =    &
                                  MIN (Atmos_input%rel_hum(:,:,k), 1.0)
      end do

!------------------------------------------------------------------
!    define the relative humidity seen by the aerosol code.
!------------------------------------------------------------------
      call lookup_es (Atmos_input%aerosoltemp(:,:,1:kmax), esat)
      do k=1,kmax
        qv(:,:,k) = Atmos_input%aerosolvapor(:,:,k) /    &
                                (1.0 + Atmos_input%aerosolvapor(:,:,k))
        psat(:,:,k) = Atmos_input%aerosolpress(:,:,k) - D378*esat(:,:,k)
        psat(:,:,k) = MAX (psat(:,:,k), esat(:,:,k))
        Atmos_input%aerosolrelhum(:,:,k) = qv(:,:,k) /   &
                                        (D622*esat(:,:,k) / psat(:,:,k))
         Atmos_input%aerosolrelhum(:,:,k) =    &
                            MIN (Atmos_input%aerosolrelhum(:,:,k), 1.0)
      end do
 
!----------------------------------------------------------------------


end subroutine calculate_auxiliary_variables



!#######################################################################





                 end module radiation_driver_mod

