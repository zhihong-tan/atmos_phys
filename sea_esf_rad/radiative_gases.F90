                  module radiative_gases_mod


use time_manager_mod,    only:   &
!                              time_manager_init, &
                               time_type
use rad_utilities_mod,   only: longwave_control_type, Lw_control, &
                               radiation_control_type, Rad_control,&
                               environment_type, Environment, &
                               atmos_input_type, &
                               rad_utilities_init, &
                               longwave_parameter_type, Lw_parameters, &
                               radiative_gases_type
use utilities_mod,       only: open_file, file_exist,    &
                               check_nml_error, error_mesg, &
                               utilities_init, &
                               print_version_number, FATAL, NOTE, &
!                              get_root_pe, &
                               WARNING, get_my_pe, close_file 

!    modules for individual radiative_gases:

use ozone_mod,           only: ozone_driver, ozone_init, ozone_end


implicit none
private

!---------------------------------------------------------------------
!    radiative_gases_mod defines mixing ratios of radiatively-active 
!    gases to be used in the calculation of longwave and shortwave
!    radiative fluxes and heating rates in the sea_esf_rad radiation
!    package.
!---------------------------------------------------------------------
!

!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: radiative_gases.F90,v 1.3 2002/07/16 22:36:32 fms Exp $'
character(len=128)  :: tag     =  '$Name: inchon $'

!---------------------------------------------------------------------
!-------  interfaces --------

public   radiative_gases_init, define_radiative_gases,  &
         radiative_gases_end


private  read_restart_radiative_gases, &
         define_ch4, define_n2o, define_co2, define_f11, &
         define_f12, define_f113, define_f22, &
         write_restart_radiative_gases


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)    :: ch4_data_source  = '   '
character(len=16)    :: n2o_data_source  = '   '
character(len=16)    :: f11_data_source  = '   '
character(len=16)    :: f12_data_source  = '   '
character(len=16)    :: f113_data_source = '   '
character(len=16)    :: f22_data_source  = '   '
character(len=16)    :: co2_data_source  = '   '

 

namelist /radiative_gases_nml/                    &
                                ch4_data_source,  &
                                n2o_data_source,  &
                                f11_data_source,  &
                                f12_data_source,  &
                                f113_data_source, &
                                f22_data_source,  &
                                co2_data_source


!---------------------------------------------------------------------
!------- public data ------



!---------------------------------------------------------------------
!------- private data ------

!---------------------------------------------------------------------
!    list of restart versions of radiation_driver.res readable by this 
!    module.
!---------------------------------------------------------------------
integer, dimension(2)    ::  restart_versions = (/ 1, 2 /)


!--------------------------------------------------------------------
!    initial mixing ratios of the various radiative gases. if a gas 
!    is not active, its mixing ratio is set to zero.
!--------------------------------------------------------------------
real         ::  rch4, rn2o, rf11, rf12, rf113, rf22, rco2

!--------------------------------------------------------------------
!    these variables contain the mixing ratios of the radiative gases
!    at the current time. if the gases are fixed, this is the same as
!    the initial value; otherwise it will be time-varying in either
!    a prescribed way, or as predicted from within the model.
!--------------------------------------------------------------------
real         :: rrvco2, rrvf11, rrvf12, rrvf113, rrvf22, rrvch4, rrvn2o 

!--------------------------------------------------------------------
!    logical variables indicating whether the gas is allowed to vary
!    in time.
!--------------------------------------------------------------------
logical      :: time_varying_ch4  = .false.,    &
                time_varying_n2o  = .false.,    &
                time_varying_co2  = .false.,    &
                time_varying_f11  = .false.,    &
                time_varying_f12  = .false.,    &
                time_varying_f113 = .false.,    &
                time_varying_f22  = .false.

!---------------------------------------------------------------------
logical      :: restart_present            =    &
                                        .false. ! restart file present ?
logical      :: radiative_gases_initialized=    &
                                        .false. ! module initialized ?



!---------------------------------------------------------------------
!---------------------------------------------------------------------



                          contains




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PUBLIC SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!####################################################################

subroutine radiative_gases_init (pref, latb, lonb)

!---------------------------------------------------------------------
!    radiative_gases_init is the constructor for radiative_gases_mod.
!---------------------------------------------------------------------

real, dimension(:,:), intent(in) :: pref
real, dimension(:),   intent(in) :: latb, lonb

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 for use in defining transmission functions [pascals]
!       latb      array of model latitudes at cell boundaries [radians]
!       lonb      array of model longitudes at cell boundaries [radians]
!
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer              :: unit, ierr, io


!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (radiative_gases_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call utilities_init
      call rad_utilities_init
!! routine not yet existent:
!     call time_manager_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
          read (unit, nml=radiative_gases_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'radiative_gases_nml')
        enddo
10      call close_file (unit)
      endif

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      unit = open_file ('logfile.out', action='append')
      if ( get_my_pe() == 0 ) then
!     if ( get_my_pe() == get_root_pe() ) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit,nml=radiative_gases_nml)
      endif
      call close_file (unit)

!---------------------------------------------------------------------
!    enforce necessary conditions on data source for ch4 and n2o. in
!    the future this constraint may be removed.
!----------------------------------------------------------------------
      if (trim(ch4_data_source) /= trim(n2o_data_source)) then
        call error_mesg ( 'radiative_gases_mod', &
          ' currently ch4 and n2o must have the same data source', &
                                                                 FATAL)
      endif

!---------------------------------------------------------------------
!   define the number of individual bands in the 1200 - 1400 cm-1 
!   region when ch4 and n2o are active.
!---------------------------------------------------------------------
      if (trim(ch4_data_source) /= 'absent' .or. &
          trim(n2o_data_source) /= 'absent') then
        Lw_control%do_ch4_n2o = .true.
      else
        Lw_control%do_ch4_n2o = .false.
      endif

!---------------------------------------------------------------------
!   if any of the cfcs are activated, set a flag indicating that cfcs
!   are active.
!---------------------------------------------------------------------
      if (trim(f11_data_source) /= 'absent' .or. &
          trim(f12_data_source) /= 'absent' .or. &
          trim(f113_data_source) /= 'absent' .or. &
          trim(f22_data_source) /= 'absent' ) then 
        Lw_control%do_cfc = .true.
      else
        Lw_control%do_cfc = .false.
      endif

!---------------------------------------------------------------------
!   define a logical variable indicating whether co2 is to be activated.
!   currently co2 must be activated. 
!---------------------------------------------------------------------
      if (trim(co2_data_source) == 'absent') then
        call error_mesg ('radiative_gases_mod', &
         ' currently MUST include co2 as a radiative gas', FATAL)
        Lw_control%do_co2 = .false.
      else
        Lw_control%do_co2 = .true.
      endif

!---------------------------------------------------------------------
!   if present, read the radiative gases restart file. set a flag 
!   indicating the presence of the file.
!---------------------------------------------------------------------
      if (file_exist ('INPUT/radiative_gases.res')) then
        call read_restart_radiative_gases
        restart_present = .true.
      endif

!---------------------------------------------------------------------
!   call a routine for each gas to initialize its mixing ratio
!   and set a flag indicating whether it is fixed in time or time-
!   varying. the values of time-varying gases will come from a restart 
!   file (if present) or an input file, REGARDLESS OF WHAT IS SPECIFIED 
!   IN THE NAMELIST. fixed-in-time gases will be defined from the source
!   specified in the namelist.
!---------------------------------------------------------------------
      call define_ch4 (ch4_data_source)
      call define_n2o (n2o_data_source)
      call define_f11 (f11_data_source)
      call define_f12 (f12_data_source)
      call define_f113(f113_data_source)
      call define_f22 (f22_data_source)
      call define_co2 (co2_data_source)

!--------------------------------------------------------------------
!    define variable which will be contain mixing ratio of each gas
!    as model is integrated in time.
!--------------------------------------------------------------------
      rrvch4  = rch4
      rrvn2o  = rn2o
      rrvf11  = rf11
      rrvf12  = rf12
      rrvf113 = rf113
      rrvf22  = rf22
      rrvco2  = rco2

!--------------------------------------------------------------------- 
!   call ozone_init to initialize the ozone field.
!---------------------------------------------------------------------
     call ozone_init (pref, latb, lonb)

!---------------------------------------------------------------------
     radiative_gases_initialized = .true.


end subroutine radiative_gases_init




!##################################################################

 subroutine define_radiative_gases (is, ie, js, je, Rad_time, lat, &
                                    Atmos_input, Time_next, Rad_gases)

!-------------------------------------------------------------------
!    define_radiative_gases returns the current values of the radiative 
!    gas mixing ratios to radiation_driver in Rad_gases.
!-------------------------------------------------------------------

integer,                    intent(in)  :: is, ie, js, je
type(time_type),            intent(in)  :: Rad_time, Time_next
real, dimension(:,:),       intent(in)  :: lat
type(atmos_input_type),     intent(in)  :: Atmos_input
type(radiative_gases_type), intent(out) :: Rad_gases

!---------------------------------------------------------------------
!
!  intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Rad_Time     time at which radiation is to be calculated
!                   [ time_type (days, seconds) ] 
!      Time_next    time on next timestep, used as stamp for diagnostic 
!                   output  [ time_type  (days, seconds) ]  
!      lat          latitude of model points  [ radians ]
!      Atmos_input  atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!
!  intent(out) variables:
!
!      Rad_gases    radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!!! NOTE :
!
!    this is where the time-variation of radiative gases will be added.
!    they may either be prognostically predicted or prescribed, as 
!    indicated in the xxx_data_source namelist variable. in either case 
!    a subroutine call can be made here, which will return the "new" 
!    value of radiative gas mixing ratio.
!
!    if prescribed, how the prescription of the gas value is made is 
!    indeterminate at this time -- the necessary mods may be made to 
!    accomodate whatever technique or techniques are selected, be they 
!    a linear trend, interpolation into a time series, or any other 
!    method.
!
!    if predicted, a call to the module predicting it (and storing it
!    as a module variable) may be made, or if available it could be 
!    passed in through an argument to radiation_driver.
!!! END NOTE 
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    call routine to return values of rrvch4 and rrvn2o at current time.
!----------------------------------------------------------------------
      if (time_varying_ch4 .or. time_varying_n2o ) then
        call error_mesg ( 'radiative_gases_mod', &
            'time-varying ch4 and / or n2o not yet implemented', FATAL)
      endif

!--------------------------------------------------------------------
!    call routine to return values of cfcs at current time.
!----------------------------------------------------------------------
      if (time_varying_f11 .or. &
          time_varying_f12 .or. &
          time_varying_f113 .or. &
          time_varying_f22 ) then 
        call error_mesg ( 'radiative_gases_mod', &
                 'time-varying cfcs not yet implemented', FATAL)
      endif

!--------------------------------------------------------------------
!    call routine to return values of co2 at current time.
!----------------------------------------------------------------------
      if (time_varying_co2) then 
        call error_mesg ( 'radiative_gases_mod', &
                  'time-varying co2 not yet implemented', FATAL)
      endif

!--------------------------------------------------------------------
!    fill the contents of the radiative_gases_type variable which
!    will be returned to the calling routine. if a gas is time-varying,
!    its newly defined value will be returned. if it is fixed in time,
!    the value defined in the _init routine will be returned.
!---------------------------------------------------------------------
      Rad_gases%rrvch4  = rrvch4
      Rad_gases%rrvn2o  = rrvn2o
      Rad_gases%rrvf11  = rrvf11
      Rad_gases%rrvf12  = rrvf12
      Rad_gases%rrvf113 = rrvf113
      Rad_gases%rrvf22  = rrvf22
      Rad_gases%rrvco2  = rrvco2
      Rad_gases%time_varying_co2  = time_varying_co2
      Rad_gases%time_varying_ch4  = time_varying_ch4
      Rad_gases%time_varying_n2o  = time_varying_n2o
      Rad_gases%time_varying_f11  = time_varying_f11
      Rad_gases%time_varying_f12  = time_varying_f12
      Rad_gases%time_varying_f113 = time_varying_f113
      Rad_gases%time_varying_f22  = time_varying_f22

!--------------------------------------------------------------------
!    call ozone_driver to define the ozone field to be used for the 
!    radiation calculation.
!--------------------------------------------------------------------
      call ozone_driver (is, ie, js, je, lat, Rad_time, Atmos_input, &
                         Rad_gases)


end subroutine define_radiative_gases


!####################################################################


subroutine radiative_gases_end

!---------------------------------------------------------------------
!    radiative_gases_end is the destructor for radiative_gases_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    write out the radiative_gases restart file.
!---------------------------------------------------------------------
      if (Environment%running_gcm) then 
        call write_restart_radiative_gases
      endif

!---------------------------------------------------------------------
!    call the destructors for the component gas module(s).
!---------------------------------------------------------------------
      call ozone_end

!--------------------------------------------------------------------
      radiative_gases_initialized = .false.


end subroutine radiative_gases_end


!####################################################################
 



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





subroutine read_restart_radiative_gases 

!---------------------------------------------------------------------
!    read_restart_radiative_gases reads the radiative_gases.res file.
!---------------------------------------------------------------------

      integer  ::  vers   ! version number of restart file 
      integer  ::  unit

!--------------------------------------------------------------------
!    determine if  a radiative_gases.parameters.res file is present.
!    this file is only present in restart version 1.
!--------------------------------------------------------------------
      if (file_exist('INPUT/radiative_gases.parameters.res' ) ) then

!---------------------------------------------------------------------
!    read radiative gas restart file, version 1.
!---------------------------------------------------------------------
        if (file_exist('INPUT/radiative_gases.res' ) ) then
          unit = open_file ('INPUT/radiative_gases.res',    &
                            form= 'unformatted', action= 'read')
          read (unit) rco2
          read (unit) rf11, rf12, rf113, rf22
          read (unit) rch4, rn2o
          call close_file (unit)
        endif ! (file_exist(.res))

      else 
!---------------------------------------------------------------------
!    read radiative gas restart file. version number will be the first
!    file record.
!---------------------------------------------------------------------
        if (file_exist('INPUT/radiative_gases.res' ) ) then
          unit = open_file ('INPUT/radiative_gases.res',   &
	                    form= 'unformatted',  action= 'read')
	  read (unit) vers

!---------------------------------------------------------------------
!    verify that this restart file version is readable by the current
!    code. if not, print a message.
!---------------------------------------------------------------------
          if ( .not. any(vers == restart_versions) ) then
            call error_mesg ('radiative_gases_mod', &
                      'radiative_gases restart problem --  may be &
		      &attempting to read version 1 file &
		      &w/o parameters.res file being present.',  FATAL)
          endif
          read (unit) rco2
          read (unit) rf11, rf12, rf113, rf22
          read (unit) rch4, rn2o
          call close_file (unit)
        endif 
      endif 

!--------------------------------------------------------------------


end subroutine read_restart_radiative_gases


!###################################################################

subroutine define_ch4 (data_source) 

!---------------------------------------------------------------------
!    define_ch4 provides initial values for ch4 mixing ratio. if ch4
!    is fixed in time, the value is given by the namelist specification.
!    if ch4 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define ch4 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rch4_ipcc_80  = 1.56900E-06
      real       ::  rch4_ipcc_92  = 1.71400E-06
      real       ::  rch4_icrccm   = 1.75000E-06
      real       ::  rch4_ipcc_98  = 1.82120E-06

      integer    :: unit, ierr, io, inrad


!---------------------------------------------------------------------
!    define initial ch4 mixing ratios to be used.
!    'icrccm'     --> rch4_icrccm
!    'ipcc80'     --> rch4_ipcc_80
!    'ipcc92'     --> rch4_ipcc_92
!    'ipcc98'     --> rch4_ipcc_98
!    'input'      --> file INPUT/id1ch4n2o, record 1
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!    'absent'     --> 0.0; ch4 not active in this experiment
!---------------------------------------------------------------------
      if (trim(data_source)      == 'icrccm') then
        rch4   = rch4_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rch4   = rch4_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rch4   = rch4_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rch4   = rch4_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (file_exist ('INPUT/id1ch4n2o') ) then
          inrad = open_file ('INPUT/id1ch4n2o', form= 'formatted', &
                             action= 'read')
          read (inrad, FMT = '(5e18.10)')  rch4
          call close_file (inrad)
        else
          call error_mesg ( 'radiative_gases_mod', &
                   'desired ch4_n2o input file is not present', FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ('radiative_gases_init', &
           'cannot use restart ch4 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1ch4n2o') ) then
            inrad = open_file ('INPUT/id1ch4n2o', form= 'formatted', &
                             action= 'read')
            read (inrad, FMT = '(5e18.10)') rch4
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor ch4 input file is present. one &
	          &of these is required when the data_source is  &
	          &"prescribed" ', FATAL)
          endif
        endif
        time_varying_ch4 = .true.

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1ch4n2o') ) then
            inrad = open_file ('INPUT/id1ch4n2o', form= 'formatted', &
                               action= 'read')
            read (inrad, FMT = '(5e18.10)') rch4
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor ch4 input file is present. one &
	          &of these is required when the data_source is  &
	          &"predicted" ', FATAL)
          endif
        endif
        time_varying_ch4 = .true.

!--------------------------------------------------------------------
!    if data_source is 'absent', set the mixing ratio to 0.0. ch4 is
!    inactive in the current experiment.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'absent') then
         rch4 = 0.0

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for ch4' , FATAL)
      endif

!---------------------------------------------------------------------


end subroutine define_ch4



!#####################################################################

subroutine define_n2o (data_source) 

!---------------------------------------------------------------------
!    define_n2o provides initial values for n2o mixing ratio. if n2o
!    is fixed in time, the value is given by the namelist specification.
!    if n2o is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define n2o initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rn2o_ipcc_80  = 3.02620E-07
      real       ::  rn2o_ipcc_92  = 3.11000E-07
      real       ::  rn2o_icrccm   = 2.80000E-07
      real       ::  rn2o_ipcc_98  = 3.16000E-07

      integer    :: unit, ierr, io, inrad


!---------------------------------------------------------------------
!    define initial n2o mixing ratios to be used.
!    'icrccm'     --> rn2o_icrccm
!    'ipcc80'     --> rn2o_ipcc_80
!    'ipcc92'     --> rn2o_ipcc_92
!    'ipcc98'     --> rn2o_ipcc_98
!    'input'      --> file INPUT/id1ch4n2o, record 2
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!    'absent'     --> 0.0; n2o not active in this experiment
!---------------------------------------------------------------------
      if (trim(data_source)      == 'icrccm') then
        rn2o   = rn2o_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rn2o   = rn2o_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rn2o   = rn2o_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rn2o   = rn2o_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (file_exist ('INPUT/id1ch4n2o') ) then
          inrad = open_file ('INPUT/id1ch4n2o', form= 'formatted', &
                             action= 'read')
          read (inrad, FMT = '(5e18.10)')  
          read (inrad, FMT = '(5e18.10)')  rn2o
          call close_file (inrad)
        else
          call error_mesg ( 'radiative_gases_mod', &
                   'desired ch4_n2o input file is not present', FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ('radiative_gases_init', &
           'cannot use restart n2o values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1ch4n2o') ) then
            inrad = open_file ('INPUT/id1ch4n2o', form= 'formatted', &
                             action= 'read')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rn2o
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor n2o input file is present. one &
	          &of these is required when the data_source is  &
	          &"prescribed" ', FATAL)
          endif
        endif
        time_varying_n2o = .true.

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1ch4n2o') ) then
            inrad = open_file ('INPUT/id1ch4n2o', form= 'formatted', &
                               action= 'read')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rn2o
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor n2o input file is present. one &
	          &of these is required when the data_source is  &
	          &"predicted" ', FATAL)
          endif
        endif
        time_varying_n2o = .true.

!--------------------------------------------------------------------
!    if data_source is 'absent', set the mixing ratio to 0.0. n2o is
!    inactive in the current experiment.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'absent') then
         rn2o = 0.0

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for n2o ', FATAL)
      endif

!---------------------------------------------------------------------


end subroutine define_n2o



!#####################################################################

subroutine define_f11 (data_source) 

!---------------------------------------------------------------------
!    define_f11 provides initial values for f11 mixing ratio. if f11
!    is fixed in time, the value is given by the namelist specification.
!    if f11 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define f11 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rf11_icrccm   = 1.00000E-09
      real       ::  rf11_ipcc_80  = 1.57500E-10
      real       ::  rf11_ipcc_92  = 2.68000E-10
      real       ::  rf11_ipcc_98  = 2.68960E-10

      integer    :: unit, ierr, io, inrad


!---------------------------------------------------------------------
!    define initial f11 mixing ratios to be used.
!    'icrccm'     --> rf11_icrccm
!    'ipcc80'     --> rf11_ipcc_80
!    'ipcc92'     --> rf11_ipcc_92
!    'ipcc98'     --> rf11_ipcc_98
!    'input'      --> file INPUT/id1cfc, record 1
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!    'absent'     --> 0.0; f11 not active in this experiment
!---------------------------------------------------------------------
      if (trim(data_source)      == 'icrccm') then
        rf11   = rf11_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rf11   = rf11_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rf11   = rf11_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rf11   = rf11_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (file_exist ('INPUT/id1cfc') ) then
          inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                             action= 'read')
          read (inrad, FMT = '(5e18.10)')  rf11
          call close_file (inrad)
        else
          call error_mesg ( 'radiative_gases_mod', &
                   'desired cfc input file is not present', FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ('radiative_gases_init', &
           'cannot use restart f11 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                             action= 'read')
            read (inrad, FMT = '(5e18.10)') rf11
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor f11 input file is present. one &
	          &of these is required when the data_source is  &
	          &"prescribed" ', FATAL)
          endif
        endif
        time_varying_f11 = .true.

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                               action= 'read')
            read (inrad, FMT = '(5e18.10)') rf11
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor f11 input file is present. one &
	          &of these is required when the data_source is  &
	          &"predicted" ', FATAL)
          endif
        endif
        time_varying_f11 = .true.

!--------------------------------------------------------------------
!    if data_source is 'absent', set the mixing ratio to 0.0. f11 is
!    inactive in the current experiment.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'absent') then
         rf11 = 0.0

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for f11 ', FATAL)
      endif

!---------------------------------------------------------------------


end subroutine define_f11




!#####################################################################

subroutine define_f12 (data_source) 

!---------------------------------------------------------------------
!    define_f12 provides initial values for f12 mixing ratio. if f12
!    is fixed in time, the value is given by the namelist specification.
!    if f12 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define f12 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rf12_icrccm   = 1.00000E-09
      real       ::  rf12_ipcc_80  = 2.72500E-10
      real       ::  rf12_ipcc_92  = 5.03000E-10
      real       ::  rf12_ipcc_98  = 5.31510E-10

      integer    :: unit, ierr, io, inrad


!---------------------------------------------------------------------
!    define initial f12 mixing ratios to be used.
!    'icrccm'     --> rf12_icrccm
!    'ipcc80'     --> rf12_ipcc_80
!    'ipcc92'     --> rf12_ipcc_92
!    'ipcc98'     --> rf12_ipcc_98
!    'input'      --> file INPUT/id1cfc, record 2
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!    'absent'     --> 0.0; f12 not active in this experiment
!---------------------------------------------------------------------
      if (trim(data_source)      == 'icrccm') then
        rf12   = rf12_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rf12   = rf12_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rf12   = rf12_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rf12   = rf12_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (file_exist ('INPUT/id1cfc') ) then
          inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                             action= 'read')
          read (inrad, FMT = '(5e18.10)')  
          read (inrad, FMT = '(5e18.10)')  rf12
          call close_file (inrad)
        else
          call error_mesg ( 'radiative_gases_mod', &
                   'desired cfc input file is not present', FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ('radiative_gases_init', &
           'cannot use restart f12 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                             action= 'read')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf12
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor f12 input file is present. one &
	          &of these is required when the data_source is  &
	          &"prescribed" ', FATAL)
          endif
        endif
        time_varying_f12 = .true.

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                               action= 'read')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf12
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor f12 input file is present. one &
	          &of these is required when the data_source is  &
	          &"predicted" ', FATAL)
          endif
        endif
        time_varying_f12 = .true.

!--------------------------------------------------------------------
!    if data_source is 'absent', set the mixing ratio to 0.0. f12 is
!    inactive in the current experiment.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'absent') then
         rf12 = 0.0

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for f12', FATAL)
      endif

!---------------------------------------------------------------------


end subroutine define_f12




!#####################################################################

subroutine define_f113 (data_source) 

!---------------------------------------------------------------------
!    define_f113 provides initial values for f113 mixing ratio. if f113
!    is fixed in time, the value is given by the namelist specification.
!    if f113 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define f113 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rf113_icrccm  = 1.00000E-09
      real       ::  rf113_ipcc_80 = 2.31400E-11
      real       ::  rf113_ipcc_92 = 8.20000E-11
      real       ::  rf113_ipcc_98 = 8.58100E-11

      integer    :: unit, ierr, io, inrad


!---------------------------------------------------------------------
!    define initial f113 mixing ratios to be used.
!    'icrccm'     --> rf113_icrccm
!    'ipcc80'     --> rf113_ipcc_80
!    'ipcc92'     --> rf113_ipcc_92
!    'ipcc98'     --> rf113_ipcc_98
!    'input'      --> file INPUT/id1cfc, record 3
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!    'absent'     --> 0.0; f113 not active in this experiment
!---------------------------------------------------------------------
      if (trim(data_source)      == 'icrccm') then
        rf113   = rf113_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rf113   = rf113_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rf113   = rf113_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rf113   = rf113_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (file_exist ('INPUT/id1cfc') ) then
          inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                             action= 'read')
          read (inrad, FMT = '(5e18.10)')  
          read (inrad, FMT = '(5e18.10)')  
          read (inrad, FMT = '(5e18.10)')  rf113
          call close_file (inrad)
        else
          call error_mesg ( 'radiative_gases_mod', &
                   'desired f113 input file is not present', FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ('radiative_gases_init', &
           'cannot use restart f113 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_file ('INPUT/id1c2o', form= 'formatted', &
                             action= 'read')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf113
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor f113 input file is present. one &
	          &of these is required when the data_source is  &
	          &"prescribed" ', FATAL)
          endif
        endif
        time_varying_f113 = .true.

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                               action= 'read')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf113
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor f113 input file is present. one &
	          &of these is required when the data_source is  &
	          &"predicted" ', FATAL)
          endif
        endif
        time_varying_f113 = .true.

!--------------------------------------------------------------------
!    if data_source is 'absent', set the mixing ratio to 0.0. f113 is
!    inactive in the current experiment.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'absent') then
         rf113 = 0.0

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for f113', FATAL)
      endif

!---------------------------------------------------------------------


end subroutine define_f113



!#####################################################################

subroutine define_f22 (data_source) 

!---------------------------------------------------------------------
!    define_f22 provides initial values for f22 mixing ratio. if f22
!    is fixed in time, the value is given by the namelist specification.
!    if f22 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define f22 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rf22_icrccm   = 1.00000E-09
      real       ::  rf22_ipcc_80  = 6.20200E-11
      real       ::  rf22_ipcc_92  = 1.05000E-10
      real       ::  rf22_ipcc_98  = 1.26520E-10

      integer    :: unit, ierr, io, inrad


!---------------------------------------------------------------------
!    define initial f22 mixing ratios to be used.
!    'icrccm'     --> rf22_icrccm
!    'ipcc80'     --> rf22_ipcc_80
!    'ipcc92'     --> rf22_ipcc_92
!    'ipcc98'     --> rf22_ipcc_98
!    'input'      --> file INPUT/id1cfc, record 4
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!    'absent'     --> 0.0; f22 not active in this experiment
!---------------------------------------------------------------------
      if (trim(data_source)      == 'icrccm') then
        rf22   = rf22_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rf22   = rf22_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rf22   = rf22_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rf22   = rf22_ipcc_98

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (file_exist ('INPUT/id1cfc') ) then
          inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                             action= 'read')
          read (inrad, FMT = '(5e18.10)')  
          read (inrad, FMT = '(5e18.10)')  
          read (inrad, FMT = '(5e18.10)')  
          read (inrad, FMT = '(5e18.10)')  rf22
          call close_file (inrad)
        else
          call error_mesg ( 'radiative_gases_mod', &
                   'desired cfc input file is not present', FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ('radiative_gases_init', &
           'cannot use restart f22 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                             action= 'read')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf22
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor f22 input file is present. one &
	          &of these is required when the data_source is  &
	          &"prescribed" ', FATAL)
          endif
        endif
        time_varying_f22 = .true.

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1cfc') ) then
            inrad = open_file ('INPUT/id1c2o', form= 'formatted', &
                               action= 'read')
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') 
            read (inrad, FMT = '(5e18.10)') rf22
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor f22 input file is present. one &
	          &of these is required when the data_source is  &
	          &"predicted" ', FATAL)
          endif
        endif
        time_varying_f22 = .true.

!--------------------------------------------------------------------
!    if data_source is 'absent', set the mixing ratio to 0.0. f22 is
!    inactive in the current experiment.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'absent') then
         rf22 = 0.0

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for f22', FATAL)
      endif

!---------------------------------------------------------------------


end subroutine define_f22



!#####################################################################

subroutine define_co2 (data_source) 

!---------------------------------------------------------------------
!    define_co2 provides initial values for co2 mixing ratio. if co2
!    is fixed in time, the value is given by the namelist specification.
!    if co2 is time-varying, the values are obtained from either a
!    restart or input data file.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       data_source     character string defining source to use to
!                       define co2 initial values
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!---------------------------------------------------------------------
!    initial trace gas volume mixing ratios in (no./no.) from various
!    sources.
!---------------------------------------------------------------------
      real       ::  rco2_icrccm   = 3.00000E-04
      real       ::  rco2_ipcc_92  = 3.56000E-04
      real       ::  rco2_ipcc_80  = 3.37320E-04
      real       ::  rco2_ipcc_98  = 3.69400E-04
      real       ::  rco2_330ppm   = 3.30000E-04
      real       ::  rco2_660ppm   = 6.60000E-04

      integer    :: unit, ierr, io, inrad


!---------------------------------------------------------------------
!    define initial co2 mixing ratios to be used.
!    'icrccm'     --> rco2_icrccm
!    'ipcc80'     --> rco2_ipcc_80
!    'ipcc92'     --> rco2_ipcc_92
!    'ipcc98'     --> rco2_ipcc_98
!    '330ppm'     --> rco2_330ppm 
!    '660ppm'     --> rco2_660ppm 
!    'input'      --> file INPUT/id1co2, record 2
!    'restart'    --> values read from restart file
!    'prescribed' --> from restart file; if restart not present, 
!                     from input file
!    'predicted'  --> from restart file; if restart not present, 
!                     from input file
!    'absent'     --> 0.0; co2 not active in this experiment
!---------------------------------------------------------------------
      if (trim(data_source)      == 'icrccm') then
        rco2   = rco2_icrccm

      else if (trim(data_source) == 'ipcc_80') then
        rco2   = rco2_ipcc_80

      else if (trim(data_source) == 'ipcc_92') then
        rco2   = rco2_ipcc_92  

      else if (trim(data_source) == 'ipcc_98') then
        rco2   = rco2_ipcc_98

      else if (trim(data_source) == '330ppm') then
        rco2   = rco2_330ppm

      else if (trim(data_source) == '660ppm') then
        rco2   = rco2_660ppm

!--------------------------------------------------------------------
!    if data_source is an input file, determine if it is present. if so,
!    open and read.  if not present, write an error message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'input') then
        if (file_exist ('INPUT/id1co2') ) then
          inrad = open_file ('INPUT/id1co2', form= 'formatted', &
                             action= 'read')
          read (inrad, FMT = '(5e18.10)')  rco2
          call close_file (inrad)
        else
          call error_mesg ( 'radiative_gases_mod', &
                   'desired co2 input file is not present', FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is a restart file and it is present, the value to be
!    used has been previously read. if not present, write an error 
!    message and stop.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'restart') then
        if (.not. restart_present ) then
          call error_mesg ('radiative_gases_init', &
           'cannot use restart co2 values without a restart file', &
                                  FATAL)
        endif

!--------------------------------------------------------------------
!    if data_source is 'prescribed' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'prescribed') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1co2') ) then
            inrad = open_file ('INPUT/id1co2', form= 'formatted', &
                             action= 'read')
            read (inrad, FMT = '(5e18.10)') rco2
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor co2 input file is present. one &
	          &of these is required when the data_source is  &
	          &"prescribed" ', FATAL)
          endif
        endif
        time_varying_co2 = .true.

!--------------------------------------------------------------------
!    if data_source is 'predicted' and a restart file is present, the 
!    value to be used has been previously read. if a restart file is
!    not present, check for an input file. if it is present, read the
!    file; if not, write an error message and stop. set the time_vary-
!    ing flag to .true.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'predicted') then
        if (.not. restart_present) then
          if (file_exist ('INPUT/id1co2') ) then
            inrad = open_file ('INPUT/id1co2', form= 'formatted', &
                               action= 'read')
            read (inrad, FMT = '(5e18.10)') rco2
            call close_file (inrad)
          else
            call error_mesg ( 'radiative_gases_init', &
              'neither restart nor co2 input file is present. one &
	          &of these is required when the data_source is  &
	          &"predicted" ', FATAL)
          endif
        endif
        time_varying_co2 = .true.

!--------------------------------------------------------------------
!    if data_source is 'absent', write an error message and stop. co2
!    currently must be activated.
!--------------------------------------------------------------------
      else if (trim(data_source) == 'absent') then
        call error_mesg ('radiative_gases_mod', &
             ' currently MUST include co2 as a radiative gas', FATAL)

!--------------------------------------------------------------------
!    write an error message if the data_source is invalid.
!--------------------------------------------------------------------
      else
        call error_mesg ('radiative_gases_mod', &
         'no valid data source was specified for co2 input', FATAL)
      endif

!---------------------------------------------------------------------


end subroutine define_co2



!####################################################################

subroutine write_restart_radiative_gases

!---------------------------------------------------------------------
!    write_restart_radiative_gases writes the radiative_gases.res file.
!---------------------------------------------------------------------

      integer                 ::  ierr, io, unit

!---------------------------------------------------------------------
!    open unit and write radiative gas restart file.
!---------------------------------------------------------------------
      if (get_my_pe() == 0) then
        unit = open_file ('RESTART/radiative_gases.res',   &
                          form= 'unformatted', action= 'write')
        write (unit) restart_versions(size(restart_versions))
        write (unit) rrvco2
        write (unit) rrvf11, rrvf12, rrvf113, rrvf22
        write (unit) rrvch4, rrvn2o
        call close_file (unit)
      endif

!----------------------------------------------------------------------
      


end subroutine write_restart_radiative_gases


!####################################################################


		  end module radiative_gases_mod

