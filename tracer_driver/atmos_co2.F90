module atmos_co2_mod
! <CONTACT EMAIL="rdslater@splash.princeton.edu">
!   Richard D. Slater
! </CONTACT>

! <REVIEWER EMAIL="none@nowhere.dot">
!   none yet
! </REVIEWER>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
! </OVERVIEW>

! <DESCRIPTION>
! </DESCRIPTION>


use              fms_mod, only : stdlog, stdout, write_version_number
use   tracer_manager_mod, only : get_tracer_index, tracer_manager_init
use    field_manager_mod, only : MODEL_ATMOS

implicit none

private

!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_co2_sourcesink
public  atmos_co2_gather_data
public  atmos_co2_flux_init
public  atmos_co2_init
public  atmos_co2_end

integer, save   :: ind_co2_flux = 0
integer, save   :: ind_co2  = 0

character(len=48), parameter    :: mod_name = 'atmos_co2_mod'

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------
!
!  When initializing additional tracers, the user needs to make the
!  following changes.
!
!  Add an integer variable below for each additional tracer. This should
!  be initialized to zero. 
!
!  Add id_tracername for each additional tracer. These are used in
!  initializing and outputting the tracer fields.
!
!-----------------------------------------------------------------------

! tracer numbers for CO2

logical :: module_is_initialized = .FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$$'
character(len=128) :: tagname = '$$'
!-----------------------------------------------------------------------

contains

!#######################################################################

!<SUBROUTINE NAME ="atmos_co2_sourcesink">
!<OVERVIEW>
!  A subroutine to calculate the internal sources and sinks of carbon dioxide.
!</OVERVIEW>
!

subroutine atmos_co2_sourcesink()

implicit none

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_sourcesink'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

return

end subroutine atmos_co2_sourcesink
!</SUBROUTINE >

!#######################################################################

!<SUBROUTINE NAME ="atmos_co2_gather_data">
!<OVERVIEW>
!  A subroutine to gather fields needed for calculating the CO2 gas flux
!</OVERVIEW>
!

subroutine atmos_co2_gather_data (gas_fields, tr_bot)

use coupler_types_mod, only: coupler_2d_bc_type, ind_pcair

implicit none

!-----------------------------------------------------------------------

type(coupler_2d_bc_type), intent(inout) :: gas_fields
real, dimension(:,:,:), intent(in)      :: tr_bot

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_gather_data'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

if (ind_co2_flux .gt. 0) then
  gas_fields%bc(ind_co2_flux)%field(ind_pcair)%values(:,:) = tr_bot(:,:,ind_co2)
endif

end subroutine atmos_co2_gather_data
!</SUBROUTINE >

!#######################################################################


!#######################################################################

!<SUBROUTINE NAME ="atmos_co2_flux_init">

!<OVERVIEW>
! Subroutine to initialize the carbon dioxide flux
!</OVERVIEW>

 subroutine atmos_co2_flux_init

use atmos_ocean_fluxes_mod, only: aof_set_coupler_flux

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: n

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_flux_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'


if ( .not. module_is_initialized) then

!----- set initial value of carbon ------------

  call tracer_manager_init      ! need to call here since the ocean pes never call it
  n = get_tracer_index(MODEL_ATMOS,'co2')
  if (n > 0) then
    ind_co2 = n
    if (ind_co2 > 0) then
      write (stdout(),*) trim(note_header), ' CO2 was initialized as tracer number ', ind_co2
      write (stdlog(),*) trim(note_header), ' CO2 was initialized as tracer number ', ind_co2
    endif
  endif
  module_is_initialized = .TRUE.
endif

!
!       initialize coupler flux
!

if (ind_co2 > 0) then
  ind_co2_flux = aof_set_coupler_flux('co2_flux',                       &
       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',       &
       atm_tr_index = ind_co2, param = (/ 9.36e-07, 9.7561e-06 /),      &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
endif

!-----------------------------------------------------------------------

end subroutine atmos_co2_flux_init
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME ="atmos_co2_init">

!<OVERVIEW>
! Subroutine to initialize the carbon dioxide module.
!</OVERVIEW>

 subroutine atmos_co2_init

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: n

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'


if (module_is_initialized) return

!----- set initial value of carbon ------------

n = get_tracer_index(MODEL_ATMOS,'co2')
if (n > 0) then
  ind_co2 = n
  if (ind_co2 > 0) then
    write (stdout(),*) trim(note_header), ' CO2 was initialized as tracer number ', ind_co2
    write (stdlog(),*) trim(note_header), ' CO2 was initialized as tracer number ', ind_co2
  endif

endif

call write_version_number (version, tagname)
module_is_initialized = .TRUE.


!-----------------------------------------------------------------------

end subroutine atmos_co2_init
!</SUBROUTINE>


!<SUBROUTINE NAME ="atmos_co2_end">
subroutine atmos_co2_end

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_co2_end'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

   module_is_initialized = .FALSE.

end subroutine atmos_co2_end
!</SUBROUTINE>

end module atmos_co2_mod
