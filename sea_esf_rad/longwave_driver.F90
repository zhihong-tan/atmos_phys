                     module longwave_driver_mod

use  utilities_mod,     only: open_file, file_exist,    &
                              check_nml_error, error_mesg, &
                              FATAL, NOTE, WARNING, get_my_pe, &
			      close_file,  utilities_init
use rad_utilities_mod,  only: rad_utilities_init, &
                              cldrad_properties_type, &
                              lw_output_type, &
                              atmos_input_type, &
                              radiative_gases_type, &
                              lw_table_type, &
                              lw_diagnostics_type, &
                              radiation_control_type, Rad_control

!   simplified-exchange-approximation longwave package:

use sealw99_mod,        only: sealw99_init, sealw99, sealw99_end

!------------------------------------------------------------------


implicit none
private

!-------------------------------------------------------------------
!    longwave_driver_mod is the driver for the longwave radiation
!    component of the sea_esf_rad radiation package.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: longwave_driver.F90,v 1.3 2002/07/16 22:35:35 fms Exp $'
character(len=128)  :: tag     =  '$Name: havana $'

!---------------------------------------------------------------------
!-------  interfaces --------

public   longwave_driver_init, longwave_driver,   &
         longwave_driver_end

private  longwave_driver_alloc


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16) :: lwform= 'sealw99'
 

namelist / longwave_driver_nml /    &
                                 lwform


!---------------------------------------------------------------------
!------- public data ------




!---------------------------------------------------------------------
!------- private data ------

logical :: longwave_driver_initialized =   &
                                      .false.   ! module initialized ?
logical :: do_sealw99 = .false.                 ! sealw99 parameter-
                                                ! ization active ?

!---------------------------------------------------------------------
!---------------------------------------------------------------------




                          contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!                    PUBLIC SUBROUTINES
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



subroutine longwave_driver_init (latb, lonb, pref, Lw_tables)
 
!---------------------------------------------------------------------
!    longwave_driver_init is the constructor for longwave_driver_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension(:),     intent(in)  :: latb, lonb
real, dimension(:,:),   intent(in)  :: pref
type(lw_table_type),    intent(out) :: Lw_tables

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       latb      array of model latitudes at cell boundaries [radians]
!       lonb      array of model longitudes at cell boundaries [radians]
!       pref      array containing two reference pressure profiles 
!                 [pascals]
!
!  intent(out) variables:
!
!       Lw_tables lw_tables_type variable containing various longwave
!                 table specifiers needed by radiation_diag_mod.
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!  local variables

      integer     :: unit, ierr, io


!--------------------------------------------------------------------
!    if routine has already been executed, return.
!--------------------------------------------------------------------
      if (longwave_driver_initialized) return

!-------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!-------------------------------------------------------------------
      call utilities_init
      call rad_utilities_init

!----------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if (file_exist('input.nml')) then
        unit = open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=longwave_driver_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'longwave_driver_nml')
        enddo
10      call close_file (unit)
      endif

!----------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      unit = open_file ('logfile.out', action='append')
      if (get_my_pe() == 0) then
!     if (get_my_pe() == get_root_pe() ) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit,nml=longwave_driver_nml)
      endif
      call close_file (unit)


!--------------------------------------------------------------------
!    determine if valid specification of lw radiation has been made.
!    if optional packages are provided at some later time, this is where
!    the choice of package will be made.
!---------------------------------------------------------------------
      if (trim(lwform) == 'sealw99') then
        do_sealw99 = .true.
        call sealw99_init ( latb, lonb, pref, Lw_tables)
      else
        call error_mesg ( 'longwave_driver_mod', &
	  'invalid longwave radiation form specified', FATAL)
      endif

!---------------------------------------------------------------------
!    set flag indicating successful initialization of module.
!---------------------------------------------------------------------
      longwave_driver_initialized = .true.



end  subroutine longwave_driver_init




!#####################################################################

subroutine longwave_driver (is, ie, js, je, Atmos_input, Rad_gases, &
                            Cldrad_props, Lw_output, Lw_diagnostics)

!--------------------------------------------------------------------
!    longwave_driver allocates and initializes longwave radiation out-
!    put variables and selects an available longwave radiation param-
!    eterization, executes it, and then returns the output fields to 
!    sea_esf_rad_mod.
!--------------------------------------------------------------------

integer,                      intent(in)  :: is, ie, js, je
type(atmos_input_type),       intent(in)  :: Atmos_input  
type(radiative_gases_type),   intent(in)  :: Rad_gases   
type(cldrad_properties_type), intent(in)  :: Cldrad_props
type(lw_output_type),         intent(out) :: Lw_output   
type(lw_diagnostics_type),    intent(out) :: Lw_diagnostics

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!      Atmos_input  atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!      Rad_gases    radiative_gases_type variable containing the radi-
!                   ative gas input fields needed by the radiation 
!                   package
!      Cldrad_props cldrad_properties_type variable containing the 
!                   cloud radiative property input fields needed by the 
!                   radiation package
!
!   intent(out) variables:
!
!      Lw_output    lw_output_type variable containing longwave 
!                   radiation output data 
!      Lw_diagnostics
!                   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!  
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables

      integer  :: ix, jx, kx

!--------------------------------------------------------------------
!    call longwave_driver_alloc to allocate component arrays of a
!    lw_output_type variable.
!----------------------------------------------------------------------
      ix = ie - is + 1
      jx = je - js + 1
      kx = size (Atmos_input%press,3) - 1
      call longwave_driver_alloc (ix, jx, kx, Lw_output)

!--------------------------------------------------------------------
!    calculate the longwave radiative heating rates and fluxes.
!--------------------------------------------------------------------
      if (do_sealw99) then

!--------------------------------------------------------------------
!    call sealw99 to use the simplified-exchange-approximation (sea)
!    parameterization.
!----------------------------------------------------------------------
        call sealw99 (is, ie, js, je, Atmos_input, Rad_gases, &
                      Cldrad_props, Lw_output, Lw_diagnostics)
      else

!--------------------------------------------------------------------
!    at the current time sealw99 is the only longwave parameterization 
!    available.
!----------------------------------------------------------------------
        call error_mesg ('longwave_driver_mod', &
         'invalid longwave radiation parameterization selected', FATAL)
      endif

!---------------------------------------------------------------------


end subroutine longwave_driver


!#####################################################################

subroutine longwave_driver_end                  

!--------------------------------------------------------------------
!    longwave_driver_end is the destructor for longwave_driver_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    call sealw99_end to close sealw99_mod.
!-------------------------------------------------------------------
      if (do_sealw99) then
        call sealw99_end
      endif

!--------------------------------------------------------------------
      longwave_driver_initialized = .false.



end subroutine longwave_driver_end                  


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################

subroutine longwave_driver_alloc (ix, jx, kx, Lw_output)

!--------------------------------------------------------------------
!    longwave_driver_alloc allocates and initializes the components
!    of the lw_output_type variable Lw_output which holds the longwave
!    output needed by radiation_driver_mod.
!--------------------------------------------------------------------

integer,                   intent(in)  :: ix, jx, kx
type(lw_output_type),      intent(out) :: Lw_output

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      ix,jx,kx     (i,j,k) lengths of radiation arrays to be allocated
!
!
!   intent(out) variables:
!
!      Lw_output    lw_output_type variable containing longwave 
!                   radiation output data 
!  
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    allocate and initialize arrays to hold net longwave fluxes and 
!    the longwave heating rate at each gridpoint. if the
!    cloud-forcing calculation is to be done, also allocate and init-
!    ialize arrays for fluxes and heating rates without clouds.
!-------------------------------------------------------------------
      allocate (Lw_output%flxnet( ix, jx, kx+1) )
      allocate (Lw_output%heatra( ix, jx, kx  ) )
      Lw_output%flxnet(:,:,:) = 0.0
      Lw_output%heatra(:,:,:) = 0.0
      if (Rad_control%do_totcld_forcing)  then
        allocate (Lw_output%flxnetcf( ix, jx, kx+1) )
        allocate (Lw_output%heatracf( ix, jx, kx  ) )
        Lw_output%flxnetcf(:,:,:) = 0.0
        Lw_output%heatracf(:,:,:) = 0.0
      endif
    
!--------------------------------------------------------------------

end subroutine longwave_driver_alloc



!###################################################################



		  end module longwave_driver_mod
