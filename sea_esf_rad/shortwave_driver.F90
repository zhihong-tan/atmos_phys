                     module shortwave_driver_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
!
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!  smf
! </REVIEWER>
! 
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
!
! <OVERVIEW>
!  Code to carry out shortwave calculation.
! </OVERVIEW>
! <DESCRIPTION>
!  This code initializes, prepares, and ends shortwave radiation calculation.
!  This code is called by sea_esf_rad.f90 and calls shortwave subroutines
!  to do shortwave flux calculation.
! </DESCRIPTION>
!

!   shared modules:

use fms_mod,              only:  open_namelist_file, fms_init, &
                                 mpp_pe, mpp_root_pe, stdlog, &
                                 file_exist, write_version_number, &
                                 check_nml_error, error_mesg, &
                                 FATAL, NOTE, WARNING, close_file

!   shared radiation package modules:
 
use rad_utilities_mod,    only:  rad_utilities_init, Rad_control,  &
                                 Environment, cldrad_properties_type, &
                                 cld_specification_type, Sw_control, &
                                 radiative_gases_type,   &
                                 aerosol_type, aerosol_properties_type,&
                                 atmos_input_type, surface_type, &
                                 astronomy_type, sw_output_type, &
                                 cld_space_properties_type
use esfsw_parameters_mod, only:  esfsw_parameters_init

!  radiation package modules:

use lhsw_driver_mod,      only:  lhsw_driver_init, swrad
use esfsw_driver_mod,     only:  esfsw_driver_init, swresf,   &
                                 esfsw_driver_end

!-------------------------------------------------------------------

implicit none
private

!------------------------------------------------------------------
!    shortwave_driver_mod is the driver for shortwave radiation 
!    component of the sea_esf_rad radiation package.
!-----------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module  -------------------------

character(len=128)  :: version =  '$Id: shortwave_driver.F90,v 10.0 2003/10/24 22:00:47 fms Exp $'
character(len=128)  :: tagname =  '$Name: jakarta $'


!---------------------------------------------------------------------
!-------  interfaces --------

public        &
          shortwave_driver_init , shortwave_driver,    & 
          shortwave_driver_end


private       &

!   called from shortwave_driver:
          shortwave_driver_alloc


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)   :: swform = '    '
  
 
namelist / shortwave_driver_nml /             &
                                     swform

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

logical :: module_is_initialized = .false.  ! module initialized ?


!-------------------------------------------------------------------
!-------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! <SUBROUTINE NAME="shortwave_driver_init">
!  <OVERVIEW>
!   Code that initializes shortwave radiation calculation.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Initialize utilities and radiation utilities if necessary. They
!   should have been initialized in the radiation initialiation subroutine
!   in the sea_esf_rad.f90. The code then reads in input.nml namelist
!   and logs input parameters to logfile. It uses lhsw or esfsw package
!   depends on namelist parameter. Initializes apropriate shortwave
!   package subroutines and set up the initialize parameter.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call shortwave_driver_init (latb, pref)
!  </TEMPLATE>
!  <IN NAME="latb" TYPE="real">
!   An array of model latitudes at cell boundaries [radians]
!  </IN>
!  <IN NAME="pref" TYPE="real">
!   An array containing two reference pressure profiles [pascals]
!  </IN>
! </SUBROUTINE>
subroutine shortwave_driver_init (latb, pref)

!---------------------------------------------------------------------
!    shortwave_driver_init is the constructor for shortwave_driver_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension(:),   intent(in) :: latb
real, dimension(:,:), intent(in) :: pref
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       latb      array of model latitudes at cell boundaries 
!                 [ radians ]
!                                
!       pref      array containing two reference pressure profiles 
!                 [ Pa ]
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer   :: unit, io, ierr

!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!                                
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!-------------------------------------------------------------------
      if (module_is_initialized) return

!-------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!-------------------------------------------------------------------
      call fms_init
      call rad_utilities_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=shortwave_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'shortwave_driver_nml')
        end do
10      call close_file (unit)
      endif
 
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
                       write (stdlog(), nml=shortwave_driver_nml)

!--------------------------------------------------------------------
!    define logicals specifying the sw package in use as an element
!    of a shortwave_control_type variable usable by other radiation-
!    related modules. initialize the modules associated with the chosen
!    sw package.
!---------------------------------------------------------------------
      if (trim(swform) == 'lhsw') then
        Sw_control%do_lhsw  = .true.
        if (Environment%running_sa_model) then
          call error_mesg ( 'shortwave_driver_mod', &
            'must use esfsw parameterization with sa_gcm',  FATAL)
        endif
        call lhsw_driver_init (pref)
      else if (trim(swform) == 'esfsw99') then
        Sw_control%do_esfsw = .true.
        call esfsw_parameters_init
        call esfsw_driver_init
      else
        call error_mesg ( 'shortwave_driver_mod',   &
        'improper specification of desired shortwave parameterization',&
                                                               FATAL)
      endif

!---------------------------------------------------------------------
!    mark the just-defined logicals as defined.
!---------------------------------------------------------------------
      Sw_control%do_lhsw_iz  = .true.
      Sw_control%do_esfsw_iz = .true.

!-------------------------------------------------------------------
!    set flag indicating successful initialization of module.
!-------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------


end subroutine shortwave_driver_init



!###########################################################
! <SUBROUTINE NAME="shortwave_driver">
!  <OVERVIEW>
!   Code that deploys shortwave radiation calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!    shortwave_driver initializes shortwave radiation output variables, 
!    determines if shortwave radiation is present in the current physics
!    window, selects one of the available shortwave parameterizations,
!    executes it, and returns the output fields to sea_esf_rad_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call shortwave_driver (is, ie, js, je, Atmos_input, Surface, Astro, &
!                           Rad_gases, Cldrad_props, Cld_spec, Sw_output, &
!                           Cldspace_rad) 
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!    starting subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!    ending subdomain i indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="js" TYPE="integer">
!    starting subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="je" TYPE="integer">
!    ending subdomain j indices of data in 
!                   the physics_window being integrated
!  </IN>
!   <IN NAME="Atmos_input" TYPE="atmos_input_type">
!     Atmos_input_type variable containing the atmospheric
!     input fields on the radiation grid 
!   </IN>
!   <IN NAME="Astro" TYPE="astronomy_type">
!     Astronomy_type variable containing the astronomical
!     input fields on the radiation grid  
!   </IN>
!   <IN NAME="Aerosol" TYPE="aerosol_type">
!    Aerosol input data of shortwave radiation calculation
!   </IN>
!   <IN NAME="Aerosol_props" TYPE="aerosol_properties_type">
!    Aerosol radiative property input data 
!   </IN>
!   <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!     Radiative_gases_type variable containing the radiative 
!     gas input fields on the radiation grid 
!   </IN>
!   <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!     The cloud radiative property input fields on the
!     radiation grid
!   </IN>
!   <INOUT NAME="Sw_output" TYPE="sw_output_type">
!     The shortwave radiation calculation result
!   </INOUT>
!   <INOUT NAME="Cldspace_rad" TYPE="cld_space_properties_type">
!     Optional cloud radiative forcing output used in lacis-hansen
!     formulation.
!   </INOUT>
!   <IN NAME="Surface" TYPE="surface_type">
!    Surface data as boundary condition to radiation
!   </IN>
!   <IN NAME="Cld_spec" TYPE="cld_specification_type">
!    Cloud specification data as initial condition to radiation
!   </IN>
! </SUBROUTINE>
subroutine shortwave_driver (is, ie, js, je, Atmos_input, Surface,  &
                             Astro, Aerosol, Aerosol_props, Rad_gases, &
                             Cldrad_props,  Cld_spec, Sw_output,   &
                             Cldspace_rad) 

!---------------------------------------------------------------------
!    shortwave_driver initializes shortwave radiation output variables, 
!    determines if shortwave radiation is present in the current physics
!    window, selects one of the available shortwave parameterizations,
!    executes it, and returns the output fields to sea_esf_rad_mod.
!---------------------------------------------------------------------

integer,                         intent(in)    :: is, ie, js, je
type(atmos_input_type),          intent(in)    :: Atmos_input     
type(surface_type),              intent(in)    :: Surface     
type(astronomy_type),            intent(in)    :: Astro           
type(radiative_gases_type),      intent(in)    :: Rad_gases   
type(aerosol_type),              intent(in)    :: Aerosol     
type(aerosol_properties_type),   intent(in)    :: Aerosol_props
type(cldrad_properties_type),    intent(in)    :: Cldrad_props
type(cld_specification_type),    intent(in)    :: Cld_spec
type(sw_output_type),            intent(inout) :: Sw_output
type(cld_space_properties_type), intent(inout) :: Cldspace_rad

!--------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Atmos_input    atmos_input_type variable containing the atmos-
!                     pheric input fields needed by the radiation 
!                     package
!      Surface        surface_type variable containing the surface input
!                     fields needed by the radiation package
!      Astro          astronomy_type variable containing the astronom-
!                     ical input fields needed by the radiation package
!      Rad_gases      radiative_gases_type variable containing the radi-
!                     ative gas input fields needed by the radiation 
!                     package
!      Aerosol        aerosol_type variable containing the aerosol input
!                     data needed by the radiation package
!      Aerosol_props  aerosol_properties_type variable containing the
!                     aerosol radiative properties input data needed by
!                     the radiation package
!      Cldrad_props   cldrad_properties_type variable containing the 
!                     cloud radiative property input fields needed by 
!                     the radiation package
!      Cld_spec       cld_specification_type variable containing the 
!                     cloud specification input fields needed by the 
!                     radiation package
!
!   intent(out) variables:
!
!      Sw_output      sw_output_type variable containing shortwave 
!                     radiation output data 
!      Cldspace_rad   cld_space_properties_type variable containing the
!                     sw output fields obtained from the lacis-hansen
!                     shortwave parameterization 
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:

      logical  :: skipswrad
      logical  :: with_clouds
      integer  :: i, j       
      integer  :: ix, jx, kx

!---------------------------------------------------------------------
!   local variables:
!
!      skipswrad    bypass calling sw package because sun is not 
!                   shining any where in current physics window ?
!      with_clouds  are clouds to be considered in determining
!                   the sw fluxes and heating rates ?
!      ix,jx,kx     dimensions of current physics window
!      i,j          do-loop indices
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('shortwave_driver_mod',   &
             'module has not been initialized', FATAL )
      endif

!----------------------------------------------------------------------
!    call shortwave_driver_alloc to initialize shortwave fluxes and 
!    heating rates.
!--------------------------------------------------------------------
      ix = ie - is + 1
      jx = je - js + 1
      kx = size (Atmos_input%press,3) - 1
      call shortwave_driver_alloc (ix, jx, kx, Sw_output) 

!--------------------------------------------------------------------
!    determine when the no-sun case exists at all points within the 
!    physics window and bypass the sw radiation calculations for that 
!    window. for do_annual_mean or do_daily_mean, only one cosz in a
!    model row need be tested, since all points in i have the same 
!    zenith angle.
!--------------------------------------------------------------------
      skipswrad = .true.
      do j=1,jx        
        if ( Astro%cosz(1,j) > 0.0 ) skipswrad = .false.
        if (Sw_control%do_diurnal) then
          do i = 2,ix         
            if (Astro%cosz(i,j) > 0.0 )  then
              skipswrad = .false.
              exit
            endif
          end do
        endif
      end do

!--------------------------------------------------------------------
!    if the sun is shining nowhere in the physics window allocate
!    output fields which will be needed later, set them to a flag
!    value and return.
!--------------------------------------------------------------------
      if (skipswrad)  then
        allocate ( Cldspace_rad%camtswkc(ie-is+1, je-js+1, 1 ))
        allocate ( Cldspace_rad%cirabswkc(ie-is+1, je-js+1, 1 ))
        allocate ( Cldspace_rad%cirrfswkc(ie-is+1, je-js+1, 1 ))
        allocate ( Cldspace_rad%cvisrfswkc(ie-is+1, je-js+1, 1 ))
        allocate ( Cldspace_rad%ktopswkc(ie-is+1, je-js+1,  1 ))
        allocate ( Cldspace_rad%kbtmswkc(ie-is+1, je-js+1,  1 ))
        Cldspace_rad%camtswkc = -99.0
        Cldspace_rad%cirabswkc = -99.0
        Cldspace_rad%cirrfswkc = -99.0
        Cldspace_rad%cvisrfswkc = -99.0
        Cldspace_rad%ktopswkc = -99.0
        Cldspace_rad%kbtmswkc = -99.0

!---------------------------------------------------------------------
!    calculate shortwave radiative forcing and fluxes using the 
!    exponential-sum-fit parameterization.
!---------------------------------------------------------------------
      else if (Sw_control%do_esfsw) then
        call swresf (is, ie, js, je, Atmos_input, Surface, Rad_gases, &
                     Aerosol, Aerosol_props, Astro, Cldrad_props,  &
                     Cld_spec, Sw_output)

!---------------------------------------------------------------------
!    calculate shortwave radiative forcing and fluxes using the 
!    lacis-hansen parameterization.
!---------------------------------------------------------------------
      else if (Sw_control%do_lhsw) then
        with_clouds = .true.
        call swrad (is, ie, js, je, Astro, with_clouds, Atmos_input, &
                    Surface, Rad_gases, Cldrad_props, Cld_spec,  &
                    Sw_output, Cldspace_rad)

!---------------------------------------------------------------------
!    lacis-hansen requires a second call to produce the cloud-free
!    fluxes.
!---------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
          with_clouds = .false.
          call swrad (is, ie, js, je, Astro, with_clouds, Atmos_input, &
                      Surface, Rad_gases, Cldrad_props, Cld_spec,  &
                      Sw_output, Cldspace_rad)  
        endif
      endif  

!--------------------------------------------------------------------


end subroutine shortwave_driver



!###################################################################
! <SUBROUTINE NAME="shortwave_driver_end">
!  <OVERVIEW>
!   Code that ends shortwave radiation calculation
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine that simply reset shortwave_driver_initialized to false
!  </DESCRIPTION>
!  <TEMPLATE>
!   call shortwave_driver_end
!  </TEMPLATE>
! </SUBROUTINE>
subroutine shortwave_driver_end

!---------------------------------------------------------------------
!    shortwave_driver_end is the destructor for shortwave_driver_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('shortwave_driver_mod',   &
             'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    close out the modules initialized by this module.
!--------------------------------------------------------------------
      if (Sw_control%do_esfsw) then
        call esfsw_driver_end
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!-------------------------------------------------------------------

end subroutine shortwave_driver_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!###################################################################
! <SUBROUTINE NAME="shortwave_driver_alloc">
!  <OVERVIEW>
!   Code that allocates and initializes shortwave output variables
!  </OVERVIEW>
!  <DESCRIPTION>
!   Shortwave_driver_alloc allocates and initializes the components
!   of the sw_output_type variable Sw_output, which is used to hold
!   output data from shortwave_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call shortwave_driver_alloc (ix, jx, kx, Sw_output)
!  </TEMPLATE>
!  <IN NAME="ix" TYPE="integer">
!   x dimention of the radiation grid where shortwave output is desired
!  </IN>
!  <IN NAME="jx" TYPE="integer">
!   y dimention of the radiation grid where shortwave output is desired
!  </IN>
!  <IN NAME="kx" TYPE="integer">
!   z dimention of the radiation grid where shortwave output is desired
!  </IN>
!  <INOUT NAME="Sw_output" TYPE="sw_output_type">
!   shortwave radiation output variable
!  </INOUT>
! </SUBROUTINE>
!
subroutine shortwave_driver_alloc (ix, jx, kx, Sw_output) 

!--------------------------------------------------------------------
!    shortwave_driver_alloc allocates and initializes the components
!    of the sw_output_type variable Sw_output, which is used to hold
!    output data from shortwave_driver_mod.
!--------------------------------------------------------------------

integer,              intent(in)     ::  ix, jx, kx 
type(sw_output_type), intent(inout)  ::  Sw_output 

!-------------------------------------------------------------------
!  intent(in) variables:
!
!    ix, jx, kx   dimensions of the radiation grid on which output 
!                 will be produced
!
!  intent(inout) variables:
!
!      Sw_output  sw_output_type variable containing shortwave 
!                 radiation output data 
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    allocate and initialize fields to contain net(up-down) sw flux 
!    (fsw), upward sw flux (ufsw), downward sw flux(dfsw) at flux 
!    levels and sw heating in model layers (hsw).
!--------------------------------------------------------------------
      allocate (Sw_output%fsw  (ix, jx, kx+1) )
      allocate (Sw_output%ufsw (ix, jx, kx+1) )
      allocate (Sw_output%dfsw (ix, jx, kx+1) )
      allocate (Sw_output%hsw  (ix, jx, kx  ) )

      Sw_output%fsw   (:,:,:) = 0.0
      Sw_output%dfsw  (:,:,:) = 0.0
      Sw_output%ufsw  (:,:,:) = 0.0
      Sw_output%hsw   (:,:,:) = 0.0

!---------------------------------------------------------------------
!    if the cloud-free values are desired, allocate and initialize 
!    arrays for the fluxes and heating rate in the absence of clouds.
!----------------------------------------------------------------------
      if (Rad_control%do_totcld_forcing) then
        allocate (Sw_output%fswcf  (ix, jx, kx+1) )
        allocate (Sw_output%dfswcf (ix, jx, kx+1) )
        allocate (Sw_output%ufswcf (ix, jx, kx+1) )
        allocate (Sw_output%hswcf  (ix, jx, kx  ) )

        Sw_output%fswcf (:,:,:) = 0.0
        Sw_output%dfswcf(:,:,:) = 0.0
        Sw_output%ufswcf(:,:,:) = 0.0
        Sw_output%hswcf (:,:,:) = 0.0
      endif

!--------------------------------------------------------------------

end  subroutine shortwave_driver_alloc



!####################################################################


                end module shortwave_driver_mod

