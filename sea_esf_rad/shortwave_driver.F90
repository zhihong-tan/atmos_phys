                     module shortwave_driver_mod

use  utilities_mod,       only:  open_file, file_exist,      &
                                 utilities_init, &
                                 check_nml_error, error_mesg, & 
                                 print_version_number, FATAL, NOTE, &
                                 WARNING, get_my_pe, close_file
use rad_utilities_mod,    only:  Rad_control, radiation_control_type, &
                                 cldrad_properties_type, &
                                 rad_utilities_init, &
                                 radiative_gases_type,   &
                                 atmos_input_type, &
                                 astronomy_type, &
                                 cld_space_properties_type, &
                                 sw_output_type, &
                                 Sw_control, shortwave_control_type, &
                                 Environment, environment_type

!  lacis-hansen shortwave package:

use lhsw_driver_mod,      only:  lhsw_driver_init, swrad

!  exponential-sum-fit shortwave package:

use esfsw_driver_mod,     only:  esfsw_driver_init, swresf
use esfsw_scattering_mod, only:  esfsw_scattering_init
use esfsw_parameters_mod, only:  esfsw_parameters_init


implicit none
private

!------------------------------------------------------------------
!    shortwave_driver_mod is the driver for shortwave radiation 
!    component of the sea_esf_rad radiation package.
!-----------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module  -------------------------

character(len=128)  :: version =  '$Id: shortwave_driver.F90,v 1.3 2002/07/16 22:37:02 fms Exp $'
character(len=128)  :: tag     =  '$Name: havana $'


!---------------------------------------------------------------------
!-------  interfaces --------

public   shortwave_driver_init , shortwave_driver,              &
         shortwave_driver_end


private  shortwave_driver_alloc


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)   :: swform = '    '
  
 
namelist / shortwave_driver_nml /             &
                                     swform


!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

logical    :: shortwave_driver_initialized =    &
                                    .false.  ! module initialized ?



!-------------------------------------------------------------------
!-------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
!       latb      array of model latitudes at cell boundaries [radians]
!                                
!       pref      array containing two reference pressure profiles 
!                 [pascals]
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables

      integer   :: unit, io, ierr


!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!-------------------------------------------------------------------
      if (shortwave_driver_initialized) return

!-------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!-------------------------------------------------------------------
      call utilities_init
      call rad_utilities_init

!----------------------------------------------------------------
!    read namelist.
!----------------------------------------------------------------
 
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=shortwave_driver_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'shortwave_driver_nml')
        enddo
10      call close_file (unit)
      endif
  
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      unit = open_file ('logfile.out', action='append')
      if (get_my_pe() == 0) then
!     if (get_my_pe() == get_root_pe() ) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit,nml=shortwave_driver_nml)
      endif
      call close_file (unit)

!--------------------------------------------------------------------
!    define logicals specifying the sw package in use as an element
!    of a shortwave_control_type variable usable by other radiation-
!    related modules. initialize the modules associated with the chosen
!    sw package.
!---------------------------------------------------------------------
      if (trim(swform) == 'lhsw') then
        Sw_control%do_lhsw  = .true.
        call lhsw_driver_init (pref)
      else if (trim(swform) == 'esfsw99') then
        Sw_control%do_esfsw = .true.
        call esfsw_parameters_init
        call esfsw_driver_init
        call esfsw_scattering_init (latb)
      else
        call error_mesg ( 'shortwave_driver_init',   &
        'improper specification of desired shortwave parameterization',&
                                                               FATAL)
      endif

!-------------------------------------------------------------------
!    set flag indicating successful initialization of module.
!-------------------------------------------------------------------
      shortwave_driver_initialized = .true.



end subroutine shortwave_driver_init



!###########################################################

subroutine shortwave_driver (is, ie, js, je, Atmos_input, Astro, &
                             Rad_gases, Cldrad_props, Sw_output, &
                             Cldspace_rad) 

!---------------------------------------------------------------------
!    shortwave_driver initializes shortwave radiation output variables, 
!    determines if shortwave radiation is present in the current physics
!    window, selects one of the available shortwave parameterizations,
!    executes it, and returns the output fields to sea_esf_rad_mod.
!---------------------------------------------------------------------

integer,                      intent(in)    :: is, ie, js, je
type(atmos_input_type),       intent(in)    :: Atmos_input     
type(astronomy_type),         intent(in)    :: Astro           
type(radiative_gases_type),   intent(in)    :: Rad_gases   
type(cldrad_properties_type), intent(in)    :: Cldrad_props
type(sw_output_type),         intent(out)   :: Sw_output
type(cld_space_properties_type), intent(out) :: Cldspace_rad

!--------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!      Atmos_input  atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!      Astro        astronomy_type variable containing the astronomical
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
!      Sw_output    sw_output_type variable containing shortwave 
!                   radiation output data 
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables
!--------------------------------------------------------------------
      logical  :: skipswrad
      logical  :: with_clouds
      integer  ::   j, i   
      integer  :: ix, jx, kx

!---------------------------------------------------------------------
!
!   local variables:
!
!      skipswrad    bypass calling sw package because sun is not 
!                   shining any where in current physics window ?
!      with_clouds  are clouds to be considered in determining
!                   the sw fluxes and heating rates ?
!
!---------------------------------------------------------------------


!----------------------------------------------------------------------
!    initialize fluxes and heating rates which are kept in a 
!    shortwave_output_type variable.
!--------------------------------------------------------------------
      ix = ie - is + 1
      jx = je - js + 1
      kx = size (Atmos_input%press,3) - 1
      call shortwave_driver_alloc (ix, jx, kx, Sw_output) 

!--------------------------------------------------------------------
!    determine when the no-sun case exists at all points within the 
!    physics window and bypass the radiation calculations for that 
!    window. for do_annual_mean or do_daily_mean, only one cosz in a
!    model row need be tested, since all points in i have the same 
!    zenith angle.
!--------------------------------------------------------------------
      skipswrad = .true.
      do j=1,jx        
        if ( Astro%cosz(1,j) > 0.0 ) skipswrad = .false.
        if (Astro%do_diurnal) then
          do i = 2,ix         
            if (Astro%cosz(i,j) > 0.0 )  then
              skipswrad = .false.
              exit
            endif
          end do
        endif
      end do

!--------------------------------------------------------------------
!    if the sun is shining nowhere in the physics window, return.
!--------------------------------------------------------------------
      if (skipswrad)  then

!---------------------------------------------------------------------
!    calculate shortwave radiative forcing and fluxes using the 
!    exponential-sum-fit parameterization.
!---------------------------------------------------------------------
      else if (Sw_control%do_esfsw) then
        call swresf (is, ie, js, je, Atmos_input, Rad_gases,    &
                     Astro, Cldrad_props, Sw_output)

!---------------------------------------------------------------------
!    calculate shortwave radiative forcing and fluxes using the 
!    lacis-hansen parameterization.
!---------------------------------------------------------------------
      else if (Sw_control%do_lhsw) then
        with_clouds = .true.
        call swrad (is, ie, js, je, Astro, with_clouds,  &
                    Atmos_input, Rad_gases, Cldrad_props, Sw_output, &
                    Cldspace_rad)

!---------------------------------------------------------------------
!    lacis-hansen requires a second call to produce the cloud-free
!    fluxes.
!---------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
          with_clouds = .false.
          call swrad (is, ie, js, je, Astro, with_clouds,  &
                      Atmos_input, Rad_gases, Cldrad_props, Sw_output, &
                      Cldspace_rad)  
        endif
      endif  

!--------------------------------------------------------------------

end subroutine shortwave_driver



!###################################################################

subroutine shortwave_driver_end

!---------------------------------------------------------------------
!    shortwave_driver_end is the destructor for shortwave_driver_mod.
!---------------------------------------------------------------------

      shortwave_driver_initialized = .false.


end subroutine shortwave_driver_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!###################################################################

subroutine shortwave_driver_alloc (ix, jx, kx, Sw_output) 

!--------------------------------------------------------------------
!    shortwave_driver_alloc allocates and initializes the components
!    of the sw_output_type variable Sw_output, which is used to hold
!    output data from shortwave_driver_mod.
!--------------------------------------------------------------------
integer,              intent(in)   ::  ix, jx, kx 
type(sw_output_type), intent(out)  ::  Sw_output 

!-------------------------------------------------------------------
!  intent(in) variables:
!
!    ix, jx, kx   dimensions of the radiation grid on which output 
!                 will be produced
!
!  intent(out) variables:
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

