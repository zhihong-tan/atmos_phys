                    module sea_esf_rad_mod


use time_manager_mod,     only:  &
!                                time_manager_init, &
                                 time_type
use utilities_mod,        only:  get_my_pe, FATAL,  WARNING, NOTE,   &
                                 close_file, open_file, &
                                 print_version_number,  &
                                 get_num_pes, &
!                                get_root_pe, &
                                 check_nml_error, file_exist,  &
                                 utilities_init, error_mesg
use rad_utilities_mod,    only:  rad_utilities_init, &
                                 radiation_control_type, Rad_control, &
                                 define_environment, &
                                 radiative_gases_type, & 
                                 cldrad_properties_type, &
				 cld_space_properties_type, &
                                 astronomy_type,  atmos_input_type, &
                                 lw_diagnostics_type, &
				 cld_diagnostics_type, &
				 lw_table_type, &
                                 aerosol_type, &
                                 sw_output_type, lw_output_type, &
                                 rad_output_type, fsrad_output_type, &
                                 environment_type, Environment, &
                                 shortwave_control_type, Sw_control, &
				 longwave_control_type, Lw_control, &
                                 cloudrad_control_type, Cldrad_control
use fms_mod,              only:  mpp_clock_id, mpp_clock_begin, &
                                 mpp_clock_end, CLOCK_MODULE,  &
                                 MPP_CLOCK_SYNC

use radiation_diag_mod,   only:  radiation_diag_init,   &
                                 radiation_diag_driver, &
                                 radiation_diag_end
use longwave_driver_mod,  only:  longwave_driver_init,   &
                                 longwave_driver, &
                                 longwave_driver_end
use shortwave_driver_mod, only:  shortwave_driver_init,  &
                                 shortwave_driver, shortwave_driver_end


implicit none 
private 

!-----------------------------------------------------------------------
!     sea_esf_rad_mod is the driver for the sea_esf_rad 
!     radiation package.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!------------ version number for this module ---------------------------

character(len=128) :: version = '$Id: sea_esf_rad.F90,v 1.4 2003/04/09 21:01:38 fms Exp $'
character(len=128) :: tag = '$Name: inchon $'


!--------------------------------------------------------------------
!-- interfaces -----

public       &
            sea_esf_rad_init, sea_esf_rad, sea_esf_rad_end


private      deallocate_arrays


!-----------------------------------------------------------------------




!---------------------------------------------------------------------
!--- namelist ---

logical :: dummy


namelist /sea_esf_rad_nml/  dummy


!---------------------------------------------------------------------
!---- public data ----





!---------------------------------------------------------------------
!---- private data ----


logical :: sea_esf_rad_initialized = .false.    ! module initialized ?

integer :: longwave_clock, shortwave_clock      ! timing clocks

!---------------------------------------------------------------------
!---------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!######################################################################

subroutine sea_esf_rad_init (lonb, latb, pref_r)

!---------------------------------------------------------------------
!   sea_esf_rad_init is the constructor for sea_esf_rad_mod.
!---------------------------------------------------------------------

real, dimension(:),      intent(in)  :: lonb, latb
real, dimension(:,:),    intent(in)  :: pref_r

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       lonb      array of model longitudes on cell boundaries 
!                 [radians]
!       latb      array of model latitudes at cell boundaries 
!                 [radians]
!       pref_r    array containing two reference pressure profiles 
!                 on the radiation grid for use in defining 
!                 transmission functions 
!                 [pascals]
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables

      integer                           :: unit, io, ierr
      logical                           :: end
      type(lw_table_type)               :: Lw_tables

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (sea_esf_rad_initialized) return

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
      if ( file_exist('input.nml')) then
        unit = open_file (file='input.nml', action='read')
        ierr=1; do while (ierr /= 0)
          read  (unit, nml=sea_esf_rad_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'sea_esf_rad_nml')
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
        write (unit, nml=sea_esf_rad_nml)
      endif
      call close_file (unit)

!---------------------------------------------------------------------
!    initialize the modules called by this module.
!---------------------------------------------------------------------
      call longwave_driver_init  (latb, lonb, pref_r, Lw_tables)
      call shortwave_driver_init (latb, pref_r)
      call radiation_diag_init   (latb, lonb, Lw_tables)

!---------------------------------------------------------------------
!    initialize clocks to time various modules called by this module.
!---------------------------------------------------------------------
      longwave_clock =      &
                  mpp_clock_id ('   Physics_down: Radiation: lw', &
                   grain=CLOCK_MODULE, flags = MPP_CLOCK_SYNC)
      shortwave_clock =     &
                  mpp_clock_id ('   Physics_down: Radiation: sw', &
                  grain=CLOCK_MODULE, flags = MPP_CLOCK_SYNC)

!-------------------------------------------------------------------
      sea_esf_rad_initialized = .true.



!-------------------------------------------------------------------


end subroutine sea_esf_rad_init



!#####################################################################

subroutine sea_esf_rad (is, ie, js, je, Atmos_input, Astro, Rad_gases, &
                        Aerosol, Cldrad_props, Cld_diagnostics,   &
                        Lw_output, Sw_output)

!-----------------------------------------------------------------------
!     sea_esf_rad calls the modules which calculate the long- and short-
!     wave radiational heating terms and fluxes and the radiation diag-
!     nostics module which provides radiation package diagnostics.
!-----------------------------------------------------------------------

integer,                      intent(in)   :: is, ie, js, je
type(atmos_input_type),       intent(in)   :: Atmos_input
type(astronomy_type),         intent(in)   :: Astro
type(radiative_gases_type),   intent(in)   :: Rad_gases
type(aerosol_type),           intent(in)   :: Aerosol      
type(cldrad_properties_type), intent(in)   :: Cldrad_props
type(cld_diagnostics_type),   intent(in)   :: Cld_diagnostics
type(lw_output_type),         intent(inout)  :: Lw_output 
type(sw_output_type),         intent(inout)  :: Sw_output 

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!      Atmos_input  atmos_input_type variable containing the atmospheric
!                   input fields on the radiation grid          
!      Astro        astronomy_type variable containing the astronomical
!                   input fields on the radiation grid           
!      Rad_gases    radiative_gases_type variable containing the radi-
!                   ative gas input fields on the radiation grid   
!      Cldrad_props cldrad_properties_type variable containing the 
!                   cloud radiative property input fields on the
!                   radiation grid
!      Cld_diagnostics
!
!  intent(out) variables:
!
!      Lw_output    lw_output_type variable containing longwave 
!                   radiation output data from the sea_esf_rad
!                   radiation package on the radiation grid
!      Sw_output    sw_output_type variable containing shortwave 
!                   radiation output data from the sea_esf_rad
!                   radiation package on the radiation grid 
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables

      type(lw_diagnostics_type)         :: Lw_diagnostics
      type(cld_space_properties_type)   :: Cldspace_rad

!---------------------------------------------------------------------
!   local variables
!
!         Lw_diagnostics      lw_diagnostics_type variable to hold
!                             desired diagnostics from longwave_driver
!                             so they may be passed to 
!                             radiation_diag_mod
!
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!    compute longwave radiation.
!----------------------------------------------------------------------
      call mpp_clock_begin (longwave_clock)
      call longwave_driver (is, ie, js, je, Atmos_input, Rad_gases,  &
                            Aerosol, Cldrad_props, Lw_output,   &
                            Lw_diagnostics)
      call mpp_clock_end (longwave_clock)

!----------------------------------------------------------------------
!    compute shortwave radiation.
!----------------------------------------------------------------------
      call mpp_clock_begin (shortwave_clock)
      call shortwave_driver (is, ie, js, je, Atmos_input, Astro,   &
                             Aerosol, Rad_gases, Cldrad_props,    &
                             Sw_output, Cldspace_rad)
      call mpp_clock_end (shortwave_clock)

!--------------------------------------------------------------------
!    call radiation_diag_driver to compute radiation diagnostics at 
!    desired points.
!--------------------------------------------------------------------
      call radiation_diag_driver (is, ie, js, je, Atmos_input, Astro, &
                                  Rad_gases, Cldrad_props,    &
                                  Cld_diagnostics, Sw_output, &
                                  Lw_output, Lw_diagnostics, &
                                  Cldspace_rad)

!---------------------------------------------------------------------
!    call subroutine to deallocate the array components of the local 
!    derived -type variables.
!---------------------------------------------------------------------
      call deallocate_arrays (Lw_diagnostics, Cldspace_rad)


end subroutine sea_esf_rad




!###################################################################
      
subroutine sea_esf_rad_end
 
!-------------------------------------------------------------------
!    sea_esf_rad_end is the destructor for the sea_esf_rad module.
!-------------------------------------------------------------------

      call longwave_driver_end
      call shortwave_driver_end
      call radiation_diag_end
 
!--------------------------------------------------------------------
      sea_esf_rad_initialized = .false.



end subroutine sea_esf_rad_end


!####################################################################
      

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################

subroutine deallocate_arrays (Lw_diagnostics, Cldspace_rad)

!---------------------------------------------------------------------
!    deallocate_arrays deallocates the array cpomponents of local
!    derived-type variables.
!---------------------------------------------------------------------

type(lw_diagnostics_type),       intent(in)   :: Lw_diagnostics
type(cld_space_properties_type), intent(in)   :: Cldspace_rad

!---------------------------------------------------------------------
!  intent(in) variables:
!
!         Lw_diagnostics      lw_diagnostics_type variable to hold
!                             desired diagnostics from longwave_driver
!                             so they may be passed to 
!                             radiation_diag_mod
!         Cldspace_rad        cld_space_properties_type variable which
!                             holds lacis-hansen sw cloud-radiation
!                             variables in cloud-space, rather than 
!                             k-space, as the third dimension.
!
!----------------------------------------------------------------------


!--------------------------------------------------------------------
!    deallocate the components of Lw_diagnostics.
!--------------------------------------------------------------------
      deallocate (Lw_diagnostics%flx1e1)
      deallocate (Lw_diagnostics%fluxn )
      deallocate (Lw_diagnostics%cts_out)
      deallocate (Lw_diagnostics%cts_outcf)
      deallocate (Lw_diagnostics%gxcts )
      deallocate (Lw_diagnostics%excts )
      deallocate (Lw_diagnostics%exctsn)
      deallocate (Lw_diagnostics%fctsg )
      if (Lw_control%do_ch4_n2o) then
        deallocate (Lw_diagnostics%flx1e1f)
      endif
      if (Rad_control%do_totcld_forcing) then
        deallocate (Lw_diagnostics%fluxncf)
      endif

!--------------------------------------------------------------------
!    deallocate the components of Cldspace_rad. these arrays are only
!    allocated when the lh sw code is called with clouds present; 
!    therefore one must test for pointer association before deallo-
!    cating the memory.
!--------------------------------------------------------------------
      if (Sw_control%do_lhsw) then
        if (associated ( Cldspace_rad%camtswkc) ) then
          deallocate (Cldspace_rad%camtswkc )
          deallocate (Cldspace_rad%cirabswkc )
          deallocate (Cldspace_rad%cirrfswkc )
          deallocate (Cldspace_rad%cvisrfswkc )
          deallocate (Cldspace_rad%ktopswkc )
          deallocate (Cldspace_rad%kbtmswkc )
        endif
      endif

!--------------------------------------------------------------------


end subroutine deallocate_arrays 


!####################################################################



                 end module sea_esf_rad_mod

