module physics_radiation_exch_mod
#include <fms_platform.h>
!---------------------------------------------------------------------

!---- module data ----
use block_control_mod,  only: block_control_type

!---- public data ----
!---
!---Exch_ctrl
 public exchange_control_type
 type  exchange_control_type
     logical           :: doing_donner
     logical           :: doing_uw_conv
     logical           :: do_cosp
     logical           :: do_modis_yim
     logical           :: donner_meso_is_largescale
     integer           :: ncol
     character(len=16) :: cloud_type_form  ! indicator of radiatively active clouds
 end type  exchange_control_type

!---
!--- Radiation Flux block type
 public radiation_flux_block_type
 type radiation_flux_block_type
     real, dimension(:,:,:), _ALLOCATABLE :: tdt_rad           _NULL, &
                                             tdt_lw            _NULL, &
                                             extinction        _NULL
     real, dimension(:,:),   _ALLOCATABLE :: flux_sw           _NULL, &
                                             flux_sw_dir            _NULL, &
                                             flux_sw_dif            _NULL, &
                                             flux_sw_down_vis_dir   _NULL, &
                                             flux_sw_down_vis_dif   _NULL, &
                                             flux_sw_down_total_dir _NULL, &
                                             flux_sw_down_total_dif _NULL, &
                                             flux_sw_vis            _NULL, &
                                             flux_sw_vis_dir        _NULL, & 
                                             flux_sw_vis_dif        _NULL, &
                                             flux_lw                _NULL, &
                                             coszen                 _NULL
 end type radiation_flux_block_type

!--- Radiation Flux control type
 public radiation_flux_control_type
 type radiation_flux_control_type
     logical :: do_rad
 end type radiation_flux_control_type

!--- Rad_flux type definition
 public radiation_flux_type
 type radiation_flux_type
     type (radiation_flux_control_type) :: control
     type (radiation_flux_block_type), dimension(:), _ALLOCATABLE :: block _NULL
 end type radiation_flux_type


!---
!--- Moist Clouds block type
 public clouds_from_moist_block_type
 type clouds_from_moist_block_type
     real, dimension(:,:,:), _ALLOCATABLE :: cell_cld_frac           _NULL, &
                                             cell_liq_amt            _NULL, &
                                             cell_liq_size           _NULL, &
                                             cell_ice_amt            _NULL, &
                                             cell_ice_size           _NULL, &
                                             cell_droplet_number     _NULL, &
                                             meso_cld_frac           _NULL, &
                                             meso_liq_amt            _NULL, &
                                             meso_liq_size           _NULL, &
                                             meso_ice_amt            _NULL, &
                                             meso_ice_size           _NULL, &
                                             meso_droplet_number     _NULL, &
                                             lsc_cloud_area          _NULL, &
                                             lsc_liquid              _NULL, &
                                             lsc_ice                 _NULL, &
                                             lsc_droplet_number      _NULL, &
                                             lsc_ice_number          _NULL, &
                                             lsc_rain                _NULL, &
                                             lsc_snow                _NULL, &
                                             lsc_rain_size           _NULL, &
                                             lsc_snow_size           _NULL, &
                                             shallow_cloud_area      _NULL, &
                                             shallow_liquid          _NULL, &
                                             shallow_ice             _NULL, &
                                             shallow_droplet_number  _NULL, &
                                             shallow_ice_number      _NULL
    integer, dimension(:,:), _ALLOCATABLE :: nsum_out                _NULL
 end type clouds_from_moist_block_type

!--- Moist_clouds type definition
 public clouds_from_moist_type
 type clouds_from_moist_type
    type (clouds_from_moist_block_type), dimension(:), _ALLOCATABLE :: block _NULL
 end type clouds_from_moist_type


!---
!--- Cosp Rad block type
 public cosp_from_rad_block_type
 type cosp_from_rad_block_type
     real, dimension(:,:,:,:), _ALLOCATABLE :: tau_stoch        _NULL, &
                                               lwem_stoch       _NULL, &
                                               stoch_cloud_type _NULL, &
                                               stoch_conc_drop  _NULL, &
                                               stoch_conc_ice   _NULL, &
                                               stoch_size_drop  _NULL, &
                                               stoch_size_ice   _NULL
     real, dimension(:,:,:),   _ALLOCATABLE :: mr_ozone         _NULL
     real, dimension(:,:),     _ALLOCATABLE :: daytime          _NULL
 end type cosp_from_rad_block_type

!--- Cosp Rad control type
 public cosp_from_rad_control_type
 type cosp_from_rad_control_type
     logical :: step_to_call_cosp
 end type cosp_from_rad_control_type

!--- Cosp_rad type definition
 public cosp_from_rad_type
 type cosp_from_rad_type
     type (cosp_from_rad_control_type) :: control
     type (cosp_from_rad_block_type), dimension(:), _ALLOCATABLE :: block _NULL
 end type cosp_from_rad_type

public :: exch_rad_phys_state, &
          alloc_clouds_from_moist_type, dealloc_clouds_from_moist_type, &
          alloc_cosp_from_rad_type, dealloc_cosp_from_rad_type, &
          alloc_radiation_flux_type, dealloc_radiation_flux_type

contains


 subroutine alloc_cosp_from_rad_type (Cosp_rad, Exch_ctrl, Atm_block)
   type (cosp_from_rad_type),    intent(inout) :: Cosp_rad(:)
   type (exchange_control_type), intent(in)    :: Exch_ctrl
   type (block_control_type),    intent(in)    :: Atm_block
!--- local variables
   integer :: n, nb, ix, jx, npz

     npz = Atm_block%npz
     ncol = Exch_ctrl%ncol
     do n = 1, size(Cosp_rad,1)
       allocate (Cosp_rad(n)%block(Atm_block%nblks))
       if (Exch_ctrl%do_cosp .or. Exch_ctrl%do_modis_yim) then
         do nb = 1,Atm_block%nblks
           ix = Atm_block%ibe(nb)-Atm_block%ibs(nb)+1
           jx = Atm_block%jbe(nb)-Atm_block%jbs(nb)+1
           allocate ( Cosp_rad(n)%block(nb)%tau_stoch       (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%lwem_stoch      (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%stoch_cloud_type(ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%stoch_conc_drop (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%stoch_conc_ice  (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%stoch_size_drop (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%stoch_size_ice  (ix,jx,npz,ncol), &
                      Cosp_rad(n)%block(nb)%mr_ozone        (ix,jx,npz),      &
                      Cosp_rad(n)%block(nb)%daytime         (ix,jx)          )
      ! initial values
           Cosp_rad(n)%block(nb)%tau_stoch        = 0.
           Cosp_rad(n)%block(nb)%lwem_stoch       = 0.
           Cosp_rad(n)%block(nb)%stoch_cloud_type = 0.
           Cosp_rad(n)%block(nb)%stoch_conc_drop  = 0.
           Cosp_rad(n)%block(nb)%stoch_conc_ice   = 0.
           Cosp_rad(n)%block(nb)%stoch_size_drop  = 0.
           Cosp_rad(n)%block(nb)%stoch_size_ice   = 0.
           Cosp_rad(n)%block(nb)%mr_ozone         = 0.
           Cosp_rad(n)%block(nb)%daytime          = 0.
         end do
       endif
     end do

 end subroutine alloc_cosp_from_rad_type



 subroutine dealloc_cosp_from_rad_type (Cosp_rad, Exch_ctrl)
   type (cosp_from_rad_type), intent(inout) :: Cosp_rad(:)
   type (exchange_control_type), intent(in) :: Exch_ctrl
!--- local variables
   integer :: n, nb

    do n = 1, size(Cosp_rad,1)
      if (Exch_ctrl%do_cosp .or. Exch_ctrl%do_modis_yim) then
        do nb = 1, size(Cosp_rad(n)%block,1)
          deallocate (Cosp_rad(n)%block(nb)%tau_stoch,        &
                      Cosp_rad(n)%block(nb)%lwem_stoch,       &
                      Cosp_rad(n)%block(nb)%stoch_cloud_type, &
                      Cosp_rad(n)%block(nb)%stoch_conc_drop,  &
                      Cosp_rad(n)%block(nb)%stoch_conc_ice,   &
                      Cosp_rad(n)%block(nb)%stoch_size_drop,  &
                      Cosp_rad(n)%block(nb)%stoch_size_ice,   &
                      Cosp_rad(n)%block(nb)%mr_ozone,         &
                      Cosp_rad(n)%block(nb)%daytime)
        end do
      endif
      deallocate(Cosp_rad(n)%block)
    end do

 end subroutine dealloc_cosp_from_rad_type



 subroutine alloc_radiation_flux_type (Rad_flux, nonzero_init, Atm_block)
   type (radiation_flux_type), intent(inout) :: Rad_flux(:)
   logical,                    intent(in)    :: nonzero_init
   type (block_control_type),  intent(in)    :: Atm_block
!--- local variables
   integer :: n, nb
   integer :: ix, jx, npz

    npz = Atm_block%npz
    do n = 1, size(Rad_flux,1)
      allocate (Rad_flux(n)%block(Atm_block%nblks))
      do nb = 1, Atm_block%nblks
        ix = Atm_block%ibe(nb) - Atm_block%ibs(nb) + 1
        jx = Atm_block%jbe(nb) - Atm_block%jbs(nb) + 1
        allocate (Rad_flux(n)%block(nb)%tdt_rad               (ix,jx,npz), &
                  Rad_flux(n)%block(nb)%tdt_lw                (ix,jx,npz), &
                  Rad_flux(n)%block(nb)%flux_sw               (ix,jx),     &
                  Rad_flux(n)%block(nb)%flux_sw_dir           (ix,jx),     &
                  Rad_flux(n)%block(nb)%flux_sw_dif           (ix,jx),     &
                  Rad_flux(n)%block(nb)%flux_sw_down_vis_dir  (ix,jx),     &
                  Rad_flux(n)%block(nb)%flux_sw_down_vis_dif  (ix,jx),     &
                  Rad_flux(n)%block(nb)%flux_sw_down_total_dir(ix,jx),     &
                  Rad_flux(n)%block(nb)%flux_sw_down_total_dif(ix,jx),     &
                  Rad_flux(n)%block(nb)%flux_sw_vis           (ix,jx),     &
                  Rad_flux(n)%block(nb)%flux_sw_vis_dir       (ix,jx),     &
                  Rad_flux(n)%block(nb)%flux_sw_vis_dif       (ix,jx),     &
                  Rad_flux(n)%block(nb)%flux_lw               (ix,jx),     &
                  Rad_flux(n)%block(nb)%coszen                (ix,jx),     &
                  Rad_flux(n)%block(nb)%extinction            (ix,jx,npz)  )
      ! initial values - only used at t=0 when there is no restart
        if (nonzero_init) then
          Rad_flux(n)%block(nb)%tdt_rad                = 0.0
          Rad_flux(n)%block(nb)%tdt_lw                 = 0.0
          Rad_flux(n)%block(nb)%flux_sw                = 50.0
          Rad_flux(n)%block(nb)%flux_sw_dir            = 25.0
          Rad_flux(n)%block(nb)%flux_sw_dif            = 25.0
          Rad_flux(n)%block(nb)%flux_sw_down_vis_dir   = 12.5 * 1.1
          Rad_flux(n)%block(nb)%flux_sw_down_vis_dif   = 12.5 * 1.1
          Rad_flux(n)%block(nb)%flux_sw_down_total_dir = 25.0 * 1.1
          Rad_flux(n)%block(nb)%flux_sw_down_total_dif = 25.0 * 1.1
          Rad_flux(n)%block(nb)%flux_sw_vis            = 25.0
          Rad_flux(n)%block(nb)%flux_sw_vis_dir        = 12.5
          Rad_flux(n)%block(nb)%flux_sw_vis_dif        = 12.5
          Rad_flux(n)%block(nb)%flux_lw                = 50.0
          Rad_flux(n)%block(nb)%coszen                 = 0.50
          Rad_flux(n)%block(nb)%extinction             = 0.0
        else
          Rad_flux(n)%block(nb)%tdt_rad                = 0.0
          Rad_flux(n)%block(nb)%tdt_lw                 = 0.0
          Rad_flux(n)%block(nb)%flux_sw                = 0.0 !50.0
          Rad_flux(n)%block(nb)%flux_sw_dir            = 0.0 !25.0
          Rad_flux(n)%block(nb)%flux_sw_dif            = 0.0 !25.0
          Rad_flux(n)%block(nb)%flux_sw_down_vis_dir   = 0.0 !12.5 * 1.1
          Rad_flux(n)%block(nb)%flux_sw_down_vis_dif   = 0.0 !12.5 * 1.1
          Rad_flux(n)%block(nb)%flux_sw_down_total_dir = 0.0 !25.0 * 1.1
          Rad_flux(n)%block(nb)%flux_sw_down_total_dif = 0.0 !25.0 * 1.1
          Rad_flux(n)%block(nb)%flux_sw_vis            = 0.0 !25.0
          Rad_flux(n)%block(nb)%flux_sw_vis_dir        = 0.0 !12.5
          Rad_flux(n)%block(nb)%flux_sw_vis_dif        = 0.0 !12.5
          Rad_flux(n)%block(nb)%flux_lw                = 0.0 !50.0
          Rad_flux(n)%block(nb)%coszen                 = 0.0 !0.50
          Rad_flux(n)%block(nb)%extinction             = 0.0 !0.0
        endif
      end do
    end do

 end subroutine alloc_radiation_flux_type



 subroutine dealloc_radiation_flux_type (Rad_flux)
   type (radiation_flux_type), intent(inout) :: Rad_flux(:)
!--- local variables
   integer :: n, nb
!---------------------------------------------------------------------
!    deallocate the variables
!---------------------------------------------------------------------
      do n = 1, size(Rad_flux,1)
        do nb = 1, size(Rad_flux(n)%block,1)
          deallocate (Rad_flux(n)%block(nb)%tdt_rad,                &
                      Rad_flux(n)%block(nb)%tdt_lw,                 &
                      Rad_flux(n)%block(nb)%flux_sw,                &
                      Rad_flux(n)%block(nb)%flux_sw_dir,            &
                      Rad_flux(n)%block(nb)%flux_sw_dif,            &
                      Rad_flux(n)%block(nb)%flux_sw_down_vis_dir,   &
                      Rad_flux(n)%block(nb)%flux_sw_down_vis_dif,   &
                      Rad_flux(n)%block(nb)%flux_sw_down_total_dir, &
                      Rad_flux(n)%block(nb)%flux_sw_down_total_dif, &
                      Rad_flux(n)%block(nb)%flux_sw_vis,            &
                      Rad_flux(n)%block(nb)%flux_sw_vis_dir,        &
                      Rad_flux(n)%block(nb)%flux_sw_vis_dif,        &
                      Rad_flux(n)%block(nb)%flux_lw,                &
                      Rad_flux(n)%block(nb)%coszen,                 &
                      Rad_flux(n)%block(nb)%extinction              )
        enddo
        deallocate (Rad_flux(n)%block)
      end do
 end subroutine dealloc_radiation_flux_type


 subroutine alloc_clouds_from_moist_type (Moist_clouds, Exch_ctrl, Atm_block)
   type (clouds_from_moist_type), intent(inout) :: Moist_clouds(:)
   type (exchange_control_type),  intent(in)    :: Exch_ctrl
   type (block_control_type),     intent(in)    :: Atm_block
!--- local variables
   integer :: n, nb, npz
   integer :: ix, jx

   npz = Atm_block%npz
!-----------------------------------------------------------------------
!    allocate derived-type that stores cloud properties
!    return from moist processes
!-----------------------------------------------------------------------
    do n=1,size(Moist_clouds,1)
      allocate(Moist_clouds(n)%block(Atm_block%nblks))
      do nb = 1, Atm_block%nblks
        ix = Atm_block%ibe(nb) - Atm_block%ibs(nb) + 1
        jx = Atm_block%jbe(nb) - Atm_block%jbs(nb) + 1
        allocate (Moist_clouds(n)%block(nb)%lsc_cloud_area     (ix,jx,npz), &
                  Moist_clouds(n)%block(nb)%lsc_liquid         (ix,jx,npz), &
                  Moist_clouds(n)%block(nb)%lsc_ice            (ix,jx,npz), &
                  Moist_clouds(n)%block(nb)%lsc_droplet_number (ix,jx,npz), &
                  Moist_clouds(n)%block(nb)%lsc_ice_number     (ix,jx,npz), &
                  Moist_clouds(n)%block(nb)%lsc_rain           (ix,jx,npz), &
                  Moist_clouds(n)%block(nb)%lsc_snow           (ix,jx,npz), &
                  Moist_clouds(n)%block(nb)%lsc_rain_size      (ix,jx,npz), &
                  Moist_clouds(n)%block(nb)%lsc_snow_size      (ix,jx,npz)  )
        Moist_clouds(n)%block(nb)%lsc_cloud_area      = -99.
        Moist_clouds(n)%block(nb)%lsc_liquid          = -99.
        Moist_clouds(n)%block(nb)%lsc_ice             = -99.
        Moist_clouds(n)%block(nb)%lsc_droplet_number  = -99.
        Moist_clouds(n)%block(nb)%lsc_ice_number      = -99.
        Moist_clouds(n)%block(nb)%lsc_rain = 0.
        Moist_clouds(n)%block(nb)%lsc_snow = 0.
        Moist_clouds(n)%block(nb)%lsc_rain_size = 0.
        Moist_clouds(n)%block(nb)%lsc_snow_size = 0.

        if (Exch_ctrl%doing_donner) then
           allocate (Moist_clouds(n)%block(nb)%cell_cld_frac       (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%cell_liq_amt        (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%cell_liq_size       (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%cell_ice_amt        (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%cell_ice_size       (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%cell_droplet_number (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%meso_cld_frac       (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%meso_liq_amt        (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%meso_liq_size       (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%meso_ice_amt        (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%meso_ice_size       (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%meso_droplet_number (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%nsum_out            (ix,jx)      )
           Moist_clouds(n)%block(nb)%cell_cld_frac = 0.
           Moist_clouds(n)%block(nb)%cell_liq_amt  = 0.
           Moist_clouds(n)%block(nb)%cell_liq_size = 0.
           Moist_clouds(n)%block(nb)%cell_ice_amt  = 0.
           Moist_clouds(n)%block(nb)%cell_ice_size = 0.
           Moist_clouds(n)%block(nb)%cell_droplet_number = 0.
           Moist_clouds(n)%block(nb)%meso_cld_frac = 0.
           Moist_clouds(n)%block(nb)%meso_liq_amt  = 0.
           Moist_clouds(n)%block(nb)%meso_liq_size = 0.
           Moist_clouds(n)%block(nb)%meso_ice_amt  = 0.
           Moist_clouds(n)%block(nb)%meso_ice_size = 0.
           Moist_clouds(n)%block(nb)%meso_droplet_number = 0.
           Moist_clouds(n)%block(nb)%nsum_out = 1
        endif

        if (Exch_ctrl%doing_uw_conv) then
           allocate (Moist_clouds(n)%block(nb)%shallow_cloud_area     (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%shallow_liquid         (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%shallow_ice            (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%shallow_droplet_number (ix,jx,npz), &
                     Moist_clouds(n)%block(nb)%shallow_ice_number     (ix,jx,npz)  )
           Moist_clouds(n)%block(nb)%shallow_cloud_area      = 0.
           Moist_clouds(n)%block(nb)%shallow_liquid          = 0.
           Moist_clouds(n)%block(nb)%shallow_ice             = 0.
           Moist_clouds(n)%block(nb)%shallow_droplet_number  = 0.
           Moist_clouds(n)%block(nb)%shallow_ice_number      = 0.
        endif
      enddo
    enddo

 end subroutine alloc_clouds_from_moist_type



 subroutine dealloc_clouds_from_moist_type (Moist_clouds, Exch_ctrl)
   type (clouds_from_moist_type), intent(inout) :: Moist_clouds(:)
   type (exchange_control_type),  intent(in)    :: Exch_ctrl
!--- local variables
   integer :: n, nb
!--------------------------------------------------------------------
!    deallocate variables
!--------------------------------------------------------------------
    do n=1,size(Moist_clouds,1)
      do nb = 1, size(Moist_clouds(n)%block,1)
        deallocate (Moist_clouds(n)%block(nb)%lsc_cloud_area,     &
                    Moist_clouds(n)%block(nb)%lsc_liquid,         &
                    Moist_clouds(n)%block(nb)%lsc_ice,            &
                    Moist_clouds(n)%block(nb)%lsc_droplet_number, &
                    Moist_clouds(n)%block(nb)%lsc_ice_number,     &
                    Moist_clouds(n)%block(nb)%lsc_snow,           &
                    Moist_clouds(n)%block(nb)%lsc_rain,           &
                    Moist_clouds(n)%block(nb)%lsc_snow_size,      &
                    Moist_clouds(n)%block(nb)%lsc_rain_size       )
        if (Exch_ctrl%doing_donner) then
          deallocate (Moist_clouds(n)%block(nb)%cell_cld_frac,       &
                      Moist_clouds(n)%block(nb)%cell_liq_amt,        &
                      Moist_clouds(n)%block(nb)%cell_liq_size,       &
                      Moist_clouds(n)%block(nb)%cell_ice_amt,        &
                      Moist_clouds(n)%block(nb)%cell_ice_size,       &
                      Moist_clouds(n)%block(nb)%cell_droplet_number, &
                      Moist_clouds(n)%block(nb)%meso_cld_frac,       &
                      Moist_clouds(n)%block(nb)%meso_liq_amt,        &
                      Moist_clouds(n)%block(nb)%meso_liq_size,       &
                      Moist_clouds(n)%block(nb)%meso_ice_amt,        &
                      Moist_clouds(n)%block(nb)%meso_ice_size,       &
                      Moist_clouds(n)%block(nb)%meso_droplet_number, &
                      Moist_clouds(n)%block(nb)%nsum_out             )
        endif
        if (Exch_ctrl%doing_uw_conv) then
          deallocate (Moist_clouds(n)%block(nb)%shallow_cloud_area,     &
                      Moist_clouds(n)%block(nb)%shallow_liquid,         &
                      Moist_clouds(n)%block(nb)%shallow_ice,            &
                      Moist_clouds(n)%block(nb)%shallow_droplet_number, &
                      Moist_clouds(n)%block(nb)%shallow_ice_number      )
        endif
      enddo
      deallocate (Moist_clouds(n)%block)
    enddo

 end subroutine dealloc_clouds_from_moist_type




! update radiation flux states
 subroutine exch_rad_phys_state (Moist_clouds, Cosp_rad, Rad_flux, Exch_ctrl)
   type (clouds_from_moist_type), intent(inout) :: Moist_clouds(2)
   type (cosp_from_rad_type),     intent(inout) :: Cosp_rad(2)
   type (radiation_flux_type),    intent(inout) :: Rad_flux(2)
   type (exchange_control_type),  intent(in)    :: Exch_ctrl
!--- local variables
   integer :: nb

    Cosp_rad(1)%control%step_to_call_cosp = Cosp_rad(2)%control%step_to_call_cosp
    Rad_flux(1)%control%do_rad = Rad_flux(2)%control%do_rad
!$OMP parallel do default(shared) private(nb)
    do nb = 1, size(Rad_flux(1)%block)
       Rad_flux(1)%block(nb)%tdt_rad                = Rad_flux(2)%block(nb)%tdt_rad
       Rad_flux(1)%block(nb)%tdt_lw                 = Rad_flux(2)%block(nb)%tdt_lw
       Rad_flux(1)%block(nb)%flux_sw                = Rad_flux(2)%block(nb)%flux_sw
       Rad_flux(1)%block(nb)%flux_sw_dir            = Rad_flux(2)%block(nb)%flux_sw_dir
       Rad_flux(1)%block(nb)%flux_sw_dif            = Rad_flux(2)%block(nb)%flux_sw_dif
       Rad_flux(1)%block(nb)%flux_sw_down_vis_dir   = Rad_flux(2)%block(nb)%flux_sw_down_vis_dir
       Rad_flux(1)%block(nb)%flux_sw_down_vis_dif   = Rad_flux(2)%block(nb)%flux_sw_down_vis_dif
       Rad_flux(1)%block(nb)%flux_sw_down_total_dir = Rad_flux(2)%block(nb)%flux_sw_down_total_dir
       Rad_flux(1)%block(nb)%flux_sw_down_total_dif = Rad_flux(2)%block(nb)%flux_sw_down_total_dif
       Rad_flux(1)%block(nb)%flux_sw_vis            = Rad_flux(2)%block(nb)%flux_sw_vis
       Rad_flux(1)%block(nb)%flux_sw_vis_dir        = Rad_flux(2)%block(nb)%flux_sw_vis_dir
       Rad_flux(1)%block(nb)%flux_sw_vis_dif        = Rad_flux(2)%block(nb)%flux_sw_vis_dif
       Rad_flux(1)%block(nb)%flux_lw                = Rad_flux(2)%block(nb)%flux_lw
       Rad_flux(1)%block(nb)%coszen                 = Rad_flux(2)%block(nb)%coszen
       Rad_flux(1)%block(nb)%extinction             = Rad_flux(2)%block(nb)%extinction
       if (Exch_ctrl%do_cosp .or. Exch_ctrl%do_modis_yim) then
         Cosp_rad(1)%block(nb)%tau_stoch        = Cosp_rad(2)%block(nb)%tau_stoch
         Cosp_rad(1)%block(nb)%lwem_stoch       = Cosp_rad(2)%block(nb)%lwem_stoch
         Cosp_rad(1)%block(nb)%stoch_cloud_type = Cosp_rad(2)%block(nb)%stoch_cloud_type
         Cosp_rad(1)%block(nb)%stoch_conc_drop  = Cosp_rad(2)%block(nb)%stoch_conc_drop
         Cosp_rad(1)%block(nb)%stoch_conc_ice   = Cosp_rad(2)%block(nb)%stoch_conc_ice
         Cosp_rad(1)%block(nb)%stoch_size_drop  = Cosp_rad(2)%block(nb)%stoch_size_drop
         Cosp_rad(1)%block(nb)%stoch_size_ice   = Cosp_rad(2)%block(nb)%stoch_size_ice 
         Cosp_rad(1)%block(nb)%mr_ozone         = Cosp_rad(2)%block(nb)%mr_ozone
         Cosp_rad(1)%block(nb)%daytime          = Cosp_rad(2)%block(nb)%daytime 
       endif

       Moist_clouds(2)%block(nb)%lsc_cloud_area     = Moist_clouds(1)%block(nb)%lsc_cloud_area
       Moist_clouds(2)%block(nb)%lsc_liquid         = Moist_clouds(1)%block(nb)%lsc_liquid   
       Moist_clouds(2)%block(nb)%lsc_ice            = Moist_clouds(1)%block(nb)%lsc_ice     
       Moist_clouds(2)%block(nb)%lsc_droplet_number = Moist_clouds(1)%block(nb)%lsc_droplet_number
       Moist_clouds(2)%block(nb)%lsc_ice_number     = Moist_clouds(1)%block(nb)%lsc_ice_number   
       Moist_clouds(2)%block(nb)%lsc_rain           = Moist_clouds(1)%block(nb)%lsc_rain        
       Moist_clouds(2)%block(nb)%lsc_snow           = Moist_clouds(1)%block(nb)%lsc_snow       
       Moist_clouds(2)%block(nb)%lsc_rain_size      = Moist_clouds(1)%block(nb)%lsc_rain_size 
       Moist_clouds(2)%block(nb)%lsc_snow_size      = Moist_clouds(1)%block(nb)%lsc_snow_size

       if (Exch_ctrl%doing_donner) then
          Moist_clouds(2)%block(nb)%cell_cld_frac       = Moist_clouds(1)%block(nb)%cell_cld_frac
          Moist_clouds(2)%block(nb)%cell_liq_amt        = Moist_clouds(1)%block(nb)%cell_liq_amt
          Moist_clouds(2)%block(nb)%cell_liq_size       = Moist_clouds(1)%block(nb)%cell_liq_size
          Moist_clouds(2)%block(nb)%cell_ice_amt        = Moist_clouds(1)%block(nb)%cell_ice_amt
          Moist_clouds(2)%block(nb)%cell_ice_size       = Moist_clouds(1)%block(nb)%cell_ice_size
          Moist_clouds(2)%block(nb)%cell_droplet_number = Moist_clouds(1)%block(nb)%cell_droplet_number
          Moist_clouds(2)%block(nb)%meso_cld_frac       = Moist_clouds(1)%block(nb)%meso_cld_frac     
          Moist_clouds(2)%block(nb)%meso_liq_amt        = Moist_clouds(1)%block(nb)%meso_liq_amt     
          Moist_clouds(2)%block(nb)%meso_liq_size       = Moist_clouds(1)%block(nb)%meso_liq_size   
          Moist_clouds(2)%block(nb)%meso_ice_amt        = Moist_clouds(1)%block(nb)%meso_ice_amt   
          Moist_clouds(2)%block(nb)%meso_ice_size       = Moist_clouds(1)%block(nb)%meso_ice_size 
          Moist_clouds(2)%block(nb)%meso_droplet_number = Moist_clouds(1)%block(nb)%meso_droplet_number 
          Moist_clouds(2)%block(nb)%nsum_out            = Moist_clouds(1)%block(nb)%nsum_out     
       endif

        if (Exch_ctrl%doing_uw_conv) then
          Moist_clouds(2)%block(nb)%shallow_cloud_area     = Moist_clouds(1)%block(nb)%shallow_cloud_area   
          Moist_clouds(2)%block(nb)%shallow_liquid         = Moist_clouds(1)%block(nb)%shallow_liquid      
          Moist_clouds(2)%block(nb)%shallow_ice            = Moist_clouds(1)%block(nb)%shallow_ice        
          Moist_clouds(2)%block(nb)%shallow_droplet_number = Moist_clouds(1)%block(nb)%shallow_droplet_number
          Moist_clouds(2)%block(nb)%shallow_ice_number     = Moist_clouds(1)%block(nb)%shallow_ice_number
        endif
     enddo

 end subroutine exch_rad_phys_state

end module physics_radiation_exch_mod
