
module vert_diff_driver_mod

!-----------------------------------------------------------------------
!   module performs vertical diffusion of atmospheric variables
!-----------------------------------------------------------------------

use    vert_diff_mod, only:  surf_diff_type,     &
                             gcm_vert_diff_init, &
                             gcm_vert_diff_down, &
                             gcm_vert_diff_up,   &
                             gcm_vert_diff_end

use  strat_cloud_mod, only:  add_strat_tend, subtract_strat_tend

use diag_manager_mod, only:  register_diag_field, send_data

use time_manager_mod, only:  time_type

use    utilities_mod, only:  file_exist, open_file, error_mesg,  &
                             check_nml_error, FATAL, get_my_pe,  &
                             close_file

use    constants_mod, only:  cp, grav

!-----------------------------------------------------------------------

implicit none
private

public :: vert_diff_driver_down, vert_diff_driver_up,  &
          vert_diff_driver_init, vert_diff_driver_end
public :: surf_diff_type


!-----------------------------------------------------------------------
!---- namelist ----

integer :: ntp = 0

!    ntp  = number of fully prognostic tracers
!           tracers (1:ntp) will be vertically diffused

namelist /vert_diff_driver_nml/  ntp

!-----------------------------------------------------------------------
!-------------------- diagnostics fields -------------------------------

integer :: id_tdt_vdif, id_qdt_vdif, id_udt_vdif, id_vdt_vdif,  &
           id_sens_vdif, id_evap_vdif

real :: missing_value = -999.

character(len=9), parameter :: mod_name = 'vert_diff'

!-----------------------------------------------------------------------
!---- version number ----

character(len=128) :: version = '$Id: vert_diff_driver.F90,v 1.2 2000/07/28 20:16:47 fms Exp $'
character(len=128) :: tag = '$Name: eugene $'

logical :: do_init = .true.


contains

!#######################################################################

 subroutine vert_diff_driver_down (is, js, Time, delt, p_half, z_full, &
                                   diff_mom, diff_heat,                &
                                   u, v, t, q, trs,                    &
                                   dtau_dv, tau_x, tau_y,              &
                                   dt_u, dt_v, dt_t, dt_q, dt_trs,     &
                                   Surf_diff,  mask, kbot              )

integer, intent(in)                     :: is, js
type(time_type),   intent(in)           :: Time
real, intent(in)                        :: delt
real, intent(in)   , dimension(:,:,:)   :: p_half, z_full,  &
                                           diff_mom, diff_heat
real, intent(in),    dimension(:,:,:)   :: u, v, t, q
real, intent(in),    dimension(:,:,:,:) :: trs
real, intent(in),    dimension(:,:)     :: dtau_dv

real, intent(inout), dimension(:,:)     :: tau_x, tau_y
real, intent(inout), dimension(:,:,:)   :: dt_u, dt_v, dt_t, dt_q
real, intent(inout), dimension(:,:,:,:) :: dt_trs

type(surf_diff_type), intent(inout)     :: Surf_diff

real   , intent(in), dimension(:,:,:), optional :: mask
integer, intent(in), dimension(:,:),   optional :: kbot

real, dimension(size(trs,1),size(trs,2),1:ntp) :: flux_trs
real, dimension(size(t,1),size(t,2),size(t,3)) :: tt, dpg
integer :: k
logical :: used
real, dimension(size(t,1),size(t,2)) :: diag2


!-----------------------------------------------------------------------

  if (do_init) call error_mesg       &
                  ('vert_diff_driver_mod',  &
                   'vert_diff_driver_init must be called first', FATAL)

!-----------------------------------------------------------------------

  if (size(trs,4) < ntp) call error_mesg ('vert_diff_driver_mod',  &
              'Number of tracers .lt. namelist value of ntp', FATAL)

!-----------------------------------------------------------------------
!---- save temperature tendency for stratiform cloud scheme ----
!   note if strat scheme is not operating this routine does nothing

    call subtract_strat_tend (is, is+size(dt_t,1)-1,     &
                              js, js+size(dt_t,2)-1, dt_t)

!-----------------------------------------------------------------------
!---- to do diagnostics on dt_t, dt_q, dt_u, and dt_v at this point add 
!-----in the negative value of the field.  Note that the multiplication
!---- by 2 is necessary to get the correct value because sum_diag_phys
!---- is called twice for this single field
!-----------------------------------------------------------------------
!------- diagnostics for dt/dt_diff -------

    if ( id_tdt_vdif > 0 ) then
       used = send_data ( id_tdt_vdif, -2.*dt_t, Time, is, js, 1, &
                           rmask=mask )
    endif
!------- diagnostics for dq/dt_diff -------
    if ( id_qdt_vdif > 0 ) then
       used = send_data ( id_qdt_vdif, -2.*dt_q, Time, is, js, 1, &
                          rmask=mask )
    endif

!------- diagnostics for uwnd_diff -------
    if ( id_udt_vdif > 0 ) then
       used = send_data ( id_udt_vdif, -2.*dt_u, Time, is, js, 1, &
                          rmask=mask )
    endif

!------- diagnostics for vwnd_diff -------
    if ( id_vdt_vdif > 0 ) then
       used = send_data ( id_vdt_vdif, -2.*dt_v, Time, is, js, 1, &
                          rmask=mask )
    endif

!------ preliminary calculations for vert integrals -----
    if ( id_sens_vdif > 0 .or. id_evap_vdif > 0 ) then
            do k = 1, size(p_half,3)-1
               dpg(:,:,k) = (p_half(:,:,k+1)-p_half(:,:,k))/grav
            enddo
            if (present(mask)) dpg = dpg*mask
    endif

!------- diagnostics for sens_diff -------
    if ( id_sens_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_t*dpg, 3 )
          used = send_data ( id_sens_vdif, -2.*cp*diag2, Time, is, js )
    endif

!------- diagnostics for evap_diff -------
    if ( id_evap_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_q*dpg, 3 )
          used = send_data ( id_evap_vdif, -2.*diag2, Time, is, js )
    endif

!-----------------------------------------------------------------------
!---- local temperature ----
!     (avoids divid by zero when computing virtual temp)

   tt = t
   if (present(mask)) where (mask < 0.5) tt = 200.

!-----------------------------------------------------------------------
!---- momentum diffusion ----
!---- heat/moisture diffusion (down only) ----
!---- tracer diffusion (no surface flux) ----

 flux_trs = 0.0

 call gcm_vert_diff_down (is, js, delt, u, v, tt, q, trs(:,:,:,1:ntp), &
                          diff_mom, diff_heat, p_half, z_full,         &
                          tau_x, tau_y, dtau_dv,  flux_trs,            &
                          dt_u, dt_v, dt_t, dt_q, dt_trs(:,:,:,1:ntp), &
                          Surf_diff,  kbot                             )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!---- to do diagnostics on dt_u, and dt_v at this point add 
!-----in the value of the field.  Note that the multiplication
!---- by 2 is necessary to get the correct value because sum_diag_phys
!---- is called twice for this single field
!-----------------------------------------------------------------------

!------- diagnostics for uwnd_diff -------
    if ( id_udt_vdif > 0 ) then
       used = send_data ( id_udt_vdif, 2.*dt_u, Time, is, js, 1, &
                          rmask=mask )
    endif

!------- diagnostics for vwnd_diff -------
    if ( id_vdt_vdif > 0 ) then
       used = send_data ( id_vdt_vdif, 2.*dt_v, Time, is, js, 1, &
                          rmask=mask )
    endif

!-----------------------------------------------------------------------

 end subroutine vert_diff_driver_down

!#######################################################################

 subroutine vert_diff_driver_up (is, js, Time, delt, p_half, &
                                 Surf_diff, dt_t, dt_q, mask, kbot)

 integer,           intent(in)            :: is, js
 type(time_type),   intent(in)            :: Time
 real,    intent(in)                      :: delt
 real,    intent(in),    dimension(:,:,:) :: p_half
 type(surf_diff_type),   intent(in)       :: Surf_diff
 real,    intent(inout), dimension(:,:,:) :: dt_t, dt_q
 real   , intent(in), dimension(:,:,:), optional :: mask
 integer, intent(in),    dimension(:,:), optional :: kbot

 integer :: k
 logical :: used
 real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: dpg
 real, dimension(size(p_half,1),size(p_half,2)) :: diag2

!-----------------------------------------------------------------------

    call gcm_vert_diff_up (is, js, delt, Surf_diff, dt_t, dt_q, kbot)

!-----------------------------------------------------------------------
!---- compute temperature tendency for stratiform cloud scheme ----
!   note if strat scheme is not operating this routine does nothing

    call add_strat_tend (is, is+size(dt_t,1)-1,     &
                         js, js+size(dt_t,2)-1, dt_t)

!-----------------------------------------------------------------------
!---- to do diagnostics on dt_t and dt_q at this point add in the
!---- the postive value of the field.  Note that the multiplication
!---- by 2 is necessary to get the correct value because sum_diag_phys
!---- is called twice for this single field
!-----------------------------------------------------------------------
!------- diagnostics for dt/dt_diff -------

    if ( id_tdt_vdif > 0 ) then
       used = send_data ( id_tdt_vdif, 2.*dt_t, Time, is, js, 1, &
                          rmask=mask )
    endif

!------- diagnostics for dq/dt_diff -------
    if ( id_qdt_vdif > 0 ) then
       used = send_data ( id_qdt_vdif, 2.*dt_q, Time, is, js, 1, &
                          rmask=mask )
    endif

!------ preliminary calculations for vert integrals -----
    if ( id_sens_vdif > 0 .or. id_evap_vdif > 0 ) then
            do k = 1, size(p_half,3)-1
               dpg(:,:,k) = (p_half(:,:,k+1)-p_half(:,:,k))/grav
            enddo
            if (present(mask)) dpg = dpg*mask
    endif

!------- diagnostics for sens_diff -------
    if ( id_sens_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_t*dpg, 3 )
          used = send_data ( id_sens_vdif, 2.*cp*diag2, Time, is, js )
    endif

!------- diagnostics for evap_diff -------
    if ( id_evap_vdif > 0 ) then
!         --- compute column changes ---
          diag2 = sum( dt_q*dpg, 3 )
          used = send_data ( id_evap_vdif, 2.*diag2, Time, is, js )
    endif


!-----------------------------------------------------------------------

 end subroutine vert_diff_driver_up

!#######################################################################

 subroutine vert_diff_driver_init ( Surf_diff, idim, jdim, kdim,  &
                                    axes, Time )

 type(surf_diff_type), intent(inout) :: Surf_diff
 integer             , intent(in)    :: idim, jdim, kdim, axes(4)
 type     (time_type), intent(in)    :: Time

 integer :: unit, io, ierr

!-----------------------------------------------------------------------
!------ read namelist ------

   if ( file_exist('input.nml')) then
      unit = open_file ('input.nml', action='read')
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=vert_diff_driver_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'vert_diff_driver_nml')
      enddo
 10   call close_file (unit)
   endif

!--------- write version number and namelist ------------------

   unit = open_file ('logfile.out', action='append')
   if ( get_my_pe() == 0 ) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit, nml=vert_diff_driver_nml)
   endif
   call close_file (unit)

!-------- initialize gcm vertical diffusion ------

   call gcm_vert_diff_init (Surf_diff, idim, jdim, kdim)

!-----------------------------------------------------------------------
!--------------- initialize diagnostic fields --------------------

   id_tdt_vdif = &
   register_diag_field ( mod_name, 'tdt_vdif', axes(1:3), Time, &
                        'Temperature tendency from vert diff',  &
                        'deg_K/s', missing_value=missing_value  )

   id_qdt_vdif = &
   register_diag_field ( mod_name, 'qdt_vdif', axes(1:3), Time, &
                        'Spec humidity tendency from vert diff',&
                        'kg/kg/s', missing_value=missing_value  )

   id_udt_vdif = &
   register_diag_field ( mod_name, 'udt_vdif', axes(1:3), Time, &
                        'Zonal wind tendency from vert diff',   &
                        'm/s2', missing_value=missing_value     )

   id_vdt_vdif = &
   register_diag_field ( mod_name, 'vdt_vdif', axes(1:3), Time,    &
                        'Meridional wind tendency from vert diff', &
                        'm/s2', missing_value=missing_value        )

   id_sens_vdif = &
   register_diag_field ( mod_name, 'sens_vdif', axes(1:2), Time,  &
                        'Integrated heat flux from vert diff',    &
                        'W/m2' )

   id_evap_vdif = &
   register_diag_field ( mod_name, 'evap_vdif', axes(1:2), Time,    &
                        'Integrated moisture flux from vert diff',  &
                        'kg/m2/s' )

!-----------------------------------------------------------------------

   do_init = .false.

!-----------------------------------------------------------------------

 end subroutine vert_diff_driver_init

!#######################################################################

 subroutine vert_diff_driver_end

   call gcm_vert_diff_end

 end subroutine vert_diff_driver_end

!#######################################################################

end module vert_diff_driver_mod

