
module vert_turb_driver_mod

!-----------------------------------------------------------------------
!
!       driver for compuing vertical diffusion coefficients
!
!         choose either:
!              1) mellor-yamada 2.5 (with tke)
!              2) non-local K scheme
!
!-----------------------------------------------------------------------
!---------------- modules ---------------------


use      my25_turb_mod, only: my25_turb_init, my25_turb_end,  &
                              my25_turb, tke_surf, tke

use    diffusivity_mod, only: diffusivity

use   shallow_conv_mod, only: shallow_conv_init, shallow_conv

use   diag_manager_mod, only: register_diag_field, send_data

use   time_manager_mod, only: time_type

use      constants_mod, only: grav, rdgas, rvgas, cp, kappa, p00, hlv

use      utilities_mod, only: error_mesg, open_file, file_exist, &
                              get_my_pe, close_file, FATAL

implicit none
private

!---------------- interfaces ---------------------

public   vert_turb_driver_init, vert_turb_driver_end, vert_turb_driver


!-----------------------------------------------------------------------
!--------------------- version number ----------------------------------

character(len=128) :: version = '$Id: vert_turb_driver.F90,v 1.2 2000/08/04 19:04:10 fms Exp $'
character(len=128) :: tag = '$Name: bombay $'

!-----------------------------------------------------------------------

 real, parameter :: p00inv = 1./p00
 real, parameter :: d622   = rdgas/rvgas
 real, parameter :: d378   = 1.-d622
 real, parameter :: d608   = d378/d622

!---------------- private data -------------------

 logical :: do_init = .true.

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

 logical :: do_shallow_conv  = .false.
 logical :: do_mellor_yamada = .true.

 namelist /vert_turb_driver_nml/ do_shallow_conv, do_mellor_yamada

!-------------------- diagnostics fields -------------------------------

integer :: id_tke,    id_lscale, id_lscale_0, id_z_pbl,  &
           id_diff_t, id_diff_m, id_diff_sc

real :: missing_value = -999.

character(len=9) :: mod_name = 'vert_turb'

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine vert_turb_driver (is, js, Time, dt, frac_land,      &
                             p_half, p_full, z_half, z_full,   &
                             u_star, b_star, rough,            &
                             u, v, t, q, diff_t, diff_m, gust, &
                             mask, kbot)

!-----------------------------------------------------------------------
integer,         intent(in)         :: is, js
type(time_type), intent(in)         :: Time
   real,         intent(in)         :: dt
   real, intent(in), dimension(:,:) :: frac_land,   &
                                       u_star, b_star, rough
   real, intent(in), dimension(:,:,:) :: p_half, p_full, &
                                         z_half, z_full, &
                                         u, v, t, q
   real, intent(out),   dimension(:,:,:) :: diff_t, diff_m
   real, intent(out),   dimension(:,:)   :: gust
   real, intent(in),optional, dimension(:,:,:) :: mask
integer, intent(in),optional, dimension(:,:) :: kbot
!-----------------------------------------------------------------------
real   , dimension(size(t,1),size(t,2),size(t,3))   :: ape, thv
logical, dimension(size(t,1),size(t,2),size(t,3)+1) :: lmask
real   , dimension(size(t,1),size(t,2),size(t,3)+1) :: el, diag3
real   , dimension(size(t,1),size(t,2))             :: el0, z_pbl
real   , dimension(size(diff_t,1),size(diff_t,2), &
                                  size(diff_t,3))   :: diff_sc
integer :: ie, je, nlev
logical :: used
!-----------------------------------------------------------------------
!----------------------- vertical turbulence ---------------------------
!-----------------------------------------------------------------------

      if (do_init)  call error_mesg  &
                     ('vert_turb_driver in vert_turb_driver_mod',  &
                      'initialization has not been called', FATAL)

     nlev = size(p_full,3)
     ie = is + size(p_full,1) - 1
     je = js + size(p_full,2) - 1

!---------------------------
 if (do_mellor_yamada) then
!---------------------------

!    ----- virtual temp ----------
     ape(:,:,:)=(p_full(:,:,:)*p00inv)**(-kappa)
     thv(:,:,:)=t(:,:,:)*(q(:,:,:)*d608+1.0)*ape(:,:,:)
     if (present(mask)) where (mask < 0.5) thv = 200.

!    --------------------- update tke-----------------------------------
!    ---- compute surface tke --------
!    ---- compute tke, master length scale (el0),  -------------
!    ---- length scale (el), and vert mix coeffs (diff_t,diff_m) ----

     call tke_surf  (u_star, tke(is:ie,js:je,:), kbot=kbot)



     if ( id_z_pbl > 0 ) then
     !------ compute pbl depth from k_profile if diagnostic needed -----
     call my25_turb (dt, frac_land, p_half, p_full, thv, u, v,  &
                     z_half, z_full, rough, tke(is:ie,js:je,:), &
                     el0, el, diff_m, diff_t, &
                     mask=mask, kbot=kbot, &
                     ustar=u_star,bstar=b_star,h=z_pbl)
     else
     call my25_turb (dt, frac_land, p_half, p_full, thv, u, v,  &
                     z_half, z_full, rough, tke(is:ie,js:je,:), &
                     el0, el, diff_m, diff_t, &
                     mask=mask, kbot=kbot)
     end if

!---------------------------
 else
!---------------------------
!------------------- non-local K scheme --------------


    call diffusivity ( t, q, u, v, z_full, z_half,   &
                       u_star, b_star, z_pbl, diff_m, diff_t, &
                       kbot = kbot)

!---------------------------
 endif
!-----------------------------------------------------------------------
!------------------ shallow convection ???? ----------------------------

   if (do_shallow_conv) then
        call shallow_conv (t, q, p_full, p_half, diff_sc, kbot)
        diff_t = diff_t + diff_sc
   endif

!-----------------------------------------------------------------------
!------------- define gustiness ------------

     gust = 1.0

!-----------------------------------------------------------------------
!------------------------ diagnostics section --------------------------

if (do_mellor_yamada) then

!     --- set up local mask for fields with surface data ---
      if ( present(mask) ) then
         lmask(:,:,1)        = .true.
         lmask(:,:,2:nlev+1) = mask(:,:,1:nlev) > 0.5
      else
         lmask = .true.
      endif

!------- tke --------------------------------
      if ( id_tke > 0 ) then
         used = send_data ( id_tke, tke(is:ie,js:je,:), Time,  &
                            is, js, 1, mask=lmask )
      endif

!------- length scale (at half levels) ------
      if ( id_lscale > 0 ) then
         used = send_data ( id_lscale, el, Time, is, js, 1,  &
                            mask=lmask )
      endif

!------- master length scale -------
      if ( id_lscale_0 > 0 ) then
         used = send_data ( id_lscale_0, el0, Time, is, js )
      endif

end if

!------- boundary layer depth -------
      if ( id_z_pbl > 0 ) then
         used = send_data ( id_z_pbl, z_pbl, Time, is, js )
      endif


!------- output diffusion coefficients ---------

  if ( id_diff_t > 0 .or. id_diff_m > 0 .or. id_diff_sc > 0 ) then
!       --- set up local mask for fields without surface data ---
        if (present(mask)) then
            lmask(:,:,1:nlev) = mask(:,:,1:nlev) > 0.5
            lmask(:,:,nlev+1) = .false.
        else
            lmask(:,:,1:nlev) = .true.
            lmask(:,:,nlev+1) = .false.
        endif
!       -- dummy data at surface --
        diag3(:,:,nlev+1)=0.0
  endif

!------- diffusion coefficient for heat/moisture -------
   if ( id_diff_t > 0 ) then
      diag3(:,:,1:nlev) = diff_t(:,:,1:nlev)
      used = send_data ( id_diff_t, diag3, Time, is, js, 1, mask=lmask )
   endif

!------- diffusion coefficient for momentum -------
   if ( id_diff_m > 0 ) then
      diag3(:,:,1:nlev) = diff_m(:,:,1:nlev)
      used = send_data ( id_diff_m, diag3, Time, is, js, 1, mask=lmask )
   endif

!------- diffusion coefficient for shallow conv -------
 if (do_shallow_conv) then
   if ( id_diff_sc > 0 ) then
      diag3(:,:,1:nlev) = diff_sc(:,:,1:nlev)
      used = send_data ( id_diff_sc, diag3, Time, is, js, 1, mask=lmask)
   endif
 endif

!-----------------------------------------------------------------------

end subroutine vert_turb_driver

!#######################################################################

subroutine vert_turb_driver_init (id, jd, kd, axes, Time)

!-----------------------------------------------------------------------
   integer,         intent(in) :: id, jd, kd, axes(4)
   type(time_type), intent(in) :: Time
!-----------------------------------------------------------------------
   integer, dimension(3) :: full = (/1,2,3/), half = (/1,2,4/)
   integer :: ierr, unit, io

      if (.not.do_init)  &
          call error_mesg  &
                   ('vert_turb_driver_init in vert_turb_driver_mod',  &
                    'attempting to call initialization twice', FATAL)

!-----------------------------------------------------------------------
!--------------- read namelist ------------------

      if (file_exist('input.nml')) then
         unit = open_file (file='input.nml', action='read')
         io=1; do while (io .ne. 0)
            read  (unit, nml=vert_turb_driver_nml, iostat=io, end=10)
         enddo
  10     call close_file (unit)
      endif

!---------- output namelist --------------------------------------------

      unit = open_file (file='logfile.out', action='append')
      if ( get_my_pe() == 0 ) then
           write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
           write (unit,nml=vert_turb_driver_nml)
      endif
      call close_file (unit)

!-----------------------------------------------------------------------

      if (do_mellor_yamada) call my25_turb_init (id, jd, kd)

      if (do_shallow_conv)  &
         call shallow_conv_init (kd, p00, hlv, cp, rdgas, rvgas,   &
                                 kappa, grav, d622, d378           )

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----

if (do_mellor_yamada) then

   id_tke = &
   register_diag_field ( mod_name, 'tke', axes(half), Time,      &
                        'turbulent kinetic energy',  'm2/s2'   , &
                        missing_value=missing_value               )

   id_lscale = &
   register_diag_field ( mod_name, 'lscale', axes(half), Time,    &
                        'turbulent length scale',  'm'   ,        &
                        missing_value=missing_value               )

   id_lscale_0 = &
   register_diag_field ( mod_name, 'lscale_0', axes(1:2), Time,   &
                        'master length scale',  'm'               )
endif

   id_z_pbl = &
   register_diag_field ( mod_name, 'z_pbl', axes(1:2), Time,       &
                        'depth of planetary boundary layer',  'm'  )

   id_diff_t = &
   register_diag_field ( mod_name, 'diff_t', axes(half), Time,    &
                        'vert diff coeff for temp',  'm2/s'   ,   &
                        missing_value=missing_value               )

   id_diff_m = &
   register_diag_field ( mod_name, 'diff_m', axes(half), Time,      &
                        'vert diff coeff for momentum',  'm2/s'   , &
                        missing_value=missing_value               )

if (do_shallow_conv) then

   id_diff_sc = &
   register_diag_field ( mod_name, 'diff_sc', axes(half), Time,      &
                        'vert diff coeff for shallow conv', 'm2/s' , &
                        missing_value=missing_value               )
endif

!-----------------------------------------------------------------------

   do_init=.false.

!-----------------------------------------------------------------------

end subroutine vert_turb_driver_init


!#######################################################################

subroutine vert_turb_driver_end

!-----------------------------------------------------------------------
      if (do_mellor_yamada) call my25_turb_end
!-----------------------------------------------------------------------

end subroutine vert_turb_driver_end

!#######################################################################

end module vert_turb_driver_mod

