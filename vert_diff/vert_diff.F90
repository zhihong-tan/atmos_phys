
module vert_diff_mod

!=======================================================================
!
!                         VERTICAL DIFFUSION MODULE
!
!      Routines for computing the tendencies due to vertical diffusion
!
!=======================================================================

use   constants_mod, only:  grav, rdgas, rvgas, cp

use   utilities_mod, only:  error_mesg, FATAL, file_exist, &
                            check_nml_error, open_file,    &
                            get_my_pe, close_file

implicit none
private


! public interfaces
!=======================================================================
public :: gcm_vert_diff_init,          &
          gcm_vert_diff_end,           &
          gcm_vert_diff,               &
          gcm_vert_diff_down,          &
          gcm_vert_diff_surf_down,     &
          gcm_vert_diff_surf_up,       &
          gcm_vert_diff_up,            &
          vert_diff,                   &
          surf_diff_type,              &
          switch_surf_diff_type_order, &
          alloc_surf_diff_type,        &
          dealloc_surf_diff_type
!=======================================================================

! form of interfaces
!=======================================================================



interface alloc_surf_diff_type
   module procedure alloc_surf_diff_type_1d, alloc_surf_diff_type_2d
end interface



type surf_diff_type
  integer :: array_order

  real, pointer, dimension(:,:) :: mu_delt_n,    &
                                   nu_n,         &
                                   e_n1,         &
                                   f_t_delt_n1,  &
                                   f_q_delt_n1,  &
                                   delta_t_n,    &
                                   delta_q_n

  real, pointer, dimension(:) :: x_mu_delt_n,    &
                                 x_nu_n,         &
                                 x_e_n1,         &
                                 x_f_t_delt_n1,  &
                                 x_f_q_delt_n1,  &
                                 x_delta_t_n,    &
                                 x_delta_q_n    

  real, pointer, dimension(:) :: x_e_t_n,        &
                                 x_e_q_n,        &
                                 x_f_t_delt_n,   &
                                 x_f_q_delt_n    
end type surf_diff_type


real,    allocatable, dimension(:,:,:) :: e_global, f_t_global, f_q_global 
integer, allocatable, dimension(:,:)   :: order

      
logical :: do_init = .true.

!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id: vert_diff.F90,v 1.3 2000/12/05 13:25:07 fms Exp $'
character(len=128) :: tag = '$Name: damascus $'

real, parameter :: d608 = (rvgas-rdgas)/rdgas

contains

!#######################################################################

subroutine gcm_vert_diff_init (Tri_surf, idim, jdim, kdim)

 type(surf_diff_type), intent(inout) :: Tri_surf
 integer,              intent(in)    :: idim, jdim, kdim

 integer :: unit

    unit = open_file ('logfile.out', action='append')
    if (get_my_pe() == 0) &
    write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
    call close_file (unit)

 if (do_init) then

    if (allocated(  order    ))    deallocate (  order    )
    if (allocated(  e_global ))    deallocate (  e_global )
    if (allocated(f_t_global ))    deallocate (f_t_global )
    if (allocated(f_q_global ))    deallocate (f_q_global )

    allocate(  order    (idim, jdim)        )
    allocate(  e_global (idim, jdim, kdim-1))
    allocate(f_t_global (idim, jdim, kdim-1))
    allocate(f_q_global (idim, jdim, kdim-1))

    order = 0
    do_init = .false.

 endif

 call alloc_surf_diff_type ( Tri_surf, idim, jdim )

end subroutine gcm_vert_diff_init

!#######################################################################

subroutine alloc_surf_diff_type_1d ( Tri_surf, idim )

type(surf_diff_type), intent(inout) :: Tri_surf
integer,              intent(in)    :: idim

    Tri_surf%array_order = 1

    allocate( Tri_surf%x_mu_delt_n   (idim) )
    allocate( Tri_surf%x_nu_n        (idim) )
    allocate( Tri_surf%x_e_n1        (idim) )
    allocate( Tri_surf%x_f_t_delt_n1 (idim) )
    allocate( Tri_surf%x_f_q_delt_n1 (idim) )
    allocate( Tri_surf%x_delta_t_n   (idim) )
    allocate( Tri_surf%x_delta_q_n   (idim) )

    allocate( Tri_surf%x_e_t_n       (idim) )
    allocate( Tri_surf%x_e_q_n       (idim) )
    allocate( Tri_surf%x_f_t_delt_n  (idim) )
    allocate( Tri_surf%x_f_q_delt_n  (idim) )

end subroutine alloc_surf_diff_type_1d

!#######################################################################

subroutine alloc_surf_diff_type_2d ( Tri_surf, idim, jdim )

type(surf_diff_type), intent(inout) :: Tri_surf
integer,              intent(in)    :: idim, jdim

    Tri_surf%array_order = 2

    allocate( Tri_surf%mu_delt_n   (idim, jdim) )
    allocate( Tri_surf%nu_n        (idim, jdim) )
    allocate( Tri_surf%e_n1        (idim, jdim) )
    allocate( Tri_surf%f_t_delt_n1 (idim, jdim) )
    allocate( Tri_surf%f_q_delt_n1 (idim, jdim) )
    allocate( Tri_surf%delta_t_n   (idim, jdim) )
    allocate( Tri_surf%delta_q_n   (idim, jdim) )

end subroutine alloc_surf_diff_type_2d

!#######################################################################

subroutine dealloc_surf_diff_type ( Tri_surf )

type(surf_diff_type), intent(inout) :: Tri_surf

  if ( Tri_surf%array_order == 1 ) then

      deallocate( Tri_surf%x_mu_delt_n   )
      deallocate( Tri_surf%x_nu_n        )
      deallocate( Tri_surf%x_e_n1        )
      deallocate( Tri_surf%x_f_t_delt_n1 )
      deallocate( Tri_surf%x_f_q_delt_n1 )
      deallocate( Tri_surf%x_delta_t_n   )
      deallocate( Tri_surf%x_delta_q_n   )

      deallocate( Tri_surf%x_e_t_n       )
      deallocate( Tri_surf%x_e_q_n       )
      deallocate( Tri_surf%x_f_t_delt_n  )
      deallocate( Tri_surf%x_f_q_delt_n  )

  else if ( Tri_surf%array_order == 2 ) then

      deallocate( Tri_surf%mu_delt_n   )
      deallocate( Tri_surf%nu_n        )
      deallocate( Tri_surf%e_n1        )
      deallocate( Tri_surf%f_t_delt_n1 )
      deallocate( Tri_surf%f_q_delt_n1 )
      deallocate( Tri_surf%delta_t_n   )
      deallocate( Tri_surf%delta_q_n   )

  endif

end subroutine dealloc_surf_diff_type

!#######################################################################

subroutine gcm_vert_diff_end

  if (.not.do_init) then

    if (allocated(   e_global ))    deallocate (   e_global)
    if (allocated( f_t_global ))    deallocate ( f_t_global)
    if (allocated( f_q_global ))    deallocate ( f_q_global)

  endif

end subroutine gcm_vert_diff_end

!#######################################################################

subroutine gcm_vert_diff_down (is, js, delt,                &
                          u, v, t, q, tr,                   &
                          diff_m, diff_t, p_half, z_full,   &
                          tau_u, tau_v, dtau_datmos,        &
                          flux_tr,                          &
                          dt_u, dt_v, dt_t, dt_q, dt_tr,    &
                          Tri_surf,                         &
                          kbot                              )

integer, intent(in)                        :: is, js
real,    intent(in)                        :: delt
real,    intent(in)   , dimension(:,:,:)   :: u, v, t, q,     &
                                              diff_m, diff_t, &
                                              p_half, z_full
real,    intent(in)   , dimension(:,:,:,:) :: tr
real,    intent(in)   , dimension(:,:)     :: dtau_datmos
real,    intent(inout), dimension(:,:)     :: tau_u, tau_v
real,    intent(inout), dimension(:,:,:)   :: flux_tr
real,    intent(inout), dimension(:,:,:)   :: dt_u, dt_v
real,    intent(in),    dimension(:,:,:)   :: dt_t, dt_q
real,    intent(inout), dimension(:,:,:,:) :: dt_tr
type(surf_diff_type), intent(inout)        :: Tri_surf

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(u,1),size(u,2),size(u,3)) :: tt, mu, nu
real    :: gcp
integer :: i, j, kb, ie, je

!-----------------------------------------------------------------------

  if(do_init) call error_mesg ('gcm_vert_diff_down in vert_diff_mod',  &
      'the initialization routine gcm_vert_diff_init has not been called', &
       FATAL)
    
 ie = is + size(t,1) -1
 je = js + size(t,2) -1

 order(is:ie,js:je) = order(is:ie,js:je) +1
 if (count(order(is:ie,js:je) /= 1) > 0) call error_mesg (  &
          'gcm_vert_diff_down in vert_diff_mod', &
          'subroutine called out of order', FATAL)

 gcp = grav/cp
 tt  = t + z_full*gcp   ! the vertical gradient of tt determines the
                        ! diffusive flux of temperature

 call compute_mu (p_half, mu)
 call compute_nu (diff_m, p_half, z_full, t, q, nu) 

!  diffuse u-momentum and v_momentum
 call uv_vert_diff (delt, mu, nu, u, v, dtau_datmos, tau_u, tau_v,  &
                    dt_u, dt_v, kbot)

!  recompute nu for a different diffusivity
 call compute_nu   (diff_t, p_half, z_full, t, q, nu)

!  diffuse tracers 
 call tr_vert_diff (delt, mu, nu, tr, flux_tr, dt_tr, kbot )

!  downward sweep of tridiagonal solver for temperature and specific humidity
 call vert_diff_down_2                            & 
         (delt, mu, nu, tt, q, dt_t, dt_q,        &  
         e_global             (is:ie,js:je,:),    &
         f_t_global           (is:ie,js:je,:),    &
         f_q_global           (is:ie,js:je,:),    &
         Tri_surf%mu_delt_n   (is:ie,js:je),      &
         Tri_surf%nu_n        (is:ie,js:je),      &
         Tri_surf%e_n1        (is:ie,js:je),      &
         Tri_surf%f_t_delt_n1 (is:ie,js:je),      &
         Tri_surf%f_q_delt_n1 (is:ie,js:je),      &
         Tri_surf%delta_t_n   (is:ie,js:je),      &
         Tri_surf%delta_q_n   (is:ie,js:je), kbot)

!-----------------------------------------------------------------------

end subroutine gcm_vert_diff_down

!#######################################################################

subroutine gcm_vert_diff_surf_down (avail, dsens_datmos, dsens_dsurf,  &
                                           devap_datmos, devap_dsurf,  &
                                           drad_dsurf,                 &
                                           sens, evap, Tri_surf        )

!-----------------------------------------------------------------------

logical, intent(in), dimension(:) :: avail
real   , intent(in), dimension(:) :: dsens_datmos, devap_datmos, drad_dsurf
real, intent(inout), dimension(:) :: sens, dsens_dsurf, evap, devap_dsurf
type(surf_diff_type), intent(inout) :: Tri_surf

real, dimension(size(sens)) :: zero

!-----------------------------------------------------------------------

! assume this is called for all points
 order = order +1
 if (count(order /= 2) > 0) call error_mesg (  &
          'gcm_vert_diff_surf_down in vert_diff_mod', &
          'subroutine called out of order', FATAL     )

! temperature
  call diff_surface_down_1d (avail, Tri_surf%x_mu_delt_n, Tri_surf%x_nu_n,   &
                                    Tri_surf%x_e_n1, Tri_surf%x_f_t_delt_n1, &
                                    Tri_surf%x_delta_t_n,                    &
                             dsens_datmos, dsens_dsurf, drad_dsurf, sens,    &
                             cp, Tri_surf%x_e_t_n, Tri_surf%x_f_t_delt_n     )

! moisture
  zero = 0.0
  call diff_surface_down_1d (avail, Tri_surf%x_mu_delt_n, Tri_surf%x_nu_n,   &
                                    Tri_surf%x_e_n1, Tri_surf%x_f_q_delt_n1, &
                                    Tri_surf%x_delta_q_n,                    &
                             devap_datmos, devap_dsurf, zero, evap,          &
                             1.0, Tri_surf%x_e_q_n, Tri_surf%x_f_q_delt_n    )

!-----------------------------------------------------------------------

end subroutine gcm_vert_diff_surf_down

!#######################################################################

subroutine gcm_vert_diff_surf_up (avail, delta_t_surf, &
                                  dsens_dsurf, devap_dsurf, drad_dsurf, &
                                  sens, evap, rad, Tri_surf)

!-----------------------------------------------------------------------

logical, intent(in), dimension(:) :: avail
real   , intent(in), dimension(:) :: delta_t_surf, dsens_dsurf,  &
                                     drad_dsurf, devap_dsurf
real, intent(inout), dimension(:) :: sens, evap, rad
type(surf_diff_type), intent(inout) :: Tri_surf

real, dimension(size(sens)) :: zero
!-----------------------------------------------------------------------

! assume this is called for all points
 order = order +1
 if (count(order /= 3) > 0) call error_mesg (  &
          'gcm_vert_diff_surf_up in vert_diff_mod', &
          'subroutine called out of order', FATAL   )

 zero = 0.0

! temperature
 call diff_surface_up_1d (avail, Tri_surf%x_e_t_n, Tri_surf%x_f_t_delt_n,   &
                          delta_t_surf, dsens_dsurf, drad_dsurf, sens, rad, &
                          Tri_surf%x_delta_t_n                              )
! water vapor
 call diff_surface_up_1d (avail, Tri_surf%x_e_q_n, Tri_surf%x_f_q_delt_n,   &
                          delta_t_surf, devap_dsurf, zero, evap, zero,      &
                          Tri_surf%x_delta_q_n                              )

!-----------------------------------------------------------------------

end subroutine gcm_vert_diff_surf_up

!#######################################################################

subroutine gcm_vert_diff_up (is, js, delt, Tri_surf, dt_t, dt_q, kbot)

integer, intent(in)                      :: is, js
real,    intent(in)                      :: delt
type(surf_diff_type), intent(in)         :: Tri_surf
real,    intent(out),   dimension(:,:,:) :: dt_t, dt_q
integer, intent(in),    dimension(:,:), optional :: kbot

integer :: ie, je

 ie = is + size(dt_t,1) -1
 je = js + size(dt_t,2) -1

 order(is:ie,js:je) = order(is:ie,js:je) +1
 if (count(order(is:ie,js:je) /= 4) > 0) call error_mesg (  &
          'gcm_vert_diff_up in vert_diff_mod',   &
          'subroutine called out of order', FATAL)

 call vert_diff_up (delt ,                               &
                    e_global           (is:ie,js:je,:) , &
                    f_t_global         (is:ie,js:je,:) , &
                    Tri_surf%delta_t_n (is:ie,js:je) ,   &
                    dt_t, kbot )

 call vert_diff_up (delt ,                               &
                    e_global           (is:ie,js:je,:) , &
                    f_q_global         (is:ie,js:je,:) , &
                    Tri_surf%delta_q_n (is:ie,js:je) ,   &
                    dt_q, kbot )

 order(is:ie,js:je) = 0

end subroutine gcm_vert_diff_up

!#######################################################################

subroutine gcm_vert_diff (delt, u, v, t, q, tr,                    &
                          diff_m, diff_t, p_half, z_full,          &
                          dtau_datmos, dsens_datmos, devap_datmos, &
                          sens, evap, tau_u, tau_v, flux_tr,       &
                          dt_u, dt_v, dt_t, dt_q, dt_tr, kbot      )

!  one-step diffusion call for gcm in which there is no implicit dependence of 
!    surface fluxes on surface temperature

real,    intent(in)                          :: delt
real,    intent(in)   , dimension(:,:,:)     :: u, v, t, q, p_half, z_full, &
                                                diff_m, diff_t
real,    intent(in)   , dimension(:,:,:,:)   :: tr
real,    intent(in)   , dimension(:,:)       :: dtau_datmos, dsens_datmos, &
                                                devap_datmos
real,    intent(inout), dimension(:,:,:)     :: flux_tr
real,    intent(inout), dimension(:,:)       :: tau_u, tau_v, sens, evap
real,    intent(inout), dimension(:,:,:)     :: dt_u, dt_v, dt_t, dt_q
real,    intent(inout), dimension(:,:,:,:)   :: dt_tr

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(u,1),size(u,2),size(u,3)) :: mu, nu


!-----------------------------------------------------------------------

 call compute_mu (p_half, mu)

 call compute_nu (diff_m, p_half, z_full, t, q, nu) 

 call uv_vert_diff (delt, mu, nu, u, v, dtau_datmos, tau_u, tau_v, &
                    dt_u, dt_v, kbot )

 call compute_nu   (diff_t, p_half, z_full, t, q, nu)

 call tq_vert_diff (delt, mu, nu, t, q, z_full,  &
                    dsens_datmos, devap_datmos,  &
                    sens, evap, dt_t, dt_q, kbot )

 call tr_vert_diff (delt, mu, nu, tr, flux_tr, dt_tr, kbot )

end subroutine gcm_vert_diff

!#######################################################################

subroutine vert_diff (delt, xi, t, q, diff, p_half, z_full, &
                      flux, dflux_datmos, factor, dt_xi, kbot)

! one-step diffusion of a single field 

real,    intent(in)                          :: delt
real,    intent(in)   , dimension(:,:,:)     :: xi, t, q, diff, p_half, z_full
real,    intent(inout), dimension(:,:)       :: flux
real,    intent(in)   , dimension(:,:)       :: dflux_datmos
real,    intent(in)                          :: factor
real,    intent(inout), dimension(:,:,:)     :: dt_xi

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(xi,1),size(xi,2),size(xi,3)  ) :: mu, nu
real, dimension(size(xi,1),size(xi,2),size(xi,3)-1) :: e, f

real, dimension(size(xi,1),size(xi,2))  :: mu_delt_n, nu_n, e_n1,  &
                                           f_delt_n1, delta_xi_n


!-----------------------------------------------------------------------

 call compute_mu    (p_half, mu)

 call compute_nu    (diff, p_half, z_full, t, q, nu) 

 call vert_diff_down &
     (delt, mu, nu, xi, dt_xi, e, f, mu_delt_n, nu_n, e_n1,  &
      f_delt_n1, delta_xi_n, kbot)

 call diff_surface (mu_delt_n, nu_n, e_n1, f_delt_n1,     &
                    dflux_datmos, flux, factor, delta_xi_n)

 call vert_diff_up (delt, e, f, delta_xi_n, dt_xi, kbot)

end subroutine vert_diff


!#######################################################################

subroutine uv_vert_diff (delt, mu, nu, u, v,  &
                         dtau_datmos, tau_u, tau_v, dt_u, dt_v, kbot )

real,    intent(in)                        :: delt
real,    intent(in)   , dimension(:,:,:)   :: u, v, mu, nu
real,    intent(in)   , dimension(:,:)     :: dtau_datmos
real,    intent(inout), dimension(:,:)     :: tau_u, tau_v
real,    intent(inout), dimension(:,:,:)   :: dt_u, dt_v

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(u,1),size(u,2)) :: mu_delt_n, nu_n, e_n1,    &
                                        f_u_delt_n1, f_v_delt_n1, &
                                        delta_u_n, delta_v_n

real, dimension(size(u,1),size(u,2),size(u,3)-1) :: e, f_u, f_v
integer :: i, j, kb

!-----------------------------------------------------------------------
  
 call vert_diff_down_2 &
     (delt, mu, nu, u, v, dt_u, dt_v, e, f_u, f_v, &
      mu_delt_n, nu_n, e_n1, f_u_delt_n1, f_v_delt_n1,  &
      delta_u_n, delta_v_n, kbot)

 call diff_surface (mu_delt_n, nu_n, e_n1, f_u_delt_n1, &
                    dtau_datmos, tau_u, 1.0, delta_u_n)
 call diff_surface (mu_delt_n, nu_n, e_n1, f_v_delt_n1, &
                    dtau_datmos, tau_v, 1.0, delta_v_n)

 call vert_diff_up (delt, e, f_u, delta_u_n, dt_u, kbot)
 call vert_diff_up (delt, e, f_v, delta_v_n, dt_v, kbot)


!-----------------------------------------------------------------------

end subroutine uv_vert_diff

!#######################################################################

subroutine tq_vert_diff (delt, mu, nu, t, q,  z_full, &
                         dsens_datmos, devap_datmos, sens, evap, &
                         dt_t, dt_q, kbot)
                         

real,    intent(in)                        :: delt
real,    intent(in)   , dimension(:,:,:)   :: t, q, z_full, mu, nu
real,    intent(in)   , dimension(:,:)     :: dsens_datmos, devap_datmos
real,    intent(inout), dimension(:,:)     :: sens, evap
real,    intent(inout), dimension(:,:,:)   :: dt_t, dt_q

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(t,1),size(t,2)) :: mu_delt_n, nu_n,          &
                                        e_n1, f_t_delt_n1, f_q_delt_n1, &
                                        delta_t_n, delta_q_n,      &
                                        flux1, dflux1_dsurf

real, dimension(size(t,1),size(t,2),size(t,3)-1) :: e, f_t, f_q
real, dimension(size(t,1),size(t,2),size(t,3)  ) :: tt

integer :: i, j, kb
real    :: gcp
!-----------------------------------------------------------------------

 gcp = grav/cp
 tt  = t + z_full*gcp
  
 call vert_diff_down_2 &
     (delt, mu, nu, tt, q, dt_t, dt_q, e, f_t, f_q,    &
      mu_delt_n, nu_n, e_n1, f_t_delt_n1, f_q_delt_n1, &
      delta_t_n, delta_q_n, kbot)


 call diff_surface (mu_delt_n, nu_n, e_n1, f_t_delt_n1,  &
                    dsens_datmos, sens, cp, delta_t_n)

 call diff_surface (mu_delt_n, nu_n, e_n1, f_q_delt_n1,  &
                    devap_datmos, evap, 1.0, delta_q_n)

 call vert_diff_up (delt, e, f_t, delta_t_n, dt_t, kbot)
 call vert_diff_up (delt, e, f_q, delta_q_n, dt_q, kbot)


!-----------------------------------------------------------------------

end subroutine tq_vert_diff

!#######################################################################

subroutine tr_vert_diff (delt, mu, nu, tr, flux, dt_tr, kbot )

real,    intent(in)                        :: delt
real,    intent(in)   , dimension(:,:,:)   :: mu, nu
real,    intent(in)   , dimension(:,:,:,:) :: tr
real,    intent(inout), dimension(:,:,:)   :: flux
real,    intent(inout), dimension(:,:,:,:) :: dt_tr

integer, intent(in)   , dimension(:,:), optional :: kbot

real, dimension(size(tr,1),size(tr,2)) :: mu_delt_n, nu_n, e_n1

real, dimension(size(tr,1),size(tr,2),size(tr,4)) :: f_delt_n1, delta_tr_n
real, dimension(size(tr,1),size(tr,2)) :: dflux_dtr

real, dimension(size(tr,1),size(tr,2),size(tr,3)-1,size(tr,4)) :: f

real, dimension(size(tr,1),size(tr,2),size(tr,3)-1) :: e
integer :: i, j, n, kb, ntr
!-----------------------------------------------------------------------

 ntr  = size(tr,4)

 dflux_dtr = 0.0
  
 call vert_diff_down_n &
     (delt, mu, nu, tr, dt_tr, e, f, mu_delt_n, nu_n,  &
      e_n1, f_delt_n1, delta_tr_n, kbot)

 do n = 1, ntr
   call diff_surface (mu_delt_n, nu_n, e_n1, f_delt_n1(:,:,n),  &
                      dflux_dtr, flux(:,:,n), 1.0, delta_tr_n(:,:,n))

   call vert_diff_up (delt, e, f(:,:,:,n), delta_tr_n(:,:,n),  &
                      dt_tr(:,:,:,n), kbot)
 end do

!-----------------------------------------------------------------------

end subroutine tr_vert_diff

!#######################################################################

subroutine vert_diff_down &
      (delt, mu, nu, tr, dt_tr, e, f, mu_delt_n, nu_n,  &
       e_n1, f_delt_n1, delta_tr_n, kbot)

!-----------------------------------------------------------------------

real,    intent(in)                         :: delt
real,    intent(in)    , dimension(:,:,:)   :: mu, nu
real,    intent(in)    , dimension(:,:,:)   :: tr
real,    intent(inout) , dimension(:,:,:)   :: dt_tr
real,    intent(out)   , dimension(:,:,:)   :: e
real,    intent(out)   , dimension(:,:,:)   :: f
real,    intent(out)   , dimension(:,:)     :: mu_delt_n, nu_n, e_n1
real,    intent(out)   , dimension(:,:)     :: f_delt_n1, delta_tr_n

integer, intent(in),    dimension(:,:), optional :: kbot

real, dimension(size(tr,1),size(tr,2),size(tr,3)) :: a, b, c, g

integer :: i, j, k, kb, nlev

!-----------------------------------------------------------------------

 call explicit_tend (mu, nu, tr, dt_tr)

 call compute_e  (delt, mu, nu, e, a, b, c, g)

 call compute_f (dt_tr, b, c, g, f)


 if (present(kbot)) then
    do j=1,size(tr,2)
    do i=1,size(tr,1)
        kb = kbot(i,j)
        mu_delt_n(i,j) =  mu(i,j,kb  )*delt
             nu_n(i,j) =  nu(i,j,kb  )
             e_n1(i,j) =   e(i,j,kb-1)
    enddo
    enddo
    do j=1,size(tr,2)
    do i=1,size(tr,1)
        kb = kbot(i,j)
         f_delt_n1(i,j) =     f(i,j,kb-1)*delt
        delta_tr_n(i,j) = dt_tr(i,j,kb  )*delt
    enddo
    enddo
 else
        nlev = size(mu,3)
        mu_delt_n(:,:) =       mu(:,:,nlev  )*delt
             nu_n(:,:) =       nu(:,:,nlev  )
             e_n1(:,:) =        e(:,:,nlev-1)
        f_delt_n1(:,:) =        f(:,:,nlev-1)*delt
       delta_tr_n(:,:) =    dt_tr(:,:,nlev  )*delt
 endif



!-----------------------------------------------------------------------

end subroutine vert_diff_down

!#######################################################################

subroutine vert_diff_down_2 &
      (delt, mu, nu, xi_1, xi_2, dt_xi_1, dt_xi_2, e, f_1, f_2, &
       mu_delt_n, nu_n, e_n1, f_1_delt_n1, f_2_delt_n1,         &
       delta_1_n, delta_2_n, kbot)

!-----------------------------------------------------------------------

real,    intent(in)                       :: delt
real,    intent(in)    , dimension(:,:,:) :: mu, nu, xi_1, xi_2
real,    intent(in)    , dimension(:,:,:) :: dt_xi_1, dt_xi_2
real,    intent(out)   , dimension(:,:,:) :: e, f_1, f_2
real,    intent(out)   , dimension(:,:)   :: mu_delt_n, nu_n, e_n1,    &
                                             f_1_delt_n1, f_2_delt_n1, &
                                             delta_1_n, delta_2_n

integer, intent(in),    dimension(:,:), optional :: kbot

real, dimension(size(xi_1,1),size(xi_1,2),size(xi_1,3)) :: a, b, c, g, &
                                                      dt_xi_11, dt_xi_22

integer :: i, j, k, kb, nlev

!-----------------------------------------------------------------------

! local copy of input 
  dt_xi_11 = dt_xi_1
  dt_xi_22 = dt_xi_2

 call explicit_tend (mu, nu, xi_1, dt_xi_11)
 call explicit_tend (mu, nu, xi_2, dt_xi_22)

 call compute_e (delt, mu, nu, e, a, b, c, g)

 call compute_f (dt_xi_11, b, c, g, f_1)
 call compute_f (dt_xi_22, b, c, g, f_2)

 if (present(kbot)) then
    do j=1,size(xi_1,2)
    do i=1,size(xi_1,1)
        kb = kbot(i,j)
        mu_delt_n(i,j)  =      mu(i,j,kb  )*delt
             nu_n(i,j)  =      nu(i,j,kb  )
            e_n1(i,j)  =       e(i,j,kb-1)
     f_1_delt_n1(i,j)  =     f_1(i,j,kb-1)*delt
     f_2_delt_n1(i,j)  =     f_2(i,j,kb-1)*delt
        delta_1_n(i,j)  = dt_xi_11(i,j,kb  )*delt
        delta_2_n(i,j)  = dt_xi_22(i,j,kb  )*delt
    enddo
    enddo
 else
        nlev = size(mu,3)
        mu_delt_n(:,:)  =      mu(:,:,nlev  )*delt
             nu_n(:,:)  =      nu(:,:,nlev  )
            e_n1(:,:)  =       e(:,:,nlev-1)
     f_1_delt_n1(:,:)  =     f_1(:,:,nlev-1)*delt
     f_2_delt_n1(:,:)  =     f_2(:,:,nlev-1)*delt
        delta_1_n(:,:)  = dt_xi_11(:,:,nlev  )*delt
        delta_2_n(:,:)  = dt_xi_22(:,:,nlev  )*delt
 endif



!-----------------------------------------------------------------------

end subroutine vert_diff_down_2

!#######################################################################

subroutine vert_diff_down_n &
      (delt, mu, nu, tr, dt_tr, e, f, mu_delt_n, nu_n,  &
       e_n1, f_delt_n1, delta_tr_n, kbot)

!-----------------------------------------------------------------------

real,    intent(in)                         :: delt
real,    intent(in)    , dimension(:,:,:)   :: mu, nu
real,    intent(in)    , dimension(:,:,:,:) :: tr
real,    intent(inout) , dimension(:,:,:,:) :: dt_tr
real,    intent(out)   , dimension(:,:,:)   :: e
real,    intent(out)   , dimension(:,:,:,:) :: f
real,    intent(out)   , dimension(:,:)     :: mu_delt_n, nu_n, e_n1
real,    intent(out)   , dimension(:,:,:)   :: f_delt_n1, delta_tr_n

integer, intent(in),    dimension(:,:), optional :: kbot

real, dimension(size(tr,1),size(tr,2),size(tr,3)) :: a, b, c, g

integer :: i, j, k, n, kb, nlev, ntr

!-----------------------------------------------------------------------

  ntr = size(tr,4)

 call compute_e  (delt, mu, nu, e, a, b, c, g)

 do n = 1, ntr
   call explicit_tend (mu, nu, tr(:,:,:,n), dt_tr(:,:,:,n))
   call compute_f (dt_tr(:,:,:,n), b, c, g, f(:,:,:,n))
 end do


 if (present(kbot)) then
    do j=1,size(tr,2)
    do i=1,size(tr,1)
        kb = kbot(i,j)
        mu_delt_n(i,j) =  mu(i,j,kb  )*delt
             nu_n(i,j) =  nu(i,j,kb  )
            e_n1(i,j) =   e(i,j,kb-1)
    enddo
    enddo
    do n=1,size(tr,4)
    do j=1,size(tr,2)
    do i=1,size(tr,1)
        kb = kbot(i,j)
       f_delt_n1(i,j,n) =     f(i,j,kb-1,n)*delt
       delta_tr_n(i,j,n) = dt_tr(i,j,kb  ,n)*delt
    enddo
    enddo
    enddo
 else
    nlev = size(mu,3)
    mu_delt_n(:,:)   =    mu(:,:,nlev    )*delt
         nu_n(:,:)   =    nu(:,:,nlev    )
        e_n1(:,:)   =     e(:,:,nlev-1  )
   f_delt_n1(:,:,:) =     f(:,:,nlev-1,:)*delt
   delta_tr_n(:,:,:) = dt_tr(:,:,nlev  ,:)*delt
 endif



!-----------------------------------------------------------------------

end subroutine vert_diff_down_n

!#######################################################################

subroutine diff_surface (mu_delt, nu, e_n1, f_n1_delt,  &
                         dflux_datmos, flux, factor, delta_xi)

!-----------------------------------------------------------------------

real, intent(in)   , dimension(:,:) :: mu_delt, nu, e_n1, f_n1_delt,  &
                                       dflux_datmos
real, intent(inout), dimension(:,:) :: flux, delta_xi
real, intent(in) :: factor

!-----------------------------------------------------------------------

 real, dimension(size(flux,1),size(flux,2)) :: b, c, g, e_n, f_n_delt, zero

 zero = 0.0

 call diff_surface_down_2d (mu_delt, nu, e_n1, f_n1_delt, delta_xi, &
                 dflux_datmos, zero, zero, flux, factor, e_n, f_n_delt)

 call diff_surface_up_2d (e_n, f_n_delt, zero, zero, zero,  &
                          flux, zero, delta_xi)

!-----------------------------------------------------------------------

end subroutine diff_surface

!#######################################################################

subroutine diff_surface_down_2d (mu_delt, nu, e_n1, f_n1_delt, delta_xi,  &
                                 dflux_datmos, dflux_dsurf, dflux1_dsurf, &
                                 flux, factor, e_n, f_n_delt)

!-----------------------------------------------------------------------

real, intent(in)   , dimension(:,:) :: mu_delt, nu, e_n1, f_n1_delt,  &
                                       dflux_datmos, delta_xi, dflux1_dsurf
real, intent(inout), dimension(:,:) :: flux, dflux_dsurf
real, intent(out)  , dimension(:,:) :: e_n, f_n_delt
real, intent(in) :: factor

!-----------------------------------------------------------------------

 real, dimension(size(flux,1),size(flux,2)) :: a, b, c, g
 real :: fff

 fff = 1.0/factor

 c     = - mu_delt * nu
 b     =   1.0 - c - mu_delt * dflux_datmos*fff
 a     = - mu_delt *dflux_dsurf*fff
 g     =   1./ (b + c * e_n1)

 e_n      = - a * g
 f_n_delt = (delta_xi + mu_delt * flux * fff- c * f_n1_delt) * g    

 flux        =  flux        + dflux_datmos * f_n_delt 
 dflux_dsurf =  dflux_dsurf + dflux_datmos * e_n  

!-----------------------------------------------------------------------

end subroutine diff_surface_down_2d

!#######################################################################

subroutine diff_surface_down_1d (avail, mu_delt, nu, e_n1, f_n1_delt, &
                                 delta_xi, dflux_datmos, dflux_dsurf, &
                                 dflux1_dsurf, flux, factor, e_n, f_n_delt)

!-----------------------------------------------------------------------

logical, intent(in)   , dimension(:) :: avail
real, intent(in)   , dimension(:) :: mu_delt, nu, e_n1, f_n1_delt,  &
                                     dflux_datmos, delta_xi, dflux1_dsurf
real, intent(inout), dimension(:) :: flux, dflux_dsurf
real, intent(out)  , dimension(:) :: e_n, f_n_delt
real, intent(in) :: factor

!-----------------------------------------------------------------------

 real, dimension(size(flux)) :: a, b, c, g
 real :: fff

 fff = 1.0/factor

 where(avail)

   c     = - mu_delt * nu
   b     =   1.0 - c - mu_delt * dflux_datmos*fff
   a     = - mu_delt *dflux_dsurf*fff
   g     =   1./ (b + c * e_n1)

   e_n      = - a * g
   f_n_delt = (delta_xi + mu_delt * flux * fff- c * f_n1_delt) * g    

   flux        =  flux        + dflux_datmos * f_n_delt 
   dflux_dsurf =  dflux_dsurf + dflux_datmos * e_n   

 endwhere

!-----------------------------------------------------------------------

end subroutine diff_surface_down_1d

!#######################################################################

subroutine diff_surface_up_2d (e_n, f_n_delt, delta_surf, dflux_dsurf, &
                               dflux1_dsurf, flux, flux1, delta_xi)

!-----------------------------------------------------------------------

real, intent(in)   , dimension(:,:) :: e_n, f_n_delt, delta_surf,  &
                                       dflux_dsurf, dflux1_dsurf
real, intent(inout), dimension(:,:) :: flux, flux1
real, intent(out)  , dimension(:,:) :: delta_xi

!-----------------------------------------------------------------------

 flux     = flux       + delta_surf * dflux_dsurf
 flux1    = flux1      - delta_surf * dflux1_dsurf
 delta_xi = f_n_delt   + delta_surf * e_n

!-----------------------------------------------------------------------

end subroutine diff_surface_up_2d

!#######################################################################

subroutine diff_surface_up_1d (avail, e_n, f_n_delt, delta_surf, dflux_dsurf, &
                               dflux1_dsurf, flux, flux1, delta_xi)

!-----------------------------------------------------------------------

logical, intent(in), dimension(:) :: avail
real, intent(in)   , dimension(:) :: e_n, f_n_delt, delta_surf,  &
                                     dflux_dsurf, dflux1_dsurf
real, intent(inout), dimension(:) :: flux, flux1
real, intent(out)  , dimension(:) :: delta_xi

!-----------------------------------------------------------------------

 delta_xi = 0.0

 where(avail)
   flux     = flux       + delta_surf * dflux_dsurf
   flux1    = flux1      - delta_surf * dflux1_dsurf
   delta_xi = f_n_delt   + delta_surf * e_n
 endwhere

!-----------------------------------------------------------------------

end subroutine diff_surface_up_1d

!#######################################################################

subroutine vert_diff_up (delt, e, f, delta_xi_n, dt_xi, kbot)

!-----------------------------------------------------------------------

real,    intent(in)                      :: delt
real,    intent(in),    dimension(:,:,:) :: e, f
real,    intent(in) ,   dimension(:,:)   :: delta_xi_n
real,    intent(out),   dimension(:,:,:) :: dt_xi
integer, intent(in),    dimension(:,:), optional :: kbot

integer :: i, j, k, kb, nlev
!-----------------------------------------------------------------------

 if (present(kbot)) then
     do j = 1, size(dt_xi,2)
     do i = 1, size(dt_xi,1)
         kb = kbot(i,j)
         dt_xi(i,j,kb) = delta_xi_n(i,j)/delt
         do k = kb -1, 1, -1
           dt_xi(i,j,k) = e(i,j,k)*dt_xi(i,j,k+1) + f(i,j,k)
         end do
     end do
     end do
 else
    nlev = size(dt_xi,3)
    dt_xi(:,:,nlev) = delta_xi_n/delt
    do k = size(dt_xi,3)-1, 1, -1
      dt_xi(:,:,k) = e(:,:,k)*dt_xi(:,:,k+1) + f(:,:,k)
    end do
 endif

!-----------------------------------------------------------------------

end subroutine vert_diff_up

!#######################################################################

subroutine compute_e (delt, mu, nu, e, a, b, c, g)

!-----------------------------------------------------------------------

real,    intent(in)                       :: delt
real,    intent(in)    , dimension(:,:,:) :: mu, nu
real,    intent(out)   , dimension(:,:,:) :: e, a, b, c, g

integer :: k, nlev

!-----------------------------------------------------------------------

 nlev = size(mu,3)

 a(:,:,1:nlev-1) = - mu(:,:,1:nlev-1)*nu(:,:,2:nlev)*delt
 a(:,:,nlev    ) =   0.0
 c(:,:,2:nlev  ) = - mu(:,:,2:nlev  )*nu(:,:,2:nlev)*delt
 c(:,:,1       ) =   0.0

 b = 1.0 - a - c

 e(:,:,1)   =   - a(:,:,1)/b(:,:,1)
 do  k= 2, nlev - 1
    g(:,:,k) = 1.0/(b(:,:,k) + c(:,:,k)*e(:,:,k-1))
    e(:,:,k) = - a(:,:,k)*g(:,:,k)
 enddo

!-----------------------------------------------------------------------

end subroutine compute_e

!#######################################################################

subroutine compute_f (dt_xi, b, c, g, f)

!-----------------------------------------------------------------------
real,    intent(in)    , dimension(:,:,:) :: dt_xi, b, c, g
real,    intent(out)   , dimension(:,:,:) :: f

integer :: k
!-----------------------------------------------------------------------

 f(:,:,1) =   dt_xi(:,:,1)/b(:,:,1)

 do  k = 2, size(b,3)-1
    f(:,:,k) = (dt_xi(:,:,k) - c(:,:,k)*f(:,:,k-1))*g(:,:,k)
 enddo

!-----------------------------------------------------------------------

end subroutine compute_f

!#######################################################################

subroutine explicit_tend (mu, nu, xi, dt_xi)

!-----------------------------------------------------------------------

real,    intent(in)    , dimension(:,:,:) :: mu, nu, xi
real,    intent(inout) , dimension(:,:,:) :: dt_xi

real, dimension(size(xi,1),size(xi,2),size(xi,3)) :: fluxx

integer :: k, nlev

!-----------------------------------------------------------------------

 nlev = size(mu,3)

 fluxx(:,:,1)      = 0.0
 fluxx(:,:,2:nlev) = nu(:,:,2:nlev)*(xi(:,:,2:nlev) - xi(:,:,1:nlev-1))

 dt_xi(:,:,1:nlev-1) = dt_xi(:,:,1:nlev-1) +  &
    mu(:,:,1:nlev-1)*(fluxx(:,:,2:nlev) - fluxx(:,:,1:nlev-1))
 dt_xi(:,:,nlev)     = dt_xi(:,:,nlev) - mu(:,:,nlev)*fluxx(:,:,nlev)

!-----------------------------------------------------------------------

end subroutine explicit_tend

!#######################################################################

subroutine compute_mu (p_half, mu)

!-----------------------------------------------------------------------
real,    intent(in)    , dimension(:,:,:) :: p_half
real,    intent(out)   , dimension(:,:,:) :: mu

integer :: nlev
!-----------------------------------------------------------------------

nlev = size(mu,3)

mu(:,:,1:nlev) = grav / (p_half(:,:,2:nlev+1) -p_half(:,:,1:nlev))

!-----------------------------------------------------------------------

end subroutine compute_mu


!#######################################################################

subroutine compute_nu (diff, p_half, z_full, t, q, nu)

!-----------------------------------------------------------------------
real,    intent(in)    , dimension(:,:,:) :: diff, p_half, z_full, t, q
real,    intent(out)   , dimension(:,:,:) :: nu

real, dimension(size(t,1),size(t,2),size(t,3)) :: rho_half, tt
integer :: nlev
!-----------------------------------------------------------------------

nlev = size(nu,3)

tt = t * (1.0 + d608*q)           ! virtual temperature
rho_half(:,:,2:nlev) =  &         ! density at half levels
      2.0*p_half(:,:,2:nlev)/(rdgas*(tt(:,:,2:nlev)+tt(:,:,1:nlev-1)))

nu(:,:,2:nlev) = rho_half(:,:,2:nlev)*diff(:,:,2:nlev) /  &
                  (z_full(:,:,1:nlev-1)-z_full(:,:,2:nlev))
!-----------------------------------------------------------------------

end subroutine compute_nu

!#######################################################################

subroutine switch_surf_diff_type_order ( Tri_surf_in, Tri_surf_out )

type(surf_diff_type), intent(in)    :: Tri_surf_in
type(surf_diff_type), intent(inout) :: Tri_surf_out

integer :: j, is, ie, nlon, nlat

  if ( Tri_surf_in%array_order == 2 .and. Tri_surf_out%array_order == 1 ) then

     nlon = size( Tri_surf_in%nu_n, 1 )
     nlat = size( Tri_surf_in%nu_n, 2 )

     if ( size( Tri_surf_out%x_nu_n ) /= nlat ) call error_mesg (   &
                    'switch_surf_diff_type_order in vert_diff_mod', &
                    'inconsistent dimensions', FATAL )

     do j = 1, nlat

       is = nlon*(j-1)+1
       ie = nlon*j

       Tri_surf_out%x_mu_delt_n   (is:ie) = Tri_surf_in%mu_delt_n   (:,j)
       Tri_surf_out%x_nu_n        (is:ie) = Tri_surf_in%nu_n        (:,j)
       Tri_surf_out%x_e_n1        (is:ie) = Tri_surf_in%e_n1        (:,j)
       Tri_surf_out%x_f_t_delt_n1 (is:ie) = Tri_surf_in%f_t_delt_n1 (:,j)
       Tri_surf_out%x_f_q_delt_n1 (is:ie) = Tri_surf_in%f_q_delt_n1 (:,j)
       Tri_surf_out%x_delta_t_n   (is:ie) = Tri_surf_in%delta_t_n   (:,j)
       Tri_surf_out%x_delta_q_n   (is:ie) = Tri_surf_in%delta_q_n   (:,j)

     enddo

  else if ( Tri_surf_in %array_order == 1 .and.  &
            Tri_surf_out%array_order == 2 ) then

     nlon = size( Tri_surf_out%nu_n, 1 )
     nlat = size( Tri_surf_out%nu_n, 2 )

     if ( size( Tri_surf_out%x_nu_n ) /= nlat ) call error_mesg (   &
                    'switch_surf_diff_type_order in vert_diff_mod', &
                    'inconsistent dimensions', FATAL )

     do j = 1, nlat

       is = nlon*(j-1)+1
       ie = nlon*j

       Tri_surf_out%delta_t_n (:,j) = Tri_surf_in%x_delta_t_n (is:ie)
       Tri_surf_out%delta_q_n (:,j) = Tri_surf_in%x_delta_q_n (is:ie)

     enddo

  else

     if ( Tri_surf_in%array_order + Tri_surf_out%array_order /= 3 )  &
      call error_mesg ('switch_surf_diff_type_order in vert_diff_mod', &
                       'arguments have the wrong order', FATAL )

  endif

end subroutine switch_surf_diff_type_order

!#######################################################################

end module vert_diff_mod

