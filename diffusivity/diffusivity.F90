
module diffusivity_mod

!=======================================================================
!
!                          DIFFUSIVITY MODULE
!
!     Routines for computing atmospheric diffusivities in the 
!       planetary boundary layer and in the free atmosphere
!
!=======================================================================


use constants_mod, only : grav, vonkarm, cp, rdgas, rvgas

use utilities_mod, only:  error_mesg, FATAL, file_exist,   &
                          check_nml_error, open_file,      &
                          get_my_pe, close_file

use monin_obukhov_mod, only : mo_diff

implicit none
private

! public interfaces
!=======================================================================

 public diffusivity, pbl_depth

!=======================================================================

! form of iterfaces

!=======================================================================
! subroutine diffusivity (t, q, u, v, z_full, z_half, u_star, b_star, 
!                         h, k_m, k_t)

! input:
  
!        t     : real, dimension(:,:,:) -- (:,:,pressure), third index running
!                          from top of atmosphere to bottom
!                 temperature (K)
!
!        q     : real, dimension(:,:,:)
!                 water vapor specific humidity (nondimensional)
!
!        u     : real, dimension(:,:)
!                 zonal wind (m/s)
!
!        v     : real, dimension(:,:,:) 
!                 meridional wind (m/s) 
!
!        z_full  : real, dimension(:,:,: 
!                 height of full levels (m)
!                 1 = top of atmosphere; size(p_half,3) = surface
!                 size(z_full,3) = size(t,3)
!
!        z_half  : real, dimension(:,:,:)
!                 height of  half levels (m)
!                 size(z_half,3) = size(t,3) +1
!              z_half(:,:,size(z_half,3)) must be height of surface!
!                                  (if you are not using eta-model)
!
!        u_star: real, dimension(:,:)
!                friction velocity (m/s)
!
!        b_star: real, dimension(:,:)
!                buoyancy scale (m/s**2)

!   (u_star and b_star can be obtained by calling 
!     mo_drag in monin_obukhov_mod)

! output:

!        h     : real, dimension(:,:,) 
!                 depth of planetary boundary layer (m)
!
!        k_m   : real, dimension(:,:,:)
!                diffusivity for momentum (m**2/s)
!
!                defined at half-levels
!                size(k_m,3) should be at least as large as size(t,3)
!                only the returned values at 
!                      levels 2 to size(t,3) are meaningful
!                other values will be returned as zero
!
!        k_t   : real, dimension(:,:,:)
!                diffusivity for temperature and scalars (m**2/s)
!
!
!=======================================================================


!--------------------- version number ----------------------------------

character(len=128) :: version = '$Id: diffusivity.F90,v 1.2 2000/08/04 18:44:23 fms Exp $'
character(len=128) :: tag = '$Name: calgary $'

!=======================================================================

!  DEFAULT VALUES OF NAMELIST PARAMETERS:

logical :: fixed_depth         = .false.
real    :: depth_0             =  5000.0
real    :: frac_inner          =  0.1
real    :: rich_crit_pbl       =  1.0
real    :: entr_ratio          =  0.2
real    :: parcel_buoy         =  2.0
real    :: znom                =  1000.0
logical :: free_atm_diff       = .false.
logical :: free_atm_skyhi_diff = .false.
logical :: pbl_supersource     = .false.
real    :: rich_crit_diff      =  0.25
real    :: mix_len             = 30.
real    :: rich_prandtl        =  1.00
real    :: background_m        =  0.0
real    :: background_t        =  0.0

namelist /diffusivity_nml/ fixed_depth, depth_0, frac_inner,& 
                           rich_crit_pbl, entr_ratio, parcel_buoy,&
                           znom, free_atm_diff, free_atm_skyhi_diff,&
                           pbl_supersource,&
                           rich_crit_diff, mix_len, rich_prandtl,&
                           background_m, background_t
                          
!=======================================================================

!  OTHER MODULE VARIABLES

real    :: small  = 1.e-04
real    :: gcp    = grav/cp
logical :: init   = .false.

real, parameter :: d608 = (rvgas-rdgas)/rdgas


contains

!=======================================================================

subroutine diffusivity_init

integer :: unit, ierr, io

!------------------- read namelist input -------------------------------

      if (file_exist('input.nml')) then
         unit = open_file ('input.nml', action='read')
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=diffusivity_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'diffusivity_nml')
         enddo
  10     call close_file (unit)

!------------------- dummy checks --------------------------------------
         if (frac_inner .le. 0. .or. frac_inner .ge. 1.) &
            call error_mesg ('diffusivity_init',  &
            'frac_inner must be between 0 and 1', FATAL) 
         if (rich_crit_pbl .lt. 0.) &
            call error_mesg ('diffusivity_init',  &
           'rich_crit_pbl must be greater than or equal to zero', FATAL)
         if (entr_ratio .lt. 0.) &
            call error_mesg ('diffusivity_init',  &
            'entr_ratio must be greater than or equal to zero', FATAL)
         if (znom .le. 0.) &
            call error_mesg ('diffusivity_init',  &
            'znom must be greater than zero', FATAL)
         if (.not.free_atm_diff .and. free_atm_skyhi_diff)&
            call error_mesg ('diffusivity_init',  &
            'free_atm_diff must be set to true if '//&
            'free_atm_skyhi_diff = .true.', FATAL) 
         if (rich_crit_diff .le. 0.) &
            call error_mesg ('diffusivity_init',  &
            'rich_crit_diff must be greater than zero', FATAL)
         if (mix_len .lt. 0.) &
            call error_mesg ('diffusivity_init',  &
            'mix_len must be greater than or equal to zero', FATAL)
         if (rich_prandtl .lt. 0.) &
            call error_mesg ('diffusivity_init',  &
            'rich_prandtl must be greater than or equal to zero', FATAL)
         if (background_m .lt. 0.) &
            call error_mesg ('diffusivity_init',  &
            'background_m must be greater than or equal to zero', FATAL)
         if (background_t .lt. 0.) &
            call error_mesg ('diffusivity_init',  &
            'background_t must be greater than or equal to zero', FATAL)
      
      endif  !end of reading input.nml

!---------- output namelist to log-------------------------------------

      unit = open_file ('logfile.out', action='append')
      if ( get_my_pe() == 0 ) then
           write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
           write (unit, nml=diffusivity_nml)
      endif
      call close_file (unit)

init = .true.

return
end subroutine diffusivity_init

!=======================================================================

subroutine diffusivity(t, q, u, v, z_full, z_half, u_star, b_star,&
                       h, k_m, k_t, kbot)


real,    intent(in),           dimension(:,:,:) :: t, q, u, v
real,    intent(in),           dimension(:,:,:) :: z_full, z_half
real,    intent(in),           dimension(:,:)   :: u_star, b_star
real,    intent(out),          dimension(:,:,:) :: k_m, k_t
real,    intent(out),          dimension(:,:)   :: h
integer, intent(in), optional, dimension(:,:)   :: kbot

real, dimension(size(t,1),size(t,2),size(t,3))  :: svcp,z_full_ag
real, dimension(size(t,1),size(t,2),size(t,3)+1):: z_half_ag
real, dimension(size(t,1),size(t,2))            :: z_surf
integer                                         :: i,j,k,nlev,nlat,nlon

if(.not.init) call diffusivity_init

nlev = size(t,3)

!compute height of surface
if (present(kbot)) then
   nlat = size(t,2)
   nlon = size(t,1)
   do j=1,nlat
   do i=1,nlon
          z_surf(i,j) = z_half(i,j,kbot(i,j)+1)
   enddo
   enddo
else
   z_surf(:,:) = z_half(:,:,nlev+1)
end if


!compute density profile, and heights relative to surface
do k = 1, nlev
  z_full_ag(:,:,k) = z_full(:,:,k) - z_surf(:,:)
  z_half_ag(:,:,k) = z_half(:,:,k) - z_surf(:,:)
  svcp(:,:,k)  =   t(:,:,k)*(1. + d608*q(:,:,k)) + gcp*(z_full_ag(:,:,k))
end do
z_half_ag(:,:,nlev+1) = z_half(:,:,nlev+1) - z_surf(:,:)


if(fixed_depth)  then
   h = depth_0
else 
   call pbl_depth(svcp,u,v,z_full_ag,u_star,b_star,h,kbot=kbot)
end if

if(pbl_supersource) then
   call diffusivity_pbl_supersource (u,v,z_full_ag,z_half_ag,h,k_m,k_t)
else
   call diffusivity_pbl  (svcp, u, v, z_half_ag, h, u_star, b_star,&
                       k_m, k_t, kbot=kbot)
end if

if(free_atm_diff) &
   call diffusivity_free (svcp, u, v, z_full_ag, z_half_ag, h, k_m, k_t)

!NOTE THAT THIS LINE MUST FOLLOW DIFFUSIVITY_FREE SO THAT ENTRAINMENT
!K's DO NOT GET OVERWRITTEN IN DIFFUSIVITY_FREE SUBROUTINE
if(entr_ratio .gt. 0. .and. .not. fixed_depth) &
    call diffusivity_entr(svcp,z_full_ag,h,u_star,b_star,k_m,k_t)

!set background diffusivities
if(background_m.gt.0.0) k_m = max(k_m,background_m)
if(background_t.gt.0.0) k_t = max(k_t,background_t)


return
end subroutine diffusivity

!=======================================================================

subroutine pbl_depth(t, u, v, z, u_star, b_star, h, kbot)


real,   intent(in) ,           dimension(:,:,:) :: t, u, v, z
real,   intent(in) ,           dimension(:,:)   :: u_star,b_star
real,   intent(out),           dimension(:,:)   :: h
integer,intent(in) , optional, dimension(:,:)   :: kbot

real,    dimension(size(t,1),size(t,2),size(t,3))  :: rich
real,    dimension(size(t,1),size(t,2))            :: ws,k_t_ref,&
                                                      h_inner,tbot
real                                               :: rich1, rich2,&
                                                      h1,h2,svp,t1,t2
integer, dimension(size(t,1),size(t,2))            :: ibot
integer                                            :: i,j,k,nlon,&
                                                      nlat, nlev

nlev = size(t,3)
nlat = size(t,2)
nlon = size(t,1)

!assign ibot, compute tbot (virtual temperature at lowest level)
if (present(kbot)) then
    ibot(:,:) = kbot
    do j = 1,nlat
    do i = 1,nlon
          tbot(i,j) = t(i,j,ibot(i,j))
    enddo
    enddo
else
    ibot(:,:) = nlev
    tbot(:,:) = t(:,:,nlev)
end if


!compute richardson number for use in pbl depth of neutral/stable side
do k = 1,nlev
  rich(:,:,k) =  z(:,:,k)*grav*(t(:,:,k)-tbot(:,:))/tbot(:,:)&
                /(u(:,:,k)*u(:,:,k) + v(:,:,k)*v(:,:,k) + small )
end do

!compute ws to be used in evaluating parcel buoyancy
!ws = u_star / phi(h_inner,u_star,b_star)  .  To find phi
!a call to mo_diff is made.

h_inner(:,:)=frac_inner*znom
call mo_diff(h_inner, u_star, b_star, ws, k_t_ref)
ws = max(small,ws/vonkarm/h_inner)


do j = 1, nlat
 do i = 1, nlon

        !do neutral or stable case 
        if (b_star(i,j).le.0.) then    
              
              h1     = z(i,j,ibot(i,j))
              h(i,j) = h1
              rich1  = rich(i,j,ibot(i,j))
              do k = ibot(i,j)-1, 1, -1
                       rich2 = rich(i,j,k)
                       h2    = z(i,j,k)
                       if(rich2.gt.rich_crit_pbl) then
                             h(i,j) = h2 + (h1 - h2)*(rich2 - rich_crit_pbl)&
                                                    /(rich2 - rich1        )
                             go to 10
                       endif
                       rich1 = rich2
                       h1    = h2
              enddo 

        !do unstable case
        else

              svp    = tbot(i,j)*(1.+ &
                       (parcel_buoy*u_star(i,j)*b_star(i,j)/grav/ws(i,j)) )
              h1     = z(i,j,ibot(i,j))
              h(i,j) = h1
              t1     = tbot(i,j)
              do k = ibot(i,j)-1 , 1, -1
                       h2 = z(i,j,k)
                       t2 = t(i,j,k)
                       if (t2.gt.svp) then
                             h(i,j) = h2 + (h1 - h2)*(t2 - svp)/(t2 - t1 )
                             go to 10
                       end if
                       h1 = h2
                       t1 = t2
              enddo

        end if
10 continue
  enddo
enddo

return
end subroutine pbl_depth

!=======================================================================

subroutine diffusivity_pbl(t, u, v, z_half, h, u_star, b_star, &
                           k_m, k_t, kbot)

real,    intent(in)  ,           dimension(:,:,:) :: t, u, v, z_half
real,    intent(in)  ,           dimension(:,:)   :: h, u_star, b_star
real,    intent(out) ,           dimension(:,:,:) :: k_m, k_t
integer, intent(in)  , optional, dimension(:,:)   :: kbot

real, dimension(size(t,1),size(t,2))              :: h_inner, k_m_ref,&
                                                     k_t_ref, factor
real, dimension(size(t,1),size(t,2),size(t,3)+1)  :: zm
real                                              :: h_inner_max
integer                                           :: i,j, k, kk, nlev


nlev = size(t,3)

!assign z_half to zm, and set to zero any values of zm < 0.
!the setting to zero is necessary so that when using eta model
!below ground half levels will have zero k_m and k_t
zm = z_half
if (present(kbot)) then
   where(zm < 0.)
        zm = 0.
   end where
end if

h_inner    = frac_inner*h
h_inner_max = maxval(h_inner)

kk = nlev
do k = 2, nlev
  if( minval(zm(:,:,k)) < h_inner_max) then
      kk = k
      exit
  end if
end do

k_m = 0.0
k_t = 0.0

call mo_diff(h_inner        , u_star, b_star, k_m_ref         , k_t_ref)
call mo_diff(zm(:,:,kk:nlev), u_star, b_star, k_m(:,:,kk:nlev), k_t(:,:,kk:nlev))

do k = 2, nlev
  where(zm(:,:,k) >= h_inner .and. zm(:,:,k) < h) 
    factor = (zm(:,:,k)/h_inner)* &
             (1.0 - (zm(:,:,k) - h_inner)/(h - h_inner))**2
    k_m(:,:,k) = k_m_ref*factor
    k_t(:,:,k) = k_t_ref*factor
  end where
end do

return
end subroutine diffusivity_pbl

!=======================================================================

subroutine diffusivity_pbl_supersource(u, v, z_full, z_half, h, k_m, k_t)

real, intent(in)  , dimension(:,:,:) :: u, v, z_full, z_half
real, intent(in)  , dimension(:,:)   :: h
real, intent(out) , dimension(:,:,:) :: k_m, k_t

integer                                        :: k, nlev
real, dimension(size(z_full,1),size(z_full,2)) :: elmix, htcrit
real, dimension(size(z_full,1),size(z_full,2)) :: delta_u, delta_v, delta_z

nlev = size(z_full,3)

k_m = 0.
htcrit = frac_inner*h

do k = 2, nlev

   !compute mixing length
   elmix(:,:) = 0.

   where(z_half(:,:,k) < htcrit .and. z_half(:,:,k) > 0.)
         elmix(:,:) = vonkarm*z_half(:,:,k)
   end where

   where(z_half(:,:,k) >= htcrit .and. z_half(:,:,k) < h)
         elmix(:,:) = vonkarm*htcrit*(h(:,:)-z_half(:,:,k))/(h(:,:)-htcrit)
   end where

   !compute k_m
   delta_z = z_full(:,:,k-1) - z_full(:,:,k)
   delta_u =      u(:,:,k-1) -      u(:,:,k)
   delta_v =      v(:,:,k-1) -      v(:,:,k)
   k_m(:,:,k) =   elmix(:,:)*elmix(:,:)*&
                  sqrt(delta_u*delta_u + delta_v*delta_v)/delta_z

end do

!set k_t = k_m
k_t = k_m

return
end subroutine diffusivity_pbl_supersource

!=======================================================================

subroutine diffusivity_free(t, u, v, z, zz, h, k_m, k_t)

real, intent(in)    , dimension(:,:,:) :: t, u, v, z, zz
real, intent(in)    , dimension(:,:)   :: h
real, intent(inout) , dimension(:,:,:) :: k_m, k_t

real, dimension(size(t,1),size(t,2))   :: dz, b, speed2, rich, fri
integer                                :: k

do k = 2, size(t,3)

  dz     = z(:,:,k-1) - z(:,:,k)
  b      = grav*(t(:,:,k-1)-t(:,:,k))/t(:,:,k)
  speed2 = (u(:,:,k-1) - u(:,:,k))**2 + (v(:,:,k-1) - v(:,:,k))**2 
  rich= b*dz/(speed2+small)
  rich = max(rich, 0.0)
  if (free_atm_skyhi_diff) &
        rich(:,:) = rich(:,:) / (1.  + 1.e-04*(dz(:,:)**1.5))
  
  fri(:,:) = (1.0 - rich/rich_crit_diff)**2
  
  if (free_atm_skyhi_diff) then
         where (rich < rich_crit_diff .and. zz(:,:,k) > h) 
                k_m(:,:,k) = mix_len*mix_len*sqrt(speed2)*fri(:,:) &
                           * (1.  + 1.e-04*(dz(:,:)**1.5))/dz
                k_t(:,:,k) = k_m(:,:,k)* (0.1 + 0.9*fri(:,:))
         end where
  else
         where (rich < rich_crit_diff .and. zz(:,:,k) > h) 
                k_t(:,:,k) = mix_len*mix_len*sqrt(speed2)*fri(:,:)/dz
                k_m(:,:,k) = k_t(:,:,k)*rich_prandtl
         end where
  end if

end do

return
end subroutine diffusivity_free

!=======================================================================

subroutine diffusivity_entr(t, z,  h, u_star, b_star, k_m, k_t)

real, intent(in)    , dimension(:,:,:) :: t, z
real, intent(in)    , dimension(:,:)   :: h, u_star, b_star
real, intent(inout) , dimension(:,:,:) :: k_m, k_t

integer                                :: k, nlev

nlev=size(t,3)

do k = 2,nlev
    where (b_star .gt. 0. .and. z(:,:,k-1) .gt. h .and. &
                                z(:,:,k)   .le. h) 
        k_t(:,:,k) = (z(:,:,k-1)-z(:,:,k))*entr_ratio*t(:,:,k)* &
                      u_star*b_star/grav/max(small,t(:,:,k-1)-t(:,:,k))
        k_m(:,:,k) = k_t(:,:,k)
    end where
enddo
end subroutine diffusivity_entr

!=======================================================================

end module diffusivity_mod

!=======================================================================
