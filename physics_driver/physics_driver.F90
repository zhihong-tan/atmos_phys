
module physics_driver_mod

!-----------------------------------------------------------------------

use  moist_processes_mod, only: moist_processes,moist_processes_init,  &
                                moist_processes_end

use vert_turb_driver_mod, only: vert_turb_driver,  &
                                vert_turb_driver_init,  &
                                vert_turb_driver_end

use vert_diff_driver_mod, only: vert_diff_driver_down,  &
                                vert_diff_driver_up,    &
                                vert_diff_driver_init,  &
                                vert_diff_driver_end,   &
                                surf_diff_type

use radiation_driver_mod, only: radiation_driver,       &
                                radiation_driver_init,  &
                                radiation_driver_end

use   damping_driver_mod, only: damping_driver,      &
                                damping_driver_init, &
                                damping_driver_end

use        constants_mod, only: tfreeze,hlv,hlf,hls,kappa,  &
                                rdgas,rvgas,cp,grav

use        utilities_mod, only: file_exist, error_mesg, FATAL, NOTE,  &
                                open_file, close_file, get_my_pe, &
                                check_nml_error

use     time_manager_mod, only: time_type, get_time, operator (-)

use atmos_tracer_driver_mod, only: atmos_tracer_driver_init, atmos_tracer_driver,  &
                                   atmos_tracer_driver_end


implicit none
private

!-----------------------------------------------------------------------
!------- namelist ------

logical :: do_mcm_moist_process = .false.

namelist / physics_driver_nml / do_mcm_moist_process 


!-----------------------------------------------------------------------
!  -------- public interfaces ---------

public  physics_driver_down, physics_driver_up,   &
        physics_driver_init, physics_driver_end
public  surf_diff_type, mcm_moist_process

!-----------------------------------------------------------------------
!  -------- private interfaces/data ----------

interface check_dim
     module procedure check_dim_2d, check_dim_3d, check_dim_4d
end interface

!--------------------- version number ----------------------------------
character(len=256) :: version = '$Id: physics_driver.F90,v 1.8 2002/07/16 22:33:47 fms Exp $'
character(len=256) :: tag = '$Name: havana $'
!-----------------------------------------------------------------------

      logical :: do_init = .true., do_check_args = .true.

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine physics_driver_down (is, ie, js, je,                       &
                                 Time_prev, Time, Time_next,           &
                                 lat, lon, area,                       &
                                 p_half, p_full, z_half, z_full,       &
                                 u, v, t, q, r, um, vm, tm, qm, rm, rd,&
                                 frac_land, rough_mom,                 &
                                 albedo,    t_surf_rad,                &
                                 u_star,    b_star,                    &
                                 dtau_dv,  tau_x,  tau_y,              &
                                 udt, vdt, tdt, qdt, rdt,              &
                                 flux_sw,  flux_lw,  coszen,  gust,    &
                                 Surf_diff,                            &
                                 mask, kbot                            )

!-----------------------------------------------------------------------
        integer,          intent(in)       :: is, ie, js, je
type(time_type),          intent(in)       :: Time_prev, Time, Time_next
      real, intent(in)   ,dimension(:,:)   :: lat, lon, area
      real, intent(in)   ,dimension(:,:,:) :: p_half, p_full,   &
                                              z_half, z_full,   &
                                              u , v , t , q ,   &
                                              um, vm, tm, qm
      real, intent(in)   ,dimension(:,:,:,:) :: r,rm

      real, intent(in),    dimension(:,:) :: frac_land, rough_mom, &
                                             albedo, t_surf_rad,   &
                                             u_star, b_star,       &
                                             dtau_dv
      real, intent(inout), dimension(:,:) :: tau_x,  tau_y

      real, intent(inout),dimension(:,:,:)   :: udt,vdt,tdt,qdt
      real, intent(inout),dimension(:,:,:,:) :: rd, rdt

      real, intent(out),   dimension(:,:) :: flux_sw, flux_lw,  &
                                             coszen,  gust

      type(surf_diff_type), intent(inout) :: Surf_diff

      real, intent(in),dimension(:,:,:),optional :: mask
   integer, intent(in),dimension(:,:)  ,optional :: kbot
!-----------------------------------------------------------------------
!
!  input
!  -----
!  is,ie,js,je  starting and ending global indices
!  Time_prev    previous time, for variables um,vm,tm,qm,rm (time_type)
!  Time         current time, for variables u,v,t,q,r  (time_type)
!  Time_next    next time, used for diagnostics   (time_type)
!  lat          latitude in radians
!  lon          longitude in radians
!  area         grid box area (in m2) - currently not used
!  p_half       pressure in pascals at half levels (offset from
!                   t,q,u,v,r)
!  p_full       pressure in pascals at full levels
!  z_half       height in meters at half levels
!  z_full       height in meters at full levels
!  u            zonal wind at current time step in m/s
!  v            meridional wind at current time step in m/s
!  t            temperature at current time step in deg k
!  q            specific humidity at current time step in
!                  kg vapor / kg air
!  r            multiple 3d tracer fields at current time step
!  um,vm        zonal and meridional wind at previous time step
!  tm,qm        temperature and specific humidity at previous time step
!  rm           multiple 3d tracer fields at previous time step
!
!  inout
!  -----
!  udt          zonal wind tendency in m/s2
!  vdt          meridional wind tendency in m/s2
!  tdt          temperature tendency in deg k per sec
!  qdt          specific humidity tendency in kg vapor/kg air/sec
!  rdt          multiple tracer tendencies (in unit/unit/sec ???)
!  rd           multiple 3d diagnostic tracer fields ((in unit/unit/sec ???)
!
!  optional input
!  --------------
!  mask         mask that designates which levels do not have data
!                 present (i.e., below ground); 0.=no data, 1.=data
!  kbot         lowest level which has data
!               note:  both mask and kbot must be present together.
!
!-----------------------------------------------------------------------

   real, dimension(size(u,1),size(u,2),size(u,3)) :: diff_t,  diff_m

   integer          :: sec, day
   real             :: dt


!-----------------------------------------------------------------------
      if (do_init) call error_mesg ('physics_driver',  &
                     'physics_driver_init must be called first.', FATAL)

!-------------------- argument checking --------------------------------
      if (do_check_args) call check_args  &
                      (lat, lon, area, p_half, p_full, z_half, z_full, &
                       u, v, t, q, r, um, vm, tm, qm, rm,              &
                       udt, vdt, tdt, qdt, rdt)

!     --- compute the time step (from tau-1 to tau+1) -----

      call get_time (Time_next-Time_prev, sec, day)
      dt = real(sec + day*86400)

!-----------------------------------------------------------------------
!-------------------------- radiation ----------------------------------

      call radiation_driver (is, ie, js, je, Time, Time_next, lat,  &
                             lon, p_full, p_half, t, q, t_surf_rad,&
                             frac_land, albedo, tdt, flux_sw, &
                             flux_lw,  coszen,     &
                             mask=mask, kbot=kbot)

!-----------------------------------------------------------------------
!------------------------ damping --------------------------------------

     call damping_driver (is, js, Time_next, dt,                   &
                          p_full, p_half, z_full, z_half,          &
                          um, vm, tm, qm, rm, udt, vdt, tdt, qdt, rdt, &
                          mask=mask, kbot=kbot)

!-----------------------------------------------------------------------
!-------- need to modify vert_turb_driver (remove t_surf) --------

     call vert_turb_driver (is, js, Time, Time_next, dt, frac_land, &
                            p_half, p_full, z_half, z_full,         &
                            u_star, b_star, rough_mom,              &
                            u, v, t, q, um, vm, tm, qm,             &
                            udt, vdt, tdt, qdt,                     &
                            diff_t, diff_m, gust,                   &
                            mask=mask, kbot=kbot                    )

!-----------------------------------------------------------------------
!----------------------- process tracers fields ------------------------

  call atmos_tracer_driver (is, ie, js, je, Time, lon, lat, frac_land, &
                            p_half, p_full, r, u, v, t, q, u_star,     &
                            rdt, rm, rd, kbot)

!-----------------------------------------------------------------------
!----------------------- vertical diffusion ----------------------------

  call vert_diff_driver_down (is, js, Time_next, dt, p_half, p_full, &
                              z_full, diff_m, diff_t,                &
                              um ,vm ,tm ,qm ,rm,                    &
                              dtau_dv, tau_x, tau_y,                 &
                              udt, vdt, tdt, qdt, rdt,               &
                              Surf_diff,                             &
                              mask=mask, kbot=kbot                   )

!-----------------------------------------------------------------------

 end subroutine physics_driver_down

!#######################################################################

 subroutine physics_driver_up (is, ie, js, je,                    &
                               Time_prev, Time, Time_next,        &
                               lat, lon, area,                    &
                               p_half, p_full, omega,             &
                               u, v, t, q, r, um, vm, tm, qm, rm, &
                               frac_land,                         &
                               udt, vdt, tdt, qdt, rdt,           &
                               Surf_diff,                         &
                               lprec,   fprec,                    &
                               mask, kbot                         )

!-----------------------------------------------------------------------
        integer,          intent(in)       :: is, ie, js, je
type(time_type),          intent(in)       :: Time_prev, Time, Time_next
      real, intent(in)   ,dimension(:,:)   :: lat, lon, area
      real, intent(in)   ,dimension(:,:,:) :: p_half, p_full, omega,  &
                                              u , v , t , q ,         &
                                              um, vm, tm, qm
      real, intent(in)   ,dimension(:,:,:,:) :: r,rm

      real, intent(in),    dimension(:,:) :: frac_land

      real, intent(inout),dimension(:,:,:)   :: udt,vdt,tdt,qdt
      real, intent(inout),dimension(:,:,:,:) :: rdt

      type(surf_diff_type), intent(inout) :: Surf_diff

      real, intent(out),   dimension(:,:) :: lprec,   fprec

      real, intent(in),dimension(:,:,:),optional :: mask
   integer, intent(in),dimension(:,:)  ,optional :: kbot
!-----------------------------------------------------------------------
!
!  input
!  -----
!  is,ie,js,je  starting and ending global indices
!  Time_prev    previous time, for variables um,vm,tm,qm,rm (time_type)
!  Time         current time, for variables u,v,t,q,r  (time_type)
!  Time_next    next time, used for diagnostics   (time_type)
!  lat          latitude in radians
!  lon          longitude in radians
!  area         grid box area (in m2) - currently not used
!  p_half       pressure in pascals at half levels (offset from
!                   t,q,u,v,r)
!  p_full       pressure in pascals at full levels
!  u            zonal wind at current time step in m/s
!  v            meridional wind at current time step in m/s
!  t            temperature at current time step in deg k
!  q            specific humidity at current time step in
!                  kg vapor / kg air
!  r            multiple 3d tracer fields at current time step
!  um,vm        zonal and meridional wind at previous time step
!  tm,qm        temperature and specific humidity at previous time step
!  rm           multiple 3d tracer fields at previous time step
!
!---------------------------------------------------------------------
   real    :: dt
   integer :: day, sec

!     --- compute the time step (from tau-1 to tau+1) ---

      call get_time (Time_next-Time_prev, sec, day)
      dt = real(sec+day*86400)

!----------------------- vertical diffusion ----------------------------

    call vert_diff_driver_up (is, js, Time_next, dt, p_half, Surf_diff,&
                              tdt, qdt,  mask=mask, kbot=kbot)

!-----------------------------------------------------------------------
!-------------------------- moist processes ----------------------------

    if(.not.do_mcm_moist_process) then
      call moist_processes (is, ie, js, je, Time_next, dt, frac_land, &
                            p_half, p_full, omega,                    &
                            t, q, r, u, v, tm, qm, rm, um, vm,        &
                            tdt, qdt, rdt, udt, vdt,                  &
                            lprec, fprec, mask=mask, kbot=kbot        )
    endif

!-----------------------------------------------------------------------

 end subroutine physics_driver_up

!#######################################################################

 subroutine physics_driver_init ( Time, lonb, latb, axes, pref, &
                                  trs, Surf_diff,  mask, kbot  )

!-----------------------------------------------------------------------
type(time_type),      intent(in)    :: Time
      real,           intent(in)    :: lonb(:), latb(:)
   integer,           intent(in)    :: axes(4)
      real,           intent(in)    :: pref(:,:)
      real,           intent(inout) :: trs(:,:,:,:)
type(surf_diff_type), intent(inout) :: Surf_diff
      real, optional, intent(in)    :: mask(:,:,:)
   integer, optional, intent(in)    :: kbot(:,:)
!-----------------------------------------------------------------------
!
!     Time       = current time (time_type)
!     lonb       = longitude in radians of the grid box edges
!     latb       = latitude  in radians of the grid box edges
!     axes       = axis indices, (/x,y,pf,ph/)
!                    (returned from diag axis manager)
!     pref       = two reference profiles of pressure at nlev+1 levels
!                    pref(nlev+1,1)=101325. and pref(nlev+1,2)=81060.
!
!-----------------------------------------------------------------------
      integer  unit, id, jd, kd
      integer  :: ierr, io, unit

      id = size(lonb)-1; jd = size(latb)-1; kd = size(trs,3)

!--------------------------------------------------------------------
!----- read namelist -----
 
   if ( file_exist('input.nml')) then
    unit = open_file (file='input.nml', action='read')
    ierr=1; do while (ierr /= 0)
    read  (unit, nml=physics_driver_nml, iostat=io, end=10)
    ierr = check_nml_error(io,'physics_driver_nml')
    enddo
10   call close_file (unit)
   endif

!!      ----- write namelist -----

    unit = open_file ('logfile.out', action='append')
    if (get_my_pe() == 0) then
      write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
      write (unit, nml=physics_driver_nml)
    endif
    call close_file (unit)

!-----------------------------------------------------------------------
!------------- initialization for various schemes ----------------------

      if(.not.do_mcm_moist_process) then
        call  moist_processes_init ( id, jd, kd, axes, Time )
      endif

      call vert_turb_driver_init ( id, jd, kd, axes, Time )

      call vert_diff_driver_init ( Surf_diff, id, jd, kd, axes, Time )

      call   damping_driver_init ( lonb, latb, axes, Time )

      call radiation_driver_init ( lonb, latb, pref, axes, time)

      call atmos_tracer_driver_init ( lonb, latb, trs, axes, time, mask)

!-----------------------------------------------------------------------

      do_init = .false.

!-----------------------------------------------------------------------

 end subroutine physics_driver_init

!#######################################################################

 subroutine physics_driver_end (Time)

   type(time_type), intent(in) :: Time
!-----------------------------------------------------------------------

      call vert_turb_driver_end
      call vert_diff_driver_end
      call radiation_driver_end

      if(.not.do_mcm_moist_process) then
        call  moist_processes_end
      endif

      call atmos_tracer_driver_end
      call   damping_driver_end

!-----------------------------------------------------------------------

 end subroutine physics_driver_end

!#######################################################################
 function mcm_moist_process()
 logical :: mcm_moist_process

 if(do_init) then
   call error_mesg('function mcm_moist_process', &
                   'physics_driver not initialized', FATAL)
 endif
 mcm_moist_process = do_mcm_moist_process

 return
 end function mcm_moist_process
!#######################################################################

 subroutine check_args (lat, lon, area, p_half, p_full, z_half, z_full,&
                        u, v, t, q, r, um, vm, tm, qm, rm,             &
                        udt, vdt, tdt, qdt, rdt, mask, kbot)
!-----------------------------------------------------------------------
  real,    intent(in), dimension(:,:)     :: lat, lon, area
  real,    intent(in), dimension(:,:,:)   :: p_half, p_full,   &
                                             z_half, z_full,   &
                                             u, v, t, q, um, vm, tm, qm
  real,    intent(in), dimension(:,:,:,:) :: r, rm
  real,    intent(in), dimension(:,:,:)   :: udt, vdt, tdt, qdt
  real,    intent(in), dimension(:,:,:,:) :: rdt

  real,    intent(in), dimension(:,:,:),optional :: mask
  integer, intent(in), dimension(:,:)  ,optional :: kbot
!-----------------------------------------------------------------------
  integer :: id,jd,kd,nt,ierr

      id=size(u,1); jd=size(u,2); kd=size(u,3); nt=size(r,4)

      ierr = 0
      ierr = ierr + check_dim (lat, 'lat',  id,jd)
      ierr = ierr + check_dim (lon, 'lon',  id,jd)
      ierr = ierr + check_dim (area,'area', id,jd)

      ierr = ierr + check_dim (p_half,'p_half', id,jd,kd+1)
      ierr = ierr + check_dim (p_full,'p_full', id,jd,kd)
      ierr = ierr + check_dim (z_half,'z_half', id,jd,kd+1)
      ierr = ierr + check_dim (z_full,'z_full', id,jd,kd)

      ierr = ierr + check_dim (u, 'u',  id,jd,kd)
      ierr = ierr + check_dim (v, 'v',  id,jd,kd)
      ierr = ierr + check_dim (t, 't',  id,jd,kd)
      ierr = ierr + check_dim (q, 'q',  id,jd,kd)
      ierr = ierr + check_dim (um,'um', id,jd,kd)
      ierr = ierr + check_dim (vm,'vm', id,jd,kd)
      ierr = ierr + check_dim (tm,'tm', id,jd,kd)
      ierr = ierr + check_dim (qm,'qm', id,jd,kd)

      ierr = ierr + check_dim (udt,'udt', id,jd,kd)
      ierr = ierr + check_dim (vdt,'vdt', id,jd,kd)
      ierr = ierr + check_dim (tdt,'tdt', id,jd,kd)
      ierr = ierr + check_dim (qdt,'qdt', id,jd,kd)

!------------------------
     if (nt > 0) then
!------------------------
      ierr = ierr + check_dim (r,  'r',   id,jd,kd,nt)
      ierr = ierr + check_dim (rm, 'rm',  id,jd,kd,nt)
      ierr = ierr + check_dim (rdt,'rdt', id,jd,kd,nt)
!------------------------
     endif
!------------------------

      if (ierr > 0) then
          call error_mesg ('physics_driver', 'bad dimensions', FATAL)
      endif

      do_check_args = .false.

!-----------------------------------------------------------------------

      end subroutine check_args

!#######################################################################

      function check_dim_2d (data,name,id,jd) result (ierr)

      real,    intent(in), dimension(:,:) :: data
      character(len=*), intent(in)        :: name
      integer, intent(in)                 :: id,jd
      integer  ierr

           ierr = 0
      if (size(data,1) /= id) then
           call error_mesg ('physics_driver',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
           call error_mesg ('physics_driver',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif

      end function check_dim_2d

!#######################################################################

      function check_dim_3d (data,name,id,jd,kd) result (ierr)

      real,    intent(in), dimension(:,:,:) :: data
      character(len=*), intent(in)          :: name
      integer, intent(in)                   :: id,jd,kd
      integer  ierr

           ierr = 0
      if (size(data,1) /= id) then
            call error_mesg ('physics_driver',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
           call error_mesg ('physics_driver',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
           call error_mesg ('physics_driver',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif

      end function check_dim_3d

!#######################################################################

      function check_dim_4d (data,name,id,jd,kd,nt) result (ierr)

      real,    intent(in), dimension(:,:,:,:) :: data
      character(len=*), intent(in)            :: name
      integer, intent(in)                     :: id,jd,kd,nt
      integer  ierr

           ierr = 0
      if (size(data,1) /= id) then
            call error_mesg ('physics_driver',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
           call error_mesg ('physics_driver',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
           call error_mesg ('physics_driver',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif
      if (size(data,4) /= nt) then
           call error_mesg ('physics_driver',  &
                'dimension 4 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif

      end function check_dim_4d

!#######################################################################
 
 
 
end module physics_driver_mod
 
 
