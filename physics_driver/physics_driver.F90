
module physics_driver_mod

!-----------------------------------------------------------------------

use constants_new_mod,    only: constants_new_init

use sea_esf_rad_mod,      only: sea_esf_rad_init,  &
                                sea_esf_rad, &
                                sea_esf_rad_time_vary, &
                                sea_esf_rad_end, &
                                get_lrad

use rad_utilities_mod,     only: Environment, environment_type

use astronomy_package_mod, only: astronomy_driver, astronomy_dealloc, &
                                 astronomy

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

use     time_manager_mod, only: time_type, get_date_julian,  &
                                set_date_julian, get_time, &
                                set_time, operator (-), operator (/=)

use    tracer_driver_mod, only: tracer_driver_init, tracer_driver,  &
                                tracer_driver_end


implicit none
private

!-----------------------------------------------------------------------
!------- namelist ------

character(len=16)    :: rad_step_physics='default'


namelist / physics_driver_nml /    &
                                 rad_step_physics



!-----------------------------------------------------------------------
!  -------- public interfaces ---------

public  physics_driver_down, physics_driver_up,   &
        physics_driver_init, physics_driver_end
public  surf_diff_type

!-----------------------------------------------------------------------
!  -------- private interfaces/data ----------

interface check_dim
     module procedure check_dim_2d, check_dim_3d, check_dim_4d
end interface

!--------------------- version number ----------------------------------
character(len=256) :: version = '$Id: physics_driver.F90,v 1.6 2001/10/25 17:48:12 fms Exp $'
character(len=256) :: tag = '$Name: fez $'
!-----------------------------------------------------------------------

      logical :: do_init = .true., do_check_args = .true.
      logical :: do_sea_esf_physics

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine physics_driver_down (is, ie, js, je,                     &
                                 Time_prev, Time, Time_next,         &
                                 lat, lon, area,                     & 
                                 p_half, p_full, z_half, z_full,     &
                                 u, v, t, q, r, um, vm, tm, qm, rm,  &
                                 frac_land, rough_mom,               &
                                 albedo,    t_surf_rad,              &
                                 u_star,    b_star,                  &
                                 dtau_dv,  tau_x,  tau_y,            &
                                 udt, vdt, tdt, qdt, rdt,            &
                                 flux_sw,  flux_lw,  coszen,  gust,  &
                                 Surf_diff,                          &
                                 mask, kbot                          )

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
      real, intent(inout),dimension(:,:,:,:) :: rdt

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

   type(time_type)  :: Rad_time, Set_dat5, Delta_t, Tim_minus
   logical          :: do_rad, do_average
   integer          :: dt_rad
   integer          :: year, month, day, hour, min2, sec, jld, jld2
   integer          :: sec_d, idt
   real             :: dt
   real             :: fjd, fjulan


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

!******************************************************************

!--------------------------------------------------------------------
!  the following is needed to implement sea_esf radiation-step physics.
!--------------------------------------------------------------------
      if (do_sea_esf_physics) then

!--------------------------------------------------------------------
!  call get_lrad to : a)determine whether this is a radiation time step,
!  b) determine if averaging of input data is desired, c) determine 
!  the radiation time-step, and d) determine the next calendar time at 
!  which radiation is to be calculated.
!--------------------------------------------------------------------
        call get_lrad (Time, do_rad, do_average, dt_rad, Rad_time)

!-------------------------------------------------------------------
!  define the julian date (julian day number and fraction).
!--------------------------------------------------------------------
        if (Environment%using_sky_periphs) then
          call get_date_julian (Rad_time, year, month, day, hour,  &
                                min2, sec)
          Set_dat5 = set_date_julian (year, month, day, hour, min2, sec)
          call get_time(Set_dat5, sec_d  , jld )
          fjd = float(sec_d)/86400.0
!-------------------------------------------------------------------
!  define the time of the last physics time step.
!--------------------------------------------------------------------
          idt =nint(dt)
          Delta_t = set_time(idt, 0)
          Tim_minus = Rad_time - Delta_t

!-------------------------------------------------------------------
!  define the julian date of the previous timestep(fjulan).
!--------------------------------------------------------------------
          call get_date_julian (Tim_minus, year, month, day, hour, &
                                min2, sec)
          Set_dat5 = set_date_julian (year, month, day, hour, min2, sec)
          call get_time(Set_dat5, sec_d, jld2 )
          fjulan = float(sec_d)/86400.0

!--------------------------------------------------------------------
!  call astronomy to calculate solar orbital parameters for the rad-
!  iation time desired. convert fms julian day to skyhi julian day
!  number.
!--------------------------------------------------------------------
          jld = jld + 1721411
          fjd = fjd - 0.5
          if (fjd < 0.0) then
            fjd = fjd + 1.0
            jld = jld - 1
          endif
          fjulan = fjulan - 0.5
          if (fjulan < 0.0) then
            fjulan = fjulan + 1.0
          endif
          call astronomy (do_rad, fjd, jld)
        endif

!-------------------------------------------------------------------
!  call sea_esf_rad_time_vary to handle the temporal variation of
!  any other radiative inputs (radiative gases, solar constant).
!-------------------------------------------------------------------
        call sea_esf_rad_time_vary (Rad_time)
 
!------------------------------------------------------------------
!  call astronomy driver to calculate zenith angles and daylight frac-
!  tion at the current time, if a radiation step, or, if not, at the 
!  next radiation time, and then call sea_esf_rad to do the calculations
!  and/or update the temperature tendency, then deallocate the astron-
!  omy variables in astronomy_dealloc.
!------------------------------------------------------------------
        if (Environment%using_sky_periphs) then
          call astronomy_driver (is, ie, js, je, dt, dt_rad, do_rad, &
		                 Time=Time, lat=lat, lon=lon,   &
			         Rad_time=Rad_time,   &
			         do_average=do_average, fjd=fjd,  &
			         fjulan=fjulan )
	else
          call astronomy_driver (is, ie, js, je, dt, dt_rad, do_rad, &
	 	                 Time=Time, lat=lat, lon=lon,   &
			         Rad_time=Rad_time,   &
				 do_average=do_average)
        endif

        call sea_esf_rad (is, ie, js, je,   &
		 	  p_full, p_half, t, q, t_surf_rad, Time=Time, &
			  Time_next=Time_next,  lat=lat, lon=lon,  &
			  land=frac_land, albedo_in=albedo, &
			  tdt=tdt, flux_sw=flux_sw, &
			  flux_lw=flux_lw, mask=mask, kbot=kbot, &
			  coszen = coszen)
        call astronomy_dealloc

!-------------------------------------------------------------------
!  if fms radiation-step physics is to be done, call radiation_driver.
!-------------------------------------------------------------------

      else

	call radiation_driver (is, ie, js, je, Time, Time_next, lat,  &
			       lon, p_full, p_half, t, q, t_surf_rad,&
			       frac_land, albedo, tdt, flux_sw, &
			       flux_lw,  coszen,     &
			       mask=mask, kbot=kbot)
      endif

!******************************************************************

!-----------------------------------------------------------------------
!------------------------ damping --------------------------------------

     call damping_driver (is, js, Time_next, dt,                   &
                          p_full, p_half, z_full, z_half,          &
                          u, v, t, q, r,  udt, vdt, tdt, qdt, rdt, &
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

  call tracer_driver (lon, lat, (frac_land > 0.5), p_half, r, rdt, kbot)

!-----------------------------------------------------------------------
!----------------------- vertical diffusion ----------------------------

  call vert_diff_driver_down (is, js, Time_next, dt, p_half, z_full, &
                              diff_m, diff_t, um ,vm ,tm ,qm ,rm,    &
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

    call moist_processes (is, ie, js, je, Time_next, dt, frac_land, &
                          p_half, p_full, omega,                    &
                          t, q, r, u, v, tm, qm, rm, um, vm,        &
                          tdt, qdt, rdt, udt, vdt,                  &
                          lprec, fprec, mask=mask, kbot=kbot        )

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


    if (rad_step_physics == 'default') then
      do_sea_esf_physics = .false.
    else if (rad_step_physics == 'sea_esf') then
      do_sea_esf_physics = .true.
    else 
      call error_mesg ('sea_esf_rad_init',  &
	'rad_step_physics is not a valid value', FATAL)
    endif

!-----------------------------------------------------------------------
!------------- initialization for various schemes ----------------------

      call  moist_processes_init ( id, jd, kd, axes, Time )

      call vert_turb_driver_init ( id, jd, kd, axes, Time )

      call vert_diff_driver_init ( Surf_diff, id, jd, kd, axes, Time )

      call   damping_driver_init ( lonb, latb, axes, Time )

      if (do_sea_esf_physics) then
        call constants_new_init
        call sea_esf_rad_init (kd, axes, lonb=lonb, latb=latb,    &
		      	       pref=pref, timesince=time)
      else
        call radiation_driver_init (lonb, latb, pref, axes, time)
      endif

      call    tracer_driver_init ( grav, trs, mask)

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

      if (do_sea_esf_physics) then
        call sea_esf_rad_end
      else
        call radiation_driver_end
      endif

      call  moist_processes_end
      call    tracer_driver_end
      call   damping_driver_end

!-----------------------------------------------------------------------

 end subroutine physics_driver_end

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
 
 
