
                   module radiation_driver_mod

!-----------------------------------------------------------------------
!
!             interface module for radiation
!
!-----------------------------------------------------------------------

   use       fsrad_mod, only: rdparm_init, fsrad, co2_data
   use   astronomy_mod, only: diurnal_solar, daily_mean_solar,  &
                              annual_mean_solar
   use zonal_ozone_mod, only: zonal_ozone

   use       constants_mod, only: p00
   use       utilities_mod, only: file_exist, error_mesg, FATAL, NOTE, &
                                  check_nml_error, open_file,          &
                                  read_data, write_data, close_file,   &
                                  get_my_pe, get_domain_decomp

   use          clouds_mod, only:  clouds, clouds_init, clouds_end

   use    diag_manager_mod, only:  register_diag_field, send_data

   use    time_manager_mod, only:  time_type, set_date, set_time,  &
                                   get_time,    operator(+),       &
                                   operator(-), operator(/=), get_date

   use     strat_cloud_mod, only:  add_strat_tend

   use   diag_integral_mod, only:     diag_integral_field_init, &
                                  sum_diag_integral_field

   implicit none
   private

!----------- public interfaces in this module -----------------------

   public    radiation_driver, radiation_driver_init, &
             radiation_driver_end

!-----------------------------------------------------------------------
!------------ version number for this module ---------------------------
character(len=128) :: version = '$Id: radiation_driver.F90,v 1.5 2001/07/05 17:26:45 fms Exp $'
character(len=128) :: tag = '$Name: eugene $'

!   ---- list of restart versions readable by this module ----
!   (sorry, but restart version 1 will not be readable by this module)
      integer, dimension(1) :: restart_versions = (/ 2 /)
!-----------------------------------------------------------------------

     real, parameter :: p00inv=1./p00

     logical :: do_init=.true.

!-----------------------------------------------------------------------
!   the following code converts rvco2, the volume mixing ratio of co2
!   (used by lwrad) into rco2, the mass mixing ratio (used by swrad).

    real, parameter :: ratco2mw=1.519449738

    real, parameter :: rvco2 = 3.3e-4
    real, parameter :: rco2 = rvco2*ratco2mw


!--------- flags -------------------------------------------------------

    logical ::  do_rad
    logical ::  use_rad_date

!-----------------------------------------------------------------------
!  ------- allocatable global data saved by this module -------
!
!   tdt_rad      = radiative (sw + lw) heating rate
!   flux_sw_surf = net (down-up) sw flux at surface
!   flux_lw_surf = downward lw flux at surface
!   coszen_angle = cosine of the zenith angle 
!                  (used for the last radiation calculation)

    real, allocatable, dimension(:,:,:) :: tdt_rad
    real, allocatable, dimension(:,:)   :: flux_sw_surf, &
                                           flux_lw_surf, &
                                           coszen_angle

!-----------------------------------------------------------------------
!   ------------------- time step constant --------------------------
!
!   rad_alarm     = time interval until the next radiation calculation
!                   (integer in seconds)
!   rad_time_step = radiation time step (integer in seconds)
!   total_pts = number of grid boxes to be processed every time step
!               (note: all grid boxes must be processed every time step)
!   num_pts   = counter for current number of grid boxes processed
!             (when num_pts=0 or num_pts=total_pts certain things happen)

           integer :: num_pts, total_pts
           integer :: rad_alarm, rad_time_step

!-----------------------------------------------------------------------
!------- private allocatable array for time averaged input data --------

      real, allocatable, dimension(:,:,:) :: psum,tsum,qsum
      real, allocatable, dimension(:,:)   :: asum,csum,ssum
   integer, allocatable, dimension(:,:)   :: nsum

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------
!
!   dt_rad  =  radiative time step in seconds
!
!   offset  =  offset for radiative time step (in seconds) 
!              note if offset=1 radition will be done on the first
!              time step after t=0, dt_rad, 2*dt_rad, etc.
!              for now do not change this value
!
!   solar_constant = solar constant in watts/m2
!
!   rad_date       = fixed date (year,month,day,hour,min,sec)
!                    for radiation (solar info, ozone, clouds)
!
!   co2std = standard co2 vol. mixing ratio (either 300 or 330 ppmv)
!
!   ratio  = co2 vol. mixing ratio in units of the standard vol.
!            mixing ratio (must lie between 0.5 and 4.0)


   integer :: dt_rad=43200, offset=1
   logical :: do_diurnal=.false., do_annual=.false.,  &
              do_clear_sky_diag=.false., do_average=.false.
   logical :: rsd=.false.
!!!   real :: solar_constant = 1.96  !(1367.44w/m2)
      real :: solar_constant = 1365.
   integer, dimension(6) :: rad_date = (/ 0, 0, 0, 0, 0, 0 /)
      real :: co2std = 330., ratio = 1.0
      integer :: jpt = -35, ipt = -35

namelist /radiation_driver_nml/ dt_rad, offset, do_average, &
                                do_diurnal, do_annual,      &
                                do_clear_sky_diag,          &
                                solar_constant, rad_date,   &
                                co2std, ratio, jpt, ipt, &
				rsd

!-----------------------------------------------------------------------
!-------------------- diagnostics fields -------------------------------

integer :: id_alb_sfc, id_coszen

integer, dimension(2) :: id_tdt_sw,   id_tdt_lw,  &
                         id_swdn_toa, id_swup_toa, id_olr, &
                         id_swup_sfc, id_swdn_sfc,         &
                         id_lwup_sfc, id_lwdn_sfc

character(len=9), parameter :: mod_name = 'radiation'

real :: missing_value = -999.
integer  :: x(4), y(4)

!-----------------------------------------------------------------------

contains

!#######################################################################
!#######################################################################

   subroutine radiation_driver (is, ie, js, je, Time, Time_next,      &
                                lat, lon, pfull, phalf, t, q, ts,     &
                                land, albedo, tdt, flux_sw, flux_lw,  &
                                coszen, mask, kbot)

!-----------------------------------------------------------------------
!
!    in:  is,ie,js,je  starting/ending i,j in global storage arrays;
!         Time       time for determining whether radiation should be done
!         Time_next  time used for diagnostic output
!         lat        latitude,in radians
!         lon        longitude,in radians
!         pfull      pressure at full levels
!         phalf      pressure at half levels
!         t          temperature in model layers (full layers)
!         q          specific humidity of water vapor in model layers
!         ts         surface temperature in deg k
!         land       fraction of surface which covered by land (real)
!         albedo     surface albedo (current only output)
!
! inout:  tdt        temperature tendency
!
!   out:  flux_sw    net shortwave surface flux (down-up) (w/m2)
!         flux_lw    downward longwave surface flux (w/m2)
!         coszen     cosine of the zenith angle used for the most
!                    recent radiation calculation
!
! optional in:
!         mask
!         kbot
!
!  NOTE:   "dt" the time step in seconds, is the interval that 
!          radiation_driver is called (not the interval that radiation
!          calculations are done.
!          It may be calculated as: dt = Time_next - Time
!-----------------------------------------------------------------------
        integer, intent(in)            :: is, ie, js, je
type(time_type), intent(in)            :: Time, Time_next
   real, intent(in), dimension(:,:)    :: lat, lon
   real, intent(in), dimension(:,:,:)  :: pfull, phalf, t, q
   real, intent(in), dimension(:,:)    :: ts
   real, intent(in), dimension(:,:)    :: land
   real, intent(in), dimension(:,:)    :: albedo
   real, intent(inout), dimension(:,:,:) :: tdt
   real, intent(out),   dimension(:,:)   :: flux_sw, flux_lw, coszen

   real, intent(in), dimension(:,:,:),optional   :: mask
integer, intent(in), dimension(:,:),optional     :: kbot
!-----------------------------------------------------------------------
integer :: day, sec, dt
 
!-----------------------------------------------------------------------

      if (do_init) call error_mesg ('radiation_driver',  &
                   'radiation_driver_init must first be called.', FATAL)

!-----------------------------------------------------------------------
!------------ is this a radiative time step ??? ------------------------

      call get_time (Time_next-Time,sec,day)
      dt = day*86400 + sec

      if (dt <= 0) call error_mesg ('radiation_driver_mod', &
	         	              'Time_next <= Time', FATAL)

      if ( num_pts == 0 ) rad_alarm = rad_alarm - dt
      num_pts = num_pts + size(t,1) * size(t,2)

      if ( rad_alarm <= 0 ) do_rad = .true.

!-----------------------------------------------------------------------
!--------- set up radiative input data and/or --------------------------
!          compute new radiative tendencies

   if (do_rad .or. do_average) then
       call radiation_calc (is, ie, js, je, Time, Time_next, dt,   &
			      lat, lon,  pfull, phalf, t, q, ts,     &
                            land, albedo, mask, kbot)
   endif

!-----------------------------------------------------------------------
!------------ update radiative tendency and fluxes ---------------------

         tdt  (:,:,:) = tdt(:,:,:) + tdt_rad(is:ie,js:je,:)

         flux_sw(:,:) =         flux_sw_surf(is:ie,js:je)
         flux_lw(:,:) =         flux_lw_surf(is:ie,js:je)
         coszen (:,:) =         coszen_angle(is:ie,js:je)

         
!-----------------------------------------------------------------------
!------------ save radiative tendency for use in -----------------------
!------------ strat cloud scheme... Note that if -----------------------
!------------ strat cloud scheme is not operating ----------------------
!------------ then nothing is done by this routine ---------------------

         call add_strat_tend(is,ie,js,je,tdt_rad(is:ie,js:je,:))
         

!-----------------------------------------------------------------------
!--- reset alarm when all points have been processed ---

    if ( num_pts == total_pts ) then
        num_pts = 0
        if ( do_rad ) then
             rad_alarm = rad_alarm + rad_time_step
             do_rad = .false.
        endif
    endif

!-----------------------------------------------------------------------

      end subroutine radiation_driver

!#######################################################################

  subroutine radiation_calc (is, ie, js, je, Time, Time_diag, dt, lat, lon, &
                             pfull, phalf, t, q, ts, land, albedo,   &
                             mask, kbot)

!-----------------------------------------------------------------------
        integer, intent(in)               :: is, ie, js, je
type(time_type), intent(in)               :: Time, Time_diag
        integer, intent(in)               :: dt
   real, intent(in),  dimension(:,:)    :: lat,lon
   real, intent(in),  dimension(:,:,:)  :: pfull, phalf, t, q
   real, intent(in),  dimension(:,:)    :: ts
   real, intent(in),  dimension(:,:)    :: land
   real, intent(in),  dimension(:,:)    :: albedo

   real, intent(in), dimension(:,:,:),optional   :: mask
integer, intent(in), dimension(:,:),optional     :: kbot
!-----------------------------------------------------------------------
!---------------- local arrays & variables -----------------------------

   real,dimension(size(t,1),size(t,2), size(t,3)+1) :: press,temp,  &
                                                       cldamt,emcld,  &
                                                       cuvrf,cirrf,  &
                                                       cirab, phaf,  &
                                                       p_int
   real,dimension(size(t,1),size(t,2), size(t,3))   :: rh2o,qo3
integer,dimension(size(t,1),size(t,2), size(t,3)+1) :: ktop,kbtm,  &
                                                       ktopsw,kbtmsw

   real,dimension(size(t,1),size(t,2),size(t,3))  :: tdtsw,tdtlw
   real,dimension(size(t,1),size(t,2))   :: swin,swout,olr, &
                                            swups,swdns,lwups,lwdns
integer,dimension(size(t,1),size(t,2))   :: nclds
   real,dimension(size(t,1),size(t,2))   :: hang,cosz,solar,  &
                                            psfc,tsfc,asfc,cosz1

integer :: i, j, k, n, id, jd, kd, ipass, npass, ip, jp, kb
!dummy variable and time of day
integer :: dum, tod(3)
logical :: no_clouds, used
type(time_type) :: Rad_time, Dt_zen

!-----------------------------------------------------------------------

      id = size(t,1);  jd = size(t,2);  kd = size(t,3)

!-----------------------------------------------------------------------
!------------- date related stuff --------------------------------------
  
      if (rsd) then

        if (.not. use_rad_date) call error_mesg ('radiation_calc', &  
              'if (rsd), must set rad_date(1:3) to valid date', FATAL)

        call get_date(Time,dum,dum,dum,tod(1),tod(2),tod(3))

        Rad_time = set_date (rad_date(1), rad_date(2), rad_date(3),   &
                tod(1),tod(2),tod(3))
	
      else if (use_rad_date) then
          Rad_time = set_date (rad_date(1), rad_date(2), rad_date(3),  &
                               rad_date(4), rad_date(5), rad_date(6))
      else
          Rad_time = Time
      endif

      
!------------- solar stuff ---------------------------------------------

      if (do_diurnal) then
         if (do_average) then
            Dt_zen = set_time(dt,0)
         else
            Dt_zen = set_time(rad_time_step,0)
         endif
         call diurnal_solar (cosz1, solar, lat,lon, Rad_time+Dt_zen, Dt_zen)
         call diurnal_solar (cosz , solar, lat,lon, Rad_time       , Dt_zen)
      else if (do_annual) then
         call annual_mean_solar (cosz, solar, lat)
	   cosz1 = cosz
      else
         call daily_mean_solar  (cosz, solar, lat, Rad_time)
	   cosz1 = cosz
      endif

      solar = solar * solar_constant

!-----------------------------------------------------------------------
!     ---- surface data ----

   if (present(kbot)) then
      do j=1,jd
      do i=1,id
         kb=kbot(i,j)
         psfc(i,j)=phalf(i,j,kb+1)
      enddo
      enddo
   else
         psfc(:,:)=phalf(:,:,kd+1)
   endif
         tsfc(:,:)=ts(:,:)
         asfc(:,:)=albedo(:,:)

!-----------------------------------------------------------------------

!     ---- data -----
      do k=1,kd
         press(:,:,k)=pfull(:,:,k)
         temp (:,:,k)=t(:,:,k)
         rh2o (:,:,k)=q(:,:,k)/(1.0-q(:,:,k))
      enddo
         press(:,:,kd+1)=phalf(:,:,kd+1)
         temp (:,:,kd+1)=tsfc(:,:)

!  ---- if kbot is present, fill in underground ----
!          temperatures with surface value

   if (present(kbot)) then
      do j=1,jd
      do i=1,id
         kb=kbot(i,j)
         if (kb < kd) then
            do k=kb+1,kd
               temp(i,j,k)=temp(i,j,kd+1)
            enddo
         endif
      enddo
      enddo
   endif

!-----------------------------------------------------------------------
!---------------- set up for averaging data ----------------------------

      if (do_average) then
          call compute_average (is,ie,js,je, &
                                press,temp,rh2o,asfc,cosz,solar)

          if (.not.do_rad) then
              return
          else
              call return_average (is,ie,js,je,  &
                                   press,temp,rh2o,asfc,cosz,solar)
          endif
      endif

!---------------- limit mixing ratio of water vapor --------------------
         where (rh2o(:,:,:) < 3.e-6) rh2o(:,:,:)=3.e-6

!-------- compute pres at layer interface as done in radiation code ----
         p_int(:,:,1)   =0.0
         p_int(:,:,2:kd)=0.5*(press(:,:,1:kd-1)+press(:,:,2:kd))
         p_int(:,:,kd+1)=press(:,:,kd+1)

!-----------------------------------------------------------------------
!     ---- data for ozone ---

      do k=1,kd+1
         phaf(:,:,k)=p_int(:,:,k)
         phaf(:,:,k)=101325.*phaf(:,:,k)/p_int(:,:,kd+1)
      enddo

         qo3  (:,:,:)=1.e-10

         call zonal_ozone (Rad_time, lat, phaf, qo3)

!-----------------------------------------------------------------------
!--------------- loop for clear sky diagnostics option -----------------

      if (do_clear_sky_diag) then
         npass=2
      else
         npass=1
      endif
	

   do ipass = 1, npass

!-----------------------------------------------------------------------
!     ---- data for clouds ---

                          no_clouds = .true.
      if (ipass == npass) no_clouds = .false.
      
      call clouds (is,js, no_clouds, Rad_time, Time_diag, lat, land, tsfc, &
		      press(:,:,1:kd), p_int, temp(:,:,1:kd), q,  &
                    cosz,nclds, ktopsw, kbtmsw, ktop, kbtm,   &
                    cldamt, cuvrf, cirrf, cirab, emcld, mask, kbot)
 
!-----------------------------------------------------------------------
!----------------------------- radiation -------------------------------

!----------------------------------------------------------------------
!   determine if radiation diagnostics column is present in current 
!   physics window. if so, determine its coordinates in the 
!   physics_window space. if not present, set ip and jp to 0.
!----------------------------------------------------------------------
     if ( (x(3) <= ipt .and. x(4) >= ipt) .and.    &
	  (y(3) <= jpt .and. y(4) >= jpt) ) then
       ip = ipt -x(3) + 1
       jp = jpt -y(3) + 1
       if (js > jp .or. je < jp ) then
         ip = 0 ; jp = 0
       else if (is > ip .or. ie < ip ) then
         ip = 0 ; jp = 0
       else
         ip = ip - is + 1; jp = jp - js + 1
       endif
     else
       ip =0 ; jp = 0
     endif

     if (present(kbot)) then
         call fsrad (ip,jp,press,temp,rh2o,qo3,  &
                     nclds,ktopsw,kbtmsw,ktop,kbtm,cldamt,  &
                     emcld,cuvrf,cirrf,cirab,asfc,rvco2,  &
                     cosz,solar,                            &
                     swin,swout,olr,swups,swdns,lwups,lwdns,  &
                     tdtsw,tdtlw, kbot,psfc)
     else
         call fsrad (ip,jp,press,temp,rh2o,qo3,  &
                     nclds,ktopsw,kbtmsw,ktop,kbtm,cldamt,  &
                     emcld,cuvrf,cirrf,cirab,asfc,rvco2,  &
                     cosz,solar,                            &
                     swin,swout,olr,swups,swdns,lwups,lwdns,  &
                     tdtsw,tdtlw)
     endif

!-----------------------------------------------------------------------
!------------------------ save global data -----------------------------

     if (present(mask)) then
          tdt_rad(is:ie,js:je,:)=(tdtsw(:,:,:)+tdtlw(:,:,:))*mask(:,:,:)
     else
          tdt_rad(is:ie,js:je,:)=(tdtsw(:,:,:)+tdtlw(:,:,:))
     endif

       flux_sw_surf(is:ie,js:je) = swdns(:,:)-swups(:,:)
       flux_lw_surf(is:ie,js:je) = lwdns(:,:)
       coszen_angle(is:ie,js:je) = cosz1(:,:)

!-----------------------------------------------------------------------
!------------------------ diagnostics section --------------------------

!------- albedo (only do once) -------------------------
      if ( id_alb_sfc > 0 .and. .not. no_clouds ) then
         used = send_data ( id_alb_sfc, asfc, Time_diag, is, js )
      endif

!------- cos zenith angle (only do once) -------------------------
      if ( id_coszen > 0 .and. .not. no_clouds ) then
         used = send_data ( id_coszen, cosz, Time_diag, is, js )
      endif

!------- sw tendency -----------
      if ( id_tdt_sw(ipass) > 0 ) then
         used = send_data ( id_tdt_sw(ipass), tdtsw, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- lw tendency -----------
      if ( id_tdt_lw(ipass) > 0 ) then
         used = send_data ( id_tdt_lw(ipass), tdtlw, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- incoming sw flux toa -------
      if ( id_swdn_toa(ipass) > 0 ) then
         used = send_data ( id_swdn_toa(ipass), swin, Time_diag, is, js )
      endif

!------- outgoing sw flux toa -------
      if ( id_swup_toa(ipass) > 0 ) then
         used = send_data ( id_swup_toa(ipass), swout, Time_diag, is, js )
      endif

!------- outgoing lw flux toa (olr) -------
      if ( id_olr(ipass) > 0 ) then
         used = send_data ( id_olr(ipass), olr, Time_diag, is, js )
      endif

!------- upward sw flux surface -------
      if ( id_swup_sfc(ipass) > 0 ) then
         used = send_data ( id_swup_sfc(ipass), swups, Time_diag, is, js )
      endif

!------- downward sw flux surface -------
      if ( id_swdn_sfc(ipass) > 0 ) then
         used = send_data ( id_swdn_sfc(ipass), swdns, Time_diag, is, js )
      endif

!------- upward lw flux surface -------
      if ( id_lwup_sfc(ipass) > 0 ) then
         used = send_data ( id_lwup_sfc(ipass), lwups, Time_diag, is, js )
      endif

!------- downward lw flux surface -------
      if ( id_lwdn_sfc(ipass) > 0 ) then
         used = send_data ( id_lwdn_sfc(ipass), lwdns, Time_diag, is, js )
      endif

!-----------------------------------------------------------------------

   enddo

!--------------- end of clear sky diagnostics loop ---------------------
!-----------------------------------------------------------------------

!     ---- accumulate global integral quantities -----

      call sum_diag_integral_field ('olr',    olr,        is, js)
      call sum_diag_integral_field ('abs_sw', swin-swout, is, js)

!-----------------------------------------------------------------------

   end subroutine radiation_calc

!#######################################################################

   subroutine radiation_driver_init ( lonb, latb, pref, axes, Time )

!-----------------------------------------------------------------------
           real, intent(in), dimension(:)   :: lonb, latb
           real, intent(in), dimension(:,:) :: pref
        integer, intent(in), dimension(4)   :: axes
type(time_type), intent(in)                 :: Time
!-----------------------------------------------------------------------
      integer :: unit, io, ierr, id, jd, kd, vers
      logical :: do_avg_init, end
      integer :: old_time_step, new_rad_time
      character(len=4) :: chvers

      id=size(lonb,1)-1; jd=size(latb,1)-1; kd=size(pref,1)-1

      if (size(pref,2) /= 2) call error_mesg ('radiation_driver_init', &
                    'must have two reference pressure profiles.', FATAL)

!       ----- read namelist -----

      if ( file_exist('input.nml')) then
         unit = open_file (file='input.nml', action='read')
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=radiation_driver_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'radiation_driver_nml')
         enddo
  10     call close_file (unit)
      endif

!      ----- write namelist -----

    unit = open_file ('logfile.out', action='append')
    if ( get_my_pe() == 0 ) then
         write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
         write (unit, nml=radiation_driver_nml)
    endif
    call close_file (unit)

    if (do_annual .and. do_diurnal) then
      call error_mesg ( 'radiation_driver_init', &
         'Cannot have both do_annual and do_diurnal as .true.', FATAL)
    endif

    
!--------------------------------------------------------------------
!  retrieve global and subdomain indices
!--------------------------------------------------------------------
      call get_domain_decomp (x,y)

!-------- set default alarm info --------

     if (offset > 0) then
        rad_alarm = offset
     else
        rad_alarm = dt_rad
     endif
        rad_time_step = dt_rad
        old_time_step = dt_rad

!-----------------------------------------------------------------------

         call rdparm_init (kd)
         call co2_data (co2std, ratio, pref)

         allocate (tdt_rad     (id,jd,kd))
         allocate (flux_sw_surf(id,jd))
         allocate (flux_lw_surf(id,jd))
         allocate (coszen_angle(id,jd))

         call clouds_init ( lonb, latb, axes, Time )

!-----------------------------------------------------------------------
!-------------- initialize time averaged input data --------------------

         if (do_average) then
             allocate (psum(id,jd,kd+1), tsum(id,jd,kd+1),    &
                       qsum(id,jd,kd), asum(id,jd), csum(id,jd),  &
                       ssum(id,jd), nsum(id,jd))
         endif
         do_avg_init = .false.

!-----------------------------------------------------------------------
!-------------- restart data -------------------------------------------

         if (file_exist('INPUT/radiation_driver.res')) then
            unit = open_file (file='INPUT/radiation_driver.res',  &
                              form='native', action='read')

!           --- read and check restart version number ---
            read (unit) vers
            if ( .not. any(vers == restart_versions) ) then
              write (chvers,'(i4)') vers
              call error_mesg ('radiation_driver_init', &
                    'restart version '//chvers//' cannot be read &
                     &by this module version', FATAL)
            endif

!           --- reading alarm info ---
!           ---- override alarm with restart values ----
            read (unit) rad_alarm, old_time_step

            call read_data (unit, tdt_rad)
            call read_data (unit, flux_sw_surf)
            call read_data (unit, flux_lw_surf)
            call read_data (unit, coszen_angle)
                if (do_average) then
                   call read_data (unit, nsum, end)
                   if (end) go to 11
                   call read_data (unit, psum)
                   call read_data (unit, tsum)
                   call read_data (unit, qsum)
                   call read_data (unit, asum)
                   call read_data (unit, csum)
                   call read_data (unit, ssum)
                   goto 12
  11               do_avg_init = .true.
  12               continue
                endif
            call close_file (unit)
         else
!-------- initial surface flux set to 100 wm-2 ---------
!---    was only used for initial guess of sea ice temp ?  -----
!--- in current model framework the ice model is called after radiation
            tdt_rad = 0.0
            flux_sw_surf = 50.0
            flux_lw_surf = 50.0
            coszen_angle = 0.50
            if (do_average) do_avg_init = .true.
         endif

!-----------------------------------------------------------------------
!------------- initialize input data averages --------------------------

         if (do_avg_init) then
            nsum=0; psum=0.; tsum=0.; qsum=0.; asum=0.; csum=0.; ssum=0.
         endif

!-----------------------------------------------------------------------
!--------- adjust radiation alarm if alarm interval has changed --------

     if ( rad_time_step /= old_time_step ) then
         new_rad_time = rad_alarm - old_time_step + rad_time_step
         if ( new_rad_time > 0 ) then
            if (get_my_pe() == 0) call error_mesg      &
                    ('radiation_driver_init',          &
                     'radiation time step has changed, &
                     &next radiation time also changed', NOTE)
            rad_alarm = new_rad_time
         endif
     endif
        
            total_pts=id*jd
            num_pts=0

!-----------------------------------------------------------------------
!---------- check if optional radiative date should be used ------------

      if (rad_date(1) > 1900 .and.                        &
          rad_date(2) >   0  .and. rad_date(2) < 13 .and. &
          rad_date(3) >   0  .and. rad_date(3) < 32 ) then
                  use_rad_date = .true.
      else
                  use_rad_date = .false.
      endif

!-----------------------------------------------------------------------
!------- initialize quantities for integral package -------

      call diag_integral_field_init ('olr',    'f8.3')
      call diag_integral_field_init ('abs_sw', 'f8.3')

!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------

      call diag_field_init ( Time, axes )

!-----------------------------------------------------------------------
         do_init=.false.
!-----------------------------------------------------------------------

   end subroutine radiation_driver_init

!#######################################################################

 subroutine radiation_driver_end

    integer :: unit

!-----------------------------------------------------------------------

    unit = open_file (file='RESTART/radiation_driver.res', &
                      form='native', action='write')

!   --- single threading of restart output ---
    if (get_my_pe() == 0) then
!      --- write last value in list of restart versions ---
       write (unit) restart_versions(size(restart_versions))

!      --- write alarm info ---
       write (unit) rad_alarm, rad_time_step
    endif

!   --- data ---
    call write_data (unit, tdt_rad)
    call write_data (unit, flux_sw_surf)
    call write_data (unit, flux_lw_surf)
    call write_data (unit, coszen_angle)

!   --- optional time avg data ---
  if (do_average) then
    call write_data (unit, nsum)
    call write_data (unit, psum)
    call write_data (unit, tsum)
    call write_data (unit, qsum)
    call write_data (unit, asum)
    call write_data (unit, csum)
    call write_data (unit, ssum)
  endif
    call close_file (unit)

!   --- release space for allocatable data ---

      deallocate (tdt_rad, flux_sw_surf, flux_lw_surf, coszen_angle)
  if (do_average) then
      deallocate (nsum, psum, tsum, qsum, asum, csum, ssum)
  endif

!   --- terminate other modules ---

    call clouds_end

!-----------------------------------------------------------------------

  end subroutine radiation_driver_end

!#######################################################################

   subroutine compute_average (is,ie,js,je,  &
                               press,temp,rh2o,albedo,cosz,solar)

!-----------------------------------------------------------------------
    integer, intent(in)                   :: is,ie,js,je
    real,    intent(in), dimension(:,:,:) :: press,temp,rh2o
    real,    intent(in), dimension(:,:)   :: albedo,cosz,solar
!-----------------------------------------------------------------------

      nsum(is:ie,js:je)   = nsum(is:ie,js:je)   + 1
      psum(is:ie,js:je,:) = psum(is:ie,js:je,:) + press(:,:,:)
      tsum(is:ie,js:je,:) = tsum(is:ie,js:je,:) + temp (:,:,:)
      qsum(is:ie,js:je,:) = qsum(is:ie,js:je,:) + rh2o (:,:,:)
      asum(is:ie,js:je)   = asum(is:ie,js:je)   + albedo(:,:)
      csum(is:ie,js:je)   = csum(is:ie,js:je)   + cosz  (:,:)
      ssum(is:ie,js:je)   = ssum(is:ie,js:je)   + solar (:,:)

!-----------------------------------------------------------------------

   end subroutine compute_average

!#######################################################################

   subroutine return_average (is,ie,js,je,  &
                                 press,temp,rh2o,albedo,cosz,solar)

!-----------------------------------------------------------------------
    integer, intent(in)                    :: is,ie,js,je
    real,    intent(out), dimension(:,:,:) :: press,temp,rh2o
    real,    intent(out), dimension(:,:)   :: albedo,cosz,solar
!-----------------------------------------------------------------------
    real, dimension(size(press,1),size(press,2)) :: fsum
    integer  n, k

!--------- check for zero divid ---------
      n = count (nsum(is:ie,js:je) <= 0)
      if ( n > 0 ) then
          call error_mesg ('return_average in radiation_driver_mod',  &
                           'dividing average by zero.',  FATAL)
      endif

!--------------- compute averages --------------------------------------

      fsum(:,:) = 1.0 / float(nsum(is:ie,js:je))

      do k=1,size(press,3)
         press(:,:,k) = psum(is:ie,js:je,k) * fsum(:,:)
         temp (:,:,k) = tsum(is:ie,js:je,k) * fsum(:,:)
      enddo

      do k=1,size(rh2o,3)
         rh2o (:,:,k) = qsum(is:ie,js:je,k) * fsum(:,:)
      enddo

      albedo(:,:) = asum(is:ie,js:je) * fsum(:,:)
      cosz  (:,:) = csum(is:ie,js:je) * fsum(:,:)
      solar (:,:) = ssum(is:ie,js:je) * fsum(:,:)

!  ----- zero out sums -----
      nsum(is:ie,js:je)   = 0
      psum(is:ie,js:je,:) = 0.0
      tsum(is:ie,js:je,:) = 0.0
      qsum(is:ie,js:je,:) = 0.0
      asum(is:ie,js:je)   = 0.0
      csum(is:ie,js:je)   = 0.0
      ssum(is:ie,js:je)   = 0.0

!-----------------------------------------------------------------------

   end subroutine return_average

!#######################################################################

   subroutine diag_field_init ( Time, axes )

     type(time_type), intent(in) :: Time
     integer        , intent(in) :: axes(4)

     character(len=4) :: clr
     character(len=9) :: clr2
     integer :: i, n

!------------ initialize diagnostic fields in this module --------------

!---- generate names for clear sky diagnostic fields -----
!---- may want to not generate names when clear sky not used? ----

                       n= 1
if (do_clear_sky_diag) n= 2

do i = 1, n

    if ( i == n) then
       clr  = "    "
       clr2 = "          "
    else
       clr  = "_clr"
       clr2 = "clear sky "
    endif

    id_tdt_sw(i) = &
    register_diag_field ( mod_name, 'tdt_sw'//trim(clr), axes(1:3), Time, &
                     trim(clr2)//'temperature tendency for SW radiation', &
                     'deg_K/sec', missing_value=missing_value             )

    id_tdt_lw(i) = &
    register_diag_field ( mod_name, 'tdt_lw'//trim(clr), axes(1:3), Time, &
                     trim(clr2)//'temperature tendency for LW radiation', &
                     'deg_K/sec', missing_value=missing_value             )

    id_swdn_toa(i) = &
    register_diag_field ( mod_name, 'swdn_toa'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'SW flux down at TOA', &
                     'watts/m2', missing_value=missing_value               )

    id_swup_toa(i) = &
    register_diag_field ( mod_name, 'swup_toa'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'SW flux up at TOA', &
                     'watts/m2', missing_value=missing_value               )

    id_olr(i) = &
    register_diag_field ( mod_name, 'olr'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'outgoing longwave radiation', &
                     'watts/m2', missing_value=missing_value               )

    id_swup_sfc(i) = &
    register_diag_field ( mod_name, 'swup_sfc'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'SW flux up at surface', &
                     'watts/m2', missing_value=missing_value               )

    id_swdn_sfc(i) = &
    register_diag_field ( mod_name, 'swdn_sfc'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'SW flux down at surface', &
                     'watts/m2', missing_value=missing_value               )

    id_lwup_sfc(i) = &
    register_diag_field ( mod_name, 'lwup_sfc'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'LW flux up at surface', &
                     'watts/m2', missing_value=missing_value               )

    id_lwdn_sfc(i) = &
    register_diag_field ( mod_name, 'lwdn_sfc'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'LW flux down at surface', &
                     'watts/m2', missing_value=missing_value               )

enddo


    id_alb_sfc = &
    register_diag_field ( mod_name, 'alb_sfc', axes(1:2), Time, &
                         'surface albedo', 'percent'            )

    id_coszen = &
    register_diag_field ( mod_name, 'coszen', axes(1:2), Time, &
                         'cosine of zenith angle', 'none'       )


!-----------------------------------------------------------------------

   end subroutine diag_field_init

!#######################################################################

end module radiation_driver_mod

