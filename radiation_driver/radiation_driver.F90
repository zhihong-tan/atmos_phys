
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
                                  print_version_number, get_my_pe,     &
                                  read_data, write_data, close_file

   use          clouds_mod, only:  clouds, clouds_init, clouds_end

   use    diag_manager_mod, only:  register_diag_field, send_data

   use    time_manager_mod, only:  time_type, increment_time,      &
                                   set_date, set_time, get_time,   &
                                   get_calendar_type, operator(+), &
                                   operator(>), operator(>=),      &
                                   operator(-), operator(/=)

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
      character(len=4), parameter :: version_number = 'v2.1'

!   ---- list of restart versions readable by this module ----
      integer, dimension(1) :: restart_versions = (/ 1 /)
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
!   tdt_rad = the current radiative (sw + lw) heating rate

    real, allocatable, dimension(:,:,:) :: tdt_rad
    real, allocatable, dimension(:,:)   :: flux_sw_surf, flux_lw_surf

!-----------------------------------------------------------------------
!   ------------------- time step constant --------------------------

           integer :: num_pts, total_pts
   type(time_type) :: Next_rad_time, Rad_time_step

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
!!!   real :: solar_constant = 1.96  !(1367.44w/m2)
      real :: solar_constant = 1365.
   integer, dimension(6) :: rad_date = (/ 0, 0, 0, 0, 0, 0 /)
      real :: co2std = 330., ratio = 1.0

namelist /radiation_driver_nml/ dt_rad, offset, do_average, &
                                do_diurnal, do_annual,      &
                                do_clear_sky_diag,          &
                                solar_constant, rad_date,   &
                                co2std, ratio

!-----------------------------------------------------------------------
!-------------------- diagnostics fields -------------------------------

integer :: id_alb_sfc, id_tot_cld_amt, id_cld_amt, id_em_cld,  &
           id_alb_uv_cld, id_alb_nir_cld, id_abs_uv_cld, id_abs_nir_cld

integer, dimension(2) :: id_tdt_sw,   id_tdt_lw,  &
                         id_swdn_toa, id_swup_toa, id_olr, &
                         id_swup_sfc, id_swdn_sfc,         &
                         id_lwup_sfc, id_lwdn_sfc

character(len=9), parameter :: mod_name = 'radiation'

real :: missing_value = -999.

!-----------------------------------------------------------------------

contains

!#######################################################################
!#######################################################################

   subroutine radiation_driver (is, ie, js, je, Time, dt,             &
                                lat, lon, pfull, phalf, t, q, ts,     &
                                land, albedo, tdt, flux_sw, flux_lw,  &
                                mask, kbot)

!-----------------------------------------------------------------------
!
!    in:  is,ie,js,je  starting/ending i,j in global storage arrays;
!         Time       time for determining whether radiation should be
!                    done 
!         dt         time step in seconds, interval that 
!                      radiation_driver is called (not the interval
!                      that radiation calculations are done)
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
!
! optional in:
!         mask
!         kbot
!
!-----------------------------------------------------------------------
        integer, intent(in)            :: is, ie, js, je
type(time_type), intent(in)            :: Time
   real, intent(in)                    :: dt
   real, intent(in), dimension(:,:)    :: lat, lon
   real, intent(in), dimension(:,:,:)  :: pfull, phalf, t, q
   real, intent(in), dimension(:,:)    :: ts
   real, intent(in), dimension(:,:)    :: land
   real, intent(in), dimension(:,:)    :: albedo
   real, intent(inout), dimension(:,:,:) :: tdt
   real, intent(out),   dimension(:,:)   :: flux_sw, flux_lw

   real, intent(in), dimension(:,:,:),optional   :: mask
integer, intent(in), dimension(:,:),optional     :: kbot
  
!-----------------------------------------------------------------------

      if (do_init) call error_mesg ('radiation_driver',  &
                   'radiation_driver_init must first be called.', FATAL)

!-----------------------------------------------------------------------
!------------ is this a radiative time step ??? ------------------------

      if ( Time >= Next_rad_time ) then
            do_rad = .true.
            num_pts = num_pts + size(t,1) * size(t,2)
!           --- increment alarm when all points have been processed ---
            if (num_pts == total_pts) then
                Next_rad_time = Next_rad_time + Rad_time_step
                num_pts = 0
            endif
      else
            do_rad = .false.
      endif

!-----------------------------------------------------------------------
!--------- set up radiative input data and/or --------------------------
!          compute new radiative tendencies

   if (do_rad .or. do_average) then
       call radiation_calc (is, ie, js, je, Time, dt, lat, lon,    &
                            pfull, phalf, t, q, ts,                &
                            land, albedo, mask, kbot)
   endif

!-----------------------------------------------------------------------
!------------ update radiative tendency and fluxes ---------------------

         tdt  (:,:,:) = tdt(:,:,:) + tdt_rad(is:ie,js:je,:)

         flux_sw(:,:) =         flux_sw_surf(is:ie,js:je)
         flux_lw(:,:) =         flux_lw_surf(is:ie,js:je)

         
!-----------------------------------------------------------------------
!------------ save radiative tendency for use in -----------------------
!------------ strat cloud scheme... Note that if -----------------------
!------------ strat cloud scheme is not operating ----------------------
!------------ then nothing is done by this routine ---------------------

         call add_strat_tend(is,ie,js,je,tdt_rad(is:ie,js:je,:))
         

!-----------------------------------------------------------------------

      end subroutine radiation_driver

!#######################################################################

  subroutine radiation_calc (is, ie, js, je, Time, dt, lat, lon,     &
                             pfull, phalf, t, q, ts, land, albedo,   &
                             mask, kbot)

!-----------------------------------------------------------------------
        integer, intent(in)               :: is, ie, js, je
type(time_type), intent(in)               :: Time
           real, intent(in)               :: dt
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
   real,dimension(size(t,1),size(t,2),size(t,3))  :: cloud
   real,dimension(size(t,1),size(t,2))   :: swin,swout,olr, &
                                            swups,swdns,lwups,lwdns
integer,dimension(size(t,1),size(t,2))   :: nclds
   real,dimension(size(t,1),size(t,2))   :: tca
   real,dimension(size(t,1),size(t,2))   :: hang,cosz,solar,  &
                                            psfc,tsfc,asfc

integer :: i, j, k, n, id, jd, kd, ipass, npass, ip, jp, kb
logical :: no_clouds, used
type(time_type) :: Rad_time, Dt_zen

!-----------------------------------------------------------------------

      id = size(t,1);  jd = size(t,2);  kd = size(t,3)

!-----------------------------------------------------------------------
!------------- date related stuff --------------------------------------
  
      if (use_rad_date) then
          Rad_time = set_date (rad_date(1), rad_date(2), rad_date(3),  &
                               rad_date(4), rad_date(5), rad_date(6))
      else
          Rad_time = Time
      endif

      
!------------- solar stuff ---------------------------------------------

      if (do_diurnal) then
         if (do_average) then
            Dt_zen = set_time(int(dt+0.5),0)
         else
            Dt_zen = Rad_time_step
         endif
         call diurnal_solar     (cosz, solar, lat,lon, Rad_time, Dt_zen)
      else if (do_annual) then
         call annual_mean_solar (cosz, solar, lat)
      else
         call daily_mean_solar  (cosz, solar, lat, Rad_time)
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
      
      if ( id_tot_cld_amt > 0 ) then


           call clouds (is,js, no_clouds, Rad_time, lat, land,   &
                        press(:,:,1:kd), p_int, temp(:,:,1:kd),   &
                        cosz,nclds, ktopsw, kbtmsw, ktop, kbtm,   &
                        cldamt, cuvrf, cirrf, cirab, emcld, kbot, tca)
           
      else

          call clouds (is,js, no_clouds, Rad_time, lat, land,   &
                        press(:,:,1:kd), p_int, temp(:,:,1:kd),   &
                        cosz,nclds, ktopsw, kbtmsw, ktop, kbtm,   &
                        cldamt, cuvrf, cirrf, cirab, emcld, kbot)
 
      endif
      
!-----------------------------------------------------------------------
!----------------------------- radiation -------------------------------

     if (js == -35) then
        ip=20; jp=1
     else
        ip=0; jp=0
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

!del      radflux(is:ie,js:je)=(swups(:,:)-swdns(:,:))+(-lwdns(:,:))

!-----------------------------------------------------------------------
!------------------------ diagnostics section --------------------------

!------- albedo (only do once) -------------------------
      if ( id_alb_sfc > 0 .and. .not. no_clouds ) then
         used = send_data ( id_alb_sfc, asfc, Time, is, js )
      endif

!------- TOTAL CLOUD AMOUNT (only do once) -----------------
      if ( id_tot_cld_amt > 0 .and. .not. no_clouds ) then
         tca = tca*100.
         used = send_data ( id_tot_cld_amt, tca, Time, is, js )
      endif

!------- cloud amount (only do once ?) -------------------------
      if ( id_cld_amt > 0 .and. .not. no_clouds ) then
!        -- insert clouds (use n+1 offset) --
         cloud=0.0
         do j=1,jd; do i=1,id
            do n=1,nclds(i,j)
              cloud(i,j,ktop(i,j,n+1):kbtm(i,j,n+1)) = cldamt(i,j,n+1)
            enddo
         enddo; enddo
         used = send_data ( id_cld_amt, cloud, Time, is, js, 1, &
                            rmask=mask )
      endif

!------- cloud emissivity ---------------------------------------
      if ( id_em_cld > 0 .and. .not. no_clouds ) then
!        -- insert emissivities (use n+1 offset) --
         cloud=0.0
         do j=1,jd; do i=1,id
            do n=1,nclds(i,j)
              cloud(i,j,ktop(i,j,n+1):kbtm(i,j,n+1)) = emcld(i,j,n+1)
            enddo
         enddo; enddo
         used = send_data ( id_em_cld, cloud, Time, is, js, 1, &
                            rmask=mask )
      endif

!------- ultra-violet reflected by cloud -----------------------------
      if ( id_alb_uv_cld > 0 .and. .not. no_clouds ) then
!        ---- insert emissivities (use n+1 offset) ----
         cloud=0.0
         do j=1,jd; do i=1,id
            do n=1,nclds(i,j)
              cloud(i,j,ktop(i,j,n+1):kbtm(i,j,n+1)) = cuvrf(i,j,n+1)
            enddo
         enddo; enddo
         used = send_data ( id_alb_uv_cld, cloud, Time, is, js, 1, &
                            rmask=mask )
      endif

!------- infra-red reflected by cloud -----------------------------
      if ( id_alb_nir_cld > 0 .and. .not. no_clouds ) then
!        ---- insert emissivities (use n+1 offset) ----
         cloud=0.0
         do j=1,jd; do i=1,id
            do n=1,nclds(i,j)
              cloud(i,j,ktop(i,j,n+1):kbtm(i,j,n+1)) = cirrf(i,j,n+1)
            enddo
         enddo; enddo
         used = send_data ( id_alb_nir_cld, cloud, Time, is, js, 1, &
                            rmask=mask )
      endif

!------- ultra-violet absorbed by cloud (not implemented)------------
!     if ( id_abs_uv_cld > 0 .and. .not. no_clouds ) then
!        ---- insert emissivities (use n+1 offset) ----
!        cloud=0.0
!        do j=1,jd; do i=1,id
!           do n=1,nclds(i,j)
!             cloud(i,j,ktop(i,j,n+1):kbtm(i,j,n+1)) = cuvab(i,j,n+1)
!           enddo
!        enddo; enddo
!        used = send_data ( id_abs_uv_cld, cloud, Time, is, js, 1, &
!                           rmask=mask )
!     endif

!------- infra-red absorbed by cloud -----------------------------
      if ( id_abs_nir_cld > 0 .and. .not. no_clouds ) then
!        ---- insert emissivities (use n+1 offset) ----
         cloud=0.0
         do j=1,jd; do i=1,id
            do n=1,nclds(i,j)
              cloud(i,j,ktop(i,j,n+1):kbtm(i,j,n+1)) = cirab(i,j,n+1)
            enddo
         enddo; enddo
         used = send_data ( id_abs_nir_cld, cloud, Time, is, js, 1, &
                            rmask=mask )
      endif

!------- sw tendency -----------
      if ( id_tdt_sw(ipass) > 0 ) then
         used = send_data ( id_tdt_sw(ipass), tdtsw, Time, is, js, 1, &
                            rmask=mask )
      endif

!------- lw tendency -----------
      if ( id_tdt_lw(ipass) > 0 ) then
         used = send_data ( id_tdt_lw(ipass), tdtlw, Time, is, js, 1, &
                            rmask=mask )
      endif

!------- incoming sw flux toa -------
      if ( id_swdn_toa(ipass) > 0 ) then
         used = send_data ( id_swdn_toa(ipass), swin, Time, is, js )
      endif

!------- outgoing sw flux toa -------
      if ( id_swup_toa(ipass) > 0 ) then
         used = send_data ( id_swup_toa(ipass), swout, Time, is, js )
      endif

!------- outgoing lw flux toa (olr) -------
      if ( id_olr(ipass) > 0 ) then
         used = send_data ( id_olr(ipass), olr, Time, is, js )
      endif

!------- upward sw flux surface -------
      if ( id_swup_sfc(ipass) > 0 ) then
         used = send_data ( id_swup_sfc(ipass), swups, Time, is, js )
      endif

!------- downward sw flux surface -------
      if ( id_swdn_sfc(ipass) > 0 ) then
         used = send_data ( id_swdn_sfc(ipass), swdns, Time, is, js )
      endif

!------- upward lw flux surface -------
      if ( id_lwup_sfc(ipass) > 0 ) then
         used = send_data ( id_lwup_sfc(ipass), lwups, Time, is, js )
      endif

!------- downward lw flux surface -------
      if ( id_lwdn_sfc(ipass) > 0 ) then
         used = send_data ( id_lwdn_sfc(ipass), lwdns, Time, is, js )
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
      integer :: unit, io, ierr, id, jd, kd, next(2), dt(2), cal, vers
      logical :: do_avg_init, end
      type(time_type) :: Old_time_step, New_rad_time
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
    call print_version_number (unit, 'radiation_driver', version_number)
    if ( get_my_pe() == 0 ) write (unit, nml=radiation_driver_nml)
    call close_file (unit)


!-------- set default alarm info --------

     if (offset > 0) then
        Next_rad_time = increment_time (Time, offset, 0)
     else
        Next_rad_time = increment_time (Time, dt_rad, 0)
     endif
        Rad_time_step = set_time (dt_rad,0)
        Old_time_step = set_time (dt_rad,0)

!-----------------------------------------------------------------------

         call rdparm_init (kd)
         call co2_data (co2std, ratio, pref)

         allocate (tdt_rad     (id,jd,kd))
         allocate (flux_sw_surf(id,jd))
         allocate (flux_lw_surf(id,jd))

         call clouds_init (lonb,latb)

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
            read (unit) next,dt,cal
!           ---- override alarm with restart values ----
               Old_time_step = set_time (dt(1),dt(2))
            if ( cal == get_calendar_type() ) then
               Next_rad_time = set_time (next(1),next(2))
            else
!           ---- calendar changed, will reset alarm ----
               call error_mesg ('radiation_driver_init', &
                   'current calendar not equal restart calendar', NOTE)
            endif

            call read_data (unit, tdt_rad)
            call read_data (unit, flux_sw_surf)
            call read_data (unit, flux_lw_surf)
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
!--------- only used for initial guess of sea ice temp ------
            tdt_rad = 0.0
            flux_sw_surf = 50.0
            flux_lw_surf = 50.0
            if (do_average) do_avg_init = .true.
         endif

!-----------------------------------------------------------------------
!------------- initialize input data averages --------------------------

         if (do_avg_init) then
            nsum=0; psum=0.; tsum=0.; qsum=0.; asum=0.; csum=0.; ssum=0.
         endif

!-----------------------------------------------------------------------
!--------- adjust radiation alarm if alarm interval has changed --------

     if ( Rad_time_step /= Old_time_step ) then
         New_rad_time = Next_rad_time - Old_time_step + Rad_time_step
         if ( New_rad_time > Time ) then
            call error_mesg ('radiation_driver_init',  &
                     'radiation time step has changed, &
                     &next radiation time also changed', NOTE)
            Next_rad_time = New_rad_time
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

    integer :: unit, next(2), dt(2), cal

    unit = open_file (file='RESTART/radiation_driver.res', &
                      form='native', action='write')

!   --- single threading of restart output ---
    if (get_my_pe() == 0) then
!      --- write last value in list of restart versions ---
       write (unit) restart_versions(size(restart_versions))

!      --- write alarm info ---
       call get_time (Next_rad_time,next(1),next(2))
       call get_time (Rad_time_step, dt (1), dt (2))
       cal = get_calendar_type()
       write (unit) next,dt,cal
    endif

!   --- data ---
    call write_data (unit, tdt_rad)
    call write_data (unit, flux_sw_surf)
    call write_data (unit, flux_lw_surf)

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

    call clouds_end

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

    id_tot_cld_amt = &
    register_diag_field ( mod_name, 'tot_cld_amt', axes(1:2), Time, &
                         'total cloud amount', 'percent'            )

    id_cld_amt = &
    register_diag_field ( mod_name, 'cld_amt', axes(1:3), Time, &
                         'cloud amount', 'percent',             &
                         missing_value=missing_value            )

    id_em_cld = &
    register_diag_field ( mod_name, 'em_cld', axes(1:3), Time, &
                         'cloud emissivity', 'percent',        &
                          missing_value=missing_value          )

    id_alb_uv_cld = &
    register_diag_field ( mod_name, 'alb_uv_cld', axes(1:3), Time, &
                         'UV reflected by cloud', 'percent',       &
                          missing_value=missing_value              )

    id_alb_nir_cld = &
    register_diag_field ( mod_name, 'alb_nir_cld', axes(1:3), Time, &
                         'IR reflected by cloud', 'percent',        &
                          missing_value=missing_value               )

!   --- do not output this field ---
!   id_abs_uv_cld = &
!   register_diag_field ( mod_name, 'abs_uv_cld', axes(1:3), Time, &
!                        'UV absorbed by cloud', 'percent',        &
!                         missing_value=missing_value              )

    id_abs_nir_cld = &
    register_diag_field ( mod_name, 'abs_nir_cld', axes(1:3), Time, &
                         'IR absorbed by cloud', 'percent',         &
                          missing_value=missing_value               )


!-----------------------------------------------------------------------

   end subroutine diag_field_init

!#######################################################################

end module radiation_driver_mod

