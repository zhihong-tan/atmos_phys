
                    module moist_processes_mod

!-----------------------------------------------------------------------
!
!         interface module for moisture processes
!         ---------------------------------------
!             moist convective adjustment
!             relaxed arakawa-schubert
!             donner deep convection
!             large-scale condensation
!             stratiform prognostic cloud scheme 
!             rel humidity cloud scheme 
!             Diagnostic cloud scheme 
!
!-----------------------------------------------------------------------

use      donner_deep_mod,only: donner_deep_init,               &
                               donner_deep, donner_deep_end
use      moist_conv_mod, only: moist_conv, moist_conv_init
use     lscale_cond_mod, only: lscale_cond, lscale_cond_init
use  sat_vapor_pres_mod, only: lookup_es

use    time_manager_mod, only: time_type

use    diag_manager_mod, only: register_diag_field, send_data

use             fms_mod, only: file_exist, check_nml_error,    &
                               open_namelist_file, close_file, &
                               write_version_number,           &
                               mpp_pe, mpp_root_pe, stdlog,    &
                               error_mesg, FATAL, NOTE

use             ras_mod, only: ras, ras_end, ras_init

use         dry_adj_mod, only: dry_adj, dry_adj_init

use     strat_cloud_mod, only: strat_cloud_init, strat_driv, strat_cloud_end, &
                               strat_cloud_sum

use       rh_clouds_mod, only: rh_clouds_init, rh_clouds_end, &
                               rh_clouds_sum

use      diag_cloud_mod, only: diag_cloud_init, diag_cloud_end, &
                               diag_cloud_sum

use diag_integral_mod, only:     diag_integral_field_init, &
                             sum_diag_integral_field

use       constants_mod, only: CP_AIR, GRAV, RDGAS, RVGAS

use     cu_mo_trans_mod, only: cu_mo_trans_init, cu_mo_trans, cu_mo_trans_end

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index,&
                              get_number_tracers, &
                              get_tracer_names, &
                              NO_TRACER
use atmos_tracer_utilities_mod, only : wet_deposition
use            fms_mod, only : mpp_clock_id, mpp_clock_begin, &
                               mpp_clock_end, CLOCK_MODULE, &
                               MPP_CLOCK_SYNC

! start chemistry modules
!use    chem_interface, only : id_wet
!use           mo_hook, only : moz_hook
! end chemistry modules
implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

   public   moist_processes, moist_processes_init, moist_processes_end, &
            doing_strat

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------


      real,parameter :: epst=200.

   real, parameter :: d622 = RDGAS/RVGAS
   real, parameter :: d378 = 1.-d622

   integer :: nsphum, nql, nqi, nqa   ! tracer indices for stratiform clouds

!--------------------- version number ----------------------------------
   character(len=128) :: version = '$Id: moist_processes.F90,v 10.0 2003/10/24 22:00:35 fms Exp $'
   character(len=128) :: tagname = '$Name: jakarta $'
   logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------
!-------------------- namelist data (private) --------------------------


   logical :: do_mca=.true., do_lsc=.true., do_ras=.false.,  &
              do_strat=.false., do_dryadj=.false., &
              do_rh_clouds=.false., do_diag_clouds=.false., &
              do_donner_deep=.false., do_cmt=.false., &
              use_tau=.false., do_gust_cv = .false.


   real :: pdepth = 150.e2
   real :: tfreeze = 273.16
   real :: gustmax = 3.             ! maximum gustiness wind (m/s)
   real :: gustconst = 10./86400.   ! constant in kg/m2/sec, default =
                                    ! 1 cm/day = 10 mm/day
                                    
!---------------- namelist variable definitions ------------------------
!
!   do_mca   = switch to turn on/off moist convective adjustment;
!                [logical, default: do_mca=true ]
!   do_lsc   = switch to turn on/off large scale condensation
!                [logical, default: do_lsc=true ]
!   do_ras   = switch to turn on/off relaxed arakawa shubert
!                [logical, default: do_ras=false ]
! do_donner_deep = switch to turn on/off donner deep convection scheme
!                [logical, default: do_donner_deep=false ]
!   do_strat = switch to turn on/off stratiform cloud scheme
!                [logical, default: do_strat=false ]
! do_rh_clouds = switch to turn on/off simple relative humidity cloud scheme
!                [logical, default: do_rh_clouds=false ]
! do_diag_clouds = switch to turn on/off (Gordon's) diagnostic cloud scheme
!                [logical, default: do_diag_clouds=false ]
!  do_dryadj = switch to turn on/off dry adjustment scheme
!                [logical, default: do_dryadj=false ]
!   use_tau  = switch to determine whether current time level (tau)
!                will be used or else future time level (tau+1).
!                if use_tau = true then the input values for t,q, and r
!                are used; if use_tau = false then input values
!                tm+tdt*dt, etc. are used.
!                [logical, default: use_tau=false ]
!
!   pdepth   = boundary layer depth in pascals for determining mean
!                temperature tfreeze (used for snowfall determination)
!   tfreeze  = mean temperature used for snowfall determination (deg k)
!                [real, default: tfreeze=273.16]
!
!  do_gust_cv = switch to use convective gustiness (default = false)
!  gustmax    = maximum convective gustiness (m/s)
!  gustconst  = precip rate which defines precip rate which begins to
!               matter for convective gustiness (kg/m2/sec)
!
!
!   notes: 1) do_lsc and do_strat cannot both be true
!          2) pdepth and tfreeze are used to determine liquid vs. solid
!             precipitation for mca, lsc, and ras schemes, the 
!             stratiform scheme determines it's own precipitation type.
!          3) if do_strat=true then stratiform cloud tracers: liq_wat,
!             ice_wat, cld_amt must be present 
!          4) do_donner_deep and do_rh_clouds cannot both be true
!             (pending revision of code flow)
!
!-----------------------------------------------------------------------

namelist /moist_processes_nml/ do_mca, do_lsc, do_ras, do_strat,  &
                               do_dryadj, pdepth, tfreeze,        &
                               use_tau, do_rh_clouds, do_diag_clouds, &
                               do_donner_deep, do_cmt, do_gust_cv, &
                               gustmax, gustconst

!-----------------------------------------------------------------------
!-------------------- diagnostics fields -------------------------------

integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_tdt_ls  , id_qdt_ls  , id_prec_ls  , id_snow_ls  , &
           id_precip  , id_WVP, id_LWP, id_IWP, id_AWP, id_gust_conv, &
           id_tdt_dadj, id_rh,  id_mc, id_mc_donner, id_mc_full, &
           id_tdt_deep_donner, id_tdt_mca_donner, id_qdt_deep_donner,  &
           id_qdt_mca_donner, id_prec_deep_donner, id_prec_mca_donner,&
           id_snow_deep_donner, id_snow_mca_donner, &
           id_qldt_ls , id_qidt_ls , id_qldt_conv, id_qidt_conv, &
           id_qadt_ls , id_qadt_conv,id_ql_ls_col, id_qi_ls_col, &
           id_ql_conv_col, id_qi_conv_col, id_qa_ls_col, id_qa_conv_col,&
           id_q_conv_col, id_q_ls_col, id_t_conv_col, id_t_ls_col, &
           id_tracerdt_conv(100), id_tracer_conv_col(100), id_prod_no
! The dimension for id_tracerdt_conv and id_tracer_conv_col is the maximum 
! number of tracers allowed. Due to predecessor issues we cannot do a use 
! tracer_manager, only : max_tracers however. If the number of tracers 
! allowed is increased then this value should also be increased.

character(len=5) :: mod_name = 'moist'

real :: missing_value = -999.
integer :: convection_clock, largescale_clock

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine moist_processes (is, ie, js, je, Time, dt, land,            &
                            phalf, pfull, zhalf, zfull, omega, diff_t, &
                            radturbten,                                &
                            t, q, r, u, v, tm, qm, rm, um, vm,         &
                            tdt, qdt, rdt, udt, vdt, diff_cu_mo,       &
                            convect, lprec, fprec, gust_cv, area,      &
                            lat, mask, kbot)

!-----------------------------------------------------------------------
!
!    in:  is,ie      starting and ending i indices for window
!
!         js,je      starting and ending j indices for window
!
!         Time       time used for diagnostics [time_type]
!
!         dt         time step (from t(n-1) to t(n+1) if leapfrog)
!                    in seconds   [real]
!
!         land        fraction of surface covered by land
!                     [real, dimension(nlon,nlat)]
!
!         phalf      pressure at half levels in pascals
!                      [real, dimension(nlon,nlat,nlev+1)]
!
!         pfull      pressure at full levels in pascals
!                      [real, dimension(nlon,nlat,nlev)]
!
!         omega      omega (vertical velocity) at full levels
!                    in pascals per second
!                      [real, dimension(nlon,nlat,nlev)]
!
!         diff_t     vertical diffusion coefficient for temperature
!                    and tracer (m*m/sec) on half levels
!                      [real, dimension(nlon,nlat,nlev)]
!
!         t, q       temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         r          tracer fields at full model levels,
!                    at the current time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         u, v,      zonal and meridional wind [m/s] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
! 
!         tm, qm     temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the previous time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rm         tracer fields at full model levels,
!                    at the previous time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         um, vm     zonal and meridional wind [m/s] at full model levels,
!                    at the previous time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev)]
!
!         area       grid box area (in m2)
!                      [real, dimension(nlon,nlat)]
!
!         lat        latitude in radians
!                      [real, dimension(nlon,nlat)]
!  
! inout:  tdt, qdt   temperature (tdt) [deg k/sec] and specific
!                    humidity of water vapor (qdt) tendency [1/sec]
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rdt        tracer tendencies 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         udt, vdt   zonal and meridional wind tendencies [m/s/s]
! 
!   out:  convect    is moist convection occurring in this grid box?
!                   [logical, dimension(nlon,nlat)]
!
!         lprec      liquid precipitiaton rate (rain) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
!
!         fprec      frozen precipitation rate (snow) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
! 
!         gust_cv    gustiness from convection  in m/s
!                      [real, dimension(nlon,nlat)]
!
!       optional
!  -----------------
! 
!    in:  mask       mask (1. or 0.) for grid boxes above or below
!                    the ground   [real, dimension(nlon,nlat,nlev)]
!
!         kbot       index of the lowest model level
!                      [integer, dimension(nlon,nlat)]
!
!
!-----------------------------------------------------------------------
integer,         intent(in)              :: is,ie,js,je
type(time_type), intent(in)              :: Time
   real, intent(in)                      :: dt
   real, intent(in) , dimension(:,:)     :: land
   real, intent(in) , dimension(:,:,:)   :: phalf, pfull, zhalf, zfull,&
                                            omega, diff_t,             &
                                            t, q, u, v, tm, qm, um, vm
   real, dimension(:,:,:), intent(in)    :: radturbten
   real, intent(in) , dimension(:,:,:,:) :: r, rm
   real, intent(inout),dimension(:,:,:)  :: tdt, qdt, udt, vdt
   real, intent(inout),dimension(:,:,:,:):: rdt
logical, intent(out), dimension(:,:)     :: convect
   real, intent(out), dimension(:,:)     :: lprec, fprec, gust_cv
   real, intent(out), dimension(:,:,:)   :: diff_cu_mo
   real, intent(in) , dimension(:,:)     :: area
   real, intent(in) , dimension(:,:)     :: lat

   
   real, intent(in) , dimension(:,:,:), optional :: mask
integer, intent(in) , dimension(:,:),   optional :: kbot

!-----------------------------------------------------------------------
real, dimension(size(t,1),size(t,2),size(t,3)) :: tin,qin,ttnd,qtnd, &
                 rlin, riin, rnew, rlnew, rinew, qnew, qinew, qlnew, &
                 rltnd, ritnd, rin, rtnd, ahuco, qrat
real, dimension(size(t,1),size(t,2),size(t,3)) :: ttnd_save,qtnd_save, &
                                                  qitnd_save,   &
                                                  qltnd_save, qatnd_save
real, dimension(size(t,1),size(t,2))           :: tsnow,snow
real, dimension(size(t,1),size(t,2))           :: snow_save, rain_save
logical,dimension(size(t,1),size(t,2))         :: coldT
real, dimension(size(t,1),size(t,2),size(t,3)) :: utnd,vtnd,uin,vin

real, dimension(size(t,1),size(t,2),size(t,3)) :: qltnd,qitnd,qatnd, &
                                                  qlin, qiin, qain
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: tracer, tracertnd
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: mc,mask3
real, dimension(size(t,1),size(t,2),size(t,3)) :: mc_full, mc_donner
real, dimension(size(t,1),size(t,2),size(t,3)) :: RH, pmass, wetdeptnd
real, dimension(size(t,1),size(t,2))           :: rain, precip, precip_dd
real, dimension(size(t,1),size(t,2))           :: wvp,lwp,iwp
real, dimension(size(t,1),size(t,2))           :: tempdiag
real, dimension(size(t,1),size(t,2),size(t,3)) :: tempdiag1
integer n

integer :: i, j, k, ix, jx, kx, nt, ip, tr
real    :: dtinv
logical :: use_mask, do_adjust, used

real, dimension(size(t,1),size(t,2),size(t,3)+1) :: press
real, dimension(size(t,1),size(t,2),size(t,3)) :: tin_in

real :: tinsave

!chemistry start
real, dimension(size(rdt,1),size(rdt,2),size(rdt,3),size(rdt,4)) :: wet_data
real, dimension(size(rdt,1),size(rdt,2),size(rdt,3)) :: prod_no
integer, dimension(size(rdt,1),size(rdt,2)) :: cldtop, cldbot
real, parameter :: boltz = 1.38044e-16
real, dimension(size(t,1),size(t,2),size(t,3)) :: conc_air
!chemistry end
!-----------------------------------------------------------------------

! The following local quantitities are used exclusively for diagnostic clouds
!      LGSCLDELQ  Averaged rate of change in mix ratio due to lg scale precip 
!               at full model levels  
!               (dimensioned IX x JX x KX)
!      CNVCNTQ  Accumulated count of change in mix ratio due to conv precip 
!               at full model levels  
!               (dimensioned IX x JX x KX)
!      CONVPRC  Accumulated conv precip rate summed over all
!               full model levels (mm/day )
!               (dimensioned IX x JX)

real, dimension(size(t,1),size(t,2),size(t,3)) :: lgscldelq,cnvcntq
real, dimension(size(t,1),size(t,2)) :: convprc


!-----------------------------------------------------------------------

      if (.not. module_is_initialized) call error_mesg ('moist_processes',  &
                     'moist_processes_init has not been called.', FATAL)

   ! --- if do_ras is false, diff_cu_mo is never been initialized, which cause
   ! --- the model crash.
   diff_cu_mo = 0.0

!-----------------------------------------------------------------------
      conc_air = 10. * pfull / (boltz * t)

!-------- input array size and position in global storage --------------

      ix=size(t,1); jx=size(t,2); kx=size(t,3); nt=size(rdt,4)

!-----------------------------------------------------------------------

      use_mask=.false.
      if (present(mask) .and. present(kbot))  use_mask=.true.

      dtinv=1./dt
      lprec=0.0; fprec=0.0; precip=0.0

!------------------ setup input data -----------------------------------

   if (use_tau) then
      tin (:,:,:)=t(:,:,:)
      qin (:,:,:)=q(:,:,:)
      uin (:,:,:)=u(:,:,:)
      vin (:,:,:)=v(:,:,:)
   else
      tin (:,:,:)=tm(:,:,:)+tdt(:,:,:)*dt
      qin (:,:,:)=qm(:,:,:)+qdt(:,:,:)*dt
      uin (:,:,:)=um(:,:,:)+udt(:,:,:)*dt
      vin (:,:,:)=vm(:,:,:)+vdt(:,:,:)*dt
   endif

   if (use_mask) then
      tin (:,:,:)=mask(:,:,:)*tin(:,:,:)+(1.0-mask(:,:,:))*epst
      qin (:,:,:)=mask(:,:,:)*qin(:,:,:)
      uin (:,:,:)=mask(:,:,:)*uin(:,:,:)
      vin (:,:,:)=mask(:,:,:)*vin(:,:,:)
   endif
   
! Make local copies of all the tracers 
      if (use_tau) then
        do tr = 1, size(r,4)
          tracer (:,:,:,tr)=r(:,:,:,tr)
        enddo  
      else
        do tr = 1, size(r,4)
          tracer (:,:,:,tr)=rm(:,:,:,tr)+rdt(:,:,:,tr)*dt
        enddo  
      endif

      if (use_mask) then
        do tr = 1, size(r,4)
          tracer (:,:,:,tr)=mask(:,:,:)*tracer(:,:,:,tr)
        enddo  
      endif

!----------------   reset mc  ------------------------------------------
 
   mc(:,:,:) = 0.0
   mc_donner(:,:,:) = 0.0

!   initialize qrat and ahuco 

   qrat = 1.0
   ahuco = 0.0

!---compute mass in each layer if needed by any of the diagnostics ----- 

      if ( id_q_conv_col  > 0 .or. id_t_conv_col  > 0 .or. &
           id_q_ls_col    > 0 .or. id_t_ls_col    > 0 .or. &
           id_ql_conv_col > 0 .or. id_qi_conv_col > 0 .or. &
           id_qa_conv_col > 0 .or. id_ql_ls_col   > 0 .or. &
           id_qi_ls_col   > 0 .or. id_qa_ls_col   > 0 .or. &
           id_WVP         > 0 .or. id_LWP         > 0 .or. &
           id_IWP         > 0 .or. id_AWP         > 0 ) then
        do k=1,kx
          pmass(:,:,k) = (phalf(:,:,k+1)-phalf(:,:,k))/GRAV
        end do
      end if

!----------------- mean temp in lower atmosphere -----------------------
!----------------- used for determination of rain vs. snow -------------
!----------------- use input temp ? ------------------------------------

   call tempavg (pdepth, phalf, t, tsnow, mask)

   where (tsnow <= tfreeze)
          coldT=.TRUE.
   elsewhere
          coldT=.FALSE.
   endwhere
      
!-----------------------------------------------------------------------
!***********************************************************************
!----------------- dry adjustment scheme -------------------------------

call mpp_clock_begin( convection_clock )
if (do_dryadj) then

         call dry_adj (tin, pfull, phalf, ttnd, mask)
         tin=tin+ttnd
         ttnd=ttnd*dtinv
         tdt = tdt + ttnd
!------------- diagnostics for dt/dt_dry_adj ---------------------------
     if ( id_tdt_dadj > 0 ) then
        used = send_data ( id_tdt_dadj, ttnd, Time, is, js, 1, &
                           rmask=mask )
     endif
! ----------------------------------------------------------------------
end if

!---------------------------------------------------------------------
!   activate donner convection scheme
!---------------------------------------------------------------------

if (do_donner_deep) then


!--------------------------------------------------------------------
! convert specific humidity to mixing ratio for input to donner_deep
!--------------------------------------------------------------------
   if (do_strat) then

!  convert specific humidity to mixing ratio
     rin = qin/(1.0 - qin)
     rlin = tracer(:,:,:,nql)/(1.0 - qin)
     riin = tracer(:,:,:,nqi)/(1.0 - qin)
     call donner_deep (is, ie, js, je, tin, rin, pfull, phalf,   &
                       omega, dt, land, Time, ttnd, rtnd, precip_dd,   &
                       ahuco, qrat, &
                       kbot=kbot, cf=tracer(:,:,:,nqa),qlin=rlin , qiin=riin, &
                       delta_qa=tracertnd(:,:,:,nqa), delta_ql=rltnd, &
                       delta_qi=ritnd, mtot=mc_donner)
     where (coldT)
       snow = precip_dd
       rain = 0.0
     elsewhere
       snow = 0.0
       rain = precip_dd
     end where
                    
!   define updated mixing ratios
     rnew = rin + rtnd
     rlnew = rlin + rltnd
     rinew = riin + ritnd
!   define updated specific humidities
     qnew = rnew / (1.0 + rnew)
     qlnew = rlnew / (1.0 + rnew)
     qinew = rinew / (1.0 + rnew)
!   define specific humidity time tendencies
     qtnd = qnew - qin
     tracertnd(:,:,:,nql) = qlnew - tracer(:,:,:,nql)
     tracertnd(:,:,:,nqi) = qinew - tracer(:,:,:,nqi)

!   Update local tracer quantity for input to RAS later.
     tracer(:,:,:,nql) = tracer(:,:,:,nql) + tracertnd(:,:,:,nql)
     tracer(:,:,:,nqi) = tracer(:,:,:,nqi) + tracertnd(:,:,:,nqi)
     tracer(:,:,:,nqa) = tracer(:,:,:,nqa) + tracertnd(:,:,:,nqa)

   else

!  convert specific humidity to mixing ratio
     rin = qin/(1.0 - qin)
     call donner_deep (is, ie, js, je, tin, rin, pfull, phalf,    &
                       omega, dt, Land, Time, ttnd, rtnd, precip_dd,   &
                       ahuco, qrat,    &
                       kbot=kbot)
     where (coldT)
       snow = precip_dd
       rain = 0.0
     elsewhere
       snow = 0.0
       rain = precip_dd
     end where
                    
!   define updated mixing ratio
     rnew = rin + rtnd
!   define updated specific humidity
     qnew = rnew / (1.0 + rnew)
!   define specific humidity time tendency
     qtnd = qnew - qin
   endif

!ljd
!     do k=1,kx
!        if (qrat(is,js,k) .ne. 1.) then
!       if (js .eq. 1) then
!          write(6,*) 'is,js,k,qrat,ahuco= ',is,js,k,qrat(is,js,k), &
!          ahuco(is,js,k)
!          stop
!        end if
!        end if
!     end do
!ljd

!------- update input values and compute tendency -------
                    
      tin=tin+ttnd;    qin=qin+qtnd
      
      ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
      rain=rain*dtinv; snow=snow*dtinv

!------- add on tendency ----------

     tdt=tdt+ttnd; qdt=qdt+qtnd

!----------------------------------------------------------------
!  save tendencies and precip forms due to deep convection as diag-
!  nosed in donner_deep
!----------------------------------------------------------------
     if (id_tdt_deep_donner > 0) then
       used = send_data ( id_tdt_deep_donner, ttnd, Time, is, js, 1, &
                        rmask=mask )
     endif

     if (id_qdt_deep_donner > 0) then
       used = send_data ( id_qdt_deep_donner, qtnd, Time, is, js, 1, &
                        rmask=mask )
     endif

     if (id_prec_deep_donner > 0) then
       used = send_data ( id_prec_deep_donner, rain+snow, Time, is, js )
     endif

     if (id_snow_deep_donner > 0) then
       used = send_data ( id_snow_deep_donner, snow, Time, is, js)
     endif

     if ( id_mc_donner > 0 ) then
       used = send_data ( id_mc_donner, mc_donner, Time, is, js, 1,&
                         rmask=mask )
     endif


!------- update input values , compute and add on tendency -----------
!-------              in the case of strat                 -----------

      if (do_strat) then
         qlin=qlin+tracertnd(:,:,:,nql) 
         qiin=qiin+tracertnd(:,:,:,nqi) 
         qain=qain+tracertnd(:,:,:,nqa)

         tracertnd(:,:,:,nql)=tracertnd(:,:,:,nql)*dtinv
         tracertnd(:,:,:,nqi)=tracertnd(:,:,:,nqi)*dtinv
         tracertnd(:,:,:,nqa)=tracertnd(:,:,:,nqa)*dtinv

         rdt(:,:,:,nql)=rdt(:,:,:,nql)+tracertnd(:,:,:,nql)
         rdt(:,:,:,nqi)=rdt(:,:,:,nqi)+tracertnd(:,:,:,nqi)
         rdt(:,:,:,nqa)=rdt(:,:,:,nqa)+tracertnd(:,:,:,nqa)

      end if

!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

! save deep component of these diagnostics
      rain_save = rain
      snow_save = snow
      ttnd_save = ttnd
      qtnd_save = qtnd
      if (do_strat .and. do_ras) then
        qltnd_save = tracertnd(:,:,:,nql)
        qitnd_save = tracertnd(:,:,:,nqi)
        qatnd_save = tracertnd(:,:,:,nqa)
      endif

!--------------------------------------------------------------------
!  call moist convective adjustment to handle any shallow convection
!--------------------------------------------------------------------

      call moist_conv (tin,qin,pfull,phalf,coldT,&
                       ttnd,qtnd,rain,snow,kbot,&
                       .FALSE., tracer, tracertnd, dtinv, Time, mask, is, js )

!------- add on tendency ----------
      tdt=tdt+ttnd; qdt=qdt+qtnd

!---------------------------------------------------------------------
!   save the tendencies and precip from the moist convective 
!   adjustment pass associated with donner_deep
!---------------------------------------------------------------------
      if (id_tdt_mca_donner > 0) then
        used = send_data ( id_tdt_mca_donner, ttnd, Time, is, js, 1, &
                         rmask=mask )
      endif

      if (id_qdt_mca_donner > 0) then
        used = send_data ( id_qdt_mca_donner, qtnd, Time, is, js, 1, &
                         rmask=mask )
      endif

      if (id_prec_mca_donner > 0) then
        used = send_data ( id_prec_mca_donner, rain+snow, Time, is, js ) 
      endif

      if (id_snow_mca_donner > 0) then
        used = send_data ( id_snow_mca_donner, snow, Time, is, js)
      endif


!------- update input values , compute and add on tendency -----------
!-------              in the case of strat                 -----------

      if (do_strat) then
         
!         qlin=qlin+tracertnd(:,:,:,nql) 
!         qiin=qiin+tracertnd(:,:,:,nqi) 
!         qain=qain+tracertnd(:,:,:,nqa)

!         tracertnd(:,:,:,nql)=tracertnd(:,:,:,nql)*dtinv
!         tracertnd(:,:,:,nqi)=tracertnd(:,:,:,nqi)*dtinv
!         tracertnd(:,:,:,nqa)=tracertnd(:,:,:,nqa)*dtinv

!         rdt(:,:,:,nql)=rdt(:,:,:,nql)+tracertnd(:,:,:,nql)
!         rdt(:,:,:,nqi)=rdt(:,:,:,nqi)+tracertnd(:,:,:,nqi)
!         rdt(:,:,:,nqa)=rdt(:,:,:,nqa)+tracertnd(:,:,:,nqa)

      end if

!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

!--------------------------------------------------------------------
!  define sum of contributions from the deep pass and the moist conv
!  pass. store in temporary variables ..._save if ras will also be 
!  activated.
!--------------------------------------------------------------------
   if (.not. do_ras) then
      rain = rain_save + rain
      snow = snow_save + snow
      ttnd = ttnd_save + ttnd
      qtnd = qtnd_save + qtnd
    else
      rain_save = rain_save + rain
      snow_save = snow_save + snow
      ttnd_save = ttnd_save + ttnd
      qtnd_save = qtnd_save + qtnd
    endif

endif  ! do_donner_deep

!-----------------------------------------------------------------------
!***********************************************************************
!----------------- moist convective adjustment -------------------------

                        if (do_mca) then
!-----------------------------------------------------------------------

         call moist_conv (tin,qin,pfull,phalf,coldT,&
                          ttnd,qtnd,rain,snow,kbot,&
                          do_strat, tracer, tracertnd, &
                          dtinv, Time, mask, is, js )

!------- add on tendency ----------
     tdt=tdt+ttnd; qdt=qdt+qtnd

!------- update input values , compute and add on tendency -----------
!-------              in the case of strat                 -----------

      if (do_strat) then
         
         rdt(:,:,:,nql)=rdt(:,:,:,nql)+tracertnd(:,:,:,nql)
         rdt(:,:,:,nqi)=rdt(:,:,:,nqi)+tracertnd(:,:,:,nqi)
         rdt(:,:,:,nqa)=rdt(:,:,:,nqa)+tracertnd(:,:,:,nqa)

      end if

!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

!-----------------------------------------------------------------------
                           endif

!-----------------------------------------------------------------------
!***********************************************************************
!----------- relaxed arakawa/schubert cumulus param scheme -------------

                        if (do_ras) then
!-----------------------------------------------------------------------
      
      call ras (is,   js,     Time,     tin,   qin,   &
                uin,  vin,    pfull,    phalf, zhalf, coldT, &
                dt,   ttnd,   qtnd,     utnd,  vtnd,  &
                rain, snow,   do_strat, mask,  kbot,  &
                mc,   tracer, tracertnd )

!------- add on tendency ----------
      tdt=tdt+ttnd; qdt=qdt+qtnd
      udt=udt+utnd; vdt=vdt+vtnd

!------- update input values , compute and add on tendency -----------
!-------              in the case of strat                 -----------

      if (do_strat) then
        do tr = 1, size(r,4)
          if ( tr .ne. nsphum ) &
            rdt (:,:,:,tr) = rdt(:,:,:,tr) + tracertnd(:,:,:,tr)
        enddo  
      end if

!------- save total precip and snow ---------
      
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

! update convection diagnostics if  donner_deep was also active
      if (do_donner_deep) then
        ttnd = ttnd + ttnd_save
        qtnd = qtnd + qtnd_save
        rain = rain + rain_save
        snow = snow + snow_save
        if (do_strat) then
          tracertnd(:,:,:,nql) = tracertnd(:,:,:,nql) + qltnd_save
          tracertnd(:,:,:,nqi) = tracertnd(:,:,:,nqi) + qitnd_save
          tracertnd(:,:,:,nqa) = tracertnd(:,:,:,nqa) + qatnd_save
        endif
      endif
          
! do diffusive cumulus momentum transport
      if ( do_cmt ) then
        call cu_mo_trans ( is, js, Time, mc, tin, &
                           phalf, pfull, zhalf, zfull, diff_cu_mo )
      endif

!-----------------------------------------------------------------------
                           endif
! Do wet deposition for the convective routine here
! qdt should give the vertical distribution here
      wet_data = 0.0
      do n = 1, size(tracertnd,4)
        wetdeptnd = 0.0
        call wet_deposition(n, t, pfull, phalf, rain, snow, qtnd, tracer(:,:,:,n), &
                            wetdeptnd, Time, 'convect', is, js, dt)
        tracertnd(:,:,:,n) = tracertnd(:,:,:,n) - wetdeptnd
        wet_data(:,:,:,n) = wetdeptnd
      enddo

!-----------------------------------------------------------------------
! do convective gustiness

  
  gust_cv = 0.0
  if (do_gust_cv) then
     where((rain+snow).gt.0.)
        gust_cv = gustmax * sqrt( (rain+snow)/(gustconst+(rain+snow)) )
     endwhere
  end if

!-----------------------------------------------------------------------
! save diagnostic of convection

  convect = .false.
  where ( (rain+snow) .gt. 0. ) convect = .true.

!-----------------------------------------------------------------------
!***********************************************************************
!--------------- DIAGNOSTICS FOR CONVECTIVE SCHEME ---------------------
!-----------------------------------------------------------------------

!------- diagnostics for dt/dt_conv -------
      if ( id_tdt_conv > 0 ) then
        used = send_data ( id_tdt_conv, ttnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for dq/dt_conv -------
      if ( id_qdt_conv > 0 ) then
        used = send_data ( id_qdt_conv, qtnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for precip_conv -------
      if ( id_prec_conv > 0 ) then
        used = send_data ( id_prec_conv, rain+snow, Time, is, js )
      endif
!------- diagnostics for snow_conv -------
      if ( id_snow_conv > 0 ) then
        used = send_data ( id_snow_conv, snow, Time, is, js )
      endif

!------- diagnostics for gust_conv -------
      if ( id_gust_conv > 0 ) then
        used = send_data ( id_gust_conv, gust_cv, Time, is, js )
      endif

!------- diagnostics for water vapor path tendency ----------
      if ( id_q_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + qtnd(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_q_conv_col, tempdiag, Time, is, js )
      end if        
   
!------- diagnostics for dry static energy tendency ---------
      if ( id_t_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + ttnd(:,:,k)*CP_AIR*pmass(:,:,k)
        end do
        used = send_data ( id_t_conv_col, tempdiag, Time, is, js )
      end if        
   
   !------- stratiform cloud tendencies from cumulus convection ------------
   if ( do_strat ) then

      !------- diagnostics for dql/dt from convection -------
      if ( id_qldt_conv > 0 ) then
        used = send_data ( id_qldt_conv, tracertnd(:,:,:,nql), Time, is, js, 1, &
                           rmask=mask )
      endif
      
      !------- diagnostics for dqi/dt from convection -------
      if ( id_qidt_conv > 0 ) then
        used = send_data ( id_qidt_conv, tracertnd(:,:,:,nqi), Time, is, js, 1, &
                           rmask=mask )
      endif
      
      !------- diagnostics for dqa/dt from convection -------
      if ( id_qadt_conv > 0 ) then
        used = send_data ( id_qadt_conv, tracertnd(:,:,:,nqa), Time, is, js, 1, &
                           rmask=mask )
      endif

      !------- diagnostics for liquid water path tendency ------
      if ( id_ql_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + tracertnd(:,:,k,nql)*pmass(:,:,k)
        end do
        used = send_data ( id_ql_conv_col, tempdiag, Time, is, js )
      end if        
      
      !------- diagnostics for ice water path tendency ---------
      if ( id_qi_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + tracertnd(:,:,k,nqi)*pmass(:,:,k)
        end do
        used = send_data ( id_qi_conv_col, tempdiag, Time, is, js )
      end if        
      
      !---- diagnostics for column integrated cloud mass tendency ---
      if ( id_qa_conv_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + tracertnd(:,:,k,nqa)*pmass(:,:,k)
        end do
        used = send_data ( id_qa_conv_col, tempdiag, Time, is, js )
      end if        
         
      do n = 1, size(tracertnd,4)
      !------- diagnostics for tracers from convection -------
        if ( id_tracerdt_conv(n) > 0 ) then
          used = send_data ( id_tracerdt_conv(n), tracertnd(:,:,:,n), Time, is, js, 1, &
                             rmask=mask )
        endif

      !------- diagnostics for tracers column integral tendency ------
        if ( id_tracer_conv_col(n) > 0 ) then
          tempdiag(:,:)=0.
          do k=1,kx
            tempdiag(:,:) = tempdiag(:,:) + tracertnd(:,:,k,n)*pmass(:,:,k)
          end do
          used = send_data ( id_tracer_conv_col(n), tempdiag, Time, is, js )
        end if        
      enddo

   end if !end do strat if


! convection diagnostics
  call mpp_clock_end ( convection_clock )
!-----------------------------------------------------------------------
  call mpp_clock_begin ( largescale_clock )

                      if (do_diag_clouds) then
!  capture convective precip and convective spec hum changes (and calculate, 
!  convective spec hum counter) which are needed as predictors 
!  for Gordon's diagnostic clouds
      where (qtnd(:,:,:) < 0.0)
            cnvcntq (:,:,:) = 1.0
      else where
            cnvcntq (:,:,:) = 0.0
      end where
      convprc = precip
                      endif
!-----------------------------------------------------------------------
!***********************************************************************
!----------------- large-scale condensation ----------------------------

                      if (do_lsc) then
!-----------------------------------------------------------------------

   call lscale_cond (tin,qin,pfull,phalf,coldT,rain,snow,ttnd,qtnd,&
                     mask=mask)

!------- (update input values and) compute tendency -----
      tin=tin+ttnd;    qin=qin+qtnd
      
      ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
      rain=rain*dtinv; snow=snow*dtinv

!------- add on tendency ----------
     tdt=tdt+ttnd; qdt=qdt+qtnd

!------- compute rh clouds if desired ------
     if (do_rh_clouds) then
           
           !calculate relative humidity
           call rh_calc(pfull,tin,qin,RH,mask)
           
           !pass RH to rh_clouds_sum
           call rh_clouds_sum (is, js, RH)
           
     end if

!------- compute diagnostic clouds if desired ------
     if (do_diag_clouds) then
           
           !calculate relative humidity
           call rh_calc(pfull,tin,qin,RH,mask)

!  capture tendency of spec. hum. due to lg scale condensation
           lgscldelq = qtnd
           
! pass cloud predictors to diag_cloud_sum
           call diag_cloud_sum (is,js, &
                    tin,qin,rh,omega,lgscldelq,cnvcntq,convprc,kbot)
           
     end if

!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

!-----------------------------------------------------------------------
                           endif
!-----------------------------------------------------------------------
!***********************************************************************
!----------------- stratiform precip/cloud scheme ----------------------
                     if (do_strat) then
     
!-----------------------------------------------------------------------
         do k=1,kx
           mc_full(:,:,k) = 0.5*(mc(:,:,k) + mc(:,:,k+1)) +   &
                                 mc_donner(:,:,k)
         end do
         call strat_driv (Time,is,ie,js,je,dt,pfull,phalf,      &
                          radturbten,                      &
                          tin,qin,tracer(:,:,:,nql),tracer(:,:,:,nqi),&
                          tracer(:,:,:,nqa),omega,mc_full,diff_t,land,&
                          ttnd,qtnd,tracertnd(:,:,:,nql),             &
                          tracertnd(:,:,:,nqi),tracertnd(:,:,:,nqa),  &
                          rain,snow,qrat,ahuco,mask=mask)
     
!
!------- update fields before passing to the clouds
      tracer(:,:,:,nql)=tracer(:,:,:,nql)+tracertnd(:,:,:,nql)
      tracer(:,:,:,nqi)=tracer(:,:,:,nqi)+tracertnd(:,:,:,nqi)
      tracer(:,:,:,nqa)=tracer(:,:,:,nqa)+tracertnd(:,:,:,nqa)
      tin=tin+ttnd; qin=qin+qtnd
      call strat_cloud_sum (is,js,tracer(:,:,:,nql),tracer(:,:,:,nqi),tracer(:,:,:,nqa))
!     call average_clouds (is,js,qlin,qiin,qain)
!     call average_clouds (is,ie,js,je,qlin,qiin,qain)

!------- compute and add on tendencies ----------

      !note mask multiplication is not necessary since the
      !strat scheme already does this
      ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
      tracertnd(:,:,:,nql)=tracertnd(:,:,:,nql)*dtinv
      tracertnd(:,:,:,nqi)=tracertnd(:,:,:,nqi)*dtinv
      tracertnd(:,:,:,nqa)=tracertnd(:,:,:,nqa)*dtinv
      rain=rain*dtinv; snow=snow*dtinv
    
      tdt=tdt+ttnd; qdt=qdt+qtnd
      rdt(:,:,:,nql)=rdt(:,:,:,nql)+tracertnd(:,:,:,nql)
      rdt(:,:,:,nqi)=rdt(:,:,:,nqi)+tracertnd(:,:,:,nqi)
      rdt(:,:,:,nqa)=rdt(:,:,:,nqa)+tracertnd(:,:,:,nqa)

!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

!-----------------------------------------------------------------------
                           endif
! Do wet deposition for the large scale routine here
! qdt should give the vertical distribution here
do n = 1, size(tracertnd,4)
wetdeptnd = 0.0
call wet_deposition(n, t, pfull, phalf, rain, snow, qtnd, tracer(:,:,:,n), &
                    wetdeptnd, Time, 'lscale', is, js, dt)
tracertnd(:,:,:,n) = tracertnd(:,:,:,n) - wetdeptnd
wet_data(:,:,:,n) = wet_data(:,:,:,n) + wetdeptnd

! start chemistry
!if (id_wet(n) /= 0 ) then
!! --------send wet deposition data to diag ----------
!used = send_data(id_wet(n), wet_data(:,:,:,n), Time,is_in=is,js_in=js)
!endif
! end chemistry
enddo

! start chemistry
!cldbot = 0
!cldtop = 0
!if ( get_tracer_index(MODEL_ATMOS,'no') .ne. NO_TRACER ) then
!  do i = 1,size(t,1)
!    do j = 1, size(t,2)
!      do k = 1, size(t,3)
!        if (mc_full(i,j,k) /= 0 ) then
!          cldtop(i,j) = k
!          exit
!        endif
!      enddo
!      do k = size(r,3),1,-1
!        if (mc_full(i,j,k) /= 0 ) then
!          cldbot(i,j) = k
!          exit
!        endif
!      enddo
!    enddo
!  enddo
!  call MOZ_HOOK(cldtop, cldbot, land, zfull, zhalf, t, prod_no, area, lat, &
!                Time, is, js)
!  rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) = &
!    rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) + &
!    prod_no/conc_air
!  if ( id_prod_no > 0 ) then
!    used = send_data(id_prod_no,prod_no, Time, is_in=is, js_in=js)
!  endif
!endif
! end chemistry
!-----------------------------------------------------------------------
!***********************************************************************
!--------------- DIAGNOSTICS FOR LARGE-SCALE SCHEME --------------------
!-----------------------------------------------------------------------

 if ( do_lsc .or. do_strat ) then
!------- diagnostics for dt/dt_strat -------
      if ( id_tdt_ls > 0 ) then
        used = send_data ( id_tdt_ls, ttnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for dq/dt_strat -------
      if ( id_qdt_ls > 0 ) then
        used = send_data ( id_qdt_ls, qtnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for precip_strat -------
      if ( id_prec_ls > 0 ) then
        used = send_data ( id_prec_ls, rain+snow, Time, is, js )
      endif
!------- diagnostics for snow_strat -------
      if ( id_snow_ls > 0 ) then
        used = send_data ( id_snow_ls, snow, Time, is, js )
      endif

!------- diagnostics for water vapor path tendency ----------
      if ( id_q_ls_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + qtnd(:,:,k)*pmass(:,:,k)
        end do
        used = send_data ( id_q_ls_col, tempdiag, Time, is, js )
      end if        
   
!------- diagnostics for dry static energy tendency ---------
      if ( id_t_ls_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + ttnd(:,:,k)*CP_AIR*pmass(:,:,k)
        end do
        used = send_data ( id_t_ls_col, tempdiag, Time, is, js )
      end if        

!------- stratiform cloud tendencies from strat cloud -------
   if ( do_strat ) then

!------- diagnostics for total cumulus mass flux ------
      if ( id_mc_full > 0 ) then
          used = send_data ( id_mc_full, mc_full, Time, is, js, 1, &
                             rmask=mask )
      endif

      !------- diagnostics for dql/dt from strat_cloud ------
      if ( id_qldt_ls > 0 ) then
        used = send_data ( id_qldt_ls, tracertnd(:,:,:,nql), Time, is, js, 1, &
                           rmask=mask )
      endif
      
      !------- diagnostics for dqi/dt from strat_cloud ------
      if ( id_qidt_ls > 0 ) then
        used = send_data ( id_qidt_ls, tracertnd(:,:,:,nqi), Time, is, js, 1, &
                           rmask=mask )
      endif
      
      !------- diagnostics for dqa/dt from strat_cloud ------
      if ( id_qadt_ls > 0 ) then
        used = send_data ( id_qadt_ls, tracertnd(:,:,:,nqa), Time, is, js, 1, &
                           rmask=mask )
      endif

      !------- diagnostics for liquid water path tendency ------
      if ( id_ql_ls_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + tracertnd(:,:,k,nql)*pmass(:,:,k)
        end do
        used = send_data ( id_ql_ls_col, tempdiag, Time, is, js )
      end if        
      
      !------- diagnostics for ice water path tendency ---------
      if ( id_qi_ls_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + tracertnd(:,:,k,nqi)*pmass(:,:,k)
        end do
        used = send_data ( id_qi_ls_col, tempdiag, Time, is, js )
      end if        
      
      !---- diagnostics for stratiform cloud volume tendency ---
      if ( id_qa_ls_col > 0 ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + tracertnd(:,:,k,nqa)*pmass(:,:,k)
        end do
        used = send_data ( id_qa_ls_col, tempdiag, Time, is, js )
      end if        
         
   end if !end do strat if

 endif  !end large-scale or strat diagnostics
  call mpp_clock_end ( largescale_clock )
 
!-----------------------------------------------------------------------
!***********************************************************************
!--------------------- GENERAL DIAGNOSTICS -----------------------------
!-----------------------------------------------------------------------

!------- diagnostics for total precip -------
   if ( id_precip > 0 ) then
        used = send_data ( id_precip, precip, Time, is, js )
   endif


!-----------------------------------------------------------------------
!------- diagnostics for column water vapor, liquid water path and
!------- ice water path


!-- compute and write out water vapor path --
   if ( id_WVP > 0 ) then
        wvp(:,:)=0.
        do k=1,kx
          wvp(:,:) = wvp(:,:) + qin(:,:,k)*pmass(:,:,k)
        end do

        used = send_data ( id_WVP, wvp, Time, is, js )
   end if

!-- compute and write out liquid and ice water path --
   if ( id_LWP > 0 .and. do_strat ) then
        lwp(:,:)=0.
        do k=1,kx
          lwp(:,:) = lwp(:,:) + tracer(:,:,k,nql)*pmass(:,:,k)
        end do

        used = send_data ( id_LWP, lwp, Time, is, js )
   end if

!-- compute and write out liquid and ice water path --
   if ( id_IWP > 0 .and. do_strat ) then
        iwp(:,:)=0.
        do k=1,kx
          iwp(:,:) = iwp(:,:) + tracer(:,:,k,nqi)*pmass(:,:,k)
        end do

        used = send_data ( id_IWP, iwp, Time, is, js )
   end if

!-- compute and write out column integrated cloud mass --
   if ( id_AWP > 0 .and. do_strat ) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + tracer(:,:,k,nqa)*pmass(:,:,k)
        end do

        used = send_data ( id_AWP, tempdiag, Time, is, js )
   end if

!-----------------------------------------------------------------------
!------- diagnostics for relative humidity -------

   if ( id_rh > 0 ) then
      if (.not.(do_rh_clouds.or.do_diag_clouds)) &
          call rh_calc (pfull,tin,qin,RH,mask)
      used = send_data ( id_rh, RH*100., Time, is, js, 1, rmask=mask )
   endif

!-----------------------------------------------------------------------
!---- accumulate global integral of precipiation (mm/day) -----

call sum_diag_integral_field ('prec', precip*86400., is, js)

!-----------------------------------------------------------------------

end subroutine moist_processes

!#######################################################################

subroutine moist_processes_init ( id, jd, kd, lonb, latb, pref, &
                                  axes, Time)

!-----------------------------------------------------------------------
integer,         intent(in) :: id, jd, kd, axes(4)
real, dimension(:), intent(in) :: lonb, latb, pref
type(time_type), intent(in) :: Time
!-----------------------------------------------------------------------
!
!      input
!     --------
!
!      id, jd        number of horizontal grid points in the global
!                    fields along the x and y axis, repectively.
!                      [integer]
!
!      kd            number of vertical points in a column of atmosphere
!-----------------------------------------------------------------------

integer :: unit,io,ierr, n, nt, ntprog
character(len=32) :: tracer_units, tracer_name
!-----------------------------------------------------------------------

       if ( module_is_initialized ) return

       if ( file_exist('input.nml')) then

         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=moist_processes_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'moist_processes_nml')
         enddo
  10     call close_file (unit)

!--------- write version and namelist to standard log ------------

      call write_version_number ( version, tagname )
      if ( mpp_pe() == mpp_root_pe() ) &
      write ( stdlog(), nml=moist_processes_nml )

!------------------- dummy checks --------------------------------------

         if ( do_mca .and. do_ras ) call error_mesg   &
                   ('moist_processes_init',  &
                    'both do_mca and do_ras cannot be specified', FATAL)

         if ( do_lsc .and. do_strat ) call error_mesg   &
                 ('moist_processes_init',  &
                  'both do_lsc and do_strat cannot be specified', FATAL)

         if ( (do_rh_clouds.or.do_diag_clouds) .and. do_strat .and. &
             mpp_pe() == mpp_root_pe() ) call error_mesg ('moist_processes_init', &
     'do_rh_clouds or do_diag_clouds + do_strat should not be specified', NOTE)

         if ( do_rh_clouds .and. do_diag_clouds .and. mpp_pe() == mpp_root_pe() ) &
            call error_mesg ('moist_processes_init',  &
       'do_rh_clouds and do_diag_clouds should not be specified', NOTE)

         if (do_mca .and. do_donner_deep) call error_mesg &
                 ('moist_processes_init',  &
            'both do_donner_deep and do_mca cannot be specified', FATAL)

         if (do_donner_deep .and. do_rh_clouds) then
           call error_mesg ('moist_processes_init',  &
            'Cannot currently activate donner_deep_mod with rh_clouds', FATAL)
         endif   

!        if (do_ras .and. do_donner_deep) call error_mesg &
!                ('moist_processes_init',  &
!           'both do_donner_deep and do_ras cannot be specified', FATAL)

      endif

!---------------------------------------------------------------------
! --- Find the tracer indices 
!---------------------------------------------------------------------

      if (do_strat) then
        ! get tracer indices for stratiform cloud variables
        nsphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
        nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
        nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
        nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
        if (min(nql,nqi,nqa) <= 0) call error_mesg ('moist_processes', &
                                                    'stratiform cloud tracer(s) not found', FATAL)
        if (nql == nqi .or. nqa == nqi .or. nql == nqa) call error_mesg ('moist_processes',  &
                                 'tracers indices cannot be the same (i.e., nql=nqi=nqa).', FATAL)
        if (mpp_pe() == mpp_root_pe()) &
            write (stdlog(),'(a,3i4)') 'Stratiform cloud tracer indices: nql,nqi,nqa =',nql,nqi,nqa
      endif

!------------ initialize various schemes ----------
      if (do_mca .or. do_donner_deep)    call  moist_conv_init (axes,Time)
      if (do_lsc) then
                     call lscale_cond_init ()
                     if (do_rh_clouds) call rh_clouds_init (id,jd,kd)
                     if (do_diag_clouds) call diag_cloud_init (id,jd,kd,ierr)
      endif
      if (do_ras)    call         ras_init (do_strat, axes,Time)
      if (do_strat)  call strat_cloud_init (axes,Time,id,jd,kd)
      if (do_dryadj) call     dry_adj_init ()
      if (do_donner_deep) call donner_deep_init (lonb, latb, pref, &
                                                 axes, Time)
      if (do_cmt)    call cu_mo_trans_init (axes,Time)

      call get_number_tracers( MODEL_ATMOS, num_tracers=nt, num_prog = ntprog)
      do n = 1,ntprog
        call get_tracer_names( MODEL_ATMOS, n, name = tracer_name, units = tracer_units)
      
         id_tracerdt_conv(n) = register_diag_field ( mod_name, &
           TRIM(tracer_name)//'dt_conv', axes(1:3), Time, &
           'Tracer tendency from moist convection ',  TRIM(tracer_units)//'/s',  &
                              missing_value=missing_value               )
      
         id_tracer_conv_col(n) = register_diag_field ( mod_name, &
           TRIM(tracer_name)//'dt_conv_col', axes(1:2), Time, &
          'Column integrated tracer tendency from moist convection', &
           TRIM(tracer_units)//'/s',  &
           missing_value=missing_value               )
      enddo   
  
!----- initialize quantities for global integral package -----

   call diag_integral_field_init ('prec', 'f6.3')
   

!----- initialize quantities for diagnostics output -----

   call diag_field_init ( axes, Time )

!----- initialize clocks -----

   convection_clock =     &
       mpp_clock_id( '   Physics_up: Moist Proc: Conv', &
           grain=CLOCK_MODULE, flags = MPP_CLOCK_SYNC )
   largescale_clock =     &
       mpp_clock_id( '   Physics_up: Moist Proc: LS', &
           grain=CLOCK_MODULE, flags = MPP_CLOCK_SYNC )

   module_is_initialized = .true.

!-----------------------------------------------------------------------

end subroutine moist_processes_init

!#######################################################################

subroutine moist_processes_end

      if( .not.module_is_initialized ) return

!----------------close various schemes-----------------

      if (do_strat)       call strat_cloud_end
      if (do_rh_clouds)   call   rh_clouds_end
      if (do_diag_clouds) call  diag_cloud_end
      if (do_donner_deep) call donner_deep_end
      if (do_cmt        ) call cu_mo_trans_end
      if (do_ras        ) call         ras_end
      module_is_initialized = .false.
      
!-----------------------------------------------------------------------

end subroutine moist_processes_end

!#######################################################################
function doing_strat()
logical :: doing_strat

  if (.not. module_is_initialized) call error_mesg ('doing_strat',  &
                     'moist_processes_init has not been called.', FATAL)

  doing_strat = do_strat

end function doing_strat
!#######################################################################

      subroutine tempavg (pdepth,phalf,temp,tsnow,mask)

!-----------------------------------------------------------------------
!
!    computes a mean atmospheric temperature for the bottom
!    "pdepth" pascals of the atmosphere.
!
!   input:  pdepth     atmospheric layer in pa.
!           phalf      pressure at model layer interfaces
!           temp       temperature at model layers
!           mask       data mask at model layers (0.0 or 1.0)
!
!   output:  tsnow     mean model temperature in the lowest
!                      "pdepth" pascals of the atmosphere
!
!-----------------------------------------------------------------------
      real, intent(in)  :: pdepth
      real, intent(in) , dimension(:,:,:) :: phalf,temp
      real, intent(out), dimension(:,:)   :: tsnow
      real, intent(in) , dimension(:,:,:), optional :: mask
!-----------------------------------------------------------------------
 real, dimension(size(temp,1),size(temp,2)) :: prsum, done, pdel, pdep
 real  sumdone
 integer  k
!-----------------------------------------------------------------------

      tsnow=0.0; prsum=0.0; done=1.0; pdep=pdepth

      do k=size(temp,3),1,-1

         if (present(mask)) then
           pdel(:,:)=(phalf(:,:,k+1)-phalf(:,:,k))*mask(:,:,k)*done(:,:)
         else
           pdel(:,:)=(phalf(:,:,k+1)-phalf(:,:,k))*done(:,:)
         endif

         where ((prsum(:,:)+pdel(:,:))  >  pdep(:,:))
            pdel(:,:)=pdepth-prsum(:,:)
            done(:,:)=0.0
            pdep(:,:)=0.0
         endwhere

         tsnow(:,:)=tsnow(:,:)+pdel(:,:)*temp(:,:,k)
         prsum(:,:)=prsum(:,:)+pdel(:,:)

         sumdone=sum(done(:,:))
         if (sumdone < 1.e-4) exit

      enddo

         tsnow(:,:)=tsnow(:,:)/prsum(:,:)

!-----------------------------------------------------------------------

      end subroutine tempavg

!#######################################################################

      subroutine rh_calc(pfull,T,qv,RH,MASK)

        IMPLICIT NONE


        REAL, INTENT (IN),    DIMENSION(:,:,:) :: pfull,T,qv
        REAL, INTENT (OUT),   DIMENSION(:,:,:) :: RH
        REAL, INTENT (IN), OPTIONAL, DIMENSION(:,:,:) :: MASK

        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: esat
        
!-----------------------------------------------------------------------
!       Calculate RELATIVE humidity.
!       This is calculated according to the formula:
!
!       RH   = qv / (epsilon*esat/ [pfull  -  (1.-epsilon)*esat])
!
!       Where epsilon = RDGAS/RVGAS = d622
!
!       and where 1- epsilon = d378
!
!       Note that RH does not have its proper value
!       until all of the following code has been executed.  That
!       is, RH is used to store intermediary results
!       in forming the full solution.

        !calculate water saturated vapor pressure from table
        !and store temporarily in the variable esat
        CALL LOOKUP_ES(T,esat)
        
        !calculate denominator in qsat formula
        RH(:,:,:) = pfull(:,:,:)-d378*esat(:,:,:)
     
        !limit denominator to esat, and thus qs to epsilon
        !this is done to avoid blow up in the upper stratosphere
        !where pfull ~ esat
        RH(:,:,:) = MAX(RH(:,:,:),esat(:,:,:)) 
        
        !calculate RH
        RH(:,:,:)=qv(:,:,:)/(d622*esat(:,:,:)/RH(:,:,:))
      
        !IF MASK is present set RH to zero
        IF (present(MASK)) RH(:,:,:)=MASK(:,:,:)*RH(:,:,:)


END SUBROUTINE rh_calc



!#######################################################################

subroutine diag_field_init ( axes, Time )

  integer,         intent(in) :: axes(4)
  type(time_type), intent(in) :: Time

  integer, dimension(3) :: half = (/1,2,4/)

!------------ initializes diagnostic fields in this module -------------

!if ( do_donner_deep .and. do_ras ) then

   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency from donner_deep plus RAS ',    'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Spec humidity tendency from donner_deep plus RAS ',  'kg/kg/s',  &
                        missing_value=missing_value               )

   id_q_conv_col = register_diag_field ( mod_name, &
     'q_conv_col', axes(1:2), Time, &
    'Water vapor path tendency from donner_deep plus RAS ',   'kg/m2/s' )
   
   id_t_conv_col = register_diag_field ( mod_name, &
     't_conv_col', axes(1:2), Time, &
    'Column static energy tendency from donner_deep plus RAS ','W/m2' )
   
   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate from donner_deep plus RAS ',       'kg/m2/s' )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate from donner_deep plus RAS ',       'kg/m2/s' )

   id_gust_conv = register_diag_field ( mod_name, &
     'gust_conv', axes(1:2), Time, &
    'Gustiness from deep convection ',       'm/s' )

   


!end if


if (( do_donner_deep ) .and. (.not. do_ras)) then

   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency from donner_deep   ',    'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Spec humidity tendency from donner_deep   ',  'kg/kg/s',  &
                        missing_value=missing_value               )

   id_q_conv_col = register_diag_field ( mod_name, &
     'q_conv_col', axes(1:2), Time, &
    'Water vapor path tendency from donner_deep',   'kg/m2/s' )
   
   id_t_conv_col = register_diag_field ( mod_name, &
     't_conv_col', axes(1:2), Time, &
    'Column static energy tendency from donner_deep','W/m2' )
   
   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate from donner_deep   ',       'kg/m2/s' )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate from donner_deep   ',       'kg/m2/s' )

end if

if ( do_donner_deep .and. do_strat ) then

   id_qldt_conv = register_diag_field ( mod_name, &
     'qldt_conv', axes(1:3), Time, &
     'Liquid water tendency from donner_deep',      'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qidt_conv = register_diag_field ( mod_name, &
     'qidt_conv', axes(1:3), Time, &
     'Ice water tendency from donner_deep',         'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qadt_conv = register_diag_field ( mod_name, &
     'qadt_conv', axes(1:3), Time, &
     'Cloud fraction tendency from donner_deep',    '1/sec',    &
                        missing_value=missing_value               )

   id_ql_conv_col = register_diag_field ( mod_name, &
     'ql_conv_col', axes(1:2), Time, &
    'Liquid water path tendency from donner_deep',  'kg/m2/s' )
   
   id_qi_conv_col = register_diag_field ( mod_name, &
     'qi_conv_col', axes(1:2), Time, &
    'Ice water path tendency from donner_deep',     'kg/m2/s' )
   
   id_qa_conv_col = register_diag_field ( mod_name, &
     'qa_conv_col', axes(1:2), Time, &
    'Cloud mass tendency from donner_deep',         'kg/m2/s' )
      
endif

if ( do_lsc ) then

   id_tdt_ls = register_diag_field ( mod_name, &
     'tdt_ls', axes(1:3), Time, &
       'Temperature tendency from large-scale cond',   'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_ls = register_diag_field ( mod_name, &
     'qdt_ls', axes(1:3), Time, &
     'Spec humidity tendency from large-scale cond', 'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_ls = register_diag_field ( mod_name, &
     'prec_ls', axes(1:2), Time, &
    'Precipitation rate from large-scale cond',     'kg/m2/s' )

   id_snow_ls = register_diag_field ( mod_name, &
     'snow_ls', axes(1:2), Time, &
    'Frozen precip rate from large-scale cond',     'kg/m2/s' )

   id_q_ls_col = register_diag_field ( mod_name, &
     'q_ls_col', axes(1:2), Time, &
    'Water vapor path tendency from large-scale cond','kg/m2/s' )
   
   id_t_ls_col = register_diag_field ( mod_name, &
     't_ls_col', axes(1:2), Time, &
    'Column static energy tendency from large-scale cond','W/m2' )
   
endif

if ( do_strat ) then

   id_mc_full = register_diag_field ( mod_name, &
     'mc_full', axes(1:3), Time, &
     'Net Mass Flux from donner deep plus RAS',   'kg/m2/s', &
                       missing_value=missing_value               )

   id_tdt_ls = register_diag_field ( mod_name, &
     'tdt_ls', axes(1:3), Time, &
     'Temperature tendency from strat cloud',        'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_ls = register_diag_field ( mod_name, &
     'qdt_ls', axes(1:3), Time, &
     'Spec humidity tendency from strat cloud',      'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_ls = register_diag_field ( mod_name, &
     'prec_ls', axes(1:2), Time, &
    'Precipitation rate from strat cloud',          'kg/m2/s' )

   id_snow_ls = register_diag_field ( mod_name, &
     'snow_ls', axes(1:2), Time, &
    'Frozen precip rate from strat cloud',          'kg/m2/s' )

   id_q_ls_col = register_diag_field ( mod_name, &
     'q_ls_col', axes(1:2), Time, &
    'Water vapor path tendency from strat cloud',   'kg/m2/s' )
   
   id_t_ls_col = register_diag_field ( mod_name, &
     't_ls_col', axes(1:2), Time, &
    'Column static energy tendency from strat cloud','W/m2' )
   
   id_qldt_ls = register_diag_field ( mod_name, &
     'qldt_ls', axes(1:3), Time, &
     'Liquid water tendency from strat cloud',       'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qidt_ls = register_diag_field ( mod_name, &
     'qidt_ls', axes(1:3), Time, &
     'Ice water tendency from strat cloud',          'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qadt_ls = register_diag_field ( mod_name, &
     'qadt_ls', axes(1:3), Time, &
     'Cloud fraction tendency from strat cloud',     '1/sec',    &
                        missing_value=missing_value               )

   id_ql_ls_col = register_diag_field ( mod_name, &
     'ql_ls_col', axes(1:2), Time, &
    'Liquid water path tendency from strat cloud',   'kg/m2/s' )
   
   id_qi_ls_col = register_diag_field ( mod_name, &
     'qi_ls_col', axes(1:2), Time, &
    'Ice water path tendency from strat cloud',      'kg/m2/s' )
   
   id_qa_ls_col = register_diag_field ( mod_name, &
     'qa_ls_col', axes(1:2), Time, &
    'Cloud mass tendency from strat cloud',          'kg/m2/s' )
      
endif

   id_precip = register_diag_field ( mod_name, &
     'precip', axes(1:2), Time, &
     'Total precipitation rate',                     'kg/m2/s' )

   id_WVP = register_diag_field ( mod_name, &
     'WVP', axes(1:2), Time, &
        'Column integrated water vapor',                'kg/m2'  )

if ( do_strat ) then

   id_LWP = register_diag_field ( mod_name, &
     'LWP', axes(1:2), Time, &
        'Liquid water path',                            'kg/m2'   )

   id_IWP = register_diag_field ( mod_name, &
     'IWP', axes(1:2), Time, &
        'Ice water path',                               'kg/m2'   )

   id_AWP = register_diag_field ( mod_name, &
     'AWP', axes(1:2), Time, &
        'Column integrated cloud mass ',                'kg/m2'   )

endif

   id_tdt_dadj = register_diag_field ( mod_name, &
     'tdt_dadj', axes(1:3), Time, &
   'Temperature tendency from dry conv adj',       'deg_K/s',  &
                        missing_value=missing_value               )

   id_rh = register_diag_field ( mod_name, &
     'rh', axes(1:3), Time, &
         'relative humidity',                            'percent',  & 
                        missing_value=missing_value               )

   
if (do_donner_deep) then

   id_tdt_deep_donner= register_diag_field ( mod_name, &
           'tdt_deep_donner', axes(1:3), Time, &
           ' heating rate - deep portion', 'deg K/s', &
                        missing_value=missing_value               )

   id_tdt_mca_donner = register_diag_field ( mod_name, &
           'tdt_mca_donner', axes(1:3), Time, &
           ' heating rate - mca  portion', 'deg K/s', &
                        missing_value=missing_value               )

   id_qdt_deep_donner = register_diag_field ( mod_name, &
           'qdt_deep_donner', axes(1:3), Time, &
           ' moistening rate - deep portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qdt_mca_donner = register_diag_field ( mod_name, &
           'qdt_mca_donner', axes(1:3), Time, &
           ' moistening rate - mca  portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_prec_deep_donner = register_diag_field ( mod_name, &
           'prc_deep_donner', axes(1:2), Time, &
           ' total precip rate - deep portion', 'kg/m2/s', &
                        missing_value=missing_value               )

   id_prec_mca_donner = register_diag_field ( mod_name, &
           'prc_mca_donner', axes(1:2), Time, &
           ' total precip rate - mca  portion', 'kg/m2/s', &
                        missing_value=missing_value               )

   id_snow_deep_donner = register_diag_field ( mod_name, &
           'snow_deep_donner', axes(1:2), Time, &
           ' frozen precip rate - deep portion', 'kg/m2/s', &
                        missing_value=missing_value               )

   id_snow_mca_donner = register_diag_field ( mod_name, &
           'snow_mca_donner', axes(1:2), Time, &
           ' frozen precip rate -  mca portion', 'kg/m2/s', &
                        missing_value=missing_value               )

   id_mc_donner = register_diag_field ( mod_name, &
           'mc_donner', axes(1:3), Time, &
           'Net Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

! start chemistry
   if (get_tracer_index(MODEL_ATMOS,'no') > 0) &
     id_prod_no = register_diag_field ( 'tracers', &
             'hook_no', axes(1:3), Time, &
             'hook_no',   '1/cm3/s')
!end chemistry
endif
!-----------------------------------------------------------------------

end subroutine diag_field_init

!#######################################################################


                 end module moist_processes_mod
