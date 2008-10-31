 
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
!             diagnostic cloud scheme 
!             lin cloud microphysics
!             betts-miller convective adjustment
!
!-----------------------------------------------------------------------

use    betts_miller_mod, only: betts_miller, betts_miller_init
use     bm_massflux_mod, only: bm_massflux, bm_massflux_init
use          bm_omp_mod, only: bm_omp, bm_omp_init

use     donner_deep_mod, only: donner_deep_init,               &
                               donner_deep, donner_deep_end,   &
                               donner_deep_restart
use      moist_conv_mod, only: moist_conv, moist_conv_init
use     lscale_cond_mod, only: lscale_cond, lscale_cond_init
use  sat_vapor_pres_mod, only: compute_qs, lookup_es

use         uw_conv_mod, only: uw_conv, uw_conv_end, uw_conv_init

use lin_cld_microphys_mod, only:  lin_cld_microphys_driver, lin_cld_microphys_init, &
                                  lin_cld_microphys_end

use    time_manager_mod, only: time_type, get_time

use    diag_manager_mod, only: register_diag_field, send_data

use             fms_mod, only: file_exist, check_nml_error,    &
                               open_namelist_file, close_file, &
                               write_version_number,           &
                               mpp_pe, mpp_root_pe, stdlog,    &
                               error_mesg, FATAL, NOTE,        &
                               mpp_clock_id, mpp_clock_begin,  &
                               mpp_clock_end, CLOCK_MODULE,    &
                               CLOCK_LOOP,  MPP_CLOCK_SYNC,    &
                               read_data, write_data

use             ras_mod, only: ras, ras_end, ras_init

use         dry_adj_mod, only: dry_adj, dry_adj_init

use     strat_cloud_mod, only: strat_cloud_init, strat_cloud, strat_cloud_end, &
                               strat_cloud_sum, strat_cloud_restart

use       rh_clouds_mod, only: rh_clouds_init, rh_clouds_end, &
                               rh_clouds_sum

use      diag_cloud_mod, only: diag_cloud_init, diag_cloud_end, &
                               diag_cloud_sum, diag_cloud_restart

use   diag_integral_mod, only: diag_integral_field_init, &
                               sum_diag_integral_field

use       constants_mod, only: CP_AIR, GRAV, HLV, HLS, HLF, &
                               RDGAS, RVGAS, TFREEZE, &
                               SECONDS_PER_DAY, KAPPA

use     cu_mo_trans_mod, only: cu_mo_trans_init, cu_mo_trans, cu_mo_trans_end

use   field_manager_mod, only: MODEL_ATMOS
use  tracer_manager_mod, only: get_tracer_index,&
                               get_number_tracers, &
                               get_tracer_names, &
                               query_method, &
                               NO_TRACER
use atmos_tracer_utilities_mod, only : wet_deposition

use       moz_hook_mod, only : moz_hook

use rad_utilities_mod,       only: aerosol_type

implicit none
private

  integer :: excl_pre_conv, excl_donner_deep, excl_donner_sph, excl_donner_moist_cons, &
             excl_donner_deep_2, excl_moist_conv_adj, excl_ras, excl_cmt, excl_tracer_dens, &
             excl_uw_conv, excl_shallow, excl_conv_diag, excl_lsc, excl_strat, excl_wet_dep, &
             excl_lsc_diag, excl_gen_diag

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

   public   moist_processes, moist_processes_init, moist_processes_end, &
            doing_strat, moist_processes_restart

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

   real,parameter :: EPST=200.

   integer :: nsphum, nql, nqi, nqa, nqn  ! tracer indices for stratiform clouds
   integer :: nqr, nqs, nqg               ! additional tracer indices for Lin Micro-Physics
   integer :: ktop                        ! top layer index for Lin Micro-Physics
!--------------------- version number ----------------------------------
   character(len=128) :: &
   version = '$Id: moist_processes.F90,v 16.0.6.4 2008/09/18 19:43:14 rab Exp $'
   character(len=128) :: tagname = '$Name: perth_2008_10 $'
   logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------
!-------------------- namelist data (private) --------------------------


   logical :: do_mca=.true., do_lsc=.true., do_ras=.false.,  &
              do_strat=.false., do_dryadj=.false., do_uw_conv=.false., &
              do_rh_clouds=.false., do_diag_clouds=.false., &
              do_donner_deep=.false., do_cmt=.false., &
              use_tau=.false., do_gust_cv = .false., &
              do_lin_cld_microphys=.false., do_liq_num = .false., &
              do_donner_mca=.true., do_bm=.false., do_bmmass =.false., & 
              do_bmomp  =.false., do_simple =.false. 

   logical :: force_donner_moist_conserv = .false.
   logical :: do_unified_convective_closure = .false.
   logical :: do_limit_donner = .false. ! .false. produces previous 
                                        ! behavior (cjg)
   logical :: do_limit_uw = .false.     ! .false. produces previous
                                        ! behavior (cjg )
   logical :: do_donner_conservation_checks = .false.
   logical :: using_fms = .true.

   integer :: tau_sg = 0
   integer :: k_sg = 2

   character(len=64)  :: cmt_mass_flux_source = 'ras'
   real :: pdepth = 150.e2
   real :: gustmax = 3.             ! maximum gustiness wind (m/s)
   real :: gustconst = 10./SECONDS_PER_DAY   ! constant in kg/m2/sec, default =
                                    ! 1 cm/day = 10 mm/day
                                    
!---------------- namelist variable definitions ------------------------
!
!   do_limit_donner = limit Donner deeo tendencies to prevent the
!               formation of grid points with negative water vapor,
!               liquid or ice.
!
!   do_limit_uw = limit UW shallow tendencies to prevent the formation
!               of grid points with negative total water specific 
!               humidities. This situation can occur because both
!               shallow and deep convection operate on the same
!               soundings without knowledge of what the other is doing
!
!   do_unified_convective_closure = use cloud base mass flux calculated
!               in uw_conv module as value for donner deep parameter-
!               ization; adjust cbmf available for uw shallow appropr-
!               iately. only available when uw shallow and donner deep
!               are the active convective schemes
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
!  do_lin_cld_microphys = switch to turn on/off the Lin Cloud Micro-Physics scheme
!                [logical, default: do_lin_cld_microphys=false ]
!  do_liq_num = switch to turn on/off the prognostic droplet number scheme.
!                [logical, default: do_liq_num=false ]
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
!  cmt_mass_flux_source = parameterization(s) being used to supply the 
!               mass flux profiles seen by the cumulus momentum transport
!               module; currently either 'ras', 'donner', 'uw', 
!               'donner_and_ras', 'donner_and_uw', 'ras_and_uw', 
!               'donner_and_ras_and_uw' or 'all'
!
!   do_bm    = switch to turn on/off betts-miller scheme
!                [logical, default: do_bm=false ]
!   do_bmmass  = switch to turn on/off betts-miller massflux scheme
!                [logical, default: do_bmmass=false ]
!   do_bmomp  = switch to turn on/off olivier's version of the betts-miller 
!               scheme (with separated boundary layer)
!                [logical, default: do_bmomp=false ]
!   do_simple = switch to turn on alternative definition of specific humidity.
!               When true, specific humidity = (rdgas/rvgas)*esat/pressure
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
   logical :: do_donner_before_uw = .false.
   logical :: use_updated_profiles_for_uw = .false.
   logical :: only_one_conv_scheme_per_column = .false.
   logical :: limit_conv_cloud_frac = .false.

namelist /moist_processes_nml/ do_mca, do_lsc, do_ras, do_uw_conv, do_strat,  &
                               do_donner_before_uw, &
                               use_updated_profiles_for_uw, &
                               only_one_conv_scheme_per_column, &
                               limit_conv_cloud_frac, &
                               do_unified_convective_closure, &
                               do_dryadj, pdepth, do_lin_cld_microphys, tau_sg, k_sg,    &
                               cmt_mass_flux_source, &
                               use_tau, do_rh_clouds, do_diag_clouds, &
                               do_donner_deep, do_cmt, do_gust_cv, &
                               gustmax, gustconst, do_liq_num, &
                               force_donner_moist_conserv, &
                               do_donner_conservation_checks, &
                               do_donner_mca, do_limit_uw, &
                               do_limit_donner, using_fms, &
                               do_bm, do_bmmass, do_bmomp, do_simple

!-----------------------------------------------------------------------
!-------------------- diagnostics fields -------------------------------

integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_tdt_ls  , id_qdt_ls  , id_prec_ls  , id_snow_ls  , &
           id_precip  , id_WVP, id_LWP, id_IWP, id_AWP, id_gust_conv, &
           id_tdt_dadj, id_rh,  id_qs, id_mc, id_mc_donner, id_mc_full, &
           id_tdt_deep_donner, id_qdt_deep_donner, &
           id_qadt_deep_donner, id_qldt_deep_donner, &
           id_qidt_deep_donner, &
           id_tdt_mca_donner, id_qdt_mca_donner, &
           id_prec_deep_donner, id_prec_mca_donner,&
           id_tdt_uw, id_qdt_uw, &
           id_qadt_uw, id_qldt_uw, id_qidt_uw, id_qndt_uw, &
           id_prec1_deep_donner, &
           id_snow_deep_donner, id_snow_mca_donner, &
           id_qadt_ls, id_qldt_ls, id_qndt_ls, id_qidt_ls, &
           id_qadt_conv, id_qldt_conv, id_qndt_conv, id_qidt_conv, &
           id_qa_ls_col, id_ql_ls_col, id_qn_ls_col, id_qi_ls_col, &
           id_qa_conv_col, id_ql_conv_col, id_qn_conv_col,  &
           id_qi_conv_col, &
           id_bmflag, id_klzbs, id_invtaubmt, id_invtaubmq, &
           id_massflux, id_entrop_ls, &
           id_cape, id_cin, id_tref, id_qref, &
           id_q_conv_col, id_q_ls_col, id_t_conv_col, id_t_ls_col, &

           id_enth_moist_col, id_wat_moist_col, &
           id_enth_ls_col, id_wat_ls_col, &
           id_enth_conv_col, id_wat_conv_col, &
           id_enth_donner_col, id_wat_donner_col, &
           id_enth_donner_col2,  &
           id_enth_donner_col3,  &
           id_enth_donner_col4,  &
           id_enth_donner_col5,  &
           id_enth_donner_col6,  &
           id_enth_donner_col7,  &
           id_enth_mca_donner_col, id_wat_mca_donner_col, &
           id_enth_uw_col, id_wat_uw_col, &
           id_scale_donner, id_scale_uw, &
           id_ras_precip, id_ras_freq, id_don_precip, id_don_freq, &
           id_lsc_precip, id_lsc_freq, id_uw_precip, id_uw_freq, &
           id_prod_no, id_m_cdet_donner, id_m_cellup, &
           id_conv_rain3d, id_conv_snow3d, id_lscale_rain3d, id_lscale_snow3d
 
integer :: id_qvout, id_qaout, id_qlout, id_qiout

integer    :: id_vaporint, id_condensint, id_precipint, id_diffint
integer    :: id_vertmotion
integer    :: id_max_enthalpy_imbal_don, id_max_water_imbal_don
integer    :: id_max_enthalpy_imbal, id_max_water_imbal
integer    :: id_enthint, id_lprcp, id_lcondensint, id_enthdiffint

integer, dimension(:), allocatable :: id_tracerdt_conv,  &
                                      id_tracerdt_conv_col, &
                                      id_conv_tracer,  &
                                      id_conv_tracer_col, &
                                      id_tracerdt_mcadon, &
                                      id_tracerdt_mcadon_col, &
                                      id_wet_deposition
character(len=5) :: mod_name = 'moist'

real :: missing_value = -999.
integer :: convection_clock, largescale_clock, &
           donner_clock, mca_clock, ras_clock, cmt_clock, &
           closure_clock, &
           lscalecond_clock, stratcloud_clock, shallowcu_clock

logical :: do_tracers_in_donner =.false.
logical :: do_tracers_in_mca = .false.
logical :: do_tracers_in_ras = .false.
logical :: do_tracers_in_uw = .false.
logical, dimension(:), allocatable :: tracers_in_donner,   &
                                      tracers_in_mca, tracers_in_ras, &
                                      tracers_in_uw
integer :: num_donner_tracers=0
integer :: num_mca_tracers=0
integer :: num_ras_tracers=0
integer :: num_uw_tracers=0
integer :: num_tracers=0

logical :: cmt_uses_donner = .false.
logical :: cmt_uses_ras = .false.
logical :: cmt_uses_uw  = .false.

logical :: doing_diffusive

real, dimension(:,:), allocatable :: max_enthalpy_imbal, max_water_imbal
real, dimension(:,:), allocatable :: max_enthalpy_imbal_don,    &
                                      max_water_imbal_don
!-----------------------------------------------------------------------




                             contains




!#######################################################################

subroutine moist_processes (is, ie, js, je, Time, dt, land,            &
                            phalf, pfull, zhalf, zfull, omega, diff_t, &
                            radturbten, cush, cbmf,                    &
                            pblht, ustar, bstar, qstar,                &
                            t, q, r, u, v, tm, qm, rm, um, vm,         &
                            tdt, qdt, rdt, udt, vdt, diff_cu_mo,       &
                            convect, lprec, fprec, gust_cv, area,      &
                            lat,   &
                            lsc_cloud_area, lsc_liquid, lsc_ice, &
                            lsc_droplet_number, &
                            Aerosol, mask, kbot, &
                            shallow_cloud_area, shallow_liquid,  &
                            shallow_ice, shallow_droplet_number, &
                            cell_cld_frac, cell_liq_amt, cell_liq_size, &
                            cell_ice_amt, cell_ice_size, &
                            cell_droplet_number, &
                            meso_cld_frac, meso_liq_amt, meso_liq_size, &
                            meso_ice_amt, meso_ice_size,  &
                            meso_droplet_number, nsum_out)

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
   real, intent(in) , dimension(:,:)     :: land, pblht, ustar, bstar, qstar
   real, intent(inout), dimension(:,:)   :: cush, cbmf
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

   real, intent(out) , dimension(:,:,:) ::  lsc_cloud_area, lsc_liquid,&
                                            lsc_ice, lsc_droplet_number

   type(aerosol_type),intent(in),       optional :: Aerosol
   real, intent(in) , dimension(:,:,:), optional :: mask
integer, intent(in) , dimension(:,:),   optional :: kbot

   real, intent(inout) , dimension(:,:,:), optional :: &      
                          shallow_cloud_area, shallow_liquid,   &
                          shallow_ice, shallow_droplet_number, &
                          cell_cld_frac, cell_liq_amt, cell_liq_size, &
                          cell_ice_amt, cell_ice_size, &
                          cell_droplet_number, &
                          meso_cld_frac, meso_liq_amt, meso_liq_size, &
                          meso_ice_amt, meso_ice_size, &
                          meso_droplet_number
integer, intent(inout) , dimension(:,:), optional ::  nsum_out

!-----------------------------------------------------------------------
real, dimension(size(t,1),size(t,2),size(t,3)) :: tin,qin,ttnd,qtnd, &
                                                  cf, delta_temp,&
                                                  delta_vapor,  &
                                                  delta_q, rin, rtnd,   &
                                                  donner_humidity_area, &
                                                donner_humidity_factor,&
                                         convective_humidity_area,&
                                          convective_humidity_ratio
real, dimension(size(t,1),size(t,2),size(t,3)) :: tin_orig,qin_orig
real, dimension(size(t,1),size(t,2),size(t,3)) :: ttnd_uw,qtnd_uw, qltnd_uw,&
                                                  qitnd_uw, qatnd_uw, &
                                                  qntnd_uw, utnd_uw, vtnd_uw
           
real, dimension(size(t,1),size(t,2)) ::  cbmf_clo        
real, dimension(size(t,1),size(t,2),size(phalf,3)) :: conv_rain3d, conv_snow3d
real, dimension(size(t,1),size(t,2),size(phalf,3)) :: lscale_rain3d, lscale_snow3d
real, dimension(size(t,1),size(t,2),size(t,3)) :: qtnd_wet, & ! specific humidity tendency (kg/kg/s)
                                                  cloud_wet, &! cloud liquid+ice (kg/kg)
                                                  cloud_frac  ! cloud area fraction
real, dimension(size(t,1),size(t,3)) :: dp, temp
real, dimension(size(t,1),size(t,2),size(t,3)) :: ttnd_conv,qtnd_conv
real, dimension(size(t,1),size(t,2))           :: tsnow,snow
logical,dimension(size(t,1),size(t,2))         :: coldT
logical,dimension(size(t,1),size(t,2))         :: tmplmask
real, dimension(size(t,1),size(t,2))           :: freq_count
real, dimension(size(t,1),size(t,2))           :: lheat_precip, vert_motion, &
                                                  total_precip
real, dimension(size(t,1),size(t,2),size(t,3)) :: liquid_precip, & 
                                                  frozen_precip
real, dimension(size(t,1),size(t,2),size(t,3)) :: utnd,vtnd,uin,vin
real, dimension(size(t,1),size(t,2),size(t,3)) :: qsat

real, dimension(size(t,1),size(t,2),size(t,3)) :: qltnd,qitnd,qatnd, qntnd,&
                                                  delta_ql, delta_qi, delta_qa, &
                                                  qlin, qiin, qain
real, dimension(size(t,1),size(t,2))           :: tkemiz !miz
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: tracer, rdt_init , qtrcumo
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: tracer_orig
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: mc,mask3, &
                                                    m_cellup, mc_cmt
real, dimension(size(t,1),size(t,2),size(t,3)) :: det0, det_cmt       
real, dimension(size(t,1),size(t,2),size(t,3)) :: tdt_init, qdt_init
real, dimension(size(t,1),size(t,2),size(t,3)) :: mc_full, mc_donner, &
                                                  m_cdet_donner, massflux
real, dimension(size(t,1),size(t,2),size(t,3)) :: RH, pmass, wetdeptnd, q_ref, t_ref
real, dimension(size(t,1),size(t,2))           :: rain, precip,  &
                                                  precip_returned, &
                                                  precip_adjustment, cape, cin
real, dimension(size(t,1),size(t,2))           :: rain_uw, snow_uw, &
                                                  rain_don, snow_don, &
                                                  rain_ras, snow_ras, &
                                                  rain_donmca, &
                                                  snow_donmca
real, dimension(size(t,1),size(t,2),size(t,3)) :: ttnd_don,qtnd_don, &
                                                  tin_pass, qin_pass, &
                                                  ttnd_donmca, &
                                                  qtnd_donmca 
real, dimension(size(t,1),size(t,2))           :: wvp,lwp,iwp
real, dimension(size(t,1),size(t,2))           :: bmflag, & 
                                                  klzbs, invtaubmt, invtaubmq


!     sfc_sh_flux      sensible heat flux across the surface
!                      [ watts / m**2 ]
!     sfc_vapor_flux   water vapor flux across the surface
!                      [ kg(h2o) / (m**2 sec) ]
real, dimension(size(t,1),size(t,2))           :: sfc_sh_flux, sfc_vapor_flux


!
! additional Lin Micro-Physics properties
real, dimension(size(t,1),size(t,2))           :: ice_lin, graupel_lin
real, dimension(size(t,1),size(t,2),size(t,3)) :: delp, delz
real, dimension(size(t,1),size(t,2),size(t,3)) :: qrtnd, qstnd, qgtnd


real, dimension(size(t,1),size(t,2))           :: tempdiag
real, dimension(size(t,1),size(t,2))           :: tempdiag4
real, dimension(size(t,1),size(t,2))           :: tempdiag5
real, dimension(size(t,1),size(t,2))           :: tempdiag6
real, dimension(size(t,1),size(t,2),size(t,3)) :: tempdiag1
integer n

real, dimension(size(t,1),size(t,2))           :: scale_donner, scale_uw
! for uw_conv
real, dimension(size(t,1),size(t,2),size(t,3)) :: cmf,thlflx,qtflx,precflx

integer :: i, j, k, ix, jx, kx, nt, ip, tr
real    :: dtinv
logical ::           do_adjust, used, avgbl

! tracer code:
integer :: nn
real, dimension (size(t,1), size(t,2), &
                 size(t,3), num_donner_tracers) :: qtrtnd , &
                                                   donner_tracers
real, dimension (size(t,1), size(t,2), &
                 size(t,3), num_ras_tracers) :: qtrras , &
                                                    ras_tracers
real, dimension (size(t,1), size(t,2), &
                 size(t,3), num_uw_tracers) :: qtruw , &
                                                    uw_tracers
real, dimension (size(t,1), size(t,2), &
                 size(t,3), num_mca_tracers) :: qtrmca , &
                                                mca_tracers
!     tr_flux        tracer fux across the surface
!                    [ kg(tracer) / (m**2 sec) ]
real, dimension (size(t,1), size(t,2), &
                 num_donner_tracers) :: tr_flux  

! end tracer code

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

real, dimension(size(t,1),size(t,2),size(t,3)) ::  cnvcntq
real, dimension(size(t,1),size(t,2)          ) ::  adjust_frac      
real, dimension(size(t,1),size(t,2),size(t,3)) ::  ttnd_adjustment
real, dimension(size(t,1),size(t,2),size(t,3)) ::  available_cf_for_uw

character(len=32) :: tracer_units, tracer_name
integer           :: secs, days
logical, dimension(size(t,1),size(t,2)) ::  conv_calc_completed
logical, dimension(size(t,1),size(t,2)) ::  ltemp

real,    dimension (size(t,1), size(t,2)) ::            &
                              enthint, lcondensint, enthdiffint,  &
                              vaporint, condensint, precipint, diffint
real :: temp_1, temp_2, temp_3, qnew 
real :: qrf, environmental_fraction, environmental_qv
real :: temp_don1, temp_don2, temp_don3, temp_uw_c
       
call mpp_clock_begin (excl_pre_conv)
!---------------------------------------------------------------------
!    verify that the module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('moist_processes_mod',  &
                 'moist_processes_init has not been called.', FATAL)
      endif

!-----------------------------------------------------------------------
      conc_air = 10. * pfull / (boltz * t)

!-------- input array size and position in global storage --------------

      ix=size(t,1); jx=size(t,2); kx=size(t,3); nt=size(rdt,4)

      conv_calc_completed = .false.
      available_cf_for_uw = 1.0

!--------------------------------------------------------------------
!    define the inverse of the time step.
!--------------------------------------------------------------------
      dtinv = 1.0/dt
!--------------------------------------------------------------------
!    initialize the arrays which will be output from this subroutine.
!--------------------------------------------------------------------

      conv_rain3d = 0.0
      conv_snow3d = 0.0
      do j=1,jx
       do i=1,ix
         lprec(i,j)         = 0.0  
         fprec(i,j)         = 0.0
         convect(i,j)       = .false.
         gust_cv(i,j)       = 0.0
         precip(i,j)        = 0.0 
       enddo
      enddo
      do k=1,kx
       do j=1,jx
        do i=1,ix
          diff_cu_mo(i,j,k)  = 0.0

!---------------------------------------------------------------------
!    initialize local arrays which will hold sums.
!---------------------------------------------------------------------
          tdt_init(i,j,k) = tdt(i,j,k)
          qdt_init(i,j,k) = qdt(i,j,k)
          ttnd_conv(i,j,k) = 0.
          qtnd_conv(i,j,k) = 0.
          qtnd(i,j,k) = 0.
          qltnd(i,j,k) = 0.
          qitnd(i,j,k) = 0.
          qatnd(i,j,k) = 0.
          if (do_liq_num) qntnd(i,j,k) = 0.

        enddo
       enddo
      enddo

!---------------------------------------------------------------------
!    define input fields to be used, either the tau time level fields,
!    or the tau - 1 time level values updated with the time tendencies
!    thus far calculated on the current step. control is through nml
!    variable use_tau.
!---------------------------------------------------------------------
      if (use_tau) then
       do k=1,kx
        do j=1,jx
         do i=1,ix
           tin(i,j,k) = t(i,j,k)
           qin(i,j,k) = q(i,j,k)
           uin(i,j,k) = u(i,j,k)
           vin(i,j,k) = v(i,j,k)
         enddo
        enddo
       enddo
       do tr=1,size(r,4)
        do k=1,kx
         do j=1,jx
          do i=1,ix

            tracer(i,j,k,tr) = r(i,j,k,tr)
            if (tr .le. nt) then
!---------------------------------------------------------------------
!    initialize local array rdt_init that will hold a sum
!---------------------------------------------------------------------
              rdt_init(i,j,k,tr) = rdt(i,j,k,tr)
            else
              rdt_init(i,j,k,tr) = 0.0
            endif
          enddo
         enddo
        enddo
       enddo  
      else
       do k=1,kx
        do j=1,jx
         do i=1,ix
           tin(i,j,k) = tm(i,j,k)+tdt(i,j,k)*dt
           qin(i,j,k) = qm(i,j,k)+qdt(i,j,k)*dt
           uin(i,j,k) = um(i,j,k)+udt(i,j,k)*dt
           vin(i,j,k) = vm(i,j,k)+vdt(i,j,k)*dt
         enddo
        enddo
       enddo
       do tr=1,size(r,4)
        do k=1,kx
         do j=1,jx
          do i=1,ix
 
            if (tr .le. nt) then
!---------------------------------------------------------------------
!    initialize local array rdt_init that will hold a sum
!---------------------------------------------------------------------
              rdt_init(i,j,k,tr) = rdt(i,j,k,tr)
              tracer(i,j,k,tr) = rm(i,j,k,tr)+rdt(i,j,k,tr)*dt
            else
              rdt_init(i,j,k,tr) = 0.0
              tracer(i,j,k,tr) = r(i,j,k,tr)
            endif
          enddo
         enddo
        enddo
       enddo  
      endif


!--------------------------------------------------------------------
!    if using eta vertical coordinate, define the appropriate values 
!    for any points located below the ground. values of 0.0 are given
!    to u, v and q, and a temperature value of EPST (=200. K) is given 
!    to sub-surface  points.
!--------------------------------------------------------------------
      if (present(mask) .and. present(kbot))  then
        do k=1,kx
         do j=1,jx
          do i=1,ix
           tin(i,j,k) = mask(i,j,k)*tin(i,j,k) + (1.0 - mask(i,j,k))*EPST 
           qin(i,j,k) = mask(i,j,k)*qin(i,j,k)
           uin(i,j,k) = mask(i,j,k)*uin(i,j,k)
           vin(i,j,k) = mask(i,j,k)*vin(i,j,k)
          enddo
         enddo
        enddo
        do tr=1,size(r,4)
          tracer(:,:,:,tr) = mask(:,:,:)*tracer(:,:,:,tr)
        end do  
      endif
   
!----------------------------------------------------------------------
!    compute the mass in each model layer.
!----------------------------------------------------------------------
      do k=1,kx
        pmass(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/GRAV
      end do

!----------------------------------------------------------------------
!    output any requested convectively-transported tracer fields 
!    and / or their column sums before convective transport.
!----------------------------------------------------------------------
      do n=1,num_tracers
        used = send_data (id_conv_tracer(n), tracer(:,:,:,n),   &
                          Time, is, js, 1, rmask=mask)
        if (id_conv_tracer_col(n) > 0) then
          tempdiag(:,:)=0.
          do k=1,kx
            tempdiag(:,:) = tempdiag(:,:) + tracer(:,:,k,n)*pmass(:,:,k)
          end do
          used = send_data (id_conv_tracer_col(n), tempdiag, Time, is, js)
        endif
      end do

!----------------------------------------------------------------------
!    compute the mean temperature in the lower atmosphere (the lowest
!    pdepth Pa), to be used to determine whether rain or snow reaches
!    the surface. define a logical variable coldT indicating whether
!    snow or rain falls in the column.
!    ????    SHOULD TIN BE USED RATHER THAN t ??
!----------------------------------------------------------------------
call mpp_clock_end   (excl_pre_conv)
      call tempavg (pdepth, phalf, t, tsnow, mask)
call mpp_clock_begin (excl_pre_conv)
      where (tsnow <= TFREEZE)
        coldT = .true.
      elsewhere
        coldT = .false.
      endwhere
      
!---------------------------------------------------------------------
!    begin the clock timing the dry and moist convection parameter-
!    izations.
!---------------------------------------------------------------------
      call mpp_clock_begin (convection_clock)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                   DRY CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    if dry adjustment is desired call subroutine dry_adj to obtain
!    the temperature tendencie swhich must be applied to adjust each
!    column to a non-superadiabatic lapse rate. 
!---------------------------------------------------------------------
      if (do_dryadj) then
call mpp_clock_end   (excl_pre_conv)
        call dry_adj (tin, pfull, phalf, delta_temp, mask)
call mpp_clock_begin (excl_pre_conv)

!-------------------------------------------------------------------
!    add the temperature change due to dry adjustment to the current
!    temperature. convert the temperature change to a heating rate and
!    add that to the temperature temndency array accumulating the ten-
!    dencies due to all physics processes.
!-------------------------------------------------------------------
        do k=1,kx
         do j=1,jx
          do i=1,ix
            tin(i,j,k)  = tin(i,j,k) + delta_temp(i,j,k)
            ttnd(i,j,k) = delta_temp(i,j,k)*dtinv
            tdt(i,j,k)  = tdt(i,j,k) + ttnd(i,j,k)
!---------------------------------------------------------------------
!    add the temperature time tendency from dry adjustment to the array
!    accumulating the total temperature time tendency from convection.
!---------------------------------------------------------------------
            ttnd_conv(i,j,k) = ttnd_conv(i,j,k) + ttnd(i,j,k)
          enddo
         enddo
        enddo

!---------------------------------------------------------------------
!    output the temperature tendency from dry adjustment, if desired.
!---------------------------------------------------------------------
        used = send_data (id_tdt_dadj, ttnd, Time, is, js, 1, rmask=mask )

      endif


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                  MOIST CONVECTION PARAMETERIZATIONS
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                0. UW SHALLOW CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

call mpp_clock_end   (excl_pre_conv)
call mpp_clock_begin (excl_shallow)
         call mpp_clock_begin (shallowcu_clock)
       cmf = 0.
     tracer_orig = tracer
  if (.not. do_donner_before_uw) then
   if (do_uw_conv) then
!---------------------------------------------------------------------
!    be sure all optional arguments associated with the uw_conv param-
!    eterization are present.
!---------------------------------------------------------------------
        if    &
           (present (shallow_cloud_area) .and.   &
            present (shallow_liquid) .and.   &
            present (shallow_ice) .and.  &
            present ( shallow_droplet_number) ) then
        else
         call error_mesg ('moist_processes_mod', 'moist_processes: &
                &not all 4 optional arguments needed for uw_conv &
              &output are present', FATAL)
       endif

!----------------------------------------------------------------------
!    if any tracers are to be transported by UW convection, check each
!    active tracer to find those to be transported and fill the 
!    ras_tracers array with these fields.
!---------------------------------------------------------------------
       nn = 1
       do n=1, num_tracers
         if (tracers_in_uw(n)) then
            uw_tracers(:,:,:,nn) = tracer(:,:,:,n)
           nn = nn + 1
         endif
      end do

call mpp_clock_end   (excl_shallow)
         call uw_conv (is, js, Time, tin, qin, uin, vin, pfull, phalf,zfull,       & !input
              zhalf, tracer, omega, dt, pblht, ustar, bstar, qstar, land, coldT,   & !input
!              Aerosol, cush, do_strat,                                                     & !input
              Aerosol, cush, do_strat,  conv_calc_completed,  & !input
              available_cf_for_uw, &
              ttnd_uw, qtnd_uw, qltnd_uw, qitnd_uw, qatnd_uw, qntnd_uw, utnd_uw, vtnd_uw, rain_uw, snow_uw,      & !output
              cmf, thlflx, qtflx, precflx, shallow_liquid, shallow_ice,&
             shallow_cloud_area, shallow_droplet_number, cbmf,       &  !output
!!5 miz does not wanty cbmf_clo as argument -- it is cbmf (intent in).
!            cbmf_clo, &
!++++yim
              uw_tracers, qtruw)                           !output
call mpp_clock_begin (excl_shallow)

    if (.not. do_limit_uw) then

        do k=1,kx
         do j=1,jx
          do i=1,ix
            tdt(i,j,k)=tdt(i,j,k)+ttnd_uw(i,j,k)
            qdt(i,j,k)=qdt(i,j,k)+qtnd_uw(i,j,k)
            udt(i,j,k)=udt(i,j,k)+utnd_uw(i,j,k)  
            vdt(i,j,k)=vdt(i,j,k)+vtnd_uw(i,j,k)

            ttnd_conv(i,j,k) = ttnd_conv(i,j,k) + ttnd_uw(i,j,k)
            qtnd_conv(i,j,k) = qtnd_conv(i,j,k) + qtnd_uw(i,j,k)
            if (do_strat) then
              rdt(i,j,k,nql) = rdt(i,j,k,nql) + qltnd_uw(i,j,k)
              rdt(i,j,k,nqi) = rdt(i,j,k,nqi) + qitnd_uw(i,j,k)
              rdt(i,j,k,nqa) = rdt(i,j,k,nqa) + qatnd_uw(i,j,k)
              if (do_liq_num) rdt(i,j,k,nqn) = rdt(i,j,k,nqn) + qntnd_uw(i,j,k)
            endif
          enddo
         enddo
        enddo

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from uw transport.
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_tracers
          if (tracers_in_uw(n)) then
            rdt(:,:,:,n) = rdt(:,:,:,n) + qtruw (:,:,:,nn)
             nn = nn + 1
          endif
        end do

        do j=1,jx
         do i=1,ix
           lprec(i,j)=lprec(i,j)+rain_uw(i,j)
           fprec(i,j)=fprec(i,j)+snow_uw(i,j)
           precip(i,j)=precip(i,j)+rain_uw(i,j)+snow_uw(i,j)
         enddo
        enddo
    endif  !(.not. do_limit_uw)

   endif  !(do_uw_conv)

  else
       tin_orig = tin
       qin_orig = qin
!     tracer_orig = tracer

  endif  ! (.not do_donner_before_uw)

call mpp_clock_end   (excl_shallow)
         call mpp_clock_end   (shallowcu_clock)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                A. DONNER DEEP CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!IF (.not. do_unified_convective_closure) then
!     cbmf_clo = 0.0
!ENDIF

!---------------------------------------------------------------------
!    if donner_deep convection is activated, execute the following code.
!---------------------------------------------------------------------
      if (do_donner_deep) then
call mpp_clock_begin (excl_donner_deep)
        
!---------------------------------------------------------------------
!    be sure all optional arguments associated with the donner param-
!    eterization are present.
!---------------------------------------------------------------------
        if    &
           (present (cell_cld_frac) .and.   &
            present (cell_liq_amt) .and. present ( cell_liq_size) .and. &
            present (cell_ice_amt) .and. present ( cell_ice_size) .and. &
            present (cell_droplet_number) .and. &
            present (meso_cld_frac) .and.   &
            present (meso_liq_amt) .and. present ( meso_liq_size) .and. &
            present (meso_ice_amt) .and. present ( meso_ice_size) .and. &
            present (meso_droplet_number) .and. &
            present (nsum_out) ) then
        else
          call error_mesg ('moist_processes_mod', 'moist_processes: &
               &not all 13 optional arguments needed for donner_deep &
              &output are present', FATAL)
        endif

        do k=1,kx
         do j=1,jx
          do i=1,ix
!--------------------------------------------------------------------
!    if strat_cloud_mod is activated, define the cloud liquid and 
!    cloud ice specific humidities and cloud area associated with 
!    strat_cloud_mod, so that they may be input to donner_deep_mod. 
!    if strat_cloud_mod is not activated, define these arrays to be 
!    zero. 
!--------------------------------------------------------------------
            if (do_strat) then
              qlin(i,j,k) = tracer(i,j,k,nql)
              qiin(i,j,k) = tracer(i,j,k,nqi)
              qain(i,j,k) = tracer(i,j,k,nqa)
           endif

!--------------------------------------------------------------------
!    convert vapor specific humidity to vapor mixing ratio so it may
!    be input to donner_deep_mod.
!--------------------------------------------------------------------
            rin(i,j,k) = qin(i,j,k)/(1.0 - qin(i,j,k))
          enddo
         enddo
        enddo

!---------------------------------------------------------------------
!    if any tracers are to be transported by donner convection, 
!    check each active tracer to find those to be transported and fill 
!    the donner_tracers array with these fields.
!---------------------------------------------------------------------
        donner_tracers(:,:,:,:) = 0.0
        nn = 1
        do n=1,num_tracers
          if (tracers_in_donner(n)) then
            donner_tracers(:,:,:,nn) = tracer(:,:,:,n)
            nn = nn + 1
          endif
        end do

!---------------------------------------------------------------------
!  NOTE 1: sfc_sh_flux, sfc_vapor_flux, tr_flux are the surface fluxes
!          that will have been obtained from the flux exchange module
!          and passed on to moist_processes and then to donner_deep.
!          FOR NOW, these values are defined herein, and given
!          values of 0.0
!---------------------------------------------------------------------
!       sfc_sh_flux(:,:) = INPUT_SFC_SH_FLUX_FROM_COUPLER(:,:)
!       sfc_vapor_flux(:,:) = INPUT_SFC_VAPOR_FLUX_FROM_COUPLER(:,:)
        sfc_sh_flux = 0.0
        sfc_vapor_flux = 0.0
        tr_flux(:,:,:) = 0.0
        nn = 1
        do n=1,num_tracers
          if (tracers_in_donner(n)) then
!           tr_flux(:,:,nn) = INPUT_SFC_FLUX_FROM_COUPLER(:,:,n)
            tr_flux(:,:,nn) = 0.0                                 
            nn = nn + 1
          endif
        end do
        do j=1,jx
         do i=1,ix
          if (pblht(i,j).lt.0.) then
           temp_1=0.0
          elseif (pblht(i,j).gt.5000.) then
           temp_1=5000.
          else
           temp_1=pblht(i,j)
          endif
          temp_2=ustar(i,j)**3.+0.6*ustar(i,j)*bstar(i,j)*temp_1
          if (temp_2 .gt. 0.) temp_2 = temp_2**(2./3.)
          tkemiz(i,j) = MAX (1.e-6, temp_2)
         enddo
        enddo

!---------------------------------------------------------------------
!    call donner_deep to compute the effects of deep convection on the 
!    temperature, vapor mixing ratio, tracers, cloud liquid, cloud ice
!    cloud area and precipitation fields.
!---------------------------------------------------------------------
        call mpp_clock_begin (donner_clock)
        call get_time (Time, secs, days)
call mpp_clock_end   (excl_donner_deep)
        if (do_strat) then
          call donner_deep (is, ie, js, je, dt, tin, rin, pfull,       &
                            phalf, zfull, zhalf, omega, pblht, tkemiz, &
                            qstar, cush, coldT, land, sfc_sh_flux,   &!miz
                            sfc_vapor_flux, tr_flux, donner_tracers, &
!!5 miz replaces cbmf_clo with cbmf
!                           Time, cbmf_clo,    &
                            secs, days, cbmf,    &
                            cell_cld_frac, cell_liq_amt, cell_liq_size, &
                            cell_ice_amt, cell_ice_size, &
                            cell_droplet_number, &
                            meso_cld_frac, meso_liq_amt, meso_liq_size, &
                            meso_ice_amt, meso_ice_size,  &
                            meso_droplet_number, nsum_out,  &
                            precip_returned, delta_temp, delta_vapor,  &
                            m_cdet_donner, m_cellup, mc_donner,  &
                            donner_humidity_area,  &
                            donner_humidity_factor, qtrtnd,    &
                            lheat_precip, vert_motion,        &
                            total_precip, liquid_precip, frozen_precip,&
                            qlin, qiin, qain, delta_ql, &   ! optional
                            delta_qi, delta_qa)             ! optional
        else
          call donner_deep (is, ie, js, je, dt, tin, rin, pfull,      &
                            phalf, zfull, zhalf, omega, pblht, tkemiz, &
                            qstar, cush, coldT, land, sfc_sh_flux,   &!miz
                            sfc_vapor_flux, tr_flux, donner_tracers, &
!!5 miz replaces cbmf_clo with cbmf
!                           Time,  cbmf_clo,                &
                            secs, days,  cbmf,                &
                            cell_cld_frac, cell_liq_amt, cell_liq_size, &
                            cell_ice_amt, cell_ice_size, &
                            cell_droplet_number, &
                            meso_cld_frac, meso_liq_amt, meso_liq_size, &
                            meso_ice_amt, meso_ice_size,  &
                            meso_droplet_number, nsum_out,  &
                            precip_returned, delta_temp, delta_vapor,  &
                            m_cdet_donner, m_cellup, mc_donner,  &
                            donner_humidity_area,  &
                            donner_humidity_factor, qtrtnd, &
                            lheat_precip, vert_motion,        &
                            total_precip, liquid_precip, frozen_precip)
        endif
        call mpp_clock_end (donner_clock)
call mpp_clock_begin (excl_donner_deep)


!---------------------------------------------------------------------
!    update the current timestep tracer changes with the contributions 
!    just obtained from donner transport.
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_tracers
          if (tracers_in_donner(n)) then
            rdt(:,:,:,n) = rdt(:,:,:,n) + qtrtnd(:,:,:,nn)
            nn = nn + 1
          endif
        end do

      if (do_donner_conservation_checks) then
       do j=1,jx
        do i=1,ix
          vaporint(i,j) = 0.
          enthint(i,j) = 0.
          condensint(i,j) = 0.
          lcondensint(i,j) = 0.
        end do
       end do
        
       do k=1,kx
        do j=1,jx
         do i=1,ix
          vaporint(i,j) = vaporint(i,j) + pmass(i,j,k)*delta_vapor(i,j,k)
          enthint(i,j) = enthint(i,j) + CP_AIR*pmass(i,j,k)*delta_temp(i,j,k)
          condensint(i,j) = condensint(i,j) + pmass(i,j,k) *  &
                            (delta_ql(i,j,k) + delta_qi(i,j,k))
          lcondensint(i,j) = lcondensint(i,j) + pmass(i,j,k) *  &
                             (HLV*delta_ql(i,j,k) + HLS*delta_qi(i,j,k))
         end do
        end do
       end do
       do j=1,jx
        do i=1,ix
          precipint(i,j) = total_precip(i,j)/seconds_per_day
          diffint(i,j) = (vaporint(i,j) + condensint(i,j))*dtinv  + precipint(i,j)
          enthdiffint(i,j) = (enthint(i,j) - lcondensint(i,j))*dtinv -    &
                   (lheat_precip(i,j) + vert_motion(i,j))/seconds_per_day 
          if (abs(enthdiffint(i,j)) > max_enthalpy_imbal_don(i,j)) then
            max_enthalpy_imbal_don(i,j) = abs (enthdiffint(i,j))
          endif
          if (abs(diffint(i,j)) > max_water_imbal_don(i,j)) then
            max_water_imbal_don(i,j) = abs (diffint(i,j))
          endif
        end do
       end do

 
       used = send_data(id_max_enthalpy_imbal_don, max_enthalpy_imbal_don, Time, is, js)
       used = send_data(id_max_water_imbal_don, max_water_imbal_don, Time, is, js)
       used = send_data(id_vaporint, vaporint*dtinv, Time, is, js)
       used = send_data(id_condensint, condensint*dtinv, Time, is, js)
       used = send_data(id_vertmotion, vert_motion/seconds_per_day, Time, is, js)
       used = send_data(id_precipint, precipint, Time, is, js)
       used = send_data(id_diffint, diffint, Time, is, js)
       used = send_data(id_enthint, enthint*dtinv, Time, is, js)
       used = send_data(id_lcondensint, lcondensint*dtinv, Time, is, js)
       used = send_data(id_lprcp, lheat_precip/seconds_per_day, Time, is, js)
       used = send_data(id_enthdiffint, enthdiffint, Time, is, js)
      endif
call mpp_clock_end   (excl_donner_deep)
call mpp_clock_begin (excl_donner_sph)

!--------------------------------------------------------------------
!    obtain updated vapor specific humidity (qnew) resulting from deep 
!    convection. define the vapor specific humidity change due to deep 
!    convection (qtnd).
!--------------------------------------------------------------------
        do k=1,kx
         do j=1,jx
          do i=1,ix
            if (delta_vapor(i,j,k) /= 0.0) then
              qnew = (rin(i,j,k) + delta_vapor(i,j,k))/   &
                            (1.0 + (rin(i,j,k) + delta_vapor(i,j,k)))
              delta_q(i,j,k) = qnew - qin(i,j,k)
            else
!              qnew  = qin(i,j,k)
              delta_q(i,j,k) = 0.
            endif
          end do
         end do
        end do

!---------------------------------------------------------------------
!    scale Donner tendencies to prevent the formation of negative
!    total water specific humidities
!---------------------------------------------------------------------

        if (do_strat .and. do_limit_donner) then

          scale_donner = HUGE(1.0)
          do k = 1,kx
           do j = 1,jx
            do i = 1,ix

!         Tendencies coming out of Donner deep are adjusted to prevent
!         the formation of negative water vapor, liquid or ice.

!         (1) Prevent negative liquid and ice specific humidities after
!             tendencies are applied

             if ((qlin(i,j,k)+delta_ql(i,j,k)) .lt. 0.) then
              delta_temp(i,j,k)  = delta_temp (i,j,k) - (qlin(i,j,k)+delta_ql(i,j,k))*HLV/CP_AIR
              delta_q    (i,j,k) = delta_q    (i,j,k) + (qlin(i,j,k)+delta_ql(i,j,k))
              delta_ql(i,j,k)    = delta_ql   (i,j,k) - (qlin(i,j,k)+delta_ql(i,j,k))
             end if

             if ((qiin(i,j,k)+delta_qi(i,j,k)) .lt. 0.) then
              delta_temp(i,j,k)  = delta_temp (i,j,k) - (qiin(i,j,k)+delta_qi(i,j,k))*HLS/CP_AIR
              delta_q    (i,j,k) = delta_q    (i,j,k) + (qiin(i,j,k)+delta_qi(i,j,k))
              delta_qi(i,j,k)    = delta_qi   (i,j,k) - (qiin(i,j,k)+delta_qi(i,j,k))
             end if

             if (abs(delta_ql(i,j,k) + delta_qi(i,j,k)) .lt. 1.e-10 ) then
              delta_qa(i,j,k) = 0.0
             end if

!         (2) Compute limit on Donner tendencies to prevent water vapor
!         from going below 1.e-10. The value of 1.e-10 is consistent with qmin
!         in strat_cloud.F90

!         scaling factor for each grid point
             if ( delta_q(i,j,k).lt.0 .and. (qin(i,j,k)+delta_q(i,j,k)).lt.1.e-10 ) then
              temp_1 = max( 0.0, -(qin(i,j,k)-1.e-10)/delta_q(i,j,k) )
             else
              temp_1 = 1.0
             end if

!         scaling factor for each column is the minimum value within that column
             scale_donner(i,j) = min ( temp_1, scale_donner(i,j) )

            enddo
           enddo
          enddo

!         scale tendencies

          do k=1,kx
           do j=1,jx
            do i=1,ix
              delta_temp(i,j,k)  = scale_donner(i,j) * delta_temp(i,j,k)
              delta_q (i,j,k)    = scale_donner(i,j) * delta_q (i,j,k)
              delta_qa(i,j,k)    = scale_donner(i,j) * delta_qa(i,j,k)
              delta_ql(i,j,k)    = scale_donner(i,j) * delta_ql(i,j,k)
              delta_qi(i,j,k)    = scale_donner(i,j) * delta_qi(i,j,k)
              liquid_precip(i,j,k) = scale_donner(i,j)*liquid_precip(i,j,k)
              frozen_precip(i,j,k) = scale_donner(i,j)*frozen_precip(i,j,k)
            end do
           end do
          end do

          nn = 1
          do n=1, num_tracers
            if (tracers_in_donner(n)) then
              do k=1,kx
                qtrtnd(:,:,k,nn) = scale_donner(:,:) * qtrtnd(:,:,k,nn)
              end do
              nn = nn + 1
            endif
          end do

          do j=1,jx
           do i=1,ix
             precip_returned(i,j) = scale_donner(i,j)*precip_returned(i,j)
             total_precip(i,j) = scale_donner(i,j)*total_precip(i,j) 
             lheat_precip(i,j) = scale_donner(i,j)*lheat_precip(i,j)
           end do
          end do

        else

          scale_donner = 1.0

        end if ! (do_limit_donner)

call mpp_clock_end   (excl_donner_sph)
call mpp_clock_begin (excl_donner_moist_cons)

!---------------------------------------------------------------------
!    recalculate the precip using the delta specific humidity tenden-
!    cies. define precip_adjustment as the change in precipitation 
!    resulting from the recalculation.
!---------------------------------------------------------------------
    if (force_donner_moist_conserv) then

!---------------------------------------------------------------------
!    calculate the precipitation needed to balance the change in water
!    content in the column.
!---------------------------------------------------------------------
     tempdiag = 0.
     do k=1,kx
      do j=1,jx
       do i=1,ix
         tempdiag(i,j) = tempdiag(i,j) -  &
                 (delta_q(i,j,k) + delta_ql(i,j,k) + delta_qi(i,j,k))* pmass(i,j,k)
         if (k.eq.kx) then
           temp_2 = (tempdiag(i,j) - precip_returned(i,j))
           if (ABS(temp_2) < 1.0e-10) then
             precip_adjustment (i,j) = 0.0
           else
             precip_adjustment (i,j) = temp_2
           endif
!----------------------------------------------------------------------
!    now adjust the temperature change to balance the precip adjustment
!    and so conserve enthalpy in the column.
!--------------------------------------------------------------------- 
           if (precip_returned(i,j) > 0.0) then
             temp_1 = precip_adjustment(i,j)/precip_returned(i,j)
           else
             temp_1 = 0.
           endif
           adjust_frac(i,j) = temp_1
         endif
       end do
      end do
     end do

     do k=1,kx
      do j=1,jx
       do i=1,ix
         temp_2 = liquid_precip(i,j,k)*adjust_frac(i,j)
         temp_3 = frozen_precip(i,j,k)*adjust_frac(i,j)
         ttnd_adjustment(i,j,k) = ((HLV*temp_2 + HLS*temp_3) * &
                                dt/seconds_per_day)/CP_AIR
         liquid_precip(i,j,k) = liquid_precip(i,j,k) + temp_2
         frozen_precip(i,j,k) = frozen_precip(i,j,k) + temp_3
       end do
      end do
     end do
  else ! (force_donner_moist_conserv)
     precip_adjustment(:,:) = 0.0
     adjust_frac      (:,:) = 0.0
      ttnd_adjustment(:,:,:) = 0.
  endif  ! (force_donner_moist_conserv)

  rain_don = 0.
  snow_don = 0.
  do k=1,kx
   do j=1,jx
    do i=1,ix
      rain_don(i,j) = rain_don(i,j) + liquid_precip(i,j,k)*pmass(i,j,k) &
                      /seconds_per_day
      snow_don(i,j) = snow_don(i,j) + frozen_precip(i,j,k)*pmass(i,j,k) &
                      /seconds_per_day

!---------------------------------------------------------------------
!    convert the changes in temperature, vapor specific humidity and 
!    precipitation resulting from deep convection to time tendencies 
!    of these quantities.
!---------------------------------------------------------------------
      ttnd_don(i,j,k) = (delta_temp(i,j,k) + ttnd_adjustment(i,j,k))*dtinv
      qtnd_don(i,j,k) = delta_q(i,j,k)*dtinv
    end do
   end do
  end do

!--------------------------------------------------------------------
!    output the time tendencies of temperature, vapor specific humid-
!    ity, precipitation and mass flux due to deep convection.
!--------------------------------------------------------------------
       used = send_data (id_tdt_deep_donner, ttnd_don, Time, is, js, 1, rmask=mask )
       used = send_data (id_qdt_deep_donner, qtnd_don, Time, is, js, 1, rmask=mask )
       used = send_data (id_qadt_deep_donner, delta_qa*dtinv, Time, is, js, 1, rmask=mask )
       used = send_data (id_qldt_deep_donner, delta_ql*dtinv, Time, is, js, 1, rmask=mask )
       used = send_data (id_qidt_deep_donner, delta_qi*dtinv, Time, is, js, 1, rmask=mask )
       used = send_data (id_prec1_deep_donner, precip_adjustment, &
                         Time, is, js, mask = precip_returned > 0.0)
       used = send_data (id_snow_deep_donner, snow_don, Time, is, js)
       used = send_data (id_mc_donner, mc_donner, Time, is, js, 1, rmask=mask )
       used = send_data (id_m_cdet_donner, m_cdet_donner, Time, is, js, 1, rmask=mask )
       used = send_data (id_m_cellup, m_cellup, Time, is, js, 1, rmask=mask )

       if (do_donner_conservation_checks) then
          tempdiag(:,:) = 0.
          tempdiag4(:,:) = 0.
          tempdiag5(:,:) = 0.
          tempdiag6(:,:) = 0.
        
          do k=1,kx
           do j=1,jx
            do i=1,ix

              temp_1 = (HLV*delta_ql(i,j,k) + HLS*delta_qi(i,j,k)) &
                       *dtinv*pmass(i,j,k)
              temp_2 = CP_AIR*ttnd_don(i,j,k) * pmass(i,j,k)
           
              tempdiag(i,j) = tempdiag(i,j) + temp_2 - temp_1 &
                    - (HLV*liquid_precip(i,j,k) + HLS*frozen_precip(i,j,k))  &
                    *pmass(i,j,k)/seconds_per_day
              tempdiag4(i,j) = tempdiag4(i,j) + temp_2
              tempdiag5(i,j) = tempdiag5(i,j) - temp_1
              tempdiag6(i,j) = tempdiag6(i,j) + CP_AIR*  &
                            ttnd_adjustment(i,j,k)*dtinv
            end do
           end do
          end do
 
          used = send_data (id_enth_donner_col, tempdiag, Time, is, js)
          used = send_data (id_enth_donner_col2, -hlv*rain_don, Time, is, js)
          used = send_data (id_enth_donner_col3, -hls*snow_don, Time, is, js)
          used = send_data (id_enth_donner_col4, tempdiag4, Time, is, js)
          used = send_data (id_enth_donner_col5, tempdiag5, Time, is, js)
          used = send_data (id_enth_donner_col6, tempdiag6, Time, is, js)
          used = send_data (id_enth_donner_col7, adjust_frac, Time, is, js)

          if (id_wat_donner_col > 0) then
            tempdiag(:,:) = rain_don + snow_don
            do k=1,kx
              tempdiag(:,:) = tempdiag(:,:)  + (qtnd_don(:,:,k)  + &
                              delta_ql(:,:,k)*dtinv +  &
                              delta_qi(:,:,k)*dtinv)*pmass(:,:,k)
            end do
            used = send_data (id_wat_donner_col, tempdiag, Time, is, js)
          endif
 
       endif ! (donner_conservation_checks)

       used = send_data (id_prec_deep_donner, rain_don + snow_don, Time, is, js )
 
!--------------------------------------------------------------------
!    save the tendencies of temperature and specific humidity resulting
!    from the deep convection component of the donner parameterization. 
!--------------------------------------------------------------------
        ttnd_conv = ttnd_conv + ttnd_don
        qtnd_conv = qtnd_conv + qtnd_don

         if (do_donner_mca) then
!--------------------------------------------------------------------
!    call subroutine moist_conv to handle any shallow convection 
!    present in the grid. in this call do_strat is always set to .false.
!    so that no convective detrainment (and corresponding change in
!    large-scale cloud amount and area) from moist convective adjustment
!    is allowed, consistent with this call being constrained to handle
!    shallow convection.
!--------------------------------------------------------------------
        call mpp_clock_begin (mca_clock)
        tin_pass = tin+delta_temp
        qin_pass = qin+delta_q
call mpp_clock_end   (excl_donner_moist_cons)
        call moist_conv (tin_pass, qin_pass, pfull, phalf, coldT, &
                         ttnd_donmca, qtnd_donmca, rain_donmca,  &
                         snow_donmca, dtinv, Time, is, js,    &
                         donner_tracers, qtrtnd, Lbot=kbot, mask=mask)           
        call mpp_clock_end (mca_clock)
call mpp_clock_begin (excl_donner_moist_cons)

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from moist convective adjustment. currently there
!    is no tracer transport by this process.
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_tracers
          if (tracers_in_donner(n)) then
            rdt(:,:,:,n) = rdt(:,:,:,n) + qtrtnd(:,:,:,nn)
            nn = nn + 1
          endif
        end do

!--------------------------------------------------------------------
!    output the time tendencies of temperature, vapor specific humid-
!    ity, precipitation and snow due to the moist convective 
!    adjustment pass of the donner parameterization.
!--------------------------------------------------------------------
        used = send_data (id_tdt_mca_donner, ttnd_donmca, Time, is, js, 1, &
                          rmask=mask)
        used = send_data (id_qdt_mca_donner, qtnd_donmca, Time, is, js, 1, &
                          rmask=mask)
        used = send_data (id_prec_mca_donner, rain_donmca+snow_donmca, Time,  &
                          is, js) 
        used = send_data (id_snow_mca_donner, snow_donmca, Time, is, js)

        if (id_enth_mca_donner_col > 0) then
          tempdiag(:,:) = -HLV*rain_donmca -HLS*snow_donmca
          do k=1,kx
            tempdiag(:,:) = tempdiag(:,:) + CP_AIR*ttnd_donmca(:,:,k)*pmass(:,:,k)
          end do
          used = send_data (id_enth_mca_donner_col, tempdiag, Time, is, js)
        endif

        if (id_wat_mca_donner_col > 0) then
          tempdiag(:,:) = rain_donmca + snow_donmca
          do k=1,kx
            tempdiag(:,:) = tempdiag(:,:) + qtnd_donmca(:,:,k)*pmass(:,:,k)
          end do
          used = send_data (id_wat_mca_donner_col, tempdiag, Time, is, js)
        endif

!--------------------------------------------------------------------
!------- diagnostics for tracers from convection -------
!  allow any tracer to be activated here (allows control cases)
!--------------------------------------------------------------------
      do n=1,num_tracers
         used = send_data (id_conv_tracer(n), tracer(:,:,:,n), Time, is, js, 1, &
                           rmask=mask )
!------- diagnostics for tracers column integral tendency ------
         if ( id_conv_tracer_col(n) > 0 ) then
           tempdiag(:,:)=0.
           do k=1,kx
             tempdiag(:,:) = tempdiag(:,:) + tracer   (:,:,k,n)*pmass(:,:,k)
           end do
           used = send_data ( id_conv_tracer_col(n), tempdiag, Time, is, js )
        end if
      enddo

 !--------------------------------------------------------------------
!    output the time tendencies of tracer and of column tracer 
!    due to the moist convective adjustment pass of the donner 
!    parameterization. currently moist convective adjustment does not
!    affect the tracer fields, so these fields are always 0.0.
!--------------------------------------------------------------------
        do n = 1, num_donner_tracers
          used = send_data (id_tracerdt_mcadon(n), qtrtnd(:,:,:,n), &
                            Time, is, js, 1, rmask=mask )
          if (id_tracerdt_mcadon_col(n) > 0 ) then
            tempdiag(:,:)=0.
            do k=1,kx
              tempdiag(:,:) = tempdiag(:,:) + qtrtnd(:,:,k,n)* &
                                                            pmass(:,:,k)
            end do
            used = send_data (id_tracerdt_mcadon_col(n), tempdiag,  &
                              Time, is, js )
          endif
        enddo

!--------------------------------------------------------------------
!    define the heating, moistening and precipitation rates as the sum 
!    of the contributions from the deep convection pass and the moist 
!    convective adjustment pass of the donner parameterization. if 
!    ras_mod is also activated, store these values in temporary arrays
!    until the contributions from ras_mod is calculated.
!--------------------------------------------------------------------
        ttnd_conv = ttnd_conv + ttnd_donmca
        qtnd_conv = qtnd_conv + qtnd_donmca
      endif
call mpp_clock_end   (excl_donner_moist_cons)


!---------------------------------------------------------------------
!    if donner_deep_mod is not active, define input fields normally 
!    produced by donner_deep_mod and needed by strat_cloud_mod 
!    appropriately.
!---------------------------------------------------------------------
      else   ! (do_donner_deep)
        mc_donner = 0.0
        m_cdet_donner = 0.0
        m_cellup = 0.0
        donner_humidity_area = 0.
        donner_humidity_factor = 0.
      endif  ! (do_donner_deep)

call mpp_clock_begin (excl_donner_deep_2)
! ADD TENDENCIES HERE, IN SAME AORDER AS ORIGINAL:
        if (do_donner_deep) then

!--------------------------------------------------------------------
!    add the contributions to the temperature and vapor specific 
!    humidity tendencies from the moist convective adjustment pass of
!    donner_deep_mod to the arrays accumulating the total tendencies 
!    due to all physics processes.
!--------------------------------------------------------------------
          if (do_donner_mca) then
           temp_1 = 1.0
           temp_2 = 0.0
          else
           temp_1 = 0.0
           temp_2 = 1.0
          endif

          do k=1,kx
           do j=1,jx
            do i=1,ix
              if (do_strat) then
                tracer(i,j,k,nql) = qlin(i,j,k) + delta_ql(i,j,k)
                tracer(i,j,k,nqi) = qiin(i,j,k) + delta_qi(i,j,k)
                tracer(i,j,k,nqa) = qain(i,j,k) + delta_qa(i,j,k)
                rdt(i,j,k,nql) = rdt(i,j,k,nql) + delta_ql(i,j,k)*dtinv
                rdt(i,j,k,nqi) = rdt(i,j,k,nqi) + delta_qi(i,j,k)*dtinv
                rdt(i,j,k,nqa) = rdt(i,j,k,nqa) + delta_qa(i,j,k)*dtinv
              endif
!--------------------------------------------------------------------
!    add the contributions to the temperature and vapor specific 
!    humidity tendencies from donner_deep mod to the arrays accumulating
!    the total tendencies due to all physics processes.
!--------------------------------------------------------------------
              tdt(i,j,k) = tdt(i,j,k) + ttnd_don(i,j,k) + temp_1*ttnd_donmca(i,j,k) 
              qdt(i,j,k) = qdt(i,j,k) + qtnd_don(i,j,k) + temp_1*qtnd_donmca(i,j,k)
!---------------------------------------------------------------------
!    update the values of temperature and vapor specific humidity to
!    include the effects of deep convection.
!---------------------------------------------------------------------
              tin(i,j,k) = temp_1*tin_pass(i,j,k) + temp_2*(tin(i,j,k)+delta_temp(i,j,k))
              qin(i,j,k) = temp_1*qin_pass(i,j,k) + temp_2*(qin(i,j,k)+delta_q(i,j,k))
            enddo
           enddo
          enddo

!--------------------------------------------------------------------
!    add the liquid (rain) and frozen (snow) precipitation generated by
!    deep convection on this step to the arrays accumulating precip-
!    itation from all sources (lprec, fprec).
!--------------------------------------------------------------------
          do j=1,jx
           do i=1,ix
             lprec(i,j)  = lprec(i,j) + rain_don(i,j) + temp_1 * rain_donmca(i,j)
             fprec(i,j)  = fprec(i,j) + snow_don(i,j) + temp_1 * snow_donmca(i,j)
           enddo
          enddo

          if (only_one_conv_scheme_per_column) then
             conv_calc_completed = (rain_don + snow_don) > 0.0
          endif
          if (limit_conv_cloud_frac) then
            ltemp = ANY(donner_humidity_area >= 0.999, dim = 3)
            where (ltemp(:,:)) conv_calc_completed(:,:) = .true.
            available_cf_for_uw = MAX(0.999 - donner_humidity_area, 0.0)
          endif

        endif   !(do_donner_deep)
call mpp_clock_end   (excl_donner_deep_2)
call mpp_clock_begin (excl_moist_conv_adj)

   if (do_donner_before_uw) then
     if (do_uw_conv) then
!---------------------------------------------------------------------
!    be sure all optional arguments associated with the uw_conv param-
!    eterization are present.
!---------------------------------------------------------------------
        if    &
           (present (shallow_cloud_area) .and.   &
            present (shallow_liquid) .and.   &
            present (shallow_ice) .and.  &
            present ( shallow_droplet_number) ) then
        else
         call error_mesg ('moist_processes_mod', 'moist_processes: &
                &not all 4 optional arguments needed for uw_conv &
                &output are present', FATAL)
        endif

        if (use_updated_profiles_for_uw) then 
!---------------------------------------------------------------------
!    update tracer fields with tendencies due to donner convection and 
!    wet deposition by donner deep precipitation.
!---------------------------------------------------------------------
          do n=1,size(rdt,4)
            if (n /= nsphum) then
              if (.not. do_strat .or. ( n /= nql .and. n /= nqi .and.   &
                   n /= nqa .and. n /= nqn) ) then
                tracer(:,:,:,n) = tracer_orig(:,:,:,n) +   &
                                    (rdt(:,:,:,n) - rdt_init(:,:,:,n)) *dt
              endif
            endif
          end do
        endif

!----------------------------------------------------------------------
!    if any tracers are to be transported by UW convection, check each
!    active tracer to find those to be transported and fill the 
!    uw_tracers array with these fields.
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_tracers
          if (tracers_in_uw(n)) then
             uw_tracers(:,:,:,nn) = tracer(:,:,:,n)
            nn = nn + 1
          endif
        end do

        if (use_updated_profiles_for_uw) then 
          call uw_conv (is, js, Time, tin, qin, uin, vin, pfull, phalf,zfull,       & !input
               zhalf, tracer, omega, dt, pblht, ustar, bstar, qstar, land, coldT,   & !input
               Aerosol, cush, do_strat, conv_calc_completed,                        & !input
               available_cf_for_uw, ttnd_uw, qtnd_uw, qltnd_uw, qitnd_uw, qatnd_uw, & !output
               qntnd_uw, utnd_uw, vtnd_uw, rain_uw, snow_uw,                        & !output
               cmf, thlflx, qtflx, precflx, shallow_liquid, shallow_ice,            &
               shallow_cloud_area, shallow_droplet_number, cbmf,                    & !output
!!5 miz does not wanty cbmf_clo as argument -- it is cbmf (intent in).
!              cbmf_clo, &
!++++yim
               uw_tracers, qtruw)                                                     !output
        else ! (use_updated_profiles_for_uw)
          call uw_conv (is, js, Time, tin_orig, qin_orig, uin, vin, pfull, phalf,zfull,  & !input
               zhalf, tracer_orig, omega, dt, pblht, ustar, bstar, qstar, land, coldT,   & !input
               Aerosol, cush, do_strat,  conv_calc_completed,                            & !input
               available_cf_for_uw, ttnd_uw, qtnd_uw, qltnd_uw, qitnd_uw, qatnd_uw,      & !output
               qntnd_uw, utnd_uw, vtnd_uw, rain_uw, snow_uw,                             & !output
               cmf, thlflx, qtflx, precflx, shallow_liquid, shallow_ice,                 &
               shallow_cloud_area, shallow_droplet_number, cbmf,                         & !output
!!5 miz does not wanty cbmf_clo as argument -- it is cbmf (intent in).
!              cbmf_clo, &
!++++yim
               uw_tracers, qtruw)                                                          !output
        endif ! (use_updated_profiles_for_uw)

        if (.not. do_limit_uw) then
          do k=1,kx
           do j=1,jx
            do i=1,ix
              tdt(i,j,k)=tdt(i,j,k)+ttnd_uw(i,j,k)
              qdt(i,j,k)=qdt(i,j,k)+qtnd_uw(i,j,k)
              udt(i,j,k)=udt(i,j,k)+utnd_uw(i,j,k)
              vdt(i,j,k)=vdt(i,j,k)+vtnd_uw(i,j,k)

              ttnd_conv(i,j,k) = ttnd_conv(i,j,k) + ttnd_uw(i,j,k)
              qtnd_conv(i,j,k) = qtnd_conv(i,j,k) + qtnd_uw(i,j,k)
              if (do_strat) then
                rdt(i,j,k,nql) = rdt(i,j,k,nql) + qltnd_uw(i,j,k)
                rdt(i,j,k,nqi) = rdt(i,j,k,nqi) + qitnd_uw(i,j,k)
                rdt(i,j,k,nqa) = rdt(i,j,k,nqa) + qatnd_uw(i,j,k)
                if (do_liq_num) rdt(i,j,k,nqn) = rdt(i,j,k,nqn) + qntnd_uw(i,j,k)
              endif
            enddo
           enddo
          enddo

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions
!    just obtained from uw transport.
!---------------------------------------------------------------------
          nn = 1
          do n=1, num_tracers
            if (tracers_in_uw(n)) then
              rdt(:,:,:,n) = rdt(:,:,:,n) + qtruw (:,:,:,nn)
               nn = nn + 1
            endif
          end do

          do j=1,jx
           do i=1,ix
             lprec(i,j)=lprec(i,j)+rain_uw(i,j)
             fprec(i,j)=fprec(i,j)+snow_uw(i,j)
             precip(i,j)=precip(i,j)+rain_uw(i,j)+snow_uw(i,j)
           enddo
          enddo
        endif  !(.not. do_limit_uw)

     endif   !(do_uw_conv)

   endif ! (do_donner_before_uw)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                B. MOIST CONVECTIVE ADJUSTMENT             
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!    execute moist convective adjustment if desired.
!-----------------------------------------------------------------------
      if (do_mca) then

!---------------------------------------------------------------------
!    check each active tracer to find any that are to be transported 
!    by moist convective adjustment and fill the mca_tracers array with
!    these fields.
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_tracers
          if (tracers_in_mca(n)) then
            mca_tracers(:,:,:,nn) = tracer(:,:,:,n)
            nn = nn + 1
          endif
        end do

!---------------------------------------------------------------------
!    call subroutine moist_conv to obtain the temperature, moisture
!    precipitation and tracer tendencies due to the moist convective
!    adjustment parameterization. currently there is no tracer tendency
!    due to this parameterization.
!---------------------------------------------------------------------
!++++yim Should also account for change in qn dut to moist convective adjustment.
call mpp_clock_end   (excl_moist_conv_adj)
        call mpp_clock_begin (mca_clock)
        if (do_strat) then
          call moist_conv (tin, qin, pfull, phalf, coldT, ttnd, qtnd,  &
                           rain, snow, dtinv, Time, is, js, mca_tracers,&
                           qtrmca,Lbot= kbot, mask=mask,  &
                           ql=tracer(:,:,:,nql), qi=tracer(:,:,:,nqi), &
                           cf=tracer(:,:,:,nqa), qldel=qltnd,  &
                           qidel=qitnd, cfdel=qatnd)
        else
          call moist_conv (tin, qin, pfull, phalf, coldT, ttnd, qtnd,  &
                           rain, snow, dtinv, Time, is, js, mca_tracers,&
                           qtrmca, Lbot=kbot, mask=mask)
        endif
        call mpp_clock_end (mca_clock)
call mpp_clock_begin (excl_moist_conv_adj)

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from moist convective adjustment. currently there
!    is no tracer transport by this process.
!    NOTE : the stratcloud tracers are updated within moist_conv.
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_tracers
          if (tracers_in_mca(n)) then
            rdt(:,:,:,n) = rdt(:,:,:,n) + qtrmca(:,:,:,nn)
            nn = nn + 1
          endif
        end do

        do k=1,kx
         do j=1,jx
          do i=1,ix
!----------------------------------------------------------------------
!    add the temperature and specific humidity tendencies from moist
!    convective adjustment (ttnd, qtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt).
!----------------------------------------------------------------------
            tdt(i,j,k) = tdt(i,j,k) + ttnd(i,j,k)
            qdt(i,j,k) = qdt(i,j,k) + qtnd(i,j,k)
            ttnd_conv(i,j,k) = ttnd_conv(i,j,k) + ttnd(i,j,k)
            qtnd_conv(i,j,k) = qtnd_conv(i,j,k) + qtnd(i,j,k)

!----------------------------------------------------------------------
!    if strat_cloud_mod is activated, add the cloud liquid, ice and area
!    tendencies from moist convective adjustment to the 
!    arrays accumulating these tendencies from all physics processes 
!    (rdt).
!----------------------------------------------------------------------
            if (do_strat) then
              rdt(i,j,k,nql) = rdt(i,j,k,nql) + qltnd(i,j,k)
              rdt(i,j,k,nqi) = rdt(i,j,k,nqi) + qitnd(i,j,k)
              rdt(i,j,k,nqa) = rdt(i,j,k,nqa) + qatnd(i,j,k)
            endif
          end do
         end do
        end do

!----------------------------------------------------------------------
!    increment the liquid, solid and total precipitation fields with 
!    the contribution from moist convective adjustment.
!----------------------------------------------------------------------
        do j=1,jx
         do i=1,ix
           lprec(i,j)  = lprec(i,j)  + rain(i,j)
           fprec(i,j)  = fprec(i,j)  + snow(i,j)
         end do
        end do
      endif ! (do_mca)
call mpp_clock_end   (excl_moist_conv_adj)
call mpp_clock_begin (excl_ras)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!           X. BETTS-MILLER CONVECTION SCHEME 
!			
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  if ( any((/do_bm,do_bmmass,do_bmomp/)) ) then

    if (do_bm) then ! betts-miller cumulus param scheme
      call betts_miller (dt,tin,qin,pfull,phalf,coldT,rain,snow,ttnd,qtnd,&
                        q_ref,bmflag,klzbs,cape,cin,t_ref,invtaubmt,&
                        invtaubmq, mask=mask)
    endif

    if (do_bmmass) then ! betts-miller-style massflux cumulus param scheme
      call bm_massflux (dt,tin,qin,pfull,phalf,coldT,rain,snow,ttnd,qtnd,&
                     q_ref,bmflag,klzbs,t_ref,massflux,&
                     mask=mask)
    endif

    if (do_bmomp) then ! olivier's betts-miller cumulus param scheme
      call bm_omp (dt,tin,qin,pfull,phalf,coldT,rain,snow,ttnd,qtnd,&
                      q_ref,bmflag,klzbs,t_ref, mask=mask)
    endif

!------- (update input values and) compute tendency -----
    tin=tin+ttnd;    qin=qin+qtnd
    ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
    rain=rain*dtinv; snow=snow*dtinv
                                                                                   
!-------- add on tendency ----------
    tdt=tdt+ttnd; qdt=qdt+qtnd
                                                                         
!------- compute rh clouds if desired ------
    if (do_rh_clouds) then
                                                                          
!calculate relative humidity
      call rh_calc(pfull,tin,qin,RH,mask)
                                                                            
!pass RH to rh_clouds_sum
      call rh_clouds_sum (is, js, RH) ! XXX  RH is not relative humidity when do_simple=.true.
                                                                     
    end if
!------- save total precip and snow ---------
    lprec=lprec+rain
    fprec=fprec+snow
    precip=precip+rain+snow
                                                                                           
!-----------------------------------------------------------------------
  endif ! if ( any((/do_bm,do_bmmass,do_bmomp/)) )


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!           C. RELAXED ARAKAWA-SCHUBERT PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!    execute relaxed arakawa/schubert cumulus parameterization scheme,
!    if desired.
!-----------------------------------------------------------------------
      if (do_ras) then

!----------------------------------------------------------------------
!    if any tracers are to be transported by ras convection, check each
!    active tracer to find those to be transported and fill the 
!    ras_tracers array with these fields.
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_tracers
          if (tracers_in_ras(n)) then
            ras_tracers(:,:,:,nn) = tracer(:,:,:,n)
            nn = nn + 1
          endif
        end do

!----------------------------------------------------------------------
!    call subroutine ras to obtain the temperature, specific humidity,
!    velocity, precipitation and tracer tendencies and mass flux 
!    associated with the relaxed arakawa-schubert parameterization.
!----------------------------------------------------------------------
        call mpp_clock_begin (ras_clock)
call mpp_clock_end   (excl_ras)
        if (do_strat .and. .not.(do_liq_num)) then
          call ras (is,   js,     Time,     tin,   qin,   &
                    uin,  vin,    pfull,    phalf, zhalf, coldT, &
                    dt,   ttnd,   qtnd,     utnd,  vtnd,  &
                    conv_rain3d, conv_snow3d, &
                    rain_ras, snow_ras,   ras_tracers, qtrras,    &
                    mask,  kbot, mc, det0,     &
                    tracer(:,:,:,nql), tracer(:,:,:,nqi), &
                    tracer(:,:,:,nqa), qltnd(:,:,:),&
                    qitnd(:,:,:), qatnd(:,:,:))       

       elseif(do_strat .and. do_liq_num) then
          call ras (is,   js,     Time,     tin,   qin,   &
                    uin,  vin,    pfull,    phalf, zhalf, coldT, &
                    dt,   ttnd,   qtnd,     utnd,  vtnd,  &
                    conv_rain3d, conv_snow3d, &
                    rain_ras, snow_ras,   ras_tracers, qtrras, &
                    mask,  kbot, mc, det0,     &
                    tracer(:,:,:,nql), tracer(:,:,:,nqi), tracer(:,:,:,nqa), &
                    qltnd(:,:,:), qitnd(:,:,:), qatnd(:,:,:),    &
                    tracer(:,:,:,nqn), qntnd(:,:,:), &
                    do_strat, Aerosol)
        else
          call ras (is,   js,     Time,     tin,   qin,          &
                    uin,  vin,    pfull,    phalf, zhalf, coldT, &
                    dt,   ttnd,   qtnd,     utnd,  vtnd,         &
                    conv_rain3d, conv_snow3d, &
                    rain_ras, snow_ras,   ras_tracers, qtrras,           &
                    mask,  kbot,  mc, det0)
        endif
        call mpp_clock_end (ras_clock)
call mpp_clock_begin (excl_ras)

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from ras transport.
!    NOTE : the stratcloud tracers are updated within ras.        
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_tracers
          if (tracers_in_ras(n)) then
            rdt(:,:,:,n) = rdt(:,:,:,n) + qtrras (:,:,:,nn)
            nn = nn + 1
          endif
        end do

        do k=1,kx
         do j=1,jx
          do i=1,ix
!----------------------------------------------------------------------
!    add the temperature, specific humidity and momentum tendencies 
!    from ras (ttnd, qtnd, utnd, vtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt, udt, vdt).
!----------------------------------------------------------------------
            tdt(i,j,k) = tdt(i,j,k) + ttnd(i,j,k)
            qdt(i,j,k) = qdt(i,j,k) + qtnd(i,j,k)
            udt(i,j,k) = udt(i,j,k) + utnd(i,j,k)
            vdt(i,j,k) = vdt(i,j,k) + vtnd(i,j,k)
            ttnd_conv(i,j,k) = ttnd_conv(i,j,k) + ttnd(i,j,k)
            qtnd_conv(i,j,k) = qtnd_conv(i,j,k) + qtnd(i,j,k)

!----------------------------------------------------------------------
!    if strat_cloud_mod is activated, add the cloud liquid, ice and area
!    tendencies from ras to the arrays accumulating these tendencies 
!    from all physics processes (rdt).
!----------------------------------------------------------------------
            if (do_strat) then
              rdt(i,j,k,nql) = rdt(i,j,k,nql) + qltnd(i,j,k)
              rdt(i,j,k,nqi) = rdt(i,j,k,nqi) + qitnd(i,j,k)
              rdt(i,j,k,nqa) = rdt(i,j,k,nqa) + qatnd(i,j,k)
              if (do_liq_num) rdt(i,j,k,nqn) = rdt(i,j,k,nqn) + qntnd(i,j,k)
            endif
          enddo
         enddo
        enddo

!----------------------------------------------------------------------
!    increment the liquid, solid and total precipitation fields with 
!    the contribution from ras.
!----------------------------------------------------------------------
        do j=1,jx
         do i=1,ix
           lprec(i,j)  = lprec(i,j)  + rain_ras(i,j)
           fprec(i,j)  = fprec(i,j)  + snow_ras(i,j)
         enddo
        enddo

!---------------------------------------------------------------------
!    if donner_deep_mod is also active, define the total time tendency 
!    due to all moist convective processes (donner (including its mca 
!    part), and ras) of temperature, specific humidity, rain, snow and,
!    if strat_cloud_mod is activated, the cloud liquid, cloud ice and 
!    cloud area. 
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!    if ras_mod is not activated, set the ras mass flux field to be 0.0.
!---------------------------------------------------------------------
      else
        mc(:,:,:) = 0.0
        det0(:,:,:) = 0.0
        rain_ras = 0.0
        snow_ras = 0.0
      endif  ! (do_ras)
call mpp_clock_end   (excl_ras)
call mpp_clock_begin (excl_cmt)

!---------------------------------------------------------------------
!    call subroutine cu_mo_trans if diffusive cumulus momentum 
!    transport is desired. 
!---------------------------------------------------------------------
        if (do_cmt) then
!  if doing nonlocal cmt, call cu_mo_trans for each convective scheme
!  separately
          if (.not. doing_diffusive) then
          if (cmt_uses_ras) then
            mc_cmt = mc
            det_cmt = det0
            call mpp_clock_begin (cmt_clock)
call mpp_clock_end   (excl_cmt)
            call cu_mo_trans (is, js, Time, mc_cmt, tin, phalf, pfull, &
                              zhalf, zfull, dt, uin, vin, tracer,   &
                              pmass, det_cmt, utnd, vtnd, ttnd,  &
                              qtrcumo, diff_cu_mo  )
            call mpp_clock_end (cmt_clock)
call mpp_clock_begin (excl_cmt)

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from cu_mo_trans.
!---------------------------------------------------------------------
            do n=1, num_tracers
              rdt(:,:,:,n) = rdt(:,:,:,n) + qtrcumo   (:,:,:,n)
            end do

!----------------------------------------------------------------------
!    add the temperature, specific humidity and momentum tendencies 
!    from ras (ttnd, qtnd, utnd, vtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt, udt, vdt).
!----------------------------------------------------------------------
          do k=1,kx
           do j=1,jx
            do i=1,ix
              tdt(i,j,k) = tdt(i,j,k) + ttnd(i,j,k) 
              udt(i,j,k) = udt(i,j,k) + utnd(i,j,k)
              vdt(i,j,k) = vdt(i,j,k) + vtnd(i,j,k)

!---------------------------------------------------------------------
!    if donner_deep_mod is also active, define the total time tendency 
!    due to all moist convective processes (donner (including its mca 
!    part), and ras) of temperature, specific humidity, rain, snow and,
!    if strat_cloud_mod is activated, the cloud liquid, cloud ice and 
!    cloud area. 
!---------------------------------------------------------------------
              ttnd_conv(i,j,k) = ttnd_conv(i,j,k) + ttnd(i,j,k)
            enddo
           enddo
          enddo
          endif
          if (cmt_uses_donner) then
            mc_cmt = m_cellup 
            det_cmt = m_cdet_donner 
            call mpp_clock_begin (cmt_clock)
call mpp_clock_end   (excl_cmt)
            call cu_mo_trans (is, js, Time, mc_cmt, tin, phalf, pfull, &
                              zhalf, zfull, dt, uin, vin, tracer,  &
                              pmass, det_cmt, utnd, vtnd, ttnd,   &
                              qtrcumo, diff_cu_mo  )
            call mpp_clock_end (cmt_clock)
call mpp_clock_begin (excl_cmt)
 
!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from cu_mo_trans.
!---------------------------------------------------------------------
            do n=1, num_tracers
              rdt(:,:,:,n) = rdt(:,:,:,n) + qtrcumo   (:,:,:,n)
            end do

!----------------------------------------------------------------------
!    add the temperature, specific humidity and momentum tendencies 
!    from ras (ttnd, qtnd, utnd, vtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt, udt, vdt).
!----------------------------------------------------------------------
          do k=1,kx
           do j=1,jx
            do i=1,ix
              tdt(i,j,k) = tdt(i,j,k) + ttnd(i,j,k)
              udt(i,j,k) = udt(i,j,k) + utnd(i,j,k)
              vdt(i,j,k) = vdt(i,j,k) + vtnd(i,j,k)
 
!---------------------------------------------------------------------
!    if donner_deep_mod is also active, define the total time tendency 
!    due to all moist convective processes (donner (including its mca 
!    part), and ras) of temperature, specific humidity, rain, snow and,
!    if strat_cloud_mod is activated, the cloud liquid, cloud ice and 
!    cloud area. 
!---------------------------------------------------------------------
              ttnd_conv(i,j,k) = ttnd_conv(i,j,k) + ttnd(i,j,k)
            enddo
           enddo
          enddo
          endif

          if (cmt_uses_uw) then
            mc_cmt(:,:,1) = 0.
            mc_cmt(:,:,kx+1) = 0.
            do k=2,kx
              mc_cmt(:,:,k) = cmf(:,:,k-1)
            end do
!   CURRENTLY no detrained mass flux provided from uw_conv; should only
!   use with 'diffusive' cmt scheme, not the non-local. (attempt to
!   use non-local will cause FATAL in _init routine.)
            det_cmt = 0.0   
call mpp_clock_end   (excl_cmt)
            call mpp_clock_begin (cmt_clock)
            call cu_mo_trans (is, js, Time, mc_cmt, tin, phalf, pfull, &
                              zhalf, zfull, dt, uin, vin, tracer,   &
                              pmass, det_cmt, utnd, vtnd, ttnd,  &
                              qtrcumo, diff_cu_mo  )
            call mpp_clock_end (cmt_clock)
call mpp_clock_begin (excl_cmt)

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from cu_mo_trans.
!---------------------------------------------------------------------
            do n=1, num_tracers
              rdt(:,:,:,n) = rdt(:,:,:,n) + qtrcumo   (:,:,:,n)
            end do

!----------------------------------------------------------------------
!    add the temperature, specific humidity and momentum tendencies 
!    from cmt due to uw (ttnd, qtnd, utnd, vtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt, udt, vdt).
!----------------------------------------------------------------------
          do k=1,kx
           do j=1,jx
            do i=1,ix
              tdt(i,j,k) = tdt(i,j,k) + ttnd(i,j,k) 
              udt(i,j,k) = udt(i,j,k) + utnd(i,j,k)
              vdt(i,j,k) = vdt(i,j,k) + vtnd(i,j,k)
              ttnd_conv(i,j,k) = ttnd_conv(i,j,k) + ttnd(i,j,k)
            enddo
           enddo
          enddo
          endif

         else ! (.not. doing_diffusive)
!  if using diffusive cmt, call cu_mo_trans once with combined mass
!  fluxes from all desired convective schemes.
          mc_cmt = 0.
          det_cmt = 0.
          if (cmt_uses_ras) then
            mc_cmt = mc_cmt + mc
          endif
          if (cmt_uses_donner) then
            mc_cmt = mc_cmt + m_cellup 
          endif
          if (cmt_uses_uw) then
            do k=2,kx
              mc_cmt(:,:,k) = mc_cmt(:,:,k) + cmf(:,:,k-1)
            end do
          endif
            call mpp_clock_begin (cmt_clock)
call mpp_clock_end   (excl_cmt)
            call cu_mo_trans (is, js, Time, mc_cmt, tin, phalf, pfull, &
                              zhalf, zfull, dt, uin, vin, tracer,   &
                              pmass, det_cmt, utnd, vtnd, ttnd,  &
                              qtrcumo, diff_cu_mo  )
            call mpp_clock_end (cmt_clock)
call mpp_clock_begin (excl_cmt)

!---------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    just obtained from cu_mo_trans.
!---------------------------------------------------------------------
            do n=1, num_tracers
              rdt(:,:,:,n) = rdt(:,:,:,n) + qtrcumo   (:,:,:,n)
            end do

!----------------------------------------------------------------------
!    add the temperature, specific humidity and momentum tendencies 
!    from ras (ttnd, qtnd, utnd, vtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt, udt, vdt).
!----------------------------------------------------------------------
          do k=1,kx
           do j=1,jx
            do i=1,ix
              tdt(i,j,k) = tdt(i,j,k) + ttnd(i,j,k) 
              udt(i,j,k) = udt(i,j,k) + utnd(i,j,k)
              vdt(i,j,k) = vdt(i,j,k) + vtnd(i,j,k)

!---------------------------------------------------------------------
!    if donner_deep_mod is also active, define the total time tendency 
!    due to all moist convective processes (donner (including its mca 
!    part), and ras) of temperature, specific humidity, rain, snow and,
!    if strat_cloud_mod is activated, the cloud liquid, cloud ice and 
!    cloud area. 
!---------------------------------------------------------------------
              ttnd_conv(i,j,k) = ttnd_conv(i,j,k) + ttnd(i,j,k)
            enddo
           enddo
          enddo

         endif ! (.not. doing_diffusive)
        endif  ! (do_cmt)
call mpp_clock_end   (excl_cmt)
call mpp_clock_begin (excl_tracer_dens)

!---------------------------------------------------------------------
!    calculate the tracer tendency due to wet deposition (wetdeptnd)
!    caused by the convectively generated precipitation (rain, snow) for
!    any tracers for which wet deposition has been activated. add this 
!    tendency to the tracer tendency due to all physics (rdt). save it 
!    also in an array which will be combined with any wet deposition 
!    resulting from large-scale precip producing the total wet deposition
!    for the tracer (wet_data).
!---------------------------------------------------------------------
      wet_data = 0.0
      cloud_wet = 1.e-3
      if (do_strat) then
         qtnd_wet = qtnd + qltnd + qitnd
      else
        qtnd_wet = qtnd
      end if
      cloud_frac = 0.1
      do n=1,size(rdt,4)
        if ( n /= nsphum ) then
          if ( .not. do_strat .or. (n /= nql .and. n /= nqi .and. n /= nqa .and. n /= nqn) ) then
            wetdeptnd = 0.0

call mpp_clock_end   (excl_tracer_dens)
!                   call wet_deposition( n, t, pfull, phalf, zfull, zhalf, rain, snow, &
                    call wet_deposition( n, t, pfull, phalf, zfull, zhalf, rain_ras, snow_ras, &
                                 qtnd_wet, cloud_wet, cloud_frac, &
                                 conv_rain3d, conv_snow3d, &
                                 tracer(:,:,:,n), wetdeptnd, &
                                 Time, 'convect', is, js, dt )
call mpp_clock_begin (excl_tracer_dens)
            rdt (:,:,:,n) = rdt(:,:,:,n) - wetdeptnd(:,:,:)
            wet_data(:,:,:,n) = wetdeptnd(:,:,:)
          endif
        endif  
      end do

      mc_full=0.; 
      do k=2,kx   
        mc_full(:,:,k) = 0.5*(mc(:,:,k) + mc(:,:,k+1)) +   &
                         0.5*(cmf(:,:,k)+cmf(:,:,k-1)) +   &
                              mc_donner(:,:,k)
      end do
!-----------------------------------------------------------------------
! lightning NOx parameterization
!-----------------------------------------------------------------------
if ( get_tracer_index(MODEL_ATMOS,'no') .ne. NO_TRACER ) then
   cldbot = 0
   cldtop = 0
   do i = 1,ix
   do j = 1,jx
      do k = 1,kx
         if (mc_full(i,j,k) /= 0 ) then
            cldtop(i,j) = k
            exit
         endif
      enddo
      do k = size(r,3),1,-1
         if (mc_full(i,j,k) /= 0 ) then
            cldbot(i,j) = k
            exit
         endif
      enddo
   enddo
   enddo
call mpp_clock_end   (excl_tracer_dens)
   call moz_hook(cldtop, cldbot, land, zfull, zhalf, t, prod_no, area, lat, &
                 Time, is, js)
call mpp_clock_begin (excl_tracer_dens)
   rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) = &
      rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) + &
      prod_no/conc_air
      used = send_data(id_prod_no,prod_no, Time, is_in=is, js_in=js)
endif

!-----------------------------------------------------------------------
!    define the total precipitation rate (precip).
!-----------------------------------------------------------------------
      precip(:,:) = lprec(:,:) + fprec(:,:)

!-----------------------------------------------------------------------
!    calculate convective gustiness, if desired.
!-----------------------------------------------------------------------
      if (do_gust_cv) then
        where((precip) > 0.0)
          gust_cv = gustmax*sqrt(   &
                           precip/(gustconst + precip) )
        end where
      end if

!---------------------------------------------------------------------
!    save a diagnostic indicating whether or not convection has occurred
!    within the column.
!---------------------------------------------------------------------
      where (precip > 0.) convect = .true.
call mpp_clock_end   (excl_tracer_dens)
call mpp_clock_begin (excl_uw_conv)


!---------------------------------------------------------------------
!    apply changes resulting from uw_conv
!---------------------------------------------------------------------
 
      if (do_uw_conv) then

       if (do_limit_uw) then

        scale_uw=HUGE(1.0)
        do k=1,kx
         do j=1,jx
          do i=1,ix

!       Tendencies coming out of UW shallow are adjusted to prevent
!       the formation of negative water vapor, liquid or ice.

!       (1) Prevent negative liquid and ice specific humidities after tendencies are applied

           temp_1 = tracer(i,j,k,nql)/dt + qltnd_uw(i,j,k)
           if (temp_1 .lt. 0.) then
            ttnd_uw (i,j,k) = ttnd_uw (i,j,k) - temp_1*HLV/CP_AIR
            qtnd_uw (i,j,k) = qtnd_uw (i,j,k) + temp_1
            qltnd_uw(i,j,k) = qltnd_uw(i,j,k) - temp_1
           end if

           temp_1 = tracer(i,j,k,nqi)/dt + qitnd_uw(i,j,k)
           if (temp_1 .lt. 0.) then
            ttnd_uw (i,j,k) = ttnd_uw (i,j,k) - temp_1*HLS/CP_AIR
            qtnd_uw (i,j,k) = qtnd_uw (i,j,k) + temp_1
            qitnd_uw(i,j,k) = qitnd_uw(i,j,k) - temp_1
           end if

           if (abs(qltnd_uw(i,j,k)+qitnd_uw(i,j,k))*dt .lt. 1.e-10 ) then
            qatnd_uw(i,j,k) = 0.0
           end if

!       (2) Compute limit on UW tendencies to prevent water vapor
!       from going below 1.e-10. The value of 1.e-10 is consistent with qmin
!       in strat_cloud.F90

!       scaling factor for each grid point
             temp_2 = qin(i,j,k) + tracer(i,j,k,nql) + tracer(i,j,k,nqi)
             temp_3  = ( qtnd_uw(i,j,k) + qltnd_uw(i,j,k) + qitnd_uw(i,j,k) )*dt
             if ( temp_3.lt.0 .and. temp_2+temp_3.lt.1.e-10 ) then
               temp_1 = max( 0.0, -(temp_2-1.e-10)/temp_3 )
             else
               temp_1 = 1.0
             end if

!       scaling factor for each column is the minimum value within that column
             scale_uw(i,j) = min( temp_1, scale_uw(i,j))

          end do
         end do
        end do

!       scale tendencies

        do k=1,kx
         do j=1,jx
          do i=1,ix
            utnd_uw(i,j,k)  = scale_uw(i,j) * utnd_uw(i,j,k)
            udt(i,j,k) = udt(i,j,k) + utnd_uw(i,j,k)
            uin(i,j,k) = uin(i,j,k) + utnd_uw(i,j,k)*dt

            vtnd_uw(i,j,k)  = scale_uw(i,j) * vtnd_uw(i,j,k)
            vdt(i,j,k) = vdt(i,j,k) + vtnd_uw(i,j,k)
            vin(i,j,k) = vin(i,j,k) + vtnd_uw(i,j,k)*dt

            ttnd_uw(i,j,k)  = scale_uw(i,j) * ttnd_uw(i,j,k)
            tdt(i,j,k) = tdt(i,j,k) + ttnd_uw(i,j,k)
            ttnd_conv(i,j,k) = ttnd_conv(i,j,k) + ttnd_uw(i,j,k)
            tin(i,j,k) = tin(i,j,k) + ttnd_uw(i,j,k)*dt

            qtnd_uw(i,j,k)  = scale_uw(i,j) * qtnd_uw(i,j,k)
            qdt(i,j,k) = qdt(i,j,k) + qtnd_uw(i,j,k)
            qtnd_conv(i,j,k) = qtnd_conv(i,j,k) + qtnd_uw(i,j,k)
            qin(i,j,k) = qin(i,j,k) + qtnd_uw(i,j,k)*dt

            qltnd_uw(i,j,k) = scale_uw(i,j) * qltnd_uw(i,j,k)
            tracer(i,j,k,nql) = tracer(i,j,k,nql) + qltnd_uw(i,j,k)*dt

            qitnd_uw(i,j,k) = scale_uw(i,j) * qitnd_uw(i,j,k)
            tracer(i,j,k,nqi) = tracer(i,j,k,nqi) + qitnd_uw(i,j,k)*dt

            qatnd_uw(i,j,k) = scale_uw(i,j) * qatnd_uw(i,j,k)
            tracer(i,j,k,nqa) = tracer(i,j,k,nqa) + qatnd_uw(i,j,k)*dt

            if (do_liq_num) then
              qntnd_uw(i,j,k) = scale_uw(i,j) * qntnd_uw(i,j,k)
              tracer(i,j,k,nqn) = tracer(i,j,k,nqn) + qntnd_uw(i,j,k)*dt
            endif

            if (do_strat) then
              rdt(i,j,k,nql) = rdt(i,j,k,nql) + qltnd_uw(i,j,k)
              rdt(i,j,k,nqi) = rdt(i,j,k,nqi) + qitnd_uw(i,j,k)
              rdt(i,j,k,nqa) = rdt(i,j,k,nqa) + qatnd_uw(i,j,k)
              if (do_liq_num) rdt(i,j,k,nqn) = rdt(i,j,k,nqn) + qntnd_uw(i,j,k)
            endif

          end do
         end do
        end do

        do j=1,jx
         do i=1,ix
           rain_uw(i,j) = scale_uw(i,j) * rain_uw(i,j)
           snow_uw(i,j) = scale_uw(i,j) * snow_uw(i,j)
!       update precipitation
           lprec(i,j) = lprec(i,j) + rain_uw(i,j)
           fprec(i,j) = fprec(i,j) + snow_uw(i,j)
           precip(i,j) = precip(i,j) + rain_uw(i,j) + snow_uw(i,j)
         end do
        end do

!       update the current tracer tendencies with the contributions 
!       obtained from uw transport.
        nn = 1
        do n=1, num_tracers
          if (tracers_in_uw(n)) then
             rdt(:,:,:,n) = rdt(:,:,:,n) + qtruw (:,:,:,nn)
              nn = nn + 1
          endif
        end do

       else

        scale_uw = 1.0

!       update input fields with changes from uw_conv
        do k=1,kx
         do j=1,jx
          do i=1,ix
            tin(i,j,k) = tin(i,j,k) + ttnd_uw(i,j,k)*dt
            qin(i,j,k) = qin(i,j,k) + qtnd_uw(i,j,k)*dt
            uin(i,j,k) = uin(i,j,k) + utnd_uw(i,j,k)*dt
            vin(i,j,k) = vin(i,j,k) + vtnd_uw(i,j,k)*dt
            tracer(i,j,k,nql) = tracer(i,j,k,nql) + qltnd_uw(i,j,k)*dt
            tracer(i,j,k,nqi) = tracer(i,j,k,nqi) + qitnd_uw(i,j,k)*dt
            tracer(i,j,k,nqa) = tracer(i,j,k,nqa) + qatnd_uw(i,j,k)*dt
            if (do_liq_num) tracer(i,j,k,nqn) = tracer(i,j,k,nqn) + qntnd_uw(i,j,k)*dt
          end do
         end do
        end do
       end if ! (do_limit_uw)

     endif ! (uw_conv)
call mpp_clock_end   (excl_uw_conv)
call mpp_clock_begin (excl_conv_diag)
 
!---------------------------------------------------------------------
!    update tracer fields with tendencies due to convection and wet 
!    deposition by convective precipitation.
!---------------------------------------------------------------------
     do n=1,size(rdt,4)
       if (n /= nsphum) then
         if (.not. do_strat .or. ( n /= nql .and. n /= nqi .and.   &
              n /= nqa .and. n /= nqn) ) then
!           tracer(:,:,:,n) = tracer(:,:,:,n) +   &
           tracer(:,:,:,n) = tracer_orig(:,:,:,n) +   &
                               (rdt(:,:,:,n) - rdt_init(:,:,:,n)) *dt
         endif
       endif
     end do

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                   CONVECTION DIAGNOSTICS      
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if ( any((/do_bm,do_bmmass,do_bmomp/)) ) then 
        used = send_data ( id_tref, t_ref, Time, is, js, 1, rmask=mask )
        used = send_data ( id_qref, q_ref, Time, is, js, 1, rmask=mask )
        used = send_data ( id_bmflag, bmflag, Time, is, js)
        used = send_data ( id_klzbs, klzbs, Time, is, js)
      endif

      if ( do_bm ) then
        used = send_data (id_invtaubmt, invtaubmt, Time, is, js)
        used = send_data (id_invtaubmq, invtaubmq, Time, is, js)
      end if

      if (do_bmmass) then
        used = send_data ( id_massflux, massflux, Time, is, js, 1, rmask=mask)
      end if 

      used = send_data (id_ras_precip, rain_ras + snow_ras, Time, is, js)
      used = send_data (id_don_precip, rain_don + snow_don + &
                           rain_donmca + snow_donmca, Time, is, js)
      used = send_data (id_scale_donner, scale_donner, Time, is, js )
      used = send_data (id_scale_uw, scale_uw, Time, is, js )
      used = send_data (id_uw_precip, rain_uw + snow_uw, Time, is, js)
        
      if (id_ras_freq > 0) then
        tmplmask = rain_ras > 0. .or. snow_ras > 0.0
        where (tmplmask) 
          freq_count = 1.
        elsewhere
          freq_count = 0.
        end where
        used = send_data (id_ras_freq, freq_count, Time, is, js)
      endif
        
      if (id_don_freq > 0) then
        tmplmask = rain_don > 0. .or. snow_don > 0.0 .or. &
                   rain_donmca > 0. .or. snow_donmca > 0.0
        where (tmplmask) 
          freq_count = 1.
        elsewhere
          freq_count = 0.
        end where
        used = send_data (id_don_freq, freq_count, Time, is, js)
      endif

      if (id_enth_uw_col > 0) then
        tempdiag(:,:) = -HLV*rain_uw -HLS*snow_uw
        do k=1,kx
          tempdiag(:,:)  &
              = tempdiag(:,:)  &
              + ( CP_AIR*ttnd_uw(:,:,k)  &
                  -HLV*qltnd_uw(:,:,k)  &
                  -HLS*qitnd_uw(:,:,k)  &
                )*pmass(:,:,k)
        end do
        used = send_data (id_enth_uw_col, tempdiag, Time, is, js)
      endif

      if (id_wat_uw_col > 0) then
        tempdiag(:,:) = rain_uw + snow_uw
        do k=1,kx
          tempdiag(:,:)  &
            = tempdiag(:,:)  &
            + ( qtnd_uw(:,:,k) + qltnd_uw(:,:,k) + qitnd_uw(:,:,k)  &
              )*pmass(:,:,k)
        end do
        used = send_data (id_wat_uw_col, tempdiag, Time, is, js)
      endif

        
      if (id_uw_freq > 0) then
        tmplmask = rain_uw > 0. .or. snow_uw > 0.0
        where (tmplmask) 
          freq_count = 1.
        elsewhere
          freq_count = 0.
        end where
        used = send_data (id_uw_freq, freq_count, Time, is, js)
      endif

      used = send_data (id_tdt_uw, ttnd_uw, Time, is, js, 1, rmask=mask)
      used = send_data (id_qdt_uw, qtnd_uw, Time, is, js, 1, rmask=mask)
      used = send_data (id_qadt_uw, qatnd_uw, Time, is, js, 1, rmask=mask)
      used = send_data (id_qldt_uw, qltnd_uw, Time, is, js, 1, rmask=mask)
      used = send_data (id_qidt_uw, qitnd_uw, Time, is, js, 1, rmask=mask)
      used = send_data (id_qndt_uw, qntnd_uw, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    temperature change due to dry and moist convection:
!---------------------------------------------------------------------
      used = send_data (id_tdt_conv, ttnd_conv, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    vapor specific humidity change due to convection:
!---------------------------------------------------------------------
      used = send_data (id_qdt_conv, qtnd_conv, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    total precipitation due to convection:
!---------------------------------------------------------------------
        used = send_data (id_prec_conv, precip, Time, is, js)

!---------------------------------------------------------------------
!    frozen precipitation (snow) due to convection:
!---------------------------------------------------------------------
        used = send_data (id_snow_conv, fprec, Time, is, js)

!---------------------------------------------------------------------
!------- diagnostics for 3D precip_conv -------
!---------------------------------------------------------------------
        used = send_data ( id_conv_rain3d, conv_rain3d, Time, is, js, 1 )

!---------------------------------------------------------------------
!------- diagnostics for 3D snow_conv -------
!---------------------------------------------------------------------
        used = send_data ( id_conv_snow3d, conv_snow3d, Time, is, js, 1 )

!---------------------------------------------------------------------
!    surface wind gustiness due to convection:
!---------------------------------------------------------------------
        used = send_data (id_gust_conv, gust_cv, Time, is, js)

!---------------------------------------------------------------------
!    water vapor path tendency due to convection:
!---------------------------------------------------------------------
      if (id_q_conv_col > 0) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + qtnd_conv(:,:,k)*pmass(:,:,k)
        end do
        used = send_data (id_q_conv_col, tempdiag, Time, is, js)
      endif        
   
!---------------------------------------------------------------------
!    dry static energy tendency due to dry and moist convection:
!---------------------------------------------------------------------
      if (id_t_conv_col > 0) then
        tempdiag(:,:)=0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + ttnd_conv(:,:,k)*CP_AIR* &
                                                            pmass(:,:,k)
        end do
        used = send_data (id_t_conv_col, tempdiag, Time, is, js)
      endif        
   
!---------------------------------------------------------------------
!    cloud liquid, ice and area tendencies due to convection:
!---------------------------------------------------------------------
      if (do_strat) then

!---------------------------------------------------------------------
!    if cloud liquid diagnostics requested:
!---------------------------------------------------------------------
        if (id_qldt_conv > 0 .or. id_ql_conv_col > 0) then
          tempdiag1 = rdt(:,:,:,nql) - rdt_init(:,:,:,nql)

!---------------------------------------------------------------------
!    cloud liquid tendency due to convection:
!---------------------------------------------------------------------
           used = send_data (id_qldt_conv, tempdiag1, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    cloud liquid water path tendency due to convection:
!---------------------------------------------------------------------
          if (id_ql_conv_col > 0) then
            tempdiag(:,:)=0.
            do k=1,kx
              tempdiag(:,:) = tempdiag(:,:) + tempdiag1(:,:,k)* &
                                                           pmass(:,:,k)
            end do
            used = send_data (id_ql_conv_col, tempdiag, Time, is, js)
          endif        
        endif

!---------------------------------------------------------------------
!    if cloud drop diagnostics requested:
!---------------------------------------------------------------------
        if (id_qndt_conv > 0 .or. id_qn_conv_col > 0) then
          tempdiag1 = rdt(:,:,:,nqn) - rdt_init(:,:,:,nqn)

!---------------------------------------------------------------------
!    cloud drop tendency due to convection:
!---------------------------------------------------------------------
          used = send_data (id_qndt_conv, tempdiag1, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    cloud drop water path tendency due to convection:
!---------------------------------------------------------------------
          if (id_qn_conv_col > 0) then
            tempdiag(:,:)=0.
            do k=1,kx
              tempdiag(:,:) = tempdiag(:,:) + tempdiag1(:,:,k)* &
                                                           pmass(:,:,k)
            end do
            used = send_data (id_qn_conv_col, tempdiag, Time, is, js)
          endif        
        endif

!---------------------------------------------------------------------
!    if cloud ice diagnostics requested:
!---------------------------------------------------------------------
        if (id_qidt_conv > 0 .or. id_qi_conv_col > 0) then
          tempdiag1 = rdt(:,:,:,nqi) - rdt_init(:,:,:,nqi)

!---------------------------------------------------------------------
!    cloud ice tendency due to convection:
!---------------------------------------------------------------------
          used = send_data (id_qidt_conv, tempdiag1, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    cloud ice water path tendency due to convection:
!---------------------------------------------------------------------
          if (id_qi_conv_col > 0) then
            tempdiag(:,:)=0.
            do k=1,kx
              tempdiag(:,:) = tempdiag(:,:) + tempdiag1(:,:,k)*  &
                                                            pmass(:,:,k)
            end do
            used = send_data (id_qi_conv_col, tempdiag, Time, is, js)
          endif        
        endif        

!---------------------------------------------------------------------
!    if cloud area diagnostics requested:
!---------------------------------------------------------------------
        if (id_qadt_conv > 0 .or.  id_qa_conv_col > 0 ) then
          tempdiag1 = rdt(:,:,:,nqa) - rdt_init(:,:,:,nqa)

!---------------------------------------------------------------------
!    cloud area tendency due to convection:
!---------------------------------------------------------------------
          used = send_data (id_qadt_conv, tempdiag1, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    column integrated cloud mass tendency due to convection:
!---------------------------------------------------------------------
          if ( id_qa_conv_col > 0 ) then
            tempdiag(:,:)=0.
            do k=1,kx
              tempdiag(:,:) = tempdiag(:,:) + tempdiag1(:,:,k)*  &
                                                            pmass(:,:,k)
            end do
            used = send_data (id_qa_conv_col, tempdiag, Time, is, js)
          endif        
        endif
      endif
         
!---------------------------------------------------------------------
!    column integrated enthalpy and total water tendencies due to 
!    convection parameterization:
!---------------------------------------------------------------------
      if (id_enth_conv_col > 0) then
        tempdiag(:,:) = -HLV*precip -HLF*fprec
        do k=1,kx
          tempdiag(:,:)  &
           = tempdiag(:,:)  &
            + ( CP_AIR*ttnd_conv(:,:,k)  &
                -HLV*(rdt(:,:,k,nql) - rdt_init(:,:,k,nql))  &
                -HLS*(rdt(:,:,k,nqi) - rdt_init(:,:,k,nqi))  &
              )*pmass(:,:,k)
        end do
        used = send_data (id_enth_conv_col, tempdiag, Time, is, js)
      endif

      if (id_wat_conv_col > 0) then
        tempdiag(:,:) = precip
        do k=1,kx
          tempdiag(:,:)  &
           = tempdiag(:,:)  &
           + ( qtnd_conv(:,:,k)  &
              + rdt(:,:,k,nql) - rdt_init(:,:,k,nql)  &
              + rdt(:,:,k,nqi) - rdt_init(:,:,k,nqi)  &
            )*pmass(:,:,k)
        end do
        used = send_data (id_wat_conv_col, tempdiag, Time, is, js)
      endif

!---------------------------------------------------------------------
!    tracer tendencies due to convection:
!---------------------------------------------------------------------
      do n=1,size(rdt,4)
        if (tracers_in_donner(n) .or. &
            tracers_in_ras(n)    .or.  &
            tracers_in_mca(n)    .or.  &
            tracers_in_uw(n))    then
          if (id_tracerdt_conv(n) > 0 .or. id_tracerdt_conv_col(n) > 0) then
            tempdiag1 = rdt(:,:,:,n) - rdt_init(:,:,:,n)
            used = send_data (id_tracerdt_conv(n), tempdiag1, &
                              Time, is, js, 1, rmask=mask )

!---------------------------------------------------------------------
!    tracer column tendencies due to convection:
!---------------------------------------------------------------------
            if (id_tracerdt_conv_col(n) > 0) then
              tempdiag(:,:)=0.
              do k=1,kx
                tempdiag(:,:) = tempdiag(:,:) + tempdiag1(:,:,k)*  &
                                                            pmass(:,:,k)
              end do
              used = send_data (id_tracerdt_conv_col(n), tempdiag,  &
                                Time, is, js)
            endif        
          endif        
        endif
      end do

!---------------------------------------------------------------------
!    end the timing of the convection code section.
!---------------------------------------------------------------------
call mpp_clock_end   (excl_conv_diag)
      call mpp_clock_end (convection_clock)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!              LARGE-SCALE CONDENSATION PARAMETERIZATIONS
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!---------------------------------------------------------------------
!    begin the timing of the large-scale condensation code section.
!---------------------------------------------------------------------
      call mpp_clock_begin (largescale_clock)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!         A. NON-PROGNOSTIC CONDENSATION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!-----------------------------------------------------------------------
!    if a non-prognostic cloud scheme is active, then call lscale_cond 
!    to calculate the temperature and specific humidity tendencies 
!    related to the latent heat release associated with the large-scale 
!    supersaturation.
!-----------------------------------------------------------------------
      if (do_lsc) then
        call mpp_clock_begin (lscalecond_clock)
        call lscale_cond (tin, qin, pfull, phalf, coldT, rain, snow,  &
                          ttnd, qtnd, mask=mask)
call mpp_clock_begin (excl_lsc)
        call mpp_clock_end (lscalecond_clock)

!-----------------------------------------------------------------------
!    add the temperature and specific humidity increments to the updated
!    temperature and specific humidity fields (tin, qin). convert these
!    increments and the precipitation increments to rates and add to 
!    the arrays accumulating the total rates for all physical processes
!    (tdt, qdt, lprec, fprec).
!-----------------------------------------------------------------------
        do k=1,kx
         do j=1,jx
          do i=1,ix
            tin(i,j,k) = tin(i,j,k) + ttnd(i,j,k) 
            qin(i,j,k) = qin(i,j,k) + qtnd(i,j,k)
            tdt(i,j,k) = tdt(i,j,k) + ttnd(i,j,k)*dtinv 
            qdt(i,j,k) = qdt(i,j,k) + qtnd(i,j,k)*dtinv
          enddo
         enddo
        enddo
        do j=1,jx
         do i=1,ix
           lprec(i,j) = lprec(i,j) + rain(i,j)*dtinv
           fprec(i,j) = fprec(i,j) + snow(i,j)*dtinv
         enddo
        enddo

!--------------------------------------------------------------------
!    if rh_clouds is active, call rh_calc to determine the grid box
!    relative humidity. call rh_clouds_sum to pass this field to 
!    rh_clouds_mod so it may be used to determine the grid boxes which
!    will contain clouds for the radiation package.
!---------------------------------------------------------------------
call mpp_clock_end   (excl_lsc)
        if (do_rh_clouds) then
          call rh_calc (pfull, tin, qin, rh, mask)
          call rh_clouds_sum (is, js, rh)
        endif
call mpp_clock_begin (excl_lsc)

!--------------------------------------------------------------------
!    if the gordon diagnostic cloud parameterization is active, set a 
!    flag to indicate those grid points where drying has resulted from 
!    convective activity (cnvcntq). call rh_calc to determine the grid 
!    box relative humidity. call diag_cloud_sum to define the cloud 
!    field that will be seen by the radiation package.
!---------------------------------------------------------------------
        if (do_diag_clouds) then
          do k=1,kx
           do j=1,jx
            do i=1,ix
              if (qtnd_conv(i,j,k) < 0.0) then
                cnvcntq (i,j,k) = 1.0
              else
                cnvcntq (i,j,k) = 0.0
              endif
            enddo
           enddo
          enddo
call mpp_clock_end   (excl_lsc)
          call rh_calc (pfull, tin, qin, rh, mask)
          call diag_cloud_sum (is, js, tin, qin, rh, omega, qtnd,  &
                               cnvcntq, precip, kbot)
call mpp_clock_begin (excl_lsc)
        endif


call mpp_clock_end   (excl_lsc)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!       B. TIEDTKE / ROTSTAYN / KLEIN PROGNOSTIC CLOUD SCHEME  
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!--------------------------------------------------------------------
!    define the model-full-level mass flux profile, composed of cont-
!    ributions from the model omega field and from donner_deep_mod 
!    (if active).
!--------------------------------------------------------------------
      else if (do_strat) then
!        do k=1,kx
!          mc_full(:,:,k) = 0.5*(mc(:,:,k) + mc(:,:,k+1)) +   &
!                                                      mc_donner(:,:,k)
!        end do
call mpp_clock_begin (excl_strat)

!----------------------------------------------------------------------
!    define the grid box specific humidity and saturation specific 
!    humidity.
!------------------------------------------------------------------
 
!----------------------------------------------------------------------
!    define the grid box area whose humidity is affected by the 
!    convective clouds and the environmental fraction and environmental
!    rh.
!-------------------------------------------------------------------
     if (do_uw_conv .or. do_donner_deep) then

       call compute_qs (tin, pfull, qsat)

       temp_don1 = 0.0
       temp_don2 = 0.0
       temp_don3 = 0.0
       temp_uw_c = 0.0

       do k=1,kx
        do j=1,jx
         do i=1,ix
           qrf = MAX (qin(i,j,k), 0.0)

           if (do_donner_deep) then
             temp_don1 = donner_humidity_area(i,j,k) 
             temp_don2 = cell_cld_frac(i,j,k) + meso_cld_frac(i,j,k)
             temp_don3 = cell_cld_frac(i,j,k) + donner_humidity_factor(i,j,k)
           endif

           if (do_uw_conv) then
             temp_uw_c = shallow_cloud_area(i,j,k)
           endif

           convective_humidity_area(i,j,k) = temp_don1 + temp_uw_c
           environmental_fraction = 1.0 - temp_don2 - temp_uw_c
           environmental_qv = qrf - qsat(i,j,k)*(temp_don3 + temp_uw_c)

!---------------------------------------------------------------------
!    define the ratio of the grid-box relative humidity to the humidity
!    in the environment of the convective clouds.
!----------------------------------------------------------------------
 
!----------------------------------------------------------------------
!    grid box has vapor and there is vapor outside of the convective a
!    clouds available for condensation.
!----------------------------------------------------------------
           if (qrf /= 0.0 .and. environmental_qv > 0.0) then
 
!--------------------------------------------------------------------
!    there is grid box area not filled with convective clouds
!--------------------------------------------------------------------  
             if (environmental_fraction > 0.0) then
               convective_humidity_ratio(i,j,k) =    &
                        MAX (qrf*environmental_fraction/   &
                                 environmental_qv, 1.0)
 
!---------------------------------------------------------------------
!    grid box is filled with convective clouds.
!----------------------------------------------------------------------
             else
               convective_humidity_ratio(i,j,k) = -10.0
             endif

!--------------------------------------------------------------------
!    either no vapor or all vapor taken up in convective clouds so 
!    none left for large-scale cd.
!---------------------------------------------------------------------
           else
             convective_humidity_ratio(i,j,k) = 1.0
           endif
         end do
        end do
       end do
     else
       convective_humidity_area = 0.0
       convective_humidity_ratio = 1.0
     endif
        
!-----------------------------------------------------------------------
!    call strat_cloud to integrate the prognostic cloud equations. 
!-----------------------------------------------------------------------
call mpp_clock_end   (excl_strat)
        call mpp_clock_begin (stratcloud_clock)
        if (do_liq_num) then 
          call strat_cloud (Time, is, ie, js, je, dt, pfull, phalf,    & 
                         radturbten, tin, qin, tracer(:,:,:,nql), &
                         tracer(:,:,:,nqi), tracer(:,:,:,nqa),  &
                         omega, mc_full, diff_t, land, ttnd, qtnd,  &
                         qltnd, qitnd, qatnd, lscale_rain3d, lscale_snow3d, &
                         rain, snow, convective_humidity_ratio, convective_humidity_area, &
                         limit_conv_cloud_frac, &
                         mask=mask, qn=tracer(:,:,:,nqn), Aerosol=Aerosol, SN=qntnd)
        else 
          if ( do_lin_cld_microphys ) then
              do k=1,kx
               do j=1,jx
                do i=1,ix
                  delp(i,j,k) =  phalf(i,j,k+1)-phalf(i,j,k)
                  delz(i,j,k) = (zhalf(i,j,k+1)-zhalf(i,j,k))*tin(i,j,k)/tm(i,j,k)
                  qrtnd(i,j,k) = 0.
                  qstnd(i,j,k) = 0.
                  qgtnd(i,j,k) = 0.
                enddo
               enddo
              enddo

              call lin_cld_microphys_driver(qin, tracer(:,:,:,nql), tracer(:,:,:,nqr), &
                       tracer(:,:,:,nqi), tracer(:,:,:,nqs), tracer(:,:,:,nqg), &
                       tracer(:,:,:,nqa),                                       &
                       qtnd, qltnd, qrtnd, qitnd, qstnd, qgtnd, qatnd,          &
                       ttnd, tin,   pfull, delz, delp, area,                    &
                       dt, land, rain, snow, ice_lin, graupel_lin,      &
                       is, ie, js, je, 1, kx, ktop, kx, Time)

              do k=1,kx
               do j=1,jx
                do i=1,ix
                  if (k.eq.kx) then
! Add all "solid" form of precipitation into surf_snow
                    snow(i,j) = (snow(i,j) + ice_lin(i,j) + graupel_lin(i,j)) * dt/86400.
                    rain(i,j) =  rain(i,j) * dt/86400.
                  endif

! Update tendencies:
                  rdt(i,j,k,nqr) = rdt(i,j,k,nqr) + qrtnd(i,j,k)
                  rdt(i,j,k,nqs) = rdt(i,j,k,nqs) + qstnd(i,j,k)
                  rdt(i,j,k,nqg) = rdt(i,j,k,nqg) + qgtnd(i,j,k)

                   ttnd(i,j,k) =  ttnd(i,j,k) * dt
                   qtnd(i,j,k) =  qtnd(i,j,k) * dt
                  qltnd(i,j,k) = qltnd(i,j,k) * dt
                  qrtnd(i,j,k) = qrtnd(i,j,k) * dt
                  qitnd(i,j,k) = qitnd(i,j,k) * dt
                  qstnd(i,j,k) = qstnd(i,j,k) * dt
                  qgtnd(i,j,k) = qgtnd(i,j,k) * dt
                  qatnd(i,j,k) = qatnd(i,j,k) * dt

! Update rain_wat, snow_wat, graupel_wat
                  tracer(i,j,k,nqr) = tracer(i,j,k,nqr) + qrtnd(i,j,k)
                  tracer(i,j,k,nqs) = tracer(i,j,k,nqs) + qstnd(i,j,k)
                  tracer(i,j,k,nqg) = tracer(i,j,k,nqg) + qgtnd(i,j,k)
                enddo
               enddo
              enddo
          else
              call strat_cloud (Time, is, ie, js, je, dt, pfull, phalf,    & 
                                radturbten, tin, qin, tracer(:,:,:,nql), &
                                tracer(:,:,:,nqi), tracer(:,:,:,nqa),  &
                                omega, mc_full, diff_t, land, ttnd, qtnd,  &
                                qltnd, qitnd, qatnd, lscale_rain3d, lscale_snow3d, &
                                rain, snow,convective_humidity_ratio, &
                                convective_humidity_area, limit_conv_cloud_frac, mask=mask)
          endif
        endif
        call mpp_clock_end (stratcloud_clock)
call mpp_clock_begin (excl_strat)
    
!----------------------------------------------------------------------
!    upon return from strat_cloud, update the cloud liquid, ice and area.
!    update the temperature and specific humidity fields.
!----------------------------------------------------------------------
        do k=1,kx
         do j=1,jx
          do i=1,ix
            tracer(i,j,k,nql) = tracer(i,j,k,nql) + qltnd(i,j,k)
            tracer(i,j,k,nqi) = tracer(i,j,k,nqi) + qitnd(i,j,k)
            tracer(i,j,k,nqa) = tracer(i,j,k,nqa) + qatnd(i,j,k)

!   save the lsc fields for use in radiation package.
            lsc_liquid(i,j,k) = tracer(i,j,k,nql)
            lsc_ice(i,j,k) = tracer(i,j,k,nqi)
            lsc_cloud_area(i,j,k) = tracer(i,j,k,nqa)

            if (do_liq_num) then 
              tracer(i,j,k,nqn) = tracer(i,j,k,nqn) + qntnd(i,j,k)
!   save the lsc fields for use in radiation package.
              lsc_droplet_number(i,j,k) = tracer(i,j,k,nqn)
            endif

            tin(i,j,k) = tin(i,j,k) + ttnd(i,j,k)
            qin(i,j,k) = qin(i,j,k) + qtnd(i,j,k)

          enddo
         enddo
        enddo

        
        used = send_data (id_qvout, qin, Time, is, js, 1, rmask=mask)
        used = send_data (id_qaout, tracer(:,:,:,nqa), Time, is, js, 1, rmask=mask)
        used = send_data (id_qlout, tracer(:,:,:,nql), Time, is, js, 1, rmask=mask)
        used = send_data (id_qiout, tracer(:,:,:,nqi), Time, is, js, 1, rmask=mask)

!----------------------------------------------------------------------
!    call strat_cloud_sum to make the cloud variables available for 
!    access by the radiation package. NOTE: this is no longer necessary,
!    and can be judiciously removed (provided other affiliated code 
!    and options are nullified).
!----------------------------------------------------------------------
call mpp_clock_end   (excl_strat)
        call strat_cloud_sum (is, js, tracer(:,:,:,nql),  &
                              tracer(:,:,:,nqi), tracer(:,:,:,nqa))
call mpp_clock_begin (excl_strat)

!----------------------------------------------------------------------
!    convert increments to tendencies.
!----------------------------------------------------------------------
        do k=1,kx
         do j=1,jx
          do i=1,ix
            ttnd(i,j,k) = ttnd(i,j,k)*dtinv 
            tdt(i,j,k) = tdt(i,j,k) + ttnd(i,j,k) 
            qtnd(i,j,k) = qtnd(i,j,k)*dtinv
            qdt(i,j,k) = qdt(i,j,k) + qtnd(i,j,k)

            qltnd(i,j,k) = qltnd(i,j,k)*dtinv
            rdt(i,j,k,nql) = rdt(i,j,k,nql) + qltnd(i,j,k)

            qitnd(i,j,k) = qitnd(i,j,k)*dtinv
            rdt(i,j,k,nqi) = rdt(i,j,k,nqi) + qitnd(i,j,k)

            qatnd(i,j,k) = qatnd(i,j,k)*dtinv
            rdt(i,j,k,nqa) = rdt(i,j,k,nqa) + qatnd(i,j,k)

            if (do_liq_num) then
              qntnd(i,j,k) = qntnd(i,j,k)*dtinv
              rdt(i,j,k,nqn) = rdt(i,j,k,nqn) + qntnd(i,j,k)
            endif

          enddo
         enddo
        enddo
        do j=1,jx
         do i=1,ix
           rain(i,j) = rain(i,j)*dtinv 
           snow(i,j) = snow(i,j)*dtinv
           lprec(i,j) = lprec(i,j) + rain(i,j)
           fprec(i,j) = fprec(i,j) + snow(i,j)
         enddo
        enddo

!----------------------------------------------------------------------
!    update the total tendency terms (temperature, vapor specific 
!    humidity, cloud liquid, cloud ice, cloud area, liquid precip,
!    frozen precip) with the contributions from the strat_cloud scheme.
!----------------------------------------------------------------------
   
call mpp_clock_end   (excl_strat)
      endif  ! (do_lsc)

call mpp_clock_begin (excl_wet_dep)
!---------------------------------------------------------------------
!    calculate the wet deposition associated with the large scale 
!    condensation. 
!---------------------------------------------------------------------
      if (do_strat) then
        do k=1,kx
         do j=1,jx
          do i=1,ix
            qtnd_wet(i,j,k) = qtnd(i,j,k) + qltnd(i,j,k) + qitnd(i,j,k)
! Count precipitation formed over timestep plus cloud amount at end of timestep
            if (do_lin_cld_microphys) then
              temp_1 = tracer(i,j,k,nqr) + tracer(i,j,k,nqs) + tracer(i,j,k,nqg)
            else
              temp_1 = (lscale_rain3d(i,j,k+1) - lscale_rain3d(i,j,k) &
                     +  lscale_snow3d(i,j,k+1) - lscale_snow3d(i,j,k)) &
                     *  dt / pmass(i,j,k)       ! convert from kg/m2/s to kg/kg
            endif
            cloud_wet(i,j,k) = temp_1 + tracer(i,j,k,nql) + tracer(i,j,k,nqi)
            cloud_frac(i,j,k) = max( min( tracer(i,j,k,nqa), 1. ), 0. )
          enddo
         enddo
        enddo
      else
         qtnd_wet = qtnd
         cloud_wet = 0.5e-3
         cloud_frac = 1.
      end if
      do n=1,size(rdt,4)
        if ( n /= nsphum ) then
          if ( .not. do_strat .or. (n /= nql .and. n /= nqi .and. n /= nqa .and. n /= nqn) ) then
            wetdeptnd = 0.0
call mpp_clock_end   (excl_wet_dep)
            call wet_deposition( n, t, pfull, phalf, zfull, zhalf, rain, snow, &
                                 qtnd_wet, cloud_wet, cloud_frac, &
                                 lscale_rain3d, lscale_snow3d, &
                                 tracer(:,:,:,n), wetdeptnd, &
                                 Time, 'lscale', is, js, dt )
call mpp_clock_begin (excl_wet_dep)
            do k=1,kx
             do j=1,jx
              do i=1,ix
                rdt (i,j,k,n) = rdt(i,j,k,n) - wetdeptnd(i,j,k)
                wet_data(i,j,k,n) = wet_data(i,j,k,n) + wetdeptnd(i,j,k)
              enddo
             enddo
            enddo

            used = send_data( id_wet_deposition(n), wet_data(:,:,:,n), &
                              Time,is_in=is,js_in=js )
          end if
        end if
      end do
call mpp_clock_end   (excl_wet_dep)
call mpp_clock_begin (excl_lsc_diag)


!---------------------------------------------------------------------
!    output diagnostics associated with the large-scale condensation
!    scheme.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    temperature change due to large-scale condensation:
!---------------------------------------------------------------------
      used = send_data (id_tdt_ls, ttnd, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    specific humidity change due to large-scale condensation:
!---------------------------------------------------------------------
      used = send_data (id_qdt_ls, qtnd, Time, is, js, 1, rmask=mask)
      used = send_data (id_lsc_precip, rain + snow, Time, is, js)
        
      if (id_lsc_freq > 0) then
        tmplmask = rain > 0. .or. snow > 0.0
        where (tmplmask) 
          freq_count = 1.
        elsewhere
          freq_count = 0.
        end where
        used = send_data (id_lsc_freq, freq_count, Time, is, js)
      endif

!---------------------------------------------------------------------
!    total precipitation rate due to large-scale condensation:
!---------------------------------------------------------------------
      used = send_data (id_prec_ls, rain+snow, Time, is, js)

!---------------------------------------------------------------------
!    snowfall rate due to large-scale condensation:
!---------------------------------------------------------------------
      used = send_data (id_snow_ls, snow, Time, is, js)

!---------------------------------------------------------------------
!    water vapor path tendency due to large-scale condensation:
!---------------------------------------------------------------------
      if (id_q_ls_col > 0) then
        tempdiag(:,:) = 0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + qtnd(:,:,k)*pmass(:,:,k)
        end do
        used = send_data (id_q_ls_col, tempdiag, Time, is, js)
      end if        
   
!---------------------------------------------------------------------
!    dry static energy tendency due to large-scale condensation:
!---------------------------------------------------------------------
      if (id_t_ls_col > 0) then
        tempdiag(:,:) = 0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + ttnd(:,:,k)*  &
                                                    CP_AIR*pmass(:,:,k)
        end do
        used = send_data (id_t_ls_col, tempdiag, Time, is, js)
      end if        

!---------------------------------------------------------------------
!    define diagnostics specific to the strat_cloud formulation:
!---------------------------------------------------------------------
      if (do_strat) then

!---------------------------------------------------------------------
!    total cumulus mass flux due to strat_cloud parameterization:
!---------------------------------------------------------------------
        used = send_data (id_mc_full, mc_full, Time, is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    cloud liquid, ice and area tendencies due to strat_cloud 
!    parameterization:
!---------------------------------------------------------------------
        used = send_data (id_qldt_ls, qltnd, Time, is, js, 1, rmask=mask)
        used = send_data (id_qidt_ls, qitnd, Time, is, js, 1, rmask=mask)
        used = send_data (id_qadt_ls, qatnd, Time, is, js, 1, rmask=mask)
        if (do_liq_num) used = send_data (id_qndt_ls, qntnd, Time, &
                                          is, js, 1, rmask=mask)

!---------------------------------------------------------------------
!    cloud liquid and ice water path tendencies due to strat_cloud 
!    parameterization:
!---------------------------------------------------------------------
        if (id_ql_ls_col > 0) then
          tempdiag(:,:) = 0.
          do k=1,kx
            tempdiag(:,:) = tempdiag(:,:) + qltnd(:,:,k)*pmass(:,:,k)
          end do
          used = send_data (id_ql_ls_col, tempdiag, Time, is, js)
        endif        
        if (do_liq_num .and. id_qn_ls_col > 0) then
          tempdiag(:,:) = 0.
          do k=1,kx
            tempdiag(:,:) = tempdiag(:,:) + qntnd(:,:,k)*pmass(:,:,k)
          end do
          used = send_data (id_qn_ls_col, tempdiag, Time, is, js)
        endif
        if (id_qi_ls_col > 0) then
          tempdiag(:,:) = 0.
          do k=1,kx
            tempdiag(:,:) = tempdiag(:,:) + qitnd(:,:,k)*pmass(:,:,k)
          end do
          used = send_data (id_qi_ls_col, tempdiag, Time, is, js)
        endif        
      
!---------------------------------------------------------------------
!    column integrated enthalpy and total water tendencies due to 
!    strat_cloud  parameterization:
!---------------------------------------------------------------------
       if (id_enth_ls_col > 0) then
         tempdiag(:,:) = -HLV*rain -HLS*snow
         do k=1,kx
           tempdiag(:,:)  &
            = tempdiag(:,:)  &
            + ( CP_AIR*ttnd(:,:,k)  &
               -HLV*qltnd(:,:,k) -HLS*qitnd(:,:,k)  &
              )*pmass(:,:,k)
         end do
         used = send_data (id_enth_ls_col, tempdiag, Time, is, js)
       endif
 
      if (id_wat_ls_col > 0) then
         tempdiag(:,:) = rain+snow
         do k=1,kx
           tempdiag(:,:)  &
            = tempdiag(:,:)  &
             + ( qtnd(:,:,k)+qltnd(:,:,k)+qitnd(:,:,k) )*pmass(:,:,k)
         end do
         used = send_data (id_wat_ls_col, tempdiag, Time, is, js)
         endif

!---------------------------------------------------------------------
!    stratiform cloud volume tendency due to strat_cloud 
!    parameterization:
!---------------------------------------------------------------------
        if (id_qa_ls_col > 0) then
          tempdiag(:,:) = 0.
          do k=1,kx
            tempdiag(:,:) = tempdiag(:,:) + qatnd(:,:,k)*pmass(:,:,k)
          end do
          used = send_data (id_qa_ls_col, tempdiag, Time, is, js)
        endif        

!---------------------------------------------------------------------
!---- diagnostics for large scale precip -----------
!---------------------------------------------------------------------
        used = send_data(id_lscale_rain3d, lscale_rain3d, Time, is, js, 1)

!---------------------------------------------------------------------
!---- diagnostics for large scale snow -------------
!---------------------------------------------------------------------
        used = send_data(id_lscale_snow3d, lscale_snow3d, Time, is, js, 1)

      endif ! (do_strat)

!---------------------------------------------------------------------
!    end the timing of the large-scale condensation code section.
!---------------------------------------------------------------------
      call mpp_clock_end (largescale_clock)
 
call mpp_clock_end   (excl_lsc_diag)
call mpp_clock_begin (excl_gen_diag)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                  GENERAL MOISTURE DIAGNOSTICS 
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!--------------------------------------------------------------------
!    output diagnostics obtained from the combination of convective and
!    large-scale parameterizations.  
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    total precipitation (all sources):
!---------------------------------------------------------------------
      precip = fprec + lprec
      used = send_data (id_precip, precip, Time, is, js)

!---------------------------------------------------------------------
!    column integrated enthalpy and total water tendencies due to 
!    moist processes:
!---------------------------------------------------------------------
      if (id_enth_moist_col > 0 .or. id_max_enthalpy_imbal > 0) then
        tempdiag(:,:) = -HLV*precip -HLF*fprec
        do k=1,kx
          tempdiag(:,:)  = tempdiag(:,:)  +   &
             (CP_AIR*(tdt(:,:,k) - tdt_init (:,:,k))     - &
               (HLV*(rdt(:,:,k,nql) - rdt_init(:,:,k,nql)) +  &
                HLS*(rdt(:,:,k,nqi) - rdt_init(:,:,k,nqi))))*  &
                                                          pmass(:,:,k)
        end do
      endif
      used = send_data (id_enth_moist_col, tempdiag, Time, is, js)

      if (id_max_enthalpy_imbal > 0) then
        do j=js,je
         do i=is,ie
           if (abs(tempdiag   (i,j)) > max_enthalpy_imbal(i,j)) then
              max_enthalpy_imbal(i,j) = abs (tempdiag  (i,j))
           endif
         end do
        end do
        used = send_data(id_max_enthalpy_imbal, max_enthalpy_imbal, &
                        Time, is, js)
      endif
  
      if (id_wat_moist_col > 0 .or. id_max_water_imbal > 0) then
        tempdiag(:,:) = precip
        do k=1,kx
          tempdiag(:,:)  = tempdiag(:,:) +  &
              (qdt(:,:,k) - qdt_init(:,:,k)  +  &
               rdt(:,:,k,nql) - rdt_init(:,:,k,nql) +  &
               rdt(:,:,k,nqi) - rdt_init(:,:,k,nqi))*pmass(:,:,k)
        end do
      endif
      used = send_data (id_wat_moist_col, tempdiag, Time, is, js)

      if (id_max_water_imbal > 0) then
        do j=js,je
         do i=is,ie
           if (abs(tempdiag(i,j)) > max_water_imbal(i,j)) then
               max_water_imbal(i,j) = abs (tempdiag(i,j))
           endif
         end do
        end do
        used = send_data(id_max_water_imbal, max_water_imbal, Time, is, js)
      endif

!---------------------------------------------------------------------
!    water vapor, liquid water and ice water column paths:
!---------------------------------------------------------------------
      if (id_WVP > 0) then
        wvp(:,:) = 0.
        do k=1,kx
          wvp(:,:) = wvp(:,:) + qin(:,:,k)*pmass(:,:,k)
        end do
        used = send_data (id_WVP, wvp, Time, is, js)
      endif
      if (id_LWP > 0 .and. do_strat) then
        lwp(:,:) = 0.
        do k=1,kx
          lwp(:,:) = lwp(:,:) + tracer(:,:,k,nql)*pmass(:,:,k)
        end do
        used = send_data (id_LWP, lwp, Time, is, js)
      endif
      if (id_IWP > 0 .and. do_strat) then
        iwp(:,:) = 0.
        do k=1,kx
          iwp(:,:) = iwp(:,:) + tracer(:,:,k,nqi)*pmass(:,:,k)
        end do
        used = send_data (id_IWP, iwp, Time, is, js)
      endif

!---------------------------------------------------------------------
!    column integrated cloud mass:
!---------------------------------------------------------------------
      if (id_AWP > 0 .and. do_strat) then
        tempdiag(:,:) = 0.
        do k=1,kx
          tempdiag(:,:) = tempdiag(:,:) + tracer(:,:,k,nqa)*pmass(:,:,k)
        end do
        used = send_data (id_AWP, tempdiag, Time, is, js)
      endif

!---------------------------------------------------------------------
!    relative humidity:         
!---------------------------------------------------------------------
      if (id_rh > 0) then
        if (.not. (do_rh_clouds .or. do_diag_clouds)) then 
          call rh_calc (pfull, tin, qin, RH, mask)
        endif
        used = send_data (id_rh, rh*100., Time, is, js, 1, rmask=mask)
      endif

!---------------------------------------------------------------------
!    saturation specific humidity:         
!---------------------------------------------------------------------
      if (id_qs > 0) then
        call compute_qs (tin, pfull, qsat, q=qin)
        used = send_data (id_qs, qsat, Time, is, js, 1, rmask=mask)
      endif

!------- diagnostics for CAPE and CIN, 

!!-- compute and write out CAPE and CIN--
   if ( id_cape > 0 .or. id_cin > 0) then
!! calculate r
         rin = qin/(1.0 - qin) ! XXX rin is not mixing ratio when do_simple=.true.
         avgbl = .false.
         do j=js,je
            do i=is,ie
               call capecalcnew( kx, pfull(i,j,:), phalf(i,j,:), CP_AIR, RDGAS, RVGAS, &
                         HLV, KAPPA, tin(i,j,:), rin(i,j,:), avgbl, cape(i,j), cin(i,j))
            end do
         end do
        if (id_cape > 0) then 
             used = send_data ( id_cape, cape, Time, is, js )
        end if

        if ( id_cin > 0 ) then
             used = send_data ( id_cin, cin, Time, is, js )
        end if
   end if

!---------------------------------------------------------------------
!    output the global integral of precipitation in units of mm/day.
!---------------------------------------------------------------------
      call sum_diag_integral_field ('prec', precip*SECONDS_PER_DAY,  &
                                                                is, js)

!-----------------------------------------------------------------------

call mpp_clock_end   (excl_gen_diag)

end subroutine moist_processes




!#######################################################################

subroutine moist_processes_init ( id, jd, kd, lonb, latb, pref, &
                                  axes, Time, doing_donner, &
                                  doing_uw_conv)

!-----------------------------------------------------------------------
integer,              intent(in)  :: id, jd, kd, axes(4)
real, dimension(:,:), intent(in)  :: lonb, latb
real, dimension(:),   intent(in)  :: pref
type(time_type),      intent(in)  :: Time
logical,              intent(out) :: doing_donner, doing_uw_conv
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
character(len=80)  :: scheme
integer            :: secs, days
integer            :: k
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

         if ( do_mca .and. do_bm ) call error_mesg   &
                   ('moist_processes_init',  &
                    'both do_mca and do_bm cannot be specified', FATAL)
         if ( do_ras .and. do_bm ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bm and do_ras cannot be specified', FATAL)
         if ( do_bm .and. do_bmmass ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bm and do_bmmass cannot be specified', FATAL)
         if ( do_bm .and. do_bmomp ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bm and do_bmomp cannot be specified', FATAL)
         if ( do_bmomp .and. do_bmmass ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bmomp and do_bmmass cannot be specified', FATAL)
         if ( do_bmmass .and. do_mca ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bmmass and do_mca cannot be specified', FATAL)
         if ( do_bmmass .and. do_ras ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bmmass and do_ras cannot be specified', FATAL)
         if ( do_bmomp .and. do_mca ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bmomp and do_mca cannot be specified', FATAL)
         if ( do_bmomp .and. do_ras ) call error_mesg   &
                    ('moist_processes_init',  &
                     'both do_bmomp and do_ras cannot be specified', FATAL)

         if ( do_lsc .and. do_strat ) call error_mesg   &
                 ('moist_processes_init',  &
                  'both do_lsc and do_strat cannot be specified', FATAL)
         if (.not. do_lsc .and. .not. do_strat) then
           call error_mesg ('moist_processes_mod', &
              'must activate either do_lsc or do_strat in order to &
                             &include large-scale condensation', FATAL)
         endif

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
         
         if (force_donner_moist_conserv .and. &
               .not. do_donner_conservation_checks) then
           call error_mesg ('moist_processes', &
              'when force_donner_moist_conserv is .true., &
                &do_donner_conservation_checks must be .true.', FATAL)
         endif

         if (use_updated_profiles_for_uw .and.   &
             .not. (do_donner_before_uw) ) then
           call error_mesg ('moist_processes_init', &
            'use_updated_profiles_for_uw is only meaningful when &
                               &do_donner_before_uw is true', FATAL)
         endif

         if (only_one_conv_scheme_per_column .and.   &
             .not. (do_donner_before_uw) ) then
           call error_mesg ('moist_processes_init', &
            'only_one_conv_scheme_per_column is only meaningful when &
                               &do_donner_before_uw is true', FATAL)
         endif

         if (limit_conv_cloud_frac .and.   &
                 .not. do_donner_before_uw) then
           call error_mesg ('moist_processes', &
              'when limit_conv_cloud_frac is .true., &
                 &do_donner_before_uw must be .true.', FATAL)
         endif

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

      nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      if (nqn == NO_TRACER .and. do_liq_num ) &
        call error_mesg ('moist_processes', &
             'prognostic droplet number scheme requested but tracer not found', FATAL)

!------------ initialize various schemes ----------

      if (do_lsc) then
                     call lscale_cond_init ()
                     if (do_rh_clouds) call rh_clouds_init (id,jd,kd)
                     if (do_diag_clouds) call diag_cloud_init (id,jd,kd,ierr)
      endif

      if (do_strat)  call strat_cloud_init (axes,Time,id,jd,kd)
      if (do_dryadj) call dry_adj_init ()
      if (do_cmt)    call cu_mo_trans_init (axes,Time, doing_diffusive)
      if (do_bm)     call betts_miller_init () 

      if (do_cmt) then
        if ( .not. do_ras .and. .not. do_donner_deep  .and. &
             .not. do_uw_conv) then
          call error_mesg ( 'moist_processes_mod', &
                'do_cmt specified but no cumulus schemes activated', &
                                                              FATAL)
        endif
        if (trim(cmt_mass_flux_source) == 'ras') then
          cmt_uses_ras = .true.
          cmt_uses_donner = .false.
          cmt_uses_uw = .false.
          if (.not. do_ras) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_ras then ras_mod must be active', FATAL)
           endif

        else if (trim(cmt_mass_flux_source) == 'donner') then
          cmt_uses_ras = .false.
          cmt_uses_donner = .true.
          cmt_uses_uw = .false.
          if (.not. do_donner_deep)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_donner then donner_deep_mod must&
                                               & be active', FATAL)
           endif

        else if (trim(cmt_mass_flux_source) == 'uw') then
          cmt_uses_ras = .false.
          cmt_uses_donner = .false.
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_uw then uw_conv_mod must&
                                               & be active', FATAL)
           endif

        else if (trim(cmt_mass_flux_source) == 'donner_and_ras') then
          cmt_uses_ras = .true.
          if (.not. do_ras) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_ras then ras_mod must be active', FATAL)
           endif
          cmt_uses_donner = .true.
          if (.not. do_donner_deep)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_donner then donner_deep_mod must&
                                               & be active', FATAL)
           endif
          cmt_uses_uw = .false.

        else if (trim(cmt_mass_flux_source) == 'donner_and_uw') then
          cmt_uses_uw = .true.
          if (.not. do_uw_conv) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_uw then uw_conv_mod must be active', FATAL)
           endif
          cmt_uses_donner = .true.
          if (.not. do_donner_deep)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_donner then donner_deep_mod must&
                                               & be active', FATAL)
           endif
          cmt_uses_ras = .false.

        else if (trim(cmt_mass_flux_source) == 'ras_and_uw') then
          cmt_uses_ras = .true.
          if (.not. do_ras) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_ras then ras_mod must be active', FATAL)
           endif
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_uw then uw_conv_mod must&
                                               & be active', FATAL)
           endif
          cmt_uses_donner = .false.

        else if (trim(cmt_mass_flux_source) == 'donner_and_ras_and_uw') then
          cmt_uses_ras = .true.
          if (.not. do_ras) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_ras then ras_mod must be active', FATAL)
           endif
          cmt_uses_donner = .true.
          if (.not. do_donner_deep)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_donner then donner_deep_mod must&
                                               & be active', FATAL)
           endif
          cmt_uses_uw = .true.
          if (.not. do_uw_conv)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_uw then uw_conv_mod must&
                                               & be active', FATAL)
           endif
        else if (trim(cmt_mass_flux_source) == 'all') then
          if (do_ras) then
            cmt_uses_ras = .true.
          else
            cmt_uses_ras = .false.
          endif
          if (do_donner_deep)  then
            cmt_uses_donner = .true.
          else
            cmt_uses_donner = .false.
          endif
          if (do_uw_conv)  then
            cmt_uses_uw = .true.
          else
            cmt_uses_uw = .false.
          endif
        else
          call error_mesg ('moist_processes_mod', &
             'invalid specification of cmt_mass_flux_source', FATAL)
        endif

        if (cmt_uses_uw .and. .not. doing_diffusive) then
          call error_mesg ('moist_processes_mod', &
             'currently cannot do non-local cmt with uw as mass &
                                                &flux_source', FATAL)
        endif
          

      endif

  
!----- initialize quantities for global integral package -----

   call diag_integral_field_init ('prec', 'f6.3')
   

!----- initialize clocks -----

   excl_pre_conv        =  mpp_clock_id( ' M_Proc excl: pre_conv      ',grain=CLOCK_LOOP)
   excl_shallow         =  mpp_clock_id( ' M_Proc excl: shallow       ',grain=CLOCK_LOOP)
   excl_donner_deep     =  mpp_clock_id( ' M_Proc excl: donner_deep   ',grain=CLOCK_LOOP)
   excl_donner_sph      =  mpp_clock_id( ' M_Proc excl: donner_sph    ',grain=CLOCK_LOOP)
   excl_donner_moist_cons  =  mpp_clock_id( ' M_Proc excl: donner_moist_conserv',grain=CLOCK_LOOP)
   excl_donner_deep_2   =  mpp_clock_id( ' M_Proc excl: donner_deep_2 ',grain=CLOCK_LOOP)
   excl_moist_conv_adj  =  mpp_clock_id( ' M_Proc excl: moist_conv_adj',grain=CLOCK_LOOP)
   excl_ras             =  mpp_clock_id( ' M_Proc excl: ras           ',grain=CLOCK_LOOP)
   excl_cmt             =  mpp_clock_id( ' M_Proc excl: cmt           ',grain=CLOCK_LOOP)
   excl_tracer_dens     =  mpp_clock_id( ' M_Proc excl: tracer_dens   ',grain=CLOCK_LOOP)
   excl_uw_conv         =  mpp_clock_id( ' M_Proc excl: uw_conv       ',grain=CLOCK_LOOP)
   excl_conv_diag       =  mpp_clock_id( ' M_Proc excl: conv_diag     ',grain=CLOCK_LOOP)
   excl_lsc             =  mpp_clock_id( ' M_Proc excl: lsc           ',grain=CLOCK_LOOP)
   excl_strat           =  mpp_clock_id( ' M_Proc excl: strat         ',grain=CLOCK_LOOP)
   excl_wet_dep         =  mpp_clock_id( ' M_Proc excl: wet_dep       ',grain=CLOCK_LOOP)
   excl_lsc_diag        =  mpp_clock_id( ' M_Proc excl: lsc_diag      ',grain=CLOCK_LOOP)
   excl_gen_diag        =  mpp_clock_id( ' M_Proc excl: gen_diag      ',grain=CLOCK_LOOP)


   convection_clock =     &
       mpp_clock_id( '   Physics_up: Moist Proc: Conv', &
           grain=CLOCK_MODULE                         )
   largescale_clock =     &
       mpp_clock_id( '   Physics_up: Moist Proc: LS', &
           grain=CLOCK_MODULE                         )

!---------------------------------------------------------------------
!    define output variables indicating whether certain convection 
!    schemes have been activated.
!---------------------------------------------------------------------
     doing_donner = do_donner_deep
     doing_uw_conv = do_uw_conv

   donner_clock =     &
       mpp_clock_id( '   Moist Processes: Donner_deep', &
           grain=CLOCK_MODULE                         )

   mca_clock =     &
       mpp_clock_id( '   Moist Processes: MCA', &
           grain=CLOCK_MODULE                         )

   ras_clock =     &
       mpp_clock_id( '   Moist Processes: RAS', &
           grain=CLOCK_MODULE                         )

   closure_clock =     &
       mpp_clock_id( '   Moist Processes: conv_closure', &
           grain=CLOCK_MODULE                         )

   shallowcu_clock =     &
       mpp_clock_id( '   Moist Processes: Shallow_cu', &
           grain=CLOCK_MODULE                         )

   cmt_clock =     &
       mpp_clock_id( '   Moist Processes: CMT', &
           grain=CLOCK_MODULE                         )

   lscalecond_clock =     &
       mpp_clock_id( '   Moist Processes: lscale_cond', &
           grain=CLOCK_MODULE                         )

   stratcloud_clock =     &
       mpp_clock_id( '   Moist Processes: Strat_cloud', &
           grain=CLOCK_MODULE                         )

!---------------------------------------------------------------------
!    retrieve the number of registered tracers in order to determine 
!    which tracers are to be convectively transported.
!---------------------------------------------------------------------
      call get_number_tracers (MODEL_ATMOS, num_prog= num_tracers)
 
!---------------------------------------------------------------------
!    allocate logical arrays to indicate the tracers which are to be
!    transported by the various available convective schemes. 
!    initialize these arrays to .false..
!---------------------------------------------------------------------
      allocate (tracers_in_donner(num_tracers))
      allocate (tracers_in_mca(num_tracers))
      allocate (tracers_in_ras(num_tracers))
      allocate (tracers_in_uw(num_tracers))
      tracers_in_donner = .false.
      tracers_in_mca = .false.
      tracers_in_ras = .false.
      tracers_in_uw = .false.

!----------------------------------------------------------------------
!    for each tracer, determine if it is to be transported by convect-
!    ion, and the convection schemes that are to transport it. set a 
!    logical flag to .true. for each tracer that is to be transported by
!    each scheme and increment the count of tracers to be transported
!    by that scheme.
!----------------------------------------------------------------------
      do n=1, num_tracers
        if (query_method ('convection', MODEL_ATMOS, n, scheme)) then
          select case (scheme)
            case ("none")
            case ("donner")
               num_donner_tracers = num_donner_tracers + 1
               tracers_in_donner(n) = .true.
            case ("mca")
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
            case ("ras")
               num_ras_tracers = num_ras_tracers + 1
               tracers_in_ras(n) = .true.
            case ("uw")
               num_uw_tracers = num_uw_tracers + 1
               tracers_in_uw(n) = .true.
            case ("donner_and_ras")
               num_donner_tracers = num_donner_tracers + 1
               tracers_in_donner(n) = .true.
               num_ras_tracers = num_ras_tracers + 1
               tracers_in_ras(n) = .true.
            case ("donner_and_mca")
               num_donner_tracers = num_donner_tracers + 1
               tracers_in_donner(n) = .true.
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
            case ("mca_and_ras")
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
               num_ras_tracers = num_ras_tracers + 1
               tracers_in_ras(n) = .true.
            case ("all")
               num_donner_tracers = num_donner_tracers + 1
               tracers_in_donner(n) = .true.
               num_mca_tracers = num_mca_tracers + 1
               tracers_in_mca(n) = .true.
               num_ras_tracers = num_ras_tracers + 1
               tracers_in_ras(n) = .true.
               num_uw_tracers = num_uw_tracers + 1
               tracers_in_uw(n) = .true.
            case default  ! corresponds to "none"
          end select
        endif
      end do

!--------------------------------------------------------------------
!    set a logical indicating if any tracers are to be transported by
!    each of the available convection parameterizations.
!--------------------------------------------------------------------
      if (num_donner_tracers > 0) then
        do_tracers_in_donner = .true.
      else
        do_tracers_in_donner = .false.
      endif
      if (num_mca_tracers > 0) then
        do_tracers_in_mca = .true.
      else
        do_tracers_in_mca = .false.
      endif
      if (num_ras_tracers > 0) then
        do_tracers_in_ras = .true.
      else
        do_tracers_in_ras = .false.
      endif
      if (num_uw_tracers > 0) then
        do_tracers_in_uw = .true.
      else
        do_tracers_in_uw = .false.
      endif     
     
!---------------------------------------------------------------------
!    check for proper use of do_unified_convective_closure.
!---------------------------------------------------------------------
      if (do_unified_convective_closure) then
        call error_mesg ('moist_processes_init', &
         'do_unified_convective_closure is currently not allowed &
               & - see rsh', FATAL)
      endif
      if (do_unified_convective_closure) then
        if (.not. (do_donner_deep) .or. .not. (do_uw_conv)   &
            .or. do_ras .or. do_mca ) then
          call error_mesg ('moist_processes_init',  &
             'must have only donner_deep and uw shallow activated &
                &when do_unified_convective_closure is .true.', FATAL)
         endif
      endif
        
!---------------------------------------------------------------------
!    allocate and initialize arrays to hold maximum enthalpy and water
!    imbalances in each column.
!---------------------------------------------------------------------
      allocate (max_enthalpy_imbal (id, jd))
      allocate (max_water_imbal (id, jd))
      max_enthalpy_imbal = 0.
      max_water_imbal = 0.


!--------------------------------------------------------------------
!    initialize the convection scheme modules.
!--------------------------------------------------------------------
      if (do_donner_deep) then
        call get_time (Time, secs, days)
        call donner_deep_init (lonb, latb, pref, axes, secs, days,  &
                               tracers_in_donner,  &
                               do_donner_conservation_checks, &
                               do_unified_convective_closure, using_fms)
        if (do_donner_conservation_checks) then
          allocate (max_enthalpy_imbal_don (id, jd))
          allocate (max_water_imbal_don (id, jd))
          max_enthalpy_imbal_don = 0.
          max_water_imbal_don = 0.
        endif
      endif ! (do_donner_deep)
 
      if (do_ras)  then
        call ras_init (do_strat, do_liq_num, axes,Time, tracers_in_ras)
      endif

      if (do_uw_conv) call uw_conv_init (do_strat, axes, Time, kd, &
                                          tracers_in_uw)

      if (do_mca .or. do_donner_deep)  then
        call  moist_conv_init (axes,Time, tracers_in_mca)
      endif
  
 
!----- initialize quantities for diagnostics output -----
 
      call diag_field_init ( axes, Time )

      if (do_lin_cld_microphys) then
         if (.not. do_strat) call error_mesg ('moist_processes_init',  &
                    'must also activate do_strat when do_lin_cld_microphys is active', FATAL)
         if (do_liq_num) call error_mesg ('moist_processes_init',  &
                    'do_lin_cld_microphys cannot be active with prognostic droplet &
                   & scheme (do_liq_num)', FATAL)
         nqr = get_tracer_index (MODEL_ATMOS, 'rainwat')
         nqs = get_tracer_index (MODEL_ATMOS, 'snowwat')
         nqg = get_tracer_index (MODEL_ATMOS, 'graupel')
         call lin_cld_microphys_init (axes, Time)
         ktop = 1
         do k = 1, kd
            if (pref(k) > 10.E2) then
              ktop=k
              exit
            endif
         enddo
         if (mpp_pe() == mpp_root_pe()) &
               write(*,*) 'Top layer for lin_cld_microphys=', ktop, pref(ktop)
      endif
 
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
      if (do_uw_conv    ) call     uw_conv_end
      if (do_lin_cld_microphys) call lin_cld_microphys_end

      deallocate (max_water_imbal)
      deallocate (max_enthalpy_imbal)
      if (do_donner_conservation_checks) then
        deallocate (max_water_imbal_don)
        deallocate (max_enthalpy_imbal_don)
      endif

      module_is_initialized = .false.
      
!-----------------------------------------------------------------------

end subroutine moist_processes_end


!#######################################################################
! <SUBROUTINE NAME="moist_processes_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine moist_processes_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp

  if (do_strat)       call strat_cloud_restart(timestamp)
  if (do_diag_clouds) call diag_cloud_restart(timestamp)
  if (do_donner_deep) call donner_deep_restart(timestamp)

end subroutine moist_processes_restart
! </SUBROUTINE> NAME="moist_processes_restart"


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

        IMPLICIT NONE

        REAL, INTENT (IN),    DIMENSION(:,:,:) :: pfull,T,qv
        REAL, INTENT (OUT),   DIMENSION(:,:,:) :: RH
        REAL, INTENT (IN), OPTIONAL, DIMENSION(:,:,:) :: MASK
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: esat

        real, parameter :: d622 = RDGAS/RVGAS
        real, parameter :: d378 = 1.-d622

! because Betts-Miller uses a simplified scheme for calculating the relative humidity
        if (do_simple) then
          call lookup_es(T, esat)
          RH(:,:,:) = pfull(:,:,:)
          RH(:,:,:) = MAX(RH(:,:,:),esat(:,:,:))  !limit where pfull ~ esat
          RH(:,:,:)=qv(:,:,:)/(d622*esat(:,:,:)/RH(:,:,:))
        else
          call compute_qs (T, pfull, rh, q=qv)
          RH(:,:,:)=qv(:,:,:)/RH(:,:,:)
        endif
      
        !IF MASK is present set RH to zero
        IF (present(MASK)) RH(:,:,:)=MASK(:,:,:)*RH(:,:,:)

END SUBROUTINE rh_calc



!#######################################################################
!cape calculation.
                                                                                
subroutine capecalcnew(kx,p,phalf,cp,rdgas,rvgas,hlv,kappa,tin,rin,&
                                avgbl,cape,cin)
                                                                                
!
!    Input:
!
!    kx          number of levels
!    p           pressure (index 1 refers to TOA, index kx refers to surface)
!    phalf       pressure at half levels
!    cp          specific heat of dry air
!    rdgas       gas constant for dry air
!    rvgas       gas constant for water vapor (used in Clausius-Clapeyron,
!                not for virtual temperature effects, which are not considered)
!    hlv         latent heat of vaporization
!    kappa       the constant kappa
!    tin         temperature of the environment
!    rin         specific humidity of the environment
!    avgbl       if true, the parcel is averaged in theta and r up to its LCL
!
!    Output:
!    cape        Convective available potential energy
!    cin         Convective inhibition (if there's no LFC, then this is set
!                to zero)
!
!    Algorithm:
!    Start with surface parcel.
!    Calculate the lifting condensation level (uses an analytic formula and a
!       lookup table).
!    Average under the LCL if desired, if this is done, then a new LCL must
!       be calculated.
!    Calculate parcel ascent up to LZB.
!    Calculate CAPE and CIN.
      implicit none
      integer, intent(in)                    :: kx
      logical, intent(in)                    :: avgbl
      real, intent(in), dimension(:)         :: p, phalf, tin, rin
      real, intent(in)                       :: rdgas, rvgas, hlv, kappa, cp
      real, intent(out)                      :: cape, cin
                                                                                
      integer            :: k, klcl, klfc, klzb, klcl2
      logical            :: nocape
      real, dimension(kx)   :: theta, tp, rp
      real                  :: t0, r0, es, rs, theta0, pstar, value, tlcl, &
                               a, b, dtdlnp, d2tdlnp2, thetam, rm, tlcl2, &
                               plcl2, plcl, plzb
                                                                                
      pstar = 1.e5
                                                                                
      nocape = .true.
      cape = 0.
      cin = 0.
      plcl = 0.
      plzb = 0.
      klfc = 0
      klcl = 0
      klzb = 0
      tp(1:kx) = tin(1:kx)
      rp(1:kx) = rin(1:kx)
                                                                                
! start with surface parcel
      t0 = tin(kx)
      r0 = rin(kx)
! calculate the lifting condensation level by the following:
! are you saturated to begin with?
      call lookup_es(t0,es)
      rs = rdgas/rvgas*es/p(kx)
      if (r0.ge.rs) then
! if you're already saturated, set lcl to be the surface value.
         plcl = p(kx)
! the first level where you're completely saturated.
         klcl = kx
! saturate out to get the parcel temp and humidity at this level
! first order (in delta T) accurate expression for change in temp
         tp(kx) = t0 + (r0 - rs)/(cp/hlv + hlv*rs/rvgas/t0**2.)
         call lookup_es(tp(kx),es)
         rp(kx) = rdgas/rvgas*es/p(kx)
      else
! if not saturated to begin with, use the analytic expression to calculate the
! exact pressure and temperature where you?re saturated.
         theta0 = tin(kx)*(pstar/p(kx))**kappa
! the expression that we utilize is 
! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00) = log(es/T**(1/kappa))
! The right hand side of this is only a function of temperature, therefore
! this is put into a lookup table to solve for temperature.
         if (r0.gt.0.) then
            value = log(theta0**(-1/kappa)*r0*pstar*rvgas/rdgas) 
            call lcltabl(value,tlcl)
            plcl = pstar*(tlcl/theta0)**(1/kappa)
! just in case plcl is very high up
            if (plcl.lt.p(1)) then
               plcl = p(1)
               tlcl = theta0*(plcl/pstar)**kappa
               write (*,*) 'hi lcl'
            end if
            k = kx
         else
! if the parcel sp hum is zero or negative, set lcl to 2nd to top level
            plcl = p(2)
            tlcl = theta0*(plcl/pstar)**kappa
!            write (*,*) 'zero r0', r0
            do k=2,kx
               tp(k) = theta0*(p(k)/pstar)**kappa
               rp(k) = 0.
! this definition of CIN contains everything below the LCL
               cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
            end do
            go to 11
         end if
! calculate the parcel temperature (adiabatic ascent) below the LCL.
! the mixing ratio stays the same
         do while (p(k).gt.plcl)
            tp(k) = theta0*(p(k)/pstar)**kappa
            call lookup_es(tp(k),es)
            rp(k) = rdgas/rvgas*es/p(k)
! this definition of CIN contains everything below the LCL
            cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
            k = k-1
         end do
! first level where you're saturated at the level
         klcl = k
	 if (klcl.eq.1) klcl = 2
! do a saturated ascent to get the parcel temp at the LCL.
! use your 2nd order equation up to the pressure above.
! moist adaibat derivatives: (use the lcl values for temp, humid, and
! pressure)
         a = kappa*tlcl + hlv/cp*r0
         b = hlv**2.*r0/cp/rvgas/tlcl**2.
         dtdlnp = a/(1. + b)
! first order in p
!         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
! second order in p (RK2)
! first get temp halfway up
         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)/2.
         if ((tp(klcl).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(klcl),es)
         rp(klcl) = rdgas/rvgas*es/(p(klcl) + plcl)*2.
         a = kappa*tp(klcl) + hlv/cp*rp(klcl)
         b = hlv**2./cp/rvgas*rp(klcl)/tp(klcl)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl)
!         d2tdlnp2 = (kappa + b - 1. - b/tlcl*(hlv/rvgas/tlcl - &
!                   2.)*dtdlnp)/ (1. + b)*dtdlnp - hlv*r0/cp/ &
!                   (1. + b)
! second order in p
!         tp(klcl) = tlcl + dtdlnp*log(p(klcl)/plcl) + .5*d2tdlnp2*(log(&
!             p(klcl)/plcl))**2.
         if ((tp(klcl).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(klcl),es)
         rp(klcl) = rdgas/rvgas*es/p(klcl)
!         write (*,*) 'tp, rp klcl:kx, new', tp(klcl:kx), rp(klcl:kx)
! CAPE/CIN stuff
         if ((tp(klcl).lt.tin(klcl)).and.nocape) then
! if you're not yet buoyant, then add to the CIN and continue
            cin = cin + rdgas*(tin(klcl) - &
                 tp(klcl))*log(phalf(klcl+1)/phalf(klcl))
         else
! if you're buoyant, then add to cape
            cape = cape + rdgas*(tp(klcl) - &
                  tin(klcl))*log(phalf(klcl+1)/phalf(klcl))
! if it's the first time buoyant, then set the level of free convection to k
            if (nocape) then
               nocape = .false.
               klfc = klcl
            endif
         end if
      end if
! then, start at the LCL, and do moist adiabatic ascent by the first order
! scheme -- 2nd order as well
      do k=klcl-1,1,-1
         a = kappa*tp(k+1) + hlv/cp*rp(k+1)
         b = hlv**2./cp/rvgas*rp(k+1)/tp(k+1)**2.
         dtdlnp = a/(1. + b)
! first order in p
!         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
! second order in p (RK2)
! first get temp halfway up
         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))/2.
         if ((tp(k).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(k),es)
         rp(k) = rdgas/rvgas*es/(p(k) + p(k+1))*2.
         a = kappa*tp(k) + hlv/cp*rp(k)
         b = hlv**2./cp/rvgas*rp(k)/tp(k)**2.
         dtdlnp = a/(1. + b)
! second half of RK2
         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1))
!         d2tdlnp2 = (kappa + b - 1. - b/tp(k+1)*(hlv/rvgas/tp(k+1) - &
!               2.)*dtdlnp)/(1. + b)*dtdlnp - hlv/cp*rp(k+1)/(1. + b)
! second order in p

!         tp(k) = tp(k+1) + dtdlnp*log(p(k)/p(k+1)) + .5*d2tdlnp2*(log( &
!             p(k)/p(k+1)))**2.
! if you're below the lookup table value, just presume that there's no way
! you could have cape and call it quits
         if ((tp(k).lt.173.16).and.nocape) go to 11
         call lookup_es(tp(k),es)
         rp(k) = rdgas/rvgas*es/p(k)
         if ((tp(k).lt.tin(k)).and.nocape) then
! if you're not yet buoyant, then add to the CIN and continue
            cin = cin + rdgas*(tin(k) - tp(k))*log(phalf(k+1)/phalf(k))
         elseif((tp(k).lt.tin(k)).and.(.not.nocape)) then
! if you have CAPE, and it's your first time being negatively buoyant,
! then set the level of zero buoyancy to k+1, and stop the moist ascent
            klzb = k+1
            go to 11
         else
! if you're buoyant, then add to cape
            cape = cape + rdgas*(tp(k) - tin(k))*log(phalf(k+1)/phalf(k))
! if it's the first time buoyant, then set the level of free convection to k
            if (nocape) then
               nocape = .false.
               klfc = k
            endif
         end if
      end do
 11   if(nocape) then
! this is if you made it through without having a LZB
! set LZB to be the top level.
         plzb = p(1)
         klzb = 0
         klfc = 0
         cin = 0.
         tp(1:kx) = tin(1:kx)
         rp(1:kx) = rin(1:kx)
      end if
!      write (*,*) 'plcl, klcl, tlcl, r0 new', plcl, klcl, tlcl, r0
!      write (*,*) 'tp, rp new', tp, rp
!       write (*,*) 'tp, new', tp
!       write (*,*) 'tin new', tin
!       write (*,*) 'klcl, klfc, klzb new', klcl, klfc, klzb
      end subroutine capecalcnew

!#######################################################################

! lookup table for the analytic evaluation of LCL
      subroutine lcltabl(value,tlcl)
!
! Table of values used to compute the temperature of the lifting condensation
! level.
!
! the expression that we utilize is 
! log(r/theta**(1/kappa)*pstar*rvgas/rdgas/es00) = log(es/T**(1/kappa))
! the RHS is tabulated for the control amount of moisture, hence the 
! division by es00 on the LHS

! Gives the values of the temperature for the following range:
!   starts with -23, is uniformly distributed up to -10.4.  There are a
! total of 127 values, and the increment is .1.
!
      implicit none
      real, intent(in)     :: value
      real, intent(out)    :: tlcl
                                                                                
      integer              :: ival
      real, dimension(127) :: lcltable
      real                 :: v1, v2
                                                                                
      data lcltable/   1.7364512e+02,   1.7427449e+02,   1.7490874e+02, &
      1.7554791e+02,   1.7619208e+02,   1.7684130e+02,   1.7749563e+02, &
      1.7815514e+02,   1.7881989e+02,   1.7948995e+02,   1.8016539e+02, &
      1.8084626e+02,   1.8153265e+02,   1.8222461e+02,   1.8292223e+02, &
      1.8362557e+02,   1.8433471e+02,   1.8504972e+02,   1.8577068e+02, &
      1.8649767e+02,   1.8723077e+02,   1.8797006e+02,   1.8871561e+02, &
      1.8946752e+02,   1.9022587e+02,   1.9099074e+02,   1.9176222e+02, &
      1.9254042e+02,   1.9332540e+02,   1.9411728e+02,   1.9491614e+02, &
      1.9572209e+02,   1.9653521e+02,   1.9735562e+02,   1.9818341e+02, &
      1.9901870e+02,   1.9986158e+02,   2.0071216e+02,   2.0157057e+02, &
      2.0243690e+02,   2.0331128e+02,   2.0419383e+02,   2.0508466e+02, &
      2.0598391e+02,   2.0689168e+02,   2.0780812e+02,   2.0873335e+02, &
      2.0966751e+02,   2.1061074e+02,   2.1156316e+02,   2.1252493e+02, &
      2.1349619e+02,   2.1447709e+02,   2.1546778e+02,   2.1646842e+02, &
      2.1747916e+02,   2.1850016e+02,   2.1953160e+02,   2.2057364e+02, &
      2.2162645e+02,   2.2269022e+02,   2.2376511e+02,   2.2485133e+02, &
      2.2594905e+02,   2.2705847e+02,   2.2817979e+02,   2.2931322e+02, &
      2.3045895e+02,   2.3161721e+02,   2.3278821e+02,   2.3397218e+02, &
      2.3516935e+02,   2.3637994e+02,   2.3760420e+02,   2.3884238e+02, &
      2.4009473e+02,   2.4136150e+02,   2.4264297e+02,   2.4393941e+02, &
      2.4525110e+02,   2.4657831e+02,   2.4792136e+02,   2.4928053e+02, &
      2.5065615e+02,   2.5204853e+02,   2.5345799e+02,   2.5488487e+02, &
      2.5632953e+02,   2.5779231e+02,   2.5927358e+02,   2.6077372e+02, &
      2.6229310e+02,   2.6383214e+02,   2.6539124e+02,   2.6697081e+02, &
      2.6857130e+02,   2.7019315e+02,   2.7183682e+02,   2.7350278e+02, &
      2.7519152e+02,   2.7690354e+02,   2.7863937e+02,   2.8039954e+02, &
      2.8218459e+02,   2.8399511e+02,   2.8583167e+02,   2.8769489e+02, &
      2.8958539e+02,   2.9150383e+02,   2.9345086e+02,   2.9542719e+02, &
      2.9743353e+02,   2.9947061e+02,   3.0153922e+02,   3.0364014e+02, &
      3.0577420e+02,   3.0794224e+02,   3.1014515e+02,   3.1238386e+02, &
      3.1465930e+02,   3.1697246e+02,   3.1932437e+02,   3.2171609e+02, &
      3.2414873e+02,   3.2662343e+02,   3.2914139e+02,   3.3170385e+02 /
                                                                                
      v1 = value
      if (value.lt.-23.0) v1 = -23.0
      if (value.gt.-10.4) v1 = -10.4
      ival = floor(10.*(v1 + 23.0))
      v2 = -230. + ival
      v1 = 10.*v1
      tlcl = (v2 + 1.0 - v1)*lcltable(ival+1) + (v1 - v2)*lcltable(ival+2)
                                                                                
      end subroutine lcltabl




!#######################################################################
subroutine diag_field_init ( axes, Time )

  integer,         intent(in) :: axes(4)
  type(time_type), intent(in) :: Time

  character(len=32) :: tracer_units, tracer_name
  character(len=128) :: diaglname
  integer, dimension(3) :: half = (/1,2,4/)
  integer   :: n, nn

!------------ initializes diagnostic fields in this module -------------

   if ( any((/do_bm,do_bmmass,do_bmomp/)) ) then
      id_qref = register_diag_field ( mod_name, &
        'qref', axes(1:3), Time, &
        'Adjustment reference specific humidity profile', &
        'kg/kg',  missing_value=missing_value               )

      id_tref = register_diag_field ( mod_name, &
        'tref', axes(1:3), Time, &
        'Adjustment reference temperature profile', &
        'K',  missing_value=missing_value                   )

      id_bmflag = register_diag_field (mod_name, &
         'bmflag', axes(1:2), Time, &
         'Betts-Miller flag', &
         'no units', missing_value=missing_value            )

      id_klzbs  = register_diag_field  (mod_name, &
         'klzbs', axes(1:2), Time, &
         'klzb', &
         'no units', missing_value=missing_value            )

      id_cape = register_diag_field ( mod_name, & 
        'cape', axes(1:2), Time, &
        'Convectively available potential energy',      'J/Kg')

      id_cin = register_diag_field ( mod_name, &
        'cin', axes(1:2), Time, &
        'Convective inhibition',                        'J/Kg')
   endif

   if ( do_bm ) then
      id_invtaubmt  = register_diag_field  (mod_name, &
         'invtaubmt', axes(1:2), Time, &
         'Inverse temperature relaxation time', &
         '1/s', missing_value=missing_value            )

      id_invtaubmq = register_diag_field  (mod_name, &
         'invtaubmq', axes(1:2), Time, &
         'Inverse humidity relaxation time', &
         '1/s', missing_value=missing_value            )
   end if  ! if ( do_bm )

   if (do_bmmass) then
      id_massflux = register_diag_field (mod_name, &
         'massflux', axes(1:3), Time, &
         'Massflux implied by temperature adjustment', &
         'm/s', missing_value=missing_value                 )
   end if  ! if ( do_bmmass )

   id_ras_precip = register_diag_field ( mod_name, &
     'ras_precip', axes(1:2), Time, &
    'Precipitation rate from ras ',       'kg/m2/s' )

   id_ras_freq = register_diag_field ( mod_name, &
     'ras_freq', axes(1:2), Time, &
    'frequency of precip from ras ',       'number' , &
         missing_value = missing_value                       )

   id_don_precip = register_diag_field ( mod_name, &
     'don_precip', axes(1:2), Time, &
    'Precipitation rate from donner ',       'kg/m2/s' )

   id_don_freq = register_diag_field ( mod_name, &
     'don_freq', axes(1:2), Time, &
    'frequency of precip from donner ',       'number', &
         missing_value = missing_value                       )

   id_lsc_precip = register_diag_field ( mod_name, &
     'lsc_precip', axes(1:2), Time, &
    'Precipitation rate from lsc ',       'kg/m2/s' )

   id_lsc_freq = register_diag_field ( mod_name, &
     'lsc_freq', axes(1:2), Time, &
    'frequency of precip from lsc ',       'number' , &
         missing_value = missing_value                       )

   id_uw_precip = register_diag_field ( mod_name, &
     'uw_precip', axes(1:2), Time, &
    'Precipitation rate from uw shallow',       'kg/m2/s', &
     interp_method = "conserve_order1" )

   id_uw_freq = register_diag_field ( mod_name, &
     'uw_freq', axes(1:2), Time, &
    'frequency of precip from uw shallow ',       'number' , &
         missing_value = missing_value                       )


   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency from convection ',    'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Spec humidity tendency from convection ',  'kg/kg/s',  &
                        missing_value=missing_value               )

   id_q_conv_col = register_diag_field ( mod_name, &
     'q_conv_col', axes(1:2), Time, &
    'Water vapor path tendency from convection ',   'kg/m2/s' )
   
   id_t_conv_col = register_diag_field ( mod_name, &
     't_conv_col', axes(1:2), Time, &
    'Column static energy tendency from convection ','W/m2' )
   
   id_enth_conv_col = register_diag_field ( mod_name, &
     'enth_conv_col', axes(1:2), Time, &
     'Column enthalpy tendency from convection','W/m2' )
 
   id_wat_conv_col = register_diag_field ( mod_name, &
     'wat_conv_col', axes(1:2), Time, &
     'Column total water tendency from convection','kg(h2o)/m2/s' )

   id_enth_donner_col2 = register_diag_field ( mod_name, &
     'enth_donner_col2', axes(1:2), Time, &
     'column enthalpy tendency from Donner liq precip','W/m2' )
 
   id_enth_donner_col3 = register_diag_field ( mod_name, &
     'enth_donner_col3', axes(1:2), Time, &
      'Column enthalpy tendency from Donner frzn precip','W/m2' )
 
   id_enth_donner_col4 = register_diag_field ( mod_name, &
      'enth_donner_col4', axes(1:2), Time, &
     'Atmospheric column enthalpy tendency from Donner convection', &
                                                            'W/m2' )
 
   id_enth_donner_col5 = register_diag_field ( mod_name, &
      'enth_donner_col5', axes(1:2), Time, &
      'Column enthalpy tendency due to condensate xfer from Donner &
                                                &to lsc','W/m2' )

  id_enth_donner_col6 = register_diag_field ( mod_name, &
     'enth_donner_col6', axes(1:2), Time, &
      'Column enthalpy tendency from donner moisture  &
                     &conservation  adjustment','W/m2' )
 
   id_enth_donner_col7 = register_diag_field ( mod_name, &
      'enth_donner_col7', axes(1:2), Time, &
      'Precip adjustment needed to balance donner moisture  &
                                           &adjustment','kg(h2o)/m2/s' )

   id_enth_donner_col = register_diag_field ( mod_name, &
     'enth_donner_col', axes(1:2), Time, &
     'Column enthalpy imbalance from Donner convection','W/m2' )

   id_wat_donner_col = register_diag_field ( mod_name, &
     'wat_donner_col', axes(1:2), Time, &
  'Column total water tendency from Donner convection','kg(h2o)/m2/s' )

   id_enth_mca_donner_col = register_diag_field ( mod_name, &
     'enth_mca_donner_col', axes(1:2), Time, &
    'Column enthalpy imbalance from Donner MCA convection','W/m2' )

   id_wat_mca_donner_col = register_diag_field ( mod_name, &
     'wat_mca_donner_col', axes(1:2), Time, &
     'Column total water imbalance from Donner MCA convection', &
                                                'kg(h2o)/m2/s' )

   id_enth_uw_col = register_diag_field ( mod_name, &
     'enth_uw_col', axes(1:2), Time, &
     'Column enthalpy tendency from UW convection','W/m2' )
 
   id_wat_uw_col = register_diag_field ( mod_name, &
     'wat_uw_col', axes(1:2), Time, &
      'Column total water tendency from UW convection','kg(h2o)/m2/s' )

   id_scale_uw = register_diag_field ( mod_name, &
     'scale_uw', axes(1:2), Time, &
     'Scaling factor applied to UW convection tendencies','1' )
          
   id_scale_donner = register_diag_field ( mod_name, &
     'scale_donner', axes(1:2), Time, &
     'Scaling factor applied to UW convection tendencies','1' )

   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate from convection ',       'kg(h2o)/m2/s', &
     interp_method = "conserve_order1" )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate from convection ',       'kg(h2o)/m2/s', &
     interp_method = "conserve_order1" )

   id_gust_conv = register_diag_field ( mod_name, &
     'gust_conv', axes(1:2), Time, &
    'Gustiness resulting from convection ',       'm/s' )

  id_conv_rain3d= register_diag_field ( mod_name, &
     'conv_rain3d', axes(half), Time, &
    'Rain fall rate from convection -3D ',       'kg(h2o)/m2/s' )

   id_conv_snow3d= register_diag_field ( mod_name, &
     'conv_snow3d', axes(half), Time, &
    'Snow fall rate from convection -3D',       'kg(h2o)/m2/s' )

   id_lscale_rain3d= register_diag_field ( mod_name, &
     'lscale_rain3d', axes(half), Time, &
    'Rain fall rate from lscale  -3D ',   'kg(h2o)/m2/s' )

   id_lscale_snow3d= register_diag_field ( mod_name, &
     'lscale_snow3d', axes(half), Time, &
    'Snow fall rate from lscale -3D',       'kg(h2o)/m2/s' )
   

    id_max_enthalpy_imbal    = register_diag_field    &
       (mod_name, 'max_enth_imbal', axes(1:2), Time,  &
        'max enthalpy  imbalance from moist_processes  ', 'W/m2',   &
              missing_value=missing_value)
    id_max_water_imbal    = register_diag_field    &
         (mod_name, 'max_water_imbal', axes(1:2), Time,   &
      'max water  imbalance from moist_processes  ', 'kg(h2o)/m2/s',  &
              missing_value=missing_value)

    id_enth_moist_col = register_diag_field ( mod_name, &
     'enth_moist_col', axes(1:2), Time, &
     'Column enthalpy imbalance from moist processes','W/m2' )
  
    id_wat_moist_col = register_diag_field ( mod_name, &
      'wat_moist_col', axes(1:2), Time, &
      'Column total water imbalance from moist processes','kg/m2/s' )

    if (do_donner_conservation_checks) then
      id_enthint    = register_diag_field    &
            (mod_name, 'enthint_don', axes(1:2), Time,  &
          'atmospheric column enthalpy change from donner', 'W/m2',  &
          missing_value=missing_value)
     id_lcondensint    = register_diag_field    &
         (mod_name, 'lcondensint_don', axes(1:2), Time, &
         'enthalpy transferred by condensate from donner to lscale', &
            'W/m2',  missing_value=missing_value)
     id_lprcp    = register_diag_field    &
             (mod_name, 'lprcpint_don', axes(1:2),   &
              Time, 'enthalpy removed by donner precip', 'W/m2',   &
             missing_value=missing_value)
      id_vertmotion    = register_diag_field    &
             (mod_name, 'vertmotion_don', axes(1:2), Time,  &
           'enthalpy change due to cell and meso motion in donner',  &
             'W/m2', missing_value=missing_value)
     id_enthdiffint    = register_diag_field    &
            (mod_name, 'enthdiffint_don', axes(1:2),   &
             Time, 'enthalpy  imbalance due to donner', 'W/m2',   &
             missing_value=missing_value)
     id_vaporint    = register_diag_field    &
           (mod_name, 'vaporint_don', axes(1:2),   &
            Time, 'column water vapor change', 'kg(h2o)/m2/s',   &
            missing_value=missing_value)
     id_max_enthalpy_imbal_don    = register_diag_field    &
            (mod_name, 'max_enth_imbal_don', axes(1:2),   &
              Time, 'max enthalpy  imbalance from donner', 'W/m**2',  &
              missing_value=missing_value)
     id_max_water_imbal_don    = register_diag_field    &
            (mod_name, 'max_water_imbal_don', axes(1:2),   &
              Time, 'max water imbalance from donner', 'kg(h2o)/m2/s', &
         missing_value=missing_value)
     id_condensint    = register_diag_field    &
           (mod_name, 'condensint_don', axes(1:2), Time,  &
         'column condensate exported from donner to lscale', &
                         'kg(h2o)/m2/s',  missing_value=missing_value )
     id_precipint    = register_diag_field    &
            (mod_name, 'precipint_don', axes(1:2),   &
             Time, 'column precip from donner', 'kg(h2o)/m2/s',   &
              missing_value=missing_value)
     id_diffint    = register_diag_field    &
          (mod_name, 'diffint_don', axes(1:2),   &
            Time, 'water imbalance due to donner', 'kg(h2o)/m2/s',   &
              missing_value=missing_value)
  endif



if (do_strat ) then

   id_qldt_conv = register_diag_field ( mod_name, &
     'qldt_conv', axes(1:3), Time, &
     'Liquid water tendency from convection',      'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qndt_conv = register_diag_field ( mod_name, &
     'qndt_conv', axes(1:3), Time, &
     'Liquid drop tendency from convection',      '#/kg/s',  &
                         missing_value=missing_value               )

   id_qidt_conv = register_diag_field ( mod_name, &
     'qidt_conv', axes(1:3), Time, &
     'Ice water tendency from convection',         'kg/kg/s',  &
                        missing_value=missing_value               )

   id_qadt_conv = register_diag_field ( mod_name, &
     'qadt_conv', axes(1:3), Time, &
     'Cloud fraction tendency from convection',    '1/sec',    &
                        missing_value=missing_value               )

   id_ql_conv_col = register_diag_field ( mod_name, &
     'ql_conv_col', axes(1:2), Time, &
    'Liquid water path tendency from convection',  'kg/m2/s' )
   
   id_qn_conv_col = register_diag_field ( mod_name, &
     'qn_conv_col', axes(1:2), Time, &
     'Liquid drp tendency from convection',  'kg/m2/s' )
 
   id_qi_conv_col = register_diag_field ( mod_name, &
     'qi_conv_col', axes(1:2), Time, &
    'Ice water path tendency from convection',     'kg/m2/s' )
   
   id_qa_conv_col = register_diag_field ( mod_name, &
     'qa_conv_col', axes(1:2), Time, &
    'Cloud mass tendency from convection',         'kg/m2/s' )
      
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
    'Precipitation rate from large-scale cond',     'kg/m2/s', &
     interp_method = "conserve_order1" )

   id_snow_ls = register_diag_field ( mod_name, &
     'snow_ls', axes(1:2), Time, &
    'Frozen precip rate from large-scale cond',     'kg/m2/s', &
     interp_method = "conserve_order1" )

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

   id_qndt_ls = register_diag_field ( mod_name, &
     'qndt_ls', axes(1:3), Time, &
     'Drop number tendency from strat cloud',        '#/kg/s',  &
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
   
   id_qn_ls_col = register_diag_field ( mod_name, &
     'qn_ls_col', axes(1:2), Time, &
     'Column drop number tendency from strat cloud',  '#/m2/s' )

   id_qi_ls_col = register_diag_field ( mod_name, &
     'qi_ls_col', axes(1:2), Time, &
    'Ice water path tendency from strat cloud',      'kg/m2/s' )
   
   id_qa_ls_col = register_diag_field ( mod_name, &
     'qa_ls_col', axes(1:2), Time, &
    'Cloud mass tendency from strat cloud',          'kg/m2/s' )
      
   id_enth_ls_col = register_diag_field ( mod_name, &
     'enth_ls_col', axes(1:2), Time, &
     'Column enthalpy tendency from strat cloud','W/m2' )
 
   id_wat_ls_col = register_diag_field ( mod_name, &
     'wat_ls_col', axes(1:2), Time, &
     'Column total water tendency from strat cloud','kg/m2/s' )

endif

   id_precip = register_diag_field ( mod_name, &
     'precip', axes(1:2), Time, &
     'Total precipitation rate',                     'kg/m2/s', &
      interp_method = "conserve_order1" )

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

   id_qs = register_diag_field ( mod_name, &
     'qs', axes(1:3), Time, &
         'saturation specific humidity',                 'kg/kg',    & 
                        missing_value=missing_value               )
   
if (do_donner_deep) then

   id_tdt_deep_donner= register_diag_field ( mod_name, &
           'tdt_deep_donner', axes(1:3), Time, &
           ' heating rate - deep portion', 'deg K/s', &
                        missing_value=missing_value               )

   id_qdt_deep_donner = register_diag_field ( mod_name, &
           'qdt_deep_donner', axes(1:3), Time, &
           ' moistening rate - deep portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qadt_deep_donner = register_diag_field ( mod_name, &
     'qadt_deep_donner', axes(1:3), Time, &
     ' cloud amount tendency - deep portion', '1/s', &
                        missing_value=missing_value               )

   id_qldt_deep_donner = register_diag_field ( mod_name, &
     'qldt_deep_donner', axes(1:3), Time, &
     ' cloud liquid tendency - deep portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qidt_deep_donner = register_diag_field ( mod_name, &
     'qidt_deep_donner', axes(1:3), Time, &
     ' ice water tendency - deep portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_tdt_mca_donner = register_diag_field ( mod_name, &
     'tdt_mca_donner', axes(1:3), Time, &
     ' heating rate - mca  portion', 'deg K/s', &
                        missing_value=missing_value               )

   id_qdt_mca_donner = register_diag_field ( mod_name, &
           'qdt_mca_donner', axes(1:3), Time, &
           ' moistening rate - mca  portion', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_prec_deep_donner = register_diag_field ( mod_name, &
           'prc_deep_donner', axes(1:2), Time, &
           ' total precip rate - deep portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_prec1_deep_donner = register_diag_field ( mod_name, &
           'prc1_deep_donner', axes(1:2), Time, &
           ' change in precip for conservation in donner', 'kg/m2/s ', &
              missing_value=missing_value, mask_variant = .true., &
             interp_method = "conserve_order1"  )

   id_prec_mca_donner = register_diag_field ( mod_name, &
           'prc_mca_donner', axes(1:2), Time, &
           ' total precip rate - mca  portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_snow_deep_donner = register_diag_field ( mod_name, &
           'snow_deep_donner', axes(1:2), Time, &
           ' frozen precip rate - deep portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_snow_mca_donner = register_diag_field ( mod_name, &
           'snow_mca_donner', axes(1:2), Time, &
           ' frozen precip rate -  mca portion', 'kg/m2/s', &
                        missing_value=missing_value, &
             interp_method = "conserve_order1"               )

   id_mc_donner = register_diag_field ( mod_name, &
           'mc_donner', axes(1:3), Time, &
           'Net Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

   id_m_cdet_donner = register_diag_field ( mod_name, &
           'm_cdet_donner', axes(1:3), Time, &
           'Detrained Cell Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

   id_m_cellup = register_diag_field ( mod_name, &
           'm_cellup', axes(half), Time, &
           'Upward Cell Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

endif


if (do_uw_conv) then

   id_tdt_uw = register_diag_field ( mod_name, &
           'tdt_uw', axes(1:3), Time, &
           'UW convection heating rate', 'deg K/s', &
                        missing_value=missing_value               )

   id_qdt_uw = register_diag_field ( mod_name, &
           'qdt_uw', axes(1:3), Time, &
           'UW convection moistening rate', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qadt_uw = register_diag_field ( mod_name, &
           'qadt_uw', axes(1:3), Time, &
           'UW convection cloud amount tendency', '1/s', &
                        missing_value=missing_value               )

   id_qldt_uw = register_diag_field ( mod_name, &
           'qldt_uw', axes(1:3), Time, &
           'UW convection cloud liquid tendency', 'kg/kg/s', &
                        missing_value=missing_value               )

   id_qidt_uw = register_diag_field ( mod_name, &
           'qidt_uw', axes(1:3), Time, &
           'UW convection ice water tendency', 'kg/kg/s', &
                        missing_value=missing_value               )

   if (do_liq_num) &
    id_qndt_uw = register_diag_field ( mod_name, &
           'qndt_uw', axes(1:3), Time, &
           'UW convection cloud drop tendency', '#/kg/s', &
                        missing_value=missing_value               )

endif

   id_qvout = register_diag_field ( mod_name, &
           'qvout', axes(1:3), Time, 'qv after strat_cloud', 'kg/kg', &
                        missing_value=missing_value               )

   id_qaout = register_diag_field ( mod_name, &
           'qaout', axes(1:3), Time, 'qa after strat_cloud', 'none', &
                        missing_value=missing_value               )

   id_qlout = register_diag_field ( mod_name, &
           'qlout', axes(1:3), Time, 'ql after strat_cloud', 'kg/kg', &
                        missing_value=missing_value               )

   id_qiout = register_diag_field ( mod_name, &
           'qiout', axes(1:3), Time, 'qi after strat_cloud', 'kg/kg', &
                        missing_value=missing_value               )

!---------------------------------------------------------------------
!    register diagnostics for lightning NOx
!---------------------------------------------------------------------

   if (get_tracer_index(MODEL_ATMOS,'no') > 0) &
     id_prod_no = register_diag_field ( 'tracers', &
             'hook_no', axes(1:3), Time, &
             'hook_no',   'molec/cm3/s')

!-----------------------------------------------------------------------
!---------------------------------------------------------------------
!    register the diagnostics associated with convective tracer 
!    transport.
!---------------------------------------------------------------------
      allocate (id_tracerdt_conv    (num_tracers))
      allocate (id_tracerdt_conv_col(num_tracers))
      allocate (id_wet_deposition(num_tracers))
      allocate (id_conv_tracer           (num_tracers))
      allocate (id_conv_tracer_col(num_tracers))
 
      do n = 1,num_tracers
        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                               units = tracer_units)
        if (tracers_in_donner(n) .or. &
            tracers_in_ras(n)      .or.  &
            tracers_in_mca(n)      .or.  &
            tracers_in_uw(n)) then
          diaglname = trim(tracer_name)//  &
                        ' total tendency from moist convection'
          id_tracerdt_conv(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

          diaglname = trim(tracer_name)//  &
                       ' total path tendency from moist convection'
          id_tracerdt_conv_col(n) =  &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv_col', &
                         axes(1:2), Time, trim(diaglname), &
                         TRIM(tracer_units)//'*(kg/m2)/s',   &
                         missing_value=missing_value)
         endif
 
         diaglname = trim(tracer_name)
         id_conv_tracer(n) =    &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name),  &
                        axes(1:3), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,  &
                        missing_value=missing_value)
         diaglname =  ' column integrated' // trim(tracer_name)
         id_conv_tracer_col(n) =  &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name)//'_col', &
                        axes(1:2), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,   &
                        missing_value=missing_value)
         id_wet_deposition(n) = register_diag_field( mod_name, &
           trim(tracer_name)//'_wetdep', axes(1:3), Time, &
           trim(tracer_name)//' tendency from wet deposition',TRIM(tracer_units)//'/sec', &
           missing_value=missing_value )
      end do

!------------------------------------------------------------------
!    register the variables associated with the mca component of 
!    donner_deep transport.
!------------------------------------------------------------------
     if (do_donner_deep) then
       allocate (id_tracerdt_mcadon  (num_donner_tracers))
       allocate (id_tracerdt_mcadon_col(num_donner_tracers))
 
       nn = 1
       do n = 1,num_tracers
         call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                                units = tracer_units)
         if (tracers_in_donner(n) ) then
           diaglname = trim(tracer_name)//  &
                       ' tendency from donner-mca'
           id_tracerdt_mcadon(nn) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_donmca',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                        missing_value=missing_value)

           diaglname = trim(tracer_name)//  &
                       ' total path tendency from donner-mca'
           id_tracerdt_mcadon_col(nn) =  &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name)//'_donmca_col', &
                        axes(1:2), Time, trim(diaglname), &
                          TRIM(tracer_units)//'*(kg/m2)/s',   &
                        missing_value=missing_value)
           nn = nn + 1
         endif
       end do
 

     endif

end subroutine diag_field_init

!#######################################################################



                 end module moist_processes_mod

