
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

use     donner_deep_mod, only: donner_deep_init,               &
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

use     strat_cloud_mod, only: strat_cloud_init, strat_cloud, strat_cloud_end, &
                               strat_cloud_sum

use       rh_clouds_mod, only: rh_clouds_init, rh_clouds_end, &
                               rh_clouds_sum

use      diag_cloud_mod, only: diag_cloud_init, diag_cloud_end, &
                               diag_cloud_sum

use   diag_integral_mod, only: diag_integral_field_init, &
                               sum_diag_integral_field

use       constants_mod, only: CP_AIR, GRAV, RDGAS, RVGAS, TFREEZE, &
                               SECONDS_PER_DAY

use     cu_mo_trans_mod, only: cu_mo_trans_init, cu_mo_trans, cu_mo_trans_end

use   field_manager_mod, only: MODEL_ATMOS
use  tracer_manager_mod, only: get_tracer_index,&
                               get_number_tracers, &
                               get_tracer_names, &
                               query_method, &
                               NO_TRACER
use atmos_tracer_utilities_mod, only : wet_deposition
use             fms_mod, only: mpp_clock_id, mpp_clock_begin, &
                               mpp_clock_end, CLOCK_MODULE, &
                               MPP_CLOCK_SYNC

use       moz_hook_mod, only : moz_hook

use rad_utilities_mod,       only: aerosol_type

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

   integer :: nsphum, nql, nqi, nqa, nqn  ! tracer indices for stratiform clouds
!--------------------- version number ----------------------------------
   character(len=128) :: &
   version = '$Id: moist_processes.F90,v 13.0.2.1 2006/07/11 17:23:51 wfc Exp $'
   character(len=128) :: tagname = '$Name: memphis_2006_08 $'
   logical            :: module_is_initialized = .false.
!-----------------------------------------------------------------------
!-------------------- namelist data (private) --------------------------


   logical :: do_mca=.true., do_lsc=.true., do_ras=.false.,  &
              do_strat=.false., do_dryadj=.false., &
              do_rh_clouds=.false., do_diag_clouds=.false., &
              do_donner_deep=.false., do_cmt=.false., &
              use_tau=.false., do_gust_cv = .false., &
              do_liq_num = .false.


   character(len=64)  :: cmt_mass_flux_source = 'ras'
   real :: pdepth = 150.e2
   real :: gustmax = 3.             ! maximum gustiness wind (m/s)
   real :: gustconst = 10./SECONDS_PER_DAY   ! constant in kg/m2/sec, default =
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
!               module; currently either 'ras', 'donner' or 
!               'donner_and_ras'
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
                               do_dryadj, pdepth,                 &
                               cmt_mass_flux_source, &
                               use_tau, do_rh_clouds, do_diag_clouds, &
                               do_donner_deep, do_cmt, do_gust_cv, &
                               gustmax, gustconst, do_liq_num

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
           id_prod_no, id_m_cdet_donner, id_m_cellup, &
           id_conv_rain3d, id_conv_snow3d, id_lscale_rain3d, id_lscale_snow3d
 
integer, dimension(:), allocatable :: id_tracerdt_conv,  &
                                      id_tracerdt_conv_col, &
                                      id_conv_tracer,  &
                                      id_conv_tracer_col, &
                                      id_tracerdt_mcadon, &
                                      id_tracerdt_mcadon_col, &
                                      id_wet_deposition
character(len=5) :: mod_name = 'moist'

real :: missing_value = -999.
integer :: convection_clock, largescale_clock

logical :: do_tracers_in_donner =.false.
logical :: do_tracers_in_mca = .false.
logical :: do_tracers_in_ras = .false.
logical, dimension(:), allocatable :: tracers_in_donner,   &
                                      tracers_in_mca, tracers_in_ras
integer :: num_donner_tracers=0
integer :: num_mca_tracers=0
integer :: num_ras_tracers=0
integer :: num_tracers=0

logical :: cmt_uses_donner = .false.
logical :: cmt_uses_ras = .false.

!-----------------------------------------------------------------------




                             contains




!#######################################################################

subroutine moist_processes (is, ie, js, je, Time, dt, land,            &
                            phalf, pfull, zhalf, zfull, omega, diff_t, &
                            radturbten,                                &
                            t, q, r, u, v, tm, qm, rm, um, vm,         &
                            tdt, qdt, rdt, udt, vdt, diff_cu_mo,       &
                            convect, lprec, fprec, gust_cv, area,      &
                            lat, Aerosol, mask, kbot, &
                            cell_cld_frac, cell_liq_amt, cell_liq_size, &
                            cell_ice_amt, cell_ice_size, &
                            meso_cld_frac, meso_liq_amt, meso_liq_size, &
                            meso_ice_amt, meso_ice_size, nsum_out)

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


   type(aerosol_type),intent(in),       optional :: Aerosol
   real, intent(in) , dimension(:,:,:), optional :: mask
integer, intent(in) , dimension(:,:),   optional :: kbot

   real, intent(inout) , dimension(:,:,:), optional :: &      
                          cell_cld_frac, cell_liq_amt, cell_liq_size, &
                          cell_ice_amt, cell_ice_size, &
                          meso_cld_frac, meso_liq_amt, meso_liq_size, &
                          meso_ice_amt, meso_ice_size
integer, intent(inout) , dimension(:,:), optional ::  nsum_out

!-----------------------------------------------------------------------
real, dimension(size(t,1),size(t,2),size(t,3)) :: tin,qin,ttnd,qtnd, &
                                                  cf, rlin, riin, &
                                                  qnew, qinew, qlnew, delta_temp, &
                                                  delta_vapor,  &
                                                  delta_q, rin, rtnd,   &
                                                  donner_humidity_area, &
                                                  donner_humidity_ratio
real, dimension(size(t,1),size(t,2),size(phalf,3)) :: conv_rain3d, conv_snow3d
real, dimension(size(t,1),size(t,2),size(phalf,3)) :: lscale_rain3d, lscale_snow3d
real, dimension(size(t,1),size(t,2),size(t,3)) :: qtnd_wet, cloud_wet
real, dimension(size(t,1),size(t,3)) :: dp, temp
real, dimension(size(t,1),size(t,2),size(t,3)) :: ttnd_conv,qtnd_conv
real, dimension(size(t,1),size(t,2))           :: tsnow,snow
logical,dimension(size(t,1),size(t,2))         :: coldT
real, dimension(size(t,1),size(t,2),size(t,3)) :: utnd,vtnd,uin,vin

real, dimension(size(t,1),size(t,2),size(t,3)) :: qltnd,qitnd,qatnd, qntnd,&
                                        delta_ql, delta_qi, delta_qa, &
                                                  qlin, qiin, qain
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: tracer, rdt_init , qtrcumo
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: mc,mask3, &
                                                    m_cellup, mc_cmt
real, dimension(size(t,1),size(t,2),size(t,3)) :: det0, det_cmt       
real, dimension(size(t,1),size(t,2),size(t,3)) :: mc_full, mc_donner, &
                                                  m_cdet_donner
real, dimension(size(t,1),size(t,2),size(t,3)) :: RH, pmass, wetdeptnd
real, dimension(size(t,1),size(t,2))           :: rain, precip, precip_dd
real, dimension(size(t,1),size(t,2))           :: wvp,lwp,iwp
real, dimension(size(t,1),size(t,2))           :: sfc_sh_flux,  &
                                                  sfc_vapor_flux
!     sfc_sh_flux      sensible heat flux across the surface
!                      [ watts / m**2 ]
!     sfc_vapor_flux   water vapor flux across the surface
!                      [ kg(h2o) / (m**2 sec) ]

real, dimension(size(t,1),size(t,2))           :: tempdiag
real, dimension(size(t,1),size(t,2),size(t,3)) :: tempdiag1
integer n

integer :: i, j, k, ix, jx, kx, nt, ip, tr
real    :: dtinv
logical ::           do_adjust, used

real, dimension(size(t,1),size(t,2),size(t,3)+1) :: press
real, dimension(size(t,1),size(t,2),size(t,3)) :: tin_in

real :: tinsave

! tracer code:
integer :: nn
real, dimension (size(t,1), size(t,2), &
                 size(t,3), num_donner_tracers) :: qtrtnd , &
                                                   donner_tracers
real, dimension (size(t,1), size(t,2), &
                 size(t,3), num_ras_tracers) :: qtrras , &
                                                    ras_tracers
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

real, dimension(size(t,1),size(t,2),size(t,3)) ::           cnvcntq

character(len=32) :: tracer_units, tracer_name

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

!--------------------------------------------------------------------
!    define the inverse of the time step.
!--------------------------------------------------------------------
      dtinv = 1.0/dt

!--------------------------------------------------------------------
!    initialize the arrays which will be output from this subroutine.
!--------------------------------------------------------------------
      lprec(:,:)  = 0.0  
      fprec(:,:)  = 0.0
      convect(:,:)  = .false.
      gust_cv(:,:) = 0.0
      diff_cu_mo(:,:,:) = 0.0

!---------------------------------------------------------------------
!    define input fields to be used, either the tau time level fields,
!    or the tau - 1 time level values updated with the time tendencies
!    thus far calculated on the current step. control is through nml
!    variable use_tau.
!---------------------------------------------------------------------
      if (use_tau) then
        tin(:,:,:) = t(:,:,:)
        qin(:,:,:) = q(:,:,:)
        uin(:,:,:) = u(:,:,:)
        vin(:,:,:) = v(:,:,:)
        do tr=1,size(r,4)
          tracer(:,:,:,tr) = r(:,:,:,tr)
        end do  
      else
        tin(:,:,:) = tm(:,:,:) + tdt(:,:,:)*dt
        qin(:,:,:) = qm(:,:,:) + qdt(:,:,:)*dt
        uin(:,:,:) = um(:,:,:) + udt(:,:,:)*dt
        vin(:,:,:) = vm(:,:,:) + vdt(:,:,:)*dt
        do tr=1,size(rdt,4)
          tracer(:,:,:,tr) = rm(:,:,:,tr) + rdt(:,:,:,tr)*dt
        end do  
        do tr=size(rdt,4) +1, size(r,4)
          tracer(:,:,:,tr) = r(:,:,:,tr)
        end do  
      endif

!--------------------------------------------------------------------
!    if using eta vertical coordinate, define the appropriate values 
!    for any points located below the ground. values of 0.0 are given
!    to u, v and q, and a temperature value of EPST (=200. K) is given 
!    to sub-surface  points.
!--------------------------------------------------------------------
      if (present(mask) .and. present(kbot))  then
        tin(:,:,:) = mask(:,:,:)*tin(:,:,:) + (1.0 - mask(:,:,:))*EPST 
        qin(:,:,:) = mask(:,:,:)*qin(:,:,:)
        uin(:,:,:) = mask(:,:,:)*uin(:,:,:)
        vin(:,:,:) = mask(:,:,:)*vin(:,:,:)
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
        if (id_conv_tracer(n) > 0) then
          used = send_data (id_conv_tracer(n), tracer(:,:,:,n),   &
                            Time, is, js, 1, rmask=mask)
        endif
        if (id_conv_tracer_col(n) > 0) then
          tempdiag(:,:)=0.
          do k=1,kx
            tempdiag(:,:) = tempdiag(:,:) + tracer(:,:,k,n)*pmass(:,:,k)
          end do
          used = send_data (id_conv_tracer_col(n), tempdiag, Time,   &
                            is, js)
        endif
      end do

!---------------------------------------------------------------------
!    initialize local arrays which will hold sums.
!---------------------------------------------------------------------
      rdt_init(:,:,:,:) = rdt(:,:,:,:)
      ttnd_conv(:,:,:) = 0.
      qtnd_conv(:,:,:) = 0.
!-->cjg: addition to initialize arrays that were causing model to crash
!        when do_cma = do_ras = do_donner_deep = .false.
      qtnd(:,:,:) = 0.
      qltnd(:,:,:) = 0.
      qitnd(:,:,:) = 0.
      conv_rain3d(:,:,:) = 0.
      conv_snow3d(:,:,:) = 0.
!<--cjg

!----------------------------------------------------------------------
!    compute the mean temperature in the lower atmosphere (the lowest
!    pdepth Pa), to be used to determine whether rain or snow reaches
!    the surface. define a logical variable coldT indicating whether
!    snow or rain falls in the column.
!    ????    SHOULD TIN BE USED RATHER THAN t ??
!----------------------------------------------------------------------
      call tempavg (pdepth, phalf, t, tsnow, mask)
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
        call dry_adj (tin, pfull, phalf, delta_temp, mask)

!-------------------------------------------------------------------
!    add the temperature change due to dry adjustment to the current
!    temperature. convert the temperature change to a heating rate and
!    add that to the temperature temndency array accumulating the ten-
!    dencies due to all physics processes.
!-------------------------------------------------------------------
        tin  = tin + delta_temp
        ttnd = delta_temp*dtinv
        tdt  = tdt + ttnd

!---------------------------------------------------------------------
!    output the temperature tendency from dry adjustment, if desired.
!---------------------------------------------------------------------
        if (id_tdt_dadj > 0) then
          used = send_data (id_tdt_dadj, ttnd, Time, is, js, 1, &
                            rmask=mask )
        endif

!---------------------------------------------------------------------
!    add the temperature time tendency from dry adjustment to the array
!    accumulating the total temperature time tendency from convection.
!---------------------------------------------------------------------
        ttnd_conv = ttnd_conv + ttnd
      endif


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                  MOIST CONVECTION PARAMETERIZATIONS
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                A. DONNER DEEP CONVECTION PARAMETERIZATION
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    if donner_deep convection is activated, execute the following code.
!---------------------------------------------------------------------
      if (do_donner_deep) then
        
!---------------------------------------------------------------------
!    be sure all optional arguments associated with the donner param-
!    eterization are present.
!---------------------------------------------------------------------
        if    &
           (present (cell_cld_frac) .and.   &
            present (cell_liq_amt) .and. present ( cell_liq_size) .and. &
            present (cell_ice_amt) .and. present ( cell_ice_size) .and. &
            present (meso_cld_frac) .and.   &
            present (meso_liq_amt) .and. present ( meso_liq_size) .and. &
            present (meso_ice_amt) .and. present ( meso_ice_size) .and. &
            present (nsum_out) ) then
        else
          call error_mesg ('moist_processes_mod', 'moist_processes: &
               &not all 11 optional arguments needed for donner_deep &
              &output are present', FATAL)
        endif

!--------------------------------------------------------------------
!    if strat_cloud_mod is activated, define the cloud liquid and 
!    cloud ice specific humidities and cloud area associated with 
!    strat_cloud_mod, so that they may be input to donner_deep_mod. 
!    if strat_cloud_mod is not activated, define these arrays to be 
!    zero. 
!--------------------------------------------------------------------
        if (do_strat) then
          qlin = tracer(:,:,:,nql)
          qiin = tracer(:,:,:,nqi)
          qain = tracer(:,:,:,nqa)
        endif

!--------------------------------------------------------------------
!    convert vapor specific humidity to vapor mixing ratio so it may
!    be input to donner_deep_mod.
!--------------------------------------------------------------------
        rin = qin/(1.0 - qin)

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

!---------------------------------------------------------------------
!    call donner_deep to compute the effects of deep convection on the 
!    temperature, vapor mixing ratio, tracers, cloud liquid, cloud ice
!    cloud area and precipitation fields.
!---------------------------------------------------------------------
        if (do_strat) then
          call donner_deep (is, ie, js, je, dt, tin, rin, pfull,       &
                            phalf, omega, land, sfc_sh_flux,   &
                            sfc_vapor_flux, tr_flux, donner_tracers, &
                            Time,     &
                            cell_cld_frac, cell_liq_amt, cell_liq_size, &
                            cell_ice_amt, cell_ice_size, &
                            meso_cld_frac, meso_liq_amt, meso_liq_size, &
                            meso_ice_amt, meso_ice_size, nsum_out,  &
                            precip_dd, delta_temp, delta_vapor, &
                            m_cdet_donner, m_cellup, mc_donner,  &
                            donner_humidity_area,  &
                            donner_humidity_ratio, qtrtnd,    &
                            qlin, qiin, qain, delta_ql, &   ! optional
                            delta_qi, delta_qa)             ! optional
        else
          call donner_deep (is, ie, js, je, dt, tin, rin, pfull,      &
                            phalf, omega, land, sfc_sh_flux,    &
                            sfc_vapor_flux, tr_flux, donner_tracers, &
                            Time,                           &
                            cell_cld_frac, cell_liq_amt, cell_liq_size, &
                            cell_ice_amt, cell_ice_size, &
                            meso_cld_frac, meso_liq_amt, meso_liq_size, &
                            meso_ice_amt, meso_ice_size, nsum_out,  &
                            precip_dd, delta_temp, delta_vapor, &
                            m_cdet_donner, m_cellup, mc_donner,  &
                            donner_humidity_area,  &
                            donner_humidity_ratio, qtrtnd)
        endif

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

!--------------------------------------------------------------------
!    define the predicted deep convective precipitation rate for  
!    either rain or snow based on the mean temperature in the bottom 
!    pdepth pascals of the atmosphere, as indicated by the logical 
!    variable coldT.
!--------------------------------------------------------------------
        where (coldT)
          snow = precip_dd*dtinv
          rain = 0.
        elsewhere
          rain = precip_dd*dtinv
          snow = 0.
        end where
                    
!--------------------------------------------------------------------
!    obtain updated vapor specific humidity (qnew) resulting from deep 
!    convection. define the vapor specific humidity change due to deep 
!    convection (qtnd).
!--------------------------------------------------------------------
        qnew = (rin + delta_vapor)/ (1.0 + (rin + delta_vapor))
        delta_q = qnew - qin

!--------------------------------------------------------------------
!    if strat_cloud_mod is activated, update the local tracer arrays 
!    (tracer) of cloud liquid, cloud ice and cloud area. add the tenden-
!    cies to the arrays accumulating tracer tendencies from all physics 
!    processes (rdt).
!--------------------------------------------------------------------
        if (do_strat) then
          qlnew = (rlin + delta_ql)/ (1.0 + (rin + delta_vapor))
          qinew = (riin + delta_qi)/ (1.0 + (rin + delta_vapor))
          tracer(:,:,:,nql) = qlin(:,:,:) + delta_ql(:,:,:)
          tracer(:,:,:,nqi) = qiin(:,:,:) + delta_qi(:,:,:)
          tracer(:,:,:,nqa) = qain(:,:,:) + delta_qa(:,:,:)
          rdt(:,:,:,nql) = rdt(:,:,:,nql) + delta_ql(:,:,:)*dtinv
          rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + delta_qi(:,:,:)*dtinv
          rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + delta_qa(:,:,:)*dtinv
        endif

!---------------------------------------------------------------------
!    update the values of temperature and vapor specific humidity to
!    include the effects of deep convection.
!---------------------------------------------------------------------
        tin = tin + delta_temp
        qin = qin + delta_q
      
!---------------------------------------------------------------------
!    convert the changes in temperature, vapor specific humidity and 
!    precipitation resulting from deep convection to time tendencies 
!    of these quantities.
!---------------------------------------------------------------------
        ttnd = delta_temp*dtinv 
        qtnd = delta_q*dtinv

!--------------------------------------------------------------------
!    output the time tendencies of temperature, vapor specific humid-
!    ity, precipitation and mass flux due to deep convection.
!--------------------------------------------------------------------
        if (id_tdt_deep_donner > 0) then
          used = send_data (id_tdt_deep_donner, ttnd, Time,   &
                            is, js, 1, rmask=mask )
        endif

        if (id_qdt_deep_donner > 0) then
          used = send_data (id_qdt_deep_donner, qtnd, Time,  &
                            is, js, 1, rmask=mask )
        endif

        if (id_prec_deep_donner > 0) then
          used = send_data (id_prec_deep_donner, rain + snow, Time,  &
                            is, js )
        endif

        if (id_snow_deep_donner > 0) then
          used = send_data (id_snow_deep_donner, snow, Time, is, js)
        endif

        if ( id_mc_donner > 0 ) then
          used = send_data (id_mc_donner, mc_donner, Time,  &
                            is, js, 1, rmask=mask )
        endif

        if ( id_m_cdet_donner > 0 ) then
          used = send_data (id_m_cdet_donner, m_cdet_donner, Time,  &
                            is, js, 1, rmask=mask )
        endif

        if ( id_m_cellup > 0 ) then
          used = send_data (id_m_cellup, m_cellup, Time,  &
                            is, js, 1, rmask=mask )
        endif

!--------------------------------------------------------------------
!    add the contributions to the temperature and vapor specific 
!    humidity tendencies from donner_deep mod to the arrays accumulating
!    the total tendencies due to all physics processes.
!--------------------------------------------------------------------
        tdt = tdt + ttnd 
        qdt = qdt + qtnd

!--------------------------------------------------------------------
!    add the liquid (rain) and frozen (snow) precipitation generated by
!    deep convection on this step to the arrays accumulating precip-
!    itation from all sources (lprec, fprec).
!--------------------------------------------------------------------
        lprec  = lprec + rain
        fprec  = fprec + snow

!--------------------------------------------------------------------
!    save the tendencies of temperature and specific humidity resulting
!    from the deep convection component of the donner parameterization. 
!--------------------------------------------------------------------
        ttnd_conv = ttnd_conv + ttnd
        qtnd_conv = qtnd_conv + qtnd

!--------------------------------------------------------------------
!    call subroutine moist_conv to handle any shallow convection 
!    present in the grid. in this call do_strat is always set to .false.
!    so that no convective detrainment (and corresponding change in
!    large-scale cloud amount and area) from moist convective adjustment
!    is allowed, consistent with this call being constrained to handle
!    shallow convection.
!--------------------------------------------------------------------
        call moist_conv (tin, qin, pfull, phalf, coldT, ttnd, qtnd,   &
                         rain, snow, dtinv, Time, is, js,    &
                         donner_tracers, qtrtnd, Lbot=kbot, mask=mask)           

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
        if (id_tdt_mca_donner > 0) then
          used = send_data (id_tdt_mca_donner, ttnd, Time, is, js, 1, &
                            rmask=mask)
        endif

        if (id_qdt_mca_donner > 0) then
          used = send_data (id_qdt_mca_donner, qtnd, Time, is, js, 1, &
                            rmask=mask)
        endif

        if (id_prec_mca_donner > 0) then
          used = send_data (id_prec_mca_donner, rain+snow, Time,  &
                            is, js) 
        endif

        if (id_snow_mca_donner > 0) then
          used = send_data (id_snow_mca_donner, snow, Time, is, js)
        endif

!--------------------------------------------------------------------
!------- diagnostics for tracers from convection -------
!  allow any tracer to be activated here (allows control cases)
!--------------------------------------------------------------------
      do n=1,num_tracers
        if ( id_conv_tracer(n) > 0 ) then
          used = send_data ( id_conv_tracer(n), tracer(:,:,:,n), Time, is, js, 1, &
                            rmask=mask )
         endif
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
          if ( id_tracerdt_mcadon(n) > 0 ) then
            used = send_data (id_tracerdt_mcadon(n), qtrtnd(:,:,:,n), &
                              Time, is, js, 1, rmask=mask )
          endif
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
!    add the contributions to the temperature and vapor specific 
!    humidity tendencies from the moist convective adjustment pass of
!    donner_deep_mod to the arrays accumulating the total tendencies 
!    due to all physics processes.
!--------------------------------------------------------------------
        tdt = tdt + ttnd 
        qdt = qdt + qtnd

!--------------------------------------------------------------------
!    add the liquid (rain) and frozen (snow) precipitation generated by
!    the moist convective adjustment pass of the donner parameterization
!    on this step to the arrays accumulating precipitation from all 
!    sources (lprec, fprec).
!--------------------------------------------------------------------
        lprec  = lprec + rain
        fprec  = fprec + snow

!--------------------------------------------------------------------
!    define the heating, moistening and precipitation rates as the sum 
!    of the contributions from the deep convection pass and the moist 
!    convective adjustment pass of the donner parameterization. if 
!    ras_mod is also activated, store these values in temporary arrays
!    until the contributions from ras_mod is calculated.
!--------------------------------------------------------------------
        ttnd_conv = ttnd_conv + ttnd
        qtnd_conv = qtnd_conv + qtnd

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
        donner_humidity_ratio = 1.
      endif  ! (do_donner_deep)


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

!----------------------------------------------------------------------
!    add the temperature and specific humidity tendencies from moist
!    convective adjustment (ttnd, qtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt).
!----------------------------------------------------------------------
        tdt = tdt + ttnd 
        qdt = qdt + qtnd

!----------------------------------------------------------------------
!    if strat_cloud_mod is activated, add the cloud liquid, ice and area
!    tendencies from moist convective adjustment to the 
!    arrays accumulating these tendencies from all physics processes 
!    (rdt).
!----------------------------------------------------------------------
        if (do_strat) then
          rdt(:,:,:,nql) = rdt(:,:,:,nql) + qltnd(:,:,:)
          rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + qitnd(:,:,:)
          rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + qatnd(:,:,:)
        endif

!----------------------------------------------------------------------
!    increment the liquid, solid and total precipitation fields with 
!    the contribution from moist convective adjustment.
!----------------------------------------------------------------------
        lprec  = lprec  + rain
        fprec  = fprec  + snow
        ttnd_conv = ttnd_conv + ttnd
        qtnd_conv = qtnd_conv + qtnd
      endif ! (do_mca)


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
        if (do_strat .and. .not.(do_liq_num)) then
          call ras (is,   js,     Time,     tin,   qin,   &
                    uin,  vin,    pfull,    phalf, zhalf, coldT, &
                    dt,   ttnd,   qtnd,     utnd,  vtnd,  &
                    conv_rain3d, conv_snow3d, &
                    rain, snow,   ras_tracers, qtrras,    &
                    mask,  kbot, mc, det0,     &
                    tracer(:,:,:,nql), tracer(:,:,:,nqi), &
                    tracer(:,:,:,nqa), qltnd(:,:,:),&
                    qitnd(:,:,:), qatnd(:,:,:))       

       elseif(do_strat .and. do_liq_num) then
          call ras (is,   js,     Time,     tin,   qin,   &
                    uin,  vin,    pfull,    phalf, zhalf, coldT, &
                    dt,   ttnd,   qtnd,     utnd,  vtnd,  &
                    conv_rain3d, conv_snow3d, &
                    rain, snow,   ras_tracers, qtrras, &
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
                    rain, snow,   ras_tracers, qtrras,           &
                    mask,  kbot,  mc, det0)
        endif

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
!----------------------------------------------------------------------
!    add the temperature, specific humidity and momentum tendencies 
!    from ras (ttnd, qtnd, utnd, vtnd) to the arrays accumulating 
!    these tendencies from all physics processes (tdt, qdt, udt, vdt).
!----------------------------------------------------------------------
        tdt = tdt + ttnd 
        qdt = qdt + qtnd
        udt = udt + utnd
        vdt = vdt + vtnd

!----------------------------------------------------------------------
!    if strat_cloud_mod is activated, add the cloud liquid, ice and area
!    tendencies from ras to the arrays accumulating these tendencies 
!    from all physics processes (rdt).
!----------------------------------------------------------------------
        if (do_strat) then
          rdt(:,:,:,nql) = rdt(:,:,:,nql) + qltnd(:,:,:)
          rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + qitnd(:,:,:)
          rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + qatnd(:,:,:)
          if (do_liq_num) rdt(:,:,:,nqn) = rdt(:,:,:,nqn) + qntnd(:,:,:)
        endif

!----------------------------------------------------------------------
!    increment the liquid, solid and total precipitation fields with 
!    the contribution from ras.
!----------------------------------------------------------------------
        lprec  = lprec  + rain
        fprec  = fprec  + snow

!---------------------------------------------------------------------
!    if donner_deep_mod is also active, define the total time tendency 
!    due to all moist convective processes (donner (including its mca 
!    part), and ras) of temperature, specific humidity, rain, snow and,
!    if strat_cloud_mod is activated, the cloud liquid, cloud ice and 
!    cloud area. 
!---------------------------------------------------------------------
        ttnd_conv = ttnd_conv + ttnd
        qtnd_conv = qtnd_conv + qtnd
!---------------------------------------------------------------------
!    if ras_mod is not activated, set the ras mass flux field to be 0.0.
!---------------------------------------------------------------------
      else
        mc(:,:,:) = 0.0
        det0(:,:,:) = 0.0
      endif  ! (do_ras)

!---------------------------------------------------------------------
!    call subroutine cu_mo_trans if diffusive cumulus momentum 
!    transport is desired. 
!---------------------------------------------------------------------
        if (do_cmt) then
          mc_cmt = 0.
          det_cmt = 0.
          if (cmt_uses_ras) then
            mc_cmt = mc_cmt + mc
            det_cmt = det_cmt + det0
          endif
          if (cmt_uses_donner) then
            mc_cmt = mc_cmt + m_cellup 
            det_cmt = det_cmt + m_cdet_donner 
          endif
          call cu_mo_trans (is, js, Time, mc_cmt, tin, phalf, pfull,  &
                            zhalf, zfull, dt, uin, vin, tracer, pmass, &
                            det_cmt, utnd, vtnd, ttnd, qtrcumo, &
                            diff_cu_mo  )

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
        tdt = tdt + ttnd 
        udt = udt + utnd
        vdt = vdt + vtnd

!---------------------------------------------------------------------
!    if donner_deep_mod is also active, define the total time tendency 
!    due to all moist convective processes (donner (including its mca 
!    part), and ras) of temperature, specific humidity, rain, snow and,
!    if strat_cloud_mod is activated, the cloud liquid, cloud ice and 
!    cloud area. 
!---------------------------------------------------------------------
          ttnd_conv = ttnd_conv + ttnd
        endif  ! (do_cmt)

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
      qtnd_wet = qtnd
      if (do_strat) then
         qtnd_wet = qtnd_wet + qltnd + qitnd
         cloud_wet = conv_rain3d(:,:,2:kx+1)-conv_rain3d(:,:,1:kx) &
                     + conv_snow3d(:,:,2:kx+1)-conv_snow3d(:,:,1:kx)
         cloud_wet = MAX(cloud_wet,1.e-20) * dt / pmass(:,:,:) ! kg/kg
      else
         cloud_wet = qtnd_wet * dt
      end if
      do n=1,size(rdt,4)
        if ( n /= nsphum ) then
          if ( .not. do_strat .or. (n /= nql .and. n /= nqi .and. n /= nqa) ) then
            wetdeptnd = 0.0

                    call wet_deposition( n, t, pfull, phalf, zfull, zhalf, rain, snow, &
                                 qtnd_wet, cloud_wet, &
                                 conv_rain3d, conv_snow3d, &
                                 tracer(:,:,:,n), wetdeptnd, &
                                 Time, 'convect', is, js, dt )
            rdt (:,:,:,n) = rdt(:,:,:,n) - wetdeptnd(:,:,:)
            wet_data(:,:,:,n) = wetdeptnd(:,:,:)
          endif
        endif  
      end do
      do k=1,kx
        mc_full(:,:,k) = 0.5*(mc(:,:,k) + mc(:,:,k+1)) +   &
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
   call moz_hook(cldtop, cldbot, land, zfull, zhalf, t, prod_no, area, lat, &
                 Time, is, js)
   rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) = &
      rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) + &
      prod_no/conc_air
   if ( id_prod_no > 0 ) then
      used = send_data(id_prod_no,prod_no, Time, is_in=is, js_in=js)
   endif
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


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!                   CONVECTION DIAGNOSTICS      
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    temperature change due to dry and moist convection:
!---------------------------------------------------------------------
      if (id_tdt_conv > 0) then
        used = send_data (id_tdt_conv, ttnd_conv, Time, is, js, 1, &
                          rmask=mask)
      endif

!---------------------------------------------------------------------
!    vapor specific humidity change due to convection:
!---------------------------------------------------------------------
      if (id_qdt_conv > 0) then
        used = send_data (id_qdt_conv, qtnd_conv, Time, is, js, 1, &
                          rmask=mask)
      endif

!---------------------------------------------------------------------
!    total precipitation due to convection:
!---------------------------------------------------------------------
      if (id_prec_conv > 0) then
        used = send_data (id_prec_conv, precip, Time, is, js)
      endif

!---------------------------------------------------------------------
!    frozen precipitation (snow) due to convection:
!---------------------------------------------------------------------
      if (id_snow_conv > 0) then
        used = send_data (id_snow_conv, fprec, Time, is, js)
      endif

!---------------------------------------------------------------------
!------- diagnostics for 3D precip_conv -------
!---------------------------------------------------------------------
      if ( id_conv_rain3d > 0 ) then
        used = send_data ( id_conv_rain3d, conv_rain3d, Time, is, js, 1 )
      endif

!---------------------------------------------------------------------
!------- diagnostics for 3D snow_conv -------
!---------------------------------------------------------------------
      if ( id_conv_snow3d > 0 ) then
        used = send_data ( id_conv_snow3d, conv_snow3d, Time, is, js, 1 )
      endif

!---------------------------------------------------------------------
!    surface wind gustiness due to convection:
!---------------------------------------------------------------------
      if (id_gust_conv > 0) then
        used = send_data (id_gust_conv, gust_cv, Time, is, js)
      endif

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
          if (id_qldt_conv > 0) then
            used = send_data (id_qldt_conv, tempdiag1, Time, &
                              is, js, 1, rmask=mask)
          endif

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
!    if cloud ice diagnostics requested:
!---------------------------------------------------------------------
        if (id_qidt_conv > 0 .or. id_qi_conv_col > 0) then
          tempdiag1 = rdt(:,:,:,nqi) - rdt_init(:,:,:,nqi)

!---------------------------------------------------------------------
!    cloud ice tendency due to convection:
!---------------------------------------------------------------------
          if (id_qidt_conv > 0) then
            used = send_data (id_qidt_conv, tempdiag1, Time, &
                              is, js, 1, rmask=mask)
          endif

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

!---------------------------------------------------------------------
!    cloud area tendency due to convection:
!---------------------------------------------------------------------
          tempdiag1 = rdt(:,:,:,nqa) - rdt_init(:,:,:,nqa)
          if (id_qadt_conv > 0 ) then          
            used = send_data (id_qadt_conv, tempdiag1, Time, &
                              is, js, 1, rmask=mask)
          endif

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
!    tracer tendencies due to convection:
!---------------------------------------------------------------------
      do n=1,size(rdt,4)
        if (tracers_in_donner(n) .or. &
            tracers_in_ras(n)    .or.  &
            tracers_in_mca(n))    then
          if (id_tracerdt_conv(n) > 0 .or.       &
              id_tracerdt_conv_col(n) > 0) then
            tempdiag1 = rdt(:,:,:,n) - rdt_init(:,:,:,n)
            if (id_tracerdt_conv(n) > 0) then
              used = send_data (id_tracerdt_conv(n), tempdiag1, &
                                Time, is, js, 1, rmask=mask )
            endif

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
        call lscale_cond (tin, qin, pfull, phalf, coldT, rain, snow,  &
                          ttnd, qtnd, mask=mask)

!-----------------------------------------------------------------------
!    add the temperature and specific humidity increments to the updated
!    temperature and specific humidity fields (tin, qin). convert these
!    increments and the precipitation increments to rates and add to 
!    the arrays accumulating the total rates for all physical processes
!    (tdt, qdt, lprec, fprec).
!-----------------------------------------------------------------------
        tin   = tin   + ttnd 
        qin   = qin   + qtnd
        tdt   = tdt   + ttnd*dtinv 
        qdt   = qdt   + qtnd*dtinv
        lprec = lprec + rain*dtinv
        fprec = fprec + snow*dtinv

!--------------------------------------------------------------------
!    if rh_clouds is active, call rh_calc to determine the grid box
!    relative humidity. call rh_clouds_sum to pass this field to 
!    rh_clouds_mod so it may be used to determine the grid boxes which
!    will contain clouds for the radiation package.
!---------------------------------------------------------------------
        if (do_rh_clouds) then
          call rh_calc (pfull, tin, qin, rh, mask)
          call rh_clouds_sum (is, js, rh)
        endif

!--------------------------------------------------------------------
!    if the gordon diagnostic cloud parameterization is active, set a 
!    flag to indicate those grid points where drying has resulted from 
!    convective activity (cnvcntq). call rh_calc to determine the grid 
!    box relative humidity. call diag_cloud_sum to define the cloud 
!    field that will be seen by the radiation package.
!---------------------------------------------------------------------
        if (do_diag_clouds) then
          where (qtnd_conv(:,:,:) < 0.0)
            cnvcntq (:,:,:) = 1.0
          elsewhere
            cnvcntq (:,:,:) = 0.0
          end where
          call rh_calc (pfull, tin, qin, rh, mask)
          call diag_cloud_sum (is, js, tin, qin, rh, omega, qtnd,  &
                               cnvcntq, precip, kbot)
        endif



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

!--------------------------------------------------------------------
!    call strat_cloud to integrate the prognostic cloud equations. 
!-----------------------------------------------------------------------
        if (do_liq_num) then 
          call strat_cloud (Time, is, ie, js, je, dt, pfull, phalf,    & 
                         radturbten, tin, qin, tracer(:,:,:,nql), &
                         tracer(:,:,:,nqi), tracer(:,:,:,nqa),  &
                         omega, mc_full, diff_t, land, ttnd, qtnd,  &
                         qltnd, qitnd, qatnd, lscale_rain3d, lscale_snow3d, &
                         rain, snow, &
                         donner_humidity_ratio, donner_humidity_area, &
                         mask=mask, qn=tracer(:,:,:,nqn), Aerosol=Aerosol, SN=qntnd)
        else 
          call strat_cloud (Time, is, ie, js, je, dt, pfull, phalf,    & 
                         radturbten, tin, qin, tracer(:,:,:,nql), &
                         tracer(:,:,:,nqi), tracer(:,:,:,nqa),  &
                         omega, mc_full, diff_t, land, ttnd, qtnd,  &
                         qltnd, qitnd, qatnd, lscale_rain3d, lscale_snow3d, &
                         rain, snow, &
                         donner_humidity_ratio, donner_humidity_area, &
                         mask=mask)
        endif
    
!----------------------------------------------------------------------
!    upon return from strat_cloud, update the cloud liquid, ice and area.
!    update the temperature and specific humidity fields.
!----------------------------------------------------------------------
        tracer(:,:,:,nql) = tracer(:,:,:,nql) + qltnd               
        tracer(:,:,:,nqi) = tracer(:,:,:,nqi) + qitnd                 
        tracer(:,:,:,nqa) = tracer(:,:,:,nqa) + qatnd                 
        if (do_liq_num) tracer(:,:,:,nqn) = tracer(:,:,:,nqn) + qntnd                 
        tin = tin + ttnd 
        qin = qin + qtnd

!----------------------------------------------------------------------
!    call strat_cloud_sum to make the cloud variables available for 
!    access by the radiation package. NOTE: this is no longer necessary,
!    and can be judiciously removed (provided other affiliated code 
!    and options are nullified).
!----------------------------------------------------------------------
        call strat_cloud_sum (is, js, tracer(:,:,:,nql),  &
                              tracer(:,:,:,nqi), tracer(:,:,:,nqa))

!----------------------------------------------------------------------
!    convert increments to tendencies.
!----------------------------------------------------------------------
        ttnd = ttnd*dtinv 
        qtnd = qtnd*dtinv
        rain = rain*dtinv 
        snow = snow*dtinv
        qltnd = qltnd*dtinv
        qitnd = qitnd*dtinv
        qatnd = qatnd*dtinv
        if (do_liq_num) qntnd = qntnd*dtinv
   
!----------------------------------------------------------------------
!    update the total tendency terms (temperature, vapor specific 
!    humidity, cloud liquid, cloud ice, cloud area, liquid precip,
!    frozen precip) with the contributions from the strat_cloud scheme.
!----------------------------------------------------------------------
        tdt = tdt + ttnd 
        qdt = qdt + qtnd
        rdt(:,:,:,nql) = rdt(:,:,:,nql) + qltnd                 
        rdt(:,:,:,nqi) = rdt(:,:,:,nqi) + qitnd                 
        rdt(:,:,:,nqa) = rdt(:,:,:,nqa) + qatnd                  
        if (do_liq_num)         rdt(:,:,:,nqn) = rdt(:,:,:,nqn) + qntnd                  
        lprec = lprec + rain
        fprec = fprec + snow
      endif  ! (do_lsc)

!---------------------------------------------------------------------
!    calculate the wet deposition associated with the large scale 
!    condensation. 
!---------------------------------------------------------------------
      qtnd_wet = qtnd
      if (do_strat) then
         qtnd_wet = qtnd_wet + qltnd + qitnd
         cloud_wet = tracer(:,:,:,nql) + tracer(:,:,:,nqi)
      else
         cloud_wet = qtnd_wet * dt
      end if
      do n=1,size(rdt,4)
        if ( n /= nsphum ) then
          if ( .not. do_strat .or. (n /= nql .and. n /= nqi .and. n /= nqa .and. n /= nqn) ) then
            wetdeptnd = 0.0
            call wet_deposition( n, t, pfull, phalf, zfull, zhalf, rain, snow, &
                                 qtnd_wet, cloud_wet, &
                                 lscale_rain3d, lscale_snow3d, &
                                 tracer(:,:,:,n), wetdeptnd, &
                                 Time, 'lscale', is, js, dt )
            rdt (:,:,:,n) = rdt(:,:,:,n) - wetdeptnd
            wet_data(:,:,:,n) = wet_data(:,:,:,n) + wetdeptnd

            if (id_wet_deposition(n) /= 0 ) then
              used = send_data( id_wet_deposition(n), wet_data(:,:,:,n), &
                                Time,is_in=is,js_in=js )
            endif
          end if
        end if
      end do


!---------------------------------------------------------------------
!    output diagnostics associated with the large-scale condensation
!    scheme.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    temperature change due to large-scale condensation:
!---------------------------------------------------------------------
      if (id_tdt_ls > 0) then
        used = send_data (id_tdt_ls, ttnd, Time, is, js, 1, &
                          rmask=mask)
      endif

!---------------------------------------------------------------------
!    specific humidity change due to large-scale condensation:
!---------------------------------------------------------------------
      if (id_qdt_ls > 0) then
        used = send_data (id_qdt_ls, qtnd, Time, is, js, 1, &
                          rmask=mask)
      endif

!---------------------------------------------------------------------
!    total precipitation rate due to large-scale condensation:
!---------------------------------------------------------------------
      if (id_prec_ls > 0) then
        used = send_data (id_prec_ls, rain+snow, Time, is, js)
      endif

!---------------------------------------------------------------------
!    snowfall rate due to large-scale condensation:
!---------------------------------------------------------------------
      if (id_snow_ls > 0) then
        used = send_data (id_snow_ls, snow, Time, is, js)
      endif

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
        if (id_mc_full > 0) then
          used = send_data (id_mc_full, mc_full, Time, is, js, 1, &
                             rmask=mask)
        endif

!---------------------------------------------------------------------
!    cloud liquid, ice and area tendencies due to strat_cloud 
!    parameterization:
!---------------------------------------------------------------------
        if (id_qldt_ls > 0) then
          used = send_data (id_qldt_ls, qltnd, Time,  &
                            is, js, 1, rmask=mask)
        endif
        if (id_qidt_ls > 0) then
          used = send_data (id_qidt_ls, qitnd, Time, &
                            is, js, 1, rmask=mask)
        endif
        if (id_qadt_ls > 0) then
          used = send_data (id_qadt_ls, qatnd, Time, &
                            is, js, 1, rmask=mask)
        endif

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
        if (id_qi_ls_col > 0) then
          tempdiag(:,:) = 0.
          do k=1,kx
            tempdiag(:,:) = tempdiag(:,:) + qitnd(:,:,k)*pmass(:,:,k)
          end do
          used = send_data (id_qi_ls_col, tempdiag, Time, is, js)
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
        if ( id_lscale_rain3d > 0) then
          used = send_data(id_lscale_rain3d, lscale_rain3d, Time, is, js, 1)
        endif

!---------------------------------------------------------------------
!---- diagnostics for large scale snow -------------
!---------------------------------------------------------------------
        if ( id_lscale_snow3d > 0 ) then
          used = send_data(id_lscale_snow3d, lscale_snow3d, Time, is, js, 1)
        endif

      endif ! (do_strat)

!---------------------------------------------------------------------
!    end the timing of the large-scale condensation code section.
!---------------------------------------------------------------------
      call mpp_clock_end (largescale_clock)
 


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
      if (id_precip > 0) then
        used = send_data (id_precip, precip, Time, is, js)
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
!    output the global integral of precipitation in units of mm/day.
!---------------------------------------------------------------------
      call sum_diag_integral_field ('prec', precip*SECONDS_PER_DAY,  &
                                                                is, js)

!-----------------------------------------------------------------------



end subroutine moist_processes




!#######################################################################

subroutine moist_processes_init ( id, jd, kd, lonb, latb, pref, &
                                  axes, Time, doing_donner)

!-----------------------------------------------------------------------
integer,         intent(in) :: id, jd, kd, axes(4)
real, dimension(:), intent(in) :: lonb, latb, pref
type(time_type), intent(in) :: Time
logical,         intent(out) :: doing_donner
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
      if (do_dryadj) call     dry_adj_init ()
      if (do_cmt)    call cu_mo_trans_init (axes,Time)

      if (do_cmt) then
        if (trim(cmt_mass_flux_source) == 'ras') then
          cmt_uses_ras = .true.
          cmt_uses_donner = .false.
          if (.not. do_ras) then
            call error_mesg ('moist_processes_mod', &
              'if cmt_uses_ras then ras_mod must be active', FATAL)
           endif
        else if (trim(cmt_mass_flux_source) == 'donner') then
          cmt_uses_ras = .false.
          cmt_uses_donner = .true.
          if (.not. do_donner_deep)  then
            call error_mesg ('moist_processes_mod', &
                'if cmt_uses_donner then donner_deep_mod must&
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
        else
          call error_mesg ('moist_processes_mod', &
             'invalid specification of cmt_mass_flux_source', FATAL)
        endif

      endif

  
!----- initialize quantities for global integral package -----

   call diag_integral_field_init ('prec', 'f6.3')
   

!----- initialize clocks -----

   convection_clock =     &
       mpp_clock_id( '   Physics_up: Moist Proc: Conv', &
           grain=CLOCK_MODULE, flags = MPP_CLOCK_SYNC )
   largescale_clock =     &
       mpp_clock_id( '   Physics_up: Moist Proc: LS', &
           grain=CLOCK_MODULE, flags = MPP_CLOCK_SYNC )

     doing_donner = do_donner_deep

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
      tracers_in_donner = .false.
      tracers_in_mca = .false.
      tracers_in_ras = .false.

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

!--------------------------------------------------------------------
!    initialize the convection scheme modules.
!--------------------------------------------------------------------
      if (do_donner_deep) then
        call donner_deep_init (lonb, latb, pref, axes, Time,  &
                               tracers_in_donner)
      endif ! (do_donner_deep)
 
      if (do_ras)  then
        call ras_init (do_strat, do_liq_num, axes,Time, tracers_in_ras)
      endif

      if (do_mca .or. do_donner_deep)  then
        call  moist_conv_init (axes,Time, tracers_in_mca)
      endif
  
 
!----- initialize quantities for diagnostics output -----
 
      call diag_field_init ( axes, Time )
 

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

  character(len=32) :: tracer_units, tracer_name
  character(len=128) :: diaglname
  integer, dimension(3) :: half = (/1,2,4/)
  integer   :: n, nn

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

  id_conv_rain3d= register_diag_field ( mod_name, &
     'conv_rain3d', axes(half), Time, &
    'Rain fall rate from donner_deep plus RAS -3D ',       'kg/m2/s' )

   id_conv_snow3d= register_diag_field ( mod_name, &
     'conv_snow3d', axes(half), Time, &
    'Snow fall rate from donner_deep plus RAS -3D',       'kg/m2/s' )

   id_lscale_rain3d= register_diag_field ( mod_name, &
     'lscale_rain3d', axes(half), Time, &
    'Rain fall rate from lscale  -3D ',   'kg/m2/s' )

   id_lscale_snow3d= register_diag_field ( mod_name, &
     'lscale_snow3d', axes(half), Time, &
    'Snow fall rate from lscale -3D',       'kg/m2/s' )
   


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

   id_m_cdet_donner = register_diag_field ( mod_name, &
           'm_cdet_donner', axes(1:3), Time, &
           'Detrained Cell Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

   id_m_cellup = register_diag_field ( mod_name, &
           'm_cellup', axes(half), Time, &
           'Upward Cell Mass Flux from donner',   'kg/m2/s', &
                        missing_value=missing_value               )

endif

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
            tracers_in_mca(n))    then
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
                         TRIM(tracer_units)//'/s',   &
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
                        TRIM(tracer_units)//'/s',   &
                        missing_value=missing_value)
           nn = nn + 1
         endif
       end do
 

     endif

end subroutine diag_field_init

!#######################################################################


                 end module moist_processes_mod
