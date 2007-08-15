
MODULE UW_CONV_MOD

  use           mpp_mod, only : mpp_pe, mpp_root_pe, stdlog
  use      Constants_Mod, ONLY: tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas
  use   Diag_Manager_Mod, ONLY: register_diag_field, send_data
  use   Time_Manager_Mod, ONLY: time_type, get_time 
  use           fms_mod, only : write_version_number, open_namelist_file, &
                                FILE_EXIST, ERROR_MESG,  &
                                lowercase, &
                                CLOSE_FILE, FATAL
  use  field_manager_mod, only: MODEL_ATMOS
  use  tracer_manager_mod, only: get_tracer_names, query_method, &
                                 get_tracer_index, NO_TRACER
  use  sat_vapor_pres_mod,only : sat_vapor_pres_init
  use atmos_tracer_utilities_mod, only : get_wetdep_param

  use  rad_utilities_mod, only : aerosol_type
  
  use  conv_utilities_mod,only :   uw_params_init
  use  conv_utilities_k_mod,only : sd_init_k, sd_copy_k, sd_end_k,  &
                                   ac_init_k, ac_clear_k, ac_end_k, &
                                   pack_sd_k, adi_cloud_k, extend_sd_k,&
                                   sd_save_k, ac_save_k, &
                                   sd_retrieve_k, ac_retrieve_k, &
                                   exn_init_k, exn_end_k, findt_init_k,&
                                   findt_end_k, &
                                   check_tracer_realizability, &
                                   qt_parcel_k, &
                                   adicloud, sounding, uw_params

  use  conv_plumes_k_mod,only    : cp_init_k, cp_end_k, cp_clear_k, &
                                   ct_init_k, ct_end_k, ct_clear_k, &
                                   cumulus_tend_k, cumulus_downdraft_k,&
                                   cumulus_plume_k, &
                                   cplume, ctend, cpnlist

  use  conv_closures_mod,only    : cclosure_bretherton,   &
                                   cclosure_relaxcbmf, &
                                   cclosure_relaxwfn,  &
!                                  cclosure_unified, &
                                   cclosure_implicit, cclosure

!---------------------------------------------------------------------
  implicit none
  private
!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: uw_conv.F90,v 15.0 2007/08/14 03:56:10 fms Exp $'
  character(len=128) :: tagname = '$Name: omsk $'

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: uw_conv, uw_conv_init, uw_conv_end, calculate_uw_closure

  real, parameter :: mv = -999.
  logical         :: module_is_initialized = .false.

  character(len=7) :: mod_name = 'uw_conv'

  !namelist parameters for UW convection scheme
  integer :: iclosure = 0      ! 0: Bretherton UWShCu orginal / -CIN/TKE based
                               ! 1: Emanuel-Rayment: quasiequilibrium PBL
  real    :: rkm_dp   = 1.0    ! fractional lateral mixing rate for deep
  real    :: rkm_sh   = 16.0   ! fractional lateral mixing rate for shallow
  real    :: cldhgt_max   = 4.e3
  real    :: cbmf_sh_frac = 1.0
  logical :: do_deep = .false.
  logical :: do_relaxcape = .false.
  logical :: do_relaxwfn  = .false.
  logical :: do_coldT = .true.
  logical :: do_lands = .false.
  logical :: do_uwcmt = .false.   
  logical :: do_fast  = .false.
  logical :: do_ice   = .true.
  logical :: do_ppen  = .true.
  logical :: do_micro = .false.   
  logical :: do_edplume = .true.
  logical :: do_forcedlifting = .false.
  real    :: atopevap = 0.
  logical :: apply_tendency = .true.
  logical :: prevent_unreasonable = .true.
  real    :: aerol = 1.e-12
  real    :: gama     = 1.0    ! 
  real    :: tkemin   = 1.e-6
  real    :: wmin_ratio = 0.
  logical :: use_online_aerosol = .false.
  logical :: do_auto_aero = .false.

  NAMELIST / uw_conv_nml / iclosure, rkm_dp, rkm_sh, cldhgt_max, cbmf_sh_frac, &
       do_deep, do_relaxcape, do_relaxwfn, do_coldT, do_lands, do_uwcmt,       &
       do_fast, do_ice, do_ppen, do_micro, do_edplume, do_forcedlifting,       &
       atopevap, apply_tendency, prevent_unreasonable, aerol, gama, tkemin,    &
       wmin_ratio, use_online_aerosol, do_auto_aero
       

  !namelist parameters for UW convective plume
  real    :: rle      = 0.10   ! for critical stopping distance for entrainment
  real    :: rpen     = 5.0    ! for entrainment efficiency
  real    :: rmaxfrac = 0.05   ! maximum allowable updraft fraction
  real    :: wmin     = 0.0    ! maximum allowable updraft fraction
  real    :: rbuoy    = 1.0    ! for nonhydrostatic pressure effects on updraft
  real    :: rdrag    = 1.0 
  real    :: frac_drs = 1.0    ! 
  real    :: bigc     = 0.7    ! for momentum transfer
  real    :: auto_th0 = 0.5e-3 ! threshold for precipitation
  real    :: auto_rate= 1.e-3
  real    :: tcrit    = -45.0  ! critical temperature 
  real    :: rad_crit = 14.0   ! critical droplet radius

  NAMELIST / uw_plume_nml / rle, rpen, rmaxfrac, wmin, rbuoy, rdrag, frac_drs, bigc, &
       auto_th0, auto_rate, tcrit, rad_crit
 
  !namelist parameters for UW convective closure
  integer :: igauss   = 1      ! options for cloudbase massflux closure
                               ! 1: cin/gaussian closure, using TKE to compute CIN.
                               ! 2: cin/gaussian closure, using W* to compute CIN.
                               ! 0: cin/tke mapse-style closure; 
  real    :: rkfre    = 0.05   ! vertical velocity variance as fraction of tke
  real    :: tau_dp   = 7200.  ! 
  real    :: tau_sh   = 7200.  ! 
  real    :: wcrit_min= 0.

  NAMELIST / uw_closure_nml / igauss, rkfre, tau_dp, tau_sh, wcrit_min

!------------------------------------------------------------------------

  integer :: nqv, nql, nqi, nqa ,nqn
  logical :: do_qn = .false.    ! use droplet number tracer field ?

  integer :: id_tdt_uwc, id_qdt_uwc, id_prec_uwc, id_snow_uwc,               &
       id_cin_uwc, id_cbmf_uwc, id_tke_uwc, id_plcl_uwc, id_zinv_uwc,  &
       id_cush_uwc, id_pct_uwc, id_pcb_uwc, id_plfc_uwc, id_enth_uwc,  &
       id_qldt_uwc, id_qidt_uwc, id_qadt_uwc, id_qndt_uwc, id_cmf_uwc, id_wu_uwc,   &
       id_fer_uwc,  id_fdr_uwc, id_fdrs_uwc, id_cqa_uwc, id_cql_uwc,   &
       id_cqi_uwc,  id_cqn_uwc, id_hlflx_uwc, id_qtflx_uwc,           &
       id_cape_uwc, id_dcin_uwc, id_dcape_uwc, id_dwfn_uwc,            &
       id_ocode_uwc, id_plnb_uwc, id_wrel_uwc, id_ufrc_uwc, id_qtmp_uwc

  integer, allocatable :: id_tracerdt_uwc(:), id_tracerdt_uwc_col(:), &
                          id_tracerdtwet_uwc(:), id_tracerdtwet_uwc_col(:)

  real, dimension(:,:), allocatable :: wrel_mod, ufrc_mod, scaleh_mod, &
                                       dcin_mod, dcape_mod, dwfn_mod, &
                                      wfn_mod, cbmfo_mod, &
                                      cbmfo_unmodified

  logical   ::  do_unified_convective_closure

  type(sounding),   save  :: sd, sd1
  type(adicloud),   save  :: ac, ac1
  type(cclosure),   save  :: cc, cc1
  type(cplume)  ,   save  :: cp, cp1
  type(ctend)   ,   save  :: ct, ct1
  type(cpnlist),    save  :: cpn
  type(uw_params),  save  :: Uw_p

  type(sounding), dimension(:,:), allocatable :: sdij
  type(adicloud), dimension(:,:), allocatable :: acij
  logical, dimension(:,:), allocatable :: exit_flag
contains

!#####################################################################
!#####################################################################

  SUBROUTINE UW_CONV_INIT( do_strat, axes, Time, kd, tracers_in_uw, &
                           doing_unified_closure)
    logical,         intent(in) :: do_strat
    integer,         intent(in) :: axes(4), kd
    type(time_type), intent(in) :: Time
    logical,         intent(in) :: tracers_in_uw(:)
    logical,         intent(in) :: doing_unified_closure
    
!---------------------------------------------------------------------
!  intent(in) variables:
!
!      tracers_in_uw 
!                   logical array indicating which of the activated 
!                   tracers are to be transported by UW convection
!
!-------------------------------------------------------------------

    integer   :: unit, io
    
    integer   :: ntracers, n, nn
    logical   :: flag
    character(len=200) :: text_in_scheme, control
 
    ntracers = count(tracers_in_uw)

    call sd_init_k(kd,ntracers,sd);
    call sd_init_k(kd,ntracers,sd1);
    call ac_init_k(kd,ac);
    call ac_init_k(kd,ac1);
    call cp_init_k(kd,ntracers,cp)
    call cp_init_k(kd,ntracers,cp1)
    call ct_init_k(kd,ntracers,ct)
    call ct_init_k(kd,ntracers,ct1)
    call uw_params_init   (Uw_p)


!   Initialize lookup tables needed for findt and exn
!   sat_vapor_pres needs to be initialized if not already done
    call sat_vapor_pres_init     
    call exn_init_k (Uw_p)
    call findt_init_k (Uw_p)

    if( FILE_EXIST( 'input.nml' ) ) then
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_closure_nml, iostat = io, end = 10 )
       end do
10     call close_file ( unit )

       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_conv_nml, iostat = io, end = 20 )
       end do
20     call close_file ( unit )
       
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_plume_nml, iostat = io, end = 30 )
       end do
30     call close_file ( unit )
    end if
    call write_version_number (version, tagname)
    WRITE( stdlog(), nml = uw_closure_nml )
    WRITE( stdlog(), nml = uw_conv_nml )
    WRITE( stdlog(), nml = uw_plume_nml )

    do_unified_convective_closure = doing_unified_closure

    !pack namelist parameters into plume and closure structure
    cpn % rle       = rle
    cpn % rpen      = rpen
    cpn % rmaxfrac  = rmaxfrac
    cpn % wmin      = wmin
    cpn % rbuoy     = rbuoy
    cpn % rdrag     = rdrag  
    cpn % frac_drs  = frac_drs
    cpn % bigc      = bigc    
    cpn % auto_th0  = auto_th0
    cpn % auto_rate = auto_rate
    cpn % tcrit     = tcrit  
    cpn % cldhgt_max= cldhgt_max
    cpn % do_ice    = do_ice
    cpn % do_ppen   = do_ppen
    cpn % do_edplume= do_edplume
    cpn % do_micro  = do_micro
    cpn % do_forcedlifting= do_forcedlifting
    cpn % atopevap  = atopevap
    cpn % wtwmin_ratio = wmin_ratio*wmin_ratio
    cpn % do_auto_aero = do_auto_aero
    cpn % rad_crit = rad_crit

    nqv = get_tracer_index ( MODEL_ATMOS, 'sphum' )
    nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
    nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
    nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
    nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
    if (nqn /= NO_TRACER) do_qn = .true.
    if (ntracers > 0) then
      allocate ( cpn%tracername   (ntracers) )
      allocate ( cpn%tracer_units (ntracers) )
      allocate ( cpn%wetdep       (ntracers) )
      nn = 1
      do n=1,size(tracers_in_uw(:))
         if (tracers_in_uw(n)) then
             call get_tracer_names (MODEL_ATMOS, n,  &
                                    name = cpn%tracername(nn), &
                                    units = cpn%tracer_units(nn))
             flag = query_method( 'wet_deposition', MODEL_ATMOS, n, &
                                  text_in_scheme, control )
             call get_wetdep_param( text_in_scheme, control, &
                                    cpn%wetdep(nn)%scheme, &
                                    cpn%wetdep(nn)%Henry_constant, &
                                    cpn%wetdep(nn)%Henry_variable, &
                                    cpn%wetdep(nn)%frac_in_cloud, &
                                    cpn%wetdep(nn)%alpha_r, &
                                    cpn%wetdep(nn)%alpha_s )
             cpn%wetdep(nn)%scheme = lowercase( cpn%wetdep(nn)%scheme )
             nn = nn + 1
          endif
       end do
    endif
    cc  % igauss    = igauss
    cc  % rkfre     = rkfre
    cc  % rmaxfrac  = rmaxfrac
    cc  % wcrit_min = wcrit_min
    cc  % rbuoy     = rbuoy
    cc  % tau_dp    = tau_dp
    cc  % tau_sh    = tau_sh
 
    id_tdt_uwc = register_diag_field ( mod_name, 'tdt_uwc', axes(1:3), Time, &
         'Temperature tendency from uw_conv', 'K/day', missing_value=mv )
    id_qdt_uwc = register_diag_field ( mod_name, 'qdt_uwc', axes(1:3), Time, &
         'Spec. humidity tendency from uw_conv', 'kg/kg/day', missing_value=mv)
    id_cmf_uwc = register_diag_field ( mod_name, 'cmf_uwc', axes(1:3), Time, &
         'Cloud vert. mass flux from uw_conv', 'kg/m2/s', missing_value=mv)
    id_wu_uwc = register_diag_field ( mod_name, 'wu_uwc', axes(1:3), Time,   &
         'Updraft vert. velocity from uw_conv', 'm/s', missing_value=mv)
    id_fer_uwc = register_diag_field ( mod_name, 'fer_uwc', axes(1:3), Time, &
         'Fractional entrainment rate from uw_conv', '1/Pa', missing_value=mv)
    id_fdr_uwc = register_diag_field ( mod_name, 'fdr_uwc', axes(1:3), Time, &
         'Fractional detrainment rate from uw_conv', '1/Pa', missing_value=mv)
    id_fdrs_uwc = register_diag_field (mod_name,'fdrs_uwc', axes(1:3), Time, &
         'Detrainment rate for sat. air from uw_conv', '1/Pa', missing_value=mv)
    id_cqa_uwc = register_diag_field ( mod_name, 'cqa_uwc', axes(1:3), Time, &
         'Updraft fraction from uw_conv', 'none', missing_value=mv)
    id_cql_uwc = register_diag_field ( mod_name, 'cql_uwc', axes(1:3), Time, &
         'Updraft liquid from uw_conv', 'kg/kg', missing_value=mv)
    id_cqi_uwc = register_diag_field ( mod_name, 'cqi_uwc', axes(1:3), Time, &
         'Updraft ice from uw_conv', 'kg/kg', missing_value=mv)
    id_cqn_uwc = register_diag_field ( mod_name, 'cqn_uwc', axes(1:3), Time, &
         'Updraft liquid drop from uw_conv', '/kg', missing_value=mv)
    id_hlflx_uwc=register_diag_field (mod_name,'hlflx_uwc',axes(1:3),Time, &
         'Liq.wat.pot.temp. flux from uw_conv', 'W/m2', missing_value=mv)
    id_qtflx_uwc = register_diag_field (mod_name,'qtflx_uwc',axes(1:3),Time, &
         'Total water flux from uw_conv', 'W/m2', missing_value=mv)
    id_prec_uwc = register_diag_field (mod_name,'prec_uwc', axes(1:2), Time, &
         'Precipitation rate from uw_conv', 'mm/day' )
    id_snow_uwc = register_diag_field (mod_name,'snow_uwc', axes(1:2), Time, &
         'Frozen precip. rate from uw_conv', 'mm/day' )
    id_cin_uwc = register_diag_field ( mod_name, 'cin_uwc', axes(1:2), Time, &
         'CIN from uw_conv', 'm2/s2' )
    id_cape_uwc= register_diag_field ( mod_name,'cape_uwc', axes(1:2), Time, &
         'CAPE from uw_conv', 'm2/s2' )
    id_cbmf_uwc = register_diag_field (mod_name,'cbmf_uwc', axes(1:2), Time, &
         'Cloud-base mass flux from uw_conv', 'kg/m2/s' )
    id_wrel_uwc = register_diag_field (mod_name,'wrel_uwc', axes(1:2), Time, &
         'Release level vertical velocity from uw_conv', 'm/s' )
    id_ufrc_uwc = register_diag_field (mod_name,'ufrc_uwc', axes(1:2), Time, &
         'Release level updraft fraction from uw_conv', 'none' )
    id_tke_uwc = register_diag_field ( mod_name, 'tke_uwc', axes(1:2), Time, &
         'PBL mean TKE from uw_conv', 'm2/s2' )
    id_plcl_uwc = register_diag_field (mod_name,'plcl_uwc', axes(1:2), Time, &
         'LCL pressure from uw_conv', 'hPa' )
    id_plfc_uwc = register_diag_field (mod_name,'plfc_uwc', axes(1:2), Time, &
         'LFC pressure from uw_conv', 'hPa' )
    id_plnb_uwc = register_diag_field (mod_name,'plnb_uwc', axes(1:2), Time, &
         'LNB pressure from uw_conv', 'hPa' )
    id_zinv_uwc = register_diag_field (mod_name,'zinv_uwc', axes(1:2), Time, &
         'Inversion pressure from uw_conv', 'm' )
    id_pct_uwc = register_diag_field ( mod_name, 'pct_uwc', axes(1:2), Time, &
         'Cloud-top pressure from uw_conv', 'hPa' )
    id_pcb_uwc = register_diag_field ( mod_name, 'pcb_uwc', axes(1:2), Time, &
         'Cloud-base pressure from uw_conv', 'hPa' )
    id_cush_uwc = register_diag_field (mod_name,'cush_uwc', axes(1:2), Time, &
         'Convective scale height from uw_conv', 'm' )
    id_dcin_uwc = register_diag_field (mod_name, 'dcin_uwc', axes(1:2), Time, &
         'dCIN/cbmf from uw_conv', 'm2/s2/(kg/m2/s)' )
    id_dcape_uwc= register_diag_field (mod_name, 'dcape_uwc', axes(1:2), Time, &
         'dCAPE/cbmf from uw_conv', 'm2/s2/(kg/m2/s)' )
    id_dwfn_uwc = register_diag_field (mod_name, 'dwfn_uwc',  axes(1:2), Time, &
         'dwfn/cbmf from uw_conv', '(m2/s2)/(kg/m2/s)' )
    id_enth_uwc = register_diag_field (mod_name,'enth_uwc', axes(1:2), Time, &
         'Column-integrated enthalpy tendency from uw_conv', 'W/m2' )
    id_qtmp_uwc = register_diag_field (mod_name,'qtmp_uwc', axes(1:2), Time, &
         'Column-integrated water tendency from uw_conv', 'kg/m2/s' )
    id_ocode_uwc = register_diag_field (mod_name,'ocode_uwc', axes(1:2), Time, &
         'Out code from uw_conv', 'none' )
    if ( do_strat ) then
       id_qldt_uwc= register_diag_field (mod_name,'qldt_uwc',axes(1:3),Time, &
            'Liquid water tendency from uw_conv', 'kg/kg/day', missing_value=mv)
       id_qidt_uwc= register_diag_field (mod_name,'qidt_uwc',axes(1:3),Time, &
            'Ice water tendency from uw_conv', 'kg/kg/day', missing_value=mv)
       id_qadt_uwc= register_diag_field (mod_name,'qadt_uwc',axes(1:3),Time, &
            'CLD fraction tendency from uw_conv', '1/day', missing_value=mv )
       id_qndt_uwc= register_diag_field (mod_name,'qndt_uwc',axes(1:3),Time, &
            'Cloud droplet number fraction tendency from uw_conv', '#/kg/day', missing_value=mv )
    end if
    if ( ntracers>0 ) then
      allocate(id_tracerdt_uwc(ntracers), id_tracerdt_uwc_col(ntracers) )
       allocate(id_tracerdtwet_uwc(ntracers), id_tracerdtwet_uwc_col(ntracers))
      do nn = 1,ntracers
         id_tracerdt_uwc(nn) = &
            register_diag_field (mod_name, trim(cpn%tracername(nn))//'d t_uwc', &
                                    axes(1:3), Time, &
                                  trim(cpn%tracername(nn)) //' tendency from uw_conv', &
                                  trim(cpn%tracer_units(nn))//'/s', missing_value=mv)
            id_tracerdt_uwc_col(nn) = &
              register_diag_field (mod_name, trim(cpn%tracername(nn))//'dt_uwc_col', &
                                     axes(1:2), Time, &
                                   trim(cpn%tracername(nn)) //' column tendency from uw_conv', &
                                   trim(cpn%tracer_units(nn))//'*(kg/m2)/s', missing_value=mv)
           id_tracerdtwet_uwc(nn) = &
              register_diag_field (mod_name, trim(cpn%tracername(nn))//'dt_uwc_wet', &
                                    axes(1:3), Time, &
                                   trim(cpn%tracername(nn)) //' tendency from uw_conv wetdep', &
                                   trim(cpn%tracer_units(nn))//'/s', missing_value=mv)
            id_tracerdtwet_uwc_col(nn) = &
              register_diag_field (mod_name, trim(cpn%tracername(nn))//'dt_uwc_wet_col', &
                                   axes(1:2), Time, &
                                   trim(cpn%tracername(nn)) //' column tendency from uw_conv wetdep', &
                                   trim(cpn%tracer_units(nn))//'*(kg/m2)/s', missing_value=mv)
        end do
     end if

    module_is_initialized = .true.
    
  end SUBROUTINE UW_CONV_INIT

!#####################################################################
!#####################################################################

  subroutine uw_conv_end
    call exn_end_k
    call findt_end_k
    call sd_end_k(sd)
    call ac_end_k(ac)
    call cp_end_k(cp)
    call cp_end_k(cp1)
    call ct_end_k(ct)
    call ct_end_k(ct1)
    module_is_initialized = .FALSE.
  end subroutine uw_conv_end

!#####################################################################
!#####################################################################

  SUBROUTINE uw_conv(is, js, Time, tb, qv, ub, vb, pmid, pint,zmid,  & !input
       zint, q, omega, delt, pblht, ustar, bstar, qstar, land, coldT,& !input
       asol,                                                         & !input
       cush, do_strat,                                               & !input
       tten, qvten, qlten, qiten, qaten, qnten,                      & !output
       uten, vten, rain, snow,                                       & !output
       cmf, hlflx, qtflx, pflx, cldql, cldqi, cldqa, cldqn, cbmfo,  & !output
        cbmf_clo,  &
        tracers, trtend)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     SHALLOW CONVECTION SCHEME
!     Described in Bretherton et. al (MWR, April 2004)
!     For info contact Ming Zhao: ming.zhao@noaa.gov
!
!     Inputs: see below
!
!     Outputs: see below
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    implicit none

    type(time_type), intent(in)  :: Time
    integer,         intent(in)  :: is, js
    real,            intent(in)  :: delt 

    real, intent(in), dimension(:,:,:)   :: ub,vb !wind profile (m/s)
    real, intent(in), dimension(:,:,:)   :: zint  !height@model interfaces(m)
    real, intent(in), dimension(:,:,:)   :: pint  !pressure@model interfaces(pa)
    real, intent(in), dimension(:,:,:)   :: tb    !temperature profile (K)
    real, intent(in), dimension(:,:,:)   :: qv    !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:,:) :: q     !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:)   :: pmid  !pressure@model mid-levels (pa)
    real, intent(in), dimension(:,:,:)   :: zmid  !height@model mid-levels (m)
    real, intent(in), dimension(:,:,:)   :: omega !omega (Pa/s)
    real, intent(in), dimension(:,:)     :: land  !land fraction
    logical,intent(in)                   :: do_strat !logical flag
    logical,intent(in), dimension(:,:)   :: coldT    !logical flag

    real, intent(in),    dimension(:,:)  :: pblht, ustar, bstar, qstar !pbl height...
    real, intent(inout), dimension(:,:)  :: cush  ! convective scale height (m) 

    type(aerosol_type),  intent (in)     :: asol
   
    real, intent(out), dimension(:,:,:)  :: tten,qvten              ! T,qv tendencies
    real, intent(out), dimension(:,:,:)  :: qlten,qiten,qaten,qnten ! q tendencies
    real, intent(out), dimension(:,:,:)  :: uten,vten               ! u,v tendencies
   
    real, intent(out), dimension(:,:,:)  :: cldql,cldqi,cldqa, cldqn!in-updraft q
    real, intent(out), dimension(:,:,:)  :: cmf    ! mass flux at level above layer (kg/m2/s)
    real, intent(out), dimension(:,:,:)  :: pflx   ! precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: hlflx ! theta_l flux
    real, intent(out), dimension(:,:,:)  :: qtflx  ! qt  flux
    real, intent(out), dimension(:,:)    :: rain, snow
    real, intent(in), dimension(:,:)  :: cbmf_clo  ! cloud-base mass flux
    real, intent(inout), dimension(:,:)  :: cbmfo  ! cloud-base mass flux
!   real, intent(in), dimension(:,:)  :: wrelo, ufrco  ! cloud-base mass flux
    real, intent(in),  dimension(:,:,:,:)  :: tracers         ! env. tracers
    real, intent(out), dimension(:,:,:,:)  :: trtend          ! calculated tracer tendencies

    integer i, j, k, kl, klm, km1, nk, ks, m, naer, na, n

    real rhos0j, thj, qvj, qlj, qij, qse, thvj
    real thvuinv, hlsrc, thcsrc, qctsrc, plcltmp, tmp
    real zsrc, psrc, cbmftmp, cbmf_old, rkm_sh1

    real, dimension(size(tb,1),size(tb,2)) :: &
         plcl,       &     ! pressure of lifting condensation level (Pa)
         plfc,       &     ! pressure of level of free convection (Pa)
         plnb,       &     ! pressure of level of neutral buoyancy (Pa)
         cino,       &     ! cin (m2/s2)
         capeo,      &     ! cape(m2/s2)
         tkeo,       &     ! tke (m2/s2)
         wrelo,      &     ! release level vertical velocity (m/s) ! cjg wrel_mod bug fix
         ufrco,      &     ! cloud-base updraft fraction           ! cjg wrel_mod bug fix
         zinvo,      &     ! surface driven mixed-layer height
         denth,      &     
         dqtmp,      &
         dcino,      &     ! dcin (m2/s2)
         dcapeo,     &     ! dcape(m2/s2)
         dwfno,      &     ! dwfn(m2/s2)
         ocode

    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: wuo,fero,fdro,fdrso

    real, dimension(size(tb,3)) :: am1, am2, am3, qntmp
    
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: pmass    ! layer mass (kg/m2)
    real, dimension(size(tb,1),size(tb,2))            :: tempdiag ! temporary diagnostic variable
    real, dimension(size(tracers,1),size(tracers,2),size(tracers,3),size(tracers,4))  :: trwet          ! calculated tracer wet deposition tendencies

    integer imax, jmax, kmax
    
    logical used

    imax  = size( tb, 1 )
    jmax  = size( tb, 2 )
    kmax  = size( tb, 3 )
    sd % kmax=kmax

    kl=kmax-1
    klm=kl-1

   !initialize 3D variables outside the loop

    tten=0.; qvten=0.; qlten=0.; qiten=0.; qaten=0.; qnten=0.;
    uten=0.; vten =0.; rain =0.; snow =0.; plcl =0.; plfc=0.; plnb=0.;  
    cldqa=0.; cldql=0.; cldqi=0.; cldqn=0.;
    hlflx=0.; qtflx=0.; pflx=0.; am1=0.; am2=0.; am3=0.;

    cino=0.; capeo=0.; tkeo=0.;                     zinvo=0.; wuo=0.; 
    fero=0.; fdro=0.; fdrso=0.; cmf=0.; denth=0.;  dqtmp=0.; ocode=0;
    dcapeo=0.; dcino=0.; dwfno=0.;
    trtend=0.
    trwet = 0.

    naer = size(asol%aerosol,4)

    do j = 1, jmax
      do i=1, imax

        call clearit(ac, cc, cp, ct, cp1, ct1);

        if (.not. do_unified_convective_closure) then

          !relaxation TKE back to 0 with time-scale of disscale
          !tkeavg = ustar(i,j)*bstar(i,j)*disscale 
          !dissipate tke with length-scale of disscale
          !tkeavg=(ustar(i,j)*bstar(i,j)*disscale)**(2./3.)
          !below following Holtslag and Boville 1993
          tkeo(i,j) = (ustar(i,j)**3.+0.6*ustar(i,j)*bstar(i,j)*pblht(i,j))**(2./3.)
          tkeo(i,j) = max(tkemin,tkeo(i,j))

          cc%scaleh = cush(i,j); 
          cush(i,j) = -1.;
          if(cc%scaleh.le.0.0) cc%scaleh=1000.

          am1(:) = 0.; am2(:) = 0.; am3(:) = 0.;

          do k=1,kmax
            pmass(i,j,k) = (pint(i,j,k+1) - pint(i,j,k))/GRAV
            tmp=1. / (zint(i,j,k)-zint(i,j,k+1)) * 1.0e9 * 1.0e-12
            if(use_online_aerosol) then
              do na = 1,naer
                if(asol%aerosol_names(na) == 'so4' .or. &
                asol%aerosol_names(na) == 'so4_anthro' .or. asol%aerosol_names(na) == 'so4_natural') then
                         am1(k)=am1(k)+asol%aerosol(i,j,k,na)*tmp
                else if(asol%aerosol_names(na) == 'omphilic' .or. &
                 asol%aerosol_names(na) == 'omphobic') then
                       am3(k)=am3(k)+asol%aerosol(i,j,k,na)*tmp
                else if(asol%aerosol_names(na) == 'seasalt1' .or. &
                asol%aerosol_names(na) == 'seasalt2') then
                       am2(k)=am2(k)+asol%aerosol(i,j,k,na)*tmp
                end if
              end do
            else
              am1(k)=(asol%aerosol(i,j,k,1)+asol%aerosol(i,j,k,2))*tmp
              am2(k)= asol%aerosol(i,j,k,5)*tmp
              am3(k)= asol%aerosol(i,j,k,3)*tmp
            endif
          end do

!========Pack column properties into a sounding structure====================

          if (do_qn) then
            qntmp(:)=q(i,j,:,nqn)
          else
            qntmp(:)=0.
          end if

          call pack_sd_k(land(i,j), coldT(i,j), delt, pmid(i,j,:), pint(i,j,:),     &
               zmid(i,j,:), zint(i,j,:), ub(i,j,:), vb(i,j,:), tb(i,j,:),   &
               q(i,j,:,nqv), q(i,j,:,nql), q(i,j,:,nqi), q(i,j,:,nqa), qntmp,       &
               am1(:), am2(:), am3(:), tracers(i,j,:,:), sd, Uw_p)

!========Finite volume intepolation==========================================

          call extend_sd_k(sd,  pblht(i,j),do_ice, Uw_p)
          zinvo (i,j) = sd%zinv


!========Find source air, and do adiabatic cloud lifting======================

          zsrc  =sd%zs (1)
          psrc  =sd%ps (1)
          thcsrc=sd%thc(1)
          qctsrc=sd%qct(1)
          hlsrc =sd%hl (1)
          rkm_sh1=rkm_sh
          if (do_lands) then
            rkm_sh1 = rkm_sh   * (1. - 0.5 * sd%land)
            rkm_sh1 = rkm_sh   * ( 500. / max(pblht(i,j), 500.))
            qctsrc = qt_parcel_k (sd%qct(1), sd%qs(1), qstar(i,j), &
                                  tkeo(i,j), sd%land, gama)
          endif
          
          call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, do_fast, do_ice, ac)
          ac % usrc = sd%u(sd%ktoppbl)
          ac % vsrc = sd%v(sd%ktoppbl)
          plcl (i,j) = ac%plcl
          plfc (i,j) = ac%plfc
          plnb (i,j) = ac%plnb
          cino (i,j) = ac%cin
          capeo(i,j) = ac%cape
          if (do_fast) then
            if (ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.) then
              ocode(i,j)=1; goto 100 !cycle;
            end if
            if (ac%plfc.lt.500.) then
              ocode(i,j)=2; goto 100 !cycle;
            end if
          end if

        else  !(.not. do_unified_convective_closure)
          cush(i,j) = -1.0
          if (.not. exit_flag(i,j)) then
            call sd_retrieve_k (sd, sdij(i,j))
            call ac_retrieve_k (ac, acij(i,j))
            zinvo (i,j) = sd%zinv
            plcl (i,j) = ac%plcl
            plfc (i,j) = ac%plfc
            plnb (i,j) = ac%plnb
            cino (i,j) = ac%cin
            capeo(i,j) = ac%cape
            rkm_sh1 = rkm_sh
            if (do_lands) then
              rkm_sh1 = rkm_sh   * (1. - 0.5 * sd%land)
              rkm_sh1 = rkm_sh   * ( 500. / max(pblht(i,j), 500.))
            endif
          else
            ocode(i,j) = 7
            go to 100
          endif
        endif  !(.not. do_unified_convective_closure)


!========Cumulus closure to determine cloud base mass flux===================


        cbmf_old=cbmfo(i,j) 

        if (do_unified_convective_closure) then
          if (cbmf_clo(i,j) == cbmfo_unmodified(i,j)) then
            cbmfo(i,j) = cbmfo_mod(i,j)
            cc%cbmf = cbmfo_mod(i,j)
          else if (cbmf_clo(i,j) > cbmfo_mod(i,j)) then
            cbmfo(i,j) = cbmfo_mod(i,j)
            cc%cbmf = cbmfo_mod(i,j)
          else
            cbmfo(i,j) = cbmf_clo(i,j)
            cc%cbmf = cbmf_clo(i,j)
          endif


          cc%wrel = wrel_mod(i,j) 
          cc%ufrc = ufrc_mod(i,j)
          cc%scaleh = scaleh_mod(i,j)
          cc%dcin = dcin_mod(i,j)
          cc%dcape = dcape_mod(i,j)
          cc%dwfn = dwfn_mod(i,j)
          cc%wfn = wfn_mod(i,j)

!--> cjg wrel_mod bug fix
          wrelo(i,j) = cc%wrel
          ufrco(i,j) = cc%ufrc
!<-- cjg
        else

          cc%cbmf=cbmf_old;

          if (iclosure.eq.0) then
             call cclosure_bretherton(tkeo(i,j), cpn, sd, Uw_p, ac, cc)
          elseif (iclosure.eq.1) then
             call cclosure_implicit (tkeo(i,j), cpn, sd, Uw_p, ac, cc, &
                                     delt, rkm_sh1, &
                  do_coldT, sd1, ac1, cc1, cp, ct) 
          elseif (iclosure.eq.2) then
             call cclosure_relaxwfn(tkeo(i,j), cpn, sd, Uw_p, ac, cc, cp, ct, delt,  &
                  rkm_sh1, do_coldT, sd1, ac1, cc1, cp1, ct1)
          end if
          

          cbmfo(i,j) = cc%cbmf
!--> cjg wrel_mod bug fix
          wrelo(i,j) = cc%wrel
          ufrco(i,j) = cc%ufrc
!<-- cjg
        endif  ! (do_unified_conv)

        if (.not.do_fast) then
          if (ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.) then
            ocode(i,j)=1; goto 100 !cycle;
          end if
          if (ac%plfc.lt.500.) then
            ocode(i,j)=2; goto 100 !cycle;
          end if
        end if

        if(cc%cbmf.lt.1.e-6 .or. cc%wrel.eq.0.) then
          ocode(i,j)=3; goto 100 !cycle;
        end if

!========Do shallow cumulus plume calculation================================

        cbmftmp=cc%cbmf * cbmf_sh_frac
        call cumulus_plume_k (cpn, sd, ac, cp, rkm_sh1, cbmftmp, &
                              cc%wrel, cc%scaleh, Uw_p)
        if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
          ocode(i,j)=4; goto 100 !cycle;
        end if
        if(cp%cldhgt.ge.cldhgt_max) then
          ocode(i,j)=5; goto 100 !cycle;
        end if
        cush(i,j)=cp%cush

!========Calculate cumulus produced tendencies===============================

        call cumulus_tend_k(cpn, sd, Uw_p, cp, ct, do_coldT)

!========Unpack convective tendencies========================================
          do k = 1,cp%ltop
             nk = kmax+1-k
             uten  (i,j,nk) = ct%uten (k)
             vten  (i,j,nk) = ct%vten (k)
             qlten (i,j,nk) = ct%qlten(k)
             qiten (i,j,nk) = ct%qiten(k)
             qaten (i,j,nk) = ct%qaten(k)
             qnten (i,j,nk) = ct%qnten(k)
             qvten (i,j,nk) = ct%qvten(k)
             pflx  (i,j,nk) = ct%pflx (k)
             tten  (i,j,nk) = ct%tten (k)
             rhos0j = sd%ps(k)/(rdgas*0.5*(cp%thvbot(k+1)+cp%thvtop(k))*sd%exners(k))
             hlflx(i,j,nk) = ct%hlflx(k)
             qtflx (i,j,nk) = rhos0j*HLv*ct%qctflx(k)
             
             cldqa (i,j,nk) = cp%ufrc(k)
             cldql (i,j,nk) = cp%qlu(k)
             cldqi (i,j,nk) = cp%qiu(k)
             cldqn (i,j,nk) = cp%qnu(k)
             cmf   (i,j,nk) = cp%umf(k)
             wuo   (i,j,nk) = cp%wu (k)
             fero  (i,j,nk) = cp%fer(k)
             fdro  (i,j,nk) = cp%fdr(k)
             fdrso (i,j,nk) = cp%fdrsat(k)*cp%umf(k)
             
!            do n = 1, size(trtend,4)
!              trtend(i,j,nk,n) = ct%trten(k,n) + ct%trwet(k,n)
!              trwet(i,j,nk,n)  = ct%trwet(k,n)
!            enddo
          enddo

! make sure the predicted tracer tendencies do not produce negative
! tracers due to convective tendencies. if necessary, adjust the 
! tendencies.
          call check_tracer_realizability (kmax, size(trtend,4), delt, &
                                           cp%tr, ct%trten, ct%trwet) 

          do k = 1,cp%ltop
             nk = kmax+1-k
             do n = 1, size(trtend,4)
               trtend(i,j,nk,n) = ct%trten(k,n) + ct%trwet(k,n)
               trwet(i,j,nk,n)  = ct%trwet(k,n)
             enddo
          enddo
          snow  (i,j)  = ct%snow
          rain  (i,j)  = ct%rain
          denth (i,j)  = ct%denth
          dqtmp (i,j)  = ct%dqtmp

!========Option for deep convection=======================================
100       if (do_deep) then
             cc%cbmf=cbmf_old
             cpn%do_forcedlifting=.true.
             call  deepconv(cpn, sd, Uw_p, ac, cc, cp, ct, delt, do_coldT, &
                  sd1, ac1, cc1, cp1, ct1, cush(i,j), ocode(i,j))
             cpn%do_forcedlifting=.false.
             if(ocode(i,j).eq.6) cycle;

             do k = 1, cp1%ltop
                nk = kmax+1-k
                uten  (i,j,nk) = uten  (i,j,nk) + ct1%uten (k)
                vten  (i,j,nk) = vten  (i,j,nk) + ct1%vten (k)
                qlten (i,j,nk) = qlten (i,j,nk) + ct1%qlten(k)
                qiten (i,j,nk) = qiten (i,j,nk) + ct1%qiten(k)
                qaten (i,j,nk) = qaten (i,j,nk) + ct1%qaten(k) 
                qnten (i,j,nk) = qnten (i,j,nk) + ct1%qnten(k) 
                qvten (i,j,nk) = qvten (i,j,nk) + ct1%qvten(k)
                pflx  (i,j,nk) = pflx  (i,j,nk) + ct1%pflx (k)
                tten  (i,j,nk) = tten  (i,j,nk) + ct1%tten (k)
                rhos0j = sd%ps(k)/(rdgas*0.5*(cp1%thvbot(k+1)+cp1%thvtop(k))*sd%exners(k))
                hlflx(i,j,nk) = hlflx(i,j,nk) + ct1%hlflx(k)
                qtflx (i,j,nk) = qtflx (i,j,nk) + rhos0j*HLv*ct1%qctflx(k)
                cmf   (i,j,nk) = cmf   (i,j,nk) +  cp1%umf(k)
                wuo   (i,j,nk) = wuo   (i,j,nk) +  cp1%wu (k)
                fero  (i,j,nk) = fero  (i,j,nk) +  cp1%fer(k)
                fdro  (i,j,nk) = fdro  (i,j,nk) +  cp1%fdr(k) 
                fdrso (i,j,nk) = fdrso (i,j,nk) + cp1%fdrsat(k)*cp1%umf(k)
!               do n = 1, size(trtend,4)
!                 trtend(i,j,nk,n) = trtend(i,j,nk,n) + ct1%trten(k,n)
!               enddo
             enddo

! make sure the predicted tracer tendencies do not produce negative
! tracers due to convective tendencies. if necessary, adjust the 
! tendencies.
          call check_tracer_realizability (kmax, size(trtend,4), delt, &
                                           cp%tr, ct1%trten, ct1%trwet) 

          do k = 1,cp%ltop
             nk = kmax+1-k
             do n = 1, size(trtend,4)
               trtend(i,j,nk,n) = trtend(i,j,nk,n) + ct1%trten(k,n) + &
                                                         ct1%trwet(k,n)
               trwet(i,j,nk,n)  = trwet(i,j,nk,n) + ct1%trwet(k,n)
             enddo
          enddo
             snow  (i,j)  = snow  (i,j) + ct1%snow
             rain  (i,j)  = rain  (i,j) + ct1%rain
             denth (i,j)  = denth (i,j) + ct1%denth
             dqtmp (i,j)  = dqtmp (i,j) + ct1%dqtmp
             cbmfo (i,j)  = cc%cbmf
             dcapeo(i,j)  = cc%dcape
             dwfno (i,j)  = cc%dwfn
          end if

       enddo
    enddo

    if (.not.do_uwcmt) then
       uten=0.;
       vten=0.;
    end if

    if ( prevent_unreasonable ) then
       where ((q(:,:,:,nqa) + qaten*delt) .lt. 0.)
          qaten(:,:,:) = -1.*q(:,:,:,nqa)/delt
       end where
       where ((q(:,:,:,nqa) + qaten*delt) .gt. 1.)
          qaten(:,:,:)= (1. - q(:,:,:,nqa))/delt
       end where
 
       where ((q(:,:,:,nql) + qlten*delt) .lt. 0.)
          tten (:,:,:) = tten(:,:,:) -(q(:,:,:,nql)/delt+qlten(:,:,:))*HLv/Cp_Air
!-->cjg
          qvten(:,:,:) = qvten(:,:,:)+(q(:,:,:,nql)/delt+qlten(:,:,:))
!<--cjg

          qlten(:,:,:) = qlten(:,:,:)-(q(:,:,:,nql)/delt+qlten(:,:,:))
      end where
 
      where ((q(:,:,:,nqi) + qiten*delt) .lt. 0.)
          tten (:,:,:) = tten(:,:,:) -(q(:,:,:,nqi)/delt+qiten(:,:,:))*HLs/Cp_Air
!-->cjg
          qvten(:,:,:) = qvten(:,:,:)+(q(:,:,:,nqi)/delt+qiten(:,:,:))
!<--cjg
    
         qiten(:,:,:) = qiten(:,:,:)-(q(:,:,:,nqi)/delt+qiten(:,:,:))
      end where

      if (do_qn) then
      where ((q(:,:,:,nqn) + qnten*delt) .lt. 0.)
         qnten(:,:,:) = qnten(:,:,:)-(q(:,:,:,nqn)/delt+qnten(:,:,:))
      end where
      endif

!     where ((tracers(:,:,:,:) + trtend(:,:,:,:)*delt) .lt. 0.)
!        trtend(:,:,:,:) = -tracers(:,:,:,:)/delt
       
!     end where
 
    endif


    !diagnostic output
    if ( id_tdt_uwc   > 0 ) &
         used = send_data( id_tdt_uwc,    tten*86400. , Time, is, js, 1)
    if ( id_qdt_uwc   > 0 ) &
         used = send_data( id_qdt_uwc,    qvten*86400., Time, is, js, 1)
    if ( id_cmf_uwc   > 0 ) &
         used = send_data( id_cmf_uwc,    cmf,          Time, is, js, 1)
    if ( id_wu_uwc    > 0 ) &
         used = send_data( id_wu_uwc,     wuo,          Time, is, js, 1)
    if ( id_fer_uwc   > 0 ) &
         used = send_data( id_fer_uwc,    fero,         Time, is, js, 1)
    if ( id_fdr_uwc   > 0 ) &
         used = send_data( id_fdr_uwc,    fdro,         Time, is, js, 1)
    if ( id_fdrs_uwc  > 0 ) &
         used = send_data( id_fdrs_uwc,   fdrso,        Time, is, js, 1)
    if ( id_cqa_uwc   > 0 ) &
         used = send_data( id_cqa_uwc,    cldqa,        Time, is, js, 1)
    if ( id_cql_uwc   > 0 ) &
         used = send_data( id_cql_uwc,    cldql,        Time, is, js, 1)
    if ( id_cqi_uwc   > 0 ) &
         used = send_data( id_cqi_uwc,    cldqi,        Time, is, js, 1)
    if ( id_cqn_uwc   > 0 ) &
         used = send_data( id_cqn_uwc,    cldqn,        Time, is, js, 1)
    if ( id_hlflx_uwc> 0 ) &
         used = send_data( id_hlflx_uwc, hlflx,       Time, is, js, 1)
    if ( id_qtflx_uwc > 0 ) &
         used = send_data( id_qtflx_uwc,  qtflx,        Time, is, js, 1)
   
    if ( id_prec_uwc > 0 ) &
         used = send_data( id_prec_uwc, (rain+snow)*86400., Time, is, js )
    if ( id_snow_uwc > 0 ) &
         used = send_data( id_snow_uwc, (snow)*86400.,      Time, is, js )
    if ( id_cin_uwc > 0 )  &
         used = send_data( id_cin_uwc,  (cino),             Time, is, js )
    if ( id_cape_uwc > 0 )  &
         used = send_data( id_cape_uwc, (capeo),            Time, is, js )
    if ( id_tke_uwc > 0 )  &
         used = send_data( id_tke_uwc,  (tkeo),             Time, is, js )
    if ( id_cbmf_uwc > 0 ) &
         used = send_data( id_cbmf_uwc, (cbmfo),            Time, is, js )
!-->cjg wrel_mod bug fix
!    if ( id_wrel_uwc > 0 ) &
!         used = send_data( id_wrel_uwc, (wrel_mod),            Time, is, js )
!    if ( id_ufrc_uwc > 0 ) &
!         used = send_data( id_ufrc_uwc, (ufrc_mod),            Time, is, js )

    if ( id_wrel_uwc > 0 ) &
         used = send_data( id_wrel_uwc, (wrelo),            Time, is, js )
    if ( id_ufrc_uwc > 0 ) &
         used = send_data( id_ufrc_uwc, (ufrco),            Time, is, js )
!<--cjg
    if ( id_plcl_uwc > 0 ) &
         used = send_data( id_plcl_uwc, (plcl*0.01),        Time, is, js )
    if ( id_plfc_uwc > 0 ) &
         used = send_data( id_plfc_uwc, (plfc*0.01),        Time, is, js )
    if ( id_plnb_uwc > 0 ) &
         used = send_data( id_plnb_uwc, (plnb*0.01),        Time, is, js )
    if ( id_zinv_uwc > 0 ) &
         used = send_data( id_zinv_uwc, (zinvo),            Time, is, js )
    if ( id_cush_uwc > 0 ) &
         used = send_data( id_cush_uwc, (cush),             Time, is, js )
    if ( id_dcin_uwc > 0 ) &
         used = send_data( id_dcin_uwc, (dcino),            Time, is, js )
    if ( id_dcape_uwc> 0 ) &
         used = send_data( id_dcape_uwc,(dcapeo),           Time, is, js )
    if ( id_dwfn_uwc > 0 ) &
         used = send_data( id_dwfn_uwc, (dwfno),            Time, is, js )
    if ( id_enth_uwc > 0 ) &
         used = send_data( id_enth_uwc, (denth),            Time, is, js )
    if ( id_qtmp_uwc > 0 ) &
         used = send_data( id_qtmp_uwc, (dqtmp),            Time, is, js )
    if ( id_ocode_uwc > 0 ) &
         used = send_data( id_ocode_uwc,(ocode),            Time, is, js )
   
    if ( do_strat ) then
       if ( id_qldt_uwc > 0 ) &
            used = send_data( id_qldt_uwc, qlten*86400.,    Time, is, js, 1)
       if ( id_qidt_uwc > 0 ) &
            used = send_data( id_qidt_uwc, qiten*86400.,    Time, is, js, 1)
       if ( id_qadt_uwc > 0 ) &
            used = send_data( id_qadt_uwc, qaten*86400.,    Time, is, js, 1)
       if ( id_qndt_uwc > 0 ) &
            used = send_data( id_qndt_uwc, qnten*86400.,    Time, is, js, 1)
    end if

    if ( allocated(id_tracerdt_uwc) ) then
       do n = 1,size(id_tracerdt_uwc)
          if ( id_tracerdt_uwc(n) > 0 ) &
            used = send_data( id_tracerdt_uwc(n), trtend(:,:,:,n), Time, is, js, 1)
       end do
    end if
    if ( allocated(id_tracerdt_uwc_col) ) then
       do n = 1,size(id_tracerdt_uwc_col)
          if ( id_tracerdt_uwc_col(n) > 0 ) then
            tempdiag = 0.
            do k = 1,kmax
               tempdiag(:,:) = tempdiag(:,:) + trtend(:,:,k,n) * pmass(:,:,k)
            end do
            used = send_data( id_tracerdt_uwc_col(n), tempdiag(:,:), Time, is, js)
          end if
       end do
    end if
    if ( allocated(id_tracerdtwet_uwc) ) then
       do n = 1,size(id_tracerdtwet_uwc)
          if ( id_tracerdtwet_uwc(n) > 0 ) &
            used = send_data( id_tracerdtwet_uwc(n), trwet(:,:,:,n), Time, is, js, 1)
       end do
    end if
    if ( allocated(id_tracerdtwet_uwc_col) ) then
       do n = 1,size(id_tracerdtwet_uwc_col)
          if ( id_tracerdtwet_uwc_col(n) > 0 ) then
             tempdiag = 0.
             do k = 1,kmax
               tempdiag(:,:) = tempdiag(:,:) + trwet(:,:,k,n) * pmass(:,:,k)
            end do
            used = send_data( id_tracerdtwet_uwc_col(n), tempdiag(:,:), Time, is, js)
          end if
       end do
    end if

    if (.not.apply_tendency) then
       uten=0.; vten=0.; tten=0.; qvten=0.; cmf=0.; rain=0.; snow=0.;
       qlten=0.; qiten=0.; qaten=0.; qnten=0.;
    end if

    if (do_unified_convective_closure) then
      do j=1,jmax
       do i=1,imax
         if (associated (acij(i,j)%t)) then
           call ac_end_k (acij(i,j))
           call sd_end_k (sdij(i,j)) 
         endif
       end do
     end do
    deallocate (cbmfo_mod )
    deallocate (cbmfo_unmodified)
    deallocate (wrel_mod )
    deallocate (ufrc_mod )
    deallocate (scaleh_mod )
    deallocate (dcin_mod )
    deallocate (dcape_mod )
    deallocate (dwfn_mod )
    deallocate (wfn_mod )
    deallocate (acij)
    deallocate (sdij)
    deallocate (exit_flag)
    endif
   
  END SUBROUTINE UW_CONV

!#####################################################################
!#####################################################################

  subroutine clearit(ac, cc, cp, ct, cp1, ct1)

    type(adicloud), intent(inout) :: ac
    type(cclosure), intent(inout) :: cc
    type(cplume),   intent(inout) :: cp,cp1
    type(ctend),    intent(inout) :: ct,ct1

    call ac_clear_k(ac); 
    ac%klcl =0;  ac%klfc =0;  ac%klnb =0; 

    cc%wrel=0.; cc%ufrc=0.; cc%scaleh=0.;

    call cp_clear_k(cp)
    call ct_clear_k(ct);
    call cp_clear_k(cp1);
    call ct_clear_k(ct1);

  end subroutine clearit

!#####################################################################
!#####################################################################

  
  subroutine deepconv(cpn, sd, Uw_p, ac, cc, cp, ct, delt, do_coldT, sd1, ac1, &
       cc1, cp1, ct1, cush, ocode)
    implicit none

    type(cpnlist),  intent(in)    :: cpn
    type(sounding), intent(in)    :: sd
    type(uw_params), intent(inout)    :: Uw_p
    type(adicloud), intent(in)    :: ac
    real,           intent(in)    :: delt
    logical,        intent(in)    :: do_coldT
    type(sounding), intent(inout) :: sd1
    type(adicloud), intent(inout) :: ac1
    type(cclosure), intent(inout) :: cc,cc1
    type(cplume),   intent(inout) :: cp,cp1
    type(ctend),    intent(inout) :: ct,ct1
    real,           intent(inout) :: cush, ocode

    integer :: k
    real    :: rkm, cbmf0=0.0001, delp, cbmf_old, tmp, cbmfs, cbmf_max
    cbmf_old= cc%cbmf

    call cumulus_plume_k(cpn, sd, ac, cp, rkm_dp, cbmf0, cc%wrel, cc%scaleh, Uw_p)
    if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
       cc % dcape=0.; cc % cbmf=0.;
       ocode=6; return
    else
       call cumulus_tend_k(cpn, sd, Uw_p, cp, ct, do_coldT)
       call sd_copy_k(sd, sd1)
       sd1 % t  = sd1 % t  + ct%tten  * delt
       sd1 % qv = sd1 % qv + ct%qvten * delt
       sd1 % ql = sd1 % ql + ct%qlten * delt
       sd1 % qi = sd1 % qi + ct%qiten * delt
       sd1 % qa = sd1 % qa + ct%qaten * delt
       sd1 % qn = sd1 % qn + ct%qnten * delt
       sd1 % u  = sd1 % u  + ct%uten  * delt
       sd1 % v  = sd1 % v  + ct%vten  * delt

       call extend_sd_k(sd1, sd%pblht, do_ice, Uw_p)

       call adi_cloud_k(sd1%zs(1), sd1%ps(1), sd1%hl(1), sd1%thc(1), sd1%qct(1), sd1, Uw_p, do_fast, do_ice, ac1)
       cc % dcape=(ac1%cape-ac%cape)/cbmf0

       call cumulus_plume_k(cpn, sd1, ac1, cp1, rkm_dp, cbmf0, cc%wrel, cc%scaleh, Uw_p)
       cc%dwfn=0.; cc%wfn=0.; delp=0.;
       do k=cp1%krel, cp1%let
          cc % wfn  = cc % wfn  + 0.5*(cp %wu(k)*cp %wu(k)) * cp%dp(k)
          cc % dwfn = cc % dwfn + 0.5*(cp1%wu(k)*cp1%wu(k) - cp%wu(k)*cp%wu(k)) * cp%dp(k)
          delp      = delp + cp%dp(k)
       end do
       cc % wfn  = cc % wfn  / delp 
       cc % dwfn = cc % dwfn / delp / cbmf0

       if (do_relaxcape) then
          cbmfs = - ac%cape / cc % dcape / (cc%tau_dp/delt)
        elseif (do_relaxwfn) then
          cbmfs = - cc%wfn  / cc % dwfn  / (cc%tau_dp/delt)
       else
          cbmfs = - cc%wfn  / cc % dwfn
          tmp   = delt/cc%tau_dp
          cbmfs = (cbmf_old+tmp*cbmfs)/(1.+tmp)
       end if


       cbmf_max=(sd%ps(0) - sd%ps(cp1%krel))*(0.25/delt)/Grav
       cc % cbmf = max(min(cbmfs, cbmf_max), 0.)

    end if

    if(cc%cbmf.lt.1.e-6 .or. cc%wrel.eq.0.) then
      cc % dcape=0.; cc % cbmf=0.;
      ocode=6; return
    end if
    call cumulus_plume_k(cpn, sd, ac, cp1, rkm_dp, cc%cbmf, cc%wrel, cc%scaleh, Uw_p)
    cush=cp1%cush; cc%scaleh=cp1%cush;
    call cumulus_tend_k(cpn, sd, Uw_p, cp1, ct1, do_coldT)

    !test
    !call cumulus_downdraft_k(sd, cp, Uw_p)


  end subroutine deepconv

!#####################################################################
!#####################################################################
 
  subroutine deepconv_test(cpn, sd, ac, cc, cp, ct, delt, do_coldT, sd1, ac1, &
       cc1, cp1, ct1, cush, ocode)
    implicit none

    type(cpnlist),  intent(in)    :: cpn
    type(sounding), intent(in)    :: sd
    type(adicloud), intent(in)    :: ac
    real,           intent(in)    :: delt
    logical,        intent(in)    :: do_coldT
    type(sounding), intent(inout) :: sd1
    type(adicloud), intent(inout) :: ac1
    type(cclosure), intent(inout) :: cc,cc1
    type(cplume),   intent(inout) :: cp,cp1
    type(ctend),    intent(inout) :: ct,ct1
    real,           intent(inout) :: cush, ocode
 
    integer :: k
    real    :: rkm, cbmf0=0.0001, delp, cbmf_old, tmp, cbmfs, cbmf_max
    cbmf_old= cc%cbmf
    
    call sd_copy_k(sd, sd1)
    sd1 % qv = sd % qs
    call extend_sd_k(sd1, sd%pblht, do_ice, Uw_p)
   cc%wrel=0.5
    call cumulus_plume_k(cpn, sd1, ac, cp, rkm_dp, cbmf0, cc%wrel, cc%scaleh, Uw_p)
    if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
       cc % dcape=0.; cc % cbmf=0.;
       ocode=6; return
    else
      call cumulus_tend_k(cpn, sd, Uw_p, cp, ct, do_coldT)
    end if
 
   end subroutine deepconv_test

!#####################################################################
!#####################################################################

  SUBROUTINE calculate_uw_closure(is, js, Time, tb, qv, ub, vb, pmid, pint,zmid,  & !input
       zint, q, omega, delt, pblht, ustar, bstar, qstar, land, coldT,& !input
       asol, cush,                                                         & !input
!                                                           cbmfo)   !output
                                                            cbmfo,  & !output
        tracers        )

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     SHALLOW CONVECTION SCHEME
!     Described in Bretherton et. al (MWR, April 2004)
!     For info contact Ming Zhao: ming.zhao@noaa.gov
!
!     Inputs: see below
!
!     Outputs: see below
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    implicit none

    type(time_type), intent(in)  :: Time
    integer,         intent(in)  :: is, js
    real,            intent(in)  :: delt 

    real, intent(in), dimension(:,:,:)   :: ub,vb !wind profile (m/s)
    real, intent(in), dimension(:,:,:)   :: zint  !height@model interfaces(m)
    real, intent(in), dimension(:,:,:)   :: pint  !pressure@model interfaces(pa)
    real, intent(in), dimension(:,:,:)   :: tb    !temperature profile (K)
    real, intent(in), dimension(:,:,:)   :: qv    !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:,:) :: q     !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:)   :: pmid  !pressure@model mid-levels (pa)
    real, intent(in), dimension(:,:,:)   :: zmid  !height@model mid-levels (m)
    real, intent(in), dimension(:,:,:)   :: omega !omega (Pa/s)
    real, intent(in), dimension(:,:)     :: land  !land fraction
    logical,intent(in), dimension(:,:)   :: coldT    !logical flag

    real, intent(in),    dimension(:,:)  :: pblht, ustar, bstar, qstar !pbl height...
    type(aerosol_type),  intent (in)     :: asol
    real, intent(in   ), dimension(:,:)  :: cush  ! convective scale height (m) 

   
   
    real, intent(inout), dimension(:,:)  :: cbmfo  ! cloud-base mass flux
    real, intent(in),  dimension(:,:,:,:)  :: tracers         ! env. tracers

    integer i, j, k, kl, klm, km1, nk, ks, m, naer, na, n

    real rhos0j, thj, qvj, qlj, qij, qse, thvj
    real thvuinv, hlsrc, thcsrc, qctsrc, plcltmp, tmp
    real zsrc, psrc, cbmftmp, cbmf_old, rkm_sh1

    real, dimension(size(tb,1),size(tb,2)) :: &
         plcl,       &     ! pressure of lifting condensation level (Pa)
         plfc,       &     ! pressure of level of free convection (Pa)
         plnb,       &     ! pressure of level of neutral buoyancy (Pa)
         cino,       &     ! cin (m2/s2)
         capeo,      &     ! cape(m2/s2)
         tkeo,       &     ! tke (m2/s2)
         wrelo,      &     ! release level vertical velocity (m/s)
         ufrco,      &     ! cloud-base updraft fraction
         zinvo,      &     ! surface driven mixed-layer height
         denth,      &     
         dqtmp,      &
         dcino,      &     ! dcin (m2/s2)
         dcapeo,     &     ! dcape(m2/s2)
         dwfno             ! dwfn(m2/s2)
!        ocode

    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: wuo,fero,fdro,fdrso

    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: pmass    ! layer mass (kg/m2)
    real, dimension(size(tb,3)) :: am1, am2, am3, qntmp
    real :: cbmf_unmod
    

    integer imax, jmax, kmax, ntr
    
    logical used

    imax  = size( tb, 1 )
    jmax  = size( tb, 2 )
    kmax  = size( tb, 3 )
    ntr = size (tracers, 4)
    sd % kmax=kmax

    naer = size(asol%aerosol,4)

    allocate (cbmfo_mod (imax, jmax))
    allocate (cbmfo_unmodified (imax, jmax))
    allocate (wrel_mod (imax, jmax))
    allocate (ufrc_mod (imax, jmax))
    allocate (scaleh_mod (imax, jmax))
    allocate (dcin_mod (imax, jmax))
    allocate (dcape_mod (imax, jmax))
    allocate (dwfn_mod (imax, jmax))
    allocate (wfn_mod (imax, jmax))
    allocate (acij(imax,jmax))
    allocate (sdij(imax,jmax))
    allocate (exit_flag(imax, jmax))

    exit_flag = .false.
    cbmfo_mod = 0.
    cbmfo_unmodified = 0.
    wrel_mod = 0.
    ufrc_mod = 0.
    scaleh_mod = 0.
    dcin_mod = 0.
    dcape_mod = 0.
    dwfn_mod = 0.
    wfn_mod = 0.


    do j = 1, jmax
      do i=1, imax

        call clearit(ac, cc, cp, ct, cp1, ct1);

          !relaxation TKE back to 0 with time-scale of disscale
          !tkeavg = ustar(i,j)*bstar(i,j)*disscale 
          !dissipate tke with length-scale of disscale
          !tkeavg=(ustar(i,j)*bstar(i,j)*disscale)**(2./3.)
          !below following Holtslag and Boville 1993
        tkeo(i,j) = (ustar(i,j)**3.+0.6*ustar(i,j)*bstar(i,j)*pblht(i,j))**(2./3.)
        tkeo(i,j) = max(tkemin,tkeo(i,j))

        cc%scaleh = cush(i,j); 
        if(cc%scaleh.le.0.0) cc%scaleh=1000.

        am1(:) = 0.; am2(:) = 0.; am3(:) = 0.;

        do k=1,kmax
          pmass(i,j,k) = (pint(i,j,k+1) - pint(i,j,k))/GRAV
          tmp=1. / (zint(i,j,k)-zint(i,j,k+1)) * 1.0e9 * 1.0e-12
          if(use_online_aerosol) then
            do na = 1,naer
              if(asol%aerosol_names(na) == 'so4' .or. &
                asol%aerosol_names(na) == 'so4_anthro' .or. asol%aerosol_names(na) == 'so4_natural') then
                         am1(k)=am1(k)+asol%aerosol(i,j,k,na)*tmp
              else if(asol%aerosol_names(na) == 'omphilic' .or. &
                 asol%aerosol_names(na) == 'omphobic') then
                       am3(k)=am3(k)+asol%aerosol(i,j,k,na)*tmp
              else if(asol%aerosol_names(na) == 'seasalt1' .or. &
                asol%aerosol_names(na) == 'seasalt2') then
                       am2(k)=am2(k)+asol%aerosol(i,j,k,na)*tmp
              end if
            end do

          else
            am1(k)=(asol%aerosol(i,j,k,1)+asol%aerosol(i,j,k,2))*tmp
            am2(k)= asol%aerosol(i,j,k,5)*tmp
            am3(k)= asol%aerosol(i,j,k,3)*tmp
          endif
        end do

!========Pack column properties into a sounding structure====================

        if (do_qn) then
          qntmp(:)=q(i,j,:,nqn)
        else
          qntmp(:)=0.
        end if
          call pack_sd_k(land(i,j), coldT(i,j), delt, pmid(i,j,:), pint(i,j,:),     &
               zmid(i,j,:), zint(i,j,:), ub(i,j,:), vb(i,j,:), tb(i,j,:),   &
               q(i,j,:,nqv), q(i,j,:,nql), q(i,j,:,nqi), q(i,j,:,nqa), qntmp,       &
               am1(:), am2(:), am3(:), tracers(i,j,:,:), sd, Uw_p)

!========Finite volume intepolation==========================================

          call extend_sd_k(sd,  pblht(i,j),do_ice, Uw_p)


!========Find source air, and do adiabatic cloud lifting======================

          zsrc  =sd%zs (1)
          psrc  =sd%ps (1)
          thcsrc=sd%thc(1)
          qctsrc=sd%qct(1)
          hlsrc =sd%hl (1)
          rkm_sh1=rkm_sh
          if (do_lands) then
            rkm_sh1 = rkm_sh   * (1. - 0.5 * sd%land)
            qctsrc = qt_parcel_k (sd%qct(1), sd%qs(1), qstar(i,j), &
                                  tkeo(i,j), sd%land, gama)
          endif

          call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, do_fast, do_ice, ac)
          ac % usrc = sd%u(sd%ktoppbl)
          ac % vsrc = sd%v(sd%ktoppbl)

          if (do_fast) then
            if (ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.) then
              exit_flag(i,j) = .true.
              cycle    !cycle;
            end if
            if (ac%plfc.lt.500.) then
              exit_flag(i,j) = .true.
              cycle    !cycle;
            end if
          end if


!========Cumulus closure to determine cloud base mass flux===================

          cbmf_old=cbmfo(i,j); cc%cbmf=cbmf_old;

          call cclosure_bretherton(tkeo(i,j), cpn, sd, Uw_p, ac, cc, &
                                    cbmf_unmod=cbmf_unmod)

          cbmfo_unmodified(i,j) = cbmf_unmod
          cbmfo_mod(i,j) = cc%cbmf
          cbmfo(i,j) = cbmf_unmod
          wrel_mod(i,j) = cc%wrel
          ufrc_mod(i,j) = cc%ufrc
          scaleh_mod(i,j) = cc%scaleh
          dcin_mod(i,j) = cc%dcin
          dcape_mod(i,j) = cc%dcape
          dwfn_mod(i,j) = cc%dwfn
          wfn_mod(i,j) = cc%wfn 

          call sd_save_k (sd, sdij(i,j))
          call ac_save_k (ac, acij(i,j))
        end do
      end do

   
  END SUBROUTINE calculate_uw_closure

!#####################################################################

end MODULE UW_CONV_MOD
