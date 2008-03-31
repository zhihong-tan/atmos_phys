
MODULE UW_CONV_MOD

  use           mpp_mod, only : mpp_pe, mpp_root_pe, stdlog
  use      Constants_Mod, ONLY: tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas
  use   Diag_Manager_Mod, ONLY: register_diag_field, send_data
  use   Time_Manager_Mod, ONLY: time_type, get_time 
  use           fms_mod, only : write_version_number, open_namelist_file, check_nml_error,&
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
                                   exn_init_k, exn_end_k, findt_init_k,&
                                   findt_end_k, &
                                   check_tracer_realizability, &
                                   qt_parcel_k, &
                                   adicloud, sounding, uw_params

  use  conv_plumes_k_mod,only    : cp_init_k, cp_end_k, cp_clear_k, &
                                   ct_init_k, ct_end_k, ct_clear_k, &
                                   cumulus_tend_k, cumulus_plume_k, &
                                   cplume, ctend, cpnlist

  use  conv_closures_mod,only    : cclosure_bretherton,   &
                                   cclosure_relaxcbmf, &
                                   cclosure_relaxwfn,  &
                                   cclosure_implicit, cclosure

!---------------------------------------------------------------------
  implicit none
  private
!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: uw_conv.F90,v 15.0.2.1.2.1.2.1.2.1.2.1 2008/02/02 10:25:31 rsh Exp $'
  character(len=128) :: tagname = '$Name: omsk_2008_03 $'

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: uw_conv, uw_conv_init, uw_conv_end

  real, parameter :: mv = -999.
  logical         :: module_is_initialized = .false.

  character(len=7) :: mod_name = 'uw_conv'

  !namelist parameters for UW convection scheme
  integer :: iclosure = 0      ! 0: Bretherton UWShCu orginal / -CIN/TKE based
                               ! 1: Emanuel-Rayment: quasiequilibrium PBL
  integer :: mixing_assumption = 0
  real    :: rkm_sh   = 16.0   ! fractional lateral mixing rate for shallow
  real    :: cldhgt_max   = 4.e3
  real    :: cbmf_dp_frac = 0.0
  real    :: landfact_m   = 0.0
  integer :: idpchoice = 0
  logical :: do_deep = .false.
  logical :: do_relaxcape = .false.
  logical :: do_relaxwfn  = .false.
  logical :: do_coldT = .true.
  logical :: do_lands = .false.
  logical :: do_uwcmt = .false.   
  logical :: do_fast  = .false.
  logical :: do_ice   = .true.
  logical :: do_ppen  = .true.
  logical :: do_pevap = .false.
  logical :: do_micro = .false.   
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
  logical :: do_rescale   = .false.
  real    :: pblht0 = 500.
  real    :: tke0 = 1.
  real    :: lofactor0 = 1.
  integer :: lochoice  = 0
  real    :: wrel_min = 0.
  real    :: om_to_oc = 1.67
  real    :: sea_salt_scale = 0.1

  NAMELIST / uw_conv_nml / iclosure, rkm_sh, cldhgt_max, cbmf_dp_frac, &
       do_deep, idpchoice, do_relaxcape, do_relaxwfn, do_coldT, do_lands, do_uwcmt,       &
       do_fast, do_ice, do_ppen, do_pevap, do_micro, mixing_assumption, do_forcedlifting, &
       atopevap, apply_tendency, prevent_unreasonable, aerol, gama, tkemin,    &
       wmin_ratio, use_online_aerosol, landfact_m, pblht0, tke0, lofactor0, lochoice, &
       do_auto_aero, do_rescale, wrel_min, om_to_oc, sea_salt_scale
       

  !namelist parameters for UW convective plume
  real    :: rle      = 0.10   ! for critical stopping distance for entrainment
  real    :: rpen     = 5.0    ! for entrainment efficiency
  real    :: rmaxfrac = 0.05   ! maximum allowable updraft fraction
  real    :: wmin     = 0.0    ! maximum allowable updraft fraction
  real    :: rbuoy    = 1.0    ! for nonhydrostatic pressure effects on updraft
  real    :: rdrag    = 1.0 
  real    :: frac_drs = 0.0    ! 
  real    :: bigc     = 0.7    ! for momentum transfer
  real    :: auto_th0 = 0.5e-3 ! threshold for precipitation
  real    :: auto_rate= 1.e-3
  real    :: tcrit    = -45.0  ! critical temperature 
  real    :: deltaqc0 = 0.5e-3
  logical :: do_pdfpcp= .false.
  logical :: do_pmadjt= .false.
  logical :: do_emmax = .false.
  logical :: do_pnqv = .false.
  real    :: rad_crit = 14.0   ! critical droplet radius
  real    :: emfrac_max = 1.0

  NAMELIST / uw_plume_nml / rle, rpen, rmaxfrac, wmin, rbuoy, rdrag, frac_drs, bigc, &
       auto_th0, auto_rate, tcrit, deltaqc0, do_pdfpcp, do_pmadjt, do_emmax, do_pnqv, rad_crit, emfrac_max
 
  !namelist parameters for UW convective closure
  integer :: igauss   = 1      ! options for cloudbase massflux closure
                               ! 1: cin/gaussian closure, using TKE to compute CIN.
                               ! 2: cin/gaussian closure, using W* to compute CIN.
                               ! 0: cin/tke mapse-style closure; 
  real    :: rkfre    = 0.05   ! vertical velocity variance as fraction of tke
  real    :: tau_sh   = 7200.  ! 
  real    :: wcrit_min= 0.

  NAMELIST / uw_closure_nml / igauss, rkfre, tau_sh, wcrit_min

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
       id_ocode_uwc, id_plnb_uwc, id_wrel_uwc, id_ufrc_uwc, id_qtmp_uwc,&
       id_tdt_pevap_uwc, id_qdt_pevap_uwc, id_xhlsrc_uwc, id_xqtsrc_uwc,&
       id_qldet_uwc, id_qidet_uwc, id_qadet_uwc, id_qtdt_uwc

  integer, allocatable :: id_tracerdt_uwc(:), id_tracerdt_uwc_col(:), &
                          id_tracerdtwet_uwc(:), id_tracerdtwet_uwc_col(:)

  logical   ::  do_unified_convective_closure

  type(sounding),   save  :: sd, sd1
  type(adicloud),   save  :: ac, ac1
  type(cclosure),   save  :: cc, cc1
  type(cplume)  ,   save  :: cp, cp1
  type(ctend)   ,   save  :: ct, ct1
  type(cpnlist),    save  :: cpn
  type(uw_params),  save  :: Uw_p

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
    
    integer   :: ntracers, n, nn, ierr
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
          ierr = check_nml_error(io,'uw_closure_nml')
       end do
10     call close_file ( unit )

       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_conv_nml, iostat = io, end = 20 )
          ierr = check_nml_error(io,'uw_conv_nml')
       end do
20     call close_file ( unit )
       
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_plume_nml, iostat = io, end = 30 )
          ierr = check_nml_error(io,'uw_plume_nml')
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
    cpn % deltaqc0  = deltaqc0
    cpn % do_pdfpcp = do_pdfpcp
    cpn % do_pmadjt = do_pmadjt
    cpn % do_emmax  = do_emmax
    cpn % do_pnqv   = do_pnqv 
    cpn % emfrac_max= emfrac_max
    cpn % auto_rate = auto_rate
    cpn % tcrit     = tcrit  
    cpn % cldhgt_max= cldhgt_max
    cpn % do_ice    = do_ice
    cpn % do_ppen   = do_ppen
    cpn % do_pevap  = do_pevap
    cpn % mixing_assumption= mixing_assumption
    cpn % do_micro  = do_micro
    cpn % do_forcedlifting= do_forcedlifting
    cpn % atopevap  = atopevap
    cpn % wtwmin_ratio = wmin_ratio*wmin_ratio
    cpn % do_auto_aero = do_auto_aero
    cpn % rad_crit = rad_crit
    cpn % wrel_min = wrel_min

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
    cc  % tau_sh    = tau_sh
 
    id_xhlsrc_uwc = register_diag_field (mod_name,'xhlsrc_uwc', axes(1:2), Time, &
         'xhlsrc', 'J/kg' )
    id_xqtsrc_uwc = register_diag_field (mod_name,'xqtsrc_uwc', axes(1:2), Time, &
         'xqtsrc', 'kg/kg' )
 
    id_tdt_pevap_uwc = register_diag_field ( mod_name, 'tdt_pevap_uwc', axes(1:3), Time, &
         'Temperature tendency due to pevap from uw_conv', 'K/day', missing_value=mv )
    id_qdt_pevap_uwc = register_diag_field ( mod_name, 'qdt_pevap_uwc', axes(1:3), Time, &
         'Spec. humidity tendency due to pevap from uw_conv', 'kg/kg/day', missing_value=mv)

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
       id_qldet_uwc = register_diag_field (mod_name,'qldet_uwc',axes(1:3),Time, &
            'ql detrainment', 'kg/kg/day', missing_value=mv)
       id_qidet_uwc = register_diag_field (mod_name,'qidet_uwc',axes(1:3),Time, &
            'qi detrainment', 'kg/kg/day', missing_value=mv)
       id_qadet_uwc = register_diag_field (mod_name,'qadet_uwc',axes(1:3),Time, &
            'qa detrainment', '1/day', missing_value=mv)
       id_qtdt_uwc= register_diag_field (mod_name,'qtdt_uwc',axes(1:3),Time, &
            'Total water tendency from uw_conv', 'kg/kg/day', missing_value=mv)
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
!!5 miz drops:
!       cbmf_clo,  &
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
!!5 miz drops
!    real, intent(in), dimension(:,:)  :: cbmf_clo  ! cloud-base mass flux
    real, intent(inout), dimension(:,:)  :: cbmfo  ! cloud-base mass flux
    real, intent(in),  dimension(:,:,:,:)  :: tracers         ! env. tracers
    real, intent(out), dimension(:,:,:,:)  :: trtend          ! calculated tracer tendencies

    integer i, j, k, kl, klm, km1, nk, ks, m, naer, na, n

    real rhos0j, thj, qvj, qlj, qij, qse, thvj
    real thvuinv, hlsrc, thcsrc, qctsrc, plcltmp, tmp, lofactor
    real zsrc, psrc, cbmf_shallow, cbmf_old, cbmf_deep, rkm_sh1, omeg_avg, dpsum

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
         dwfno,      &     ! dwfn(m2/s2)
!!5 miz drops as input arg; define as local
         cbmf_clo,   &
         ocode,      &
         xhlsrc,     &
         xqtsrc

    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: wuo,fero,fdro,fdrso, tten_pevap, qvten_pevap
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: qldet, qidet, qadet

    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: tempdiag1
    real, dimension(size(tb,1),size(tb,2))            :: qtin, dqt, scale_uw

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
    tten_pevap=0.; qvten_pevap=0.;

    cino=0.; capeo=0.; tkeo=0.; wrelo=0.; ufrco=0.; zinvo=0.; wuo=0.; 
    fero=0.; fdro=0.; fdrso=0.; cmf=0.; denth=0.;  dqtmp=0.; ocode=0;
    dcapeo=0.; dcino=0.; dwfno=0.; xhlsrc=0.; xqtsrc=0.;
    trtend=0.; qldet=0.; qidet=0.; qadet=0.;
    trwet = 0.

    naer = size(asol%aerosol,4)

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
              am2(k)= sea_salt_scale*asol%aerosol(i,j,k,5)*tmp
              am3(k)= om_to_oc*asol%aerosol(i,j,k,3)*tmp
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
            !wstar   = (ustar(i,j)*bstar(i,j)*pblht(i,j))**(1./3.)
             cpn % auto_th0 = auto_th0 * (1. + landfact_m * sd%land)
             call qt_parcel_k (sd%qs(1), qstar(i,j), pblht(i,j), tkeo(i,j), sd%land, gama, &
                  pblht0, tke0, lofactor0, lochoice, qctsrc, lofactor)
             rkm_sh1 = rkm_sh   * lofactor
          endif
          
          call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, do_fast, do_ice, ac)
          ac % usrc = sd%u(sd%ktoppbl)
          ac % vsrc = sd%v(sd%ktoppbl)
          plcl (i,j) = ac%plcl
          plfc (i,j) = ac%plfc
          plnb (i,j) = ac%plnb
          cino (i,j) = ac%cin
          capeo(i,j) = ac%cape
          xhlsrc(i,j)= ac%hlsrc;
          xqtsrc(i,j)= ac%qctsrc;
          if (do_fast) then
            if (ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.) then
              ocode(i,j)=1; goto 100 !cycle;
            end if
            if (ac%plfc.lt.500.) then
              ocode(i,j)=2; goto 100 !cycle;
            end if
          end if



!========Cumulus closure to determine cloud base mass flux===================


        cbmf_old=cbmfo(i,j) 


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
          wrelo(i,j) = cc%wrel
          ufrco(i,j) = cc%ufrc

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

!       cbmftmp=cc%cbmf * cbmf_sh_frac
!       call cumulus_plume_k (cpn, sd, ac, cp, rkm_sh1, cbmftmp, &
!                             cc%wrel, cc%scaleh, Uw_p)
          cbmf_deep    = min(cbmf_dp_frac * cc%cbmf * cc%cbmf / cc%wrel,cc%cbmf*0.9)
          cbmf_shallow = cc%cbmf - cbmf_deep
          call cumulus_plume_k(cpn, sd, ac, cp, rkm_sh1, cbmf_shallow, cc%wrel, cc%scaleh, Uw_p)
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
             qldet (i,j,nk) = ct%qldet(k)
             qidet (i,j,nk) = ct%qidet(k)
             qadet (i,j,nk) = ct%qadet(k)
             qvten (i,j,nk) = ct%qvten(k)
             pflx  (i,j,nk) = ct%pflx (k)
             tten  (i,j,nk) = ct%tten (k)
             rhos0j = sd%ps(k)/(rdgas*0.5*(cp%thvbot(k+1)+cp%thvtop(k))*sd%exners(k))
             hlflx(i,j,nk) = ct%hlflx(k)
             qtflx (i,j,nk) = rhos0j*HLv*ct%qctflx(k)
             tten_pevap (i,j,nk) = ct%tevap (k)
             qvten_pevap(i,j,nk) = ct%qevap (k)
             
             cldqa (i,j,nk) = cp%ufrc(k)
             cldql (i,j,nk) = cp%qlu(k)
             cldqi (i,j,nk) = cp%qiu(k)
             cldqn (i,j,nk) = cp%qnu(k)
             cmf   (i,j,nk) = cp%umf(k) !+ cp%emf(k)
             wuo   (i,j,nk) = cp%wu (k)
             fero  (i,j,nk) = cp%fer(k)
             fdro  (i,j,nk) = cp%fdr(k)
             fdrso (i,j,nk) = cp%fdrsat(k)*cp%umf(k)
             
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
          end if
!========Option for deep convection=======================================

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
          qvten(:,:,:) = qvten(:,:,:)+(q(:,:,:,nql)/delt+qlten(:,:,:))
          qlten(:,:,:) = qlten(:,:,:)-(q(:,:,:,nql)/delt+qlten(:,:,:))
      end where
 
      where ((q(:,:,:,nqi) + qiten*delt) .lt. 0.)
          tten (:,:,:) = tten(:,:,:) -(q(:,:,:,nqi)/delt+qiten(:,:,:))*HLs/Cp_Air
          qvten(:,:,:) = qvten(:,:,:)+(q(:,:,:,nqi)/delt+qiten(:,:,:))
    
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

   if (do_rescale) then
      !rescaling to prevent negative specific humidity for each grid point
      tempdiag1 = 1.0
      do k=1,kmax
         qtin =  q(:,:,k,nqv)
         dqt  =  qvten(:,:,k) * delt
         where ( dqt.lt.0 .and. qtin+dqt.lt.1.e-10 )
            tempdiag1(:,:,k) = max( 0.0, -(qtin-1.e-10)/dqt )
         end where
      end do

      !scaling factor for each column is the minimum value within that column
      scale_uw = minval( tempdiag1, dim=3 )

      !scale tendencies
      do k=1,kmax
         uten (:,:,k)  = scale_uw(:,:) * uten (:,:,k)
         vten (:,:,k)  = scale_uw(:,:) * vten (:,:,k)
         tten (:,:,k)  = scale_uw(:,:) * tten (:,:,k)
         qvten(:,:,k)  = scale_uw(:,:) * qvten(:,:,k)
         qlten(:,:,k)  = scale_uw(:,:) * qlten(:,:,k)
         qiten(:,:,k)  = scale_uw(:,:) * qiten(:,:,k)
         qaten(:,:,k)  = scale_uw(:,:) * qaten(:,:,k)
      end do

      if (do_qn) then
         do k=1,kmax
            qnten(:,:,k) = scale_uw(:,:) * qnten(:,:,k)
         end do
      end if

      rain(:,:) = scale_uw(:,:) * rain(:,:)
      snow(:,:) = scale_uw(:,:) * snow(:,:)
 end if
 
    endif


    !diagnostic output
    if ( id_xhlsrc_uwc   > 0 ) &
         used = send_data( id_xhlsrc_uwc,       xhlsrc,             Time, is, js)
    if ( id_xqtsrc_uwc   > 0 ) &
         used = send_data( id_xqtsrc_uwc,       xqtsrc,             Time, is, js)
    if ( id_tdt_pevap_uwc   > 0 ) &
         used = send_data( id_tdt_pevap_uwc,    tten_pevap*86400. , Time, is, js, 1)
    if ( id_qdt_pevap_uwc   > 0 ) &
         used = send_data( id_qdt_pevap_uwc,    qvten_pevap*86400., Time, is, js, 1)

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

    if ( id_wrel_uwc > 0 ) &
         used = send_data( id_wrel_uwc, (wrelo),            Time, is, js )
    if ( id_ufrc_uwc > 0 ) &
         used = send_data( id_ufrc_uwc, (ufrco),            Time, is, js )
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
       if ( id_qldet_uwc > 0 ) &
            used = send_data( id_qldet_uwc,  qldet*86400.,  Time, is, js, 1)
       if ( id_qidet_uwc > 0 ) &
            used = send_data( id_qidet_uwc,  qidet*86400.,  Time, is, js, 1)
       if ( id_qadet_uwc > 0 ) &
            used = send_data( id_qadet_uwc,  qadet*86400.,  Time, is, js, 1)
       if ( id_qtdt_uwc  > 0 ) &
            used = send_data( id_qtdt_uwc,(qvten+qlten+qiten)*86400.,Time, is, js, 1)
       
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

end MODULE UW_CONV_MOD
