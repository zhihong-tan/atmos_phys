                       module donner_deep_mod

use time_manager_mod,       only: time_type, set_time, &
                                  set_date, get_time,   &
                                  get_calendar_type, &
                                  operator(-), &
                                  operator(>=), operator (<)
use diag_manager_mod,       only: register_diag_field, send_data
use field_manager_mod,      only: MODEL_ATMOS, field_manager_init, &
                                  fm_query_method, get_field_info, &
                                  parse
use tracer_manager_mod,     only: get_tracer_names,get_number_tracers, &
                                  get_tracer_indices, &
!++lwh
                                  query_method
use atmos_tracer_utilities_mod, only : get_wetdep_param
use  sat_vapor_pres_mod,only : sat_vapor_pres_init
!--lwh
use fms_mod,                only: mpp_pe, mpp_root_pe,  &
                                  file_exist,  check_nml_error,  &
                                  error_mesg, FATAL, WARNING, NOTE,  &
                                  close_file, open_namelist_file,    &
                                  stdlog, write_version_number,  &
                                  field_size, &
                                  read_data, write_data, lowercase,    &
                                  open_restart_file
use mpp_io_mod,             only: mpp_open, mpp_close, fieldtype,  &
                                  mpp_read_meta, mpp_get_info, &
                                  mpp_get_fields, mpp_read, &
                                  MPP_NETCDF, MPP_SINGLE,   &
                                  MPP_SEQUENTIAL, MPP_RDONLY, MPP_NATIVE, &
                                  mpp_get_field_name
use constants_mod,          only: DENS_H2O, RDGAS, GRAV, CP_AIR,  &
                                  pie=>PI, KAPPA, RVGAS, &
                                  SECONDS_PER_DAY, HLV, HLF, HLS, KELVIN
use column_diagnostics_mod, only: initialize_diagnostic_columns, &
                                  column_diagnostics_header, &
                                  close_column_diagnostics_units
use donner_types_mod,       only: donner_initialized_type, &
                                  donner_save_type, donner_rad_type, &
                                  donner_nml_type, donner_param_type, &
                                  donner_budgets_type, &
                                  donner_column_diag_type, &
                                  MAXMAG, MAXVAL, MINMAG, MINVAL, &
                                  DET_MASS_FLUX, MASS_FLUX,  &
                                  CELL_UPWARD_MASS_FLUX, TEMP_FORCING, &
                                  MOIST_FORCING, PRECIP,  FREEZING, &
                                  RADON_TEND, &
                                  donner_conv_type, donner_cape_type
use  conv_utilities_k_mod,  only: sd_init_k, sd_end_k, ac_init_k,  &
                                  ac_end_k, uw_params_init_k, &
                                  exn_init_k, exn_end_k, findt_init_k, &
                                  findt_end_k, &
                                  adicloud, sounding, uw_params
use  conv_plumes_k_mod,     only: cp_init_k, cp_end_k, ct_init_k,  &
                                  ct_end_k, cplume, ctend

implicit none
private

!--------------------------------------------------------------------
!        donner_deep_mod diagnoses the location and computes the 
!        effects of deep convection on the model atmosphere
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------


character(len=128)  :: version =  '$Id: donner_deep.F90,v 15.0.2.1.2.1.2.1.2.1.2.1 2007/11/13 11:30:02 rsh Exp $'
character(len=128)  :: tagname =  '$Name: omsk_2007_12 $'


!--------------------------------------------------------------------
!---interfaces------

public   &
        donner_deep_init, donner_deep, donner_deep_end

private   &
!  module subroutines called by donner_deep_init:
        register_fields, read_restart, read_restart_nc,  &
        process_coldstart,&
!  module subroutines called by donner_deep:
        donner_deep_netcdf, donner_column_control,     &
!  module subroutines called from donner_deep_end:
        write_restart, write_restart_nc, deallocate_variables


!---------------------------------------------------------------------
!---namelist----

!----------------------------------------------------------------------
!  the following nml variables are stored in a donner_nml_type derived
!  type variable so they may be conveniently passed to kernel subroutines
!  as needed:
!----------------------------------------------------------------------

integer, parameter  :: MAX_ENSEMBLE_MEMBERS = 7
                             ! maximum number of cumulus ensemble 
                             ! members
integer             :: model_levels_in_sfcbl = 2 
                             ! number of levels at which the temperature 
                             ! and vapor profiles are not allowed to 
                             ! change from lag value when calculating the
                             ! time tendency of cape
integer             :: parcel_launch_level = 2 
                             ! large-scale model level from which a par-
                             ! cel is launched to determine the lifting 
                             ! condensation level (level 1 nearest the 
                             ! surface)
logical             :: allow_mesoscale_circulation = .true.
                             ! a mesoscale circulation will be included 
                             ! in those columns which satisfy the 
                             ! required conditions ?
 real               ::  CDEEP_CV = 100.   
                        ! maximum value of convective inhibition (J/kg) 
                        ! that allows convection. Value of 10 suggested 
                        ! by Table 2 in Thompson et al. (1979, JAS).
logical             :: do_freezing_for_cape = .false.
                        ! include freezing in cape parcel calculation
real                :: tfre_for_cape = 258.0
                        ! temperature at which freezing begins for   
                        ! parcel in cape parcel calculation [ deg K ]
real                :: dfre_for_cape =  10.
                        ! all liquid freezes between tfre_for_cape and 
                        ! tfre_for_cape + dfre_for_cape  in cape parcel
                        ! calculation [ deg K ]
real                :: rmuz_for_cape = 0.0       
                        ! parcel entrainment factor used in cape parcel
                        ! calculation
logical             :: do_freezing_for_closure = .false.
                        ! include freezing in closure calculation
real                :: tfre_for_closure = 258.0
                        ! temperature at which freezing begins for
                        ! parcel in closure calculation [ deg K ]
real                :: dfre_for_closure =  10.
                        ! all liquid freezes tfre_for_closure and 
                        ! tfre_for_closure + dfre_for_closure in  
                        ! closure calculation [ deg K ]
real                :: rmuz_for_closure = 0.0     
                        ! parcel entrainment factor used in cape closure
                        ! calculation
logical             :: do_donner_cape    = .true.
logical             :: do_donner_plume   = .true.
logical             :: do_donner_closure = .true.
logical             :: do_dcape          = .true.
logical             :: do_lands          = .false.
real                :: tau               = 28800.
real                :: cape0             = 0.
real                :: rhavg0            = 1.
real                :: plev0             = 50000.
logical             :: do_rh_trig        = .false.
logical             :: do_capetau_land   = .false.
real                :: pblht0            = 500.
real                :: tke0              = 1.
real                :: lofactor0         = 1.
real                :: deephgt0          = 4000.
integer             :: lochoice          = 0
integer             :: deep_closure      = 0
real                :: gama              = 0.0001
logical             :: do_ice            = .false.
real                :: atopevap          = 0.
logical             :: do_donner_lscloud = .true.
logical             :: use_llift_criteria =.true.
logical             :: use_pdeep_cv       =.true.
real                :: auto_rate = 1.0e-3
real                :: auto_th   = 0.5e-3
real                :: frac      = 1.
real                :: ttend_max = 1000.
real                :: EVAP_IN_DOWNDRAFTS  = 0.25
                        ! fraction of condensate available to the meso-
                        ! scale which is evaporated in convective down-
                        ! drafts (from Leary and Louze, 1980)
                        ! [ dimensionless ]
real                :: EVAP_IN_ENVIRON     = 0.13
                        ! fraction of condensate available to the meso-
                        ! scale which is evaporated in the cell environ-
                        ! ment (from Leary and Louze, 1980)
                        ! [ dimensionless ]
real                :: ENTRAINED_INTO_MESO = 0.62
                        ! fraction of condensate available to the meso-
                        ! scale which is entrained into the mesoscale 
                        ! circulation (from Leary and Louze, 1980)
                        ! [ dimensionless ]
real                :: ANVIL_PRECIP_EFFICIENCY = 0.5
                        ! fraction of total condensate in anvil (trans-
                        ! fer from cell plus in situ condensation) that 
                        ! precipitates out (from Leary and Louze, 1980)
                        ! [ dimensionless ]
real                :: MESO_DOWN_EVAP_FRACTION = 0.4
                        ! fraction of total anvil condensate assumed 
                        ! evaporated in the mesoscale downdraft
                        ! [ fraction ]
real                :: MESO_UP_EVAP_FRACTION   = 0.1
                        ! fraction of total anvil condensate assumed 
                        ! evaporated in outflow from the mesoscale 
                        ! updraft [ fraction ]    
real,dimension(MAX_ENSEMBLE_MEMBERS)   :: arat
data  arat / 1.0, 0.26, 0.35, 0.32, 0.3, 0.54, 0.66 /
                        ! ratio at cloud base of the fractional area of 
                        ! ensemble member i relative to ensemble member 
                        ! 1. (taken from GATE data).
real,dimension(MAX_ENSEMBLE_MEMBERS)  :: erat
data  erat / 1.0, 1.30, 1.80, 2.50, 3.3, 4.50, 10.0 /
                        ! ratio of entrainment constant between ensemble
                        ! member 1 and ensemble member i for gate-based
                        ! ensemble

integer             :: donner_deep_freq = 1800 
                             ! frequency of calling donner_deep [ sec ]; 
                             ! must be <= 86400 
character(len=16)   :: cell_liquid_size_type = 'bower' 
                             ! choose either 'input' or 'bower' 
real                :: cell_liquid_eff_diam_input = -1.0 
                             ! input cell droplet effective diameter 
                             ! [ microns ];
                             ! needed when cell_liquid_size_type == 
                             ! 'input'
character(len=16)   :: cell_ice_size_type = 'default' 
                             ! choose either 'input' or 'default'
real                :: cell_ice_geneff_diam_input = -1.0 
                             ! input cell ice generalized effective diam-
                             ! eter [ microns ]; needed when 
                             ! cell_ice_size_type == 'input'
real                :: meso_liquid_eff_diam_input = -1.0 
                             ! input mesoscale droplet effective diameter
                             ! [ microns ]; currently no liquid allowed 
                             ! in mesoscale cloud
logical             :: do_average = .false.   
                             ! time-average donner cloud properties for 
                             ! use by radiation package?
character(len=32)   :: entrainment_constant_source = 'gate'
                             ! source of cloud entrainment constants for
                             !  cumulus ensemble; either 'gate' or 'kep'
logical             :: use_memphis_size_limits = .false.
                             ! this should be set to .true. only to
                             ! eliminate the answer changes resulting
                             ! from the change to drop size limits
                             ! introduced in  modset
                             ! memphis_cloud_diagnostics_rsh. in any 
                             ! new or ongoing production, this variable
                             ! should be .false.
real                :: wmin_ratio = 0. 
                             ! ratio of current updraft velocity to 
                             ! maximum updraft velocity of plume at 
                             ! which time plume is assumed to stop
                             ! rising; used in donner lite 

logical             :: do_budget_analysis = .false.
                             ! save arrays containing terms involved in
                             ! enthalpy and water budgets for netcdf
                             ! output ?

logical             :: force_internal_enthalpy_conservation = .false.
                             ! modify the temperature tendency so that
                             ! enthalpy is conserved by the cell and
                             ! mesoscale motions, rather than forcing
                             ! entropy conservation

!----------------------------------------------------------------------
!   The following nml variables are not needed in any kernel subroutines
!   and so are not included in the donner_nml_type variable:
!----------------------------------------------------------------------

logical             :: do_netcdf_restart = .true.
                             ! restart file written in netcdf format 
                             ! (option is native mode) 
                             ! NOTE: current code will produce ONLY
                             ! a netcdf restart; if native mode restart 
                             ! is desired, user must update the source 
                             ! code appropriately.
logical             :: write_reduced_restart_file = .false.
                             ! by setting this variable to .true., the 
                             ! user is asserting that the donner deep 
                             ! calculation will be made on the first step
                             ! of any job reading the restart file, so 
                             ! that those variables needed on steps when
                             ! donner_deep is not calculated may be 
                             ! omitted from the file, thus saving archive
                             ! space. a check is provided in the code 
                             ! which will force the writing of a full 
                             ! file (even with this variable set to 
                             !.true.), if it is determined that the cur-
                             ! rent job does not end just prior to a 
                             ! donner step. 
                             ! user should set this variable to .false. 
                             ! if it is known that donner_deep_freq is to
                             ! be changed at the next restart to avoid a
                             ! fatal error in that job.
integer, parameter  :: MAX_PTS = 20
                             ! max nunber of diagnostic columns
real                :: diagnostics_pressure_cutoff =  50.E+02   
                             ! column data will be output on model levels
                             ! with pressures greater than this value 
                             ! [ Pa ]
integer,                        &
    dimension(6)    :: diagnostics_start_time= (/ 0,0,0,0,0,0 /)
                             ! integer specification of time at which 
                             ! column diagnostics is to be activated
                             ! [year, month, day, hour, min, sec ]
integer             :: num_diag_pts_ij = 0   
                             ! number of diagnostic columns specified by 
                             ! global(i,j) coordinates
integer             :: num_diag_pts_latlon = 0 
                             ! number of diagnostic columns specified by
                             ! lat-lon coordinates 
integer,                       &
 dimension(MAX_PTS) :: i_coords_gl = -100
                             ! global i coordinates for ij diagnostic 
                             ! columns
integer,                        &
 dimension(MAX_PTS) :: j_coords_gl = -100
                             ! global j coordinates for ij diagnostic 
                             ! columns
real,                           &
 dimension(MAX_PTS) :: lat_coords_gl = -999.
                             ! latitudes for lat-lon  diagnostic columns 
                             ! [ degrees, -90. -> 90. ]
real,                            &
 dimension(MAX_PTS) :: lon_coords_gl = -999.
                             ! longitudes for lat-lon  diagnostic columns
                             ! [ degrees, 0. -> 360. ]

namelist / donner_deep_nml /      &

! contained in donner_rad_type variable:
                            model_levels_in_sfcbl, &
                            parcel_launch_level, &
                            allow_mesoscale_circulation, &
                            cdeep_cv,         &
                            do_freezing_for_cape, tfre_for_cape, &
                            dfre_for_cape, rmuz_for_cape, &
                            do_freezing_for_closure, tfre_for_closure, &
                            dfre_for_closure, rmuz_for_closure, &
                            do_donner_cape,   &!miz
                            do_donner_plume,  &!miz
                            do_donner_closure,&!miz
                            do_dcape,         &!miz
                            do_lands,         &!miz
                            tau,              &!miz
                            cape0,            &!miz
                            rhavg0,           &!miz
                            plev0,            &!miz
                            do_rh_trig,       &!miz
                            do_capetau_land,  &!miz
                            pblht0,           &!miz
                            tke0,             &!miz
                            lofactor0,        &!miz
                            deephgt0,         &!miz
                            lochoice,         &!miz
                            deep_closure,     &!miz
                            gama,             &!miz
                            do_ice,           &!miz
                            atopevap,         &!miz
                            do_donner_lscloud,&!miz
                            use_llift_criteria,&!miz
                            use_pdeep_cv,     &!miz
                            auto_rate,        &!miz
                            auto_th,          &!miz
                            frac,             &!miz
                            ttend_max,        &!miz
                            EVAP_IN_DOWNDRAFTS,      &!miz
                            EVAP_IN_ENVIRON,         &!miz
                            ENTRAINED_INTO_MESO,     &!miz
                            ANVIL_PRECIP_EFFICIENCY, &!miz
                            MESO_DOWN_EVAP_FRACTION, &!miz
                            MESO_UP_EVAP_FRACTION,   &!miz
                            arat,                    &!miz
                            erat,                    &
                            donner_deep_freq, &
                            cell_liquid_size_type,   &
                            cell_liquid_eff_diam_input, &
                            cell_ice_size_type, &
                            cell_ice_geneff_diam_input, &
                            meso_liquid_eff_diam_input, &
                            do_average,  &
                            entrainment_constant_source, &
                            use_memphis_size_limits, &
                            wmin_ratio, &
                            do_budget_analysis, &
                            force_internal_enthalpy_conservation, &

! not contained in donner_nml_type variable:
                            do_netcdf_restart, &
                            write_reduced_restart_file, &
                            diagnostics_pressure_cutoff, &
                            diagnostics_start_time, &
                            num_diag_pts_ij, num_diag_pts_latlon, &
                            i_coords_gl, j_coords_gl, &
                            lat_coords_gl, lon_coords_gl


!--------------------------------------------------------------------
!--- public data ----------




!--------------------------------------------------------------------
!----private data-----------



!---------------------------------------------------------------------
!  parameters stored in the donner_param derived type variable to facili-
!  tate passage to kernel subroutines:
!

real,                          &
  parameter                     &
             ::  CP_VAPOR= 4.0*RVGAS  
                       ! specific heat of water vapor at constant pres-
                       ! sure [ J/(kg K) ]
real,                             &
  parameter                    &
             ::  D622 = RDGAS/RVGAS 
                       ! ratio of molecular weights of water vapor and 
                       ! dry air
real,                        &
  parameter                 &
             ::  D608 = (RVGAS/RDGAS - 1.0)
                        ! factor in virtual temperature definition
integer,                   &
  parameter                  &
             ::  KPAR=7         
                        ! number of members in cumulus ensemble
integer,                   &
  parameter                  &
             ::  NLEV_HIRES=100       
                        ! number of levels in cloud model
real,                       &
  parameter                 &
             ::  PDEEP_CV = 500.e02 
                        ! minimum pressure difference between level of 
                        ! free convection and level of zero buoyancy 
                        ! needed for deep convection to occur [ Pa ].
real,                       &
  parameter                 &
             ::  MAX_ENTRAINMENT_CONSTANT_GATE = 0.0915
                        ! entrainment constant based on gate data for 
                        ! most entraining ensemble member
real,                       &
  parameter                 &
             ::  MAX_ENTRAINMENT_CONSTANT_KEP  = 0.0915
                        ! entrainment constant based on kep data for most
                        ! entraining ensemble member
real,                       &
  parameter,                 &
  dimension(KPAR)              &
             ::  ENSEMBLE_ENTRAIN_FACTORS_KEP  = (/ 1.0, 1.22, 1.56,  &
                                                    2.05, 2.6, 3.21,   &
                                                    7.84 /)
                        ! ratio of entrainment constant between ensemble
                        ! member 1 and ensemble member i for kep-based
                        ! ensemble
real,                       &
  parameter                 &
             ::  CLD_BASE_VERT_VEL = 0.5                          
                        ! vertical velocity assumed present at cloud base
                        ! [ m / sec ]
real,                       &
  parameter                 &
             ::  PSTOP = 40.0e02   
                        ! lowest possible pressure to which a cloud may 
                        ! extend in the cloud model [ Pa ]
real,                       &
  parameter                 &
             ::  PARCEL_DP = -1.0e02
                        ! pressure increment used for parcel calculations
                        ! [ Pa ]
real,                       &
  parameter                 &
             ::  UPPER_LIMIT_FOR_LFC = 500.e02 
                        ! lowest pressure allowed for level of free conv-
                        ! ection [ Pa ]
real,                       &
  parameter                 &
             ::  DP_OF_CLOUD_MODEL = -10.e02
                        ! pressure thickness (Pa) of the layers in the
                        ! donner parameterization's cloud model.
real,                       &
  parameter                 &
             ::  CLOUD_BASE_RADIUS = 1000.
                        ! radius assumed for cloud ensemble member #1 at
                        ! cloud base [ m ]
real,                       &
  parameter                 &
             ::  WDET = .1   
                        ! vertical velocity at which detrainment from the
                        ! clouds begins [ m/s ]
real,                       &
  parameter                 &
             ::  RBOUND = 0.01    
                        ! value of cumulus radius at which cloud effect-
                        ! ively disappears and cloud model calculation 
                        ! stops [ m ]
real,                       &
  parameter                 &
             ::  WBOUND = 0.01  
                        ! value of cumulus vertical velocity at which 
                        ! cloud model calculation stops [ m / sec ]
real,                       &
  parameter                 &
             ::  FREEZE_FRACTION = 0.52
                        ! fraction of liquid in cloud updraft which may 
                        ! be frozen. (Leary and Houze (JAS,1980)) 
                        ! [ dimensionless ]
real,                       &
  parameter                 &
             ::  VIRT_MASS_CO = 0.5
                        ! virtual mass coefficient [ dimensionless ]
real,                       &
  parameter                 &
             ::  PDEEP_MC = 200.e02 
                        ! pressure thickness [ Pa ] required for meso-
                        ! scale circulation. It refers to the least
                        ! penetrative ensemble member. For this check 
                        ! to function properly, the entrainment coeffic-
                        ! ient in cloud_model for kou=1 must be the 
                        ! largest entrainment coefficient.
real,                       &
  parameter                 &
             ::  TR_INSERT_TIME = 0.0
                        ! fractional point (based on mass increase) 
                        ! during a timestep at which an entraining parcel
                        ! takes on internally-generated tracer 
                        ! [ dimensionless, value between 0.0 and 1.0 ]
real,                       &
  parameter                 &
             ::  AUTOCONV_RATE = 1.0e-03
                        ! rate of autoconversion of cloud to rainwater 
                        ! [ sec**(-1) ]
real,                       &
  parameter                 &
             ::  AUTOCONV_THRESHOLD =  0.5    
                        ! threshold of cloud water at which autoconver-
                        ! sion of cloud to rainwater begins  [ g / m**3 ]
real,                       &
  parameter                 &
             ::  TFRE = 258.  
                        ! temperature at which cloud liquid begins to 
                        ! freeze [ deg K ]
real,                       &
  parameter                 &
             ::  DFRE = 10.   
                        ! range of temperature between the onset and 
                        ! completion of freezing  [ deg K ]
real,                       &
  parameter                 &
             ::  UPPER_LIMIT_FOR_LCL = 500.0E02
                        ! lowest pressure allowable for lifting condens-
                        ! ation level; deep convection will not be pres-
                        ! ent if lcl not reached before this pressure 
                        ! [ Pa ]
integer,                   &
  parameter         &
             ::  ISTART = 1    
                        ! index of level in cape grid from which the 
                        ! parcel originates for the cape calculations
real,                       &
  parameter                 &
             ::  TMIN = 154.       
                        ! cape calculations are terminated when parcel 
                        ! temperature goes below TMIN [ deg K ]
real,                       &
  parameter                 &
             ::  MESO_LIFETIME = 64800.
                        ! assumed lifetime of mesoscale circulation 
                        ! (from Leary and Louze, 1980) [ sec ]
real,                       &
  parameter                 &
             ::  MESO_REF_OMEGA = -0.463
                        ! assumed reference omega for mesoscale updraft 
                        ! (from Leary and Louze, 1980) [ Pa / sec ]
real,                       &
  parameter                 &
             ::  TPRIME_MESO_UPDRFT = 1.0    
                        ! assumed temperature excess of mesoscale updraft
                        ! over its environment [ deg K ]
real,                       &
  parameter                 &
             ::  MESO_SEP = 200.0E+02
                        ! pressure separation between base of mesoscale
                        ! updraft and top of mesoscale downdraft [ Pa ]
real,                       &
  parameter                 &
             ::  REF_PRESS = 1.0E05
                        ! reference pressure used in calculation of exner
                        ! fumction [ Pa ]
real,                       &
  parameter                 &
             ::  R_CONV_LAND  = 10.0    
                        ! assumed convective cloud droplet radius over 
                        ! land [ microns ]   
real,                       &
  parameter                 &
             ::  R_CONV_OCEAN = 16.0  
                        ! assumed convective cloud droplet radius over 
                        ! ocean [ microns ]   
real,                       &
  parameter                 &
             ::  N_LAND = 600*1.0e6 
                        ! assumed droplet number conc over land (m**-3)
real,                       &
  parameter                 &
             ::  N_OCEAN = 150*1.0e6 
                        ! assumed droplet number conc over ocean (m**-3)
real,                       &
  parameter                 &
             ::  DELZ_LAND = 500.0   
                        ! assumed cloud depth over land (m) 
real,                       &
  parameter                 &
             ::  DELZ_OCEAN = 1500.0   
                        ! assumed cloud depth over ocean (m)
real,                       &
  parameter                 &
             ::  CELL_LIQUID_EFF_DIAM_DEF = 15.0    
                        ! default cell liquid eff diameter [ microns ]
real,                       &
  parameter                 &
             ::  CELL_ICE_GENEFF_DIAM_DEF = 18.6   
                        ! default cell ice generalized effective diameter
                        ! [ microns ]
integer,                   &
  parameter         &
             ::  ANVIL_LEVELS = 6  
                        ! number of levels assumed to be in anvil clouds
real,                       &
  parameter,                &
  dimension(ANVIL_LEVELS)   &
             ::  DGEICE  = (/ 38.5, 30.72, 28.28, 25.62, 24.8, 13.3 /)
                        ! generalized effective size of hexagonal ice 
                        ! crystals, defined as in Fu (1996, J. Clim.) 
                        ! values from Table 2 of McFarquhar et al. 
                        ! (1999, JGR) are averaged over all grid boxes 
                        ! for which D_ge is defined for all altitudes 
                        ! between 9.9 and 13.2 km. index 1 at bottom of 
                        ! anvil
real,                       &
  parameter,                &
  dimension(ANVIL_LEVELS)   &
             ::  RELHT  =  (/0.0, 0.3, 0.45, 0.64, 0.76, 1.0/)
                        ! distance from anvil base, normalized by total 
                        ! anvil thickness. from Table 2 of McFarquhar et
                        ! al. (1999, JGR) for grid boxes with data 
                        ! between 9.9 and 13.2 km. index 1 at anvil 
                        ! bottom

integer,                      &
  parameter                   &
             ::  N_WATER_BUDGET = 9
                        ! number of terms in vapor budget

integer,                      &
  parameter                   &
             ::  N_ENTHALPY_BUDGET =  19 
                        ! number of terms in enthalpy budget

integer,               &
  parameter            &
             ::  N_PRECIP_PATHS = 5
                        ! number of paths precip may take from 
                        ! condensing until it reaches the ground
                        ! (liquid; liquid which freezes; liquid which
                        ! freezes and then remelts; ice; ice which 
                        ! melts)

integer,               &
  parameter            &
             ::  N_PRECIP_TYPES = 3
                        ! number of precip types (cell, cell condensate
                        ! tranmsferred to mesoscale circulation,
                        ! mesoscale condensation and deposition)



!--------------------------------------------------------------------
!   list of native mode restart versions usable by this module:
!
!   NOTE: none of the earlier versions of restart files can be used to
!         initiate an experiment with this code version due to a change 
!         in the calculation algorithm. experiments begun with this code
!         must be coldstarted, or use a native mode restart file gener-
!         ated by an experiment using this code version (restart version
!         #8), or a netcdf restart file.
!          
!   version 8 has the lag temp, vapor and pressure fields needed to cal-
!             culate the lag time value of cape. tempbl and ratpbl
!             removed. 
!
!   version 9 is reserved for the native mode restart file version cor-
!             responding to the current netcdf restart file. it is up to 
!             the user to generate the code needed to read and write this
!             version, if needed, using the subroutines read_restart and 
!             write_restart that are provided as starting points, since 
!             only netcdf restarts are currently supported.
!
!   version 10 contains donner_humidity_factor rather than 
!             donner_humidity_ratio, a change necessitated by the intro-
!             duction of the uw_conv shallow convection scheme.

integer, dimension(3)  :: restart_versions = (/ 8, 9, 10 /)


!--------------------------------------------------------------------
!   variables associated with netcdf diagnostic output from this module:
!
!   id_xxxx         indices associated with each potential netcdf 
!                   diagnostic field:
!   missing value   value used by netcdf routines if data not present
!   mod_name        module name associated with these diagnostics; used
!                   to connect these diagnostics to the diag_table
!

integer    :: id_leff
integer    :: id_cemetf_deep, id_ceefc_deep, id_cecon_deep, &
              id_cemfc_deep, id_cememf_deep, id_cememf_mod_deep, &
              id_cual_deep, id_fre_deep, id_elt_deep, &
              id_cmus_deep, id_ecds_deep, id_eces_deep, &
              id_emds_deep, id_emes_deep, id_qmes_deep,&
              id_wmps_deep, id_wmms_deep, id_tmes_deep,&
              id_dmeml_deep, id_uceml_deep, id_umeml_deep, &
              id_xice_deep,  id_dgeice_deep, id_dgeliq_deep,  &
              id_xliq_deep,    &
              id_cuqi_deep, id_cuql_deep, &
              id_plcl_deep, id_plfc_deep, id_plzb_deep, &
              id_xcape_deep, id_coin_deep,  &
              id_dcape_deep, id_qint_deep, id_a1_deep, &
              id_amax_deep, id_amos_deep, &
              id_tprea1_deep, id_ampta1_deep, &
              id_omint_deep, id_rcoa1_deep, id_detmfl_deep

integer, dimension(:), allocatable :: id_qtren1, id_qtmes1, &
                                      id_wtp1, id_qtceme, &
                                      id_total_wet_dep, &
                                      id_meso_wet_dep, id_cell_wet_dep
integer, dimension(:), allocatable :: id_qtren1_col, id_qtmes1_col, &
                                      id_wtp1_col, id_qtceme_col, &
                                      id_total_wet_dep_col, &
                                      id_meso_wet_dep_col,   &
                                      id_cell_wet_dep_col
integer, dimension(:), allocatable :: id_extremes, id_hits
integer, dimension(N_WATER_BUDGET)    :: id_water_budget, &
                                         id_ci_water_budget        
integer, dimension(N_ENTHALPY_BUDGET) :: id_enthalpy_budget,   &
                                         id_ci_enthalpy_budget
integer, dimension (N_PRECIP_PATHS, N_PRECIP_TYPES) ::            &
                                         id_precip_budget, &
                                         id_ci_precip_budget
integer   :: id_ci_prcp_heat_liq_cell, id_ci_prcp_heat_frz_cell, &
             id_ci_prcp_heat_liq_meso, id_ci_prcp_heat_frz_meso, &
             id_ci_prcp_heat_total, id_ci_prcp_total

real              :: missing_value = -999.
character(len=16) :: mod_name = 'donner_deep'

!--------------------------------------------------------------------
!   variables for column diagnostics option
!
!   arrays containing information for all requested diagnostic columns
!   (1:num_diag_pts):
!    col_diag_unit         unit numbers for each column's output file 
!    col_diag_lon          each column's longitude 
!                          [ degrees, 0 < lon < 360 ]
!    col_diag_lat          each column's latitude 
!                          [degrees, -90 < lat < 90 ]
!    col_diag_j            each column's j index (processor coordinates)
!    col_diag_i            each column's i index (processor coordinates) 
!
!    Time_col_diagnostics  time in model simulation at which to activate
!                          column diagnostics 
!

integer, dimension(:), allocatable :: col_diag_unit
real   , dimension(:), allocatable :: col_diag_lon, col_diag_lat   
integer, dimension(:), allocatable :: col_diag_j, col_diag_i        
type(time_type)                    :: Time_col_diagnostics  

!---------------------------------------------------------------------
!    derived type variables present for duration of job:
!    (see donner_types.h for documentation of their contents)
!

type(donner_param_type),       save :: Param
type(donner_column_diag_type), save :: Col_diag
type(donner_nml_type),         save :: Nml
type(donner_save_type),        save :: Don_save
type(donner_initialized_type), save :: Initialized
 
type(sounding),                save :: sd
type(uw_params),                save :: Uw_p
type(adicloud),                save :: ac
type(cplume),                  save :: cp
type(ctend ),                  save :: ct


!-----------------------------------------------------------------------
!   miscellaneous variables
!
!     module_is_initialized       module has been initialized ?
!

logical :: module_is_initialized = .false. 


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                   PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################

subroutine donner_deep_init (lonb, latb, pref, axes, Time,  &
                             tracers_in_donner, do_conservation_checks,&
                             using_unified_closure)

!---------------------------------------------------------------------
!    donner_deep_init is the constructor for donner_deep_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
real,            dimension(:,:), intent(in)   :: lonb, latb
real,            dimension(:),   intent(in)   :: pref
integer,         dimension(4),   intent(in)   :: axes
type(time_type),                 intent(in)   :: Time
logical,         dimension(:),   intent(in)   :: tracers_in_donner
logical,                         intent(in)   :: do_conservation_checks
logical,                         intent(in)   :: using_unified_closure

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      lonb         array of model longitudes on cell corners     
!                   [ radians ]
!      latb         array of model latitudes on cell corners   
!                   [ radians ]
!      pref         array of reference pressures at full levels (plus 
!                   surface value at nlev+1), based on 1013.25 hPa pstar
!                   [ Pa ]
!      axes         data axes for diagnostics
!      Time         current time [ time_type ]
!      tracers_in_donner 
!                   logical array indicating which of the activated 
!                   tracers are to be transported by donner_deep_mod
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables:

      integer                             :: unit, ierr, io
      integer                             :: idf, jdf, nlev, ntracers
      integer                             :: secs, days
      logical, dimension(size(latb,2)-1) :: do_column_diagnostics
      integer                             :: k, n, nn
!++lwh
      logical                             :: flag
      character(len=200)                  :: method_name, method_control
!--lwh
  
!-------------------------------------------------------------------
!  local variables:
!
!     unit                   unit number for nml file
!     ierr                   error return flag
!     io                     error return code
!     idf                    number of columns in the x dimension on the
!                            processors domain
!     jdf                    number of columns in the y dimension on the
!                            processors domain
!     nlev                   number of model layers 
!     ntracers               number of tracers to be transported by
!                            the donner deep convection parameterization
!     secs                   seconds component of time_type variable Time
!     days                   days component of time_type variable Time
!     do_column_diagnostics  logical array indicating which latitude rows
!                            in the processor domain contain diagnostic
!                            columns
!     k, n                   do-loop indices
!     nn                     counter of tracers transported by 
!                            donner_deep_mod
!                         
!-------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    1. READ NAMELIST AND WRITE IT TO LOG FILE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if (file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read (unit, nml=donner_deep_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'donner_deep_nml')
        enddo
10      call close_file (unit)
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() )    &
                                 write (stdlog(), nml=donner_deep_nml)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    2. DO CONSISTENCY / VALIDITY TESTS ON NML AND PARAMETER VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

       if (do_donner_cape .and. gama /= 0.0) then
         call error_mesg ('donner_deep_mod;', 'donner_deep_init: &
            & gama must be 0.0 if do_donner_cape is .true.; code for &
            & gama /=  0.0 not yet implemented', FATAL)
       endif
       if (deep_closure /= 0 ) then
         call error_mesg ('donner_deep_mod;', 'donner_deep_init: &
              & deep_closure must be 0; code for &
              & deep_closure /=  0 not yet implemented', FATAL)
       endif
       if (do_rh_trig .and. do_donner_closure) then
         call error_mesg ('donner_deep_mod;', 'donner_deep_init: &
             & do_rh_trig must be .false. for donner full  &
               &parameterization; its use not yet implemented', FATAL)
        endif

!---------------------------------------------------------------------
!    check for a valid value of donner_deep_freq. 
!---------------------------------------------------------------------
      if (donner_deep_freq > 86400) then
        call error_mesg ( 'donner_deep_mod', 'donner_deep_init: &
         & donner convection must be called at least once per day', &
                                                           FATAL)
      else if (donner_deep_freq <= 0) then
        call error_mesg ( 'donner_deep_mod', 'donner_deep_init: &
          & a positive value must be assigned to donner_deep_freq', &
                                                            FATAL)
      endif

!---------------------------------------------------------------------
!    check for valid value of entrainment_constant_source.
!---------------------------------------------------------------------
      if (trim(entrainment_constant_source) == 'gate' .or. &
          trim(entrainment_constant_source) == 'kep' ) then
      else
        call error_mesg ('donner_deep_mod', 'donner_deep_init: &
        & invalid string for nml variable entrainment_constant_source', &
                                                             FATAL)
      endif

!---------------------------------------------------------------------
!    test that PSTOP is smaller than UPPER_LIMIT_FOR_LFC.
!---------------------------------------------------------------------
      if (PSTOP > UPPER_LIMIT_FOR_LFC) then
        call error_mesg ('donner_deep_mod', 'donner_deep_init: &
           & pstop must be above the upper limit of &
                                &the level of free convection', FATAL)
      endif

!---------------------------------------------------------------------
!    test that cell_liquid_size_type has been validly specified, and if
!    it is specified as 'input', an appropriate input value has been
!    supplied.
!---------------------------------------------------------------------
      if (trim(cell_liquid_size_type) == 'input') then
        Initialized%do_input_cell_liquid_size = .true.
        Initialized%do_bower_cell_liquid_size = .false.
        if (cell_liquid_eff_diam_input < 0.0) then
          call error_mesg ('donner_deep_mod', 'donner_deep_init: &
            & cell liquid size must be input, but no value supplied', &
                                                               FATAL)
        endif
      else if (trim(cell_liquid_size_type) == 'bower') then
        Initialized%do_input_cell_liquid_size = .false.
        Initialized%do_bower_cell_liquid_size = .true.
      else
        call error_mesg ( 'donner_deep_mod', 'donner_deep_init: &
           & cell_liquid_size_type must be either input or bower', &
                                                                FATAL)
      endif

!---------------------------------------------------------------------
!    test that cell_ice_size_type has been validly specified, and if
!    specified as 'input', that cell_ice_geneff_diam_input has also 
!    been appropriately defined.
!---------------------------------------------------------------------
      if (trim(cell_ice_size_type) == 'input') then
        Initialized%do_input_cell_ice_size = .true.
        Initialized%do_default_cell_ice_size = .false.
        if (cell_ice_geneff_diam_input <= 0.0) then
          call error_mesg ('donner_deep_mod', 'donner_deep_init: &
               & must define a nonnegative generalized effective '//&
                'diameter for ice when cell_ice_size_type is input', &
                                                                 FATAL)
        endif
      else if (trim(cell_ice_size_type) == 'default') then
        Initialized%do_input_cell_ice_size = .false.
        Initialized%do_default_cell_ice_size = .true.
      else
        call error_mesg ( 'donner_deep_init', 'donner_deep_init: &
             & cell_ice_size_type must be input or default',  FATAL)
      endif


!---------------------------------------------------------------------
!    place the logical input argument indicating whether the cloud 
!    base mass flux calculated by uw_conv_mod is also to be used 
!    in defining the closure for donner deep convection in the 
!    donner_initialized_type variable Initialized.
!    place the conservation check flag in the Initialized variable.
!---------------------------------------------------------------------
      Initialized%using_unified_closure = using_unified_closure
      Initialized%do_conservation_checks = do_conservation_checks

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    3. PROCESS TRACERS THAT ARE TO BE TRANSPORTED BY THE DONNER DEEP
!       CONVECTION PARAMETERIZATION.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    determine how many tracers are to be transported by donner_deep 
!    convection. allocate arrays to contain their names and units for use
!    with diagnostics and restarts. define a logical variable indicating
!    if any tracers are to be so transported. obtain the tracer names and
!    units.
!---------------------------------------------------------------------
      ntracers = count(tracers_in_donner)
      allocate ( Don_save%tracername   (ntracers) )
      allocate ( Don_save%tracer_units (ntracers) )
!++lwh
      allocate ( Initialized%wetdep(ntracers) )
!--lwh
      if (ntracers > 0) then
        Initialized%do_donner_tracer = .true.
        nn = 1
        do n=1,size(tracers_in_donner(:))
          if (tracers_in_donner(n)) then
            call get_tracer_names (MODEL_ATMOS, n,  &
                                   name = Don_save%tracername(nn), &
                                   units = Don_save%tracer_units(nn))
!++lwh
            Initialized%wetdep(nn)%units = Don_save%tracer_units(nn)
            flag = query_method( 'wet_deposition', MODEL_ATMOS, n, &
                                 method_name, method_control )
            call get_wetdep_param( method_name, method_control, &
                                   Initialized%wetdep(nn)%scheme, &
                                   Initialized%wetdep(nn)%Henry_constant, &
                                   Initialized%wetdep(nn)%Henry_variable, &
                                   Initialized%wetdep(nn)%frac_in_cloud, &
                                   Initialized%wetdep(nn)%alpha_r, &
                                   Initialized%wetdep(nn)%alpha_s )
            Initialized%wetdep(nn)%scheme = lowercase( Initialized%wetdep(nn)%scheme )
!-lwh
            nn = nn + 1
          endif
        end do
      else
        Initialized%do_donner_tracer = .false.
      endif
      

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    4. DEFINE PROCESSOR DIMENSIONS AND ALLOCATE SPACE FOR MODULE 
!       VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-------------------------------------------------------------------
!    define the grid dimensions. idf and jdf are the (i,j) dimensions of
!    the domain on this processor, nlev is the number of model layers.
!-------------------------------------------------------------------
      nlev = size(pref(:)) - 1
      idf  = size(lonb,1) - 1
      jdf  = size(latb,2) - 1

!---------------------------------------------------------------------
!    initialize the points processed counter. define the total number 
!    of columns present on the processor. 
!---------------------------------------------------------------------
      Initialized%pts_processed_conv   = 0
      Initialized%total_pts = idf*jdf

!--------------------------------------------------------------------
!    allocate module variables that will be saved across timesteps.
!    these are stored in the derived-type variable Don_save. see 
!    donner_types.h for description of these variables.
!--------------------------------------------------------------------
      allocate ( Don_save%cemetf             (idf, jdf, nlev ) )
      allocate ( Don_save%lag_temp           (idf, jdf, nlev ) )
      allocate ( Don_save%lag_vapor          (idf, jdf, nlev ) )
      allocate ( Don_save%lag_press          (idf, jdf, nlev ) )
      allocate ( Don_save%cememf             (idf, jdf, nlev ) )
      allocate ( Don_save%mass_flux          (idf, jdf, nlev ) )
      allocate ( Don_save%cell_up_mass_flux  (idf, jdf, nlev+1 ) )
      allocate ( Don_save%det_mass_flux      (idf, jdf, nlev ) )
      allocate ( Don_save%dql_strat          (idf, jdf, nlev ) )
      allocate ( Don_save%dqi_strat          (idf, jdf, nlev ) )
      allocate ( Don_save%dqa_strat          (idf, jdf, nlev ) )
      allocate ( Don_save%humidity_area      (idf, jdf, nlev ) )
      allocate ( Don_save%humidity_factor    (idf, jdf, nlev ) )
      allocate ( Don_save%tracer_tends       (idf, jdf, nlev, ntracers) )
      allocate ( Don_save%parcel_disp        (idf, jdf ) )
      allocate ( Don_save%tprea1             (idf, jdf ) )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    4. INITIALIZE THE NETCDF OUTPUT VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!--------------------------------------------------------------------
!    activate the netcdf diagnostic fields.
!-------------------------------------------------------------------
      call register_fields (Time, axes)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    5. PROCESS THE RESTART FILE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!--------------------------------------------------------------------
!    if a netcdf restart file is present, call read_restart_nc to read 
!    it.
!--------------------------------------------------------------------
      if (file_exist ('INPUT/donner_deep.res.nc') ) then
        Initialized%coldstart= .false.
        call read_restart_nc (ntracers)

!--------------------------------------------------------------------
!    if a native mode restart file is present, call read_restart 
!    to read it.
!--------------------------------------------------------------------
      else if (file_exist ('INPUT/donner_deep.res') ) then
        Initialized%coldstart= .false.
        call read_restart (ntracers, Time)

!--------------------------------------------------------------------
!    if no restart file is present, call subroutine process_coldstart
!    to define the needed variables.
!--------------------------------------------------------------------
      else
        call process_coldstart (Time)
      endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    6. INITIALIZE VARIABLES NEEDED FOR COLUMN_DIAGNOSTICS_MOD OUTPUT.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    define the total number of columns for which diagnostics
!    are desired.
!---------------------------------------------------------------------
      Col_diag%num_diag_pts = num_diag_pts_ij + num_diag_pts_latlon

!---------------------------------------------------------------------
!    initialize the value of the k index associated with diagnostics
!    cutoff.
!---------------------------------------------------------------------
      Col_diag%kstart = -99

!---------------------------------------------------------------------
!    if any diagnostics are requested, perform various consistency
!    checks.
!---------------------------------------------------------------------
      if (Col_diag%num_diag_pts > 0) then

!---------------------------------------------------------------------
!    check that array dimensions are sufficiently large for the number 
!    of columns requested.
!---------------------------------------------------------------------
        if (Col_diag%num_diag_pts > MAX_PTS) then
          call error_mesg ('donner_deep_mod', 'donner_deep_init: &
         &must reset MAX_PTS or reduce number of diagnostic points', &
                                                           FATAL)  
        endif

!---------------------------------------------------------------------
!    check that the specified time at which diagnostics are to be 
!    activated has been specified.
!---------------------------------------------------------------------
        do n=1,3
          if (diagnostics_start_time(n) == 0) then
            call error_mesg ('donner_deep_mod', 'donner_deep_init:&
             &year, month and/or day invalidly specified for column '//&
                  'diagnostics starting time', FATAL)
          endif
        end do

!---------------------------------------------------------------------
!    define a time_type variable indicating the requested time to begin
!    outputting diagnostics.
!---------------------------------------------------------------------
        Time_col_diagnostics = set_date (diagnostics_start_time(1), &
                                         diagnostics_start_time(2), &   
                                         diagnostics_start_time(3), &   
                                         diagnostics_start_time(4), &   
                                         diagnostics_start_time(5), &   
                                         diagnostics_start_time(6) )    

!---------------------------------------------------------------------
!    allocate space for the arrays used to specify the diagnostics 
!    columns and the output units. initialize the arrays with bogus
!    values.
!---------------------------------------------------------------------
        allocate (col_diag_unit    (Col_diag%num_diag_pts) )
        allocate (col_diag_lon     (Col_diag%num_diag_pts) )
        allocate (col_diag_lat     (Col_diag%num_diag_pts) )
        allocate (col_diag_i       (Col_diag%num_diag_pts) )
        allocate (col_diag_j       (Col_diag%num_diag_pts) )
        col_diag_unit  = -1
        col_diag_lon   = -1.0
        col_diag_lat   = -1.0
        col_diag_i     = -1
        col_diag_j     = -1

!---------------------------------------------------------------------
!    call initialize_diagnostic_columns to determine the locations 
!    (i,j,lat and lon) of any diagnostic columns in this processor's
!    space and to open output files for the diagnostics.
!---------------------------------------------------------------------
        call initialize_diagnostic_columns   &
                     (mod_name, num_diag_pts_latlon, num_diag_pts_ij, &
                      i_coords_gl, j_coords_gl, lat_coords_gl, &
                      lon_coords_gl, lonb(:,1), latb(1,:),  &
                      do_column_diagnostics, &
                      col_diag_lon, col_diag_lat, col_diag_i,  &
                      col_diag_j, col_diag_unit)

!---------------------------------------------------------------------
!    verify that requested pressure cutoff for column diagnostics output
!    is valid. define the model k index which corresponds (kstart).
!---------------------------------------------------------------------
        do k=1,size(pref(:))
          if (pref(k) >= diagnostics_pressure_cutoff) then
            Col_diag%kstart = k
            exit
          endif
        end do

!----------------------------------------------------------------------
!    if the specified pressure is larger than any pressure level in the
!    model grid, write an error message.
!----------------------------------------------------------------------
        if (Col_diag%kstart == -99) then
          call error_mesg ( 'donner_deep_mod', 'donner_deep_init: &
           &diagnostics_pressure_cutoff is higher than pressure at '//&
                                     'any model level', FATAL)
        endif

!----------------------------------------------------------------------
!   if column diagnostics is not requested, define the components of
!   Col_diag that will be needed.
!----------------------------------------------------------------------
      else
        Col_diag%in_diagnostics_window = .false.
        Col_diag%ncols_in_window = 0
      endif

!----------------------------------------------------------------------
!    allocate space for the array elements of the donner_column_diag_type
!    variable Col_diag. These arrays remain for the life of the job and
!    will be defined for each physics window as it is entered.
!----------------------------------------------------------------------
      allocate (Col_diag%i_dc(Col_diag%num_diag_pts))
      allocate (Col_diag%j_dc(Col_diag%num_diag_pts))
      allocate (Col_diag%unit_dc(Col_diag%num_diag_pts))
      allocate (Col_diag%igl_dc(Col_diag%num_diag_pts))
      allocate (Col_diag%jgl_dc(Col_diag%num_diag_pts))

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    7. FILL THE DONNER_PARAM_TYPE VARIABLE WITH VALUES THAT HAVE BEEN 
!       DEFINED HERE.                  
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------
!    define the components of Param that come from constants_mod. see 
!    donner_types.h for their definitions.
!----------------------------------------------------------------------
      Param%dens_h2o        = DENS_H2O
      Param%rdgas           = RDGAS
      Param%grav            = GRAV
      Param%cp_air          = CP_AIR  
      Param%pie             = PIE
      Param%kappa           = KAPPA
      Param%rvgas           = RVGAS
      Param%seconds_per_day = SECONDS_PER_DAY
      Param%hlv             = HLV
      Param%hlf             = HLF
      Param%hls             = HLS
      Param%kelvin          = KELVIN

!----------------------------------------------------------------------
!    store the parameters defined in this module into the 
!    donner_parameter_type variables Param. these variables are defined
!    above.
!----------------------------------------------------------------------
      Param%cp_vapor                = CP_VAPOR
      Param%parcel_dp               = PARCEL_DP
      Param%upper_limit_for_lfc     = UPPER_LIMIT_FOR_LFC
      Param%pstop                   = PSTOP
      Param%cld_base_vert_vel       = CLD_BASE_VERT_VEL
      Param%dp_of_cloud_model       = DP_OF_CLOUD_MODEL
      Param%cloud_base_radius       = CLOUD_BASE_RADIUS
      Param%wdet                    = WDET
      Param%rbound                  = RBOUND
      Param%wbound                  = WBOUND
      Param%freeze_fraction         = FREEZE_FRACTION
      Param%virt_mass_co            = VIRT_MASS_CO
      Param%pdeep_mc                = PDEEP_MC
      Param%tr_insert_time          = TR_INSERT_TIME
      Param%autoconv_rate           = AUTOCONV_RATE
      Param%autoconv_threshold      = AUTOCONV_THRESHOLD
      Param%tfre                    = TFRE
      Param%dfre                    = DFRE
      Param%evap_in_downdrafts      = EVAP_IN_DOWNDRAFTS
      Param%evap_in_environ         = EVAP_IN_ENVIRON      
      Param%entrained_into_meso     = ENTRAINED_INTO_MESO
      Param%d622                    = D622
      Param%d608                    = D608
      Param%upper_limit_for_lcl     = UPPER_LIMIT_FOR_LCL
      Param%tmin                    = TMIN
      Param%anvil_precip_efficiency = ANVIL_PRECIP_EFFICIENCY
      Param%meso_lifetime           = MESO_LIFETIME
      Param%meso_ref_omega          = MESO_REF_OMEGA
      Param%tprime_meso_updrft      = TPRIME_MESO_UPDRFT
      Param%meso_sep                = MESO_SEP
      Param%ref_press               = REF_PRESS
      Param%meso_down_evap_fraction = MESO_DOWN_EVAP_FRACTION
      Param%meso_up_evap_fraction   = MESO_UP_EVAP_FRACTION
      Param%istart                  = ISTART

      Param%max_entrainment_constant_gate = MAX_ENTRAINMENT_CONSTANT_GATE
      Param%max_entrainment_constant_kep  = MAX_ENTRAINMENT_CONSTANT_KEP
      Param%pdeep_cv                      = PDEEP_CV
      Param%cdeep_cv                      = CDEEP_CV
      Param%kpar                          = KPAR
      Param%r_conv_land                   = R_CONV_LAND
      Param%r_conv_ocean                  = R_CONV_OCEAN 
      Param%n_land                        = N_LAND
      Param%n_ocean                       = N_OCEAN
      Param%delz_land                     = DELZ_LAND
      Param%delz_ocean                    = DELZ_OCEAN
      Param%cell_liquid_eff_diam_def      = CELL_LIQUID_EFF_DIAM_DEF 
      Param%cell_ice_geneff_diam_def      = CELL_ICE_GENEFF_DIAM_DEF
      Param%anvil_levels                  = ANVIL_LEVELS 

      allocate (Param%arat(kpar))
      allocate (Param%ensemble_entrain_factors_gate(kpar))
      allocate (Param%ensemble_entrain_factors_kep(kpar))
      Param%arat                          = ARAT
      Param%ensemble_entrain_factors_gate = erat
      Param%ensemble_entrain_factors_kep  = ENSEMBLE_ENTRAIN_FACTORS_KEP

      allocate (Param%dgeice (ANVIL_LEVELS))
      allocate (Param%relht  (ANVIL_LEVELS))
      Param%dgeice  = DGEICE               
      Param%relht   = RELHT                

      call sd_init_k (nlev, ntracers, sd);
      call uw_params_init_k (Param%hlv, Param%hls, Param%hlf, &
          Param%cp_air, Param%grav, Param%kappa, Param%rdgas,  &
          Param%ref_press, Param%d622,Param%d608, Param%kelvin -160., & 
          Param%kelvin + 100. ,Uw_p)
      call ac_init_k (nlev, ac);
      call cp_init_k (nlev, ntracers, cp);
      call ct_init_k (nlev, ntracers, ct)
      call sat_vapor_pres_init
      call exn_init_k (Uw_p)
      call findt_init_k (Uw_p)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    8. STORE THE NAMELIST VARIABLES THAT NEED TO BE MADE AVAILABLE 
!       OUTSIDE OF THIS MODULE INTO THE DONNER_NML_TYPE VARIABLE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      Nml%parcel_launch_level         = parcel_launch_level
      Nml%allow_mesoscale_circulation = allow_mesoscale_circulation
      Nml%do_donner_cape              = do_donner_cape    !miz
      Nml%do_donner_plume             = do_donner_plume   !miz
      Nml%do_donner_closure           = do_donner_closure !miz
      Nml%do_dcape                    = do_dcape          !miz
      Nml%do_lands                    = do_lands          !miz
      Nml%tau                         = tau               !miz
      Nml%cape0                       = cape0             !miz
      Nml%rhavg0                      = rhavg0            !miz
      Nml%plev0                       = plev0             !miz
      Nml%do_rh_trig                  = do_rh_trig        !miz
      Nml%do_capetau_land             = do_capetau_land   !miz
      Nml%pblht0                      = pblht0            !miz
      Nml%tke0                        = tke0              !miz
      Nml%lofactor0                   = lofactor0         !miz
      Nml%deephgt0                    = deephgt0          !miz
      Nml%lochoice                    = lochoice          !miz
      Nml%deep_closure                = deep_closure      !miz
      Nml%gama                        = gama              !miz
      Nml%do_ice                      = do_ice            !miz
      Nml%atopevap                    = atopevap          !miz
      Nml%do_donner_lscloud           = do_donner_lscloud !miz
      Nml%auto_rate                   = auto_rate         !miz
      Nml%auto_th                     = auto_th           !miz
      Nml%frac                        = frac              !miz
      Nml%ttend_max                   = ttend_max         !miz
      Nml%use_llift_criteria          = use_llift_criteria
      Nml%use_pdeep_cv                = use_pdeep_cv
      Nml%entrainment_constant_source = entrainment_constant_source
      Nml%donner_deep_freq            = donner_deep_freq             
      Nml%model_levels_in_sfcbl       = model_levels_in_sfcbl        
      Nml%cell_liquid_size_type       = cell_liquid_size_type 
      Nml%cell_ice_size_type          = cell_ice_size_type
      Nml%cell_liquid_eff_diam_input  = cell_liquid_eff_diam_input
      Nml%cell_ice_geneff_diam_input  = cell_ice_geneff_diam_input
      Nml%meso_liquid_eff_diam_input  = meso_liquid_eff_diam_input
      Nml%do_average                  = do_average
      Nml%use_memphis_size_limits     = use_memphis_size_limits
      Nml%wmin_ratio                  = wmin_ratio
      Nml%do_freezing_for_cape         = do_freezing_for_cape
      Nml%tfre_for_cape               = tfre_for_cape
      Nml%dfre_for_cape               = dfre_for_cape
      Nml%rmuz_for_cape               = rmuz_for_cape
      Nml%do_freezing_for_closure     = do_freezing_for_closure
      Nml%tfre_for_closure            = tfre_for_closure
      Nml%dfre_for_closure            = dfre_for_closure
      Nml%rmuz_for_closure            = rmuz_for_closure
      Nml%do_budget_analysis          = do_budget_analysis
      Nml%force_internal_enthalpy_conservation =  &
                                 force_internal_enthalpy_conservation


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    9. SET UP CODE TO MONITOR SELECTED OUTPUT VARIABLES.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

       call process_monitors (idf, jdf, nlev, ntracers, axes, Time)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!   10. END OF SUBROUTINE.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!--------------------------------------------------------------------
!    set flag to indicate that donner_deep_mod has been initialized.
!--------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------



end subroutine donner_deep_init



!###################################################################

subroutine donner_deep (is, ie, js, je, dt, temp, mixing_ratio, pfull, &
                        phalf, zfull, zhalf, omega, pblht, tkemiz, &
                        qstar, cush, coldT, land, sfc_sh_flux,  &
                        sfc_vapor_flux,&               !miz
                        tr_flux, tracers, Time, cbmf, cell_cld_frac,  &
                        cell_liq_amt, cell_liq_size, cell_ice_amt,   &
                        cell_ice_size, cell_droplet_number, &
                        meso_cld_frac, meso_liq_amt, &
                        meso_liq_size, meso_ice_amt, meso_ice_size,  &
                        meso_droplet_number, &
                        nsum, precip, delta_temp, delta_vapor, detf, &
                        uceml_inter, mtot, donner_humidity_area,    &
                        donner_humidity_factor, qtrtnd, &
                        lheat_precip, vert_motion,        &
                        total_precip, liquid_precip, frozen_precip, &
                        qlin, qiin, qain,              &      ! optional
                        delta_ql, delta_qi, delta_qa)         ! optional
                        
!-------------------------------------------------------------------
!    donner_deep is the prognostic driver subroutine of donner_deep_mod.
!    it takes as input the temperature (temp), vapor mixing ratio 
!    (mixing_ratio), pressure at full and half-levels (pfull, phalf),
!    vertical velocity at full levels (omega), the large scale cloud 
!    variables (qlin, qiin, qain), the land fraction (land),  the heat 
!    (sfc_sh_flux) , moisture (sfc_vapor_flux) and tracer (tr_flux) 
!    fluxes across the surface that are to be seen by this parameter-
!    ization, the tracers to be transported by the donner convection
!    parameterization (tracers), and the current time (as time_type 
!    variable Time). the routine returns the precipitation (precip),
!    increments to the temperature (delta_temp) and mixing ratio 
!    (delta_vapor), the detrained mass flux (detf), upward cell mass 
!    flux at interface levels  (uceml_inter) and total mass flux at full
!    levels (mtot), two arrays needed to connect the donner convection 
!    and strat cloud parameterizations (donner_humidity_area, 
!    donner_humidity_ratio), increments to the cloudwater (delta_ql), 
!    cloudice (delta_qi) and cloud area (delta_qa) fields and tendencies
!    for those tracers that are to be transported by the donner convect-
!    ion parameterization (qtrtnd). there are an additional eleven arrays
!    defining the donner scheme cloud characteristics needed by the rad-
!    iation package, which are passed in and updated on donner calcul-
!    ation steps.
!-------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                      intent(in)    :: is, ie, js, je
real,                         intent(in)    :: dt
real, dimension(:,:,:),       intent(in)    :: temp, mixing_ratio, &
                                               pfull, phalf, zfull, zhalf, omega
real, dimension(:,:),         intent(in)    :: pblht, tkemiz, qstar,cush
real, dimension(:,:),         intent(in)    :: land
logical, dimension(:,:),      intent(in)    :: coldT
real, dimension(:,:),         intent(in)    :: sfc_sh_flux, &
                                               sfc_vapor_flux
real, dimension(:,:,:),       intent(in)    :: tr_flux 
real, dimension(:,:,:,:),     intent(in)    :: tracers 
type(time_type),              intent(in)    :: Time
real, dimension(:,:),         intent(inout)    :: cbmf              
real, dimension(:,:,:),       intent(inout) :: cell_cld_frac,  &
                                               cell_liq_amt,  &
                                               cell_liq_size, &
                                               cell_ice_amt,  &
                                               cell_ice_size, &
                                           cell_droplet_number, &
                                               meso_cld_frac,  &
                                               meso_liq_amt, &
                                               meso_liq_size, &
                                               meso_ice_amt,   &
                                               meso_ice_size, &
                                           meso_droplet_number
integer, dimension(:,:),      intent(inout) :: nsum
real, dimension(:,:),         intent(out)   :: precip, &
                                               lheat_precip, &
                                               vert_motion, &
                                               total_precip
real, dimension(:,:,:),       intent(out)   :: delta_temp, delta_vapor,&
                                               detf, uceml_inter, mtot, &
                                               donner_humidity_area,&
                                               donner_humidity_factor, &
                                               liquid_precip, &
                                               frozen_precip
real, dimension(:,:,:,:),     intent(out)   :: qtrtnd 
real, dimension(:,:,:),       intent(in),                &
                                   optional :: qlin, qiin, qain
real, dimension(:,:,:),       intent(out),               &
                                   optional :: delta_ql, delta_qi, &
                                               delta_qa

!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     dt             physics time step [ sec ]
!     temp           temperature field at model levels [ deg K ]
!     mixing_ratio   vapor mixing ratio field at model levels 
!                    [ kg(h20) / kg(dry air) ]
!     pfull          pressure field on model full levels [ Pa ]
!     phalf          pressure field at half-levels 1:nlev+1  [ Pa ]
!     omega          model omega field at model full levels [ Pa / sec ]
!     qlin           large-scale cloud liquid specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qiin           large-scale cloud ice specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qain           large-scale cloud fraction  
!                    [ fraction ]
!     land           fraction of grid box covered by land
!                    [ fraction ]
!     sfc_sh_flux    sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     tr_flux        surface flux of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     tracers        tracer mixing ratios
!                    [ kg(tracer) / kg (dry air) ]
!     Time           current time (time_type)
!
!   intent(out) variables:
!
!     precip         precipitation generated by deep convection
!                    [ kg(h2o) / m**2 ]
!     delta_temp     temperature increment due to deep convection 
!                    [ deg K ]
!     delta_vapor    water vapor mixing ratio increment due to deep 
!                    convection [ kg(h2o) / kg (dry air) ]
!     detf           detrained cell mass flux at model levels 
!                    [ (kg / (m**2 sec) ) ]
!     uceml_inter    upward cell mass flux at interface levels 
!                    [ (kg / (m**2 sec) ) ]
!     mtot           mass flux at model full levels, convective plus 
!                    mesoscale, due to donner_deep_mod 
!                    [ (kg / (m**2 sec) ) ]
!     donner_humidity_area
!                    fraction of grid box in which humidity is affected
!                    by the deep convection, defined as 0.0 below cloud
!                    base and above the mesoscale updraft, and as the
!                    sum of the cell and mesoscale cloud areas in 
!                    between. it is used in strat_cloud_mod to determine
!                    the large-scale specific humidity field for the
!                    grid box. DO NOT use for radiation calculation,
!                    since not all of this area includes condensate.
!                    [ fraction ]
!     donner_humidity_ratio
!                    ratio of large-scale specific humidity to specific 
!                    humidity in environment outside convective system
!                    [ dimensionless ]
!     delta_ql       cloud water specific humidity increment due to 
!                    deep convection over the timestep
!                    [ kg (h2o) / kg (moist air) ]
!     delta_qi       cloud ice specific humidity increment due to deep 
!                    convection over the timestep 
!                    [ kg (h2o) / kg (moist air) ]
!     delta_qa       cloud area increment due to deep convection
!                    over the time step [ fraction ]
!     qtrtnd         tracer time tendencies due to deep convection
!                    during the time step
!                    [ kg(tracer) / (kg (dry air) sec) ]
!
!   intent(inout) variables:
!
!     cell_cld_frac  fractional coverage of convective cells in
!                    grid box [ dimensionless ]
!     cell_liq_amt   liquid water content of convective cells
!                    [ kg(h2o) / kg(air) ]
!     cell_liq_size  assumed effective size of cell liquid drops
!                    [ microns ]
!     cell_ice_amt   ice water content of cells
!                    [ kg(h2o) / kg(air) ]
!     cell_ice_size  generalized effective diameter for ice in
!                    convective cells [ microns ]
!     meso_cld_frac  fractional area of mesoscale clouds in grid
!                    box [ dimensionless ]
!     meso_liq_amt   liquid water content in mesoscale clouds
!                    [ kg(h2o) / kg(air) ]
!     meso_liq_size  assumed effective size of mesoscale drops
!                    [ microns ]
!     meso_ice_amt   ice water content of mesoscale elements
!                    [ kg(h2o) / kg(air) ]
!     meso_ice_size  generalized ice effective size for anvil ice
!                    [ microns ]
!     nsum           number of time levels over which the above variables
!                    have so far been summed
!
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!    local variables:

      real,    dimension (size(temp,1), size(temp,2), size(temp,3)) :: &
                       temperature_forcing, moisture_forcing, pmass, &
                       qlin_arg, qiin_arg, qain_arg, delta_ql_arg, & 
                       delta_qi_arg, delta_qa_arg

      real,    dimension (size(temp,1), size(temp,2)) ::   parcel_rise

      type(donner_conv_type)            :: Don_conv
      type(donner_budgets_type)         :: Don_budgets
      type(donner_cape_type)            :: Don_cape
      type(donner_rad_type)             :: Don_rad
      character(len=128)                :: ermesg
      integer                           :: isize, jsize, nlev_lsm
      integer                           :: ntr, me
      logical                           :: calc_conv_on_this_step 
      logical                           :: cloud_tracers_present
      integer                           :: num_cld_tracers
      integer                           :: i, j, k, n   
      logical                           :: used

!--------------------------------------------------------------------
!   local variables:
!
!     temperature_forcing  temperature tendency due to donner convection
!                          [ deg K / sec ]
!     moisture_forcing     vapor mixing ratio tendency due to donner 
!                          convection [ kg(h2o) / (kg(dry air) sec ) ]
!     pmass                mass per unit area within the grid box
!                          [ kg (air) / (m**2) ]
!     parcel_rise          accumulated vertical displacement of a 
!                          near-surface parcel as a result of the lowest
!                          model level omega field [ Pa ]
!     total_precip         total precipitation rate produced by the
!                          donner parameterization [ mm / day ]
!     exit_flag            logical array indicating whether deep conv-
!                          ection exists in a column
!     Don_conv             donner_convection_type derived type variable 
!                          containing diagnostics and intermediate
!                          results describing the nature of the convec-
!                          tion produced by the donner parameterization
!     Don_cape             donner_cape type derived type variable con-
!                          taining diagnostics and intermediate results
!                          related to the cape calculation associated 
!                          with the donner convection parameterization
!     Don_rad              donner_rad_type derived type variable used
!                          to hold those fields needed to connect the
!                          donner deep convection parameterization and
!                          the model radiation package
!     ermesg               character string containing any error message
!                          that is returned from a kernel subroutine
!     isize                x-direction size of the current physics window
!     isize, jsize         y-direction size of the current physics window
!     nlev_lsm             number of model layers in large-scale model
!     ntr                  number of tracers to be transported by donner
!                          convection 
!     me                   local pe number
!     calc_conv_on_this_step 
!                          is this a step on which to calculate 
!                          convection ?
!     k                    do-loop index
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    check that the module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('donner_deep_mod', 'donner_deep: &
             &donner_deep_init was not called before subroutine   &
                                                  &donner_deep', FATAL)
      endif

!----------------------------------------------------------------------
!    determine if the arguments needed when run with the strat_cloud_mod 
!    are present; set cloud_tracers_present appropriately.
!----------------------------------------------------------------------
      num_cld_tracers = count( (/present(qlin), present(qiin),   &
                                 present(qain), present(delta_ql), &
                                 present(delta_qi),present(delta_qa)/) )
      if (num_cld_tracers == 0) then
        cloud_tracers_present = .false.
        qlin_arg = 0.
        qiin_arg = 0.
        qain_arg = 0.
      else if (num_cld_tracers == 6) then
        cloud_tracers_present = .true.
        qlin_arg = qlin 
        qiin_arg = qiin
        qain_arg = qain
      else
        call error_mesg ('donner_deep_mod','donner_deep: &
                        &Either none or all of the cloud tracers '// &
                         'and their tendencies must be present',FATAL)
      endif

!--------------------------------------------------------------------
!    if column diagnostics have been requested for any column, call 
!    donner_column_control to define the components of the 
!    donner_column_diag_type variable for the diagnostic columns in this 
!    window. if column diagnostics have not been requested, the needed
!    variables so indicating have already been set.
!--------------------------------------------------------------------
      if (Col_diag%num_diag_pts > 0) then
        call donner_column_control (is, ie, js, je, Time)
      endif

!-------------------------------------------------------------------
!    define the dimensions for the variables in this physics window.
!    define the pe number of the current pe.
!-------------------------------------------------------------------
      isize     = ie - is + 1
      jsize     = je - js + 1
      nlev_lsm  = size(temp,3)
      ntr       = size(tracers,4) 
      me        = mpp_pe()
      Don_budgets%n_water_budget      = N_WATER_BUDGET
      Don_budgets%n_enthalpy_budget   = N_ENTHALPY_BUDGET
      Don_budgets%n_precip_paths      = N_PRECIP_PATHS     
      Don_budgets%n_precip_types      = N_PRECIP_TYPES     

!-----------------------------------------------------------------------
!    call the kernel subroutine don_d_donner_deep_k to obtain the
!    output fields resulting from the donner deep convection parameter-
!    ization.
!-----------------------------------------------------------------------
      call don_d_donner_deep_k   &
           (is, ie, js, je, isize, jsize, nlev_lsm, NLEV_HIRES, ntr, me,&
            cloud_tracers_present,  cbmf,    &
            dt, Param, Nml, temp, mixing_ratio, pfull,    &
            phalf, zfull, zhalf, omega, pblht, tkemiz, qstar, cush, coldT,&
!           qlin, qiin, qain, land, sfc_sh_flux, sfc_vapor_flux,    &
            qlin_arg, qiin_arg, qain_arg, land, sfc_sh_flux,  &
            sfc_vapor_flux,    &
            tr_flux, tracers, cell_cld_frac, cell_liq_amt,      &
            cell_liq_size, cell_ice_amt, cell_ice_size,   &
            cell_droplet_number, meso_cld_frac,  &
            meso_liq_amt, meso_liq_size, meso_ice_amt, meso_ice_size,  &
            meso_droplet_number, &
            nsum, precip, delta_temp, delta_vapor, detf, uceml_inter,  &
            mtot, donner_humidity_area, donner_humidity_factor, &
            total_precip, temperature_forcing, moisture_forcing,    &
!           parcel_rise, delta_ql, delta_qi, delta_qa, qtrtnd,         &
            parcel_rise, delta_ql_arg, delta_qi_arg, delta_qa_arg,   &
            qtrtnd,         &
            calc_conv_on_this_step, ermesg, Initialized, Col_diag,   &
            Don_rad, Don_conv, Don_cape, Don_save, &!miz
            sd, Uw_p, ac, cp, ct,  Don_budgets)

!----------------------------------------------------------------------
!    if strat_cloud is active, move the output arguments into the proper
!    locations.
!----------------------------------------------------------------------
      if (cloud_tracers_present) then
        delta_ql = delta_ql_arg
        delta_qi = delta_qi_arg
        delta_qa = delta_qa_arg
      endif

      if (Initialized%do_conservation_checks .or.   &
                                          Nml%do_budget_analysis) then
        lheat_precip = Don_budgets%lheat_precip
        vert_motion = Don_budgets%vert_motion
      else
        lheat_precip = 0.
        vert_motion = 0.
      endif
      liquid_precip = Don_budgets%liq_prcp
      frozen_precip = Don_budgets%frz_prcp

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, process the error message.
!!!  HOW TO DISTINGUISH FATAL, WARNING, NOTE ??
!    FOR NOW, ALL messages considered FATAL.
!----------------------------------------------------------------------
      if (trim(ermesg) /= ' ') then
        print *, 'ermesg', ermesg, me
        call error_mesg ('donner_deep_mod', ermesg, FATAL)
      endif

!---------------------------------------------------------------------
!    if this is a calculation step for donner_deep, define a mass
!    weighting factor (mass per unit area) needed for some of the netcdf
!    diagnostics (pmass). call donner_deep_netcdf to send the requested 
!    diagnostic data to the diag_manager for output.
!---------------------------------------------------------------------
      if (calc_conv_on_this_step) then
        do k=1,nlev_lsm
          pmass(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/Param%GRAV   
        end do

        call donner_deep_netcdf (is, ie, js, je, Nml, Time, Don_conv,  &
                                 Don_cape, parcel_rise, pmass, &
                                 total_precip,  Don_budgets, &
                                 temperature_forcing, &
                                 moisture_forcing)

!----------------------------------------------------------------------
!    on calculation steps, update the values of the cell and
!    mesoscale cloud variables to be returned to moist_processes_mod. 
!    (on non-calculation steps, the values that were passed in are 
!    simply passed back.)
!----------------------------------------------------------------------
        cell_cld_frac = Don_rad%cell_cloud_frac
        cell_liq_amt  = Don_rad%cell_liquid_amt
        cell_liq_size = Don_rad%cell_liquid_size
        cell_ice_amt  = Don_rad%cell_ice_amt
        cell_ice_size = Don_rad%cell_ice_size
        cell_droplet_number = Don_rad%cell_droplet_number
        meso_cld_frac = Don_rad%meso_cloud_frac
        meso_liq_amt  = Don_rad%meso_liquid_amt
        meso_liq_size = Don_rad%meso_liquid_size
        meso_ice_amt  = Don_rad%meso_ice_amt
        meso_ice_size = Don_rad%meso_ice_size
        meso_droplet_number = Don_rad%meso_droplet_number
        nsum          = Don_rad%nsum

!--------------------------------------------------------------------
!    call deallocate_local_variables to deallocate space used by the
!    local derived-type variables.
!--------------------------------------------------------------------
        call don_d_dealloc_loc_vars_k   &
               (Don_conv, Don_cape, Don_rad, Don_budgets, Nml, &
                                               Initialized, ermesg)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, process the error message.
!----------------------------------------------------------------------
        if (trim(ermesg) /= ' ') then
          call error_mesg ('donner_deep_mod', ermesg, FATAL)
        endif
      endif  ! (calc_conv_on_this_step)

!--------------------------------------------------------------------


end subroutine donner_deep




!####################################################################

subroutine donner_deep_end

!---------------------------------------------------------------------
!   donner_deep_end is the destructor for donner_deep_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variable

      integer  :: ntracers     ! number of tracers transported by the
                               ! donner deep convection parameterization

!-------------------------------------------------------------------
!    if module has not been initialized, return.
!-------------------------------------------------------------------
      if (.not. module_is_initialized) return

!-------------------------------------------------------------------
!    define the number of tracers that have been transported by the 
!    donner deep convection parameterization.
!-------------------------------------------------------------------
      ntracers = size(Don_save%tracername(:))

!-------------------------------------------------------------------
!    call subroutine to write restart file. NOTE: only the netcdf 
!    restart file is currently supported.
!-------------------------------------------------------------------
      if (do_netcdf_restart) then
        call write_restart_nc (ntracers)
      else
        call write_restart (ntracers)
      endif

!-------------------------------------------------------------------
!    close any column diagnostics units which are open.
!------------------------------------------------------------------
      if (Col_diag%num_diag_pts > 0) then
        call close_column_diagnostics_units (col_diag_unit)
      endif

!----------------------------------------------------------------------
!    call deallocate_variables to deallocate the module variables.
!----------------------------------------------------------------------
      call deallocate_variables 

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.


!---------------------------------------------------------------------

end subroutine donner_deep_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                   PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      1. ROUTINES CALLED BY DONNER_DEEP_INIT
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
 
subroutine register_fields (Time, axes)

!----------------------------------------------------------------------
!    subroutine register_fields registers all of the potential diagnos-
!    tics written by this module with diag_manager_mod.
!----------------------------------------------------------------------

type(time_type),               intent(in)   :: Time
integer,         dimension(4), intent(in)   :: axes

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      Time         current time [ time_type ]
!      axes         data axes for diagnostics
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer  :: ntracers     ! number of tracers transported by the
                               ! donner deep convection parameterization
      integer :: nn            ! do-loop index

!----------------------------------------------------------------------
!    define the number of tracers that are to be transported by the 
!    donner deep convection parameterization.
!-------------------------------------------------------------------
      ntracers = size(Don_save%tracername(:))

!---------------------------------------------------------------------
!    register the various diagnostic fields.
!---------------------------------------------------------------------

    if (do_budget_analysis) then
      id_water_budget(1)    = register_diag_field    &
            (mod_name, 'vapor_net_tend', axes(1:3),   &
             Time, 'net water vapor tendency', &
             'g(h2o) / kg(air) / day',    &
             missing_value=missing_value)
      
      id_water_budget(2)    = register_diag_field    &
            (mod_name, 'vapor_cell_dynam', axes(1:3),   &
             Time, 'vapor tendency due to cell dynamics', &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(3)    = register_diag_field    &
            (mod_name, 'vapor_meso_depo', axes(1:3),   &
             Time, 'vapor tendency from mesoscale deposition', &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(4)    = register_diag_field    &
            (mod_name, 'vapor_meso_cd', axes(1:3),   &
             Time, 'vapor tendency from mesoscale condensation',  &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(5)    = register_diag_field    &
            (mod_name, 'vapor_cell_evap', axes(1:3),   &
             Time, 'vapor tendency from cell evaporation',  &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(6)    = register_diag_field    &
            (mod_name, 'vapor_cell_meso_trans', axes(1:3),   &
             Time, 'vapor tendency from cell to mesoscale transfer',  &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(7)    = register_diag_field    &
            (mod_name, 'vapor_meso_evap', axes(1:3),   &
             Time, 'vapor tendency from mesoscale evaporation', &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(8)    = register_diag_field    &
            (mod_name, 'vapor_meso_dynam_up', axes(1:3),   &
             Time, 'vapor tendency from mesoscale updrafts',  &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_water_budget(9)    = register_diag_field    &
            (mod_name, 'vapor_meso_dynam_dn',  axes(1:3),   &
             Time, 'vapor tendency from mesoscale downdrafts',  &
             ' g(h2o) / kg(air) / day', &
             missing_value=missing_value)
      
      id_enthalpy_budget(1)    = register_diag_field    &
            (mod_name, 'enth_net_tend', axes(1:3),   &
             Time, 'net temp tendency', 'deg K  /day',    &
             missing_value=missing_value)

      id_enthalpy_budget(2)    = register_diag_field    &
            (mod_name, 'enth_cell_dynam', axes(1:3),   &
             Time, 'temp tendency due to cell dynamics', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(3)    = register_diag_field    &
            (mod_name, 'enth_meso_depo_liq', axes(1:3), Time, &
             'temp tendency from mesoscale deposition on liquid&
                    & condensate',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(4)    = register_diag_field    &
            (mod_name, 'enth_meso_cd_liq', axes(1:3), Time, &
             ' temp tendency from mesoscale liquid condensation', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(5)    = register_diag_field    &
            (mod_name, 'enth_cell_evap_liq', axes(1:3),   &
             Time, 'temp tendency from evap of liquid condensate', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(6)    = register_diag_field    &
            (mod_name, 'enth_meso_evap_liq_up', axes(1:3),   &
             Time, 'temp tendency from evaporation of liquid &
              &condensate in mesoscale updrafts',  &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(7)    = register_diag_field    &
            (mod_name, 'enth_meso_evap_liq_dn', axes(1:3),   &
             Time, 'temp tendency from evaporation of liquid &
              &condensate in mesoscale downdrafts',  &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(8)    = register_diag_field    &
            (mod_name, 'enth_meso_depo_ice', axes(1:3),   &
             Time, ' temp tendency from mesoscale deposition on &
              &ice condensate',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(9)    = register_diag_field    &
            (mod_name, 'enth_meso_cd_ice', axes(1:3),   &
             Time, 'temp tendency from mesoscale ice condensation', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(10)    = register_diag_field    &
            (mod_name, 'enth_cell_evap_ice', axes(1:3),   &
             Time, 'temp tendency from evap of ice condensate', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(11)    = register_diag_field    &
            (mod_name, 'enth_meso_evap_ice_up', axes(1:3),   &
             Time, 'temp tendency from evaporation of ice condensate &
              &in mesoscale updrafts',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(12)    = register_diag_field    &
            (mod_name, 'enth_meso_evap_ice_dn', axes(1:3),   &
             Time, 'temp tendency from evaporation of ice &
               &condensate in mesoscale downdrafts',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(13)    = register_diag_field    &
            (mod_name, 'enth_meso_freeze', axes(1:3),   &
             Time, 'temp tendency from the freezing of liquid &
              &condensate when it enters the mesoscale circulation',  &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(14)    = register_diag_field    &
            (mod_name, 'enth_cell_freeze', axes(1:3),   &
             Time, 'temp tendency from the freezing of liquid &
                &cell condensate',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(15)    = register_diag_field    &
            (mod_name, 'enth_cell_precip_melt', axes(1:3),   &
             Time, 'temp tendency from the melting of cell frozen &
             liquid and ice that is precipitating out', 'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(16)    = register_diag_field    &
            (mod_name, 'enth_meso_melt', axes(1:3), Time, &
             'temp tendency from melting bogus frozen condensate',  &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(17)    = register_diag_field    &
            (mod_name, 'enth_meso_precip_melt', axes(1:3),   &
             Time, 'temp tendency from the melting of frozen &
               &mesoscale precipitation',  'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(18)    = register_diag_field    &
            (mod_name, 'enth_meso_dynam_up', axes(1:3),   &
             Time, 'temp tendency from mesoscale updraft', &
             'deg K / day', &
             missing_value=missing_value)

      id_enthalpy_budget(19)    = register_diag_field    &
            (mod_name, 'enth_meso_dynam_dn', axes(1:3),   &
             Time, 'temp tendency from mesoscale downdraft', &
             'deg K / day', &
             missing_value=missing_value)

      id_precip_budget(1,1)    = register_diag_field    &
            (mod_name, 'precip_cell_liq', axes(1:3),   &
             Time, 'precip from cell liquid condensate', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(2,1)    = register_diag_field    &
            (mod_name, 'precip_cell_liq_frz', axes(1:3),   &
             Time, 'precip from cell liquid condensate which froze', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(3,1)    = register_diag_field    &
            (mod_name, 'precip_cell_liq_frz_melt', axes(1:3), Time, &
              'precip from cell liquid condensate which froze &
               &and remelted', 'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(4,1)    = register_diag_field    &
            (mod_name, 'precip_cell_ice', axes(1:3),   &
             Time, 'precip from cell ice condensate', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(5,1)    = register_diag_field    &
            (mod_name, 'precip_cell_ice_melt', axes(1:3),   &
             Time, 'precip from cell ice condensate which melted', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(1,2)    = register_diag_field    &
            (mod_name, 'precip_trans_liq', axes(1:3),   &
             Time, 'precip from cell liquid transferred to meso', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(2,2)    = register_diag_field    &
            (mod_name, 'precip_trans_liq_frz', axes(1:3),   &
             Time, 'precip from cell liquid transferred to meso &
              which froze', 'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(3,2)    = register_diag_field    &
            (mod_name, 'precip_trans_liq_frz_melt', axes(1:3), Time, &
             'precip from cell liquid transferred to meso which &
              &froze and remelted', 'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(4,2)    = register_diag_field    &
            (mod_name, 'precip_trans_ice', axes(1:3),   &
             Time, 'precip from cell ice transferred to meso', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(5,2)    = register_diag_field    &
            (mod_name, 'precip_trans_ice_melt', axes(1:3),   &
             Time, 'precip from cell ice transferred to meso &
              &which melted', 'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(1,3)    = register_diag_field    &
            (mod_name, 'precip_meso_liq', axes(1:3),   &
             Time, 'precip from meso liq condensate', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(2,3)    = register_diag_field    &
            (mod_name, 'precip_meso_liq_frz', axes(1:3),   &
             Time, 'precip from meso liq  condensate which froze', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(3,3)    = register_diag_field    &
            (mod_name, 'precip_meso_liq_frz_melt', axes(1:3), Time, &
            'precip from meso condensate liq which froze and &
             &remelted', 'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(4,3)    = register_diag_field    &
            (mod_name, 'precip_meso_ice', axes(1:3),   &
             Time, 'precip from meso ice condensate', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_precip_budget(5,3)    = register_diag_field    &
            (mod_name, 'precip_meso_ice_melt', axes(1:3),   &
             Time, 'precip from meso ice condensate which melted', &
              'kg(h2o) / kg(air) / day', &
             missing_value=missing_value)

      id_ci_precip_budget(1,1)    = register_diag_field    &
            (mod_name, 'ci_precip_cell_liq', axes(1:2),   &
             Time, 'col intg precip from cell liquid condensate', &
             'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(2,1)    = register_diag_field    &
            (mod_name, 'ci_precip_cell_liq_frz', axes(1:2),   &
             Time, 'col intg precip from cell liquid condensate &
             which froze',  'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(3,1)    = register_diag_field    &
            (mod_name, 'ci_precip_cell_liq_frz_melt', axes(1:2), Time, &
             'col intg precip from cell liquid condensate which &
              &froze and remelted',  'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(4,1)    = register_diag_field    &
            (mod_name, 'ci_precip_cell_ice', axes(1:2),   &
             Time, 'col intg precip from cell ice condensate', &
             'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(5,1)    = register_diag_field    &
            (mod_name, 'ci_precip_cell_ice_melt', axes(1:2),   &
             Time, 'col intg precip from cell ice condensate &
             &which melted',  'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(1,2)    = register_diag_field    &
            (mod_name, 'ci_precip_trans_liq', axes(1:2),   &
             Time, 'col intg precip from cell liquid transferred &
             &to meso',  'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(2,2)    = register_diag_field    &
            (mod_name, 'ci_precip_trans_liq_frz', axes(1:2),   &
             Time, 'col intg precip from cell liquid transferred &
              to meso  which froze', 'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(3,2)    = register_diag_field    &
            (mod_name, 'ci_precip_trans_liq_frz_melt', axes(1:2), &
             Time, 'col intg precip from cell liquid transferred &
             &to meso which froze and remelted', 'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(4,2)    = register_diag_field    &
            (mod_name, 'ci_precip_trans_ice', axes(1:2),   &
             Time, 'col intg precip from cell ice transferred &
              &to meso', 'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(5,2)    = register_diag_field    &
            (mod_name, 'ci_precip_trans_ice_melt', axes(1:2),   &
             Time, 'col intg precip from cell ice transferred to &
             &meso which melted', 'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(1,3)    = register_diag_field    &
            (mod_name, 'ci_precip_meso_liq', axes(1:2),   &
             Time, 'col intg precip from meso liq condensate', &
             'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(2,3)    = register_diag_field    &
            (mod_name, 'ci_precip_meso_liq_frz', axes(1:2),   &
             Time, 'col intg precip from meso liq  condensate &
             &which froze',  'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(3,3)    = register_diag_field    &
            (mod_name, 'ci_precip_meso_liq_frz_melt', axes(1:2), Time, &
             'col intg precip from meso condensate liq which froze &
               &and remelted', 'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(4,3)    = register_diag_field    &
            (mod_name, 'ci_precip_meso_ice', axes(1:2),   &
             Time, 'col intg precip from meso ice condensate', &
             'mm / day', &
             missing_value=missing_value)

      id_ci_precip_budget(5,3)    = register_diag_field    &
            (mod_name, 'ci_precip_meso_ice_melt', axes(1:2),   &
             Time, 'col intg precip from meso ice condensate &
              &which melted', 'mm / day', &
             missing_value=missing_value)

      id_ci_water_budget(1)    = register_diag_field    &
            (mod_name, 'ci_vapor_net_tend', axes(1:2),   &
             Time, 'col intg net water vapor tendency', 'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(2)    = register_diag_field    &
            (mod_name, 'ci_vapor_cell_dynam', axes(1:2),   &
             Time, 'col intg vapor tendency due to cell dynamics', &
              'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(3)    = register_diag_field    &
            (mod_name, 'ci_vapor_meso_depo', axes(1:2),   &
             Time, 'col intg vapor tendency from mesoscale deposition',&
              'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(4)    = register_diag_field    &
            (mod_name, 'ci_vapor_meso_cd', axes(1:2),   &
             Time, 'col intg vapor tendency from mesoscale &
              &condensation',  'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(5)    = register_diag_field    &
            (mod_name, 'ci_vapor_cell_evap', axes(1:2),   &
             Time, 'col intg vapor tendency from cell evaporation', &
              'mm / day', missing_value=missing_value)
      
      id_ci_water_budget(6)    = register_diag_field    &
            (mod_name, 'ci_vapor_cell_meso_trans', axes(1:2),   &
             Time, 'col intg vapor tendency from cell to mesoscale &
              &transfer',  'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(7)    = register_diag_field    &
            (mod_name, 'ci_vapor_meso_evap', axes(1:2),   &
             Time, 'col intg vapor tendency from mesoscale &
              &evaporation', 'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(8)    = register_diag_field    &
            (mod_name, 'ci_vapor_meso_dynam_up', axes(1:2),   &
             Time, 'col intg vapor tendency from mesoscale updrafts',  &
              'mm / day',    &
             missing_value=missing_value)
      
      id_ci_water_budget(9)    = register_diag_field    &
            (mod_name, 'ci_vapor_meso_dynam_dn',  axes(1:2),   &
             Time, 'col intg vapor tendency from mesoscale downdrafts',&
              'mm / day',    &
             missing_value=missing_value)
      
      id_ci_enthalpy_budget(1)    = register_diag_field    &
            (mod_name, 'ci_enth_net_tend', axes(1:2),   &
             Time, 'col intg net enthalpy tendency', 'J/m**2 / day',   &
             missing_value=missing_value)

      id_ci_enthalpy_budget(2)    = register_diag_field    &
            (mod_name, 'ci_enth_cell_dynam', axes(1:2),   &
             Time, 'col intg enthalpy tendency due to cell dynamics', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(3)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_depo_liq', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale &
             deposition on liquid condensate',  'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(4)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_cd_liq', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale &
             liquid condensation', 'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(5)    = register_diag_field    &
            (mod_name, 'ci_enth_cell_evap_liq', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evap of liquid &
             &condensate',  'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(6)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_evap_liq_up', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evaporation of &
             &liquid condensate in mesoscale updrafts',  &
             'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(7)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_evap_liq_dn', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evaporation &
             &of liquid condensate in mesoscale downdrafts',  &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(8)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_depo_ice', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale &
              &deposition on ice condensate',  &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(9)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_cd_ice', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale ice &
             condensation', 'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(10)    = register_diag_field    &
            (mod_name, 'ci_enth_cell_evap_ice', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evap of ice &
              &condensate', 'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(11)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_evap_ice_up', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evaporation of &
             &ice condensate in mesoscale updrafts',  'J/m**2 / day',  &
             missing_value=missing_value)

      id_ci_enthalpy_budget(12)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_evap_ice_dn', axes(1:2),   &
             Time, 'col intg enthalpy tendency from evaporation of &
             &ice condensate in mesoscale downdrafts',  &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(13)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_freeze', axes(1:2),   &
             Time, 'col intg enthalpy tendency from the freezing of &
             &liquid condensate when it enters the mesoscale &
             &circulation',  'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(14)    = register_diag_field    &
            (mod_name, 'ci_enth_cell_freeze', axes(1:2),   &
             Time, 'col intg enthalpy tendency from the freezing of &
             liquid cell condensate', 'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(15)    = register_diag_field    &
            (mod_name, 'ci_enth_cell_precip_melt', axes(1:2),   &
             Time, 'col intg enthalpy tendency from the melting of &
             &cell frozen liquid and ice that is precipitating out',  &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(16)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_melt', axes(1:2),   &
             Time, 'col intg enthalpy tendency from melting bogus &
              &frozen condensate',  'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(17)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_precip_melt', axes(1:2),   &
             Time, 'col intg enthalpy tendency from the melting of &
              frozen mesoscale precipitation',  &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(18)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_dynam_up', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale updraft',&
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_enthalpy_budget(19)    = register_diag_field    &
            (mod_name, 'ci_enth_meso_dynam_dn', axes(1:2),   &
             Time, 'col intg enthalpy tendency from mesoscale &
             &downdraft',  'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_heat_frz_cell =  register_diag_field & 
            (mod_name, 'ci_prcp_heat_frz_cell', axes(1:2),   &
             Time, 'col intg heat removed by frozen cell precip', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_heat_liq_cell =  register_diag_field & 
            (mod_name, 'ci_prcp_heat_liq_cell', axes(1:2),   &
             Time, 'col intg heat removed by liquid cell precip', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_heat_frz_meso =  register_diag_field & 
            (mod_name, 'ci_prcp_heat_frz_meso', axes(1:2),   &
             Time, 'col intg heat removed by frozen meso precip', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_heat_liq_meso =  register_diag_field & 
            (mod_name, 'ci_prcp_heat_liq_meso', axes(1:2),   &
             Time, 'col intg heat removed by liquid meso precip', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_heat_total =  register_diag_field & 
            (mod_name, 'ci_prcp_heat_total', axes(1:2),   &
             Time, 'col intg total heat removed by precip', &
              'J/m**2 / day',    &
             missing_value=missing_value)

      id_ci_prcp_total =  register_diag_field & 
            (mod_name, 'ci_prcp_total', axes(1:2),   &
             Time, 'col intg total precip', &
              'mm / day',    &
             missing_value=missing_value)

    endif

      id_leff          = register_diag_field    &
            (mod_name, 'leff_don', axes(1:2),   &
             Time, 'effective latent heat with donner precip ',  &
             'J/kg(h2o)',  missing_value=missing_value)

!    heating rate:
      id_cemetf_deep = register_diag_field    &
            (mod_name, 'cemetf_deep', axes(1:3),   &
             Time, 'heating rate, c + m ', 'K/s',   &
             missing_value=missing_value)

!    cell entropy flux convergence:
      id_ceefc_deep = register_diag_field   &
            (mod_name, 'ceefc_deep', axes(1:3),   &
             Time, 'cell entrpy flx cnvrgnc', 'K/s',   &
             missing_value=missing_value)

!    cell condensation / evaporation:
      id_cecon_deep = register_diag_field      &
            (mod_name, 'cecon_deep', axes(1:3),   &
             Time, 'cell cond/evap ', 'K/s',   &
             missing_value=missing_value)

!    cell moisture flux convergence:
      id_cemfc_deep = register_diag_field       &
            (mod_name, 'cemfc_deep', axes(1:3),   &
             Time, 'cell moist flx cnvgnc', 'kg(h2o)/kg/s',   &
             missing_value=missing_value)

!    moistening rate:
      id_cememf_deep = register_diag_field        &
            (mod_name, 'cememf_deep', axes(1:3),   &
             Time, 'moistening rate, c + m ', 'kg(h2o)/kg/s',   &
             missing_value=missing_value)

!    moistening rate after adjustment for negative vapor mixing ratio:
      id_cememf_mod_deep = register_diag_field       &
            (mod_name, 'cememf_mod_deep', axes(1:3),&
             Time, 'mod cememf due to negative q ', 'kg(h2o)/kg/s',   &
             missing_value=missing_value)

!    cell + mesoscale cloud fraction:
      id_cual_deep = register_diag_field      &
            (mod_name, 'cual_deep', axes(1:3),   &
             Time, 'c + m cld frac ', 'percent',   &
             missing_value=missing_value)

!    heating rate due to freezing:
      id_fre_deep = register_diag_field     &
            (mod_name, 'fre_deep', axes(1:3),   &
             Time, 'freezing ', 'K/sec',   &
             missing_value=missing_value)

!    heating rate due to melting:
      id_elt_deep = register_diag_field         &
            (mod_name, 'elt_deep', axes(1:3),   &
             Time, 'melting', 'K/sec',   &
             missing_value=missing_value)

!    deposition in mesoscale updraft:
      id_cmus_deep = register_diag_field        &
            (mod_name, 'cmus_deep', axes(1:3),   &
             Time, 'meso-up deposition', 'kg(h2o)/kg/sec)',   &
             missing_value=missing_value)

!    evaporation in convective downdraft:
      id_ecds_deep = register_diag_field    &
            (mod_name, 'ecds_deep', axes(1:3),   &
             Time, 'convective dwndrft evap ', 'kg(h2o)/kg/sec', &
             missing_value=missing_value)

!    evaporation / sublimation in convective updraft:
      id_eces_deep = register_diag_field       &
            (mod_name, 'eces_deep', axes(1:3),   &
             Time, 'convective updrft evap/subl ', 'kg(h2o)/kg/sec',  &
             missing_value=missing_value)

!    sublimation in mesoscale downdraft:
      id_emds_deep = register_diag_field     &
            (mod_name, 'emds_deep', axes(1:3),   &
             Time, 'meso-dwn subl ', 'kg(h2o)/kg/sec',   &
             missing_value=missing_value)

!    sublimation in mesoscale updraft:
      id_emes_deep = register_diag_field        &
            (mod_name, 'emes_deep', axes(1:3),   &
             Time, 'meso-up subl ', 'kg(h2o)/kg/sec',   &
             missing_value=missing_value)

!    mesoscale moisture flux convergence:
      id_qmes_deep = register_diag_field     &
             (mod_name, 'qmes_deep', axes(1:3),   &
              Time, 'meso moist flux conv', 'kg(h2o)/kg/sec',   &
              missing_value=missing_value)

!    transfer of vapor from cells to mesoscale:
      id_wmps_deep = register_diag_field      &
             (mod_name, 'wmps_deep', axes(1:3),   &
              Time, 'meso redistrib of vapor from cells',  &
              'kg(h2o)/kg/sec', missing_value=missing_value)

!    deposition of vapor from cells to mesoscale:
      id_wmms_deep = register_diag_field         &
             (mod_name, 'wmms_deep', axes(1:3),   &
              Time, 'meso depo of vapor from cells',    &
              'kg(h2o)/kg/sec',  missing_value=missing_value)

!    mesoscale entropy flux convergesnce:
      id_tmes_deep = register_diag_field         &
            (mod_name, 'tmes_deep', axes(1:3),   &
             Time, 'meso entropy flux conv',  'K/sec',   &
              missing_value=missing_value)
 
!    mass flux in mesoscale downdrafts:    
      id_dmeml_deep = register_diag_field      &
            (mod_name, 'dmeml_deep', axes(1:3), &  
             Time, 'mass flux meso dwndrfts', 'kg/((m**2) s)',   &
             missing_value=missing_value)

!    mass flux in cell updrafts:
      id_uceml_deep = register_diag_field     &
            (mod_name, 'uceml_deep', axes(1:3), &
             Time, 'mass flux cell updrfts', 'kg/((m**2) s)',   &
             missing_value=missing_value)

!    mass flux in mesoscale updrafts:
      id_umeml_deep = register_diag_field       &
            (mod_name, 'umeml_deep', axes(1:3), &
             Time, 'mass flux meso updrfts', 'kg/((m**2) s)',   &
             missing_value=missing_value)

!    mesoscale ice mass mixing ratio:
      id_xice_deep = register_diag_field     &
            (mod_name, 'xice_deep', axes(1:3),  &
             Time, 'meso ice mass mixing ratio ', 'kg(ice)/kg',   &
             missing_value=missing_value)

!    mesoscale liquid mass mixing ratio:
      id_xliq_deep = register_diag_field       &
            (mod_name, 'xliq_deep', axes(1:3),  &
             Time, 'meso liq mass mixing ratio ', 'kg(liq)/kg',   &
             missing_value=missing_value)

!    detrained mass flux:
      id_detmfl_deep = register_diag_field       &
            (mod_name, 'detmfl_deep', axes(1:3),  &
             Time, 'detrained mass flux ', 'kg/((m**2) s)',   &
             missing_value=missing_value)

!---------------------------------------------------------------------
!    if tracers are being transported by donner_deep_mod, allocate diag-
!    nostic indices for each tracer and register their diagnostics.
!---------------------------------------------------------------------
      if (ntracers > 0) then
        allocate (id_qtren1 (ntracers))
        allocate (id_qtmes1 (ntracers))
        allocate (id_wtp1   (ntracers))
        allocate (id_qtceme (ntracers))
        allocate (id_total_wet_dep (ntracers))
        allocate (id_meso_wet_dep (ntracers))
        allocate (id_cell_wet_dep (ntracers))
        allocate (id_qtren1_col (ntracers))
        allocate (id_qtmes1_col (ntracers))
        allocate (id_wtp1_col   (ntracers))
        allocate (id_qtceme_col (ntracers))
        allocate (id_total_wet_dep_col (ntracers))
        allocate (id_meso_wet_dep_col (ntracers))
        allocate (id_cell_wet_dep_col (ntracers))
        do nn=1,ntracers

!    tracer tendency due to cells:
          id_qtren1(nn) = register_diag_field     &
                (mod_name, trim(Don_save%tracername(nn)) // '_qtren1',  &
                 axes(1:3), Time,  &
                 trim(Don_save%tracername(nn)) // ' cell tendency ', &
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    tracer tendency due to mesoscale circulation:
          id_qtmes1(nn) = register_diag_field    &
                (mod_name, trim(Don_save%tracername(nn)) // '_qtmes1', &
                 axes(1:3), Time,   &
                 trim(Don_save%tracername(nn)) //' mesoscale tendency',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    tracer tendency due to mesoscale redistribution:
          id_wtp1(nn) = register_diag_field         &
                (mod_name, trim(Don_save%tracername(nn)) // '_wtp1',  &
                 axes(1:3), Time,  &
                 trim(Don_save%tracername(nn)) //' mesoscale redist',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

 !    tracer tendency due to deep convective wet deposition:
          id_total_wet_dep(nn) = register_diag_field         &
              (mod_name, trim(Don_save%tracername(nn)) // '_totwdep',  &
                 axes(1:3), Time,  &
                 trim(Don_save%tracername(nn)) //' deep conv wet depo',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    tracer tendency due to wet deposition in mesoscale updrafts:
         id_meso_wet_dep(nn) = register_diag_field         &
                 (mod_name, trim(Don_save%tracername(nn)) // '_mwdep', &
                  axes(1:3), Time,   &
                 trim(Don_save%tracername(nn)) //' mesoscale wet depo',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    tracer tendency due to wet deposition in cells:
          id_cell_wet_dep(nn) = register_diag_field         &
                (mod_name, trim(Don_save%tracername(nn)) // '_cwdep', &
                 axes(1:3), Time,  &
                 trim(Don_save%tracername(nn)) //' cell wet depo',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    total tracer tendency:
          id_qtceme(nn) = register_diag_field     &
                (mod_name, trim(Don_save%tracername(nn)) // '_qtceme', &
                 axes(1:3), Time,  &
                 trim(Don_save%tracername(nn)) // ' total tendency ',&
                 trim(Don_save%tracer_units(nn))//'/s', &
                 missing_value=missing_value)

!    column-integrated tracer tendency due to cells:
          id_qtren1_col(nn) = register_diag_field      &
                (mod_name,       &
                 trim(Don_save%tracername(nn)) // '_qtren1_col',  &
                 axes(1:2), Time,  & 
                 'column integrated ' //trim(Don_save%tracername(nn)) //&
                 ' cell tendency ', &
                 trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                 missing_value=missing_value)

!    column-integrated tracer tendency due to mesoscale circulation:
          id_qtmes1_col(nn) = register_diag_field    &
                (mod_name,          &
                 trim(Don_save%tracername(nn)) // '_qtmes1_col',  &
                 axes(1:2), Time,   &
                 'column integrated ' //trim(Don_save%tracername(nn)) //&
                 ' mesoscale tendency',&
                trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                 missing_value=missing_value)

!    column-integrated tracer tendency due to mesoscale redistribution:
          id_wtp1_col(nn) = register_diag_field     &
                (mod_name,  &
                 trim(Don_save%tracername(nn)) // '_wtp1_col',   &
                 axes(1:2), Time,  &
                 'column integrated '//trim(Don_save%tracername(nn)) // &
                 ' mesoscale redist',&
               trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                missing_value=missing_value)

!    column-integrated tracer tendency due to deep convective wet 
!    deposition: 
          id_total_wet_dep_col(nn) = register_diag_field     &
                 (mod_name,  &
                  trim(Don_save%tracername(nn)) // '_totwdep_col',   &
                  axes(1:2), Time,  &
                'column integrated '//trim(Don_save%tracername(nn)) // &
                ' deep convective wet depo',&
                 trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                 missing_value=missing_value)

!    column-integrated tracer tendency due to mesocscale updraft  wet 
!    deposition: 
          id_meso_wet_dep_col(nn) = register_diag_field     &
                (mod_name,  &
                  trim(Don_save%tracername(nn)) // '_mwdep_col',   &
                  axes(1:2), Time,  &
                'column integrated '//trim(Don_save%tracername(nn)) // &
                 ' meso updraft wet depo',&
                  trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                 missing_value=missing_value)

!    column-integrated tracer tendency due to wet deposition in cells:
          id_cell_wet_dep_col(nn) = register_diag_field     &
                (mod_name,  &
                 trim(Don_save%tracername(nn)) // '_cwdep_col',   &
                  axes(1:2), Time,  &
                'column integrated '//trim(Don_save%tracername(nn)) // &
                  ' cell wet depo',&
                 trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                  missing_value=missing_value)

!    column-integrated total tracer tendency:
          id_qtceme_col(nn) = register_diag_field     &
                (mod_name,  &
                 trim(Don_save%tracername(nn)) // '_qtceme_col',  &
                 axes(1:2), Time,  &
                 'column integrated ' //trim(Don_save%tracername(nn)) //&
                 ' total tendency ', &
                  trim(Don_save%tracer_units(nn)) // '* kg/(m**2 s) ', &
                 missing_value=missing_value)
        end do
      endif

!    mesoscale ice generalized effective size:
      id_dgeice_deep = register_diag_field    &
            (mod_name, 'dgeice_deep', axes(1:3), &
             Time, 'meso ice gen eff size ', 'micrometers',   &
             missing_value=missing_value)

!    cell ice mixing ratio:
      id_cuqi_deep = register_diag_field         &
            (mod_name, 'cuqi_deep', axes(1:3),  &
             Time, 'cell ice ', 'kg(H2O)/kg',   &
             missing_value=missing_value)

!    cell liquid mixing ratio:
      id_cuql_deep = register_diag_field     &
            (mod_name, 'cuql_deep', axes(1:3),  &
             Time, 'cell liquid ', 'kg(H2O)/kg',   &
             missing_value=missing_value)

!    cell liquid generalized effective size:
      id_dgeliq_deep = register_diag_field    &
            (mod_name, 'dgeliq_deep', axes(1:3), &
             Time, 'cell liq gen eff size ', 'micrometers',   &
             missing_value=missing_value)

!    pressure at lifting condensation level:
      id_plcl_deep = register_diag_field       &
            (mod_name, 'plcl_deep', axes(1:2),   &
             Time, 'pressure at lcl ', 'Pa ',   &
             missing_value=missing_value)

!    pressure at level of free convection:
      id_plfc_deep = register_diag_field     &
            (mod_name, 'plfc_deep', axes(1:2),   &
             Time, 'pressure at lfc ', 'Pa ',   &
             missing_value=missing_value)

!    pressure at level of zero buoyancy:  
      id_plzb_deep = register_diag_field      &
            (mod_name, 'plzb_deep', axes(1:2),   &
             Time, 'pressure at lzb ', 'Pa ',   &
             missing_value=missing_value)

!    convective available potential energy (cape):
      id_xcape_deep = register_diag_field      &
            (mod_name, 'xcape_deep', axes(1:2),  &
             Time, 'cape', 'J/kg',   &
             missing_value=missing_value)

!    convective inhibition:
      id_coin_deep = register_diag_field      &
            (mod_name, 'coin_deep', axes(1:2),   &
             Time, 'convective inhibition ', 'J/kg',   &
             missing_value=missing_value)

!    time tendency of cape:
      id_dcape_deep = register_diag_field      &
            (mod_name, 'dcape_deep', axes(1:2), &
             Time, 'time tendency of cape ', 'J/kg/sec',   &
             missing_value=missing_value)

!    column integrated water vapor:
      id_qint_deep = register_diag_field    &
            (mod_name, 'qint_deep', axes(1:2),   &
             Time, 'column moisture ', 'kg(h2o)/m**2',   &
             missing_value=missing_value)

!    fractional area of cumulus ensemble member:
      id_a1_deep = register_diag_field           &
            (mod_name, 'a1_deep', axes(1:2),   &
             Time, 'fractional area of cu subensemble ', 'percent',   &
             missing_value=missing_value)

!    fractional area of largest cumulus ensemble member:
      id_amax_deep = register_diag_field      &
            (mod_name, 'amax_deep', axes(1:2),   &
             Time, 'fractional area of largest cu subensemble ',  &
             'percent',  missing_value=missing_value)

!    upper limit onfractional area based on moisture constraint:
      id_amos_deep = register_diag_field      &
            (mod_name, 'amos_deep', axes(1:2),   &
             Time, 'uppr lmt on frac area from moisture', 'percent', &
             missing_value=missing_value)

!    area-weighted total precipitation:
      id_tprea1_deep = register_diag_field         &
            (mod_name, 'tprea1_deep', axes(1:2), &
             Time, 'area wtd total precip ', 'mm/day',   &
             missing_value=missing_value)

!    mesoscale cloud fraction:
      id_ampta1_deep = register_diag_field       &
            (mod_name, 'ampta1_deep', axes(1:2), &
             Time, 'meso cld frac', 'percent',   &
             missing_value=missing_value)

!    accumulated low-level vertical displacement:
      id_omint_deep = register_diag_field      &
            (mod_name, 'omint_deep', axes(1:2), &
             Time, 'accumulated low-lvl displ', 'Pa ',   &
             missing_value=missing_value)

!    area-weighted convective precipitation:
      id_rcoa1_deep = register_diag_field     &
            (mod_name, 'rcoa1_deep', axes(1:2),  &
             Time, 'area wtd cnvctv precip ', 'mm/day',   &
             missing_value=missing_value)

!----------------------------------------------------------------------


end subroutine register_fields 



!####################################################################

subroutine read_restart (ntracers, Time)

!---------------------------------------------------------------------
!    subroutine read_restart reads a native mode restart file, which are
!    not written by this code version. currently only restart version #8 
!    may be read to provide initial conditions for an experiment run with
!    this code version. this routine remains as a template for any user 
!    who is unable to process the current standard netcdf restart file, 
!    and must modify the current code to write a native mode file. 
!---------------------------------------------------------------------

integer, intent(in)         :: ntracers
type(time_type), intent(in) :: Time

!----------------------------------------------------------------------
!   intent(in) variables:
!
!     ntracers               number of tracers to be transported by
!                            the donner deep convection parameterization
!     Time                   current time [ time_type ]
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!   local variables:

      logical, dimension(ntracers)  :: success
      integer                       :: old_freq
      integer                       :: unit, vers
      character(len=8)              :: chvers
      character(len=32)             :: tracername_in
      integer                       :: ntracers_in
      integer                       :: n, nn, k

!-----------------------------------------------------------------------
!   local variables:
!
!     success      logical array indicating whether data for each trans-
!                  ported tracer is present in restart file
!     old_freq     donner_Deep_freq used in job which wrote the restart
!                  file, used in versions 5 and higher [ seconds ]
!     unit         io unit number assigned to restart file
!     vers         restart version number of file being read
!     chvers       character representation of restart version of file
!                  being read
!     tracername_in
!                  tracer name read from restart file, used in versions
!                  6, 7 and 8
!     ntracers_in  number of tracers contained in restart file, used in
!                  versions 6, 7 and 8.
!     n, nn, k     do-loop indices
!
!--------------------------------------------------------------------


!-------------------------------------------------------------------- 
!    open the restart file.
!--------------------------------------------------------------------- 
      unit = open_restart_file ('INPUT/donner_deep.res', 'read')

!--------------------------------------------------------------------- 
!    read and check restart version number. 
!-------------------------------------------------------------------- 
      read (unit) vers 
      if ( .not. any(vers == restart_versions) ) then 
        write (chvers,'(i4)') vers 
        call error_mesg ('donner_deep_mod', 'read_restart: &  
            &restart version '//chvers//' cannot be used'//& 
            'as a restart file for the current code release; &
            & a COLDSTART will be initiated', NOTE)
         call process_coldstart (Time)
         return
      endif 
      if (vers >= 9) then 
        call error_mesg ('donner_deep_mod', 'read_restart: & 
         &native mode restart versions above #8 are totally the &
         &responsibility of the user; be sure you process it properly!',&
                                                                   NOTE) 
      endif 

!-------------------------------------------------------------------- 
!    read the time remaining before the next calculation call ( which
!    becomes Initialized%conv_alarm, in seconds) and the donner deep 
!    frequency used in the job writing the file, also in seconds 
!    (old_freq).
!---------------------------------------------------------------------
      read (unit) Initialized%conv_alarm, old_freq

!--------------------------------------------------------------------
!    determine if it is desired to change the donner_deep_freq from that
!    used in the previous job. if so, modify the alarm as read from the 
!    restart file.
!--------------------------------------------------------------------
      if (donner_deep_freq /= old_freq ) then
        Initialized%conv_alarm = Initialized%conv_alarm - old_freq + &
                                 donner_deep_freq
        if (mpp_pe() == mpp_root_pe()) then
          call error_mesg ('donner_deep_mod', 'read_restart:  &
            &donner_deep time step has changed', NOTE)
        endif
      endif

!---------------------------------------------------------------------
!    read the total heating and moistening rates produced by the donner
!    deep convection parameterization from the restart file.
!---------------------------------------------------------------------
      call read_data (unit, Don_save%cemetf)
      call read_data (unit, Don_save%cememf)

!----------------------------------------------------------------------
!    read the mass flux and large-scale cloud tendencies needed by 
!    strat_cloud_mod. if this is an earlier file, set these values to 
!    0.0.
!----------------------------------------------------------------------
      call read_data (unit, Don_save%mass_flux)
      call read_data (unit, Don_save%dql_strat )
      call read_data (unit, Don_save%dqi_strat )
      call read_data (unit, Don_save%dqa_strat )

!----------------------------------------------------------------------
!    read the accumulated vertical displacement of a boundary layer 
!    parcel.
!----------------------------------------------------------------------
      call read_data (unit, Don_save%parcel_disp)

!----------------------------------------------------------------------
!    read the total precipitation produced by the donner parameteriz-
!    ation.
!----------------------------------------------------------------------
      call read_data (unit, Don_save%tprea1)

!----------------------------------------------------------------------
!    read the temperature, mixing ratio and pressure fields at the lag 
!    time step from the restart file.        
!----------------------------------------------------------------------
      call read_data (unit, Don_save%lag_temp)
      call read_data (unit, Don_save%lag_vapor)
      call read_data (unit, Don_save%lag_press)

!----------------------------------------------------------------------
!    two fields which are needed by strat_cloud_mod are available and 
!    are read in. 
!----------------------------------------------------------------------
      call read_data (unit, Don_save%humidity_area)
      if (vers == 9) then
        call error_mesg ('donner_deep_mod', &
          'version 9 not acceptable restart -- needs to have humidity_factor&
            & rather than humidity_ratio', FATAL)
      else
        call read_data (unit, Don_save%humidity_factor)
      endif

!------------------------------------------------------------------
!    if tracers are to be transported by the donner parameterization,
!    determine if the current tendencies are available on the restart.
!------------------------------------------------------------------
      if (Initialized%do_donner_tracer) then

!------------------------------------------------------------------
!    read the number of tracers whose tendencies are included in 
!    this file. tracer tendencies are available only in version #6 and
!    higher.
!-------------------------------------------------------------------
        success = .false.
        read (unit) ntracers_in 

!--------------------------------------------------------------------
!    read each restart file tracer's name and see if it is to be 
!    transported in the current job.
!--------------------------------------------------------------------
        do n=1,ntracers_in
          read (unit) tracername_in
          do nn=1,ntracers

!--------------------------------------------------------------------
!    if the tracer is needed in the current job, read its data and
!    store it in the appropriate array. write a note indicating that 
!    the data has bben found and set a logical variable to also 
!    indicate such. exit this loop and process the next tracer present
!    in the restart file.
!--------------------------------------------------------------------
            if (trim(tracername_in) ==     &
                trim(Don_save%tracername(nn))) then
              call read_data(unit, Don_save%tracer_tends(:,:,:,nn))
              if (mpp_pe() == mpp_root_pe() ) then
                call error_mesg ('donner_deep_mod', 'read_restart: &
                         &found tracer restart data for ' // &
                         trim(Don_save%tracername(nn)), NOTE)
              endif
              success(nn) = .true.
              exit 

!---------------------------------------------------------------------
!    if the tracer in the restart file is not needed by the current
!    job, do a dummy read to get to the next record.
!---------------------------------------------------------------------
            else 
              if (nn == ntracers) then
                read (unit)
              endif
            endif
          end do
        end do

!---------------------------------------------------------------------
!    after having completely read the file, initialize the time ten-
!    dencies to 0.0 for any tracers whose tinme tendencies were not
!    found on the restart file and enter a message in the output file.
!---------------------------------------------------------------------
        do nn=1,ntracers
          if (success(nn) ) then
          else
            call error_mesg ('donner_deep_mod', 'read_restart: &
                  &did not find tracer restart data for ' //  &
                  trim(Don_save%tracername(nn)) //  &
                  '; am initializing tendency to 0.0', NOTE)
            Don_save%tracer_tends(:,:,:,nn) = 0.0
          endif   
        end do
      endif  ! (do_donner_tracer)

!-------------------------------------------------------------------- 
!    close the restart file.
!--------------------------------------------------------------------- 
      call close_file (unit)

!--------------------------------------------------------------------- 



end subroutine read_restart

!#####################################################################

subroutine process_coldstart (Time)

!-----------------------------------------------------------------------
!    subroutine process_coldstart provides initialization that is needed
!    when the job is a donner_deep coldstart, or if the user-supplied 
!    restart file is not usable for a restart with the current code 
!    version.
!-----------------------------------------------------------------------

type(time_type), intent(in) :: Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!        Time      current time [ time_type, secs and days ]
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      integer  :: days, secs   ! components of current time

!---------------------------------------------------------------------
!    set the coldstart flag to .true.. set the time until the first cal-
!    culation call to donner_deep_mod, donner_deep calculation calls will
!    be every donner_deep_freq seconds after the start of the day.
!---------------------------------------------------------------------
      Initialized%coldstart = .true.
      call get_time (Time, secs, days)
      if (secs == 0) then    ! i.e., 00Z
        Initialized%conv_alarm = donner_deep_freq
      else 
        Initialized%conv_alarm = donner_deep_freq -   &
                                 MOD (secs, donner_deep_freq)
      endif

!----------------------------------------------------------------------
!    initialize the variables which must be returned from donner_deep_mod
!    on the first step when coldstarting.
!----------------------------------------------------------------------
      Don_save%cemetf            = 0.
      Don_save%cememf            = 0.
      Don_save%tracer_tends      = 0.
      Don_save%mass_flux         = 0.
      Don_save%cell_up_mass_flux = 0.
      Don_save%det_mass_flux     = 0.
      Don_save%dql_strat         = 0.
      Don_save%dqi_strat         = 0.
      Don_save%dqa_strat         = 0.
      Don_save%humidity_area     = 0.
      Don_save%humidity_factor   = 0.
      Don_save%tprea1            = 0.
      Don_save%parcel_disp       = 0.

!----------------------------------------------------------------------


end subroutine process_coldstart



!#####################################################################
! <SUBROUTINE NAME="read_restart_nc">
!  <OVERVIEW>
!    read_restart_nc reads a netcdf restart file containing donner_deep
!    restart information.
!  </OVERVIEW>
!  <DESCRIPTION>
!    read_restart_nc reads a netcdf restart file containing donner_deep
!    restart information.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_restart_nc
!  </TEMPLATE>
! </SUBROUTINE>
!


subroutine read_restart_nc (ntracers)

!-----------------------------------------------------------------------
!    subroutine read_restart_nc reads a netcdf restart file to obtain 
!    the variables needed upon experiment restart. 
!-----------------------------------------------------------------------

integer, intent(in) :: ntracers

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      ntracers    number of tracers being transported by the
!                  donner deep convection parameterization in this job
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      logical,         dimension(ntracers)  :: success
      integer,         dimension(:), allocatable :: ntindices
      type(fieldtype), dimension(:), allocatable :: tracer_fields

      character(len=64)     :: fname2='INPUT/donner_deep.res'
      character(len=64)     :: fname='INPUT/donner_deep.res.nc'
      character(len=128)    :: tname
      integer               :: ndim, natt, nvar, ntime
      integer               :: old_freq
      integer               :: n_alltracers, iuic, n
      logical               :: is_tracer_in_restart_file
      integer, dimension(4) :: siz
      logical               :: field_found, field_found2, field_found3,&
                               field_found4
      integer               :: it, jn, nn

!---------------------------------------------------------------------
!   local variables:
!
!        success          logical indicating if needed data for tracer n 
!                         was obtained from restart file
!        ntindices        array of all tracer indices
!        tracer_fields    field_type variable containing information on
!                         all restart file variables
!        fname2           restart file name without ".nc" appended, 
!                         needed as argument in call to mpp_open
!        fname            restart file name
!        tname            contains successive variable names from 
!                         restart file
!        ndim             number of dimensions in restart file
!        natt             number of attributes in restart file
!        nvar             number of variables in restart file
!        ntime            number of time levels in restart file
!        old_freq         donner_deep_freq as read from restart file;
!                         value used during previous job
!        n_alltracers     number of tracers registered with 
!                         tracer_manager_mod
!        iuic             unit number assigned to restart file
!        is_tracer_in_restart_file  
!                         should we stop searching the restart file 
!                         for the current tracer name because it has 
!                         been found ?
!        siz              sizes (each dimension) of netcdf variable 
!        field_found      is the requested variable in the restart file ?
!                         if it is not, then this is a reduced restart
!                         file
!        field_found2     is the requested variable in the restart file ?
!                         if it is not, then Don_save%det_mass_flux and
!                         Don_save%cell_up_mass_flux must be initialized
!        it, jn, nn       do-loop indices
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!    output a message indicating entrance into this routine.
!--------------------------------------------------------------------
      if (mpp_pe() == mpp_root_pe() ) then
        call error_mesg ('donner_deep_mod',  'read_restart_nc:&
             &Reading netCDF formatted restart file: &
                                 &INPUT/donner_deep.res.nc', NOTE)
      endif

!-------------------------------------------------------------------
!    read the values of conv_alarm when the restart file was written and
!    the frequency of calculating donner deep convection effects in the
!    job which wrote the file.
!-------------------------------------------------------------------
      call read_data(fname, 'conv_alarm', Initialized%conv_alarm,   &
                                                       no_domain=.true.)
      call read_data(fname, 'donner_deep_freq', old_freq,   &
                                                       no_domain=.true.)
  
!----------------------------------------------------------------------
!    call field_size to determine if variable cemetf is present in the
!    restart file.
!----------------------------------------------------------------------
      call field_size(fname, 'cemetf', siz, field_found=field_found)

!---------------------------------------------------------------------
!    if the frequency of calculating deep convection has changed, 
!    redefine the time remaining until the next calculation.
!---------------------------------------------------------------------
      if (donner_deep_freq /= old_freq) then
        Initialized%conv_alarm = Initialized%conv_alarm - old_freq +  &
                                 donner_deep_freq
        if (mpp_pe() == mpp_root_pe()) then
          call error_mesg ('donner_deep_mod', 'read_restart_nc:  &
                   &donner_deep time step has changed', NOTE)
        endif

!----------------------------------------------------------------------
!    if cemetf is not present, then this is a reduced restart file. it 
!    is not safe to change the frequency of calculating donner 
!    effects when reading a reduced restart file, so a fatal error is
!    generated.
!----------------------------------------------------------------------
        if (.not. field_found) then
          call error_mesg ('donner_deep_mod', 'read_restart_nc: &
           & cannot use reduced restart file and change donner_deep_freq&
           & within experiment and guarantee restart reproducibility', &
                                                                  FATAL)
        endif
      endif  !(donner_deep_freq /= old_freq)

!---------------------------------------------------------------------
!    read the restart data that is present in a full restart but absent
!    in a reduced restart.
!---------------------------------------------------------------------
      if (field_found) then
        call read_data (fname, 'cemetf',  Don_save%cemetf)
        call read_data (fname, 'cememf',  Don_save%cememf)            
        call read_data (fname, 'mass_flux', Don_save%mass_flux)
        call read_data (fname, 'dql_strat', Don_save%dql_strat)
        call read_data (fname, 'dqi_strat', Don_save%dqi_strat)
        call read_data (fname, 'dqa_strat', Don_save%dqa_strat)
        call read_data (fname, 'tprea1', Don_save%tprea1)       
        call read_data (fname, 'humidity_area', Don_save%humidity_area) 

!---------------------------------------------------------------------
!  determine if humidity_factor is in file. if it is, read the values 
!  into Don_Save%humidity_factor. if it is not (it is an older file), 
!  it is only required if donner_deep will not be called on the first 
!  step of this job.
!  if that is the case, stop with a fatal error; otherwise, continue on,
!  since humidity_factor will be calculated before it is used.
!---------------------------------------------------------------------
        call field_size(fname, 'humidity_factor', siz,   &
                                              field_found=field_found4)
        if (field_found4) then
          call read_data (fname, 'humidity_factor',  &
                                              Don_save%humidity_factor)
        else if (Initialized%conv_alarm > 0.0) then
          call error_mesg ('donner_deep_mod', &
             'cannot restart with this restart file unless donner_deep &
                &calculated on first step', FATAL)
        endif

!----------------------------------------------------------------------
!    determine if det_mass_flux is present in the file.
!----------------------------------------------------------------------
        call field_size(fname, 'det_mass_flux', siz,    &
                                               field_found=field_found2)

!----------------------------------------------------------------------
!    if it is present, then read det_mass_flux and cell_up_mass_flux.
!----------------------------------------------------------------------
        if (field_found2) then
          call read_data (fname, 'det_mass_flux', Don_save%det_mass_flux)
          call read_data (fname, 'cell_up_mass_flux',    &
                                              Don_save%cell_up_mass_flux)

!----------------------------------------------------------------------
!    if it is not present (an earlier version of this file), set 
!    det_mass_flux and cell_up_mass_flux to default values.
!----------------------------------------------------------------------
        else
          Don_save%det_mass_flux     = 0.0
          Don_save%cell_up_mass_flux = 0.0
        endif

!------------------------------------------------------------------
!    if tracers are to be transported, see if tendencies are available
!    in the restart file.
!------------------------------------------------------------------
        if (Initialized%do_donner_tracer) then

!---------------------------------------------------------------------
!    initialize a logical array indicating whether the data for each
!    tracer is available.
!---------------------------------------------------------------------
          success = .false.

!---------------------------------------------------------------------
!    open the restart file with mpp_open so that the unit number is 
!    available. obtain needed file characteristics by calling 
!    mpp_read_meta and  mpp_get_info. 
!---------------------------------------------------------------------
          call mpp_open(iuic, fname2, &
               action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_SINGLE )
          call mpp_read_meta (iuic)
          call mpp_get_info (iuic, ndim, nvar, natt, ntime)

!---------------------------------------------------------------------
!    obtain information on the file variables by calling mpp_get_fields.
!    it is returned in a field_type variable tracer_fields; the specific
!    information needed is the variable name.
!---------------------------------------------------------------------
          allocate (tracer_fields(nvar))
          if (mpp_pe() == mpp_root_pe()) then
            call mpp_get_fields (iuic, tracer_fields)
          endif

!---------------------------------------------------------------------
!    call get_number_tracers to determine how many tracers are registered
!    with tracer manager. allocate an array to hold their tracer indices.
!    call get_tracer_indices to retrieve the tracer indices. 
!---------------------------------------------------------------------
          call get_number_tracers (MODEL_ATMOS, num_tracers=n_alltracers)
          allocate (ntindices(n_alltracers))
          call get_tracer_indices (MODEL_ATMOS, ind=ntindices)

!----------------------------------------------------------------------
!    loop over the tracers, obtaining their names via a call to
!    get_tracer_names. bypass those tracers known to not be transported
!    by donner convection.
!----------------------------------------------------------------------
          do it=1,n_alltracers
            call get_tracer_names (MODEL_ATMOS, ntindices(it), tname)
            if (tname == "sphum"  ) cycle
            if (tname == "liq_wat") cycle
            if (tname == "ice_wat") cycle
            if (tname == "cld_amt") cycle

!--------------------------------------------------------------------
!    initialize a logical indicating whether this tracer is in the 
!    restart file.
!--------------------------------------------------------------------
            is_tracer_in_restart_file = .FALSE.

!---------------------------------------------------------------------
!    loop over the variables in the restart file to determine if the
!    current tracer's time tendency field is present.
!---------------------------------------------------------------------
            do jn=1,nvar 
              if (lowercase (trim(mpp_get_field_name(tracer_fields(jn)))) ==   &
                  lowercase ('tracer_tends_' // trim(tname)) ) then 

!---------------------------------------------------------------------
!    if tracer tendency is in restart file, write a message. set the 
!    logical flag indicating such to .true..
!---------------------------------------------------------------------
                if (mpp_pe() == mpp_root_pe() )  then
                  print *,'tracer_tends_' // trim(tname), ' found!'
                endif
                is_tracer_in_restart_file = .TRUE.

!---------------------------------------------------------------------
!    loop over the tracers being transported by donner convection in this
!    job to determine if this tracer is one of those being transported.
!    determine the tracer index in tracername array corresponding to 
!    this tracer.
!---------------------------------------------------------------------
                do nn=1,ntracers
                  if (lowercase( 'tracer_tends_' // trim(tname) ) == &
                      'tracer_tends_' // Don_save%tracername(nn) )  then
                  
!---------------------------------------------------------------------
!    if data for this tracer is needed, read data into proper section of
!    array tracer_tends. set the logical flag for this tracer indicating 
!    successful retrieval. exit this loop.
!---------------------------------------------------------------------
                    call read_data (fname,   &
                                  'tracer_tends_' // trim(tname),   &
                                   Don_save%tracer_tends(:,:,:,nn))
                    success(nn) = .true.
                    exit
                  endif 
                end do  ! (nn)
              endif

!---------------------------------------------------------------------
!    if desired tracer has been found, stop searching the restart file
!    variables for this tracer and cycle to begin searching the restart
!    file for the next field_table tracer.
!---------------------------------------------------------------------
              if (is_tracer_in_restart_file) exit
            end do !  (jn)
          end do ! (it)

!---------------------------------------------------------------------
!    initialize the time tendencies to 0.0 for any tracers that are to
!    be transported and whose time tendencies were not found on the 
!    restart file.  enter a message in the output file.
!---------------------------------------------------------------------
          do nn=1,ntracers
            if (success(nn) ) then
            else
              call error_mesg ('donner_deep_mod', 'read_restart_nc: &
                  &did not find tracer restart data for ' //  &
                  trim(Don_save%tracername(nn)) //  &
                  '; am initializing tendency to 0.0', NOTE)
              Don_save%tracer_tends(:,:,:,nn) = 0.0
            endif   
          end do

!----------------------------------------------------------------------
!    deallocate local variables.
!----------------------------------------------------------------------
          deallocate (ntindices)
          deallocate (tracer_fields)
        endif  ! (do_donner_tracer)
      endif  ! (field_found)

!---------------------------------------------------------------------
!    read the restart data that is present in both full and reduced
!    restart files.
!---------------------------------------------------------------------
      call read_data (fname, 'parcel_disp', Don_save%parcel_disp)
      call read_data (fname, 'lag_temp',    Don_save%lag_temp)     
      call read_data (fname, 'lag_vapor',   Don_save%lag_vapor)     
      call read_data (fname, 'lag_press',   Don_save%lag_press)     

!---------------------------------------------------------------------




end subroutine read_restart_nc



!#####################################################################

subroutine process_monitors (idf, jdf, nlev, ntracers, axes, Time)

integer,                       intent(in)  :: idf, jdf, nlev, ntracers
integer,         dimension(4), intent(in)  :: axes
type(time_type),               intent(in)  :: Time

!-------------------------------------------------------------------
!  local variables:

      integer             :: k, n, nn, nx, nc
      logical             :: flag, success
      integer             :: nfields, model, num_methods
      character(len=200)  :: method_name, field_type, method_control,&
                             field_name, list_name
      character(len=32)   :: path_name = '/atmos_mod/don_deep_monitor/'

!---------------------------------------------------------------------
!    determine if and how many output variables are to be monitored. 
!    set a flag indicating if monitoring is activated.
!---------------------------------------------------------------------
      call field_manager_init (nfields)
      nx = 0
      do n=1,nfields
        call get_field_info (n, field_type, field_name, model, &
                             num_methods)
        if (trim(field_type) == 'don_deep_monitor') then
          nx = nx + 1
        endif
      end do
      if (nx > 0) then
        Initialized%monitor_output = .true.
      else
        Initialized%monitor_output = .false.
      endif

!---------------------------------------------------------------------
!    allocate arrays needed for each monitored variable. 
!---------------------------------------------------------------------
      if (Initialized%monitor_output) then
        allocate (Initialized%Don_monitor(nx))
        allocate (id_extremes(nx))
        allocate (id_hits(nx))

!---------------------------------------------------------------------
!    read the field_table to determine the nature of the monitors
!    requested.
!---------------------------------------------------------------------
        nx = 1
        do n = 1,nfields
          call get_field_info (n, field_type, field_name, model, &
                               num_methods)

!---------------------------------------------------------------------
!    define the list name used by field_manager_mod to point to 
!    monitored variables.
!---------------------------------------------------------------------
          if (trim(field_type) == 'don_deep_monitor') then
            list_name = trim(path_name) // trim(field_name) // '/'

!--------------------------------------------------------------------
!    place name of field in don_monitor_type variable.
!--------------------------------------------------------------------
            Initialized%Don_monitor(nx)%name = trim(field_name)

!--------------------------------------------------------------------
!    map the field name to the list of acceptable field names. store
!    the index of this field name in the don_monitor_type variable.
!    note that any tracer variables need to have 'tr_' as the first
!    three characters in their name to allow proper processing. store
!    the appropriate tracer index for any tracer arrays.
!--------------------------------------------------------------------
            if (trim(field_name(1:3)) == 'tr_') then
              select case (trim(field_name(4:9)))
                case ('rn_ten')
                  Initialized%Don_monitor(nx)%index = RADON_TEND
                  success = .false.
                  do nc=1,ntracers
                    if (trim(Don_save%tracername(nc)) == 'radon') then
                      Initialized%Don_monitor(nx)%tracer_index = nc
                      success = .true.
                      exit
                    endif
                  end do
                  if (.not. success) then
                    call error_mesg ('donner_deep_mod', &
                     'not able to find "radon" tracer index', FATAL)
                  endif
                case default
                  call error_mesg ('donner_deep_mod', &
                 'tracer variable name in field_table don_deep_monitor &
                                             &type is invalid', FATAL)
              end select

!---------------------------------------------------------------------
!    for non-tracer variables, set the tracer index to an arbitrary 
!    value.
!---------------------------------------------------------------------
            else
              Initialized%Don_monitor(nx)%tracer_index = 0
              select case (trim(field_name(1:6)))
                case ('det_ma')
                  Initialized%Don_monitor(nx)%index = DET_MASS_FLUX
                case ('mass_f')
                  Initialized%Don_monitor(nx)%index = MASS_FLUX
                case ('cell_u')
                  Initialized%Don_monitor(nx)%index =   &
                                                  CELL_UPWARD_MASS_FLUX
                case ('temp_f')
                  Initialized%Don_monitor(nx)%index = TEMP_FORCING
                case ('moistu')
                  Initialized%Don_monitor(nx)%index = MOIST_FORCING
                case ('precip')
                  Initialized%Don_monitor(nx)%index = PRECIP
                case ('freeze')
                  Initialized%Don_monitor(nx)%index = FREEZING
                case default
                  call error_mesg ('donner_deep_mod', &
                      'variable name in field_table don_deep_monitor &
                                              &type is invalid', FATAL)
              end select
            endif

!---------------------------------------------------------------------
!    read the units for this variable from the field_table entry.
!    if the units method is missing, set units to be 'missing'.
!---------------------------------------------------------------------
            flag = fm_query_method (trim(list_name) //  'units',    &
                                    method_name, method_control)
            if (flag) then
              Initialized%Don_monitor(nx)%units = trim(method_name)
            else
              Initialized%Don_monitor(nx)%units = 'missing'
            endif

!---------------------------------------------------------------------
!    determine the type of limit being imposed for this variable from 
!    the field_table entry.
!---------------------------------------------------------------------
            flag = fm_query_method (trim(list_name) // 'limit_type',  &
                                    method_name, method_control)

!----------------------------------------------------------------------
!    include the limit_type for this variable in its don_monitor type
!    variable.
!    register diagnostics associated with the monitored output fields
!    (extreme values and number of times threshold was exceeeded).
!----------------------------------------------------------------------
            if ( flag) then
              if (trim(method_name) == 'maxmag') then
                Initialized%Don_monitor(nx)%initial_value = 0.0
                Initialized%Don_monitor(nx)%limit_type =   MAXMAG
                id_extremes(nx) = register_diag_field (mod_name,   &
                  'maxmag_'// trim(Initialized%Don_monitor(nx)%name),  &
                   axes(1:3),  Time,  'maxmag values of ' // &
                            trim(Initialized%Don_monitor(nx)%name),  &
                   Initialized%Don_monitor(nx)%units,   &
                   mask_variant = .true., missing_value=missing_value)
                id_hits(nx) = register_diag_field (mod_name,   &
                  'num_maxmag_'// &
                              trim(Initialized%Don_monitor(nx)%name) , &
                   axes(1:3),  Time,    &
                   '# of times that magnitude of '&
                     // trim(Initialized%Don_monitor(nx)%name) //  &
                 ' > ' // trim(method_control(2:)) // ' ' // &
                    trim(Initialized%Don_monitor(nx)%units) ,  &
                   'number', mask_variant = .true., & 
                   missing_value=missing_value)
              else if (trim(method_name) == 'minmag') then
                Initialized%Don_monitor(nx)%initial_value = 1.0e30
                Initialized%Don_monitor(nx)%limit_type =   MINMAG
                id_extremes(nx) = register_diag_field (mod_name,   &
                  'minmag_'// trim(Initialized%Don_monitor(nx)%name),  &
                   axes(1:3),  Time,  'minmag values of ' // &
                            trim(Initialized%Don_monitor(nx)%name),  &
                   Initialized%Don_monitor(nx)%units,   &
                   mask_variant = .true., missing_value=missing_value)
                id_hits(nx) = register_diag_field (mod_name,   &
                  'num_minmag_'//     &
                             trim(Initialized%Don_monitor(nx)%name) , &
                  axes(1:3),  Time,    &
                   '# of times that magnitude of '&
                     // trim(Initialized%Don_monitor(nx)%name) //  &
                ' < ' // trim(method_control(2:)) // ' ' // &
                    trim(Initialized%Don_monitor(nx)%units) ,  &
                  'number', mask_variant = .true., & 
                  missing_value=missing_value)
              else if (trim(method_name) == 'minval') then
                Initialized%Don_monitor(nx)%initial_value = 1.0e30
                Initialized%Don_monitor(nx)%limit_type =   MINVAL
                id_extremes(nx) = register_diag_field (mod_name,   &
                  'minval_'// trim(Initialized%Don_monitor(nx)%name),  &
                   axes(1:3),  Time,  'minimum values of ' // &
                            trim(Initialized%Don_monitor(nx)%name),  &
                  Initialized%Don_monitor(nx)%units,   &
                  mask_variant = .true., missing_value=missing_value)
                id_hits(nx) = register_diag_field (mod_name,   &
                  'num_minval_'//   &
                             trim(Initialized%Don_monitor(nx)%name) , &
                   axes(1:3),  Time,    &
                   '# of times that value of '&
                     // trim(Initialized%Don_monitor(nx)%name) //  &
                ' < ' // trim(method_control(2:)) // ' ' // &
                    trim(Initialized%Don_monitor(nx)%units) ,  &
                  'number', mask_variant = .true., & 
                  missing_value=missing_value)
              else if (trim(method_name) == 'maxval') then
                Initialized%Don_monitor(nx)%initial_value = -1.0e30
                Initialized%Don_monitor(nx)%limit_type = MAXVAL 
                id_extremes(nx) = register_diag_field (mod_name,   &
                  'maxval_'// trim(Initialized%Don_monitor(nx)%name),  &
                  axes(1:3),  Time,  'maximum values of ' // &
                            trim(Initialized%Don_monitor(nx)%name),  &
                  Initialized%Don_monitor(nx)%units,  &
                  mask_variant = .true., missing_value=missing_value)
                id_hits(nx) = register_diag_field (mod_name,   &
                  'num_maxval_'//    &
                             trim(Initialized%Don_monitor(nx)%name) , &
                  axes(1:3),  Time,    &
                   '# of times that value of '&
                     // trim(Initialized%Don_monitor(nx)%name) //  &
                    ' > ' // trim(method_control(2:)) // ' ' // &
                    trim(Initialized%Don_monitor(nx)%units) ,  &
                  'number', mask_variant = .true., & 
                  missing_value=missing_value)
              else
                call error_mesg ('donner_deep_mod', &
                    'invalid limit_type for monitored variable', FATAL)
              endif

!----------------------------------------------------------------------
!    if limit_type not in field_table, set it to look for maximum
!    magnitude.
!----------------------------------------------------------------------
            else
              Initialized%Don_monitor(nx)%initial_value = 0.0
              Initialized%Don_monitor(nx)%limit_type =   MAXMAG
              id_extremes(nx) = register_diag_field (mod_name,   &
                  'maxmag_'// trim(Initialized%Don_monitor(nx)%name),  &
                   axes(1:3),  Time,  'maxmag values of ' // &
                            trim(Initialized%Don_monitor(nx)%name),  &
                   Initialized%Don_monitor(nx)%units,   &
                   mask_variant = .true., missing_value=missing_value)
              id_hits(nx) = register_diag_field (mod_name,   &
                  'num_maxmag_'// &
                              trim(Initialized%Don_monitor(nx)%name) , &
                   axes(1:3),  Time,    &
                   '# of times that magnitude of '&
                     // trim(Initialized%Don_monitor(nx)%name) //  &
                ' > ' // trim(method_control(2:)) // ' ' // &
                    trim(Initialized%Don_monitor(nx)%units) ,  &
                   'number', mask_variant = .true., & 
                   missing_value=missing_value)
            endif

!----------------------------------------------------------------------
!    obtain the magnitude of the limit being monitored for this 
!    variable from the field_table. 
!----------------------------------------------------------------------
            flag = parse (method_control, 'value',   &
                            Initialized%Don_monitor(nx)%threshold ) > 0

!----------------------------------------------------------------------
!    if no limit_type and / or value has been given, the
!    field will be flagged for magnitudes  > 0.0, i.e., if deep 
!    convection has affected the point.
!----------------------------------------------------------------------
            if ( .not. flag) then
              Initialized%Don_monitor(nx)%threshold = 0.0
            endif

!-------------------------------------------------------------------
!    allocate and initialize arrays to hold the extrema and a count of 
!    times the threshold was exceeded at each point.
!-------------------------------------------------------------------
            allocate (Initialized%Don_monitor(nx)%extrema(idf,jdf,nlev))
            Initialized%Don_monitor(nx)%extrema(:,:,:) =  &
                         Initialized%Don_monitor(nx)%initial_value
            allocate (Initialized%Don_monitor(nx)%hits(idf,jdf,nlev))
            Initialized%Don_monitor(nx)%hits(:,:,:) = 0.0
            nx = nx + 1
          endif
        end do
      endif 

end subroutine process_monitors



!#####################################################################


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      2. ROUTINES CALLED BY DONNER_DEEP
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#######################################################################

subroutine donner_column_control (is, ie, js, je, Time)                

!---------------------------------------------------------------------
!    subroutine donner_column_control returns the number, location
!    (processor and window indices) and output units associated with 
!    any diagnostic columns requested within the current physics window.
!---------------------------------------------------------------------

integer,                       intent(in)   :: is, ie, js, je
type(time_type),               intent(in)   :: Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points
!                    in this physics window (processor coordinates)
!     Time           current model time [ time_type, days, seconds ]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer  :: isize      !   i-dimension of physics window
      integer  :: jsize      !   j-dimension of physics window
      integer  :: nn, j, i   !   do-loop indices

!--------------------------------------------------------------------
!    define the sizes of the current physics window's horizontal
!    dimensions.
!--------------------------------------------------------------------
      isize = ie - is + 1
      jsize = je - js + 1

!-------------------------------------------------------------------
!    initialize the output variables.
!-------------------------------------------------------------------
      Col_diag%i_dc(:) = -99
      Col_diag%j_dc(:) = -99
      Col_diag%unit_dc(:) = -1
      Col_diag%jgl_dc(:) = -99
      Col_diag%igl_dc(:) = -99
      Col_diag%ncols_in_window = 0

!--------------------------------------------------------------------
!    if any requested diagnostic columns are present within the current
!    physics window, and if it is at or past the time to start output-
!    ting column diagnostics, save the relevant variables describing
!    those diagnostic columns in arrays to be returned to the calling
!    routine. call column_diagnostics_header to write the file header
!    for the diagnostic columns in this window. 
!--------------------------------------------------------------------
      if (Col_diag%num_diag_pts > 0) then
        if (Time >= Time_col_diagnostics) then
          do nn=1,Col_diag%num_diag_pts
            do j=1,jsize      
              if (js + j - 1 == col_diag_j(nn)) then
                do i=1,isize       
                  if (is + i - 1 == col_diag_i(nn)) then
                    Col_diag%ncols_in_window =   &
                                           Col_diag%ncols_in_window + 1
                    Col_diag%i_dc(Col_diag%ncols_in_window) = i
                    Col_diag%j_dc(Col_diag%ncols_in_window) = j
                    Col_diag%igl_dc(COl_diag%ncols_in_window) =  &
                                                          col_diag_i(nn)
                    Col_diag%jgl_dc(Col_diag%ncols_in_window) =   &
                                                           col_diag_j(nn)
                    Col_diag%unit_dc(Col_diag%ncols_in_window) =   &
                                                        col_diag_unit(nn)
                    call column_diagnostics_header &
                            (mod_name, col_diag_unit(nn), Time, nn,  &
                             col_diag_lon, col_diag_lat, col_diag_i,  &
                             col_diag_j)
                  endif
                end do  ! (i loop)
              endif
            end do  ! (j loop) 
          end do  ! (num_diag_pts loop)
        endif  ! (Time >= starting time)
      endif ! (num_diag_pts > 0)

!---------------------------------------------------------------------

end subroutine donner_column_control



!######################################################################

subroutine donner_deep_netcdf (is, ie, js, je, Nml, Time, Don_conv, Don_cape,&
                               parcel_rise, pmass, total_precip, &
                               Don_budgets, &
                               temperature_forcing, moisture_forcing)  

!---------------------------------------------------------------------
!    subroutine donner_deep_netcdf sends the fields requested in the
!    diag_table to diag_manager_mod so that they may be appropriately
!    processed for output.
!---------------------------------------------------------------------

integer,                intent(in) :: is, ie, js, je
type(time_type),        intent(in) :: Time
type(donner_nml_type), intent(in) :: Nml   
type(donner_conv_type), intent(in) :: Don_conv
type(donner_budgets_type), intent(in) :: Don_budgets
type(donner_cape_type), intent(in) :: Don_cape
real, dimension(:,:,:), intent(in) :: pmass, temperature_forcing,&
                                      moisture_forcing
real, dimension(:,:),   intent(in) :: parcel_rise, total_precip

!---------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     Time           current time (time_type)
!     Don_conv       donner_convection_type derived type variable con-
!                    taining diagnostics describing the nature of the 
!                    convection produced by the donner parameterization
!     Don_cape       donner_cape type derived type variable containing
!                    diagnostics related to the cape calculation assoc-
!                    iated with the donner convection parameterization
!     temperature_forcing  
!                    temperature tendency due to donner convection
!                    [ deg K / sec ]
!     moisture_forcing  
!                    vapor mixing ratio tendency due to donner 
!                    convection [ kg(h2o) / (kg(dry air) sec ) ]
!     pmass          mass per unit area within the grid box
!                    [ kg (air) / (m**2) ]
!     parcel_rise    accumulated vertical displacement of a near-surface
!                    parcel as a result of the lowest model level omega 
!                    field [ Pa ]
!     total_precip   total precipitation rate produced by the
!                    donner parameterization [ mm / day ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      real, dimension (ie-is+1, je-js+1)  :: tempdiag, tempdiag2, tempdiag3  
                           ! array used to hold various data fields being
                           ! sent to diag_manager_mod
      logical :: used      ! logical indicating data has been received 
                           ! by diag_manager_mod 
      integer :: nlev      ! number of large-scale model layers
      integer :: ntr       ! number of tracers transported by the
                           ! donner deep convection parameterization
      integer :: k, n, nn  ! do-loop indices

!----------------------------------------------------------------------
!    define the number of model layers (nlev) and number of transported
!    tracers (ntr).
!----------------------------------------------------------------------
      nlev = size (pmass,3)
      ntr  = size (Don_conv%qtren1,4)

!---------------------------------------------------------------------
!    send the 3D convective output variables to diag_manager_mod.
!!   NOTE: effective with code mod lima_donnermod3_rsh (7-19-05) the
!!         temperature and moisture forcing fields passed to diag_manager
!!         (id_cemetf_deep, id_cememf_deep) are the total convective
!!         forcings calculated by the donner parameterization. Previous
!!         code versions run in models in which strat_cloud_mod was 
!!         activated output the forcing fields less the terms related to 
!!         the flux convergence of the large-scale condensate and the 
!!         mesoscale detrainment.
!---------------------------------------------------------------------

!   total convective temperature forcing:
      used = send_data (id_cemetf_deep, Don_conv%conv_temp_forcing,  &
                        Time, is, js, 1)

!   cell entropy flux convergence:
      used = send_data (id_ceefc_deep, Don_conv%ceefc, Time, is, js, 1)

!   cell condensation / evaporation:
      used = send_data (id_cecon_deep, Don_conv%cecon, Time, is, js, 1)

!   cell moisture flux convergence:
      used = send_data (id_cemfc_deep, Don_conv%cemfc, Time, is, js, 1)

!   total convective moistening forcing:
      used = send_data (id_cememf_deep, Don_conv%conv_moist_forcing,  &
                        Time, is, js, 1)

!   total convective moistening rate after adjustnment for negative 
!   vapor mixing ratio:
      used = send_data (id_cememf_mod_deep, Don_conv%cememf_mod,   &
                        Time, is, js, 1)

!   cell + mesoscale cloud fraction:
      used = send_data (id_cual_deep, Don_conv%cual, Time, is, js, 1)

!   heating rate due to freezing:
      used = send_data (id_fre_deep, Don_conv%fre, Time, is, js, 1)

!   heating rate due to melting:
      used = send_data (id_elt_deep, Don_conv%elt, Time, is, js, 1)

!   deposition in mesoscale updraft:
      used = send_data (id_cmus_deep, Don_conv%cmus, Time, is, js, 1)

!   evaporation in convective downdrafts:
      used = send_data (id_ecds_deep, Don_conv%ecds, Time, is, js, 1)

!   evaporation / sublimation in convective updrafts:
      used = send_data (id_eces_deep, Don_conv%eces, Time, is, js, 1)

!   sublimation in mesoscale downdrafts:
      used = send_data (id_emds_deep, Don_conv%emds, Time, is, js, 1)

!   sublimation in mesoscale updrafts:
      used = send_data (id_emes_deep, Don_conv%emes, Time, is, js, 1)

!   mesoscale moisture flux convergence:
      used = send_data (id_qmes_deep, Don_conv%mrmes, Time, is, js, 1)

!   transfer of vapor from cells to mesoscale:
      used = send_data (id_wmps_deep, Don_conv%wmps, Time, is, js, 1)

!   deposition of vapor from cells to mesoscale:
      used = send_data (id_wmms_deep, Don_conv%wmms, Time, is, js, 1)

!   mesoscale entropy flux convergence:
      used = send_data (id_tmes_deep, Don_conv%tmes, Time, is, js, 1)

!   mass flux in mesoscale downdrafts:
      used = send_data (id_dmeml_deep, Don_conv%dmeml, Time, is, js, 1)

!   mass flux in cell updrafts:
      used = send_data (id_uceml_deep, Don_conv%uceml, Time, is, js, 1)

!   detrained mass flux:
      used = send_data (id_detmfl_deep, Don_conv%detmfl, Time, is, js, 1)

!   mass flux in mesoscale updrafts:
      used = send_data (id_umeml_deep, Don_conv%umeml, Time, is, js, 1)

!   mesoscale ice mixing ratio:
      used = send_data (id_xice_deep, Don_conv%xice, Time, is, js, 1)

!   mesoscale liquid mass mixing ratio
      used = send_data (id_xliq_deep, Don_conv%xliq, Time, is, js, 1)

!   mesoscale ice generalized effective size:
      used = send_data (id_dgeice_deep, Don_conv%dgeice,      &
                        Time, is, js, 1)

!   cell ice mixing ratio:
      used = send_data (id_cuqi_deep, Don_conv%cuqi, Time, is, js, 1)

!   cell liquid mixing ratio:
      used = send_data (id_cuql_deep, Don_conv%cuql, Time, is, js, 1)

!   cell liquid generalized effective size:
      used = send_data (id_dgeliq_deep, Don_conv%cell_liquid_eff_diam, &
                        Time, is, js, 1)

     if (Nml%do_budget_analysis) then
       do n=1,Don_budgets%N_WATER_BUDGET
         if (id_water_budget(n) > 0) then
            used = send_data (id_water_budget(n), &
                              Don_budgets%water_budget(:,:,:,n), &
                              Time, is, js, 1)
         endif
       end do
       do n=1,Don_budgets%N_PRECIP_TYPES
         do nn=1,Don_budgets%N_PRECIP_PATHS
           if (id_precip_budget(nn,n) > 0) then
             used = send_data (id_precip_budget(nn,n), &
                               Don_budgets%precip_budget(:,:,:,nn,n), &
                               Time, is, js, 1)
           endif
         end do
       end do
       do n=1,Don_budgets%N_ENTHALPY_BUDGET
         if (id_enthalpy_budget(n) > 0) then
           used = send_data (id_enthalpy_budget(n),   &
                             Don_budgets%enthalpy_budget(:,:,:,n), &
                             Time, is, js, 1)
         endif
       end do
       do n=1,Don_budgets%N_WATER_BUDGET
         tempdiag(:,:) = 0.
         do k=1,nlev
           tempdiag(:,:) = tempdiag(:,:) + &
                           Don_budgets%water_budget(:,:,k,n)* &
                                                     pmass(:,:,k)/1000.
         end do
         if (id_ci_water_budget(n) > 0) then
           used = send_data (id_ci_water_budget(n), tempdiag, &
                             Time, is, js)
         endif
       end do
       tempdiag3(:,:) = 0.
       do n=1,Don_budgets%N_PRECIP_TYPES
         do nn=1,Don_budgets%N_PRECIP_PATHS
           tempdiag(:,:) = 0.
           do k=1,nlev
             tempdiag(:,:) = tempdiag(:,:) + &
                             Don_budgets%precip_budget(:,:,k,nn,n)* &
                                                           pmass(:,:,k)
           end do
           if (id_ci_precip_budget(nn,n) > 0) then
             used = send_data (id_ci_precip_budget(nn,n), tempdiag, &
                               Time, is, js)
           endif
           tempdiag3(:,:) = tempdiag3(:,:) + tempdiag(:,:)
         end do
       end do
       do n=1,Don_budgets%N_ENTHALPY_BUDGET
         tempdiag(:,:) = 0.
         do k=1,nlev
           tempdiag(:,:) = tempdiag(:,:) +  &
                           Don_budgets%enthalpy_budget(:,:,k,n)* &
                                                    pmass(:,:,k)*CP_AIR
         end do
         if (id_ci_enthalpy_budget(n) > 0) then
           used = send_data (id_ci_enthalpy_budget(n), tempdiag, &
                             Time, is, js)
         endif
       end do
           
        
       tempdiag2(:,:) = 0.
       tempdiag(:,:) = 0.
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) +  &
                         (Don_budgets%precip_budget(:,:,k,2,1) +  &
                          Don_budgets%precip_budget(:,:,k,4,1))* &
                                                 Param%hls*pmass(:,:,k)
       end do
       if (id_ci_prcp_heat_frz_cell > 0) then
         used = send_data (id_ci_prcp_heat_frz_cell, tempdiag, &
                           Time, is, js)
       endif
       tempdiag2 = tempdiag2 + tempdiag
           
       tempdiag(:,:) = 0.
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) +  &
                         (Don_budgets%precip_budget(:,:,k,1,1) +   &
                          Don_budgets%precip_budget(:,:,k,3,1) + &
                          Don_budgets%precip_budget(:,:,k,5,1))* &
                                                 Param%hlv*pmass(:,:,k)
       end do
       if (id_ci_prcp_heat_liq_cell > 0) then
         used = send_data (id_ci_prcp_heat_liq_cell, tempdiag, &
                           Time, is, js)
       endif
       tempdiag2 = tempdiag2 + tempdiag
           
       tempdiag(:,:) = 0.
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) +  &
                         (Don_budgets%precip_budget(:,:,k,2,2) + &
                          Don_budgets%precip_budget(:,:,k,4,2) + &
                          Don_budgets%precip_budget(:,:,k,2,3) + &
                          Don_budgets%precip_budget(:,:,k,4,3))* &
                                                 Param%hls*pmass(:,:,k)
       end do
       if (id_ci_prcp_heat_frz_meso > 0) then
         used = send_data (id_ci_prcp_heat_frz_meso, tempdiag, &
                           Time, is, js)
       endif
       tempdiag2 = tempdiag2 + tempdiag
           
       tempdiag(:,:) = 0.
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) +  &
                         (Don_budgets%precip_budget(:,:,k,1,2) +   &
                          Don_budgets%precip_budget(:,:,k,3,2) +  &
                          Don_budgets%precip_budget(:,:,k,5,2) + &
                          Don_budgets%precip_budget(:,:,k,1,3) +   &
                          Don_budgets%precip_budget(:,:,k,3,3) +  &
                          Don_budgets%precip_budget(:,:,k,5,3))* &
                                                  Param%hlv*pmass(:,:,k)
       end do
       if (id_ci_prcp_heat_liq_meso > 0) then
         used = send_data (id_ci_prcp_heat_liq_meso, tempdiag, &
                           Time, is, js)
       endif
       tempdiag2 = tempdiag2 + tempdiag
       if ( id_ci_prcp_heat_total > 0) then
         used = send_data (id_ci_prcp_heat_total, tempdiag2, &
                           Time, is, js)
       endif
       if (id_ci_prcp_total > 0) then
         used = send_data (id_ci_prcp_total, tempdiag3, &
                           Time, is, js)
       endif
       if ( id_leff > 0) then
         used = send_data(id_leff, tempdiag2/(tempdiag3+1.0e-40), &
                           Time, is, js)
       endif
           
     endif

!--------------------------------------------------------------------
!    send the tracer-related arrays to diag_manager_mod.
!--------------------------------------------------------------------
      do n=1,ntr    

!   tracer tendency due to cells:
        if (id_qtren1(n) > 0) then
        used = send_data (id_qtren1(n), Don_conv%qtren1(:,:,:,n), &
                          Time, is, js, 1)
        endif

!   tracer tendency due to mesoscale:
         if (id_qtmes1(n) > 0) then
        used = send_data (id_qtmes1(n), Don_conv%qtmes1(:,:,:,n),   &
                          Time, is, js, 1)
        endif

!   tracer tendency due to mesoscale redistribution:
        if (id_wtp1(n) > 0) then
        used = send_data (id_wtp1(n), Don_conv%wtp1(:,:,:,n),     &
                          Time, is, js, 1)
        endif

!   tracer tendency due to deep convective wet deposition:
       if (id_total_wet_dep(n) > 0) then
     used = send_data (id_total_wet_dep(n), Don_conv%wetdept(:,:,:,n), &
                            Time, is, js, 1)
        endif
!   tracer tendency due to wet deposition in mesoscale updrafts:
       if ( id_meso_wet_dep(n) > 0) then
     used = send_data (id_meso_wet_dep(n), Don_conv%wetdepm(:,:,:,n), &
                            Time, is, js, 1)
      endif
 
!   tracer tendency due to wet deposition in cells:
      if (id_cell_wet_dep(n) > 0) then
     used = send_data (id_cell_wet_dep(n), Don_conv%wetdepc(:,:,:,n), &
                           Time, is, js, 1)
      endif

!   total tracer tendency:
      if (id_qtceme(n) > 0) then
        used = send_data (id_qtceme(n), Don_conv%qtceme(:,:,:,n), &
                          Time, is, js, 1)
      endif

!---------------------------------------------------------------------
!    define the column-integrated tracer tendency due to convective
!    cells, in units of kg (tracer) / (m**2 sec). send it to 
!    diag_manager_mod.
!---------------------------------------------------------------------
        tempdiag = 0.0
        do k=1,nlev
          tempdiag(:,:) = tempdiag(:,:) + Don_conv%qtren1(:,:,k,n)* &
                          pmass(:,:,k)
        end do
        if (id_qtren1_col(n) > 0) then
        used = send_data (id_qtren1_col(n), tempdiag, Time, is, js)
        endif

!---------------------------------------------------------------------
!    define the column-integrated tracer tendency due to mesoscale circ-
!    ulation, in units of kg (tracer) / (m**2 sec). send it to 
!    diag_manager_mod.
!---------------------------------------------------------------------
        tempdiag = 0.0
        do k=1,nlev
          tempdiag(:,:) = tempdiag(:,:) + Don_conv%qtmes1(:,:,k,n)* &
                          pmass(:,:,k)
        end do
        if (id_qtmes1_col(n) > 0) then
        used = send_data (id_qtmes1_col(n), tempdiag, Time, is, js)
        endif

!---------------------------------------------------------------------
!    define the column-integrated tracer redistribution due to meso-
!    scale circulation, in units of kg (tracer) / (m**2 sec). send it 
!    to diag_manager_mod.
!---------------------------------------------------------------------
        tempdiag = 0.0
        do k=1,nlev
          tempdiag(:,:) = tempdiag(:,:) + Don_conv%wtp1(:,:,k,n)*   &
                          pmass(:,:,k)
        end do
        if (id_wtp1_col(n) > 0) then
        used = send_data (id_wtp1_col(n), tempdiag, Time, is, js)
        endif

!---------------------------------------------------------------------
!    define the column-integrated tracer change due to wet deposition in
!    deep convection (cells and mesoscale) in units of kg (tracer) / 
!    (m**2 sec). send it to diag_manager_mod.
!---------------------------------------------------------------------
        tempdiag = 0.0
        do k=1,nlev
          tempdiag(:,:) = tempdiag(:,:) + Don_conv%wetdept(:,:,k,n)*   &
                          pmass(:,:,k)
        end do
        if (id_total_wet_dep_col(n) > 0) then
        used = send_data (id_total_wet_dep_col(n), tempdiag, Time,  &
                                                                 is, js)
        endif

!---------------------------------------------------------------------
!    define the column-integrated tracer change due to wet deposition in
!    mesoscale updrafts, in units of kg (tracer) / (m**2 sec). send it 
!    to diag_manager_mod.
!---------------------------------------------------------------------
       tempdiag = 0.0
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) + Don_conv%wetdepm(:,:,k,n)*   &
                         pmass(:,:,k)
       end do
       if (id_meso_wet_dep_col(n) > 0) then
       used = send_data (id_meso_wet_dep_col(n), tempdiag, Time,  &
                                                                 is, js)
       endif

!---------------------------------------------------------------------
!    define the column-integrated tracer change due to wet deposition 
!    by convective cells, in units of kg (tracer) / (m**2 sec). send it 
!    to diag_manager_mod.
!---------------------------------------------------------------------
       tempdiag = 0.0
       do k=1,nlev
         tempdiag(:,:) = tempdiag(:,:) + Don_conv%wetdepc(:,:,k,n)*   &
                         pmass(:,:,k)
       end do
        if (id_cell_wet_dep_col(n) > 0) then
       used = send_data (id_cell_wet_dep_col(n), tempdiag, Time,  &
                                                                 is, js)
        endif

!-----------------------------------------------------------------
!    define the column-integrated total tracer tendency, in units of 
!    kg (tracer) / (m**2 sec). send it to diag_manager_mod.
!---------------------------------------------------------------------
        tempdiag = 0.0
        do k=1,nlev
          tempdiag(:,:) = tempdiag(:,:) + Don_conv%qtceme(:,:,k,n)* &
                          pmass(:,:,k)
        end do
         if (id_qtceme_col(n) > 0) then
        used = send_data (id_qtceme_col(n), tempdiag, Time, is, js)
        endif
      end do

!---------------------------------------------------------------------
!    send the 2D convection-related diagnostics to diag_manager_mod.
!---------------------------------------------------------------------

!   pressure at lifting condensation level:
       if (id_plcl_deep > 0) then
      used = send_data (id_plcl_deep, Don_cape%plcl, Time, is, js)
       endif

!   pressure at level of free convection:
       if (id_plfc_deep > 0) then
      used = send_data (id_plfc_deep, Don_cape%plfc, Time, is, js)
       endif

!   pressure at level of zero buoyancy:
       if (id_plzb_deep > 0) then
      used = send_data (id_plzb_deep, Don_cape%plzb, Time, is, js)
       endif

!   convective available potential energy:
      if (id_xcape_deep > 0) then
      used = send_data (id_xcape_deep, Don_cape%xcape_lag, Time, is, js)
       endif

!   convective inhibition:
      if (id_coin_deep > 0) then
      used = send_data (id_coin_deep, Don_cape%coin, Time, is, js)
       endif

!   time tendency of cape:
      if (id_dcape_deep > 0) then
      used = send_data (id_dcape_deep, Don_conv%dcape, Time, is, js)
       endif

!   column integrated water vapor:
      if (id_qint_deep > 0) then
      used = send_data (id_qint_deep, Don_cape%qint_lag, Time, is, js)
       endif

!   fractional area of cumulus ensemble members:
      if (id_a1_deep > 0) then
      used = send_data (id_a1_deep, Don_conv%a1, Time, is, js)
       endif

!   fractional area of largest cumulus ensemble member:
      if (id_amax_deep > 0) then
      used = send_data (id_amax_deep, Don_conv%amax, Time, is, js)
       endif

!   upper limit of fractional area based on moisture constraint:
      if (id_amos_deep > 0) then
      used = send_data (id_amos_deep, Don_conv%amos, Time, is, js)
       endif

!   area-weighted total precipitation:
      if (id_tprea1_deep > 0) then
      used = send_data (id_tprea1_deep, total_precip, Time, is, js)
       endif

!   mesoscale cloud fraction:
       if (id_ampta1_deep > 0) then
      used = send_data (id_ampta1_deep, Don_conv%ampta1, Time, is, js)
       endif

!   accumulated low-level parcel displacement:
       if (id_omint_deep > 0) then
         used = send_data (id_omint_deep, parcel_rise, Time, is, js)
       endif

!   area weighted convective precipitation:
       if (id_rcoa1_deep > 0) then
      used = send_data (id_rcoa1_deep, Don_conv%cell_precip,    &
                        Time, is, js)
       endif

!----------------------------------------------------------------------
!    send diagnostics associated with the monitored output fields.
!----------------------------------------------------------------------
      if (Initialized%monitor_output) then
        do n=1,size(Initialized%Don_monitor,1)
          if (id_extremes(n) > 0) then
            used = send_data (id_extremes(n),   &
                    Initialized%Don_monitor(n)%extrema(is:ie,js:je,:), &
                    Time, is, js,1, mask =    &
                  Initialized%Don_monitor(n)%extrema(is:ie,js:je,:) /= &
                              Initialized%Don_monitor(n)%initial_value )
          endif
          if (id_hits(n) > 0) then
            used = send_data (id_hits(n),  &
                       Initialized%Don_monitor(n)%hits(is:ie,js:je,:), &
                       Time, is, js,1, mask =   &
                 Initialized%Don_monitor(n)%extrema(is:ie,js:je,:) /= &
                              Initialized%Don_monitor(n)%initial_value )
          endif
        end do
      endif

!----------------------------------------------------------------------


end subroutine donner_deep_netcdf


!######################################################################


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      3. ROUTINES CALLED BY DONNER_DEEP_END
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################

subroutine write_restart_nc (ntracers) 

!----------------------------------------------------------------------
!    subroutine write_restart_nc writes a netcdf restart file, either
!    a full file needed when the first step after restart will not be
!    a calculation step, or a reduced file when it is known (expected)
!    that it will be.
!----------------------------------------------------------------------

integer, intent(in) :: ntracers

!----------------------------------------------------------------------
!   intent(in) variables:
!
!     ntracers               number of tracers to be transported by
!                            the donner deep convection parameterization
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      character(len=65)  :: fname = 'RESTART/donner_deep.res.nc'
                                    ! name of restart file to be written
      integer            ::  n      ! do-loop index

!---------------------------------------------------------------------
!    write a message indicating the type of restart to be written and
!    the reason for it.
!---------------------------------------------------------------------
      if (mpp_pe() == mpp_root_pe() ) then
        if (.not. (write_reduced_restart_file) ) then
          call error_mesg ('donner_deep_mod', 'wrote_restart_nc: &
            &Writing FULL netCDF formatted restart file as requested: &
                 &RESTART/donner_deep.res.nc', NOTE)
        else
          if (Initialized%conv_alarm >= Initialized%physics_dt)  then
            call error_mesg ('donner_deep_mod', 'write_restart_nc: &
            &Writing FULL netCDF formatted restart file; it is needed &
             &to allow seamless restart because next step is not a &
             &donner calculation step: RESTART/donner_deep.res.nc', NOTE)
          else
            call error_mesg ('donner_deep_mod', 'write_restart_nc: &
              &Writing REDUCED netCDF formatted restart file as  &
                &requested: RESTART/donner_deep.res.nc', NOTE)
          endif 
        endif 
      endif

!---------------------------------------------------------------------
!    write out the alarm information -- the time remaining before the
!    next calculation and the current donner calculation frequency.
!---------------------------------------------------------------------
      call write_data(fname, 'conv_alarm', Initialized%conv_alarm, &
                                                       no_domain=.true.)
      call write_data(fname, 'donner_deep_freq', donner_deep_freq, &
                                                       no_domain=.true.)

!---------------------------------------------------------------------
!    write out the restart data that is present in a full restart, but
!    which is not needed if donner_deep is calculated on the first
!    step of the restarted job.
!---------------------------------------------------------------------
      if (.not. (write_reduced_restart_file) .or. &
          Initialized%conv_alarm >= Initialized%physics_dt)  then
        call write_data (fname, 'cemetf',            Don_save%cemetf)
        call write_data (fname, 'cememf',            Don_save%cememf)   
        call write_data (fname, 'mass_flux',         Don_save%mass_flux)
        call write_data (fname, 'cell_up_mass_flux',    &
                                              Don_save%cell_up_mass_flux)
        call write_data (fname, 'det_mass_flux',        &
                                                  Don_save%det_mass_flux)
        call write_data (fname, 'dql_strat',         Don_save%dql_strat)
        call write_data (fname, 'dqi_strat',         Don_save%dqi_strat)
        call write_data (fname, 'dqa_strat',         Don_save%dqa_strat)
        call write_data (fname, 'tprea1',            Don_save%tprea1)   
        call write_data (fname, 'humidity_area',        &
                                                  Don_save%humidity_area)
        call write_data (fname, 'humidity_factor',       &
                                               Don_save%humidity_factor)
        if (Initialized%do_donner_tracer) then
          do n=1,ntracers
            call write_data (fname,   &
                  'tracer_tends_'// trim(Don_save%tracername(n)),  &
                                          Don_save%tracer_tends(:,:,:,n))
          end do
        endif
      endif

!----------------------------------------------------------------------
!    write out the restart data which is always needed, regardless of
!    when the first donner calculation step is after restart.
!----------------------------------------------------------------------
      call write_data (fname, 'parcel_disp', Don_save%parcel_disp)
      call write_data (fname, 'lag_temp',    Don_save%lag_temp)     
      call write_data (fname, 'lag_vapor',   Don_save%lag_vapor)     
      call write_data (fname, 'lag_press',   Don_save%lag_press)     

!----------------------------------------------------------------------



end subroutine write_restart_nc




!#####################################################################

subroutine write_restart (ntracers)          

!--------------------------------------------------------------------
!    subroutine write_restart is a template to be used if a native mode
!    restart file MUST be generated. currently, if a native mode file is
!    requested, a netcdf file will be witten instead, and an informative
!    message provided.
!--------------------------------------------------------------------
 
integer, intent(in) :: ntracers

!----------------------------------------------------------------------
!   intent(in) variables:
!
!     ntracers               number of tracers to be transported by
!                            the donner deep convection parameterization
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

!     integer :: unit          ! unit number for restart file
!     integer :: n             ! do-loop index

!-------------------------------------------------------------------
!    currently code is provided only for writing netcdf restart files.
!    if a non-netcdf restart file has been requested, this routine will 
!    issue a message, and then call the routine to write the netcdf file.
!    if the user is insistent on a native mode restart file, the code to
!    read and write such files (subroutines write_restart and 
!    read_restart_file) must be updated to be compatible with  the cur-
!    rent versions of write_restart_nc and read_restart_nc, and the 
!    code immediately below eliminated. the commented code below repres-
!    ents a starting point for the write_restart routine; it is not 
!    kept up-to-date as far as the variables which must be written.
!-------------------------------------------------------------------
      call error_mesg ('donner_deep_mod', 'write_restart: &
          &writing a netcdf restart despite request for native &
           &format (not currently supported); if you must have native &
           &mode, then you must update the source code and remove &
                                               &this if loop.', NOTE)
      call write_restart_nc (ntracers) 

!-------------------------------------------------------------------
!    open unit for restart file.
!-------------------------------------------------------------------
!      unit = open_restart_file ('RESTART/donner_deep.res', 'write')

!-------------------------------------------------------------------
!    file writing is currently single-threaded. write out restart
!    version, time remaining until next call to donner_deep_mod and
!    the frequency of calculating donner_deep convection.
!-------------------------------------------------------------------
!     if (mpp_pe() == mpp_root_pe()) then
!       write (unit) restart_versions(size(restart_versions(:)))
!       write (unit) Initialized%conv_alarm, donner_deep_freq
!     endif

!-------------------------------------------------------------------
!    write out the donner_deep restart variables.
!    cemetf    - heating rate due to donner_deep
!    cememf    - moistening rate due to donner_deep
!    xcape_lag - cape value which will be used on next step in
!                calculation od dcape/dt
!-------------------------------------------------------------------
!     call write_data (unit, Don_save%cemetf)
!     call write_data (unit, Don_save%cememf)
      
!--------------------------------------------------------------------
!    the following variables are needed when a prognostic cloud scheme
!    is being used. they are always present in the restart file, having
!    been initialized to zero, if prognostic clouds are not active.
!--------------------------------------------------------------------
!     call write_data (unit, Don_save%mass_flux)
!     call write_data (unit, Don_save%dql_strat )
!     call write_data (unit, Don_save%dqi_strat )
!     call write_data (unit, Don_save%dqa_strat )

!----------------------------------------------------------------------
!    
!-------------------------------------------------------------------
!    write out more donner_deep restart variables.
!    qint_lag   - column integrated water vapor mixing ratio
!    parcel_disp  - time-integrated low-level vertical displacement
!    tprea1     - precipitation due to donner_deep_mod
!----------------------------------------------------------------------
!     call write_data (unit, Don_save%parcel_disp)
!     call write_data (unit, Don_save%tprea1)
!     call write_data (unit, Don_save%lag_temp)
!     call write_data (unit, Don_save%lag_vapor)
!     call write_data (unit, Don_save%lag_press)
!     call write_data (unit, Don_save%humidity_area)
!     call write_data (unit, Don_save%humidity_ratio)

!---------------------------------------------------------------------
!    write out the number of tracers that are being transported by
!    donner_deep_mod.
!---------------------------------------------------------------------
!     if (mpp_pe() == mpp_root_pe()) then
!       write (unit) ntracers
!     endif

!----------------------------------------------------------------------
!    if tracers are being transported, write out their names and 
!    current time tendencies.
!----------------------------------------------------------------------
!     if (Initialized%do_donner_tracer) then
!       do n=1,ntracers
!         if (mpp_pe() == mpp_root_pe()) then
!           write (unit) Don_save%tracername(n)         
!         endif
!         call write_data(unit, Don_save%tracer_tends(:,:,:,n))
!       end do
!     endif

!-------------------------------------------------------------------
!    close restart file unit.
!------------------------------------------------------------------
!     call close_file (unit)

!---------------------------------------------------------------------


end subroutine write_restart




!######################################################################

subroutine deallocate_variables 

!---------------------------------------------------------------------
!    subroutine deallocate_variables deallocates the space used by the
!    module variables.
!---------------------------------------------------------------------
 
      deallocate ( Don_save%cemetf              )
      deallocate ( Don_save%lag_temp            )
      deallocate ( Don_save%lag_vapor           )
      deallocate ( Don_save%lag_press           )
      deallocate ( Don_save%cememf              )
      deallocate ( Don_save%mass_flux           )
      deallocate ( Don_save%cell_up_mass_flux   )
      deallocate ( Don_save%det_mass_flux       )
      deallocate ( Don_save%dql_strat           )
      deallocate ( Don_save%dqi_strat           )
      deallocate ( Don_save%dqa_strat           )
      deallocate ( Don_save%humidity_area       )
      deallocate ( Don_save%humidity_factor     )
      deallocate ( Don_save%tracer_tends        )
      deallocate ( Don_save%parcel_disp         )
      deallocate ( Don_save%tprea1              )
      deallocate ( Don_save%tracername          )
      deallocate ( Don_save%tracer_units        )

      deallocate (Param%arat)
      deallocate (Param%ensemble_entrain_factors_gate)
      deallocate (Param%ensemble_entrain_factors_kep )
      deallocate (Param%dgeice)
      deallocate (Param%relht )


      deallocate (Col_diag%i_dc    )
      deallocate (Col_diag%j_dc    ) 
      deallocate (Col_diag%unit_dc )
      deallocate (Col_diag%igl_dc  )
      deallocate (Col_diag%jgl_dc  )

      if (allocated(col_diag_unit)) then
        deallocate (col_diag_unit  )
        deallocate (col_diag_lon   )
        deallocate (col_diag_lat   )
        deallocate (col_diag_i     )
        deallocate (col_diag_j     )
      endif
     
      if (allocated (id_qtren1)) then
      deallocate (id_qtren1)
      deallocate (id_qtmes1)
      deallocate (id_wtp1  )
      deallocate (id_qtceme)
      deallocate (id_total_wet_dep)
      deallocate (id_meso_wet_dep)
      deallocate (id_cell_wet_dep)
      deallocate (id_qtren1_col)
      deallocate (id_qtmes1_col)
      deallocate (id_wtp1_col  )
      deallocate (id_qtceme_col)
      deallocate (id_total_wet_dep_col)
      deallocate (id_meso_wet_dep_col)
      deallocate (id_cell_wet_dep_col)
      endif

      call sd_end_k(sd);
      call ac_end_k(ac);
      call cp_end_k(cp);
      call ct_end_k(ct);
      call exn_end_k
      call findt_end_k

      if (Initialized%monitor_output) then
        deallocate (id_extremes)
        deallocate (id_hits)
      endif

!----------------------------------------------------------------------


end subroutine deallocate_variables 




!######################################################################



                     end module donner_deep_mod

