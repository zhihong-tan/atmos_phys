!  $Id: donner_nml.h,v 16.0 2008/07/30 22:06:58 fms Exp $

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

logical             :: do_ensemble_diagnostics = .false.
                             ! include netcdf diagnostics for selected
                             ! variables for each ensemble member ?

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
                            do_ensemble_diagnostics, &

! not contained in donner_nml_type variable:
                            do_netcdf_restart, &
                            write_reduced_restart_file, &
                            diagnostics_pressure_cutoff, &
                            diagnostics_start_time, &
                            num_diag_pts_ij, num_diag_pts_latlon, &
                            i_coords_gl, j_coords_gl, &
                            lat_coords_gl, lon_coords_gl

