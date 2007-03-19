!FDOC_TAG_GFDL
module strat_cloud_mod
  ! <CONTACT EMAIL="Stephen.Klein@noaa.gov">
  !   Stephen Klein
  ! </CONTACT>
  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
  ! <OVERVIEW>
  !   Code to compute time tendencies of stratiform clouds and diagnoses
  !   rain and snow flux with prognostic scheme.
  !   
  ! </OVERVIEW>
  ! <DESCRIPTION>
  !
  !
  !       The prognostic scheme returns the time tendencies of liquid,
  !       ice, and saturated volume fraction that are suspended in 
  !       stratiform clouds.  The scheme diagnoses the fluxes of rain
  !       and snow in saturated and unsaturated areas.
  !
  !       The prognostic cloud scheme is responsible for determing
  !       cloud volume fractions, condensed water tendencies, and
  !       the stratiform precipitation rate.  It includes processes
  !       for evaporation, condensation, deposition, and sublimation
  !       of cloud water, conversion of cloud water to precipitation,
  !       evaporation of falling precipitation, the bergeron-findeisan 
  !       process, freezing of cloud liquid, accretion of cloud water 
  !       by precipitation, and melting of falling precipitation.
  !
  !       This scheme is based on the experience the author had 
  !       at the ECMWF in 1997. The saturated volume fraction formalism 
  !       and type of solution follow directly from the scheme of Tiedtke
  !       (1993): Monthly Weather Review, Volume 121, pages 3040-3061.
  !       The form of most of the microphysics follows Rotstayn , 1997:
  !       Quart. J. Roy. Met. Soc. vol 123, pages 1227-1282. The partial
  !       precipitation area formulism follows Jakob and Klein, 2000:
  !       Quart. J. Roy. Met. Soc. vol 126, pages 2525-2544. 
  !
  !       The statistical cloud scheme treatment, which is used as
  !       a replacement for the Tiedtke cloud fraction scheme, is based
  !       on a number of publications: Tompkins, A., 2002: J. Atmos. 
  !       Sci., 59, 1917-1942, Klein et al., 2005: J. Geophys. Res., 
  !       110, D15S06, doi:10.1029/2004JD005017. 
  ! </DESCRIPTION>
  !

! <DATASET NAME="strat_cloud.res">
!   native format of the restart file
! </DATASET>
! <DATASET NAME="strat_cloud.res.nc">
!   netcdf format of the restart file
! </DATASET>


! <INFO>
!   <REFERENCE>           
!The saturation volume fraction formalism comes from:
!
!Tiedtke, M., 1993: Representation of clouds in large-scale models. Mon. Wea. Rev., 121, 3040-3061.
!
! </REFERENCE>
!   <REFERENCE>           
!The form of most of the microphysics follows:
!
!Rotstayn, L., 1997: A physically based scheme for the treatment of stratiform clouds and precipitation in large-scale models. I: Description and evaluation of microphysical processes. Quart. J. Roy. Met. Soc. 123, 1227-1282. 
! </REFERENCE>
!   <COMPILER NAME="">     </COMPILER>
!   <PRECOMP FLAG="">      </PRECOMP>
!   <LOADER FLAG="">       </LOADER>
!   <TESTPROGRAM NAME="">  </TESTPROGRAM>
!   <BUG>                  </BUG>
!   <NOTE> 
!1. qmin should be chosen such that the range of {qmin, max(qa,ql,qi)} is resolved by the precision of the numbers used. (default = 1.E-10)
!   </NOTE>

!   <NOTE> 
!2. Dmin will be MACHINE DEPENDENT and occur when
!   </NOTE>

!   <NOTE> 
!a. 1. -exp(-Dmin) = 0. instead of Dmin in the limit of very small Dmin
!   </NOTE>

!AND

!   <NOTE> 
!b. 1. - exp(-D) < D for all D > Dmin
!   </NOTE>
!   <FUTURE>               </FUTURE>

! </INFO>

  use  sat_vapor_pres_mod, only :  lookup_es,lookup_des
  use             fms_mod, only :  file_exist, open_namelist_file,  &
       error_mesg, FATAL, NOTE,         &
       mpp_pe, mpp_root_pe, close_file, &
       read_data, write_data,           &
       check_nml_error, &
       write_version_number, stdlog, &
       open_restart_file, open_ieee32_file, &
       mpp_error
  use  fms_io_mod,         only :  get_restart_io_mode
  use  constants_mod,      only :  rdgas,rvgas,hlv,hlf,hls,      &
       cp_air,grav,tfreeze,dens_h2o
  use  cloud_rad_mod,      only :  cloud_rad_init
  use  diag_manager_mod,   only :  register_diag_field, send_data
  use  time_manager_mod,   only :  time_type, get_date
  use cloud_generator_mod, only :  do_cloud_generator,           &
       compute_overlap_weighting
  use beta_dist_mod,        only: beta_dist_init, beta_dist_end, &
                                  incomplete_beta
  use  rad_utilities_mod,  only : aerosol_type
  use  aer_ccn_act_mod,    only : aer_ccn_act

  implicit none

  public  strat_cloud_init,    &
       strat_cloud,         &
       strat_cloud_end,     &
       strat_cloud_sum,     &
       strat_cloud_avg,     &
       do_strat_cloud,      &
       strat_cloud_on

  !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !              GLOBAL STORAGE VARIABLES
  !
  !     radturbten    The sum of radiation and turbulent tendencies
  !                   for each grid box. (K/sec)
  !


  !
  !     ------ data for cloud averaging code ------
  !

  real,    allocatable, dimension (:,:,:) :: qlsum, qisum, cfsum
  integer, allocatable, dimension (:,:)   :: nsum
  !
  !     ------ constants used by the scheme -------
  !

  real, parameter :: d608 = (rvgas-rdgas) / rdgas
  real, parameter :: d622 = rdgas / rvgas
  real, parameter :: d378 = 1. - d622
  logical         :: do_netcdf_restart = .true.
  !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !       DECLARE CONSTANTS AND SET DEFAULT VALUES FOR PARAMETERS OF 
  !       THE SCHEME
  !
  !
  !
  !                  PHYSICAL CONSTANTS USED IN THE SCHEME
  !
  !
  !         constant              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !           grav       gravitational acceleration      m/(s*s)
  !
  !           hlv        latent heat of vaporization     J/kg condensate
  !
  !           hlf        latent heat of fusion           J/kg condensate
  !
  !           hls        latent heat of sublimation      J/kg condensate
  !
  !           rdgas      gas constant of dry air         J/kg air/K
  !
  !           rvgas      gas constant of water vapor     J/kg air/K
  !
  !           cp_air     specific heat of air at         J/kg air/K
  !                      constant pressure
  !
  !           d622       rdgas divided by rvgas          dimensionless
  !
  !           d378       One minus d622                  dimensionless
  !
  !           tfreeze    Triple point of water           K
  !
  !           dens_h2o   Density of pure liquid          kg/(m*m*m)
  !
  !
  !
  !
  !                          PARAMETERS OF THE SCHEME 
  !
  !
  !         parameter              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !            U00       threshold relative humidity     fraction
  !                      for cloud formation 
  !
  !       u00_profile    should low-level u00 ECMWF profile be applied?
  !
  !       rthresh        liquid cloud drop radius        microns
  !                      threshold for autoconversion
  !
  !       use_kk_auto    Should one use the Khairoutdinv and Kogan (2000)
  !                      autoconversion formula?
  !
  !       N_land         fixed number of cloud drops     1/(m*m*m)
  !                      per unit volume in liquid 
  !                      clouds on land
  !
  !       N_ocean        fixed number of cloud drops     1/(m*m*m)
  !                      per unit volume in liquid 
  !                      clouds over ocean
  !
  !       rho_ice        mass density of ice crystals    kg/(m*m*m)
  !
  !       ELI            collection efficiency of        dimensionless
  !                      cloud liquid by falling ice
  !
  !       U_evap         critical relative humidity      fraction
  !                      above which rain does not
  !                      evaporate
  !
  !       eros_scale     normal erosion rate constant    1/sec
  !                      cloud destruction
  !
  !       eros_choice    should enhanced erosion in      logical
  !                      turbulent conditions be done?
  !
  !       eros_scale_c   erosion rate constant for       1/sec
  !                      convective conditions
  !
  !       eros_scale_t   erosion rate constant for       1/sec
  !                      cloud destruction for
  !                      turbulent conditions
  !
  !       mc_thresh      Convective mass-flux            kg/m2/sec
  !                      threshold for enhanced
  !                      erosion to turn on.
  !
  !       diff_thresh    Diffusion coefficient           m2/s
  !                      threshold for enhanced
  !                      erosion to turn on.
  !
  !       super_choice   Should should excess vapor      logical
  !                      in supersaturated conditions    
  !                      be put into cloud water (true) 
  !                      or precipitation fluxes (false)
  !                 
  !       tracer_advec   Are cloud liquid,ice and        logical
  !                      fraction advected by the
  !                      grid resolved motion?
  !
  !       qmin           minimum permissible value of    kg condensate/
  !                      cloud liquid, cloud ice,        kg air
  !                      saturated volume fraction,
  !                      or rain and snow areas
  !
  !                      NOTE: qmin should be chosen
  !                      such that the range of
  !                      {qmin, max(qa,ql,qi)} is
  !                      resolved by the precision 
  !                      of the numbers used.
  !
  !       Dmin           minimum permissible             dimensionless
  !                      dissipation in analytic 
  !                      integration of qa, ql, qi
  !                      equations. This constant
  !                      only affects the method by
  !                      which the prognostic equations
  !                      are integrated.
  !
  !                      NOTE: Dmin will be MACHINE 
  !                      DEPENDENT and occur when 
  !                      a. 1. -exp(-Dmin)  = 0. 
  !                         instead of Dmin in the 
  !                         limit of very small Dmin
  !
  !                      AND 
  !
  !                      b. 1. - exp(-D) < D for
  !                         all D > Dmin
  !
  !       do_average     Average stratiform cloud properties
  !                      before computing clouds used by radiation?
  !
  !       strat_cloud_on Is the stratiform cloud scheme
  !                      operating? 
  !
  !       do_budget_diag Are any of the budget diagnostics
  !                      requested from this run?
  !                          
  !       num_strat_pts  number of grid points where 
  !                      instantaneous output will be 
  !                      saved to file strat.data
  !
  !                      num_strat_pts <= max_strat_pts
  !       
  !       max_strat_pts  maximum number of strat pts
  !                      for instantaneous output
  !
  !       strat_pts      "num_strat_pts" pairs of grid
  !                      indices, e.g., the global 
  !                      indices for i,j.
  !
  !       overlap        value of the overlap parameter
  !                      from cloud rad
  !                      overlap = 1 is maximum-random
  !                      overlap = 2 is random
  !
  !     do_old_snowmelt  Should the cloud scheme be run with
  !                      the snowmelt bug?
  !
  !       do_liq_num     Should the prognostic droplet number
  !                      concentration be used?
  ! 
  !     do_pdf_clouds    Should the statistical cloud scheme
  !                      be used?  
  ! 
  !         betaP        the p parameter in the beta distribution
  !
  !      qthalfwidth     The fraction of qtbar (mean total water in the
  !                      grid box) that the maximum and minimum of the 
  !                      distribution differ from qtbar. That is, total 
  !                      water at the sub-grid scale may take on values 
  !                      anywhere between (1.-qthalfwidth)*qtbar and 
  !                      (1.+qthalfwidth)*qtbar
  !
  !      nsublevels      This is the number of sub-levels to be used
  !                      for sub-grid scale vertical structure to
  !                      clouds. If equal to 1, then no vertical
  !                      sub-grid scale structure is calculated.
  !
  !      kmap, kord      Quantities related to the PPM vertical inter-
  !                      polation calculation.
  !
  !                      THIS IS ONLY USED WITH DIAGNOSTIC VARIANCE
  
  real              :: U00            =  0.80
  logical           :: u00_profile    =  .false.
  real              :: rthresh        =  10.
  logical           :: use_kk_auto    =  .false.
  real              :: N_land         =  250.E+06
  real              :: N_ocean        =  100.E+06
  real,   parameter :: rho_ice        =  100.
  real,   parameter :: ELI            =  0.7
  real              :: U_evap         =  1.0
  real              :: eros_scale     =  1.E-06
  logical           :: eros_choice    =  .false.
  real              :: eros_scale_c   =  8.E-06
  real              :: eros_scale_t   =  5.E-05
  real              :: mc_thresh      =  0.001
  real              :: diff_thresh    =  1.0
  logical           :: super_choice   =  .false.
  logical           :: tracer_advec   =  .false.
  real              :: qmin           =  1.E-10
  real              :: Dmin           =  1.E-08
  logical           :: do_average     =  .false.
  logical           :: strat_cloud_on =  .false.
  logical           :: do_budget_diag =  .false.
  integer,parameter :: max_strat_pts  =  5
  integer           :: num_strat_pts  =  0
  integer,dimension(2,max_strat_pts) :: strat_pts = 0
  integer           :: overlap        =  2
  real              :: efact          = 0.0
  logical           :: do_old_snowmelt= .false.
  logical           :: do_liq_num   = .false.
  real              :: N_min          = 1.E6
  logical           :: do_pdf_clouds  = .false.
  real              :: qthalfwidth    = 0.1
  integer           :: nsublevels     = 1
  integer           :: kmap           = 1
  integer           :: kord           = 7
  integer           :: betaP          = 5

  !
  !-----------------------------------------------------------------------
  !-------------------- diagnostics fields -------------------------------

  integer :: id_droplets 
  integer :: id_aliq,         id_aice,            id_aall,       &
       id_rvolume,      id_autocv,          id_vfall 
  integer :: id_qldt_cond,    id_qldt_eros,       id_qldt_fill,  &
       id_qldt_accr,    id_qldt_evap,       id_qldt_freez, &
       id_qldt_berg,    id_qldt_destr,      id_qldt_rime,  &
       id_qldt_auto,    id_qndt_cond,       id_qndt_evap

  integer :: id_rain_clr,     id_rain_cld,        id_a_rain_clr, &
       id_a_rain_cld,   id_rain_evap,       id_liq_adj
  integer :: id_qidt_fall,    id_qidt_fill,       id_qidt_melt,  &
       id_qidt_dep,     id_qidt_subl,       id_qidt_eros,  &
       id_qidt_destr
  integer :: id_snow_clr,     id_snow_cld,        id_a_snow_clr, &
       id_a_snow_cld,   id_snow_subl,       id_snow_melt,  &
       id_ice_adj
  integer :: id_ql_eros_col,  id_ql_cond_col,   id_ql_evap_col,  &
       id_ql_accr_col,  id_ql_auto_col,   id_ql_fill_col,  &
       id_ql_berg_col,  id_ql_destr_col,  id_ql_rime_col,  &       
       id_ql_freez_col    
  integer :: id_rain_evap_col,id_liq_adj_col
  integer :: id_qi_fall_col,  id_qi_fill_col,   id_qi_subl_col,  &
       id_qi_melt_col,  id_qi_destr_col,  id_qi_eros_col,  &
       id_qi_dep_col
  integer :: id_snow_subl_col,id_snow_melt_col, id_ice_adj_col
  integer :: id_qadt_lsform,  id_qadt_eros,     id_qadt_fill,    &
       id_qadt_rhred,   id_qadt_destr,   &
       id_qadt_lsdiss,  id_qadt_super   
  integer :: id_qa_lsform_col,id_qa_eros_col,   id_qa_fill_col,  &
       id_qa_rhred_col, id_qa_destr_col, &
       id_qa_lsdiss_col,id_qa_super_col
  integer :: id_a_precip_cld, id_a_precip_clr

  character(len=5) :: mod_name = 'strat'
  real :: missing_value = -999.

  !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !       CREATE NAMELIST
  !

! <NAMELIST NAME="strat_cloud_nml">
!  <DATA NAME="do_netcdf_restart" UNITS="fraction" TYPE="logical" DIM="" DEFAULT="">
!   netcdf/native restart format
!  </DATA>
!  <DATA NAME="U00" UNITS="" TYPE="real" DIM="" DEFAULT="">
!     Threshold relative humidity for cloud formation by large-scale condensation. (default = 0.80) 
!  </DATA>
!  <DATA NAME="u00_profile" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! Should low-level u00 ECMWF profile be applied? (default = .false.) 
!  </DATA>
!  <DATA NAME="rthresh" UNITS="microns" TYPE="real" DIM="" DEFAULT="">
!  Liquid cloud drop radius threshold for autoconversion. (default = 10.)
!  </DATA>
!  <DATA NAME="use_kk_auto" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! Should the Khairoutdinov and Kogan (2000) autoconversion be used? ( default = .false.)
!  </DATA>
!  <DATA NAME="N_land" UNITS="1/(m*m*m)" TYPE="real" DIM="" DEFAULT="">
! Fixed number of cloud drops per unit volume in liquid clouds on land. ( default = 250.E+06)
!  </DATA>
!  <DATA NAME="N_ocean" UNITS="1/(m*m*m)" TYPE="real" DIM="" DEFAULT="">
!  Fixed number of cloud drops per unit volume in liquid clouds over ocean. ( default = 100.E+06)
!  </DATA>
!  <DATA NAME="U_evap" UNITS="fraction" TYPE="real" DIM="" DEFAULT="">
!    Critical relative humidity above which rain does not evaporate. (default = 1.0) 
!  </DATA>
!  <DATA NAME="eros_scale" UNITS="1/sec" TYPE="real" DIM="" DEFAULT="">
! Normal erosion rate constant cloud destruction (default = 1.E-06) 
!  </DATA>
!  <DATA NAME="eros_choice" UNITS="" TYPE="real" DIM="" DEFAULT="">
! Should enhanced erosion in turbulent conditions be done? (default = .false.)
!  </DATA>
!  <DATA NAME="eros_scale_c" UNITS="1/sec" TYPE="real" DIM="" DEFAULT="">
!  Erosion rate constant for convective conditions. (default = 8.E-05)
!  </DATA>
!  <DATA NAME="eros_scale_t" UNITS="1/sec" TYPE="real" DIM="" DEFAULT="">
! Erosion rate constant for cloud destruction for turbulent conditions. (default = 5.E-05)
!  </DATA>
!  <DATA NAME="mc_thresh" UNITS="kg/m2/sec" TYPE="real" DIM="" DEFAULT="">
!  Convective mass-flux threshold for enhanced erosion to turn on. (default = 0.001) 
!  </DATA>
!  <DATA NAME="diff_thresh" UNITS="m2/s" TYPE="real" DIM="" DEFAULT="">
!  Diffusion coefficient threshold for enhanced erosion to turn on. (default = 1.0) 
!  </DATA>
!  <DATA NAME="super_choice" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! Should should excess vapor in supersaturated conditions be put into cloud water (true) or precipitation fluxes (false)? (default = .false.) 
!  </DATA>
!  <DATA NAME="tracer_advec" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! Are cloud liquid,ice and fraction advected by the grid resolved motion? (default = .false.) 
!  </DATA>
!  <DATA NAME="qmin" UNITS="kg condensate/kg air" TYPE="real" DIM="" DEFAULT="">
!  Minimum permissible value of cloud liquid, cloud ice, saturated volume fraction, or rain and snow areas.

! NOTE: qmin should be chosen such that the range of {qmin, max(qa,ql,qi)} is resolved by the precision of the numbers used. (default = 1.E-10) 
!  </DATA>
!  <DATA NAME="Dmin" UNITS="Dimensionless" TYPE="real" DIM="" DEFAULT="">
! Minimum permissible dissipation in analytic integration of qa, ql, qi equations. This constant only affects the method by which the prognostic equations are integrated.

!NOTE: Dmin will be MACHINE DEPENDENT and occur when

!a. 1. -exp(-Dmin) = 0. instead of Dmin in the limit of very small Dmin

!AND

!b. 1. - exp(-D) < D for all D > Dmin

!(default = 1.E-08) 
!  </DATA>
!  <DATA NAME="num_strat_pts" UNITS="" TYPE="integer" DIM="" DEFAULT="">
! Number of grid points where instantaneous output will be saved to file strat.data

!num_strat_pts <= max_strat_pts

!(default = 0)
!  </DATA>
!  <DATA NAME="strat_pts" UNITS="" TYPE="integer" DIM="" DEFAULT="">
!num_strat_pts" pairs of grid indices, e.g., the global indices for i,j. (default = 0) 
!  </DATA>
!  <DATA NAME="efact" UNITS="" TYPE="real" DIM="" DEFAULT="">
! (default = 0.0) 
!  </DATA>
!  <DATA NAME="do_old_snowmelt" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! Should the old version of snow melting, which has a bug,be run? (default = .false.) 
!  </DATA>
!  <DATA NAME="do_liq_num" UNITS="" TYPE="logical" DIM="" DEFAULT="">
!  Should the prognostic droplet number code be run? (default = .false.) 
!  </DATA>
!  <DATA NAME="do_pdf_clouds" UNITS="" TYPE="logical" DIM="" DEFAULT="">
! Should the statistical cloud scheme be run? (default = .false.) 
!  </DATA>
!  <DATA NAME="qthalfwidth" UNITS="none" TYPE="real" DIM="" DEFAULT="">
! Half-width to the qt PDF - used only if do_pdf_clouds is true and diagnostic variance(default = 0.1) 
!  </DATA>
!  <DATA NAME="nsublevels" UNITS="none" TYPE="integer" DIM="" DEFAULT="">
! Number of sublevels to vertical sub-grid cloud structure - used only if do_pdf_cloud is true (default = 1) 
!  </DATA>
!  <DATA NAME="kmap" UNITS="none" TYPE="integer" DIM="" DEFAULT="">
! PPM partial remap integer - used only if do_pdf_cloud is true and if vertical subgrid structure is used(default = 1) 
!  </DATA>
!  <DATA NAME="kord" UNITS="none" TYPE="integer" DIM="" DEFAULT="">
! PPM method number - used only if do_pdf_cloud is true and if vertical subgrid structure is used (default = 7) 
!  </DATA>
!  <DATA NAME="betaP" UNITS="none" TYPE="integer" DIM="" DEFAULT="">
! p-parameter to the beta distribution - used only if do_pdf_clouds is true (default = 5) 
!  </DATA>
! </NAMELIST>

  NAMELIST /strat_cloud_nml/ do_netcdf_restart,   &
       U00,u00_profile,rthresh,use_kk_auto,N_land,&
       N_ocean,U_evap,eros_scale,eros_choice,   &
       eros_scale_c,eros_scale_t,mc_thresh,     &
       diff_thresh,super_choice,tracer_advec,   &
       qmin,Dmin,num_strat_pts,strat_pts,efact, &
       do_old_snowmelt, do_pdf_clouds, betaP,   &
       qthalfwidth,nsublevels,kmap,kord, do_liq_num
       
  !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !       DECLARE VERSION NUMBER OF SCHEME
  !

  Character(len=128) :: Version = '$Id: strat_cloud.F90,v 13.1 2006/03/30 17:24:03 fms Exp $'
  Character(len=128) :: Tagname = '$Name: nalanda $'
  logical            :: module_is_initialized = .false.
  integer, dimension(1) :: restart_versions = (/ 1 /)
  !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !       The module contains the following subroutines:
  !
  !
  !       strat_cloud_init    read namelist file, open logfile, initialize
  !                     constants and fields, read restart file
  !
  !       diag_field_init 
  !                     initializes diagnostic fields
  !
  !       strat_cloud   calculations of the cloud scheme are performed 
  !                     here
  !
  !       add_strat_tend  
  !                     Adds a field to radturbten. This subroutine is 
  !                     needed because of the method to calculate the 
  !                     radiative and turbulent tendencies.
  !
  !       subtract_strat_tend 
  !                     Subtracts a field from radturbten.
  !
  !       strat_cloud_end     writes out restart data to a restart file.
  !
  !       strat_cloud_sum
  !                     sum cloud scheme variables
  !
  !       strat_cloud_avg
  !                     return average of summed cloud scheme variables
  !
  !       do_strat_cloud
  !                     logical flag, is the scheme on?
  !


CONTAINS



  !#######################################################################
  !#######################################################################


  ! <SUBROUTINE NAME="strat_cloud_init">
  !  <OVERVIEW>
  !   
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !       Initializes strat_cloud.  Reads namelist, calls cloud_rad_init,
  !       reads restart (if present), initializes netcdf output.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call strat_cloud_init(axes,Time,idim,jdim,kdim)
  !                
  !  </TEMPLATE>
  !  <IN NAME="axes" TYPE="integer">
  !       Axes integer vector used for netcdf initialization.
  !  </IN>
  !  <IN NAME="Time" TYPE="time_type">
  !       Time type variable used for netcdf.
  !  </IN>
  !  <IN NAME="idim" TYPE="integer">
  !       Size of first array (usually longitude) dimension.
  !  </IN>
  !  <IN NAME="jdim" TYPE="integer">
  !       Size of second array (usually latitude) dimension.
  !  </IN>
  !  <IN NAME="kdim" TYPE="integer">
  !       Size of vertical array (usually height) dimension.
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine strat_cloud_init(axes,Time,idim,jdim,kdim)


    !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !       This subroutine reads the namelist file, opens a logfile, 
    !       and initializes the physical constants of the routine.
    !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !
    !
    !       VARIABLES
    !
    !
    !       ------
    !       INPUT:
    !       ------
    !
    !         variable              definition                  unit
    !       ------------   -----------------------------   ---------------
    !
    !       axes           integers corresponding to the
    !                      x,y,z,z_half axes types
    !
    !       Time           time type variable
    !
    !       idim,jdim      number of points in first 
    !                      and second dimensions
    !
    !       kdim           number of points in vertical
    !                      dimension
    !
    !
    !       -------------------
    !       INTERNAL VARIABLES:
    !       -------------------
    !
    !         variable              definition                  unit
    !       ------------   -----------------------------   ---------------
    !
    !       unit           unit number for namelist and
    !                      restart file
    !
    !       io             internal variable for reading
    !                      of namelist file
    !
    !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !


    !        
    !       user Interface variables
    !       ------------------------
    !

    integer, intent (in)                   :: idim,jdim,kdim,axes(4)
    type(time_type), intent(in)            :: Time

    !
    !       Internal variables
    !       ------------------
    !

    integer                                :: unit,io,ierr
    integer                                :: vers, vers2
    character(len=4)                       :: chvers
    character(len=5)                       :: chio
    character(len=64)                      :: fname='INPUT/strat_cloud.res.nc'
    integer, dimension(4)                  :: siz
    !-----------------------------------------------------------------------
    !       
    !
    if(module_is_initialized) then
      return
    else
      module_is_initialized = .true.
    endif
    !-----------------------------------------------------------------------
    !
    !       Namelist functions

    !       ----- read namelist -----

    if ( file_exist('input.nml')) then
       unit = open_namelist_file ()
       ierr=1; do while (ierr /= 0)
       read  (unit, nml=strat_cloud_nml, iostat=io, end=10)
       ierr = check_nml_error(io,'strat_cloud_nml')
    enddo
10  call close_file (unit)
 endif
 call get_restart_io_mode(do_netcdf_restart)

 !write namelist variables to logfile
 if ( mpp_pe() == mpp_root_pe() ) then
    call write_version_number(Version, Tagname)
    write (stdlog(),nml=strat_cloud_nml)
 endif

!-----------------------------------------------------------------------
!    qthalfwidth must be greater than 0.001
!-----------------------------------------------------------------------
        if (qthalfwidth .lt. 1.e-03) then
          call error_mesg ( 'strat_cloud_mod', &
           'qthalfwidth must be greater than 0.001', FATAL)
        endif

!-----------------------------------------------------------------------
!    nsublevels must be greater than 0
!-----------------------------------------------------------------------
        if (nsublevels .lt. 1) then
          call error_mesg ( 'strat_cloud_mod', &
           'nsublevels must be greater than 0', FATAL)
        endif

 !-----------------------------------------------------------------------
 !
 !       initialize qmin, N_land, N_ocean and selected physical constants
 !       in cloud_rad_mod

 if (do_liq_num) then 
   call cloud_rad_init(axes,Time,qmin_in=qmin,N_land_in=N_land,&
        N_ocean_in=N_ocean,prog_droplet_in=do_liq_num,overlap_out=overlap)
 else
   call cloud_rad_init(axes,Time,qmin_in=qmin,N_land_in=N_land,&
        N_ocean_in=N_ocean,overlap_out=overlap)
 endif
 !-----------------------------------------------------------------------
 !
 !       initialize strat_cloud_on to true

     strat_cloud_on = .TRUE.

 !-----------------------------------------------------------------------
 !
 !       initialize the beta distribution module if needed
 
 if (do_pdf_clouds) call beta_dist_init
 
 !-----------------------------------------------------------------------
 !
 !       Read Restart file


 !set up stratiform cloud storage
 !           PRINT *, idim, jdim
 allocate(nsum(idim, jdim),      &
      qlsum(idim,jdim,kdim), &
      qisum(idim,jdim,kdim), &
      cfsum(idim,jdim,kdim)  )

 !see if restart file exists
 if (file_exist('INPUT/strat_cloud.res.nc') ) then
    if(mpp_pe() == mpp_root_pe() ) call mpp_error ('strat_cloud_mod', &
         'Reading netCDF formatted restart file: INPUT/strat_cloud.res.nc', NOTE)
    call read_data(fname, 'vers', vers, no_domain=.true.)
    call read_data(fname, 'nsum',  nsum)
    call read_data(fname, 'qlsum', qlsum)
    call read_data(fname, 'qisum', qisum)
    call read_data(fname, 'cfsum', cfsum)
 else
    If (file_exist('INPUT/strat_cloud.res')) Then
       unit = open_restart_file (FILE='INPUT/strat_cloud.res', &
            ACTION='read')
       if(mpp_pe() == mpp_root_pe() ) call mpp_error ('strat_cloud_mod', &
            'Reading native formatted restart file.', NOTE)
       read (unit, iostat=io, err=142) vers, vers2

142    continue
       if (io == 0) then

          !--------------------------------------------------------------------
          !    if eor is not encountered, then the file includes radturbten.
          !    that data is not needed, simply continue by reading next record.
          !--------------------------------------------------------------------
          call error_mesg ('strat_cloud_mod',  &
               'reading pre-version number strat_cloud.res file, &
               &ignoring radturbten', NOTE)

          !--------------------------------------------------------------------
          !    the file is a newer one with a version number included. read the 
          !    version number. if it is not a valid version, stop execution with
          !    a message.
          !--------------------------------------------------------------------
       else
          if (.not. any(vers == restart_versions) ) then
             write (chvers, '(i4)') vers
             call error_mesg ('strat_cloud_mod',  &
                  'restart version ' // chvers//' cannot be read &
                  &by this version of strat_cloud_mod.', FATAL)
          endif
       endif
       call read_data (unit, nsum)
       call read_data (unit, qlsum)
       call read_data (unit, qisum)
       call read_data (unit, cfsum)
       call close_file (unit)
    else
       qlsum=0.0; qisum=0.0; cfsum=0.0; nsum=0
    endif
 endif
 !-----------------------------------------------------------------------
 !
 !       Setup Diagnostics

 call diag_field_init(axes,Time)

 !-----------------------------------------------------------------------

end subroutine strat_cloud_init


!#######################################################################
!#######################################################################


! <SUBROUTINE NAME="diag_field_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!         Initializes netcdf diagnostics.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_field_init (axes,Time)
!
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="integer">
!         Integer array containing axes integers.
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!         Time 
!  </IN>
! </SUBROUTINE>
!
subroutine diag_field_init (axes,Time)


 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time

 integer, dimension(3) :: half = (/1,2,4/)


 ! assorted items

 id_droplets = register_diag_field ( mod_name, 'droplets', &
      axes(1:3), Time, 'Droplet number conentration',     &
      '/m3', missing_value=missing_value )

 id_aall = register_diag_field ( mod_name, 'aall', axes(1:3),   &
      Time, 'Cloud fraction for all clouds at midtimestep',     &
      'dimensionless', missing_value=missing_value )

 id_aliq = register_diag_field ( mod_name, 'aliq', axes(1:3),   &
      Time, 'Cloud fraction for liquid clouds', 'dimensionless',&
      missing_value=missing_value )

 id_aice = register_diag_field ( mod_name, 'aice', axes(1:3),   &
      Time, 'Cloud fraction for ice clouds', 'dimensionless',   &
      missing_value=missing_value )

 id_rvolume = register_diag_field ( mod_name, 'rv', axes(1:3),  &
      Time, 'Cloud liquid mean volume radius', 'microns',       &
      missing_value=missing_value )

 id_autocv = register_diag_field ( mod_name, 'aauto', axes(1:3),&
      Time, 'Cloud fraction where autoconversion is occurring', &
      'dimensionless', missing_value=missing_value )

 id_vfall = register_diag_field ( mod_name, 'vfall', axes(1:3), &
      Time, 'Ice crystal fall speed', 'meters/second',          &
      missing_value=missing_value )


 !liquid water tendencies


 id_qldt_cond = register_diag_field ( mod_name, &
      'qldt_cond', axes(1:3), Time, &
      'Liquid water specific humidity tendency from LS condensation', &
      'kg/kg/sec', missing_value=missing_value               )

 id_qldt_evap = register_diag_field ( mod_name, &
      'qldt_evap', axes(1:3), Time, &
      'Liquid water specific humidity tendency from LS evaporation', &
      'kg/kg/sec', missing_value=missing_value               )

 id_qldt_eros = register_diag_field ( mod_name, &
      'qldt_eros', axes(1:3), Time, &
      'Liquid water specific humidity tendency from erosion',   &
      'kg/kg/sec', missing_value=missing_value               )

 id_qldt_berg = register_diag_field ( mod_name, &
      'qldt_berg', axes(1:3), Time, &
      'Liquid water specific humidity tendency from Bergeron process',&
      'kg/kg/sec', missing_value=missing_value               )

 id_qldt_freez = register_diag_field ( mod_name, &
      'qldt_freez', axes(1:3), Time, &
      'Liquid water specific humidity tendency from homogenous freezing',&
      'kg/kg/sec', missing_value=missing_value               )

 id_qldt_rime = register_diag_field ( mod_name, &
      'qldt_rime', axes(1:3), Time, &
      'Liquid water specific humidity tendency from riming',    &
      'kg/kg/sec', missing_value=missing_value               )

 id_qldt_accr = register_diag_field ( mod_name, &
      'qldt_accr', axes(1:3), Time, &
      'Liquid water specific humidity tendency from accretion', &
      'kg/kg/sec', missing_value=missing_value               )

 id_qldt_auto = register_diag_field ( mod_name, &
      'qldt_auto', axes(1:3), Time, &
      'Liquid water specific humidity tendency from autoconversion',&
      'kg/kg/sec', missing_value=missing_value               )

 id_qldt_fill = register_diag_field ( mod_name, &
      'qldt_fill', axes(1:3), Time, &
      'Liquid water specific humidity tendency from filler',    &
      'kg/kg/sec', missing_value=missing_value               )

 id_qldt_destr = register_diag_field ( mod_name, &
      'qldt_destr', axes(1:3), Time, &
      'Liquid water specific humidity tendency from cloud destruction',&
      'kg/kg/sec', missing_value=missing_value               )


!cloud droplet number tendencies

 id_qndt_cond = register_diag_field ( mod_name, &
       'qndt_cond', axes(1:3), Time, &
       'Cloud droplet tendency from LS condensation', &
       '#/kg/sec', missing_value=missing_value               )

 id_qndt_evap = register_diag_field ( mod_name, &
       'qndt_evap', axes(1:3), Time, &
       'Cloud droplet tendency from LS evaporation', &
       '#/kg/sec', missing_value=missing_value               )

 !rain stuff


 id_liq_adj = register_diag_field ( mod_name, &
      'liq_adj', axes(1:3), Time, &
      'Liquid condensation rate from removal of supersaturation',&
      'kg/kg/sec', missing_value=missing_value               )

 id_rain_clr = register_diag_field ( mod_name, &
      'rain_clr', axes(half), Time, &
      'Clear sky rain rate averaged to grid box mean',          &
      'kg/m2/s', missing_value=missing_value                 )

 id_rain_cld = register_diag_field ( mod_name, &
      'rain_cld', axes(half), Time, &
      'cloudy sky rain rate averaged to grid box mean',         &
      'kg/m2/s', missing_value=missing_value                 )

 id_a_rain_clr = register_diag_field ( mod_name, &
      'a_rain_clr', axes(half), Time, &
      'Clear sky rain fractional coverage',                     &
      'fraction', missing_value=missing_value                 )

 id_a_rain_cld = register_diag_field ( mod_name, &
      'a_rain_cld', axes(half), Time, &
      'cloudy sky rain fractional coverage',                    &
      'fraction', missing_value=missing_value                 )

 id_rain_evap = register_diag_field ( mod_name, &
      'rain_evap', axes(1:3), Time, &
      'Water vapor tendency from rain evaporation',             &
      'kg/kg/sec', missing_value=missing_value               )

 id_a_precip_clr = register_diag_field ( mod_name, &
      'a_precip_clr', axes(half), Time, &
      'Clear sky precip fractional coverage',                   &
      'fraction', missing_value=missing_value                 )

 id_a_precip_cld = register_diag_field ( mod_name, &
      'a_precip_cld', axes(half), Time, &
      'cloudy sky precip fractional coverage',                  &
      'fraction', missing_value=missing_value                 )


 !ice water tendencies


 id_qidt_dep = register_diag_field ( mod_name, &
      'qidt_dep', axes(1:3), Time, &
      'Ice water specific humidity tendency from LS deposition', &
      'kg/kg/sec', missing_value=missing_value               )

 id_qidt_subl = register_diag_field ( mod_name, &
      'qidt_subl', axes(1:3), Time, &
      'Ice water specific humidity tendency from LS sublimation', &
      'kg/kg/sec', missing_value=missing_value               )

 id_qidt_fall = register_diag_field ( mod_name, &
      'qidt_fall', axes(1:3), Time, &
      'Ice water specific humidity tendency from ice settling', &
      'kg/kg/sec', missing_value=missing_value               )

 id_qidt_eros = register_diag_field ( mod_name, &
      'qidt_eros', axes(1:3), Time, &
      'Ice water specific humidity tendency from erosion',      &
      'kg/kg/sec', missing_value=missing_value               )

 id_qidt_melt = register_diag_field ( mod_name, &
      'qidt_melt', axes(1:3), Time, &
      'Ice water specific humidity tendency from melting to rain',&
      'kg/kg/sec', missing_value=missing_value               )

 id_qidt_fill = register_diag_field ( mod_name, &
      'qidt_fill', axes(1:3), Time, &
      'Ice water specific humidity tendency from filler',       &
      'kg/kg/sec', missing_value=missing_value               )

 id_qidt_destr = register_diag_field ( mod_name, &
      'qidt_destr', axes(1:3), Time, &
      'Ice water specific humidity tendency from cloud destruction',&
      'kg/kg/sec', missing_value=missing_value               )



 !snow stuff


 id_ice_adj = register_diag_field ( mod_name, &
      'ice_adj', axes(1:3), Time, &
      'Frozen condensation rate from removal of supersaturation', &
      'kg/kg/sec', missing_value=missing_value               )

 id_snow_clr = register_diag_field ( mod_name, &
      'snow_clr', axes(half), Time, &
      'Clear sky snow rate averaged to grid box mean',          &
      'kg/m2/s', missing_value=missing_value                 )

 id_snow_cld = register_diag_field ( mod_name, &
      'snow_cld', axes(half), Time, &
      'cloudy sky snow rate averaged to grid box mean',         &
      'kg/m2/s', missing_value=missing_value                 )

 id_a_snow_clr = register_diag_field ( mod_name, &
      'a_snow_clr', axes(half), Time, &
      'Clear sky snow fractional coverage',                     &
      'fraction', missing_value=missing_value                 )

 id_a_snow_cld = register_diag_field ( mod_name, &
      'a_snow_cld', axes(half), Time, &
      'cloudy sky snow fractional coverage',                    &
      'fraction', missing_value=missing_value                 )

 id_snow_subl = register_diag_field ( mod_name, &
      'snow_subl', axes(1:3), Time, &
      'Water vapor tendency from snow sublimation',             &
      'kg/kg/sec', missing_value=missing_value               )

 id_snow_melt = register_diag_field ( mod_name, &
      'snow_melt', axes(1:3), Time, &
      'Rain water tendency from snow melting',                  &
      'kg/kg/sec', missing_value=missing_value               )


 !cloud fraction tendencies

 id_qadt_lsform = register_diag_field ( mod_name, &
      'qadt_lsform', axes(1:3), Time, &
      'cloud fraction tendency from LS condensation',                 &
      '1/sec', missing_value=missing_value               )

 id_qadt_lsdiss = register_diag_field ( mod_name, &
      'qadt_lsdiss', axes(1:3), Time, &
      'cloud fraction tendency from LS evaporation',                 &
      '1/sec', missing_value=missing_value               )

 id_qadt_rhred = register_diag_field ( mod_name, &
      'qadt_rhred', axes(1:3), Time, &
      'cloud fraction tendency from RH limiter',                      &
      '1/sec', missing_value=missing_value               )

 id_qadt_eros = register_diag_field ( mod_name, &
      'qadt_eros', axes(1:3), Time, &
      'cloud fraction tendency from erosion',                   &
      '1/sec', missing_value=missing_value               )

 id_qadt_fill = register_diag_field ( mod_name, &
      'qadt_fill', axes(1:3), Time, &
      'cloud fraction tendency from filler',                    &
      '1/sec', missing_value=missing_value               )

 id_qadt_super = register_diag_field ( mod_name, &
      'qadt_super', axes(1:3), Time, &
      'cloud fraction tendency from supersaturation formation', &
      '1/sec', missing_value=missing_value               )

 id_qadt_destr = register_diag_field ( mod_name, &
      'qadt_destr', axes(1:3), Time, &
      'cloud fraction tendency from cloud destruction',             &
      '1/sec', missing_value=missing_value               )


 !column integrated liquid tendencies


 id_ql_cond_col = register_diag_field ( mod_name, &
      'ql_cond_col', axes(1:2), Time, &
      'Column integrated condensation',                         &
      'kg/m2/sec', missing_value=missing_value               )

 id_ql_evap_col = register_diag_field ( mod_name, &
      'ql_evap_col', axes(1:2), Time, &
      'Column integrated evaporation',                          &
      'kg/m2/sec', missing_value=missing_value               )

 id_ql_eros_col = register_diag_field ( mod_name, &
      'ql_eros_col', axes(1:2), Time, &
      'Column integrated liquid erosion',                       &
      'kg/m2/sec', missing_value=missing_value               )

 id_ql_accr_col = register_diag_field ( mod_name, &
      'ql_accr_col', axes(1:2), Time, &
      'Column integrated accretion',                            &
      'kg/m2/sec', missing_value=missing_value               )

 id_ql_berg_col = register_diag_field ( mod_name, &
      'ql_berg_col', axes(1:2), Time, &
      'Column integrated Bergeron process',                     &
      'kg/m2/sec', missing_value=missing_value               )

 id_ql_freez_col = register_diag_field ( mod_name, &
      'ql_freez_col', axes(1:2), Time, &
      'Column integrated homogeneous freezing',                 &
      'kg/m2/sec', missing_value=missing_value               )

 id_ql_destr_col = register_diag_field ( mod_name, &
      'ql_destr_col', axes(1:2), Time, &
      'Column integrated liquid destruction',                   &
      'kg/m2/sec', missing_value=missing_value               )

 id_ql_rime_col = register_diag_field ( mod_name, &
      'ql_rime_col', axes(1:2), Time, &
      'Column integrated riming',                               &
      'kg/m2/sec', missing_value=missing_value               )

 id_ql_auto_col = register_diag_field ( mod_name, &
      'ql_auto_col', axes(1:2), Time, &
      'Column integrated autoconversion',                       &
      'kg/m2/sec', missing_value=missing_value               )

 id_ql_fill_col = register_diag_field ( mod_name, &
      'ql_fill_col', axes(1:2), Time, &
      'Column integrated liquid filler',                        &
      'kg/m2/sec', missing_value=missing_value               )

 id_liq_adj_col = register_diag_field ( mod_name, &
      'liq_adj_col', axes(1:2), Time, &
      'Column integrated liquid condensation by adjustment',    &
      'kg/m2/sec', missing_value=missing_value               )

 id_rain_evap_col = register_diag_field ( mod_name, &
      'rain_evap_col', axes(1:2), Time, &
      'Column integrated rain evaporation',                     &
      'kg/m2/sec', missing_value=missing_value               )



 !column integrated ice tendencies


 id_qi_fall_col = register_diag_field ( mod_name, &
      'qi_fall_col', axes(1:2), Time, &
      'Column integrated ice settling',                         &
      'kg/m2/sec', missing_value=missing_value               )

 id_qi_fill_col = register_diag_field ( mod_name, &
      'qi_fill_col', axes(1:2), Time, &
      'Column integrated ice filler',                           &
      'kg/m2/sec', missing_value=missing_value               )

 id_qi_eros_col = register_diag_field ( mod_name, &
      'qi_eros_col', axes(1:2), Time, &
      'Column integrated ice erosion',                          &
      'kg/m2/sec', missing_value=missing_value               )

 id_qi_dep_col = register_diag_field ( mod_name, &
      'qi_dep_col', axes(1:2), Time, &
      'Column integrated large-scale deposition',               &
      'kg/m2/sec', missing_value=missing_value               )

 id_qi_subl_col = register_diag_field ( mod_name, &
      'qi_subl_col', axes(1:2), Time, &
      'Column integrated large-scale sublimation',              &
      'kg/m2/sec', missing_value=missing_value               )

 id_qi_destr_col = register_diag_field ( mod_name, &
      'qi_destr_col', axes(1:2), Time, &
      'Column integrated ice destruction',                      &
      'kg/m2/sec', missing_value=missing_value               )

 id_qi_melt_col = register_diag_field ( mod_name, &
      'qi_melt_col', axes(1:2), Time, &
      'Column integrated ice melting',                          &
      'kg/m2/sec', missing_value=missing_value               )

 id_ice_adj_col = register_diag_field ( mod_name, &
      'ice_adj_col', axes(1:2), Time, &
      'Column integrated frozen condesation by adjustment',     &
      'kg/m2/sec', missing_value=missing_value               )

 id_snow_subl_col = register_diag_field ( mod_name, &
      'snow_subl_col', axes(1:2), Time, &
      'Column integrated snow sublimation',                     &
      'kg/m2/sec', missing_value=missing_value               )

 id_snow_melt_col = register_diag_field ( mod_name, &
      'snow_melt_col', axes(1:2), Time, &
      'Column integrated snow melting',                         &
      'kg/m2/sec', missing_value=missing_value               )


 !column integrated cloud fraction tendencies


 id_qa_lsform_col = register_diag_field ( mod_name, &
      'qa_lsform_col', axes(1:2), Time, &
      'Column integrated large-scale formation',                &
      'kg/m2/sec', missing_value=missing_value               )

 id_qa_lsdiss_col = register_diag_field ( mod_name, &
      'qa_lsdiss_col', axes(1:2), Time, &
      'Column integrated large-scale dissipation',                &
      'kg/m2/sec', missing_value=missing_value               )

 id_qa_rhred_col = register_diag_field ( mod_name, &
      'qa_rhred_col', axes(1:2), Time, &
      'Column integrated RH reduction',                         &
      'kg/m2/sec', missing_value=missing_value               )

 id_qa_eros_col = register_diag_field ( mod_name, &
      'qa_eros_col', axes(1:2), Time, &
      'Column integrated cloud fraction erosion',               &
      'kg/m2/sec', missing_value=missing_value               )

 id_qa_fill_col = register_diag_field ( mod_name, &
      'qa_fill_col', axes(1:2), Time, &
      'Column integrated cloud fraction filler',                &
      'kg/m2/sec', missing_value=missing_value               )

 id_qa_super_col = register_diag_field ( mod_name, &
      'qa_super_col', axes(1:2), Time, &
      'Column integrated cloud fraction supersaturation'//      &
      ' formation', 'kg/m2/sec', missing_value=missing_value    )

 id_qa_destr_col = register_diag_field ( mod_name, &
      'qa_destr_col', axes(1:2), Time, &
      'Column integrated cloud fraction destruction',           &
      'kg/m2/sec', missing_value=missing_value               )


 !-----------------------------------------------------------------------
 !
 !       set diagnostic flag

 do_budget_diag = .false.
 if ( id_qldt_cond     > 0 .or. id_qldt_eros     > 0 .or. &
      id_qldt_fill     > 0 .or. id_qldt_accr     > 0 .or. &
      id_qldt_evap     > 0 .or. id_qldt_freez    > 0 .or. &
      id_qldt_berg     > 0 .or. id_qldt_destr    > 0 .or. &
      id_qldt_rime     > 0 .or. id_qldt_auto     > 0 .or. &
      id_qndt_cond     > 0 .or. id_qndt_evap     > 0) then
    do_budget_diag = .true.
 end if
 if ( id_rain_clr      > 0 .or. id_rain_cld      > 0 .or. &
      id_a_rain_clr    > 0 .or. id_a_rain_cld    > 0 .or. &
      id_rain_evap     > 0 .or. id_liq_adj       > 0) then
    do_budget_diag = .true.
 end if
 if ( id_qidt_fall     > 0 .or. id_qidt_fill     > 0 .or. &
      id_qidt_melt     > 0 .or. id_qidt_dep      > 0 .or. &
      id_qidt_subl     > 0 .or. id_qidt_eros     > 0 .or. &
      id_qidt_destr    > 0) then
    do_budget_diag = .true.
 end if
 if ( id_snow_clr      > 0 .or. id_snow_cld      > 0 .or. &
      id_a_snow_clr    > 0 .or. id_a_snow_cld    > 0 .or. &
      id_snow_subl     > 0 .or. id_snow_melt     > 0 .or. &
      id_ice_adj       > 0 ) then
    do_budget_diag = .true.
 end if
 if ( id_ql_eros_col   > 0 .or. id_ql_cond_col   > 0 .or. &
      id_ql_evap_col   > 0 .or. id_ql_accr_col   > 0 .or. &
      id_ql_auto_col   > 0 .or. id_ql_fill_col   > 0 .or. &
      id_ql_berg_col   > 0 .or. id_ql_destr_col  > 0 .or. &
      id_ql_rime_col   > 0 .or. id_ql_freez_col  > 0) then
    do_budget_diag = .true.
 end if
 if ( id_rain_evap_col > 0 .or. id_liq_adj_col   > 0 ) then
    do_budget_diag = .true. 
 end if
 if ( id_qi_fall_col   > 0 .or. id_qi_fill_col   > 0 .or. &
      id_qi_subl_col   > 0 .or. id_qi_melt_col   > 0 .or. &
      id_qi_destr_col  > 0 .or. id_qi_eros_col   > 0 .or. &
      id_qi_dep_col    > 0) then
    do_budget_diag = .true.
 end if
 if ( id_snow_subl_col > 0 .or. id_snow_melt_col > 0 .or. &
      id_ice_adj_col   > 0 ) then
    do_budget_diag = .true.
 end if
 if ( id_qadt_lsform   > 0 .or. id_qadt_eros     > 0 .or. &
      id_qadt_fill     > 0 .or. id_qadt_rhred    > 0 .or. &
      id_qadt_destr    > 0 .or. id_qa_lsdiss_col > 0 .or. &
      id_qadt_lsdiss   > 0 .or. id_qadt_super    > 0 .or. &
      id_qa_lsform_col > 0 .or. id_qa_super_col  > 0 .or. &
      id_qa_eros_col   > 0 .or. id_qa_fill_col   > 0 .or. &
      id_qa_rhred_col  > 0 .or. id_qa_destr_col  > 0) then
    do_budget_diag = .true.
 end if
 if ( id_a_precip_cld  > 0 .or. id_a_precip_clr  > 0 ) then
    do_budget_diag = .true.
 end if

 !-----------------------------------------------------------------------



end subroutine diag_field_init


!#######################################################################
!#######################################################################


! <SUBROUTINE NAME="strat_cloud">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!       
!  </DESCRIPTION>
!  <TEMPLATE>
!   call strat_cloud(Time,is,ie,js,je,dtcloud,pfull,phalf,radturbten2, 
!                   T,qv,ql,qi,qa,omega,Mc,diff_t,LAND,              
!                   ST,SQ,SL,SI,SA,surfrain,                         
!                   surfsnow,qrat,ahuco,MASK, 
!                    qn, SN, Aerosol)
!
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!         Time
!  </IN>
!  <IN NAME="is" TYPE="integer">
!         Indice of starting point in the longitude direction of the slab being passed to strat_cloud
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!         Indice of ending point in the longitude direction of the slab being passed 
!  </IN>
!  <IN NAME="js" TYPE="integer">
!         Indice of starting point in the latitude direction of the slab being passed
!  </IN>
!  <IN NAME="je" TYPE="integer">
!         Indice of ending point in the latitude direction of the slab being passed 
!  </IN>
!  <IN NAME="dtcloud" TYPE="real">
!         Physics time step (sec)
!  </IN>
!  <IN NAME="pfull" TYPE="real">
!         Pressure on model full levels (Pa)
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!         Pressure on model half levels (Pa)
!  </IN>
!  <IN NAME="radturbten2" TYPE="real">
!         Sum of the tendencies of temperature from turbulence and radiation schemes (K/s)
!  </IN>
!  <IN NAME="T" TYPE="real">
!         Temperature (K)         
!  </IN>
!  <IN NAME="qv" TYPE="real">
!         Water vapor specific humidity (kg vapor/kg air)
!  </IN>
!  <IN NAME="ql" TYPE="real">
!         Grid-box mean liquid water specific humidity (kg liquid/kg air)
!  </IN>
!  <IN NAME="qi" TYPE="real">
!         Grid-box mean ice water specific humidity (kg ice/kg air)
!  </IN>
!  <IN NAME="qa" TYPE="real">
!         Cloud fraction (3d array and a prognostic variable) (fraction)
!  </IN>
!  <IN NAME="qn" TYPE="real">
!         Cloud droplet number (3d array and a prognostic variable) (#/kg air)
!  </IN>
!  <IN NAME="omega" TYPE="real">
!         Vertical pressure velocity (Pa/sec)
!  </IN>
!  <IN NAME="Mc" TYPE="real">
!         Cumulus mass flux (defined positive as upward) (kg air/m2/sec)
!  </IN>
!  <IN NAME="diff_t" TYPE="real">
!         Vertical diffusion coefficient for temperature and tracer from vertical diffusion scheme (m2/sec) 
!  </IN>
!  <IN NAME="LAND" TYPE="real">
!         Fraction of surface that contains land (fraction)
!  </IN>
!  <OUT NAME="ST" TYPE="real">
!         Change in temperature due to strat_cloud (K) 
!  </OUT>
!  <OUT NAME="SQ" TYPE="real">
!         Change in water vapor due to strat_cloud (kg vapor/kg air) 
!  </OUT>
!  <OUT NAME="SL" TYPE="real">
!         Change in cloud liquid due to strat_cloud (kg liquid/kg air)
!  </OUT>
!  <OUT NAME="SI" TYPE="real">
!         Change in cloud ice due to strat_cloud (kg ice/kg air)
!  </OUT>
!  <OUT NAME="SA" TYPE="real">
!         Change in cloud fraction due to strat_cloud (fraction)
!  </OUT>
!  <OUT NAME="SN" TYPE="real">
!         Change in cloud droplet number due to strat_cloud (fraction)
!  </OUT>
!  <OUT NAME="surfrain" TYPE="real">
!         Surface rain fall over time step dtcloud (kg liquid/m2)
!  </OUT>
!  <OUT NAME="surfsnow" TYPE="real">
!         Surface snow fall over time step dtcloud (kg ice/m2)
!  </OUT>
!  <OUT NAME="rain3d" TYPE="real">
!         3D rain fall over time step dtcloud (kg liquid/m2)
!  </OUT>
!  <OUT NAME="snow3d" TYPE="real">
!         3D snow fall over time step dtcloud (kg ice/m2)
!  </OUT>
!  <IN NAME="qrat" TYPE="real">
!         Ratio of large-scale specific humidity to specific humidity in 
!         environment outside convective system (from donner_deep) 
!         
!         Will be equal to 1 for all normal AM2 operations (i.e. donner_deep is not activated)              
!         
!         Note that index 1 is nearest ground
!
!  </IN>
!  <IN NAME="ahuco" TYPE="real">
!         The fraction of the grid box containing either cumulus cells or the mesoscale circulation (from donner_deep).
!
!         Will be equal to 0 for all normal AM2 operations (i.e. donner_deep is not activated)              
!         
!         Note that index 1 is nearest ground
!
!  </IN>
!  <IN NAME="MASK" TYPE="real">
!         Optional input real array indicating the point is above the surface
!         if equal to 1.0 and indicating the point is below the surface if 
!         equal to 0.
!
!         Used only in eta vertical coordinate model.
!  </IN>
! </SUBROUTINE>
!


subroutine strat_cloud(Time,is,ie,js,je,dtcloud,pfull,phalf,radturbten2,&
    T,qv,ql,qi,qa,omega,Mc,diff_t,LAND,              &
    ST,SQ,SL,SI,SA,rain3d,snow3d,surfrain,                         &
    surfsnow,qrat,ahuco,MASK, &
    qn, Aerosol, SN)


  !        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !
  !       VARIABLES
  !
  !
  !
  !       ------
  !       INPUT:
  !       ------
  !
  !         variable              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !       Time           time type variable 
  !
  !       is,ie          starting and ending i indices 
  !                      for data window
  !
  !       js,je          starting and ending j indices 
  !                      for data window
  !
  !       dtcloud        time between this call and      s
  !                      the next call to strat_cloud
  !
  !       pfull          pressure at full model levels   Pa
  !                      IMPORTANT NOTE: p(j)<p(j+1)
  !
  !       phalf          pressure at half model levels   Pa
  !                      phalf(j)<pfull(j)<phalf(j+1)
  !
  !       T              temperature                     K
  !
  !       qv             specific humidity of water      kg vapor/kg air
  !                      vapor
  !
  !       ql             specific humidity of cloud      kg condensate/
  !                      liquid                          kg air
  !
  !       qi             specific humidity of cloud      kg condensate/
  !                      ice                             kg air
  !
  !       qa             saturated volume fraction       fraction
  !
  !       qn             cloud droplet number            #/kg air
  !
  !       qrat           ratio of large-scale spec       fraction
  !                      humidity to spec humidity
  !                      in environment outside
  !                      convective system (from
  !                      donner_deep) 
  !                      index 1 nearest ground
  !
  !       ahuco          fraction, cell+meso, from       fraction
  !                      donner_deep
  !                      index 1 nearest ground
  !
  !       omega          vertical pressure velocity      Pa/s
  !
  !       Mc             Cumulus mass flux defined       kg/(m*m)/s
  !                      on full levels
  !
  !       diff_t         Vertical diffusion coefficient  (m*m)/s
  !                      for temperature and tracer
  !
  !       LAND           the fraction of the grid box    fraction
  !                      covered by land
  !                               
  !
  !       -------
  !       OUTPUT:
  !       -------
  !
  !         variable              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !       ST             temperature change due to       K
  !                      all stratiform processes
  !
  !       SQ             water vapor change due to       kg vapor/kg air
  !                      all stratiform processes
  !
  !       SL             cloud liquid change due to      kg condensate/
  !                      all stratiform processes        kg air
  !
  !       SI             cloud ice change due to         kg condensate/
  !                      all stratiform processes        kg air
  !
  !       SA             saturated volume fraction       fraction
  !                      change due to all stratiform 
  !                      processes
  !
  !       SN             cloud droplet number            #/kg air
  !                      change due to all stratiform 
  !                      processes
  !
  !       surfrain       rain that falls through the     kg condensate/
  !                      bottom of the column over       (m*m)
  !                      the time dtcloud
  !
  !       surfsnow       snow that falls through the     kg condensate/
  !                      bottom of the column over       (m*m)
  !                      the time dtcloud
  !
  !       rain3d         rain that falls through the     kg condensate/
  !                      each of the model layer         (m*m)/sec
  !
  !       snow3d         snow that falls through the     kg condensate/
  !                      each of the model layer         (m*m)/sec
  !
  !
  !       ---------------
  !       optional INPUT:
  !       ---------------
  !
  !         variable              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !
  !       MASK           real array indicating the 
  !                      point is above the surface
  !                      if equal to 1.0 and 
  !                      indicating the point is below
  !                      the surface if equal to 0.
  !
  !       -------------------
  !       INTERNAL VARIABLES:
  !       -------------------
  !
  !         variable              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !       kdim           number of vertical levels
  !
  !       j              model vertical level being 
  !                      processed
  !
  !       ipt,jpt        i and j point indice used only
  !                      in instantaneous diag-
  !                      nostic output
  !
  !       i,unit,nn      temporary integers used in
  !                      instantaneous diagnostic
  !                      output.
  !
  !       inv_dtcloud    1 / dtcloud                     1/sec
  !
  !       airdens        air density                     kg air/(m*m*m)
  !
  !       qs             saturation specific humidity    kg vapor/kg air
  !
  !       dqsdT          T derivative of qs              kg vapor/kg air/K
  !
  !       gamma          (L/cp)*dqsdT                    dimensionless
  !
  !       rain_clr       grid mean flux of rain enter-   kg condensate/
  !                      ing the grid box from above     (m*m)/s
  !                      and entering the unsaturated 
!                      portion of the grid box
!
!       rain_cld       grid mean flux of rain enter-   kg condensate/
!                      ing the grid box from above     (m*m)/s
!                      and entering the saturated 
!                      portion of the grid box
!
!       a_rain_clr     fraction of grid box occupied   fraction
!                      by rain_clr
!
!       a_rain_cld     fraction of grid box occupied   fraction
!                      by rain_cld
!
!       snow_cld       flux of ice entering the        kg condensate/
!                      saturated portion of the        (m*m)/s
!                      grid box from above by means
!                      of gravitational settling 
!
!       snow_clr       flux of ice outside of cloud    kg condensate/
!                      entering the unsaturated        (m*m)/s
!                      portion of the grid box from      
!                      above
!
!       a_snow_clr     area fraction of grid box       fraction
!                      covered by snow flux in
!                      unsaturated air
!
!       a_snow_cld     area fraction of grid box       fraction
!                      covered by the snow flux in
!                      saturated air
!
!       deltpg         pressure thickness of grid box  kg air/(m*m)
!                      divided by gravity
!
!       U              grid box relative humidity      fraction
!
!       U00p           critical relative humidity      fraction
!                      which is a function of pressure 
!
!       dqs_ls         change in saturation specific   kg vapor/kg air
!                      due to large-scale processes,
!                      such as large-scale vertical
!                      motion, compensating convective
!                      mass flux, or radiative cooling
!
!       da_ls          change in saturated volume      fraction
!                      fraction due to large-scale
!                      processes
!
!       C_dt           product of A and dtcloud in     dimensionless in 
!                      in the analytic integration     qa integration
!                      of the qa equation, or C and
!                      dtcloud in the analytic         kg condensate/
!                      integration of the ql and qi    kg air in ql or 
!                      equations.                      qi integration
!
!       D_dt           product of B and dtcloud in     dimensionless in
!                      in the analytic integration     qa, ql, and qi
!                      of the qa equation, or D and    integration
!                      dtcloud in the analytic         
!                      integration of the ql and qi    
!                      equations.                      
!
!       qceq           equilibrium value of cloud      dimensionless or
!                      fraction or cloud condensate    kg condensate /
!                      that the analytic integration   kg air
!                      approaches                      
!
!       qcbar          mean value of cloud fraction    dimensionless or
!                      or cloud condensate over the    kg condensate /
!                      t0 to t0 + dtcloud interval     kg air
!
!       qc0            value of cloud fraction or      dimensionless or
!                      cloud condensate at the         kg condensate /
!                      initial time                    kg air
!        
!       qc1            value of cloud fraction or      dimensionless or
!                      cloud condensate at the final   kg condensate /
!                      time                            kg air
!       
!       D1_dt          first sink in either ql or qi   dimensionless
!                      equation. This is analogous to
!                      D_dt.  In ql equation, this 
!                      sink represents the conversion
!                      of cloud liquid to rain. In the
!                      qi equation it represents the
!                      settling of ice crystals.
!
!       D2_dt          second sink in ql or qi         dimensionless
!                      equation. This is analogous 
!                      to D_dt. In ql equation this
!                      sink represents the conversion
!                      of cloud liquid to ice. In the
!                      qi equation this sink 
!                      represents the melting of 
!                      cloud ice into rain.
!
!       D_eros         Sink in ql, qi and qa equation  dimensionless
!                      due to turbulent erosion of
!                      cloud sides
! 
!       ql_upd         updated value of ql             kg condensate/
!                                                      kg air
!       
!       qi_upd         updated value of qi             kg condensate/
!                                                      kg air
!
!       qa_upd         updated value of qa             fraction
!
!       qa_mean        qa + SA; semi-implicit          fraction
!                      saturated volume fraction
!
!       qa_mean_lst    qa_mean of the level above      fraction
!
!       ql_mean        ql + positive increment         kg condensate/
!                      of ql; i.e. a sort of           kg air
!                      intermediate ql
!
!       qi_mean        ql + positive increment         kg condensate/
!                      of qi; i.e. a sort of           kg air
!                      intermediate qi
!
!       dcond_ls       change in condensate due to     kg condensate/
!                      non-convective condensation.    kg air
!                      After phase determination,
!                      this variable refers only to
!                      liquid condensation.
!
!       dcond_ls_ice   change in ice due to            kg condensate/
!                      non-convective condensation.    kg air
!
!       da_cld2clr     fraction of the area in which   fraction
!                      rain/snow in saturated volume 
!                      above falls into unsaturated 
!                      volume in the current layer.
!
!       da_clr2cld     as in da_cld2clr except for     fraction
!                      the transfer from unsaturated
!                      to saturated volume
!
!       dprec_cld2clr  grid mean flux that is trans-   kg condensate/
!                      ferred from rain/snow in        (m*m)/s
!                      saturated volume to rain/snow 
!                      in unsaturated volume at layer 
!                      interfaces.
!
!       dprec_clr2cld  as in dprec_cld2clr except for  kg condensate/
!                      the transfer from unsaturated   (m*m)/s
!                      to saturated volume.
!
!       N              fixed number of cloud drops     1/(m*m*m)
!                      per unit volume in liquid
!                      clouds
!     
!       rad_liq        mean volume radius of liquid    microns
!                      cloud drops
!
!       A_plus_B       sum of vapor diffusion factor   m*s/kg
!                      and thermal conductivity factor
!                      which is used in various 
!                      microphysical formula for the 
!                      evaporation of rain and snow
!
!       Vfall          fall speed of ice crystals      m/s
!
!       lamda_f        slope factor in the SIZE        1/m
!                      distribution of ice crystals
!
!       U_clr          relative humidity in the clear  fraction
!                      portion of the grid box.
!
!       tmp1,tmp2,tmp3 temporary numbers used at 
!                      several points within the
!                      subroutine
!
!
!                STATISTICAL CLOUD SCHEME VARIABLES
!
!       qag            equilibrium value of cloud      dimensionless 
!                      fraction for statistical 
!                      cloud scheme           
!
!       qcg            equilibrium value of cloud      kg condensate /
!                      condensate that PDF clouds      kg air
!                      wants                          
!
!       qcg_ice        equilibrium value of cloud      kg condensate /
!                      ice condensate that PDF clouds  kg air
!                      wants
!
!       qvg            equilibrium value of water      kg vapor /
!                      vapor in the clear portion      kg air
!                      of the grid box that PDF 
!                      clouds wants                          
!
!       qtbar          total water specific humidity   kg water /
!                      which is equal to the sum of    kg air
!                      liquid water, ice water, and
!                      water vapor
!
!       deltaQ         the width of the total water    kg water /
!                      subgrid distribution (= qtmax   kg air
!                      minus qtmin)
!
!       qtmin          the minimum value to the total  kg water /
!                      sub-grid scale distribution     kg air
!
!       qs_norm        the difference between the      dimensionless
!                      saturation specific humidity    
!                      and qtmin normalized by deltaQ
!
!       icbp           the value of the incomplete     dimensionless
!                      beta function evaluated with
!                      x=qs_norm, p=betaP, and q=betaP
!
!       icbp1          the value of the incomplete     dimensionless
!                      beta function evaluated with                  
!                      x=qs_norm, p=betaP+1, and 
!                      q=betaP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


!        
!       user Interface variables
!       ------------------------
!

        type(time_type), intent (in)           :: Time
        integer, intent (in)                   :: is,ie,js,je
        real, intent (in)                      :: dtcloud
        real, intent (in),    dimension(:,:,:) :: pfull,phalf
        real, intent (in),    dimension(:,:,:) :: T,qv,ql,qi,qa,omega
        real, intent (in),    dimension(:,:,:) :: Mc, diff_t
        real, intent (in),    dimension(:,:,:) :: qrat,ahuco
        real, intent (in),    dimension(:,:)   :: LAND
        real, intent (in),    dimension(:,:,:) :: radturbten2
        real, intent (out),   dimension(:,:,:) :: ST,SQ,SL,SI,SA
        real, intent (out),   dimension(:,:)   :: surfrain,surfsnow
        real, intent (in), optional, dimension(:,:,:) :: MASK
        real, intent (in),  optional, dimension(:,:,:) :: qn
        real, intent (out),   dimension(:,:,:) :: rain3d
        real, intent (out),   dimension(:,:,:) :: snow3d
        type(aerosol_type), intent (in), optional      :: Aerosol  
        real, intent (out), optional, dimension(:,:,:) :: SN

!
!       Internal variables
!       ------------------
!

        integer                                        :: idim,jdim,kdim
        integer                                        :: id,jd,ns
        integer                                        :: j,ipt,jpt
        integer                                        :: i,unit,nn        
!ljd
        integer                                        :: itest,jtest
!ljd
        real                                           :: inv_dtcloud
        real                                           :: icbp, icbp1, pnorm
        real, dimension(size(T,1),size(T,2),size(T,3)) :: airdens
        real, dimension(size(T,1),size(T,2),size(T,3)) :: qs,dqsdT
        real, dimension(size(T,1),size(T,2),size(T,3)) :: gamma
        real, dimension(size(T,1),size(T,2),size(T,3)) :: A_plus_B
        real, dimension(4,size(T,1),size(T,2),size(T,3)) :: qta4,qtqsa4
        real, dimension(size(T,1),size(T,2),size(T,3)) :: delp
        real, dimension(size(T,1),size(T,2))   :: rain_clr,rain_cld
        real, dimension(size(T,1),size(T,2))   :: a_rain_clr,a_rain_cld
        real, dimension(size(T,1),size(T,2))   :: snow_clr,snow_cld
        real, dimension(size(T,1),size(T,2))   :: a_snow_clr,a_snow_cld
        real, dimension(size(T,1),size(T,2))   :: deltpg,U,U00p
!ljd
        real, dimension(size(T,1),size(T,2))   :: U00pr
!ljd
        real, dimension(size(T,1),size(T,2))   :: dqs_ls,da_ls
        real, dimension(size(T,1),size(T,2))   :: C_dt, D_dt
        real, dimension(size(T,1),size(T,2))   :: D1_dt,D2_dt,D_eros
        real, dimension(size(T,1),size(T,2))   :: qcg_ice, qvg
        real, dimension(size(T,1),size(T,2))   :: qceq, qcbar, qcg, qag
        real, dimension(size(T,1),size(T,2))   :: qagtmp,qcgtmp,qvgtmp
        real, dimension(size(T,1),size(T,2))   :: qc1, qc0
        real, dimension(size(T,1),size(T,2))   :: Dterm
        real, dimension(size(T,1),size(T,2))   :: ql_upd,qi_upd,qa_upd, qn_upd
        real, dimension(size(T,1),size(T,2))   :: ql_mean,qi_mean, qn_mean
        real, dimension(size(T,1),size(T,2))   :: qa_mean,qa_mean_lst
        real, dimension(size(T,1),size(T,2))   :: dcond_ls,dcond_ls_ice
        real, dimension(size(T,1),size(T,2))   :: da_cld2clr,da_clr2cld
        real, dimension(size(T,1),size(T,2))   :: dprec_clr2cld
        real, dimension(size(T,1),size(T,2))   :: dprec_cld2clr
        real, dimension(size(T,1),size(T,2))   :: N,rad_liq
        real, dimension(size(T,1),size(T,2),size(T,3)) :: N3D,concen, &
                                                concen_na, concen_an
        real, dimension(size(T,1),size(T,2))   :: Vfall,lamda_f
        real, dimension(size(T,1),size(T,2))   :: U_clr
        real, dimension(size(T,1),size(T,2))   :: tmp1,tmp2,tmp3,tmp5,drop1
        real, dimension(size(T,1),size(T,2))   :: qtbar,deltaQ
        real, dimension(size(T,1),size(T,2))   :: qtmin,qs_norm          

        !diagnostic variables
        real, allocatable, dimension(:,:,:) :: &
             qldt_cond,  qldt_evap, qldt_berg, qldt_freez,&
             qldt_rime,  qldt_accr, qldt_auto, qldt_fill, &
             qldt_destr, qldt_eros, liq_adj,     rain_evap,            &
             qidt_dep,   qidt_subl, qidt_fill, &
             qidt_melt,  qidt_fall, qidt_destr,  qidt_eros, ice_adj,   &
             snow_subl,  snow_melt, qadt_lsform, qadt_eros, qadt_rhred,&
             qadt_destr, qadt_fill, qadt_lsdiss, qadt_super           ,&
             qndt_cond, qndt_evap
        real, allocatable, dimension(:,:,:) :: &
             rain_clr_diag,     rain_cld_diag,    a_rain_clr_diag, &
             a_rain_cld_diag,   snow_clr_diag,    snow_cld_diag,   &
             a_snow_clr_diag,   a_snow_cld_diag,  mask3,           &
             a_precip_clr_diag, a_precip_cld_diag
        real, allocatable, dimension(:,:,:) :: areaall, arealiq,   &
             areaice, areaautocv, rvolume, vfalldiag
        logical :: used, cloud_generator_on

        integer :: year, month, day, hour, minute, second

        integer            :: k,s,ix,jx,kx,sx
        real               :: thickness, up_strat, drop
        real, dimension(3) :: totalmass, mass_ratio
        real, dimension(size(T,1),size(T,2),size(T,3),3) :: totalmass1

!-----------------------------------------------------------------------
!       
!
     if (.not.module_is_initialized) then
       call error_mesg('strat_cloud','strat_cloud is not initialized',FATAL)
     endif

!-----------------------------------------------------------------------
!       
!       allocate space for diagnostic variables if needed
!

     if (id_aall > 0) then
             if (allocated(areaall)) deallocate (areaall)
             allocate(areaall(size(T,1),size(T,2),size(T,3)))
             areaall(:,:,:) = 0.0
     end if
     if (id_aliq > 0 .or. id_rvolume > 0) then
             if (allocated(arealiq)) deallocate (arealiq)
             allocate(arealiq(size(T,1),size(T,2),size(T,3)))
             arealiq(:,:,:) = 0.0
     end if
     if (id_aice > 0 .or. id_vfall > 0) then
             if (allocated(areaice)) deallocate (areaice)
             allocate(areaice(size(T,1),size(T,2),size(T,3)))
             areaice(:,:,:) = 0.0
     end if
     if (id_rvolume > 0) then
             if (allocated(rvolume)) deallocate (rvolume)
             allocate(rvolume(size(T,1),size(T,2),size(T,3)))
             rvolume(:,:,:) = 0.0
     end if
     if (id_autocv > 0) then
             if (allocated(areaautocv)) deallocate (areaautocv)
             allocate(areaautocv(size(T,1),size(T,2),size(T,3)))
             areaautocv(:,:,:) = 0.0
     end if
     if (id_vfall > 0) then
             if (allocated(vfalldiag)) deallocate (vfalldiag)
             allocate(vfalldiag(size(T,1),size(T,2),size(T,3)))
             vfalldiag(:,:,:) = 0.0
     end if

     if (do_budget_diag) then

     if (do_liq_num) then
       if (allocated(qndt_cond)) deallocate (qndt_cond)
               allocate(qndt_cond(size(T,1),size(T,2),size(T,3)))
       if (allocated(qndt_evap)) deallocate (qndt_evap)
               allocate(qndt_evap(size(T,1),size(T,2),size(T,3)))
     endif
     if (allocated(qldt_cond)) deallocate (qldt_cond)
             allocate(qldt_cond(size(T,1),size(T,2),size(T,3)))
     if (allocated(qldt_evap)) deallocate (qldt_evap)
             allocate(qldt_evap(size(T,1),size(T,2),size(T,3)))
     if (allocated(qldt_eros)) deallocate (qldt_eros)
             allocate(qldt_eros(size(T,1),size(T,2),size(T,3)))
     if (allocated(qldt_berg)) deallocate (qldt_berg)
             allocate(qldt_berg(size(T,1),size(T,2),size(T,3)))
     if (allocated(qldt_freez)) deallocate (qldt_freez)
             allocate(qldt_freez(size(T,1),size(T,2),size(T,3)))
     if (allocated(qldt_rime)) deallocate (qldt_rime)
             allocate(qldt_rime(size(T,1),size(T,2),size(T,3)))
     if (allocated(qldt_accr)) deallocate (qldt_accr)
             allocate(qldt_accr(size(T,1),size(T,2),size(T,3)))
     if (allocated(qldt_auto)) deallocate (qldt_auto)
             allocate(qldt_auto(size(T,1),size(T,2),size(T,3)))
     if (allocated(qldt_fill)) deallocate (qldt_fill)
             allocate(qldt_fill(size(T,1),size(T,2),size(T,3)))
     if (allocated(qldt_destr)) deallocate (qldt_destr)
             allocate(qldt_destr(size(T,1),size(T,2),size(T,3)))
     if (allocated(rain_evap)) deallocate (rain_evap)
             allocate(rain_evap(size(T,1),size(T,2),size(T,3)))
     if (allocated(liq_adj)) deallocate (liq_adj)
             allocate(liq_adj(size(T,1),size(T,2),size(T,3)))
     if (allocated(qidt_dep)) deallocate (qidt_dep)
             allocate(qidt_dep(size(T,1),size(T,2),size(T,3)))
     if (allocated(qidt_eros)) deallocate (qidt_eros)
             allocate(qidt_eros(size(T,1),size(T,2),size(T,3)))
     if (allocated(qidt_fall)) deallocate (qidt_fall)
             allocate(qidt_fall(size(T,1),size(T,2),size(T,3)))
     if (allocated(qidt_fill)) deallocate (qidt_fill)
             allocate(qidt_fill(size(T,1),size(T,2),size(T,3)))
     if (allocated(qidt_subl)) deallocate (qidt_subl)
             allocate(qidt_subl(size(T,1),size(T,2),size(T,3)))
     if (allocated(qidt_melt)) deallocate (qidt_melt)
             allocate(qidt_melt(size(T,1),size(T,2),size(T,3)))
     if (allocated(qidt_destr)) deallocate (qidt_destr)
             allocate(qidt_destr(size(T,1),size(T,2),size(T,3)))
     if (allocated(snow_melt)) deallocate (snow_melt)
             allocate(snow_melt(size(T,1),size(T,2),size(T,3)))
     if (allocated(ice_adj)) deallocate (ice_adj)
             allocate(ice_adj(size(T,1),size(T,2),size(T,3)))
     if (allocated(snow_subl)) deallocate (snow_subl)
             allocate(snow_subl(size(T,1),size(T,2),size(T,3)))
     if (allocated(qadt_lsform)) deallocate (qadt_lsform)
             allocate(qadt_lsform(size(T,1),size(T,2),size(T,3)))
     if (allocated(qadt_lsdiss)) deallocate (qadt_lsdiss)
             allocate(qadt_lsdiss(size(T,1),size(T,2),size(T,3)))
     if (allocated(qadt_eros)) deallocate (qadt_eros)
             allocate(qadt_eros(size(T,1),size(T,2),size(T,3)))
     if (allocated(qadt_rhred)) deallocate (qadt_rhred)
             allocate(qadt_rhred(size(T,1),size(T,2),size(T,3)))
     if (allocated(qadt_destr)) deallocate (qadt_destr)
             allocate(qadt_destr(size(T,1),size(T,2),size(T,3)))
     if (allocated(qadt_fill)) deallocate (qadt_fill)
             allocate(qadt_fill(size(T,1),size(T,2),size(T,3)))
     if (allocated(qadt_super)) deallocate (qadt_super)
             allocate(qadt_super(size(T,1),size(T,2),size(T,3)))
     
     if (allocated(rain_cld_diag)) deallocate (rain_cld_diag)
             allocate(rain_cld_diag(size(T,1),size(T,2),size(T,3)+1))
     if (allocated(rain_clr_diag)) deallocate (rain_clr_diag)
             allocate(rain_clr_diag(size(T,1),size(T,2),size(T,3)+1))
     if (allocated(a_rain_cld_diag)) deallocate (a_rain_cld_diag)
             allocate(a_rain_cld_diag(size(T,1),size(T,2),size(T,3)+1))
     if (allocated(a_rain_clr_diag)) deallocate (a_rain_clr_diag)
             allocate(a_rain_clr_diag(size(T,1),size(T,2),size(T,3)+1))
     if (allocated(snow_cld_diag)) deallocate (snow_cld_diag)
             allocate(snow_cld_diag(size(T,1),size(T,2),size(T,3)+1))
     if (allocated(snow_clr_diag)) deallocate (snow_clr_diag)
             allocate(snow_clr_diag(size(T,1),size(T,2),size(T,3)+1))
     if (allocated(a_snow_cld_diag)) deallocate (a_snow_cld_diag)
             allocate(a_snow_cld_diag(size(T,1),size(T,2),size(T,3)+1))
     if (allocated(a_snow_clr_diag)) deallocate (a_snow_clr_diag)
             allocate(a_snow_clr_diag(size(T,1),size(T,2),size(T,3)+1))
     if (allocated(a_precip_cld_diag)) deallocate (a_precip_cld_diag)
             allocate(a_precip_cld_diag(size(T,1),size(T,2),size(T,3)+1))
     if (allocated(a_precip_clr_diag)) deallocate (a_precip_clr_diag)
             allocate(a_precip_clr_diag(size(T,1),size(T,2),size(T,3)+1))
     
     if (allocated(mask3)) deallocate (mask3)
             allocate(mask3(size(T,1),size(T,2),size(T,3)+1))
                         
     end if

!-----------------------------------------------------------------------
!
!       initialize select variables to zero. The variables reset
!       are:
!
!       (1) changes of prognostic variables
!
!       (2) variables dealing with the rain/snow fluxes. 
!
!       (3) qa_mean of the level above the top level.
!        
!       (4) diagnostic output fields

        ST = 0.
        SQ = 0.
        SL = 0.
        SI = 0.
        SA = 0.
        if (present(SN)) SN = 0.
   
        rain_cld   = 0.
        rain_clr   = 0.
        a_rain_cld = 0.
        a_rain_clr = 0.
        snow_cld   = 0.
        snow_clr   = 0.
        a_snow_clr = 0.
        a_snow_cld = 0.
        
        qa_mean_lst= 0.

        dcond_ls      = 0.
        dcond_ls_ice  = 0.
        qcg           = 0.
        qcg_ice       = 0.
        
        if (do_budget_diag) then
             if(do_liq_num) then 
               qndt_cond       = 0.
               qndt_evap       = 0.
             endif
             qldt_cond         = 0.
             qldt_evap         = 0.
             qldt_eros         = 0.
             qldt_accr         = 0.
             qldt_auto         = 0.
             qldt_berg         = 0.
             qldt_freez        = 0.
             qldt_rime         = 0.
             qldt_fill         = 0.
             qldt_destr        = 0.
             liq_adj           = 0.
             rain_evap         = 0.
             qidt_dep          = 0.
             qidt_subl         = 0.
             qidt_eros         = 0.
             qidt_fall         = 0.
             qidt_melt         = 0.
             qidt_destr        = 0.
             qidt_fill         = 0.
             snow_subl         = 0.
             snow_melt         = 0.
             ice_adj           = 0.
             rain_clr_diag     = 0.
             rain_cld_diag     = 0.
             a_rain_clr_diag   = 0.
             a_rain_cld_diag   = 0.
             a_precip_clr_diag = 0.
             a_precip_cld_diag = 0.
             snow_clr_diag     = 0.
             snow_cld_diag     = 0.
             a_snow_clr_diag   = 0.
             a_snow_cld_diag   = 0.
             qadt_lsform       = 0.
             qadt_lsdiss       = 0.
             qadt_rhred        = 0.
             qadt_eros         = 0.
             qadt_fill         = 0.
             qadt_super        = 0.
             qadt_destr        = 0.
        end if
        
        cloud_generator_on = do_cloud_generator()
                     
!-----------------------------------------------------------------------
!
!       Determine dimensions of slab

        idim = SIZE(T,1)
        jdim = SIZE(T,2)
        kdim = SIZE(T,3)
        
!-----------------------------------------------------------------------
!
!       compute inverse time step

        inv_dtcloud = 1.0 / dtcloud

!-----------------------------------------------------------------------
!
!       Calculate saturation specific humidity and its temperature 
!       derivative, and thermal conductivity plus vapor diffusivity
!       factor.
!
!       These are calculated according to the formulas:
!
!   (1)  qs   = d622*esat/ [pfull  -  (1.-d622)*esat]
!
!   (2) dqsdT = d622*pfull*(desat/dT)/[pfull-(1.-d622)*esat]**2.
!
!   (3) gamma = (L/cp) * dqsdT
!       
!       where d622 = rdgas/rvgas; esat = saturation vapor pressure;
!       and desat/dT is the temperature derivative of esat.
!       Note that in the calculation of gamma, 
!
!            {             hlv          for T > tfreeze             }
!       L =  { 0.05*(T-tfreeze+20.)*hlv + 0.05*(tfreeze-T)*hls      }
!            {                          for tfreeze-20.< T < tfreeze}
!            {             hls          for T < tfreeze-20.         }
!
!       This linear form is chosen because at tfreeze-20. es = esi, and
!       at tfreeze, es = esl, with linear interpolation in between.
!
!       The conductivity/diffusivity factor, A_plus_B is given by:
!
!   (4) A_plus_B =   { (hlv/Ka/T)*((hlv/rvgas/T)-1.) } + 
!
!                    { (rvgas*T/chi*esat) }
!
!       where Ka is the thermal conductivity of air = 0.024 J/m/s/K
!       and chi is the diffusitivy of water vapor in air which is
!       given by
!
!   (5) chi = 2.21 E-05 (m*m)/s  * (1.E+05)/pfull
!
!       where p is the pressure in Pascals.
!    
!
!       Note that qs, dqsdT, and gamma do not have their proper values
!       until all of the following code has been executed.  That
!       is qs and dqsdT are used to store intermediary results
!       in forming the full solution.

        !calculate water saturated vapor pressure from table
        !and store temporarily in the variable gamma
        if (do_pdf_clouds) then 
             call lookup_es(T-((hlv*ql+hls*qi)/cp_air),gamma)
        else
             call lookup_es(T,gamma)
        end if
                     
        !compute A_plus_B
        A_plus_B = ( (hlv/0.024/T) * ((hlv/rvgas/T)-1.) ) +            &
           (rvgas*T*pfull/2.21/gamma)  
        
        !calculate denominator in qsat formula
        qs = pfull - d378*gamma
     
        !limit denominator to esat, and thus qs to d622
        !this is done to avoid blow up in the upper stratosphere
        !where pfull ~ esat
        qs = max(qs,gamma) 
        
        !compute temperature derivative of esat 
        if (do_pdf_clouds) then
             call lookup_des(T-((hlv*ql+hls*qi)/cp_air),dqsdT)
        else
             call lookup_des(T,dqsdT)
        end if
        
        !calculate dqsdT
        dqsdT = d622*pfull*dqsdT/qs/qs

        !calculate qs
        qs=d622*gamma/qs
         
        !calculate gamma
        if (do_pdf_clouds) then
        gamma = dqsdT *(min(1.,max(0.,0.05*(T-((hlv*ql+hls*qi)/cp_air) &
                                             -tfreeze+20.)))*hlv +     &
                        min(1.,max(0.,0.05*(tfreeze -T+((hlv*ql+hls*qi)&
                                              /cp_air)   )))*hls)/cp_air
        else
        gamma = dqsdT *(min(1.,max(0.,0.05*(T-tfreeze+20.)))*hlv +     &
                        min(1.,max(0.,0.05*(tfreeze -T   )))*hls)/cp_air
        end if             
!-----------------------------------------------------------------------
!
!       Calculate air density

!       airdens = pfull / (rdgas * T * (1. + d608*qv  - ql - qi) )
        where (qrat .gt. 0.) 
             airdens = pfull / (rdgas * T *(1.+(d608*qv/qrat)-ql-qi) )
        elsewhere
             airdens = pfull / (rdgas * T * (1.   - ql - qi) )
        end where

!-----------------------------------------------------------------------
!
!       Assign cloud droplet number based on land or ocean point.
        
        N = N_land*LAND + N_ocean*(1.-LAND)
        if (do_liq_num) then
        ix = size(ql,1)
        jx = size(ql,2)
        kx = size(ql,3)
      
        N3D(:,:,:) = 0.
        concen_an(:,:,:) = 0.
        concen_na(:,:,:) = 0.
        concen(:,:,:) = 0.
        totalmass1(:,:,:,:) = 0.
        qn_mean(:,:) = 0.
        qn_upd(:,:) = 0.
        
        do k = 1,kx
          do j = 1,jx
            do i = 1,ix
              do s = 1,2
                concen(i,j,k)=concen(i,j,k)+ &
                Aerosol%aerosol(i,j,k,s)
              end do
      
      
              if(phalf(i,j,k)<1.) then
                thickness=log(phalf(i,j,k+1)/1.)* &
                8.314*T(i,j,k)/(9.8*0.02888)
              else
                thickness=log(phalf(i,j,k+1)/ &
                phalf(i,j,k))*8.314*T(i,j,k)/(9.8*0.02888)
              end if

! anthro. and natural sulfate concentration (ug so4/m3)
              concen_an(i,j,k)=0.7273*Aerosol%aerosol(i,j,k,1)/thickness*1.0e9
              concen_na(i,j,k)=0.7273*Aerosol%aerosol(i,j,k,2)/thickness*1.0e9

! NO. 1 Sulfate; NO. 2 Sea Salt; NO. 3 Organics
              totalmass1(i,j,k,1)=(Aerosol%aerosol(i,j,k,1)+Aerosol%aerosol(i,j,k,2))  &
              /thickness*1.0e9*1.0e-12
      
              totalmass1(i,j,k,2)=Aerosol%aerosol(i,j,k,5)  &
              /thickness*1.0e9*1.0e-12
      
              totalmass1(i,j,k,3)=Aerosol%aerosol(i,j,k,3)  &
              /thickness*1.0e9*1.0e-12
            end do
          end do
        end do
        endif

!-----------------------------------------------------------------------
!
!       Is a sub-vertical grid scale distribution going to be neededed?  
!       If yes, then do ppm fits
!
        if (do_pdf_clouds) then
        
        !initialize quantities
        do j = 1, kdim
             delp(:,:,j) = phalf(:,:,j+1)-phalf(:,:,j)
        enddo     
        qta4(1,:,:,:) = max(qmin,qv+ql+qi)
        qtqsa4(1,:,:,:) = qta4(1,:,:,:)-qs
        
        if (nsublevels.gt.1) then
            do id = 1,idim
                call ppm2m_sak(qta4(:,id,:,:),delp(id,:,:),kdim,kmap,1,jdim,0,kord)
                call ppm2m_sak(qtqsa4(:,id,:,:),delp(id,:,:),kdim,kmap,1,jdim,0,kord)
            enddo                
        else
            qta4(2,:,:,:) = qta4(1,:,:,:)
            qta4(3,:,:,:) = qta4(1,:,:,:)
            qta4(4,:,:,:) = 0.
            qtqsa4(2,:,:,:) = qtqsa4(1,:,:,:)
            qtqsa4(3,:,:,:) = qtqsa4(1,:,:,:)
            qtqsa4(4,:,:,:) = 0.   
        end if

        end if  !end of do_pdf_clouds section

!-----------------------------------------------------------------------
!
!       Enter the large loop over vertical levels.  Level 1 is the top
!       level of the column and level kdim is the bottom of the model.
!       If MASK is present, each column may not have kdim valid levels.
!
        rain3d = 0.
        snow3d = 0.

        DO j = 1, kdim

       
!-----------------------------------------------------------------------
!
!       Calculate pressure thickness of level and relative humidity 

        !calculate difference in pressure across the grid box divided
        !by gravity
        deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
         
        !calculate GRID box mean relative humidity 
        !       U = min(max(0.,qv(:,:,j)/qs(:,:,j)),1.)
 
        where (qrat(:,:,j) .gt. 0.)
             U = min(max(0.,(qv(:,:,j)/(qrat(:,:,j)*qs(:,:,j)))),1.)
        elsewhere
             U = 0.
        end where
        
!-----------------------------------------------------------------------
!
!       Account for the fact that other processes may have created
!       negative tracer or extremely small values of tracer fields.
!       The general reason for the extremely small values of the 
!       tracer fields is due to vertical diffusion, advection of 
!       condensate or cumulus induced subsidence (also a form of 
!       advection) of condensate.
!
!       In this step any values of the prognostic variables which are 
!       less than qmin are reset to zero, while conserving total 
!       moisture.
!
!       Note that this is done slightly different for the Tiedtke
!       cloud fraction than it is for pdf clouds. In the former, 
!       the filling requires that cloud liquid, cloud ice, and 
!       cloud fraction are greater than qmin. For PDF clouds, 
!       cloud fraction need not be considered since it is diagnosed
!       below from the PDF clouds
!

    if (.not. do_pdf_clouds) then

        where (qa(:,:,j) .le. qmin)
             SA(:,:,j)   = SA(:,:,j) - qa(:,:,j)
             qa_upd = 0.
        elsewhere
             qa_upd = qa(:,:,j)      
        end where
        
        if (do_budget_diag) then
             where (qa(:,:,j) .le. qmin)
                    qadt_fill(:,:,j)  =    -qa(:,:,j) * inv_dtcloud
             endwhere
        end if

!       Correct for qa > RH, which is not permitted under the 
!       assumption that the cloudy air is saturated and the temperature 
!       inside and outside of the cloud are about the same.
     
        if (do_budget_diag) then
             where (qa_upd .gt. U)
                    qadt_rhred(:,:,j) =   (qa_upd-U)  * inv_dtcloud
             endwhere
        end if
        
        where (qa_upd .gt. U)
             SA(:,:,j)   = SA(:,:,j) + U - qa_upd
             qa_upd = U      
        end where
        
        if (.not. do_liq_num) then
          where (ql(:,:,j) .le. qmin .or. qa(:,:,j) .le. qmin)
             SL(:,:,j)   = SL(:,:,j) - ql(:,:,j)
             SQ(:,:,j)   = SQ(:,:,j) + ql(:,:,j)
             ST(:,:,j)   = ST(:,:,j) - hlv*ql(:,:,j)/cp_air
             ql_upd = 0.
          elsewhere
             ql_upd = ql(:,:,j)
          end where
        else
          do k = 1,jx
            do i = 1,ix
              if (ql(i,k,j) .le. qmin .or. qa(i,k,j) .le. qmin .or. qn(i,k,j) .le. qmin) then
                SL(i,k,j)   = SL(i,k,j) - ql(i,k,j)
                SQ(i,k,j)   = SQ(i,k,j) + ql(i,k,j)
                ST(i,k,j)   = ST(i,k,j) - hlv*ql(i,k,j)/cp_air
                SN(i,k,j)   = SN(i,k,j) - qn(i,k,j)
                ql_upd(i,k) = 0.
                qn_upd(i,k) = 0.
                N3D(i,k,j)    = 0.
              else
                ql_upd(i,k) = ql(i,k,j)
                qn_upd(i,k) = qn(i,k,j)
                N3D(i,k,j) = qn(i,k,j)*airdens(i,k,j)*1.e-6
              endif
            end do
          end do
        endif

        if (do_budget_diag) then
             where (ql(:,:,j) .le. qmin .or. qa(:,:,j) .le. qmin)
                  qldt_fill(:,:,j) = -1.*ql(:,:,j)*inv_dtcloud
             endwhere
        end if

        where (qi(:,:,j) .le. qmin .or. qa(:,:,j) .le. qmin)
             SI(:,:,j)   = SI(:,:,j) - qi(:,:,j)
             SQ(:,:,j)   = SQ(:,:,j) + qi(:,:,j)
             ST(:,:,j)   = ST(:,:,j) - hls*qi(:,:,j)/cp_air
             qi_upd = 0.
        elsewhere
             qi_upd = qi(:,:,j)
        end where
        
        if (do_budget_diag) then
             where (qi(:,:,j) .le. qmin .or. qa(:,:,j) .le. qmin)
                  qidt_fill(:,:,j) = -1.*qi(:,:,j)*inv_dtcloud
             endwhere
        end if
        
    else
        
        where (ql(:,:,j) .le. qmin)
             SL(:,:,j)   = SL(:,:,j) - ql(:,:,j)
             SQ(:,:,j)   = SQ(:,:,j) + ql(:,:,j)
             ST(:,:,j)   = ST(:,:,j) - hlv*ql(:,:,j)/cp_air
             ql_upd = 0.
        elsewhere
             ql_upd = ql(:,:,j)
        end where

        if (do_budget_diag) then
             where (ql(:,:,j) .le. qmin)
                  qldt_fill(:,:,j) = -1.*ql(:,:,j)*inv_dtcloud
             endwhere
        end if

        where (qi(:,:,j) .le. qmin)
             SI(:,:,j)   = SI(:,:,j) - qi(:,:,j)
             SQ(:,:,j)   = SQ(:,:,j) + qi(:,:,j)
             ST(:,:,j)   = ST(:,:,j) - hls*qi(:,:,j)/cp_air
             qi_upd = 0.
        elsewhere
             qi_upd = qi(:,:,j)
        end where
        
        if (do_budget_diag) then
             where (qi(:,:,j) .le. qmin)
                  qidt_fill(:,:,j) = -1.*qi(:,:,j)*inv_dtcloud
             endwhere
        end if
        
    end if !for do_pdf_clouds
        
        
!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                 NON-CONVECTIVE CONDENSATION                          !
!                                                                      !
!                                                                      !
!                                                                      !
!                         METHOD 1                                     !
!                                                                      !
!                                                                      !
!                  TIEDTKE (1993) CLOUD FRACTION                       !
!                                                                      !
!       ANALYTIC INTEGRATION OF SATURATED VOLUME FRACTION EQUATION     !
!                                                                      !
!                                                                      !
!
!       Do non-convective condensation following Tiedtke, pages 3044-5.
!       In this formulation stratiform clouds are only formed/destroyed 
!       when there is upward or downward motion to support/destroy it. 
!
!       The first step is to compute the change in qs due to large-
!       scale processes, dqs_ls.   In Tiedtke, it has contributions from 
!       large-scale uplift, convection induced compensating subsidence,
!       turbulence cooling and radiative cooling.  dqs_ls has the form:
!
!               (((omega+ grav*Mc)/airdens/cp)+radturbten)*dqsdT*dtcloud
!   (6) dqs_ls= --------------------------------------------------------
!                  1.  +   ( qa +  (da_ls/2.) ) * gamma
!
!       Here da_ls is the increase in cloud fraction due to non-
!       convective processes.  Because this increase is also a function
!       of dqs_ls, a quadratic equation must be solved for dqs_ls in
!       the case that da_ls is not equal to zero.
!
!       Note that if the PDF cloud scheme is active the Tiedtke large-
!       scale condensation is bypassed.

    if (.not.do_pdf_clouds) then

        dqs_ls =(((omega(:,:,j)+grav*Mc(:,:,j))/airdens(:,:,j)/cp_air)+&
                 radturbten2(:,:,j))*dtcloud*dqsdT(:,:,j)

        !compute pressure dependent U00 following ECMWF formula if 
        !desired
        U00p = U00
        if (u00_profile) then
             where (pfull(:,:,j) .gt. 0.8*phalf(:,:,KDIM+1)) 
                    U00p = U00 + (1.-U00)* &
                         (((pfull(:,:,j)-(0.8*phalf(:,:,KDIM+1))) &
                                    /    (0.2*phalf(:,:,KDIM+1)) )**2.)
             end where
        end if       

!ljd
!       modify u00p to account for humidity in convective system
!       See "Tiedtke u00 adjustment" notes, 10/22/02
!ljd
 
        u00pr=u00p+(1.-u00p)*ahuco(:,:,j)
        u00p=u00pr

        where (dqs_ls .le. 0. .and. U .ge. U00p .and. qa_upd .lt. 1.)
             tmp1 = sqrt( (1.+qa_upd*gamma(:,:,j))**2. - (1.-qa_upd) * &
                    (1.-qa_upd)*gamma(:,:,j)*dqs_ls/qs(:,:,j)/         &
                    max(1.-U,qmin) ) - (1.+qa_upd*gamma(:,:,j))
             tmp1 = -1. * tmp1 / ((1.-qa_upd)*(1.-qa_upd)*gamma(:,:,j)/&
                    qs(:,:,j)/max(1.-U,qmin)/2.)
             dqs_ls = min(tmp1,dqs_ls/(1.+0.5*(1.+qa_upd)*gamma(:,:,j)))
        elsewhere
             dqs_ls = dqs_ls/(1.+qa_upd*gamma(:,:,j))
        endwhere
      
!       The next step is to compute the change in saturated volume
!       fraction due to non-convective condensation, da_ls.   This 
!       occurs in two conditions:
!
!       (a)  dqs_ls < 0. and U00 < U < 1., where U00 is the threshold
!            relative humidity for non-convective condensation. Note 
!            that if U is greater than or equal to 1., ideally qa = 1,
!            and da_ls = 0.  However this may not be the case for 
!            numerical reasons so this must be assured after analytic 
!            integration of the qa equation.
!
!            For these cases the change in saturated volume fraction is:
!
!   (7)      da_ls = - (1.-qa)*(1.-qa)*dqs_ls/2./qs/(1.-U)
!
!            This formula arises from the assumption that vapor is uni-
!            formly distributed in the range [qv_clr - (qs - qv_clr),qs]
!            where qv_clr is the amount of vapor in the unsaturated 
!            volume and is given from the following equation:
!
!   (8)      qv  =   qa * qs      +   (1.-qa) * qv_clr
!          
!            Implicit in equation (7) is the following assumption:
!            As qsat changes, the distribution of qv+ql+qi 
!            remains constant.  That is as qsat rises, portions where
!            qv+ql+qi > qsat+dqsat remain saturated.  This can only
!            occur if it is assumed that ql+qi evaporate-sublimate or
!            condense-deposit to keep qv = qsat. 
!
!       (b)  dqs_ls > 0.  Ideally some portion of the cloud should
!            evaporate however this is not accounted for at present.
!            

        !compute formula for da_ls
!       where (dqs_ls .le. 0. .and. U .ge. U00p)
!            da_ls = -0.5 * (1.-qa_upd) * (1.-qa_upd) * dqs_ls /       &
!             qs(:,:,j) / max(1.-U,qmin)
!       elsewhere
!            da_ls = 0.
!       end where 
 
        where ((dqs_ls .le. 0. .and. U .ge. U00p) .and. &
              (qa_upd+ahuco(:,:,j) .le. 1.))
             da_ls = -0.5 * (1.-qa_upd-ahuco(:,:,j)) * (1.-qa_upd-    &
                ahuco(:,:,j))    * dqs_ls/   qs(:,:,j) / max(1.-U,qmin)
        elsewhere
             da_ls = 0.
        end where 

!       Turbulent erosion of clouds
!
!       As in Tiedtke (1993) this is calculated using the eros_scale
!       parameter as:
!
!   (9) dql/dt    =  - qa * eros_scale * (qs - qv) * (ql/ ql+qi )
!
!  (10) dqi/dt    =  - qa * eros_scale * (qs - qv) * (qi/ ql+qi )
!
!  (11) dqa/dt    =  - qa * eros_scale * (qs - qv) * (qa/ ql+qi )
!
!       for which the erosion sink term (B in equation 13) is
!
!  (12) B = qa * eros_scale * (qs - qv) / (ql + qi)  
!
!
!       Theory for eros_scale
!
!       If eros_choice equals false, then a single erosion time scale
!       is used in all conditions (eros_scale).  If eros_choice equals
!       true then it is assumed that the timescale for turbulent 
!       evaporation is a function of the conditions in the grid box.  
!       Specifically, if the flow is highly turbulent then the scale is 
!       short, and eros_scale is large.  Likewise if convection is 
!       occurring, then it is assumed that the erosion term is larger 
!       than backround conditions. 
!
!       Here are the typical values for the timescales and the 
!       switches used:
!
!         Mixing type      eros_scale (sec-1)          Indicator
!       ----------------   ------------------     --------------------
!
!       Background            1.e-06              always present
!       Convective layers     5.e-06              Mc > Mc_thresh
!       Turbulent  layers     5.e-05              diff_t > diff_thresh
!

        !Background erosion scale
        tmp2 = eros_scale

        !Do enhanced erosion in convective or turbulent layers?
        !
        !              IMPORTANT NOTE
        !                
        !Note that convection is considered first, so that if 
        !turbulence and convection occur in the same layer, the
        !erosion rate for turbulence is selected.                
        !
        
        if (eros_choice) then
                
             !Enhanced erosion in convective layers
             where (Mc(:,:,j) .gt. mc_thresh)
                  tmp2 = eros_scale_c
             endwhere
        
             !Enhanced erosion in turbulent layers
             where ( (diff_t(:,:,j)            .gt.diff_thresh)  .or.  &
                     (diff_t(:,:,min(j+1,KDIM)).gt.diff_thresh) )
                  tmp2 = eros_scale_t
             endwhere
        
        end if   !for erosion choice

        where (ql_upd .gt. qmin .or. qi_upd .gt. qmin)
             D_eros = qa_upd * tmp2 * dtcloud * qs(:,:,j) *      &
                      (1.-U) / (qi_upd + ql_upd)
             where (pfull(:,:,j) .gt. 400.e02)
                 D_eros=D_eros+efact*D_eros*((pfull(:,:,kdim)-  &
                       pfull(:,:,j))/(pfull(:,:,kdim)-400.e02))
             elsewhere
                 D_eros=D_eros+efact*D_eros
             endwhere
        elsewhere
             D_eros = 0.
        end where
     
!    
!       The next step is to analytically integrate the saturated volume
!       fraction equation.  This follows the Tiedtke approach
!
!       The qa equation is written in the form:
!
!  (13) dqa/dt    =   (1.-qa) * A   -  qa * B 
!
!       Note that over the physics time step, A, B are assumed to be 
!       constants.
!
!       Defining qa(t) = qa0 and qa(t+dtcloud) = qa1, the analytic
!       solution of the above equation is:
!
!  (14) qa1 = qaeq -  (qaeq - qa0) * exp (-(A+B)*dtcloud)
! 
!       where qaeq is the equilibrium cloud fraction that is approached
!       with an time scale of 1/(A+B),
!
!  (15) qaeq  =  A/(A+B)
!
!
!       To diagnose the magnitude of each of the right hand terms of
!       (13) integrated over the time step, define the average cloud
!       fraction in the interval t to t + dtcloud qabar as:
!
!  (16) qabar  = qaeq - [ (qa1-qa0) / ( dtcloud * (A+B) ) ]
! 
!       from which the magnitudes of the A and B terms integrated
!       over the time step are:
!
!       A * (1-qabar)    and    -B * (qabar)
!
!       Additional notes on this analytic integration:
!
!       1.   For large-scale cloud formation or destruction from 
!            the dqs_ls term the contributions to A or B are defined
!            from:
!
!  (19)      A_ls * (1. - qa) = da_ls / dtcloud      if da_ls >= 0.
! 
!  (20)      B_ls * qa        = da_ls / dtcloud      if da_ls < 0.
!
!       2.   Note that to reduce the number of variables, the following
!            equivalency exists:
!
!               Ql or Qi equation              Qa equation
!             --------------------         -------------------
! 
!                     C_dt                        A_dt
!                     D_dt                        B_dt
!                     qceq                        qaeq
!                     qcbar                       qabar
!                     qc1                         qa1
!                     qc0                         qa0
!
!       3.   Qa goes to zero only in the case of ql and qi less than or
!            equal to qmin; see 'cloud destruction code' near the end of 
!            this loop over levels.
!
        
        !compute C_dt; This is assigned to the large-scale source term
        !following (18). Reset D_dt.
        C_dt = da_ls/max((1.-qa_upd),qmin)
        D_dt = D_eros
  
        !do analytic integration      
        where ( (C_dt.gt.Dmin) .or. (D_dt.gt.Dmin) ) 
             qc0   = qa_upd
             qceq  = C_dt  / (C_dt + D_dt)
             qc1   = qceq - (qceq - qc0) * exp ( -1.*(C_dt+D_dt) )
             qcbar = qceq - ((qc1 - qc0)/ (C_dt + D_dt))
        elsewhere
             qc0   = qa_upd
             qceq  = qc0   
             qc1   = qc0   
             qcbar = qc0  
        end where

        !set total tendency term and update cloud fraction    
        SA(:,:,j)  = SA(:,:,j) + qc1 - qc0
        qa_upd     = qc1
        
        if (do_budget_diag) then
             qadt_lsform(:,:,j) =  C_dt * (1.-qcbar) * inv_dtcloud 
             qadt_eros  (:,:,j) =  D_dt *     qcbar  * inv_dtcloud
        end if
        tmp5 = C_dt * (1.-qcbar)

!       The next step is to calculate the change in condensate
!       due to non-convective condensation, dcond_ls. Note that this is
!       not the final change but is used only to apportion condensate
!       change between phases. According to Tiedtke 1993 this takes the
!       form:
!
!  (21) dcond_ls = -1. * (qa +  0.5*da_ls) * dqs_ls
!
!       Here the 0.5*da_ls represents using a midpoint cloud fraction.
!       This is accomplished by using the variable qcbar.

        dcond_ls = -1. * qcbar * dqs_ls
            

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                 NON-CONVECTIVE CONDENSATION                          !
!                                                                      !
!                                                                      !
!                                                                      !
!                         METHOD 2                                     !
!                                                                      !
!                                                                      !
!                STATISTICAL CLOUD FRACTION                            !
!                                                                      !
!                                                                      !

    else                   !for doing PDF cloud scheme
    

        !set Tiedtke erosion term to zero            
        D_eros = 0.
        
        !compute pdf cloud fraction and condensate
        ! 
        !Note that the SYMMETRIC beta distribution is used here.
        !
        !
        ! Initialize grid-box mean values of cloud fraction (qag),
        ! cloud condensate(qcg), and clear sky water vapor (qvg)

        qcg = 0.
        qvg = 0.
        qag = 0.
        
        !Create loop over sub-levels within a grid box
        do ns = 1, nsublevels
        
             !calculate normalized vertical level
             ! 0. = top of gridbox
             ! 1. = bottom of gridbox
        
             pnorm =  (real(ns) - 0.5 )/real(nsublevels)
        
             !First step is to calculating the minimum (qtmin)
             !of the total water distribution and 
             !the width of the qt distribution (deltaQ)
             !
             !For diagnostic variance this is set to (1.-qthalfwidth)*qtbar
             !and 2*qthalfwidth*qtbar, respectively, where qtbar is the
             !mean total water in the grid box.        
             !
             !

             qtbar = qta4(2,:,:,j)+pnorm*( (qta4(3,:,:,j)-qta4(2,:,:,j)) + &
                                       qta4(4,:,:,j)*(1-pnorm) )
             
             qtbar = max(qmin,qtbar)
             deltaQ = 2.*qthalfwidth*qtbar
             qtmin = (1.-qthalfwidth)*qtbar
        
             !From this the variable normalized saturation specific
             !humidity qs_norm is calculated.
             !
             !  qs_norm = (qs(Tl) - qtmin)/(qtmax-qtmin)
             !
             !          = 0.5  - (qtbar - qs(Tl))/deltaQ
             !
             !Note that if qs_norm > 1., the grid box is fully clear.
             !If qs_norm < 0., the grid box is fully cloudy.
        
             qs_norm = qtqsa4(2,:,:,j)+  &
                       pnorm*( (qtqsa4(3,:,:,j)-qtqsa4(2,:,:,j)) + &
                       qtqsa4(4,:,:,j)*(1-pnorm) )
      
             qs_norm = 0.5 - ( qs_norm/deltaQ )
             
             !Calculation of cloud fraction (qagtmp), cloud condensate 
             !(qcgtmp), and water vapor in clear air part of the grid 
             !box (qvgtmp)
             !
             !Formulas (from Tompkins, and personal derivations):
             !
             !  Define icbp  = incomplete_beta(qs_norm,p,q)
             !         icbp1 = incomplete_beta(qs_norm,p+1,q)
             !
             !  qagtmp = 1. - icbp
             !
             !  qcgtmp = aThermo * {  (qtbar-qtmin)*(1.-icbp1) - 
             !                       qs_norm*deltaQ*(1.-icbp ) }
             !
             !
             !  qvgtmp = qtmin + (p/(p+q))*(icbp1/icbp)*deltaQ
             !
             !  
             ! where aThermo = 1./(1.+(L/cp)*dqsdT)
             !
             ! note that in the qvg formula below the factor of 0.5
             ! is equal to (p/(p+q)).
             !

             do jd = 1,jdim
             do id = 1,idim
        
             if (qs_norm(id,jd).le.1.) then
                 
                 icbp = incomplete_beta(max(0.,qs_norm(id,jd)), &
                                      p = betaP    , q = betaP)
                 icbp1= incomplete_beta(max(0.,qs_norm(id,jd)), &
                                      p = betaP + 1, q = betaP)
                 qagtmp(id,jd) = 1.-icbp
                 qcgtmp(id,jd) = (qtbar(id,jd)-qtmin(id,jd))*(1.-icbp1)&
                               - qs_norm(id,jd)*deltaQ(id,jd)*(1.-icbp)    
                 qcgtmp(id,jd) = qcgtmp(id,jd)/(1.+gamma(id,jd,j))
                 qvgtmp(id,jd) = qtmin(id,jd) + &
                               0.5*(icbp1/max(icbp,qmin))*deltaQ(id,jd)
             
                 !bound very very small cloud fractions which may
                 !cause negative cloud condensates due to roundoff 
                 !errors or similar errors in the beta table lookup.
                 if((qagtmp(id,jd).lt.0.).or.(qcgtmp(id,jd).le.0.))then
                      qagtmp(id,jd) = 0.
                      qcgtmp(id,jd) = 0.
                      qvgtmp(id,jd) = qtbar(id,jd)
                 end if
                 
             else             
                 qagtmp(id,jd) = 0.
                 qcgtmp(id,jd) = 0.
                 qvgtmp(id,jd) = qtbar(id,jd)             
             end if
             
             enddo
             enddo
    
             !sum vertically
             !
             !note special averaging of clear-sky water vapor
             !this is weighting clear-sky relative humidity by the 
             !clear-sky fraction
         
             qag = qag + qagtmp
             qcg = qcg + qcgtmp
             qvg = qvg+(1.-qagtmp)*min(max(qvgtmp/max(qmin, &
                       (qtbar+((qs_norm-0.5)*deltaQ))),0.),1.)
             
        enddo !for number of sublevels loop
        
        !compute grid-box average cloud fraction, cloud condensate
        !and water vapor
        
        if (nsublevels.gt.1) then
             qag = qag / real(nsublevels)
             qcg = qcg / real(nsublevels)
             
             !note special averaging of clear-sky water vapor
             do jd = 1,jdim
             do id = 1,idim
                  if ((1.-qag(id,jd)).gt.qmin) then
                      qvg(id,jd) =qvg(id,jd)/real(nsublevels)/&
                                  (1.-qag(id,jd))
                      qvg(id,jd) =qvg(id,jd) * qs(id,jd,j)
                  else
                      qvg(id,jd) = qs(id,jd,j)
                  end if
             enddo
             enddo
        else
             ! for nsublevels = 1, qag and qcg already hold their
             ! final values
             qvg = qvgtmp
        end if
             
        !do adjustment of cloud fraction
        qc0 = qa(:,:,j)
        qc1 = qag

        !set total tendency term and update cloud fraction    
        SA(:,:,j)  = SA(:,:,j) + qc1 - qc0
        qa_upd     = qc1

        if (do_budget_diag) then
             qadt_lsform(:,:,j) =  max(qc1-qc0,0.) * inv_dtcloud 
             qadt_lsdiss(:,:,j) =  max(qc0-qc1,0.) * inv_dtcloud
        end if

        !compute large-scale condensation / evaporation
        dcond_ls = qcg - (ql_upd + qi_upd)

    end if    !for doing PDF clouds

!       The next step is the apportionment on the non-convective 
!       condensation between liquid and ice phases. Following the
!       suggestion of Rostayn (2000), all condensation at temperatures
!       greater than -40C is in liquid form as ice nuclei are generally 
!       limited in the atmosphere. The droplets may subsequently be 
!       converted to ice by the Bergeron-Findeisan mechanism.  
!
!       One problem with this formulation is that the proper saturation
!       vapor pressure is not used for cold clouds as it should be
!       liquid saturation in the case of first forming liquid, but
!       change to ice saturation as the cloud glaciates.  The current
!       use of ice saturation beneath -20C thus crudely mimics the
!       result that nearly all stratiform clouds are glaciated for
!       temperatures less than -15C.
!
!       In the case of large-scale evaporation (dcond_ls<0.), it is
!       assumed that cloud liquid will evaporate faster than cloud
!       ice because if both are present in the same volume the
!       saturation vapor pressure over the droplet is higher than 
!       that over the ice crystal.
!
!       The fraction of large-scale condensation that is liquid
!       is stored in the temporary variable tmp1.   

        !assume liquid fractionation
        tmp1 = 1.
        tmp2 = 0.

        !For cases of cloud condensation where temperatures are
        !less than -40C create only ice
        where (dcond_ls .ge. 0. .and. T(:,:,j) .lt. tfreeze-40.)
             tmp1 = 0.
             tmp2 = 1.
        endwhere

        !For cases of cloud evaporation of mixed phase clouds
        !set liquid evaporation to preferentially occur first
        where ( (dcond_ls.lt.0.) .and. (ql_upd.gt.qmin)                &
                         .and. (qi_upd.gt.qmin) )   
             tmp1 = min(-1.*dcond_ls,ql_upd)/max(-1.*dcond_ls,qmin)
             tmp2 = 1.-tmp1
        end where

        !do evaporation of pure ice cloud
        where ( (dcond_ls.lt.0.) .and. (ql_upd.le.qmin)                &
                                 .and. (qi_upd.gt.qmin) )
             tmp1 = 0.
             tmp2 = 1.
        end where
        
        !calculate partitioning among liquid and ice to dcond_ls
        dcond_ls_ice = tmp2 * dcond_ls
        dcond_ls     = tmp1 * dcond_ls      

!       The next step is to compute semi-implicit qa,ql,qi which are 
!       used in many of the formulas below.  This gives a somewhat 
!       implicitness to the scheme. In this calculation an estimate 
!       is made of what the cloud fields would be in the absence of 
!       cloud microphysics and cloud erosion.
!
!       In the case of the Tiedtke cloud scheme, the mean cloud 
!       condensate is incremented if large-scale condensation is 
!       occuring. For cloud fraction, the value from the analytic 
!       integration above is used.
!
!       For the statistical cloud scheme these are set equal to the
!       values diagnosed from the beta-distribution apart from the
!       corrections for mixed phase clouds.

        qa_mean = qa_upd
        if (.not. do_pdf_clouds) then
             ql_mean = ql_upd + max(dcond_ls    ,0.)        
             qi_mean = qi_upd + max(dcond_ls_ice,0.)
        else
             ql_mean = max(ql_upd + dcond_ls    ,qmin)
             qi_mean = max(qi_upd + dcond_ls_ice,qmin)  
        end if
        if(do_liq_num) then 
!yim's CCN activation
          do k = 1,jx
            do i = 1,ix
              drop = 0.
              drop1(i,k) = 0.                
              if ( da_ls(i,k) > 0.0 ) then
                up_strat=-1.*(((omega(i,k,j)+grav*Mc(i,k,j))/airdens(i,k,j)/grav)+&
                        radturbten2(i,k,j)*cp_air/grav)

                totalmass(1) =totalmass1(i,k,j,1)
                totalmass(2) =0.1*totalmass1(i,k,j,2)
                totalmass(3) =1.67*totalmass1(i,k,j,3)
                          
                call aer_ccn_act(T(i,k,j),pfull(i,k,j),up_strat,totalmass,drop)
                drop1(i,k)=drop
              endif
            end do
          end do        
  
          qn_mean = qn_upd + max(tmp5,0.)*drop1*1.e6/airdens(:,:,j)
        endif


    
        !compute diagnostics for cloud fraction
        if (id_aall > 0) areaall(:,:,j) = qa_mean
        if (id_aliq > 0 .or. id_rvolume > 0) then
              where (ql_mean .gt. qmin) arealiq(:,:,j) = qa_mean
        end if
        if (id_aice > 0 .or. id_vfall > 0) then
              where (qi_mean .gt. qmin) areaice(:,:,j) = qa_mean
        end if              

  
!-----                                                            -----! 
!                                                                      !
!                                  END OF                              !
!                                                                      !
!                        NON-CONVECTIVE CONDENSATION                   !
!                                                                      !
!                                 SECTION                              !
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!
!
!            TRANSFER OF RAIN FLUXES AND SNOW FLUXES BETWEEN 
!      
!           SATURATED AND UNSATURATED PORTIONS OF THE GRID BOX
!
!                           AT LAYER INTERFACES
!                           
!
!       Do transfer of rain and snow fluxes between clear and cloud 
!       portions of the grid box. This formulism follows the rain/snow
!       parameterization developed by Jakob and Klein (2000).  It 
!       treats the grid mean flux of rain and snow separately according 
!       to the portion that exists in unsaturated air and the portion 
!       that exists in saturated air.  At the top of each level, some 
!       precipitation that was in unsaturated air in the levels above 
!       may enter saturated air and vice versa. These transfers are 
!       calculated here under an assumption for the overlap between 
!       precipitation area and saturated area at the same and adjacent 
!       levels.
!
!       For rain, the area of the grid box in which rain in saturated 
!       air of the level above enters unsaturated air in the current 
!       level is called da_cld2clr and is given by 
!       maximum(a_rain_cld - qa , 0.) in the case of maximum overlap of
!       rain_cld and qa, but equal to a_rain_cld*(1.-qa) in the case of
!       random overlap of cloud and precipitation. The grid mean flux of 
!       rain which is transfered from saturated air to unsaturated air 
!       is given by rain_cld*da_cld2clr/a_rain_cld.
!
!       The area of the grid box where rain in unsaturated air of the
!       level above enters saturated air in the current level is
!       called da_clr2cld and is given by
!       maximum(0.,mininum(a_rain_clr,qa(j)-qa(j-1))) in the case of
!       maximum overlap of cloud and rain_clr, but equal to a_rain_clr*
!       qa for the case of random overlap of cloud and precipitation.  
!       The grid mean flux of rain transfered from unsaturated air to 
!       saturated air is given by rain_clr*da_clr2cld/a_rain_clr.
!
!       NOTE: the overlap assumption used is set by the namelist 
!             variable in cloud_rad
!
!       Overlap values are:      1  = maximum
!                                2  = random
!
!       If cloud_generator is on, the overlap choice is taken from
!       there, by computing a quantity alpha, which is weighting 
!       between maximum and random overlap solutions.
!
!                alpha = 1 ---> maximum,
!                alpha = 0 ---> random
!
!       alpha is stored in tmp3.
!       tmp1 has the maximum overlap solution
!       tmp2 has the random overlap solution
       
!
!
!       Rain transfers are done first
!
!

        !-------------------------------
        !compute cloud to clear transfer
        if (overlap .eq. 1)                                            &
             da_cld2clr= min(a_rain_cld,max(0.,a_rain_cld - qa_mean)   )
        
        if (overlap .eq. 2)                                            &
             da_cld2clr= min(a_rain_cld,max(0.,a_rain_cld*(1.-qa_mean)))
       
        if (cloud_generator_on) then
             tmp3 = 0.0
             if (j.gt.1) then
                  tmp3 = compute_overlap_weighting(qa_mean_lst,qa_mean,&
                         pfull(:,:,j-1),pfull(:,:,j))
                  tmp3 = min(1.,max(0.,tmp3))       
             end if
             tmp1 =      min(a_rain_cld,max(0.,a_rain_cld - qa_mean)   )             
             tmp2 =      min(a_rain_cld,max(0.,a_rain_cld*(1.-qa_mean)))
             da_cld2clr=min(a_rain_cld,max(0.,tmp3*tmp1+(1.-tmp3)*tmp2))             
        end if      

        !-------------------------------
        !compute clear to cloud transfer
        if (overlap .eq. 1)                                            &
             da_clr2cld = min(max(qa_mean-qa_mean_lst,0.),a_rain_clr)
        if (overlap .eq. 2)                                            &
             da_clr2cld = min(max( a_rain_clr*qa_mean,0.),a_rain_clr)

        if (cloud_generator_on) then
             tmp1 =       min(max(qa_mean-qa_mean_lst,0.),a_rain_clr)             
             tmp2 =       min(max( a_rain_clr*qa_mean,0.),a_rain_clr)
             da_clr2cld=min(a_rain_clr,max(0.,tmp3*tmp1+(1.-tmp3)*tmp2))
        end if      
        
        !---------------------------------
        !calculate precipitation transfers
        dprec_cld2clr = rain_cld*(da_cld2clr/max(a_rain_cld,qmin))
        dprec_clr2cld = rain_clr*(da_clr2cld/max(a_rain_clr,qmin))
        
        !----------------
        !add in transfers
        a_rain_clr = a_rain_clr + da_cld2clr - da_clr2cld
        a_rain_cld = a_rain_cld - da_cld2clr + da_clr2cld
        rain_clr   = rain_clr + dprec_cld2clr - dprec_clr2cld        
        rain_cld   = rain_cld - dprec_cld2clr + dprec_clr2cld

!
!
!       Snow transfers are done second, in a manner exactly like that
!       done for the rain fluxes
!
        

        !-------------------------------
        !compute cloud to clear transfer
        if (overlap .eq. 1)                                            &
             da_cld2clr= min(a_snow_cld,max(0.,a_snow_cld-qa_mean)     ) 
        if (overlap .eq. 2)                                            &
             da_cld2clr= min(a_snow_cld,max(0.,a_snow_cld*(1.-qa_mean)))
      
        if (cloud_generator_on) then
             tmp1 =      min(a_snow_cld,max(0.,a_snow_cld-qa_mean)     )             
             tmp2 =      min(a_snow_cld,max(0.,a_snow_cld*(1.-qa_mean)))
             da_cld2clr=min(a_snow_cld,max(0.,tmp3*tmp1+(1.-tmp3)*tmp2))
        end if      

        !-------------------------------
        !compute clear to cloud transfer
        if (overlap .eq. 1)                                            &
             da_clr2cld = min(max(qa_mean-qa_mean_lst,0.),a_snow_clr)
        if (overlap .eq. 2)                                            &
             da_clr2cld = min(max(a_snow_clr*qa_mean,0.) ,a_snow_clr)

        if (cloud_generator_on) then
             tmp1 =       min(max(qa_mean-qa_mean_lst,0.),a_snow_clr)             
             tmp2 =       min(max(a_snow_clr*qa_mean,0.) ,a_snow_clr)
             da_clr2cld=min(a_snow_clr,max(0.,tmp3*tmp1+(1.-tmp3)*tmp2))
        end if      
        
        !---------------------------------
        !calculate precipitation transfers
        dprec_cld2clr = snow_cld*(da_cld2clr/max(a_snow_cld,qmin))
        dprec_clr2cld = snow_clr*(da_clr2cld/max(a_snow_clr,qmin))
        
        !----------------
        !add in transfers
        a_snow_clr = a_snow_clr + da_cld2clr - da_clr2cld
        a_snow_cld = a_snow_cld - da_cld2clr + da_clr2cld
        snow_clr = snow_clr + dprec_cld2clr - dprec_clr2cld
        snow_cld = snow_cld - dprec_cld2clr + dprec_clr2cld
   
               
!-----------------------------------------------------------------------
!
!
!                        MELTING OF CLEAR SKY SNOW FLUX
!
!
!       Melting of falling ice to rain occurs when T > tfreeze. The 
!       amount of melting is limited to the melted amount that would 
!       cool the temperature to tfreeze.
!
!       In the snowmelt bug version, the temperature of melting was 
!       tfreeze + 2. like the original Tiedtke (1993) paper, instead of 
!       tfreeze.

        !compute grid mean change in snow flux to cool the
        !grid box to tfreeze and store in temporary variable tmp1
        if (do_old_snowmelt) then
             tmp1 = cp_air*(T(:,:,j)-tfreeze-2.)*deltpg*inv_dtcloud/hlf
        else
             tmp1 = cp_air*(T(:,:,j)-tfreeze)*deltpg*inv_dtcloud/hlf
        end if
        
        ! If snow_clr > tmp1, then the amount of snow melted is
        ! limited to tmp1, otherwise melt snow_clr.  The amount
        ! melted is stored in tmp2
        tmp2 = max(min(snow_clr,tmp1),0.)     

        ST(:,:,j) = ST(:,:,j) - hlf*tmp2*dtcloud/deltpg/cp_air                
        rain_clr  = rain_clr + tmp2
        
        !raise a_rain_clr to a_snow_clr IF AND only IF melting occurs
        !and a_rain_clr < a_snow_clr
        where (tmp2 .gt. 0. .and. a_snow_clr .gt. qmin)
             a_rain_clr = max(a_rain_clr,a_snow_clr)
        end where

        ! If all of the snow has melted, then zero out a_snow_clr
        where (snow_clr.lt.tmp1 .and. a_snow_clr.gt.qmin)
             snow_clr = 0.
             a_snow_clr = 0.
        elsewhere
             snow_clr = snow_clr - tmp2          
        end where

        if (do_budget_diag) snow_melt(:,:,j) = tmp2/deltpg             
             
!-----------------------------------------------------------------------
!
!
!                        MELTING OF CLOUDY SKY SNOW FLUX
!
!
!       Melting of falling ice to rain occurs when T > tfreeze. The 
!       amount of melting is limited to the melted amount that would 
!       cool the temperature to tfreeze.
!

        if (.not.do_old_snowmelt) then

        !compute grid mean change in snow flux to cool the
        !grid box to tfreeze and store in temporary variable tmp1
        !
        !note that tmp1 already has the value of this variable 
        !from the clear-sky melt calculation, so one does not need
        !to repeat the calculation here.
        !
        !However, note that clear-sky snow melt may have already 
        !reduced the temperature of the grid box - this snow melt is in 
        !variable tmp2 from lines above. Thus the amount that one
        !can melt is less.
        
        tmp1 = tmp1 - tmp2
        
        ! If snow_cld > tmp1, then the amount of snow melted is
        ! limited to tmp1, otherwise melt snow_cld.  The amount
        ! melted is stored in tmp2
        tmp2 = max(min(snow_cld,tmp1),0.)     

        ST(:,:,j) = ST(:,:,j) - hlf*tmp2*dtcloud/deltpg/cp_air                
        rain_cld  = rain_cld + tmp2
        
        !raise a_rain_cld to a_snow_cld IF AND only IF melting occurs
        !and a_rain_cld < a_snow_cld
        where (tmp2 .gt. 0. .and. a_snow_cld .gt. qmin)
             a_rain_cld = max(a_rain_cld,a_snow_cld)
        end where

        ! If all of the snow has melted, then zero out a_snow_cld
        where (snow_cld.lt.tmp1 .and. a_snow_cld.gt.qmin)
             snow_cld = 0.
             a_snow_cld = 0.
        elsewhere
             snow_cld = snow_cld - tmp2          
        end where

        if (do_budget_diag) snow_melt(:,:,j) =  snow_melt(:,:,j) + &
                                                tmp2/deltpg             

        end if  !for snowmelt bugfix
                            
!----------------------------------------------------------------------!
!
!
!              COMPUTE SLOPE FACTOR FOR ICE MICROPHYSICS                  
!
!       [The following microphysics follows that of Rotstayn (1997)]
!       The number concentration of ice crystals of diameter D in the
!       SIZE interval D to D+dD is assumed to be 
!       distributed as in a Marshall Palmer distribution :
!
!  (22) N(D)dD = Nof * Exp( - lamda_f * D)
!
!       The slope factor and intercept are not assumed to be constant,
!       but the slope factor is
!
!  (23) lamda_f = 1.6X10^(3.+0.023(tfreeze-T))
!
!       Integration of (22) over all particle sizes with a constant
!       density of ice crystals , rho_ice, and assumed spherical shape
!       yields a relationship between the intercept parameter and qi
!
!  (24) Nof = airdens*qi_local*(lamda_f^4)/pi*rho_ice
!       
!
!       For the calculation of riming and sublimation of snow, 
!       lamda_f is needed, so it is calculated here.
!
!       Also qi_mean is updated here with the flux of ice that falls 
!       into the cloudy portion of the grid box from above. This permits
!       the Bergeron and riming process to know about the ice that falls
!       into the grid box in the same time step.

        !Increment qi_mean by the ice flux entering the
        !the grid box. To convert ice_flux to units of condensate by
        !multiply by dtcloud and dividing by the mass per unit area 
        !of the grid box. Implicit here is the assumption that the
        !ice entering the cloud will be spread instantaneously over
        !all of the cloudy area.
        qi_mean = qi_mean + snow_cld*dtcloud/deltpg        

        !snow falling into cloud reduces the amount that
        !falls out of cloud: a loss of cloud ice from settling
        !is defined to be positive
        if (do_budget_diag) qidt_fall(:,:,j)= -1.*snow_cld/deltpg
         
        !compute lamda_f
        lamda_f = 1.6 * 10**(3.+0.023*(tfreeze-T(:,:,j)))
        

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                       LIQUID PHASE MICROPHYSICS                      !
!                                                                      !
!                                 AND                                  !
!                                                                      !
!                  ANALYTIC INTEGRATION OF QL EQUATION                 !
!                                                                      !
!                                                                      !
!                                                                      !
!       Accretion
!
!       The parameterization of collection of cloud liquid drops by
!       rain drops follows the parameterization of Rotstayn (1997).
!       The parameterization starts with the continous-collection 
!       equation of drops by a rain drop of diameter D, falling with
!       speed V(D).   The fall speed of rain drops is taken from
!       Gunn and Kinzer(1949):
!
!  (25) V(D) = 141.4 (m**0.5,s-1) * sqrt(D) * (rho_ref/airdens)**0.5
!
!       where D is the radius of the rain drops in meters. Here
!       rho_ref is a reference air density.  This formula is generally
!       good for 1 mm < D < 4 mm.
!
!       The distribution of rain drops by SIZE follows a Marshall-
!       Palmer distribution:
!
!  (26) N(D) dD = Nor * Exp (- lamda *D)
!
!       where N(D)dD is the number of rain drops with SIZE D in the
!       interval D to D+dD, Nor is the intercept (assumed fixed at
!       8E+04 (1/m*m*m*m).  lamda is the slope intercept parameter
!       and with (21) it can be shown that lamda is a function of
!       the rain rate.
!
!       With these assumptions the local rate of accretion of cloud
!       liquid reduces to:
!
!  (27) dl/dt_local = - CB*Eco*((rain_rate_local/dens_h2o)**(7/9))*
!                        
!                       (rho_ref/airdens)**(1/9)   * ql_local
!
!       where CB is the accretion constant:
! 
!       CB = 65.772565 [(m)**(-7/9)] * [(s)**(-2/9)]
!
!       AND Eco is the collection efficiency of a cloud droplet by a 
!       rain droplet.   A good fit to the Table 8.2 of Rogers and Yau 
!       (1988) for rain drops between SIZE 1mm and 3mm is:
!
!  (28) Eco = rad_liq**2 / (rad_liq**2 + 20.5 microns**2) .
!
!       In generalizing to grid mean conditions (27) becomes:
!
!  (29) dl/dt = - (arain_cld/qa_mean) * CB * Eco * 
!
!                 [(rain_cld/a_rain_cld/dens_h2o)**(7/9)] * ql
!        
!       Note that the very weak dependence on air density is
!       neglected at this point.
!
!       The particle sizes are computed from the following equation
!
!  (30) rad_liq = (3*airdens*ql/(4*pi*liq_dens*qa*N)^(1/3)
!
!       
!       For numerical treatment we write (25) as:
!
!       dl/dt = - D_acc * l 
!
!       and if we do so:
!
!  (31) D_acc   =  (arain_cld/qa_mean) * CB * Eco * 
!
!                 [(rain_cld/a_rain_cld/dens_h2o)**(7/9)] 
!
!       In the work below, D_acc is added to D1_dt, the first processes
!       contributing to the depletion of cloud liquid in the analytic
!       integration.  D1_dt represents the conversion of liquid to rain.

        if (.not. do_liq_num) then
        !compute rad_liq.  The constant below is equal to  
        !1.E+06 * (3/4*pi)^(1/3), where the 1E+06 is
        !the factor to convert meters to microns.
        
          rad_liq= 620350.49 *( (airdens(:,:,j)*ql_mean/ &
                        max(qa_mean,qmin)/N/dens_h2o)**(1./3.))

        !do not let very small cloud fractions contribution to
        !autoconversion or accretion
          where (qa_mean .le. qmin) rad_liq = 0.
        else        
!yim The 1st place droplet number is used
          rad_liq= 620350.49 *( (ql_mean/max(qn_mean,qmin)/dens_h2o)**(1./3.))
          !do not let very small cloud fractions contribution to
          !autoconversion or accretion
          where (qa_mean .le. qmin .or. qn_upd .le.qmin) rad_liq = 0.
        endif
        
        if (id_rvolume > 0) rvolume(:,:,j) = rad_liq*arealiq(:,:,j)

        !compute accretion D term
        D1_dt =  dtcloud * 65.772565 * (a_rain_cld/max(qa_mean,qmin))* &
                 ( rad_liq*rad_liq / (rad_liq*rad_liq+20.5) ) *        &
                 ((rain_cld/max(a_rain_cld,qmin)/dens_h2o)**(7./9.))
            
        if (do_budget_diag) qldt_accr(:,:,j) = D1_dt
    
!       Autoconversion
!
!       The autoconversion parameterization follow that of Manton
!       and Cotton (1977).  This formula has been used in Chen and
!       Cotton (1987) and is used in the CSIRO GCM (Rotstayn 1997)
!       and the LMD GCM (Boucher, Le Treut and Baker 1995).  In this
!       formulation the time rate of change of grid mean liquid is
!
!  (32) dl/dt= -CA * qa * [(ql/qa)^(7/3)] * [(N*dens_h2o)^(-1/3)] 
!
!               * H(rad_liq - rthresh)
!
!       where N is the number of cloud droplets per cubic metre,
!       rthresh is a particle radius threshold needed to for autoconv-
!       ersion to occur, H is the Heaviside function, and CA is
!       a constant which is:
!
!  (33) CA =  0.104 * grav * Ec * (airdens)^(4/3) / mu 
!        
!       where grav is gravitational acceleration, Ec is the collection
!       efficiency, airdens is the density of air, and mu is the 
!       dynamic viscosity of air.   This constant is evaluated 
!       ignoring the temperature dependence of mu and with a fixed 
!       airdens of 1 kg/m3.
!
!       With   Ec = 0.55        (standard choice - see references)
!            grav = 9.81        m/(s*s)
!         airdens = 1.00        kg air/(m*m*m)
!              mu = 1.717  E-05 kg condensate/(m*s) 
!
!              CA = 32681. [(kg air)^(4/3)]/kg liq/m/s
!
!
!       For numerical treatment we write (32) as:
!
!       dl/dt = - D_aut * l 
!
!       and if we do so:
!
!  (34) D_aut   =   CA * [(N*dens_h2o)^(-1/3)] * [(ql/qa)^(4/3)] * &
!                   H(r-rthresh)
!
!       In the work below, D_aut is temporarily stored in the variable
!       tmp1 before being added to D1_dt.  D1_dt represents the 
!       conversion of liquid to rain.
!
!       Following Rotstayn, autoconversion is limited to the amount that
!       would reduce the local liquid cloud condensate to the critical
!       value at which autoconversion begins. This limiter is likely to
!       be invoked frequently and is computed from
!
!  (35) D_dt = log( (rad_liq/rthresh)**3. )
!
!       This limiter is stored in tmp2.
!
!
!
!       -------------------------------------------------
!
!       Khairoutdinov and Kogan (2000) Autoconversion
!
!       Reference: Khairoutdinov, M. and Y. Kogan, 2000: A new cloud 
!                  physics parameterization in a large-eddy simulation
!                  model of marine stratocumulus. Mon. Wea. Rev., 128,
!                  229-243.
!
!       If the namelist parameter use_kk_auto = true, then the 
!       Khairoutdinov and Kogan (KK) autoconversion parameterization
!       is used in place of the Manton and Cotton formula described
!       above.
!
!       In SI units this formula is:
!
!  (32A) dl/dt= -CA * qa * [(ql/qa)^(2.47)] * [(N)^(-1.79)] 
!
!
!       where N is the number of cloud droplets per cubic metre
!       and CA is a constant which is:
!
!  (33A) CA =  7.4188E+13 (kg condensate/kg air)**(-1.47)
!                         (# drops/meters**3)**(1.79)
!                         seconds**(-1) 
!        
!       For numerical treatment we write (32A) as:
!
!       dl/dt = - D_aut * l 
!
!       and if we do so:
!
!  (34A) D_aut   =   CA * [(N)^(-1.79)] * [(ql/qa)^(1.47)]
!

        if (.not. use_kk_auto ) then

             !compute autoconversion sink as in (34)
             tmp1 = 32681. * dtcloud * ((N*dens_h2o)**(-1./3.))*       &
                    ((ql_mean/max(qa_mean,qmin))**(4./3.))
        
             !compute limiter as in (35)
             tmp2 =max(3*log(max(rad_liq,qmin)/rthresh),0.)

             !limit autoconversion to the limiter
             tmp1 = min(tmp1,tmp2)
        
        else
        
          if (do_liq_num) then 
!*************************yim's version based on Khairoutdinov and Kogan (2000)
!The second place N is used
          tmp1 = dtcloud * 1350. *  &
                 (1.e-6*max(qn_mean,N_min)*airdens(:,:,j))**(-1.79)*  &
                 (ql_mean)**(1.47)*max(qa_mean,qmin)**(0.32)
!**************************
          else

             tmp1 = 0.0
             
             !compute autoconversion sink as in (34A)
             where (ql_mean.gt.qmin)
                  tmp1 = 7.4188E+13 * dtcloud *  (N**(-1.79))*         &
                    ((ql_mean/max(qa_mean,qmin))**(1.47))
             endwhere
               
          endif
        end if

        !add autoconversion to D1_dt
        D1_dt = D1_dt + tmp1

        !auto conversion will change a_rain_cld upto area of cloud
        where (tmp1 .gt. Dmin) a_rain_cld = qa_mean

        if (do_budget_diag) qldt_auto(:,:,j) = tmp1        

        if ( id_autocv > 0 ) then
             where ( rad_liq .gt. rthresh ) areaautocv(:,:,j) = qa_mean       
        end if
        

!       Bergeron-Findeisan Process 
!
!       where ice and liquid coexist, the differential saturation
!       vapor pressure between liquid and ice phases encourages
!       the growth of ice crystals at the expense of liquid droplets.
!
!       Rotstayn (2000) derive an equation for the growth of ice by
!       starting with the vapor deposition equation for an ice crystal
!       and write it in terms of ice specific humidity as:
!
!                 {(Ni/airdens)**2/3}*7.8
!  (36) dqi/dt =  ------------------------  X [(esl-esi)/esi] X
!                 [rhoice**1/3]* A_plus_B
!
!                 ((max(qi,Mio*Ni/airdens))**(1/3))*
!
!       Here Ni is the ice crystal number which is taken from the 
!       parameterization of Meyers et al. :
!
!  (37) Ni = 1000 * exp( (12.96* [(esl-esi)/esi]) - 0.639 )
!
!       The use of the maximum operator assumed that there is a 
!       background ice crystal always present on which deposition can 
!       occur.  Mio is an initial ice crystal mass taken to be 10-12kg.
!
!       Figure 9.3 of Rogers and Yau (1998) shows the nearly linear
!       variation of [(esl-esi)/esi] from 0. at 273.16K to 0.5 at 
!       233.16K.  Analytically this is parameterized as (tfreeze-T)/80.
!
!
!       Generalizing (36) to grid mean conditions and writing it in 
!       terms of a loss of cloud liquid yields:
!
!                  (1000*exp((12.96*(tfreeze-T)/80)-0.639)/airdens)**2/3
!  (38) dql/dt = - ----------------------------------------------------- 
!                           [rhoice**1/3]* A_plus_B * ql
!
!       *qa*7.8*((max(qi/qa,Mio*Ni/airdens))**(1/3))*[(tfreeze-T)/80]*ql
!
!       Note that the density of ice is set to 700 kg m3 the value 
!       appropriate for pristine ice crystals.  This value is 
!       necessarily different than the value of ice used in the riming 
!       and sublimation part of the code which is for larger particles 
!       which have much lower densities.
        
        !do Bergeron process
        where ( (T(:,:,j) .lt. tfreeze) .and. (ql_mean .gt. qmin)      &
                                        .and. (qa_mean .gt. qmin))           
             D2_dt =  dtcloud * qa_mean * ((1000.*exp((12.96*0.0125*   &
                      (tfreeze-T(:,:,j)))-0.639)/airdens(:,:,j))**(2./ &
                      3.))* 7.8* ((max(qi_mean/qa_mean,1.E-12*1000.*   &
                      exp((12.96*0.0125*(tfreeze-T(:,:,j)))-0.639)     &
                      /airdens(:,:,j)))**(1./3.))*0.0125*              &
                      (tfreeze-T(:,:,j))/((700.**(1./3.))*             &
                      A_plus_B(:,:,j)*ql_mean)
        elsewhere
             D2_dt = 0.0        
        end where
       
        if (do_budget_diag) qldt_berg(:,:,j) = D2_dt
       
!       Accretion of cloud liquid by ice ('Riming')
!       
!       [The below follows Rotstayn]
!       Accretion of cloud liquid by ice ('Riming') is derived in
!       the same way as the accretion of cloud liquid by rain. That
!       is the continous-collection equation for the growth of an
!       ice crystal is integrated over all ice crystal sizes. This
!       calculation assumes all crystals fall at the mass weighted
!       fall speed, Vfall.  This yields the following equation after
!       accounting for the area of the interaction
!
!  (39) dql/dt = -  ( a_snow_cld / qa/a_snow_cld ) *
!                   ( ELI * lamda_f * snow_cld / 2/ rho_ice ) * ql 
!
!
!       Note that in the version with the snowmelt bug, riming was
!       prevented when temperatures were in excess of freezing.

        !add in accretion of cloud liquid by ice
        tmp1 = 0.0
        if (do_old_snowmelt) then
             where ((a_snow_cld.gt.qmin) .and. (ql_mean.gt.qmin) .and. &
                    (   qa_mean.gt.qmin) .and. (T(:,:,j) .lt. tfreeze) )            
                 tmp1 = dtcloud*0.5*ELI*lamda_f*snow_cld/qa_mean/rho_ice              
             end where
        else
             where ((a_snow_cld.gt.qmin) .and. (ql_mean.gt.qmin) .and. &
                    (   qa_mean.gt.qmin) )            
                 tmp1 = dtcloud*0.5*ELI*lamda_f*snow_cld/qa_mean/rho_ice              
             end where
        end if        
        
        D2_dt = D2_dt + tmp1

        if (do_budget_diag) qldt_rime(:,:,j) = tmp1

!       Freezing of cloud liquid to cloud ice occurs when
!       the temperature is less than -40C. At these very cold temper-
!       atures it is assumed that homogenous freezing of cloud liquid
!       droplets will occur.   To accomplish this numerically in one 
!       time step:
!
!  (40) D*dtcloud =  ln( ql / qmin ).
!
!       With this form it is guaranteed that if this is the only
!       process acting that ql = qmin after one integration.
!
               
        !do homogeneous freezing
        where ( T(:,:,j).lt.tfreeze-40..and.(ql_mean.gt.qmin).and.     &
               (qa_mean.gt.qmin))
             D2_dt = log ( ql_mean / qmin )
        end where
        
        if (do_budget_diag) then
             where (T(:,:,j).lt.(tfreeze-40.).and.(ql_mean.gt.qmin)    &
                  .and.(qa_mean.gt.qmin))
                  qldt_freez(:,:,j) = D2_dt
                  qldt_rime (:,:,j) = 0.
                  qldt_berg (:,:,j) = 0.     
             end where
        end if
  
!       Analytic integration of ql equation
!
!
!       The next step is to analytically integrate the cloud liquid
!       condensate equation.  This follows the Tiedtke approach.
!
!       The qc equation is written in the form:
!
!  (41) dqc/dt    =   C   -  qc * D   
!
!       Note that over the physics time step, C and D are assumed to 
!       be constants.
!
!       Defining qc(t) = qc0 and qc(t+dtcloud) = qc1, the analytic
!       solution of the above equation is:
!
!  (42) qc1 = qceq -  (qceq - qc0) * exp (-D*dtcloud)
! 
!       where qceq is the equilibrium cloud condensate that is approached
!       with an time scale of 1/D,
!
!  (43) qceq  =   C / D 
!
!
!       To diagnose the magnitude of each of the right hand terms of
!       (41) integrated over the time step, define the average cloud
!       condensate in the interval t to t + dtcloud qcbar as:
!
!  (44) qcbar  = qceq - [ (qc1-qc0) / ( dtcloud * D ) ]
! 
!       from which the magnitudes of the C and D terms integrated
!       over the time step are:
!
!       C   and   -D * (qcbar)
!   
!
!       Additional notes on this analytic integration:
!
!       1.   Because of finite machine precision it is required that
!            D*dt is greater than a minimum value.  This minimum
!            alue occurs where 1. - exp(-D*dt) = 0. instead of 
!            D*dt.  This value will be machine dependent. See discussion
!            at top of code for Dmin.
!

        !C_dt is set to large-scale condensation. Sink of cloud liquid 
        !is set to the sum of D1 (liquid to rain component), and D2 
        !(liquid to ice component), D_eros (erosion), and large-scale 
        !evaporation (note use of ql mean).        

        C_dt = max(dcond_ls,0.)
        D_dt = D1_dt + D2_dt + D_eros +                                &
               (max(-1.*dcond_ls,0.)/max(ql_mean,qmin)) 
                             
        !do analytic integration      
        where ( D_dt.gt.Dmin ) 
             qc0   = ql_upd
             qceq  = C_dt   /  D_dt
             qc1   = qceq - (qceq - qc0) * exp ( -1.* D_dt )
             qcbar = qceq - ((qc1 - qc0)/ D_dt)
        elsewhere
             qc0   = ql_upd
             qceq  = qc0 + C_dt   
             qc1   = qc0 + C_dt
             qcbar = qc0 + 0.5*C_dt
        end where

        !set total tendency term and update cloud
        !Note that the amount of SL calculated here is stored in tmp1.
        SL(:,:,j)  = SL(:,:,j) + qc1 - qc0
        tmp1       = qc1 - qc0        
        ql_upd     = qc1

        !compute the amount each term contributes to the change     
        Dterm  = -D_dt *      qcbar

!       Apportion SL between various processes.  This is necessary to
!       account for how much the temperature changes due to various
!       phase changes.   For example:
!
!       liquid to ice   = (D2/D)*(-Dterm)
!
!       liquid to rain  = (D1/D)*(-Dterm) 
!
!       (no phase change but needed to know how much to increment 
!        rainflux)
!
!       vapor to liquid = - { ((-dcond_ls/ql_mean)+D_eros)/D}*(-Dterm) 
!                        where dcond_ls < 0 
!                                  
!                         but
!
!                        dcond_ls  -(D_eros/D)*(-Dterm)
!                        where dcond_ls > 0
!

        !initialize tmp2 to hold (-Dterm)/D
        tmp2 = -Dterm/max(D_dt,Dmin)
        
        !do phase changes from large-scale processes and boundary
        !layer condensation/evaporation
 
        ST(:,:,j) = ST(:,:,j) + (hlv*max(dcond_ls,0.)/cp_air) -          &
             (hlv*(max(-1.*dcond_ls,0.) /max(ql_mean,qmin))*tmp2/cp_air)
   
        SQ(:,:,j) = SQ(:,:,j) -      max(dcond_ls,0.)     +            &
                  (max(-1.*dcond_ls,0.) /max(ql_mean,qmin))*tmp2
            
        !add in liquid to ice and cloud erosion to temperature tendency
        ST(:,:,j) = ST(:,:,j) + (hlf*D2_dt-hlv*D_eros)*tmp2/cp_air

        !cloud evaporation adds to water vapor
        SQ(:,:,j) = SQ(:,:,j) + D_eros*tmp2  
             
        !add conversion of liquid to rain to the rainflux
        rain_cld = rain_cld +D1_dt*tmp2*deltpg*inv_dtcloud
     
        !save liquid converted to ice into tmp3 and increment qi_mean
        tmp3    = tmp2*D2_dt
        qi_mean = qi_mean + tmp3
     
        
        if (do_liq_num) then
!******************************************************************

!       Analytic integration of qn equation
!
!       The qn equation is written in the form:
!
!  (m1) dqc/dt    =   (1 - qabar) * A * qc^   -  qc * D
!                 =   C - qc * D
!
!       where  qc^ is the large-scale qn calculated from the
!       activation parameterization.   Note that over the physics
!       time step, C and D are assumed to be constants.
!
!       Defining qc(t) = qc0 and qc(t+dtcloud) = qc1, the analytic
!       solution of the above equation is:
!
!  (m2) qc1 = qceq -  (qceq - qc0) * exp (-D*dtcloud)
! 
!       where qceq is the equilibrium cloud droplet number that is approached
!       with an time scale of 1/D,
!
!  (m3) qceq  = C / D
!
!
!       To diagnose the magnitude of each of the right hand terms of
!       (m1) integrated over the time step, define the average cloud
!       condensate in the interval t to t + dtcloud qcbar as:
!
!  (m4) qcbar  = qceq - [ (qc1-qc0) / ( dtcloud * D ) ]
! 
!       from which the magnitudes of the C and D terms integrated
!       over the time step are:
!
!       C and -D * (qcbar)
!

!Calculate C_dt
        C_dt=max(tmp5,0.)*drop1*1.e6/airdens(:,:,j)
!For replying the review, substract autoconversion
         D_dt = D1_dt + D2_dt + D_eros 

        !do analytic integration      
        where ( (D_dt.gt.Dmin) ) 
             qc0   = qn_upd
             qceq  = C_dt  / max(D_dt, Dmin)
             qc1   = qceq - (qceq - qc0) * exp ( -1.* D_dt )
             qcbar = qceq - ((qc1 - qc0)/ max(D_dt, Dmin))
        elsewhere
             qc0   = qn_upd
             qceq  = qc0 + C_dt   
             qc1   = qc0 + C_dt
             qcbar = qc0 + 0.5*C_dt
        end where

        !set total tendency term and update cloud
        !Note that the amount of SN calculated here is stored in tmp1.
        SN(:,:,j)  = SN(:,:,j) + qc1 - qc0
        qn_upd     = qc1


        !compute the amount each term contributes to the change 
!        where ( C_dt .gt. 0 )
!                Cterm  =  C_dt             
!                Dterm  =  D_dt *      qcbar 
!        elsewhere
!                Cterm  =  0.             
!                Dterm  =  D_dt *      qcbar 
!        end where

        if (do_budget_diag) then
          qndt_cond(:,:,j) = 0.
          where ( C_dt .gt. 0 )
            qndt_cond(:,:,j)  =  C_dt 
          endwhere

          qndt_evap(:,:,j) = D_dt *      qcbar !Dterm
        end if
        endif

!****************************************************************************


!
!       diagnostics for cloud liquid tendencies
!       

        if (do_budget_diag) then
             qldt_cond(:,:,j)  = max(dcond_ls,0.) *inv_dtcloud
             qldt_evap(:,:,j)  = (max(0.,-1.*dcond_ls )/max(ql_mean,   &
                                 qmin))           *tmp2*inv_dtcloud
             qldt_accr(:,:,j)  = qldt_accr (:,:,j)*tmp2*inv_dtcloud
             qldt_auto(:,:,j)  = qldt_auto (:,:,j)*tmp2*inv_dtcloud
             qldt_eros(:,:,j)  = D_eros           *tmp2*inv_dtcloud 
             qldt_berg(:,:,j)  = qldt_berg (:,:,j)*tmp2*inv_dtcloud
             qldt_rime(:,:,j)  = qldt_rime (:,:,j)*tmp2*inv_dtcloud
             qldt_freez(:,:,j) = qldt_freez(:,:,j)*tmp2*inv_dtcloud
           if(do_liq_num) then
             qndt_cond(:,:,j)  = qndt_cond(:,:,j)*inv_dtcloud 
             qndt_evap(:,:,j)  = qndt_evap(:,:,j)*inv_dtcloud 
           endif  
        end if
        


!-----                                                            -----! 
!                                                                      !
!                                END OF                                !
!                                                                      !
!                       LIQUID PHASE MICROPHYSICS                      !
!                                                                      !
!                                 AND                                  !
!                                                                      !
!                  ANALYTIC INTEGRATION OF QL EQUATION                 !
!                                                                      !
!                               SECTION                                !
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                        ICE PHASE MICROPHYSICS                        !
!                                                                      !
!                                 AND                                  !
!                                                                      !
!                  ANALYTIC INTEGRATION OF QI EQUATION                 !
!                                                                      !
!                                                                      !
!                                                                      !
!       Ice settling
!
!       Ice settling is treated as in Heymsfield Donner 1990. 
!       The mass weighted fall speed is parameterized as in equation
!       #49.
!
!       In terms of the analytic integration with respect qi of Tiedtke,
!       the flux in from the top of the grid layer is equated to the
!       source term, and the flux out of the bottom of the layer is 
!       equated to the sink term:
!
!  (47) C_dt =  snow_cld * dtcloud * grav / deltp
!
!  (48) D_dt =  airdens * grav * Vfall * dtcloud / deltp
!
!       All ice crystals are assumed to fall with the same fall speed
!       which is given as in Heymsfield and Donner (1990) as:
!
!  (49) Vfall = 3.29 * ( (airdens*qi_mean/qa_mean)**0.16)
!
!       which is the formula in Heymsfield and Donner.  Note however
!       that because this is uncertain, sensitivity runs will be made
!       with different formulations. Note that when Vfall is computed 
!       the source incremented qi, qi_mean, is used.  This gives some 
!       implicitness to the scheme.

        !compute Vfall
        Vfall = 3.29*((airdens(:,:,j)*qi_mean/max(qa_mean,qmin))**0.16)

        if (id_vfall > 0) vfalldiag(:,:,j) = Vfall(:,:)*areaice(:,:,j)

        !add to ice source the settling ice flux from above
        !also note that tmp3 contains the source
        !of liquid converted to ice from above
        C_dt = tmp3 + snow_cld*dtcloud/deltpg
        
        !Compute settling of ice. The result is multiplied by 
        !dtcloud/deltp to convert to units of D_dt.  
        !Note that if tracers are not advected then this is done
        !relative to the local vertical motion.
        if (tracer_advec) then
             tmp1 = 0.
        else
             tmp1 = omega(:,:,j)
        end if 
        
        where (qi_mean .gt. qmin .and. qa_mean .gt. qmin)
             D1_dt      = max(0.,((airdens(:,:,j)*Vfall)+(tmp1/grav))* &
                          dtcloud/deltpg )
             a_snow_cld = qa_mean     
        elsewhere
             D1_dt      = 0.
             snow_cld   = 0.
             a_snow_cld = 0.    
        end where 

        
!       Melting of in-cloud ice
!
!       Melting occurs where the temperature is greater than Tfreezing. 
!       This is an instaneous process such that no stratiform ice will
!       remain in the grid at the end of the timestep. 
!       No ice settles out of the grid box (i.e. D1 is set to zero) when
!       the amount of ice to be melted is less than that that would
!       bring the grid box to the freezing point.
!
!       The ice that melts becomes rain.  This is because if an ice
!       crystal of dimension ~100 microns and mass density of 100 kg/m2
!       melts it will become a droplet of SIZE 40 microns which is 
!       clearly a drizzle SIZE drop.  Ice crystals at temperatures near 
!       freezing are assumed to be this large, consistent with the 
!       assumption of particle SIZE temperature dependence.
!

        !compute grid mean change in cloud ice to cool the
        !grid box to 0C and store in temporary variable tmp1
        tmp1 = cp_air*(T(:,:,j)-tfreeze)/hlf

        ! If qi_mean > tmp1, then the amount of ice melted is
        ! limited to tmp1, otherwise melt all qi_mean.  The amount
        ! melted is stored in tmp2
        tmp2  = max(min(qi_mean,tmp1),0.)     
        D2_dt = max(0.,log(max(qi_mean,qmin)/max(qi_mean-tmp2,qmin)))
            
        !melting of ice creates area to a_rain_cld
        where (D2_dt .gt. Dmin) 
             a_rain_cld = qa_mean
        endwhere
                  
        !If all of the ice can melt, then don't permit any ice to fall
        !out of the grid box and set a_snow_cld to zero.
        where (qi_mean .lt. tmp1 .and. qi_mean .gt. qmin) 
             D1_dt = 0.
             snow_cld = 0.
             a_snow_cld = 0.
        end where
         
!       Analytic integration of qi equation
!
!       This repeats the procedure done for the ql equation.
!       See above notes for detail.

        !At this point C_dt already includes the source of cloud ice 
        !falling from above as well as liquid converted to ice. 
        !Therefore add in large_scale deposition.
        !
        !Compute D_dt which has contributions from D1_dt (ice settling)
        !and D2_dt (ice melting), D_eros (cloud erosion), and large-
        !scale sublimation (note use of qi mean).

        C_dt = C_dt + max(dcond_ls_ice,0.)
        D_dt =  D1_dt + D2_dt + D_eros +                               &
                (max(-1.*dcond_ls_ice,0.)/max(qi_mean,qmin))
        
        !do analytic integration      
        where ( D_dt.gt.Dmin ) 
             qc0   = qi_upd
             qceq  = C_dt / D_dt
             qc1   = qceq - (qceq - qc0) * exp ( -1.* D_dt )
             qcbar = qceq - ((qc1 - qc0)/D_dt)
        elsewhere
             qc0   = qi_upd
             qceq  = qc0 + C_dt   
             qc1   = qc0 + C_dt
             qcbar = qc0 + 0.5*C_dt
        end where

        !set total tendency term and update cloud
        !Note that the amount of SL calculated here is stored in tmp1.
        SI(:,:,j)  = SI(:,:,j) + qc1 - qc0
        tmp1       = qc1 - qc0        
        qi_upd     = qc1

        !compute the amount each term contributes to the change     
        Dterm  = -D_dt *          qcbar 
      
!       Apportion SI between various processes.  This is necessary to
!       account for how much the temperature and water vapor changes 
!       due to various phase changes.   For example:
!
!       ice settling = (D1/D)*(-Dterm)*deltp/grav/dtcloud
!                     
!       vapor to ice =
!           -{ ((-dcond_ls_ice/qi_mean)+D_eros)/ D }*(-Dterm) 
!           where dcond_ls_ice  < 0. 
!
!           but
!       
!           dcond_ls_ice -(D_eros/D)* (-Dterm)
!
!           where dcond_ls_ice > 0.
!       
!       melting of ice = (D2/D)*(-Dterm)*deltp/grav/dtcloud
!

        !initialize tmp2 to hold (-Dterm)/D
        tmp2 = -Dterm/max(D_dt,Dmin)
        
        !do phase changes from large-scale processes 
        ST(:,:,j) = ST(:,:,j) +  hls*max(dcond_ls_ice,0.)/cp_air -    &
         hls*(max(-1.*dcond_ls_ice,0.)/max(qi_mean,qmin))*tmp2/cp_air
       
        SQ(:,:,j) = SQ(:,:,j) -      max(dcond_ls_ice,0.)    +         &
             (max(-1.*dcond_ls_ice,0.)/max(qi_mean,qmin))*tmp2
     
        !cloud erosion changes temperature and vapor
        ST(:,:,j) = ST(:,:,j) - hls*D_eros* tmp2/cp_air
        SQ(:,:,j) = SQ(:,:,j) +     D_eros* tmp2

        !add settling ice flux to snow_cld 
        snow_cld = D1_dt*tmp2*deltpg*inv_dtcloud
       
        !add melting of ice to temperature tendency
        ST(:,:,j) = ST(:,:,j) - hlf*D2_dt*tmp2/cp_air

        !add melting of ice to the rainflux
        rain_cld = rain_cld + D2_dt*tmp2*deltpg*inv_dtcloud

!
!       diagnostics for cloud ice tendencies
!       
        
        if (do_budget_diag) then
             qidt_dep (:,:,j) = max(dcond_ls_ice,0.)*inv_dtcloud
             qidt_subl(:,:,j) = (max(0.,-1.*dcond_ls_ice)/max(qi_mean, &
                        qmin))*tmp2*inv_dtcloud
             qidt_melt(:,:,j) = D2_dt *tmp2*inv_dtcloud
             qidt_eros(:,:,j) = D_eros*tmp2*inv_dtcloud       
        end if


        
!-----                                                            -----! 
!                                                                      !
!                                END OF                                !
!                                                                      !
!                        ICE PHASE MICROPHYSICS                        !
!                                                                      !
!                                 AND                                  !
!                                                                      !
!                  ANALYTIC INTEGRATION OF QI EQUATION                 !
!                                                                      !
!                                 SECTION                              !
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!
!
!                       RAIN EVAPORATION
!                           
!
!       Rain evaporation is derived by integration of the growth 
!       equation of a droplet over the assumed Marshall-Palmer 
!       distribution of rain drops (equation #22).  This leads to the 
!       following formula:
!
!  (50) dqv/dt_local =  56788.636 * {rain_rate/dens_h2o}^(11/18) *(1-U)/
!
!                      ( SQRT(airdens)* A_plus_B)
!
!       Numerically this equation integrated by use of time-centered
!       values for qs and qv.   This leads to the solution:
!
!  (51) qv_clr(t+1)-qv_clr(t) = K3 *[qs(t)-qv_clr(t)]/
!                               {1.+0.5*K3*(1+gamma)}
!       where 
!
!       K3= 56788.636 * dtcloud * {rain_rate_local/dens_h2o}^(11/18) /
!           ( SQRT(airdens)* A_plus_B * qs)
!
!       and gamma is given by (3). Note that in (51), it is made 
!       explicit that it is the vapor concentration in the unsaturated 
!       part of the grid box that is used in the rain evaporation 
!       formula.
!
!       Now there are several limiters to this formula. First,
!       you cannot evaporate more than is available in a time step.
!       The amount available for evaporation locally is 
!       (rain_clr/a_rain_clr)*(grav*dtcloud/deltp).   Second, to
!       avoid supersaturating the box or exceeding the critical
!       relative humidity above which rain does not evaporate, 
!       the amount of evaporation is limited.
!
!       Finally rain evaporation occurs only if the relative humidity
!       in the unsaturated portion of the grid box, U_clr, is less
!       then a threshold, U_evap.   U_evap, will not necessarily be
!       one.   For example, stratiform precipitation in convective
!       regions rarely saturates subcloud air because of the presence
!       of downdrafts.  If the convection scheme does not have down-
!       drafts then it doesn't make sense to allow the sub-cloud layer
!       to saturate. U_clr may be solved from (8) as:
!
!  (52) U_clr = ( U - qa ) / (1. - qa)
!
!       Some variables are temporarily stored in tmp1.
!
!       Note that for pdf clouds the relative humidity in the clear part
!       of the grid box can be calculated exactly from the beta distr-
!       ibution. 

        !compute U_clr
        if (.not. do_pdf_clouds) then 
             U_clr =  (U-qa_mean)/max((1.-qa_mean),qmin)
        else
             U_clr = qvg/qs(:,:,j)
        end if
        
        !keep U_clr > 0. and U_clr < 1.
        U_clr = min(max(U_clr,0.),1.)
        
        !compute K3
        tmp1 = 56788.636 * dtcloud * ((rain_clr/max(a_rain_clr,qmin)/  &
             dens_h2o)**(11./18.))/SQRT(airdens(:,:,j))/A_plus_B(:,:,j)&
             /qs(:,:,j)

        !compute local change in vapor mixing ratio due to 
        !rain evaporation
        tmp1 = tmp1*qs(:,:,j)*(1.-U_clr)/(1.+0.5*tmp1*(1.+gamma(:,:,j)))

        !limit change in qv to the amount that would raise the relative
        !humidity to U_evap in the clear portion of the grid box
        tmp1 = min(tmp1,((1.-qa_mean)/max(a_rain_clr,qmin))*qs(:,:,j)* &
               max(0.,U_evap-U_clr)/(1.+(U_evap*(1.-qa_mean)+qa_mean)* &
               gamma(:,:,j)) )
        
        !do limiter by amount available
        tmp1= tmp1*a_rain_clr*deltpg*inv_dtcloud
        tmp2= max(min(rain_clr,tmp1),0.)
    
        SQ(:,:,j) = SQ(:,:,j) +     tmp2*dtcloud/deltpg
        ST(:,:,j) = ST(:,:,j) - hlv*tmp2*dtcloud/deltpg/cp_air
        
        !if all of the rain evaporates set things to zero.    
        where (tmp1.gt.rain_clr.and.a_rain_clr.gt.qmin)         
             rain_clr = 0.
             a_rain_clr = 0.
        elsewhere
             rain_clr = rain_clr - tmp2     
        endwhere
        
        if (do_budget_diag) rain_evap(:,:,j) = tmp2/deltpg


!-----------------------------------------------------------------------
!
!
!
!                              SNOW SUBLIMATION
!                           
!
!       Sublimation of cloud ice
!
!       [The following follows Rotstayn (1997)]
!       Given the assumptions of the Marshall-Palmer distribution of
!       ice crystals (18), the crystal growth equation as a function
!       of the humidity of the air and the diffusivity of water vapor
!       and thermal conductivity of air is integrated over all crystal
!       sizes.   This yields:
!
!  (53) dqi/dt_local = - a_snow_clr* K3 * (qs - qv_clr)
!
!       where the leading factor of a_snow_clr is the portion of the
!       grid box undergoing sublimation. K3 is given by
!
!  (54) K3 = (4/(pi*rho_air*qs*rho_ice*A_plus_B))*
!            ((snow_clr/a_snow_clr/3.29)**1.16 ) *
!           [ 0.65*lamda_f^2 + 
!             198.92227 * (airdens)^0.5 * 
!             ((snow_clr/a_snow_clr)**(1/14.5)) * lamda_f^(3/2) ]
!
!       Note that equation (53) is identical to equation (30) of 
!       Rotstayn.
!
!       Numerically this is integrated as in rain evaporation.


        !compute K3
        tmp1 = dtcloud * (4./3.14159/rho_ice/airdens(:,:,j)/           &
               A_plus_B(:,:,j)/qs(:,:,j))*((snow_clr/max(a_snow_clr,   &
               qmin)/3.29)**(1./1.16))*(0.65*lamda_f*lamda_f +         &
               198.92227*lamda_f*SQRT(airdens(:,:,j)*lamda_f)*         &
               ( (snow_clr/max(a_snow_clr,qmin))**(1./14.5) )  )

        !compute local change in vapor mixing ratio due to 
        !snow sublimation
        tmp1 = tmp1*qs(:,:,j)*(1.-U_clr)/(1.+0.5*tmp1*(1.+gamma(:,:,j)))

        !limit change in qv to the amount that would raise the relative
        !humidity to U_evap in the clear portion of the grid box
        tmp1 = min(tmp1,((1.-qa_mean)/max(a_snow_clr,qmin))*qs(:,:,j)* &
               max(0.,U_evap-U_clr)/(1.+(U_evap*(1.-qa_mean)+qa_mean)* &
               gamma(:,:,j)) )
        
        !do limiter by amount available
        tmp1= tmp1*a_snow_clr*deltpg*inv_dtcloud
        tmp2= max(min(snow_clr,tmp1),0.)
    
        SQ(:,:,j) = SQ(:,:,j) +     tmp2*dtcloud/deltpg
        ST(:,:,j) = ST(:,:,j) - hls*tmp2*dtcloud/deltpg/cp_air
        
        !if all of the snow sublimates set things to zero.    
        where (tmp1.gt.snow_clr.and.a_snow_clr.gt.qmin)         
             snow_clr = 0.
             a_snow_clr = 0.
        elsewhere
             snow_clr = snow_clr - tmp2     
        endwhere
         
        if (do_budget_diag) snow_subl(:,:,j) = tmp2/deltpg        


!-----------------------------------------------------------------------
!
!       Adjustment.  Due to numerical errors in detrainment or advection
!       sometimes the current state of the grid box may be super-
!       saturated. Under the assumption that the temperature is constant
!       in the grid box and that q <= qs, the excess vapor is condensed. 
!       
!       What happens to the condensed water vapor is determined by the
!       namelist parameter super_choice.
!
!       If super_choice = .false. (default), then the condensed water is
!       is added to the precipitation fluxes.  If super_choice = .true.,
!       then the condensed water is added to the cloud condensate field.
!       Note that in this case the cloud fraction is raised to one.
!
!       The phase partitioning depends on super_choice; if super_choice
!       is false then at T < -20C, snow is produced.  If super_choice
!       is true, then at T < -40C, ice is produced.  The latter choice 
!       is consistent with that done in the large-scale condensation
!       section above.        
!
!       If pdf clouds are operating then this section is bypassed - 
!       as statistical clouds should remove supersaturation according
!       to the beta distribution used.
             
        if (.not.do_pdf_clouds) then
               
        !estimate current qs
        tmp2 = qs(:,:,j)+dqsdT(:,:,j)*ST(:,:,j)

        !compute excess over saturation
        tmp1 = max(0.,qv(:,:,j)+SQ(:,:,j)-tmp2)/(1.+gamma(:,:,j))
        
        !change vapor content
        SQ(:,:,j)=SQ(:,:,j)-tmp1

        if (super_choice) then
        
             ! Put supersaturation into cloud

             !cloud fraction source diagnostic
             if (do_budget_diag) then
                  where (tmp1 .gt. 0.)
                       qadt_super(:,:,j)  = (1.-qa_upd) * inv_dtcloud
                  endwhere
             end if

             !add in excess to cloud condensate, change cloud area and 
             !increment temperature
             where (T(:,:,j) .le. tfreeze-40. .and. tmp1 .gt. 0.)
                  qi_upd   = qi_upd + tmp1
                  SI(:,:,j) = SI(:,:,j) + tmp1
                  SA(:,:,j) = SA(:,:,j) + (1.-qa_upd)
                  qa_upd   = 1.
                  ST(:,:,j)  = ST(:,:,j) + hls*tmp1/cp_air
             end where
             where (T(:,:,j) .gt. tfreeze-40. .and. tmp1 .gt. 0.)        
                   ql_upd   = ql_upd + tmp1
                   SL(:,:,j) = SL(:,:,j) + tmp1
                   SA(:,:,j) = SA(:,:,j) + (1.-qa_upd)
                   qa_upd   = 1.
                   ST(:,:,j)  = ST(:,:,j) + hlv*tmp1/cp_air              
             end where
             if (do_liq_num) then 
               where (T(:,:,j) .gt. tfreeze-40. .and. tmp1 .gt. 0.)        
                   qn_upd = qn_upd + drop1*1.e6/airdens(:,:,j)*(1.-qa_upd)
                   SN(:,:,j) = SN(:,:,j) + drop1*1.e6/airdens(:,:,j)*(1.-qa_upd)
               end where
             endif
                
             if (do_budget_diag) then             
               where (T(:,:,j) .le. tfreeze-40.)
                    ice_adj(:,:,j) = tmp1*inv_dtcloud
               elsewhere
                    liq_adj(:,:,j) = tmp1*inv_dtcloud
               endwhere
               where (T(:,:,j) .gt. tfreeze-40. .and. tmp1 .gt. 0.)
                    qndt_cond(:,:,j) = qndt_cond(:,:,j) + drop1*1.e6 / &
                                       airdens(:,:,j)*(1.-qa_upd)*inv_dtcloud
               endwhere                       
             end if

        else

             !Put supersaturation into precip

             !add in excess to precipitation fluxes, change their area 
             !and increment temperature
             where (T(:,:,j) .le. tfreeze-20. .and. tmp1 .gt. 0.)
                  snow_cld = snow_cld + qa_mean *tmp1*deltpg*inv_dtcloud
                  snow_clr = snow_clr + (1.-qa_mean)*tmp1*deltpg*      &
                                                             inv_dtcloud
                  a_snow_cld = qa_mean
                  a_snow_clr = 1.-qa_mean
                  ST(:,:,j)  = ST(:,:,j) + hls*tmp1/cp_air
             end where
             where (T(:,:,j) .gt. tfreeze-20. .and. tmp1 .gt. 0.)        
                  rain_cld = rain_cld + qa_mean *tmp1*deltpg*inv_dtcloud
                  rain_clr = rain_clr + (1.-qa_mean)*tmp1*deltpg*      &
                                                             inv_dtcloud
                  a_rain_cld = qa_mean
                  a_rain_clr = 1.-qa_mean
                  ST(:,:,j)  = ST(:,:,j) + hlv*tmp1/cp_air
             end where
        
             if (do_budget_diag) then     
                  where (T(:,:,j) .le. tfreeze-20.)
                       ice_adj(:,:,j) = tmp1*inv_dtcloud
                  elsewhere
                       liq_adj(:,:,j) = tmp1*inv_dtcloud
                  endwhere
             end if

        end if !super choice
        
        end if !for do_pdf_clouds
                      
!-----------------------------------------------------------------------
!
!       Cloud Destruction occurs where both ql and qi are .le. qmin, 
!       or if qa is .le. qmin. In this case, ql, qi, and qa are set to 
!       zero conserving moisture and energy.

        where ((ql_upd .le. qmin .and. qi_upd .le. qmin)               &
               .or. (qa_upd .le. qmin))
             SL(:,:,j) = SL(:,:,j) - ql_upd
             SI(:,:,j) = SI(:,:,j) - qi_upd
             SQ(:,:,j) = SQ(:,:,j) + ql_upd + qi_upd
             ST(:,:,j) = ST(:,:,j) - (hlv*ql_upd + hls*qi_upd)/cp_air
             SA(:,:,j) = SA(:,:,j) - qa_upd
        end where
        if(do_liq_num) then 
          where ((ql_upd .le. qmin .and. qi_upd .le. qmin)               &
               .or. (qa_upd .le. qmin))
             SN(:,:,j) = SN(:,:,j) - qn_upd
          end where
        endif  

        if (do_budget_diag) then

             where ((ql_upd .le. qmin .and. qi_upd .le. qmin)          &
                    .or. (qa_upd .le. qmin))
                  qldt_destr(:,:,j) = ql_upd*inv_dtcloud
                  qidt_destr(:,:,j) = qi_upd*inv_dtcloud
                  qadt_destr(:,:,j) = qa_upd*inv_dtcloud
             endwhere

        end if

!-----------------------------------------------------------------------
!
!       Save qa_mean of current level into qa_mean_lst.   This is used
!       in transferring rain and snow fluxes between levels.

        qa_mean_lst = qa_mean
        
!-----------------------------------------------------------------------
!
!       add the ice falling out from cloud to qidt_fall
        
        if (do_budget_diag) qidt_fall(:,:,j) = qidt_fall(:,:,j) +    & 
                                         (snow_cld/deltpg)
        
!-----------------------------------------------------------------------
!
!       save profiles of rain and snow 
!
        rain3d(:,:,j+1) = rain_clr(:,:) + rain_cld(:,:)
        snow3d(:,:,j+1) = snow_clr(:,:) + snow_cld(:,:)

!-----------------------------------------------------------------------
!
!       Save rain and snow diagnostics

        if (do_budget_diag) then
             rain_clr_diag(:,:,j+1)     = rain_clr
             rain_cld_diag(:,:,j+1)     = rain_cld
             a_rain_clr_diag(:,:,j+1)   = a_rain_clr
             a_rain_cld_diag(:,:,j+1)   = a_rain_cld
             snow_clr_diag(:,:,j+1)     = snow_clr
             snow_cld_diag(:,:,j+1)     = snow_cld
             a_snow_clr_diag(:,:,j+1)   = a_snow_clr
             a_snow_cld_diag(:,:,j+1)   = a_snow_cld
             a_precip_clr_diag(:,:,j+1) = max(a_rain_clr,a_snow_clr)
             a_precip_cld_diag(:,:,j+1) = max(a_rain_cld,a_snow_cld)
        end if

!-----------------------------------------------------------------------
!
!
!       Put rain and ice fluxes into surfrain and surfsnow if the
!       grid point is at the bottom of a column.   If MASK is not
!       present then this code is executed only if j .eq. kdim.
!       IF MASK is present some grid points may be beneath ground. 
!       If a given grid point is at the bottom of the column then
!       the surface values of rain and snow must be created.
!       Also if the MASK is present then the code forces all tenden-
!       cies below ground to be zero. Note that MASK = 1. equals above
!       ground point, MASK = 0. equals below ground point.

        if (present(MASK)) then

             !zero out all tendencies below ground
             ST(:,:,j)=MASK(:,:,j)*ST(:,:,j)
             SQ(:,:,j)=MASK(:,:,j)*SQ(:,:,j)
             SL(:,:,j)=MASK(:,:,j)*SL(:,:,j)
             SI(:,:,j)=MASK(:,:,j)*SI(:,:,j)
             SA(:,:,j)=MASK(:,:,j)*SA(:,:,j)
             if (do_liq_num) then
               SN(:,:,j)=MASK(:,:,j)*SN(:,:,j)
             endif
 
             if (j .lt. kdim) then
                  
                  !bottom of true points in columns which contain some
                  !dummy points
                  where(MASK(:,:,j) .eq. 1. .and. MASK(:,:,j+1) .eq. 0.)
                       surfrain = dtcloud*(rain_clr+rain_cld)
                       surfsnow = dtcloud*(snow_clr+snow_cld)
                       rain_clr = 0.
                       rain_cld = 0.
                       snow_clr = 0.
                       snow_cld = 0.
                       a_rain_clr = 0.
                       a_rain_cld = 0.
                       a_snow_clr = 0.
                       a_snow_cld = 0.
                  end where

             else

                  !bottom of column for those columns which contain no
                  !dummy points
                  where(MASK(:,:,j) .eq. 1.)
                       surfrain = dtcloud*(rain_clr+rain_cld)
                       surfsnow = dtcloud*(snow_clr+snow_cld)
                       rain_clr = 0.
                       rain_cld = 0.
                       snow_clr = 0.
                       snow_cld = 0.
                       a_rain_clr = 0.
                       a_rain_cld = 0.
                       a_snow_clr = 0.
                       a_snow_cld = 0.                  
                  end where

             end if

        else

             !do code if we are at bottom of column
             if (j .eq. kdim) then
                  surfrain = dtcloud*(rain_clr+rain_cld)
                  surfsnow = dtcloud*(snow_clr+snow_cld)
             end if

        end if 
                  

!-----------------------------------------------------------------------
!
!       END LOOP OVER VERTICAL LEVELS
!

        enddo




!-----------------------------------------------------------------------
!
!       INSTANTANEOUS OUTPUT DIAGNOSTICS
!
     
        if (num_strat_pts > 0) then
         do nn=1,num_strat_pts
          if (strat_pts(1,nn) >= is .and. strat_pts(1,nn) <= ie .and.  &
             strat_pts(2,nn) >= js .and. strat_pts(2,nn) <= je) then
                ipt=strat_pts(1,nn); jpt=strat_pts(2,nn)
                i=ipt-is+1; j=jpt-js+1
                unit = open_ieee32_file ('strat.data', action='append')
                write (unit) ipt,jpt,     ql(i,j,:)+SL(i,j,:)
                write (unit) ipt,jpt,     qi(i,j,:)+SI(i,j,:)
                write (unit) ipt,jpt,     qa(i,j,:)+SA(i,j,:)
                write (unit) ipt,jpt,      T(i,j,:)+ST(i,j,:) 
                write (unit) ipt,jpt,     qv(i,j,:)+SQ(i,j,:)
                write (unit) ipt,jpt,     pfull(i,j,:)
                call close_file(unit)
          endif
         enddo
        endif


!-----------------------------------------------------------------------
!
!       DIAGNOSTICS
!

        if (id_droplets > 0) then
             used = send_data ( id_droplets, N3D, Time, is, js, 1,&
                                  rmask=mask )
        end if

        if (id_aall > 0) then
             used = send_data ( id_aall, areaall, Time, is, js, 1,     &
                                  rmask=mask )
             deallocate(areaall)
        end if
        if (id_aliq > 0) then
             used = send_data ( id_aliq, arealiq, Time, is, js, 1,     &
                                  rmask=mask )
             deallocate(arealiq)
        end if
        if (id_aice > 0) then
             used = send_data ( id_aice, areaice, Time, is, js, 1,     &
                                  rmask=mask )
             deallocate(areaice)
        end if
        if (id_rvolume > 0) then
             used = send_data ( id_rvolume, rvolume, Time, is, js, 1,  &
                                  rmask=mask )
             deallocate(rvolume)
        end if
        if (id_autocv > 0) then
             used = send_data ( id_autocv, areaautocv, Time, is, js, 1,&
                                  rmask=mask )
             deallocate(areaautocv)
        end if
        if (id_vfall > 0) then
             used = send_data ( id_vfall, vfalldiag, Time, is, js, 1,  &
                                  rmask=mask )  
             deallocate(vfalldiag)
        end if

        if (do_budget_diag) then

        !------- set up half level mask --------
        mask3(:,:,1:(kdim+1)) = 1.
        if (present(mask)) then
             where (mask(:,:,1:kdim) <= 0.5)
                  mask3(:,:,2:(kdim+1)) = 0.
             end where
        endif
        

        !cloud liquid, droplet number and rain diagnostics
        if ( id_qldt_cond > 0 ) then
             used = send_data ( id_qldt_cond, qldt_cond, &
                               Time, is, js, 1, rmask=mask )
        endif

        if ( id_qldt_evap > 0 ) then
             used = send_data ( id_qldt_evap, qldt_evap, &
                               Time, is, js, 1, rmask=mask )
        endif

        if ( id_qldt_eros > 0 ) then
             used = send_data ( id_qldt_eros, qldt_eros, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qldt_accr > 0 ) then
             used = send_data ( id_qldt_accr, qldt_accr, &
                               Time, is, js, 1, rmask=mask )
        endif

        if ( id_qldt_auto > 0 ) then
             used = send_data ( id_qldt_auto, qldt_auto, &
                               Time, is, js, 1, rmask=mask )
        endif

        if ( id_liq_adj > 0 ) then
             used = send_data ( id_liq_adj, liq_adj, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qldt_fill > 0 ) then
             used = send_data ( id_qldt_fill, qldt_fill, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qldt_berg > 0 ) then
             used = send_data ( id_qldt_berg, qldt_berg, &
                                Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qldt_freez > 0 ) then
             used = send_data ( id_qldt_freez, qldt_freez, &
                                Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qldt_rime > 0 ) then
             used = send_data ( id_qldt_rime, qldt_rime, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qldt_destr > 0 ) then
             used = send_data ( id_qldt_destr, qldt_destr, &
                               Time, is, js, 1, rmask=mask )
        endif

        if ( id_qndt_cond > 0 ) then
             used = send_data ( id_qndt_cond, qndt_cond, &
                               Time, is, js, 1, rmask=mask )
        endif

        if ( id_qndt_evap > 0 ) then
             used = send_data ( id_qndt_evap, qndt_evap, &
                               Time, is, js, 1, rmask=mask )
        endif

        if ( id_rain_evap > 0 ) then
             used = send_data ( id_rain_evap, rain_evap, &
                               Time, is, js, 1, rmask=mask )
        endif

        if ( id_rain_clr > 0 ) then
             used = send_data ( id_rain_clr, rain_clr_diag, &
                               Time, is, js, 1, rmask=mask3 )
        end if
         
        if ( id_a_rain_clr > 0 ) then
             used = send_data ( id_a_rain_clr, a_rain_clr_diag, &
                               Time, is, js, 1, rmask=mask3 )
        end if
         
        if ( id_rain_cld > 0 ) then
             used = send_data ( id_rain_cld, rain_cld_diag, &
                               Time, is, js, 1, rmask=mask3 )
        end if
         
        if ( id_a_rain_cld > 0 ) then
             used = send_data ( id_a_rain_cld, a_rain_cld_diag, &
                               Time, is, js, 1,rmask=mask3 )
        end if
       
        if ( id_a_precip_clr > 0 ) then
             used = send_data ( id_a_precip_clr, a_precip_clr_diag, &
                               Time, is, js, 1, rmask=mask3 )
        end if
         
        if ( id_a_precip_cld > 0 ) then
             used = send_data ( id_a_precip_cld, a_precip_cld_diag, &
                               Time, is, js, 1,rmask=mask3 )
        end if
       
       
        !ice and snow diagnostics        
        if ( id_qidt_dep > 0 ) then
             used = send_data ( id_qidt_dep, qidt_dep, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qidt_subl> 0 ) then
             used = send_data ( id_qidt_subl, qidt_subl, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qidt_eros > 0 ) then
             used = send_data ( id_qidt_eros, qidt_eros, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qidt_fall > 0 ) then
             used = send_data ( id_qidt_fall, qidt_fall, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qidt_melt > 0 ) then
             used = send_data ( id_qidt_melt, qidt_melt, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_ice_adj > 0 ) then
             used = send_data ( id_ice_adj, ice_adj, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_qidt_destr > 0 ) then
             used = send_data ( id_qidt_destr, qidt_destr, &
                               Time, is, js, 1, rmask=mask )
        endif

        if ( id_qidt_fill > 0 ) then
             used = send_data ( id_qidt_fill, qidt_fill, &
                               Time, is, js, 1, rmask=mask )
        endif
      
        if ( id_snow_clr > 0 ) then
             used = send_data ( id_snow_clr, snow_clr_diag, &
                               Time, is, js, 1, rmask=mask3 )
        end if
         
        if ( id_a_snow_clr > 0 ) then
             used = send_data ( id_a_snow_clr, a_snow_clr_diag, &
                               Time, is, js, 1, rmask=mask3 )
        end if
         
        if ( id_snow_cld > 0 ) then
             used = send_data ( id_snow_cld, snow_cld_diag, &
                               Time, is, js, 1, rmask=mask3 )
        end if
         
        if ( id_a_snow_cld > 0 ) then
             used = send_data ( id_a_snow_cld, a_snow_cld_diag, &
                               Time, is, js, 1, rmask=mask3 )
        end if
         
        if ( id_snow_subl > 0 ) then
             used = send_data ( id_snow_subl, snow_subl, &
                               Time, is, js, 1, rmask=mask )
        endif
   
        if ( id_snow_melt > 0 ) then
             used = send_data ( id_snow_melt, snow_melt, &
                               Time, is, js, 1, rmask=mask )
        endif


        !cloud fraction diagnostics        
        if ( id_qadt_lsform > 0 ) then
             used = send_data ( id_qadt_lsform, qadt_lsform, &
                               Time, is, js, 1, rmask=mask )
        endif
        
        if ( id_qadt_lsdiss > 0 ) then
             used = send_data ( id_qadt_lsdiss, qadt_lsdiss, &
                               Time, is, js, 1, rmask=mask )
        endif
        
        if ( id_qadt_rhred > 0 ) then
             used = send_data ( id_qadt_rhred, qadt_rhred, &
                               Time, is, js, 1, rmask=mask )
        endif
        
        if ( id_qadt_eros > 0 ) then
             used = send_data ( id_qadt_eros, qadt_eros, &
                               Time, is, js, 1, rmask=mask )
        endif
        
        if ( id_qadt_fill > 0 ) then
             used = send_data ( id_qadt_fill, qadt_fill, &
                               Time, is, js, 1, rmask=mask )
        endif
        
        if ( id_qadt_super > 0 ) then
             used = send_data ( id_qadt_super, qadt_super, &
                               Time, is, js, 1, rmask=mask )
        endif
        
        if ( id_qadt_destr > 0 ) then
             used = send_data ( id_qadt_destr, qadt_destr, &
                               Time, is, js, 1, rmask=mask )
        endif
        
           
        !-------write out column integrated diagnostics------!

        !liquid and rain diagnostics
        if ( id_ql_cond_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qldt_cond    (:,:,j) = qldt_cond   (:,:,j)*deltpg        
             enddo
             do j = kdim-1, 1, -1
                  qldt_cond    (:,:,kdim) = qldt_cond    (:,:,kdim) &
                                          + qldt_cond    (:,:,j)
             enddo
             used = send_data ( id_ql_cond_col, qldt_cond(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_ql_evap_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qldt_evap    (:,:,j) = qldt_evap   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qldt_evap    (:,:,kdim) = qldt_evap    (:,:,kdim) &
                                          + qldt_evap    (:,:,j)             
             enddo
             used = send_data ( id_ql_evap_col, qldt_evap(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_ql_eros_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qldt_eros    (:,:,j) = qldt_eros   (:,:,j)*deltpg             
             enddo
             do j = kdim-1, 1, -1
                  qldt_eros    (:,:,kdim) = qldt_eros    (:,:,kdim) &
                                          + qldt_eros    (:,:,j)
             enddo
             used = send_data ( id_ql_eros_col, qldt_eros(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_ql_accr_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qldt_accr    (:,:,j) = qldt_accr   (:,:,j)*deltpg             
             enddo
             do j = kdim-1, 1, -1
                  qldt_accr    (:,:,kdim) = qldt_accr    (:,:,kdim) &
                                          + qldt_accr    (:,:,j)
             enddo
             used = send_data ( id_ql_accr_col, qldt_accr(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_ql_auto_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qldt_auto    (:,:,j) = qldt_auto   (:,:,j)*deltpg             
             enddo
             do j = kdim-1, 1, -1
                  qldt_auto    (:,:,kdim) = qldt_auto    (:,:,kdim) &
                                          + qldt_auto    (:,:,j)
             enddo
             used = send_data ( id_ql_auto_col, qldt_auto(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_ql_berg_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qldt_berg    (:,:,j) = qldt_berg   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qldt_berg    (:,:,kdim) = qldt_berg    (:,:,kdim) &
                                          + qldt_berg    (:,:,j)
             enddo
             used = send_data ( id_ql_berg_col, qldt_berg(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_ql_freez_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qldt_freez   (:,:,j) = qldt_freez  (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qldt_freez   (:,:,kdim) = qldt_freez   (:,:,kdim) &
                                          + qldt_freez   (:,:,j)
             enddo
             used = send_data ( id_ql_freez_col, qldt_freez(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_ql_destr_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qldt_destr   (:,:,j) = qldt_destr  (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qldt_destr   (:,:,kdim) = qldt_destr   (:,:,kdim) &
                                          + qldt_destr   (:,:,j)
             enddo
             used = send_data ( id_ql_destr_col, qldt_destr(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_ql_rime_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qldt_rime    (:,:,j) = qldt_rime   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qldt_rime    (:,:,kdim) = qldt_rime    (:,:,kdim) &
                                          + qldt_rime    (:,:,j)
             enddo
             used = send_data ( id_ql_rime_col, qldt_rime(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_ql_fill_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qldt_fill    (:,:,j) = qldt_fill   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qldt_fill    (:,:,kdim) = qldt_fill    (:,:,kdim) &
                                          + qldt_fill    (:,:,j)
             enddo
             used = send_data ( id_ql_fill_col, qldt_fill(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_liq_adj_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  liq_adj      (:,:,j) = liq_adj     (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  liq_adj      (:,:,kdim) = liq_adj      (:,:,kdim) &
                                          + liq_adj      (:,:,j)
             enddo
             used = send_data ( id_liq_adj_col, liq_adj(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_rain_evap_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  rain_evap    (:,:,j) = rain_evap   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  rain_evap    (:,:,kdim) = rain_evap    (:,:,kdim) &
                                          + rain_evap    (:,:,j)
             enddo
             used = send_data ( id_rain_evap_col, rain_evap(:,:,kdim), &
                                  Time, is, js )
        endif
   
   
        !ice and snow diagnostics   
        if ( id_qi_fall_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qidt_fall    (:,:,j) = qidt_fall   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qidt_fall    (:,:,kdim) = qidt_fall    (:,:,kdim) &
                                          + qidt_fall    (:,:,j)
             enddo
             used = send_data ( id_qi_fall_col, qidt_fall(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_qi_fill_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qidt_fill    (:,:,j) = qidt_fill   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qidt_fill    (:,:,kdim) = qidt_fill    (:,:,kdim) &
                                          + qidt_fill    (:,:,j)
             enddo
             used = send_data ( id_qi_fill_col, qidt_fill(:,:,kdim), &
                                          Time, is, js )
        endif
   
        if ( id_qi_dep_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qidt_dep     (:,:,j) = qidt_dep    (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qidt_dep     (:,:,kdim) = qidt_dep     (:,:,kdim)    &
                                          + qidt_dep     (:,:,j)
             enddo
             used = send_data ( id_qi_dep_col, qidt_dep(:,:,kdim),     &
                                  Time, is, js )
        endif
   
        if ( id_qi_subl_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qidt_subl    (:,:,j) = qidt_subl   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qidt_subl    (:,:,kdim) = qidt_subl    (:,:,kdim) &
                                          + qidt_subl    (:,:,j)
             enddo
             used = send_data ( id_qi_subl_col, qidt_subl(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_qi_eros_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qidt_eros    (:,:,j) = qidt_eros   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qidt_eros    (:,:,kdim) = qidt_eros    (:,:,kdim) &
                                          + qidt_eros    (:,:,j)
             enddo
             used = send_data ( id_qi_eros_col, qidt_eros(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_qi_destr_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qidt_destr   (:,:,j) = qidt_destr  (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qidt_destr   (:,:,kdim) = qidt_destr   (:,:,kdim) &
                                          + qidt_destr   (:,:,j)
             enddo
             used = send_data ( id_qi_destr_col, qidt_destr(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_qi_melt_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qidt_melt    (:,:,j) = qidt_melt   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qidt_melt    (:,:,kdim) = qidt_melt    (:,:,kdim) &
                                          + qidt_melt    (:,:,j)
             enddo
             used = send_data ( id_qi_melt_col, qidt_melt(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_ice_adj_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  ice_adj     (:,:,j)  = ice_adj     (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  ice_adj      (:,:,kdim) = ice_adj      (:,:,kdim) &
                                          + ice_adj      (:,:,j)
             enddo
             used = send_data ( id_ice_adj_col,  ice_adj(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_snow_melt_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  snow_melt    (:,:,j) = snow_melt   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  snow_melt    (:,:,kdim) = snow_melt    (:,:,kdim) &
                                          + snow_melt    (:,:,j)
             enddo
             used = send_data ( id_snow_melt_col, snow_melt(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_snow_subl_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  snow_subl    (:,:,j) = snow_subl   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  snow_subl    (:,:,kdim) = snow_subl    (:,:,kdim) &
                                          + snow_subl    (:,:,j)
             enddo
             used = send_data ( id_snow_subl_col, snow_subl(:,:,kdim), &
                                  Time, is, js )
        endif


        !cloud fraction and volume diagnostics
        if ( id_qa_lsform_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qadt_lsform  (:,:,j) = qadt_lsform (:,:,j)*deltpg             
             enddo
             do j = kdim-1, 1, -1
                  qadt_lsform  (:,:,kdim) = qadt_lsform  (:,:,kdim) &
                                          + qadt_lsform  (:,:,j)
             enddo
             used = send_data (id_qa_lsform_col, qadt_lsform(:,:,kdim),&
                                  Time, is, js )
        endif

        if ( id_qa_lsdiss_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qadt_lsdiss  (:,:,j) = qadt_lsdiss (:,:,j)*deltpg             
             enddo
             do j = kdim-1, 1, -1
                  qadt_lsdiss  (:,:,kdim) = qadt_lsdiss  (:,:,kdim) &
                                          + qadt_lsdiss  (:,:,j)
             enddo
             used = send_data (id_qa_lsdiss_col, qadt_lsdiss(:,:,kdim),&
                                  Time, is, js )
        endif

        if ( id_qa_rhred_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qadt_rhred   (:,:,j) = qadt_rhred  (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qadt_rhred   (:,:,kdim) = qadt_rhred   (:,:,kdim) &
                                          + qadt_rhred   (:,:,j)
             enddo
             used = send_data ( id_qa_rhred_col, qadt_rhred(:,:,kdim), &
                                  Time, is, js )
        endif

        if ( id_qa_eros_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qadt_eros    (:,:,j) = qadt_eros   (:,:,j)*deltpg             
             enddo
             do j = kdim-1, 1, -1
                  qadt_eros    (:,:,kdim) = qadt_eros    (:,:,kdim) &
                                          + qadt_eros    (:,:,j)
             enddo
             used = send_data ( id_qa_eros_col, qadt_eros(:,:,kdim), &
                                  Time, is, js )
        endif
   
        if ( id_qa_fill_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qadt_fill    (:,:,j) = qadt_fill   (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qadt_fill    (:,:,kdim) = qadt_fill    (:,:,kdim) &
                                          + qadt_fill    (:,:,j)
             enddo
             used = send_data ( id_qa_fill_col, qadt_fill(:,:,kdim), &
                                  Time, is, js )
        endif
           
        if ( id_qa_super_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qadt_super   (:,:,j) = qadt_super  (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qadt_super   (:,:,kdim) = qadt_super   (:,:,kdim) &
                                          + qadt_super   (:,:,j)
             enddo
             used = send_data ( id_qa_super_col, qadt_super(:,:,kdim), &
                                  Time, is, js )
        endif
                 
        if ( id_qa_destr_col > 0 ) then
             do j = 1, kdim
                  deltpg = (phalf(:,:,j+1)-phalf(:,:,j))/grav
                  if (present(MASK)) deltpg=deltpg*MASK(:,:,j)
                  qadt_destr   (:,:,j) = qadt_destr  (:,:,j)*deltpg
             enddo
             do j = kdim-1, 1, -1
                  qadt_destr   (:,:,kdim) = qadt_destr   (:,:,kdim) &
                                          + qadt_destr   (:,:,j)
             enddo
             used = send_data ( id_qa_destr_col, qadt_destr(:,:,kdim), &
                                  Time, is, js )
        endif


        !other stuff

        
   
        !---------------------------------
        !deallocate space for diagnostics

        if (do_liq_num) then
          deallocate (qndt_cond)
          deallocate (qndt_evap)
        endif  
        deallocate (qldt_cond)
        deallocate (qldt_evap)
        deallocate (qldt_eros)
        deallocate (qldt_fill)
        deallocate (qldt_accr)
        deallocate (qldt_auto)
        deallocate (qldt_freez)
        deallocate (qldt_berg)
        deallocate (qldt_destr)
        deallocate (qldt_rime)
        deallocate (rain_clr_diag)
        deallocate (rain_cld_diag)
        deallocate (a_rain_clr_diag)
        deallocate (a_rain_cld_diag)
        deallocate (liq_adj)
        deallocate (rain_evap)
        deallocate (qidt_fall)
        deallocate (qidt_fill)
        deallocate (qidt_melt)
        deallocate (qidt_dep)
        deallocate (qidt_subl)
        deallocate (qidt_eros)
        deallocate (qidt_destr)
        deallocate (snow_clr_diag)
        deallocate (snow_cld_diag)
        deallocate (a_snow_clr_diag)
        deallocate (a_snow_cld_diag)
        deallocate (snow_subl)
        deallocate (snow_melt)
        deallocate (ice_adj)
        deallocate (qadt_lsform)
        deallocate (qadt_lsdiss)
        deallocate (qadt_eros)
        deallocate (qadt_fill)
        deallocate (qadt_super)
        deallocate (qadt_rhred)
        deallocate (qadt_destr)
        deallocate (a_precip_cld_diag)
        deallocate (a_precip_clr_diag)
        deallocate (mask3)
     
        end if  !for budget_diag

        
!-----------------------------------------------------------------------
!
!
!       end of subroutine



end subroutine strat_cloud


!#######################################################################
!#######################################################################


! <SUBROUTINE NAME="strat_cloud_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   This writes out a restart (if needed).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call strat_cloud_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine strat_cloud_end()


!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This subroutine writes out radturbten to a restart file.
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!       VARIABLES
!
!
!       -------------------
!       INTERNAL VARIABLES:
!       -------------------
!
!         variable              definition                  unit
!       ------------   -----------------------------   ---------------
!
!       unit           unit number for restart file
!
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!
!       Internal variables
!       ------------------
!

        integer                                :: unit
        character(len=64)                      :: fname='RESTART/strat_cloud.res.nc'

!-----------------------------------------------------------------------
!       
       if(.not.module_is_initialized) then
          return
        else
          module_is_initialized = .false.
        endif
!
        if( do_netcdf_restart) then
           if (mpp_pe() == mpp_root_pe()) then
              call mpp_error ('strat_cloud_mod', 'Writing netCDF formatted restart file: RESTART/strat_cloud.res.nc', NOTE)
           endif
           call write_data (fname, 'vers', restart_versions(size(restart_versions(:))) , no_domain=.true.)
           call write_data (fname, 'nsum', nsum)
           call write_data (fname, 'qlsum', qlsum)
           call write_data (fname, 'qisum', qisum)
           call write_data (fname, 'cfsum', cfsum)
        else
           if (mpp_pe() == mpp_root_pe()) then
              call mpp_error ('strat_cloud_mod', 'Writing native formatted restart file.', NOTE)
           endif
           unit = open_restart_file ('RESTART/strat_cloud.res', ACTION='write')
           if (mpp_pe() == mpp_root_pe()) then
              write (unit) restart_versions(size(restart_versions(:)))
           endif
           call write_data (unit, nsum)
           call write_data (unit, qlsum)
           call write_data (unit, qisum)
           call write_data (unit, cfsum)
           call close_file (unit)
        endif

    !-------------------------------------------
    ! end beta distribution module if used

    if (do_pdf_clouds) call beta_dist_end

end subroutine strat_cloud_end


!#######################################################################
!#######################################################################


! <SUBROUTINE NAME="strat_cloud_sum">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!     This increments cloud variables for passing to radiation.
!     It is expected that this will become obsolete soon.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  strat_cloud_sum (is, js, ql, qi, cf)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!        Starting integer for longitude window.
!  </IN>
!  <IN NAME="js" TYPE="integer">
!        Starting integer for latitude window.
!  </IN>
!  <IN NAME="ql" TYPE="real">
!        Cloud liquid water specific humidity (kg/kg)
!  </IN>
!  <IN NAME="qi" TYPE="real">
!        Cloud ice water specific humidity (kg/kg)
!  </IN>
!  <IN NAME="cf" TYPE="real">
!        Cloud fraction (fraction, 0-1)
!  </IN>
! </SUBROUTINE>
!
 subroutine strat_cloud_sum (is, js, ql, qi, cf)


!-----------------------------------------------------------------------
   integer, intent(in)                   :: is, js
      real, intent(in), dimension(:,:,:) :: ql, qi, cf
!-----------------------------------------------------------------------
   integer :: ie, je

   if(.not.module_is_initialized) then
     call error_mesg('strat_cloud_sum','strat_cloud is not initialized',FATAL)
   endif

   ie = is + SIZE(ql,1) - 1
   je = js + SIZE(ql,2) - 1
     
!--------- use time-averaged or instantaneous clouds -----------

  if (do_average) then
       nsum(is:ie,js:je)   =  nsum(is:ie,js:je)   +  1
       qlsum(is:ie,js:je,:) = qlsum(is:ie,js:je,:) + ql
       qisum(is:ie,js:je,:) = qisum(is:ie,js:je,:) + qi
       cfsum(is:ie,js:je,:) = cfsum(is:ie,js:je,:) + cf
  else
       nsum(is:ie,js:je)   =  1
       qlsum(is:ie,js:je,:) = ql
       qisum(is:ie,js:je,:) = qi
       cfsum(is:ie,js:je,:) = cf
  endif

!-----------------------------------------------------------------------


 end subroutine strat_cloud_sum


!#######################################################################
!#######################################################################


! <SUBROUTINE NAME="strat_cloud_avg">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!      Averaging routine for cloud variables to be passed to radiation.
!      Expected to be removed shortly.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call  strat_cloud_avg (is, js, ql, qi, cf, ierr)
!
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
!      Starting integer for longitude window.
!  </IN>
!  <IN NAME="js" TYPE="integer">
!      Starting integer for latitude window.
!  </IN>
!  <OUT NAME="ql" TYPE="real">
!      Cloud liquid water specific humidity (kg/kg)
!  </OUT>
!  <OUT NAME="qi" TYPE="real">
!      Cloud ice water specific humidity (kg/kg)
!  </OUT>
!  <OUT NAME="cf" TYPE="real">
!      Cloud fraction (0-1)
!  </OUT>
!  <OUT NAME="ierr" TYPE="integer">
!      Error integer.
!  </OUT>
! </SUBROUTINE>
!
 subroutine strat_cloud_avg (is, js, ql, qi, cf, ierr)


!-----------------------------------------------------------------------
   integer, intent(in)                    :: is, js
   real, intent(out), dimension(:,:,:) :: ql, qi, cf
   integer, intent(out)                   :: ierr
!-----------------------------------------------------------------------
   integer :: ie, je, num, k
!-----------------------------------------------------------------------
      
   if(.not.module_is_initialized) then
     call error_mesg('strat_cloud_avg','strat_cloud is not initialized',FATAL)
   endif

   if (SIZE(ql,3) /= SIZE(qlsum,3)) then
     call error_mesg ('strat_cloud_avg in strat_cloud_mod',  &
                              'input argument has the wrong SIZE',FATAL)
   endif

   ie = is + SIZE(ql,1) - 1
   je = js + SIZE(ql,2) - 1
   num = count(nsum(is:ie,js:je) == 0)

   if (num > 0) then

!     ----- no average, return error flag -----

      ierr = 1

   else

!     ----- compute average -----

      do k = 1, SIZE(ql,3)
        ql(:,:,k) = qlsum(is:ie,js:je,k) / float(nsum(is:ie,js:je))
        qi(:,:,k) = qisum(is:ie,js:je,k) / float(nsum(is:ie,js:je))
        cf(:,:,k) = cfsum(is:ie,js:je,k) / float(nsum(is:ie,js:je))
      enddo
      ierr = 0

   endif

    nsum(is:ie,js:je)   = 0
   qlsum(is:ie,js:je,:) = 0.0
   qisum(is:ie,js:je,:) = 0.0
   cfsum(is:ie,js:je,:) = 0.0

!-----------------------------------------------------------------------


 end subroutine strat_cloud_avg


!#######################################################################
!#######################################################################


! <FUNCTION NAME="do_strat_cloud">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!     Logical function to indicate whether or not strat_cloud is running.
!  </DESCRIPTION>
!  <TEMPLATE>
!   result =  do_strat_cloud ( ) result (answer)
!
!  </TEMPLATE>
! </FUNCTION>
!
 function do_strat_cloud ( ) result (answer)


   logical :: answer
   answer = strat_cloud_on


 end function do_strat_cloud


!#######################################################################
!#######################################################################


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppm2m_sak --- Piecewise parabolic method for fields
!
! !INTERFACE:
 subroutine ppm2m_sak(a4, delp, km, kmap, i1, i2, iv, kord)

 implicit none

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kmap    ! partial remap to start
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real, intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! !DESCRIPTION:
!
!   Perform the piecewise parabolic method 
! 
! !REVISION HISTORY: 
!   ??.??.??    Lin        Creation
!   02.04.04    Sawyer     Newest release from FVGCM
! 
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! local arrays:
      real   dc(i1:i2,km)
      real   h2(i1:i2,km)
      real delq(i1:i2,km)
      real  df2(i1:i2,km)
      real   d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt
      integer it
      real fac
      real a1, a2, c1, c2, c3, d1, d2
      real qmax, qmin, cmax, cmin
      real qm, dq, tmp
      real qmp, pmp
      real lac

      km1 = km - 1
       it = i2 - i1 + 1

      do k=max(2,kmap-2),km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo
 
      do k=max(2,kmap-2),km1
         do i=i1,i2
            c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
            c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
             dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
            df2(i,k) = tmp
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=max(3,kmap), km1
      do i=i1,i2
        c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
        a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
        a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
        a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                  ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

      if(km>8 .and. kord>3) call steepz_sak(i1, i2, km, kmap, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      if ( kmap <= 2 ) then
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
         dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
         cmax = max(a4(1,i,1), a4(1,i,2))
         cmin = min(a4(1,i,1), a4(1,i,2))
         a4(2,i,2) = max(cmin,a4(2,i,2))
         a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
! No over- and under-shoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
      enddo

      do k=max(1,kmap),km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo

! Enforce monotonicity of the "slope" within the top layer
      if ( kmap <= 2 ) then
      do i=i1,i2
         if ( a4(2,i,1) * a4(1,i,1) <= 0. ) then 
              a4(2,i,1) = 0.
                dc(i,1) = a4(1,i,1)
         endif
         if ( dc(i,1) * (a4(2,i,2) - a4(1,i,1)) <= 0. ) then
! Setting DC==0 will force piecewise constant distribution after
! calling kmppm_sak
              dc(i,1) = 0.
         endif
      enddo
      endif

! Enforce constraint on the "slope" at the surface

      do i=i1,i2
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) then
!            a4(3,i,km) = 0.
!              dc(i,km) =  -a4(1,i,km)
               dc(i,km) = 0.
         endif
         if( dc(i,km) * (a4(1,i,km) - a4(2,i,km)) <= 0. ) then
             dc(i,km) = 0.
         endif
      enddo
 
!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      if ( kmap <= 2 ) then
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
            call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
      endif

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=max(2,kmap-1), km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2
!           h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
!    &               / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )
!    &               * delp(i,k)**2
! Method#3
            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      if( kord == 7 ) then
         fac = 1.5           ! original quasi-monotone
      else
         fac = 0.125         ! full monotone
      endif

      do k=max(3,kmap), km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord == 7) then
             call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 2)
         endif
      enddo

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

      do k=max(3,kmap), km-2
      if( kord /= 4) then
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
      endif
         call kmppm_sak(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm_sak(dc(i1,k), a4(1,i1,k), it, 0)
      enddo
!EOC
 end subroutine ppm2m_sak
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  kmppm_sak --- Perform piecewise parabolic method in vertical
!
! !INTERFACE:
 subroutine kmppm_sak(dm, a4, itot, lmt)

 implicit none

! !INPUT PARAMETERS:
      real, intent(in):: dm(*)     ! the linear slope
      integer, intent(in) :: itot      ! Total Longitudes
      integer, intent(in) :: lmt       ! 0: Standard PPM constraint
                                       ! 1: Improved full monotonicity constraint (Lin)
                                       ! 2: Positive definite constraint
                                       ! 3: do nothing (return immediately)
! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout) :: a4(4,*)   ! PPM array
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)

! !DESCRIPTION:
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    Incorporated newest FVGCM version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      real, parameter:: r12 = 1./12.
      real qmp
      real da1, da2, a6da
      real fmin
      integer i

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2003)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

!EOC
 end subroutine kmppm_sak
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  steepz_sak --- Calculate attributes for PPM
!
! !INTERFACE:
 subroutine steepz_sak(i1, i2, km, kmap, a4, df2, dm, dq, dp, d4)

   implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: km                   ! Total levels
      integer, intent(in) :: kmap                 ! 
      integer, intent(in) :: i1                   ! Starting longitude
      integer, intent(in) :: i2                   ! Finishing longitude
      real, intent(in) ::  dp(i1:i2,km)       ! grid size
      real, intent(in) ::  dq(i1:i2,km)       ! backward diff of q
      real, intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
      real, intent(in) :: df2(i1:i2,km)       ! first guess mismatch
      real, intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch

! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened

!
! !DESCRIPTION:
!   This is complicated stuff related to the Piecewise Parabolic Method
!   and I need to read the Collela/Woodward paper before documenting
!   thoroughly.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, k
      real alfa(i1:i2,km)
      real    f(i1:i2,km)
      real  rat(i1:i2,km)
      real  dg2

! Compute ratio of dq/dp
      do k=max(2,kmap-1),km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=max(2,kmap-1),km-1
         do i=i1,i2
            f(i,k) =   (rat(i,k+1) - rat(i,k))                          &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=max(3,kmap),km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<0. .and. df2(i,k)/=0.) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0., min(0.5, -0.1875*dg2/df2(i,k))) 
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=max(4,kmap+1),km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

!EOC
 end subroutine steepz_sak
!----------------------------------------------------------------------- 


!#######################################################################
!#######################################################################



end module strat_cloud_mod
