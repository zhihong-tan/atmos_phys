 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               module rad_utilities_mod

use fms_mod,          only : error_mesg, FATAL

use time_manager_mod, only : time_type

!--------------------------------------------------------------------

implicit none
private

character(len=128)  :: version =  '$Id: rad_utilities.F90,v 20.0 2013/12/13 23:21:23 fms Exp $'
character(len=128)  :: tagname =  '$Name: testing $'

!---------------------------------------------------------------------
!-------  interfaces --------

public rad_utilities_init, check_derived_types, & 
       locate_in_table,  &
       looktab, table_alloc, &
       thickavg, thinavg,   &
       rad_utilities_end,   &
       get_radiative_param, assignment(=)

interface looktab
    module procedure  looktab_type1, looktab_type2, looktab_type3
end interface

interface table_alloc
   module procedure    table1_alloc, table2_alloc, table3_alloc
end interface

interface thickavg
   module procedure thickavg_3d
   module procedure thickavg_0d
   module procedure thickavg_1band
   module procedure thickavg_isccp
end interface

interface assignment(=)
  module procedure lw_output_type_eq
  module procedure sw_output_type_eq
  module procedure aerosol_props_type_eq
end interface

!------------------------------------------------------------------

public   aerosol_type              
  
!    aerosol
!    aerosol_names
!    family_members

type aerosol_type
     real, dimension(:,:,:,:),        pointer :: aerosol=>NULL()
     logical, dimension(:,:),        pointer :: family_members=>NULL()
     character(len=64), dimension(:), pointer :: aerosol_names=>NULL()
end type aerosol_type              

!------------------------------------------------------------------
 
public   aerosol_diagnostics_type
 
!    extopdep
!    absopdep
!    extopdep_vlcno
!    absopdep_vlcno
!    lw_extopdep_vlcno
!    lw_absopdep_vlcno
!    sw_heating_vlcno
 
type aerosol_diagnostics_type
!    real, dimension(:,:,:),   pointer  :: sw_heating_vlcno=>NULL()
     real, dimension(:,:,:,:),   pointer  :: sw_heating_vlcno=>NULL()
     real, dimension(:,:,:,:,:), pointer  :: extopdep=>NULL(), &
                                             absopdep=>NULL()
     real, dimension(:,:,:,:,:), pointer  :: asymdep=>NULL()

     real, dimension(:,:,:,:), pointer  :: extopdep_vlcno=>NULL(), &
                                           absopdep_vlcno=>NULL(), &
                                           lw_extopdep_vlcno=>NULL(), &
                                           lw_absopdep_vlcno=>NULL()
 
end type aerosol_diagnostics_type
 
!------------------------------------------------------------------

public   aerosol_properties_type

!    aerextband
!    aerssalbband
!    aerasymmband
!    aerextbandlw
!    aerssalbbandlw
!    sulfate_index
!    optical_index

type aerosol_properties_type
     integer, dimension(:,:,:), pointer  :: ivol=>NULL()
     real, dimension(:,:), pointer  :: aerextband=>NULL(),   &
                                       aerssalbband=>NULL(), &
                                       aerasymmband=>NULL(), &
                                       aerextbandlw=>NULL(), &
                                       aerssalbbandlw=>NULL(), &
                                       aerextbandlw_cn=>NULL(), &
                                       aerssalbbandlw_cn=>NULL()
     real, dimension(:,:,:,:), pointer :: sw_ext=>NULL(), &
                                          sw_ssa=>NULL(), &
                                          sw_asy=>NULL(), &
                                          lw_ext=>NULL(), &
                                          lw_ssa=>NULL(), &
                                          lw_asy=>NULL()
!yim
     integer, dimension(:,:), pointer :: sulfate_index=>NULL()
     integer, dimension(:), pointer :: optical_index=>NULL()
     integer, dimension(:), pointer :: omphilic_index=>NULL()
     integer, dimension(:), pointer :: bcphilic_index=>NULL()
     integer, dimension(:), pointer :: seasalt1_index=>NULL()
     integer, dimension(:), pointer :: seasalt2_index=>NULL()
     integer, dimension(:), pointer :: seasalt3_index=>NULL()
     integer, dimension(:), pointer :: seasalt4_index=>NULL()
     integer, dimension(:), pointer :: seasalt5_index=>NULL()
     integer                        :: sulfate_flag
     integer                        :: omphilic_flag
     integer                        :: bcphilic_flag
     integer                        :: seasalt1_flag
     integer                        :: seasalt2_flag
     integer                        :: seasalt3_flag
     integer                        :: seasalt4_flag
     integer                        :: seasalt5_flag
!yim
     integer                        :: bc_flag
end type aerosol_properties_type

!------------------------------------------------------------------

public astronomy_type
    
!    solar
!    cosz
!    fracday
!    rrsun

type astronomy_type
     real, dimension(:,:), pointer  :: solar=>NULL(),   &
                                       cosz=>NULL(),  &
                                       fracday=>NULL()
     real, dimension(:,:,:), pointer  :: solar_p=>NULL(),   &
                                       cosz_p=>NULL(),  &
                                       fracday_p=>NULL()
     real    :: rrsun
end type astronomy_type

!--------------------------------------------------------------------
 
public astronomy_inp_type
 
!    zenith_angle   specified zenith angles [ degrees ]
!    fracday        specified daylight fraction [ fraction ]
!    rrsun          specified earth-sun distance, normalized by mean
!                   distance 

type astronomy_inp_type
     real, dimension(:,:), pointer  :: zenith_angle=>NULL()
     real, dimension(:,:), pointer  :: fracday=>NULL()
     real                           :: rrsun
end type astronomy_inp_type

!--------------------------------------------------------------------

public atmos_input_type

!    press
!    temp
!    rh2o
!    zfull
!    pflux
!    tflux
!    deltaz
!    phalf
!    rel_hum
!    cloudtemp
!    clouddeltaz
!    cloudvapor
!    aerosoltemp
!    aerosolvapor
!    aerosolpress
!    aerosolrelhum
!    tracer_co2
!    g_rrvco2
!    tsfc
!    psfc

type atmos_input_type
     real, dimension(:,:,:), pointer :: press=>NULL(),   &
                                        temp=>NULL(), &
                                        rh2o=>NULL(),  &
                                        zfull=>NULL(),  &
                                        pflux=>NULL(), &
                                        tflux=>NULL(),  &
                                        deltaz=>NULL(),  &
                                        phalf=>NULL(),   &
                                        rel_hum=>NULL(), &
                                        cloudtemp=>NULL(),   &
                                        clouddeltaz=>NULL(), &
                                        cloudvapor=>NULL(), &
                                        aerosoltemp=>NULL(), &
                                        aerosolvapor=>NULL(), &
                                        aerosolpress=>NULL(), &
                                        aerosolrelhum=>NULL(), &
                                        tracer_co2 => NULL()
     real, dimension(:,:),   pointer :: tsfc=>NULL(),   &
                                        psfc=>NULL()              
     real                            :: g_rrvco2
end type atmos_input_type

!-------------------------------------------------------------------

public cldrad_properties_type

!    cldext
!    cldasymm
!    cldsct
!    emmxolw
!    emrndlw
!    abscoeff
!    cldemiss
!    cirabsw
!    cirrfsw
!    cvisrfsw

type cldrad_properties_type
     real, dimension(:,:,:,:,:), pointer :: cldext=>NULL(),   &
                                          cldasymm=>NULL(), &
                                          cldsct=>NULL()
     real, dimension(:,:,:,:,:), pointer ::                   &
                                          emmxolw=>NULL(), &
                                          emrndlw=>NULL(),  &
                                          abscoeff=>NULL(),  &
                                          cldemiss=>NULL()
     real, dimension(:,:,:), pointer ::   cirabsw=>NULL(), &
                                          cirrfsw=>NULL(),   &
                                          cvisrfsw=>NULL()
end type cldrad_properties_type

!------------------------------------------------------------------

public cld_space_properties_type

!    camtswkc
!    cirabswkc
!    cirrfswkc
!    cvisrfswkc
!    ktopswkc
!    kbtmswkc

type cld_space_properties_type
     real, dimension(:,:,:),    pointer :: camtswkc=>NULL()        
     real, dimension(:,:,:),    pointer :: cirabswkc=>NULL(),  &
                                           cirrfswkc=>NULL(), &
                                           cvisrfswkc=>NULL()
     integer, dimension(:,:,:), pointer :: ktopswkc=>NULL(),   &
                                           kbtmswkc=>NULL()
end type cld_space_properties_type

!------------------------------------------------------------------

public cld_specification_type

!    tau
!    lwp
!    iwp
!    reff_liq
!    reff_ice
!    liq_frac
!    cloud_water
!    cloud_ice
!    cloud_area
!    cloud_droplet
!    reff_liq_micro
!    reff_ice_micro
!    camtsw
!    cmxolw
!    crndlw
!    cld_thickness
!    ncldsw
!    nmxolw
!    nrndlw
!    hi_cloud
!    mid_cloud
!    low_cloud
!    ice_cloud

type cld_specification_type
   real, dimension(:,:,:,:),  pointer :: tau=>NULL(),  &
                                         camtsw_band=>NULL(), &
                                         crndlw_band=>NULL(), &
                                         lwp_lw_band=>NULL(), &
                                         iwp_lw_band=>NULL(), &
                                         lwp_sw_band=>NULL(), &
                                         iwp_sw_band=>NULL(), &
                                         reff_liq_lw_band=>NULL(),   &
                                         reff_ice_lw_band=>NULL(), &
                                         reff_liq_sw_band=>NULL(),   &
                                         reff_ice_sw_band=>NULL()
   real, dimension(:,:,:),    pointer :: lwp=>NULL(),   &
                                         iwp=>NULL(),  &
                                         reff_liq=>NULL(),   &
                                         reff_ice=>NULL(), &
                                         reff_liq_lim=>NULL(),   &
                                         reff_ice_lim=>NULL(), &
                                         liq_frac=>NULL(), &
                                         cloud_water=>NULL(), &
                                         cloud_ice=>NULL(),  &
                                         cloud_area=>NULL(), &
                                         cloud_droplet=>NULL(), &
                                         cloud_ice_num=>NULL(), &
  ! snow, rain
                                         rain =>NULL(), &
                                         snow =>NULL(), &
                                         rain_size =>NULL(), &
                                         snow_size =>NULL(), &

                                         reff_liq_micro=>NULL(),   &
                                         reff_ice_micro=>NULL(),&
                                         camtsw=>NULL(),   &
                                         cmxolw=>NULL(),  &
                                         crndlw=>NULL()
   integer, dimension(:,:,:), pointer :: cld_thickness=>NULL()
   integer, dimension(:,:,:,:), pointer :: stoch_cloud_type=>NULL()
   integer, dimension(:,:,:,:), pointer :: cld_thickness_lw_band=>NULL()
   integer, dimension(:,:,:,:), pointer :: cld_thickness_sw_band=>NULL()
   integer, dimension(:,:),   pointer :: ncldsw=>NULL(),   &
                                         nmxolw=>NULL(),&
                                         nrndlw=>NULL()
   integer, dimension(:,:,:), pointer :: ncldsw_band=>NULL(),   &
                                         nrndlw_band=>NULL()
   logical, dimension(:,:,:), pointer :: hi_cloud=>NULL(),   &
                                         mid_cloud=>NULL(),  &
                                         low_cloud=>NULL(),   &
                                         ice_cloud=>NULL()
end type cld_specification_type

!------------------------------------------------------------------

public cloudrad_control_type

type cloudrad_control_type
    logical :: do_pred_cld_microphys
    logical :: do_presc_cld_microphys
    logical :: do_bulk_microphys
    logical :: do_sw_micro
    logical :: do_lw_micro
    logical :: do_rh_clouds        
    logical :: do_strat_clouds        
    logical :: do_zonal_clouds        
    logical :: do_mgroup_prescribed
    logical :: do_obs_clouds        
    logical :: do_no_clouds        
    logical :: do_diag_clouds        
    logical :: do_specified_clouds        
    logical :: do_donner_deep_clouds
    logical :: do_uw_clouds
    logical :: do_zetac_clouds
    logical :: do_random_overlap
    logical :: do_max_random_overlap
    logical :: do_stochastic_clouds
    logical :: use_temp_for_seed
    logical :: do_specified_strat_clouds
    logical :: do_ica_calcs
    logical :: do_liq_num
    logical :: do_ice_num
    logical :: using_fu2007
    integer :: nlwcldb                   !   number of frequency bands 
                                         !   for which lw cloud emissiv-
                                         !   ities are defined.
    integer :: cloud_data_points
    integer :: ich
    integer :: icm
    integer :: ict
    integer :: icb
    logical :: do_pred_cld_microphys_iz
    logical :: do_presc_cld_microphys_iz
    logical :: do_bulk_microphys_iz
    logical :: do_sw_micro_iz
    logical :: do_lw_micro_iz
    logical :: do_rh_clouds_iz
    logical :: do_strat_clouds_iz
    logical :: do_zonal_clouds_iz
    logical :: do_mgroup_prescribed_iz
    logical :: do_obs_clouds_iz
    logical :: do_no_clouds_iz
    logical :: do_diag_clouds_iz
    logical :: do_specified_clouds_iz
    logical :: do_donner_deep_clouds_iz
    logical :: do_uw_clouds_iz
    logical :: do_zetac_clouds_iz
    logical :: do_random_overlap_iz
    logical :: do_max_random_overlap_iz
    logical :: do_stochastic_clouds_iz
    logical :: use_temp_for_seed_iz
    logical :: do_specified_strat_clouds_iz
    logical :: do_ica_calcs_iz
    logical :: do_liq_num_iz
    logical :: do_ice_num_iz
    logical :: using_fu2007_iz
end type cloudrad_control_type

!------------------------------------------------------------------

public fsrad_output_type

!    tdtsw
!    tdtlw
!    tdtsw_clr
!    tdtlw_clr
!    swdns
!    swups
!    lwdns
!    lwups
!    swin
!    swout
!    olr
!    swdns_clr
!    swups_clr
!    lwdns_clr
!    lwups_clr`
!    swin_clr
!    swout_clr
!    olr_clr
!    npass

type fsrad_output_type
     real, dimension(:,:,:), pointer :: tdtsw=>NULL(), &
                                        tdtlw=>NULL(),  &
                                        tdtsw_clr=>NULL(),  &
                                        tdtlw_clr=>NULL()
     real, dimension(:,:),   pointer :: swdns=>NULL(),   &
                                        swups=>NULL(),  &
                                        lwups=>NULL(), &
                                        lwdns=>NULL(), &
                                        swin=>NULL(), &
                                        swout=>NULL(), &
                                        olr=>NULL(), &
                                        swdns_clr=>NULL(),  &
                                        swups_clr=>NULL(),  &
                                        lwups_clr=>NULL(),&
                                        lwdns_clr=>NULL(),   &
                                        swin_clr=>NULL(),  &
                                        swout_clr=>NULL(), &
                                        olr_clr=>NULL()
     integer      :: npass
end type fsrad_output_type

!-------------------------------------------------------------------

public gas_tf_type

!    tdav
!    tlsqu
!    tmpdiff
!    tstdav
!    co2nbl
!    n2o9c
!    tn2o17
!    co2spnb
!    a1
!    a2

type gas_tf_type
     real, dimension(:,:,:),   pointer :: tdav=>NULL(),   &
                                          tlsqu=>NULL(),   &
                                          tmpdiff=>NULL(),   &
                                          tstdav=>NULL(),  &
                                          co2nbl=>NULL(),   &
                                          n2o9c=>NULL(),   &
                                          tn2o17=>NULL()
     real, dimension(:,:,:,:), pointer :: co2spnb=>NULL()
     real, dimension(:,:),     pointer :: a1=>NULL(),    &
                                          a2=>NULL()
end type gas_tf_type

!------------------------------------------------------------------

public longwave_control_type 

type longwave_control_type
    character(len=16) :: lw_form
    character(len=16) :: continuum_form
    character(len=16) :: linecatalog_form
    logical           :: do_cfc
    logical           :: do_lwaerosol
    logical           :: do_ch4
    logical           :: do_n2o
    logical           :: do_ch4lbltmpint
    logical           :: do_n2olbltmpint
    logical           :: do_co2
    logical           :: do_lwcldemiss
    logical           :: do_h2o
    logical           :: do_o3 
    logical           :: do_cfc_iz
    logical           :: do_lwaerosol_iz
    logical           :: do_ch4_iz
    logical           :: do_n2o_iz
    logical           :: do_ch4lbltmpint_iz
    logical           :: do_n2olbltmpint_iz
    logical           :: do_co2_iz
    logical           :: do_lwcldemiss_iz
    logical           :: do_h2o_iz
    logical           :: do_o3_iz 
end type longwave_control_type

!---------------------------------------------------------------------

public longwave_parameter_type

type longwave_parameter_type
     integer   :: offset
     integer   :: NBTRG
     integer   :: NBTRGE
     integer   :: NBLY
     integer   :: n_lwaerosol_bands
     real      :: lw_band_resolution
     logical   :: offset_iz
     logical   :: NBTRG_iz
     logical   :: NBTRGE_iz
     logical   :: NBLY_iz
     logical   :: n_lwaerosol_bands_iz
     logical   :: lw_band_resolution_iz
end type longwave_parameter_type

!--------------------------------------------------------------------

public longwave_tables1_type

!    vae
!    td
!    md
!    cd

type longwave_tables1_type
    real, dimension(:,:), pointer  ::  vae=>NULL(),   &
                                       td=>NULL(), &
                                       md=>NULL(), &
                                       cd=>NULL()
end type longwave_tables1_type

!--------------------------------------------------------------------

public longwave_tables2_type

!    vae
!    td
!    md
!    cd

type longwave_tables2_type
    real, dimension(:,:,:), pointer  ::  vae=>NULL(),  &
                                         td=>NULL(),  &
                                         md=>NULL(),   &
                                         cd=>NULL()
end type longwave_tables2_type

!---------------------------------------------------------------------

public longwave_tables3_type

!    vae
!    td

type longwave_tables3_type
     real,  dimension(:,:), pointer    ::  vae=>NULL(),   &
                                           td=>NULL()          
end type longwave_tables3_type

!---------------------------------------------------------------------

public lw_clouds_type

!    taucld_rndlw
!    taucld_mxolw
!    taunbl_mxolw

type lw_clouds_type
     real, dimension(:,:,:,:),   pointer :: taucld_rndlw=>NULL(), &
                                            taucld_mxolw=>NULL(), &
                                            taunbl_mxolw=>NULL()
end type lw_clouds_type

!------------------------------------------------------------------

public lw_diagnostics_type

!    flx1e1
!    gxcts
!    flx1e1f
!    excts
!    fctsg
!    fluxn
!    fluxncf
!    exctsn
!    cts_out
!    cts_outcf

type lw_diagnostics_type
     real, dimension(:,:),   pointer   :: flx1e1=>NULL(),  &
                                          gxcts=>NULL()
     real, dimension(:,:,:), pointer   :: flx1e1f=>NULL(),  &
                                          excts=>NULL(),&
                                          fctsg=>NULL()
     real, dimension(:,:,:,:), pointer :: fluxn=>NULL(),   &
                                          fluxncf=>NULL(),   &
                                          exctsn=>NULL(),  &
                                          cts_out=>NULL(), &
                                          cts_outcf=>NULL()
end type lw_diagnostics_type

!-------------------------------------------------------------------

public lw_output_type

!    heatra
!    flxnet
!    heatracf
!    flxnetcf
!    netlw_special
!    netlw_special_clr

type lw_output_type
     real, dimension(:,:,:), pointer :: heatra=>NULL(), &
                                        flxnet=>NULL(),  &
                                        heatracf=>NULL(), &
                                        flxnetcf=>NULL()
     real, dimension(:,:,:), pointer   :: netlw_special=>NULL(), &
                                          netlw_special_clr=>NULL(), &
                                          bdy_flx=>NULL(), &
                                          bdy_flx_clr=>NULL()
end type lw_output_type

!------------------------------------------------------------------

public lw_table_type

!    bdlocm
!    bdhicm
!    bandlo
!    bandhi
!    iband

type lw_table_type
     real, dimension(:),    pointer :: bdlocm=>NULL(),   &
                                       bdhicm=>NULL(),  &
                                       bandlo=>NULL(),  &
                                       bandhi=>NULL()
     integer, dimension(:), pointer :: iband=>NULL()
end type lw_table_type

!------------------------------------------------------------------

public microphysics_type
 
!    conc_ice
!    conc_drop
!    size_ice
!    size_drop
!    size_snow
!    conc_snow
!    size_rain
!    conc_rain
!    cldamt
!    stoch_conc_ice
!    stoch_conc_drop
!    stoch_size_ice
!    stoch_size_drop
!    stoch_cldamt
!    stoch_cloud_type
!    lw_stoch_conc_ice
!    lw_stoch_conc_drop
!    lw_stoch_size_ice
!    lw_stoch_size_drop
!    lw_stoch_cldamt
!    sw_stoch_conc_ice
!    sw_stoch_conc_drop
!    sw_stoch_size_ice
!    sw_stoch_size_drop
!    sw_stoch_cldamt

type microphysics_type
   real, dimension(:,:,:), pointer :: conc_ice=>NULL(),   &
                                      conc_drop=>NULL(),      &
                                      size_ice=>NULL(),   &
                                      size_drop=>NULL(),     &
                                      size_snow=>NULL(),   &
                                      conc_snow=>NULL(),     &
                                      size_rain=>NULL(),     &
                                      conc_rain=>NULL(),   &
                                      cldamt=>NULL(),      &
                                      droplet_number=>NULL(), &
                                      ice_number=>NULL()
real, dimension(:,:,:,:), pointer :: stoch_conc_ice=>NULL(),   &
                                     stoch_conc_drop=>NULL(),  &
                                     stoch_size_ice=>NULL(),   &
                                     stoch_size_drop=>NULL(),  &
                                     stoch_cldamt=>NULL(),     &
                                     stoch_droplet_number=>NULL(), &
                                     stoch_ice_number=>NULL()
integer, dimension(:,:,:,:), pointer ::  stoch_cloud_type=>NULL()
!
! In practice, we allocate a single set of columns for the stochastic
!   clouds, then point to sections of the larger array with the 
!   lw_ and sw_type. 
!   I.e. lw_stoch_conc_ice => stoch_conc_ice(:, :, :, 1:numLwBands)
!
real, dimension(:,:,:,:), pointer :: lw_stoch_conc_ice=>NULL(),   &
                                     lw_stoch_conc_drop=>NULL(),  &
                                     lw_stoch_size_ice=>NULL(),   &
                                     lw_stoch_size_drop=>NULL(),  &
                                     lw_stoch_cldamt=>NULL(),     &
                                     lw_stoch_droplet_number=>NULL(), &
                                     lw_stoch_ice_number=>NULL(), &
                                     sw_stoch_conc_ice=>NULL(),   &
                                     sw_stoch_conc_drop=>NULL(),  &
                                     sw_stoch_size_ice=>NULL(),   &
                                     sw_stoch_size_drop=>NULL(),  &
                                     sw_stoch_cldamt=>NULL(),     &
                                     sw_stoch_droplet_number=>NULL(), &
                                     sw_stoch_ice_number=>NULL()
end type microphysics_type

!-------------------------------------------------------------------

public microrad_properties_type
 
!    cldext
!    cldsct
!    cldasymm
!    abscoeff

type microrad_properties_type
   real, dimension(:,:,:,:), pointer :: cldext=>NULL(),  &
                                        cldsct=>NULL(), &
                                        cldasymm=>NULL(),    &
                                        abscoeff=>NULL()
end type microrad_properties_type

!--------------------------------------------------------------------

public optical_path_type

!    empl1f
!    vrpfh2o
!    xch2obd
!    tphfh2o
!    avephif
!    totaerooptdep
!    empl1
!    empl2
!    var1
!    var2
!    emx1f
!    emx2f
!    totvo2
!    avephi
!    totch2obdwd
!    xch2obdwd
!    totphi
!    cntval
!    toto3
!    tphio3
!    var3
!    var4
!    wk
!    rh2os
!    rfrgn
!    tfac
!    totaerooptdep15
!    totf11
!    totf12
!    totf113
!    totf22
!    emx1
!    emx2
!    csfah2o
!    aerooptdep_KE_15

type optical_path_type
     real, dimension (:,:,:,:), pointer :: empl1f=>NULL(),  &
                                           empl2f=>NULL(),  &
                                           vrpfh2o=>NULL(), &
                                           xch2obd=>NULL(),  &
                                           tphfh2o=>NULL(), &
                                           avephif=>NULL(), &
                                           totaerooptdep=>NULL()
     real, dimension (:,:,:),   pointer :: empl1=>NULL(), &
                                           empl2=>NULL(),  &
                                           var1=>NULL(), &
                                           var2=>NULL(), &
                                           emx1f=>NULL(),   &
                                           emx2f=>NULL(),   &
                                           totvo2=>NULL(),  &
                                           avephi=>NULL(),&
                                           totch2obdwd=>NULL(), &
                                           xch2obdwd=>NULL(), &
                                           totphi=>NULL(),   &
                                           cntval=>NULL(), &
                                           toto3=>NULL(),   &
                                           tphio3=>NULL(),  &
                                           var3=>NULL(),  &
                                           var4=>NULL(),        &
                                           wk=>NULL(),         &
                                           rh2os=>NULL(),  &
                                           rfrgn=>NULL(),  &
                                           tfac=>NULL(), &
                                           totaerooptdep_15=>NULL(), &
                                           totf11=>NULL(),   &
                                           totf12=>NULL(),  &
                                           totf113=>NULL(),   &
                                           totf22=>NULL()
      real, dimension (:,:), pointer    :: emx1=>NULL(),  &
                                           emx2=>NULL(),  &
                                           csfah2o=>NULL(), &
                                           aerooptdep_KE_15=>NULL()
end type optical_path_type

!------------------------------------------------------------------

public radiation_control_type

type radiation_control_type
    logical  :: do_totcld_forcing
    logical  :: do_aerosol
    integer  :: rad_time_step
    integer  :: lw_rad_time_step
    integer  :: sw_rad_time_step
    logical  :: do_sw_rad
    logical  :: do_lw_rad
    logical  :: hires_coszen
    integer  :: nzens
    real     :: co2_tf_calc_intrvl
    logical  :: use_current_co2_for_tf
    logical  :: calc_co2_tfs_on_first_step
    logical  :: calc_co2_tfs_monthly
    real     :: co2_tf_time_displacement
    real     :: ch4_tf_calc_intrvl
    logical  :: use_current_ch4_for_tf
    logical  :: calc_ch4_tfs_on_first_step
    logical  :: calc_ch4_tfs_monthly
    real     :: ch4_tf_time_displacement
    real     :: n2o_tf_calc_intrvl
    logical  :: use_current_n2o_for_tf
    logical  :: calc_n2o_tfs_on_first_step
    logical  :: calc_n2o_tfs_monthly
    real     :: n2o_tf_time_displacement
    integer  :: mx_spec_levs
    logical  :: time_varying_solar_constant
    logical  :: volcanic_sw_aerosols
    logical  :: volcanic_lw_aerosols
    logical  :: using_solar_timeseries_data
    logical  :: do_lwaerosol_forcing
    logical  :: do_swaerosol_forcing
    integer  :: indx_swaf
    integer  :: indx_lwaf
    logical  :: using_im_bcsul
    logical  :: do_totcld_forcing_iz
    logical  :: do_aerosol_iz
    logical  :: rad_time_step_iz
    logical  :: lw_rad_time_step_iz
    logical  :: sw_rad_time_step_iz
    logical  :: do_sw_rad_iz
    logical  :: do_lw_rad_iz
    logical  :: hires_coszen_iz
    logical  :: nzens_iz  
    logical  :: co2_tf_calc_intrvl_iz
    logical  :: use_current_co2_for_tf_iz
    logical  :: calc_co2_tfs_on_first_step_iz
    logical  :: calc_co2_tfs_monthly_iz
    logical  :: ch4_tf_calc_intrvl_iz
    logical  :: use_current_ch4_for_tf_iz
    logical  :: calc_ch4_tfs_on_first_step_iz
    logical  :: calc_ch4_tfs_monthly_iz
    logical  :: n2o_tf_calc_intrvl_iz
    logical  :: use_current_n2o_for_tf_iz
    logical  :: calc_n2o_tfs_on_first_step_iz
    logical  :: calc_n2o_tfs_monthly_iz
    logical  :: co2_tf_time_displacement_iz
    logical  :: ch4_tf_time_displacement_iz
    logical  :: n2o_tf_time_displacement_iz
    logical  :: mx_spec_levs_iz
    logical  :: time_varying_solar_constant_iz
    logical  :: volcanic_sw_aerosols_iz
    logical  :: volcanic_lw_aerosols_iz
    logical  :: using_solar_timeseries_data_iz
    logical  :: do_lwaerosol_forcing_iz
    logical  :: do_swaerosol_forcing_iz
    logical  :: indx_swaf_iz
    logical  :: indx_lwaf_iz
    logical  :: using_im_bcsul_iz
end type radiation_control_type

!------------------------------------------------------------------

public   radiative_gases_type
 
!    qo3
!    rrvch4
!    rrvn2o
!    rrvco2
!    rrvf11
!    rrvf12
!    rrvf113
!    rrvf22
!    rf11air
!    rf12air
!    rf113air
!    rf22air
!    time_varying_co2
!    time_varying_f11
!    time_varying_f12
!    time_varying_f113
!    time_varying_f22
!    time_varying_ch4
!    time_varying_n2o

type radiative_gases_type
     real, dimension(:,:,:), pointer :: qo3=>NULL()
     real                            :: rrvch4, rrvn2o, rrvco2,    &
                                        rrvf11, rrvf12, rrvf113,  &
                                        rrvf22, rf11air, rf12air,  &
                                        rf113air, rf22air, &
                                        co2_for_last_tf_calc,  &
                                        co2_tf_offset, &
                                        co2_for_next_tf_calc, &
                                        ch4_for_last_tf_calc,  &
                                        ch4_tf_offset, &
                                        ch4_for_next_tf_calc, &
                                        n2o_for_last_tf_calc,  &
                                        n2o_tf_offset, &
                                        n2o_for_next_tf_calc
     logical                         :: time_varying_co2,  &
                                        time_varying_f11, &
                                        time_varying_f12,  &
                                        time_varying_f113, &
                                        time_varying_f22,  &
                                        time_varying_ch4, &
                                        time_varying_n2o, &
                                        use_model_supplied_co2
     type(time_type)                 :: Co2_time, Ch4_time, N2o_time
end type radiative_gases_type

!------------------------------------------------------------------

public rad_output_type

!    tdt_rad
!    tdt_rad_clr
!    tdtsw
!    tdtsw_clr
!    tdtlw
!    flux_sw_surf
!    flux_lw_surf
!    coszen_angle

type rad_output_type
     real, dimension(:,:,:,:), pointer :: tdt_rad=>NULL(),  &
                                        ufsw=>NULL(),  &
                                        dfsw=>NULL(),  &
                                        tdtsw=>NULL()  
     real, dimension(:,:,:,:), pointer :: tdt_rad_clr=>NULL(), &
                                        ufsw_clr=>NULL(),  &
                                        dfsw_clr=>NULL(),  &
                                        tdtsw_clr=>NULL()
                                        
     real, dimension(:,:,:), pointer :: tdtlw=>NULL()
     real, dimension(:,:,:), pointer :: flxnet=>NULL()
     real, dimension(:,:,:), pointer :: flxnetcf=>NULL()
     real, dimension(:,:,:), pointer :: tdtlw_clr=>NULL()
     real, dimension(:,:,:),   pointer :: flux_sw_surf=>NULL(), &
                                        flux_sw_surf_refl_dir=>NULL(), &
                                        flux_sw_surf_dir=>NULL(), &
                                        flux_sw_surf_dif=>NULL(), &
                                        flux_sw_down_vis_dir=>NULL(), &
                                        flux_sw_down_vis_dif=>NULL(), &
                                       flux_sw_down_total_dir=>NULL(), &
                                       flux_sw_down_total_dif=>NULL(), &
                                        flux_sw_vis=>NULL(), &
                                        flux_sw_vis_dir=>NULL(), &
                                        flux_sw_refl_vis_dir=>NULL(), &
                                        flux_sw_vis_dif=>NULL()
     real, dimension(:,:,:),   pointer :: flux_sw_down_vis_clr=>NULL(), &
                                  flux_sw_down_total_dir_clr=>NULL(), &
                                  flux_sw_down_total_dif_clr=>NULL()
     real, dimension(:,:),   pointer :: flux_lw_surf=>NULL(), &
                                        coszen_angle=>NULL()
end type rad_output_type

!-------------------------------------------------------------------

public shortwave_control_type

type shortwave_control_type
    logical  :: do_lhsw
    logical  :: do_esfsw
    logical  :: do_swaerosol
    logical  :: do_diurnal
    logical  :: do_annual
    logical  :: do_daily_mean
    logical  :: do_cmip_diagnostics
    real     :: solar_constant
    logical  :: do_lhsw_iz
    logical  :: do_esfsw_iz
    logical  :: do_swaerosol_iz
    logical  :: do_diurnal_iz
    logical  :: do_annual_iz
    logical  :: do_daily_mean_iz
    logical  :: do_cmip_diagnostics_iz
end type shortwave_control_type

!---------------------------------------------------------------------

public solar_spectrum_type

!    solarfluxtoa      highly-resolved solar flux at toa in 
!                      Sw_control%tot_wvnums 
!                      bands [         ]
!    solflxband_lean   a time series of toa solar flux in each
!                      parameterization band
!    solflxband_lean_ann_1882   1882 average toa solar flux in each
!                               parameterization band
!    solflxband_lean_ann_2000   2000 average toa solar flux in each
!                               parameterization band
!    solflxband        toa solar flux in each parameterization band
!                      [         ]
!    endwvnbands       highest wave number in each of the solar
!                      spectral parameterization bands [ cm (-1) ]
!    tot_wvnums
!    nbands
!    nfrqpts
!    nstreams
!    nh2obands
!    visible_band_indx
!    one_micron_indx
!    eight70_band_indx
!    visible_band_indx_iz
!    one_micron_indx_iz
!    eight70_band_indx_iz

type solar_spectrum_type
    real, dimension(:),    pointer   :: solarfluxtoa=>null()
    real, dimension(:),    pointer   :: solflxband=>NULL()
    real, dimension(:),    pointer   :: solflxbandref=>NULL()
    real, dimension(:),    pointer   :: solflxband_lean_ann_1882=>NULL()
    real, dimension(:),    pointer   :: solflxband_lean_ann_2000=>NULL()
    real, dimension(:,:,:),pointer   :: solflxband_lean=>NULL()
    integer, dimension(:), pointer   :: endwvnbands=>NULL()
    integer         :: tot_wvnums
    integer         :: nbands
    integer         :: nfrqpts
    integer         :: nstreams
    integer         :: nh2obands
    integer         :: visible_band_indx, one_micron_indx
    integer         :: eight70_band_indx
    logical         :: visible_band_indx_iz, one_micron_indx_iz
    logical         :: eight70_band_indx_iz
    integer         :: w340_band_indx, w380_band_indx,  &
                       w440_band_indx, w670_band_indx
    logical         :: w340_band_iz, w380_band_iz, &
                       w440_band_iz, w670_band_iz
end type solar_spectrum_type

!---------------------------------------------------------------------

public surface_type

!    asfc
!    land

type surface_type
    real, dimension(:,:),   pointer ::  asfc=>NULL(),   &
                                        land=>NULL(),  &
                                        asfc_vis_dir=>NULL(), &
                                        asfc_nir_dir=>NULL(), &
                                        asfc_vis_dif=>NULL(), &
                                        asfc_nir_dif=>NULL()
end type surface_type
 
!-------------------------------------------------------------------

public sw_output_type

!    dfsw
!    ufsw
!    fsw
!    hsw
!    dfswcf
!    ufswcf
!    fswcf
!    hswcf
!    swdn_special
!    swup_special
!    swdn_special_clr
!    swup_special_clr

type sw_output_type
     real, dimension(:,:,:,:), pointer :: dfsw=>NULL(),   &
                                        ufsw=>NULL(),  &
                                        fsw=>NULL(),   &
                                        hsw=>NULL()   
     real, dimension(:,:,:,:), pointer :: dfswcf=>NULL(),   &
                                        ufswcf=>NULL(),&
                                        fswcf=>NULL(),  &
                                        hswcf=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_vis_sfc=>NULL(),   &
                                       ufsw_vis_sfc=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_dir_sfc=>NULL()
      real, dimension(:,:,:), pointer :: ufsw_dir_sfc=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_dir_sfc_clr=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_dif_sfc=>NULL(),   &
                                       ufsw_dif_sfc=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_dif_sfc_clr=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_vis_sfc_dir=>NULL()
      real, dimension(:,:,:), pointer :: ufsw_vis_sfc_dir=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_vis_sfc_clr=>NULL()
      real, dimension(:,:,:), pointer :: dfsw_vis_sfc_dif=>NULL(),   &
                                       ufsw_vis_sfc_dif=>NULL()
      real, dimension(:,:,:,:), pointer   ::  bdy_flx=>NULL()
      real, dimension(:,:,:,:), pointer   ::  bdy_flx_clr=>NULL()
      real, dimension(:,:,:,:), pointer   ::                       &
                                        swup_special=>NULL(), &
                                        swup_special_clr=>NULL()
     real, dimension(:,:,:,:), pointer   :: swdn_special=>NULL(), &
                                          swdn_special_clr=>NULL()
end type sw_output_type

!-------------------------------------------------------------------

public table_axis_type

type table_axis_type
  integer :: first_col
  real    :: min_val
  real    :: max_val
  real    :: tab_inc
end type table_axis_type

!---------------------------------------------------------------------

!private      &


!--------------------------------------------------------------------
!-------- namelist  ---------

integer  ::  dummy = 0


namelist / rad_utilities_nml /  dummy


!---------------------------------------------------------------------
!------- public data ------


type (longwave_control_type),  public   ::    &
     Lw_control = longwave_control_type( '    ', '    ', '    ', &
                                         .false., .false., .false.,  &
                                         .false., .false., .false.,  &
                                         .false., .false.,  &
                                         .false., .false.,  &
                                         .false., .false., .false.,  &
                                         .false., .false.,  &
                                         .false., .false.,  &
                                         .false., .false., .false.   )

type (shortwave_control_type), public   ::  &
    Sw_control = shortwave_control_type( .false., .false., .false. , &
                                         .false., .false., .false., &
                                         .false., &
                                         0.0, &
                                         .false., .false., .false. , &
                                         .false., &
                                         .false., .false., .false.)

type (radiation_control_type), public   ::  &
   Rad_control = radiation_control_type( .false., .false., 0, 0, 0, &
                                         .false., .false., &
                                         .false., 1, &
                                         0.0,  .true.,  .false.,&
                                         .false.,  0.0, &
                                         0.0, .true., .false.,  &
                                         .false.,  0.0,  &
                                         0.0, .true., .false.,   &
                                         .false., 0.0, &
                                         0, .false., .false.,   &
                                         .false., .false.,&
                                         .false., .false.,      &
                                         0, 0, .false., &
! _iz variables:
                                         .false., .false., .false., &
                                         .false., .false., .true., &
                                         .true.,  &
                                         .false., .false., &
                                         .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., &
                                         .false., .false., &
                                         .false., .false., &
                                         .false.,          &
                                         .false., .false., .false.,  &
                                         .false., .false.,   &
                                         .false., .false., &
                                         .false., .false., .false.)

type (cloudrad_control_type), public    ::   &
 Cldrad_control = cloudrad_control_type( .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false.,                   &
                                         0,0,0,0,0,0 , &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false., .false., .false., &
                                         .false.   )


type (longwave_parameter_type), public  ::   &
Lw_parameters = longwave_parameter_type( 0, 0, 0, 0, 0, 10.0, &
                                         .false., .false., .false.,  &
                                         .false., .false., .true.)

type (table_axis_type),        public   ::    &
               temp_1 = table_axis_type( 1, 100.0, 370.0, 10.0), &
               mass_1 = table_axis_type( 1, -16.0,   1.9,  0.1)

!---------------------------------------------------------------------

                           contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine rad_utilities_init

call error_mesg ('rad_utilities_mod', &
  'This is the null version of rad_utilities_init and should not be executed', FATAL)

end subroutine rad_utilities_init

!#####################################################################

subroutine check_derived_types 

call error_mesg ('rad_utilities_mod', &
  'This is the null version of check_derived_types and should not be executed', FATAL)

end subroutine check_derived_types 

!#####################################################################

subroutine locate_in_table (table_axis, x, dx, ix, k_min, k_max) 

type(table_axis_type),     intent(in)  :: table_axis
real,    dimension(:,:,:), intent(in)  :: x
integer,                   intent(in)  :: k_min, k_max
real,    dimension(:,:,:), intent(out) :: dx
integer, dimension(:,:,:), intent(out) :: ix


call error_mesg ('rad_utilities_mod', &
  'This is the null version of locate_in_table and should not be executed', FATAL)

end subroutine locate_in_table

!####################################################################

subroutine looktab_type1 (tab, ix, iy, dx, dy, answer, k_min, k_max)   

type(longwave_tables1_type), intent(in)  :: tab
integer,dimension(:,:,:),    intent(in)  :: ix, iy
real,   dimension(:,:,:),    intent(in)  :: dx, dy   
real,   dimension(:,:,:),    intent(out) :: answer
integer,                     intent(in)  :: k_min, k_max

call error_mesg ('rad_utilities_mod', &
  'This is the null version of looktab and should not be executed', FATAL)

end subroutine looktab_type1

!#####################################################################

subroutine looktab_type2 (tab, ix, iy, dx, dy, answer, k_min, k_max, m)

type(longwave_tables2_type), intent(in)   :: tab
integer, dimension (:,:,:),  intent(in)   :: ix, iy
integer,                     intent(in)   :: m
real, dimension (:,:,:),     intent(in)   :: dx, dy
real, dimension (:,:,:),     intent(out)  :: answer
integer,                     intent(in)   :: k_min, k_max

call error_mesg ('rad_utilities_mod', &
  'This is the null version of looktab and should not be executed', FATAL)

end subroutine looktab_type2

!###################################################################

subroutine looktab_type3 (tab, ix, dx,  answer, k_min, k_max, n)

type(longwave_tables3_type), intent(in)  :: tab
integer, dimension (:,:,:),  intent(in)  :: ix
integer,                     intent(in)  :: n
real,    dimension(:,:,:),   intent(in)  :: dx
real,    dimension(:,:,:),   intent(out) :: answer 
integer,                     intent(in)  :: k_min, k_max

call error_mesg ('rad_utilities_mod', &
  'This is the null version of looktab and should not be executed', FATAL)

end subroutine  looktab_type3

!#####################################################################

subroutine table1_alloc (tab, dim1, dim2)

type(longwave_tables1_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2

call error_mesg ('rad_utilities_mod', &
  'This is the null version of table_alloc and should not be executed', FATAL)

end subroutine table1_alloc

!####################################################################

subroutine table2_alloc (tab, dim1, dim2, dim3)

type(longwave_tables2_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2, dim3

call error_mesg ('rad_utilities_mod', &
  'This is the null version of table_alloc and should not be executed', FATAL)

end subroutine table2_alloc

!#####################################################################

subroutine table3_alloc (tab, dim1, dim2)

type(longwave_tables3_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2

call error_mesg ('rad_utilities_mod', &
  'This is the null version of table_alloc and should not be executed', FATAL)

end subroutine table3_alloc

!##################################################################

subroutine thickavg_3d (nivl1, nivl2, nivls, nbands, extivl, ssalbivl,&
                        asymmivl, solflxivl, solflxband, mask, extband,&
                        ssalbband, asymmband)

integer, dimension(:),    intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:,:,:,:), intent(in)       :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real, dimension(:),       intent(in)       :: solflxband            
real, dimension(:,:,:,:), intent(out)      :: extband, ssalbband, asymmband
logical, dimension(:,:,:), intent(in)      :: mask

call error_mesg ('rad_utilities_mod', &
  'This is the null version of thickavg and should not be executed', FATAL)
  
end subroutine thickavg_3d

!####################################################################

subroutine thickavg_0d (nivl1, nivl2, nivls, nbands, extivl, ssalbivl,&
                        asymmivl, solflxivl, solflxband, extband,  &
                        ssalbband , asymmband)
 
integer, dimension(:),    intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:),       intent(in)       :: extivl, asymmivl
real, dimension(:),       intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real, dimension(:),       intent(in)       :: solflxband            
real, dimension(:),       intent(out)      :: extband, ssalbband, asymmband

call error_mesg ('rad_utilities_mod', &
  'This is the null version of thickavg and should not be executed', FATAL)

end subroutine thickavg_0d

!##################################################################

subroutine thickavg_isccp (nband, nivl1, nivl2, extivl,          &
                           solflxivl, solflxband, mask, extband)

integer,                  intent(in)       :: nband
integer,                  intent(in)       :: nivl1, nivl2
real, dimension(:,:,:,:), intent(in)       :: extivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real,                     intent(in)       :: solflxband            
logical, dimension(:,:,:),intent(in)       :: mask
real, dimension(:,:,:),   intent(out)      :: extband

call error_mesg ('rad_utilities_mod', &
  'This is the null version of thickavg and should not be executed', FATAL)

end subroutine thickavg_isccp

!##################################################################

subroutine thickavg_1band (nband, nivl1, nivl2, nivls, nbands, extivl, &
                           ssalbivl, asymmivl, solflxivl, solflxband, &
                           mask, extband, ssalbband, asymmband)
 
integer,                  intent(in)       :: nband
integer,                  intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:,:,:,:), intent(in)       :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real,                     intent(in)       :: solflxband            
real, dimension(:,:,:  ), intent(inout)      :: extband, ssalbband,   &
                                              asymmband
logical, dimension(:,:,:), intent(in)      :: mask

call error_mesg ('rad_utilities_mod', &
  'This is the null version of thickavg and should not be executed', FATAL)
  
end subroutine thickavg_1band

!####################################################################

subroutine thinavg (nivl1, nivl2, nivls, nbands, extivl, ssalbivl, &
                    asymmivl,  solflxivl, solflxband, extband,   &
                    ssalbband , asymmband)
 
integer, dimension(:),    intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
integer,                  intent(in)       :: nbands
real, dimension(:,:,:,:), intent(in)       :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real, dimension(:),       intent(in)       :: solflxband            
real, dimension(:,:,:,:), intent(out)      :: extband, ssalbband,   &
                                              asymmband
call error_mesg ('rad_utilities_mod', &
  'This is the null version of thinavg and should not be executed', FATAL)

end subroutine thinavg 

!#########################################################################
subroutine get_radiative_param(text_in_scheme,text_in_param, &
                               rad_forc_online, tr_rad_name,  &
                               tr_clim_name, tr_rad_scale_factor)

character(len=*), intent(in)    :: text_in_scheme, text_in_param
logical, intent(out)            :: rad_forc_online
character(len=*), intent(out)   :: tr_rad_name,tr_clim_name
real,             intent(out)   :: tr_rad_scale_factor

call error_mesg ('rad_utilities_mod', &
  'This is the null version of get_radiative_param and should not be executed', FATAL)

end subroutine get_radiative_param

!#####################################################################

subroutine rad_utilities_end


end subroutine rad_utilities_end

!#####################################################################

subroutine aerosol_props_type_eq(aerosol_props_out,aerosol_props_in)

   type(aerosol_properties_type), intent(inout) :: aerosol_props_out
   type(aerosol_properties_type), intent(in)    :: aerosol_props_in

call error_mesg ('rad_utilities_mod', &
  'This is the null version of aerosol_props_type_eq and should not be executed', FATAL)
   
end subroutine aerosol_props_type_eq

!#####################################################################

subroutine lw_output_type_eq(lw_output_out,lw_output_in)

   type(lw_output_type), intent(inout) :: lw_output_out
   type(lw_output_type), intent(in)    :: lw_output_in

call error_mesg ('rad_utilities_mod', &
  'This is the null version of lw_output_type_eq and should not be executed', FATAL)

end subroutine lw_output_type_eq

!#####################################################################

subroutine sw_output_type_eq(sw_output_out,sw_output_in)

   type(sw_output_type), intent(inout) :: sw_output_out
   type(sw_output_type), intent(in)    :: sw_output_in

call error_mesg ('rad_utilities_mod', &
  'This is the null version of sw_output_type_eq and should not be executed', FATAL)

end subroutine sw_output_type_eq

!####################################################################

                     end module rad_utilities_mod
