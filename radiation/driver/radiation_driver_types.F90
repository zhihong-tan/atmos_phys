 
module radiation_driver_types_mod

!
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!  smf
! </REVIEWER>
! 
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Code to define the derived data types.
! </OVERVIEW>
! <DESCRIPTION>
!  This code is used in the radiation code package as a helper module.
!  It defines many derived data types used in radiation calculation.
! </DESCRIPTION>
!
!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    radiation_driver_types_mod contains definitions for
!    derived-type variables used by the radiation/driver routines.
!---------------------------------------------------------------------

public astronomy_type

!    solar     shortwave flux factor
!    cosz      cosine of zenith angle
!    fracday   fraction of timestep during which the sun is shining
!    rrsun     inverse of square of earth-sun distance
!              relative to the mean square of earth-sun distance

type astronomy_type
     real, dimension(:,:), pointer  :: solar=>NULL(),   &
                                       cosz=>NULL(),  &
                                       fracday=>NULL()
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
!    relhum
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
                                        relhum=>NULL(), &
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

public radiation_control_type

type radiation_control_type
    logical  :: do_totcld_forcing
    logical  :: using_restart_file
    logical  :: do_sw_rad
    logical  :: do_lw_rad
    logical  :: renormalize_sw_fluxes
  ! from shortwave_control_type
    logical  :: do_esfsw
    logical  :: do_diurnal
    logical  :: do_annual
    logical  :: do_daily_mean
    logical  :: do_cmip_sw_diagnostics
end type radiation_control_type

!------------------------------------------------------------------

public rad_output_type

!
!  rad_output_type variable with the following components:
!
!      tdt_rad                radiative (sw + lw) heating rate
!      flux_sw_surf           net (down-up) sw flux at surface
!      flux_sw_surf_dir       net (down-up) sw flux at surface
!      flux_sw_surf_refl_dir       dir sw flux reflected at surface
!      flux_sw_surf_dif            net (down-up) sw flux at surface
!      flux_sw_down_vis_dir        downward visible sw flux at surface
!      flux_sw_down_vis_dif        downward visible sw flux at surface
!      flux_sw_down_total_dir      downward total sw flux at surface
!      flux_sw_down_total_dif      downward total sw flux at surface
!      flux_sw_down_total_dir_clr  downward total direct sw flux at surface  (clear sky)
!      flux_sw_down_total_dif_clr  downward total diffuse sw flux at surface (clear sky)
!      flux_sw_down_vis_clr        downward visible sw flux at surface  (clear sky)
!      flux_sw_vis            net visible sw flux at surface
!      flux_sw_vis_dir        net visible sw flux at surface
!      flux_sw_refl_vis_dir   reflected direct visible sw flux at surface
!      flux_sw_vis_dif        net visible sw flux at surface
!      flux_lw_surf           downward lw flux at surface
!      coszen_angle    cosine of the zenith angle (used for the last radiation calculation)
!      tdt_rad_clr     net radiative heating rate in the absence of cloud
!      tdtsw           shortwave heating rate
!      tdtsw_clr       shortwave heating rate in he absence of cloud
!      tdtlw_clr       longwave heating rate in he absence of cloud
!      tdtlw           longwave heating rate
!      ufsw            upward sw flux
!      dfsw            downward sw flux
!      ufsw_clr        upward sw flux
!      dfsw_clr        downward sw flux
!      flxnet          net lw flux
!      flxnetcf        net lw flux, cloud free
!      extinction      SW extinction (band 4 centered on 1 micron) for volcanoes

type rad_output_type
     real, dimension(:,:,:), pointer :: tdt_rad=>NULL(),  &
                                        ufsw=>NULL(),  &
                                        dfsw=>NULL(),  &
                                        tdtsw=>NULL()  
     real, dimension(:,:,:), pointer :: tdt_rad_clr=>NULL(), &
                                        ufsw_clr=>NULL(),  &
                                        dfsw_clr=>NULL(),  &
                                        tdtsw_clr=>NULL()
                                        
     real, dimension(:,:,:), pointer :: tdtlw=>NULL()
     real, dimension(:,:,:), pointer :: flxnet=>NULL()
     real, dimension(:,:,:), pointer :: flxnetcf=>NULL()
     real, dimension(:,:,:), pointer :: tdtlw_clr=>NULL()
     real, dimension(:,:),   pointer :: flux_sw_surf=>NULL(), &
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
     real, dimension(:,:),   pointer :: flux_sw_down_vis_clr=>NULL(), &
                                  flux_sw_down_total_dir_clr=>NULL(), &
                                  flux_sw_down_total_dif_clr=>NULL()
     real, dimension(:,:),   pointer :: flux_lw_surf=>NULL(), &
                                        coszen_angle=>NULL()
     real, dimension(:,:,:), pointer :: extinction=>NULL()
end type rad_output_type

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

!####################################################################


end module radiation_driver_types_mod


