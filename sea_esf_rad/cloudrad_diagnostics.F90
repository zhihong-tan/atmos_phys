                 module cloudrad_diagnostics_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    cloudrad_diagnostics_mod generates any desired netcdf output
!    fields involving the cloud fields seen by the radiation package
!    or the cloud radiation interaction variables.
! </OVERVIEW>
! <DESCRIPTION>
!    cloudrad_diagnostics_mod generates any desired netcdf output
!    fields involving the cloud fields seen by the radiation package
!    or the cloud radiation interaction variables.
! </DESCRIPTION>

! shared modules:

use fms_mod,                 only: fms_init, open_namelist_file, &
                                   write_version_number, mpp_pe, &
                                   mpp_root_pe, stdlog, file_exist,  &
                                   check_nml_error, error_mesg,   &
                                   FATAL, NOTE, WARNING, close_file,  &
                                   read_data, write_data
use time_manager_mod,        only: time_type, time_manager_init
use diag_manager_mod,        only: register_diag_field, send_data, &
                                   diag_manager_init

! shared radiation package modules:

use rad_utilities_mod,       only: rad_utilities_init, Environment, &
                                   cldrad_properties_type, &
                                   cld_specification_type, &
                                   microrad_properties_type, &
                                   microphysics_type, atmos_input_type,&
                                   Cldrad_control

!  other cloud diagnostics modules

use cloud_rad_mod,           only: get_strat_cloud_diagnostics, &
                                   cloud_rad_init
use isccp_clouds_mod,        only: isccp_clouds_init, isccp_clouds_end,&
                                   isccp_output, tau_reff_diag2,  &
                                   isccp_cloudtypes

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    cloudrad_diagnostics_mod generates any desired netcdf output
!    fields involving the cloud fields seen by the radiation package
!    or the cloud radiation interaction variables.
!
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: cloudrad_diagnostics.F90,v 10.0 2003/10/24 22:00:39 fms Exp $'
character(len=128)  :: tagname =  '$Name: jakarta $'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
         cloudrad_diagnostics_init, cloudrad_netcdf, &
         cloudrad_diagnostics_end

private          &
!   called from cloudrad_diagnostics_init:
         diag_field_init, &
!   called from cloudrad_netcdf:
         isccp_diag, compute_isccp_clds,  &
!   called from isccp_diag:  
         cloud_optical_properties_diag


!---------------------------------------------------------------------
!-------- namelist  ---------

integer  ::  dummy = 0

 namelist /cloudrad_diagnostics_nml /                             &
                                       dummy


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

real, parameter     :: taumin = 1.E-06  ! minimum value allowed for 
                                        ! optical depth 
                                        ! [ dimensionless ]
real                :: qmin             ! minimum permissible cloud  
                                        ! condensate 
                                        ! [ kg condensate / kg air ]
integer             :: overlap          ! variable indicating which 
                                        ! overlap assumption to use:
                                        ! overlap = 1. means condensate
                                        ! in adjacent levels is treated
                                        ! as part of the same cloud
                                        ! i.e. maximum-random overlap
                                        ! overlap = 2. means condensate
                                        ! in adjacent levels is treated 
                                        ! as different clouds
                                        ! i.e. random overlap

!----------------------------------------------------------------------
!    diagnostics variables.     
!----------------------------------------------------------------------
character(len=8)    :: mod_name = 'cloudrad'
real                :: missing_value = -999.

integer :: id_tot_cld_amt, id_cld_amt, id_em_cld_lw, id_em_cld_10u, & 
           id_abs_lsc_cld_10u, id_abs_lsc_cld_lw,  &
           id_abs_cell_cld_10u, id_abs_cell_cld_lw,  &
           id_abs_meso_cld_10u, id_abs_meso_cld_lw,  &
           id_abs_cld_10u, id_abs_cld_lw,  &
           id_high_cld_amt, id_mid_cld_amt, id_low_cld_amt,  &
           id_lsc_cld_amt, id_cell_cld_amt, id_meso_cld_amt,  &
           id_lsc_cld_ext_uv, id_lsc_cld_ext_vis, id_lsc_cld_ext_nir, &
           id_lsc_cld_sct_uv, id_lsc_cld_sct_vis, id_lsc_cld_sct_nir, &
           id_lsc_cld_asymm_uv, id_lsc_cld_asymm_vis,    &
           id_lsc_cld_asymm_nir, &
           id_cell_cld_ext_uv, id_cell_cld_ext_vis,    &
           id_cell_cld_ext_nir, &
           id_cell_cld_sct_uv, id_cell_cld_sct_vis,    &
           id_cell_cld_sct_nir, &
           id_cell_cld_asymm_uv, id_cell_cld_asymm_vis,    &
           id_cell_cld_asymm_nir, &
           id_meso_cld_ext_uv, id_meso_cld_ext_vis,   &
           id_meso_cld_ext_nir, &
           id_meso_cld_sct_uv, id_meso_cld_sct_vis,   &
           id_meso_cld_sct_nir, &
           id_meso_cld_asymm_uv, id_meso_cld_asymm_vis,    &
           id_meso_cld_asymm_nir, &
           id_ext_cld_uv,   id_sct_cld_uv,  id_asymm_cld_uv, &
           id_ext_cld_vis,  id_sct_cld_vis, id_asymm_cld_vis, &
           id_ext_cld_nir,  id_sct_cld_nir, id_asymm_cld_nir, &
           id_alb_uv_cld, id_alb_nir_cld, id_abs_uv_cld, id_abs_nir_cld

!   ISCCP diagnostic variables:

integer :: id_pc1tau1,id_pc1tau2,id_pc1tau3,id_pc1tau4,id_pc1tau5, &
           id_pc1tau6,id_pc1tau7, &
           id_pc2tau1,id_pc2tau2,id_pc2tau3,id_pc2tau4,id_pc2tau5, &
           id_pc2tau6,id_pc2tau7, &
           id_pc3tau1,id_pc3tau2,id_pc3tau3,id_pc3tau4,id_pc3tau5, &
           id_pc3tau6,id_pc3tau7, &
           id_pc4tau1,id_pc4tau2,id_pc4tau3,id_pc4tau4,id_pc4tau5, &
           id_pc4tau6,id_pc4tau7, &
           id_pc5tau1,id_pc5tau2,id_pc5tau3,id_pc5tau4,id_pc5tau5, &
           id_pc5tau6,id_pc5tau7, &
           id_pc6tau1,id_pc6tau2,id_pc6tau3,id_pc6tau4,id_pc6tau5, &
           id_pc6tau6,id_pc6tau7, &
           id_pc7tau1,id_pc7tau2,id_pc7tau3,id_pc7tau4,id_pc7tau5, &
           id_pc7tau6,id_pc7tau7, &
           id_nisccp, id_aice, id_reffice, id_aliq, id_reffliq, &
           id_alow, id_tauicelow, id_tauliqlow, id_tlaylow, id_tcldlow

logical :: do_isccp    = .false.     ! are isccp diagnostics desired ?
logical :: do_tau_reff = .false.     ! are the isccp summary diagnostics
                                    ! desired ?
logical :: do_sunlit   = .false.    ! do isccp diagnostics only in day-
                                    ! light ?


logical :: module_is_initialized =                            &
                         .false.    ! module  initialized ?


!----------------------------------------------------------------------
!----------------------------------------------------------------------



                        contains 



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#####################################################################
! <SUBROUTINE NAME="cloudrad_diagnostics_init">
!  <OVERVIEW>
!    cloudrad_diagnostics_init is the constructor for 
!    cloudrad_diagnostics_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloudrad_diagnostics_init is the constructor for 
!    cloudrad_diagnostics_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_diagnostics_init (axes, Time)
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
!
subroutine cloudrad_diagnostics_init (axes, Time)

!---------------------------------------------------------------------
!    cloudrad_diagnostics_init is the constructor for 
!    cloudrad_diagnostics_mod.
!------------------------------------------------------------------

integer, dimension(4),   intent(in)    ::   axes
type(time_type),         intent(in)    ::   Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       axes      diagnostic variable axes
!       Time      current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer         :: unit, io, ierr

!---------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      io       error status returned from io operation  
!      ierr     error code
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call rad_utilities_init
      call time_manager_init
      if (Environment%running_gcm .or. &
          Environment%running_sa_model) then
        call diag_manager_init
      endif

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=cloudrad_diagnostics_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'cloudrad_diagnostics_nml')
        enddo
10      call close_file (unit)
      endif
 
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() )    &
                       write (stdlog(), nml=cloudrad_diagnostics_nml)
 

!-------------------------------------------------------------------
!    initialize the netcdf diagnostics provided with this module.
!-------------------------------------------------------------------
      if ( (Environment%running_gcm .or.    &
            Environment%running_sa_model)   .and. &
            .not. Cldrad_control%do_no_clouds) then
        call diag_field_init (Time, axes)
      endif

!---------------------------------------------------------------------
!    if strat_cloud_mod is active, verify that that module has been 
!    activated and retrieve the overlap parameter and the value used 
!    for the minimum cloud water amount.
!---------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then
        call cloud_rad_init
        call get_strat_cloud_diagnostics (overlap, qmin)
        if (EnvironmenT%running_gcm .or.   &
            Environment%running_sa_model) then
          call isccp_clouds_init (axes, Time, do_isccp_out=do_isccp,   &
                                  do_tau_reff_out=do_tau_reff, &
                                  do_sunlit_out=do_sunlit,   &
                                  overlap_in=overlap, qmin_in=qmin)
        endif 
      endif 

!--------------------------------------------------------------------
!    mark the module initialized.
!--------------------------------------------------------------------
      module_is_initialized= .true.

!--------------------------------------------------------------------



end subroutine cloudrad_diagnostics_init


!###################################################################
! <SUBROUTINE NAME="cloudrad_netcdf">
!  <OVERVIEW>
!    cloudrad_netcdf generates and outputs netcdf fields describing the
!    cloud radiative properties, along with isccp cloud diagnostics
!    fields. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloudrad_netcdf generates and outputs netcdf fields describing the
!    cloud radiative properties, along with isccp cloud diagnostics
!    fields. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_netcdf (is, js, Time_diag, Atmos_input, cosz, &
!                            Lsc_microphys, Meso_microphys, &
!                            Cell_microphys, Lscrad_props,   &
!                            Mesorad_props, Cellrad_props, Cldrad_props,&
!                            Cld_spec, mask)
!  </TEMPLATE>
!  <IN NAME="is,js" TYPE="integer">
!   starting subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Time_diag" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                        diagnostic output [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </IN>
!  <IN NAME="cosz" TYPE="real">
!    cosine of solar zenith angle
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </IN>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale 
!                        clouds
!  </IN>
!  <IN NAME="Meso_microphys" TYPE="microphysics_type">
!   microphysical specification for meso-scale 
!                        clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cell_microphys" TYPE="microphysics_type">
!   microphysical specification for convective cell
!                        clouds associated with donner convection
!  </IN>
!  <IN NAME="Lscrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the large-scale 
!                      clouds   
!  </IN>
!  <IN NAME="Mesorad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the meso-scale
!                      clouds assciated with donner convection
!  </IN>
!  <IN NAME="Cellrad_props" TYPE="microrad_properties_type">
!   cloud radiative properties for the convective cell
!                      clouds associated with donner convection 
!  </IN>
!  <IN NAME="Cldrad_props" TYPE="cldrad_properties_type">
!   cloud radiative properties on model grid
!  </IN>
!  <IN NAME="mask" TYPE="real">
!   OPTIONAL: present when running eta vertical coordinate,
!                        mask to remove points below ground
!  </IN>
! </SUBROUTINE>
!
subroutine cloudrad_netcdf (is, js, Time_diag, Atmos_input, cosz, &
                            Lsc_microphys, Meso_microphys, &
                            Cell_microphys, Lscrad_props,   &
                            Mesorad_props, Cellrad_props, Cldrad_props,&
                            Cld_spec, mask)

!---------------------------------------------------------------------
!    cloudrad_netcdf generates and outputs netcdf fields describing the
!    cloud radiative properties, along with isccp cloud diagnostics
!    fields. 
!---------------------------------------------------------------------

integer,                        intent(in)           :: is, js
type(time_type),                intent(in)           :: Time_diag
type(atmos_input_type),         intent(in)           :: Atmos_input
real, dimension(:,:),           intent(in)           :: cosz        
type(microphysics_type),        intent(in)           :: Lsc_microphys, &
                                                        Meso_microphys,&
                                                        Cell_microphys
type(microrad_properties_type), intent(in)           :: Lscrad_props, &
                                                        Mesorad_props, &
                                                        Cellrad_props
type(cldrad_properties_type),   intent(in)           :: Cldrad_props
type(cld_specification_type),   intent(in)           :: Cld_spec       
real, dimension(:,:,:),         intent(in), optional :: mask

!-------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting subdomain i,j indices of data 
!                      in the physics_window being integrated
!      Time_diag       time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!      Atmos_input     atmospheric input fields on model grid,
!                      [ atmos_input_type ]
!      cosz            cosine of zenith angle --  mean value over
!                      appropriate averaging interval
!                      [ non-dimensional ]
!      Lsc_microphys   microphysical specification for large-scale 
!                      clouds
!                      [ microphysics_type ]
!      Meso_microphys  microphysical specification for meso-scale 
!                      clouds assciated with donner convection
!                      [ microphysics_type ]
!      Cell_microphys  microphysical specification for convective cell
!                      clouds associated with donner convection
!                      [ microphysics_type ]
!      Lscrad_props    cloud radiative properties for the large-scale 
!                      clouds   
!                      [ microrad_properties_type ]
!      Mesorad_props   cloud radiative properties for meso-scale 
!                      clouds associated with donner convection   
!                      [ microrad_properties_type ]
!      Cellrad_props   cloud radiative properties for convective cell
!                      clouds associated with donner convection  
!                      [ microrad_properties_type ]
!      Cldrad_props    total-cloud radiative properties,
!                      [ cldrad_properties_type ]
!      Cld_spec        variables on the model grid which define the 
!                      cloud location and amount     
!                      [ cld_specification_type ]
!
!   intent(in), optional variables:
!
!      mask              present when running eta vertical coordinate,
!                        mask to remove points below ground
!
!------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2),                       &
                      size(Atmos_input%rh2o,3))     :: cloud

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2))     :: tca       

      real, dimension(size(Atmos_input%rh2o,1),                       &
                      size(Atmos_input%rh2o,2), 3)  :: hml_ca        

      logical    ::  used
      integer    ::  kx
      integer    ::  i, j, k   

!---------------------------------------------------------------------
!  local variables:
!
!      cloud                array used to hold the various netcdf 
!                           output fields as they are sent off to 
!                           diag_manager_mod
!      tca                  total column cloud amount [ dimensionless ]
!      hml_ca               total column cloud amount for isccp high, 
!                           middle and low clouds, individually
!      used                 flag returned from send_data indicating
!                           whether diag_manager_mod has received 
!                           data that was sent
!      kx                   number of model layers
!      i,j,k                do-loop indices
!
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('cloudrad_diagnostics_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
!    define the number of model layers.
!---------------------------------------------------------------------
      kx =  size(Cld_spec%camtsw,3)

!--------------------------------------------------------------------
!    when running strat clouds, call isccp_diag to generate isccp-
!    relevant diagnostics.
!---------------------------------------------------------------------
      if (Environment%running_gcm .or.     &
          Environment%running_sa_model) then
        if (Cldrad_control%do_strat_clouds) then
          call isccp_diag (is, js, Cld_spec, Atmos_input, cosz,    &
                          Time_diag)
        endif
      endif

!---------------------------------------------------------------------
!    if desired as a diagnostic, define the total cloud amount. send to
!    diag_manager_mod for netcdf output.
!---------------------------------------------------------------------
      if ( id_tot_cld_amt > 0 ) then
        tca = 1.0
        do k=1,kx        
          tca(:,:) = tca(:,:)*(1.0 - Cld_spec%camtsw(:,:,k))
        end do
        tca = 100.*(1. - tca)
        used = send_data (id_tot_cld_amt, tca, Time_diag, is, js)
      endif

!---------------------------------------------------------------------
!    if high, mid or low cloud diagnostics are desired, call 
!    compute_isccp_clds to define the amount of each. send to 
!    diag_manager_mod.
!---------------------------------------------------------------------
      if (id_high_cld_amt > 0 .or. id_mid_cld_amt > 0 .or. &
          id_low_cld_amt > 0) then
        call compute_isccp_clds (Atmos_input%pflux, Cld_spec%camtsw, &
                                 hml_ca)
        if (id_high_cld_amt > 0)  used =    &
           send_data (id_high_cld_amt, hml_ca(:,:,1), Time_diag, is, js)
        if (id_mid_cld_amt > 0)  used =     &
           send_data (id_mid_cld_amt, hml_ca(:,:,2), Time_diag, is, js)
        if (id_low_cld_amt > 0)  used =    &
           send_data (id_low_cld_amt, hml_ca(:,:,3), Time_diag, is, js)
      endif

!----------------------------------------------------------------------
!    send the 3D cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
      if (id_cld_amt > 0)  used =        &
        send_data (id_cld_amt, Cld_spec%camtsw, Time_diag, is, js, 1, &
                   rmask=mask)

!----------------------------------------------------------------------
!    the following diagnostics are meaningful only when strat_clouds
!    is active:
!----------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then

!----------------------------------------------------------------------
!    send the 3D large-scale cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
 
        if (id_lsc_cld_amt > 0)  used =      &
           send_data (id_lsc_cld_amt, Lsc_microphys%cldamt, Time_diag, &
                      is, js, 1, rmask=mask)

!----------------------------------------------------------------------
!    send various large-scale cloud shortwave radiative property fields
!    to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_lsc_cld_ext_uv > 0) then
          cloud(:,:,:) = Lscrad_props%cldext(:,:,:,22)
          used = send_data (id_lsc_cld_ext_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_ext_vis > 0) then
          cloud(:,:,:) = Lscrad_props%cldext(:,:,:,12)
          used = send_data (id_lsc_cld_ext_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_ext_nir > 0) then
          cloud(:,:,:) = Lscrad_props%cldext(:,:,:,8)
          used = send_data (id_lsc_cld_ext_nir, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_sct_uv > 0) then
          cloud(:,:,:) = Lscrad_props%cldsct(:,:,:,22)
          used = send_data (id_lsc_cld_sct_uv, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_sct_vis > 0) then
          cloud(:,:,:) = Lscrad_props%cldsct(:,:,:,12)
          used = send_data (id_lsc_cld_sct_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_sct_nir > 0) then
          cloud(:,:,:) = Lscrad_props%cldsct(:,:,:,8)
          used = send_data (id_lsc_cld_sct_nir, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_asymm_uv > 0) then
          cloud(:,:,:) = Lscrad_props%cldasymm(:,:,:,22)
          used = send_data (id_lsc_cld_asymm_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_asymm_vis > 0) then
          cloud(:,:,:) = Lscrad_props%cldasymm(:,:,:,12)
          used = send_data (id_lsc_cld_asymm_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_lsc_cld_asymm_nir > 0) then
          cloud(:,:,:) = Lscrad_props%cldasymm(:,:,:,8)
          used = send_data (id_lsc_cld_asymm_nir, cloud, Time_diag,  &
                            is, js, 1, rmask=mask)
        endif
      endif ! (do_strat_clouds)

!----------------------------------------------------------------------
!    the following diagnostics are meaningful only when 
!    donner_deep_clouds is active:
!----------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then

!----------------------------------------------------------------------
!    send the 3D cell-scale cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_cell_cld_amt > 0 ) then
          used = send_data (id_cell_cld_amt, Cell_microphys%cldamt, &
                            Time_diag, is, js, 1, rmask=mask)
        endif

!----------------------------------------------------------------------
!    send various cell-scale cloud shortwave radiative property fields
!    to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_cell_cld_ext_uv > 0) then
          cloud(:,:,:) = Cellrad_props%cldext(:,:,:,22)
          used = send_data (id_cell_cld_ext_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_ext_vis > 0) then
          cloud(:,:,:) = Cellrad_props%cldext(:,:,:,12)
          used = send_data (id_cell_cld_ext_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_ext_nir > 0) then
          cloud(:,:,:) = Cellrad_props%cldext(:,:,:,8)
          used = send_data (id_cell_cld_ext_nir, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_sct_uv > 0) then
          cloud(:,:,:) = Cellrad_props%cldsct(:,:,:,22)
          used = send_data (id_cell_cld_sct_uv, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_sct_vis > 0) then
          cloud(:,:,:) = Cellrad_props%cldsct(:,:,:,12)
          used = send_data (id_cell_cld_sct_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_sct_nir > 0) then
          cloud(:,:,:) = Cellrad_props%cldsct(:,:,:,8)
          used = send_data (id_cell_cld_sct_nir, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_cell_cld_asymm_uv > 0) then
          cloud(:,:,:) = Cellrad_props%cldasymm(:,:,:,22)
          used = send_data (id_cell_cld_asymm_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if ( id_cell_cld_asymm_vis > 0 ) then
          cloud(:,:,:) = Cellrad_props%cldasymm(:,:,:,12)
          used = send_data (id_cell_cld_asymm_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if ( id_cell_cld_asymm_nir > 0 ) then
          cloud(:,:,:) = Cellrad_props%cldasymm(:,:,:,8)
          used = send_data (id_cell_cld_asymm_nir, cloud, Time_diag,  &
                            is, js, 1, rmask=mask )
        endif

!----------------------------------------------------------------------
!    send the 3D meso-scale cloud amount field to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_meso_cld_amt > 0) then
          used = send_data (id_meso_cld_amt, Meso_microphys%cldamt,  &
                            Time_diag, is, js, 1, rmask=mask)
        endif

!----------------------------------------------------------------------
!    send various meso-scale cloud shortwave radiative property fields
!    to diag_manager_mod.
!----------------------------------------------------------------------
        if (id_meso_cld_ext_uv > 0) then
          cloud(:,:,:) = Mesorad_props%cldext(:,:,:,22)
          used = send_data (id_meso_cld_ext_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_ext_vis > 0) then
          cloud(:,:,:) = Mesorad_props%cldext(:,:,:,12)
          used = send_data (id_meso_cld_ext_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_ext_nir > 0) then
          cloud(:,:,:) = Mesorad_props%cldext(:,:,:,8)
          used = send_data (id_meso_cld_ext_nir, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_sct_uv > 0) then
          cloud(:,:,:) = Mesorad_props%cldsct(:,:,:,22)
          used = send_data (id_meso_cld_sct_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_sct_vis > 0) then
          cloud(:,:,:) = Mesorad_props%cldsct(:,:,:,12)
          used = send_data (id_meso_cld_sct_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_sct_nir > 0) then
          cloud(:,:,:) = Mesorad_props%cldsct(:,:,:,8)
          used = send_data (id_meso_cld_sct_nir, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_asymm_uv > 0) then
          cloud(:,:,:) = Mesorad_props%cldasymm(:,:,:,22)
          used = send_data (id_meso_cld_asymm_uv, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_asymm_vis > 0) then
          cloud(:,:,:) = Mesorad_props%cldasymm(:,:,:,12)
          used = send_data (id_meso_cld_asymm_vis, cloud, Time_diag,  &
                            is, js, 1, rmask=mask)
        endif
        if (id_meso_cld_asymm_nir > 0) then
          cloud(:,:,:) = Mesorad_props%cldasymm(:,:,:,8)
          used = send_data (id_meso_cld_asymm_nir, cloud, Time_diag,  &
                            is, js, 1, rmask=mask)
        endif
      endif ! (do_donner_deep_clouds)

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    a total cloud field emissivity that is the weighted average
!    of the random and max overlap emissivities, over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_em_cld_10u > 0) then
          cloud(:,:,:) =    &
               (Cld_spec%crndlw(:,:,:)*Cldrad_props%emrndlw(:,:,:,5) + &
                Cld_spec%cmxolw(:,:,:)*Cldrad_props%emmxolw(:,:,:,5))/ &
               (Cld_spec%crndlw(:,:,:) + Cld_spec%cmxolw(:,:,:) +      &
                                                               1.0E-10)
          used = send_data (id_em_cld_10u, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define a total cloud field emissivity that is the weighted average
!    of the random and max overlap emissivities, over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
      else
        if (id_em_cld_lw > 0   ) then
          cloud(:,:,:) =      &
              (Cld_spec%crndlw(:,:,:)*Cldrad_props%emrndlw(:,:,:,1) +  &
               Cld_spec%cmxolw(:,:,:)*Cldrad_props%emmxolw(:,:,:,1))/ &
              (Cld_spec%crndlw(:,:,:) + Cld_spec%cmxolw(:,:,:) +     &
                                                                1.0E-10)
          used = send_data (id_em_cld_lw, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    a large scale cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_abs_lsc_cld_10u > 0) then
          used = send_data (id_abs_lsc_cld_10u,    &
                            Lscrad_props%abscoeff(:,:,:,5), Time_diag, &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the large scale cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
      else
        if (id_abs_lsc_cld_lw > 0) then
          used = send_data (id_abs_lsc_cld_lw,      &
                            Lscrad_props%abscoeff(:,:,:,1), Time_diag, &
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    the cell scale cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_abs_cell_cld_10u > 0) then
          used = send_data (id_abs_cell_cld_10u,     &
                            Cellrad_props%abscoeff(:,:,:,5), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
      else

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the cell scale cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
        if (id_abs_cell_cld_lw > 0) then
          used = send_data (id_abs_cell_cld_lw,     &
                            Cellrad_props%abscoeff(:,:,:,1), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    a meso-scale cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_abs_meso_cld_10u > 0) then
          used = send_data (id_abs_meso_cld_10u,    &
                            Mesorad_props%abscoeff(:,:,:,5), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
      else
 
!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the meso-scale cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
        if ( id_abs_meso_cld_lw > 0  ) then
          used = send_data (id_abs_meso_cld_lw,    &
                            Mesorad_props%abscoeff(:,:,:,1), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
      endif

!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is active, define 
!    the total-cloud absorption coefficient over the 990-1070 cm-1 
!    band (band 5 of 7).
!---------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        if (id_abs_cld_10u > 0) then
          used = send_data (id_abs_cld_10u,      &
                            Cldrad_props%abscoeff(:,:,:,5), Time_diag,&
                            is, js, 1, rmask=mask)
        endif
!---------------------------------------------------------------------
!    if a multi-band lw cloud emissivity formulation is not active, 
!    define the total-cloud absorption coefficient over 1 band 
!    (0-2200 cm-1).
!---------------------------------------------------------------------
      else
        if (id_abs_cld_lw > 0) then
          used = send_data (id_abs_cld_lw,    &
                            Cldrad_props%abscoeff(:,:,:,1), Time_diag, &
                            is, js, 1, rmask=mask)
        endif
      endif


!---------------------------------------------------------------------
!    define total-cloud diagnostics that are associated with the micro-
!    physically-based cloud shortwave radiative properties.
!---------------------------------------------------------------------
      if (Cldrad_control%do_sw_micro) then
        if (id_ext_cld_uv > 0) then
          cloud(:,:,:) =  Cldrad_props%cldext(:,:,:,22)
          used = send_data (id_ext_cld_uv, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif
        if (id_sct_cld_uv > 0) then
          cloud(:,:,:) =  Cldrad_props%cldsct(:,:,:,22)
          used = send_data (id_sct_cld_uv, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif
        if (id_asymm_cld_uv > 0) then
          cloud(:,:,:) = 100.0*Cldrad_props%cldasymm(:,:,:,22)
          used = send_data (id_asymm_cld_uv, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_ext_cld_vis > 0) then
          cloud(:,:,:) =  Cldrad_props%cldext(:,:,:,12)
          used = send_data (id_ext_cld_vis, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_sct_cld_vis > 0) then
          cloud(:,:,:) =  Cldrad_props%cldsct(:,:,:,12)
          used = send_data (id_sct_cld_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_asymm_cld_vis > 0) then
          cloud(:,:,:) = 100.0*Cldrad_props%cldasymm(:,:,:,12)
          used = send_data (id_asymm_cld_vis, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
        if (id_ext_cld_nir > 0) then
          cloud(:,:,:) =  Cldrad_props%cldext(:,:,:,8)
          used = send_data (id_ext_cld_nir, cloud, Time_diag,     &
                            is, js, 1, rmask=mask)
        endif
        if (id_sct_cld_nir > 0) then
          cloud(:,:,:) =  Cldrad_props%cldsct(:,:,:,8)
          used = send_data (id_sct_cld_nir, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif
        if (id_asymm_cld_nir > 0) then
          cloud(:,:,:) = 100.0*Cldrad_props%cldasymm(:,:,:,8)
          used = send_data (id_asymm_cld_nir, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    define total-cloud diagnostics that are associated with the bulk 
!    cloud shortwave radiative properties.
!---------------------------------------------------------------------
      else

!---------------------------------------------------------------------
!    define the reflected ultra-violet.
!---------------------------------------------------------------------
        if (id_alb_uv_cld > 0) then
          cloud(:,:,:) = Cldrad_props%cvisrfsw(:,:,:)
          used = send_data (id_alb_uv_cld, cloud, Time_diag,    &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    define the reflected infra-red. 
!---------------------------------------------------------------------
        if (id_alb_nir_cld > 0) then
          cloud(:,:,:) =  Cldrad_props%cirrfsw(:,:,:)
          used = send_data (id_alb_nir_cld, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif

!---------------------------------------------------------------------
!    define the absorbed  ultra-violet (not implemented).
!---------------------------------------------------------------------
!       if ( id_abs_uv_cld > 0 ) then
!         cloud = 0.0
!         used = send_data (id_abs_uv_cld, cloud, Time_diag,    &
!                           is, js, 1, rmask=mask)
!       endif

!---------------------------------------------------------------------
!    define the absorbed  infra-red.
!---------------------------------------------------------------------
        if (id_abs_nir_cld > 0) then
          cloud(:,:,:) =  Cldrad_props%cirabsw(:,:,:)
          used = send_data (id_abs_nir_cld, cloud, Time_diag,   &
                            is, js, 1, rmask=mask)
        endif
      endif 

!------------------------------------------------------------------



end subroutine cloudrad_netcdf


!####################################################################
! <SUBROUTINE NAME="cloudrad_diagnostics_end">
!  <OVERVIEW>
!    cloudrad_diagnostics_end is the destructor for 
!    cloudrad_diagnostics_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloudrad_diagnostics_end is the destructor for 
!    cloudrad_diagnostics_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloudrad_diagnostics_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine cloudrad_diagnostics_end

!-------------------------------------------------------------------
!    cloudrad_diagnostics_end is the destructor for 
!    cloudrad_diagnostics_mod.
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    be sure module has been initialized.
!--------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('cloudrad_diagnostics_mod',  &
         'initialization routine of this module was never called', &
                                                                 FATAL)
      endif

!--------------------------------------------------------------------
!    close out the component modules.
!--------------------------------------------------------------------
      if (Environment%running_gcm .or.    &
          Environment%running_sa_model) then
        if (Cldrad_control%do_strat_clouds) then
          call isccp_clouds_end
        endif
      endif

!--------------------------------------------------------------------
!    mark the module as not initialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------



end subroutine cloudrad_diagnostics_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################
! <SUBROUTINE NAME="diag_field_init">
!  <OVERVIEW>
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call diag_field_init (axes, Time)
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
!
subroutine diag_field_init (Time, axes )

!---------------------------------------------------------------------
!    diag_field_init registers the potential netcdf output variables
!    with diag_manager_mod.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time
integer        , intent(in) :: axes(4)

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       Time      initialization time for the netcdf output fields
!       axes      diagnostic variable axes
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    register the total-cloud diagnostic fields in this module.
!---------------------------------------------------------------------
      id_tot_cld_amt = register_diag_field    &
                       (mod_name, 'tot_cld_amt', axes(1:2), Time, &
                        'total cloud amount', 'percent'            )

      id_high_cld_amt = register_diag_field   &
                        (mod_name, 'high_cld_amt', axes(1:2), Time, &
                         'high cloud amount', 'percent'            )

      id_mid_cld_amt =  register_diag_field     &
                        (mod_name, 'mid_cld_amt', axes(1:2), Time, &
                         'mid cloud amount', 'percent'            )
  
      id_low_cld_amt = register_diag_field    &
                       (mod_name, 'low_cld_amt', axes(1:2), Time, &
                        'low cloud amount', 'percent'            )

      id_cld_amt =  register_diag_field     &
                    (mod_name, 'cld_amt', axes(1:3), Time,      &
                     'cloud amount', 'percent',     &
                     missing_value=missing_value            )

      id_em_cld_lw =  register_diag_field    &
                      (mod_name, 'em_cld_lw', axes(1:3), Time, &
                       'lw cloud emissivity', 'percent',        &
                       missing_value=missing_value          )
 
      id_em_cld_10u = register_diag_field    &
                      (mod_name, 'em_cld_10u', axes(1:3), Time, &
                       'cloud emissivity 10 um band', 'percent',    &
                       missing_value=missing_value          )

      id_abs_cld_lw = register_diag_field    &
                      (mod_name, 'abs_lw', axes(1:3), Time, &
                       'cloud abs coeff lw', 'percent',        &
                       missing_value=missing_value          )

     id_abs_cld_10u = register_diag_field     &
                      (mod_name, 'abs_10u', axes(1:3), Time, &
                       'cloud abs coeff 10um band', 'percent',    &
                       missing_value=missing_value          )

!---------------------------------------------------------------------
!    register diagnostic fields associated with the bulk shortwave
!    parameterization.
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_sw_micro) then
        id_alb_uv_cld = register_diag_field     &
                        (mod_name, 'alb_uv_cld', axes(1:3), Time, &
                         'UV reflected by cloud', 'percent',       &
                         missing_value=missing_value              )

        id_alb_nir_cld = register_diag_field      &
                         (mod_name, 'alb_nir_cld', axes(1:3), Time, &
                          'IR reflected by cloud', 'percent',        &
                          missing_value=missing_value               )

!   --- do not output this field ---
!       id_abs_uv_cld =  register_diag_field    &
!                        (mod_name, 'abs_uv_cld', axes(1:3), Time, &
!                         'UV absorbed by cloud', 'percent',        &
!                         missing_value=missing_value              )

        id_abs_nir_cld = register_diag_field     &
                         (mod_name, 'abs_nir_cld', axes(1:3), Time, &
                          'IR absorbed by cloud', 'percent',         &
                          missing_value=missing_value               )

!---------------------------------------------------------------------
!    register the microphysically-based total-cloud diagnostic fields.
!---------------------------------------------------------------------
      else 
        id_ext_cld_uv = register_diag_field       &
                        (mod_name, 'ext_cld_uv', axes(1:3), Time, &
                         '.27um cloud extinction coeff', 'km-1',  &
                         missing_value=missing_value          )

        id_sct_cld_uv = register_diag_field     &
                        (mod_name, 'sct_cld_uv', axes(1:3), Time, &
                         '.27um cloud scattering coeff', 'km-1', &
                         missing_value=missing_value          )

        id_asymm_cld_uv = register_diag_field    &
                          (mod_name, 'asymm_cld_uv', axes(1:3), Time, &
                           '.27um cloud asymmetry parameter',   &
                           'percent', missing_value=missing_value  ) 

        id_ext_cld_vis =  register_diag_field     &
                          (mod_name, 'ext_cld_vis', axes(1:3), Time, &
                           '.55um cloud extinction coeff', 'km-1', &
                           missing_value=missing_value          )

        id_sct_cld_vis = register_diag_field    &
                         (mod_name, 'sct_cld_vis', axes(1:3), Time, &
                          '.55um cloud scattering coeff', 'km-1', &
                          missing_value=missing_value          )

        id_asymm_cld_vis = register_diag_field      &
                           (mod_name, 'asymm_cld_vis', axes(1:3), Time,&
                            '.55um cloud asymmetry parameter',   &
                            'percent', missing_value=missing_value  )
        id_ext_cld_nir = register_diag_field    &
                         (mod_name, 'ext_cld_nir', axes(1:3), Time, &
                          '1.4um cloud extinction coeff', 'km-1', &
                          missing_value=missing_value          )

        id_sct_cld_nir = register_diag_field    &
                         (mod_name, 'sct_cld_nir', axes(1:3), Time, &
                          '1.4um cloud scattering coeff', 'km-1', &
                          missing_value=missing_value          )
 
        id_asymm_cld_nir = register_diag_field   &
                           (mod_name, 'asymm_cld_nir', axes(1:3), Time,&
                            '1.4um cloud asymmetry parameter',   &
                            'percent', missing_value=missing_value )

!---------------------------------------------------------------------
!    register the microphysically-based large-scale cloud diagnostic 
!    fields.
!---------------------------------------------------------------------
        id_lsc_cld_ext_uv = register_diag_field    &
                            (mod_name, 'lsc_cld_ext_uv', axes(1:3),   &
                             Time, '.27um lsc cloud ext coeff', 'km-1',&
                             missing_value=missing_value            )

        id_lsc_cld_ext_vis = register_diag_field   &
                             (mod_name, 'lsc_cld_ext_vis', axes(1:3), &
                              Time, '.55um lsc cloud ext coeff',   &
                              'km-1', missing_value=missing_value  )
 
        id_lsc_cld_ext_nir = register_diag_field    &
                             (mod_name, 'lsc_cld_ext_nir', axes(1:3),  &
                              Time, '1.4um lsc cloud ext coeff',   &
                              'km-1', missing_value=missing_value    )

        id_lsc_cld_sct_uv = register_diag_field    &
                            (mod_name, 'lsc_cld_sct_uv', axes(1:3),  &
                             Time, '.27um lsc cloud sct coeff', 'km-1',&
                             missing_value=missing_value            )

        id_lsc_cld_sct_vis = register_diag_field    &
                             (mod_name, 'lsc_cld_sct_vis', axes(1:3),  &
                              Time, '.55um lsc cloud sct coeff',  &
                              'km-1', missing_value=missing_value  )

        id_lsc_cld_sct_nir = register_diag_field    &
                             (mod_name, 'lsc_cld_sct_nir', axes(1:3), &
                              Time, '1.4um lsc cloud sct coeff',  &
                              'km-1', missing_value=missing_value  )
 
        id_lsc_cld_asymm_uv = register_diag_field   &
                              (mod_name, 'lsc_cld_asymm_uv', axes(1:3),&
                               Time, '.27um lsc cloud asymm coeff',  &
                               'percent', missing_value=missing_value  )

        id_lsc_cld_asymm_vis = register_diag_field  &
                               (mod_name, 'lsc_cld_asymm_vis',  &
                                axes(1:3), Time,    &
                                '.55um lsc cloud asymm coeff',   &
                                'percent', missing_value=missing_value)

        id_lsc_cld_asymm_nir = register_diag_field   &
                               (mod_name, 'lsc_cld_sct_nir', axes(1:3),&
                                Time, '1.4um lsc cloud sct coeff', &
                                'percent', missing_value=missing_value)

      endif

!---------------------------------------------------------------------
!    register the microphysically-based large-scale cloud amount and
!    lw properties diagnostic fields.
!---------------------------------------------------------------------
      id_lsc_cld_amt = register_diag_field    &
                       (mod_name, 'lsc_cld_amt', axes(1:3), Time, &
                        'lsc cloud amount', 'percent',             &
                        missing_value=missing_value            )

      id_abs_lsc_cld_lw = register_diag_field    &
                          (mod_name, 'lsc_abs_lw', axes(1:3), Time, &
                           'lsc cloud abs coeff lw', 'percent',     &
                           missing_value=missing_value          )

      id_abs_lsc_cld_10u = register_diag_field   &
                           (mod_name, 'lsc_abs_10u', axes(1:3), Time,&
                            'lsc cloud abs coeff 10um band',   &
                            'percent', missing_value=missing_value   )

!---------------------------------------------------------------------
!    register the cell-scale cloud diagnostic fields.
!---------------------------------------------------------------------
      if (Cldrad_control%do_donner_deep_clouds) then
        id_cell_cld_amt = register_diag_field    &
                          (mod_name, 'cell_cld_amt', axes(1:3), Time, &
                           'cell cloud amount', 'percent',             &
                           missing_value=missing_value            )

        id_cell_cld_ext_uv = register_diag_field   &
                             (mod_name, 'cell_cld_ext_uv', axes(1:3), &
                              Time, '.27um cell cloud ext coeff',  &
                              'km-1', missing_value=missing_value  )

        id_cell_cld_ext_vis = register_diag_field   &
                              (mod_name, 'cell_cld_ext_vis', axes(1:3),&
                               Time, '.55um cell cloud ext coeff',  &
                               'km-1', missing_value=missing_value  )

        id_cell_cld_ext_nir = register_diag_field   &
                              (mod_name, 'cell_cld_ext_nir', axes(1:3),&
                               Time, '1.4um cell cloud ext coeff',  &
                               'km-1', missing_value=missing_value  )

        id_cell_cld_sct_uv = register_diag_field    &
                             (mod_name, 'cell_cld_sct_uv', axes(1:3),  &
                              Time, '.27um cell cloud sct coeff',   &
                              'km-1', missing_value=missing_value    )

        id_cell_cld_sct_vis = register_diag_field    &
                              (mod_name, 'cell_cld_sct_vis', axes(1:3),&
                               Time, '.55um cell cloud sct coeff',  &
                               'km-1', missing_value=missing_value   )

        id_cell_cld_sct_nir = register_diag_field    &
                              (mod_name, 'cell_cld_sct_nir', axes(1:3),&
                               Time, '1.4um cell cloud sct coeff', &
                               'km-1', missing_value=missing_value )

        id_cell_cld_asymm_uv = register_diag_field    &
                               (mod_name, 'cell_cld_asymm_uv',   &
                                axes(1:3), Time,     &
                                '.27um cell cloud asymm coeff',   &
                                'percent', missing_value=missing_value )

        id_cell_cld_asymm_vis = register_diag_field     &
                                (mod_name, 'cell_cld_asymm_vis',   &
                                 axes(1:3), Time,     &
                                 '.55um cell cloud asymm coeff',   &
                                 'percent', missing_value=missing_value)

        id_cell_cld_asymm_nir = register_diag_field    &
                                (mod_name, 'cell_cld_sct_nir',    &
                                 axes(1:3), Time,&
                                 '1.4um cell cloud sct coeff',    &
                                 'percent', missing_value=missing_value)

        id_abs_cell_cld_lw = register_diag_field    &
                             (mod_name, 'cell_abs_lw', axes(1:3), Time,&
                              'cell cloud abs coeff lw', 'percent',   &
                              missing_value=missing_value          )

        id_abs_cell_cld_10u = register_diag_field    &
                              (mod_name, 'cell_abs_10u', axes(1:3), &
                               Time, 'cell cloud abs coeff 10um band', &
                               'percent',  missing_value=missing_value )

!---------------------------------------------------------------------
!    register the meso-scale cloud diagnostic fields.
!---------------------------------------------------------------------
        id_meso_cld_amt = register_diag_field     &
                          (mod_name, 'meso_cld_amt', axes(1:3), Time, &
                           'meso cloud amount', 'percent',             &
                           missing_value=missing_value            )

        id_meso_cld_ext_uv = register_diag_field    &
                             (mod_name, 'meso_cld_ext_uv', axes(1:3), &
                              Time, '.27um meso cloud ext coeff',   &
                              'km-1', missing_value=missing_value)

        id_meso_cld_ext_vis = register_diag_field   &
                              (mod_name, 'meso_cld_ext_vis', axes(1:3),&
                               Time, '.55um meso cloud ext coeff',  &
                               'km-1', missing_value=missing_value )

        id_meso_cld_ext_nir = register_diag_field   &
                              (mod_name, 'meso_cld_ext_nir', axes(1:3),&
                               Time, '1.4um meso cloud ext coeff',  &
                               'km-1', missing_value=missing_value    )

        id_meso_cld_sct_uv = register_diag_field   &
                             (mod_name, 'meso_cld_sct_uv', axes(1:3), &
                              Time, '.27um meso cloud sct coeff',  &
                              'km-1', missing_value=missing_value )

        id_meso_cld_sct_vis = register_diag_field  &
                              (mod_name, 'meso_cld_sct_vis', axes(1:3),&
                               Time, '.55um meso cloud sct coeff',  &
                               'km-1', missing_value=missing_value   )

        id_meso_cld_sct_nir = register_diag_field  &
                              (mod_name, 'meso_cld_sct_nir', axes(1:3),&
                               Time, '1.4um meso cloud sct coeff',  &
                               'km-1', missing_value=missing_value )

        id_meso_cld_asymm_uv = register_diag_field  &
                               (mod_name, 'meso_cld_asymm_uv',   &
                                axes(1:3), Time,     &
                                '.27um meso cloud asymm coeff',   &
                                'percent', missing_value=missing_value)

        id_meso_cld_asymm_vis = register_diag_field   &
                                (mod_name, 'meso_cld_asymm_vis',  &
                                 axes(1:3), Time,      &
                                 '.55um meso cloud asymm coeff',    &
                                 'percent', missing_value=missing_value)

        id_meso_cld_asymm_nir = register_diag_field    &
                                (mod_name, 'meso_cld_sct_nir',   &
                                 axes(1:3), Time,&
                                 '1.4um meso cloud sct coeff',   &
                                 'percent', missing_value=missing_value)

        id_abs_meso_cld_lw = register_diag_field    &
                             (mod_name, 'meso_abs_lw', axes(1:3), Time,&
                              'meso cloud abs coeff lw', 'percent',   &
                              missing_value=missing_value          )

        id_abs_meso_cld_10u = register_diag_field   &
                              (mod_name, 'meso_abs_10u', axes(1:3),  &
                               Time, 'meso cloud abs coeff 10um band', &
                               'percent', missing_value=missing_value)
      endif


!---------------------------------------------------------------------

 
end subroutine diag_field_init



!#####################################################################
! <SUBROUTINE NAME="isccp_diag">
!  <OVERVIEW>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_diag (is, js, Cld_spec, Atmos_input, coszen, Time)
!  </TEMPLATE>
!  <IN NAME="is,js" TYPE="integer">
!   starting subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                        diagnostic output [ time_type (days, seconds) ]
!  </IN>
!  <IN NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </IN>
!  <IN NAME="coszen" TYPE="real">
!    cosine of solar zenith angle
!  </IN>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </IN>
! </SUBROUTINE>
!
subroutine isccp_diag (is, js, Cld_spec, Atmos_input, coszen, Time)

!--------------------------------------------------------------------
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!---------------------------------------------------------------------
 
integer,                      intent(in)   :: is,js
type(cld_specification_type), intent(in)   :: Cld_spec
type(atmos_input_type),       intent(in)   :: Atmos_input
real, dimension(:,:),         intent(in)   :: coszen
type(time_type),              intent(in)   :: Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting/ending subdomain i,j indices of data 
!                      in the physics_window being integrated
!      Cld_spec        cloud specification properties on model grid,
!                      [ cld_specification_type ]
!      Atmos_input     atmospheric input fields on model grid,
!                      [ atmos_input_type ] 
!      coszen          cosine of zenith angle [ dimensionless ]
!      Time            time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      real, dimension (size(Cld_spec%lwp,3)) :: &
                                    tau_local, em_local, cldamt_local

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                       size(Cld_spec%lwp,3) ) ::  qv, em_lw_local, rh2o

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                       size(Cld_spec%lwp,3)+1 ) ::  temp   

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                       size(Cld_spec%lwp,3), 4 ) ::  tau, tau_ice

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2), &
                      7, 7) ::       fq_isccp

      real, dimension (size(Cld_spec%lwp,1), size(Cld_spec%lwp,2)) :: &
                                     npoints

      integer      :: kdim
      integer      :: max_cld
      integer      :: i, j, k

!---------------------------------------------------------------------
!   local variables:
!
!      tau_local        optical depth in band 1 in the current column
!                       [ dimensionless ]
!      em_local         lw cloud emissivity in current column
!                       [ dimensionless ]
!      cldamt_local     cloud fraction in current column 
!                       [ dimensionless ]
!      qv               water vapor specific humidity
!                       [ kg vapor / kg air ]
!      em_lw_local      lw cloud emissivity [ dimensionless ]
!      rh2o             mixing ratio of water vapor at model full levels
!                       [ non-dimensional ]
!      temp             temperature at model levels (1:nlev), surface
!                       temperature is stored at value nlev+1; if eta
!                       coordinates, surface value stored in below 
!                       ground points
!                       [ deg K ]
!      tau              optical depth in 4 bands [ dimensionless ]
!      tau_ice          optical depth due to cloud ice (4 bands)
!                       [ dimensionless ]
!      fq_isccp         matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges
!      npoints          flag indicating whether isccp cloud is present
!                       in column (cloud + daylight needed)
!      kdim             number of model layers
!      max_cld          greatest number of clouds in any column in the
!                       current physics window
!      i,j,k            do-loop indices
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    define number of model layers.
!----------------------------------------------------------------------
      kdim = size (Cld_spec%lwp,3)

!---------------------------------------------------------------------
!    if any netcdf isccp diagnostics have been requested through the
!    diag_table, continue with this routine. if not, fall through if 
!    loop and return.
!---------------------------------------------------------------------
      if (do_isccp .or. do_tau_reff) then
        
!---------------------------------------------------------------------
!    determine if there are any clouds in this physics window.
!---------------------------------------------------------------------
        max_cld = MAXVAL(Cld_spec%ncldsw(:,:))
         
!---------------------------------------------------------------------
!    if clouds exist in the window, call cloud_optical_properties_diag
!    to define the cloud optical depth, the optical depth due to cloud
!    ice and the longwave emissivity.
!---------------------------------------------------------------------
        if (max_cld >= 1) then
          call cloud_optical_properties_diag (Cld_spec, tau, tau_ice, &
                                              em_lw_local)

!----------------------------------------------------------------------
!    if no clouds exist in the window, set all the optical depths and 
!    emissivities to zero.
!----------------------------------------------------------------------
        else
           em_lw_local = 0.
           tau = 0.
           tau_ice = 0.
        end if  

!---------------------------------------------------------------------
!    if any isccp diagnostics are desired,  define the needed input
!    fields that determine the isccp category into which different
!    clouds will fall.
!---------------------------------------------------------------------
        if (do_isccp) then
          do j=1,size(Cld_spec%lwp,2)
            do i=1,size(Cld_spec%lwp,1)

!---------------------------------------------------------------------
!    isccp clouds are only defined in sunlight.
!---------------------------------------------------------------------
              if (coszen(i,j) > 1.E-06) then
                do k=1,kdim                      

!--------------------------------------------------------------------
!    define the specific humidity from the mixing ratio which has been
!    input.
!--------------------------------------------------------------------
                  qv(i,j,k) = Atmos_input%cloudvapor(i,j,k)/   &
                              (1. + Atmos_input%cloudvapor(i,j,k))

!---------------------------------------------------------------------
!    define the column values of cloud fraction, cloud optical depth, 
!    and lw cloud emissivity. if cloud is not present, set these var-
!    iables to clear sky values.
!---------------------------------------------------------------------
                  if (Cld_spec%camtsw(i,j,k) > 0.0) then
                    cldamt_local(k) = Cld_spec%camtsw(i,j,k)  
                    tau_local(k) = tau(i,j,k,1)/ &
                                 real(Cld_spec%cld_thickness(i,j,k))
                    em_local(k) = 1. - ( (1.-em_lw_local(i,j,k))** &
                              (1./real(Cld_spec%cld_thickness(i,j,k))) )
                  else
                    cldamt_local(k) = 0.
                    tau_local(k) = 0.
                    em_local(k) = 0.
                  endif
                end do

!---------------------------------------------------------------------
!    call isccp_cloudtypes to map each model cloud to an isccp cloud
!    type, based on its optical depth and height above the surface.
!    set a flag to indicate the presence of isccp cloud in this column.
!---------------------------------------------------------------------
                call isccp_cloudtypes (Atmos_input%press(i,j,1:kdim), &
                                       Atmos_input%pflux(i,j,:),&
                                       qv(i,j,:),       &
                                       Atmos_input%cloudtemp(i,j,:),  &
                                       Atmos_input%temp(i,j,kdim+1),  &
                                       cldamt_local, tau_local,   &
                                       em_local, fq_isccp(i,j,:,:))
                npoints(i,j) = 1.

!----------------------------------------------------------------------
!    if it is not daylight, set the isccp clouds to zero.
!----------------------------------------------------------------------
              else
                npoints(i,j) = 0.
                fq_isccp(i,j,:,:) = 0.
              end if
            end do
          end do
         
!----------------------------------------------------------------------
!    send any desired diagnostics to the diag_manager_mod.
!----------------------------------------------------------------------
          call isccp_output (is, js, fq_isccp, npoints, Time)
        end if   !(do_isccp)

!---------------------------------------------------------------------
!    if any isccp summary diagnostics are desired, call tau_reff_diag2
!    to process them.
!---------------------------------------------------------------------
        if (do_tau_reff) then
          call tau_reff_diag2 (is, js, Time, coszen, Cld_spec%ncldsw, &
                               Cld_spec%cld_thickness,   &
                               Cld_spec%camtsw,  &
                               Atmos_input%cloudtemp(:,:,1:kdim),   &
                               Atmos_input%pflux, tau(:,:,:,1),  &
                               tau_ice(:,:,:,1), Cld_spec%reff_liq,  &
                               Cld_spec%reff_ice)
        end if   
      end if   !(do_isccp .or. do_tau_reff)
   
!---------------------------------------------------------------------
    
    
end subroutine isccp_diag        



!#####################################################################
! <SUBROUTINE NAME="compute_isccp_clds">
!  <OVERVIEW>
!    subroutine compute_isccp_clds maps the model clouds into isccp
!    categories (high, middle, low) and defines the cloud fraction of
!    each.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine compute_isccp_clds maps the model clouds into isccp
!    categories (high, middle, low) and defines the cloud fraction of
!    each.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call compute_isccp_clds (pflux, camtsw, hml_ca)
!  </TEMPLATE>
!  <IN NAME="pflux" TYPE="real">
!   average of pressure at adjacent model levels
!  </IN>
!  <IN NAME="camtsw" TYPE="real">
!   total cloud amount [ nondimensional ]
!  </IN>
!  <OUT NAME="hml_ca" TYPE="real">
!   cloud fraction for the 3 isccp cloud types
!  </OUT>
! </SUBROUTINE>
!
subroutine compute_isccp_clds (pflux, camtsw, hml_ca)

!---------------------------------------------------------------------
!    subroutine compute_isccp_clds maps the model clouds into isccp
!    categories (high, middle, low) and defines the cloud fraction of
!    each.
!--------------------------------------------------------------------- 

real,  dimension(:,:,:),   intent(in)  :: pflux, camtsw
real,  dimension(:,:,:),   intent(out) :: hml_ca

!---------------------------------------------------------------------
!  intent(in) variables:
!
!        pflux           average of pressure at adjacent model levels
!                        [ (kg /( m s^2) ]
!        camtsw          total cloud amount [ nondimensional ]
!
!  intent(out) variable:
!
!        hml_ca          cloud fraction for the 3 isccp cloud types
!                        [ nondimensional ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real,  parameter ::   mid_btm  = 6.8e4
      real, parameter  ::   high_btm = 4.4e4
      integer          ::   i, j, k

!---------------------------------------------------------------------
!  local variables:
!
!         mid_btm     pressure boundary between isccp middle and 
!                     isccp low clouds [ Pa ]
!         high_btm    pressure boundary between isccp middle and
!                     isccp high clouds [ Pa ]
!         i,j,k       do-loop indices
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    initialize a column integrated hi-mid-low cloud-free area array.
!--------------------------------------------------------------------
      hml_ca = 1.0
 
!---------------------------------------------------------------------
!    define arrays giving the cloud-free area in each column within the
!    pressure regionscorresponding to the ISCCP definitions of high 
!    (10-440 hPa), middle (440-680 hPa) and low (680-1000 hPa) clouds. 
!    compute high, middle and low cloud amounts assuming that independ-
!    ent clouds overlap randomly. note that model clouds above 10 hPa 
!    and below 1000 hPa are included in the totals.
!----------------------------------------------------------------------
      do k = 1,size(pflux,3)-1
        do j=1,size(pflux,2)
          do i=1,size(pflux,1)
            if (pflux(i,j,k)  <=  high_btm) then
              hml_ca(i,j,1) = hml_ca(i,j,1) * (1. - camtsw(i,j,k))
            else if ( (pflux(i,j,k) >  high_btm) .and.  &
                      (pflux(i,j,k) <=  mid_btm) ) then
              hml_ca(i,j,2) = hml_ca(i,j,2) * (1. - camtsw(i,j,k))
            else  if ( pflux(i,j,k) > mid_btm ) then
              hml_ca(i,j,3) = hml_ca(i,j,3) * (1. - camtsw(i,j,k))
            endif
          end do
        end do
      end do

!--------------------------------------------------------------------
!    convert the cloud-free area to an integrated cloud fraction in 
!    the column. express the cloud area in percent.
!--------------------------------------------------------------------
      hml_ca = 1. - hml_ca
      hml_ca = 100.*hml_ca
  
!-------------------------------------------------------------------


end subroutine compute_isccp_clds



!####################################################################
! <SUBROUTINE NAME="cloud_optical_properties_diag">
!  <OVERVIEW>
!    cloud_optical_properties_diag calculates the cloud optical depth,
!    ice cloud optical depth and longwave cloud emissivity for each
!    cloudy grid box.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloud_optical_properties_diag calculates the cloud optical depth,
!    ice cloud optical depth and longwave cloud emissivity for each
!    cloudy grid box.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_optical_properties_diag (Cld_spec, tau, tau_ice, em_lw)
!  </TEMPLATE>
!  <IN NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid
!  </IN>
!  <OUT NAME="tau" TYPE="real">
!   cloud optical depth in each of the
!                     num_slingo_bands
!  </OUT>
!  <OUT NAME="tau_ice" TYPE="real">
!   ice cloud optical depth in each  of the 
!                     num_slingo_bands
!  </OUT>
!  <OUT NAME="em_lw" TYPE="real">
!   longwave cloud emissivity
!  </OUT>
! </SUBROUTINE>
!
subroutine cloud_optical_properties_diag (Cld_spec, tau, tau_ice, em_lw)

!---------------------------------------------------------------------
!    cloud_optical_properties_diag calculates the cloud optical depth,
!    ice cloud optical depth and longwave cloud emissivity for each
!    cloudy grid box.
!---------------------------------------------------------------------
                              
type(cld_specification_type), intent(in)   :: Cld_spec
real, dimension(:,:,:,:),     intent(out)  :: tau, tau_ice
real, dimension(:,:,:),       intent(out)  :: em_lw       

!--------------------------------------------------------------------
!   intent(in) variable:
!
!      Cld_spec       cloud specification properties on model grid,
!                     [ cld_specification_type ]
!
!   intent(out) variables:
!
!      tau            cloud optical depth in each of the
!                     num_slingo_bands [ dimensionless ]
!      tau_ice        ice cloud optical depth in each  of the 
!                     num_slingo_bands [ dimensionless ]
!      em_lw          longwave cloud emissivity [ dimensionless ]  
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:


      real, dimension (size(Cld_spec%lwp,1),size(Cld_spec%lwp,2),  &
                       size(Cld_spec%lwp,3), 4) ::     &
                                                tau_liq

      real, dimension (size(Cld_spec%lwp,1),size(Cld_spec%lwp,2),   &
                       size(Cld_spec%lwp,3)) ::     &
                                                k_liq, k_ice

!--------------------------------------------------------------------
!   local variables:
!
!       tau_liq    liquid cloud optical depth [ dimensionless ]
!       k_liq      liquid cloud mass absorption coefficient for longwave
!                  portion of the spectrum [ m**2 / kg condensate ]
!       k_ice      ice cloud mass absorption coefficient for longwave
!                  portion of the spectrum [ m**2 / kg condensate ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    compute uv cloud optical depths due to liquid. the formula for 
!    optical depth comes from: 
!    Slingo (1989), J. Atmos. Sci., vol. 46, pp. 1419-1427
!---------------------------------------------------------------------
      tau_liq(:,:,:,1) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.02817 + (1.305/Cld_spec%reff_liq(:,:,:)))
      tau_liq(:,:,:,2) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.02682 + (1.346/Cld_spec%reff_liq(:,:,:)))
      tau_liq(:,:,:,3) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.02264 + (1.454/Cld_spec%reff_liq(:,:,:)))
      tau_liq(:,:,:,4) = Cld_spec%lwp(:,:,:) * 1000. * &
                         (0.01281 + (1.641/Cld_spec%reff_liq(:,:,:)))
        
!---------------------------------------------------------------------
!    compute uv cloud optical depths due to ice. the formula for 
!    optical depth comes from:
!    Ebert and Curry (1992), J. Geophys. Res., vol. 97, pp. 3831-3836
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    IMPORTANT:  WE ARE CHEATING HERE BECAUSE WE ARE FORCING THE FIVE 
!                BAND MODEL OF EBERT AND CURRY INTO THE FOUR BAND MODEL 
!                OF SLINGO. THIS IS DONE BY COMBINING BANDS 3 and 4 OF 
!                EBERT AND CURRY. EVEN SO THE EXACT BAND LIMITS DO NOT 
!                MATCH.  FOR COMPLETENESS HERE ARE THE BAND LIMITS IN 
!                MICRONS:
!            BAND               SLINGO                 EBERT AND CURRY
!
!             1               0.25-0.69                0.25 - 0.7
!             2               0.69-1.19                0.7 - 1.3
!             3               1.19-2.38                1.3 - 2.5
!             4               2.38-4.00                2.5 - 3.5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!---------------------------------------------------------------------
      tau_ice(:,:,:,1) = Cld_spec%iwp(:,:,:) * 1000. * &
                         (0.003448 + (2.431/Cld_spec%reff_ice(:,:,:)))
      tau_ice(:,:,:,2) = tau_ice(:,:,:,1)
      tau_ice(:,:,:,3) = tau_ice(:,:,:,1)
      tau_ice(:,:,:,4) = tau_ice(:,:,:,1)
        
!---------------------------------------------------------------------
!    compute total cloud optical depth. the mixed phase optical prop-
!    erties are based upon equation 14 of Rockel et al. 1991, 
!    Contributions to Atmospheric Physics, volume 64, pp.1-12. 
!    thus:

!          tau = tau_liq + tau_ice
!
!---------------------------------------------------------------------
      tau(:,:,:,:) = tau_liq(:,:,:,:) + tau_ice(:,:,:,:)
        
!----------------------------------------------------------------------
!    place a minimum value on tau.
!----------------------------------------------------------------------
      where (tau(:,:,:,:) .lt. taumin)
        tau(:,:,:,:) = taumin
      end where   
        
!----------------------------------------------------------------------
!    define the  mass absorption coefficient for longwave radiation 
!    for cloud ice and cloud liquid.
!----------------------------------------------------------------------
      k_liq(:,:,:) = 140.
      k_ice(:,:,:) = 4.83591 + 1758.511/Cld_spec%reff_ice(:,:,:)
        
!----------------------------------------------------------------------
!    compute combined lw emmisivity. the mixed phase optical properties
!    are based upon equation 14 of Rockel et al. 1991, Contributions to
!    Atmospheric Physics,  volume 64, pp.1-12.  thus:
!
!    transmivvity_lw =   transmissivity_lw_ice * transmissivity_lw_liq
!
!    which can also be written as:
!
!    em_lw =  em_lw_liq + em_lw_ice -  (em_lw_liq * em_lw_ice )
!
!    which is what is solved here.
!----------------------------------------------------------------------
      em_lw(:,:,:) =  1. - exp(-1.*(k_liq(:,:,:)*Cld_spec%lwp(:,:,:) + &
                                    k_ice(:,:,:)*Cld_spec%iwp(:,:,:)))

!---------------------------------------------------------------------


end subroutine cloud_optical_properties_diag


!######################################################################



                end module cloudrad_diagnostics_mod

