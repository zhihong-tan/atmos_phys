                 module isccp_clouds_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    isccp_clouds partitions the model cloud fields into the isccp
!    cloud categories, by cld top height and cld amount, and provides
!    diagnostic and netcdf output.
! </OVERVIEW>
! <DESCRIPTION>
!    isccp_clouds partitions the model cloud fields into the isccp
!    cloud categories, by cld top height and cld amount, and provides
!    diagnostic and netcdf output.
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

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    isccp_clouds partitions the model cloud fields into the isccp
!    cloud categories, by cld top height and cld amount, and provides
!    diagnostic and netcdf output.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: isccp_clouds.F90,v 10.0 2003/10/24 22:00:41 fms Exp $'
character(len=128)  :: tagname =  '$Name: jakarta $'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
         isccp_clouds_init, &
         isccp_output, tau_reff_diag2, isccp_cloudtypes, &
         isccp_clouds_end

private          &
!   called from isccp_clouds_init:
         diag_field_init, &
!   called from isccp_cloudtypes:
         ran0


!---------------------------------------------------------------------
!-------- namelist  ---------
!
!                         ISCCP CLOUD PROCESSING
!    
!       The following variables only are used if the ISCCP cloud view
!       processing is done.  ISCCP cloud processing is activated by
!       requesting one of the diagnostic fields produced by it from
!       the diag_table. 
!
!       top_height
!
!                      integer variable indicating whether 
!                      or not to simulate 10.5 micron brightness
!                      temperatures to adjust top heights according
!                      to the emissivity of the cloud. 
!                     
!                      1 = adjust top height using both a computed
!                          infrared brightness temperature and the
!                          visible optical depth to adjust cloud top
!                          pressure. Note that this calculation is
!                          most appropriate to compare to ISCCP data
!                          during sunlit hours.
!       
!                      2 = do not adjust top height, that is cloud top
!                          pressure is the actual cloud top pressure
!                          in the model
!       
!                      3 = adjust top height using only the computed
!                          infrared brightness temperature. Note that
!                          this calculation is most appropriate to
!                          compare to ISCCP IR only algortihm (i.e.
!                          you can compare to nighttime ISCCP data
!                          with this option)
!                      
!       ncol           number of columns used in ISCCP cloud type
!                      simulations
! 
!       isccp_taumin   minimum optical depth ISCCP can see
!                      [ dimensionless ]
!
!       emsfclw        assumed constant fraction of longwave emissivity
!                      of the surface [ dimensionless ]
! 
!       do_sunlit      should ISCCP diagnostics be done during sunlit
!                      hours only?
!       
!----------------------------------------------------------------------

integer  ::  ncol = 50
integer  ::  top_height = 1
real     ::  isccp_taumin = 0.3
real     ::  emsfclw = 0.94
logical  ::  do_sunlit = .false.

namelist /isccp_clouds_nml /                             &
                                      ncol,                      &
                                      top_height,   &
                                      isccp_taumin, emsfclw,    &
                                      do_sunlit

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
character(len=8)    :: mod_name = 'isccp_clouds'
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
! <SUBROUTINE NAME="isccp_clouds_init">
!  <OVERVIEW>
!    isccp_clouds_init is the constructor for isccp_clouds_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    isccp_clouds_init is the constructor for isccp_clouds_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_clouds_init (axes, Time, do_isccp_out,    &
!                              do_tau_reff_out, do_sunlit_out,   &
!                              overlap_in, qmin_in)
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
!  <OUT NAME="do_isccp_out" TYPE="logical">
!   OPTIONAL: processing isccp diagnostics output from model
!  </OUT>
!  <OUT NAME="do_tau_reff_out" TYPE="logical">
!   OPTIONAL: processing isccp summary diagnostics output from model
!  </OUT>
!  <OUT NAME="do_sunlit_out" TYPE="logical">
!   OPTIONAL: processing isccp diagnosics only for sunlit hours
!  </OUT>
!  <IN NAME="overlap_in" TYPE="integer">
!    if strat_cloud_mod is active, verify that that module has been 
!    activated and retrieve the overlap parameter and the value used 
!    for the minimum cloud water amount. Cloud overlap parameter.
!  </IN>
!  <IN NAME="qmin_in" TYPE="real">
!    if strat_cloud_mod is active, verify that that module has been 
!    activated and retrieve the overlap parameter and the value used 
!    for the minimum cloud water amount. Minimum cloud water amount.
!  </IN> 
! </SUBROUTINE>
!
subroutine isccp_clouds_init (axes, Time, do_isccp_out,    &
                              do_tau_reff_out, do_sunlit_out,   &
                              overlap_in, qmin_in)

!---------------------------------------------------------------------
!    isccp_clouds_init is the constructor for isccp_clouds_mod.
!------------------------------------------------------------------

integer, dimension(4),   intent(in)              :: axes
type(time_type),         intent(in)              :: Time
integer,                 intent(in),  optional   :: overlap_in
real,                    intent(in),  optional   :: qmin_in
logical,                 intent(out), optional   :: do_isccp_out,   &
                                                    do_tau_reff_out, &
                                                    do_sunlit_out

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       axes             diagnostic variable axes
!       Time             current time [time_type(days, seconds)]
!
!   intent(in),optional variables:
!
!       overlap_in       type of cloud overlap assumption
!       qmin_in          minimum value of cloud condensate
!
!   intent(out),optional variables:
!
!       do_isccp_out     isccp diagnostics are desired ?
!       do_tau_reff_out  isccp summary diagnostics are desired ?
!       do_sunlit_out    isccp diagnostics only consider daytime clouds?
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
      if (module_is_initialized)   then
        if (present(do_isccp_out)) do_isccp_out = do_isccp
        if (present(do_sunlit_out)) do_sunlit_out = do_sunlit
        if (present(do_tau_reff_out)) do_tau_reff_out = do_tau_reff
        return
      endif

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call diag_manager_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=isccp_clouds_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'isccp_clouds_nml')
        enddo
10      call close_file (unit)
      endif
 
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() )    &
                       write (stdlog(), nml=isccp_clouds_nml)
 

!-------------------------------------------------------------------
!    initialize the netcdf diagnostics provided with this module.
!-------------------------------------------------------------------
      call diag_field_init (Time, axes)

!-------------------------------------------------------------------
!    provide output variables if they are desired.
!-------------------------------------------------------------------
      if (present(do_isccp_out))    do_isccp_out = do_isccp
      if (present(do_tau_reff_out)) do_tau_reff_out = do_tau_reff
      if (present(do_sunlit_out))   do_sunlit_out = do_sunlit

!---------------------------------------------------------------------
!    if strat_cloud_mod is active, verify that that module has been 
!    activated and retrieve the overlap parameter and the value used 
!    for the minimum cloud water amount.
!---------------------------------------------------------------------
      if (present(overlap_in)) overlap = overlap_in
      if (present(qmin_in))   qmin = qmin_in

!--------------------------------------------------------------------
!    mark the module initialized.
!--------------------------------------------------------------------
      module_is_initialized= .true.

!--------------------------------------------------------------------



end subroutine isccp_clouds_init



!######################################################################
! <SUBROUTINE NAME="isccp_output">
!  <OVERVIEW>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </OVERVIEW>
!  <DESCRIPTION>
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_output (is, js, fq_isccp, npoints, Time)
!  </TEMPLATE>
!  <IN NAME="is, js" TYPE="integer">
!	starting/ending subdomain i,j indices of data 
!                      in the physics_window being integrated
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                      diagnostic output 
!  </IN>
!  <IN NAME="fq_isccp" TYPE="real">
!   matrix of fractional area covered by cloud
!                      types of a given optical depth and cloud
!                      top pressure range.  The matrix is 7x7 for
!                      7 cloud optical depths and 7 cloud top 
!                      pressure ranges
!  </IN>
!  <IN NAME="npoints" TYPE="real">
!   flag indicating whether isccp cloud is present
!                      in column (cloud + daylight needed)
!  </IN>
! </SUBROUTINE>
!
subroutine isccp_output (is, js, fq_isccp, npoints, Time)

!--------------------------------------------------------------------
!    subroutine isccp_diag maps the model cloud distribution to the
!    isccp cloud categories, and provides netcdf output if desired.
!---------------------------------------------------------------------
 
integer,                      intent(in)   :: is,js
real, dimension(:,:,:,:),     intent(in)   :: fq_isccp
real, dimension(:,:),         intent(in)   :: npoints 
type(time_type),              intent(in)   :: Time

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting/ending subdomain i,j indices of data 
!                      in the physics_window being integrated
!      fq_isccp        matrix of fractional area covered by cloud
!                      types of a given optical depth and cloud
!                      top pressure range.  The matrix is 7x7 for
!                      7 cloud optical depths and 7 cloud top 
!                      pressure ranges
!      npoints         flag indicating whether isccp cloud is present
!                      in column (cloud + daylight needed)
!      Time            time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!
!---------------------------------------------------------------------

!local variable:

     logical :: used    !  flag returned from send_data indicating
                        !  whether diag_manager_mod has received 
                        !  data that was sent

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('isccp_clouds_mod',   &
               'module has not been initialized', FATAL )
      endif

!----------------------------------------------------------------------
!    send any desired diagnostics to the diag_manager_mod.
!----------------------------------------------------------------------
      used = send_data (id_pc1tau1, fq_isccp(:,:,1,1), Time, &
                        is, js )
      used = send_data (id_pc1tau2, fq_isccp(:,:,2,1), Time, &
                        is, js )
      used = send_data (id_pc1tau3, fq_isccp(:,:,3,1), Time, &
                        is, js )
      used = send_data (id_pc1tau4, fq_isccp(:,:,4,1), Time, &
                        is, js )
      used = send_data (id_pc1tau5, fq_isccp(:,:,5,1), Time, &
                        is, js )
      used = send_data (id_pc1tau6, fq_isccp(:,:,6,1), Time, &
                        is, js )
      used = send_data (id_pc1tau7, fq_isccp(:,:,7,1), Time, &
                        is, js )
      used = send_data (id_pc2tau1, fq_isccp(:,:,1,2), Time, &
                        is, js )
      used = send_data (id_pc2tau2, fq_isccp(:,:,2,2), Time, &
                        is, js )
      used = send_data (id_pc2tau3, fq_isccp(:,:,3,2), Time, &
                        is, js )
      used = send_data (id_pc2tau4, fq_isccp(:,:,4,2), Time, &
                        is, js )
      used = send_data (id_pc2tau5, fq_isccp(:,:,5,2), Time, &
                        is, js )
      used = send_data (id_pc2tau6, fq_isccp(:,:,6,2), Time, &
                        is, js )
      used = send_data (id_pc2tau7, fq_isccp(:,:,7,2), Time, &
                        is, js )
      used = send_data (id_pc3tau1, fq_isccp(:,:,1,3), Time, &
                        is, js )
      used = send_data (id_pc3tau2, fq_isccp(:,:,2,3), Time, &
                        is, js )
      used = send_data (id_pc3tau3, fq_isccp(:,:,3,3), Time, &
                        is, js )
      used = send_data (id_pc3tau4, fq_isccp(:,:,4,3), Time, &
                        is, js )
      used = send_data (id_pc3tau5, fq_isccp(:,:,5,3), Time, &
                        is, js )
      used = send_data (id_pc3tau6, fq_isccp(:,:,6,3), Time, &
                        is, js )
      used = send_data (id_pc3tau7, fq_isccp(:,:,7,3), Time, &
                        is, js )
      used = send_data (id_pc4tau1, fq_isccp(:,:,1,4), Time, &
                        is, js )
      used = send_data (id_pc4tau2, fq_isccp(:,:,2,4), Time, &
                        is, js )
      used = send_data (id_pc4tau3, fq_isccp(:,:,3,4), Time, &
                        is, js )
      used = send_data (id_pc4tau4, fq_isccp(:,:,4,4), Time, &
                        is, js )
      used = send_data (id_pc4tau5, fq_isccp(:,:,5,4), Time, &
                        is, js )
      used = send_data (id_pc4tau6, fq_isccp(:,:,6,4), Time, &
                        is, js )
      used = send_data (id_pc4tau7, fq_isccp(:,:,7,4), Time, &
                        is, js )
      used = send_data (id_pc5tau1, fq_isccp(:,:,1,5), Time, &
                        is, js )
      used = send_data (id_pc5tau2, fq_isccp(:,:,2,5), Time, &
                        is, js )
      used = send_data (id_pc5tau3, fq_isccp(:,:,3,5), Time, &
                        is, js )
      used = send_data (id_pc5tau4, fq_isccp(:,:,4,5), Time, &
                        is, js )
      used = send_data (id_pc5tau5, fq_isccp(:,:,5,5), Time, &
                        is, js )
      used = send_data (id_pc5tau6, fq_isccp(:,:,6,5), Time, &
                        is, js )
      used = send_data (id_pc5tau7, fq_isccp(:,:,7,5), Time, &
                        is, js )
      used = send_data (id_pc6tau1, fq_isccp(:,:,1,6), Time, &
                        is, js )
      used = send_data (id_pc6tau2, fq_isccp(:,:,2,6), Time, &
                        is, js )
      used = send_data (id_pc6tau3, fq_isccp(:,:,3,6), Time, &
                        is, js )
      used = send_data (id_pc6tau4, fq_isccp(:,:,4,6), Time, &
                        is, js )
      used = send_data (id_pc6tau5, fq_isccp(:,:,5,6), Time, &
                        is, js )
      used = send_data (id_pc6tau6, fq_isccp(:,:,6,6), Time, &
                        is, js )
      used = send_data (id_pc6tau7, fq_isccp(:,:,7,6), Time, &
                        is, js )
      used = send_data (id_pc7tau1, fq_isccp(:,:,1,7), Time, &
                        is, js )
      used = send_data (id_pc7tau2, fq_isccp(:,:,2,7), Time, &
                        is, js )
      used = send_data (id_pc7tau3, fq_isccp(:,:,3,7), Time, &
                        is, js )
      used = send_data (id_pc7tau4, fq_isccp(:,:,4,7), Time, &
                        is, js )
      used = send_data (id_pc7tau5, fq_isccp(:,:,5,7), Time, &
                        is, js )
      used = send_data (id_pc7tau6, fq_isccp(:,:,6,7), Time, &
                        is, js )
      used = send_data (id_pc7tau7, fq_isccp(:,:,7,7), Time, &
                        is, js )
      used = send_data (id_nisccp, npoints, Time, is, js )

!--------------------------------------------------------------------
         

end subroutine isccp_output



!#######################################################################
! <SUBROUTINE NAME="tau_reff_diag2">
!  <OVERVIEW>
!    tau_reff_diag2 computes the typical cloud optical depth, effective
!    radius and temperature for each column of the atmosphere, and for 
!    liquid and ice clouds separately. the method is to select the 
!    compacted cloud with the greatest horizontal area coverage and use 
!    the properties of that cloud in determining  ?????. the routine 
!    also determines low liquid cloud and low ice cloud optical depths
!    and low cloud cloudtop and mean temperatures. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    tau_reff_diag2 computes the typical cloud optical depth, effective
!    radius and temperature for each column of the atmosphere, and for 
!    liquid and ice clouds separately. the method is to select the 
!    compacted cloud with the greatest horizontal area coverage and use 
!    the properties of that cloud in determining  ?????. the routine 
!    also determines low liquid cloud and low ice cloud optical depths
!    and low cloud cloudtop and mean temperatures. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call tau_reff_diag2 (is, js, Time, coszen, nclds, cld_thickness, &
!                           cldamt, tkel, phalf, tau, tau_ice,  &
!                           reff_liq, reff_ice)
!  </TEMPLATE>
!  <IN NAME="is, js" TYPE="integer">
!	starting/ending subdomain i,j indices of data 
!                      in the physics_window being integrated
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   time on next timestep, used as stamp for 
!                      diagnostic output 
!  </IN>
!  <IN NAME="coszen" TYPE="real">
!   Cosine of zenith angle --  mean value over
!                      appropriate averaging interval
!  </IN>
!  <IN NAME="nclds" TYPE="integer">
!   number of clouds in the column
!  </IN>
!  <IN NAME="cld_thickness" TYPE="real">
!    thickness of cloud at this grid point
!  </IN>
!  <IN NAME="cldamt" TYPE="real">
!   cloud fraction [ dimensionless ]
!  </IN>
!  <IN NAME="tkel" TYPE="real">
!   temperature at model levels
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!   pressure on model half levels
!  </IN>
!  <IN NAME="tau" TYPE="real">
!   visible cloud optical depth
!  </IN>
!  <IN NAME="tau_ice" TYPE="real">
!   visible ice cloud optical depth
!  </IN>
!  <IN NAME="reff_liq" TYPE="real">
!   liquid cloud effective radius
!  </IN>
!  <IN NAME="reff_ice" TYPE="real">
!   ice cloud effective radius
!  </IN>
! </SUBROUTINE>
!
subroutine tau_reff_diag2 (is, js, Time, coszen, nclds, cld_thickness, &
                           cldamt, tkel, phalf, tau, tau_ice,  &
                           reff_liq, reff_ice)

!----------------------------------------------------------------------
!    tau_reff_diag2 computes the typical cloud optical depth, effective
!    radius and temperature for each column of the atmosphere, and for 
!    liquid and ice clouds separately. the method is to select the 
!    compacted cloud with the greatest horizontal area coverage and use 
!    the properties of that cloud in determining  ?????. the routine 
!    also determines low liquid cloud and low ice cloud optical depths
!    and low cloud cloudtop and mean temperatures. 
!----------------------------------------------------------------------
 
integer,                   intent(in)    :: is, js
type(time_type),           intent(in)    :: Time
real, dimension(:,:),      intent(in)    :: coszen
integer, dimension(:,:),   intent(in)    :: nclds   
integer, dimension(:,:,:), intent(in)    :: cld_thickness
real, dimension(:,:,:),    intent(in)    :: cldamt, tkel, phalf, tau, &
                                            tau_ice, reff_liq, reff_ice

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      is,js           starting subdomain i,j indices of data 
!                      in the physics_window being integrated
!      Time            time on next timestep, used as stamp for 
!                      diagnostic output [ time_type (days, seconds) ]
!      coszen          cosine of zenith angle --  mean value over
!                      appropriate averaging interval
!                      [ non-dimensional ]
!      nclds           number of clouds in the column
!      cld_thickness   thickness of cloud at this grid point
!                      [ number of model layers ]
!      cldamt          cloud fraction [ dimensionless ]
!      tkel            temperature at model levels [ deg K ]
!      phalf           pressure on model half levels [ Pa ]
!      tau             visible cloud optical depth [ dimensionless ]
!      tau_ice         visible ice cloud optical depth [ dimensionless ]
!      reff_liq        liquid cloud effective radius [ microns ]
!      reff_ice        ice cloud effective radius [ microns ]
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real, dimension (size(tkel,1), size(tkel,2),    &
                       size(tkel,3))   :: tau_liq

      real, dimension (size(tkel,1), size(tkel,2))   ::    &
                                          tau_cld_liq, tau_cld_ice, &
                                          reff_cld_liq, reff_cld_ice, &
                                          t_cld_low, t_cld_lay, &
                                          cldamt_space_liq, &
                                          cldamt_space_ice, &
                                          tca_space, tca_space_last,&
                                          aliq_space, aice_space,&
                                          alow_space, tlay, sum_logp
      real     :: tca_space_loc, tca_space_loc_high, tca_space_loc_last
      logical  :: used
      integer  :: kdim 
      integer  :: maxclds, kcld
      integer  :: i,j,k

!---------------------------------------------------------------------
!   local variables:
!
!      tau_liq                   liquid cloud optical depth
!                                [ dimensionless ]
!      tau_cld_liq               liquid cloud optical depth of low
!                                clouds as seen from space 
!                                [ dimensionless ]
!      tau_cld_ice               ice cloud optical depth of low clouds 
!                                as seen from space [ dimensionless ]
!      reff_cld_liq              effective drop size of liquid clouds
!                                seen from space [ microns ]
!      reff_cld_ice              effective drop size of ice clouds
!                                seen from space [ microns ]
!      t_cld_low                 low cloud cloudtop temperature as
!                                seen from space [ 
!      t_cld_lay                 log pressure weighted mean atmospheric
!                                temperature below 680 hPa [ deg K ]
!      cldamt_space_liq          liquid cloud fraction added to total
!                                fraction seen from space from current 
!                                level
!                                [ dimensionless ]
!      cldamt_space_ice          ice cloud fraction added to total
!                                fraction seen from space from current 
!                                level
!                                [ dimensionless ]
!      tca_space                 total cloud fraction seen from space
!                                in region extending to current level
!                                [ dimensionless ]
!      tca_space_last            total cloud fraction seen from space
!                                in region extending to previous level
!                                [ dimensionless ]
!      aliq_space                area covered by water clouds as seen
!                                from space [ dimensionless ]
!      aice_space                area covered by ice clouds as seen 
!                                from space [ dimensionless ]
!      alow_space                area covered by low clouds as seen
!                                from space [ dimensionless ]
!      tlay                      log-pressure weighted temperature of
!                                low cloud layers [ deg K ]
!      sum_logp                  sum of log-pressure over the low cloud
!                                region 
!      tca_space_loc             variable used to accumulate total cloud
!                                area as seen from space as one marches
!                                down a column [ dimensionless ]
!      tca_space_loc_high        variable used to accumulate high cloud
!                                area as seen from space as one marches
!                                down a column [ dimensionless ]
!      tca_space_loc_last        variable used to accumulate total cloud
!                                area as seen from space from model top
!                                to the previous layer[ dimensionless ]
!      used                      flag returned from send_data indicating
!                                whether diag_manager_mod has received 
!                                data that was sent
!      kdim                      number of model layers
!      maxclds                   maximum number of clouds in any column
!                                in current physics window
!      kcld                      next k index at which to look for 
!                                another cloud
!      i,j,k                     do-loop indices
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('isccp_clouds_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    define the number of model layers.
!---------------------------------------------------------------------
      kdim = size (tkel,3)

!---------------------------------------------------------------------
!    determine if any clouds are present in the current window.
!---------------------------------------------------------------------
      maxclds = MAXVAL(nclds(:,:))
    
!--------------------------------------------------------------------
!    if any cloud is present, define the desired diagnostics.
!--------------------------------------------------------------------
      if (maxclds > 0) then    
    
!-------------------------------------------------------------------
!    compute liquid cloud optical depth as the difference between total
!    cloud optical depth and ice cloud optical depth.
!-------------------------------------------------------------------
        tau_liq = MAX (tau - tau_ice, 0.)

!---------------------------------------------------------------------
!    determine the ice and liquid cloud characteristics as seen from
!    space in each cloud column. initialize the area and effective
!    particle sizes to 0.0. 
!---------------------------------------------------------------------
        do j=1,size(tkel,2)
          do i=1,size(tkel,1)
            reff_cld_liq(i,j) = 0.
            reff_cld_ice(i,j) = 0.
            aliq_space(i,j) = 0.
            aice_space(i,j) = 0.

!--------------------------------------------------------------------
!    if it is daytime, calculate the total cloud area in a column as 
!    seen from space. this field is set to zero at model top.
!--------------------------------------------------------------------
            if (coszen(i,j) > 1.E-06) then        
              tca_space(i,j) = 0.

!---------------------------------------------------------------------
!    define kcld to be the next model level at which to check for the 
!    presence of cloud. this is used so thaty multi-layer clouds are
!    only counted once.
!---------------------------------------------------------------------
              kcld = 1
              do k=1,kdim       
                if (k >= kcld) then
                  if (cldamt(i,j,k) > 0.0) then

!---------------------------------------------------------------------
!    upon encountering the next cloud, update the total cloud amount 
!    seen from space.
!---------------------------------------------------------------------
                    tca_space_last(i,j) = tca_space(i,j)
                    tca_space(i,j) = 1. - ((1. - tca_space(i,j))*   &
                                     (1. - cldamt(i,j,k)))
!---------------------------------------------------------------------
!    if this is a liquid cloud, update the liquid cloud amount seen
!    space. otherwise set the contribution to the liquid cloud amount 
!    seen from space from this level to be zero, since the cloud is 
!    an ice cloud.
!---------------------------------------------------------------------
                    if ((tau_liq(i,j,k) > 0.) .and.   &
                        (tau_ice(i,j,k) <= 0.)) then
                      cldamt_space_liq(i,j) =  MIN(   &
                        MAX(0., tca_space(i,j) - tca_space_last(i,j)), &
                                    1.)
                    else
                      cldamt_space_liq(i,j) = 0.
                    endif

!---------------------------------------------------------------------
!    if this is an ice cloud, update the ice cloud amount seen from 
!    space. otherwise set the contribution to the ice cloud amount seen
!    from space from this level to be zero, since the cloud is a 
!    water cloud.
!---------------------------------------------------------------------
                    if ( (tau_ice(i,j,k) > 0.)  .and.   &
                         (tau_liq(i,j,k) <= 0.) ) then
                      cldamt_space_ice(i,j) = MIN(   &
                        MAX(0., tca_space(i,j) - tca_space_last(i,j)),&
                                    1.)
                    else
                      cldamt_space_ice(i,j) = 0.
                    endif

!---------------------------------------------------------------------
!    add the liquid and ice areas and area-weighted particle sizes from
!    the current level to the accumulated sum. note that particle sizes
!    are multiplied by the area of the cloud seen from space so 
!    that division by the area will yield the local cloud properties.
!    skip over any remaining layers in the current cloud by incrementing
!    kcld by the cloud thickness.
!---------------------------------------------------------------------
                    aliq_space(i,j) = aliq_space(i,j) +    &
                                      cldamt_space_liq(i,j)
                    reff_cld_liq(i,j) = reff_cld_liq(i,j) + &
                                        cldamt_space_liq(i,j)*  &
                                        reff_liq(i,j,k) 
                    aice_space(i,j) = aice_space(i,j) +     &
                                      cldamt_space_ice(i,j)
                    reff_cld_ice(i,j) = reff_cld_ice(i,j) + &
                                        cldamt_space_ice(i,j)*  &
                                        reff_ice(i,j,k) 
                    kcld = kcld + cld_thickness(i,j,k)

!--------------------------------------------------------------------
!    if cloud is not present on the current level, increment kcld by 1
!    so that the next level may be checked.
!--------------------------------------------------------------------
                  else  ! (cldamt > 0)
                    kcld = kcld + 1 
                  endif ! (cldamt > 0)
                endif ! (k >= kcld)
              end do  !loop over levels
            endif ! (coszen)
          end do
        end do

!---------------------------------------------------------------------
!    compute low cloud temperature as the log pressure weighted
!    mean temperature of layers where the pressure is greater than
!    680 mb.  This choice conforms to isccp definitions of low cloud.
!---------------------------------------------------------------------
        sum_logp(:,:) = 0.
        tlay(:,:) = 0.
        do k=1,kdim
          where (((phalf(:,:,k+1) + phalf(:,:,k))/2.) >= 68000.)
               sum_logp(:,:) = sum_logp(:,:) + &
                             (log(phalf(:,:,k+1)) - log(phalf(:,:,k)))
               tlay(:,:) = tlay(:,:) + TKel(:,:,k)* &
                         (log(phalf(:,:,k+1)) - log(phalf(:,:,k)))
          end where
        end do
        where (sum_logp(:,:) > 0.)
          tlay(:,:) = tlay(:,:)/sum_logp(:,:)
        elsewhere 
          tlay(:,:) = 0.
        end where
    
!---------------------------------------------------------------------
!    loop over columns to determine the low liquid cloud and low ice 
!    cloud optical depths and the low cloud cloudtop temperature in a
!    each column. 
!---------------------------------------------------------------------
        do j=1,size(tkel,2)
          do i=1,size(tkel,1)

!--------------------------------------------------------------------
!    initialize the arrays to hold low liquid cloud optical depth, low
!    ice cloud optical depth, low cloud cloudtop temperature,  mean 
!    temperature of the low clouds and the clow cloud area as seen
!    from space.
!---------------------------------------------------------------------
            tau_cld_liq(i,j)  = 0.
            tau_cld_ice(i,j)  = 0.
            T_cld_low(i,j)  = 0.
            T_cld_lay(i,j)  = 0.
            alow_space(i,j) = 0.

!---------------------------------------------------------------------
!    set the total cloud amount and total non-low cloud amount as seen
!    from space to 0.0 at model top.
!---------------------------------------------------------------------
            if (coszen(i,j) > 1.E-06) then        
              tca_space_loc = 0.
              tca_space_loc_high = 0.

!---------------------------------------------------------------------
!    define kcld to be the next model level at which to check for the 
!    presence of cloud. this is used so that multi-layer clouds are
!    only counted once.
!---------------------------------------------------------------------
              kcld = 1

!--------------------------------------------------------------------
!    march downwards, looking for clouds.
!--------------------------------------------------------------------
              do k=1,kdim        

!---------------------------------------------------------------------
!    if k .lt. kcld, then skip this level. if cloud is present it is
!    part of a multi-layer cloud which has been previously processed.
!---------------------------------------------------------------------
                if (k >= kcld) then

!---------------------------------------------------------------------
!    if cloud is present, add the contribution of this cloud layer to
!    the total cloud amount as seen from space. if there is no cloud in
!    this layer, skip to the next layer.
!---------------------------------------------------------------------
                  if (cldamt(i,j,k) > 0.0) then
                    tca_space_loc_last = tca_space_loc
                    tca_space_loc = 1. - ((1. - tca_space_loc)*   &
                                    (1. - cldamt(i,j,k)))

!--------------------------------------------------------------------
!    the total non-low cloud amount as seen from space is the sum of
!    all contributions of clouds above 680 hPa. 
!--------------------------------------------------------------------
                    if (((phalf(i,j,k+1) +     &
                          phalf(i,j,k))/2.) < 68000.) then
                      tca_space_loc_high = tca_space_loc
                    endif
        
!--------------------------------------------------------------------
!    generate the low ice and low liquid cloud optical depths as the
!    area-weighted sum of the optical depths of each appropriate low
!    cloud. the low cloud top temperature is also area-weighted.
!--------------------------------------------------------------------
                    if (((phalf(i,j,k+1) + &
                          phalf(i,j,k))/2.) >= 68000. ) then      
                      tau_cld_ice(i,j) = tau_cld_ice(i,j) + &
                                         (cldamt(i,j,k)*    &
                                         (1. - tca_space_loc_high)*&
                                         tau_ice(i,j,k))
                      tau_cld_liq(i,j) = tau_cld_liq(i,j) + &
                                         (cldamt(i,j,k)*    &
                                         (1. - tca_space_loc_high)*&
                                         tau_liq(i,j,k))
                      T_cld_low(i,j) = T_cld_low(i,j) + &
                                       tKel(i,j,k)* &
                                       (tca_space_loc -    &
                                                    tca_space_loc_last)
                    endif

!---------------------------------------------------------------------
!    set kcld to point to the next level to check for cloudiness. if
!    multi-layer clouds are present, the remaining layers in that cloud 
!    are skipped over. 
!----------------------------------------------------------------------
                    kcld = kcld + cld_thickness(i,j,k)
                  else
            kcld = kcld + 1 
                  end if
                end if
              end do ! (loop over k)
  
!---------------------------------------------------------------------
!    define low cloud area as the difference between the total and the
!    high cloud area. define the mean low cloud temperature as the 
!    low cloud area times the log pressure-weighted low-level 
!    temperature.
!---------------------------------------------------------------------
              alow_space(i,j) = tca_space_loc - tca_space_loc_high
              T_cld_lay(i,j) = (tca_space_loc - tca_space_loc_high)*   &
                               Tlay(i,j)
           endif  ! (coszen > 1.0e-6)
          end do  !loop over i
        end do  !loop over j

!--------------------------------------------------------------------
!    if there are no clouds in the window, zero out the diagnostics
!    arrays.
!--------------------------------------------------------------------
      else
        tau_cld_liq(:,:)  = 0.
        tau_cld_ice(:,:)  = 0.
        reff_cld_liq(:,:) = 0.
        reff_cld_ice(:,:) = 0.
        T_cld_lay(:,:)  = 0.
        T_cld_low(:,:)  = 0.
        aliq_space(:,:) = 0.
        aice_space(:,:) = 0.
        alow_space(:,:) = 0.
      end if  
    
!---------------------------------------------------------------------
!    send any desired netcdf diagnostics to diag_manager_mod.
!---------------------------------------------------------------------
      if (id_aice > 0) then
        used = send_data ( id_aice, aice_space, Time, is, js )
      end if
      if (id_reffice > 0) then
        used = send_data ( id_reffice, reff_cld_ice, Time, is, js )
      end if
      if (id_aliq > 0) then
        used = send_data ( id_aliq, aliq_space, Time, is, js )
      end if
      if (id_reffliq > 0) then
        used = send_data ( id_reffliq, reff_cld_liq, Time, is, js )
      end if
      if (id_alow > 0) then
        used = send_data ( id_alow, alow_space, Time, is, js )
      end if
      if (id_tauicelow > 0) then
        used = send_data ( id_tauicelow, tau_cld_ice, Time, is, js )
      end if
      if (id_tauliqlow > 0) then
        used = send_data ( id_tauliqlow, tau_cld_liq, Time, is, js )
      end if
      if (id_tlaylow > 0) then
        used = send_data ( id_tlaylow, T_cld_lay, Time, is, js )
      end if
      if (id_tcldlow > 0) then
        used = send_data ( id_tcldlow, T_cld_low, Time, is, js )
      end if

!---------------------------------------------------------------------



end subroutine tau_reff_diag2



!######################################################################
! <SUBROUTINE NAME="isccp_cloudtypes">
!  <OVERVIEW>
!    isccp_cloudtypes calculates the fraction of each model grid box 
!    covered by each of the 49 ISCCP D level cloud types 
!    (i.e. stratified by optical depth and cloud top pressure) by 
!    accounting for model overlap. 
!  </OVERVIEW>
!  <DESCRIPTION>
!    isccp_cloudtypes calculates the fraction of each model grid box 
!    covered by each of the 49 ISCCP D level cloud types 
!    (i.e. stratified by optical depth and cloud top pressure) by 
!    accounting for model overlap. For further explanation see Klein 
!    and Jakob, Monthly Weather Review, (2000), vol x, pp. .
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_cloudtypes (pfull, phalf, qv, at, skt, cc, dtau_s,  &
!                             dem_s, fq_isccp)
!  </TEMPLATE>
!  <IN NAME="pfull" TYPE="real">
!   pressure of full model levels, pfull(1) is top
!                        level of model, pfull(nlev) is bottom level of
!                        model
!  </IN>
!  <IN NAME="phalf" TYPE="real">
!   pressure of half model levels, phalf(1) is top
!                        of model, phalf(nlev+1) is the surface pressure
!  </IN>
!  <IN NAME="qv" TYPE="real">
!   water vapor specific humidity on model levels.
!                        used only if top_height = 1 or top_height = 3.
!  </IN>
!  <IN NAME="at" TYPE="real">
!   temperature in each model level [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!  </IN>
!  <IN NAME="skt" TYPE="real">
!   skin temperature [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!  </IN>
!  <INOUT NAME="cc" TYPE="real">
!   cloud cover in each model layer [ fraction ]
!                        this includes convective clouds if any
!                        NOTE:  This is the HORIZONTAL area of each
!                               grid box covered by clouds
!  </INOUT>
!  <INOUT NAME="dtau_s" TYPE="real">
!   mean 0.67 micron optical depth of stratiform
!                        clouds in each model level [ dimensionless ]
!                        NOTE:  this the cloud optical depth of only the
!                               cloudy part of the grid box, it is not 
!                               weighted with the 0 cloud optical depth 
!                               of the clear part of the grid box
!  </INOUT>
!  <INOUT NAME="dem_s" TYPE="real">
!   10.5 micron longwave emissivity of stratiform
!                        clouds in each model level. 
!                        used only if top_height = 1 or top_height = 3.
!                        Same note applies as in dtau. [ dimensionless ]
!  </INOUT>
!  <OUT NAME="fq_isccp" TYPE="real">
!   matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges. [ fraction ]
!  </OUT>
! </SUBROUTINE>
!
subroutine isccp_cloudtypes (pfull, phalf, qv, at, skt, cc, dtau_s,  &
                             dem_s, fq_isccp)

!---------------------------------------------------------------------
!    isccp_cloudtypes calculates the fraction of each model grid box 
!    covered by each of the 49 ISCCP D level cloud types 
!    (i.e. stratified by optical depth and cloud top pressure) by 
!    accounting for model overlap. For further explanation see Klein 
!    and Jakob, Monthly Weather Review, (2000), vol x, pp. .
!
!---------------------------------------------------------------------

real, dimension(:),   intent(in)      :: pfull, phalf, qv, at
real,                 intent(in)      :: skt
real, dimension(:),   intent(inout)   :: cc, dtau_s, dem_s
real, dimension(7,7), intent(out)     :: fq_isccp

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       pfull            pressure of full model levels, pfull(1) is top
!                        level of model, pfull(nlev) is bottom level of
!                        model [ Pa ]
!       phalf            pressure of half model levels, phalf(1) is top
!                        of model, phalf(nlev+1) is the surface pressure
!                        [ Pa ]
!       qv               water vapor specific humidity on model levels.
!                        used only if top_height = 1 or top_height = 3.
!                        [ kg vapor / kg air ]
!       at               temperature in each model level [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!       skt              skin temperature [ deg K ]
!                        used only if top_height = 1 or top_height = 3.
!
!   intent(inout) variables:
!
!       cc               cloud cover in each model layer [ fraction ]
!                        this includes convective clouds if any
!                        NOTE:  This is the HORIZONTAL area of each
!                               grid box covered by clouds
!       dtau_s           mean 0.67 micron optical depth of stratiform
!                        clouds in each model level [ dimensionless ]
!                        NOTE:  this the cloud optical depth of only the
!                               cloudy part of the grid box, it is not 
!                               weighted with the 0 cloud optical depth 
!                               of the clear part of the grid box
!       dem_s            10.5 micron longwave emissivity of stratiform
!                        clouds in each model level. 
!                        used only if top_height = 1 or top_height = 3.
!                        Same note applies as in dtau. [ dimensionless ]
!
!       NOTE :  OPTION TO RUN WITH CONVECTIVE CLOUDS IS NOT
!               IMPLEMENTED YET
!
!       conv             convective cloud cover in each model 
!                        level (fraction) this includes convective 
!                        clouds if any
!  
!                        NOTE:  This is the HORIZONTAL area of each
!                               grid box covered by clouds
!                         
!       dtau_c           mean 0.67 micron optical depth of convective
!                        clouds in each model level
!
!                        NOTE:  this the cloud optical depth of only the
!                               cloudy part of the grid box, it is not 
!                               weighted with the 0 cloud optical depth 
!                               of the clear part of the grid box
!
!   intent(out) variable:
!
!       fq_isccp        matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges. [ fraction ]
!
!---------------------------------------------------------------------
    
!---------------------------------------------------------------------
!   local variables:

      real,    dimension (size(qv,1))        :: demwv
      integer, dimension (size(qv,1)-1)      :: match
      integer, dimension (ncol)              :: levmatch
      real,    dimension (ncol)              :: tau,tb,ptop
      real,    dimension (ncol,size(qv,1))   :: frac_out
      real,    dimension (ncol,0:size(qv,1)) :: tca
      real,    dimension (ncol)              :: threshold,maxosc,boxpos
      real,    dimension (ncol)              :: threshold_min
!convective cloud stuff:
!     real,    dimension (ncol,size(qv,1))  :: cca
!     real,    dimension (ncol)             :: maxocc

      real        ::      ptrop, attrop, atmax, atmin, btcmin,  &
                          transmax, fluxtop, fluxtopclrsky, trans,  &
                          transclrsky, wtmair, wtmh20, Navo, grav,  &
                          pstd, t0, press, dpress, atmden, rvh20,   &
                          wk, rhoave, rh20s, rfrgn, tmpexp, tauwv,  &
                          emcld, taumintmp, bb, dem, fluxtopinit,   &
                          tauir, boxarea
      integer     ::      ilev, j, ibox, itrop, ipres, itau, nmatch, &
                          nlev, seed, icycle

!---------------------------------------------------------------------
!   local variables:
!
!       demwv
!       match
!       levmatch
!       tau
!       tb
!       ptop
!       frac_out
!       tca
!       threshold
!       maxosc
!       boxpos
!       threshold_min
!       cca
!       maxocc
!       ptrop          pressure at tropopause [ pa ]
!       attrop         temperature at troopause [ deg K ]
!       atmax          maximum temperature in column [ deg K ]
!       atmin          minimum temperature in column [ deg K ]
!       btcmin
!       transmax
!       fluxtop
!       fluxtopclrsky
!       trans
!       transclrsky
!       wtmair
!       wtmh20
!       Navo
!       grav
!       pstd
!       t0
!       press
!       dpress
!       atmden
!       rvh20
!       wk
!       rhoave
!       rh20s
!       rfrgn
!       tmpexp
!       tauwv
!       emcld
!       taumintmp
!       bb
!       dem
!       fluxtopinit
!       tauir
!       boxarea
!       ilev
!       j
!       ibox
!       itrop          k index of tropopause
!       ipres
!       itau
!       nmatch
!       nlev           number of model layers
!       seed
!       icycle
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('isccp_clouds_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    define the number of model layers.
!---------------------------------------------------------------------
      nlev = size(qv,1) 
     
!---------------------------------------------------------------------
!    find tropopause pressure and minimum and maximum temperature in the
!    column. the tropopause pressure is used only in the assignment of
!    cloud top pressure in the case where the infrared brightness temp-
!    erature is calculated  (top_height = 1 or 3). Clouds whose infrared
!    brightness temperatures are less than the tropopause temperature 
!    are assigned to a cloud top pressure of the tropopause pressure 
!    (see ISCCP documentation)
!----------------------------------------------------------------------
      if (top_height .eq. 1 .or. top_height .eq. 3) then
        ptrop = 5000.
        atmin = 400.
        atmax = 0.
        do ilev=1,nlev-1
           if ((pfull(ilev)/phalf(nlev+1)) .lt. 0.4 .and. &
                at(ilev) .gt. at(ilev+1)) then
             ptrop = pfull(ilev+1)
             attrop = at(ilev+1)
             itrop=ilev+1
           end if
           if (at(ilev) > atmax) atmax = at(ilev)
           if (at(ilev) < atmin) atmin = at(ilev)
        end do
      end if

!---------------------------------------------------------------------
!    remove any non-physical input values of cloud fraction (must be 
!    between 0 and 1), optical depth (non-negative), emissivity (between
!    0 and 1).
!---------------------------------------------------------------------
      where (cc(:) .lt. 0.) 
             cc(:) = 0.
      end where
      where (cc(:) .gt. 1.) 
             cc(:) = 1.
      end where
      where (dtau_s(:) .lt. 0.) 
             dtau_s(:) = 0.
      end where
      where (dem_s(:) .lt. 0.) 
             dem_s(:) = 0.
      end where
      where (dem_s(:) .gt. 1.) 
             dem_s(:) = 1.
      end where
     
      !where (dtau_c(:) .lt. 0.) 
      !       dtau_c(:) = 0.
      !end where
      !where (dem_c(:) .lt. 0.) 
      !       dem_c(:) = 0.
      !end where
      !where (dem_c(:) .gt. 1.) 
      !       dem_c(:) = 1.
      !end where
     

!     ---------------------------------------------------!
!     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
!

      !initialize variables
      
      !assign 2d tca array using 1d input array cc
      !assign 2d cca array using 1d input array conv
      tca(:,0) = 0.
      do ilev = 1,nlev
             tca(:,ilev) = cc(ilev)
             !cca(:,ilev  ) = conv(ilev)
      enddo
            
      seed = (pfull(nlev)-int(pfull(nlev)))*100+1
      frac_out(:,:)=0.

      do ibox=1,ncol
           boxpos(ibox)=(real(ibox)-.5)/real(ncol)
      enddo

      !frac_out is the array that contains the information 
      !where 0 is no cloud, 1 is a stratiform cloud and 2 is a
      !convective cloud
      do ilev=1, nlev


           !Initialise threshold
           IF (ilev.eq.1) then
               DO ibox=1,ncol
            ! If max overlap 
            ! select pixels spread evenly
            ! across the gridbox
                    !threshold(ibox)=boxpos(ibox)
          
                    ! for non-maximum overlap options
            ! select random pixels from the non-convective
            ! part the gridbox ( some will be converted into
            ! convective pixels below )
                    threshold(ibox)= ran0(seed)
!                   threshold(ibox)= cca(ibox,ilev)+ &
!                                    (1-cca(ibox,ilev))*ran0(seed)            
               ENDDO
           ENDIF

           DO ibox=1,ncol

                   !All versions
                   !if (boxpos(ibox).le.cca(ibox,ilev)) then
                   !     maxocc(ibox) = 1
                   !else
                   !     maxocc(ibox) = 0
                   !end if

                   !Max overlap
                   !      threshold_min(ibox)=cca(ibox,ilev)
                   !       maxosc(ibox)=1
                   
                   !Random overlap
                   if (overlap .eq. 2) then
                   !      threshold_min(ibox)=cca(ibox,ilev)
                          threshold_min(ibox)=0.
                          maxosc(ibox)=0
                   end if

                   !Max/Random overlap
                   if (overlap .eq. 1) then
                     
                         threshold_min(ibox)= &
                               min(tca(ibox,ilev-1),tca(ibox,ilev))
                              !max(cca(ibox,ilev), &
                              !min(tca(ibox,ilev-1),tca(ibox,ilev)))
                         if (threshold(ibox).lt.  &
                               min(tca(ibox,ilev-1),tca(ibox,ilev))) then
                              !.and.(threshold(ibox).gt.cca(ibox,ilev))          
                               maxosc(ibox)=1.
                         else
                               maxosc(ibox)=0.
                         end if 

                   end if 
                         
    
                   !Reset threshold 
                   threshold(ibox)=     &                                
                           (maxosc(ibox) * threshold(ibox)) + &
                           ((1-maxosc(ibox)) * &
                            ( threshold_min(ibox)+ &
                             (1-threshold_min(ibox))*ran0(seed) ) )

                   !original code...................
                   !threshold(ibox)=     &                                
                   !if max overlapped conv cloud
                   !maxocc(ibox) * (  &                                     
                   ! boxpos(ibox)     &                                         
                   !        ) +       &                                               
                   !else
                   !(1-maxocc(ibox)) * ( &                                   
                   !if max overlapped strat cloud
                   ! (maxosc(ibox)) * (  &                                
                   !threshold=boxpos
                   !      threshold(ibox) &                                        
                   !            ) +       &                                          
                   !else
                   !    (1-maxosc(ibox)) * (   &                            
                   !threshold_min=random[thrmin,1]
                   !     threshold_min(ibox)+  &
                   !  (1-threshold_min(ibox))*ran0(seed)  
                   !           ) 
                   !    )

           ENDDO


           !Fill frac_out with 1's where tca is greater than the threshold
           WHERE(tca(:,ilev).gt.threshold(:))
                 frac_out(:,ilev)=1
           ENDWHERE

           !Code to partition boxes into startiform and convective parts
           !goes here
           !
           !DO ibox=1,ncol
           !    if (threshold(ibox).le.cca(ibox,ilev)) then
           !    ! = 2 IF threshold le cca(ibox)
           !        frac_out(ibox,ilev) = 2
           !    else
           !    ! = the same IF NOT threshold le cca(ibox) 
           !        frac_out(ibox,ilev) = frac_out(ibox,ilev)
           !    end if
           !ENDDO

      enddo      !loop over levels



!     ---------------------------------------------------!

      
!
!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau
 
      !compute total cloud optical depth for each column     
      tau(:)=0.
      do ilev=1,nlev
            tau(:)=tau(:)+ &
                 real(frac_out(:,ilev))*dtau_s(ilev)
                 !if (frac_out(ibox,ilev).eq.1) then
                 !       dtautmp(ibox)= dtau_s(ilev)
                 !else if (frac_out(ibox,ilev).eq.2) then
                 !       dtautmp(ibox)= dtau_c(ilev)
                 !else
                 !       dtautmp(ibox)= 0.
                 !end if
                 !tau(ibox)=tau(ibox)+dtautmp(ibox)
      end do

!
!     ---------------------------------------------------!



!     
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done if top_height = 1 or 3
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     fluxtopclrsky is the 10.5 micron radiance at the top of
!              the atmosphere under clear skies
!     trans is the transmissivity from the top of level to
!              the top of the atmosphere
!     transclrsky is the clear sky transmissivity from the top
!              of level to the top of the atmosphere


      if (top_height .eq. 1 .or. top_height .eq. 3) then




      !compute water vapor continuum emissivity
      !this treatment follows Schwarkzopf and Ramaswamy
      !JGR 1999,vol 104, pages 9467-9499.
      !the emissivity is calculated at a wavenumber of 955 cm-1, 
      !or 10.47 microns 
      do ilev=1,nlev
               wtmair = 28.9644
               wtmh20 = 18.01534
               Navo = 6.023E+23
               grav = 9.806650E+02
               pstd = 1.013250E+06
               t0 = 296.
               !press and dpress are dyne/cm2 = Pascals *10
               press = pfull(ilev)*10.
               dpress = (phalf(ilev+1)-phalf(ilev))*10
               !atmden = g/cm2 = kg/m2 / 10 
               atmden = dpress/grav
               rvh20 = qv(ilev)*wtmair/wtmh20
               wk = rvh20*Navo*atmden/wtmair
               rhoave = (press/pstd)*(t0/at(ilev))
               rh20s = rvh20*rhoave
               rfrgn = rhoave-rh20s
               tmpexp = exp(-0.02*(at(ilev)-t0))
               tauwv = wk*1.e-20*( (0.0224697*rh20s*tmpexp) + &
                      (3.41817e-7*rfrgn)         )*0.98
               demwv(ilev) = 1. - exp( -1. * tauwv)
               
      end do

      !loop over columns 
      do ibox=1,ncol
      
            fluxtop=0.
            fluxtopclrsky=0.
            trans=1.
            transclrsky=1.
          
            do ilev=1,nlev
 
                   !blackbody emission
                   bb=1 / ( exp(1307.27/at(ilev)) - 1. )
         
                   !compute emissivity for this layer
                   dem= 1. -  &
                           (1. - demwv(ilev)) *  &
                           (1. - frac_out(ibox,ilev)*dem_s(ilev))

                            !if (frac_out(ibox,ilev).eq.1) then
                            ! dem(ibox)= 1. -  ( (1. - demwv(ilev)) *
                            ! (1. -  dem_s(ilev)) )
                            !else if (frac_out(ibox,ilev).eq.2) then
                            !  dem(ibox)= 1. - ( (1. - demwv(ilev)) *
                            ! (1. -  dem_c(ilev)) )
                            ! else
                            !dem(ibox)=  dem_wv(ilev)
                            !end if

                   ! increase TOA flux by flux emitted from layer
                   ! times total transmittance in layers above
                   fluxtop = fluxtop   + dem * bb  * trans 
                   fluxtopclrsky = fluxtopclrsky   +  &
                           demwv(ilev) * bb  * transclrsky 
                   
                   ! update trans_layers_above with transmissivity
                   ! from this layer for next time around loop
                   trans=trans*(1.-dem)
                   transclrsky=transclrsky*(1.-demwv(ilev))

            end do

            !add in surface emission
            fluxtop = fluxtop + trans * emsfclw/(exp(1307.27/skt)-1.)
            fluxtopclrsky = fluxtopclrsky +                        &
                      (transclrsky * emsfclw/(exp(1307.27/skt)-1.))
            
            
            !now that you have the top of atmosphere radiance account
            !for ISCCP procedures to determine cloud top temperature

            !account for partially transmitting cloud recompute flux 
            !ISCCP would see assuming a single layer cloud
            !note choice here of 2.13, as it is primarily ice
            !clouds which have partial emissivity and need the 
            !adjustment performed in this section
            !
            !If it turns out that the cloud brightness temperature
            !is greater than 260K, then the liquid cloud conversion
            !factor of 2.56 is used.
            !
            !Note that this is discussed on pages 85-87 of 
            !the ISCCP D level documentation (Rossow et al. 1996)
           
            
            !compute minimum brightness temperature and optical depth
            btcmin = 1. /  ( exp(1307.27/(attrop-5.)) - 1. ) 
            transmax = (fluxtop-btcmin)/(fluxtopclrsky-btcmin)
            taumintmp = -1. * log(max(min(transmax,(1.-qmin)),0.001))
            
            !note that the initial setting of tauir is needed so that
            !tauir has a realistic value should the next if block be
            !bypassed
            tauir = tau(ibox) / 2.13
            
            if (top_height .eq. 1 .and. transmax .gt. 0.001 .and. &
                transmax .le. (1.-qmin)) then
                    icycle = 1
                    fluxtopinit = fluxtop
                    tauir = tau(ibox) / 2.13
243                 emcld = 1. - exp(-1. * tauir  )
                    fluxtop = fluxtopinit -                            &
                              ((1.-emcld)*fluxtopclrsky)
                    fluxtop = max(1.E-06,(fluxtop/emcld))
                    tb(ibox)= 1307.27/ (log(1. + (1./fluxtop)))
                    if (icycle .eq. 1 .and. tb(ibox) .gt. 260.) then
                         tauir = tau(ibox) / 2.56
                         icycle = 2
                         go to 243
                     end if
            end if

            if (tau(ibox) .gt. (-1.*log((1.-qmin)))) then 
                
                !cloudy box
                tb(ibox)= 1307.27/ (log(1. + (1./fluxtop)))
                if (top_height .eq. 1 .and. tauir .lt. taumintmp) then
                         tb(ibox) = attrop
                         tau(ibox) = 2.13*taumintmp
                end if

            else

                !clear sky brightness temperature
                tb(ibox) = 1307.27/(log(1.+(1./fluxtopclrsky)))

            end if
            
      end do

      end if
!
!     ---------------------------------------------------!



!     
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use   the physical cloud top pressure
!     (top_height = 2)
!     or the radiatively determined cloud top pressure 
!     (top_height = 1 or 3)
!



      !compute cloud top pressure
      do ibox=1,ncol
      
           !segregate according to optical thickness
           if (tau(ibox).le. (-1.*log((1.-qmin)))) then

                ptop(ibox)=0.
                levmatch(ibox)=0      

           else 

                if (top_height .eq. 1 .or. top_height .eq. 3) then
                                               
                     !find level whose temperature
                     !most closely matches brightness temperature
                     nmatch=0
                     do ilev=1,nlev-1
                        
                          if ((at(ilev)   .ge. tb(ibox) .and.       &
                               at(ilev+1) .lt. tb(ibox)) .or.       &
                              (at(ilev) .le. tb(ibox) .and.         &
                               at(ilev+1) .gt. tb(ibox))) then 
   
                               nmatch=nmatch+1
                               if(abs(at(ilev)-tb(ibox)) .lt.      &
                                  abs(at(ilev+1)-tb(ibox))) then
                                    match(nmatch)=ilev
                               else
                                    match(nmatch)=ilev+1
                               end if
                          end if                        
                     end do

                     if (nmatch .ge. 1) then
                                 
                          ptop(ibox)=pfull(match(nmatch))
                          levmatch(ibox)=match(nmatch)   
                     else
                                                        
                          if (tb(ibox) .lt. atmin) then
                               ptop(ibox)=ptrop
                               levmatch(ibox)=itrop
                          end if
                          if (tb(ibox) .gt. atmax) then
                               ptop(ibox)=pfull(nlev)
                               levmatch(ibox)=nlev
                          end if                                
                                                                
                     end if
                                                               
                else             

                     ptop(ibox)=0.
                     ilev=1
                     do while(ptop(ibox) .eq. 0. .and. ilev .lt. nlev+1)
                          if (frac_out(ibox,ilev) .ne. 0) then
                               ptop(ibox)=pfull(ilev)
                               levmatch(ibox)=ilev
                          end if
                          ilev=ilev+1              
                     end do
                end if                            
           end if
      
      end do
              

!
!     ---------------------------------------------------!



!     
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYP FREQUENCIES
!
!     Now that ptop and tau have been determined, 
!     determine amount of each of the 49 ISCCP cloud
!     types
!



      !compute isccp frequencies
      fq_isccp(:,:)=0. 
      
      !compute boxarea
      boxarea = 1./real(ncol)

      do ibox=1,ncol
      
            !convert ptop to millibars
            ptop(ibox)=ptop(ibox) / 100.
            
      if (tau(ibox) .gt. (-1.*log((1.-qmin))) .and. ptop(ibox) .gt. 0.) then


            !reset itau, ipres
            itau = 0
            ipres = 0


            !determine optical depth category
            if (tau(ibox) .lt. isccp_taumin) then
                itau=1
            else if (tau(ibox).ge.isccp_taumin .and.tau(ibox) .lt. 1.3) then
                itau=2
            else if (tau(ibox) .ge. 1.3 .and. tau(ibox) .lt. 3.6) then
                itau=3
            else if (tau(ibox) .ge. 3.6 .and. tau(ibox) .lt. 9.4) then
                itau=4
            else if (tau(ibox) .ge. 9.4 .and. tau(ibox) .lt. 23.) then
                itau=5
            else if (tau(ibox) .ge. 23. .and. tau(ibox) .lt. 60.) then
                itau=6
            else if (tau(ibox) .ge. 60.) then
                itau=7
            end if

            !determine cloud top pressure category
            if (ptop(ibox) .gt. 0. .and. ptop(ibox) .lt. 180.) then
                ipres=1
            else if (ptop(ibox) .ge. 180..and.ptop(ibox) .lt. 310.) then
                ipres=2
            else if (ptop(ibox) .ge. 310..and.ptop(ibox) .lt. 440.) then
                ipres=3
            else if (ptop(ibox) .ge. 440..and.ptop(ibox) .lt. 560.) then
                ipres=4
            else if(ptop(ibox) .ge. 560..and.ptop(ibox) .lt. 680.) then
                ipres=5
            else if (ptop(ibox) .ge. 680..and.ptop(ibox) .lt. 800.) then
                ipres=6
            else if (ptop(ibox) .ge. 800.) then
                ipres=7
            end if 

            !update frequencies
            if(ipres .gt. 0.and.itau .gt. 0) then
            fq_isccp(itau,ipres)=fq_isccp(itau,ipres)+boxarea 
            end if

      end if
                       
      end do

!---------------------------------------------------------------------
      

end subroutine isccp_cloudtypes

!###################################################################

! <SUBROUTINE NAME="isccp_clouds_end">
!  <OVERVIEW>
!    isccp_clouds_end is the destructor for isccp_clouds_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    isccp_clouds_end is the destructor for isccp_clouds_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call isccp_clouds_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine isccp_clouds_end

!-------------------------------------------------------------------
!    isccp_clouds_end is the destructor for isccp_clouds_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('isccp_clouds_mod',   &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    mark the module as not initialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------



end subroutine isccp_clouds_end


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
!   call diag_field_init (Time, axes )
!  </TEMPLATE>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   initialization time for the netcdf output fields
!  </IN>
! </SUBROUTINE>
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


!-----------------------------------------------------------------------
!    register the isccp diagnostic fields with diag_manager_mod.
!-----------------------------------------------------------------------
        id_pc1tau1 = register_diag_field ( mod_name, &
                             'pc1tau1', axes(1:2), Time, &
                      '    pc<180;     0<tau<taumin', 'fraction' )
        id_pc1tau2 = register_diag_field ( mod_name, &
                             'pc1tau2', axes(1:2), Time, &
                      '    pc<180;taumin<tau<1.3   ', 'fraction' )
        id_pc1tau3 = register_diag_field ( mod_name, &
                             'pc1tau3', axes(1:2), Time, &
                      '    pc<180;   1.3<tau<3.6   ', 'fraction' )
        id_pc1tau4 = register_diag_field ( mod_name, &
                             'pc1tau4', axes(1:2), Time, &
                     '    pc<180;   3.6<tau<9.4   ', 'fraction' )
        id_pc1tau5 = register_diag_field ( mod_name, &
                             'pc1tau5', axes(1:2), Time, &
                      '    pc<180;   9.4<tau<23    ', 'fraction' )
        id_pc1tau6 = register_diag_field ( mod_name, &
                             'pc1tau6', axes(1:2), Time, &
                      '    pc<180;    23<tau<60    ', 'fraction' )
        id_pc1tau7 = register_diag_field ( mod_name, &
                             'pc1tau7', axes(1:2), Time, &
                      '    pc<180;    60<tau       ', 'fraction' )
        id_pc2tau1 = register_diag_field ( mod_name, &
                             'pc2tau1', axes(1:2), Time, &
                      '180<pc<310;     0<tau<taumin', 'fraction' )
        id_pc2tau2 = register_diag_field ( mod_name, &
                             'pc2tau2', axes(1:2), Time, &
                      '180<pc<310;taumin<tau<1.3   ', 'fraction' )
        id_pc2tau3 = register_diag_field ( mod_name, &
                             'pc2tau3', axes(1:2), Time, &
                      '180<pc<310;   1.3<tau<3.6   ', 'fraction' )
        id_pc2tau4 = register_diag_field ( mod_name, &
                             'pc2tau4', axes(1:2), Time, &
                      '180<pc<310;   3.6<tau<9.4   ', 'fraction' )
        id_pc2tau5 = register_diag_field ( mod_name, &
                             'pc2tau5', axes(1:2), Time, &
                      '180<pc<310;   9.4<tau<23    ', 'fraction' )
        id_pc2tau6 = register_diag_field ( mod_name, &
                             'pc2tau6', axes(1:2), Time, &
                      '180<pc<310;    23<tau<60    ', 'fraction' )
        id_pc2tau7 = register_diag_field ( mod_name, &
                             'pc2tau7', axes(1:2), Time, &
                      '180<pc<310;    60<tau       ', 'fraction' )
        id_pc3tau1 = register_diag_field ( mod_name, &
                             'pc3tau1', axes(1:2), Time, &
                      '310<pc<440;     0<tau<taumin', 'fraction' )
        id_pc3tau2 = register_diag_field ( mod_name, &
                             'pc3tau2', axes(1:2), Time, &
                      '310<pc<440;taumin<tau<1.3   ', 'fraction' )
        id_pc3tau3 = register_diag_field ( mod_name, &
                             'pc3tau3', axes(1:2), Time, &
                      '310<pc<440;   1.3<tau<3.6   ', 'fraction' )
        id_pc3tau4 = register_diag_field ( mod_name, &
                             'pc3tau4', axes(1:2), Time, &
                      '310<pc<440;   3.6<tau<9.4   ', 'fraction' )
        id_pc3tau5 = register_diag_field ( mod_name, &
                             'pc3tau5', axes(1:2), Time, &
                      '310<pc<440;   9.4<tau<23    ', 'fraction' )
        id_pc3tau6 = register_diag_field ( mod_name, &
                             'pc3tau6', axes(1:2), Time, &
                      '310<pc<440;    23<tau<60    ', 'fraction' )
        id_pc3tau7 = register_diag_field ( mod_name, &
                             'pc3tau7', axes(1:2), Time, &
                      '310<pc<440;    60<tau       ', 'fraction' )
        id_pc4tau1 = register_diag_field ( mod_name, &
                             'pc4tau1', axes(1:2), Time, &
                      '440<pc<560;     0<tau<taumin', 'fraction' )
        id_pc4tau2 = register_diag_field ( mod_name, &
                             'pc4tau2', axes(1:2), Time, &
                      '440<pc<560;taumin<tau<1.3   ', 'fraction' )
        id_pc4tau3 = register_diag_field ( mod_name, &
                             'pc4tau3', axes(1:2), Time, &
                      '440<pc<560;   1.3<tau<3.6   ', 'fraction' )
        id_pc4tau4 = register_diag_field ( mod_name, &
                             'pc4tau4', axes(1:2), Time, &
                      '440<pc<560;   3.6<tau<9.4   ', 'fraction' )
        id_pc4tau5 = register_diag_field ( mod_name, &
                             'pc4tau5', axes(1:2), Time, &
                      '440<pc<560;   9.4<tau<23    ', 'fraction' )
        id_pc4tau6 = register_diag_field ( mod_name, &
                             'pc4tau6', axes(1:2), Time, &
                      '440<pc<560;    23<tau<60    ', 'fraction' )
        id_pc4tau7 = register_diag_field ( mod_name, &
                             'pc4tau7', axes(1:2), Time, &
                      '440<pc<560;    60<tau       ', 'fraction' )
        id_pc5tau1 = register_diag_field ( mod_name, &
                             'pc5tau1', axes(1:2), Time, &
                      '560<pc<680;     0<tau<taumin', 'fraction' )
        id_pc5tau2 = register_diag_field ( mod_name, &
                             'pc5tau2', axes(1:2), Time, &
                      '560<pc<680;taumin<tau<1.3   ', 'fraction' )
        id_pc5tau3 = register_diag_field ( mod_name, &
                             'pc5tau3', axes(1:2), Time, &
                      '560<pc<680;   1.3<tau<3.6   ', 'fraction' )
        id_pc5tau4 = register_diag_field ( mod_name, &
                             'pc5tau4', axes(1:2), Time, &
                      '560<pc<680;   3.6<tau<9.4   ', 'fraction' )
        id_pc5tau5 = register_diag_field ( mod_name, &
                             'pc5tau5', axes(1:2), Time, &
                      '560<pc<680;   9.4<tau<23    ', 'fraction' )
        id_pc5tau6 = register_diag_field ( mod_name, &
                             'pc5tau6', axes(1:2), Time, &
                      '560<pc<680;    23<tau<60    ', 'fraction' )
        id_pc5tau7 = register_diag_field ( mod_name, &
                            'pc5tau7', axes(1:2), Time, &
                      '560<pc<680;    60<tau       ', 'fraction' )
        id_pc6tau1 = register_diag_field ( mod_name, &
                             'pc6tau1', axes(1:2), Time, &
                      '680<pc<800;     0<tau<taumin', 'fraction' )
        id_pc6tau2 = register_diag_field ( mod_name, &
                             'pc6tau2', axes(1:2), Time, &
                     '680<pc<800;taumin<tau<1.3   ', 'fraction' )
        id_pc6tau3 = register_diag_field ( mod_name, &
                             'pc6tau3', axes(1:2), Time, &
                      '680<pc<800;   1.3<tau<3.6   ', 'fraction' )
        id_pc6tau4 = register_diag_field ( mod_name, &
                             'pc6tau4', axes(1:2), Time, &
                      '680<pc<800;   3.6<tau<9.4   ', 'fraction' )
        id_pc6tau5 = register_diag_field ( mod_name, &
                             'pc6tau5', axes(1:2), Time, &
                      '680<pc<800;   9.4<tau<23    ', 'fraction' )
        id_pc6tau6 = register_diag_field ( mod_name, &
                             'pc6tau6', axes(1:2), Time, &
                      '680<pc<800;    23<tau<60    ', 'fraction' )
        id_pc6tau7 = register_diag_field ( mod_name, &
                             'pc6tau7', axes(1:2), Time, &
                      '680<pc<800;    60<tau       ', 'fraction' )
        id_pc7tau1 = register_diag_field ( mod_name, &
                             'pc7tau1', axes(1:2), Time, &
                      '680<pc<800;     0<tau<taumin', 'fraction' )
        id_pc7tau2 = register_diag_field ( mod_name, &
                             'pc7tau2', axes(1:2), Time, &
                      '800<pc    ;taumin<tau<1.3   ', 'fraction' )
        id_pc7tau3 = register_diag_field ( mod_name, &
                             'pc7tau3', axes(1:2), Time, &
                      '800<pc    ;   1.3<tau<3.6   ', 'fraction' )
        id_pc7tau4 = register_diag_field ( mod_name, &
                             'pc7tau4', axes(1:2), Time, &
                      '800<pc    ;   3.6<tau<9.4   ', 'fraction' )
        id_pc7tau5 = register_diag_field ( mod_name, &
                             'pc7tau5', axes(1:2), Time, &
                      '800<pc    ;   9.4<tau<23    ', 'fraction' )
        id_pc7tau6 = register_diag_field ( mod_name, &
                             'pc7tau6', axes(1:2), Time, &
                      '800<pc    ;    23<tau<60    ', 'fraction' )
        id_pc7tau7 = register_diag_field ( mod_name, &
                             'pc7tau7', axes(1:2), Time, &
                      '800<pc    ;    60<tau       ', 'fraction' )
        id_nisccp = register_diag_field ( mod_name, &
                             'nisccp', axes(1:2), Time, &
                      'frequency of ISCCP calculations', 'fraction' )
       
!----------------------------------------------------------------------
!    if any of the isccp cloudtype variables are desired, set a flag
!    indicating that the model clouds are to be mapped to isccp cloud 
!    categories.
!----------------------------------------------------------------------
        if (id_pc1tau1>0 .or. id_pc1tau2>0 .or.    &
            id_pc1tau3>0 .or. id_pc1tau4>0 .or.    &       
            id_pc1tau5>0 .or. id_pc1tau6>0 .or.    &  
            id_pc1tau7>0) then
          do_isccp=.true.
        else if (id_pc2tau1>0 .or. id_pc2tau2>0 .or.    &
                 id_pc2tau3>0 .or. id_pc2tau4>0 .or.    &       
                 id_pc2tau5>0 .or. id_pc2tau6>0 .or.    &  
                 id_pc2tau7>0) then              
          do_isccp=.true.
        else if (id_pc3tau1>0 .or. id_pc3tau2>0 .or.    &
                 id_pc3tau3>0 .or. id_pc3tau4>0 .or.    &       
                 id_pc3tau5>0 .or. id_pc3tau6>0 .or.    &  
                 id_pc3tau7>0) then              
          do_isccp=.true.
        else if (id_pc4tau1>0 .or. id_pc4tau2>0 .or.    &
                 id_pc4tau3>0 .or. id_pc4tau4>0 .or.    &       
                 id_pc4tau5>0 .or. id_pc4tau6>0 .or.    &  
                 id_pc4tau7>0) then              
          do_isccp=.true.
        else if (id_pc5tau1>0 .or. id_pc5tau2>0 .or.    &
                 id_pc5tau3>0 .or. id_pc5tau4>0 .or.    &       
                 id_pc5tau5>0 .or. id_pc5tau6>0 .or.    &  
                 id_pc5tau7>0) then              
          do_isccp=.true.
        else if (id_pc6tau1>0 .or. id_pc6tau2>0 .or.    &
                 id_pc6tau3>0 .or. id_pc6tau4>0 .or.    &       
                 id_pc6tau5>0 .or. id_pc6tau6>0 .or.    &  
                 id_pc6tau7>0) then              
          do_isccp=.true.
        else if (id_pc7tau1>0 .or. id_pc7tau2>0 .or.    &
                 id_pc7tau3>0 .or. id_pc7tau4>0 .or.    &       
                 id_pc7tau5>0 .or. id_pc7tau6>0 .or.    &  
                 id_pc7tau7>0) then              
          do_isccp=.true.
        endif

!---------------------------------------------------------------------
!    register the isccp cloud summary diagnostics.
!---------------------------------------------------------------------
        id_aice = register_diag_field ( mod_name, &
                             'aice', axes(1:2), Time, &
                      'area of ice clouds seen from space', 'fraction' )
        id_reffice = register_diag_field ( mod_name, &
                         'reffice', axes(1:2), Time, &
                      'ice cloud effective radius', 'microns' )
        id_aliq = register_diag_field ( mod_name, &
                             'aliq', axes(1:2), Time, &
                   'area of liquid clouds seen from space', 'fraction' )
        id_reffliq = register_diag_field ( mod_name, &
                         'reffliq', axes(1:2), Time, &
                      'liquid cloud effective radius', 'microns' )
        id_alow = register_diag_field ( mod_name, &
                             'alow', axes(1:2), Time, &
                      'area of low clouds seen from space', 'fraction' )
        id_tauicelow = register_diag_field ( mod_name, &
                             'tauicelow', axes(1:2), Time, &
                      'low ice cloud optical depth', 'fraction' )
        id_tauliqlow = register_diag_field ( mod_name, &
                             'tauliqlow', axes(1:2), Time, &
                      'low liquid cloud optical depth', 'fraction' )
        id_tlaylow = register_diag_field ( mod_name, &
                             'tlaylow', axes(1:2), Time, &
                      'low layer temperature', 'fraction' )
        id_tcldlow = register_diag_field ( mod_name, &
                             'tcldlow', axes(1:2), Time, &
                      'low cloud top temperature', 'fraction' )
      
!---------------------------------------------------------------------
!    if any of the isccp cloud summary diagnostics are desired, set a
!    logical flag to so indicate.
!----------------------------------------------------------------------
        if (id_aice      > 0  .or. id_reffice   > 0   .or.   &
            id_aliq      > 0  .or. id_reffliq   > 0  .or.    &
            id_alow      > 0  .or. id_tauicelow > 0  .or.    &
            id_tauliqlow > 0  .or. id_tlaylow   > 0  .or.    &
            id_tcldlow   > 0 ) then
           do_tau_reff=.true.
        endif

!---------------------------------------------------------------------

 
end subroutine diag_field_init




!######################################################################
! <FUNCTION NAME="ran0">
!  <OVERVIEW>
!    ran0 is a platform-independent random number generator from
!    Numerical Recipes -- Mark Webb July 1999
!  </OVERVIEW>
!  <DESCRIPTION>
!    ran0 is a platform-independent random number generator from
!    Numerical Recipes -- Mark Webb July 1999
!  </DESCRIPTION>
!  <TEMPLATE>
!   x = ran0 (idum)
!  </TEMPLATE>
!  <IN NAME="idum" TYPE="real">
!   seed for random number generator
!  </IN>
! </FUNCTION>
function ran0 (idum)

!--------------------------------------------------------------------- 
!    ran0 is a platform-independent random number generator from
!    Numerical Recipes -- Mark Webb July 1999
!---------------------------------------------------------------------
       
real                    :: ran0               
integer, intent(inout)  :: idum
                                 
!--------------------------------------------------------------------
!  intent(out) variable:
!
!      ran0           random number generated by this function
!
!  intent(inout) variable:
!
!      idum           seed for random number generator
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer :: ia = 16807        ! constant in random number generator
      integer :: im = 2147483647   ! constant in random number generator
      integer :: iq = 127773       ! constant in random number generator
      integer :: ir = 2836         ! constant in random number generator
      real    :: am                ! constant in random number generator
      integer :: k                 ! work variable              

!---------------------------------------------------------------------
!    define a needed  variable.
!---------------------------------------------------------------------
      am = 1.0/im

!---------------------------------------------------------------------
!    verify that the seed is valid.
!---------------------------------------------------------------------
      if (idum == 0) then
        call error_mesg ('isccp_clouds_mod', &
         'ZERO seed not allowed in ran0', FATAL)
      endif
 
!---------------------------------------------------------------------
!    compute next random number in sequence, using input value of
!    idum. return a new value of idum for use on next call.
!---------------------------------------------------------------------
      k = idum/iq
      idum = ia*(idum - k*iq) - ir*k
      if (idum < 0) idum = idum + im
      ran0 = am*idum

!---------------------------------------------------------------------


end function ran0



!#####################################################################

                end module isccp_clouds_mod


