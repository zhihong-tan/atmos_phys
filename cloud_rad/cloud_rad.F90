
MODULE CLOUD_RAD_MOD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	CLOUD RADIATIVE PROPERTIES
!
!       February 2001
!       Contact person: Steve Klein
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This module does one of two things...
!
!       (a)  If using the "old" Fels-Schwarzkopf LW / Lacis-Hansen SW
!       radiation code this module solves for the radiative properties 
!       of every cloud.  In particular it uses either the
!       two stream approach or the delta-Eddington approach
!       to solve for the longwave emmissivity, the ultra-violet-
!       visible reflectivities and absorptions, and the
!       near-infrared reflectivities and absorptions.
!
!       (b)  If using the "new sea-esf" radiation code this returns
!       the cloud amounts, layers in which these occur, and concentrations
!       of microphysical properties, and the effective diameters
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        The module contains the following:
!
!	SUBROUTINES
!
!
!            CLOUD_RAD_INIT
!                        Initializes values of qmin, N_land, and 
!                        N_ocean using values from strat_cloud namelist
!                        as well as reads its own namelist variables. 
!                        In addition, it registed diagnostic fields
!                        if needed, and returns the value of the
!                        cloud overlap to strat_cloud.
!            CLOUD_SUMMARY
!                        This is the main driver program of the module
!            CLOUD_ORGANIZE
!                        for each cloud this computes the cloud amount,
!                        the liquid and ice water paths, the effective
!                        particle sizes, the tops and bottoms of clouds.
!    	     CLOUD_RAD   
!                        this solves for the radiative properties of the
!                        clouds given the cloud optical properties
!                        (tau,w0,gg) for each cloud using either a 
!                        Delta-Eddington solution (default) or the
!                        two stream approximation.
!            CLOUD_OPTICAL_PROPERTIES
!                        for each cloud this calculates the mixed phase
!                        values of the optical depth, the single scattering
!                        albedo, and the asymmetry parameter (tau, w0,and g)
!                        for the visible and near infrared bands.  It also
!                        computes the longwave emmissivity of each cloud.
!            ISCCP_CLOUDTYPES
!                        computes the ISCCP view of model clouds
!            TAU_REFF_DIAG
!                        computes diagnostics on cloud optical depth,
!                        effective radius, and cloud temperature.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!       Declare outside modules used, make implicit nones, declare
!       public routines/variables.
!

        USE  Constants_Mod,      ONLY :  Rdgas,Grav,Tfreeze,Dens_h2o
        USE  Utilities_Mod,      ONLY :  File_Exist, Open_File,  &
                                         print_version_number,  &
                                         error_mesg, FATAL,     &
                                         Close_File, get_my_pe, &
					 check_nml_error
        use  diag_manager_mod,    only:  register_diag_field, &
                                         send_data
        use  time_manager_mod,    only:  time_type
  
        IMPLICIT NONE
        PRIVATE

        PUBLIC   CLOUD_RAD_INIT, &
                 CLOUD_SUMMARY, &
                 CLOUD_RAD

     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!                           PARAMETERS OF THE SCHEME
!
!
!       taumin         minimum permissible tau         dimensionless
!
!       qmin           minimum permissible cloud       kg condensate/kgair                  
!                      condensate
!              
!       N_land         number of cloud droplets        1/(m*m*m)
!                      in liquid clouds over land
!
!       N_ocean        number of cloud droplets        1/(m*m*m)
!                      in liquid clouds over ocean
!
!       k_land         ratio of effective radius to    dimensionless
!                      volume radius for continental
!                      air masses
!
!       k_ocean        ratio of effective radius to    dimensionless
!                      volume radius for maritime
!                      air masses
!
!       IMPORTANT NOTE qmin, N_land, N_ocean are initialized with a 
!       call from strat_cloud_init to cloud_rad_init.  This guarantees
!       that both strat_cloud and cloud_rad have the exact same values
!       for these parameters.
!
!      
!                         NAMELIST VARIABLES OF THE SCHEME
!
!       overlap        integer variable indicating which 
!                      overlap assumption to use:
!
!                      overlap = 1. means condensate in 
!                                   adjacent levels is 
!                                   treated as part of 
!                                   the same cloud
!
!                      i.e. maximum-random overlap
!                          
!                      overlap = 2. means condensate in 
!                                   adjacent levels is 
!                                   treated as different 
!                                   clouds
!
!                      i.e. random overlap
!
!
!       l2strem        logical variable indicating whether
!                      what solution for cloud radiative
!                      properties is being used.
!          
!                      l2strem = T  2 stream solution
!                      l2strem = F  Delta-Eddington solution
!
!                      Note that IF l2strem = T then the 
!                      solution does not depend on solar 
!                      zenith angle
!
!                      used in the radiative transfer 
!                      solution module (cloud_rad)
!
!       taucrit        critical optical depth for switching 
!                      direct beam to diffuse beam for 
!                      use in Delta-Eddington solution
!
!                      used in the radiative transfer
!                      solution module (cloud_rad)
!                         
!       NOTE: l2strem and taucrit are only used in subroutine cloud_rad!
!
!
!
!       adjust_top     logical variable indicating whether 
!                      or not to use the code which places 
!                      the top and bottom of the cloud at 
!                      the faces which are most in view from
!                      the top and bottom of the cloud block.  
!                      This is done to avoid undue influence 
!                      to very small cloud fractions.  If 
!                      true this adjustment of tops is 
!                      performed; If false this is not 
!                      performed.
!
!       scale_factor   Factor which multiplies actual cloud
!                      optical depths to account for the
!                      plane-parallel homogenous cloud bias
!                      (e.g. Cahalan effect).
!
!                      used only by strat_cloud_mod
!
!       qamin          minimum permissible cloud       dimensionless                  
!                      fraction 
!
!
!       ------------------------------------------------------------------
!
!                         ISCCP CLOUD PROCESSING
!
!    
!       The following variables only apply if the ISCCP cloud view
!       processing is done.  To turn on ISCCP cloud processing you
!       must require one of the diagnostic fields produced by it from
!       the diag_table. 
!
!
!       adjust_isccp_top_height
!
!                      logical variable indicating whether 
!                      or not to simulate 10.5 micron brightness
!                      temperatures to adjust top heights according
!                      to the emissivity of the cloud. 
!                     
!                      If TRUE:
!
!                      then infrared brightness temperature are 
!                      computed and cloud top pressures are adjusted 
!                      
!                      If FALSE:
!                
!                      reported cloud top pressures is the actual cloud 
!                      top pressure in the model
!
!       ncol           number of columns used in ISCCP cloud type
!                      simulations
! 
!       isccp_taumin   minimum optical depth ISCCP can see
!
!       emsfclw        assumed constant of longwave emissivity of
!                      the surface (fraction)
!
!
!                  PHYSICAL CONSTANTS USED IN THE SCHEME
!
!	
!         constant              definition                  unit
!       ------------   -----------------------------   ---------------
!
!	Grav	       gravitational acceleration      m/(s*s)
!
!	Rdgas	       gas constant of dry air         J/kg air/K
!
!       Tfreeze        Triple point of water           K
!
!       Dens_h2o       density of pure liquid          kg/(m*m*m)
!



REAL,    PRIVATE, PARAMETER :: taumin = 1.E-06
REAL,    PRIVATE            :: qmin = 1.E-10
REAL,    PRIVATE            :: N_land = 3.E+08
REAL,    PRIVATE            :: N_ocean = 1.E+08
REAL,    PRIVATE, PARAMETER :: k_land = 1.143
REAL,    PRIVATE, PARAMETER :: k_ocean = 1.077
LOGICAL, PRIVATE            :: do_tau_reff = .FALSE.
LOGICAL, PRIVATE            :: do_isccp = .FALSE.

INTEGER, PRIVATE            :: overlap = 1
LOGICAL, PRIVATE            :: l2strem = .FALSE.
REAL,    PRIVATE            :: taucrit = 1.
LOGICAL, PRIVATE            :: adjust_top = .TRUE.
REAL,    PRIVATE            :: scale_factor = 0.8
REAL,    PRIVATE            :: qamin = 1.E-2
LOGICAL, PRIVATE            :: adjust_isccp_top_height = .TRUE.
INTEGER, PRIVATE            :: ncol = 50
REAL,    PRIVATE            :: isccp_taumin = 0.1
REAL,    PRIVATE            :: emsfclw = 0.94


!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CREATE NAMELIST
!

        NAMELIST /CLOUD_RAD_NML/ overlap,l2strem,taucrit,adjust_top, &
                                 scale_factor,qamin,&
                                 adjust_isccp_top_height,&
                                 ncol, isccp_taumin, emsfclw

!-----------------------------------------------------------------------
!-------------------- diagnostics fields -------------------------------

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

character(len=5) :: mod_name = 'cloud'

!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       DECLARE VERSION NUMBER OF SCHEME
!
        
        character(len=128) :: version = '$Id: cloud_rad.F90,v 1.6 2001/07/05 17:18:16 fms Exp $'
        character(len=128) :: tag = '$Name: havana $'

! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CONTAINS


!#######################################################################
!#######################################################################


SUBROUTINE CLOUD_RAD_INIT(axes,Time,qmin_in,N_land_in,N_ocean_in,overlap_out)
                               

!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This subroutine initializes values of qmin, N_land, and 
!       N_ocean using values from the strat_cloud_module, 
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!	VARIABLES
!
!
!       --------------
!       OPTIONAL INPUT
!       --------------
!
!	
!         variable              definition                  unit
!       ------------   -----------------------------   ---------------
!
!       axes           axis integers for diagnostics
!
!       Time           time type variable for 
!                      diagnostics
!     
!       qmin_in        input value of minimum per-     kg condensate/ 
!                      missible cloud liquid, ice,     kg air
!                      or fraction                     or fraction
!
!       N_land_in      input value of number of        #/(m*m*m)
!                      of cloud drop per cubic meter
!                      over land
!
!       N_ocean_in     input value of number of        #/(m*m*m)
!                      of cloud drop per cubic meter
!                      over ocean
!
!       ---------------
!       OPTIONAL OUTPUT
!       ---------------
!
!       overlap_out    value of the namelist variable overlap
!
!
!       -------------------
!	INTERNAL VARIABLES:
!       -------------------
!
!       unit,io        namelist integers
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

integer,         intent(in), optional     :: axes(4)
type(time_type), intent(in), optional     :: Time
REAL,     INTENT (IN),  OPTIONAL          :: qmin_in,N_land_in,&
                                             N_ocean_in
INTEGER,  INTENT (OUT), OPTIONAL          :: overlap_out

!  Internal variables
!  ------------------


INTEGER                                  :: unit,io,ierr

!-----------------------------------------------------------------------
!       
!       Code
!
    
        
!-----------------------------------------------------------------------
!
!       Namelist functions

	!read namelist file if it exists
        if ( file_exist('input.nml')) then
        unit = open_file (file='input.nml', action='read')
        ierr=1; do while (ierr /= 0)
           read  (unit, nml=cloud_rad_nml, iostat=io, end=10)
           ierr = check_nml_error(io,'cloud_rad_nml')
        enddo
10      call close_file (unit)
        endif
	
	!write namelist variables to logfile
        unit = Open_File ('logfile.out', action='APPEND')
        if ( get_my_pe() == 0 ) then
             Write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
             Write (unit,nml=CLOUD_RAD_NML)
        endif
        Call Close_File (unit)
       
!-----------------------------------------------------------------------
!
!       Prevent unreasonable values

        if (overlap.ne.1 .and. overlap.ne.2) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                                'overlap must be either 1 or 2 ', FATAL)
        if (taucrit .lt. 0.) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                  'taucrit must be greater than or equal to 0. ', FATAL)
        if (scale_factor .lt. 0.) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                         'scale_factor must be greater than 0. ', FATAL)
        if (qamin .le. 0. .or. qamin .ge. 1.) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                               'qamin must be between 0. and 1.', FATAL)
        
!-----------------------------------------------------------------------
!
!       Assign values

        if (present(qmin_in)) then
              qmin = qmin_in
        end if
        if (present(N_land_in)) then
              N_land = N_land_in
        end if
        if (present(N_ocean_in)) then
              N_ocean = N_ocean_in
        end if
        if (present(overlap_out)) then
              overlap_out = overlap
        end if
        
!-----------------------------------------------------------------------
!
!
!       NETCDF
!
!

   if(PRESENT(axes).and.PRESENT(Time)) then

!      ISCCP
!      -----

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
                      'frequency of sunlit times', 'fraction' )
       
       !set do_isccp
       do_isccp = .false.
       if (id_pc1tau1>0) do_isccp=.true.; if (id_pc1tau2>0) do_isccp=.true.
       if (id_pc1tau3>0) do_isccp=.true.; if (id_pc1tau4>0) do_isccp=.true.
       if (id_pc1tau5>0) do_isccp=.true.; if (id_pc1tau6>0) do_isccp=.true.
       if (id_pc1tau7>0) do_isccp=.true.
       if (id_pc2tau1>0) do_isccp=.true.; if (id_pc2tau2>0) do_isccp=.true.
       if (id_pc2tau3>0) do_isccp=.true.; if (id_pc2tau4>0) do_isccp=.true.
       if (id_pc2tau5>0) do_isccp=.true.; if (id_pc2tau6>0) do_isccp=.true.
       if (id_pc2tau7>0) do_isccp=.true.
       if (id_pc3tau1>0) do_isccp=.true.; if (id_pc3tau2>0) do_isccp=.true.
       if (id_pc3tau3>0) do_isccp=.true.; if (id_pc3tau4>0) do_isccp=.true.
       if (id_pc3tau5>0) do_isccp=.true.; if (id_pc3tau6>0) do_isccp=.true.
       if (id_pc3tau7>0) do_isccp=.true.
       if (id_pc4tau1>0) do_isccp=.true.; if (id_pc4tau2>0) do_isccp=.true.
       if (id_pc4tau3>0) do_isccp=.true.; if (id_pc4tau4>0) do_isccp=.true.
       if (id_pc4tau5>0) do_isccp=.true.; if (id_pc4tau6>0) do_isccp=.true.
       if (id_pc4tau7>0) do_isccp=.true.
       if (id_pc5tau1>0) do_isccp=.true.; if (id_pc5tau2>0) do_isccp=.true.
       if (id_pc5tau3>0) do_isccp=.true.; if (id_pc5tau4>0) do_isccp=.true.
       if (id_pc5tau5>0) do_isccp=.true.; if (id_pc5tau6>0) do_isccp=.true.
       if (id_pc5tau7>0) do_isccp=.true.
       if (id_pc5tau1>0) do_isccp=.true.; if (id_pc6tau2>0) do_isccp=.true.
       if (id_pc6tau3>0) do_isccp=.true.; if (id_pc6tau4>0) do_isccp=.true.
       if (id_pc6tau5>0) do_isccp=.true.; if (id_pc6tau6>0) do_isccp=.true.
       if (id_pc6tau7>0) do_isccp=.true.
       if (id_pc7tau1>0) do_isccp=.true.; if (id_pc7tau2>0) do_isccp=.true.
       if (id_pc7tau3>0) do_isccp=.true.; if (id_pc7tau4>0) do_isccp=.true.
       if (id_pc7tau5>0) do_isccp=.true.; if (id_pc7tau6>0) do_isccp=.true.
       if (id_pc7tau7>0) do_isccp=.true.
      

!      TAU_REFF
!      --------

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
       
       if (id_aice>0) do_tau_reff=.true.; if (id_reffice>0) do_tau_reff=.true.
       if (id_aliq>0) do_tau_reff=.true.; if (id_reffliq>0) do_tau_reff=.true.
       if (id_alow>0) do_tau_reff=.true.
       if (id_tauicelow>0) do_tau_reff=.true.
       if (id_tauliqlow>0) do_tau_reff=.true.
       if (id_tlaylow>0) do_tau_reff=.true.
       if (id_tcldlow>0) do_tau_reff=.true.

   end if !for setting up diagnostics

!-----------------------------------------------------------------------
!
!
!       END OF SUBROUTINE
!
!


END SUBROUTINE CLOUD_RAD_INIT


!########################################################################
!########################################################################

SUBROUTINE CLOUD_SUMMARY(is,js,&
                  LAND,ql,qi,qa,qv,pfull,phalf,TKel,coszen,skt,&
                  nclds,ktop,kbot,cldamt,Time,&
                  r_uv,r_nir,ab_uv,ab_nir,em_lw,&
                  conc_drop,conc_ice,size_drop,size_ice)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine returns the following properties of clouds
!
!               1. nclds: # of clouds
!               2. ktop : integer level for top of cloud
!               3. kbot : integer level for bottom of cloud
!               4. cldamt:horizontal cloud amount of every cloud
!
!      Optional arguments
!               5. r_uv : cloud reflectance in uv band
!               6. r_nir: cloud reflectance in nir band
!               7. ab_uv: cloud absorption in uv band
!               8. ab_nir:cloud absorption in nir band
!               9. em_lw :longwave cloud emmissivity
!              10. conc_drop : liquid cloud droplet mass concentration
!              11. conc_ice  : ice cloud mass concentration
!              12. size_drop : effective diameter of liquid cloud droplets
!              13. size_ice  : effective diameter of ice cloud 
!
!      given inputs of ql and qi (liquid and ice condensate),
!      cloud volume fraction, and pressure at the half and full levels
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	VARIABLES
!
!       ------
!	INPUT:
!       ------
!
!       is,js        indices for model slab
!       LAND         fraction of the grid box covered by LAND
!       ql           cloud liquid condensate (kg condensate/kg air)
!       qi           cloud ice condensate (kg condensate/kg air)
!       qa           cloud volume fraction (fraction)
!       qv           water vapor specific humidity (kg vapor/kg air)
!       pfull        pressure at full levels (Pascals)
!       phalf        pressure at half levels (Pascals)
!                     NOTE: it is assumed that phalf(j+1) > phalf(j)
!       TKel            temperature (Kelvin)
!       coszen       cosine of the zenith angle
!       skt          surface skin temperature (Kelvin)
!
!       -------------
!	INPUT/OUTPUT:
!       -------------
!
!       nclds        number of (random overlapping) clouds in column and also
!                        the current # for clouds to be operating on
!       ktop         level of the top of the cloud
!       kbot         level of the bottom of the cloud
!       cldamt       cloud amount of condensed cloud
!
!       ---------------------
!       OPTIONAL INPUT/OUTPUT
!       ---------------------
!
!       r_uv         cloud reflectance in uv band
!       r_nir        cloud reflectance in nir band
!       ab_uv        cloud absorption in uv band
!       ab_nir       cloud absorption in nir band
!       em_lw        longwave cloud emmissivity
!       conc_drop    liquid cloud droplet mass concentration (g /m3)
!       conc_ice     ice cloud mass concentration (g /m3)
!       size_drop    effective diameter of liquid cloud droplets (microns)
!       size_ice   : effective diameter of ice clouds (microns)
!
!       -------------------
!	INTERNAL VARIABLES:
!       -------------------
!
!       i,j,k,t      looping variables
!       IDIM         number of first dimension points
!       JDIM         number of second dimension points
!       KDIM         number of third dimension points
!       N_drop       number of cloud droplets per cubic meter
!       k_ratio      ratio of effective radius to mean volume radius
!       max_cld      maximum number of clouds in whole array
!       qa_local     local value of qa (fraction)
!       ql_local     local value of ql (kg condensate / kg air)
!       qi_local     local value of qi (kg condensate / kg air)
!       LWP          cloud liquid water path (kg condensate per square meter)
!       IWP          cloud ice path (kg condensate per square meter)
!       Reff_liq     effective radius for liquid clouds (microns)
!       Reff_ice     effective particle size for ice clouds (microns)
!       tau          optical depth in 4 bands (dimensionless)
!       w0           single scattering albedo in 4 bands (dimensionless)
!       gg           asymmetry parameter in 4 bands (dimensionless)
!       rad_prop     logical indicating if you are requesting the
!                    radiative properties of the clouds
!       wat_prop     logical determining if you are requesting the
!                    concentrations and particle sizes of clouds
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

INTEGER,  INTENT (IN)                    :: is,js
REAL,     INTENT (IN), DIMENSION(:,:)    :: LAND,skt
REAL,     INTENT (IN), DIMENSION(:,:)    :: coszen
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: ql,qi,qa,qv,pfull,TKel
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: phalf
INTEGER,  INTENT (INOUT),DIMENSION(:,:)  :: nclds
INTEGER,  INTENT (INOUT),DIMENSION(:,:,:):: ktop,kbot
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: cldamt
type(time_type), intent(in), optional    :: Time
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: r_uv,r_nir,ab_uv,ab_nir,em_lw
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: conc_drop,conc_ice
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: size_drop,size_ice



!  Internal variables
!  ------------------

INTEGER                                           :: i,j,k,IDIM,JDIM,KDIM,max_cld
LOGICAL                                           :: rad_prop, wat_prop
REAL, DIMENSION(SIZE(ql,1),SIZE(ql,2),SIZE(ql,3)) :: qa_local,ql_local,qi_local
REAL, DIMENSION(SIZE(ql,1),SIZE(ql,2))            :: N_drop, k_ratio
REAL, DIMENSION(:,:,:), allocatable               :: r_uv_local, r_nir_local
REAL, DIMENSION(:,:,:), allocatable               :: ab_uv_local, ab_nir_local
REAL, DIMENSION(:,:,:), allocatable               :: em_lw_local
REAL, DIMENSION(:,:,:), allocatable               :: conc_drop_local,conc_ice_local
REAL, DIMENSION(:,:,:), allocatable               :: size_drop_local,size_ice_local
REAL, DIMENSION(:,:,:), allocatable               :: LWP,IWP,Reff_liq,Reff_ice
REAL, DIMENSION(:,:,:,:), allocatable             :: tau,tau_ice,w0,gg
REAL, DIMENSION(:,:,:,:), allocatable             :: fq_isccp
REAL, DIMENSION(:,:), allocatable                 :: npoints
REAL, DIMENSION(:), allocatable                   :: tau_local,em_local,cldamt_local
LOGICAL                                           :: used

!
! Code
! ----

    
    ! reinitialize variables
    IDIM=SIZE(ql,1)
    JDIM=SIZE(ql,2)
    KDIM=SIZE(ql,3)
    nclds(:,:)    = 0
    ktop(:,:,:)   = 0
    kbot(:,:,:)   = 0
    cldamt(:,:,:) = 0.

    rad_prop = .FALSE.
    wat_prop = .FALSE.
    if (PRESENT(r_uv).or.PRESENT(r_nir).or.PRESENT(ab_uv).or.PRESENT(ab_nir).or.&
        PRESENT(em_lw)) then
        rad_prop = .TRUE.
    end if
    if (PRESENT(conc_drop).or.PRESENT(conc_ice).or.PRESENT(size_drop).or.&
        PRESENT(size_ice)) then
        wat_prop = .TRUE.
    end if
    if ((.not.rad_prop).and.(.not.wat_prop)) then
        rad_prop = .TRUE.
    end if

    !create local values of ql and qi
    !this step is necessary to remove the values of (qi,ql) which are
    !           0 < (qi,ql) < qmin   or
    !               (qi,ql) > qmin and qa <= qamin

    ql_local(:,:,:) = 0.
    qi_local(:,:,:) = 0.
    qa_local(:,:,:) = 0.
    WHERE ( (qa(:,:,:) .gt. qamin) .and. (ql(:,:,:) .gt. qmin) )
                  ql_local(:,:,:) = ql(:,:,:)
                  qa_local(:,:,:) = qa(:,:,:)
    END WHERE
    WHERE ( (qa(:,:,:) .gt. qamin) .and. (qi(:,:,:) .gt. qmin) )
                  qi_local(:,:,:) = qi(:,:,:)
                  qa_local(:,:,:) = qa(:,:,:)
    END WHERE

    !compute N_drop and k_ratio
    N_drop(:,:)=N_land*LAND(:,:) + N_ocean*(1.-LAND(:,:))
    k_ratio(:,:)=k_land*LAND(:,:) + k_ocean*(1.-LAND(:,:))



    !do solution for new radiation code
    if (wat_prop) then

         ALLOCATE(conc_drop_local(IDIM,JDIM,KDIM))
         ALLOCATE(conc_ice_local(IDIM,JDIM,KDIM))
         ALLOCATE(size_drop_local(IDIM,JDIM,KDIM))
         ALLOCATE(size_ice_local(IDIM,JDIM,KDIM))
         
         conc_drop_local(:,:,:) = 0.
         conc_ice_local(:,:,:)  = 0.
         size_drop_local(:,:,:) = 20.
         size_ice_local(:,:,:)  = 60.

         call  cloud_organize(ql_local,qi_local,qa_local,&
                  pfull,phalf,TKel,coszen,N_drop,&
                  k_ratio,nclds,ktop,kbot,cldamt,&
                  conc_drop_org=conc_drop_local,&
                  conc_ice_org =conc_ice_local,&
                  size_drop_org=size_drop_local,&
                  size_ice_org =size_ice_local)

         !assign to output
         if (PRESENT(conc_drop)) conc_drop = conc_drop_local
         if (PRESENT(conc_ice)) conc_ice = conc_ice_local
         if (PRESENT(size_drop)) size_drop = size_drop_local
         if (PRESENT(size_ice)) size_ice = size_ice_local
    
         DEALLOCATE(conc_drop_local)
         DEALLOCATE(conc_ice_local)
         DEALLOCATE(size_drop_local)
         DEALLOCATE(size_ice_local)
         
    end if
     
         
    !do solution for old radiation code
    if (rad_prop) then

    
         ALLOCATE(r_uv_local(IDIM,JDIM,KDIM))
         ALLOCATE(r_nir_local(IDIM,JDIM,KDIM))
         ALLOCATE(ab_uv_local(IDIM,JDIM,KDIM))
         ALLOCATE(ab_nir_local(IDIM,JDIM,KDIM))
         ALLOCATE(em_lw_local(IDIM,JDIM,KDIM))
         ALLOCATE(LWP(IDIM,JDIM,KDIM))
         ALLOCATE(IWP(IDIM,JDIM,KDIM))
         ALLOCATE(Reff_liq(IDIM,JDIM,KDIM))
         ALLOCATE(Reff_ice(IDIM,JDIM,KDIM))
         ALLOCATE(tau(IDIM,JDIM,KDIM,4))
         ALLOCATE(w0(IDIM,JDIM,KDIM,4))
         ALLOCATE(gg(IDIM,JDIM,KDIM,4))
         
         r_uv_local(:,:,:)   = 0.
         r_nir_local(:,:,:)  = 0.
         ab_uv_local(:,:,:)  = 0.
         ab_nir_local(:,:,:) = 0.
         em_lw_local(:,:,:)  = 0.
         LWP(:,:,:)    = 0.
         IWP(:,:,:)    = 0.
         Reff_liq(:,:,:) = 10.
         Reff_ice(:,:,:) = 30.
         tau(:,:,:,:)    = 0.
         w0(:,:,:,:)     = 0.
         gg(:,:,:,:)     = 0.

    
         call  cloud_organize(ql_local,qi_local,qa_local,&
                  pfull,phalf,TKel,coszen,N_drop,&
                  k_ratio,nclds,ktop,kbot,cldamt,&
                  LWP=LWP,IWP=IWP,Reff_liq=Reff_liq,Reff_ice=Reff_ice)
         
         !find maximum number of clouds
         max_cld  = MAXVAL(nclds(:,:))
         
         !compute cloud radiative properties
         IF (max_cld .gt. 0) then

              CALL CLOUD_OPTICAL_PROPERTIES(LWP(:,:,1:max_cld),IWP(:,:,1:max_cld),&
                      Reff_liq(:,:,1:max_cld),Reff_ice(:,:,1:max_cld),&
                      tau(:,:,1:max_cld,:),w0(:,:,1:max_cld,:),gg(:,:,1:max_cld,:),&
                      em_lw_local(:,:,1:max_cld))
              
              !Account for plane-parallel homogenous cloud bias
              tau(:,:,:,:) = scale_factor * tau(:,:,:,:)

              !compute cloud radiative properties
              CALL CLOUD_RAD(tau(:,:,1:max_cld,:),w0(:,:,1:max_cld,:),&
                       gg(:,:,1:max_cld,:),coszen,&
                       r_uv_local(:,:,1:max_cld),r_nir_local(:,:,1:max_cld),&
                       ab_uv_local(:,:,1:max_cld),ab_nir_local(:,:,1:max_cld))
              
              !assure that zero clouds have properties of zero clouds
              DO i = 1, IDIM
              DO j = 1, JDIM
                       if (nclds(i,j).lt.max_cld) then
                               r_uv_local(i,j,nclds(i,j)+1:max_cld)   = 0.
                               r_nir_local(i,j,nclds(i,j)+1:max_cld)  = 0.
                               ab_uv_local(i,j,nclds(i,j)+1:max_cld)  = 0.
                               ab_nir_local(i,j,nclds(i,j)+1:max_cld) = 0.
                               em_lw_local(i,j,nclds(i,j)+1:max_cld)  = 0.
                       end if
              ENDDO
              ENDDO

         END IF

         if (PRESENT(r_uv)) r_uv = r_uv_local
         if (PRESENT(r_nir)) r_nir = r_nir_local
         if (PRESENT(ab_uv)) ab_uv = ab_uv_local
         if (PRESENT(ab_nir)) ab_nir = ab_nir_local
         if (PRESENT(em_lw)) em_lw = em_lw_local
          
         DEALLOCATE(r_uv_local)
         DEALLOCATE(r_nir_local)
         DEALLOCATE(ab_uv_local)
         DEALLOCATE(ab_nir_local)
         DEALLOCATE(em_lw_local)
         DEALLOCATE(LWP)
         DEALLOCATE(IWP)
         DEALLOCATE(Reff_liq)
         DEALLOCATE(Reff_ice)
         DEALLOCATE(tau)
         DEALLOCATE(w0)
         DEALLOCATE(gg)
         
    end if

    !do solution for isccp clouds
    if (do_isccp .or. do_tau_reff) then

         !make sure Time is present
         if (.not.present(Time)) then
              call error_mesg  ('Time has not been passed to cloud',&
                                'summary, needed for diagnostics', FATAL)
         end if
        
         !initialize temporary variables
         ALLOCATE(em_lw_local(IDIM,JDIM,KDIM))
         ALLOCATE(LWP(IDIM,JDIM,KDIM))
         ALLOCATE(IWP(IDIM,JDIM,KDIM))
         ALLOCATE(Reff_liq(IDIM,JDIM,KDIM))
         ALLOCATE(Reff_ice(IDIM,JDIM,KDIM))
         ALLOCATE(tau(IDIM,JDIM,KDIM,4))
         ALLOCATE(tau_ice(IDIM,JDIM,KDIM,4))
         ALLOCATE(w0(IDIM,JDIM,KDIM,4))
         ALLOCATE(gg(IDIM,JDIM,KDIM,4))


         em_lw_local(:,:,:)  = 0.
         LWP(:,:,:)    = 0.
         IWP(:,:,:)    = 0.
         Reff_liq(:,:,:) = 10.
         Reff_ice(:,:,:) = 30.
         tau(:,:,:,:)    = 0.
         tau_ice(:,:,:,:) = 0.
         w0(:,:,:,:)     = 0.
         gg(:,:,:,:)     = 0.

         
         call cloud_organize(ql_local,qi_local,qa_local,&
              pfull,phalf,TKel,coszen,N_drop,&
              k_ratio,nclds,ktop,kbot,cldamt,&
              LWP=LWP,IWP=IWP,Reff_liq=Reff_liq,Reff_ice=Reff_ice)
         
         !find maximum number of clouds
         max_cld  = MAXVAL(nclds(:,:))
         
         if (max_cld .ge. 1) then
         
         CALL CLOUD_OPTICAL_PROPERTIES(          LWP(:,:,1:max_cld)  ,&
                       IWP(:,:,1:max_cld)  ,Reff_liq(:,:,1:max_cld)  ,&
                  Reff_ice(:,:,1:max_cld)  ,     tau(:,:,1:max_cld,:),&
                        w0(:,:,1:max_cld,:),      gg(:,:,1:max_cld,:),&
               em_lw_local(:,:,1:max_cld)  ,  &
                  tau_ice_diag = tau_ice(:,:,1:max_cld,:))
                          
         end if  !max_cld

         if (do_isccp) then

         !allocate temporary arrays
         ALLOCATE(fq_isccp(IDIM,JDIM,7,7))
         ALLOCATE(npoints(IDIM,JDIM))
         ALLOCATE(cldamt_local(KDIM))
         ALLOCATE(tau_local(KDIM))
         ALLOCATE(em_local(KDIM))
         fq_isccp(:,:,:,:) = 0.
         npoints(:,:) = 0.
                
         !---loop over points and nclds----!
         DO i = 1, IDIM
         DO j = 1, JDIM

              !move from cloud space to model level space
              cldamt_local(:) = 0.
              tau_local(:) = 0.
              em_local(:) = 0.
              do k = 1, nclds(i,j)
                    cldamt_local(ktop(i,j,k):kbot(i,j,k)) = cldamt(i,j,k)                  
                    tau_local(ktop(i,j,k):kbot(i,j,k)) = tau(i,j,k,1)/ &
                                          real(kbot(i,j,k)-ktop(i,j,k)+1)
                    em_local(ktop(i,j,k):kbot(i,j,k)) = 1. - &
                         ( (1.-em_lw_local(i,j,k))** &
                           (1./real(kbot(i,j,k)-ktop(i,j,k)+1)) )
              enddo

              if (coszen(i,j) .gt. 1.E-06) then
                    call ISCCP_CLOUDTYPES(pfull(i,j,:),phalf(i,j,:),&
                         qv(i,j,:),TKel(i,j,:),skt(i,j),cldamt_local(:),&
                         tau_local(:),em_local(:),fq_isccp(i,j,:,:))
                    npoints(i,j) = 1.
              end if

         END DO
         END DO
         
         !5. netcdf diagnostics....
         
         used = send_data ( id_pc1tau1, fq_isccp(:,:,1,1), Time, &
                           is, js )
         used = send_data ( id_pc1tau2, fq_isccp(:,:,2,1), Time, &
                           is, js )
         used = send_data ( id_pc1tau3, fq_isccp(:,:,3,1), Time, &
                           is, js )
         used = send_data ( id_pc1tau4, fq_isccp(:,:,4,1), Time, &
                           is, js )
         used = send_data ( id_pc1tau5, fq_isccp(:,:,5,1), Time, &
                           is, js )
         used = send_data ( id_pc1tau6, fq_isccp(:,:,6,1), Time, &
                           is, js )
         used = send_data ( id_pc1tau7, fq_isccp(:,:,7,1), Time, &
                           is, js )
         used = send_data ( id_pc2tau1, fq_isccp(:,:,1,2), Time, &
                           is, js )
         used = send_data ( id_pc2tau2, fq_isccp(:,:,2,2), Time, &
                           is, js )
         used = send_data ( id_pc2tau3, fq_isccp(:,:,3,2), Time, &
                           is, js )
         used = send_data ( id_pc2tau4, fq_isccp(:,:,4,2), Time, &
                           is, js )
         used = send_data ( id_pc2tau5, fq_isccp(:,:,5,2), Time, &
                           is, js )
         used = send_data ( id_pc2tau6, fq_isccp(:,:,6,2), Time, &
                           is, js )
         used = send_data ( id_pc2tau7, fq_isccp(:,:,7,2), Time, &
                           is, js )
         used = send_data ( id_pc3tau1, fq_isccp(:,:,1,3), Time, &
                           is, js )
         used = send_data ( id_pc3tau2, fq_isccp(:,:,2,3), Time, &
                           is, js )
         used = send_data ( id_pc3tau3, fq_isccp(:,:,3,3), Time, &
                           is, js )
         used = send_data ( id_pc3tau4, fq_isccp(:,:,4,3), Time, &
                           is, js )
         used = send_data ( id_pc3tau5, fq_isccp(:,:,5,3), Time, &
                           is, js )
         used = send_data ( id_pc3tau6, fq_isccp(:,:,6,3), Time, &
                           is, js )
         used = send_data ( id_pc3tau7, fq_isccp(:,:,7,3), Time, &
                           is, js )
         used = send_data ( id_pc4tau1, fq_isccp(:,:,1,4), Time, &
                           is, js )
         used = send_data ( id_pc4tau2, fq_isccp(:,:,2,4), Time, &
                           is, js )
         used = send_data ( id_pc4tau3, fq_isccp(:,:,3,4), Time, &
                           is, js )
         used = send_data ( id_pc4tau4, fq_isccp(:,:,4,4), Time, &
                           is, js )
         used = send_data ( id_pc4tau5, fq_isccp(:,:,5,4), Time, &
                           is, js )
         used = send_data ( id_pc4tau6, fq_isccp(:,:,6,4), Time, &
                           is, js )
         used = send_data ( id_pc4tau7, fq_isccp(:,:,7,4), Time, &
                           is, js )
         used = send_data ( id_pc5tau1, fq_isccp(:,:,1,5), Time, &
                           is, js )
         used = send_data ( id_pc5tau2, fq_isccp(:,:,2,5), Time, &
                           is, js )
         used = send_data ( id_pc5tau3, fq_isccp(:,:,3,5), Time, &
                           is, js )
         used = send_data ( id_pc5tau4, fq_isccp(:,:,4,5), Time, &
                           is, js )
         used = send_data ( id_pc5tau5, fq_isccp(:,:,5,5), Time, &
                           is, js )
         used = send_data ( id_pc5tau6, fq_isccp(:,:,6,5), Time, &
                           is, js )
         used = send_data ( id_pc5tau7, fq_isccp(:,:,7,5), Time, &
                           is, js )
         used = send_data ( id_pc6tau1, fq_isccp(:,:,1,6), Time, &
                           is, js )
         used = send_data ( id_pc6tau2, fq_isccp(:,:,2,6), Time, &
                           is, js )
         used = send_data ( id_pc6tau3, fq_isccp(:,:,3,6), Time, &
                           is, js )
         used = send_data ( id_pc6tau4, fq_isccp(:,:,4,6), Time, &
                           is, js )
         used = send_data ( id_pc6tau5, fq_isccp(:,:,5,6), Time, &
                           is, js )
         used = send_data ( id_pc6tau6, fq_isccp(:,:,6,6), Time, &
                           is, js )
         used = send_data ( id_pc6tau7, fq_isccp(:,:,7,6), Time, &
                           is, js )
         used = send_data ( id_pc7tau1, fq_isccp(:,:,1,7), Time, &
                           is, js )
         used = send_data ( id_pc7tau2, fq_isccp(:,:,2,7), Time, &
                           is, js )
         used = send_data ( id_pc7tau3, fq_isccp(:,:,3,7), Time, &
                           is, js )
         used = send_data ( id_pc7tau4, fq_isccp(:,:,4,7), Time, &
                           is, js )
         used = send_data ( id_pc7tau5, fq_isccp(:,:,5,7), Time, &
                           is, js )
         used = send_data ( id_pc7tau6, fq_isccp(:,:,6,7), Time, &
                           is, js )
         used = send_data ( id_pc7tau7, fq_isccp(:,:,7,7), Time, &
                           is, js )
   
         used = send_data ( id_nisccp, npoints, Time, is, js )

         !deallocate temporary arrays
         DEALLOCATE(fq_isccp)
         DEALLOCATE(npoints)
         DEALLOCATE(cldamt_local)
         DEALLOCATE(tau_local)
         DEALLOCATE(em_local)
         
         end if   !for do_isccp

         if (do_tau_reff) then

         call TAU_REFF_DIAG(is,js,Time,coszen,nclds,ktop,&
             kbot,cldamt,TKel,phalf,tau(:,:,:,1),tau_ice(:,:,:,1),&
             Reff_liq,Reff_ice)
         
         end if   !for do_tau_reff
        
         DEALLOCATE(em_lw_local)
         DEALLOCATE(LWP)
         DEALLOCATE(IWP)
         DEALLOCATE(Reff_liq)
         DEALLOCATE(Reff_ice)
         DEALLOCATE(tau)
         DEALLOCATE(tau_ice)
         DEALLOCATE(w0)
         DEALLOCATE(gg)
         
    end if   !for do_isccp .or. do_tau_reff
   
    
    
END SUBROUTINE CLOUD_SUMMARY




!########################################################################
!########################################################################

SUBROUTINE CLOUD_ORGANIZE(ql,qi,qa,pfull,phalf,TKel,coszen,N_drop,&
                  k_ratio,nclds,ktop,kbot,cldamt,&
                  LWP,IWP,Reff_liq,Reff_ice, &
                  conc_drop_org,conc_ice_org,size_drop_org,size_ice_org)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine returns the following properties of clouds
!
!               1. nclds: # of clouds
!               2. ktop : integer level for top of cloud
!               3. kbot : integer level for bottom of cloud
!               4. cldamt:horizontal cloud amount of every cloud
!               5. LWP :
!
!      given inputs of ql and qi (liquid and ice condensate),
!      cloud volume fraction, and pressure at the half and full levels
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	VARIABLES
!
!       ------
!	INPUT:
!       ------
!
!       LAND         fraction of the grid box covered by LAND
!       ql           cloud liquid condensate (kg condensate/kg air)
!       qi           cloud ice condensate (kg condensate/kg air)
!       qa           cloud volume fraction (fraction)
!       pfull        pressure at full levels (Pascals)
!       phalf        pressure at half levels (Pascals)
!                     NOTE: it is assumed that phalf(j+1) > phalf(j)
!       TKel            temperature (Kelvin)
!       coszen       cosine of the zenith angle
!       N_drop       number of cloud droplets per cubic meter
!       k_ratio      ratio of effective radius to mean volume radius
!
!       -------------
!	INPUT/OUTPUT:
!       -------------
!
!       nclds        number of (random overlapping) clouds in column and also
!                        the current # for clouds to be operating on
!       ktop         level of the top of the cloud
!       kbot         level of the bottom of the cloud
!       cldamt       cloud amount of condensed cloud
!
!       ---------------------
!       OPTIONAL INPUT/OUTPUT
!       ---------------------
!
!       LWP          cloud liquid water path (kg condensate per square meter)
!       IWP          cloud ice path (kg condensate per square meter)
!       Reff_liq     effective radius for liquid clouds (microns)
!       Reff_ice     effective particle size for ice clouds (microns)
!       conc_drop_org liquid cloud droplet mass concentration (g /m3)
!       conc_ice_org  ice cloud mass concentration (g /m3)
!       size_drop_org effective diameter of liquid cloud droplets (microns)
!       size_ice_org  effective diameter of ice clouds (microns)
!
!       -------------------
!	INTERNAL VARIABLES:
!       -------------------
!
!       i,j,k,t      looping variables
!       IDIM         number of first dimension points
!       JDIM         number of second dimension points
!       KDIM         number of vertical levels
!       nlev         number of levels in the cloud
!       reff_liq_local   reff of liquid clouds used locally (microns)
!       sum_reff_liq  a sum of reff_liq_local
!       reff_ice_local   reff of ice clouds used locally (microns)
!       sum_reff_ice  a sum of reff_liq_local
!       sum_liq      sum of liquid in cloud (kg condensate per square meter)
!       sum_ice      sum of ice in cloud (kg condensate per square meter)
!       maxcldfrac   maximum cloud fraction in cloud block (fraction)
!       top_t,bot_t  temporary integers used to identify cloud edges
!       totcld_bot   total cloud fraction from bottom view
!       max_bot      largest cloud fraction face from bottom view
!       totcld_top   total cloud fraction from top view
!       max_top      largest cloud fraction face from top view
!       tmp_val      temporary number used in the assigning of top and bottoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      NOTE ON THE FORMULAS FOR EFFECTIVE RADIUS OF LIQUID AND ICE
!      CLOUDS:
!
!
!      FOR LIQUID CLOUDS THE FOLLOWING FORMULA IS USED:
!
!      THIS FORMULA IS THE RECOMMENDATION OF
!      Martin et al., J. Atmos. Sci, vol 51, pp. 1823-1842
!
!
!        reff (in microns) =  k * 1.E+06 *
!                    (3*airdens*(ql/qa)/(4*pi*Dens_h2o*N_liq))**(1/3)
!
!       where airdens = density of air in kg air/m3
!                  ql = liquid condensate in kg cond/kg air
!                  qa = cloud fraction
!                  pi = 3.14159
!            Dens_h2o = density of pure liquid water (kg liq/m3) 
!               N_liq = density of cloud droplets (number per cubic meter)
!                   k = factor to account for difference between 
!                       mean volume radius and effective radius
!
!        IN THIS PROGRAM reff_liq is limited to be between 4.2 microns
!        and 16.6 microns, which is the range of validity for the
!        Slingo (1989) radiation.
!
!     For single layer liquid or mixed phase clouds it is assumed that
!     cloud liquid is vertically stratified within the cloud.  Under
!     such situations for observed stratocumulus clouds it is found
!     that the cloud mean effective radius is between 80 and 100% of
!     the cloud top effective radius. (Brenguier et al., Journal of
!     Atmospheric Sciences, vol. 57, pp. 803-821 (2000))  For linearly 
!     stratified cloud in liquid specific humidity, the cloud top 
!     effective radius is greater than the effective radius of the 
!     cloud mean specific humidity by a factor of 2**(1./3.).
!
!     This correction, 0.9*(2**(1./3.)) = 1.134, is applied only to 
!     single layer liquid or mixed phase clouds.
!
!
!     FOR ICE CLOUDS THE EFFECTIVE RADIUS IS TAKEN FROM THE FORMULATION
!     IN DONNER (1997, J. Geophys. Res., 102, pp. 21745-21768) WHICH IS
!     BASED ON HEYMSFIELD AND PLATT (1984) WITH ENHANCEMENT FOR PARTICLES
!     SMALLER THAN 20 MICRONS.  
!
!              T Range (K)               Reff (microns)   Deff (microns)
!     -------------------------------    --------------   --------------
!
!     Tfreeze-25. < T                       92.46298         100.6
!     Tfreeze-30. < T <= Tfreeze-25.        72.35392          80.8
!     Tfreeze-35. < T <= Tfreeze-30.        85.19071          93.5
!     Tfreeze-40. < T <= Tfreeze-35.        55.65818          63.9
!     Tfreeze-45. < T <= Tfreeze-40.        35.29989          42.5
!     Tfreeze-50. < T <= Tfreeze-45.        32.89967          39.9
!     Tfreeze-55  < T <= Tfreeze-50         16.60895          21.6
!                   T <= Tfreeze-55.        15.41627          20.2
!
!        IN THIS PROGRAM, reff_ice is limited to be between 10 microns
!        and 130 microns, which is the range of validity for the Ebert
!        and Curry (1992) radiation.
!
!        IN THIS PROGRAM, size_ice (i.e. Deff) is limited to be between
!        18.6 microns and 130.2 microns, which is the range of validity 
!        for the Fu Liou JAS 1993 radiation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

REAL,     INTENT (IN), DIMENSION(:,:)    :: coszen, N_drop, k_ratio
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: ql,qi,qa,pfull,TKel
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: phalf
INTEGER,  INTENT (INOUT),DIMENSION(:,:)  :: nclds
INTEGER,  INTENT (INOUT),DIMENSION(:,:,:):: ktop,kbot
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: cldamt
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: LWP,IWP,Reff_liq,Reff_ice
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: conc_drop_org,conc_ice_org
REAL,     INTENT (INOUT),OPTIONAL,DIMENSION(:,:,:):: size_drop_org,size_ice_org


!  Internal variables
!  ------------------

INTEGER                                  :: i,j,k,IDIM,JDIM,KDIM
INTEGER                                  :: t,top_t,bot_t
INTEGER                                  :: tmp_top,tmp_bot,nlev
LOGICAL                                  :: add_cld,lhsw,sea_esf
REAL                                     :: sum_liq,sum_ice,maxcldfrac
REAL                                     :: totcld_bot,max_bot
REAL                                     :: totcld_top,max_top,tmp_val
REAL                                     :: reff_liq_local,sum_reff_liq
REAL                                     :: reff_ice_local,sum_reff_ice

!
! Code
! ----


    ! reinitialize variables
    IDIM=SIZE(ql,1)
    JDIM=SIZE(ql,2)
    KDIM=SIZE(ql,3)
    nclds(:,:)    = 0
    ktop(:,:,:)   = 0
    kbot(:,:,:)   = 0
    cldamt(:,:,:) = 0.

    !decide which type of output is necessary
    lhsw = .FALSE.
    sea_esf = .FALSE.
    if (PRESENT(LWP).or.PRESENT(IWP).or.PRESENT(Reff_liq).or. &
        PRESENT(Reff_ice)) then
        lhsw = .TRUE.
    end if
    if (PRESENT(conc_drop_org).or.PRESENT(conc_ice_org).or.  &
        PRESENT(size_drop_org).or.PRESENT(size_ice_org)) then
        sea_esf = .TRUE.
    end if
    if ((.not.lhsw).and.(.not.sea_esf)) then
        lhsw = .TRUE.
    end if

    !initialize output fields
    if (lhsw) then
         LWP(:,:,:)    = 0.
         IWP(:,:,:)    = 0.
         Reff_liq(:,:,:) = 10.
         Reff_ice(:,:,:) = 30.
    end if

    if (sea_esf) then
         conc_drop_org(:,:,:) = 0.
         conc_ice_org (:,:,:) = 0.
         size_drop_org(:,:,:) = 20.
         size_ice_org (:,:,:) = 60.
    end if
        


    !-----------  DETERMINE CLOUD AMOUNT, LWP, IWP FOR EACH CLOUD ------!

    if (overlap .eq. 1) then


         !-----------  prevent user from attempting to do maximum
         !-----------  -random overlap with new radiation code

         if (sea_esf) then

              call error_mesg  ('cloud_rad_organize in cloud_rad module',&
                                'maximum random overlap is not currently '//&
                                'available with sea_esf radiation', FATAL)

        end if

        !---- DO CONDENSING OF CLOUDS ----!

        
        !---loop over vertical levels----!
        DO i = 1, IDIM
        DO j = 1, JDIM
        
        add_cld  = .FALSE.

        DO k = 1, KDIM

                 !identify new cloud tops
                 IF ( ( (ql(i,j,k) .gt. qmin) .or. &
                        (qi(i,j,k) .gt. qmin) ) .and. &
                           (.NOT. add_cld) ) then
                       nclds(i,j) = nclds(i,j) + 1
                       add_cld = .TRUE.
                       ktop(i,j,nclds(i,j)) = k
                       sum_liq          = 0.
                       sum_ice          = 0.
                       maxcldfrac       = 0.
                       sum_reff_liq     = 0.
                       sum_reff_ice     = 0.        
                 END IF

                 !increment sums where cloud
                 IF (   (ql(i,j,k) .gt. qmin) .or. &
                        (qi(i,j,k) .gt. qmin)  ) then

                       !compute reff
                       if (ql(i,j,k) .gt. qmin) then
                          reff_liq_local = k_ratio(i,j)* 620350.49 *    &
                             (pfull(i,j,k)*ql(i,j,k)/qa(i,j,k)/Rdgas/   &
                             TKel(i,j,k)/Dens_h2o/N_drop(i,j))**(1./3.)
                       else
                          reff_liq_local = 0.
                       end if

                       if ( (k .eq. 1    .and. qa(i,j,2)      .lt. qamin) &
                          .or. &
                          (k .eq. KDIM .and. qa(i,j,KDIM-1) .lt. qamin) &
                          .or. &
                          (k .gt. 1 .and. k .lt. KDIM .and. &
                          qa(i,j,k-1) .lt. qamin .and. &
                          qa(i,j,k+1) .lt. qamin) ) then
                          reff_liq_local = 1.134 * reff_liq_local
                       end if

                       !limit reff_liq_local to values for which
                       !Slingo radiation is valid :
                       ! 4.2 microns < reff < 16.6 microns
                       reff_liq_local = MIN(16.6,reff_liq_local)
                       reff_liq_local = MAX(4.2, reff_liq_local)

                       if (qi(i,j,k) .gt. qmin) then
                          
                          if (TKel(i,j,k) .gt. Tfreeze-25.) then
                             reff_ice_local = 92.46298
                          else if (TKel(i,j,k) .gt. Tfreeze-30. .and. &
                                   TKel(i,j,k) .le. Tfreeze-25.) then
                             reff_ice_local = 72.35392
                          else if (TKel(i,j,k) .gt. Tfreeze-35. .and. &
                                   TKel(i,j,k) .le. Tfreeze-30.) then
                             reff_ice_local = 85.19071 
                          else if (TKel(i,j,k) .gt. Tfreeze-40. .and. &
                                   TKel(i,j,k) .le. Tfreeze-35.) then
                             reff_ice_local = 55.65818
                          else if (TKel(i,j,k) .gt. Tfreeze-45. .and. &
                                   TKel(i,j,k) .le. Tfreeze-40.) then
                             reff_ice_local = 35.29989
                          else if (TKel(i,j,k) .gt. Tfreeze-50. .and. &
                                   TKel(i,j,k) .le. Tfreeze-45.) then
                             reff_ice_local = 32.89967
                          else if (TKel(i,j,k) .gt. Tfreeze-55. .and. &
                                   TKel(i,j,k) .le. Tfreeze-50.) then
                             reff_ice_local = 16.60895
                          else
                             reff_ice_local = 15.41627
                          end if
                          !limit values to that for which Ebert and
                          !Curry radiation is valid :
                          !  10 microns < reff < 130 microns
                          !
                          reff_ice_local = MIN(130.,reff_ice_local)
                          reff_ice_local = MAX(10.,reff_ice_local)

                       else
                          reff_ice_local = 0.
                       end if   !end if for qi > qmin

                       !increment sums
                       sum_liq = sum_liq + ql(i,j,k)* &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav
                       sum_ice = sum_ice + qi(i,j,k)* &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav
                       maxcldfrac = MAX(maxcldfrac,qa(i,j,k))
                       sum_reff_liq  = sum_reff_liq + &
                          ( reff_liq_local * ql(i,j,k) * &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav )
                       sum_reff_ice  = sum_reff_ice + &
                          ( reff_ice_local * qi(i,j,k) * &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav )

                 END IF


                 !where the first cloud gap exists after a cloud
                 ! or bottom level is reached compute kbot, cldamt,
                 ! LWP, IWP, Reff_liq, and Reff_ice
                 IF (  ( (ql(i,j,k) .le. qmin) .and. &
                         (qi(i,j,k) .le. qmin) .and. &
                         (add_cld) ) .or. &
                         (add_cld .and. k .eq. KDIM)) then

                    !reset add_cld
                    add_cld     = .FALSE.

                    !determine kbot
                    kbot(i,j,nclds(i,j))= k-1
                    if ((ql(i,j,k) .gt. qmin) .or. &
                        (qi(i,j,k) .gt. qmin)) then
                       kbot(i,j,nclds(i,j)) = k
                    end if

                    cldamt(i,j,nclds(i,j)) = maxcldfrac
                    LWP(i,j,nclds(i,j)) = sum_liq / cldamt(i,j,nclds(i,j))
                    IWP(i,j,nclds(i,j)) = sum_ice / cldamt(i,j,nclds(i,j))
                    if (sum_liq .gt. 0.) then
                       Reff_liq(i,j,nclds(i,j)) = sum_reff_liq / sum_liq
                    end if
                    if (sum_ice .gt. 0.) then
                       Reff_ice(i,j,nclds(i,j)) = sum_reff_ice / sum_ice
                    end if

                    ! If adjust_top is T then
                    !change top and bottom indices to those that
                    !are at the most exposed to top and bottom
                    !view
                    nlev = kbot(i,j,nclds(i,j))-ktop(i,j,nclds(i,j))+1
                    if (adjust_top .and. nlev .gt. 1) then

                       !reset tmp_top,tmp_bot
                       tmp_top = ktop(i,j,nclds(i,j))
                       tmp_bot = kbot(i,j,nclds(i,j))

                       !find top and base of cloud
                       totcld_bot=0.
                       totcld_top=0.
                       max_bot=0.
                       max_top=0.
          
                       DO t = 1,nlev

                          top_t = ktop(i,j,nclds(i,j))+t-1
                          bot_t = kbot(i,j,nclds(i,j))-t+1
                          
                          tmp_val = MAX(0.,qa(i,j,top_t)-totcld_top)
                          if (tmp_val .gt. max_top) then
                             max_top = tmp_val
                             tmp_top = top_t
                          end if
                          totcld_top = totcld_top+tmp_val         
                              
                          tmp_val = MAX(0.,qa(i,j,bot_t)-totcld_bot)
                          if (tmp_val .gt. max_bot) then
                             max_bot = tmp_val
                             tmp_bot = bot_t
                          end if
                          totcld_bot = totcld_bot+tmp_val         
                               
                       END DO
                       
                       !assign tmp_top and tmp_bot to ktop and kbot
                       ktop(i,j,nclds(i,j)) = tmp_top
                       kbot(i,j,nclds(i,j)) = tmp_bot

                    end if  !for adjust_top

                 END IF  !for end of cloud

        END DO
        END DO
        END DO

    else if (overlap .eq. 2) then

           
        !---loop over vertical levels----!
        DO i = 1, IDIM
        DO j = 1, JDIM
        DO k = 1, KDIM
               
                 !where cloud exists compute ktop,kbot, cldamt and LWP and IWP
                 IF ( (ql(i,j,k) .gt. qmin) .or. &
                      (qi(i,j,k) .gt. qmin)  ) then

                    nclds(i,j) = nclds(i,j) + 1
                    ktop(i,j,nclds(i,j)) = k
                    kbot(i,j,nclds(i,j)) = k

                    if (lhsw) then
                       cldamt(i,j,nclds(i,j)) = qa(i,j,k)
                       LWP(i,j,nclds(i,j))    = ql(i,j,k)*    &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav/ &
                          cldamt(i,j,nclds(i,j))
                       IWP(i,j,nclds(i,j))    = qi(i,j,k)*    &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Grav/ &
                          cldamt(i,j,nclds(i,j))
                    end if  !lhsw if

                    if (sea_esf) then
                       cldamt(i,j,k) = qa(i,j,k)
                       !Note units are in g/m3!
                       conc_drop_org(i,j,k) = 1000.*ql(i,j,k)*                &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Rdgas/TKel(i,j,k)/    &
                          log(phalf(i,j,k+1)/MAX(phalf(i,j,k),pfull(i,j,1)))/ &
                          cldamt(i,j,k)
                       conc_ice_org (i,j,k) = 1000.*qi(i,j,k)*                &
                          (phalf(i,j,k+1)-phalf(i,j,k))/Rdgas/TKel(i,j,k)/    &
                          log(phalf(i,j,k+1)/MAX(phalf(i,j,k),pfull(i,j,1)))/ &
                          cldamt(i,j,k)
                    end if  !sea_esf if

                    !compute reff_liquid
                    if (ql(i,j,k) .gt. qmin) then

                       reff_liq_local = k_ratio(i,j)* 620350.49 *    &
                          (pfull(i,j,k)*ql(i,j,k)/qa(i,j,k)/Rdgas/   &
                          TKel(i,j,k)/Dens_h2o/N_drop(i,j))**(1./3.)
                       
                       if ( (k .eq. 1    .and. qa(i,j,2)      .lt. qamin) &
                            .or. &
                            (k .eq. KDIM .and. qa(i,j,KDIM-1) .lt. qamin) &
                            .or. &
                            (k .gt. 1 .and. k .lt. KDIM .and. &
                            qa(i,j,k-1) .lt. qamin .and. &
                            qa(i,j,k+1) .lt. qamin) ) then
                          reff_liq_local = 1.134 * reff_liq_local
                       end if

                       !limit reff_liq_local to values for which
                       !Slingo radiation is valid :
                       ! 4.2 microns < reff < 16.6 microns
                       reff_liq_local = MIN(16.6,reff_liq_local)
                       reff_liq_local = MAX(4.2, reff_liq_local)

                       if (lhsw) Reff_liq(i,j,nclds(i,j)) =  reff_liq_local
                       if (sea_esf) size_drop_org(i,j,k) = 2. * reff_liq_local

                    end if  !ql calculation

                    !compute reff_ice
                    if (qi(i,j,k) .gt. qmin) then
                          
                       if (lhsw) then
                       
                          if (TKel(i,j,k) .gt. Tfreeze-25.) then
                              reff_ice_local = 92.46298
                          else if (TKel(i,j,k) .gt. Tfreeze-30. .and. &
                                   TKel(i,j,k) .le. Tfreeze-25.) then
                              reff_ice_local = 72.35392
                          else if (TKel(i,j,k) .gt. Tfreeze-35. .and. &
                                   TKel(i,j,k) .le. Tfreeze-30.) then
                              reff_ice_local = 85.19071 
                          else if (TKel(i,j,k) .gt. Tfreeze-40. .and. &
                                   TKel(i,j,k) .le. Tfreeze-35.) then
                              reff_ice_local = 55.65818
                          else if (TKel(i,j,k) .gt. Tfreeze-45. .and. &
                                   TKel(i,j,k) .le. Tfreeze-40.) then
                              reff_ice_local = 35.29989
                          else if (TKel(i,j,k) .gt. Tfreeze-50. .and. &
                                   TKel(i,j,k) .le. Tfreeze-45.) then
                              reff_ice_local = 32.89967
                          else if (TKel(i,j,k) .gt. Tfreeze-55. .and. &
                                   TKel(i,j,k) .le. Tfreeze-50.) then
                              reff_ice_local = 16.60895
                          else
                              reff_ice_local = 15.41627
                          end if

                          !limit values to that for which Ebert and
                          !Curry radiation is valid :
                          !  10 microns < reff < 130 microns
                          !
                          reff_ice_local = MIN(130.,reff_ice_local)
                          Reff_ice(i,j,nclds(i,j)) = MAX(10.,reff_ice_local)                  
    
                       end if  !end of lhsw if

                       if (sea_esf) then

                          if (TKel(i,j,k) .gt. Tfreeze-25.) then
                              reff_ice_local = 100.6
                          else if (TKel(i,j,k) .gt. Tfreeze-30. .and. &
                                   TKel(i,j,k) .le. Tfreeze-25.) then
                              reff_ice_local = 80.8
                          else if (TKel(i,j,k) .gt. Tfreeze-35. .and. &
                                   TKel(i,j,k) .le. Tfreeze-30.) then
                              reff_ice_local = 93.5 
                          else if (TKel(i,j,k) .gt. Tfreeze-40. .and. &
                                   TKel(i,j,k) .le. Tfreeze-35.) then
                              reff_ice_local = 63.9
                          else if (TKel(i,j,k) .gt. Tfreeze-45. .and. &
                                   TKel(i,j,k) .le. Tfreeze-40.) then
                              reff_ice_local = 42.5
                          else if (TKel(i,j,k) .gt. Tfreeze-50. .and. &
                                   TKel(i,j,k) .le. Tfreeze-45.) then
                              reff_ice_local = 39.9
                          else if (TKel(i,j,k) .gt. Tfreeze-55. .and. &
                                   TKel(i,j,k) .le. Tfreeze-50.) then
                              reff_ice_local = 21.6
                          else
                              reff_ice_local = 20.2
                          end if

                          !the ice crystal effective size can          
                          !only be 18.6 <= D^sub^e <= 130.2 microns. 
                          !for Fu Liou JAS 1993 code                    
                          reff_ice_local = MIN(130.2,reff_ice_local)
                          size_ice_org(i,j,k) = &
                                           MAX(18.6,reff_ice_local)

                       end if  !end of sea_esf if
                             
                    end if !qi loop                       
                          
                 END IF  !cloud exist if

           END DO
           END DO
           END DO



    end if   !overlap = 2  if

    
END SUBROUTINE CLOUD_ORGANIZE


!########################################################################
!########################################################################

SUBROUTINE CLOUD_RAD(tau,w0,gg,coszen,r_uv,r_nir,ab_uv,ab_nir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine calculates the following radiative properties
!      for each cloud:
!
!               1. r_uv : cloud reflectance in uv band
!               2. r_nir: cloud reflectance in nir band
!               3. ab_uv: cloud absorption in uv band
!               4. ab_nir:cloud absorption in nir band
!               
!
!      These quantities are computed by dividing the shortwave
!      spectrum into 4 bands and then computing the reflectance
!      and absorption for each band individually and then setting
!      the uv reflectance and absorption equal to that of band
!      1 and the nir reflectance and absorption equal to the
!      spectrum weighted results of bands 2,3,and 4.  The limits
!      of bands are described in CLOUD_OPTICAL_PROPERTIES.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	VARIABLES
!
!       ------
!	INPUT:
!       ------
!
!       tau          optical depth in 4 bands (dimensionless)
!       w0           single scattering albedo in 4 bands (dimensionless)
!       gg           asymmetry parameter in 4 bands (dimensionless)
!       coszen       cosine of the zenith angle
!
!       ------
!	INPUT/OUTPUT:
!       ------
!
!       r_uv         cloud reflectance in uv band
!       r_nir        cloud reflectance in nir band
!       ab_uv        cloud absorption in uv band
!       ab_nir       cloud absorption in nir band
!
!       -------------------
!	INTERNAL VARIABLES:
!       -------------------
!
!      direct       logical variable for each cloud indicating whether
!                       or not to use the direct beam solution for the
!                       delta-eddington radiation or the diffuse beam
!                       radiation solution.
!      tau_local    optical depth for the band being solved
!      w0_local     single scattering albedo for the band being solved
!      g_local      asymmetry parameter for the band being solved
!      coszen_3d    3d version of coszen
!      I            looping variable
!      iband        looping variables over band number
!      taucum       cumulative sum of visible optical depth
!      g_prime      scaled g
!      w0_prime     scaled w0
!      tau_prime    scaled tau
!      crit         variable equal to 1./(4 - 3g')
!      AL           variable equal to sqrt(3*(1-w0')*(1-w0'*g'))
!      ALPHV        temporary work variable
!      GAMV         temporary work variable
!      T1V          exp( -1.*AL * tau')
!      trans_dir    direct radiation beam transmittance
!      U            1.5 * (1. - w0'*g')/AL
!      r_diffus     diffuse beam reflection
!      trans_diffus diffuse beam transmission
!
!      r            cloud reflectance for each cloud in each band
!      ab           cloud absorption for each cloud in each band
!      r_dir_uv       direct beam reflection for uv band
!      r_dir_nir      direct beam reflection for nir band
!      trans_dir_uv   direct beam transmission for uv band
!      trans_dir_nir  direct beam transmission for uv band
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

REAL,     INTENT (IN), DIMENSION(:,:,:,:):: tau,w0,gg
REAL,     INTENT (IN), DIMENSION(:,:)    :: coszen
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: r_uv,r_nir,ab_uv,ab_nir

!  Internal variables
!  ------------------

INTEGER                                  :: I,iband,J,K
REAL,    DIMENSION(SIZE(tau,1),SIZE(tau,2)) :: taucum
LOGICAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: direct
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: coszen_3d
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: tau_local
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: w0_local,g_local
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: g_prime,w0_prime
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: tau_prime,crit,AL
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: ALPHV,GAMV,T1V,U
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: r_diffus
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: trans_diffus
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3)) :: r_dir,trans_dir
REAL, DIMENSION(SIZE(tau,1),SIZE(tau,2),SIZE(tau,3),SIZE(tau,4)) :: r,ab

!
! Code
! ----

        ! reinitialize variables
        r_uv(:,:,:) = 0.
        r_nir(:,:,:)= 0.
        ab_uv(:,:,:)= 0.
        ab_nir(:,:,:)=0.

        !create 3d zenith angle
        DO I = 1, SIZE(tau,3)
               coszen_3d(:,:,I)=coszen(:,:)
        END DO
        WHERE (coszen_3d(:,:,:) .lt. 1.E-06)
                coszen_3d(:,:,:) = 1.E-06
        END WHERE

        
        !------------------------------------------------------------------
        !do logical variable to determine where total cloud optical depth
        !at uv wavelengths exceeds taucrit
        IF (.NOT. l2strem) THEN

                !---- initialize taucum and direct
                taucum(:,:)=0.
                direct(:,:,:)=.TRUE.

                DO I = 1, SIZE(tau,3)

                      !find if taucum to levels above has exceeded taucrit
                      WHERE (taucum(:,:) .gt. taucrit)
                            direct(:,:,I)=.FALSE.
                      END WHERE

                      !increment cumulative tau
                      taucum(:,:)=taucum(:,:)+tau(:,:,I,1)

                END DO
        END IF

    !----------------- LOOP OVER BAND -----------------------------!

    DO iband = 1, SIZE(tau,4)

        !-----------------------------------------------------------
        !  assign w0, g, tau to the value appropriate for the band

        w0_local(:,:,:) = w0(:,:,:,iband)
        tau_local(:,:,:)= tau(:,:,:,iband)
        g_local(:,:,:) =  gg(:,:,:,iband)

        !-------------------------------------------------------------------
        ! for delta-Eddington scaled ('prime') g, w0, tau where:
        !
        !               g' = g / (1 + g)
        !              w0' = (1 - g*g) * w0 / (1 - w*g*g)
        !             tau' = (1 - w*g*g) * tau
        !

                 tau_prime(:,:,:) = 1. - &
                      (w0_local(:,:,:)*g_local(:,:,:)*g_local(:,:,:))
                 w0_prime(:,:,:) = w0_local(:,:,:) * &
                   (1. - (g_local(:,:,:)*g_local(:,:,:)))/tau_prime(:,:,:)
                 tau_prime(:,:,:) = tau_prime(:,:,:) * tau_local(:,:,:)
                 g_prime(:,:,:) = g_local(:,:,:) / (1. + g_local(:,:,:))

        !-------------------------------------------------------------------
        ! create other variables
        !
        !        crit = 1./(4 - 3g')
        !
        !      and where w0' < crit set w0' = crit
        !
        !        AL = sqrt( 3. * (1. - w0') * (1. - w0'*g') )
        !

                 crit(:,:,:) = 1./(4.- 3.*g_prime(:,:,:))

                 WHERE (w0_prime(:,:,:) .lt. crit(:,:,:) )
                           w0_prime(:,:,:) = crit(:,:,:)
                 END WHERE

                 AL(:,:,:) =  ( 3. * (1. - w0_prime(:,:,:) ) &
                    * (1. - (w0_prime(:,:,:)*g_prime(:,:,:)))  )**0.5

                 !set up a minimum to AL
                 WHERE (AL(:,:,:) .lt. 1.E-06)
                        AL(:,:,:) = 1.E-06
                 END WHERE


        !-------------------------------------------------------------------
        ! simplifications if not two stream
        !
        !        ALPHV = 0.75*w0'*coszen*(1.+g'(1.-w0'))/
        !                          (1.-(AL*coszen)**2.)
        !        GAMV = 0.5*w0'*(3.*g'*(1.-w0')*coszen*coszen + 1.)/
        !                          (1.-(AL*coszen)**2.)
        !

        IF (.NOT. l2strem) THEN


                ALPHV(:,:,:) = 0.75 * w0_prime(:,:,:)*coszen_3d(:,:,:) * &
                 (1. + (g_prime(:,:,:)*(1. - w0_prime(:,:,:)))) / &
                 (1. - (AL(:,:,:)*coszen_3d(:,:,:))**2.0)

                GAMV(:,:,:) =  0.50 * w0_prime(:,:,:) * &
                (  (3.* g_prime(:,:,:) * (1. - w0_prime(:,:,:)) * &
                    coszen_3d(:,:,:) * coszen_3d(:,:,:)) + 1. ) / &
                 (1. - (AL(:,:,:)*coszen_3d(:,:,:))**2.0)

        END IF


        !-------------------------------------------------------------------
        ! calculate T1V
        !
        !    T1V = exp (-1* AL * tau' )


                  T1V(:,:,:) = exp( -1.*AL(:,:,:) * tau_prime(:,:,:) )


        !-------------------------------------------------------------------
        !calculate diffuse beam reflection and transmission
        !

        !first calculate U  = 1.5 * (1. - w0'*g')/AL
        U(:,:,:) = 1.5 *(1. - w0_prime(:,:,:)*g_prime(:,:,:))/AL(:,:,:)

        !initialize variables
        r_diffus(:,:,:)= 0.
        trans_diffus(:,:,:) = 1.



        trans_diffus(:,:,:) = 4. * U(:,:,:) * T1V(:,:,:) / &
            ( ( (U(:,:,:)+1.) * (U(:,:,:)+1.)  ) - &
              ( (U(:,:,:)-1.) * (U(:,:,:)-1.) * &
                   T1V(:,:,:) *   T1V(:,:,:)   )    )

        r_diffus(:,:,:) =     ((U(:,:,:)*U(:,:,:))-1.) * &
                   ( 1. -   (T1V(:,:,:)*T1V(:,:,:)) ) / &
             ( ( (U(:,:,:)+1.) * (U(:,:,:)+1.)  ) - &
               ( (U(:,:,:)-1.) * (U(:,:,:)-1.) * &
                    T1V(:,:,:) *   T1V(:,:,:)   )    )



        !-------------------------------------------------------------------
        ! calculate direct bean transmission
        !
        !
        IF (.NOT. l2strem) THEN


            !initialize variables
            trans_dir(:,:,:) = 1.
            r_dir(:,:,:) = 0.

            r_dir(:,:,:) = ( (ALPHV(:,:,:) - GAMV(:,:,:)) * &
               exp(-1.*tau_prime(:,:,:)/coszen_3d(:,:,:)) * &
               trans_diffus(:,:,:) ) +  &
              ( (ALPHV(:,:,:) + GAMV(:,:,:)) * &
              r_diffus(:,:,:) )  -  (ALPHV(:,:,:) - GAMV(:,:,:))

            trans_dir(:,:,:) = &
              ( (ALPHV(:,:,:)+GAMV(:,:,:))*trans_diffus(:,:,:) ) + &
              ( exp(-1.*tau_prime(:,:,:)/coszen_3d(:,:,:)) * &
              ( ( (ALPHV(:,:,:)-GAMV(:,:,:))*r_diffus(:,:,:) ) - &
                (ALPHV(:,:,:)+GAMV(:,:,:)) + 1. )   )

        END IF


        !-------------------------------------------------------------------
        ! patch together final solution
        !
        !


        IF (l2strem) THEN

             !two-stream solution
             r(:,:,:,iband) = r_diffus(:,:,:)
             ab(:,:,:,iband) = 1. - trans_diffus(:,:,:) - r_diffus(:,:,:)

        ELSE

             !delta-Eddington solution
             WHERE (.not. direct)

                   r(:,:,:,iband) = r_diffus(:,:,:)
                   ab(:,:,:,iband) = 1. - trans_diffus(:,:,:) &
                                     - r_diffus(:,:,:)


             END WHERE

             WHERE (direct)

                   r(:,:,:,iband) = r_dir(:,:,:)
                   ab(:,:,:,iband) = 1. - trans_dir(:,:,:) &
                                     - r_dir(:,:,:)

             END WHERE

        END IF

    !----------------- END LOOP OVER BAND -----------------------------!

    END DO


    !----------------- CREATE SUM OVER BAND ---------------------------!

    r_uv(:,:,:) = r(:,:,:,1)
    ab_uv(:,:,:) = ab(:,:,:,1)

    r_nir(:,:,:) =  (  0.326158 * r(:,:,:,2) + &
                       0.180608 * r(:,:,:,3) + &
                       0.033474 * r(:,:,:,4) ) / 0.540240

    ab_nir(:,:,:) =  (  0.326158 * ab(:,:,:,2) + &
                        0.180608 * ab(:,:,:,3) + &
                        0.033474 * ab(:,:,:,4) ) / 0.540240


        !-------------------------------------------------------------------
        ! guarantee that clouds for tau = 0. have the properties
        ! of no cloud

        
        WHERE(tau(:,:,:,1) .le. 0.)
             r_uv(:,:,:) = 0.
             ab_uv(:,:,:)= 0.                       
        END WHERE
        WHERE((tau(:,:,:,2)+tau(:,:,:,3)+tau(:,:,:,4)) .le. 0.)
             r_nir(:,:,:)= 0.
             ab_nir(:,:,:)=0.
        END WHERE       

        !-------------------------------------------------------------------
        ! guarantee that for coszen lt. or equal to zero that solar
        ! reflectances and absorptances are equal to zero.
        DO I = 1, SIZE(tau,3)
               WHERE (coszen(:,:) .lt. 1.E-06)
                    r_uv(:,:,I) = 0. 
                    ab_uv(:,:,I) = 0.
                    r_nir(:,:,I) = 0.
                    ab_nir(:,:,I) = 0.
               END WHERE
        END DO
        
        !-------------------------------------------------------------------
        ! guarantee that each cloud has some transmission by reducing
        ! the actual cloud reflectance in uv and nir band
        ! this break is necessary to avoid the rest of the
        ! radiation code from breaking up.
        !

        WHERE ( (1. - r_uv(:,:,:) - ab_uv(:,:,:)) .lt. 0.01)
                      r_uv(:,:,:) = r_uv(:,:,:) - 0.01
        END WHERE
        WHERE ( (1. - r_nir(:,:,:) - ab_nir(:,:,:)) .lt. 0.01)
                      r_nir(:,:,:) = r_nir(:,:,:) - 0.01
        END WHERE

        !-------------------------------------------------------------------
        ! guarantee that cloud reflectance and absorption are greater than
        ! or equal to zero

        WHERE (r_uv(:,:,:) .lt. 0.)
               r_uv(:,:,:) = 0.
        END WHERE
        WHERE (r_nir(:,:,:) .lt. 0.)
               r_nir(:,:,:) = 0.
        END WHERE
        WHERE (ab_uv(:,:,:) .lt. 0.)
               ab_uv(:,:,:) = 0.
        END WHERE
        WHERE (ab_nir(:,:,:) .lt. 0.)
               ab_nir(:,:,:) = 0.
        END WHERE


END SUBROUTINE CLOUD_RAD

!########################################################################
!########################################################################


SUBROUTINE CLOUD_OPTICAL_PROPERTIES(LWP,IWP,Reff_liq,Reff_ice,&
                        tau,w0,gg,em_lw,tau_ice_diag)

                              

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine calculates the following optical properties
!      for each cloud:
!
!               1. tau   :optical depth in each band
!               2. w0    :single scattering albedo for each band
!               3. gg     :asymmetry parameter for each band
!               4. em_lw    :longwave cloud emmissivity
!
!   The formulas for optical depth come from Slingo (1989) for liquid
!   clouds and from Ebert and Curry (1992) for ice clouds.
!
!   Slingo (1989) is at J. Atmos. Sci., vol. 46, pp. 1419-1427
!   Ebert and Curry (1992) is at J. Geophys. Res., vol. 97, pp. 3831-3836
!
!                    IMPORTANT!!!
!
!    NOTE WE ARE CHEATING HERE BECAUSE WE ARE FORCING THE FIVE BAND
!    MODEL OF EBERT AND CURRY INTO THE FOUR BAND MODEL OF SLINGO
!
!    THIS IS DONE BY COMBINING BANDS 3 and 4 OF EBERT AND CURRY TOGETHER
!
!   EVEN SO THE EXACT BAND LIMITS DO NOT MATCH.  FOR COMPLETENESS
!   HERE ARE THE BAND LIMITS IN MICRONS
!
!            BAND               SLINGO                 EBERT AND CURRY
!
!             1               0.25-0.69                0.25 - 0.7
!             2               0.69-1.19                0.7 - 1.3
!             3               1.19-2.38                1.3 - 2.5
!             4               2.38-4.00                2.5 - 3.5
!
!
!   The mixed phase optical properties are based upon equation 14
!   of Rockel et al. 1991, Contributions to Atmospheric Physics,
!   volume 64, pp.1-12.   These equations are:
!
!   (1)    tau = tau_liq + tau_ice
!
!   (2)    w0  =   ( w0_liq * tau_liq  +  w0_ice * tau_ice ) /
!                  (          tau_liq  +           tau_ice )
!
!   (3)     g  = ( g_liq * w0_liq * tau_liq +  g_ice * w0_ice * tau_ice ) /
!                (         w0_liq * tau_liq +          w0_ice * tau_ice )
!
!
!   (4) transmivvity_lw =   transmissivity_lw_ice * transmissivity_lw_liq
!
!    The last equation can be rewritten as:
!
!   (5)  em_lw =  em_lw_liq + em_lw_ice -  (em_lw_liq * em_lw_ice )
!
!   Which is what is solved here.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	VARIABLES
!
!       ------
!	INPUT:
!       ------
!
!       LWP          cloud liquid water path (kg condensate per square meter)
!       IWP          cloud ice path (kg condensate per square meter)
!       Reff_liq     effective radius for liquid clouds (microns)
!       Reff_ice     effective particle size for ice clouds (microns)
!
!       ------
!	INPUT/OUTPUT:
!       ------
!
!      tau          optical depth in each band
!      w0           single scattering albedo for each band
!      gg           asymmetry parameter for each band
!      em_lw        longwave cloud emmissivity
!
!       ---------------------
!       OPTIONAL INPUT/OUTPUT
!       ---------------------
!
!       tau_ice_diag    optical depth in each band
!
!
!       -------------------
!	INTERNAL VARIABLES:
!       -------------------
!
!       tau_liq   optical depth            at each band for cloud liquid
!       tau_ice   optical depth            at each band for cloud ice
!       w0_liq    single scattering albedo at each band for cloud liquid
!       w0_ice    single scattering albedo at each band for cloud ice
!       g_liq     asymmetry parameter      at each band for cloud liquid
!       g_ice     asymmetry parameter      at each band for cloud ice
!       k_liq        liquid cloud mass absorption coefficient for longwave
!                         portion of the spectrum (meters**2./kg condensate)
!       k_ice           ice cloud mass absorption coefficient for longwave
!                         portion of the spectrum (meters**2./kg condensate)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------


REAL,     INTENT (IN)   ,DIMENSION(:,:,:)   :: LWP,IWP,Reff_liq,Reff_ice
REAL,     INTENT (INOUT),DIMENSION(:,:,:,:) :: tau,w0,gg
REAL,     INTENT (INOUT),DIMENSION(:,:,:)   :: em_lw
REAL,     INTENT (INOUT), OPTIONAL, DIMENSION(:,:,:,:) :: tau_ice_diag
        
!  Internal variables
!  ------------------
REAL, DIMENSION(SIZE(LWP,1),SIZE(LWP,2),SIZE(LWP,3),4) :: tau_liq,tau_ice
REAL, DIMENSION(SIZE(LWP,1),SIZE(LWP,2),SIZE(LWP,3),4) :: w0_liq,w0_ice
REAL, DIMENSION(SIZE(LWP,1),SIZE(LWP,2),SIZE(LWP,3),4) :: g_liq,g_ice
REAL, DIMENSION(SIZE(LWP,1),SIZE(LWP,2),SIZE(LWP,3))   :: k_liq,k_ice

!
! Code
! ----
        
        ! reinitialize output variables to default values
        ! (not usually used)
        tau(:,:,:,:)=0.
        gg(:,:,:,:) = 0.85
        w0(:,:,:,:) = 0.95
        em_lw(:,:,:) = 0.
        
        ! reinitialize internal variables (not usually used)
        w0_liq(:,:,:,:) = 0.95
        w0_ice(:,:,:,:) = 0.95
        g_liq(:,:,:,:)  = 0.85
        g_ice(:,:,:,:)  = 0.85
        tau_liq(:,:,:,:)= 0.
        tau_ice(:,:,:,:)= 0.



   !---------------   COMPUTE OPTICAL DEPTH ---------------------------!

        ! compute uv cloud optical depths due to liquid
        ! and ice phase separately

        tau_liq(:,:,:,1) = LWP(:,:,:) * 1000. * &
                           (0.02817 + (1.305/Reff_liq(:,:,:)))
        tau_liq(:,:,:,2) = LWP(:,:,:) * 1000. * &
                           (0.02682 + (1.346/Reff_liq(:,:,:)))
        tau_liq(:,:,:,3) = LWP(:,:,:) * 1000. * &
                           (0.02264 + (1.454/Reff_liq(:,:,:)))
        tau_liq(:,:,:,4) = LWP(:,:,:) * 1000. * &
                           (0.01281 + (1.641/Reff_liq(:,:,:)))
        
        tau_ice(:,:,:,1) = IWP(:,:,:) * 1000. * &
                           (0.003448 + (2.431/Reff_ice(:,:,:)))
        tau_ice(:,:,:,2) = tau_ice(:,:,:,1)
        tau_ice(:,:,:,3) = tau_ice(:,:,:,1)
        tau_ice(:,:,:,4) = tau_ice(:,:,:,1)
        

        ! compute total cloud optical depth
        tau(:,:,:,:) = tau_liq(:,:,:,:) + tau_ice(:,:,:,:)
        

   !---------------   COMPUTE SINGLE SCATTERING ALBEDO ----------------!

        w0_liq(:,:,:,1) =  5.62E-08   - 1.63E-07*Reff_liq(:,:,:)
        w0_liq(:,:,:,2) =  6.94E-06   - 2.35E-05*Reff_liq(:,:,:)
        w0_liq(:,:,:,3) = -4.64E-04   - 1.24E-03*Reff_liq(:,:,:)
        w0_liq(:,:,:,4) = -2.01E-01   - 7.56E-03*Reff_liq(:,:,:)
        w0_liq(:,:,:,:) = w0_liq(:,:,:,:) + 1.

        w0_ice(:,:,:,1) = -1.00E-05
        w0_ice(:,:,:,2) = -1.10E-04   - 1.41E-05*Reff_ice(:,:,:)
        w0_ice(:,:,:,3) = -1.86E-02   - 8.33E-04*Reff_ice(:,:,:)
        w0_ice(:,:,:,4) = -4.67E-01   - 2.05E-05*Reff_ice(:,:,:)
        w0_ice(:,:,:,:) = w0_ice(:,:,:,:) + 1.


        ! compute total single scattering albedo
        WHERE (tau(:,:,:,:) .gt. 0.)
               w0(:,:,:,:) = ( w0_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                               w0_ice(:,:,:,:) * tau_ice(:,:,:,:) )  /&
                             tau(:,:,:,:)
        END WHERE
        
   !---------------   COMPUTE ASYMMETRY PARAMETER --------------------!


       g_liq(:,:,:,1) = 0.829 + 2.482E-03*Reff_liq(:,:,:)
       g_liq(:,:,:,2) = 0.794 + 4.226E-03*Reff_liq(:,:,:)
       g_liq(:,:,:,3) = 0.754 + 6.560E-03*Reff_liq(:,:,:)
       g_liq(:,:,:,4) = 0.826 + 4.353E-03*Reff_liq(:,:,:)

       g_ice(:,:,:,1) = 0.7661+ 5.851E-04*Reff_ice(:,:,:)
       g_ice(:,:,:,2) = 0.7730+ 5.665E-04*Reff_ice(:,:,:)
       g_ice(:,:,:,3) = 0.7940+ 7.267E-04*Reff_ice(:,:,:)
       g_ice(:,:,:,4) = 0.9595+ 1.076E-04*Reff_ice(:,:,:)

        ! compute  asymmetry parameter
        WHERE (tau(:,:,:,:) .gt. 0. )
              gg(:,:,:,:) = ( &
                 w0_liq(:,:,:,:) * g_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                 w0_ice(:,:,:,:) * g_ice(:,:,:,:) * tau_ice(:,:,:,:) ) &
                       /          (w0_liq(:,:,:,:) * tau_liq(:,:,:,:) + &
                                   w0_ice(:,:,:,:) * tau_ice(:,:,:,:) )
        END WHERE

        
   !---------------   COMPUTE LONGWAVE EMMISSIVITY --------------------!


        k_liq(:,:,:) = 140.
        k_ice(:,:,:) = 4.83591 + 1758.511/Reff_ice(:,:,:)
        
        ! compute combined emmisivity
        em_lw(:,:,:) =  1. - exp( -1. * ( k_liq(:,:,:) * LWP(:,:,:) + &
                                          k_ice(:,:,:) * IWP(:,:,:) ) )

        
   !--------------    RANGE LIMIT QUANTITIES --------------------------!

        WHERE (tau(:,:,:,:) .lt. taumin)
               tau(:,:,:,:) = taumin
        END WHERE

   !----- ---------    EVALUATE TAU_ICE_DIAG (OPTIONAL) ---------------!

        
        if (present(tau_ice_diag)) then
               tau_ice_diag(:,:,:,:) = 0.
               tau_ice_diag(:,:,:,:) = tau_ice(:,:,:,:)
        end if
 

END SUBROUTINE CLOUD_OPTICAL_PROPERTIES

!########################################################################
!########################################################################

SUBROUTINE ISCCP_CLOUDTYPES(pfull,phalf,qv,at,skt,cc,dtau_s,dem_s,&
                            fq_isccp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine calculates the fraction of each model 
!      grid box covered by each of the 49 ISCCP D level cloud 
!      types (i.e. stratified by optical depth and cloud top
!      pressure) by accounting for model overlap. For further 
!      explanation see Klein and Jakob, Monthly Weather Review,
!      (2000), vol x, pp. .
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	VARIABLES
!
!       ------
!	INPUT:
!       ------
!  
!       pfull            pressure of full model levels (Pascals)
!                        pfull(1)    is    top level of model
!                        pfull(nlev) is bottom level of model
!       phalf            pressure of half model levels (Pascals)
!                        phalf(1)    is    top       of model
!                        phalf(nlev+1) is the surface pressure
!       cc               cloud cover in each model level (fraction)
!                        this includes convective clouds if any
!  
!                        NOTE:  This is the HORIZONTAL area of each
!                               grid box covered by clouds
!                         
!       dtau_s           mean 0.67 micron optical depth of stratiform
!                        clouds in each model level
!
!                        NOTE:  this the cloud optical depth of only the
!                               cloudy part of the grid box, it is not 
!                               weighted with the 0 cloud optical depth 
!                               of the clear part of the grid box
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
!       The following input variables are used only if 
!       adjust_isccp_top_height = .TRUE.
!
!       qv               water vapor specific humidity on model levels
!                        (kg vapor / kg air)
!       at               temperature in each model level (K)
!       dem_s            10.5 micron longwave emissivity of stratiform
!                        clouds in each model level.  Same note applies 
!                        as in dtau.
!       skt              skin Temperature (K)
!
!
!       -------------------
!	OUTPUT
!       -------------------
!
!       fq_isccp        matrix of fractional area covered by cloud
!                       types of a given optical depth and cloud
!                       top pressure range.  The matrix is 7x7 for
!                       7 cloud optical depths and 7 cloud top 
!                       pressure ranges
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

REAL,     INTENT (IN), DIMENSION(:)      :: pfull,phalf,at,qv
REAL,     INTENT (IN)                    :: skt
REAL,     INTENT (INOUT), DIMENSION(:)   :: cc,dtau_s,dem_s
!REAL,     INTENT (INOUT), DIMENSION(:)   :: conv,dtau_c,dem_c

REAL,     INTENT (OUT), DIMENSION(7,7)   :: fq_isccp

!
!
!  Internal variables
!  ------------------
!

  
REAL                        :: ptrop,attrop,atmax,atmin,btcmin,transmax, &
                               fluxtop,fluxtopclrsky,trans,transclrsky, &
                               wtmair, wtmh20, Navo, grav, pstd, t0, &
                               press, dpress, atmden, rvh20, wk, rhoave, &
                               rh20s, rfrgn, tmpexp,tauwv, emcld, &
                               taumintmp
INTEGER                     :: ilev,j,ibox,itrop,ipres,itau,nmatch,nlev

REAL,    DIMENSION(size(qv,1))      :: demwv
INTEGER, DIMENSION(size(qv,1)-1)    :: match
INTEGER, DIMENSION(ncol)            :: levmatch
REAL,    DIMENSION(ncol)            :: tau,tb,ptop

REAL,    DIMENSION(ncol,size(qv,1))  :: frac_out
REAL,    DIMENSION(ncol,size(qv,1)+1):: tca
REAL,    DIMENSION(ncol)             :: threshold,maxosc,boxpos
REAL,    DIMENSION(ncol)             :: threshold_min
REAL                                 :: bb,dem
INTEGER                              :: seed

!convective cloud stuff
!REAL,    DIMENSION(ncol,size(qv,1))  :: cca
!REAL,    DIMENSION(ncol)             :: maxocc

!
! Code
! ----



      nlev = size(qv,1) 
     
      
!     ---------------------------------------------------!
!     find tropopause pressure and temperature min and max
!
!     The tropopause pressure is used only in the assignment
!     of cloud top pressure in the case where the infrared bright-
!     ness temperature is calculated (adjust_isccp_top_height = 
!     .TRUE.). Cloud swho infrared brightness temperatures are 
!     less than the tropopause temperature are assigned to a 
!     cloud top pressure of the tropopause pressure (see 
!     ISCCP documentation)
!

      if (adjust_isccp_top_height) then

      ptrop=5000.
      atmin = 400.
      atmax = 0.
      do  ilev=1,nlev-1
           if ((pfull(ilev)/phalf(nlev+1)) .lt. 0.4 .and. &
                at(ilev) .gt. at(ilev+1)) then
                ptrop = pfull(ilev+1)
                attrop = at(ilev+1)
                itrop=ilev+1
           end if
           if (at(ilev) .gt. atmax) atmax=at(ilev)
           if (at(ilev) .lt. atmin) atmin=at(ilev)
      end do

      end if

!     -----------------------------------------------------!

!     ---------------------------------------------------!
!     find unpermittable data.....
!
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


!     ---------------------------------------------------!
!     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
!

      !initialize variables
      
      !assign 2d tca array using 1d input array cc
      !assign 2d cca array using 1d input array conv
      tca(:,1) = 0.
      do ilev = 1,nlev
             tca(:,ilev+1) = cc(ilev)
             !cca(:,ilev  ) = conv(ilev)
      enddo
            
      seed = ncol
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
                   !maxocc(ibox)=(boxpos(ibox).le.cca(ibox,ilev))

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
                               min(tca(ibox,ilev),tca(ibox,ilev+1))
                              !max(cca(ibox,ilev), &
                              !min(tca(ibox,ilev),tca(ibox,ilev+1)))
                         if (threshold(ibox).lt.  &
                               min(tca(ibox,ilev),tca(ibox,ilev+1))) then
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
           WHERE(tca(:,ilev+1).gt.threshold(:))
                 frac_out(:,ilev)=1
           ENDWHERE

           !Code to partition boxes into startiform and convective parts
           !goes here
           !
           !DO ibox=1,ncol
           !    frac_out(ibox,ilev)=
           !    ! = the same IF NOT threshold le cca(ibox) 
           !         frac_out(ibox,ilev) + 
           !     ! = 2 IF threshold le cca(ibox)
           !        (threshold(ibox).le.cca(ibox,ilev)) 
           !      * ( 2 - frac_out(ibox,ilev) ) 
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
                 !(frac_out(:,ilev).eq.1)*dtau_s(ilev)+
                 !(frac_out(:,ilev).eq.2)*dtau_c(ilev)
      end do

!
!     ---------------------------------------------------!



!     
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done if adjust_isccp_top_height = T
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     fluxtopclrsky is the 10.5 micron radiance at the top of
!              the atmosphere under clear skies
!     trans is the transmissivity from the top of level to
!              the top of the atmosphere
!     transclrsky is the clear sky transmissivity from the top
!              of level to the top of the atmosphere


      if (adjust_isccp_top_height) then




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
                           !(1. - &
                           !((frac_out(ibox,ilev).eq.1)*dem_s(ilev) &
                           !+(frac_out(ibox,ilev).eq.2)*dem_c(ilev)))

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
            !Note that this is discussed on pages 85-87 of 
            !the ISCCP D level documentation (Rossow et al. 1996)
           
            
            !compute minimum brightness temperature and optical depth
            btcmin = 1. /  ( exp(1307.27/(attrop-5.)) - 1. ) 
            transmax = (fluxtop-btcmin)/(fluxtopclrsky-btcmin)
            taumintmp = -1. * log(max(min(transmax,(1.-qmin)),0.001))
            
            if (transmax .gt. 0.001 .and. transmax .le. (1.-qmin)) then
                    emcld = 1. - exp(-1. * tau(ibox) / 2.13 )
                  
                    fluxtop = fluxtop -   ((1.-emcld)*fluxtopclrsky)
                    fluxtop = max(qmin,fluxtop / emcld)

                    
            end if

            if (tau(ibox) .gt. (-1.*log((1.-qmin)))) then 
                
                !cloudy box
                tb(ibox)= 1307.27/ (log(1. + (1./fluxtop)))
                
                if (tau(ibox) .lt. taumintmp) then
                         tb(ibox) = attrop - 5.
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
!     (adjust_isccp_top_height = F)
!     or the radiatively determined cloud top pressure 
!     (adjust_isccp_top_height = T)
!



      !compute cloud top pressure
      do ibox=1,ncol
      
               !segregate according to optical thickness
               if (tau(ibox).le. (-1.*log((1.-qmin)))) then

                         ptop(ibox)=0.
                         levmatch(ibox)=0      

               else 

                     if (adjust_isccp_top_height) then  
                                               
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
                          do while(ptop(ibox) .eq. 0.        &
                                    .and. ilev .lt. nlev+1)
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
            fq_isccp(itau,ipres)=fq_isccp(itau,ipres)+(1./real(ncol))
            end if

      end if
                       
      end do
      
END SUBROUTINE ISCCP_CLOUDTYPES

!########################################################################
!########################################################################

SUBROUTINE TAU_REFF_DIAG(is,js,Time,coszen,nclds,ktop_phys,&
          kbot_phys,cldamt,TKel,phalf,tau,tau_ice,Reff_liq,Reff_ice)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine computes the typical cloud optical
!      depth, effective radius and temperature for each
!      column of the atmosphere, and for liquid and ice 
!      clouds separately.
!
!      The method is to select the compacted cloud with the
!      greatest horizontal area coverage and use the properties
!      of that cloud in determining 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	VARIABLES
!
!       ------
!	INPUT:
!       ------
!
!       is,js        (is/js) indices for model slab
!       coszen       cosine of zenith angle
!       nclds        number of clouds in column
!       ktop_phys    physical tops of compacted clouds
!       kbot_phys    physical bottoms of compacted clouds
!       cldamt       amount of clouds in matrix (fraction)
!       TKel            temperature (Kelvin)
!       phalf        pressure on model half levels (Pa)
!       tau          visible optical depth (dimensionless)
!       tau_ice      visible ice cloud optical depth (dimensionless)
!       Reff_liq     liquid cloud effective radius (microns)
!       Reff_ice     ice cloud effective radius (microns)
!
!       -------------------
!	INTERNAL VARIABLES:
!       -------------------
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

INTEGER,  INTENT (IN)                    :: is, js
type(time_type), intent(in)              :: Time
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: TKel,phalf
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: cldamt,tau,tau_ice
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: Reff_liq,Reff_ice
REAL,     INTENT (IN), DIMENSION(:,:)    :: coszen
INTEGER,  INTENT (IN), DIMENSION(:,:)    :: nclds
INTEGER,  INTENT (IN), DIMENSION(:,:,:)  :: ktop_phys,kbot_phys

!  Internal variables
!  ------------------

INTEGER                                  :: i,j,k,IDIM,JDIM,KDIM,ifld,maxclds
REAL, DIMENSION(SIZE(TKel,1),SIZE(TKel,2))   :: &
                                            tau_cld_liq, tau_cld_ice, &
                                            reff_cld_liq, reff_cld_ice, &
                                            T_cld_low,T_cld_lay
REAL, DIMENSION(SIZE(TKel,1),SIZE(TKel,2),SIZE(TKel,3))   :: tau_liq

REAL, DIMENSION(SIZE(TKel,1),SIZE(TKel,2))   :: cldamt_space_liq, &
                                                cldamt_space_ice, &
                                                cldamt_space_low, &
                                                tca_space,tca_space_last,&
                                                aliq_space,&
                                                aice_space,&
                                                alow_space,&
                                                Tlay,sum_logp
REAL :: tca_space_loc, tca_space_loc_high, tca_space_loc_last
INTEGER, DIMENSION(SIZE(TKel,1),SIZE(TKel,2))   :: ncld_liq,ncld_ice
LOGICAL :: used

!
! Code
! ----

    


    !1. initialize some variables
    IDIM=SIZE(TKel,1)
    JDIM=SIZE(TKel,2)
    KDIM=SIZE(TKel,3)
    tau_cld_liq(:,:)  = 0.
    tau_cld_ice(:,:)  = 0.
    reff_cld_liq(:,:) = 0.
    reff_cld_ice(:,:) = 0.
    T_cld_lay(:,:)  = 0.
    T_cld_low(:,:)  = 0.
    aliq_space(:,:) = 0.
    aice_space(:,:) = 0.
    alow_space(:,:) = 0.
    maxclds = MAXVAL(nclds(:,:))
    
    if (maxclds .gt. 0) then    
    
    
    !compute liquid cloud optical depth
    tau_liq = MAX(tau - tau_ice,0.)

    !compute properties of liquid and ice cloud reff
    if ( (id_aice > 0) .or. (id_aliq > 0)) then
  
       tca_space(:,:) = 0.
       do k = 1, maxclds

            !update cld amount from space calculation
            tca_space_last = tca_space
            tca_space(:,:) = 1.-((1.-tca_space(:,:))*(1.-cldamt(:,:,k)))
            WHERE ((tau_liq(:,:,k)>0.).and.(coszen(:,:)>1.E-06).and.&
                   (tau_ice(:,:,k)<=0.))
            cldamt_space_liq(:,:) = MIN(MAX(0.,tca_space-tca_space_last),1.)
            ELSEWHERE
            cldamt_space_liq(:,:) = 0.
            END WHERE
            WHERE ((tau_ice(:,:,k)>0.).and.(coszen(:,:)>1.E-06).and.&
                   (tau_liq(:,:,k)<=0.))
            cldamt_space_ice(:,:) = MIN(MAX(0.,tca_space-tca_space_last),1.)
            ELSEWHERE
            cldamt_space_ice(:,:) = 0.
            END WHERE
           
            !Sum up areas and particle sizes. Note that particle sizes
            !are multiplied by the area of the cloud seen from space so 
            !that division by the area will yield the local cloud properties
            aliq_space(:,:) = aliq_space(:,:) + cldamt_space_liq(:,:)
            reff_cld_liq(:,:) = reff_cld_liq(:,:) + &
                               cldamt_space_liq(:,:)*Reff_liq(:,:,k) 
            aice_space(:,:) = aice_space(:,:) + cldamt_space_ice(:,:)
            reff_cld_ice(:,:) = reff_cld_ice(:,:) + &
                               cldamt_space_ice(:,:)*Reff_ice(:,:,k) 
                                 
       end do  !loop over clouds

    end if   !for doing liquid and ice cloud properties


    if (id_alow > 0) then

       !compute low cloud temperature as the log pressure weighted
       !mean temperature of layers where the pressure is greater than
       !680 mb.  This choice conforms to Isccp definitions of low cloud
       sum_logp(:,:)=0.
       Tlay(:,:)=0.
       do k = 1, KDIM
          WHERE (((phalf(:,:,k+1)+phalf(:,:,k))/2.).ge. 68000.)
                  sum_logp(:,:)=sum_logp(:,:)+ &
                        (log(phalf(:,:,k+1))-log(phalf(:,:,k)))
                  Tlay(:,:)=Tlay(:,:)+ TKel(:,:,k)* &
                        (log(phalf(:,:,k+1))-log(phalf(:,:,k)))
          ENDWHERE
       enddo
       WHERE(sum_logp(:,:) .gt. 0.)
          Tlay(:,:)=Tlay(:,:)/sum_logp(:,:)
       ELSEWHERE
          Tlay(:,:)=0.
       ENDWHERE
    
       !loop over column becase algorithm cannot be made more
       !vectorized
        
       do i = 1,IDIM
       do j = 1,JDIM
             
            tca_space_loc = 0.
            tca_space_loc_high = 0.
            do k = 1,nclds(i,j)

                  !update cld amount from space calculation
                  tca_space_loc_last = tca_space_loc
                  tca_space_loc = 1.-((1.-tca_space_loc)*(1.-cldamt(i,j,k)))
                  
                  !should high cloud amount be set
                  if (((phalf(i,j,ktop_phys(i,j,k)+1)+ &
                        phalf(i,j,ktop_phys(i,j,k)))/2.) .lt. 68000.) then
                           tca_space_loc_high = tca_space_loc
                  end if
        
                  !is this a low cloud visible under sunlit conditions
                  if (((phalf(i,j,ktop_phys(i,j,k)+1)+ &
                        phalf(i,j,ktop_phys(i,j,k)))/2.) .ge. 68000. &
                        .and. (coszen(i,j)>1.E-06)) then

                       tau_cld_ice(i,j)=tau_cld_ice(i,j)+ &
                            (cldamt(i,j,k)*(1.-tca_space_loc_high)*&
                             tau_ice(i,j,k))
                       tau_cld_liq(i,j)=tau_cld_liq(i,j)+ &
                            (cldamt(i,j,k)*(1.-tca_space_loc_high)*&
                             tau_liq(i,j,k))
                       T_cld_low(i,j) = T_cld_low(i,j)+ &
                             TKel(i,j,ktop_phys(i,j,k))* &
                             (tca_space_loc-tca_space_loc_last)
                       
                  end if

            end do !loop over k

            !record low cloud area and temperature
            if (coszen(i,j)>1.E-06) then
            alow_space(i,j) = tca_space_loc - tca_space_loc_high
            T_cld_lay(i,j) = (tca_space_loc - tca_space_loc_high)* Tlay(i,j)
            end if

       enddo  !loop over j
       enddo  !loop over i

    end if   !for low cloud diagnostics

    end if   !for any clouds in the column
    
    
    !netcdf diagnostics....

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

END SUBROUTINE TAU_REFF_DIAG

!########################################################################
!########################################################################
 
FUNCTION ran0(idum)


!     $Id: cloud_rad.F90,v 1.6 2001/07/05 17:18:16 fms Exp $
!     Platform independent random number generator from
!     Numerical Recipies
!     Mark Webb July 1999
      
      REAL :: ran0
      INTEGER, INTENT (INOUT) :: idum
      
      INTEGER :: IA,IM,IQ,IR,k
      REAL :: AM

      IA=16807; IM=2147483647; AM=1.0/IM; IQ=127773; IR=2836
      
      if (idum.eq.0) then
        write(6,*) 'idum=',idum
        write(6,*) 'ZERO seed not allowed'
        stop
      endif

      k=idum/IQ
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      ran0=am*idum

END FUNCTION ran0
      

!########################################################################
!########################################################################
 

END MODULE CLOUD_RAD_MOD


      

