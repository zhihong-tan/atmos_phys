MODULE STRAT_CLOUD_MOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!       PROGNOSTIC CLOUD SCHEME
!
!
!       April 2000
!       Contact person: Steve Klein
!
!
!	The prognostic scheme returns the time tendencies of liquid,
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
!	This scheme is based on the experience the author had 
!       at the ECMWF in 1997. The saturated volume fraction formalism 
!       and type of solution follow directly from the scheme of Tiedtke
!       (1993): Monthly Weather Review, Volume 121, pages 3040-3061.
!       The form of most of the microphysics follows Rotstayn , 1997:
!       Quart. J. Roy. Met. Soc. vol 123, pages 1227-1282.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!       Outside modules used include:
!
!
!       (1) ESCOMP, DESCOMP: these return the saturation vapor presure
!                            and its temperature derivative given the
!                            temperature
!
!       (2) Utilities_Mod  : File handling utility which verifies the
!                            existence of a file, and returns a unit 
!                            number
!       (3) Constants_Mod
!
!       (4) cloud_rad_mod    need to initialize values of qmin, N_land
!                            and N_ocean

        USE  SAT_VAPOR_PRES_MOD, ONLY :  ESCOMP,DESCOMP
        USE  Utilities_Mod,      ONLY :  File_Exist, Open_File,  &
                                         error_mesg, FATAL, &
                                         get_my_pe, Close_File, &
                                         read_data, write_data
        USE  Constants_Mod,      ONLY :  RDgas,RVgas,HLv,HLf,HLs,Cp,&
                                         Grav,Tfreeze,Dens_h2o
        USE  CLOUD_RAD_MOD,      ONLY :  CLOUD_RAD_INIT

        IMPLICIT NONE
        PRIVATE

        PUBLIC  STRAT_INIT,          &
                STRAT_DRIV,          &
                ADD_STRAT_TEND,      &
                SUBTRACT_STRAT_TEND, &
                STRAT_END,           &
                STRAT_CLOUD_SUM,     &
                STRAT_CLOUD_AVG,     &
                DO_STRAT_CLOUD
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!              GLOBAL STORAGE VARIABLES
!
!	radturbten    The sum of radiation and turbulent tendencies
!                     for each grid box. (K/sec)
!

        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: radturbten


!
!     ------ data for cloud averaging code ------
!

      real,    allocatable, dimension (:,:,:) :: qlsum, qisum, cfsum
      integer, allocatable, dimension (:,:)   :: nsum

!
!     ------ constants used by the scheme -------
!
      real, parameter :: d622 = RDgas / RVgas
      real, parameter :: d378 = 1. - d622
                                      
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
!	Grav	       gravitational acceleration      m/(s*s)
!
!       HLv            latent heat of vaporization     J/kg condensate
!
!	HLf	       latent heat of fusion           J/kg condensate
!
!	HLs	       latent heat of sublimation      J/kg condensate
!
!	RDgas	       gas constant of dry air         J/kg air/K
!
!	RVgas	       gas constant of water vapor     J/kg air/K
!
!	Cp	       specific heat of air at         J/kg air/K
!                      constant pressure
!
!       d622           RDgas divided by RVgas          dimensionless
!
!       d378           One minus d622                  dimensionless
!
!       Tfreeze        Triple point of water           K
!
!       Dens_h2o       Density of pure liquid          kg/(m*m*m)
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
!       U00	       threshold relative humidity     fraction
!                      for cloud formation 
!
!       rthresh        liquid cloud drop radius        microns
!                      threshold for autoconversion
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
!       eros_scale     erosion rate constant for       1/sec
!                      cloud destruction
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
!                      {qmin, MAX(qa,ql,qi)} is
!                      resolved by the precision 
!                      of the numbers used.
!
!       Dmin           minimum permissible             dimensionless
!                      dissipation in analytic 
!                      integration of qa, ql, qi
!                      equations. This constant
!                      ONLY affects the method by
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

  
        REAL              :: U00            =  0.825
        REAL              :: rthresh        =  10.
        REAL              :: N_land         =  6.E+08
        REAL              :: N_ocean        =  1.E+08
        REAL              :: rho_ice        =  100.
        REAL              :: ELI            =  0.7
        REAL              :: U_evap         =  1.0
        REAL              :: eros_scale     =  1.E-06
        LOGICAL           :: tracer_advec   =  .false.
        REAL              :: qmin           =  1.E-10
        REAL              :: Dmin           =  1.E-08
        LOGICAL           :: do_average     =  .false.
        LOGICAL           :: strat_cloud_on =  .false.
        INTEGER,PARAMETER :: max_strat_pts  = 5
        INTEGER           :: num_strat_pts = 0
        INTEGER,DIMENSION(2,max_strat_pts) :: strat_pts = 0

!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!       CREATE NAMELIST
!

        NAMELIST /STRAT_CLOUD_NML/ U00,rthresh,N_land,N_ocean,&
                              rho_ice,ELI,U_evap,eros_scale,tracer_advec,&
                              qmin,Dmin,num_strat_pts,strat_pts          

!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!       DECLARE VERSION NUMBER OF SCHEME
!
        
        Character(len=128) :: Version = '$Id: strat_cloud.F90,v 1.2 2000/08/04 17:24:20 fms Exp $'
        Character(len=128) :: Tag = '$Name: bombay $'
        
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!       The module contains the following SUBROUTINES:
!
!
!       STRAT_INIT    read namelist file, open logfile, initialize
!                     constants and fields, read restart file
!
!       STRAT_DRIV    calculations of the cloud scheme are performed 
!                     here
!
!       ADD_STRAT_TEND  
!                     Adds a field to radturbten. This subroutine
!                     is needed because of the method to calculate
!                     the radiative and turbulent tendencies
!
!       SUBTRACT_STRAT_TEND 
!                     Subtracts a field from radturbten.
!
!       STRAT_END     writes out radturbten to restart file
!
!       STRAT_CLOUD_SUM
!                     sum cloud scheme variables
!
!       STRAT_CLOUD_AVG
!                     return average of summed cloud scheme variables
!
!       DO_STRAT_CLOUD
!                     logical flag, is the scheme on?
!


CONTAINS


!#######################################################################
!#######################################################################



SUBROUTINE STRAT_INIT(IDIM,JDIM,KDIM)
                               

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
!	VARIABLES
!
!
!       ------
!	INPUT:
!       ------
!	
!         variable              definition                  unit
!       ------------   -----------------------------   ---------------
!
!       IDIM,JDIM      number of points in first 
!                      and second dimensions
!
!       KDIM           number of points in vertical
!                      dimension
!
!
!       -------------------
!	INTERNAL VARIABLES:
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
!       User Interface variables
!       ------------------------
!

        INTEGER, INTENT (IN)                   :: IDIM,JDIM,KDIM
        
!
!       Internal variables
!       ------------------
!

        INTEGER                                :: unit,io


!-----------------------------------------------------------------------
!       
!       Code
!

!-----------------------------------------------------------------------
!
!       Namelist functions

        !read namelist if it exists
        If (File_Exist('input.nml')) Then
             unit = Open_File ('input.nml', action='read')
             io=1
             Do While (io .ne. 0)
                   Read  (unit, nml=STRAT_CLOUD_NML, iostat=io, End=10)
             EndDo
  10         Call Close_File (unit)
        EndIf

        !write namelist variables to logfile
        unit = Open_File ('logfile.out', action='append')
        if ( get_my_pe() == 0 ) then
           Write (unit,'(/,80("="),/(a))') trim(Version), trim(Tag)
           Write (unit,nml=STRAT_CLOUD_NML)
        endif
        Call Close_File (unit)

!-----------------------------------------------------------------------
!
!       INITIALIZE qmin, N_land, N_ocean and selected physical constants
!       in cloud_rad_mod

        CALL CLOUD_RAD_INIT(qmin,N_land,N_ocean)

!-----------------------------------------------------------------------
!
!       INITIALIZE strat_cloud_on to true

        strat_cloud_on = .TRUE.

!-----------------------------------------------------------------------
!
!       Read Restart file

        
       !set up stratiform cloud storage
       if (ALLOCATED(radturbten)) DEALLOCATE (radturbten)
           ALLOCATE(radturbten(IDIM,JDIM,KDIM))

           ALLOCATE( nsum(IDIM,JDIM),      &
                    qlsum(IDIM,JDIM,KDIM), &
                    qisum(IDIM,JDIM,KDIM), &
                    cfsum(IDIM,JDIM,KDIM)  )
       
       !see if restart file exists for radturbten
       If (File_Exist('INPUT/strat_cloud.res')) Then
                unit = Open_File (FILE='INPUT/strat_cloud.res', &
                                  FORM='native', ACTION='read')
                call read_data (unit, radturbten)
                call read_data (unit, nsum)
                call read_data (unit, qlsum)
                call read_data (unit, qisum)
                call read_data (unit, cfsum)
                Call Close_File (unit)
       ELSE
                radturbten(:,:,:)=0.
                qlsum=0.0; qisum=0.0; cfsum=0.0; nsum=0
       EndIf


!-----------------------------------------------------------------------
!
!
!       END OF SUBROUTINE
!
!


END SUBROUTINE STRAT_INIT



!#######################################################################
!#######################################################################



SUBROUTINE STRAT_DRIV(is,ie,js,je,dtcloud,pfull,phalf,T,qv,ql,qi,qa,&
                      omega,Mc,LAND,ST,SQ,SL,SI,SA,surfrain,surfsnow,&
                      MASK)
                               

!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This subroutine is the main driver of the cloud scheme.
!       The subroutine returns the time tendencies of temperature,
!       water vapor, cloud condensate, saturated volume
!       fraction, and surface values of rain and snow.
!
!       An analytic integration of the equations of cloud liquid and
!       ice and saturated volume fraction is performed following  
!       Tiedtke (1993): Monthly Weather Review, Volume 121, page 3046.
!       The difference of the result of the analytic integration with
!       the input values are the increments (changes) of the cloud 
!       prognostic variables, water vapor, and temperature which are
!       returned in the fields ST,SQ,SL,SI,SA.
!
!       The formulas for individual microphysics and large-scale 
!       processes are detailed below in comment statements.
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!	VARIABLES
!
!
!
!       ------
!	INPUT:
!       ------
!	
!         variable              definition                  unit
!       ------------   -----------------------------   ---------------
!
!       is,ie          starting and ending i indices 
!                      for data window
!
!       js,je          starting and ending j indices 
!                      for data window
!
!       dtcloud        time between this call and      s
!                      the next call to STRAT_DRIV
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
!       omega          vertical pressure velocity      Pa/s
!
!       Mc             Cumulus mass flux defined       kg/(m*m)/s
!                      on same levels as phalf
!
!       LAND           the fraction of the grid box    fraction
!                      covered by land
!                               
!
!
!       -------
!	OUTPUT:
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
!       surfrain       rain that falls through the     kg condensate/
!                      bottom of the column over       (m*m)
!                      the time dtcloud
!
!       surfsnow       snow that falls through the     kg condensate/
!                      bottom of the column over       (m*m)
!                      the time dtcloud
!
!
!       ---------------
!	OPTIONAL INPUT:
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
!
!       -------------------
!	INTERNAL VARIABLES:
!       -------------------
!	
!         variable              definition                  unit
!       ------------   -----------------------------   ---------------
!
!	KDIM	       number of vertical levels
!
!       j	       model vertical level being 
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
!       airdens        air density                     kg air/(m*m*m)
!
!       qs             saturation specific humidity    kg vapor/kg air
!
!       dqsdT          T derivative of qs              kg vapor/kg air/K
!
!       gamma          1. + (L/Cp)*dqsdT               dimensionless
!
!	rain_clr       grid mean flux of rain enter-   kg condensate/
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
!       deltp	       pressure thickness of grid box  Pa
!
!	U	       grid box relative humidity      fraction
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
!       lamda_f        slope factor in the size        1/m
!                      distribution of ice crystals
!
!       U_clr          relative humidity in the clear  fraction
!                      portion of the grid box.
!
!       tmp1,tmp2      temporary numbers used at 
!                      several points within the
!                      subroutine
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


!        
!       User Interface variables
!       ------------------------
!

        INTEGER, INTENT (IN)                   :: is,ie,js,je
        REAL, INTENT (IN)                      :: dtcloud
        REAL, INTENT (IN),    DIMENSION(:,:,:) :: pfull,phalf
        REAL, INTENT (IN),    DIMENSION(:,:,:) :: T,qv,ql,qi,qa,omega
        REAL, INTENT (IN),    DIMENSION(:,:,:) :: Mc
        REAL, INTENT (IN),    DIMENSION(:,:)   :: LAND
        REAL, INTENT (OUT),   DIMENSION(:,:,:) :: ST,SQ,SL,SI,SA
        REAL, INTENT (OUT),   DIMENSION(:,:)   :: surfrain,surfsnow
        REAL, INTENT (IN), OPTIONAL, DIMENSION(:,:,:) :: MASK
        
!
!       Internal variables
!       ------------------
!

        INTEGER                                :: KDIM,j,ipt,jpt
        INTEGER                                :: i,unit,nn
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: airdens
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: qs,dqsdT
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: gamma
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: A_plus_B
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: rain_clr,rain_cld
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: a_rain_clr,a_rain_cld
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: snow_clr,snow_cld
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: a_snow_clr,a_snow_cld
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: deltp,U
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: dqs_ls,da_ls
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: C_dt, D_dt
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: D1_dt,D2_dt,D_eros
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: ql_upd,qi_upd,qa_upd
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: ql_mean,qi_mean
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: qa_mean
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: dcond_ls,dcond_ls_ice        
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: da_cld2clr,da_clr2cld
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: dprec_clr2cld
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: dprec_cld2clr
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: N,rad_liq
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: Vfall,lamda_f
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: U_clr
        REAL, DIMENSION(SIZE(T,1),SIZE(T,2))   :: tmp1,tmp2
        

!-----------------------------------------------------------------------
!       
!       Code
!

!-----------------------------------------------------------------------
!
!       Determine number of vertical levels

        KDIM = SIZE(T,3)


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
!   (3) gamma = 1.    +   (L/Cp) * dqsdT
!       
!       where d622 = RDgas/RVgas; esat = saturation vapor pressure;
!       and desat/dT is the temperature derivative of esat.
!       Note that in the calculation of gamma, 
!
!            {             HLv          for T > Tfreeze             }
!       L =  { 0.05*(T-Tfreeze+20.)*HLv + 0.05*(Tfreeze-T)*HLs      }
!            {                          for Tfreeze-20.< T < Tfreeze}
!            {             HLs          for T < Tfreeze-20.         }
!
!       This linear form is chosen because at Tfreeze-20. es = esi, and
!       at Tfreeze, es = esl, with linear interpolation in between.
!
!       The conductivity/diffusivity factor, A_plus_B is given by:
!
!   (4) A_plus_B =   { (HLv/Ka/T)*((HLv/RVgas/T)-1.) } + 
!
!                    { (RVgas*T/chi*esat) }
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
        CALL ESCOMP(T(:,:,:),gamma(:,:,:))
        
        !compute A_plus_B
        A_plus_B(:,:,:) = ( (HLv/0.024/T(:,:,:)) *     &
                           ((HLv/RVgas/T(:,:,:)) -1.) ) &
                   +  (RVgas*T(:,:,:)*pfull(:,:,:)/2.21/gamma(:,:,:))  
    
        !calculate denominator in qsat formula
        qs(:,:,:) = pfull(:,:,:)-(d378*gamma(:,:,:))
     
        !limit denominator to esat, and thus qs to d622
        !this is done to avoid blow up in the upper stratosphere
        !where pfull ~ esat
        qs(:,:,:) = MAX(qs(:,:,:),gamma(:,:,:)) 
        
        !compute temperature derivative of esat 
        CALL DESCOMP(T(:,:,:),dqsdT(:,:,:))
        
        !calculate dqsdT
        dqsdT(:,:,:)=d622*pfull(:,:,:)*dqsdT(:,:,:)/ &
                       qs(:,:,:)/qs(:,:,:)

        !calculate qs
        qs(:,:,:)=d622*gamma(:,:,:)/qs(:,:,:)
         
        !calculate gamma
        gamma(:,:,:) =  1.  + (HLv/Cp)*dqsdT(:,:,:)
        
        WHERE (T(:,:,:) .ge. Tfreeze-20. .and. &
               T(:,:,:) .lt. Tfreeze)
             gamma(:,:,:) =  1.  +   (    &
                       (0.05*(T(:,:,:)-Tfreeze+20.)*HLv + &
                        0.05*(Tfreeze-T(:,:,:))*HLs)*dqsdT(:,:,:)/Cp)
        END WHERE

        WHERE (T(:,:,:) .lt. Tfreeze-20.)
             gamma(:,:,:) =  1.  + (HLs/Cp)*dqsdT(:,:,:)
        END WHERE
        
        
!-----------------------------------------------------------------------
!
!       Calculate air density ignoring contributions from vapor and
!       and condensed phase.  This can be done because this is only
!       used a scaling term.
!
!   (6) airdens  =   pfull / RDgas / T

        airdens(:,:,:) = pfull(:,:,:)/RDgas/T(:,:,:)

!-----------------------------------------------------------------------
!
!       Assign cloud droplet number based on land or ocean point.

        
        N(:,:)=N_land*LAND(:,:) + N_ocean*(1.-LAND(:,:))

!-----------------------------------------------------------------------
!
!       Initialize select variables to zero. The variables reset
!       are:
!
!       (1) changes of prognostic variables
!
!       (2) variables dealing with the rain/snow fluxes. 
!
!       (3) qa_mean of the level above the top level.

        ST(:,:,:) = 0.
        SQ(:,:,:) = 0.
        SL(:,:,:) = 0.
        SI(:,:,:) = 0.
        SA(:,:,:) = 0.
   
        rain_cld(:,:)   = 0.
        rain_clr(:,:)   = 0.
        a_rain_cld(:,:) = 0.
        a_rain_clr(:,:) = 0.
        snow_cld(:,:)   = 0.
        snow_clr(:,:)   = 0.
        a_snow_clr(:,:) = 0.
        a_snow_cld(:,:) = 0.
        
        
!-----------------------------------------------------------------------
!
!       Enter the large loop over vertical levels.  Level 1 is the top
!       level of the column and level KDIM is the bottom of the model.
!       If MASK is present, each column may not have KDIM valid levels.

        DO j = 1, KDIM

       
!-----------------------------------------------------------------------
!
!       Calculate pressure thickness of level and relative humidity 

        !calculate difference in pressure across the grid box
        deltp(:,:) = phalf(:,:,j+1)-phalf(:,:,j)
         
        !calculate GRID box mean relative humidity 
        U(:,:) = MIN(MAX(0.,qv(:,:,j)/qs(:,:,j)),1.)
        
!-----------------------------------------------------------------------
!
!       Account for the fact that other processes may have created
!       negative tracer or extremely small values of tracer fields.
!       The general reason for the extremely small values of the 
!       tracer fields is due to advection of condensate or cumulus
!       induced subsidence (also a form of advection) of condensate.  
!       In this step any values of prognostic variables which are 
!       less than qmin are reset to zero, while conserving total 
!       moisture.
!
!       Also correct for qa > RH, which is not permitted under the
!       assumption that the cloudy air is saturated and the temperature
!       inside and outside of the cloud are about the same.

        WHERE (qa(:,:,j) .le. qmin)
               SA(:,:,j)   = SA(:,:,j) - qa(:,:,j)
               qa_upd(:,:) = 0.
        ELSEWHERE
               qa_upd(:,:) = qa(:,:,j)      
        END WHERE
        
        WHERE (qa_upd(:,:) .gt. U(:,:))
               SA(:,:,j)   = SA(:,:,j) + U(:,:) - qa(:,:,j)
               qa_upd(:,:) = U(:,:)      
        END WHERE
                
        WHERE (ql(:,:,j) .le. qmin .or. qa(:,:,j) .le. qmin)
               SL(:,:,j)   = SL(:,:,j) - ql(:,:,j)
               SQ(:,:,j)   = SQ(:,:,j) + ql(:,:,j)
               ST(:,:,j)   = ST(:,:,j) - HLv*ql(:,:,j)/Cp
               ql_upd(:,:) = 0.
        ELSEWHERE
               ql_upd(:,:) = ql(:,:,j)
        END WHERE
        
        WHERE (qi(:,:,j) .le. qmin .or. qa(:,:,j) .le. qmin)
               SI(:,:,j)   = SI(:,:,j) - qi(:,:,j)
               SQ(:,:,j)   = SQ(:,:,j) + qi(:,:,j)
               ST(:,:,j)   = ST(:,:,j) - HLs*qi(:,:,j)/Cp
               qi_upd(:,:) = 0.
        ELSEWHERE
               qi_upd(:,:) = qi(:,:,j)
        END WHERE
        
        
!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!                      NON-CONVECTIVE CONDENSATION                     !
!                                                                      !
!                                   AND                                !
!                                                                      !
!          ANALYTIC INTEGRATION OF SATURATED VOLUME FRACTION EQUATION  !
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
!       and radiative cooling.  For the moment the radiative cooling
!       term is not included.  dqs_ls has the form:
!
!   (7) dqs_ls   =   (omega+ Grav*Mc) * dtcloud * (dqs/dp)_ma
!
!
!       where (dqs/dp)_ma is the change of qs following a moist adiabat
!       and is given by the form:
!
!   (8) (dqs/dp)_ma = dqsdT/Cp/airdens/gamma  
!

        
        dqs_ls(:,:) = (omega(:,:,j)+0.5*Grav*(Mc(:,:,j)+Mc(:,:,j+1))) &
               *dtcloud*dqsdT(:,:,j)/airdens(:,:,j)/gamma(:,:,j)/Cp

!       dqs_ls is augmented by radiation and turbulent tendencies
!
!   (9) dqs_ls = dqs_ls + radturbten * dtcloud * dqsdT
!

        dqs_ls(:,:) = dqs_ls(:,:) + &
                      radturbten(is:ie,js:je,j)*dtcloud*dqsdT(:,:,j)
       
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
!   (9)      da_ls = - (1.-qa)*(1.-qa)*dqs_ls/2./qs/(1.-U)
!
!            This formula arises from the assumption that vapor is uni-
!            formly distributed in the range [qv_clr - (qs - qv_clr),qs]
!            where qv_clr is the amount of vapor in the unsaturated 
!            volume and is given from the following equation:
!
!  (10)      qv  =   qa * qs      +   (1.-qa) * qv_clr
!          
!            Implicit in equation (9) is the following assumption:
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
        WHERE (dqs_ls(:,:) .le. 0. .and. U(:,:) .ge. U00)
             da_ls(:,:) = -0.5 * (1.-qa_upd(:,:))* (1.-qa_upd(:,:))* &
                    dqs_ls(:,:) / qs(:,:,j) / MAX(1.-U(:,:),qmin)
        ELSEWHERE
             da_ls(:,:) = 0.
        END WHERE 

!       Turbulent erosion of clouds
!
!       As in Tiedtke (19883) this is calculated using the eros_scale
!       parameter as:
!
!  (12) dql/dt    =  - qa * eros_scale * (qs - qv) * (ql/ ql+qi )
!
!  (13) dqi/dt    =  - qa * eros_scale * (qs - qv) * (qi/ ql+qi )
!
!  (14) dqa/dt    =  - qa * eros_scale * (qs - qv) * (qa/ ql+qi )
!
!       for which the erosion sink term (B in equation 16) is
!
!  (15) B = qa * eros_scale * (qs - qv) / (ql + qi)  
!

        WHERE (ql_upd(:,:) .gt. qmin .or. qi_upd(:,:) .gt. qmin)
             D_eros(:,:) = qa_upd(:,:) * eros_scale * dtcloud *&
                           qs(:,:,j) * (1.-U(:,:)) / &
                           (qi_upd(:,:) + ql_upd(:,:))
        ELSEWHERE
             D_eros(:,:) = 0.
        END WHERE
      
!       The next step is to analytically integrate the saturated volume
!       fraction equation. Following Tiedtke, the qa equation is written 
!       in the form:
!
!  (12) dqa/dt    =   (1.-qa) * A    -   qa * B
! 
!       For which the analtyic solution is given by:
!
!  (13) qa(t+dtcloud) =   qa(t)  *       exp (-(A+B)*dtcloud)    +
!                       (A/(A+B))* (1. - exp (-(A+B)*dtcloud) )
!
!       Or that 
!
!  (14) SA = qa(t+dtcloud)-qa(t) =  
!                           
!       SA = (  (A/(A+B)) - qa(t) ) * (1. - exp (-(A+B)*dtcloud) )
!                        
!
!       Note that convective detrainment contributes to the source
!       of saturated volume.
!
!       Note that in the below, A and B are computed as A*dtcloud
!       and B*dtcloud for ease of solution.  Following equations
!       arise from:
!
!  (15) A * (1. - qa) = da_ls / dtcloud      if da_ls >= 0.
!     
!       and 
!
!  (16) B * qa        = da_ls / dtcloud      if da_ls < 0.
!
!       Note that to reduce the number of variables, A is called C here,
!       and B is called D here. At present, B = 0 because no large-scale
!       or turbulent cloud evaporation is included.  qa goes to zero only
!       in the case of ql and qi less than or equal to qmin; see code 
!       near the end of this loop over levels.
        
        !compute C_dt; This is assigned to the large-scale source term
        !following (15). Reset D_dt.
        C_dt(:,:) = da_ls(:,:)/MAX((1.-qa_upd(:,:)),qmin)
        D_dt(:,:) = D_eros(:,:)

        !compute formula for total change in da
        WHERE (D_dt(:,:) .gt. Dmin)
             tmp1(:,:) =  &
               ( (C_dt(:,:)/(C_dt(:,:)+D_dt(:,:))) - qa_upd(:,:))  *   &
               (  1.   -      exp( -1.* (C_dt(:,:)+D_dt(:,:))))
        ELSEWHERE
             tmp1(:,:) = (1.-qa_upd(:,:))*(1.-exp(-1.*C_dt(:,:)))
        END WHERE
          
!       The next step is to calculate the change in condensate
!       due to non-convective condensation, dcond_ls. Note that this is
!       not the final change but is used only to apportion condensate
!       change between phases. According to Tiedtke 1993 this takes the
!       form:
!
!  (17) dcond_ls =  -1. (qa +  0.5*da_ls) dqs_ls
!
!       For ease of calculation, da_ls as in (9) is used instead
!       of an exact calculation using information of SA. That is
!       to exactly calculate da_ls at this point would require
!       subtracting the portions of the SA which are not due
!       large-scale condensation.  For example if turbulent 
!       cloud evaporation were present one would use:
!
!       tmp1  + qa*(1.-exp(-B_turb*dt))
! 
!       as the approximation to da_ls in (17).  However to reduce
!       the number of exponentials calculated in the program,
!       this is not done.
        
        dcond_ls(:,:) = -1. * dqs_ls(:,:) * (qa_upd(:,:) + &
                              0.5*MIN(da_ls(:,:),1.-qa_upd(:,:)))
             
        !note that non-convective condensation/deposition can only 
        !occur where U > U00
        WHERE (dqs_ls(:,:) .lt. 0. .and. U(:,:) .lt. U00)
             dcond_ls(:,:) = 0.
        END WHERE
        
        !also set non-convective evaporation/sublimation to zero
        !where no condensate is present
        WHERE (dcond_ls(:,:) .lt. 0. .and. &
               (ql_upd(:,:)+qi_upd(:,:)) .le. qmin)
             dcond_ls(:,:) = 0.
        END WHERE
        
!       Now that qa_upd has been used to calculate (17), we can now
!       update it and SA.
        SA(:,:,j) = SA(:,:,j) + tmp1(:,:)
        qa_upd(:,:) = qa_upd(:,:) + tmp1(:,:) 

!       The next step is the apportionment on the non-convective 
!       condensation between liquid and ice phases. Following the
!       suggestion of Rostayn (2000), all condensation is in liquid
!       form as ice nuclei are generally limited in the atmosphere.
!       The droplets may subsequently be nucleated by homogeneous
!       freezing or be converted to ice by the Bergeron-Findeisan
!       mechanism.  
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
!       ice because if both are present in the same volume the saturation
!       vapor pressure over the droplet is higher than that over the
!       ice crystal
!
!       The fraction of large-scale condensation that is liquid
!       is stored in the temporary variable tmp1.   

        !assume liquid fractionation
        tmp1(:,:)=1.
        tmp2(:,:)=0.

        !For cases of cloud evaporation, set liquid evaporation
        !to preferentially occur first
        WHERE (dcond_ls(:,:) .lt. 0. .and. ql_upd(:,:) .gt. qmin)
             tmp1(:,:) = MIN(-1.*dcond_ls(:,:),ql_upd(:,:))/ &
                         MAX(-1.*dcond_ls(:,:),qmin)
             tmp2(:,:) = 1.-tmp1(:,:)
        END WHERE

        !prevent evaporation of zero liquid /ice clouds
        WHERE (dcond_ls(:,:) .lt. 0. .and. ql_upd(:,:) .le. qmin)
             tmp1(:,:) = 0.
        END WHERE
        WHERE (dcond_ls(:,:) .lt. 0. .and. qi_upd(:,:) .le. qmin)
             tmp2(:,:) = 0.
        END WHERE

        !calculate partitioning among liquid and ice to dcond_ls
        dcond_ls_ice(:,:)= tmp2(:,:)*dcond_ls(:,:)
        dcond_ls(:,:)    = tmp1(:,:)*dcond_ls(:,:)
        
!       The next step is to compute semi-implicit qa,ql,qi which are 
!       used in many of the formulas below.  This gives a somewhat 
!       implicitness to the scheme. Note that ql and qi are 
!       incremented only when dcond_ls > 0.
  
        qa_mean(:,:) = qa_upd(:,:)
        ql_mean(:,:) = ql_upd(:,:) + MAX(dcond_ls(:,:)    ,0.)        
        qi_mean(:,:) = qi_upd(:,:) + MAX(dcond_ls_ice(:,:),0.)
        
!-----                                                            -----! 
!                                                                      !
!                                  END OF                              !
!                                                                      !
!                        NON-CONVECTIVE CONDENSATION                   !
!                                                                      !
!                                   AND                                !
!                                                                      !
!           ANALYTIC INTEGRATION OF SATURATED VOLUME FRACTION EQUATION !
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
!       parameterization developed by Jakob and Klein in 1997.  It 
!       treats the grid mean flux of rain and snow separately according 
!       to the portion that exists in unsaturated air and the portion 
!       that exists in saturated air.  At the top of each level, some 
!       precipitation that was in unsaturated air in the levels above 
!       may enter saturated air and vice versa. These transfers are 
!       calculated here under the assumption of maximum overlap between 
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
!       NOTE: Random overlap of cloud and precipitation is assumed in
!             the equations below.
!
!       For snow fluxes the transfers are computed in exactly the same
!       manner.
!
         
       
!
!
!       Rain transfers are done first
!
!

        !-------------------------------
        !compute cloud to clear transfer
        da_cld2clr(:,:) = MIN(   a_rain_cld(:,:), &
                          MAX(0.,a_rain_cld(:,:)*(1. - qa_mean(:,:))))
        
        !-------------------------------
        !compute clear to cloud transfer
        da_clr2cld(:,:) = MIN(MAX(a_rain_clr(:,:)*qa_mean(:,:),0.),&
                                  a_rain_clr(:,:))
        
        !---------------------------------
        !calculate precipitation transfers
        dprec_cld2clr(:,:) = rain_cld(:,:)*(da_cld2clr(:,:)/&
                             MAX(a_rain_cld(:,:),qmin))

        dprec_clr2cld(:,:) = rain_clr(:,:)*(da_clr2cld(:,:)/ &
                             MAX(a_rain_clr(:,:),qmin))
        
        !----------------
        !add in transfers
        a_rain_clr(:,:)=a_rain_clr(:,:)+da_cld2clr(:,:)-da_clr2cld(:,:)
        a_rain_cld(:,:)=a_rain_cld(:,:)-da_cld2clr(:,:)+da_clr2cld(:,:)
        rain_clr(:,:)=rain_clr(:,:)+dprec_cld2clr(:,:)&
                                   -dprec_clr2cld(:,:)        
        rain_cld(:,:)=rain_cld(:,:)-dprec_cld2clr(:,:)&
                                   +dprec_clr2cld(:,:)

        

!
!
!       Snow transfers are done second
!
!
        

        !-------------------------------
        !compute cloud to clear transfer
        da_cld2clr(:,:) = MIN(   a_snow_cld(:,:), &
                          MAX(0.,a_snow_cld(:,:)*(1. - qa_mean(:,:))))
        
        !-------------------------------
        !compute clear to cloud transfer
        da_clr2cld(:,:) = MIN(MAX(a_snow_clr(:,:)*qa_mean(:,:),0.), &
                                  a_snow_clr(:,:))

        !---------------------------------
        !calculate precipitation transfers
        dprec_cld2clr(:,:) = snow_cld(:,:)*(da_cld2clr(:,:)/ &
                             MAX(a_snow_cld(:,:),qmin))

        dprec_clr2cld(:,:) = snow_clr(:,:)*(da_clr2cld(:,:)/ &
                             MAX(a_snow_clr(:,:),qmin))
        
        !----------------
        !add in transfers
        a_snow_clr(:,:)=a_snow_clr(:,:)+da_cld2clr(:,:)-da_clr2cld(:,:)
        a_snow_cld(:,:)=a_snow_cld(:,:)-da_cld2clr(:,:)+da_clr2cld(:,:)
        snow_clr(:,:)=snow_clr(:,:)+dprec_cld2clr(:,:)&
                                   -dprec_clr2cld(:,:)
        snow_cld(:,:)=snow_cld(:,:)-dprec_cld2clr(:,:)&
                                   +dprec_clr2cld(:,:)
        
!-----------------------------------------------------------------------
!
!
!                           MELTING OF CLEAR SKY SNOW FLUX
!
!
!       Melting of falling ice to rain occurs when T 0C. The amount 
!       of melting is limited to the melted amount that would cool the 
!       temperature to Tfreeze.   In cloud melting of ice occurs
!       in the ice microphysics section.
!

        !compute grid mean change in snow flux to cool portion
        !of grid box with snow_clr flux to 2C and store in
        !temporary variable tmp1
        tmp1(:,:) = MAX((Cp*(T(:,:,j)-Tfreeze)*a_snow_clr(:,:)* & 
                         deltp(:,:) / Grav / dtcloud / HLf) ,0.)
        
        !change a_rain_clr
        WHERE (T(:,:,j) .gt. Tfreeze .and. &
               a_snow_clr(:,:) .gt. qmin .and. &
               a_rain_clr(:,:) .lt. a_snow_clr(:,:))
             a_rain_clr(:,:) = a_snow_clr(:,:)
        END WHERE

        !do melting in non-limited case first
        WHERE (snow_clr(:,:) .lt. tmp1(:,:))
             ST(:,:,j) = ST(:,:,j) - HLf*snow_clr(:,:)* &
                                Grav*dtcloud/deltp(:,:)/Cp                
             rain_clr(:,:) = rain_clr(:,:) + snow_clr(:,:)
             snow_clr(:,:) = 0.
             a_snow_clr(:,:) = 0.
        ELSEWHERE  !do flux limited case 
             ST(:,:,j) = ST(:,:,j) - HLf*tmp1(:,:)* &
                                Grav*dtcloud/deltp(:,:)/Cp                
             rain_clr(:,:) = rain_clr(:,:) + tmp1(:,:)
             snow_clr(:,:) = snow_clr(:,:) - tmp1(:,:)
        END WHERE
        
!----------------------------------------------------------------------!
!
!
!              COMPUTE SLOPE FACTOR FOR ICE MICROPHYSICS                  
!
!       [The following microphysics follows that of Rotstayn (1997)]
!       The number concentration of ice crystals of diameter D in the
!       size interval D to D+dD is assumed to be 
!       distributed as in a Marshall Palmer distribution :
!
!  (18) N(D)dD = Nof * Exp( - lamda_f * D)
!
!       The slope factor and intercept are not assumed to be constant,
!       but the slope factor is
!
!  (19) lamda_f = 1.6X10^(3.+0.023(Tfreeze-T))
!
!       Integration of (18) over all particle sizes with a constant
!       density of ice crystals , rho_ice, and assumed spherical shape
!       yields a relationship between the intercept parameter and qi
!
!  (20) Nof = airdens*qi_local*(lamda_f^4)/pi*rho_ice
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
        qi_mean(:,:) = qi_mean(:,:) + &
                       snow_cld(:,:)*Grav*dtcloud/deltp(:,:)        

        !compute lamda_f
        lamda_f(:,:) = 1.6 * 10**(3.+0.023*(Tfreeze-T(:,:,j)))
        

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
!  (21) V(D) = 141.4 (m**0.5,s-1) * sqrt(D) * (rho_ref/airdens)**0.5
!
!       where D is the radius of the rain drops in meters. Here
!       rho_ref is a reference air density.  This formula is generally
!       good for 1 mm < D < 4 mm.
!
!       The distribution of rain drops by size follows a Marshall-
!       Palmer distribution:
!
!  (22) N(D) dD = Nor * Exp (- lamda *D)
!
!       where N(D)dD is the number of rain drops with size D in the
!       interval D to D+dD, Nor is the intercept (assumed fixed at
!       8E+04 (1/m*m*m*m).  lamda is the slope intercept parameter
!       and with (21) it can be shown that lamda is a function of
!       the rain rate.
!
!       With these assumptions the local rate of accretion of cloud
!       liquid reduces to:
!
!  (23) dl/dt_local = - CB*Eco*((rain_rate_local/Dens_h2o)**(7/9))*
!                        
!                       (rho_ref/airdens)**(1/9)   * ql_local
!
!       where CB is the accretion constant:
! 
!       CB = 65.772565 [(m)**(-7/9)] * [(s)**(-2/9)]
!
!       AND Eco is the collection efficiency of a cloud droplet by a 
!       rain droplet.   A good fit to the Table 8.2 of Rogers and Yau 
!       (1988) for rain drops between size 1mm and 3mm is:
!
!  (24) Eco = rad_liq**2 / (rad_liq**2 + 20.5 microns**2) .
!
!       In generalizing to grid mean conditions (23) becomes:
!
!  (25) dl/dt = - (arain_cld/qa_mean) * CB * Eco * 
!
!                 [(rain_cld/a_rain_cld/Dens_h2o)**(7/9)] * ql
!        
!       Note that the very weak dependence on air density is
!       neglected at this point.
!
!       The particle sizes are computed from the following equation
!
!  (26) rad_liq = (3*airdens*ql/(4*pi*liq_dens*qa*N)^(1/3)
!
!       
!       For numerical treatment we write (25) as:
!
!       dl/dt = - D_acc * l 
!
!       and if we do so:
!
!  (27) D_acc   =  (arain_cld/qa_mean) * CB * Eco * 
!
!                 [(rain_cld/a_rain_cld/Dens_h2o)**(7/9)] 
!
!       In the work below, D_acc is added to D1_dt, the first processes
!       contributing to the depletion of cloud liquid in the analytic
!       integration.  D1_dt represents the conversion of liquid to rain.

        !compute rad_liq.  The constant below is equal to  
        !1.E+06 * (3/4*pi)^(1/3), where the 1E+06 is
        !the factor to convert meters to microns.
        
        rad_liq(:,:)= 620350.49 *( (airdens(:,:,j)*ql_mean(:,:)/ &
                      MAX(qa_mean(:,:),qmin)/N(:,:)/Dens_h2o)**(1./3.))
        
        !do not let very small cloud fractions contribution to
        !autoconversion or accretion
        WHERE (qa_mean(:,:) .le. qmin)
                      rad_liq(:,:) = 0.
        ENDWHERE

        !compute accretion D term
        D1_dt(:,:) =  dtcloud * 65.772565 * &
                     (a_rain_cld(:,:)/MAX(qa_mean(:,:),qmin)) * &
                     ( rad_liq(:,:)*rad_liq(:,:) / &
                      (rad_liq(:,:)*rad_liq+20.5) ) * &
                     ((rain_cld(:,:)/MAX(a_rain_cld(:,:),qmin)&
                      /Dens_h2o)**(7./9.))
                
!       Autoconversion
!
!       The autoconversion parameterization follow that of Manton
!       and Cotton (1977).  This formula has been used in Chen and
!       Cotton (1987) and is used in the CSIRO GCM (Rotstayn 1997)
!       and the LMD GCM (Boucher, Le Treut and Baker 1995).  In this
!       formulation the time rate of change of grid mean liquid is
!
!  (28) dl/dt= -CA * qa * [(ql/qa)^(7/3)] * [(N*Dens_h2o)^(-1/3)] 
!
!               * H(rad_liq - rthresh)
!
!       Where N is the number of cloud droplets per cubic metre,
!       rthresh is a particle radius threshold needed to for autoconv-
!       ersion to occur, H is the Heaviside function, and CA is
!       a constant which is:
!
!  (29) CA =  0.104 * Grav * Ec * (airdens)^(4/3) / mu 
!        
!       where Grav is gravitational acceleration, Ec is the collection
!       efficiency, airdens is the density of air, and mu is the 
!       dynamic viscosity of air.   This constant is evaluated 
!       ignoring the temperature dependence of mu and with a fixed 
!       airdens of 1 kg/m3.
!
!       With   Ec = 0.55        (standard choice - see references)
!            Grav = 9.81        m/(s*s)
!         airdens = 1.00        kg air/(m*m*m)
!              mu = 1.717  E-05 kg condensate/(m*s) 
!
!              CA = 32681. [(kg air)^(4/3)]/kg liq/m/s
!
!
!       For numerical treatment we write (28) as:
!
!       dl/dt = - D_aut * l 
!
!       and if we do so:
!
!  (30) D_aut   =   CA * [(N*Dens_h2o)^(-1/3)] * [(ql/qa)^(4/3)] * &
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
!  (31) D_dt = log( (rad_liq/rthresh)**3. )
!
!       This limiter is stored in tmp2.

        !compute autoconversion sink as in (30)
        tmp1(:,:)=    32681. * dtcloud * &
                      ((N(:,:)*Dens_h2o)**(-1./3.))* &
                      ((ql_mean(:,:)/MAX(qa_mean(:,:),qmin))**(4./3.))

        !compute limiter as in (31)
        tmp2(:,:) =MAX(3*log(MAX(rad_liq(:,:),qmin)/rthresh),0.)

        !limit autoconversion to the limiter
        tmp1(:,:) = MIN(tmp1(:,:),tmp2(:,:))
        
        !add autoconversion to D1_dt
        D1_dt(:,:) = D1_dt(:,:) + tmp1(:,:)
        
        
!       Bergeron-Findeisan Process 
!
!       Where ice and liquid coexist, the differential saturation
!       vapor pressure between liquid and ice phases encourages
!       the growth of ice crystals at the expense of liquid droplets.
!
!       Rotstayn (2000) show derive an equation for the growth of
!       ice by starting with the vapor deposition equation for an
!       ice crystal and write it in terms of ice specific humidity
!       as:
!
!                 {(Ni/airdens)**2/3}*7.8*((Max(qi,Mio*Ni/airdens))**(1/3))
!  (32) dqi/dt =  -------------------------------------------------------- X
!                          [rhoice**1/3]* A_plus_B
!
!                          [(esl-esi)/esi]
!
!       Here Ni is the ice crystal number which is taken from the 
!       parameterization of Meyers et al. :
!
!  (33) Ni = 1000 * exp( (12.96* [(esl-esi)/esi]) - 0.639 )
!
!       The use of the maximum operator assumed that there is a background
!       ice crystal always present on which deposition can occur.  Mio
!       is an initial ice crystal mass taken to be 10-12kg.
!
!       Figure 9.3 of Rogers and Yau (1998) shows the nearly linear
!       variation of [(esl-esi)/esi] from 0. at 273.16K to 0.5 at 
!       233.16K.  Analytically this is parameterized as (Tfreeze-T)/80.
!
!
!       Generalizing (32) to grid mean conditions and writing it in terms
!       of a loss of cloud liquid yields:
!
!                   qa {(1000 * exp((12.96* (Tfreeze-T)/80)-0.639)/airdens)**2/3}
!  (34) dql/dt = -  -------------------------------------------------------- 
!                            [rhoice**1/3]* A_plus_B * ql
!
!                  *7.8*((Max(qi/qa,Mio*Ni/airdens))**(1/3))*[(Tfreeze-T)/80]*ql
!
!       Note that the density of ice is set to 700 kg m3 the value appropriate
!       for pristine ice crystals.  This value is necessarily different than
!       the value of ice used in the riming and sublimation part of the code
!       which is for larger particles which have much lower densities.
        
        !reset D2_dt
        D2_dt(:,:)=0.
        
        !do Bergeron process
        WHERE (T(:,:,j) .lt. Tfreeze .and. ql_mean(:,:) .gt. qmin &
                    .and. qa_mean(:,:) .gt. qmin)
           
             D2_dt(:,:) =  dtcloud * qa_mean(:,:) * &
                         ((1000.*exp((12.96*0.0125*(Tfreeze-T(:,:,j)))-0.639)&
                           / airdens(:,:,j))**(2./3.))* &
                           7.8* ((MAX(qi_mean(:,:)/qa_mean(:,:), &
                           1.E-12*1000.* &
                           exp((12.96*0.0125*(Tfreeze-T(:,:,j)))-0.639)/&
                           airdens(:,:,j)))**(1./3.))*&
                           0.0125*(Tfreeze-T(:,:,j))/ &
                           ((700.**(1./3.))* A_plus_B(:,:,j)* ql_mean(:,:))
        END WHERE
       
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
!  (35) dql/dt = -  ( a_snow_cld / qa/a_snow_cld ) *
!                   ( ELI * lamda_f * snow_cld /
!                     2/ rho_ice ) * ql. 
!
!
!       Note that the temperature range of riming is limited to 
!       temperatures less than Tfreeze.  Freezing cannot take place
!       of rimed condensate for temperatures greater than this.       

        !add in accretion of cloud liquid by ice
        WHERE (a_snow_cld(:,:) .gt. qmin .and. &
                  ql_mean(:,:) .gt. qmin .and. &
                  qa_mean(:,:) .gt. qmin .and. &
                      T(:,:,j) .lt. Tfreeze)
             D2_dt(:,:) = D2_dt(:,:) + (dtcloud*0.5*ELI* &
                          lamda_f(:,:)*snow_cld(:,:)/ &
                          qa_mean(:,:)/rho_ice)
        END WHERE
        
!       Freezing of cloud liquid to cloud ice occurs when
!       the temperature is less than -40C. At these very cold temper-
!       atures it is assumed that homogenous freezing of cloud liquid
!       droplets will occur.   To accomplish this numerically in one 
!       time step:
!
!  (36) D*dtcloud =  ln( ql/qmin).
!
!       With this form it is guaranteed that if this is the only
!       process acting that ql = qmin after one integration.

        !do homogenous freezing
        WHERE (T(:,:,j)     .lt. (Tfreeze-40.) .and. &
               ql_mean(:,:) .gt. qmin .and. &
               qa_mean(:,:) .gt. qmin)
             D2_dt(:,:) = log ( ql_mean(:,:) / qmin )
        END WHERE
       
        
!       Analytic integration of ql equation
!
!       Following Tiedtke, the ql equation is written in the form
!
!  (37) dl/dt =    C  -  Dl
!
!       for which the analytic solution is
!
!  (38) l(t+dt) = l(t)*exp(-D*dt) + (C*dt/D*dt) * (1. - exp(-D*dt))
!
!       or
!
!  (39) l(t+dt)-l(t) = (1.-exp(-D*dt)) * ( (C*dt/D*dt) - l(t) )
!
!       Because of finite machine precision it is required that D*dt is
!       greater than a minimum value.  This minimum value occurs where
!       1. - exp(-D*dt) = 0. instead of D*dt.  This value will be
!       machine dependent. See discussion at top of code for Dmin.
!
!       Note that the amount of SL calculated here is stored in tmp1.

        !reset C_dt and D_dt from previous use; Source of
        !cloud liquid is initially set to 0. Sink of cloud
        !liquid is set to the sum of D1 (liquid to rain component), 
        !and D2 (liquid to ice component), and D_eros (erosion)
        C_dt(:,:) = 0.
        D_dt(:,:) = D1_dt(:,:) + D2_dt(:,:) + D_eros(:,:)
        
        !add in large_scale term,note use of ql mean
        WHERE (dcond_ls(:,:).ge.0.)
             C_dt(:,:)=C_dt(:,:)+dcond_ls(:,:)
        ELSEWHERE
             D_dt(:,:)=D_dt(:,:)-(dcond_ls(:,:)/ql_mean(:,:))
        END WHERE
        
        !do analytic integration
        WHERE (D_dt(:,:) .gt. Dmin)
             tmp1(:,:) = ( (C_dt(:,:)/D_dt(:,:)) - ql_upd(:,:) ) * &
                         ( 1.  - exp( -1.*D_dt(:,:) ) )
        ELSEWHERE
             tmp1 = C_dt(:,:)
        END WHERE
        
        !update SL and ql
        ql_upd(:,:) = ql_upd(:,:) + tmp1(:,:)
        SL(:,:,j) = SL(:,:,j) + tmp1(:,:)        
        
!       Apportion SL between various processes.  This is necessary to
!       account for how much the temperature changes due to various
!       phase changes.   For example:
!
!       liquid to ice =  (D2/D)*(C*dt-SL)
!
!       liquid to rain = (D1/D)*(C*dt-SL) (no phase change but needed to
!                                   know how much to increment rainflux)
!
!       vapor to liquid =-{ ((-dcond_ls/ql_mean)+D_eros)/ D }*(C*dt-SL)
!                         where dcond_ls < 0. but
!
!                         (dcond_ls/C_dt)*(C_dt) -(D_eros/D)*(C*dt-SL)                      
!                         where dcond_ls > 0.

        
        !do phase changes due to large-scale condensation
        WHERE (dcond_ls(:,:).ge.0.)
             ST(:,:,j) = ST(:,:,j) + HLv*dcond_ls(:,:)/Cp
             SQ(:,:,j) = SQ(:,:,j) - dcond_ls(:,:)
        END WHERE
        
        WHERE (dcond_ls(:,:) .lt. 0. .and. D_dt(:,:) .gt. Dmin)
             ST(:,:,j) = ST(:,:,j) + &
                         HLv*(dcond_ls(:,:)/ql_mean(:,:)/D_dt(:,:))* &
                         (C_dt(:,:) - tmp1(:,:))/Cp
             SQ(:,:,j) = SQ(:,:,j) - &
                         (dcond_ls(:,:)/ql_mean(:,:)/D_dt(:,:))* &
                         (C_dt(:,:) - tmp1(:,:))
        END WHERE
        
        
        WHERE (D_dt(:,:) .gt. Dmin)

             !add in liquid to ice and cloud erosion to temperature 
             !tendency
             ST(:,:,j) = ST(:,:,j) + &
                         (HLf*D2_dt(:,:)-HLv*D_eros(:,:))* &
                         (C_dt(:,:)-tmp1(:,:))/Cp/D_dt(:,:)

             !add conversion of liquid to rain to the rainflux
             rain_cld(:,:) = rain_cld(:,:) + &
                        (D1_dt(:,:)/D_dt(:,:))*(C_dt(:,:)-tmp1(:,:))* &
                        deltp(:,:)/Grav/dtcloud
        
             !cloud evaporation adds to water vapor
             SQ(:,:,j) = SQ(:,:,j) + &
                         D_eros(:,:)*(C_dt(:,:)-tmp1(:,:))/D_dt(:,:)
        
             !add liquid converted to ice into ice equation source
             C_dt(:,:) = (C_dt(:,:)-tmp1(:,:))*D2_dt(:,:)/D_dt(:,:)

             !increment qi_mean by the frozen liquid
             qi_mean(:,:) = qi_mean(:,:) + C_dt(:,:)

        ELSEWHERE  !needed to initialize ice source with zero condensate

             C_dt(:,:) = 0.     
             
        END WHERE 
        
        
        !auto conversion will change a_rain_cld upto area of cloud
        WHERE (D_dt(:,:) .gt. Dmin .and. D1_dt(:,:) .gt. 0.)
             a_rain_cld(:,:) = qa_mean(:,:)
        END WHERE
       
        
        
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
!       #42.
!
!       In terms of the analytic integration with respect qi of Tiedtke,
!       the flux in from the top of the grid layer is equated to the
!       source term, and the flux out of the bottom of the layer is 
!       equated to the sink term:
!
!  (40) C_dt =  snow_cld * dtcloud *Grav/deltp
!
!  (41) D_dt =  airdens*Grav*Vfall* dtcloud /deltp
!
!       All ice crystals are assumed to fall with the same fall speed
!       which is given as in Heymsfield and Donner (1990) as:
!
!  (42) Vfall = 3.29 * ( (airdens*qi_mean/qa_mean)**0.16)
!
!       which is the formula in Heymsfield and Donner.  Note however
!       that because this is uncertain, sensitivity runs will be made
!       with different formulations. Note that when Vfall is computed 
!       the source incremented qi, qi_mean, is used.  This gives some 
!       implicitness to the scheme.

        !compute Vfall
        Vfall(:,:) = 3.29*((airdens(:,:,j)*qi_mean(:,:)/ &
                              MAX(qa_mean(:,:),qmin))**0.16)

        !add to ice source the settling ice flux from above
        !note that because C_dt already contains the source
        !of liquid converted to ice from above
        C_dt(:,:) = C_dt(:,:) + snow_cld(:,:)*Grav*dtcloud/deltp(:,:)
        
        !Compute settling of ice. The result is multiplied by 
        !dtcloud/deltp to convert to units of D_dt.  
        !Note that if tracers are not advected then this is done
        !relative to the local vertical motion.
        if (tracer_advec) then
             tmp1(:,:) = 0.
        else
             tmp1(:,:) = omega(:,:,j)
        end if 
        
        WHERE (qi_mean(:,:) .gt. qmin .and. qa_mean(:,:) .gt. qmin)
             D1_dt(:,:)= MAX( 0., &
                         ((airdens(:,:,j)*Grav*Vfall(:,:))+tmp1(:,:)) &
                          * dtcloud/deltp(:,:)) 
        ELSEWHERE
             D1_dt(:,:)= 0.
        END WHERE 

        
!       Melting of in-cloud ice
!
!       Melting occurs where the temperature is greater than Tfreezing. 
!       This is an instaneous process such that no stratiform ice will
!       remain in the grid at the end of the timestep. 
!       No ice settles out of the grid box (i.e. D1 is set to zero) when
!       the amount of ice to be melted is less than that that would
!       bring the cloudy part of the grid box to the freezing point.
!
!       The ice that melts becomes rain.  This is because if an ice
!       crystal of dimension ~100 microns melts it will become a
!       droplet of size 50-100 microns which is clearly drizzle size
!       drop.  Ice crystals at temperatures near freezing are assumed
!       to be this large, consistent with the assumption of particle
!       size temperature dependence.

        !reset D2
        D2_dt(:,:) = 0.0
               
        !do melting in non-limited case
        WHERE ( T(:,:,j) .gt. Tfreeze .and. &
               qi_mean(:,:) .gt. qmin .and. &
               qa_mean(:,:) .gt. qmin .and. &
               (HLf*qi_mean(:,:)) .lt. (qa_mean(:,:)*Cp*(T(:,:,j)-Tfreeze)))
             D2_dt(:,:) = log ( qi_mean(:,:) / qmin )
             D1_dt(:,:) = 0.
        END WHERE
               
        !do melting in limited case
        WHERE ( T(:,:,j) .gt. Tfreeze .and. &
               qi_mean(:,:) .gt. qmin .and. &
               qa_mean(:,:) .gt. qmin .and. &
               (HLf*qi_mean(:,:)) .ge. (qa_mean(:,:)*Cp*(T(:,:,j)-Tfreeze)))
             D2_dt(:,:) = log ( qi_mean(:,:) / &
                (qi_mean(:,:)-(qa_mean(:,:)*Cp*(T(:,:,j)-Tfreeze)/HLf)) )
        END WHERE
        


!       Analytic integration of qi equation
!
!       Following Tiedtke, the qi equation is written in the form
!
!  (43) dl/dt =    C  -  Dl
!
!       for which the analytic solution is
!
!  (44) l(t+dt) = l(t)*exp(-D*dt) + (C*dt/D*dt) * (1. - exp(-D*dt))
!
!       or
!
!  (45) l(t+dt)-l(t) = (1.-exp(-D*dt)) * ( (C*dt/D*dt) - l(t) )
!
!       Because of finite machine precision it is required that D*dt is
!       greater than a minimum value.  This minimum value occurs where
!       1. - exp(-D*dt) = 0. instead of D*dt.  This value will be
!       machine dependent. See discussion at top of code for Dmin.


        !Compute D_dt which has contributions from D1_dt (ice settling)
        !and D2_dt(:,:) (ice melting), and D_eros(:,:) (cloud erosion). 
        !There is no calculation of C_dt done here because C_dt already 
        !includes the source of cloud ice falling from above as well as 
        !liquid converted to ice. 
        
        D_dt(:,:) =  D1_dt(:,:) + D2_dt(:,:) + D_eros(:,:)
        
        !add in large_scale term,note use of qi mean
        WHERE (dcond_ls_ice(:,:).ge.0.)
             C_dt(:,:)=C_dt(:,:)+dcond_ls_ice(:,:)
        ELSEWHERE
             D_dt(:,:)=D_dt(:,:)-(dcond_ls_ice(:,:)/qi_mean(:,:))
        END WHERE
        
        !do analytic integration
        WHERE (D_dt(:,:) .gt. Dmin)
            tmp1(:,:) = ( (C_dt(:,:)/D_dt(:,:)) - qi_upd(:,:) ) * &
                        ( 1.  - exp( -1.*D_dt(:,:) ) )
        ELSEWHERE
            tmp1(:,:) = C_dt(:,:)
        END WHERE
        

        !update SI and qi
        qi_upd(:,:) = qi_upd(:,:) + tmp1(:,:)
        SI(:,:,j) = SI(:,:,j) + tmp1(:,:)
        
!       Apportion SI between various processes.  This is necessary to
!       account for how much the temperature and water vapor changes 
!       due to various phase changes.   For example:
!
!       ice settling = (D1/D)*(C*dt-SI)*deltp/Grav/dtcloud
!                     
!       vapor to ice =
!           -{ ((-dcond_ls_ice/qi_mean) + D_eros)/ D }*(C*dt-SI)
!                         where dcond_ls < 0. but
!
!              (dcond_ls_ice/C_dt)*(C_dt) - (D_eros/D)*(C*dt-SI)
!                         where dcond_ls > 0.
!       
!       melting of ice = (D2/D)*(C*dt-SI)*deltp/Grav/dtcloud
!

        
        !do phase changes due to large-scale condensation
        WHERE (dcond_ls_ice(:,:) .ge. 0.)
             ST(:,:,j) = ST(:,:,j) + HLs*dcond_ls_ice(:,:)/Cp
             SQ(:,:,j) = SQ(:,:,j) - dcond_ls_ice(:,:)
        END WHERE
        
        WHERE (dcond_ls_ice(:,:) .lt. 0. .and. D_dt(:,:) .gt. Dmin)
             ST(:,:,j) = ST(:,:,j) + &
                        HLs*(dcond_ls_ice(:,:)/qi_mean(:,:)/D_dt(:,:))*&
                         (C_dt(:,:) - tmp1(:,:))/Cp
             SQ(:,:,j) = SQ(:,:,j) - &
                         (dcond_ls_ice(:,:)/qi_mean(:,:)/D_dt(:,:))* &
                         (C_dt(:,:) - tmp1(:,:))
        END WHERE
        
        !add settling ice flux to snow_cld and change a_snow_cld     
        WHERE (D_dt(:,:) .gt. Dmin .and. D1_dt(:,:) .gt. 0.)        

            !add settling ice flux to snow_cld 
            snow_cld(:,:) = D1_dt(:,:)* &
                            (C_dt(:,:)-tmp1(:,:))*deltp(:,:) / &
                            D_dt(:,:)/ Grav/ dtcloud
       
            a_snow_cld(:,:) = qa_mean(:,:)

        ELSEWHERE

            snow_cld(:,:) = 0.
            a_snow_cld(:,:) = 0.

        END WHERE
        
        WHERE (D_dt(:,:) .gt. Dmin .and. D2_dt(:,:) .gt. 0.)

             !add melting of ice to temperature tendency
             ST(:,:,j) = ST(:,:,j) - &
                         HLf*D2_dt(:,:)* &
                         (C_dt(:,:)-tmp1(:,:))/Cp/D_dt(:,:)

             !add melting of ice to the rainflux
             rain_cld(:,:) = rain_cld(:,:) + &
                        (D2_dt(:,:)/D_dt(:,:))*(C_dt(:,:)-tmp1(:,:))* &
                        deltp(:,:)/Grav/dtcloud
              
             !melting of ice creates area to a_rain_cld
             a_rain_cld(:,:) = qa_mean(:,:)

        END WHERE

        !Cloud evaporation changes temperature and vapor
        WHERE (D_dt(:,:) .gt. Dmin .and. D_eros(:,:) .gt. 0.) 
        
             ST(:,:,j) = ST(:,:,j) - HLs*D_eros(:,:)* &
                         (C_dt(:,:)-tmp1(:,:))/Cp/D_dt(:,:)
        
             SQ(:,:,j) = SQ(:,:,j) + D_eros(:,:)* &
                         (C_dt(:,:)-tmp1(:,:))/D_dt(:,:)
         
        END WHERE



        
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
!                              RAIN EVAPORATION
!                           
!
!       Rain evaporation is derived by integration of the growth 
!       equation of a droplet over the assumed Marshall-Palmer 
!       distribution of rain drops (equation #22).  This leads to the 
!       following formula:
!
!  (46) dqv/dt_local =  56788.636 * {rain_rate/Dens_h2o}^(11/18) *(1-U)/
!
!                      ( SQRT(airdens)* A_plus_B)
!
!       Numerically this equation integrated by use of time-centered
!       values for qs and qv.   This leads to the solution:
!
!  (47) qv_clr(t+1)-qv_clr(t) = K3 *[qs(t)-qv_clr(t)]/{1.+0.5*K3*gamma}
!
!       where 
!
!       K3= 56788.636 * dtcloud * {rain_rate_local/Dens_h2o}^(11/18) /
!           ( SQRT(airdens)* A_plus_B * qs)
!
!       and gamma is given by (3). Note that in (47), it is made 
!       explicit that it is the vapor concentration in the unsaturated 
!       part of the grid box that is used in the rain evaporation 
!       formula.
!
!       Now there are several limiters to this formula. First,
!       you cannot evaporate more than is available in a time step.
!       The amount available for evaporation locally is 
!       (rain_clr/a_rain_clr)*(Grav*dtcloud/deltp).   Second, to
!       avoid supersaturating the box or exceeding the critical
!       relative humidity above which rain does not evaporate, 
!       the amount of evaporation is limited to U_evap*qs(t)-qv_clr(t)/
!       gamma, where the factor of gamma is the appropriate reduction 
!       accounting for the cooling of the grid box while evaporation
!       occurs. 
!
!       Finally rain evaporation occurs ONLY if the relative humidity
!       in the unsaturated portion of the grid box, U_clr, is less
!       then a threshold, U_evap.   U_evap, will not necessarily be
!       one.   For example, stratiform precipitation in convective
!       regions rarely saturates subcloud air because of the presence
!       of downdrafts.  If the convection scheme does not have down-
!       drafts then it doesn't make sense to allow the sub-cloud layer
!       to saturate. U_clr may be solved from (10) as:
!
!  (48) U_clr = ( U - qa ) / (1. - qa)
!
!       Some variables are temporarily stored in tmp1.

        !compute U_clr
        U_clr(:,:) =  (U(:,:)-qa_mean(:,:))/MAX((1.-qa_mean(:,:)),qmin)
        
        !keep U_clr > 0. and U_clr < 1.
        U_clr(:,:) = MIN(MAX(U_clr(:,:),0.),1.)
        
        !compute K3
        tmp1(:,:) = 56788.636 * dtcloud * &
                  ((rain_clr(:,:)/MAX(a_rain_clr(:,:),qmin)/Dens_h2o)&
                    **(11./18.)) / &
                  SQRT(airdens(:,:,j))/A_plus_B(:,:,j)/qs(:,:,j)

        !compute local change in vapor mixing ratio due to 
        !rain evaporation
        tmp1(:,:) = tmp1(:,:) * qs(:,:,j) * (1. - U_clr(:,:)) / &
                  (1. + 0.5*tmp1(:,:)*gamma(:,:,j))

        !limit change in qv to the amount that would raise the relative
        !humidity to U_evap in the portion of the grid box that rain 
        !evaporation occurs in
        tmp1(:,:) = MIN(tmp1(:,:),&
                    qs(:,:,j)*MAX(0.,U_evap-U_clr(:,:))/gamma(:,:,j))
        
        !do limiter by amount available
        WHERE ((tmp1(:,:)*deltp(:,:)*a_rain_clr(:,:)/Grav/dtcloud).ge. &
               rain_clr(:,:) )
             
             SQ(:,:,j) = SQ(:,:,j)+rain_clr(:,:)*dtcloud*Grav/deltp(:,:)
             ST(:,:,j) = ST(:,:,j)-HLv*rain_clr(:,:)*dtcloud*Grav &
                                    /deltp(:,:)/Cp
             rain_clr(:,:) = 0.
             a_rain_clr(:,:) = 0.

        ELSEWHERE
          
             SQ(:,:,j) = SQ(:,:,j) + a_rain_clr(:,:)*tmp1(:,:)
             ST(:,:,j) = ST(:,:,j) - HLv*a_rain_clr(:,:)*tmp1(:,:)/Cp        
             rain_clr(:,:) = rain_clr(:,:) - &
                       a_rain_clr(:,:)*tmp1(:,:)*deltp(:,:)/Grav/dtcloud
        END WHERE
        
        

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
!  (49) dqi/dt_local = - a_snow_clr* K3 * (qs - qv_clr)
!
!       where the leading factor of a_snow_clr is the portion of the
!       grid box undergoing sublimation. K3 is given by
!
!  (50) K3 = (4/(pi*rho_air*qs*rho_ice*A_plus_B))*
!            ((snow_clr/a_snow_clr/3.29)**1.16 ) *
!           [ 0.65*lamda_f^2 + 
!             198.92227 * (airdens)^0.5 * 
!             ((snow_clr/a_snow_clr)**(1/14.5)) * lamda_f^(3/2) ]
!
!       Note that equation (49) is identical to equation (30) of 
!       Rotstayn.
!
!       Numerically this is integrated as in rain evaporation.


        !compute K3
        tmp1(:,:) = dtcloud * &
                 (4./3.14159/rho_ice/airdens(:,:,j)/&
                  A_plus_B(:,:,j)/qs(:,:,j))* &
                  ((snow_clr(:,:)/MAX(a_snow_clr(:,:),qmin)/3.29)&
                    **(1./1.16))* &
                  (0.65*lamda_f(:,:)*lamda_f(:,:) + &
                   198.92227*lamda_f(:,:)* &
                   SQRT(airdens(:,:,j)*lamda_f(:,:))*&
                  ( (snow_clr(:,:)/MAX(a_snow_clr(:,:),qmin))&
                    **(1./14.5) )  )

        !compute local change in vapor mixing ratio due to 
        !snow evaporation
        tmp1(:,:) = tmp1(:,:) * qs(:,:,j) * (1. - U_clr(:,:)) / &
                  (1. + 0.5*tmp1(:,:)*gamma(:,:,j))

        !limit change in qv to the amount that would raise the relative
        !humidity to U_evap in the portion of the grid box that snow 
        !evaporation occurs in
        tmp1(:,:) = MIN(tmp1(:,:),&
                    qs(:,:,j)*MAX(0.,U_evap-U_clr(:,:))/gamma(:,:,j))
        
        !do limiter by amount available
        WHERE ((tmp1(:,:)*deltp(:,:)*a_snow_clr(:,:)/Grav/dtcloud).ge. &
               snow_clr(:,:) )
             
             SQ(:,:,j) = SQ(:,:,j)+snow_clr(:,:)*dtcloud*Grav/deltp(:,:)
             ST(:,:,j) = ST(:,:,j)-HLs*snow_clr(:,:)*dtcloud*Grav &
                                    /deltp(:,:)/Cp
             snow_clr(:,:) = 0.
             a_snow_clr(:,:) = 0.

        ELSEWHERE
          
             SQ(:,:,j) = SQ(:,:,j) + a_snow_clr(:,:)*tmp1(:,:)
             ST(:,:,j) = ST(:,:,j) - HLs*a_snow_clr(:,:)*tmp1(:,:)/Cp        
             snow_clr(:,:) = snow_clr(:,:) - &
                       a_snow_clr(:,:)*tmp1(:,:)*deltp(:,:)/Grav/dtcloud
        END WHERE
        

!-----------------------------------------------------------------------
!
!       Adjustment.  Due to numerical errors in detrainment or advection
!       sometimes the current state of the grid box may be supersaturated.
!       Under the assumption that the temperature is constant in the grid
!       box and that q <= qs, the excess vapor is condensed and will be
!       put into the precipitation fluxes at the end of the loop
        
        !estimate current qs
        tmp2(:,:) = qs(:,:,j)+dqsdT(:,:,j)*ST(:,:,j)

        !compute excess over saturation
        WHERE ((qv(:,:,j)+SQ(:,:,j)) .gt. tmp2(:,:))
              tmp1(:,:) = (qv(:,:,j)+SQ(:,:,j)-tmp2(:,:))/gamma(:,:,j)
        ELSEWHERE
              tmp1(:,:) = 0.
        END WHERE

        !change vapor content
        SQ(:,:,j)=SQ(:,:,j)-tmp1(:,:)

        !add in excess to precipitation fluxes, change their area and 
        !increment temperature
        WHERE (T(:,:,j) .le. Tfreeze-20. .and. tmp1(:,:) .gt. 0.)
              snow_cld(:,:) = snow_cld(:,:) + &
                        qa_mean(:,:) *tmp1(:,:)*deltp(:,:)/Grav/dtcloud
              snow_clr(:,:) = snow_clr(:,:) + &
                    (1.-qa_mean(:,:))*tmp1(:,:)*deltp(:,:)/Grav/dtcloud
              a_snow_cld(:,:) = qa_mean(:,:)
              a_snow_clr(:,:) = 1.-qa_mean(:,:)
              ST(:,:,j) = ST(:,:,j) + HLs*tmp1(:,:)/Cp
        END WHERE
        WHERE (T(:,:,j) .gt. Tfreeze-20. .and. tmp1(:,:) .gt. 0.)
        
              rain_cld(:,:) = rain_cld(:,:) + &
                        qa_mean(:,:) *tmp1(:,:)*deltp(:,:)/Grav/dtcloud
              rain_clr(:,:) = rain_clr(:,:) + &
                    (1.-qa_mean(:,:))*tmp1(:,:)*deltp(:,:)/Grav/dtcloud
              a_rain_cld(:,:) = qa_mean(:,:)
              a_rain_clr(:,:) = 1.-qa_mean(:,:)
              ST(:,:,j) = ST(:,:,j) + HLv*tmp1(:,:)/Cp
        END WHERE
        
!-----------------------------------------------------------------------
!
!       Cloud Destruction occurs where both ql and qi are .le. qmin, 
!       or if qa is .le. qmin.
!       In this case all tracers are set to zero conserving moisture
!       and energy.

        WHERE ((ql_upd(:,:) .le. qmin .and. qi_upd(:,:) .le. qmin)&
               .or. (qa_upd(:,:) .le. qmin))
             SL(:,:,j) = SL(:,:,j) - ql_upd(:,:)
             SI(:,:,j) = SI(:,:,j) - qi_upd(:,:)
             SQ(:,:,j) = SQ(:,:,j) + ql_upd(:,:) + qi_upd(:,:)
             ST(:,:,j) = ST(:,:,j) - &
                         (HLv*ql_upd(:,:) + HLs*qi_upd(:,:))/Cp
             SA(:,:,j) = SA(:,:,j) - qa_upd(:,:)
        END WHERE
                
        
!-----------------------------------------------------------------------
!
!
!       Put rain and ice fluxes into surfrain and surfsnow if the
!       grid point is at the bottom of a column.   If MASK is not
!       present then this code is executed only if j .eq. KDIM.
!       IF MASK is present some grid points may be beneath ground. 
!       If a given grid point is at the bottom of the column then
!       the surface values of rain and snow must be created.
!       Also if the MASK is present then the code forces all tenden-
!       cies below ground to be zero. Note that MASK = 1. equals above
!       ground point, MASK = 0. equals below ground point.

        IF (present(MASK)) THEN

           !zero out all tendencies below ground
           ST(:,:,j)=MASK(:,:,j)*ST(:,:,j)
           SQ(:,:,j)=MASK(:,:,j)*SQ(:,:,j)
           SL(:,:,j)=MASK(:,:,j)*SL(:,:,j)
           SI(:,:,j)=MASK(:,:,j)*SI(:,:,j)
           SA(:,:,j)=MASK(:,:,j)*SA(:,:,j)

           IF (j .lt. KDIM) THEN
                  
                  !bottom of true points in columns which contain some
                  !dummy points
                  WHERE(MASK(:,:,j) .eq. 1. .and. MASK(:,:,j+1) .eq. 0.)
                          surfrain(:,:) = dtcloud* &
                                           (rain_clr(:,:)+rain_cld(:,:))
                          surfsnow(:,:) = dtcloud* &
                                           (snow_clr(:,:)+snow_cld(:,:))
                          rain_clr(:,:) = 0.
                          rain_cld(:,:) = 0.
                          snow_clr(:,:) = 0.
                          snow_cld(:,:) = 0.
                          a_rain_clr(:,:) = 0.
                          a_rain_cld(:,:) = 0.
                          a_snow_clr(:,:) = 0.
                          a_snow_cld(:,:) = 0.
                  END WHERE

           ELSE

                  !bottom of column for those columns which contain no
                  !dummy points
                  WHERE(MASK(:,:,j) .eq. 1.)
                          surfrain(:,:) = dtcloud* &
                                           (rain_clr(:,:)+rain_cld(:,:))
                          surfsnow(:,:) = dtcloud* &
                                           (snow_clr(:,:)+snow_cld(:,:))
                          rain_clr(:,:) = 0.
                          rain_cld(:,:) = 0.
                          snow_clr(:,:) = 0.
                          snow_cld(:,:) = 0.
                          a_rain_clr(:,:) = 0.
                          a_rain_cld(:,:) = 0.
                          a_snow_clr(:,:) = 0.
                          a_snow_cld(:,:) = 0.                  
                  END WHERE

           END IF

       ELSE

           !do code if we are at bottom of column
           IF (j .eq. KDIM) THEN
                  surfrain(:,:) = dtcloud*(rain_clr(:,:)+rain_cld(:,:))
                  surfsnow(:,:) = dtcloud*(snow_clr(:,:)+snow_cld(:,:))
           END IF

       END IF 
                  

!-----------------------------------------------------------------------
!
!       END LOOP OVER VERTICAL LEVELS
!

        ENDDO


!-----------------------------------------------------------------------
!
!       CLEAR OUT RADTURBTEN FOR FUTURE USE
!

        radturbten(is:ie,js:je,:) = 0.


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
                unit = Open_File ('strat.data', form='ieee32',&
                                  action='append')
                write (unit) ipt,jpt,     ql(i,j,:)+SL(i,j,:)
                write (unit) ipt,jpt,     qi(i,j,:)+SI(i,j,:)
                write (unit) ipt,jpt,     qa(i,j,:)+SA(i,j,:)
                write (unit) ipt,jpt,      T(i,j,:)+ST(i,j,:) 
                write (unit) ipt,jpt,     qv(i,j,:)+SQ(i,j,:)
                write (unit) ipt,jpt,     pfull(i,j,:)
                close (unit)
         endif
       enddo
     endif


!-----------------------------------------------------------------------
!
!
!       END OF SUBROUTINE
!
!


END SUBROUTINE STRAT_DRIV



!#######################################################################
!#######################################################################


SUBROUTINE ADD_STRAT_TEND(is,ie,js,je,tend)


!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This subroutine adds the fields tend to the radturbten
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!	VARIABLES
!
!
!       -------------------
!	INTERNAL VARIABLES:
!       -------------------
!	
!         variable              definition                  unit
!       ------------   -----------------------------   ---------------
!
!       is,ie          starting and ending i indices 
!                      for data window
!
!       js,je          starting and ending j indices 
!                      for data window
!
!       tend           tendency of physical field      K/sec
!
!       
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!        
!       User Interface variables
!       ------------------------
!

        INTEGER, INTENT (IN)                   :: is,ie,js,je
        REAL, INTENT(IN), DIMENSION(:,:,:)     :: tend
        


!-----------------------------------------------------------------------
!       
!       Code
!

        if (.not. strat_cloud_on) return

        radturbten(is:ie,js:je,:)=radturbten(is:ie,js:je,:)+tend(:,:,:)

       
END SUBROUTINE ADD_STRAT_TEND


!#######################################################################
!#######################################################################


SUBROUTINE SUBTRACT_STRAT_TEND(is,ie,js,je,tend)


!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This subroutine subtracts the field tend from the radturbten
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!	VARIABLES
!
!
!       -------------------
!	INTERNAL VARIABLES:
!       -------------------
!	
!         variable              definition                  unit
!       ------------   -----------------------------   ---------------
!
!       is,ie          starting and ending i indices 
!                      for data window
!
!       js,je          starting and ending j indices 
!                      for data window
!
!       tend           tendency of physical field      K/sec
!
!       
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!        
!       User Interface variables
!       ------------------------
!

        INTEGER, INTENT (IN)                   :: is,ie,js,je
        REAL, INTENT(IN), DIMENSION(:,:,:)     :: tend
        


!-----------------------------------------------------------------------
!       
!       Code
!

        if (.not. strat_cloud_on) return

        radturbten(is:ie,js:je,:)=radturbten(is:ie,js:je,:)-tend(:,:,:)

       
END SUBROUTINE SUBTRACT_STRAT_TEND


!#######################################################################
!#######################################################################


SUBROUTINE STRAT_END()


!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This subroutine writes out radturbten to a restart file.
!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!	VARIABLES
!
!
!       -------------------
!	INTERNAL VARIABLES:
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

        INTEGER                                :: unit


!-----------------------------------------------------------------------
!       
!       Code
!

        unit = Open_File ('RESTART/strat_cloud.res', &
                          FORM='native', ACTION='write')
        call write_data (unit, radturbten)
        call write_data (unit, nsum)
        call write_data (unit, qlsum)
        call write_data (unit, qisum)
        call write_data (unit, cfsum)
        Call Close_File (unit)
       

END SUBROUTINE STRAT_END


!#######################################################################

 SUBROUTINE STRAT_CLOUD_SUM (is, js, ql, qi, cf)

!-----------------------------------------------------------------------
   integer, intent(in)                   :: is, js
      real, intent(in), dimension(:,:,:) :: ql, qi, cf
!-----------------------------------------------------------------------
   integer :: ie, je

   ie = is + size(ql,1) - 1
   je = js + size(ql,2) - 1
     
!--------- use time-averaged or instantaneous clouds -----------

  if (do_average) then
       nsum(is:ie,js:je)   =  nsum(is:ie,js:je)   +  1
      qlsum(is:ie,js:je,:) = qlsum(is:ie,js:je,:) + ql(:,:,:)
      qisum(is:ie,js:je,:) = qisum(is:ie,js:je,:) + qi(:,:,:)
      cfsum(is:ie,js:je,:) = cfsum(is:ie,js:je,:) + cf(:,:,:)
  else
       nsum(is:ie,js:je)   =  1
      qlsum(is:ie,js:je,:) = ql(:,:,:)
      qisum(is:ie,js:je,:) = qi(:,:,:)
      cfsum(is:ie,js:je,:) = cf(:,:,:)
  endif

!-----------------------------------------------------------------------

 END SUBROUTINE STRAT_CLOUD_SUM

!#######################################################################

 SUBROUTINE STRAT_CLOUD_AVG (is, js, ql, qi, cf, ierr)

!-----------------------------------------------------------------------
   integer, intent(in)                    :: is, js
      real, intent(out), dimension(:,:,:) :: ql, qi, cf
   integer, intent(out)                   :: ierr
!-----------------------------------------------------------------------
   integer :: ie, je, num, k
!-----------------------------------------------------------------------
      
   if (size(ql,3) .ne. size(qlsum,3)) call error_mesg ( &
                              'strat_cloud_avg in strat_cloud_mod',  &
                              'input argument has the wrong size',FATAL)

   ie = is + size(ql,1) - 1
   je = js + size(ql,2) - 1
   num = count(nsum(is:ie,js:je) == 0)

   if (num > 0) then

!     ----- no average, return error flag -----

!!!   call error_mesg ('strat_cloud_avg in strat_cloud_mod',  &
!!!                    'dividing by a zero counter.', FATAL)
      ierr = 1

   else

!     ----- compute average -----

      do k = 1, size(ql,3)
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

 END SUBROUTINE STRAT_CLOUD_AVG

!#######################################################################

 FUNCTION DO_STRAT_CLOUD ( ) RESULT (answer)
  logical :: answer

  answer = strat_cloud_on

 END FUNCTION DO_STRAT_CLOUD

!#######################################################################
!#######################################################################


END MODULE STRAT_CLOUD_MOD

                
