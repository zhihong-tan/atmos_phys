
MODULE CLOUD_RAD_MOD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	CLOUD RADIATIVE PROPERTIES
!
!       April 1999
!       Contact person: Steve Klein
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This module solves for the radiative properties of
!       every cloud.  In particular it uses either the
!       two stream approach or the delta-Eddington approach
!       to solve for the longwave emmissivity, the ultra-violet-
!       visible reflectivities and absorptions, and the
!       near-infrared reflectivities and absorptions.
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
!            CLOUD_SUMMARY
!                        for each cloud it computes the cloud amount,
!                        the liquid and ice water paths, the effective
!                        particle sizes, the tops and bottoms of clouds
!                        and then calls CLOUD_OPTICAL_PROPERTIES to 
!                        retrieve the optical properties of each cloud.
!                        Finally it calls CLOUD_RAD to determine the 
!                        radiative properties of each cloud given its
!                        optical properties.
!    	     CLOUD_RAD   solves for the radiative properties of the
!                        clouds given the cloud optical properties
!                        (tau,w0,gg) for each cloud using either a 
!                        Delta-Eddington solution (default) or the
!                        two stream approximation.
!            CLOUD_OPTICAL_PROPERTIES
!                        for each cloud it calculates the mixed phase
!                        values of the optical depth, the single scattering
!                        albedo, and the asymmetry parameter (tau, w0,and g)
!                        for the visible and near infrared bands.  It also
!                        computes the longwave emmissivity of each cloud.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!       Declare outside modules used, make implicit nones, declare
!       public routines/variables.
!

        USE  Constants_Mod,      ONLY :  Rdgas,Grav,Tfreeze,Dens_h2o
        USE  Utilities_Mod,      ONLY :  File_Exist, Open_File,  &
                                         error_mesg, FATAL,     &
                                         Close_File, get_my_pe
          
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
!       qmin           minimum permissible cloud       kg condensate/kg air                  
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
!                      overlap = 2. means condensate in 
!                                   adjacent levels is 
!                                   treated as different 
!                                   clouds
!
!                      overlap is used only by the stratiform
!                      cloud scheme.
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
REAL,    PRIVATE            :: N_land = 6.E+08
REAL,    PRIVATE            :: N_ocean = 1.E+08
REAL,    PRIVATE, PARAMETER :: k_land = 1.143
REAL,    PRIVATE, PARAMETER :: k_ocean = 1.077

INTEGER, PRIVATE            :: overlap = 1
LOGICAL, PRIVATE            :: l2strem = .FALSE.
REAL,    PRIVATE            :: taucrit = 1.
LOGICAL, PRIVATE            :: adjust_top = .TRUE.
REAL,    PRIVATE            :: scale_factor = 0.8
REAL,    PRIVATE            :: qamin = 1.E-2


!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CREATE NAMELIST
!

        NAMELIST /CLOUD_RAD_NML/ overlap,l2strem,taucrit,adjust_top, &
                                 scale_factor,qamin

!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       DECLARE VERSION NUMBER OF SCHEME
!
        
        character(len=128) :: version = '$Id: cloud_rad.F90,v 1.2 2000/08/04 19:30:28 fms Exp $'
        character(len=128) :: tag = '$Name: bombay $'

! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CONTAINS


!#######################################################################
!#######################################################################


SUBROUTINE CLOUD_RAD_INIT(qmin_in,N_land_in,N_ocean_in)
                               

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
!       -----
!       INPUT
!       -----
!
!	
!         variable              definition                  unit
!       ------------   -----------------------------   ---------------
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

REAL,     INTENT (IN)                    :: qmin_in,N_land_in,&
                                            N_ocean_in

!  Internal variables
!  ------------------


INTEGER                                  :: unit,io

!-----------------------------------------------------------------------
!       
!       Code
!

!-----------------------------------------------------------------------
!
!       Assign values

        qmin = qmin_in
        N_land = N_land_in
        N_ocean = N_ocean_in
        
!-----------------------------------------------------------------------
!
!       Namelist functions

        !read namelist if it exists
        If (File_Exist('input.nml')) Then
             unit = Open_File ('input.nml', action='read')
             io=1
             Do While (io .ne. 0)
                   Read  (unit, nml=CLOUD_RAD_NML, iostat=io, End=10)
             EndDo
  10         Call Close_File (unit)
        EndIf

        !write namelist variables to logfile
        unit = Open_File ('logfile.out', action='APPEND')
        if ( get_my_pe() == 0 ) then
             Write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
             Write (unit,nml=CLOUD_RAD_NML)
        endif
        Call Close_File (unit)

!-----------------------------------------------------------------------
!
!
!       END OF SUBROUTINE
!
!


END SUBROUTINE CLOUD_RAD_INIT


!########################################################################
!########################################################################

SUBROUTINE CLOUD_SUMMARY(LAND,ql,qi,qa,pfull,phalf,TKel,coszen,&
                  nclds,ktop,kbot,cldamt,&
                  r_uv,r_nir,ab_uv,ab_nir,em_lw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine returns the following properties of clouds
!
!               1. nclds: # of clouds
!               2. ktop : integer level for top of cloud
!               3. kbot : integer level for bottom of cloud
!               4. cldamt:horizontal cloud amount of every cloud
!               5. r_uv : cloud reflectance in uv band
!               6. r_nir: cloud reflectance in nir band
!               7. ab_uv: cloud absorption in uv band
!               8. ab_nir:cloud absorption in nir band
!               9. em_lw :longwave cloud emmissivity
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
!       r_uv         cloud reflectance in uv band
!       r_nir        cloud reflectance in nir band
!       ab_uv        cloud absorption in uv band
!       ab_nir       cloud absorption in nir band
!       em_lw        longwave cloud emmissivity
!
!       -------------------
!	INTERNAL VARIABLES:
!       -------------------
!
!       i,j,k,t      looping variables
!       IDIM         number of first dimension points
!       JDIM         number of second dimension points
!       KDIM         number of vertical levels
!
!       N_drop       number of cloud droplets per unit volume
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
!
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
!     FOR ICE CLOUDS THE EFFECTIVE RADIUS IS TAKEN FROM THE FORMULATION
!     IN DONNER (1997, J. Geophys. Res., 102, pp. 21745-21768) WHICH IS
!     BASED ON HEYMSFIELD AND PLATT (1984) WITH ENHANCEMENT FOR PARTICLES
!     SMALLER THAN 20 MICRONS.  
!
!              T Range (K)                            Reff (microns)
!     ---------------------------------------       -----------------
!
!     Tfreeze-25. < T                                    92.46298
!     Tfreeze-30. < T <= Tfreeze-25.                     72.35392 
!     Tfreeze-35. < T <= Tfreeze-30.                     85.19071 
!     Tfreeze-40. < T <= Tfreeze-35.                     55.65818
!     Tfreeze-45. < T <= Tfreeze-40.                     35.29989
!     Tfreeze-50. < T <= Tfreeze-45.                     32.89967
!     Tfreeze-55  < T <= Tfreeze-50                      16.60895
!                   T <= Tfreeze-55.                     15.41627
!
!        IN THIS PROGRAM, reff_ice is limited to be between 10 microns
!        and 130 microns, which is the range of validity for the Ebert
!        and Curry (1992) radiation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

REAL,     INTENT (IN), DIMENSION(:,:)    :: LAND
REAL,     INTENT (IN), DIMENSION(:,:)    :: coszen
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: ql,qi,qa,pfull,TKel
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: phalf
INTEGER,  INTENT (INOUT),DIMENSION(:,:)  :: nclds
INTEGER,  INTENT (INOUT),DIMENSION(:,:,:):: ktop,kbot
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: cldamt
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: r_uv,r_nir,ab_uv,ab_nir,em_lw

!  Internal variables
!  ------------------

INTEGER                                  :: i,j,k,IDIM,JDIM,KDIM
INTEGER                                  :: t,top_t,bot_t
INTEGER                                  :: tmp_top,tmp_bot,nlev,max_cld
LOGICAL                                  :: add_cld
REAL                                     :: sum_liq,sum_ice,maxcldfrac
REAL                                     :: totcld_bot,max_bot
REAL                                     :: totcld_top,max_top,tmp_val
REAL                                     :: reff_liq_local,sum_reff_liq
REAL                                     :: reff_ice_local,sum_reff_ice
REAL, DIMENSION(SIZE(ql,1),SIZE(ql,2),SIZE(ql,3)) :: qa_local,ql_local,qi_local
REAL, DIMENSION(SIZE(ql,1),SIZE(ql,2),SIZE(ql,3)) :: LWP,IWP,Reff_liq,Reff_ice
REAL, DIMENSION(SIZE(ql,1),SIZE(ql,2))   :: N_drop, k_ratio
REAL, DIMENSION(SIZE(ql,1),SIZE(ql,2),SIZE(ql,3),4) :: tau,w0,gg

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
        r_uv(:,:,:)   = 0.
        r_nir(:,:,:)  = 0.
        ab_uv(:,:,:)  = 0.
        ab_nir(:,:,:) = 0.
        em_lw(:,:,:)  = 0.
        LWP(:,:,:)    = 0.
        IWP(:,:,:)    = 0.
        Reff_liq(:,:,:) = 10.
        Reff_ice(:,:,:) = 30.
        tau(:,:,:,:)    = 0.
        w0(:,:,:,:)     = 0.
        gg(:,:,:,:)     = 0.

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


    !-----------  DETERMINE CLOUD AMOUNT, LWP, IWP FOR EACH CLOUD ------!

    if (overlap .eq. 1) then


        !---- DO CONDENSING OF CLOUDS ----!

        !initialize variables
        sum_liq          = 0.
        sum_ice          = 0.
        maxcldfrac       = 0.
        sum_reff_liq     = 0.
        sum_reff_ice     = 0.
        reff_liq_local   = 10.
        reff_ice_local   = 50.
        add_cld          = .FALSE.

        !---loop over vertical levels----!
        DO i = 1, IDIM
        DO j = 1, JDIM
        DO k = 1, KDIM

                 !identify new cloud tops
                 IF ( ( (ql_local(i,j,k) .gt. qmin) .or. &
                           (qi_local(i,j,k) .gt. qmin) ) .and. &
                           (.NOT. add_cld) ) then
                       nclds(i,j) = nclds(i,j) + 1
                       add_cld = .TRUE.
                       ktop(i,j,nclds(i,j)) = k
                 END IF

                 !increment sums where cloud
                 IF (   (ql_local(i,j,k) .gt. qmin) .or. &
                        (qi_local(i,j,k) .gt. qmin)  ) then

                       !compute reff
                       if (ql_local(i,j,k) .gt. qmin) then
                           reff_liq_local = k_ratio(i,j)* 620350.49 * &
                                 (  pfull(i,j,k)*  &
                        ql_local(i,j,k)/qa_local(i,j,k)/Rdgas/&
                           TKel(i,j,k)/Dens_h2o/N_drop(i,j))**(1./3.)
                       end if
                       !limit reff_liq_local to values for which
                       !Slingo radiation is valid :
                       ! 4.2 microns < reff < 16.6 microns
                       reff_liq_local = MIN(16.6,reff_liq_local)
                       reff_liq_local = MAX(4.2, reff_liq_local)

                       if (qi_local(i,j,k) .gt. qmin) then
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
                       end if

                       !increment sums
                       sum_liq = sum_liq + ql_local(i,j,k)* &
                                      (phalf(i,j,k+1)-phalf(i,j,k))/Grav
                       sum_ice = sum_ice + qi_local(i,j,k)* &
                                      (phalf(i,j,k+1)-phalf(i,j,k))/Grav
                       maxcldfrac = MAX(maxcldfrac,qa_local(i,j,k))
                       sum_reff_liq  = sum_reff_liq + &
                           ( reff_liq_local * ql_local(i,j,k) * &
                           (phalf(i,j,k+1)-phalf(i,j,k))/Grav )
                       sum_reff_ice  = sum_reff_ice + &
                           ( reff_ice_local * qi_local(i,j,k) * &
                           (phalf(i,j,k+1)-phalf(i,j,k))/Grav )

                 END IF


                 !where the first cloud gap exists after a cloud
                 ! or bottom level is reached compute kbot, cldamt,
                 ! LWP, IWP, Reff_liq, and Reff_ice
                 IF (  ( (ql_local(i,j,k) .le. qmin) .and. &
                           (qi_local(i,j,k) .le. qmin) .and. &
                           (add_cld) ) .or. &
                        (add_cld .and. k .eq. KDIM)) then

                        !determine kbot
                        kbot(i,j,nclds(i,j))= k-1
                        if ((ql_local(i,j,k) .gt. qmin) .or. &
                            (qi_local(i,j,k) .gt. qmin)) then
                              kbot(i,j,nclds(i,j)) = k
                        end if

                        cldamt(i,j,nclds(i,j)) = maxcldfrac
                        LWP(i,j,nclds(i,j)) = sum_liq / &
                               cldamt(i,j,nclds(i,j))
                        IWP(i,j,nclds(i,j)) = sum_ice / &
                               cldamt(i,j,nclds(i,j))
                        if (sum_liq .gt. 0.) then
                             Reff_liq(i,j,nclds(i,j)) = &
                             sum_reff_liq / sum_liq
                        end if
                        if (sum_ice .gt. 0.) then
                             Reff_ice(i,j,nclds(i,j)) = &
                             sum_reff_ice / sum_ice
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
                             
                              tmp_val = MAX(0.,qa_local(i,j,top_t)-totcld_top)
                              if (tmp_val .gt. max_top) then
                                max_top = tmp_val
                                tmp_top = top_t
                              end if
                              totcld_top = totcld_top+tmp_val         
                              
                              tmp_val = MAX(0.,qa_local(i,j,bot_t)-totcld_bot)
                              if (tmp_val .gt. max_bot) then
                                max_bot = tmp_val
                                tmp_bot = bot_t
                              end if
                              totcld_bot = totcld_bot+tmp_val         
                               
                              END DO
                       
                              !assign tmp_top and tmp_bot to ktop and kbot
                              ktop(i,j,nclds(i,j)) = tmp_top
                              kbot(i,j,nclds(i,j)) = tmp_bot

                        end if


                        !reset values
                        add_cld     = .FALSE.
                        sum_liq     = 0.
                        sum_ice     = 0.
                        maxcldfrac  = 0.
                        sum_reff_liq = 0.
                        sum_reff_ice = 0.
                        reff_liq_local = 10.
                        reff_ice_local = 50.


                 END IF


        END DO
        END DO
        END DO

    else if (overlap .eq. 2) then

           !---loop over vertical levels----!
           DO i = 1, IDIM
           DO j = 1, JDIM
           DO k = 1, KDIM

               !where cloud exists compute ktop,kbot, cldamt and LWP and IWP
                 IF (   (ql_local(i,j,k) .gt. qmin) .or. &
                           (qi_local(i,j,k) .gt. qmin)  ) then
                        nclds(i,j) = nclds(i,j) + 1
                        ktop(i,j,nclds(i,j)) = k
                        kbot(i,j,nclds(i,j)) = k

                        cldamt(i,j,nclds(i,j)) = qa_local(i,j,k)

                        LWP(i,j,nclds(i,j)) = ql_local(i,j,k)* &
                                      (phalf(i,j,k+1)-phalf(i,j,k))/Grav/&
                                       cldamt(i,j,nclds(i,j))
                        IWP(i,j,nclds(i,j)) = qi_local(i,j,k)* &
                                      (phalf(i,j,k+1)-phalf(i,j,k))/Grav/&
                                       cldamt(i,j,nclds(i,j))

                       !compute reff_liqiud
                       if (ql_local(i,j,k) .gt. qmin) then
                           reff_liq_local = k_ratio(i,j)* 620350.49 * &
                                 ( pfull(i,j,k)*  &
                        ql_local(i,j,k)/qa_local(i,j,k)/Rdgas/&
                           TKel(i,j,k)/Dens_h2o/N_drop(i,j))**(1./3.)
                       end if
                       !limit reff_liq_local to values for which
                       !Slingo radiation is valid :
                       ! 4.2 microns < reff < 16.6 microns
                       if (ql_local(i,j,k) .gt. qmin) then
                             reff_liq_local = MIN(16.6,reff_liq_local)
                             Reff_liq(i,j,nclds(i,j)) = &
                                              MAX(4.2, reff_liq_local)
                       end if

                       !compute reff_ice
                       if (qi_local(i,j,k) .gt. qmin) then
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
                          Reff_ice(i,j,nclds(i,j)) = &
                                           MAX(10.,reff_ice_local)
                       end if

                 END IF

           END DO
           END DO
           END DO

    end if

    !find maximum number of clouds
    max_cld  = MAXVAL(nclds(:,:))

    !compute cloud radiative properties
    IF (max_cld .gt. 0) then

         ! compute cloud optical properties
         CALL CLOUD_OPTICAL_PROPERTIES(LWP(:,:,1:max_cld),IWP(:,:,1:max_cld),&
                  Reff_liq(:,:,1:max_cld),Reff_ice(:,:,1:max_cld),&
                  tau(:,:,1:max_cld,:),w0(:,:,1:max_cld,:),gg(:,:,1:max_cld,:),&
                  em_lw(:,:,1:max_cld))
         
         !Account for plane-parallel homogenous cloud bias
         tau(:,:,:,:) = scale_factor * tau(:,:,:,:)

         !compute cloud radiative properties
         CALL CLOUD_RAD(tau(:,:,1:max_cld,:),w0(:,:,1:max_cld,:),&
                  gg(:,:,1:max_cld,:),coszen,&
                  r_uv(:,:,1:max_cld),r_nir(:,:,1:max_cld),&
                  ab_uv(:,:,1:max_cld),ab_nir(:,:,1:max_cld))

         !assure that zero clouds have properties of zero clouds
         DO i = 1, IDIM
         DO j = 1, JDIM
                  if (nclds(i,j).lt.max_cld) then
                          r_uv(i,j,nclds(i,j)+1:max_cld)   = 0.
                          r_nir(i,j,nclds(i,j)+1:max_cld)  = 0.
                          ab_uv(i,j,nclds(i,j)+1:max_cld)  = 0.
                          ab_nir(i,j,nclds(i,j)+1:max_cld) = 0.
                          em_lw(i,j,nclds(i,j)+1:max_cld)  = 0.
                  end if
         ENDDO
         ENDDO

    END IF
    
END SUBROUTINE CLOUD_SUMMARY


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
!       l2strem      logical variable indicating 2 stream operating or not
!                          l2strem = T  2 stream solution to radiation
!                          l2strem = F  Delta-Eddington solution to radiation
!
!                            IF l2strem = T then the solution doesnot
!                            depend on solar zenith angle
!
!       taucrit      critical tau for switching direct beam to diffuse beam
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
                        tau,w0,gg,em_lw)


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

END SUBROUTINE CLOUD_OPTICAL_PROPERTIES


!########################################################################
!########################################################################


END MODULE CLOUD_RAD_MOD

