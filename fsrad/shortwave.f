
                        MODULE SHORTWAVE_MOD

C-----------------------------------------------------------------------

      USE  RDPARM_MOD, ONLY: LMAX,LP1,LLP1,LP2,LLP2,NB
      USE  HCONST_MOD, ONLY: DIFFCTR,GINV,O3DIFCTR,RADCON

      Use  Utilities_Mod, ONLY: Error_Mesg, FATAL

CCC   PRIVATE   LMAX,LP1,LLP1,NB
      PRIVATE   LMAX,LP1,LLP1,LP2,LLP2,NB
      PRIVATE   DIFFCTR,GINV,O3DIFCTR,RADCON
      Private   Error_Mesg

      INTEGER, PRIVATE :: IMAX

C------- interfaces -------
      PUBLIC  SWRAD
      PRIVATE SWRAD_ORIG


      CONTAINS

C#######################################################################

      SUBROUTINE SWRAD (NCLDS,KTOPSW,KBTMSW,PRESS,RH2O,QO3,CAMT,
     &                  CUVRF,CIRRF,CIRAB,RRCO2,COSZRO,SSOLAR,
     &                  SALB, FSW,DFSW,UFSW,HSW, LSFC,PSFC)

C-----------------------------------------------------------------------
C              WRAPPER FOR  SHORT WAVE RADIATION CODE
C     inserts surface albedo into appropriate cloud property arrays
C-----------------------------------------------------------------------

      INTEGER, INTENT(IN), DIMENSION(:,:)   :: NCLDS
      INTEGER, INTENT(IN), DIMENSION(:,:,:) :: KTOPSW,KBTMSW
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: PRESS,RH2O,QO3
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: CAMT,CUVRF,CIRRF,CIRAB
      REAL,    INTENT(IN)                   :: RRCO2
      REAL,    INTENT(IN), DIMENSION(:,:)   :: COSZRO,SSOLAR
      REAL,    INTENT(IN), DIMENSION(:,:)   :: SALB

      REAL,   INTENT(OUT), DIMENSION(:,:,:) :: FSW,DFSW,UFSW,HSW

      INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:,:)   :: LSFC
         REAL, INTENT(IN), OPTIONAL, DIMENSION(:,:)   :: PSFC

C-----------------------------------------------------------------------
      REAL,DIMENSION(SIZE(CUVRF,1),SIZE(CUVRF,3)) :: CUVRF2,CIRRF2
      REAL,DIMENSION(SIZE(CAMT ,1),SIZE(CAMT ,3)) :: CAMT2

      INTEGER  i,j
C-----------------------------------------------------------------------

      IMAX=SIZE(PRESS,1)

      DO j=1,SIZE(PRESS,2)

      CUVRF2(:,:)=CUVRF(:,j,:)
      CIRRF2(:,:)=CIRRF(:,j,:)
      CAMT2 (:,:)=CAMT (:,j,:)

      DO i=1,IMAX
         IF (COSZRO(i,j) > 0.0) THEN
            CUVRF2(i,NCLDS(i,j)+2)=SALB(i,j)
            CIRRF2(i,NCLDS(i,j)+2)=SALB(i,j)
            CAMT2 (i,NCLDS(i,j)+2)=1.0
         ENDIF
      ENDDO

C----- check usage of optional arguments -------
      IOPT=0
      IF (PRESENT(LSFC)) IOPT=IOPT+1
      IF (PRESENT(PSFC)) IOPT=IOPT+2

C     ------------------
      SELECT CASE (IOPT)
C     ------------------
          CASE (0)
C     ------------------
      CALL SWRAD_ORIG (NCLDS(:,j),KTOPSW(:,j,:),KBTMSW(:,j,:),
     &                 PRESS(:,j,:),RH2O(:,j,:),QO3(:,j,:),CAMT2(:,:),
     &                 CUVRF2,CIRRF2,CIRAB(:,j,:),RRCO2,
     &                 COSZRO(:,j),SSOLAR(:,j),
     &                 FSW(:,j,:),DFSW(:,j,:),UFSW(:,j,:),HSW(:,j,:))
C     ------------------
          CASE (3)
C     ------------------
      CALL SWRAD_ORIG (NCLDS(:,j),KTOPSW(:,j,:),KBTMSW(:,j,:),
     &                 PRESS(:,j,:),RH2O(:,j,:),QO3(:,j,:),CAMT2(:,:),
     &                 CUVRF2,CIRRF2,CIRAB(:,j,:),RRCO2,
     &                 COSZRO(:,j),SSOLAR(:,j),
     &                 FSW(:,j,:),DFSW(:,j,:),UFSW(:,j,:),HSW(:,j,:),
     &                 LSFC(:,j),PSFC(:,j))
C     ------------------
          CASE DEFAULT
C     ------------------
             Call Error_Mesg ('SWRAD in SHORTWAVE_MOD',
     &                    'LSFC and PSFC must be used together.', FATAL)
C     ------------------
      END SELECT
C     ------------------

      ENDDO
C-----------------------------------------------------------------------

      END SUBROUTINE SWRAD
      
C#######################################################################

      SUBROUTINE SWRAD_ORIG (NCLDS,KTOPSW,KBTMSW,PRESS,RH2O,QO3,CAMT,
     &                  CUVRF,CIRRF,CIRAB,RRCO2,COSZRO,SSOLAR,
     &                  FSW,DFSW,UFSW,HSW,LSFC,PSFC)

C***********************************************************************
C                    SHORT WAVE RADIATION CODE
C***********************************************************************
C-----------------------------------------------------------------------
C                        INPUT PARAMETERS
C                        ----------------
C
C      NCLDS   =  NO. CLOUDS AT EACH GRID PT. 
C      KTOPSW  =  INDEX OF (FLUX LEVEL) PRESSURE OF CLOUD TOP, USED 
C                    IN THE SHORTWAVE PROGRAM 
C      KBTMSW  =  INDEX OF (FLUX LEVEL) PRESSURE OF CLOUD BOTTOM, 
C                    USED IN THE SHORTWAVE PROGRAM
C      PRESS   =  PRESSURE (CGS UNITS) AT DATA LEVELS OF MODEL
C      RH2O    =  MASS MIXING RATIO (G/G) OF H2O AT MODEL DATA LVLS.
C      QO3     =  MASS MIXING RATIO (G/G) OF O3 AT MODEL DATA LVLS. 
C      CAMT    =  CLOUD AMOUNTS OF CLOUDS (THEIR LOCATIONS ARE
C                    SPECIFIED IN THE KTOP/KBTM INDICES)
C      CUVRF   =  REFLECTIVITY OF CLOUDS IN THE VISIBLE FREQ. BAND
C                    USED IN SHORTWAVE CALCS. ONLY
C      CIRRF   =  REFLECTIVITY OF CLOUDS IN THE INFRARED FREQ. BAND 
C                    USED IN SHORTWAVE CALCS. ONLY
C      CIRAB   =  ABSORPTIVITY OF CLOUDS IN THE INFRARED FREQ. BAND 
C                    USED IN SHORTWAVE CALCS. ONLY
C      RRCO2   =  MASS MIXING RATIO (G/G) OF CO2,USED IN SHORTWAVE
C                    CALCS. ONLY (scalar)
C      COSZRO  =  ZENITH ANGLE AT GRID PT. USED ON SHORTWAVE CALCS. 
C      SSOLAR  =  TOTAL SOLAR FLUX EITHER FOR THE TIMESTEP OR AVERAGED
C                 OVER THE DAY OR YEAR,
C                 EQUALS THE SOLAR CONSTANT x NORMALIZED SOLAR FLUX
C                 (INCLUDES THE COS ZENITH ANGLE AND DAY FRACTION)
C                 (AT PRESENT,IN LY/MIN).
C
C      LSFC    =  Vertical index of the lowest model level,
C                    dimensioned by IMAX.
C      PSFC    =  Surface pressure
C
C-----------------------------------------------------------------------
C                        OUTPUT PARAMETERS
C                        -----------------
C
C      FSW     = NET RADIATION (UP-DOWN) IN CGS UNITS AT ALL
C                PRESSURE LEVELS
C     DFSW     = DOWNWARD RADIATION AT ALL PRESSURE LEVELS
C     UFSW     = UPWARD RADIATION AT ALL PRESSURE LEVELS
C      HSW     = SHORTWAVE HEATING RATES IN K/DAY FOR PRESSURE
C                LAYERS.
C
C-----------------------------------------------------------------------

      INTEGER, INTENT(IN), DIMENSION(:)   :: NCLDS
      INTEGER, INTENT(IN), DIMENSION(:,:) :: KTOPSW,KBTMSW
      REAL,    INTENT(IN), DIMENSION(:,:) :: PRESS,RH2O,QO3
      REAL,    INTENT(IN), DIMENSION(:,:) :: CAMT,CUVRF,CIRRF,CIRAB
      REAL,    INTENT(IN)                 :: RRCO2
      REAL,    INTENT(IN), DIMENSION(:)   :: COSZRO,SSOLAR

      REAL,   INTENT(OUT), DIMENSION(:,:) :: FSW,DFSW,UFSW,HSW

      INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:)   :: LSFC
      REAL,    INTENT(IN), OPTIONAL, DIMENSION(:)   :: PSFC
       
C-----------------------------------------------------------------------
C                        D I M E N S I O N
C    &  NCLDS(IMAX),KTOPSW(IMAX,LP1),KBTMSW(IMAX,LP1)
C    2 ,PRESS(IMAX,LP1),RH2O(IMAX,LMAX),QO3(IMAX,LMAX)
C    3 ,CAMT(IMAX,LP1),CUVRF(IMAX,LP1),CIRRF(IMAX,LP1),CIRAB(IMAX,LP1)
C    4 ,COSZRO(IMAX),SSOLAR(IMAX),LSFC(IMAX),PSFC(IMAX)
C                        D I M E N S I O N
C    &  FSW(IMAX,LP1),DFSW(IMAX,LP1),UFSW(IMAX,LP1),HSW(IMAX,LMAX)
C-----------------------------------------------------------------------
C***********************************************************************
      LOGICAL BCLDS,BJTOP
C-----------------------------------------------------------------------
                         D I M E N S I O N
     &  BCLDS(IMAX,LP1),BJTOP(IMAX,LP1),
     1  ICNT(LP1),ICNT1(LP1),IINCL(LP1),INDX4(IMAX,LP1),INDXK(IMAX),
     2  IBETCL(LP1) 
C-----------------------------------------------------------------------
                         D I M E N S I O N
     &  DFN(IMAX,LP1,NB),UFN(IMAX,LP1,NB),
     1  TTD(IMAX,LP1,NB),TTU(IMAX,LP1,NB),
     2  PP    (IMAX,LP1),PPTOP (IMAX,LP1),
     &  DP    (IMAX,LP1),DPCLD (IMAX,LP1),
     3  CR    (IMAX,LP1),CT    (IMAX,LP1),
     4  TDCL1 (IMAX,LP1),TDCL2 (IMAX,LP1),TUCL1 (IMAX,LP1),
     5  TUCL1I(IMAX,LP1),TDCL1I(IMAX,LP1),TDCL2I(IMAX,LP1),
     6  TCLU  (IMAX,LP1),TCLD  (IMAX,LP1),ALFAU (IMAX,LP1),
     7  UFNCLU(IMAX,LP1),UFNCLD(IMAX,LP1),
     8  DFNCLU(IMAX,LP1),DFNCLD(IMAX,LP1),
     9  TEMP1(IMAX),TEMP2(IMAX),TEMP3(IMAX),TEMP4(IMAX),
     A  TEMP5(IMAX),TEMP6(IMAX),ALFA (IMAX),
     B  VV   (IMAX),REFL (IMAX),SECZ (IMAX),RRAY (IMAX)
C***********************************************************************
C-------DIMENSION OF VARIABLES EQUIVALENCED TO THOSE PREVIOUSLY---------
C***********************************************************************
C-----------------------------------------------------------------------
C                        D I M E N S I O N
      DIMENSION PPBOT(IMAX,LP1)
      DIMENSION DFNTOP(IMAX,NB)
      DIMENSION UD(IMAX,LP1),UR(IMAX,LP1)
C     DIMENSION UCO2(IMAX,LLP2),UO3(IMAX,LLP2),ACO2(IMAX,LLP2)
C     DIMENSION AO3(IMAX,LLP2)
      DIMENSION VTROW1(IMAX,LP1)
      DIMENSION VTROW2(IMAX,LP1),VTROW3(IMAX,LP1)
      DIMENSION FF(IMAX,LP1),FFCO2(IMAX,LP1),FFO3(IMAX,LP1)
      DIMENSION PR2(IMAX,LP1)
      DIMENSION DU(IMAX,LP1),DUCO2(IMAX,LP1),DUO3(IMAX,LP1)
      DIMENSION ADCO2(IMAX,LP1),AUCO2(IMAX,LP1),UDCO2(IMAX,LP1),
     &          URCO2(IMAX,LP1)
      DIMENSION ABSDO3(IMAX,LP1),ABSUO3(IMAX,LP1),UDO3(IMAX,LP1),
     &          URO3(IMAX,LP1)
      DIMENSION CR1D(IMAX*LP1),ALFU1D(IMAX*LP1),TCLU1D(IMAX*LP1),
     &          TCLD1D(IMAX*LP1)
C--DIMENSIONS OF LOCAL DATA VARIABLES---
      DIMENSION ABCFF(NB),PWTS(NB)
C***********************************************************************
C-----------------------------------------------------------------------
C     EQUIVALENCE (ADCO2(1,1),ACO2(1,1)),(AUCO2(1,1),ACO2(1,LP2))
C     EQUIVALENCE (UDCO2(1,1),UCO2(1,1)),(URCO2(1,1),UCO2(1,LP2))
C     EQUIVALENCE (ABSDO3(1,1),AO3(1,1)),(ABSUO3(1,1),AO3(1,LP2))
C     EQUIVALENCE (UDO3(1,1),UO3(1,1)),(URO3(1,1),UO3(1,LP2))

C     EQUIVALENCE (CR,UD),(CT,UR)
C     EQUIVALENCE (TDCL1I,ABSDO3,PPBOT),(TDCL2I,ABSUO3)
C     EQUIVALENCE (UFNCLU,UDO3),(UFNCLD,URO3)
C     EQUIVALENCE (DFNCLU,UDCO2),(DFNCLD,URCO2)
C     EQUIVALENCE (TCLU,ADCO2),(TCLD,AUCO2)
C     EQUIVALENCE (TTD(1,1,1),FF(1,1))
C     EQUIVALENCE (TTD(1,1,2),FFCO2(1,1))
C     EQUIVALENCE (TTD(1,1,3),PR2(1,1))
C     EQUIVALENCE (TTD(1,1,4),DU(1,1))
C     EQUIVALENCE (TTD(1,1,5),DUCO2(1,1))
C     EQUIVALENCE (TTD(1,1,6),DUO3(1,1))
C     EQUIVALENCE (TTD(1,1,7),VTROW1(1,1))
C     EQUIVALENCE (TTD(1,1,8),VTROW2(1,1))
C     EQUIVALENCE (TTD(1,1,9),VTROW3(1,1))
C     EQUIVALENCE (TTU(1,1,1),DFNTOP(1,1))
C     EQUIVALENCE (CR1D,CR),(ALFU1D,ALFAU),(TCLD1D,TCLD),(TCLU1D,TCLU)
C-----------------------------------------------------------------------
      DATA ABCFF /2*4.0E-5,.002,.035,.377,1.95,9.40,44.6,190./
      DATA PWTS /.5000,.1470,.0698,.1443,.0584,.0335,.0225,.0158,.0087/
      DATA CFCO2,CFO3 /508.96,466.64/
      DATA REFLO3 /1.9/
      DATA RRAYAV /0.144/
C-----------------------------------------------------------------------

      DO K=1,LP1
      DO I=1,IMAX
         FF   (I,K)=DIFFCTR
         FFCO2(I,K)=DIFFCTR
         FFO3 (I,K)=O3DIFCTR
      ENDDO
      ENDDO

C------ NOTE: converting pressures (PRESS) to CGS units -------

      DO K=2,LMAX
      DO I=1,IMAX
CCCC     PP(I,K)=0.50*(PRESS(I,K)+PRESS(I,K-1))
         PP(I,K)=5.0*(PRESS(I,K)+PRESS(I,K-1))
      ENDDO
      ENDDO

      DO I=1,IMAX
         SECZ(I)=35./SQRT(1224.*COSZRO(I)*COSZRO(I)+1.0)
C        SECZ(I)=1./COSZRO(I)
      ENDDO

      IF (.not.PRESENT(LSFC)) THEN
         DO I=1,IMAX
            PP(I,1)=0.00
            PP(I,LP1)=10.*PRESS(I,LP1)
         ENDDO
      ELSE
         DO I=1,IMAX
            PP(I,1)=0.00
            PP(I,LP1)=10.*PRESS(I,LP1)
            PP(I,LSFC(I)+1)=10.*PSFC(I)
         ENDDO
      ENDIF

      DO K=1,LMAX
      DO I=1,IMAX
         DP(I,K)=PP(I,K+1)-PP(I,K)
      ENDDO
      ENDDO

      IF (.not.PRESENT(LSFC)) THEN
         DO K=1,LMAX
         DO I=1,IMAX
CCCC        PR2(I,K)=0.50*(PP(I,K)+PP(I,K+1))/PRESS(I,LP1)
            PR2(I,K)=0.050*(PP(I,K)+PP(I,K+1))/PRESS(I,LP1)
         ENDDO
         ENDDO
      ELSE
         DO K=1,LMAX
         DO I=1,IMAX
CCCC        PR2(I,K)=0.50*(PP(I,K)+PP(I,K+1))/PSFC(I)
            PR2(I,K)=0.050*(PP(I,K)+PP(I,K+1))/PSFC(I)
         ENDDO
         ENDDO
      ENDIF

      DO K=1,LMAX
      DO I=1,IMAX
         DUO3 (I,K)=QO3 (I,K)*DP(I,K)*GINV
         DUCO2(I,K)=RRCO2    *DP(I,K)*GINV
         DU   (I,K)=RH2O(I,K)*DP(I,K)*GINV
      ENDDO
      ENDDO

      DO I=1,IMAX
         JTOP=KTOPSW(I,2) 
         DO K=1,JTOP 
            FFO3 (I,K)=SECZ(I) 
            FFCO2(I,K)=SECZ(I)
            FF   (I,K)=SECZ(I) 
         ENDDO
      ENDDO

C-----------------------------------------------------------------------
C     CALCULATE PRESSURE-WEIGHTED OPTICAL PATHS IN UNITS OF G/CM2.
C     PRESSURE WEIGHTING IS BY P**0.5 
C     UD IS THE DOWNWARD PATH,UR THE UPWARD PATH, 
C     AND THE CALCULATION IS MADE BY TAKING A PATH WITH AN ANGLE
C     OF (SECZ) FROM THE TOP OF THE ATMOSPHERE TO THE TOPMOST CLOUD,
C     THEN USING THE DIFFUSIVITY FACTOR (1.66) TO THE SURFACE AND FOR 
C     REFLECTED RADIATION. THE CODE BELOW REFLECTS THIS.
C-----------------------------------------------------------------------
C*****************************************
      IF (.not.PRESENT(LSFC)) THEN
         DO K=1,LMAX
         DO I=1,IMAX
            VTROW1(I,K)=DU   (I,K)*PR2(I,K)
            VTROW2(I,K)=DUCO2(I,K)*PR2(I,K)*CFCO2
            VTROW3(I,K)=DUO3 (I,K)*CFO3
         ENDDO
         ENDDO
      ELSE
         DO K=1,LMAX
         DO I=1,IMAX
            VTROW1(I,K)=0.00
            VTROW2(I,K)=0.00
            VTROW3(I,K)=0.00
         ENDDO
         ENDDO
         DO I=1,IMAX
            LMA=LSFC(I)
            DO K=1,LMA
               VTROW1(I,K)=DU   (I,K)*PR2(I,K)
               VTROW2(I,K)=DUCO2(I,K)*PR2(I,K)*CFCO2
               VTROW3(I,K)=DUO3 (I,K)*CFO3
            ENDDO
         ENDDO
      ENDIF
C*****************************************

      DO I=1,IMAX
         UD   (I,1)=0.00 
         UDCO2(I,1)=0.00
         UDO3 (I,1)=0.00 
      ENDDO

      DO K=2,LP1
C     DO I=1,IMAX
         UD   (:,K)=UD   (:,K-1)+VTROW1(:,K-1)*FF   (:,K) 
         UDCO2(:,K)=UDCO2(:,K-1)+VTROW2(:,K-1)*FFCO2(:,K)
         UDO3 (:,K)=UDO3 (:,K-1)+VTROW3(:,K-1)*FFO3 (:,K) 
C     ENDDO
      ENDDO

C   UDO3,URO3 ARE IN UNITS OF CM 
C   CFCO2,CFO3 IS THE CONVERSION FACTOR FROM GM/CM2 TO CM

      DO I=1,IMAX
         UR   (I,LP1)=UD   (I,LP1) 
         URCO2(I,LP1)=UDCO2(I,LP1) 
         URO3 (I,LP1)=UDO3 (I,LP1) 
      ENDDO

      DO K=LMAX,1,-1 
      DO I=1,IMAX
         UR   (I,K)=UR   (I,K+1)+VTROW1(I,K)*DIFFCTR
         URCO2(I,K)=URCO2(I,K+1)+VTROW2(I,K)*DIFFCTR
         URO3 (I,K)=URO3 (I,K+1)+VTROW3(I,K)*REFLO3 
      ENDDO
      ENDDO

C   CALCULATE WATER VAPOR TRANSMISSION FUNCTIONS FOR BANDS 2-9;
C   T.F. FOR BAND 1= T.F FOR BAND 2

      DO N=2,NB 
      DO K=1,LP1
      DO I=1,IMAX
         TTD(I,K,N)=ABCFF(N)*UD(I,K) 
         TTU(I,K,N)=ABCFF(N)*UR(I,K) 
      ENDDO
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LP1
      DO I=1,IMAX
         IF (TTD(I,K,N).GE.50.)  TTD(I,K,N)=50.
         IF (TTU(I,K,N).GE.50.)  TTU(I,K,N)=50.
      ENDDO
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LP1
      DO I=1,IMAX
         TTD(I,K,N)=-1.0*TTD(I,K,N)
         TTU(I,K,N)=-1.0*TTU(I,K,N)
      ENDDO
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LP1
      DO I=1,IMAX
         TTD(I,K,N)=EXP(TTD(I,K,N))
         TTU(I,K,N)=EXP(TTU(I,K,N))
      ENDDO
      ENDDO
      ENDDO

C   CALCULATE CO2 ABSORPTIONS . THEY WILL BE USED IN BANDS 2-9. 
C   SINCE THESE OCCUPY 50 PERCENT OF THE SOLAR SPECTRUM THE 
C   ABSORPTIONS WILL BE MULTIPLIED BY 2 

      DO I=1,IMAX
         ADCO2(I,1)=0.00
      ENDDO

      DO K=2,LP1
      DO I=1,IMAX
         ADCO2(I,K)=UDCO2(I,K)+0.0129
         ADCO2(I,K)=0.26*LOG(ADCO2(I,K))
         ADCO2(I,K)=EXP(ADCO2(I,K))
         ADCO2(I,K)=2.35E-3*ADCO2(I,K)-7.58265E-4 
      ENDDO
      ENDDO

      DO K=1,LMAX
      DO I=1,IMAX
         AUCO2(I,K)=URCO2(I,K)+0.0129
         AUCO2(I,K)=0.26*LOG(AUCO2(I,K))
         AUCO2(I,K)=EXP(AUCO2(I,K))
         AUCO2(I,K)=2.35E-3*AUCO2(I,K)-7.58265E-4 
      ENDDO
      ENDDO

C     DO K=2,LLP1
C     DO I=1,IMAX
C        ACO2(I,K)=UCO2(I,K)+0.0129
C     ENDDO
C     ENDDO

C     DO K=2,LLP1
C     DO I=1,IMAX
C        ACO2(I,K)=LOG(ACO2(I,K))
C     ENDDO
C     ENDDO

C     DO K=2,LLP1
C     DO I=1,IMAX
C        ACO2(I,K)=0.26*ACO2(I,K)
C     ENDDO
C     ENDDO

C     DO K=2,LLP1
C     DO I=1,IMAX
C        ACO2(I,K)=EXP(ACO2(I,K))
C     ENDDO
C     ENDDO

C     DO K=2,LLP1
C     DO I=1,IMAX
C        ACO2(I,K)=2.35E-3*ACO2(I,K)-7.58265E-4 
C     ENDDO
C     ENDDO

      DO I=1,IMAX
         AUCO2(I,LP1)=ADCO2(I,LP1) 
      ENDDO

      DO K=1,LP1 
      DO I=1,IMAX 
         ADCO2(I,K)=2.0*ADCO2(I,K) 
         AUCO2(I,K)=2.0*AUCO2(I,K) 
      ENDDO
      ENDDO

C   NOW CALCULATE OZONE ABSORPTIONS. THESE WILL BE USED IN
C   BAND 1. AS THIS OCCUPIES 50 PERCENT OF THE SOLAR SPECTRUM 
C   THE OZONE ABSORPTIONS WILL BE MULTIPLIED BY 2

      DO I=1,IMAX
         ABSDO3(I,1)=0.00 
      ENDDO

      H103P6=103.6*103.6*103.6

      DO K=2,LP1
      DO I=1,IMAX
         ABSDO3(I,K)=1.0+138.6*UDO3(I,K) 
         ABSDO3(I,K)=-0.805*LOG(ABSDO3(I,K))
         ABSDO3(I,K)=EXP(ABSDO3(I,K))
         ABSDO3(I,K)=1.082*UDO3(I,K)*ABSDO3(I,K)+ 
     &      0.0658*UDO3(I,K)/(1.0+H103P6*UDO3(I,K)*UDO3(I,K)*UDO3(I,K))+
     &      0.02118*UDO3(I,K)/(1.0+0.042*UDO3(I,K)+ 
     &      0.000323*UDO3(I,K)*UDO3(I,K))
      ENDDO
      ENDDO

      DO K=1,LMAX
      DO I=1,IMAX
         ABSUO3(I,K)=1.0+138.6*URO3(I,K) 
         ABSUO3(I,K)=-0.805*LOG(ABSUO3(I,K))
         ABSUO3(I,K)=EXP(ABSUO3(I,K))
         ABSUO3(I,K)=1.082*URO3(I,K)*ABSUO3(I,K)+ 
     &      0.0658*URO3(I,K)/(1.0+H103P6*URO3(I,K)*URO3(I,K)*URO3(I,K))+
     &      0.02118*URO3(I,K)/(1.0+0.042*URO3(I,K)+ 
     &      0.000323*URO3(I,K)*URO3(I,K))
      ENDDO
      ENDDO

C     DO K=2,LLP1
C     DO I=1,IMAX
C        AO3(I,K)=1.0+138.6*UO3(I,K) 
C     ENDDO
C     ENDDO

C     DO K=2,LLP1
C     DO I=1,IMAX
C        AO3(I,K)=LOG(AO3(I,K))
C     ENDDO
C     ENDDO

C     DO K=2,LLP1
C     DO I=1,IMAX
C        AO3(I,K)=-0.805*AO3(I,K)
C     ENDDO
C     ENDDO

C     DO K=2,LLP1
C     DO I=1,IMAX
C        AO3(I,K)=EXP(AO3(I,K)) 
C     ENDDO
C     ENDDO

C     H103P6=103.6*103.6*103.6
C     DO K=2,LLP1
C     DO I=1,IMAX 
C        AO3(I,K)=1.082*UO3(I,K)*AO3(I,K)+ 
C    &        0.0658*UO3(I,K)/(1.0+H103P6*UO3(I,K)*UO3(I,K)*UO3(I,K))+
C    &        0.02118*UO3(I,K)/(1.0+0.042*UO3(I,K)+ 
C    &        0.000323*UO3(I,K)*UO3(I,K))
C     ENDDO
C     ENDDO

      DO I=1,IMAX
         ABSUO3(I,LP1)=ABSDO3(I,LP1) 
      ENDDO

C     DO K=1,LLP2
      DO K=1,LP1
      DO I=1,IMAX 
C        AO3(I,K)=2.0*AO3(I,K)
         ABSDO3(I,K)=2.0*ABSDO3(I,K) 
         ABSUO3(I,K)=2.0*ABSUO3(I,K) 
      ENDDO
      ENDDO

C     WRITE (*,101) ((K,UD(IP,K),UR(IP,K),K=1,LP1),IP=1,IMAX) 
C     WRITE (*,105) ((K,UDO3(IP,K),URO3(IP,K),ABSDO3(IP,K), 
C    1 ABSUO3(IP,K),K=1,LP1),IP=1,IMAX) 
C     WRITE (*,105) ((K,UDCO2(IP,K),URCO2(IP,K),ADCO2(IP,K),
C    1 AUCO2(IP,K),K=1,LP1),IP=1,IMAX)

C   COMBINE ABSORPTIONS AND TRANSMISSIONS TO OBTAIN A 
C   TRANSMISSION FUNCTION FOR EACH OF THE 9 BANDS.

      DO K=1,LP1 
      DO I=1,IMAX 
         TTD(I,K,1)=TTD(I,K,2)*(1.0-ABSDO3(I,K))
      ENDDO
      ENDDO

      DO K=1,LMAX 
      DO I=1,IMAX 
         TTU(I,K,1)=TTU(I,K,2)*(1.0-ABSUO3(I,K))
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LP1 
      DO I=1,IMAX 
         TTD(I,K,N)=TTD(I,K,N)*(1.0-ADCO2(I,K)) 
      ENDDO
      ENDDO
      ENDDO

      DO N=1,NB
      DO I=1,IMAX 
         TTU(I,LP1,N)=TTD(I,LP1,N) 
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LMAX 
      DO I=1,IMAX 
         TTU(I,K,N)=TTU(I,K,N)*(1.0-AUCO2(I,K)) 
      ENDDO
      ENDDO
      ENDDO

C   IN THE 850 LOOP BELOW AND THE 855 LOOP,WE WILL SCALE DFN(IP,1,N)
C   TO UNITY. AFTER THE 850 LOOP,WE MULTIPLY BY DFN(IP,1,N) FROM THE
C   6 LOOP TO GET THE ACTUAL DFN'S AND UFN'S.
C   THE 855 LOOP=SCALED DOWNWARD FLUX ADOVE TOPMOST CLOUD(REDUN-
C   DANTLY OBTAINED TO THE GROUND) .ALSO,UFN IS INITIALIZED TO 0

      DO N=1,NB
      DO K=1,LP1
      DO I=1,IMAX
         DFN(I,K,N)=TTD(I,K,N)
         UFN(I,K,N)=0.00
      ENDDO
      ENDDO
      ENDDO
 
C***EVALUATE LOGICAL ARRAYS USED FOR CLOUD CALCULATIONS
C----BCLDS : TRUE IF KK IS < OR = TO NCLDS (NO. CLOUDS AT GRID PT I)
C----BJTOP : TRUE IF K IS < OR = TO KTOPSW(I,2) (INDEX OF TOP CLOUD,
C                    IF ANY; LP1 IS NO CLOUD AT GRID PT)
 
      DO KK=1,LP1
      DO I=1,IMAX
         BCLDS(I,KK)=KK.LE.NCLDS(I)
      ENDDO
      ENDDO

      DO K=1,LP1
      DO I=1,IMAX
         BJTOP(I,K)=K.LE.KTOPSW(I,2)
      ENDDO
      ENDDO

C---COUNT NO. OF PTS IN EACH ROW FOR WHICH BCLDS IS TRUE (ICNT)
C   AND FOR WHICH BJTOP IS TRUE (ICNT1)
 
      DO K=1,LP1
         ICNT(K)=0
         ICNT1(K)=0
         DO I=1,IMAX
            IF (BCLDS(I,K)) THEN
               ICNT(K)=ICNT(K)+1
            ENDIF
            IF (BJTOP(I,K)) THEN
               ICNT1(K)=ICNT1(K)+1
            ENDIF
         ENDDO
      ENDDO

C---FIND NO. OF CLOUD LEVELS WITH NONZERO VALUES OF ICNT
C---FIND NO. OF PRESSURE LEVELS WITH NONZERO VALUES OF ICNT1
 
      KCLDS=0
      KJTOP=0
      DO K=1,LP1
         IF (ICNT(K).GT.0) THEN
            KCLDS=KCLDS+1 
         ENDIF
         IF (ICNT1(K).GT.0) THEN
            KJTOP=KJTOP+1
         ENDIF
      ENDDO

C***IF NO CLOUDS AT ALL EXIST IN THE ROW, THE CALCULATIONS ARE
C   DRASTICALLY SIMPLIFIED.

CMPP  IF (KCLDS.EQ.0) THEN
      IF (KCLDS.EQ.-1) THEN
         DO N=1,NB

            IF (N.EQ.1) THEN
               DO I=1,IMAX
                  REFL(I)=CUVRF(I,2)
CCCCC             REFL(I)=SALB(I)
                  RRAY(I)=0.219/(1.0+0.816*COSZRO(I))
                  REFL(I)=RRAY(I)+(1.0-RRAY(I))*(1.0-RRAYAV)*REFL(I)/
     &                    (1.0-REFL(I)*RRAYAV)
                  ALFA(I)=REFL(I)
               ENDDO
            ELSE
               DO I=1,IMAX
                  ALFA(I)=CIRRF(I,2)
CCCCC             ALFA(I)=SALB(I)
               ENDDO
            ENDIF

            DO I=1,IMAX
               VV(I)=ALFA(I)*DFN(I,LP1,N)/TTU(I,LP1,N)
            ENDDO

            DO K=1,LP1
            DO I=1,IMAX
               UFN(I,K,N)=VV(I)*TTU(I,K,N)
            ENDDO
            ENDDO

         ENDDO
      ENDIF

C***********************************************************************
C     ****** COMPUTE NORMAL CASE: AT LEAST 1 PT HAS A CLOUD ******

CMPP                  IF (KCLDS.NE.0) THEN
                      IF (KCLDS.GE.0) THEN

C***********************************************************************

C---FIND  HIGHEST PRESSURE LEVEL WITH AT LEAST 1 PT BELOW TOP CLOUD
      KCLDS2=0
      DO K=1,LP1
         IF (ICNT1(K).EQ.IMAX) THEN
            KCLDS2=KCLDS2+1
         ENDIF
      ENDDO
      KCLDS2=KCLDS2+1

C    -------------------
      DO 2105 KK=1,KCLDS
C    -------------------
C---DETERMINE WHETHER A CLOUD LAYER KK HAS AT LEAST 1 GRID PT WITH A
C   "THICK CLOUD"
C---DETERMINE WHETHER THERE IS AT LEAST 1 GRID PT WHERE THE BOTTOM
C  OF CLOUD LAYER KK DOES NOT COINCIDE WITH THE TOP OF CLOUD LAYER
C  KK+1

      IINCL(KK)=0
      IBETCL(KK)=0
      DO K=KCLDS2,LP1
      DO I=1,IMAX
         IF (BCLDS(I,KK) .AND.
     &     (K.GT.KTOPSW(I,KK+1) .AND. K.LT.KBTMSW(I,KK+1))) THEN
               IINCL(KK)=IINCL(KK)+1
         ENDIF
         IF (BCLDS(I,KK) .AND. 
     &     (K.GE.KBTMSW(I,KK+1) .AND. K.LE.KTOPSW(I,KK+2))) THEN
               IBETCL(KK)=IBETCL(KK)+1
         ENDIF
      ENDDO
      ENDDO
C    -------------------
2105  CONTINUE
C    -------------------

C***COMPUTE VISIBLE BAND GROUND REFLECTIVITY USING LACIS-HANSEN 
C   PARAMETERIZATION

      DO I=1,IMAX
         REFL(I)=CUVRF(I,NCLDS(I)+2)
CCCCC    REFL(I)=SALB(I)
      ENDDO

      DO IP=1,IMAX
         RRAY(IP)=0.219/(1.0+0.816*COSZRO(IP))
         REFL(IP)=RRAY(IP)+(1.0-RRAY(IP))*(1.0-RRAYAV)*REFL(IP)/
     &                     (1.0-REFL(IP)*RRAYAV)
      ENDDO

      DO KK=1,KCLDS
      DO I=1,IMAX
         PPTOP(I,KK)=PP(I,KTOPSW(I,KK+1)) 
         PPBOT(I,KK)=PP(I,KBTMSW(I,KK+1))
      ENDDO
      ENDDO

      DO KK=1,KCLDS
      DO I=1,IMAX
         IF (PPTOP(I,KK).NE.PPBOT(I,KK)) THEN
            DPCLD(I,KK)=1.0/(PPTOP(I,KK)-PPBOT(I,KK))
         ELSE
            DPCLD(I,KK)=0.00
         ENDIF
      ENDDO
      ENDDO

C***WE NOW OBTAIN AN INDEX FOR (I,NCLDS(I)+1-KK).WE FORCE THIS
C   INDEX TO HAVE A MINIMUM VALUE IN THE (0,IMAX) RANGE.

      DO KK=1,KCLDS+1
      DO I=1,IMAX
         IF (BCLDS(I,KK)) THEN
            INDXK(I)=KK
         ELSE
            INDXK(I)=NCLDS(I)
         ENDIF
         INDX4(I,KK)=(NCLDS(I)-INDXK(I))*IMAX+I
      ENDDO
      ENDDO

C-----------------------------------------------------------------------
C   THE REST OF THE CLOUD CALCULATION IS PERFORMED INSIDE A
C   BAND (FREQUENCY) LOOP OVER N, RUNNING FROM 1 TO NB
C-----------------------------------------------------------------------

                         DO 2301 N=1,NB

C     print *, 'NB,N=',NB,N
C-----------------------------------------------------------------------

C***INITIALIZE CR TO ZERO AND CT TO ONE***
      DO K=1,LP1
      DO I=1,IMAX
         CR(I,K)=0.00
         CT(I,K)=1.0
      ENDDO
      ENDDO

C***OBTAIN CLOUD REFLECTION AND TRANSMISSION COEFFICIENTS, FIRST FOR
C   VISIBLE BAND (N=1) THEN FOR NEAR IR BANDS (N=2-NB) 
C---FIRST, THE VISIBLE BAND:

C               ----------------
                IF (N == 1) THEN
C               ----------------

      DO KK=1,KCLDS
      DO I=1,IMAX
         IF (BCLDS(I,KK)) THEN
            CR(I,KK+1)=CUVRF(I,KK+1)*CAMT(I,KK+1)
            CT(I,KK+1)=1.0-CR(I,KK+1)
         ENDIF
      ENDDO
      ENDDO

C---USE THIS INDEX FOR SPECIFYING VISIBLE BAND GROUND ALBEDO 
C   AS REFL(I):

      DO I=1,IMAX
         CR(I,NCLDS(I)+2)=REFL(I)
      ENDDO

C               ----------------
                     ENDIF
C               ----------------

C---NOW, THE NEAR IR BANDS; HERE THE GROUND CAN BE HANDLED AS PART
C   OF THE CLOUD LOOP

C               ----------------
                IF (N > 1) THEN
C               ----------------

      DO I=1,IMAX 
          CR(I,2)=CIRRF(I,2)*CAMT(I,2)
          CT(I,2)=1.0-CAMT(I,2)*(CIRRF(I,2)+CIRAB(I,2))
      ENDDO

      DO KK=1,KCLDS
      DO I=1,IMAX
C       print *, 'I,KK,BCLDS,CIRRF,CAMT,CIRAB=',I,KK,
C    &           BCLDS(I,KK),CIRRF(I,KK+2),CAMT(I,KK+2),CIRAB(I,KK+2)
         IF (BCLDS(I,KK)) THEN
            CR(I,KK+2)=CIRRF(I,KK+2)*CAMT(I,KK+2)
            CT(I,KK+2)=1.0-CAMT(I,KK+2)*
     &                    (CIRRF(I,KK+2)+CIRAB(I,KK+2))
         ENDIF
      ENDDO
      ENDDO

C               ----------------
                     ENDIF
C               ----------------

      DO K=1,LP1
      DO I=1,IMAX
         ALFAU(I,K)=0.00
      ENDDO
      ENDDO


C***FOR EXECUTION OF THE CLOUD LOOP, IT IS NECESSARY TO SEPARATE OUT
C   TRANSMISSION FCTNS AT THE TOP AND BOTTOM OF THE CLOUDS, FOR
C   EACH BAND N. THE REQUIRED QUANTITIES ARE:
C      TTD(I,KTOPSW(I,K),N)  K RUNS FROM 2 TO NCLDS(I)+2: 
C      TTD(I,KBTMSW(I,K),N)  K RUNS FROM 2 TO NCLDS(I)+1: 
C      TTU(I,KTOPSW(I,K),N)  K RUNS FROM 2 TO NCLDS(I)+2:
C      AND INVERSES OF THE ABOVE. THE ABOVE QUANTITIES ARE STORED 
C      IN TDCL1,TDCL2,TUCL1,TDCL1I,TDCL2I,TUCLI,RESPECTIVELY, AS
C      THEY HAVE MULTIPLE USE IN THE PGM.

C---COMPUTE GATHERS
      DO KK=1,KCLDS+1
      DO I=1,IMAX
         TDCL1(I,KK)=TTD(I,KTOPSW(I,KK+1),N)
         TUCL1(I,KK)=TTU(I,KTOPSW(I,KK+1),N)
      ENDDO
      ENDDO

      DO KK=1,KCLDS
      DO I=1,IMAX
         TDCL2(I,KK)=TTD(I,KBTMSW(I,KK+1),N)
      ENDDO
      ENDDO

C---COMPUTE INVERSES
      DO KK=1,KCLDS
      DO I=1,IMAX
         TDCL2I(I,KK)=1.0/TDCL2(I,KK)
      ENDDO
      ENDDO

      DO KK=1,KCLDS+1
      DO I=1,IMAX
        TDCL1I(I,KK)=1.0/TDCL1(I,KK)
        TUCL1I(I,KK)=1.0/TUCL1(I,KK)
      ENDDO
      ENDDO


C   TCLU(LL) IS TRANSMISSION FUNCTION FROM TOP OF NEXT LOWER
C   CLOUD TO TOP OF UPPER CLOUD. TCLD(LL) IS T.F. FROM TOP
C   OF NEXT LOWER CLOUD TO BOTTOM OF UPPER CLOUD. LL=NC1 IS
C   THE LOWEST BOTTOM CLOUD (THE GROUND) ; LL=1 IS THE
C   HIGHEST UPPER CLOUD. 

         TCLU = 0.0
         TCLD = 0.0
      DO KK=1,KCLDS
      DO I=1,IMAX
         TCLU(I,KK+1)=TDCL1(I,KK+1)*TDCL1I(I,KK)*CT(I,KK+1)
         TCLD(I,KK+1)=TDCL1(I,KK+1)*TDCL2I(I,KK) 
      ENDDO
      ENDDO

C***WE DEFINE TCLD (I,1) AS TTD(I,KTOPSW(I,2),N)
      DO I=1,IMAX
         TCLD(I,1)=TDCL1(I,1)
      ENDDO

      DO I=1,IMAX
         DFNCLU(I,1)=TCLD(I,1)
      ENDDO


C   THE FOLLOWING CALCULATION IS FOR ALFAT: THE RATIO BETWEEN
C   THE DOWNWARD FLUX AT THE TOP OF THE HIGHEST CLOUD (IF
C   PRESENT) OR THE GROUND TO THE UPWARD FLUX AT THE SAME
C   LEVEL, TAKING INTO ACCOUNT MULTIPLE REFLECTIONS FROM
C   CLOUDS, IF PRESENT

C  --- Reshape 2-D arrays to 1-D ---
      DO K=1,LP1
      DO I=1,IMAX
        I1=(K-1)*IMAX+I
        CR1D(I1)=CR(I,K)
        ALFU1D(I1)=ALFAU(I,K)
        TCLD1D(I1)=TCLD(I,K)
        TCLU1D(I1)=TCLU(I,K)
      ENDDO
      ENDDO

C     print *, 'KCLDS=',KCLDS

      DO KK=1,KCLDS
C     -------------

      DO I=1,IMAX
         TEMP1(I)=CR1D(INDX4(I,KK)+2*IMAX)
         TEMP2(I)=ALFU1D(INDX4(I,KK)+IMAX)
         TEMP3(I)=TCLU1D(INDX4(I,KK)+IMAX)
         TEMP4(I)=CR1D(INDX4(I,KK)+IMAX)
         TEMP5(I)=TCLD1D(INDX4(I,KK)+IMAX)
      ENDDO
C     print *, 'TEMP1=',TEMP1

      DO I=1,IMAX
         TEMP6(I)=(TEMP1(I)+TEMP2(I))*TEMP3(I)*TEMP3(I) /
     &            (1.0-(TEMP1(I)+TEMP2(I))*TEMP4(I)*TEMP5(I)*TEMP5(I))
      ENDDO
C     print *, 'TEMP6=',TEMP6

      DO I=1,IMAX
         ALFU1D(INDX4(I,KK))=TEMP6(I)
      ENDDO

C  --- Reshape 1-D array into 2-D array -----
      DO K=1,LP1
      DO I=1,IMAX
         I1=(K-1)*IMAX+I
         ALFAU(I,K)=ALFU1D(I1)
      ENDDO
      ENDDO

C     print *, 'ALFAU=',ALFAU
C     -------------
      ENDDO

C***DEFINE ALFA FROM ALFAU(I,1) AND CR(I,2):
C***ALFA IS THE SYSTEM REFLECTION COEFFICIENT ABOVE THE TOP CLOUD
C    (OR GROUND, IF NO CLOUD AT GRID PT I )

      DO I=1,IMAX
         ALFA(I)=ALFAU(I,1)+CR(I,2)
      ENDDO

C     UPWARD FLUX ABOVE TOPMOST CLOUD
      DO I=1,IMAX
         UFNCLU(I,1)=TCLD(I,1)*ALFA(I)
      ENDDO

      DO I=1,IMAX
         TEMP2(I)=TUCL1I(I,1)*UFNCLU(I,1)
      ENDDO
C     print *, 'TEMP2=',TEMP2

      DO K=1,KJTOP      
      DO I=1,IMAX
         IF (BJTOP(I,K)) THEN
            UFN(I,K,N)=TEMP2(I)*TTU(I,K,N)
         ENDIF
      ENDDO
      ENDDO
C     print *, 'UFN=',UFN

C   CALCULATE UFN AT CLOUD TOPS AND DFN AT CLOUD BOTTOMS

      DO I=1,IMAX
         VV(I)=1.0
      ENDDO

      DO KK=1,KCLDS 
C     -------------
C     print *, 'KK,KCLDS=',KK,KCLDS
C     print *, 'TCLU=',TCLU

      DO I=1,IMAX
         IF (BCLDS(I,KK)) THEN
            UFNCLU(I,KK+1)=ALFAU(I,KK)*VV(I)*TCLD(I,KK)/TCLU(I,KK+1)
            DFNCLD(I,KK)=VV(I)*TCLD(I,KK)*
     &                   TCLU(I,KK+1)*TDCL2(I,KK)*TDCL1I(I,KK+1) +
     &                   UFNCLU(I,KK+1)*TCLD(I,KK+1)*CR(I,KK+1)
         ELSE
            UFNCLU(I,KK+1)=UFN(I,LP1,N)
            DFNCLD(I,KK)=DFN(I,LP1,N)
         ENDIF
      ENDDO
C     print *, 'UFNCLU=',UFNCLU
C     print *, 'DFNCLD=',DFNCLD

      DO I=1,IMAX
         VV(I)=DFNCLD(I,KK)
      ENDDO

C     print *, 'VV=',VV
C     print *, 'KTOPSW(KK+2)=',KTOPSW(:,KK+2)
C     print *, 'KBTMSW(KK+1)=',KBTMSW(:,KK+1)

      DO I=1,IMAX
         UFN(I,KTOPSW(I,KK+2),N)=UFNCLU(I,KK+1)
         DFN(I,KBTMSW(I,KK+1),N)=DFNCLD(I,KK)
      ENDDO

C     -------------
      ENDDO


C     NOW OBTAIN DFN AND UFN FOR LEVELS BETWEEN THE CLOUDS
      DO 2401 KK=1,KCLDS
C---SKIP IF THERE ARE NO SPACES BETWEEN CLOUD LAYERS KK AND KK+1,
C    FOR ANY GRID PT:
      IF (IBETCL(KK).EQ.0) GO TO 2401
         DO K=KCLDS2,LP1
         DO I=1,IMAX
            IF (BCLDS(I,KK) .AND. 
     1         (K.GE.KBTMSW(I,KK+1) .AND. K.LE.KTOPSW(I,KK+2))) THEN
                  UFN(I,K,N)=UFNCLU(I,KK+1)*TTU(I,K,N)*TUCL1I(I,KK+1)
                  DFN(I,K,N)=DFNCLD(I,KK)*TTD(I,K,N)*TDCL2I(I,KK)
            ENDIF
         ENDDO
         ENDDO
2401  CONTINUE


C     NOW OBTAIN DOWNWARD AND UPWARD FLUXES FOR LEVELS,IF ANY,
C     BETWEEN THE TOPS AND BOTTOMS OF CLOUDS. THE ASSUMPTION OF
C     CONSTANT HEATING RATE IN THESE REGIONS IS USED.

C***OBTAIN FLUXES AT TOP AND BOTTOM OF CLOUDS
      DO 2501 KK=1,KCLDS
C---SKIP IF THERE ARE NO "THICK CLOUDS" AT ALL IN CLOUD LEVEL KK
      IF (IINCL(KK).EQ.0) GO TO 2501
C
C***OBTAIN DOWNWARD FLUXES AT CLOUD TOPS AND UPWARD FLUXES AT
C   CLOUD BOTTOMS

      IF (KK.GT.1) THEN
         DO I=1,IMAX
            DFNCLU(I,KK)=DFN(I,KTOPSW(I,KK+1),N)
         ENDDO
      ENDIF

      DO I=1,IMAX
         UFNCLD(I,KK)=UFN(I,KBTMSW(I,KK+1),N)
      ENDDO

      DO I=1,IMAX
         TEMP1(I)=(UFNCLU(I,KK)-UFNCLD(I,KK))*DPCLD(I,KK)
         TEMP2(I)=(DFNCLU(I,KK)-DFNCLD(I,KK))*DPCLD(I,KK)
      ENDDO


      DO K=KCLDS2,LP1
      DO I=1,IMAX
         IF (BCLDS(I,KK) .AND.
     &       (K.GT.KTOPSW(I,KK+1) .AND. K.LT.KBTMSW(I,KK+1))) THEN
                UFN(I,K,N)=UFNCLU(I,KK)+TEMP1(I)*(PP(I,K)-PPTOP(I,KK))
                DFN(I,K,N)=DFNCLU(I,KK)+TEMP2(I)*(PP(I,K)-PPTOP(I,KK))
         ENDIF
      ENDDO
      ENDDO

2501  CONTINUE

C-----------------------------------------------------------------------

2301                           CONTINUE

C-----------------------------------------------------------------------

                                ENDIF

C***********************************************************************

C   CALCULATE ENTERING FLUX AT THE TOP
C   LOOP 860 SCALES THE DFN'S AND UFN'S TO THE CORRECT DFN(I,1,N)

      DO N=1,NB
      DO I=1,IMAX
         DFNTOP(I,N)=SSOLAR(I)*6.97667E5*PWTS(N)
COLD     DFNTOP(I,N)=SSOLAR*6.97667E5*COSZRO(I)*TAUDAR(I)*PWTS(N)
      ENDDO
      ENDDO

      DO N=1,NB
      DO K=1,LP1
      DO I=1,IMAX
         DFN(I,K,N)=DFN(I,K,N)*DFNTOP(I,N)
         UFN(I,K,N)=UFN(I,K,N)*DFNTOP(I,N)
      ENDDO
      ENDDO
      ENDDO

C   SUM OVER BANDS

      DO K=1,LP1
      DO I=1,IMAX
         DFSW(I,K)=DFN(I,K,1)
         UFSW(I,K)=UFN(I,K,1)
      ENDDO
      ENDDO

      DO N=2,NB
      DO K=1,LP1
      DO I=1,IMAX
         DFSW(I,K)=DFSW(I,K)+DFN(I,K,N)
         UFSW(I,K)=UFSW(I,K)+UFN(I,K,N)
      ENDDO
      ENDDO
      ENDDO

      DO K=1,LP1
      DO I=1,IMAX
         FSW(I,K)=UFSW(I,K)-DFSW(I,K)
      ENDDO
      ENDDO

      DO K=1,LMAX
      DO I=1,IMAX
         HSW(I,K)=RADCON*(FSW(I,K+1)-FSW(I,K))/DP(I,K)
      ENDDO
      ENDDO

C-----------------------------------------------------------------------

      END SUBROUTINE SWRAD_ORIG

C#######################################################################

      END MODULE SHORTWAVE_MOD
