
                        MODULE RDPARM_MOD

!-----------------------------------------------------------------------
!
!   PARAMETER SETTINGS FOR THE LONGWAVE AND SHORTWAVE RADIATION CODE: 
!   ----------------------------------------------------------------- 
!
!          IMAX   =  NO. POINTS ALONG THE LAT. CIRCLE USED IN CALCS.
!          JMAX   =  NO. POINTS ALONG THE MERIDIONAL AXIS
!          LMAX   =  NO. VERTICAL LEVELS (ALSO LAYERS) IN MODEL 
!
!      *** NOTE: THE USER NORMALLY WILL MODIFY ONLY THE
!                IMAX AND LMAX VARIABLES 
!
!          NBLW   =  NO. FREQ. BANDS FOR APPROX COMPUTATIONS. SEE 
!                      BANDTA FOR DEFINITION
!          NBLX   =  NO. FREQ BANDS FOR APPROX CTS COMPUTATIONS 
!          NBLY   =  NO. FREQ. BANDS FOR EXACT CTS COMPUTATIONS. SEE
!                      BDCOMB FOR DEFINITION
!          INLTE  =  NO. LEVELS USED FOR NLTE CALCS.
!          NNLTE  =  INDEX NO. OF FREQ. BAND IN NLTE CALCS. 
!
!          NB,KO2 ARE SHORTWAVE PARAMETERS; OTHER QUANTITIES ARE
!                    DERIVED FROM THE ABOVE PARAMETERS. 
!
!-----------------------------------------------------------------------

!!!!  INTEGER, PUBLIC, SAVE :: IMAX,JMAX,LMAX
      INTEGER, PUBLIC, SAVE :: LMAX=0
      INTEGER, PUBLIC, SAVE :: LP1,LP2,LP3,LM1,LM2,LM3
      INTEGER, PUBLIC, SAVE :: LL,LLP1,LLP2,LLP3,LLM1,LLM2,LLM3
      INTEGER, PUBLIC, SAVE :: LP1M,LP1M1,LP1V,LP121,LL3P
      INTEGER, PUBLIC, SAVE :: LP1I,LLP1I,LL3PI

      INTEGER, PARAMETER :: NBLW=163,NBLX=47,NBLY=15,NBLM=NBLY-1
      INTEGER, PARAMETER :: NB=9,NB1=NB-1
      INTEGER, PARAMETER :: INLTE=3,INLTEP=INLTE+1
      INTEGER, PARAMETER :: NNLTE=56
      INTEGER, PARAMETER :: KO2=12,KO21=KO2+1,KO2M=KO2-1


      CONTAINS

!#######################################################################

      SUBROUTINE RDPARM_INIT (KDIM)

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: KDIM

         LMAX=KDIM

         LP1=LMAX+1; LP2=LMAX+2; LP3=LMAX+3 
         LM1=LMAX-1; LM2=LMAX-2; LM3=LMAX-3 
         LL=2*LMAX; LLP1=LL+1; LLP2=LL+2; LLP3=LL+3
         LLM1=LL-1; LLM2=LL-2; LLM3=LL-3 
         LP1M=LP1*LP1; LP1M1=LP1M-1 
         LP1V=LP1*(1+2*LMAX/2)
         LP121=LP1*NBLY
         LL3P=3*LMAX+2

!!!!  Not Used ?????
!!!!     LP1I=IMAX*LP1; LLP1I=IMAX*LLP1; LL3PI=IMAX*LL3P 

      END SUBROUTINE RDPARM_INIT

!#######################################################################

                        END MODULE RDPARM_MOD

