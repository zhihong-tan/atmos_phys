
module moist_conv_mod

!-----------------------------------------------------------------------

use  sat_vapor_pres_mod, ONLY: EsComp, DEsComp
use       utilities_mod, ONLY:  error_mesg, file_exist, open_file,  &
                                check_nml_error, close_file,        &
                                FATAL, WARNING, NOTE, get_my_pe
use       constants_mod, ONLY: HLv, HLs, cp, grav, rdgas, rvgas

implicit none
private

!------- interfaces in this module ------------

public :: moist_conv, moist_conv_Init

!-----------------------------------------------------------------------
!---- namelist ----

 real :: HC   = 1.00
 real :: beta = 0.0
 real :: TOLmin=.02, TOLmax=.10
 integer :: ITSMOD=30

!----- note beta is the fraction of convective condensation that is
!----- detrained into a stratiform cloud

 namelist /moist_conv_nml/  HC, beta, TOLmin, TOLmax, ITSMOD

!-----------------------------------------------------------------------
!---- VERSION NUMBER -----

 character(len=128) :: version = '$Id: moist_conv.F90,v 1.2 2000/07/28 20:16:40 fms Exp $'
 character(len=128) :: tag = '$Name: damascus $'

!---------- initialize constants used by this module -------------------

 real, parameter :: d622 = rdgas/rvgas
 real, parameter :: d378 = 1.0-d622
 real, parameter :: grav_inv = 1.0/grav
 real, parameter :: rocp = rdgas/cp

 logical :: do_init=.true.

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CONTAINS

!#######################################################################

 subroutine moist_conv ( Tin, Qin, Pfull, Phalf, coldT,    &
                         Tdel, Qdel, Rain, Snow, Lbot, cf, &
                         Conv, cfdel, qldel, qidel)

!-----------------------------------------------------------------------
!
!                       MOIST CONVECTIVE ADJUSTMENT
!
!-----------------------------------------------------------------------
!
!   INPUT:   Tin     temperature at full model levels
!            Qin     specific humidity of water vapor at full
!                      model levels
!            Pfull   pressure at full model levels
!            Phalf   pressure at half model levels
!            coldT   Should MCA produce snow in this column?
!
!   OUTPUT:  Tdel    temperature adjustment at full model levels (deg k)
!            Qdel    specific humidity adjustment of water vapor at
!                       full model levels
!            Rain    liquid precipitiation (in Kg m-2)
!            Snow    ice phase precipitation (kg m-2)
!  OPTIONAL
!
!   INPUT:   Lbot    integer index of the lowest model level,
!                      Lbot is always <= size(Tin,3)
!              cf    stratiform cloud fraction (used only when
!                    operating with stratiform cloud scheme) (fraction)
!
!  OUTPUT:   Conv    logical flag; TRUE then moist convective
!                       adjustment was performed at that model level.
!            cfdel   change in stratiform cloud fraction (fraction)
!            qldel   change in liquid water condensate due to
!                    convective detrainment (kg condensate /kg air)
!            qidel   change in ice condensate due to
!                    convective detrainment (kg condensate /kg air)
!
!-----------------------------------------------------------------------
!----------------------PUBLIC INTERFACE ARRAYS--------------------------
    real, intent(IN) , dimension(:,:,:) :: Tin, Qin, Pfull, Phalf
 logical, intent(IN) , dimension(:,:)   :: coldT
    real, intent(OUT), dimension(:,:,:) :: Tdel, Qdel
    real, intent(OUT), dimension(:,:)   :: Rain, Snow
 integer, intent(IN) , optional, dimension(:,:)   :: Lbot
    real, intent(IN) , optional, dimension(:,:,:) :: cf
 logical, intent(OUT), optional, dimension(:,:,:) :: Conv
    real, intent(OUT), optional, dimension(:,:,:) :: cfdel, qldel, qidel
      
!-----------------------------------------------------------------------
!----------------------PRIVATE (LOCAL) ARRAYS---------------------------
! logical, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) :: DO_ADJUST
!    real, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) ::  &
!-----------------------------------------------------------------------
!----------------------PRIVATE (LOCAL) ARRAYS---------------------------
integer, dimension(size(Tin,1),size(Tin,2)) :: ISMVD,ISMVF
integer, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) :: IVF

real, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) ::   &
  Qdif,Temp,Qmix,Esat,Qsat,Test1,Test2

real, dimension(size(Tin,1),size(Tin,2),size(Tin,3)-1) ::   &
  Thalf,DelPoP,Esm,Esd,ALRM

real, dimension(size(Tin,3)) :: C,Ta,Qa
real, dimension(size(Tin,1),size(Tin,2)) :: HL

integer :: i,j,k,kk,KX,ITER,MXLEV,MXLEV1,kstart,KTOP,KBOT,KBOTM1
real    :: ALTOL,Sum0,Sum1,Sum2,EsDiff,EsVal,Thaf,Pdelta,DTinv
!-----------------------------------------------------------------------

      if (do_init) call moist_conv_init ( )

      KX=size(Tin,3)

!------ compute Proper HL
      WHERE (coldT)
            HL = HLs
      ELSEWHERE
            HL = HLv
      END WHERE

!------ convert spec hum to mixing ratio ------
      Temp(:,:,:)=Tin(:,:,:)
      Qmix(:,:,:)=Qin(:,:,:)

      do k=1,KX-1
         DelPoP(:,:,k)=(Pfull(:,:,k+1)-Pfull(:,:,k))/Phalf(:,:,k+1)
      enddo

!-------------SATURATION VAPOR PRESSURE FROM ETABL----------------------

      call EsComp (Temp,Esat)

      Esat(:,:,:)=Esat(:,:,:)*HC
      Qsat(:,:,:)=Pfull(:,:,:)-d378*Esat(:,:,:)
      Qsat(:,:,:)=Max(0.0,d622*Esat(:,:,:)/Qsat(:,:,:))
      Qdif(:,:,:)=Max(0.0,Qmix(:,:,:)-Qsat(:,:,:))

!-----------------------------------------------------------------------
!                  MOIST CONVECTIVE ADJUSTMENT
!-----------------------------------------------------------------------

!  *** Set initial tolerance ***

           ALTOL = TOLmin

      do k=1,KX-1
         Thalf(:,:,k)=0.50*(Temp(:,:,k)+Temp(:,:,k+1))
      enddo

      call  EsComp (Thalf,Esm)
      call DEsComp (Thalf,Esd)

      do k=1,KX-1
         ALRM(:,:,k)=rocp*DelPoP(:,:,k)*Thalf(:,:,k)  &
         *(Phalf(:,:,k+1)+d622*HL(:,:)*Esm(:,:,k)/Thalf(:,:,k)/rdgas)  &
         /(Phalf(:,:,k+1)+d622*HL(:,:)*Esd(:,:,k)/cp)
      enddo

      IVF  (:,:,KX)=0
      Test1(:,:,KX)=0.0
      Test2(:,:,KX)=0.0

      do k=1,KX-1
         Test1(:,:,k)=Temp(:,:,k+1)-Temp(:,:,k)
         Test2(:,:,k)=ALRM(:,:,k)+ALTOL-Test1(:,:,k)
      enddo

!!!!! Test1(:,:,:)=0.0-Qdif(:,:,:)
      Test1(:,:,:)=(0.0-Qdif(:,:,:))*Qsat(:,:,:)

!-------IVF=1 in unstable layers where both levels are saturated--------

      do k=1,KX-1
         where (Test1(:,:,k) < 0.0 .and. Test1(:,:,k+1) < 0.0 .and. &
                Test2(:,:,k) < 0.0)
                         IVF(:,:,k)=1
         elsewhere
                         IVF(:,:,k)=0
         endwhere
      enddo

!  ------ Set convection flag (for optional output only) --------

      if (Present(Conv)) then
         Conv(:,:,1)=(IVF(:,:,1) == 1)
         do k=1,KX-1
            Conv(:,:,k+1)=(IVF(:,:,k) == 1 .or. IVF(:,:,k+1) == 1)
         enddo
      endif

!  ----- Set counter for each column -----

         ISMVF(:,:)=0
      do k=1,KX-1
         ISMVF(:,:)=ISMVF(:,:)+IVF(:,:,k)
      enddo

!-----------------------------------------------------------------------
!---------------LOOP OVER EACH VERTICAL COLUMN--------------------------
                       do j=1,size(Tin,2)
           OUTER_LOOP: do i=1,size(Tin,1)
!-----------------------------------------------------------------------
      if (ISMVF(i,j) == 0)  CYCLE

      if (Present(Lbot)) then
         MXLEV=Lbot(i,j)
      else
         MXLEV=KX
      endif
      MXLEV1=MXLEV-1

!  *** Re-set initial tolerance ***
           ALTOL = TOLmin

!----------(return here after increasing tolerance)--------------------
1450  CONTINUE

!--------------Iterations at the same tolerance-------------------------
                     do 1740 ITER=1,ITSMOD
!-----------------------------------------------------------------------
      kstart=1
 1500 CONTINUE
!-------------TEST TO DETERMINE UNSTABLE LAYER BLOCKS-------------------
!-------Find top (KTOP) and bottom (KBOT) of unstable layers------------
      do k=kstart,MXLEV1
          if (IVF(i,j,k) == 1)  GO TO 1505
      enddo
      CYCLE OUTER_LOOP
1505  KTOP=k

      do k=KTOP,MXLEV1
          if (IVF(i,j,k+1) == 0) then
             KBOT=k+1
             GO TO 1510
          endif
      enddo
      KBOT=MXLEV
1510  CONTINUE
!-----------------------------------------------------------------------

      KBOTM1=KBOT - 1
      Sum1=0.0
      Sum2=0.0
!-----------------------------------------------------------------------
                      do 1630 k=KTOP,KBOT
!-----------------------------------------------------------------------
         C(k)=Pfull(i,j,k)-d378*Esat(i,j,k)
      if (C(k) <= 0.0) then
         C(k)=0.0
      else
!DIR$ INLINE
         call DEsComp (Temp(i,j,k),EsDiff)
!DIR$ NOINLINE
         C(k)=d622*Pfull(i,j,k)*HC*EsDiff/  &
                        ((Pfull(i,j,k)-d378*Esat(i,j,k))**2)
      endif

      Sum0=0.0
      if (k == KBOT) GO TO 1625
      kk=k
1620  if (kk > KBOTM1) GO TO 1625
      Sum0=Sum0+ALRM(i,j,kk)
      kk=kk+1
      GO TO 1620

1625  CONTINUE
      Pdelta=Phalf(i,j,k+1)-Phalf(i,j,k)
      Sum1=Sum1 + Pdelta*((cp+HL(i,j)*C(k))*(Temp(i,j,k)+Sum0)+  &
                           HL(i,j)*(Qmix(i,j,k)-Qsat(i,j,k)))
      Sum2=Sum2 + Pdelta*(cp+HL(i,j)*C(k))
!-----------------------------------------------------------------------
1630                   CONTINUE
!-----------------------------------------------------------------------
      Ta(KBOT)=Sum1/Sum2
      k=KTOP
1645  if (k > KBOTM1) GO TO 1641
      Sum0=0.0
      kk=k
1640  if (kk > KBOTM1) GO TO 1642
      Sum0=Sum0+ALRM(i,j,kk)
      kk=kk+1
      GO TO 1640
1642  Ta(k)=Ta(KBOT)-Sum0
      k=k+1
      GO TO 1645

!---------UPDATE T,R,ES,Esm,Esd & Qsat FOR THE ADJUSTED POINTS----------

1641  do k=KTOP,KBOT
        Qa(k)=Qsat(i,j,k)+C(k)*(Ta(k)-Temp(i,j,k))
        Temp(i,j,k)=Ta(k)
        Qmix(i,j,k)=Qa(k)
!DIR$ INLINE
        call EsComp (Temp(i,j,k),EsVal)
!DIR$ NOINLINE
        Esat(i,j,k)=HC*EsVal
        Qsat(i,j,k)=Pfull(i,j,k)-d378*Esat(i,j,k)
        Qsat(i,j,k)=Max(0.0,d622*Esat(i,j,k)/Qsat(i,j,k))
        Qdif(i,j,k)=Max(0.0,Qmix(i,j,k)-Qsat(i,j,k))
      enddo

      do k=KTOP,KBOTM1
        Thaf=0.50*(Temp(i,j,k)+Temp(i,j,k+1))
!DIR$ INLINE
        call  EsComp (Thaf,EsVal)
        call DEsComp (Thaf,EsDiff)
!DIR$ NOINLINE
        Esm (i,j,k)=HC*EsVal
        Esd (i,j,k)=HC*EsDiff
        ALRM(i,j,k)=rocp*DelPoP(i,j,k)*Thaf*  &
               (Phalf(i,j,k+1)+d622*HL(i,j)*Esm(i,j,k)/Thaf/rdgas)/  &
               (Phalf(i,j,k+1)+d622*HL(i,j)*Esd(i,j,k)/cp)
      enddo

!------------Is this the bottom of the current column ???---------------
      kstart=KBOT+1
      if (kstart <= MXLEV1) GO TO 1500
!-----------------------------------------------------------------------
                if (ITER == ITSMOD) GO TO 1740
!-----------------------------------------------------------------------

      do k=1,MXLEV1
        Thaf=0.50*(Temp(i,j,k)+Temp(i,j,k+1))
!DIR$ INLINE
        call  EsComp (Thaf,EsVal)
        call DEsComp (Thaf,EsDiff)
!DIR$ NOINLINE
        Esm (i,j,k)=HC*EsVal
        Esd (i,j,k)=HC*EsDiff
        ALRM(i,j,k)=rocp*DelPoP(i,j,k)*Thaf*  &
               (Phalf(i,j,k+1)+d622*HL(i,j)*Esm(i,j,k)/Thaf/rdgas)/  &
               (Phalf(i,j,k+1)+d622*HL(i,j)*Esd(i,j,k)/cp)
      enddo

      do k=1,MXLEV1
        IVF(i,j,k)=0
!!!!    if (Qdif(i,j,k) > 0.0 .and. Qdif(i,j,k+1) > 0.0 .and.  &
        if (Qdif(i,j,k)*Qsat(i,j,k) > 0.0 .and.     &
             Qdif(i,j,k+1)*Qsat(i,j,k+1) > 0.0 .and.  &
              (Temp(i,j,k+1)-Temp(i,j,k)) > (ALRM(i,j,k)+ALTOL)) then
                       IVF(i,j,k) = 1
        endif
      enddo

!   ------ reset optional convection flag ------

      if (Present(Conv)) then
         Conv(i,j,1)=(IVF(i,j,1) == 1)
         do k=1,MXLEV1
            Conv(i,j,k+1)=(IVF(i,j,k) == 1 .or. IVF(i,j,k+1) == 1)
         enddo
      endif

!   ------ Are all layers sufficiently stable ??? ------

              ISMVF(i,j)=0
           do k=1,MXLEV1
              ISMVF(i,j)=ISMVF(i,j)+IVF(i,j,k)
           enddo
              if (ISMVF(i,j) == 0) CYCLE OUTER_LOOP

!-----------------------------------------------------------------------
1740                       CONTINUE
!-----------------------------------------------------------------------

!---------Maximum iterations reached: Increase tolerance (ALTOL)--------
      ALTOL=2.0*ALTOL
!del  WRITE (*,9902) I,ALTOL
      call error_mesg ('moist_conv', 'Tolerence (ALTOL) doubled', NOTE)
      if (ALTOL <= TOLmax)  GO TO 1450

!     WRITE (*,9903)
!     WRITE (*,9904) (k,Temp(i,j,k),Qmix(i,j,k),Qsat(i,j,k),  &
!                       Qdif(i,j,k),ALRM(i,j,k),k=1,MXLEV1)
!     WRITE (*,9904) (k,Temp(i,j,k),Qmix(i,j,k),Qsat(i,j,k),  &
!                       Qdif(i,j,k)            ,k=MXLEV,MXLEV)

   call error_mesg ('moist_conv', 'maximum iterations reached', WARNING)
!-----------------------------------------------------------------------
                            enddo OUTER_LOOP
                            enddo
!-----------------------------------------------------------------------
!---------------------- END OF i,j LOOP --------------------------------
!-----------------------------------------------------------------------

!---- call Convective Detrainment subroutine -----

     if (Present(cf)) then

          !reset quantities
          cfdel = 0.
          qldel = 0.
          qidel = 0.
     
          CALL CONV_DETR(Qmix,Qin,Phalf,Temp,cf,coldT,cfdel,qldel,qidel)
          
     endif

!----- compute adjustments to temp and spec hum ----

      Tdel(:,:,:)=Temp(:,:,:)-Tin(:,:,:)
      Qdel(:,:,:)=Qmix(:,:,:)-Qin(:,:,:)

!----- integrate precip -----

      Rain(:,:)=0.0
      Snow(:,:)=0.0
   do k =1,KX

      WHERE(coldT(:,:)) 
      Snow(:,:)=Snow(:,:)+(Phalf(:,:,k)-Phalf(:,:,k+1))*  &
                               Qdel(:,:,k)*grav_inv
      ELSEWHERE
      Rain(:,:)=Rain(:,:)+(Phalf(:,:,k)-Phalf(:,:,k+1))*  &
                               Qdel(:,:,k)*grav_inv
      END WHERE

      !subtract off detrained condensate from surface precip
      if (present(cf)) then
      WHERE(coldT(:,:)) 
      Snow(:,:)=Snow(:,:)+(Phalf(:,:,k)-Phalf(:,:,k+1))*  &
                               qidel(:,:,k)*grav_inv
      ELSEWHERE
      Rain(:,:)=Rain(:,:)+(Phalf(:,:,k)-Phalf(:,:,k+1))*  &
                               qldel(:,:,k)*grav_inv
      END WHERE      
      end if

   enddo
      Rain(:,:)=Max(Rain(:,:),0.0)
      Snow(:,:)=Max(Snow(:,:),0.0)
!-----------------------------------------------------------------------
!-----------------   PRINT FORMATS   -----------------------------------

 9902 FORMAT(    ' *** ALTOL DOUBLED IN CONVAD AT I=',  &
                 I5,' ,ALTOL=', F10.4 )
 9903 FORMAT(/,' *** DIVERGENCE IN MOIST CONVECTIVE ADJUSTMENT ',/,  &
           4X,'K',14X,'T',14X,'R',13X,'Qsat',14X,'Qdif',12X,'ALRM',/)
 9904 FORMAT (I5,5E15.7)
!-----------------------------------------------------------------------

 end subroutine moist_conv

!#######################################################################
!#######################################################################

SUBROUTINE CONV_DETR(qvout,qvin,phalf,T,cf,coldT,cfdel,qldel,qidel)


IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This subroutine takes a fraction of the water condensed
!      by the convection scheme and detrains in the top level
!      undergoing convective adjustment.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	VARIABLES
!
!   ------
!   INPUT:
!   ------
!
!       qvout    water vapor specific humidity AFTER adjustment
!                (kg vapor/kg air)
!       qvin     water  vapor specific humidity BEFORE adjustment
!                (kg vapor/kg air)
!       phalf    pressure at model half levels (Pascals)
!       T        Temperature (Kelvin)
!       cf       cloud fraction (fraction)
!       coldT    is condensation of ice nature?
!
!   -------------
!   INPUT/OUTPUT:
!   -------------
!
!       cfdel    Change in cloud fraction due to detrainment (fraction)
!       qldel    Increase in liquid water due to detrainment
!                (kg condensate/kg air)
!       qidel    Increase in ice due to detrainment
!                (kg condensate/kg air)
!
!   -------------------
!	INTERNAL VARIABLES:
!   -------------------
!
!       precipsource  accumulated source of precipitation
!                     (kg condensate /meter/ (seconds*squared))
!       ktop          integer of top level undergoing convection
!       accum         logical variable indicating whether or not to
!                     add precip
!       fT            fraction of condensate that is liquid
!       i,j,k         looping variables
!       IDIM,JDIM,KDIM dimensions of input arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  User Interface variables
!  ------------------------

REAL,     INTENT (IN), DIMENSION(:,:,:)  :: qvout,qvin,T,cf
REAL,     INTENT (IN), DIMENSION(:,:,:)  :: phalf
LOGICAL,  INTENT (IN), DIMENSION(:,:)    :: coldT
REAL,     INTENT (INOUT),DIMENSION(:,:,:):: cfdel,qldel,qidel

!  Internal variables
!  ------------------

INTEGER                                  :: i,j,k,IDIM,JDIM,KDIM,ktop
LOGICAL                                  :: accum
REAL                                     :: precipsource

!
! Code
! ----

        ! reinitialize variables
        cfdel(:,:,:)   = 0.
        qidel(:,:,:)   = 0.
        qldel(:,:,:)   = 0.
        IDIM           = SIZE(qvout,1)
        JDIM           = SIZE(qvout,2)
        KDIM           = SIZE(qvout,3)

        !---loop over grid columns----!
        DO i = 1, IDIM
        DO j = 1, JDIM

             !reset variables
             precipsource     = 0.
             accum            = .FALSE.

             DO k = 1, KDIM

                 !begin new convective event
                 IF ((qvout(i,j,k) .ne. qvin(i,j,k)) .and. &
                     (.NOT. accum)) THEN
                     ktop  = k
                     accum = .TRUE.
                 END IF

                 !if convective event is over compute detrainment
                 IF ( (accum) .and. (qvout(i,j,k) .eq. qvin(i,j,k)) &
                      .and. (precipsource .gt. 0.)) THEN                      
                      if (coldT(i,j)) then
                      qidel(i,j,ktop) = beta * precipsource / &
                                    (phalf(i,j,ktop+1)-phalf(i,j,ktop))
                      else
                      qldel(i,j,ktop) = beta * precipsource / &
                                    (phalf(i,j,ktop+1)-phalf(i,j,ktop))
                      end if
                      cfdel(i,j,ktop) = MAX(0.,HC-cf(i,j,ktop))
                      accum        = .FALSE.
                      precipsource = 0.
                 END IF

                 !accumulate precip
                 IF (accum) THEN
                      precipsource = precipsource + &
                                   ( qvin(i,j,k)  -qvout(i,j,k))* &
                                   (phalf(i,j,k+1)-phalf(i,j,k))
                 END IF

             END DO    !---end k loop over vertical column

            !---clear any remaining precip
            IF ( (precipsource .gt. 0.) .and. (accum) ) THEN
                 if (coldT(i,j)) then
                 qidel(i,j,ktop) = beta * precipsource / &
                                    (phalf(i,j,ktop+1)-phalf(i,j,ktop))
                 else
                 qldel(i,j,ktop) = beta * precipsource / &
                                    (phalf(i,j,ktop+1)-phalf(i,j,ktop))
                 end if
            END IF

        END DO    !---end j loop
        END DO    !---end i loop

END SUBROUTINE CONV_DETR

!#######################################################################

 subroutine moist_conv_init ()

!-----------------------------------------------------------------------
      
 integer :: unit, io, ierr

!-----------------------------------------------------------------------

    if (file_exist('input.nml')) then
        unit = open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
            read  (unit, nml=moist_conv_nml, iostat=io, end=10)
            ierr = check_nml_error (io,'moist_conv_nml')
        enddo
 10     call close_file (unit)
    endif

!---------- output namelist --------------------------------------------

      unit = open_file ('logfile.out', action='append')
      if ( get_my_pe() == 0 ) then
           write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
           write (unit,nml=moist_conv_nml)
      endif
      call close_file (unit)


      do_init = .false.

!-----------------------------------------------------------------------

 end subroutine moist_conv_init

!#######################################################################

end module moist_conv_mod

