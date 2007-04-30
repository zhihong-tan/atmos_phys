
MODULE CONV_UTILITIES_MOD
  
  use Sat_Vapor_Pres_Mod, ONLY: ESCOMP, DESCOMP
  use      Constants_Mod, ONLY: tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas

!---------------------------------------------------------------------
  implicit none
  private

!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: conv_utilities.F90,v 14.0 2007/03/15 22:08:41 fms Exp $'
  character(len=128) :: tagname = '$Name: nalanda_2007_04 $'

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: sd_init, sd_copy, sd_end, ac_init, ac_clear, ac_end, qsat, qses, exn, &
       conden, findt, pack_sd, pack_sd_lsm, extend_sd, adi_cloud

  real, parameter :: p00   = 1.E5
  real, parameter :: epsilo= rdgas/rvgas      !ratio of h2o to dry air molecular weights 
  real, parameter :: zvir  = rvgas/rdgas - 1. !rh2o/rair - 1

  character(len=7) :: mod_name = 'conv_utilities'

  public sounding
  type sounding
     logical  :: coldT
     integer  :: kmax, kinv, ktoppbl
     real     :: psfc, pinv, zinv, thvinv, land, pblht, qint
     real, pointer :: t     (:)=>NULL(), qv   (:)=>NULL(), u     (:)=>NULL()
     real, pointer :: v     (:)=>NULL(), ql   (:)=>NULL(), qi    (:)=>NULL()
     real, pointer :: qa    (:)=>NULL(), thc  (:)=>NULL(), qct   (:)=>NULL()
     real, pointer :: thv   (:)=>NULL(), rh   (:)=>NULL(), p     (:)=>NULL()
     real, pointer :: z     (:)=>NULL(), dp   (:)=>NULL(), dz    (:)=>NULL()
     real, pointer :: rho   (:)=>NULL(), nu   (:)=>NULL(), leff  (:)=>NULL()
     real, pointer :: exner (:)=>NULL(), ps   (:)=>NULL(), exners(:)=>NULL()
     real, pointer :: zs    (:)=>NULL(), ssthc(:)=>NULL(), ssqct (:)=>NULL()
     real, pointer :: dudp  (:)=>NULL(), dvdp (:)=>NULL(), thvbot(:)=>NULL()
     real, pointer :: thvtop(:)=>NULL(), qn   (:)=>NULL(), qs    (:)=>NULL()
     real, pointer :: am1   (:)=>NULL(), am2  (:)=>NULL(), am3   (:)=>NULL()
     real, pointer :: hl    (:)=>NULL(), sshl (:)=>NULL()
  end type sounding

  public adicloud
  type adicloud
     real     :: usrc, vsrc, hlsrc, thcsrc, qctsrc
     integer  :: klcl, klfc, klnb
     real     :: plcl, zlcl, thvlcl, thv0lcl, rho0lcl
     real     :: plfc, plnb, cape, cin
     real, pointer :: t  (:)=>NULL(), qv  (:)=>NULL(), ql  (:)=>NULL()
     real, pointer :: qi (:)=>NULL(), thc (:)=>NULL(), qct (:)=>NULL()
     real, pointer :: thv(:)=>NULL(), nu  (:)=>NULL(), leff(:)=>NULL()
     real, pointer :: hl (:)=>NULL()
  end type adicloud

contains

!#####################################################################
!#####################################################################

  subroutine sd_init(kd, sd)
    implicit none
    integer, intent(in) :: kd
    type(sounding), intent(inout) :: sd
    
    allocate ( sd%t     (1:kd)); sd%t     =0.;
    allocate ( sd%qv    (1:kd)); sd%qv    =0.;
    allocate ( sd%u     (1:kd)); sd%u     =0.;
    allocate ( sd%v     (1:kd)); sd%v     =0.;
    allocate ( sd%qs    (1:kd)); sd%qs    =0.;
    allocate ( sd%ql    (1:kd)); sd%ql    =0.;
    allocate ( sd%qi    (1:kd)); sd%qi    =0.;
    allocate ( sd%qa    (1:kd)); sd%qa    =0.;
    allocate ( sd%qn    (1:kd)); sd%qn    =0.;
    allocate ( sd%thc   (1:kd)); sd%thc   =0.;
    allocate ( sd%qct   (1:kd)); sd%qct   =0.;
    allocate ( sd%thv   (1:kd)); sd%thv   =0.;
    allocate ( sd%rh    (1:kd)); sd%rh    =0.;
    allocate ( sd%p     (1:kd)); sd%p     =0.;
    allocate ( sd%z     (1:kd)); sd%z     =0.;
    allocate ( sd%dp    (1:kd)); sd%dp    =0.;
    allocate ( sd%dz    (1:kd)); sd%dz    =0.;
    allocate ( sd%rho   (1:kd)); sd%rho   =0.;
    allocate ( sd%nu    (1:kd)); sd%nu    =0.;
    allocate ( sd%leff  (1:kd)); sd%leff  =0.;
    allocate ( sd%exner (1:kd)); sd%exner =0.;
    allocate ( sd%ps    (0:kd)); sd%ps    =0.;
    allocate ( sd%zs    (0:kd)); sd%zs    =0.;
    allocate ( sd%exners(0:kd)); sd%exners=0.;
    allocate ( sd%ssthc (1:kd)); sd%ssthc =0.;
    allocate ( sd%ssqct (1:kd)); sd%ssqct =0.;
    allocate ( sd%dudp  (1:kd)); sd%dudp  =0.;
    allocate ( sd%dvdp  (1:kd)); sd%dvdp  =0.;
    allocate ( sd%thvbot(1:kd)); sd%thvbot=0.;
    allocate ( sd%thvtop(1:kd)); sd%thvtop=0.;
    allocate ( sd%am1   (1:kd)); sd%am1   =0.;
    allocate ( sd%am2   (1:kd)); sd%am2   =0.;
    allocate ( sd%am3   (1:kd)); sd%am3   =0.;
    allocate ( sd%hl    (1:kd)); sd%hl    =0.;
    allocate ( sd%sshl  (1:kd)); sd%sshl  =0.;
  end subroutine sd_init

!#####################################################################
!#####################################################################

  subroutine sd_copy(sd, sd1)
    implicit none
    type(sounding), intent(in)    :: sd
    type(sounding), intent(inout) :: sd1
    
    sd1% kmax = sd % kmax
    sd1% land = sd % land
    sd1% coldT= sd % coldT
    sd1%p     = sd%p;    sd1%z     =sd%z;
    sd1%ps    = sd%ps;   sd1%zs    =sd%zs;
    sd1%t     = sd%t;    sd1%qv    =sd%qv;
    sd1%u     = sd%u;    sd1%v     =sd%v;
    sd1%ql    = sd%ql;   sd1%qi    =sd%qi;
    sd1%qa    = sd%qa;   sd1%qn    =sd%qn;    
    sd1%am1   = sd%am1;  sd1%am2   =sd%am2;  sd1%am3   =sd%am3;
    sd1%hl    = sd%hl;   sd1%sshl  =sd%sshl;
  end subroutine sd_copy

!#####################################################################
!#####################################################################

  subroutine sd_end(sd)
    type(sounding), intent(inout) :: sd
    deallocate ( sd%t, sd%qv, sd%u, sd%v, sd%ql, sd%qi, sd%qa, sd%thc, sd%qct,&
         sd%thv, sd%rh, sd%p, sd%z, sd%dp, sd%dz, sd%rho, sd%nu, sd%leff,     &
         sd%exner, sd%ps, sd%exners, sd%zs, sd%ssthc, sd%ssqct, sd%dudp,      &
         sd%dvdp, sd%thvbot, sd%thvtop, sd%qn, sd%am1, sd%am2, sd%am3, sd%qs, &
         sd%hl, sd%sshl)
  end subroutine sd_end

!#####################################################################
!#####################################################################

  subroutine ac_init(kd, ac)
    implicit none
    integer, intent(in) :: kd
    type(adicloud), intent(inout) :: ac
    
    allocate ( ac%t     (1:kd)); ac%t    =0.;
    allocate ( ac%qv    (1:kd)); ac%qv   =0.;
    allocate ( ac%ql    (1:kd)); ac%ql   =0.;
    allocate ( ac%qi    (1:kd)); ac%qi   =0.;
    allocate ( ac%thc   (1:kd)); ac%thc  =0.;
    allocate ( ac%qct   (1:kd)); ac%qct  =0.;
    allocate ( ac%thv   (1:kd)); ac%thv  =0.;
    allocate ( ac%nu    (1:kd)); ac%nu   =0.;
    allocate ( ac%leff  (1:kd)); ac%leff =0.;
    allocate ( ac%hl    (1:kd)); ac%hl   =0.;
  end subroutine ac_init

!#####################################################################
!#####################################################################

  subroutine ac_clear(ac)
    implicit none
    type(adicloud), intent(inout) :: ac
    ac%t    =0.;    ac%qv   =0.;    ac%ql   =0.;
    ac%qi   =0.;    ac%thc  =0.;    ac%qct  =0.;
    ac%thv  =0.;    ac%nu   =0.;    ac%leff =0.; ac%hl   =0.;
  end subroutine ac_clear

!#####################################################################
!#####################################################################

  subroutine ac_end(ac)
    type(adicloud), intent(inout) :: ac
    deallocate ( ac%t, ac%qv, ac%ql, ac%qi, ac%thc, ac%qct, ac%thv, ac%nu, ac%leff, ac%hl )
  end subroutine ac_end

!#####################################################################
!#####################################################################

  subroutine pack_sd(land, coldT, pmid, pint, zmid, zint, u, v, t, qv, ql, qi, qa, &
                    qn, am1, am2, am3, sd)
    implicit none

    real,    intent(in)              :: land
    logical, intent(in)              :: coldT
    real, intent(in), dimension(:)   :: pmid, zmid !pressure&height@mid level
    real, intent(in), dimension(:)   :: pint, zint !pressure&height@ interface level
    real, intent(in), dimension(:)   :: u, v       !wind profile (m/s)
    real, intent(in), dimension(:)   :: t, qv      !temperature and specific humidity
    real, intent(in), dimension(:)   :: ql, qi, qa, qn !cloud tracers
    real, intent(in), dimension(:)   :: am1, am2, am3  !3 aerosal species
    type(sounding), intent(inout)    :: sd

    integer :: k, nk, kmax

    !Pack environmental sounding; layers are numbered from bottom up!=
    sd % kmax   = size(t)
    sd % land   = land
    sd % coldT  = coldT
    sd % ps(0)  = pint(sd%kmax+1);
    sd % zs(0)  = zint(sd%kmax+1);

    do k=1, sd%kmax
       nk=sd%kmax-k+1
       sd % p     (k) = pmid(nk)
       sd % z     (k) = zmid(nk);
       sd % ps    (k) = pint(nk); 
       sd % zs    (k) = zint(nk); 
       sd % t     (k) = t   (nk)
       !prevent negative values for qv,ql,qi,qa,qn
       sd % qv    (k) = max(qv(nk), 4.e-10) 
       sd % ql    (k) = max(ql(nk), 0.)
       sd % qi    (k) = max(qi(nk), 0.)
       sd % qa    (k) = max(qa(nk), 0.)
       sd % qn    (k) = max(qn(nk), 0.)
       sd % u     (k) = u(nk)
       sd % v     (k) = v(nk)
       !yim's aerosol
       !sd % am1  (k) = am1(nk)
       !sd % am2  (k) = am2(nk)
       !sd % am3  (k) = am3(nk)
    end do
  end subroutine pack_sd

!#####################################################################
!#####################################################################

  subroutine extend_sd(sd, pblht, doice)
    implicit none
    type(sounding), intent(inout) :: sd
    real, intent(in)              :: pblht
    logical, intent(in)           :: doice

    integer :: k, kl, ktoppbl, id_check
    real    :: sshl0a, sshl0b, ssthc0a, ssthc0b, ssqct0a, ssqct0b
    real    :: hl0bot, thc0bot, qct0bot, hl0top, thc0top, qct0top
    real    :: thj, qvj, qlj, qij, qse, thvj

    sd % exners(0) = exn(sd%ps(0));
    if (doice) then
       sd%nu(:)= max(min((268. - sd % t(:))/20.,1.0),0.0);
    else
       sd%nu(:)=0.
    end if
    sd % leff(:) = (1-sd%nu(:))*HLv + sd%nu(:)*HLs
    sd % qct (:) = sd%qv(:)+sd%ql(:)+sd%qi(:)
    sd % hl  (:) = cp_air*sd%t(:)+grav*sd%z(:)-sd%leff(:)*(sd%ql(:)+sd%qi(:))
    sd % qint = 0.
    do k=1, sd%kmax
       sd % dp    (k) = sd%ps(k-1)-sd%ps(k)
       sd % dz    (k) = sd%zs(k)  -sd%zs(k-1)
       sd % exner (k) = exn(sd%p (k))
       sd % exners(k) = exn(sd%ps(k))
       !sd % thc   (k) = (sd%t(k) - sd%leff(k)*(sd%ql(k)+sd%qi(k))/cp_air)/sd%exner(k) 
       sd % thc   (k) = sd%t(k) / sd%exner(k) 
       sd % qs    (k) = qsat(sd%t(k), sd%p(k))
       sd % rh    (k) = min(sd%qv(k)/sd%qs(k),1.)
       sd % thv   (k) = sd%t(k)/sd%exner(k) * (1.+zvir*sd%qv(k)-sd%ql(k)-sd%qi(k))
       sd % rho   (k) = sd % p(k)/(rdgas * sd % thv(k) * sd % exner(k))
       sd % qint      =sd % qint + sd%qct(k)*sd%dp(k)
    end do
    sd % qint = sd % qint / grav

   !Finite-Volume intepolation
    kl=sd%kmax-1
    sshl0b  = (sd%hl (2)-sd%hl (1))/(sd%p(2)-sd%p(1))
    ssthc0b = (sd%thc(2)-sd%thc(1))/(sd%p(2)-sd%p(1))
    ssqct0b = (sd%qct(2)-sd%qct(1))/(sd%p(2)-sd%p(1))

    do k=2,kl
       sshl0a  = (sd%hl (k)-sd%hl (k-1))/(sd%p(k)-sd%p(k-1))
       if(sshl0a.gt.0)then
          sd%sshl (k-1) = max(0.,min(sshl0a,sshl0b))
       else
          sd%sshl (k-1) = min(0.,max(sshl0a,sshl0b))
       endif
       sshl0b = sshl0a
       ssthc0a = (sd%thc(k)-sd%thc(k-1))/(sd%p(k)-sd%p(k-1))
       if(ssthc0a.gt.0)then
          sd%ssthc(k-1) = max(0.,min(ssthc0a,ssthc0b))
       else
          sd%ssthc(k-1) = min(0.,max(ssthc0a,ssthc0b))
       endif
       ssthc0b = ssthc0a
       ssqct0a = (sd%qct(k)-sd%qct(k-1))/(sd%p(k)-sd%p(k-1))
       if(ssqct0a.gt.0)then
          sd%ssqct(k-1) = max(0.,min(ssqct0a,ssqct0b))
       else
          sd%ssqct(k-1) = min(0.,max(ssqct0a,ssqct0b))
       endif
       ssqct0b = ssqct0a
    enddo
    do k = 2,kl-1 !wind shear
       sd%dudp(k) = (sd%u(k+1)-sd%u(k-1))/(sd%p(k+1)-sd%p(k-1))
       sd%dvdp(k) = (sd%v(k+1)-sd%v(k-1))/(sd%p(k+1)-sd%p(k-1))
    end do
    sd%ssthc(kl)=sd%ssthc(kl-1)
    sd%ssqct(kl)=sd%ssqct(kl-1)

    do k = 1,kl !kl cannot be pver since ps0(pver)=0
       hl0bot  = sd%hl (k)+sd%sshl (k)*(sd%ps(k-1)-sd%p(k))
       thc0bot = sd%thc(k)+sd%ssthc(k)*(sd%ps(k-1)-sd%p(k))
       qct0bot = sd%qct(k)+sd%ssqct(k)*(sd%ps(k-1)-sd%p(k))
       !call conden(sd%ps(k-1),thc0bot,qct0bot,thj,qvj,qlj,qij,qse,sd%thvbot(k))
       call findt(sd%zs(k-1),sd%ps(k-1),hl0bot,qct0bot,thj,qvj,qlj,qij,qse,sd%thvbot(k),doice)
       
       hl0top  = sd%hl (k)+sd%sshl (k)*(sd%ps(k)-sd%p(k))
       thc0top = sd%thc(k)+sd%ssthc(k)*(sd%ps(k)-sd%p(k))
       qct0top = sd%qct(k)+sd%ssqct(k)*(sd%ps(k)-sd%p(k))
       !call conden(sd%ps(k),thc0top,qct0top,thj,qvj,qlj,qij,qse,sd%thvtop(k))
       call findt(sd%zs(k),sd%ps(k),hl0top,qct0top,thj,qvj,qlj,qij,qse,sd%thvtop(k),doice)
    enddo

    ktoppbl=1;
    do k = 2, kl-1
       if ((pblht+1.-sd%zs(k))*(pblht+1.-sd%zs(k+1)).lt.0.) then 
          ktoppbl=k; exit;
       endif
    end do
    !given a layer index k (here k=kinv); !its bottom interface 
    !level pressure is ps0(k-1) [here ps0(kinv-1)] and its bottom
    !interface level virt. pot. temperature is thv(k) [thv0bot(kinv)]
    sd % ktoppbl = ktoppbl
    sd % kinv    = ktoppbl+1 
    sd % pinv    = sd % ps    (sd % kinv-1) 
    sd % zinv    = sd % zs    (sd % kinv-1) 
    sd % thvinv  = sd % thvbot(sd % kinv)
    sd % pblht   = pblht
  end subroutine extend_sd

!#####################################################################
!#####################################################################

  subroutine adi_cloud(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, dofast, doice, ac)
  
    implicit none
    real, intent(inout) :: zsrc, psrc, hlsrc, thcsrc, qctsrc
    type(sounding), intent(in)    :: sd 
    logical,        intent(in)    :: dofast, doice
    type(adicloud), intent(inout) :: ac
    
    integer :: k, kl, klcl, klfc, klnb, id_check
    real    :: esrc,qs,tdsrc,temsrc,zlcl,tlcl
    real    :: hl0lcl, thc0lcl, qct0lcl, thv0lcl, rho0lcl
    real    :: cin, cinlcl, plfc, thvubot, thvutop
    real    :: thj, qvj, qlj, qij, qse, thvj
    real    :: cape, plnb, chi, tmp, rhtmp

    call ac_clear(ac);
    ac%klcl=0; ac%klfc=0; ac%klnb=0; 
    ac%plcl=0.; ac%zlcl=0.; ac%thvlcl=0.; ac%thv0lcl=0; ac%rho0lcl=0.;
    ac%plfc=0.; ac%plnb=0.; ac%cape=0.; ac%cin=0.;

!!$!below is pzhu's version now commented out
!!$    esrc=qctsrc*psrc/100./(qctsrc+epsilo)             ! water vapor pressure
!!$    tdsrc=tfreeze/(1-tfreeze*rvgas*log(esrc/6.11)/HLv)! dew-point of source air
!!$    temsrc=thcsrc*exn(psrc)                           ! temperature of source air
!!$    zlcl=123.5*(temsrc-tdsrc)+zsrc                    ! from sea-level
!!$    tlcl=temsrc-0.0098*(zlcl-zsrc)
!!$    ac % zlcl =zlcl
!!$    ac % plcl =psrc*(tlcl/temsrc)**(1./Kappa)


    !calculate lifted condensation level of air at parcel origin level
    !(within 0.2% of formula of Bolton, mon. wea. rev.,1980)
    !call conden(psrc, thcsrc, qctsrc, thj, qvj, qlj, qij, qse, thvj)
    call findt(zsrc, psrc, hlsrc, qctsrc, thj, qvj, qlj, qij, qse, thvj, doice)
    tmp=thj*exn(psrc)
    rhtmp=min(qctsrc/qse,1.)
    chi=tmp/(1669.0-122.0*rhtmp-tmp)
    ac%plcl=psrc*(rhtmp**chi); !Emanuel's calculation, results nearly identical to RAS

    klcl=0;  !klcl is the layer containing the LCL, i.e., ps0(klcl)<=plcl(i,j)
    do k=1,sd % kmax
       if(sd%ps(k).le.ac%plcl) then
          klcl=k; 
          ac%zlcl=sd%zs(k)-(ac%plcl-sd%ps(k))/sd%dp(k)*sd%dz(k);
          exit
       end if
    end do
    if (sd%ps(1).le.ac%plcl) then
       klcl=2; ac%plcl=sd%ps(1); ac%zlcl=sd%zs(1); 
    end if
    ac % klcl=klcl; 

    if (dofast.and.(ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.)) return;

    !call conden(ac % plcl, thcsrc, qctsrc, thj, qvj, qlj, qij, qse, ac%thvlcl)
    call findt(ac%zlcl, ac%plcl, hlsrc, qctsrc, thj, qvj, qlj, qij, qse, ac%thvlcl, doice)
    ac % hlsrc  = hlsrc
    ac % thcsrc = thcsrc
    ac % qctsrc = qctsrc

    hl0lcl  = sd%hl (klcl)+sd%sshl (klcl)*(ac%plcl-sd%p(klcl))
    thc0lcl = sd%thc(klcl)+sd%ssthc(klcl)*(ac%plcl-sd%p(klcl))
    qct0lcl = sd%qct(klcl)+sd%ssqct(klcl)*(ac%plcl-sd%p(klcl))
    !call conden(ac%plcl,thc0lcl,qct0lcl,thj,qvj,qlj,qij,qse,thv0lcl)
    call findt(ac%zlcl,ac%plcl,hl0lcl,qct0lcl,thj,qvj,qlj,qij,qse,thv0lcl, doice)
    rho0lcl = ac%plcl/(rdgas*thv0lcl*exn(ac%plcl))
    ac % thv0lcl= thv0lcl
    ac % rho0lcl= rho0lcl


    kl=sd % kmax-3
    do k=1,kl
       !call conden(sd%ps(k), thcsrc, qctsrc, thj, ac%qv(k), ac%ql(k), ac%qi(k), qs, ac%thv(k))
       call findt(sd%zs(k),sd%ps(k), hlsrc, qctsrc, thj, ac%qv(k), ac%ql(k), ac%qi(k), qs, ac%thv(k), doice)
       ac%t(k) = thj*exn(sd%ps(k))
!       qctsrc = qctsrc - (ac%ql(k) + ac%qi(k))
!       hlsrc  = hlsrc  + HLv*(ac%ql(k)+ac%qi(k))
    end do


    !Determine the convective inhibition (CIN)
    CIN = 0.
    cinlcl = 0.
    plfc   = 0.

    !define CIN based on LFC  
    do k = sd % kinv, kl-1
       if(k.eq.klcl-1) then !klcl-1 < layer < klcl
          thvubot=ac % thv (k); thvutop=ac % thvlcl
          call getcin(sd%ps(k),sd%thvtop(k),ac%plcl,thv0lcl,thvubot,thvutop,plfc,cin)
          cinlcl = cin
          thvubot=thvutop; thvutop=ac % thv (k+1)
          call getcin(ac%plcl,thv0lcl,sd%ps(k+1),sd%thvtop(k+1),thvubot,thvutop,plfc,cin)
          if(plfc.gt.0. .and. plfc.lt.ac%plcl) exit
       else
          thvubot=ac % thv (k); thvutop=ac % thv (k+1)
          call getcin(sd%ps(k),sd%thvtop(k),sd%ps(k+1),sd%thvtop(k+1),thvubot,thvutop,plfc,cin)
          if(plfc.gt.0. .and. plfc.lt.ac%plcl) exit
       endif
    enddo
 
    ac % cin =cin; !CIN has been estimated
    ac % plfc=plfc;

    if (dofast .and. (ac%plfc .lt. 500.) ) return; !miz

    !calculate cape=================
    if (ac%plfc.eq.0.0) then
       ac%cape=0.0;
       ac%plnb=0.0;
    else
       ac % klfc=0; !klfc is the layer containing the plfc, i.e., ps0(klfc)<=plfc(i,j)
       do k=1,kl 
          if(sd%ps(k).le.ac%plfc) then
             ac % klfc=max(k,2); exit
          end if
       end do
       plnb = 0.; cape=0.0;
       do k = ac % klfc-1, sd%kmax
          thvubot=ac % thv (k); thvutop=ac % thv (k+1)
          call getcape(sd%ps(k),sd%thvtop(k),sd%ps(k+1),sd%thvtop(k+1),thvubot,thvutop,plnb,cape)
          if(plnb.gt.0.) exit
       enddo
       ac%cape=cape
       ac%plnb=plnb
    end if

    return
  end subroutine adi_cloud

!#####################################################################
!#####################################################################

  function qsat(temp, p)
    real, intent(in)  :: temp, p
    real :: qsat, es, rs,t
    real :: tmin = 139.01, tmax = 399.99   

    t=max(temp,35.)
    if (t.le.tmin) then !accurate to within 0.14% in the range -80<T<0
       es=exp(23.33086-6111.72784/t+0.15215*log(t))*100.
       !print*, t,'Temperature <',tmin,'K in UW-ShCu'
    elseif (t.ge.tmax) then !accurate to within 0.3% in the range -35<T<35
       es=6.112*exp(17.67*(t-273.15)/(t-273.15+243.5))*100.
       !print*, t,'Temperature >',tmax,'K in UW-ShCu'
    else
       call escomp(t,es)
    end if
    rs=epsilo/(p/es-1.)
    qsat=rs/(1.+rs)
   
  end function qsat

!#####################################################################
!#####################################################################

  subroutine qses(temp, p, qs, es)
    real, intent(in)    :: temp, p
    real, intent(inout) :: qs, es
    real :: rs, t
    real :: tmin = 139.01, tmax = 399.99   

    t=max(temp,35.)
    if (t.le.tmin) then !accurate to within 0.14% in the range -80<T<0
       es=exp(23.33086-6111.72784/t+0.15215*log(t))*100.
       !print*, t,'Temperature <',tmin,'K in UW-ShCu'
    elseif (t.ge.tmax) then !accurate to within 0.3% in the range -35<T<35
       es=6.112*exp(17.67*(t-273.15)/(t-273.15+243.5))*100.
       !print*, t,'Temperature >',tmax,'K in UW-ShCu'
    else
       call escomp(t,es)
    end if
    rs=epsilo/(p/es-1.)
    qs=rs/(1.+rs)
    return
  end subroutine qses

!#####################################################################
!#####################################################################

  function exn(p)
    real :: exn
    real, intent(in)  :: p
    exn = (p/p00) ** Kappa
  end function exn

!#####################################################################
!#####################################################################

  subroutine conden(p,thc,qt,th,qv,ql,qi,qs,thv)
    implicit none
    real,     intent(in)  :: p, thc, qt
    real,     intent(out) :: th, qv, ql, qi, qs, thv
    real      tc, exn, leff, nu, qc, temps, tc1
    integer   iteration, id_check
    integer :: niteration = 5

    exn = (p/p00)**Kappa
    tc = thc * exn
    nu = max(min((268.-tc)/20.,1.0),0.0);
    leff = (1-nu)*HLv + nu*HLs
  
    temps = tc
    qs=qsat(temps,p)
    
    if(qs.gt.qt) then
       id_check=0
    else
       do iteration = 1,niteration
          temps = temps + ((tc-temps)*cp_air/leff + (qt -qs))/        &
               (cp_air/leff + epsilo*leff*qs/rdgas/temps/temps)
          qs = qsat(temps,p)
       enddo
       tc1=temps-leff/cp_air*(qt-qs)
       if(abs(tc1-tc).lt.1.0) then
          id_check=0
       else
          id_check=1; print*,'ID_CHECK=11111111111111111111111111111111'
       endif
    endif
    qc = max(qt-qs, 0.)
    qv = qt - qc
    ql = (1-nu)*qc
    qi = nu*qc !temps=tc+leff/cp_air*qc
    th = temps/exn
    thv=th*(1.+zvir*qv-ql-qi)
    return
  end subroutine conden

!#####################################################################
!#####################################################################
  subroutine findt(z,p,hl,qt,th,qv,ql,qi,qs,thv,doice)
    implicit none
    real,     intent(in)  :: z, p, hl, qt
    logical,  intent(in)  :: doice
    real,     intent(out) :: th, qv, ql, qi, qs, thv
    real      tc, leff, nu, qc, temps, hl1, es, dqsdt, dherror
    integer   iteration, id_check
    integer,  save :: niteration = 10
 
    tc = (hl-grav*z)/cp_air
    if (doice) then
       nu = max(min((268.-tc)/20.,1.0),0.0);
    else
       nu = 0.;
    end if
    leff = (1.-nu)*HLv + nu*HLs
    temps = tc
    call qses(temps,p,qs,es)
      id_check = 1
    if(qs.gt.qt) then
       id_check=0
    else
       do iteration = 1,niteration
          !dqsdt=epsilo*leff*qs/rdgas/temps/temps
          dqsdt=epsilo*leff*es/(rdgas*temps*temps)
          dqsdt=dqsdt*epsilo*p/(p-es*(1.-epsilo))**2.
          temps=temps+((tc-temps)*cp_air/leff+(qt-qs))/(cp_air/leff+dqsdt)
          call qses(temps,p,qs,es)
 
          hl1=cp_air*temps-leff*(qt-qs)+grav*z
          dherror=abs(hl1-hl)
          if(dherror.lt.1.0) then
             id_check=0; exit
          end if
       enddo
       if (id_check.ne.0 .and. dherror .gt. 1000)  then
          print*,'ID_CHECK=111111',z,p,hl,qt,temps,qs,hl1,leff,niteration
       end if
    endif
    qc = max(qt-qs, 0.)
    qv = qt - qc
    ql = (1.-nu)*qc
    qi = nu*qc
    th = temps/exn(p)
    thv=th*(1.+zvir*qv-ql-qi)
    return
  end subroutine findt

!#####################################################################
!#####################################################################

  subroutine getcin(pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin)
    
    implicit none
    real,    intent(in)    :: pbot,thv0bot,ptop,thv0top,thvubot,thvutop
    real,    intent(inout) :: plfc,cin
    real                   :: frc, rhom, delp

    delp=(pbot-ptop)
    rhom = pbot/(rdgas*thv0bot*exn(pbot))+ptop/(rdgas*thv0top*exn(ptop))

    if(thvubot.gt.thv0bot.and.thvutop.gt.thv0top)then
       !Both top and bottom positively buoyant
       plfc = pbot
    elseif(thvubot.le.thv0bot.and.thvutop.le.thv0top)then 
       !Both top and bottom negatively buoyant
       cin  = cin  - ((thvubot/thv0bot-1.)+(thvutop/thv0top-1.)) * delp / rhom
    elseif(thvutop.le.thv0top.and.thvubot.gt.thv0bot)then
       !Top negatively buoyant; Bottom positively buoyant
       frc  = (thvutop/thv0top-1.)/((thvutop/thv0top-1.)-(thvubot/thv0bot-1.))
       delp = (ptop+frc*delp) - ptop
       cin  = cin - (thvutop/thv0top-1.) * delp / rhom
    else                                                  
       !Top positively buoyant; Bottom negatively buoyant
       frc = (thvubot/thv0bot-1.)/((thvubot/thv0bot-1.)-(thvutop/thv0top-1.))
       plfc = pbot - frc * (pbot-ptop)
       delp = pbot - plfc
       cin  = cin - (thvubot/thv0bot-1.) * delp / rhom
    endif
    return
  end subroutine getcin

!#####################################################################
!#####################################################################

  subroutine getcape(pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plnb,cape)
   
    implicit none
    real,    intent(in)    :: pbot,thv0bot,ptop,thv0top,thvubot,thvutop
    real,    intent(inout) :: plnb,cape
    real                   :: frc, rhom, delp

    delp=(pbot-ptop)
    rhom = pbot/(rdgas*thv0bot*exn(pbot))+ptop/(rdgas*thv0top*exn(ptop))

    if(thvubot.gt.thv0bot.and.thvutop.gt.thv0top)then
       !Both top and bottom positively buoyant
       cape = cape + ((thvubot/thv0bot-1.)+(thvutop/thv0top-1.)) * delp / rhom
    elseif(thvubot.le.thv0bot.and.thvutop.le.thv0top)then 
       !Both top and bottom negatively buoyant
       plnb = pbot
    elseif(thvutop.le.thv0top.and.thvubot.gt.thv0bot)then
       !Top negatively buoyant; Bottom positively buoyant
       frc  = (thvubot/thv0bot-1.)/((thvubot/thv0bot-1.)-(thvutop/thv0top-1.))
       plnb = pbot - frc * (pbot-ptop)
       delp = pbot - plnb
       cape = cape + (thvubot/thv0bot-1.) * delp / rhom
    else                                                  
       !Top positively buoyant; Bottom negatively buoyant
       frc  = (thvutop/thv0top-1.)/((thvutop/thv0top-1.)-(thvubot/thv0bot-1.))
       delp = (ptop+frc*delp) - ptop
       cape = cape + (thvutop/thv0top-1.) * delp / rhom
    endif

    return
  end subroutine getcape

!#####################################################################
!#####################################################################


!###################################################################
!###################################################################

subroutine pack_sd_lsm(pf, ph, zf, zh, t, qv, sd)

  implicit none
  real, dimension(:), intent(in)    :: pf, ph, zf, zh, t, qv
  type(sounding),     intent(inout) :: sd

  integer :: k, nk, kmax

  kmax=size(t)
  sd % kmax   = kmax
  sd % land   = 0.
  sd % coldT  = .false.
  sd % ps(0)  = ph(kmax+1)
  sd % zs(0)  = zh(kmax+1)

  do k=1, kmax
     nk=kmax-k+1
     sd % p (k) = pf(nk)
     sd % z (k) = zf(nk)
     sd % ps(k) = ph(nk)
     sd % zs(k) = zh(nk)
     sd % t (k) = t (nk)
     sd % qv(k) = max(qv(nk), 4.e-10) 
     sd % ql(k) = 0.
     sd % qi(k) = 0.
     sd % qa(k) = 0.
     sd % qn(k) = 0.
     sd % u (k) = 0.
     sd % v (k) = 0.
  end do

end subroutine pack_sd_lsm

!###################################################################
!###################################################################


end MODULE CONV_UTILITIES_MOD
