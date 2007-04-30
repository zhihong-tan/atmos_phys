
MODULE CONV_PLUMES_MOD

  use Sat_Vapor_Pres_Mod, ONLY: ESCOMP, DESCOMP
  use      Constants_Mod, ONLY: tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas
  use  aer_ccn_act_mod,   only: aer_ccn_act, aer_ccn_act2
  use  conv_utilities_mod,only: qsat, exn, findt, conden, adicloud, sounding

!---------------------------------------------------------------------
  implicit none
  private

!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id: conv_plumes.F90,v 14.0 2007/03/15 22:08:38 fms Exp $'
  character(len=128) :: tagname = '$Name: nalanda_2007_04 $'

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: cp_init, cp_end, cp_clear, ct_init, ct_end, ct_clear, &
       cumulus_plume, cumulus_tend, cumulus_downdraft

  character(len=11) :: mod_name = 'conv_plumes'

  real, parameter :: zvir  = rvgas/rdgas - 1. !rh2o/rair - 1

  public cpnlist
  type cpnlist
     real :: rle, rpen, rmaxfrac, rbuoy, rdrag, frac_drs, bigc
     real :: auto_th0, auto_rate, tcrit, cldhgt_max
     logical :: do_ice, do_ppen, do_edplume, do_micro, do_forcedlifting, do_topevap
  end type cpnlist

  public cplume
  type cplume
     integer :: ltop, let, krel
     real    :: cush, cldhgt, prel, zrel
     real, pointer :: thcu  (:)=>NULL(), qctu  (:)=>NULL(), uu    (:)=>NULL()
     real, pointer :: vu    (:)=>NULL(), qlu   (:)=>NULL(), qiu   (:)=>NULL()
     real, pointer :: pptr  (:)=>NULL(), ppti  (:)=>NULL(), wu    (:)=>NULL()
     real, pointer :: umf   (:)=>NULL(), emf   (:)=>NULL(), thvu  (:)=>NULL()
     real, pointer :: rei   (:)=>NULL(), fer   (:)=>NULL(), fdr   (:)=>NULL()
     real, pointer :: dp    (:)=>NULL(), thc   (:)=>NULL(), qct   (:)=>NULL()
     real, pointer :: ql    (:)=>NULL(), qi    (:)=>NULL(), qa    (:)=>NULL()
     real, pointer :: u     (:)=>NULL(), v     (:)=>NULL(), p     (:)=>NULL()
     real, pointer :: ps    (:)=>NULL(), ufrc  (:)=>NULL(), thvtop(:)=>NULL()
     real, pointer :: thvbot(:)=>NULL(), fdrsat(:)=>NULL(), z     (:)=>NULL()
     real, pointer :: qn    (:)=>NULL(), qnu   (:)=>NULL(), zs    (:)=>NULL()
     real, pointer :: hlu   (:)=>NULL(), hl    (:)=>NULL(), clu   (:)=>NULL()
     real, pointer :: ciu   (:)=>NULL(), buo   (:)=>NULL(), t     (:)=>NULL()
     real, pointer :: crate (:)=>NULL(), prate (:)=>NULL()
  end type cplume

  public ctend
  type ctend
     integer :: botlev, toplev
     real    :: rain, snow, denth, uav, vav, conint, freint, dtint, dqint
     real, pointer :: uten  (:)=>NULL(), vten  (:)=>NULL(), tten  (:)=>NULL()
     real, pointer :: qvten (:)=>NULL(), qlten (:)=>NULL(), qiten (:)=>NULL()
     real, pointer :: qaten (:)=>NULL(), thcten(:)=>NULL(), qctten(:)=>NULL()
     real, pointer :: qvdiv (:)=>NULL(), qldiv (:)=>NULL(), qidiv (:)=>NULL()
     real, pointer :: thcflx(:)=>NULL(), qctflx(:)=>NULL()
     real, pointer :: umflx (:)=>NULL(), vmflx (:)=>NULL(), qvflx (:)=>NULL()
     real, pointer :: qlflx (:)=>NULL(), qiflx (:)=>NULL(), qaflx (:)=>NULL()
     real, pointer :: qnflx (:)=>NULL(), qnten (:)=>NULL(), pflx  (:)=>NULL()
     real, pointer :: hlflx (:)=>NULL(), hlten (:)=>NULL()
  end type ctend

contains

!#####################################################################
!#####################################################################

  subroutine cp_init(kd, cp)
    implicit none
    integer, intent(in) :: kd
    type(cplume), intent(inout) :: cp
    
    allocate ( cp%hlu   (0:kd)); cp%hlu   =0.;
    allocate ( cp%thcu  (0:kd)); cp%thcu  =0.;
    allocate ( cp%qctu  (0:kd)); cp%qctu  =0.;
    allocate ( cp%uu    (0:kd)); cp%uu    =0.;
    allocate ( cp%vu    (0:kd)); cp%vu    =0.;
    allocate ( cp%qlu   (0:kd)); cp%qlu   =0.;
    allocate ( cp%qiu   (0:kd)); cp%qiu   =0.;
    allocate ( cp%clu   (0:kd)); cp%clu   =0.;
    allocate ( cp%ciu   (0:kd)); cp%ciu   =0.;
    allocate ( cp%buo   (0:kd)); cp%buo   =0.;
    allocate ( cp%t     (0:kd)); cp%t     =0.;
    allocate ( cp%crate (0:kd)); cp%crate =0.;
    allocate ( cp%prate (0:kd)); cp%prate =0.;
    allocate ( cp%qnu   (0:kd)); cp%qnu   =0.;
    allocate ( cp%pptr  (1:kd)); cp%pptr  =0.;
    allocate ( cp%ppti  (1:kd)); cp%ppti  =0.;
    allocate ( cp%wu    (0:kd)); cp%wu    =0.;
    allocate ( cp%umf   (0:kd)); cp%umf   =0.;
    allocate ( cp%emf   (0:kd)); cp%emf   =0.;
    allocate ( cp%thvu  (1:kd)); cp%thvu  =0.;
    allocate ( cp%rei   (1:kd)); cp%rei   =0.;
    allocate ( cp%fer   (1:kd)); cp%fer   =0.;
    allocate ( cp%fdr   (1:kd)); cp%fdr   =0.;
    allocate ( cp%dp    (1:kd)); cp%dp    =0.;
    allocate ( cp%hl    (1:kd)); cp%hl    =0.;
    allocate ( cp%thc   (1:kd)); cp%thc   =0.;
    allocate ( cp%qct   (1:kd)); cp%qct   =0.;
    allocate ( cp%u     (1:kd)); cp%u     =0.;
    allocate ( cp%v     (1:kd)); cp%v     =0.;
    allocate ( cp%ql    (1:kd)); cp%ql    =0.;
    allocate ( cp%qi    (1:kd)); cp%qi    =0.;
    allocate ( cp%qa    (1:kd)); cp%qa    =0.;
    allocate ( cp%qn    (1:kd)); cp%qn    =0.;
    allocate ( cp%p     (1:kd)); cp%p     =0.;
    allocate ( cp%ps    (0:kd)); cp%ps    =0.;
    allocate ( cp%z     (1:kd)); cp%z     =0.;
    allocate ( cp%zs    (0:kd)); cp%zs    =0.;
    allocate ( cp%ufrc  (1:kd)); cp%ufrc  =0.;
    allocate ( cp%thvbot(1:kd)); cp%thvbot=0.;
    allocate ( cp%thvtop(1:kd)); cp%thvtop=0.;
    allocate ( cp%fdrsat(1:kd)); cp%fdrsat=0.;
  end subroutine cp_init

!#####################################################################
!#####################################################################

  subroutine cp_end(cp)
    implicit none
    type(cplume), intent(inout) :: cp
    deallocate ( cp%thcu, cp%qctu, cp%uu, cp%vu, cp%qlu, cp%qiu, cp%pptr,     &
         cp%ppti, cp%wu, cp%umf, cp%emf, cp%thvu, cp%rei, cp%fer, cp%fdr,     &
         cp%dp, cp%thc, cp%qct, cp%u, cp%v, cp%p, cp%ps, cp%ufrc, cp%thvbot,  &
         cp%thvtop,  cp%fdrsat, cp%qnu, cp%ql, cp%qi, cp%qa, cp%qn, cp%z,     &
         cp%zs, cp%hl, cp%hlu, cp%clu, cp%ciu, cp%buo, cp%t, cp%crate, cp%prate)
  end subroutine cp_end

!#####################################################################
!#####################################################################

  subroutine cp_clear(cp)
    implicit none
    type(cplume), intent(inout) :: cp
    cp%thcu  =0.;    cp%qctu  =0.;    cp%uu    =0.;    cp%vu    =0.;
    cp%qlu   =0.;    cp%qiu   =0.;    cp%qnu   =0.;    cp%pptr  =0.;
    cp%ppti  =0.;    cp%wu    =0.;    cp%umf   =0.;    cp%emf   =0.;
    cp%thvu  =0.;    cp%rei   =0.;    cp%fer   =0.;    cp%fdr   =0.;
    cp%dp    =0.;    cp%thc   =0.;    cp%qct   =0.;    cp%u     =0.;
    cp%v     =0.;    cp%ql    =0.;    cp%qi    =0.;    cp%qa    =0.; 
    cp%qn    =0.;    cp%p     =0.;    cp%ps    =0.;
    cp%ufrc  =0.;    cp%thvbot=0.;    cp%thvtop=0.;    cp%hlu   =0.;
    cp%fdrsat=0.;    cp%z     =0.;    cp%zs    =0.;    cp%hl    =0.;
    cp%clu   =0.;    cp%ciu   =0.;    cp%buo   =0.;    cp%t     =0.;
    cp%crate =0.;    cp%prate =0.;
  end subroutine cp_clear

!#####################################################################
!#####################################################################

  subroutine ct_init(kd, ct)
    implicit none
    integer, intent(in) :: kd
    type(ctend), intent(inout) :: ct
    allocate ( ct%uten  (1:kd)); ct%uten  =0.;
    allocate ( ct%vten  (1:kd)); ct%vten  =0.;
    allocate ( ct%tten  (1:kd)); ct%tten  =0.;
    allocate ( ct%qvten (1:kd)); ct%qvten =0.;
    allocate ( ct%qlten (1:kd)); ct%qlten =0.;
    allocate ( ct%qiten (1:kd)); ct%qiten =0.;
    allocate ( ct%qaten (1:kd)); ct%qaten =0.;
    allocate ( ct%qnten (1:kd)); ct%qnten =0.;
    allocate ( ct%hlten (1:kd)); ct%hlten =0.;
    allocate ( ct%thcten(1:kd)); ct%thcten=0.;
    allocate ( ct%qctten(1:kd)); ct%qctten=0.;
    allocate ( ct%qvdiv (1:kd)); ct%qvdiv =0.;
    allocate ( ct%qldiv (1:kd)); ct%qldiv =0.;
    allocate ( ct%qidiv (1:kd)); ct%qidiv =0.;
    allocate ( ct%hlflx (0:kd)); ct%hlflx =0.;
    allocate ( ct%thcflx(0:kd)); ct%thcflx=0.;
    allocate ( ct%qctflx(0:kd)); ct%qctflx=0.;
    allocate ( ct%qvflx (0:kd)); ct%qvflx =0.;
    allocate ( ct%qlflx (0:kd)); ct%qlflx =0.;
    allocate ( ct%qiflx (0:kd)); ct%qiflx =0.;
    allocate ( ct%qaflx (0:kd)); ct%qaflx =0.;
    allocate ( ct%qnflx (0:kd)); ct%qnflx =0.;
    allocate ( ct%umflx (0:kd)); ct%umflx =0.;
    allocate ( ct%vmflx (0:kd)); ct%vmflx =0.;
    allocate ( ct%pflx  (1:kd)); ct%pflx  =0.;
  end subroutine ct_init

!#####################################################################
!#####################################################################

  subroutine ct_end(ct)
    implicit none
    type(ctend), intent(inout) :: ct
    deallocate ( ct%uten, ct%vten, ct%tten, ct%qvten, ct%qlten, ct%qiten,     &
         ct%qaten, ct%qnten, ct%hlten, ct%thcten, ct%qctten, ct%qvdiv, ct%qldiv, ct%qidiv, &
         ct%hlflx, ct%thcflx, ct%qctflx, ct%qvflx, ct%qlflx, ct%qiflx, ct%qaflx, ct%qnflx, &
         ct%umflx, ct%vmflx, ct%pflx)
  end subroutine ct_end

!#####################################################################
!#####################################################################

  subroutine ct_clear(ct)
    implicit none
    type(ctend), intent(inout) :: ct
    
    ct%uten  =0.;    ct%vten  =0.;    ct%tten  =0.;
    ct%qvten =0.;    ct%qlten =0.;    ct%qiten =0.;
    ct%qaten =0.;    ct%qnten =0.;    ct%thcten=0.;
    ct%qctten=0.;    
    ct%qvdiv =0.;    ct%qldiv =0.;    ct%qidiv =0.;
    ct%thcflx=0.;    ct%qctflx=0.;    ct%qvflx =0.;
    ct%qlflx =0.;    ct%qiflx =0.;    ct%qaflx =0.;    ct%qnflx =0.;
    ct%umflx =0.;    ct%vmflx =0.;    ct%pflx  =0.;
    ct%hlflx =0.;    ct%hlten =0.;
    ct%denth =0.;    ct%rain  =0.;    ct%snow  =0.;    

  end subroutine ct_clear

!#####################################################################
!#####################################################################

  subroutine mixing(cpn, z0, p0, hl0, thc0, qct0, hlu, thcu, qctu, wu, &
       scaleh, rei, fer, fdr, fdrsat, rho0j, rkm)
  
    implicit none
    type(cpnlist),  intent(in)    :: cpn
    real,           intent(in)    :: z0, p0, hl0, thc0, qct0 !envirn. properties at level k
    real,           intent(in)    :: hlu, thcu, qctu, wu !updraft properties at level k-1
    real,           intent(in)    :: scaleh, rkm
    real,           intent(inout) :: rei, fer, fdr, fdrsat, rho0j
    
    integer :: id_check
    real    :: excessu, excess0, hlfs, thlfs, qtfs, thvfs, xbuo0, xsat, xs, xs1, xs2
    real    :: thj, qvj, qlj, qij, qse, thvj, thv0j
    real    :: aquad, bquad, cquad, ee2, ud2

!-----A.  Entrainment and Detrainment
!     first, to determine fraction (xsat) of mixture that is to be detrained out 
!     of clouds, i.e., the mixture with negative buoyancy. We consider a thin 
!     layer between two interfaces, so using mid-point value to represent the 
!     mean value of the layer. The properties of updraft at midpoint is assumed
!     to be undiluted from the lower interface.  

!-----calculate fraction of mixture that is just saturated

    !excessu = qctu - qsat(thcu*exn(p0), p0); excessu = max(excessu,0.0)
    !excess0 = qct0 - qsat(thc0*exn(p0), p0); 
    excessu = qctu - qsat((hlu-grav*z0)/cp_air, p0); excessu = max(excessu,0.0)
    excess0 = qct0 - qsat((hl0-grav*z0)/cp_air, p0); 

    if(excessu*excess0.le.0)then
       xsat = -excessu/(excess0-excessu)
    else
       xsat = 1.0
    endif
    hlfs =(1.-xsat)*hlu  + xsat*hl0
    !thlfs=(1.-xsat)*thcu + xsat*thc0
    qtfs =(1.-xsat)*qctu + xsat*qct0
    !thvfs=thlfs*(1+zvir*qtfs)
    call findt (z0,p0,hlfs,qtfs, thj, qvj, qlj, qij, qse, thvfs, cpn%do_ice)

    !call conden(p0, thcu, qctu, thj, qvj, qlj, qij, qse, thvj)
    call findt(z0,p0,hlu, qctu, thj, qvj, qlj, qij, qse, thvj, cpn%do_ice)
    !call conden(p0, thc0, qct0, thj, qvj, qlj, qij, qse, thv0j)
    call findt(z0,p0,hl0, qct0, thj, qvj, qlj, qij, qse, thv0j, cpn%do_ice)
    rho0j = p0/(rdgas*thv0j*exn(p0))

!-----calculate fraction of mixture with zero buoyancy
    if(thvfs.ge.thv0j) then
       xbuo0=xsat
    elseif(thvj.le.thv0j) then
       xbuo0=0.
    else
       !xbuo0=xsat*(thv0j-thvfs)/(thvj-thvfs)
       xbuo0=xsat*(thvj-thv0j)/(thvj-thvfs)
    endif

    !-----calculate fraction of mixture with negative buoyancy but can 
    !     penetrate a critical distance lc=rle*scaleh
    if(thvfs.ge.thv0j.or.xsat.le.0.05) then
       xs=xsat !mixture has to be saturated
    else
       aquad = wu**2.
       bquad = -(2.*wu**2. + 2.*cpn%rbuoy*grav*cpn%rle*scaleh*(thvj-thvfs)/thv0j/xsat)
       cquad = wu**2.      - 2.*cpn%rbuoy*grav*cpn%rle*scaleh*(1-thvj/thv0j)
       call roots(aquad,bquad,cquad,xs1,xs2)
       xs=min(xs1,xs2)
    endif
    xs=min(xs,xsat)
    xs=max(xbuo0,xs)
    xs=min(1.0,xs)
    
    ee2     = xs**2.
    ud2     = 1. - 2.*xs + xs**2.
    rei     = rkm/scaleh/grav/rho0j  !make entrainment rate in unit of 1/Pa
    fer     = rei * ee2
    fdr     = rei * ud2
    fdrsat = rei * (ud2-(1. - 2.*xsat + xsat**2.))

  end subroutine mixing

!#####################################################################
!#####################################################################

  subroutine cumulus_plume(cpn, sd, ac, cp, rkm, cbmf, wrel, scaleh)
  
    implicit none

    type(cpnlist),  intent(in)    :: cpn
    type(sounding), intent(in)    :: sd
    type(adicloud), intent(in)    :: ac
    real,           intent(in)    :: rkm, cbmf, wrel, scaleh
    type(cplume),   intent(inout) :: cp

    real, dimension(1:size(sd%p)) :: p0, dp
    real, dimension(0:size(sd%p)) :: ps0

    real, dimension(3)            :: totalmass
    real                          :: thickness, drop=0
   
    integer :: k, klm, km1, krel, let, ltop, id_check
    real    :: thv0rel, thv0t1, wexp, wtw
    real    :: thj, qvj, qlj, qij, qse, thvj, thv0j, rhos0j, rho0j
    real    :: aquad, bquad, cquad, xs1, xs2, ppen
    real    :: bogtop, bogbot, delbog, drage, expfac
    real    :: zrel, prel, exnj, nu, leff, qrj, qsj, tmp, temp
    real    :: qctu_new, hlu_new, qlu_new, qiu_new, clu_new, ciu_new
    real    :: auto_th

    call cp_clear(cp)
    cp%p=sd%p; cp%ps=sd%ps; cp%dp=sd%dp; cp%u=sd%u; cp%v=sd%v;
    cp%hl=sd%hl; cp%thc=sd%thc; cp%qct=sd%qct; 
    cp%ql=sd%ql; cp%qi=sd%qi; cp%qa=sd%qa; cp%qn=sd%qn;
    cp%thvbot=sd%thvbot; cp%thvtop=sd%thvtop;
    cp%z=sd%z; cp%zs=sd%zs;

    wtw  = wrel*wrel

    !determine release height and parcel properties (krel, prel, thv0rel, thvurel)
    if(ac % plcl .gt. sd % pinv)then
       krel    = sd % kinv
       prel    = sd % pinv
       zrel    = sd % zinv
       thv0rel = sd % thvinv
    else
       krel    = ac % klcl
       prel    = ac % plcl
       zrel    = ac % zlcl
       thv0rel = ac % thv0lcl
    endif

    cp%krel=krel
    cp%prel=prel
    cp%zrel=zrel

    !(krel-1) represents the bottom of the updraft
    !call conden(prel, ac%thcsrc, ac%qctsrc, thj, qvj, qlj, qij, qse, cp%thvu(krel-1))
    call findt(zrel,prel,ac%hlsrc,ac%qctsrc, thj, qvj, qlj, qij, qse, cp%thvu(krel-1), cpn%do_ice)
    cp%ps   (krel-1) = prel
    cp%hlu  (krel-1) = ac % hlsrc
    cp%thcu (krel-1) = ac % thcsrc
    cp%qctu (krel-1) = ac % qctsrc
    cp%uu   (krel-1) = ac % usrc
    cp%vu   (krel-1) = ac % vsrc
    cp%umf  (krel-1) = cbmf
    cp%wu   (krel-1) = wrel
    cp%ufrc (krel-1) = cp%umf(krel-1)/(sd%rho(krel-1)*cp%wu(krel-1))
    !==================================================
    !     yim's CONVECTIVE NUCLEATION
    !==================================================
    totalmass(1)=     sd%am1(krel-1); !totalmass(1)=aerol;
    totalmass(2)=0.1 *sd%am2(krel-1); !totalmass(2)=0.;
    totalmass(3)=1.67*sd%am3(krel-1); !totalmass(3)=0.;
    !call aer_ccn_act(thj*exn(prel), prel, wrel, totalmass, drop)
    cp%qnu(krel-1) = drop * 1.0e6 / (prel / (rdgas*cp%thvu(krel-1)*exn(prel)))


    !(krel) represents the first partial updraft layer
    cp%z      (krel) = (cp%zs(krel-1) + cp%zs(krel))*0.5
    cp%p      (krel) = (cp%ps(krel-1) + cp%ps(krel))*0.5
    cp%dp     (krel) =  cp%ps(krel-1) - cp%ps(krel)   
    cp%thvbot (krel) = thv0rel
    if(krel.ne. sd % kinv) then
       cp%hl  (krel) = cp%hl (krel)+sd%sshl (krel)*(cp%p(krel)-sd%p(krel))
       cp%thc (krel) = cp%thc(krel)+sd%ssthc(krel)*(cp%p(krel)-sd%p(krel))
       cp%qct (krel) = cp%qct(krel)+sd%ssqct(krel)*(cp%p(krel)-sd%p(krel))
       !call conden(cp%p(krel), cp%thc(krel), cp%qct(krel), thj, qvj, qlj, qij, qse, cp%thvu(krel))
       call findt(cp%z(krel),cp%p(krel),cp%hl(krel), cp%qct(krel), thj, qvj, &
            qlj, qij, qse, cp%thvu(krel), cpn%do_ice)
    endif

    !Compute updraft properties above the LCL
    let=krel
    klm=sd%kmax-1
    do k=krel,klm
       km1=k-1

       !Calculation entrainment and detrainment rate
       if (cpn%do_edplume) then
          call mixing(cpn, cp%z(k), cp%p(k), cp%hl(k), cp%thc(k), cp%qct(k), &
               cp%hlu(km1), cp%thcu(km1), cp%qctu(km1), cp%wu(km1), scaleh,  &
               cp%rei(k), cp%fer(k), cp%fdr(k), cp%fdrsat(k), rho0j, rkm)
       else
          temp         = sqrt(cp%ufrc(km1)) !scaleh for fixed length scale
          rho0j        = sd%rho(k)
          cp%rei(k)    = rkm/temp/grav/rho0j
          cp%fer(k)    = cp%rei(k)
          cp%fdr(k)    = 0.
          cp%fdrsat(k) = 0.
       end if

       !Calculate the mass flux
       cp%umf(k)=cp%umf(km1)*exp(cp%dp(k)*(cp%fer(k)-cp%fdr(k)))
       cp%emf(k)=0.0

       !Thermodynamics for the dilute plume
       cp%hlu (k)=cp%hl (k)-(cp%hl (k)-cp%hlu (km1))*exp(-cp%fer(k)*cp%dp(k))
       !cp%thcu(k)=cp%thc(k)-(cp%thc(k)-cp%thcu(km1))*exp(-cp%fer(k)*cp%dp(k))
       cp%qctu(k)=cp%qct(k)-(cp%qct(k)-cp%qctu(km1))*exp(-cp%fer(k)*cp%dp(k))
       cp%qnu (k)=cp%qn (k)-(cp%qn (k)-cp%qnu (km1))*exp(-cp%fer(k)*cp%dp(k))
       if(cp%fer(k)*cp%dp(k).lt.1.e-4)then
          cp%uu(k)=cp%uu(km1) - sd%dudp(k)*cp%dp(k)
          cp%vu(k)=cp%vu(km1) - sd%dvdp(k)*cp%dp(k)
       else
          cp%uu(k)=cp%u(k)-cpn%bigc*sd%dudp(k)/cp%fer(k)-exp(-cp%fer(k)*cp%dp(k))* &
                  (cp%u(k)-cpn%bigc*sd%dudp(k)/cp%fer(k) - cp%uu(km1))
          cp%vu(k)=cp%v(k)-cpn%bigc*sd%dvdp(k)/cp%fer(k)-exp(-cp%fer(k)*cp%dp(k))* &
                  (cp%v(k)-cpn%bigc*sd%dvdp(k)/cp%fer(k) - cp%vu(km1))
       endif

       if (cpn%do_micro) then
          call micro_donner (cpn, cp%zs(k), cp%ps(k), cp%hlu(k), cp%qctu(k),               &
               cp%zs(km1), cp%qlu(km1), cp%clu(km1), cp%qiu(km1), cp%ciu(km1), cp%wu(km1), &
               cp%crate(k), cp%prate(k),                                                   &
               qrj, qsj, qlu_new, clu_new, qiu_new, ciu_new, hlu_new, qctu_new, temp, cpn%do_ice)
       else
          call precipitation(cp%zs(k), cp%ps(k), cp%hlu(k), cp%qctu(k), &
               cpn, qrj, qsj, hlu_new, qctu_new, qlu_new, qiu_new, temp, cpn%do_ice)
       clu_new = 0.0 ! Need to initialize these variables to zero for now.
       ciu_new = 0.0 ! A new version of precipitation will be included after this release.
       end if

       cp%qctu(k)=qctu_new
       cp%hlu (k)=hlu_new
       cp%qlu (k)=qlu_new
       cp%qiu (k)=qiu_new
       cp%clu (k)=clu_new
       cp%ciu (k)=ciu_new
       cp%thvu(k)=temp/exn(cp%ps(k))*(1.+zvir*(cp%qctu(k)-cp%qlu(k)-cp%qiu(k))-cp%qlu(k)-cp%qiu(k))
       cp%buo (k)=cp%thvu(k)-cp%thvtop(k)
       cp%t   (k)=temp
       nu = max(min((268. - temp)/20.,1.0),0.0)
       leff = (1-nu)*HLv + nu*HLs
       !cp%thcu(k)=cp%thcu(k) + (qrj + qsj)*leff/cp_air/exn(cp%ps(k))
       cp%thcu(k)=temp/exn(cp%ps(k))

       !Calculate vertical velocity
!!$       bogbot = (cp%thvu(km1)/cp%thvbot(k) - 1.)
!!$       if(bogbot.gt.0.)then
!!$          bogbot =  bogbot /cpn%rbuoy
!!$       else
!!$          bogbot =  bogbot*cpn%rbuoy
!!$       endif
!!$       bogtop = (cp%thvu(k)/cp%thvtop(k) - 1.)
!!$       if(bogtop.gt.0.)then
!!$          bogtop =  bogtop /cpn%rbuoy
!!$       else
!!$          bogtop =  bogtop*cpn%rbuoy
!!$       endif

       bogbot = (cp%thvu(km1)/cp%thvbot(k) - 1.)
       bogbot =  bogbot*cpn%rbuoy
       bogtop = (cp%thvu(k)/cp%thvtop(k) - 1.)
       bogtop =  bogtop*cpn%rbuoy

       if(bogbot.gt.0.and.bogtop.gt.0) let = k
       delbog = bogtop - bogbot
       drage = cp%fer(k) * ( 1. + cpn%rdrag )
       expfac = exp(-2.* drage * cp%dp(k))
       if(drage * cp%dp(k).gt.1.e-3) then
          wtw = wtw*expfac + (delbog + (1.-expfac) *               &
               (bogbot+delbog/(-2.*drage*cp%dp(k))))/(rho0j*drage)
       else
          wtw = wtw + cp%dp(k) * (bogbot+bogtop)/rho0j
       endif

       if (cpn%do_forcedlifting) then
          if(cp%buo(k).le.0. .and. k <= ac%klfc) then
             wtw=wrel*wrel
          elseif (wtw.le.0.) then
             exit
          end if
       else
          if(wtw.le.0.) exit
       end if

       cp%wu(k) = sqrt(wtw)
       if(cp%wu(k).gt.100.)then
          print *, 'Very big wu in UW-ShCu',bogbot,bogtop,expfac,cp%fer(k)
          stop
       endif
       
       rhos0j     = cp%ps(k)/(rdgas*0.5*(cp%thvbot(k+1)+cp%thvtop(k))*exn(cp%ps(k)))
       cp%ufrc(k) = cp%umf(k)/(rhos0j*cp%wu(k))
       if(cp%ufrc(k).gt.cpn%rmaxfrac)then
          cp%ufrc(k) = cpn%rmaxfrac
          cp%fdr (k) = cp%fer(k)-log(cpn%rmaxfrac*rhos0j*cp%wu(k)/cp%umf(km1))/cp%dp(k)
          cp%umf (k) = cpn%rmaxfrac * rhos0j * cp%wu(k)
       endif
       cp%pptr(k) = qrj*cp%umf(k)
       cp%ppti(k) = qsj*cp%umf(k)
       
    enddo !End of Updraft Loop

    cp%let =let
    cp%ltop=k
    ltop = k !ltop is the 1st level with negative vertical velocity at top
    cp%umf(ltop)=0.0
    cp%cldhgt=sd%z(ltop)-ac%zlcl
    if (qrj+qsj .gt. 0) then
       cp%pptr(ltop) = 0.
       cp%ppti(ltop) = 0.
       cp%qlu (ltop) = cp%qlu(ltop)+qrj
       cp%qiu (ltop) = cp%qiu(ltop)+qsj
       cp%qctu(ltop) = cp%qctu(ltop) + qrj + qsj
       cp%hlu (ltop) = cp%hlu (ltop) - leff*(qrj+qsj)
       !cp%thcu(ltop) = cp%thcu(ltop) - leff*(qrj+qsj)/cp_air/exn(cp%ps(k))
   end if

    !Restriction of convection too deep or too shallow
    if(cp%cldhgt.ge.cpn%cldhgt_max .or. ltop.lt.krel+2 .or. let.le.krel+1) then
       return
    endif

    !convective scale height
    cp % cush=sd%z(ltop)

    if (cpn%do_ppen) then !Calculate penetrative entrainment
       call penetrative_mixing(cpn, sd, cp)
    end if

  end subroutine cumulus_plume


  subroutine precipitation(zs, ps, hlu, qctu, cpn, qrj, qsj, hlu_new, qctu_new, qlu_new, qiu_new, temp, doice)
    implicit none
    type(cpnlist),  intent(in)    :: cpn
    real,           intent(in)    :: zs, ps, hlu, qctu
    real,           intent(inout) :: qrj, qsj, hlu_new, qctu_new, qlu_new, qiu_new, temp
    logical,        intent(in)    :: doice

    real    :: thj, qvj, qlj, qij, qse, thvj, thv0j, nu, exnj, auto_th, leff

    !Precip at the flux level
    !call conden(ps,thcu,qctu,thj,qvj,qlj,qij,qse,thvj)
    call findt(zs,ps,hlu,qctu,thj,qvj,qlj,qij,qse,thvj,doice)
    exnj=exn(ps)
    temp=thj*exnj-273.15
    if (temp.ge.0.0) then
       auto_th=cpn%auto_th0
    else
       auto_th=cpn%auto_th0*(1.0-temp/cpn%tcrit)
    end if
    auto_th=max(auto_th,0.0)

    temp=temp+273.15
    if((qlj+qij).gt.auto_th)then
       qrj = (qlj+qij-auto_th)*qlj/(qlj+qij)
       qsj = (qlj+qij-auto_th)*qij/(qlj+qij)
       nu = max(min((268. - temp)/20.,1.0),0.0)
     else
       qrj = 0.0
       qsj = 0.0
       nu  = 0.0
    endif
    leff     = (1-nu)*HLv + nu*HLs
    qctu_new = qctu - (qrj + qsj)
    hlu_new  = hlu  + (qrj + qsj)*leff
!    thcu_new = thcu + (qrj + qsj)*leff/cp_air/exnj
    qlu_new  = qlj - qrj
    qiu_new  = qij - qsj

    return
    
  end subroutine precipitation



  subroutine micro_donner(cpn, zs, ps, hlu, qctu, zs1, qlu1, clu1, qiu1, ciu1, w1, cr12, pr12, &
       qrj, qsj, qlu_new, clu_new, qiu_new, ciu_new, hlu_new, qctu_new, temp, doice)
    implicit none
    type(cpnlist),  intent(in)    :: cpn
    real,           intent(in)    :: zs, ps, hlu, qctu
    real,           intent(in)    :: zs1, qlu1, clu1, qiu1, ciu1, w1
    real,           intent(inout) :: qrj, qsj, cr12, pr12
    real,           intent(inout) :: qlu_new, clu_new, qiu_new, ciu_new, hlu_new, qctu_new, temp
    logical,        intent(in)    :: doice

    real    :: thj, qvj, qlj, qij, qse, thvj, thv0j, exnj, nu, leff
    real    :: dt_micro, rw1, cw1, drwa, drwb, flw, rw2, cw2, pw2, dcw

    call findt(zs,ps,hlu,qctu,thj,qvj,qlj,qij,qse,thvj,doice)
    temp = thj*exn(ps)
    nu   = max(min((268. - temp)/20.,1.0),0.0)
    leff = (1.-nu)*HLv + nu*HLs
    if (qlj+qij .gt. 0.0) then
       flw  = qlj/(qlj+qij)
    else
       qrj      = 0.
       qsj      = 0.
       qlu_new  = 0.
       qiu_new  = 0.
       qctu_new = qctu
       hlu_new  = hlu
       clu_new  = 0.
       ciu_new  = 0.
       return
    end if

    cw1 = clu1 + ciu1
    rw1 = qlu1 + qiu1 - cw1

    cw2 = qlj + qij - rw1

    dcw = cw2 - cw1; 

    dt_micro = (zs - zs1) / w1

    cr12 = dcw/dt_micro

    drwa = cpn%auto_rate * (cw2 - cpn%auto_th0) * dt_micro

    drwa=min(max(drwa, 0.0), cw2-cpn%auto_th0)

    cw2 = cw2 - drwa

    drwb = 5.26e-03 * cw2 * (rw1**0.875) * dt_micro

    drwb=min(max(drwb, 0.0), cw2)

    cw2 = cw2 - drwb

    rw2 = rw1 + drwa + drwb

    rw2 = max(rw2, 0.0)

    pw2 = 5.1*(rw2**1.125)*dt_micro/100.

    pw2 =min(pw2, rw2)

!    pw2 =min(pw2, max(dcw,0.))

    pr12=pw2/dt_micro
   
    rw2 =rw2-pw2

    qrj = pw2*flw
    qsj = pw2*(1.-flw)

    qlu_new  = qlj - qrj
    qiu_new  = qij - qsj
    qctu_new = qctu - (qrj + qsj)
    hlu_new  = hlu  + (qrj + qsj)*leff

    clu_new  = (qlu_new-rw2)*flw
    ciu_new  = (qiu_new-rw2)*(1.-flw)

    if (qctu_new .le. 0.) then
       print*, qctu_new, qrj, qsj, clu_new, ciu_new,'??????????????????'
    end if
    return
    
  end subroutine micro_donner



!#####################################################################
!#####################################################################

  subroutine penetrative_mixing(cpn, sd, cp)
    implicit none
    type(cpnlist),  intent(in)    :: cpn
    type(sounding), intent(in)    :: sd
    type(cplume),   intent(inout) :: cp

    integer :: k, ltop, let
    real    :: rhos0j, bogtop, bogbot
    real    :: aquad, bquad, cquad, xs1, xs2, ppen

    ltop=cp%ltop
    let =cp%let

    cp % emf(ltop)=0.0
    do k=ltop-1,let,-1
       rhos0j = cp%ps(k) /(rdgas*0.5*(cp%thvbot(k+1)+cp%thvtop(k))*exn(cp%ps(k)))
       if(k.eq.ltop-1)then
          !Calculate ppen
!!$          bogbot = (cp%thvu(k)/cp%thvbot(ltop) - 1.)/cpn%rbuoy
!!$          if(bogbot.gt.0.)then
!!$             bogbot =  bogbot /cpn%rbuoy
!!$          else
!!$             bogbot =  bogbot*cpn%rbuoy
!!$          endif
!!$          bogtop = (cp%thvu(ltop)/cp%thvtop(ltop) - 1.)/cpn%rbuoy
!!$          if(bogtop.gt.0.)then
!!$             bogtop =  bogtop /cpn%rbuoy
!!$          else
!!$             bogtop =  bogtop*cpn%rbuoy
!!$          endif

          bogbot = (cp%thvu(k)/cp%thvbot(ltop) - 1.)
          bogbot =  bogbot*cpn%rbuoy
          bogtop = (cp%thvu(ltop)/cp%thvtop(ltop) - 1.)
          bogtop =  bogtop*cpn%rbuoy

          aquad = (bogtop - bogbot) / (cp%ps(ltop)-cp%ps(k))
          bquad = 2*bogbot
          cquad = -cp%wu(k) * cp%ps(k) / (rdgas*cp%thvbot(ltop)*exn(cp%ps(k)))
          call roots(aquad,bquad,cquad,xs1,xs2)
          if(xs1.le.0..and.xs2.le.0.)then
             ppen = max(xs1,xs2)
          else
             ppen = min(xs1,xs2)
          endif
          ppen = min(0.,max(-cp%dp(k+1),ppen))
          if(xs1.eq.-9.99e33.or.xs2.eq.-9.99e33) ppen=0.
          !Calculate returning mass flux
          cp%emf (k)=max(cp%umf(k)*ppen*cp%rei(ltop)*cpn%rpen,-0.1*rhos0j)
          cp%hlu (k)=cp%hl (ltop)+sd%sshl (ltop)*(cp%ps(k)-cp%p(ltop))
          !cp%thcu(k)=cp%thc(ltop)+sd%ssthc(ltop)*(cp%ps(k)-cp%p(ltop))
          cp%qctu(k)=cp%qct(ltop)+sd%ssqct(ltop)*(cp%ps(k)-cp%p(ltop))
       else
          cp%emf (k)=max(cp%emf(k+1)-cp%umf(k)*cp%dp(k+1)*cp%rei(k+1)*cpn%rpen,-0.1*rhos0j)
          cp%hlu (k)=(cp%hlu (k+1)*cp%emf(k+1)+cp%hl (k+1)*(cp%emf(k)-cp%emf(k+1)))/cp%emf(k)
          !cp%thcu(k)=(cp%thcu(k+1)*cp%emf(k+1)+cp%thc(k+1)*(cp%emf(k)-cp%emf(k+1)))/cp%emf(k)
          cp%qctu(k)=(cp%qctu(k+1)*cp%emf(k+1)+cp%qct(k+1)*(cp%emf(k)-cp%emf(k+1)))/cp%emf(k)
       endif
       cp%umf(k)=0.0 !the line is commented out by pzhu
    enddo
  end subroutine penetrative_mixing

!#####################################################################
!#####################################################################

  subroutine cumulus_tend(cpn, sd, cp, ct, do_coldT)
  
    implicit none
    type(cpnlist),  intent(in)    :: cpn
    type(sounding), intent(in)    :: sd
    type(cplume),   intent(inout) :: cp
    type(ctend),    intent(inout) :: ct
    logical,        intent(in)    :: do_coldT

    integer :: k, krel, ltop, kp1, km1, ktop
    real    :: dpsum, qtdef, umftmp, qlutmp, qiutmp, fdrtmp

    call ct_clear(ct);

    krel=cp%krel
    ltop=cp%ltop

    ! Calculate Fluxes of heat, moisture, momentum
    dpsum = 0.0
    do k = 1, krel-1
       dpsum = dpsum + sd%dp(k)
    enddo

    qtdef = max(0.,cp%umf(krel)*(cp%qctu(krel) - cp%qct(krel)))
    !yy1  = min(0.,umf(krel)*(thcu(krel) - thc0(krel)))
    do k=1,krel-1
       ct%hlflx (k)=0.0; 
       ct%thcflx(k)=0.0; !thcflx(k)=thcflx(k-1) + yy1*dp(k)/dpsum
       ct%qctflx(k)=ct%qctflx(k-1) + qtdef*sd%dp(k)/dpsum;
       ct%qlflx(k)=0.0;
       ct%qiflx(k)=0.0;
       ct%umflx(k)=0.0;
       ct%vmflx(k)=0.0;
       cp%pptr (k)=0.0; 
       cp%ppti (k)=0.0;
    enddo

    do k = krel,ltop-1 !pzhu do k = krel,ltop
       kp1 = k+1

       ct%hlflx (k)= cp%umf(k)*(cp%hlu (k)-(cp%hl (kp1)+sd%sshl (kp1)*(sd%ps(k)-sd%p(kp1)))) + &
            cp%emf(k) * (cp%hlu (k)-(cp%hl (k)+sd%sshl (k)*(sd%ps(k)-sd%p(k))))
       ct%thcflx(k)= cp%umf(k)*(cp%thcu(k)-(cp%thc(kp1)+sd%ssthc(kp1)*(sd%ps(k)-sd%p(kp1)))) + &
            cp%emf(k) * (cp%thcu(k)-(cp%thc(k)+sd%ssthc(k)*(sd%ps(k)-sd%p(k))))
       ct%qctflx(k)= cp%umf(k)*(cp%qctu(k)-(cp%qct(kp1)+sd%ssqct(kp1)*(sd%ps(k)-sd%p(kp1)))) + &
            cp%emf(k) * (cp%qctu(k)-(cp%qct(k)+sd%ssqct(k)*(sd%ps(k)-sd%p(k))))
       
       ct%umflx(k) =cp%umf(k) * (cp%uu(k) - cp%u(kp1))  + cp%emf(k) * (cp%u(kp1) -cp%u(k))
       ct%vmflx(k) =cp%umf(k) * (cp%vu(k) - cp%v(kp1))  + cp%emf(k) * (cp%v(kp1) -cp%v(k))

       ct%qlflx(k) =cp%umf(k) * (cp%qlu(k)- cp%ql(kp1)) + cp%emf(k) * (cp%ql(kp1)-cp%ql(k))
       ct%qiflx(k) =cp%umf(k) * (cp%qiu(k)- cp%qi(kp1)) + cp%emf(k) * (cp%qi(kp1)-cp%qi(k))
       ct%qaflx(k) =cp%umf(k) * (1.       - cp%qa(kp1)) + cp%emf(k) * (cp%qa(kp1)-cp%qa(k))
       ct%qnflx(k) =cp%umf(k) * (cp%qnu(k)- cp%qn(kp1)) + cp%emf(k) * (cp%qn(kp1)-cp%qn(k))

       ct%qvflx(k) =ct%qctflx(k)-ct%qlflx(k)-ct%qiflx(k)
    enddo

    do k = 0,krel
       ct%thcflx(k)=ct%thcflx(krel)
       ct%qlflx (k)=ct%qlflx (krel)
       ct%qiflx (k)=ct%qiflx (krel)
       ct%qaflx (k)=ct%qaflx (krel)
       ct%qnflx (k)=ct%qnflx (krel)
    end do

    ! Calculate model tendencies
    ct%denth=0.; ct%uav=0.;ct%vav=0.; dpsum=0.;
    do k = 1,ltop
       km1 = k-1
       kp1 = k+1
       ct%uten (k) = (ct%umflx(km1)-ct%umflx(k))*grav/sd%dp(k)
       ct%vten (k) = (ct%vmflx(km1)-ct%vmflx(k))*grav/sd%dp(k)

       umftmp = cp%umf(k) !(umf(km1)*(p0(k)-ps0(k))+umf(k)*(ps0(km1)-p0(k)))/(ps0(km1)-ps0(k))
       qlutmp = (cp%qlu(km1)*(sd%p(k)-sd%ps(k))+cp%qlu(k)*(sd%ps(km1)-sd%p(k)))/(sd%ps(km1)-sd%ps(k))
       qiutmp = (cp%qiu(km1)*(sd%p(k)-sd%ps(k))+cp%qiu(k)*(sd%ps(km1)-sd%p(k)))/(sd%ps(km1)-sd%ps(k))
       fdrtmp = min(cp%fdrsat(k),cp%fdr(k))*cpn%frac_drs
       
       ct%qlten(k) = grav*umftmp*(fdrtmp*(qlutmp-sd%ql(k)) + (sd%ql(k)-sd%ql(kp1))/sd%dp(k))
       ct%qiten(k) = grav*umftmp*(fdrtmp*(qiutmp-sd%qi(k)) + (sd%qi(k)-sd%qi(kp1))/sd%dp(k))
       ct%qaten(k) = grav*umftmp*(fdrtmp*(1.    -sd%qa(k)) + (sd%qa(k)-sd%qa(kp1))/sd%dp(k))
       
       ct%qlten(k) = max(ct%qlten(k),0.)
       ct%qiten(k) = max(ct%qiten(k),0.)
       ct%qaten(k) = max(ct%qaten(k),0.)

       ct%hlten (k) = (ct%hlflx (km1) - ct%hlflx (k))* grav/sd%dp(k)
       ct%thcten(k) = (ct%thcflx(km1) - ct%thcflx(k))* grav/sd%dp(k)
       ct%qctten(k) = (ct%qctflx(km1) - ct%qctflx(k) - cp%pptr(k) - cp%ppti(k)) *grav/sd%dp(k)
       ct%qvten (k) = ct%qctten(k)-ct%qlten(k)-ct%qiten(k)
       ct%pflx  (k) = cp%pptr(k) + cp%ppti(k)

!!$       ct%tten  (k) = ct%thcten(k)+sd%leff(k)* &
!!$                     (ct%qlten(k)+ct%qiten(k)+(cp%pptr(k)+cp%ppti(k))*grav/sd%dp(k))/cp_air

       ct%tten  (k) = ct%hlten(k)/cp_air+sd%leff(k)* &
                     (ct%qlten(k)+ct%qiten(k)+(cp%pptr(k)+cp%ppti(k))*grav/sd%dp(k))/cp_air

       ct%qnten(k) = (ct%qnflx(km1)-ct%qnflx(k))*grav/sd%dp(k)

       ct%denth = ct%denth + (cp_air*ct%tten(k)+sd%leff(k)*ct%qvten(k))*sd%dp(k)/cp_air
       ct%uav   = ct%uav   + ct%uten(k)*sd%dp(k)
       ct%vav   = ct%vav   + ct%vten(k)*sd%dp(k)
       dpsum    = dpsum    + sd%dp(k)
    enddo
    ct%denth = ct%denth/dpsum
    ct%uav   = ct%uav  /dpsum
    ct%vav   = ct%vav  /dpsum


    ct%dtint=0.; ct%dqint=0.; ct%conint=0.; ct%freint=0.; !dpsum=0.;

    if (cpn%do_topevap) then
       ktop = ltop
    else
       ktop = ltop-1
    end if
    do k = 1,ktop
       km1 = k-1
       ct%qvdiv(k) = -(ct%qvflx(k) - ct%qvflx(km1))* grav/sd%dp(k)
       ct%qldiv(k) = -(ct%qlflx(k) - ct%qlflx(km1) + cp%pptr(k))* grav/sd%dp(k)
       ct%qidiv(k) = -(ct%qiflx(k) - ct%qiflx(km1) + cp%ppti(k))* grav/sd%dp(k)
       ct%dtint    = ct%dtint    + ct%tten (k)*cp_air/sd%leff(k)* sd%dp(k)/grav
       ct%dqint    = ct%dqint    + ct%qvten(k)                  * sd%dp(k)/grav
       ct%conint   = ct%conint   + ct%qldiv(k)                  * sd%dp(k)/grav
       ct%freint   = ct%freint   + ct%qidiv(k)                  * sd%dp(k)/grav
       !ct%freint   = ct%freint   + cp%crate(k)*cp%ufrc(k)*sd%dp(k)/grav
       !dpsum       = dpsum       + sd%dp(k)
    end do

    !slightly adjust tendencies to force exact enthalpy and momentum conservation
    do k = 1,ltop
       ct%tten(k)=ct%tten(k)-ct%denth
       ct%uten(k)=ct%uten(k)-ct%uav
       ct%vten(k)=ct%vten(k)-ct%vav
    end do

    if (do_coldT) then
       do k = 1,ltop
          if (sd%coldT) then
             cp%ppti(k)=cp%ppti(k)+cp%pptr(k)
             cp%pptr(k)=0.
             ct%tten(k)=ct%tten(k)+cp%pptr(k)*HLf*grav/sd%dp(k)/cp_air
          else
             cp%pptr(k)=cp%pptr(k)+cp%ppti(k)
             cp%ppti(k)=0.
             ct%tten(k)=ct%tten(k)-cp%ppti(k)*HLf*grav/sd%dp(k)/cp_air
          end if
          ct%snow  = ct%snow  + cp%ppti(k)
          ct%rain  = ct%rain  + cp%pptr(k)
       end do
    end if

  end subroutine cumulus_tend

!#####################################################################
!#####################################################################

 subroutine roots(a,b,c,r1,r2)
   implicit none
   real a,b,c,r1,r2,q
   
   if(a.eq.0)then            ! form b*x + c = 0
      if(b.eq.0)then         ! failure: c = 0
         r1 = -9.99e33
      else                   ! b*x + c = 0
         r1 = -c / b
      endif
      r2 = r1
   else
      if(b.eq.0.)then        ! form a*x**2 + c = 0
         if(a*c.gt.0.)then   ! failure: x**2 = -c/a < 0
            r1 =  -9.99e33
         else                ! x**2 = -c/a
            r1 = sqrt(-c/a)
         endif
         r2 = -r1
      else 
         if((b**2. - 4.*a*c).lt.0.)then ! failure, no real(r8) roots
            r1 =  -9.99e33
            r2 = -r1
         else
            q = - 0.5 * ( b + sign(1.0,b) * sqrt(b**2. - 4.*a*c) )
            r1 = q/a
            r2 = c/q
         endif
      endif
   endif
   return
 end subroutine roots

!#####################################################################
!#####################################################################

  subroutine cumulus_downdraft(sd, cp)
  
    implicit none

    type(sounding), intent(in)    :: sd
    type(cplume),   intent(inout) :: cp

    integer :: k, klm, km1, krel, let, ltop, jtt
    real    :: coeff, rat, hk, hkm1
    real    :: wdetrain, qsm, qstm, afac, B6, C6, revap, dhdp, fac, precip
    real, dimension(1:size(sd%p)) :: wt, evap, water, qp, mp, up, vp

    real :: coeffr=1.0, coeffs=0.8, omtsnow=5.5, omtrain=50.0, sigd=0.05, sigs=0.12
    
    qp=0.; water=0.; evap =0.; mp=0.; up=0.; vp=0.; qp=0.; wt=omtrain
    precip=0.;
    do k = cp%ltop-1, 1, -1
       if(sd%coldT)then
          wdetrain=cp%ppti(k)
          coeff=coeffs
          wt(k)=omtsnow
       else
          wdetrain=cp%pptr(k)
          coeff=coeffr
          wt(k)=omtrain
       end if

       !temporary evaluation of specific humdity of precipitation downdraft
       qsm=0.5*(cp%qct(k)-cp%qlu(k)-cp%qiu(k)+qp(k+1)) 

       afac=coeff*sd%ps(k)*(sd%qs(k)-qsm)/(1.0e4+2.0e3*sd%ps(k)*sd%qs(k))
       afac=max(afac,0.0)
       B6=100.*(sd%ps(k)-sd%ps(k+1))*sigs*afac/wt(k)
       C6=(water(k+1)*wt(k+1)+wdetrain/sigd)/wt(k)
       revap=0.5*(-B6+SQRT(B6*B6+4.*C6))
       evap (k)=sigs*afac*revap
       water(k)=revap*revap

       !CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER
       !HYDROSTATIC APPROXIMATION

       if(k.eq.1) goto 360
       hk  =CP_air*sd%t(k)  +GRAV*sd%z(k)  +HLv*sd%qv(k)
       hkm1=CP_air*sd%t(k-1)+GRAV*sd%z(k-1)+HLv*sd%qv(k-1)
       dhdp =(hk-hkm1)/sd%dp(k)
       dhdp =max(dhdp,10.0)
       mp(k)=100./grav*HLv*sigd*evap(k)/dhdp
       mp(k)=max(mp(k),0.0)

       !ADD SMALL AMOUNT OF INERTIA TO DOWNDRAFT 
       fac  =20.0/sd%dp(k)
       mp(k)=(fac*mp(k+1)+mp(k))/(1.+fac)

       !FORCE MP TO DECREASE LINEARLY TO ZERO
       !BETWEEN ABOUT 950 MB AND THE SURFACE

       if(sd%p(k).gt.(0.949*sd%p(1))) then
          jtt  =max(jtt,k)
          mp(k)=mp(jtt)*(sd%p(1)-sd%p(k))/(sd%p(1)-sd%p(jtt))
       end if
360    continue

       !FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT

       if(k.eq.cp%ltop) goto 400
       if(k.eq.1) then
          qstm=sd%qs(1)
       else
          qstm=sd%qs(k-1)
       end if
       if(mp(k).gt.mp(k+1)) then
          rat=mp(k+1)/mp(k)
          qp(k)=qp(k+1)*rat+sd%qv(k)*(1.0-rat)+100./grav* &
               SIGD*(sd%ps(k)-sd%ps(k+1))*(evap(k)/mp(k))
          up(k)=up(k+1)*rat+sd%u(k)*(1.-rat)
          vp(k)=vp(k+1)*rat+sd%v(k)*(1.-rat)
       else
          if(mp(k+1).gt.0.0) then
             qp(k)=qp(k+1)
             !(grav*sd%z(k+1)-grav*sd%z(k) + qp(k+1)*(HLv+sd%t(k+1)
             !*(CL-CPD))+CPD*(sd%t(k+1)-sd%t(k)))/(HLv+sd%t(k)*(CL-CPD))
             up(k)=up(k+1)
             vp(k)=vp(k+1)
          end if
       end if
       qp(k)=min(qp(k),qstm)
       qp(k)=max(qp(k),0.0)
400    continue

       precip=precip+wt(1)*SIGD*water(1)*3600.*24000./(1000.*GRAV)
405    continue
       
!       WD=BETA*ABS(mp(ICB))*0.01*RD*T(ICB)/(SIGD*P(ICB))
!       QPRIME=0.5*(QP(1)-Q(1))
!       TPRIME=LV0*QPRIME/CPD
    end do

  end subroutine cumulus_downdraft

!#####################################################################

!#####################################################################


end MODULE CONV_PLUMES_MOD
