  MODULE RAS_MOD

!=======================================================================
!  RELAXED ARAKAWA/SCHUBERT (RAS) CUMULUS PARAM SCHEME MODULE          !
!=======================================================================

 use  Sat_Vapor_Pres_Mod, ONLY: ESCOMP, DESCOMP
 use Utilities_Mod,       ONLY: FILE_EXIST, OPEN_FILE, ERROR_MESG,  &
                                get_my_pe, CLOSE_FILE, FATAL
 USE  Constants_Mod,      ONLY:  HLv,HLs,Cp,Grav,Kappa,rdgas,rvgas

!---------------------------------------------------------------------
 implicit none
 private
!---------------------------------------------------------------------

 public  :: ras, ras_init, ras_bdgt

!---------------------------------------------------------------------

!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 character(len=128) :: version = '$Id: ras.F90,v 1.4 2001/03/06 18:51:50 fms Exp $'
 character(len=128) :: tag = '$Name: damascus $'
!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 real :: cp_div_grav
 real :: one_plus_kappa, one_minus_kappa
 real :: onebcp, onebg
 real :: rn_pfac

 logical :: do_init = .true.

 real, parameter :: d622 = rdgas/rvgas
 real, parameter :: d378 = 1.0-d622

!---------------------------------------------------------------------
! --- Climatological cloud work function data
!---------------------------------------------------------------------
 
 real                :: actop 
 real, dimension(15) :: ac, ad
 real, dimension(15) :: ph, a

 data ph / 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0,   &
           550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0 /

 data a / 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677,            &
          0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664,            &
          0.0553, 0.0445, 0.0633  /

!---------------------------------------------------------------------
! --- NAMELIST
!---------------------------------------------------------------------
!     fracs  - Fraction of PBL mass allowed to be used 
!              by a cloud-type in time DT
!     rasal0 - Base value for cloud type relaxation parameter
!     puplim - Upper limit for cloud tops of deepest clouds
!     aratio - Ratio of critical cloud work function to standard
!              value of cloud work function
!     cufric - Should Cumulus friction (momentum transport) occur?
! --- CLOUD ORDER SPECIFICATION ---
!     ncrnd  - Number of random cloud-types between krmin and krmax to be
!              invoked in a single call to ras
!     iseed  - Integer seed used in generating random numbers
!              -- used only when ncrnd > 0
!     krmin  - Index of the top most level to which random clouds may
!              be invoked
!     krmax  - Index of the bottom most level to which random clouds 
!              may be invoked. krmin should be <= krmax. 
!              If ncrnd is specified as zero, then all cloud-types 
!              below the level krmax will be called sequentially.
!     botop  - A logical variable -- .true. if sequential clouds are 
!              called from bottom to top and .false. if top to bottom.
!              Level krmax will be called sequentially.
! --- PARTION CLOUD LIQUID WATER INTO PRECIP & DETRAINED VAPOR AND LIQUID --- 
!     rn_ptop - rn_frac_top of parcel liquid water converted to precip 
!               for cloud top pressures above rn_ptop
!     rn_pbot - rn_frac_bot of parcel liquid water converted to precip 
!               for cloud top pressures below rn_pbot (linear profile in 
!               between)
!     rn_frac_bot - Fraction of parcel liquid water converted to 
!                   precip for cloud top pressures below rn_pbot
!     rn_frac_top - Fraction of liquid water converted to
!                   precip for cloud top pressures above rn_ptop
! --- EVAPORATION OF CONVECTIVE SCALE PRECIP ---
!     evap_on - Turn on evaporation if true
!     cfrac   - Fraction of grid box assumed to be covered by convection
!     hcevap  - Evap allowed while q <= hcevap * qsat
!---------------------------------------------------------------------

 real :: fracs  = 0.25
 real :: rasal0 = 0.25
 real :: puplim = 20.0E2
 real :: aratio = 1.4
 logical :: cufric = .false.

! --- cloud order specification ---
 integer :: ncrnd = 0
 integer :: iseed = 123
 integer :: krmin = 2
 integer :: krmax = 2
 logical :: botop = .true.

! --- partion cloud liquid water into precip & detrained water --- 
 real :: rn_ptop = 500.0E2
 real :: rn_pbot = 800.0E2
 real :: rn_frac_bot = 0.8
 real :: rn_frac_top = 1.0

! --- evaporation of convective scale precip ---
 logical :: evap_on = .true.   
 real    :: cfrac   = 0.05
 real    :: hcevap  = 0.80    


    NAMELIST / ras_nml /                         &
      fracs,   rasal0,  puplim, aratio, cufric,  &
      ncrnd,   iseed,   krmin,   krmax, botop,   &
      rn_ptop, rn_pbot, rn_frac_top, rn_frac_bot,&
      evap_on, cfrac,   hcevap

!---------------------------------------------------------------------

 contains

!#####################################################################
!#####################################################################

  SUBROUTINE RAS( temp0, qvap0, uwnd0, vwnd0, pres0, pres0_int, coldT0, &
                  dtime, dtemp0, dqvap0, duwnd0, dvwnd0, rain0, snow0, &
                  kbot, ql0,  qi0,    qa0,   mc0,     Dl0,    Di0,   Da0)

!=======================================================================
! ***** DRIVER FOR RAS
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!     dtime     - Size of time step in seconds
!     pres0     - Pressure
!     pres0_int - Pressure at layer interface
!     temp0     - Temperature
!     qvap0     - Water vapor 
!     uwnd0     - U component of wind
!     vwnd0     - V component of wind
!     coldT0    - should the precipitation assume to be frozen?
!     kbot      - OPTIONAL;lowest model level index (integer)
!     ql0       - OPTIONAL;cloud liquid
!     qi0       - OPTIONAL;cloud ice
!     qa0       - OPTIONAL;cloud/saturated volume fraction
!---------------------------------------------------------------------
  real, intent(in)                   :: dtime
  real, intent(in), dimension(:,:,:) :: pres0, pres0_int
  real, intent(in), dimension(:,:,:) :: temp0, qvap0, uwnd0, vwnd0
  logical, intent(in), dimension(:,:):: coldT0
  integer, intent(in), OPTIONAL, dimension(:,:)   :: kbot
  real,    intent(in), OPTIONAL, dimension(:,:,:) :: ql0,qi0,qa0
!---------------------------------------------------------------------
! Arguments (Intent out)
!       rain0  - surface rain
!       snow0  - surface snow
!       dtemp0 - Temperature change 
!       dqvap0 - Water vapor change 
!       duwnd0 - U wind      change 
!       dvwnd0 - V wind      change 
!       mc0    - OPTIONAL; cumulus mass flux
!       Dl0    - OPTIONAL; cloud liquid change
!       Di0    - OPTIONAL; cloud ice change
!       Da0    - OPTIONAL; cloud fraction change
!---------------------------------------------------------------------
  real, intent(out), dimension(:,:,:) :: dtemp0, dqvap0, duwnd0, dvwnd0
  real, intent(out), dimension(:,:)   :: rain0,snow0
  
  real, intent(out), OPTIONAL, dimension(:,:,:) :: mc0,Dl0,Di0,Da0

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

 real, parameter :: p00 = 1000.0E2

 real, allocatable, dimension(:,:) :: t_parc, q_parc, p_parc 

 integer, allocatable, dimension(:) :: idex, jdex,   ic,   ksfc
 real,    allocatable, dimension(:) :: dpcu, precip, psfc, hl
 logical, allocatable, dimension(:) :: coldT    

 real, allocatable, dimension(:,:) ::                         &
   temp,   qvap,    uwnd,   vwnd,      pres,      pres_int,   &
   pi,     pi_int,  theta,  cp_by_dp,  qvap_sat,  dqvap_sat,  &
   alpha,  beta,    gamma,  mass,      dtcu,      dqcu,       &
   ducu,   dvcu,    dtemp,  dqvap,     duwnd,     dvwnd,      &
   mccu,   Dlcu,    Dicu,   Dacu,      mc,        ql,         &
   qi,     qa,      Dl,     Di,        Da
       
!---------------------------------------------------------------------

 logical, dimension(size(temp0,1),size(temp0,2)) :: lcl_mask
 integer, dimension(size(temp0,1),size(temp0,2)) :: kcbase
 real,    dimension(size(temp0,1),size(temp0,2)) :: psfc0    

 logical :: setras
 integer :: i, imax, j, jmax, k, kmax
 integer :: ncmax, nc, ib
 integer :: nn, klcl, nlcl, kfull, klcl_min, klcl_max
 real    :: rasal, frac

!=====================================================================

! --- Check to see if ras has been initialized
  if( do_init ) CALL ERROR_MESG( 'RAS',  &
                                 'ras_init has not been called', FATAL )
  
! --- Set dimensions
  imax  = size( temp0, 1 )
  jmax  = size( temp0, 2 )
  kmax  = size( temp0, 3 )
  kfull = kmax / 2

! --- Initalize
   dtemp0 = 0.0                               
   dqvap0 = 0.0                               
   duwnd0 = 0.0                               
   dvwnd0 = 0.0                                                                        
    rain0 = 0.0   
    snow0 = 0.0

    frac  = fracs  / dtime
    rasal = rasal0 / dtime  

  if (PRESENT(Da0)) then
      Da0 = 0.0
      Dl0 = 0.0
      Di0 = 0.0
  end if
  if (PRESENT(mc0)) then
      mc0 = 0.0
  end if
       
! --- Find LCL ---> cloud base
  allocate (   t_parc(imax,jmax) )
  allocate (   q_parc(imax,jmax) )
  allocate (   p_parc(imax,jmax) )
  allocate ( qvap_sat(imax,jmax) )

  if( PRESENT( kbot ) ) then
     do j = 1,jmax
     do i = 1,imax
        k = kbot(i,j)
          t_parc(i,j) = temp0(i,j,k)
          q_parc(i,j) = qvap0(i,j,k)
          p_parc(i,j) = pres0(i,j,k)
     end do
     end do
  else
          t_parc(:,:) = temp0(:,:,kmax)
          q_parc(:,:) = qvap0(:,:,kmax)
          p_parc(:,:) = pres0(:,:,kmax)
  end if

  CALL ESCOMP( t_parc, qvap_sat )

  qvap_sat = d622 * qvap_sat / MAX(qvap_sat,(p_parc-d378*qvap_sat))

  q_parc = MAX( q_parc, 1.0E-6   )
  q_parc = MIN( q_parc, qvap_sat )

  CALL COMP_LCL( t_parc, q_parc, p_parc, pres0, kcbase )

  deallocate (   t_parc )
  deallocate (   q_parc )
  deallocate (   p_parc )
  deallocate ( qvap_sat )

  klcl_min = MINVAL( kcbase )
  klcl_max = MAXVAL( kcbase )

! --- Set surface pressure
  if( PRESENT( kbot ) ) then
     do j = 1,jmax
     do i = 1,imax
        k = kbot(i,j) + 1
        psfc0(i,j) = pres0_int(i,j,k)
     end do
     end do
  else
        psfc0(:,:) = pres0_int(:,:,kmax+1)
  end if

!---------------------------------------------------------------------
! Cloud base loop starts
!---------------------------------------------------------------------
  do klcl = klcl_min, klcl_max
!---------------------------------------------------------------------

! --- Identify points with cloud base at klcl
  lcl_mask = klcl .eq. kcbase
  nlcl     = COUNT( lcl_mask ) 

! --- If there are no points - go to end of loop
  if ( nlcl .eq. 0 ) go to 9000

! --- Allocate storage for packed input/output variables 
  allocate ( temp    (nlcl,kmax  ) )
  allocate ( qvap    (nlcl,kmax  ) )
  allocate ( uwnd    (nlcl,kmax  ) )
  allocate ( vwnd    (nlcl,kmax  ) )
  allocate ( pres    (nlcl,kmax  ) )
  allocate ( pres_int(nlcl,kmax+1) )
  allocate ( dtemp   (nlcl,kmax)   )
  allocate ( dqvap   (nlcl,kmax)   )
  allocate ( duwnd   (nlcl,kmax)   )
  allocate ( dvwnd   (nlcl,kmax)   )
  allocate ( precip  (nlcl)        )

  if (PRESENT(Da0)) then
  allocate ( qa      (nlcl,kmax  ) )
  allocate ( ql      (nlcl,kmax  ) )
  allocate ( qi      (nlcl,kmax  ) )
  allocate ( Da      (nlcl,kmax  ) )
  allocate ( Dl      (nlcl,kmax  ) )
  allocate ( Di      (nlcl,kmax  ) )
  end if
  if (PRESENT(mc0)) then
  allocate ( mc      (nlcl,kmax+1) )
  end if

! --- Allocate storage for index arrays
  allocate ( idex (nlcl) )
  allocate ( jdex (nlcl) )
  allocate ( ksfc (nlcl) )
  allocate ( ic   (kmax) )

! --- Allocate storage for local arrays
  allocate ( pi       (nlcl,kmax)   )
  allocate ( pi_int   (nlcl,kmax+1) )
  allocate ( theta    (nlcl,kmax)   )
  allocate ( cp_by_dp (nlcl,kmax)   )
  allocate ( qvap_sat (nlcl,kmax)   )
  allocate ( dqvap_sat(nlcl,kmax)   )
  allocate ( alpha    (nlcl,kmax)   )
  allocate ( beta     (nlcl,kmax)   )
  allocate ( gamma    (nlcl,kmax)   )
  allocate ( dtcu     (nlcl,kmax)   )
  allocate ( dqcu     (nlcl,kmax)   )
  allocate ( ducu     (nlcl,kmax)   )
  allocate ( dvcu     (nlcl,kmax)   )
  allocate ( dpcu     (nlcl)        )
  allocate ( psfc     (nlcl)        )
  allocate ( mass     (nlcl,kmax)   )
  allocate ( hl       (nlcl)        )
  allocate ( coldT    (nlcl)        )
  
  if (PRESENT(Da0)) then
  allocate ( Dacu     (nlcl,kmax  ) )
  allocate ( Dlcu     (nlcl,kmax  ) )
  allocate ( Dicu     (nlcl,kmax  ) )
  end if
  if (PRESENT(mc0)) then
  allocate ( mccu     (nlcl,kmax+1) )
  end if

! --- Set index arrays
    nn = 0
  do j = 1,jmax
  do i = 1,imax
  if ( lcl_mask(i,j) ) then
    nn = nn + 1
    idex(nn) = i
    jdex(nn) = j
  end if
  end do
  end do

! --- Pack input variables
  do i = 1,nn
       psfc(i) = psfc0(idex(i),jdex(i))
  end do
  do i = 1,nn
  do k = 1,kmax
       temp(i,k) = temp0(idex(i),jdex(i),k)
       qvap(i,k) = qvap0(idex(i),jdex(i),k)
       uwnd(i,k) = uwnd0(idex(i),jdex(i),k)
       vwnd(i,k) = vwnd0(idex(i),jdex(i),k)
       pres(i,k) = pres0(idex(i),jdex(i),k)
  end do
  end do

  do i = 1,nn
       coldT(i)=coldT0(idex(i),jdex(i))
       if(coldT(i)) then
          hl(i)=HLs
       else  
          hl(i)=HLv
       end if
  end do
 
  if( PRESENT( kbot ) ) then
     do i = 1,nn
       ksfc(i) = kbot(idex(i),jdex(i)) + 1
     end do
  else
       ksfc(1:nn) = kmax + 1
  endif

  if (PRESENT(Da0)) then
  do i = 1,nn
  do k = 1,kmax
       qa(i,k)   = qa0(idex(i),jdex(i),k)
       ql(i,k)   = ql0(idex(i),jdex(i),k)
       qi(i,k)   = qi0(idex(i),jdex(i),k)
  end do
  end do
  end if

  do i = 1,nn
  do k = 1,kmax+1
       pres_int(i,k) = pres0_int(idex(i),jdex(i),k)
  end do
  end do

! --- Initalize
   dtemp = 0.0                               
   dqvap = 0.0                               
   duwnd = 0.0                               
   dvwnd = 0.0                                                                        
  precip = 0.0     

  if (PRESENT(Da0)) then
      Da  = 0.0
      Dl  = 0.0
      Di  = 0.0
  end if
  if (PRESENT(mc0)) then
      mc  = 0.0
  end if
   
! --- Thickness of layers 
  mass(:,1:kmax) = pres_int(:,2:kmax+1) - pres_int(:,1:kmax) 
  mass(:,1:kmax) = MAX( mass(:,1:kmax), 1.0e-5 )

! --- Compute exner functions 
! --- at layer interfaces
  pi_int = ( pres_int / p00  ) ** Kappa                                   
! --- at full levels
  pi(:,1:kmax) = ( pi_int(:,2:kmax+1) * pres_int(:,2:kmax+1)                 &
                 - pi_int(:,1:kmax  ) * pres_int(:,1:kmax  ) )               &
                 / ( mass(:,1:kmax  ) * one_plus_kappa )
  pi(:,1:kmax) = MAX( pi(:,1:kmax), 1.0e-5 )

! --- Compute potential temperature
  theta = temp / pi                                 

! --- Compute Cp divided by dpres                  
  cp_by_dp(:,1:kmax) = Cp / mass(:,1:kmax)

! --- Compute mass of layers              
  mass = mass / Grav

! --- Set order in which clouds are to be done                   
  CALL RAS_CLOUD_INDEX ( ic, ncmax )

  setras  = .true.

!---------------------------------------------------------------------
! Cloud top loop starts
!---------------------------------------------------------------------
  do nc = 1,ncmax
!---------------------------------------------------------------------

 ib = ic(nc)
 if (ib .ge. klcl) CYCLE

 if ( setras ) then
! --- Compute saturation value of water vapor & its derivative
  CALL  ESCOMP( temp,  qvap_sat )
  CALL DESCOMP( temp, dqvap_sat )
     dqvap_sat = d622 * pres * dqvap_sat / &
                    ((MAX(qvap_sat,pres-d378*qvap_sat))**2.)
     qvap_sat = d622 *  qvap_sat / MAX(qvap_sat,pres-d378*qvap_sat) 
     
! --- Compute some other stuff
     alpha =  qvap_sat - dqvap_sat * temp
     beta  = dqvap_sat * pi
     do k=1,kmax
        gamma(:,k) = 1.0 /((1.0+(hl(:)*dqvap_sat(:,k)/ Cp))*pi(:,k))
     enddo
 endif

! --- Do adjustment
  if (PRESENT(Da0)) then
  CALL RAS_CLOUD(                                                    &
       klcl, ib, rasal, frac, hl, coldT,                             &
       theta, qvap, uwnd, vwnd, pres_int, pi_int, pi, psfc,          &
       alpha, beta, gamma, cp_by_dp,                                 &
       dtcu, dqcu, ducu, dvcu, dpcu, &
       ql, qi, qa, mccu, Dlcu, Dicu, Dacu )
  else if (PRESENT(mc0).and..not.PRESENT(Da0)) then
  CALL RAS_CLOUD(                                                    &
       klcl, ib, rasal, frac, hl, coldT,                             &
       theta, qvap, uwnd, vwnd, pres_int, pi_int, pi, psfc,          &
       alpha, beta, gamma, cp_by_dp,                                 &
       dtcu, dqcu, ducu, dvcu, dpcu,mccu=mccu)
  else
  CALL RAS_CLOUD(                                                    &
       klcl, ib, rasal, frac, hl, coldT,                             &
       theta, qvap, uwnd, vwnd, pres_int, pi_int, pi, psfc,          &
       alpha, beta, gamma, cp_by_dp,                                 &
       dtcu, dqcu, ducu, dvcu, dpcu)
  end if

! --- Multiply tendencies by size of time step
  dtcu =  dtcu * dtime
  dqcu =  dqcu * dtime
  ducu =  ducu * dtime
  dvcu =  dvcu * dtime
  dpcu =  dpcu * dtime

  if (PRESENT(Da0)) then
  Dacu =  Dacu * dtime
  Dlcu =  Dlcu * dtime
  Dicu =  Dicu * dtime
  end if

! --- Evaporate some precip
  if( evap_on ) then
  CALL RAS_CEVAP ( ib, temp, qvap, pres, mass, pi, qvap_sat, &
                   dqvap_sat, psfc, hl, dtime, ksfc, dtcu, dqcu, dpcu )
  endif

! --- Update potential temperature, water vapor, winds
  theta = theta + dtcu
  qvap  = qvap  + dqcu
  uwnd  = uwnd  + ducu
  vwnd  = vwnd  + dvcu

  if (PRESENT(Da0)) then
    ql = ql + Dlcu
    qi = qi + Dicu
    qa = qa + Dacu
  end if

! --- Recover temperature 
  temp = theta * pi

! --- Accumulate precip 
  precip = precip + dpcu 

! --- Accumulate tendencies 
  dtemp = dtemp + dtcu * pi
  dqvap = dqvap + dqcu
  duwnd = duwnd + ducu 
  dvwnd = dvwnd + dvcu

  if (PRESENT(Da0)) then
  Da = Da + Dacu
  Dl = Dl + Dlcu
  Di = Di + Dicu
  end if
  if (PRESENT(mc0)) then
  mc = mc + mccu
  end if

  setras = .false.

!---------------------------------------------------------------------
  end do
!---------------------------------------------------------------------
! Cloud top loop ends
!---------------------------------------------------------------------

! --- Unpack variables
  do i = 1,nn
    do k = 1,kmax
      dtemp0(idex(i),jdex(i),k) = dtemp(i,k)
      dqvap0(idex(i),jdex(i),k) = dqvap(i,k) 
      duwnd0(idex(i),jdex(i),k) = duwnd(i,k)
      dvwnd0(idex(i),jdex(i),k) = dvwnd(i,k)
    end do
    if (coldT(i)) then
      snow0(idex(i),jdex(i)) =  precip(i)
    else
      rain0(idex(i),jdex(i)) =  precip(i)
    end if
  end do

  if (PRESENT(Da0)) then
  do i = 1,nn
    do k = 1,kmax
      Da0(idex(i),jdex(i),k) = Da(i,k) 
      Dl0(idex(i),jdex(i),k) = Dl(i,k) 
      Di0(idex(i),jdex(i),k) = Di(i,k)
    end do
  end do
  end if
  if (PRESENT(mc0)) then
  do i = 1,nn
    do k = 1,kmax
      mc0(idex(i),jdex(i),k) = mc(i,k)
    end do
      mc0(idex(i),jdex(i),kmax+1) =  mc(i,kmax+1)
  end do
  end if

! --- Deallocate storage  
  deallocate ( temp     )
  deallocate ( qvap     )
  deallocate ( uwnd     )
  deallocate ( vwnd     )
  deallocate ( pres     )
  deallocate ( pres_int )
  deallocate ( dtemp    )
  deallocate ( dqvap    )
  deallocate ( duwnd    )
  deallocate ( dvwnd    )
  deallocate ( precip   )
  deallocate ( idex     )
  deallocate ( jdex     )
  deallocate ( ic       )
  deallocate ( ksfc     )
  deallocate ( pi       )
  deallocate ( pi_int   )
  deallocate ( theta    )
  deallocate ( cp_by_dp )
  deallocate ( qvap_sat )
  deallocate ( dqvap_sat)
  deallocate ( alpha    )
  deallocate ( beta     )
  deallocate ( gamma    )
  deallocate ( dtcu     )
  deallocate ( dqcu     )
  deallocate ( ducu     )
  deallocate ( dvcu     )
  deallocate ( dpcu     )
  deallocate ( psfc     )
  deallocate ( mass     )
  deallocate ( hl       )
  deallocate ( coldT    )

  if (PRESENT(Da0)) then
  deallocate ( qa       )
  deallocate ( qi       )
  deallocate ( ql       )
  deallocate ( Da       )
  deallocate ( Dl       )
  deallocate ( Di       )
  deallocate ( Dacu     )
  deallocate ( Dlcu     )
  deallocate ( Dicu     )
  end if
  if (PRESENT(mc0)) then
  deallocate ( mc       )
  deallocate ( mccu     )
  end if

!---------------------------------------------------------------------
  9000 continue
  end do
!---------------------------------------------------------------------
! Cloud base loop ends
!---------------------------------------------------------------------

!=====================================================================
  end SUBROUTINE RAS


!#####################################################################
!#####################################################################

  SUBROUTINE RAS_INIT()

!=======================================================================
! ***** INITIALIZE RAS
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
! none
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

 integer             :: unit, io
 real                :: actp, facm 
 real, dimension(15) :: au,   tem

!=====================================================================

!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------

  if( FILE_EXIST( 'input.nml' ) ) then
      unit = OPEN_FILE ( file = 'input.nml', action = 'read' )
      io = 1
  do while ( io .ne. 0 )
      READ( unit, nml = ras_nml, iostat = io, end = 10 )
  end do
  10  CALL CLOSE_FILE ( unit )
  end if

!---------------------------------------------------------------------
! --- Write namelist
!---------------------------------------------------------------------

  unit = OPEN_FILE ( file = 'logfile.out', action = 'APPEND')
  if ( get_my_pe() == 0 ) then
       WRITE( unit, '(/,80("="),/(a))') trim(version),trim(tag)
       WRITE( unit, nml = ras_nml ) 
  endif
  CALL CLOSE_FILE ( unit )

!---------------------------------------------------------------------
! --- Initialize constants
!---------------------------------------------------------------------

  
  cp_div_grav = Cp / Grav

  one_plus_kappa  = 1.0 + Kappa
  one_minus_kappa = 1.0 - Kappa

  onebcp = 1.0 / Cp
  onebg  = 1.0 / Grav
  
  rn_pfac = ( rn_frac_top -  rn_frac_bot)  / ( rn_pbot - rn_ptop ) 
  
! --- Climatological cloud work function

 actp  = 1.7
 facm  = 0.01                               
 actop = actp * facm

  a  = facm * a
  ph = 100.0*ph
   
  ac = 0.0     
  ad = 0.0     
  au = 0.0     
 tem = 0.0     

 tem(2:15) = ph(2:15) -  ph(1:14)
  au(2:15) =  a(1:14) / tem(2:15) 
  ad(2:15) =  a(2:15) / tem(2:15)

  ac(2:15) = ph(2:15) * au(2:15) - ph(1:14) * ad(2:15)
  ad(2:15) = ad(2:15) - au(2:15)

!-------------------------------------------------------------------

  do_init = .false.

!=====================================================================
  end SUBROUTINE RAS_INIT


!#####################################################################
!#####################################################################

  SUBROUTINE COMP_LCL( t_parc, q_parc, p_parc, pres, k_lcl )

!=======================================================================
! ***** COMPUTE LCL ( CLOUD BASE )
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!       t_parc   Initial parcel temperature
!       q_parc   Initial parcel mixing ratio
!       p_parc   Initial parcel pressure
!       pres     Pressure in colunm
!---------------------------------------------------------------------
  real, intent(in), dimension(:,:)   :: t_parc, q_parc, p_parc
  real, intent(in), dimension(:,:,:) :: pres

!---------------------------------------------------------------------
! Arguments (Intent out)
!       k_lcl   Index of LCL in column
!---------------------------------------------------------------------
  integer, intent(out), dimension(:,:) :: k_lcl

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

  real, dimension(size(t_parc,1),size(t_parc,2)) :: &
        esat, qsat, rhum, chi, p_lcl

  integer :: k, kmax, k_lcl_min

!=====================================================================

! --- Index of lowest level
  kmax = size( pres, 3 )

! --- Compute relative humidity
  CALL  ESCOMP ( t_parc, esat )
  qsat(:,:) = d622 * esat(:,:) / p_parc(:,:) 
  rhum(:,:) = q_parc(:,:) / qsat(:,:)
  rhum(:,:) = MIN( rhum(:,:), 1.0 )

! --- Compute exponent
   chi(:,:) = t_parc(:,:) / ( 1669.0 - 122.0*rhum(:,:) - t_parc(:,:) )

! --- Compute pressure at LCL
    rhum(:,:) =    chi(:,:) *  LOG( rhum(:,:) )
   p_lcl(:,:) = p_parc(:,:) *  EXP( rhum(:,:) )

! --- Bound p_lcl 
  p_lcl(:,:) = MAX( p_lcl(:,:), pres(:,:,1) )
  p_lcl(:,:) = MIN( p_lcl(:,:), p_parc(:,:) )

! --- Find index of LCL
  do k = 2,kmax
  where ( ( p_lcl(:,:) >= pres(:,:,k-1) ) .and. &
          ( p_lcl(:,:) <= pres(:,:,k) ) )
     k_lcl(:,:) = k
  end where
  end do

! --- Bound k_lcl
  k_lcl_min  = kmax / 2
  k_lcl(:,:) = MAX( k_lcl(:,:), k_lcl_min )

!=====================================================================
  end SUBROUTINE COMP_LCL

!#####################################################################
!#####################################################################

 SUBROUTINE RAS_ACRITN ( len, pl, pl_base, acrit )

!=======================================================================
! ***** GET VALUE OF THE CLIMATOLOGICAL CLOUD WORKFUNCTION          
!=======================================================================
!---------------------------------------------------------------------
!  Arguments (Intent in)
!      len     - size of arrays
!      pl      - Pressure
!      pl_base - Pressure at cloud base
!---------------------------------------------------------------------
  integer, intent(in)                 :: len
  real   , intent(in), dimension(len) :: pl, pl_base

!---------------------------------------------------------------------
!  Arguments (Intent out)
!      acrit - Climatological cloud workfunction
!---------------------------------------------------------------------
  real, intent(out), dimension(len) :: acrit 

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

 integer, dimension(len) :: iwk
 integer                 :: i

!=======================================================================

 do i = 1,len
   iwk(i) = pl(i) * 0.02E-2 - 0.999999999
 end do

 do i = 1,len
 if ( ( iwk(i) > 1 ) .and. ( iwk(i) <= 15 ) ) then
         acrit(i) = ac(iwk(i)) + pl(i) * ad(iwk(i))
 else if ( iwk(i) <= 1 ) then
         acrit(i) = actop
 else if ( iwk(i) > 15 ) then
         acrit(i) = a(15)
 else
         CALL ERROR_MESG ('RAS_ACRITN in RAS_MOD', 'STOP?', FATAL)
!!!!     print *, 'ERROR IN ras_acritn - STOP'
!!!!     STOP
 end if
 end do

 do i = 1,len
    acrit(i) = aratio * acrit(i) * ( pl_base(i) - pl(i) )
 end do

!=======================================================================
 end SUBROUTINE RAS_ACRITN

!#####################################################################
!#####################################################################

  SUBROUTINE RAS_BDGT ( precip, coldT, dtemp, dqvap, duwnd, dvwnd, &
                        pres_int, dql, dqi)

!=====================================================================
! ***** BUDGET CHECK FOR RAS - A DEBUGGING TOOL
!=====================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!       precip   - either rain or snow
!       coldT    - is the precip snow?
!       dtemp    - Temperature change 
!       dqvap    - Water vapor change 
!       duwnd    - U wind change 
!       dvwnd    - V wind change 
!       pres_int - Pressure at layer interface
!       dql      - OPTIONAL, liquid change
!       dqi      - OPTIONAL, ice change
!---------------------------------------------------------------------
  real, intent(in), dimension(:,:,:) :: dtemp, dqvap, duwnd, dvwnd, pres_int
  real, intent(in), dimension(:,:)   :: precip
  logical, intent(in), dimension(:,:):: coldT
  real, intent(in), optional, dimension(:,:,:) :: dql, dqi
!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

 integer :: imax, jmax, kmax, i, j, k
 real    :: sum_dtemp, sum_dqvap, sum_duwnd, sum_dvwnd
 real    :: dqvap_prec, dqvap_dtemp
 real , dimension(size(dtemp,1),size(dtemp,2),size(dtemp,3)) :: mass
 
!=====================================================================

  imax = size ( dtemp, 1 )
  jmax = size ( dtemp, 2 )
  kmax = size ( dtemp, 3 )

  mass(:,:,1:kmax) = pres_int(:,:,2:kmax+1) - pres_int(:,:,1:kmax)
  mass = mass / Grav

  do j = 1,jmax
  do i = 1,imax

    sum_dtemp = 0.                                                          
    sum_dqvap = 0. 
    sum_duwnd = 0. 
    sum_dvwnd = 0. 

  do k = 1,kmax
    sum_dtemp = sum_dtemp + dtemp(i,j,k)*mass(i,j,k)                                   
    sum_dqvap = sum_dqvap + dqvap(i,j,k)*mass(i,j,k)  
    sum_duwnd = sum_duwnd + duwnd(i,j,k)*mass(i,j,k)  
    sum_dvwnd = sum_dvwnd + dvwnd(i,j,k)*mass(i,j,k)
    if (present(dql)) then
    sum_dqvap = sum_dqvap + (dql(i,j,k)+dqi(i,j,k))*mass(i,j,k)      
    end if  
  end do

    dqvap_prec  = sum_dqvap + precip(i,j)
     
     if (coldT(i,j)) then
     dqvap_dtemp = sum_dqvap + Cp*sum_dtemp/HLs
     else
     dqvap_dtemp = sum_dqvap + Cp*sum_dtemp/HLv
     end if

   if ( ( abs( dqvap_prec  ) > 1.0E-4 ) .or.     &
        ( abs( dqvap_dtemp ) > 1.0E-4 ) .or.     &
        ( abs( sum_duwnd   ) > 1.0E-1 ) .or.     &      
        ( abs( sum_dvwnd   ) > 1.0E-1 ) )  then
    print *
    print *, ' RAS BUDGET CHECK AT i,j = ', i,j
    print *, ' dqvap + (Cp/hl)*dtemp = ',  dqvap_dtemp                                                                     
    print *, ' dqvap + precip        = ',  dqvap_prec                                                                     
    print *, ' duwnd                 = ',  sum_duwnd                                                                    
    print *, ' dvwnd                 = ',  sum_dvwnd                                                                     
    print *, 'STOP'
    STOP
   endif

  end do
  end do

!=====================================================================
  end SUBROUTINE RAS_BDGT

!#######################################################################
!#######################################################################

 SUBROUTINE RAS_CEVAP ( type, temp, qvap, pres, mass, pi, qvap_sat, &
                        dqvap_sat,psfc,hl,dtime,ksfc,dtcu, dqcu, dpcu )

!=======================================================================
! EVAPORATION OF CONVECTIVE SCALE PRECIP         
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!     type  - cloud type index
!     temp  - Temperature
!     qvap  - Water vapor
!     pres  - Pressure
!     mass  - Mass of layers
!     pi    - Exner function
!     qvap_sat - Saturation value of qvap
!     dqvap_sat - Temperature derivative of qvap_sat
!     psfc  - surface presure
!     hl    - proper latent heat for the column
!     dtime - size of time step in seconds
!     ksfc  - index of lowest model level
!---------------------------------------------------------------------

  integer :: type
  real    :: dtime

  real,    intent(in), dimension(:,:) :: temp, qvap, pres, mass, pi
  real, intent(in), dimension(:,:)    :: qvap_sat, dqvap_sat
  real,    intent(in), dimension(:)   :: psfc,hl
  integer, intent(in), dimension(:)   :: ksfc

!---------------------------------------------------------------------
! Arguments (Intent inout)
!     dtcu - Temperature change due to convection
!     dqcu - Water vapor change due to convection
!     dpcu - Precip in mm
!---------------------------------------------------------------------
  real, intent(inout), dimension(:,:) :: dtcu, dqcu
  real, intent(inout), dimension(:)   :: dpcu

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

  real, parameter :: cem  = 0.054
  real, parameter :: ceta = -544.0E-6

  integer :: itopp1, kmax, k

  real, dimension(size(temp,1),size(temp,2)) ::  temp_new, qvap_new

  real, dimension(size(temp,1)) :: pevap, prec, def, evef
  real, dimension(size(temp,1)) :: prec_mmph, pfac, emx

!=======================================================================

  kmax   = size(temp,2)
  itopp1 = type + 1

! --- Initalize
  pevap    = 0.0
  qvap_new = qvap
  temp_new = temp

!=======================================================================

  do k = itopp1,kmax
!-----------------------------------------

! --- Compute whats available for evaporation 
   prec = MAX(dpcu - pevap, 0.0 )

! --- Compute precipitation efficiency factor
  prec_mmph = prec * 3600.0 / dtime
  pfac      = SQRT( pres(:,k) / psfc )
  emx       = SQRT( cem * cfrac * prec_mmph * pfac )   
  evef      = 1.0 - EXP( ceta * dtime * emx ) 

! --- Evaporate precip where needed
  where ( ( hcevap*qvap_sat(:,k) >= qvap(:,k) ) .and.  &
          ( prec(:)              > 0.0        ) .and.  &
          ( ksfc(:)              > k          ) )
            def = ( hcevap*qvap_sat(:,k) - qvap(:,k) ) /    &
                  ( 1.0 + (hl(:)* hcevap* dqvap_sat(:,k)/Cp) )
            def = evef*def
            def = MIN( def, prec/mass(:,k) )
  qvap_new(:,k) = qvap(:,k) + def
  temp_new(:,k) = temp(:,k) - (def * hl(:)/Cp)
          pevap = pevap + def * mass(:,k)
  end where

!-----------------------------------------
  end do

! --- Add on changes to temperature and water vapor from evaporation
  dtcu = dtcu + ( temp_new - temp ) / pi
  dqcu = dqcu + ( qvap_new - qvap )

! --- Update precip 
  dpcu = MAX(dpcu - pevap,0.)

!=======================================================================
  end SUBROUTINE RAS_CEVAP

!#######################################################################
!#######################################################################

 SUBROUTINE RAS_CLOUD(                                           &
            k, ic, rasal, frac, hl, coldT,                       &
            theta, qvap, uwnd, vwnd, pres_int, pi_int, pi, psfc, &
            alf, bet, gam, cp_by_dp,                             &
            dtcu, dqcu, ducu,  dvcu, dpcu,                       &
            ql, qi, qa, mccu, Dlcu, Dicu, Dacu )
!=======================================================================
! RAS Cu Parameterization 
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!     k        : Cloud base index.
!     ic       : Cloud top index.
!     hl       : proper latent heat for the column
!     coldT    : is the precipitation frozen?
!     theta    : Potential temperature
!     qvap     : Water vapor 
!     uwnd     : U component of wind
!     vwnd     : V component of wind
!     pres_int : Pressure       at layer interface
!     pi_int   : Exner function at layer interface
!     pi       : Exner function 
!     psfc     : Surface pressure
!     rasal    : Relaxation parameter for cloud type ic
!     ql       : OPTIONAL, cloud liquid
!     qi       : OPTIONAL, cloud ice
!     qa       : OPTIONAL, cloud fraction or saturated volume fraction
!---------------------------------------------------------------------

 real,    intent(in) :: rasal, frac
 integer, intent(in) :: ic, k
 real, intent(in), dimension(:)   :: hl, psfc
 logical, intent(in), dimension(:):: coldT
 real, intent(in), dimension(:,:) :: theta, qvap, uwnd, vwnd, pres_int, pi_int
 real, intent(in), dimension(:,:) :: alf, bet, gam, pi, cp_by_dp
 real, intent(in), OPTIONAL, dimension(:,:) :: ql,qi,qa
!---------------------------------------------------------------------
! Arguments (Intent out)
!     dpcu    : Precip for cloud type ic.
!     dtcu    : Potential temperature change
!     dqcu    : Water vapor change
!     ducu    : Cu friction - u component.
!     dvcu    : Cu friction - v component.
!     mccu    : OPTIONAL ; Cumulus Mass Flux
!     Dacu    : OPTIONAL ; Detrained saturated mass fraction tendency
!     Dlcu    : OPTIONAL ; Detrained cloud liquid tendency
!     Dicu    : OPTIONAL ; Detrained cloud ice tendency
!---------------------------------------------------------------------

 real, intent(out), dimension(:)   :: dpcu
 real, intent(out), dimension(:,:) :: dtcu, dqcu, ducu, dvcu
 real, intent(out), OPTIONAL, dimension(:,:) :: mccu, Dacu, Dlcu, Dicu

!---------------------------------------------------------------------
!    (Intent local)
!---------------------------------------------------------------------

 real, parameter :: rhmax  = 0.9999 

 integer :: len

 real,    dimension(size(theta,1)) :: wfn, akm, qs1, uht, vht, wlq, alm, pcu 
 real,    dimension(size(theta,1)) :: tx1, tx2, tx3, tx4, tx5, tx6, tx7, tx8
 integer, dimension(size(theta,1)) :: ia, i1, i2
 real,    dimension(size(theta,1)) :: wll,wli
 
 real, dimension(size(theta,1),size(theta,2)) :: gmh, eta, hol, hst, qol

 integer :: km1, ic1, len1, len11, len2, lena, lenb, lena1
 integer :: i, l, ii
 real    :: tem1, tem, tem2, rasalf

!=====================================================================

! Initialize
  dtcu = 0.0
  dqcu = 0.0
  ducu = 0.0
  dvcu = 0.0
  dpcu = 0.0
  if (PRESENT(Dacu)) then
  Dacu = 0.0
  Dlcu = 0.0
  Dicu = 0.0
  end if
  if (PRESENT(mccu)) then
  mccu = 0.0
  end if

! Initialize
  wfn=0.0
  akm=0.0
  qs1=0.0
  uht=0.0
  vht=0.0
  wlq=0.0
  alm=0.0
  pcu=0.0 
  tx1=0.0
  tx2=0.0
  tx3=0.0
  tx4=0.0
  tx5=0.0
  tx6=0.0
  tx7=0.0
  tx8=0.0
  ia=0
  i1=0
  i2=0
  wll=0.0
  wli=0.0
  gmh=0.0
  eta=0.0
  hol=0.0
  hst=0.0
  qol=0.0
  km1=0
  ic1=0
  len=0
  len1=0
  len11=0
  len2=0
  lena=0
  lenb=0
  lena1=0
  i=0
  l=0
  ii=0
  tem1=0.0
  tem=0.0
  tem2=0.0
  rasalf=0.0


 len = size( theta, 1 )

  km1 = k  - 1
  ic1 = ic + 1

!---------------------------------------------------------------------

! --- at cloud base
 do i = 1,len
    tx1(i)   = pi_int(i,k+1) * theta(i,k)
    qs1(i)   = alf(i,k) + bet(i,k)*theta(i,k)
    qol(i,k) = MIN( qs1(i)*rhmax, qvap(i,k) )
    hol(i,k) = tx1(i)*Cp + qol(i,k)*hl(i)
    eta(i,k) = 0.0
    tx2(i)   = (pi_int(i,k+1) - pi_int(i,k)) * theta(i,k) * Cp
 end do

! --- between cloud base & cloud top
 if ( ic < km1 ) then
 do l = km1,ic1,-1
 do i = 1,len
    qs1(i)   = alf(i,l) + bet(i,l)*theta(i,l)
    qol(i,l) = MIN( qs1(i)*rhmax, qvap(i,l) )
        tem1 = tx2(i) + pi_int(i,l+1) * theta(i,l) * Cp 
    hol(i,l) = tem1 + qol(i,l )* hl(i)
    hst(i,l) = tem1 + qs1(i)   * hl(i)
    tx1(i)   = (pi_int(i,l+1) - pi_int(i,l)) * theta(i,l)
    eta(i,l) = eta(i,l+1) + tx1(i)*cp_div_grav
    tx2(i)   = tx2(i)     + tx1(i)*Cp    
 end do
 end do
 end if

! --- at cloud top
 do i = 1,len
    hol(i,ic) = tx2(i)
    qs1(i)    = alf(i,ic) + bet(i,ic)*theta(i,ic)
    qol(i,ic) = MIN( qs1(i)*rhmax, qvap(i,ic) ) 
         tem1 = tx2(i) + pi_int(i,ic1) * theta(i,ic) * Cp 
    hol(i,ic) = tem1 + qol(i,ic) * hl(i)
    hst(i,ic) = tem1 + qs1(i)    * hl(i)
    tx3(i)    = (pi_int(i,ic1) - pi(i,ic)) * theta(i,ic)
    eta(i,ic) = eta(i,ic1) + cp_div_grav * tx3(i)
 end do

!---------------------------------------------------------------
!     ENTRAINMENT PARAMETER
!---------------------------------------------------------------

 do i = 1,len
    tx2(i) = hol(i,k)  - hst(i,ic)
    tx1(i) = 0.0
 end do

 do l = ic,km1
 do i = 1,len
    tx1(i) = tx1(i) + (hst(i,ic) - hol(i,l)) * (eta(i,l) - eta(i,l+1))
 end do
 end do

 len1 = 0
 len2 = 0

 do i = 1,len
 if ( (tx1(i) > 0.0) .and. (tx2(i) > 0.0) ) then
        len1 = len1 + 1
    ia(len1) = i
   alm(len1) = tx2(i) / tx1(i)
 end if
 end do

 len2 = len1

 if ( ic1 < k ) then
 do i = 1,len
 if ( (tx2(i) <= 0.0) .and. (hol(i,k) > hst(i,ic1))) then
         len2 = len2 + 1
     ia(len2) = i
    alm(len2) = 0.0
 end if
 end do
 end if
 
 if ( len2 .eq. 0 )  RETURN

 len11 = len1 + 1
 
!---------------------------------------------------------------
!    NORMALIZED MASSFLUX
!---------------------------------------------------------------

 do i = 1,len2
   ii = ia(i)
      eta(i,k) = 1.0
      tx2(i)   = 0.5 * (pres_int(ii,ic) + pres_int(ii,ic1))
      tx4(i)   = pres_int(ii,k)
 end do

 do i = len11,len2
   ii = ia(i)
      wfn(i) = 0.0
 if ( hst(ii,ic1) < hst(ii,ic) ) then
      tx6(i) = (hst(ii,ic1)-hol(ii,k))/(hst(ii,ic1)-hst(ii,ic))
 else
      tx6(i) = 0.0
 endif
      tx2(i) = 0.5*(pres_int(ii,ic1) + pres_int(ii,ic1+1)) * (1.0-tx6(i))   &
                              + tx2(i) * tx6(i)
 end do

  CALL RAS_ACRITN ( len2, tx2, tx4, tx3 )
 do l = km1,ic,-1
 do i = 1,len2
    tx1(i) = eta(ia(i),l)
 end do
 do i = 1,len2
    eta(i,l) = 1.0 + alm(i) * tx1(i)
 end do
 end do

!---------------------------------------------------------------
!     CLOUD WORKFUNCTION
!---------------------------------------------------------------

 if ( len1 > 0 ) then
 do i = 1,len1
   ii = ia(i)
     wfn(i) = - gam(ii,ic) * (pi_int(ii,ic1) - pi(ii,ic))     &
                           *  hst(ii,ic) * eta(i,ic1)
 end do
 end if

 do i = 1,len2
   ii = ia(i)
     tx1(i) = hol(ii,k)
 end do

 if (ic1 <= km1) then
! ----------------------
 do l = km1,ic1,-1
 do i = 1,len2
   ii = ia(i)
       tem = tx1(i) + (eta(i,l) - eta(i,l+1)) * hol(ii,l)
    pcu(i) = pi_int(ii,l+1) - pi(ii,l)
      tem1 = eta(i,l+1) * pcu(i)
    tx1(i) = tx1(i)*pcu(i)
    pcu(i) = pi(ii,l) - pi_int(ii,l)
      tem1 = (tem1 + eta(i,l) * pcu(i)) * hst(ii,l)
    tx1(i) = tx1(i) + tem*pcu(i)
    wfn(i) = wfn(i) + (tx1(i) - tem1) * gam(ii,l)
    tx1(i) = tem
 end do
 end do
! ----------------------
 endif

 lena = 0

 if (len1 > 0) then
! ----------------------
 do i = 1,len1
   ii = ia(i)
      wfn(i) = wfn(i) + tx1(i) * gam(ii,ic)*(pi_int(ii,ic1)-pi(ii,ic))  &
                      - tx3(i)
 if ( wfn(i) > 0.0 ) then
           lena = lena + 1
       i1(lena) = ia(i)
       i2(lena) = i
      tx1(lena) = wfn(i)
      tx2(lena) = qs1(ia(i))
      tx6(lena) = 1.0
 end if
 end do
! ----------------------
 end if

 lenb = lena

 do i = len11,len2
       wfn(i) = wfn(i) - tx3(i)
 if ( (wfn(i) > 0.0) .and. (tx6(i) > 0.0) ) then
         lenb = lenb + 1
     i1(lenb) = ia(i)
     i2(lenb) = i
    tx1(lenb) = wfn(i)
    tx2(lenb) = qs1(ia(i))
    tx4(lenb) = tx6(i)
 end if
 end do
 
 if ( lenb <= 0 )  RETURN

 do i = 1,lenb
    wfn(i) = tx1(i)
    qs1(i) = tx2(i)
 end do

 do l = ic,k
 do i = 1,lenb
    tx1(i) = eta(i2(i),l)
 end do
 do i = 1,lenb
    eta(i,l) = tx1(i)
 end do
 end do

  lena1 = lena + 1

 do i = 1,lena
   ii = i1(i)
     tx8(i) = hst(ii,ic) - hol(ii,ic)
 end do

 do i = lena1,lenb
   ii = i1(i)
     tx6(i) = tx4(i)
        tem = tx6(i) * (hol(ii,ic)-hol(ii,ic1)) + hol(ii,ic1)
     tx8(i) = hol(ii,k) - tem
       tem1 = tx6(i) * (qol(ii,ic)-qol(ii,ic1)) + qol(ii,ic1)
     tx5(i) = tem    - tem1 * hl(ii)
     qs1(i) = tem1   + tx8(i)*(1.0/hl(ii))
     tx3(i) = hol(ii,ic)

 end do

 do i = 1,lenb
   ii = i1(i)
     wlq(i) = qol(ii,k) - qs1(i)     * eta(i,ic)
     uht(i) = uwnd(ii,k) - uwnd(ii,ic) * eta(i,ic)
     vht(i) = vwnd(ii,k) - vwnd(ii,ic) * eta(i,ic)
     tx7(i) = hol(ii,k)
 end do

 do l = km1,ic,-1
 do i = 1,lenb
   ii = i1(i)
        tem = eta(i,l) - eta(i,l+1)
     wlq(i) = wlq(i) + tem * qol(ii,l)
     uht(i) = uht(i) + tem * uwnd(ii,l)
     vht(i) = vht(i) + tem * vwnd(ii,l)
 end do
 end do

!---------------------------------------------------------------
!         CALCULATE TRACER UPDRAFT PROPERTIES
!            AND CONVECTIVE SOURCE OF TRACER
!---------------------------------------------------------------
 
 if (PRESENT(Dacu)) then

 do i = 1,lenb
   ii = i1(i)
     wll(i) = ql(ii,k)
     wli(i) = qi(ii,k)
 end do

 do l = km1,ic,-1
 do i = 1,lenb
   ii = i1(i)
        tem = eta(i,l) - eta(i,l+1)
     wll(i) = wll(i) + tem * ql(ii,l)
     wli(i) = wli(i) + tem * qi(ii,l)
 end do
 end do
 
 !do cloud base level
 do i = 1,lenb
   ii = i1(i)
     Dlcu(ii,k) = (ql(ii,km1) - ql(ii,k))*0.5*cp_by_dp(ii,k)/Cp/onebg
     Dicu(ii,k) = (qi(ii,km1) - qi(ii,k))*0.5*cp_by_dp(ii,k)/Cp/onebg
     Dacu(ii,k) = (qa(ii,km1) - qa(ii,k))*0.5*cp_by_dp(ii,k)/Cp/onebg
     mccu(ii,k) = eta(i,k)
 end do

  if (ic1 <= km1) then
! ---------------------------------------
 do l = km1,ic1,-1
 do i = 1,lenb
   ii = i1(i)
     Dlcu(ii,l) = (eta(i,l+1)*(ql(ii,l  ) - ql(ii,l+1)) + &
                   eta(i,l  )*(ql(ii,l-1) - ql(ii,l  ))   ) * &
                  0.5*cp_by_dp(ii,l)/Cp/onebg
     Dicu(ii,l) = (eta(i,l+1)*(qi(ii,l  ) - qi(ii,l+1)) + &
                   eta(i,l  )*(qi(ii,l-1) - qi(ii,l  ))   ) * &
                  0.5*cp_by_dp(ii,l)/Cp/onebg
     Dacu(ii,l) = (eta(i,l+1)*(qa(ii,l  ) - qa(ii,l+1)) + &
                   eta(i,l  )*(qa(ii,l-1) - qa(ii,l  ))   ) * &
                  0.5*cp_by_dp(ii,l)/Cp/onebg
     mccu(ii,l) = eta(i,l)
 end do
 end do
! ---------------------------------------
  endif

 !do cloud top level
 do i = 1,lenb
   ii = i1(i)
     Dlcu(ii,ic) = (eta(i,ic1)*(ql(ii,ic) - ql(ii,ic1))* &
                   0.5*cp_by_dp(ii,ic)/Cp/onebg ) + &
                   (wll(i) - eta(i,ic)*ql(ii,ic)) &
                      *cp_by_dp(ii,ic)/Cp/onebg
     Dicu(ii,ic) = (eta(i,ic1)*(qi(ii,ic) - qi(ii,ic1))* &
                   0.5*cp_by_dp(ii,ic)/Cp/onebg ) + &
                   (wli(i) - eta(i,ic)*qi(ii,ic)) &
                      *cp_by_dp(ii,ic)/Cp/onebg
     Dacu(ii,ic) = (eta(i,ic1)*(qa(ii,ic) - qa(ii,ic1))* &
                   0.5*cp_by_dp(ii,ic)/Cp/onebg ) + &
                   (eta(i,ic) - qa(ii,ic)*eta(i,ic)) &
                      *cp_by_dp(ii,ic)/Cp/onebg
 end do

 end if


 if (PRESENT(mccu).and..not.PRESENT(Dacu)) then

 !do cloud base level
 do i = 1,lenb
   ii = i1(i)
     mccu(ii,k) = eta(i,k)
 end do

  if (ic1 <= km1) then
! ---------------------------------------
 do l = km1,ic1,-1
 do i = 1,lenb
   ii = i1(i)
   mccu(ii,l) = eta(i,l)
 end do
 end do
! ---------------------------------------
  endif

 end if

!---------------------------------------------------------------
!     CALCULATE GS AND PART OF AKM (THAT REQUIRES ETA)
!---------------------------------------------------------------

 do i = 1,lenb
   ii = i1(i)
            tem = (theta(ii,km1) - theta(ii,k)) / (pi(ii,k) - pi(ii,km1))
     hol(i,k)   = tem * (pi_int(ii,k)-pi(ii,km1))*pi(ii,k)*cp_by_dp(ii,k)
     hol(i,km1) = tem * (pi(ii,k)-pi_int(ii,k))*pi(ii,km1)*cp_by_dp(ii,km1)
     akm(i)     = 0.0
     tx2(i)     = 0.5 * (pres_int(ii,ic) + pres_int(ii,ic1))
 end do

 if ( ic1 <= km1 ) then
! ---------------------------------------
 do l = km1,ic1,-1
 do i = 1,lenb
   ii = i1(i)
            tem = ( theta(ii,l-1) - theta(ii,l) ) * eta(i,l)            &
                                 / ( pi(ii,l) - pi(ii,l-1) )
     hol(i,l)   = tem * ( pi_int(ii,l)-pi(ii,l-1) ) * pi(ii,l)          &
                      *   cp_by_dp(ii,l)  + hol(i,l)  
     hol(i,l-1) = tem * ( pi(ii,l)-pi_int(ii,l) ) * pi(ii,l-1)          &
                                                * cp_by_dp(ii,l-1)
     akm(i)     = akm(i) - hol(i,l)                                     &
                * ( eta(i,l)   * ( pi(ii,l  ) - pi_int(ii,l) ) +        &
                    eta(i,l+1) * ( pi_int(ii,l+1) - pi(ii,l) ) ) / pi(ii,l)
 end do
 end do
! ---------------------------------------
 end if


  if (PRESENT(Dacu)) then !detrain non-precipitated fraction of condensed vapor
    
    CALL RAS_RNCL( lenb, tx2, tx1)
    
    do i = 1,lenb
         ii = i1(i)
         if (coldT(ii)) then
         Dicu(ii,ic)= Dicu(ii,ic) + &
                      (1.0 - tx1(i)) * wlq(i) * cp_by_dp(ii,ic)/Cp/onebg
         else
         Dlcu(ii,ic)= Dlcu(ii,ic) + &
                      (1.0 - tx1(i)) * wlq(i) * cp_by_dp(ii,ic)/Cp/onebg
         end if
         tx2(i) = 0.
         wlq(i) = tx1(i) * wlq(i)
         tx1(i) = hol(i,ic)
    end do    

  else  !do original code

    CALL RAS_RNCL( lenb, tx2, tx1 )

    do i = 1,lenb
         tx2(i) = (1.0 - tx1(i)) * wlq(i)
         wlq(i) = tx1(i) * wlq(i)
         tx1(i) = hol(i,ic)
    end do

  end if


 do i = lena1, lenb
   ii = i1(i)
   tx1(i) = tx1(i) + (tx5(i)-tx3(i)+qol(ii,ic)*hl(ii))*(cp_by_dp(ii,ic)/Cp)
 end do

 do i = 1,lenb
    hol(i,ic) = tx1(i) - (tx2(i) * hl(i1(i)) * cp_by_dp(i1(i),ic) / Cp)
 end do

 if ( lena > 0 ) then
 do i = 1,lena
   ii = i1(i)
   akm(i) = akm(i) - eta(i,ic1) * (pi_int(ii,ic1) - pi(ii,ic))     & 
                                      * tx1(i) / pi(ii,ic)
 end do
 end if

!---------------------------------------------------------------
!     CALCULATE GH
!---------------------------------------------------------------
 
 do i = 1,lenb
   ii = i1(i)
     tx3(i)   = qol(ii,km1) - qol(ii,k)
     gmh(i,k) = hol(i,k) + (tx3(i) * cp_by_dp(ii,k) * hl(ii) * 0.5 /Cp)
     akm(i)   = akm(i) + gam(ii,km1)*(pi_int(ii,k)-pi(ii,km1)) * gmh(i,k)
 end do

  if (ic1 <= km1) then
! ---------------------------------------
 do l = km1,ic1,-1
 do i = 1,lenb
   ii = i1(i)
     tx2(i)   = tx3(i)
     tx3(i)   = (qol(ii,l-1) - qol(ii,l)) * eta(i,l)
     tx2(i)   = tx2(i) + tx3(i)
     gmh(i,l) = hol(i,l) + (tx2(i)   * cp_by_dp(ii,l) * hl(ii)*0.5 /Cp)
 end do
 end do
! ---------------------------------------
  endif

 do i = lena1,lenb
    tx3(i) = tx3(i) + &
            (2.*(tx7(i)-tx8(i)-tx5(i)-qol(i1(i),ic)*hl(i1(i)))/hl(i1(i)))
 end do

 do i = 1,lenb
     gmh(i,ic) = tx1(i) + cp_by_dp(i1(i),ic) * onebcp          &
              * (tx3(i)*hl(i1(i))*0.5 + eta(i,ic) * tx8(i))
 end do

!---------------------------------------------------------------
!     CALCULATE HC PART OF AKM
!---------------------------------------------------------------

 if ( ic1 <= km1 ) then
!---------------------------------------

 do i = 1,lenb
    tx1(i) = gmh(i,k)
 end do

 do  l = km1,ic1,-1
!---------------------------------------

 do i = 1,lenb
   ii = i1(i)
     tx1(i) = tx1(i) + (eta(i,l) - eta(i,l+1)) * gmh(i,l)
     tx2(i) = gam(ii,l-1) * (pi_int(ii,l) - pi(ii,l-1))
 end do

  if ( l .eq. ic1 ) then
  do i = lena1,lenb
     tx2(i) = 0.0
  end do
  end if

 do i = 1,lenb
   ii = i1(i)
     akm(i) = akm(i) + tx1(i) *                               &
            ( tx2(i) + gam(ii,l)*( pi(ii,l)-pi_int(ii,l) ) )
 end do

!---------------------------------------
 end do

!---------------------------------------
 end if

 do i = lena1,lenb
   ii = i1(i)
         tx2(i) = 0.5*( pres_int(ii,ic  ) + pres_int(ii,ic1) )              &
                + 0.5*( pres_int(ii,ic+2) - pres_int(ii,ic ) )*( 1.0-tx6(i) )
         tx1(i) =       pres_int(ii,ic1 )
         tx5(i) = 0.5*( pres_int(ii,ic1 ) + pres_int(ii,ic+2) )
 if ( (tx2(i) >= tx1(i)) .and. (tx2(i) < tx5(i)) ) then
         tx6(i) = 1.0 - (tx2(i) - tx1(i)) / (tx5(i) - tx1(i))
            tem = cp_by_dp(ii,ic1) / cp_by_dp(ii,ic)
     hol(i,ic1) = hol(i,ic1) + hol(i,ic) * tem
     gmh(i,ic1) = gmh(i,ic1) + gmh(i,ic) * tem
     hol(i,ic)  = 0.0
     gmh(i,ic)  = 0.0
 else if ( tx2(i) < tx1(i) ) then
     tx6(i) = 1.0
 else
     tx6(i) = 0.0
 end if
 end do

!---------------------------------------------------------------
!   MASS FLUX
!---------------------------------------------------------------
 
 do i = 1,lenb
   ii = i1(i)
 if ( ( akm(i) < 0.0 ) .and. ( wlq(i) >= 0.0 ) ) then
!jjs
  rasalf = rasal * ( pres_int(ii,ic+1) - puplim ) /      &
                   (     psfc(ii)      - puplim )
  rasalf = MAX( 0.0, rasalf )
!jjs
     wfn(i) = - tx6(i) * wfn(i) * rasalf / akm(i)
 else
     wfn(i) = 0.0
 end if
        tem = (pres_int(ii,k+1) - pres_int(ii,k))*frac
     wfn(i) = MIN( wfn(i), tem )
 end do


 if (PRESENT(Dacu)) then
     do l = ic,k
     do i = 1,lenb
          ii = i1(i)
          Dacu(ii,l) = Dacu(ii,l) * wfn(i) * onebg
          Dlcu(ii,l) = Dlcu(ii,l) * wfn(i) * onebg
          Dicu(ii,l) = Dicu(ii,l) * wfn(i) * onebg
          mccu(ii,l) = mccu(ii,l) * wfn(i) * onebg
     end do
     end do
 end if
 if (PRESENT(mccu).and..not.PRESENT(Dacu)) then
     do l = ic,k
     do i = 1,lenb
          ii = i1(i)
          mccu(ii,l) = mccu(ii,l) * wfn(i) * onebg
     end do
     end do
 end if
 
!---------------------------------------------------------------
!     PRECIPITATION
!---------------------------------------------------------------

 do i = 1,lenb
   ii = i1(i)
   dpcu(ii) =  wlq(i) * wfn(i) * onebg
 end do

!---------------------------------------------------------------
!     THETA AND Q CHANGE DUE TO CLOUD TYPE IC
!---------------------------------------------------------------
 
 do l = ic,k
 do i = 1,lenb
   ii = i1(i)
    tx4(i)    = wfn(i) * (1.0/hl(ii))
    tx5(i)    = wfn(i) * onebcp
          tem = ( gmh(i,l) - hol(i,l) ) * tx4(i)
         tem1 =  hol(i,l) * tx5(i)
   dtcu(ii,l) = tem1 / pi(ii,l)
   dqcu(ii,l) = tem
 end do
 end do

!---------------------------------------------------------------
!  CUMULUS FRICTION 
!---------------------------------------------------------------
 if (cufric) then

 do i = 1,lenb
      tx5(i) = tx5(i) * 0.5
 end do

! --- At cloud base
 do i = 1,lenb
   ii = i1(i)
          tem = tx5(i) * cp_by_dp(ii,k)
    tx1(i)    = ( uwnd(ii,km1) - uwnd(ii,k) )
    tx2(i)    = ( vwnd(ii,km1) - vwnd(ii,k) )
   ducu(ii,k) = tem * tx1(i)
   dvcu(ii,k) = tem * tx2(i)
 end do

! --- Between cloud base & cloud top
 do l = km1,ic1,-1
 do i = 1,lenb
   ii = i1(i)
         tem  = tx5(i) * cp_by_dp(ii,l)
         tem1 = tx1(i)
         tem2 = tx2(i)
    tx1(i)    = ( uwnd(ii,l-1) - uwnd(ii,l) ) * eta(i,l)
    tx2(i)    = ( vwnd(ii,l-1) - vwnd(ii,l) ) * eta(i,l)
   ducu(ii,l) = ( tx1(i) + tem1 ) * tem
   dvcu(ii,l) = ( tx2(i) + tem2 ) * tem
 end do
 end do

! --- At cloud top
  do i = 1,lenb
    ii = i1(i)
            tem =   tx5(i) * cp_by_dp(ii,ic)
    ducu(ii,ic) = ( tx1(i) + uht(i) + uht(i) ) * tem
    dvcu(ii,ic) = ( tx2(i) + vht(i) + vht(i) ) * tem
  end do

  end if

!=======================================================================
  end SUBROUTINE RAS_CLOUD

!#######################################################################
!#######################################################################

 SUBROUTINE RAS_CLOUD_INDEX ( ic, ncmax )

!=======================================================================
! Set order in which clouds are to be done                   
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent out)
!       ic      Cloud order index
!       ncmax   Max number of clouds
!---------------------------------------------------------------------
 integer, dimension(:) :: ic
 integer :: ncmax

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

 integer :: kmax
 integer :: km1, kcr, kfx, nc, irnd

!=======================================================================

 kmax = size( ic )

 km1   = kmax  - 1
 kcr   = MIN( km1, krmax )
 kfx   = km1 - kcr
 ncmax = kfx + ncrnd

! --- Sequential clouds
  if ( kfx > 0) then
     if ( botop ) then
! --- bottom to top
       do nc = 1,kfx
          ic(nc) = kmax - nc
       end do
     else
! --- top to bottom
       do nc = kfx,1,-1
          ic(nc) = kmax - nc
       end do
     endif
  endif

! --- set non-used clouds to zero
  ic((kfx+1):kmax) = 0

! $$$$$$ COMMENTED OUT FOR WORKSTATION
! --- Random clouds
!if ( ncrnd > 0 ) then
!      CALL RANSET( iseed )
!    do nc = 1,ncrnd
!      irnd      = ( RANF() - 0.0005 ) * ( kcr - krmin + 1 )
!      ic(kfx+nc) = irnd + krmin
!    end do
! endif
! $$$$$$ COMMENTED OUT FOR WORKSTATION

!=======================================================================
  end SUBROUTINE RAS_CLOUD_INDEX 

!#######################################################################
!#######################################################################

  SUBROUTINE RAS_RNCL( len, pl, rno )

!====================================================================
! PARTION CLOUD LIQUID WATER INTO PRECIP & DETRAINED WATER  
!====================================================================
!---------------------------------------------------------------------
!  Arguments (Intent in)
!      len - array size
!      pl  - Pressure at detrainment level
!---------------------------------------------------------------------
  integer, intent(in)                 :: len
  real,    intent(in), dimension(len) :: pl

!---------------------------------------------------------------------
!  Arguments (Intent out)
!      rno - fraction of liquid water converted to precip
!---------------------------------------------------------------------
  real, intent(out), dimension(len) :: rno
!=====================================================================

 rno = rn_frac_top
 
 where ( ( pl >= rn_ptop ) .and. ( pl <= rn_pbot ) )
    rno = ( rn_pbot - pl ) * rn_pfac +  rn_frac_bot
 end where

 where ( pl > rn_pbot )
    rno =   rn_frac_bot
 end where
 

!=====================================================================
  end SUBROUTINE RAS_RNCL

!#######################################################################
!#######################################################################
  end MODULE RAS_MOD
