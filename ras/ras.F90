  MODULE RAS_MOD

!=======================================================================
!  RELAXED ARAKAWA/SCHUBERT (RAS) CUMULUS PARAM SCHEME MODULE          !
!=======================================================================
!  SUBROUTINE RAS_INIT        - INITIALIZE RAS
!  SUBROUTINE RAS             - DRIVER
!  SUBROUTINE RAS_CLOUD       - RAS CU PARAMETERIZATION
!  SUBROUTINE COMP_LCL        - COMPUTE LCL (CLOUD BASE)
!  SUBROUTINE RAS_CEVAP       - EVAPORATION OF CONVECTIVE SCALE PRECIP
!  SUBROUTINE RAS_CLOUD_INDEX - SET ORDER IN WHICH CLOUDS ARE TO BE DONE 
!  SUBROUTINE RAS_CLOUD_EXIST - TEST FOR INSTABILITY IN COLUMN
!  SUBROUTINE RAS_BDGT        - BUDGET CHECK FOR RAS
!=======================================================================

 use Sat_Vapor_Pres_Mod, ONLY: ESCOMP, DESCOMP
 use      Constants_Mod, ONLY:  HLv, HLs, Cp, Grav, Kappa, rdgas, rvgas
 use   Diag_Manager_Mod, ONLY: register_diag_field, send_data
 use   Time_Manager_Mod, ONLY: time_type
 use      Utilities_Mod, ONLY: FILE_EXIST, OPEN_FILE, ERROR_MESG,  &
                                get_my_pe, CLOSE_FILE, FATAL

!---------------------------------------------------------------------
 implicit none
 private
 public  :: ras, ras_init, ras_bdgt
!---------------------------------------------------------------------

!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 character(len=128) :: version = '$Id: ras.F90,v 1.5 2001/10/25 17:48:21 fms Exp $'
 character(len=128) :: tag = '$Name: fez $'
!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 real :: cp_div_grav
 real :: one_plus_kappa, one_minus_kappa
 real :: onebcp, onebg
 real :: rn_pfac

 logical :: do_init = .true.

 real, parameter :: d622        = rdgas/rvgas
 real, parameter :: d378        = 1.0-d622
 real, parameter :: Tokioka_con = 0.05 

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
!    rh_trig - Convection takes place only if the relative humidity
!              of the lowest model level exceeds rh_trig
!    alm_min - Min value for entrainment parameter.
! Tokioka_on - If true, alm_min computed using Tokioka formulation
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

 real    :: fracs      = 0.25
 real    :: rasal0     = 0.25
 real    :: puplim     = 20.0E2
 real    :: aratio     = 1.4
 logical :: cufric     = .false.
 real    :: rh_trig    = 0.0      
 real    :: alm_min    = 0.0  
 logical :: Tokioka_on = .false.

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


    NAMELIST / ras_nml /                          &
      fracs,   rasal0,  puplim, aratio, cufric,   &
      ncrnd,   iseed,   krmin,   krmax, botop,    &
      rn_ptop, rn_pbot, rn_frac_top, rn_frac_bot, &
      evap_on, cfrac,   hcevap, rh_trig, alm_min, &
      Tokioka_on

!---------------------------------------------------------------------
! DIAGNOSTICS FIELDS 
!---------------------------------------------------------------------

integer :: id_tdt_revap,  id_qdt_revap,    id_prec_revap,  &
           id_snow_revap, id_prec_conv_3d, id_pcldb

character(len=3) :: mod_name = 'ras'

real :: missing_value = -999.

!---------------------------------------------------------------------

 contains

!#####################################################################
!#####################################################################

  SUBROUTINE RAS_INIT( axes, Time )

!=======================================================================
! ***** INITIALIZE RAS
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!---------------------------------------------------------------------

 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time

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

!---------------------------------------------------------------------
! --- initialize quantities for diagnostics output
!---------------------------------------------------------------------

   id_tdt_revap = register_diag_field ( mod_name, &
     'tdt_revap', axes(1:3), Time, &
     'Temperature tendency from RAS prec evap',      'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_revap = register_diag_field ( mod_name, &
     'qdt_revap', axes(1:3), Time, &
     'Spec humidity tendency from RAS prec evap',    'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_revap = register_diag_field ( mod_name, &
     'prec_revap', axes(1:2), Time, &
     'Evap Precip rate from RAS',                   'kg/m2/s' )

   id_snow_revap = register_diag_field ( mod_name, &
     'snow_revap', axes(1:2), Time, &
     'Evap Frozen precip rate from RAS',             'kg/m2/s' )

   id_prec_conv_3d = register_diag_field ( mod_name, &
     'prec_conv_3d', axes(1:3), Time, &
     '3D Precipitation rate from RAS',      'kg/m2/s',  &
                        missing_value=missing_value               )

   id_pcldb = register_diag_field ( mod_name, &
     'pcldb', axes(1:2), Time, &
     'Pressure at cloud base from RAS',              'hPa' )

!---------------------------------------------------------------------

  do_init = .false.

!=====================================================================
  end SUBROUTINE RAS_INIT

!#####################################################################
!#####################################################################

  SUBROUTINE RAS( is,     js,      Time,    temp0,      qvap0,     &
                  uwnd0,  vwnd0,   pres0,   pres0_int,  coldT0,    &
                  dtime,  dtemp0,  dqvap0,  duwnd0,     dvwnd0,    &
                  rain0,  snow0,   kbot,    ql0,        qi0,       &
                  qa0,    mc0,     Dl0,     Di0,        Da0  )

!=======================================================================
! ***** DRIVER FOR RAS
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!     is, js    - Starting indices for window
!     Time      - Time used for diagnostics [time_type]
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
  type(time_type), intent(in)                   :: Time
  integer,         intent(in)                   :: is, js
  real,            intent(in)                   :: dtime
  real,            intent(in), dimension(:,:,:) :: pres0, pres0_int
  real,            intent(in), dimension(:,:,:) :: temp0, qvap0, uwnd0, vwnd0
  logical,         intent(in), dimension(:,:)   :: coldT0
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
  real, intent(out), dimension(:,:)   :: rain0,  snow0
  
  real, intent(out), OPTIONAL, dimension(:,:,:) :: mc0,Dl0,Di0,Da0

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

 real, parameter :: p00 = 1000.0E2

 logical :: coldT,  exist    
 real    :: precip, Hl, psfc, dpcu
 integer :: ksfc,   klcl

 integer, dimension(SIZE(temp0,3)) :: ic

 real, dimension(SIZE(temp0,3)) :: &
       temp, qvap, uwnd, vwnd, pres, dtemp, dqvap, duwnd, dvwnd, &    
       ql,   qi,   qa,   Dl,   Di,   Da,    mass,  pi,    theta, &
       cp_by_dp,   dqvap_sat,  qvap_sat,    alpha, beta,  gamma, &
       dtcu, dqcu, ducu, dvcu, Dlcu, Dicu,  Dacu

 real, dimension(SIZE(temp0,3)+1) ::  pres_int, mc, pi_int, mccu

 logical, dimension(size(temp0,1),size(temp0,2)) :: rhtrig_mask
 integer, dimension(size(temp0,1),size(temp0,2)) :: kcbase

 real,    dimension(size(temp0,1),size(temp0,2)) :: &
          psfc0, t_parc, q_parc, p_parc, qs_parc    

 real,    dimension(size(temp0,1),size(temp0,2),size(temp0,3)) :: &
          qvap_sat0, dqvap_sat0

 logical :: setras, Lkbot, Lda0, Lmc0
 integer :: i, imax, j, jmax, k, kmax
 integer :: ncmax, nc, ib
 real    :: rasal, frac

!--- For extra diagnostics

 logical :: used
 real    :: dpevap, precip_ev

 real, dimension(SIZE(temp0,3)) :: &
       cup, dtevap, dqevap, dtemp_ev, dqvap_ev

 real, dimension(size(temp0,1),size(temp0,2)) :: &
       pcldb0, rain_ev0, snow_ev0

 real, dimension(size(temp0,1),size(temp0,2),size(temp0,3)) :: &
       cuprc3d, dtemp_ev0, dqvap_ev0

!=====================================================================

! --- Check to see if ras has been initialized
  if( do_init ) CALL ERROR_MESG( 'RAS',  &
                                 'ras_init has not been called', FATAL )
  
! --- Check for presence of optional arguments
  Lkbot  = PRESENT( kbot )
  Lda0   = PRESENT( Da0  )
  Lmc0   = PRESENT( mc0  )

! --- Set dimensions
  imax  = size( temp0, 1 )
  jmax  = size( temp0, 2 )
  kmax  = size( temp0, 3 )

! --- Initalize

   dtemp0 = 0.0                               
   dqvap0 = 0.0                               
   duwnd0 = 0.0                               
   dvwnd0 = 0.0                                                                        
    rain0 = 0.0   
    snow0 = 0.0

  dtemp_ev0 = 0.0
  dqvap_ev0 = 0.0
   rain_ev0 = 0.0
   snow_ev0 = 0.0
    cuprc3d = 0.0

  if ( Lda0 ) then
      Da0 = 0.0
      Dl0 = 0.0
      Di0 = 0.0
  end if
  if ( Lmc0 ) then
      mc0 = 0.0
  end if

    frac  = fracs  / dtime
    rasal = rasal0 / dtime 
       
! --- Compute saturation value of water vapor & its derivative

  CALL  ESCOMP( temp0,  qvap_sat0 )
  CALL DESCOMP( temp0, dqvap_sat0 )

  dqvap_sat0(:,:,:) = d622 * dqvap_sat0(:,:,:) * pres0(:,:,:) /        &
           ( MAX( qvap_sat0(:,:,:), pres0(:,:,:)-d378*qvap_sat0(:,:,:) ) **2 )

   qvap_sat0(:,:,:) = d622 * qvap_sat0(:,:,:) /                       &
             MAX( qvap_sat0(:,:,:), pres0(:,:,:)-d378*qvap_sat0(:,:,:) )


! --- Find LCL ---> cloud base

  if( Lkbot ) then
     do j = 1,jmax
     do i = 1,imax
        k = kbot(i,j)
          t_parc(i,j) =     temp0(i,j,k)
          q_parc(i,j) =     qvap0(i,j,k)
          p_parc(i,j) =     pres0(i,j,k)
         qs_parc(i,j) = qvap_sat0(i,j,k)
     end do
     end do
  else
          t_parc(:,:) =     temp0(:,:,kmax)
          q_parc(:,:) =     qvap0(:,:,kmax)
          p_parc(:,:) =     pres0(:,:,kmax)
         qs_parc(:,:) = qvap_sat0(:,:,kmax)
  end if

  q_parc  = MAX( q_parc, 1.0E-6   )
  q_parc  = MIN( q_parc, qs_parc )

  CALL COMP_LCL( t_parc, q_parc, p_parc, pres0, kcbase )

! --- set rh trigger
  rhtrig_mask(:,:) = q_parc(:,:) >= rh_trig*qs_parc(:,:) 

! --- Set surface pressure
  if( Lkbot ) then
     do j = 1,jmax
     do i = 1,imax
        k = kbot(i,j) + 1
        psfc0(i,j) = pres0_int(i,j,k)
     end do
     end do
  else
        psfc0(:,:) = pres0_int(:,:,kmax+1)
  end if

! --- Save cloud base pressure
  if ( id_pcldb > 0 ) then
     do j = 1,jmax
     do i = 1,imax
     klcl = kcbase(i,j)
            pcldb0(i,j) = pres0_int(i,j,klcl)
     end do
     end do
  end if

!---------------------------------------------------------------------

! --- Set order in which clouds are to be done                   
  CALL RAS_CLOUD_INDEX ( ic, ncmax )

! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  do j = 1,jmax
  do i = 1,imax
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

! --- rh trigger
  if ( .not.rhtrig_mask(i,j) ) CYCLE

! --- Pack input variables / Initalize output variables 

      temp(:) =      temp0(i,j,:)
      qvap(:) =      qvap0(i,j,:)
      uwnd(:) =      uwnd0(i,j,:)
      vwnd(:) =      vwnd0(i,j,:)
      pres(:) =      pres0(i,j,:)
  pres_int(:) =  pres0_int(i,j,:)
  qvap_sat(:) =  qvap_sat0(i,j,:) 
 dqvap_sat(:) = dqvap_sat0(i,j,:) 
      psfc    =      psfc0(i,j)
     coldT    =     coldT0(i,j)
      klcl    =     kcbase(i,j)

     dtemp(:) = 0.0                               
     dqvap(:) = 0.0                               
     duwnd(:) = 0.0                               
     dvwnd(:) = 0.0  
    precip    = 0.0   

          cup(:) = 0.0
     dtemp_ev(:) = 0.0 
     dqvap_ev(:) = 0.0 
    precip_ev    = 0.0 

  if ( Lda0 ) then
        qa(:) = qa0(i,j,:)
        ql(:) = ql0(i,j,:)
        qi(:) = qi0(i,j,:)
        Da(:) = 0.0
        Dl(:) = 0.0
        Di(:) = 0.0
  end if
  if ( Lmc0 ) then
        mc(:) = 0.0
  end if

    if( coldT ) then
     Hl = HLs
  else  
     Hl = HLv
  end if

  if( Lkbot ) then
       ksfc = kbot(i,j) + 1
  else
       ksfc = kmax + 1
  endif

! --- Thickness of layers 
  mass(1:kmax) =  pres_int(2:kmax+1) - pres_int(1:kmax) 
  mass(1:kmax) = MAX( mass(1:kmax), 1.0e-5 )

! --- Compute exner functions 
! --- at layer interfaces
  pi_int(:) = ( pres_int(:) / p00  ) ** Kappa                                   
! --- at full levels
  pi(1:kmax) = ( pi_int(2:kmax+1) * pres_int(2:kmax+1)     &
               - pi_int(1:kmax  ) * pres_int(1:kmax  ) )   &
               / ( mass(1:kmax  ) * one_plus_kappa )
  pi(1:kmax) =  MAX( pi(1:kmax), 1.0e-5 )

! --- Compute potential temperature
  theta(:) = temp(:) / pi(:)                                 

! --- Compute Cp divided by dpres                  
  cp_by_dp(1:kmax) = Cp / mass(1:kmax)

! --- Compute mass of layers              
  mass(:) = mass(:) / Grav

  setras  = .true.

! --- Test for convective instability in column             
 CALL RAS_CLOUD_EXIST( klcl, qvap,   qvap_sat, theta,    &
                       pi,   pi_int, ic,       ncmax,    &
                       Hl,   exist  )
 if ( .not. exist ) CYCLE

!---------------------------------------------------------------------
! Cloud top loop starts
!---------------------------------------------------------------------
  do nc = 1,ncmax
!---------------------------------------------------------------------

     ib  = ic(nc)
 if( ib >= klcl) CYCLE

 if ( setras ) then
! --- Compute some stuff
     alpha(:) =  qvap_sat(:) - dqvap_sat(:) * temp(:)
      beta(:) = dqvap_sat(:) * pi(:)
     gamma(:) = 1.0 / ( ( 1.0 + ( Hl * dqvap_sat(:) / Cp ) ) * pi(:) )
 endif

! --- Do adjustment
  if ( Lda0 ) then
  CALL RAS_CLOUD(                                                    &
       klcl,  ib,   rasal, frac, Hl, coldT,                          &
       theta, qvap, uwnd,  vwnd, pres_int, pi_int, pi, psfc,         &
       alpha, beta, gamma, cp_by_dp,                                 &
       dtcu,  dqcu, ducu,  dvcu, dpcu,                               &
       ql,    qi,   qa,    mccu, Dlcu, Dicu, Dacu )
  else if ( Lmc0 .and. .not.Lda0 ) then
  CALL RAS_CLOUD(                                                    &
       klcl,  ib,   rasal, frac, Hl,  coldT,                         &
       theta, qvap, uwnd,  vwnd, pres_int, pi_int, pi, psfc,         &
       alpha, beta, gamma, cp_by_dp,                                 &
       dtcu,  dqcu, ducu,  dvcu, dpcu, mccu=mccu )
  else
  CALL RAS_CLOUD(                                                    &
       klcl,  ib,   rasal, frac, Hl,  coldT,                         &
       theta, qvap, uwnd,  vwnd, pres_int, pi_int, pi, psfc,         &
       alpha, beta, gamma, cp_by_dp,                                 &
       dtcu,  dqcu, ducu,  dvcu, dpcu )
  end if

! --- For optional diagnostic output
  cup(ib) = dpcu

! --- Multiply tendencies by size of time step
  dtcu(:) =  dtcu(:) * dtime
  dqcu(:) =  dqcu(:) * dtime
  ducu(:) =  ducu(:) * dtime
  dvcu(:) =  dvcu(:) * dtime
  dpcu    =  dpcu    * dtime

  if ( Lda0 ) then
  Dacu(:) =  Dacu(:) * dtime
  Dlcu(:) =  Dlcu(:) * dtime
  Dicu(:) =  Dicu(:) * dtime
  end if

! --- Evaporate some precip

      dtevap(:) = 0.0
      dqevap(:) = 0.0
      dpevap    = 0.0

  if( evap_on .and. ( dpcu > 0.0 ) ) then

  CALL RAS_CEVAP ( ib, temp, qvap, pres, mass, qvap_sat,       &
                   dqvap_sat, psfc, Hl, dtime, ksfc, dpcu,     &
                   dtevap, dqevap,  dpevap )

      dtcu(:) =  dtcu(:) + dtevap(:) / pi(:)
      dqcu(:) =  dqcu(:) + dqevap(:)
      dpcu    = MAX(dpcu - dpevap, 0.)

  endif

! --- Update potential temperature, water vapor, winds
  theta(:) = theta(:) + dtcu(:)
   qvap(:) =  qvap(:) + dqcu(:)
   uwnd(:) =  uwnd(:) + ducu(:)
   vwnd(:) =  vwnd(:) + dvcu(:)

  if ( Lda0 ) then
    ql(:) = ql(:) + Dlcu(:)
    qi(:) = qi(:) + Dicu(:)
    qa(:) = qa(:) + Dacu(:)
  end if

! --- Recover temperature 
  temp(:) = theta(:) * pi(:)

! --- Accumulate precip 
  precip    = precip    + dpcu 
  precip_ev = precip_ev + dpevap

! --- Accumulate tendencies 
     dtemp(:) =    dtemp(:) +   dtcu(:) * pi(:)
     dqvap(:) =    dqvap(:) +   dqcu(:)
     duwnd(:) =    duwnd(:) +   ducu(:) 
     dvwnd(:) =    dvwnd(:) +   dvcu(:)
  dtemp_ev(:) = dtemp_ev(:) + dtevap(:)
  dqvap_ev(:) = dqvap_ev(:) + dqevap(:)

  if ( Lda0 ) then
    Da(:) = Da(:) + Dacu(:)
    Dl(:) = Dl(:) + Dlcu(:)
    Di(:) = Di(:) + Dicu(:)
  end if
  if ( Lmc0 ) then
    mc(:) = mc(:) + mccu(:)
  end if

  setras = .false.

 if ( setras ) then
! --- Re-Compute saturation value of water vapor & its derivative
  CALL  ESCOMP( temp,  qvap_sat )
  CALL DESCOMP( temp, dqvap_sat )
         alpha(:) = MAX( qvap_sat(:), pres(:) - d378*qvap_sat(:) )
     dqvap_sat(:) = d622 * dqvap_sat(:) * pres(:) / ( alpha(:) * alpha(:) )
      qvap_sat(:) = d622 *  qvap_sat(:) / alpha(:)
  end if
     
!---------------------------------------------------------------------
  end do
!---------------------------------------------------------------------
! Cloud top loop ends
!---------------------------------------------------------------------

! --- Unpack variables
      dtemp0(i,j,:) = dtemp(:)
      dqvap0(i,j,:) = dqvap(:) 
      duwnd0(i,j,:) = duwnd(:)
      dvwnd0(i,j,:) = dvwnd(:)

   if ( coldT ) then
      snow0(i,j) = precip
   else
      rain0(i,j) = precip
   end if

  if ( Lda0 ) then
        Da0(i,j,:) = Da(:) 
        Dl0(i,j,:) = Dl(:) 
        Di0(i,j,:) = Di(:)
  end if
  if ( Lmc0 ) then
        mc0(i,j,:) = mc(:)
  end if

! --- For extra diagnostics
 
     cuprc3d(i,j,:) =       cup(:) 
   dtemp_ev0(i,j,:) =  dtemp_ev(:)
   dqvap_ev0(i,j,:) =  dqvap_ev(:)

    if ( coldT ) then
       snow_ev0(i,j) = precip_ev
    else
       rain_ev0(i,j) = precip_ev
    end if

! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  end do
  end do
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

!---------------------------------------------------------------------
! --- Extra diagnostics
!---------------------------------------------------------------------
 
   dtemp_ev0(:,:,:) = dtemp_ev0(:,:,:) / dtime
   dqvap_ev0(:,:,:) = dqvap_ev0(:,:,:) / dtime
    snow_ev0(:,:)   =  snow_ev0(:,:)   / dtime
    rain_ev0(:,:)   =  rain_ev0(:,:)   / dtime

     if ( id_tdt_revap > 0 ) then
        used = send_data ( id_tdt_revap, dtemp_ev0, Time, is, js, 1 )
     endif
     if ( id_qdt_revap > 0 ) then
        used = send_data ( id_qdt_revap, dqvap_ev0, Time, is, js, 1 )
     endif
     if ( id_prec_revap > 0 ) then
        used = send_data ( id_prec_revap, rain_ev0+snow_ev0, Time, is, js )
     endif
     if ( id_snow_revap > 0 ) then
        used = send_data ( id_snow_revap, snow_ev0, Time, is, js )
     endif
     if ( id_prec_conv_3d > 0 ) then
        used = send_data ( id_prec_conv_3d, cuprc3d, Time, is, js, 1 )
     endif
     if ( id_pcldb > 0 ) then
        used = send_data ( id_pcldb, pcldb0, Time, is, js )
     endif

!=====================================================================
  end SUBROUTINE RAS


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
 real,    intent(in) :: hl,    psfc
 integer, intent(in) :: ic,    k
 logical, intent(in) :: coldT
 real, intent(in), dimension(:) :: theta, qvap, uwnd, vwnd, pres_int, pi_int
 real, intent(in), dimension(:) :: alf,   bet,  gam,  pi,   cp_by_dp
 real, intent(in), OPTIONAL, dimension(:) :: ql,qi,qa
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

 real, intent(out)  :: dpcu
 real, intent(out), dimension(:) :: dtcu, dqcu, ducu, dvcu
 real, intent(out), OPTIONAL, dimension(:) :: mccu, Dacu, Dlcu, Dicu

!---------------------------------------------------------------------
!    (Intent local)
!---------------------------------------------------------------------

 real, parameter :: rhmax  = 0.9999 

  logical :: lcase1,    lcase2
  real    :: wfn_crit,  ftop, rn_frac
  real    :: wfn,  akm, qs1,  uht, vht, wlq, alm 
  real    :: rasalf
  real    :: wll,  wli
  integer :: km1,  ic1, l,    iwk
  real    :: xx1,  xx2, xx3,  xx4
  real    :: ssl, dtemp, zzl, hccp, hcc,  dpib, dpit, zbase 
  real    :: dut, dvt,   dub, dvb,  wdet, dhic, hkb, hic, sic
 
 real, dimension(size(theta,1)) :: gmh, eta, hol, hst, qol

 logical :: Ldacu, Lmccu

!=====================================================================

! --- Check for presence of optional arguments
  Ldacu = PRESENT( Dacu )
  Lmccu = PRESENT( mccu ) 

! Initialize
  dtcu = 0.0
  dqcu = 0.0
  ducu = 0.0
  dvcu = 0.0
  dpcu = 0.0
  if ( Ldacu ) then
  Dacu = 0.0
  Dlcu = 0.0
  Dicu = 0.0
  end if
  if ( Lmccu ) then
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
  wll=0.0
  wli=0.0
  gmh=0.0
  eta=0.0
  hol=0.0
  hst=0.0
  qol=0.0
  km1=0
  ic1=0
  l=0
  rasalf=0.0

  km1 = k  - 1
  ic1 = ic + 1

!=====================================================================
! --- RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
!=====================================================================

! --- at cloud base
    qs1    = alf(k) + bet(k) * theta(k)
    qol(k) = MIN( qs1 * rhmax, qvap(k) )
    hol(k) = pi_int(k+1) * theta(k) * Cp + qol(k) * Hl
    eta(k) = 0.0
    zzl    = ( pi_int(k+1) - pi_int(k) ) * theta(k) * Cp
    zbase  = zzl

! --- between cloud base & cloud top
 if ( ic < km1 ) then
 do l = km1,ic1,-1
    qs1    = alf(l) + bet(l) * theta(l)
    qol(l) = MIN( qs1 * rhmax, qvap(l) )
    ssl    = zzl + pi_int(l+1) * theta(l) * Cp 
    hol(l) = ssl + qol(l) * Hl
    hst(l) = ssl + qs1    * Hl
    dtemp  = ( pi_int(l+1) - pi_int(l) ) * theta(l)
    eta(l) = eta(l+1) + dtemp * cp_div_grav
    zzl    = zzl      + dtemp * Cp    
 end do
 end if

! --- at cloud top
    qs1     = alf(ic) + bet(ic) * theta(ic)
    qol(ic) = MIN( qs1 * rhmax, qvap(ic) ) 
    ssl     = zzl + pi_int(ic1) * theta(ic) * Cp 
    hol(ic) = ssl + qol(ic) * Hl
    hst(ic) = ssl + qs1     * Hl
    dtemp   = ( pi_int(ic1) - pi(ic) ) * theta(ic)
    eta(ic) = eta(ic1) + dtemp * cp_div_grav

!=====================================================================
!     ENTRAINMENT PARAMETER
!=====================================================================

    xx1 = hol(k) - hst(ic)

    xx2 = 0.0
 do l = ic,km1
    xx2 = xx2 + ( hst(ic) - hol(l) ) * ( eta(l) - eta(l+1) )
 end do

   lcase1 = ( xx2    >  0.0      ) .and.  &
            ( xx1    >  0.0      )

   lcase2 = ( xx1    <= 0.0      ) .and. &
            ( hol(k) >  hst(ic1) ) .and. &
            ( ic1    <  k ) 

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 if ( .not.lcase1 .and. .not.lcase2 )  RETURN
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  if ( lcase1 ) then
     alm = xx1 / xx2
  else 
     alm = 0.0
  end if

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  if( Tokioka_on ) alm_min = Tokioka_con / zbase
  if( alm < alm_min ) RETURN
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!=====================================================================
!    NORMALIZED MASSFLUX
!=====================================================================

 do l = ic,km1
    eta(l) = 1.0 + alm * eta(l)
 end do
    eta(k) = 1.0

!=====================================================================
!     CLOUD WORKFUNCTION
!=====================================================================

    wfn  = 0.0
    hccp = hol(k)

 if ( ic1 <= km1 ) then
 do l = km1,ic1,-1
    hcc = hccp + ( eta(l) - eta(l+1) ) * hol(l)
    dpib = pi_int(l+1) - pi(l)
    dpit = pi(l)       - pi_int(l)
    xx1  = ( eta(l+1) * dpib + eta(l) * dpit ) * hst(l)
    xx2  =       hccp * dpib +    hcc * dpit
    wfn  = wfn + ( xx2 - xx1 ) * gam(l)
    hccp = hcc
 end do
 end if

 if ( lcase1 ) then
    wfn = wfn + gam(ic) * ( pi_int(ic1) - pi(ic) ) * &
               ( hccp - hst(ic) * eta(ic1) )
 end if

!=====================================================================
!    CRITICAL CLOUD WORK FUNCTION
!=====================================================================

      xx1  = 0.5 * ( pres_int(ic) + pres_int(ic1) )
      xx2  = pres_int(k)
      ftop = 1.0

 if ( lcase2 ) then
   if ( hst(ic1) < hst(ic) ) then
      ftop = ( hst(ic1) - hol(k) ) / ( hst(ic1) - hst(ic) )
   else
      ftop = 0.0
   endif
      xx3 = 0.5 * ( pres_int(ic1) + pres_int(ic1+1) )
      xx1 = xx3 * ( 1.0 - ftop ) + xx1 * ftop
 endif

! $$$  CALL RAS_ACRITN ( xx1, xx2, wfn_crit )

        iwk = xx1 * 0.02E-2 - 0.999999999
 if ( ( iwk > 1 ) .and. ( iwk <= 15 ) ) then
         wfn_crit = ac(iwk) + xx1 * ad(iwk)
 else if ( iwk <= 1 ) then
         wfn_crit = actop
 else if ( iwk > 15 ) then
         wfn_crit = a(15)
 end if
         wfn_crit = aratio * wfn_crit * ( xx2 - xx1 )

  wfn = wfn - wfn_crit

 lcase1 = lcase1 .and. ( wfn > 0.0 )
 lcase2 = lcase2 .and. ( wfn > 0.0 ) .and. ( ftop > 0.0 )

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  if ( .not.lcase1 .and. .not.lcase2 )  RETURN
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!=====================================================================

  if ( lcase1 ) then
      dhic = hst(ic) - hol(ic)
  else
      dhic = ( hol(k)  - hol(ic1) ) &
           - ( hol(ic) - hol(ic1) ) * ftop
  end if

  if ( lcase2 ) then
     xx1 = ftop * ( hol(ic) - hol(ic1) ) + hol(ic1)
     xx2 = ftop * ( qol(ic) - qol(ic1) ) + qol(ic1)
     sic = xx1 - xx2 * Hl
     qs1 = xx2 + dhic * ( 1.0 / Hl )
  end if

     hic = hol(ic)
     hkb = hol(k)

!=====================================================================

     wlq =  qol(k) - qs1      * eta(ic)
     uht = uwnd(k) - uwnd(ic) * eta(ic)
     vht = vwnd(k) - vwnd(ic) * eta(ic)

 do l = km1,ic,-1
     xx1 = eta(l) - eta(l+1)
     wlq = wlq + xx1 *  qol(l)
     uht = uht + xx1 * uwnd(l)
     vht = vht + xx1 * vwnd(l)
 end do

!=====================================================================
!         CALCULATE TRACER UPDRAFT PROPERTIES
!            AND CONVECTIVE SOURCE OF TRACER
!=====================================================================
 
 if ( Ldacu ) then

     wll = ql(k)
     wli = qi(k)

 do l = km1,ic,-1
     xx1 = eta(l) - eta(l+1)
     wll = wll + xx1 * ql(l)
     wli = wli + xx1 * qi(l)
 end do
 
!do cloud base level
     xx1     = 0.5 * cp_by_dp(k) / Cp / onebg
     Dlcu(k) = ( ql(km1) - ql(k) ) * xx1
     Dicu(k) = ( qi(km1) - qi(k) ) * xx1
     Dacu(k) = ( qa(km1) - qa(k) ) * xx1
     mccu(k) = eta(k)

 if ( ic1 <= km1 ) then
 do l = km1,ic1,-1
     xx1     = 0.5 * cp_by_dp(l) / Cp / onebg
     Dlcu(l) = ( eta(l+1) * ( ql(l  ) - ql(l+1) ) + &
                 eta(l  ) * ( ql(l-1) - ql(l  ) ) ) * xx1
     Dicu(l) = ( eta(l+1) * ( qi(l  ) - qi(l+1) ) + &
                 eta(l  ) * ( qi(l-1) - qi(l  ) ) ) * xx1
     Dacu(l) = ( eta(l+1) * ( qa(l  ) - qa(l+1) ) + &
                 eta(l  ) * ( qa(l-1) - qa(l  ) ) ) * xx1
     mccu(l) = eta(l)
 end do
 end if

 !do cloud top level
     xx1      = cp_by_dp(ic) / Cp / onebg
     xx2      = 0.5 * xx1
     Dlcu(ic) = ( eta(ic1) * ( ql(ic) - ql(ic1) ) * xx2 ) + &
                ( wll       - eta(ic) * ql(ic)  ) * xx1
     Dicu(ic) = ( eta(ic1) * ( qi(ic) - qi(ic1) ) * xx2 ) + &
                ( wli       - eta(ic) * qi(ic)  ) * xx1
     Dacu(ic) = ( eta(ic1) * ( qa(ic) - qa(ic1) ) * xx2 ) + &
                ( eta(ic)   - eta(ic) * qa(ic)  ) * xx1

 end if

 if ( Lmccu .and. .not.Ldacu ) then

 !do cloud base level
     mccu(k) = eta(k)

 if ( ic1 <= km1 ) then
 do l = km1,ic1,-1
   mccu(l) = eta(l)
 end do
 end if

 end if

!=======================================================================
!     CALCULATE GS AND PART OF AKM (THAT REQUIRES ETA)
!=======================================================================

     xx1      = ( theta(km1) - theta(k) ) / ( pi(k) - pi(km1) )
     hol(k)   = xx1 * ( pi_int(k) - pi(km1) ) * pi(k)   * cp_by_dp(k)
     hol(km1) = xx1 * ( pi(k) - pi_int(k)   ) * pi(km1) * cp_by_dp(km1)
     akm      = 0.0

 if ( ic1 <= km1 ) then
 do l = km1,ic1,-1
     xx1      = ( theta(l-1) - theta(l) ) * eta(l) / ( pi(l) - pi(l-1) )
     hol(l)   = xx1 * ( pi_int(l) - pi(l-1) ) * pi(l)   * cp_by_dp(l) + hol(l)  
     hol(l-1) = xx1 * ( pi(l) - pi_int(l)   ) * pi(l-1) * cp_by_dp(l-1)
     akm      = akm - hol(l)  * ( eta(l)   * ( pi(l  ) - pi_int(l) ) +  &
                                  eta(l+1) * ( pi_int(l+1) - pi(l) ) ) / pi(l)
 end do
 end if

!=======================================================================
!  PARTION CLOUD LIQUID WATER INTO PRECIP & DETRAINED WATER 
!=======================================================================

     xx1 = 0.5 * ( pres_int(ic) + pres_int(ic1) )

! $$$ CALL RAS_RNCL( xx1, rn_frac)
       rn_frac = rn_frac_top
   if ( ( xx1 >= rn_ptop ) .and. ( xx1 <= rn_pbot ) ) then
       rn_frac = ( rn_pbot - xx1 ) * rn_pfac +  rn_frac_bot
   end if
   if ( xx1 > rn_pbot ) then
       rn_frac = rn_frac_bot
   end if

       wdet = ( 1.0 - rn_frac ) * wlq
       wlq  =         rn_frac   * wlq

  if ( Ldacu ) then !detrain non-precipitated fraction of condensed vapor
    
      if (coldT) then
         Dicu(ic) = Dicu(ic) + wdet * cp_by_dp(ic)/Cp/onebg
      else
         Dlcu(ic) = Dlcu(ic) + wdet * cp_by_dp(ic)/Cp/onebg
      end if
          wdet = 0.0

  end if

!=======================================================================
!     CALCULATE GH
!=======================================================================

   xx1 = hol(ic)
 if ( lcase2 ) then
   xx1 = xx1 + ( sic - hic + qol(ic) * Hl ) * ( cp_by_dp(ic) / Cp )
 end if

   hol(ic) = xx1 - ( wdet * Hl * cp_by_dp(ic) / Cp )

 if ( lcase1 ) then
   akm = akm - eta(ic1) * ( pi_int(ic1) - pi(ic) ) * xx1 / pi(ic)
 end if

    xx2    = qol(km1) - qol(k)
    gmh(k) = hol(k) + ( xx2 * cp_by_dp(k) * Hl * 0.5 / Cp )
    akm    = akm + gam(km1) * ( pi_int(k) - pi(km1) ) * gmh(k)

 if ( ic1 <= km1 ) then
 do l = km1,ic1,-1
     xx3    = xx2
     xx2    = ( qol(l-1) - qol(l) ) * eta(l)
     xx3    = xx3 + xx2
     gmh(l) = hol(l) + ( xx3 * cp_by_dp(l) * Hl * 0.5 / Cp )
 end do
 end if

 if ( lcase2 ) then
    xx2 = xx2 + 2.0 * ( hkb - dhic - sic - qol(ic) * Hl ) / Hl
 end if

  gmh(ic) = xx1 + cp_by_dp(ic) * onebcp * ( xx2 * Hl * 0.5 + eta(ic) * dhic )

!=======================================================================
!     CALCULATE HC PART OF AKM
!=======================================================================

 if ( ic1 <= km1 ) then
    xx1 = gmh(k)
 do  l = km1,ic1,-1
   xx1 = xx1 + ( eta(l) - eta(l+1) ) * gmh(l)
   xx2 = gam(l-1) * ( pi_int(l) - pi(l-1) )
   if ( lcase2 .and. ( l == ic1 ) ) xx2 = 0.0
   akm = akm + xx1 * ( xx2 + gam(l) * ( pi(l) - pi_int(l) ) )
 end do
 end if

!=======================================================================

 if ( lcase2 ) then

      xx1 = 0.5*( pres_int(ic  ) + pres_int(ic1) )   &
          + 0.5*( pres_int(ic+2) - pres_int(ic ) ) * ( 1.0 - ftop )
      xx2 =       pres_int(ic1 )
      xx3 = 0.5*( pres_int(ic1 ) + pres_int(ic+2) )

 if ( ( xx1 >= xx2 ) .and. ( xx1 < xx3 ) ) then
        ftop = 1.0 - ( xx1 - xx2 ) / ( xx3 - xx2 )
         xx4 = cp_by_dp(ic1) / cp_by_dp(ic)
    hol(ic1) = hol(ic1) + hol(ic) * xx4
    gmh(ic1) = gmh(ic1) + gmh(ic) * xx4
    hol(ic)  = 0.0
    gmh(ic)  = 0.0
 else if ( xx1 < xx2 ) then
     ftop = 1.0
 else
     ftop = 0.0
 end if

 end if

!=======================================================================
!   MASS FLUX
!=======================================================================
 
 if ( ( akm < 0.0 ) .and. ( wlq >= 0.0 ) ) then
!jjs
  rasalf = rasal * ( pres_int(ic+1) - puplim ) /      &
                   (     psfc       - puplim )
  if (puplim > psfc) rasalf=0.
  rasalf = MAX( 0.0, rasalf )
!jjs
     wfn = -ftop * wfn * rasalf / akm
 else
     wfn = 0.0
 end if

     xx1 = ( pres_int(k+1) - pres_int(k) ) * frac
     wfn = MIN( wfn, xx1 )

!=======================================================================
!     FOR SAK CLOUDS
!=======================================================================

 if ( Ldacu ) then
      xx1 = wfn * onebg
     do l = ic,k
          Dacu(l) = Dacu(l) * xx1
          Dlcu(l) = Dlcu(l) * xx1
          Dicu(l) = Dicu(l) * xx1
          mccu(l) = mccu(l) * xx1
     end do
 end if

 if ( Lmccu .and. .not.Ldacu ) then
      xx1 = wfn * onebg
     do l = ic,k
          mccu(l) = mccu(l) * xx1
     end do
 end if
 
!=======================================================================
!     PRECIPITATION
!=======================================================================

   dpcu =  wlq * wfn * onebg

!=======================================================================
!     THETA AND Q CHANGE DUE TO CLOUD TYPE IC
!=======================================================================
 
    xx1 = wfn * onebcp
    xx2 = wfn * ( 1.0 / Hl )

 do l = ic,k
   dtcu(l) = xx1 * hol(l) / pi(l)
   dqcu(l) = xx2 * ( gmh(l) - hol(l) )
 end do

!=======================================================================
!  CUMULUS FRICTION 
!=======================================================================
 if (cufric) then

    xx1 = 0.5 * wfn * onebcp

! --- At cloud base
    xx2    = xx1 * cp_by_dp(k)
    dut    = ( uwnd(km1) - uwnd(k) )
    dvt    = ( vwnd(km1) - vwnd(k) )
   ducu(k) = dut * xx2
   dvcu(k) = dvt * xx2

! --- Between cloud base & cloud top
 do l = km1,ic1,-1
    xx2    = xx1 * cp_by_dp(l)
    dub    = dut
    dvb    = dvt
    dut    = ( uwnd(l-1) - uwnd(l) ) * eta(l)
    dvt    = ( vwnd(l-1) - vwnd(l) ) * eta(l)
   ducu(l) = ( dut + dub ) * xx2
   dvcu(l) = ( dvt + dvb ) * xx2
 end do

! --- At cloud top
    xx2      =   xx1 * cp_by_dp(ic)
    ducu(ic) = ( dut + uht + uht ) * xx2
    dvcu(ic) = ( dvt + vht + vht ) * xx2

  end if

!=======================================================================
  end SUBROUTINE RAS_CLOUD

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
    rhum(:,:) =    chi(:,:) * LOG( rhum(:,:) )
   p_lcl(:,:) = p_parc(:,:) * EXP( rhum(:,:) )

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

!#######################################################################
!#######################################################################

 SUBROUTINE RAS_CEVAP ( type,     temp,      qvap,   pres,   mass,  &
                        qvap_sat, dqvap_sat, psfc,   hl,     dtime, &
                        ksfc,     dpcu,      dtevap, dqevap, dpevap )

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
!     qvap_sat - Saturation value of qvap
!     dqvap_sat - Temperature derivative of qvap_sat
!     psfc  - surface presure
!     hl    - proper latent heat for the column
!     dtime - size of time step in seconds
!     ksfc  - index of lowest model level
!     dpcu - Precip in mm
!---------------------------------------------------------------------

  integer :: type
  real    :: dtime

  real,    intent(in), dimension(:) :: temp, qvap, pres, mass
  real,    intent(in), dimension(:) :: qvap_sat, dqvap_sat
  real,    intent(in)               :: psfc, hl, dpcu
  integer, intent(in)               :: ksfc

!---------------------------------------------------------------------
! Arguments (Intent out)
!     dtevap - Temperature change due to precip evap
!     dqevap - Water vapor change due to precip evap
!     dpevap - Amount of precip evaporated
!---------------------------------------------------------------------
  real, intent(out), dimension(:) :: dtevap, dqevap
  real, intent(out)               :: dpevap

!---------------------------------------------------------------------
!  (Intent local)
!---------------------------------------------------------------------

  real, parameter :: cem  = 0.054
  real, parameter :: ceta = -544.0E-6

  real, dimension(size(temp,1)) ::  temp_new, qvap_new

  real    :: prec, def, evef
  real    :: prec_mmph, pfac, emx
  integer :: itopp1, kmax, k

!=======================================================================

  kmax   = size(temp)
  itopp1 = type + 1

! --- Initalize
  dpevap   = 0.0
  qvap_new = qvap
  temp_new = temp

!=======================================================================

  do k = itopp1,kmax
!-----------------------------------------

! --- Compute whats available for evaporation 
   prec = MAX(dpcu - dpevap, 0.0 )

! --- Compute precipitation efficiency factor
  prec_mmph = prec * 3600.0 / dtime
  pfac      = SQRT( pres(k) / psfc )
  emx       = SQRT( cem * cfrac * prec_mmph * pfac )   
  evef      = 1.0 - EXP( ceta * dtime * emx ) 

! --- Evaporate precip where needed
  if ( ( hcevap*qvap_sat(k) >= qvap(k) ) .and.  &
          ( prec            > 0.0      ) .and.  &
          ( ksfc            > k        ) ) then
            def = ( hcevap*qvap_sat(k) - qvap(k) ) /    &
                  ( 1.0 + (hl * hcevap * dqvap_sat(k) / Cp ) )
            def = evef*def
            def = MIN( def, prec/mass(k) )
    qvap_new(k) = qvap(k) + def
    temp_new(k) = temp(k) - (def * hl/Cp)
         dpevap = dpevap + def * mass(k)
  end if

!-----------------------------------------
  end do

! --- Changes to temperature and water vapor from evaporation
  dtevap(:) = temp_new(:) - temp(:) 
  dqevap(:) = qvap_new(:) - qvap(:) 

!=======================================================================
  end SUBROUTINE RAS_CEVAP

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

!#####################################################################
!#####################################################################


 SUBROUTINE RAS_CLOUD_EXIST( k,      qvap,   qsat, theta,    &
                             pifull, pihalf, nc,   ncmax,    &
                             Hl,     exist  )
!=======================================================================
 implicit none
!=======================================================================

 integer, intent(in)               :: k, ncmax
 real,    intent(in)               :: Hl
 integer, intent(in), dimension(:) :: nc 
 real,    intent(in), dimension(:) :: qvap, qsat, theta, pifull, pihalf 

 logical, intent(out) :: exist

!---------------------------------------------------------------------
 real, parameter :: rhmax  = 0.9999 

 real, dimension(size(theta)) :: ssl, hst
 real                         :: zzl, hol_k, qol_k
 integer                      :: km1, ic, ic1, l

!=======================================================================

 ic  = MINVAL( nc(1:ncmax) )

 km1 = k  - 1
 ic1 = ic + 1

! --- at cloud base
    qol_k = MIN( qsat(k) * rhmax, qvap(k) )
    ssl(k) = pihalf(k+1) * theta(k) * Cp 
    hol_k  = ssl(k) +  qol_k  * Hl
    hst(k) = ssl(k) + qsat(k) * Hl
    zzl    = ( pihalf(k+1) - pihalf(k) ) * theta(k) * Cp

! --- between cloud base & cloud top
 do l = km1,ic1,-1
    ssl(l) = zzl + pihalf(l+1) * theta(l) * Cp 
    hst(l) = ssl(l) + qsat(l) * Hl
    zzl    = zzl + ( pihalf(l+1) - pihalf(l) ) * theta(l) * Cp    
 end do

! --- at cloud top
    ssl(ic) = zzl  + pihalf(ic1) * theta(ic) * Cp 
    hst(ic) = ssl(ic) + qsat(ic) * Hl

! --- test
    exist = hol_k > MINVAL( hst(ic:k) )

!=======================================================================
 end SUBROUTINE RAS_CLOUD_EXIST

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
  end do

    dqvap_prec  = sum_dqvap + precip(i,j)
     
    if (present(dql)) then
  do k = 1,kmax
    dqvap_prec = dqvap_prec + (dql(i,j,k)+dqi(i,j,k))*mass(i,j,k)      
  end do
    end if  

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
  end MODULE RAS_MOD
