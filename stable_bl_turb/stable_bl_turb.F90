  MODULE STABLE_BL_TURB_MOD

!=======================================================================
  use    Utilities_Mod, ONLY: FILE_EXIST, OPEN_FILE, ERROR_MESG, FATAL,     &
                              get_my_pe, read_data, write_data, CLOSE_FILE, &
                              check_nml_error
 use  Diag_Manager_Mod, ONLY: register_diag_field, send_data
 use  Time_Manager_Mod, ONLY: time_type
 use     Constants_Mod, ONLY: grav, vonkarm, omega
 use Monin_Obukhov_Mod, ONLY: stable_mix

 implicit none
 private
 public :: STABLE_BL_TURB, STABLE_BL_TURB_INIT


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  character(len=128) :: version = '$Id: stable_bl_turb.F90,v 1.2 2003/04/09 21:02:16 fms Exp $'
  character(len=128) :: tag = '$Name: inchon $'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 logical :: do_init = .true.

!---------------------------------------------------------------------
! --- CONSTANTS
!---------------------------------------------------------------------

  real :: oalsm, oalsh

!---------------------------------------------------------------------
! --- NAMELIST
!---------------------------------------------------------------------

  real    :: alpha    = 0.5
  real    :: alsm     = 150.0
  real    :: alsh     = 150.0
  real    :: fmin     = 5.0e-5
  real    :: hpbl_cap = 2000.
  logical :: taper_ak = .true.

  NAMELIST / stable_bl_turb_nml / alpha, alsm, alsh, fmin, &
                                  hpbl_cap, taper_ak

   real, dimension(:,:), allocatable :: sgsmtn

!---------------------------------------------------------------------
! DIAGNOSTICS FIELDS 
!---------------------------------------------------------------------

integer :: id_z_sbl, id_f_sbl

character(len=14) :: mod_name = 'stable_bl_turb'

real :: missing_value = -999.

!---------------------------------------------------------------------
 contains

!#######################################################################

 SUBROUTINE STABLE_BL_TURB( is, js, Time, theta,  um,  vm,  zhalf, zfull, &
                            u_star, b_star, lat, akm, akh,  mask,  &
                            edttop  )

!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!       is, js   -  Starting indices for window
!       Time     -  Time used for diagnostics [time_type]
!       theta    -  Potential temperature
!       um, vm   -  Wind components
!       zhalf    -  Height at half levels
!       zfull    -  Height at full levels
!       u_star   -  surface friction velocity (m/s)
!       b_star   -  surface buoyancy
!       lat      -  latitude in radians
!       mask     -  OPTIONAL; floating point mask (0. or 1.) designating
!                   where data is present
!       edttop   -  OPTIONAL; maximum altitude for stable bl to operate
!                   as determined by the EDT turbulence scheme
!---------------------------------------------------------------------
! Arguments (Intent out)
!       akm  -  mixing coefficient for momentum
!       akh  -  mixing coefficient for heat and moisture
!---------------------------------------------------------------------

  type(time_type), intent(in)                    :: Time
  integer,         intent(in)                    :: is, js
  real,            intent(in),  dimension(:,:)   :: u_star, b_star, lat
  real,            intent(in),  dimension(:,:,:) :: theta,  um,     vm 
  real,            intent(in),  dimension(:,:,:) :: zhalf,  zfull
  real,            intent(out), dimension(:,:,:) :: akm,    akh

  real, intent(in), OPTIONAL, dimension(:,:,:) :: mask
  real, intent(in), OPTIONAL, dimension(:,:)   ::         edttop

!---------------------------------------------------------------------

  real, dimension(SIZE(um,1),SIZE(um,2),SIZE(um,3)-1) ::     &
        dsdzh, shear, buoync, Ri, fm, fh, lm, lh, xxm1, xxm2, phi

  real, dimension(SIZE(um,1),SIZE(um,2),SIZE(um,3)) ::       &
        zfunc

  real, dimension(SIZE(um,1),SIZE(um,2)) ::       &
        fcor, hpbl, z_sbl, f_sbl

  integer :: ix, jx, kx, i, j, k
  integer :: kxp, kxm
  integer :: shape1(1), shape3(3)
  logical :: used

  real, dimension(SIZE(um,1)*SIZE(um,2)*(SIZE(um,3)-1)) ::     &
        Ri_1d, phi_1d

!=======================================================================
! --- Initalize
!=======================================================================

! --- Check to see if STABLE_BL_TURB has been initialized
  if( do_init ) CALL ERROR_MESG( ' STABLE_BL_TURB',                          &
                                 ' STABLE_BL_TURB_INIT has not been called', &
                                   FATAL)

! --- Zero out output arrays
    akm(:,:,:) = 0.0
    akh(:,:,:) = 0.0
  z_sbl(:,:)   = 0.0
  f_sbl(:,:)   = 0.0

! --- Go home if no stable points at surface
  if( COUNT( b_star < 0.0 ) == 0 )  GO TO 100

! --- Set dimensions etc
  ix  = SIZE( um, 1 )
  jx  = SIZE( um, 2 )
  kx  = SIZE( um, 3 )
  kxp = kx + 1
  kxm = kx - 1

  shape1 =    ix * jx * kxm
  shape3 = (/ ix,  jx,  kxm /)

!====================================================================
! --- DYNAMIC HEIGHT      
!====================================================================

  fcor(:,:) = 2.0 * omega * SIN( lat(:,:) )
  fcor(:,:) = ABS( fcor(:,:) )
  fcor(:,:) = MAX( fcor(:,:), fmin )
  hpbl(:,:) = alpha * u_star(:,:) / fcor(:,:)

! --- impact of topography 
     hpbl(:,:) = MAX(hpbl(:,:), sgsmtn(is:ix, js:jx))

! --- limit from edt turbulence
  if( PRESENT( edttop ) ) then
      hpbl(:,:) = min ( hpbl(:,:), edttop(:,:) )
  end if

! --- bound
  hpbl(:,:) = MIN( hpbl(:,:), hpbl_cap )

!====================================================================
! --- COMPUTE RICHARDSON NUMBER                 
!====================================================================

! --- D( )/DZ OPERATOR  
  
  dsdzh(:,:,1:kxm) = 1.0 / ( zfull(:,:,1:kxm) - zfull(:,:,2:kx) )

! --- WIND SHEAR SQUARED

  xxm1(:,:,1:kxm) = dsdzh(:,:,1:kxm)*( um(:,:,1:kxm) - um(:,:,2:kx) )
  xxm2(:,:,1:kxm) = dsdzh(:,:,1:kxm)*( vm(:,:,1:kxm) - vm(:,:,2:kx) )

  shear(:,:,:) = xxm1(:,:,:)*xxm1(:,:,:) + xxm2(:,:,:)*xxm2(:,:,:)

! --- BUOYANCY                 

  xxm1(:,:,1:kxm) =       theta(:,:,1:kxm) - theta(:,:,2:kx)
  xxm2(:,:,1:kxm) = 0.5*( theta(:,:,1:kxm) + theta(:,:,2:kx) )
 
 
  buoync(:,:,:) = grav * dsdzh(:,:,:) * xxm1(:,:,:) / xxm2(:,:,:)

! --- RICHARDSON NUMBER

  Ri(:,:,:) = buoync(:,:,:) / shear(:,:,:)   

!====================================================================
! --- MASK OUT UNDERGROUND VALUES FOR ETA COORDINATE
!====================================================================

  if( PRESENT( mask ) ) then
     shear(:,:,1:kxm) =  shear(:,:,1:kxm) * mask(:,:,2:kx) 
    buoync(:,:,1:kxm) = buoync(:,:,1:kxm) * mask(:,:,2:kx) 
        Ri(:,:,1:kxm) =     Ri(:,:,1:kxm) * mask(:,:,2:kx) 
  endif

!====================================================================
! --- STABILITY FUNCTIONS                 
!====================================================================

  Ri_1d = RESHAPE( Ri, shape1 )

  CALL STABLE_MIX( Ri_1d, phi_1d )

  phi = RESHAPE( phi_1d, shape3 )

  fm(:,:,:) = 0.0
  fh(:,:,:) = 0.0

  where ( Ri(:,:,:) > 0.0 )
          fm(:,:,:) = phi(:,:,:)
          fh(:,:,:) =  fm(:,:,:) 
  endwhere

!====================================================================
! --- MIXING LENGTHS                 
!====================================================================

 do k = 1,kxm
   xxm1(:,:,k) = 1.0 / ( vonkarm*( zhalf(:,:,k+1) - zhalf(:,:,kxp) ) )
 end do 

  lm(:,:,:) = 1.0 / ( xxm1(:,:,1:kxm) + oalsm )
  lh(:,:,:) = 1.0 / ( xxm1(:,:,1:kxm) + oalsh )

!====================================================================
! --- MIXING COEFFICENTS                 
!====================================================================

  shear(:,:,:) = SQRT( shear(:,:,:) )

! --- Momentum
  xxm1(:,:,:)    = lm(:,:,:) * lm(:,:,:) * fm(:,:,:)
   akm(:,:,2:kx) = xxm1(:,:,1:kxm) * shear(:,:,1:kxm) 

! --- Heat and Moisture
  xxm1(:,:,:)    = lh(:,:,:) * lh(:,:,:) * fh(:,:,:)
   akh(:,:,2:kx) = xxm1(:,:,1:kxm) * shear(:,:,1:kxm)

!====================================================================
! --- BOUNDS                 
!====================================================================

! --- Use only below pbl top
 do k = 1,kx
   zfunc(:,:,k) = ( 1.0 - ( zhalf(:,:,k) - zhalf(:,:,kxp) ) / max(1.,hpbl(:,:)) )
 end do   
   zfunc(:,:,:) = MAX( zfunc(:,:,:), 0.0 )

  if( .not.taper_ak ) then
  where( zfunc(:,:,:) > 0.0 )
         zfunc(:,:,:) = 1.0
  endwhere
  endif

     akm(:,:,:) = zfunc(:,:,:) * akm(:,:,:)
     akh(:,:,:) = zfunc(:,:,:) * akh(:,:,:)

! --- Use only where b_star < 0
  do k = 1,kx
  where( b_star(:,:)  >= 0.0 )
            akm(:,:,k) = 0.0
            akh(:,:,k) = 0.0
  endwhere
  end do

!====================================================================
! --- Extra diagnostics
!====================================================================

  where( b_star(:,:)  < 0.0 .and. hpbl (:,:) > 0.0 )
          z_sbl(:,:) = hpbl(:,:)
          f_sbl(:,:) = 1.0
  endwhere
  

  100 continue
 
  if ( id_z_sbl > 0 ) then
     used = send_data ( id_z_sbl, z_sbl, Time, is, js )
  endif
  if ( id_f_sbl > 0 ) then
     used = send_data ( id_f_sbl, f_sbl, Time, is, js )
  endif

!=======================================================================
  end SUBROUTINE STABLE_BL_TURB

!#######################################################################

! SUBROUTINE STABLE_BL_TURB_INIT ( axes, Time )
 SUBROUTINE STABLE_BL_TURB_INIT ( axes, Time, sgsmtn_in )
!=======================================================================
                   
 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time
 real, dimension(:,:), intent(in) :: sgsmtn_in

!       sgsmtn_in   -            sub grid scale topography variance
 integer :: unit, io, ierr

!=======================================================================

!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------

  if( FILE_EXIST( 'input.nml' ) ) then
! -------------------------------------
   unit = OPEN_FILE ( file = 'input.nml', action = 'read' )
   ierr = 1
   do while( ierr .ne. 0 )
   READ ( unit,  nml = stable_bl_turb_nml, iostat = io, end = 10 ) 
   ierr = check_nml_error (io, 'stable_bl_turb_nml')
   end do
10 continue
   CALL CLOSE_FILE( unit )
! -------------------------------------
  end if

!---------------------------------------------------------------------
! --- Output version
!---------------------------------------------------------------------

  unit = OPEN_FILE ( file = 'logfile.out', action = 'append' )
  if ( get_my_pe() == 0 ) then
       WRITE( unit,'(/,80("="),/(a))') trim(version), trim(tag)
       WRITE( unit, nml = stable_bl_turb_nml ) 
  endif
  CALL CLOSE_FILE( unit )

!--------------------------------------------------------------------
!    save the sgs topography variance for later use.
!--------------------------------------------------------------------
   allocate (sgsmtn (size(sgsmtn_in,1), size(sgsmtn_in,2)) )
   sgsmtn(:,:) = sgsmtn_in(:,:) 

!---------------------------------------------------------------------
! --- CONSTANTS
!---------------------------------------------------------------------

  oalsm = 1.0 / alsm
  oalsh = 1.0 / alsh

!---------------------------------------------------------------------
! --- initialize quantities for diagnostics output
!---------------------------------------------------------------------

   id_z_sbl = register_diag_field ( mod_name, &
     'z_sbl', axes(1:2), Time, &
     'Depth of stable boundary layer',              'm' )

   id_f_sbl = register_diag_field ( mod_name, &
     'f_sbl', axes(1:2), Time, &
     'Frequency of stable boundary layer',          ' ' )

!---------------------------------------------------------------------
 do_init = .false.
!=======================================================================
 end SUBROUTINE STABLE_BL_TURB_INIT

!#######################################################################

  end MODULE STABLE_BL_TURB_MOD
