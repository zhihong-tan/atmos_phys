module topo_drag_mod

!==========================================================================
! TOPOGRAPHIC DRAG CLOSURE FOR GENERAL CIRCULATION MODELS -- Garner (2001)
!==========================================================================

!--------------------------------------------------------------------------
!  Calculates horizontal velocity tendency due to topographic drag
!--------------------------------------------------------------------------

  use       Fms_Mod, only: FILE_EXIST, OPEN_NAMELIST_FILE, ERROR_MESG, FATAL, NOTE, &
                           READ_DATA, WRITE_DATA, CLOSE_FILE, mpp_pe, mpp_root_pe, &
                           write_version_number, stdlog, open_restart_file
  use fms_io_mod,    only: get_restart_io_mode, write_data
  use Constants_Mod, only: Grav,Cp_Air,Rdgas,Radius,Pi,Radian

  implicit none

  private

  logical :: do_init = .true.
  logical :: do_restart_write = .true.

  character(len=128) :: version = '$Id: topo_drag.F90,v 11.0 2004/09/28 19:25:02 fms Exp $'
  character(len=128) :: tagname = '$Name: khartoum $'
  logical            :: module_is_initialized = .false.

! horizontal array size

  integer :: nlon,nlat

! arrays defined by topo_drag_init:

  real,allocatable,dimension(:)     :: flat
  real,allocatable,dimension(:,:)   :: t11,t21,t12,t22        ! drag tensor
  real,allocatable,dimension(:,:)   :: umin,umax,vmin,vmax

! parameters:

!RSH ADD:
  real,parameter :: u0=1.        ! arbitrary velocity scale for diagnostics
  real,parameter :: xl=80.e3     ! arbitrary horiz length scale for diagnostics
  real,parameter :: ro=1.2       ! arbitrary density scale for diagnostics
  real,parameter :: resfac=3.    ! residual flux is cut off when h(residual) = h/resfac

! parameters in namelist (topo_drag_nml):

  real :: &
    frcrit=1.0   &     ! critical value of Froude # for nonlinear flow
   ,anonlin=7.0  &     ! amplitude of nonpropagating drag
!RSH ADD:
!  ,gamma=1.2    &     ! exponent in aspect ratio power law
   ,gamma=1.8    &     ! exponent in aspect ratio power law
!RSH ADD:
   ,epsi=0.0     &     ! exponent in distribution power law
   ,beta=1.0     &     ! bluntness of topo features
!RSH  ,gamma=1.8    &     ! exponent in aspect ratio power law
   ,zref_fac=1.0 &     ! to adjust level separating breaking/laminar flow
   ,no_drag_frac=0.15  ! fraction of atmosphere where wave breaking is disallowed (PBL)
  logical :: &
   calculate_pbl_top=.true. ! calculate pbl top in this module, rather
                            ! than using input array
  logical :: do_netcdf_restart = .true. ! use netCDF version of the restart file

  NAMELIST /topo_drag_nml/ do_netcdf_restart, &
     & frcrit,anonlin,beta,gamma, &
!    & zref_fac,no_drag_frac
!RSH ADD:
     & epsi,                 &
     & zref_fac,no_drag_frac, &
     & calculate_pbl_top

  public topo_drag, topo_drag_init, topo_drag_end

contains

!#############################################################################      

  subroutine topo_drag(is,js,uwnd,vwnd,atmp,pfull,phalf,zfull,zhalf,    &
!   taux,tauy,dtaux,dtauy,taus)
    z_pbl, taux,tauy,dtaux,dtauy,taus)

    integer,intent(in) :: is,js

!   INPUT
!   -----

!   UWND     Zonal wind (dimensioned IDIM x JDIM x KDIM)
!   VWND     Meridional wind (dimensioned IDIM x JDIM x KDIM)
!   ATMP     Temperature at full model levels (IDIM x JDIM x KDIM)
!   PFULL    Pressure at full model levels (IDIM x JDIM x KDIM)
!   PHALF    Pressure at half model levels (IDIM x JDIM x KDIM+1)
!   ZFULL    Height at full model levels (IDIM x JDIM x KDIM)
!   ZHALF    Height at half model levels (IDIM x JDIM x KDIM+1)

    real,intent(in),dimension(:,:,:) :: uwnd,vwnd,atmp
    real,intent(in),dimension(:,:,:) :: pfull,phalf,zfull,zhalf
    real,intent(in),dimension(:,:  ) :: z_pbl                    

!   OUTPUT
!   ------

!   TAUX,TAUY    Base momentum flux in kg/m/s^2 (IDIM x JDIM) for diagnostics
!   DTAUX,DTAUY  Tendency of the vector wind in m/s^2 (IDIM x JDIM x KDIM)
!   TAUS         normalized, "clipped" saturation momentum flux at 1/2 levels

    real,intent(out),dimension(:,:)   :: taux,tauy
    real,intent(out),dimension(:,:,:) :: dtaux,dtauy,taus

    integer,dimension(size(pfull,1),size(pfull,2)) :: ktop
    integer,dimension(size(pfull,1),size(pfull,2)) :: ktop2
    real,dimension(size(pfull,1),size(pfull,2)) :: taulin,taup,taun,frulo,fruhi,ptop
    real,dimension(size(phalf,1),size(phalf,2),size(phalf,3)) :: tausat
    real,dimension(size(taux,1),size(taux,2)) :: taub

    integer :: i,idim,jdim,kdim,k,ibeg,iend,jbeg,jend
    integer :: j
    integer  :: nbad

    idim = SIZE(uwnd,1)
    jdim = SIZE(uwnd,2)
    kdim = SIZE(uwnd,3)

!   get k-index at top of estimated boundary layer

! if (calculate_pbl_top) then
    nbad = COUNT (z_pbl == -999.)
!   print *, 'pe=,js =,  nbad = ', get_my_pe(), js, nbad
  if (calculate_pbl_top .or. nbad == idim*jdim) then
    call get_pbl(is,js,uwnd,vwnd,atmp,zfull,zhalf,pfull,ktop)
  else
    call compute_pbl_top_index (is, js, zfull, z_pbl, ktop)
  endif

!   if (get_my_pe() == 15) then
!   do j=1,jdim
!   do i=1,idim
!   if (ktop(i,j) /= ktop2(i,j) ) then
!   print *, i,j, 'pe = ', get_my_pe(), 'top indices:', ktop(i,j), ktop2(i,j), z_pbl(i,j)
!   endif
!   end do
!   end do
!   endif

!   calculate base flux

    call base_flux (                                                    &
                                              is, js, uwnd, vwnd, atmp, &
                    taux, tauy, taub, taulin, taup, taun, frulo, fruhi, &
                               dtaux, dtauy, zfull, zhalf, pfull, ktop )

!   calculate saturation flux profile

    call satur_flux (                                                   &
                                                      uwnd, vwnd, atmp, &
                          taux, tauy, taub, taup, tausat, frulo, fruhi, &
                               dtaux, dtauy, zfull, zhalf, phalf, ktop )

!   put saturation flux profile into taus for diagnostics

    do k=1,kdim
      taus(:,:,k) = .5*(tausat(:,:,k) + tausat(:,:,k+1))
    enddo

!   calculate momentum tendency

    call topo_drag_tend (                                               &
                                                      uwnd, vwnd, atmp, &
                                taux, tauy, tausat, taulin, taup, taun, &
                              dtaux, dtauy, zfull, zhalf, pfull, phalf )

!   put total drag into tau for diagnostics

    taux = taux*(taup + taun)/taulin
    tauy = tauy*(taup + taun)/taulin
!    taux = taux*taup/taulin
!    tauy = tauy*taup/taulin

    return
  endsubroutine topo_drag

  !=====================================================================
                                  
  subroutine base_flux (                                                &
                                              is, js, uwnd, vwnd, atmp, &
                    taux, tauy, taub, taulin, taup, taun, frulo, fruhi, &
                               dtaux, dtauy, zfull, zhalf, pfull, ktop )

    integer,intent(in) :: is,js
    integer,intent(in),dimension(:,:)  :: ktop
    real,intent(in), dimension(:,:,:)  :: uwnd,vwnd,atmp,zfull,zhalf,pfull
    real,intent(out),dimension(:,:)    :: taux,tauy,taub,taulin,taup,taun,frulo,fruhi
    real,intent(out),dimension(:,:,:)  :: dtaux,dtauy

    integer :: idim,jdim,kdim
    integer :: i,j,k,id,jd,kr,kbl
    real :: dzfull,dzhalf,dzhalf1,dzhalf2
    real :: rnormal,gterm,rfac,bfac
    real :: frmin,frmax,frumin,frumax,fruclp,frumin1,frumax1,fruclp1,frusat
    real :: usat,bfreq2,vvtau

    idim = size(uwnd,1)
    jdim = size(uwnd,2)
    kdim = size(uwnd,3)

!   get curvature of wind at full levels (for WKB correction of wavelength)

    dtaux = 0.
    dtauy = 0.

    do j=1,jdim
      do i=1,idim
        kbl = ktop(i,j)
        do k=2,kbl
          dzfull = zhalf(i,j,k) - zhalf(i,j,k+1)
          dzhalf1 = zfull(i,j,k-1) - zfull(i,j,k)
          dzhalf2 = zfull(i,j,k) - zfull(i,j,k+1)
          dtaux(i,j,k) = ((uwnd(i,j,k-1) - uwnd(i,j,k  ))/dzhalf1 &
                        - (uwnd(i,j,k  ) - uwnd(i,j,k+1))/dzhalf2)/dzfull
          dtauy(i,j,k) = ((vwnd(i,j,k-1) - vwnd(i,j,k  ))/dzhalf1 &
                        - (vwnd(i,j,k  ) - vwnd(i,j,k+1))/dzhalf2)/dzfull
        enddo
      enddo
    enddo

!   get base flux at k = ktop

    do j=1,jdim
      jd=j+js-1
      do i=1,idim
        id=i+is-1

!       get low-level buoyancy frequency and density

        kbl = ktop(i,j)
        dzhalf = zfull(i,j,kbl) - zfull(i,j,kbl+1)
        bfreq2 = Grav*((atmp(i,j,kbl) - atmp(i,j,kbl+1))/dzhalf + Grav/Cp_Air) & 
                 / (.5*(atmp(i,j,kbl) + atmp(i,j,kbl+1)))
        rfac = pfull(i,j,kbl)/(Rdgas*atmp(i,j,kbl))
        bfac = sqrt(max(1.e-6, bfreq2))

!       get maximum drag

        taux(i,j) = (uwnd(i,j,kbl)*t11(id,jd) + vwnd(i,j,kbl)*t21(id,jd)) &
            * rfac * bfac
        tauy(i,j) = (uwnd(i,j,kbl)*t12(id,jd) + vwnd(i,j,kbl)*t22(id,jd)) &
            * rfac * bfac
        taub(i,j) = max(1.e-8, sqrt(taux(i,j)**2 + tauy(i,j)**2))

!       get min/max Froude numbers based on surface flow

        vvtau = max(1.e-8, -(uwnd(i,j,kbl)*taux(i,j) + vwnd(i,j,kbl)*tauy(i,j)))
        frmin = abs(umin(id,jd)*taux(i,j) + vmin(id,jd)*tauy(i,j)) / vvtau * bfac
        frmax = abs(umax(id,jd)*taux(i,j) + vmax(id,jd)*tauy(i,j)) / vvtau * bfac

!       get linear momentum flux associated with min/max Froude numbers

        vvtau = vvtau/taub(i,j)

        usat = max(1.e-6, sqrt(rfac/ro*vvtau**3/bfac/xl))
        frusat = frcrit*usat

        frumin = frmin*usat      ! linear momentum flux = constant as per EP theorem
        frumax = frmax*usat
        frumax = max(frumax,frumin+1.e-4)
        fruclp = min(frumax,max(frumin,frusat))

        frumin1 = min(frumin,resfac*frusat)      ! in order to cut off residual flux
        frumax1 = min(frumax,resfac*frusat)
        frumax1 = max(frumax1,frumin1+1.e-4)
        fruclp1 = min(frumax1,max(frumin1,frusat))

!       get total drag in linear limit

!RSH     rnormal = usat**gamma*2.*gamma/(frumax**(2.*gamma) - frumin**(2.*gamma))
!RSH     gterm = (frumax1**(gamma-beta) - fruclp1**(gamma-beta))/(gamma-beta)
!RSH2    rnormal = u0**gamma*(2.*gamma-epsi)/(frumax**(2.*gamma-epsi) -&
         rnormal = usat**gamma*(2.*gamma-epsi)/(frumax**(2.*gamma-epsi) -&
                   frumin**(2.*gamma-epsi))
         gterm = (frumax1**(gamma-epsi-beta) - fruclp1**   &
                   (gamma-epsi-beta))/(gamma-epsi-beta)

!RSH    taulin(i,j) = rnormal*(frumax**(gamma+2.) - frumin**(gamma+2.))/(gamma+2.)
        taulin(i,j) = rnormal*(frumax**(gamma-epsi+2.) - frumin**(gamma-epsi+2.))/(gamma-epsi+2.)

!       get propagating and nonpropagating parts of total drag

        taup(i,j) = rnormal &
!RSH      *( (fruclp**(gamma+2.) - frumin**(gamma+2.))/(gamma+2.) &
!RSH      + frusat**(2.+beta)*gterm )
           *( (fruclp**(gamma-epsi+2.) - frumin**(gamma-epsi+2.))/(gamma-epsi+2.) &
            + frusat**(2.+beta)*gterm )

        taun(i,j) = anonlin*usat/(1.+beta)*rnormal &
!RSH      *( (frumax**(gamma+1.) - fruclp**(gamma+1.))/(gamma+1.) &
          *( (frumax**(gamma-epsi+1.) - fruclp**(gamma-epsi+1.))/(gamma-epsi+1.) &
            - frusat**(1.+beta)*gterm )

!       save square root of min/max momentum flux

        frulo(i,j) = frumin
        fruhi(i,j) = frumax

      enddo
    enddo

    return
  endsubroutine base_flux

  !=====================================================================

  subroutine satur_flux (                                               &
                                                      uwnd, vwnd, atmp, &
                          taux, tauy, taub, taup, tausat, frulo, fruhi, &
                               dtaux, dtauy, zfull, zhalf, phalf, ktop )

    integer,intent(in),dimension (:,:) :: ktop
    real,intent(in), dimension (:,:,:) :: uwnd,vwnd,atmp,zfull,zhalf,phalf
    real,intent(in), dimension (:,:,:) :: dtaux,dtauy
    real,intent(in), dimension (:,:)   :: frulo,fruhi,taux,tauy,taub,taup
    real,intent(out),dimension (:,:,:) :: tausat

    integer :: i,j,k,kt
    integer :: idim,jdim,kdim
    real :: dzhalf,rnormal,gterm,rfac,bfac,xl1
    real :: frumin,frumax,fruclp,frumin1,frumax1,fruclp1,frusat
    real :: usat,bfreq2,vvtau,dv2dz2

    idim = SIZE(uwnd,1)
    jdim = SIZE(uwnd,2)
    kdim = SIZE(uwnd,3)

!   get vertical profile of propagating part of momentum flux

    do k=kdim,2,-1
      do j=1,jdim
        do i=1,idim

!         get buoyancy frequency, velocity and density at half levels

          dzhalf = zfull(i,j,k-1) - zfull(i,j,k)
          bfreq2 = Grav*((atmp(i,j,k-1) - atmp(i,j,k))/dzhalf + Grav/Cp_Air) & 
                    /(.5*(atmp(i,j,k-1) + atmp(i,j,k)))
          bfac = sqrt(max(bfreq2, 1.e-6))
          rfac = phalf(i,j,k)/(Rdgas*.5*(atmp(i,j,k-1) + atmp(i,j,k)))
          vvtau = max(0.,-.5*((uwnd(i,j,k-1) + uwnd(i,j,k))*taux(i,j) &
                            + (vwnd(i,j,k-1) + vwnd(i,j,k))*tauy(i,j)))/taub(i,j)
          dv2dz2 = -.5*((dtaux(i,j,k-1) + dtaux(i,j,k))*taux(i,j) &
                      + (dtauy(i,j,k-1) + dtauy(i,j,k))*tauy(i,j))/taub(i,j)
          xl1 = xl*max(0.5, 1.0 - vvtau*dv2dz2/(bfac*bfac))

!         get min/max and critical momentum flux values at 1/2 levels

          usat = max(1.e-6, sqrt(rfac/ro*vvtau**3/bfac/xl1))
          frusat = frcrit*usat

          frumin = frulo(i,j)     ! get momentum flux and clip to saturation value
          frumax = fruhi(i,j)
          fruclp = min(frumax,max(frumin,frusat))

          frumin1 = min(frumin,resfac*frusat)  ! in order to cut off residual flux
          frumax1 = min(frumax,resfac*frusat)
          frumax1 = max(frumax1,frumin1+1.e-4)
          fruclp1 = min(frumax1,max(frumin1,frusat))

!         get propagating part of momentum flux (from WKB or EP)

!RSH      rnormal = usat**gamma*2.*gamma/(frumax**(2.*gamma) - frumin**(2.*gamma))
!RSH      gterm = (frumax1**(gamma-beta) - fruclp1**(gamma-beta))/(gamma-beta)
!RSH2    rnormal = u0**gamma*(2.*gamma-epsi)/(frumax**(2.*gamma-epsi) - frumin**(2.*gamma-epsi))
          rnormal = usat**gamma*(2.*gamma-epsi)/(frumax**(2.*gamma-epsi) - frumin**(2.*gamma-epsi))
          gterm = (frumax1**(gamma-epsi-beta) - fruclp1**(gamma-epsi-beta))/(gamma-epsi-beta)

          tausat(i,j,k) = rnormal &
!RSH        *( (fruclp**(gamma+2.) - frumin**(gamma+2.))/(gamma+2.) &
!RSH        + frusat**(2.+beta)*gterm )
            *( (fruclp**(gamma-epsi+2.) - frumin**(gamma-epsi+2.))/(gamma-epsi+2.) &
            + frusat**(2.+beta)*gterm )
        enddo
      enddo
    enddo

!   make propagating flux constant with height in zero-drag layer

    do k=kdim+1,1,-1
      where (k >= ktop(:,:))
        tausat(:,:,k) = taup(:,:)
      endwhere
    enddo

!   clip momentum flux again if incident value is smaller

    do k=kdim,2,-1
      tausat(:,:,k) = min(tausat(:,:,k),tausat(:,:,k+1))
    enddo

!   linear momentum flux at highest model level

    tausat(:,:,1) = 0.             ! use all forcing
!    tausat(:,:,1) = tausat(:,:,2)  ! let remaining flux escape thru upper boundary

    return
  endsubroutine satur_flux

  !=====================================================================

  subroutine topo_drag_tend (                                           &
                                                      uwnd, vwnd, atmp, &
                                taux, tauy, tausat, taulin, taup, taun, &
                              dtaux, dtauy, zfull, zhalf, pfull, phalf )

    real,intent(in), dimension(:,:,:) :: uwnd,vwnd,atmp,zfull,zhalf,pfull,phalf
    real,intent(in), dimension(:,:,:) :: tausat
    real,intent(in), dimension(:,:)   :: taux,tauy,taulin,taun
    real,intent(out),dimension(:,:)   :: taup
    real,intent(out),dimension(:,:,:) :: dtaux,dtauy

    integer,dimension(size(pfull,1),size(pfull,2)) :: kref
    real,   dimension(size(pfull,1),size(pfull,2)) :: delp

    integer idim,jdim,kdim
    integer i,j,k,kr
    real :: dzhalf,zlast,rscale,phase,bfreq2,vvtau2
    real :: gfac,gfac1,dp,weight,wtsum

    real,parameter :: bfmin=.7e-2,bfmax=1.7e-2  ! min, max buoy. freq's (1/s)
    real,parameter :: vvmin=1.                  ! minimum surface wind (m/s)

    idim = SIZE(dtaux,1)
    jdim = SIZE(dtaux,2)
    kdim = SIZE(dtaux,3)

!   CALCULATE DECELERATION DUE TO LINEAR DRAG (~-rho^-1 dtau/dz)

    taup = 0.

    do k=1,kdim
      do j=1,jdim
        do i=1,idim
          dp = phalf(i,j,k+1) - phalf(i,j,k)
          gfac = (tausat(i,j,k+1) - tausat(i,j,k))
          taup(i,j) = taup(i,j) + gfac
          gfac1 = gfac*Grav/(dp*taulin(i,j))
          dtaux(i,j,k) = gfac1*taux(i,j)
          dtauy(i,j,k) = gfac1*tauy(i,j)
        enddo
      enddo
    enddo

!   find reference level (z ~ pi U/N)

    do j=1,jdim
      do i=1,idim
        kr = kdim
        phase = 0.
        zlast = zhalf(i,j,kdim)
        do while (phase <= Pi*zref_fac .and. kr > 1)
          kr = kr-1
          vvtau2 = (.5*((taux(i,j)*(uwnd(i,j,kr) + uwnd(i,j,kr-1)) &
                       + tauy(i,j)*(vwnd(i,j,kr) + vwnd(i,j,kr-1)))))**2 &
                 /max(1.e-8,(taux(i,j)*taux(i,j) + tauy(i,j)*tauy(i,j)))
          dzhalf = zfull(i,j,kr-1) - zfull(i,j,kr)
          bfreq2 = Grav*((atmp(i,j,kr-1) - atmp(i,j,kr))/dzhalf + Grav/Cp_Air) &
                    /(.5*(atmp(i,j,kr-1) + atmp(i,j,kr)))
          rscale = sqrt(max(bfmin*bfmin,min(bfmax*bfmax,bfreq2))/ &
                        max(vvmin*vvmin,vvtau2))
          dzhalf = zfull(i,j,kr-1) - zlast
          phase = phase + dzhalf*rscale
          zlast = zfull(i,j,kr-1)
        enddo
        kref(i,j) = kr
      enddo
    enddo

!   CALCULATE DECELERATION DUE TO NONLINEAR DRAG

    do j=1,jdim
      do i=1,idim
        kr=kref(i,j)
        dp = phalf(i,j,kdim+1) - phalf(i,j,kr)
        gfac = Grav/dp*taun(i,j)/taulin(i,j)
        wtsum = 0.
        do k=kr,kdim
          weight = pfull(i,j,k) - phalf(i,j,kr)
          wtsum = wtsum + weight
        enddo
        do k=kr,kdim
          weight = pfull(i,j,k) - phalf(i,j,kr)
          gfac1 = gfac*weight/wtsum
          dtaux(i,j,k) = dtaux(i,j,k) + gfac1*taux(i,j)
          dtauy(i,j,k) = dtauy(i,j,k) + gfac1*tauy(i,j)
        enddo
      enddo
    enddo

    return
  endsubroutine topo_drag_tend

  !=====================================================================

  subroutine topo_drag_init(lonb,latb,ierr)

    integer :: io,err,unit=19
    integer,intent(out) :: ierr
    integer :: i,j
    real,intent(in),dimension(:) :: lonb,latb

    if (module_is_initialized) return

    nlon = size(lonb(:))-1
    nlat = size(latb(:))-1

    !   Read namelist

    if (FILE_EXIST('input.nml')) then
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while (io /= 0)
          read (unit, nml = topo_drag_nml, iostat = io, end = 10) 
       enddo
10     continue
       call CLOSE_FILE (unit)
    endif
    call get_restart_io_mode(do_netcdf_restart)

    !   Output version details

    if ( mpp_pe() == mpp_root_pe() ) then
       call write_version_number(version, tagname)
       write(stdlog(), nml = topo_drag_nml) 
    endif

    !RSHif (gamma == beta) gamma = beta + 1.e-6
    if (gamma == beta + epsi) gamma = beta + epsi + 1.e-6


    allocate (flat(nlat))
    do j=1,nlat
       flat(j) = 1./sqrt(max(1., 2.*sin(.5*Radian*(latb(j) + latb(j+1)))))
    enddo

    module_is_initialized = .true.

    !   Read and interpolate mountain drag dataset
    if (file_exist('INPUT/topo_drag.res.nc')) then

       if(mpp_pe() == mpp_root_pe()) call error_mesg('topo_drag_mod', &
            'Reading NetCDF formatted restart file : INPUT/topo_drag.res.nc', NOTE)

       allocate (t11(nlon,nlat),t21(nlon,nlat),t12(nlon,nlat),t22(nlon,nlat))
       allocate (umin(nlon,nlat),vmin(nlon,nlat),umax(nlon,nlat),vmax(nlon,nlat))

       call read_data('INPUT/topo_drag.res.nc', 't11', t11, no_domain=.true.)
       call read_data('INPUT/topo_drag.res.nc', 't21', t21, no_domain=.true.)
       call read_data('INPUT/topo_drag.res.nc', 't12', t12, no_domain=.true.)
       call read_data('INPUT/topo_drag.res.nc', 't22', t22, no_domain=.true.)
       call read_data('INPUT/topo_drag.res.nc', 'umin', umin, no_domain=.true.)
       call read_data('INPUT/topo_drag.res.nc', 'vmin', vmin, no_domain=.true.)
       call read_data('INPUT/topo_drag.res.nc', 'umax', umax, no_domain=.true.)
       call read_data('INPUT/topo_drag.res.nc', 'vmax', vmax, no_domain=.true.)
    else       

       if (FILE_EXIST('INPUT/topo_drag.res')) then

          unit = OPEN_RESTART_FILE (file = 'INPUT/topo_drag.res', action = 'READ')

       else if (FILE_EXIST('INPUT/topo_drag')) then

          unit = OPEN_RESTART_FILE (file = 'INPUT/topo_drag', action = 'READ')

       else

          call ERROR_MESG ('topo_drag_init',  &
               'No sub-grid orography specified in topo_drag', &
               FATAL)
          do_restart_write = .false.
          return

       endif

       if(mpp_pe() == mpp_root_pe()) call error_mesg('topo_drag_mod', &
            'Reading native formatted restart file : INPUT/topo_drag.res', NOTE)

       allocate (t11(nlon,nlat),t21(nlon,nlat),t12(nlon,nlat),t22(nlon,nlat))
       allocate (umin(nlon,nlat),vmin(nlon,nlat),umax(nlon,nlat),vmax(nlon,nlat))

       call READ_DATA (unit,t11)
       call READ_DATA (unit,t21)
       call READ_DATA (unit,t12)
       call READ_DATA (unit,t22)
       call READ_DATA (unit,umin)
       call READ_DATA (unit,vmin)
       call READ_DATA (unit,umax)
       call READ_DATA (unit,vmax)

       call close_file (unit)

       return
    endif
  end subroutine topo_drag_init

  !=====================================================================

  subroutine get_pbl(is,js,uwnd,vwnd,atmp,zfull,zhalf,pfull,ktop)

    integer,intent(in) :: is,js
    integer,intent(out),dimension(:,:) :: ktop
    real,intent(in),dimension(:,:,:)   :: uwnd,vwnd,atmp,zfull,zhalf,pfull
    real,dimension(size(pfull,1),size(pfull,2))  :: ptop

    integer :: i,j,k,idim,jdim,kdim,jd,kbl
    real :: dzhalf,bfreq2

    real,parameter :: bf2 = 1.e-4

    idim = size(uwnd,1)
    jdim = size(uwnd,2)
    kdim = size(uwnd,3)

!   find highest model level in no-drag layer

    do j=1,jdim
      jd=j+js-1
      ptop(:,j) = (1. - no_drag_frac)*pfull(:,j,kdim)*flat(jd)
    enddo

    do k=kdim,1,-1 
      where (pfull(:,:,k) >= ptop(:,:)) 
        ktop(:,:) = k
      endwhere
    enddo

    ktop(:,:) = min(ktop(:,:),kdim-1)

    do j=1,jdim
      do i=1,idim
        kbl = ktop(i,j)
        dzhalf = zfull(i,j,kbl) - zfull(i,j,kbl+1)
        bfreq2 = Grav*((atmp(i,j,kbl) - atmp(i,j,kbl+1))/dzhalf + Grav/Cp_Air) & 
                 / (.5*(atmp(i,j,kbl) + atmp(i,j,kbl+1)))
        ptop(i,j) = pfull(i,j,kdim) &
                 - (pfull(i,j,kdim) - ptop(i,j))*2./(1. + bfreq2/bf2)
      enddo
    enddo

    do k=kdim,1,-1 
      where (pfull(:,:,k) >= ptop(:,:)) 
        ktop(:,:) = k
      endwhere
    enddo

    ktop(:,:) = min(ktop(:,:),kdim-1)

    return
  endsubroutine get_pbl

  !=====================================================================

  subroutine compute_pbl_top_index(is,js,zfull,z_pbl,ktop)

    integer,intent(in) :: is,js
    integer,intent(out),dimension(:,:) :: ktop
    real,intent(in),dimension(:,:,:)   :: zfull
    real,intent(in),dimension(:,:  )   :: z_pbl       

    integer :: i,j,k,idim,jdim,kdim,jd,kbl
    real :: dzhalf,bfreq2

    real,parameter :: bf2 = 1.e-4

    idim = size(zfull,1)
    jdim = size(zfull,2)
    kdim = size(zfull,3)

    do j=1,jdim
    do i=1,idim
    do k=kdim,1, -1
      if (zfull(i,j,k) > z_pbl(i,j)) then
        ktop(i,j) = k+1
        exit
      endif
    end do
    end do
    end do

    ktop(:,:) = MIN(ktop(:,:), kdim-1)

  end subroutine compute_pbl_top_index

  !=====================================================================

  subroutine topo_drag_end

    integer :: unit=99

    if (.not. do_restart_write) return

!   write out global arrays to restart file
    if (do_netcdf_restart) then

       if(mpp_pe() == mpp_root_pe()) call error_mesg('topo_drag_mod', &
            'Writing NetCDF formatted restart file : RESTART/topo_drag.res.nc', NOTE)
       call write_data('RESTART/topo_drag.res.nc', 't11', t11, no_domain=.true.)
       call write_data('RESTART/topo_drag.res.nc', 't21', t21, no_domain=.true.)
       call write_data('RESTART/topo_drag.res.nc', 't12', t12, no_domain=.true.)
       call write_data('RESTART/topo_drag.res.nc', 't22', t22, no_domain=.true.)
       call write_data('RESTART/topo_drag.res.nc', 'umin', umin, no_domain=.true.)
       call write_data('RESTART/topo_drag.res.nc', 'vmin', vmin, no_domain=.true.)
       call write_data('RESTART/topo_drag.res.nc', 'umax', umax, no_domain=.true.)
       call write_data('RESTART/topo_drag.res.nc', 'vmax', vmax, no_domain=.true.)
    else    
       
       if(mpp_pe() == mpp_root_pe()) call error_mesg('topo_drag_mod', &
            'Writing native formatted restart file : RESTART/topo_drag.res', NOTE)
       unit = OPEN_RESTART_FILE (file = 'RESTART/topo_drag.res', action = 'WRITE')
       
       call WRITE_DATA (unit,t11)
       call WRITE_DATA (unit,t21)
       call WRITE_DATA (unit,t12)
       call WRITE_DATA (unit,t22)
       call WRITE_DATA (unit,umin)
       call WRITE_DATA (unit,vmin)
       call WRITE_DATA (unit,umax)
       call WRITE_DATA (unit,vmax)

       call CLOSE_FILE (unit)
    endif
 
      module_is_initialized = .false.
 
    return
  end subroutine topo_drag_end

!#############################################################################      

endmodule topo_drag_mod
