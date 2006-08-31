module atmos_soa_mod
! <DESCRIPTION>
!   This module is an implementation of Secondary organic aerosols (SOA)
!   from anthropogenic activities, and is based on Tie et al. (JGR, 2003).
!   The only souce of SOA is due to the oxydation of C4H10 by OH.
!   The concentrations of these 2 gas species are read as input.
! </DESCRIPTION>
! <WARNING>
!  To save space only the actual month of input files are kept in memory. 
!  This implies that the "atmos_SOA_init" should be executed at the begining 
!  of each month. In other words, the script should not run more than 1 month
!  without a restart.
! </WARNING>
! <CONTACT EMAIL="Paul.Ginouxe@noaa.gov">
!   Paul Ginoux
! </CONTACT>
!-----------------------------------------------------------------------

use                    fms_mod, only : file_exist,              &
                                       write_version_number,    &
                                       mpp_pe,                  &
                                       mpp_root_pE,             &
                                       close_file,              &
                                       stdlog
use           time_manager_mod, only : time_type
use           diag_manager_mod, only : send_data,               &
                                       register_diag_field,     &
                                       register_static_field
use         tracer_manager_mod, only : get_tracer_index,        &
                                       set_tracer_atts
use          field_manager_mod, only : MODEL_ATMOS
use              constants_mod, only : PI, GRAV, RDGAS, WTMAIR
use atmos_tracer_utilities_mod, only : interp_emiss

implicit none

private
!-----------------------------------------------------------------------
!----- interfaces -------
!
public  atmos_SOA_init, atmos_SOA_end, SOA_source_input, chem_SOA

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------

!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

integer :: nSOA = 0  ! tracer number for Secondary Organic Aerosol 

!--- identification numbers for  diagnostic fields and axes ----
integer ::   id_OH_conc            = 0
integer ::   id_C4H10_conc         = 0
integer ::   id_SOA_chem           = 0

real, allocatable, dimension(:,:,:) :: OH_conc
real, allocatable, dimension(:,:,:) :: C4H10_conc

logical :: module_is_initialized=.FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id: atmos_soa.F90,v 13.0 2006/03/28 21:15:37 fms Exp $'
character(len=128) :: tagname = '$Name: memphis_2006_08 $'
!-----------------------------------------------------------------------

contains


!#######################################################################

!<SUBROUTINE NAME="atmos_SOA_init">
!<OVERVIEW>
! The constructor routine for the soa module.
!</OVERVIEW>
 subroutine atmos_SOA_init ( lonb, latb, nlev, axes, Time, mask)
!-----------------------------------------------------------------------
real,             intent(in), dimension(:)          :: lonb, latb
integer,          intent(in)                        :: nlev
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
real, intent(in), dimension(:,:,:), optional        :: mask
character(len=7), parameter :: mod_name = 'tracers'
logical :: flag
integer :: n, m
!
!-----------------------------------------------------------------------
!
      integer  log_unit,unit,io,index,ntr,nt
      character*3 :: SOA_tracer
!
!     1. C4H10     = nButane
!                                                                      
      data SOA_tracer/     'SOA'/
      
      if (module_is_initialized) return

!---- write namelist ------------------

      call write_version_number (version, tagname)

!----- set initial value of soa ------------

      nSOA = get_tracer_index(MODEL_ATMOS,'SOA')
      if (nSOA > 0) then
         call set_tracer_atts(MODEL_ATMOS,'SOA','SOA','mmr')
         if (nSOA > 0 .and. mpp_pe() == mpp_root_pe()) &
                 write (*,30) SOA_tracer,nsoa
         if (nSOA > 0 .and. mpp_pe() == mpp_root_pe()) &
                 write (stdlog(),30) SOA_tracer,nsoa
      endif


  30   format (A,' was initialized as tracer number ',i2)
      if (.not.allocated(OH_conc)) then
        allocate (OH_conc(size(lonb)-1,size(latb)-1,nlev)  )
      endif
      if (id_OH_conc .eq. 0 ) &
        id_OH_conc    = register_diag_field ( mod_name,           &
                      'OH_SOA_conc',axes(1:3),Time,                        &
                      'Hydroxyl radical concentration',           &
                      'molec.cm-3')

      if (.not.allocated(C4H10_conc)) & 
        allocate (C4H10_conc(size(lonb)-1,size(latb)-1,nlev)  )

      id_C4H10_conc    = register_diag_field ( mod_name,           &
                      'C4H10_mmr',axes(1:3),Time,                        &
                      'nButane concentration',           &
                      'mmr')

      id_SOA_chem    = register_diag_field ( mod_name,       &
                      'SOA_chem',axes(1:3),Time,            &
                      'SOA production by C4H10 + OH',        &
                      'kgC/m2/s')


      call write_version_number (version, tagname)

      module_is_initialized = .TRUE.

!-----------------------------------------------------------------------
 end subroutine atmos_SOA_init
!</SUBROUTINE>

!#######################################################################
!<SUBROUTINE NAME="atmos_SOA_end">
!<OVERVIEW>
!  The destructor routine for the soa module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_SOx_end
!</TEMPLATE>
 subroutine atmos_SOA_end

      if (allocated(OH_conc)) deallocate(OH_conc)
      if (allocated(C4H10_conc)) deallocate(C4H10_conc)

      module_is_initialized = .FALSE.

 end subroutine atmos_SOA_end
!</SUBROUTINE>
!#######################################################################
!<SUBROUTINE NAME="SOA_source_input">
!<OVERVIEW>
!  This subroutine read the monthly mean concentrations of OH and C4H10
!  *****WARNING:
!  To save space only the actual month is kept in memory which implies
!  that the "atmos_SOA_init" should be executed at the begining of each
!  month. In other words, the script should not run more than 1 month
!  without a restart
!</OVERVIEW>

 subroutine SOA_source_input( lon, lat, imonth, Time, is, ie, js, je, kbot)
!--- Input variables
        integer, intent(in)              :: imonth
        real, dimension(:,:), intent(in) :: lon, lat
        type(time_type),intent(in)       :: Time
        integer, intent(in), dimension(:,:), optional :: kbot
        integer, intent(in)              :: is, ie, js, je
!--- Working variables
        real                        :: dtr, lat_S, lon_W, dlon, dlat
        real                        :: variable
        integer                     :: i, j, l, im, unit, ios
        logical                     :: opened
        real, dimension(144, 90)    :: data2D
        real, dimension(144, 90,24) :: data3D
        character (len=3)           :: month(12)
!--- Molecular weight [g/mole]
        real, parameter             :: wtm_C   = 12.
        real, parameter             :: wtm_C4H10 = 64.

!--- Input filenames
        character (len=7 ) :: FNMOH     = 'OH_AM2_'
        character (len=10) :: FNMC4H10  = 'C4H10_AM2_'

        data month/'jan','feb','mar','apr','may','jun','jul', &
                   'aug','sep','oct','nov','dec'/
!---
        dtr = PI/180.
        lat_S= -90.*dtr; lon_W= -180.*dtr
        dlon = 2.5*dtr; dlat = 2.*dtr;
!
!------------------------------------------------------------------------
! Read 12 monthly C4H10 concentration
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) write(*,*) 'Reading C4H10 concentration'
!
        c4h10_conc(:,:,:) = 0.  ![molec/cm3]
        data3D(:,:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo

        open (unit,file='INPUT/'//FNMC4H10//month(imonth)//'.txt', &
              form='formatted',action='read')
        do
          read (unit, '(3i4,e12.4)',iostat=ios) i,j,l,variable
          data3D(i,j,l)=variable
          if ( ios.ne.0 ) exit
        enddo
        call close_file (unit)
      
!------------------------------------------------------------------------
! End reading C4H10 concentration
!------------------------------------------------------------------------
! --- Interpolate data
      do l=1,24
        call interp_emiss ( data3D(:,:,l), lon_W, lat_S, dlon, dlat, &
                            C4H10_conc(:,:,l))
      enddo
! Send the C4H10 data to the diag_manager for output.
         if (id_C4H10_conc > 0 ) &
           used = send_data ( id_C4H10_conc, &
                C4H10_conc, Time, is_in=is, js_in=js )
!------------------------------------------------------------------------
! Read 12 monthly OH concentration
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) write(*,*) 'Reading OH concentration'
!
        OH_conc(:,:,:) = 0.  ![molec/cm3]
        data3D(:,:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo

        open (unit,file='INPUT/'//FNMOH//month(imonth)//'.txt', &
              form='formatted',action='read')
        do
          read (unit, '(3i4,e12.4)',iostat=ios) i,j,l,variable
          data3D(i,j,l)=variable
          if ( ios.ne.0 ) exit
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading OH concentration
!------------------------------------------------------------------------
! --- Interpolate data
      do l=1,24
        call interp_emiss ( data3D(:,:,l), lon_W, lat_S, dlon, dlat, &
                            OH_conc(:,:,l))
      enddo
! Send the OH data to the diag_manager for output.
         if (id_OH_conc > 0 ) &
           used = send_data ( id_OH_conc, &
                OH_conc, Time, is_in=is, js_in=js , ks_in=1 )

end subroutine SOA_source_input
!</SUBROUTINE>
!-----------------------------------------------------------------------
      SUBROUTINE CHEM_SOA(pwt,temp,pfull, dt, &
                          jday,hour,minute,second,lat,lon, &
                          SOA, SOA_dt, Time,is,ie,js,je,kbot)

! ****************************************************************************
      real, intent(in),    dimension(:,:,:)          :: pwt
      real, intent(in),    dimension(:,:,:)          :: temp,pfull
      real, intent(in)                               :: dt
      integer, intent(in)                            :: jday, hour,minute,second
      real, intent(in),  dimension(:,:)              :: lat, lon  ! [radian]
      real, intent(in),    dimension(:,:,:)          :: SOA
      real, intent(out),   dimension(:,:,:)          :: SOA_dt
      type(time_type), intent(in)                    :: Time
      integer, intent(in),  dimension(:,:), optional :: kbot
      integer, intent(in)                            :: is,ie,js,je
! Working vectors
      real, dimension(size(SOA,1),size(SOA,2),size(SOA,3)) :: SOA_chem
      real, dimension(size(SOA,1),size(SOA,2)) :: &
               xu, dayl, h, hl, hc, hred, fac_OH, fact_OH
      real                                       :: oh, c4h10
      real, parameter                            :: wtm_C = 12.
      real, parameter                            :: yield = 0.1
      real, parameter                            :: small_value=1.e-21
      real, parameter                            :: A0 = 0.006918
      real, parameter                            :: A1 = 0.399912
      real, parameter                            :: A2 = 0.006758
      real, parameter                            :: A3 = 0.002697
      real, parameter                            :: B1 = 0.070257
      real, parameter                            :: B2 = 0.000907
      real, parameter                            :: B3 = 0.000148
      real                                       :: decl, hd, x
      integer :: i,j,k,id,jd,kd
      integer                                    :: istep, nstep
! Local grid sizes
      id=size(SOA,1); jd=size(SOA,2); kd=size(SOA,3)

      x = 2. *pi *float(jday-1)/365.
      decl = A0 - A1*cos(  X) + B1*sin(  X) - A2*cos(2.*X) + B2*sin(2.*X) &
           - A3*cos(3.*X) + B3*sin(3.*X)
      xu(:,:) = -tan(lat(:,:))*tan(decl)
      where ( xu > -1 .and. xu < 1 ) dayl=acos(xu)/pi
      where ( xu <= -1 ) dayl = 1.
      where ( xu >= 1 ) dayl = 0.
!   Calculate normalization factors for OH and NO3 such that
!   the diurnal average respect the monthly input values.
      hd=0.
      fact_OH(:,:)  = 0.
      nstep = int(24.*3600./dt)
      do istep=1,nstep
        hd=hd+dt/3600./24.
        hl(:,:) = pi*(1.-dayl(:,:))
        hc(:,:) = pi*(1.+dayl(:,:))
        h(:,:)=2.*pi*mod(hd+lon(:,:)/2./pi,1.)
        where ( h.ge.hl .and. h.lt.hc )
! Daytime
          hred=(h-hl)/(hc-hl)
          fact_OH  = fact_OH + amax1(0.,sin(pi*hred)/2.)/nstep
        endwhere
      enddo


      hd=amax1(0.,amin1(1.,(hour+minute/60.+second/3600.)/24.))
      hl(:,:) = pi*(1.-dayl(:,:))
      hc(:,:) = pi*(1.+dayl(:,:))
      h(:,:)=2.*pi*mod(hd+lon(:,:)/2./pi,1.)
      fac_OH(:,:)  = 0.
      where ( h.ge.hl .and. h.lt.hc )
! Daytime
          hred=(h-hl)/(hc-hl)
          fac_OH  = amax1(0.,sin(pi*hred)/2.)/fact_OH
      elsewhere
! Nightime
          fac_OH  = 0.
      endwhere

      do i=1,id
        do j=1,jd
          do k=1,kd 
            c4h10=max(small_value,C4H10_conc(i+is-1,j+js-1,k))
            oh   =max(small_value,OH_conc(i+is-1,j+js-1,k))*fac_oh(i,j)
            SOA_dt(i,j,k) = 1.55E-11 * exp( -540./temp(i,j,k) ) *yield &
                * c4h10 * oh
          enddo
        enddo
      enddo

      SOA_chem(:,:,:)=SOA_dt(:,:,:)*pwt(:,:,:)*wtm_C/WTMAIR

      if (id_SOA_chem > 0) then
        used = send_data ( id_SOA_chem, &
              SOA_chem, Time,is_in=is,js_in=js,ks_in=1)
      endif


end subroutine chem_soa


end module atmos_soa_mod
