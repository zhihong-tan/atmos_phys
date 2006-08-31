module atmos_dust_mod
! <DESCRIPTION>
!   This module evaluates the change of mass mixing ratio for mineral dust
!   particles due to their emission from preferential sources, and the removal
!   by gravitational settling. The dust particles are transported as dry
!   particles. No hygroscopic growth is considered.
!   The size distribution of sea salt ranges from 0.1 to 10 um (dry radius)
!   and is divided into 5 bins. For each bin, the volume size distribution
!   dV/dlnr is considered constant.
! </DESCRIPTION>
! <CONTACT EMAIL="Paul.Ginouxe@noaa.gov">
!   Paul Ginoux
! </CONTACT>
!-----------------------------------------------------------------------

use              fms_mod, only : file_exist, &
                                 write_version_number, &
                                 mpp_pe, &
                                 mpp_root_pE, &
                                 close_file,           &
                                 stdlog
use     time_manager_mod, only : time_type
use     diag_manager_mod, only : send_data,            &
                                 register_diag_field
use   tracer_manager_mod, only : get_tracer_index, &
                                 set_tracer_atts
use    field_manager_mod, only : MODEL_ATMOS
use atmos_tracer_utilities_mod, only : wet_deposition,       &
                                 dry_deposition, &
                                  interp_emiss
use     constants_mod, only : PI, GRAV, RDGAS, DENS_H2O, PSTD_MKS, WTMAIR


implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_dust_sourcesink, atmos_dust_init, atmos_dust_end

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------

!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

integer :: ndust=0  ! tracer number for dust
!--- identification numbers for  diagnostic fields and axes ----

integer :: id_DU_emis(5), id_DU_setl(5)
integer :: id_DU_source

!--- Arrays to help calculate tracer sources/sinks ---
real, allocatable, dimension(:,:) :: DU_source

logical :: module_is_initialized=.FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id: atmos_dust.F90,v 13.0 2006/03/28 21:15:27 fms Exp $'
character(len=128) :: tagname = '$Name: memphis_2006_08 $'
!-----------------------------------------------------------------------

contains


!#######################################################################
!<SUBROUTINE NAME="atmos_dust_sourcesink">
!<OVERVIEW>
! The routine that calculate the sources and sinks of dust.
!</OVERVIEW>
 subroutine atmos_dust_sourcesink (i_DU,ra,rb,dustref,dustden, &
       lon, lat, frac_land, pwt, &
       zhalf, pfull, w10m, t, rh, &
       dust, dust_dt, Time, is,ie,js,je,kbot)

!-----------------------------------------------------------------------
   integer, intent(in)                 :: i_DU
   real, intent(in)                    :: ra
   real, intent(in)                    :: rb
   real, intent(in)                    :: dustref
   real, intent(in)                    :: dustden
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:)   :: frac_land
   real, intent(in),  dimension(:,:)   :: w10m
   real, intent(in),  dimension(:,:,:) :: pwt, dust
   real, intent(in),  dimension(:,:,:) :: zhalf, pfull, t, rh
   real, intent(out), dimension(:,:,:) :: dust_dt
   type(time_type), intent(in) :: Time     
   integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
 real, dimension(size(dust,1),size(dust,2)) :: DU_setl, DU_emis
 logical :: flag
integer  i, j, k, m, id, jd, kd, kb, ir
integer, intent(in)                    :: is, ie, js, je
!----------------------------------------------
!     Dust parameters
!----------------------------------------------
      real, dimension(5) ::   frac_s
      real :: u_ts0, diam, den

      real, dimension(size(dust,3)) :: setl
      real, dimension(size(dust,1),size(dust,2)) :: u_ts, source, maxgw

      real       ::  CH           ! THE tuning factor for dust emission

      real, parameter :: small_value = 1.e-20
      real, parameter :: mtcm = 100.            ! meter to cm
      real, parameter :: mtv  = 1. ! factor conversion for mixing ratio of dust
      real, parameter :: ptmb = 0.01     ! pascal to mb

      real :: rhb, rcm, g0
      real :: ratio_r, rho_wet_dust, viscosity, free_path, C_c, vdep
      real :: rho_air
      real :: rmid, rwet
!-----------------------------------
!    SET-Up  DATA
!-----------------------------------

      data frac_s/0.1,0.225,0.225,0.225,0.225/
      data CH/0.5e-9/

!-----------------------------------------------------------------------

      id=size(dust,1); jd=size(dust,2); kd=size(dust,3)

!----------- dust sources on local grid
      do i=1,id
        do j=1,jd
          source(i,j)=DU_source(i+is-1,j+js-1)
        enddo
      enddo

!----------- compute dust emission ------------
      DU_emis(:,:)   = 0.0
      DU_setl(:,:)   = 0.0
      dust_dt(:,:,:) = 0.0

      den=dustden*1.e-3
      diam=2.*dustref*1.e2
      g0=GRAV*1.e2

!       rho_air = pfull(i,j,kd)/t(i,j,kd)/RDGAS   ! Air density [kg/m3]
!        rho_air=1.25
!        u_ts=0.13*1.e-2*sqrt(den*g0*diam/(rho_air*1.e-3))* &
!              sqrt(1.+0.006/den/g0/(diam)**2.5)/ &
!              sqrt(1.928*(1331*(diam)**1.56+0.38)**0.092-1) 
        where ( frac_land.gt.0.1 )
          u_ts=0.
          DU_emis = CH * frac_s(i_DU)*source * frac_land &
             * w10m**2 * (w10m - u_ts)
        endwhere
        dust_dt(:,:,kd)=dust_dt(:,:,kd)+DU_emis(:,:)/pwt(:,:,kd)*mtv

! Send the emission data to the diag_manager for output.
      if (id_DU_emis(i_DU) > 0 ) then
        used = send_data ( id_DU_emis(i_DU), DU_emis, Time, &
              is_in=is,js_in=js )
      endif

         rcm=dustref*mtcm            ! Particles radius in centimeters
!------------------------------------------
!       Solve at the model TOP (layer plev-10)
!------------------------------------------
      do j=1,jd
        do i=1,id
          setl(:)=0.
          if (present(kbot)) then
              kb=kbot(i,j)
          else
             kb=kd
          endif
          do k=1,kb
              rhb=amin1(0.99,rh(i,j,k))
              rhb=amax1(0.01,rhb)
!----------------------------------------------------------
!     Aerosol growth with relative humidity
!----------------------------------------------------------

            rwet=dustref  ! Add any particle growth here
            ratio_r=(dustref/rwet)**3.   ! Ratio dry over wet radius cubic power
            rho_wet_dust=ratio_r*dustden+(1.-ratio_r)*DENS_H2O     ! Density of wet aerosol [kg/m3]
            viscosity = 1.458E-6 * t(i,j,k)**1.5/(t(i,j,k)+110.4)     ! Dynamic viscosity
            free_path=6.6e-8*t(i,j,k)/293.15*(PSTD_MKS/pfull(i,j,k))
            C_c=1. + free_path/dustref* &          ! Slip correction [none]
                  (1.257+0.4*exp(-1.1*dustref/free_path))
            Vdep=2./9.*C_c*GRAV*rho_wet_dust*rwet**2./viscosity   ! Settling velocity [m/s]
            rho_air = pfull(i,j,k)/t(i,j,k)/RDGAS      ! Air density [kg/m3]
            if (dust(i,j,k).gt.0.) then
              setl(k)=dust(i,j,k)*rho_air/mtv*vdep    ! settling flux [kg/m2/s]
            endif
          enddo
          dust_dt(i,j,1)=dust_dt(i,j,1)-setl(1)/pwt(i,j,1)*mtv
          dust_dt(i,j,2:kb)=dust_dt(i,j,2:kb) &
             + ( setl(1:kb-1) - setl(2:kb) )/pwt(i,j,2:kb)*mtv
          DU_setl(i,j)=setl(kb)
        enddo
      enddo 

! Send the settling data to the diag_manager for output.
      if (id_DU_setl(i_DU) > 0 ) then
        used = send_data ( id_DU_setl(i_DU), DU_setl, Time, &
              is_in=is,js_in=js )
      endif


!-----------------------------------------------------------------------

 end subroutine atmos_dust_sourcesink
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_dust_init">
!<OVERVIEW>
! The constructor routine for the dust module.
!</OVERVIEW>
 subroutine atmos_dust_init (lonb, latb, axes, Time, mask)
!-----------------------------------------------------------------------
real, intent(in),    dimension(:)               :: lonb, latb
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
      character(len=16) ::  fld
      character*1 :: numb(5)
      data numb/'1','2','3','4','5'/


      if (module_is_initialized) return

!---- write namelist ------------------

      call write_version_number (version, tagname)

!----- set initial value of dust ------------
    do m=1,5

       n = get_tracer_index(MODEL_ATMOS,'dust'//numb(m))
       if (n>0) then
         ndust=n
       call set_tracer_atts(MODEL_ATMOS,'dust'//numb(m),'dust'//numb(m),'mmr')
         if (ndust > 0 .and. mpp_pe() == mpp_root_pe()) &
                write (*,30) 'dust'//numb(m),ndust
         if (ndust > 0 .and. mpp_pe() == mpp_root_pe()) &
                write (stdlog(),30) 'dust '//numb(m),ndust
       endif


  30        format (A,' was initialized as tracer number ',i2)
! Register a diagnostic field : emission of dust
     id_DU_emis(m) = register_diag_field ( mod_name,            &
                     'dust'//numb(m)//'_emis', axes(1:2),Time,  &
                     'dust'//numb(m)//'_emis', 'kg/m2/s',       &
                     missing_value=-999.  )

! Register a diagnostic field : settling of dust
     id_DU_setl(m) = register_diag_field ( mod_name,            &
                     'dust'//numb(m)//'_setl', axes(1:2),Time,  &
                     'dust'//numb(m)//'_setl', 'kg/m2/s',       &
                     missing_value=-999.  )
enddo
!
     if (.not.allocated(DU_source)) allocate (DU_source(size(lonb)-1,size(latb)-1))
     id_DU_source  = register_diag_field ( mod_name,             &
                      'DU_source',axes(1:2),Time,                    &
                      'Dust_source', 'none')

     call dust_source_input(lonb, latb, Time)

     call write_version_number (version, tagname)

      module_is_initialized = .TRUE.


!-----------------------------------------------------------------------

 end subroutine atmos_dust_init
!</SUBROUTINE>

!#######################################################################
!<SUBROUTINE NAME="atmos_dust_end">
!<OVERVIEW>
!  The destructor routine for the dust module.
!</OVERVIEW>
 subroutine atmos_dust_end

      if (allocated(DU_source)) deallocate(DU_source) 
      module_is_initialized = .FALSE.

 end subroutine atmos_dust_end
!</SUBROUTINE>
!#######################################################################
!<SUBROUTINE NAME="dust_source_input">
!<OVERVIEW>
!  This subroutine read the fraction of dust source at each grid cell
!  based on Ginoux et al. (JGR, 2001). This source function is
!  considered invariant in time (no change in topography and vegetation).
!</OVERVIEW>
 subroutine dust_source_input(lonb, latb, Time)
real, dimension(:),    intent(in) :: lonb, latb
type(time_type),intent(in) :: Time

integer      :: i, j, unit, io, ios
real         :: dtr, lat_S, lon_W, dlon, dlat
real         :: DU_source1(360,180)
logical :: opened
!
!
      dtr= PI/180.
      lat_S= -89.5*dtr; lon_W= 0.5*dtr
      dlon = 1.*dtr; dlat = 1.*dtr;

do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         DU_source1 = 0.0

         open (unit,file='INPUT/dust_source_1x1.txt', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
           read  (unit,FMT=2000, end=11) ((DU_source1(i,j),i=1,360),j=1,180)
         else
           write(*,*) '***ERROR: Opening dust source file'
         endif
  11    call close_file (unit)
 2000 format(10f8.5)
! Interpolate the R30 emission field to the resolution of the model.
        call interp_emiss ( DU_source1, lon_W, lat_S, dlon, dlat, &
                     DU_source)


if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading Dust source'
endif

! Send the dust source data to the diag_manager for output.
         if (id_DU_source > 0 ) &
           used = send_data ( id_DU_source, DU_source , Time )

end subroutine dust_source_input
!</SUBROUTINE>

end module atmos_dust_mod
