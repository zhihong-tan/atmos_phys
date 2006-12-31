module atmos_age_tracer_mod
! <CONTACT EMAIL="William.Cooke@noaa.gov">
!   William Cooke
! </CONTACT>

! <REVIEWER EMAIL="Larry.Horowitz@noaa.gov">
!   Larry Horowitz
! </REVIEWER>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!     This code implements an age-of-air tracer.
! </OVERVIEW>

! <DESCRIPTION>
!     This code implements an age-of-air tracer, based on that in
!     the stratospheric chemistry code.
! </DESCRIPTION>

!-----------------------------------------------------------------------

use              fms_mod, only : file_exist, &
                                 write_version_number, &
                                 mpp_pe, &
                                 mpp_root_pe, &
                                 error_mesg, &
                                 FATAL,WARNING, NOTE, &
                                 stdlog
use     time_manager_mod, only : time_type
use     diag_manager_mod, only : send_data,            &
                                 register_static_field
use   tracer_manager_mod, only : get_tracer_index
use    field_manager_mod, only : MODEL_ATMOS


implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_age_tracer, atmos_age_tracer_init, atmos_age_tracer_end

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------
! namelist /atmos_age_tracer_nml/  


!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

!--- identification numbers for  diagnostic fields and axes ----

integer :: id_emiss

logical :: module_is_initialized=.FALSE.
real, parameter :: trop_age_cutoff = 0.1, trop_age_sq = trop_age_cutoff**2
real, parameter :: sec_per_day = 86400., &
                   age_relax_time = 10., & ! timescale for relaxation to zero in trop (days)
                   k_relax = 1./(age_relax_time*sec_per_day), & ! (1/sec)
                   days_per_year = 365.25, &
                   k_aging = 1./(days_per_year*sec_per_day) ! increase age at 1 yr/yr (convert to yr/sec)

!---- version number -----
character(len=128) :: version = ''
character(len=128) :: tagname = ''
!-----------------------------------------------------------------------

contains


!#######################################################################
!<SUBROUTINE NAME="atmos_age_tracer">
!<OVERVIEW>
! The routine that calculate the sources and sinks of age tracer.
!</OVERVIEW>
!<DESCRIPTION>
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_age_tracer (lon, lat, pwt, age, age_dt, Time, kbot)
!</TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     Longitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     Latitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     The pressure weighting array. = dP/grav
!   </IN>
!   <IN NAME="age" TYPE="real" DIM="(:,:,:)">
!     The array of the age tracer
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

!   <OUT NAME="age_dt" TYPE="real" DIM="(:,:,:)">
!     The array of the tendency of the age tracer
!   </OUT>
 subroutine atmos_age_tracer (lon, lat, pwt, age, age_dt, Time, kbot)

!-----------------------------------------------------------------------
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:,:) :: pwt, age
   real, intent(out), dimension(:,:,:) :: age_dt
   type(time_type), intent(in)         :: Time     
   integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(age,1),size(age,2),size(age,3)) ::  &
         source, sink
   integer :: i,j,k,id,jd,kd
   real :: dagesq(size(age,1))
!-----------------------------------------------------------------------

!<ERROR MSG="tropchem_driver_init must be called first." STATUS="FATAL">
!   Tropchem_driver_init needs to be called before tracer_driver.
!</ERROR>
   if (.not. module_is_initialized)  &
      call error_mesg ('atmos_age_tracer','atmos_age_tracer_init must be called first.', FATAL)

   id=size(age,1); jd=size(age,2); kd=size(age,3)

!----------- compute age tracer source and sink------------
!
!  Increase at the rate of 1 per second outside the troposphere.
!  Results expressed in years. Relax towards zero with a 10 
!  day timescale in the troposphere, denoted by DAGESQ less than 0.01 
!
   sink = 0.
   source = 0.
   do k = 1,kd
   do j = 1,jd
      dagesq(:) = (age(:,j,k) - age(:,j,kd))**2

      where (dagesq(:) < trop_age_sq) 
           sink(:,j,k) = -age(:,j,k)*k_relax
      elsewhere
           source(:,j,k) = k_aging
      end where

   end do
   end do


!------- tendency ------------------

   age_dt=source+sink
      

!-----------------------------------------------------------------------

 end subroutine atmos_age_tracer
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_age_tracer_init">
!<OVERVIEW>
! The constructor routine for the age tracer module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to initialize the age tracer module.
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_age_tracer_init (r, mask, axes, Time)
!</TEMPLATE>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace). 
!   </INOUT>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>
 subroutine atmos_age_tracer_init (r, axes, Time, nage, mask)

!-----------------------------------------------------------------------
!
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
real,             intent(inout), dimension(:,:,:,:) :: r
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
integer,          intent(out)                       :: nage
real, intent(in), dimension(:,:,:), optional        :: mask

!
!-----------------------------------------------------------------------
!
      integer :: n

      if (module_is_initialized) return

!---- write namelist ------------------

      call write_version_number (version, tagname)
!     if ( mpp_pe() == mpp_root_pe() ) &
!       write ( stdlog(), nml=atmos_age_tracer_nml )
 
      nage = -1
      n = get_tracer_index(MODEL_ATMOS,'AGE' )
      if (n>0) nage = n

      module_is_initialized = .TRUE.

!-----------------------------------------------------------------------

 end subroutine atmos_age_tracer_init
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_age_tracer_end">
!<OVERVIEW>
!  The destructor routine for the age tracer module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_age_tracer_end
!</TEMPLATE>
 subroutine atmos_age_tracer_end
 
      module_is_initialized = .FALSE.

 end subroutine atmos_age_tracer_end
!</SUBROUTINE>


end module atmos_age_tracer_mod



