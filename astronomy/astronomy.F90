module astronomy_mod

!----------------------------------------------------------------------

use utilities_mod,    only: error_mesg,                             &
                            open_file, file_exist, check_nml_error, &
                            close_file, FATAL, get_my_pe

use time_manager_mod, only: time_type, set_time, get_time, set_date, &
                            operator(-), operator(//), operator(<),  &
                            length_of_year

!----------------------------------------------------------------------

implicit none
private

public :: set_orbital_parameters, &
          get_orbital_parameters, &
          set_ref_date_of_ae,     &
          get_ref_date_of_ae,     &
          diurnal_solar,          &
          daily_mean_solar,       &
          annual_mean_solar,      &
          astronomy_init,         &
          set_period,             &
          get_period

character(len=128) :: version = '$Id: astronomy.F90,v 1.2 2000/08/04 18:05:13 fms Exp $'
character(len=128) :: tag = '$Name: calgary $'

! description of public interfaces

!----------------------------------------------------------------
! CALL  SET_ORBITAL_PARAMETERS(ECC, OBLIQ, PER)
!   real, intent(in) :: ecc, obliq, per
!
!   ecc = eccentricity of orbital ellipse (non-dimensional)
!   obliq = obliquity (DEGREES)
!   per = longitude of perihelion with respect to autumnal equinox 
!     in northern hemisphere (DEGREES)
!
!  default values are 
!                 ecc   = 0.01671
!                 per   = 102.932  (perihelion near winter soltice in NH)
!                 obliq = 23.439
!
!----------------------------------------------------------------
! CALL  GET_ORBITAL_PARAMETERS(ECC, OBLIQ, PER)
!   real, intent(out) :: ecc, obliq, per
!
!----------------------------------------------------------------
! CALL  SET_PERIOD(PERIOD)
!
! integer, intent(in) :: period
!    or
! type(time_type), intent(in) :: period
!
!    default value is obtained from the time
!    manager via function length_of_year
!
!----------------------------------------------------------------
! CALL  GET_PERIOD(PERIOD)
!
! integer, intent(out) :: period
!    or
! type(time_type), intent(out) :: period
!
!----------------------------------------------------------------
! CALL SET_REF_DATE_OF_AE(DAY, MONTH, YEAR) 
! integer, intent(in) :: day, month, year  
!    or
! CALL SET_REF_DATE_OF_AE(DAY, MONTH, YEAR, SECOND, MINUTE, HOUR)
! integer, intent(in) :: day, month, year, second, minute, hour
!
! only used if calls are made to the calandar versions of the routines
!     diurnal_solar and   daily_mean_solar
!
!  sets the calandar date of autumnal equinox in the Northern hemisphere
!   for a particular year
!  if the NO_LEAP calandar is used, then the date of autumnal equinox will
!    be the same every year
!  if JULIAN is used, then the date of autumnal equinox will return 
!    to the same value every 4th year
!
! defaults :: day    = 23
!             month  = 9
!             year   = 1998
!             hour   = 5
!             minute = 37
!             second = 0
! 
!----------------------------------------------------------------
! CALL GET_REF_DATE_OF_AE(DAY, MONTH, YEAR, SECOND, MINUTE, HOUR)
! integer, intent(out) :: day, month, year, second, minute, hour
!
!----------------------------------------------------------------
! CALL DIURNAL_SOLAR(COSZ, SOLAR, LAT, LON, TIME, DT_TIME)
!   or 
! CALL DIURNAL_SOLAR(COSZ, SOLAR, LAT, LON, GMT, TIME_SINCE_AE, DT)
!
!the first option is for those using time_manager_mod for time-keeping
!
! the output are the cosine of the zenith angle and the normalized solar
! flux, solar = cosz*((a/r)**2), where 
! r is the earth-sun distance and a is the semi-major axis of the orbital 
!  ellipse
! cosz and solar are both set to 0.0 in polar night
!  (should we output (a/r)**2 or r instead ??)
!  these can be 2d, 1d, or 0d fields:
!    real, intent(in), dimension(:,:) :: cosz, solar
! or real, intent(in), dimension(:)   :: cosz, solar
! or real, intent(in)                 :: cosz, solar
!
! lat and lon are the input latitudes and longitudes (in RADIANS) at
! which the solar flux or zenith angle are required, they may be either
! 2d, 1d, or 0d (consistent with the desired output)
!    real, intent(in), dimension(:,:) :: lat, lon
! or real, intent(in), dimension(:)   :: lat, lon
! or real, intent(in)                 :: lat, lon
!
! if NOT using the calendars defined in time_manager_mod, the time of day
! is set by 
!    real, intent(in) :: gmt
! gmt is measured in RADIANS (1 DAY = 2* pi) 
!          and = 0.0 at midnight at longitude = 0.0
! the time of year is set by 
!    real, intent(in) :: time_since_ae
! and is also measured in radians (1 ORBITAL YEAR = 2*pi)
!          and = 0.0 at autumnal equinox
! One also has the option of averaging the flux and cos(z) over the 
! time interval between gmt and gmt + dt by including the argument
!   real, intent(in), optional :: dt
! This average is computed analytically, and should be exact
!   except for the fact that changes in earth-sun distance over
!   the time interval dt are ignored
! NOTE:  dt is also measured in radians and must be less than pi = 1/2 day
! In the context of a diurnal GCM, this option should always be employed 
!   to insure that the total flux at the top of the atmosphere is not 
!   modified by time truncation error
!
! if using the calendars defined in time_manager_mod, then the time 
! is simply set by 
!    type(time_type), intent(in) :: time
! and the optional time interval for averaging is also of this type 
!    type(time_type), optional :: dt_time
!  (see test.90 for examples of the use of these types)
!-------------------------------------------------------------------------
! CALL DAILY_MEAN_SOLAR(COSZ, SOLAR, LAT, TIME)
!   or 
! CALL DAILY_MEAN_SOLAR(COSZ, SOLAR, LAT, TIME_SINCE_AE)
!
! as above, but only latitude is input and there is no option for 
!    time averaging
! input and output can be 2d, 1d, or 0d
! the relation between cosz and solar is now
! solar = cosz*((a/r)**2)*h/pi
!   where h is the "half-day-length" in radians 
!  (so h/pi is the fraction of the day in which the sun is shining )
!
!  Therefore, cosz is an effective zenith angle that would produce 
!  the correct flux if the sun where fixed at that position for the 
!  period of daylight.  
!
! Should one also, or instead, compute cosz weighted by the 
!   instantaneous flux, averaged over the day ??
!-------------------------------------------------------------------------
! CALL ANNUAL_MEAN_SOLAR(COSZ, SOLAR, LAT)
!  real, intent(in), dimension(:) :: lat
!  real, intent(out), dimension(:) :: cosz, solar
!
!   annual mean values are obtained by averaging output from daily_mean_solar
!   over the year
!   in this average, the daily mean effective cosz is weighted by the 
!    daily mean solar flux
!-------------------------------------------------------------------------

interface diurnal_solar
   module procedure diurnal_solar_2d
   module procedure diurnal_solar_1d
   module procedure diurnal_solar_0d
   module procedure diurnal_solar_cal_2d
   module procedure diurnal_solar_cal_1d
   module procedure diurnal_solar_cal_0d
end interface

interface daily_mean_solar
   module procedure daily_mean_solar_2d
   module procedure daily_mean_solar_1d
   module procedure daily_mean_solar_0d
   module procedure daily_mean_solar_cal_2d
   module procedure daily_mean_solar_cal_1d
   module procedure daily_mean_solar_cal_0d
end interface

interface annual_mean_solar
   module procedure annual_mean_solar_2d
   module procedure annual_mean_solar_1d
end interface

interface get_period
   module procedure get_period_time_type, get_period_integer
end interface

interface set_period
   module procedure set_period_time_type, set_period_integer
end interface

real :: ecc   = 0.01671
real :: per   = 102.932
real :: obliq = 23.439

integer :: day_ae    = 23
integer :: month_ae  = 9
integer :: year_ae   = 1998
integer :: hour_ae   = 5
integer :: minute_ae = 37
integer :: second_ae = 0

! time_type is defined in time_manager_mod

type(time_type) :: autumnal_eq_ref

!------------------------------------------------------------------
! period_time_type = period of one absolute orbit.
!                    default value is obtained from the time
!                    manager via function length_of_year.

type(time_type) :: period_time_type

! To override default value of period_time_type it must be
! specified in the namelist as an integer number of seconds.

integer :: period=0
!------------------------------------------------------------------

!   num_angles = number of intervals into which the year is divided
!     for the purpose of computing orbital positions
!   orb_angle = table of orbital positions (0 to 2*pi) from which
!    position is determined by interpolation in time
!   init = logical flag set to .true. after module is initialized

integer, parameter :: num_angles = 3600
real               :: orb_angle(0:num_angles)
logical            :: init = .false.
real               :: pi, twopi, deg_to_rad

namelist /astronomy_nml/ ecc, obliq, per, period, &
                         year_ae, month_ae,  day_ae,         & 
                         hour_ae, minute_ae, second_ae

contains

!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine astronomy_init()
integer :: unit, ierr, io, seconds, days

if(file_exist('input.nml')) then
  unit = open_file('input.nml', action='read')
  ierr = 1
  do while (ierr /= 0)
    read(unit,nml=astronomy_nml,iostat=io,end=10)
    ierr = check_nml_error(io,'astronomy_nml')
  enddo
10 continue
  call close_file (unit)
else
  call error_mesg('astronomy_init','input.nml does not exist.', FATAL)
endif

autumnal_eq_ref = set_date(year_ae,month_ae,day_ae, &
                           hour_ae,minute_ae,second_ae)

if(period == 0) then
  period_time_type = length_of_year()
  call get_time(period_time_type, seconds, days)
  period = 86400*days + seconds
else
  period_time_type = set_time(period,0)
endif

unit = open_file ('logfile.out', action='append')
if ( get_my_pe() == 0 ) then
     write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
     write (unit, nml=astronomy_nml)
endif
call close_file (unit)

pi    = 4*atan(1.)
twopi = 2*pi
deg_to_rad = twopi/360

call orbit
init = .true.

return
end subroutine astronomy_init
!-----------------------------------------------------------
function orbital_time(time) result(t)

type(time_type), intent(in) :: time
real :: t

!  the orbital time is the time (1 year = 2*pi) since autumnal equinox;
!  autumnal_eq_ref is a module variable of time_type and will have been 
!  defined by default 
!  or by a call to set_ref_date_of_ae; length_of_year is available
!  through the time manager and is set at the value approriate for 
!  the calandar being used

t = real((time - autumnal_eq_ref)//period_time_type)
t = twopi*(t - floor(t))
if(time < autumnal_eq_ref) t = twopi - t 

return
end function orbital_time

!------------------------------------------------------------

function universal_time(time) result(t)

type(time_type), intent(in) :: time
integer ::  seconds, days
real :: t

! univseral time is the time of day at longitude = 0.0 (1 day = 2*pi) 

call get_time(time, seconds, days)
t = twopi*real(seconds)/86400.0

return
end function universal_time

!--------------------------------------------------

subroutine get_orbital_parameters(ecc_out, obliq_out, per_out)

real, intent(out) :: ecc_out, obliq_out, per_out

ecc_out = ecc
obliq_out = obliq
per_out = per

return
end subroutine get_orbital_parameters

!---------------------------------------------------
subroutine set_orbital_parameters(ecc_in, obliq_in, per_in)

real, intent(in) :: ecc_in, obliq_in, per_in

if(ecc.lt.0.0.or.ecc.gt.0.99) &
    call error_mesg('SET_ORBITAL_PARAMTERS in ASTRONOMY_MOD', &
                              'ecc must be between 0 and 0.99', FATAL)

ecc = ecc_in
obliq = obliq_in
per = per_in

call orbit

return
end subroutine set_orbital_parameters

!---------------------------------------------------

subroutine get_period_integer(period_out)
integer, intent(out) :: period_out
integer :: seconds, days

call get_time(period_time_type, seconds, days)
period_out = 86400*days + seconds

return
end subroutine get_period_integer

!---------------------------------------------------

subroutine set_period_integer(period_in)
integer, intent(in) :: period_in

period_time_type = set_time(period_in, 0)

return
end subroutine set_period_integer
!---------------------------------------------------

subroutine get_period_time_type(period_out)
type(time_type), intent(out) :: period_out

period_out = period_time_type

return
end subroutine get_period_time_type

!---------------------------------------------------
subroutine set_period_time_type(period_in)
type(time_type), intent(in) :: period_in

period_time_type = period_in

return
end subroutine set_period_time_type

!---------------------------------------------------

subroutine get_ref_date_of_ae(day_out,month_out,year_out,&
                          second_out,minute_out,hour_out)

integer, intent(out) :: day_out, month_out, year_out,  &
                     second_out, minute_out, hour_out

   day_out =    day_ae
 month_out =  month_ae
  year_out =   year_ae
second_out = second_ae
minute_out = minute_ae
  hour_out =   hour_ae

return
end subroutine get_ref_date_of_ae

!---------------------------------------------------
subroutine set_ref_date_of_ae(day_in,month_in,year_in, &
                          second_in,minute_in,hour_in)

integer, intent(in) :: day_in, month_in, year_in
integer, optional :: second_in, minute_in, hour_in


   day_ae =    day_in
 month_ae =  month_in
  year_ae =   year_in

second_ae = 0
minute_ae = 0
hour_ae   = 0

if(present(second_in)) then
  second_ae = second_in
  minute_ae = minute_in
  hour_ae =   hour_in
endif

autumnal_eq_ref = set_date(year_ae,month_ae,day_ae, &
                            hour_ae,minute_ae,second_ae)

return
end subroutine set_ref_date_of_ae

!---------------------------------------------------

function angle(t)

real, intent(in) :: t
real :: angle, time_since_per, norm_time, x
integer :: int, int_1

norm_time = t*float(num_angles)/twopi
int = floor(norm_time)
int = modulo(int,num_angles)
int_1 = int+1

x = norm_time - floor(norm_time)
angle = (1.0 -x)*orb_angle(int) + x*orb_angle(int_1)
angle = modulo(angle, twopi)

return
end function angle
!----------------------------------------------------------

subroutine orbit

!  computes and stores a table of value of orbital angles as a function
!    of orbital time (both the angle and time are zero at autumnal
!    equinox in the NH, and range from 0 tp 2*pi) 
!    num_angles is a namelist parameter 

integer :: j
real :: d1, d2, d3, d4, d5, dt, norm

orb_angle(0) = 0.0

dt = twopi/float(num_angles)
norm = sqrt(1.0 - ecc**2)
dt = dt*norm

do j = 1, num_angles
  d1 = dt*r_inv_squared(orb_angle(j-1))
  d2 = dt*r_inv_squared(orb_angle(j-1)+0.5*d1)
  d3 = dt*r_inv_squared(orb_angle(j-1)+0.5*d2)
  d4 = dt*r_inv_squared(orb_angle(j-1)+d3)
  d5 = d1/6.0 + d2/3.0 +d3/3.0 +d4/6.0
  orb_angle(j) = orb_angle(j-1) + d5
end do
  
return
end subroutine orbit
!-----------------------------------------------------------

function r_inv_squared(ang)

real, intent(in) :: ang
real :: r_inv_squared, r, rad_per

rad_per = per*deg_to_rad
r = (1 - ecc**2)/(1. + ecc*cos(ang - rad_per))
r_inv_squared = r**(-2)

return
end function r_inv_squared

!----------------------------------------------------------

function declination(ang)

real, intent(in) :: ang
real :: declination

real :: rad_obliq, sin_dec

rad_obliq   =   obliq*deg_to_rad
sin_dec     = - sin(rad_obliq)*sin(ang)
declination =   asin(sin_dec)

return
end function declination

!------------------------------------------------------------

subroutine diurnal_solar_2d(cosz, solar, lat, lon, gmt, time_since_ae, dt) 

real, intent(in), dimension(:,:) :: lat, lon
real, intent(in) :: gmt, time_since_ae
real, intent(in), optional :: dt
real, intent(out), dimension(size(lat,1),size(lat,2)) :: cosz, solar


real, dimension(size(lat,1),size(lat,2)) :: t, tt, h, aa, bb, st, stt, sh
real :: ang, dec, rr

if(.not.init) call astronomy_init

if(time_since_ae.lt.0.0.or.time_since_ae.gt.twopi) &
     call error_mesg('DIURNAL_SOLAR in ASTRONOMY_MOD', &
                           'time_since_ae not between 0 and 2pi', FATAL)

ang = angle(time_since_ae)
dec = declination(ang)
rr  = r_inv_squared(ang)

aa = sin(lat)*sin(dec)
bb = cos(lat)*cos(dec)

t = gmt + lon - pi
where(t >= pi) t = t - twopi  
where(t < -pi) t = t + twopi   

! t = time of day with respect to local noon (2*pi = 1 day)

if(.not.present(dt)) then
  cosz = aa + bb*cos(t)

else  ! averaging from t to t+dt 
  tt = t + dt
  h   = half_day_2d(lat,dec)
  st  = sin(t)
  stt = sin(tt)
  sh  = sin(h)

  where(      t < -h .and.      tt < -h) cosz = 0.0
  where(      t < -h .and. abs(tt) <= h) cosz = aa*(tt + h ) + bb*(stt + sh)
  where(      t < -h .and.       h < tt) cosz = aa*( h + h ) + bb*( sh + sh)
  where( abs(t) <= h .and. abs(tt) <= h) cosz = aa*(tt - t ) + bb*(stt - st)
  where( abs(t) <= h .and.       h < tt) cosz = aa*( h - t ) + bb*( sh - st)
  where(      h <  t                   ) cosz = 0.0
  where( twopi - h < tt) cosz = cosz + aa*(tt + h - twopi) + bb*(stt + sh)
  cosz = cosz/dt
endif

cosz = max(0.0, cosz)
solar = cosz*rr
  
return
end subroutine diurnal_solar_2d
!------------------------------------------------------------

subroutine diurnal_solar_1d(cosz, solar, lat, lon, gmt, time_since_ae, dt)

real, intent(in), dimension(:) :: lat, lon
real, intent(in) :: gmt, time_since_ae
real, optional :: dt
real, intent(out), dimension(size(lat)) :: cosz, solar

real, dimension(size(lat),1) :: lat_2d, lon_2d, cosz_2d, solar_2d

lat_2d(:,1) = lat
lon_2d(:,1) = lon

if(present(dt)) then
  call diurnal_solar_2d &
            (cosz_2d, solar_2d, lat_2d, lon_2d, gmt, time_since_ae, dt)
else
    call diurnal_solar_2d &
            (cosz_2d, solar_2d, lat_2d, lon_2d, gmt, time_since_ae)
end if
cosz  = cosz_2d (:,1)
solar = solar_2d(:,1)

return
end subroutine diurnal_solar_1d

!------------------------------------------------------------

subroutine diurnal_solar_0d(cosz, solar, lat, lon, gmt, time_since_ae, dt)

real, intent(in) :: lat, lon, gmt, time_since_ae
real, optional :: dt
real, intent(out) :: cosz, solar

real, dimension(1,1) :: lat_2d, lon_2d, cosz_2d, solar_2d

lat_2d = lat
lon_2d = lon

if(present(dt)) then
  call diurnal_solar_2d &
            (cosz_2d, solar_2d, lat_2d, lon_2d, gmt, time_since_ae, dt)
else
    call diurnal_solar_2d &
            (cosz_2d, solar_2d, lat_2d, lon_2d, gmt, time_since_ae)
end if
cosz = cosz_2d(1,1)
solar = solar_2d(1,1)

return
end subroutine diurnal_solar_0d

!-------------------------------------------------------------

subroutine diurnal_solar_cal_2d(cosz, solar, lat, lon, time, dt_time) 

real, intent(in), dimension(:,:) :: lat, lon
type(time_type), intent(in) :: time
type(time_type), optional:: dt_time
real, intent(out), dimension(size(lat,1),size(lat,2)) :: cosz, solar

real :: dt
real :: gmt, time_since_ae

if(.not.init) call astronomy_init

gmt           = universal_time(time)
dt = 0.0
if(present(dt_time)) dt = universal_time(dt_time)
time_since_ae = orbital_time(time)

if(dt.gt.0.0) then
  call diurnal_solar_2d(cosz, solar, lat, lon, gmt, time_since_ae, dt)
else
  call diurnal_solar_2d(cosz, solar, lat, lon, gmt, time_since_ae)
end if

return
end subroutine diurnal_solar_cal_2d
!------------------------------------------------------------

subroutine diurnal_solar_cal_1d(cosz, solar, lat, lon, time, dt_time)

real, intent(in), dimension(:) :: lat, lon
type(time_type), intent(in) :: time
type(time_type), optional   :: dt_time
real, intent(out), dimension(size(lat)) :: cosz, solar


real, dimension(size(lat),1) :: lat_2d, lon_2d, cosz_2d, solar_2d

lat_2d(:,1) = lat
lon_2d(:,1) = lon

if(present(dt_time)) then
  call diurnal_solar_cal_2d(cosz_2d, solar_2d, lat_2d, lon_2d, time, dt_time)
else
  call diurnal_solar_cal_2d(cosz_2d, solar_2d, lat_2d, lon_2d, time)
end if
cosz  = cosz_2d (:,1)
solar = solar_2d(:,1)

return
end subroutine diurnal_solar_cal_1d
!------------------------------------------------------------
subroutine diurnal_solar_cal_0d(cosz, solar, lat, lon, time, dt_time)

real, intent(in) :: lat, lon
type(time_type), intent(in) :: time
type(time_type), optional   :: dt_time
real, intent(out) :: cosz, solar

real, dimension(1,1) :: lat_2d, lon_2d, cosz_2d, solar_2d

lat_2d = lat
lon_2d = lon

if(present(dt_time)) then
  call diurnal_solar_cal_2d(cosz_2d, solar_2d, lat_2d, lon_2d, time, dt_time)
else
  call diurnal_solar_cal_2d(cosz_2d, solar_2d, lat_2d, lon_2d, time)
end if

cosz = cosz_2d(1,1)
solar = solar_2d(1,1)

return
end subroutine diurnal_solar_cal_0d
!-------------------------------------------------------------

function half_day_2d(latitude, dec) result(h)

real, dimension(:,:), intent(in) :: latitude
real, intent(in) :: dec
real, dimension(size(latitude,1),size(latitude,2)):: &
   h, cos_half_day, lat
real :: tan_dec, eps

eps = 1.0e-05
tan_dec = tan(dec)

lat = latitude
where(latitude.eq. 0.5*pi) lat= latitude - eps
where(latitude.eq.-0.5*pi) lat= latitude + eps

cos_half_day = - tan(lat)*tan_dec
where(cos_half_day.le.-1.0)  h = pi
where(cos_half_day.ge.+1.0)  h = 0.0
where(cos_half_day.gt.-1.0.and.cos_half_day.lt.1.0) &
  h = acos(cos_half_day)

return
end function half_day_2d

!---------------------------------------------------------------

function half_day_0d(latitude, dec) result(h)

real, intent(in) :: latitude, dec
real :: h
real, dimension(1,1) :: lat_2d, h_2d

lat_2d = latitude
h_2d = half_day_2d(lat_2d,dec)
h = h_2d(1,1)

return
end function half_day_0d
!---------------------------------------------------------------

subroutine daily_mean_solar_2d(cosz, solar, lat, time_since_ae) 

real, intent(in), dimension(:,:) :: lat
real, intent(in) :: time_since_ae
real, intent(out), dimension(size(lat,1),size(lat,2)) :: solar, cosz
real, dimension(size(lat,1),size(lat,2)) :: h
real :: ang, dec, rr

if(.not.init) call astronomy_init

if(time_since_ae.lt.0.0.or.time_since_ae.gt.twopi) &
     call error_mesg('DAILY_MEAN_SOLAR in ASTRONOMY_MOD', &
                           'time_since_ae not between 0 and 2pi', FATAL)

ang = angle(time_since_ae)
dec = declination(ang)
rr  = r_inv_squared(ang)

h = half_day_2d(lat, dec)

where (h.eq.0.0)
  cosz = 0.0
elsewhere
  cosz = sin(lat)*sin(dec) + cos(lat)*cos(dec)*sin(h)/h
end where

solar = cosz*rr*h/pi



return
end subroutine daily_mean_solar_2d

!--------------------------------------------------------------
subroutine daily_mean_solar_1d(cosz, solar, lat, time_since_ae) 

real, intent(in), dimension(:) :: lat
real, intent(in) :: time_since_ae
real, intent(out), dimension(size(lat)) :: solar, cosz
real, dimension(size(lat),1) :: lat_2d, solar_2d, cosz_2d

lat_2d(:,1) = lat
call daily_mean_solar_2d(cosz_2d, solar_2d, lat_2d, time_since_ae)
solar = solar_2d(:,1)
cosz  = cosz_2d(:,1)

return
end subroutine daily_mean_solar_1d

!--------------------------------------------------------------
subroutine daily_mean_solar_0d(cosz, solar, lat, time_since_ae) 

real, intent(in) :: lat, time_since_ae
real, intent(out) :: solar, cosz
real, dimension(1,1) :: lat_2d, solar_2d, cosz_2d

lat_2d = lat
call daily_mean_solar_2d(cosz_2d, solar_2d, lat_2d, time_since_ae)
solar = solar_2d(1,1)
cosz  = cosz_2d(1,1)

return
end subroutine daily_mean_solar_0d

!--------------------------------------------------------------
subroutine daily_mean_solar_cal_2d(cosz, solar, lat, time) 

real, intent(in), dimension(:,:) :: lat
type(time_type), intent(in) :: time
real, intent(out), dimension(size(lat,1),size(lat,2)) :: solar, cosz

real :: time_since_ae

if(.not.init) call astronomy_init

time_since_ae = orbital_time(time)

if(time_since_ae.lt.0.0.or.time_since_ae.gt.twopi) &
     call error_mesg('DAILY_MEAN_SOLAR in ASTRONOMY_MOD', &
                           'time_since_ae not between 0 and 2pi', FATAL)

call daily_mean_solar_2d(cosz, solar, lat, time_since_ae) 

return
end subroutine daily_mean_solar_cal_2d

!--------------------------------------------------------------
subroutine daily_mean_solar_cal_1d(cosz, solar, lat, time) 

real, intent(in), dimension(:) :: lat
type(time_type), intent(in) :: time
real, intent(out), dimension(size(lat)) :: solar, cosz
real, dimension(size(lat),1) :: lat_2d, solar_2d, cosz_2d

lat_2d(:,1) = lat
call daily_mean_solar_cal_2d(cosz_2d, solar_2d, lat_2d, time)
solar = solar_2d(:,1)
cosz  = cosz_2d(:,1)

return
end subroutine daily_mean_solar_cal_1d

!--------------------------------------------------------------
subroutine daily_mean_solar_cal_0d(cosz, solar, lat, time) 

real, intent(in) :: lat
type(time_type), intent(in) :: time
real, intent(out) :: solar, cosz
real, dimension(1,1) :: lat_2d, solar_2d, cosz_2d

lat_2d = lat
call daily_mean_solar_cal_2d(cosz_2d, solar_2d, lat_2d, time)
solar = solar_2d(1,1)
cosz  = cosz_2d(1,1)

return
end subroutine daily_mean_solar_cal_0d

!--------------------------------------------------------------

subroutine annual_mean_solar_2d(cosz,solar,lat)

real, intent(in) :: lat(:,:)
real, intent(out), dimension(size(lat,1),size(lat,2)) :: solar, cosz
real, dimension(size(lat,1),size(lat,2)) :: s,z
real :: t
integer :: i

if(.not.init) call astronomy_init

solar = 0.0
cosz = 0.0
do i =1, num_angles
  t = float(i-1)*twopi/float(num_angles)
  call daily_mean_solar(z,s,lat,t)
  solar = solar + s
  cosz  = cosz  + z*s
end do
solar = solar/float(num_angles)
cosz  = cosz/float(num_angles)

where(solar.eq.0.0) 
  cosz = 0.0
elsewhere
  cosz = cosz/solar
end where

return
end subroutine annual_mean_solar_2d

!--------------------------------------------------------------

subroutine annual_mean_solar_1d(cosz,solar,lat)

real, intent(in) :: lat(:)
real, intent(out), dimension(size(lat)) :: solar, cosz

real, dimension(size(lat),1) :: lat_2d, solar_2d, cosz_2d

lat_2d(:,1) = lat
call annual_mean_solar_2d(cosz_2d, solar_2d, lat_2d)
solar = solar_2d(:,1)
cosz  =  cosz_2d(:,1)


return
end subroutine annual_mean_solar_1d



!----------------------------------------------------
end module astronomy_mod

