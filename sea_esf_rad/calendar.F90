
                  module calendar_mod

use time_manager_mod,        only: time_type, days_in_year, &
				   get_time, get_date

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    this module  is needed by the compiler to mimic a skyhi
!    module when running FMS and is used in the standalone program.
!
!--------------------------------------------------------------------



!--------------------------------------------------------------------
!------------- ****** VERSION NUMBER ****** -------------------------

!	  character(len=5), parameter  :: version_number = 'v0.09'
          character(len=128)  :: version =  '$Id: calendar.F90,v 1.2 2001/08/30 15:12:43 fms Exp $'
          character(len=128)  :: tag     =  '$Name: galway $'




!--------------------------------------------------------------------
!---- interfaces ------

public calendar_4seasons_harmonic, &
       calendar_o3monloss

!--------------------------------------------------------------------
!--------------------------------------------------------------------


                  contains

 

subroutine calendar_4seasons_harmonic (lag, xday, arg, Time)

type(time_type), intent(in), optional    :: Time
real,            intent(out)             :: xday, arg
real,            intent(in)              :: lag

	integer             :: ndiy, year, month, day, &
			       hour, minute, second
	integer             :: m
	real                :: po4, fndiy, lag2, rday
        real, dimension(12) :: daypmn  

        data daypmn/        &
           31.0E+00, 28.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, &
           31.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, 31.0E+00 /

!-------------------------------------------------------------------
!   define position of current time within "year".
!---------------------------------------------------------------------
	ndiy = days_in_year (Time)
	if (ndiy == 366) daypmn(2) = 29.0

        po4 = 0.25E+00*ndiy
	fndiy = ndiy
	if (ndiy == 366.0) then
	  lag2 = lag + 1.0
        endif

        call get_date (Time, year, month, day, hour, minute, second)
	
	do m=1, month-1
	day = day + daypmn(m) 
	end do
	rday = day + FLOAT(hour)/24.0 + FLOAT(minute)/1440.0 + &
              FLOAT(second)/86400.0
	xday = MOD(rday + (po4 - lag2) + ndiy, fndiy)
	arg = 2.0E+00*acos(-1.0)*xday/ndiy

end subroutine calendar_4seasons_harmonic 



!####################################################################

subroutine calendar_o3monloss (o3monf, im_out, im1_out, imp_out, Time)

type(time_type), intent(in), optional    :: Time
real,            intent(out)   :: o3monf
integer,         intent(out)   :: im_out, im1_out, imp_out

!----------------------------------------------------------------------
!    define the cumulative days of the year at the end of each month,
!     beginning with december and ending with december (dy). also define
!     the number of days per month for a non leap year (daypmn) and the
!     number of days from 15th of the month to 15th of the next month
!     (delday).
!----------------------------------------------------------------------
         real, dimension(12)   ::  daypmn, delday 
	 integer               ::  im, ndayom, ihour, imp, im1, year, &
				   minute, second
	 real                  ::  deld, dpm

         data daypmn/        &
           31.0E+00, 28.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, &
           31.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, 31.0E+00 /

         data delday/   &
           31.0E+00, 28.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, &
           31.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, 31.0E+00 /

	 integer                 :: ndiy, year, month, day, &
				    hour, minute, second
!--------------------------------------------------------------------- 
!    compute ozone depletion profile factors for appropriate day and
!    month using above day/month values.
!---------------------------------------------------------------------
 
	 ndiy = days_in_year (Time)
	 call get_date (Time, year, im, ndayom, ihour, minute, second)
         if (ndayom .GE. 15) then
           imp = MOD(im,12) + 1
           im1 = im

!---------------------------------------------------------------------
!     if the date is past the 15th of the month, find the fraction of
!    the month between the 15th of this month and the 15th of next
!     month which has elapsed (o3monf). account for leap year when nec-
!     essary. deld is the number of days between the 15th of this month
!    and the 15th of next month.
!---------------------------------------------------------------------
           if (ndiy .EQ. 366 .AND. im .EQ. 2) then
             deld = delday(2) + 1
           else
             deld = delday(im)
           endif
           o3monf = (ndayom + ihour/24.0E+00 - 15.0E+00)/deld
         else
!---------------------------------------------------------------------
!     if the date precedes the 15th of the month, find the fraction of
!     the month between the 15th of last month and the 15th of this
!     month which has elapsed (o3monf). account for leap year when nec-
!     essary. deld is the number of days between the 15th of last month
!     and the 15th of this month.
!---------------------------------------------------------------------
            if (im .EQ. 1) then
              imp = 12
            else
              imp = im - 1
           endif
           im1 = imp
           if (ndiy .EQ. 366 .AND. imp .EQ. 2) then
             deld = delday(2) + 1
             dpm  = daypmn(2) + 1
           else
             deld = delday(imp)
             dpm  = daypmn(imp)
            endif
            o3monf = -(dpm + ndayom + ihour/24.0E+00 - 15.0E+00)/deld
          endif
         im_out = im
	 im1_out = im1
         imp_out = imp
  


end subroutine calendar_o3monloss 




             end module calendar_mod
