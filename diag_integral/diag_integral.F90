
module diag_integral_mod

!-----------------------------------------------------------------------
!
!    computes and outputs global integrals
!
!-----------------------------------------------------------------------

use     time_manager_mod, only:  time_type, get_time, set_time,  &
                                 operator(+),  operator(-),      &
                                 operator(==), operator(>=),     &
                                 operator(/=)
use        utilities_mod, only: file_exist, error_mesg, open_file,  &
                                check_nml_error, get_my_pe, FATAL,  &
                                close_file
use        constants_mod, only: radius

use              mpp_mod, only: mpp_sum

implicit none
private

!-----------------------------------------------------------------------
!------ interfaces ------

public   diag_integral_output, diag_integral_init, diag_integral_end,  &
         sum_diag_integral_field, diag_integral_field_init

!-----------------------------------------------------------------------

interface sum_diag_integral_field
   module procedure sum_field_2d,   &
                    sum_field_3d,   &
                    sum_field_wght_3d
end interface

!-----------------------------------------------------------------------
!------ namelist -------

   integer, parameter :: mxch = 64
   integer, parameter :: max_num_field = 10

   real                :: output_interval = -1.0
   character(len=8)    :: time_units = 'hours'
   character(len=mxch) :: file_name = ' '
   logical             :: print_header = .true.

   namelist /diag_integral_nml/ output_interval, time_units,  &
                                file_name, print_header

!-----------------------------------------------------------------------
!----- version number -----

   character(len=128) :: version = '$Id: diag_integral.F90,v 1.2 2000/08/04 18:43:06 fms Exp $'
   character(len=128) :: tag = '$Name: bombay $'

   type (time_type) :: Next_alarm_time, Alarm_interval,  &
                       Zero_time
   type (time_type) :: Time_init_save

   integer :: diag_unit = 0
   logical :: do_init   = .true.

!-----------------------------------------------------------------------
!---- global storage ----

   real, allocatable, dimension(:,:) :: area
   integer                           :: idim, jdim, field_size
   real                              :: sum_area

!-----------------------------------------------------------------------

   integer, parameter :: max_len_name   = 8
   integer, parameter :: max_field_name = 10

   integer                     :: num_field = 0
   character(len=max_len_name) :: field_name   (max_num_field)
   character(len=16)           :: field_format (max_num_field)
   real                        :: field_sum    (max_num_field)
   integer                     :: field_count  (max_num_field)

   integer            :: nt,             nd
   character(len=160) :: format_text,    format_data
   logical            ::              do_format_data = .true.
!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine diag_integral_output (Time)

   type (time_type), intent(in) :: Time

!-----------------------------------------------------------------------
!----- check alarm -------

      if ( diag_unit == 0 )                  return
      if ( .not. diag_integral_alarm(Time) ) return

!-----------------------------------------------------------------------
!---- output -----

      call write_field_averages (Time)

!---- increment diagnostics alarm -----

      Next_alarm_time = Next_alarm_time + Alarm_interval

!-----------------------------------------------------------------------

   end subroutine diag_integral_output

!#######################################################################

 subroutine diag_integral_end (Time)

   type (time_type), intent(in) :: Time

!-----------------------------------------------------------------------
!----- output here if alarm not turned on -----

      if ( diag_unit == 0 )              return
      if ( Alarm_interval /= Zero_time ) return

!-----------------------------------------------------------------------
!---- output -----

      call write_field_averages (Time)

!-----------------------------------------------------------------------

   end subroutine diag_integral_end

!#######################################################################

 subroutine write_field_averages (Time)

  type (time_type), intent(in) :: Time

  integer :: i, kount
  real    :: xtime, field_avg(max_num_field)

  real :: rcount
!-----------------------------------------------------------------------
!------ compute averages and output -------

      if (print_header)   call format_text_init
      if (do_format_data) call format_data_init

      do i = 1, num_field
         rcount = real(field_count(i))
         call mpp_sum (rcount, 1)
         call mpp_sum (field_sum(i), 1)
         field_count(i) = nint(rcount)

         if ( field_count(i) == 0 ) call error_mesg &
                     ('write_field_averages in diag_integral_mod',  &
                      'field_count equals zero for field_name ' //  &
                       field_name(i)(1:len_trim(field_name(i))), FATAL )
         kount = field_count(i) / field_size
         if ( field_size * kount /= field_count(i) ) call error_mesg &
                    ('write_field_averages in diag_integral_mod',  &
                     'field_count not a multiple of field_size', FATAL )

         field_avg  (i) = field_sum(i) / (sum_area * float(kount))
         field_sum  (i) = 0.0
         field_count(i) = 0
      enddo

      if ( get_my_pe() /= 0 ) return

      xtime = get_axis_time (Time-Time_init_save, time_units)

      write (diag_unit,format_data(1:nd)) &
                       xtime, (field_avg(i),i=1,num_field)

!-----------------------------------------------------------------------

   end subroutine write_field_averages

!#######################################################################

   subroutine diag_integral_init (Time_init, Time, blon, blat)

      type (time_type), intent(in) :: Time_init, Time
      real            , intent(in) :: blon(:), blat(:)
      
      integer :: unit, io, ierr, seconds, nc, i, j
      real, dimension(size(blat)) :: slat
      real :: r2
      real :: rsize

      Time_init_save = Time_init

!       ----- read namelist -----
!      ----- write namelist (to standard output) -----

      if ( file_exist('input.nml')) then
         unit = open_file ('input.nml', action='read')
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=diag_integral_nml, iostat=io, end=10)
            ierr = check_nml_error (io, 'diag_integral_nml')
         enddo
  10     call close_file (unit)
      endif

!      ----- write namelist -----
      unit = open_file ('logfile.out', action='append')
      if ( get_my_pe() == 0 ) then
           write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
           write (unit, nml=diag_integral_nml)
      endif
      call close_file (unit)

!----- grid setup ----

      idim = size(blon)-1
      jdim = size(blat)-1
      allocate (area(idim,jdim))
      field_size = idim*jdim

      r2 = radius*radius
      slat = sin(blat)
      do j = 1, jdim
      do i = 1, idim
         area(i,j) = r2 * (blon(i+1)-blon(i)) * (slat(j+1)-slat(j))
      enddo
      enddo

      sum_area = sum(area)

      rsize = real(field_size)
      call mpp_sum (rsize, 1)
      call mpp_sum (sum_area, 1)
      field_size = nint(rsize)

!     print *, 'from diag_integral_init: sum_area  = ', sum_area
!     print *, 'from diag_integral_init: 4*pi*r**2 = ', 8.*acos(0.0)*r2

!--- initialize diagnostics output unit/file ? ----

      if ( file_name(1:1) /= ' ' ) then
           nc = len_trim(file_name)
           diag_unit = open_file (file_name(1:nc), action='write')
      endif

      do_init = .false.

!-----------------------------------------------------------------------
!--------------------- initialize alarm --------------------------------

      Zero_time = set_time (0,0)

!---- set output interval -----

   if ( output_interval >= -0.01) then
      Alarm_interval = set_axis_time (output_interval, time_units)
      Next_alarm_time = Time + Alarm_interval
   else
!     --- zero time means the alarm is not set ---
      Alarm_interval = Zero_time
   endif

!---- set alarm to current time + interval ----

      Next_alarm_time = Time + Alarm_interval

!-----------------------------------------------------------------------

 end subroutine diag_integral_init

!#######################################################################

 function diag_integral_alarm (Time) result (answer)

 type (time_type), intent(in) :: Time
 logical                      :: answer

!-----------------------------------------------------------------------

    if (do_init) call error_mesg ('diag_integral_alarm',  &
                                  'must call diag_integral_init', FATAL)

!----- check the alarm -----

    answer = .false.
    if (Alarm_interval == Zero_time) return

!----- sound the alarm -----

    if (Time >= Next_alarm_time) answer = .true.

 end function diag_integral_alarm

!#######################################################################

   subroutine format_text_init

   integer :: i, nc

     if (diag_unit   == 0) return
     if (get_my_pe() /= 0) return

      nt = 11
      format_text(1:nt) = "('#    time"

      do i = 1, num_field
!        ---- text format ----
         nc = len_trim(field_name(i))
         format_text(nt+1:nt+nc+5) =  '     ' // field_name(i)(1:nc)
         nt = nt+nc+5
      enddo
         format_text(nt+1:nt+2) = "')"
         nt = nt+2

      write (diag_unit, format_text(1:nt))

      print_header = .false.

   end subroutine format_text_init

!#######################################################################

   subroutine format_data_init

   integer :: i, nc

      nd = 9
      format_data(1:nd) = '(1x,f10.2'

      do i = 1, num_field
!        ---- data format ----
         nc = len_trim(field_format(i))
         format_data(nd+1:nd+nc+5) =  ',1x,' // field_format(i)(1:nc)
         nd = nd+nc+5
      enddo
      format_data(nd+1:nd+1) = ')'
      nd=nd+1

      do_format_data = .false.

   end subroutine format_data_init

!#######################################################################

 function get_axis_time (Time, units) result (atime)

   type(time_type),  intent(in) :: Time
   character(len=*), intent(in) :: units
   real                         :: atime
   integer                      :: sec, day

!---- returns real time in the time axis units ----
!---- convert time type to appropriate real units ----

      call get_time (Time, sec, day)

      if (units(1:3) == 'sec') then
         atime = float(sec) + 86400.*float(day)
      else if (units(1:3) == 'min') then
         atime = float(sec)/60. + 1440.*float(day)
      else if (units(1:3) == 'hou') then
         atime = float(sec)/3600. + 24.*float(day)
      else if (units(1:3) == 'day') then
         atime = float(sec)/86400. + float(day)
      endif

 end function get_axis_time

!#######################################################################

 function set_axis_time (atime, units) result (Time)

   real,             intent(in) :: atime
   character(len=*), intent(in) :: units
   type(time_type)  :: Time
   integer          :: sec, day = 0

!---- returns time type given real time in axis units ----
!---- convert real time units to time type ----

      if (units(1:3) == 'sec') then
         sec = int(atime + 0.5)
      else if (units(1:3) == 'min') then
         sec = int(atime*60. + 0.5)
      else if (units(1:3) == 'hou') then
         sec = int(atime*3600. + 0.5)
      else if (units(1:3) == 'day') then
         sec = int(atime*86400. + 0.5)
      endif

      Time = set_time (sec, day)

 end function set_axis_time

!#######################################################################

 subroutine diag_integral_field_init (name, format)

   character(len=*), intent(in) :: name, format

   integer :: field

   if (len(name) > max_len_name ) call error_mesg  &
                    ('diag_integral_field_init', 'name too long', FATAL)

   field = get_field_index (name)

   if (field /= 0) call error_mesg ('diag_integral_field_init', &
                                    'name already exists', FATAL)


   num_field = num_field + 1

   if (num_field > max_num_field) call error_mesg  &
                                  ('diag_integral_field_init', &
                                  'too many fields initialized', FATAL)

   field_name   (num_field) = name
   field_format (num_field) = format
   field_sum    (num_field) = 0.0
   field_count  (num_field) = 0


 end subroutine diag_integral_field_init

!#######################################################################

 subroutine sum_field_2d (name, data, is, js)

   character(len=*),  intent(in) :: name
   real,              intent(in) :: data(:,:)
   integer, optional, intent(in) :: is, js

   integer :: field, i1, j1, i2, j2

   field = get_field_index (name)
   if (field == 0) call error_mesg ('sum_diag_integral_field', &
                                    'field does not exist', FATAL)

   i1 = 1;  if (present(is)) i1 = is
   j1 = 1;  if (present(js)) j1 = js
   i2 = i1 + size(data,1) - 1
   j2 = j1 + size(data,2) - 1

   field_count (field) = field_count (field) + size(data,1)*size(data,2)
   field_sum   (field) = field_sum   (field) +  &
                               sum (data * area(i1:i2,j1:j2))

 end subroutine sum_field_2d

!#######################################################################

 subroutine sum_field_3d (name, data, is, js)

   character(len=*),  intent(in) :: name
   real,              intent(in) :: data(:,:,:)
   integer, optional, intent(in) :: is, js

   real, dimension(size(data,1),size(data,2)) :: data2
   integer :: field, i1, j1, i2, j2

   field = get_field_index (name)
   if (field == 0) call error_mesg ('sum_diag_integral_field', &
                                    'field does not exist', FATAL)

   i1 = 1;  if (present(is)) i1 = is
   j1 = 1;  if (present(js)) j1 = js
   i2 = i1 + size(data,1) - 1
   j2 = j1 + size(data,2) - 1

   field_count (field) = field_count (field) + size(data,1)*size(data,2)

   data2 = sum(data,3)
   field_sum   (field) = field_sum   (field) +  &
                               sum (data2 * area(i1:i2,j1:j2))

 end subroutine sum_field_3d

!#######################################################################

 subroutine sum_field_wght_3d (name, data, wt, is, js)

   character(len=*),  intent(in) :: name
   real,              intent(in) :: data(:,:,:), wt(:,:,:)
   integer, optional, intent(in) :: is, js

   real, dimension(size(data,1),size(data,2)) :: data2
   integer :: field, i1, j1, i2, j2

   field = get_field_index (name)
   if (field == 0) call error_mesg ('sum_diag_integral_field', &
                                    'field does not exist', FATAL)

   i1 = 1;  if (present(is)) i1 = is
   j1 = 1;  if (present(js)) j1 = js
   i2 = i1 + size(data,1) - 1
   j2 = j1 + size(data,2) - 1

   field_count (field) = field_count (field) + size(data,1)*size(data,2)

   data2 = vert_diag_integral (data, wt) 
   field_sum   (field) = field_sum   (field) +  &
                               sum (data2 * area(i1:i2,j1:j2))

 end subroutine sum_field_wght_3d

!#######################################################################

 function vert_diag_integral (data, wt) result (data2)

   real, intent(in), dimension(:,:,:) :: data, wt
   real, dimension(size(data,1),size(data,2)) :: data2, wt2

     wt2 = sum(wt,3)
     if (count(wt2 == 0.) > 0) call error_mesg ('vert_diag_integral',  &
                               'vert sum of weights equals zero', FATAL)

     data2 = sum(data*wt,3) / wt2

 end function vert_diag_integral

!#######################################################################

 function get_field_index (name) result (index)
   character(len=*),  intent(in) :: name

   integer :: i, nc, index
   character(len=max_len_name) :: fname

   nc = len_trim (name)
   if (nc > max_len_name) call error_mesg ('get_field_index',  &
                                           'name too long', FATAL)
   index = 0

   do i = 1, num_field
      if ( name(1:nc) == field_name(i)(1:len_trim(field_name(i))) ) then
           index = i
           exit
      endif
   enddo

 end function get_field_index

!#######################################################################

end module diag_integral_mod

