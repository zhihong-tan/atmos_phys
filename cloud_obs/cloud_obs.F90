
                    module cloud_obs_mod

!-----------------------------------------------------------------------
!
!           sets up observed (climatological) clouds
!
!-----------------------------------------------------------------------

use horiz_interp_mod, only: horiz_interp
use    utilities_mod, only: file_exist, error_mesg, FATAL,  &
                            open_file, close_file,          &
                            check_nml_error, get_my_pe
use time_manager_mod, only: time_type, get_date
use  time_interp_mod, only: time_interp

implicit none
private

!---------- public interfaces ----------

public  cloud_obs, cloud_obs_init

!-----------------------------------------------------------------------
!   ---------- private data ------------

   character(len=128) :: version = '$Id: cloud_obs.F90,v 1.2 2000/08/04 18:18:22 fms Exp $'
   character(len=128) :: tag = '$Name: eugene $'

      real, allocatable, dimension(:,:,:) :: clda,cldb
      real, allocatable, dimension(:)     :: londat,latdat
   integer :: yrclda=-99,moclda=-99, yrcldb=-99,mocldb=-99
   logical :: do_init=.true.

!   ---------- namelist ---------------

   logical :: use_climo = .true.
   integer :: verbose = 0

   namelist /cloud_obs_nml/ use_climo, verbose


!------------ input grid parameters ----------
     integer, parameter :: mobs=144, nobs=72
        real :: sb, wb, dx, dy
!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine cloud_obs ( is, js, Time, cldamt )

!-----------------------------------------------------------------------
!    routine that reads monthly records of climatological
!    isccp cloud amount and then linearly interpolates between them
!-----------------------------------------------------------------------
!     input
!     -----
!     is, js   starting i,j indices (dimension(2))
!     Time     current time (time_type)
!
!     output
!     ------
!     cldamt    cloud amount data on horizontal grid,
!               dimensioned ix x jx x 3, for high,med, & low clouds.
!-----------------------------------------------------------------------
        integer, intent(in)                    :: is, js
type(time_type), intent(in)                    :: Time
           real, intent(out), dimension(:,:,:) :: cldamt
!-----------------------------------------------------------------------
      real, dimension(mobs,nobs,3) :: obs

   integer  day,month,year, second,minute,hour
   integer  month1,month2,mo,year1,year2,yr,unit,irec,n
   integer  ie,je,ix,jx,pe
      real  dmonth,dif
   logical,save :: useclimo1,useclimo2
   logical      :: unit_opened
!-----------------------------------------------------------------------

   if (do_init) call error_mesg ('cloud_obs',  &
                          'cloud_obs_init has not been called.',FATAL)

   if (size(cldamt,3) < 3) call error_mesg ('cloud_obs',  &
                                'dimension 3 of cldamt is < 3', FATAL)

   pe = get_my_pe()

!------------ size & position of this window ---------------------------
      ix=size(cldamt,1); jx=size(cldamt,2)
      ie=is+ix-1;        je=js+jx-1

!  --- check existence of cloud data set --------

      if (.not.file_exist('INPUT/cloud_obs.data')) then
        call error_mesg ('observed_cloud',  &
                    'file INPUT/cloud_obs.data does not exist.', FATAL)
      endif

!-----------------------------------------------------------------------
! ---- time interpolation for months -----

      call time_interp (Time, dmonth, year1, year2, month1, month2)
      
! ---- force climatology ----

      if (use_climo) then
          year1 = 0; year2 = 0
      endif

!-----------------------------------------------------------------------

      unit_opened=.false.

!    assumption is being made that the record for (year1,month1)
!    precedes the record for (year2,month2)

      if (year1 .ne. yrclda .or. month1 .ne. moclda) then
          unit_opened=.true.
          unit = open_file ( 'INPUT/cloud_obs.data', form='ieee32', &
                             action='read' )
          irec=0
          do
!!!!               read (unit,end=380)  yr,mo,obs
             yr=0; read (unit,end=380)     mo,obs
             irec=irec+1
             dif=12*(year1-yr)+month1-mo
             if (dif == 0) then
                yrclda=yr
                moclda=mo
                useclimo1=.false.
                if (yr == 0) useclimo1=.true.
                exit
             endif
!           --- otherwise use climo ---
             if (yr == 0 .and. month1 == mo) then
                yrclda=yr
                moclda=mo
                useclimo1=.true.
                exit
             endif
          enddo
          do n=1,3
            call horiz_interp (obs(:,:,n),wb,sb,dx,dy,  &
                               londat,latdat,clda(:,:,n),  &
                               verbose=verbose)
          enddo
      endif

      if (year2 .ne. yrcldb .or. month2 .ne. mocldb) then
          if (.not.unit_opened) then
             unit_opened=.true.
             unit = open_file ( 'INPUT/cloud_obs.data', form='ieee32', &
                                action='read' )
          endif
          if (useclimo1 .and. month2 <= month1 ) then
             if (verbose > 1 .and. pe == 0)  &
                       print *, ' rewinding INPUT/cloud_obs.data'
             rewind unit
          endif
          irec=0
          do
!!!!               read (unit,end=380)  yr,mo,obs
             yr=0; read (unit,end=380)     mo,obs
             irec=irec+1
             dif=12*(year2-yr)+month2-mo
             if (dif == 0) then
                yrcldb=yr
                mocldb=mo
                useclimo2=.false.
                if (yr == 0) useclimo2=.true.
                exit
             endif
!           --- climo ---
             if (yr == 0 .and. month2 == mo) then
                yrcldb=yr
                mocldb=mo
                useclimo2=.true.
                exit
             endif
          enddo
          do n=1,3
            call horiz_interp (obs(:,:,n),wb,sb,dx,dy,  &
                               londat,latdat,cldb(:,:,n),  &
                               verbose=verbose)
          enddo
      endif
          goto 381

 380  if (pe == 0) print *, ' month1,month2=',month1,month2
      if (pe == 0) print *, ' useclimo1,useclimo2=',useclimo1,useclimo2
      call error_mesg ('observed_cloud',  &
                       'eof reading file=INPUT/cloud_obs.data', FATAL)

 381  continue

   if (unit_opened) then
      call close_file (unit)
      if (verbose > 0 .and. pe == 0) then
         call get_date (Time, year, month, day, hour, minute, second)
         write (*,600) year,month,day, hour,minute,second
 600     format (/,'from cloud_obs:',   &
                 /,' date(y/m/d h:m:s) = ', &
                     i4,2('/',i2.2),1x,2(i2.2,':'),i2.2)
         print *, ' dmonth=',dmonth
         print *, ' year1,month1, yrclda,moclda, useclimo1=',  &
                    year1,month1, yrclda,moclda, useclimo1
         print *, ' year2,month2, yrcldb,mocldb, useclimo2=',  &
                    year2,month2, yrcldb,mocldb, useclimo2
         print *, ' '
      endif
   endif

!------------ time interpolation ---------------------------------------

      do n=1,3
         cldamt(:,:,n)=clda(is:ie,js:je,n)+  &
                       dmonth*(cldb(is:ie,js:je,n)-clda(is:ie,js:je,n))
      enddo

!-----------------------------------------------------------------------

 end subroutine cloud_obs

!#######################################################################

 subroutine cloud_obs_init (lonb,latb)

!-----------------------------------------------------------------------
!  lonb  =   longitude in radians of the grid box edges
!  latb  =   longitude in radians of the grid box edges
!-----------------------------------------------------------------------
   real, intent(in), dimension(:) :: lonb,latb
!-----------------------------------------------------------------------
   real    :: hpie
   integer :: i,j,in,jn, unit, ierr, io

!------- read namelist --------

      if (file_exist('input.nml')) then
          unit = open_file ('input.nml', action='read')
          ierr=1; do while (ierr /= 0)
             read  (unit, nml=cloud_obs_nml, iostat=io, end=10)
             ierr = check_nml_error(io,'cloud_obs_nml')
          enddo
  10      call close_file (unit)
      endif

!------- write version number and namelist ---------

      unit = open_file ('logfile.out', action='append')
      if ( get_my_pe() == 0 ) then
           write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
           write (unit, nml=cloud_obs_nml)
      endif
      call close_file (unit)

!------- setup for observed grid -------

      hpie=acos(0.0)
      sb=-hpie; wb=0.0; dx=4.0*hpie/float(mobs); dy=2.0*hpie/float(nobs)

!------- setup for data grid -------

      in=size(lonb,1); jn=size(latb,1)
      allocate (londat(in), latdat(jn))
      allocate (clda(in-1,jn-1,3), cldb(in-1,jn-1,3))

      latdat = latb; londat = lonb

      do_init = .false.

!-----------------------------------------------------------------------

 end subroutine cloud_obs_init

!#######################################################################

end module cloud_obs_mod

