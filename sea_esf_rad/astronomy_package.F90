                module astronomy_package_mod

use time_manager_mod,    only:  time_type, set_time, get_date_julian, &
			        set_date_julian, get_time, operator(-)
use astronomy_mod,       only:  daily_mean_solar, diurnal_solar, &
			        annual_mean_solar
use radiation_diag_mod,  only:  radiag_from_astronomy
use constants_new_mod,   only:  pie, calyear, radians_to_degrees, &
                                secday
use rad_utilities_mod,   only:  radiation_control_type, Rad_control, &
				environment_type, Environment
use stratchem_mod,       only:  inquire_stratchem
use utilities_mod,       only:  open_file, file_exist,    &
                                check_nml_error, error_mesg, &
                                print_version_number, get_my_pe, &
				close_file, FATAL, NOTE,   &
				get_domain_decomp

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                      astronomy module
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!         character(len=5), parameter  ::  version_number = 'v0.09'
	  character(len=128)  :: version =  '$Id: astronomy_package.F90,v 1.2 2001/07/05 17:27:17 fms Exp $'
          character(len=128)  :: tag     =  '$Name: eugene $'



!---------------------------------------------------------------------
!-------  interfaces --------


public       &
              astronomy_package_init, astronomy,  &
	      astronomy_driver, astronomy_dealloc,  &
	      get_astronomy_for_sfcalbedo,   &
	      get_astronomy_for_swrad_init, &
	      get_astronomy_for_swrad,  &
	      get_astronomy_for_clouds_init, &
              get_astronomy_for_clouds,  &
	      return_cosz, get_solar_distance

interface get_astronomy_for_swrad
     module procedure    get_astronomy_for_swrad_3d, &
			 get_astronomy_for_swrad_2d
end interface

private       &
	       zenitha, getastronomy, sol88, sol88zen, &
	       sol88zeng, zenithb, gptinit

!---------------------------------------------------------------------
!-------- namelist  ---------


logical                     :: do_annual=.false.
logical                     :: ldiurn=.true.
logical                     :: lswg=.false.
integer                     :: nswg=0
logical                     :: do_oldrightasc=.false.
integer                     :: verbose = 0

namelist /astronomy_package_nml /      &
                                   ldiurn, lswg,     &
                                   do_oldrightasc,    &
                                   nswg, do_annual, &
				   verbose

!--------------------------------------------------------------------
!------   public data ----------


!--------------------------------------------------------------------
!------   private data ----------

real, dimension(:),     allocatable :: zhangle1, zfracday1, zrlat
real, dimension(:,:),   allocatable :: zrlong
real, dimension(:,:),   allocatable :: hangle, fracday, zenith, szenij,&
				       rlat, rlong, ssolwd
real, dimension(:,:,:), allocatable :: cosangsolar         
real, dimension(:),     allocatable :: gpt

integer                             :: nsolwg
real                                :: rsun, declin, rightasc, slag
real                                :: rrsun
integer                             :: IMINP, IMAXP, JMINP, JMAXP
integer                             :: x(4), y(4), idf, jdf



!--------------------------------------------------------------------
!--------------------------------------------------------------------




                       contains




subroutine astronomy_package_init (zrlat_in, zrlong_in  )

!--------------------------------------------------------------------
real, dimension(:),   intent(in)  ::  zrlat_in
real, dimension(:,:), intent(in)  ::  zrlong_in
!-------------------------------------------------------------------

      integer    ::  unit, ierr, io
      integer    :: x(4), y(4), j

!---------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=astronomy_package_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'astronomy_package_nml')
        enddo
10     call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'astronomy_package',  &
!				 version_number)
      if  (get_my_pe() == 0)  then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=astronomy_package_nml)
      endif
      call close_file (unit)

!-------------------------------------------------------------------

      call get_domain_decomp (x, y)
      jdf = y(4) -y(3) + 1 
      idf = x(4) -x(3) + 1 

!--------------------------------------------------------------------
!    allocate module variables
!-------------------------------------------------------------------
      allocate  (zhangle1  (jdf) )
      allocate  (zfracday1 (jdf) )
      allocate  (zrlat     (jdf) )
      allocate  (zrlong    (idf, jdf) )

!--------------------------------------------------------------------
!    save model latitudes and longitudes for use within module
!-------------------------------------------------------------------
      zrlong(:,:) = zrlong_in(:,:)
      zrlat(:) = zrlat_in(:)

 
!---------------------------------------------------------------------
!     check for a single, consistent zenith angle specification. if lwsg
!     is true, then nsolwg should be either 1,2,4, or 8.  otherwise
!     there is an error, and the model is stopped. 
!---------------------------------------------------------------------
 
      if (lswg .AND. (nswg .NE. 1) .AND. (nswg .NE. 2) .AND.  &
         (nswg .NE. 4) .AND. (nswg .NE. 8))  then
        call error_mesg ( 'astronomy_package_init', &
                      ' nsolwg does NOT have an acceptable value.', &
                                                              FATAL)
      endif
      if (ldiurn .AND. lswg) then
        call error_mesg ( 'astronomy_package_init', &
            'both lswg and ldiurn are true: ONLY ONE can be true !', &
                                                                FATAL)
      endif

      if (ldiurn .AND. do_annual) then
        call error_mesg ( 'astronomy_package_init', &
         'both do_annual and ldiurn are true: ONLY ONE can be true !', &
                                                                FATAL)
      endif

      if (lswg .AND. do_annual) then
        call error_mesg ( 'astronomy_package_init', &
         'both do_annual and lswg are true: ONLY ONE can be true !', &
                                                                FATAL)
      endif

      if (Environment%using_sky_periphs) then
	if (do_annual) then
          call error_mesg ( 'astronomy_package_init', &
         'do_annual is not valid when using skyhi peripherals !', FATAL)
        endif
      endif

      if (Environment%using_fms_periphs) then
	if (lswg) then
          call error_mesg ( 'astronomy_package_init', &
         'lswg is not valid when using fms peripherals !', FATAL)
        endif
      endif

!--------------------------------------------------------------------
!     define the number of zenith angles that are being used.
!--------------------------------------------------------------------
      if (lswg) then
	nsolwg = nswg
      else
	nsolwg = 1
      endif

!--------------------------------------------------------------------
!    call gptinit to define gaussian quadrature points.
!--------------------------------------------------------------------
      call gptinit


end subroutine astronomy_package_init


!####################################################################

subroutine astronomy  (lrad, fjd, jld)

!---------------------------------------------------------------------
!     astronomy is an interface between the astronomy package and the
!     GCM.  it determines characteristics of the relation of the earth  
!     to the sun, including declination, julian day, earth-sun distance,
!     hour angle and fraction of daylight.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!       fjd =  current fraction of julian day.
!       jld =  current julian day number (integral portion)
!--------------------------------------------------------------------
logical,                 intent(in) :: lrad
real,                    intent(in) :: fjd
integer,                 intent(in) :: jld
!--------------------------------------------------------------------

      real, dimension(jdf) ::  cc, ss
      real                 ::  ddlt, fid, xmin, fix, xsec, halp, fih, &
			       ymin, fiy,  asec
      integer              ::  idlt, ix, ihalp, iy, j 
      character(len=1)     ::  dsig,sign,sigb
      data sign /'-'/
      data sigb /' '/

!-----------------------------------------------------------------------
!     call sol88 to define time-dependent astronomical quantities (dec-
!     lination, right ascension, slag, earth-sun distance).
!-----------------------------------------------------------------------
      call sol88 (fjd, jld)

!---------------------------------------------------------------------
!     write out the current astronomical parameters (radius vector, the 
!     right ascension, and the declination of the sun, in both 
!     hours/min/sec and in fractions of hours or degrees, respectively. 
!---------------------------------------------------------------------
      ddlt  = 180.0E+00*declin/pie 
      idlt  = ddlt
      fid   = idlt
      if ((ddlt .LT. 0.0E+00) .AND. (idlt .EQ. 0)) then
        dsig  = sign
      else
        dsig  = sigb
      endif
      xmin  = ABS(ddlt - fid)*60.0E+00
      ix    = xmin
      fix   = ix
      xsec  = ABS(xmin - fix)*60.0E+00
      halp  = 12.0E+00*rightasc/pie 
      if (halp .LT. 0.0E+00) halp = halp + 24.0E+00
      ihalp = halp
      fih   = ihalp
      ymin  = ABS(halp - fih)*60.0E+00
      iy    = ymin
      fiy   = iy
      asec  = ABS(ymin - fiy)*60.0E+00
      if (lrad .and. verbose > 0  .and. get_my_pe() == 0) then
        write (*, 90060) rsun, ihalp, iy, asec, dsig,  &
			 idlt, ix, xsec
        write (*, 90070) halp, ddlt
      endif
 
!--------------------------------------------------------------------
!     compute hour angle of sunrise/sunset and daylight fraction at 
!     each model latitude.
!----------------------------------------------------------------------
      do j=1,jdf
        ss(j) = SIN(zrlat(j))*SIN(declin)
        cc(j) = COS(zrlat(j))*COS(declin)
        if (ABS(cc(j)) .GE. ABS(ss(j))) then
          zhangle1(j) = ACOS(-ss(j)/cc(j))
          zfracday1(j) = zhangle1(j)/pie
        elseif(ABS(cc(j)) .LT. ABS(ss(j)) .AND. ss(j) .GT. 0.0E+00) then
	  zhangle1(j) = pie
          zfracday1(j) = 1.0E+00
        elseif(ABS(cc(j)) .LT. ABS(ss(j)) .AND. ss(j) .LE. 0.0E+00) then
	  zhangle1(j) = 0.0E+00
          zfracday1(j) = 0.0E+00
        end if
      end do  

!--------------------------------------------------------------------
90060 format(/' radius vector=',1pe19.11,/,   &
              ' r.a. of sun=',i3,'hr',i3,'min',0pf6.2,'sec',5x,   &
              ' declination=',a1,i3,'deg',i3,'min',0pf5.1,'sec')
90070 format( '          or ',f11.7,' hrs',19x,' or ',f12.7,' degs')
!--------------------------------------------------------------------


end  subroutine astronomy



!#####################################################################

subroutine astronomy_driver (is, ie, js, je, dt, dt_rad, do_rad, &
			     Time, lat, lon, Rad_time,  &
			     do_average, fjd, fjulan)

!-------------------------------------------------------------------
integer, intent(in)                         :: is, ie, js, je, dt_rad
real, intent(in)                            :: dt
logical, intent(in)                         :: do_rad
type(time_type), intent(in),optional        :: Time, Rad_time
real, dimension(:,:), intent(in), optional  :: lat, lon
logical, intent(in),optional                :: do_average
real, intent(in), optional                  :: fjd, fjulan

!------------------------------------------------------------------
     integer, dimension(:), allocatable :: iabsp, jabsp
     integer                            :: j, i

!------------------------------------------------------------------

     jminp = 1
     jmaxp = je-js+1
     allocate (jabsp(jminp:jmaxp))

     do j=jminp,jmaxp
       jabsp(j) = js+j-1
     end do

     iminp = 1
     imaxp = ie-is+1
     allocate (iabsp(iminp:imaxp) )

     do i=iminp,imaxp
       iabsp(i) = is+i-1
     end do

     if (Environment%running_gcm) then
       if (Environment%running_fms) then
	 if (present(fjd) ) then   ! this implies using_sky_periphs
           call getastronomy (iabsp, jabsp, do_rad, dt_rad, dt,  &
	  	 	      Rad_time_sv=Rad_time, lat_sv=lat,   &
			      do_average=do_average,    &
			      fjulan=fjulan, fjd=fjd)
	 else    ! this implies using_fms_periphs
           call getastronomy (iabsp, jabsp, do_rad, dt_rad, dt,  &
	         	      Rad_time_sv=Rad_time, lat_sv=lat,  &
			      do_average=do_average)
	 endif
       else if (Environment%running_skyhi) then
         call getastronomy (iabsp, jabsp, do_rad, dt_rad, dt,  &
			    fjulan=fjulan, fjd=fjd)
       endif
     else if (Environment%running_standalone) then
       call getastronomy (iabsp, jabsp, do_rad, dt_rad, dt, fjd=fjd, &
		          fjulan=fjulan)
     endif

     deallocate (iabsp)
     deallocate (jabsp)


end subroutine astronomy_driver 
 

!####################################################################

subroutine astronomy_dealloc

    deallocate (cosangsolar)
    deallocate (zenith) 
    deallocate (fracday)
    deallocate (ssolwd)
    if (allocated (szenij)) then
      deallocate (szenij)
    endif


end subroutine astronomy_dealloc



!###################################################################

subroutine get_astronomy_for_sfcalbedo (zenith_out)

real, dimension(:,:), intent(out)   ::  zenith_out
 
      zenith_out(:,:) = zenith(:,:)

end subroutine get_astronomy_for_sfcalbedo



!##################################################################

subroutine get_astronomy_for_swrad_init (ldiurn_out, lswg_out, &
					 nsolwg_out, do_annual_out)

!-------------------------------------------------------------------
      logical lswg_out, ldiurn_out, do_annual_out
      integer nsolwg_out
!-------------------------------------------------------------------

      ldiurn_out = ldiurn
      lswg_out = lswg
      nsolwg_out = nsolwg
      do_annual_out = do_annual

end subroutine get_astronomy_for_swrad_init


!###################################################################

subroutine get_astronomy_for_swrad_3d (cosangsolar_out, fracday_out)

!------------------------------------------------------------------
      real, dimension(:, :,:), intent(out) :: cosangsolar_out
      real, dimension(:, :),   intent(out) ::  fracday_out

      fracday_out(:,:) = fracday(:,:)
      cosangsolar_out(:,:,:) = cosangsolar(:,:,:)

end subroutine get_astronomy_for_swrad_3d




!##################################################################

subroutine get_astronomy_for_swrad_2d (cosangsolar_out, fracday_out)

!------------------------------------------------------------------
      real, dimension(:, :),   intent(out) :: cosangsolar_out
      real, dimension(:, :),   intent(out) :: fracday_out

      fracday_out(:,:) = fracday(:,:)
      cosangsolar_out(:,:) = cosangsolar(:,:,1)

end subroutine get_astronomy_for_swrad_2d



!######################################################################

subroutine get_astronomy_for_clouds_init (nsolwg_out)

integer,              intent(out)  :: nsolwg_out

      nsolwg_out = nsolwg


end  subroutine get_astronomy_for_clouds_init


!##################################################################

subroutine get_astronomy_for_clouds (cosangsolar_out)

real, dimension(:,:,:), intent(out) :: cosangsolar_out


      cosangsolar_out(:,:,:) = cosangsolar(:,:,:)


end  subroutine get_astronomy_for_clouds


!######################################################################

subroutine return_cosz (cosz_in, fracday_in, rrin)

real, dimension(:,:), intent(in) :: cosz_in, fracday_in
real,                 intent(in) :: rrin
  
   cosangsolar(:,:,1) = cosz_in(:,:)
   fracday(:,:) = fracday_in(:,:)

   if (Environment%using_fms_periphs) then
     rrsun = rrin
   else if (Environment%using_sky_periphs) then
     rsun = rrin
   endif

end subroutine return_cosz 


!####################################################################


subroutine get_solar_distance  (rsun_out)

!---------------------------------------------------------------------
!     get_solar_distance returns the earth-sun distance.
!---------------------------------------------------------------------
real,             intent(out)  :: rsun_out

   if (Environment%using_fms_periphs) then
     rsun_out = rrsun
   else if (Environment%using_sky_periphs) then
     rsun_out = rsun
   endif

end subroutine get_solar_distance



!#####################################################################

subroutine zenitha (irad, fjulan)

!---------------------------------------------------------------------
!     zenitha computes the effective mean cosine of zenith angle and the
!     daylight fraction of the averaging period.
!     author: c. h. goldberg
!     revised: 2/20/93
!---------------------------------------------------------------------
 
integer,           intent(in) :: irad
real,              intent(in) :: fjulan
!-------------------------------------------------------------------

!----------------------------------------------------------------
!       arg     = integration period in radians.
!       denom   = cumulative angle of daylight integrated over.
!       gha     = hour angle of sun at greenwich (west of meridian is
!                 plus).
!       he      = MIN(next sunset, end of averaging period).
!       hlrise  = hour angle of sunrise.
!       hlset   = hour angle of sunset.
!       hlstart = hour angle of start of averaging period.
!       hlend   = hour angle of end of averaging period.
!       hs      = MIN(next sunrise, end of averaging period).
!       rnum    = cumulative contribution to integral of cos(h).
!---------------------------------------------------------------------

      real, dimension(:,:), allocatable :: denom, he, hlrise, hlset, &
                                           hlstart, hlend, rnum, dflag
      real                              :: dtsec, gha, arg
      integer                           :: i, j
!---------------------------------------------------------------------

      allocate ( denom   (IMINP:IMAXP, JMINP:JMAXP) )
      allocate ( he      (IMINP:IMAXP, JMINP:JMAXP) )
      allocate ( hlrise  (IMINP:IMAXP, JMINP:JMAXP) )
      allocate ( hlset   (IMINP:IMAXP, JMINP:JMAXP) )
      allocate ( hlstart (IMINP:IMAXP, JMINP:JMAXP) )
      allocate ( hlend   (IMINP:IMAXP, JMINP:JMAXP) )
      allocate ( rnum    (IMINP:IMAXP, JMINP:JMAXP) )
      allocate ( dflag   (IMINP:IMAXP, JMINP:JMAXP) )

!---------------------------------------------------------------------
!     calculations common to all grid points
!---------------------------------------------------------------------
      dtsec  = REAL(irad)
      gha    = fjulan*2.0E+00*pie + slag
      arg    = dtsec*2.0E+00*pie/secday

      do j=JMINP,JMAXP
        do i=IMINP,IMAXP

!---------------------------------------------------------------------
!     compute local angle of start and end of the averaging 
!     period.
!---------------------------------------------------------------------

          hlstart(i,j) = gha + rlong(i,j)
          hlend  (i,j) = hlstart(i,j) + arg

!---------------------------------------------------------------------
!     compute next sunrise and sunset after start of averaging period.
!     beware: MOD returns negative answers for negative arguments!
!     ensure that hlrise and hlset follow hlstart.
!---------------------------------------------------------------------

          hlrise(i,j) = hlstart(i,j) +    &
                        MOD(-hangle(i,j) - hlstart(i,j), 2.0E+00*pie)
          if (hlrise(i,j) .LT. hlstart(i,j))    &
            hlrise(i,j) = hlrise(i,j) + 2.0E+00*pie
            hlset (i,j) = hlstart(i,j) +      &
                          MOD(hangle(i,j) - hlstart(i,j), 2.0E+00*pie)
          if (hlset(i,j) .LT. hlstart(i,j))    &
            hlset (i,j) = hlset(i,j) + 2.0E+00*pie

!---------------------------------------------------------------------
!     if the averaging period starts with the sun up, then compute
!     integral over period to first sunset or end of averaging period,
!     whichever comes first.
!     note:  if hangle = pie then sun is always up.
!---------------------------------------------------------------------

          if ((hangle(i,j) .EQ. pie) .OR.     &
              (hlset(i,j) .LT. hlrise(i,j))) then
            he   (i,j) = MIN(hlset(i,j), hlend(i,j))
	    if ((hangle(i,j) .EQ. pie) .OR.     &
                ((hlset(i,j) .LT. hlrise(i,j)) .AND.    &
                 (he(i,j) .EQ. hlend(i,j))) ) then
              dflag(i,j) = 1.0
            else 
              dflag(i,j) = 0.0
            endif
            hlset(i,j) = hlset(i,j) + 2.0E+00*pie
            rnum (i,j) = SIN(he(i,j)) - SIN(hlstart(i,j))
            denom(i,j) = he(i,j) - hlstart(i,j)

!---------------------------------------------------------------------
!     if the averaging period starts in darkness, then there is no
!     contribution to the integral before the main loop.
!---------------------------------------------------------------------
          else
            rnum (i,j) = 0.0E+00
            denom(i,j) = 0.0E+00
	    dflag(i,j) = 0.0E+00
          endif
        end do
      end do

!---------------------------------------------------------------------
!     main loop:  in each iteration, integrate from next sunrise to
!                 next sunset or end of averaging period.
!---------------------------------------------------------------------
      do j=JMINP,JMAXP
        do i=IMINP,IMAXP
          do while (hlrise(i,j) .LT. hlend(i,j))
            he    (i,j) = MIN(hlset(i,j), hlend(i,j))
            rnum  (i,j) = rnum  (i,j) + SIN(he(i,j)) - SIN(hlrise(i,j))
            denom (i,j) = denom (i,j) + he(i,j) - hlrise(i,j)
            hlrise(i,j) = hlrise(i,j) + 2.0E+00*pie
            hlset (i,j) = hlset (i,j) + 2.0E+00*pie
          end do
        end do
      end do

!--------------------------------------------------------------------
!         denom(i,j) is the sum of the angles of daylight integrated
!         over.
!--------------------------------------------------------------------
      do j=JMINP,JMAXP
        do i=IMINP,IMAXP
          fracday(i,j)= denom(i,j)/arg

!--------------------------------------------------------------------
!         eliminate roundoff 1 as the daylight fraction -- set it to
!         exactly 1 when the sun is always up. similar problem does
!         not exist with the fracday .eq. 0 case.
!--------------------------------------------------------------------
	  if (dflag(i,j) .EQ. 1.0) then
	    fracday(i,j) = 1.0E+00
          endif
          if (denom(i,j) .GT. 0.0E+00) then
            zenith(i,j) = SIN(rlat(i,j))*SIN(declin) +    &
                     COS(rlat(i,j))*COS(declin) * rnum(i,j)/denom(i,j)
          else
            zenith(i,j) = 0.0E+00
          end if
        end do
      end do

!-------------------------------------------------------------------
   
      deallocate (denom    )
      deallocate (he       )
      deallocate (hlrise   )
      deallocate (hlset    )
      deallocate (hlstart  )
      deallocate (hlend    )
      deallocate (rnum     )
      deallocate (dflag    )

!---------------------------------------------------------------------

end  subroutine zenitha



!###################################################################

subroutine getastronomy  (iabsp, jabsp, lrad, irad, tdt, &
		          Rad_time_sv, lat_sv, do_average, fjulan, fjd)

!-----------------------------------------------------------------------
real,                           intent(in) ::   tdt
integer, dimension(:),          intent(in) ::   iabsp
integer, dimension(:),          intent(in) ::   jabsp
integer,                        intent(in) ::   irad
logical,                        intent(in) ::   lrad
real,                optional,  intent(in) ::   fjulan, fjd
logical,             optional,  intent(in) ::   do_average
real, dimension(:,:),optional,  intent(in) ::   lat_sv
type(time_type),     optional,  intent(in) ::   Rad_time_sv
!---------------------------------------------------------------------

      integer            :: j, i
      logical            :: do_avgt           
      type(time_type)    :: Dt_zen
      logical            :: do_stratchem

!--------------------------------------------------------------------
!     allocate space for zenith angle arrays and daylight fraction 
!     array. 
!--------------------------------------------------------------------
      allocate ( cosangsolar(IMINP:IMAXP, JMINP:JMAXP, nsolwg) )
      allocate ( zenith     (IMINP:IMAXP, JMINP:JMAXP        ) )
      allocate ( fracday    (IMINP:IMAXP, JMINP:JMAXP        ) )
      allocate ( ssolwd     (IMINP:IMAXP, JMINP:JMAXP        ) )

!-----------------------------------------------------------------------
!     if stratchem is activated, allocate space for its zenith angle.
!-----------------------------------------------------------------------
      call inquire_stratchem (do_stratchem)
      if (do_stratchem) then
        allocate ( szenij (IMINP:IMAXP, JMINP:JMAXP) )
      endif

!--------------------------------------------------------------------
!     allocate space for some locally used arrays. 
!--------------------------------------------------------------------
      allocate (hangle      (IMINP:IMAXP, JMINP:JMAXP) )
      allocate (rlat        (IMINP:IMAXP, JMINP:JMAXP) )
      allocate (rlong       (IMINP:IMAXP, JMINP:JMAXP) )

!-----------------------------------------------------------------------
!     define (i,j) arrays containing the hour angle and daylight frac-
!     tion and the grid point latitude and longitudes to be used to cal-
!     culate the zenith angle and the incoming solar radiation.
!-----------------------------------------------------------------------
      do j=JMINP,JMAXP
	do i=IMINP,IMAXP
          hangle(i,j)    = zhangle1(jabsp(j))
          fracday(i,j)   = zfracday1(jabsp(j))
	  rlat(i,j)      = zrlat(jabsp(j))
	  rlong(i,j)     = zrlong(iabsp(i), jabsp(j))
	end do
      end do

!-----------------------------------------------------------------------
!     define the zenith angle(s) needed for the desired radiation im-
!     plementation.
!-----------------------------------------------------------------------
      if (Environment%using_fms_periphs) then
        if (ldiurn) then
          if (do_average) then
            Dt_zen = set_time(int(tdt+0.5),0)
          else
            Dt_zen = set_time(irad,0)
          endif
          call diurnal_solar   (cosangsolar(:,:,1), ssolwd, lat_sv, &
		   		rlong, Rad_time_sv, Dt_zen, &
				fracday=fracday, rrsun=rrsun) 
        else if (do_annual) then
          call annual_mean_solar (cosangsolar(:,:,1), ssolwd,  &
			          lat_sv, fracday=fracday,   &
				  rrsun=rrsun)
        else
          call daily_mean_solar (cosangsolar(:,:,1), ssolwd,  &
			         lat_sv, Rad_time_sv,    &
				 fracday=fracday, rrsun=rrsun)
        endif
        fracday = MIN (fracday, 1.00)
      else if (Environment%using_sky_periphs) then
        if (ldiurn) then
          call zenitha  (irad, fjulan)
        else if (lswg) then
          call sol88zeng  
        else
          call sol88zen 
        endif
!---------------------------------------------------------------------
!    define an array containing these (this) zenith angle(s) returned
!    from zenitha, sol88zeng or sol88zen.  
!---------------------------------------------------------------------
        if (.not. lswg) then
          cosangsolar(:,:,1) = zenith(:,:)
        endif
      endif

!---------------------------------------------------------------------
!    determine whether averaged astronomical values are to be used in
!    the radiation calculation (only available in FMS).
!---------------------------------------------------------------------
      if (Environment%running_gcm .and. Environment%running_fms) then
        do_avgt = do_average
      else if (Environment%running_gcm .and.    &
		     Environment%running_skyhi) then
        do_avgt = .false.
      else if (Environment%running_standalone) then
        do_avgt = .false.
      endif

!---------------------------------------------------------------------
!    send astronomy data to radiation_diag_mod if the data is to be
!    used and not averaged.  
!---------------------------------------------------------------------
      if (lrad .and. .not. do_avgt) then
        do j=JMINP,JMAXP
          if (Rad_control%do_raddg(jabsp(j)  )) then
            call radiag_from_astronomy(cosangsolar, fracday, nsolwg, & 
                                       jabsp(j), j, iabsp(IMINP),    &
                                       iabsp(IMAXP))
          endif
        end do
      endif

!-----------------------------------------------------------------------
!     if stratchem is activated, determine the zenith angle for the 
!     chemistry calculations.
!-----------------------------------------------------------------------
      if (do_stratchem) then
        call zenithb (tdt, fjd)
      endif

!------------------------------------------------------------------
!     deallocate local arrays
!------------------------------------------------------------------
      deallocate (hangle)
      deallocate (rlong)
      deallocate (rlat)


end subroutine getastronomy




!###################################################################

subroutine sol88 (fjd, jld)
!--------------------------------------------------------------------
!     sol88 computes radius vector, declination and right ascension of
!     sun, equation of time, and hour angle of sun at sunset, and 
!     daylight fraction for nlat specified latitudes given the current 
!     julian day and fraction.
!     perihelion-to-equinox days, deleqn, are calculated more directly
!     from ayearfd and tyearfd and t.  proposed, but commented-out 
!     formula takes dependence of ayearfd and tyearfd on t into
!     account.
!     author: c. h. goldberg
!     revised: 1/1/93
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!       fjd    = fractional part of current julian day.
!       jld    = integral part of current julian day.
!--------------------------------------------------------------------
real,                  intent(in)  :: fjd
integer,               intent(in)  :: jld

!--------------------------------------------------------------------
!     jdor  = jld of epoch which is january 0, 1900 at 12 hours ut.
!     tpp   = days between epoch and perihelion passage of 1900.
!     svt6  = days between perihelion passage and march equinox
!             of 1900.
!     small = convergence criterion.
!--------------------------------------------------------------------

      integer :: jdor= 2415020
      real    :: tpp=1.5500E+00
      real    :: svt6=7.8035E+01
      real    :: small=1.0E-07
      real    :: ayearfd0=0.25964134E+00
      real    :: ayearfd1=0.304E-05
      real    :: tyearfd0=0.24219879E+00
      real    :: tyearfd1=-0.614E-05
      real    :: dat, t, ayearfd, tyearfd, ec, angin, angincr, &
                 ador, deleqn, ayear, er, qq, e0, e, diff, he, eq, &
		 date, daterad, w, sind, alp, tst, dateang
      integer :: jdoe

!-------------------------------------------------------------------

      dat = jld - jdor - tpp + fjd

!-------------------------------------------------------------------
!     compute time in julian centuries after epoch.
!--------------------------------------------------------------------
      t = (jld - jdor)/(1.0E+02*calyear)

!--------------------------------------------------------------------
!     compute length of anomalistic and tropical ayearfds (minus 365
!     days).
!--------------------------------------------------------------------
      ayearfd = ayearfd0 + ayearfd1*t
      tyearfd = tyearfd0 + tyearfd1*t

!--------------------------------------------------------------------
!     compute orbit eccentricity and angle of earth's inclination
!     from julian centuries ater epoch.
!--------------------------------------------------------------------
      ec    = 0.1675104E-01 - (0.418000E-04 + 0.126E-06*t)*t
      angin = 2.3452294E+01 - (0.130125E-01 + 0.164E-05*t)*t
      angincr = angin / radians_to_degrees

!--------------------------------------------------------------------
!     deleqn=updated svt6 for current date. (days to equinox)
!--------------------------------------------------------------------
!     Three formulae appear here:
!      (1) the original formula.
!      (2) a clearer, mathematically equivalent formula
!           [this results in some last digit differences in the arrays
!            calculated in Sol, but not in the overall answers.]
!      (3) a better formula that integrates the linear expression
!          for anomalous year and tropical year lengths from the epoch
!          to the given day instead of using their values only on the
!          given day.
!           [this results in calculated values that differ in the 5th
!            significant digit here, and .... in the final answers.]
!--------------------------------------------------------------------
!     original formula.
!--------------------------------------------------------------------
      ador = jdor
      jdoe = ador + (svt6*calyear)/(ayearfd - tyearfd)
      deleqn = (jdoe -jld)*(ayearfd - tyearfd)/calyear

!--------------------------------------------------------------------
!     equivalent formula.
!--------------------------------------------------------------------
!     deleqn = svt6 - t*1.0E+02*(ayearfd - tyearfd)
!--------------------------------------------------------------------
!     better formula.
!     deleqn = svt6 - t*1.0E+02*
!    &        (ayearfd0 - tyearfd0 + 0.5E+00*t*(ayearfd1-tyearfd1))
!--------------------------------------------------------------------
      ayear  = ayearfd + 3.65E+02
      er     = SQRT((1.0E+00 + ec)/(1.E+00 - ec))

!--------------------------------------------------------------------
!     determine true anomaly at equinox.
!     solve:  e - ec*SIN(e) = deleqn*(2*pi/ayear) by newton's method
!---------------------------------------------------------------------
      qq = deleqn*2.0E+00*pie/ayear
      e0 = 1.0E+00

!---------------------------------------------------------------------
!     repeat until ABS(e-e0) <= small
!---------------------------------------------------------------------
101   continue
      e    = e0 - (e0 - ec*SIN(e0) - qq)/(1.0E+00 - ec*COS(e0))
      diff = ABS(e-e0)
      e0   = e
      if (diff .GT. small) goto 101
      he = 0.5E+00*e
      eq = 2.0E+00*ATAN(er*TAN(he))

!--------------------------------------------------------------------
!     solve orbit equations by Newton's method
!           e - ec*SIN(e) = date*(2*pi/ayear)
!     date=days since last perihelion passage.
!--------------------------------------------------------------------
      date = MOD(dat, ayear)
      daterad   = 2.0E+00*pie*date/ayear
      e0   = 1.0E+00

!--------------------------------------------------------------------
!     repeat until ABS(e-e0) <= small
!--------------------------------------------------------------------
111   continue
         e = e0 - (e0 - ec*SIN(e0) -daterad)/(1.0E+00 - ec*COS(e0))
	 diff = ABS(e-e0)
	 e0 = e
	 if (diff .GT. small) goto 111

!--------------------------------------------------------------------
!     define the earth-sun distance. 
!--------------------------------------------------------------------
      rsun   = 1.0E+00 - ec*COS(e)

!--------------------------------------------------------------------
!     define the declination of the earth.
!--------------------------------------------------------------------
      he     = 0.5E+00*e
      w      = 2.0E+00*ATAN(er*TAN(he))
      sind   = SIN(angincr)*SIN(w - eq)
      declin = ASIN(sind)

!--------------------------------------------------------------------
!     define the right ascension of the sun.
!--------------------------------------------------------------------
      if (do_oldrightasc) then
        alp    = ASIN(TAN(declin)/TAN(angincr))
        tst  = COS(w - eq)
        if(tst .LT. 0.0E+00) alp = pie - alp
        if(alp .LT. 0.0E+00) alp = alp + 2.0E+00*pie
        rightasc = alp
      else
        rightasc = ACOS(COS(w-eq)/COS(declin))
        if (SIN(w - eq) .LT. 0.0E+00) rightasc = 2.0E+00*pie - rightasc
      endif

      dateang  = 2.0E+00*pie*(date - deleqn)/ayear
      if(dateang .LT. 0.0E+00) then
        dateang = dateang + 2.0E+00*pie
      endif

!---------------------------------------------------------------------
!     define the apparent sun lag angle.    
!---------------------------------------------------------------------
      slag = dateang - rightasc - 0.3255E-01

!---------------------------------------------------------------------



end subroutine sol88


!####################################################################

subroutine sol88zen 

!--------------------------------------------------------------------
!     sol88zen computes the zenith cosines at each latitude, averaged
!     over the complete period of daylight containing the julian day
!     and fractional day specified.
!     author: c. h. goldberg
!     revised: 1/1/93
!--------------------------------------------------------------------

      real, dimension(:, :), allocatable   ::  cc, ss
      integer                              ::  i, j

!--------------------------------------------------------------------
      allocate ( ss (IMINP:IMAXP, JMINP:JMAXP) )
      allocate ( cc (IMINP:IMAXP, JMINP:JMAXP) )

!--------------------------------------------------------------------
      do j=JMINP,JMAXP
	do i=IMINP,IMAXP
          if(hangle(i,j) .EQ. 0.0E+00) then
  	    zenith(i,j) = 0.0E+00
          else
	    ss(i,j) = SIN(rlat(i,j))*SIN(declin)
	    cc(i,j) = COS(rlat(i,j))*COS(declin)
            zenith(i,j) = ss(i,j) + cc(i,j)*SIN(hangle(i,j))/hangle(i,j)
          endif
	enddo
      enddo

!---------------------------------------------------------------------
!  deallocate local arrays
!---------------------------------------------------------------------
      deallocate (ss)
      deallocate (cc)


end subroutine sol88zen


!####################################################################

subroutine sol88zeng

!-------------------------------------------------------------------
!     sol88zeng computes the zenith cosines at each latitude at 
!     gaussian times over the half day from local noon to sunset.
!     author: c. h. goldberg
!     revised: 1/1/93
!--------------------------------------------------------------------

      real, dimension(:,:), allocatable :: cc, ss
      real, dimension(:,:,:), allocatable :: hmu
      integer                             :: i, j, ngp

!---------------------------------------------------------------------
      allocate (hmu  (IMINP:IMAXP, JMINP:JMAXP, nsolwg) )
      allocate (ss   (IMINP:IMAXP, JMINP:JMAXP        ) )
      allocate (cc   (IMINP:IMAXP, JMINP:JMAXP        ) )

!---------------------------------------------------------------------
!     calculation of gaussian coszwgs.
!     calculation of gaussian cosangsolars.
!---------------------------------------------------------------------
      do j=JMINP,JMAXP
	do i=IMINP,IMAXP
	  if (hangle(i,j) .EQ. 0.0E+00) then
	    zenith(i,j) = 0.0E+00
          else
	    ss(i,j) = SIN(rlat(i,j))*SIN(declin)
	    cc(i,j) = COS(rlat(i,j))*COS(declin)
            zenith(i,j) = ss(i,j) + cc(i,j)*SIN(hangle(i,j))/hangle(i,j)
	  endif
        enddo
      enddo

      do ngp=1,nsolwg
        do j=JMINP,JMAXP
	  do i=IMINP,IMAXP
	    if (hangle(i,j) .EQ. 0.0E+00) then
	        cosangsolar(i,j,ngp) = 0.0E+00
            else
              hmu(i,j,ngp) = gpt(ngp)*hangle(i,j)
	      cosangsolar(i,j,ngp) =     &
                   MAX(ss(i,j) + cc(i,j)*COS(hmu(i,j,ngp)), 0.0E+00)
	    endif
          enddo
        enddo
      enddo

!---------------------------------------------------------------------
!  deallocate local arrays
!---------------------------------------------------------------------
      deallocate (hmu)
      deallocate (ss )
      deallocate (cc )

!-------------------------------------------------------------------

end  subroutine sol88zeng


!#####################################################################

      subroutine zenithb (tdt, fjd)

!---------------------------------------------------------------------
!     zenithb computes the effective mean cosine of zenith angle over
!     the given model time step.
!---------------------------------------------------------------------

real,           intent(in) :: tdt, fjd

!---------------------------------------------------------------------
      real, dimension(:,:), allocatable :: hloc, sinsin, coscos
      integer                           :: i, j
      real                              :: gha, arg, piet2

!-----------------------------------------------------------------
!  allocate local arrays.
!-----------------------------------------------------------------
      allocate (hloc  (IMINP:IMAXP, JMINP:JMAXP) )
      allocate (sinsin(IMINP:IMAXP, JMINP:JMAXP) )
      allocate (coscos(IMINP:IMAXP, JMINP:JMAXP) )

!---------------------------------------------------------------------
!     define the increase in time from the mid time step to the end of 
!     the time step ( delta t/2, in radians, where 1 day is 2 pi 
!     radians). 
!---------------------------------------------------------------------
      piet2 = 2.0E+00*pie
      arg = 0.5E+00*tdt/3600.E+00*(piet2/24.0E+00)

!---------------------------------------------------------------------
!     define the local hour angle at the end of the time step.
!---------------------------------------------------------------------
      gha = fjd*piet2 + slag + arg

!---------------------------------------------------------------------
!     define the local zenith angle. if the current values of declin 
!     and slag are desired, call Astronomy which will call Sol88 to 
!     generate these values before calling Getastronomy which calls 
!     zenithb. otherwise, the value of declin from the last radiation 
!     step is used.
!---------------------------------------------------------------------
      do j=JMINP,JMAXP
        do i=IMINP,IMAXP
          sinsin(i,j) = SIN(rlat(i,j))*SIN(declin)
          coscos(i,j) = COS(rlat(i,j))*COS(declin)
          hloc(i,j) = gha + rlong(i,j)
          hloc(i,j) = MOD(hloc(i,j), piet2)
          if (hloc(i,j) .GT. pie) then
            hloc(i, j) = hloc(i,j) - piet2
          endif
          szenij(i,j) = ACOS(sinsin(i,j) + coscos(i,j)*COS(hloc(i,j)))
        end do
      end do

!-------------------------------------------------------------------
!  deallocate local arrays.
!-------------------------------------------------------------------
      deallocate (hloc   )
      deallocate (sinsin )
      deallocate (coscos )

!--------------------------------------------------------------------



end  subroutine zenithb



!####################################################################

subroutine gptinit

!-----------------------------------------------------------------------
!     gptinit initializes data for gaussian quadrature.
!     author: c. h. goldberg
!     revised: 1/1/93
!--------------------------------------------------------------------

      real, dimension(1) ::  gpt1
      real, dimension(2) ::  gpt2
      real, dimension(4) ::  gpt4
      real, dimension(8) ::  gpt8

!-----------------------------------------------------------------------
!     define gaussian quadrature points for one, two, four and eight
!     point quadratures.
!-----------------------------------------------------------------------

      data gpt1 /0.5000000000E+00/
      data gpt2 /0.2113248654E+00, 0.7886751346E+00/
      data gpt4 /0.0694318442E+00, 0.3300094782E+00,   &
                 0.6699905218E+00, 0.9305681558E+00/
      data gpt8/ 0.0198550718E+00, 0.1016667613E+00,   &
                 0.2372337950E+00, 0.4082826788E+00,   &
                 0.5917173212E+00, 0.7627662050E+00,   &
                 0.8983332387E+00, 0.9801449282E+00/

!----------------------------------------------------------------------
!     fill array with desired quadrature points. 
!----------------------------------------------------------------------

      allocate ( gpt (nsolwg)  )

      if (nsolwg == 1) then
	gpt(:) = gpt1(:)
      else if (nsolwg == 2) then
	gpt(:) = gpt2(:)
      else if (nsolwg == 4) then
	gpt(:) = gpt4(:)
      else if (nsolwg == 8) then
	gpt(:) = gpt8(:)
      endif
!------------------------------------------------------------------
	
	

end  subroutine gptinit



!####################################################################



                   end module astronomy_package_mod


