		module std_pressures_mod

use  utilities_mod,       only: open_file, file_exist,   & 
			        check_nml_error, error_mesg, &
			        print_version_number, FATAL, NOTE, &
				WARNING, get_my_pe, close_file
use rad_utilities_mod,    only: Environment, environment_type
use constants_new_mod,    only: rgas, rearth, grav, wtmair, pstd

!---------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!               module to hold standard pressure levels used 
!                       by radiation package modules
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: std_pressures.F90,v 1.2 2001/08/30 15:12:37 fms Exp $'
   character(len=128)  :: tag     =  '$Name: eugene $'



!---------------------------------------------------------------------
!-------  interfaces --------

public  std_pressures_init,     &
	skyp_std, sigp_std, lanz_std,       &
	set_std_pressures, get_std_pressures

!---------------------------------------------------------------------
!-------- namelist  ---------

integer                    :: LUS = 0


namelist / std_pressures_nml /            &
	         	          LUS


!---------------------------------------------------------------------
!------ public data -----





!---------------------------------------------------------------------
!------ private data -----


real, allocatable, dimension (:)      :: q, qi, qmh
real, allocatable, dimension (:)      :: zh, zm, hm
real, allocatable, dimension (:)      :: pd, pd8, plm, plm8

real, dimension(40)             :: q40
real, dimension(41)             :: q40h
real, dimension(80)             :: q80
real, dimension(81)             :: q80h
real, dimension(9)              :: q9e
real, dimension(9)              :: q9g
real, dimension(18)             :: q18e
real, dimension(18)             :: q18m
real, dimension(12)             :: q12
real, dimension(30)             :: q30
real, dimension(14)             :: q14

real, dimension(46)             :: z45h
real, dimension(45)             :: z45
real, dimension(7)              :: dtdz
real, dimension(8)              :: hlyr_top
real, dimension(8)              :: plyr_top
real, dimension(8)              :: tlyr_top
real                            :: plyr_top_2km

integer                         :: k
character(len=10)               :: type_form
real                            :: pstd_08

!----------------------------------------------------------------------
!      specified values for some cases follow
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!                          skyhi l40
!---------------------------------------------------------------------

data q40 /                                                       &
     0.9447012E-05, 0.3042544E-04, 0.6876343E-04, 0.1295448E-03, &
     0.2187045E-03, 0.3446506E-03, 0.5182561E-03, 0.7525096E-03, &
     0.1065881E-02, 0.1484524E-02, 0.2044372E-02, 0.2791420E-02, &
     0.3783706E-02, 0.5096454E-02, 0.6827019E-02, 0.9101167E-02, &
     0.1208109E-01, 0.1597555E-01, 0.2105252E-01, 0.2765502E-01, &
     0.3622044E-01, 0.4730383E-01, 0.6154377E-01, 0.7960552E-01, &
     0.1022674E+00, 0.1304765E+00, 0.1652889E+00, 0.2078405E+00, &
     0.2584953E+00, 0.3166243E+00, 0.3816890E+00, 0.4530773E+00, &
     0.5297009E+00, 0.6099063E+00, 0.6913949E+00, 0.7711706E+00, &
     0.8455313E+00, 0.9101312E+00, 0.9601436E+00, 0.9905502E+00/
 
data q40h /                                                      &
     0.0          , 0.1911410E-04, 0.4843058E-04, 0.9763264E-04, &
     0.1718878E-03, 0.2782727E-03, 0.4268617E-03, 0.6292185E-03, &
     0.8999584E-03, 0.1262393E-02, 0.1745741E-02, 0.2394086E-02, &
     0.3254696E-02, 0.4398700E-02, 0.5904891E-02, 0.7893145E-02, &
     0.1049407E-01, 0.1390813E-01, 0.1835028E-01, 0.2415268E-01, &
     0.3166522E-01, 0.4143094E-01, 0.5400919E-01, 0.7012945E-01, &
     0.9036189E-01, 0.1157414E+00, 0.1470876E+00, 0.1857426E+00, &
     0.2325674E+00, 0.2873136E+00, 0.3489250E+00, 0.4175294E+00, &
     0.4916515E+00, 0.5706949E+00, 0.6518117E+00, 0.7333819E+00, &
     0.8109064E+00, 0.8816345E+00, 0.9395488E+00, 0.9811897E+00, &
     0.1000000E+01/

!---------------------------------------------------------------------
!                           skyhi l80
!---------------------------------------------------------------------

data (q80(k),k=1,40) /                                           &
     3.0048389E-07, 7.7225698E-07, 1.6748142E-06, 3.2701594E-06, &
     5.9000549E-06, 1.0061819E-05, 1.6268606E-05, 2.4933648E-05, &
     3.6552699E-05, 5.1701686E-05, 7.1004114E-05, 9.5131055E-05, &
     1.2480116E-04, 1.6078072E-04, 2.0392604E-04, 2.5549930E-04, &
     3.1701189E-04, 3.8987631E-04, 4.7564623E-04, 5.7602575E-04, &
     6.9316870E-04, 8.3003441E-04, 9.8990388E-04, 1.1760961E-03, &
     1.3933353E-03, 1.6476632E-03, 1.9453214E-03, 2.2933711E-03, &
     2.6999937E-03, 3.1746642E-03, 3.7283529E-03, 4.3737555E-03, &
     5.1255595E-03, 6.0007491E-03, 7.0189556E-03, 8.2028578E-03, &
     9.5786412E-03, 1.1176520E-02, 1.3031334E-02, 1.5183222E-02/
 
data (q80(k),k=41,80) /                                          &
     1.7678386E-02, 2.0569954E-02, 2.3918942E-02, 2.7795333E-02, &
     3.2279265E-02, 3.7462350E-02, 4.3449095E-02, 5.0358442E-02, &
     5.8306312E-02, 6.7394856E-02, 7.7744210E-02, 8.9501313E-02, &
     1.0282456E-01, 1.1788332E-01, 1.3485707E-01, 1.5393386E-01, &
     1.7530822E-01, 1.9917813E-01, 2.2571259E-01, 2.5473484E-01, &
     2.8597155E-01, 3.1939434E-01, 3.5494172E-01, 3.9251494E-01, &
     4.3197377E-01, 4.7313236E-01, 5.1575521E-01, 5.5955356E-01, &
     6.0418221E-01, 6.4923710E-01, 6.9425381E-01, 7.3870736E-01, &
     7.8201346E-01, 8.2353171E-01, 8.6257094E-01, 8.9839709E-01, &
     9.3024409E-01, 9.5732772E-01, 9.7886309E-01, 9.9408549E-01/

data (q80h(k),k=1,40) /                                          &
     0.0,           5.0947389E-07, 1.1705817E-06, 2.3962467E-06, &
     4.4627886E-06, 7.8002008E-06, 1.2979178E-05, 2.0391704E-05, &
     3.0487240E-05, 4.3824885E-05, 6.0994211E-05, 8.2656765E-05, &
     1.0948792E-04, 1.4225615E-04, 1.8171755E-04, 2.2884873E-04, &
     2.8525346E-04, 3.5230612E-04, 4.3145301E-04, 5.2436611E-04, &
     6.3277482E-04, 7.5932675E-04, 9.0732629E-04, 1.0799970E-03, &
     1.2807462E-03, 1.5158220E-03, 1.7909715E-03, 2.1129734E-03, &
     2.4891705E-03, 2.9286727E-03, 3.4413176E-03, 4.0393293E-03, &
     4.7358696E-03, 5.5473149E-03, 6.4912468E-03, 7.5895646E-03, &
     8.8657096E-03, 1.0348903E-02, 1.2070324E-02, 1.4068858E-02/
 
data (q80h(k),k=41,81) /                                         &
     1.6385852E-02, 1.9072877E-02, 2.2184541E-02, 2.5788941E-02, &
     2.9957823E-02, 3.4780597E-02, 4.0350880E-02, 4.6785197E-02, &
     5.4204596E-02, 6.2718409E-02, 7.2419991E-02, 8.3459859E-02, &
     9.5980093E-02, 1.1015711E-01, 1.2615144E-01, 1.4416346E-01, &
     1.6436642E-01, 1.8697840E-01, 2.1217386E-01, 2.4011521E-01, &
     2.7024459E-01, 3.0261376E-01, 3.3710545E-01, 3.7372171E-01, &
     4.1225322E-01, 4.5263768E-01, 4.9455500E-01, 5.3786422E-01, &
     5.8211753E-01, 6.2708324E-01, 6.7217362E-01, 7.1705931E-01, &
     7.6100896E-01, 8.0359771E-01, 8.4396020E-01, 8.8159207E-01, &
     9.1552246E-01, 9.4520244E-01, 9.6960855E-01, 9.8820596E-01, &
     1.0/

!-------------------------------------------------------------------
!                   9 level sigmas for exp pred gp
!--------------------------------------------------------------------

data q9e/.0089163236,.074074073,.18861453,.33607680,.49999997,     &
       .66392314,.81138542,.92592590,.99108367/
 
!-----------------------------------------------------------------
!                   9 level sigmas for manabe model
!-----------------------------------------------------------------

data q9g/.025,.095,.205,.350,.515,.680,.830,.940,.990/
 
!-----------------------------------------------------------------
!                14 level sigmas for wetherald r30l14 model
!-----------------------------------------------------------------

data q14 /.0150,.05035,.1009,.17065,.2569,.3549,.4600,.5682,     &
             .6755,.77695,.86605,.9353,.97865,.99665/
 
!-----------------------------------------------------------------
!                  30 level sigmas for manabe model
!-----------------------------------------------------------------

data q30/                                                        &
       .00334,.01534,.02963,.04656,.06649,.08979,.11680,.14781,   &
       .18301,.22246,.26606,.31348,.36419,.41744,.47230,.52770,   &
       .58256,.63581,.68652,.73394,.77754,.81699,.85219,.88320,   &
       .91021,.93351,.95344,.97037,.98466,.99666/
 
!-----------------------------------------------------------------
!                  18 level sigmas for exp pred gp
!-----------------------------------------------------------------

data q18e/.0022719,.0196759,.0525120,.0987226,.1562500,           &
       .2230367,.2970250,.3761574,.4583762,.5416238,.6238425,   &
       .7029749,.7769632,.8437499,.9012774,.9474880,.9803240,   &
       .9977280/
 
!-----------------------------------------------------------------
!                  18 level sigmas for nmc mrf model
!-----------------------------------------------------------------

data q18m/                                                         &
       .0207469, .0739862, .1244004, .1745733, .2246687, .2747291,   &
       .3247711, .3748014, .4248250, .4974484, .5935378, .6881255,   &
       .7772229, .8563145, .9204018, .9604809, .9814907, .9949968/
 
!-----------------------------------------------------------------
!                12 level sigmas for ncar ccm model
!-----------------------------------------------------------------

data q12/.009,.025,.060,.110,.165,.245,.355,.500,.664,               &
       .811,.926,.991/

!-------------------------------------------------------------------
!    dtdz is the temperature gradient of the US Standard Atmospheres
!    (1977). this quantity is in (degrees/geopotential km).
!-------------------------------------------------------------------

data dtdz / -6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0/

!-------------------------------------------------------------------
!    hlyr_top gives altitude boundaries (in geopotential km) of 
!    regions in the US Standard Atmospheres with different temperature
!    gradients.
!-------------------------------------------------------------------

data hlyr_top /                                          &
         0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.8520/
 
!-------------------------------------------------------------------
!    plyr_top gives pressures (in mb) tabulated in the US Standard
!   Atmospheres for the altitudes given in hlyr_top.
!-------------------------------------------------------------------

data plyr_top /                                          &
       1.01325E+03, 2.2632E+02, 5.4748E+01, 8.6801E+00,        &
       1.1090E+00,  6.6938E-01, 3.9564E-02, 3.7338E-03/
 
!-------------------------------------------------------------------
!    plyr_top_2km gives the pressure (in mb) of the US Standard
!   Atmospheres profile at 2 km
!-------------------------------------------------------------------

data plyr_top_2km /                                          &
       7.9495E+02/
 
!-------------------------------------------------------------------
!    tlyr_top gives temperatures (in K) tabulated in the US Standard
!   Atmospheres for the altitudes given in hlyr_top.
!-------------------------------------------------------------------

data tlyr_top /                                          &
       288.150, 216.650, 216.650, 228.650, 270.650, 270.650,   &
       214.650, 186.870/
 
!-------------------------------------------------------------------
!     z45h are data for layer boundaries of the 45 layer version of
!     the lan model. one (dummy) layer is added between the highest
!     model layer and the top of the atmosphere.
!                 (in meters, first value at the surface)
!-------------------------------------------------------------------

data z45h/                                                    &
     0.0000000E+00, 5.0000000E+02, 1.0000000E+03, 1.5000000E+03,   &
     2.0000000E+03, 2.5000000E+03, 3.0000000E+03, 3.5000000E+03,   &
     4.0000000E+03, 4.5000000E+03, 5.0000000E+03, 5.5000000E+03,   &
     6.0000000E+03, 6.5000000E+03, 7.0000000E+03, 7.5000000E+03,   &
     8.0000000E+03, 8.5000000E+03, 9.0000000E+03, 9.5000000E+03,   &
     1.0000000E+04, 1.0500000E+04, 1.1000000E+04, 1.1500000E+04,   &
     1.2000000E+04, 1.2500000E+04, 1.3000000E+04, 1.3500000E+04,   &
     1.4000000E+04, 1.4500000E+04, 1.5000000E+04, 1.5500000E+04,   &
     1.6000000E+04, 1.6500000E+04, 1.7000000E+04, 1.7500000E+04,   &
     1.8000000E+04, 1.8500000E+04, 1.9000000E+04, 1.9500000E+04,   &
     2.0000000E+04, 2.0500000E+04, 2.1000000E+04, 2.1500000E+04,   &
     2.2000000E+04, 1.0000000E+05/
 
!-------------------------------------------------------------------
!     z45 are data for layer-mean altitudes of the 45 layer version of
!     the lan model. one (dummy) layer is added between the highest
!     model layer and the top of the atmosphere.
!                 (in meters, first value at the surface)
!-------------------------------------------------------------------

data z45/                                                  &
     2.5000000E+02, 7.5000000E+02, 1.2500000E+03, 1.7500000E+03,   &
     2.2500000E+03, 2.7500000E+03, 3.2500000E+03, 3.7500000E+03,   &
     4.2500000E+03, 4.7500000E+03, 5.2500000E+03, 5.7500000E+03,   &
     6.2500000E+03, 6.7500000E+03, 7.2500000E+03, 7.7500000E+03,   &
     8.2500000E+03, 8.7500000E+03, 9.2500000E+03, 9.7500000E+03,   &
     1.0250000E+04, 1.0750000E+04, 1.1250000E+04, 1.1750000E+04,   &
     1.2250000E+04, 1.2750000E+04, 1.3250000E+04, 1.3750000E+04,   &
     1.4250000E+04, 1.4750000E+04, 1.5250000E+04, 1.5750000E+04,   &
     1.6250000E+04, 1.6750000E+04, 1.7250000E+04, 1.7750000E+04,   &
     1.8250000E+04, 1.8750000E+04, 1.9250000E+04, 1.9750000E+04,   &
     2.0250000E+04, 2.0750000E+04, 2.1250000E+04, 2.1750000E+04,   &
     2.2250000E+04/

!--------------------------------------------------------------------

integer :: kmin, kmax



!--------------------------------------------------------------------
!--------------------------------------------------------------------




                         contains




subroutine std_pressures_init (kmin_in, kmax_in)

integer, intent(in)    :: kmin_in, kmax_in

!--------------------------------------------------------------------
      integer            :: unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=std_pressures_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'std_pressures_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'std_pressures', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=std_pressures_nml)
      endif
      call close_file (unit)

      kmin = kmin_in
      kmax = kmax_in

      pstd_08 = 0.8*pstd

      allocate (pd   (kmin:kmax+1))
      allocate (pd8  (kmin:kmax+1))
      allocate (plm  (kmin:kmax+1))
      allocate (plm8 (kmin:kmax+1))

!-------------------------------------------------------------------


end subroutine std_pressures_init




!####################################################################

subroutine skyp_std(qlevel_in, qmh_in)
 
!-------------------------------------------------------------------
real, dimension (:), intent(in), optional     :: qlevel_in, qmh_in
!-------------------------------------------------------------------

!--------------------------------------------------------------------
!      The general formulas for skyhi pressures are: 
!
!             for the upper sandwich-k=1,LUS-----
!
!                qi(i,k)=pstd*q(k)/pss(i)
!
!          for the lower sandwich k= LUP,KMAX
!
!    qi(i,k)=pc*(1.-q(k))+pss(i)*(q(k)-qmh(LUP))
!            ---------------------------------
!             pss(i)*(1.-qmh(LUP))
!
!    where pc=pstd*qmh(LUP),pss(i)=sfc pressure at point i,
!    and q,qmh are specified levels
!--------------------------------------------------------------------
      integer    :: LUP, k
      real       :: pss


      LUP = LUS + 1
!--------------------------------------------------------------------
!  allocate temporary arrays
!--------------------------------------------------------------------
      allocate ( q(1:kmax))
      allocate ( qi(1:kmax))
      allocate ( qmh(1:kmax+1))

!---------------------------------------------------------------------
!    define sigmas (q,qmh) for the appropriate skyhi model structure
!---------------------------------------------------------------------
      if (present(qlevel_in)) then
        do k=1,kmax
	  q(k) = qlevel_in(k)
        enddo
        do k=1,kmax+1
	  qmh(k) = qmh_in(k)
        enddo
      else
        type_form = Environment%column_type
        if (trim(type_form) == 'skyl40') then
	  if (kmax == 40) then
            do k=1,kmax
              q(k) = q40(k)
            enddo
            do k=1,kmax+1
              qmh(k) = q40h(k)
            enddo
	  else
	    call error_mesg ( 'skyp_std', &
	      'KMAX of model disagrees with designation selected', &
							       FATAL)
	  endif
        else if (trim(type_form) == 'skyl80') then
	  if (KMAX == 80) then
            do k=1,KMAX
              q(k) = q80(k)
            enddo
            do k=1,KMAX+1
              qmh(k) = q80h(k)
            enddo
	  else
	    call error_mesg ( 'skyp_std', &
	      'KMAX of model disagrees with designation selected', &
							       FATAL)
	  endif
        endif
      endif

!---------------------------------------------------------------------
!    for the case pss=pstd (standard pressure) both formulas
!     reduce to qi(i,k)=q(k) for all k.
!---------------------------------------------------------------------
      pss= pstd

!---------------------------------------------------------------------
!   define layer-mean pressures (pd) from the appropriate q's and pss
!---------------------------------------------------------------------
      do k=1, LUS
         qi(k)=pstd*q(k)/pss
         pd(k)=qi(k)*pss
      enddo
      do k= LUP,KMAX
         qi(k)=(pstd*qmh(LUP)*(1.-q(k))+pss*(q(k)-qmh(LUP)))/   &
                   (pss*(1.-qmh(LUP)))
         pd(k)=qi(k)*pss
      enddo
      pd(KMAX+1)=pss
 
!---------------------------------------------------------------------
!   define all layer-boundary pressures (plm) as the mean of adjacent
!   layer-mean pressures. this is necessary for co2 interpolation.
!---------------------------------------------------------------------
      plm(1)=0.
      do k=2,KMAX
        plm(k)=0.5*(pd(k-1)+pd(k))
      enddo
      plm(KMAX+1)=pss
!---------------------------------------------------------------------
!  the following puts p-data into mb
!---------------------------------------------------------------------
      do k=1,KMAX+1
        pd(k)=pd(k)*1.0e-3
        plm(k)=plm(k)*1.0e-3
      enddo

!---------------------------------------------------------------------
!***second pass: pss = ~810mb
!---------------------------------------------------------------------
      pss = pstd_08

!---------------------------------------------------------------------
!   define layer-mean pressures (pd) from the appropriate q's and pss
!---------------------------------------------------------------------
      do k=1, LUS
        qi(k)=pstd*q(k)/pss
        pd8(k)=qi(k)*pss
      enddo
      do k= LUP,KMAX
        qi(k)=(pstd*qmh(LUP)*(1.-q(k))+pss*(q(k)-qmh(LUP)))/  &
                  (pss*(1.-qmh(LUP)))
        pd8(k)=qi(k)*pss
      enddo
      pd8(KMAX+1)=pss
 
!---------------------------------------------------------------------
!   define all layer-boundary pressures (plm) as the mean of adjacent
!   layer-mean pressures. this is necessary for co2 interpolation.
!---------------------------------------------------------------------
      plm8(1)=0.
      do k=1,KMAX-1
        plm8(k+1)=0.5*(pd8(k)+pd8(k+1))
      enddo
      plm8(KMAX+1)=pss

!---------------------------------------------------------------------
!  the following puts p-data into mb
!---------------------------------------------------------------------
      do k=1,KMAX+1
        pd8(k)=pd8(k)*1.0e-3
        plm8(k)=plm8(k)*1.0e-3
      enddo

      deallocate (q)
      deallocate (qi)
      deallocate (qmh)
!-------------------------------------------------------------------



end subroutine skyp_std




!####################################################################

subroutine sigp_std(qlevel_in, qmh_in)
 
!--------------------------------------------------------------------
real, dimension (:), intent(in), optional     :: qlevel_in, qmh_in
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!      The general formulas for sigma pressures are: 
!
!      qi(i,k) = q(k)*pss(i)
!    where pss(i)=sfc pressure at point i,
!    and q are specified values (either by formula as a function
!    of KMAX or specified in a DATA statement)
!--------------------------------------------------------------------

      real          :: pss

!-------------------------------------------------------------------
!  allocate local variables
!-------------------------------------------------------------------
      allocate ( q   (KMIN:KMAX)  )
      allocate ( qi  (KMIN:KMAX)  )
      allocate ( qmh (KMIN:KMAX+1))

!-------------------------------------------------------------------
!    define sigmas (q) for the appropriate sigma model structure
!-------------------------------------------------------------------
      if (present(qlevel_in)) then
        do k=1,KMAX
	  q(k) = qlevel_in(k)
        enddo
        do k=1,KMAX+1
	  qmh(k) = qmh_in(k)
        enddo
      else
        type_form = Environment%column_type
        if (trim(type_form) == 'egrpl9') then
	  if (KMAX == 9) then
            do k=1,KMAX
              q(k)=q9e(k)
            enddo
	  else
	    call error_mesg ( 'sigp_std', &
	      'KMAX of model disagrees with designation selected', &
							       FATAL)
	  endif
        else if (trim(type_form) == 'climl9') then
	  if (KMAX == 9) then
            do k=1,KMAX
              q(k)=q9g(k)
            enddo
	  else
	    call error_mesg ( 'sigp_std', &
	      'KMAX of model disagrees with designation selected', &
							       FATAL)
	  endif
        else if (trim(type_form) == 'ccml12') then
	  if (KMAX == 12) then
            do k=1,KMAX
              q(k)=q12(k)
            enddo
	  else
	    call error_mesg ( 'sigp_std', &
	      'KMAX of model disagrees with designation selected', &
							       FATAL)
	  endif
        else if (trim(type_form) == 'climl14') then
	  if (KMAX == 14) then
            do k=1,KMAX
              q(k) = q14(k)
            enddo
	  else
	    call error_mesg ( 'sigp_std', &
	      'KMAX of model disagrees with designation selected', &
							       FATAL)
	  endif
        else if (trim(type_form) == 'egrpl18') then
	  if (KMAX == 18) then
            do k=1,KMAX
              q(k)=q18e(k)
            enddo
	  else
	    call error_mesg ( 'sigp_std', &
	      'KMAX of model disagrees with designation selected', &
							       FATAL)
	  endif
        else if (trim(type_form) == 'climl30') then
	  if (KMAX == 30) then
            do k=1,KMAX
              q(k)=q30(k)
            enddo
	  else
	    call error_mesg ( 'sigp_std', &
	      'KMAX of model disagrees with designation selected', &
							       FATAL)
	  endif
        else if (trim(type_form) == 'nmcl18') then
	  if (KMAX == 18) then
            do k=1,KMAX
              q(k)=q18m(k)
            enddo
	  else
	    call error_mesg ( 'sigp_std', &
	      'KMAX of model disagrees with designation selected', &
							       FATAL)
	  endif
        endif
      endif

!---------------------------------------------------------------------
!    for the case pss=pstd (standard pressure) both formulas
!     reduce to qi(i,k)=q(k) for all k.
!---------------------------------------------------------------------
      pss= pstd

!---------------------------------------------------------------------
!   define layer-mean pressures (pd) from the appropriate q's and pss
!---------------------------------------------------------------------
      do k=1, KMAX
        qi(k)=q(k)
        pd(k)=qi(k)*pss
      enddo
      pd(KMAX+1)=pss
 
!---------------------------------------------------------------------
!   define all layer-boundary pressures (plm) as the mean of adjacent
!   layer-mean pressures. this is necessary for co2 interpolation.
!---------------------------------------------------------------------
      plm(1)=0.
      do k=2,KMAX
        plm(k)=0.5*(pd(k-1)+pd(k))
      enddo
      plm(KMAX+1)=pss

!---------------------------------------------------------------------
!  the following puts p-data into mb
!---------------------------------------------------------------------
      do k=1,KMAX+1
        pd(k)=pd(k)*1.0e-3
        plm(k)=plm(k)*1.0e-3
      enddo

!---------------------------------------------------------------------
!***second pass: pss = ~810mb
!---------------------------------------------------------------------
      pss = pstd_08

!---------------------------------------------------------------------
!   define layer-mean pressures (pd) from the appropriate q's and pss
!---------------------------------------------------------------------
      do k=1, KMAX
        qi(k)=q(k)
        pd8(k)=qi(k)*pss
      enddo
      pd8(KMAX+1)=pss
 
!---------------------------------------------------------------------
!   define all layer-boundary pressures (plm) as the mean of adjacent
!   layer-mean pressures. this is necessary for co2 interpolation.
!---------------------------------------------------------------------
      plm8(1)=0.
      do k=1,KMAX-1
        plm8(k+1)=0.5*(pd8(k)+pd8(k+1))
      enddo
      plm8(KMAX+1)=pss

!---------------------------------------------------------------------
!  the following puts p-data into mb
!---------------------------------------------------------------------
      do k=1,KMAX+1
        pd8(k)=pd8(k)*1.0e-3
        plm8(k)=plm8(k)*1.0e-3
      enddo

!---------------------------------------------------------------------
      deallocate (q)
      deallocate (qi)
      deallocate (qmh)



end subroutine sigp_std






!####################################################################

subroutine lanz_std

     integer      :: hlyr, k, kh

!------------------------------------------------------------------
!   allocate local arrays
!------------------------------------------------------------------
     allocate (zh(KMIN:KMAX+1))
     allocate (zm(KMIN:KMAX))
     allocate (hm(KMIN:KMAX))

!------------------------------------------------------------------
!     pd are data for layer-mean pressures (in mb) for standard lan
!     layers assuming the US Standard atmospheres (1977) profiles at
!     the lan altitudes (in z45) and a surface altitude of 0 m.
!     values are in mb, first value at the top of the atmosphere
!
!------------------------------------------------------------------
     type_form = Environment%column_type
     if (trim(type_form) == 'lanl45') then
       if (KMAX == 45) then	  
         do k=1,KMAX+1
           zh(k) = z45h(k)*1.0E-03
         enddo
         do k=1,KMAX
           zm(k) = z45(k)*1.0E-03
         enddo
       else
         call error_mesg ( 'lanz_std', &
	      'KMAX of model disagrees with designation selected', &
							       FATAL)
       endif
     endif

!------------------------------------------------------------------
!   convert geometric layer-mean altitude (zm) to geopotential
!   altitude (hm) using Eq. (18) in the US standard atmospheres
!   (hm is in km)
!------------------------------------------------------------------
     do k=1,KMAX
       hm(k) = rearth*1.0E-05*zm(k)/(rearth*1.0E-05 + zm(k))
     enddo
  
!------------------------------------------------------------------
!  establish the altitude range to do the integration
!------------------------------------------------------------------
     do k=1,KMAX
       do kh = 1,7
         if (hm(k) .ge. hlyr_top(kh) .and.                   &
      	     hm(k) .lt. hlyr_top(kh+1)) then
	   hlyr = kh
         endif
       enddo

!------------------------------------------------------------------
!  perform the integration, storing results in (pd), with index
!  1 for top layer and (KMAX) for bottom layer
!------------------------------------------------------------------
       if (dtdz(hlyr) .eq. 0.0) then
         pd(KMAX+1-k) = plyr_top(hlyr)*EXP(-grav*wtmair*          &
                     (hm(k)-hlyr_top(hlyr))*1.0e5/(rgas*tlyr_top(hlyr)))
       else
         pd(KMAX+1-k) = plyr_top(hlyr)*                           &
                  EXP(-grav*wtmair/(rgas*dtdz(hlyr)*1.0e-5)*          &
            LOG((tlyr_top(hlyr) + dtdz(hlyr)*(hm(k)-hlyr_top(hlyr))) / &
                tlyr_top(hlyr)))
       endif
     enddo
     pd(KMAX+1) = plyr_top(1)

!------------------------------------------------------------------
!   define all layer-boundary pressures as the mean of adjacent layer-
!   mean pressures. this is necessary for co2 interpolation.
!------------------------------------------------------------------
     plm(1)=0.
     do k=1,KMAX-1
       plm(k+1)=0.5*(pd(k)+pd(k+1))
     enddo
     plm(KMAX+1)=plyr_top(1)

!------------------------------------------------------------------
!  second pass: repeat process, but for geometric altitudes 2 km higher
!------------------------------------------------------------------
     if (trim(type_form) == 'lanl45') then
       do k=1,KMAX+1
         zh(k) = z45h(k)*1.0E-03 + 2.0
       enddo
       do k=1,KMAX
         zm(k) = z45(k)*1.0E-03 + 2.0
       enddo
     endif

!------------------------------------------------------------------
!   convert geometric layer-mean altitude (zm) to geopotential
!   altitude (hm) using Eq. (18) in the US standard atmospheres
!   (hm is in km)
!------------------------------------------------------------------
     do k=1,KMAX
       hm(k) = rearth*1.0E-05*zm(k)/(rearth*1.0E-05 + zm(k))
     enddo
  
!------------------------------------------------------------------
!  establish the altitude range to do the integration
!------------------------------------------------------------------
     do k=1,KMAX
       do kh = 1,7
         if (hm(k) .ge. hlyr_top(kh) .and.              &
             hm(k) .lt. hlyr_top(kh+1)) then
           hlyr = kh
         endif
       enddo

!------------------------------------------------------------------
!  perform the integration, storing results in (pd8), with index
!  1 for top layer and (KMAX) for bottom layer
!------------------------------------------------------------------
       if (dtdz(hlyr) .eq. 0.0) then
         pd8(KMAX+1-k) = plyr_top(hlyr)*EXP(-grav*wtmair*             &
                  (hm(k)-hlyr_top(hlyr))*1.0e5/(rgas*tlyr_top(hlyr)))
       else
         pd8(KMAX+1-k) = plyr_top(hlyr)*                             &
             EXP(-grav*wtmair/(rgas*dtdz(hlyr)*1.0e-5)*             &
            LOG((tlyr_top(hlyr) + dtdz(hlyr)*(hm(k)-hlyr_top(hlyr))) / &
                tlyr_top(hlyr)))
       endif
     enddo
     pd8(KMAX+1) = plyr_top_2km

!------------------------------------------------------------------
!   define all layer-boundary pressures as the mean of adjacent layer-
!   mean pressures. this is necessary for co2 interpolation.
!------------------------------------------------------------------
     plm8(1)=0.
     do k=1,KMAX-1
       plm8(k+1)=0.5*(pd8(k)+pd8(k+1))
     enddo
     plm8(KMAX+1)=plyr_top_2km

!---------------------------------------------------------------------
     deallocate (zh)
     deallocate (zm)
     deallocate (hm)



end subroutine lanz_std




!####################################################################

subroutine set_std_pressures (plm_in, pd_in, plm8_in, pd8_in)

!------------------------------------------------------------------
real,intent(in), dimension(:)  :: pd_in, plm_in, plm8_in, pd8_in
!------------------------------------------------------------------

      plm(:)  = plm_in(:)
      pd(:)   = pd_in(:)
      plm8(:) = plm8_in(:)
      pd8(:)  = pd8_in(:)
!------------------------------------------------------------------


end subroutine set_std_pressures





!####################################################################

subroutine get_std_pressures (plm_out, pd_out, plm8_out, pd8_out)

!------------------------------------------------------------------
real,intent(out), dimension(:), optional  :: pd_out, plm_out,    &
					     plm8_out, pd8_out
!------------------------------------------------------------------


      if (present (  pd_out) ) pd_out(:)   = pd(:)
      if (present ( plm_out) ) plm_out(:)  = plm(:)
      if (present ( pd8_out) ) pd8_out(:)  = pd8(:)
      if (present (plm8_out) ) plm8_out(:) = plm8(:)

!------------------------------------------------------------------


end subroutine get_std_pressures




!####################################################################



                  end module std_pressures_mod
