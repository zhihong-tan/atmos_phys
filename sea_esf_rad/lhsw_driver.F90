		      module lhsw_driver_mod

use rad_step_setup_mod,    only: KS=>KSRAD, KE=>KERAD, ISRAD, IERAD, &
			         JSRAD, JERAD
use rad_utilities_mod,     only: Environment, environment_type
use constants_new_mod,     only: diffac, grav, radcon, alogmin
use std_pressures_mod,     only: get_std_pressures
use co2_mod,               only: rrco2
use cloudrad_package_mod,  only: get_ncldsw, get_clouds_for_lhsw
use  utilities_mod,        only: open_file, file_exist,     &
                                 check_nml_error, error_mesg, &  
                                 print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file

!--------------------------------------------------------------------

implicit none
private


!--------------------------------------------------------------------
!            lacis-hansen shortwave parameterization
!
!------------------------------------------------------------------


!--------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.09'
    character(len=128)  :: version =  '$Id: lhsw_driver.F90,v 1.2 2001/08/30 15:10:24 fms Exp $'
    character(len=128)  :: tag     =  '$Name: galway $'



!---------------------------------------------------------------------
!-------  interfaces --------
 
public  lhsw_driver_init, swrad



!-------------------------------------------------------------------
!-------- namelist ----------


real     :: dummy = 0.0

namelist /lhsw_driver_nml /         &
                            dummy


!--------------------------------------------------------------------
!------ public data -------





!--------------------------------------------------------------------
!------ private data -------


!---------------------------------------------------------------------
!       abcff   = absorption coefficients for bands in k-distribution
!                 that were originally given by lacis and hansen and
!                 revised by ramaswamy.
!       pwts    = the corresponding weights assigned to the 
!                 shortwave radiation k-distribution.
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!   absorption coefficients and weights for bands as revised by 
!   ramaswamy.
!-------------------------------------------------------------------

real, dimension(12)   :: abcff_12, pwts_12

data abcff_12 /4.0000E-05, 4.0000E-05, 0.0020E+00, 0.0350E+00,   &
               0.3770E+00, 1.9500E+00, 9.4000E+00, 4.4600E+01,   &
               1.9000E+02, 9.8900E+02, 2.7060E+03, 3.9011E+04/
  
data pwts_12 /0.50000E+00, 0.121416E+00, 0.06980E+00, 0.15580E+00, &
              0.06310E+00, 0.036200E+00, 0.02430E+00, 0.01580E+00,&
              0.00870E+00, 0.001467E+00, 0.23420E-02, 0.10750E-02/

!---------------------------------------------------------------------
!     the original 9-band lacis-hansen coefficients and weights.
!---------------------------------------------------------------------
 
real, dimension(9)    :: abcff_9, pwts_9

data abcff_9 /4.000E-05, 4.000E-05, 0.002E+00, 0.035E+00,  &
              0.377E-00, 1.950E+00, 9.400E+00, 4.460E+01,  &
              1.900E+02/

data pwts_9 /0.5000E+00, 0.1470E+00, 0.0698E+00, 0.1443E+00,  &
             0.0584E+00, 0.0335E+00, 0.0225E+00, 0.0158E+00, &
             0.0087E+00/

 
!--------------------------------------------------------------------
!   cfco2, cfo3 = conversion factors from gm/cm**2 to cm-atm (standard
!                 temperature and pressure).
!reflo3, rrayav = reflection coefficients given by lacis and hansen to 
!                 account for effects of rayleigh scattering in the 
!                 visible band one frequencies.         
!         prko2 = pressure (mb) below which o2 absorption affects short-
!                 wave transmissivities. 
!       efmago3 = ??????????
!--------------------------------------------------------------------
real           :: cfco2=5.0896E+02
real           :: cfo3=4.6664E+02
real           :: reflo3=1.9000E+00
real           :: rrayav=1.4400E-01
real           :: prko2=2.8
real           :: efmago3=1.900
!---------------------------------------------------------------------
 
integer        :: NB
integer        :: ixprko2


real, dimension(:), allocatable   :: abcff, pwts







!------------------------------------------------------------------
!------------------------------------------------------------------




	contains


subroutine lhsw_driver_init (kmin, kmax)

integer, intent(in)    :: kmin, kmax

      integer :: unit, ierr, io
      integer :: k, KSL, KEL

      real, dimension(:), allocatable :: plm

!---------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=lhsw_driver_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'lhsw_driver_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'lhsw_driver', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=lhsw_driver_nml)
      endif
      call close_file (unit)

!---------------------------------------------------------------------
      allocate ( plm(kmin:kmax+1) )
      call get_std_pressures(plm_out=plm)
      KSL = lbound(plm,1)
      KEL = ubound(plm,1) - 1
 
      if (Environment%using_fms_periphs) then
	NB = 9
	allocate (abcff (9) )
	allocate (pwts  (9) )
	abcff(:) = abcff_9(:)
	pwts(:)  = pwts_9(:)
      else if (Environment%using_sky_periphs) then
        NB = 12
        allocate (abcff (12) )
        allocate (pwts  (12) )
        abcff(:) = abcff_12(:)
        pwts(:)  = pwts_12(:)
      endif

!-------------------------------------------------------------------
!   convert pressure specification for bottom (flux) pressure level
!   for o2 calculation into an index (ixprko2)
!------------------------------------------------------------------- 

      ixprko2 = KSL 
      do k=KSL+1,KEL
        if ((plm(k) - prko2) .LT. 0.0) then
          ixprko2 = k + KSL - 1
        else
          exit
        endif
      enddo

      deallocate (plm)

!------------------------------------------------------------------- 


end subroutine lhsw_driver_init





!######################################################################
 
subroutine swrad (cosangsolar, fracday, with_clouds, cirabgd, cirrfgd, &
	  cvisrfgd, press, qo3, rh2o, ssolar, nsolwg, dfsw,  &
		  fsw, hsw, ufsw, gwt)
       
!-----------------------------------------------------------------------
!
!     Swrad solves for shortwave radiation.
!
!     references:
!
!     (1)  lacis, a. a. and j. e. hansen, "a parameterization for the
!          absorption of solar radiation in the earth's atmosphere," 
!          journal of the atmospheric sciences, 31 (1974), 118-133.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!
!-----------------------------------------------------------------------
!     intent in:
!
!     fracday =  fraction of day (or timestep) that sun is above 
!                horizon.
! 
!     press   =  pressure at data levels of model.
!
!     qo3     =  mass mixing ratio of o3 at model data levels.
!
!     rh2o    =  mass mixing ratio of h2o at model data levels.
!
!     ssolar  =  solar constant (may vary over one year). units: Wm-2.
!
! cosangsolar =  zenith angle at grid point.
!-----------------------------------------------------------------------

logical,                  intent(in)     :: with_clouds
real,                     intent(in)     :: ssolar
integer,                  intent(in)     :: nsolwg
real, dimension(:,:,:),   intent(in)     :: press, qo3, rh2o,  &
					    cosangsolar
real, dimension(:,:),     intent(in)     :: fracday, cirabgd, &
					    cirrfgd, cvisrfgd
real, dimension(:), optional, intent(in) :: gwt

!-----------------------------------------------------------------------
!     intent out:
!
!     dfsw    =  downward radiation at all pressure levels.
!
!     fsw     =  net radiation (up-down) at all pressure levels.
!
!     hsw     =  radiation heating rates at all pressure layers.
!
!     ufsw    =  upward radiation at all pressure levels.
!-----------------------------------------------------------------------

real, dimension(:,:,:), intent(out)   :: dfsw, fsw, hsw, ufsw

!-----------------------------------------------------------------------
!     intent local:
!
!     absdo3  =  absorption of o3 down the atmosphere.
!
!     absuo3  =  absorption of o3 up the atmosphere.
!
!     alfa    =  effective reflective coefficient of a cloud and all
!                clouds below it.
! 
!     alfau   =  effective reflective coefficient excluding current
!                cloud.
!
!     cr      =  coefficient of reflection of a cloud.
!
!     ct      =  coefficient of transmission of a cloud.
!
!     dfn     =  downwards flux for current band n.
!
!     dfnclb  =  down flux at cloud bottom.
!
!     dfnclt  =  down flux at cloud top.
!
!     dfntop  =  down flux at atmosphere top. (in cgs units)
!
!     dfntrn  =  down flux partial product at cloud.
!
!     dp      =  derivative of pressure.
!
!     dpcldi  =  inverse of pressure difference between cloud top and
!                bottom.
!
!     du      =  delta optical path.
!
!     duco2   =  delta optical path for co2.
!
!     duo3    =  delta optical path for o3.
!
!     ff      =  angular optical factor.
!
!     ffco2   =  angular optical for co2.
!
!     ffo3    =  angular optical for o3.
!
!     pp      =  pressure.
!
!     ppbtm   =  pressure at bottom of cloud.
!
!     pptop   =  pressure at top of cloud.
!
!     pr2     =  normalized pressure.
!
!     refl    =  ground coefficent of reflection.
!
!     rray    =  parameterization for rayleigh scrattering.
!
!     secz    =  average secant of solar angle.
!
!     tdclbt  =  transmission down from cloud top to next cloud top
!
!     tdcltt  =  transmission down from cloud top to to next cloud top.
!
!     tdcltop =  transmission down at cloud top.
!
!     tdcltopi=  inverse of tdcltop.
!
!     tdclbtm =  transmission down at cloud bottom.
!
!     tdco2   =  transmission down co2.
!
!     ttd     =  transmission down h2o.
!
!     ttdb1   =  band one transmission quantity.
!
!     ttu     =  transmission up h2o.
!
!     ttub1   =  band one transmission quantity.
!
!     tuclbtm =  transmission up at cloud bottom.
!
!     tucltop =  transmission up at cloud top.
!
!     tucltop =  inverse of tucltop.
!
!     tuco2   =  transmission up co2.
!
!     ud      =  optical path down h2o.
!
!     udco2   =  optical path down co2.
!
!     udo3    =  optical path down o3.
!
!     ufn     =  upward flux for current band n.
!
!     ufnclb  =  upward flux at cloud bottom.
!
!     ufnclt  =  upward flux cloud top.
!
!     ufntrn  =  up flux partial product at cloud.
!
!     uu      =  optical path up h2o.
!
!     uuco2   =  optical path up co2.
!
!     uuo3    =  optical path up o3.
!
!     absdo2  =  absorption down o2.
!     udo2    =  optical path o2 (up-down).
!     uo2     =  optical path o2 (up-down).
!-----------------------------------------------------------------------
      real, dimension(:,:),   allocatable :: refl, rray, seczen
      real, dimension(:,:,:), allocatable ::   &
		  absdo3,        absuo3,          alfa,   &
		  alfau,         alogtt,          cr,       &
		  ct,            dfn,             dfncltop, &
		  dfnlbc,        dfntop,          dp,       &
		  dpcldi,        du,              duco2,    &
		  duo3,          ff,              ffco2,    &
		  ffo3,          pp,              pptop,    &
		  pr2,           tdcltt,          tdclbt,   &
		  tdcltop,       tdclbtm,         tdco2,    &
		  ttd,           ttdb1,           ttu,      &
		  ttub1,         tucltop,         tuco2,    &
		  ud,            udco2,           udo3,     &
		  ufn,           ufncltop,        ufnlbc,   &
		  uu,            uuco2,           uuo3,     &
		  absdo2,        uo2,             dfswg,    &
		  fswg,          hswg,            ufswg
     integer :: k,j,i,kc, kc1, nband, ko, ngp, kcldsw
     real    :: tempu, tempd
     real, dimension(:,:,:,:),  allocatable   :: cirabsw, cirrfsw, &
 				                 cvisrfsw
     real, dimension(:,:,:),    allocatable   :: camtsw
     integer, dimension(:,:),   allocatable   :: ncldsw            
     integer, dimension(:,:,:), allocatable   :: ktopsw, kbtmsw    

!--------------------------------------------------------------------
!     allocate space for and then retrieve the number of clouds in 
!     in each model column.
!------------------------------------------------------------------
      allocate( ncldsw (ISRAD:IERAD, JSRAD:JERAD) )
      call get_ncldsw (ncldsw )

!---------------------------------------------------------------------
!     define the maximum number of clouds in any column of the chunk
!     and the distribution of clouds with grid column to be used in 
!     this execution of the shortwave radiation code.
!---------------------------------------------------------------------
      if (with_clouds) then
        kcldsw = MAXVAL(ncldsw)
      else
        kcldsw = 0
        ncldsw(:,:) = 0
      endif

!--------------------------------------------------------------------
!     allocate space for and then retrieve the shortwave radiative 
!     cloud properties in each model column.
!------------------------------------------------------------------
      if (kcldsw /= 0) then
        allocate( camtsw (ISRAD:IERAD, JSRAD:JERAD, 1:kcldsw) )
        allocate( ktopsw (ISRAD:IERAD, JSRAD:JERAD, 1:kcldsw) )
        allocate( kbtmsw (ISRAD:IERAD, JSRAD:JERAD, 1:kcldsw+1) )
        allocate( cirabsw(ISRAD:IERAD, JSRAD:JERAD, 1:kcldsw, nsolwg) )
        allocate( cirrfsw(ISRAD:IERAD, JSRAD:JERAD, 1:kcldsw, nsolwg) )
        allocate( cvisrfsw(ISRAD:IERAD, JSRAD:JERAD,1:kcldsw, nsolwg) )

        call get_clouds_for_lhsw (cirabsw, cirrfsw, cvisrfsw,  &
	 			  ktopsw, kbtmsw, camtsw)
      endif
  
!--------------------------------------------------------------------
!    allocate local arrays
!--------------------------------------------------------------------

      allocate( absdo3  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( absuo3  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( alfa    (ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( alfau   (ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( alogtt  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( cr      (ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( ct      (ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( dfn     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( dfncltop(ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( dfnlbc  (ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( dfntop  (ISRAD:IERAD, JSRAD:JERAD, NB) )
      allocate( dp      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( dpcldi  (ISRAD:IERAD, JSRAD:JERAD, kcldsw) )
      allocate( du      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( duco2   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( duo3    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( ff      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( ffco2   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( ffo3    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( pp      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( pptop   (ISRAD:IERAD, JSRAD:JERAD, kcldsw) )
      allocate( pr2     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( refl    (ISRAD:IERAD, JSRAD:JERAD) )
      allocate( rray    (ISRAD:IERAD, JSRAD:JERAD) )
      allocate( seczen  (ISRAD:IERAD, JSRAD:JERAD) )
      allocate( tdcltt  (ISRAD:IERAD, JSRAD:JERAD, kcldsw  ))
      allocate( tdclbt  (ISRAD:IERAD, JSRAD:JERAD, kcldsw  ) ) 
      allocate( tdcltop (ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( tdclbtm (ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( tdco2   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( ttd     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( ttdb1   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( ttu     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( ttub1   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( tucltop (ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( tuco2   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( ud      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( udco2   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( udo3    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( ufn     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( ufncltop(ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( ufnlbc  (ISRAD:IERAD, JSRAD:JERAD, kcldsw+1) )
      allocate( uu      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( uuco2   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) ) 
      allocate( uuo3    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )

      allocate( absdo2  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate( uo2     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )

      if (present(gwt)) then
        allocate( dfswg   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
        allocate( fswg    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
        allocate( hswg    (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
        allocate( ufswg   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      endif
!----------------------------------------------------------------------


    do ngp=1,nsolwg

!-----------------------------------------------------------------------
!     calculate secant of zenith angle seczen, flux pressures pp, layer
!     width dp, and pressure scaling factor pr2.  see reference (1).
!-----------------------------------------------------------------------

      seczen(:,:        ) = 3.5E+01/SQRT(1.224E+03*  &
			    cosangsolar(:,:,ngp)*  &
                            cosangsolar(:,:,ngp) + 1.0E+00)
      pp    (:,:,KS     ) = 0.0E+00
      pp    (:,:,KE+1   ) = press(:,:,KE+1)
      pp    (:,:,KS+1:KE) = (press(:,:,KS+1:KE) + press(:,:,KS:KE-1))* &
                             0.5E+00
      dp    (:,:,KS:KE  ) = (pp(:,:,KS+1:KE+1) - pp(:,:,KS:KE  ))
      pr2   (:,:,KS:KE  ) = (pp(:,:,KS:KE  )+ pp(:,:,KS+1:KE+1))*0.5E+00
      do k=KS,KE
        pr2(:,:,k) = pr2(:,:,k)/press(:,:,KE+1)
      end do

!-----------------------------------------------------------------------
!     set up angular factor to be multiplied by the optical factor.
!     above the highest cloud, this is the secant of the zenith angle 
!     seczen (modified slightly for refractive effects).  below the 
!     highest cloud, this is diffac (efmago3 for ozone, in accordance 
!     with lacis-hansen parameterization).  this factor is used 
!     regardless of cloud amount and for direct and diffuse radiation 
!     (this is not a 2-1/2 stream model).
!-----------------------------------------------------------------------

      ff   (:,:,KS:KE+1) = diffac
      ffco2(:,:,KS:KE+1) = diffac
      ffo3 (:,:,KS:KE+1) = efmago3
      do j=JSRAD,JERAD
        do i=ISRAD,IERAD
	  if (ncldsw(i,j) .NE. 0) then
            kc1=ktopsw(i,j,ncldsw(i,j))
	  else
	    kc1 = KE+1
	  endif
          do kc=KS,kc1 
            ffo3 (i,j,kc) = seczen(i,j) 
            ffco2(i,j,kc) = seczen(i,j)
            ff   (i,j,kc) = seczen(i,j) 
          end do
        end do
      end do

!-----------------------------------------------------------------------
!     calculate pressure-weighted optical paths for each layer.
!     pressure weighting uses pr2.  du equals value for h2o; duco2 for
!     co2; duo3 for o3.
!-----------------------------------------------------------------------
      duo3 (:,:,KS:KE) = (cfo3/grav)*qo3(:,:,KS:KE)*dp(:,:,KS:KE)
      duco2(:,:,KS:KE) = (rrco2/grav*cfco2)*dp(:,:,KS:KE)*  &
                          pr2(:,:,KS:KE)
      du   (:,:,KS:KE) = rh2o(:,:,KS:KE)*dp(:,:,KS:KE)* &
                          pr2(:,:,KS:KE)/grav

!-----------------------------------------------------------------------
!     obtain the optical path from the top of the atmosphere to the
!     flux pressure.  angular factors are now included.  ud equals 
!     downward path for h2o, with uu the upward path for h2o. 
!     corresponding quantities for co2 and o3 are udco2/uuco2 and 
!     udo3/uuo3 respectively.
!-----------------------------------------------------------------------
      ud   (:,:,KS) = 0.0E+00 
      udco2(:,:,KS) = 0.0E+00
      udo3 (:,:,KS) = 0.0E+00
      do k=KS+1,KE+1
        ud   (:,:,k) = ud   (:,:,k-1) + du   (:,:,k-1)*ff   (:,:,k) 
        udco2(:,:,k) = udco2(:,:,k-1) + duco2(:,:,k-1)*ffco2(:,:,k)
        udo3 (:,:,k) = udo3 (:,:,k-1) + duo3 (:,:,k-1)*ffo3 (:,:,k) 
      end do
      uu   (:,:,KE+1) = ud   (:,:,KE+1) 
      uuco2(:,:,KE+1) = udco2(:,:,KE+1) 
      uuo3 (:,:,KE+1) = udo3 (:,:,KE+1)
      do k=KE,KS,-1
        uu   (:,:,k) = uu   (:,:,k+1) + du   (:,:,k)*diffac
        uuco2(:,:,k) = uuco2(:,:,k+1) + duco2(:,:,k)*diffac
        uuo3 (:,:,k) = uuo3 (:,:,k+1) + duo3 (:,:,k)*reflo3 
      end do

!-----------------------------------------------------------------------
!     for the skyhi model only, obtain the oxygen optical path, using
!     the o3 angular integration.
!-----------------------------------------------------------------------
      uo2(:,:,KS+1:ixprKO2) = pp(:,:,KS+1:ixprKO2)*    &
                              ffo3(:,:,KS+1:ixprKO2)

!-----------------------------------------------------------------------
!     calculate co2 absorptions.  they will be used in near infrared
!     bands.  since the absorption amount is given (in the formula used
!     below, derived from sasamori) in terms of the total solar flux,
!     and the absorption is only included in the near infrared 
!     (50 percent of the solar spectrum), the absorptions are multiplied
!     by two.
!-----------------------------------------------------------------------
      tdco2(:,:,KS:KE+1) =    &
!                     1.0E+00 - 2.0E+00*(2.35E-03* EXP(0.26E+00*ALOG( &
!                      (udco2(:,:,Ks   :Ke   +1) + &
!                     1.29E-02))) - 7.58265E-04)
                     1.0E+00 - 2.0E+00*(2.35E-03*(udco2(:,:,KS:KE+1) + &
                     1.29E-02)**2.6E-01 - 7.58265E-04)
      tuco2(:,:,KS:KE+1) =    &
!          1.0E+00 - 2.0E+00*(2.35E-03* EXP(0.26E+00*ALOG(  &
!                     (uuco2(:,:,Ks   :Ke   +1) +  &
!                     1.29E-02))) - 7.58265E-04)
                     1.0E+00 - 2.0E+00*(2.35E-03*(uuco2(:,:,KS:KE+1) + &
                     1.29E-02)**2.6E-01 - 7.58265E-04)

!-----------------------------------------------------------------------
!     now calculate ozone absorptions.  these will be used in the
!     visible band one.  just as in the co2 case, since this band is
!     50 percent of the solar spectrum, the absorptions are multiplied
!     by two.  see reference (1).
!-----------------------------------------------------------------------
       absdo3(:,:,KS:KE+1) =   &
!     &   2.0E+00*udo3(:,:,Ks :Ke   +1)*(1.082E+00*EXP(-.805E+00*ALOG( &
!     &    (1.0E+00 + 1.386E+02*udo3(:,:,Ks   :Ke   +1)))) + &
!     &    6.58E-02/(1.0E+00 + (1.036E+02)**3* &
!     &    udo3(:,:,Ks   :Ke   +1)**3) +  &
!     &    2.118E-02/(1.0E+00 + udo3(:,:,Ks   :Ke   +1)*(4.2E-02 + &
!     &    3.23E-04*   &
!     &    udo3(:,:,KS   :Ke   +1))))
          2.0E+00*udo3(:,:,KS:KE+1)*(1.082E+00*    &
          (1.0E+00 + 1.386E+02*udo3(:,:,KS:KE+1))**(-8.05E-01) +   &
          6.58E-02/(1.0E+00 + (1.036E+02)**3*udo3(:,:,KS:KE+1)**3) + &
          2.118E-02/(1.0E+00 + udo3(:,:,KS:KE+1)*(4.2E-02 + 3.23E-04* &
          udo3(:,:,KS:KE+1))))
       absuo3(:,:,KS:KE+1) =   &

!         2.0E+00*uuo3(:,:,Ks :Ke   +1)*(1.082E+00*EXP(-.805E+00*ALOG( &
!         (1.0E+00 + 1.386E+02*uuo3(:,:,Ks   :Ke   +1)))) +  &
!         6.58E-02/(1.0E+00 + (1.036E+02)**3*  &
!         uuo3(:,:,Ks   :Ke   +1)**3) +  &
!         2.118E-02/(1.0E+00 + uuo3(:,:,Ks   :Ke   +1)*(4.2E-02 + &
!         3.23E-04*   &
!         uuo3(:,:,Ks   :Ke   +1)))) 
          2.0E+00*uuo3(:,:,KS:KE+1)*(1.082E+00*  &
          (1.0E+00 + 1.386E+02*uuo3(:,:,KS:KE+1))**(-8.05E-01) + &
          6.58E-02/(1.0E+00 + (1.036E+02)**3*uuo3(:,:,KS:KE+1)**3) + &
          2.118E-02/(1.0E+00 + uuo3(:,:,KS:KE+1)*(4.2E-02 + 3.23E-04* &
          uuo3(:,:,KS:KE+1))))

!---------------------------------------------------------------------
!     calculate o2 absorptions (in visible band one).
!
!     formula is: abs=1.02E-05*uo2(k)**0.3795 for uo2 < 673.9057
!
!                 abs=4.51E-06*uo2(k)**0.5048 for uo2 > 673.9057
!
!     the absorption is constant below 35 km (level ixprKO2).
!---------------------------------------------------------------------
    if (Environment%using_sky_periphs) then
      absdo2(:,:,KS) = 0.0E+00
      do ko=KS+1,ixprKO2
        do j=JSRAD,JERAD
          do i=ISRAD,IERAD
            if(uo2(i,j,ko) .LE. 6.739057E+02) then
!	    absdo2(i,j,ko) = 1.02E-05*EXP(.3795E+00*ALOG(uo2(i,j,ko)))
               absdo2(i,j,ko) = 1.02E-05*uo2(i,j,ko)**0.3795E+00
            else
!	    absdo2(i,j,ko) = 4.51E-06*EXP(.5048E+00*ALOG(uo2(i,j,ko)))
               absdo2(i,j,ko) = 4.51E-06*uo2(i,j,ko)**0.5048E+00
            endif
          end do
        end do
      end do

!-----------------------------------------------------------------------
!     add o2 absorption to o3 absorption.
!-----------------------------------------------------------------------
      do k=KS,KE+1
        if(k .LE. ixprKO2) then
          absdo3(:,:,k) = absdo3(:,:,k) + absdo2(:,:,k  )
        else
          absdo3(:,:,k) = absdo3(:,:,k) + absdo2(:,:,ixprKO2)
        endif
        absuo3(:,:,k) = absuo3(:,:,k) + absdo2(:,:,ixprKO2)
      end do
    endif

!-----------------------------------------------------------------------
!     computations for bands have been divided into: band one, band two,
!     and band three through NB.  begin by calculating flux entering
!     at the top for each band.
!-----------------------------------------------------------------------
      do nband=1,NB
        dfntop(:,:,nband) = ssolar*1.0E+03*cosangsolar(:,:,ngp)*   &
                                fracday(:,:)*pwts(nband)
      end do

!-----------------------------------------------------------------------
!     initialize dfsw and ufsw for the bands.
!-----------------------------------------------------------------------
      if (present(gwt)) then
        dfswg(:,:,KS:KE+1) = 0.0E+00
        ufswg(:,:,KS:KE+1) = 0.0E+00
      endif

!-----------------------------------------------------------------------
!     execute the reflectivity parameterization for the visible band 
!     one.  see reference (1)
!-----------------------------------------------------------------------
      rray(:,:) = 2.19E-01/(1.0E+00 + 8.16E-01*cosangsolar(:,:,ngp))
      refl(:,:) = rray(:,:) + (1.0E+00 - rray(:,:))*(1.0E+00 -    &
                  rrayav)*cvisrfgd(:,:)/(1.0E+00 - cvisrfgd(:,:  )* &
                  rrayav)

!-----------------------------------------------------------------------
!     where clouds exist.
!-----------------------------------------------------------------------
      if(kcldsw .NE. 0) then
!-----------------------------------------------------------------------
!     the first cloud is the ground; its properties are given by refl. 
!     the transmission ct(:,:,1) is irrelevant for now.
!-----------------------------------------------------------------------
        cr(:,:,1) = refl(:,:)
!-----------------------------------------------------------------------
!     obtain cloud reflection and transmission coefficients for
!     remaining clouds in the visible band one.  the maximum number of 
!     clouds kcldsw is used.  this creates extra work and may be removed
!     in a subsequent update.
!-----------------------------------------------------------------------
        cr(:,:,2:kcldsw+1) = cvisrfsw(:,:,1:kcldsw,ngp)*   &
 		     camtsw(:,:,1:kcldsw)
        ct(:,:,2:kcldsw+1) = 1.0E+00 - cr(:,:,2:kcldsw+1)
      end if

!-----------------------------------------------------------------------
!     begin computations for visible band one, near infrared band two,
!     and near infrared bands three thru NB.
!-----------------------------------------------------------------------
      do nband=1,NB
!-----------------------------------------------------------------------
!     calculations for visible band one which includes o3, o2, and 
!     (negligible) h2o absorption.
!-----------------------------------------------------------------------
        if(nband .EQ. 1) then
          alogtt(:,:,KS:KE+1) = -abcff(1)*ud(:,:,KS:KE+1)
          ttdb1 (:,:,KS:KE+1) = EXP(MAX(alogmin, alogtt(:,:,KS:KE+1)))
          ttd   (:,:,KS:KE+1) = ttdb1 (:,:,KS:KE+1)*(1.0E+00 -   &
                                absdo3(:,:,KS:KE+1))
          alogtt(:,:,KS:KE+1) = -abcff(1)*uu(:,:,KS:KE+1)
          ttub1 (:,:,KS:KE+1) = EXP(MAX(alogmin, alogtt(:,:,KS:KE+1)))
          ttu   (:,:,KS:KE+1) = ttub1 (:,:,KS:KE+1)*(1.0E+00 -   &
                                absuo3(:,:,KS:KE+1))

!-----------------------------------------------------------------------
!     calculations for the near infrared band two where the water vapor
!     transmission function ttdb1 and ttub1 for band two is equal to
!     that of band one.
!-----------------------------------------------------------------------
	else if(nband .EQ. 2) then
          alogtt(:,:,KS:KE+1) = -abcff(1)*ud(:,:,KS:KE+1)
          ttdb1 (:,:,KS:KE+1) = EXP(MAX(alogmin,alogtt(:,:,KS:KE+1)))
          ttd   (:,:,KS:KE+1) = ttdb1(:,:,KS:KE+1)*tdco2(:,:,KS:KE+1)
          alogtt(:,:,KS:KE+1) = -abcff(1)*uu(:,:,KS:KE+1)
          ttub1 (:,:,KS:KE+1) = EXP(MAX(alogmin,alogtt(:,:,KS:KE+1)))
          ttu   (:,:,KS:KE+1) = ttub1(:,:,KS:KE+1)*tuco2(:,:,KS:KE+1)

!-----------------------------------------------------------------------
!     calculate water vapor transmission functions for near infrared
!     bands.  include co2 transmission tdco2/tuco2, which is the
!     same for all infrared bands.
!-----------------------------------------------------------------------
	else if(nband .GE. 3) then
          alogtt(:,:,KS:KE+1) = -abcff(nband)*ud(:,:,KS:KE+1)
          ttd   (:,:,KS:KE+1) = EXP(MAX(alogmin,alogtt(:,:,KS:KE+1)))* &
                                tdco2(:,:,KS:KE+1)
          alogtt(:,:,KS:KE+1) = -abcff(nband)*uu(:,:,KS:KE+1)
          ttu   (:,:,KS:KE+1) = EXP(MAX(alogmin,alogtt(:,:,KS:KE+1)))* &
                                tuco2(:,:,KS:KE+1)
        endif

!-----------------------------------------------------------------------
!     at this point, include ttd(KS), ttu(KE+1), noting that ttd(KS)=1 
!     for all bands, and that ttu(KE+1)=ttd(KE+1) for all bands.
!-----------------------------------------------------------------------
        ttd(:,:,KS  ) = 1.0E+00
        ttu(:,:,KE+1) = ttd(:,:,KE+1)
!-----------------------------------------------------------------------
!     where no clouds exist.
!-----------------------------------------------------------------------
        if (kcldsw .EQ. 0) then
          if(nband .EQ. 1) then
            do k=KS,KE+1
              ufn(:,:,k) = refl(:,:)*ttu(:,:,k)
              dfn(:,:,k) = ttd(:,:,k)
            end do
          else
            do k=KS,KE+1
              ufn(:,:,k) = cirrfgd(:,:)*ttu(:,:,k)
              dfn(:,:,k) = ttd(:,:,k)
            end do      
          endif
!-----------------------------------------------------------------------
!     where clouds exist.
!-----------------------------------------------------------------------
        else
!-----------------------------------------------------------------------
!     for execution of the cloud loop, it is necessary to separate the
!     transmission functions at the top and bottom of the clouds, for
!     each band.  the required quantities are:
!
!       ttd(:,:,ktopsw(:,:,kc),nband) kc=1,ncldsw(:,:)+1 
!       ttd(:,:,kbtmsw(:,:,kc),nband) kc=1,ncldsw(:,:)+1 
!       ttu(:,:,ktopsw(:,:,kc),nband) kc=1,ncldsw(:,:)+1.
!
!     the above quantities are stored in tdcltop, tdclbtm, and tucltop 
!     respectively, as they have multiple use in the program.
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!     for first cloud layer (i.e. the ground) tdcltop and tucltop are
!     known.
!-----------------------------------------------------------------------
          tdcltop (:,:,1) = ttd    (:,:,KE+1)
          tucltop (:,:,1) = ttu    (:,:,KE+1)
          tdclbtm (:,:,1) = tdcltop(:,:,1 )

!-----------------------------------------------------------------------
!     use gathers for remaining ones.
!-----------------------------------------------------------------------
              do kc=1,kcldsw
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
                tdcltop(i,j,kc+1) = ttd(i,j,ktopsw(i,j,kc))
                tucltop(i,j,kc+1) = ttu(i,j,ktopsw(i,j,kc))
                tdclbtm(i,j,kc+1) = ttd(i,j,kbtmsw(i,j,kc))
              end do
            end do
          end do
!-----------------------------------------------------------------------
!     compute the transmissivity from the top of cloud kc+1 to the top
!     of cloud kc.  the cloud transmission ct is included.  this
!     quantity is called tdclbt(kc)-transmission down cloud bottom-top.
!     also, obtain the transmissivity from the bottom of cloud kc+1
!     to the top of cloud kc (a path entirely outside clouds).  this
!     quantity is called tdcltt(kc)-transmission down cloud top-top.
!-----------------------------------------------------------------------
          tdclbt(:,:,1:kcldsw) = tdcltop(:,:,1:kcldsw)/  &
                              tdcltop(:,:,2:kcldsw+1)*ct(:,:,2:kcldsw+1)
          tdcltt(:,:,1:kcldsw) = tdcltop(:,:,1:kcldsw)/   &
                              tdclbtm(:,:,2:kcldsw+1)

!-----------------------------------------------------------------------
!     the following is the recursion relation for alfa. the reflection
!     coefficient for a system including the cloud in question and the
!     flux coming out of the cloud system including all clouds below
!     the cloud in question.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     alfau is alfa without the reflection of the cloud in question.
!----------------------------------------------------------------------!
          alfa (:,:,1) = cr(:,:,1)
          alfau(:,:,1) = 0.0E+00
          do kc=1,kcldsw
!-----------------------------------------------------------------------
!     excessive calculations.
!-----------------------------------------------------------------------
            alfau(:,:,kc+1) = (tdclbt(:,:,kc)**2)*alfa(:,:,kc)/ &
                            (1.0E+00 - (tdcltt(:,:,kc)**2)* &
                            alfa(:,:,kc)*cr(:,:,kc+1))
            alfa (:,:,kc+1) = alfau(:,:,kc+1) + cr(:,:,kc+1)
          end do
!-----------------------------------------------------------------------
!     the following calculation is done for visible band one case only.
!     obtain the pressure at the top, bottom and the thickness of thick 
!     clouds (those at least 2 layers thick).  this is used later is 
!     obtaining fluxes inside the thick clouds, for all frequency bands.
!-----------------------------------------------------------------------
          if(nband .EQ. 1) then
            do kc=1,kcldsw
              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
                  if(ktopsw(i,j,kc  ) .NE. KS) then 
                  if((kbtmsw(i,j,kc  ) - ktopsw(i,j,kc  )) .GT. 1) then
                      pptop (i,j,kc) = pp(i,j,ktopsw(i,j,kc  )) 
                      dpcldi(i,j,kc) = 1.0E+00/(pptop(i,j,kc) -  &
                                       pp(i,j,kbtmsw(i,j,kc  )))
                    endif
                  endif
                end do
              end do
            end do
            cr(:,:,1) = cirrfgd(:,:)
            ct(:,:,1) = 1.0E+00 - (cirrfgd(:,:) + cirabgd(:,:))
            cr(:,:,2:kcldsw+1) = cirrfsw(:,:,1:kcldsw,ngp)*  &
                                 camtsw(:,:,1:kcldsw)
            ct(:,:,2:kcldsw+1) = 1.0E+00 - camtsw(:,:,1:kcldsw  )*  &
                                 (cirrfsw(:,:,1:kcldsw,ngp  ) +    &
                                  cirabsw(:,:,1:kcldsw,ngp  ))
          endif

!-----------------------------------------------------------------------
!     calculate ufn at cloud tops and dfn at cloud bottoms note that 
!     ufncltop(:,:,kcldsw+1) gives the upward flux at the top of the 
!     highest real cloud (if ncldsw(:,:)=kcldsw).  it gives the flux at 
!     the top of the atmosphere if ncldsw(:,:) < kcldsw.  it the first 
!     case, tdcltop equals the transmission function to the top of the 
!     highest cloud, as we want.  in the second case, tdcltop=1, so
!     ufncltop equals alfa. this is also correct.
!-----------------------------------------------------------------------
          ufncltop(:,:,kcldsw+1) = alfa(:,:,kcldsw+1)*   &
                                   tdcltop(:,:,kcldsw+1)
          dfncltop(:,:,kcldsw+1) = tdcltop(:,:,kcldsw+1)
!-----------------------------------------------------------------------
!     this calculation is the reverse of the recursion relation used
!     above.
!-----------------------------------------------------------------------
          do kc=kcldsw,1,-1
            ufncltop(:,:,kc) = ((ufncltop(:,:,kc+1)*alfau(:,:,kc+1))/ &
                               alfa(:,:,kc+1))/tdclbt(:,:,kc)
            dfncltop(:,:,kc) = ufncltop(:,:,kc)/alfa (:,:,kc)
          end do

!-----------------------------------------------------------------------
!     now obtain dfn and ufn for levels between the clouds.
!-----------------------------------------------------------------------
          ufnlbc(:,:,1:kcldsw+1) = ufncltop(:,:,1:kcldsw+1)/  &
                                   tucltop (:,:,1:kcldsw+1)
          dfnlbc(:,:,1:kcldsw+1) = dfncltop(:,:,1:kcldsw+1)/  &
                                   tdcltop (:,:,1:kcldsw+1)
!-----------------------------------------------------------------------
!     case of kc=1 (from the ground to the bottom of the lowest cloud).
!-----------------------------------------------------------------------
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
              do k=kbtmsw(i,j,1),KE+1
                ufn(i,j,k) = ufnlbc(i,j,1)*ttu(i,j,k)
                dfn(i,j,k) = dfnlbc(i,j,1)*ttd(i,j,k)
              end do
            end do
          end do
!-----------------------------------------------------------------------
!     remaining levels.
!-----------------------------------------------------------------------
            do j=JSRAD,JERAD
              do i=ISRAD,IERAD
		do kc=1,ncldsw(i,j)
                  do k=kbtmsw(i,j,kc+1),ktopsw(i,j,kc  )
                    ufn(i,j,k) = ufnlbc(i,j,kc+1)*ttu(i,j,k) 
                    dfn(i,j,k) = dfnlbc(i,j,kc+1)*ttd(i,j,k)
                  end do
!-----------------------------------------------------------------------
!     for the thick clouds, the flux divergence through the cloud layer
!     is assumed to be constant.  the flux derivative is given by
!     tempu (for the upward flux) and tempd (for the downward flux).
!-----------------------------------------------------------------------
                  if((kbtmsw(i,j,kc  ) - ktopsw(i,j,kc  )) .GT. 1) then
               tempu=(ufncltop(i,j,kc+1) - ufn(i,j,kbtmsw(i,j,kc  )))*&
                           dpcldi(i,j,kc  )
                tempd=(dfncltop(i,j,kc+1) - dfn(i,j,kbtmsw(i,j,kc  )))*&
                           dpcldi(i,j,kc  )
                    do k=ktopsw(i,j,kc  )+1,kbtmsw(i,j,kc  )-1
                      ufn(i,j,k) = ufncltop(i,j,kc+1) +   &
                                   tempu*(pp(i,j,k) - pptop(i,j,kc  ))
                      dfn(i,j,k) = dfncltop(i,j,kc+1) +   &
                                   tempd*(pp(i,j,k) - pptop(i,j,kc  ))
                    end do
                  endif
              end do
            end do
          end do
        endif  ! (end of kcldsw loop)

!-----------------------------------------------------------------------
!     scale the previously computed fluxes by the flux at the top of the
!     atmosphere.
!-----------------------------------------------------------------------
        do k=KS,KE+1
          dfn(:,:,k) = dfn(:,:,k)*dfntop(:,:,nband)
          ufn(:,:,k) = ufn(:,:,k)*dfntop(:,:,nband)
        end do

!-----------------------------------------------------------------------
!     sum over bands.
!-----------------------------------------------------------------------
        if (present(gwt)) then
          dfswg(:,:,KS:KE+1) = dfswg(:,:,KS:KE+1) + dfn(:,:,KS:KE+1)
          ufswg(:,:,KS:KE+1) = ufswg(:,:,KS:KE+1) + ufn(:,:,KS:KE+1)
        else
          dfsw (:,:,KS:KE+1) = dfsw (:,:,KS:KE+1) + dfn(:,:,KS:KE+1)
          ufsw (:,:,KS:KE+1) = ufsw (:,:,KS:KE+1) + ufn(:,:,KS:KE+1)
        endif
      end do

!-----------------------------------------------------------------------
!     obtain the net flux and the shortwave heating rate for all bands.
!-----------------------------------------------------------------------
      if (present(gwt)) then
        fswg(:,:,KS:KE+1) = ufswg(:,:,KS:KE+1) - dfswg(:,:,KS:KE+1)
        hswg(:,:,KS:KE  ) = radcon*(fswg(:,:,KS+1:KE+1) -    &
	  		    fswg(:,:,KS:KE))/dp(:,:,KS:KE)
        fsw  (:,:,:) = fsw  (:,:,:) + gwt(ngp)*fswg (:,:,:)
        dfsw (:,:,:) = dfsw (:,:,:) + gwt(ngp)*dfswg(:,:,:)
        ufsw (:,:,:) = ufsw (:,:,:) + gwt(ngp)*ufswg(:,:,:)
        hsw  (:,:,:) = hsw  (:,:,:) + gwt(ngp)*hswg (:,:,:)
      else
        fsw (:,:,KS:KE+1) = ufsw (:,:,KS:KE+1) - dfsw (:,:,KS:KE+1)
        hsw (:,:,KS:KE  ) = radcon*(fsw (:,:,KS+1:KE+1) -    &
	  		    fsw (:,:,KS:KE))/dp(:,:,KS:KE)
      endif
    end do ! (ngp loop)

!--------------------------------------------------------------------
!   deallocate local arrays
!-------------------------------------------------------------------
      if (present(gwt)) then
        deallocate (ufswg      )
        deallocate (hswg       )
        deallocate (fswg       )
        deallocate (dfswg      )
      endif
      deallocate ( absdo3    )
      deallocate ( absuo3    )
      deallocate ( alfa      )
      deallocate ( alfau     )
      deallocate ( alogtt    )
      deallocate ( cr        )
      deallocate ( ct        )
      deallocate ( dfn       )
      deallocate ( dfncltop  )
      deallocate ( dfnlbc    )
      deallocate ( dfntop    )
      deallocate ( dp        )
      deallocate ( dpcldi    )
      deallocate ( du        )
      deallocate ( duco2     )
      deallocate ( duo3      )
      deallocate ( ff        )
      deallocate ( ffco2     )
      deallocate ( ffo3      )
      deallocate ( pp        )
      deallocate ( pptop     )
      deallocate ( pr2       )
      deallocate ( refl      )
      deallocate ( rray      )
      deallocate ( seczen    )
      deallocate ( tdcltt    )
      deallocate ( tdclbt    )
      deallocate ( tdcltop   )
      deallocate ( tdclbtm   )
      deallocate ( tdco2     )
      deallocate ( ttd       )
      deallocate ( ttdb1     )
      deallocate ( ttu       )
      deallocate ( ttub1     )
      deallocate ( tucltop   )
      deallocate ( tuco2     )
      deallocate ( ud        )
      deallocate ( udco2     )
      deallocate ( udo3      )
      deallocate ( ufn       )
      deallocate ( ufncltop  )
      deallocate ( ufnlbc    )
      deallocate ( uu        )
      deallocate ( uuco2     )
      deallocate ( uuo3      )
      deallocate ( absdo2    )
      deallocate ( uo2       )
  
      if (kcldsw /= 0) then
        deallocate ( camtsw )
        deallocate ( cirabsw    )
        deallocate ( cirrfsw    )
        deallocate ( cvisrfsw     )
        deallocate ( ktopsw )
        deallocate ( kbtmsw )
      endif

      deallocate ( ncldsw )

!-------------------------------------------------------------------

      end subroutine swrad 



!##################################################################


             end module lhsw_driver_mod
