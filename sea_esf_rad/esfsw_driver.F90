
		       module esfsw_driver_mod


use rad_step_setup_mod,   only:  ISRAD, IERAD, JSRAD, JERAD, &
				 KSRAD, KERAD, &
				 pflux, deltaz
use  utilities_mod,       only:  open_file, file_exist,    &
                                 check_nml_error, error_mesg, &
                                 print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file
use radiative_gases_mod,  only : rrvco2
use constants_new_mod,    only:  radcon_mks, o2mixrat, gasconst, &
				 rhoair, pie, grav_mks, pstd_mks
use esfsw_parameters_mod, only:  nbands, nfrqpts, put_solarfluxes, &
				 TOT_WVNUMS, nh2obands, nstreams
use esfsw_scattering_mod, only:  get_swaerprops_from_scattering, &
                                 scatpar
use cloudrad_package_mod, only:  get_swcldprops_from_cloudrad
use rad_utilities_mod,    only:  Rad_control, radiation_control_type, &
                                 Sw_control, shortwave_control_type

!---------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!                    esf shortwave driver module
!
!------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.08'
   character(len=128)  :: version =  '$Id: esfsw_driver.F90,v 1.3 2001/10/25 17:48:29 fms Exp $'
   character(len=128)  :: tag     =  '$Name: fez $'



!---------------------------------------------------------------------
!-------  interfaces --------

public   esfsw_driver_init, swresf


private  adding, deledd


!---------------------------------------------------------------------
!-------- namelist  ---------

real              :: dummy = 1.0
  


namelist / esfsw_driver_nml /    &
                                dummy




!---------------------------------------------------------------------
!------- public data ------



!---------------------------------------------------------------------
!------- private data ------


real, dimension (:), allocatable    :: betaddensitymol
real, dimension (:), allocatable :: c1co2, c1co2str, c1o2, c1o2str, &
	                       c2co2, c2co2str, c2o2, c2o2str, &
 			       c3co2, c3co2str, c3o2, c3o2str, &
 			       c4co2, c4co2str, c4o2, c4o2str, &
 			       powph2o, p0h2o
real                        :: c1o2strschrun, c2o2strschrun, &
		               c3o2strschrun, c4o2strschrun, &
			       refquanray, solflxtotal
integer                     :: FIRSTRAYBAND, NIRBANDS
integer, dimension (:), allocatable :: nfreqpts
real,    dimension(:), allocatable  :: solflxband
real, dimension (:), allocatable   :: kh2o, ko3, wtfreq
logical, dimension(:), allocatable :: strterm

real, dimension(:), allocatable   :: wtstr, cosangstr

real, dimension(4)          :: wtstr_4 =(/0.347854845, 0.652145155,&
		                        0.347854845, 0.652145155/)

integer                     :: IJK
logical                     :: do_esfsw_band_diagnostics = .false.
logical                     :: nstr4 = .false.

!------------------------------------------------------------------
!  variables used in subroutine deledd  
!------------------------------------------------------------------
integer, dimension(:,:,:), allocatable ::  cld_index
!real,    dimension(:,:),   allocatable ::  q, r, t
real,    dimension(:,:),   allocatable ::  q
real,    dimension(:  ),   allocatable ::  r, t
real,    dimension(:),     allocatable ::  gstr2, taustr2, omegastr2, &
			  	           cosangzk2, rlayerdir2,    &
				           tlayerde2, tlayerdir2, &
 		         	           sumr, sumt

!------------------------------------------------------------------
!  variables used in subroutine adding  
!------------------------------------------------------------------
real, dimension(:,:,:), allocatable ::  alpp, dm1tl, dm2tl, dm3, dm3r, &
					dm3r1p, radddowndif,    &
					raddupdif, raddupdir, rdm2tl, &
                                        tadddowndir, tlevel

!---------------------------------------------------------------------
!---------------------------------------------------------------------
 




                          contains
 

subroutine esfsw_driver_init
 
!---------------------------------------------------------------------- 
! define the time-independent quantities associated with the incoming 
! shortwave radiation in the multiple-band solar radiation parameter-
! ization.                                    
!---------------------------------------------------------------------- 
!                                                                       
! solflxband   = the solar flux in each parameterization band           
!                                                                       
! endwvnbands  = the wavenumber value for the band limits               
!                                                                       
!                                                                       
! nintsolar    = the number of wavenumber regions where the solar flux  
!                is constant.                                           
!                                                                       
! nwvnsolar    = the number of wavenumbers in each region where the     
!                solar flux is constant                                 
!                                                                       
! solint       = the solar flux in watts per meter**2 in each           
!                wavenumber region where it is constant                 
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
! define parameters for use in determining the absorption due to co2,   
! h2o, o2 and o3 in each band.                                          
!                                                                       
! nfreqpts     = the number of pseudo-monochromatic frequencies         
!                                                                       
! firstrayband = the first band number where the contribution by        
!                rayleigh scattering is included in the solar           
!                calculations                                           
!                                                                       
! nirbands     = the number of bands in the near-infrared (used in      
!                assigning the value of the surface albedo for the      
!                near-infrared, and the visible and ultraviolet         
!                regions, separately)                                   
!                                                                       
! powph2o      = the scaling factor used in the fit of the h2o          
!                transmission function                                  
!                                                                       
! p0h2o        = the reference pressure (mb) used in the fit of the     
!                h2o transmission function                              
!                                                                       
! c(n)co2(str) = coefficients for the absorptivity expression for co2   
!                for the pressure-scaled and non-scaled, respectively,  
!                portions of the fit (=1.0e-99 if no absorption)        
!                                                                       
! c(n)o2(str)  = coefficients for the absorptivity expression for o2    
!                for the pressure-scaled and non-scaled, respectively,  
!                portions of the fit (=1.0e-99 if no absorption)        
!                                                                       
! ""(schrun)   = coefficients for the absorptivity expression for the   
!                Schuman-Runge o2 band (non-scaled only)                
!                                                                       
! wtfreq       = the weight associated with each exponential term       
!                                                                       
! kh2o         =  the psuedo-absorption coefficients in cm2/gm for h2o  
!                                                                       
! ko3          = the absorption coefficients in cm2/gm for o3           
!                                                                       
! strterm      = logical flag to indicate whether or not a h2o pseudo-  
!                absorption coefficient is assigned a non-scaled        
!                (true) or pressure-scaled (false) gas amount           
! ptstr        = gaussian points and weights for evaluation of the
!                diffuse beam.
!----------------------------------------------------------------------c

      real, dimension(4)          ::  ptstr_4 = (/-0.861136312,&
                                               -0.339981044, &
					        0.861136312,  &
						0.339981044 /)
      real                        ::  ptstr_1 = 0.2
      real, dimension(:), allocatable     ::  freqnu, ptstr

!----------------------------------------------------------------------c
! parameters for determining rayleigh optical depth.                   c
!----------------------------------------------------------------------c
 
      real      :: temprefray  = 288.15
      real      :: pressrefray = 101325.       ! MKS units
      real      :: densmolref  = 2.54743E+19
      real      :: convfac     = 1.0E+18

!---------------------------------------------------------------------
      integer   :: iounit, nband, nf, ni, nw, nw1, nw2, nintsolar
      integer   :: unit, io, ierr
      integer   :: i
      real      :: corrfac, gamma, f1, f2, f3, pinteg, &
	           twopiesq, densmolrefsqt3, wavelength, freqsq, &
		   ristdm1, ri

character(len=64)                   :: file_name
integer, dimension(:), allocatable  :: nwvnsolar
real   , dimension(:), allocatable  :: solint   
real   , dimension(TOT_WVNUMS)      :: solarfluxtoa        

integer, dimension(:), allocatable        :: endwvnbands

!--------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=esfsw_driver_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'esfsw_driver_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'esfsw_driver', version_number)
      if (get_my_pe() == 0) then 
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=esfsw_driver_nml)
      endif
      call close_file (unit)

!---------------------------------------------------------------------
!     allocate module variables
!---------------------------------------------------------------------
      allocate ( betaddensitymol (nbands) )
      allocate ( c1co2   (nh2obands),  &
                 c1co2str(nh2obands),  &
                 c1o2    (nh2obands),  &
                 c1o2str (nh2obands),  &
                 c2co2   (nh2obands),  &
                 c2co2str(nh2obands),  &
                 c2o2    (nh2obands),  &
                 c2o2str (nh2obands),  &
                 c3co2   (nh2obands),  &
                 c3co2str(nh2obands),  &
                 c3o2    (nh2obands),  &
                 c3o2str (nh2obands),  &
                 c4co2   (nh2obands),  &
                 c4co2str(nh2obands),  &
                 c4o2    (nh2obands),  &
                 c4o2str (nh2obands),  &
                 powph2o (nh2obands),  &
                 p0h2o   (nh2obands)    )
      allocate ( nfreqpts        (nbands) )
      allocate ( solflxband      (nbands) )
      allocate ( kh2o            (nfrqpts),   & 
		 ko3             (nfrqpts),   &
		 wtfreq          (nfrqpts),   &
		 strterm         (nfrqpts)   )
      allocate ( wtstr           (nstreams),   & 
		 cosangstr       (nstreams)  )


!---------------------------------------------------------------------
!   allocate local variables.
!---------------------------------------------------------------------
  
      allocate (endwvnbands (0:nbands) )
      allocate (freqnu      (nbands) )
      allocate (ptstr       (nstreams) )

      if (nstreams == 4) then
        ptstr(:) = ptstr_4(:)
        wtstr(:) = wtstr_4(:)
        nstr4 = .true.
      else if (nstreams == 1) then
        ptstr(1) = ptstr_1
        wtstr(1) = 1.0
        nstr4 = .false.
      endif

!---------------------------------------------------------------------
!  read input file for band positions, solar fluxes and band
!  strengths
!---------------------------------------------------------------------
      if (nbands == 25 .and. nfrqpts == 72) then 
	file_name = 'INPUT/esf_sw_input_data_n72b25'
      else if (nbands == 18 .and. nfrqpts == 38) then
	file_name = 'INPUT/esf_sw_input_data_n38b18'
      else
	call error_mesg ( 'esfsw_driver_init', &
	 'input file for desired bands and frqs is not available', &
	 FATAL)
      endif
      iounit = open_file (file_name, action='read',  &
			   form='formatted' )
      read(iounit,101) ( solflxband(nband), nband=1,NBANDS )
      read(iounit,102) ( nfreqpts(nband), nband=1,NBANDS )
      read(iounit,103) ( endwvnbands(nband), nband=1,NBANDS )
      read(iounit,103) FIRSTRAYBAND,NIRBANDS
      read(iounit,104) ( powph2o(nband), nband=1,NH2OBANDS )
      read(iounit,104) ( p0h2o(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1co2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1co2str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2co2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2co2str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3co2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3co2str(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1o2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c1o2str(nband), nband=1,NH2OBANDS ),  &
                         c1o2strschrun
      read(iounit,105) ( c2o2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c2o2str(nband), nband=1,NH2OBANDS ), &
                         c2o2strschrun
      read(iounit,105) ( c3o2(nband), nband=1,NH2OBANDS )
      read(iounit,105) ( c3o2str(nband), nband=1,NH2OBANDS ),  &
                         c3o2strschrun
      do nf = 1,nfrqpts
        read(iounit,106) wtfreq(nf),kh2o(nf),ko3(nf),strterm (nf)
      end do

      read(iounit,107) nintsolar

      allocate ( nwvnsolar (nintsolar) )
      allocate ( solint    (nintsolar) )
 
      do ni = 1,nintsolar
        read(iounit,107) nwvnsolar (ni),solint(ni)
      end do
 
      call close_file (iounit)
 
!----------------------------------------------------------------------c
! define the one wavenumber solar fluxes.                              c
!----------------------------------------------------------------------c
 
      do ni = 1,nintsolar
 
	if ( ni.eq.1 ) then
	  nw1 = 1
        else
          nw1 = nw1 + nwvnsolar(ni-1)
        end if
 
	nw2 = nw1 + nwvnsolar(ni) - 1
 
	do nw = nw1,nw2
	  solarfluxtoa(nw) = solint(ni)
        end do
 
      end do

 
      deallocate  (solint    )
      deallocate  (nwvnsolar )
 
!----------------------------------------------------------------------
!   send the solar flux data to be stored in esfsw_parameters_mod.
!--------------------------------------------------------------------

      call put_solarfluxes(solarfluxtoa, solflxband, endwvnbands)

!----------------------------------------------------------------------
!   convert some values to mks to match model units
!--------------------------------------------------------------------

!     p0h2o = 1.0E+2*p0h2o  !  mb to mks
      p0h2o = 1.0E-2/p0h2o  ! invert, and convert mb to mks
      kh2o = kh2o *1.0E-01   ! cgs to mks
      ko3  = ko3 *1.0E-01    ! cgs to mks


!----------------------------------------------------------------------
 
      c4co2(:) = c1co2(:) * c2co2(:) ** c3co2(:)
      c4co2str(:) = c1co2str(:) * c2co2str(:) ** c3co2str(:)
      c4o2(:) = c1o2(:) * c2o2(:) ** c3o2(:)
      c4o2str(:) = c1o2str(:) * c2o2str(:) ** c3o2str(:)
      c4o2strschrun = c1o2strschrun * c2o2strschrun ** c3o2strschrun
 
      solflxtotal = 0.0
      do nband = 1,NBANDS
	solflxtotal = solflxtotal + solflxband(nband)
      end do
 
!----------------------------------------------------------------------c
! define the wavenumbers to evaluate rayleigh optical depth.           c
!----------------------------------------------------------------------c
 
      endwvnbands(0) = 0
      do nband = 1,NBANDS
	freqnu(nband) = 0.5 * ( endwvnbands(nband-1) +   &
                                endwvnbands(nband) )
      end do
 
!----------------------------------------------------------------------c
! define quantities used to determine the rayleigh optical depth.      c
!                                                                      c
! notes: refquanray is the quantity which multiplies pressure /        c
!        temperature to yield the molecular density.                   c
!                                                                      c
!        betaddensitymol is the quantity which multiples the           c
!        molecular density to yield the rayleigh scattering            c
!        coefficient.                                                  c
!                                                                      c
!        1.39E-02 is the depolorization factor.                        c
!----------------------------------------------------------------------c
 
      refquanray = densmolref * temprefray / pressrefray 
 
      corrfac = ( 6.0E+00 + 3.0E+00 * 1.39E-02 )/( 6.0E+00 - 7.0E+00 * &
                  1.39E-02 )
      gamma = 1.39E-02 / ( 2.0E+00 - 1.39E-02 )
      f1 = 7.5E-01 / ( gamma * 2.0E+00 + 1.0E+00 )
      f2 = gamma * 3.0E+00 + 1.0E+00 
      f3 = 1.0E+00 - gamma 
      pinteg = 2.0E+00 * pie * ( 2.0E+00 * f1 * f2 * ( 1.0E+00 + f3 / &
               f2 / 3.0E+00 ) )
      twopiesq = 2.0E+00 * pie ** 2 
      densmolrefsqt3 = 3.0E+00 * densmolref ** 2
 
      do nband = 1,NBANDS
        wavelength = 1.0E+04 / freqnu(nband)
        freqsq = 1.0 / ( wavelength ) ** 2
        ristdm1 = ( 6.4328E+03 + 2.94981E+06 / ( 1.46E+02 - freqsq ) + &
                    2.554E+04 / ( 4.1E+01 - freqsq ) ) * 1.0E-08
        ri = ristdm1 + 1
        betaddensitymol(nband) = twopiesq*( ri ** 2 - 1.0E+00 ) ** 2 * &
                                 corrfac / ( densmolrefsqt3 *  &
                                 wavelength ** 4 ) * pinteg * convfac 
      end do
 
!----------------------------------------------------------------------c
! define the gaussian angles for evaluation of the diffuse beam.       c
!----------------------------------------------------------------------c
 
      do i = 1,nstreams
        cosangstr(i) = ( ptstr(i) + 1. ) * 5.0E-01
      end do

!----------------------------------------------------------------------c
! deallocate local arrays.       
!----------------------------------------------------------------------c
      deallocate ( freqnu)
      deallocate (ptstr)
      deallocate (endwvnbands)

!--------------------------------------------------------------------
 101  format( 12f10.4 )
 102  format( 32i4 )
 103  format( 20i6 )
 104  format( 12f10.2 )
 105  format( 1p,16e8.1 )
 106  format( 1p,3e16.6,l16 )
 107  format( i5,1p,e14.5 )
 
!---------------------------------------------------------------------

end subroutine esfsw_driver_init 



!#################################################################

subroutine swresf (cloudfrac, cirrfgd, cvisrfgd, fracday,  &
                   press, qo3, rh2o, ssolar, temp,         &
                   dfsw, fsw, hsw, ufsw,                   &
                   cosangsolar, gausswt, nsolwg,           &
                   dfswclr, fswclr, hswclr, ufswclr)

!----------------------------------------------------------------------c
! swresf uses the delta-eddington technique in conjunction with a      c
! multiple-band parameterization for h2o+co2+o2+o3 absorption to       c
! derive solar fluxes and heating rates.                               c
!                                                                      c
! notes: drops are assumed if temp>273.15K, ice crystals otherwise.    c
!----------------------------------------------------------------------c
!                                                                      c
! intent in:                                                           c
!                                                                      c
! cloudfrac   = cloud amount                                           c
!                                                                      c
! cirrfgd     = the ground albedo for infrared wavelengths             c
!                                                                      c
! cosangsolar = cosine of the solar zenith angle, or gaussian values of 
!	        the cosine of the zenith angle, if lswg = true    
!                                                                      c
! cvisrfgd    = the ground albedo for visible and ultraviolet          c
!               wavelengths                                            c
!                                                                      c
! fracday     = fraction of day (or timestep) that sun is above the    c
!               horizon                                                c
!                                                                      c
! gausswt     = gaussian weights for zenith angles (if lswg is true); 
!               otherwise = 1.0
!
! nsolwg      = number of gaussian angles used; if lswg is false, nsolwg
!               is equal to 1
!
! press       = pressure at data levels of model                       c
!                                                                      c
! qo3         = mass mixing ratio of o3 at model data levels           c
!                                                                      c
! rh2o        = mass mixing ratio of h2o at model data levels          c
!                                                                      c
! ssolar      = solar constant in watts/meter**2 (may vary over one    c
!               year)                                                  c
!                                                                      c
! temp        = temperature in degrees K                               c
!                                                                      c
!        optional arguments
!
! sw_with_clouds = present if this is call made with clouds allowed;
!                  i.e., it is not a cloud - forcing diagnostic call
!
!----------------------------------------------------------------------c

real,dimension(:,:,:), intent(in)  :: cloudfrac, press, qo3, &
				      rh2o, temp,  cosangsolar
real, dimension(:,:  ), intent(in) :: cirrfgd, cvisrfgd, fracday
real, dimension(:),     intent(in) :: gausswt
real,                   intent(in) :: ssolar
integer,                intent(in) :: nsolwg
!logical,      intent(in), optional :: sw_with_clouds

!----------------------------------------------------------------------
! intent out:                                                          
!                                                                      
! dfsw =  downward flux at all pressure levels                    
!                                                                      
! fsw  =  net flux (up-down) at all pressure levels               
!                                                                      
! hsw  =  shortwave heating rates at all pressure layers               
!                                                                     
! ufsw =  upward flux at all pressure levels                      
!
!         optional arguments
!
! dfswcf =  clear-sky downward flux at all pressure levels
!                                                                       
! fswcf  =  clear-sky net flux (up-down) at all pressure levels
!                                                                       
! hswcf  =  clear-sky shortwave heating rates at all pressure layers
!                                                                       
! ufswcf =  clear-sky upward flux at all pressure levels
!
!----------------------------------------------------------------------c
 
real, dimension(:,:,:), intent(inout)            :: dfsw, fsw, hsw, ufsw
real, dimension(:,:,:), intent(inout) , optional :: dfswclr, fswclr,   &
                                                  hswclr, ufswclr
 
!-----------------------------------------------------------------------
!     local variables:
!-----------------------------------------------------------------------
 
   logical, dimension(:,:,:), allocatable :: cloud
   logical, dimension(:,:  ), allocatable :: cloud_in_column
   logical, dimension(:,:), allocatable   :: daylight

   real,    dimension(:,:,:), allocatable ::     &
         aeroasymfac,     aerosctopdep,   aeroextopdep,     &
         alphaco2,        alphaco2str,    alphao2,          &
         alphao2str,      cloudasymfac,   cloudextopdep,    &
         cloudsctopdep,   delpdig,        deltap,           &
                          densitymol,     extopdepclr,      &
	 extopdepovc,     fclr,           fovc,             &
         gocpdp,                                            &
         gasopdep,        gclr,           gstrclr,          &
         gstrovc,         govc,           omegastrclr,      &
         omegastrovc,     opdep,                            &
         rayopdep,        reflectance,    rlayerdif,        &
         rlayerdifclr,    rlayerdifovc,   rlayerdir,        &
         rlayerdirclr,    rlayerdirovc,   scale,            &
         scalestr,        sctopdepclr,    sctopdepovc,      &
         ssalbclr,        ssalbovc,       sumre,            &
         sumtr,           taustrclr,      taustrovc,        &
         tco2,            tlayerde,       tlayerdeclr,      &
         tlayerdeovc,     tlayerdif,      tlayerdifclr,     &
         tlayerdifovc,    tlayerdir,      tlayerdirclr,     &
         tlayerdirovc,    to2,            transmittance,    &
         wh2o,            wh2ostr,        wo3              

   real, dimension (:,:,:), allocatable ::      &
         reflectanceclr,  transmittanceclr,                 &
         sumtrclr,        sumreclr,  &
         dfswband, fswband, hswband, ufswband

   real, dimension (:,:,:,:), allocatable ::                &
                                          efftauo2,         &
	 efftauco2,                                         &
	 totco2,          totco2str,      toto2,            &
	 toto2str

   real, dimension (:,:,:), allocatable ::                &
         dfswbandclr,     fswbandclr,     hswbandclr,       &
         ufswbandclr

   real, dimension (:,:,:), allocatable ::                &
	 cldext,          cldsct,         cldasymm

   real, dimension (:,:,:), allocatable ::                &
	 aerext,          aerssalb,       aerasymm
   real, dimension (:,:,:), allocatable ::                  &
	 aeramtsc

   real, dimension (:,:,:), allocatable   :: pflux_mks

   real, dimension(:,:), allocatable       ::               &
	 sfcalb,          solarflux,      wtfac, denom

   integer            :: j, i, k, ng, ncld, ncldlyrs, nhiclds,    &
			 nmiclds, nloclds, np, nband, nf, ns
   character(len=15)  :: swaer_form


!--------------------------------------------------------------------
!    define limits and dimensions 
!--------------------------------------------------------------------

      IJK   = IERAD-ISRAD+1

!---------------------------------------------------------------------

if (Rad_control%do_totcld_forcing) then
  allocate ( reflectanceclr  (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
  allocate ( sumreclr        (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
  allocate ( sumtrclr        (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
  allocate ( transmittanceclr(ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
endif

allocate ( opdep           (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD) )
allocate ( reflectance     (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
allocate ( sfcalb          (ISRAD:IERAD,JSRAD:JERAD) )
allocate ( solarflux       (ISRAD:IERAD, JSRAD:JERAD) )
allocate ( sumre           (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
allocate ( sumtr           (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
allocate ( transmittance   (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD+1))
allocate ( wh2o            (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD) )
allocate ( wh2ostr         (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD))
allocate ( wo3             (ISRAD:IERAD, JSRAD:JERAD,KSRAD:KERAD) )
allocate ( wtfac           (ISRAD:IERAD, JSRAD:JERAD)  )

allocate ( densitymol      (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( cloud           (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( cloud_in_column (ISRAD:IERAD, JSRAD:JERAD               ) )
allocate ( daylight        (ISRAD:IERAD, JSRAD:JERAD               ) )
allocate ( deltap          (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( gocpdp          (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( denom           (ISRAD:IERAD, JSRAD:JERAD               ) )
allocate ( delpdig         (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( totco2    (ISRAD:IERAD, JSRAD:JERAD, KSRAD+1:KERAD+1, NSOLWG) )
allocate ( totco2str (ISRAD:IERAD, JSRAD:JERAD, KSRAD+1:KERAD+1, NSOLWG) )
allocate ( toto2     (ISRAD:IERAD, JSRAD:JERAD, KSRAD+1:KERAD+1, NSOLWG) )
allocate ( toto2str  (ISRAD:IERAD, JSRAD:JERAD, KSRAD+1:KERAD+1, NSOLWG) )
allocate ( rlayerdif       (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( rlayerdifclr    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( rlayerdifovc    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( rlayerdir       (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( rlayerdirclr    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( rlayerdirovc    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( tlayerdif       (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( tlayerdifclr    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( tlayerdifovc    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( tlayerdir       (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( tlayerdirclr    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( tlayerdirovc    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( tlayerde        (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( tlayerdeovc     (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
allocate ( tlayerdeclr     (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )

allocate ( cldext    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
allocate ( cldsct    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
allocate ( cldasymm  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )

allocate ( aerext    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
allocate ( aerssalb  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
allocate ( aerasymm  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
allocate ( aeramtsc  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )

allocate ( aeroextopdep    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
allocate ( aerosctopdep    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
allocate ( aeroasymfac     (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
allocate ( cloudextopdep   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
allocate ( cloudsctopdep   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
allocate ( cloudasymfac    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
allocate ( rayopdep        (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
allocate ( alphaco2        (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
allocate ( alphaco2str     (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
allocate ( alphao2         (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
allocate ( alphao2str      (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
allocate ( efftauco2 (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NSOLWG))
allocate ( efftauo2  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NSOLWG))
allocate ( tco2            (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
allocate ( to2             (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
allocate ( pflux_mks       (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))

!---------------------------------------------------------------------
! initialize local variables.                                        
!--------------------------------------------------------------------

        alphaco2   (:,:,1) = 0.0
        alphaco2str(:,:,1) = 0.0
        alphao2    (:,:,1) = 0.0
        alphao2str (:,:,1) = 0.0
	 
!----------------------------------------------------------------------c
! define a flag indicating columns in which there is sunshine during
! this radiation time step. define a flag indicating points with both 
! sunlight and cloud.      
!----------------------------------------------------------------------c

            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
                if ( fracday(i,j) /= 0.0) then
                  daylight(i,j) = .true.                 
                else
                  daylight(i,j) = .false.                
                endif     
                cloud_in_column(i,j) = .false.
            end do
          end do
        
           do k=KSRAD,KERAD
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
           if (cloudfrac(i,j,k) > 0.0 .and. daylight(i,j) )  then
	      cloud_in_column(i,j) = .true.
              cloud(i,j,k) = .true.
            else
              cloud(i,j,k) = .false.
              endif
              end do
            end do
          end do

!--------------------------------------------------------------------
! if aerosol parameters are to be determined by the model, not
! prescribed by observations, call scatpar to place into form usable
! by the shortwave algorithm
!-------------------------------------------------------------------
       swaer_form = Sw_control%swaerosol_form
       if (trim(swaer_form) == 'predicted') then 
         call scatpar 
       endif

!----------------------------------------------------------------------c
! define pressure related quantities, pressure is in mks units. 
!----------------------------------------------------------------------c
 
	pflux_mks = pflux*1.0E-1

	allocate ( scale      (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )
	allocate ( scalestr   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )

        do k = KSRAD+1,KERAD+1
          deltap(:,:,k-1) = pflux_mks(:,:,k) - pflux_mks(:,:,k-1)
          gocpdp(:,:,k-1) = radcon_mks/deltap(:,:,k-1)
          delpdig(:,:,k-1) = deltap(:,:,k-1)/ grav_mks
          scalestr(:,:,k) = pflux_mks(:,:,k) 
          scale(:,:,k) = scalestr(:,:,k)*pflux_mks(:,:,k)/pstd_mks
        end do
 
!----------------------------------------------------------------------c
! define the scaled and unscaled co2 and o2 pathlengths in centimeter- c
! atm, and the unscaled h2o and o3 amounts in kgrams/meter**2. 
! cm-atm needed as units because of c2co2 having those units.
!----------------------------------------------------------------------c
  
        do ng = 1,nsolwg
          denom(:,:) = 1.0/(grav_mks*rhoair*cosangsolar(:,:,ng)*2.0)
          do k = KSRAD+1,KERAD+1
            totco2(:,:,k,ng) = 1.0E+02*rrvco2*scale(:,:,k)*denom(:,:)       
            totco2str(:,:,k,ng) = 2.0E+02 * rrvco2 * scalestr(:,:,k) * &
                                  denom(:,:)
            toto2(:,:,k,ng) = 1.0E+02 * o2mixrat * scale(:,:,k) *  &
                              denom(:,:)
            toto2str(:,:,k,ng) = 2.0E+02*o2mixrat * scalestr(:,:,k) *  &
                                 denom(:,:)
          end do
        end do

        do k = KSRAD,KERAD
          wh2ostr(:,:,k) = rh2o(:,:,k) * delpdig(:,:,k)
          wo3(:,:,k)     = qo3(:,:,k) * delpdig(:,:,k)
        end do
  
!----------------------------------------------------------------------c
! define the molecular density for use in calculating the              c
! rayleigh optical depth (deltaz is in meters).                        c
!----------------------------------------------------------------------c
 
        do k = KSRAD,KERAD
          densitymol(:,:,k) = refquanray * press(:,:,k) / temp(:,:,k)
        end do
 
	deallocate ( scalestr     )
	deallocate ( scale        )
	deallocate (pflux_mks     )
 
if (Rad_control%do_totcld_forcing) then
  allocate ( dfswbandclr(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
  allocate (  fswbandclr(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1)) 
  allocate (  hswbandclr(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  )) 
  allocate ( ufswbandclr(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))   
endif

  allocate ( dfswband   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
  allocate (  fswband   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
  allocate (  hswband   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ))
  allocate ( ufswband   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
  allocate ( sctopdepclr    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( gclr           (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( fclr           (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( extopdepovc    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( govc           (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( fovc           (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( gasopdep       (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( extopdepclr    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( ssalbclr       (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( sctopdepovc    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( ssalbovc       (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( taustrclr      (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( omegastrclr    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( taustrovc      (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( omegastrovc    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( gstrclr        (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))
  allocate ( gstrovc        (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD))

!--------------------------------------------------------------------
!   allocate variables for use in deledd
!---------------------------------------------------------------------

  allocate (cld_index  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
  allocate ( gstr2          (IJK             ))
  allocate ( taustr2        (IJK             ))
  allocate ( omegastr2      (IJK             ))
  allocate ( cosangzk2      (IJK             ))
  allocate ( rlayerdir2     (IJK             ))
  allocate ( tlayerde2      (IJK             ))
  allocate ( tlayerdir2     (IJK             ))
  allocate ( sumr           (IJK             ))
  allocate ( sumt           (IJK             ))

!--------------------------------------------------------------------
!   allocate variables for use in adding
!---------------------------------------------------------------------

  allocate( radddowndif (ISRAD:IERAD, JSRAD:JERAD,KSRAD  :KERAD  ) )
  allocate( raddupdif   (ISRAD:IERAD, JSRAD:JERAD,KSRAD  :KERAD+1) )
  allocate( raddupdir   (ISRAD:IERAD, JSRAD:JERAD,KSRAD  :KERAD+1) )
  allocate( tlevel      (ISRAD:IERAD, JSRAD:JERAD,KSRAD  :KERAD+1) )
  allocate( tadddowndir (ISRAD:IERAD, JSRAD:JERAD,KSRAD  :KERAD  ) )
  allocate( dm1tl       (ISRAD:IERAD, JSRAD:JERAD,KSRAD+1:KERAD  ) )
  allocate( dm2tl       (ISRAD:IERAD, JSRAD:JERAD,KSRAD  :KERAD  ) )
  allocate( rdm2tl      (ISRAD:IERAD, JSRAD:JERAD,KSRAD  :KERAD  ) )
  allocate( alpp        (ISRAD:IERAD, JSRAD:JERAD,KSRAD+1:KERAD+1) )
  allocate( dm3         (ISRAD:IERAD, JSRAD:JERAD,KSRAD+1:KERAD+1) )
  allocate( dm3r        (ISRAD:IERAD, JSRAD:JERAD,KSRAD+1:KERAD+1) )
  allocate( dm3r1p      (ISRAD:IERAD, JSRAD:JERAD,KSRAD+1:KERAD+1) )

!----------------------------------------------------------------------c
! np is a counter for the pseudo-monochromatic frequency point number. c
!----------------------------------------------------------------------c
 
        np = 0
 
!----------------------------------------------------------------------c
! begin band loop                                                      c
!----------------------------------------------------------------------c
 
        do nband = 1,NBANDS
 
!----------------------------------------------------------------------c
! define the surface albedo (infrared value for infrared bands,        c
! visible value for the remaining bands).                              c
!----------------------------------------------------------------------c
 
          if ( nband.le.NIRBANDS ) then
              sfcalb(:,:) = cirrfgd(:,:)
          else
              sfcalb(:,:) = cvisrfgd(:,:)
          end if
 
!----------------------------------------------------------------------c
! obtain aerosol properties from esfsw_scattering
!----------------------------------------------------------------------c

          call get_swaerprops_from_scattering ( nband,   &
	                         aerext, aerssalb, aerasymm, aeramtsc)

!----------------------------------------------------------------------c
! define the local variables for the band values of aerosol and cloud  c
! single scattering parameters.                                        c
!                                                                      c
! note: the unit for the aerosol extinction is kilometer**(-1).        c
!----------------------------------------------------------------------c
 
          do k = KSRAD,KERAD
            do  j=JSRAD,JERAD
               do i=ISRAD,IERAD
            if (daylight(i,j) ) then
            aeroextopdep(i,j,k) = 1.0e-03*aerext(i,j,k)*   &
                               aeramtsc(i,j,k)*deltaz(i,j,k)
            aerosctopdep(i,j,k) = aerssalb(i,j,k) * &
                                  aeroextopdep(i,j,k)
            aeroasymfac(i,j,k) = aerasymm(i,j,k)
           endif
          end do
          end do
          end do
 

!----------------------------------------------------------------------c
! obtain cloud properties from cloudrad_package
!----------------------------------------------------------------------c

            call get_swcldprops_from_cloudrad (nband, cldext, cldsct, &
					       cldasymm)

            do k = KSRAD,KERAD
              do j=JSRAD,JERAD
                do i=ISRAD, IERAD
                 if (cloud(i,j,k) ) then
                   cloudextopdep(i,j,k) = 1.0E-03*cldext(i,j,k) *    &
                                   deltaz(i,j,k) 
                cloudsctopdep(i,j,k) = 1.0E-03*cldsct(i,j,k) *    &
                                   deltaz(i,j,k) 
                cloudasymfac(i,j,k) = cldasymm(i,j,k)

                endif
            end do
            end do
            end do
!----------------------------------------------------------------------c
! define the rayleigh optical depths.                                  c
!----------------------------------------------------------------------c
 
          rayopdep(:,:,:) = betaddensitymol(nband)*densitymol(:,:,:) * &
                            deltaz(:,:,:)
 
          if ( nband.le.NH2OBANDS ) then
 
!----------------------------------------------------------------------c
! define the h2o scaled gas amounts in kgrams/meter**2             c
!----------------------------------------------------------------------c
 
            do k = KSRAD,KERAD
              wh2o(:,:,k) = rh2o(:,:,k) * delpdig(:,:,k) *   &
                            exp( powph2o(nband)*alog(press(:,:,k)* &
                                 p0h2o(nband) ) )
            end do
 
!----------------------------------------------------------------------c
! calculate the "effective" co2 and o2 gas optical depths for the      c
! appropriate absorbing bands.                                         c
!                                                                      c
! note: a correction is applied to the determined transmissions for    c
!       co2 and o2 if t < 0, which can occur for large zenith angles.  c
!----------------------------------------------------------------------c
 
            if ( c1co2(nband).ne.1.0E-99 ) then
	      do ng = 1,NSOLWG
                do k = KSRAD+1,KERAD+1
		  do j=JSRAD,JERAD
		    do i = ISRAD,IERAD
		      if (daylight(i,j) ) then
                      alphaco2(i,j,k) = c1co2(nband)*exp(c3co2(nband)* &
                                        alog( ( totco2(i,j,k,ng) +  &
                                                c2co2(nband) ) ) )  - &
				        c4co2(nband)
                      alphaco2str(i,j,k) = c1co2str(nband) *   &
                                           exp( c3co2str(nband) *   &
                                           alog((totco2str(i,j,k,ng) + &
                                           c2co2str(nband) ) ) ) -   &
                                           c4co2str(nband)
                      tco2(i,j,k-1) = ( 1.0 - alphaco2(i,j,k) ) *   &
                                      ( 1.0 - alphaco2str(i,j,k) ) /  &
                                    ( ( 1.0 - alphaco2(i,j,k-1) ) *   &
                                    ( 1.0 - alphaco2str(i,j,k-1) ) )
                      if ( tco2(i,j,k-1).le.0.0 )     &
						tco2(i,j,k-1) = 1.0E-09
                      efftauco2(i,j,k-1,ng) = -cosangsolar(i,j,ng) *   &
                                               alog( tco2(i,j,k-1) )
                    else
		      efftauco2(i,j,k-1,ng) = 0.0
                    endif
	  	    end do
                  end do
                end do
              end do
	    else
              efftauco2(:,:,:,:) = 0.0
	    end if
 
            if ( c1o2(nband).ne.1.0E-99 ) then
              do ng = 1,NSOLWG
                do k = KSRAD+1,KERAD+1
		  do j=JSRAD,JERAD
                    do i = ISRAD,IERAD
		      if (daylight(i,j) ) then
                      alphao2(i,j,k) = c1o2(nband)*exp( c3o2(nband) * &
                                       alog( ( toto2(i,j,k,ng) +  &
                                               c2o2(nband) ) ) ) -  &
				       c4o2(nband)
                      alphao2str(i,j,k) = c1o2str(nband) *  &
                                          exp( c3o2str(nband) *  &
                                          alog( ( toto2str(i,j,k,ng) + &
                                          c2o2str(nband) ) ) ) -  &
                                          c4o2str(nband)
                      to2(i,j,k-1) = ( 1.0 - alphao2(i,j,k) ) *   &
                                     ( 1.0 - alphao2str(i,j,k) ) /  &
                                   ( ( 1.0 - alphao2(i,j,k-1) ) *  &
                                   ( 1.0 - alphao2str(i,j,k-1) ) )
                      if ( to2(i,j,k-1).le.0.0 ) to2(i,j,k-1) = 1.0E-09
                      efftauo2(i,j,k-1,ng) = -cosangsolar(i,j,ng) *  &
                                              alog( to2(i,j,k-1) )
                    else
		      efftauo2(i,j,k-1,ng) = 0.0
                    endif
                    end do
                  end do
                end do
	      end do
            else
              efftauo2(:,:,:,:) = 0.0
	    end if
 
          end if
 
!----------------------------------------------------------------------c
! calculate the "effective" o2 gas optical depths for the Schuman-     c
! Runge band.                                                          c
!----------------------------------------------------------------------c

	  if ( nband.EQ.NBANDS ) then
            do ng = 1,NSOLWG
              do k = KSRAD+1,KERAD+1
	        do j=JSRAD,JERAD
                  do i = ISRAD,IERAD
		      if (daylight(i,j) ) then
                    alphao2str(i,j,k) = c1o2strschrun*exp  &
					(c3o2strschrun* &
                                         alog ( ( toto2str(i,j,k,ng) + &
                                         c2o2strschrun ) ) ) -    &
					 c4o2strschrun
                    to2(i,j,k-1) = ( 1.0 - alphao2str(i,j,k  ) ) /  &
                                   ( 1.0 - alphao2str(i,j,k-1) )
	            if ( to2(i,j,k-1).le.0.0 ) to2(i,j,k-1) = 1.0E-09
                    efftauo2(i,j,k-1,ng) = -cosangsolar(i,j,ng) * &
                                           alog( to2(i,j,k-1) )
                    else
		      efftauo2(i,j,k-1,ng) = 0.0
                    endif
                  end do
                end do
              end do
	    end do
	  end if
 
 !-----------------------------------------------------------------
 !    initialize summing arrays
 !-----------------------------------------------------------------

          sumtr(:,:,:) = 0.0
           sumre(:,:,:) = 0.0

          if (Rad_control%do_totcld_forcing) then
            sumtrclr(:,:,:) = 0.0
            sumreclr(:,:,:) = 0.0
          endif

 !-----------------------------------------------------------------
 !    define clear sky arrays
 !-----------------------------------------------------------------

          if (nband >= FIRSTRAYBAND .or. Sw_control%do_swaerosol ) then
            do k=ksrad,kerad
            do j=jsrad,jerad
            do i=israd,ierad
              if (daylight(i,j) ) then
            sctopdepclr(i,j,k) = rayopdep(i,j,k) + aerosctopdep(i,j,k)
            gclr(i,j,k) = aeroasymfac(i,j,k)*aerosctopdep(i,j,k) /  &
                          sctopdepclr(i,j,k)
            fclr(i,j,k) = aeroasymfac(i,j,k)*gclr(i,j,k)
            gstrclr(i,j,k) = ( gclr(i,j,k)  - fclr(i,j,k) ) /   &
                             ( 1.0 - fclr(i,j,k) )
	      endif
            end do
            end do
            end do
          endif

 !-----------------------------------------------------------------
 !    define cloudy sky arrays
 !-----------------------------------------------------------------

            do k=KSRAD,KERAD
              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
                  if (cloud(i,j,k)  )then
                    sctopdepovc(i,j,k) = rayopdep(i,j,k) +    &
					 aerosctopdep(i,j,k) + &
                                         cloudsctopdep(i,j,k) 
                    govc(i,j,k) = ( ( cloudasymfac(i,j,k) *   &
                                   cloudsctopdep(i,j,k) ) +  &
                                  ( aeroasymfac(i,j,k) *   &
                                   aerosctopdep(i,j,k)))/   &
				   sctopdepovc(i,j,k)
                    fovc(i,j,k) = ( ( cloudasymfac(i,j,k) ** 2 *  &
                                  cloudsctopdep(i,j,k) ) + &
                                  ( aeroasymfac(i,j,k) ** 2 *  &
                                  aerosctopdep(i,j,k) ))/  &
				  sctopdepovc(i,j,k)
                    gstrovc(i,j,k) = ( govc(i,j,k)  - fovc(i,j,k))/  &
                                     ( 1.0 - fovc(i,j,k) )
                  endif
                end do
              end do
            end do

!----------------------------------------------------------------------c
! begin frequency points in the band loop                              c
!----------------------------------------------------------------------c

          do nf = 1,nfreqpts(nband)
            np = np + 1
 
!----------------------------------------------------------------------c
! define the h2o + o3 gas optical depths.                              c
!----------------------------------------------------------------------c
 
            if ( strterm(np) ) then
              opdep(:,:,:) = kh2o(np)*wh2ostr(:,:,:) +    &
			     ko3(np)*wo3(:,:,:)
	    else
              opdep(:,:,:) = kh2o(np) * wh2o(:,:,:) +    &
			     ko3(np) * wo3(:,:,:)
            end if

!----------------------------------------------------------------------c
! begin gaussian angle loop (ng > 1 only when lswg = true).            c
!----------------------------------------------------------------------c
 
            do ng = 1,NSOLWG
              gasopdep(:,:,:) = opdep(:,:,:) + efftauco2(:,:,:,ng) +   &
                                efftauo2(:,:,:,ng) 
!----------------------------------------------------------------------c
! clear sky mode                                                       c
!                                                                      c
! note: in this mode, the delta-eddington method is performed for all  c
!       spatial columns experiencing sunlight.              
!----------------------------------------------------------------------c

!----------------------------------------------------------------------c
! calculate the scaled single-scattering quantities for use in the     c
! delta-eddington routine.                                             c
!----------------------------------------------------------------------c
              if (nband >= FIRSTRAYBAND .or.     &
		                       Sw_control%do_swaerosol )  then

                do k=ksrad,kerad
                do j=jsrad,jerad
                do i=israd,ierad
                  if (daylight(i,j) ) then
                extopdepclr(i,j,k) = gasopdep(i,j,k) +    &
				     rayopdep(i,j,k) +    &
				     aeroextopdep(i,j,k)
                ssalbclr(i,j,k) = sctopdepclr(i,j,k)/extopdepclr(i,j,k)
 
 
                taustrclr(i,j,k) = extopdepclr(i,j,k) * ( 1.0 -    &
                                   ssalbclr(i,j,k) * fclr(i,j,k) )
                omegastrclr(i,j,k) = ssalbclr(i,j,k)*((1.0 -     &
				     fclr(i,j,k))/(1.0 -   &
				     ssalbclr(i,j,k)*fclr(i,j,k)))
                  endif
                end do
                end do
                end do
!----------------------------------------------------------------------c
! calculate the reflection and transmission in the scattering layers   c
! using the delta-eddington method.                                    c
!----------------------------------------------------------------------c
                call deledd (taustrclr, omegastrclr, gstrclr,   &
			     cosangsolar, ng , daylight,rlayerdirclr,  &
			     tlayerdirclr, rlayerdifclr,  &
                             tlayerdifclr, tlayerdeclr)

                if (ng /= 1) then
                  tlayerdifclr(:,:,:) = 0.0       
		  if (NSTREAMS == 1) then
                    tlayerdifclr(:,:,:) = exp( -gasopdep(:,:,:)/ &
                                               cosangstr(1) )
                  else
                    do ns = 1,NSTREAMS
                      tlayerdifclr(:,:,:) = tlayerdifclr(:,:,:) +    &
	    				    exp( -gasopdep(:,:,:)/&
                                            cosangstr(ns) )*wtstr(ns)* &
                                            cosangstr(ns)
                    end do
		  endif
                  rlayerdifclr(:,:,:) = 0.0
                endif
		  
!----------------------------------------------------------------------c
! initialize the layer reflection and transmission arrays with the     c
! non-scattering case.                                                 c
!----------------------------------------------------------------------c
              else
                tlayerdifclr(:,:,:) = 0.0       
		if (NSTREAMS == 1) then
		  tlayerdifclr(:,:,:) = exp( -gasopdep(:,:,:)/ &
		                             cosangstr(1) )
		else
                  do ns = 1,NSTREAMS
                    tlayerdifclr(:,:,:) = tlayerdifclr(:,:,:) +    &
			   		  exp( -gasopdep(:,:,:)/&
                                          cosangstr(ns) )*wtstr(ns)*  & 
                                          cosangstr(ns)
                  end do
		endif
                rlayerdirclr(:,:,:) = 0.0
	        do k=KSRAD,KERAD
                  tlayerdirclr(:,:,k) = exp( -gasopdep(:,:,k) /   &
                                            cosangsolar(:,:,ng) )
	        end do				  
                tlayerdeclr(:,:,:) = tlayerdirclr(:,:,:)
                rlayerdifclr(:,:,:) = 0.0
              endif

!----------------------------------------------------------------------c
! overcast sky mode                                                    c
!                                                                      c
! note: in this mode, the delta-eddington method is performed only for c
!       spatial columns containing a cloud and experiencing sunlight. 
!----------------------------------------------------------------------c

!----------------------------------------------------------------------c
! calculate the scaled single-scattering quantities for use in the     c
! delta-eddington routine.                                             c
!----------------------------------------------------------------------c

                do k=KSRAD,KERAD
                  do j=JSRAD,JERAD
                    do i=ISRAD,IERAD
                      if (cloud(i,j,k) ) then
                        extopdepovc(i,j,k) = gasopdep(i,j,k) +    &
					     rayopdep(i,j,k) +  &
                                             aeroextopdep(i,j,k) +  &
					     cloudextopdep(i,j,k)
                        ssalbovc(i,j,k) = sctopdepovc(i,j,k) /    &
					  extopdepovc(i,j,k)
                        taustrovc(i,j,k) = extopdepovc(i,j,k)*( 1.0 - &
                                           ssalbovc(i,j,k)*fovc(i,j,k) )
                        omegastrovc(i,j,k) = ssalbovc(i,j,k)*( ( 1.0 - &
					     fovc(i,j,k) )/( 1.0 -   &
					     ssalbovc(i,j,k) *   &
					     fovc(i,j,k) ) )
 
                      endif
                    end do
                  end do
                end do
!----------------------------------------------------------------------c
! calculate the reflection and transmission in the scattering layers   c
! using the delta-eddington method.                                    c
!----------------------------------------------------------------------c
                call deledd (taustrovc, omegastrovc, gstrovc,   &
			     cosangsolar, ng, daylight, rlayerdirovc, &
			     tlayerdirovc, rlayerdifovc,  & 
                             tlayerdifovc, tlayerdeovc, cloud)
                if (ng /= 1) then
		  tlayerdifovc(:,:,:) = tlayerdifclr(:,:,:)
		  rlayerdifovc(:,:,:) = rlayerdifclr(:,:,:)
		endif
 
!----------------------------------------------------------------------c
! weight the reflection and transmission arrays for clear and          c
! overcast sky conditions by the cloud fraction, to calculate the      c
! resultant values.                                                    c
!----------------------------------------------------------------------c
 
              do k=KSRAD,KERAD
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
		    if ( cloud(i,j,k) ) then
                      rlayerdir(i,j,k) = cloudfrac(i,j,k)*   &
	  		                 rlayerdirovc(i,j,k) +  &
                                         (1.0 - cloudfrac(i,j,k) ) *   &
                                         rlayerdirclr(i,j,k)
                      rlayerdif(i,j,k) = cloudfrac(i,j,k) *  &
					 rlayerdifovc(i,j,k) +  &
                                         ( 1.0 - cloudfrac(i,j,k) ) *  &
                                         rlayerdifclr(i,j,k)
                      tlayerdir(i,j,k) = cloudfrac(i,j,k) *   &
					 tlayerdirovc(i,j,k) +  &
                                         ( 1.0 - cloudfrac(i,j,k) ) *  &
                                         tlayerdirclr(i,j,k)
                      tlayerdif(i,j,k) = cloudfrac(i,j,k) *   &
				 	 tlayerdifovc(i,j,k) +  &
                                         ( 1.0 - cloudfrac(i,j,k) ) *  &
                                         tlayerdifclr(i,j,k)
                      tlayerde(i,j,k) = cloudfrac(i,j,k) *   &
				        tlayerdeovc(i,j,k) +  &
                                        (1.0 - cloudfrac(i,j,k) ) *  &
                                        tlayerdeclr(i,j,k)
                    else if (daylight(i,j)) then
	              rlayerdir(i,j,k) = rlayerdirclr(i,j,k)
	              tlayerdir(i,j,k) = tlayerdirclr(i,j,k)
	              rlayerdif(i,j,k) = rlayerdifclr(i,j,k)
	              tlayerdif(i,j,k) = tlayerdifclr(i,j,k)
	              tlayerde (i,j,k) = tlayerdeclr (i,j,k)
                    endif
          	  end do
                end do
	      end do
 
!----------------------------------------------------------------------c
! calculate the reflection and transmission at flux levels from the    c
! direct and diffuse values of reflection and transmission in the      c
! corresponding layers using the adding method.                        c
!----------------------------------------------------------------------c
 
              call adding (rlayerdir, tlayerdir, rlayerdif, tlayerdif, &
                           tlayerde , sfcalb   , daylight, &
                           reflectance, transmittance)    
 
              if (Rad_control%do_totcld_forcing) then

                call adding (                                    &
                             rlayerdirclr, tlayerdirclr,         &
                             rlayerdifclr, tlayerdifclr,         &
                             tlayerdeclr,  sfcalb, cloud_in_column, &
                             reflectanceclr, transmittanceclr)

              endif

!----------------------------------------------------------------------c
! weight and sum the reflectance and transmittance to calculate the    c
! band values.                                                         c
!----------------------------------------------------------------------c
 
              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
              wtfac(i,j) = wtfreq(np)*gausswt(ng)*cosangsolar(i,j,ng)
                   end do
                end do
 
              do k = KSRAD,KERAD+1
                do j=JSRAD,JERAD
                 do i=ISRAD,IERAD
                   if (daylight(i,j) ) then
                sumtr(i,j,k) = sumtr(i,j,k) + transmittance(i,j,k)*   &
                              wtfac(i,j)
                sumre(i,j,k) = sumre(i,j,k) + reflectance(i,j,k)* &
                           wtfac(i,j)
                 endif
              end do
              end do
              end do

              if (Rad_control%do_totcld_forcing) then

                do k = KSRAD,KERAD+1
                do j=JSRAD,JERAD
                 do i=ISRAD,IERAD
                  if (cloud_in_column(i,j)       ) then
                  sumtrclr(i,j,k) = sumtrclr(i,j,k) +         &
                                    transmittanceclr(i,j,k)*wtfac(i,j) 
                  sumreclr(i,j,k) = sumreclr(i,j,k) +         &
                                    reflectanceclr(i,j,k)*wtfac(i,j)
                   else if (daylight(i,j) ) then
          sumtrclr(i,j,k) = sumtrclr(i,j,k) + transmittance(i,j,k)*   &
                              wtfac(i,j)
             sumreclr(i,j,k) = sumreclr(i,j,k) + reflectance(i,j,k)* &
                           wtfac(i,j)
                 endif
              end do
              end do
                end do
              endif

!--------------------------------------------------------------------

            end do    ! end of gaussian loop
          end do  ! end of frequency points in the band loop
 
!----------------------------------------------------------------------c
! normalize the solar flux in the band to the appropriate value for    c
! the given total solar insolation.                                    c
!----------------------------------------------------------------------c
 
          solarflux(:,:) = fracday(:,:) * solflxband(nband) * ssolar / &
                           solflxtotal
 
!----------------------------------------------------------------------c
! calculate the band fluxes and heating rates.                         c
!----------------------------------------------------------------------c
 
         if (do_esfsw_band_diagnostics) then
          do k = KSRAD,KERAD+1
           do j=JSRAD,JERAD
             do i=ISRAD,IERAD
            dfswband(i,j,k      ) = sumtr(i,j,k) * solarflux(i,j) 
            ufswband(i,j,k      ) = sumre(i,j,k) * solarflux(i,j)
          end do
          end do
          end do
         endif
 
!----------------------------------------------------------------------c
! sum the band fluxes and heating rates to calculate the total         c
! spectral values.                                                     c
!----------------------------------------------------------------------c
          do k = KSRAD,KERAD+1
           do j=JSRAD,JERAD
             do i=ISRAD,IERAD
              if (daylight(i,j) ) then
          dfsw (i,j,k) = dfsw(i,j,k) + sumtr(i,j,k)*solarflux(i,j)
          ufsw (i,j,k) = ufsw(i,j,k) + sumre(i,j,k)*solarflux(i,j)
          fswband(i,j,k      ) = (sumre(i,j,k)*solarflux(i,j)) - &
                                 (sumtr(i,j,k)*solarflux(i,j))
          fsw(i,j,k) = fsw(i,j,k) + fswband(i,j,k)
            endif
          end do
          end do
          end do
 
          do k = KSRAD,KERAD
           do j=JSRAD,JERAD
             do i=ISRAD,IERAD
              if (daylight(i,j) ) then
            hswband(i,j,k      ) = (fswband(i,j,k+1      ) -    &
                                    fswband(i,j,k      ) )*gocpdp(i,j,k)
            hsw(i,j,k) = hsw(i,j,k) + hswband(i,j,k      )
            endif
          end do
          end do
          end do

          if (Rad_control%do_totcld_forcing) then

!----------------------------------------------------------------------c
! calculate the band fluxes and heating rates.                         c
!----------------------------------------------------------------------c
         if (do_esfsw_band_diagnostics) then
            do k = KSRAD,KERAD+1
           do j=JSRAD,JERAD
             do i=ISRAD,IERAD
              dfswbandclr(i,j,k      ) = sumtrclr(i,j,k)*solarflux(i,j)
              ufswbandclr(i,j,k      ) = sumreclr(i,j,k)*solarflux(i,j)
          end do
          end do
            end do
	    endif

!----------------------------------------------------------------------c
! sum the band fluxes and heating rates to calculate the total         c
! spectral values.                                                     c
!----------------------------------------------------------------------c
            do k = KSRAD,KERAD+1
           do j=JSRAD,JERAD
             do i=ISRAD,IERAD
              dfswclr(i,j,k) = dfswclr(i,j,k) +sumtrclr(i,j,k)* &
                               solarflux(i,j)
              ufswclr(i,j,k) = ufswclr(i,j,k) + sumreclr   (i,j,k)*  &
                                 solarflux(i,j)
              fswbandclr(i,j,k  ) = (sumreclr(i,j,k)*solarflux(i,j))- &
                                      (sumtrclr(i,j,k)*solarflux(i,j))
               fswclr(i,j,k) = fswclr(i,j,k) + fswbandclr(i,j,k      )
          end do
          end do
            end do
           

            do k = KSRAD,KERAD
           do j=JSRAD,JERAD
             do i=ISRAD,IERAD
              hswbandclr(i,j,k      ) = (fswbandclr(i,j,k+1      ) -  &
                                         fswbandclr(i,j,k      ) ) *  &
                                         gocpdp(i,j,k)
               hswclr(i,j,k) = hswclr(i,j,k) + hswbandclr(i,j,k      )
          end do
          end do
            end do
          endif

!----------------------------------------------------------------------c
! sum the band fluxes and heating rates to calculate the total         c
! spectral values.                                                     c
!----------------------------------------------------------------------c
 

        end do      ! end of band loop
 
!-----------------------------------------------------------------
!   deallocate arrays used in adding
!-----------------------------------------------------------------

      deallocate  (  dm3r1p     )
      deallocate  (  dm3r       )
      deallocate  (  dm3        )
      deallocate  (  alpp       )
      deallocate  (  rdm2tl     )
      deallocate  (  dm2tl      )
      deallocate  (  dm1tl      )
      deallocate  (  tadddowndir)
      deallocate  (  tlevel     )
      deallocate  (  raddupdir  )
      deallocate  (  raddupdif  )
      deallocate  (  radddowndif)

!-----------------------------------------------------------------
!   deallocate arrays used in deledd
!-----------------------------------------------------------------

      deallocate  (  sumt      )
      deallocate  (  sumr      )
      deallocate  (  tlayerdir2)
      deallocate  (  tlayerde2 )
      deallocate  (  rlayerdir2)
      deallocate  (  cosangzk2 )
      deallocate  (  omegastr2 )
      deallocate  (  taustr2   )
      deallocate  (  gstr2     )
      deallocate  (cld_index   )

!-----------------------------------------------------------------
!   deallocate swresf arrays 
!-----------------------------------------------------------------

      deallocate (gstrovc           )
      deallocate (gstrclr           )
      deallocate (omegastrovc       )
      deallocate (taustrovc         )
      deallocate (omegastrclr       )
      deallocate (taustrclr         )
      deallocate (ssalbovc          )
      deallocate (sctopdepovc       )
      deallocate (ssalbclr          )
      deallocate (extopdepclr       )
      deallocate (gasopdep          )
      deallocate (fovc              )
      deallocate (govc              )
      deallocate (extopdepovc       )
      deallocate (fclr              )
      deallocate (gclr              )
      deallocate (sctopdepclr       )

      deallocate (dfswband          )
      deallocate ( fswband          )
      deallocate ( hswband          )
      deallocate (ufswband          )
 
      if (Rad_control%do_totcld_forcing) then
        deallocate (dfswbandclr       )
        deallocate ( fswbandclr       )
        deallocate ( hswbandclr       )
        deallocate (ufswbandclr       )
      endif

      deallocate (to2               )
      deallocate (tco2              )
      deallocate (efftauo2          )
      deallocate (efftauco2          )
      deallocate (alphao2str         )
      deallocate (alphao2            )
      deallocate (alphaco2str         )
      deallocate (alphaco2            )

      deallocate (rayopdep          )
      deallocate (cloudasymfac      )
      deallocate (cloudsctopdep     )
      deallocate (cloudextopdep     )
      deallocate (aeroasymfac       )
      deallocate (aerosctopdep      )
      deallocate (aeroextopdep      )
      deallocate (tlayerdeovc       )
      deallocate (tlayerdeclr       )
      deallocate (tlayerde          )
      deallocate (tlayerdir         )
      deallocate (tlayerdirclr      )
      deallocate (tlayerdirovc      )
      deallocate (tlayerdif         )
      deallocate (tlayerdifclr      )
      deallocate (tlayerdifovc      )
      deallocate (rlayerdir         )
      deallocate (rlayerdirclr      )
      deallocate (rlayerdirovc      )
      deallocate (rlayerdif         )
      deallocate (rlayerdifclr      )
      deallocate (rlayerdifovc      )
      deallocate (totco2      )
      deallocate (totco2str   )
      deallocate (toto2      )
      deallocate (toto2str   )
      deallocate ( delpdig        )
      deallocate ( deltap         )
      deallocate ( gocpdp         )
      deallocate ( denom          )
      deallocate ( densitymol     )
      deallocate ( cloud          )
      deallocate ( cloud_in_column)
      deallocate ( daylight       )

      deallocate (cldext          )
      deallocate (cldasymm        )
      deallocate (cldsct          )

      deallocate (aerext          )
      deallocate (aerasymm        )
      deallocate (aerssalb        )
      deallocate (aeramtsc        )

      deallocate (  wtfac         )
      deallocate (  wo3           )
      deallocate (  wh2ostr       )
      deallocate (  wh2o          )
      deallocate (  transmittance )
      deallocate (  sumtr         )
      deallocate (  sumre         )
      deallocate (  solarflux     )
      deallocate (  sfcalb        )
      deallocate (  reflectance   )
      deallocate (  opdep         )

      if (Rad_control%do_totcld_forcing) then
        deallocate (  transmittanceclr )
        deallocate (  sumtrclr         )
        deallocate (  sumreclr         )
        deallocate (  reflectanceclr   )
      endif

!---------------------------------------------------------------------


end  subroutine swresf



!#####################################################################

subroutine adding (rlayerdir, tlayerdir, rlayerdif, tlayerdif,  &
                   tlayerde, sfcalb, calc_flag, reflectance,   &
                   transmittance)
 
!----------------------------------------------------------------------c
! calculate the reflection and transmission at flux levels from the    c
! direct and diffuse values of reflection and transmission in the      c
! corresponding layers using the adding method.                        c
!                                                                      c
! references:                                                          c
!                                                                      c
! bowen, m.m., and v. ramaswamy, effects of changes in radiatively   
!        active species upon the lower stratospheric temperatures.,    c
!        j. geophys. res., 18909-18921, 1994.                          c
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!                                                                      c
! intent in:                                                           c
!                                                                      c
! rlayerdir = the layer reflectivity to a direct incident beam         c
!                                                                      c
! tlayerdir = the layer transmissivity to a direct incident beam       c
!                                                                      c
! rlayerdif = the layer reflectivity to a diffuse incident beam        c
!                                                                      c
! tlayerdif = the layer transmissivity to a diffuse incident beam      c
!                                                                      c
! tlayerde  = the layer transmissivity (non-scattered) to the direct   c
!             incident beam                                            c
!                                                                      c
! sfcalb    = surface albedo                                           c
!
! calc_flag = flag to indicate columns where adding is to be done.
!             calculations not done in "dark" columns and on clr sky
!             pass in columns without any clouds
!
!----------------------------------------------------------------------c
 
real, dimension(:,:,:),  intent(in)  :: rlayerdir, rlayerdif, &
                                        tlayerdir, tlayerdif, & 
                                        tlayerde
real, dimension (:,:) ,  intent(in)  :: sfcalb
logical, dimension (:,:) ,  intent(in)  :: calc_flag
!----------------------------------------------------------------------c
!                                                                      c
! intent out:                                                          c
!                                                                      c
! reflectance   = the reflectance of the scattered radiation at a      c
!                 level                                                c
!                                                                      c
! transmittance = the transmittance of the scattered radiation at a    c
!                 level                                                c
!----------------------------------------------------------------------c
 
real, dimension(:,:,:)  , intent(out)  :: reflectance, transmittance

!----------------------------------------------------------------------c
! local variables:                                                     c
!----------------------------------------------------------------------c
 
      integer     ::  k, j, i
      real, dimension (lbound(rlayerdir,3):ubound(rlayerdir,3)+1) ::  &
                                 raddupdif2, raddupdir2,  tlevel2
                                              

      real, dimension (lbound(rlayerdir,3):ubound(rlayerdir,3)  ) ::  &
                                      radddowndif2,  tadddowndir2, &
				      tlayerdif2, tlayerdir2, &
				      rlayerdif2, rlayerdir2, &
				      tlayerde2

      real :: dm1tl2, dm2tl2, rdm2tl2, dm32, dm3r2, dm3r1p2, alpp2, &
              raddupdif2p, raddupdir2p, tlevel2p, radddowndifm, &
	      tadddowndirm

!----------------------------------------------------------------------c
! initialization for the surface layer.                                c
!----------------------------------------------------------------------c
 
      do j=jsrad,jerad
      do i=israd,ierad

        if (calc_flag(i,j) ) then
      

 
!----------------------------------------------------------------------c
! add the inhomogeneous layers upward from the surface to the top of   c
! the atmosphere.                                                      c
!                                                                      c
! radiation incident from above for diffuse beam, reflection of        c
! direct beam and conversion to diffuse.                               c
!----------------------------------------------------------------------c
 
        raddupdif2p = sfcalb(i,j)
        raddupdir2p = sfcalb(i,j)
      do k = KERAD,KSRAD,-1
        dm2tl2    = tlayerdif(i,j,k) / ( 1.0 - rlayerdif(i,j,k) *  &
                     raddupdif2p )
        rdm2tl2    = dm2tl2    * raddupdif2p     
        raddupdif2(k) = rlayerdif(i,j,k) + tlayerdif(i,j,k)*   &
                        rdm2tl2    
 
        raddupdir2(k) = rlayerdir(i,j,k) + tlayerde(i,j,k) *   &
                        raddupdir2p     * dm2tl2    +   &     
                        (tlayerdir(i,j,k) - tlayerde(i,j,k)) *  &
                        rdm2tl2   
        raddupdir2p = raddupdir2(k)
        raddupdif2p = raddupdif2(k)
      end do
 
!----------------------------------------------------------------------c
! define the direct transmittance.                                     c
! add the inhomogeneous layers downward from the second layer to the   c
! surface.                                                             c
!                                                                      c
! radiation incident from below for diffuse beam, transmission of      c
! direct beam and conversion to diffuse.                               c
!----------------------------------------------------------------------c
 
!----------------------------------------------------------------------c
! initialization for the first scattering layer.                       c
!----------------------------------------------------------------------c
 
      tlevel2p         = tlayerde(i,j,KSRAD)
      radddowndifm    =  rlayerdif(i,j,KSRAD)
      tadddowndirm    =  tlayerdir(i,j,KSRAD)


      do k = KSRAD+1,KERAD
        dm1tl2       = tlayerdif(i,j,k) / ( 1.0 - rlayerdif(i,j,k) *  &
                     radddowndifm      )
        radddowndif2(k) = rlayerdif(i,j,k) + radddowndifm      * &
                           tlayerdif(i,j,k) * dm1tl2      
        tadddowndir2(k) = tlevel2p   * ( tlayerdir(i,j,k) + &
                           rlayerdir(i,j,k) * radddowndifm      * &
                           dm1tl2    ) + ( tadddowndirm      -  &
                           tlevel2p   ) * dm1tl2                    

!----------------------------------------------------------------------c
! add downward to calculate the resultant reflectances and             c
! transmittances at flux levels.                                       c
!----------------------------------------------------------------------c
        dm32    = 1.0/(1.0 - raddupdif2(k)*radddowndifm     )
        dm3r2    = dm32    * radddowndifm      
        dm3r1p2    = 1.0 + raddupdif2(k) * dm3r2   
        alpp2    = (tadddowndirm      - tlevel2p   )*dm32   
        transmittance(i,j,k) = (tlevel2p  *(1.0 + raddupdir2(k)* &
                                dm3r2   ) + alpp2   )
        reflectance(i,j,k) = ( tlevel2p   * raddupdir2(k) *   &
                             dm3r1p2    + alpp2    *   &
                             raddupdif2(k) )
        tlevel2p     = tlevel2p   * tlayerde (i,j,k) 
        radddowndifm = radddowndif2(k)
        tadddowndirm = tadddowndir2(k)
      end do
 
        dm32          = 1.0/(1.0 - sfcalb(i,j)*   &
                              radddowndifm       )
        dm3r2          = dm32          * radddowndifm       
        dm3r1p2          = 1.0 + sfcalb(i,j) * dm3r2         
        alpp2          = (tadddowndirm        - tlevel2p         )* &
                           dm32          
        transmittance(i,j,kERAD+1) = (tlevel2p        *(1.0 +   &
	                           sfcalb(i,j)* &
                                dm3r2         ) + alpp2         )
        reflectance(i,j,kERAD+1) = ( tlevel2p         *    &
	                              sfcalb(i,j) *   &
                             dm3r1p2          + alpp2          *   &
                             sfcalb(i,j) )

      reflectance(i,j,KSRAD) = raddupdir2p         
      transmittance(i,j,KSRAD) = 1.0


        endif
      end do
      end do
!------------------------------------------------------------------

end subroutine adding 



!####################################################################

subroutine deledd (taustr, omegastr, gstr, cosang, ng , daylight,  &
                    rlayerdir, tlayerdir, rlayerdif, tlayerdif,   &
		    tlayerde,  cloud)
 
!----------------------------------------------------------------------c
! calculate the reflection and transmission in the scattering layers   c
! using the delta-eddington method.                                    c
!                                                                      c
! references:                                                          c
!                                                                      c
! joseph, j.h., w. wiscombe, and j.a. weinman, the delta-eddington    
!      approximation for radiative flux transfer.,j. atmos. sci.,33,   c
!      2452-2459, 1976.                                                c
!----------------------------------------------------------------------c
!                                                                      c
! intent in:                                                           c
!                                                                      c
! taustr   = the scaled extinction optical depth                       c
!                                                                      c
! omegastr = the scaled single-scattering albedo                       c
!                                                                      c
! gstr     = the scaled asymmetry factor                               c
!                                                                      c
! cosang   = the cosine of the solar zenith angle                      c
!                                                                      c
! ng       = the number of gaussian angles to compute the diurnally    c
!            averaged solar radiation (=1 unless lswg = true)          c
!                                                                      c
! cloud    = flag for existence of a cloud (used only in 'ovc' mode)   c
!                                                                      c
!----------------------------------------------------------------------c
 
real, dimension(:,:,:),    intent(in)             :: gstr, cosang
logical, dimension(:,:),   intent(in)             :: daylight
real, dimension(:,:,:),    intent(inout)          :: taustr, omegastr
integer,                   intent(in)             :: ng
logical, dimension(:,:,:), intent(in), optional   :: cloud          

!----------------------------------------------------------------------c
! intent out:                                                          c
!                                                                      c
! rlayerdir = the layer reflectivity to a direct incident beam         c
!                                                                      c
! tlayerdir = the layer transmissivity to a direct incident beam       c
!                                                                      c
! rlayerdif = the layer reflectivity to a diffuse incident beam        c
!                                                                      c
! tlayerdif = the layer transmissivity to a diffuse incident beam      c
!                                                                      c
! tlayerde  = the layer transmissivity (non-scattered) to the direct   c
!             incident beam                                            c
!----------------------------------------------------------------------c
 
real, dimension(:,:,:), intent(out)    :: rlayerdir, rlayerdif,   &
                                          tlayerdir, tlayerdif,   &
	   				  tlayerde
 
!----------------------------------------------------------------------c
! local variables:                                                     c
!----------------------------------------------------------------------c

    real                     :: onedi3 = 1.0/3.0           
    real                     :: twodi3 = 2.0/3.0             
    integer                  :: k, i, ns, j, nn, ntot
    real  :: rsum, tsum, alpha
     real :: qq(7), rr(5), ss(8), tt(8),        ww(4)
 
       do k=KSRAD,KERAD
         do j=JSRAD,JERAD
!----------------------------------------------------------------------c
! overcast sky mode                                                    c
!                                                                      c
! note: in this mode, the delta-eddington method is performed only for c
!       spatial points containing a cloud.                             c
!----------------------------------------------------------------------c

     nn = 0
     if (present(cloud)) then
           do i=ISRAD,IERAD
             if (cloud(i,j,k) ) then
               nn = nn + 1
	       cld_index(i,j,k) = nn
               gstr2(nn) = gstr(i,j,k)
               taustr2(nn) = taustr(i,j,k)
               omegastr2(nn) = omegastr(i,j,k)
               cosangzk2(nn) = cosang(i,j,ng)
!----------------------------------------------------------------------c
! note: the following are done to avoid the conservative scattering    c
!       case, and to eliminate floating point errors in the            c
!       exponential calculations, respectively.                        c
!----------------------------------------------------------------------c
      if ( omegastr2(nn   ).ge.1.0 ) omegastr2(nn   ) = 9.9999999E-01
      if ( taustr2(nn   ).ge.1.0E+02 ) taustr2(nn   ) = 1.0E+02
             endif
           end do

!----------------------------------------------------------------------c
! clear sky mode                                                       c
!                                                                      c
! note: in this mode, the delta-eddington method is performed for all  c
!       spatial points.                                                c
!----------------------------------------------------------------------c

     else
           do i=ISRAD,IERAD
            if (daylight(i,j) ) then
             nn = nn + 1
	     cld_index(i,j,k) = nn
             gstr2(nn) = gstr(i,j,k)
             taustr2(nn) = taustr(i,j,k)
             omegastr2(nn) = omegastr(i,j,k)
             cosangzk2(nn) = cosang(i,j,ng)
!----------------------------------------------------------------------c
! note: the following are done to avoid the conservative scattering    c
!       case, and to eliminate floating point errors in the            c
!       exponential calculations, respectively.                        c
!----------------------------------------------------------------------c
      if ( omegastr2(nn   ).ge.1.0 ) omegastr2(nn   ) = 9.9999999E-01
      if ( taustr2(nn   ).ge.1.0E+02 ) taustr2(nn   ) = 1.0E+02
	    endif
           end do
    endif

    ntot = nn


    do nn = 1,ntot      

 
!----------------------------------------------------------------------c
! direct quantities                                                    c
!----------------------------------------------------------------------c
 
      ww(1) = omegastr2(nn)
      ww(2) = gstr2(nn)
      ww(3) = taustr2(nn)
      ww(4) = cosangzk2(nn)

      qq(1)     = 3.0 * ( 1.0 - ww(1) )
      qq(2)         = 1.0 - ww(1) * ww(2)
      qq(3)     = qq(1)/qq(2)
      qq(4) = sqrt( qq(1) * qq(2) )
      qq(5) = sqrt (qq(3))
      qq(6) = 1.0 + twodi3 * qq(5)         
      qq(7) = 1.0 - twodi3 * qq(5)       

      rr(1) = 1./qq(6)
      rr(2) = qq(7)*rr(1)
      rr(3) = exp( -ww(3)          * qq(4) )
      rr(4) = 1.0/rr(3)
      rr(5) = 1.0/(qq(6) * rr(4) - qq(7) * rr(3) *   &
                      rr(2) )



      ss(1) = 0.75 * ww(1)            / ( 1.0 - ( qq(4) *  &
                  ww(4)               ) ** 2 )
      ss(2) = ss(1) * ww(4)              *( 1.0 + ww(2)       *&
                  qq(1) * onedi3 )
      ss(3) = ss(1) * ( 1.0 + ww(2)        * qq(1) * &
                  ww(4)               ** 2 )
      ss(4) = ss(2) - twodi3 * ss(3)     
      ss(5) = ss(2) + twodi3 * ss(3)     
       ss(6) = exp( -ww(3)          / ww(4) )
      ss(7) = ( ss(4) * ss(6) - ss(5) * rr(3) * &
                  rr(2)      ) * rr(5)
      ss(8) = ( ss(5) - qq(7) * ss(7) ) * rr(1)


!----------------------------------------------------------------------c
! diffuse quantities                                                   c
!                                                                      c
! notes: the number of streams for the diffuse beam is fixed at 4.     c
!                                                                      c
!        this calculation is done only for ng=1.                       c
!----------------------------------------------------------------------c
 
       if ( ng.eq.1 ) then   
      rsum = 0.0
      tsum = 0.0
      if (nstr4) then
      do ns = 1,NSTREAMS


          tt(1) = 0.75 * ww(1)            / ( 1.0 - ( qq(4) * &
                       cosangstr(ns) ) ** 2 )
          tt(2) = tt(1) * cosangstr(ns) * ( 1.0 +  &
                       ww(2)        * qq(1) * onedi3 )
          tt(3) = tt(1) * ( 1.0 + ww(2)        * qq(1)*&
                       cosangstr(ns) ** 2 )
          tt(4) = tt(2) - twodi3 * tt(3)
          tt(5) = tt(2) + twodi3 * tt(3)
          tt(6) = exp( -ww(3)          / cosangstr(ns) )
          tt(7) = ( tt(4) * tt(6) - tt(5) *  &
                    rr(3) * rr(2)   ) * rr(5)
          tt(8) = ( tt(5) - qq(7) * tt(7) )*rr(1)



          rsum = rsum + (qq(7)*tt(8) + qq(6)*tt(7) - tt(4))*  &
                   wtstr(ns)*cosangstr(ns)
          tsum = tsum + ( (rr(3)*qq(6)*tt(8) + qq(7)*rr(4)*tt(7)    &
                   -tt(5)*tt(6)) + tt(6))*  &
                       wtstr(ns)*cosangstr(ns)
            
        end do
        else 
      do ns = 1,NSTREAMS


          tt(1) = 0.75 * ww(1)            / ( 1.0 - ( qq(4) * &
                       cosangstr(ns) ) ** 2 )
          tt(2) = tt(1) * cosangstr(ns) * ( 1.0 +  &
                       ww(2)        * qq(1) * onedi3 )
          tt(3) = tt(1) * ( 1.0 + ww(2)        * qq(1)*&
                       cosangstr(ns) ** 2 )
          tt(4) = tt(2) - twodi3 * tt(3)
          tt(5) = tt(2) + twodi3 * tt(3)
          tt(6) = exp( -ww(3)          / cosangstr(ns) )
          tt(7) = ( tt(4) * tt(6) - tt(5) *  &
                    rr(3) * rr(2)   ) * rr(5)
          tt(8) = ( tt(5) - qq(7) * tt(7) )*rr(1)



          rsum = rsum + (qq(7)*tt(8) + qq(6)*tt(7) - tt(4))
          tsum = tsum + ( (rr(3)*qq(6)*tt(8) + qq(7)*rr(4)*tt(7)    &
                   -tt(5)*tt(6)) + tt(6))
            
        end do
         endif
           sumr(nn) = rsum
           sumt(nn) = tsum


       endif  !  ng == 1

      rlayerdir2(nn) = qq(7) * ss(8) + qq(6)*ss(7) - &
                       ss(4)
      tlayerdir2(nn) = ((rr(3) * qq(6) * ss(8) + &
                       qq(7) * rr(4) * ss(7) -  &
                       ss(5) * ss(6) ) + ss(6) )
      tlayerde2(nn) = ss(6)

        end do  ! ntot loop

!---------------------------------------------------------------------
!     return results in proper locations in (i,j,k) arrays
!---------------------------------------------------------------------

   if (present(cloud)) then
         do i=ISRAD,IERAD
           if (cloud(i,j,k) ) then
             rlayerdir(i,j,k) = rlayerdir2(cld_index(i,j,k))
             tlayerdir(i,j,k) = tlayerdir2(cld_index(i,j,k))
             tlayerde(i,j,k) = tlayerde2(cld_index(i,j,k))
             if (ng .eq. 1) then
               rlayerdif(i,j,k) = sumr(cld_index(i,j,k))
               tlayerdif(i,j,k) = sumt(cld_index(i,j,k))
             endif
           endif
         end do
   else
         do i=ISRAD,IERAD
	   if (daylight(i,j) ) then
             rlayerdir(i,j,k) = rlayerdir2(cld_index(i,j,k))
             tlayerdir(i,j,k) = tlayerdir2(cld_index(i,j,k))
             tlayerde(i,j,k) = tlayerde2(cld_index(i,j,k))
             if (ng .eq. 1) then
               rlayerdif(i,j,k) = sumr(cld_index(i,j,k))
               tlayerdif(i,j,k) = sumt(cld_index(i,j,k))
             endif
           endif
         end do
   endif


     end do
     end do
 
end  subroutine deledd




!#####################################################################

			 end module esfsw_driver_mod
