
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
!use rad_utilities_mod,    only:  shortwave_control_type, Sw_control

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
   character(len=128)  :: version =  '$Id: esfsw_driver.F90,v 1.2 2001/07/05 17:31:09 fms Exp $'
   character(len=128)  :: tag     =  '$Name: eugene $'



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

!------------------------------------------------------------------
!  variables used in subroutine deledd  
!------------------------------------------------------------------
integer, dimension(:,:,:), allocatable ::  cld_index
real,    dimension(:,:),   allocatable ::  q, r, t
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
      else if (nstreams == 1) then
        ptstr(1) = ptstr_1
        wtstr(1) = 1.0
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

      p0h2o = 1.0E+2*p0h2o  ! mb to mks
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
 
real, dimension(:,:,:), intent(out)            :: dfsw, fsw, hsw, ufsw
real, dimension(:,:,:), intent(out) , optional :: dfswclr, fswclr,   &
                                                  hswclr, ufswclr
 
!-----------------------------------------------------------------------
!     local variables:
!-----------------------------------------------------------------------
 
   logical, dimension(:,:,:), allocatable :: cloud

   real,    dimension(:,:,:), allocatable ::     &
         aeroasymfac,     aerosctopdep,   aeroextopdep,     &
         alphaco2,        alphaco2str,    alphao2,          &
         alphao2str,      cloudasymfac,   cloudextopdep,    &
         cloudsctopdep,   delpdig,        deltap,           &
                          densitymol,     extopdepclr,      &
	 extopdepovc,     fclr,           fovc,             &
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
         sumtrclr,        sumreclr

   real, dimension (:,:,:,:), allocatable ::                &
                          dfswband,       efftauo2,         &
	 efftauco2,       fswband,        hswband,          &
	 totco2,          totco2str,      toto2,            &
	 toto2str,        ufswband

   real, dimension (:,:,:,:), allocatable ::                &
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
	 sfcalb,          solarflux,      wtfac

   integer            :: j, i, k, ng, ncld, ncldlyrs, nhiclds,    &
			 nmiclds, nloclds, np, nband, nf, ns
   character(len=15)  :: swaer_form

!--------------------------------------------------------------------
!    define limits and dimensions 
!--------------------------------------------------------------------

      IJK   = (IERAD-ISRAD+1)*(JERAD-JSRAD+1)*(KERAD-KSRAD+1)

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
allocate ( deltap          (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
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
! define a cloud existence variable.      
!----------------------------------------------------------------------c

	  do k=KSRAD,KERAD
            do j = JSRAD,JERAD
              do i = ISRAD,IERAD
	        if (cloudfrac(i,j,k) > 0.0)  then
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
          do k = KSRAD+1,KERAD+1
            totco2(:,:,k,ng) = 1.0E+02*rrvco2*scale(:,:,k)/(grav_mks* &
                             rhoair * cosangsolar(:,:,ng) * 2.0 )
            totco2str(:,:,k,ng) = 1.0E+02 * rrvco2 * scalestr(:,:,k) / &
                              (grav_mks * rhoair * cosangsolar(:,:,ng) )
            toto2(:,:,k,ng) = 1.0E+02 * o2mixrat * scale(:,:,k) /  &
                          (grav_mks * rhoair * cosangsolar(:,:,ng)*2.0 )
            toto2str(:,:,k,ng) = 1.0E+02*o2mixrat * scalestr(:,:,k) /  &
                             (grav_mks * rhoair * cosangsolar(:,:,ng) )
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
  allocate ( dfswbandclr(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1, NBANDS))
  allocate (  fswbandclr(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1, NBANDS)) 
  allocate (  hswbandclr(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  , NBANDS)) 
  allocate ( ufswbandclr(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1, NBANDS))   
endif

  allocate ( dfswband   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1, NBANDS))
  allocate (  fswband   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1, NBANDS))
  allocate (  hswband   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  , NBANDS))
  allocate ( ufswband   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1, NBANDS))
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
  allocate ( q              (IJK, 15         ))
  allocate ( r              (IJK, NSTREAMS   ))
  allocate ( t              (IJK, NSTREAMS   ))
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
            aeroextopdep(:,:,k) = aerext(:,:,k)*aeramtsc(:,:,k)*&
                                  deltaz(:,:,k) / 1.0E+03
            aerosctopdep(:,:,k) = aerssalb(:,:,k) * &
                                  aeroextopdep(:,:,k)
            aeroasymfac(:,:,k) = aerasymm(:,:,k)
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
		cloudextopdep(i,j,k) = cldext(i,j,k) *    &
                                   deltaz(i,j,k) / 1.0E+03
                cloudsctopdep(i,j,k) = cldsct(i,j,k) *    &
                                   deltaz(i,j,k) / 1.0E+03
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
                            exp( powph2o(nband)*alog(press(:,:,k)/ &
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
		      if (fracday(i,j) /= 0.0) then
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
		      if (fracday(i,j) /= 0.0) then
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
		      if (fracday(i,j) /= 0.0) then
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
            sctopdepclr(:,:,:) = rayopdep(:,:,:) + aerosctopdep(:,:,:)
            gclr(:,:,:) = aeroasymfac(:,:,:)*aerosctopdep(:,:,:) /  &
                          sctopdepclr(:,:,:)
            fclr(:,:,:) = aeroasymfac(:,:,:) ** 2*aerosctopdep(:,:,:)/&
                          sctopdepclr(:,:,:)
            gstrclr(:,:,:) = ( gclr(:,:,:)  - fclr(:,:,:) ) /   &
                             ( 1.0 - fclr(:,:,:) )
          endif

 !-----------------------------------------------------------------
 !    define cloudy sky arrays
 !-----------------------------------------------------------------

            do k=KSRAD,KERAD
              do j=JSRAD,JERAD
                do i=ISRAD,IERAD
                  if (cloud(i,j,k) ) then
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
                extopdepclr(:,:,:) = gasopdep(:,:,:) +    &
				     rayopdep(:,:,:) +    &
				     aeroextopdep(:,:,:)
                ssalbclr(:,:,:) = sctopdepclr(:,:,:)/extopdepclr(:,:,:)
 
 
                taustrclr(:,:,:) = extopdepclr(:,:,:) * ( 1.0 -    &
                                   ssalbclr(:,:,:) * fclr(:,:,:) )
                omegastrclr(:,:,:) = ssalbclr(:,:,:)*((1.0 -     &
				     fclr(:,:,:))/(1.0 -   &
				     ssalbclr(:,:,:)*fclr(:,:,:)))
!----------------------------------------------------------------------c
! calculate the reflection and transmission in the scattering layers   c
! using the delta-eddington method.                                    c
!----------------------------------------------------------------------c
                call deledd (taustrclr, omegastrclr, gstrclr,   &
			     cosangsolar, ng , fracday,rlayerdirclr,   &
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
			     cosangsolar, ng, fracday, rlayerdirovc,   &
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
                    if ( .not. cloud(i,j,k) ) then
	              rlayerdir(i,j,k) = rlayerdirclr(i,j,k)
	              tlayerdir(i,j,k) = tlayerdirclr(i,j,k)
	              rlayerdif(i,j,k) = rlayerdifclr(i,j,k)
	              tlayerdif(i,j,k) = tlayerdifclr(i,j,k)
	              tlayerde (i,j,k) = tlayerdeclr (i,j,k)
                    else
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
                           tlayerde , sfcalb   ,  &
                           reflectance, transmittance)    
 
              if (Rad_control%do_totcld_forcing) then

                call adding (                                    &
                             rlayerdirclr, tlayerdirclr,         &
                             rlayerdifclr, tlayerdifclr,         &
                             tlayerdeclr,  sfcalb,               &
                             reflectanceclr, transmittanceclr)

              endif

!----------------------------------------------------------------------c
! weight and sum the reflectance and transmittance to calculate the    c
! band values.                                                         c
!----------------------------------------------------------------------c
 
              wtfac(:,:) = wtfreq(np)*gausswt(ng)*cosangsolar(:,:,ng)
 
              do k = KSRAD,KERAD+1
                sumtr(:,:,k) = sumtr(:,:,k) + transmittance(:,:,k)*   &
		               wtfac(:,:)
                sumre(:,:,k) = sumre(:,:,k) + reflectance(:,:,k)* &
			       wtfac(:,:)
              end do

              if (Rad_control%do_totcld_forcing) then

                do k = KSRAD,KERAD+1
                  sumtrclr(:,:,k) = sumtrclr(:,:,k) +         &
                                    transmittanceclr(:,:,k)*wtfac(:,:) 
                  sumreclr(:,:,k) = sumreclr(:,:,k) +         &
                                    reflectanceclr(:,:,k)*wtfac(:,:)
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
 
          do k = KSRAD,KERAD+1
            dfswband(:,:,k,nband) = sumtr(:,:,k) * solarflux(:,:) 
            ufswband(:,:,k,nband) = sumre(:,:,k) * solarflux(:,:)
          end do
 
          fswband(:,:,:,nband) = ufswband(:,:,:,nband) -   &
				 dfswband(:,:,:,nband)
 
          do k = KSRAD,KERAD
            hswband(:,:,k,nband) = (fswband(:,:,k+1,nband) -    &
                                    fswband(:,:,k,nband) )*radcon_mks/ &
				    deltap(:,:,k) 
          end do

	  if (Rad_control%do_totcld_forcing) then

            do k = KSRAD,KERAD+1
              dfswbandclr(:,:,k,nband) = sumtrclr(:,:,k)*solarflux(:,:)
              ufswbandclr(:,:,k,nband) = sumreclr(:,:,k)*solarflux(:,:)
            end do
  
            fswbandclr(:,:,:,nband) = ufswbandclr(:,:,:,nband) -   &
                                      dfswbandclr(:,:,:,nband)

            do k = KSRAD,KERAD
              hswbandclr(:,:,k,nband) = (fswbandclr(:,:,k+1,nband) -  &
                                         fswbandclr(:,:,k,nband) ) *  &
                                         radcon_mks / deltap(:,:,k)
            end do
          endif

!--------------------------------------------------------------------
 
        end do      ! end of band loop

!----------------------------------------------------------------------c
! sum the band fluxes and heating rates to calculate the total         c
! spectral values.                                                     c
!----------------------------------------------------------------------c
 
        do nband = 1,NBANDS
          do k = KSRAD,KERAD+1
            dfsw(:,:,k) = dfsw(:,:,k) + dfswband(:,:,k,nband) 
            ufsw(:,:,k) = ufsw(:,:,k) + ufswband(:,:,k,nband)
            fsw(:,:,k) = fsw(:,:,k) + fswband(:,:,k,nband) 
          end do 
          do k = KSRAD,KERAD
            hsw(:,:,k) = hsw(:,:,k) + hswband(:,:,k,nband)
          end do 
        end do
      
        if (Rad_control%do_totcld_forcing) then
          do nband = 1,NBANDS
            do k = KSRAD,KERAD+1
              dfswclr(:,:,k) = dfswclr(:,:,k) + dfswbandclr(:,:,k,nband)
	      ufswclr(:,:,k) = ufswclr(:,:,k) + ufswbandclr(:,:,k,nband)
	      fswclr(:,:,k) = fswclr(:,:,k) + fswbandclr(:,:,k,nband)
	    end do
	    do k = KSRAD,KERAD
	      hswclr(:,:,k) = hswclr(:,:,k) + hswbandclr(:,:,k,nband)
	    end do
	  end do
	endif

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
      deallocate  (  t         )
      deallocate  (  r         )
      deallocate  (  q         )
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
      deallocate ( densitymol     )
      deallocate ( cloud          )

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
                   tlayerde, sfcalb, reflectance, transmittance)
 
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
!----------------------------------------------------------------------c
 
real, dimension(:,:,:),  intent(in)  :: rlayerdir, rlayerdif, &
                                        tlayerdir, tlayerdif, & 
                                        tlayerde
real, dimension (:,:) ,  intent(in)  :: sfcalb
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

!----------------------------------------------------------------------c
! initialization for the surface layer.                                c
!----------------------------------------------------------------------c
 
      raddupdif(:,:,KERAD+1) = sfcalb(:,:)
      raddupdir(:,:,KERAD+1) = sfcalb(:,:)
 
!----------------------------------------------------------------------c
! initialization for the first scattering layer.                       c
!----------------------------------------------------------------------c
 
      tadddowndir(:,:,KSRAD) = tlayerdir(:,:,KSRAD)
      radddowndif(:,:,KSRAD) = rlayerdif(:,:,KSRAD)
 
!----------------------------------------------------------------------c
! define the direct transmittance.                                     c
!----------------------------------------------------------------------c
 
      tlevel(:,:,KSRAD) = 1.0
      tlevel(:,:,KSRAD+1) = tlayerde(:,:,KSRAD)
 
      do k = KSRAD+1,KERAD
        tlevel(:,:,k+1) = tlevel(:,:,k) * tlayerde(:,:,k) 
      end do
 
!----------------------------------------------------------------------c
! add the inhomogeneous layers downward from the second layer to the   c
! surface.                                                             c
!                                                                      c
! radiation incident from below for diffuse beam, transmission of      c
! direct beam and conversion to diffuse.                               c
!----------------------------------------------------------------------c
 
      do k = KSRAD+1,KERAD
        dm1tl(:,:,k) = tlayerdif(:,:,k) / ( 1.0 - rlayerdif(:,:,k) *  &
                     radddowndif(:,:,k-1) )
        radddowndif(:,:,k) = rlayerdif(:,:,k) + radddowndif(:,:,k-1) * &
                           tlayerdif(:,:,k) * dm1tl(:,:,k)
        tadddowndir(:,:,k) = tlevel(:,:,k) * ( tlayerdir(:,:,k) + &
                           rlayerdir(:,:,k) * radddowndif(:,:,k-1) * &
                           dm1tl(:,:,k) ) + ( tadddowndir(:,:,k-1) -  &
                           tlevel(:,:,k) ) * dm1tl(:,:,k)              
      end do
 
!----------------------------------------------------------------------c
! add the inhomogeneous layers upward from the surface to the top of   c
! the atmosphere.                                                      c
!                                                                      c
! radiation incident from above for diffuse beam, reflection of        c
! direct beam and conversion to diffuse.                               c
!----------------------------------------------------------------------c
 
      do k = KERAD,KSRAD,-1
        dm2tl(:,:,k) = tlayerdif(:,:,k) / ( 1.0 - rlayerdif(:,:,k) *  &
                     raddupdif(:,:,k+1) )
        rdm2tl(:,:,k) = dm2tl(:,:,k) * raddupdif(:,:,k+1)
        raddupdif(:,:,k) = rlayerdif(:,:,k) + tlayerdif(:,:,k)*   &
			   rdm2tl(:,:,k) 
      end do
 
      do k = KERAD,KSRAD,-1
        raddupdir(:,:,k) = rlayerdir(:,:,k) + tlayerde(:,:,k) *   &
                           raddupdir(:,:,k+1) * dm2tl(:,:,k) +   &     
                           (tlayerdir(:,:,k) - tlayerde(:,:,k)) *  &
			   rdm2tl(:,:,k)
      end do
 

!----------------------------------------------------------------------c
! add downward to calculate the resultant reflectances and             c
! transmittances at flux levels.                                       c
!----------------------------------------------------------------------c
 
      do k = KSRAD+1,KERAD+1
        dm3(:,:,k) = 1.0/(1.0 - raddupdif(:,:,k)*radddowndif(:,:,k-1))
        dm3r(:,:,k) = dm3(:,:,k) * radddowndif(:,:,k-1)
        dm3r1p(:,:,k) = 1.0 + raddupdif(:,:,k) * dm3r(:,:,k)
        alpp(:,:,k) = (tadddowndir(:,:,k-1) - tlevel(:,:,k) )*dm3(:,:,k)
        transmittance(:,:,k) = (tlevel(:,:,k)*(1.0 + raddupdir(:,:,k)* &
                                dm3r(:,:,k)) + alpp(:,:,k))
        reflectance(:,:,k) = ( tlevel(:,:,k) * raddupdir(:,:,k) *   &
                             dm3r1p(:,:,k) + alpp(:,:,k) *   &
                             raddupdif(:,:,k) )
      end do
 

      reflectance(:,:,KSRAD) = raddupdir(:,:,KSRAD)
      transmittance(:,:,KSRAD) = 1.0
!------------------------------------------------------------------

end subroutine adding 



!####################################################################

subroutine deledd (taustr, omegastr, gstr, cosang, ng , fracday,  &
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
real, dimension(:,:),      intent(in)             :: fracday
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
 
!----------------------------------------------------------------------c
! overcast sky mode                                                    c
!                                                                      c
! note: in this mode, the delta-eddington method is performed only for c
!       spatial points containing a cloud.                             c
!----------------------------------------------------------------------c

     if (present(cloud)) then
       nn = 0
       do k=KSRAD,KERAD
         do j=JSRAD,JERAD
           do i=ISRAD,IERAD
             if (cloud(i,j,k) .and. fracday(i,j) /= 0.0) then
               nn = nn + 1
	       cld_index(i,j,k) = nn
               gstr2(nn) = gstr(i,j,k)
               taustr2(nn) = taustr(i,j,k)
               omegastr2(nn) = omegastr(i,j,k)
               cosangzk2(nn) = cosang(i,j,ng)
             endif
           end do
         end do
       end do

!----------------------------------------------------------------------c
! clear sky mode                                                       c
!                                                                      c
! note: in this mode, the delta-eddington method is performed for all  c
!       spatial points.                                                c
!----------------------------------------------------------------------c

     else
       nn = 0
       do k=KSRAD,KERAD
         do j=JSRAD,JERAD
           do i=ISRAD,IERAD
            if (fracday(i,j) /= 0.0) then
             nn = nn + 1
	     cld_index(i,j,k) = nn
             gstr2(nn) = gstr(i,j,k)
             taustr2(nn) = taustr(i,j,k)
             omegastr2(nn) = omegastr(i,j,k)
             cosangzk2(nn) = cosang(i,j,ng)
	    endif
           end do
         end do
       end do
    endif

    ntot = nn

!----------------------------------------------------------------------c
! note: the following are done to avoid the conservative scattering    c
!       case, and to eliminate floating point errors in the            c
!       exponential calculations, respectively.                        c
!----------------------------------------------------------------------c

    do nn = 1,ntot      
      if ( omegastr2(nn   ).ge.1.0 ) omegastr2(nn   ) = 9.9999999E-01
      if ( taustr2(nn   ).ge.1.0E+02 ) taustr2(nn   ) = 1.0E+02
    end do
 
!----------------------------------------------------------------------c
! direct quantities                                                    c
!----------------------------------------------------------------------c
 
    do nn = 1,ntot      
      q (nn,1 ) = 3.0 * ( 1.0 - omegastr2(nn   ) )
      q (nn,2 ) = 1.0 - omegastr2(nn   ) * gstr2(nn   )
      q (nn,3 ) = sqrt( q (nn,1 ) * q (nn,2 ) )
      q (nn,4 ) = 1.0 + twodi3 * sqrt( q (nn,1 ) / q (nn,2 ) ) 
      q (nn,5 ) = 1.0 - twodi3 * sqrt( q (nn,1 ) / q (nn,2 ) )
      q (nn,6 ) = exp( -taustr2(nn   ) * q (nn,3 ) )
      q (nn,7 ) = q (nn,4 ) / q (nn,6 ) - q (nn,5 ) * q (nn,6 ) *   &
                  q (nn,5 ) / q (nn,4 )
      q (nn,8 ) = 0.75 * omegastr2(nn   ) / ( 1.0 - ( q (nn,3 ) *  &
                  cosangzk2(nn      ) ) ** 2 )
      q (nn,9 ) = q (nn,8 ) * cosangzk2(nn      )*( 1.0 + gstr2(nn   )*&
                  q (nn,1 ) * onedi3 )
      q (nn,10) = q (nn ,8) * ( 1.0 + gstr2(nn   ) * q (nn,1 ) * &
                  cosangzk2(nn      ) ** 2 )
      q (nn,11) = q (nn,9 ) - twodi3 * q  (nn,10)     
      q (nn,12) = q (nn,9 ) + twodi3 * q  (nn,10)     
      q (nn,13) = exp( -taustr2(nn   ) / cosangzk2(nn      ) )
      q (nn,14) = ( q  (nn,11) * q  (nn,13) - q  (nn,12) * q (nn,6 ) / &
                  q (nn,4 ) * q (nn,5 ) ) / q (nn,7 )
      q (nn,15) = ( q  (nn,12) - q (nn,5 ) * q  (nn,14) ) / q (nn,4 )
      rlayerdir2(nn) = q (nn,5 ) * q  (nn,15) + q (nn,4 )*q  (nn,14) - &
                       q  (nn,11)
      tlayerdir2(nn) = ((q (nn,6 ) * q (nn,4 ) * q  (nn,15) + &
                       q (nn,5 ) / q (nn,6 ) * q  (nn,14) -  &
                       q  (nn,12) * q  (nn,13) ) + q  (nn,13) )
      tlayerde2(nn) = q (nn,13)
    end do

!----------------------------------------------------------------------c
! diffuse quantities                                                   c
!                                                                      c
! notes: the number of streams for the diffuse beam is fixed at 4.     c
!                                                                      c
!        this calculation is done only for ng=1.                       c
!----------------------------------------------------------------------c
 
    if ( ng.eq.1 ) then   
      do ns = 1,NSTREAMS
	do nn=1,ntot
          q  (nn,8 ) = 0.75 * omegastr2(nn   ) / ( 1.0 - ( q (nn,3 ) * &
                       cosangstr(ns) ) ** 2 )
          q  (nn,9 ) = q  (nn,8 ) * cosangstr(ns) * ( 1.0 +  &
                       gstr2(nn   ) * q (nn,1 ) * onedi3 )
          q  (nn,10) = q  (nn,8 ) * ( 1.0 + gstr2(nn   ) * q (nn,1 )*&
                       cosangstr(ns) ** 2 )
          q  (nn,11) = q  (nn,9 ) - twodi3 * q   (nn,10)
          q  (nn,12) = q  (nn,9 ) + twodi3 * q   (nn,10)
          q  (nn,13) = exp( -taustr2(nn   ) / cosangstr(ns) )
          q  (nn,14) = ( q   (nn,11) * q   (nn,13) - q   (nn,12) *  &
                       q (nn,6 ) / q (nn,4 ) * q (nn,5 ) ) / q (nn,7 )
          q  (nn,15) = ( q (nn,12) - q (nn,5 ) * q (nn,14) ) / q (nn,4 )
          r  (nn,ns) = q (nn,5 ) * q(nn,15) + q(nn,4 ) * q (nn,14) -&
                       q   (nn,11)
          t  (nn,ns) = (q (nn,6 ) * q (nn,4 ) * q(nn,15) + q (nn,5 ) /&
                       q (nn,6 ) * q   (nn,14) - q   (nn,12) * &
                       q   (nn,13) ) + q   (nn,13)
        end do
      end do

      sumr(:) = 0.0
      sumt(:) = 0.0
 
      if (NSTREAMS == 1) then
        do nn=1,ntot
          sumr(nn) = r(nn,1)
          sumt(nn) = t(nn,1)
        end do
      else
        do ns = 1,NSTREAMS
          do nn=1,ntot
            sumr(nn) = sumr(nn) + r(nn,ns) * wtstr(ns) * cosangstr(ns)
            sumt(nn) = sumt(nn) + t(nn,ns) * wtstr(ns) * cosangstr(ns)
          end do
        end do
      endif

   endif

!---------------------------------------------------------------------
!     return results in proper locations in (i,j,k) arrays
!---------------------------------------------------------------------

   if (present(cloud)) then
     nn = 0
     do k=KSRAD,KERAD
       do j=JSRAD,JERAD
         do i=ISRAD,IERAD
           if (cloud(i,j,k) .and. fracday(i,j) /= 0.0 ) then
             nn = nn + 1
             rlayerdir(i,j,k) = rlayerdir2(cld_index(i,j,k))
             tlayerdir(i,j,k) = tlayerdir2(cld_index(i,j,k))
             tlayerde(i,j,k) = tlayerde2(cld_index(i,j,k))
             if (ng .eq. 1) then
               rlayerdif(i,j,k) = sumr(cld_index(i,j,k))
               tlayerdif(i,j,k) = sumt(cld_index(i,j,k))
             endif
	   else
             rlayerdir(i,j,k) = 0.0
             tlayerdir(i,j,k) = 0.0
             tlayerde(i,j,k) = 0.0
             if (ng .eq. 1) then
               rlayerdif(i,j,k) = 0.0
               tlayerdif(i,j,k) = 0.0
             endif
           endif
         end do
       end do
     end do
   else
     nn = 0
     do k=KSRAD,KERAD
       do j=JSRAD,JERAD
         do i=ISRAD,IERAD
	   if (fracday(i,j) /= 0.0) then
             nn = nn + 1
             rlayerdir(i,j,k) = rlayerdir2(cld_index(i,j,k))
             tlayerdir(i,j,k) = tlayerdir2(cld_index(i,j,k))
             tlayerde(i,j,k) = tlayerde2(cld_index(i,j,k))
             if (ng .eq. 1) then
               rlayerdif(i,j,k) = sumr(cld_index(i,j,k))
               tlayerdif(i,j,k) = sumt(cld_index(i,j,k))
             endif
	   else
             rlayerdir(i,j,k) = 0.0
             tlayerdir(i,j,k) = 0.0
             tlayerde(i,j,k) = 0.0
             if (ng .eq. 1) then
               rlayerdif(i,j,k) = 0.0
               tlayerdif(i,j,k) = 0.0
             endif
           endif
         end do
       end do
     end do
   endif

!-------------------------------------------------------------------
   if  ( nn /= ntot) then
     call error_mesg ('deledd', &
          'final number of deledd points differs from original', FATAL)
   endif

!--------------------------------------------------------------------

 
end  subroutine deledd




!#####################################################################

			 end module esfsw_driver_mod
