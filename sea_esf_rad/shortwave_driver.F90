
                     module shortwave_driver_mod

use rad_output_file_mod,  only:  hold_sw, hold_sfcalbedo
use rad_step_setup_mod,   only:  temp, press, rh2o, jabs, &
                                 ISRAD,IERAD, JSRAD,JERAD, KSRAD, &
			         KERAD, albedo_sv, sealp
use rad_utilities_mod,    only:  Rad_control, radiation_control_type, &
				 Sw_control, shortwave_control_type, &
				 Environment, environment_type
use  utilities_mod,       only:  open_file, file_exist,      &
                                 check_nml_error, error_mesg, & 
                                 print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file
use esfsw_driver_mod,     only:  swresf, esfsw_driver_init
use esfsw_scattering_mod, only:  esfsw_scattering_init
use lhsw_driver_mod,      only:  swrad, lhsw_driver_init
use radiation_diag_mod,   only:  radiag_from_sw_driver_init, &
			         radiag_from_sw_driver
use ozone_mod,            only:  get_ozone
use cloudrad_package_mod, only:  get_clouds_for_esfsw
use astronomy_package_mod,only:  get_astronomy_for_swrad,   &
			         get_solar_distance,  &
			         get_astronomy_for_swrad_init
use surface_albedo_mod,   only:  get_surfacealbedo_for_swrad
use esfsw_parameters_mod, only:  esfsw_parameters_init



!-------------------------------------------------------------------

implicit none
private

!------------------------------------------------------------------
!              driver for shortwave radiation calculation
!
!-----------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.09'
    character(len=128)  :: version =  '$Id: shortwave_driver.F90,v 1.2 2001/07/05 17:33:21 fms Exp $'
    character(len=128)  :: tag     =  '$Name: fez $'



!---------------------------------------------------------------------
!-------  interfaces --------

public   shortwave_driver_init , shortwave, solar_constant, &
	 shortwave_driver_alloc, shortwave_driver_dealloc, &
	 get_gwt, get_sw_output, get_swflux_toa


private  gwtinit


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=10)   :: swform = '    '
integer            :: verbose= 0
  
  
 
namelist / shortwave_driver_nml /             &
                                     swform, verbose



!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


!------------------------------------------------------------------
!the following constants relate to the solar constant (watts/m**2):
!---------------------------------------------------------------------

! standard specification 

real,         parameter :: solar_mgroup =  1370.0   
real,         parameter :: solar_fms    =  1365.0   

!#ifdef oldsolarcnst (previous skyhi - conversion to and from ly/min)
! (1370. / 1.434e-3 * 697.667)
real,         parameter :: solar_oldsky =  1370.62263486 

!ifdef offlinesolarcnst (1.96 ly/min)
 real,       parameter :: solar_offline =  1367.42732      

!#ifdef rcemnsolarcnst (specify as 1370. / 697.66 * 697.67)
 real,         parameter :: solar_rcemn =  1370.0137459507           
!----------------------------------------------------------------------

real           :: solar
logical        :: do_lhsw  = .false.
logical        :: do_esfsw = .false.
logical        :: ldiurn, lswg, do_annual
integer        :: nsolwg
real           :: ssolar

!--------------------------------------------------------------------
!     fsw     =  net shortwave flux (up-down) at all model flux levels.
!     hsw     =  shortwave heating rates in model layers.
!     dfsw    =  downward short-wave radiation at pressure levels
!     fswcf, hswcf = values for corresponding variables (without ...cf)
!                    computed for cloud-free case
!--------------------------------------------------------------------
real, dimension(:,:,:), allocatable     :: dfsw, fsw, ufsw, hsw, &
                                           dfswcf, fswcf, ufswcf, &
				           hswcf



real, dimension (:), allocatable        :: gwt



!-------------------------------------------------------------------
!-------------------------------------------------------------------




                         contains



subroutine shortwave_driver_init (kmin, kmax)

integer, intent(in)    :: kmin, kmax

!---------------------------------------------------------------------
      integer   :: unit, io, ierr

!----------------------------------------------------------------
!-----  read namelist  ------
 
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=shortwave_driver_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'shortwave_driver_nml')
        enddo
10      call close_file (unit)
      endif
  
      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'shortwave_driver',   &
!						version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=shortwave_driver_nml)
      endif
      call close_file (unit)

      if (Environment%using_fms_periphs) then
	solar = solar_fms
      else if (Environment%using_sky_periphs) then
        solar = solar_mgroup
      endif

      Sw_control%sw_form = swform

      if (trim(swform) == 'lhsw') then
        do_lhsw  = .true.
        do_esfsw = .false.
      else if (trim(swform) == 'esfsw99') then
        do_lhsw  = .false.
        do_esfsw = .true.
      else
        call error_mesg ( 'shortwave_driver_init',   &
        'improper specification of desired shortwave parameterization',&
                                                               FATAL)
      endif


      call get_astronomy_for_swrad_init (ldiurn, lswg, nsolwg,   &
					 do_annual)

!---------------------------------------------------------------------
!     check for acceptable zenith angle specification.
!---------------------------------------------------------------------
      if (Environment%running_skyhi .and. do_annual) then
	call error_mesg ( 'shortwave_driver_init', &
	               ' cannot use do_annual in skyhi.',  FATAL)
      endif

!---------------------------------------------------------------  
!     define gaussian weights for gaussian calculations
!----------------------------------------------------------------

      call gwtinit

!---------------------------------------------------------------  
!     pass needed values to the diagnostics package    
!----------------------------------------------------------------
      call radiag_from_sw_driver_init (ldiurn,lswg, nsolwg, gwt, &
                                       do_lhsw, do_esfsw, do_annual)

!----------------------------------------------------------------------
      if (do_lhsw) then
	call lhsw_driver_init (kmin, kmax)
      else if (do_esfsw) then
	call esfsw_parameters_init
	call esfsw_driver_init
	call esfsw_scattering_init (kmin, kmax)
      endif
!-------------------------------------------------------------------



end subroutine shortwave_driver_init





!###########################################################

subroutine shortwave                                               

!-----------------------------------------------------------------------
!   shortwave handles the calculation of the shortwave radiation.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!
!-----------------------------------------------------------------------
!     intent in:
!     camtsw  =  shortwave cloud amounts. their locations are specified
!                in the ktopsw/kbtmsw indices. when no clouds are 
!                present camtsw should be given a value of zero.
!
!     cirabsw =  absorptivity of clouds in the infrared frequency band.
!                when no clouds are present cirabsw should be given a
!                value of zero. may be zenith angle dependent.
!
!     cirrfsw =  reflectivity of clouds in the infrared frequency band.
!                when no clouds are present cirrfsw should be given a
!                value of zero. may be zenith angle dependent.
!
!     cvisrfsw=  reflectivity of clouds in the visible frequency band.
!                when no clouds are present cvisrfsw should be given a
!                value of zero. may be zenith angle dependent.
!
!     kbtmsw  =  index of flux level pressure of cloud bottom.  when no
!                clouds are present kbtmsw should be given a value of
!                KS.
!
!     ktopsw  =  index of flux level pressure of cloud top.  when no
!                clouds are present ktopsw should be given a value of
!                KS.
!
!     ncldsw  =  number of clouds at each grid point.
!
!     press   =  pressure at data levels of model.
!
!     qo3     =  mass mixing ratio of o3 at model data levels.
!
!     rh2o    =  mass mixing ratio of h2o at model data levels.
!
!     ssolar  =  solar constant.
!     rsun    =  distance earth to sun.
!
!     temp    =  temperature at data levels of model.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
      real, dimension(:,:,:), allocatable  :: dfswg, fswg, hswg,    &
					      ufswg, dfswgcf, fswgcf,  &
					      hswgcf, ufswgcf,  &
					      camtswcf, press_mks
      integer, dimension(:,:), allocatable :: ncldswcf
      real, dimension(:,:,:),  allocatable :: cosangsolar
      real, dimension(:,:),    allocatable :: fracday              
      real, dimension(:,:),    allocatable :: cvisrfgd, cirabgd, cirrfgd
      real, dimension(:,:,:),  allocatable :: camtsw             
      real, dimension (:),     allocatable :: gausswt
      real, dimension(:,:,:),  allocatable :: qo3

      logical  :: sw_with_clouds=.true.
      logical  :: skipswrad, with_clouds
      integer  :: k,j, i, n
      real     :: rsun

!----------------------------------------------------------------------
!    initialize fluxes and heating rates
!--------------------------------------------------------------------
      dfsw(:,:,:) = 0.0
      ufsw(:,:,:) = 0.0
      fsw(:,:,:) = 0.0
      hsw(:,:,:) = 0.0

      if (Rad_control%do_totcld_forcing) then
        fswcf (:,:,:) = 0.0E+00
        dfswcf(:,:,:) = 0.0E+00
        ufswcf(:,:,:) = 0.0E+00
        hswcf (:,:,:) = 0.0E+00
      endif

!------------------------------------------------------------------
!  retrieve the zenith angles, daylight fraction and solar distance.
!------------------------------------------------------------------
      allocate (cosangsolar (ISRAD:IERAD, JSRAD:JERAD, nsolwg) )
      allocate (fracday     (ISRAD:IERAD, JSRAD:JERAD        ) )
      call get_astronomy_for_swrad(cosangsolar, fracday)

      if (Environment%using_fms_periphs) then
        call get_solar_distance(rsun)
        ssolar = solar*rsun
        if (get_my_pe() == 0 .and.jabs(JSRAD) == 1 .and.   &
						    verbose > 0) then
          write (*,90070) ssolar
        endif
      endif
90070 format (/ ' solar flux = ', f12.6, ' watts/m**2'/)

!--------------------------------------------------------------------
! determine when the no-sun case exists at all points within a chunk   
! and bypass the radiation calculations for that chunk. for the gaus-
! sian quadrature case, test the first angle , which is the smallest of 
! the nswg zenith angles.
!--------------------------------------------------------------------
      skipswrad = .true.

      do j=JSRAD,JERAD
        if ( cosangsolar(ISRAD,j,1)  >  0.0 ) skipswrad = .false.
        if (ldiurn) then
          do i = ISRAD+1,IERAD
            if (cosangsolar(i,j,1)  >  0.0 )  then
              skipswrad      = .false.
              exit
            endif
          end do
        endif
      end do

!--------------------------------------------------------------------
!  if the sun is shining nowhere in the chunk, send needed information
!  to various modules, clean up the allocated arrays in this routine, 
!  and return.
!--------------------------------------------------------------------

      if (skipswrad)  then
        allocate( camtsw  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
        call get_clouds_for_esfsw (camtsw)
        if (Rad_control%do_totcld_forcing) then
          call hold_sw  (fsw, ufsw, fswcf, ufswcf)
        else
          call hold_sw (fsw, ufsw)
        endif
	allocate ( cirrfgd (ISRAD:IERAD, JSRAD:JERAD) )
	allocate ( cvisrfgd(ISRAD:IERAD, JSRAD:JERAD) )
	cirrfgd(:,:) = 0.0
	cvisrfgd(:,:) = 0.0
        call hold_sfcalbedo (cirrfgd, cvisrfgd)
        deallocate (cirrfgd     )
        deallocate (cvisrfgd     )
        deallocate (camtsw   )
        deallocate (cosangsolar)
        deallocate (fracday    )
        return       
      endif

!---------------------------------------------------------------------
! allocate space for and then retrieve the ozone mixing ratio from the
! ozone module.
!--------------------------------------------------------------------
      allocate (qo3 (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
      call get_ozone (qo3)

!---------------------------------------------------------------------
! allocate space for and then retrieve the surface radiative properties
! from either the surface albedo module (standalone) or use the values
! input from the land_surface module now resident in rad_step_setup_mod
! (gcm).
!--------------------------------------------------------------------
      allocate ( cirabgd (ISRAD:IERAD, JSRAD:JERAD) )
      allocate ( cirrfgd (ISRAD:IERAD, JSRAD:JERAD) )
      allocate ( cvisrfgd(ISRAD:IERAD, JSRAD:JERAD) )

      if (Environment%running_gcm) then
        do j=JSRAD,JERAD
          cirrfgd(:,j) = albedo_sv(:,j)
          cvisrfgd(:,j) = albedo_sv(:,j)
          cirabgd(:,j) = 0.0
        end do
      else if (Environment%running_standalone) then
        call get_surfacealbedo_for_swrad (cirabgd, cirrfgd, cvisrfgd)
      endif

!--------------------------------------------------------------------
!  send the surface albedos to be saved to the output files
!-------------------------------------------------------------------
      call hold_sfcalbedo (cirrfgd, cvisrfgd)

!---------------------------------------------------------------------
! calculate shortwave radiative forcing and fluxes using the 
! exponential-sum-fit parameterization.
!---------------------------------------------------------------------
      if (do_esfsw) then

!---------------------------------------------------------------------
! allocate space for and initialize variables for the esfsw code.
!---------------------------------------------------------------------
        allocate( press_mks(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1)  )
        press_mks(:,:,:) = 1.0E-01*press(:,:,:)

        allocate (gausswt(nsolwg))

        allocate( camtsw  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
        call get_clouds_for_esfsw (camtsw)

        if (Rad_control%do_totcld_forcing) then
          allocate( camtswcf(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
          camtswcf(:,:,:) = 0.0
        endif

!----------------------------------------------------------------------
! note: the cosine of the solar zenith angle is set to 1 when it is 0  
!       to avoid a divide by zero in swresf. since the daylight 
!       fraction is 0 for these points, solar irradiance will still
!       properly be computed as zero.
!----------------------------------------------------------------------
        do n=1,nsolwg
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
              if (cosangsolar(i,j,n) == 0.0) then
	        cosangsolar(i,j,n) = 1.0
              endif
            end do
          end do
        end do

!-----------------------------------------------------------------------
! compute diurnally-varying shortwave radiation.
!-----------------------------------------------------------------------
        if (ldiurn) then
          gausswt(:) = 1.0
	  if (Rad_control%do_totcld_forcing) then
            call swresf (camtsw, cirrfgd, cvisrfgd, fracday,  &
                         press_mks, qo3, rh2o, ssolar, temp,  &
                         dfsw, fsw, hsw, ufsw,                &
                         cosangsolar, gausswt, nsolwg,        &
                         dfswclr=dfswcf, fswclr=fswcf,        &
                         hswclr=hswcf, ufswclr=ufswcf)
          else
            call swresf (camtsw, cirrfgd, cvisrfgd, fracday,  &
                         press_mks, qo3, rh2o, ssolar, temp,  &
                         dfsw, fsw, hsw, ufsw,                &
                         cosangsolar, gausswt, nsolwg)           
          endif

!-----------------------------------------------------------------------
!     compute shortwave radiation with gaussian integration.
!-----------------------------------------------------------------------
        else if (lswg) then
	  gausswt(:) = gwt(:) 
	  if (Rad_control%do_totcld_forcing) then
            call swresf (camtsw, cirrfgd, cvisrfgd, fracday,  &
                         press_mks, qo3, rh2o, ssolar, temp,  &
                         dfsw, fsw, hsw, ufsw,                &
                         cosangsolar, gausswt, nsolwg,        &
                         dfswclr=dfswcf, fswclr=fswcf,        &
                         hswclr=hswcf, ufswclr=ufswcf)
          else
            call swresf (camtsw, cirrfgd, cvisrfgd, fracday,  &
                         press_mks, qo3, rh2o, ssolar, temp,  &
                         dfsw, fsw, hsw, ufsw,                &
                         cosangsolar, gausswt, nsolwg)
          endif

!-----------------------------------------------------------------------
!     compute diurnally-averaged shortwave radiation.
!-----------------------------------------------------------------------
        else 
          gausswt(:) = 1.0
	  if (Rad_control%do_totcld_forcing) then
            call swresf (camtsw, cirrfgd, cvisrfgd, fracday,  &
                         press_mks, qo3, rh2o, ssolar, temp,  &
                         dfsw, fsw, hsw, ufsw,                &
                         cosangsolar, gausswt, nsolwg,        &
                         dfswclr=dfswcf, fswclr=fswcf,        &
                         hswclr=hswcf, ufswclr=ufswcf)
          else
            call swresf (camtsw, cirrfgd, cvisrfgd, fracday,  &
                         press_mks, qo3, rh2o, ssolar, temp,  &
                         dfsw, fsw, hsw, ufsw,                &
                         cosangsolar, gausswt, nsolwg)
          endif
        endif

!---------------------------------------------------------------------
! calculate shortwave radiative forcing and fluxes using the 
! lacis-hansen parameterization.
!---------------------------------------------------------------------
      else if (do_lhsw) then

!---------------------------------------------------------------------
! allocate space for and initialize variables for the lhsw code.
!---------------------------------------------------------------------
        if (lswg) then
          allocate( dfswg   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
          allocate( fswg    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
          allocate( hswg    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
          allocate( ufswg   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )
          if (Rad_control%do_totcld_forcing) then
            allocate( dfswgcf (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
            allocate( fswgcf  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
            allocate( hswgcf  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ))
            allocate( ufswgcf (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1))
          endif
        endif

        with_clouds = .true.
!-----------------------------------------------------------------------
!     compute shortwave radiation with gaussian integration.
!-----------------------------------------------------------------------
        if (lswg) then
          call swrad (cosangsolar, fracday, with_clouds, cirabgd,    &
	 	      cirrfgd, cvisrfgd, press, qo3, rh2o, ssolar,   &
		      nsolwg, dfsw, fsw, hsw, ufsw, gwt=gwt)      
          if (Rad_control%do_totcld_forcing) then
	    with_clouds = .false.
            call swrad (cosangsolar, fracday, with_clouds, cirabgd,  &
	     	        cirrfgd, cvisrfgd, press, qo3, rh2o, ssolar,  &
		        nsolwg, dfswcf, fswcf, hswcf, ufswcf, gwt=gwt)
          endif

!---------------------------------------------------------------------
!   integrate either diurnal averaged or diurnally varying case
!--------------------------------------------------------------------
        else  
          call swrad (cosangsolar, fracday, with_clouds, cirabgd,    &
		      cirrfgd, cvisrfgd, press, qo3, rh2o, ssolar,  &
		      nsolwg, dfsw, fsw, hsw, ufsw)
          if (Rad_control%do_totcld_forcing) then
            with_clouds = .false.
            call swrad (cosangsolar, fracday, with_clouds, cirabgd,  &
	     	        cirrfgd, cvisrfgd, press, qo3, rh2o, ssolar, &
		        nsolwg, dfswcf, fswcf, hswcf, ufswcf)
          endif
        endif
      endif  

!--------------------------------------------------------------------
!  for SKYHI the fluxes are converted to dynes per centimeter**2  
!        from watts per meter**2 using the factor 1.0E+03.       
!--------------------------------------------------------------------
      if (do_esfsw) then 
        dfsw(:,:,:) = dfsw(:,:,:) *1.0E+3
        ufsw(:,:,:) = ufsw(:,:,:) *1.0E+3
        fsw(:,:,:) =  fsw(:,:,:) *1.0E+3
        if (Rad_control%do_totcld_forcing) then
          dfswcf(:,:,:) = dfswcf(:,:,:) *1.0E+3
          ufswcf(:,:,:) = ufswcf(:,:,:) *1.0E+3
          fswcf(:,:,:) =  fswcf(:,:,:) *1.0E+3
        endif
      endif

!--------------------------------------------------------------------
!  deallocate arrays of this module
!-------------------------------------------------------------------
      if (do_lhsw) then
        if (lswg) then
          if (Rad_control%do_totcld_forcing) then
            deallocate ( ufswgcf  )
            deallocate ( hswgcf   )
            deallocate ( fswgcf   )
            deallocate ( dfswgcf  )
          endif
          deallocate ( ufswg  )
          deallocate ( hswg   )
          deallocate ( fswg   )
          deallocate ( dfswg  )
        endif
      else 
        if (Rad_control%do_totcld_forcing) then
          deallocate ( camtswcf  )
        endif
        deallocate (camtsw   )
        deallocate (gausswt )
        deallocate ( press_mks)
      endif
      deallocate (cvisrfgd )
      deallocate (cirrfgd  )
      deallocate (cirabgd  )
      deallocate (qo3)
      deallocate (fracday    )
      deallocate (cosangsolar)

!-------------------------------------------------------------------
!  send diagnostics to radiation_diag_mod
!-------------------------------------------------------------------
      if (Rad_control%do_diagnostics) then
        if (Rad_control%do_totcld_forcing) then
	  do j=JSRAD,JERAD
            call radiag_from_sw_driver (j,    &
				        ssolar, dfsw, fsw, hsw, &
				        ufsw, dfswcf, fswcf, hswcf,  &
				        ufswcf)
          end do
        else
	  do j=JSRAD,JERAD
            call radiag_from_sw_driver (j,     &
				        ssolar, dfsw, fsw, hsw, &
				        ufsw)
          end do
        endif
      endif

!--------------------------------------------------------------------
!  send the shortwve fluxes to be held for later use by the archive 
!  package
!--------------------------------------------------------------------
      if (Rad_control%do_totcld_forcing) then
        call hold_sw  (fsw, ufsw, fswcf, ufswcf)
      else
        call hold_sw (fsw, ufsw)
      endif
!--------------------------------------------------------------------


end subroutine shortwave




!###################################################################

subroutine solar_constant 

!--------------------------------------------------------------
      real      :: rsun
!-----------------------------------------------------------------

      call get_solar_distance (rsun)
      if (Environment%using_sky_periphs) then
        ssolar = solar/(rsun**2)
	if (verbose > 0) then
          write (*,90070) ssolar
        endif
      endif
!--------------------------------------------------------------------


90070 format (/ ' solar flux = ', f12.6, ' watts/m**2'/)


end subroutine solar_constant 


!###################################################################

subroutine shortwave_driver_alloc


   allocate (fsw    ( ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )
   allocate (dfsw   ( ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )
   allocate (ufsw   ( ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )
   allocate (hsw    ( ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )

!  initialize needed so that avg history files are valid

   fsw   (:,:,:) = 0.0
   dfsw  (:,:,:) = 0.0
   ufsw  (:,:,:) = 0.0
   hsw   (:,:,:) = 0.0

   if (Rad_control%do_totcld_forcing) then
     allocate (fswcf  ( ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )
     allocate (dfswcf ( ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )
     allocate (ufswcf ( ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )
     allocate (hswcf  ( ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )

!  initialize needed so that avg history files are valid

     fswcf (:,:,:) = 0.0
     dfswcf(:,:,:) = 0.0
     ufswcf(:,:,:) = 0.0
     hswcf (:,:,:) = 0.0
   endif
!--------------------------------------------------------------------


end  subroutine shortwave_driver_alloc




!####################################################################

subroutine shortwave_driver_dealloc

    if (Rad_control%do_totcld_forcing) then
      deallocate (hswcf )
      deallocate (ufswcf)
      deallocate (dfswcf)
      deallocate (fswcf )
    endif

      deallocate (hsw )
      deallocate (ufsw)
      deallocate (dfsw)
      deallocate (fsw )


end subroutine shortwave_driver_dealloc




!####################################################################

subroutine get_gwt (gwt_out)

real, dimension(:), intent(out)  :: gwt_out

     gwt_out(:) = gwt(:)


end subroutine get_gwt 

!###################################################################

subroutine get_sw_output (dfsw_out, ufsw_out, fsw_out, hsw_out,    &
			  dfswcf_out, ufswcf_out, fswcf_out, hswcf_out)

real, dimension(:,:,:), intent(out)            :: dfsw_out, ufsw_out,&
                                                  fsw_out, hsw_out
real, dimension(:,:,:), intent(out), optional  :: dfswcf_out, &
					          ufswcf_out, &
					          fswcf_out,  &
					          hswcf_out

!----------------------------------------------------------------------
!  here the "..._out" arrays have model dimensions while the right-side
!  arrays are dimensioned by the radiation limits (ISRAD, etc.). if
!  these differ, the proper mapping must be made here.
!----------------------------------------------------------------------
  
       dfsw_out(:,:,:) = dfsw(:,:,:)
       ufsw_out(:,:,:) = ufsw(:,:,:)
       fsw_out (:,:,:) = fsw (:,:,:)
       hsw_out (:,:,:) = hsw (:,:,:)
  
       if (Rad_control%do_totcld_forcing) then
         dfswcf_out(:,:,:) = dfswcf(:,:,:)
         ufswcf_out(:,:,:) = ufswcf(:,:,:)
         fswcf_out (:,:,:) = fswcf (:,:,:)
         hswcf_out (:,:,:) = hswcf (:,:,:)
       endif
  
end subroutine get_sw_output
 

!####################################################################
  
subroutine get_swflux_toa (dfsw_out, ufsw_out, fsw_out)
 
real, dimension(:,:), intent(out)        :: dfsw_out, ufsw_out, &
					    fsw_out

!----------------------------------------------------------------------
!  here the "..._out" arrays have model dimensions while the right-side
!  arrays are dimensioned by the radiation limits (ISRAD, etc.). if
!  these differ, the proper mapping must be made here.
!----------------------------------------------------------------------
 
     dfsw_out(:,:) = dfsw  (:,:,KSRAD)
     ufsw_out(:,:) = ufsw  (:,:,KSRAD)
     fsw_out (:,:) = fsw   (:,:,KSRAD)

  
end subroutine get_swflux_toa
 


!####################################################################

subroutine gwtinit

!-----------------------------------------------------------------------
!     gwtinit initializes data for gaussian quadrature.
!
!     author: c. h. goldberg
!
!     revised: 1/1/93
!-----------------------------------------------------------------------
!     intent in:
!     nswg = number of quadrature points.
!-----------------------------------------------------------------------
!     intent out:
!     gpt  = gaussian quadrature points
!     gwt  = gaussian weights
!-----------------------------------------------------------------------
!     intent local:
!     gpt1-gpt8  = gaussian quadrature points for one, two, four, and
!                  eight-point quadratures.
!     gwt1-gwt8  = gaussian weights for one, two, four, and eight-point
!                  quadratures.
!-----------------------------------------------------------------------

      real, dimension(1) ::   gwt1
      real, dimension(2) ::   gwt2
      real, dimension(4) ::   gwt4
      real, dimension(8) ::   gwt8

!-----------------------------------------------------------------------
!
!     define gaussian quadrature points and weights.
!
!-----------------------------------------------------------------------

      data gwt1 /1.0E+00/
       
      data gwt2 /0.5E+00, 0.5E+00/

      data gwt4 /0.1739274226E+00, 0.3260725774E+00, 0.3260725774E+00, &
                 0.1739274226E+00/

      data gwt8 /0.0506142681E+00, 0.1111905172E+00, 0.1568533229E+00, &
                 0.1813418917E+00, 0.1813418917E+00, 0.1568533229E+00, &
                 0.1111905172E+00, 0.0506142681E+00/

!----------------------------------------------------------------------
!     fill desired quadrature points and weights arrays. 
!----------------------------------------------------------------------
      allocate ( gwt (nsolwg)  )

      if (nsolwg == 1) then
	gwt(:) = gwt1(:)
      else if (nsolwg == 2) then
	gwt(:) = gwt2(:)
      else if (nsolwg == 4) then
	gwt(:) = gwt4(:)
      else if (nsolwg == 8) then
	gwt(:) = gwt8(:)
      endif
!------------------------------------------------------------------
	
end  subroutine gwtinit



!####################################################################


                end module shortwave_driver_mod

