
                 module rh_based_clouds_mod

use rh_clouds_mod,          only: get_global_clouds
use microphys_rad_mod,      only: microphys_rad_driver,   &
				  microphys_rad_init
use utilities_mod,          only: open_file, file_exist,   &
                                  check_nml_error, error_mesg,   &
                                  print_version_number, FATAL, NOTE, &
				  WARNING, get_my_pe, close_file, &
				  get_domain_decomp
use rad_step_setup_mod,     only: press, jabs, iabs,  &
			          IMINP, IMAXP, JMINP, JMAXP, &
                                  ISRAD, IERAD, JSRAD, JERAD, & 
                                  KSRAD, KERAD
use rad_utilities_mod,      only: Environment, environment_type, &
                                  shortwave_control_type, Sw_control, &
                                  longwave_control_type, Lw_control
use longwave_setup_mod,     only: longwave_parameter_type, &    
                                  Lw_parameters
use astronomy_package_mod,  only: get_astronomy_for_clouds,  &
			          get_astronomy_for_clouds_init
use constants_new_mod,      only: radians_to_degrees
use donner_deep_mod,        only: get_cemetf, inquire_donner_deep

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!	           module which defines cloud locations
!                     based on model relative humidity
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

! character(len=5), parameter  ::  version_number = 'v0.09'
  character(len=128)  :: version =  '$Id: rh_based_clouds.F90,v 1.2 2001/07/05 17:33:04 fms Exp $'
  character(len=128)  :: tag     =  '$Name: fez $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          rh_clouds_init, rh_clouds_calc


private         &
	  cldalb, albcld_lw, albcld_sw


!---------------------------------------------------------------------
!-------- namelist  ---------


character(len=5)             :: cldht_type_form       = '     '
character(len=4)             :: cirrus_cld_prop_form  = '    '

!    logical variables derived from namelist input

logical                      :: do_part_black_cirrus=.false.
logical                      :: do_full_black_cirrus=.false.

logical                      :: do_cldht60 = .false.
logical                      :: do_cldht93 = .false.
logical                      :: do_cldhtskyhi = .false.





namelist /rh_based_clouds_nml /     &
			       cldht_type_form, &
			       cirrus_cld_prop_form


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!--------------------------------------------------------------------
!     define radiative properties for low, middle and high clouds
!     indices 1, 2 and 3, respectively).
!     cldem     : infrared emissivity
!     crfvis    : visible band reflectivity
!     crfir     : near-ir band reflectivity
!     cabir     : near-ir band absorptivity
!     crz       : cloud fraction
!--------------------------------------------------------------------
 
integer, parameter             :: NOFCLDS_SP=3  
integer, parameter             :: LATOBS=19     

real, dimension(NOFCLDS_SP)    :: crfvis_m,   crfir_m,   cabir_m,  &
                                  crfvis_fms, crfir_fms, cabir_fms, &
				  crfvis,     crfir,     cabir,  &
                                  cldem
real                           :: crz


data cldem / 1.0E+00, 1.0E+00, 1.0E+00 /  
data crz   / 1.00 /                     

!---------------------------------------------------------------------
!  these are m group values :
!---------------------------------------------------------------------

data crfvis_m  / 0.69E+00, 0.48E+00, 0.21E+00 / 
data crfir_m   / 0.69E+00, 0.48E+00, 0.21E+00 / 
data cabir_m   / 0.30E+00, 0.30E+00, 0.04E+00 /

!---------------------------------------------------------------------
!  these are fms values :
!---------------------------------------------------------------------
  
data crfvis_fms/ 0.59E+00, 0.45E+00, 0.21E+00 /
data crfir_fms / 0.59E+00, 0.45E+00, 0.21E+00 /
data cabir_fms / 0.40E+00, 0.30E+00, 0.04E+00 /
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!     these arrays define the cloud reflectivities as a function of 
!     zenith angle and the radiation band (visible and infrared) for 
!     high, middle and low clouds.
!     NREFL_BDS = number of radiative bands over which reflectivities
!              are provided.
!     NANGS     = number of zenith angles at which reflectivity values
!              are given.
!--------------------------------------------------------------------
 
integer, parameter                    ::  NANGS=17
integer, parameter                    ::  NREFL_BDS=2
real, dimension(NANGS,NREFL_BDS)      ::  albch, albcm, albcl 
 

!---------------------------------------------------------------------
!     albedos for high clouds at zenith angles from 0-80 deg. at 5 deg.
!     intervals for 1) visible and 2) infrared radiation.
!     for original skyhi values use:
!    &             .04,.05,.05,.05,.06,.06,.07,.07,.08,.11,.13,.16,.21,
!    &             .28,.39,.48,.61,
!    &             .04,.05,.05,.05,.06,.06,.07,.07,.08,.10,.11,.14,.21,
!    &             .26,.35,.44,.55 /
!     and specify a zenith angle of 60.000001.
!---------------------------------------------------------------------

data albch /      &
                 .04,.05,.05,.05,.06,.06,.07,.07,.08,.11,.13,.16,.21,  &
                 .28,.39,.48,.61,                                     &
                 .04,.05,.05,.05,.06,.06,.07,.07,.08,.10,.11,.14,.19, &
                 .26,.35,.44,.55 /

!---------------------------------------------------------------------
!     albedos for middle clouds at zenith angles from 0-80 deg. at 5 deg
!     intervals for 1) visible and 2) infrared radiation.
!     for original skyhi values use:
!    &            .18,.18,.19,.20,.21,.23,.24,.26,.29,.33,.37,.42,.48,
!    &            .55,.64,.71,.79,
!    &            .14,.14,.15,.16,.17,.18,.18,.20,.23,.25,.29,.32,.48,
!    &            .43,.50,.55,.61 /
!     and specify a zenith angle of 60.000001.
!----------------------------------------------------------------------

data albcm /     &
               .18,.18,.19,.20,.21,.23,.24,.26,.29,.33,.37,.42,.47,  &
               .55,.64,.71,.79,                                      &
               .14,.14,.15,.16,.17,.18,.18,.20,.23,.25,.29,.32,.37, &
               .43,.50,.55,.61 /

!-----------------------------------------------------------------------
!     albedos for low clouds at zenith angles from 0-80 deg. at 5 deg
!     intervals for 1) visible and 2) infrared radiation.
!     for original skyhi values use:
!    &             .50,.50,.51,.51,.52,.53,.54,.56,.58,.62,.65,.67,.69,
!    &             .73,.78,.82,.86,
!    &             .42,.42,.43,.43,.44,.45,.46,.48,.50,.52,.55,.57,.69,
!    &             .63,.66,.70,.74 /
!-----------------------------------------------------------------------

data albcl /     &
                .50,.50,.51,.51,.52,.53,.54,.56,.58,.62,.65,.67,.69, &
                .73,.78,.82,.86,                                     &
                .42,.42,.43,.43,.44,.45,.46,.48,.50,.52,.55,.57,.59, &
                .63,.66,.70,.74 /


 
!-------------------------------------------------------------------
!     this array defines the zenith angle dependent albedo for each of 
!     the different cloud types (NOFCLDS_SP) for each of the radiative 
!     bands (NREFL_BDS). currently NREFL_BDS are the visible and the
!     infrared.
!-------------------------------------------------------------------
 
real, dimension (:,:), allocatable   ::  zza 
 
!--------------------------------------------------------------------
!     these variables define the boundaries (in sigma coordinates) 
!     between high and middle and middle and low clouds. 
!--------------------------------------------------------------------

real, dimension(:), allocatable   :: cldhm_abs, cldml_abs,   &
				     cldhm_abs_gl, cldml_abs_gl
real                              :: cldhp, cldhe, cldmp, cldme
!-----------------------------------------------------------------


real, dimension(:), allocatable   :: qlevel

!-------------------------------------------------------------------
!    NLWCLDB is the actual number of frequency bands for which lw
!    emissitivies are defined. 
!--------------------------------------------------------------------
integer               :: NLWCLDB 
integer               :: nsolwg
character(len=10)     :: swform
logical               :: do_lwcldemiss
integer               :: x(4), y(4), id, jd, jdf


!----------------------------------------------------------------------
!----------------------------------------------------------------------




	 contains 


subroutine rh_clouds_init (th, qlevel_in)


real, dimension(:), intent(in)    :: th, qlevel_in

!--------------------------------------------------------------------
     real                              :: alatn, fl
     integer                           :: j,li, jj
     integer                           :: unit, ierr, io
     real, dimension(:), allocatable   :: cldhm, cldml

!---------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=rh_based_clouds_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'rh_based_clouds_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'rh_based_clouds',    &
!					       version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=rh_based_clouds_nml)
      endif
      call close_file (unit)

      call get_domain_decomp (x, y)
      jdf = y(4) -y(3) + 1
      id  = x(2) -x(1) + 1
      jd  = y(2) -y(1) + 1
      call get_astronomy_for_clouds_init (nsolwg)

      allocate ( qlevel(KSRAD:KERAD) )
      qlevel = qlevel_in

      NLWCLDB = Lw_parameters%NLWCLDB
      swform = Sw_control%sw_form
      do_lwcldemiss = Lw_control%do_lwcldemiss

      if (trim(swform) /= 'esfsw99' ) then
        allocate ( zza(NOFCLDS_SP, NREFL_BDS) )
      endif

!--------------------------------------------------------------------
!  define the formulation to be used for defining the interfaces 
!  between high, middle and low clouds
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!  use cloud height (low, middle, high) specification as in 1960
!  definition in SUPERSOURCE model
!--------------------------------------------------------------------
      if (trim(cldht_type_form) == '60') then
        do_cldht60 = .true.

!--------------------------------------------------------------------
!  use cloud height (low, middle, high) specification as in 1993
!  definition in SUPERSOURCE model
!--------------------------------------------------------------------
      else if (trim(cldht_type_form) == '93') then
        do_cldht93 = .true.

!--------------------------------------------------------------------
!  use cloud height (low, middle, high) specification as in 1993
!  definition in SKYHI model
!--------------------------------------------------------------------
      else if (trim(cldht_type_form) == 'skyhi') then
        do_cldhtskyhi = .true.

!--------------------------------------------------------------------
! error condition
!--------------------------------------------------------------------
      else
        call error_mesg( 'cloudrad_package_init',  &
                ' cldht_type_form is not an acceptable value.', FATAL)
      endif

!--------------------------------------------------------------------
! define the "blackness" of the cirrus clouds. cirrus clouds are either 
! "part black" with emissivity of 0.6 or are "full black" with emis-
! sivity of 1.0
!--------------------------------------------------------------------
      if ( .not. do_lwcldemiss) then
        if (trim(cirrus_cld_prop_form) == 'part') then
          do_part_black_cirrus = .true.
        else if (trim(cirrus_cld_prop_form) == 'full') then
          do_full_black_cirrus = .true.
        else
          call error_mesg( 'cloudrad_package_init',  &
                ' cirrus_cld_prop_form is not an acceptable value.', & 
								FATAL)
        endif
      endif

!--------------------------------------------------------------------
!   allocate arrays for the sigmas of the low-mid and mid-high cloud 
!   interfaces
!---------------------------------------------------------------------
      allocate ( cldhm_abs(jdf) )
      allocate ( cldml_abs(jdf) )
      allocate ( cldhm    (LATOBS) )
      allocate ( cldml    (LATOBS) )

!---------------------------------------------------------------------
!    select the ir absorptivity, ir reflectivity and visible reflect-
!    ivity values to use. 
!---------------------------------------------------------------------
      if (Environment%using_fms_periphs) then
        crfvis(:) = crfvis_fms(:)
        crfir (:) = crfir_fms (:)
        cabir(:) = cabir_fms(:)
      else if (Environment%using_sky_periphs) then
        crfvis(:) = crfvis_m(:)
        crfir (:) = crfir_m (:)
        cabir(:) = cabir_m(:)
      endif
!-----------------------------------------------------------------------
!     define the cloud height boundaries that distinguish low, middle 
!     and high clouds. values are defined in terms of the standard 
!     sigma coordinate at full levels.
!---------------------------------------------------------------------
        if (do_cldhtskyhi) then

!---------------------------------------------------------------------
!      this block supplies the standard skyhi values for these param-
!      eters, to be used for checkout when needed.
!----------------------------------------------------------------------
          do jj=1,LATOBS
            cldhm(jj) = 0.52
            cldml(jj) = 0.70         
          end do
          cldml(1) = 0.78
          cldml(2) = 0.78
          cldml(18) = 0.78
          cldml(19) = 0.78

!-----------------------------------------------------------------------
!     this definition was determined by Wetherald and Manabe circa 1960.
!-----------------------------------------------------------------------
        else  if (do_cldht60) then
          do j=1,LATOBS
            cldhm(j) = 0.4    
            cldml(j) = 0.7      
          end do

!-----------------------------------------------------------------------
!     this definition was determined by Wetherald, Manabe, Spelman and
!     Stouffer in 1993. here the cloud boundaries vary with latitude.
!-----------------------------------------------------------------------
        else if (do_cldht93) then
          cldhp = 0.7E+00
          cldhe = 0.4E+00
          cldmp = 0.85E+00
          cldme = 0.7E+00
          do jj=1,LATOBS
            alatn = ABS(-90.0E+00 + (jj-1)*10.0E+00)
            cldhm(jj) = cldhp + (90.0E+00-alatn)*(cldhe-cldhp)/90.0E+00
            cldml(jj) = cldmp + (90.0E+00-alatn)*(cldme-cldmp)/90.0E+00
          end do
        endif
!---------------------------------------------------------------------
!    perform latitude "interpolation" to the (JD) model latitudes
!    in default case, no interpolation is actually done; the nearest
!    latitude available (using NINT function) is used.
!---------------------------------------------------------------------
	allocate ( cldhm_abs_gl (jd) )
	allocate ( cldml_abs_gl (jd) )
        do j=1,jd
	  if (Environment%using_sky_periphs) then
            fl = 9.0E+00 - 9.0E+00*(FLOAT(JD+1 -j - JD/2) -    &
                 0.5E+00)/FLOAT(JD/2)
            li = NINT(fl) + 1
            cldhm_abs_gl(j) = cldhm(li)
            cldml_abs_gl(j) = cldml(li)
	  endif
        end do
        do j=1,jdf
	  if (Environment%using_sky_periphs) then
	    cldhm_abs(j) =  cldhm_abs_gl(j+y(3)-1)
	    cldml_abs(j) =  cldml_abs_gl(j+y(3)-1)
	  else if (Environment%using_fms_periphs) then
            cldhm_abs(j) = cldhp + (90.0E+00-abs(th(j)*   &
			   radians_to_degrees))*(cldhe-cldhp)/90.0E+00
            cldml_abs(j) = cldmp + (90.0E+00-abs(th(j)*    &
			   radians_to_degrees))*(cldme-cldmp)/90.0E+00
	  endif
        end do
	deallocate (cldhm_abs_gl)
	deallocate (cldml_abs_gl)

!---------------------------------------------------------------------
!    if a cloud microphysics scheme is to be employed with the cloud
!    scheme, initialize the microphysics_rad module.
!--------------------------------------------------------------------
        if (trim(swform) == 'esfsw99' .or. do_lwcldemiss) then
          call microphys_rad_init (cldhm_abs, cldml_abs)
        endif

!---------------------------------------------------------------------
!    deallocate local variables
!---------------------------------------------------------------------
        deallocate (cldhm) 
        deallocate (cldml) 


end subroutine rh_clouds_init




!######################################################################

subroutine rh_clouds_calc (camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
			   nrndlw, emmxolw, emrndlw, cirabsw,   &
			   cvisrfsw, cirrfsw, cldext, cldsct, &
			   cldasymm)

real,    dimension(:,:,:),   intent(inout) :: camtsw, cmxolw, crndlw
integer, dimension(:,:),     intent(inout) :: ncldsw, nrndlw, nmxolw
real,    dimension(:,:,:,:), intent(inout) :: emmxolw, emrndlw
real,    dimension(:,:,:,:), intent(inout), optional ::       &
				              cirabsw, cirrfsw, &
				              cvisrfsw, cldext, &
					      cldsct, cldasymm

!--------------------------------------------------------------------
!   these arrays define the basic radiative properties of the clouds.
!
!     cmxolw  =  amounts of maximally overlapped longwave clouds in
!                layers from KSRAD to KERAD.
!     crndlw  =  amounts of randomly overlapped longwave clouds in
!                layers from KSRAD to KERAD.
!     nmxolw  =  number of maximally overlapped longwave clouds
!                in each grid column.
!     nrndlw  =  number of maximally overlapped longwave clouds
!                in each grid column.
!     emmxolw =  longwave cloud emissivity for maximally overlapped
!                clouds through layers from KSRAD to KERAD.
!     emrndlw =  longwave cloud emissivity for randomly overlapped
!                clouds through layers from KSRAD to KERAD.
!     camtsw  =  shortwave cloud amounts. the sum of the maximally
!                overlapped and randomly overlapped longwave
!                cloud amounts.
!     ncldsw  =  number of shortwave clouds in each grid column.
!                                                                  
!    cldext   =  the parameterization band values of the cloud      
!                extinction coefficient in kilometer**(-1)          
!                                                                  
!    cldsct   =  the parameterization band values of the cloud      
!                scattering coefficient in kilometer**(-1)          
!                                                                  
!    cldasymm =  the parameterization band values of the asymmetry  
!                factor                                             
!     cirabsw =  absorptivity of clouds in the infrared frequency band.
!                may be zenith angle dependent.
!     cirrfsw =  reflectivity of clouds in the infrared frequency band.
!                may be zenith angle dependent.
!    cvisrfsw =  reflectivity of clouds in the visible frequency band.
!                may be zenith angle dependent.
!---------------------------------------------------------------------
 
      real, dimension(:,:,:), allocatable    :: ccover, cosangsolar, &
						cemetf
      real, dimension(:), allocatable        :: qlsig
      integer, dimension(:,:,:), allocatable :: iflagglbl, ifcd,&
                                                ifcv
      logical, dimension(:,:,:), allocatable :: hi_cloud, mid_cloud, &
						low_cloud, cvflag
      integer                                :: k, j, i, ngp
      real                                   :: cld
      logical                                :: do_donner_deep

!---------------------------------------------------------------------

      allocate (cvflag(IMINP:IMAXP, JMINP:JMAXP, ksrad:kerad) )

!--------------------------------------------------------------------
!   if donner_deep is active, retrieve the array of heating rates
!   due to convection. wherever the heating rate from convection is
!   non-zero, set a flag which will be used to indicate to the radiat-
!   ion package that cloud is present in that grid box.
!---------------------------------------------------------------------
      call inquire_donner_deep (do_donner_deep)
      if (do_donner_deep) then
        allocate (cemetf(IMINP:IMAXP, JMINP:JMAXP, ksrad:kerad))
        call get_cemetf(jabs(jminp), cemetf)
        do k=ksrad,kerad
          do j=jminp,jmaxp
            do i=iminp,imaxp
              if (cemetf(i,j,k) /=  0.0)  then
                cvflag(i,j,k) = .true.
              else
                cvflag(i,j,k) = .false.
              endif
            end do
          end do
        end do
        deallocate (cemetf)
      else
        do k=ksrad,kerad
          do j=jminp,jmaxp
            do i=iminp,imaxp
              cvflag(i,j,k) = .false.
            end do
          end do
        end do
      endif
  


      allocate (cosangsolar (ISRAD:IERAD, JSRAD:JERAD, nsolwg) )
      allocate (ifcd   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
      allocate (ifcv   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
      allocate (ccover (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
      allocate (qlsig  (                          KSRAD:KERAD) )
      allocate (hi_cloud(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
      allocate (mid_cloud(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
      allocate (low_cloud(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )

      if (Environment%running_fms) then
        allocate (iflagglbl (IMINP:IMAXP, JMINP:JMAXP, KSRAD:KERAD) )
        iflagglbl(:,:,:) = 0
      else if (Environment%running_skyhi) then
        allocate (iflagglbl (id, jd, KSRAD:KERAD) )
      endif

!---------------------------------------------------------------------
!     obtain the appropriate zenith angles that are to be used here.
!---------------------------------------------------------------------
      call get_astronomy_for_clouds (cosangsolar)

!--------------------------------------------------------------------- 
!     for the predicted clouds based on the rh distribution, obtain the
!     cloud flag array
!-------------------------------------------------------------
      call get_global_clouds (iabs(IMINP), jabs(JMINP), press,  &
                              iflagglbl)
!-----------------------------------------------------------------------
!     locate those points with cloudiness. set a flag to indicate if
!     it is stable or convective cloudiness. do not allow clouds at
!     model level 1. 
!-----------------------------------------------------------------------
      if (Environment%running_fms) then
        do k=KSRAD+1,KERAD
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
              if (iflagglbl(i,j,k) .EQ. 3 .or. cvflag(i,j,k) ) then
                ifcv(i,j,k) = 1
                ifcd(i,j,k) = 0
                ccover(i,j,k) = crz
              else if (iflagglbl(i,j,k) .EQ. 1) then
                ifcv(i,j,k) = 0
                ifcd(i,j,k) = 1
                ccover(i,j,k) = crz
              else
                ifcv(i,j,k) = 0
                ifcd(i,j,k) = 0
                ccover(i,j,k) = 0.0E+00
              endif
            end do
          end do
        end do
      else if (Environment%running_skyhi) then
        do k=KSRAD+1,KERAD
          do j=JSRAD,JERAD
            do i=ISRAD,IERAD
              if (iflagglbl(iabs (i),jabs (j),k) .EQ. 3 .or.  &
		  cvflag(i,j,k) ) then
	        ifcv(i,j,k) = 1
	        ifcd(i,j,k) = 0
	        ccover(i,j,k) = crz 
	      else if (iflagglbl(iabs (i),jabs (j),k) .EQ. 1) then
	        ifcv(i,j,k) = 0
	        ifcd(i,j,k) = 1
	        ccover(i,j,k) = crz 
	      else
	        ifcv(i,j,k) = 0
	        ifcd(i,j,k) = 0
	        ccover(i,j,k) = 0.0E+00
	      endif
	    end do
	  end do
        end do
      endif

!---------------------------------------------------------------------
!     define the number of each type of cloud in each column and the 
!     amount of each cloud type present in each grid box. allow for 
!     thick clouds if cloud is present at two adjacent levels.
!---------------------------------------------------------------------
      do j=JSRAD,JERAD
        do i=ISRAD,IERAD
          cld = 0.
          do k=KERAD,KSRAD+1,-1
            if (cld .EQ. 0.0E+00) then
              if (ccover(i,j,k) .NE. 0.0E+00) then
                if (ccover(i,j,k-1) .EQ. 0.0E+00 ) then
                  nrndlw(i,j) = nrndlw(i,j) + 1
                  crndlw(i,j,k) = ccover(i,j,k)
                else
                  cld = 1.0E+00
                  cmxolw(i,j,k) = ccover(i,j,k)
	        endif
              endif
            else
              if (ccover(i,j,k) .NE. 0.0E+00) then
                cmxolw(i,j,k) = ccover(i,j,k)
                if (k .EQ. KSRAD+1) then
                  nmxolw(i,j) = nmxolw(i,j) + 1
                endif
                if (ccover(i,j,k-1) .EQ. 0.0E+00) then
                  nmxolw(i,j) = nmxolw(i,j) + 1
                  cld = 0.0E+00
                endif
              endif
            endif
            camtsw(i,j,k) = cmxolw(i,j,k) + crndlw(i,j,k)

            if (camtsw(i,j,k) > 0.0) then

!-------------------------------------------------------------------
!      if cloud is present, define it to be either high, middle or low.
!      define a sigma variable to be used for this purpose.
!-------------------------------------------------------------------
              if (Environment%using_fms_periphs) then
                qlsig(k) = press(i,j,k)/press(i,j,KERAD+1)
              else if (Environment%using_sky_periphs) then
                qlsig(k) = qlevel(k)
              endif

	      if (qlsig(k) .LE. cldhm_abs(jabs(j)) ) then
	        hi_cloud(i,j,k) = .true.
	        mid_cloud(i,j,k) = .false.
	        low_cloud(i,j,k) = .false.
	      else if (qlsig(k) .GT. cldhm_abs(jabs(j))  .and.  &
	               qlsig(k) .LT. cldml_abs(jabs(j)) ) then
	        hi_cloud(i,j,k) = .false.
	        mid_cloud(i,j,k) = .true.
	        low_cloud(i,j,k) = .false.
	      else if (qlsig(k) .GE. cldml_abs(jabs(j)) ) then
	        hi_cloud(i,j,k) = .false.
	        mid_cloud(i,j,k) = .false.
	        low_cloud(i,j,k) = .true.
              else
	        call error_mesg ('rh_clouds_calc',  &
                     'model level is not mapped to a cloud type', FATAL)
              endif
!-------------------------------------------------------------------
!      if no cloud is present, define all cloud types as false.
!-------------------------------------------------------------------
	    else
	      hi_cloud(i,j,k) = .false.
	      mid_cloud(i,j,k) = .false.
	      low_cloud(i,j,k) = .false.
            endif
          end do
        end do
      end do

!-------------------------------------------------------------------
!      define the number of clouds present in each column.
!-------------------------------------------------------------------
      do j=JSRAD,JERAD
        do i=ISRAD,IERAD
          ncldsw(i,j) = nmxolw(i,j) + nrndlw(i,j)
        end do
      end do

!-------------------------------------------------------------------
!     if microphysics not being used for lw calculation, call albcld_lw
!     to define long-wave emissivities for the model clouds.
!-------------------------------------------------------------------
      if (.not. do_lwcldemiss) then
        call albcld_lw (hi_cloud, mid_cloud,  &
	 	       low_cloud, cmxolw, crndlw,  emmxolw,  emrndlw)
      endif

!-------------------------------------------------------------------
!     if microphysics not being used for sw calculation, calculate the
!     needed sw reflectivities and absorptivities.
!-------------------------------------------------------------------
      if (trim(swform) /= 'esfsw99' ) then
        do j=JSRAD,JERAD
          do i=ISRAD,IERAD
	    if (ncldsw(i,j) > 0 ) then
	      do ngp=1,nsolwg
!-----------------------------------------------------------------------
!     call cldalb to define the zenith angle dependent cloud visible
!     and infrared reflectivities.
!-----------------------------------------------------------------------
	        call cldalb (cosangsolar(i,j,ngp))

!-----------------------------------------------------------------------
!     call albcld_sw to assign the zenith-angle-dependent visible and 
!     infrared reflectivities and infrared absorptivities to the model 
!     clouds.
!-----------------------------------------------------------------------
                call albcld_sw (i, j, hi_cloud, mid_cloud, low_cloud,  &
		 	        ngp, camtsw, cmxolw, crndlw, cvisrfsw, &
			        cirrfsw, cirabsw)
	      end do
            endif
	  end do
        end do
      endif

!--------------------------------------------------------------------
!  if microphysics is being used for either sw or lw calculation, call 
!  microphys_rad_driver to obtain cloud radiative properties.
!--------------------------------------------------------------------
      if (trim(swform) == 'esfsw99' ) then
        call microphys_rad_driver (camtsw, emmxolw, emrndlw, &
				   cldext=cldext, cldsct=cldsct,    &
				   cldasymm=cldasymm) 
      else if (do_lwcldemiss) then
        call microphys_rad_driver (camtsw, emmxolw, emrndlw)
      endif

!-------------------------------------------------------------------
!    deallocate  local arrays
!-------------------------------------------------------------------
      deallocate (cvflag)
      deallocate (cosangsolar  )
      deallocate (qlsig   )
      deallocate (ifcd    )
      deallocate (ifcv    )
      deallocate (ccover )
      deallocate (iflagglbl  )
      deallocate (hi_cloud)
      deallocate (mid_cloud)
      deallocate (low_cloud)



end subroutine rh_clouds_calc


!####################################################################


subroutine cldalb (zenith)

!---------------------------------------------------------------------
!     cldalb calculates a zenith angle dependency for the cloud albedos.
!     the cloud albedos are interpolated using data adapted from fritz 
!     (1954).  the solar zenith angle is the only input required.
!-----------------------------------------------------------------------

real, intent(in)           ::  zenith

!--------------------------------------------------------------------
      real                 :: zangle, remain
      integer              :: nband, indx

!--------------------------------------------------------------------
!     define zenith angle in degrees. for original Skyhi results, use 
!     a zenith angle specified as 60.00001 and the skyhi albedo values.
!-----------------------------------------------------------------------

        zangle = ACOS(zenith)*radians_to_degrees
!       zangle = 60.0000001

!-----------------------------------------------------------------------
!     define reflectivities for each cloud level.
!-----------------------------------------------------------------------
        if (zangle .GE. 80.0E+00) then
!-----------------------------------------------------------------------
!     if zenith angle is greater than 80 degrees, define the reflect-
!     ivities as those values in the table at 80 degrees (last entry).
!-----------------------------------------------------------------------
          do nband=1,NREFL_BDS
            zza(3,nband) = albch(NANGS,nband)
            zza(2,nband) = albcm(NANGS,nband)
            zza(1,nband) = albcl(NANGS,nband)
          end do  
        else
!-----------------------------------------------------------------------
!     if zenith angle is less than 80 degrees, interpolate albedos from 
!     tables for each cloud level.
!-----------------------------------------------------------------------
          indx   = IFIX(zangle/5.0E+00) + 1
          remain = AMOD(zangle, 5.0E+00)
          do nband=1,NREFL_BDS
            zza(3,nband) = albch(indx,nband) + (remain/5.0E+00)*  &
                           (albch(indx+1,nband) - albch(indx,nband))
            zza(2,nband) = albcm(indx,nband) + (remain/5.0E+00)*  &
                           (albcm(indx+1,nband) - albcm(indx,nband))
            zza(1,nband) = albcl(indx,nband) + (remain/5.0E+00)*   &
                           (albcl(indx+1,nband) - albcl(indx,nband))
          end do 
        endif



end subroutine cldalb




!##################################################################

subroutine albcld_lw(hi_cloud, mid_cloud, low_cloud,       &
	             cmxolw, crndlw, emmxolw, emrndlw)

!-----------------------------------------------------------------------
!     albcld_lw computes the lw cloud emissivities. This calculation is 
!     based on sigma and cloud thickness in the old scheme (cldht60) 
!     and sigma, cloud thickness and latitude in the new scheme 
!     (cldht93).
!-----------------------------------------------------------------------

real, dimension(:,:,:),   intent(in)    :: cmxolw, crndlw
real, dimension(:,:,:,:), intent(inout) :: emmxolw, emrndlw
logical, dimension(:,:,:),intent(in)    :: hi_cloud, mid_cloud,   &
					   low_cloud

!---------------------------------------------------------------------
!    local variables
!---------------------------------------------------------------------

       integer  ::  i, j, k, kk

!-----------------------------------------------------------------------
!     compute the emissivities for each cloud in the column. 
!-----------------------------------------------------------------------

       if (do_cldht60) then
         do k=KSRAD,KERAD
	   do j=JSRAD,JERAD
	     do i=ISRAD,IERAD
               if ((cmxolw(i,j,k) + crndlw(i,j,k) ) > 0.0) then
!-----------------------------------------------------------------------
!     case of a thick cloud. note that all thick clouds are given the 
!     properties of low clouds.
!-----------------------------------------------------------------------
                 if (cmxolw(i,j,k) .NE. 0.0E+00) then
	           emmxolw(i,j,k,:) = cldem(1)
                 endif
                 if (crndlw(i,j,k) .NE. 0.0) then
!-----------------------------------------------------------------------
!     case of a thin high cloud.
!-----------------------------------------------------------------------
		   if (hi_cloud(i,j,k)) then
	             emrndlw(i,j,k,:) = cldem(3)
!-----------------------------------------------------------------------
!     case of a thin middle cloud.
!-----------------------------------------------------------------------
                   else if (mid_cloud(i,j,k))  then
	             emrndlw(i,j,k,:) = cldem(2)
!-----------------------------------------------------------------------
!     case of a thin low cloud. note that thick clouds and low clouds 
!     have the same radiative properties. 
!-----------------------------------------------------------------------
                   else if (low_cloud(i,j,k)) then
	             emrndlw(i,j,k,:) = cldem(1)
		   endif
                 endif
	       endif
             end do
           end do
         end do
       endif

       if (do_cldht93 .or. do_cldhtskyhi) then
         do k=KSRAD,KERAD
	   do j=JSRAD,JERAD
	     do i=ISRAD,IERAD
               if ((cmxolw(i,j,k) + crndlw(i,j,k) ) > 0.0) then
!-----------------------------------------------------------------------
!     case of a thick cloud. note that thick cloud properties are deter-
!     mined by the height of the base of the thick cloud, so that if 
!     there are two adjacent cirrus level clouds, they are assigned
!     cirrus cloud properties, incontrast to the cldht60 treatment.
!-----------------------------------------------------------------------
                 if (cmxolw(i,j,k) .NE. 0.0E+00) then
!-----------------------------------------------------------------------
!     case of a thick high cloud.
!-----------------------------------------------------------------------
                   if (hi_cloud(i,j,k)) then
         	     if (Environment%using_fms_periphs) then
                       emmxolw(i,j,k,:) = cldem(3)*0.6
                     else if (Environment%using_sky_periphs) then
                       emmxolw(i,j,k,:) = cldem(3)
                     endif
!-----------------------------------------------------------------------
!     case of a thick middle cloud.
!-----------------------------------------------------------------------
                   else if (mid_cloud(i,j,k)) then
	             emmxolw(i,j,k,:) = cldem(2)
!-----------------------------------------------------------------------
!     case of a thick low cloud.
!-----------------------------------------------------------------------
                   else if (low_cloud (i,j,k)) then
	             emmxolw(i,j,k,:) = cldem(1)
                   endif
                 endif
	         if (crndlw(i,j,k) .NE. 0.0E+00) then
!-----------------------------------------------------------------------
!     case of a thin high cloud.
!-----------------------------------------------------------------------
                   if (hi_cloud(i,j,k)) then
	             if (do_full_black_cirrus) then
                       emrndlw(i,j,k,:) = cldem(3)
	             else if (do_part_black_cirrus) then
	               emrndlw(i,j,k,:) = 0.6E+00*cldem(3)
	             endif
!-----------------------------------------------------------------------
!     case of a thin middle cloud.
!-----------------------------------------------------------------------
                   else if (mid_cloud(i,j,k)) then 
	             emrndlw(i,j,k,:) = cldem(2)
!-----------------------------------------------------------------------
!     case of a thin low cloud.
!-----------------------------------------------------------------------
                   else if (low_cloud(i,j,k)) then
	             emrndlw(i,j,k,:) = cldem(1)
	           endif
                 endif
	       endif
             end do     
           end do
         end do
!------------------------------------------------------------------
!! for fms formulation, set thick cloud properties based on cloud base,
!  even if some cloud layers extend out of the cloud base type region
!-------------------------------------------------------------------

         if (Environment%using_fms_periphs) then 
           do k=KERAD,KSRAD+1,-1
	     do j=JSRAD,JERAD
	       do i=ISRAD,IERAD
                 if (cmxolw(i,j,k) /= 0.0) then
                   kk = k-1
                   if (cmxolw(i,j,kk) /= 0.0) then
                     emmxolw(i,j,kk,:) = emmxolw(i,j,k,:)
                   endif
                 endif
               end do
             end do
           end do
         endif
       endif
     
 
end subroutine albcld_lw


!####################################################################

subroutine albcld_sw(i,j, hi_cloud, mid_cloud, low_cloud, ngp,    &
		     camtsw, cmxolw, crndlw, cvisrfsw, cirrfsw, cirabsw)

!-----------------------------------------------------------------------
!     albcld_sw computes the cloud albedos. This calculation is based on
!     sigma and cloud thickness in the old scheme (cldht60) and sigma, 
!     cloud thickness  and latitude in the new scheme (cldht93).
!-----------------------------------------------------------------------

real, dimension(:,:,:),   intent(inout) :: camtsw, cmxolw, crndlw
real, dimension(:,:,:,:), intent(inout) :: cvisrfsw, cirrfsw, cirabsw
logical, dimension(:,:,:),intent(in)    :: hi_cloud, mid_cloud,   &
					   low_cloud
integer,                  intent(in)    :: i, j, ngp
!---------------------------------------------------------------------

       integer  ::  k, n, kk

!-----------------------------------------------------------------------
!     compute the reflectivities and absorptivities for each cloud in 
!     the column. cldhm and cldml are sigma levels which serve as 
!     boundaries between low, middle and high clouds. 
!-----------------------------------------------------------------------
       if (do_cldht60) then
         do k=KSRAD+1,KERAD
	   if (camtsw(i,j,k) .NE. 0.0) then
!-----------------------------------------------------------------------
!     case of a thick cloud. note that all thick clouds are given the 
!     properties of low clouds.
!-----------------------------------------------------------------------
             if (cmxolw(i,j,k) .NE. 0.0E+00) then
               cvisrfsw(i,j,k,ngp) = zza(1,1)
               cirrfsw(i,j,k,ngp) = zza(1,2)
               cirabsw(i,j,k,ngp) = cabir(1)
             endif
             if (crndlw(i,j,k) .NE. 0.0 ) then    

!-----------------------------------------------------------------------
!     case of a thin high cloud.
!-----------------------------------------------------------------------
               if ( hi_cloud(i,j,k)) then
                 cvisrfsw(i,j,k,ngp) = zza(3,1)
                 cirrfsw(i,j,k,ngp) = zza(3,2)
                 cirabsw(i,j,k,ngp) = cabir(3)

!-----------------------------------------------------------------------
!     case of a thin middle cloud.
!-----------------------------------------------------------------------
               else if (mid_cloud(i,j,k)) then
                 cvisrfsw(i,j,k,ngp) = zza(2,1)
                 cirrfsw(i,j,k,ngp) = zza(2,2)
                 cirabsw(i,j,k,ngp) = cabir(2)

!-----------------------------------------------------------------------
!     case of a thin low cloud. note that thick clouds and low clouds 
!     have the same radiative properties. 
!-----------------------------------------------------------------------
               else if (low_cloud(i,j,k)) then
                 cvisrfsw(i,j,k,ngp) = zza(1,1)
                 cirrfsw(i,j,k,ngp) = zza(1,2)
                 cirabsw(i,j,k,ngp) = cabir(1)
               endif
	     endif
           endif
         end do
       endif

       if (do_cldht93 .or. do_cldhtskyhi) then
         do k=KSRAD+1,KERAD
           if (camtsw(i,j,k) .NE. 0.0E+00) then
!-----------------------------------------------------------------------
!     case of a thick cloud. note that thick cloud properties are deter-
!     mined by the height of the base of the thick cloud, so that if 
!     there are two adjacent cirrus level clouds, they are assigned
!     cirrus cloud properties.
!-----------------------------------------------------------------------
             if (cmxolw(i,j,k) .NE. 0.0E+00) then
!-----------------------------------------------------------------------
!     case of a thick high cloud.
!-----------------------------------------------------------------------
               if (hi_cloud(i,j,k)) then
                 cvisrfsw(i,j,k,ngp) = zza(3,1)
                 cirrfsw(i,j,k,ngp) = zza(3,2)
	         if (Environment%using_fms_periphs) then
                   cirabsw(i,j,k,ngp) = MIN(0.99-cirrfsw(i,j,k,ngp), &
                                            cabir(3) )
                 else if (Environment%using_sky_periphs) then
                   cirabsw(i,j,k,ngp) = cabir(3)
                 endif

!-----------------------------------------------------------------------
!     case of a thick middle cloud.
!-----------------------------------------------------------------------
               else if (mid_cloud(i,j,k)) then
                 cvisrfsw(i,j,k,ngp) = zza(2,1)
                 cirrfsw(i,j,k,ngp) = zza(2,2)
	         if (Environment%using_fms_periphs) then
                   cirabsw(i,j,k,ngp) = MIN(0.99-cirrfsw(i,j,k,ngp), &
                                            cabir(2) )
                 else if (Environment%using_sky_periphs) then
                   cirabsw(i,j,k,ngp) = cabir(2)
                 endif

!-----------------------------------------------------------------------
!     case of a thick low cloud.
!-----------------------------------------------------------------------
               else if (low_cloud(i,j,k)) then
                 cvisrfsw(i,j,k,ngp) = zza(1,1)
                 cirrfsw(i,j,k,ngp) = zza(1,2)
	         if (Environment%using_fms_periphs) then
                   cirabsw(i,j,k,ngp) = MIN(0.99-cirrfsw(i,j,k,ngp), &
	                                    cabir(1) )
                 else if (Environment%using_sky_periphs) then
                   cirabsw(i,j,k,ngp) = cabir(1)
                 endif
               endif
             endif
	     if (crndlw(i,j,k) .NE. 0.0E+00) then

!-----------------------------------------------------------------------
!     case of a thin high cloud.
!-----------------------------------------------------------------------
               if (hi_cloud(i,j,k)) then
                 cvisrfsw(i,j,k,ngp) = zza(3,1)
                 cirrfsw(i,j,k,ngp) = zza(3,2)
	         if (Environment%using_fms_periphs) then
                   cirabsw(i,j,k,ngp) = MIN(0.99-cirrfsw(i,j,k,ngp), &
                                            cabir(3) )
                 else if (Environment%using_sky_periphs) then
                   cirabsw(i,j,k,ngp) = cabir(3)
                 endif

!-----------------------------------------------------------------------
!     case of a thin middle cloud.
!-----------------------------------------------------------------------
               else if (mid_cloud(i,j,k))  then
                 cvisrfsw(i,j,k,ngp) = zza(2,1)
                 cirrfsw(i,j,k,ngp) = zza(2,2)
	         if (Environment%using_fms_periphs) then
                   cirabsw(i,j,k,ngp) = MIN(0.99-cirrfsw(i,j,k,ngp), &
                                            cabir(2) )
                 else if (Environment%using_sky_periphs) then
                   cirabsw(i,j,k,ngp) = cabir(2)
                 endif

!-----------------------------------------------------------------------
!     case of a thin low cloud.
!-----------------------------------------------------------------------
               else if (low_cloud(i,j,k)) then
                 cvisrfsw(i,j,k,ngp) = zza(1,1)
                 cirrfsw(i,j,k,ngp) = zza(1,2)
	         if (Environment%using_fms_periphs) then
                   cirabsw(i,j,k,ngp) = MIN(0.99-cirrfsw(i,j,k,ngp), &
                                            cabir(1) )
                 else if (Environment%using_sky_periphs) then
                   cirabsw(i,j,k,ngp) = cabir(1)
                 endif
               endif
             endif
           endif
         end do     

!------------------------------------------------------------------
!! for fms formulation, set thick cloud properties based on cloud base,
!  even if some cloud layers extend out of the cloud base type region
!-------------------------------------------------------------------
         if (Environment%using_fms_periphs) then
           do k=KERAD,KSRAD+1,-1
             if (cmxolw(i,j,k) /= 0.0) then
               kk = k-1
               if (cmxolw(i,j,kk) /= 0.0) then
                 cvisrfsw(i,j,kk,ngp) = cvisrfsw(i,j,k,ngp)
                 cirrfsw(i,j,kk,ngp) = cirrfsw(i,j,k,ngp)
                 cirabsw(i,j,kk,ngp) = cirabsw(i,j,k,ngp)
               endif
             endif
           end do
         endif
       endif
     
 
end subroutine albcld_sw



	       end module rh_based_clouds_mod

