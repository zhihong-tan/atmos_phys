
                 module standalone_clouds_mod

use time_manager_mod,    only: time_type
use microphys_rad_mod,   only: microphys_rad_driver, microphys_rad_init
use utilities_mod,       only: open_file, file_exist,   &
                               check_nml_error, error_mesg,   &
                               print_version_number, FATAL, NOTE, &
			       WARNING, get_my_pe, close_file, &
                               get_my_pe, get_domain_decomp
use rad_step_setup_mod,  only: ISRAD, IERAD, JSRAD, JERAD, & 
                               KSRAD, KERAD
use rad_utilities_mod,   only: Environment, environment_type, &
			       longwave_control_type, Lw_control, &
			       shortwave_control_type, Sw_control
use longwave_setup_mod,  only: longwave_parameter_type, Lw_parameters
use constants_new_mod,   only: radians_to_degrees
use strat_clouds_W_mod,  only: strat_clouds_calc
use cloud_rad_mod,       only: cloud_rad_init
!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!	   standalone cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

! character(len=5), parameter  ::  version_number = 'v0.09'
  character(len=128)  :: version =  '$Id: standalone_clouds.F90,v 1.2 2001/07/05 17:33:29 fms Exp $'
  character(len=128)  :: tag     =  '$Name: eugene $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          standalone_clouds_init, standalone_clouds_driver

!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=10) :: cldht_type_form='   '
character(len=10) :: cloud_data_form='   '
character(len=16) :: cloud_overlap_form='   '
character(len=10) :: lhsw_cld_prop_form='   '
character(len=10) :: lw_cld_prop_form='   '
integer           :: cloud_data_points = 0

!  logical variables derived from namelist variables

logical     :: do_cldhtskyhi=.true.
logical     :: do_cldht60=.false.
logical     :: do_cldht93=.false.
logical     :: do_random_overlap = .false.
logical     :: do_max_random_overlap = .false.
logical     :: do_cloud_data_spec = .false.
logical     :: do_cloud_data_input = .false.
logical     :: do_lhsw_cld_prop_spec = .false.
logical     :: do_lhsw_cld_prop_input = .false.
logical     :: do_lw_cld_prop_spec = .false.
logical     :: do_lw_cld_prop_input = .false.

namelist /standalone_clouds_nml /     &
	  	          cldht_type_form, cloud_data_form, &
                          cloud_data_points, cloud_overlap_form, &
                          lhsw_cld_prop_form,  &
                          lw_cld_prop_form


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

character(len=10)      :: swform

!------------------------------------------------------------------
!    test (singlecolumn) values for cloud amount, height index,
!    infrared emissivity and shortwave reflectivity and absorptivity
!------------------------------------------------------------------

integer                               :: ich, icm, ict, icb
real                                  :: ch, cm, cl
real, dimension(:,:), allocatable     :: cvis_rf_in, cir_rf_in,  &
                                         cir_abs_in
real, dimension(:,:), allocatable     :: cloud_amount_in
real, dimension(:,:), allocatable     :: max_cloud_amount_in
real, dimension(:,:), allocatable     :: rnd_cloud_amount_in
integer, dimension(:), allocatable    :: nmax_cloud_in
integer, dimension(:), allocatable    :: nrnd_cloud_in
real, dimension(:,:,:), allocatable   :: emlw_band_in
real, dimension(:,:), allocatable     :: emlw_in

!---------------------------------------------------------------------
!    NLWCLDB is the actual number of frequency bands for which lw
!    emissitivies are defined. 
!---------------------------------------------------------------------

integer   :: NLWCLDB 
 
integer, parameter   :: LATOBS=19
integer              :: x(4), y(4), jd, jdf

!----------------------------------------------------------------------
!   define default values for shortwave cloud absorptivity and
!   reflectivity for use only in the lacis-hansen implementation
!   of shortwave radiative transfer
!      (these are the values previously used in SKYHI)
!----------------------------------------------------------------------
real                 :: lowcloud_refl_visband = 0.66E+00
real                 :: midcloud_refl_visband = 0.54E+00
real                 :: highcloud_refl_visband = 0.21E+00
real                 :: lowcloud_refl_nearirband = 0.50E+00
real                 :: midcloud_refl_nearirband = 0.46E+00
real                 :: highcloud_refl_nearirband = 0.19E+00
real                 :: lowcloud_abs_visband = 0.0E+00
real                 :: midcloud_abs_visband = 0.0E+00
real                 :: highcloud_abs_visband = 0.0E+00
real                 :: lowcloud_abs_nearirband = 0.30E+00
real                 :: midcloud_abs_nearirband = 0.20E+00
real                 :: highcloud_abs_nearirband = 0.04E+00

!----------------------------------------------------------------------
!   define default (grey) values for longwave emissivity for low, mid,
!   high clouds. these are used if no microphysics parameterizations
!   (or assumptions) are used in the longwave radiative transfer
!      (these are the values previously used in SKYHI)
!----------------------------------------------------------------------
real                 :: lowcloud_emiss = 1.00E+00
real                 :: midcloud_emiss = 1.00E+00
real                 :: highcloud_emiss = 1.00E+00

!----------------------------------------------------------------------
!----------------------------------------------------------------------



	 contains 



subroutine standalone_clouds_init (theta, do_strat_clouds,   &
                                          do_specified_clouds)

real, dimension(:), intent(in)    :: theta
logical,            intent(in)    :: do_strat_clouds
logical,            intent(in)    :: do_specified_clouds

!--------------------------------------------------------------------
!
!   do_strat_clouds (if true) indicates that cloud microphysical
!   properties are to be obtained from the Klein parameterization
!
!   do_specified_clouds (if true) indicates that the cloud microphysical
!   properties are specified.
!
!   one of these quantities must be true. neither quantity pertains
!   to the specification of cloud fractions and altitudes, although
!   the altitudes of the specified cloud microphysical properties
!   must be consistent with the altitude of the specified cloud(s).
!
!--------------------------------------------------------------------

      integer            :: unit, ierr, io
      integer            :: k, j, i, n
      integer            :: ni
      character(len=10)  :: type_form
      real               :: alatn, fl
      integer            :: li, jj
      integer            :: kk,ktop,kbot
      real               :: max_cld_calc

!--------------------------------------------------------------------
!     these variables define the boundaries (in sigma coordinates) 
!     between high and middle and middle and low clouds. 
!--------------------------------------------------------------------
      real, dimension(:), allocatable   :: cldhm, cldml
      real, dimension(:), allocatable   :: cldhm_abs, cldml_abs, &
					   cldhm_abs_gl, cldml_abs_gl
      real                              :: cldhp, cldhe, cldmp, cldme

!---------------------------------------------------------------------
!-----  read namelist  ------
  
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=standalone_clouds_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'standalone_clouds_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'standalone_clouds',   &
!						       version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=standalone_clouds_nml)
      endif
      call close_file (unit)

!-------------------------------------------------------------------
!  ensure that cloud_data_points (from namelist) is at least one
!-------------------------------------------------------------------
      if (cloud_data_points < 1) then
        call error_mesg( 'standalone_clouds_init',  &
                 ' cloud_data_points must be greater than zero', FATAL)
      endif
!--------------------------------------------------------------------
!  define the formulation to be used for defining the interfaces 
!  between high, middle and low clouds (using cldht_type_form)
!--------------------------------------------------------------------
      if (trim(cldht_type_form) == '60') then

!--------------------------------------------------------------------
!  use cloud height (low, middle, high) specification as in 1960
!  definition in SUPERSOURCE model
!--------------------------------------------------------------------
        do_cldht60 = .true.
      else if (trim(cldht_type_form) == '93') then

!--------------------------------------------------------------------
!  use cloud height (low, middle, high) specification as in 1993
!  definition in SUPERSOURCE model
!--------------------------------------------------------------------
        do_cldht93 = .true.
      else if (trim(cldht_type_form) == 'skyhi') then

!--------------------------------------------------------------------
!  use cloud height (low, middle, high) specification as in 1993
!  definition in SKYHI model
!--------------------------------------------------------------------
        do_cldhtskyhi = .true.
      else

!--------------------------------------------------------------------
! error condition
!--------------------------------------------------------------------
        call error_mesg( 'standalone_clouds_init',  &
                   ' cldht_type_form is not an acceptable value.', &  
							      FATAL)
      endif

!-------------------------------------------------------------------
!  define cloud overlap type (using cloud_overlap_form)
!-------------------------------------------------------------------
      if (trim(cloud_overlap_form) == 'random') then
        do_random_overlap = .true.
      else if (trim(cloud_overlap_form) == 'max_random') then
        do_max_random_overlap = .true.
      else
        call error_mesg( 'standalone_clouds_init',  &
                  ' cloud_overlap_form is not an acceptable value.', &
                                                              FATAL)
      endif
 
!-------------------------------------------------------------------
!  define cloud input data type
!-------------------------------------------------------------------
      if (trim(cloud_data_form) == 'input') then
        do_cloud_data_input = .true.
      else if (trim(cloud_data_form) == 'specified') then
        do_cloud_data_spec = .true.
      else
        call error_mesg( 'standalone_clouds_init',  &
                   ' cloud_data_form is not an acceptable value.', &
                                                              FATAL)
      endif

!-------------------------------------------------------------------
!  define cloud property data type
!-------------------------------------------------------------------
      if (trim(lhsw_cld_prop_form) == 'input') then
        do_lhsw_cld_prop_input = .true.
      else if (trim(lhsw_cld_prop_form) == 'specified') then
        do_lhsw_cld_prop_spec = .true.
      else
        call error_mesg( 'standalone_clouds_init',  &
               ' lhsw_cld_prop_form is not an acceptable value.', &
						              FATAL)
      endif
 
      if (trim(lw_cld_prop_form) == 'input') then
        do_lw_cld_prop_input = .true.
      else if (trim(lw_cld_prop_form) == 'specified') then
        do_lw_cld_prop_spec = .true.
      else
        call error_mesg( 'standalone_package_init',  &
               ' lw_cld_prop_form is not an acceptable value.', &
                                                              FATAL)
      endif

!---------------------------------------------------------------------
!  retrieve model size parameters from temporary module
!---------------------------------------------------------------------
      call get_domain_decomp (x, y)
      jd = y(2)
      jdf = y(4) - y(3) + 1

!---------------------------------------------------------------------
!    define the number of cloud emissivity bands for use in this module.
!    define the shortwave parameterization being used.
!---------------------------------------------------------------------

      NLWCLDB = Lw_parameters%NLWCLDB
      swform = Sw_control%sw_form

      allocate (cloud_amount_in (1:cloud_data_points,KSRAD:KERAD ))
      allocate (max_cloud_amount_in (1:cloud_data_points,KSRAD:KERAD ))
      allocate (rnd_cloud_amount_in (1:cloud_data_points,KSRAD:KERAD ))
      allocate (nmax_cloud_in (1:cloud_data_points))
      allocate (nrnd_cloud_in (1:cloud_data_points))
!   default cloud amounts are zero
      cloud_amount_in = 0.0E+00
      max_cloud_amount_in = 0.0E+00
      rnd_cloud_amount_in = 0.0E+00
      nmax_cloud_in = 0
      nrnd_cloud_in = 0
      if (do_cloud_data_input) then

        unit = open_file ('INPUT/cld_amt_file', action='read', &
                          form='formatted')
        do i=1,cloud_data_points
          read (unit,FMT = '(5e18.10)')            &
                     (cloud_amount_in(i,k),k=KSRAD,KERAD)
        enddo
        call close_file (unit)
      
!    compute max_random and random overlap cloud numbers, depending on
!    whether the maximum-random or the random overlap assumption
!    is made (and on the cloud fraction data)
        if (do_random_overlap) then
!   random = count of cloud number at gridpt i, max_random = 0
          nrnd_cloud_in = COUNT(cloud_amount_in > 0.0E+00,DIM=2)
          nmax_cloud_in = 0
          max_cloud_amount_in = 0.0E+00
          rnd_cloud_amount_in = cloud_amount_in
        else if (do_max_random_overlap) then
          do i=1,cloud_data_points
            ktop = KSRAD
            kbot = KSRAD
            do k=KSRAD,KERAD-1
!   the cloud code below accounts for clouds in the (ktop, kbot) layers.
!    cycle k to (kbot+1) to continue processing
             if (k  > KSRAD .AND. k <= kbot) CYCLE
             if (cloud_amount_in(i,k) > 0.0E+00) then
!   check how many layers below (ie, higher k-index) have clouds
                ktop = k
                kbot = k
                do kk=ktop+1,KERAD
!   exit the loop if a layer with zero clouds is hit; otherwise
!   increment kbot index
                  if ( cloud_amount_in(i,kk) == 0.0E+00) EXIT
                  kbot = kk
                enddo
                if (ktop == kbot) then
!    single-layer cloud (always random)
                  max_cloud_amount_in(i,k) = 0.0E+00
                  rnd_cloud_amount_in(i,k) = cloud_amount_in(i,k)
                  nrnd_cloud_in(i) = nrnd_cloud_in(i) + 1
                else
!   TEMPORARY:
!    max cloud amount = max(cloud amounts from ktop to kbot) ;
!    rnd cloud amount = zero
                  max_cld_calc = MAXVAL          &
                                     (cloud_amount_in(i,ktop:kbot))
                  max_cloud_amount_in(i,ktop:kbot) = max_cld_calc
                  rnd_cloud_amount_in(i,ktop:kbot) = 0.0E+00
                  nmax_cloud_in(i) = nmax_cloud_in(i) + 1
                endif
              else
!    give values for (nonexistent) clouds
                ktop = k
                kbot = k
              endif
            enddo
!   deal with special case of cloud in lowest layer, no cloud in
!   next lowest layer (should never happen)
            if (cloud_amount_in(i,KERAD) > 0.0E+00 .AND.     &
                cloud_amount_in(i,KERAD-1) == 0.0E+00   ) then
              rnd_cloud_amount_in(i,KERAD) = cloud_amount_in(i,KERAD)
              max_cloud_amount_in(i,KERAD) = 0.0E+00
              nrnd_cloud_in(i) = nrnd_cloud_in(i) + 1
            endif
          enddo
        else
!   some other cloud scheme (not defined now)
        endif
!
      else if (do_cloud_data_spec) then
!-------------------------------------------------------------------
!   define remaining cloud properties based as specified here.
!-------------------------------------------------------------------
      cloud_amount_in(:,:) = 0.0   ! default cloudiness is zero
      type_form = Environment%column_type
      if (trim(type_form) == 'skyl40') then
        ch   = 0.159E+00
        cm   = 0.070E+00
        cl   = 0.269E+00
        ich  = 29 + KSRAD - 1
        icm  = 34 + KSRAD - 1
        ict  = 35 + KSRAD - 1
        icb  = 37 + KSRAD - 1
      else if (trim(type_form) == 'nmcl18') then
        ch   = 0.159E+00
        cm   = 0.070E+00
        cl   = 0.269E+00
        ich  =  5 + KSRAD - 1
        icm  = 11 + KSRAD - 1
        ict  = 12 + KSRAD - 1
        icb  = 14 + KSRAD - 1
      else if (trim(type_form) == 'r30l14') then
        ch   = 0.159E+00
        cm   = 0.070E+00
        cl   = 0.269E+00
        ich  = 6 + KSRAD - 1
        icm  = 9 + KSRAD - 1
        ict  = 10 + KSRAD - 1
        icb  = 12 + KSRAD - 1
      else
	call error_mesg ('standalone_clouds_init', &
                 'cloud properties have not been specified for&
		 & this column type', FATAL)
      endif

!   point 1 contains 3 clouds (max-random) 
        rnd_cloud_amount_in(1,ich) = ch
	rnd_cloud_amount_in(1,icm) = cm
	do k=ict,icb
	max_cloud_amount_in(1,k) = cl
	enddo
	nmax_cloud_in(1) = 1
	nrnd_cloud_in(1) = 2
      if (cloud_data_points > 1) then
!   point 2 contains zero clouds
        max_cloud_amount_in(2,:) = 0.0
        rnd_cloud_amount_in(2,:) = 0.0
	nmax_cloud_in(2) = 0
	nrnd_cloud_in(2) = 0
      endif
      if (cloud_data_points > 2) then
!   point 3 contains 5 clouds (random)
        rnd_cloud_amount_in(3,ich) = ch
	rnd_cloud_amount_in(3,icm) = cm
	do k=ict,icb
	rnd_cloud_amount_in(3,k) = cl
	enddo
	nmax_cloud_in(3) = 0
	nrnd_cloud_in(3) = 5
      endif
    else
	call error_mesg ('standalone_clouds_init', &
                 'cloud data type has not been specified properly', &
		  FATAL)
    endif

!-------------------------------------------------------------------
!   allocate lw emissivity array and fill with either input or
!   specified values
!-------------------------------------------------------------------
      allocate (emlw_band_in(1:cloud_data_points, KSRAD:KERAD, NLWCLDB))

! default emissivity value is unity
      emlw_band_in(:,:,:) = 1.00E+00

      if (do_lw_cld_prop_input) then
	allocate (emlw_in (1:cloud_data_points,KSRAD:KERAD ))
        unit = open_file ('INPUT/lw_cld_prop_file', action='read', &
		         form='formatted')
        do i=1,cloud_data_points
        read (unit,FMT = '(5e18.10)')            &
                 (emlw_in(i,k),k=KSRAD,KERAD)
        enddo
	do n=1,NLWCLDB
	  emlw_band_in(:,:,n) = emlw_in(:,:)
	enddo
!   note that there must be consistency between the cloud fraction/ 
!        altitude values and the emissivity values. will have to 
!        handle this later

        call close_file (unit)
      
      else if (do_lw_cld_prop_spec) then
!   lw emissivities are given for the specified clouds and cloud
!   heights given above
        do n=1,NLWCLDB
!   point 1 contains 3 clouds (max-random) 
          emlw_band_in(1,ich,n) = highcloud_emiss
          emlw_band_in(1,icm,n) = midcloud_emiss
	  do k=ict,icb
            emlw_band_in(1,k,n) = lowcloud_emiss
	  enddo
	enddo
        if (cloud_data_points > 1) then
!   point 2 contains zero clouds
          emlw_band_in(2,:,:) = 1.00E+00
        endif
        if (cloud_data_points > 2) then
  	  do n=1,NLWCLDB
!   point 3 contains 5 clouds (random)
          emlw_band_in(3,ich,n) = highcloud_emiss
          emlw_band_in(3,icm,n) = midcloud_emiss
	    do k=ict,icb
              emlw_band_in(3,k,n) = lowcloud_emiss
	    enddo
          enddo
        endif
      else
	call error_mesg ('standalone_clouds_init', &
             'lw cld properties have not been specified correctly', &
	       FATAL)
      endif

      if (trim(swform) /= 'esfsw99') then
!-------------------------------------------------------------------
!   allocate lhsw cloud property arrays and either read input file
!   containing them or use specified values
!-------------------------------------------------------------------
	allocate (cvis_rf_in (1:cloud_data_points,KSRAD:KERAD+1 ))
	allocate (cir_rf_in (1:cloud_data_points,KSRAD:KERAD+1 ))
	allocate (cir_abs_in (1:cloud_data_points,KSRAD:KERAD+1 ))


! default lhsw cloud property value is unity reflection, zero absorption
      cvis_rf_in(:,:) = 1.00E+00
      cir_rf_in(:,:) = 1.00E+00
      cir_abs_in(:,:) = 0.0

      if (do_lhsw_cld_prop_input) then
        unit = open_file ('INPUT/lhsw_cld_prop_file', action='read', &
		         form='formatted')
        do i=1,cloud_data_points
        read (unit,FMT = '(5e18.10)')            &
                 (cvis_rf_in(i,k),k=KSRAD,KERAD)
        read (unit,FMT = '(5e18.10)')            &
                 (cir_rf_in(i,k),k=KSRAD,KERAD)
        read (unit,FMT = '(5e18.10)')            &
                 (cir_abs_in(i,k),k=KSRAD,KERAD)
        enddo
        call close_file (unit)
      
      else if (do_lhsw_cld_prop_spec) then
!   lhsw properties are given for the specified clouds and cloud
!   heights given above
!   point 1 contains 3 clouds (max-random) 
          cvis_rf_in(1,ich) = highcloud_refl_visband
          cvis_rf_in(1,icm) = midcloud_refl_visband
          cir_rf_in(1,ich) = highcloud_refl_nearirband
          cir_rf_in(1,icm) = midcloud_refl_nearirband
          cir_abs_in(1,ich) = highcloud_abs_nearirband
          cir_abs_in(1,icm) = midcloud_abs_nearirband
	  do k=ict,icb
            cvis_rf_in(1,k) = lowcloud_refl_visband
            cir_rf_in(1,k) = lowcloud_refl_nearirband
            cir_abs_in(1,k) = lowcloud_abs_nearirband
	  enddo
        if (cloud_data_points > 1) then
!   point 2 contains zero clouds
          cvis_rf_in(2,:) = 1.00E+00
          cir_rf_in(2,:) = 1.00E+00
          cir_abs_in(2,:) = 0.0
        endif
        if (cloud_data_points > 2) then
!   point 3 contains 5 clouds (random)
          cvis_rf_in(3,ich) = highcloud_refl_visband
          cvis_rf_in(3,icm) = midcloud_refl_visband
          cir_rf_in(3,ich) = highcloud_refl_nearirband
          cir_rf_in(3,icm) = midcloud_refl_nearirband
          cir_abs_in(3,ich) = highcloud_abs_nearirband
          cir_abs_in(3,icm) = midcloud_abs_nearirband
	    do k=ict,icb
              cvis_rf_in(3,k) = lowcloud_refl_visband
              cir_rf_in(3,k) = lowcloud_refl_nearirband
            cir_abs_in(3,k) = lowcloud_abs_nearirband
	    enddo
        endif
      else
	call error_mesg ('standalone_clouds_init', &
             'lhsw cld properties have not been specified correctly', &
	       FATAL)
      endif

      endif !  lhsw code

!-----------------------------------------------------------------------
!     define the cloud height boundaries that distinguish low, middle 
!     and high clouds. values are defined in terms of the standard 
!     sigma coordinate at full levels.
!---------------------------------------------------------------------
      allocate (cldhm (LATOBS) )
      allocate (cldml (LATOBS) )
      allocate (cldhm_abs (jdf) )
      allocate (cldml_abs (jdf) )
      allocate (cldhm_abs_gl (jd) )
      allocate (cldml_abs_gl (jd) )

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
      do j=1,jd
        if (Environment%using_sky_periphs) then
	  if ( jd > 1) then
            fl = 9.0E+00 - 9.0E+00*(FLOAT(jd+1 -j - jd/2) -    &
                 0.5E+00)/FLOAT(jd/2)
          else
	    fl = 9.0E+00
	  endif
          li = NINT(fl) + 1
          cldhm_abs_gl(j) = cldhm(li)
          cldml_abs_gl(j) = cldml(li)
        endif
      end do
      do j=1,jdf
        if (Environment%using_sky_periphs) then
	  cldhm_abs(j) = cldhm_abs_gl(j+y(3)-1)
	  cldml_abs(j) = cldml_abs_gl(j+y(3)-1)
        else if (Environment%using_fms_periphs) then
          cldhm_abs(j) = cldhp + (90.0E+00-abs(theta(j)*   &
	                 radians_to_degrees))*(cldhe-cldhp)/90.0E+00
          cldml_abs(j) = cldmp + (90.0E+00-abs(theta(j)*    &
	                 radians_to_degrees))*(cldme-cldmp)/90.0E+00
        endif
      end do
      deallocate (cldhm_abs_gl)
      deallocate (cldml_abs_gl)

!-------------------------------------------------------------------
!     if the strat clouds module is to be used, initialize the
!     cloud_rad module
!-------------------------------------------------------------------
      if (do_strat_clouds) then
        call cloud_rad_init
      endif

!-------------------------------------------------------------------
!    if a cloud microphysics scheme is to be employed with the cloud
!    scheme, initialize the microphysics_rad module.
!--------------------------------------------------------------------
      if (trim(swform) == 'esfsw99'  .or.     &
                                       Lw_control%do_lwcldemiss) then
	call microphys_rad_init                                  &
             (cldhm_abs_in=cldhm_abs, cldml_abs_in=cldml_abs)
      endif

!-------------------------------------------------------------------- 
!  deallocate arrays
!-------------------------------------------------------------------
      deallocate (cldhm)
      deallocate (cldml)
      deallocate (cldhm_abs)
      deallocate (cldml_abs)

!--------------------------------------------------------------------



end subroutine standalone_clouds_init




!#################################################################

subroutine standalone_clouds_driver (do_strat_clouds, &
				     camtsw, cmxolw, crndlw, ncldsw, &
				     nmxolw, nrndlw, emmxolw, emrndlw, &
				     Time_next,   &
				     cirabsw, cvisrfsw, cirrfsw, &
				     cldext, cldsct, cldasymm )

real,    dimension(:,:,:),   intent(inout) :: camtsw, cmxolw, crndlw
integer, dimension(:,:),     intent(inout) :: ncldsw, nrndlw, nmxolw
real,    dimension(:,:,:,:), intent(inout) :: emmxolw, emrndlw
type(time_type),                intent(in), optional ::  Time_next
real,    dimension(:,:,:,:), intent(inout), optional ::       &
                                              cirabsw, cirrfsw, &
                                              cvisrfsw, cldext, &
                                              cldsct, cldasymm
logical                   , intent(in)     :: do_strat_clouds

!----------------------------------------------------------------------
!   local variables
!----------------------------------------------------------------------

      integer       :: i, j, k, n

!-----------------------------------------------------------------------
!    subroutine default_clouds (called previously from clouddrvr)
!    has initialized all cloud property fields
!-----------------------------------------------------------------------

!---------------------------------------------------------------------- 
!     define the cloud and cloud index fields used by both longwave
!     and shortwave radiation.
!----------------------------------------------------------------------
 
      do j=JSRAD,JERAD
        do i=1,cloud_data_points
          crndlw(i,j,:) = rnd_cloud_amount_in(i,:)
          cmxolw(i,j,:) = max_cloud_amount_in(i,:)
          nmxolw(i,j) = nmax_cloud_in(i)
          nrndlw(i,j) = nrnd_cloud_in(i)
        enddo
      enddo

      camtsw(:,:,:) = crndlw(:,:,:) + cmxolw(:,:,:)
      ncldsw(:,:) = nmxolw(:,:) + nrndlw(:,:) 

      do j=JSRAD,JERAD
        do i=1,cloud_data_points
          emmxolw(i,j,:,:) = emlw_band_in(i,:,:)
          emrndlw(i,j,:,:) = emlw_band_in(i,:,:)
          if (present (cirabsw) ) then
            do k=KSRAD,KERAD+1
              cirabsw(i,j,k,:) = cir_abs_in(i,k)
              cirrfsw(i,j,k,:) = cir_rf_in(i,k)
              cvisrfsw(i,j,k,:) = cvis_rf_in(i,k)
            enddo
          endif
        enddo
      enddo

!--------------------------------------------------------------------
!  call cloudrad_driver to obtain cloud radiative properties from
!  specified or computed microphysical cloud properties
!--------------------------------------------------------------------
      if (do_strat_clouds) then
        call strat_clouds_calc  &
                             (camtsw, cmxolw, crndlw, ncldsw, nmxolw, &
                              nrndlw, emmxolw, emrndlw,           &
                              Time_next,                  &
                              cldext=cldext, &
                              cldsct=cldsct, cldasymm=cldasymm)

      else    ! (if not strat_clouds)
        if (trim(swform) == 'esfsw99' ) then
          call microphys_rad_driver (camtsw, emmxolw, emrndlw,   &
                                     cldext=cldext, cldsct=cldsct,   &
                                     cldasymm=cldasymm )
         else if (Lw_control%do_lwcldemiss) then
           call microphys_rad_driver ( camtsw, emmxolw, emrndlw)
         endif
      endif

!--------------------------------------------------------------------


end subroutine standalone_clouds_driver



!###################################################################



	       end module standalone_clouds_mod

