	       module donner_deep_mod
 
use  utilities_mod,       only:  open_file, file_exist, read_data, &
                                 write_data, check_nml_error,    &
                                 error_mesg,   &
                                 print_version_number, FATAL, NOTE,  &
                                 WARNING, get_my_pe, close_file, &
                                 get_domain_decomp
use time_manager_mod,     only:  time_type, increment_time, &
			         set_time, operator(>=), operator(+),&
	  		         get_time, get_calendar_type, &
                                 operator(/=), operator(-), &
				 operator(>), operator(==)
use diag_manager_mod,     only:  register_diag_field, send_data
use sat_vapor_pres_mod,   only:  lookup_es


!implicit none
private

!--------------------------------------------------------------------
!         module to compute the effects of deep convection
!
!         Primary remaining work:
!               optimization / cleanup of code -- current code adds 
!               ~ 40% to total N30L40 FMS model time, after donner_deep
!               is spunup.
!     
!
!--------------------------------------------------------------------





!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------


character(len=128)  :: version =  '$Id: donner_deep.F90,v 1.4 2002/02/22 18:59:37 fms Exp $'
character(len=128)  :: tag     =  '$Name: galway $'


!--------------------------------------------------------------------
!---interfaces------

public   &
	 donner_deep_init, donner_deep, fill_donner_c_variables, &
	 get_cemetf, get_tprea1_donner_deep,   &
	 store_donner_deep_variables, &
	 fill_donner_deep_variables, donner_deep_time_vary,     &
	 write_restart_donner_deep, inquire_donner_deep, donner_deep_end

private   &
	cupar_vect, precu_vect, polat_vect, cape_vect, ieq_x, satad, &
	lcl, mulsub_vect, cloudm_vect, simult, meens, mesub, micro,  &
	prean, andge, write_diagnostics, STRAT_CLOUD_DONNER_TEND


!---------------------------------------------------------------------
!---namelist----

logical :: initialize_donner_deep = .false.    ! should the scheme be 
			                       ! initialized now ?
logical :: save_donner_deep_diagnostics=.true. ! should module diag-
					       ! nostics be archived ? 

integer :: donner_deep_freq   = 1080  ! frequency of calling donner_deep
                                      ! assumption made that <= 86400 s
integer :: donner_deep_offset = 720   ! offset of first daily calcul-
				      ! ation time from 00Z

logical :: debug = .false.          ! should debug data be written ?
integer :: kttest = -1              ! number of "standard" steps after 
				    ! start of job to do debug. for 
				    ! constant delta t, kttest = 
				    ! (time interval from start at
			            ! which debug is desired)/ delta t.
integer :: itest = -1               ! i index of debug column
integer :: jtest = -1               ! global j index of debug column
integer :: ktest_model = 100        ! debug data will be printed from
				    ! model level ktest_model to kma

integer :: ntracers_in_usrtr_file =  43   ! number of tracers in rstrt
					  ! file - applicable only in
					  ! SKYHI


namelist / donner_deep_nml /      &
                              initialize_donner_deep, &
			      save_donner_deep_diagnostics, &
			      donner_deep_freq, donner_deep_offset, &
			      debug, kttest, itest, jtest, ktest_model,&
			      ntracers_in_usrtr_file




!--------------------------------------------------------------------
!--- public data ----------




!--------------------------------------------------------------------
!----private data-----------

!--------------------------------------------------------------------
!---list of restart versions readable by this module--

integer, dimension(2)  :: restart_versions = (/ 1, 2 /)
!   version 1 does not have the tprea1_3d array included

!--------------------------------------------------------------------
!--- 3d module arrays (k index 1 is closest to the ground) ---

       
!     cecon          normalized cell condensation/deposition
!                    [ kg(H2O)/kg/sec ]
!     ceefc          normalized cell entropy-flux convergence [ K/sec ]
!                    (excludes convergence of surface flux) Entropy-flux
!                    convergence divided by (p0/p)**(rd/cp).
!     cememf         normalized moisture forcing, cells+meso 
!                    [ kg(H2O)/kg/sec ]
!     cememf_mod     normalized moisture forcing, cells+meso, after 
!                    any modifications needed to prevent the creation
!                    of negative values [ kg(H2O)/kg/sec ]
!     cemetf         normalized thermal forcing, cells+meso [ K/sec ]
!                    (excludes convergence of surface heat flux)
!                    Cumulus thermal forcing defined as in Fig. 3 of 
!                    Donner (1993, JAS).
!     cemfc          normalized cell moisture-flux convergence
!                    (excludes convergence of surface moisture flux)
!                    [ kg(H2O)/kg/sec ]
!     cmus_3d        normalized mesoscale-updraft deposition
!                    [ kg(H2O)/kg/sec]
!     cual_3d        cloud fraction, cells+meso, normalized by a(1,p_b)
!     dmeml_3d       mass flux in mesoscale downdraft [ kg/((m**2) s) ]
!                    (normalized by a(1,p_b)) 
!     ecds_3d        normalized convective downdraft evaporation
!                    [ kg(H2O)/kg/sec ]
!     eces_3d        normalzed convective-updraft evporation/sublimation
!                    [ kg(H2O)/kg/sec ]
!     elt_3d         normalized melting [ K/sec ]
!     emds_3d        normalized mesoscale-downdraft sublimation
!                    [ kg(H2O)/kg/sec ]
!     emes_3d        normalized mesoscale-updraft sublimation
!                    [ kg(H2O)/kg/sec ]
!     fre_3d         normalized freezing [ K/sec ]
!     qmes_3d        normalized mesoscale moisture-flux convergence
!                    [ kg(H2O)/kg/sec ]
!     tmes_3d        normalized mesoscale entropy-flux convergence
!                    [ K/sec ]
!                    Entropy-flux convergence is mesoscale component
!                    of second term in expression for cumulus thermal
!                    forcing in Fig. 3 of Donner (1993, JAS).
!     uceml_3d       normalized mass fluxes in cell updrafts
!                    [ kg/((m**2)*s ] 
!     umeml_3d       mass flux in mesoscale updraft [ kg/((m**2) s) ]
!                    (normalized by a(1,p_b)) 
!     wmms_3d        normalized mesoscale deposition of water vapor from
!                    cells [ kg(H2O)/kg/sec ]
!     wmps_3d        normalized mesoscale redistribution of water vapor
!                    from cells [ kg(H2O)/kg/sec ]
!     xice_3d        mesoscale ice mass mixing ratio (kg(ice)/kg)
!     dgeice_3d      mesoscale ice generalized effective size, defined 
!                    as in Fu  (1996, J. Clim.) (micrometers)
!                    index 1 at model top
!     cuqi_3d        cell ice content (kg(ice)/kg)
!                    index 1 at model top
!     cuql_3d        cell liquid content (kg{water)/kg)
!                    index 1 at model top

real, dimension(:,:,:), allocatable  :: cemetf        !  ooo(24)
real, dimension(:,:,:), allocatable  :: ceefc         !  ooo(25)
real, dimension(:,:,:), allocatable  :: cecon         !  ooo(26)
real, dimension(:,:,:), allocatable  :: cemfc         !  ooo(27)
real, dimension(:,:,:), allocatable  :: cememf        !  ooo(28)
real, dimension(:,:,:), allocatable  :: cememf_mod              
real, dimension(:,:,:), allocatable  :: cual_3d       !  ooo(29)
real, dimension(:,:,:), allocatable  :: fre_3d        !  ooo(30)
real, dimension(:,:,:), allocatable  :: elt_3d        !  ooo(31)
real, dimension(:,:,:), allocatable  :: cmus_3d       !  ooo(32)
real, dimension(:,:,:), allocatable  :: ecds_3d       !  ooo(33)
real, dimension(:,:,:), allocatable  :: eces_3d       !  ooo(34)
real, dimension(:,:,:), allocatable  :: emds_3d       !  ooo(35)
real, dimension(:,:,:), allocatable  :: emes_3d       !  ooo(36)
real, dimension(:,:,:), allocatable  :: qmes_3d       !  ooo(37)
real, dimension(:,:,:), allocatable  :: wmps_3d       !  ooo(38)
real, dimension(:,:,:), allocatable  :: wmms_3d       !  ooo(39)
real, dimension(:,:,:), allocatable  :: tmes_3d       !  ooo(40)
real, dimension(:,:,:), allocatable  :: dmeml_3d      !  ooo(41)
real, dimension(:,:,:), allocatable  :: uceml_3d      !  ooo(42)
real, dimension(:,:,:), allocatable  :: umeml_3d      !  ooo(43)
real, dimension(:,:,:), allocatable  :: xice_3d       !
real, dimension(:,:,:), allocatable  :: dgeice_3d       !
real, dimension(:,:,:), allocatable  :: cuqi_3d       !
real, dimension(:,:,:), allocatable  :: cuql_3d       !

!--------------------------------------------------------------------
!--- 2d module arrays (k index 1 is closest to the ground) ---

!     a1_3d            fractional area of index-1 cu subensemble
!     amax_3d          maximum value for a_1(p_b)
!                      See "a Bounds 6/7/97" notes
!     amos_3d          upper limit on cloud fractional area based on
!                      moisture constraint See "Moisture Constraint," 
!                      8/8/97.
!     ampta1_3d        area weighted mesoscale cloud fraction, normal-
!                      ized by a(1,p_b)
!     coin_3d          convective inhibition 
!                      energy required to lift parcel from level istart
!                      to level of free convection. [ J/kg ]
!     contot_v         ratio of convective to total precipitation
!     dcape_3d         time rate of change of xcape_3d [ J/(kg s) ]
!     emdi_v           vertical integral of mesoscale-downdraft 
!                      sublimation
!     omint_3d         time-integrated low level dispalcement [ Pa ]
!     plcl_3d          pressure at lifting condensation level [ Pa ]
!     plfc_3d          pressure at level of free convection [ Pa ]
!                      height of plfc .le. height of plcl. if parcel 
!                      becomes buoyant below plcl, cin can be .lt. 0
!     plzb_3d          pressure at level of zero buoyancy [ Pa ]
!     qint_3d          vertically integrated column moisture
!                      [ kg (h20)/(m**2) ]
!     rcoa1_3d         area weighted convective precipitation rate
!                      [ mm/day ]
!     tprea1_3d        area weighted total normalized precipitation 
!                      [ mm/day ]
!     xcape_3d         convective available potential energy. energy 
!                      released as parcel moves from level of free 
!                      convection to level of zero buoyancy [ J/kg ]

real, dimension(:,:), allocatable  :: plcl_3d      !  spare(20)
real, dimension(:,:), allocatable  :: plfc_3d      !  spare(21)
real, dimension(:,:), allocatable  :: plzb_3d      !  spare(22)
real, dimension(:,:), allocatable  :: xcape_3d     !  spare(23)
real, dimension(:,:), allocatable  :: coin_3d      !  spare(24)
real, dimension(:,:), allocatable  :: dcape_3d     !  spare(25)
real, dimension(:,:), allocatable  :: qint_3d      !  spare(26)
real, dimension(:,:), allocatable  :: a1_3d        !  spare(27)
real, dimension(:,:), allocatable  :: amax_3d      !  spare(28)
real, dimension(:,:), allocatable  :: amos_3d      !  spare(29)
real, dimension(:,:), allocatable  :: tprea1_3d    !  spare(30)
real, dimension(:,:), allocatable  :: ampta1_3d    !  spare(31)
real, dimension(:,:), allocatable  :: omint_3d     !  spare(32)
real, dimension(:,:), allocatable  :: rcoa1_3d     !  spare(33)
real, dimension(:,:), allocatable  :: contot_v     
real, dimension(:,:), allocatable  :: emdi_v       

!------------------------------------------------------------------
!   module loop and dimension variables

integer :: iminp=1       ! initial physics window i index
integer :: jminp=1       ! initial physics window j index
integer :: imaxp, jmaxp  ! physics_window sizes
integer :: idf, jdf      ! subdomain sizes
integer :: id, jd        ! global domain sizes
integer :: nlev          ! number of model half-levels (layers)

!------------------------------------------------------------------
!   module control variables

logical  :: do_donner_deep=.false. ! is this module active ?
logical  :: do_init         !  is initialization necessary ?
integer  :: tstep_counter   !  number of timesteps executed in this job
integer  :: total_pts       !  total number of points in the subdomain
integer  :: num_pts         !  number of pts processed in the subdomain
			    !  so far on this time step
integer  :: num_pts2        !  number of pts processed in the subdomain
			    !  so far on this time step
logical  :: ldeep           !  is this a step to calculate donner deep ?
logical  :: ldeep_next      !  is deep calculated on the next step ?
logical  :: standard_step   !  is this a "standard" time step ? (cur-
			    !  rently only meaningful in SKYHI)

!---------------------------------------------------------------------
!  constants required for deep cumulus parameterization. these will
!  ultimately be "used" from a constants module.
 
real, parameter   ::  rgas = 8.314320E+03    ! universal gas constant
					     ! [ J/(kg K) ]
real, parameter   ::  wtmair = 2.896440E+01  ! molecular wt of air 
real, parameter   ::  wtmh2o = 1.801534E+01  ! molecular wt of water
real, parameter   ::  cpi    = 1.004840E+03  ! specific heat of dry air 
					     ! at constant pressure 
					     ! [ J/(kg K) ]
real, parameter   ::  gravm  = 9.806650      ! acceleration of gravity
                                             ! [ m/(sec**2) ]
real, parameter   ::  cpv    = 1850.         ! specific heat of water
					     ! vapor at constant pres-
					     ! sure [ J/(kg K) ]
real              ::  rair   = rgas/wtmair   ! gas constant for dry air
					     ! [ J/(kg K) ]
real, parameter   ::  rvap   = rgas/wtmh2o   ! gas constant for water
					     ! vapor [ J/(kg K) ]
real              ::  rh2o   = 461.          ! gas constant for water
					     ! vapor [ J/(kg K) ]
real, parameter   ::  latvap = 2.5104e06     ! latent heat of 
					     ! vaporization [ J/kg ]
real, parameter   ::  rocp   = 0.622         ! ratio of molecular 
					     ! weights of water vapor 
					     ! and dry air
real              ::  epsilo = 0.622         ! ratio of molecular 
					     ! weights of water vapor 
					     ! and dry air
real, parameter   ::  frezdk = 273.16        ! melting temperature [ K ]

!--------------------------------------------------------------------
!  module time and time-step variables

real            :: dt                    !  model time step [ sec ]
type(time_type) :: Next_donner_deep_time !  next time to calculate 
					 !  donner_deep
type(time_type) :: Donner_deep_timestep  !  frequency of calculating
					 !  donner_deep
type(time_type) :: Old_time_step         !  frequency of calculating
					 !  donner_deep, as read from
					 !  restart file (subject to
					 !  modification to namelist
					 !  supplied value)

!--------------------------------------------------------------------
!   module internal parameters

integer, parameter :: kpar=7             !  number of cumulus    
					 !  subensembles
integer, parameter :: ncap=100           !  number of levels in cloud 
					 !  model
real, parameter    :: pdeep_cv = 500.e02 !  required pressure difference
					 !  between GCM level closest to
					 !  ground and pressure at level
					 !  of zero buoyancy for deep 
                                         !  convection  to occur (Pa).
real, parameter    :: cdeep_cv = 10.     !  maximum value of convective 
					 !  inhibition (J/kg) that 
					 !  allows convection. Value of 
					 !  10 suggested by Table 2 in 
					 !  Thompson et al. (1979, JAS).

real, parameter    :: pdeep_mc = 200.e02 !  pressure thickness (Pa) 
					 !  required for mesoscale circ-
					 !  ulation. It refers to the 
					 !  least penetrative ensemble 
					 !  member. For this check to 
					 !  function properly, the en-
					 !  trainment coefficient 
                                         !  in Cloudm for kou=1 must be
					 !  the largest entrainment
					 !  coefficient. 


!--------------------------------------------------------------------
!   variables for debugging option

logical       :: in_debug_window     ! indicates if debug column is
				     ! in current window
integer       :: jdebug              ! window j index of debug column

!--------------------------------------------------------------------
!----variables for diagnostics:-----

integer          :: id_cemetf_deep, id_ceefc_deep, id_cecon_deep, &
		    id_cemfc_deep, id_cememf_deep, id_cememf_mod_deep, &
		    id_cual_deep, id_fre_deep, id_elt_deep, &
		    id_cmus_deep, id_ecds_deep, id_eces_deep, &
		    id_emds_deep, id_emes_deep, id_qmes_deep,&
		    id_wmps_deep, id_wmms_deep, id_tmes_deep,&
		    id_dmeml_deep, id_uceml_deep, &
		    id_umeml_deep, id_xice_deep, id_dgeice_deep, &
                    id_cuqi_deep, id_cuql_deep, &
		    id_plcl_deep, id_plfc_deep, id_plzb_deep, &
		    id_xcape_deep, id_coin_deep,  &
		    id_dcape_deep, id_qint_deep, id_a1_deep, &
		    id_amax_deep, id_amos_deep, &
		    id_tprea1_deep, id_ampta1_deep, &
		    id_omint_deep, id_rcoa1_deep

real             :: Missing_value = -999.
character(len=8) :: mod_name = 'don_deep'






!---------------------------------------------------------------------
!---------------------------------------------------------------------




	          contains




subroutine donner_deep_init (kmax_in, Time, axes)

!---------------------------------------------------------------------
!     initialization subroutine for module
!---------------------------------------------------------------------

integer,               intent(in)           :: kmax_in
type(time_type),       intent(in), optional :: Time
integer, dimension(4), intent(in), optional :: axes

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      kmax_in        number of model levels
!
!  intent(in), optional variables:
!
!    these variables are present when running in FMS, not in SKYHI:
!
!      Time           current Time
!      axes           data axes for diagnostics
!
!-------------------------------------------------------------------


      integer, dimension(4)   :: x, y
      integer                 :: unit, ierr, io
      integer                 :: secs_from_start, time_to_donner_deep, &
				 time_remaining, secs, days
      character(len=4)        :: chvers
      type(time_type)         :: New_donner_deep_time
  
!---------------------------------------------------------------------
!-----  read namelist  ------
  
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=donner_deep_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'donner_deep_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
      if (get_my_pe() == 0)  then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit,nml=donner_deep_nml)
      endif
      call close_file (unit)

!--------------------------------------------------------------------
!  execute the remainder of this subroutine only if do_donner_deep is
!  true.
!--------------------------------------------------------------------
!     if (do_donner_deep) then

!---------------------------------------------------------------------
!  test that namelist values are valid
!---------------------------------------------------------------------
        if (donner_deep_freq > 86400.) then
	  call error_mesg ( 'donner_deep_init', &
	    ' donner convection must be called at least once per day', &
						       FATAL)
        endif
        if (donner_deep_freq <= 0.) then
	  call error_mesg ( 'donner_deep_init', &
	    ' a positive value must be assigned to donner_deep_freq', &
						       FATAL)
        endif
        if (donner_deep_offset < 0.) then
	  call error_mesg ( 'donner_deep_init', &
	    ' a non-negative value must be given donner_deep_offset', &
						       FATAL)
        endif
        if (donner_deep_offset >= 86400.) then
	  call error_mesg ( 'donner_deep_init', &
	    ' donner_deep_offset must be less than a full day', &
						       FATAL)
        endif

!-------------------------------------------------------------------
!  define module variables specifying the grid dimensions
!   idf, jdf = dimensions of subdomain on this processor
!   id, jd   = total model dimensions
!   x(1), x(2), y(1), y(2) = beginning and ending full model dimensions
!   x(3), x(4), y(3), y(4) = beginning and ending domain indices present
!                            on this processor
!   NOTE: on a unitasked PVP machine with a slab decomposition, 
!   (1) = (3) and (2) = (4), since only 1 processor is employed
!-------------------------------------------------------------------
        nlev = kmax_in
        call get_domain_decomp (x, y)
        idf = x(4) - x(3) + 1
        jdf = y(4) - y(3) + 1
        id  = x(2) - x(1) + 1
        jd  = y(2) - y(1) + 1

!---------------------------------------------------------------------
!  verify that requested debug indices are valid
!---------------------------------------------------------------------
        if (debug) then
	  if (itest < x(1) .or. itest > x(2)) &
	    call error_mesg ('donner_deep_init',   &
	    'desired column for debug diagnostics is not in domain',  &
	 					     FATAL)
	  if (jtest < y(1) .or. jtest > y(2)) &
	    call error_mesg ('donner_deep_init',   &
	    'desired row for debug diagnostics is not in domain', FATAL)

      	  if (ktest_model < 1 .or. ktest_model > nlev) &
            call error_mesg ('donner_deep_init',   &
   	   'cutoff level for debug diagnostics is not in domain', FATAL)
        endif

!--------------------------------------------------------------------
!  initialize the es table to be used with this module. AS OF NOW, this
!  has been imported from SKYHI; the ultimate plan is to use the sat-
!  uration vapor pressure module in FMS.
!--------------------------------------------------------------------
!       call esinit

!---------------------------------------------------------------------
!  initialize various points processed counters. define the total num-
!  ber of columns present in this subdomain. initialize the timestep
!  counter (needed only in SKYHI).
!---------------------------------------------------------------------
        num_pts   = 0
        num_pts2  = 0
        total_pts = idf*jdf
        tstep_counter = 0

!---------------------------------------------------------------------
!  define the next time at which this module is to be called, using
!  the namelist supplied frequency and offset from 00Z at which 
!  donner_deep is to be called. assumption is made that donner_deep is
!  called at least once per day. define a time_type version of the
!  donner_deep_freq, and initialize Old_time_step to this value, if
!  there is no restart file. if a restart file exists, the value found
!  there will be used, provided the donner_deep_freq has not changed.
!---------------------------------------------------------------------
!!! CHECK OUT ALL POSSIBLE CASES HERE TO VERIFY THIS CODE IS
!!  GENERAL ENOUGH for various relations between time, offset and
!!  frequency.  3/6/01 : appears OK, as modified. handles all tested
!!  cases properly. 3/7/01 : additional slight mods made in response
!!  to running in FMS.
        if (present (Time) ) then
          call get_time (Time, secs, days)
	  if (secs > donner_deep_offset) then
	    secs_from_start = MOD (secs, donner_deep_freq)
	    time_to_donner_deep = donner_deep_freq - secs_from_start
	    time_remaining = MOD (time_to_donner_deep, donner_deep_freq)
	  else
            time_remaining = donner_deep_offset - secs
	    if (time_remaining == 0 .and. initialize_donner_deep) then
	      time_remaining = donner_deep_freq
            endif
	  endif
          Next_donner_deep_time = increment_time (Time,    &
				   time_remaining, 0)
	  Donner_deep_timestep = set_time (donner_deep_freq,0)
	  if (.not. file_exist ('INPUT/donner_deep.res') ) then
	    Old_time_step = set_time (donner_deep_freq, 0)
          endif
        endif

!--------------------------------------------------------------------
!   allocate module variables that must be saved across timesteps.
!--------------------------------------------------------------------
        call allocate_variables

!--------------------------------------------------------------------
!   set a flag to indicate whether module being coldstarted. Send
!   message if it is being coldstarted even though a restart file is
!   in place.to advise user.
!--------------------------------------------------------------------
        if (initialize_donner_deep) then
	  do_init = .true.
          if (file_exist ('INPUT/donner_deep.res') ) then
	    if (get_my_pe() == 0) then
	      call error_mesg ('donner_deep_init', &
	          ' BE ADVISED: donner_deep_mod is being coldstarted, &
	             &even though a restart file is present. Is&
	             & this what you want ?', NOTE)
	    endif
	  endif
	else 
	  do_init = .false.
	endif

!--------------------------------------------------------------------
!   provide initial values for the module variables that are saved
!   across timesteps. these will come from either an initialization
!   subroutine, a module restart file or an old SKYHI usrtr file. if
!   none of these options is indicated, a message is written and the 
!   job aborts.
!--------------------------------------------------------------------
        if (do_init) then
          call initialize_variables
        else if (file_exist ('INPUT/donner_deep.res') ) then
	  call read_restart_donner_deep
        else if (file_exist ('INPUT/usrtr_rstrt') ) then
	  call read_usrtr_file
        else
	   call error_mesg ('donner_deep_init', &
          'either provide restart file (old or new) or request&
	  & initialization', FATAL)
        endif

!--------------------------------------------------------------------
!    initialize netcdf diagnostic fields only when running in FMS.
!-------------------------------------------------------------------
        if (present (Time) ) then
          if (save_donner_deep_diagnostics) then
            id_cemetf_deep =  &
           register_diag_field (mod_name, 'cemetf_deep', axes(1:3),   &
 	   Time, 'nrmlzd heating rate, c + m ', 'K/s',   &
           missing_value=missing_value)
            id_ceefc_deep =  &
           register_diag_field (mod_name, 'ceefc_deep', axes(1:3),   &
 	   Time, 'nrmlzd cell entrpy flx cnvrgnc', 'K/s',   &
           missing_value=missing_value)
            id_cecon_deep =  &
           register_diag_field (mod_name, 'cecon_deep', axes(1:3),   &
 	   Time, 'nrmlzd cell cond/evap ', 'kg(h2o)/kg/s',   &
           missing_value=missing_value)
            id_cemfc_deep =  &
           register_diag_field (mod_name, 'cemfc_deep', axes(1:3),   &
 	   Time, 'nrmlzd cell moist flx cnvgnc', 'kg(h2o)/kg/s',   &
           missing_value=missing_value)
            id_cememf_deep =  &
           register_diag_field (mod_name, 'cememf_deep', axes(1:3),   &
 	   Time, 'nrmlzd moistening rate, c + m ', 'kg(h2o)/kg/s',   &
           missing_value=missing_value)
            id_cememf_mod_deep =  &
           register_diag_field (mod_name, 'cememf_mod_deep', axes(1:3),&
 	   Time, 'mod cememf due to negative q ', 'kg(h2o)/kg/s',   &
           missing_value=missing_value)
            id_cual_deep =  &
           register_diag_field (mod_name, 'cual_deep', axes(1:3),   &
 	   Time, 'nrmlzd c + m cld frac ', 'percent',   &
           missing_value=missing_value)
            id_fre_deep =  &
           register_diag_field (mod_name, 'fre_deep', axes(1:3),   &
 	   Time, 'nrmlzd freezing ', 'K/sec',   &
           missing_value=missing_value)
            id_elt_deep =  &
           register_diag_field (mod_name, 'elt_deep', axes(1:3),   &
 	   Time, 'nrmlzd melting', 'K/sec',   &
           missing_value=missing_value)
            id_cmus_deep =  &
           register_diag_field (mod_name, 'cmus_deep', axes(1:3),   &
 	   Time, 'nrmlzd meso-up deposition', 'kg(h2o)/((m**2)sec)',   &
           missing_value=missing_value)
            id_ecds_deep =  &
           register_diag_field (mod_name, 'ecds_deep', axes(1:3),   &
 	   Time, 'nrmlzd convective dwndrft evap ', 'kg(h2o)/kg/sec', &
           missing_value=missing_value)
            id_eces_deep =  &
           register_diag_field (mod_name, 'eces_deep', axes(1:3),   &
 	   Time, 'nrmlzd convective updrft evap/subl ',    &
	   'kg(h2o)/kg/sec',   missing_value=missing_value)
            id_emds_deep =  &
           register_diag_field (mod_name, 'emds_deep', axes(1:3),   &
 	   Time, 'nrmlzd meso-dwn subl ', 'kg(h2o)/kg/sec',   &
           missing_value=missing_value)
            id_emes_deep =  &
           register_diag_field (mod_name, 'emes_deep', axes(1:3),   &
 	   Time, 'nrmlzd meso-up subl ', 'kg(h2o)/kg/sec',   &
           missing_value=missing_value)
            id_qmes_deep =  &
           register_diag_field (mod_name, 'qmes_deep', axes(1:3),   &
 	   Time, 'nrmlzd meso moist flux conv', 'kg(h2o)/kg/sec',   &
           missing_value=missing_value)
            id_wmps_deep =  &
           register_diag_field (mod_name, 'wmps_deep', axes(1:3),   &
 	   Time, 'nrmlzd meso redistrib of vapor from cells',    &
	   'kg(h2o)/kg/sec', missing_value=missing_value)
            id_wmms_deep =  &
           register_diag_field (mod_name, 'wmms_deep', axes(1:3),   &
 	   Time, 'nrmlzd meso depo of vapor from cells',    &
	   'kg(h2o)/kg/sec',  missing_value=missing_value)
            id_tmes_deep =  &
           register_diag_field (mod_name, 'tmes_deep', axes(1:3),   &
 	   Time, 'nrmlzd meso entropy flux conv',  'K/sec',   &
           missing_value=missing_value)
            id_dmeml_deep =  &
           register_diag_field (mod_name, 'dmeml_deep', axes(1:3), &
 	   Time, 'nrmlzd mass flux meso dwndrfts', 'kg/((m**2) s)',   &
           missing_value=missing_value)
            id_uceml_deep =  &
           register_diag_field (mod_name, 'uceml_deep', axes(1:3), &
 	   Time, 'nrmlzd mass flux cell updrfts', 'kg/((m**2) s)',   &
           missing_value=missing_value)
            id_umeml_deep =  &
           register_diag_field (mod_name, 'umeml_deep', axes(1:3), &
 	   Time, 'nrmlzd mass flux meso updrfts', 'kg/((m**2) s)',   &
           missing_value=missing_value)
	    id_xice_deep =  &
           register_diag_field (mod_name, 'xice_deep', axes(1:3),  &
           Time, 'meso ice mass mixing ratio ', 'kg(ice)/kg',   &
           missing_value=missing_value)
            id_dgeice_deep =  &
           register_diag_field (mod_name, 'dgeice_deep', axes(1:3), &
           Time, 'meso ice gen eff size ', 'micrometers',   &
           missing_value=missing_value)
            id_cuqi_deep =  &
           register_diag_field (mod_name, 'cuqi_deep', axes(1:3),  &
           Time, 'cell ice ', 'kg(H2O)/kg',   &
           missing_value=missing_value)
            id_cuql_deep =  &
           register_diag_field (mod_name, 'cuql_deep', axes(1:3),  &
           Time, 'cell liquid ', 'kg(H2O)/kg',   &
           missing_value=missing_value)

            id_plcl_deep =  &
           register_diag_field (mod_name, 'plcl_deep', axes(1:2),   &
 	   Time, 'pressure at lcl ', 'Pa ',   &
           missing_value=missing_value)
            id_plfc_deep =  &
           register_diag_field (mod_name, 'plfc_deep', axes(1:2),   &
 	   Time, 'pressure at lfc ', 'Pa ',   &
           missing_value=missing_value)
            id_plzb_deep =  &
           register_diag_field (mod_name, 'plzb_deep', axes(1:2),   &
 	   Time, 'pressure at lzb ', 'Pa ',   &
           missing_value=missing_value)
            id_xcape_deep =  &
           register_diag_field (mod_name, 'xcape_deep', axes(1:2),  &
 	   Time, 'cape', 'J/kg',   &
           missing_value=missing_value)
            id_coin_deep =  &
           register_diag_field (mod_name, 'coin_deep', axes(1:2),   &
 	   Time, 'convective inhibition ', 'J/kg',   &
           missing_value=missing_value)
            id_dcape_deep =  &
           register_diag_field (mod_name, 'dcape_deep', axes(1:2), &
 	   Time, 'time tendency of cape ', 'J/kg/sec',   &
           missing_value=missing_value)
            id_qint_deep =  &
           register_diag_field (mod_name, 'qint_deep', axes(1:2),   &
 	   Time, 'column moisture ', 'kg(h2o)/m**2',   &
           missing_value=missing_value)
            id_a1_deep =  &
           register_diag_field (mod_name, 'a1_deep', axes(1:2),   &
 	   Time, 'fractional area of cu subensemble ', 'percent',   &
           missing_value=missing_value)
            id_amax_deep =  &
           register_diag_field (mod_name, 'amax_deep', axes(1:2),   &
 	   Time, 'fractional area of largest cu subensemble ',  &
	   'percent',  missing_value=missing_value)
            id_amos_deep =  &
           register_diag_field (mod_name, 'amos_deep', axes(1:2),   &
 	   Time, 'uppr lmt on frac area from moisture', 'percent',   &
           missing_value=missing_value)
            id_tprea1_deep =  &
           register_diag_field (mod_name, 'tprea1_deep', axes(1:2), &
 	   Time, 'area wtd total precip ', 'mm/day',   &
           missing_value=missing_value)
            id_ampta1_deep =  &
           register_diag_field (mod_name, 'ampta1_deep', axes(1:2), &
 	   Time, 'meso cld frac', 'percent',   &
           missing_value=missing_value)
            id_omint_deep =  &
           register_diag_field (mod_name, 'omint_deep', axes(1:2), &
 	   Time, 'integrated low-lvl displ', 'Pa ',   &
           missing_value=missing_value)
            id_rcoa1_deep =  &
           register_diag_field (mod_name, 'rcoa1_deep', axes(1:2),  &
 	   Time, 'area wtd cnvctv precip ', 'mm/day',   &
           missing_value=missing_value)
          endif

!--------------------------------------------------------------------
!  override previously calculated next time to call donner_deep_mod if 
!  timestep as input from namelist differs from that in restart file.
!--------------------------------------------------------------------
          if (Donner_deep_timestep /= Old_time_step) then
	    New_donner_deep_time = Next_donner_deep_time -    &
				   Old_time_step + Donner_deep_timestep
	    if (New_donner_deep_time > Time) then
	      if (get_my_pe() == 0) then
	        call error_mesg ('donner_deep_init',   &
	             'donner_deep time step has changed, &
	             &next donner_deep time also changed', NOTE)
	      end if
              Next_donner_deep_time = New_donner_deep_time
            endif
          endif

!--------------------------------------------------------------------
!   if running in SKYHI, open a diagnostics file, if desired.
!--------------------------------------------------------------------
        else ! present(time)
          if (save_donner_deep_diagnostics) then
            unit = open_file (        &
		    file='history/donner_deep_diagnostics.his', &
!                   form='native', action='write')
                    form='unformatted', action='write')
            call close_file (unit)
          endif
        endif  !  present(time)
!     endif     !  if do_donner_deep

!--------------------------------------------------------------------
!  set flag to indicate this module has been activated
!--------------------------------------------------------------------
      do_donner_deep = .true.

!     do k=1,100 
!     tttt = 253.0 + (k-1)*0.01
!     call establ (es, tttt)
!     call lookup_es ( tttt, ess)
!     print *, 'temp, es(SKY), es (FMS) : ', tttt, es, ess
!     end do
!     stop

!------------------------------------------------------------------




end subroutine donner_deep_init

 

!#####################################################################

subroutine inquire_donner_deep (do_donner_deep_out)

!--------------------------------------------------------------------
!   inquire_donner_deep passes out a flag indicating whether this 
!   module has been activated.
!--------------------------------------------------------------------

logical, intent(out) :: do_donner_deep_out

!---------------------------------------------------------------------
!  intent(out) variables:
!
!    do_donner_deep_out    flag indicating whether donner_deep has 
!                          been activated
!
!---------------------------------------------------------------------

do_donner_deep_out = do_donner_deep


end subroutine inquire_donner_deep


!####################################################################

subroutine donner_deep_time_vary (ifcst, idt_in, nflip, dt_in, Time_in)

!---------------------------------------------------------------------
!     donner_deep_time_vary defines time and time-step dependent 
!     quantities, primarily related to controlling module execution.
!---------------------------------------------------------------------

integer,         intent(in), optional    :: ifcst, idt_in, nflip
real,            intent(in), optional    :: dt_in
type(time_type), intent(in), optional    :: Time_in

!---------------------------------------------------------------------
!   intent(in), optional variables:
!
!     the following variables are present when running in SKYHI:
!
!        ifcst     elapsed time in seconds from beginning of run
!        idt_in    model time step in seconds
!        nflip     indicator of current time step type (euler backward
!                  or not)
!
!     the following variables are present when running in FMS:
!     
!        dt_in     atmospheric time step in seconds
!        Time_in   current time 
!
!---------------------------------------------------------------------

       type(time_type)    :: Time_next, Next_check_time

!--------------------------------------------------------------------
!  this code is executed when in SKYHI:
!--------------------------------------------------------------------
       if (present (ifcst)) then

!-------------------------------------------------------------------
!   save the model timestep. verify that donner_deep_freq is an integral
!   multiple of it.
!--------------------------------------------------------------------
	 dt = real (idt_in)
	 if (MOD ( donner_deep_freq, idt_in) /= 0) then
	   call error_mesg ('donner_deep_time_vary',  &
		'donner_deep timestep NOT an integral multiple of &
		 &physics timestep', FATAL)
	 endif

!---------------------------------------------------------------------
!   define a flag indicating if donner_deep is to be calculated on this
!   timestep. define a flag indicating if it is to be calculated on
!   the next timestep. define a flag indicating if this is a leapfrog
!   step, as opposed to an euler backward step. increment the timestep
!   counter on "normal" steps.
!---------------------------------------------------------------------
	 if (do_init) then
	   ldeep = .false.
         else
	   ldeep      = (MOD(ifcst, donner_deep_freq) ==       &
						   donner_deep_offset )
	 endif
	 ldeep_next = (MOD(ifcst + idt_in, donner_deep_freq) ==     &
					           donner_deep_offset )
	 if (nflip .eq. 0 .or. nflip .eq. 2) then
           standard_step = .true.
	   tstep_counter = tstep_counter + 1
         else
           standard_step = .false.
         endif

!--------------------------------------------------------------------
!  this code is executed when in FMS:
!--------------------------------------------------------------------
       else if (present (Time_in) ) then

!-------------------------------------------------------------------
!   save the model timestep. verify that donner_deep_freq is an integral
!   multiple of it.
!--------------------------------------------------------------------
         dt = dt_in
	 idt = int (dt_in)
         if (MOD ( donner_deep_freq, idt) /= 0) then
	     call error_mesg ('donner_deep_time_vary',  &
		'donner_deep timestep NOT an integral multiple of &
		 &physics timestep', FATAL)
         endif

         if (idt == donner_deep_freq) then
	   if (do_init) then

!--------------------------------------------------------------------
!    if the model timestep equals the donner_deep timestep, and the
!    donner_deep module is being coldstarted, do not execute the module
!    code on this step, since the time tendency of cape is not avail-
!    able. set the next time to execute this code to be one model 
!    (= donner_deep timestep) from now. set the time to check for the 
!    step before the next donner_deep step to also be one model step 
!    from now.
!---------------------------------------------------------------------
	     ldeep = .false.
	     Next_check_time = Time_in + Donner_deep_timestep
	     Next_donner_deep_time  = Time_in + Donner_deep_timestep
           else if (Time_in == Next_donner_deep_time) then

!--------------------------------------------------------------------
!    if the model timestep equals the donner_deep timestep, and the
!    donner_deep module is NOT being coldstarted, and the current time 
!    is the Next_donner_deep_time, set the flag to indicate execution 
!    on this step and reset the Next_time to be 1 donner_deep step from
!    now.
!---------------------------------------------------------------------
	       ldeep = .true.
	       Next_check_time = Next_donner_deep_time +    &
				  Donner_deep_timestep 
           else

!--------------------------------------------------------------------
!    if the model timestep equals the donner_deep timestep, and the
!    donner_deep module is NOT being coldstarted, and the current time 
!    is NOT the Next_donner_deep_time, write an error message and abort.
!---------------------------------------------------------------------
	       call error_mesg ( 'donner_deep_time_vary' ,  &
		  'idt=donner_dt, not do_init and Times are /=', FATAL)
           endif
         else

!--------------------------------------------------------------------
!    if the model timestep does not equal the donner_deep timestep, 
!    set the flag to indicate either execution or non-execution on this
!    step, depending on whether or not the current time equals the next
!    donner_deep time. set the time to be checked against the step 
!    preceding the donner_deep step to be the next_donner_deep_time.
!---------------------------------------------------------------------
	   if (do_init) then
	     ldeep = .false.
           else
             if (Time_in == Next_donner_deep_time) then
	       ldeep = .true.
             else
	       ldeep = .false.
             endif
	   endif
           Next_check_time = Next_donner_deep_time 
         endif

!--------------------------------------------------------------------
!    determine if donner_deep is to be activated on the next model
!    time step. set the ldeep_next flag appropriately.
!--------------------------------------------------------------------
         Time_next = increment_time( Time_in, idt, 0)
         if (Time_next == Next_check_time) then
           ldeep_next = .true.
         else
           ldeep_next = .false.
         endif
          
!------------------------------------------------------------------
!    define standard_step and the time step counter appropriately.
!    these are currently only needed in SKYHI.
!------------------------------------------------------------------
         standard_step = .true.
	 tstep_counter = tstep_counter + 1
	 
!--------------------------------------------------------------------
!  if this code is executed, then an error is present, since neither
!  ifcst nor Time is present:
!--------------------------------------------------------------------
       else 
         call error_mesg ('donner_deep_time_vary', &
       ' either Time or ifcst must be passed as input argument', FATAL)
       endif


!--------------------------------------------------------------------
	   

end subroutine donner_deep_time_vary




!###################################################################

subroutine donner_deep (is, js, temold, ratold, pfull, phalf, coldT,   &
		        omega_btm, iflag, ttnd, qtnd, rain, snow, &
			Time, Lbot, cf, qlin,qiin, qatend, qltend, qitend)

!-------------------------------------------------------------------
!   donner_deep is the prognostic interface between the model and the 
!   donner_deep module and adds the tendencies due to deep convection
!   to the temperature and moisture fields.
!-------------------------------------------------------------------

!--------------------------------------------------------------------
integer,                    intent(in)           :: is, js
!real, dimension(:,:,:),     intent(in)           :: press
real, dimension(:,:,:),     intent(in)           :: pfull, phalf
logical, dimension(:,:),    intent(in)           :: coldT
real, dimension(:,:),       intent(in)           :: omega_btm
!real, dimension(:,:,:),     intent(inout)        :: temold, ratold
real, dimension(:,:,:),     intent(in)        :: temold, ratold
real, dimension(:,:,:),     intent(in),optional        :: qlin,qiin
real, dimension(:,:,:),     intent(out)       :: ttnd, qtnd
real, dimension(:,:),       intent(out)       :: rain, snow
integer,dimension(:,:,:),   intent(inout)        :: iflag
type(time_type),            intent(in), optional :: Time
real, dimension(:,:,:),     intent(out), optional:: qatend, qltend, qitend
!real, dimension(:,:,:),     intent(out), optional:: cfdel, qldel, qidel
real, dimension(:,:,:),     intent(in), optional :: cf
integer, dimension(:,:), intent(in), optional   :: Lbot
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!     is              lowest i subdomain index in the physics window 
!                     being integrated
!     js              lowest j subdomain index in the physics window 
!                     being integrated
!     pfull           pressure array, k dimension 1:nlev; contains
!                     model full level pressures 
!                     (i,j) dimensions are size of physics_window
!     phalf           pressure array, k dimension 1:nlev+1; contains
!                     model half level pressures for k indices 1:nlev+1
!                     (i,j) dimensions are size of physics_window
!     omega_btm       omega field at lowest model level
!                     (i,j) dimensions are size of physics_window
!     qliN             on entry, large-scale cloud liquid specific
!                     humidity (kg(liquid)/kg)
!     qiin            on entry, large-scale cloud ice specific
!                     humidity (kg(ice)/kg)
!     cf              on entry, large-scale cloud fraction (0-1)
!
!   intent(inout) variables:
!
!     temold          on entry, the current tau+1 temperature field;
!                     on exit, the tau+1 temperature field after the
!                     effects of donner_deep have been included
!                     (i,j) dimensions are size of physics_window, 
!                     k dimension is nlev
!     ratold          on entry, the current tau+1 moisture field;
!                     on exit, the tau+1 moisture field after the
!                     effects of donner_deep have been included
!                     (i,j) dimensions are size of physics_window, 
!                     k dimension is nlev
!     iflag           on entry, array is zeroed out;
!                     on exit contains flag indicating grid boxes
!                     where donner_deep occurs.
!                     (i,j) dimensions are size of physics_window, 
!                     k dimension is nlev
!
!   intent(inout), optional variables:
!
!    the following variable is present when running in FMS:
!
!     Time            current time (time_type)
!
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!   local variables:
       real, dimension(size(temold,1), size(temold,2)) ::    &
			         coin_vect, plcl_vect, plfc_vect,  &
				 plzb_vect, xcape_vect, qint_save,  &
				 qint_vect, qlsd_vect, dqls_vect, &
	                         xcapet, omint2_vect

       real, dimension(size(temold,1), size(temold,2), ncap) ::    &
			         rpc_vect, tpca_vect, pcape_vect,   &
				 rcape_vect, tcape_vect

       real, dimension (size(temold,1), size(temold,2),    &
						 size(temold,3) ) ::  &
	                         prini_vect, rini_vect, tini_vect,  &
				 rsave, dmeso_3d, xliq_3d

       real, dimension (size(temold,1), size(temold,2),    &
                        size(temold,3)+1)    :: mhalf_3d

       real ::   ampu, contot, emdi, tprei, przm, prent, prztm, dgeicer 
       real ::   dnna, dnnb
       real, dimension(size(temold,3))  :: prinp, rmuf, trf, xice  

       integer, dimension (size(temold,1) )     :: iabs
       integer, dimension (size(temold,2) )     :: jabs
       logical                                  :: used

!--------------------------------------------------------------------
!   local variables:
!
!      ampu       fractional area of mesoscale circulation
!      contot     ratio of convective to total precipitation
!      dnna       weight for averaging full-level mass fluxes to half levels
!      dnnb       weight for averaging full-level mass fluxes to half levels
!      emdi       vertical integral of  mesoscale-downdraft 
!                 sublimation (mm/d)
!      prinp      pressures at full levels (index 1 at model top) (Pa)
!      rmuf       mesoscale updraft mass flux  (kg/[(m**2) s])
!      trf        temperature at full GCM levels (K) 
!                 index 1 at model top
!      tprei      total convective-system precipitation (mm/d)
!      xice       anvil ice (kg(ice)/kg)
!      przm       pressure at anvil base (Pa), inferred from presence 
!                 of ice at full GCM levels
!      prztm      pressure at anvil top (Pa), inferred from presence 
!                 of ice at full GCM levels
!      dgeicer    generalized effective size of ice crystals in anvil 
!                 (micrometers) defined as in Fu (1996, J. Clim.)
!
!---------------------------------------------------------------------





!       call tic ( 'cvct', '1')

!-------------------------------------------------------------------
!    define i and j dimensions of physics window and allocate 
!    diagnostic arrays so dimensioned.
!-------------------------------------------------------------------
       imaxp = size (temold,1)
       jmaxp = size (temold,2)
  
!-------------------------------------------------------------------
!    define arrays which map physics window indices to the subdomain
!    indices. 
!-------------------------------------------------------------------
       do i=iminp,imaxp
         iabs(i) = is + i - 1
       end do
       do j=jminp,jmaxp
         jabs(j) = js + j - 1
       end do

!-------------------------------------------------------------------
!    allocate those module variables dimensioned by the size of the 
!    physics window. once allocated, the space for these variables 
!    remains for the duration of the model run.  initialize these
!    variables after allocating them.
!-------------------------------------------------------------------
       if (.not. allocated (plcl_3d) ) then
         call allocate_local_variables
         call initialize_local_variables
       endif

!-------------------------------------------------------------------
!    local variables are initialized upon each entry into the module
!    whem running in FMS. In Skyhi, initialization is not done, since
!    values have been retrieved from disk and are needed for diagnostics
!    performed on non-donner_deep steps. initialization here would prevent 
!    the calculation of those diagnostics. it has been verified that
!    obtaining these values from disk results in no changes to the 
!    prognostic integration of the model, indicating that these var-
!    iables are not needed from the previous timestep.
!-------------------------------------------------------------------
       if (present (Time) ) then
         call initialize_local_variables
       endif

!--------------------------------------------------------------------
!    initialize debug control in this subdomain. if debug is desired, 
!    this is the requested debug timestep, and data from this window
!    is desired as output, in_debug_window will be true and jdebug
!    is set to the debug column's j index. 
!-------------------------------------------------------------------
       in_debug_window = .false.
       jdebug = 0
       if (debug) then
         if (tstep_counter .ge. kttest) then
           do j=jminp,jmaxp
             if (jabs(j) == jtest) then
               print *, 'DONNER_DEEP/donner_deep: time, standard_step',    &
			             tstep_counter, standard_step
	       in_debug_window = .true.
	       jdebug = j
	       exit
             endif
           end do
         endif
       endif
!       call toc ( 'cvct', '1')

       ttnd(:,:,:) = 0.0
       qtnd(:,:,:) = 0.0
       qatend(:,:,:) = 0.0
       qltend(:,:,:) = 0.0
       qitend(:,:,:) = 0.0

!---------------------------------------------------------------------
!   calculate time integrated low-level displacement (Cu Closure E 
!   notes, 8 feb 98). omint_3d is the time integrated value, omint2_vect
!   is the value to be used on this step, which is the time integrated
!   value, unless the current displacement is downward. print out debug 
!   quantities, if desired.
!---------------------------------------------------------------------
!       call tic ( 'cvct', '2')
       if (standard_step) then
         do j=jminp,jmaxp
           do i=iminp,imaxp
             omint_3d(i,jabs(j)) = omint_3d(i,jabs(j)) +   &
	 			   omega_btm(i,j)*dt
             omint_3d(i,jabs(j)) = MIN (0.0, omint_3d(i,jabs(j)) )
	     omint_3d(i,jabs(j)) = MAX (omint_3d(i,jabs(j)),  &
!				  -press(i,j,nlev+1) )
					  -phalf(i,j,nlev+1) )
             omint2_vect(i,j)    = omint_3d(i,jabs(j) )
             if (omega_btm(i,j) > 0.) omint2_vect(i,j) = 0. 
           end do
         end do

         if (in_debug_window) then
           print *, 'DONNER_DEEP/donner_deep: dt, itest, jtest= ', dt, itest, jtest
	   print *, 'DONNER_DEEP/donner_deep: sfcprs, omint, omega_btm= ',   &
!                 press(itest,jdebug,nlev+1),   &
	                 phalf(itest,jdebug,nlev+1),   &
			 omint_3d(itest,jtest),&
 	                 omega_btm(itest,jdebug) 
         endif
       end if  !  standard_step
!       call toc ( 'cvct', '2')

!---------------------------------------------------------------------
!   DONNER_DEEP CALCULATION TIMESTEP CODE
!   execute the following code only on steps where convective forcing
!   is being calculated.
!---------------------------------------------------------------------
       if (ldeep) then

!--------------------------------------------------------------------
!   call precu to prepare single-column arrays for Donner (1993, JAS)
!   deep cumulus parameterization
!--------------------------------------------------------------------
!         call tic ( 'cvct', '4')
!        call precu_vect (temold, ratold, press, pcape_vect, &
         call precu_vect (temold, ratold, pfull, pcape_vect, &
                          prini_vect, jabs, rcape_vect, rini_vect,   &
		          tcape_vect, tini_vect)
!         call toc ( 'cvct', '4')

!--------------------------------------------------------------------
!   call cape to evaluate the column soundings for Donner (1993, JAS)
!   deep cumulus parameterization
!--------------------------------------------------------------------
!         call tic ( 'cvct', '5')
         call cape_vect (pcape_vect, rcape_vect, tcape_vect,   &
			 coin_vect, plcl_vect, plfc_vect, plzb_vect, &
			 rpc_vect, tpca_vect, xcape_vect, jabs)
!         call toc ( 'cvct', '5')

!---------------------------------------------------------------------
!     calculate integrated vapor.
!--------------------------------------------------------------------
!         call tic ( 'cvct', '6')
         do i=iminp,imaxp
           do j=jminp,jmaxp
             qint_vect(i,j) = rcape_vect(i,j,1)*(pcape_vect(i,j,1) -  &
	 		      pcape_vect(i,j,2))
             do k=2,ncap-1
               qint_vect(i,j) = qint_vect(i,j) + rcape_vect(i,j,k)*  &
			        0.5*(pcape_vect(i,j,k-1) -  &
				     pcape_vect(i,j,k+1))
             end do
             qint_vect(i,j) = qint_vect(i,j) + rcape_vect(i,j,ncap)*  &
			      (pcape_vect(i,j,ncap-1) -  &
		               pcape_vect(i,j,ncap))
             qint_vect(i,j) = qint_vect(i,j)/gravm
           end do
         end do

!---------------------------------------------------------------------
!    if desired, print out debug quantities.
!---------------------------------------------------------------------
         if (in_debug_window) then
           print *, 'DONNER_DEEP/donner_deep: plfc,plzb= ',  &
	         	    plfc_vect(itest,jdebug),  & 
			    plzb_vect(itest,jdebug)
           print *, 'DONNER_DEEP/donner_deep: plcl,coin,xcape= ',   &
                            plcl_vect(itest,jdebug),  &
			    coin_vect(itest,jdebug),   &
			    xcape_vect(itest,jdebug)
           print *,  'DONNER_DEEP/donner_deep: qint= ',qint_vect(itest,j)  
           do k=1,ncap      
             print *,   'DONNER_DEEP/donner_deep: k,pcape= ',k,  &
			     pcape_vect(itest,jdebug,k)
             print   *, 'DONNER_DEEP/donner_deep: k,tcape,tpca= ',k,   &
	 	  	     tcape_vect(itest,jdebug,k),   &
			     tpca_vect(itest,jdebug,k)
	     print *,   'DONNER_DEEP/donner_deep: k,rcape,rpca= ',k,   &
			     rcape_vect(itest,jdebug,k),    &
			     rpc_vect(itest,jdebug,k)
             if (pcape_vect(itest,jdebug,k) <     &
				 plzb_vect(itest,jdebug))  exit
           end do
         endif
!         call toc ( 'cvct', '6')

!---------------------------------------------------------------------
!   STANDARD TIMESTEP CODE
!   execute the following code on "standard" timesteps.
!---------------------------------------------------------------------
         if (standard_step) then
!           call tic ( 'cvct', '7')

!---------------------------------------------------------------------
!    save variables from local to module variables.
!---------------------------------------------------------------------
           do i=iminp,imaxp
             do j=jminp,jmaxp
 	       plcl_3d(i,j) = plcl_vect(i,j)
   	       plfc_3d(i,j) = plfc_vect(i,j)
 	       plzb_3d(i,j) = plzb_vect(i,j)
 	       coin_3d(i,j) = coin_vect(i,j)
 	       qint_save(i,j) = qint_3d(i,jabs(j))
	       dqls_vect(i,j) = (qint_vect(i,j) - qint_save(i,j ))/dt
	       qlsd_vect(i,j) = qint_vect(i,j)/donner_deep_freq
	       dcape_3d(i,j) = (xcape_vect(i,j) -    &
					       xcape_3d(i,jabs(j)))/dt

!  ASK LEO as to utility of retaining:
!ljd
!              if (dcape_3d(i,j) .gt. 0.) then
!	         xcapet(i,j) = xcape_vect(i,j)
!                if (xcapet(i,j) .gt. 500.) xcapet(i,j)=500.
!	         dcape_3d(i,j)=7.8e-6*xcapet(i,j)
!	       end if
!ljd
             end do
           end do

!---------------------------------------------------------------------
!    if desired, print out debug quantities.
!---------------------------------------------------------------------
           if (in_debug_window) then
             if (dcape_3d(itest,jdebug) .gt. 20.) then
               print *, 'DONNER_DEEP/donner_deep: dcape = ',dcape_3d(itest,jdebug)
  	       print *, 'DONNER_DEEP/donner_deep: xcape(-1)= ',xcape_3d(itest,jtest)
             endif
           endif
!           call toc ( 'cvct', '7')

!---------------------------------------------------------------------
!   call cupar_vect to calculate normalized deep convective forcing
!---------------------------------------------------------------------
!           call tic ( 'cvct', '9')
           call cupar_vect (xcape_vect, coin_vect, dcape_3d, dqls_vect,&
                            omint2_vect, pcape_vect, plfc_vect,   &
			    plzb_vect, prini_vect, rini_vect,    &
			    qlsd_vect, rcape_vect, rpc_vect, tini_vect,&
			    tcape_vect, tpca_vect, jabs)
!           call toc ( 'cvct', '9')
            
	    do ilon=iminp,imaxp
	      do jlat=jminp,jmaxp
	        jgl=jabs(jlat)
	        ampu=ampta1_3d(ilon,jlat)
	        emdi=emdi_v(ilon,jlat)
	        contot=contot_v(ilon,jlat)
	        tprei=tprea1_3d(ilon,jgl)
	        do k=1,nlev
	          prinp(k)=pfull(ilon,jlat,k)
	          rmuf(k)=umeml_3d(ilon,jlat,k)
	          trf(k)=temold(ilon,jlat,k)
	        end do
	        przm=prinp(nlev)
	        prztm=-10.
!	        if (ampu .ne. 0.) then
!                 write(6,*) 'pre prean i,j,jgl= ',ilon,jlat,jgl
!                 write(6,*) 'ampu,emdi,contot= ',ampu,emdi,contot
!                 write(6,*) 'tprei= ',tprei
!                 write(6,*) 'ampta1= ',ampta1_3d(ilon,jlat)
!                 write(6,*) 'emdi_v= ',emdi_v(ilon,jlat)
!                 write(6,*) 'contot_v= ',contot_v(ilon,jlat)
!                 write(6,*) 'tprea1= ',tprea1_3d(ilon,jgl)
!               end if
                call prean(ampu,contot,emdi,prinp,rmuf,trf,tprei,xice)
                do k=1,nlev
                  xice_3d(ilon,jlat,k)=xice(k)
                  if ((xice(k) .ge. 1.e-10) .and. (prztm .le. 0.)) then 
                    prztm=prinp(k)
                    cycle
                  end if
                end do
                do k=1,nlev
                  if ((xice(k) .lt. 1.e-11) .and. (prinp(k) .ge. &
                     prztm)) then
                    przm=prinp(k-1)
                    exit
                  end if
                end do

!               if (prztm .ge. 0.) then
!                 do k=1,nlev
!                   write(6,*) 'k,prinp,xice= ',k,prinp(k),   &
!		                xice_3d(ilon,jlat,k)
!                 end do
!               end if
                do k=1,nlev
                  prent=prinp(k)
                  dgeicer=0.
                  if (xice(k) .ge. 1.e-10) call andge(prent,przm,prztm,&
                     dgeicer)
                  dgeice_3d(ilon,jlat,k)=dgeicer
                end do
              end do
            end do

!
!        Code to link Donner convection to Tiedtke/Rotstayn cloud
!        fraction/microphysics
!
!        Calculate dmeso_3d, the detrainment rate from the convective system
!        dmeso_3d has units of s**-1 and is (1/rho) dM/dz, where M is the
!        detraining mass flux of the convective system (kg/(m**2 s))
!
         do ilon=iminp,imaxp
           do jlat=jminp,jmaxp
             do k=1,nlev
               xliq_3d(ilon,jlat,k)=0.
               if (xice_3d(ilon,jlat,k) .ge. 1.0e-10) then
                 dmeso_3d(ilon,jlat,k)=emes_3d(ilon,jlat,k)/ &
                                       xice_3d(ilon,jlat,k)
               else
                 dmeso_3d(ilon,jlat,k)=0.
               end if
             end do
             do k=1,nlev-1
               if ((uceml_3d(ilon,jlat,k) .le. 1.0e-10)  .and.   &
                   (umeml_3d(ilon,jlat,k) .le. 1.0e-10)) then
                 mhalf_3d(ilon,jlat,k)=0.
                 mhalf_3d(ilon,jlat,k+1)=0.
                 cycle
               end if
               dnna=phalf(ilon,jlat,k+1)-pfull(ilon,jlat,k)
               dnnb=pfull(ilon,jlat,k+1)-phalf(ilon,jlat,k+1)
               mhalf_3d(ilon,jlat,k+1)=(dnnb*uceml_3d(ilon,jlat,k) + &
                                      dnna*uceml_3d(ilon,jlat,k+1)   &
                                      +dnnb*umeml_3d(ilon,jlat,k)    &
                              +dnna*umeml_3d(ilon,jlat,k+1))/(dnna+dnnb)
             end do
             mhalf_3d(ilon,jlat,nlev+1)=0.
             mhalf_3d(ilon,jlat,1)=0.
           end do
         end do
      call STRAT_CLOUD_DONNER_TEND(dmeso_3d,xliq_3d,xice_3d,mhalf_3d, &
              phalf,qlin,qiin,cf,qltend,qitend,qatend)

!
!     Convert tendencies to increments.
!
       qltend(:,:,:)=qltend(:,:,:)*dt
       qitend(:,:,:)=qitend(:,:,:)*dt
       qatend(:,:,:)=qatend(:,:,:)*dt


!---------------------------------------------------------------------
!   END OF STANDARD TIMESTEP CODE
!---------------------------------------------------------------------
         endif ! standard_step

!---------------------------------------------------------------------
!   END OF DONNER_DEEP CALCULATION STEP CODE
!---------------------------------------------------------------------
       endif ! ldeep

!---------------------------------------------------------------------
!   STANDARD TIMESTEP CODE
!   execute the following code on "standard" timesteps.
!---------------------------------------------------------------------
       if (standard_step) then
!----------------------------------------------------------------------
!    if this is is a standard timestep, add the
!    previously obtained temperature and moisture tendencies from deep
!    convection to those fields. if the moisture tendency results in
!    the production of a negative value, reduce the tendency to produce
!    a zero value for moisture, and save the modified tendency for 
!    analysis purposes.
!----------------------------------------------------------------------
!         call tic ( 'cvct', '10')
         do k=1,nlev
           do j=jminp,jmaxp
             do i=iminp,imaxp
!              temold(i,j,k) = temold(i,j,k) + cemetf(i,jabs(j),k)*dt
               ttnd  (i,j,k) = cemetf(i,jabs(j),k)*dt
               rsave(i,j,k)  = ratold(i,j,k)
!              ratold(i,j,k) = ratold(i,j,k) + cememf(i,jabs(j),k)*dt
               qtnd  (i,j,k) = cememf(i,jabs(j),k)*dt
               if ((ratold(i,j,k) + qtnd(i,j,k)) .lt. 0.) then
                 if (rsave(i,j,k) .gt. 0.) then
!           ratold(i,j,k) = 0.
	           qtnd  (i,j,k) = -rsave(i,j,k)/dt
                   cememf_mod(i,j,k) = -rsave(i,j,k)/dt
                 else 
                   cememf_mod(i,j,k) = 0.
                   qtnd(i,j,k) = 0.0
!           ratold(i,j,k) = rsave(i,j,k)
                 end if
               else
	         cememf_mod(i,j,k) = cememf(i,jabs(j),k)
               end if
             end do
           end do
         end do

!----------------------------------------------------------------------
!  define the rainfall and snowfall accrued on the current timestep
!  from deep convection. note: tprea1_3d [mm/day] * 1/86400 [day/sec] *
!  1/1000 [ m/mm] * 1000 [kg(h2o)/m**3] * dt [sec] = kg/m2, as desired. 
!----------------------------------------------------------------------
       do j=jminp,jmaxp
	 do i=iminp,imaxp
	   if (coldT(i,j)) then
	     snow(i,j) = tprea1_3d(i,jabs(j))*dt/86400.
	     rain(i,j) = 0.0
           else
	     rain(i,j) = tprea1_3d(i,jabs(j))*dt/86400.
	     snow(i,j) = 0.0
	   endif
         end do
       end do

!---------------------------------------------------------------------
!    if desired, print out debug quantities.
!---------------------------------------------------------------------
	 if (in_debug_window) then
           do k=ktest_model,nlev
  	     if (abs(cemetf(itest,jtest,k     )) .gt. .002) then
  	       print *, 'DONNER_DEEP/donner_deep: k, cemetf= ',k,   &
		 			     cemetf(itest,jtest,k)
      	       print *, 'DONNER_DEEP/donner_deep: k, cual= ',k,    &
					   cual_3d(itest,jdebug,k )
  	       print *, 'DONNER_DEEP/donner_deep: k, xcape= ',k,    &
					       xcape_3d(itest,jtest) 
  	       print *, 'DONNER_DEEP/donner_deep: k, dcape = ',k,    &
					  dcape_3d(itest,jdebug)
     	       print *, 'DONNER_DEEP/donner_deep: k,a1    = ',k,    &
					    a1_3d (itest,jdebug)
               print *, 'DONNER_DEEP/donner_deep: k, amax  = ',k,   &
					  amax_3d(itest,jdebug)
             end if
	   end do
         endif
       end if      ! (standard_step)

!---------------------------------------------------------------------
!   END OF STANDARD TIMESTEP CODE
!---------------------------------------------------------------------
!       call toc ( 'cvct', '10')

!---------------------------------------------------------------------
!   DONNER_DEEP CALCULATION ON NEXT TIMESTEP CODE
!   execute this code when donner_deep is to be called on the next timestep.
!---------------------------------------------------------------------
       if (ldeep_next) then

!--------------------------------------------------------------------
!   call precu to prepare single-column arrays for Donner (1993, JAS)
!   deep cumulus parameterization. note that the fields have been
!   updated from their entry into this subroutine by the current 
!   tendencies.
!--------------------------------------------------------------------
!         call tic ( 'cvct', '11')
!        call precu_vect (temold, ratold, press, pcape_vect,   &
!        call precu_vect (temold+ttnd, ratold+qtnd, press, pcape_vect, &
         call precu_vect (temold+ttnd, ratold+qtnd, pfull, pcape_vect, &
                          prini_vect, jabs, rcape_vect, rini_vect,    &
		          tcape_vect,  tini_vect)
!         call toc ( 'cvct', '11')

!--------------------------------------------------------------------
!   call cape to evaluate the column soundings for Donner (1993, JAS)
!   deep cumulus parameterization
!--------------------------------------------------------------------
!         call tic ( 'cvct', '12')
         call cape_vect (pcape_vect, rcape_vect, tcape_vect,   &
			 coin_vect, plcl_vect, plfc_vect, plzb_vect,  &
			 rpc_vect, tpca_vect, xcape_vect, jabs)
!         call toc ( 'cvct', '12')

!---------------------------------------------------------------------
!   save cape returned above into a module variable for use on next
!   timestep. calculate integrated vapor; also save it as a module 
!   variable for use on the next timestep.
!--------------------------------------------------------------------
!         call tic ( 'cvct', '13')
         do j=jminp,jmaxp
           do i=iminp,imaxp
             xcape_3d(i,jabs(j)) = xcape_vect(i,j)

             qint_vect(i,j) = rcape_vect(i,j,1)*(pcape_vect(i,j,1) -  &
			      pcape_vect(i,j,2))
             do k=2,ncap-1
               qint_vect(i,j) = qint_vect(i,j) + rcape_vect(i,j,k)*  &
			        0.5*(pcape_vect(i,j,k-1) -  &
	                             pcape_vect(i,j,k+1))
             end do
             qint_vect(i,j) = qint_vect(i,j) + rcape_vect(i,j,ncap)*  &
			      (pcape_vect(i,j,ncap-1) -  &
		               pcape_vect(i,j,ncap))
             qint_vect(i,j) = qint_vect(i,j)/gravm
             qint_3d(i,jabs(j)) = qint_vect(i,j)
           end do
         end do
!         call toc ( 'cvct', '13')
       end if  ! ldeep_next
 
!---------------------------------------------------------------------
!   END OF DONNER_DEEP CALCULATION ON NEXT TIMESTEP CODE
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if desired, print out debug quantities.
!---------------------------------------------------------------------
!       call tic ( 'cvct', '14')
       if (in_debug_window) then
         print *, 'DONNER_DEEP/donner_deep: plcl ', plcl_3d(itest,jdebug)
         print *, 'DONNER_DEEP/donner_deep: plfc ', plfc_3d(itest,jdebug)
         print *, 'DONNER_DEEP/donner_deep: plzb ', plzb_3d(itest,jdebug)
         print *, 'DONNER_DEEP/donner_deep: xcape ', xcape_3d(itest,jtest)
         print *, 'DONNER_DEEP/donner_deep: coin ', coin_3d(itest,jdebug)
         print *, 'DONNER_DEEP/donner_deep: dcape ', dcape_3d(itest,jdebug)
         print *, 'DONNER_DEEP/donner_deep: qint ', qint_3d(itest,jtest)
         print *, 'DONNER_DEEP/donner_deep: a1   ', a1_3d  (itest,jdebug)
         print *, 'DONNER_DEEP/donner_deep: amax ', amax_3d(itest,jdebug)
         print *, 'DONNER_DEEP/donner_deep: amos ', amos_3d(itest,jdebug)
         print *, 'DONNER_DEEP/donner_deep: tprea1 ', tprea1_3d(itest,jtest)
         print *, 'DONNER_DEEP/donner_deep: ampta1 ', ampta1_3d(itest,jdebug)
         print *, 'DONNER_DEEP/donner_deep: omint', omint_3d(itest,jtest)
         print *, 'DONNER_DEEP/donner_deep: rcoa1 ', rcoa1_3d(itest,jdebug)

         do k=ktest_model,nlev
           print *, 'DONNER_DEEP/donner_deep: k = ', k
           print *, 'DONNER_DEEP/donner_deep: cemetf', cemetf (itest,jtest,k)
           print *, 'DONNER_DEEP/donner_deep: ceefc ', ceefc  (itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: cecon ', cecon  (itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: cemfc ', cemfc  (itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: cememf', cememf (itest,jtest,k)
           print *, 'DONNER_DEEP/donner_deep: cememf_mod',    &
					  cememf_mod (itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: cual  ', cual_3d(itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: fre   ', fre_3d (itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: elt   ', elt_3d (itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: cmus  ', cmus_3d(itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep ', ecds_3d(itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: eces  ', eces_3d(itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: emds  ', emds_3d(itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: emes  ', emes_3d(itest,jdebug,k)
           print *, 'DONNER_DEEP/donner_deep: qmes  ', qmes_3d(itest,jdebug,k)
	   print *, 'DONNER_DEEP/donner_deep: wmps  ', wmps_3d(itest,jdebug,k)
	   print *, 'DONNER_DEEP/donner_deep: wmms  ', wmms_3d(itest,jdebug,k)
	   print *, 'DONNER_DEEP/donner_deep: tmes  ', tmes_3d(itest,jdebug,k)
	   print *, 'DONNER_DEEP/donner_deep: dmeml ', dmeml_3d(itest,jdebug,k)
	   print *, 'DONNER_DEEP/donner_deep: uceml ', uceml_3d(itest,jdebug,k)
	   print *, 'DONNER_DEEP/donner_deep: umeml ', umeml_3d(itest,jdebug,k)
	 end do
       endif

!---------------------------------------------------------------------
!  define a flag array to indicate those points where donner deep
!  convection has occurred.
!--------------------------------------------------------------------
       do k=1,nlev
         do j=jminp,jmaxp
           do i=iminp,imaxp
             if (tprea1_3d(i,jabs(j)) .ne. 0. .and.  &
                 cemetf(i,jabs(j),k) .ne. 0.) then   
	       iflag(i,j,k) = 5
             else
	       iflag(i,j,k) = 0
             end if
           end do
         end do
       end do


!---------------------------------------------------------------------
!  if running in FMS, determine if all points in this subdomain have 
!  been integrated on this timestep. if so, set the point counter to
!  zero and define the next time at which this module is to be executed.
!----------------------------------------------------------------------
       if (present (Time) ) then
         if (ldeep) then
           num_pts2 = num_pts2 + imaxp*jmaxp
           if (num_pts2 >= total_pts) then
             Next_donner_deep_time = Next_donner_deep_time +    &
				     Donner_deep_timestep
             num_pts2 = 0 
           endif

!---------------------------------------------------------------------
!  if running in FMS and saving diagnostics and this is a calculation
!  time step, send the diagnostic data to the diagnostics handler 
!  module.
!---------------------------------------------------------------------
           if (save_donner_deep_diagnostics) then
	     used = send_data (id_cemetf_deep,    &
			       cemetf(:,jabs(jminp):jabs(jmaxp),:),   &
			       Time, is, js, 1)
             used = send_data (id_ceefc_deep,   ceefc(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_cecon_deep,  cecon(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_cemfc_deep, cemfc(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_cememf_deep,    &
           	               cememf(:,jabs(jminp):jabs(jmaxp),:),  &
			       Time, is, js, 1)
             used = send_data (id_cememf_mod_deep, cememf_mod(:,:,:), &
			       Time, is, js, 1)
             used = send_data (id_cual_deep, cual_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_fre_deep, fre_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_elt_deep, elt_3d(:,:,:), Time,&
    	                       is, js, 1)
             used = send_data (id_cmus_deep, cmus_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_ecds_deep, ecds_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_eces_deep, eces_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_emds_deep, emds_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_emes_deep, emes_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_qmes_deep, qmes_3d(:,:,:), Time,&
          	               is, js, 1)
             used = send_data (id_wmps_deep, wmps_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_wmms_deep, wmms_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_tmes_deep, tmes_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_dmeml_deep, dmeml_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_uceml_deep, uceml_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_umeml_deep, umeml_3d(:,:,:), Time,&
 	                       is, js, 1)
             used = send_data (id_xice_deep, xice_3d(:,:,:), Time,&
                               is, js, 1)
             used = send_data (id_dgeice_deep, dgeice_3d(:,:,:), Time,&
                               is, js, 1)
             used = send_data (id_cuqi_deep, cuqi_3d(:,:,:), Time,&
                               is, js, 1)
             used = send_data (id_cuql_deep, cuql_3d(:,:,:), Time,&
                               is, js, 1)

             used = send_data (id_plcl_deep, plcl_3d(:,:), Time,&
                	       is, js)
             used = send_data (id_plfc_deep, plfc_3d(:,:), Time,&
 	                       is, js)
             used = send_data (id_plzb_deep, plzb_3d(:,:), Time,&
 	                       is, js)
             used = send_data (id_xcape_deep,    &
                               xcape_3d(:,jabs(jminp):jabs(jmaxp)),  &
			       Time, is, js)
             used = send_data (id_coin_deep, coin_3d(:,:), Time,&
 	                       is, js)
             used = send_data (id_dcape_deep, dcape_3d(:,:), Time,&
 	                       is, js)
             used = send_data (id_qint_deep,     &
                               qint_3d(:,jabs(jminp):jabs(jmaxp)),   &
			       Time, is, js)
             used = send_data (id_a1_deep, a1_3d(:,:), Time,&
 	                       is, js)
             used = send_data (id_amax_deep, amax_3d(:,:), Time,&
 	                       is, js)
             used = send_data (id_amos_deep, amos_3d(:,:), Time,&
 	                       is, js)
             used = send_data (id_tprea1_deep,    &
                               tprea1_3d(:,jabs(jminp):jabs(jmaxp)),  &
			       Time, is, js)
             used = send_data (id_ampta1_deep, ampta1_3d(:,:), Time,&
 	                       is, js)
             used = send_data (id_omint_deep,    &
                               omint_3d(:,jabs(jminp):jabs(jmaxp)),  &
			       Time, is, js)
             used = send_data (id_rcoa1_deep, rcoa1_3d(:,:), Time,&
 	                       is, js)
           endif
	 endif  ! ldeep

!---------------------------------------------------------------------
!     if running in Skyhi and if archived diagnostics from this module 
!     are desired on this timestep, call the routine to produce them.
!---------------------------------------------------------------------
       else  !  (present(Time)
         if (save_donner_deep_diagnostics .and. ldeep) then
           call write_diagnostics (jabs(1))
         endif
       endif  ! present(Time)

!---------------------------------------------------------------------
!    keep count of the number of points in the subdomain that have been 
!    initialized, so that initialization can be turned off once all 
!    points have been treated.
!---------------------------------------------------------------------
       if (do_init) then
         num_pts = num_pts + size(temold,1)*size(temold,2)
         if (num_pts == total_pts) do_init = .false.
       endif
!       call toc ( 'cvct', '14')

!--------------------------------------------------------------------


end subroutine donner_deep






!####################################################################

subroutine read_restart_donner_deep

!---------------------------------------------------------------------
!   this subroutine reads a restart file previously written by this 
!   module.
!---------------------------------------------------------------------

      integer             :: next(2), dt(2), cal
      integer             :: unit, vers
      character(len=4)    :: chvers

!-------------------------------------------------------------------- 
!   open the restart file.
!--------------------------------------------------------------------- 
      unit = open_file (file='INPUT/donner_deep.res',  &
!                       form='native', action='read')
                        form='unformatted', action='read')

!--------------------------------------------------------------------- 
!   read and check restart version number
!-------------------------------------------------------------------- 
      read (unit) vers
      if ( .not. any(vers == restart_versions) ) then
        write (chvers,'(i4)') vers
        call error_mesg ('read_restart_donner_deep', &
            'restart version '//chvers//' cannot be read &
            &by this module version', FATAL)
      endif

!--------------------------------------------------------------------
!   read time step and next occurrence from restart file. override 
!   provisional values, unless calendar type has changed.
!---------------------------------------------------------------------
      read (unit) next, dt, cal
      Old_time_step = set_time (dt(1), dt(2))
      if (cal == get_calendar_type() ) then
        Next_donner_deep_time = set_time( next(1), next(2) ) 	
      else
	if (get_my_pe() == 0) then
          call error_mesg ('read_restart_donner_deep',  &
	        'current calendar not same as restart calendar', NOTE)
	endif
      endif

!---------------------------------------------------------------------
!   read the restart data fields.
!---------------------------------------------------------------------
      call read_data (unit, cemetf)
      call read_data (unit, cememf)
      call read_data (unit, xcape_3d)
      call read_data (unit, qint_3d )
      call read_data (unit, omint_3d)
      if (vers /= 1) then
        call read_data (unit, tprea1_3d)
      else
        tprea1_3d = 0.0
      endif

!-------------------------------------------------------------------- 
!   close the restart file.
!--------------------------------------------------------------------- 
      call close_file (unit)

!--------------------------------------------------------------------



end subroutine read_restart_donner_deep




!##################################################################

subroutine read_usrtr_file

!---------------------------------------------------------------------
!    this subroutine reads an old usrtr restart file written in 
!    skyhi and extracts the donner_deep variables, placing them
!    in the proper storage locations.
!--------------------------------------------------------------------

      logical            :: lskp
      character(len=8)   :: string
      character(len=16)  :: string2
      integer            :: iskp, nxblck, ktbb, jjnn, jread, i, iounit
      real               :: skip

      real, dimension(:,:,:,:), allocatable  :: usmid
      real, dimension(:,:,:)  , allocatable  :: dymid

!---------------------------------------------------------------------
!    read usrtr file created as in SKYHI runs (these are needed only
!    as transition from SKYHI to new f90 donner_deep module structure -
!    new donner_deep.res file will be written by this module.
!---------------------------------------------------------------------
      iounit = open_file ('INPUT/usrtr_rstrt', form= 'unformatted', &
	 	           action='read')

!---------------------------------------------------------------------
!    read over file control block(s)
!---------------------------------------------------------------------
      read (iounit) iskp    , string  , string2  ,  iskp    , iskp  , &
                    string  , iskp    , iskp     ,  iskp    , iskp  , &
                    iskp    , iskp    , string   ,  string  , iskp  , &
                    iskp    , iskp    , iskp     ,  skip    , skip  , &
                    skip    , nxblck
      if (nxblck /= 0) then
	read (iounit)
      endif

!---------------------------------------------------------------------
!   allocate arrray to read tracer data into
!---------------------------------------------------------------------
      allocate ( usmid (idf, 1, nlev, ntracers_in_usrtr_file) ) 

!---------------------------------------------------------------------
!   read file and assign data to proper module variable
!---------------------------------------------------------------------
      do jread=1,jdf
	read (iounit)  ktbb, jjnn
	read (iounit) usmid
        cemetf     (:,jjnn,:  ) = usmid(:,1,:,24)
!!! NOTE ERROR
!       cememf     (:,jjnn,:  ) = usmid(:,1,:,24)
        cememf     (:,jjnn,:  ) = usmid(:,1,:,28)

	read (iounit) usmid

      end do
      deallocate (usmid)

      call close_file (iounit)

!---------------------------------------------------------------------
!    read dynamics file created as in SKYHI runs (these are needed only
!    as transition from SKYHI to new f90 donner_deep module
!    structure -- new donner_deep.res file will be written by
!    this code.
!---------------------------------------------------------------------
      iounit = open_file ('INPUT/dyn_rstrt', form= 'unformatted', &
	 	           action='read')

!---------------------------------------------------------------------
!   read over file control block(s)
!---------------------------------------------------------------------
      read (iounit) iskp    , string  , string2  ,  iskp    , iskp  , &
                    string  , iskp    , iskp     ,  iskp    , iskp  , &
                    iskp    , iskp    , string   ,  string  , iskp  , &
                    iskp    , iskp    , iskp     ,  skip    , skip  , &
                    skip    , nxblck
      if (nxblck /= 0) then
	read (iounit)
      endif

!---------------------------------------------------------------------
!   allocate arrray to read dynamics data into
!---------------------------------------------------------------------
      allocate ( dymid (idf, 1, 65+1+6*nlev                 ) ) 

!---------------------------------------------------------------------
!   read file and assign data to proper module variable
!---------------------------------------------------------------------
      do jread=1,jdf
	read (iounit)  ktbb, jjnn
	read (iounit) dymid

	xcape_3d(:,jjnn) = dymid(:,1,20+23)
	qint_3d(:,jjnn) = dymid(:,1,20+26)
	omint_3d(:,jjnn) = dymid(:,1,20+32)
	tprea1_3d(:,jjnn) = dymid(:,1,20+30)

	read (iounit) dymid

      end do
      deallocate (dymid)

      call close_file (iounit)

!---------------------------------------------------------------------


end subroutine read_usrtr_file



!###################################################################

subroutine initialize_variables

!--------------------------------------------------------------------
!      this subroutine initializes the subdomain-dimensioned module 
!      variables when the donner_deep scheme is being coldstarted.
!--------------------------------------------------------------------
     
      cemetf       = 0.0
      cememf       = 0.0

      xcape_3d     = 0.0
      qint_3d      = 0.0
      omint_3d     = 0.0
      tprea1_3d    = 0.0

end subroutine initialize_variables


!###################################################################

subroutine allocate_variables

!--------------------------------------------------------------------
!   this subroutine allocates space for the subdomain-dimensioned 
!   module variables.
!--------------------------------------------------------------------

     allocate ( cemetf         (idf, jdf, nlev ) )
     allocate ( cememf         (idf, jdf, nlev ) )

     allocate ( xcape_3d       (idf, jdf ) )
     allocate ( qint_3d        (idf, jdf ) )
     allocate ( omint_3d       (idf, jdf ) )
     allocate ( tprea1_3d      (idf, jdf ) )

!-------------------------------------------------------------------


end subroutine allocate_variables


!###################################################################

subroutine allocate_local_variables

!--------------------------------------------------------------------
!   this subroutine allocates space for the window-dimensioned module 
!   variables.
!--------------------------------------------------------------------

     allocate ( ceefc          (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( cecon          (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( cemfc          (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( cememf_mod     (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( cual_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( fre_3d         (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( elt_3d         (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( cmus_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( ecds_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( eces_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( emds_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( emes_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( qmes_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( wmps_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( wmms_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( tmes_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( dmeml_3d       (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( uceml_3d       (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( umeml_3d       (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( xice_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( dgeice_3d      (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( cuqi_3d        (iminp:imaxp, jminp:jmaxp, nlev) )
     allocate ( cuql_3d        (iminp:imaxp, jminp:jmaxp, nlev) )

     allocate ( plcl_3d        (iminp:imaxp, jminp:jmaxp ) )
     allocate ( plfc_3d        (iminp:imaxp, jminp:jmaxp ) )
     allocate ( plzb_3d        (iminp:imaxp, jminp:jmaxp ) )
     allocate ( coin_3d        (iminp:imaxp, jminp:jmaxp ) )
     allocate ( dcape_3d       (iminp:imaxp, jminp:jmaxp ) )
     allocate ( a1_3d          (iminp:imaxp, jminp:jmaxp ) )
     allocate ( amax_3d        (iminp:imaxp, jminp:jmaxp ) )
     allocate ( amos_3d        (iminp:imaxp, jminp:jmaxp ) )
     allocate ( ampta1_3d      (iminp:imaxp, jminp:jmaxp ) )
     allocate ( rcoa1_3d       (iminp:imaxp, jminp:jmaxp ) )
     allocate ( contot_v       (iminp:imaxp, jminp:jmaxp ) )
     allocate ( emdi_v         (iminp:imaxp, jminp:jmaxp ) )

!-------------------------------------------------------------------


end subroutine allocate_local_variables



!####################################################################

subroutine initialize_local_variables

!--------------------------------------------------------------------
!   this subroutine initializes the window-dimensioned module variables.
!--------------------------------------------------------------------

     ceefc      = 0.0
     cecon      = 0.0
     cemfc      = 0.0
     cual_3d    = 0.0
     fre_3d     = 0.0
     elt_3d     = 0.0
     cmus_3d    = 0.0
     ecds_3d    = 0.0
     eces_3d    = 0.0
     emds_3d    = 0.0
     emes_3d    = 0.0
     qmes_3d    = 0.0
     wmps_3d    = 0.0
     wmms_3d    = 0.0
     tmes_3d    = 0.0
     dmeml_3d   = 0.0
     uceml_3d   = 0.0
     umeml_3d   = 0.0
     xice_3d    = 0.0
     dgeice_3d  = 0.0
     cuqi_3d    = 0.0
     cuql_3d    = 0.0

     plcl_3d    = 0.0
     plfc_3d    = 0.0
     plzb_3d    = 0.0
     coin_3d    = 0.0
     dcape_3d   = 0.0
     a1_3d      = 0.0
     amax_3d    = 0.0
     amos_3d    = 0.0
     ampta1_3d  = 0.0
     rcoa1_3d   = 0.0

!-------------------------------------------------------------------


end subroutine initialize_local_variables



!####################################################################

subroutine fill_donner_c_variables (j, oom, spare)

!--------------------------------------------------------------------
!   collects the desired module variables for collective processing
!   THIS SUBROUTINE IS USED ONLY IN SKYHI
!---------------------------------------------------------------------

integer,                  intent(in)              :: j
real, dimension(:,:,:,:), intent(out), optional   :: oom
real, dimension(:,:,:)  , intent(out), optional   :: spare

!--------------------------------------------------------------------
!   intent(in) variables:
! 
!      j                  subdomain starting index of the segment
!                         being processed
!
!   intent(out), optional variables:
!
!      oom                collected 3D module diagnostic variables
!      spare              collected 2D module diagnostci variables
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!   if the 3D diagnostic variables are desired, collect them.
!-------------------------------------------------------------------
    if (present (oom)) then
      oom(:,jminp:jmaxp,:, 1) =  cemetf(:,j+jminp-1:j+jmaxp-1,:) 
      oom(:,:,:, 2) =  ceefc          (:,:,:  ) 
      oom(:,:,:, 3) =  cecon          (:,:,:  ) 
      oom(:,:,:, 4) =  cemfc          (:,:,:  ) 
      oom(:,:,:, 5) =  cememf_mod     (:,:,:  ) 
      oom(:,:,:, 6) =  cual_3d        (:,:,:  ) 
      oom(:,:,:, 7) =  fre_3d         (:,:,:  ) 
      oom(:,:,:, 8) =  elt_3d         (:,:,:  ) 
      oom(:,:,:, 9) =  cmus_3d        (:,:,:  ) 
      oom(:,:,:,10) =  ecds_3d        (:,:,:  ) 
      oom(:,:,:,11) =  eces_3d        (:,:,:  ) 
      oom(:,:,:,12) =  emds_3d        (:,:,:  ) 
      oom(:,:,:,13) =  emes_3d        (:,:,:  ) 
      oom(:,:,:,14) =  qmes_3d        (:,:,:  ) 
      oom(:,:,:,15) =  wmps_3d        (:,:,:  ) 
      oom(:,:,:,16) =  wmms_3d        (:,:,:  ) 
      oom(:,:,:,17) =  tmes_3d        (:,:,:  ) 
      oom(:,:,:,18) =  dmeml_3d       (:,:,:  ) 
      oom(:,:,:,19) =  uceml_3d       (:,:,:  ) 
      oom(:,:,:,20) =  umeml_3d       (:,:,:  ) 
     endif

!-------------------------------------------------------------------
!   if the 2D diagnostic variables are desired, collect them.
!-------------------------------------------------------------------
     if (present (spare) ) then
       spare (:,:, 1) = plcl_3d          (:,:)
       spare (:,:, 2) = plfc_3d          (:,:)
       spare (:,:, 3) = plzb_3d          (:,:)
       spare (:,jminp:jmaxp, 4) = xcape_3d(:,j+jminp-1:j+jmaxp-1)
       spare (:,:, 5) = coin_3d          (:,:)
       spare (:,:, 6) = dcape_3d          (:,:)
       spare (:,jminp:jmaxp, 7) = qint_3d (:,j+jminp-1:j+jmaxp-1)
       spare (:,:, 8) = a1_3d            (:,:)
       spare (:,:, 9) = amax_3d            (:,:)
       spare (:,:,10) = amos_3d            (:,:)
       spare (:,:,11) = tprea1_3d            (:,j+jminp-1:j+jmaxp-1)
       spare (:,:,12) = ampta1_3d            (:,:)
       spare (:,jminp:jmaxp,13) = omint_3d  (:,j+jminp-1:j+jmaxp-1)
       spare (:,:,14) = rcoa1_3d            (:,:)
     endif

!--------------------------------------------------------------------



end subroutine fill_donner_c_variables 




!###################################################################

subroutine store_donner_deep_variables (j, sparep, oom)

!--------------------------------------------------------------------
!   collects the diagnostic module variables at end of timestep to be 
!   stored on disk so that they may be retrieved on the next timestep.
!   needed with slab models when donner_deep not always called on every 
!   step. THIS SUBROUTINE IS USED ONLY IN SKYHI.
!---------------------------------------------------------------------

integer,                  intent(in)              :: j
real, dimension(:,:,:,:), intent(out)             :: oom
real, dimension(:,:,:)  , intent(out)             :: sparep

!--------------------------------------------------------------------
!   intent(in) variables:
! 
!      j                  subdomain starting index of the segment
!                         being processed
!
!   intent(out) variables:
!
!      oom                collected 3D module diagnostic variables
!      sparep             collected 2D module diagnostci variables
!
!-------------------------------------------------------------------

      oom(:,jminp:jmaxp,:, 1) =  cemetf     (:,j+jminp-1:j+jmaxp-1,:  ) 
      oom(:,:,:, 2) =  ceefc          (:,:,:  ) 
      oom(:,:,:, 3) =  cecon          (:,:,:  ) 
      oom(:,:,:, 4) =  cemfc          (:,:,:  ) 
      oom(:,:,:, 5) =  cememf         (:,j+jminp-1:j+jmaxp-1,:  ) 
      oom(:,:,:, 6) =  cual_3d           (:,:,:  ) 
      oom(:,:,:, 7) =  fre_3d         (:,:,:  ) 
      oom(:,:,:, 8) =  elt_3d         (:,:,:  ) 
      oom(:,:,:, 9) =  cmus_3d        (:,:,:  ) 
      oom(:,:,:,10) =  ecds_3d        (:,:,:  ) 
      oom(:,:,:,11) =  eces_3d        (:,:,:  ) 
      oom(:,:,:,12) =  emds_3d        (:,:,:  ) 
      oom(:,:,:,13) =  emes_3d        (:,:,:  ) 
      oom(:,:,:,14) =  qmes_3d        (:,:,:  ) 
      oom(:,:,:,15) =  wmps_3d        (:,:,:  ) 
      oom(:,:,:,16) =  wmms_3d        (:,:,:  ) 
      oom(:,:,:,17) =  tmes_3d        (:,:,:  ) 
      oom(:,:,:,18) =  dmeml_3d       (:,:,:  ) 
      oom(:,:,:,19) =  uceml_3d       (:,:,:  ) 
      oom(:,:,:,20) =  umeml_3d       (:,:,:  ) 

      sparep(:,:,20) = plcl_3d          (:,:)
      sparep(:,:,21) = plfc_3d          (:,:)
      sparep(:,:,22) = plzb_3d          (:,:)
      sparep(:,jminp:jmaxp,23) = xcape_3d(:,j+jminp-1:j+jmaxp-1)
      sparep(:,:,24) = coin_3d          (:,:)
      sparep(:,:,25) = dcape_3d          (:,:)
      sparep(:,jminp:jmaxp,26) = qint_3d (:,j+jminp-1:j+jmaxp-1)
      sparep(:,:,27) = a1_3d            (:,:)
      sparep(:,:,28) = amax_3d            (:,:)
      sparep(:,:,29) = amos_3d            (:,:)
      sparep(:,jminp:jmaxp,30) = tprea1_3d (:,j+jminp-1:j+jmaxp-1)
      sparep(:,:,31) = ampta1_3d            (:,:)
      sparep(:,jminp:jmaxp,32) = omint_3d  (:,j+jminp-1:j+jmaxp-1)
      sparep(:,:,33) = rcoa1_3d            (:,:)


!--------------------------------------------------------------------



end subroutine store_donner_deep_variables 



!###################################################################

subroutine fill_donner_deep_variables (j, sparep, oom)

!--------------------------------------------------------------------
!   returns the collected module variables to the module at start of
!   timestep. needed with slab models when donner_deep not always 
!   called on every step. THIS SUBROUTINE IS USED ONLY IN SKYHI.
!---------------------------------------------------------------------

integer,                  intent(in)              :: j
real, dimension(:,:,:,:), intent(out)             :: oom
real, dimension(:,:,:)  , intent(out)             :: sparep

!--------------------------------------------------------------------
!   intent(in) variables:
! 
!      j                  global starting index of the domain segment
!                         being processed
!
!   intent(out) variables:
!
!      oom                collected 3D module diagnostic variables
!      spare              collected 2D module diagnostci variables
!
!-------------------------------------------------------------------

      cemetf     (:,j+jminp-1:j+jmaxp-1,:  )  = oom(:,jminp:jmaxp,:,1)
      ceefc          (:,:,:  )   =   oom(:,:,:,2)
      cecon          (:,:,:  )   = oom(:,:,:,3) 
      cemfc          (:,:,:  )  = oom(:,:,:,4)
      cememf     (:,j+jminp-1:j+jmaxp-1,:  )  = oom(:,jminp:jmaxp,:,5)
      cual_3d        (:,:,:  )  = oom(:,:,:,6)
      fre_3d         (:,:,:  )  = oom(:,:,:,7)
      elt_3d         (:,:,:  )  = oom(:,:,:,8)
      cmus_3d        (:,:,:  )  = oom(:,:,:,9)
      ecds_3d        (:,:,:  )  = oom(:,:,:,10)
      eces_3d        (:,:,:  )  = oom(:,:,:,11)
      emds_3d        (:,:,:  )  = oom(:,:,:,12)
      emes_3d        (:,:,:  ) = oom(:,:,:,13)
      qmes_3d        (:,:,:  )   = oom(:,:,:,14)
      wmps_3d        (:,:,:  )  = oom(:,:,:,15)
      wmms_3d        (:,:,:  )  = oom(:,:,:,16)
      tmes_3d        (:,:,:  ) = oom(:,:,:,17)
      dmeml_3d       (:,:,:  ) = oom(:,:,:,18)
      uceml_3d       (:,:,:  ) = oom(:,:,:,19)
      umeml_3d       (:,:,:  ) = oom(:,:,:,20)

      plcl_3d          (:,:)  = sparep (:,:,1)
      plfc_3d          (:,:)  = sparep (:,:,2)
      plzb_3d          (:,:)  = sparep (:,:,3)
      xcape_3d(:,j+jminp-1:j+jmaxp-1) = sparep (:,jminp:jmaxp, 4)
      coin_3d          (:,:) = sparep (:,:,5)
      dcape_3d          (:,:) = sparep (:,:,6)
      qint_3d (:,j+jminp-1:j+jmaxp-1) = sparep(:,jminp:jmaxp, 7)
      a1_3d            (:,:)  = sparep(:,:,8)
      amax_3d            (:,:) = sparep(:,:,9)
      amos_3d            (:,:) = sparep(:,:,10)
      tprea1_3d            (:,j+jminp-1:j+jmaxp-1) = sparep(:,:,11)
      ampta1_3d            (:,:) = sparep(:,:,12)
      omint_3d  (:,j+jminp-1:j+jmaxp-1) = sparep(:,jminp:jmaxp,13)
      rcoa1_3d            (:,:) = sparep(:,:,14)

!--------------------------------------------------------------------



end subroutine fill_donner_deep_variables 



!###################################################################


subroutine get_cemetf (j, oom)

!--------------------------------------------------------------------
!   sends the heating rate due to deep convection to another module.
!---------------------------------------------------------------------

integer,                  intent(in)    :: j
real, dimension(:,:,:),   intent(out)   :: oom

!--------------------------------------------------------------------
!   intent(in) variables:
! 
!      j                  subdomain starting index of the  segment
!                         being processed
!
!   intent(out) variables:
!
!      oom                convective heating rate at the requested
!                         grid points
!
!-------------------------------------------------------------------

      oom(:,jminp:jmaxp,:) =  cemetf(:,j+jminp-1:j+jmaxp-1,:) 

!-------------------------------------------------------------------



end subroutine get_cemetf



!###################################################################

subroutine get_tprea1_donner_deep (j, spare30)

!--------------------------------------------------------------------
!   sends the total precip due to deep convection to another module.
!   THIS SUBROUTINE IS ONLY USED IN SKYHI
!---------------------------------------------------------------------

integer,              intent(in)     :: j
real, dimension(:,:), intent(out)    :: spare30

!--------------------------------------------------------------------
!   intent(in) variables:
! 
!      j                  subdomain starting index of the segment
!                         being processed
!
!   intent(out) variables:
!
!      spare30            total precip from deep convection
!
!-------------------------------------------------------------------

       spare30(:,jminp:jmaxp) = tprea1_3d(:,j:j+jmaxp-1)

!-------------------------------------------------------------------



end subroutine get_tprea1_donner_deep


!###################################################################

subroutine prean(ampu,contot,emdi,prf,rmuf,trf,tprei,xice)

! Calculates ice content for mesoscale circulation.
!
!      Leo Donner
!      GFDL
!      17 May 2001
!
!   Input arguments
!
real,  intent(in)  ::  ampu  ! fractional area of mesoscale
                             !   circulation
real,  intent(in)  ::  contot ! ratio of convective to total
!   precipitation
real,  intent(in)  ::  emdi  ! vertical integral of
!    mesoscale-downdraft sublimation  (mm/d)
real, dimension(:), intent(in)  ::  prf ! pressure at full levels
!   (Pa) index 1 at top of model
real, dimension(:), intent(in)  ::  rmuf ! mesoscale updraft mass flux
!  (kg/[(m**2) s]) Index 1 at model top
real, dimension(:), intent(in)  ::  trf  ! temperature at full GCM levels
!   (K)  index 1 at model top
real,  intent(in)  ::  tprei ! total convective-system
!   precipitation (mm/d)
!
!   Input arguments as module variables:
!
!   nlev                       number of GCM levels
!   rair  (module variable)    gas constant for dry air (J/(kg K))
!
!   Output arguments
!
real, dimension(:), intent(out) :: xice ! anvil ice (kg(ice)/kg)
! index 1 at model top
!
!      local workspace
!
  integer :: k                    ! vertical index
  integer :: kou                  ! counter
  real    :: endiw
 real    :: rho                  ! anvil (ht. avg.) air density (kg/(m**3))
 real    :: rhol
  real    :: tprew
 real    :: xicet                ! anvil ice work variable (kg(ice)/kg)
!
!   Initialize
!
 xicet=0.
 do k=1,nlev
    xice(k)=0.
   end do
!
!    Return if no mesoscale anvil exists.
!
  if (ampu .le. 0.) return
   if (contot .ge. 1.) return
 if ((tprei .le. 0.) .and. (emdi .le. 0.)) return
!
!    Convert precipitation and sublimation integral to kg/[(m**2) s].
!
   tprew=tprei/86400.
   emdiw=emdi/86400.
!
!     Calculate average air density over vertical extent of anvil.
!
   rho=0.
   kou=0
   do k=1,nlev
    if (rmuf(k) .ne. 0.) then
         kou=kou+1
         rhol=prf(k)/(rair*trf(k))
         rho=rho+rhol
        cycle

 end if
 end do
 if (kou .eq. 0) return
   rho=rho/kou
!
!      Calculate mesoscale ice content by balancing fall at anvil base with
!      mesoscale precipitation and sublimation in mesoscale downdraft.
!
  xicet=(1.-contot)*tprew
   xicet=xicet+emdiw
   xicet=xicet/(3.29*ampu)
   xicet=xicet**.862
  xicet=xicet/rho
!
!      Assign anvil ice to layers with postive meso updraft mass flux
!
  do k=1,nlev
    if (rmuf(k) .le. 0.) xice(k)=0.
      if (rmuf(k) .gt. 0.) xice(k)=xicet
  end do
! ljd
!   do k=1,nlev
! if (xice(k) .gt. .01) then
! write(6,*) 'prean ampu,contot,emdi= ',ampu,contot,emdi
! write(6,*) 'rho,rhol,kou= ',rho,rhol,kou
! write(6,*) 'tprei,rair= ',tprei, rair
!    do kk=1,nlev
!    write(6,*) 'k,prf,xice= ',kk,prf(kk),xice(kk)
!    write(6,*) 'k,prf,trf= ',kk,prf(kk),trf(kk)
!    write(6,*) 'k,prf,rmuf= ',kk,prf(kk),rmuf(kk)
!    end do
!    stop
!    end if
! end do
! ljd


end subroutine prean

!###################################################################
subroutine andge(p,pzm,pztm,dgeicer)
!
!  test routine for calculating anvil d_ge
!
!----------------------------------------------------------------------
       implicit none
!      Input Arguments
!
      integer, parameter          :: nmclev=6        ! number of levels
! in anvil
      real, intent(in)            :: p               ! pressure
      real, intent(in)            :: pzm             ! anvil base pressure
     real, intent(in)            :: pztm            ! anvil top pressure
     real, intent(out)           :: dgeicer         ! generalized effective
!  size at pressure p (micrometers)
! p, pzm, and pztm must have identical units
!
!      Local Variables
!
     real, dimension(nmclev)     :: dgeice  ! generalized effective
!  size of hexagonal ice crystals, defined as in Fu (1996, J. Clim.)
!  values from Table 2 of McFarquhar et al. (1999, JGR) are averaged over
!  all grid boxes for which D_ge is defined for all altitudes between
!  9.9 and 13.2 km. index 1 at bottom of anvil
        real, dimension(nmclev)     :: relht   ! distance from anvil
!  base, normalized by total anvil thickness. from Table 2 of McFarquhar
!  et al. (1999, JGR) for grid boxes with data between 9.9 and 13.2 km.
!  index 1 at anvil bottom
!
        integer                     :: k               ! vertical index
       real                        :: znor            ! normalized distance

!  from anvil base
      data dgeice/38.5,30.72,28.28,25.62,24.8,13.3/
      data relht/0.,.3,.45,.64,.76,1./
!      write(6,*) 'p,pzm,pztm= ',p,pzm,pztm
        if (pzm .lt. pztm) then
          write(6,*) 'error in angde pzm,pztm= ',pzm,pztm
          stop
        end if
       if (pzm .eq. pztm) then
         znor=.5
         go to 12
       end if
       znor=(pzm-p)/(pzm-pztm)
12   continue
!      write(6,*) 'znor= ',znor
        do k=2,nmclev
!       write(6,*) 'relhts= ',relht(k-1),relht(k)
          if ((znor .ge. relht(k-1)) .and. (znor .le. relht(k))) then
             dgeicer=dgeice(k-1)+( (znor-relht(k-1))*(dgeice(k)- &
               dgeice(k-1))/(relht(k)-relht(k-1)) )
!       write(6,*) 'dgeicer= ',dgeicer
        go to 11
       end if
      end do
11   continue
!      do k=1,nmclev
!      write(6,*) 'k,dgeice,relht= ',k,dgeice(k),relht(k)
!      end do

 end subroutine andge

!###################################################################

subroutine cupar_vect (cape_v, cin_v,           dcape_v, dqls_v,   &
                                          omint_v, pcape_v, plfc_v,  &
                  plzb_v, pr_v, q_v, qlsd_v, r_v, rpc_v,     t_v, &
                  tcape_v, tpc_v,            jabs)

!---------------------------------------------------------------------
real, dimension(:,:), intent(in)   :: cape_v, cin_v, dcape_v, dqls_v, &
				      omint_v, plfc_v, plzb_v, qlsd_v

real, dimension(:,:,:), intent(in) :: pcape_v, pr_v, q_v, r_v, rpc_v, &
				      t_v, tcape_v, tpc_v

integer, dimension(:), intent(in)  :: jabs


!-----------------------------------------------------------------------
!      Cupar drives the parameterization for deep cumulus convection
!      based on Donner (1993, J. Atmos. Sci.).
!
!-----------------------------------------------------------------------
!     On Input:
!
!        cape     convective available potential energy (J/kg)
!        cin      convective inhibtion (J/kg)
!        cpd      specific heat of dry air at constant pressure (J/(kg K))
!        cpv      specific heat of water vapor [J/(kg K)]
!        dcape    local rate of CAPE change by all processes
!                 other than deep convection [J/(kg s)]
!        dqls     local rate of change in column-integrated vapor
!                 by all processes other than deep convection
!                 {kg(H2O)/[(m**2) s]}
!        epsilo   ratio of molecular weights of water vapor to dry air
!        gravm    gravity constant [m/(s**2)]
!        ilon     longitude index
!        jlat     latitude index
!        mcu      frequency (in time steps) of deep cumulus
!        omint    integrated low-level displacement (Pa)
!        pcape    pressure at Cape.F resolution (Pa)
!                 Index 1 at bottom of model.
!        plfc     pressure at level of free convection (Pa)
!        plzb     pressure at level of zero buoyancy (Pa)
!        pr       pressure at Skyhi vertical resolution (Pa)
!                 Index 1 nearest ground  
!        q        large-scale vapor mixing ratio at Skyhi vertical resolution
!                 [kg(h2O)/kg]
!                 Index 1 nearest ground 
!        qlsd     column-integrated vapor divided by timestep for cumulus
!                 parameterization {kg(H2O)/[(m**2) s]}
!        r        large-scale vapor mixing ratio at Cape.F resolution
!                 [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rpc      parcel vapor mixing ratio from Cape.F [kg(h2O)/kg]
!                 Index 1 at bottom of model.
!        rd       gas constant for dry air (J/(kg K))
!        rlat     latent heat of vaporization (J/kg)
!        rv       gas constant for water vapor (J/(kg K))
!        t        large-scale temperature at Skyhi vertical resolution (K)
!                 Index 1 nearest ground
!        tcape    large-scale temperature at Cape.F resolution (K)
!                 Index 1 at bottom of model.
!        tpc      parcel temperature from from Cape.F (K)
!                 Index 1 at bottom of model.
!
!     On Input as Parameters:
!
!        kmax     number of vertical levels at Skyhi resolution
!        kpar     number of cumulus sub-ensembles
!        ncap     number of vertical levels in Cape.F resolution
!

      real, dimension(kpar) :: arat
      real, dimension (size(t_v,1), size(t_v,2), kpar    ) :: arat_v

      real, dimension(size(t_v,3)) :: cmus, ecds, eces, emds, pr, q, &
				      t, cual, emes, disa, disb, disc, &
				      disd, dise, elt, fre, qmes,   &
				      tmes, wmms, wmps, sig, dmeml, &
				      uceml, umeml
      real, dimension (ncap) ::qli0, qli1, qr, qt, r, ri, rl, tc,   &
			       pcape, tcape, tpc, rpc
      real, dimension (size(t_v,1), size(t_v,2), size(t_v,3) ) ::   &
			       sig_v, cmus_v, cual_v, ecds_v, &
				   eces_v, emds_v, emes_v, disa_v, &
				   disb_v, disc_v, disd_v, dise_v, &
				   dmeml_v, elt_v, fre_v, qmes_v, &
				   tmes_v, uceml_v, umeml_v, wmms_v, &
				   wmps_v, qtest, a1_vk, cuq_v,cuql_v

      real, dimension (size(t_v,1), size(t_v,2), ncap)  ::     &
			       qli0_v, qli1_v, qr_v, qt_v, ri_v, &
			       rl_v

      real, dimension (size(t_v,1), size(t_v,2)        ) ::     &
	       ampt_v, sfcqf_v, sfcsf_v, amax_v, &
	       cmui_v,  tpre_v, emei_v, &
                            amos_v, rco_v, tfint, a1_v, pdeet1, &
			    pdeet2
     logical, dimension (size(t_v,1), size(t_v,2)) :: exit_flag,    &
				     nine_flag
			       
      logical, dimension (size(t_v,1), size(t_v,2), size(t_v,3) ) ::   &
	        exit_3d


!--------------------------------------------------------------------
!
!     test set to true for graphics
!     must also dimension sig for graphics
!
      logical test
      test=.false.
!
!     arat(i) is the ratio at cloud base of the fractional area
!     of ensemble i to ensemble 1.
!
!     GATE:
!
      data arat/1.,.26,.35,.32,.3,.54,.66/



!        do j=jminp,jmaxp
!  if (debug_jt(j)) then
	  if (in_debug_window) then
            print *, 'DONNER_DEEP/cupar: entering cupar with itest, jtest,&
	                          & kttest:', itest, jtest, kttest
          endif
!       end do

      do jlat=jminp,jmaxp
        do ilon=iminp,imaxp

          nine_flag(ilon,jlat) = .false.

!---------------------------------------------------------------------
!  check that criteria for deep convection are satisfifed.
!---------------------------------------------------------------------
          pdeet1(ilon,jlat) = plfc_v(ilon,jlat) - plzb_v(ilon,jlat)
          pdeet2(ilon,jlat) = plfc_v(ilon,jlat) - pr_v(ilon,jlat,1)
          if ((cape_v (ilon,jlat) <= 0.) .or.  &
	      (dcape_v(ilon,jlat) <= 0.) .or. & 
	      (pdeet1(ilon,jlat) < pdeep_cv )       .or. &
	      (pdeet2(ilon,jlat) < omint_v(ilon,jlat)) .or.  &
              (cin_v(ilon,jlat) > cdeep_cv) ) then
            exit_flag(ilon,jlat) = .true.
	  else
            exit_flag(ilon,jlat) = .false.
          endif
        end do
      end do

!	  do j=jminp,jmaxp
!            if (debug_jt(j)) then
	    if (in_debug_window) then
	      print *, 'DONNER_DEEP/cupar: omint,cape,dcape, exit_flag',    &
			omint_v(itest,jdebug), cape_v(itest,jdebug),   &
			dcape_v(itest,jdebug), exit_flag(itest,jdebug)
            endif
!          end do

      do jlat=jminp,jmaxp
        do ilon=iminp,imaxp

	 jgl = jabs(jlat) 



!-------------------------------------------------------------------
!   initialize fields.
!--------------------------------------------------------------------
!     do j=1,kmax
      do j=1,nlev
          cemetf(ilon,jgl ,j          )=0.
          ceefc (ilon,jlat,j          )=0.
	  cecon(ilon,jlat,j          )=0.
	  cemfc(ilon,jlat,j          )=0.
	  cememf(ilon,jgl ,j          )=0.
	  cememf_mod(ilon,jlat,j          )=0.
	  cual_3d(ilon,jlat,j          )=0.
	  fre_3d(ilon,jlat,j          )=0.


          elt_3d  (ilon,jlat,j) = 0.0
          cmus_3d (ilon,jlat,j) = 0.0
          ecds_3d (ilon,jlat,j) = 0.0
          eces_3d (ilon,jlat,j) = 0.0
          emds_3d (ilon,jlat,j) = 0.0
          emes_3d (ilon,jlat,j) = 0.0
          qmes_3d (ilon,jlat,j) = 0.0
          wmps_3d (ilon,jlat,j) = 0.0
          wmms_3d (ilon,jlat,j) = 0.0
          tmes_3d (ilon,jlat,j) = 0.0
          dmeml_3d(ilon,jlat,j) = 0.0
          uceml_3d(ilon,jlat,j) = 0.0
          umeml_3d(ilon,jlat,j) = 0.0
	  xice_3d(ilon,jlat,j) = 0.0
	  dgeice_3d(ilon,jlat,j) = 0.0
          cuqi_3d(ilon,jlat,j) = 0.0
          cuql_3d(ilon,jlat,j) = 0.0
          cuq_v(ilon,jlat,j) = 0.0
          cuql_v(ilon,jlat,j) = 0.0
          a1_3d(ilon,jlat) = 0.0
          amax_3d(ilon,jlat) = 0.0
          amos_3d(ilon,jlat) = 0.0
          tprea1_3d(ilon,jgl ) = 0.0
          ampta1_3d(ilon,jlat) = 0.0
          rcoa1_3d(ilon,jlat) = 0.0
      end do
      end do
      end do

!  do j=jminp,jmaxp
!	    if (debug_jt(j)) then
	    if (in_debug_window) then
              print *, 'DONNER_DEEP/cupar: cpi,cpv= ',cpi,cpv
              print *, 'DONNER_DEEP/cupar: rocp,rair,latvap= ',rocp,rair, &
			latvap
              print *, 'DONNER_DEEP/cupar: rvap= ',rvap
              print *, 'DONNER_DEEP/cupar: cape,cin,plzb= ',  &
			cape_v(itest,jdebug), &
			cin_v(itest,jdebug),plzb_v(itest,jdebug)
              print *, 'DONNER_DEEP/cupar: plfc= ',plfc_v(itest,jdebug)
              print *, 'DONNER_DEEP/cupar: dqls,qlsd= ',dqls_v(itest,jdebug),&
			qlsd_v(itest,jdebug)
              do k=1,nlev-ktest_model+1
                print *, 'DONNER_DEEP/cupar: k,pr,t,q= ',k,   &
			  pr_v(itest,jdebug,k),   &
			 t_v(itest,jdebug,k), q_v(itest,jdebug,k)
              end do
            endif
!          end do

      
!---------------------------------------------------------------------
!  call routine to calculate non-normalized cumulus forcings.
!---------------------------------------------------------------------
      sfcqf_v(:,:      )=0.
      sfcsf_v(:,:      )=0.
       do k=1,kpar
       arat_v(:,:,k) = arat(k)
       end do


       if (test) then
!        call opngks
       endif

       call mulsub_vect (ampt_v,arat_v,pr_v,q_v,sfcqf_v,sfcsf_v,t_v, &
			  amax_v,cmui_v,cmus_v, cual_v,cuq_v, cuql_v, &
			  ecds_v,eces_v,emds_v,emei_v,emes_v,disa_v, &
             disb_v, disc_v,disd_v,dise_v,dmeml_v, &
       elt_v,fre_v,qmes_v,tmes_v,tpre_v,uceml_v,umeml_v,wmms_v,wmps_v, &
       jabs, exit_flag)


!ljd
!         do jlat=jminp,jmaxp
!        do ilon=iminp,imaxp
!           if (ampt_v(ilon,jlat) .ne. 0.) then
!           write(6,*) 'cupar 1 ilon,jlat,jgl= ',ilon,jlat,jgl
!             write(6,*) 'cupar 1 contot= ',contot_v(ilon,jlat)
!           end if
!        end do
!        end do
!ljd


!	  do j=jminp,jmaxp
!	    if (debug_jt(j)) then
            if (in_debug_window) then
	  if (.not. exit_flag(itest,jdebug)) then
             print *,   'DONNER_DEEP/cupar:  ampt,tpre= ',  &
		      ampt_v(itest,jdebug),  tpre_v(itest,jdebug)
      	     do k=1,nlev-ktest_model+1    
 	       print *,  'DONNER_DEEP/cupar: k,cual= ',k,cual_v(itest,jdebug,k)
             end do
            endif
            endif
!         end do


      do jlat=jminp,jmaxp
        do ilon=iminp,imaxp
	  if (.not. exit_flag(ilon,jlat)) then
 	    amos_v(ilon,jlat)=0.
	    rco_v(ilon,jlat)=tpre_v(ilon,jlat)*contot_v(ilon,jlat)
            if (tpre_v(ilon,jlat) .eq. 0.) then
	      a1_v(ilon,jlat)=0.
	      amax_v(ilon,jlat)=0.
	      nine_flag(ilon,jlat) = .true.
	    end if
          endif
        end do
      end do




!          do j=jminp,jmaxp
!     if (debug_jt(j)) then
       if (in_debug_window) then
       if (exit_flag(itest,jdebug) ) then
 	print *, 'DONNER_DEEP/cupar: deep convection not present in&
		  &test column. exit_flag is set to .true.'
       else
 	print *, 'DONNER_DEEP/cupar: deep convection is possible&
		& -- exit_flag is .false.'
              endif
            endif
!          end do

!---------------------------------------------------------------------
!  interpolate forcings to resolution of cuclo.
!--------------------------------------------------------------------
      call polat_vect (dise_v, pr_v, pcape_v, qr_v) 
      call polat_vect (disa_v, pr_v, pcape_v, qt_v) 

!--------------------------------------------------------------------
!   initialize variables.
!--------------------------------------------------------------------
      do k=1,ncap
	do j=jminp,jmaxp
	do i=iminp,imaxp
	 qli0_v(i,j,k)=0.
	 qli1_v(i,j,k)=0.
	 rl_v(i,j,k)=0.
	 ri_v(i,j,k)=0.
      end do
      end do
      end do

!	  do j=jminp,jmaxp
!	    if (debug_jt(j)) then
            if (in_debug_window) then
          if (.not. exit_flag(itest,jdebug) ) then
    	      do k=1,ncap
                print *,   'DONNER_DEEP/cupar: k,qr,qt= ',k,   &
			    qr_v(itest,jdebug,k), &
			    qt_v(itest,jdebug,k)
     	      end do
 	      do k=1,nlev
 	        print *,   'DONNER_DEEP/cupar: k,dise,disa= ',k,   &
			    dise_v(itest,jdebug,k),   &
			    disa_v(itest,jdebug,k)
 	      end do
            endif
            endif
!         end do

!--------------------------------------------------------------------
!   call routine to close deep-cumulus parameterization.
!--------------------------------------------------------------------
      call cuclo_vect(        dcape_v,       pcape_v,plfc_v,plzb_v, &
		qli0_v,qli1_v,qr_v, &
          qt_v,r_v,     ri_v,rl_v,       rpc_v,     tcape_v,tpc_v,a1_v,&
	   exit_flag, nine_flag, jabs)


!--------------------------------------------------------------------
!   vertical integral of normalized moisture forcing:
!-------------------------------------------------------------------
      do j=jminp,jmaxp
        do i=iminp,imaxp
	  tfint(i,j) = 0.0
        end do
      end do


!     do k=2,kmax
      do k=2,nlev
        do j=jminp,jmaxp
          do i=iminp,imaxp
          if (.not. exit_flag(i,j) ) then
	    disbar = 0.5*(dise_v(i,j,k-1)+dise_v(i,j,k))
	    tfint(i,j) = tfint(i,j) - disbar*   &
			 (pr_v(i,j,k-1) - pr_v(i,j,k))
           endif
          end do
        end do
      end do

!      do k=1,kmax
       do k=1,nlev
       exit_3d(:,:,k) = exit_flag(:,:)
       end do
	

      do jlat=jminp,jmaxp
        do ilon=iminp,imaxp
          if (.not. exit_flag(ilon,jlat) ) then

!---------------------------------------------------------------------
!       restrict fractional area of first ensemble. see 
!       "a Bounds 6/7/97" notes.
!---------------------------------------------------------------------
 	    a1_v(ilon,jlat) = min(amax_v(ilon,jlat), a1_v(ilon,jlat))

!---------------------------------------------------------------------
!      fractional area further restricted by moisture constraint.
!      See "Moisture Constraint," 8/8/97.
!---------------------------------------------------------------------
            if ((tfint(ilon,jlat) .eq. 0.) .or.    &
	        (tpre_v(ilon,jlat) .eq. 0.)) then
	      a1_v(ilon,jlat) = 0.
	      nine_flag(ilon,jlat) = .true.
            end if
	  else
	      a1_v(ilon,jlat) = 0.
          endif
        end do
      end do

!          do j=jminp,jmaxp
!	    if (debug_jt(j)) then
             if (in_debug_window) then
       if (.not. exit_flag(itest,jdebug) .and.   &
	    .not. nine_flag(itest,jdebug))  then
                print   *, 'DONNER_DEEP/cupar: tfint= ',tfint(itest,jdebug)
                print   *, 'DONNER_DEEP/cupar: a1_v = ',a1_v (itest,jdebug)
            endif
            endif
!          end do

      do jlat=jminp,jmaxp
        do ilon=iminp,imaxp
          if (.not. exit_flag(ilon,jlat) .and.   &
	      .not. nine_flag(ilon,jlat))  then
            tfint(ilon,jlat)=tfint(ilon,jlat)/gravm
            amos_v(ilon,jlat) = (dqls_v(ilon,jlat) +     &
				 qlsd_v(ilon,jlat))/tfint(ilon,jlat)
            if (a1_v(ilon,jlat) .gt. amos_v(ilon,jlat))     &
	                a1_v(ilon,jlat) = max(amos_v(ilon,jlat),0.)

          endif
        end do
      end do

!-----------------------------------------------------------------------

!         do j=jminp,jmaxp
!	    if (debug_jt(j)) then
	    if (in_debug_window) then
       if (.not. exit_flag(itest,jdebug) .and.   &
	    .not. nine_flag(itest,jdebug))  then
                print   *, 'DONNER_DEEP/cupar: tfint,amos,a1,gravm= ',  &
			   tfint(itest,jdebug), amos_v(itest,jdebug), &
			   a1_v(itest,jdebug), gravm 
            endif
            endif
!          end do

      do k=1,nlev
      a1_vk(:,:,k) = a1_v(:,:)
      end do

!---------------------------------------------------------------------
!   if a1 allows negative mixing ratio, reset its value.
!--------------------------------------------------------------------
!     do k=1,kmax
      do k=1,nlev
        do j=jminp,jmaxp
          do i=iminp,imaxp
            if (.not. exit_flag(i,j) .and.   &
	        .not. nine_flag(i,j))  then
	      qtest(i,j,k) = q_v(i,j,k) + a1_vk(i,j,k)* &
		                   donner_deep_freq*dise_v(i,j,k)
              if (qtest(i,j,k) .lt. 0.) then
                a1_vk(i,j,k) = -q_v(i,j,k)/   &
		                  (dise_v(i,j,k)*donner_deep_freq)
               endif
            endif
          end do
        end do
      end do

!      do k=2,kmax
       do k=2,nlev
       do j=jminp,jmaxp
       do i=iminp,imaxp
       a1_v(i,j) = min(a1_vk(i,j,k), a1_v(i,j))
       end do
       end do
       end do

!         do j=jminp,jmaxp
!	    if (debug_jt(j)) then
	    if (in_debug_window) then
          if (.not. exit_flag(itest,jdebug) ) then
                print   *, 'DONNER_DEEP/cupar: a1= ',a1_v (itest,jdebug)
            endif
            endif
!         end do

!-----------------------------------------------------------------------
!      remove normalization from cumulus forcings.
!-----------------------------------------------------------------------
!      do j=1,kmax
       do j=1,nlev
      do jlat=jminp,jmaxp
      do ilon=iminp,imaxp
       if (.not. exit_3d  (ilon,jlat,j) ) then


	  disa_v(ilon,jlat,j)=disa_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  disb_v(ilon,jlat,j)=disb_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  disc_v(ilon,jlat,j)=disc_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  disd_v(ilon,jlat,j)=disd_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  dise_v(ilon,jlat,j)=dise_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  fre_v(ilon,jlat,j)=fre_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  elt_v(ilon,jlat,j)=elt_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  cmus_v(ilon,jlat,j)=cmus_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  ecds_v(ilon,jlat,j)=ecds_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  eces_v(ilon,jlat,j)=eces_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  emds_v(ilon,jlat,j)=emds_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  emes_v(ilon,jlat,j)=emes_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  wmms_v(ilon,jlat,j)=wmms_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  wmps_v(ilon,jlat,j)=wmps_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  tmes_v(ilon,jlat,j)=tmes_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  qmes_v(ilon,jlat,j)=qmes_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  cual_v(ilon,jlat,j)=cual_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  uceml_v(ilon,jlat,j)=uceml_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  umeml_v(ilon,jlat,j)=umeml_v(ilon,jlat,j)*a1_v(ilon,jlat)
	  dmeml_v(ilon,jlat,j)=dmeml_v(ilon,jlat,j)*a1_v(ilon,jlat)
	 endif
	 end do
	 end do
	 end do


!do k=kmax-4,kmax
	do k=nlev-4,nlev
!         do j=jminp,jmaxp
            do i=iminp,imaxp
!	      if (debug_jt(j)) then
               if (in_debug_window) then
                  if (.not. exit_flag(i,jdebug) ) then
                    if (disa_v(i,jdebug,k) .ne. 0.) then
 	              print *, 'DONNER_DEEP/cupar:kttest,k, jtest,i, disa= ',&
			 kttest, k, jtest, i,  disa_v(i,jdebug,k)
	              call error_mesg ('cupar_vect', &
		       'disa found to be /= 0.0 in top 5 levels', FATAL)
                    endif
                    if (fre_v(i,jdebug,k) .lt. 0.) then
 	              print *, 'DONNER_DEEP/cupar: kttest, k,jtest,i,fre,a1=',&
			      kttest, k,jtest,i, fre_v(i,jdebug,k),  &
					 a1_v(i,jdebug)
		      call error_mesg ('cupar_vect',  &
		                    ' fre_v found to be < 0.0', FATAL)
 	            end if
                    if (cmus_v(i,jdebug,k) .lt. 0.) then
 	              print *, 'DONNER_DEEP/cupar: kttest,k,jtest,i,cmus= ', &
				kttest, k,jtest,i, cmus_v(i,jdebug,k)
		      call error_mesg ('cupar_vect',  &
		                  ' cmus_v found to be < 0.0', FATAL)
                    end if	    
                endif
              endif
            end do
!         end do
        end do



!      do j=1,kmax
       do j=1,nlev
!      jinv=kmax+1-j
       jinv=nlev+1-j
      do jlat=jminp,jmaxp
!! THIS WAS THE LIKELY ERROR  !!
!         jgl = jabs(1)
 	  jgl = jabs(jlat)
      do ilon=iminp,imaxp
       if (.not. exit_3d  (ilon,jlat,j) ) then
	  cemetf(ilon,jgl ,jinv          )=disa_v(ilon,jlat,j)
	  ceefc(ilon,jlat,jinv          )=disb_v(ilon,jlat,j)
	  cecon(ilon,jlat,jinv          )=disc_v(ilon,jlat,j)
	  cemfc(ilon,jlat,jinv          )=disd_v(ilon,jlat,j)
	  cememf(ilon,jgl ,jinv          )=dise_v(ilon,jlat,j)
!  cememf_mod(ilon,jlat,jinv          )=dise_v(ilon,jlat,j)
	  cual_3d(ilon,jlat,jinv          )=cual_v(ilon,jlat,j)
	  fre_3d(ilon,jlat,jinv          )=fre_v(ilon,jlat,j)

          elt_3d(ilon,jlat,jinv) = elt_v(ilon,jlat,j)
	  cmus_3d(ilon,jlat,jinv          )=cmus_v(ilon,jlat,j)
	  ecds_3d(ilon,jlat,jinv          )=ecds_v(ilon,jlat,j)
	  eces_3d(ilon,jlat,jinv          )=eces_v(ilon,jlat,j)
	  emds_3d(ilon,jlat,jinv          )=emds_v(ilon,jlat,j)
	  emes_3d(ilon,jlat,jinv          )=emes_v(ilon,jlat,j)
	  qmes_3d(ilon,jlat,jinv          )=qmes_v(ilon,jlat,j)
	 wmps_3d(ilon,jlat,jinv          )=wmps_v(ilon,jlat,j)
	  wmms_3d(ilon,jlat,jinv           )=wmms_v(ilon,jlat,j)
	  tmes_3d(ilon,jlat,jinv          )=tmes_v(ilon,jlat,j)
	  dmeml_3d(ilon,jlat,jinv          )=dmeml_v(ilon,jlat,j)
	  uceml_3d(ilon,jlat,jinv          )=uceml_v(ilon,jlat,j)
	  umeml_3d(ilon,jlat,jinv          )=umeml_v(ilon,jlat,j)
	  cuqi_3d(ilon,jlat,jinv )  = cuq_v(ilon,jlat,j)
          cuql_3d(ilon,jlat,jinv )  = cuql_v(ilon,jlat,j)
      endif
      end do
      end do
       end do

!     if (get_my_pe() == 2) then
!print *, 'cemetf in cupar', (cemetf(33,3,k),k=1,nlev)
!     endif



      do jlat=jminp,jmaxp
      do ilon=iminp,imaxp
       if (.not. exit_flag(ilon,jlat) ) then
          a1_3d(ilon,jlat) = a1_v(ilon,jlat)
          amax_3d(ilon,jlat) = amax_v(ilon,jlat)
	  amos_3d(ilon,jlat         )=amos_v(ilon,jlat)
!	  tprea1_3d(ilon,jlat     )=tpre_v(ilon,jlat)*a1_v(ilon,jlat)
 	  tprea1_3d(ilon,jgl      )=tpre_v(ilon,jlat)*a1_v(ilon,jlat)
	  ampta1_3d(ilon,jlat    )=ampt_v(ilon,jlat)*a1_v(ilon,jlat)
	  rcoa1_3d(ilon,jlat)=rco_v(ilon,jlat)*a1_v(ilon,jlat)
	  emdi_v(ilon,jlat)=emdi_v(ilon,jlat)*a1_v(ilon,jlat)


      endif
      end do
      end do

!ljd
!         do jlat=jminp,jmaxp
!        do ilon=iminp,imaxp
!           if (ampt_v(ilon,jlat) .ne. 0.) then
!           write(6,*) 'cupar 3 ilon,jlat,jgl= ',ilon,jlat,jgl
!             write(6,*) 'cupar 3 ampt,tpre= ',ampta1_3d(ilon,jlat), &
!              tprea1_3d(ilon,jgl)
!             do k=1,nlev
!             write(6,*) 'k,cual= ',k,cual_3d(ilon,jlat,k)
!             end do
!           end if
!        end do
!        end do
!ljd
!ljd
!      do jlat=jminp,jmaxp
!      do ilon=iminp,imaxp
!              if (ampta1_3d(ilon,jlat) .ne. 0.) then
!             write(6,*) ' Cupar 4 i,j,jgl= ',ilon,jlat,jgl
!             write(6,*) 'emdi= ',emdi_v(ilon,jlat)
!             do k=1,nlev
!               write(6,*) 'k,umem= ',k,umeml_3d(ilon,jlat,k)
!             end do
!             end if
!             end do
!             end do
!ljd

!         do j=jminp,jmaxp
!	    if (debug_jt(j)) then
	    if ( in_debug_window) then
              do i=iminp,imaxp
                if (.not. exit_flag(i,jdebug) ) then
                  if (rcoa1_3d(i,jdebug ) .gt. 1.) then
 	            call error_mesg ('cupar_vect',  &
             ' the following columns contain deep convection', NOTE)
    	            print *, 'DONNER_DEEP/cupar: kttest,i,jtest,meso a= ',  &
			      kttest,i, jtest,rcoa1_3d(i,jdebug)
           print *, 'DONNER_DEEP/cupar: amax,ampt= ',amax_3d(i,jdebug), &
			      ampta1_3d(i,jdebug)
 	          end if
                endif
              end do
            endif
!	  end do

!          do j=jminp,jmaxp
!	    if (debug_jt(j)) then
	    if (in_debug_window) then
              if (.not. exit_flag(itest,jdebug) ) then
                if (tprea1_3d(itest,jtest ) .ne. 0.) then
	          do k=1,nlev
                    if ((pr_v(itest,jdebug,k) .gt. 100.e02) .and.   &
		        (pr_v(itest,jdebug,k) .lt. 500.e02)) then
                      if (disa_v(itest,jdebug,k) .ne. 0.) then
                        print *,'DONNER_DEEP/cupar: kttest,j,itest,k,t= ',  &
				kttest, j, itest,k,t_v(itest,jdebug,k)
                        print *,'DONNER_DEEP/cupar: tprea1,k,pr,cemetf= ',  &
			         tprea1_3d(itest,jtest), k,    &
				 pr_v(itest,jdebug,k),   &
				 cemetf(itest,jtest ,k )
	              endif
	            endif
	          end do
	        endif
	      endif
	    endif
!          end do

!         do j=jminp,jmaxp
!	    if (debug_jt(j)) then
	    if (in_debug_window) then
	      print *, 'DONNER_DEEP/cupar: contot,tpre=',   &
		       contot_v(itest,jdebug), tpre_v(itest,jdebug)
	      print *, 'DONNER_DEEP/cupar: a1,ampt =', a1_3d (itest,jdebug), &
			 ampta1_3d(itest,jdebug)
	    endif
!         end do

!         do j=jminp,jmaxp
!	    if (debug_jt(j)) then
	    if (in_debug_window) then
              print *,   'DONNER_DEEP/cupar: amax= ',amax_3d(itest,jdebug)
              do k=1,nlev
                print *, 'DONNER_DEEP/cupar: k,pr,uceml,dmeml,umeml= ',  &
			  k,pr_v(itest,jdebug,k),    &
			  uceml_3d(itest,jdebug,k), &
			  dmeml_3d(itest,jdebug,k),  &
			  umeml_3d(itest,jdebug,k)
              end do
            end if
!         end do


      if (test) then

        do j=jminp,jmaxp
        do i=iminp,imaxp
          if (.not. exit_flag(i,j) ) then
!    do k=1,kmax
	    do k=1,nlev
      	      sig(k)  = pr_v(ilon,jlat,k)/pr_v(ilon,jlat,1)
 	      disa(k) = disa_v(i,j,k)*86400.
 	      disb(k) = disb_v(i,j,k)*86400.
 	      disc(k) = disc_v(i,j,k)*86400.
       	      disd(k) = disd_v(i,j,k)*8.64e07
 	      dise(k) = dise_v(i,j,k)*8.64e07
 	      fre(k)  = fre_v (i,j,k)*86400.
 	      elt(k)  =elt_v (i,j,k)*86400.
 	      cmus(k) =cmus_v(i,j,k)*8.64e07
 	      ecds(k) =ecds_v(i,j,k)*8.64e07
 	      eces(k) = eces_v(i,j,k)*8.64e07
 	      emds(k)= emds_v(i,j,k)*8.64e07
 	      emes(k)=emes_v(i,j,k)*8.64e07
 	      wmms(k)=wmms_v(i,j,k)*8.64e07
 	      wmps(k)=wmps_v(i,j,k)*8.64e07
 	      tmes(k)=tmes_v(i,j,k)*86400.
 	      qmes(k)=qmes_v(i,j,k)*8.64e07
	    end do


!      CALL AGSETF('Y/ORDER.',1.)
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISB,SIG ,nlev,'ENTROPY FLUX CONVERGENCE')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISD,SIG ,nlev,'MOISTURE FLUX CONVERGENCE')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISC,SIG,nlev,'LATENT HEAT RELEASE')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(FRE ,SIG,nlev,'FREEZE ')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(ELT, SIG,nlev,'MELT ')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISA ,SIG,nlev,'CUMULUS THERMAL FORCING')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISE, SIG,nlev,'CUMULUS MOISTURE FORCING')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(CMUS, SIG,nlev,'MESOSCALE UPDRAFT CONDENSATION')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(ECDS, SIG,nlev,'CONVECTIVE DOWNDRAFT EVAPORATION')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(ECES, SIG,nlev,'CONVECTIVE UPDRAFT EVAPORATION')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(EMDS, SIG,nlev,'MESOSCALE DOWNDRAFT EVAPORATION')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(EMES, SIG,nlev,'MESOSCALE UPDRAFT EVAPORATION')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(WMMS, SIG,nlev,'MESOSCALE CELL CONDENS')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(WMPS, SIG,nlev,'MESOSCALE VAPOR Redistribution')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(tmes,SIG ,nlev,'Meso EFC')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(qmes,SIG ,nlev,'Meso MFC')
!      call clsgks
     end if

    end do
    end do

    endif


end subroutine cupar_vect



subroutine precu_vect( temold, ratold, press, pcape, prini, &
			      jabs, &
			     rcape,  &
                             rini, tcape, tini)

integer, dimension(:),  intent(in)  :: jabs
real, dimension(:,:,:), intent(in) :: temold, ratold, press
real, dimension(:,:,:), intent(out) :: pcape, prini, rcape, rini, &
				       tcape, tini

!-------------------------------------------------------------------
!     precu_vect defines pressure, temperature and moisture values on
!     an enhanced vertical grid that will be used in the deep cumulus
!     convection parameterization (Donner, 1993, JAS).
!-------------------------------------------------------------------
!     On Input:
!
!        ilon          longitude index
!        jlat          latitude index
!        nlev          Cape.F resolution (parameter)
!        temold        temperautre at GCM full levels (deg C)
!                      For third array index, 1 nearest surface.
!
!     On Output:
!
!        pcape         pressure at Cape.F resolution  (Pa)
!                      Index 1 nearest surface.
!        rcape         vapor mixing ratio at Cape.F resolution (kg(H2O)/kg)
!                      Index 1 nearest surface.
!        rini          vapor mixing ratio at Skyhi resolution (kg(H2O)/kg)
!                      Index 1 nearest surface.
!        tcape         temperautre at Cape.F resolution (K)
!        tini          temperaute at Skyhi resolution (K)
!                      Index 1 nearest surface.
!        prini         pressure at Skyhi resolution (Pa)
!                      Index 1 nearest surface.
!!
!

      real, dimension (size(temold,1), size(temold,2) ) :: dp

!--------------------------------------------------------------------
!   create arrays of moisture, temperature (K) and pressure with 
!   index 1 being closest to the surface. require that the moisture 
!   value be  non-negative.
!--------------------------------------------------------------------
!     do k=1,kmax
      do k=1,nlev
        do j=jminp,jmaxp
          do i=iminp,imaxp
	    rini (i,j,nlev+1-k) = amax1 (ratold(i,j,k), 0.0e00)
!           tini (i,j,nlev+1-k) = temold(i,j,k) + frezdk
            tini (i,j,nlev+1-k) = temold(i,j,k)
	    prini(i,j,nlev+1-k) = press(i,j,k)
          end do
        end do
      end do

!-------------------------------------------------------------------
!   define the vertical resolution of the convection parameterization
!   grid. define the top level pressure in that grid to be zero.
!   interpolate to define the pressure levels of that grid.
!-------------------------------------------------------------------
      do j=jminp,jmaxp
        do i=iminp,imaxp
          dp(i,j) = (prini(i,j,1) - prini(i,j,nlev))/(ncap-1)
          pcape(i,j,ncap) = 0.
        end do
      end do
      do k=1,ncap-1
        do j=jminp,jmaxp
          do i=iminp,imaxp
	    pcape(i,j,k) = prini(i,j,1) - (k-1)*dp(i,j)
          end do
        end do
      end do

!--------------------------------------------------------------------
!   call polat_vect to define values of temperature and moisture on
!   the enhanced convection grid by interpolating from the model grid
!   values. insure that the moisture field is positive-definite.
!--------------------------------------------------------------------
      call polat_vect (tini, prini, pcape, tcape)
      call polat_vect (rini, prini, pcape, rcape)
      rcape(:,:,:) = MAX (rcape(:,:,:), 0.0)

!---------------------------------------------------------------------
!   if debugging is activated, print out the input pressure, moisture 
!   and temperature fields in the  desired column.
!---------------------------------------------------------------------
!         do j=jminp,jmaxp
!    if (debug_jt(j)) then
	    if (in_debug_window) then
	      do k=1,nlev-ktest_model+1
		print *, 'DONNER_DEEP/precu: k, prini,tini,rini = ', k,   &
			    prini(itest,jdebug,k),   &
			    tini(itest,jdebug,k), &
			    rini(itest,jdebug,k)
              end do
            endif
!         end do

!--------------------------------------------------------------------

end subroutine precu_vect


!###################################################################

subroutine polat_vect (xv, pv, p, x)

!------------------------------------------------------------------
!    polat_vect interpolates the field xv on a  pressure grid pv to
!    a field x on a pressure grid p.
!------------------------------------------------------------------
 
real, dimension(:,:,:), intent(in)   :: xv, pv, p
real, dimension(:,:,:), intent(out)  :: x

!     ON INPUT:
!
!         XV(N)   DATA AT RESOLUTION N
!         PV(N)   PRESSURE AT N LEVELS
!
!     ON OUTPUT:
!
!         X       DATA AT PRESSURE P
!
!       PARAMETER(N=40      ,NM1=N-1)
!       DIMENSION XV(iminp:imaxp,jminp:jmaxp,N),    &
!		 PV(iminp:imaxp,jminp:jmaxp,N)
!       dimension x (iminp:imaxp, jminp:jmaxp,ncap)
!       dimension p (iminp:imaxp, jminp:jmaxp,ncap)


!       call tic ( 'polat_vect', '1')
       n = size(pv,3)
       nm1 = n - 1

       do k=1,ncap
       do j=jminp,jmaxp
       do i=iminp,imaxp
       IF (P(i,j,k) .GE. PV(i,j,1)) THEN
          X(i,j,k)=(XV(i,j,2)-XV(i,j,1))/(PV(i,j,2)-PV(i,j,1))
          X(i,j,k)=X(i,j,k)*(P(i,j,k)-PV(i,j,1))+XV(i,j,1)
       else IF (P(i,j,k) .LE. PV(i,j,n)) THEN
          X(i,j,k)=(XV(i,j,n)-XV(i,j,nm1))/(PV(i,j,n)-PV(i,j,nm1))
          X(i,j,k)=X(i,j,k)*(P(i,j,k)-PV(i,j,n))+XV(i,j,n)
       ENDIF
       end do
       end do
       end do

!       call toc ( 'polat_vect', '1')
!       call tic ( 'polat_vect', '2')

       do k=1,ncap
       do j=jminp,jmaxp
       do i=iminp,imaxp
       DO  kk=1,NM1
       IF ((P(i,j,k) .GE. PV(i,j,kk+1)) .AND.     &
	   (P(i,j,k) .LE. PV(i,j,kk))) THEN
          X(i,j,k)=(XV(i,j,kk+1)-XV(i,j,kk))/(PV(i,j,kk+1)-PV(i,j,kk))
          X(i,j,k)=X(i,j,k)*(P(i,j,k)-PV(i,j,kk+1))+XV(i,j,kk+1)
	  cycle
       ENDIF
       end do

       end do
       end do
       end do

!       call toc ( 'polat_vect', '2')
!       call tic ( 'polat_vect', '3')

!      do j=jminp,jmaxp
! if (debug_jt(j)) then
	 if (in_debug_window) then
	   do k=1,ncap
             print *, 'DONNER_DEEP/polat: k,p,x=', k, p(itest,jdebug,k),    &
						    x(itest,jdebug,k)
           end do
         endif
!      end do

!       call toc ( 'polat_vect', '3')



end subroutine polat_vect



!###################################################################

subroutine cape_vect (p_v, r_v, t_v, cin_v, plcl_v, plfc_v, plzb_v, &
		      rpc_v, tpc_v, xcape_v, jabs)
 
!----------------------------------------------------------------------
!    cape_vect calculates convective available potential energy for a 
!    cloud whose temperature follows a saturated adiabat.
!---------------------------------------------------------------------

integer, dimension(:),   intent(in)  :: jabs
real, dimension (:,:,:), intent(in)  :: p_v, r_v, t_v
real, dimension (:,:)  , intent(out) :: cin_v, plcl_v, plfc_v, &
					plzb_v, xcape_v
real, dimension (:,:,:), intent(out) :: rpc_v, tpc_v

!--------------------------------------------------------------------
!  intent(in) variables:
!
!      p_v     pressure (Pa)
!              Index 1 refers to level nearest earth's surface.
!      r_v     mixing ratio (kg(H2O)/kg)
!              Index 1 refers to level nearest earth's surface.
!      t_v     temperature (K)
!              Index 1 refers to level nearest earth's surface.
!
!   intent(out) variables:
!
!     cin_v    convective inhibition (J/kg)
!              energy required to lift parcel from level istart to
!              level of free convection
!     plcl_v   pressure at lifting condensation level (Pa)
!     plfc_v   pressure at level of free convection (Pa)
!              height of plfc .le. height of plcl
!              If parcel becomes buoyant below plcl, cin can be .lt. 0
!     plzb_v   pressure at level of zero buoyancy (Pa)
!     rpc_v    parcel mixing ratio
!              Set to environment below istart.
!              Index 1 refers to level nearest earth's surface.
!     tpc_v    parcel temperature (K)
!              Set to environment below istart.
!              Index 1 refers to level nearest earth's surface.
!     xcape_v  convective available potential energy (J/kg)
!              energy released as parcel moves from level of free
!              convection to level of zero buoyancy
!
!     Calculated:
!   
!     tot      xcape+cin (J/kg)
!--------------------------------------------------------------------
!

   integer, dimension(size(t_v,1), size(t_v,2) ) :: klcl, klzb, klfc, &
                        ieqv      

   logical, dimension(size(t_v,1), size(t_v,2) ) :: capepos, cape_exit,&
				lcl_found, lzb_found, skip_search, &
				cin_done, cape_skip

      real, dimension(size(t_v,1), size(t_v,2) ) :: tot, ro, tc_1, tp, &
				  q_ve, cp_v, rs_v, es_v, tlcl, rlcl, &
				  dt_v, dtdp_v, tc_v, qe_v, qs_v, &
				  tve_v, tv_v, tvc_v, tve_v, delt_v

     integer :: s_found, z_found, cin_counter, error_flag



!--------------------------------------------------------------------
!     call tic ('cape', '1')
!         do j=jminp,jmaxp
!	    if (debug_jt(j)) then
	    if (in_debug_window) then
	      print *, 'DONNER_DEEP/cape: cpi= ', cpi
	      print *, 'DONNER_DEEP/cape: cpv= ', cpv
	      print *, 'DONNER_DEEP/cape: rocp= ', rocp
	      print *, 'DONNER_DEEP/cape: rair= ', rair
	      print *, 'DONNER_DEEP/cape: latvap= ', latvap
	      print *, 'DONNER_DEEP/cape: rvap= ',   rvap
	      do k=1,ncap
		print *, 'DONNER_DEEP/cape: k, p = ', k, p_v  (itest,jdebug,k)
              end do
	      do k=1,ncap
		print *, 'DONNER_DEEP/cape: k,  r= ', k,  r_v (itest,jdebug,k)
              end do
	      do k=1,ncap
		print *, 'DONNER_DEEP/cape: k,  t = ', k,  t_v(itest,jdebug,k)
              end do
            endif
!         end do

!--------------------------------------------------------------------
!     Stop CAPE calculation when pressure falls to pstop or
!     parcel temperature falls below tmin.
!     istart-index of level whose mixing ratio is conserved as a parcel 
!            leaves it undergoing dry adiabatic ascent
!--------------------------------------------------------------------
      pstop=4.0e03
      tmin=154.
      istart=1

!-------------------------------------------------------------------
!     initialize variables.
!------------------------------------------------------------------
      do j=jminp,jmaxp
        do i=iminp,imaxp
          capepos(i,j)=.false.
          plfc_v(i,j)=pstop
          plzb_v(i,j)=pstop
          plcl_v(i,j)=pstop
          klfc(i,j)=ncap-1
          klcl(i,j)=ncap-1
          cin_v(i,j)=0.
          xcape_v(i,j)=0.
          tot(i,j)=0.
        end do
      end do
      do k=1,ncap
        do j=jminp,jmaxp
          do i=iminp,imaxp
	    tpc_v(i,j,k)=t_v(i,j,k)
	    rpc_v(i,j,k)=r_v(i,j,k)
          end do
        end do
      end do
!     call toc ('cape', '1')
!     call tic ('cape', '2')
!--------------------------------------------------------------------
!  define parcel departure point values. if temperature at departure
!  point is too cold, set flag to suppress calculations in that column.
!  convert mixing ratio to specific humidity.
!--------------------------------------------------------------------

      s_found = 0
      z_found = 0
      cin_counter = 0

      do j=jminp,jmaxp
        do i=iminp,imaxp
          ro(i,j)=r_v(i,j,istart)
          tc_1(i,j)=t_v(i,j,istart)
          tp(i,j)=tc_1(i,j)
          if (tp(i,j) .lt. tmin)  then
	    cape_exit(i,j) = .true.
	    s_found = s_found + 1
	    lcl_found(i,j) = .true.
	    z_found = z_found + 1
	    lzb_found(i,j) = .true.
	  else
	    cape_exit(i,j) = .false.
	    lcl_found(i,j) = .false.
	    lzb_found(i,j) = .false.
          endif
          q_ve(i,j) = ro(i,j)/(1.+ro(i,j))
          cp_v(i,j) = cpi*(1.+((cpv/cpi)-1.)*q_ve(i,j))
        end do
      end do

!     call toc ('cape', '2')
!--------------------------------------------------------------------
!  integrate in the vertical to find the lcl in each column. only 
!  calculate in those columns which are still searching.
!--------------------------------------------------------------------
!     call tic ('cape', '3')
      do k=istart,ncap
!if (s_found < 100) then
	if (s_found < imaxp*jmaxp) then

!--------------------------------------------------------------------
!  determine saturation mixing ratio for parcels at this level.
!---------------------------------------------------------------------
!         call establ_vect(es_v, tp, lcl_found)
!         call lookup_es     (tp, es_v, lcl_found)
          call lookup_es     (tp, es_v)
          do j=jminp,jmaxp
            do i=iminp,imaxp
	      if (.not. lcl_found(i,j) ) then
	        rs_v(i,j) = rocp*es_v(i,j)/(p_v(i,j,k) +  &
			    (rocp - 1.)*es_v(i,j))
              endif
	    end do
	  end do

!--------------------------------------------------------------------
!  check for parcel saturation at this level, if it has not yet been
!  reached and level is below upper limit for cape calculation.
!---------------------------------------------------------------------
          call ieq_z (rs_v, ro, ieqv, lcl_found)
          do j=jminp,jmaxp
            do i=iminp,imaxp
	      if (.not. lcl_found(i,j) ) then
 	        if (p_v(i,j,k) .ge. pstop) then    

!--------------------------------------------------------------------
!   saturation reached exactly; save pressure, temp, moisture and 
!   vertical index for cloud base.  set flag and update counter.
!--------------------------------------------------------------------
	          if (ieqv(i,j) .eq. 0) then
	            plcl_v(i,j) = p_v(i,j,k)
	            rlcl(i,j) = r_v(i,j,k)
	            tlcl(i,j) = t_v(i,j,k)
	            klcl(i,j) = k
		    s_found = s_found + 1
                    lcl_found(i,j) = .true.

!--------------------------------------------------------------------
!   saturation exceeded; save pressure, temp, moisture and vertical
!   index for cloud base. set flag and update counter.
!--------------------------------------------------------------------
                  else if (ieqv(i,j) .lt. 0) then
	            if (k .eq. istart) then
	              plcl_v(i,j)=p_v(i,j,istart)
	              tlcl(i,j) =t_v(i,j,istart)
	              rlcl(i,j) =r_v(i,j,istart)
	              klcl(i,j)=istart
		      s_found = s_found + 1
                      lcl_found(i,j) = .true.
	            else                         
	              plcl_v(i,j) = (p_v(i,j,k)+p_v(i,j,k-1))/2.
	              tlcl(i,j) = (t_v(i,j,k)+t_v(i,j,k-1))/2.
	              rlcl(i,j) = (r_v(i,j,k)+r_v(i,j,k-1))/2.
	              klcl(i,j) = k
		      s_found = s_found + 1
                      lcl_found(i,j) = .true.
	   	    endif

!---------------------------------------------------------------------
!    parcel still unsaturated; move parcel to next pressure level 
!    along the dry adiabat. define temperature at this level; verify
!    that it is within bounds of scheme. update parcel values at the
!    k+1 level to be those thus calculated.
!---------------------------------------------------------------------
	   	  else  ! (ieqv(i,j) .gt. 0) 
	   	    if (k .lt. ncap) then
                      dtdp_v(i,j)=rair*tp(i,j)/cp_v(i,j)
                      dt_v(i,j)=dtdp_v(i,j)*    &
				    alog((p_v(i,j,k+1)   )/p_v(i,j,k))
                      tp(i,j)=tp(i,j)+dt_v(i,j)
      	              if (tp(i,j) .lt. tmin)  then
	                cape_exit(i,j) = .true.
		        lcl_found(i,j) = .true.
		        s_found = s_found + 1
			z_found = z_found + 1
			lzb_found(i,j) = .true.
                      else  
	                tpc_v(i,j,k+1)=tp(i,j)
	                rpc_v(i,j,k+1)=ro(i,j)
                      endif
	            endif
		  endif ! (ieqv(i,j) .gt. 0)

!--------------------------------------------------------------------
!    if parcel is now at pressure below cutoff for cape calculation,
!    set flag to stop calculation.
!--------------------------------------------------------------------
	        else
		  lcl_found(i,j) = .true.
		  s_found = s_found + 1
 	        endif      ! (p(k) .lt. pstop )
	      endif  ! (lcl_found)
            end do
          end do
        endif
      end do   ! k loop

!     call toc ('cape', '3')
!--------------------------------------------------------------------
!  verify that the lcl calculated is below the limiting pressure level.
!--------------------------------------------------------------------
!     call tic ('cape', '4')
      call ieq_y (plcl_v, pstop, ieqv)
      do j=jminp,jmaxp
        do i=iminp,imaxp
          if (ieqv(i,j) .le. 0) then
            cape_exit(i,j) = .true.
	    if ( .not. (lzb_found(i,j))) then
	      z_found = z_found + 1
	      lzb_found(i,j) = .true.
	    endif
          end if
        end do
      end do

!     call toc ('cape', '4')
!     call tic ('cape', '5')


!        do j=jminp,jmaxp
!	   if (debug_jt(j)) then
           if (in_debug_window) then
             print *,  'DONNER_DEEP/cape: plcl,klcl,tlcl,rlcl= ',   &
			 plcl_v(itest,jdebug), klcl(itest,jdebug),    &
			 tlcl(itest,jdebug),rlcl(itest,jdebug)
             print *,  'DONNER_DEEP/cape: p(klcl)= ',   &
			    p_v(itest,jdebug,klcl(itest,j))
           endif
!        end do

!-------------------------------------------------------------------
!   calculate temperature along saturated adiabat, starting at p(klcl)
!   and a temperature tp.
!-------------------------------------------------------------------
      do j=jminp,jmaxp
        do i=iminp,imaxp
            tc_v(i,j)=tp(i,j)
      end do
      end do
!     call toc ('cape', '5')

!--------------------------------------------------------------------
!   integrate in the vertical to find the level of free convection and
!   the level of zero buoyancy. only integrate in those columns which
!   are still active.
!--------------------------------------------------------------------
      do k=istart,ncap-1
!	if (z_found < 100) then
	if (z_found < imaxp*jmaxp) then

!--------------------------------------------------------------------
!   define a flag indicating whether one should search for either the
!   lfc or lzb at this level in each column. 
!--------------------------------------------------------------------
!	  call tic ('cape', '5a')
          do j=jminp,jmaxp
            do i=iminp,imaxp
               if (p_v(i,j,k+1) >= pstop ) then
                 if ( (.not. lzb_found(i,j)) .and.     &
		      (k .ge. klcl(i,j)) ) then
                  skip_search(i,j) = .false.
                else
                  skip_search(i,j) = .true.
                endif
	      else
                skip_search(i,j) = .true.
              endif  ! p >= pstop
            end do
          end do

!	  call toc ('cape', '5a')
!	  call tic ('cape', '5b')
!--------------------------------------------------------------------
!   define saturation vapor pressure for the parcel in active columns.
!--------------------------------------------------------------------
!         call establ_vect (es_v, tc_v, skip_search)
!         call lookup_es      (tc_v, es_v, skip_search)
          call lookup_es      (tc_v, es_v)

!	  call toc ('cape', '5b')
!	  call tic ('cape', '5c')
!--------------------------------------------------------------------
!   define the environmental and parcel virtual temperature and specific
!   humidity where needed.
!--------------------------------------------------------------------
          do j=jminp,jmaxp
            do i=iminp,imaxp
	      if ( .not. skip_search(i,j)) then
	        qe_v(i,j) = r_v(i,j,k)/(1.+r_v(i,j,k))
	        tve_v(i,j) = t_v(i,j,k)*(1.+.61*qe_v(i,j))
     	        rs_v(i,j) = rocp  *es_v(i,j)/(p_v(i,j,k)+   &
		      	    (rocp  -1.)*es_v(i,j))
                qs_v(i,j)=rs_v(i,j)/(1.+rs_v(i,j))
	        tv_v(i,j)=tc_v(i,j)*(1.+.61*qs_v(i,j))
	      endif
	    end do
	  end do
!	  call toc ('cape', '5c')

!--------------------------------------------------------------------
!   determine whether the parcel temperature is cooler or warmer than 
!   the environment.
!--------------------------------------------------------------------
!	  call tic ('cape', '5d')
          call ieq_z (tv_v, tve_v, ieqv, skip_search  )
!	  call toc ('cape', '5d')

!---------------------------------------------------------------------
!   integrate parcel upward, finding level of free convection and 
!   level of zero buoyancy.
!---------------------------------------------------------------------
!	  call tic ('cape', '5e')
          do j=jminp,jmaxp
            do i=iminp,imaxp
	      if ( .not. skip_search(i,j)) then

!!!COMMENT: ieqv == 0 should be included below:
!!!             if ((ieqv(i,j) .gt. 0) .and. (.not. capepos(i,j))) then

!-------------------------------------------------------------------
!   determine if the level of free convection has been reached. 
!-------------------------------------------------------------------
                if ((ieqv(i,j) .ge. 0) .and. (.not. capepos(i,j))) then
                  capepos(i,j) = .true.
                  plfc_v(i,j) = p_v(i,j,k)
                  klfc(i,j) = k
                end if

!-------------------------------------------------------------------
!   determine if the level of zero buoyancy has been reached.  if so,
!   set flag so that calculation will be ended in this column.
!-------------------------------------------------------------------
	        if ((ieqv(i,j) .lt. 0) .and. (capepos(i,j))) then
	          klzb(i,j) = k
	          plzb_v(i,j) = (p_v(i,j,k)+p_v(i,j,k-1))/2.
                  z_found = z_found + 1
	          lzb_found(i,j) = .true.
!-------------------------------------------------------------------
!   if not, continue moving parcel up pseudo-adiabat to next cape-
!   calculation pressure level. define new parcel temperature and
!   mixing ratio at this level; if temperature is colder than allowed,
!   end integration.
!-------------------------------------------------------------------
                else   !  (cape is pos, parcel warmer than env)
	          rc = (1.-qs_v(i,j))*rair+qs_v(i,j)*rvap
	          pb = 0.5*(p_v(i,j,k)+p_v(i,j,k+1))
	          fact1 = rair/cpi
	          fact2 = tv_v(i,j)+(latvap*qs_v(i,j)/rc)
	          fact1 = fact1*fact2
	          fact3 = rocp  *(latvap**2)*es_v(i,j)/    &
		         (cpi*pb*rvap*(tv_v(i,j)**2))
	          fact3 = 1.+fact3
	          dtdp = fact1/fact3
	          tc_v(i,j) = tc_v(i,j)+dtdp*     &
			      alog(p_v(i,j,k+1)/p_v(i,j,k))
           	  if (tc_v(i,j) .lt. tmin        )  then
	            cape_exit(i,j) = .true.
	            z_found = z_found + 1
	            lzb_found(i,j) = .true.
	          else
	            tpc_v(i,j,k+1) = tc_v(i,j)
	            rpc_v(i,j,k+1) = rs_v(i,j)
                  endif
                endif   !  (ieq < 0, capepos)

!--------------------------------------------------------------------
!    if pressure has gone below the minimum at which deep convection 
!    is allowed, set flag to end calculation in this column.
!--------------------------------------------------------------------
              else  ! (skip_search)
                if (p_v(i,j,k+1)  < pstop ) then
		 if (.not. lzb_found(i,j) ) then
	            lzb_found(i,j) = .true.
	            z_found = z_found + 1
		  endif
                endif                    
              endif    !   (skip_search)
            end do  ! i loop
          end do  ! j loop
!	  call toc ('cape', '5e')
        endif   !  (z_found)
      end do   ! k loop

!      call tic ('cape', '5f')
!         do j=jminp,jmaxp
!    if (debug_jt(j)) then
	    if (in_debug_window) then
              do k=klcl(itest,j),ncap-1
             print *,   'DONNER_DEEP/cape: k,tv,tve= ',k,tv_v(itest,jdebug),  &
				    tve_v(itest,jdebug)
 	        print *,   'DONNER_DEEP/cape: klzb,plzb,p(klzb)= ',  &
			     klzb(itest,jdebug),   &
 		            plzb_v(itest,jdebug),   &
			    p_v(itest,jdebug,klzb(itest,jdebug))
 	        print *,   'DONNER_DEEP/cape: fact1,fact2,rc= ',fact1,fact2,rc
 	        print *,   'DONNER_DEEP/cape: fact1,fact3= ',fact1,fact3
 	        print *,   'DONNER_DEEP/cape: dtdp= ',dtdp
 	        print *,   'DONNER_DEEP/cape: tc,t= ',tc_v(itest,jdebug),  &
					   t_v(itest,jdebug,k+1)
       print *,   'DONNER_DEEP/cape: p,r,rs= ',p_v(itest,jdebug,k+1),   &
			    r_v(itest,jdebug,k+1),  rs_v(itest,jdebug)
	      end do  ! k loop
            endif
!         end do
!      call toc ('cape', '5f')

!-------------------------------------------------------------------
!   define flag to indicate where to calculate convective inhibition.
!-------------------------------------------------------------------
!      call tic ('cape', '6')
      do j=jminp,jmaxp
        do i=iminp,imaxp
	  if (.not. cape_exit(i,j)) then
	    cin_done(i,j) = .false.
          else
	    cin_done(i,j) = .true.
	    cin_counter = cin_counter + 1
	  endif
        end do
      end do

!-------------------------------------------------------------------
!   calculate convective inhibition.
!--------------------------------------------------------------------
      do k=istart,ncap-1
! 	if (cin_counter < 100) then
 	if (cin_counter < imaxp*jmaxp) then

!------------------------------------------------------------------
!   determine if sounding fails to produce a level of free convection.
!   if so, set flag to avoid cape calculation. If desired, print out
!   columns where lcl exists, but no lfc.
!------------------------------------------------------------------
        do j=jminp,jmaxp
          do i=iminp,imaxp
             if (.not. (cin_done (i,j)) ) then
              if ( p_v(i,j,k+1) <= pstop              ) then
		cin_done(i,j) = .true.
	        cape_exit(i,j) = .true.
		cin_counter = cin_counter + 1
                if (debug) then
                  print *, 'cape = 0 (NO LFC): i, jrow, cin= ',   &
			   i, jabs(j), cin_v(i,j)
		endif
	      end if
            endif
          end do
	end do

!--------------------------------------------------------------------
!    define the specific humidity and virtual temperature of the
!    parcel and environment.
!--------------------------------------------------------------------
        do j=jminp,jmaxp
          do i=iminp,imaxp
 	    if ( .not. cin_done(i,j) ) then
	      rbc =     (rpc_v(i,j,k) + rpc_v(i,j,k+1))/2.
	      rbe =     (r_v(i,j,k) + r_v(i,j,k+1))/2.
	      qc = rbc/(1. + rbc)
	      qe = rbe/(1. + rbe)
	      tvc_v(i,j) =     (tpc_v(i,j,k) + tpc_v(i,j,k+1))/2.
	      tve_v(i,j) =     (t_v(i,j,k) + t_v(i,j,k+1))/2.
	      tvc_v(i,j) = tvc_v(i,j)*(1.+.61*qc)
	      tve_v(i,j) = tve_v(i,j)*(1.+.61*qe)
             endif
          end do
	end do

!           do j=jminp,jmaxp
!	      if (debug_jt(j)) then
               if (in_debug_window) then
	      if (.not. cin_done(itest,jdebug)) then
                print *, 'DONNER_DEEP/cape: k,tvc,tve= ', k,  &
			   tvc_v(itest,jdebug),  &
			  tve_v(itest,jdebug)
              endif
              endif
!           end do

!---------------------------------------------------------------------
!   determine whether the parcel temperature is cooler or warmer than 
!   the environment.
!--------------------------------------------------------------------
        call ieq_z(tvc_v, tve_v, ieqv, cin_done)

!---------------------------------------------------------------------
!   add the contribution to cin from this pressure layer.
!---------------------------------------------------------------------
        do j=jminp,jmaxp
          do i=iminp,imaxp
	    if ( .not. cin_done(i,j) ) then
	      if ((ieqv(i,j) .lt. 0) .or.      &
 		  (p_v(i,j,k) .gt. plfc_v(i,j)))  then
                delt_v(i,j) = rair*(tvc_v(i,j)-tve_v(i,j))*   &
			 alog(p_v(i,j,k)/p_v(i,j,k+1))
                cin_v(i,j) = cin_v(i,j) - delt_v(i,j)
	      else
                  cin_done(i,j) = .true.
	  	  cin_counter = cin_counter + 1
              end if
            endif
          end do
        end do
      endif
      end do
!      call toc ('cape', '6')

!      call tic ('cape', '7')
!-------------------------------------------------------------------
!   if desired, print out lfc k index and pressure.
!-------------------------------------------------------------------
!         do j=jminp,jmaxp
!    if (debug_jt(j)) then
	    if (in_debug_window) then
           print *, 'DONNER_DEEP/cape: klfc, p(klfc)= ', klfc(itest,jdebug),  &
	               p_v(itest,jdebug,klfc(itest,jdebug))
            endif
!          end do

!--------------------------------------------------------------------
!  calculate convective available potential energy.
!--------------------------------------------------------------------
      do k=istart,ncap-1

!--------------------------------------------------------------------
!  define flag to indicate which columns are actively computing cape.
!-------------------------------------------------------------------
        do j=jminp,jmaxp
          do i=iminp,imaxp
            if (.not. (cape_exit(i,j))  .and.  &
	         ( k .ge. klfc(i,j) ) ) then
              if ( p_v(i,j,k+1) .gt. plzb_v(i,j) ) then
	        cape_skip(i,j) = .false.
	      else
	        cape_skip(i,j) = .true.
              endif
	    else
	      cape_skip(i,j) = .true.
            endif
	  end do 
	end do 

!--------------------------------------------------------------------
!  define virtual temperature and specific humidity of parcel and 
!  environment.
!-------------------------------------------------------------------
        do j=jminp,jmaxp
          do i=iminp,imaxp
	    if (.not. cape_skip(i,j) ) then
              rbc = (rpc_v(i,j,k)+rpc_v(i,j,k+1))/2.
	      rbe = (r_v(i,j,k)+r_v(i,j,k+1))/2.
	      qc = rbc/(1.+rbc)
	      qe = rbe/(1.+rbe)
	      tvc_v(i,j) = (tpc_v(i,j,k)+tpc_v(i,j,k+1))/2.
	      tve_v(i,j) = (t_v(i,j,k)+t_v(i,j,k+1))/2.
	      tvc_v(i,j) = tvc_v(i,j)*(1.+.61*qc)
	      tve_v(i,j) = tve_v(i,j)*(1.+.61*qe)
	    endif
	  end do
        end do

!--------------------------------------------------------------------
!   determine whether the parcel temperature is cooler or warmer than 
!   the environment.
!--------------------------------------------------------------------
        call ieq_z(tvc_v, tve_v, ieqv, cape_skip)

!---------------------------------------------------------------------
!   add the contribution to column cape from this pressure layer.
!---------------------------------------------------------------------
        do j=jminp,jmaxp
          do i=iminp,imaxp
	    if (.not. cape_skip(i,j) ) then
	      if (ieqv(i,j) .ge. 0) then
                delt_v(i,j) = rair*(tvc_v(i,j)-tve_v(i,j))*    &
			     alog(p_v(i,j,k)/p_v(i,j,k+1))
	        xcape_v(i,j)=xcape_v(i,j)+delt_v(i,j)
              end if
            endif
	  end do
	end do

!---------------------------------------------------------------------
!   print out cape and cape contribution from this level.
!---------------------------------------------------------------------
!           do j=jminp,jmaxp
!	      if (debug_jt(j)) then
	      if (in_debug_window) then
	        if (.not. cape_skip(itest,jdebug) ) then
	          if (ieqv(itest,jdebug) .ge. 0) then
    		    print *,   'DONNER_DEEP/cape: k,delt,xcape= ',k,   &
				delt_v(itest,jdebug),    &
				xcape_v(itest,jdebug)
                  endif
                endif
              endif
!           end do

      end do  ! end of k loop

!--------------------------------------------------------------------
!  print out diagnostics (cape, cin, tot), if desired.
!--------------------------------------------------------------------
!       do j=jminp,jmaxp
!  if (debug_jt(j)) then
          if (in_debug_window) then
            if (.not. cape_skip(itest,jdebug) ) then
              tot(itest,jdebug) = xcape_v(itest,jdebug) -   &
			       cin_v(itest,jdebug)
              print *,   'DONNER_DEEP/cape: cin= ',cin_v(itest,jdebug),' J/kg'
              print *,   'DONNER_DEEP/cape: xcape= ',xcape_v(itest,jdebug),  &
					     ' J/kg'
              print *,   'DONNER_DEEP/cape: tot= ',tot(itest,jdebug),' J/kg'
            endif
          endif
!       end do

!--------------------------------------------------------------------
!  check for error in cape calculation. stop execution if present.
!--------------------------------------------------------------------
      error_flag = 0
      do j=jminp,jmaxp
        do i=iminp,imaxp
          if (.not. (cape_exit(i,j)) )  then
            if (xcape_v(i,j) .lt. 0.) error_flag = 1
            exit
          endif
        end do
      end do
      if (error_flag == 1)  &
        call error_mesg ( 'cape_vect',  &
	                  ' xcape error -- value < 0.0 ', FATAL)
!      call toc ('cape', '7')

!---------------------------------------------------------------------






end subroutine cape_vect




!#####################################################################

function iequ(x,y)

!
!     Checks for equality of two variable, within a tolerance eps.
!
!     On Input:
!
!        x    first variable
!        y    second variable
!
!     On Output:
!
!        equ  flag, equal to zero if x=y within eps
!                   equal to 10 if x greater than y
!                   equal to -10, if x less than y
!
      iequ=0
!     eps=1.e-14
      eps=1.e-13
      epsm=-eps
      d=x-y
      if (d .gt. eps) iequ=10
      if (d .lt. epsm) iequ=-10
!      return



end function iequ

subroutine  ieq_y (x, y, ieq)

!-------------------------------------------------------------------
real, dimension(:,:),   intent(in)   :: x
real,                   intent(in)   :: y
integer,dimension(:,:), intent(out)  :: ieq
!-------------------------------------------------------------------
!
!     Checks for equality of two variable, within a tolerance eps.
!
!     On Input:
!
!        x    first variable
!        y    second variable
!
!     On Output:
!
!        equ  flag, equal to zero if x=y within eps
!                   equal to 10 if x greater than y
!                   equal to -10, if x less than y
!
!      iequ=0
!!     eps=1.e-14
     real, dimension(size(x,1), size(x,2)) :: d

      eps=1.e-13
      epsm=-eps

      d(:,:) = x(:,:) - y

      do j=1,size(x,2)
      do i=1,size(x,1)

        if (d(i,j) .gt. eps) then
	 ieq(i,j) = 10
      else if (d (i,j).lt. epsm) then
	 ieq(i,j) = -10
      else 
	 ieq(i,j) = 0
      endif
      end do
      end do



!     if (d .gt. eps) iequ=10
!     if (d .lt. epsm) iequ=-10
!      return



end subroutine ieq_y



subroutine  ieq_z (x, y, ieq, flag)

!-------------------------------------------------------------------
real, dimension(:,:),   intent(in)   :: x, y
integer,dimension(:,:), intent(out)  :: ieq
logical, dimension(:,:), intent(in)  :: flag
!-------------------------------------------------------------------
!
!     Checks for equality of two variable, within a tolerance eps.
!
!     On Input:
!
!        x    first variable
!        y    second variable
!
!     On Output:
!
!        equ  flag, equal to zero if x=y within eps
!                   equal to 10 if x greater than y
!                   equal to -10, if x less than y
!
!      iequ=0
!!     eps=1.e-14
     real, dimension(size(x,1), size(x,2)) :: d

      eps=1.e-13
      epsm=-eps
      

!      d(:,:) = x(:,:) - y(:,:)

      do j=1,size(x,2)
      do i=1,size(x,1)

	if (.not. flag(i,j)) then
        d(i,j) = x(i,j) - y(i,j)
        if (d(i,j) .gt. eps) then
	 ieq(i,j) = 10
      else if (d (i,j).lt. epsm) then
	 ieq(i,j) = -10
      else 
	 ieq(i,j) = 0
      endif
      endif
      end do
      end do



end subroutine ieq_z






subroutine  ieq_x (x, y, ieq)

!
!     Checks for equality of two variable, within a tolerance eps.
!
!     On Input:
!
!        x    first variable
!        y    second variable
!
!     On Output:
!
!        equ  flag, equal to zero if x=y within eps
!                   equal to 10 if x greater than y
!                   equal to -10, if x less than y
!
!      iequ=0
!!     eps=1.e-14

      eps=1.e-13
      epsm=-eps
      d=x-y
      if (d .gt. eps) then
	 ieq = 10
      else if (d .lt. epsm) then
	 ieq = -10
      else 
	 ieq = 0
      endif
!     if (d .gt. eps) iequ=10
!     if (d .lt. epsm) iequ=-10
!      return



end subroutine ieq_x



!----------------------------------------------------------------------

subroutine cuclo_vect (dcape_v,         pcape_v, plfc_v, plzb_v,  &
   qli0_v, qli1_v, qr_v, qt_v,  &
	  r_v,     ri_v, rl_v, rpc_v,     tcape_v, tpc_v, a1_v,  &
		  exit_flag, nine_flag, jabs)

!---------------------------------------------------------------------
integer, dimension(:), intent(in)   :: jabs
real, dimension(:,:  ), intent(in) :: dcape_v, plfc_v, plzb_v
logical, dimension(:,:), intent(in) ::  exit_flag, nine_flag
real, dimension(:,:,:), intent(in) :: pcape_v, qli0_v, qli1_v,    &
				      qr_v, qt_v, &
				      r_v, ri_v, rl_v, rpc_v, tcape_v, &
				      tpc_v
real, dimension(:,:), intent(out)  :: a1_v
!---------------------------------------------------------------------

!
!     Calculates a_1(p_b) for closing cumulus parameterization.
!     See LJD notes, "Cu Closure D," 6/11/97
!
!
!     On Input:
!
!     Parameters:
! 
!        nc        number of levels at cloud-model resolution
!
!     Variables:
!
!        cpd       specific heat of dry air at constant pressure [J/(kg K)]
!        cpv       specific heat of water vapor
!        dcape     rate of change of convective available potential
!                  energy due to large-scale processes [J/(kg s)]
!        epsilo    ratio of molecular weights of water vapor to air
!        p         pressure (Pa)
!                  level 1 nearest ground.
!        plfc      pressure at level of free convection (Pa)
!        plzb      pressure at level of zero buoyancy (Pa)
!        qli0      normalized component of cumulus condensate forcing
!                  [kg(H2O)/(kg sec)]
!                  level 1 nearest ground.
!                  defined in "Cu Closure D," p. 4.
!        qli1      un-normalized component of condensate forcing
!                  [kg(h2O)/(kg sec)]
!                  level 1 nearest ground.
!                  defined in "Cu Closure D," p. 4.
!        qr        normalized cumulus moisture forcing [kg(H2O)/(kg sec)]
!                  level 1 nearest ground.
!                  defined in "Cu Closure D," p. 1.
!        qt        normalized cumulus thermal forcing (K/sec)
!                  level 1 nearest ground.
!                  defined in "Cu Closure D," p. 1.
!        r         large-scale water-vapor mixing ratio (kg(H2O)/kg)
!                  level 1 nearest ground.
!        rd        gas constant for dry air [J/(kg K)]
!        ri        large-scale ice mixing ratio (kg(H2O)/kg)
!                  level 1 nearest ground.
!        rl        large-scale liquid mixing ratio (kg(H2O)/kg)
!                  level 1 nearest ground.
!        rlat      latent heat of vaporization (J/kg)
!        rpc       parcel vapor mixing ratio from Cape.F [kg(H2O)/kg]
!                  Index 1 at bottom of model.
!        rv        gas constant for water vapor (J/(kg K))
!        t         large-scale temperature (K)
!                  level 1 nearest ground.
!        tpc       parcel temperature from Cape.F (K)
!                  Index 1 at bottom of model.
!   
!     On Output:
! 
!         a1       fractional area of index-1 cu subensemble
!        

      real, dimension (ncap) :: p, qli0, qli1, qr, qt, r, ri, rl, t, &
				rt, tpc, tpca, ta, ra, tden, tdena, &
				rpc, rpca, dtpdta

      real, dimension (size(tcape_v,1), size(tcape_v,2), ncap) :: &
		       ra_v, ta_v, rpca_v, tpca_v

   logical, dimension (size(tcape_v,1), size(tcape_v,2), ncap) :: &
		       exit_3d

      real, dimension (size(tcape_v,1), size(tcape_v,2) ) :: &
                     cin_v2, plcl_v2, plfc_v2, plzb_v2, xcape_v2

      logical  :: debug_ijt



!--------------------------------------------------------------------
!  define 3d execution flag.
!---------------------------------------------------------------------
      do k=1,ncap
        do j=jminp,jmaxp
          do i=iminp,imaxp
	    if ( (.not. exit_flag(i,j))   .and. &
	         ( .not. nine_flag(i,j)) ) then
              exit_3d(i,j,k) = .false.
            else
              exit_3d(i,j,k) = .true.
            endif
          end do
        end do
      end do

!--------------------------------------------------------------------
!  initialize perturbed parcel profile and perturbed parcel 
!  environmental profiles where calculation is to be done.
!--------------------------------------------------------------------
      do k=1,ncap
        do j=jminp,jmaxp
          do i=iminp,imaxp
	    if ( .not. exit_3d(i,j,k) ) then

!--------------------------------------------------------------------
!  initialize perturbed parcel profiles to actual parcel profiles.
!---------------------------------------------------------------------
   	      tpca_v(i,j,k) = tcape_v(i,j,k)
              rpca_v(i,j,k) = r_v(i,j,k)

!-------------------------------------------------------------------
!  define environmental profiles for perturbed parcel as the actual 
!  parcel soundings.
!-------------------------------------------------------------------
	      ra_v(i,j,k) = r_v(i,j,k)
	      ta_v(i,j,k) = tcape_v(i,j,k)

!--------------------------------------------------------------------
!  if calculation not needed in column, set temperature to ~0.0 to
!  prevent calculation in cape_vect subroutine.
!--------------------------------------------------------------------
            else
              ta_v(i,j,k) = 0.007
	      ra_v(i,j,k) = 0.0          
            endif
          end do
        end do
      end do

!--------------------------------------------------------------------
!   perturb surface mixing ratio and temperature to calculate 
!   derivative of parcel density temperature w.r.t. surface
!   large-scale density temperature.
!--------------------------------------------------------------------
      do j=jminp,jmaxp
        do i=iminp,imaxp
	  if ( (.not. exit_flag(i,j) )    .and. &
	       (.not. nine_flag(i,j) ) )  then
            ra_v(i,j,1) = ra_v(i,j,1) - 0.01*ra_v(i,j,1)
	    ra_v(i,j,1) = max(ra_v(i,j,1), 0.0)
!            if (ra_v(i,j,1) .lt. 0.) ra_v(i,j,1) = 0.
            ta_v(i,j,1) = tcape_v(i,j,1) - 1.0
          endif
        end do
      end do

!--------------------------------------------------------------------
!   call cape_vect to calculate t and r profiles for perturbed parcel.
!   only tpca_v and rpca_v are needed as output from this call.
!--------------------------------------------------------------------
      call cape_vect (pcape_v, ra_v, ta_v, cin_v2, plcl_v2, plfc_v2,  &
		      plzb_v2, rpca_v, tpca_v, xcape_v2, jabs)

!---------------------------------------------------------------------

      do j=jminp,jmaxp
      do i=iminp,imaxp

!        if (debug_jt(j) .and. i == itest) then
        if (in_debug_window .and. j == jdebug .and. i == itest) then
	  debug_ijt = .true.
        else
	  debug_ijt = .false.
        endif

	dcape = dcape_v(i,j)
	plfc = plfc_v(i,j)
	plzb = plzb_v(i,j)
	do k=1,ncap

	rpca(k) = rpca_v(i,j,k)
	tpca(k) = tpca_v(i,j,k)

	p(k) = pcape_v(i,j,k)
	qli0(k) = qli0_v(i,j,k)
	qli1(k) = qli1_v(i,j,k)
	qr(k) = qr_v(i,j,k)
	qt(k) = qt_v(i,j,k)
	 r(k) = r_v(i,j,k)
	 t(k) = tcape_v(i,j,k)
	 ri(k) = ri_v(i,j,k)
	 rl(k) = rl_v(i,j,k)
	 rpc(k) = rpc_v(i,j,k)
	 tpc(k) = tpc_v(i,j,k)
	 end do
	if ( (.not. exit_flag(i,j) ).and. &
	     ( .not. nine_flag(i,j) )) then

          if (debug_ijt) then
	    print *, 'DONNER_DEEP/cuclo: cpi,cpv,dcape=', cpi, cpv, dcape
	    print *, 'DONNER_DEEP/cuclo: rocp, rair, latvap,rvap=', rocp, &
					rair, latvap, rvap
	    do k=1,ncap
	    print *, 'DONNER_DEEP/cuclo: k,p,qr,qt,r,t  =', k, p(k), qr(k), &
					  qt(k), r(k), t(k)
            end do
	    do k=1,ncap
	    print *, 'DONNER_DEEP/cuclo: k,p,tpc, rpc   =', k, p(k), tpc(k),&
					              rpc(k)            
            end do
	    do k=1,ncap
	    print *, 'DONNER_DEEP/cuclo: k,p,qli0,qli1,ri,rl =', k, p(k),   &
				     qli0(k), qli1(k), ri(k), rl(k)
            end do
            print *, 'DONNER_DEEP/cuclo: plfc,plzb= ',plfc,plzb
         endif
!---------------------------------------------------------------------
!     calculate total-water mixing ratio. 
!---------------------------------------------------------------------
      do k=1,ncap
        rt(k)=r(k)+ri(k)+rl(k)
      end do



!
!     Calculate density temperatures. No liquid water in cape calculation.
!
      do k=1,ncap
	tden(k)=tpc(k)*(1.+(rpc(k)/rocp  )) 
	tdena(k)=tpca(k)*(1.+(rpca(k)/rocp  ))
      end do
      tdens=t(1)*(1.+(r(1)/rocp  ))
      tdensa=ta_v(i,j,1)*(1.+(ra_v(i,j,1)/rocp  ))

!
!     Evaluate derivative of parcel density temperature w.r.t. surface
!     large-scale density temperature.
!

      do k=1,ncap
	 dtpdta(k)=(tdena(k)-tden(k))/(tdensa-tdens)

! if (debug_ijt) then
! print *, 'DONNER_DEEP/cuclo: k,tden(k),tdena(k),dtpdta(k)= ',   &
!	   k,tden(k), tdena(k),dtpdta(k)
!        endif

!	 write(6,*) 'k,tden(k),tdena(k),dtpdta(k)= ',k,tden(k),
!     1    tdena(k),dtpdta(k)

      end do
       do k=1,ncap

          if (debug_ijt) then
             print *, 'DONNER_DEEP/cuclo: k,tden(k),tdena(k),dtpdta(k)= ',   &
 	                 k,tden(k), tdena(k),dtpdta(k)
           endif
	 if (debug_ijt) then
	 print *, 'DONNER_DEEP/cuclo: k,tpca,rpca= ', k, tpca(k),rpca(k)
         endif

!	 write(6,*) 'k,tpca,rpca= ',k,tpca(k),rpca(k)

       end do
!
!     Calculate I1 and I2 integrals from p. 5 of "Cu Closure D" notes.
!
      ri1=0.
      ri2=0.
      k=1
      rild=qt(k)*(rocp  +r(k))/(rocp  *(1.+rt(k)))
      rile=t(k)*(1.+rl(k)+ri(k)-rocp  )*qr(k)
      rile=rile/(rocp  *((1.+rt(k))**2))
      rilf=-t(k)*(rocp  +r(k))*qli0(k)
      rilf=rilf/(rocp  *((1.+rt(k))**2))
      ri2b=t(k)*(rocp  +r(k))/(rocp  *((1.+rt(k))**2))
      ri2b=ri2b*qli1(k)
      sum2=rild+rile+rilf


      do k=2,ncap
	 if (p(k) .eq. 0.) go to 3
	 rilak=-qt(k)*(rocp  +r(k))/(rocp  *(1.+rt(k)))
         rilbk=-t(k)*(1.+rl(k)+ri(k)-rocp  )*qr(k)
	 rilbk=rilbk/(rocp  *((1.+rt(k))**2))
	 rilck=t(k)*(rocp  +r(k))*qli0(k)
	 rilck=rilck/(rocp  *((1.+rt(k))**2))
	 rilakm=-qt(k-1)*(rocp  +r(k-1))/(rocp  *(1.+rt(k-1)))
         rilbkm=-t(k-1)*(1.+rl(k-1)+ri(k-1)-rocp  )*qr(k-1)
	 rilbkm=rilbkm/(rocp  *((1.+rt(k-1))**2))
	 rilckm=t(k-1)*(rocp  +r(k-1))*qli0(k-1)
	 rilckm=rilckm/(rocp  *((1.+rt(k-1))**2))
	 rila=.5*(rilak+rilakm)
	 rilb=.5*(rilbk+rilbkm)
	 rilc=.5*(rilck+rilckm)
	 ri2ak=t(k)*(rocp  +r(k))/(rocp  *((1.+rt(k))**2))
	 ri2ak=ri2ak*qli1(k)
	 ri2akm=t(k-1)*(rocp  +r(k-1))/(rocp  *((1.+rt(k-1))**2))
	 ri2akm=ri2akm*qli1(k-1)
	 ri2a=.5*(ri2ak+ri2akm)
	 sum1=rila+rilb+rilc
	 ri1=ri1+(alog(p(k-1)/p(k)))*(sum1+dtpdta(k)*sum2)
	 ri2=ri2+(alog(p(k-1)/p(k)))*(ri2a-dtpdta(k)*ri2b)
	 if (debug_ijt) then
	   print *, 'DONNER_DEEP/cuclo: k,dtpdta(k)= ',k,dtpdta(k)
 	   print *, 'DONNER_DEEP/cuclo: rila,rilb,rilc= ', rila,rilb,rilc
           print *, 'DONNER_DEEP/cuclo: ri1,ri2= ',ri1,ri2
           print *, 'DONNER_DEEP/cuclo: sum1,sum2= ',sum1,sum2
         endif
      end do

      if (debug_ijt) then
	print *, 'DONNER_DEEP/cuclo: rild,rile,rilf= ', rild, rile, rilf
      endif

 3    continue
      if (ri1 .ge. 0) then
	a1_v(i,j) = 0.
	cycle 
      end if
      ri1=rair*ri1
      ri2=rair*ri2
 2    continue
      a1_v(i,j)=-(ri2+dcape)/ri1


     endif
     end do
     end do



end  subroutine cuclo_vect



!#####################################################################

subroutine mulsub_vect (ampt_v, arat_v, pr_v, q_v, sfcqf_v, sfcsf_v,&
			t_v, amax_v, &
			cmui_v, cmus_v, cual_v,cuq_v, cuql_v, ecds_v,&
			 eces_v, emds_v, &
      			emei_v, emes_v, disa_v, disb_v, disc_v,    &
			disd_v, dise_v,  &
			dmeml_v, elt_v, fre_v, qmes_v, tmes_v,    &
			tpre_v, uceml_v,  &
                     umeml_v, wmms_v, wmps_v, jabs, exit_flag)

!---------------------------------------------------------------------
integer, dimension(:), intent(in)   :: jabs
real, dimension(:,:,:),intent(in)   :: arat_v
real, dimension(:,:), intent(in)    :: sfcqf_v, sfcsf_v
real, dimension(:,:), intent(out)   :: ampt_v, amax_v,   &
				       cmui_v, emei_v, tpre_v

real, dimension(:,:,:), intent(in)  ::  pr_v, q_v, t_v
real, dimension(:,:,:), intent(out) ::  cmus_v, cual_v, ecds_v, &
		 eces_v, emds_v, emes_v, disa_v, disb_v, disc_v, &
		 disd_v, dise_v, dmeml_v, elt_v, fre_v, qmes_v, &
		 tmes_v, uceml_v, umeml_v, wmms_v, wmps_v, cuq_v, &
                 cuql_v
logical, dimension(:,:), intent(in) :: exit_flag
!--------------------------------------------------------------------




!
!     Calculates thermal and moisture forcing by an ensemble of cumulus
!     elements, following Donner (1993, JAS). See also LJD notes, "Cu 
!     Closure A," 2/97. Normalized forcings by a_1(p_b).
!     L. Donner  GFDL 27 Apr 97
!
!     On Input:
!
!     arat(kpar)       a_i(p_b)/a_1(p_b) for each of kpar subensembles
!                      Subsensemble kpar should have least entrainment.
!     pr(nlev)         pressure at GCM resolution (Pa)
!                      index 1 at ground.
!     q(nlev)          mixing ratio at GCM resolution (kg(H2O)/kg)          
!                      index 1 at ground.
!     sfcqf            surface sensible heat flux (W/((m**2))
!                      set to zero if treated by large-scale eddy-fluxes
!     sfcsf            surface moisture flux (kg(H2O)/((m**2) sec))
!                      set to zero if treated by large-scale moisture
!                      fluxes
!     t(nlev)          temperature at GCM resolution (K)
!                      index 1 at ground.
!
!     On Input as Parameters:
!
!     kpar             number of cumulus subensembles
!     nlev             number of levels at coarse resolution
!     ncm              number of levels at cloud-model resolution
!
!     On Output:
!     
!     amax             maximum value for a_1(p_b)
!                      See "a Bounds 6/7/97" notes
!     ampt             mesoscale cloud fraction, normalized by a(1,p_b)
!     contot           ratio of convective to total precipitation
!     cmui             normalized vertical integral of mesoscale-updraft
!                      deposition (kg(H2O)/((m**2) sec)
!     cmus(nlev)       normalized mesoscale-updraft deposition
!                      (kg(H2O)/kg/sec)
!     cual(nlev)       cloud fraction, cells+meso, normalized by a(1,p_b)
!     cuq(nlev)        ice content in cells, weighted by cell area,
!                      (kg(H2O)/kg)
!                      index 1 at model bottom
!     cuqll(nlev)      liquid content in cells, weighted by cell area,
!                      (kg(H2O)/kg)
!                      index 1 at model bottom
!     ecds(nlev)       normalized convective downdraft evaporation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     eces(nlev)       normalzed convective-updraft evporation/sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emds(nlev)       normalized mesoscale-downdraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     emei             normalized vertical integral of mesoscale-updraft
!                      sublimation (kg(h2O)/((m**2) sec)
!     emes(nlev)       normalized mesoscale-updraft sublimation
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     disa(nlev)       normalized thermal forcing, cells+meso (K/sec)
!                      (excludes convergence of surface heat flux)
!                      index 1 at ground. Cumulus thermal forcing defined
!                      as in Fig. 3 of Donner (1993, JAS).
!     disb(nlev)       normalized cell entropy-flux convergence (K/sec)
!                      (excludes convergence of surface flux)
!                      index 1 at ground. Entropy-flux convergence divided
!                      by (p0/p)**(rd/cp).
!     disc(nlev)       normalized cell condensation/deposition
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     disd(nlev)       normalized cell moisture-flux convergence
!                      (excludes convergence of surface moisture flux)
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dise(nlev)       normalized moisture forcing, cells+meso (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     dmeml(nlev)      mass flux in mesoscale downdraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!     elt(nlev)        normalized melting (K/sec)
!                      index 1 at ground.
!     fre(nlev)        normalized freezing (K/sec)
!                      index 1 at ground.
!     qmes(nlev)       normalized mesoscale moisture-flux convergence
!                      (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     sfcq(nlev)       boundary-layer mixing-ratio tendency due to surface
!                      moisture flux (kg(H2O)/kg/sec)
!     sfch(nlev)       boundary-layer heating due to surface heat flux
!                      (K/sec)
!     tmes(nlev)       normalized mesoscale entropy-flux convergence
!                      (K/sec)
!                      Entropy-flux convergence is mesoscale component
!                      of second term in expression for cumulus thermal
!                      forcing in Fig. 3 of Donner (1993, JAS).
!                      index 1 at ground.
!     tpre             total normalized precipitation (mm/day)
!     uceml(nlev)      normalized mass fluxes in cell updrafts
!                      (kg/((m**2)*s) 
!     umeml(nlev)      mass flux in mesoscale updraft (kg/((m**2) s))
!                      (normalized by a(1,p_b)) (index 1 at atmosphere
!                      bottom)
!                      index 1 at ground.
!     wmms(nlev)       normalized mesoscale deposition of water vapor from
!                      cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!     wmps(nlev)       normalized mesoscale redistribution of water vapor
!                      from cells (kg(H2O)/kg/sec)
!                      index 1 at ground.
!--------------------------------------------------------------------

      logical, dimension (size(t_v,1), size(t_v,2)) :: exit_mulsub
      real   , dimension (size(t_v,1), size(t_v,2)) :: tb_v, pb_v, &
						       qb_v, pt_v, &
						       rr_v, ps_v, &
						       precip_v, &
						       conint_v, &
						       dint_v
      real   , dimension (size(t_v,1), size(t_v,2), size(t_v,3) ) :: &
                                  sig_v         
      real   , dimension (size(t_v,1), size(t_v,2), ncap ) :: &
			    tcc_v, wv_v, rcl_v, te_v, qe_v, dpf_v, &
			    dfr_v, flux_v, qlw_v
      real, dimension (ncap) :: rlsm, emsm, efchr, emfhr, rlhr, ctfhr, &
			        cmfhr, qllw, rsc, tcc, wv, te, qe, &
                                rcl, cuah, disp, dis, dpf, qlw, dfr, &
                                alp, flux, ucemh, cuql, cuqli
      real, dimension (size(t_v,3))   :: t, q, pr, h1, q1, sig, efc, &
					 em, rlh, ctf, cmf, fre_vk, &
					 elt_vk, frea, elta, cual_vk, &
					 cuml, evap, disf, disg, dish, &
					 disl, dism, disn, diso, cmu, &
					 ecd, ece, emd, eme, wmm, wmp, &
                                         thlr, qlr, enctf, encmf, enev,&
                                         fres, elts, tmes_vk, qmes_vk, &
                                         sfcq, sfch, dmeml_vk,   &
					 uceml_vk, umeml_vk, disg_sv, &
					 cuq_vk, cuqll_vk

      real, dimension (size(t_v,3)+1) :: phr                   
      real, dimension (kpar) :: cuto, preto, pbma, ptma

!     logical :: test2, constab, thetl, lmeso, debug_ijt
    logical :: test2, constab, thetl, lmeso, debug_ijt, in_debug_column
      real    ::         latsub, latice, rair_mul, cpair_mul, emdi




!---------------------------------------------------------------------
!      call tic ('mulsub', '1')

!           do j=jminp,jmaxp
!      if (debug_jt(j)) then
	      if (in_debug_window) then
                print *, 'DONNER_DEEP/mulsub: kttest,itest,jtest= ',kttest, &
                          itest, jtest
              endif
!           end do

!---------------------------------------------------------------------
!   define constants. 
!---------------------------------------------------------------------
       DP=-1000.
      RR=1000.

!     EPSILO=.622
!     rh2o=461.
      GRAVIT=9.80616
      cpair_mul=1004.64
      rair_mul=287.04
      CAPPA=rair_mul/cpair_mul
      LATICE=3.336E05
      LATSUB=latvap+LATICE
      tmel=273.15
      tfre=258.
      dfre= 10. 

       do j=jminp,jmaxp
       do i=iminp,imaxp

      tpre_v(i,j)=0.
      ampt_v(i,j)   = 0.
      amax_v(i,j)   = 0.
      contot_v(i,j) = 1.
      cmui_v(i,j)   = 0.
      emei_v(i,j)   = 0.
!
      DO k=1,nlev
      cual_v(i,j,k)=0.
      uceml_v(i,j,k)=0.
      umeml_v(i,j,k)=0.
      dmeml_v(i,j,k)=0.
      cmus_v(i,j,k)=0.
      emds_v(i,j,k)=0.
      emes_v(i,j,k)=0.
      wmms_v(i,j,k)=0.
      wmps_v(i,j,k)=0.
      fre_v(i,j,k)=0.
      elt_v(i,j,k)=0.
      tmes_v(i,j,k)=0.
      disd_v(i,j,k)=0.
      disa_v(i,j,k)=0.
      disb_v(i,j,k)=0.
      disc_v(i,j,k)=0.
      dise_v(i,j,k)=0.
      qmes_v(i,j,k)=0.
      ecds_v(i,j,k)=0.
      eces_v(i,j,k)=0.
!     SIG(k)=PR(k)/PS
     end do
       if (.not. exit_flag(i,j)) then

!	if (debug_jt(j) .and. i == itest) then
!	  debug_ijt =.true.
!	else
!	  debug_ijt = .false.
!	endif

        if (in_debug_window .and. j == jdebug .and. i == itest) then
	  in_debug_column = .true.
	  debug_ijt = .true.
        else
	  in_debug_column = .false.
	  debug_ijt = .false.
	endif



	exit_mulsub(i,j) = .false.

! lmeso = .true. for mesoscale, .false. for no mesoscale
!  initialize to .true.
!        lmeso = .true.

	pr(:)   = pr_v(i,j,:)
	q (:)   =  q_v(i,j,:)
	t (:)   =  t_v(i,j,:)
      
!---------------------------------------------------------------------
!   define constants. 
!---------------------------------------------------------------------

      PS=PR(1)
      CONSTAB=.FALSE.
!!! NOTE THIS IS NEEDED TO FIX BUG PRESENT IN LCL
      pb = 0.0
!     tpre_v(i,j)=0.
!     ampt_v(i,j)   = 0.
!     amax_v(i,j)   = 0.
!     contot_v(i,j) = 1.
!     cmui_v(i,j)   = 0.
!     emei_v(i,j)   = 0.
!
      DO 11 k=1,nlev
!     cual_v(i,j,k)=0.
!     uceml_v(i,j,k)=0.
!     umeml_v(i,j,k)=0.
!     dmeml_v(i,j,k)=0.
!     cmus_v(i,j,k)=0.
!     emds_v(i,j,k)=0.
!     emes_v(i,j,k)=0.
!     wmms_v(i,j,k)=0.
!     wmps_v(i,j,k)=0.
!     fre_v(i,j,k)=0.
!     elt_v(i,j,k)=0.
!     tmes_v(i,j,k)=0.
!     disd_v(i,j,k)=0.
!     disa_v(i,j,k)=0.
!     disb_v(i,j,k)=0.
!     disc_v(i,j,k)=0.
!     dise_v(i,j,k)=0.
!     qmes_v(i,j,k)=0.
!     ecds_v(i,j,k)=0.
!     eces_v(i,j,k)=0.
      SIG(k)=PR(k)/PS
 11   CONTINUE


             if (debug_ijt) then
!	do k=1,kmax-ktest_model+1
		do k=1,nlev-ktest_model+1
                  print *, 'DONNER_DEEP/mulsub: k,T,Q,P= ',k,T(k),Q(k),PR(k)
		end do
              endif








      CALL LCL(TB,CONSTAB,PB,QB,PS,SIG,T,Q)

      tb_v(i,j) = tb
      pb_v(i,j) = pb
      qb_v(i,j) = qb

              if (debug_ijt) then
                print *, 'DONNER_DEEP/mulsub: tb,pb,qb= ',TB,PB,QB
                if (constab) then
                  print *, 'DONNER_DEEP/mulsub: constab true'
                endif
              endif

      if (constab) then
	exit_mulsub(i,j) = .true.
      endif 
 

     endif
     end do
     end do

!      call toc ('mulsub', '1')
!      call tic ('mulsub', '2')


!     if (get_my_pe() == 3) then
!       print *, 'in section2'
!     endif


       do j=jminp,jmaxp
       do i=iminp,imaxp

       if (.not. exit_flag(i,j)) then

!if (debug_jt(j) .and. i == itest) then
     if (in_debug_window .and. i == itest .and. j == jdebug) then
	  debug_ijt =.true.
	else
	  debug_ijt = .false.
	endif

! lmeso = .true. for mesoscale, .false. for no mesoscale
!  initialize to .true.
        lmeso = .true.

	pr(:)   = pr_v(i,j,:)
	q (:)   =  q_v(i,j,:)
	t (:)   =  t_v(i,j,:)
      
!---------------------------------------------------------------------
!   define constants. 
!
!---------------------------------------------------------------------

      PS=PR(1)
       AL=.5
      THETL=.FALSE.
      TEST2=.false.
      PRETOT=0.
      APTSUM=0.
      CUTOT=0.
      CATOT=0.
      ca=0.
      APT=0.
!
      DO  k=1,nlev
      sfcq(k)=0.
      sfch(k)=0.
      cuml(k)=0.
      DISG(k)=0.
      DISN(k)=0.
      FRES(k)=0.
      frea(k)=0.
      elta(k)=0.
      ELTS(k)=0.
      WMM(k)=0.
      wmp(k)=0.
      ecd(k)=0.
      ece(k)=0.
      cmu(k)=0.
      emd(k)=0.
      eme(k)=0.
      ENEV(k)=0.
      ENCTF(k)=0.
      ENCMF(k)=0.
      SIG(k)=PR(k)/PS
      end do

      PHR(1)=PS
      PHR(nlev+1)=0.
!     do k=2,kmax
      do k=2,nlev
       PHR(k)=(PR(k)+PR(k-1))/2.
      end do

      DO k=1,ncap
      rcl(k)=0.
      EMSM(k)=0.
      RLSM(k)=0.
      cuah(k)=0.
      cuql(k)=0.
      cuqli(k)=0.
      ucemh(k)=0.
      flux(k)=0.
      end do

       if ( .not. exit_mulsub(i,j)) then

      do k=1,ncap
	 alp(k)=0.
      end do

      tb = tb_v(i,j)
      qb = qb_v(i,j)
      pb = pb_v(i,j)



      DO 31 KOU=1,KPAR
!     if (get_my_pe() == 3) then
!       print *, 'in      31 loop, kou = ', kou
!     endif
      DO 578 k=1,ncap
      EFCHR(k)=0.
      EMFHR(k)=0.
      RLHR(k)=0.
      CTFHR(k)=0.
      CMFHR(k)=0.
      QLLW(k)=0.
 578  RSC(k)=0.


	rr_v(i,j) = rr
	ps_v(i,j) = ps
	sig_v(i,j,:) = sig(:)

!     if (get_my_pe() == 3) then
!       print *, 'before cloudm'
!     endif
      call cloudm_vect(i, j, tb_v, pb_v, tcc_v, pt_v, wv_v,   &
!  rr_v,rcl_v, te_v,qe_v,ps_v,t_v,q_v, sig_v,rh2o, cappa, &
	  rr_v,rcl_v, te_v,qe_v,ps_v,t_v,q_v, sig_v,      cappa, &
	  rair_mul,  &
!   EPSILO,cpair_mul,GRAVIT,LATICE,KOU,TFRE,precip_v,conint_v,dpf_v,   &
     cpair_mul,GRAVIT,LATICE,KOU,dfre,TFRE,precip_v,conint_v,dpf_v,   &
!          cpair_mul,GRAVIT,LATICE,KOU,TFRE,precip_v,conint_v,dpf_v,   &
!          cpair_mul,       LATICE,KOU,TFRE,precip_v,conint_v,dpf_v,   &
       dfr_v,dint_v,flux_v,qlw_v, debug_ijt)

       pmel = pb_v(i,j)
       
       tcc(:) = tcc_v(i,j,:)
       pt = pt_v(i,j)
       wv(:) = wv_v(i,j,:)
       rcl (:) = rcl_v(i,j,:) 
       te(:) = te_v(i,j,:)
       qe(:) = qe_v(i,j,:)
       precip = precip_v(i,j)
       conint = conint_v(i,j)
       dpf(:) = dpf_v(i,j,:)
       dfr(:) = dfr_v(i,j,:)
       dint = dint_v(i,j)
       flux(:) = flux_v(i,j,:)
       qlw(:) = qlw_v(i,j,:)
!     if (get_my_pe() == 3) then
!       print *, 'after  cloudm'
!     endif
!
!     accumulate normalized convective cloud fraction and
!     mass flux in cell updrafts
!
      do k=1,ncap
	 cuah(k)=cuah(k)+arat_v(i,j,kou)*(rcl(k)/rr)**2
	 cfracl=1.
         cfraci=0.
         cfract=tfre-dfre
         if (tcc(k) .lt. cfract) then
           cfracl=0.
           cfraci=1.
         end if
         if ((tcc(k) .ge. cfract) .and. (tcc(k) .le. tfre)) then
           cfraci=(tfre-tcc(k))/dfre
           cfracl=(tcc(k)-tfre+dfre)/dfre
         end if
         cuql(k)=cuql(k)+arat_v(i,j,kou)*cfraci*qlw(k)*(rcl(k)/rr)**2
         cuqli(k)=cuqli(k)+arat_v(i,j,kou)*cfracl*qlw(k)*(rcl(k)/rr)**2
	 ucemh(k)=ucemh(k)+arat_v(i,j,kou)*flux(k)/(rr**2)
	 if ((kou .eq. kpar) .and. (cuah(k) .gt. 0.)) then
           cuql(k)=cuql(k)/cuah(k)
           cuqli(k)=cuqli(k)/cuah(k)
         end if

      if (debug_ijt) then
        print *, 'DONNER_DEEP/mulsub: kou,k,ucemh= ',kou,k,ucemh(k)
      endif

!	 write(6,*) 'kou,k,ucemh= ',kou,k,ucemh(k)

      end do
!
!     If cloud thickness less than pdeep_mc, de-activate mesoscale
!     circulation.
!
      pdeet=pb-pt
      if (pdeet .lt. pdeep_mc) lmeso=.false.
!
!
!     call establ(esc,tb)
      call lookup_es(tb, esc)
      qcc=epsilo*esc/(pb-esc)
!
      conint=conint/(rr**2)
      precip=precip/(rr**2)
      dint=dint/(rr**2)

      if (debug_ijt) then
        print *, 'DONNER_DEEP/mulsub: conint,precip,dint= ',conint,precip,dint
      endif


      do k=1,ncap

      if (debug_ijt) then
        print *, 'DONNER_DEEP/mulsub: k,dfr,dpr,rcl= ',k,dfr(k),dpf(k),rcl(k)
        print *, 'DONNER_DEEP/mulsub: cuah(k)= ',cuah(k)
      endif


	 dfr(k)=dfr(k)/(rr**2)
	 dpf(k)=dpf(k)/(rr**2)

      if (debug_ijt) then
        print *, 'DONNER_DEEP/mulsub: k,dfr,dpr,rcl= ',k,dfr(k),dpf(k),rcl(k)
      endif


      end do
!
      IF (KOU .EQ. 1) THEN
         PBMU=PB
         PTMU=PT
      END IF
!
!
      IF (PB .EQ. PT) THEN

       if (debug_ijt) then
         print *,  'DONNER_DEEP/mulsub: kou,pb,pt= ',kou,pb,pt
       endif


      go to 165
      END IF
      IF (PRECIP .EQ. 0.) THEN

	if (debug_ijt) then
	  print *, 'DONNER_DEEP/mulsub: PRECIP=0 AFTER CLOUD MODEL'
        endif


        go to 165
      END IF
      CONPRE=CONINT
      CPRE=PRECIP
      CU=CONPRE*86400.
      RC=CPRE*86400.
      PRETOT=PRETOT+RC*arat_v(i,j,kou)
      CUTOT=CUTOT+CU*arat_v(i,j,kou)
!
       if (debug_ijt) then
        print *, 'DONNER_DEEP/mulsub: CONPRE, CPRE= ', conpre, cpre, & 
						 ' KG/(M**2)/SEC'
        print *,  'DONNER_DEEP/mulsub: CONDENSATION PRE= ',CU,' MM/DAY'
        print *,  'DONNER_DEEP/mulsub: CLOUD MODEL PRE= ',RC,' MM/DAY'
       endif
 4    CONTINUE
 5    CONTINUE
          DISP(1)=PB
      DO 22 k=1,ncap-1
      P=PB+k*DP
      DISP(k+1)=P
      NCC=k
      IF (P .LT. PT) GO TO 21
  22  CONTINUE
 21   CONTINUE
      NCCM=NCC-1
!
!     Calculate quantities required for realizability check on
!     cloud fraction. See a bounds notes (7/6/97).
!
      do k=1,ncap
	 alp(k)=alp(k)+((arat_v(i,j,kou)*(rcl(k)**2))/rcl(1)**2)
      end do
!
!     CALCULATE CUMULUS THERMAL FORCING AND MOISTURE FORCING AT
!     CLOUD-MODEL RESOLUTION
!

    if (debug_ijt) then
       print *,  'DONNER_DEEP/mulsub: PB,PT= ',PB,PT
    endif

!     if (get_my_pe() == 3) then
!       print *, 'before 511 loop'
!     endif

      SUMEMF=0.
      SUMLHR=0.
      SUMEFC=0.
         SUMTHET=0.
      summel=0.
      DO 511 IT=1,ncap-1
!     if (get_my_pe() == 3) then
!print *, 'in 511 loop, it = ', it
!     endif
      IF (IT.EQ.1) THEN
         PL=PB
         DPP=DP/2.
      END IF
      IF (IT.GT. 1) THEN
         PL=PB+(IT-2)*DP
         DPP=DP
      END IF
      PH=PB+(IT  )*DP
      IF (PH .GE. PT) ITH=IT+1
      IF (PH .LT. PT) THEN
         ITH=IT
         PH=PT
         P=PT
      END IF
!     if (get_my_pe() == 3) then
!print *, 'pl, ph, p, pb, pt = ', pl, ph, p, pb, pt
!     endif
      IF (PL .LE. PT) GO TO 502
      IF (IT.EQ. 1) ITL=1
      IF (IT.NE. 1) ITL =IT-1

      if (debug_ijt) then
        print *,  'DONNER_DEEP/mulsub: IT,ITL,ITH= ',IT,ITL,ITH
        print *,  'DONNER_DEEP/mulsub: TCC = ',TCC(IT),TCC(ITL),TCC(ITH)
      endif


!     CALL ESTABL(ESH,TCC(ITH))
!     CALL ESTABL(ESL,TCC(ITL))
      CALL lookup_es(TCC(ITH), ESH)
      CALL lookup_es(TCC(ITL), ESL)
      targ=tcc(it)
!     call establ(es,targ)   
      call lookup_es(targ, es)   
      rh=epsilo*esh/(ph+(epsilo-1.)*esh)
      rl=epsilo*esl/(pl+(epsilo-1.)*esl)
      pit=pb+(it-1)*dp
      rsc(it)=epsilo*es/(pit+(epsilo-1.)*es)
      RGH=287.05*(1.+.608*RH)
      RGL=287.05*(1.+.608*RL)

      if (debug_ijt) then
        print *,  'DONNER_DEEP/mulsub: QE= ',QE(IT),QE(ITL),QE(ITH)
        print *,  'DONNER_DEEP/mulsub: TE= ',TE(IT),TE(ITL),TE(ITH)
      endif


      TVEH=TE(ITH)*(1.+.61*QE(ITH))
      TVCH=TCC(ITH)*(1.+.61*RH)
! COMPUTE DP/DZ (DPDZH,DPDZL) AS EQUATION (B3) IN DONNER (1987),
! J. ATM. SCI.,??,PP ???-????.
      DPDZH=-GRAVIT*PH*(AL *TVEH+TVCH)/(rair_mul*(1.+AL )*TVEH*TVCH)
! COMPUTE VERTICAL EDDY TRANSPORT OF T (EHFH,EHFL) USING EQUATION (3).
      EHFH=((rcl(ith)/rr)**2)*WV(ITH)*DPDZH*(TCC(ITH)-TE(ITH))
      TVEL=TE(ITL)*(1.+.61*QE(ITL))
      TVCL=TCC(ITL)*(1.+.61*RL)
      DPDZL=-GRAVIT*PL*(AL *TVEL+TVCL)/(rair_mul*(1.+AL )*TVEL*TVCL)
      EHFL=((rcl(itl)/rr)**2)*WV(ITL)*DPDZL*(TCC(ITL)-TE(ITL))
! COMPUTE THE EDDY FLUX CONVERGENCE (EHF) BY FLUX DIFFERENCING
! ACROSS THE LAYER.
      IF (IT .EQ. ITH) THEN
         EFCHR(IT+1)=EHFH/DP
         SUMEFC=SUMEFC+EFCHR(IT+1)*DP/2.
         PTT=PT+DP
         SUMTHET=SUMTHET+EFCHR(IT+1)*( (1.0E05/PTT)**CAPPA )*DP/2.

        if (debug_ijt) then 
           print *, 'DONNER_DEEP/mulsub: EFCHR(IT+1)= ',EFCHR(IT+1)
	endif


         EHFH=0.
      END IF
      EHF=(EHFL-EHFH)/(2.*DPP)
      TVE=TE(IT)*(1.+.609*QE(IT))
      TVC=TCC(IT)*(1.+.609*RSC(IT))
      IF (ITH .NE. IT) P=(PL+PH)/2.
      dtv=tvc-tve
      DPDZ=-GRAVIT*P*(AL *TVE +TVC )/(rair_mul*(1.+AL )*TVE*TVC)
! COMPUTE THE CUMULUS VERTICAL-FLUX CONVERGENCE OF ENTROPY (EXF)
! AS GIVEN IN UN-NUMBERED EQUATION ON P. 2163.
      EXF=rair_mul*WV(IT)*DPDZ* (TCC(IT)-TE(IT))*((rcl(it)/rr)**2)/   &
       (cpair_mul*p )
      EHF=EHF+EXF
      PI= (1.0E05/P)**CAPPA
      EFCHR(IT)       =EHF
      SUMTHET=SUMTHET+EFCHR(IT)* ( (1.0E05/P)**CAPPA )*DPP
      RLHR(IT)=-DPF(IT)
      if (tcc(it) .ge. tfre) then
      convrat=latvap/cpair_mul
      end if
      if (tcc(it) .lt. tfre) then
      convrat=LATSUB/cpair_mul
      end if
!     if (get_my_pe() == 3) then
!print *, 'it, ith, tcc(it), tcc(ith), tmel', it, ith, tcc(it), &
!	  tcc(ith), tmel
!     endif
      if ((tcc(it) .ge. tmel) .and. (tcc(ith) .le. tmel)) pmel=p
!
!     NO MESOSCALE-NEXT LINE
!
      if (.not. lmeso)    &
            qlw(it)=-dpf(it)*(1.-(rc/cu))
      IF (RLHR(IT) .LT. 0.)  THEN

      if (debug_ijt) then
	print *, 'DONNER_DEEP/mulsub: RLHR .LT. 0.'
      endif


	  go to 165
      END IF
      CTFHR(IT)=RLHR(IT)*convrat +EHF
      QLHR  =DPF(IT)
      IF (QLHR   .GT. 0.) QLHR  =0.
! COMPUTE EDDY FLUX CONVERGENCE OF MOISTURE AS IN EQUATION (3).
      EMFH=((rcl(ith)/rr)**2)*WV(ITH)*DPDZH*(RH-QE(ITH))
      EMFL=((rcl(itl)/rr)**2)*WV(ITL)*DPDZL*(RL-QE(ITL))
      IF (IT .EQ. 1) THEN
         THETF=EHFL
         THETF=THETF*( (1.0E05/PL)**CAPPA )
         EMFF=EMFL

       if (debug_ijt) then
          print *,  'DONNER_DEEP/mulsub: THETF,EMFF= ',THETF,EMFF
       endif


      END IF
      IF (IT .EQ. ITH) THEN
         EMFHR(IT+1)=EMFH/DP

	 if (debug_ijt) then
          print *,  'DONNER_DEEP/mulsub: EMFHRIT=ITH ',EMFHR(IT+1)
	 endif


         SUMEMF=SUMEMF+EMFHR(IT+1)*DP/2.
                       EMFH=0.
      END IF
      EMF=(EMFL-EMFH)/(2.*DPP)
      SUMEMF=SUMEMF+(EMF*(DPP  )   )
      SUMLHR=SUMLHR+(DPF(IT)*(DPP/GRAVIT) )

      if (debug_ijt) then
        print *,  'DONNER_DEEP/mulsub: IT,SUMLHR= ',IT,SUMLHR
      endif


      if (tcc(it) .le. tfre) then
      summel=summel+dpp*dpf(it)/gravit
      end if
      SUMEFC=SUMEFC+(EHF*(DPP  )   )

      if (debug_ijt) then
        print *,  'DONNER_DEEP/mulsub: SUMTHET,SUMEFC= ',SUMTHET,SUMEFC
      endif


      EMFHR(IT)=EMF
      CMFHR(IT)=QLHR  +EMF
      IF (IT .EQ. ITH) THEN
         RLHR(IT+1)=0.
         CMFHR(IT+1)=EMFHR(IT+1)
         CTFHR(IT+1)=EFCHR(IT+1)
      END IF

      if (debug_ijt) then
       IF (kou .eq. kpar) THEN
          print *,   'DONNER_DEEP/mulsub: IT,PL,PH= ',IT,PL,PH
          print *,   'DONNER_DEEP/mulsub: EHFH,EHFL,EXF,EHF= ',EHFH,EHFL,   &
							  EXF,EHF
         print *,   'DONNER_DEEP/mulsub: EMFH,EMFL,EMF= ',EMFH,EMFL,EMF
         print *,    'DONNER_DEEP/mulsub: WV,RH,QE= ',WV(ITH),RH,QE(ITH)
         print *,   'DONNER_DEEP/mulsub: RLHR,DPF= ',RLHR(IT),DPF(IT)
         print *,   'DONNER_DEEP/mulsub: WV,RL,QE= ',WV(ITL),RL,QE(ITL)
        END IF
      endif


      IF (P .EQ. PT) THEN
         CTFHR(IT)=0.
         CMFHR(IT)=0.
         APT=(rcl(it)/rcl(1))**2
      END IF
 511  CONTINUE
!     if (get_my_pe() == 3) then
!print *, 'after 511 loop  '
!     endif
 502  CONTINUE
!     if (get_my_pe() == 3) then
!print *, 'after 502     '
!     endif

      if (debug_ijt) then
         print *, 'DONNER_DEEP/mulsub: SUMLHR,SUMEMF= ',SUMLHR,SUMEMF
         print *, 'DONNER_DEEP/mulsub: SUMTHET=',SUMTHET
      endif


      SBL=0.
!

      if (debug_ijt) then
         print *, 'DONNER_DEEP/mulsub: VERAV DFR'
      endif


      CALL VERAV(DFR,PB,PT,THETL,FREA,SBL,PS,PR,PHR,CAPPA, debug_ijt)

      if (debug_ijt) then
        do jk=1,nlev
           print *,   'DONNER_DEEP/mulsub: jk,fres,frea= ',jk,fres(jk),frea(jk)
          end do
       print *,   'summel= ',summel
      endif

!     if (get_my_pe() == 3) then
!print *, 'before summel '
!     endif

! if (get_my_pe() == 3) then
! print *, 'summel, rc, cu', summel, rc, cu, dint
! endif
      summel=summel*rc/cu
      dints=dint+summel
!
!     No MESOSCALE, Next 6 Lines
!
! if (get_my_pe() == 3) then
! print *, 'lmeso', lmeso
! endif
      if (.not. lmeso) then 
! if (get_my_pe() == 3) then
! print *, 'dint ', dint 
! endif
      if (dint .ne. 0.) then
! if (get_my_pe() == 3) then
! print *, 'pmel,p1, pb', pmel, p1, pb
! endif
      if (pb > pmel) then
      p1=pmel
	  dmela=(summel+dint)*gravit/(p1-pb)
      dmela=-dmela*8.64e07
!     endif

!     if ( pb > p1 ) then
      call ver(dmela,pb,p1,pr,phr,elta, debug_ijt)
!     call ver(dmela,pb,p1,   phr,elta, debug_ijt)
      endif

      end if
      end if
	  summel=-summel*86400.

      if (debug_ijt) then
       print *,   'DONNER_DEEP/mulsub: summel,rc,cu= ',summel,rc,cu,' mm/day'
       print *,   'DONNER_DEEP/mulsub: VERAV QLW'
      endif

!     if (get_my_pe() == 3) then
!print *, 'before verav - qlw'
!     endif

      CALL VERAV(QLW,PB,PT,THETL,EVAP,SBL,PS,PR,PHR,CAPPA, debug_ijt)
      
      if (debug_ijt) then
       print *,   'DONNER_DEEP/mulsub: VERAV RLHR'
      endif

!     if (get_my_pe() == 3) then
!print *, 'before verav - rlhr'
!     endif

      CALL VERAV(RLHR,PB,PT,THETL,H1  ,SBL,PS,PR ,PHR ,CAPPA, debug_ijt)
!     if (get_my_pe() == 2) then
!if (i == 33 ) then
!print *, 'kou, rlhr', kou, (rlhr(k), k=1,100)
!print *, 'kou, h1', kou, (h1(k), k=1,40)
!endif
!     endif
!
!
      SSBL=0.
!
!     CALCULATE SUBCLOUD MOISTURE-FLUX CONVERGENCE (KG(H2O)/KG/SEC)
!
      SBL=(SSBL*GRAVIT+SUMEMF)/(PS-PB)
!
      if (debug_ijt) then
       print *,   'DONNER_DEEP/mulsub: VERAV EMFHR'
      endif

!     if (get_my_pe() == 3) then
!print *, 'before verav - emfhr'
!     endif

      CALL VERAV(EMFHR,PB,PT,   THETL,Q1,SBL,PS,PR,PHR,CAPPA, debug_ijt)

!     if (get_my_pe() == 3) then
!print *, 'after  verav - emfhr'
!     endif

      DO 517 JK=1,nlev
         fre_v(i,j,JK)=0.
         EM(JK)=Q1(JK)*8.64E07
	 disd_v(i,j,jk)=disd_v(i,j,jk)+em(jk)*arat_v(i,j,kou)
         CMF(JK)=(-H1(JK)+Q1(JK))*8.64E07
         if (t(jk) .ge. tfre) then
         convrat=latvap/cpair_mul
         end if
         if (t(jk) .lt. tfre) then
         convrat=latsub/cpair_mul
         end if
         RLH(JK)=H1(JK)*86400.*convrat
 517  CONTINUE
      THETL=.TRUE.
!
!
      SSBL=0.
!
!     CALCULATE SUBCLOUD ENTROPY-FLUX CONVERGENCE (K/S) DUE ONLY
!     TO SURFACE HEAT FLUX
!
      SBL=GRAVIT*SSBL/((PS-PB)*cpair_mul)
!

      if (debug_ijt) then
       print *,   'DONNER_DEEP/mulsub: VERAV EFCHR'
      endif

!     if (get_my_pe() == 3) then
!print *, 'before mesub  '
!     endif

      CALL VERAV(EFCHR,PB,PT,THETL,H1,SBL,PS,PR,PHR,CAPPA, debug_ijt)
      SBL=0.
      THETL=.FALSE.
      PBMA(KOU)=PB
      PTMA(KOU)=PT
      CUTO(KOU)=CU
      PRETO(KOU)=RC

!     if (get_my_pe() == 2) then
! if (i == 33) then
!print *, 'elt_v pre   mesub   ', (elt_v(33,1,k), k=1,40)
!print *, 'elta  pre   mesub   ', (elta (     k), k=1,40)
!endif
!endif

      if (lmeso) then
      elt_vk(:) = elt_v(i,j,:)
      fre_vk(:) = fre_v(i,j,:)
      CALL MESub(CU,RC,DINTS,PR,PHR,PS,PB,PT,GRAVIT,    &
!     CALL MESub(CU,RC,DINTS,PR,PHR,PS,PB,PT,           &
!       CAPPA,EPSILO,rair_mul,T,CA,                         &
        CAPPA,       rair_mul,T,CA,                         &
        tmel,ECD,ECE,fre_vk,elt_vk, debug_ijt)
	elt_v(i,j,:) = elt_vk(:)
	fre_v(i,j,:) = fre_vk(:)
      APTSUM=APTSUM+APT*arat_v(i,j,kou)
      else
      ca=0.
      do jk=1,nlev
	 ecd(jk)=0.
	 ece(jk)=0.
      end do
      end if
      CATOT=CATOT+CA*arat_v(i,j,kou)

!     if (get_my_pe() == 2) then
! if (i == 33) then
!print *, 'elt_v pre   518 loop', (elt_v(33,1,k), k=1,40)
!print *, 'elta  pre   518 loop', (elta (     k), k=1,40)
!endif
!endif



      DO 518 JK=1,nlev
!
!     MESOSCALE, NEXT 2 LINES
!
      if (lmeso) then
      FRES(JK)=FRES(JK)+fre_v(i,j,JK)*arat_v(i,j,kou)

      if (debug_ijt) then
        print *,  'DONNER_DEEP/mulsub: jk,fres,fre= ',jk,fres(jk),fre_v(i,j,jk)
      endif

            ELTS(JK)=ELTS(JK)+elt_v(i,j,JK)*arat_v(i,j,kou)
      end if
            FRES(JK)=FRES(JK)+FREA(JK)*arat_v(i,j,kou)

      if (debug_ijt) then
        print *,  'DONNER_DEEP/mulsub: jk,fres,frea= ',jk,fres(jk),frea(jk)
      endif


       elts(jk)=elts(jk)+arat_v(i,j,kou)*elta(jk)
         EFC(JK)=H1(JK)*86400.
         CTF(JK)=EFC(JK)+RLH(JK)
	 disb_v(i,j,jk)=disb_v(i,j,jk)+efc(jk)*arat_v(i,j,kou)
	 disc_v(i,j,jk)=disc_v(i,j,jk)+rlh(jk)*arat_v(i,j,kou)
         DISN(JK)=DISN(JK)+RLH(JK)
         ecds_v(i,j,JK)=ecds_v(i,j,JK)+ECD(JK)*arat_v(i,j,kou)
         eces_v(i,j,JK)=eces_v(i,j,JK)+ECE(JK)*arat_v(i,j,kou)
         ENCTF(JK)=ENCTF(JK)+arat_v(i,j,kou)*CTF(JK)
         ENCMF(JK)=ENCMF(JK)+arat_v(i,j,kou)*CMF(JK)
         ENEV(JK)=ENEV(JK)+arat_v(i,j,kou)*EVAP(JK)
         IF (DINT .EQ. 0.) DISG(JK)=DISG(JK)-arat_v(i,j,kou)*((ECD(JK) &
         +ECE(JK))*latvap/(cpair_mul*1000.))
         IF (DINT .NE. 0) THEN
   DISG(JK)=DISG(JK)-arat_v(i,j,kou)*(ECE(JK)*LATSUB/(cpair_mul*1000.))
   DISG(JK)=DISG(JK)-arat_v(i,j,kou)*(ECD(JK)*latvap/(cpair_mul*1000.))
         END IF
 518  CONTINUE
      DO 519 k=1,ncap
         RLSM(k)=(RLHR(k)*arat_v(i,j,kou))+RLSM(k)
         EMSM(k)=(EMFHR(k)*arat_v(i,j,kou))+EMSM(k)
      DIS(k)=PB+(k-1)*DP
 519  CONTINUE
      NDIA=1

      if (debug_ijt) then
       print *,   'DONNER_DEEP/mulsub: P & QLW'
      endif

!     if (get_my_pe() == 2) then
! if (i == 33) then
!print *, 'elt_v after 518 loop', (elt_v(33,1,k), k=1,40)
!endif
!endif

      DO 716 k=1,ncap-1
      if (debug_ijt) then
        print *, 'DONNER_DEEP/mulsub: k, P & QLW', k, DIS(k), QLW(k)
      endif


      IF (WV(k+1) .LE. 0.) GO TO 707
      NDIA=NDIA+1
 716  CONTINUE
 707  CONTINUE
      if (test2) then
!       CALL AGSETF('Y/ORDER.',1.)
      imi=2
!      CALL ANOTAT('M/S','PA',0,0,-1,0)
!      CALL EZXY(WV,DIS,NDIA,'CLOUD VERT VEL')
       end if
!      CALL EZXY(FLUX,DISP,NCC,'CU FLUX')
!      CALL AGSETF('Y/ORDER.',1.)


       if (debug_ijt) then
	 do k=1,nlev
	 print *, 'DONNER_DEEP/mulsub: k, p & ctf', k, pr(k), ctf(k)
	 end do
	 do k=1,nlev
	 print *, 'DONNER_DEEP/mulsub: k, p & cmf', k, pr(k), cmf(k)
	 end do
       endif


       if (debug_ijt) then
         print *,  'DONNER_DEEP/mulsub: kou,pb,pt= ',kou,pb,pt
       endif

 31   CONTINUE

     
!     if (get_my_pe() == 3) then
!       print *, 'after   31 loop'
!     endif


!
!     Select pt so it refers to the most penetrative sub-ensemble.
!     This is frequently, but not always, the sub-ensemble with
!     the lowest entrainment.
!
      pt=min(ptma(1),ptma(2),ptma(3),ptma(4),ptma(5),ptma(6),ptma(7))
!
!     Convert normalized convective cloud fraction and mass fluxes 
!     to GCM resolution.
!
      sbl=0.

      cual_vk(:) = cual_v(i,j,:)
     call verav(cuah,pb,pt,thetl,cual_vk,sbl,ps,pr,phr,cappa, debug_ijt)
      cual_v(i,j,:) = cual_vk(:)
      call verav(cuql,pb,pt,thetl,cuq_vk,sbl,ps,pr,phr,cappa,debug_ijt)
      cuq_v(i,j,:) = cuq_vk(:)
      call verav(cuqli,pb,pt,thetl,cuqll_vk,sbl,ps,pr,phr,cappa,debug_ijt)
      cuql_v(i,j,:) = cuqll_vk(:)

      uceml_vk(:) = uceml_v(i,j,:)
   call verav(ucemh,pb,pt,thetl,uceml_vk,sbl,ps,pr,phr,cappa, debug_ijt)
      uceml_v(i,j,:) = uceml_vk(:)

      apt=aptsum
      ampt_v(i,j)=5.*apt
!
!     Calculate quantities required for realizability check on
!     cloud fraction. See a bounds notes (7/6/97).
!
      al=maxval(alp)
      aal=0.
      aalm=0.
      if (al .gt. 0.) then
      aal=1./al
      aalm=1./(al+ampt_v(i,j))
      end if
      if (lmeso) then
         amax_v(i,j)=aalm
      else
	 amax_v(i,j)=aal
      end if

      if (debug_ijt) then
        print *,  'DONNER_DEEP/mulsub: CUTOT=',CUTOT,' PRETOT=',PRETOT
         print *, 'DONNER_DEEP/mulsub: CATOT=',CATOT
      endif


      CATOT=-CATOT
      PB=PBMA(1)

      if (debug_ijt) then     
       print *, 'DONNER_DEEP/mulsub: ps,pb,pt,lmeso= ',ps,pb,pt,lmeso
      endif


      if (lmeso) then

        if (debug_ijt) then
         print *, 'DONNER_DEEP/mulsub: ampt= ',ampt_v(i,j)
	endif


      dmeml_vk(:) = dmeml_v(i,j,:)
      umeml_vk(:) = umeml_v(i,j,:)
      elt_vk(:) = elt_v(i,j,:)
      qmes_vk(:) = qmes_v(i,j,:)
      tmes_vk(:) = tmes_v(i,j,:)
      CALL MEens(CUTOT,PRETOT,PR,PHR,PS,PB,PT,GRAVIT,       &
!       CAPPA,EPSILO,rair_mul,cpair_mul,LATVAP,rh2o,T,Q,APT,RLSM,   &
!        CAPPA,EPSILO,rair_mul,cpair_mul,latvap, rh2o,T,Q,APT,RLSM,   &
         CAPPA,EPSILO,rair_mul,cpair_mul,latvap,      T,Q,APT,RLSM,   &
!        CAPPA,       rair_mul,cpair_mul,latvap,      T,Q,APT,RLSM,   &
        EMSM,   &
!       catot,tmel,cuml,CMU,cmui_v(i,j),dmeml_vk,EMD,EME,emei_v(i,j), &
!       catot,tmel,cuml,CMU,cmuxxx     ,dmeml_vk,EMD,EME,emexxx     , &
	catot,tmel,cuml,CMU,cmuxxx     ,dmeml_vk,EMD,emdi,EME,emexxx, &
	  WMM,WMP,  &
!       elt_vk,contot_v(i,j),tmes_vk,qmes_vk,umeml_vk, debug_ijt)
        elt_vk,contotxx     ,tmes_vk,qmes_vk,umeml_vk, debug_ijt)

	cmui_v(i,j) = cmuxxx
	emei_v(i,j) = emexxx
	contot_v(i,j) = contotxx
	emdi_v(i,j) = emdi
	dmeml_v(i,j,:) = dmeml_vk(:)
	umeml_v(i,j,:) = umeml_vk(:)
	elt_v(i,j,:) = elt_vk(:)
	qmes_v(i,j,:) = qmes_vk(:)
	tmes_v(i,j,:) = tmes_vk(:)

       end if
	 if (pretot .ne. 0.) then
	    if (lmeso) tpre_v(i,j)=pretot/contot_v(i,j)
	    if (.not. lmeso) tpre_v(i,j)=pretot
         end if
         sumwmp=0.

!     if (get_my_pe() == 2) then
! if (i == 33) then
!print *, 'elt_v after meens', (elt_v(33,1,k), k=1,40)
!endif
!endif

         DO 132 JK=1,nlev

	 if (debug_ijt) then
            print *,'DONNER_DEEP/mulsub: jk, pr,cual,cuml= ',jk,pr(jk),  &
		    cual_v(i,j,jk), cuml(jk)
	 endif


            cual_v(i,j,jk)=cual_v(i,j,jk)+cuml(jk)

          if (debug_ijt) then
        print *, 'DONNER_DEEP/mulsub: jk,cuml,cual= ',jk,cuml(jk),   &
					      cual_v(i,j,jk)
        print *, 'DONNER_DEEP/mulsub: jk,cmu,emd,eme= ',jk,cmu(jk),emd(jk), &
		      eme(jk)
       print *, 'DONNER_DEEP/mulsub: jk,wmm,wmp,elt= ',jk,wmm(jk),wmp(jk),  &
			       elt_v(i,j,jk)
       print *, 'DONNER_DEEP/mulsub: jk,tmes,qmes= ',jk,tmes_v(i,j,jk),  &
			    qmes_v(i,j,jk)
	  endif


            cmus_v(i,j,JK)=CMU(JK)-WMM(JK)
            emds_v(i,j,JK)=EMD(JK)
            emes_v(i,j,JK)=EME(JK)
            wmms_v(i,j,JK)=wmms_v(i,j,JK)+WMM(JK)
            wmps_v(i,j,JK)=wmps_v(i,j,JK)+WMP(JK)
            sumwmp=sumwmp+wmps_v(i,j,jk)*(phr(jk)-phr(jk+1))
       if ((phr(jk+1) .le. pb) .and.       &
        (phr(jk) .ge. pb)) psmx=phr(jk+1)
!
!      MESOSCALE, Next Line
!
      if (lmeso) then
            ELTS(JK)=ELTS(JK)+elt_v(i,j,JK)
      end if
 132      CONTINUE
          sumwmp=sumwmp/(gravit*1000.)

	  if (debug_ijt) then
            print *,  'DONNER_DEEP/mulsub:sumwmp= ',sumwmp,' mm/day'
	  endif


 131  CONTINUE
 134  CONTINUE

	  if (debug_ijt) then
          print *, 'DONNER_DEEP/mulsub: CATOT= ',CATOT,' contot=',contot_v(i,j)
	  endif


!
!     surface heat flux (W/(m**2))
!
!     if (get_my_pe() == 3) then
!print *, 'ps, psmx', ps, psmx
!     endif
      ssbl=sfcsf_v(i,j)
      sbl=gravit*ssbl/((ps-psmx)*cpair_mul)
      if (ps > psmx) then
      call ver(sbl,ps,psmx,pr,phr,sfch, debug_ijt)
!     call ver(sbl,ps,psmx,   phr,sfch, debug_ijt)
      endif

!
!     surface moisture flux (kg(H2O)/((m**2) sec)
!
      ssbl=sfcqf_v(i,j)
      sbl=(ssbl*gravit)/(ps-psmx)
      if (ps > psmx) then
      call ver(sbl,ps,psmx,pr,phr,sfcq, debug_ijt)
!     call ver(sbl,ps,psmx,   phr,sfcq, debug_ijt)
      endif
      ESUMB=0.
      ESUMC=0.
      SUMF=0.
      SUMM=0.
      sumqme=0.
      DO 513 JK=1,nlev
      EVAP(JK)=ENEV(JK)*8.64E07
      DISH(JK)=-cmus_v(i,j,JK)
      DISF(JK)=DISH(JK)+ecds_v(i,j,JK)+eces_v(i,j,JK)+wmps_v(i,j,JK)
      disf(jk)=disf(jk)+emes_v(i,j,jk)+emds_v(i,j,jk)
      disf(jk)=disf(jk)+qmes_v(i,j,jk)
      CMF(JK)=ENCMF(JK)
!
      DISM(JK)=CMF(JK)
!      DISM(JK)=DISM(JK)+ecds_v(i,j,JK)+eces_v(i,j,JK)
!
!     MESOSCALE, NEXT LINE
!
      if (lmeso)    &
            CMF(JK)=CMF(JK)+DISF(JK)
!
!     NO MESOSCALE, NEXT 2 LINES
!
      if (.not. lmeso) then
      CMF(JK)=CMF(JK)+EVAP(JK)
      DISF(JK)=EVAP(JK)
      end if
!
      SUMF=SUMF+DISF(JK)*(PHR(JK)-PHR(JK+1))
      sumqme=sumqme+qmes_v(i,j,jk)*(phr(jk)-phr(jk+1))
      SUMM=SUMM+DISM(JK)*(PHR(JK)-PHR(JK+1))
      ESUMB=ESUMB+CMF(JK)*(PHR(JK)-PHR(JK+1))
 513  CONTINUE
!
      ESUMB=ESUMB/(GRAVIT*1000.)
      SUMF=SUMF/(GRAVIT*1000.)
      sumqme=sumqme/(gravit*1000.)
      SUMM=SUMM/(GRAVIT*1000.)

      if (debug_ijt) then
        print *,  'DONNER_DEEP/mulsub: SUMF= ',SUMF, ' MM/DAY'
        print *,  'DONNER_DEEP/mulsub: SUMM= ',SUMM, ' MM/DAY'
        print *,  'DONNER_DEEP/mulsub: ESUMB=',ESUMB, ' MM/DAY'
        print *,  'DONNER_DEEP/mulsub: sumqme= ',sumqme, ' mm/day'
      endif


      SUMG=0.
      SUMN=0.
	  sumelt=0.
	  sumfre=0.
	  summes=0.
      DO 514 JK=1,nlev

      if (debug_ijt) then
         print *, 'DONNER_DEEP/mulsub: JK,H1= ',JK,H1(JK)
      endif



      CTF(JK)=ENCTF(JK)
!
      ESUMC=ESUMC+CTF(JK)*(PHR(JK)-PHR(JK+1))
      DISL(JK)=cmus_v(i,j,JK)*LATSUB/(cpair_mul*1000.)
	  sumelt=sumelt+elts(jk)*(phr(jk)-phr(jk+1))
	  sumfre=sumfre+fres(jk)*(phr(jk)-phr(jk+1))
      fre_v(i,j,JK)=FRES(JK)*LATICE/(cpair_mul*1000.)

      if (debug_ijt) then
          print *, 'DONNER_DEEP/mulsub: jk,fres= ', jk, fres(jk)
      endif


      elt_v(i,j,JK)=-ELTS(JK)*LATICE/(cpair_mul*1000.)
      DISN(JK)=CTF(JK)+DISG(JK)+fre_v(i,j,JK)+elt_v(i,j,JK)
      disga=((emes_v(i,j,jk)+emds_v(i,j,jk))*latsub/(cpair_mul    &
         *1000.))

      if (debug_ijt) then
 	  print *, 'DONNER_DEEP/mulsub: jk,pr,emds,disga= ',jk,pr(jk),   &
		    emds_v(i,j,jk), disga
      endif


       disg_sv(jk) = disg(jk)
      DISG(JK)=DISL(JK)+DISG(JK)+fre_v(i,j,JK)+elt_v(i,j,jk)-disga+   &
        tmes_v(i,j,jk)
      if (t(jk) .ge. tfre) then
      DISO(JK)=EVAP(JK)*latvap/(cpair_mul*1000.)
      end if
      if (t(jk) .le. tfre) then
      DISO(JK)=EVAP(JK)*LATSUB/(cpair_mul*1000.)
      end if
!
!     MESOSCALE, NEXT LINE
!
      if (lmeso)     &
            CTF(JK)=CTF(JK)+DISG(JK)
!
!     NO MESOSCALE, NEXT 2 LINES
!
       if (.not. lmeso) then
      CTF(JK)=CTF(JK)-DISO(JK)
      DISG(JK)=-DISO(JK)
      end if
!
      disa_v(i,j,JK)=CTF(JK)
      dise_v(i,j,JK)=CMF(JK)
!
      SUMG=SUMG+DISG(JK)*(PHR(JK)-PHR(JK+1))
      SUMN=SUMN+DISN(JK)*(PHR(JK)-PHR(JK+1))
      summes=summes+tmes_v(i,j,jk)*(phr(jk)-phr(jk+1))
 514  CONTINUE

!     if (get_my_pe() == 2) then
! if (i == 33) then
!print *, 'disa_v in mulsub', (disa_v(33,1,k), k=1,40)
!print *, 'lmeso', lmeso
!print *, 'ctf    in mulsub', (ctf   (     k), k=1,40)
!print *, 'enctf  in mulsub', (enctf (     k), k=1,40)
!print *, 'diso   in mulsub', (diso  (     k), k=1,40)
!print *, 'disg   in mulsub', (disg  (     k), k=1,40)
!print *, 'disg_svin mulsub', (disg_sv(     k), k=1,40)
!print *, 'disl   in mulsub', (disl  (     k), k=1,40)
!print *, 'disga  in mulsub',  disga                  
!print *, 'fre_v  in mulsub', (fre_v (33,1,k), k=1,40)
!print *, 'elt_v  in mulsub', (elt_v (33,1,k), k=1,40)
!print *, 'tmes_v in mulsub', (tmes_v(33,1,k), k=1,40)
!  endif
!      endif
	  sumelt=sumelt/(gravit*1000.)
	  sumfre=sumfre/(gravit*1000.)
	  


      ESUMC=(ESUMC*cpair_mul)/(GRAVIT*latvap)


      summes=summes*cpair_mul/(gravit*latvap)


      SUMG=SUMG*cpair_mul/(GRAVIT*latvap)
      SUMN=SUMN*cpair_mul/(GRAVIT*latvap)


	  if (debug_ijt) then
 	    print *, 'DONNER_DEEP/mulsub: sumelt,sumfre= ',sumelt,sumfre,  &
							    ' mm/day'
            print *, 'DONNER_DEEP/mulsub: ESUMC=',ESUMC
            print *, 'DONNER_DEEP/mulsub: summes=',summes,' mm/day'
            print *, 'DONNER_DEEP/mulsub: SUMG=',SUMG,' MM/DAY'
            print *, 'DONNER_DEEP/mulsub: SUMN= ',SUMN,' MM/DAY'
	  endif

      ESUM=0.
      SUMEV=0.
      ESUMA=0.
      DO 516 k=1,nlev
      ESUMA=ESUMA+CTF(k)*(PHR(k)-PHR(k+1))
      ESUM=ESUM+CMF(k)*(PHR(k)-PHR(k+1))
      SUMEV=SUMEV+EVAP(k)*(PHR(k)-PHR(k+1))
 516  CONTINUE
      ESUM=ESUM/(GRAVIT*1000.)
      SUMEV=SUMEV/(GRAVIT*1000.)
      ESUMA=(ESUMA*cpair_mul)/(GRAVIT*latvap)

      if (debug_ijt) then
       print *,   'DONNER_DEEP/mulsub: ESUM= ',ESUM,' ESUMA=',ESUMA,    &
		   '  SUMEV=',SUMEV
       do k=1,nlev
	 if (disc_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & LHR =',  k, PR(k),disc_v(i,j,k)
	 endif
       end do
       do k=1,nlev
	 if (disb_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & EFC =',  k, PR(k),disb_v(i,j,k)
	 endif
       end do
       do k=1,nlev
	 if (disd_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & EMF =',  k, PR(k),disd_v(i,j,k)
	 endif
       end do
       do k=1,nlev
	 if (disn(k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & cell thermal forcing =',    &
				      k, PR(k),disn(k)
	 endif
       end do
       do k=1,nlev
	 if (dism(k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & cell moisture forcing =',    &
				      k, PR(k),dism(k)
	 endif
       end do
       do k=1,nlev
	 if ( fre_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & meso up freeze        =',    &
				      k, PR(k),fre_v(i,j,k)
	 endif
       end do
       do k=1,nlev
	 if ( elt_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & meso down melt        =',    &
				      k, PR(k),elt_v(i,j,k)
	 endif
       end do
       do k=1,ncc 
         print *, 'DONNER_DEEP/mulsub: k, P & cond/efc              =',    &
				      k, dis(k),ctfhr(k)
       IF (WV(k+1) .LE. 0.) exit      
       end do
       do k=1,nlev
         print *, 'DONNER_DEEP/mulsub: k, P & cond/mfc              =',    &
				      k, dis(k),cmfhr(k)
       IF (WV(k+1) .LE. 0.) exit      
       end do
       do k=1,nlev
	 if (cmus_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & meso up con           =',    &
				      k, PR(k),cmus_v(i,j,k)
         endif
       end do
       do k=1,nlev
	 if (emds_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & meso down evap        =',    &
				      k, PR(k),emds_v(i,j,k)
         endif
       end do
       do k=1,nlev
	 if (emes_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & meso up evap        =',    &
				      k, PR(k),emes_v(i,j,k)
         endif
       end do
       do k=1,nlev
	 if (wmms_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & meso cell con       =',    &
				      k, PR(k),wmms_v(i,j,k)
         endif
       end do
       do k=1,nlev
	 if (wmps_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & meso vap redist     =',    &
				      k, PR(k),wmps_v(i,j,k)
         endif
       end do
       do k=1,nlev
	 if (tmes_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & meso efc            =',    &
				      k, PR(k),tmes_v(i,j,k)
         endif
       end do
       do k=1,nlev
	 if (tmes_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & meso mfc            =',    &
				      k, PR(k),qmes_v(i,j,k)
         endif
       end do
       do k=1,nlev
	 if (eces_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & up con evap         =',    &
				      k, PR(k),eces_v(i,j,k)
         endif
       end do
       do k=1,nlev
	 if (ecds_v(i,j,k) /= 0.00 ) then
         print *, 'DONNER_DEEP/mulsub: k, P & down con evap         =',    &
				      k, PR(k),ecds_v(i,j,k)
         endif
       end do
      endif











!      GO TO 135
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISB,SIG ,NLEV,'ENTROPY FLUX CONVERGENCE')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISD,SIG ,NLEV,'MOISTURE FLUX CONVERGENCE')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISC,SIG,NLEV,'LATENT HEAT RELEASE')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(FRE ,SIG,NLEV,'FREEZE ')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(ELT, SIG,NLEV,'MELT ')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISA ,SIG,NLEV,'CUMULUS THERMAL FORCING')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISE, SIG,NLEV,'CUMULUS MOISTURE FORCING')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(CMUS, SIG,NLEV,'MESOSCALE UPDRAFT CONDENSATION')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(ECDS, SIG,NLEV,'CONVECTIVE DOWNDRAFT EVAPORATION')

 135  CONTINUE



!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(ECES, SIG,NLEV,'CONVECTIVE UPDRAFT EVAPORATION')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(EMDS, SIG,NLEV,'MESOSCALE DOWNDRAFT EVAPORATION')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(EMES, SIG,NLEV,'MESOSCALE UPDRAFT EVAPORATION')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(WMMS, SIG,NLEV,'MESOSCALE CELL CONDENS')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(WMPS, SIG,NLEV,'MESOSCALE VAPOR Redistribution')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISH, SIG,NLEV,'MESOSCALE MOISTURE TENDENCY')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISL, SIG,NLEV,'MESOSCALE THERMAL TENDENCY')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISF, SIG,NLEV,'MOISTURE MODIFICATIONS ')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISG, SIG,NLEV,'THERMAL MODIFICATIONS ')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISM, SIG,NLEV,'CONVECTIVE MOISTURE FORCING')
!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(DISN, SIG,NLEV,'CONVECTIVE THERMAL FORCING')
!      CALL ANOTAT('K/SEC','KPA',0,0,-1,0)
!      CALL EZXY(CTFHR,DISP,NCCM,'COND/EFC')
!      CALL ANOTAT('K/SEC','KPA',0,0,-1,0)
!      CALL EZXY(RLHR,DISP,NCCM,'RLHR')
!      CALL ANOTAT('K/SEC','KPA',0,0,-1,0)
!      CALL EZXY(EFCHR,DISP,NCCM,'EFCHR')
!      CALL ANOTAT('KG/KG/SEC','KPA',0,0,-1,0)
!      CALL EZXY(CMFHR,DISP,NCCM,'COND/MFC')


       if (debug_ijt) then
	  do k=1,nlev
	    if (ctf(k) /= 0.0) then
	    print *, 'DONNER_DEEP/mulsub: k, p & ens thermal forc', &
			k, pr(k), ctf(k)
            endif
	  end do
	  do k=1,nlev
	    if (cmf(k) /= 0.0) then
	    print *, 'DONNER_DEEP/mulsub: k, p & ens moisture forc', &
			k, pr(k), cmf(k)
            endif
	  end do
	  do k=1,nlev
	    if (disg(k) /= 0.0) then
	    print *, 'DONNER_DEEP/mulsub: k, p & thermal modifications', &
			k, pr(k), disg(k)
            endif
	  end do
	  do k=1,nlev
	    if (disf(k) /= 0.0) then
	    print *, 'DONNER_DEEP/mulsub: k, p & moisture modifications', &
			k, pr(k), disf(k)
            endif
	  end do
       endif





!      CALL ANOTAT('K/DAY','SI',0,0,-1,0)
!      CALL EZXY(tmes,SIG ,NLEV,'Meso EFC')
!      CALL ANOTAT('G/KG/DAY','SI',0,0,-1,0)
!      CALL EZXY(qmes,SIG ,NLEV,'Meso MFC')
! 

!     CLOSE GKS GRAPHICS.
!
!      call frame
!	  CALL CLSGKS
!

!     Convert to MKS units.
!
       do k=1,nlev
           disa_v(i,j,k)=disa_v(i,j,k)/86400.
	   disb_v(i,j,k)=disb_v(i,j,k)/86400.
	   disc_v(i,j,k)=disc_v(i,j,k)/86400.
	   disd_v(i,j,k)=disd_v(i,j,k)/8.64e07
	   dise_v(i,j,k)=dise_v(i,j,k)/8.64e07
	   fre_v(i,j,k)=fre_v(i,j,k)/86400.
	   elt_v(i,j,k)=elt_v(i,j,k)/86400.
	   cmus_v(i,j,k)=cmus_v(i,j,k)/8.64e07
	   ecds_v(i,j,k)=ecds_v(i,j,k)/8.64e07
	   eces_v(i,j,k)=eces_v(i,j,k)/8.64e07
	   emds_v(i,j,k)=emds_v(i,j,k)/8.64e07
	   emes_v(i,j,k)=emes_v(i,j,k)/8.64e07
	   wmms_v(i,j,k)=wmms_v(i,j,k)/8.64e07
	   wmps_v(i,j,k)=wmps_v(i,j,k)/8.64e07
	   tmes_v(i,j,k)=tmes_v(i,j,k)/86400.
	   qmes_v(i,j,k)=qmes_v(i,j,k)/8.64e07
       end do

       endif



165   continue





endif




end do
end do


!       call toc ('mulsub', '2')

end subroutine mulsub_vect

!####################################################################


subroutine cloudm_vect(i, j, TB_v,PB_v,TCC_v,PT_v,WV_v,RR_v,rcl_v,    &
!      TE_v,QE_v,PRESS_v,T_v,Q_v,SIG_v,   rh2o, cappa,  &
       TE_v,QE_v,PRESS_v,T_v,Q_v,SIG_v,         cappa,  &
!rair_mul,EPSILO,cpair_mul,GRAVIT,LATICE,KOU,TFRE,PRECIP_v,CONINT_v,  &
!rair_mul,       cpair_mul,GRAVIT,LATICE,KOU,TFRE,PRECIP_v,CONINT_v,  &
 rair_mul,       cpair_mul,GRAVIT,LATICE,KOU,dfre,TFRE,PRECIP_v,CONINT_v,  &
!rair_mul,       cpair_mul,       LATICE,KOU,TFRE,PRECIP_v,CONINT_v,  &
	DPF_v,  &
       DFR_v,DINT_v,flux_v,QLWA_v, debug_ijt)

logical, intent(in)                :: debug_ijt
integer, intent(in)                :: i, j, kou
!real,    intent(in)           :: rair_mul, epsilo, cpair_mul, gravit, &
 real,    intent(in)           :: rair_mul,         cpair_mul, gravit, &
!real,    intent(in)           :: rair_mul,         cpair_mul,         &
!			      latice, tfre, rh2o
				      latice, dfre, tfre
real, dimension(:,:), intent(in)   ::  tb_v, pb_v, rr_v, press_v
real, dimension(:,:), intent(out)  ::  pt_v, precip_v, conint_v,dint_v
real, dimension(:,:,:), intent(in)  :: t_v, q_v, sig_v
real, dimension(:,:,:), intent(out) :: tcc_v, wv_v, rcl_v, te_v, qe_v, &
                                       dpf_v, dfr_v, flux_v, qlwa_v

!
!     ONE-DIMENSIONAL CLOUD MODEL 
!     L. DONNER     NCAR     3 OCT 1984
!     Modified to eliminate dependence of cloud properties on
!     cloud fractional area.
!
!
!     CLOUD MODEL
!
!     ON INPUT:
!
!     TB     cloud temperature at cloud base (K)
!     PB     pressure at cloud base (Pa)
!     rr     cloud radius at base (m)
!     press  surface pressure (Pa)
!     t      large-scale temperature at GCM resolution (K)
!            Index 1 at physical bottom.
!     q      large-scale mixing ratio at GCM resolution (K)
!            Index 1 at physical bottom.
!     sig    (pressure/surface pressure) GCM vertical coordinate
!            Index 1 at physical bottom.
!     rair   gas constant for dry air (J/(kg K))
!     epsilo ratio of molecular weights, water vapor to dry air
!     cpair  specific heat for dry air (J/(kg K))
!     gravit gravity constant (m/(s**2))
!     LATICE LATENT HEAT OF FUSION (J/KG)
!     KOU    INDICATOR FOR ENTRAINMENT COEFFICEINT IN ENSEMBLE
!     dfre   freezing range (K)
!     TFRE   freezing temperature (K)
!
!     ON OUTPUT:
!     
!     DPF    CONDENSATION RATE*(CLOUD RADIUS**2)  ((M**2)*KG(H2O)/KG/S)
!     precip precipiation integral*(cloud radius**2) (kg/s)
!     CONINT CONDENSATION INTEGRAL*(CLOUD RADIUS**2)  (KG/S)
!     QLWA   CLOUD LIQUID WATER CONTENT (KG(H2O)/KG)
!     DFR    MOISTURE TENDENCY DUE TO FREEZING IN
!            CONVECTIVE UPDRAFT*(CLOUD RADIUS**2) ((M**2)*G(H2O)/(KG DAY))
!     DINT   INTEGRATED WATER MASS FROZEN IN CONVECTIVE
!            UPDRAFT*(CLOUD RADIUS**2) (KG(WATER) /SEC)
!     flux   updraft density*(radius**2)*vert vel (kg/s)
!            Index 1 at physical base of cloud.
!     tcc    cloud temperature (K)
!            Index 1 at physical base of cloud.
!     pt     pressure at cloud top (Pa)
!     wv     cloud vertical velocity (m/s)
!            Index 1 at physical base of cloud.
!     rcl    cloud radius (m)
!            Index 1 at physical base of cloud.
!     te     large-scale temperature at cloud-model resolution (K)
!            Index 1 at physical base of cloud.
!     qe     large-scale mixing ratio at cloud-model resolution (kg/kg)
!            Index 1 at physical base of cloud.
!
!     PARAMETER(NCM=100,NCM1=NCM-1)
!     PARAMETER(NLEV=40      )

      real, dimension (ncap) :: sub1, sub2, sub3, sub4, sub5, sub6, &
				sub7, pf, dpf, qllw, qe, te, dfr, &
				tcc, wv, rsc, qlwa, disc, disd, disp, &
				dis, disb, dadp, rcl, flux
      real, dimension ( size(t_v,3) ) :: sig, t, q

      logical :: testlc, test2
!     real    :: latvap


      save sub1sum, pcsave




  
       
        tb = tb_v(i,j)	
	pb  = pb_v(i,j)
	rr  = rr_v(i,j)
	press = press_v(i,j)
	  t(:) = t_v(i,j,:)
	  q(:) = q_v(i,j,:)
	  sig(:) = sig_v(i,j,:)
	  
      TESTLC=.TRUE.
      test2=.false.

      if (debug_ijt) then
          print *,   'DONNER_DEEP/cloudm: kou= ',kou
      endif

!
!     Set pstop (in Pa) as lowest presure to which cloud can extend.
!
      pstop=4.e03
!  LJD fix 3/21/01, see email

!     IF (PB .EQ. 0.) PT=0.
!     IF (PB .EQ. 0.) go to 307
      IF (PB .le. pstop) PT=pstop
      IF (PB .le. pstop) go to 307

      ALP=.5
      DINT=0. 
      DTFR=0.
!     LATVAP=2.5104E06
      FP=0.
      PF(1)=0.
      DP=-1000.
      if (kou .eq. 1) sub1sum=0.
      if (kou .eq. 1) pcsave=press
!
!     GATE:
!     These values are used for the GATE case in Donner (1993, JAS).
!     The most penetrative sub-ensemble must have KOU=KPAR.
!
      ALPP=.0915
      IF (KOU .EQ. 2) ALPP=ALPP/(1.3)
      IF (KOU .EQ. 3) ALPP=ALPP/(1.8)
      IF (KOU .EQ. 4) ALPP=ALPP/(2.5)
      IF (KOU .EQ. 5) ALPP=ALPP/(3.3)
      IF (KOU .EQ. 6) ALPP=ALPP/(4.5)
      IF (KOU .EQ. 7) ALPP=ALPP/(40.)
!
!     KEP:
!
!      ALPP=.0915
!      IF (KOU .EQ. 2) ALPP=ALPP/(1.22)
!      IF (KOU .EQ. 3) ALPP=ALPP/(1.56)
!      IF (KOU .EQ. 4) ALPP=ALPP/(2.05)
!      IF (KOU .EQ. 5) ALPP=ALPP/(2.6)
!      IF (KOU .EQ. 6) ALPP=ALPP/(3.21)
!      IF (KOU .EQ. 7) ALPP=ALPP/(7.84)
!      if (kou .eq. 1) alpp=alpp/(5.5)
!
      DO 3 k=1,ncap
      flux(k)=0.
      DPF(k)=0.
      DFR(k)=0.
      DISC(k)=0.
      DISD(k)=0.
      DIS(k)=0.
      QLWA(k)=0.
      SUB1(k)=0.
      SUB2(k)=0.
      SUB3(k)=0.
      SUB4(k)=0.
      SUB5(k)=0.
      SUB6(k)=0.
      SUB7(k)=0.
      PF(k)=0.
      rcl(k)=0.
      DPF(k)=0.
      DISB(k)=PB+(k-1)*DP
      TCC(k)=0.
      QLLW(k)=0.
      QE(k)=0.
      TE(k)=0.
 3    WV(k)=0.
      WV(1)=0.5
      rcl(1)=rr
      QCW=0.
      QLW=0.
      QRW=0.
      TCC(1)=TB
!     CALL ESTABL(ES,TB)
      CALL lookup_es(TB, ES)
      RSC(1)=ES*EPSILO/(PB-ES)
      PRECIP=0.
      CONINT=0.
      ILVM=NLEV-1
      N=0

      if (debug_ijt) then
       print *,   'DONNER_DEEP/cloudm: RR,PB,TB= ',RR,PB,TB
      endif


      ACCOND=0.
      ACPRE=0.
      dtupa=0.
      dfrac=0.
      sumfrea=0.

      if (debug_ijt) then
      do k=1,nlev-ktest_model+1
       print *,   'DONNER_DEEP/cloudm: k,sig   = ',k,sig(k)    
      end do
      endif

      DO 1 k=1,ncap-1
      IF (k .NE. 1) GO TO 10
!     CALL ESTABL(ES,TCC(k))
      CALL lookup_es(TCC(k), ES)
      RSC(k)=EPSILO*ES/(PB-ES)
      IH=0

      if (debug_ijt) then
       print *,  'DONNER_DEEP/cloudm: ILVM=',ILVM
      endif


      DO 2 JK=1,ILVM


      IF ( (PRESS*SIG(JK) .GE. PB) .AND. (PRESS*SIG(JK+1) .LE. PB) )  &
        IH=JK
 2    CONTINUE

      IF (IH .EQ. 0) THEN
         QE(1)=Q(1)
         TE(1)=T(1)
         GO TO 10
      END IF
      QE(1)=Q(IH)+ ( (Q(IH+1)-Q(IH) )*ALOG(PB/(SIG(IH)*PRESS))/  &
        ALOG(SIG(IH+1)/SIG(IH)) )
      TE(1)=T(IH)+ ( (T(IH+1)-T(IH))*ALOG(PB/(SIG(IH)*PRESS))/   &
        ALOG(SIG(IH+1)/SIG(IH)))

       if (debug_ijt) then
         print *,   'DONNER_DEEP/cloudm: QE,TE= ',QE(1),TE(1)
       endif
	

 10   IH=0

       if (debug_ijt) then
         print *,   'DONNER_DEEP/cloudm: TE,TCC= ',TE(1),TCC(1)
       endif


      P=PB+k*DP
!
      if ( p .lt. pstop) then

       if (debug_ijt) then
         print *,   'DONNER_DEEP/cloudm: pstop in Clouda= ',pstop
       endif


!      if (p .lt. pstop) test2=.true.
      if (p .lt. pstop) go to 4
      end if
!
      IF ( (P .LE. 50.E03) .AND. (TESTLC) ) THEN
         N=0
         GO TO 4
      END IF
      DO 11 JK=1,ILVM
      IF ( (PRESS*SIG(JK) .GE. P) .AND. (PRESS*SIG(JK+1) .LE. P) )  &
       IH=JK
  11  CONTINUE
      IF (IH .EQ. 0) THEN
         QE(k+1)=Q(1)
         TE(k+1)=T(1)
         GO TO 12
      END IF
      QE(k+1)=Q(IH)+ ( (Q(IH+1)-Q(IH))*ALOG(P/(PRESS*SIG(IH)))/  &
       ALOG(SIG(IH+1)/SIG(IH)))
      TE(k+1)=T(IH)+( (T(IH+1)-T(IH))*ALOG(P/(PRESS*SIG(IH)))/   &
       ALOG(SIG(IH+1)/SIG(IH)))
 12   CONTINUE
!     CALL ESTABL(ES,TCC(k))
      CALL lookup_es(TCC(k), ES)
      PP=P-DP
      DISP(k)=PP

      if (debug_ijt) then
        print *, 'DONNER_DEEP/cloudm: WV,PP,TE= ',WV(k),PP,TE(k)
        print *, 'DONNER_DEEP/cloudm: TCC(k),qe(k),QLW= ',TCC(k),WV(k),QLW
      endif




      CALL SIMULT(TCC(k),rcl(k),WV(k),PP,TE(k),QE(k),QLW,    &
!    TESTLC,pcsave,DTDP,DPD ,DWDP,EPSILO,rair_mul,LATVAP,LATICE,TFRE  &
!    TESTLC,pcsave,DTDP,DPD ,DWDP,EPSILO,rair_mul,       LATICE,TFRE  &
     TESTLC,pcsave,DTDP,DPD ,DWDP,       rair_mul,       LATICE,TFRE  &
!     ,rh2o , cappa, GRAVIT,cpair_mul,RMU,ALPP, debug_ijt)
            , cappa, GRAVIT,cpair_mul,RMU,ALPP, debug_ijt)
!           , cappa,        cpair_mul,RMU,ALPP, debug_ijt)

      if (debug_ijt) then
       print *,   'DONNER_DEEP/cloudm: QE,QLW,RR= ',QE(k),QLW,RR
        print *,  'DONNER_DEEP/cloudm: DPD,DWDP= ',DPD,DWDP
      endif


      TCEST=TCC(k)+DTDP*DP
      WEST=WV(k)+DWDP*DP
      rEST=rcl(k)+DPD*DP

       if (debug_ijt) then
        print *,  'DONNER_DEEP/cloudm: rest,west= ',rest,west
       endif


          IF (REST .LE. 0.) GO TO 4
	  IF (WEST .LT. .01) GO TO 4
!     CALL ESTABL(ESTEST,TCEST)
      CALL lookup_es(TCEST, ESTEST)
	  QCEST=EPSILO*ESTEST/(P-ESTEST)
      TEST=TE(k+1)
      QEST=QE(k+1)
      QRWT=QRW
      QCWT=QCW
      QLWT=QLW
      CALL MICRO(TCC(k),TCEST,PP,P,TE(k),TEST,QE(k),QEST,WV(k),   &
       WEST,rest,RMU,QRWT,QCWT,QLWT,DCW1,DQRW3, debug_ijt)

       if (debug_ijt) then
        print *, 'DONNER_DEEP/cloudm: P,TEST,QEST= ',P,TEST,QEST
       print *,  'DONNER_DEEP/cloudm: TCEST,rest,QLWT= ',TCEST,REST,QLWT
       endif


      CALL SIMULT(TCEST,rEST,WEST,P,TEST,QEST,QLWT,TESTLC,    &
!      pcsave,DTDP2,DPD2,DWDP2,EPSILO,rair_mul,LATVAP,LATICE,TFRE  &
!      pcsave,DTDP2,DPD2,DWDP2,EPSILO,rair_mul,       LATICE,TFRE  &
       pcsave,DTDP2,DPD2,DWDP2,       rair_mul,       LATICE,TFRE  &
!     ,rh2o , cappa, GRAVIT,cpair_mul,RMU,ALPP, debug_ijt)
            , cappa, GRAVIT,cpair_mul,RMU,ALPP, debug_ijt)
!           , cappa,        cpair_mul,RMU,ALPP, debug_ijt)

      if (debug_ijt) then
       print *,   'DONNER_DEEP/cloudm: TESTLC,DTDP2,DPD2= ',TESTLC,DTDP2,DPD2
       print *,   'DONNER_DEEP/cloudm: DWDP2= ',DWDP2
      endif


      rcl(k+1)=rcl(k)+DPD2*DP
      WV(k+1)=WV(k)+DWDP2*DP
      IF ((WV(k+1) .LT. .01) .OR. (rcl(k+1) .LE. 0.))  THEN
         rcl(k+1) = 0.0
         WV(k+1)=0.
         GO TO 4
      END IF
      DTDP=(DTDP+DTDP2)/2.
      DWDP=(DWDP2+DWDP)/2.
      DPD=(DPD+DPD2)/2.
      TCC(k+1)=TCC(k)+DTDP*DP
!
!     ADD EFFECT OF FREEZING
!
      IF ((TCC(k) .GE. TFRE) .AND. (TCC(k+1) .LE. TFRE) .AND.    &
        (DTFR .EQ. 0.)) THEN
         DTFR=QLWT*LATICE/cpair_mul
!
!     Multiply by factor to take account that not all this
!     water will freeze before falling out. Use Leary and Houze
!     (JAS,1980) rc/cu ratio to estimate this ratio.
      dtfr=.52*dtfr
!
      end if
         dfraca=(tfre-tcc(k+1))/dfre
         if (dfraca .gt. 1.) dfraca=1.
         dfrac=amax1(dfrac,dfraca)
         dtupb=dtfr*dfrac 
         dfr(k+1)=dtupb-dtupa
         sumfrea=sumfrea+dfr(k+1)
         dtupa=dtupb
         TCC(k+1)=TCC(k+1)+dfr(k+1)

	 if (debug_ijt) then
          print *,   'DONNER_DEEP/cloudm: DTFR=',DTFR,' dfr=',dfr(k+1),   &
			    'LATICE=',LATICE, 'CPAIR=',cpair_mul,'k= ',k
	 endif


      WV(k+1)=WV(k)+DWDP*DP
      rcl(k+1)=rcl(k)+DPD*DP
      IF ((WV(k+1) .LT. .01) .OR. (rcl(k+1) .LE. 0.))  THEN
         rcl(k+1) = 0.0
         WV(k+1)=0.
         GO TO 4
      end if
!
!     Ensure that all ensemble members have pt .le. pressure at which
!     ensemble member 1 becomes buoyant. Simult subroutine will use
!     value of pcsave to do so.
!
      if (kou .eq. 1) then
      if ((testlc) .and. (WV(k+1) .GT. WV(k))) pcsave=p
      IF (WV(k+1) .GT. WV(k)) TESTLC=.FALSE.
      else
      IF (WV(k+1) .GT. WV(k)) TESTLC=.FALSE.
      end if
!     CALL ESTABL(ES,TCC(k+1))
      CALL lookup_es(TCC(k+1), ES)
       RSC(k+1)=EPSILO*ES/(P-ES)

       if (debug_ijt) then
       print *,   'DONNER_DEEP/cloudm: TE,QE,RSC= ',TEST,QEST,RSC(k+1)
       endif

 
      RBAR=(rcl(k)+rcl(k+1))/2.
      RMUB=2.*ALPP/RBAR
      TVB=TE(k)*(1.+.61*QE(k))
      TVB=TVB+TEST*(1.+.61*QEST)
      TVB=TVB/2.
      TVC=TCC(k)*(1.+.61*RSC(k))
      TVC=TVC+TCC(k+1)*(1.+.61*RSC(k+1))
      TVC=TVC/2.
      DZDP=-(1.+ALP)*TVC*TVB*rair_mul
      DZDP=2.*DZDP/((P+PP)*GRAVIT*(ALP*TVB+TVC))
      sub1(k)=wv(k)/dzdp
      sub1(k+1)=wv(k+1)/dzdp
      DZ=DP*DZDP
!
!     CALL MICROPHYSICS ROUTINE.
!
      CALL MICRO(TCC(k),TCC(k+1),PP,P,TE(k),TE(k+1),QE(k),QE(k+1),  &
       WV(k),WV(k+1),RBAR,RMUB,QRW,QCW,QLW,DCW1,DQRW3, debug_ijt)
      DISC(k+1)=QCW
      DISD(k+1)=QRW
      TAV=TCC(k)*(1.+.61*RSC(k))
      TAVP1=TCC(k+1)*(1.+.61*RSC(k+1))
      QLWA(k+1)=QLW

     if (debug_ijt) then
       print *,   'DONNER_DEEP/cloudm: TCC(k+1),WV(k+1),QLW= ',TCC(k+1),  &
				WV(k+1),QLW
       print *,   'DONNER_DEEP/cloudm: DCW1,DQRW3= ',DCW1,DQRW3
     endif


      RM=rair_mul*(1.+.609*RSC(k))
      RMA=rair_mul*(1.+.609*RSC(k+1))
      CPM=cpair_mul*(1.+.87*RSC(k))
      CPMA=cpair_mul*(1.+.87*RSC(k+1))
      WVB=(WV(k)+WV(k+1))/2.
!
!
      pdab=rbar**2
      pdam=(rcl(k))**2
      pdap=(rcl(k+1))**2
!
      DIS(k)=DQRW3*PDAB*.001*WVB*GRAVIT/DP
      PF(k)=.001*DCW1*WVB*PDAB*GRAVIT/DP
!
!     CALCULATE MOISTURE SUBJECT TO FREEZING
!
      IF (  (DFR(k+1) .NE. 0.) ) THEN
      DFR(k+1)=DFR(k+1)*P*WV(k+1)*pdap*GRAVIT/(rair_mul*TAVP1*DP)
      DFR(k+1)=-DFR(k+1)*cpair_mul/LATICE
      DINT=DINT-DFR(k+1)*DP/GRAVIT
      DFR(k+1)=DFR(k+1)*8.64E07
      END IF
!
      ACCOND=ACCOND+DCW1
      ACPRE=ACPRE+DQRW3
      ACTOT=ACCOND+ACPRE

      if (debug_ijt) then
       print *,   'DONNER_DEEP/cloudm: k,ACCOND,ACPRE= ',k,ACCOND,ACPRE
        print *,  'DONNER_DEEP/cloudm:  k,ACTOT= ',k,ACTOT
      endif


      IF (k .EQ. 1) FLUX(k)=(rcl(k)**2)*WV(k)*PP/(RM*TCC(k))
      FLUX(k+1)=(rcl(k+1)**2)*WV(k+1)*P /(RMA*TCC(k+1))

      if (debug_ijt) then
        print *,  'DONNER_DEEP/cloudm: k,rcl,RMU= ',k,rcl(k),RMU
      endif


      IF (QLW .lt. 0.) GO TO 4
      PI=(1.0E05/PP)**(rair_mul/cpair_mul)
      PIP=(1.0E05/P)**(rair_mul/cpair_mul)
      PIM=PP-DP
      PIM=(1.0E05/PIM)**(rair_mul/cpair_mul)
      IF (k .EQ. 1) THEN
         SUB1(k)=pdam*sub1(k)*(RSC(k+1)-RSC(k))/DP
!         SUB1(k)=-SUB1(k)*PI*latvap/cpair_mul
         SUB3(k)=-(PIP*TCC(k+1)-PI*TCC(k))/DP
      END IF
      IF (k .NE. 1) THEN
         SUB1(k)=pdam*sub1(k)*(RSC(k+1)-RSC(k-1))/(2.*DP)
!         SUB1(k)=-PI*latvap*SUB1(k)/cpair_mul
         SUB3(k)=-(PIP*TCC(k+1)-PIM*TCC(k-1))/(2.*DP)
         END IF
	 sub1(k)=sub1(k)-pf(k)
         RMUP=RMUB*DZDP
         SUB2(k)=RMUP*PI*latvap*(QE(k)-RSC(k))/cpair_mul
         SUB4(k)=RMUP*PI*(TE(k)-TCC(k))
      N=k
 1    CONTINUE
 4    CONTINUE

      if (debug_ijt) then
       print *,   'DONNER_DEEP/cloudm: DINT IN CLOUDM,sumfrea= ',DINT,sumfrea
      endif


      KK=N
      IF (KK .EQ. 0.) GO TO 23

	 if (debug_ijt) then
           print *,  'DONNER_DEEP/cloudm: n,rsc(n+2),rsc(n)= ',n,   &
						 rsc(n+1),rsc(n)
	 endif

 
	 pdanp=rcl(n+1)**2
         SUB1(n+1)=pdanp*sub1(n+1)*(RSC(n+1)-RSC(n))/(2.*DP)
 24   FORMAT(2X,'PRECIP=',E15.5,'KG/(M**2)*SEC')
!      sub1sum=sub1sum+sub1(n+1)*dp
      DO 22 K=1,KK
      P=PB+(K-1)*DP
       sub1sum=sub1sum+sub1(k)*dp
       IF (K .EQ. 1) THEN
          DPF(1)=PF(1)
       ELSE
          DPF(K)=(PF(K)+PF(K-1))/2.
       END IF
       IF (K .EQ. KK) DPF(KK+1)=PF(KK)/2.

	if (debug_ijt) then
        print *,   'DONNER_DEEP/cloudm: K,PF,DIS= ',K,PF(K),DIS(K)
	endif


       PI=(1.0E05/P)**(rair_mul/cpair_mul)
       TV=TCC(K)*(1.+.61*RSC(K))
       SUB5(K)=-PI*latvap*DPF(K)*rair_mul*TV/(cpair_mul*   &
       WV(K)*P*GRAVIT*pdam)
!       SUB6(K)=SUB1(K)+SUB2(K)+SUB5(K)
       SUB6(K)=SUB1(K)+SUB2(K)
!       SUB7(K)=SUB3(K)+SUB4(K)-SUB5(K)
       SUB7(K)=SUB3(K)+SUB4(K)
       SUB6(K)=-SUB6(K)*WV(K)*P*GRAVIT*PDAM/(rair_mul*TV)
       SUB7(K)=-SUB7(K)*WV(K)*P*GRAVIT*PDAM/(rair_mul*TV)

       if (debug_ijt) then
        print *,   'DONNER_DEEP/cloudm: K,P,SUB1= ',K,P,SUB1(K)
       endif


       if (k .eq. n) then

	 if (debug_ijt) then
           print *,  'DONNER_DEEP/cloudm: N,P,SUB1+= ',n,P,SUB1(n+1)
	 endif


       end if
       DIS(K)=-DIS(K)
       CONINT=CONINT+PF(K)*DP/GRAVIT
 22   CONTINUE
      if (kou .eq. 7) sub1sum=86400.*sub1sum/gravit

      if (debug_ijt) then
        print *,  'DONNER_DEEP/cloudm: kou= ',kou,' Vapr Adv= ',sub1sum,  &
					  ' mm/day'
        print *,  'DONNER_DEEP/cloudm: CONINT= ',CONINT,' KG/(M**2)/SEC'
      endif


      PRECIP=CONINT*ACPRE/ACCOND
      SUB6(1)=0.
      SUB7(1)=0.
 23   PT=PB+KK*DP

       if (debug_ijt) then
         print *, 'DONNER_DEEP/cloudm: kou,pcsave= ',kou,pcsave
       endif


       if (test2) then

       if (debug_ijt) then
	 do k=1,ncap
	   print *, 'DONNER_DEEP/cloudm: p & r    ',  &
		disb(k), rcl(k)
           IF (WV(k+1) .LE. 0.) exit      
	 end do
	 do k=1,ncap
	   print *, 'DONNER_DEEP/cloudm: p & te   ',  &
		disb(k),  te(k)
           IF (WV(k+1) .LE. 0.) exit      
	 end do
	 do k=1,ncap
	   print *, 'DONNER_DEEP/cloudm: p & wv   ',  &
		disb(k),  wv(k)
           IF (WV(k+1) .LE. 0.) exit      
	 end do
	 do k=1,ncap
	   print *, 'DONNER_DEEP/cloudm: p & qlw  ',  &
		disb(k), qlwa(k)
           IF (WV(k+1) .LE. 0.) exit      
	 end do
	 do k=1,ncap
	   print *, 'DONNER_DEEP/cloudm: p & tcc  ',  &
		disb(k),  tcc(k)
           IF (WV(k+1) .LE. 0.) exit      
	 end do
       endif

       end if


       if (test2) then
!      CALL AGSETF('Y/ORDER.',1.)
!      CALL ANOTAT('M','KPA',0,0,-1,0)
!      CALL EZXY(rcl,DISB,95,'RADIUS')
!      CALL ANOTAT('K/PA','KPA',0,0,-1,0)
!      CALL EZXY(SUB1,DISB,95,'SUB1')
!      CALL EZXY(SUB2,DISB,95,'SUB2')
!      CALL EZXY(SUB3,DISB,95,'SUB3')
!      CALL EZXY(SUB4,DISB,95,'SUB4')
!      CALL EZXY(SUB5,DISB,95,'SUB5')
!      CALL EZXY(SUB6,DISB,95,'SUB6')
!      CALL EZXY(SUB7,DISB,95,'SUB7')
!      CALL ANOTAT('KG(H2O)/KG/SEC','KPA',0,0,-1,0)
!      CALL EZXY(DIS,DISB,95,'PRECIP')
!      CALL ANOTAT('KG(H2O)/KG','PA',0,0,-1,0)
!      CALL EZXY(QLWA,DISB,95,'LIQ WATER')
!      CALL ANOTAT('G(H2O)/M**3','PA',0,0,-1,0)
!      CALL EZXY(DISC,DISB,95,'CLOUD WATER')
!      CALL ANOTAT('G(H2O)/M**3','PA',0,0,-1,0)
!      CALL EZXY(DISD,DISB,95,'RAIN WATER')
!      IF (TEST2 .AND. (PT .LE. 60000.) ) THEN
!      CALL AGSETF('Y/ORDER.',1.)
!      CALL ANOTAT('KG/((M**2)*S)','KPA',0,0,-1,0)
!      CALL EZXY(FLUX,DISB,95,'FLUX')
!      CALL ANOTAT('  /PA        ','KPA',0,0,-1,0)
!      CALL EZXY(DADP,DISB,95,'DADP')
!      CALL ANOTAT('             ','KPA',0,0,-1,0)
!      CALL EZXY(RCL,DISB,95,'RCL')
!     CALL ANOTAT('M/S          ','KPA',0,0,-1,0)
!     CALL EZXY(WV,DISB,95,'WV')
       END IF
!     IF (RR .EQ. 2000.)    &
!      CALL ANOTAT(14HPRESSURE(KPA)$,23HVERTICAL VELOCITY(M/S)$,0, &
!      0,-1,0)
!     IF (RR .EQ. 2000.)    &
!      CALL EZMXY(DISB,DISA,95,3,95,18HVERTICAL VELOCITY$)
      if (test2) stop
!      RETURN



        dpf_v(i,j,:) = dpf(:)
	precip_v(i,j) = precip
	conint_v(i,j) = conint
	dfr_v(i,j,:) = dfr(:)
	dint_v(i,j) = dint
	tcc_v(i,j,:) = tcc(:)
	wv_v(i,j,:) = wv(:)
	rcl_v(i,j,:) = rcl(:)
	te_v(i,j,:) = te(:)
	qe_v(i,j,:) = qe(:)
	flux_v(i,j,:) = flux(:)
	qlwa_v(i,j,:) = qlwa(:)

    307  continue

	pt_v(i,j) = pt

!      end do
!      end do

end subroutine cloudm_vect

!#####################################################################

subroutine lcl(tb,constab,pb,qb,press,sig,t,qin)

real, intent(inout)   ::  pb
real, intent(out)   :: tb, qb
logical, intent(inout) :: constab
real, dimension(:), intent(in) :: t, sig, qin
real, intent(in)     :: press


      istart=2
      dp=-100.
      q=qin(istart)
      tc=t(istart)
      do i=1,600
      p=press*sig(istart)+dp*(i-1.)
      pt=press*sig(istart)+dp*i
      GAM=.286*(1.-.26*Q)*TC/(P/      press       )
      TC=TC+GAM*DP/    press
!     CALL ESTABL(ES,TC)
      CALL lookup_es(TC, ES)
      RS=.622*ES/(PT-ES)
      IF (RS .LE. Q) PB=PT
      IF (RS .LE. Q) TB=TC
      IF (RS .LE. Q) QB=qin(istart)
      IF (RS .LE. Q) GO TO 6
      end do
      IF (PB .EQ. 0.) CONSTAB=.TRUE.
 6    CONTINUE



end subroutine lcl


!#####################################################################


subroutine satad(t,p,lat,dtp)
!!  IT APPEARS THAT THIS ROUTINE IS NOT CURRENTLY ACCESSED BY THE MODEL

real, intent(in) :: t, p,lat
real, intent(out) :: dtp


!
!      computes saturated adiabatic lapse rate
!
!      on input:
!        t    temperature (K)
!        p    pressure (Pa)
!        lat  latent heat (J/kg)
!      on output:
!        dtp  saturated adiabatic lapse rate (K/m)
!
       real lat
!
!     constants
!
      rair_satad=287.04
      gravit=9.80616
       cpair_satad=1004.64
      epsilo_satad=.622
      rh2o_satad=461.
!     write(6,*) 'satad t,p= ',t,p
!
!     calculate saturation specific humidity
!
!     call establ(es,t    )
      call lookup_es(t, es    )
      qs=epsilo_satad*es/(p-es)
      tv=t*(1.+.61*qs)
!
!     calculate saturated adiabatic lapse rate
!
      dtp=tv+(lat*qs/rair_satad)
      dtp=-dtp*gravit/(cpair_satad*t)
      desdt=lat*es/(rh2o_satad*(t**2))
      dtpd=1.+(epsilo_satad*lat*desdt/(cpair_satad*p))
      dtp=dtp/dtpd
!     return

end subroutine satad





!####################################################################


subroutine simult(tc,rda,wv,p,te,qe,qlw,testlc,pcsave,   &
!      dtdp,drdp,dwdp,epsilo,rair_mul,latvap,latice,tfre   &
!      dtdp,drdp,dwdp,epsilo,rair_mul,       latice,tfre   &
       dtdp,drdp,dwdp,       rair_mul,       latice,tfre   &
!      ,rh2o, cappa   &
            , cappa   &
       ,gravit,cpair_mul,rmu,alpp, debug_ijt)
!             ,cpair_mul,rmu,alpp, debug_ijt)

!--------------------------------------------------------------------
logical, intent(in) :: testlc, debug_ijt
real, intent(in) ::  tc, rda, wv, p, te, qe, qlw, pcsave, &
!	     epsilo, rair_mul, latvap, latice, tfre, rh2o, gravit, &
!	     epsilo, rair_mul,cappa,   latice, tfre, gravit, &
 	             rair_mul,cappa,   latice, tfre, gravit, &
!		             rair_mul,cappa,   latice, tfre,         &
		     cpair_mul, alpp
real, intent(out) :: dtdp, drdp, dwdp, rmu
!--------------------------------------------------------------------
!
!     Generates cloud profiles.
!     See LJD "Cloud Model 89" notes on Generalized mu (10/1/89)
!     and dwdz (10/2/89). The value of epsilon is taken as 1.
!     Version where cloud propeties independent of cloud area.
!
!     On input:
!
!        tc   cu temperature (K) at pressure p (Pa)
!        rda  cu radius at pressure p
!        te   environmental temperature (K) at pressure p
!        qe   environmental mixing ratio at pressure p
!        wv   cu vertical velocity (m/s) at pressure p
!        qlw  liquid water
!        testlc  indicator(logical)
!        pcsave  pressure at which cloud ensemble 1 becomes
!                buoyant (Pa)
!
!     On Output:
!
!        dtdp    temperature derivative (K/Pa)
!        drdp    cu radius derivative (m/Pa)
!        dwdp    cu vertical velocity derivative (m/s/Pa)
!        rmu  entrainment coefficient (/s)
!
!
!     logical debug_ijt
!     logical testlc
!      real latvap,lat,latice



      real lat



!
!      assign constants
!
!     rh2o=461.
      alp=.5
      dp=-1000.
      epm=0.
      lat=latvap
      if (tc .lt. tfre) lat=latvap+latice
!
!     call establ(es,tc)
      call lookup_es(tc, es)
      rsc=epsilo*es/(p+(epsilo-1.)*es)
       htve=te   *(1.+.609*qe   )
      HTV=TC    *(1.+.609*RSC   )

      if (debug_ijt) then
        print *,  'DONNER_DEEP/simult: htve,htv= ',htve,htv
      endif

      DA=rair_mul*(1.+ALP)*htv*htve    /(GRAVIT*p * (ALP*htve +htv   ))
 102  FORMAT(2X,3E15.5,I5,2E15.5)
      DER=latvap*ES/(RH2O*(TC **2   ))
      rmu=2.*alpp/rda
      fm=p*(rda**2)*wv/(rair_mul*htv)
      rmup=-da*rmu
      fmp=fm*exp(rmup*dp)

      if (debug_ijt) then
        print *,  'DONNER_DEEP/simult: p,fm,fmp= ',p,fm,fmp
      endif


      rst=rair_mul*(1.+.608*rsc)
      c2=(htv+(lat*rsc/rst))*gravit/(cpair_mul*tc)

      if (debug_ijt) then
        print *,  'DONNER_DEEP/simult: c2= ',c2
        print *,  'DONNER_DEEP/simult: te,p,qe= ',te,p,qe
      endif


!     call tae(te,p,qe,lat,epsilo,dp,rair_mul,cpair_mul, cappa,teae)
      call tae(te,p,qe,lat,       dp,rair_mul,cpair_mul, cappa,teae)

      if (debug_ijt) then
        print *,  'DONNER_DEEP/simult: teae= ',teae
      endif


      tcae=tc*exp(lat*rsc/(cpair_mul*tc))
      dz=-dp*da
      c3=(tcae-teae)*(1.-exp(rmu*dz))/(exp(rmu*dz)*dz   &
       *exp(lat*rsc/(cpair_mul*tc)))
      c4=1.+(epsilo*lat*der/(cpair_mul*p))

      if (debug_ijt) then
        print *,  'DONNER_DEEP/simult: c3,c4= ',c3,c4
      endif


      dtdp=(c2-c3)*da/c4
      test=tc+dtdp*dp
      c5=gravit*(htv-htve)/(htve*(1.+alp))
!
!     No  Liq  Remove Next Line
!
      c5=c5-gravit*qlw
      dwdz=2.*c5
      wtest=(wv**2)+dwdz*dz

      if (debug_ijt) then
        print *,  'DONNER_DEEP/simult: testlc,wtest= ',testlc,wtest
      endif


      dwdp=(wtest-(wv**2))/dp
      if ( ((.not. testlc) .and. (p .le. pcsave)) .or.   &
          (dwdp .le. 0.) ) then
      if  (wtest .le. 0.)  return
      wtest=sqrt(wtest)
      dwdp=(wtest-wv)/dp
      dwdz=wv*(1.-exp(rmu*dz))
      dwdz=dwdz/(dz*exp(rmu*dz))
      dwdp=(-dwdz*da)+dwdp
      end if
      
      if (debug_ijt) then
        print *,  'DONNER_DEEP/simult: testlc,dwdp,test= ',testlc,dwdp,test
      endif


      if ( (testlc) .and. (dwdp .gt. 0.) ) dwdp=0.
      if ( (.not. testlc) .and. (p .gt. pcsave) .and.    &
        (dwdp .gt. 0.) ) dwdp=0.
      west=wv+dwdp*dp
!     call establ(es,test)
      call lookup_es(test, es)
!      rsc=epsilo*es/(p-es)
      rsc=epsilo*es/(p+(epsilo-1.)*es)
      htv=test*(1.+.609*rsc)
      pdap=fmp*rair_mul*htv/((p+dp)*west)
      dadp=(pdap-(rda**2))/dp
      drdp=dadp/(2.*rda)
!     return

end subroutine simult




subroutine meens (cu,rc,pr,phr,ps,pb,pt,gravit     &
!       ,cappa,epsilo,rair_mul,cpair_mul,lat,rh2o,t,q,apt  &
       ,cappa,epsilo,rair_mul,cpair_mul,lat,     t,q,apt  &
!      ,cappa,       rair_mul,cpair_mul,lat,     t,q,apt  &
       ,rlhr,emfhr,ca,tmel,cuml              &
       ,cmu,cmui,dmeml,emd,emdi,eme,emei,wmm,wmp,elt,contot, &
       tmes,qmes,umeml, debug_ijt)


!------------------------------------------------------------------
logical, intent(in)          ::  debug_ijt        
real, dimension(:), intent(in) :: pr, phr, t, q, rlhr, emfhr
real, dimension(:), intent(out) :: cuml, cmu, dmeml, emd, eme, wmm, &
!			   wmp, elt, umeml, tmes, qmes
				   wmp,      umeml, tmes, qmes
real,         intent(in)  :: cu, rc, ps, pb, pt, gravit, cappa, &
!      epsilo, rair_mul, cpair_mul, lat, rh2o, apt, &
!!    epsilo, rair_mul, cpair_mul, lat,       apt, &
              rair_mul, cpair_mul, lat,       apt, &
			         tmel
real, intent(inout) ::  ca
real, dimension(:), intent(inout) ::  elt
real,   intent(out) :: cmui, emei, contot, emdi
!------------------------------------------------------------------




!     Calculate mesoscale heat and moisture sources, using
!     variation on Leary and Houze (JAS, 1980).
!     Performs calculations after all sub-ensemble calculations complete.
!     For notation, see "Cu Closure A notes," 2/97
!
!     On Input:
!       cu      condenstation integral
!               sigma(i=1,N) (a(i,p_b)/a(1,p_b))*
!               (sigma(j=1,M)((r(i,j)**2)/(r(i,p_b)**2))*c_u(i,j))
!               (mm/day)
!       rc      precipitation integral
!               sigma(i=1,N) (a(i,p_b)/a(1,p_b))*
!               (sigma(j=1,M)((r(i,j)**2)/(r(i,p_b)**2))*r_c(i,j))
!               (mm/day)
!       pr,phr  low-resolution pressure levels (Pa)
!       ps      surface pressure (Pa)
!       pb      cloud-base pressure (Pa)
!       pt      cloud-top pressure (Pa)
!       gravit  gravity constant (m/s**2)
!       cappa   ratio of gas constant to specific heat for dry air
!       epsilo  ratio of gas constants, dry air to water vapor
!       rair    gas constant for dry air (J/kg/K)
!       cpair   specific heat for dry air (J/kg/K)
!       lat     latent heat for phase change (J/kg)
!       rh2o    gas constant for water vapor (J/kg/K)
!       t       low-resolution temperature (K)
!       q       low-resolution specific humidity (kg(H2O)/kg)
!       apt     mesoscale fraction, normalized by [a(1,p_b)*5]
!       rlhr    latent heat release (kg/kg/s)-high resolution
!               sigma(i=1,N)(a(i,p_b)/a(1,p_b))*((r(i,j)**2)/(r(i,p_b)**2))*
!               (sigma(k=1,4) gamma^*(k,i,j))
!       emfhr   moisture flux convergence (kg/kg/s)-high resolution
!               sigma(i=1,N)(a(i,p_b)/a(1,p_b))*d(((r(i,j)**2)/(r(i,p_b)**2))*
!               omega*'(i,j) * q*'(i,j))/dp
!       ca      condensed water X-fer from cells to anvil (mm/day)
!               weighted as cu,rc
!       tmel    melting temperature (K)
!
!     On Output:
!       cuml    fractional mesoscale area, normalized by
!               a(1,p_b) at resolution of GCM
!       cmu     water mass condensed in mesoscale updraft
!               (g/kg/day) (normalized by a(1,p_b))
!       cmui    vertical integral of mesoscale-updraft deposition
!               (kg(H2O)/((m**2)*sec) 
!       dmeml   mass flux in mesoscale downdraft (kg/((m**2) s))
!               (normalized by a(1,p_b)) (index 1 at atmosphere bottom)
!               (resolution of GCM)
!       emd     water mass evaporated in mesoscale
!               downdraft (g/kg/day) (normalized by a(1,p_b))
!       emdi    vertical integral of mesoscale-downdraft sublimation
!               (mm/d)
!       eme     water mass evaporated from mesoscale
!               updraft (g/kg/day) (normalized by a(1,p_b))
!       emei    vertical integral of mesoscale-updraft sublimation
!               (kg(h2O)/((m**2)*sec)
!       wmm     water vapor removal by condensation of
!               cell vapor source (g/kg/day) (normalized by a(1,p_b))
!       wmp     water vapor redistributed from cell vapor source
!               (g/kg/day) (normalized by a(1,p_b))
!       elt     melting of ice in mesoscale updraft-
!               equivalent (g/kg/day)-which falls as meso sfc precip
!               (normalized by a(1,p_b))
!       contot  ratio of convective to total precipitation
!       tmes    temperature tendency due to mesoscale entropy-flux-
!               convergence (K/day) (normalized by a(1,p_b))
!       qmes    moisture tendency due to mesoscale moisture-flux
!               convergence (g/kg/day) (normalized by a(1,p_b))
!       umeml   mass flux in mesoscale updraft (kg/((m**2) s))
!               (normalized by a(1,p_b)) (index 1 at atmosphere bottom)
!               (resolution of GCM)
!

!


      real, dimension (size(t)+1) :: emt, emq
      real, dimension (size(t)) :: owm, omv, emdx, tempq, tempqa, &
				   tempt
      real, dimension (ncap) ::  wmhr, wphr, p, cumh, dmemh

      logical :: thetl







      if (debug_ijt) then
	print *, 'DONNER_DEEP/meens: entering meens'
      endif


!
!
!      define constants
!
      dp=-1000.
      ptt=pt+dp
      pztm=ptt-30.e03
      tpri=1.
      po=100.e03
!
!     Restrict pztm to .ge. 10 kPa, cf Ackerman et al (JAS,1988)
!     (unless pt .le. 10kPa)
!
      if (pztm .lt. 10.e03) pztm=10.e03
      if (ptt .lt. 10.e03) pztm=ptt+dp
      thetl=.false.
!
!      water-budget constants calculated from flux-
!      condensation parameterization
! 
       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: rc,cu= ',rc,cu
       endif


      gnu=rc/cu
!
!     maintain Leary and Houze ratios of alpha, beta, and
!     eta to their sum
!
      sum=1.-gnu
      alrat=.25
      berat=.13
      etrat=.62
      alpha=alrat*sum
      beta=berat*sum
      eta=etrat*sum
!
!     use Leary and Houze ratios for gnu-sub-m
!
      gnum=.5

       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: gnu= ',gnu
       endif

!
!     calculate mass of cu incorporated into mesoscale anvil,
!     integrated water evaporated in convective downdraft,
!     and integrated water evaporated from convective updraft
 
       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: alpha,beta,eta= ',alpha,beta,eta
       endif


      if (ca .ge. 0.) ca=eta*cu
      if (ca .lt. 0.) ca=-ca

       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: ca= ',ca
       endif

!
!     condensation of tower vapor source
!
      emt(nlev+1)=0.
      emq(nlev+1)=0.
      do 18 jk=1,nlev
         emt(jk)=0.
         emq(jk)=0.
         emd(jk)=0.
	 eme(jk)=0.
         qmes(jk)=0.
         tmes(jk)=0.
         omv(jk)=0.
         tempqa(jk)=q(jk)
         wmp(jk)=0.
         wmm(jk)=0.
         cmu(jk)=0.
         tempq(jk)=q(jk)
         tempt(jk)=t(jk)
	 cuml(jk)=0.
	 dmeml(jk)=0.
	 umeml(jk)=0.
 18   continue
      do 12 i=1,ncap
      wmhr(i)=0.
      wphr(i)=0.
      cumh(i)=0.
      dmemh(i)=0.
      p(i)=pb+(i-1)*dp
 12   continue
      pzm=0.
      do 11 i=1,ncap



      cmfhr=-rlhr(i)+emfhr(i)
      if (cmfhr .gt. 0.) then
         wmhr(i)=-cmfhr
         if (pzm .eq. 0.) then
            pzm=p(i)
         end if
      end if


 11   continue
      do i=1,ncap

       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: i,rlhr,emfhr= ',i,rlhr(i),emfhr(i)
       endif

       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: i,wmhr= ',i,wmhr(i)
       endif

      end do

      if (pzm .eq. 0.) pzm=pt
      sbl=0.
      call verav(wmhr,pb,pt,thetl,owm,sbl,ps,pr,phr,cappa, debug_ijt)

       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: pzm= ',pzm
      do jk=1,nlev
        print *, 'DONNER_DEEP/meens: jk,owm= ',jk,owm(jk)
      end do
       endif


      ome=-.463
        ampt=5.*apt
         tme=6.48e04
!
!     calculate redistribution of Cu H2O source
!
      do 13 jk=1,nlev
         if (owm(jk) .ge. 0.) go to 16
         pc1=phr(jk)
         pc2=phr(jk+1)
         omer=ome
         pctm=pc2+ome*tme
         if (pctm .le. pztm) omer=(pztm-pc2)/tme
         pctm=pc2+omer*tme
         q1=owm(jk)*(pc2-pc1)*tme/(pc1-pc2-omer*tme)
         q4=q1/2.

       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: pctm,pztm,q4= ',pctm,pztm,q4
       endif

         do 17 jj=jk,nlev



            if (phr(jj) .lt. pctm) go to 16
            tempq(jj)=tempq(jj)+(q1/ampt)
            tempqa(jj)=tempqa(jj)+(q4/ampt)
            wmpt=(q1/tme)
            if (phr(jj+1) .le. pctm) wmpt=wmpt*(phr(jj)-pctm)/   &
            (phr(jj)-phr(jj+1))
            wmp(jj)=wmpt+wmp(jj)


 17    continue

       do jj=jk,nlev
       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: jj,pr= ',jj,pr(jj)
       endif
            if (phr(jj) .lt. pctm) go to 216
       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: jj,q1,tempq,wmm= ',jj,q1,tempq(jj),wmm(jj)
        print *, 'DONNER_DEEP/meens: wmp= ',wmp(jj)
        print *, 'DONNER_DEEP/meens: jj,tempqa= ',jj,tempqa(jj)
       endif

       end do

       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: jk,q1,tempq,wmm= ',jk,q1,tempq(jk),wmm(jk)
!       print *, 'DONNER_DEEP/meens: jk,wmp,owm= ',jk,wmp(jk),owm(jk)
       endif

 216   continue
 16   continue


       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: jk,wmp,owm= ',jk,wmp(jk),owm(jk)
       endif


         wmp(jk)=wmp(jk)+owm(jk)
 13   continue
!
!     calculate condensed portion of redistributed H2O
!
	  do 23 jk=1,nlev
         if ( (phr(jk+1) .le. pzm) .and. (phr(jk) .ge. pztm) ) then
         pc2=phr(jk+1)
         omer=ome
         pctm=pc2+ome*tme
         if (pctm .le. pztm) omer=(pztm-pc2)/tme
         pctm=pc2+omer*tme
      ta=t(jk)+tpri
!  call establ(es,ta)
	  call lookup_es(ta, es)
	  qsat=epsilo*es/(pr(jk)+(epsilo-1.)*es)
	  q3=qsat-tempq(jk)
	  q5=qsat-tempqa(jk)
	  if (q3 .le. 0.) then
		 wmmt=q3*ampt/tme
		 if (phr(jk+1) .le. pctm) wmmt=wmmt*(phr(jk)-pctm)/  &
         (phr(jk)-phr(jk+1))
		 wmm(jk)=wmmt
      end if
	  if (q5 .le. 0.) tempqa(jk)=qsat
      end if
 23   continue
!
!     calculate meso up con portion due to meso lift
!
         anv=0.
         do 14 jk=1,nlev
            if (pr(jk) .gt. pzm) go to 19
               if (anv .eq. 0.) qref=tempqa(jk)
               anv=1.
            if (pr(jk) .lt. pztm) go to 21
            te=t(jk)+tpri
!           call establ(es,te)
            call lookup_es(te, es)
            qs=epsilo*es/(pr(jk)+(epsilo-1.)*es)
            jsave=jk

       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: qs,tempqa= ',qs,tempqa(jk)
       endif


            if (qref .ge. qs) go to 21
 19      continue
 14     continue
 21     continue
        alp=6.*ome/((pzm-pztm)**2)
        do 22 jk=1,nlev
	    if (jk .eq. 1) then
	       jkm=jk
            else
	       jkm=jk-1
            end if
	    if (jk .eq. nlev) then
	       jkp=jk
            else
	       jkp=jk+1
            end if
            if (pr(jk) .lt. pztm) go to 20
            if (pr(jk) .gt. pzm) go to 24
            pp=phr(jk+1)
            pm=phr(jk)
            if (phr(jk+1) .lt. pztm) pp=pztm
            if (phr(jk) .gt. pzm) pm=pzm
!
!        Calculate mesoscale entropy-flux convergence
!
         tmes(jk)=(pzm+pztm)*(rair_mul-cpair_mul)*(pp-pm)/cpair_mul
         tmes(jk)=tmes(jk)+((2.*cpair_mul-rair_mul)*((pp**2)-(pm**2))/   &
        (2.*cpair_mul))
         tmes(jk)=tmes(jk)-(rair_mul*pztm*pzm/cpair_mul)*alog(pp/pm)
         tmes(jk)=tmes(jk)/(phr(jk+1)-phr(jk))
!         tmes(jk)=tmes(jk)*ampt*tpri*alp/(1.-ampt)
         tmes(jk)=tmes(jk)*ampt*tpri*alp
!
!         Calculate mesoscale vertical velocities
!
            omv(jk)=(pzm+pztm)*((pp**2)-(pm**2))/2.
            omv(jk)=omv(jk)-(((pp**3)-(pm**3))/3.)
            omv(jk)=omv(jk)-pztm*pzm*(pp-pm)
            omv(jk)=omv(jk)/(phr(jk+1)-phr(jk))
            omv(jk)=omv(jk)*alp
          if (jk .lt. jsave) go to 24
            pre=pr(jk)
	       te=t(jk)+tpri
!       call establ(es,te)
	       call lookup_es(te, es)
               tempqa(jk)=epsilo*es/(pr(jk)+(epsilo-1.)*es)
	       if (qref .ge. tempqa(jk)) then
               tep=t(jkp)+tpri
            if (pr(jkp) .le. pztm) then
               cmu(jk)=-omv(jk)*(tempqa(jk)-tempqa(jkm))/(pr(jk)   &
              -pr(jkm))
            else if (jk .eq. jsave) then
!              call establ(es,tep)
               call lookup_es(tep, es)
               tempqa(jkp)=epsilo*es/(pr(jkp)+(epsilo-1.)*   &
              es)
              cmu(jk)=-omv(jk)*(tempqa(jkp)-tempqa(jk))/(pr(jkp)   &
               -pr(jk))
              qref=tempqa(jkp)
            else
!              call establ(es,tep)
               call lookup_es(tep, es)
               tempqa(jkp)=epsilo*es/(pr(jkp)+(epsilo-1.)*es)
               cmu(jk)=-omv(jk)*(tempqa(jkp)-tempqa(jkm))/(pr(jkp)   &
              -pr(jkm))
              qref=tempqa(jkp)
            end if
	    if (cmu(jk) .lt. 0.) cmu(jk)=0.
	    else
	    cmu(jk)=0.
	    end if
         cmu(jk)=cmu(jk)*ampt*8.64e07

       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: jk,t,omv= ',jk,t(jk),omv(jk)
       endif


 24     continue
 22     continue
 20      continue
!
!        Calculate mesoscale moisture-flux convergence
!
        sumq=0.
        do 25 jk=1,nlev 
	    if (jk .eq. 1) then
	       jkm=jk
            else
	       jkm=jk-1
            end if
	    if (jk .eq. nlev) then
	       jkp=jk
            else
	       jkp=jk+1
            end if
            if (pr(jk) .lt. pztm) go to 26
            if (pr(jk) .gt. pzm) go to 27
            pp=phr(jk+1)
            pm=phr(jk)
            if (phr(jk+1) .lt. pztm) pp=pztm
            if (phr(jk) .gt. pzm) pm=pzm
         qprip=(tempqa(jkp)+tempqa(jk)-q(jkp)-q(jk))/2.
         qprim=(tempqa(jk)+tempqa(jkm)-q(jk)-q(jkm))/2. 
         eqfp=ampt*qprip*alp*(phr(jk+1)-pztm)*(pzm-phr(jk+1))
!         eqfp=eqfp/(1.-ampt)
         eqfm=ampt*qprim*alp*(phr(jk)-pztm)*(pzm-phr(jk))
!         eqfm=eqfm/(1.-ampt)
         if ((phr(jk) .le. pzm) .and. (phr(jk+1) .ge. pztm)) then
            qmes(jk)=(eqfm-eqfp)/(phr(jk+1)-phr(jk))
         end if
         if ((pzm .lt. phr(jk)) .and. (pzm .gt. phr(jk+1))) then
         qmes(jk)=eqfp/(phr(jk)-phr(jk+1))
         end if
         if ((pztm .gt. phr(jk+1)) .and. (pztm .le. phr(jk))) then
           qmes(jk)=eqfm/(phr(jk+1)-phr(jk))
           end if
 27   continue
 25   continue
 26   continue
!
!     Calculate eddy flux of moist static energy in mesoscale
!     updraft and identify its minimum.
!
      hfmin=0.
      do 41 jk=1,nlev
         if (pr(jk) .lt. pztm) go to 42
         if (pr(jk) .gt. pzm) go to 43
         tmu=t(jk)+tpri
!        call establ(es,tmu)
         call lookup_es(tmu, es)
         qmu=epsilo*es/(pr(jk)+(epsilo-1.)*es)
         hflux=omv(jk)*(cpair_mul*tpri+lat*(qmu-q(jk)))
         if (hflux .lt. hfmin) then
            hfmin=hflux
            pfmin=pr(jk)
         end if
 43      continue
 41   continue
 42   continue

       if (debug_ijt) then
        print *, 'DONNER_DEEP/meens: hfmin,pfmin= ',hfmin,pfmin
       endif


!
!     Calculate mesoscale downdraft speed at pmd. Follow Leary
!     and Houze (1980,JAS) and set magnitude to half that in
!     mesoscale updraft; vertical pressure velocity constant with ht. 
!
	  omd=-alp*((pzm-pztm)**2)/8.
      omd=omd/2.
      pmd=pzm+20.e03
      if (pmd .gt. ps) pmd=ps
      if (pmd .gt. ps) pmd=ps
!
!     Calculate temperature and specific humidity in mesoscale
!     downdraft. 
!
      do 44 jk=1,nlev
         if (pr(jk) .lt. pmd) go to 45
         if (pr(jk) .gt. pb) go to 46
!
!     Calculate c2, the relative humidity in the mesoscale downdraft,
!     after Table 3 of Leary and Houze (1980, JAS).
!
         c2=1.-(.3*(pr(jk)-pmd)/(pb-pmd))
!
!     Calculate c3, the factor which yields the eddy flux of moist
!     static energy when multiplied by the minimum of moist static
!     energy in the mesoscale updraft. Multiply by 1.3 to take account
!     of convective downdrafts. See Fig. 7 of Leary and Houze
!     (1980,JAS).
!
         c3=(pr(jk)-pmd)/(pb-pmd)
         c3=1.3*c3
!
!     See "Moist Static Energy A, 1/26/91" notes.
!
         targ=t(jk)
!        call establ(es,targ)
         call lookup_es(targ, es)
         qs=epsilo*es/(pr(jk)+(epsilo-1.)*es)
         c1=epsilo*lat*es/(pr(jk)*rh2o*(t(jk)**2))
         tprimd=c3*hfmin/omd
         tprimd=tprimd-lat*(c2*qs-q(jk))
         tprimd=tprimd/(cpair_mul+lat*c1*c2)
         tempt(jk)=t(jk)+tprimd
         targ=tempt(jk)
!        call establ(es,targ)
         call lookup_es(targ, es)
         tempqa(jk)=c2*es*epsilo/(pr(jk)+(epsilo-1.)*es)

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: tprimd,tempqa,q,qs= ',tprimd,   &
					     tempqa(jk),q(jk),qs
         print *, 'DONNER_DEEP/meens: pr,rh,factr= ',pr(jk),c2,c3
       endif


 46      continue
 44   continue
 45   continue
!
!     Calculate eddy fluxes of potential temperature and specific
!     humidity in mesoscale downdraft.
!
      do 28 jk=2,nlev-1
         if (phr(jk) .lt. pmd) go to 47
         if (phr(jk) .gt. pb) go to 29
         if ((pr(jk-1) .le. pb) .and. (pr(jk) .ge. pmd)) then
            fjk=ampt*omd*((po/pr(jk))**(rair_mul/cpair_mul))  &
            *(tempt(jk)-t(jk))    
            fjkm=ampt*omd*((po/pr(jk-1))**(rair_mul/cpair_mul))*   &
            (tempt(jk-1)-t(jk-1))
            emt(jk)=(fjk+fjkm)/2.
            fjk=ampt*omd*(tempqa(jk)-q(jk))
            fjkm=ampt*omd*(tempqa(jk-1)-q(jk-1))
            emq(jk)=(fjk+fjkm)/2.
         end if
         if (pr(jk-1) .ge. pb) then
            fjk=ampt*omd*((po/pr(jk))**(rair_mul/cpair_mul))*   &
            (tempt(jk)-t(jk))
            call polat(q,pr,qb,pb, debug_ijt)
            call polat(t,pr,tb,pb, debug_ijt)
!           call establ(es,tb)
            call lookup_es(tb, es)
            qsb=epsilo*es/(pb+(epsilo-1.)*es)
            tprimd=hfmin/omd
            tprimd=tprimd-lat*(.7*qsb-qb)
            c1=epsilo*lat*es/(pb*rh2o*(tb**2))
            tprimd=tprimd/(cpair_mul+.7*lat*c1)
            fjkb=ampt*omd*((po/pb)**(rair_mul/cpair_mul))*tprimd
            wa=(phr(jk)-pr(jk))/(pb-pr(jk))
            wb=(pb-phr(jk))/(pb-pr(jk))
            emt(jk)=wa*fjkb+wb*fjk
            fjk=ampt*omd*(tempqa(jk)-q(jk))
            targ=tb+tprimd
!           call establ(es,targ)
            call lookup_es(targ, es)
            qbm=.7*epsilo*es/(pb+(epsilo-1.)*es)
            fjkb=ampt*omd*(qbm-qb)
            emq(jk)=wa*fjkb+wb*fjk
         end if
         if (pr(jk) .le. pmd) then
            fjkm=ampt*omd*((po/pr(jk-1))**(rair_mul/cpair_mul))*   &
            (tempt(jk-1)-t(jk-1))
            call polat(q,pr,qmd,pmd, debug_ijt)
            call polat(t,pr,tmd,pmd, debug_ijt)
!           call establ(es,tmd)
            call lookup_es(tmd, es)
            qsmd=epsilo*es/(pmd+(epsilo-1.)*es)
            c1=epsilo*lat*es/(pmd*rh2o*(tmd**2))
            tprimd=-lat*(qsmd-qmd)/(cpair_mul+lat*c1)
            fjkmd=ampt*omd*((po/pmd)**(rair_mul/cpair_mul))*tprimd
            wa=(pr(jk-1)-phr(jk))/(pr(jk-1)-pmd)
            wb=(phr(jk)-pmd)/(pr(jk-1)-pmd)
            emt(jk)=fjkmd*wa+fjkm*wb
            targ=tmd+tprimd
!           call establ(es,targ)
            call lookup_es(targ, es)
            qmmd=epsilo*es/(pmd+(epsilo-1.)*es)
            fjkm=ampt*omd*(tempqa(jk-1)-q(jk-1))
            fjkmd=ampt*omd*(qmmd-qmd)
            emq(jk)=fjkmd*wa+fjkm*wb
          end if

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: jk,phr,emt,emq= ',jk,phr(jk),   &
						 emt(jk),emq(jk)
       endif


          emt(jk)=((po/pr(jk))**(rair_mul/cpair_mul))*emt(jk)
 29       continue
 28   continue
 47   continue
!
!     Calculate temperature and specific humidity tendencies due
!     to eddy-flux convergences in mesoscale downdraft.
!
      rin=0.
!     do 31 jjk=1,nlev
!         jk=nlev+1-jjk
       do 31 jk=nlev,1, -1
         if ((phr(jk+1) .le. pzm) .and. (phr(jk) .ge. pzm))  &
         jksave=jk+1
         pi=(po/pr(jk))**(rair_mul/cpair_mul)
         if ((emt(jk+1) .ne. 0.) .and. (emt(jk) .eq. 0.)   &
         .and. (rin .eq. 0.)) then
            tten=-emt(jk+1)/(phr(jk+1)-ps)
            qten=-emq(jk+1)/(phr(jk+1)-ps)
            rin=1.
         end if
         if (rin .eq. 1.) then
            tmes(jk)=tmes(jk)+(tten/pi)
            qmes(jk)=qmes(jk)+qten
         end if
         if ((rin .eq. 0.) .and. (emt(jk+1) .ne. 0.) .and.   &
         (emt(jk) .ne. 0.)) then
            tten=(emt(jk+1)-emt(jk))/(phr(jk+1)-phr(jk))
            tten=-tten/pi
            qten=(emq(jk+1)-emq(jk))/(phr(jk+1)-phr(jk))
            qten=-qten
            tmes(jk)=tmes(jk)+tten
            qmes(jk)=qmes(jk)+qten
         end if

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: jk,pr,tmes,qmes= ',jk,pr(jk),  &
						 tmes(jk),qmes(jk)
       endif


 31   continue
!
!     Apply flux at top of mesoscale downdraft to level between
!     pzm and pmd, as flux at bottom of mesoscale downdraft is
!     applied to PBL.
!
      psa=0.
      do 32 jk=1,nlev
         if ((emt(jk) .ne. 0.) .and. (emt(jk+1) .eq. 0.)) then
            tten=emt(jk)/(phr(jksave)-phr(jk))
            qten=emq(jk)/(phr(jksave)-phr(jk))
            psa=phr(jk)
         end if
 32      continue

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: pmd,pb= ',pmd,pb
       endif

	 do jk=1,nlev
         if ((pr(jk) .le. psa) .and. (pr(jk) .ge. phr(jksave))) then

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: po,psa,phr(jksave)= ',po,psa,phr(jksave)
         print *, 'DONNER_DEEP/meens: jk,qmes,qten= ',jk,qmes(jk),qten
         print *, 'DONNER_DEEP/meens: jk,tmes,tten= ',jk,tmes(jk),tten
       endif

            qmes(jk)=qmes(jk)+qten
            pi=(po/pr(jk))**(rair_mul/cpair_mul)
            tmes(jk)=tmes(jk)+(tten/pi)
         end if
	 end do
      cmui=0.
      wmc=0.
      wpc=0.
      owms=0.
      do 15 jk=1,nlev
      tmes(jk)=tmes(jk)*86400.
	  qmes(jk)=qmes(jk)*8.64e07
      wmm(jk)=wmm(jk)*8.64e07
      owm(jk)=owm(jk)*8.64e07
      wmp(jk)=wmp(jk)*8.64e07
      wmc=wmc+wmm(jk)*(phr(jk)-phr(jk+1))
      owms=owms+owm(jk)*(phr(jk)-phr(jk+1))
      wpc=wpc+wmp(jk)*(phr(jk)-phr(jk+1))
      cmui=cmui+cmu(jk)*(phr(jk)-phr(jk+1))
 15   continue
      wmc=wmc/(gravit*1000.)
      wpc=wpc/(gravit*1000.)
      owms=owms/(gravit*1000.)
      cmui=cmui/(gravit*1000.)

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: wmc=',wmc,' mm/day',' wpc=',wpc, 'mm/day'
         print *, 'DONNER_DEEP/meens: owms= ',owms,' mm/day',' cmui= ',   &
							 cmui,'mm/day'
       endif


! 
!
!     calculate mesoscale precipitation
!
      cmui=cmui-wmc
      rm=gnum*(cmui+ca)
      contot=rc/(rm+rc)

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: cmui,ca=', cmui,ca
         print *, 'DONNER_DEEP/meens: rm= ',rm,'rc= ',rc
       endif
!
!      calculate integrated water evaporated in mesoscale  
!      downdraft using Leary and Houze coefficients
!
      a=.4
      emdi=a*(cmui+ca)
      b=.1
      emei=b*(cmui+ca)
!
!     calculate vertical tendency distributions (g/kg/day)
!
!      evaporation from mesoscale updraft
!
      p1=pzm
      emea=emei*gravit*1000./(p1-pztm)

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: emea,emei= ',emea,emei
       endif


      if (p1 > pztm) then
      call ver(emea,p1,pztm,pr,phr,eme, debug_ijt)
!     call ver(emea,p1,pztm,   phr,eme, debug_ijt)
      endif
!
!
      p3=pzm
      pst=pzm+30.e03
      pst=ps
      p1=pst
      emda=emdi*gravit*1000./(p1-p3)
!
!      evaporation in mesoscale downdraft
!
       do 35 jk=1,nlev
          if (phr(jk) .lt. p3) go to 36
          if (phr(jk+1) .gt. pst) go to 37
          pm=phr(jk)
          pp=phr(jk+1)
          if ((phr(jk) .ge. pst) .and. (phr(jk+1) .le.   &
          pst)) pm=pst
          if ((phr(jk) .ge. p3) .and. (phr(jk+1) .le. p3))  &
          pp=p3
          emd(jk)=emda*(pm-pp)*(pm+pp-2.*pst)

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: jk,emda,emd= ',jk,emda,emd(jk)
       endif

          emd(jk)=emd(jk)/((phr(jk)-phr(jk+1))*(p3-pst))
 37   continue
 35   continue
 36   continue

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: emdi,ps= ',emdi,ps
       endif


!      emda=.75*emda
!      if (p1 > p3) then
!!     call ver(emda,p1,p3,pr,phr,emd, debug_ijt)
!      call ver(emda,p1,p3,   phr,emd, debug_ijt)
!      endif
!
!      evaporation in mesoscale downdraft
!

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: emdi,ps= ',emdi,ps
       endif


!      p3=pzm   
!      p1=pzm+15.e03
!      emda=emdi*gravit*1000./(p1-p3)
!      emda=.25*emda
!      if (p1 > p3) then
!!     call ver(emda,p1,p3,pr,phr,emdx, debug_ijt)
!      call ver(emda,p1,p3,   phr,emdx, debug_ijt)
!      endif
!      delp=pb-pt
!      p1=pt+.50*delp
!      emd1=emdi*.2
!      emda=emd1*gravit*1000./(p1-pzm)
!      if (p1 > pzm) then
!!     call ver(emda,p1,pzm,pr,phr,emdx, debug_ijt)
!      call ver(emda,p1,pzm,   phr,emdx, debug_ijt)
!      endif
!      do 2 i=1,nlev
!      emd(i)=emd(i)+emdx(i)
!      emdx(i)=0.
! 2    continue
!      emd1=0.5*emdi
!      p1=pt+.25*delp

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: pzm= ',pzm
       endif


!      do 5 i=1,nlev
!      emd(i)=0.
! 5    continue
!      factsum=0.
!      do 4 i=1,nlev
!      if ( ( (phr(i) .le. pst) .and. (phr(i) .ge. p3) )
!     1 .or. ( (phr(i+1) .le. pst) .and. (phr(i+1) .ge. p3) ) ) then
!      if (phr(i) .le. pst) facta=phr(i)
!      if (phr(i) .gt. pst) facta=pst
!      if (phr(i+1) .ge. p3) factb=phr(i+1)
!      if (phr(i+1) .lt. p3) factb=p3
!      fact=(facta-factb)*(pst-.5*(facta+factb))
!      fact=2.*fact/((p3-pst)**2)
!      write(6,*) 'i,fact= ',i,fact
!      factsum=factsum+fact
!      emd1=(emdi)*fact 
!      emda=emd1*gravit*1000./(facta-factb)
!      if (facta > factb) then
!!     call ver(emda,facta,factb,pr,phr,emdx, debug_ijt)
!      call ver(emda,facta,factb,   phr,emdx, debug_ijt)
!      endif
!      do 3 j=1,nlev
!      emd(j)=emd(j)+(emdx(j)*emdi/(emdi))
!      write(6,*) 'i,pr,emd= ',i,pr(i),emd(i)
! 3    continue
!      end if
! 4    continue

!      if (debug_ijt) then
!        print *, 'DONNER_DEEP/meens: factsum= ',factsum
!      endif


!
!     equivalent for melting of mesoscale sfc precip
!     generated by deposition in mesoscale updraft
!
      p2=-10.
      do 51 j=1,nlev-1
	 if (phr(j+1) < pztm ) go to 52
         if ((t(j) .ge. tmel) .and. (t(j+1) .le. tmel))   &
          then
            p2=phr(j+1)
            go to 52
	  end if
 51   continue
 52   continue
      if (p2 .ne. -10.) then
	  rm=(rc/contot)-rc
      else
	  rm=0.
      end if
      rmm=rm
      rma=rmm*gravit*1000./(pb-p2)
      if (pb > p2) then
      call ver(rma,pb,p2,pr,phr,elt, debug_ijt)
!     call ver(rma,pb,p2,   phr,elt, debug_ijt)
      endif
!
!     Convert cmui and emei from mm/day to kg(H2O)/((m**2) sec)
      cmui=cmui/86400.
      emei=emei/86400.
      do j=1,ncap
	 if ( (p(j) .le. pzm) .and. (p(j) .ge. pztm)) cumh(j)=ampt
	 if ( (p(j) .le. pb) .and. (p(j) .ge. pmd) ) dmemh(j)=    &
         -omd*ampt/gravit
      end do
	call verav(cumh,pb,pztm,thetl,cuml,sbl,ps,pr,phr,cappa, debug_ijt)
       call verav(dmemh,pb,pmd,thetl,dmeml,sbl,ps,pr,phr,cappa, debug_ijt)
      do j=1,nlev
	 umeml(j)=-omv(j)*ampt/gravit
      end do

       if (debug_ijt) then
         print *, 'DONNER_DEEP/meens: pzm,pztm,pb,p2= ',pzm,pztm,pb,p2
          do k=1,nlev-ktest_model+1
	 if (omv(k) /= 0.00) then
         print *, 'DONNER_DEEP/meens: j,omv= ',k,omv(k)
	 endif
	 if (elt(k) /= 0.0) then 
         print *, 'DONNER_DEEP/meens: j,elt= ',k,elt(k)
	 endif
           end do
       endif



end subroutine meens





subroutine mesub(cu,rc,dint,pr,phr,ps,pb,pt,gravit    &
!subroutine mesub(cu,rc,dint,pr,phr,ps,pb,pt           &
!      ,cappa,epsilo,rair_mul,t,ca,tmel   &
       ,cappa,       rair_mul,t,ca,tmel   &
       ,ecd,ece,fre,elt, debug_ijt)

!--------------------------------------------------------------------
real,               intent(inout) :: dint
real,               intent(in) :: cu, rc,       ps, pb, pt, gravit, &
!real,               intent(in) :: cu, rc,       ps, pb, pt,         &
!                                   cappa, epsilo, rair_mul, tmel
                                    cappa,         rair_mul, tmel
real, dimension(:), intent(in) :: pr, phr, t
logical,            intent(in) :: debug_ijt
real, dimension(:), intent(inout) :: ecd, ece, fre, elt
real,               intent(inout) :: ca
!---------------------------------------------------------------------

!
!     Calculate mesoscale heat and moisture sources, using
!     variation on Leary and Houze (JAS, 1980).
!     Performs calculations related to individual subensembles only.
!     For notation, see "Cu Closure A notes," 2/97
!
!     On Input:
!
!       cu      condenstation integral
!               sigma(j=1,M)((r(i,j)**2)/(r(i,p_b)**2))*c_u(i,j)
!               for ensemble i
!               (mm/day)
!       rc      precipitation integral
!               sigma(j=1,M)((r(i,j)**2)/(r(i,p_b)**2))*r_c(i,j)
!               for ensemble i
!               (mm/day)
!       pr,phr  low-resolution pressure levels (Pa)
!       ps      surface pressure (Pa)
!       pb      cloud-base pressure (Pa)
!       pt      cloud-top pressure (Pa)
!       gravit  gravity constant (m/s**2)
!       cappa   ratio of gas constant to specific heat for dry air
!       epsilo  ratio of gas constants, dry air to water vapor
!       rair    gas constant for dry air (J/kg/K)
!       t       low-resolution temperature (K)
!       tmel    melting temperature (K)
!       dint    water mass frozen in convective updraft
!               plus ice deposited convective updraft
!               (kg(H2O)/((m**2) sec)
!               weighted as cu,rc
!
!     On Output:
!
!       ca      condensed water X-fer from cells to anvil (mm/day)
!               weighted as cu,rc
!       ecd     water mass evaporated in convective
!               downdraft (g/kg/day) (vertical integral ecdi
!               weighted as rc,cu)
!       ece     water mass evaporated from convective updraft
!               (g/kg/day) (vertical integral ecei weighted as cu,rc)
!       fre     equivalent for freezing in mesoscale updraft
!               (g/kg/day) (vertical integral caa weighted as cu,rc)
!       elt     equivalent for melting in mesoscale downdraft
!               (g/kg/day) (vertical integral elt weighted as cu,rc)
!


!
!      define constants
!
      dp=-1000.
      ptt=pt+dp
      pztm=ptt-30.e03
!     Restrict pztm to .ge. 10 kPa, cf Ackerman et al (JAS,1988)
!     (unless pt .le. 10 kPa)
!
      if (pztm .lt. 10.e03) pztm=10.e03
      if (ptt .lt. 10.e03) pztm=ptt+dp
      pzm=ptt
      if (pzm .le. pztm) pzm=pztm-dp
!
!      water-budget constants calculated from flux-
!      condensation parameterization
! 

       if (debug_ijt) then
         print *, 'DONNER_DEEP/mesub: rc,cu= ',rc,cu
       endif

      gnu=rc/cu
!
!     maintain Leary and Houze ratios of alpha, beta, and
!     eta to their sum
!
      sum=1.-gnu
      alrat=.25
      berat=.13
      etrat=.62
      alpha=alrat*sum
      beta=berat*sum
      eta=etrat*sum
!
!     use Leary and Houze ratios for gnu-sub-m
!
      gnum=.5

       if (debug_ijt) then
         print *, 'DONNER_DEEP/mesub: gnu= ',gnu
       endif


!
!     calculate mass of cu incorporated into mesoscale anvil,
!     integrated water evaporated in convective downdraft,
!     and integrated water evaporated from convective updraft
!


      ca=eta*cu


       if (debug_ijt) then
         print *, 'DONNER_DEEP/mesub: alpha,beta,eta= ',alpha,beta,eta
         print *, 'DONNER_DEEP/mesub: ca= ',ca
       endif

      ecei=beta*cu
      ecdi=alpha*cu
!
!
!     calculate vertical tendency distributions (g/kg/day)
!
!      calculate equivalent for freezing in mesoscale updraft
!
      dint=dint*8.64e07*gravit/(pzm-pztm)
      if (dint .ne. 0.)    &
       caa=(ca+ecei)*gravit*1000./(pzm-pztm)
      if (dint .eq. 0.)   &
       caa=ca*gravit*1000./(pzm-pztm)

       if (debug_ijt) then
         print *, 'DONNER_DEEP/mesub: caa,dint=',caa,dint
       endif


      caa=caa-dint

       if (debug_ijt) then
      if (caa .gt. 0.) then
         print *, 'DONNER_DEEP/mesub: caa,pzm,pztm= ',caa,pzm,pztm
          do i=1,nlev-ktest_model+1
	   if (fre(i) /= 0.0) then
         print *, 'DONNER_DEEP/mesub: i,fre= ',i,fre(i)
	   endif
          end do
       endif
       endif

      if ( pzm > pztm) then
      if (caa .gt. 0.)  call ver(caa,pzm,pztm,pr,phr,fre, debug_ijt)
!     if (caa .gt. 0.)  call ver(caa,pzm,pztm,   phr,fre, debug_ijt)
      endif
!
!     evaporation in convective downdraft
!
      pz0=ptt
      p1=ps
      ecda=ecdi*gravit*1000./(p1-pz0)

       if (debug_ijt) then
         print *, 'DONNER_DEEP/mesub: ecda,p1,pz0= ',ecda,p1,pz0
          do i=1,nlev-ktest_model+1
	 if (ecd(i) /= 0.0) then
         print *, 'DONNER_DEEP/mesub: i,ecd= ',i,ecd(i)
	 endif
          end do
       endif


      if (p1 > pz0)then
      call ver(ecda,p1,pz0,pr,phr,ecd, debug_ijt)
!     call ver(ecda,p1,pz0,   phr,ecd, debug_ijt)
      endif
!
!     calculate melting due to excess freezing, over
!     ecei and ca
!
      nlevm=nlev-1
      p2=0.
      do 51 jk=1,nlev
	 if (phr(jk+1) < pztm) exit
         if ((t(jk) .ge. tmel) .and. (t(jk+1) .le. tmel))   &
         then
         p2=phr(jk+1)
         go to 52
         end if
 51   continue
 52   continue
      elta=0.
      if (caa .le. 0.) then 
         caa=-caa*(pzm-pztm)/(pb-p2)
         elta=caa
      end if
      if (p2 .eq. 0.) elta=0.

       if (debug_ijt) then
         print *, 'DONNER_DEEP/mesub: elta,pb,p2= ',elta,pb,p2
          do i=1,nlev-ktest_model+1
	 if (elt(i) /= 0.0) then
         print *, 'DONNER_DEEP/mesub: i,elt= ',i,elt(i)
	 endif
          end do
       endif


      if (pb > p2) then
      call ver(elta,pb,p2,pr,phr,elt, debug_ijt)
!     call ver(elta,pb,p2,   phr,elt, debug_ijt)
      endif
!
!     evaporation from convective updraft
!
      p1=pt+5.0e03
      p2=ptt
      ecea=ecei*gravit*1000./(p1-p2)

       if (debug_ijt) then
         print *, 'DONNER_DEEP/mesub: ecea,ecei= ',ecea,ecei
         print *, 'DONNER_DEEP/mesub: ecea,p1,p2= ',ecea,p1,p2
          do i=1,nlev-ktest_model+1
	 if (ece(i) /= 0.0) then
         print *, 'DONNER_DEEP/mesub: i,ece= ',i,ece(i)
	 endif
          end do
       endif


      if (p1 > p2) then
      call ver(ecea,p1,p2,pr,phr,ece, debug_ijt)
!     call ver(ecea,p1,p2,   phr,ece, debug_ijt)
      endif


end subroutine mesub



!####################################################################

subroutine micro(tc1,tc2,p1,p2,te1,te2,qe1,qe2,w1,w2, rr,rmu,  &
         qrw,qcw,qlw,dcw1, dqrw3, debug_ijt)

 !--------------------------------------------------------------------
real,      intent(inout) :: qcw, qrw, rmu
real,        intent(in) :: tc1, tc2, p1, p2, te1, te2, qe1, qe2, &
                             w1, w2, rr
real,        intent(inout) :: qlw, dcw1, dqrw3
logical,      intent(in) :: debug_ijt
!--------------------------------------------------------------------




!
!     Kessler microphysics
!
!     On Input:
!        tc1    cloud temperature (K) at pressure p1 (Pa)
!        tc2    cloud temperature (K) at pressure p2 (Pa)
!        te1    environmental temperature (K) at pressure p1 (Pa)
!        te2    environmental temperature (K) at pressure p2 (Pa)
!        qe1    environmental mixing ratio (kg(H2O)/kg) at pressure p1
!        qe2    environmental mixing ratio (kg(H2O)/kg) at pressure p2
!        w1     cloud vertical velocity (m/s) at pressure p1
!        w2     cloud vertical velocity (m/s) at pressure p2
!        rr     cloud radius (m)
!        qrw    rain water (g(H2O)/m**3)
!        qcw    cloud water (g(H2O)/m**3)
!        rmu    entrainment coefficient (/s)
!
!     On Output:
!        qlw    total liquid water (kg(H2O)/kg)
!        dcw1   condensation increment (g(H2O)/m**3)
!        dqrw3  condensation increment (g(H2O)/m**3)
!
!     CONSTANTS
!

        if (debug_ijt) then
          print *, 'DONNER_DEEP/micro: qrw,qcw= ',qrw,qcw
         print *,  'DONNER_DEEP/micro: rr= ',rr
	endif

      ep=.622
      g=9.81
      rd=287.05
      alp=.5
      dp=-1000.
!
!     epm is defined in "Generalized mu 10/1/89"
!
      epm=0.
      hr=rr
!
!     call establ(es1,tc1)
!     call establ(es2,tc2)
      call lookup_es(tc1, es1)
      call lookup_es(tc2, es2)
      rs1=ep*es1/(p1-es1)
      rs2=ep*es2/(p2-es2)
      tcb=(tc1+tc2)/2.

        if (debug_ijt) then
          print *, 'DONNER_DEEP/micro: rs1,rs2,p1,p2= ',rs1,rs2,p1,p2
          print *, 'DONNER_DEEP/micro: es1,es2,tc1,tc2= ',es1,es2,tc1,tc2
          print *, 'DONNER_DEEP/micro: te1,te2= ',te1,te2
	endif


      qcb=(rs1+rs2)/2.
      d1=rd*(1.+alp)*tc1*te1/(g*p1*(alp*te1+tc1))
      d2=rd*(1.+alp)*tc2*te2/(g*p2*(alp*te2+tc2))
      dz=-(d1+d2)*dp/2
!
!     calculate condensation
!
      qeb=(qe1+qe2)/2.

        if (debug_ijt) then
          print *, 'DONNER_DEEP/micro: qeb,rmu= ',qeb,rmu
	endif


      rmusa=rmu
!      rmu=0.
      cond=rs2*exp(rmu*dz)-rs1+(1.-exp(rmu*dz))*qeb
      cond=-cond/(1.-epm+epm*exp(rmu*dz))
      if (cond .lt. 0.) cond=0.

        if (debug_ijt) then
          print *, 'DONNER_DEEP/micro: cond= ',cond
	endif


      dcw1 =(cond)*(p1+p2)*500./(tcb*rd*(1.+.61*qcb))
      w=(w1+w2)/2.
      IF (QCW .GE. .5) DCW2=DZ*(QCW-.5)/W
      IF (QCW .GE. .5) DCW2=DCW2*1.E-03
      IF (QCW .LT. .5) DCW2=0.
      DQCW3=0.
      IF ((QCW .EQ. 0.) .OR. (QRW .EQ. 0.)) GO TO 7
      DQCW3=5.26E-03*QCW*(QRW**.875)*DZ/W
!
!     calculate effect of entrainment on cloud water
!
 7    continue
      ent=rmu*dz
      qcw=qcw/exp(ent)
      qcwa=qcw+dcw1-dcw2-dqcw3
      if (qcwa .lt. 0.) then
      red=(qcw+dcw1)/(dcw2+dqcw3)
         dcw2=dcw2*red
         dqcw3=dqcw3*red
      end if
      QCW=QCW+DCW1-DCW2-DQCW3

        if (debug_ijt) then
          print *, 'DONNER_DEEP/micro: ent,qcw,dcw1= ',ent,qcw,dcw1
	endif


      IF (QCW .LT. 0.) QCW=0.
      DQRW3=0.
      IF (QRW .EQ. 0.) GO TO 8
      DQRW3=(QRW**1.125)*5.1*DZ/(W    *HR)
!
!     calculate effect of entrainment on rain water
!
 8    continue
      qrw=qrw/exp(ent)
      qrwa=qrw+dcw2+dqcw3-dqrw3
      if (qrwa .lt. 0.) then
         red=(qrw+dcw2+dqcw3)/dqrw3
         dqrw3=red*dqrw3
      end if
      QRW=QRW+DCW2+DQCW3-DQRW3

        if (debug_ijt) then
          print *, 'DONNER_DEEP/micro: ent,qrw= ',ent,qrw
	endif


      IF (QRW .LT. 0.) QRW=0.
      QLW=QRW+QCW

        if (debug_ijt) then
          print *, 'DONNER_DEEP/micro: exit micro dcw1,dqrw3,qlw= ',   &
				   dcw1,dqrw3,qlw
	endif


!     QLW IN UNITS OF G/M**3
      TV=TCB*(1.+.61*((RS1   +RS2     )/2.))
      QEB=(QE1  +QE2    )/2.
      TEB=(TE1  +TE2    )/2.
      TVE=TEB*(1.+.61*QEB)
      QLW=2.0E-03*QLW*TV*RD/(P1+p2    )
      rmu=rmusa



!-------------------------------------------------------------------

end  subroutine micro




subroutine clotr(alp, ent, gravit, p1, p2, rd, sou1, sou2, tc1, tc2, &
                       te1, te2, xe1, xe2, xc, debug_ijt)

!-----------------------------------------------------------------------
!     Calculates in-cloud tracer profile.
!
!     Leo Donner
!     GFDL
!     14 Jan 00
!
!-----------------------------------------------------------------------
!
!     IMPLICIT NONE
!
!---------------------Parameters----------------------------------------
!
      integer ncont            !  number of tracers
      parameter(ncont=1)
!
!
!---------------------Arguments-----------------------------------------
!
!     Input arguments
!
      real alp                 !  virtual mass coefficient (dimensionless)
      real ent                 !  entrainment coefficient (m**-1)
      real gravit              !  gravity constant (m/s2)
      real p1                  !  pressure at level nearer earth surface (Pa)
      real p2                  !  pressure at level farther from surface
			       !  (Pa)
      real rd                  !  gas constant (J/kg K)
      real sou1(ncont)         !  in-cloud source of x at pressure p1 (Pa)
      real sou2(ncont)         !  in-cloud source of x at pressure p2 (Pa) 
      real tc1                 !  cloud temperature (K) at pressure p1 (Pa)
      real tc2                 !  cloud temperature (K) at pressure p2 (Pa)
      real te1                 !  environmental temperature (K) at pressure  
			       !  p1 (Pa)
      real te2                 !  environmental temperature (K) at pressure  
			       !  p2 (Pa)
      real xe1(ncont)          !  large-scale tracer concentration at
			       !  pressure p1 (Pa)
      real xe2(ncont)          !  large-scale tracer concentration at
			       !  pressure p2 (Pa)

      logical debug_ijt
!
!     Input/output arguments
!
      real xc(ncont)           !  cloud-tracer concentration at pressure
			       !  p1 (Pa) on input, p2 (Pa) on output
!
!---------------------Local Workspace-----------------------------------
!
      real d1                  !  -dz/dp from p1 side  (m/Pa) 
      real d2                  !  -dz/dp from p2 side (m/Pa)
      real dz                  !  height increment (m)
      real epm                 !  defined in "Generalized mu, 10/1/89"
      integer kcont            !  counter over tracers
      real seb                 !  layer average of sou
      real xeb                 !  layer average of xe
!
!     NOTES:  xc and xe must have same units.
!             sou must have the (units of xc and xe)/sec.
!
!-----------------------------------------------------------------------
      integer ktest
      integer itest
      integer jtest
      integer ktdiag
      integer idiag
      integer jdiag
      real tediag
!#include "cudiag.H"

      if (debug_ijt) then
	print *, 'DONNER_DEEP/clotr: entering clotr'
      endif


!
!     Initialization.
!
      epm=0.
!
!     Calculate in-cloud profile of tracer.
!
        d1=rd*(1.+alp)*tc1*te1/(gravit*p1*(alp*te1    &
         +tc1))
        d2=rd*(1.+alp)*tc2*te2/(gravit*p2*(alp*te2   &
         +tc2))
	dz=(d1+d2)*(p1-p2)/2.
	do kcont=1,ncont
	   xeb=(xe1(kcont)+xe2(kcont))/2.
	   seb=(sou1(kcont)+sou2(kcont))/2.
	   xc(kcont)=xc(kcont)/exp(ent*dz)   
	   xc(kcont)=xc(kcont)+(((exp(ent*dz)-1.)*xeb)/    &
           exp(ent*dz))
	   xc(kcont)=xc(kcont)+(seb*(1.+epm*(exp(ent*dz)-1.))   &
                /exp(ent*dz))
	   if (xc(kcont) .lt. 0.) xc(kcont)=0.
	end do

      if (debug_ijt) then
	print *, 'DONNER_DEEP/clotr: xc= ',xc(ncont)
	print *, 'DONNER_DEEP/clotr: d1,d2,dz= ',d1,d2,dz
	print *, 'DONNER_DEEP/clotr: xeb= ',xeb
	print *, 'DONNER_DEEP/clotr: seb= ',seb
	print *, 'DONNER_DEEP/clotr: ent= ',ent
      endif




end subroutine clotr




!subroutine tae(t,p,q,latvap,epsilo,dp,rair_mul,cpair_mul,ta)
!subroutine tae(t,p,q, lat,  epsilo,dp,rair_mul,cpair_mul, cappa, ta)
subroutine tae(t,p,q, lat,         dp,rair_mul,cpair_mul, cappa, ta)

!real,   intent(in) :: t, p, q, latvap, epsilo, dp, rair_mul, cpair_mul
!real,    intent(in) :: t, p, q,  lat,epsilo, dp, rair_mul, cpair_mul, &
real,     intent(in) :: t, p, q,  lat,        dp, rair_mul, cpair_mul, &
			 cappa
real,     intent(out)  :: ta

!
!     calculates adiabatic equivalent temperature
!
!     On Input:
!        t     temperature (K)
!        p     pressure (Pa)
!        q     specific humidity (kg(H2O)/kg)
!              various constants
!
!     On Output:
!        ta    adiabatic equivalent temperature
!


!
!     define constants
!
!     cappa=rair_mul/cpair_mul
!
      ta=t
!     do 1 i=1,nc
      do 1 i=1,ncap
      pr=p+(i-1)*dp
      if (pr .le. 0.) go to 2
      te=t*((pr/p)**cappa)
!     call establ(es,te)
      call lookup_es(te, es)
      qe=epsilo*es/(p-es)
      if (q .ge. qe) then
!        ta=t*exp(latvap*q/(te*cpair_mul))
         ta=t*exp(lat*q/(te*cpair_mul))
         go to 2
      end if
 1    continue
 2    continue

end subroutine tae





subroutine verav(qrnh,pb,pt,thetl,qrnl,qbl,ps,pr,phr,cappa, debug_ijt)


!--------------------------------------------------------------------
logical,      intent(in) :: debug_ijt, thetl
real, dimension(:), intent(in) :: qrnh, pr, phr
real,         intent(inout)  ::  qbl
real,         intent(in)  :: pb, pt,      ps, cappa
real, dimension(:), intent(inout) :: qrnl
!--------------------------------------------------------------------

!
!      layer averaging for large-scale       source due to cumulus
!      convection
!
!      on input:
!        qrnh(ncm) large-scale       source at cloud-model resolution
!        qbl  large-scale PBL    source due to Cu
!        pb        cloud-base pressure (Pa)
!        pt        cloud-top pressure (Pa)
!        ps        surface pressure (Pa)
!        thetl     indicator for temperature calculation
!        pr        large-scale pressure levels (Pa)
!        phr       large-scale pressure half levels (Pa)
!        cappa     constant
!
!     on output:
!
!        qrnl(nlev)large-scale       source at CCM resolution
!

real, dimension (ncap) :: p
logical                :: test2








      test2=.false.

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: entering verav'
      endif


      dp=-1000.
      ptt=pt+dp
!BUG??     do 1 i=1,ncap-1
      do 1 i=1,ncap
      p(i)=pb+(i-1)*dp
 1    continue
      do 8 i=1,nlev
      qrnl(i)=0.
  8   continue
      j1=1
      rintsum=0.
      do 2 i=1,nlev
      ph=phr(i+1)
      pl=phr(i)

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: pl,ph,pb,pt',pl,ph,pb,pt
      endif


      if (ph .ge. pb) go to 30

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: thetl= ',thetl
      endif


      if (pl .le. ptt) go to 7
      rint=0.0
      rkou=0.
         do 4 j=j1,ncap-1
         phrh=(p(j)+p(j+1))/2.
         phrl= p(j)-(dp/2.)
            if (j .eq. 1 .and. (pl .ge. pb)) rkou=pl-pb
         if (phrl .gt. pl) phrl=pl
         if (phrl .gt. pb) phrl=pb
         if (phrh .lt. ph) phrh=ph
         if (phrh .le. ptt) then
            rkou=rkou+phrl-ph
            rint=rint+qrnh(j)*(phrl-ptt)
            qrnl(i)=rint/rkou
            go to 7
         end if

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: j,phrl,phrh= ',j,phrl,phrh
      endif


         rkou=rkou+phrl-phrh
         rint=rint+qrnh(j)*(phrl-phrh)
         rintsum=rintsum+qrnh(j)*(phrl-phrh)

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: rintsum= ',rintsum
      endif


         j1=j

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: j,qrnh,p,rkou',j,qrnh(j),p(j),rkou
	print *, 'DONNER_DEEP/verav: rint= ',rint
      endif


         if (phrh   .le. ph       ) then
            rint=rint/rkou
            qrnl(i)=rint
            go to 30
         end if
 4       continue
 30   continue

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: i,pl,ph,qrnl(i)',i,pl,ph,qrnl(i)
      endif


 2    continue
 7    continue
      if (thetl) then
      qlsum=0.
      qlsu=0.
      do 9 i=1,nlev

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: i,pr,thetl= ',i,pr(i),thetl
      endif


         pi=(1.0e05/pr(i) )**cappa
         qlsum=qlsum+qrnl(i)*pi*(phr(i)-phr(i+1))
         if (phr(i) .le. pb)    &
         qlsu=qlsu+qrnl(i)*pi*(phr(i)-phr(i+1))
 9    continue
      qbl0=qbl
      qbl=-qlsum/(ps-pb)

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: qbl= ',qbl
	print *, 'DONNER_DEEP/verav: qlsu=',qlsu
	print *, 'DONNER_DEEP/verav: thetl qlsum= ',qlsum
      endif


      end if

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: first qlsum= ',qlsum
      endif


      sb=pb/ps

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: ps= ',ps
      endif


      do 3 i=1,nlev

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: i,qbl,phr,pr= ',i,qbl,phr(i+1),pr(i)
      endif


      if (phr(i+1) .ge. pb) then
         qrnl(i)=qbl
         if (thetl) then
         pi=(1.0e05/pr(i))**cappa
         qrnl(i)=qbl/pi 
         qrnl(i)=qrnl(i)+qbl0
         end if
      end if
      if (phr(i+1) .lt. pb) then

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: phr,phr+,qrnl= ',phr(i),phr(i+1),qrnl(i)
      endif


         if (phr(i) .le. pb) go to 5
         wta=qrnl(i)*(phr(i)-phr(i+1))
         wtb=qbl*(phr(i)-pb)
         qrnl(i)=(wta+wtb)/(phr(i)-phr(i+1))
          if (thetl) then
            pi=(1.0e05/pr(i))**cappa
          thetsum=qlsum*(ps-phr(i))/(pb-ps)

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: qbl0,qlsum,qlsu,thetsum= ',qbl0,qlsum,  &
                                                      qlsu,thetsum
      endif


            qrnl(i)=(qlsu)+thetsum
            qrnl(i)=-qrnl(i)
            qrnl(i)=qrnl(i)/(pi*(phr(i)-phr(i+1)))
!!! BUG FOUND -- BAD INDEX !!
!           qrnl(i)=qrnl(i) + ( (qbl0/ ( (phr(i)-phr(i-1))))   &
            qrnl(i)=qrnl(i) + ( (qbl0/ ( (phr(i)-phr(i+1))))   &
            *(pb-phr(i) ))
          end if
      end if
 3    continue
 5    continue
      qlsum=0.
      do 11 i=1,nlev

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: i,phr,phr+= ',i,phr(i),phr(i+1) 
	print *, 'DONNER_DEEP/verav: pr= ',pr(i),'qrnl= ',qrnl(i)
      endif


         pi=(1.0e05/pr(i))**cappa
         if (.not. thetl) qlsum=qlsum+qrnl(i)*(phr(i)-phr(i+1))
         if (thetl) qlsum=qlsum+qrnl(i)*pi*(phr(i)-phr(i+1))
 11   continue

      if (debug_ijt) then
	print *, 'DONNER_DEEP/verav: qlsum= ',qlsum
      endif


end subroutine verav



subroutine ver (xav, p1, p2, pr, phr, x, debug_ijt)
!subroutine ver (xav, p1, p2,     phr, x, debug_ijt)

!------------------------------------------------------------------
 logical,   intent(in) :: debug_ijt
 real,     intent(in) :: xav,p1, p2
  real, dimension(:), intent(in) :: phr
!  real, dimension(:), intent(in) :: pr
!!! COMPILER BUG ?????
!!! pr must be dimensioned below to work properly ???  
 real, dimension(:), intent(out) :: x
!------------------------------------------------------------------

!
!     vertical averaging subroutine for mesoscale tendencies
!
!     On Input:
!       xav     average tendency
!       p1      high pressure boundary(Pa)
!       p2      low pressure boundary(Pa)
!       pr,phr  low-resolution pressure levels
!
!     On Output:
!       x       vertically averaged tendency (g/kg/day)
!

!      dimension         pr(nlev)               




!******************************************************************


!---- NOTICE:
!
! THE FOLLOWING LINE IS NEEDED, BUT NEED NOT BE CHANGED IF VERTICAL|
! RESOLUTION IS CHANGED
! IT LIKELY REFLECTS THE PRESENCE OF A BUG SOMEWHERE WHICH HAS YET TO
! BE LOCATED
!
!
       dimension         pr(127 )               
!
!
!
!---- END OF NOTICE



!******************************************************************


!      dimension         pr( size(phr)-1 )               

    

      if (debug_ijt) then
       print *,  'DONNER_DEEP/ver: xav,p1,p2= ',xav,p1,p2
      endif

      if ( p1 < p2 ) then
       call error_mesg ('DONNER_DEEP/ver',  &
           ' input pressure p1 is less than input pressure p2', FATAL)
      endif



!
      do 1 i=1,nlev
        if (p1 .le. phr(i+1)) x(i)=0.
        if (p2 .ge. phr(i)) x(i)=0.
        if ((p1 .ge. phr(i)) .and. (p2 .le. phr(i+1)) )   &
        x(i)=xav
        if ( (p1 .le. phr(i)) .and. (p1 .ge. phr(i+1))     &
        .and. (p2 .le. phr(i+1)) ) x(i)=xav*(p1-phr(i+1))/  &
         (phr(i)-phr(i+1))
        if ( (p2 .ge. phr(i+1)) .and. (p2 .le. phr(i))   &
        .and. (p1 .ge. phr(i)) ) x(i)=xav*(phr(i)-p2)/  &
        (phr(i)-phr(i+1))
      if ((p1 .le. phr(i)) .and. (p2 .ge. phr(i+1))) x(i)=   &
        xav*(p1-p2)/(phr(i)-phr(i+1))
 1    continue

      if (debug_ijt) then
      do i=1,nlev
       if (x(i) /= 0.0) then
         print *,  'DONNER_DEEP/ver: i,x= ',i,x(i)
       endif
      end do
      endif



end subroutine ver




subroutine polat (xv, pv,x,p, debug_ijt)

!---------------------------------------------------------------------
logical,            intent(in) :: debug_ijt
real, dimension(:), intent(in) :: xv, pv
real,              intent(in)  :: p
real,               intent(inout)  :: x
!---------------------------------------------------------------------



!
!     INTERPOLATES XV TO PRESSURE P.
!
!     ON INPUT:
!
!         XV(N)   DATA AT RESOLUTION N
!         PV(N)   PRESSURE AT N LEVELS
!
!     ON OUTPUT:
!
!         X       DATA AT PRESSURE P
!



       if (debug_ijt) then
	 print *, 'DONNER_DEEP/polat: p=', p
       endif


       IF (P .GE. PV(1)) THEN

       if (debug_ijt) then
	 print *, 'DONNER_DEEP/polat: p= ', p
       endif


          X=(XV(2)-XV(1))/(PV(2)-PV(1))
          X=X*(P-PV(1))+XV(1)
          GO TO 1
       ENDIF
       IF (P .LE. PV(nlev)) THEN

       if (debug_ijt) then
	 print *, 'DONNER_DEEP/polat: p= ', p
       endif


          X=(XV(nlev)-XV(nlev-1))/(PV(nlev)-PV(nlev-1))
          X=X*(P-PV(nlev))+XV(nlev)
          GO TO 1
       ENDIF
       DO 10 I=1,nlev-1
       IF ((P .GE. PV(i+1)) .AND. (P .LE. PV(I))) THEN

       if (debug_ijt) then
	 print *, 'DONNER_DEEP/polat: p= ', p
	 print *, 'DONNER_DEEP/polat: i,xv(i),xv(i+1)= ',i,xv(i),xv(i+1)
       endif


          X=(XV(I+1)-XV(I))/(PV(I+1)-PV(I))
          X=X*(P-PV(I+1))+XV(I+1)
          GO TO 1
       ENDIF
 10    CONTINUE
 1     CONTINUE

end subroutine polat


!####################################################################

subroutine write_diagnostics (jsg) 

!--------------------------------------------------------------------
!   this subroutine writes the donner_deep_mod history file
!--------------------------------------------------------------------

integer,                intent(in) :: jsg

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      jsg                global starting index of the domain segment
!                         being processed
!
!--------------------------------------------------------------------

      integer :: unit
     
!---------------------------------------------------------------------
!  open the history file and write the diagnostic and restart 
!  variables to it. in FMS, this will be a netCDF file.
!---------------------------------------------------------------------
      unit = open_file (file='history/donner_deep_diagnostics.his', &
!                       form='native', action='write')
                        form='unformatted', action='append')

!--------------------------------------------------------------------
!   these are the restart variables
!--------------------------------------------------------------------
      write (unit) cemetf    (:,jsg:jsg+jmaxp-jminp,:)
      write (unit) cememf    (:,jsg:jsg+jmaxp-jminp,:)
      write (unit) xcape_3d  (:,jsg:jsg+jmaxp-jminp)
      write (unit) qint_3d   (:,jsg:jsg+jmaxp-jminp)
      write (unit) omint_3d  (:,jsg:jsg+jmaxp-jminp)
      write (unit) tprea1_3d (:,jsg:jsg+jmaxp-jminp)

!--------------------------------------------------------------------
!   these are the 3D diagnostic variables
!--------------------------------------------------------------------
      write(unit) ceefc           
      write(unit) cecon          
      write(unit) cemfc        
      write(unit) cememf_mod
      write(unit) cual_3d   
      write(unit) fre_3d      
      write(unit) elt_3d     
      write(unit) cmus_3d       
      write(unit) ecds_3d      
      write(unit) eces_3d       
      write(unit) emds_3d                              
      write(unit) emes_3d 
      write(unit) qmes_3d      
      write(unit) wmps_3d          
      write(unit) wmms_3d        
      write(unit) tmes_3d    
      write(unit) dmeml_3d  
      write(unit) uceml_3d
      write(unit) umeml_3d
      write(unit) xice_3d
      write(unit) dgeice_3d
      write(unit) cuqi_3d
      write(unit) cuql_3d

!--------------------------------------------------------------------
!   these are the 2D diagnostic variables
!--------------------------------------------------------------------
      write(unit)  plcl_3d      
      write(unit)  plfc_3d       
      write(unit)  plzb_3d        
      write(unit)  coin_3d        
      write(unit)  dcape_3d      
      write(unit)  a1_3d            
      write(unit)  amax_3d      
      write(unit)  amos_3d            
      write(unit)  ampta1_3d 
      write(unit)  rcoa1_3d     

      call close_file (unit)

!------------------------------------------------------------------



end subroutine write_diagnostics


!####################################################################

subroutine donner_deep_end

!---------------------------------------------------------------------
!   this subroutine writes the donner_deep_mod restart file
!--------------------------------------------------------------------

       integer :: unit, next(2), dt(2), cal

!-------------------------------------------------------------------
!    open unit for restart file
!-------------------------------------------------------------------
       unit = open_file (file='RESTART/donner_deep.res', &
			 form='native', action='write')

!-------------------------------------------------------------------
!    file writing is currently single-threaded. write out restart
!    version, next time to call donner_deep_mod, donner_deep_timestep 
!    and calendar type being used.
!-------------------------------------------------------------------
       if (get_my_pe() == 0) then
	 write (unit) restart_versions(size(restart_versions))

	 call get_time (Next_donner_deep_time,next(1), next(2) )
	 call get_time (Donner_deep_timestep, dt(1), dt(2))
	 cal = get_calendar_type()
	 write (unit) next, dt, cal
       endif

!-------------------------------------------------------------------
!    write out the donner_deep restart variables.
!-------------------------------------------------------------------
      call write_data (unit, cemetf  )
      call write_data (unit, cememf  )
      call write_data (unit, xcape_3d)
      call write_data (unit, qint_3d )
      call write_data (unit, omint_3d)
      call write_data (unit, tprea1_3d)

!-------------------------------------------------------------------
!    close restart file unit.
!------------------------------------------------------------------
       call close_file (unit)

!---------------------------------------------------------------------

 
end subroutine donner_deep_end



!####################################################################

subroutine write_restart_donner_deep

!---------------------------------------------------------------------
!   this subroutine writes the donner_deep_mod restart file
!--------------------------------------------------------------------

    integer :: unit, next(2), dt(2), cal


    unit = open_file (file='RESTART/donner_deep.res', &
!                     form='native', action='write')
                      form='unformatted', action='write')
 
!   --- single threading of restart output ---
    if (get_my_pe() == 0) then
!      --- write last value in list of restart versions ---
      write (unit) restart_versions(size(restart_versions))

      call get_time (Next_donner_deep_time, next(1), next(2) )
      call get_time (Donner_deep_timestep, dt(1), dt(2) )
      cal = get_calendar_type()
      write (unit) next, dt, cal
    endif
 
!   --- data ---
!     call write_data (unit, cemetf  )
!     call write_data (unit, cememf  )
!     call write_data (unit, xcape_3d)
!     call write_data (unit, qint_3d )
!     call write_data (unit, omint_3d)
!     call write_data (unit, tprea1_3d)


   write (unit) cemetf  
   write (unit) cememf  
   write (unit) xcape_3d
   write (unit) qint_3d 
   write (unit) omint_3d
   write (unit) tprea1_3d
 
   call close_file (unit)
 
!---------------------------------------------------------------------

 
end subroutine write_restart_donner_deep


!##################################################################


 SUBROUTINE STRAT_CLOUD_DONNER_TEND (Dmeso, qlmeso, qimeso, &
         Mtot, phalf, ql, qi, cf, qltend, qitend, cftend)

!-----------------------------------------------------------------------
! input
!
!               vertical index 1 at model top
! Dmeso		mass detrainment rate from mesoscale region to large-scale
! 		region (sec-1)
! qlmeso	cloud liquid specific humidity (kg condensate/kg air)
! qimeso	cloud ice specific humidity (kg condensate/kg air)
! Mtot		total mass flux = mesoscale_mass_flux + convective_mass_flux
!               (kg /m2/sec) defined on level interfaces
!
!               NOTE: Regardless of what they contain, Mtot(:,:,1)
!                     Mtot(:,:,kdim+1) will be assumed to be zero.
!
! phalf		pressure on model interfaces (Pa)
! ql		large-scale cloud liquid specific humidity
! qi		large-scale cloud ice specific humidity
! cf		large-scale cloud fraction (0-1)
!
! output	
!
! qltend	large-scale cloud liquid tendency (kg cond/kg air/sec)
! qitend	large-scale cloud ice tendency (kg cond/kg air/sec)
! cftend	large-scale cloud fraction tendency (1/sec)
!
!-----------------------------------------------------------------------
   real, intent(in),  dimension(:,:,:) :: Dmeso, qlmeso, qimeso
   real, intent(in),  dimension(:,:,:) :: ql, qi, cf
   real, intent(in),  dimension(:,:,:) :: Mtot,phalf
   real, intent(out), dimension(:,:,:) :: qltend, qitend,cftend
!-----------------------------------------------------------------------
   integer kdim
!ljd
    integer kdima
    integer kdimb
!ljd
   real, dimension(size(Dmeso,1),size(Dmeso,2),size(Dmeso,3)) :: mass    ! kg
!  air/m2 of each level
   integer k
   integer kk
!-----------------------------------------------------------------------

   kdim = size(Dmeso,3)
!ljd
   kdima=size(Dmeso,1)
   kdimb=size(Dmeso,2)
!ljd
   mass(:,:,1:kdim) = (phalf(:,:,2:kdim+1)-phalf(:,:,1:kdim))/gravm

   qltend=0.
   qitend=0.
   cftend=0.
    
   !do incoming compensating subsidence fluxes from above
   qltend(:,:,2:kdim)=qltend(:,:,2:kdim) +  Mtot(:,:,2:kdim) * &
                      0.5*(ql(:,:,1:kdim-1)+ql(:,:,2:kdim))/mass(:,:,2:kdim)
   qitend(:,:,2:kdim)=qitend(:,:,2:kdim) +  Mtot(:,:,2:kdim) * &
                      0.5*(qi(:,:,1:kdim-1)+qi(:,:,2:kdim))/mass(:,:,2:kdim)
   cftend(:,:,2:kdim)=cftend(:,:,2:kdim) +  Mtot(:,:,2:kdim) * &
                      0.5*(cf(:,:,1:kdim-1)+cf(:,:,2:kdim))/mass(:,:,2:kdim)
!ljd
!  do ilon=1,kdima
! do jlat=1,kdimb
!		      do k=1,kdim
!  if (qimeso(ilon,jlat,k) .ge. 1.0e-09) then
!   kk=6
!  write(6,*) 'kk,qltend,qitend,cftend= ',kk,qltend(ilon,jlat,kk), &
!  qitend(ilon,jlat,kk), cftend(ilon,jlat,kk)
!  go to 11
!  end if
!  end do
!  end do
!  end do
! 11  continue
!ljd   
   !do outgoing compensating subsidence fluxes out the bottom
   qltend(:,:,1:kdim-1)=qltend(:,:,1:kdim-1) -  Mtot(:,:,2:kdim) * &
                      0.5*(ql(:,:,1:kdim-1)+ql(:,:,2:kdim))/mass(:,:,1:kdim-1)
   qitend(:,:,1:kdim-1)=qitend(:,:,1:kdim-1) -  Mtot(:,:,2:kdim) * &
                      0.5*(qi(:,:,1:kdim-1)+qi(:,:,2:kdim))/mass(:,:,1:kdim-1)
   cftend(:,:,1:kdim-1)=cftend(:,:,1:kdim-1) -  Mtot(:,:,2:kdim) * &
                      0.5*(cf(:,:,1:kdim-1)+cf(:,:,2:kdim))/mass(:,:,1:kdim-1)
!ljd
!  do ilon=1,kdima
! do jlat=1,kdimb
!  do k=1,kdim
!  if (qimeso(ilon,jlat,k) .ge. 1.0e-09) then
!   kk=6
!  write(6,*) 'kk,qltend,qitend,cftend= ',kk,qltend(ilon,jlat,kk), &
!  qitend(ilon,jlat,kk), cftend(ilon,jlat,kk)
!  go to 12
!  end if
!  end do
!  end do
!  end do
!  12 continue
!ljd   
   !do detrainment from meso region
   qltend(:,:,:) = qltend(:,:,:) + Dmeso(:,:,:)*qlmeso(:,:,:)
   qitend(:,:,:) = qitend(:,:,:) + Dmeso(:,:,:)*qimeso(:,:,:)
   where ((qlmeso+qimeso) .ge. 1.e-10)
   cftend(:,:,:) = cftend(:,:,:) + Dmeso(:,:,:)
   end where
!ljd
!  do ilon=1,kdima
! do jlat=1,kdimb
!  do k=1,kdim
!  if (qimeso(ilon,jlat,k) .ge. 1.0e-09) then
!   kk=6
!  write(6,*) 'kk,qltend,qitend,cftend= ',kk,qltend(ilon,jlat,kk), &
!  qitend(ilon,jlat,kk), cftend(ilon,jlat,kk)
!  go to 13
!  end if
!  end do
!  end do
!  end do
! 13 continue
!ljd
!ljd
!  do ilon=1,kdima
! do jlat=1,kdimb
! do k=1,kdim
!  if (qimeso(ilon,jlat,k) .ge. 1.0e-09) then
!  do kk=1,kdim
!  write(6,*) 'kk,dmeso,qimeso= ',kk,dmeso(ilon,jlat,kk),  &
!  qimeso(ilon,jlat,kk)
!  write(6,*) 'kk,phalf,mtot= ',kk,phalf(ilon,jlat,kk),    &
!  mtot(ilon,jlat,kk)
!  write(6,*) 'kk,ql,qi,cf= ',kk,ql(ilon,jlat,kk),qi(ilon,jlat,kk), &
!  cf(ilon,jlat,kk)
!  write(6,*) 'kk,qltend,qitend,cftend= ',qltend(ilon,jlat,kk), &
!  qitend(ilon,jlat,kk), cftend(ilon,jlat,kk)
!  end do
!  call error_mesg('strat_CLOUD_DONNER_TEND','diag stop',FATAL)
!  end if
!  end do
!  end do
!  end do
!ljd
      
!-----------------------------------------------------------------------

END SUBROUTINE STRAT_CLOUD_DONNER_TEND




	  end module donner_deep_mod


