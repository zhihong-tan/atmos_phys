                   module sea_esf_rad_mod

!-----------------------------------------------------------------------
!
!             interface module for sea_esf radiation-step physics
!
!-----------------------------------------------------------------------

   use utilities_mod,        only: get_my_pe, FATAL, NOTE, close_file,&
				   read_data, write_data, open_file, &
				   print_version_number,  &
				   check_nml_error, file_exist,  &
				   error_mesg, get_domain_decomp
   use rad_output_file_mod,  only: write_rad_output_file,    &
				   rad_output_file_end,&
			 	   rad_output_file_init, hold_tendencies
   use astronomy_mod,        only: set_ref_date_of_ae
   use store_rad_output_mod, only: store_radiation_output
   use diag_manager_mod,     only: register_diag_field, send_data
   use time_manager_mod,     only: time_type, increment_time,     &
				   get_date, &
                                   set_date, set_time, operator(<), &
				   operator(-), get_date_julian, &
				   set_date_julian, get_time, &
				   operator(+), operator(>=), &
				   operator(>),   &
				   operator(/=), get_calendar_type
   use strat_cloud_mod,      only: add_strat_tend
   use diag_integral_mod,    only: diag_integral_field_init, &
                                   sum_diag_integral_field
   use radiative_gases_mod,  only: radiative_gases_init,    &
				   radiative_gases_end
   use cool_to_space_mod,    only: cool_to_space_init
   use optical_path_mod,     only: optical_path_init
   use longwave_tables_mod,  only: longwave_tables_init
   use longwave_aerosol_mod, only: longwave_aerosol_init
   use gas_tf_mod,           only: gas_tf_init, gas_tf_end
   use co2_source_mod,       only: co2_source_init
   use longwave_setup_mod,   only: longwave_setup_init
   use longwave_clouds_mod,  only: longwave_clouds_init
   use radiation_diag_mod,   only: radiation_diag_init, &
			           radiag_from_astronomy, &
				   radiation_diag_end
   use longwave_driver_mod,  only: longwave_driver_init, &
				   get_lw_output
   use longwave_fluxes_mod,  only: longwave_fluxes_init
   use constants_new_mod,    only: define_constants,         &
				   pstd_mks, sigmasb, secday, &
				   frezdk, rh2oair
   use longwave_params_mod,  only: longwave_params_init
   use rad_utilities_mod,    only: rad_utilities_init, &
				   radiation_control_type, &
				   Rad_control, &
     			           define_environment, &
				   environment_type, Environment, &
				   cloudrad_control_type, &
                                   Cldrad_control
   use shortwave_driver_mod, only: shortwave_driver_init, &
				   get_sw_output
   use esfsw_parameters_mod, only: esfsw_parameters_init
   use lhsw_driver_mod,      only: lhsw_driver_init
   use radiation_package_mod,only: radpack_init
   use rad_step_setup_mod,   only: rad_step_setup_init, &
				   rad_step_setup_dr, &
                                   press, temp, rh2o, &
				   rad_step_setup_dealloc
   use radiation_step_mod,   only: radiation_step_init, &
				   radiation_step_alloc, &
				   radiation_step_dr, &
				   radiation_step_dealloc, &
				   radiation_step_time_vary
   use ozone_mod,            only: ozone_init
   use cloudrad_package_mod, only: cloudrad_package_init
   use astronomy_package_mod,only: astronomy_package_init, &
				   get_astronomy_for_swrad, &
				   return_cosz, get_solar_distance
   use surface_albedo_mod,   only: surface_albedo_init
   use std_pressures_mod,    only: set_std_pressures, skyp_std, &
				   sigp_std, lanz_std,   &
				   std_pressures_init
   use lw_gases_stdtf_mod,   only: lw_gases_stdtf_init

!---------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!             driver module for sea_esf radiation step physics
!
!--------------------------------------------------------------------




!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!        character(len=5), parameter  ::  version_number = 'v0.09'
         character(len=128)  :: version =  '$Id: sea_esf_rad.F90,v 1.2 2001/07/05 17:33:13 fms Exp $'
         character(len=128)  :: tag     =  '$Name: fez $'


!---------------------------------------------------------------------
!-------  interfaces --------

public     &
          sea_esf_rad_init, sea_esf_rad_time_vary, &
	  get_lrad, sea_esf_rad, sea_esf_rad_end


private    &
	  compute_average, return_average, radiation_netcdf,  &
	  update_rad_fields, define_rad_input, define_astro_input, &
	  diag_field_init


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=11)              :: model_name= ' ',   &
			          peripherals_source=' ', &
				  application_type = '  ', &
				  column_type = '    '
logical                        :: rsd = .false.

! Namelist variables:
!   rsd (repeat same day) : call radiation for the specified rad_date,
!                           running through the diurnal cycle
!*********************************************************************
! these variables are only used when running fms: 
!   dt_rad  =  radiative time step in seconds (integer)
!   offset  =  offset for radiative time step (in seconds) 
!              note if offset=1 radition will be done on the first
!              time step after t=0, dt_rad, 2*dt_rad, etc.
!              for now do not change this value
!   do_average = if true, astronomy, etc.forcing  is averaged over 
!                radiative timestep
!   rad_date  = fixed date (year,month,day,hour,min,sec)
!               for radiation (solar info, ozone, clouds)

integer                       :: dt_rad=14400
integer                       :: offset=1 
logical                       :: do_average=.false.    
integer, dimension(6)         :: rad_date = (/ 0, 0, 0, 0, 0, 0 /)
!*********************************************************************


namelist /sea_esf_rad_nml/ dt_rad, offset, do_average, &
                           rad_date, model_name, &
	        	   peripherals_source, &
			   application_type, &
			   column_type, rsd

!---------------------------------------------------------------------
!------- public data ------





!---------------------------------------------------------------------
!------- private data ------

!-------------------------------------------------------------------- 
!   list of restart versions readable by this module ----
!-------------------------------------------------------------------- 
integer, dimension(1) :: restart_versions = (/ 1 /)

!--------------------------------------------------------------------
!   tdt_rad = the current radiative (sw + lw) heating rate
!--------------------------------------------------------------------
!real,    dimension(:,:,:), allocatable   :: tdt_rad, tdt_radcf, &
real,    dimension(:,:,:), allocatable   :: tdt_rad, tdt_rad_clr, &
                                            ttt_l, rrr_l, press3d_l
real,    dimension(:,:),   allocatable   :: flux_sw_surf, flux_lw_surf,&
					    ts_l, fracday, cosz
integer, dimension(:),     allocatable   :: jabsp, iabsp
integer, dimension(:,:),   allocatable   :: nsum

real,    dimension(:,:,:), allocatable   :: psum,tsum,qsum
real,    dimension(:,:),   allocatable   :: asum,csum, usum, fsum, &
                                            albedo, snow, rrsum
real                                     :: rrsun
real,    dimension(:,:),   allocatable   :: coszen_save
real,    dimension(:,:),   allocatable   :: cosz_sv, fracday_sv

logical            ::  do_init=.true.
logical            ::  do_rad
integer            ::  num_pts, total_pts
type(time_type)    ::  Next_rad_time, Rad_time_step
integer            ::  xjd, kmin, kmax
logical            ::  use_rad_date

!---------------------------------------------------------------------
! diagnostics field informations 
!---------------------------------------------------------------------
!integer :: id_alb_sfc 
integer :: id_alb_sfc, id_cosz, id_fracday 

integer, dimension(2) :: id_tdt_sw,   id_tdt_lw,  &
                         id_swdn_toa, id_swup_toa, id_olr, &
                         id_swup_sfc, id_swdn_sfc,         &
                         id_lwup_sfc, id_lwdn_sfc

character(len=9), parameter :: mod_name = 'radiation'

real :: missing_value = -999.

integer            ::  x(4), y(4)
integer            ::  kmin, kmax





!---------------------------------------------------------------------
!---------------------------------------------------------------------






                         contains


subroutine sea_esf_rad_init (kd, axes, lonb, latb, pref,&
				  Timesince, th_in, zrlong_in, &
				  qlevel_in)

!---------------------------------------------------------------------
!  all input arguments but kd are optional. their presence depends on 
!  the  application running the radiation package
!  in SKYHI: qlevel_in, th_in, zrlong_in
!  in STANDALONE: th_in, zrlong_in
!  in FMS:  lat, lonb, latb, sfull, lon, pref, Timesince, axes
!--------------------------------------------------------------------
integer,         intent(in)                           :: kd
integer,         intent(in), dimension(4),   optional :: axes
real,            intent(in), dimension(:),   optional :: qlevel_in,  &
							 th_in,  &
							 lonb, latb
real,            intent(in), dimension(:,:), optional :: zrlong_in, &
                                                         pref
type(time_type), intent(in),                 optional :: Timesince
!-----------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables
!-------------------------------------------------------------------
    real, dimension(:) ,  allocatable   ::  pd, plm, pd8, plm8
    real, dimension(:,:), allocatable   ::  zrlong
    real, dimension(:) ,  allocatable   ::  th, qlevel, &
					    sfull, lon, lat

    integer            ::  unit,io,ierr, k, j, i, xid, next(2), dt(2), &
			   cal, vers
    logical            ::  do_avg_init, end
    logical            ::  skyhi_type, sigma_type, lanz_type
    character(len=11)  ::  type_form
    type(time_type)    ::  Old_time_step, New_rad_time
    character(len=4)   ::  chvers
    integer            ::  idf, jdf
    integer            ::  hemi
    logical            ::  extended

!--------------------------------------------------------------------
!----- read namelist -----

    if ( file_exist('input.nml')) then
      unit = open_file (file='input.nml', action='read')
      ierr=1; do while (ierr /= 0)
      read  (unit, nml=sea_esf_rad_nml, iostat=io, end=10)
      ierr = check_nml_error(io,'sea_esf_rad_nml')
      enddo
10    call close_file (unit)
    endif

!      ----- write namelist -----

    unit = open_file ('logfile.out', action='append')
    if (get_my_pe() == 0) then
      write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
      write (unit, nml=sea_esf_rad_nml)
    endif
    call close_file (unit)

!-------------------------------------------------------------------
!   define vertical index extents of model
!-------------------------------------------------------------------
    kmin = 1
    kmax = kd

!-------------------------------------------------------------------
! pass information concerning the environment in which this package
! is being run to rad_step_setup_mod. 
!-------------------------------------------------------------------
    call define_environment (model_name, peripherals_source,   &
	    		     application_type, column_type)

!----------------------------------------------------------------------
!      set time step flag (time_type) 
!----------------------------------------------------------------------
!   if (Environment%running_fms) then
    if (Environment%running_gcm .and. Environment%running_fms) then
      if (offset > 0) then
        Next_rad_time = increment_time (Timesince, offset, 0)
      else
        Next_rad_time = increment_time (Timesince, dt_rad, 0)
      endif
      Rad_time_step = set_time (dt_rad, 0)
      Old_time_step = set_time (dt_rad, 0)
    endif

!-------------------------------------------------------------------
! define constants which depend on the environment in which the
! radiation code is being run
!-------------------------------------------------------------------
    call define_constants 
   
!-----------------------------------------------------------------------
! retrieve model dimensions for use within radiation package
!-----------------------------------------------------------------------
    call get_domain_decomp (x, y)
    idf = x(4)-x(3)+1
    jdf = y(4)-y(3)+1

!-----------------------------------------------------------------------
! allocate space for and define arrays holding the model latitudes 
! and longitudes.
!-----------------------------------------------------------------------
    if (Environment%running_gcm) then
      allocate (sfull (1:kmax) )
      allocate (lon   (1:idf ) )
      allocate (lat   (1:jdf ) )
      sfull(1:kmax) = 0.01*pref(1:kmax,1)
      lon(1:idf) = 0.5*(lonb(1:idf) + lonb(2:idf+1))
      lat(1:jdf) = 0.5*(latb(1:jdf) + latb(2:jdf+1))
    endif
    allocate (zrlong (idf,jdf) )
    allocate (th     (jdf) )
    if (Environment%running_gcm) then
      if (Environment%running_skyhi) then
        zrlong(:,:) = zrlong_in(:,:) 
        th(:)       = th_in(:) 
      else if (Environment%running_fms) then
        do j=1,jdf
          zrlong(:,j) = lon(1:idf)
        end do
        th(:) = lat(1:jdf)
      endif
    else if (Environment%running_standalone) then
      zrlong(:,:) = zrlong_in(:,:) 
      th(:)       = th_in(:) 
    endif

!--------------------------------------------------------------------
!  initialize the various modules associated with the radiation code.
!  order does matter, so be wary of any order changes that are consid-
!  ered!!
!-----------------------------------------------------------------------
    call astronomy_package_init (th, zrlong) 
    call rad_utilities_init
    call longwave_params_init
    call longwave_setup_init
    call radiation_diag_init (th, zrlong)
    call std_pressures_init (kmin, kmax)

!-----------------------------------------------------------------------
! define and fill arrays holding the model sigma levels.
!-----------------------------------------------------------------------
    if (Environment%running_gcm) then
      allocate (qlevel (kmin:kmax) )
      if (Environment%running_skyhi) then
        qlevel(:)   = qlevel_in(:)

      else if (Environment%running_fms) then
        qlevel(:) = sfull(:)/ (0.01*pstd_mks)
      endif
    endif

!--------------------------------------------------------------------
!  define the pressure level profiles for use with tf calculations. 
!-----------------------------------------------------------------------
    if (present(pref)) then

!--------------------------------------------------------------------
!  if the pressure array pref is input, define the pressure arrays to be
!  stored in std_pressures_mod. make sure that the array pref is 
!  properly dimensioned.
!-----------------------------------------------------------------------
      if (size(pref,2) /= 2) then
        call error_mesg ('sea_esf_rad_init', &
                  'must have two reference pressure profiles.', FATAL)
      endif

      allocate (pd    (kmin:kmax+1) )
      allocate (plm   (kmin:kmax+1) )
      allocate (pd8   (kmin:kmax+1) )
      allocate (plm8  (kmin:kmax+1) )

      pd(:) = pref(:,1) *1.0E-02 ! convert to mb
      plm(kmin)=0.
      do k=kmin+1,kmax
        plm(k)=0.5*(pd(k-1)+pd(k))
      enddo
      plm(kmax+1)=pd(kmax+1)

      pd8(:) = pref(:,2)*1.0E-02  ! convert to mb
      plm8(1)=0.
      do k=kmin+1,kmax
        plm8(k)=0.5*(pd8(k-1)+pd8(k))
      enddo
      plm8(kmax+1)=pd8(kmax+1)
      call set_std_pressures (plm, pd, plm8, pd8)
      deallocate (pd)
      deallocate (pd8)
      deallocate (plm)
      deallocate (plm8)

!--------------------------------------------------------------------
!  if the pressure array is not input, call a routine to calculate the
!  model pressure levels, based on model sigma data specified within
!  std_pressures_mod.
!-----------------------------------------------------------------------
    else
      if (Environment%running_gcm .and. Environment%running_skyhi) then

!-----------------------------------------------------------------------
!  call skyp_std if running skyhi
!-----------------------------------------------------------------------
        call skyp_std

!--------------------------------------------------------------------
!   determine column type being used and compute appropriate vertical 
!   standard pressures.
!--------------------------------------------------------------------
      else if (Environment%running_standalone) then
        type_form = Environment%column_type
        skyhi_type =                                                &
          trim(type_form) == 'skyl40' .or.     &
          trim(type_form) == 'skyl80'
        sigma_type =                                                &
          trim(type_form) == 'egrpl9' .or.     &
          trim(type_form) == 'climl9' .or.     &
          trim(type_form) == 'ccml12' .or.     &
          trim(type_form) == 'climl14' .or.    &
          trim(type_form) == 'egrpl18' .or.    &
          trim(type_form) == 'climl30' .or.    &
          trim(type_form) == 'nmcl18'
        lanz_type =                                                 &
          trim(type_form) == 'lanl45'
  
        if (skyhi_type) then
          call skyp_std
        else if (sigma_type) then
          call sigp_std
        else if (lanz_type) then
          call lanz_std
	else
	  call error_mesg ( 'sea_esf_rad_init', &
	       'specified column type is not recognized', FATAL)
        endif

!--------------------------------------------------------------------
! error condition
!--------------------------------------------------------------------
      else
        call error_mesg ('sea_esf_rad_init',   &
	  'pref must be passed in when running fms', FATAL)
      endif
    endif

!-------------------------------------------------------------------
!  complete the initialization of the remaining radiation package
!  modules. remember, order matters !!
!-------------------------------------------------------------------
    call rad_step_setup_init (kmax) 
    call radiation_step_init
    call ozone_init (kmin, kmax)
    call radpack_init
    call radiative_gases_init 
    call longwave_aerosol_init (kmin, kmax)
    call optical_path_init
    call longwave_tables_init
!-------------------------------------------------------------------
! co2_source_init must follow rad_step_setup_init
!-------------------------------------------------------------------
    call co2_source_init (kmin, kmax)
    call gas_tf_init (kmin, kmax) 
    call lw_gases_stdtf_init (kmin, kmax)
    call longwave_clouds_init
    call cool_to_space_init
    call longwave_driver_init
    call longwave_fluxes_init
    call shortwave_driver_init (kmin, kmax)
    if (Environment%running_gcm) then
      if (Environment%running_fms) then
        call cloudrad_package_init (kmin, kmax,  &
				    th, axes=axes, Time=Timesince, &
				  qlyr=qlevel, lonb=lonb, latb=latb)
      else if (Environment%running_skyhi) then
        call cloudrad_package_init (kmin, kmax,   &
				    th, Time=Timesince, qlyr=qlevel)
      endif
    else if (Environment%running_standalone) then
      call cloudrad_package_init (kmin, kmax, th)
    endif
    if (Environment%running_gcm) then
      if (Environment%running_fms) then
        call rad_output_file_init (axes=axes, Time=Timesince)
!id_cosz = register_diag_field    &
!	                (mod_name,'cosz',axes(1:2),  &
!	                 Timesince, 'cosz',&
!                                'unknown',missing_value=missing_value)
!       id_fracday = register_diag_field    &
!			(mod_name,'fracday',axes(1:2),   &
!			 Timesince, 'fracday',&
!                                'unknown',missing_value=missing_value)
      else if (Environment%running_skyhi) then
        call rad_output_file_init 
      endif
    else if (Environment%running_standalone) then
      call rad_output_file_init 
    endif
    call surface_albedo_init 
    
!------------------------------------------------------------------
!  be sure that do_average is only true when running fms with fms
!  peripherals
!------------------------------------------------------------------
    if (do_average .and. (.not. Environment%running_fms) ) then
      call error_mesg ('sea_esf_rad_init',  &
           'must not use avgd radiation inputs except in  fms', FATAL)
    endif

!-------------------------------------------------------------------
!  deallocate arrays
!-------------------------------------------------------------------
    deallocate (zrlong)
    deallocate (th )
    if (allocated (qlevel)) then
      deallocate (qlevel)
    endif

!------------------------------------------------------------------
!  perform some initialization which is unique to fms
!------------------------------------------------------------------
    if ( Environment%running_gcm .and. Environment%running_fms) then

      allocate (tdt_rad     (idf,jdf,kmax))
      allocate (flux_sw_surf(idf,jdf))
      allocate (flux_lw_surf(idf,jdf))

      if (Rad_control%do_totcld_forcing) then
        allocate (tdt_rad_clr  (idf,jdf,kmax))
      endif

!-----------------------------------------------------------------------
!  initialize time averaged input data 
!-----------------------------------------------------------------------
      if (do_average) then
        allocate (psum(idf,jdf,kmax+1), tsum(idf,jdf,kmax),    &
                  qsum(idf,jdf,kmax), asum(idf,jdf), csum(idf,jdf),  &
                  fsum(idf,jdf), usum(idf,jdf), nsum(idf,jdf), &
		  rrsum(idf,jdf))
      endif
      allocate (coszen_save( idf, jdf) )

      do_avg_init = .false.

!-----------------------------------------------------------------------
!   read restart data                                             
!-----------------------------------------------------------------------
      if (file_exist('INPUT/sea_esf_rad.res')) then
        unit = open_file (file='INPUT/sea_esf_rad.res',  &
                          form='native', action='read')

!-----------------------------------------------------------------------
! read and check restart version number 
!-----------------------------------------------------------------------
        read (unit) vers
        if ( .not. any(vers == restart_versions) ) then
          write (chvers,'(i4)') vers
          call error_mesg ('sea_esf_rad_init', &
                    'restart version '//chvers//' cannot be read &
                     &by this module version', FATAL)
        endif

!-----------------------------------------------------------------------
!  reading alarm info 
!-----------------------------------------------------------------------
        read (unit) next,dt,cal

!-----------------------------------------------------------------------
!  override alarm with restart values 
!-----------------------------------------------------------------------
        Old_time_step = set_time (dt(1),dt(2))
        if ( cal == get_calendar_type() ) then
          Next_rad_time = set_time (next(1),next(2))
        else

!-----------------------------------------------------------------------
! calendar changed, will reset alarm 
!-----------------------------------------------------------------------
          call error_mesg ('sea_esf_rad_init', &
                   'current calendar not equal restart calendar', NOTE)
        endif

        call read_data (unit, tdt_rad)
        call read_data (unit, flux_sw_surf)
        call read_data (unit, flux_lw_surf)
        if (do_average) then
          call read_data (unit, nsum, end)
          if (end) go to 11
          call read_data (unit, psum)
          call read_data (unit, tsum)
          call read_data (unit, usum)
          call read_data (unit, fsum)
          call read_data (unit, qsum)
          call read_data (unit, asum)
          call read_data (unit, csum)
          call read_data (unit, rrsum)
          goto 12
11        do_avg_init = .true.
12        continue
        endif
        call close_file (unit)
      else
!-----------------------------------------------------------------------
!-------- initial surface flux set to 100 wm-2 ---------
!--------- only used for initial guess of sea ice temp ------
!-----------------------------------------------------------------------
        tdt_rad = 0.0
        flux_sw_surf = 50.0
        flux_lw_surf = 50.0
        if (do_average) do_avg_init = .true.
      endif

!----------------------------------------------------------------------
!     initialize input data averages 
!----------------------------------------------------------------------
      if (do_avg_init) then
        fsum=0; nsum=0; psum=0.; tsum=0.; usum=0.  
        qsum=0.; asum=0.; csum=0.; rrsum = 0.0
      endif

!----------------------------------------------------------------------
!  adjust radiation alarm if alarm interval has changed 
!----------------------------------------------------------------------
      if ( Rad_time_step /= Old_time_step ) then
        New_rad_time = Next_rad_time - Old_time_step + Rad_time_step
        if ( New_rad_time > Timesince ) then
          call error_mesg ('sea_esf_rad_init',  &
                     'radiation time step has changed, &
                     &next radiation time also changed', NOTE)
          Next_rad_time = New_rad_time
        endif
      endif

      total_pts=idf*jdf
      num_pts=0

!----------------------------------------------------------------------
!     check if optional radiative date should be used 
!----------------------------------------------------------------------
      if (rad_date(1) > 1900 .and.                        &
          rad_date(2) >   0  .and. rad_date(2) < 13 .and. &
          rad_date(3) >   0  .and. rad_date(3) < 32 ) then
        use_rad_date = .true.
      else
        use_rad_date = .false.
      endif


!----------------------------------------------------------------------
!   initialize quantities for integral package. hemi is the global
!   normalization factor used for hemispheric integrals. 
!----------------------------------------------------------------------
      hemi = 2
      extended = .true.
      call diag_integral_field_init ('olr',    'f8.3')
      call diag_integral_field_init ('abs_sw', 'f8.3')
      call diag_integral_field_init ('sntop_tot_sh',    'f16.11', extended=extended, hemi=hemi)
      call diag_integral_field_init ('lwtop_tot_sh',    'f16.11', extended=extended, hemi=hemi)
      call diag_integral_field_init ('sngrd_tot_sh',    'f16.11', extended=extended, hemi=hemi)
      call diag_integral_field_init ('lwgrd_tot_sh',    'f16.11', extended=extended, hemi=hemi)
      call diag_integral_field_init ('sntop_tot_nh',    'f16.11', extended=extended, hemi=hemi)
      call diag_integral_field_init ('lwtop_tot_nh',    'f16.11', extended=extended, hemi=hemi)
      call diag_integral_field_init ('sngrd_tot_nh',    'f16.11', extended=extended, hemi=hemi)
      call diag_integral_field_init ('lwgrd_tot_nh',    'f16.11', extended=extended, hemi=hemi)
      call diag_integral_field_init ('sntop_tot_gl',    'f16.11', extended=extended)
      call diag_integral_field_init ('lwtop_tot_gl',    'f16.11', extended=extended)
      call diag_integral_field_init ('sngrd_tot_gl',    'f16.11', extended=extended)
      call diag_integral_field_init ('lwgrd_tot_gl',    'f16.11', extended=extended)
      if (Rad_control%do_totcld_forcing) then
        call diag_integral_field_init ('sntop_clr_sh',    'f16.11', extended=extended, hemi=hemi)
        call diag_integral_field_init ('lwtop_clr_sh',    'f16.11', extended=extended, hemi=hemi)
        call diag_integral_field_init ('sngrd_clr_sh',    'f16.11', extended=extended, hemi=hemi)
        call diag_integral_field_init ('lwgrd_clr_sh',    'f16.11', extended=extended, hemi=hemi)
        call diag_integral_field_init ('sntop_clr_nh',    'f16.11', extended=extended, hemi=hemi)
        call diag_integral_field_init ('lwtop_clr_nh',    'f16.11', extended=extended, hemi=hemi)
        call diag_integral_field_init ('sngrd_clr_nh',    'f16.11', extended=extended, hemi=hemi)
        call diag_integral_field_init ('lwgrd_clr_nh',    'f16.11', extended=extended, hemi=hemi)
        call diag_integral_field_init ('sntop_clr_gl',    'f16.11', extended=extended)
        call diag_integral_field_init ('lwtop_clr_gl',    'f16.11', extended=extended)
        call diag_integral_field_init ('sngrd_clr_gl',    'f16.11', extended=extended)
        call diag_integral_field_init ('lwgrd_clr_gl',    'f16.11', extended=extended)
      endif

      call diag_field_init (timesince, axes)

    endif  ! (if running_fms)



!---------------------------------------------------------------------
!   set flag to indicate that module has been initialized.
!---------------------------------------------------------------------
    do_init=.false.

!-------------------------------------------------------------------


end subroutine sea_esf_rad_init




!#######################################################################


subroutine sea_esf_rad_time_vary (Time) 

type(time_type), intent(in), optional ::  Time
! Time will be present when running in FMS

!---------------------------------------------------------------------
! called from bgrid_physics_mod in fms and from rad_physics_driver.F 
! in SKYHI
!---------------------------------------------------------------------


   if (present (Time)) then
     call radiation_step_time_vary (Time)
   else                                      
     call radiation_step_time_vary
   endif


end subroutine sea_esf_rad_time_vary 




!####################################################################

 subroutine get_lrad (Time, do_rad, do_average_out, dt_rad_out, &
		       rad_time_out)

!------------------------------------------------------------------
type(time_type), intent(in)   :: Time
type(time_type), intent(out)  :: rad_time_out
logical,         intent(out)  :: do_rad, do_average_out
integer,         intent(out)  :: dt_rad_out
!------------------------------------------------------------------

      integer                 :: dum, tod(3)

!-------------------------------------------------------------------
! called from bgrid_physics_mod and physics_driver_mod, needed  only  
! in fms for sea_esf_rad_physics case.
!-------------------------------------------------------------------
      if (Time < Next_rad_time) then
        do_rad = .false.
      else
        do_rad = .true.
      endif

      if (rsd) then
        if ( .not. use_rad_date) call error_mesg ( 'get_lrad', &
             'if (rsd), must set rad_date(1:3) to valid date', FATAL)
        call get_date (Time, dum, dum, dum, tod(1), tod(2), &
	               tod(3))
        rad_time_out = set_date (rad_date(1), rad_date(2), rad_date(3),&
	                         tod(1), tod(2), tod(3))
      else if (use_rad_date) then
        rad_time_out = set_date (rad_date(1), rad_date(2), rad_date(3),&
                                 rad_date(4), rad_date(5), rad_date(6))
      else
        if (do_rad) then
          rad_time_out = Time
        else
          rad_time_out = Next_rad_time
        endif
      endif

      do_average_out = do_average
      dt_rad_out     = dt_rad


end subroutine get_lrad



!#######################################################################

subroutine sea_esf_rad (is, ie, js, je, pfull, phalf, t, q, ts, &
	        Time, Time_next, lat, lon, land, albedo_in, tdt, &
		        flux_sw, flux_lw, mask, kbot, coszen, &
		        snow_l_in, sealp_l_in, &
			cloud_ice_model, cloud_water_model)

!--------------------------------------------------------------------
!   the following variables are always present as arguments, regardless
!   of application:
!         is,ie,js,je  starting/ending i,j in global storage arrays;
!         pfull        pressure at full levels
!         phalf        pressure at half levels
!         t            temperature in model layers (full layers)
!         q            specific humidity of water vapor in model layers
!         ts           surface temperature in deg k
!--------------------------------------------------------------------
integer, intent(in)                         :: is, ie, js, je
real,    intent(in), dimension(:,:,:)       :: pfull, phalf, t, q
real,    intent(in), dimension(:,:)         :: ts


!-----------------------------------------------------------------------
!    the following variables are present only when sea_esf_rad is 
!    run within fms:
!         Time       time for determining whether radiation should be
!                    done 
!         Time_next  time used for diagnostics
!         lat        latitude,in radians
!         lon        longitude,in radians
!         land       fraction of surface which covered by land (real)
!         albedo_in  surface albedo (current only output)
!         tdt        temperature tendency
!         flux_sw    net shortwave surface flux (down-up) (w/m2)
!         flux_lw    downward longwave surface flux (w/m2)
!         mask
!         kbot
!
!         dt    NOTE:     time step in seconds, interval that 
!                    sea_esf_rad is called (not the interval
!                    that radiation calculations are done)
!
!   the following variables are present only when the standalone
!   code is run with predicted-cloud microphysics:
!
!         cloud_ice_model
!         cloud_water_model
!
!----------------------------------------------------------------------
type(time_type), intent(in),               optional :: Time, Time_next
real,            intent(in),    dimension(:,:),   optional ::    &
						     land, albedo_in, &
					             lat,  lon
real,            intent(inout), dimension(:,:,:), optional :: tdt
real,            intent(out),   dimension(:,:),   optional :: flux_sw, &
							      flux_lw
real,            intent(in),    dimension(:,:,:), optional :: mask, &
					            cloud_ice_model, &
                                                    cloud_water_model
integer,         intent(in),    dimension(:,:),   optional :: kbot
real,            intent(out),   dimension(:,:),   optional :: coszen
 

!-----------------------------------------------------------------------
!    the following variables are present only when sea_esf_rad is 
!    run within skyhi:
!         snow_l_in
!         sealp_l_in
!-----------------------------------------------------------------------
real, intent(in), dimension(:,:),optional     :: sealp_l_in, snow_l_in
  

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
     real, dimension(:,:,:), allocatable   :: cosangsolar
     real, dimension(:,:),   allocatable   :: sealp_l, snow_l
     type(time_type)   :: Rad_time
     logical           :: safe_to_dealloc
     integer           :: kb, k, i, j, nsolwg
     integer           :: dt, day, sec
     logical           :: used
     integer           :: dum, tod(3)

!-------------------------------------------------------------------
!   be sure this module has been initialized.
!-------------------------------------------------------------------
     if (do_init) call error_mesg ('sea_esf_rad',  &
                    'sea_esf_rad_init must first be called.', &  
							       FATAL)

!-------------------------------------------------------------------
!   compute the model timestep.
!-------------------------------------------------------------------
	if (present (Time) ) then
          call get_time ( Time_next-Time, sec, day)
	  dt = day*86400 + sec

	  if (dt <= 0) call error_mesg ('sea_esf_rad_mod',  &
		    'Time_next <= Time', FATAL)
        endif

!--------------------------------------------------------------------
!   determine if this is a radiation time step. set the logical control
!   appriopriately, and increment the counter used to determine when all
!   grid points have been updated. set do_rad to true for SKYHI case, 
!   since routine would not be called if it were not a radiation time 
!   step, and also to true for the standalone case.
!--------------------------------------------------------------------
     if (Environment%running_gcm) then
       if (Environment%running_skyhi) then
         do_rad = .true.
       else if (Environment%running_fms) then
         if (Time < Next_rad_time ) then
           do_rad = .false.
         else
           do_rad = .true.
           num_pts = num_pts + size(t,1) * size(t,2)
           if (num_pts == total_pts) then
             Next_rad_time = Next_rad_time + Rad_time_step
             num_pts = 0
           endif
         endif
       endif
     else if (Environment%running_standalone) then
       do_rad = .true.
     endif

!---------------------------------------------------------------------
!   note: this routine will be called on every time step in fms. 
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   if radiation is not to be calculated on this step and the option
!   to time-average the radiative forcing fields has not been activated,
!   all that is necessary is that one calls update_rad_fields to apply 
!   the existing radiative forcings and fluxes to those parts of the 
!   model where needed.
!--------------------------------------------------------------------
     if (.not. do_rad ) then
!  needed to provide coszen for physics_driver in bombay
         call define_astro_input (is, ie, js, je)
 	 coszen(:,:) = cosz(:,:)
 	 coszen_save(is:ie,js:je) = cosz(:,:)
!  this if included sends the cosz data every ts and eliminates the
!  wave pattern. (along with setting rad_time = time in get_lrad)
!   if ( id_cosz > 0 ) then
!            allocate (cosz_sv(is:ie,js:je))
!            cosz_sv(is:ie,js:je) = cosz(:,:)
!            used = send_data ( id_cosz, cosz_sv, Time_next, is, js )
!            deallocate (cosz_sv)
!          endif
!          if ( id_fracday > 0 ) then
!            allocate (fracday_sv(is:ie,js:je))
!            fracday_sv(is:ie,js:je) = fracday(:,:)
!            used = send_data (id_fracday, fracday_sv,    &
!		       Time_next, is, js )
!            deallocate (fracday_sv)
!          endif
! if (get_my_pe() == 0) then
!   do j=js,je
!   print *, js,js,js,js
!   print *, (coszen_save(i,j), i=1,100)
!   end do
! endif
	 deallocate (fracday)
	 deallocate (cosz)
       if (.not. do_average) then
         call update_rad_fields (is, ie, js, je, tdt, flux_sw, flux_lw)

!-----------------------------------------------------------------------
!   if time-averaged fields are desired as input to the next radiation
!   calculation, retrieve and save the fields on this step. this option
!   only available within fms.
!-----------------------------------------------------------------------
       else if (do_average) then

!---------------------------------------------------------------------
!    define_rad_input will produce the t, q, p fields 
!---------------------------------------------------------------------
	 call define_rad_input (is, ie, js, je, t, q, pfull, phalf, &
				ts)

!---------------------------------------------------------------------
!    define_astro_input will produce the zenith angles and fracday for
!    input to the averaging routine
!---------------------------------------------------------------------
          call define_astro_input (is, ie, js, je)

!--------------------------------------------------------------------
!  pass the needed variables to compute_average to be held until the
!  next radiation step.
!-------------------------------------------------------------------
         call compute_average (is, ie, js, je, press3d_l, ttt_l,  &
    	                       ts_l, rrr_l, albedo_in, cosz,  &
			       fracday, rrsun)
!-------------------------------------------------------------------
!  deallocate the arrays that have been allocated and call 
!  update_rad_fields to apply the radiative forcings and fluxes as 
!  needed.
!-------------------------------------------------------------------
! coszen(:,:) = cosz(:,:)
!???? coszen_save(is:ie,js:je) = cosz(:,:)
         deallocate (cosz     )
         deallocate (fracday  )
         deallocate (ttt_l   )
         deallocate (rrr_l   )
         deallocate (ts_l    )
         deallocate (press3d_l)
         call update_rad_fields (is, ie, js, je, tdt, flux_sw,  &
 	          	         flux_lw)
       endif  ! (do_average)
!!2     coszen(:,:) = coszen_save(is:ie, js:je)

!-----------------------------------------------------------------
!    if this is a radiation step, set up to compute radiation
!-----------------------------------------------------------------
     else if  (do_rad) then

!-----------------------------------------------------------------
!   define dimensions and mapping arrays from global to local coord-
!   inates using the input arguments provided.
!-----------------------------------------------------------------
       allocate (iabsp(1:ie-is+1))
       allocate (jabsp(1:je-js+1))
       do j=1,je-js+1
         jabsp(j) = js+j-1
       end do
       do i=1,ie-is+1
         iabsp(i) = is+i-1
       end do

!-------------------------------------------------------------------
!   define the date to be used in defining the astronomical parameters 
!   for the shortwave radiation calculation. in fms the option to
!   specify a radiation date different from the current model time is
!   available.
!-------------------------------------------------------------------
       if (Environment%running_gcm) then
         if (Environment%running_fms) then
           if (rsd) then
             if ( .not. use_rad_date) call error_mesg ( 'get_lrad', &
                'if (rsd), must set rad_date(1:3) to valid date', FATAL)
             call get_date (Time, dum, dum, dum, tod(1), tod(2), tod(3))
             Rad_time = set_date (rad_date(1), rad_date(2),& 
                                  rad_date(3), tod(1), tod(2), &
                                  tod(3))
           else if (use_rad_date) then
             Rad_time = set_date (rad_date(1), rad_date(2),    &
                                  rad_date(3), rad_date(4),    &
                                  rad_date(5), rad_date(6))
           else
             Rad_time = time
           endif
         endif
       endif

!---------------------------------------------------------------------
!    define_rad_input will produce the t, q, p fields 
!---------------------------------------------------------------------
       call define_rad_input (is, ie, js, je, t, q, pfull, phalf, &
	          	      ts)

!-------------------------------------------------------------------
!  allocate space for and define the surface type (land, sea, seaice 
!  or glacier) and snow cover. these are needed only for skyhi and the 
!  skyhi-mimic version of fms.
!------------------------------------------------------------------
       if (Environment%running_gcm) then
         if ( Environment%running_skyhi  ) then
           allocate (sealp_l   (is:ie, js:je) )
           allocate (snow_l    (is:ie, js:je) )
         endif 
!------------------------------------------------------------------
!  in skyhi, seal and snow cover are input arguments.
!------------------------------------------------------------------
         if (Environment%running_skyhi) then
           sealp_l(:,:) = sealp_l_in(:,:)
           snow_l (:,:) = snow_l_in(:,:)
         endif
       endif

!-----------------------------------------------------------------------
!   obtain variables needed when executing within fms. currently these
!   values are passed to rad_step_setup_mod as arguments in fms; in 
!   skyhi they are obtained directly from the appropriate modules 
!   within the radiation package when needed.
!-----------------------------------------------------------------------
       if (Environment%running_gcm) then
         if (Environment%running_fms) then

!-------------------------------------------------------------------
!  allocate and define the surface albedo array to be passed into the
!  rad_step_setup_mod.
!-----------------------------------------------------------------------
           allocate (albedo (is:ie, js:je)  )
           albedo(:,:) = albedo_in(:,:)

!---------------------------------------------------------------------
!    if time-averaged values are to be used as input to the radiation
!    package, call define_astro_input to produce the zenith angles,
!    fracday and earth-sun distance values.
!---------------------------------------------------------------------
!  call this to get cosz for return to physics_driver for bombay
           call define_astro_input (is, ie, js, je)
	   if ( id_cosz > 0 ) then
             allocate (cosz_sv(is:ie,js:je))
             cosz_sv(is:ie,js:je) = cosz(:,:)
             used = send_data ( id_cosz, cosz_sv, Time_next, is, js )
             deallocate (cosz_sv)
           endif
           if ( id_fracday > 0 ) then
             allocate (fracday_sv(is:ie,js:je))
             fracday_sv(is:ie,js:je) = fracday(:,:)
             used = send_data (id_fracday, fracday_sv,    &
			       Time_next, is, js )
             deallocate (fracday_sv)
           endif
           if (do_average) then

!--------------------------------------------------------------------
!  call compute_average to pass in this time step's values to be used
!  in creating the time-averaged values.
!---------------------------------------------------------------------
             call compute_average (is, ie, js, je, press3d_l, ttt_l,  &
	    	                   ts_l, rrr_l, albedo_in, cosz,  &
			           fracday, rrsun)

!-------------------------------------------------------------------
!  bring back the time-averaged values of the input quantities via a 
!  call to return_average.  send these values of zenith angle, daylight 
!  fraction and earth-sun distance to radiation_diag_mod and to 
!  astronomy_package_mod.
!-------------------------------------------------------------------
             call return_average (is,ie,js,je, press3d_l, ttt_l,   &
	          		  ts_l, rrr_l, albedo, cosz, fracday, &
				  rrsun)
             where (fracday > 1.00) fracday =  1.00
             call return_cosz (cosz, fracday, rrsun)
             nsolwg = 1
             allocate (cosangsolar(is:ie, 1:je-js+1, 1))
             cosangsolar(:,:,1) = cosz(:,:)
	     coszen_save(is:ie,js:je) = cosz(:,:)
	     coszen(:,:) = cosz(:,:)
             do j=1,je-js+1      
               if (Rad_control%do_raddg(js+j-1)) then
                 call radiag_from_astronomy (cosangsolar, fracday,  &
		                        nsolwg, jabsp(j), j, 1, ie-is+1)
               endif
             end do
             deallocate (cosangsolar)
!            deallocate (cosz)
!             deallocate (fracday)
	   else
	     coszen_save(is:ie,js:je) = cosz(:,:)
	     coszen(:,:) = cosz(:,:)
! if (get_my_pe() == 0) then
!   do j=js,je
!   print *, js,js,js,js
!   print *, (coszen_save(i,j), i=1,100)
!   end do
! endif
           endif       ! (end of do_average loop)
             deallocate (cosz)
             deallocate (fracday)
         endif     ! (end of running_fms loop)
       endif   ! (end of running_gcm loop)

!--------------------------------------------------------------------
!   call the setup routines for the radiation package. 
!--------------------------------------------------------------------
       if (Environment%running_gcm) then
         if (Environment%running_fms) then
           call rad_step_setup_dr (is, ie, js, je, &
				   ttt_l, rrr_l, ts_l, jabsp,   &
    		                   iabsp, press3d_l, Rad_time=Rad_time,&
				   lat=lat, albedo=albedo, land_in=land)
         else if (Environment%running_skyhi) then
           call rad_step_setup_dr (is, ie, js, je, &
				   ttt_l, rrr_l, ts_l, jabsp,   &
	  		           iabsp, press3d_l, snow_l=snow_l,  &
				   sealp_l=sealp_l)
         endif
       else if (Environment%running_standalone) then
	 if (Cldrad_control%do_pred_cld_microphys) then
           call rad_step_setup_dr (is, ie, js, je, &
                                   ttt_l, rrr_l, ts_l, jabsp,   &
                                   iabsp, press3d_l,            &
                                   land_in=land,                &
                                   cloud_ice_in=cloud_ice_model, &
                                   cloud_water_in=cloud_water_model )
         else
         call rad_step_setup_dr (is, ie, js, je, &
				 ttt_l, rrr_l, ts_l, jabsp,   &
  	 	                 iabsp, press3d_l)
	 endif
       endif

!----------------------------------------------------------------------
!   deallocate arrays which are no longer needed
!----------------------------------------------------------------------
       deallocate (ttt_l   )
       deallocate (rrr_l   )
       deallocate (ts_l    )
       deallocate (press3d_l)
       if (allocated (sealp_l) )then
         deallocate (sealp_l)
         deallocate (snow_l)
       endif
       deallocate (iabsp)
       deallocate (jabsp)

!---------------------------------------------------------------------
!  call routine to allocate major module variables
!---------------------------------------------------------------------
       call radiation_step_alloc
 
!---------------------------------------------------------------------
!  call routine which performs radiation calculation.
!---------------------------------------------------------------------
       if (Environment%running_fms) then
         if (present (kbot)) then
           call radiation_step_dr (is, ie, js, je,   &
			   Time_next=Time_next, kbot=kbot, mask=mask)
         else
           call radiation_step_dr (is, ie, js, je,   &
				   Time_next=Time_next)
         endif
       else 
           call radiation_step_dr (is, ie, js, je)
       endif

!-------------------------------------------------------------------
!   handle the output from the radiation package. store data where 
!   desired and calculate diagnostic quantities, both as fields and as
!   integrals.
!-------------------------------------------------------------------
       if (Environment%running_gcm) then

!-------------------------------------------------------------------
!   call the output interface routine for the radiation package output.
!   for skyhi, access store_radiation_output_mod to produce skyhi-like
!   output. in fms call radiation_output to  produce integrals and a
!   netcdf output file.
!-------------------------------------------------------------------
         if ( Environment%running_skyhi) then
           call store_radiation_output
         else if (Environment%running_fms) then 
           call radiation_netcdf (is, ie, js, je, ts, Time_next)
         endif
       endif
         
!-------------------------------------------------------------------
! deallocate arrays from radiation calculation
!-------------------------------------------------------------------
       call radiation_step_dealloc
       call rad_step_setup_dealloc

       if (Environment%running_gcm) then
         if (Environment%running_fms) then
           deallocate (albedo  )
         endif
       endif

!-------------------------------------------------------------------
!   call routine to produce new values of boundary fluxes and heating
!   rate and to write output file, if desired.
!-------------------------------------------------------------------
       if (Environment%running_gcm) then
         if (Environment%running_fms) then
           call update_rad_fields (is, ie, js, je, tdt, flux_sw,  &
				   flux_lw)
	   safe_to_dealloc = .true.
           call write_rad_output_file (safe_to_dealloc, Time, is, js)
         else if (Environment%running_skyhi) then
           safe_to_dealloc = .false.
           call write_rad_output_file (safe_to_dealloc)
         endif
       else if (Environment%running_standalone) then
         safe_to_dealloc = .true.
        call write_rad_output_file (safe_to_dealloc)
       endif
     endif  ! (not do_rad, not do_average)
!-------------------------------------------------------------------


end subroutine sea_esf_rad




!#######################################################################


subroutine sea_esf_rad_end

    integer :: unit, next(2), dt(2), cal

    if (Environment%running_gcm) then
      if (Environment%running_fms) then
    unit = open_file (file='RESTART/sea_esf_rad.res', &
                      form='native', action='write')

!   --- single threading of restart output ---
    if (get_my_pe() == 0) then
!      --- write last value in list of restart versions ---
       write (unit) restart_versions(size(restart_versions))

!      --- write alarm info ---
       call get_time (Next_rad_time,next(1),next(2))
       call get_time (Rad_time_step, dt (1), dt (2))
       cal = get_calendar_type()
       write (unit) next,dt,cal
    endif

!   --- data ---
    call write_data (unit, tdt_rad)
    call write_data (unit, flux_sw_surf)
    call write_data (unit, flux_lw_surf)

!   --- optional time avg data ---
  if (do_average) then
    call write_data (unit, nsum)
    call write_data (unit, psum)
    call write_data (unit, tsum)
    call write_data (unit, usum)
    call write_data (unit, fsum)
    call write_data (unit, qsum)
    call write_data (unit, asum)
    call write_data (unit, csum)
    call write_data (unit, rrsum)
  endif
    call close_file (unit)

      endif
    endif

    call radiative_gases_end
    call gas_tf_end
    call rad_output_file_end
    call radiation_diag_end

end subroutine sea_esf_rad_end



!#######################################################################

subroutine compute_average (is, ie, js, je, press, temp, tsur, rh2o, &
			    albedo, cosz, fracday, rrsun)

!-----------------------------------------------------------------------
real,    intent(in)                     :: rrsun
integer, intent(in)                     :: is,ie,js,je
real,    intent(in), dimension(:,:,:)   :: press,temp,rh2o
real,    intent(in), dimension(:,:)     :: albedo,cosz, fracday, tsur
!-----------------------------------------------------------------------
      integer  :: k, j, i

      nsum(is:ie,js:je)   = nsum(is:ie,js:je)   + 1
      psum(is:ie,js:je,:) = psum(is:ie,js:je,:) + press(:,:,:  )
      tsum(is:ie,js:je,:) = tsum(is:ie,js:je,:) + temp (:,:,:  )
      qsum(is:ie,js:je,:) = qsum(is:ie,js:je,:) + rh2o (:,:,:  )
      usum(is:ie,js:je  ) = usum(is:ie,js:je  ) + tsur (:,:  )
      asum(is:ie,js:je)   = asum(is:ie,js:je)   + albedo(:,:)
      csum(is:ie,js:je)   = csum(is:ie,js:je)   + cosz  (:,:)
      fsum(is:ie,js:je)   = fsum(is:ie,js:je)   + fracday(:,:)
      rrsum(is:ie,js:je)   = rrsum(is:ie,js:je)   + rrsun        

!-----------------------------------------------------------------------

end subroutine compute_average

!#######################################################################

subroutine return_average (is, ie, js, je, press, temp, tsur, rh2o, &
			   albedo, cosz, fracday, rrsun)

!-----------------------------------------------------------------------
integer, intent(in)                    :: is,ie,js,je
real,    intent(out)                   :: rrsun
real,    intent(out), dimension(:,:,:) :: press,temp,rh2o
real,    intent(out), dimension(:,:)   :: albedo,cosz, tsur, fracday
!-----------------------------------------------------------------------

      real, dimension(size(press,1),size(press,2)) :: dfsum
      integer  n, k

!--------- check for zero divid ---------
      n = count (nsum(is:ie,js:je) <= 0)
      if ( n > 0 ) then
          call error_mesg ('return_average in sea_esf_rad_mod',  &
                           'dividing average by zero.',  FATAL)
      endif

!--------------- compute averages --------------------------------------

      dfsum(:,:) = 1.0 / float(nsum(is:ie,js:je))

      do k=1,size(press,3)
         press(:,:,k) = psum(is:ie,js:je,k) * dfsum(:,:)
      enddo

      do k=1,size(rh2o,3)
         temp (:,:,k) = tsum(is:ie,js:je,k) * dfsum(:,:)
         rh2o (:,:,k) = qsum(is:ie,js:je,k) * dfsum(:,:)
      enddo

      albedo(:,:) = asum(is:ie,js:je) * dfsum(:,:)
      cosz  (:,:) = csum(is:ie,js:je) * dfsum(:,:)
      tsur  (:,:) = usum(is:ie,js:je) * dfsum(:,:)
      fracday(:,:) = fsum(is:ie,js:je) * dfsum(:,:)
      rrsun        = rrsum(is,js) * dfsum(1,1)


!  ----- zero out sums -----
      nsum(is:ie,js:je)   = 0
      psum(is:ie,js:je,:) = 0.0
      tsum(is:ie,js:je,:) = 0.0
      usum(is:ie,js:je  ) = 0.0
      fsum(is:ie,js:je  ) = 0.0
      qsum(is:ie,js:je,:) = 0.0
      asum(is:ie,js:je)   = 0.0
      csum(is:ie,js:je)   = 0.0
      rrsum(is:ie,js:je)   = 0.0

!-----------------------------------------------------------------------

end subroutine return_average

!#######################################################################

subroutine radiation_netcdf(is, ie, js, je, ts, Time_diag, mask) 

!--------------------------------------------------------------------
integer,               intent(in)           :: is, ie, js, je
real,dimension(:,:),   intent(in)           :: ts
type(time_type),       intent(in)           :: Time_diag
real,dimension(:,:,:), intent(in), optional :: mask
!--------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

     real, dimension(:,:,:), allocatable   ::   &
                                          tdtsw, tdtlw, &
                                          tdt_rad_loc,&
                                          tdt_rad_clr_loc, tdtsw_clr, &
                                          tdtlw_clr, flxnet, heatra,&
                                          flxnetcf, heatracf, &
                                          dfsw, ufsw, fsw, hsw, &
                                          dfswcf, ufswcf, fswcf, &
                                          hswcf
     real, dimension(:,:),   allocatable   :: swin, &
                                              swout, olr, swups,   &
                                              swdns, lwups,lwdns, &
                                              swin_clr, swout_clr,   &
                                              olr_clr, swups_clr, &
                                              swdns_clr, lwups_clr, &
                                              lwdns_clr, asfc
					      
     integer           :: k, idiag, ifld, i, j
     integer           :: ipass
     logical           :: used

!--------------------------------------------------------------------

!-------------------------------------------------------------------
!  retrieve lw heating rates and fluxes
!-------------------------------------------------------------------

     allocate ( flxnet (is:ie, js:je, KMIN:kmax+1) )
     allocate (heatra ( is:ie, js:je, KMIN:kmax) )

     if (Rad_control%do_totcld_forcing) then
       allocate (flxnetcf ( is:ie, js:je, KMIN:kmax+1) )
       allocate (heatracf ( is:ie, js:je, KMIN:kmax  ) )
       call get_lw_output (heatra, flxnet, heatracf, flxnetcf)
     else
       call get_lw_output (heatra, flxnet)
     endif

!-------------------------------------------------------------------
!  retrieve sw heating rates and fluxes
!-------------------------------------------------------------------
     allocate ( dfsw   (is:ie, js:je, KMIN:kmax+1) )
     allocate ( ufsw   (is:ie, js:je, KMIN:kmax+1) )
     allocate (  fsw   (is:ie, js:je, KMIN:kmax+1) )
     allocate (hsw    ( is:ie, js:je, KMIN:kmax  ) )

     if (Rad_control%do_totcld_forcing) then
       allocate (dfswcf   ( is:ie, js:je, KMIN:kmax+1) )
       allocate (ufswcf   ( is:ie, js:je, KMIN:kmax+1) )
       allocate ( fswcf   ( is:ie, js:je, KMIN:kmax+1) )
       allocate ( hswcf   ( is:ie, js:je, KMIN:kmax  ) )
       call get_sw_output (dfsw, ufsw, fsw, hsw, dfswcf, ufswcf,  &
			   fswcf, hswcf)
     else
       call get_sw_output (dfsw, ufsw, fsw, hsw)
     endif
!-------------------------------------------------------------------
!    for fms, set up diagnostic fields for netcdf files
!-------------------------------------------------------------------
     allocate ( swin   ( is:ie, js:je ) )
     allocate ( swout  ( is:ie, js:je ) )
     allocate ( olr    ( is:ie, js:je ) )
     allocate ( swups  ( is:ie, js:je ) )
     allocate ( swdns  ( is:ie, js:je ) )
     allocate ( lwups  ( is:ie, js:je ) )
     allocate ( lwdns  ( is:ie, js:je ) )
     allocate ( tdtsw  ( is:ie, js:je, kmin:kmax ) )
     allocate ( tdtlw  ( is:ie, js:je, kmin:kmax ) )

     if (Rad_control%do_totcld_forcing) then
       allocate ( swin_clr   ( is:ie, js:je ) )
       allocate ( swout_clr  ( is:ie, js:je ) )
       allocate ( olr_clr    ( is:ie, js:je ) )
       allocate ( swups_clr  ( is:ie, js:je ) )
       allocate ( swdns_clr  ( is:ie, js:je ) )
       allocate ( lwups_clr  ( is:ie, js:je ) )
       allocate ( lwdns_clr  ( is:ie, js:je ) )
       allocate ( tdtsw_clr  ( is:ie, js:je, kmin:kmax ) )
       allocate ( tdtlw_clr  ( is:ie, js:je, kmin:kmax ) )
     endif

!---------------------------------------------------------------------
!  define the netcdf diagnostic arrays.
!---------------------------------------------------------------------
     swin (:,:)   = dfsw(:,:,1)*1.0E-03
     swout(:,:)   = ufsw(:,:,1)*1.0E-03
     olr  (:,:)   = flxnet(:,:,1)*1.0E-03
     swups(:,:)   = ufsw(:,:,kmax+1)*1.0E-03
     swdns(:,:)   = dfsw(:,:,kmax+1)*1.0E-03
     lwups(:,:)   = 1.0E-03*sigmasb*ts(:,:  )**4
     lwdns(:,:)   = lwups(:,:) - flxnet(:,:,kmax+1)*1.0E-03
     tdtsw(:,:,:) = hsw(:,:,:)/secday
     tdtlw(:,:,:) = heatra(:,:,:)/ secday

     if (Rad_control%do_totcld_forcing) then
       swin_clr (:,:)   = dfswcf(:,:,1)*1.0E-03
       swout_clr(:,:)   = ufswcf(:,:,1)*1.0E-03
       olr_clr  (:,:)   = flxnetcf(:,:,1)*1.0E-03
       swups_clr(:,:)   = ufswcf(:,:,kmax+1)*1.0E-03
       swdns_clr(:,:)   = dfswcf(:,:,kmax+1)*1.0E-03
       lwups_clr(:,:)   = 1.0E-03*sigmasb*ts(:,:  )**4
       lwdns_clr(:,:)   = lwups_clr(:,:) - flxnetcf(:,:,kmax+1)*1.0E-03
       tdtsw_clr(:,:,:) = hswcf(:,:,:)/secday
       tdtlw_clr(:,:,:) = heatracf(:,:,:)/secday
     endif

!-----------------------------------------------------------------------
!  save heating rates to the global arrays.
!-----------------------------------------------------------------------
     if (present(mask)) then
       tdt_rad(is:ie,js:je,:)=(tdtsw(:,:,:)+tdtlw(:,:,:))*mask(:,:,:)
       if (Rad_control%do_totcld_forcing) then
         tdt_rad_clr(is:ie,js:je,:)=(tdtsw_clr(:,:,:)+  &
                                     tdtlw_clr(:,:,:))*mask(:,:,:)
       endif
     else
       tdt_rad(is:ie,js:je,:)=(tdtsw(:,:,:)+tdtlw(:,:,:))
       if (Rad_control%do_totcld_forcing) then
         tdt_rad_clr(is:ie,js:je,:)=(tdtsw_clr(:,:,:)+tdtlw_clr(:,:,:))
       endif
     endif

!--------------------------------------------------------------------
!   define the surface long- and short-wave fluxes.
!-----------------------------------------------------------------------
     flux_sw_surf(is:ie,js:je) = swdns(:,:) - swups(:,:)
     flux_lw_surf(is:ie,js:je) = lwdns(:,:)

!---------------------------------------------------------------------
!  define locally dimensioned heating rate arrays to be sent off for
!  later archiving.
!-----------------------------------------------------------------------
     allocate (tdt_rad_loc(is:ie, js:je, kmin:kmax) )
     tdt_rad_loc(:,:,:) = tdt_rad(is:ie,js:je,:)
     if (Rad_control%do_totcld_forcing) then
       allocate (tdt_rad_clr_loc(is:ie, js:je, kmin:kmax) )
       tdt_rad_clr_loc(:,:,:) = tdt_rad_clr(is:ie,js:je,:)
       call hold_tendencies (tdt_rad_loc, tdtsw, tdt_rad_clr_loc,  &
                             tdtsw_clr)
       deallocate (tdt_rad_clr_loc)
     else
       call hold_tendencies (tdt_rad_loc, tdtsw)
     endif
     deallocate (tdt_rad_loc)

!-----------------------------------------------------------------------
!    diagnostics section 
!-----------------------------------------------------------------------

!------- albedo (only do once) -------------------------
      if ( id_alb_sfc > 0 ) then
         allocate (asfc ( is:ie, js:je ) )
         asfc(:,:) = 100.0*albedo   (:,:)
         used = send_data ( id_alb_sfc, asfc, Time_diag, is, js )
      endif

      ipass = 1

!------- sw tendency -----------
      if ( id_tdt_sw(ipass) > 0 ) then
         used = send_data ( id_tdt_sw(ipass), tdtsw, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- lw tendency -----------
      if ( id_tdt_lw(ipass) > 0 ) then
         used = send_data ( id_tdt_lw(ipass), tdtlw, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- incoming sw flux toa -------
      if ( id_swdn_toa(ipass) > 0 ) then
         used = send_data ( id_swdn_toa(ipass), swin, Time_diag, is, js )
      endif

!------- outgoing sw flux toa -------
      if ( id_swup_toa(ipass) > 0 ) then
         used = send_data ( id_swup_toa(ipass), swout, Time_diag, is, js )
      endif

!------- outgoing lw flux toa (olr) -------
      if ( id_olr(ipass) > 0 ) then
         used = send_data ( id_olr(ipass), olr, Time_diag, is, js )
      endif

!------- upward sw flux surface -------
      if ( id_swup_sfc(ipass) > 0 ) then
         used = send_data ( id_swup_sfc(ipass), swups, Time_diag, is, js )
      endif

!------- downward sw flux surface -------
      if ( id_swdn_sfc(ipass) > 0 ) then
         used = send_data ( id_swdn_sfc(ipass), swdns, Time_diag, is, js )
      endif

!------- upward lw flux surface -------
      if ( id_lwup_sfc(ipass) > 0 ) then
         used = send_data ( id_lwup_sfc(ipass), lwups, Time_diag, is, js )
      endif

!------- downward lw flux surface -------
      if ( id_lwdn_sfc(ipass) > 0 ) then
         used = send_data ( id_lwdn_sfc(ipass), lwdns, Time_diag, is, js )
      endif

!-----------------------------------------------------------------------
   if (Rad_control%do_totcld_forcing) then

      ipass = 2

!------- sw tendency -----------
      if ( id_tdt_sw(ipass) > 0 ) then
       used = send_data ( id_tdt_sw(ipass), tdtsw_clr, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- lw tendency -----------
      if ( id_tdt_lw(ipass) > 0 ) then
      used = send_data ( id_tdt_lw(ipass), tdtlw_clr, Time_diag, is, js, 1, &
                            rmask=mask )
      endif

!------- incoming sw flux toa -------
      if ( id_swdn_toa(ipass) > 0 ) then
         used = send_data ( id_swdn_toa(ipass), swin_clr, Time_diag, is, js )
      endif

!------- outgoing sw flux toa -------
      if ( id_swup_toa(ipass) > 0 ) then
         used = send_data ( id_swup_toa(ipass), swout_clr, Time_diag, is, js )
      endif

!------- outgoing lw flux toa (olr) -------
      if ( id_olr(ipass) > 0 ) then
         used = send_data ( id_olr(ipass), olr_clr, Time_diag, is, js )
      endif

!------- upward sw flux surface -------
      if ( id_swup_sfc(ipass) > 0 ) then
!        used = send_data ( id_swup_sfc(ipass), swups_cf, Time_diag, is, js )
         used = send_data ( id_swup_sfc(ipass), swups_clr, Time_diag, is, js )
      endif

!------- downward sw flux surface -------
      if ( id_swdn_sfc(ipass) > 0 ) then
!        used = send_data ( id_swdn_sfc(ipass), swdns_cf, Time_diag, is, js )
         used = send_data ( id_swdn_sfc(ipass), swdns_clr, Time_diag, is, js )
      endif

!------- upward lw flux surface -------
      if ( id_lwup_sfc(ipass) > 0 ) then
!        used = send_data ( id_lwup_sfc(ipass), lwups_cf, Time_diag, is, js )
         used = send_data ( id_lwup_sfc(ipass), lwups_clr, Time_diag, is, js )
      endif

!------- downward lw flux surface -------
      if ( id_lwdn_sfc(ipass) > 0 ) then
!        used = send_data ( id_lwdn_sfc(ipass), lwdns_cf, Time_diag, is, js )
         used = send_data ( id_lwdn_sfc(ipass), lwdns_clr, Time_diag, is, js )
      endif

  endif
!-----------------------------------------------------------------------



!--------------------------------------------------------------------
! accumulate global integral quantities 
!--------------------------------------------------------------------
     call sum_diag_integral_field ('olr',    olr,        is, js)
     call sum_diag_integral_field ('abs_sw', swin-swout, is, js)
!--------------------------------------------------------------------
! accumulate hemispheric integral quantities 
!--------------------------------------------------------------------
     do j=js,je        
       if ( y(3)+j-1 <= y(2)/2 ) then
         call sum_diag_integral_field ('sntop_tot_sh ', swin-swout,   &
		                 is, ie, j, j)
         call sum_diag_integral_field ('lwtop_tot_sh ', olr,     &
		                 is, ie, j, j)
         call sum_diag_integral_field ('sngrd_tot_sh ', swdns-swups, &
		                 is, ie, j, j)
         call sum_diag_integral_field ('lwgrd_tot_sh ', flxnet(:,:,kmax+1)*  &
		                 1.0e-03, is, ie,  j, j)
         if (Rad_control%do_totcld_forcing) then 
           call sum_diag_integral_field ('sntop_clr_sh ', swin_clr-swout_clr, &
		                   is, ie, j, j)
           call sum_diag_integral_field ('lwtop_clr_sh ', olr_clr,     &
		                   is, ie, j, j)
           call sum_diag_integral_field ('sngrd_clr_sh ', swdns_clr-swups_clr, &
		                   is, ie, j, j)
           call sum_diag_integral_field ('lwgrd_clr_sh ', flxnetcf(:,:,kmax+1)*&
		                   1.0e-03, is, ie, j, j)
         endif
       else
         call sum_diag_integral_field ('sntop_tot_nh ', swin-swout,   &
		                 is, ie, j, j)
         call sum_diag_integral_field ('lwtop_tot_nh ', olr,     &
		                 is, ie, j, j)
         call sum_diag_integral_field ('sngrd_tot_nh ', swdns-swups,  &
		                 is, ie, j, j)
         call sum_diag_integral_field ('lwgrd_tot_nh ', flxnet(:,:,kmax+1)*&
		                 1.0e-03, is, ie, j, j)
         if (Rad_control%do_totcld_forcing) then 
           call sum_diag_integral_field ('sntop_clr_nh ', swin_clr-swout_clr,  &
		                   is, ie, j, j)
           call sum_diag_integral_field ('lwtop_clr_nh ', olr_clr,   &
		                   is, ie, j, j)
           call sum_diag_integral_field ('sngrd_clr_nh ', swdns_clr-swups_clr, &
		                   is, ie, j, j)
           call sum_diag_integral_field ('lwgrd_clr_nh ', flxnetcf(:,:,kmax+1)*&
		                   1.0e-03, is, ie, j, j)
         endif
       end if
     end do
!--------------------------------------------------------------------
! accumulate global integral quantities 
!--------------------------------------------------------------------
     call sum_diag_integral_field ('sntop_tot_gl ', swin-swout, is, js)
     call sum_diag_integral_field ('lwtop_tot_gl ', olr,        is, js)
     call sum_diag_integral_field ('sngrd_tot_gl ', swdns-swups,is, js)
     call sum_diag_integral_field ('lwgrd_tot_gl ',  &
       	            flxnet(:,:,kmax+1)*1.0e-03, is, js)
     if (Rad_control%do_totcld_forcing) then 
       call sum_diag_integral_field ('sntop_clr_gl ', swin_clr-swout_clr, &
  			             is,  js)
       call sum_diag_integral_field ('lwtop_clr_gl ', olr_clr,   &
		                     is, js)
       call sum_diag_integral_field ('sngrd_clr_gl ', swdns_clr-swups_clr,&
		                     is,  js)
       call sum_diag_integral_field ('lwgrd_clr_gl ', flxnetcf(:,:,kmax+1)*&
				     1.0e-03, is, js)
     endif

     deallocate (dfsw )
     deallocate (ufsw )
     deallocate ( fsw )
     deallocate ( hsw )
     deallocate (heatra  )
     deallocate (flxnet  )
     deallocate ( swin    )
     deallocate ( swout   )
     deallocate ( olr    )
     deallocate ( swups  )
     deallocate ( swdns  )
     deallocate ( lwups   )
     deallocate ( lwdns  )
     deallocate ( tdtsw )
     deallocate ( tdtlw )
     if (Rad_control%do_totcld_forcing) then
       deallocate (heatracf)
       deallocate (flxnetcf)
       deallocate (dfswcf )
       deallocate (ufswcf )
       deallocate ( fswcf )
       deallocate ( hswcf )
       deallocate ( swin_clr    )
       deallocate ( swout_clr   )
       deallocate ( olr_clr    )
       deallocate ( swups_clr  )
       deallocate ( swdns_clr  )
       deallocate ( lwups_clr   )
       deallocate ( lwdns_clr  )
       deallocate ( tdtsw_clr )
       deallocate ( tdtlw_clr )
     endif
!-----------------------------------------------------------------------



end subroutine radiation_netcdf




!###################################################################

subroutine update_rad_fields (is, ie, js, je, tdt, flux_sw, flux_lw)

!--------------------------------------------------------------------
integer, intent(in)                      :: is, ie, js, je
real,    intent(inout), dimension(:,:,:) :: tdt
real,    intent(out),   dimension(:,:)   :: flux_sw, flux_lw
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  update radiative tendency and fluxes
!-------------------------------------------------------------------
      tdt  (:,:,:) = tdt(:,:,:) + tdt_rad(is:ie,js:je,:)
      flux_sw(:,:) =         flux_sw_surf(is:ie,js:je)
      flux_lw(:,:) =         flux_lw_surf(is:ie,js:je)

!-------------------------------------------------------------------
!  save radiative tendency for use in strat cloud scheme. Note that 
!  if strat cloud scheme is not operating then nothing is done by 
!  this routine.
!-------------------------------------------------------------------
      call add_strat_tend(is,ie,js,je,tdt_rad(is:ie,js:je,:))


end subroutine update_rad_fields 



!####################################################################

subroutine define_rad_input (is, ie, js, je, t, q, pfull, phalf, &
			     ts, kbot)

!------------------------------------------------------------------
integer, intent(in)                           :: is, ie, js, je
real,    intent(in), dimension(:,:,:)         :: pfull, phalf, t, q
real,    intent(in), dimension(:,:)           :: ts
integer, intent(in), dimension(:,:), optional :: kbot
!------------------------------------------------------------------

         integer     :: i, j, k, kb

!--------------------------------------------------------------------
!   allocate space for the input arrays which will be sent to 
!   rad_step_setup_mod from whence they will be usable by all of 
!   the modules in the radiation package.
!--------------------------------------------------------------------
         allocate (ttt_l     (is:ie, js:je, kmin:kmax) )
         allocate (rrr_l     (is:ie, js:je, kmin:kmax) )
         allocate (press3d_l (is:ie, js:je, kmin:kmax+1) )
         allocate (ts_l      (is:ie, js:je) )

!------------------------------------------------------------------
!   put the basic model input arrays into the form expected by the 
!   radiation code.
!   temperature :  3d array(is:ie, js:je, kmin:kmax), degrees Kelvin
!   moisture : 3d array (is:ie, js:je, kmin:kmax), mixing ratio
!   pressure : 3d array (is:ie, js:je, kmin:kmax+1), pressures at model 
!              full levels plus the surface pressure in the kmax+1 slot,
!              in cgs units (FOR NOW !!)
!   surface temperature : 2d array (is:ie, js:je), degrees Kelvin
!------------------------------------------------------------------
         if (Environment%running_gcm) then
           if (Environment%running_fms) then
             ttt_l(is:ie,js:je,:) = t(:,:,:) 
             rrr_l(is:ie,js:je,:) = q(:,:,:)/(1.0-q(:, :, :)/rh2oair)
             press3d_l(:,:,kmin:kmax) = 10.0*   &
					   pfull(:, :,kmin:kmax)
             press3d_l(:,:,kmax+1  ) = 10.0*  &
					   phalf(:, :, kmax+1)

!------------------------------------------------------------------
!  for the eta coordinate case, set the kmax+1 pressure to the value
!  at the surface. set the values of pressure at the below-ground levels
!  to the surface pressure value. set the values of temperature at
!  below-ground levels to the surface temperature value.
!------------------------------------------------------------------
             if (present(kbot)) then
               do j=js,je 
                 do i=is,ie
                   kb=kbot(i,j)
                   if (kb < kmax) then
                     press3d_l(i,j,kmax+1  )=press3d_l(i,j,kb+1  )
                     do k=kb+1,kmax
                       ttt_l(i,j,k  )=ts(i,j )
	               press3d_l(i,j,k  ) = press3d_l(i,j,kmax+1  )
                     enddo
                   endif
                 enddo
               enddo
             endif

!---------------------------------------------------------------------
!   define base variables when running skyhi
!---------------------------------------------------------------------
           else if (Environment%running_skyhi) then
             ttt_l(is:ie,js:je,:) = t(:,:,:) + frezdk
             rrr_l(is:ie,js:je,:) = q(:,:,:)
             press3d_l(is:ie,js:je,kmin:kmax) =    &
					  pfull(:, :, kmin:kmax)
             press3d_l(is:ie,js:je,kmax+1  ) = phalf(:, :, kmax+1)
           endif  ! (running fms)

!---------------------------------------------------------------------
!   define base variables when running standalone code
!---------------------------------------------------------------------
         else if (Environment%running_standalone) then
           ttt_l(is:ie,js:je,:) = t(:,:,:) 
           rrr_l(is:ie,js:je,:) = q(:,:,:)
           press3d_l(is:ie,js:je,kmin:kmax  ) =    &
					 pfull(:, :, kmin:kmax)
           press3d_l(is:ie,js:je,kmax+1  ) = phalf(:, :, kmax+1)
         endif ! (running gcm)

!--------------------------------------------------------------------
!  surface temp is defined the same in all cases
!--------------------------------------------------------------------
         ts_l(is:ie,js:je    )         = ts(:,:)

!---------------------------------------------------------------------



end subroutine define_rad_input



!####################################################################

subroutine define_astro_input (is, ie, js, je)

!-------------------------------------------------------------------
integer, intent(in)           :: is, ie, js, je
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!  obtain the zenith angle and daylight fraction.
!-------------------------------------------------------------------
         allocate (cosz       (1:ie-is+1, 1:je-js+1) )
         allocate (fracday    (1:ie-is+1, 1:je-js+1) )
         call get_astronomy_for_swrad (cosz, fracday)

!-------------------------------------------------------------------
!  obtain the earth-sun distance.
!-------------------------------------------------------------------
	 call get_solar_distance (rrsun)

!----------------------------------------------------------------------


end subroutine define_astro_input


!####################################################################

   subroutine diag_field_init ( Time, axes )

     type(time_type), intent(in) :: Time
     integer        , intent(in) :: axes(4)

     character(len=4) :: clr
     character(len=9) :: clr2
     integer :: i, n

!------------ initialize diagnostic fields in this module --------------

!---- generate names for clear sky diagnostic fields -----
!---- may want to not generate names when clear sky not used? ----

                       n= 1
if (Rad_control%do_totcld_forcing) n= 2

do i = 1, n

    if ( i == 1) then
!  i = 1 is the total-sky case
!   if ( i == n) then
       clr  = "    "
       clr2 = "          "
!   else
!   i = 2 is the clear-sky case
    else if (i == 2) then
       clr  = "_clr"
       clr2 = "clear sky "
    endif

    id_tdt_sw(i) = &
    register_diag_field ( mod_name, 'tdt_sw'//trim(clr), axes(1:3), Time, &
                     trim(clr2)//'temperature tendency for SW radiation', &
                     'deg_K/sec', missing_value=missing_value             )

    id_tdt_lw(i) = &
    register_diag_field ( mod_name, 'tdt_lw'//trim(clr), axes(1:3), Time, &
                     trim(clr2)//'temperature tendency for LW radiation', &
                     'deg_K/sec', missing_value=missing_value             )

    id_swdn_toa(i) = &
    register_diag_field ( mod_name, 'swdn_toa'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'SW flux down at TOA', &
                     'watts/m2', missing_value=missing_value               )

    id_swup_toa(i) = &
    register_diag_field ( mod_name, 'swup_toa'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'SW flux up at TOA', &
                     'watts/m2', missing_value=missing_value               )

    id_olr(i) = &
    register_diag_field ( mod_name, 'olr'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'outgoing longwave radiation', &
                     'watts/m2', missing_value=missing_value               )

    id_swup_sfc(i) = &
    register_diag_field ( mod_name, 'swup_sfc'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'SW flux up at surface', &
                     'watts/m2', missing_value=missing_value               )

    id_swdn_sfc(i) = &
    register_diag_field ( mod_name, 'swdn_sfc'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'SW flux down at surface', &
                     'watts/m2', missing_value=missing_value               )

    id_lwup_sfc(i) = &
    register_diag_field ( mod_name, 'lwup_sfc'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'LW flux up at surface', &
                     'watts/m2', missing_value=missing_value               )

    id_lwdn_sfc(i) = &
    register_diag_field ( mod_name, 'lwdn_sfc'//trim(clr), axes(1:2), Time, &
                     trim(clr2)//'LW flux down at surface', &
                     'watts/m2', missing_value=missing_value               )

enddo


     id_alb_sfc = &
       register_diag_field ( mod_name, 'alb_sfc', axes(1:2), Time, &
                    'surface albedo', 'percent'            )

     if (Environment%running_gcm .and. Environment%running_fms) then
       id_cosz = register_diag_field    &
		                (mod_name,'cosz',axes(1:2),  &
		                 Time, 'cosz',&
                                 'unknown',missing_value=missing_value)
       id_fracday = register_diag_field    &
				(mod_name,'fracday',axes(1:2),   &
				 Time, 'fracday',&
                                 'unknown',missing_value=missing_value)
     endif
     
        
!-------------------------------------------------------------------

end subroutine diag_field_init
        
      
!####################################################################
      
      
        
                  end module sea_esf_rad_mod
      
      
