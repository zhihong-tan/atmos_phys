                       module ozone_mod

use time_manager_mod,    only:  time_type, days_in_year, &
!                               time_manager_init, &
                                get_time, get_date
use   time_interp_mod,   only:  fraction_of_year
!use   time_interp_mod,   only:  fraction_of_year, &
!                               time_interp_init  ! doesn't exist yet
use  utilities_mod,      only:  open_file, file_exist,    &
                                check_nml_error, error_mesg,  &
                                utilities_init, &
                                print_version_number, FATAL, NOTE, &
                                WARNING, get_my_pe, close_file
use rad_utilities_mod,   only:  Rad_control, radiation_control_type, &
                                rad_utilities_init, &
                                radiative_gases_type,   &
                                atmos_input_type, &
                                Environment, environment_type
use constants_mod,       only:  GRAV, radian
use interpolator_mod,    only:  interpolate_type, interpolator_init, &
                                interpolator, interpolator_end, &
                                CONSTANT, INTERP_WEIGHTED_P


implicit none
private

!-------------------------------------------------------------------
!    ozone_mod supplies mass mixing ratios of ozone (g/g) to the 
!    sea_esf_rad radiation_package (and the original_fms_rad package).
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

character(len=128)  :: version =  '$Id: ozone.F90,v 1.4 2003/04/09 21:01:03 fms Exp $'
character(len=128)  :: tag     =  '$Name: inchon $'



!---------------------------------------------------------------------
!-------  interfaces --------

public  ozone_init, ozone_driver, ozone_end


private   &
        obtain_annual_mean_data, &
        obtain_seasonal_data, &
        obtain_randel_data, &
        obtain_input_file_data, &
!        obtain_zonal_ozone_data, &
        obtain_fms_zonal_ozone_data, &
        obtain_clim_zonal_ozone_data, &
        obtain_photochem_adjust_data, &
        obtain_ozone_loss_data, & 
        find_nearest_lower_index, &
        obtain_current_ozone, &
        harmonic_seasonal_interpolation, &
        interpolate_monthly_data, &
!       zonal_ozone
        fms_zonal_ozone, &
        get_clim_ozone

!interface zonal_ozone
interface fms_zonal_ozone
    module procedure geto3_2d, geto3_3d
end interface


!---------------------------------------------------------------------
!-------- namelist  ---------

logical            :: o3monloss=.false.
character(len=24)  :: o3monloss_type = '    '
!character(len=24)  :: basic_ozone_type= 'zonal_ozone'
!character(len=24)  :: zonal_ozone_type= 'seasonally_varying'
character(len=24)  :: basic_ozone_type= 'fms_zonal_ozone'
character(len=24)  :: fms_zonal_ozone_type= 'seasonally_varying'
logical            :: do_ozone_adjustment = .false.
logical            :: do_mcm_o3_clim = .false.
real               :: ozone_adjustment_limit = 1000.0 
                                 ! 10 mb, pressure above which ozone 
                                 ! photochemical adjustment is applied, 
                                 ! if activated


namelist /ozone_nml/             &
                       basic_ozone_type, &
                       o3monloss, o3monloss_type, &
                       ozone_adjustment_limit, &
!                      zonal_ozone_type, &
                       fms_zonal_ozone_type, &
                       do_ozone_adjustment, &
                       do_mcm_o3_clim

!---------------------------------------------------------------------
!------- public data ------



!---------------------------------------------------------------------
!------- private data ------

!--------------------------------------------------------------------
!    xo3ann contains the annual mean ozone mass mixing ratio (10**4g/g).
!--------------------------------------------------------------------
real, dimension(:,:), allocatable    ::  xo3ann

!--------------------------------------------------------------------
!    xbar, a1, b1, b2, xbarp, a1p, b1p, b2p contain the components
!    needed for harmonic time interpolation of the seasonal ozone
!    data. units are 10**4 g/g. displacement is the fractional distance
!    of model grid points from the nearest lower latitude of the input 
!    data set.
!--------------------------------------------------------------------
real, dimension(:,:), allocatable    ::  xbar, a1, b1, b2, &
                                         xbarp, a1p, b1p, b2p
real, dimension(:),   allocatable    ::  displacement

!------------------------------------------------------------------
!    qo3_randel contains ozone mass mixing ratio in g/g as a function 
!    of latitude, height and month of year from the randel data set.
!--------------------------------------------------------------------
real, dimension(:,:,:), allocatable  ::  qo3_randel

!-------------------------------------------------------------------
!    qqo3 contains the column input ozone mass mixing ratio (g/g).
!-------------------------------------------------------------------
real, dimension (:), allocatable     ::  qqo3

!-------------------------------------------------------------------
!    o3data contains the zonal ozone mass mixing ratio (10**6 g/g)
!    at 19 specified latitudes, 81 specified pressures and up to 4
!    specified times. rstd contains the zonal ozone mass mixing ratio 
!    (10**6 g/g) at 19 specified latitudes, 81 specified pressures and 
!    at the currently specified time. ph is the array defining the 
!    interface levels in the zonal ozone data set.
!-------------------------------------------------------------------
real, dimension (82)        ::    ph
real, dimension (19,81,4)   ::    o3data
real, dimension (19,81)     ::    rstd

!--------------------------------------------------------------------
!    o3depl, o3depb and o3dept define an ozone depletion profile, 
!    either from toms1 or toms2. tomsprofile1 values assume depletion 
!    from the arbitrarily-specified tropopause (as in RSS) to 7 km 
!    above. tomsprofile2 values assume depletion from the arbitrarily-
!    specified tropopause (as in RSS) to 10 km above.
!    o3depl is the o3 column change in dobson units as a function
!    of latitude and month. 
!    o3dept is the index of top layer in which ozone change is to occur.
!    varies with season and latitude.
!    o3depb is the index of bottom layer in which ozone change is to 
!    occur. varies with latitude.
!---------------------------------------------------------------------
real,    dimension(:,:), allocatable ::  o3depl
integer, dimension(:,:), allocatable ::  o3dept
integer, dimension(:),   allocatable ::  o3depb

!----------------------------------------------------------------------
!    ztrad is the reference temperature distribution at model grid 
!    points corresponding to the input ozone distribution.
!----------------------------------------------------------------------
real, dimension(:,:), allocatable    ::  ztrad

!----------------------------------------------------------------------
!    o3 is an interpolate_type variable containing the relevant 
!    information about the clim_zonal_ozone data set.
!----------------------------------------------------------------------
type(interpolate_type)               ::  o3


integer     ::  kobeg      ! lowest vertical index at which photo-
                           ! chemical adjustment is applied
integer     ::  koend      ! highest vertical index at which photo-
                           ! chemical adjustment is applied
integer     ::  kmin=1     ! lowest index in model vertical grid
integer     ::  kmax       ! highest index in model vertical grid
integer     ::  jdf        ! number of latitudes on processor
integer     ::  iseason=-1 ! flag indicating type of zonal ozone time
                           !  variation to use
real        ::  current_fyear=-1.0  
                           ! flag to force interpolation on initial call
logical     ::  do_tomsprofile1 =.false.  
                           ! using toms 1 ozone loss profile ?
logical     ::  do_tomsprofile2 =.false.
                           ! using toms 2 ozone loss profile ?
logical     ::  do_randelprofile=.false.
                           ! using randel ozone data set ?
!logical     ::  do_zonal_ozone=.false.
!                          ! using FMS zonal ozone data set ?
logical     ::  do_fms_zonal_ozone=.false.
                           ! using FMS zonal ozone data set ?
logical     ::  do_clim_zonal_ozone=.false.
                           ! using clim zonal ozone data set ?
logical     ::  do_annual_mean_ozone=.false.
                           ! using annual mean ozone data set ?
logical     ::  do_seasonal_ozone=.false.
                           ! using seasonal ozone data set ?
logical     ::  do_column_input=.false.
                           ! using column ozone data set ?
logical     ::  ozone_initialized=.false.  ! module initialized ?


!-------------------------------------------------------------------
!-------------------------------------------------------------------



                          contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!subroutine ozone_init (pref, latb)
subroutine ozone_init (pref, latb, lonb)

!---------------------------------------------------------------------
!    ozone_init is the constructor for ozone_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
real, dimension(:,:), intent(in) :: pref
!real, dimension(:),   intent(in) :: latb
real, dimension(:),   intent(in) :: latb, lonb
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       pref      array containing two reference pressure profiles 
!                 [pascals]
!       latb      array of model latitudes at cell boundaries [radians]
!       lonb      array of model longitudes at cell boundaries [radians]
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  local variables:

       integer           ::  unit, ierr, io

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!-------------------------------------------------------------------
      if (ozone_initialized) return

!-------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!-------------------------------------------------------------------
      call utilities_init
      call rad_utilities_init
!     call time_manager_init   ! doesn't exist yet
!     call time_interp_init    ! doesn't exist yet
 
!----------------------------------------------------------------
!    read namelist.
!----------------------------------------------------------------
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=ozone_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'ozone_nml')
        enddo
10      call close_file (unit)
      endif

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      unit = open_file ('logfile.out', action='append')
      if (get_my_pe() == 0) then
!     if (get_my_pe() == get_root_pe()) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit,nml=ozone_nml)
      endif
      call close_file (unit)

!-------------------------------------------------------------------
!    define the number of layers in the model and the number of 
!    latitude rows local to the processor.
!-------------------------------------------------------------------
      kmax = size(pref, 1) - 1
      jdf = size (latb) - 1

!---------------------------------------------------------------------
!    when running standalone code, basic_ozone_type must be 'input'.
!---------------------------------------------------------------------
       if ((Environment%running_standalone  .and.       &
            .not. Environment%running_sa_model) .and.   &
             trim(basic_ozone_type) /= 'input') then
        call error_mesg ('ozone_mod', &
            ' must supply ozone via input file when running &
	        &standalone code', FATAL)
      endif

!------------------------------------------------------------------
!    use the annual mean ozone data. this reads a specific file with
!    ozone data in a lat x ht array dimensioned (19, 40), with the 
!    latitudes being every 10 degrees and the levels being the 40 
!    level skyhi grid. it is currently only usable with that vertical
!    grid, although a horizontal interpolation ability is provided.
!-------------------------------------------------------------------
      if (trim(basic_ozone_type) == 'annual_mean' ) then
        call obtain_annual_mean_data (latb)

!------------------------------------------------------------------
!    use the seasonal ozone data. this reads a specific file with
!    ozone data in lat x ht arrays dimensioned (19, 40), with the 
!    latitudes being every 10 degrees and the levels being the 40 
!    level skyhi grid. it is currently only usable with that vertical
!    grid, although a horizontal interpolation ability is provided.
!-------------------------------------------------------------------
      else if (trim(basic_ozone_type) == 'seasonally_varying' ) then
        call obtain_seasonal_data (latb)

!------------------------------------------------------------------
!    use the ozone data file for the randel profile case. the randel
!    ozone profile is created based on the depletions given by randel 
!    from a standard ozone profile, in this case taken from skyhi 
!    zonally-averaged ozone profiles, as a function of height and 
!    month of year.  it is currently only usable with an n3o skyhi
!    grid.
!-----------------------------------------------------------------
      else if (trim(basic_ozone_type) == 'randelprofile' ) then
        call obtain_randel_data (latb) 

!------------------------------------------------------------------
!    read a single column ozone profile.
!-----------------------------------------------------------------
      else if (trim(basic_ozone_type) == 'input' ) then
        call obtain_input_file_data 

!------------------------------------------------------------------
!    use the original fms ozone module.
!-----------------------------------------------------------------
!      else if (trim(basic_ozone_type) == 'zonal_ozone' ) then
      else if (trim(basic_ozone_type) == 'fms_zonal_ozone' ) then
!        do_zonal_ozone = .true.
        do_fms_zonal_ozone = .true.
!       if (trim(zonal_ozone_type)      == 'winter' ) then
        if (trim(fms_zonal_ozone_type)      == 'winter' ) then
          iseason = 1   
        else if (trim(fms_zonal_ozone_type) == 'spring' ) then
          iseason = 2   
        else if (trim(fms_zonal_ozone_type) == 'summer' ) then
          iseason = 3   
        else if (trim(fms_zonal_ozone_type) == 'autumn' ) then
          iseason = 4   
        else if (trim(fms_zonal_ozone_type) == 'annual_mean' ) then
          iseason = 0   
        else if (trim(fms_zonal_ozone_type) == 'seasonally_varying' ) then
          iseason = 5   
        else
          call error_mesg ( 'ozone_mod', &
         'improper specification of nml variable zonal_ozone_type',   &
                                                            FATAL)
        endif
        call obtain_fms_zonal_ozone_data (iseason)

!------------------------------------------------------------------
!    use the new zonal ozone climatology.
!-----------------------------------------------------------------
      else if (trim(basic_ozone_type) == 'clim_zonal_ozone' ) then
        do_clim_zonal_ozone = .true.
        call obtain_clim_zonal_ozone_data (lonb, latb)

!------------------------------------------------------------------
!    error condition -- improper basic_ozone_type.
!-----------------------------------------------------------------
      else
        call error_mesg ( 'ozone_init',    &
               ' basic_ozone_type is not an acceptable value.', FATAL)
      endif

!-----------------------------------------------------------------------
!    if adjustment of the ozone profile to attempt to account for
!    the effects of ozone photochemistry is desired, obtain the 
!    necessary input data.
!-----------------------------------------------------------------------
!     if (do_ozone_adjustment .and. do_zonal_ozone) then
      if (do_ozone_adjustment .and. (do_fms_zonal_ozone .or.    &
                                     do_clim_zonal_ozone) ) then
        call error_mesg ( 'ozone_mod', &
      'ozone adjustment is not currently available with zonal ozone', &
                                                            FATAL)
      else if (do_ozone_adjustment) then
        call obtain_photochem_adjust_data (latb, pref)  
      endif  

!---------------------------------------------------------------------
!    if a modification to the ozone profile is to be made based on 
!    toms data, call obtain_ozone_loss_data to determine the loss 
!    distribution.
!--------------------------------------------------------------------
!     if (o3monloss .and. do_zonal_ozone) then
      if (o3monloss .and. (do_fms_zonal_ozone .or. &
                           do_clim_zonal_ozone) ) then
        call error_mesg ( 'ozone_mod', &
       'ozone depletion is not currently available with zonal_ozone', &
                                                            FATAL)
      else if (o3monloss) then
        call obtain_ozone_loss_data (latb, pref)
      endif   

!---------------------------------------------------------------------
      ozone_initialized = .true.
    


end subroutine ozone_init



!######################################################################

subroutine ozone_driver (is, ie, js, je, lat, Rad_time, Atmos_input, &
                         Rad_gases )
 
!--------------------------------------------------------------------
!   ozone_driver obtains the current ozone distributions and returns 
!   them in Rad_gases%qo3.
!--------------------------------------------------------------------

integer,                    intent(in)  :: is, ie, js, je
real, dimension(:,:),       intent(in)  :: lat
type(time_type),            intent(in)  :: Rad_time
type(atmos_input_type),     intent(in)  :: Atmos_input
type(radiative_gases_type), intent(inout) :: Rad_gases

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      lat          latitude of model points  [ radians ]
!      Rad_time     time at which the climatologically-determined,
!                   time-varying ozone field should apply
!                   [ time_type (days, seconds) ] 
!      Atmos_input  atmos_input_type variable containing the atmospheric
!                   input fields needed by the radiation package
!
!   intent(out) variables:
!
!      Rad_gases    radiative_gases_type variable which will return
!                   the ozone mass mixing ratio (g/g) to the calling
!                   routine
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  local variables:
!
!      phaf         pressure at model interface levels, normalized
!                   by the mean sea level pressure (101325 N/m**2).
!                   [ N/m**2 ]
!      em           adjustment factor to ozone mixing ratio as a 
!                   result of adjusting temperature from that con-
!                   sistent with the assumed ozone distribution to the 
!                   actual model temperature field [ dimensionless ]
!      vx1          photochemical adjustment factor
!      vx2          photochemical adjustment factor
!      tcolo3       total column ozone in region over which
!                   toms-based ozone depletion is to apply
!      ztrad1       factor used in photochemical damping
!      ztrad2       factor used in photochemical damping
!      zdduo3n      ozone amounts (mass mixing ratio) at model 
!                   latitudes and levels at specified time.
!      zdoblos      ozone depletion at the specified time
!      o3ducst      conversion factor from dobson units to cgs units
!      tmean        mean temperature used in damping calculations
!      fudgek       fudge factor
!      bb1          constant used in photochemical adjustment
!      bb2          constant used in photochemical adjustment
!      con1         factor   used in photochemical adjustment
!
!---------------------------------------------------------------------
      real, dimension (size(Atmos_input%press,1), &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3))    :: phaf, em,   &
                                                        vx1, vx2

      real, dimension (size(Atmos_input%press,1), &
                       size(Atmos_input%press,2))    :: tcolo3

      real, dimension (size(ztrad,1),             &
                       size(ztrad,2))                :: ztrad1, ztrad2

      real, dimension (size(Atmos_input%press,2), &
                       size(Atmos_input%press,3)-1)  :: zdduo3n

      real, dimension (size(Atmos_input%press,2) )   :: zdoblos

      real      ::  o3ducst = 466.968 
      real      ::  tmean   = 265.0   
      real      ::  fudgek  = 4.0 
      real      ::  bb1     = 510.0 
      real      ::  bb2     = 2300.0 
      real      ::  con1    

      integer   ::  k,j,i   ! do-loop indices 

!--------------------------------------------------------------------
      if ( .not. ozone_initialized)    & 
        call error_mesg ('ozone_mod',  &         
           'module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    define the values of ozone valid at the specified time, and any
!    ozone depletion amounts which are to be applied.
!--------------------------------------------------------------------
      call obtain_current_ozone (js, je, lat, Rad_time, zdoblos,  &
                                 zdduo3n)
 
!-------------------------------------------------------------------
!    if an adjustment is to be made to the ozone distribution to 
!    reflect the effect of temperature change from the temperatures
!    for which the input field is compatible, determine this adjustment
!    and apply it. it will be applied between levels kobeg and koend.
!-------------------------------------------------------------------
      if (do_ozone_adjustment) then
        con1 = fudgek*EXP(-bb2/tmean)
        do k=kobeg,koend
          do j=1,jdf
            ztrad1(j,k) = EXP(-bb1/ztrad(j,k))
            ztrad2(j,k) = EXP(-bb2/ztrad(j,k))
          end do
        end do
        em (:,:,kmin:kmax  ) = 1.0E+00
        vx1(:,:,kobeg:koend) = EXP(bb1/   &
                                     Atmos_input%temp(:,:,kobeg:koend))
        vx2(:,:,kobeg:koend) = EXP(-bb2/  &
                                     Atmos_input%temp(:,:,kobeg:koend))
        vx2(:,:,kobeg:koend) = 1.0E+00/(con1 + vx2(:,:,kobeg:koend))
        do k=kobeg,koend
          do j=1,je-js+1      
            do i=1,ie-is+1      
              em (i,j,k) = SQRT(ztrad(j+js-1 ,k)/   &
                                Atmos_input%temp(i,j,k)*(vx1(i,j,k)*   &
                                ztrad1(j+js-1 ,k)*(con1 +    &
                                ztrad2(j+js-1 ,k))*vx2(i,j,k)))
            end do
          end do
        end do
      endif

!---------------------------------------------------------------------
!    allocate space for the ozone mixing ratio array.
!---------------------------------------------------------------------
      allocate ( Rad_gases%qo3 (ie-is+1, je-js+1, kmax) )
       Rad_gases%qo3  = 0.

!-----------------------------------------------------------------------
!    define the (i,j,k) ozone field valid at the specified time, 
!    including the effect of the photochemical adjustment, if desired.
!-----------------------------------------------------------------------
!     if ( .not. do_zonal_ozone) then 
      if ( .not. (do_fms_zonal_ozone .or. do_clim_zonal_ozone) ) then 
        if (do_ozone_adjustment) then
          do k=kmin,kmax
            do j=1,je-js+1    
              do i=1,ie-is+1    
                Rad_gases%qo3(i,j,k) = zdduo3n(j,k)*em(i,j,k)
              end do
            end do
          end do
        else 
          do k=kmin,kmax
            do j=1,je-js+1      
              do i=1,ie-is+1      
                Rad_gases%qo3(i,j,k) = zdduo3n(j      ,k)
              end do
            end do
          end do
        endif
      else 
        do k=kmin,kmax+1
          if (do_mcm_o3_clim) then
            phaf(:,:,k) = 100000.*Atmos_input%phalf(:,:,k)/Atmos_input%phalf(:,:,kmax+1)
          else
            phaf(:,:,k) = (Atmos_input%pflux(:,:,k))*101325./   &
                          (Atmos_input%pflux(:,:,kmax +1))
          endif
        end do
        if (do_fms_zonal_ozone) then
          call fms_zonal_ozone (Rad_time, lat, phaf, Rad_gases%qo3)
        else if (do_clim_zonal_ozone) then
          call get_clim_ozone (Rad_time, phaf, Rad_gases%qo3, is, js)
        endif
      endif

!----------------------------------------------------------------------
!    if ozone is to be depleted, calculate the value to be used.
!    evaluate the column ozone, in dobson units, in the skyhi layers
!    (idept-idepb) (ie, those in which ozone change is to occur)
!    the altitude (index values) of the depletion layers is set, as in 
!    the annually averaged case, at the april values. 
!----------------------------------------------------------------------
      if (o3monloss) then
        tcolo3(:,:) = 0.0E+00
        do j=1,je-js+1    
          do i=1,ie-is+1      
            do k=o3dept(j+js-1,2), o3depb(j+js-1 )
              tcolo3(i,j) = tcolo3(i,j) + o3ducst*1.0E+02*  &
                            Rad_gases%qo3(i,j,k)* &
                            0.5E+00*(Atmos_input%press(i,j,k+1) -  &
                            Atmos_input%press(i,j,k-1))/(GRAV)
            end do
            do k=o3dept(j+js-1,2), o3depb(j+js-1 )
              Rad_gases%qo3(i,j,k) = Rad_gases%qo3(i,j,k)*   &
                                    (1.0E+00 + zdoblos(j)/tcolo3(i,j))
            end do
          end do
        end do
      endif

!--------------------------------------------------------------------




end subroutine ozone_driver



!#####################################################################

subroutine ozone_end

!---------------------------------------------------------------------
!    ozone_end is the destructor for ozone_mod.
!---------------------------------------------------------------------
     
!---------------------------------------------------------------------
!    deallocate any module variables that have been allocated.
!---------------------------------------------------------------------
      if (allocated (xo3ann      ) ) deallocate ( xo3ann        )
      if (allocated (displacement) ) deallocate ( displacement  )
      if (allocated (xbar        ) ) deallocate ( xbar          )
      if (allocated (a1          ) ) deallocate ( a1            )
      if (allocated (b1          ) ) deallocate ( b1            )
      if (allocated (b2          ) ) deallocate ( b2            )
      if (allocated (xbarp       ) ) deallocate ( xbarp         )
      if (allocated (a1p         ) ) deallocate ( a1p           )
      if (allocated (b1p         ) ) deallocate ( b1p           )
      if (allocated (b2p         ) ) deallocate ( b2p           )
      if (allocated (qo3_randel  ) ) deallocate ( qo3_randel    )
      if (allocated (qqo3        ) ) deallocate ( qqo3          )
      if (allocated (ztrad       ) ) deallocate ( ztrad         )
      if (allocated (o3depb      ) ) deallocate ( o3depb        )
      if (allocated (o3dept      ) ) deallocate ( o3dept        )
      if (allocated (o3depl      ) ) deallocate ( o3depl        )

!---------------------------------------------------------------------
      if (do_clim_zonal_ozone) then
        call interpolator_end (o3)
      endif

!---------------------------------------------------------------------
      ozone_initialized = .false.



end subroutine ozone_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!######################################################################

subroutine obtain_annual_mean_data (latb)

!--------------------------------------------------------------------
!    obtain_annual_mean_data reads an input file containing an annual
!    mean ozone distribution as a function of latitude and height.
!--------------------------------------------------------------------

real, dimension(:), intent(in) :: latb

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       latb      array of model latitudes at cell boundaries [radians]
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  local variables:
!
!       num_level_records     number of data levels in input file
!       num_latitude_records  number of latitude records in data file
!       initial_lat           latitude of first data record in data file
!       delta_lat             latitudinal increment between data records
!       hemi_data             does the data set contain only 
!                             hemispheric data (in contrast to global) ?
!       xo3ann_r              array into which the input data is read
!                             [ 10**4 g / g ]
!       jindx                 latitude index in the data set closest to
!                             but less than the model latitude.
!       displacement          distance of model latitude from data lat-
!                             itude jindx relative to model grid delta
!
!---------------------------------------------------------------------
      integer, parameter   :: NUM_LEVEL_RECORDS = 40 
      integer, parameter   :: NUM_LAT_RECORDS = 19 
      real                 :: initial_lat = -90.0   
      real                 :: delta_lat = 10.0    
      logical              :: hemi_data = .false.

      real,    dimension(NUM_LAT_RECORDS,     &
                         NUM_LEVEL_RECORDS) :: xo3ann_r
      integer, dimension(size(latb)-1)      :: jindx
      real   , dimension(size(latb)-1)      :: displacement


      integer           ::  iounit ! unit number for reading input file
      integer           ::  k, j   ! various indices
      real              ::  inilat ! lowest latitude present in data set
                                   ! [ radians ]
      real              ::  dellat ! latitudinal resolution of data 
                                   ! set [ radians ]

!---------------------------------------------------------------------
!    verify that the model and input data are both on a 40 level
!    vertical grid. if not, write an error message, since interpolation
!    to a different grid is not currently available. if they both are 
!    on a 40 level grid, send a note to remind the user that this data 
!    set is applicable to the 40 level skyhi grid, and that if a 
!    different 40 level grid is being used, the data is incorrect.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!! THIS DOES NOT PREVENT A 40 level NON-SKYHI GRID FROM USING THIS 
!!! OPTION.  BEWARE !!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!---------------------------------------------------------------------
      if (kmax /= 40) then
        call error_mesg ('ozone_mod', &
           'can only use annual_mean ozone with 40 level &
           &skyhi grid -- no vertical interpolation available.', &
                                                       FATAL)
      else
        call error_mesg ('ozone_mod', &
           'ARE YOU USING THE 40 level SKYHI GRID?? IF NOT,&
           & please stop the run since annual mean ozone &
           &profiles will be incorrect!!!', NOTE)
      endif

!---------------------------------------------------------------------
!    set a flag indicating that annual mean ozone data is to be used.
!---------------------------------------------------------------------
      do_annual_mean_ozone = .true.

!---------------------------------------------------------------------
!    determine if the annual mean ozone input data file exists. if so, 
!    open the file and read the data set.  close the file upon 
!    completion.
!---------------------------------------------------------------------
      if (file_exist ( 'INPUT/annual_mean_ozone') ) then
        iounit = open_file ('INPUT/annual_mean_ozone',   &
                            form='formatted', action='read')
        do k=1,NUM_LEVEL_RECORDS
          read (iounit, FMT = '(4e13.5)')   &
	                  (xo3ann_r(j,k), j=1,NUM_LAT_RECORDS)
        end do
        call close_file (iounit)

!---------------------------------------------------------------------
!    if file is not present, write an error message.
!---------------------------------------------------------------------
      else
        call error_mesg ( 'ozone_mod', &
            'annual_mean_ozone data file is not present',FATAL)
      endif

!--------------------------------------------------------------------
!    find the location of model latitudes with respect to the latitude
!    grid that the input data is on.
!--------------------------------------------------------------------
      inilat = initial_lat/radian
      dellat    = delta_lat/radian
      call find_nearest_lower_index (latb, hemi_data, inilat, dellat,&
                                     jindx, displacement)

!--------------------------------------------------------------------
!    allocate space and then define the ozone mass mixing ratio 
!    (10**4 g/g) at the model grid points. 
!--------------------------------------------------------------------
      allocate (xo3ann (jdf, kmax))
      do j=1,jdf
        xo3ann(j,:) = xo3ann_r(jindx(j),:) + displacement(j)*  &
                      (xo3ann_r(jindx(j)+1,:) - xo3ann_r(jindx(j),:)) 
      end do

!---------------------------------------------------------------------


end subroutine obtain_annual_mean_data 



!###################################################################

subroutine obtain_seasonal_data (latb)

!---------------------------------------------------------------------
!    obtain_seasonal_data reads an input file containing seasonal
!    ozone data as a function of latitude and height and retrieves 
!    the data needed by the given processor.
!---------------------------------------------------------------------

real, dimension(:), intent(in) :: latb

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       latb      array of model latitudes at cell boundaries [radians]
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  local variables:
!
!       num_level_records     number of data levels in input file
!       num_latitude_records  number of latitude records in data file
!       initial_lat           latitude of first data record in data file
!       delta_lat             latitudinal increment between data records
!       hemi_data             does the data set contain only 
!                             hemispheric data (in contrast to global) ?
!       xo3ws                 arrays containing the winter - summer 
!                             ozone data [ 10**4 g / g ]
!       xo3sf                 arrays containing the spring - fall ozone
!                             data [ 10**4 g / g ]
!       jindx                 latitude index in the data set closest to
!                             but less than the model latitude.
!
!---------------------------------------------------------------------
      integer, parameter   :: NUM_LEVEL_RECORDS = 40 
      integer, parameter   :: NUM_LAT_RECORDS = 19 
      real                 :: initial_lat = -90.0
      real                 :: delta_lat = 10.0
      logical              :: hemi_data = .false.

      real,    dimension (NUM_LAT_RECORDS, NUM_LEVEL_RECORDS) :: xo3ws
      real,    dimension (NUM_LAT_RECORDS, NUM_LEVEL_RECORDS) :: xo3sf
      integer, dimension(size(latb)-1)                        :: jindx

      integer           ::  iounit ! unit for reading input file
      integer           ::  k, j   ! various indices
      real              ::  inilat ! lowest latitude present in data set
                                   ! [ radians ]
      real              ::  dellat ! latitudinal resolution of data 
                                   ! set [ radians ]

!---------------------------------------------------------------------
!    verify that the model and input data are both on a 40 level
!    vertical grid. if not, write an error message, since interpolation
!    to a different gridis not currently available. if they are on the
!    same grid, send a note to remind the user that this data set is 
!    applicable to the 40 level skyhi grid, and that if a different 40 
!    level grid is being used, the data is incorrect.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!! THIS DOES NOT PREVENT A 40 level NON-SKYHI GRID FROM USING THIS 
!!! OPTION.  BEWARE !!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!---------------------------------------------------------------------
      if (kmax /= 40) then
        call error_mesg ('ozone_mod', &
            'can only use annual_mean ozone with 40 level &
            &skyhi grid -- no vertical interpolation available.', &
                                                         FATAL)
      else
        call error_mesg ('ozone_mod', &
            'ARE YOU USING THE 40 level SKYHI GRID?? IF NOT,&
            & please stop the run since seasonally-varying &
	       &ozone profiles will be incorrect!!!', NOTE)
      endif

!---------------------------------------------------------------------
!    set a flag indicating that the seasonal ozone data is to be used.
!---------------------------------------------------------------------
      do_seasonal_ozone = .true.

!---------------------------------------------------------------------
!    determine if the seasonal ozone input data file exists. if so, 
!    open the file and read the data set.  close the file upon 
!    completion.
!---------------------------------------------------------------------
      if (file_exist ( 'INPUT/seasonal_ozone') ) then
        iounit = open_file ('INPUT/seasonal_ozone',   &
                            form='formatted', action='read')
        do k=1,NUM_LEVEL_RECORDS
          read (iounit, FMT = '(4e18.10)')    &
	                             (xo3sf(j,k), j=1,NUM_LAT_RECORDS)
        end do
        do k=1,NUM_LEVEL_RECORDS
          read (iounit, FMT = '(4e18.10)')     &                      
                          (xo3ws(j,k), j=1,NUM_LAT_RECORDS)
        end do
        call close_file (iounit)

!---------------------------------------------------------------------
!    if file is not present, write an error message.
!---------------------------------------------------------------------
      else
        call error_mesg ( 'ozone_init', &
            'seasonal_ozone data file is not present',FATAL)
      endif

!--------------------------------------------------------------------
!    find the location of model latitudes with respect to the latitude
!    grid that the input data is on.
!--------------------------------------------------------------------
      inilat = initial_lat/radian
      dellat = delta_lat/radian
      allocate (displacement (jdf))
      call find_nearest_lower_index (latb, hemi_data, inilat, dellat,&
                                     jindx, displacement)
!---------------------------------------------------------------------
!    allocate space for and then define the components of the ozone 
!    field needed for harmonic interpolation.
!---------------------------------------------------------------------
      allocate (xbar        (jdf, kmax))
      allocate (a1          (jdf, kmax))
      allocate (b1          (jdf, kmax))
      allocate (b2          (jdf, kmax))
      allocate (xbarp       (jdf, kmax))
      allocate (a1p         (jdf, kmax))
      allocate (b1p         (jdf, kmax))
      allocate (b2p         (jdf, kmax))

      do j=1,jdf
        xbar(j,:) = 0.25*(xo3ws(20-jindx(j),:) + &
                          xo3ws(jindx(j),:) + &
                          xo3sf(20-jindx(j),:) + &
                          xo3sf(jindx(j),:) )   
        a1(j,:) = 0.50*(xo3sf(20-jindx(j),:) - &
                        xo3sf(jindx(j),:) )   
        b1(j,:) = 0.50*(xo3ws(20-jindx(j),:) - &
                        xo3ws(jindx(j),:) )   
        b2  (j,:) = 0.25*(xo3ws(20-jindx(j),:) + &
                          xo3ws(jindx(j),:) - &
                          xo3sf(20-jindx(j),:) - &
                          xo3sf(jindx(j),:) )   
        xbarp(j,:) = 0.25*(xo3ws(20-(jindx(j)+1),:) + &
                           xo3ws(jindx(j)+1,:) + &
                           xo3sf(20-(jindx(j)+1),:) + &
                           xo3sf(jindx(j)+1,:) )   
        a1p(j,:) = 0.50*(xo3sf(20-(jindx(j)+1),:) - &
                         xo3sf(jindx(j)+1,:) )   
        b1p(j,:) = 0.50*(xo3ws(20-(jindx(j)+1),:) - &
                         xo3ws(jindx(j)+1,:) )   
        b2p (j,:) = 0.25*(xo3ws(20-(jindx(j)+1),:) + &
                          xo3ws(jindx(j)+1,:) - &
                          xo3sf(20-(jindx(j)+1),:) - &
                          xo3sf(jindx(j)+1,:) )   
      end do

!--------------------------------------------------------------------


end subroutine obtain_seasonal_data 



!######################################################################

subroutine obtain_randel_data (latb)

!--------------------------------------------------------------------
!    obtain_randel_data reads an input file to obtain the desired 
!    ozone profiles provided as a function of latitude, height and 
!    month of the year. the current input file is for an n30 skyhi
!    grid. no interpolation algorithm is currently available, meaning
!    that this data may only be used when that grid is activated 
!    within FMS.
!--------------------------------------------------------------------

real, dimension(:),   intent(in) :: latb

!--------------------------------------------------------------------
!  intent(in) variables:
!
!       latb      array of model latitudes at cell boundaries [radians]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
!
!       num_level_records     number of data levels in input file
!       num_latitude_records  number of latitude records in data file
!       initial_lat           latitude of first data record in data file
!       delta_lat             latitudinal increment between data records
!       num_time_records      number of time level records in data file
!       qo3_col               array into which ozone data is read
!---------------------------------------------------------------------
      integer, parameter   :: NUM_LEVEL_RECORDS = 40 
      integer, parameter   :: NUM_LAT_RECORDS = 60 
      real                 :: initial_lat = -88.5 
      real                 :: delta_lat = 3.0  
      integer, parameter   :: NUM_TIME_RECORDS = 12 

      real, dimension (NUM_LEVEL_RECORDS) :: qo3_col  

      integer              :: iounit    ! unit number for file reading
      integer              :: nt, k, j, jj  ! various indices
      real                 :: record_lat    ! latitude of current record

!---------------------------------------------------------------------
!    verify that the model and the randel data are on the same
!    latitudinal grid. if not, the data cannot be used without inter-
!    polation, which is currently not supplied.
!---------------------------------------------------------------------
      if (nint(acos(-1.0)/(latb(2) - latb(1))) /= 60) then
        call error_mesg ('ozone_mod' , &
          'randelprofile is on a different latitude grid than &
           &current model. currently no interpolation is  &
	   &available.', FATAL)
      endif

!---------------------------------------------------------------------
!    verify that the model grid and the randel data have the same
!    number of vertical layers. if not, the data cannot be used without
!    interpolation, which is currently not supplied. if they do, send
!    a note asking the user to be certain that the current 40 levels
!    are the same as the randel data's 40 levels. set a flag indicating
!    that the randel data is to be used.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!! THIS DOES NOT PREVENT A 40 level NON-SKYHI GRID FROM USING THIS 
!!! OPTION.  BEWARE !!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!---------------------------------------------------------------------
      if (kmax /= 40) then
        call error_mesg ('ozone_mod', &
          'can only use randelprofile with 40 level &
           &skyhi grid -- no vertical interpolation available.', FATAL)
      else
        call error_mesg ('ozone_mod', &
          'ARE YOU USING THE 40 level SKYHI GRID?? IF NOT,&
           & please stop the run since randel ozone profiles &
           &will be incorrect!!!', NOTE)
      endif

      do_randelprofile = .true.

!---------------------------------------------------------------------
!    allocate an array to hold the randel data set.
!---------------------------------------------------------------------
      allocate (qo3_randel(jdf, kmax, NUM_TIME_RECORDS))

!---------------------------------------------------------------------
!    determine if the randel data input file exists. if so, open the 
!    file and read the data set. save only the data needed by the 
!    current processor in the qo3_randel array. it is assumed that the 
!    jdf latitude rows on the processor are consecutive. close the
!    file upon completion.
!---------------------------------------------------------------------
      if (file_exist ( 'INPUT/randelo3data') ) then
        iounit = open_file ('INPUT/randelo3data', form='formatted',  &
                             action= 'read')
        do nt = 1,NUM_TIME_RECORDS
          jj = 0
          do j=1,NUM_LAT_RECORDS
            read (iounit, FMT = '(5E14.6)')    &
                                   (qo3_col(k),k=1, NUM_LEVEL_RECORDS)
            if (jj /= jdf) then
              record_lat = initial_lat + (j-1) *delta_lat
              if (record_lat >= latb(1)*radian ) then
                jj= jj + 1
                qo3_randel(jj,:,nt) = qo3_col(:)
              endif
            endif
          enddo
        enddo
        call close_file (iounit)

!---------------------------------------------------------------------
!    if file is not present, write an error message.
!---------------------------------------------------------------------
      else
        call error_mesg ( 'ozone_init', &
            'randelo3data file is not present',FATAL)
      endif

!----------------------------------------------------------------------


end subroutine obtain_randel_data




!######################################################################

subroutine obtain_input_file_data 

!---------------------------------------------------------------------
!    obtain_input_file_data reads an input file containing a single
!    column ozone profile.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer   :: iounit    ! unit to read file on
      integer   :: k         ! vertical loop index
      integer   :: kmax_file ! number of levels of data in file


!-------------------------------------------------------------------
!    set a flag indicating that the ozone data source is an input 
!    column.
!-------------------------------------------------------------------
      do_column_input = .true.

!-------------------------------------------------------------------
!    determine if the input data input file exists. if so, verify that 
!    the number of data records in the file matches the number of model
!    levels. 
!---------------------------------------------------------------------
      if (file_exist ( 'INPUT/id1o3') ) then
         iounit = open_file ('INPUT/id1o3', form='formatted',  &
                             action= 'read')
         read (iounit,FMT = '(i4)') kmax_file
         if (kmax_file /= kmax) then
           call error_mesg ('ozone_mod', &
             'size of ozone profile in input file does not match &
	      &model grid.', FATAL)
         endif

!-------------------------------------------------------------------
!    allocate space for the input data. read the data set. close the 
!    file upon completion.
!---------------------------------------------------------------------
         allocate (qqo3(kmax_file) )
         read (iounit,FMT = '(5e18.10)') (qqo3(k),k=1,kmax_file)
         call close_file (iounit)

!---------------------------------------------------------------------
!    if file is not present, write an error message.
!---------------------------------------------------------------------
       else
         call error_mesg ( 'ozone_init', &
              'desired ozone input file is not present',FATAL)
       endif

!----------------------------------------------------------------------


end subroutine obtain_input_file_data 




!######################################################################

subroutine obtain_fms_zonal_ozone_data (season)

!---------------------------------------------------------------------
!    obtain_fms_zonal_ozone_data generates data at the appropriate time
!    from the basic fms_zonal_ozone input data set, allowing the use of
!    annual mean, fixed seasonal, or seasonally-varying ozone distrib-
!    utions.
!---------------------------------------------------------------------

integer, intent(in) :: season

!--------------------------------------------------------------------
!  intent(in) variables:
!
!      season          scalar integer between 0-5, where 1-4 uses fixed
!                      data (1=winter, 2=spring, etc.), season=0 is 
!                      annual mean ozone, and season=5 is seasonally
!                      varying ozone
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables:
! 
!     ro31, ro32       ozone values for the proper time, at 41 standard
!                      levels and 10 standard latitudes (index 1 = 
!                      equator, index 10= pole, 10 deg resolution)
!     duo3n            array of ozone at proper time, at 41 standard 
!                      levels, over 19 global latitudes (10 deg 
!                      resolution), index 1 = north pole, index 10=
!                      equator, index 19 = south pole)
!     data4            array of ozone at proper time, at 81 levels, over
!                      19 global latitudes (10 deg resolution). if 
!                      seasonally-varying ozone is desired, there are
!                      4 such arrays, to which harmonic interpolation
!                      will be applied.
!     ph3              sigma levels at which zonal ozone data set data
!                      exists
!     o3hi             array containing ozone data at 25 high layers in
!                      which mixing ratios are invariant with season
!     o3lo1            array containing ozone data at 16 lower layers
!                      valid for nh spring, latitudinal indices from
!                      equator to pole
!     o3lo2            array containing ozone data at 16 lower layers
!                      valid for nh fall, latitudinal indices from
!                      equator to pole
!     o3lo3            array containing ozone data at 16 lower layers
!                      valid for nh winter, latitudinal indices from
!                      equator to pole
!     o3lo4            array containing ozone data at 16 lower layers
!                      valid for nh summer, latitudinal indices from
!                      equator to pole
!     pref             assumed surface pressure value used to convert 
!                      ph3 from sigma to pressure
!     
!--------------------------------------------------------------------
      real, dimension (10,41)             ::  ro31, ro32
      real, dimension (19,41)             ::  duo3n
      real, dimension (19,81,4)           ::  data4
      real, dimension (82)                ::  ph3
      real, dimension (10,25)             ::  o3hi
      real, dimension (10,16)             :: o3lo1, o3lo2, o3lo3, o3lo4

      real                      ::  pref = 101325. 
      character(len=48)         ::  err_string ! part of error message
      integer                   ::  iounit   ! used to read input data
      integer                   ::  j, k, kk, n   ! various indices

!---------------------------------------------------------------------
!    be sure that the input argument is valid.
!---------------------------------------------------------------------
      if (season < 0 .or. season > 5) then
          write (err_string,9001) season
 9001     format ('invalid value of season=',i10)
          call error_mesg ('obtain_fms_zonal_ozone_data', err_string, FATAL)
      endif

!---------------------------------------------------------------------
!    save input argument as a module variable.
!---------------------------------------------------------------------
      iseason = season

!---------------------------------------------------------------------
!    determine if the zonal ozone input data file exists. if so, 
!    open the file and read the data set.  close the file upon 
!    completion. if it is not present, write an error message.
!---------------------------------------------------------------------
      if (file_exist ( 'INPUT/zonal_ozone_data') ) then
        iounit = open_file ('INPUT/zonal_ozone_data',   &
                            form='unformatted', action='read')
        read (iounit) ph3
        read (iounit)    
        read (iounit) o3hi
        read (iounit) o3lo1
        read (iounit) o3lo2
        read (iounit) o3lo3
        read (iounit) o3lo4
        call close_file (iounit)
      else
        call error_mesg ( 'obtain_fms_zonal_ozone_data', &
                'zonal_ozone_data data file is not present',FATAL)
      endif


!---------------------------------------------------------------------
!    define standard pressure interfaces (ph).
!---------------------------------------------------------------------
      ph(:) = ph3(:)*pref

!-----------------------------------------------------------------------
!    define the seasonally-invariant elements of arrays ro31, ro32 from
!    the values in o3hi.
!---------------------------------------------------------------------
      do k=1,25
         ro31(:,k) = o3hi(:,k)
         ro32(:,k) = o3hi(:,k)
      end do

!---------------------------------------------------------------------
!    define the ozone values at the lower levels to be seasonally
!    varying. northern hemisphere seasons are used to define the
!    indices, with n = 1 being winter, going to n = 4 being fall.
!--------------------------------------------------------------------
      do n=1,4

!---------------------------------------------------------------------
!    for nh spring or fall, obtain the o3lo1 and o3lo2 data (ro31 and 
!    ro32). 
!----------------------------------------------------------------------
        if (n == 2 .or. n == 4) then
          do k=26,41
            ro31(:,k) = o3lo1(:,k-25)
            ro32(:,k) = o3lo2(:,k-25)
          end do
        endif

!---------------------------------------------------------------------
!    for nh winter or summer, obtain the o3lo3 and o3lo4 data (ro31 and
!    ro32). 
!----------------------------------------------------------------------
        if (n == 1 .or. n == 3) then
          do k=26,41
            ro31(:,k) = o3lo3(:,k-25)
            ro32(:,k) = o3lo4(:,k-25)
          end do
        endif

!---------------------------------------------------------------------
!    define ozone values for both hemispheres -- indices 1->9 = nh,
!    index 10 = equator, indices 11-19 = sh. for nh spring (n=2)
!    and nh winter (n=1), nh values are contained in ro31, in
!    reverse latitudinal order. sh values are contained in ro32, going
!    from equator to south pole.
!----------------------------------------------------------------------
        if (n == 2 .or. n == 1) then
          do k=1,41
            do j=1,10
              duo3n(j  ,k) = ro31(11-j,k)
              duo3n(j+9,k) = ro32(j   ,k)
            end do
            duo3n(10 ,k) = 0.50*(ro31(1,k) + ro32(1,k))
          end do

!---------------------------------------------------------------------
!    for nh summer (n=3) and nh fall (n=4), nh values are 
!    contained in ro32, in reverse latitudinal order. sh values are 
!    contained in ro31, going from equator to south pole.
!---------------------------------------------------------------------
        else if(n == 4 .or. n == 3) then
          do k=1,41
            do j=1,10
              duo3n(j  ,k) = ro32(11-j,k)
              duo3n(j+9,k) = ro31(j   ,k)
            end do
            duo3n(10 ,k) = 0.50*(ro31(1,k) + ro32(1,k))
          end do
        endif

!-----------------------------------------------------------------------
!    vertical interp between original data using bessels half-point 
!    interpolation formula.
!-----------------------------------------------------------------------
        do kk=4,78,2
          k = kk/2
          o3data(:,kk,n) = 0.50*(duo3n(:,k) + duo3n(:,k+1)) -  &
                           (duo3n(:,k+2) - duo3n(:,k+1) -   &
                            duo3n(:,k) + duo3n(:,k-1))/16.
        end do
        o3data(:, 2,n) = 0.50*(duo3n(:,2) + duo3n(:,1))
        o3data(:,80,n) = 0.50*(duo3n(:,41) + duo3n(:,40))

!---------------------------------------------------------------------
!    put intermediate (unchanged) data into new array.                
!---------------------------------------------------------------------
        do kk=1,81,2
          k = (kk + 1)/2
          o3data(:,kk,n) = duo3n(:,k)
        end do
      end do  ! n loop

!-----------------------------------------------------------------------
!    prepare seasonal interpolation.
!-----------------------------------------------------------------------
      if (iseason == 5) then
        data4(:,:,1) = 0.25*(o3data(:,:,1) + o3data(:,:,2)  &
                           + o3data(:,:,3) + o3data(:,:,4))
        data4(:,:,2) = 0.50*(o3data(:,:,2) - o3data(:,:,4))
        data4(:,:,3) = 0.50*(o3data(:,:,1) - o3data(:,:,3))
!! BUGFIX  *******************
!! the following line should be replaced:
!!!     data4(:,:,4) = 0.25*(o3data(:,:,1) - o3data(:,:,1)  &
!!!  with:
        data4(:,:,4) = 0.25*(o3data(:,:,1) - o3data(:,:,2)  &
                           + o3data(:,:,3) - o3data(:,:,4))
        o3data(:,:,1) = data4(:,:,1)
        o3data(:,:,2) = data4(:,:,2)
        o3data(:,:,3) = data4(:,:,3)
        o3data(:,:,4) = data4(:,:,4)
      endif

!-----------------------------------------------------------------------
!    prepare annual mean data.
!-----------------------------------------------------------------------
      if (iseason == 0) then
        data4(:,:,1) = 0.25*(o3data(:,:,1) + o3data(:,:,2)  &
                           + o3data(:,:,3) + o3data(:,:,4))
        o3data(:,:,1) = data4(:,:,1)
        iseason = 1
      endif

!---------------------------------------------------------------------



end subroutine obtain_fms_zonal_ozone_data




!######################################################################

subroutine obtain_clim_zonal_ozone_data (lonb, latb)

!----------------------------------------------------------------------
!    obtain_clim_zonal_ozone_data provides the necessary information 
!    to interpolator_mod so that the appropriate clim_ozone data may
!    be obtained later on when needed.
!---------------------------------------------------------------------

real, dimension(:), intent(in) :: lonb, latb

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       lonb      array of model longitudes at cell boundaries [radians]
!       latb      array of model latitudes at cell boundaries [radians]
!
!-----------------------------------------------------------------

      call interpolator_init (o3, "o3.climatology.nc", lonb, latb, &
                              data_out_of_bounds=(/CONSTANT/), &
                              vert_interp=(/INTERP_WEIGHTED_P/) )

end subroutine obtain_clim_zonal_ozone_data


!######################################################################

subroutine obtain_photochem_adjust_data (latb, pref)

!---------------------------------------------------------------------
!    obtain_photochem_adjust_data reads the temperature distribution for
!    which the ozone distribution applies. changes in temperature 
!    from these values will result in a change of ozone amount, sim-
!    ulating the photochemical adjustment of ozone to the ambient 
!    temperature. 
!---------------------------------------------------------------------

real, dimension(:),   intent(in) :: latb
real, dimension(:,:), intent(in) :: pref

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       latb      array of model latitudes at cell boundaries [radians]
!       pref      array containing two reference pressure profiles 
!                 [pascals]
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  local variables:
!
!       num_level_records     number of data levels in input file
!       num_latitude_records  number of latitude records in data file
!       hemi_data             does the data set contain only 
!                             hemispheric data (in contrast to global) ?
!       tradin2               standard temperature distribution con-
!                             sistent with ozone distribution
!       jindx                 latitude index in the data set closest to
!                             but less than the model latitude.
!       displacement          distance of model latitude from data lat-
!                             itude jindx relative to delta latitude of
!                             input data
!
!---------------------------------------------------------------------
      integer, parameter   :: NUM_LEVEL_RECORDS = 20 
      integer, parameter   :: NUM_LAT_RECORDS = 36 
      logical              :: hemi_data = .true.

      real,    dimension     &
             (NUM_LEVEL_RECORDS,NUM_LAT_RECORDS)    :: tradin2
      real,    dimension(size(latb)-1)              :: displacement
      integer, dimension(size(latb)-1)              :: jindx


      integer    ::  unit         ! unit for reading input data 
      integer    :: ll, k, j      ! various indices
      integer    :: nlats         ! number of latitudes of data in input
                                  ! file
      integer    :: nlevs         ! number of levels of data in input 
                                  ! file
      real       :: init_lat      ! lowest latitude of input file data

!-------------------------------------------------------------------
!    define the region of the model grid in which photochemical 
!    adjustment is desired (levels kobeg -> koend). if no part of the
!    grid will be adjusted, set do_ozone_adjustment to .false., write
!    a note, and return.
!-------------------------------------------------------------------
      if (0.5*(pref(1,1) + pref(2,1)) < ozone_adjustment_limit) then
        kobeg = 1
      else
        do_ozone_adjustment = .false.
        call error_mesg ('obtain_photochem_adjust_data', &
            'do_ozone_adjustment was set to .true.; however, no model &
	    &levels are located above the ozone_adjustment_limit, so &
	    &that do_ozone_adjustment is reset to .false.', NOTE)
        return
      endif

!----------------------------------------------------------------------
!    define the lower limit (highest k index) at which adjustment is to
!    be applied.
!----------------------------------------------------------------------
      do k=1,kmax
        if (0.5*(pref(k,1) + pref(k+1,1)) > ozone_adjustment_limit) then
          koend = k - 1
          exit
        endif
      end do

!---------------------------------------------------------------------
!    verify that the model grid and the photochemical data have the same
!    number of vertical layers. if not, the data cannot be used without
!    interpolation, which is currently not supplied. if they do, send
!    a note asking the user to be certain that the current 40 levels
!    are the same as the randel data's 40 levels. set a flag indicating
!    that the randel data is to be used.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!! THIS DOES NOT PREVENT A 40 level NON-SKYHI GRID FROM USING THIS 
!!! OPTION.  BEWARE !!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if (kmax /= 40) then
        call error_mesg ('obtain_photochem_adjust_data', &
           'ozone_adjustment_file has data for top 20 levels of &
           &skyhi 40 level grid -- currently can only be used &
           &with that vertical grid', FATAL)
      else 
        call error_mesg ('obtain_photochem_adjust_data', &
          'ARE YOU USING THE 40 level SKYHI GRID?? IF NOT,&
          & please turn off ozone_adjustment since the input &
          &data file is inconsistent with your grid!!!', NOTE)
      endif

!---------------------------------------------------------------------
!    read the data file containing the reference temperature profile
!    consistent with the ozone distribution. the number of latitudes 
!    and levels of data present are first read.
!---------------------------------------------------------------------
      if (file_exist ( 'INPUT/ozone_adjustment_file') ) then
        unit = open_file('INPUT/ozone_adjustment_file',&
                         form='unformatted', action = 'read')
        read  (unit) nlats,  nlevs

!---------------------------------------------------------------------
!    verify that koend is not larger than the vertical dimension of the 
!    photochemical input data. 
!---------------------------------------------------------------------
        if (koend > nlevs) then
          call error_mesg ('obtain_photochem_adjust_data', &
            ' apparent error -- data required for more levels than is &
	     &present in adjustment input file, suggesting you are not &
	     &running a  40 level skyhi grid', FATAL)
        endif

!---------------------------------------------------------------------
!    allocate space and then read the latitudes at which the input data
!    is present and the reference temperature profiles. close the unit
!    upon completion.
!---------------------------------------------------------------------
        read (unit) init_lat             
        read (unit) tradin2        
        call close_file (unit)

!---------------------------------------------------------------------
!    if file is not present, write an error message.
!---------------------------------------------------------------------
      else
        call error_mesg ( 'obtain_photochem_adjust_data', &
             'desired ozone input file is not present',FATAL)
      endif

!----------------------------------------------------------------------
!    find the location of the model latitudes within the input data
!    array. jindx is the nearest lower index to the model latitude 
!    and displacement is the displacement of the model latitude from
!    that index value, relative to the input data resolution.
!----------------------------------------------------------------------
      call find_nearest_lower_index (latb, hemi_data, init_lat,  &
                                     acos(-1.0)/nlats, jindx,    &
                                     displacement)

!----------------------------------------------------------------------
!    obtain values at model latitudes by linearly interpolating the
!    input data. save in module variable ztrad.
!----------------------------------------------------------------------
      allocate ( ztrad(jdf, KOBEG:KOEND) )
      do j=1,jdf
        ll = jindx(j)
        do k=kobeg,koend
          ztrad(j,k) = tradin2(k,ll) + (tradin2(k,ll+1) -  &
                       tradin2(k,ll))*displacement(j)      
        end do
      end do

!-------------------------------------------------------------------


end subroutine obtain_photochem_adjust_data




!######################################################################

subroutine obtain_ozone_loss_data (latb, pref)

!---------------------------------------------------------------------
!    obtain_ozone_loss_data reads an input file to define the ozone
!    loss to be applied to the input ozone profile.
!---------------------------------------------------------------------

real, dimension(:),   intent(in) :: latb
real, dimension(:,:), intent(in) :: pref

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       latb      array of model latitudes at cell boundaries [radians]
!       pref      array containing two reference pressure profiles 
!                 [pascals]
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  local variables:
!
!       modlats      latitudes at which the data in the file apply
!       top_press    pressure at top of region over which depletion
!                    occurs. 
!       btm_press    pressure at bottom of region over which depletion
!                    occurs. 
!       o3depl_prof  o3 column change in dobson units as a function
!                    of (JD=60) latitudes and 12 months. varies with
!                    season and latitude.
!       nlats        number of latitudes for which data is present in
!                    input file
!       nlevs        number of levels for which data is present in the
!                    input file
!       jind_start   global index of first latitude row on processor
!
!-----------------------------------------------------------------
      real, dimension(:,:), allocatable   :: top_press, o3depl_prof   
      real, dimension(:),   allocatable   :: btm_press, modlats        

      integer   :: nlats, nlevs, jind_start

      integer   ::  unit            !  unit used to read input data
      integer   ::  k, nn, nseason  !  various indices


!-----------------------------------------------------------------
!    determine the ozone depletion source and open the appropriate
!    input files. if file is not present, write an error message. if
!    an improper depletion source is indicated, write an error message.
!-----------------------------------------------------------------
      if (trim(o3monloss_type) == 'tomsprofile1') then
        do_tomsprofile1 = .true.
        if (file_exist ( 'INPUT/ozone_depletion_toms1') ) then
          unit = open_file('INPUT/ozone_depletion_toms1',&
                           form='unformatted', action = 'read')
        else
          call error_mesg ( 'obtain_ozone_loss_data', &
                  'desired ozone input file is not present',FATAL)
        endif
      else if (trim(o3monloss_type) == 'tomsprofile2') then
        do_tomsprofile2 = .true.
        if (file_exist ( 'INPUT/ozone_depletion_toms2') ) then
          unit = open_file('INPUT/ozone_depletion_toms2',&
                           form='unformatted', action = 'read')
        else
          call error_mesg ( 'obtain_ozone_loss_data', &
                'desired ozone input file is not present',FATAL)
        endif
      else
        call error_mesg ( 'obtain_ozone_loss_data',    &
                 ' o3monloss_type is not an acceptable value.', FATAL)
      endif   

!--------------------------------------------------------------------
!    read the number of latitudes present in the depletion data. if
!    this differs from the model resolution, write an error message
!    and stop, since latitudinal interpolation is not currently 
!    available.
!--------------------------------------------------------------------
      read  (unit) nlats
      if (nlats /= nint(acos(-1.0)/(latb(2) - latb(1))) ) then
        call error_mesg ('obtain_ozone_loss_data' , &
              ' tomsprofile input file is on a different latitude &
	       & grid than is being used. currently no interpolation &
	       &is available.', FATAL)
      endif

!--------------------------------------------------------------------
!    allocate arrays for and read the input data latitudes, the top
!    pressure level where depletion is to apply, the bottom pressure
!    level where depletion is to apply, and the magnitudes of column
!    depletion at these latitudes for each month of the year. upon
!    completion, close the input file.
!---------------------------------------------------------------------
      allocate (modlats (nlats))
      allocate (top_press(nlats,4))
      allocate (btm_press(nlats))
      allocate (o3depl_prof (nlats,12))

      read  (unit) modlats
      read  (unit) top_press      
      read  (unit) btm_press     
      read  (unit) o3depl_prof  

      call close_file (unit)

!---------------------------------------------------------------------
!    define the starting latitude index from this global data set 
!    needed by the current processor. assumption is made that processor
!    owns jdf consecutive latitude rows beginning at jind_start.
!---------------------------------------------------------------------
      do nn=1,nlats
        if (latb(1) .lt. modlats(nn) .and.   &
            latb(2) .gt. modlats(nn)) then
          jind_start = nn
          exit
        endif
      end do

!-----------------------------------------------------------------
!    allocate arrays to hold the ozone depletion and the pressure 
!    limits of region to be depleted. define the depletion amount and
!    the vertical k index region over which it is to be applied.
!    the depletion differs for each month and the top pressure varies
!    seasonally.
!-----------------------------------------------------------------
      allocate ( o3depb(jdf)    )
      allocate ( o3dept(jdf,4)  )
      allocate ( o3depl(jdf,12) )
      do nn=1,jdf
        o3depl(nn,:) = o3depl_prof (nn+jind_start-1,:)
        do k=2, kmax
          if (btm_press(nn+jind_start-1) <     &
	         0.5*(pref(k,1) + pref(k+1,1)) ) then
            o3depb (nn) = k-1
            exit
          endif
        end do
        do nseason=1,4
          do k=2, kmax
            if (top_press(nn+jind_start-1,nseason) <    &
                    0.5*(pref(k,1) + pref(k+1,1)) ) then
              o3dept (nn,nseason) = k-1
              exit
            endif
          end do
         end do
      end do  ! (nn loop)

!---------------------------------------------------------------------
!    deallocate local arrays.
!---------------------------------------------------------------------
      deallocate (o3depl_prof)
      deallocate (modlats)
      deallocate (top_press)
      deallocate (btm_press)


end subroutine obtain_ozone_loss_data



!#####################################################################

subroutine find_nearest_lower_index (latb, hemi_data, inilat, dellat,&
                                     jindx, displa)

!---------------------------------------------------------------------
!    find_nearest_lower_index finds the location of the model latitudes
!    within the input data set. it returns the index closest to but
!    smaller than the model latitude (jindx) and the fractional distance
!    of the model latitude between jindx and jindx + 1 (displa).
!---------------------------------------------------------------------
 
real,    dimension(:), intent(in)   :: latb
logical,               intent(in)   :: hemi_data
real,                  intent(in)   :: inilat, dellat
integer, dimension(:), intent(out)  :: jindx
real   , dimension(:), intent(out)  :: displa

!---------------------------------------------------------------------
!  intent(in) variables:
!
!       latb        array of model latitudes at cell boundaries 
!                   [ radians ]
!       hemi_data   is the input data set hemispheric (as opposed to 
!                   global)?
!       inilat      lowest latitude at which data is present in the 
!                   input data set
!       dellat      latitudinal resolution of input data set
!
!  intent(out) variables:
! 
!       jindx       latitude index in the data set closest to but less 
!                   than the model latitude.
!       displa      distance of model latitude from data latitude jindx 
!                   relative to delta latitude of input data
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  local variables:
!
!       lat         latitudes at model grid points [ radians ]
!       ref_lat     data set latitude closest to but less than the
!                   model grid point [ radians ]
!       final_lat   highest latitude in input data set [ radians ]
!       npts        number of latitudes in input data set
!
!------------------------------------------------------------------
      real, dimension(size(latb,1)-1) :: lat
      real                            :: ref_lat, final_lat
      integer                         :: npts

      integer                         ::  j   ! do-loop index  

!---------------------------------------------------------------------
!    define the number of latitiude points in the data set, based on
!    the initial latitude provided and the data set latitudinal reso-
!    lution.
!---------------------------------------------------------------------
      npts = nint((0.4999*acos(-1.0) - inilat)/dellat) + 1
      if (hemi_data) then
        npts = npts/2
      endif

!---------------------------------------------------------------------
!    define the latitudes of the model grid points. if the input data
!    set is hemispheric, set the model latitudes to be the colatitudes
!    in the northern hemisphere.
!---------------------------------------------------------------------
      do j = 1,jdf
        lat(j) = 0.5*(latb(j) + latb(j+1))
        if (hemi_data) then
          if (lat(j) > 0.0) then
            lat(j) = -1.0*lat(j)
          endif
        endif

!--------------------------------------------------------------------
!    define the index of the data set latitude closest to but less than
!    the model latitude. handle cases where the input data set does not
!    extend from pole to pole by extrapolating interior gradients 
!    outward. define the fractional displacement of the model grid 
!    point from this reference index.
!--------------------------------------------------------------------
        if (lat(j) > inilat) then
          jindx(j) = int((lat(j) - inilat)/dellat)+1
          final_lat = inilat + (npts-1)*dellat
          if (hemi_data .and. lat(j) > final_lat) then
            jindx(j) = npts - 1
            displa(j) = (lat(j) - (final_lat - dellat))/dellat
          else
            ref_lat = inilat + (jindx(j)-1)*dellat
            displa(j) = (lat(j) - ref_lat)/dellat
          endif
        else
          jindx(j) = 1
          ref_lat = inilat + (jindx(j)-1)*dellat
          displa(j) = (lat(j) - ref_lat)/dellat
        endif
      end do

!------------------------------------------------------------------



end subroutine find_nearest_lower_index



!##################################################################

subroutine obtain_current_ozone (js, je, lat, Time, zdoblos, zdduo3n)

!---------------------------------------------------------------------
!    obtain_current_ozone defines the ozone field  (zdduo3n) and, if
!    desired, the ozone depletion profile (zdoblos) at the specified
!    time.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
integer,              intent(in)   :: js, je
real, dimension(:,:), intent(in)   :: lat
type(time_type),      intent(in)   :: Time
real, dimension(:),   intent(out)  :: zdoblos
real, dimension(:,:), intent(out)  :: zdduo3n
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      js,je        starting/ending subdomain j indices of data in 
!                   the physics_window being integrated
!      lat          latitude of model points  [ radians ]
!      Time         current model time [ time_type (days, seconds) ] 
!
!   intent(out) variables:  
!
!      zdoblos      total column ozone loss for each column, interp-
!                   olated to the specified time [ dobson units ]
!      zdduo3n      ozone distribution at the specified time 
!                   mass mixing ratio [ g / g ]
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:
!
!      dduo3n       harmonically interpolated, seasonally varying    
!                   ozone field at the speciied time at the data
!                   set latitude closest to but less than the model
!                   latitude [ 10**4 g / g ]
!      dduo3np      harmonically interpolated, seasonally varying    
!                   ozone field at the specified time at the data
!                   set latitude closest to but greater than the model
!                   latitude [ 10**4 g / g ]
!      o3monf       fraction of month which has passed; used for monthly
!                   interpolation [ dimensionless ]
!      o3lag        displacement into calendar year at which spring -
!                   fall data is valid for seasonally-varying ozone 
!                   case (april 13 or 14) [ days ]
!      arg          displacement in year of current time from start of 
!                   ozone year (time at which winter-summer distrib-
!                   utions in seasonal ozone case are valid, ~ 
!                   january 12 / 13) [ radians ]
!      cosarg, sinarg, cos2ar
!                   components used in the harmonic interpolation of 
!                   the seasonal ozone data
!      im, imp, im1
!                   indices indicating the month entry into monthly
!                   data sets, used for time interpolation of these
!                   fields
!
!---------------------------------------------------------------------
      real, dimension (size(zdduo3n,1), size(zdduo3n,2)) :: &
                                                    dduo3n, dduo3np
      real               ::   o3monf
      real               ::   o3lag=104.0
      real               ::   arg, cosarg, sinarg, cos2ar
      integer            ::   im, imp, im1

      integer            :: j, k   ! do-loop indices

!---------------------------------------------------------------------
!    determine the ozone fields valid at the specified time (Time).
!---------------------------------------------------------------------

!-----------------------------------------------------------------
!    if ozone distribution uses monthly input data, call 
!    interpolate_monthly_data to calculate index for linear inter-
!    polation into monthly ozone data fields.
!-------------------------------------------------------------------
        if (o3monloss .or. do_randelprofile) then
          call interpolate_monthly_data (Time, o3monf, im, im1, imp)
        endif  

!--------------------------------------------------------------------
!     if the ozone field is to be seasonally varying, determine current 
!     time with respect to "ozone year". determine argument for harmonic
!     interpolation of seasonal ozone data. do an harmonic time inter-
!     polation of the seasonal (xo3sf, xo3ws) ozone data to the spec-
!     ified time. 
!--------------------------------------------------------------------
        if (do_seasonal_ozone) then
          call harmonic_seasonal_interpolation (o3lag, Time, arg)
          cosarg = COS(arg)
          sinarg = SIN(arg)
          cos2ar = COS(2.0E+00*arg)
          do k=kmin,kmax
            do j=js,je
              dduo3n(j-js+1,k) = xbar(j,k) + a1(j,k)*sinarg +  &
                                 b1(j,k)*cosarg + b2(j,k)*cos2ar
              dduo3np(j-js+1,k) = xbarp(j,k) + a1p(j,k)*sinarg +   &
                                  b1p(j,k)*cosarg + b2p(j,k)*cos2ar
            end do
          end do

!----------------------------------------------------------------------
!    convert this time-interpolated ozone field to model latitudes and 
!    levels (zdduo3n). convert units to mass mixing ratio, g/g.
!-----------------------------------------------------------------------
          do k=kmin,kmax
            do j=1,je-js+1
              zdduo3n(j,k) = 1.0E-04*(dduo3n(j,k) +   &
                             displacement(j+js-1)*&
                             (dduo3np(j,k) - dduo3n(j,k)))
            end do
          end do

!-----------------------------------------------------------------------
!     for the annual mean ozone case, retrieve the proper data for this
!     latitude and convert the ozone mixing ratio to g/g.
!-----------------------------------------------------------------------
        else if (do_annual_mean_ozone) then
          do k=kmin,kmax
            do j=1,je-js+1
              zdduo3n(j,k) = 1.0E-04*xo3ann(j+js-1,k)               
            end do
          end do

!-----------------------------------------------------------------------
!     time interpolate the monthly randel profile data. values are in
!     g/g.
!----------------------------------------------------------------------
        else if (do_randelprofile) then
          do k=kmin,kmax
            do j = 1,je-js+1
              zdduo3n(j,k) = qo3_randel(j+js-1,k,im1) +   &
                             (qo3_randel(j+js-1,k,imp) -   &
                              qo3_randel(j+js-1,k,im))*o3monf
            end do
          end do

!-----------------------------------------------------------------------
!     use single column ozone data read from input file.
!----------------------------------------------------------------------
        else if (do_column_input) then
          do k=kmin,kmax
            do j = 1,je-js+1
              zdduo3n(j,k) = qqo3(k)
            end do
          end do
        endif

!---------------------------------------------------------------------
!    obtain the time-interpolated ozone depletion in each column.
!---------------------------------------------------------------------
        if (o3monloss) then
          do j = 1,je-js+1
            zdoblos(j) = o3depl (js+j-1,im1) + (o3depl (j+js-1,imp) -  &
                         o3depl (j+js-1,im))*o3monf
          end do
        endif  

!-----------------------------------------------------------------------
!    for standalone case, use single column ozone data read from input
!    file.
!----------------------------------------------------------------------
      if (do_column_input) then
        do k=kmin,kmax
          do j=1,je-js+1
            zdduo3n(j,k) = qqo3(k)
          end do
        end do
      endif 

!-------------------------------------------------------------------



end subroutine obtain_current_ozone



!#####################################################################

subroutine harmonic_seasonal_interpolation (lag, Time, arg)

!----------------------------------------------------------------------
!    harmonic_seasonal_interpolation returns the specified time in
!    terms of an angular displacement into the ozone year.
!----------------------------------------------------------------------

real,            intent(in)              :: lag
type(time_type), intent(in)              :: Time
real,            intent(out)             :: arg

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      lag          displacement into calendar year at which spring -
!                   fall data is valid for seasonally-varying ozone 
!                   case (april 13 or 14) [ days ]
!      Time         date and time to which data is to be interpolated
!                   [ time_type (days, seconds) ]
! 
!   intent(out) variables:
!
!      arg          displacement in year of current time from start of 
!                   ozone year (time at which winter-summer distrib-
!                   utions in seasonal ozone case are valid, ~ 
!                   january 12 / 13) [ radians ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:
!
!      daypmn       number of days in each month, adjusted for leap 
!                   years
!      xday         displacement of current time from start of ozone
!                   year (time at which winter-summer distributions
!                   in seasonal ozone case are valid, ~ january 12 / 13)
!                   [ days ]
!      ndiy, fndiy 
!                   number of days in the current year
!      year, month, day, hour, minute, second
!                   current time components
!      po4          one-fourth of current year [ days ]
!      lag2         value of lag, adjusted if needed for leap year
!                   [ days ]
!      rday         specified time, measured from start of calendar
!                   year [ days ]
!
!---------------------------------------------------------------------
      real, dimension(12) ::  daypmn  
      data daypmn/        &
           31.0E+00, 28.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, &
           31.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, 31.0E+00 /


      real                ::  xday
      integer             ::  ndiy
      real                ::  fndiy
      integer             ::  year, month, day, &
                              hour, minute, second
      real                ::  po4, lag2, rday

      integer             ::  m  ! loop index

!-------------------------------------------------------------------
!    define day and month structure of current year.
!---------------------------------------------------------------------
      ndiy = days_in_year (Time)
      fndiy = ndiy
      po4 = 0.25E+00*ndiy
      if (ndiy == 366) daypmn(2) = 29.0

!--------------------------------------------------------------------
!    in leap years, define the valid time of the spring / fall obs to
!    be one day later in the year, accounting for february 29.
!---------------------------------------------------------------------
      if (ndiy == 366.0) then
        lag2 = lag + 1.0
      else
        lag2 = lag
      endif

!---------------------------------------------------------------------
!    retrieve the date and time components from the time_type specified
!    time.
!---------------------------------------------------------------------
      call get_date (Time, year, month, day, hour, minute, second)

!---------------------------------------------------------------------
!    convert the time components to be in units of days from the 
!    beginning of the calendar year.
!---------------------------------------------------------------------
      do m=1, month-1
        day = day + daypmn(m) 
      end do
      rday = day + FLOAT(hour)/24.0 + FLOAT(minute)/1440.0 + &
             FLOAT(second)/86400.0

!---------------------------------------------------------------------
!    define the time in terms of days from the beginning of the
!    ozone year, i.e., the valid time of the w/s seasonal ozone obs
!    (~ january 12). convert this to an angular measure.
!---------------------------------------------------------------------
      xday = MOD(rday + (po4 - lag2) + ndiy, fndiy)
      arg = 2.0E+00*acos(-1.0)*xday/ndiy

!---------------------------------------------------------------------


end subroutine harmonic_seasonal_interpolation 



!####################################################################

subroutine interpolate_monthly_data (Time, o3monf, im, im1, imp)

!---------------------------------------------------------------------
!    interpolate_monthly_data defines the indices and fractional part
!    of the month to use when interpolating to the specified time 
!    within a monthly data set.
!----------------------------------------------------------------------

type(time_type), intent(in)    :: Time
real,            intent(out)   :: o3monf
integer,         intent(out)   :: im, im1, imp

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      Time         specified time at which data is desired
!                   [ time_type (days, seconds) ]
!
!   intent(out) variables:
!
!      o3monf       fractional displacement of specified time from
!                   previous 15th of the month [ dimensionless ]
!      im           index of current month [ dimensionless ]
!      im1          index of previous 15th of the month 
!                   [ dimensionless ]
!      imp          index of upcoming 15th of the month 
!                   [ dimensionless ]
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:
! 
!      daypmn       number of days in each month, adjusted for leap 
!                   years
!      delday       number of days from 15th of month to 15th of 
!                   next month
!      ndiy         number of days in the current year
!      year, month, day, hour, minute, second
!                   current time components
!      deld         the number of days between the 15th of this month
!                   and the 15th of next month.
!      dpm          number of days in current month
!
!---------------------------------------------------------------------

      real, dimension(12)   ::  daypmn, delday 
      data daypmn/        &
         31.0E+00, 28.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, &
         31.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, 31.0E+00 /

      data delday/   &
         31.0E+00, 28.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, &
         31.0E+00, 31.0E+00, 30.0E+00, 31.0E+00, 30.0E+00, 31.0E+00 /

      integer       ::  ndiy, year, month, day, hour, minute, second
      real          ::  deld, dpm

 
!---------------------------------------------------------------------
!    define number of days in current year.
!---------------------------------------------------------------------
      ndiy = days_in_year (Time)

!---------------------------------------------------------------------
!    convert specified time to its components.
!---------------------------------------------------------------------
      call get_date (Time, year, im, day, hour, minute, second)

!---------------------------------------------------------------------
!    if the date is past the 15th of the month, find the fraction of
!    the period between the 15th of this month and the 15th of next
!    month which has elapsed, resolved only to the nearest hour 
!    (o3monf). account for leap year when necessary. also define
!    the month indices for these two months to be used in interpolating
!    (im, imp), in addition to the calendar month which will be the
!    base of the interpolation (im1).
!---------------------------------------------------------------------
      if (day .GE. 15) then
        imp = MOD(im,12) + 1
        im1 = im
        if (ndiy .EQ. 366 .AND. im .EQ. 2) then
          deld = delday(2) + 1
        else
          deld = delday(im)
        endif
        o3monf = (day + hour/24.0E+00 - 15.0E+00)/deld

!---------------------------------------------------------------------
!    if the date precedes the 15th of the month, find the fraction of
!    the month between the 15th of last month and the 15th of this
!    month which has elapsed (o3monf). account for leap year when nec-
!    essary.  also define the month indices for these two months to be 
!    used in interpolating (im, imp), in addition to the calendar month 
!   which will be the base of the interpolation (im1).
!---------------------------------------------------------------------
      else
        if (im .EQ. 1) then
          imp = 12
        else
          imp = im - 1
        endif
        im1 = imp
        if (ndiy .EQ. 366 .AND. imp .EQ. 2) then
          deld = delday(2) + 1
          dpm  = daypmn(2) + 1
        else
          deld = delday(imp)
          dpm  = daypmn(imp)
        endif
        o3monf = -(dpm + day + hour/24.0E+00 - 15.0E+00)/deld
      endif
  
!---------------------------------------------------------------------- 


end subroutine interpolate_monthly_data


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    INTERFACE FMS_ZONAL_OZONE
!
!
! call fms_zonal_ozone (Rad_time, lat, phaf, ozone)
!
!  separate routines exist within this interface for 3d and 2d
!  array input and output:
!
!  real, dimension(:,:),    intent(in)  :: lat
!  real, dimension(:,:,:),  intent(in)  :: phalf
!  real, dimension(:,:,:),  intent(out) :: ozone
! OR
!  real, dimension(:),    intent(in)    :: lat
!  real, dimension(:,:),  intent(in)    :: phalf
!  real, dimension(:,:),  intent(out)   :: ozone
!
!--------------------------------------------------------------------
!
!  intent(in) variables:
!
!      Time         current model time [ time_type (days, seconds) ] 
!      lat          latitude of model points  [ radians ]
!      phalf        pressure at model layer interfaces [ kg / (m s^2) ]
!
!  intent(out) variables:
!
!      ozone        ozone mass mixing ratio at model levels [ g / g ]
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#######################################################################

subroutine geto3_3d (Time, lat, phalf, ozone)

!---------------------------------------------------------------------
!    geto3_3d retrieves an (i,j,k) array of ozone mass mixing ratio 
!    (g / g ) valid at the specified time to be returned to the 
!    calling routine.
!---------------------------------------------------------------------

type(time_type),         intent(in)  :: Time
real, dimension(:,:),    intent(in)  :: lat
real, dimension(:,:,:),  intent(in)  :: phalf
real, dimension(:,:,:),  intent(out) :: ozone

!---------------------------------------------------------------------
!  local variables:
!
!      rlag         time lag of  valid time of seasonal ozone data from 
!                   start of calendar year [ years ]
!      profile      ozone profile in model grid column
!      dp           pressure difference between data layer interfaces
!      dp1, dp2, dp3, dp4          
!                   pressure differences used to assign data set layer
!                   ozone to the proper model layer
!      o3col        total ozone in a model pressure layer
!      rang         angular position of specified time from start of
!                   ozone year
!      rsin1        argument for harmonic interpolation
!      rcos1        argument for harmonic interpolation
!      rcos2        argument for harmonic interpolation
!      phd          model latitude, guaranteed to be between -90 and +90
!                   degrees
!      th           relative latitudinal distance of model grid point
!                   from nearest lower data set index, normalized by
!                   data set resolution 
!      fyear        fraction of the year which has passed at the 
!                   specified time
!      j1           nearest lower latitude index in data set to the
!                   model latitude
!      j2           j1 + 1; data from index j1 and j2 will be inter-
!                   polated to the model grid point
!
!------------------------------------------------------------------
      real                :: rlag
      real, dimension(81) :: profile, dp, dp1, dp2, dp3, dp4, o3col
      real                :: rang, rsin1, rcos1, rcos2, phd, th, fyear
      integer             :: j1,j2

      integer             :: i,j,k,l  ! various indices

!--------------------------------------------------------------------
!    when seasonally-varying ozone is desired, perform a harmonic time 
!    interpolation to obtain values at the specified time.
!--------------------------------------------------------------------
      if (iseason == 5) then

      if(do_mcm_o3_clim) then
        rlag = 12.6875/365.
      else
        rlag = 1./24.
      endif

!        TK NOTE:  I tried supersource version of rlag, which is
!              rlag = 12.6875/365.0   Source ???
!              This differs from fms rlag, which is 1./24.
!              This did not fix the discrepency in ozone
!              at the sample point.  It made it very slightly worse.
!              7/11/01

        fyear = fraction_of_year (time)
        if (fyear /= current_fyear) then
             rang = 4.0*acos(0.0)*(fyear-rlag)
             rsin1 = sin(rang)
             rcos1 = cos(rang)
             rcos2 = cos(2.0*rang)
             rstd(:,:) =       o3data(:,:,1) + rsin1*o3data(:,:,2)  &
                       + rcos1*o3data(:,:,3) + rcos2*o3data(:,:,4)
             current_fyear = fyear
        endif
!---------------------------------------------------------------------
!    otherwise, no interpolation is needed, use the specified seasonal 
!    data. if annual mean value has been specified, that data will be
!    at iseason = 1.
!---------------------------------------------------------------------
      else
        rstd(:,:) = o3data(:,:,iseason)
      endif

!---------------------------------------------------------------------
!    define the pressure increments of the standard data levels.
!---------------------------------------------------------------------
      do l=1,81
         dp (l) = ph(l+1) - ph(l)
      end do

!---------------------------------------------------------------------
!    perform horizontal interpolation into the data set. profile is
!    the vertical ozone data at the grid point.
!---------------------------------------------------------------------
      do j=1,size(lat,2)
        do i=1,size(lat,1)
          phd = max(min(lat(i,j)*radian, 90.), -90.)
          j1 = 10.000 - phd*0.10
          j1 = max(min(j1, 18), 1)
          j2 = j1 + 1
          th = (10-j1) - phd*0.10
          profile(:) = rstd(j1,:) + th*(rstd(j2,:) - rstd(j1,:))

!---------------------------------------------------------------------
!    now interpolate in the vertical to produce values at model
!    levels. calculate layer-mean ozone mixing ratio from data set for 
!    each model layer.
!---------------------------------------------------------------------
          do k=1,size(ozone,3)

!---------------------------------------------------------------------
!    loop to calculate sums over data set layers to get model layer 
!    ozone mean (o3col).
!---------------------------------------------------------------------

            if(do_mcm_o3_clim) then

              ozone(i,j,k) = 0.0
              do l = 1, 81
                 if ((ph(l+1) .ge. phalf(i,j,k)) .and. (ph(l) .le. phalf(i,j,k+1))) then
                    if ((ph(l+1) .lt. phalf(i,j,k+1)) .and. (ph(l) .lt. phalf(i,j,k)))   &
                      &  ozone(i,j,k) = ozone(i,j,k) + profile(l) * (ph(l+1)-phalf(i,j,k))

                    if ((ph(l+1) .lt. phalf(i,j,k+1)) .and. (ph(l) .ge. phalf(i,j,k)))   &
                      &  ozone(i,j,k) = ozone(i,j,k) + profile(l) * (ph(l+1)-ph(l))

                    if ((ph(l+1) .gt. phalf(i,j,k+1)) .and. (ph(l) .gt. phalf(i,j,k)))     &
                      &  ozone(i,j,k) = ozone(i,j,k) + profile(l) * (phalf(i,j,k+1)-ph(l))
                 endif
              enddo

            else if(.not.do_mcm_o3_clim) then
              do l=1,81
                dp1(l) = ph(l+1) - phalf(i,j,k)
                dp2(l) = phalf(i,j,k+1) - ph(l)
                dp3(l) = ph(l+1) - phalf(i,j,k+1)
                dp4(l) = ph(l) - phalf(i,j,k)
              end do
              where (dp1(:) < 0.0) dp1(:)=0.0
              where (dp2(:) < 0.0) dp2(:)=0.0
              do l=1,81
                o3col(l) = 0.0
                if ( dp3(l) < 0.0 ) then
                  if ( dp4(l) < 0.0 ) then
                    o3col(l) = profile(l)*dp1(l)
                  else
                    o3col(l) = profile(l)*dp(l)
                  endif
                else
                  if ( dp4(l) < 0.0 ) then
                    o3col(l) = profile(l)*(phalf(i,j,k+1) - phalf(i,j,k))
                  else
                    o3col(l) = profile(l)*dp2(l)
                  endif
                endif
              end do

!---------------------------------------------------------------------
!    sum the contributions from each data set layer, then normalize
!    by the model pressure depth to produce a mass mixing ratio.
!---------------------------------------------------------------------
              ozone(i,j,k) = sum(o3col(:))
            endif ! do_mcm_o3_clim if block

            ozone(i,j,k) = ozone(i,j,k)/(phalf(i,j,k+1)-phalf(i,j,k))

!---------------------------------------------------------------------
!    code to cover case when surface pressure is greater than pref.
!---------------------------------------------------------------------
            if (.not.do_mcm_o3_clim .and. ph(82) < phalf(i,j,k+1)) then
              ozone(i,j,k) = profile(81)
            endif

            if(do_mcm_o3_clim) then

!             code to cover case when model resolution is so fine that no value
!             of ph(l) in the ozone data array falls betwen phalf(k+1) and
!             phalf(k).   procedure is to simply grab the nearest value from
!             rdata (or profile in the fms code).

              if (ozone(i,j,k) <= 0.0) then
                do l = 1, 81
                  if(ph(l) < phalf(i,j,k) .and. ph(l+1) > phalf(i,j,k+1) ) then
                    ozone(i,j,k) = profile(l)
                  endif
                enddo
              endif
            endif ! do_mcm_o3_clim if block

          end do ! k loop
        end do   ! i loop
      end do     ! j loop

!---------------------------------------------------------------------
!    convert units from micrograms/gram to kg/kg.
!---------------------------------------------------------------------
      ozone(:,:,:) = ozone(:,:,:)*1.e-6

!-----------------------------------------------------------------------



end subroutine geto3_3d




!#######################################################################

subroutine geto3_2d (Time, lat, phalf, ozone)

!---------------------------------------------------------------------
!    geto3_2d retrieves an (i,k) array of ozone mass mixing ratio 
!    (g / g ) valid at the specified time to be returned to the 
!    calling routine.
!---------------------------------------------------------------------

type(time_type),         intent(in)  :: Time
real, dimension(:),    intent(in)  :: lat
real, dimension(:,:),  intent(in)  :: phalf
real, dimension(:,:),  intent(out) :: ozone

!---------------------------------------------------------------------
!  local variables:
!
!      lat3         2d equivalent of lat
!      phalf3       3d equivalent of phalf
!      ozone3       3d equivalent of ozone
!
!--------------------------------------------------------------------

      real,dimension (size(lat,1),1)                  :: lat3
      real,dimension (size(phalf,1),1, size(phalf,2)) :: phalf3
      real,dimension (size(ozone,1),1, size(ozone,2)) :: ozone3

!---------------------------------------------------------------------
!    add an extra dummy dimension to the input variables.
!---------------------------------------------------------------------
      lat3(:,1)     = lat(:)
      phalf3(:,1,:) = phalf(:,:)

!---------------------------------------------------------------------
!    call the 3d interface of this procedure.
!---------------------------------------------------------------------
      call geto3_3d (time, lat3, phalf3, ozone3)

!---------------------------------------------------------------------
!    remove the extra dummy dimension from the output variables.
!---------------------------------------------------------------------
      ozone(:,:) = ozone3(:,1,:)



end subroutine geto3_2d


 
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                END INTERFACE FMS_ZONAL_OZONE
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!#######################################################################

subroutine get_clim_ozone (model_time, p_half, model_data, is, js)

!--------------------------------------------------------------------
!    get_clim_ozone retrieves the clim_ozone field at the desired place 
!    and time from the o3.climatology.nc file by accessing 
!    interpolator_mod.
!---------------------------------------------------------------------

type(time_type),        intent(in)           :: model_time
real, dimension(:,:,:), intent(in)           :: p_half
real, dimension(:,:,:), intent(out)          :: model_data
integer,                intent(in), optional :: is,js

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      model_time   time at which the climatologically-determined,
!                   time-varying ozone field should apply
!                   [ time_type (days, seconds) ] 
!      phalf        pressure at model interface levels
!                   [ N / m**2 ]
!
!   intent(out) variables:
!
!      model_data   output field containing ozone field at desired time
!                   [ g / g ]
!
!   intent(in), optional variables:
!
!      is,js        starting subdomain i,j indices of data in 
!                   the physics_window being integrated
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    check for module initialization.
!---------------------------------------------------------------------
      if ( .not. ozone_initialized ) &
        call error_mesg ('ozone_mod',  &         
           'module has not been initialized', FATAL)

!--------------------------------------------------------------------
!    obtain appropriate ozone data.
!--------------------------------------------------------------------
      call interpolator (o3, model_time, p_half, model_data, "ozone", &
                         is, js)



end subroutine get_clim_ozone



!#####################################################################



                       end module ozone_mod



