                 module aerosolrad_package_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="">
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    aerosolrad_package_mod provides the radiative properties 
!    associated with the atmospheric aerosols.
! </OVERVIEW>
! <DESCRIPTION>
!    aerosolrad_package_mod provides the radiative properties 
!    associated with the atmospheric aerosols.
! </DESCRIPTION>
!    shared modules:

use fms_mod,               only: open_namelist_file, fms_init, &
                                 mpp_pe, mpp_root_pe, stdlog, &
                                 file_exist, write_version_number, &
                                 check_nml_error, error_mesg, &
                                 FATAL, NOTE, WARNING, close_file
use mpp_io_mod,            only: mpp_open, mpp_close, MPP_RDONLY,   &
                                 MPP_ASCII, MPP_SEQUENTIAL, MPP_MULTI, &
                                 MPP_SINGLE, mpp_io_init

! shared radiation package modules:
                                
use rad_utilities_mod,     only: shortwave_control_type, Sw_control, &
                                 longwave_control_type, Lw_control, &
                                 radiation_control_type, Rad_control,&
                                 aerosol_type, aerosol_properties_type,&
                                 longwave_parameter_type,  &
                                 Lw_parameters, rad_utilities_init, &
                                 solar_spectrum_type,  &
                                 thickavg, thinavg
use esfsw_parameters_mod,  only: Solar_spect, esfsw_parameters_init 
use longwave_params_mod,   only: NBLW, longwave_params_init

!-------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    aerosolrad_package_mod provides the radiative properties 
!    associated with the atmospheric aerosols.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id: aerosolrad_package.F90,v 10.0 2003/10/24 22:00:39 fms Exp $'
character(len=128)  :: tagname =  '$Name: jakarta $'


!---------------------------------------------------------------------
!-------  interfaces --------

public           &
       aerosolrad_package_init, aerosol_radiative_properties, &
       aerosolrad_package_end,  get_aerosol_optical_info, &
       get_aerosol_optical_index
     
private          &

!  called from aerosolrad_package_init: 
   assign_aerosol_opt_props, read_optical_input_file, &
   sw_aerosol_interaction,  lw_aerosol_interaction
       

!---------------------------------------------------------------------
!-------- namelist  ---------

integer, parameter   ::        &
             MAX_OPTICAL_FIELDS = 100  ! maximum number of aerosol 
                                       ! optical property types
logical              ::        &
             do_lwaerosol = .false.    ! aerosol efects included in lw
                                       ! radiation ?
logical              ::        &
             do_swaerosol = .false.    ! aerosol effects included in sw
                                       ! radiation ?
character(len=48)    ::        &
             aerosol_data_set = ' '    ! source of aerosol data; if 
                                       ! aerosols not desired remains
                                       ! ' ', otherwise is set to either
                                       ! 'shettle-fenn' or 'predicted'
                                       ! (predicted not yet implemented)
character(len=64)    ::        &
             aerosol_optical_names(MAX_OPTICAL_FIELDS) = '  '
                                       ! names associated with the 
                                       ! optical property types that
                                       ! are to be used in this 
                                       ! experiment
character(len=64)    ::        &
             optical_filename = ' '    ! name of file containing the
                                       ! aerosol optical property types


namelist / aerosolrad_package_nml /                          &
                                    do_lwaerosol, do_swaerosol, &
                                    aerosol_data_set, &
                                    aerosol_optical_names, &
                                    optical_filename   

!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


!---------------------------------------------------------------------
!    the following variables define the number and type of different
!    bands over which the radiation package calculates aerosol 
!    radiative properties.
!---------------------------------------------------------------------
integer, parameter ::    &
         N_AEROSOL_BANDS_FR = 8 ! number of non-continuum ir aerosol
                                ! emissivity bands 
integer, parameter ::     &
         N_AEROSOL_BANDS_CO = 1 ! number of continuum ir aerosol
                                ! emissivity bands  
integer, parameter ::    &
         N_AEROSOL_BANDS = N_AEROSOL_BANDS_FR + N_AEROSOL_BANDS_CO
                                ! total number of ir aerosol emissivity
                                ! bands 

!--------------------------------------------------------------------
!    num_wavenumbers is the number of wavenumber bands over which the
!    aerosol parameterization provides aerosol radiative property data.
!--------------------------------------------------------------------
integer     ::  num_wavenumbers = 0 ! number of wavenumber bands 
                                    ! present in the aerosol 
                                    ! parameterization

!----------------------------------------------------------------------
!    the following variable defines the number of aerosol property 
!    types that are active.
!----------------------------------------------------------------------
integer     :: naermodels = 0   ! number of aerosol optical properties
                                ! types that are active

!----------------------------------------------------------------------
!    Aerosol_props is an aerosol_properties_type variable which retains
!    the information needed throughout the run concerning aerosol
!    radiative properties and provides the means to transfer that
!    information to other modules.
!----------------------------------------------------------------------
type(aerosol_properties_type), save         :: Aerosol_props

!---------------------------------------------------------------------
!    the following arrays related to sw aerosol effects are allocated 
!    during initialization and retained throughout the integration.
!    here n refers to the bands of the solar parameterization, ni
!    to the bands of the aerosol parameterization, and na to the optical
!    properties type.
!
!      solivlaero(n,ni)  amount of toa incoming solar from solar
!                        spectral band n that is in aerosol parameter-
!                        ization band ni
!      nivl1aero(n)      the aerosol band index corresponding to the 
!                        lowest wave number of spectral band n
!      nivl2aero(n)      the aerosol band index corresponding to the 
!                        highest wave number of spectral band n
!      endaerwvnsf(ni)   ending wave number of aerosol parameterization
!                        band ni
!      aeroextivl(ni,na) extinction coefficient for aerosol parameter-
!                        ization band ni for aerosol optical property 
!                        type na
!      aerossalbivl(ni,na) 
!                        single-scattering albedo for aerosol band 
!                        ni and aerosol optical property type na
!      aeroasymmivl(ni,na)
!                        asymmetry factor for aerosol band ni  and 
!                        aerosol optical property type na
!
!---------------------------------------------------------------------
real,    dimension(:,:), allocatable   :: solivlaero  
integer, dimension(:),   allocatable   :: nivl1aero, nivl2aero
integer, dimension(:),   allocatable   :: endaerwvnsf
real,    dimension(:,:), allocatable   :: aeroextivl, aerossalbivl, &
                                          aeroasymmivl

!---------------------------------------------------------------------
!    sfl following arrays related to lw aerosol effects are allocated 
!    during initialization and retained throughout the integration.
!
!    sflwwts(n,ni)     the fraction of the planck function in aerosol 
!                      emissivity band n that is in aerosol param-
!                      eterization band ni
!
!----------------------------------------------------------------------
real,    dimension(:,:), allocatable   :: sflwwts

!--------------------------------------------------------------------
!    logical flags 
!--------------------------------------------------------------------
logical :: swprops_completed          = .false. ! sw properties have
                                                ! been calculated ?
logical :: band_calculation_completed = .false. ! lw properties have
                                                ! been calculated ?
logical :: module_is_initialized      = .false. ! module has been
                                                ! initialized ?
logical :: doing_predicted_aerosols   = .false. ! predicted aerosol 
                                                ! scheme being used ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------
 


                         contains
 

 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="aerosolrad_package_init">
!  <OVERVIEW>
!     aerosolrad_package_init is the constructor for 
!     aerosolrad_package_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!     aerosolrad_package_init is the constructor for 
!     aerosolrad_package_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosolrad_package_init (aerosol_names)
!  </TEMPLATE>
!  <IN NAME="aerosol_names" TYPE="character">
!   names of the activated aerosol species
!  </IN>
! </SUBROUTINE>
!
subroutine aerosolrad_package_init (aerosol_names)

!---------------------------------------------------------------------
!     aerosolrad_package_init is the constructor for 
!     aerosolrad_package_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!character(len=64), dimension(:), intent(in)  :: aerosol_names
character(len=*), dimension(:), intent(in)  :: aerosol_names
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  intent(in) variables:
!
!      aerosol_names     the names assigned to each of the activated
!                        aerosol species
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer        :: unit, ierr, io

!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call mpp_io_init
      call fms_init
      call rad_utilities_init
      call esfsw_parameters_init
      call longwave_params_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosolrad_package_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'aerosolrad_package_nml')
        end do
10      call close_file (unit)
      endif
 
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
                          write (stdlog(), nml=aerosolrad_package_nml)

!---------------------------------------------------------------------
!   exit if aerosols are desired with the lacis-hansen parameterization.
!---------------------------------------------------------------------
     if (Sw_control%do_lhsw_iz) then
       if (Sw_control%do_lhsw .and. do_swaerosol) then
         call error_mesg ('aerosolrad_package_mod', &
         ' cannot activate sw aerosols with lhsw', FATAL)
        endif
     else
       call error_mesg ('aerosolrad_package_mod', &
        'Sw_control%do_lhsw not yet initialized', FATAL)
     endif

!----------------------------------------------------------------------
!    define control variables which indicate whether the impact of 
!    aerosols on radiation is to be included in the sw and lw rad-
!    iation calculations. define a control variable which will be true 
!    if aerosols are included in either the sw or the lw radiation 
!    (Rad_control%do_aerosol).
!----------------------------------------------------------------------
      Sw_control%do_swaerosol = do_swaerosol
      Lw_control%do_lwaerosol = do_lwaerosol
      if (do_lwaerosol .or. do_swaerosol ) then   
        Rad_control%do_aerosol = .true.
      else
        Rad_control%do_aerosol = .false.
      endif

!--------------------------------------------------------------------
!    mark the just defined logicals as initialized.
!--------------------------------------------------------------------
      Sw_control%do_swaerosol_iz = .true.        
      Lw_control%do_lwaerosol_iz = .true.        
      Rad_control%do_aerosol_iz  = .true.
     
!---------------------------------------------------------------------
!    exit if an aerosol_data_set is provided when do_aerosol is 
!    .false..
!---------------------------------------------------------------------
      if ( .not. Rad_control%do_aerosol .and.    &
           trim(aerosol_data_set) /= ' ') then
        call error_mesg ('aerosolrad_package_mod', &
           'if aerosol impacts are not desired, aerosol_data_set '//&
            'must be set to "   "', FATAL)
      endif

!---------------------------------------------------------------------
!    exit if no aerosol_data_set is provided when do_aerosol  
!    is .true..
!---------------------------------------------------------------------
      if ( Rad_control%do_aerosol .and.    &
           trim(aerosol_data_set) == ' ') then
        call error_mesg ('aerosolrad_package_mod', &
           'if aerosol impacts are desired, aerosol_data_set '//&
            'must be non-blank', FATAL)
      endif

!---------------------------------------------------------------------
!    exit if aerosol effects are desired but the aerosol input file
!    provided no aerosol fields.
!---------------------------------------------------------------------
      if (Rad_control%do_aerosol .and. size(aerosol_names) == 0) then
        call error_mesg ('aerosolrad_package_mod', &
          ' aerosols desired  for radiation but no aerosol '//&
            'data_names supplied', FATAL)
      endif

!---------------------------------------------------------------------
!    if predicted aerosols are desired, exit with an error message.
!    code to handle prognostic aerosol has not yet been produced.
!---------------------------------------------------------------------
      if (Rad_control%do_aerosol .and.    &
          trim(aerosol_data_set) == 'predicted') then
        doing_predicted_aerosols = .true.
        call error_mesg ('aerosolrad_package_mod',  &
               'predicted aerosols not currently implemented', FATAL)
      endif

!----------------------------------------------------------------------
!    the only aerosol data set currently implemented is the
!    shettle-fenn data set.
!    references:                                            
!    shettle, e.p. and r.w. fenn, models for the aerosols of the lower 
!    atmosphere and the effects of humidity variations on their       
!    optical properties,afgl-tr-79-0214,1979,94pp.                    
!----------------------------------------------------------------------
      if (Rad_control%do_aerosol ) then         
        if (trim(aerosol_data_set) == 'shettle_fenn')  then

!-------------------------------------------------------------------
!  error condition -- prescribed type is not acceptable
!-------------------------------------------------------------------
        else
          call error_mesg ('aerosolrad_package_mod',  &
                       ' aerosol data set is not recognized.', FATAL)
        endif
      endif

!----------------------------------------------------------------------
!    if aerosol radiative effects are to be included, call 
!    assign_aerosol_opt_props to assign the proper aerosol 
!    properties type to each aerosol type. then call 
!    read_optical_input_file to read the optical input file contain-
!    ing the aerosol parameterization information and data.
!----------------------------------------------------------------------
      if (Rad_control%do_aerosol) then
        call assign_aerosol_opt_props (aerosol_names)
        call read_optical_input_file
      endif
 
!---------------------------------------------------------------------
!    if aerosol effects are to be included in the sw calculation,
!    call sw_aerosol_interaction to define the weights needed to prop-
!    erly map the input data from the aerosol parameterization bands to 
!    the solar parameterization bands that the model is using.
!--------------------------------------------------------------------
!     if (Rad_control%do_aerosol) then
!       call read_optical_input_file
!     endif
      if (do_swaerosol) then
        call sw_aerosol_interaction                 
      endif

!---------------------------------------------------------------------
!    if aerosol effects are to be included in the lw calculation,
!    call lw_aerosol_interaction to define the weights needed to prop-
!    erly map the input data from the aerosol parameterization bands to 
!    the solar parameterization bands that the model is using. if
!    they are not, indicate that this part of the code has been 
!    executed.
!---------------------------------------------------------------------
      if (do_lwaerosol) then
        call lw_aerosol_interaction
      else
        Lw_parameters%n_lwaerosol_bands_iz = .true.
      endif

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!----------------------------------------------------------------------


end subroutine aerosolrad_package_init




!####################################################################
! <SUBROUTINE NAME="aerosol_radiative_properties">
!  <OVERVIEW>
!    aerosol_radiative_properties defines and returns the radiative
!    properties for each aerosol properties type and for each solar
!    parameterization band in the shortwave and for each aerosol 
!    emissivity band in the longwave.
!  </OVERVIEW>
!  <DESCRIPTION>
!    aerosol_radiative_properties defines and returns the radiative
!    properties for each aerosol properties type and for each solar
!    parameterization band in the shortwave and for each aerosol 
!    emissivity band in the longwave.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_radiative_properties (is, ie, js, je, &
!                                         Aerosol, Aerosol_props_out)
!  </TEMPLATE>
!  <IN NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatology input
!  </IN>
!  <INOUT NAME="Aerosol_props_out" TYPE="aerosol_properties_type">
!   Aerosol radiative properties in radiation package
!  </INOUT>
!  <IN NAME="is, ie" TYPE="integer">
!   The longitude index of model physics window domain
!  </IN>
!  <IN NAME="js, je" TYPE="integer">
!   The latitude index of model physics window domain
!  </IN>
! </SUBROUTINE>
!
subroutine aerosol_radiative_properties (is, ie, js, je, &
                                         Aerosol, Aerosol_props_out)

!---------------------------------------------------------------------
!    aerosol_radiative_properties defines and returns the radiative
!    properties for each aerosol properties type and for each solar
!    parameterization band in the shortwave and for each aerosol 
!    emissivity band in the longwave.
!---------------------------------------------------------------------

integer,                       intent(in)    :: is, ie, js, je
type(aerosol_type),            intent(in)    :: Aerosol
type(aerosol_properties_type), intent(inout) :: Aerosol_props_out
 
!----------------------------------------------------------------------
! local variables:                                                     

      integer  :: na, nw, ni, nmodel       ! do-loop indices
     
!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosolrad_package_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    the following code applies to the non-predicted-aerosol case, 
!    where these calculations need be done only once and the 
!    results saved for access on succeeeding timesteps. if predicted
!    aerosols were being used, there may be a need to do new 
!    calculations on each timestep.
!---------------------------------------------------------------------
     if (.not. doing_predicted_aerosols) then

!--------------------------------------------------------------------
!    if sw aerosol properties are desired and have not yet been calc-
!    ulated, use the thick-averaging technique to define the single-
!    scattering properties for each solar parameterization band n 
!    from the specified properties on the aerosol parameterization 
!    bands ni for each aerosol properties type nmodel. 
! references:                                                          
!    edwards,j.m. and a. slingo, studies with a flexible new radiation  
!    code I: choosing a configuration for a large-scale model.,     
!    q.j.r. meteorological society, 122, 689-719, 1996.              
!                                                                      
! note: a thin-averaging technique (subroutine thinavg in 
!    rad_utilities_mod) is also available.   
!--------------------------------------------------------------------
        if (Sw_control%do_swaerosol .and. .not. swprops_completed) then 
          do nmodel=1,naermodels
            call thickavg (nivl1aero, nivl2aero, num_wavenumbers,   &
                           Solar_spect%nbands, aeroextivl(:,nmodel), &
                           aerossalbivl(:,nmodel),    &
                           aeroasymmivl(:,nmodel), solivlaero,   &
                           Solar_spect%solflxband,       & 
                           Aerosol_props%aerextband(:,nmodel),    &
                           Aerosol_props%aerssalbband(:,nmodel),   &
                           Aerosol_props%aerasymmband(:,nmodel))
          end do

!--------------------------------------------------------------------
!    mark the sw aerosol properties as having been calculated.
!--------------------------------------------------------------------
          swprops_completed = .true.
        endif

!---------------------------------------------------------------------
!    if longwave aerosol effects are desired, and the following cal-
!    culation has not already been done, calculate the aerosol 
!    properties for each aerosol properties type nw over each aerosol 
!    emissivity band na using the weighted contributions from each
!    aerosol parameterization band ni. mark the calculation as com-
!    pleted.
!---------------------------------------------------------------------
        if (Lw_control%do_lwaerosol .and.  &
            .not. band_calculation_completed) then
          do nw=1,naermodels    
            do na=1,N_AEROSOL_BANDS  
              do ni=1,num_wavenumbers 
                Aerosol_props%aerextbandlw(na,nw) =    &
                                  Aerosol_props%aerextbandlw(na,nw) + &
                                  aeroextivl(ni,nw)*sflwwts(na,ni)
                Aerosol_props%aerssalbbandlw(na,nw) =     &
                                  Aerosol_props%aerssalbbandlw(na,nw) +&
                                  aerossalbivl(ni,nw)*sflwwts(na,ni)
              end do
            end do
          end do
          band_calculation_completed = .true.
        endif ! (do_lwaerosol)

!---------------------------------------------------------------------
!    code for predicted aerosols will be placed here, when available.
!---------------------------------------------------------------------
      else  ! ( .not. doing_predicted_aerosols)
        call error_mesg ('aerosolrad_package_mod', &
                    'predicted aerosols not yet implenmented', FATAL)
      endif ! ( ,not. doing_predicted_aerosols)

!---------------------------------------------------------------------
!    return the aerosol_properties_type variable to the calling 
!    routine. this variable contains the aerosol radiative properties 
!    for each aerosol properties type over each solar and aerosol
!    emissivity band.
!---------------------------------------------------------------------
      Aerosol_props_out = Aerosol_props 

!--------------------------------------------------------------------


end subroutine aerosol_radiative_properties




!#####################################################################
! <SUBROUTINE NAME="aerosolrad_package_end">
!  <OVERVIEW>
!    aerosolrad_package_end is the destructor for 
!    aerosolrad_package_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    aerosolrad_package_end is the destructor for 
!    aerosolrad_package_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosolrad_package_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine aerosolrad_package_end

!--------------------------------------------------------------------
!    aerosolrad_package_end is the destructor for 
!    aerosolrad_package_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosolrad_package_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    deallocate module variables.
!---------------------------------------------------------------------
      if (do_swaerosol) then
        deallocate (solivlaero, nivl1aero, nivl2aero, endaerwvnsf, &
                    aeroextivl, aerossalbivl, aeroasymmivl)
      endif
      if (do_lwaerosol) then
        deallocate ( sflwwts)
      endif
      
!---------------------------------------------------------------------
!    deallocate elements of the aerosol_properties_type array.
!---------------------------------------------------------------------
      if (do_swaerosol) then
        deallocate ( Aerosol_props%aerextband,      &
                     Aerosol_props%aerssalbband,    &
                     Aerosol_props%aerasymmband     )
      endif
      if (do_lwaerosol) then
        deallocate ( Aerosol_props%aerextbandlw,    &
                     Aerosol_props%aerssalbbandlw   )
      endif
      if (Rad_control%do_aerosol) then
        deallocate ( Aerosol_props%sulfate_index,   &
                     Aerosol_props%optical_index    )
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.




end subroutine aerosolrad_package_end


!####################################################################
! <SUBROUTINE NAME="get_aerosol_optical_info">
!  <OVERVIEW>
!    get_aerosol_optical_info accesses data stored by this module.
!  </OVERVIEW>
!  <DESCRIPTION>
!    get_aerosol_optical_info accesses data stored by this module.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call get_aerosol_optical_info( num_categories, nwavenumbers, &
!                                     names, wavenumbers, &
!                                     aer_ext, aer_ss_alb, aer_asymm)
!  </TEMPLATE>
!  <OUT NAME="num_categories" TYPE="integer">
!   number of aerosol properties types
!  </OUT>
!  <OUT NAME="nwavenumbers" TYPE="integer">
!   number of wavenumber bands over which
!                           aerosol properties are defined
!  </OUT>
!  <OUT NAME="names" TYPE="character">
!   names assigned to the optical properties types
!  </OUT>
!  <OUT NAME="wavenumbers" TYPE="real">
!   wavenumber limits for each of the bands for
!                           which aerosol properties are defined
!  </OUT>
!  <OUT NAME="aer_ext, aer_ss_alb, aer_asymm" TYPE="real">
!   Aerosol extinction coefficient, single scattering albedo, and
!   asymmetry parameter
!  </OUT>
! </SUBROUTINE>
!
subroutine get_aerosol_optical_info( num_categories, nwavenumbers, &
                                     names, wavenumbers, &
                                     aer_ext, aer_ss_alb, aer_asymm)

!-----------------------------------------------------------------------
!    get_aerosol_optical_info accesses data stored by this module.
!-----------------------------------------------------------------------

integer,                        intent(out), optional ::       &
                                            num_categories, nwavenumbers
character(len=*), dimension(:), intent(out), optional :: names
integer, dimension(:),          intent(out), optional :: wavenumbers
real, dimension(:,:),           intent(out), optional :: aer_ext, &
                                                         aer_ss_alb, &
                                                         aer_asymm

!----------------------------------------------------------------------
!   intent(out), optional variables:
!
!      num_categories       number of aerosol properties types
!      nwavenumbers         number of wavenumber bands over which
!                           aerosol properties are defined
!      names                names assigned to the optical properties 
!                           types
!      wavenumbers          wavenumber limits for each of the bands for
!                           which aerosol properties are defined
!      aer_ext              extinction coefficient for each aerosol
!                           spectral band and each aerosol optical 
!                           property type
!      aer_ss_ab            single-scattering albedo for each aerosol 
!                           band and each aerosol optical property type 
!      aer_asymm            asymmetry factor for each aerosol band and 
!                           each aerosol optical property type
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosolrad_package_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    define the desired output variables.
!---------------------------------------------------------------------
      if( present(num_categories) ) num_categories = naermodels
      if( present(nwavenumbers))    nwavenumbers   = num_wavenumbers
      if( present(names) )          names(:naermodels) =    &
                                     aerosol_optical_names(:naermodels)
      if( present(wavenumbers) )    wavenumbers(:num_wavenumbers) =  &
                                           endaerwvnsf(:num_wavenumbers)
      if( present(aer_ext) )        aer_ext(:,:)    = aeroextivl(:,:)
      if( present(aer_ss_alb) )     aer_ss_alb(:,:) = aerossalbivl(:,:)
      if( present(aer_asymm) )      aer_asymm(:,:)  = aeroasymmivl(:,:)

!---------------------------------------------------------------------

end subroutine get_aerosol_optical_info


!######################################################################
! <SUBROUTINE NAME="get_aerosol_optical_index">
!  <OVERVIEW>
!    get_aerosol_optical_index returns the aerosol optical property
!    index for given aerosol number and relative humidity.
!  </OVERVIEW>
!  <DESCRIPTION>
!    get_aerosol_optical_index returns the aerosol optical property
!    index for given aerosol number and relative humidity.
!  </DESCRIPTION>
!  <TEMPLATE>
!   index = get_aerosol_optical_index( name, naerosol, rh )
!  </TEMPLATE>
!  <IN NAME="name" TYPE="real">
!   aerosol species name for which the optical 
!                      properties index is desired
!  </IN>
!  <IN NAME="naerosol" TYPE="integer">
!   aerosol index of the aerosol for whoch the 
!                      optical properties index is desired
!  </IN>
!  <IN NAME="rh" TYPE="real">
!    relative humidity 
!  </IN>
! </SUBROUTINE>
!
function get_aerosol_optical_index( name, naerosol, rh ) result(index)

!-----------------------------------------------------------------------
!    get_aerosol_optical_index returns the aerosol optical property
!    index for given aerosol number and relative humidity.
!-----------------------------------------------------------------------

!character(len=64),         intent(in) :: name
character(len=*),         intent(in) :: name
integer,                   intent(in) :: naerosol
real,                      intent(in) :: rh
integer                               :: index  ! function value

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      name            aerosol species name for which the optical 
!                      properties index is desired
!      naerosol        aerosol index of the aerosol for whoch the 
!                      optical properties index is desired
!      rh              relative humidity 
!
!  function value
!
!      index           returned optical properties index for aerosol
!                      name (aerosol index = naerosol) when the 
!                      relative humidity is rh
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer   :: irh     ! integer value for relative humidity, 
                           ! used as an index
      integer   :: nfields ! total number of active aerosols

!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module is initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg( 'aerosolrad_package_mod',  &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    be sure the desired aerosol index is valid.
!---------------------------------------------------------------------
      nfields = size(Aerosol_props%optical_index)
      if (naerosol > nfields) then
        call error_mesg( 'aerosolrad_package_mod', &
           'aerosol index exceeds number of aerosol fields', FATAL )
      end if

!---------------------------------------------------------------------
!    determine if the desired aerosol is a sulfate or non-sulfate. if
!    sulfate, then the optical propeties index will depend on the
!    relative humidity.
!---------------------------------------------------------------------
      if (Aerosol_props%optical_index(naerosol) == 0 ) then
        irh = MIN( 100, MAX( 0, NINT(100.*rh) ) )
        index = Aerosol_props%sulfate_index( irh )
      else
        index = Aerosol_props%optical_index(naerosol)
      endif

!---------------------------------------------------------------------
!    if no value was obtained for the optical index, stop execution.
!---------------------------------------------------------------------
      if (index == 0 ) then
        call error_mesg ('aerosolrad_package_mod', &
           'Cannot find aerosol optical properties for species = ' // &
                  trim (name), FATAL )
      endif

!----------------------------------------------------------------------


end function get_aerosol_optical_index




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                
!                    PRIVATE SUBROUTINES
!                                
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  
                                  

!#####################################################################
! <SUBROUTINE NAME="assign_aerosol_opt_props">
!  <OVERVIEW>
!    assign_aerosol_opt_props assigns an index for an available optical
!    properties type to each activated aerosol type. for sulfates, a 
!    flag is set, since the aerosol properties type is a function 
!    of model relative humidity, and will vary with time.
!  </OVERVIEW>
!  <DESCRIPTION>
!    assign_aerosol_opt_props assigns an index for an available optical
!    properties type to each activated aerosol type. for sulfates, a 
!    flag is set, since the aerosol properties type is a function 
!    of model relative humidity, and will vary with time.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call assign_aerosol_opt_props (aerosol_names)
!  </TEMPLATE>
!  <IN NAME="aerosol_names" TYPE="character">
!   names associated with each aerosol species
!  </IN>
! </SUBROUTINE>
!
subroutine assign_aerosol_opt_props (aerosol_names)

!----------------------------------------------------------------------
!    assign_aerosol_opt_props assigns an index for an available optical
!    properties type to each activated aerosol type. for sulfates, a 
!    flag is set, since the aerosol properties type is a function 
!    of model relative humidity, and will vary with time.
!---------------------------------------------------------------------

!character(len=64), dimension(:), intent(in) :: aerosol_names
character(len=*), dimension(:), intent(in) :: aerosol_names

!----------------------------------------------------------------------
!  intent(in) variables:
!
!     aerosol_names     names associated with each aerosol species
!
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    local variables:

      character(len=64) :: name_in, target_name
      integer           :: nfields
      integer           :: n, noptical

!---------------------------------------------------------------------
!   local variables:
!
!       name_in          variable to hold current aerosol name 
!                        being processed
!       target_name      aerosol_optical_name associated with a given
!                        aerosol species     
!       nfields          number of activated aerosol species
!       n, noptical      do-loop indices
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    count the number of aerosol optical property categories requested
!    via the namelist input.
!----------------------------------------------------------------------
      do n=1,MAX_OPTICAL_FIELDS
        if (aerosol_optical_names(n) /= ' '  ) then
          naermodels = n
        else
          exit
        endif
      end do

!---------------------------------------------------------------------
!    define the number of activated aerosol species.
!---------------------------------------------------------------------
      nfields = size (aerosol_names)

!---------------------------------------------------------------------
!    allocate components of the aerosol_properties_type module variable
!    which will contain the indices for the different aerosols.
!---------------------------------------------------------------------
       allocate (Aerosol_props%sulfate_index (0:100 ) )
       allocate (Aerosol_props%optical_index(nfields) )

!----------------------------------------------------------------------
!    match aerosol optical property indices with aerosol indices.
!    sulfate aerosols are handled separately (below) with RH dependence.
!----------------------------------------------------------------------
      do n=1,nfields
        name_in = aerosol_names(n)
        Aerosol_props%optical_index(n) = 0
        if (name_in /= "so4_anthro" .and.    &
            name_in /= "so4_natural" ) then
          select case( name_in )
            case( "anthro_dust_0.1" )
              target_name = "dust_0.1"
            case( "anthro_dust_0.2" )
              target_name = "dust_0.2"
            case( "anthro_dust_0.4" )
              target_name = "dust_0.4"
            case( "anthro_dust_0.8" )
              target_name = "dust_0.8"
            case( "anthro_dust_1.0" )
              target_name = "dust_1.0"
            case( "anthro_dust_2.0" )
              target_name = "dust_2.0"
            case( "anthro_dust_4.0" )
              target_name = "dust_4.0"
            case( "anthro_dust_8.0" )
              target_name = "dust_8.0"
            case( "natural_dust_0.1" )
              target_name = "dust_0.1"
            case( "natural_dust_0.2" )
              target_name = "dust_0.2"
            case( "natural_dust_0.4" )
              target_name = "dust_0.4"
            case( "natural_dust_0.8" )
              target_name = "dust_0.8"
            case( "natural_dust_1.0" )
              target_name = "dust_1.0"
            case( "natural_dust_2.0" )
              target_name = "dust_2.0"
            case( "natural_dust_4.0" )
              target_name = "dust_4.0"
            case( "natural_dust_8.0" )
              target_name = "dust_8.0"
            case( "black_carbon" )
              target_name = "soot"
            case DEFAULT
              target_name = name_in
          end select  

!--------------------------------------------------------------------
!    go through the set of aerosol properties types looking for 
!    the target_name defined above. when found, associate the
!    optical properties type index with the current aerosol species.
!--------------------------------------------------------------------
          do noptical=1,naermodels
            if (aerosol_optical_names(noptical) == target_name) then
              Aerosol_props%optical_index(n) = noptical
              exit
            end if
          end do

!--------------------------------------------------------------------
!    if the target_name is not found, exit with an error message.
!----------------------------------------------------------------------
          if (Aerosol_props%optical_index(n) == 0 ) then
            call error_mesg( 'aerosolrad_package_mod', &
                'Cannot find aerosol optical model = ' //    &
                                           TRIM( target_name ), FATAL )
          endif
        endif
      end do

!----------------------------------------------------------------------
!    set up RH-dependent sulfate aerosol optical property indices.
!-------------------------------------------------------------------
      do n=0,100

!--------------------------------------------------------------------
!    define the optical properties type for all possible values of 
!    relative humidity.
!--------------------------------------------------------------------
        select case( n )
          case ( 0:32 )
            target_name = "sulfate_30%"
          case ( 33:37 )
            target_name = "sulfate_35%"
          case ( 38:42 )
            target_name = "sulfate_40%"
          case ( 43:47 )
            target_name = "sulfate_45%"
          case ( 48:52 )
            target_name = "sulfate_50%"
          case ( 53:57 )
            target_name = "sulfate_55%"
          case ( 58:62 )
            target_name = "sulfate_60%"
          case ( 63:67 )
            target_name = "sulfate_65%"
          case ( 68:72 )
            target_name = "sulfate_70%"
          case ( 73:77 )
            target_name = "sulfate_75%"
          case ( 78:81 )
            target_name = "sulfate_80%"
          case ( 82:83 )
            target_name = "sulfate_82%"
          case ( 84:85 )
            target_name = "sulfate_84%"
          case ( 86:87 )
            target_name = "sulfate_86%"
          case ( 88:89 )
            target_name = "sulfate_88%"
          case ( 90 )
            target_name = "sulfate_90%"
          case ( 91 )
            target_name = "sulfate_91%"
          case ( 92 )
            target_name = "sulfate_92%"
          case ( 93 )
            target_name = "sulfate_93%"
          case ( 94 )
            target_name = "sulfate_94%"
          case ( 95 )
            target_name = "sulfate_95%"
          case ( 96 )
            target_name = "sulfate_96%"
          case ( 97 )
            target_name = "sulfate_97%"
          case ( 98 )
            target_name = "sulfate_98%"
          case ( 99 )
            target_name = "sulfate_99%"
          case ( 100 )
            target_name = "sulfate_100%"
        end select

!---------------------------------------------------------------------
!    associate an index value with each possible relative humidity.
!---------------------------------------------------------------------
        Aerosol_props%sulfate_index(n) = 0
        do noptical=1,naermodels
          if (aerosol_optical_names(noptical) == target_name ) then
            Aerosol_props%sulfate_index(n) = noptical
            exit
          end if
        end do

!---------------------------------------------------------------------
!    if the  aerosol_optical name_is not included in the potential
!    set listed above, exit with an error message.
!---------------------------------------------------------------------
        if (Aerosol_props%sulfate_index(n) == 0 ) then
          call error_mesg( 'aerosolrad_package_mod', &
                 'Cannot find aerosol optical model = ' // &
                                          TRIM( target_name), FATAL )
        endif
      end do

!---------------------------------------------------------------------



end subroutine assign_aerosol_opt_props



!######################################################################
! <SUBROUTINE NAME="read_optical_input_file">
!  <OVERVIEW>
!    read_optical_input_file reads the optical properties input file
!    to obtain the specified aerosol radiative properties for each 
!    aerosol in each of the aerosol parameterization spectral bands.
!  </OVERVIEW>
!  <DESCRIPTION>
!    read_optical_input_file reads the optical properties input file
!    to obtain the specified aerosol radiative properties for each 
!    aerosol in each of the aerosol parameterization spectral bands.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call read_optical_input_file
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine read_optical_input_file

!-----------------------------------------------------------------------
!    read_optical_input_file reads the optical properties input file
!    to obtain the specified aerosol radiative properties for each 
!    aerosol in each of the aerosol parameterization spectral bands.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    local variables:

      real,    dimension(:), allocatable    :: aeroext_in,   &
                                               aerossalb_in,   &
                                               aeroasymm_in
      logical, dimension(:), allocatable    :: found

      integer           :: unit, num_input_categories
      character(len=64) :: name_in
      integer           :: n, noptical

!---------------------------------------------------------------------
!   local variables:
!
!       aeroext_in       aerosol extinction coefficient read from 
!                        input file
!       aerossalb_in     aerosol single scattering albedo read from 
!                        input file
!       aeroasymm_in     aerosol asymmetry factor read from 
!                        input file
!       found            aerosol radiative property data has been
!                        obtained from input file for the given
!                        optical properties type ?
!       unit             io unit number used for optical properties file
!       num_input_categories
!                        number of optical properties types contained
!                        in optical data input file
!       name_in          name of optical properties type being processed
!       n, noptical      do-loop indices
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!    open the ASCII input file containing aerosol optical property
!    information.
!----------------------------------------------------------------------
      call mpp_open (unit, 'INPUT/'//optical_filename, MPP_RDONLY,  &
                     MPP_ASCII, MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE)

!----------------------------------------------------------------------
!    read the dimension information contained in the input file.
!----------------------------------------------------------------------
      read ( unit,* ) num_wavenumbers
      read ( unit,* ) num_input_categories

!----------------------------------------------------------------------
!    read wavenumber limits for aerosol parameterization bands from 
!    the input file.
!----------------------------------------------------------------------
       allocate (endaerwvnsf(num_wavenumbers) )
       read (unit,* )
       read (unit,* ) endaerwvnsf
 
!----------------------------------------------------------------------
!    allocate module arrays to hold the specified sw properties for 
!    each parameterization bnad and each aerosol properties type.
!----------------------------------------------------------------------
      allocate (       &
            aeroextivl   (num_wavenumbers, naermodels),&
            aerossalbivl (num_wavenumbers, naermodels), &
            aeroasymmivl (num_wavenumbers, naermodels) )

!----------------------------------------------------------------------
!    allocate local working arrays.
!----------------------------------------------------------------------
      allocate (aeroext_in   (num_wavenumbers ),             &
                aerossalb_in (num_wavenumbers ),           &
                aeroasymm_in (num_wavenumbers ),           &
                found        (naermodels ) )

!----------------------------------------------------------------------
!    match the names of optical property categories from input file with
!    those specified in the namelist, and store the following data
!    appropriately. indicate that the data has been found.
!----------------------------------------------------------------------
      found(:) = .false.
      do n=1,num_input_categories
        read( unit,* ) name_in
        read( unit,* )
        read( unit,* ) aeroext_in
        read( unit,* )
        read( unit,* ) aerossalb_in
        read( unit,* )
        read( unit,* ) aeroasymm_in
        do noptical=1,naermodels
          if (aerosol_optical_names(noptical) == name_in) then
            aeroextivl(:,noptical)   = aeroext_in
            aerossalbivl(:,noptical) = aerossalb_in
            aeroasymmivl(:,noptical) = aeroasymm_in
            found( noptical ) = .true.
            exit
          endif
        end do
      end do

!----------------------------------------------------------------------
!    close the ASCII input file.
!----------------------------------------------------------------------
      call mpp_close( unit )

!----------------------------------------------------------------------
!    check to make sure data for all aerosol optical property
!    categories specified in namelist were contained in ASCII
!    input file. if not, exit with a message.
!----------------------------------------------------------------------
      do noptical = 1,naermodels
        if (.not. found( noptical ) ) then
              call error_mesg( 'aerosolrad_package_mod', &
              'Cannot find aerosol optical properties for ' // &
                TRIM(aerosol_optical_names(noptical)),  FATAL )
        endif
      end do

!----------------------------------------------------------------------
!    deallocate local working arrays.
!----------------------------------------------------------------------
      deallocate (aeroext_in, aerossalb_in, aeroasymm_in, found)



end subroutine read_optical_input_file



!#####################################################################
! <SUBROUTINE NAME="sw_aerosol_interaction">
!  <OVERVIEW>
!    sw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the solar parameterization
!    bands being used by the model.
!  </OVERVIEW>
!  <DESCRIPTION>
!    sw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the solar parameterization
!    bands being used by the model.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call sw_aerosol_interaction
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine sw_aerosol_interaction

!-----------------------------------------------------------------------
!    sw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the solar parameterization
!    bands being used by the model.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    local variables:

      real,    dimension(:), allocatable    :: aeroext_in,   &
                                               aerossalb_in,   &
                                               aeroasymm_in
      logical, dimension(:), allocatable    :: found

      integer           :: unit, num_input_categories
      character(len=64) :: name_in
      integer           :: nbands, nband, nivl3
      real              :: sumsol3
      integer           :: n, nw, noptical

!---------------------------------------------------------------------
!   local variables:
!
!       aeroext_in       aerosol extinction coefficient read from 
!                        input file
!       aerossalb_in     aerosol single scattering albedo read from 
!                        input file
!       aeroasymm_in     aerosol asymmetry factor read from 
!                        input file
!       found            aerosol radiative property data has been
!                        obtained from input file for the given
!                        optical properties type ?
!       unit             io unit number used for optical properties file
!       num_input_categories
!                        number of optical properties types contained
!                        in optical data input file
!       name_in          name of optical properties type being processed
!       nbands           number of bands in solar spectral param-
!                        eterization
!       nband            currently active solar spectrum band 
!       nivl3            currently active aerosol parameterization band
!       sumsol3          sum of solar input in current aerosol param-
!                        eterization band
!       n, nw, noptical  do-loop indices
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    allocate and initialize the array components of an 
!    aerosol_properties_type variable which will hold the aerosol 
!    radiative properties for each solar spectral parameterization band.
!    aerextband     the solar band values of the extinction 
!                   coefficient for aerosols                           
!    aerssalbband   the solar band values of the single-     
!                   scattering albedo for aerosols                      
!    aerasymmband   the solar band values of the asymmetry   
!                   factor for aerosols                                 
!---------------------------------------------------------------------
      allocate     &
        (Aerosol_props%aerextband   (Solar_spect%nbands, naermodels), &
         Aerosol_props%aerssalbband (Solar_spect%nbands, naermodels), &
         Aerosol_props%aerasymmband (Solar_spect%nbands, naermodels) )
      Aerosol_props%aerextband   = 0.
      Aerosol_props%aerssalbband = 0.
      Aerosol_props%aerasymmband = 0.

!---------------------------------------------------------------------
!    define the number of bands in the solar spectrum parameterization.
!    allocate space for variables defining the highest and lowest 
!    aerosol parameterization wavenumber in each solar spectral band
!    and the solar flux common to solar spectral band n and aerosol
!    parameterization band ni.
!---------------------------------------------------------------------
      nbands = Solar_spect%nbands 
      allocate ( nivl1aero  (nbands) )
      allocate ( nivl2aero  (nbands) )
      allocate ( solivlaero (nbands, num_wavenumbers))

!---------------------------------------------------------------------
!    define the solar weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the solar
!    spectral intervals and so determine the single-scattering proper-
!    ties on the solar spectral intervals.
!--------------------------------------------------------------------
      nivl3 = 1
      sumsol3 = 0.0
      nband = 1
      solivlaero(:,:) = 0.0
      nivl1aero(1) = 1
      do nw = 1,Solar_spect%endwvnbands(nbands)
        sumsol3 = sumsol3 + Solar_spect%solarfluxtoa(nw)
        if (nw == endaerwvnsf(nivl3) ) then
          solivlaero(nband,nivl3) = sumsol3
          sumsol3 = 0.0
        end if
        if ( nw == Solar_spect%endwvnbands(nband) ) then
          if ( nw /= endaerwvnsf(nivl3) ) then
            solivlaero(nband,nivl3) = sumsol3 
            sumsol3 = 0.0
          end if
          nivl2aero(nband) = nivl3
          nband = nband + 1
          if ( nband <= nbands ) then
            if ( nw == endaerwvnsf(nivl3) ) then
              nivl1aero(nband) = nivl3 + 1
            else
              nivl1aero(nband) = nivl3
            end if
          end if
        end if
        if ( nw == endaerwvnsf(nivl3) ) nivl3 = nivl3 + 1
      end do

!---------------------------------------------------------------------



end subroutine sw_aerosol_interaction   



!#####################################################################
! <SUBROUTINE NAME="lw_aerosol_interaction">
!  <OVERVIEW>
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!  </OVERVIEW>
!  <DESCRIPTION>
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call lw_aerosol_interaction
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine lw_aerosol_interaction      

!----------------------------------------------------------------------
!    lw_aerosol_interaction defines the weights and interval infor-
!    mation needed to map the aerosol radiative properties from the
!    aerosol parameterization bands to the aerosol emissivity bands
!    being used by the model.
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  local variables:

!---------------------------------------------------------------------
!    the following arrays define the wavenumber ranges for the separate
!    aerosol emissivity bands in the model infrared parameterization. 
!    these may be changed only by the keeper of the radiation code.
!    the order of the frequency bands corresponds to the order used
!    in the lw radiation code.
!
!      aerbandlo_fr      low wavenumber limit for the non-continuum 
!                        aerosol emissivity bands
!      aerbandhi_fr      high wavenumber limit for the non-continuum
!                        aerosol emissivity bands
!      istartaerband_fr  starting wavenumber index for the non-continuum
!                        aerosol emissivity bands
!      iendaerband_fr    ending wavenumber index for the non-continuum
!                        aerosol emissivity bands
!      aerbandlo_co      low wavenumber limit for the continuum 
!                        aerosol emissivity bands
!      aerbandhi_co      high wavenumber limit for the continuum
!                        aerosol emissivity bands
!      istartaerband_co  starting wavenumber index for the continuum
!                        aerosol emissivity bands
!      iendaerband_co    ending wavenumber index for the continuum
!                        aerosol emissivity bands
!      aerbandlo         low wavenumber limit for the entire set of
!                        aerosol emissivity bands
!      aerbandhi         high wavenumber limit for the entire set of
!                        aerosol emissivity bands
!      istartaerband     starting wavenumber index for the entire set of
!                        aerosol emissivity bands
!      iendaerband       ending wavenumber index for the entire set of
!                        aerosol emissivity bands
!
!----------------------------------------------------------------------
      real, dimension (N_AEROSOL_BANDS_FR)     :: aerbandlo_fr =  &
      (/ 560.0, 630.0, 700.0, 800.0, 900.0,  990.0, 1070.0, 1200.0 /)

      real, dimension (N_AEROSOL_BANDS_FR)     :: aerbandhi_fr =  &
      (/ 630.0, 700.0, 800.0, 900.0, 990.0, 1070.0, 1200.0, 1400.0 /)

      integer, dimension (N_AEROSOL_BANDS_FR)  :: istartaerband_fr =  &
      (/ 57,  64,  71,  81,  91, 100, 108, 121 /)

      integer, dimension (N_AEROSOL_BANDS_FR)  :: iendaerband_fr =  &
      (/ 63,  70,  80,  90,  99, 107, 120, 140 /)

      real, dimension (N_AEROSOL_BANDS_CO)     :: aerbandlo_co =  &
      (/ 560.0 /)

      real, dimension (N_AEROSOL_BANDS_CO)     :: aerbandhi_co =  &
      (/ 800.0 /)

      integer, dimension (N_AEROSOL_BANDS_CO)  :: istartaerband_co =  &
      (/ 57  /)

      integer, dimension (N_AEROSOL_BANDS_CO)  :: iendaerband_co =  &
      (/ 80  /)

      real,    dimension(N_AEROSOL_BANDS)      :: aerbandlo, aerbandhi
      integer, dimension(N_AEROSOL_BANDS)      :: istartaerband,    &
                                                  iendaerband

!---------------------------------------------------------------------
!    the following arrays define how the ir aerosol band structure 
!    relates to the aerosol parameterization bands.
!
!      nivl1aer_fr(n)    aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber of the 
!                        non-continuum ir aerosol emissivity band n
!      nivl2aer_fr(n)    aerosol parameterization band index corres-
!                        ponding to the highest wavenumber of the 
!                        non-continuum ir aerosol emissivity band n
!      nivl1aer_co(n)    aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber of the 
!                        continuum ir aerosol emissivity band n
!      nivl2aer_co(n)    aerosol parameterization band index corres-
!                        ponding to the highest wavenumber of the 
!                        continuum ir aerosol emissivity band n
!      nivl1aer(n)       aerosol parameterization band index corres-
!                        ponding to the lowest wavenumber for the 
!                        ir aerosol emissivity band n
!      nivl2aer(n)       aerosol parameterization band index corres-
!                        ponding to the highest wavenumber for the 
!                        ir aerosol emissivity band n
!      planckaerband(n)  planck function summed over each lw param-
!                        eterization band that is contained in the 
!                        ir aerosol emissivity band n
!
!---------------------------------------------------------------------
      integer, dimension (N_AEROSOL_BANDS_FR)  :: nivl1aer_fr,   &
                                                  nivl2aer_fr
      integer, dimension (N_AEROSOL_BANDS_CO)  :: nivl1aer_co,   &
                                                  nivl2aer_co
      integer, dimension (N_AEROSOL_BANDS)     :: nivl1aer, nivl2aer
      real,    dimension (N_AEROSOL_BANDS)     :: planckaerband

!----------------------------------------------------------------------
!    the following arrays relate the ir aerosol emissivity band n to
!    either the aerosol optical properties type na or to the aerosol 
!    parameterization band ni.
!        aerextbandlw_fr(n,na)  band averaged extinction coefficient
!                               for non-continuum aerosol emissivity 
!                               band n and aerosol properties type na
!        aerssalbbandlw_fr(n,na)
!                               band averaged single-scattering
!                               coefficient for non-continuum aerosol
!                               emissivity band n and aerosol properties
!                               type na
!        aerextbandlw_co(n,na)  band averaged extinction coefficient
!                               for the continuum aerosol emissivity
!                               band n and aerosol properties type na
!        aerssalbbandlw_co(n,na)
!                               band averaged single-scattering
!                               coefficient for continuum aerosol
!                               emissivity band n and aerosol properties
!                               type na
!        planckivlaer_fr(n,ni)  planck function over the spectral range
!                               common to aerosol emissivity non-
!                               continuum band n and aerosol parameter-
!                               ization band ni
!        planckivlaer_co(n,ni)  planck function over the spectral range
!                               common to aerosol emissivity continuum 
!                               band n and aerosol parameterization 
!                               band ni
!        sflwwts_fr(n,ni)       band weights for the aerosol emissivity
!                               non-continuum band n and the aerosol 
!                               parameterization band ni 
!        sflwwts_co(n,ni)       band weights for the aerosol emissivity
!                               continuum band n and the aerosol 
!                               parameterization band ni 
!        planckivlaer(n,ni)     planck function over the spectral range
!                               common to aerosol emissivity band n and
!                               aerosol parameterization band ni
!        iendsfbands(ni)        ending wavenumber index for aerosol 
!                               parameterization band ni
!
!----------------------------------------------------------------------
      real,    dimension (N_AEROSOL_BANDS_FR, naermodels) ::   &
                                                  aerextbandlw_fr, &
                                                  aerssalbbandlw_fr
      real,    dimension (N_AEROSOL_BANDS_CO, naermodels) ::   &
                                                  aerextbandlw_co, &
                                                  aerssalbbandlw_co
      real,    dimension (N_AEROSOL_BANDS_FR, num_wavenumbers) :: &
                                                  planckivlaer_fr, &
                                                  sflwwts_fr
      real,    dimension (N_AEROSOL_BANDS_CO, num_wavenumbers) :: &
                                                  planckivlaer_co, &
                                                  sflwwts_co
      real,    dimension (N_AEROSOL_BANDS, num_wavenumbers)  ::    &
                                                  planckivlaer
      integer, dimension (num_wavenumbers)    ::  iendsfbands

!---------------------------------------------------------------------
!    variables associated with the planck function calculation.
!    the planck function is defined for each of the NBLW longwave 
!    parameterization bands.
!---------------------------------------------------------------------
      real, dimension(NBLW)  :: c1, centnb, sc, src1nb, x, x1
      real                   :: del, xtemv, sumplanck

!---------------------------------------------------------------------
!    miscellaneous variables:

     logical         :: do_band1   !  should we do special calculation 
                                   !  for band 1 ?
     integer         :: ib, nw, nivl, nband, n, ni 
                                   !  do-loop indices and counters

!--------------------------------------------------------------------
!    define arrays containing the characteristics of all the ir aerosol
!    emissivity bands, both continuum and non-continuum.
!--------------------------------------------------------------------
      do n=1,N_AEROSOL_BANDS_FR
        aerbandlo(n)     = aerbandlo_fr(n)
        aerbandhi(n)     = aerbandhi_fr(n)
        istartaerband(n) = istartaerband_fr(n)
        iendaerband(n)   = iendaerband_fr(n)
      end do
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        aerbandlo(n)     = aerbandlo_co     (n - N_AEROSOL_BANDS_FR)
        aerbandhi(n)     = aerbandhi_co     (n - N_AEROSOL_BANDS_FR)
        istartaerband(n) = istartaerband_co (n - N_AEROSOL_BANDS_FR)
        iendaerband(n)   = iendaerband_co   (n - N_AEROSOL_BANDS_FR)
      end do

!---------------------------------------------------------------------
!    define the number of aerosol ir bands to be used in other modules.
!    set the initialization flag to .true.
!---------------------------------------------------------------------
      Lw_parameters%n_lwaerosol_bands = N_AEROSOL_BANDS
      Lw_parameters%n_lwaerosol_bands_iz = .true.

!--------------------------------------------------------------------
!    allocate a module variable which will store the weighting function
!    between the aerosol emissivity bands and the aerosol parameter-
!    ization bands.
!--------------------------------------------------------------------
      allocate (sflwwts (N_AEROSOL_BANDS, num_wavenumbers))

!--------------------------------------------------------------------
!    define the ending aerosol band index for each of the aerosol
!    parameterization bands.
!--------------------------------------------------------------------
      iendsfbands(:) = INT((endaerwvnsf(:) + 0.01)/10.0)

!--------------------------------------------------------------------
!    compute the planck function at 10C over each of the longwave
!    parameterization bands to be used as the weighting function. 
!--------------------------------------------------------------------
      do n=1,NBLW 
        del  = 10.0E+00
        xtemv = 283.15
        centnb(n) = 5.0 + (n - 1)*del
        c1(n)     = (3.7412E-05)*centnb(n)**3
        x(n)      = 1.4387E+00*centnb(n)/xtemv
        x1(n)     = EXP(x(n))
        sc(n)     = c1(n)/(x1(n) - 1.0E+00)
        src1nb(n) = del*sc(n)
      end do
 
!--------------------------------------------------------------------
!    sum the weighting function calculated over the longwave param-
!    eterization bands that are contained in each of the aerosol 
!    emissivity bands. 
!--------------------------------------------------------------------
      planckaerband(:) = 0.0E+00
      do n = 1,N_AEROSOL_BANDS
        do ib = istartaerband(n),iendaerband(n)
          planckaerband(n) = planckaerband(n) + src1nb(ib)
        end do
      end do
 
!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the non-
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_fr(:,:) = 0.0
      nivl1aer_fr(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_fr(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_fr(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_fr(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_fr(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_FR ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_fr(nband) = nivl + 1
            else
              nivl1aer_fr(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband .eq. 1 .and.   &
              iendsfbands(nivl-1) >= istartaerband_fr(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_fr(1)) then
            nivl1aer_fr(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if (nw >= iendaerband_fr(N_AEROSOL_BANDS_FR) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the weights and interval counters that are needed to  
!    map the aerosol parameterization spectral intervals onto the 
!    continuum ir aerosol emissivity bands and so determine the 
!    single-scattering properties on the ir aerosol emissivity bands.
!--------------------------------------------------------------------
      nivl = 1
      sumplanck = 0.0
      nband = 1
      planckivlaer_co(:,:) = 0.0
      nivl1aer_co(1) = 1
      do_band1 = .true.
 
      do nw = 1,NBLW
        sumplanck = sumplanck + src1nb(nw)
        if ( nw == iendsfbands(nivl) ) then
          planckivlaer_co(nband,nivl) = sumplanck
          sumplanck = 0.0
        end if
        if ( nw == iendaerband_co(nband) ) then
          if ( nw /= iendsfbands(nivl) ) then
            planckivlaer_co(nband,nivl) = sumplanck 
            sumplanck = 0.0
          end if
          nivl2aer_co(nband) = nivl
          nband = nband + 1
          if ( nband <= N_AEROSOL_BANDS_CO ) then
            if ( nw == iendsfbands(nivl) ) then
              nivl1aer_co(nband) = nivl + 1
            else
              nivl1aer_co(nband) = nivl
            end if
          end if
        end if
        if ( nw == iendsfbands(nivl) ) then
          nivl = nivl + 1
          if (do_band1 .and. nband == 1 .and.  &
              iendsfbands(nivl-1) >= istartaerband_co(1) .and.  &
              iendsfbands(nivl-1) < iendaerband_co(1)) then
            nivl1aer_co(nband) = nivl-1
            do_band1 = .false.
          endif
        endif
        if ( nw >= iendaerband_co(N_AEROSOL_BANDS_CO) ) then
          exit
        endif
      end do

!--------------------------------------------------------------------
!    define the planck-function-weighted band weights for the aerosol
!    parameterization bands onto the non-continuum and continuum ir 
!    aerosol emissivity bands.
!--------------------------------------------------------------------
      sflwwts_fr(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_FR
        do ni=nivl1aer_fr(n),nivl2aer_fr(n)
          sflwwts_fr(n,ni) = planckivlaer_fr(n,ni)/planckaerband(n)
        end do
      end do
      sflwwts_co(:,:) = 0.0E+00
      do n=1,N_AEROSOL_BANDS_CO
        do ni=nivl1aer_co(n),nivl2aer_co(n)
          sflwwts_co(n,ni) = planckivlaer_co(n,ni)/     &
                             planckaerband(N_AEROSOL_BANDS_FR+n)
        end do
      end do

!--------------------------------------------------------------------
!    consolidate the continuum and non-continuum weights into an
!    array covering all ir aerosol emissivity bands.
!--------------------------------------------------------------------
      do n=1,N_AEROSOL_BANDS_FR
        do ni = 1,num_wavenumbers
          sflwwts(n,ni) = sflwwts_fr(n,ni)
        end do
      end do
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        do ni = 1,num_wavenumbers
          sflwwts(n,ni) = sflwwts_co(n-N_AEROSOL_BANDS_FR,ni)
        end do
      end do

!-----------------------------------------------------------------
!    allocate and initialize the arrays in the aerosol_properties_type 
!    variable that will contain the ir aerosol properties for each
!    aerosol optical type over each ir aerosol emissivity band.
!----------------------------------------------------------------
      allocate     &
         (Aerosol_props%aerextbandlw   (N_AEROSOL_BANDS, naermodels), &
          Aerosol_props%aerssalbbandlw (N_AEROSOL_BANDS, naermodels) )
      Aerosol_props%aerextbandlw   = 0.0E+00
      Aerosol_props%aerssalbbandlw = 0.0E+00

!----------------------------------------------------------------------



end subroutine lw_aerosol_interaction


!######################################################################



               end module aerosolrad_package_mod





