                       module rad_output_file_mod


use  utilities_mod,     only:  open_file, file_exist,    &
                               check_nml_error, error_mesg, &
                               print_version_number, FATAL, NOTE, &
                               WARNING, get_my_pe, close_file, &
                               utilities_init, &
                               get_num_pes
use time_manager_mod,   only:     &
!	                       time_manager_init, &
                               time_type
use diag_manager_mod,   only:  register_diag_field,  &
                               diag_manager_init, &
                               send_data
use rad_utilities_mod,  only:  Environment, environment_type, &
                               rad_utilities_init, &
                               radiative_gases_type, &
                               rad_output_type, &
                               cldrad_properties_type, &
                               cld_diagnostics_type, &
                               atmos_input_type, &
                               sw_output_type,  &
                               lw_output_type,  &
                               Cldrad_control,  &
                               cloudrad_control_type, &
                               Rad_control, radiation_control_type


implicit none
private

!------------------------------------------------------------------
!    rad_output_file_mod writes an output file containing an assort-
!    ment of variables related to the sea_esf_rad radiation package.
!    this is an optionally-generated file, which may be used to sup-
!    plement the standard diagnostic model output files. NOTE THAT
!    THIS FILE IS GENERATED ONLY ON RADIATION TIMESTEPS, SO THAT WHEN 
!    SW FLUXES ARE BEING RENORMALIZED ON EACH PHYSICS STEP, VARIABLES 
!    IN THIS FILE RELATED TO SW RADIATION WILL NOT REFLECT THE EFFECTS 
!    OF THE RENORMALIZATION.
!------------------------------------------------------------------


!-------------------------------------------------------------------
!----------- version number for this module ------------------------

character(len=128)  :: version =  '$Id: rad_output_file.F90,v 1.3 2002/07/16 22:36:36 fms Exp $'
character(len=128)  :: tag     =  '$Name: havana $'


!---------------------------------------------------------------------
!-------  interfaces --------

public   &
         rad_output_file_init, write_rad_output_file,   &
         rad_output_file_end

private  register_fields


!-------------------------------------------------------------------
!-------- namelist  ---------

logical :: write_data_file=.false.  ! data file to be written  ?


namelist  / rad_output_file_nml /  &
                                  write_data_file
     

!------------------------------------------------------------------
!----public data------


!------------------------------------------------------------------
!----private data------

logical        :: do_isccp=.false.     ! isccp cloud output included  ?
logical        :: rad_output_file_initialized= .false. 
                                       ! module initialized ?

!--------------------------------------------------------------------
! netcdf diagnostics field variables
!--------------------------------------------------------------------
character(len=16), parameter :: mod_name='radiation'
real                         :: missing_value = -999.
integer                      :: id_radswp, id_radp, id_temp, id_rh2o, &
                                id_qo3, id_cmxolw, id_crndlw,   &
                                id_flxnet, id_fsw, id_ufsw, id_psj, &
                                id_tmpsfc, id_cvisrfgd, id_cirrfgd, &
                                id_tot_clds, id_cld_isccp_hi, &
                                id_cld_isccp_mid, id_cld_isccp_low, &
                                id_radswpcf, id_radpcf, id_flxnetcf, &
                                id_fswcf, id_ufswcf



!---------------------------------------------------------------------
!---------------------------------------------------------------------


                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
!#####################################################################

subroutine rad_output_file_init (axes, Time)

!--------------------------------------------------------------------
!    rad_output_file_init is the constructor for rad_output_file_mod.
!--------------------------------------------------------------------

integer, dimension(4), intent(in), optional    :: axes
type(time_type),       intent(in), optional    :: Time

!--------------------------------------------------------------------
!  intent(in), optional variables:
!
!    these variables are present when running the gcm, not present
!    when running the standalone code.
!  
!       axes      diagnostic variable axes for netcdf files
!       Time      current time [ time_type(days, seconds) ]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables

integer   :: unit, io, ierr

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (rad_output_file_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call utilities_init
      call rad_utilities_init
      if (Environment%running_gcm) then
        call diag_manager_init  
      endif
!     call time_manager_init ! doesn't exist
 
!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
          read (unit, nml=rad_output_file_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'rad_output_file_nml')
        enddo
10      call close_file (unit)
      endif

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      unit = open_file ('logfile.out', action='append')
      if (get_my_pe() == 0) then
!     if ( get_my_pe() == get_root_pe() ) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit,nml=rad_output_file_nml)
      endif
      call close_file (unit)

!--------------------------------------------------------------------
!    if running gcm, continue on if data file is to be written. 
!--------------------------------------------------------------------
      if (Environment%running_gcm) then
        if (write_data_file) then

!--------------------------------------------------------------------
!    define a variable indicating whether isccp output variables will
!    be present.
!--------------------------------------------------------------------
          if (Cldrad_control%do_rh_clouds) then
            do_isccp = .true.
          endif

!---------------------------------------------------------------------
!    register the diagnostic fields for output.
!---------------------------------------------------------------------
          call register_fields (Time, axes)
        endif

!--------------------------------------------------------------------
!    if running standalone and data file desired, write error message.
!--------------------------------------------------------------------
      else
        if (write_data_file) then
          call error_mesg ('rad_output_file_mod', &
            ' can only write output file when running in gcm', &
                                                           FATAL)
        endif
      endif

!--------------------------------------------------------------------
    rad_output_file_initialized = .true.



end subroutine rad_output_file_init



!################################################################

subroutine write_rad_output_file (is, ie, js, je, Atmos_input,  &
                                  Rad_output, &
                                  Sw_output, Lw_output, Rad_gases, &
                                  Cldrad_props, Cld_diagnostics,   &
                                  Time_diag)

!----------------------------------------------------------------
!    write_rad_output_file produces a netcdf output file containing
!    the user-specified radiation-related variables.
!----------------------------------------------------------------

integer,                      intent(in)  ::  is, ie, js, je
type(atmos_input_type),       intent(in)  ::  Atmos_input
type(rad_output_type),        intent(in)  ::  Rad_output
type(sw_output_type),         intent(in)  ::  Sw_output
type(lw_output_type),         intent(in)  ::  Lw_output
type(radiative_gases_type),   intent(in)  ::  Rad_gases
type(cldrad_properties_type), intent(in)  ::  Cldrad_props
type(cld_diagnostics_type),   intent(in)  ::  Cld_diagnostics 
type(time_type),              intent(in)  ::  Time_diag

!------------------------------------------------------------------
!  intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data 
!                        in the physics_window being integrated
!      Atmos_input       atmos_input_type variable containing atmos-
!                        pheric input data for the radiation package 
!                        on the model grid
!      Rad_output        rad_output_type variable containing radiation
!                        output data needed by other modules
!      Sw_output         sw_output_type variable containing shortwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid
!      Lw_output         lw_output_type variable containing longwave 
!                        radiation output data from the sea_esf_rad
!                        radiation package on the model grid 
!      Rad_gases         radiative_gases_type variable containing rad-
!                        iative gas input data for the radiation package
!                        on the model grid
!      Cldrad_props      cldrad_properties_type variable containing 
!                        cloud radiative property input data for the 
!                        radiation package on the model grid
!      Cld_diagnostics   cld_diagnostics_type variable containing cloud
!                        microphysical and /or isccp diagnostic fields
!      Time_diag         time on next timestep, used as stamp for diag-
!                        nostic output  [ time_type  (days, seconds) ]  
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:
!
!      tmpsfc         surface temperature [ deg K ]
!      psj            surface pressure [ hPa ]
!      cvisrfgd       surface visible light albedo [ dimensionless ]
!      cirrfgd        surface ir albedo [ dimensionless ]
!      tot_clds       total column isccp clouds [ percent ]
!      cld_isccp_hi   number of isccp high clouds [ percent ]
!      cld_isccp_mid  number of isccp middle clouds [ percent ]
!      cld_isccp_low  number of isccp low clouds [ percent ]
!      fsw            net shortwave flux [ W / m**2 ]
!      ufsw           upward shortwave flux [ W / m**2 ]
!      fswcf          net sw flux in the absence of clouds [ W / m**2 ]
!      ufswcf         upward sw flux in absence of clouds [ W / m**2]
!      flxnet         net longwave flux [ W / m**2 ]
!      flxnetcf       net lw flux in the absence of clouds [ W / m**2 ]
!      temp           temperature [ deg K ]
!      rh2o           water vapor specific humidity [ g / g ]
!      qo3            ozone mixing ratio [ g / g ]
!      heatra         lw heating rate [ deg K / day ]
!      heatracf       lw heating rate without cloud [ deg K / day ]
!      cmxolw         amount of maximal overlap clouds [ percent]
!      crndlw         amount of ramndom overlap clouds [ percent]
!      radp           lw + sw heating rate [ deg K / sec ]
!      radswp         sw heating rate [ deg K / sec ]
!      radpcf         lw + sw heating rate w/o clouds [ deg K / sec ]
!      radswpcf       sw heating rate w/o clouds [ deg K / sec ]
!
!---------------------------------------------------------------------
      real, dimension (size(Atmos_input%press,1),   &
                       size(Atmos_input%press,2) )  ::  &
                    tmpsfc, psj, cvisrfgd, cirrfgd, tot_clds,   &
                    cld_isccp_hi, cld_isccp_mid, cld_isccp_low

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3)  ) ::   &
                    fsw, ufsw, fswcf, ufswcf, flxnet, flxnetcf

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
                       size(Atmos_input%press,3)-1 ) ::  &
                    temp, rh2o, qo3, cmxolw, crndlw, radp, radswp, &
                    radpcf, radswpcf

      logical   :: used  
      integer   :: kerad ! number of model layers

!--------------------------------------------------------------------
!    if the file is not to be written, do nothing.e desired fields from
!    the input derived data types.
!--------------------------------------------------------------------
      if (write_data_file) then

!--------------------------------------------------------------------
!    if the file is to be written, define the number of model layers.
!--------------------------------------------------------------------
        kerad = ubound(Atmos_input%temp,3) - 1

!--------------------------------------------------------------------
!    retrieve the desired fields from the input derived data types.
!--------------------------------------------------------------------
        tmpsfc(:,:)     = Atmos_input%temp(:,:,kerad+1)
        psj   (:,:)     = 0.01*Atmos_input%press(:,:,kerad+1)
        temp(:,:,:)     = Atmos_input%temp(:,:,1:kerad)
        rh2o(:,:,:)     = Atmos_input%rh2o(:,:,:)
        radp(:,:,:)     = Rad_output%tdt_rad(:,js:je,:)
        radswp(:,:,:)   = Rad_output%tdtsw  (:,js:je,:)
        cirrfgd(:,:)    = Atmos_input%asfc(:,:)
        cvisrfgd(:,:)   = Atmos_input%asfc(:,:)
        fsw(:,:,:)      = Sw_output% fsw(:,:,:)
        ufsw(:,:,:)     = Sw_output%ufsw(:,:,:)
        flxnet(:,:,:)   = Lw_output%flxnet(:,:,:)
        qo3(:,:,:)      = Rad_gases%qo3(:,:,:)
        cmxolw(:,:,:)   = 100.0*Cldrad_props%cmxolw(:,:,:)
        crndlw(:,:,:)   = 100.0*Cldrad_props%crndlw(:,:,:)

        if (Rad_control%do_totcld_forcing) then 
          fswcf(:,:,:)    = Sw_output%fswcf(:,:,:)
          ufswcf(:,:,:)   = Sw_output%ufswcf(:,:,:)
          flxnetcf(:,:,:) = Lw_output%flxnetcf(:,:,:)
          radpcf(:,:,:)   = Rad_output%tdt_rad_clr(:,js:je,:)
          radswpcf(:,:,:) = Rad_output%tdtsw_clr  (:,js:je,:)
        endif
        if (do_isccp) then
          tot_clds(:,:)      = 100.0*Cld_diagnostics%tot_clds(:,:)
          cld_isccp_hi (:,:) = 100.0*Cld_diagnostics%cld_isccp_hi (:,:)
          cld_isccp_mid(:,:) = 100.0*Cld_diagnostics%cld_isccp_mid(:,:)
          cld_isccp_low(:,:) = 100.0*Cld_diagnostics%cld_isccp_low(:,:)
        endif

!---------------------------------------------------------------------
!    send the user-designated data to diag_manager_mod for processing.
!---------------------------------------------------------------------
        if (id_radswp > 0 ) then
          used = send_data (id_radswp, radswp, Time_diag, is, js, 1)
        endif

        if (id_radp > 0 ) then
          used = send_data (id_radp, radp, Time_diag, is, js, 1)
        endif

        if (id_temp > 0 ) then
          used = send_data (id_temp, temp, Time_diag, is, js, 1)
        endif

        if (id_rh2o > 0 ) then
          used = send_data (id_rh2o, rh2o, Time_diag, is, js, 1)
        endif

        if (id_qo3  > 0 ) then
          used = send_data (id_qo3, qo3, Time_diag, is, js, 1)
        endif

        if (id_cmxolw > 0 ) then
          used = send_data (id_cmxolw, cmxolw, Time_diag, is, js, 1)
        endif

        if (id_crndlw > 0 ) then
          used = send_data (id_crndlw, crndlw, Time_diag, is, js, 1)
        endif

        if (id_flxnet > 0 ) then
          used = send_data (id_flxnet, flxnet, Time_diag, is, js, 1)
        endif

        if (id_fsw    > 0 ) then
          used = send_data (id_fsw   , fsw   , Time_diag, is, js, 1)
        endif

        if (id_ufsw   > 0 ) then
          used = send_data (id_ufsw  , ufsw  , Time_diag, is, js, 1)
        endif

        if (id_psj    > 0 ) then
          used = send_data (id_psj   , psj   , Time_diag, is, js)
        endif

        if (id_tmpsfc > 0 ) then
          used = send_data (id_tmpsfc, tmpsfc, Time_diag, is, js)
        endif

        if (id_cvisrfgd > 0 ) then
          used = send_data (id_cvisrfgd, cvisrfgd, Time_diag, is, js)
        endif

        if (id_cirrfgd > 0 ) then
          used = send_data (id_cirrfgd , cirrfgd , Time_diag, is, js)
        endif

        if (do_isccp) then
          if (id_tot_clds > 0 ) then
            used = send_data (id_tot_clds , tot_clds, Time_diag, is, js)
          endif

          if (id_cld_isccp_hi > 0 ) then
            used = send_data (id_cld_isccp_hi, cld_isccp_hi,   &
                              Time_diag, is, js)
          endif

          if (id_cld_isccp_mid > 0 ) then
            used = send_data (id_cld_isccp_mid, cld_isccp_mid,   &
                              Time_diag, is, js)
          endif

          if (id_cld_isccp_low > 0 ) then
            used = send_data (id_cld_isccp_low, cld_isccp_low,   &
                              Time_diag, is, js)
          endif
        endif

        if (Rad_control%do_totcld_forcing) then
          if (id_radswpcf > 0 ) then
            used = send_data (id_radswpcf, radswpcf, Time_diag,   &
                              is, js, 1)
          endif

          if (id_radpcf > 0 ) then
            used = send_data (id_radpcf, radpcf, Time_diag, is, js, 1)
          endif

          if (id_flxnetcf > 0 ) then
            used = send_data (id_flxnetcf, flxnetcf, Time_diag,   &
                              is, js, 1)
          endif

          if (id_fswcf  > 0 ) then
            used = send_data (id_fswcf , fswcf , Time_diag, is, js, 1)
          endif

          if (id_ufswcf  > 0 ) then
            used = send_data (id_ufswcf , ufswcf , Time_diag, is, js, 1)
          endif
        endif
      endif


!------------------------------------------------------------------


end subroutine write_rad_output_file



!#####################################################################

subroutine rad_output_file_end

!-------------------------------------------------------------------
!    rad_output_file_end is the destructor for rad_output_file_mod.
!-------------------------------------------------------------------

      rad_output_file_initialized= .false. 


end subroutine rad_output_file_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!#################################################################

subroutine register_fields (Time, axes)

!--------------------------------------------------------------------
!    register_fields send the relevant information concerning the 
!    user-desired output fields to diag_manager_mod.
!--------------------------------------------------------------------

type(time_type),       intent(in)          :: Time
integer, dimension(4), intent(in)          :: axes

!--------------------------------------------------------------------
!  intent(in) variables:
!
!       Time      current time [ time_type(days, seconds) ]
!       axes      diagnostic variable axes for netcdf files
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:
!
!       bxes      diagnostic variable axes with elements (1:3) valid 
!                 for variables defined at flux levels
!
!--------------------------------------------------------------------
      integer, dimension(4)    :: bxes
 
!-------------------------------------------------------------------
!    define variable axis array with elements (1:3) valid for variables
!    defined at flux levels.
!-------------------------------------------------------------------
      bxes(1:2) = axes(1:2)
      bxes(3) = axes(4)
      bxes(4) = axes(4)

!---------------------------------------------------------------------
!    register the potential diagnostic variables from this module.
!    the variables will actually be saved and then output only if
!    they are activated through the input diag_file.
!---------------------------------------------------------------------
      id_radswp = &
         register_diag_field (mod_name, 'radswp', axes(1:3), Time, &
                          'temperature tendency for SW radiation', &
                          'deg_K/sec', missing_value=missing_value)

      id_radp = &
         register_diag_field (mod_name, 'radp', axes(1:3), Time, &
                          'temperature tendency for radiation', &
                          'deg_K/sec', missing_value=missing_value)

      id_temp   = &
         register_diag_field (mod_name, 'temp', axes(1:3), Time, &
                          'temperature field', &
                          'deg_K', missing_value=missing_value)

      id_rh2o   = &
         register_diag_field (mod_name, 'rh2o', axes(1:3), Time, &
                          'water vapor mixing ratio', &
                          'kg/kg', missing_value=missing_value)

      id_qo3    = &
         register_diag_field (mod_name, 'qo3', axes(1:3), Time, &
                          'ozone mixing ratio', &
                          'kg/kg', missing_value=missing_value)

      id_cmxolw = &
         register_diag_field (mod_name, 'cmxolw', axes(1:3), Time, &
                          'maximum overlap cloud amount', &
                          'percent', missing_value=missing_value)

      id_crndlw = &
         register_diag_field (mod_name, 'crndlw', axes(1:3), Time, &
                          'random overlap cloud amount', &
                          'percent', missing_value=missing_value)

      id_flxnet = &
         register_diag_field (mod_name, 'flxnet', bxes(1:3), Time, &
                          'net longwave radiative flux', &
                          'W/m**2', missing_value=missing_value)

      id_fsw    = &
         register_diag_field (mod_name, 'fsw', bxes(1:3), Time, &
                          'net shortwave radiative flux', &
                          'W/m**2', missing_value=missing_value)

      id_ufsw   = &
         register_diag_field (mod_name, 'ufsw', bxes(1:3), Time, &
                          'upward shortwave radiative flux ', &
                          'W/m**2', missing_value=missing_value)

      id_psj    = &
         register_diag_field (mod_name, 'psj', axes(1:2), Time, &
                          'surface pressure', &
                          'hPa', missing_value=missing_value)

      id_tmpsfc = &
         register_diag_field (mod_name, 'tmpsfc', axes(1:2), Time, &
                          'surface temperature', &
                          'deg_K', missing_value=missing_value)

      id_cvisrfgd = &
         register_diag_field (mod_name, 'cvisrfgd', axes(1:2), Time, &
                          'visible surface albedo', &
                          'dimensionless', missing_value=missing_value)

      id_cirrfgd = &
         register_diag_field (mod_name, 'cirrfgd', axes(1:2), Time, &
                          'visible infra-red albedo', &
                          'dimensionless', missing_value=missing_value)

      if (do_isccp) then
        id_tot_clds = &
           register_diag_field (mod_name, 'tot_clds', axes(1:2), Time, &
                            'total isccp cloud cover', 'percent',  &
                            missing_value=missing_value)

        id_cld_isccp_hi = &
           register_diag_field (mod_name, 'cld_isccp_hi', axes(1:2),   &
                            Time, 'isccp high clouds', 'percent',&
                            missing_value=missing_value)

        id_cld_isccp_mid = &
           register_diag_field (mod_name, 'cld_isccp_mid', axes(1:2), &
                            Time, 'isccp middle clouds', 'percent',   &
                            missing_value=missing_value)

        id_cld_isccp_low = &
           register_diag_field (mod_name, 'cld_isccp_low', axes(1:2), &
                            Time, 'isccp low clouds', 'percent',   &
                            missing_value=missing_value)

      endif 

      if (Rad_control%do_totcld_forcing) then
        id_radswpcf = &
           register_diag_field (mod_name, 'radswpcf', axes(1:3), Time, &
                            'temperature forcing from sw w/o clouds', &
                            'deg_K/sec', missing_value=missing_value)

        id_radpcf = &
           register_diag_field (mod_name, 'radpcf', axes(1:3), Time, &
                            'temperature forcing w/o clouds', &
                            'deg_K/sec', missing_value=missing_value)

        id_flxnetcf = &
           register_diag_field (mod_name, 'flxnetcf', bxes(1:3), Time, &
                            'net longwave flux w/o clouds', &
                            'W/m**2', missing_value=missing_value)

        id_fswcf = &
           register_diag_field (mod_name, 'fswcf', bxes(1:3), Time, &
                            'net shortwave flux w/o clouds', &
                            'W/m**2', missing_value=missing_value)

        id_ufswcf   = &
           register_diag_field (mod_name, 'ufswcf', bxes(1:3), Time, &
                            'upward shortwave flux w/o clouds', &
                            'W/m**2', missing_value=missing_value)

      endif

!---------------------------------------------------------------------


end subroutine register_fields


!####################################################################





                  end module rad_output_file_mod
