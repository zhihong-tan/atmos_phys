module atmos_sulfate_mod
! <DESCRIPTION>
!   This module is an implementation of sulfate chemistry. It contains
!   tracer emissions and chemistry. The chemistry is partly based on MOZART.
!   The change of concentration of SO2, DMS, SO4, MSA and H2O2 are
!   calculated using monthly mean concentration of OH, HO2, jH2O2, NO3, O3,
!   pH. The emissions include:
!     - DMS from seawater
!     - SO2 by fossil fuel, biomass burning, non-eruptive volcanoes and aircraft
!     - SO4 by fossil fuel
! </DESCRIPTION>
! <WARNING>
!  To save space only the actual month of input files are kept in memory.
!  This implies that the "atmos_SOx_init" should be executed at the begining
!  of each month. In other words, the script should not run more than 1 month
!  without a restart.
! </WARNING>
! <CONTACT EMAIL="Paul.Ginouxe@noaa.gov">
!   Paul Ginoux
! </CONTACT>
!-----------------------------------------------------------------------

use                    fms_mod, only : file_exist,              &
                                       write_version_number,    &
                                       mpp_pe,                  &
                                       mpp_root_pE,             &
                                       close_file,              &
                                       stdlog,                  &
                                       check_nml_error, error_mesg, &
                                       open_namelist_file, FATAL

use           time_manager_mod, only : time_type, get_date_julian
use           diag_manager_mod, only : send_data,               &
                                       register_diag_field,     &
                                       register_static_field
use         tracer_manager_mod, only : get_tracer_index,        &
                                       set_tracer_atts
use          field_manager_mod, only : MODEL_ATMOS
use              constants_mod, only : PI, GRAV, RDGAS, WTMAIR
use atmos_tracer_utilities_mod, only : interp_emiss

implicit none

private
!-----------------------------------------------------------------------
!----- interfaces -------
!
public  atmos_SOx_init, atmos_SOx_end, &
        atmos_DMS_emission, atmos_SO2_emission, atmos_SO4_emission, &
        SOx_source_input, chem_sox, &
        get_SO2_nerup_volc_emis

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------

!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

integer :: nSO4 = 0  ! tracer number for Sulfate               = SO4=
integer :: nDMS = 0  ! tracer number for Dimethyl sulfide      = CH3SCH3
integer :: nSO2 = 0  ! tracer number for Sulfur dioxide        = SO2
integer :: nMSA = 0  ! tracer number for Methane sulfonic acid = CH3SO3H
integer :: nH2O2= 0  ! tracer number for Hydrogen peroxyde     = H2O2

real , parameter :: WTM_S     = 32.06599807
real , parameter :: WTM_SO2   = 64.06480408
real , parameter :: WTM_SO4   = 96.06359863
real , parameter :: WTM_DMS   = 62.13240051
real , parameter :: WTM_MSA   = 96.06359863

!--- identification numbers for  diagnostic fields and axes ----
integer ::   id_OH_conc             = 0
integer ::   id_HO2_conc            = 0
integer ::   id_NO3_conc            = 0
integer ::   id_OH_diurnal          = 0
integer ::   id_HO2_diurnal         = 0
integer ::   id_jH2O2_diurnal       = 0
integer ::   id_NO3_diurnal         = 0
integer ::   id_O3_vmr              = 0
integer ::   id_pH                  = 0
integer ::   id_jH2O2               = 0

integer ::   id_pph =0
integer ::   id_po3 =0
integer ::   id_DMSo                = 0
integer ::   id_sstemp              = 0
integer ::   id_SO2_anth_l1_emis    = 0
integer ::   id_SO2_anth_l2_emis    = 0
integer ::   id_SO4_anth_l1_emis    = 0
integer ::   id_SO4_anth_l2_emis    = 0
integer ::   id_SO2_aircraft_emis   = 0
integer ::   id_SO2_bioburn_emis    = 0
integer ::   id_SO2_nerup_volc_emis = 0
integer ::   id_DMS_emis            = 0
integer ::   id_SO2_emis            = 0
integer ::   id_SO4_emis            = 0
integer ::   id_DMS_chem            = 0
integer ::   id_SO2_chem            = 0
integer ::   id_SO4_chem            = 0
integer ::   id_MSA_chem            = 0
integer ::   id_H2O2_chem           = 0

!--- Arrays of SOx emission rates and gas species needed for SOx chemistry
real, allocatable, dimension(:,:)   :: DMSo        ![nM/L]
real, allocatable, dimension(:,:)   :: SO4_anth_l1_emis
real, allocatable, dimension(:,:)   :: SO4_anth_l2_emis
real, allocatable, dimension(:,:)   :: SO2_anth_l1_emis
real, allocatable, dimension(:,:)   :: SO2_anth_l2_emis
real, allocatable, dimension(:,:)   :: SO2_bioburn_emis
real, allocatable, dimension(:,:,:) :: SO2_aircraft_emis
real, allocatable, dimension(:,:,:) :: OH_conc
real, allocatable, dimension(:,:,:) :: HO2_conc
real, allocatable, dimension(:,:,:) :: jH2O2
real, allocatable, dimension(:,:,:) :: NO3_conc
real, allocatable, dimension(:,:,:) :: O3_vmr
real, allocatable, dimension(:,:,:) :: pH

!trim(runtype) 
!biomass_only; fossil_fuels_only, natural_only, anthrop
character(len=20):: runtype = "anthrop"

namelist / aerosol_emissions_nml/  &
           runtype                         

logical :: module_is_initialized=.FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id: atmos_sulfate.F90,v 13.0 2006/03/28 21:15:41 fms Exp $'
character(len=128) :: tagname = '$Name: memphis_2006_08 $'
!-----------------------------------------------------------------------

contains


!#######################################################################

!<SUBROUTINE NAME="atmos_SOx_init">
!<OVERVIEW>
! The constructor routine for the sulfate module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to initialize the sulfate module.
!</DESCRIPTION>
 subroutine atmos_SOx_init ( lonb, latb, nlev, axes, Time, mask)

!-----------------------------------------------------------------------
real, intent(in),    dimension(:)                   :: lonb, latb
integer, intent(in)                                 :: nlev
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
real, intent(in), dimension(:,:,:), optional        :: mask
character(len=7), parameter :: mod_name = 'tracers'
logical :: flag
integer :: n, m, nsulfate
!
!----------------------------------------------------------------------
!  local variables:

      integer   ::   ierr

!---------------------------------------------------------------------
!    local variables:
!
!         unit       io unit number used to read namelist file
!         ierr       error code
!         io         error status returned from io operation
!         n          do-loop index
!
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!
      integer  log_unit,unit,io,index,ntr,nt
      real :: initial_values(5)
      character*3 :: SOx_tracer(5)
!
!     1. DMS       = Dimethyl sulfide            = CH3SCH3
!     2. SO2       = Sulfur dioxide              = SO2     
!     3. SO4       = Sulfate                     = SO4=   
!     4. MSA       = Methane sulfonic acid       = CH3SO3H
!     5. H2O2      = Hydrogen peroxyde           = H2O2
!                                                                      
      data SOx_tracer/'simpleDMS', &
                      'simpleSO2', &
                      'simpleSO4', &
                      'simpleMSA', &
                      'simpleH2O2' /
      
      if (module_is_initialized) return


!---- write namelist ------------------

      call write_version_number (version, tagname)

!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosol_emissions_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'aerosol_emissions_nml')
        end do
10      call close_file (unit)
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
                          write (stdlog(), nml=aerosol_emissions_nml)


!----- set initial value of sulfate ------------

     do m=1,size(SOx_tracer)

       n = get_tracer_index(MODEL_ATMOS,SOx_tracer(m))
       if (n>0) then
         nsulfate=n
!        call set_tracer_atts(MODEL_ATMOS,SOx_tracer(m),SOx_tracer(m),'mmr')
         if (nsulfate > 0 .and. mpp_pe() == mpp_root_pe()) &
                 write (stdlog(),30) SOx_tracer(m),nsulfate
       endif
     enddo


  30   format (A,' was initialized as tracer number ',i2)
! Register a diagnostic field : emission of SOx species
      id_DMS_emis   = register_diag_field ( mod_name,             &
                      'simpleDMS_emis', axes(1:2),Time,                 &
                      'simpleDMS_emis', 'kgS/m2/s',                     &
                       missing_value=-999.  )
      id_SO2_emis   = register_diag_field ( mod_name,             &
                      'simpleSO2_emis', axes(1:3),Time,                 &
                      'simpleSO2_emis', 'kgS/m2/s',                     &
                       missing_value=-999.  )
      id_SO4_emis   = register_diag_field ( mod_name,             &
                      'simpleSO4_emis', axes(1:3),Time,                 &
                      'simpleSO4_emis', 'kgS/m2/s',                     &
                       missing_value=-999.  )
      if (.not.allocated(DMSo)) then
        allocate (DMSo(size(lonb)-1,size(latb)-1)                   )
        id_DMSo       = register_diag_field ( mod_name,             &
                      'DMSo',axes(1:2),Time,                      &
                      'Dimethylsulfide seawater concentration',   &
                      'nM/L')
      endif

      id_pph = register_diag_field ( mod_name,             &
                      'pph',axes(1:3),Time,                    &
                      'pH in chem',                  &
                      'none')

      id_po3 = register_diag_field ( mod_name,             &
                      'po3',axes(1:3),Time,                    &
                      'o3 in chem',                  &
                      'none')

      if (.not.allocated(SO2_anth_l1_emis)) then
        allocate (SO2_anth_l1_emis(size(lonb)-1,size(latb)-1)    )
        id_SO2_anth_l1_emis = register_diag_field ( mod_name,     &
                      'simpleSO2_anth_l1_emis',axes(1:2),Time,          &
                      'SO2 anthropogenic emission GEIA-85 level1',&
                      'kgS/m2/s')
      endif
      if (.not.allocated(SO2_anth_l2_emis)) then
        allocate (SO2_anth_l2_emis(size(lonb)-1,size(latb)-1)    )
        id_SO2_anth_l2_emis = register_diag_field ( mod_name,     &
                      'simpleSO2_anth_l2_emis',axes(1:2),Time,          &
                      'SO2 anthropogenic emission GEIA-85 level2',&
                      'kgS/m2/s')
      endif
      if (.not.allocated(SO2_aircraft_emis)) then
        allocate (SO2_aircraft_emis(size(lonb)-1,size(latb)-1,nlev))
        id_SO2_aircraft_emis= register_diag_field ( mod_name,     &
                      'simpleSO2_aircraft_emis',axes(1:3),Time,         &
                      'SO2 emission by aircraft',                 &
                      'kgS/m2/s')
      endif
      if (.not.allocated(SO2_bioburn_emis)) then
        allocate (SO2_bioburn_emis(size(lonb)-1,size(latb)-1)    )
        id_SO2_bioburn_emis = register_diag_field ( mod_name,     &
                      'simpleSO2_bioburn_emis',axes(1:2),Time,          &
                      'SO2 emission from biomass burning',        &
                      'kgS/m2/s')
      endif
      id_SO2_nerup_volc_emis = register_diag_field ( mod_name,     &
                      'simpleSO2_nerup_volc_emis',axes(1:3),Time,     &
                      'SO2 emission from non-eruptive volcanoes', &
                      'kgS/m2/s')
      if (.not.allocated(SO4_anth_l1_emis)) then
        allocate (SO4_anth_l1_emis(size(lonb)-1,size(latb)-1)    )
        id_SO4_anth_l1_emis = register_diag_field ( mod_name,     &
                      'simpleSO4_anth_l1_emis',axes(1:2),Time,          &
                      'SO4 anthropogenic emission GEIA-85 level1',&
                      'kgS/m2/s')
      endif
      if (.not.allocated(SO4_anth_l2_emis)) then
        allocate (SO4_anth_l2_emis(size(lonb)-1,size(latb)-1)    )
        id_SO4_anth_l2_emis = register_diag_field ( mod_name,     &
                      'simpleSO4_anth_l2_emis',axes(1:2),Time,          &
                      'SO4 anthropogenic emission GEIA-85 level2',&
                      'kgS/m2/s')
      endif
      if (.not.allocated(NO3_conc) ) &
        allocate (NO3_conc(size(lonb)-1,size(latb)-1,nlev) )
      if (id_NO3_conc.eq.0) &
        id_NO3_conc   = register_diag_field ( mod_name,           &
                      'simpleNO3_conc',axes(1:3),Time,                       &
                      'NO3 concentration',                        &
                      'molec.cm-3')
      if (id_NO3_diurnal.eq.0) &
        id_NO3_diurnal= register_diag_field ( mod_name,           &
                      'simpleNO3_diurnal',axes(1:3),Time,                       &
                      'Time varying NO3 concentration',                        &
                      'molec.cm-3')

      if (.not.allocated(OH_conc) ) &
        allocate (OH_conc(size(lonb)-1,size(latb)-1,nlev)  )
      if (id_OH_conc.eq.0) &
        id_OH_conc    = register_diag_field ( mod_name,           &
                      'simpleOH_SOx_conc',axes(1:3),Time,                        &
                      'Hydroxyl radical concentration',           &
                      'molec.cm-3')
      if (id_OH_diurnal.eq.0) &
        id_OH_diurnal = register_diag_field ( mod_name,           &
                      'simpleOH_SOx_diurnal',axes(1:3),Time,                        &
                      'Varying Hydroxyl radical concentration',           &
                      'molec.cm-3')


      if (.not.allocated(jH2O2) ) &
        allocate (jH2O2(size(lonb)-1,size(latb)-1,nlev)  )
      if (id_jH2O2.eq.0) &
        id_jH2O2   = register_diag_field ( mod_name,           &
                     'simplejH2O2',axes(1:3),Time,              &
                     'H2O2 photodissociation',           &
                     's-1')
      if (id_jH2O2_diurnal.eq.0) &
        id_jH2O2_diurnal = register_diag_field ( mod_name,           &
                      'simplejH2O2_diurnal',axes(1:3),Time,               &
                      'Varying H2O2 photodissociation',   &
                      's-1')

      if (.not.allocated(HO2_conc) ) &
        allocate (HO2_conc(size(lonb)-1,size(latb)-1,nlev)  )
      if (id_HO2_conc.eq.0) &
        id_HO2_conc   = register_diag_field ( mod_name,           &
                      'simpleHO2_SOx_conc',axes(1:3),Time,              &
                      'Hydroperoxyl radical concentration',           &
                      'molec.cm-3')
      if (id_HO2_diurnal.eq.0) &
        id_HO2_diurnal = register_diag_field ( mod_name,           &
                      'simpleHO2_SOx_diurnal',axes(1:3),Time,               &
                      'Varying Hydroperoxyl radical concentration',   &
                      'molec.cm-3')

      if (.not.allocated(O3_vmr) ) &
        allocate (O3_vmr(size(lonb)-1,size(latb)-1,nlev) )
      if (id_O3_vmr.eq.0) &
        id_O3_vmr        = register_diag_field ( mod_name,           &
                      'simpleO3_vmr',axes(1:3),Time,                       &
                      'Ozone volume mixing ratio',               &
                      'vmr')

      if (.not.allocated(pH) ) &
        allocate (pH(size(lonb)-1,size(latb)-1,nlev) )
      if (id_pH.eq.0) &
        id_pH        = register_diag_field ( mod_name,           &
                      'simplepH',axes(1:3),Time,                       &
                      'pH',      &
                      'none')

      if (id_DMS_chem.eq.0) &
      id_DMS_chem   = register_diag_field ( mod_name,           &
                      'simpleDMS_chem',axes(1:3),Time,                       &
                      'DMS chemical production',      &
                      'kgS/m2/s')

      if (id_SO2_chem.eq.0) &
      id_SO2_chem   = register_diag_field ( mod_name,           &
                      'simpleSO2_chem',axes(1:3),Time,                       &
                      'SO2 chemical production',      &
                      'kgS/m2/s')

      if (id_SO4_chem.eq.0) &
      id_SO4_chem   = register_diag_field ( mod_name,           &
                      'simpleSO4_chem',axes(1:3),Time,                       &
                      'SO4 chemical production',      &
                      'kgS/m2/s')

      if (id_MSA_chem.eq.0) &
      id_MSA_chem   = register_diag_field ( mod_name,           &
                      'simpleMSA_chem',axes(1:3),Time,                       &
                      'MSA chemical production',      &
                      'kgS/m2/s')

      if (id_H2O2_chem.eq.0) &
      id_H2O2_chem   = register_diag_field ( mod_name,           &
                      'simpleH2O2_chem',axes(1:3),Time,                       &
                      'H2O2 chemical production',      &
                      'kgH2O2/m2/s')

      call write_version_number (version, tagname)

      module_is_initialized = .TRUE.

!-----------------------------------------------------------------------
 end subroutine atmos_SOx_init
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_SOx_end">
!<OVERVIEW>
!  The destructor routine for the sulfate module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_SOx_end
!</TEMPLATE>
 subroutine atmos_SOx_end

        if (allocated(DMSo))                deallocate(DMSo)
        if (allocated(SO4_anth_l1_emis))    deallocate(SO4_anth_l1_emis)
        if (allocated(SO4_anth_l2_emis))    deallocate(SO4_anth_l2_emis)
        if (allocated(SO2_anth_l1_emis))    deallocate(SO2_anth_l1_emis)
        if (allocated(SO2_anth_l2_emis))    deallocate(SO2_anth_l2_emis)
        if (allocated(SO2_bioburn_emis))    deallocate(SO2_bioburn_emis)
        if (allocated(SO2_aircraft_emis))   deallocate(SO2_aircraft_emis)
        if (allocated(OH_conc))             deallocate(OH_conc)
        if (allocated(HO2_conc))            deallocate(HO2_conc)
        if (allocated(jH2O2))               deallocate(jH2O2)
        if (allocated(NO3_conc))            deallocate(NO3_conc)
        if (allocated(O3_vmr))              deallocate(O3_vmr)
        if (allocated(pH))                  deallocate(pH)
 
        module_is_initialized = .FALSE.

 end subroutine atmos_SOx_end
!</SUBROUTINE>
!#######################################################################
!<SUBROUTINE NAME="SOx_source_input">
!<OVERVIEW>
!  This subroutine read the monthly mean concentrations of 
!  OH, HO2, NO3, O3, and the monthly photodissociation rates jH2o2 and
!  pH, as well as the emissions for DMS, SO2, and SO4
!  *****WARNING:
!  To save space only the actual month is kept in memory which implies
!  that the "atmos_SOx_init" should be executed at the begining of each
!  month. In other words, the script should not run more than 1 month
!  without a restart
!</OVERVIEW>
 subroutine SOx_source_input( pfull, lon, lat, imonth, Time, is, ie, js, je, &
                              kbot)
!--- Input variables
        real, intent(in), dimension(:,:,:) :: pfull
        real, dimension(:,:), intent(in)   :: lon, lat
        integer, intent(in)                :: imonth
        type(time_type),intent(in)         :: Time
        integer, intent(in)                :: is, ie, js, je
        integer, intent(in), dimension(:,:), optional :: kbot
!--- Working variables
        real                        :: dtr, lat_S, lon_W, dlon, dlat
        real                        :: variable
        integer                     :: i, j, k, kd, im, unit, ios
        logical                     :: opened
        real, dimension(144, 90)    :: data2D
        real, dimension(144, 90,24) :: data3D
        character (len=3)           :: month(12)

!--- Input filenames
        character (len=8 ) :: FNMDMS    = 'DMS_AM2_'
        character (len=8 ) :: FNMNO3    = 'NO3_AM2_'
        character (len=7 ) :: FNMOH     = 'OH_AM2_'
        character (len=8 ) :: FNMHO2    = 'HO2_AM2_'
        character (len=10) :: FNMJH2O2  = 'jH2O2_AM2_'
        character (len=7 ) :: FNMO3     = 'O3_AM2_'
        character (len=7 ) :: FNMPH     = 'PH_AM2_'
        character (len=8 ) :: FNMLAI    = 'LAI_AM2_'
        character (len=17) :: FNMSO2a   = 'SO2_aircraft_AM2_'
        character (len=11) :: FNMBB     = 'biobur_AM2_'
        character (len=14) :: FNMSO2an1 = 'SO2_GEIA1_AM2_'
        character (len=14) :: FNMSO2an2 = 'SO2_GEIA2_AM2_'
        character (len=14) :: FNMSO4an1 = 'SO4_GEIA1_AM2_'
        character (len=14) :: FNMSO4an2 = 'SO4_GEIA2_AM2_'

        data month/'jan','feb','mar','apr','may','jun','jul', &
                   'aug','sep','oct','nov','dec'/
!---
        dtr = PI/180.
        kd=size(pfull,3)
!
!------------------------------------------------------------------------
! Read 12 monthly average DMS seawater concentrations [nM/L] from Andreae
! on a 2x2.5 grid
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
   write(stdlog(),*) 'Reading DMS seawater concentration, month=',imonth
!
        data2D(:,:)=0.
        DMSo(:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo
!
        open (unit,file='INPUT/'//FNMDMS//month(imonth)//'.txt', &
              form='formatted', action='read')
        do
          read (unit, '(2i4,e12.4)',iostat=ios) i,j,variable
          if ( ios.ne.0 ) exit
          data2D(i,j)=variable   ![nM/L]
        enddo

        call close_file (unit)
!------------------------------------------------------------------------
! End reading DMS seawater concentration
!------------------------------------------------------------------------

        lat_S= -90.*dtr; lon_W= -180.*dtr
        dlon = 2.5*dtr; dlat = 2.*dtr;
! --- Interpolate data
        call interp_emiss ( data2D, lon_W, lat_S, dlon, dlat, DMSo )

! --- Send the DMS data to the diag_manager for output.
        if (id_DMSo > 0 ) &
          used = send_data ( id_DMSo, DMSo, Time, is_in=is, js_in=js )
!------------------------------------------------------------------------
!    Read anthropogenic SO4 emission lower level (GEIA 1985).
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading GEIA anthropogenic SO4 emission level 1, month=',imonth
!
        SO4_anth_l1_emis(:,:) = 0.  ![kgSO4/m2/s]
        data2D(:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo
!
        open (unit,file='INPUT/'//FNMSO4an1//month(imonth)//'.txt', &
              form='formatted', action='read')
        do
          read (unit, '(2i4,e12.4)',iostat=ios) i,j,variable
          if ( ios.ne.0 ) exit
          data2D(i,j)=variable
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading SO4 emission lower level (GEIA 1985)
!------------------------------------------------------------------------
! --- Interpolate data
        call interp_emiss ( data2D, lon_W, lat_S, dlon, dlat, SO4_anth_l1_emis )
! --- Send the anthropogenic emission of SO2 to the diag_manager for output.
        if (id_SO4_anth_l1_emis > 0 ) &
          used = send_data ( id_SO4_anth_l1_emis, &
                             SO4_anth_l1_emis*wtm_S/wtm_SO4, Time, &
                              is_in=is, js_in=js )
!------------------------------------------------------------------------
!    Read anthropogenic SO4 emission upper level (GEIA 1985).
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading GEIA anthropogenic SO4 emission level 2, month=',imonth
!
        SO4_anth_l2_emis(:,:) = 0.  ![kgSO4/m2/s]
        data2D(:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo
!
        open (unit,file='INPUT/'//FNMSO4an2//month(imonth)//'.txt', &
              form='formatted', action='read')
        do
          read (unit, '(2i4,e12.4)',iostat=ios) i,j,variable
          if ( ios.ne.0 ) exit
          data2D(i,j)=variable
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading SO4 emission upper level (GEIA 1985)
!------------------------------------------------------------------------
! --- Interpolate data
        call interp_emiss ( data2D, lon_W, lat_S, dlon, dlat, SO4_anth_l2_emis )
! --- Send the anthropogenic emission of SO2 to the diag_manager for output.
        if (id_SO4_anth_l2_emis > 0 ) &
          used = send_data ( id_SO4_anth_l2_emis, &
                             SO4_anth_l2_emis*wtm_S/wtm_SO4, Time, &
                              is_in=is, js_in=js )
!------------------------------------------------------------------------
!    Read anthropogenic SO2 emission lower level (GEIA 1985).
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading GEIA anthropogenic SO2 emission level 1, month=',imonth
!
        SO2_anth_l1_emis(:,:) = 0.  ![kgSO2/grid box/s]
        data2D(:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo
!
        open (unit,file='INPUT/'//FNMSO2an1//month(imonth)//'.txt', &
              form='formatted', action='read')
        do
          read (unit, '(2i4,e12.4)',iostat=ios) i,j,variable
          if ( ios.ne.0 ) exit
          data2D(i,j)=variable
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading SO2 emission lower level (GEIA 1985)
!------------------------------------------------------------------------
! --- Interpolate data
        call interp_emiss ( data2D, lon_W, lat_S, dlon, dlat, SO2_anth_l1_emis )
! --- Send the anthropogenic emission of SO2 to the diag_manager for output.
        if (id_SO2_anth_l1_emis > 0 ) &
          used = send_data ( id_SO2_anth_l1_emis, &
                             SO2_anth_l1_emis*wtm_S/wtm_SO2, Time, &
                              is_in=is, js_in=js )
!------------------------------------------------------------------------
!    Read anthropogenic SO2 emission upper level (GEIA 1985).
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading GEIA anthropogenic SO2 emission level 2, month=',imonth
!
        SO2_anth_l2_emis(:,:) = 0.  ![kgSO2/m2/s]
        data2D(:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo
!
        open (unit,file='INPUT/'//FNMSO2an2//month(imonth)//'.txt', &
              form='formatted', action='read')
        do
          read (unit, '(2i4,e12.4)',iostat=ios) i,j,variable
          data2D(i,j)=variable
          if ( ios.ne.0 ) exit
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading SO2 emission upper level (GEIA 1985)
!------------------------------------------------------------------------
! --- Interpolate data
        call interp_emiss ( data2D, lon_W, lat_S, dlon, dlat, SO2_anth_l2_emis )
! --- Send the anthropogenic emission of SO2 to the diag_manager for output.
        if (id_SO2_anth_l2_emis > 0 ) &
          used = send_data ( id_SO2_anth_l2_emis, &
                     SO2_anth_l2_emis*wtm_S/wtm_SO2, Time, &
                              is_in=is, js_in=js )
          
!------------------------------------------------------------------------
! Read 12 monthly 3D aircraft emissions
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading Aircraft emissions, month=',imonth
!
        SO2_aircraft_emis(:,:,:) = 0.  ![kgSO2/box/sec]
        data3D(:,:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo

        open (unit,file='INPUT/'//FNMSO2a//month(imonth)//'.txt', &
              form='formatted',action='read')
        do
          read (unit, '(3i4,e12.4)',iostat=ios) i,j,k,variable
          if ( ios.ne.0 ) exit
          data3D(i,j,k)=variable
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading aircraft emission of SO2
!------------------------------------------------------------------------
! --- Interpolate data
      do k=1,kd 
        call interp_emiss ( data3D(:,:,k), lon_W, lat_S, dlon, dlat, &
                            SO2_aircraft_emis(:,:,k))
      enddo
! Send the emission data to the diag_manager for output.
         if (id_SO2_aircraft_emis > 0 ) &
           used = send_data ( id_SO2_aircraft_emis, &
                So2_aircraft_emis*wtm_S/wtm_SO2, Time, &
                is_in=is, js_in=js , ks_in=1 )
!------------------------------------------------------------------------
! Read 12 monthly 2D biomass burning emission of SO2
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading Aircraft emissions, month=',imonth
        SO2_bioburn_emis(:,:) = 0.  ![kgSO2/box/sec]
        data2D(:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo
        open (unit,file='INPUT/'//FNMBB//month(imonth)//'.txt', &
              form='formatted', action='read')
        do
          read (unit, '(2i4,e12.4)',iostat=ios) i,j,variable
          data2D(i,j)=variable
          if ( ios.ne.0 ) exit
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading biomass burning emission of SO2
!------------------------------------------------------------------------
! --- Interpolate data
        call interp_emiss ( data2D, lon_W, lat_S, dlon, dlat, SO2_bioburn_emis)
! Send the emission data to the diag_manager for output.
         if (id_SO2_bioburn_emis > 0 ) &
           used = send_data ( id_SO2_bioburn_emis, &
                SO2_bioburn_emis*wtm_S/wtm_SO2, Time, is_in=is, js_in=js )
!------------------------------------------------------------------------
! Read 12 monthly NO3 concentration
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading NO3 concentration, month=',imonth
!
        NO3_conc(:,:,:) = 0.  ![molec/cm3]
        data3D(:,:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo

        open (unit,file='INPUT/'//FNMNO3//month(imonth)//'.txt', &
              form='formatted',action='read')
        do
          read (unit, '(3i4,e12.4)',iostat=ios) i,j,k,variable
          if ( ios.ne.0 ) exit
          data3D(i,j,k)=variable
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading NO3 concentration
!------------------------------------------------------------------------
! --- Interpolate data
      do k=1,kd
        call interp_emiss ( data3D(:,:,k), lon_W, lat_S, dlon, dlat, &
                            NO3_conc(:,:,k))
      enddo
! Send the NO3 data to the diag_manager for output.
         if (id_NO3_conc > 0 ) &
           used = send_data ( id_NO3_conc, &
                NO3_conc, Time, is_in=is, js_in=js , ks_in=1 )
!------------------------------------------------------------------------
! Read 12 monthly OH concentration
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading OH concentration, month=',imonth
!
        OH_conc(:,:,:) = 0.  ![molec/cm3]
        data3D(:,:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo

        open (unit,file='INPUT/'//FNMOH//month(imonth)//'.txt', &
              form='formatted',action='read')
        do
          read (unit, '(3i4,e12.4)',iostat=ios) i,j,k,variable
          data3D(i,j,k)=variable
          if ( ios.ne.0 ) exit
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading OH concentration
!------------------------------------------------------------------------
! --- Interpolate data
      do k=1,kd
        call interp_emiss ( data3D(:,:,k), lon_W, lat_S, dlon, dlat, &
                            OH_conc(:,:,k))
      enddo
! Send the OH data to the diag_manager for output.
         if (id_OH_conc > 0 ) &
           used = send_data ( id_OH_conc, &
                OH_conc, Time, is_in=is, js_in=js , ks_in=1 )
!------------------------------------------------------------------------
! Read 12 monthly HO2 concentration
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading HO2 concentration, month=',imonth
!
        HO2_conc(:,:,:) = 0.  ![molec/cm3]
        data3D(:,:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo

        open (unit,file='INPUT/'//FNMHO2//month(imonth)//'.txt', &
               form='formatted',action='read')
        do
          read (unit, '(3i4,e12.4)',iostat=ios) i,j,k,variable
          if ( ios.ne.0 ) exit
          data3D(i,j,k)=variable
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading HO2 concentration
!------------------------------------------------------------------------
! --- Interpolate data
        do k=1,kd
          call interp_emiss ( data3D(:,:,k), lon_W, lat_S, dlon, dlat, &
                            HO2_conc(:,:,k))
        enddo
! Send the HO2 data to the diag_manager for output.
        if (id_HO2_conc > 0 ) &
           used = send_data ( id_HO2_conc, &
                HO2_conc, Time, is_in=is, js_in=js , ks_in=1 )
!------------------------------------------------------------------------
! Read 12 monthly H2O2 photodissociation rates
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading jH2O2, month=',imonth
!
        jH2O2(:,:,:) = 0.  ![molec/cm3]
        data3D(:,:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo

        open (unit,file='INPUT/'//FNMJH2O2//month(imonth)//'.txt', &
               form='formatted',action='read')
        do
          read (unit, '(3i4,e12.4)',iostat=ios) i,j,k,variable
          if ( ios.ne.0 ) exit
          data3D(i,j,k)=variable
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading jH2O2
!------------------------------------------------------------------------
! --- Interpolate data
        do k=1,kd
          call interp_emiss ( data3D(:,:,k), lon_W, lat_S, dlon, dlat, &
                            jH2O2(:,:,k))
        enddo
! Send the jH2O2 data to the diag_manager for output.
        if (id_jH2O2 > 0 ) &
           used = send_data ( id_jH2O2, &
                jH2O2, Time, is_in=is, js_in=js , ks_in=1 )
!------------------------------------------------------------------------
! Read 12 monthly O3 volume mixing ratio
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading O3, month=',imonth
!
        O3_vmr(:,:,:) = 0.  
        data3D(:,:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo

        open (unit,file='INPUT/'//FNMO3//month(imonth)//'.txt', &
               form='formatted',action='read')
        do
          read (unit, '(3i4,e12.4)',iostat=ios) i,j,k,variable
          if ( ios.ne.0 ) exit
          data3D(i,j,k)=variable
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading O3
!------------------------------------------------------------------------
! --- Interpolate data
        do k=1,kd
          call interp_emiss ( data3D(:,:,k), lon_W, lat_S, dlon, dlat, &
                            O3_vmr(:,:,k))
        enddo
! Send the O3 data to the diag_manager for output.
        if (id_O3_vmr > 0 ) &
           used = send_data ( id_O3_vmr, &
                O3_vmr, Time, is_in=is, js_in=js , ks_in=1 )
!------------------------------------------------------------------------
! Read 12 monthly O3
!------------------------------------------------------------------------
if (mpp_pe()== mpp_root_pe() ) &
  write(stdlog(),*) 'Reading pH, month=',imonth
!
        pH(:,:,:) = 0.  
        data3D(:,:,:)=0.
        do unit = 30,100
          INQUIRE(unit=unit, opened= opened)
          if (.NOT. opened) exit
        enddo

        open (unit,file='INPUT/'//FNMPH//month(imonth)//'.txt', &
               form='formatted',action='read')
        do
          read (unit, '(3i4,e12.4)',iostat=ios) i,j,k,variable
          if ( ios.ne.0 ) exit
          data3D(i,j,k)=variable
        enddo
        call close_file (unit)
!------------------------------------------------------------------------
! End reading pH
!------------------------------------------------------------------------
! --- Interpolate data
        do k=1,kd
          call interp_emiss ( data3D(:,:,k), lon_W, lat_S, dlon, dlat, &
                            pH(:,:,k))
        enddo
! Send the pH data to the diag_manager for output.
        if (id_pH > 0 ) &
           used = send_data ( id_pH, &
                pH, Time, is_in=is, js_in=js , ks_in=1 )


! select sources/emissions

   if(trim(runtype) == 'fossil_fuels_only') then
      DMSo             = 0.0        ![nM/L]
      SO2_bioburn_emis =0.0
   else if(trim(runtype) == 'biomass_only') then
      DMSo = 0.0        ![nM/L]
      SO4_anth_l1_emis =0.0
      SO4_anth_l2_emis =0.0
      SO2_anth_l1_emis =0.0
      SO2_anth_l2_emis =0.0
      SO2_aircraft_emis =0.0
   else if(trim(runtype) == 'natural_only')then
      SO4_anth_l1_emis =0.0
      SO4_anth_l2_emis =0.0
      SO2_anth_l1_emis =0.0
      SO2_anth_l2_emis =0.0
      SO2_aircraft_emis =0.0
      SO2_bioburn_emis =0.0
   else if(trim(runtype) == 'anthrop')then
   else
        call error_mesg ('atmos_sulfate_mod', &
             'runtype is not recognized', FATAL)
   endif
!
end subroutine SOx_source_input
!#######################################################################
!</SUBROUTINE>
!<SUBROUTINE NAME="atmos_DMS_emission">
!<OVERVIEW>
! The constructor routine for the sulfate module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to calculate dimethyl sulfide emission form the ocean
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_DMS_emission (r, mask, axes, Time)
!</TEMPLATE>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace). 
!   </INOUT>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>
subroutine atmos_DMS_emission (lon, lat, area, frac_land, t_surf_rad, w10m, &
       pwt, dms_dt, Time, is,ie,js,je,kbot)
!
      real, intent(in),    dimension(:,:)           :: lon, lat
      real, intent(in),    dimension(:,:)           :: frac_land
      real, intent(in),    dimension(:,:)           :: t_surf_rad
      real, intent(in),    dimension(:,:)           :: w10m
      real, intent(in),    dimension(:,:)           :: area
      real, intent(in),    dimension(:,:,:)         :: pwt
      real, intent(out),   dimension(:,:,:)         :: DMS_dt
      type(time_type), intent(in)                   :: Time    
      integer, intent(in)                           :: is, ie, js, je
      integer, intent(in), dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
      real, dimension(size(DMS_dt,1),size(DMS_dt,2)) :: DMS_emis
      integer                                        :: i, j, id, jd, kd
      real                                           :: sst, Sc, conc, w10 
      real                                           :: ScCO2, Akw
      real, parameter                                :: Sc_min=1.

      id=size(dms_dt,1); jd=size(dms_dt,2); kd=size(dms_dt,3)

      dms_dt(:,:,:) =0.0

! ****************************************************************************
! *  If frac_land < 0.5: DMS_emis = seawaterDMS * transfer velocity.         *
! *  Otherwise,  DMS_emis = 0.                                               *
! ****************************************************************************
!

      do j = 1, jd
      do i = 1, id
       SST = t_surf_rad(i,j)-273.15     ! Sea surface temperature [Celsius]
       if (frac_land(i,j).le.0.5) then

!  < Schmidt number for DMS (Saltzman et al., 1993) >
        Sc = 2674.0 - 147.12*SST + 3.726*(SST**2) - 0.038*(SST**3)
        Sc = max(Sc_min, Sc)

! ****************************************************************************
! *  Calculate transfer velocity in cm/hr  (AKw)                             *
! *                                                                          *
! *  Tans et al. transfer velocity (1990) for CO2 at 25oC (Erickson, 1993)   *
! *                                                                          *
! *  Tans et al. assumed AKW=0 when W10<=3. I modified it to let             *
! *  DMS emit at low windseeds too. Chose 3.6m/s as the threshold.           *
! *                                                                          *
! *  Schmidt number for CO2:       Sc = 600  (20oC, fresh water)             *
! *                                Sc = 660  (20oC, seawater)                *
! *                                Sc = 428  (25oC, Erickson 93)             *
! ****************************************************************************
!
        CONC = DMSo(i+is-1,j+js-1)

        W10  = W10M(i,j)

! ---  Tans et al. (1990) -----------------
!       ScCO2 = 428.
!       if (W10 .le. 3.6) then
!        AKw = 1.0667 * W10
!       else
!        AKw = 6.4 * (W10 - 3.)
!       end if

! ---  Wanninkhof (1992) ------------------
!       ScCO2 = 660.
!       AKw = 0.31 * W10**2

! ---  Liss and Merlivat (1986) -----------
        ScCO2 = 600.
        if (W10 .le. 3.6) then
         AKw = 0.17 * W10
        else if (W10 .le. 13.) then
         AKw = 2.85 * W10 - 9.65
        else
         AKw = 5.90 * W10 - 49.3
        end if
!------------------------------------------

        if (W10 .le. 3.6) then
         AKw = AKw * ((ScCO2/Sc) ** 0.667)
        else
         AKw = AKw * sqrt(ScCO2/Sc)
        end if

! ****************************************************************************
! *  Calculate emission flux in kg/m2/s                                  *
! *                                                                          *
! *   AKw is in cm/hr:             AKw/100/3600    -> m/sec.                 *
! *   CONC is in nM/L (nM/dm3):    CONC*1E-12*1000 -> kmole/m3.              *
! *   WTM_DMS          : kgDMS/kmol.                                         *
! *   DMS_EMIS         : kgDMS/m2/s.                                         *
! ****************************************************************************
!
        DMS_emis(i,j) = AKw/100./3600. * CONC*1.e-12*1000.* WTM_DMS &
            * (1.-frac_land(i,j))
!
       else                !  frac_land <> 1 (water)
        DMS_emis(i,j) = 0.

       end if              ! -- if frac_land = 1.

      end do
      end do
!--------------------------------------------------------------------
! Update DMS concentration in level kd (where emission occurs)
!--------------------------------------------------------------------
      dms_dt(:,:,kd)=DMS_emis(:,:)/pwt(:,:,kd)* WTMAIR/WTM_DMS
!------------------------------------------------------------------
! DIAGNOSTICS:      DMS surface emission in kg/m2/s     
!--------------------------------------------------------------------
      if (id_DMS_emis > 0) then
        used = send_data ( id_DMS_emis, dms_emis*WTM_S/WTM_DMS, Time, &
              is_in=is,js_in=js )
      endif

end subroutine atmos_DMS_emission
!#######################################################################
!</SUBROUTINE>
!<SUBROUTINE NAME="atmos_SO2_emission">
!<OVERVIEW>
! The constructor routine for the sulfate module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to calculate SO2 emission from volcanoes, biomass burning,
! anthropogenic sources, aircraft.
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_SO2_emission ()
!</TEMPLATE>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace).
!   </INOUT>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>
subroutine atmos_SO2_emission (lon, lat, area, frac_land, &
       z_pbl, zhalf, pwt, so2_nerup_volc_emis, SO2_dt, Time, is,ie,js,je,kbot)
!
      real, intent(in),    dimension(:,:)           :: lon, lat
      real, intent(in),    dimension(:,:)           :: frac_land
      real, intent(in),    dimension(:,:)           :: area
      real, intent(in),    dimension(:,:)           :: z_pbl
      real, intent(in),    dimension(:,:,:)         :: so2_nerup_volc_emis
      real, intent(in),    dimension(:,:,:)         :: zhalf
      real, intent(in),    dimension(:,:,:)         :: pwt
      real, intent(out),   dimension(:,:,:)         :: SO2_dt
      type(time_type), intent(in)                   :: Time
      integer, intent(in)                           :: is, ie, js, je
      integer, intent(in), dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      real, dimension(size(SO2_dt,1),size(SO2_dt,2),size(SO2_dt,3)) :: SO2_emis
      integer  :: i, j, l, id, jd, kd

      real :: fbb(size(SO2_dt,3)),fa1(size(SO2_dt,3)),fa2(size(SO2_dt,3))

      real, parameter :: ze1=100.
      real, parameter :: ze2=500.
      real :: z1, z2, bltop, fbt

      id=size(SO2_dt,1); jd=size(SO2_dt,2); kd=size(SO2_dt,3)

      SO2_dt(:,:,:) = 0.0

      do j = 1, jd
      do i = 1, id

        fbb(:) = 0.
        fa1(:) = 0.
        fa2(:) = 0.

! --- Assuming biomass burning emission within the PBL -------
        fbt=0.
        BLTOP = z_pbl(i,j)
        do l = kd,2,-1
          z1=zhalf(i,j,l+1)-zhalf(i,j,kd+1)
          z2=zhalf(i,j,l)-zhalf(i,j,kd+1)
          if (bltop.lt.z1) exit
          if (bltop.ge.z2) fbb(l)=pwt(i,j,l)
          if (bltop.gt.z1.and.bltop.lt.z2) then
            fbb(l) = pwt(i,j,l)*(bltop-z1)/(z2-z1)
          endif
          fbt=fbt+fbb(l)
        enddo
        if (fbt .gt. 0.) fbb(:)=fbb(:)/fbt
 

! --- Assuming anthropogenic source L1 emitted below Ze1, and L2
!     emitted between Ze1 and Ze2.
        do l = kd,2,-1
          Z1 = zhalf(i,j,l+1)-zhalf(i,j,kd+1)
          Z2 = zhalf(i,j,l)-zhalf(i,j,kd+1)
          if (Z2.lt.Ze1) then
            fa1(l) = (Z2-Z1)/Ze1
            fa2(l) = 0.
          endif
          if (Z2.ge.Ze1 .and. Z1.lt.Ze1) then
            fa1(l) = (Ze1-Z1)/Ze1
            if ((Ze2-Ze1).gt.0.) fa2(l) = (Z2-Ze1)/(Ze2-Ze1)
          endif
          if (Z1.gt.Ze1) then
            if (Z2.lt.Ze2.and.(Ze2-Ze1).gt.0.) &
              fa2(l) = (z2-z1)/(Ze2-Ze1)
            if (Z2.ge.Ze2 .and. Z1.lt.Ze2) then
              if ((Ze2-Ze1).gt.0.) fa2(l) = (Ze2-Z1)/(Ze2-Ze1)
            endif
          endif
          if (Z1.gt.Ze2) exit
        enddo
! --- Total SO2 source ----
! The emission rates are in units of kgSO2/box/sec
! SO2_emis: [kgSO2/m2/s]
        SO2_emis(i,j,:) =( SO2_aircraft_emis(i+is-1,j+js-1,:)           &
                        +  SO2_nerup_volc_emis(i,j,:)                   &
                        +  SO2_anth_l1_emis(i+is-1,j+js-1)*fa1(:)       &
                        +  SO2_anth_l2_emis(i+is-1,j+js-1)*fa2(:)       &
                        +  SO2_bioburn_emis(i+is-1,j+js-1)*fbb(:) )     &
                        /  area(i,j)
      end do
      end do
!
      SO2_dt(:,:,:)=SO2_emis(:,:,:)/pwt(:,:,:)*WTMAIR/WTM_SO2
!------------------------------------------------------------------
! DIAGNOSTICS:      SO2 and SO4 emission in kg/timestep
!--------------------------------------------------------------------
      if (id_SO2_emis > 0) then
        used = send_data ( id_SO2_emis, SO2_emis*WTM_S/WTM_SO2, Time, &
              is_in=is,js_in=js,ks_in=1)
      endif

end subroutine atmos_SO2_emission
!</SUBROUTINE>
!-----------------------------------------------------------------------
!#######################################################################
!<SUBROUTINE NAME="atmos_SO4_emission">
!<OVERVIEW>
! The constructor routine for the sulfate module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to calculate SO4 emission from volcanoes, biomass burning,
! anthropogenic sources, aircraft.
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_SO4_emission ()
!</TEMPLATE>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace).
!   </INOUT>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>
subroutine atmos_SO4_emission (lon, lat, area, frac_land, &
       z_pbl, zhalf, pwt, SO4_dt, &
       Time, is,ie,js,je,kbot)
!
      real, intent(in),    dimension(:,:)           :: lon, lat
      real, intent(in),    dimension(:,:)           :: frac_land
      real, intent(in),    dimension(:,:)           :: area
      real, intent(in),    dimension(:,:)           :: z_pbl
      real, intent(in),    dimension(:,:,:)         :: zhalf
      real, intent(in),    dimension(:,:,:)         :: pwt
      real, intent(out),   dimension(:,:,:)         :: SO4_dt
      type(time_type), intent(in)                   :: Time
      integer, intent(in)                           :: is, ie, js, je
      integer, intent(in), dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      real, dimension(size(SO4_dt,1),size(SO4_dt,2),size(SO4_dt,3)) :: SO4_emis
      integer  :: i, j, l, id, jd, kd

      real :: fbb(size(SO4_dt,3)),fa1(size(SO4_dt,3)),fa2(size(SO4_dt,3))

      real, parameter :: ze1=100.
      real, parameter :: ze2=500.
      real :: z1, z2, bltop

      id=size(SO4_dt,1); jd=size(SO4_dt,2); kd=size(SO4_dt,3)

      SO4_dt(:,:,:) =0.0

      do j = 1, jd
      do i = 1, id

        fbb(:) = 0.
        fa1(:) = 0.
        fa2(:) = 0.

        BLTOP = z_pbl(i,j)
        do l = kd,1,-1
          z1=zhalf(i,j,l+1)-zhalf(i,j,kd+1)
          z2=zhalf(i,j,l)-zhalf(i,j,kd+1)
          if (bltop.lt.z1) exit
          if (bltop.ge.z2) fbb(l)=(z2-z1)/z_pbl(i,j)
          if (bltop.gt.z1.and.bltop.lt.z2) then
            fbb(l) = (bltop-z1)/z_pbl(i,j)
          endif
        enddo
! --- Assuming anthropogenic source L1 emitted below Ze1, and L2
!     emitted between Ze1 and Ze2.
        do l = kd,2,-1
          Z1 = zhalf(i,j,l+1)-zhalf(i,j,kd+1)
          Z2 = zhalf(i,j,l)-zhalf(i,j,kd+1)
          if (Z2.lt.Ze1) then
            fa1(l) = (Z2-Z1)/Ze1
            fa2(l) = 0.
          endif
          if (Z2.ge.Ze1 .and. Z1.lt.Ze1) then
            fa1(l) = (Ze1-Z1)/Ze1
            if ((Ze2-Ze1).gt.0.) fa2(l) = (Z2-Ze1)/(Ze2-Ze1)
          endif
          if (Z1.gt.Ze1) then
            if (Z2.lt.Ze2.and.(Ze2-Ze1).gt.0.) &
              fa2(l) = (z2-z1)/(Ze2-Ze1)
            if (Z2.ge.Ze2 .and. Z1.lt.Ze2) then
              if ((Ze2-Ze1).gt.0.) fa2(l) = (Ze2-Z1)/(Ze2-Ze1)
            endif
          endif
          if (Z1.gt.Ze2) exit
        enddo

! --- Total SO4 source ----
!    Anthropogenic SOx emission from GEIA 1985.
!    Assuming:   Europe:      5.0% SOx emission is SO4;
!                US + Canada: 1.4% SOx emission is SO4;
!                The rest:    2.5% SOx emission is SO4.
!
!
        SO4_emis(i,j,:) =( SO4_anth_l1_emis(i+is-1,j+js-1)*fa1(:)       &
                        +  SO4_anth_l2_emis(i+is-1,j+js-1)*fa2(:) )     &
                        /  area(i,j) 
      end do
      end do
!

      SO4_dt(:,:,:)=SO4_emis(:,:,:)/pwt(:,:,:)*WTMAIR/WTM_SO4
!------------------------------------------------------------------
! DIAGNOSTICS:      SO4 emission in kg/timestep
!--------------------------------------------------------------------
      if (id_SO4_emis > 0) then
        used = send_data ( id_SO4_emis, SO4_emis*WTM_S/WTM_SO4, Time, &
              is_in=is,js_in=js,ks_in=1)
      endif

end subroutine atmos_SO4_emission
!</SUBROUTINE>
!-----------------------------------------------------------------------
!#######################################################################
      subroutine chem_sox(pwt,temp,pfull, dt, lwc, &
        jday,hour,minute,second,lat,lon, &
        SO2, SO4, DMS, MSA, H2O2, &
        SO2_dt, SO4_dt, DMS_dt, MSA_dt, H2O2_dt, &
        Time,is,ie,js,je,kbot)
!
      real, intent(in)                   :: dt
      integer, intent(in)                :: jday, hour,minute,second
      real, intent(in),  dimension(:,:)  :: lat, lon  ! [radi
      real, intent(in), dimension(:,:,:) :: pwt
      real, intent(in), dimension(:,:,:) :: lwc
      real, intent(in), dimension(:,:,:) :: temp, pfull
      real, intent(in), dimension(:,:,:) :: SO2, SO4, DMS, MSA, H2O2
      real, intent(out),dimension(:,:,:) :: SO2_dt,SO4_dt,DMS_dt,MSA_dt,H2O2_dt

      type(time_type), intent(in)                    :: Time
      integer, intent(in),  dimension(:,:), optional :: kbot
      integer, intent(in)                            :: is,ie,js,je
! Working vectors
      integer :: i,j,k,id,jd,kd, istop
      integer                                    :: istep, nstep
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: ppH
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: pO3
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: no3_diurnal
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: oh_diurnal
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) ::jh2o2_diurnal
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: ho2_diurnal
      real, dimension(size(pfull,1),size(pfull,2)) :: &
               xu, dayl, h, hl, hc, hred
      real, dimension(size(pfull,1),size(pfull,2)) :: fac_NO3, fact_NO3
      real, dimension(size(pfull,1),size(pfull,2)) :: fac_OH , fact_OH 
      real, dimension(size(pfull,1),size(pfull,2)) :: fac_HO2            
      real, parameter                            :: A0 = 0.006918
      real, parameter                            :: A1 = 0.399912
      real, parameter                            :: A2 = 0.006758
      real, parameter                            :: A3 = 0.002697
      real, parameter                            :: B1 = 0.070257
      real, parameter                            :: B2 = 0.000907
      real, parameter                            :: B3 = 0.000148
      real                                       :: decl, hd, x
      real :: f, f1, tk, rho_air
      real :: SO2_0,SO4_0,MSA_0,DMS_0,H2O2_0    ! initial concentrations
      real :: xSO2,xSO4,xMSA,xDMS,xH2O2,xno3,xo3,xoh,xho2,xjh2o2 ! update conc.
      real :: rk0, rk1, rk2, rk3, rk, rkt  ! kinetic rates
      real :: cfact, work1, xk, xe, x2, xph
      real :: heh2o2, h2o2g, rah2o2, px, heso2, so2g, heo3, o3g, rao3
      real :: pso4a, pso4b
      real :: xlwc, xhnm, ccc1, ccc2
      real :: pmsa, pso2, pso4, ph2o2    ! chemical production terms
      real :: ldms, lso2, lh2o2          ! chemical loss terms
      real :: o2
      real, parameter        :: small_value=1.e-21
      real, parameter        :: t0 = 298.
      real, parameter        :: Ra = 8314./101325.
      real, parameter        :: xkw = 1.e-14 ! water acidity
      real, parameter        :: const0 = 1.e3/6.022e23


! Local grid sizes
      id=size(pfull,1) ; jd=size(pfull,2) ; kd=size(pfull,3)

      so2_dt(:,:,:) = 0.0
      so4_dt(:,:,:) = 0.0
      dms_dt(:,:,:) = 0.0
      msa_dt(:,:,:) = 0.0
      h2o2_dt(:,:,:) = 0.0

      x = 2. *pi *float(jday-1)/365.
      decl = A0 - A1*cos(  X) + B1*sin(  X) - A2*cos(2.*X) + B2*sin(2.*X) &
           - A3*cos(3.*X) + B3*sin(3.*X)
      xu(:,:) = -tan(lat(:,:))*tan(decl)

      where ( xu > -1 .and. xu < 1 ) dayl=acos(xu)/pi
      where ( xu <= -1 ) dayl = 1.
      where ( xu >= 1 ) dayl = 0.
!   Calculate normalization factors for OH and NO3 such that
!   the diurnal average respect the monthly input values.
      hd=0.
      fact_OH(:,:)  = 0.
      fact_NO3(:,:) = 0.
      nstep = int(24.*3600./dt)
      do istep=1,nstep
        hd=hd+dt/3600./24.
        hl(:,:) = pi*(1.-dayl(:,:))
        hc(:,:) = pi*(1.+dayl(:,:))
        h(:,:)=2.*pi*mod(hd+lon(:,:)/2./pi,1.)
        where ( h.ge.hl .and. h.lt.hc )
! Daytime
          hred=(h-hl)/(hc-hl)
          fact_OH  = fact_OH + amax1(0.,sin(pi*hred)/2.)/nstep
        elsewhere
! Nightime
          fact_NO3 = fact_NO3 + amax1(0.,amin1(1.,(1.-dayl)))/nstep
        endwhere
      enddo

      hd=amax1(0.,amin1(1.,(hour+minute/60.+second/3600.)/24.))
      hl(:,:) = pi*(1.-dayl(:,:))
      hc(:,:) = pi*(1.+dayl(:,:))
      h(:,:)=2.*pi*mod(hd+lon(:,:)/2./pi,1.)
      fac_OH(:,:)  = 0.
      fac_NO3(:,:) = 0.
      fac_HO2(:,:) = 1.
      where ( h.ge.hl .and. h.lt.hc )
! Daytime
          fac_NO3 = 0.
          hred=(h-hl)/(hc-hl)
          fac_OH  = amax1(0.,sin(pi*hred)/2.)/fact_OH
      elsewhere
! Nightime
          fac_NO3 = amax1(0.,amin1(1.,(1.-dayl)))/fact_NO3
          fac_OH  = 0.
      endwhere


! < Factor to convert AIRDEN from kgair/m3 to molecules/cm3: >
      f  = 1000. / WTMAIR * 6.022e23 * 1.e-6

      do k = 1, kd
      do j = 1, jd
      do i = 1, id
       tk    = temp(i,j,k)
       rho_air = pfull(i,j,k)/tk/RDGAS             ! Air density [kg/m3]
       xhnm  = rho_air * f
       O2    = xhnm * 0.21
       xlwc  = lwc(i,j,k)*rho_air *1.e-3
       DMS_0 = max(0.,DMS(i,j,k))
       MSA_0 = max(0.,MSA(i,j,k))
       SO4_0 = max(0.,SO4(i,j,k))
       SO2_0 = max(0.,SO2(i,j,k))
       H2O2_0= max(0.,H2O2(i,j,k))
       xSO2  = SO2_0
       xSO4  = SO4_0
       xH2O2 = H2O2_0
       xDMS  = DMS_0
       xMSA  = MSA_0
       xph   = max(1.e-7,       pH(i+is-1,j+js-1,k))
       xoh   = max(0.         , OH_conc(i+is-1,j+js-1,k)  *fac_OH(i,j))
       xho2  = max(0.         , HO2_conc(i+is-1,j+js-1,k) *fac_HO2(i,j))
       xjh2o2= max(0.         , jH2O2(i+is-1,j+js-1,k)    *fac_OH(i,j))
       xno3  = max(0.         , NO3_conc(i+is-1,j+js-1,k) *fac_NO3(i,j))
       xo3   = max(small_value, O3_vmr(i+is-1,j+js-1,k))
       pph(i,j,k)=xph
       po3(i,j,k)=xo3
       oh_diurnal(i,j,k)=xoh
       no3_diurnal(i,j,k)=xno3
       ho2_diurnal(i,j,k)=xho2
       jh2o2_diurnal(i,j,k)=xjh2o2
! ****************************************************************************
! *  H2O2 production by HO2 + HO2 reactions
! ****************************************************************************
       PH2O2=(2.2e-13*exp(619./tk)+xhnm*1.9e-33*exp(980./tk))* xHO2**2 /xhnm
! ****************************************************************************
! *  H2O2 loss by OH and photodissociation
! ****************************************************************************
       LH2O2= ( 2.9e-12*exp(-160./tk)* xOH + xjH2O2 )
       if (LH2O2 .gt. 0.) then
         xH2O2= H2O2_0 * exp(-LH2O2*dt) + PH2O2*(1.-exp(-LH2O2*dt))/LH2O2
       else
         xH2O2= H2O2_0 + PH2O2 * dt
       endif
! ****************************************************************************
! *  (1) DMS + OH:  RK1 - addition channel;  RK2 - abstraction channel.      *
! ****************************************************************************
       rk1 = (1.7e-42 * exp(7810./TK) * O2) /   &
              (1. + 5.5e-31 * exp(7460./TK) * O2 ) * xoh
       rk2 = 1.2e-11*exp(-260./TK) * xoh
! ****************************************************************************
! *  (2) DMS + NO3 (only happen at night):                                   *
! ****************************************************************************
!  < XNO3 fields are in molecules/cm3.        >
        rk3 = 1.9e-13 * exp(500./TK) * xno3
! ****************************************************************************
! *  Update DMS concentration after gas phase chemistry                      *
! ****************************************************************************
       LDMS = RK1 + RK2 + RK3
       if ( LDMS .gt. 0. ) then
         xDMS = DMS_0 * exp( - LDMS*dt)
       endif
! ****************************************************************************
! *  Update MSA concentration after gas phase chemistry                      *
! ****************************************************************************
       PMSA = RK1*0.25 * xDMS
       xMSA = MSA_0 + PMSA * dt
! ****************************************************************************
! *  SO2 oxydation by OH
! ****************************************************************************
       PSO2 = ( RK1*0.75 + RK2 + RK3 ) * xDMS
       rk0 = 3.0E-31 * (300./TK)**3.3
       rk1 = rk0 * xhnm / 1.5e-12
       f1 = ( 1.+ ( log10(rk1) )**2 )**(-1)
       LSO2 = ( rk0 * xhnm / (1.+ rk1) ) * 0.6**f1 * xoh
! ****************************************************************************
! *  Update SO2 concentration after gas phase chemistry                      *
! ****************************************************************************
       xSO2 = SO2_0 + dt * ( PSO2 - LSO2 * SO2_0)
       if (xSO2 .lt. 0.) then
        xSO2 = SO2_0 * exp(-LSO2*dt) + PSO2 * (1.-exp(-LSO2*dt)) / LSO2
       end if
! ****************************************************************************
! *  Update SO4 concentration after gas phase chemistry                      *
! ****************************************************************************
       xso4 = SO4_0 + LSO2*xso2 * dt
! ****************************************************************************
! < Cloud chemistry (above 258K): >
       work1 = (t0 - tk)/(tk*t0)
!-----------------------------------------------------------------------
!         ... h2o2
!-----------------------------------------------------------------------
       xk = 7.4e4   *exp( 6621.* work1 )
       xe = 2.2e-12 *exp(-3730.* work1 )
       heh2o2  = xk*(1. + xe/xph)
       px = heh2o2 * Ra * tk * xlwc
       h2o2g = xh2o2 /(1.+px) 
!-----------------------------------------------------------------------
!         ... so2
!-----------------------------------------------------------------------
       xk = 1.23   * exp(3120. * work1 )
       xe = 1.7e-2 * exp(2090. * work1 ) 
       x2 = 6.0e-8 * exp(1120. * work1 )
       heso2 = xk*(1. + xe/xph *(1. + x2/xph) ) 
!       heso2 = 1.e2 ! xk*(1. + xe/xph *(1. + x2/xph) )
       px = heso2 * Ra * tk * xlwc
       so2g = xso2/(1.+px)
!-----------------------------------------------------------------------
!         ... o3
!-----------------------------------------------------------------------
       xk = 1.15e-2 * exp( 2560. * work1 )
       heo3 = xk
       px = heo3 * Ra * tk *xlwc
       o3g = xo3 / (1.+px) 
!-----------------------------------------------
!       ... Aqueous phase reaction rates
!           SO2 + H2O2 -> SO4
!           SO2 + O3   -> SO4
!-----------------------------------------------

!------------------------------------------------------------------------
!       ... S(IV) (HSO3) + H2O2
!------------------------------------------------------------------------
            rah2o2 = 8.e4 * EXP( -3650.*work1 )  / (.1 + xph)

!------------------------------------------------------------------------
!        ... S(IV)+ O3
!------------------------------------------------------------------------
            rao3   = 4.39e11 * EXP(-4131./tk)  &
                  + 2.56e3  * EXP(-996. /tk) /xph

!-----------------------------------------------------------------
!       ... Prediction after aqueous phase
!       so4
!       When Cloud is present
!
!       S(IV) + H2O2 = S(VI)
!       S(IV) + O3   = S(VI)
!
!       reference:
!           (1) Seinfeld
!           (2) Benkovitz
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!       ... S(IV) + H2O2 = S(VI)
!-----------------------------------------------------------------
       ccc1=0.
       ccc2=0.
       if( xlwc >= 1.e-8 ) then                    ! when cloud is present
               pso4a = rah2o2 * heh2o2*h2o2g  &
                             * heso2 *so2g            ! [M/s]
               pso4a = pso4a       &                    ! [M/s] =  [mole/L(w)/s]
                    * xlwc       &                    ! [mole/L(a)/s]
                    / const0     &                    ! [/L(a)/s]
                    / xhnm                            ! [mixing ratio/s]

          ccc1 = pso4a*dt
          ccc1 = max(min(ccc1,xso2,xh2o2), 0.)
          xso4 = xso4 + ccc1
          xh2o2 = max(xh2o2 - ccc1, small_value)
          xso2 =  max(xso2  - ccc1, small_value)
!          ccc1 = max(ccc1, 0.)
!          if( xh2o2 > xso2 ) then
!              if( ccc1 > xso2 ) then
!                  xso4  = xso4 + xso2
!                  xso2  = small_value
!                  xh2o2 = xh2o2 - xso2
!              else
!                  xso4  = xso4  + ccc1
!                  xh2o2 = xh2o2 - ccc1
!                  xso2  = xso2  - ccc1
!              end if
!          else
!               if( ccc1 > xh2o2 ) then
!                   xso4  = xso4 + xh2o2
!                   xso2  = xso2 - xh2o2
!                   xh2o2 = small_value
!               else
!                   xso4  = xso4  + ccc1
!                   xh2o2 = xh2o2 - ccc1
!                   xso2  = xso2  - ccc1
!               end if
!          end if


!-----------------------------------------------
!       ... S(IV) + O3 = S(VI)
!-----------------------------------------------
           pso4b = rao3 * heo3*o3g * heso2*so2g       ! [M/s]
           pso4b = pso4b        &        ! [M/s] =  [mole/L(w)/s]
                * xlwc        &        ! [mole/L(a)/s]
                / const0      &        ! [/L(a)/s]
                / xhnm                 ! [mixing ratio/s]
 
           ccc2 = pso4b*dt
            ccc2 = max(min(ccc2, xso2), 0.)               ! mozart2
            xso4 = xso4 + ccc2                           ! mozart2
            xso2 = max(xso2 - ccc2, small_value)         ! mozart2
!          ccc2 = max(ccc2, 0.)
!          if( ccc2 > xso2 ) then
!             xso4 = xso4 + xso2
!             xso2 = small_value
!          else
!             xso4 = xso4  + ccc2
!             xso2 = xso2  - ccc2
!             xso2 = max(xso2, small_value)
!          end if
       end if
       MSA_dt(i,j,k) = (xMSA-MSA_0)/dt
       DMS_dt(i,j,k) = (xDMS-DMS_0)/dt
       SO2_dt(i,j,k) = (xso2-SO2_0)/dt
       SO4_dt(i,j,k) = (xso4-SO4_0)/dt
       H2O2_dt(i,j,k)= (xh2o2-H2O2_0)/dt
      end do
      end do
      end do
      if ( id_NO3_diurnal > 0) then
         used = send_data ( id_NO3_diurnal, NO3_diurnal, &
           Time,is_in=is,js_in=js,ks_in=1)
      endif
      if ( id_OH_diurnal > 0) then
        used = send_data ( id_OH_diurnal, OH_diurnal, &
          Time, is_in=is, js_in=js )
      endif
      if ( id_HO2_diurnal > 0) then
        used = send_data ( id_HO2_diurnal, HO2_diurnal, &
          Time, is_in=is, js_in=js )
      endif
      if ( id_jH2O2_diurnal > 0) then
        used = send_data ( id_jH2O2_diurnal, jH2O2_diurnal, &
          Time, is_in=is, js_in=js )
      endif
      if (id_pph > 0) then
        used = send_data ( id_pph, &
              pph, Time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_po3 > 0) then
        used = send_data ( id_po3, &
              po3, Time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_SO2_chem > 0) then
        used = send_data ( id_SO2_chem, &
              SO2_dt*pwt*WTM_S/WTMAIR, Time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_SO4_chem > 0) then
        used = send_data ( id_SO4_chem, &
              SO4_dt*pwt*WTM_S/WTMAIR, Time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_DMS_chem > 0) then
        used = send_data ( id_DMS_chem, &
              DMS_dt*pwt*WTM_S/WTMAIR, Time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_MSA_chem > 0) then
        used = send_data ( id_MSA_chem, &
              MSA_dt*pwt*WTM_S/WTMAIR, Time,is_in=is,js_in=js,ks_in=1)
      endif

      if (id_H2O2_chem > 0) then
        used = send_data ( id_H2O2_chem, &
              H2O2_dt*pwt, Time,is_in=is,js_in=js,ks_in=1)
      endif
end subroutine chem_sox
!=============================================================================
!<SUBROUTINE NAME="get_SO2_nerup_volc_emis">
!<OVERVIEW>
!  This subroutine builds the emission rates of SO2 by non-eruptive volcanoes
!</OVERVIEW>
subroutine get_SO2_nerup_volc_emis(so2_nerup_volc_emis, &
        zhalf,lon,lat,Time,is,ie,js,je)
!-----------------------------------------------------------------------
        real, intent(out),dimension(:,:,:) :: so2_nerup_volc_emis
        real, intent(in), dimension(:,:,:) :: zhalf
        real, intent(in), dimension(:,:)   :: lon, lat
        integer, intent(in)                :: is,ie,js,je
        type(time_type), intent(in)        :: Time

        integer              :: iv, k, ilon,ilat
        integer              :: id, jd, kd
        real                 :: latmin,latmax,lonmin,lonmax
        type volcano
          character (len=26) :: name  ! Volcano name
          real               :: lat   ! Volcano latitude in degrees [-90,90]
          real               :: lon   ! Volcano longitude in degrees [-180,180]
          real               :: alt   ! Volcano altitude in meters
          real               :: emis  ! Volcano emission of MgSO2/day
        end type volcano
        integer, parameter               :: nerup=45
        type (volcano), dimension(nerup) :: nerup_volc

        nerup_volc(:)%name=" "
        nerup_volc(:)%lat=0.
        nerup_volc(:)%lon=0.
        nerup_volc(:)%emis=0.
        
!--------------------------------------------------------------------------
! Input values for the nerup volcanoes continuously emtting SO2
!--------------------------------------------------------------------------
        nerup_volc( 1)%name = "                 STROMBOLI"
        nerup_volc( 1)%lat  =     38.789
        nerup_volc( 1)%lon  =     15.213
        nerup_volc( 1)%alt  =    926.000
        nerup_volc( 1)%emis =   0.730E+03
        nerup_volc( 2)%name = "                   VULCANO"
        nerup_volc( 2)%lat  =     38.404
        nerup_volc( 2)%lon  =     14.962
        nerup_volc( 2)%alt  =    500.000
        nerup_volc( 2)%emis =   0.440E+02
        nerup_volc( 3)%name = "                      ETNA"
        nerup_volc( 3)%lat  =     37.734
        nerup_volc( 3)%lon  =     15.004
        nerup_volc( 3)%alt  =   3350.000
        nerup_volc( 3)%emis =   0.400E+04
        nerup_volc( 4)%name = "                      ERTA"
        nerup_volc( 4)%lat  =     13.600
        nerup_volc( 4)%lon  =     40.670
        nerup_volc( 4)%alt  =    613.000
        nerup_volc( 4)%emis =   0.210E+02
        nerup_volc( 5)%name = "                   LENGAI,"
        nerup_volc( 5)%lat  =     -2.751
        nerup_volc( 5)%lon  =     35.902
        nerup_volc( 5)%alt  =   2890.000
        nerup_volc( 5)%emis =   0.160E+02
        nerup_volc( 6)%name = "                     WHITE"
        nerup_volc( 6)%lat  =    -37.520
        nerup_volc( 6)%lon  =    177.180
        nerup_volc( 6)%alt  =    321.000
        nerup_volc( 6)%emis =   0.520E+03
        nerup_volc( 7)%name = "                     MANAM"
        nerup_volc( 7)%lat  =     -4.100
        nerup_volc( 7)%lon  =    145.061
        nerup_volc( 7)%alt  =   1807.000
        nerup_volc( 7)%emis =   0.920E+03
        nerup_volc( 8)%name = "                   LANGILA"
        nerup_volc( 8)%lat  =     -5.525
        nerup_volc( 8)%lon  =    148.420
        nerup_volc( 8)%alt  =   1330.000
        nerup_volc( 8)%emis =   0.690E+03
        nerup_volc( 9)%name = "                    ULAWUN"
        nerup_volc( 9)%lat  =     -5.050
        nerup_volc( 9)%lon  =    151.330
        nerup_volc( 9)%alt  =   2334.000
        nerup_volc( 9)%emis =   0.480E+03
        nerup_volc(10)%name = "                    BAGANA"
        nerup_volc(10)%lat  =     -6.140
        nerup_volc(10)%lon  =    155.195
        nerup_volc(10)%alt  =   1750.000
        nerup_volc(10)%emis =   0.330E+04
        nerup_volc(11)%name = "                     YASUR"
        nerup_volc(11)%lat  =    -19.520
        nerup_volc(11)%lon  =    169.425
        nerup_volc(11)%alt  =    361.000
        nerup_volc(11)%emis =   0.900E+03
        nerup_volc(12)%name = "           TANGKUBANPARAHU"
        nerup_volc(12)%lat  =     -6.770
        nerup_volc(12)%lon  =    107.600
        nerup_volc(12)%alt  =   2084.000
        nerup_volc(12)%emis =   0.750E+02
        nerup_volc(13)%name = "                    SLAMET"
        nerup_volc(13)%lat  =     -7.242
        nerup_volc(13)%lon  =    109.208
        nerup_volc(13)%alt  =   3432.000
        nerup_volc(13)%emis =   0.580E+02
        nerup_volc(14)%name = "                    MERAPI"
        nerup_volc(14)%lat  =     -7.542
        nerup_volc(14)%lon  =    110.442
        nerup_volc(14)%alt  =   2911.000
        nerup_volc(14)%emis =   0.140E+03
        nerup_volc(15)%name = "                   TENGGER"
        nerup_volc(15)%lat  =     -7.942
        nerup_volc(15)%lon  =    112.950
        nerup_volc(15)%alt  =   2329.000
        nerup_volc(15)%emis =   0.140E+02
        nerup_volc(16)%name = "                   BULUSAN"
        nerup_volc(16)%lat  =     12.770
        nerup_volc(16)%lon  =    124.050
        nerup_volc(16)%alt  =   1565.000
        nerup_volc(16)%emis =   0.370E+03
        nerup_volc(17)%name = "                     MAYON"
        nerup_volc(17)%lat  =     13.257
        nerup_volc(17)%lon  =    123.685
        nerup_volc(17)%alt  =   2462.000
        nerup_volc(17)%emis =   0.530E+03
        nerup_volc(18)%name = "                     KIKAI"
        nerup_volc(18)%lat  =     30.780
        nerup_volc(18)%lon  =    130.280
        nerup_volc(18)%alt  =    717.000
        nerup_volc(18)%emis =   0.570E+03
        nerup_volc(19)%name = "               SAKURA-JIMA"
        nerup_volc(19)%lat  =     31.580
        nerup_volc(19)%lon  =    130.670
        nerup_volc(19)%alt  =   1117.000
        nerup_volc(19)%emis =   0.190E+04
        nerup_volc(20)%name = "                     UNZEN"
        nerup_volc(20)%lat  =     32.750
        nerup_volc(20)%lon  =    130.300
        nerup_volc(20)%alt  =   1359.000
        nerup_volc(20)%emis =   0.130E+03
        nerup_volc(21)%name = "                       ASO"
        nerup_volc(21)%lat  =     32.880
        nerup_volc(21)%lon  =    131.100
        nerup_volc(21)%alt  =   1592.000
        nerup_volc(21)%emis =   0.760E+02
        nerup_volc(22)%name = "                      KUJU"
        nerup_volc(22)%lat  =     33.080
        nerup_volc(22)%lon  =    131.250
        nerup_volc(22)%alt  =   1788.000
        nerup_volc(22)%emis =   0.140E+03
        nerup_volc(23)%name = "                     ASAMA"
        nerup_volc(23)%lat  =     36.400
        nerup_volc(23)%lon  =    138.530
        nerup_volc(23)%alt  =   2560.000
        nerup_volc(23)%emis =   0.370E+03
        nerup_volc(24)%name = "                    OSHIMA"
        nerup_volc(24)%lat  =     34.730
        nerup_volc(24)%lon  =    139.380
        nerup_volc(24)%alt  =    758.000
        nerup_volc(24)%emis =   0.270E+03
        nerup_volc(25)%name = "                       USU"
        nerup_volc(25)%lat  =     42.530
        nerup_volc(25)%lon  =    140.830
        nerup_volc(25)%alt  =    731.000
        nerup_volc(25)%emis =   0.560E+02
        nerup_volc(26)%name = "                 MEDVEZHIA"
        nerup_volc(26)%lat  =     45.380
        nerup_volc(26)%lon  =    148.830
        nerup_volc(26)%alt  =   1124.000
        nerup_volc(26)%emis =   0.680E+02
        nerup_volc(27)%name = "                    MARTIN"
        nerup_volc(27)%lat  =     58.170
        nerup_volc(27)%lon  =   -155.350
        nerup_volc(27)%alt  =   1860.000
        nerup_volc(27)%emis =   0.300E+01
        nerup_volc(28)%name = "                 AUGUSTINE"
        nerup_volc(28)%lat  =     59.370
        nerup_volc(28)%lon  =   -153.420
        nerup_volc(28)%alt  =   1252.000
        nerup_volc(28)%emis =   0.480E+02
        nerup_volc(29)%name = "                   ILIAMNA"
        nerup_volc(29)%lat  =     60.030
        nerup_volc(29)%lon  =   -153.080
        nerup_volc(29)%alt  =   3053.000
        nerup_volc(29)%emis =   0.220E+02
        nerup_volc(30)%name = "                   KILAUEA"
        nerup_volc(30)%lat  =     19.425
        nerup_volc(30)%lon  =   -155.292
        nerup_volc(30)%alt  =   1222.000
        nerup_volc(30)%emis =   0.103E+04
        nerup_volc(31)%name = "                    COLIMA"
        nerup_volc(31)%lat  =     19.514
        nerup_volc(31)%lon  =   -103.620
        nerup_volc(31)%alt  =   4100.000
        nerup_volc(31)%emis =   0.140E+03
        nerup_volc(32)%name = "                     SANTA"
        nerup_volc(32)%lat  =     14.756
        nerup_volc(32)%lon  =    -91.552
        nerup_volc(32)%alt  =   3772.000
        nerup_volc(32)%emis =   0.230E+03
        nerup_volc(33)%name = "                     FUEGO"
        nerup_volc(33)%lat  =     14.473
        nerup_volc(33)%lon  =    -90.880
        nerup_volc(33)%alt  =   3763.000
        nerup_volc(33)%emis =   0.640E+03
        nerup_volc(34)%name = "                    PACAYA"
        nerup_volc(34)%lat  =     14.381
        nerup_volc(34)%lon  =    -90.601
        nerup_volc(34)%alt  =   2552.000
        nerup_volc(34)%emis =   0.510E+03
        nerup_volc(35)%name = "                     SANTA"
        nerup_volc(35)%lat  =     13.853
        nerup_volc(35)%lon  =    -89.630
        nerup_volc(35)%alt  =   2365.000
        nerup_volc(35)%emis =   0.200E+02
        nerup_volc(36)%name = "                    IZALCO"
        nerup_volc(36)%lat  =     13.813
        nerup_volc(36)%lon  =    -89.633
        nerup_volc(36)%alt  =   1950.000
        nerup_volc(36)%emis =   0.200E+02
        nerup_volc(37)%name = "                       SAN"
        nerup_volc(37)%lat  =     12.702
        nerup_volc(37)%lon  =    -87.004
        nerup_volc(37)%alt  =   1745.000
        nerup_volc(37)%emis =   0.590E+03
        nerup_volc(38)%name = "                    TELICA"
        nerup_volc(38)%lat  =     12.603
        nerup_volc(38)%lon  =    -86.845
        nerup_volc(38)%alt  =   1010.000
        nerup_volc(38)%emis =   0.157E+03
        nerup_volc(39)%name = "                    MASAYA"
        nerup_volc(39)%lat  =     11.984
        nerup_volc(39)%lon  =    -86.161
        nerup_volc(39)%alt  =    635.000
        nerup_volc(39)%emis =   0.790E+03
        nerup_volc(40)%name = "                    ARENAL"
        nerup_volc(40)%lat  =     10.463
        nerup_volc(40)%lon  =    -84.703
        nerup_volc(40)%alt  =   1657.000
        nerup_volc(40)%emis =   0.110E+03
        nerup_volc(41)%name = "                      POAS"
        nerup_volc(41)%lat  =     10.200
        nerup_volc(41)%lon  =    -84.233
        nerup_volc(41)%alt  =   2708.000
        nerup_volc(41)%emis =   0.500E+03
        nerup_volc(42)%name = "                      RUIZ"
        nerup_volc(42)%lat  =      4.895
        nerup_volc(42)%lon  =    -75.323
        nerup_volc(42)%alt  =   5321.000
        nerup_volc(42)%emis =   0.190E+04
        nerup_volc(43)%name = "                   GALERAS"
        nerup_volc(43)%lat  =      1.220
        nerup_volc(43)%lon  =    -77.370
        nerup_volc(43)%alt  =   4276.000
        nerup_volc(43)%emis =   0.650E+03
        nerup_volc(44)%name = "                    LASCAR"
        nerup_volc(44)%lat  =    -23.370
        nerup_volc(44)%lon  =    -67.730
        nerup_volc(44)%alt  =   5592.000
        nerup_volc(44)%emis =   0.240E+04
        nerup_volc(45)%name = "                KVERKFJOLL"
        nerup_volc(45)%lat  =     64.650
        nerup_volc(45)%lon  =    -16.720
        nerup_volc(45)%alt  =   1920.000
        nerup_volc(45)%emis =   0.300E+01
!--------------------------------------------------------------------------

        id=size(lon,1); jd=size(lat,2); kd=size(zhalf,3)-1

        SO2_nerup_volc_emis(:,:,:)=0.
        latmin = minval(lat)*180./pi-1.0
        latmax = maxval(lat)*180./pi+1.0
        lonmin = minval(lon)*180./pi-1.25
        lonmax = maxval(lon)*180./pi+1.25

        do iv=1,nerup
          if (nerup_volc(iv)%lon .lt.0.) nerup_volc(iv)%lon = &
            360.+nerup_volc(iv)%lon
          if (nerup_volc(iv)%lat.ge.latmin .and. &
              nerup_volc(iv)%lat.le.latmax .and. &
              nerup_volc(iv)%lon.ge.lonmin .and. &
              nerup_volc(iv)%lon.le.lonmax ) then 
            ilon=int((nerup_volc(iv)%lon-(lonmin+1.25))/2.5) + 1
            ilon=max(1,min(id,ilon))
            ilat=int((nerup_volc(iv)%lat-(latmin+1.0))/2.) + 1
            ilat=max(1,min(jd,ilat))
            do k = kd, 1, -1
              if (zhalf(ilon,ilat,k+1).le. nerup_volc(iv)%alt .and.&
                zhalf(ilon,ilat,k).gt. nerup_volc(iv)%alt) then
! Convert the unit from MgSO2/box/day to kgSO2/box/sec.                    
                SO2_nerup_volc_emis(ilon,ilat,k) =  &
                   SO2_nerup_volc_emis(ilon,ilat,k) +&
                   nerup_volc(iv)%emis * 1000. / (24.*3600.)
                exit
              endif
            enddo
          endif
        enddo

      if (id_SO2_nerup_volc_emis > 0) then
        used = send_data ( id_SO2_nerup_volc_emis, &
              SO2_nerup_volc_emis* WTM_S /WTM_SO2, &
              Time,is_in=is,js_in=js,ks_in=1)
      endif
!shm
    !if trim(runtype) is not natural or anthrop set narua_nerup_volc_emis ==0.0
     if(trim(runtype) == 'biomass_only' .or.  &
        trim(runtype) == 'fossil_fuels_only') then
       SO2_nerup_volc_emis =0.0
     else if (trim(runtype)== 'anthrop' .or.  &
        trim(runtype)== 'natural_only')then
     else 
        call error_mesg ('atmos_sulfate_mod', &
             'runtype is not recognized', FATAL)
     endif

end subroutine get_SO2_nerup_volc_emis
!</SUBROUTINE>

end module atmos_sulfate_mod
