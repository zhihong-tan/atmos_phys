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
!  This implies that the "atmos_sulfate_init" should be executed at the begining
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
use           interpolator_mod, only:  interpolate_type, interpolator_init, &
                                       interpolator, interpolator_end,     &
                                       CONSTANT, INTERP_WEIGHTED_P
use              constants_mod, only : PI, GRAV, RDGAS, WTMAIR

implicit none

private
!-----------------------------------------------------------------------
!----- interfaces -------
!
public  atmos_sulfate_init, atmos_sulfate_end, &
        atmos_DMS_emission, atmos_SOx_emission, atmos_SOx_chem

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

real , parameter :: WTM_S     = 32.0
real , parameter :: WTM_O3    = 48.0
real , parameter :: WTM_SO2   = 64.0
real , parameter :: WTM_SO4   = 96.0
real , parameter :: WTM_NH4_2SO4   = 132.00
real , parameter :: WTM_DMS   = 62.0
real , parameter :: WTM_MSA   = 96.0

!--- identification numbers for  diagnostic fields and axes ----
integer ::   id_OH                  = 0
integer ::   id_HO2                 = 0
integer ::   id_NO3                 = 0
integer ::   id_jH2O2               = 0
integer ::   id_O3                  = 0
integer ::   id_pH                  = 0

integer ::   id_DMSo                = 0
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

type(interpolate_type),save         ::  gas_conc_interp
type(interpolate_type),save         ::  aerocom_emission_interp
type(interpolate_type),save         ::  gocart_emission_interp
type(interpolate_type),save         ::  aircraft_emission_interp

character(len=20)  :: runtype = "gocart"

character(len=32)  :: gas_conc_filename = 'gas_conc_3D.nc'
character(len=32), dimension(6) :: gas_conc_name
data gas_conc_name/'OH','HO2','NO3','O3','jH2O2','pH'/

character(len=32)  :: gocart_emission_filename = 'gocart_emission.nc'
character(len=32), dimension(6) :: gocart_emission_name
data gocart_emission_name/'DMSo','SO2_GEIA1','SO2_GEIA2', &
                       'SO4_GEIA1','SO4_GEIA2','SO2_biobur'/

character(len=32)  :: aerocom_emission_filename = 'aerocom_emission.nc'
integer, parameter :: max_aerocom_emission=18
character(len=32), dimension(max_aerocom_emission)  :: aerocom_emission_name
data aerocom_emission_name/ &
                   'SO2_RoadTransport',&
                   'SO2_Off-road', &
                   'SO2_Domestic', &
                   'SO2_Industry', &
                   'SO2_International_Shipping', &
                   'SO2_Powerplants', &
                   'SO2_cont_volc', &
                   'alt_cont_volc_low', &
                   'alt_cont_volc_high', &
                   'SO2_expl_volc', &
                   'alt_expl_volc_low', &
                   'alt_expl_volc_high', &
                   'GFED_SO2_l1', &
                   'GFED_SO2_l2', &
                   'GFED_SO2_l3', &
                   'GFED_SO2_l4', &
                   'GFED_SO2_l5', &
                   'GFED_SO2_l6'/

character(len=32)  :: aircraft_emission_filename = 'aircraft_emission.nc'
character(len=32)  :: aircraft_emission_name(1)
data aircraft_emission_name/'fuel'/

namelist /simple_sulfate_nml/  &
                runtype,                         &
                gas_conc_filename, gas_conc_name,          &
                aerocom_emission_filename, aerocom_emission_name,  &
                gocart_emission_filename, gocart_emission_name,  &
                aircraft_emission_filename, aircraft_emission_name

!trim(runtype) 
!biomass_only; fossil_fuels_only, natural_only, anthrop

logical :: module_is_initialized=.FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id: atmos_sulfate.F90,v 15.0 2007/08/14 03:57:00 fms Exp $'
character(len=128) :: tagname = '$Name: omsk_2007_10 $'
!-----------------------------------------------------------------------

contains


!#######################################################################

!<SUBROUTINE NAME="atmos_sulfate_init">
!<OVERVIEW>
! The constructor routine for the sulfate module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to initialize the sulfate module.
!</DESCRIPTION>
 subroutine atmos_sulfate_init ( lonb, latb, nlev, axes, Time, mask)

!-----------------------------------------------------------------------
real, intent(in),    dimension(:,:)                 :: lonb, latb
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
      character*12 :: SOx_tracer(5)
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
        read  (unit, nml=simple_sulfate_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'simple_sulfate_nml')
        end do
10      call close_file (unit)
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
                          write (stdlog(), nml=simple_sulfate_nml)


!----- set initial value of sulfate ------------

     do m=1,size(SOx_tracer)

       n = get_tracer_index(MODEL_ATMOS,trim(SOx_tracer(m)))
       if (n>0) then
         nsulfate=n
         if (nsulfate > 0 .and. mpp_pe() == mpp_root_pe()) &
                 write (stdlog(),30) trim(SOx_tracer(m)),nsulfate
       endif
     enddo


  30   format (A,' was initialized as tracer number ',i2)
     call interpolator_init (gas_conc_interp, trim(gas_conc_filename),  &
                             lonb, latb,&
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = gas_conc_name, &
                             vert_interp=(/INTERP_WEIGHTED_P/) )
     call interpolator_init (aerocom_emission_interp, &
                             trim(aerocom_emission_filename),  &
                             lonb, latb,&
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = aerocom_emission_name, &
                             vert_interp=(/INTERP_WEIGHTED_P/) )

     call interpolator_init (gocart_emission_interp, &
                             trim(gocart_emission_filename),  &
                             lonb, latb,&
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = gocart_emission_name, &
                             vert_interp=(/INTERP_WEIGHTED_P/) )

     call interpolator_init (aircraft_emission_interp, &
                             trim(aircraft_emission_filename),  &
                             lonb, latb,&
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = aircraft_emission_name, &
                             vert_interp=(/INTERP_WEIGHTED_P/) )

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
      id_DMSo       = register_diag_field ( mod_name,             &
                      'DMSo',axes(1:2),Time,                      &
                      'Dimethylsulfide seawater concentration',   &
                      'nM/L')
      id_ph          = register_diag_field ( mod_name,             &
                      'pH_simple_sulfate',axes(1:3),Time,                    &
                      'pH in simple-sulfate',                  &
                      'none')
      id_O3           = register_diag_field ( mod_name,             &
                      'O3_simple_sulfate',axes(1:3),Time,                    &
                      'O3 in simple-sulfate',                  &
                      'none')

      id_SO2_anth_l1_emis = register_diag_field ( mod_name,     &
                      'simpleSO2_anth_l1_emis',axes(1:2),Time,          &
                      'SO2 anthropogenic emission GEIA-85 level1',&
                      'kgS/m2/s')
      id_SO2_anth_l2_emis = register_diag_field ( mod_name,     &
                      'simpleSO2_anth_l2_emis',axes(1:2),Time,          &
                      'SO2 anthropogenic emission GEIA-85 level2',&
                      'kgS/m2/s')
      id_SO2_aircraft_emis= register_diag_field ( mod_name,     &
                      'simpleSO2_aircraft_emis',axes(1:3),Time,         &
                      'SO2 emission by aircraft',                 &
                      'kgS/m2/s')
      id_SO2_bioburn_emis = register_diag_field ( mod_name,     &
                      'simpleSO2_bioburn_emis',axes(1:2),Time,          &
                      'SO2 emission from biomass burning',        &
                      'kgS/m2/s')
      id_SO2_nerup_volc_emis = register_diag_field ( mod_name,     &
                      'simpleSO2_nerup_volc_emis',axes(1:3),Time,     &
                      'SO2 emission from non-eruptive volcanoes', &
                      'kgS/m2/s')
      id_SO4_anth_l1_emis = register_diag_field ( mod_name,     &
                      'simpleSO4_anth_l1_emis',axes(1:2),Time,          &
                      'SO4 anthropogenic emission GEIA-85 level1',&
                      'kgS/m2/s')
      id_SO4_anth_l2_emis = register_diag_field ( mod_name,     &
                      'simpleSO4_anth_l2_emis',axes(1:2),Time,          &
                      'SO4 anthropogenic emission GEIA-85 level2',&
                      'kgS/m2/s')
      id_NO3        = register_diag_field ( mod_name,           &
                      'simpleNO3_diurnal',axes(1:3),Time,       &
                      'Time varying NO3 concentration',         &
                      'molec.cm-3')

      id_OH         = register_diag_field ( mod_name,         &
                      'OH_simple_sulfate',axes(1:3),Time,    &
                      'Varying Hydroxyl radical concentration',           &
                      'molec.cm-3')
      id_jH2O2         = register_diag_field ( mod_name,           &
                      'jH2O2_simple_sulfate',axes(1:3),Time,               &
                      'Varying H2O2 photodissociation',   &
                      's-1')

      id_HO2         = register_diag_field ( mod_name,           &
                      'HO2_simple_sulfate',axes(1:3),Time,               &
                      'Varying Hydroperoxyl radical concentration',   &
                      'molec.cm-3')
      id_DMS_chem   = register_diag_field ( mod_name,           &
                      'simpleDMS_chem',axes(1:3),Time,                       &
                      'DMS chemical production',      &
                      'kgS/m2/s')
      id_SO2_chem   = register_diag_field ( mod_name,           &
                      'simpleSO2_chem',axes(1:3),Time,                       &
                      'SO2 chemical production',      &
                      'kgS/m2/s')
      id_SO4_chem   = register_diag_field ( mod_name,           &
                      'simpleSO4_chem',axes(1:3),Time,                       &
                      'SO4 chemical production',      &
                      'kgS/m2/s')
      id_MSA_chem   = register_diag_field ( mod_name,           &
                      'simpleMSA_chem',axes(1:3),Time,                       &
                      'MSA chemical production',      &
                      'kgS/m2/s')
      id_H2O2_chem   = register_diag_field ( mod_name,           &
                      'simpleH2O2_chem',axes(1:3),Time,                       &
                      'H2O2 chemical production',      &
                      'kgH2O2/m2/s')

      call write_version_number (version, tagname)

      module_is_initialized = .TRUE.

!-----------------------------------------------------------------------
 end subroutine atmos_sulfate_init
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_sulfate_end">
!<OVERVIEW>
!  The destructor routine for the sulfate module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_sulfate_end
!</TEMPLATE>
 subroutine atmos_sulfate_end

        call interpolator_end (gas_conc_interp) 
        call interpolator_end (aerocom_emission_interp) 
        call interpolator_end (gocart_emission_interp) 
        call interpolator_end (aircraft_emission_interp) 
        module_is_initialized = .FALSE.

 end subroutine atmos_sulfate_end
!</SUBROUTINE>
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
       pwt, DMS_dt, Time, is,ie,js,je,kbot)
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
      real, dimension(size(DMS_dt,1),size(DMS_dt,2)) :: DMSo, DMS_emis
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
        DMSo(:,:)=0.0
        call interpolator(gocart_emission_interp, Time, DMSo, &
                       trim(gocart_emission_name(1)), is, js)
! --- Send the DMS data to the diag_manager for output.
        if (id_DMSo > 0 ) &
          used = send_data ( id_DMSo, DMSo, Time, is_in=is, js_in=js )

        CONC = DMSo(i,j)

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
subroutine atmos_SOx_emission (lon, lat, area, frac_land, &
       z_pbl, zhalf, phalf, pwt, SO2_dt, SO4_dt, Time, is,ie,js,je,kbot)
!
! This subroutine calculates the tendencies of SO2 and SO4 due to
! their emissions.
! The inventories are based from AEROCOM (cf. Dentener, ACPD, 2006)
! except the aircraft emission.
! The emission of SO4 is assumed to be fe=2.5% of all sulfur emission
! (cf. Dentener, ACPD, 2006). NB. Some authors consider 5%
!
      real, intent(in),    dimension(:,:)           :: lon, lat
      real, intent(in),    dimension(:,:)           :: frac_land
      real, intent(in),    dimension(:,:)           :: area
      real, intent(in),    dimension(:,:)           :: z_pbl
      real, intent(in),    dimension(:,:,:)         :: zhalf, phalf
      real, intent(in),    dimension(:,:,:)         :: pwt
      real, intent(out),   dimension(:,:,:)         :: SO2_dt, SO4_dt
      type(time_type), intent(in)                   :: Time
      integer, intent(in)                           :: is, ie, js, je
      integer, intent(in), dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      integer, parameter :: nlevel_fire = 6
      real, dimension(size(SO4_dt,1),size(SO4_dt,2),size(SO4_dt,3)) :: SO4_emis
      real, dimension(size(SO2_dt,1),size(SO2_dt,2),size(SO2_dt,3)) :: SO2_emis
! Input emission fields
      real, dimension(size(SO2_dt,1),size(SO2_dt,2),size(SO2_dt,3)) :: &
             SO2_aircraft
      real, dimension(size(SO2_dt,1),size(SO2_dt,2)) :: &
             SO2_GEIA1, SO2_GEIA2, SO4_GEIA1, SO4_GEIA2,&
             SO2_biobur,                                &
             SO2_RoadTransport,                         &
             SO2_Off_road,                              &
             SO2_Domestic,                              &
             SO2_Industry,                              &
             SO2_International_Shipping,                &
             SO2_Powerplants,                           &
             SO2_cont_volc,                             &
             SO2_expl_volc,                             &
             SO2_cont_volc_h1,                          &
             SO2_cont_volc_h2,                          &
             SO2_expl_volc_h1,                          &
             SO2_expl_volc_h2
      real, dimension(size(SO2_dt,1),size(SO2_dt,2),nlevel_fire) :: &
             SO2_wildfire
! Factors of vertical distribution of emissions
      real, dimension(size(SO2_dt,3)) :: fbb, fa1, fa2, fcv, fev
      real, dimension(size(SO2_dt,3),nlevel_fire) :: ff
! Lower altitude of injection of SO2 from wild fires 
! These values correspond to the AEROCOM input data (cf. Dentener, ACPD, 2006)
      real, dimension(nlevel_fire) :: &
             alt_fire_min=(/0.,100.,500.,1000.,2000.,3000./)
! Upper altitude of injection of SO2 from wild fires 
! These values correspond to the AEROCOM input data (cf. Dentener, ACPD, 2006)
      real, dimension(nlevel_fire) :: &
             alt_fire_max=(/100.,500.,1000.,2000.,3000.,6000./)
! Altitude of injection of surafce anthropogenic emissions
      real :: ze1
! Altitude of injection of SO2 from industries and power plants.   
      real :: ze2
! Emission factor for SO4
      real, parameter :: fe = 0.025

      real :: z1, z2, bltop, fbt, del
      integer  :: i, j, k, l, id, jd, kd, il, lf

      id=size(SO2_dt,1); jd=size(SO2_dt,2); kd=size(SO2_dt,3)
!
! Initialize
!
      SO2_dt(:,:,:) = 0.0
      SO4_dt(:,:,:) = 0.0
! GOCART emissions
      SO2_GEIA1(:,:)=0.0
      SO2_GEIA2(:,:)=0.0
      SO4_GEIA1(:,:)=0.0
      SO4_GEIA2(:,:)=0.0
      SO2_biobur(:,:)=0.0
! AEROCOM emissions
      SO2_RoadTransport(:,:)=0.0
      SO2_Off_road(:,:)=0.0
      SO2_Domestic(:,:)=0.0
      SO2_Industry(:,:)=0.0
      SO2_International_Shipping(:,:)=0.0
      SO2_Powerplants(:,:)=0.0
      SO2_aircraft(:,:,:)=0.0
      SO2_cont_volc(:,:)=0.0
      SO2_expl_volc(:,:)=0.0
      SO2_wildfire(:,:,:)=0.0

      SO2_cont_volc_h1(:,:)=0.0
      SO2_cont_volc_h2(:,:)=0.0
      SO2_expl_volc_h1(:,:)=0.0
      SO2_expl_volc_h2(:,:)=0.0
!
      if ( runtype .eq. 'gocart') then
        call interpolator(gocart_emission_interp, Time, SO2_GEIA1, &
                       trim(gocart_emission_name(2)), is, js)
        call interpolator(gocart_emission_interp, Time, SO2_GEIA2, &
                       trim(gocart_emission_name(3)), is, js)
        call interpolator(gocart_emission_interp, Time, SO4_GEIA1, &
                       trim(gocart_emission_name(4)), is, js)
        call interpolator(gocart_emission_interp, Time, SO4_GEIA2, &
                       trim(gocart_emission_name(5)), is, js)
        call interpolator(gocart_emission_interp, Time, SO2_biobur, &
                       trim(gocart_emission_name(6)), is, js)
      endif
      if ( runtype .eq. 'aerocom') then
        call interpolator(aerocom_emission_interp, Time, SO2_RoadTransport, &
                       trim(aerocom_emission_name(1)), is, js)
!
        call interpolator(aerocom_emission_interp, Time, SO2_Off_road, &
                       trim(aerocom_emission_name(2)), is, js)
!
        call interpolator(aerocom_emission_interp, Time, SO2_Domestic, &
                       trim(aerocom_emission_name(3)), is, js)
!
        call interpolator(aerocom_emission_interp, Time, SO2_Industry, &
                       trim(aerocom_emission_name(4)), is, js)
!
        call interpolator(aerocom_emission_interp, Time, &
                       SO2_International_Shipping, &
                       trim(aerocom_emission_name(5)), is, js)
!
        call interpolator(aerocom_emission_interp, Time, SO2_Powerplants, &
                       trim(aerocom_emission_name(6)), is, js)
! Wildfire emissions at 6 levels from 0 to 6 km
! (cf. AEROCOM web site or Dentener et al., ACPD, 2006)
        do il=1,nlevel_fire
          call interpolator(aerocom_emission_interp, Time, SO2_wildfire(:,:,il), &
                       trim(aerocom_emission_name(12+il)), is, js)
        enddo
      endif
!
! Aircraft emissions
!
      call interpolator(aircraft_emission_interp, Time, phalf, SO2_aircraft, &
                       trim(aircraft_emission_name(1)), is, js)
!
! Continuous volcanoes
!
      call interpolator(aerocom_emission_interp, Time, SO2_cont_volc, &
                       trim(aerocom_emission_name(7)), is, js)
      call interpolator(aerocom_emission_interp, Time, SO2_cont_volc_h1, &
                       trim(aerocom_emission_name(8)), is, js)
      call interpolator(aerocom_emission_interp, Time, SO2_cont_volc_h2, &
                       trim(aerocom_emission_name(9)), is, js)
!
! Explusive volcanoes
!
      call interpolator(aerocom_emission_interp, Time, SO2_expl_volc, &
                       trim(aerocom_emission_name(10)), is, js)
      call interpolator(aerocom_emission_interp, Time, SO2_expl_volc_h1, &
                       trim(aerocom_emission_name(11)), is, js)
      call interpolator(aerocom_emission_interp, Time, SO2_expl_volc_h2, &
                       trim(aerocom_emission_name(12)), is, js)

      do j = 1, jd
      do i = 1, id

! --- Assuming biomass burning emission within the PBL -------
        fbb(:) = 0.
        if (runtype.eq.'gocart') then
          ze1=100.
          ze2=500.
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
        endif
! --- Assuming anthropogenic source L1 emitted below Ze1, and L2
!     emitted between Ze1 and Ze2.
        ff(:,:)=0.
        if (runtype.eq.'aerocom') then
          ze1=100.
          ze2=300.
          do l = kd,2,-1
            Z1 = zhalf(i,j,l+1)-zhalf(i,j,kd+1)
            Z2 = zhalf(i,j,l)-zhalf(i,j,kd+1)
            do lf=1,nlevel_fire
              del=alt_fire_max(lf)-alt_fire_min(lf)
              if (del.gt.0. .and. Z1.le.alt_fire_max(lf) ) then
                if (Z2.lt.alt_fire_max(lf)) then
                  ff(l,lf) = 0.
                else
                  if (Z1.lt.alt_fire_min(lf)) then
                    ff(l,lf) = (Z2-alt_fire_min(lf))/del
                  else
                    if (Z2.lt.alt_fire_max(lf)) then
                      ff(l,lf) = (z2-z1)/del
                    else
                      ff(l,lf)=(alt_fire_max(lf)-Z1)/del
                    endif
                  endif
                endif
              endif
            enddo
          enddo
        endif
! --- For continuous volcanoes, calculate the fraction of emission for
! --- each vertical levels
        fcv(:)=0.
        do l = kd,2,-1
          Z1 = zhalf(i,j,l+1)-zhalf(i,j,kd+1)
          Z2 = zhalf(i,j,l)-zhalf(i,j,kd+1)
          del=SO2_cont_volc_h2(i,j)-SO2_cont_volc_h1(i,j)
          if (del.gt.0. .and. Z1.lt.SO2_cont_volc_h2(i,j) ) then
            if (Z2.lt.SO2_cont_volc_h1(i,j)) then
              fcv(l) = 0.
            else
              if (Z1.lt.SO2_cont_volc_h1(i,j)) then
                fcv(l) = (Z2-SO2_cont_volc_h1(i,j))/del
              else
                if (Z2.lt.SO2_cont_volc_h2(i,j)) then
                  fcv(l)=(z2-z1)/del
                else
                  fcv(l)=(SO2_cont_volc_h2(i,j)-Z1)/del
                endif
              endif
            endif
          endif
        enddo

! --- For explosive volcanoes, calculate the fraction of emission for
! --- each vertical levels
        fev(:)=0.
        do l = kd,2,-1
          Z1 = zhalf(i,j,l+1)-zhalf(i,j,kd+1)
          Z2 = zhalf(i,j,l)-zhalf(i,j,kd+1)
          del=SO2_expl_volc_h2(i,j)-SO2_expl_volc_h1(i,j)
          if (del.gt.0. .and. Z1.lt.SO2_expl_volc_h2(i,j) ) then
            if (Z2.lt.SO2_expl_volc_h1(i,j)) then
              fev(l) = 0.
            else
              if (Z1.lt.SO2_expl_volc_h1(i,j)) then
                fev(l) = (Z2-SO2_expl_volc_h1(i,j))/del
              else
                if (Z2.lt.SO2_expl_volc_h2(i,j)) then
                  fev(l)=(z2-z1)/del
                else
                  fev(l)=(SO2_expl_volc_h2(i,j)-Z1)/del
                endif
              endif
            endif
          endif
        enddo
! --- For fosil fuel emissions, calculate the fraction of emission for
! --- each vertical levels
        fa1(:) = 0.
        fa2(:) = 0.
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
! SO2_emis: [kgSO2/m2/s]
        SO2_emis(i,j,:) =( &
!              Assuming that 1g of SO2 is emitted from 1kg of fuel: 1.e-3
               1.e-3  * SO2_aircraft(i,j,:)                   &
             + fa1(:) * SO2_GEIA1(i,j)                        &
             + fa2(:) * SO2_GEIA2(i,j)                        &
             + fbb(:) * SO2_biobur(i,j)                       &
             + fa1(:) * SO2_RoadTransport(i,j)                &
             + fa1(:) * SO2_Off_road(i,j)                     &
             + fa1(:) * SO2_Domestic(i,j)                     &
             + fa1(:) * SO2_International_Shipping(i,j)       &
             + fa2(:) * SO2_Industry(i,j)                     &
             + fa2(:) * SO2_Powerplants(i,j)                  &
             + fcv(:) * SO2_cont_volc(i,j)                    &
             + fev(:) * SO2_expl_volc(i,j)                    )
        do lf = 1, nlevel_fire
          SO2_emis(i,j,:) = SO2_emis(i,j,:) + ff(:,lf) * SO2_wildfire(i,j,lf)
        enddo
!
! Aerocom assumes a constant emission index for sulfate (2.5%)
        if (runtype .eq. 'aerocom') then
          SO4_emis(i,j,:) = fe * SO2_emis(i,j,:)
          SO2_emis(i,j,:) = (1.-fe)* SO2_emis(i,j,:)
        endif
!
! GOCART assumes continent based emission index for sulfate:
!    Anthropogenic SOx emission from GEIA 1985.
!    Assuming:   Europe:      5.0% SOx emission is SO4;
!                US + Canada: 1.4% SOx emission is SO4;
!                The rest:    2.5% SOx emission is SO4.
        if (runtype .eq. 'gocart' ) then
          SO4_emis(i,j,:) = &
               fa1(:) * SO4_GEIA1(i,j) + fa2(:) * SO4_GEIA2(i,j)
        endif

      end do   ! end i loop
      end do   ! end j loop
!
      SO2_dt(:,:,:)= SO2_emis(:,:,:)/pwt(:,:,:)*WTMAIR/WTM_SO2
      SO4_dt(:,:,:)= SO4_emis(:,:,:)/pwt(:,:,:)*WTMAIR/WTM_SO4

!------------------------------------------------------------------
! DIAGNOSTICS:      SO2 and SO4 emission in kg/timestep
!--------------------------------------------------------------------
      if (id_SO2_emis > 0) then
        used = send_data ( id_SO2_emis, SO2_emis*WTM_S/WTM_SO2, Time, &
              is_in=is,js_in=js,ks_in=1)
      endif
      if (id_SO4_emis > 0) then
        used = send_data ( id_SO4_emis, SO4_emis*WTM_S/WTM_SO4, Time, &
              is_in=is,js_in=js,ks_in=1)
      endif

end subroutine atmos_SOx_emission
!</SUBROUTINE>
!-----------------------------------------------------------------------
!#######################################################################
      subroutine atmos_SOx_chem(pwt,temp,pfull, phalf, dt, lwc, &
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
      real, intent(in), dimension(:,:,:) :: temp, pfull, phalf
      real, intent(in), dimension(:,:,:) :: SO2, SO4, DMS, MSA, H2O2
      real, intent(out),dimension(:,:,:) :: SO2_dt,SO4_dt,DMS_dt,MSA_dt,H2O2_dt

      type(time_type), intent(in)                    :: Time
      integer, intent(in),  dimension(:,:), optional :: kbot
      integer, intent(in)                            :: is,ie,js,je
! Working vectors
      integer :: i,j,k,id,jd,kd, istop
      integer                                    :: istep, nstep
!!! Input fields from interpolator
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: pH
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: O3_mmr
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: no3_conc
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: oh_conc
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: jh2o2
      real, dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: ho2_conc
!!! Time varying fields
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

      OH_conc(:,:,:)=1.e5  ! molec/cm3
      call interpolator(gas_conc_interp, Time, phalf, OH_conc, &
                       trim(gas_conc_name(1)), is, js)

      HO2_conc(:,:,:)=1.e6  ! molec/cm3
      call interpolator(gas_conc_interp, Time, phalf, HO2_conc, &
                       trim(gas_conc_name(2)), is, js)

      NO3_conc(:,:,:)=0.0  ! molec/cm3
      call interpolator(gas_conc_interp, Time, phalf, NO3_conc, &
                       trim(gas_conc_name(3)), is, js)

      O3_mmr(:,:,:)=0  ! Ozone mass mixing ratio
      call interpolator(gas_conc_interp, Time, phalf, O3_mmr, &
                       trim(gas_conc_name(4)), is, js)
      O3_mmr(:,:,:)=O3_mmr(:,:,:)*WTM_O3/WTMAIR

      jH2O2(:,:,:)=1.e-6 ! s-1
      call interpolator(gas_conc_interp, Time, phalf, jH2O2, &
                       trim(gas_conc_name(5)), is, js)

      pH(:,:,:)=1.e-5
      call interpolator(gas_conc_interp, Time, phalf, pH, &
                       trim(gas_conc_name(6)), is, js)

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
       xph   = max(1.e-7,       pH(i,j,k))
       xoh   = max(0.         , OH_conc(i,j,k)  *fac_OH(i,j))
       xho2  = max(0.         , HO2_conc(i,j,k) *fac_HO2(i,j))
       xjh2o2= max(0.         , jH2O2(i,j,k)    *fac_OH(i,j))
       xno3  = max(0.         , NO3_conc(i,j,k) *fac_NO3(i,j))
       xo3   = max(small_value, O3_mmr(i,j,k))
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
      if ( id_NO3 > 0) then
        used = send_data ( id_NO3, NO3_diurnal, &
                           Time,is_in=is,js_in=js,ks_in=1)
      endif
      if ( id_OH > 0) then
        used = send_data ( id_OH, OH_diurnal, &
                           Time, is_in=is, js_in=js,ks_in=1 )
      endif
      if ( id_HO2 > 0) then
        used = send_data ( id_HO2, HO2_diurnal, &
                           Time, is_in=is, js_in=js,ks_in=1 )
      endif
      if ( id_jH2O2 > 0) then
        used = send_data ( id_jH2O2, jH2O2_diurnal, &
                           Time, is_in=is, js_in=js,ks_in=1 )
      endif
      if (id_ph > 0) then
        used = send_data ( id_ph, ph, &
                           Time,is_in=is,js_in=js,ks_in=1)
      endif
      if (id_o3 > 0) then
        used = send_data ( id_o3, o3_mmr, &
                           Time,is_in=is,js_in=js,ks_in=1)
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
end subroutine atmos_SOx_chem

end module atmos_sulfate_mod
