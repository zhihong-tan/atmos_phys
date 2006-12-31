Module atmos_carbon_aerosol_mod

! <CONTACT EMAIL="Shekar.Reddy@noaa.gov">
!   Shekar Reddy
! </CONTACT>
use fms_mod,           only: open_namelist_file, fms_init, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             file_exist, write_version_number, &
                             check_nml_error, error_mesg, &
                             FATAL, NOTE, WARNING, close_file
use     time_manager_mod, only : time_type
use     diag_manager_mod, only : send_data,            &
                                 register_diag_field,  &
                                 register_static_field

use   tracer_manager_mod, only : get_tracer_index, &
                                 set_tracer_atts
use    field_manager_mod, only : MODEL_ATMOS
use     time_manager_mod, only : time_type, get_date
use     interpolator_mod, only : interpolate_type, interpolator_init, &
                                 interpolator, interpolator_end, &
                                 CONSTANT, INTERP_WEIGHTED_P
use        constants_mod, only : PI,RDGAS

implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_carbon_aerosol_driver,   &
        atmos_carbon_aerosol_init, &
        atmos_carbon_aerosol_end

!-----------------------------------------------------------------------
! tracer number for carbonaceous aerosols
integer :: nbcphobic=0
integer :: nbcphilic=0
integer :: nomphobic=0
integer :: nomphilic=0

!--- identification numbers for  diagnostic fields and axes ----
integer :: id_bcphob_sink, id_bcphil_sink
integer :: id_omphob_sink, id_omphil_sink
!
integer :: id_bcemisff_l, id_bcemisbb_l
integer :: id_omemisff_l, id_omemisbb_l
integer :: id_bcemisff_h, id_bcemisbb_h
integer :: id_omemisff_h, id_omemisbb_h
integer :: id_omemisnat
!----------------------------------------------------------------------
type(interpolate_type),save         ::carbon_aerosol_interp
!----------------------------------------------------------------------
!-------- namelist  ---------
character(len=26)  :: emission_filename = 'carbon_aerosol_emission.nc'
character(len=12), dimension(15)  :: emission_name
   data emission_name/'bcff_cooke99','omff_cooke99', &
                      'bcbb_cw96','bcff_cw96', &
                      'bcff_bond','omff_bond','bcbf_bond','ombf_bond', &
                      'bcob_bond','omob_bond', &
                      'bc_geia1','bc_geia2','om_geia1','om_geia2', &
                      'omemisnat' /
character(len=12)  :: emission
! Default values for carbon_aerosol_nml
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FOSSIL FUEL source can be either:
! for black carbon: 'cooke_and_wilson_1996' 
!                   'cooke_1999' 
!                   'bond_2004'
character(len=9)  :: bcff_source = 'bond_2004'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fossil fuel emission 
! for organic matter: 'cooke_1999' 
!                     'bond_2004'
character(len=20)  :: omff_source = 'bond_2004'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BIOMASS BURNING emissions can be either
! for black carbon: 'cooke_and_wilson_1996' 
!                   'bond_2004' 
!                   'reddy_and_boucher_2004'
character(len=20)  :: bcbb_source = 'bond_2004'
! biomass burning emission
! for organic matter: 'cooke_and_wilson_1996' 
!                     'bond_2004' 
!                     'reddy_and_boucher_2004'
character(len=20)  :: ombb_source = 'Bond_2004'
! for natural emssion, there is no choice (and thus no selection): Gunther et al., 1996 inventory
namelist /carbon_aerosol_nml/ &
          bcff_source, &
          omff_source, &
          bcbb_source, &
          ombb_source,  &
          emission_filename, emission_name

character(len=6), parameter :: module_name = 'tracer'

logical :: module_is_initialized = .FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id: atmos_carbon_aerosol.F90,v 13.0.2.1 2006/09/12 21:11:09 wfc Exp $'
character(len=128) :: tagname = '$Name: memphis_2006_12 $'
!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine atmos_carbon_aerosol_driver(lon, lat, land, pfull,phalf,z_half, z_pbl, &
                               T, pwt, &
                               black_cphob, black_cphob_dt,  &
                               black_cphil, black_cphil_dt,  &
                               omphob, omphob_dt, &
                               omphil, omphil_dt, &
                               Time, is, ie, js, je )

!-----------------------------------------------------------------------
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:)   :: land
   real, intent(in),  dimension(:,:)   :: z_pbl
   real, intent(in),  dimension(:,:,:) :: z_half
   real, intent(in),  dimension(:,:,:) :: pwt,pfull,phalf,T
   real, intent(in),  dimension(:,:,:) :: black_cphob,black_cphil
   real, intent(in),  dimension(:,:,:) :: omphob,omphil
   real, intent(out), dimension(:,:,:) :: black_cphob_dt,black_cphil_dt
   real, intent(out), dimension(:,:,:) :: omphob_dt,omphil_dt
type(time_type), intent(in)            :: Time
integer, intent(in)                    :: is, ie, js, je
!-----------------------------------------------------------------------
real, dimension(size(black_cphob,1),size(black_cphob,2),size(black_cphob,3)) ::&
   bc_sourcephob, bc_sourcephil, om_sourcephob, om_sourcephil, sinkphob
   real  dtr,bltop,z1,z2
real, dimension(size(black_cphob,3)) :: fbb, fff
integer  i,j,l, kb,id,jd,kd,lat1,k
REAL,PARAMETER :: frac_bc_phobic = 0.8
REAL,PARAMETER :: frac_bc_philic = 0.2
REAL,PARAMETER :: frac_om_phobic = 0.5
REAL,PARAMETER :: frac_om_philic = 0.5

!
!-------------------------------------------------------------------------
!Emissions profiles as per SAFARI extinction 
!----
   real, dimension(size(black_cphob,1),size(black_cphob,2)) ::          &
              bcemisff_l, bcemisff_h, bcemisbb_l, bcemisbb_h,           &
              omemisff_l, omemisff_h, omemisbb_l, omemisbb_h, omemisnat

!
!-----------------------------------------------------------------------
!
      id=size(black_cphob,1); jd=size(black_cphob,2); kd=size(black_cphob,3)

      dtr= PI/180.

!----------- compute black carbon source ------------

    bc_sourcephob = 0.0
    bc_sourcephil = 0.0
    om_sourcephob = 0.0
    om_sourcephil = 0.0
    black_cphob_dt = 0.0
    black_cphil_dt = 0.0
    omphob_dt = 0.0
    omphil_dt = 0.0
! emission rates (kg/m2/s)
    bcemisff_l(:,:) = 0.0
    bcemisff_h(:,:) = 0.0
    omemisff_l(:,:) = 0.0
    omemisff_h(:,:) = 0.0
    bcemisbb_l(:,:) = 0.0
    bcemisbb_h(:,:) = 0.0
    omemisbb_l(:,:) = 0.0
    omemisbb_h(:,:) = 0.0
    omemisnat(:,:)  = 0.0

    select case (trim(bcff_source))
      case ('cooke_and_wilson_1996')
             emission = emission_name(4)
      case ('cooke_1999')
             emission = emission_name(1)
      case ('bond_2004') 
             emission = emission_name(5)
    end select
    call interpolator(carbon_aerosol_interp, Time, bcemisff_l, &
                      trim(emission), is, js)

    select case (trim(omff_source))
      case ('cooke_1999')
             emission = emission_name(2)
      case ('bond_2004') 
             emission = emission_name(6)
    end select
    call interpolator(carbon_aerosol_interp, Time, omemisff_l, &
                      trim(emission), is, js)

    select case (trim(bcbb_source))
      case ('cooke_and_wilson_1996')
        call interpolator(carbon_aerosol_interp, Time, bcemisbb_l, &
                         trim(emission_name(3)), is, js)
      case ('bond_2004') 
        call interpolator(carbon_aerosol_interp, Time, bcemisbb_l, &
                         trim(emission_name(7)), is, js)
        call interpolator(carbon_aerosol_interp, Time, bcemisbb_h, &
                         trim(emission_name(8)), is, js)
      case ('reddy_and_boucher_2004') 
        call interpolator(carbon_aerosol_interp, Time, bcemisbb_l, &
                         trim(emission_name(11)), is, js)
        call interpolator(carbon_aerosol_interp, Time, bcemisbb_h, &
                         trim(emission_name(12)), is, js)
    end select

    select case (trim(ombb_source))
      case ('cooke_and_wilson_1996')
        call interpolator(carbon_aerosol_interp, Time, omemisbb_l, &
                           trim(emission_name(3)), is, js)
        omemisbb_l(:,:) = omemisbb_l(:,:)*7.0
        omemisbb_h(:,:) = 0.0
      case ('bond_2004') 
        call interpolator(carbon_aerosol_interp, Time, omemisbb_l, &
                       trim(emission_name(8)), is, js)
        call interpolator(carbon_aerosol_interp, Time, omemisbb_h, &
                       trim(emission_name(10)), is, js)
      case ('reddy_and_boucher_2004') 
          call interpolator(carbon_aerosol_interp, Time, omemisbb_l, &
                       trim(emission_name(13)), is, js)
          call interpolator(carbon_aerosol_interp, Time, omemisbb_h, &
                       trim(emission_name(14)), is, js)
    end select

    call interpolator(carbon_aerosol_interp, Time, omemisnat, &
                       trim(emission_name(15)), is, js)
!
    do j = 1, jd
      do i = 1, id

        fbb(:) = 0.
        fbb(kd)= 1.
        fff(:) = 0.
        fff(kd-1)= 1.
!
! --- Assuming biomass burning emission within the PBL -------
!
        bltop = z_pbl(i,j)
        do l = kd,1,-1
          z1=z_half(i,j,l+1)-z_half(i,j,kd+1)
          z2=z_half(i,j,l)-z_half(i,j,kd+1)
          if (bltop.lt.z1) exit
          if (bltop.ge.z2) fbb(l)=(z2-z1)/z_pbl(i,j)
          if (bltop.gt.z1.and.bltop.lt.z2) then
            fbb(l) = (bltop-z1)/z_pbl(i,j)
          endif
        enddo
!
        bc_sourcephob(i,j,kd) =  frac_bc_phobic  &
          *(bcemisff_l(i,j)+ bcemisbb_l(i,j))/pwt(i,j,kd)

        bc_sourcephil(i,j,kd) =  frac_bc_philic  &
          *(bcemisff_l(i,j)+ bcemisbb_l(i,j))/pwt(i,j,kd)

        bc_sourcephob(i,j,:) =  bc_sourcephob(i,j,:) + frac_bc_phobic &
          *(fff(:)*bcemisff_h(i,j)+fbb(:)*bcemisbb_h(i,j)) /pwt(i,j,:)
        bc_sourcephil(i,j,:) =  bc_sourcephil(i,j,:) + frac_bc_philic &
          *(fff(:)*bcemisff_h(i,j)+fbb(:)*bcemisbb_h(i,j)) /pwt(i,j,:)
!
        om_sourcephob(i,j,kd) =  frac_om_phobic  &
          *(omemisff_l(i,j)+ omemisbb_l(i,j) &
            +omemisnat(i,j))/pwt(i,j,kd)

        om_sourcephil(i,j,kd) =  frac_om_philic  &
          *(omemisff_l(i,j)+ omemisbb_l(i,j) &
            +omemisnat(i,j))/pwt(i,j,kd)

        om_sourcephob(i,j,:) =  om_sourcephob(i,j,:) + frac_om_phobic &
          *(fff(:)*omemisff_h(i,j)+fbb(:)*omemisbb_h(i,j))/pwt(i,j,:)
        om_sourcephil(i,j,:) =  om_sourcephil(i,j,:) + frac_om_philic &
          *(fff(:)*omemisff_h(i,j)+fbb(:)*omemisbb_h(i,j))/pwt(i,j,:)
      enddo
    enddo

!------- compute black carbon phobic sink --------------
!
!  BCphob has a half-life time of 1.0days 
!   (corresponds to an e-folding time of 1.44 days)
!
!  sink = 1./(86400.*1.44) = 8.023e-6

    sinkphob(:,:,:) = 0.0
    where (black_cphob > 0.0)
        sinkphob = 8.038e-6*black_cphob
    elsewhere
        sinkphob = 0.0
    endwhere
!

!------- tendency ------------------

    black_cphob_dt = bc_sourcephob - sinkphob
    black_cphil_dt = bc_sourcephil + sinkphob
!
!------- compute organic carbon sink --------------
!
!  OCphob has a half-life time of 2.0days 
!   (corresponds to an e-folding time of 2.88 days)
!
!  sink = 1./(86400.*1.44) = 8.023e-6
!
    sinkphob(:,:,:) = 0.0
    where (omphob >= 0.0)
       sinkphob = 8.023e-6*omphob
    elsewhere
       sinkphob = 0.0
    endwhere

!------- tendency ------------------

      omphob_dt = om_sourcephob - sinkphob
      omphil_dt = om_sourcephil + sinkphob
!
 end subroutine atmos_carbon_aerosol_driver

!#######################################################################
!<SUBROUTINE NAME ="atmos_carbon_aerosol_init">

!<OVERVIEW>
! Subroutine to initialize the carbon aerosol module.
!</OVERVIEW>
!<DESCRIPTION>
! This subroutine querys the tracer manager to find the indices for the 
! various carbonaceous aerosol tracers. It also registers the emission 
! fields for diagnostic purposes.
!  
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_carbon_aerosol_init (lonb, latb, r, axes, Time, mask)
!</TEMPLATE>
!   <IN NAME="lonb" TYPE="real" DIM="(:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="latb" TYPE="real" DIM="(:)">
!     The latitudes for the local domain.
!   </IN>
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

 subroutine atmos_carbon_aerosol_init (lonb, latb, axes, Time, mask)

!-----------------------------------------------------------------------
real, dimension(:),    intent(in) :: lonb, latb
integer        , intent(in)                        :: axes(4)
type(time_type), intent(in)                        :: Time
real,            intent(in),    dimension(:,:,:), optional :: mask
character(len=7), parameter :: mod_name = 'tracers'
integer :: n
integer ::  unit, ierr, io


   if (module_is_initialized) return
!----------------------------------
!namelist files
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=carbon_aerosol_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'carbon_aerosol_nml')
        end do
10      call close_file (unit)
      endif
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
                          write (stdlog(), nml=carbon_aerosol_nml)

!--------------------------------------------------------
!------namelist

!----- set initial value of carbon ------------

   n = get_tracer_index(MODEL_ATMOS,'bcphob')
   if (n>0) then
      nbcphobic = n
      call set_tracer_atts(MODEL_ATMOS,'bcphob','hphobic_bc','mmr')
      if (nbcphobic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (*,30) 'Hydrophobic BC',nbcphobic
      if (nbcphobic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (stdlog(),30) 'Hydrophobic BC',nbcphobic
   endif

   n = get_tracer_index(MODEL_ATMOS,'bcphil')
   if (n>0) then
      nbcphilic=n
      call set_tracer_atts(MODEL_ATMOS,'bcphil','hphilic_bc','mmr')
      if (nbcphilic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (*,30) 'Hydrophilic BC',nbcphilic
      if (nbcphilic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (stdlog(),30) 'Hydrophilic BC',nbcphilic
   endif

   n = get_tracer_index(MODEL_ATMOS,'omphob')
   if (n>0) then
      nomphobic=n
      call set_tracer_atts(MODEL_ATMOS,'omphob','phobic_om','mmr')
      if (nomphobic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (*,30) 'Hydrophobic OC',nomphobic
      if (nomphobic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (stdlog(),30) 'Hydrophobic OC',nomphobic
   endif

   n = get_tracer_index(MODEL_ATMOS,'omphil')
   if (n>0) then
      nomphilic=n
      call set_tracer_atts(MODEL_ATMOS,'omphil','philic_om','mmr')
      if (nomphilic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (*,30) 'Hydrophilic OC',nomphilic
      if (nomphilic > 0 .and. mpp_pe() == mpp_root_pe()) &
          write (stdlog(),30) 'Hydrophilic OC',nomphilic
   endif

30        format (A,' was initialized as tracer number ',i2)
!   Register Emissions as static fields (monthly)
!
       id_bcemisff_l= register_diag_field ( mod_name,     &
                    'bcemisff_l', axes(1:2), Time,                &
                    'BC emissions from FF', 'kg/m2/sec', -999. )

       id_bcemisbb_l= register_diag_field ( mod_name,     &
                    'bcemisbb_l', axes(1:2), Time,                &
                    'BC emissions from BB', 'kg/m2/sec')

       id_omemisff_l= register_diag_field ( mod_name,     &
                    'omemisff_l', axes(1:2), Time,                &
                    'OM emissions from FF', 'kg/m2/sec')

       id_omemisbb_l= register_diag_field ( mod_name,     &
                    'omemisbb_l', axes(1:2), Time,                &
                    'OM emissions from BB', 'kg/m2/sec')
       id_bcemisff_h= register_diag_field ( mod_name,     &
                    'bcemisff_h', axes(1:2), Time,                &
                    'BC emissions from FF', 'kg/m2/sec')

       id_bcemisbb_h= register_diag_field ( mod_name,     &
                    'bcemisbb_h', axes(1:2), Time,                &
                    'BC emissions from BB', 'kg/m2/sec')

       id_omemisff_h= register_diag_field ( mod_name,     &
                    'omemisff_h', axes(1:2), Time,                &
                    'OM emissions from FF', 'kg/m2/sec')

       id_omemisbb_h= register_diag_field ( mod_name,     &
                    'omemisbb_h', axes(1:2), Time,                &
                    'OM emissions from BB', 'kg/m2/sec')

       id_omemisnat= register_diag_field ( mod_name,     &
                    'omemisnat', axes(1:2), Time,                &
                    'OM emissions from BVOCs', 'kg/m2/sec')
!
     id_bcphob_sink = register_diag_field ( mod_name,           &
                    'bcphob_sink', axes(1:3),Time,              &
                    'BC phobic loss rate', 'kg/m2/sec' )

     id_bcphil_sink = register_diag_field ( mod_name,           &
                    'bcphil_sink', axes(1:3),Time,              &
                    'BC phylic loss rate', 'kg/m2/sec' )

     id_omphob_sink = register_diag_field ( mod_name,           &
                    'omphob_sink', axes(1:3),Time,              &
                    'OM phobic loss rate', 'kg/m2/sec' )

     id_omphil_sink = register_diag_field ( mod_name,           &
                    'omphil_sink', axes(1:3),Time,              &
                    'OM phylic loss rate', 'kg/m2/sec' )
     call interpolator_init (carbon_aerosol_interp,             &
                             trim(emission_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )


   call write_version_number (version, tagname)
   module_is_initialized = .TRUE.

!-----------------------------------------------------------------------

end subroutine atmos_carbon_aerosol_init
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_carbon_aerosol_end">
!<OVERVIEW>
!  The destructor routine for the carbon module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_carbon_aero_end
!</TEMPLATE>
 subroutine atmos_carbon_aerosol_end
      call interpolator_end ( carbon_aerosol_interp)
      module_is_initialized = .FALSE.

 end subroutine atmos_carbon_aerosol_end
!</SUBROUTINE>
!#######################################################################
end module atmos_carbon_aerosol_mod
