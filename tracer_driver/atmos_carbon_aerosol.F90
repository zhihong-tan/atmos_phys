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
use atmos_tracer_utilities_mod, only : interp_emiss
use     time_manager_mod, only : time_type, get_date
use        constants_mod, only : PI,RDGAS

implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_black_carbon_driver,   &
        atmos_organic_carbon_driver,  &
        atmos_carbon_aerosol_init, &
        atmos_carbon_aerosol_end

!-----------------------------------------------------------------------
!----------- namelist -------------------
!integer      ::  bcff_source  = 1
!integer      ::  ocff_source = 2
!integer      ::  bcbb_source = 1
!integer      ::  ocbb_source = 1
!logical      ::  streets_asia = .false.
!logical      ::  rv02_india  = .false.
!biomass_only; fossil_fuels_only, natural_only, anthrop
character(len=20):: runtype = "anthrop"

namelist / aerosol_emissions_nml/  &
                                runtype
!
! tracer number for carbonaceous aerosols
integer :: nbcphobic=0
integer :: nbcphilic=0
integer :: nomphobic=0
integer :: nomphilic=0

!--- identification numbers for  diagnostic fields and axes ----
integer :: id_bcphob_sink, id_bcphil_sink
integer :: id_omphob_sink, id_omphil_sink
!
integer :: id_garea
integer :: id_bcemisff_l, id_bcemisbb_l
integer :: id_omemisff_l, id_omemisbb_l
integer :: id_bcemisff_h, id_bcemisbb_h
integer :: id_omemisff_h, id_omemisbb_h
integer :: id_omemisnat
integer :: id_faemisff_l, id_faemisff_h
integer :: id_faemisbb_l, id_faemisbb_h

! William Cooke's inventories; Cooke and Wilson, JGR, 1996; Cooke et al., 1999.
integer :: id_bcff_cw96, id_bcbb_cw96
integer :: id_bcff_cooke99, id_omff_cooke99

! Tami Bond emissions, Bond et al., JGR, 2004.
integer :: id_bcbf_bond, id_ombf_bond
integer :: id_bcff_bond, id_omff_bond
integer :: id_bcob_bond, id_omob_bond

! ATSR corrected emissions, Reddy & Boucher, JGR, 2004
integer :: id_bcbb_l_atsr, id_bcbb_h_atsr

! Asian Emissions, Streets et al., 2003;
integer :: id_bcff_l_asia, id_omff_l_asia
integer :: id_bcbb_l_asia, id_ombb_l_asia
integer :: id_bcff_h_asia, id_omff_h_asia
integer :: id_bcbb_h_asia, id_ombb_h_asia

! Reddy and Venkataraman inventories for India
!Reddy and Venkataraman, Atmos. Env.,  2002a;2002b
integer :: id_bcff_l_india, id_omff_l_india
integer :: id_bcbb_l_india, id_ombb_l_india
integer :: id_bcff_h_india, id_omff_h_india
integer :: id_bcbb_h_india, id_ombb_h_india
integer :: id_faff_l_india, id_faff_h_india
integer :: id_fabb_l_india, id_fabb_h_india

!--- Arrays to help calculate tracer sources/sinks ---
real, allocatable, dimension(:,:) :: garea
real, allocatable, dimension(:,:) :: bcemisff_l, bcemisbb_l
real, allocatable, dimension(:,:) :: omemisff_l, omemisbb_l
real, allocatable, dimension(:,:) :: bcemisff_h, bcemisbb_h
real, allocatable, dimension(:,:) :: omemisff_h, omemisbb_h
real, allocatable, dimension(:,:) :: omemisnat
real, allocatable, dimension(:,:) :: faemisff_l, faemisbb_l
real, allocatable, dimension(:,:) :: faemisff_h, faemisbb_h
!
real, allocatable, dimension(:,:) :: bcbb_l_atsr, bcbb_h_atsr

!
real, allocatable, dimension(:,:) :: bcff_cw96,bcbb_cw96
real, allocatable, dimension(:,:) :: bcff_cooke99,omff_cooke99
!
real, allocatable, dimension(:,:) :: bcbf_bond, ombf_bond
real, allocatable, dimension(:,:) :: bcff_bond, omff_bond
real, allocatable, dimension(:,:) :: bcob_bond, omob_bond
!
real, allocatable, dimension(:,:) :: bcff_l_asia, omff_l_asia
real, allocatable, dimension(:,:) :: bcff_h_asia, omff_h_asia
real, allocatable, dimension(:,:) :: bcbb_l_asia, ombb_l_asia
real, allocatable, dimension(:,:) :: bcbb_h_asia, ombb_h_asia
!
real, allocatable, dimension(:,:) :: bcff_l_india, omff_l_india
real, allocatable, dimension(:,:) :: bcff_h_india, omff_h_india
real, allocatable, dimension(:,:) :: bcbb_l_india, ombb_l_india
real, allocatable, dimension(:,:) :: bcbb_h_india, ombb_h_india
real, allocatable, dimension(:,:) :: fabb_h_india, fabb_l_india
real, allocatable, dimension(:,:) :: faff_h_india, faff_l_india

!################### CHOICE OF EMISSIONS TO THE TRANSPORT ###################
! Following choices of emission inventories are possible
! bcff_source = 1 Cooke and Wilson [1996]
!             = 2 Cooke et al. [1999]
!             = 3 Bond et al. [2004]
! ocff_source = 2 Cooke et al. [1999]
!             = 3 Bond et al. [2004]
!
! bcbb_source = 1 Cooke and Wilson [1996] with all emissions 
!                 are emitted at surface 
!             = 3 Bond et al. [2004] Open burning emissions seasonality is 
!                 inferred from Van der Werf et al. [2004]
!             = 4 Cooke and Wilson [1996] are splitted into surface/
!                 elevated sources with ATSR fire coutns
!             = 5 Van der Werf open for open burning source &
!                 Bond et al. [2004] for biofuel sources
!               
! ocbb_source = 1 Assumed OC/BC 7.0  for Cooke and Wilson [1996] 
!                 BC emissions with all emissions are emitted at surface
!             = 3 Bond et al. [2004] Open burning emissions seasonality is 
!                 inferred from Van der Werf et al. [2004]
!             = 4 Assuemed oc/bc=7 from  Cooke and Wilson [1996] and are 
!                 splitted into surface/elevated sources with ATSR fire coutns
!             = 5 Van der Werf open for open burning source [oc = bc*7.0] 
!                 Bond et al. [2004] for biofuel sources
!streets_asia= .true. | mask Asian sources with Streets et al. [2004]
!rv02_india=   .true. | mask Indian sources with Reddy and Venkataraman [2002a, b]
!##########################################################################

!Choice of emission source for bc and oc for biomass burning 
!and fossil fuels
INTEGER, PARAMETER    ::  bcff_source=3, ocff_source=3
INTEGER, PARAMETER    ::  bcbb_source=3, ocbb_source=3
LOGICAL, PARAMETER    ::  streets_asia=.false.
LOGICAL, PARAMETER    ::  rv02_india  =.false.

character(len=6), parameter :: module_name = 'tracer'

logical :: module_is_initialized = .FALSE.
logical :: used
!  local variables:

      integer        :: unit, ierr, io
      integer        :: n
      character(len=16) :: chvers

!---- version number -----
character(len=128) :: version = '$Id: atmos_carbon_aerosol.F90,v 13.0 2006/03/28 21:15:14 fms Exp $'
character(len=128) :: tagname = '$Name: memphis_2006_08 $'
!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine atmos_black_carbon_driver(lon, lat, land, pfull,phalf,T,pwt, z_half, z_pbl, &
                               black_cphob, black_cphob_dt,  &
                               black_cphil, black_cphil_dt,  &
                               Time, is, ie, js, je )

!-----------------------------------------------------------------------
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:)   :: land
   real, intent(in),  dimension(:,:)   :: z_pbl
   real, intent(in),  dimension(:,:,:) :: z_half
   real, intent(in),  dimension(:,:,:) :: pwt, black_cphob,black_cphil,pfull,phalf,T
   real, intent(out), dimension(:,:,:) :: black_cphob_dt,black_cphil_dt
type(time_type), intent(in)            :: Time
integer, intent(in)                    :: is, ie, js, je
!-----------------------------------------------------------------------
   real, dimension(size(black_cphob,1),size(black_cphob,2),size(black_cphob,3)) ::  &
         sourcephob, sinkphob, sourcephil, sinkphil
   real  dtr,bltop,z1,z2
   real, dimension(size(black_cphob,3)) :: fbb, fff
integer  i,j,l, kb,id,jd,kd,lat1,k
REAL,PARAMETER :: frac_bc_phobic = 0.8
REAL,PARAMETER :: frac_bc_philic = 0.2

!
!-------------------------------------------------------------------------
!Emissions profiles as per SAFARI extinction 
!----
   real, dimension(size(black_cphob,1),size(black_cphob,2),size(black_cphob,3)) :: dz3d,rho_air

!
!-----------------------------------------------------------------------
!
      id=size(black_cphob,1); jd=size(black_cphob,2); kd=size(black_cphob,3)

      dtr= PI/180.

!----------- compute black carbon source ------------

      sourcephob = 0.0
      sourcephil = 0.0

    rho_air(:,:,:) = pfull(:,:,:)/T(:,:,:)/RDGAS
    DO k =1, kd
      dz3d(:,:,k) = pwt(:,:,k)/rho_air(:,:,k)
    ENDDO



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
        sourcephob(i,j,kd) =  frac_bc_phobic*  &
           (bcemisff_l(i+is-1,j+js-1)+ bcemisbb_l(i+is-1,j+js-1))/pwt(i,j,kd)

        sourcephil(i,j,kd) =  frac_bc_philic*  &
           (bcemisff_l(i+is-1,j+js-1)+ bcemisbb_l(i+is-1,j+js-1))/pwt(i,j,kd)

        sourcephob(i,j,:) =  sourcephob(i,j,:) + frac_bc_phobic* &
           (fff(:)*bcemisff_h(i+is-1,j+js-1)+fbb(:)*bcemisbb_h(i+is-1,j+js-1)) &
                                                             /pwt(i,j,:)
        sourcephil(i,j,:) =  sourcephil(i,j,:) + frac_bc_philic* &
           (fff(:)*bcemisff_h(i+is-1,j+js-1)+fbb(:)*bcemisbb_h(i+is-1,j+js-1)) &
                                                             /pwt(i,j,:)
      enddo
      enddo

!------- compute black carbon phobic sink --------------
!
!  BCphob has a half-life time of 1.0days 
!   (corresponds to an e-folding time of 1.44 days)
!
!  sink = 1./(86400.*1.44) = 8.023e-6

      sinkphil(:,:,:) = 0.0
      sinkphob(:,:,:) = 0.0
      where (black_cphob > 0.0)
        sinkphob = 8.038e-6*black_cphob
      elsewhere
        sinkphob = 0.0
      endwhere
!

!------- tendency ------------------

      black_cphob_dt=sourcephob-sinkphob
      black_cphil_dt=sourcephil+sinkphob
!
 end subroutine atmos_black_carbon_driver

!#######################################################################

 subroutine atmos_organic_carbon_driver(lon, lat, land, pfull,phalf,T,pwt,z_half, z_pbl, &
                               omphob, omphob_dt,  &
                               omphil, omphil_dt,  &
                               Time, is, ie, js, je)

!-----------------------------------------------------------------------
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:)   :: land
   real, intent(in),  dimension(:,:)   :: z_pbl
   real, intent(in),  dimension(:,:,:) :: z_half
   real, intent(in),  dimension(:,:,:) :: pwt, omphob,omphil,pfull,phalf,T
   real, intent(out), dimension(:,:,:) :: omphob_dt,omphil_dt
type(time_type), intent(in)            :: Time
integer, intent(in)                    :: is, ie, js, je
!-----------------------------------------------------------------------
   real, dimension(size(omphob,1),size(omphob,2),size(omphob,3)) ::  &
         sourcephob, sinkphob, sourcephil, sinkphil
   real, dimension(size(omphob,3)) :: fbb, fff, fna
   real  dtr, bltop, z1, z2
integer  i,j,l, kb,id,jd,kd,lat1,k
REAL,PARAMETER :: frac_om_phobic = 0.5
REAL,PARAMETER :: frac_om_philic = 0.5
!
!-----------------------------------------------------------------------
!----
   real, dimension(size(omphob,1),size(omphob,2),size(omphob,3)) :: dz3d,rho_air

      id=size(omphob,1); jd=size(omphob,2); kd=size(omphob,3)

      dtr= PI/180.


    rho_air(:,:,:) = pfull(:,:,:)/T(:,:,:)/RDGAS
    DO k =1, kd
      dz3d(:,:,k) = pwt(:,:,k)/rho_air(:,:,k)
    ENDDO

!---------------------------------------------- 
      sourcephob = 0.0
      sourcephil = 0.0

!----------- compute organic carbon source ------------
! --- Assuming biomass burning emission within the PBL -------

      do j = 1, jd
      do i = 1, id

        fbb(:) = 0.
        fbb(kd)= 1.
        fff(:)= 0; 
        fff(kd-1)= 1.
        fna(:)= 0; 
        fna(kd)= 1.
!
        BLTOP = z_pbl(i,j)
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
        sourcephob(i,j,kd) =  frac_om_phobic*  &
           (omemisff_l(i+is-1,j+js-1)+ omemisbb_l(i+is-1,j+js-1) &
            +omemisnat(i+is-1,j+js-1))/pwt(i,j,kd)

        sourcephil(i,j,kd) =  frac_om_philic*  &
           (omemisff_l(i+is-1,j+js-1)+ omemisbb_l(i+is-1,j+js-1) &
            +omemisnat(i+is-1,j+js-1))/pwt(i,j,kd)

        sourcephob(i,j,:) =  sourcephob(i,j,:) + frac_om_phobic* &
           (fff(:)*omemisff_h(i+is-1,j+js-1)+fbb(:)*omemisbb_h(i+is-1,j+js-1)) &
                                                             /pwt(i,j,:)
        sourcephil(i,j,:) =  sourcephil(i,j,:) + frac_om_philic* &
           (fff(:)*omemisff_h(i+is-1,j+js-1)+fbb(:)*omemisbb_h(i+is-1,j+js-1)) &
                                                             /pwt(i,j,:)
      enddo
      enddo
!
!------- compute organic carbon sink --------------
!
!  OCphob has a half-life time of 2.0days 
!   (corresponds to an e-folding time of 2.88 days)
!
!  sink = 1./(86400.*1.44) = 8.023e-6
!
    sinkphob(:,:,:) = 0.0
    sinkphil(:,:,:) = 0.0
    where (omphob >= 0.0)
       sinkphob = 8.023e-6*omphob
    elsewhere
       sinkphob = 0.0
    endwhere

!------- tendency ------------------

      omphob_dt=sourcephob-sinkphob
      omphil_dt=sourcephil+sinkphob
      

 end subroutine atmos_organic_carbon_driver


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


   if (module_is_initialized) return
!----------------------------------
!namelist files
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

!--------------------------------------------------------
!------namelist

!----- set initial value of carbon ------------

   n = get_tracer_index(MODEL_ATMOS,'bcphob')
   if (n>0) then
      nbcphobic = n
      call set_tracer_atts(MODEL_ATMOS,'bcphob','hphobic_bc','mmr')
      if (nbcphobic > 0 .and. mpp_pe() == mpp_root_pe()) write (*,30) 'Hydrophobic BC',nbcphobic
      if (nbcphobic > 0 .and. mpp_pe() == mpp_root_pe()) write (stdlog(),30) 'Hydrophobic BC',nbcphobic
   endif

   n = get_tracer_index(MODEL_ATMOS,'bcphil')
   if (n>0) then
      nbcphilic=n
      call set_tracer_atts(MODEL_ATMOS,'bcphil','hphilic_bc','mmr')
      if (nbcphilic > 0 .and. mpp_pe() == mpp_root_pe()) write (*,30) 'Hydrophilic BC',nbcphilic
      if (nbcphilic > 0 .and. mpp_pe() == mpp_root_pe()) write (stdlog(),30) 'Hydrophilic BC',nbcphilic
   endif

   n = get_tracer_index(MODEL_ATMOS,'omphob')
   if (n>0) then
      nomphobic=n
      call set_tracer_atts(MODEL_ATMOS,'omphob','phobic_om','mmr')
      if (nomphobic > 0 .and. mpp_pe() == mpp_root_pe()) write (*,30) 'Hydrophobic OC',nomphobic
      if (nomphobic > 0 .and. mpp_pe() == mpp_root_pe()) write (stdlog(),30) 'Hydrophobic OC',nomphobic
   endif

   n = get_tracer_index(MODEL_ATMOS,'omphil')
   if (n>0) then
      nomphilic=n
      call set_tracer_atts(MODEL_ATMOS,'omphil','philic_om','mmr')
      if (nomphilic > 0 .and. mpp_pe() == mpp_root_pe()) write (*,30) 'Hydrophilic OC',nomphilic
      if (nomphilic > 0 .and. mpp_pe() == mpp_root_pe()) write (stdlog(),30) 'Hydrophilic OC',nomphilic
   endif

30        format (A,' was initialized as tracer number ',i2)
      !Read in emission files
!####################################################################
!####################Register Emissions as static fields (monthly)
!
     if (.not.allocated(garea)) then
       allocate (garea(size(lonb)-1,size(latb)-1))
       id_garea = register_static_field ( mod_name,     &
                    'garea', axes(1:2), 'Grid area', 'm2' )
     endif
!
       allocate (bcemisff_l(size(lonb)-1,size(latb)-1))
       allocate (bcemisbb_l(size(lonb)-1,size(latb)-1))
       allocate (omemisff_l(size(lonb)-1,size(latb)-1))
       allocate (omemisbb_l(size(lonb)-1,size(latb)-1))
       allocate (bcemisff_h(size(lonb)-1,size(latb)-1))
       allocate (bcemisbb_h(size(lonb)-1,size(latb)-1))
       allocate (omemisff_h(size(lonb)-1,size(latb)-1))
       allocate (omemisbb_h(size(lonb)-1,size(latb)-1))
       allocate (omemisnat(size(lonb)-1,size(latb)-1))
!
       allocate (faemisff_l(size(lonb)-1,size(latb)-1))
       allocate (faemisbb_l(size(lonb)-1,size(latb)-1))
       allocate (faemisff_h(size(lonb)-1,size(latb)-1))
       allocate (faemisbb_h(size(lonb)-1,size(latb)-1))

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

       id_faemisff_l= register_diag_field ( mod_name,     &
                    'faemisff_l', axes(1:2), Time,                &
		    'FA emissions from FF', 'kg/m2/sec')

       id_faemisbb_l= register_diag_field ( mod_name,     &
                    'faemisbb_l', axes(1:2), Time,                &
		    'FA emissions from BB', 'kg/m2/sec')

       id_faemisff_h= register_diag_field ( mod_name,     &
                    'faemisff_h', axes(1:2), Time,                &
		    'FA emissions from FF', 'kg/m2/sec')

       id_faemisbb_h= register_diag_field ( mod_name,     &
                    'faemisbb_h', axes(1:2), Time,                &
		    'FA emissions from BB', 'kg/m2/sec')

! 
! ATSR corrected emissions
       allocate (bcbb_l_atsr(size(lonb)-1,size(latb)-1))
       allocate (bcbb_h_atsr(size(lonb)-1,size(latb)-1))

       id_bcbb_l_atsr= register_diag_field ( mod_name,     &
                    'bcbb_l_atsr', axes(1:2), Time,                &
		    'BC emissions from CW96 -low', 'kg/m2/sec')

       id_bcbb_h_atsr= register_diag_field ( mod_name,     &
                    'bcbb_h_atsr', axes(1:2), Time,                &
		    'BC emissions from CW96 -low', 'kg/m2/sec')
!
! Cooke's Inventories
!
       allocate (bcff_cw96(size(lonb)-1,size(latb)-1))
       allocate (bcbb_cw96(size(lonb)-1,size(latb)-1))
       allocate (bcff_cooke99(size(lonb)-1,size(latb)-1))
       allocate (omff_cooke99(size(lonb)-1,size(latb)-1))

       id_bcff_cw96= register_diag_field ( mod_name,     &
                    'bcff_cw96', axes(1:2), Time,                &
		    'BC emissions from FF CW96', 'kg/m2/sec')
       id_bcbb_cw96= register_diag_field ( mod_name,     &
                    'bcbb_cw96', axes(1:2), Time,                &
		    'BC emissions from BB W96', 'kg/m2/sec')
       id_bcff_cooke99= register_diag_field ( mod_name,     &
                    'bcff_cooke99', axes(1:2), Time,                &
		    'BC emissions from FF Cooke et al.,1999 ', 'kg/m2/sec')
       id_omff_cooke99= register_diag_field ( mod_name,     &
                    'omff_cooke99', axes(1:2), Time,                &
		    'OM emissions from FF Cooke et al.,1999 ', 'kg/m2/sec')
!
! Tami Bond emissions
!
       allocate (bcbf_bond(size(lonb)-1,size(latb)-1))
       allocate (ombf_bond(size(lonb)-1,size(latb)-1))
       allocate (bcff_bond(size(lonb)-1,size(latb)-1))
       allocate (omff_bond(size(lonb)-1,size(latb)-1))
       allocate (bcob_bond(size(lonb)-1,size(latb)-1))
       allocate (omob_bond(size(lonb)-1,size(latb)-1))
!
       id_bcbf_bond= register_diag_field ( mod_name,     &
                    'bcbf_bond', axes(1:2), Time,                &
		    'BC emissions from Biofuels', 'kg/m2/sec')

       id_ombf_bond= register_diag_field ( mod_name,     &
                    'ombf_bond', axes(1:2), Time,                &
		    'OM emissions from Biofuels', 'kg/m2/sec')

       id_bcff_bond= register_diag_field ( mod_name,     &
                    'bcff_bond', axes(1:2), Time,                &
		    'BC emissions from FF', 'kg/m2/sec')

       id_omff_bond= register_diag_field ( mod_name,     &
                    'omff_bond', axes(1:2), Time,                &
		    'OM emissions from FF', 'kg/m2/sec')

       id_bcob_bond= register_diag_field ( mod_name,     &
                    'bcob_bond', axes(1:2), Time,                &
		    'BC emissions from OB', 'kg/m2/sec')

       id_omob_bond= register_diag_field ( mod_name,     &
                    'omob_bond', axes(1:2), Time,                &
		    'OM emissions from OB', 'kg/m2/sec')
!
       allocate (bcbb_l_asia(size(lonb)-1,size(latb)-1))
       allocate (bcbb_h_asia(size(lonb)-1,size(latb)-1))
       allocate (ombb_l_asia(size(lonb)-1,size(latb)-1))
       allocate (ombb_h_asia(size(lonb)-1,size(latb)-1))
       allocate (bcff_l_asia(size(lonb)-1,size(latb)-1))
       allocate (bcff_h_asia(size(lonb)-1,size(latb)-1))
       allocate (omff_l_asia(size(lonb)-1,size(latb)-1))
       allocate (omff_h_asia(size(lonb)-1,size(latb)-1))
!
       allocate (bcbb_l_india(size(lonb)-1,size(latb)-1))
       allocate (bcbb_h_india(size(lonb)-1,size(latb)-1))
       allocate (ombb_l_india(size(lonb)-1,size(latb)-1))
       allocate (ombb_h_india(size(lonb)-1,size(latb)-1))
       allocate (bcff_l_india(size(lonb)-1,size(latb)-1))
       allocate (bcff_h_india(size(lonb)-1,size(latb)-1))
       allocate (omff_l_india(size(lonb)-1,size(latb)-1))
       allocate (omff_h_india(size(lonb)-1,size(latb)-1))
       allocate (faff_h_india(size(lonb)-1,size(latb)-1))
       allocate (faff_l_india(size(lonb)-1,size(latb)-1))
       allocate (fabb_h_india(size(lonb)-1,size(latb)-1))
       allocate (fabb_l_india(size(lonb)-1,size(latb)-1))

       id_bcbb_l_asia= register_diag_field ( mod_name,     &
                    'bcbb_l_asia', axes(1:2), Time,                &
		    'BC emissions from low/BB ', 'kg/m2/sec')
       id_bcbb_h_asia= register_diag_field ( mod_name,     &
                    'bcbb_h_asia', axes(1:2), Time,                &
		    'BC emissions from high/BB ', 'kg/m2/sec')
       id_ombb_l_asia= register_diag_field ( mod_name,     &
                    'ombb_l_asia', axes(1:2), Time,                &
		    'BC emissions from low/BB ', 'kg/m2/sec')
       id_ombb_h_asia= register_diag_field ( mod_name,     &
                    'ombb_h_asia', axes(1:2), Time,                &
		    'BC emissions from high/BB ', 'kg/m2/sec')
       id_bcff_l_asia= register_diag_field ( mod_name,     &
                    'bcff_l_asia', axes(1:2), Time,                &
		    'BC emissions from low/FF', 'kg/m2/sec')
       id_bcff_h_asia= register_diag_field ( mod_name,     &
                    'bcff_h_asia', axes(1:2), Time,                &
		    'BC emissions from high/FF', 'kg/m2/sec')
       id_omff_l_asia= register_diag_field ( mod_name,     &
                    'omff_l_asia', axes(1:2), Time,                &
		    'OM emissions from low/FF', 'kg/m2/sec')
       id_omff_h_asia= register_diag_field ( mod_name,     &
                    'omff_h_asia', axes(1:2), Time,                &
		    'OM emissions from high/FF', 'kg/m2/sec')
!
       id_bcbb_l_india= register_diag_field ( mod_name,     &
                    'bcbb_l_india', axes(1:2), Time,                &
		    'BC emissions from low/BB ', 'kg/m2/sec')
       id_bcbb_h_india= register_diag_field ( mod_name,     &
                    'bcbb_h_india', axes(1:2), Time,                &
		    'BC emissions from high/BB ', 'kg/m2/sec')
       id_ombb_l_india= register_diag_field ( mod_name,     &
                    'ombb_l_india', axes(1:2), Time,                &
		    'BC emissions from low/BB ', 'kg/m2/sec')
       id_ombb_h_india= register_diag_field ( mod_name,     &
                    'ombb_h_india', axes(1:2), Time,                &
		    'BC emissions from high/BB ', 'kg/m2/sec')
       id_bcff_l_india= register_diag_field ( mod_name,     &
                    'bcff_l_india', axes(1:2), Time,                &
		    'BC emissions from low/FF', 'kg/m2/sec')
       id_bcff_h_india= register_diag_field ( mod_name,     &
                    'bcff_h_india', axes(1:2), Time,                &
		    'BC emissions from high/FF', 'kg/m2/sec')
       id_omff_l_india= register_diag_field ( mod_name,     &
                    'omff_l_india', axes(1:2), Time,                &
		    'OM emissions from low/FF', 'kg/m2/sec')
       id_omff_h_india= register_diag_field ( mod_name,     &
                    'omff_h_india', axes(1:2), Time,                &
!
		    'OM emissions from high/FF', 'kg/m2/sec')
       id_fabb_l_india= register_diag_field ( mod_name,     &
                    'fabb_l_india', axes(1:2), Time,                &
		    'FA emissions from low/BB ', 'kg/m2/sec')
       id_fabb_h_india= register_diag_field ( mod_name,     &
                    'fabb_h_india', axes(1:2), Time,                &
		    'FA emissions from high/BB ', 'kg/m2/sec')
       id_faff_l_india= register_diag_field ( mod_name,     &
                    'faff_l_india', axes(1:2), Time,                &
		    'FA emissions from low/FF', 'kg/m2/sec')
       id_faff_h_india= register_diag_field ( mod_name,     &
                    'faff_h_india', axes(1:2), Time,                &
		    'FA emissions from high/FF', 'kg/m2/sec')
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

   call carbon_aerosol_source_input(lonb, latb, Time)

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
 
      module_is_initialized = .FALSE.

 end subroutine atmos_carbon_aerosol_end
!</SUBROUTINE>
!#######################################################################
 subroutine carbon_aerosol_source_input(lonb, latb, Time)
real, dimension(:),    intent(in) :: lonb, latb
type(time_type),intent(in) :: Time

integer      :: i, j, unit, io, ios, itime,ii,nn
INTEGER, PARAMETER        :: imax = 360, jmax=180
INTEGER, PARAMETER        :: imax_india = 69, jmax_india=73
real         :: dtr, lat_S, lon_W, dlon, dlat,sumbc,bcsum,ocsum
real         :: lat_S_india, lon_W_india, dlon_india, dlat_india
real,DIMENSION(imax, jmax)    :: emiss1, emiss11,grid_area
real,DIMENSION(imax, jmax)    :: bcl,bch,ocl,och,fal,fah
real,DIMENSION(imax, jmax,12) :: emiss2
real,DIMENSION(12) :: rdata(12)

REAL, PARAMETER           :: sec_in_day = 86400., day_in_year = 365.
REAL, PARAMETER           :: r_globe = 6371229.
REAL, PARAMETER           :: om_oc_ff= 1.4, om_oc_bb = 1.6,om_oc_nat = 1.4
CHARACTER(len=100)        :: header,junk
INTEGER                   :: second, day,minute,year,month,hour
LOGICAL                   :: opened
LOGICAL, PARAMETER        :: tilleof = .false.
REAL, DIMENSION(12)       :: day_in_month 
real :: areasum
! TAMI BOND emissions
REAL                      :: bc,oc,so
REAL                      :: rlon,rlat,remiss
real,DIMENSION(imax, jmax)    :: bcbond,ocbond,bc_asia,oc_asia 
CHARACTER(len=20),DIMENSION(9):: rdummy

!----------------------------------------------------------------------
!  local variables:

      integer   ::   ierr, n
!---------------------------------------------------------------------
!    local variables:
!
!         unit       io unit number used to read namelist file
!         ierr       error code
!         io         error status returned from io operation
!         n          do-loop index
!
!-----------------------------------------------------------------------
    DATA day_in_month/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./
!
!-----------some definations----------------------
    dtr= PI/180.
    lon_W = 0.0*dtr; lat_S = -90.0*dtr   ! west & south corner grid
    dlon  = 1.0*dtr; dlat  = 1.0*dtr     ! lon & lat ressolution
!
    lon_W_india = 248.125*dtr;lat_S_india=6.125*dtr   ! west & south corner grid
    dlon_india  = 0.25*dtr;  dlat_india  = 0.25*dtr     ! lon & lat ressolution

!-----Calculating grid-box area (1x1)-----
!
    DO j=1, jmax
    DO i=1, imax
    grid_area(i,j) = r_globe*r_globe*2.*PI/360.* &          !m2
                     (SIN((j-90)*2.*PI/360.)-SIN((j-91)*2.*PI/360.))
    ENDDO
    ENDDO
!
    call interp_emiss ( grid_area, lon_W, lat_S, dlon, dlat, garea)
! Send the source data to the diag_manager for output.
         if (id_garea> 0 )  &
           used = send_data ( id_garea, garea, Time )
!
!  getting month and day of the year for emssions

call    get_date (Time, year, month, day, hour, minute, second)
!
!################################# COOKE's Inventories
! Emissions sare given as tonnes per grid box (1x1)
! Reading Cooke et al. [1999] emissions : Fossil fuels only!
!
! Reading Cooke and Wilson [1996]/GEIA emissions 
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         emiss1= 0.0
         open (unit,file='INPUT/BC_annual.1x1', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
            DO WHILE (.not.tilleof)
                read  (unit, FMT=2001, end=21) j,i,emiss1(i,j)
            ENDDO
         else
             write(*,*) '***ERROR: Opening BC FF emissions file'
         endif
21      call close_file (unit)
2001    FORMAT(2(I3),F11.3)

        emiss1 = emiss1*1.E3/grid_area/day_in_year/sec_in_day  !kg BC/m2/sec
        call interp_emiss ( emiss1, lon_W, lat_S, dlon, dlat, &
                     bcff_cooke99)
         if (id_bcff_cooke99 > 0 ) &
           used = send_data ( id_bcff_cooke99, bcff_cooke99, Time )

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading bcff_cooke99'
endif

! Reading OC 
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         emiss1= 0.0
         open (unit,file='INPUT/OC_annual.1x1', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
           DO WHILE (.not.tilleof)
           read  (unit, FMT=2001, end=22) j,i,emiss1(i,j)
           ENDDO
22    call close_file (unit)
         else
            write(*,*) '***ERROR: Opening OC FF emissions file'
         endif
         emiss1 = om_oc_ff*emiss1*1.E3/grid_area/day_in_year/sec_in_day  !kg OC/m2/sec
!
!NOTE: Organic Matter (OM) emissions are dervied assumng an OM to OC 
!      ratio of 1.4 for fossil fuel combustion 
!      Reddy and Boucher [2004] and Turpin and Lin [2002].
!
        call interp_emiss ( emiss1, lon_W, lat_S, dlon, dlat, &
                     omff_cooke99)
         if (id_omff_cooke99 > 0 ) &
           used = send_data ( id_omff_cooke99, omff_cooke99, Time )
!
if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading omff_cooke99'
endif
!
!#############################################################
!  Reading BC emissions for biomass burning [Cooke and Wilson, 1996]
!  Seasonsal emissions are given as tonnes per grid box per month
!
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         emiss1= 0.0
         emiss2= 0.0
         open (unit,file='INPUT/biobc87mn1.1a', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
         DO WHILE (.not.tilleof)
         read  (unit, FMT=2003, end=23) j,i, &
                (emiss2(i,j,itime),itime=1,12)
         ENDDO
         else
         write(*,*) '***ERROR: Opening BC BB emissions file'
         endif
23      call close_file (unit)
2003    FORMAT(2(I3),12(E11.4))

! assigning data corresponding to simmulation month
        emiss1(:,:) =  emiss2(:,:,month)         
        emiss1 = emiss1*1.E3/grid_area/day_in_month(month)/sec_in_day  !kg BC/m2/sec
        call interp_emiss ( emiss1, lon_W, lat_S, dlon, dlat, &
                     bcbb_cw96)
         if (id_bcbb_cw96 > 0 ) &
           used = send_data ( id_bcbb_cw96, bcbb_cw96, Time )
if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading bcbb_cw96'
endif
!
! BC  Fossil fuels

do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         emiss1= 0.0
         open (unit,file='INPUT/antbc84yr1.1a', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
         DO WHILE (.not.tilleof)
         read  (unit, FMT=2004, end=24) j,i, emiss1(i,j)
         ENDDO
         else
         write(*,*) '***ERROR: Opening BC BB emissions file'
         endif
24      call close_file (unit)
2004    FORMAT(2(I3),E11.4)

! assigning data corresponding to simmulation month
        emiss1(:,:) = emiss1(:,:)*1.E3/grid_area(:,:) &
	              /day_in_year/sec_in_day  !kg BC/m2/sec
        call interp_emiss ( emiss1, lon_W, lat_S, dlon, dlat, &
                     bcff_cw96)
         if (id_bcff_cw96 > 0 ) &
           used = send_data ( id_bcff_cw96, bcff_cw96, Time )
!
if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading bcff_cw96'
endif

!#############################################################
!  Reading BC emissions for biomass burning which are scaled to 
!   ATSR fire coutns [Reddy and Boucher, 2004]
!  Seasonsal emissions are given as g/cm2/sec
!
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo

         emiss1= 0.0
         emiss2= 0.0
         open (unit,file='INPUT/bc_geia_low_atsr.1x1', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
         DO WHILE (.not.tilleof)
         read  (unit, FMT=2005, end=25) j,i,(emiss2(i,j,itime),itime=1,12)
         ENDDO
         else
         write(*,*) '***ERROR: Opening BC BB emissions file'
         endif
25      call close_file (unit)
2005    FORMAT(2(I3),12(E12.3))
!
! assigning data corresponding to simmulation month
        emiss1(:,:)= emiss2(:,:,month)*1.e3/sec_in_day/day_in_month(month) &
	              /grid_area(:,:)       

        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat, &
                     bcbb_l_atsr)

         if (id_bcbb_l_atsr> 0 ) &
           used = send_data ( id_bcbb_l_atsr, bcbb_l_atsr, Time )

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading bcbb_l_atsr'
endif


!###########################################
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         emiss1= 0.0
         emiss2= 0.0
         open (unit,file='INPUT/bc_geia_high_atsr.1x1', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
         DO WHILE (.not.tilleof)
         read  (unit, FMT=2006, end=26) j,i,(emiss2(i,j,itime),itime=1,12)
         ENDDO
         else
         write(*,*) '***ERROR: Opening BC BB emissions file'
         endif
26      call close_file (unit)

! assigning data corresponding to simmulation month
!
        emiss1(:,:) =  emiss2(:,:,month)*1.e3/sec_in_day/day_in_month(month) &
	              /grid_area(:,:)       

        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat, &
                     bcbb_h_atsr)
         if (id_bcbb_h_atsr> 0 ) &
           used = send_data ( id_bcbb_h_atsr, bcbb_h_atsr, Time )

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading bcbb_h_atsr'
endif
!
!################################################################
! TAMI BOND EMISSION INVENTORIES
! Fossil fuels

do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         bcbond= 0.0
         ocbond= 0.0
         open (unit,file='INPUT/ff_bond.1x1', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
	   DO WHILE (.not.tilleof)
             READ(unit,FMT=111,end=11)j,i,bcbond(i,j), &
	                                  ocbond(i,j)
          ENDDO
         else
            write(*,*) '***ERROR: Opening fossil fuels emissions file: ff_bond.1x1'
         endif
11      call close_file (unit)
111   FORMAT(2(I3),2(E12.3))

        bcbond= bcbond*1.E3/grid_area/day_in_year/sec_in_day      !kg BC/m2/sec
        ocbond= om_oc_ff*ocbond*1.E3/grid_area/day_in_year/sec_in_day  !kg OM/m2/sec

        call interp_emiss (bcbond, lon_W, lat_S, dlon, dlat, &
                     bcff_bond)
        call interp_emiss (ocbond, lon_W, lat_S, dlon, dlat, &
                     omff_bond)

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading Fossil Fuel source = BOND'
endif

! BIOFUELS EMISSIONS
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         bcbond= 0.0
         ocbond= 0.0
         open (unit,file='INPUT/bf_bond.1x1', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
	   DO WHILE (.not.tilleof)
             READ(unit,FMT=111,end=12)j,i,bcbond(i,j),ocbond(i,j)
          ENDDO
         else
            write(*,*) '***ERROR: Opening biofuels emissions file: bf_bond.1x1'
         endif
12      call close_file (unit)

        bcbond= bcbond*1.E3/grid_area/day_in_year/sec_in_day  !kg BC/m2/sec
        ocbond= om_oc_bb*ocbond*1.E3/grid_area/day_in_year/sec_in_day  !kg BC/m2/sec

        call interp_emiss (bcbond, lon_W, lat_S, dlon, dlat, &
                     bcbf_bond)
        call interp_emiss (ocbond, lon_W, lat_S, dlon, dlat, &
                     ombf_bond)

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading Biofuels source = BOND'
endif
!
!-- Tami BOND Emissions  with seasonality inferred from ATSR fire counts
!
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         emiss1= 0.0
         emiss2= 0.0
         open (unit,file='INPUT/BC_BOND_m.txt', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
	  DO WHILE (.not.tilleof)
             READ(unit,FMT=2025,end=16) & 
	                   j,i,(emiss2(i,j,itime),itime=1,12) !kt/grid/month
          ENDDO
         else
             write(*,*) '***ERROR: Opening BC_OB_BOND emissions file'
         endif
16      call close_file (unit)
2025    FORMAT(2(I3),11(E12.4,','),',',E12.4)

        emiss1(:,:) = emiss2(:,:,month)*1.e6/grid_area/sec_in_day/day_in_month(month) !kg /m2/sec
        call interp_emiss ( emiss1, lon_W, lat_S, dlon, dlat, &
                     bcob_bond)
        

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading BC open burning source = BOND'
endif
!
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         emiss1= 0.0
         emiss2= 0.0
         open (unit,file='INPUT/OC_BOND_m.txt', &
               form='formatted', action='read', iostat=ios )

         if (ios.eq.0 ) then
	  DO WHILE (.not.tilleof)
             READ(unit,FMT=2025,end=17) & 
	                   j,i,(emiss2(i,j,itime),itime=1,12) !kt/grid/month
          ENDDO
         else
             write(*,*) '***ERROR: Opening OC_OB_BOND emissions file'
         endif
17      call close_file (unit)

        emiss1(:,:) = om_oc_bb*emiss2(:,:,month)*1.e6/grid_area/sec_in_day/day_in_month(month) !kg /m2/sec
        call interp_emiss ( emiss1, lon_W, lat_S, dlon, dlat, &
                     omob_bond)

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading OC open burning source = BOND'
endif
!
! write to output
!
         if (id_bcbf_bond > 0 ) &
           used = send_data ( id_bcbf_bond, bcbf_bond, Time )
         if (id_ombf_bond > 0 ) &
           used = send_data ( id_ombf_bond, ombf_bond, Time )
         if (id_bcff_bond > 0 ) &
           used = send_data ( id_bcff_bond, bcff_bond, Time )
         if (id_omff_bond > 0 ) &
           used = send_data ( id_omff_bond, omff_bond, Time )
         if (id_bcob_bond > 0 ) &
           used = send_data ( id_bcob_bond, bcob_bond, Time )
         if (id_omob_bond > 0 ) &
           used = send_data ( id_omob_bond, omob_bond, Time )
!
!#############################################################
!
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         emiss1= 0.0
         emiss2= 0.0
         open (unit,file='INPUT/ocob_atsr_asia.1x1', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
         DO WHILE (.not.tilleof)
         read  (unit, FMT=2006, end=51) j,i,(emiss2(i,j,itime),itime=1,12)
         ENDDO
         else
         write(*,*) '***ERROR: Opening BC BB emissions file'
         endif
2006    FORMAT(2(I3),12(E12.3))
51      call close_file (unit)

!
! assigning data corresponding to simmulation month
!
        emiss1(:,:)= om_oc_bb*emiss2(:,:,month)*1.e3/sec_in_day/day_in_month(month) &
	              /grid_area(:,:)       

        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat, &
                     ombb_h_asia)

         if (id_ombb_h_asia> 0 ) &
           used = send_data ( id_ombb_h_asia, ombb_h_asia, Time )

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading bcob_h_asia'
endif

!
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         emiss1= 0.0
         emiss2= 0.0
         open (unit,file='INPUT/bcob_atsr_asia.1x1', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
         DO WHILE (.not.tilleof)
         read  (unit, FMT=2006, end=45) j,i,(emiss2(i,j,itime),itime=1,12)
         ENDDO
         else
         write(*,*) '***ERROR: Opening BC BB emissions file'
         endif
45      call close_file (unit)

        emiss1(:,:)= emiss2(:,:,month)*1.e3/sec_in_day/day_in_month(month) &
	              /grid_area(:,:)       
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat, &
                     bcbb_h_asia)
         if (id_bcbb_h_asia > 0 ) &
           used = send_data ( id_bcbb_h_asia, bcbb_h_asia, Time )

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading bcob_h_asia'
endif

!###################################
         bc_asia= 0.0
         oc_asia= 0.0
         open (unit=30,file='INPUT/bc_lps_asia_2000.1x1', &
               form='formatted', action='read')
         open (unit=40,file='INPUT/oc_lps_asia_2000.1x1', &
               form='formatted', action='read')

         READ (30,*,end=53) junk, junk,junk
         DO WHILE (.not.tilleof)
           READ(30,*,end=53)rlon,rlat,remiss
           i=181+INT(rlon)
	   IF(rlat.GT.0.0)THEN 
	     j=91+INT(rlat)  
           ELSE
	     j=91-INT(-1.0*rlat)
           ENDIF
	   bc_asia(i,j) = remiss
         ENDDO
53      call close_file (30)

         READ (40,*,end=54) junk, junk,junk
         DO WHILE (.not.tilleof)
           READ(40,*,end=54)rlon,rlat,remiss
           i=181+INT(rlon)
	   IF(rlat.GT.0.0)THEN 
	     j=91+INT(rlat)  
           ELSE
	     j=91-INT(-1.0*rlat)
           ENDIF
	   oc_asia(i,j) = remiss
          ENDDO
54      call close_file (40)
!
        emiss1(:,:)= 1.2611*bc_asia(:,:)*1.e3/sec_in_day/day_in_year &
	              /grid_area(:,:)       
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat, &
                     bcff_h_asia)
!
        emiss1(:,:)= 1.2611*om_oc_ff*oc_asia(:,:)*1.e3/sec_in_day/day_in_year &
	              /grid_area(:,:)       
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat, &
                     omff_h_asia)

         if (id_bcff_h_asia > 0 ) &
           used = send_data ( id_bcff_h_asia, bcff_h_asia, Time )

         if (id_omff_h_asia > 0 ) &
           used = send_data ( id_omff_h_asia, omff_h_asia, Time )

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading Emissions: Asia for HIGH FF'
endif
!
         bc_asia= 0.0
         oc_asia= 0.0
         open (unit=30,file='INPUT/bc_are_asia_2000.1x1', &
               form='formatted', action='read')
         open (unit=40,file='INPUT/oc_are_asia_2000.1x1', &
               form='formatted', action='read')

         READ (30,*,end=55) junk, junk,junk

         DO WHILE (.not.tilleof)
           READ(30,*,end=55)rlon,rlat,remiss
           i=181+INT(rlon)
	   IF(rlat.GT.0.0)THEN 
	     j=91+INT(rlat)  
           ELSE
	     j=91-INT(-1.0*rlat)
           ENDIF
	   bc_asia(i,j) = remiss
         ENDDO
55      call close_file (30)

         READ (40,*,end=56) junk, junk,junk
         DO WHILE (.not.tilleof)
           READ(40,*,end=56)rlon,rlat,remiss
           i=181+INT(rlon)
	   IF(rlat.GT.0.0)THEN 
	     j=91+INT(rlat)  
           ELSE
	     j=91-INT(-1.0*rlat)
           ENDIF
	   oc_asia(i,j) = remiss
         ENDDO
56      call close_file (40)
!
        bc_asia = bc_asia*1.e3/sec_in_day/day_in_year/grid_area
        oc_asia = om_oc_ff*oc_asia*1.e3/sec_in_day/day_in_year/grid_area

        call interp_emiss (bc_asia, lon_W, lat_S, dlon, dlat, &
                     bcff_l_asia)
        call interp_emiss (oc_asia, lon_W, lat_S, dlon, dlat, &
                     omff_l_asia)

         if (id_bcff_l_asia > 0 ) &
           used = send_data ( id_bcff_l_asia, bcff_l_asia, Time )
         if (id_omff_l_asia > 0 ) &
           used = send_data ( id_omff_l_asia, omff_l_asia, Time )

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading Emissions: Asia for Low FF'
endif
!#
!#
         bc_asia= 0.0
         oc_asia= 0.0
         open (unit=30,file='INPUT/bb_are_asia_2000.1x1', &
               form='formatted', action='read')

         READ (30,*,end=57) junk, junk,junk,junk,junk

         DO WHILE (.not.tilleof)
           READ(30,*,end=57)rlon,rlat,so,bc,oc
           i=181+INT(rlon)
	   IF(rlat.GT.0.0)THEN 
	     j=91+INT(rlat)  
           ELSE
	     j=91-INT(-1.0*rlat)
           ENDIF
	   bc_asia(i,j) = bc
	   oc_asia(i,j) = oc
         ENDDO
57      call close_file (30)

        bc_asia = bc_asia*1.e3/sec_in_day/day_in_year/grid_area
        oc_asia = om_oc_bb*oc_asia*1.e3/sec_in_day/day_in_year/grid_area

        call interp_emiss (bc_asia, lon_W, lat_S, dlon, dlat, &
                     bcbb_l_asia)
        call interp_emiss (oc_asia, lon_W, lat_S, dlon, dlat, &
                     ombb_l_asia)

         if (id_bcbb_l_asia > 0 ) &
           used = send_data ( id_bcbb_l_asia, bcbb_l_asia, Time )
         if (id_ombb_l_asia > 0 ) &
           used = send_data ( id_ombb_l_asia, ombb_l_asia, Time )
!
!
!##############I N D I A N EMISSIONS ##############################
! Indian emissions for BC/OC are from Reddy and Venkataraman,2002a; 2002b
! Reddy et al., 2002. The emissions are for fossil fuels: area and lrage point 
! sources, biomass burning: biofuels, forest fires/field burning of crop waste
! Emisions are for base year 1999. The  seasonality for forest fires are 
! imposed from ATSR fires [Reddy  and Boucher, 2004].
!
!OPEN BIOMASS BURNING 
!
         rdata  = 0.0
	 emiss2 = 0.0
	 emiss1 = 0.0

         open (30,file='INPUT/bcob_india.0.25x0.25', &
               form='formatted', action='read')

         DO WHILE (.not.tilleof)
         read  (30, FMT=2010,end=41)rlat,rlon,(rdata(itime),itime=1,12)
           i=181+INT(rlon-0.125)
	   j= 91+INT(rlat-0.125)  
	   emiss2(i,j,:) = emiss2(i,j,:)+rdata(:)
         ENDDO
2010    FORMAT(2(F6.2),12(E12.3))
41      call close_file (30)

        emiss1(:,:)= emiss2(:,:,month)/day_in_month(month) &
                                      /sec_in_day/grid_area(:,:)
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   bcbb_h_india)
         if (id_bcbb_h_india> 0 ) &
           used = send_data ( id_bcbb_h_india, bcbb_h_india, Time )
!
         rdata  = 0.0
	 emiss2 = 0.0
	 emiss1 = 0.0

         open (30,file='INPUT/omob_india.0.25x0.25', &
               form='formatted', action='read')

         DO WHILE (.not.tilleof)
         read  (30, FMT=2010,end=52)rlat,rlon,(rdata(itime),itime=1,12)
           i=181+INT(rlon-0.125)
	   j= 91+INT(rlat-0.125)  
	   emiss2(i,j,:) = emiss2(i,j,:)+rdata(:)
         ENDDO
52      call close_file (30)

        emiss1(:,:)= 1.6/1.3*emiss2(:,:,month)/day_in_month(month) &
                                             /sec_in_day/grid_area(:,:)

        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   ombb_h_india)
         if (id_ombb_h_india> 0 ) &
           used = send_data ( id_ombb_h_india, ombb_h_india, Time )
!
!-----------------------fly ash emissions-
!
         rdata  = 0.0
	 emiss2 = 0.0
	 emiss1 = 0.0

         open (30,file='INPUT/ifob_india.0.25x0.25', &
               form='formatted', action='read')

         DO WHILE (.not.tilleof)
         read  (30, FMT=2010,end=42)rlat,rlon,(rdata(itime),itime=1,12)
           i=181+INT(rlon-0.125)
	   j= 91+INT(rlat-0.125)  
	   emiss2(i,j,:) = emiss2(i,j,:)+rdata(:)
         ENDDO
42      call close_file (30)
        emiss1(:,:)= emiss2(:,:,month)/day_in_month(month) &
                                             /sec_in_day/grid_area(:,:)

        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   fabb_h_india)
         if (id_fabb_h_india> 0 ) &
           used = send_data ( id_fabb_h_india, fabb_h_india, Time )

if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading bc/om/fa BB_H India'
endif
!
! Fossil fuels
!
         rdata  = 0.0
	 bcl=0.0; bch = 0.0
	 ocl=0.0; och = 0.0
!
         open (30,file='INPUT/ff_india.0.25x0.25', &
               form='formatted', action='read')

         DO WHILE (.not.tilleof)
         read  (30, FMT=2011,end=43)rlat,rlon, &
	                (rdata(itime),itime=1,10)
           i=181+INT(rlon-0.125)
	   j= 91+INT(rlat-0.125)  
!
	   bcl(i,j) = bcl(i,j)+rdata(3)
	   bch(i,j) = bch(i,j)+rdata(8)
	   ocl(i,j) = ocl(i,j)+rdata(4)
	   och(i,j) = och(i,j)+rdata(9)
	   fal(i,j) = fal(i,j)+rdata(5)
	   fah(i,j) = fah(i,j)+rdata(10)

         ENDDO
43      call close_file (30)
2011    FORMAT(2(F6.2),10(E12.3))
!
        emiss1(:,:)= bcl(:,:)/day_in_year/sec_in_day/grid_area(:,:)
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   bcff_l_india)
!
        emiss1(:,:)= bch(:,:)/day_in_year/sec_in_day/grid_area(:,:)
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   bcff_h_india)
!
        emiss1(:,:)= om_oc_ff/1.3*ocl(:,:)/day_in_year/sec_in_day/grid_area(:,:)
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   omff_l_india)
!
        emiss1(:,:)= om_oc_ff/1.3*och(:,:)/day_in_year/sec_in_day/grid_area(:,:)
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   omff_h_india)
!
        emiss1(:,:)= fal(:,:)/day_in_year/sec_in_day/grid_area(:,:)
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   faff_l_india)
!
        emiss1(:,:)= fah(:,:)/day_in_year/sec_in_day/grid_area(:,:)
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   faff_h_india)
!
         if (id_bcff_h_india> 0 ) &
           used = send_data ( id_bcff_h_india, bcff_h_india, Time )

         if (id_bcff_l_india> 0 ) &
           used = send_data ( id_bcff_l_india, bcff_l_india, Time )

         if (id_omff_h_india> 0 ) &
           used = send_data ( id_omff_h_india, omff_h_india, Time )

         if (id_omff_l_india> 0 ) &
           used = send_data ( id_omff_l_india, omff_l_india, Time )

         if (id_faff_h_india> 0 ) &
           used = send_data ( id_faff_h_india, faff_h_india, Time )

         if (id_faff_l_india> 0 ) &
           used = send_data ( id_faff_l_india, faff_l_india, Time )
!
         rdata = 0.0
	 emiss1 = 0.0
	 bcl=0.0; ocl=0.0

         open (30,file='INPUT/bf_india.0.25x0.25', &
               form='formatted', action='read')

         DO WHILE (.not.tilleof)
         read  (30, FMT=2012,end=44)rlat,rlon, (rdata(itime),itime=1,5)
           i=181+INT(rlon-0.125)
	   j= 91+INT(rlat-0.125)  

	   bcl(i,j) = bcl(i,j)+rdata(3)
	   ocl(i,j) = ocl(i,j)+rdata(4)
	   fal(i,j) = fal(i,j)+rdata(5)
         ENDDO
44      call close_file (30)
2012    FORMAT(2(F6.2),5(E12.3))
!
        emiss1(:,:)= bcl(:,:)/day_in_year/sec_in_day/grid_area(:,:)
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   bcbb_l_india)
!
        emiss1(:,:)= 1.6/1.3*ocl(:,:)/day_in_year/sec_in_day/grid_area(:,:)
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   ombb_l_india)
!
        emiss1(:,:)= fal(:,:)/day_in_year/sec_in_day/grid_area(:,:)
        call interp_emiss (emiss1, lon_W, lat_S, dlon, dlat,  &
	                   fabb_l_india)
!
         if (id_bcbb_l_india> 0 ) &
           used = send_data ( id_bcbb_l_india, bcbb_l_india, Time )
         if (id_ombb_h_india> 0 ) &
           used = send_data ( id_ombb_l_india, ombb_l_india, Time )
         if (id_fabb_l_india> 0 ) &
           used = send_data ( id_fabb_l_india, fabb_l_india, Time )
!
!################### CHOICE OF EMISSIONS TO THE TRANSPORT ###################
! Following choices of emission inventories are possible
! bcff_source = 1 Cooke and Wilson [1996]
!             = 2 Cooke et al. [1999]
!             = 3 Bond et al. [2004]
! ocff_source = 2 Cooke et al. [1999]
!             = 3 Bond et al. [2004]
! bcbb_source = 1 Cooke and Wilson [1996] with all emissions 
!                 are emitted at surface 
!             = 3 Bond et al. [2004] Open burning emissions seasonality is 
!                 inferred from ATSR fire counts
!             = 4 Cooke and Wilson [1996] are splitted into surface/
!                 elevated sources with ATSR fire coutns
!             = 5 BF emissions are from Tami Bond;
!                 Open burning emissions seasonalilty is from van der Werf 
! ocbb_source = 1 Assumed OC/BC 7.0  for Cooke and Wilson [1996] 
!                 BC emissions with all emissions are emitted at surface
!             = 3 Bond et al. [2004] Open burning emissions seasonality is 
!                 inferred from ATSR fire counts
!             = 4 Assuemed oc/bc=7 from  Cooke and Wilson [1996] and are 
!                 splitted into surface/elevated sources with ATSR fire coutns
!             = 5 BF emissions are from Tami Bond;
!                 Open burning emissions seasonalilty is from van der Werf 
!##########################################################################

        IF(bcff_source.EQ.1) THEN
           bcemisff_l(:,:) = bcff_cw96(:,:)
           bcemisff_h(:,:) = 0.0
        ELSEIF(bcff_source.EQ.2) THEN
           bcemisff_l(:,:) = bcff_cooke99(:,:)
           bcemisff_h(:,:) = 0.0
        ELSEIF(bcff_source.EQ.3) THEN
           bcemisff_l(:,:) = bcff_bond(:,:)
           bcemisff_h(:,:) = 0.0
        ENDIF 
!
        IF(ocff_source.EQ.2) THEN
           omemisff_l(:,:) = omff_cooke99(:,:)
           omemisff_h(:,:) = 0.0
        ELSEIF(ocff_source.EQ.3) THEN
           omemisff_l(:,:) = omff_bond(:,:)
           omemisff_h(:,:) = 0.0
        ENDIF 
!
        IF(bcbb_source.EQ.1) THEN
           bcemisbb_l(:,:) = bcbb_cw96(:,:)
           bcemisbb_h(:,:) = 0.0
        ELSEIF(bcbb_source.EQ.3) THEN
           bcemisbb_l(:,:) = bcbf_bond(:,:)
           bcemisbb_h(:,:) = bcob_bond(:,:)
        ELSEIF(bcbb_source.EQ.4) THEN
           bcemisbb_l(:,:) = bcbb_l_atsr(:,:)
           bcemisbb_h(:,:) = bcbb_h_atsr(:,:)
        ENDIF 
!
        IF(ocbb_source.EQ.1) THEN
           omemisbb_l(:,:) = bcbb_cw96(:,:)*7.0
           omemisbb_h(:,:) = 0.0
        ELSEIF(ocbb_source.EQ.3) THEN
           omemisbb_l(:,:) = ombf_bond(:,:)
           omemisbb_h(:,:) = omob_bond(:,:)
        ELSEIF(ocbb_source.EQ.4) THEN
           omemisbb_l(:,:) = bcbb_l_atsr(:,:)*7.0*1.6
           omemisbb_h(:,:) = bcbb_h_atsr(:,:)*7.0*1.6
        ENDIF 
!
!# selection of Asian/Indian emission inventories
!
!
        IF(streets_asia) THEN
	    where(bcbb_l_asia.GT.0.0) bcemisbb_l=bcbb_l_asia
	    where(bcbb_h_asia.GT.0.0) bcemisbb_h=bcbb_h_asia
	    where(bcff_l_asia.GT.0.0) bcemisff_l=bcbb_l_asia
	    where(bcff_h_asia.GT.0.0) bcemisff_h=bcbb_h_asia
!
	    where(ombb_l_asia.GT.0.0) omemisbb_l=ombb_l_asia
	    where(ombb_h_asia.GT.0.0) omemisbb_h=ombb_h_asia
	    where(omff_l_asia.GT.0.0) omemisff_l=ombb_l_asia
	    where(omff_h_asia.GT.0.0) omemisff_h=ombb_h_asia
	ENDIF

        IF(rv02_india) THEN
	    where(bcbb_l_india.GT.0.0) bcemisbb_l=bcbb_l_india
	    where(bcbb_h_india.GT.0.0) bcemisbb_h=bcbb_h_india
	    where(bcff_l_india.GT.0.0) bcemisff_l=bcbb_l_india
	    where(bcff_h_india.GT.0.0) bcemisff_h=bcbb_h_india
!
	    where(ombb_l_india.GT.0.0) omemisbb_l=ombb_l_india
	    where(ombb_h_india.GT.0.0) omemisbb_h=ombb_h_india
	    where(omff_l_india.GT.0.0) omemisff_l=ombb_l_india
	    where(omff_h_india.GT.0.0) omemisff_h=ombb_h_india
!
	ENDIF
!initialise fly-ash emissions 
!
        IF(rv02_india) THEN
	    faemisbb_l=fabb_l_india
	    faemisbb_h=fabb_h_india
	    faemisff_l=faff_l_india
	    faemisff_h=faff_h_india
	    else
	    faemisbb_l=0.0
	    faemisbb_h=0.0
	    faemisff_l=0.0
	    faemisff_h=0.0
	ENDIF
!
!
!################## NATURAL BIOGENIC EMISSIONS ###########################
!
! Reading natural Terpene emissions natural sources [Guenther et al., 1996]
! Monthly mean emissions are ginven as mg/m2/month
!
do unit = 30,100
INQUIRE(unit=unit, opened= opened)
if (.NOT. opened) exit
enddo
         emiss1= 0.0
         emiss2= 0.0
!
         open (unit,file='INPUT/terpen.nat.1x1', &
               form='formatted', action='read', iostat=ios )
         if (ios.eq.0 ) then
	 DO WHILE (.not.tilleof)
         read  (unit, FMT=3000, end=99) j,i, &
	                     (emiss2(i,j,itime),itime=1,12)
	 ENDDO
         else
         write(*,*) '***ERROR: Opening Terpene emissions file'
         endif
99      call close_file (unit)
3000    FORMAT(2(I3),12(E11.5))
        emiss1(:,:) =  emiss2(:,:,month)
        emiss1 = emiss1*1.E-6/day_in_month(month)/sec_in_day    !kgC /m2/sec
!
!      It is assumed that 11% of terpene emissions are converted to 
!      OC. Once again an OM to OC ratio of 1.4 is assumed.
       emiss1 = emiss1*0.11*om_oc_nat

       call interp_emiss ( emiss1, lon_W, lat_S, dlon, dlat, &
                     omemisnat)
if (mpp_pe()== mpp_root_pe() ) then
   write(*,*) 'Reading terpen emissions'
endif

!
!#############################################
! select source depending upon runttype
!
   if(trim(runtype) == 'fossil_fuels_only') then
           bcemisbb_l = 0.0
           bcemisbb_h = 0.0
           omemisbb_l = 0.0
           omemisbb_h = 0.0
           faemisbb_l = 0.0
           faemisbb_h = 0.0
           omemisnat  = 0.0
   else if(trim(runtype) == 'biomass_only') then
           bcemisff_l = 0.0
           bcemisff_h = 0.0
           omemisff_l = 0.0
           omemisff_h = 0.0
	   faemisff_l = 0.0
	   faemisff_h = 0.0
           omemisnat  = 0.0
   else if(trim(runtype) == 'natural_only') then
           bcemisbb_l = 0.0
           bcemisbb_h = 0.0
           omemisbb_l = 0.0
           omemisbb_h = 0.0
           faemisbb_l = 0.0
           faemisbb_h = 0.0
!
           bcemisff_l = 0.0
           bcemisff_h = 0.0
           omemisff_l = 0.0
           omemisff_h = 0.0
	   faemisff_l = 0.0
	   faemisff_h = 0.0
   else if(trim(runtype) == 'anthrop') then
   else
        call error_mesg ('atmos_carbon_aerosol_mod', &
             'runtype is not recognized', FATAL)
   endif
!#######################################!
!Write to output
!
         if (id_bcemisff_l > 0 ) &
           used = send_data ( id_bcemisff_l, bcemisff_l, Time )
         if (id_bcemisbb_l > 0 ) &
           used = send_data ( id_bcemisbb_l, bcemisbb_l, Time )
         if (id_omemisff_l > 0 ) &
           used = send_data ( id_omemisff_l, omemisff_l, Time )
         if (id_omemisbb_l > 0 ) &
           used = send_data ( id_omemisbb_l, omemisbb_l, Time )
!
         if (id_bcemisff_h > 0 ) &
           used = send_data ( id_bcemisff_h, bcemisff_h, Time )
         if (id_bcemisbb_h > 0 ) &
           used = send_data ( id_bcemisbb_h, bcemisbb_h, Time )
         if (id_omemisff_h > 0 ) &
           used = send_data ( id_omemisff_h, omemisff_h, Time )
         if (id_omemisbb_h > 0 ) &
           used = send_data ( id_omemisbb_h, omemisbb_h, Time )
!
         if (id_faemisff_h > 0 ) &
           used = send_data ( id_faemisff_h, faemisff_h, Time )
         if (id_faemisbb_h > 0 ) &
           used = send_data ( id_faemisbb_h, faemisbb_h, Time )
         if (id_faemisff_l > 0 ) &
           used = send_data ( id_faemisff_l, faemisff_l, Time )
         if (id_faemisbb_l > 0 ) &
           used = send_data ( id_faemisbb_l, faemisbb_l, Time )
         if (id_omemisnat > 0 ) &
           used = send_data ( id_omemisnat, omemisnat, Time )
!################ END OF CARBONACEOUS AEROSOL SOURCES ###################

end subroutine carbon_aerosol_source_input
end module atmos_carbon_aerosol_mod
