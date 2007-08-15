Module atmos_carbon_aerosol_mod

! <CONTACT EMAIL="Shekar.Reddy@noaa.gov">
!   Shekar Reddy
! </CONTACT>
use fms_mod,                    only : file_exist, close_file, &
                                       write_version_number, &
                                       mpp_pe, mpp_root_pE, &
                                       open_namelist_file,  &
                                       check_nml_error, error_mesg,  &
                                       stdlog, FATAL, NOTE, WARNING
use time_manager_mod,           only : time_type, &
                                       days_in_month, days_in_year, &
                                       set_date, set_time, get_date_julian, &
                                       print_date, get_date, &
                                       operator(>), operator(+), operator(-)
use time_interp_mod,            only:  fraction_of_year, &
                                       time_interp_init
use diag_manager_mod,           only : send_data, register_diag_field, &
                                       diag_manager_init, get_base_time
use tracer_manager_mod,         only : get_tracer_index, &
                                       set_tracer_atts
use field_manager_mod,          only : MODEL_ATMOS
use atmos_tracer_utilities_mod, only : wet_deposition,       &
                                       dry_deposition
use interpolator_mod,           only:  interpolate_type, interpolator_init, &
                                       interpolator, interpolator_end, &
                                       CONSTANT, INTERP_WEIGHTED_P
use constants_mod,              only : PI, GRAV, RDGAS
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
!--- Interpolate_type variable containing all the information needed to
! interpolate the emission provided in the netcdf input file.
type(interpolate_type),save  ::bcff_aerosol_interp
type(interpolate_type),save  ::bcbb_aerosol_interp
type(interpolate_type),save  ::bcob_aerosol_interp
type(interpolate_type),save  ::omff_aerosol_interp
type(interpolate_type),save  ::ombb_aerosol_interp
type(interpolate_type),save  ::omob_aerosol_interp
type(interpolate_type),save  ::omna_aerosol_interp
! Initial calendar time for model
type(time_type) :: model_init_time

! Difference between model initial time and dust source timeseries applied
! at model initial time
type(time_type), save :: bcff_offset
type(time_type), save :: bcbb_offset
type(time_type), save :: bcob_offset
type(time_type), save :: omff_offset
type(time_type), save :: ombb_offset
type(time_type), save :: omob_offset
type(time_type), save :: omna_offset

! Timein dust source timeseries which is mapped to model initial time
type(time_type), save :: bcff_entry
type(time_type), save :: bcbb_entry
type(time_type), save :: bcob_entry
type(time_type), save :: omff_entry
type(time_type), save :: ombb_entry
type(time_type), save :: omob_entry
type(time_type), save :: omna_entry

! The model initial time is later than the dust_dataset_entry time  ?
logical, save    :: bcff_negative_offset
logical, save    :: bcbb_negative_offset
logical, save    :: bcob_negative_offset
logical, save    :: omff_negative_offset
logical, save    :: ombb_negative_offset
logical, save    :: omob_negative_offset
logical, save    :: omna_negative_offset

integer, save    :: bcff_time_serie_type
integer, save    :: bcbb_time_serie_type
integer, save    :: bcob_time_serie_type
integer, save    :: omff_time_serie_type
integer, save    :: ombb_time_serie_type
integer, save    :: omob_time_serie_type
integer, save    :: omna_time_serie_type
!----------------------------------------------------------------------
!-------- namelist  ---------
character(len=26)  :: bcff_filename = 'carbon_aerosol_emission.nc'
character(len=26)  :: bcbb_filename = 'carbon_aerosol_emission.nc'
character(len=26)  :: bcob_filename = 'carbon_aerosol_emission.nc'
character(len=26)  :: omff_filename = 'carbon_aerosol_emission.nc'
character(len=26)  :: ombb_filename = 'carbon_aerosol_emission.nc'
character(len=26)  :: omob_filename = 'carbon_aerosol_emission.nc'
character(len=26)  :: omna_filename = 'carbon_aerosol_emission.nc'

character(len=12), dimension(3)  :: bcff_emission_name
data bcff_emission_name/'bcff_cw96','bcff_cooke99','bcff_bond'/
character(len=12), dimension(4)  :: bcbb_emission_name
data bcbb_emission_name/'bcbb_cw96','bcbf_bond','bc_geia1','bc_geia2'/
character(len=12), dimension(1)  :: bcob_emission_name
data bcob_emission_name/'bcob_bond'/
character(len=12), dimension(4)  :: omff_emission_name
data omff_emission_name/'omff_cooke99','omff_bond','om_geia1','om_geia2'/
character(len=12), dimension(4)  :: ombb_emission_name
data ombb_emission_name/'bcbb_cw96','ombf_bond','om_geia1','om_geia2'/
character(len=12), dimension(1)  :: omob_emission_name
data omob_emission_name/'omob_bond'/
character(len=12), dimension(1)  :: omna_emission_name
data omna_emission_name/'omemisnat'/
! Default values for carbon_aerosol_nml
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FOSSIL FUEL source can be either:
! for black carbon: 'cooke_and_wilson_1996' 
!                   'cooke_1999' 
!                   'bond_2004'
character(len=9)      :: bcff_source = 'bond_2004'
character(len=24)     :: bcff_time_dependency_type
integer, dimension(6) :: bcff_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fossil fuel emission 
! for organic matter: 'cooke_1999' 
!                     'bond_2004'
character(len=20)     :: omff_source = 'bond_2004'
character(len=24)     :: omff_time_dependency_type
integer, dimension(6) :: omff_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BIOMASS BURNING emissions can be either
! for black carbon: 'cooke_and_wilson_1996' 
!                   'bond_2004' 
!                   'reddy_and_boucher_2004'
character(len=20)     :: bcbb_source = 'bond_2004'
character(len=24)     :: bcbb_time_dependency_type
integer, dimension(6) :: bcbb_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
! biomass burning emission
! for organic matter: 'cooke_and_wilson_1996' 
!                     'bond_2004' 
!                     'reddy_and_boucher_2004'
character(len=20)     :: ombb_source = 'bond_2004'
character(len=24)     :: ombb_time_dependency_type
integer, dimension(6) :: ombb_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open burning emission
! for black carbon     'bond_2004' 
character(len=20)     :: bcob_source = 'bond_2004'
character(len=24)     :: bcob_time_dependency_type
integer, dimension(6) :: bcob_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
! for organic matter:  'bond_2004' 
character(len=20)     :: omob_source = 'bond_2004'
character(len=24)     :: omob_time_dependency_type
integer, dimension(6) :: omob_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for natural emssion: Gunther et al., 1996 inventory
character(len=20)     :: omna_source = 'gunther'
character(len=24)     :: omna_time_dependency_type
integer, dimension(6) :: omna_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
namelist /carbon_aerosol_nml/ &
 bcff_source,bcff_emission_name,bcff_time_dependency_type, bcff_dataset_entry, &
 bcbb_source,bcbb_emission_name,bcbb_time_dependency_type, bcbb_dataset_entry, &
 bcob_source,bcob_emission_name,bcob_time_dependency_type, bcob_dataset_entry, &
 omff_source,omff_emission_name,omff_time_dependency_type, omff_dataset_entry, &
 ombb_source,ombb_emission_name,ombb_time_dependency_type, ombb_dataset_entry, &
 omob_source,omob_emission_name,omob_time_dependency_type, omob_dataset_entry, &
 omna_source,omna_emission_name,omna_time_dependency_type, omna_dataset_entry

character(len=6), parameter :: module_name = 'tracer'

logical :: module_is_initialized = .FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id: atmos_carbon_aerosol.F90,v 15.0 2007/08/14 03:56:43 fms Exp $'
character(len=128) :: tagname = '$Name: omsk $'
!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine atmos_carbon_aerosol_driver(lon, lat, land, pfull,phalf,z_half, z_pbl, &
                               T, pwt, &
                               black_cphob, black_cphob_dt,  &
                               black_cphil, black_cphil_dt,  &
                               omphob, omphob_dt, &
                               omphil, omphil_dt, &
                               model_time, is, ie, js, je )

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
type(time_type), intent(in)            :: model_time
integer, intent(in)                    :: is, ie, js, je
!-----------------------------------------------------------------------
type(time_type)                        :: bcff_time
type(time_type)                        :: bcbb_time
type(time_type)                        :: bcob_time
type(time_type)                        :: omff_time
type(time_type)                        :: ombb_time
type(time_type)                        :: omob_time
type(time_type)                        :: omna_time
real, dimension(size(black_cphob,1),size(black_cphob,2),size(black_cphob,3)) ::&
   bc_sourcephob, bc_sourcephil, om_sourcephob, om_sourcephil, sinkphob
   real  dtr,bltop,z1,z2
real, dimension(size(black_cphob,3)) :: fbb, fff
integer        :: i,j,l, kb,id,jd,kd,lat1,k
integer        :: yr, mo, dy, hr, mn, sc, dum, mo_yr, dayspmn
REAL,PARAMETER :: frac_bc_phobic = 0.8
REAL,PARAMETER :: frac_bc_philic = 0.2
REAL,PARAMETER :: frac_om_phobic = 0.5
REAL,PARAMETER :: frac_om_philic = 0.5

!
!-------------------------------------------------------------------------
   real, dimension(size(black_cphob,1),size(black_cphob,2)) ::          &
     bcemisff_l, bcemisff_h, bcemisbb_l, bcemisbb_h, bcemisob, &
     omemisff_l, omemisff_h, omemisbb_l, omemisbb_h, omemisob, omemisnat
!
!-----------------------------------------------------------------------
!
      id=size(black_cphob,1); jd=size(black_cphob,2); kd=size(black_cphob,3)

      dtr= PI/180.

!----------- compute black carbon source ------------

    bc_sourcephob = 0.0
    bc_sourcephil = 0.0
    black_cphob_dt = 0.0
    black_cphil_dt = 0.0
    om_sourcephob = 0.0
    om_sourcephil = 0.0
    omphob_dt = 0.0
    omphil_dt = 0.0
! emission rates (kg/m2/s)
    bcemisff_l(:,:) = 0.0
    bcemisff_h(:,:) = 0.0
    bcemisbb_l(:,:) = 0.0
    bcemisbb_h(:,:) = 0.0
    bcemisob(:,:)  = 0.0

    omemisff_l(:,:) = 0.0
    omemisff_h(:,:) = 0.0
    omemisbb_l(:,:) = 0.0
    omemisbb_h(:,:) = 0.0
    omemisob(:,:)  = 0.0
    omemisnat(:,:)  = 0.0

!--------------------------------------------------------------------
!    define the time in the bcff data set from which data is to be 
!    taken. if bcff is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(bcff_time_serie_type .eq. 3) then
       if (bcff_negative_offset) then
         bcff_time = model_time - bcff_offset
       else
         bcff_time = model_time + bcff_offset
       endif
     else 
       if(bcff_time_serie_type .eq. 2 ) then
         call get_date (bcff_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(bcff_entry)
           if (dayspmn /= 29) then
             bcff_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             bcff_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           bcff_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         bcff_time = model_time
       endif
     endif
    select case (trim(bcff_source))
      case ('cooke_and_wilson_1996')
             call interpolator(bcff_aerosol_interp, bcff_time, bcemisff_l, &
                  trim(bcff_emission_name(1)), is, js)
      case ('cooke_1999')
             call interpolator(bcff_aerosol_interp, bcff_time, bcemisff_l, &
                  trim(bcff_emission_name(2)), is, js)
      case ('bond_2004') 
             call interpolator(bcff_aerosol_interp, bcff_time, bcemisff_l, &
                  trim(bcff_emission_name(3)), is, js)
    end select

!--------------------------------------------------------------------
!    define the time in the bcbb data set from which data is to be 
!    taken. if bcbb is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(bcbb_time_serie_type .eq. 3) then
       if (bcbb_negative_offset) then
         bcbb_time = model_time - bcbb_offset
       else
         bcbb_time = model_time + bcbb_offset
       endif
     else 
       if(bcbb_time_serie_type .eq. 2 ) then
         call get_date (bcbb_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(bcbb_entry)
           if (dayspmn /= 29) then
             bcbb_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             bcbb_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           bcbb_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         bcbb_time = model_time
       endif
     endif
    select case (trim(bcbb_source))
      case ('cooke_and_wilson_1996')
        call interpolator(bcbb_aerosol_interp, bcbb_time, bcemisbb_l, &
                         trim(bcbb_emission_name(1)), is, js)
      case ('bond_2004') 
        call interpolator(bcbb_aerosol_interp, bcbb_time, bcemisbb_l, &
                         trim(bcbb_emission_name(2)), is, js)
      case ('reddy_and_boucher_2004') 
        call interpolator(bcbb_aerosol_interp, bcbb_time, bcemisbb_l, &
                         trim(bcbb_emission_name(3)), is, js)
        call interpolator(bcbb_aerosol_interp, bcbb_time, bcemisbb_h, &
                         trim(bcbb_emission_name(4)), is, js)
    end select
!--------------------------------------------------------------------
!    define the time in the bcob data set from which data is to be 
!    taken. if bcob is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(bcob_time_serie_type .eq. 3) then
       if (bcob_negative_offset) then
         bcob_time = model_time - bcob_offset
       else
         bcob_time = model_time + bcob_offset
       endif
     else 
       if(bcob_time_serie_type .eq. 2 ) then
         call get_date (bcob_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(bcob_entry)
           if (dayspmn /= 29) then
             bcob_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             bcob_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           bcob_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         bcob_time = model_time
       endif
     endif
    select case (trim(bcob_source))
      case ('bond_2004') 
        call interpolator(bcob_aerosol_interp, bcob_time, bcemisob, &
                       trim(bcob_emission_name(1)), is, js)
    end select
!--------------------------------------------------------------------
!    define the time in the omff data set from which data is to be 
!    taken. if omff is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(omff_time_serie_type .eq. 3) then
       if (omff_negative_offset) then
         omff_time = model_time - omff_offset
       else
         omff_time = model_time + omff_offset
       endif
     else 
       if(omff_time_serie_type .eq. 2 ) then
         call get_date (omff_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(omff_entry)
           if (dayspmn /= 29) then
             omff_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             omff_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           omff_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         omff_time = model_time
       endif
     endif
    select case (trim(omff_source))
      case ('cooke_1999')
          call interpolator(omff_aerosol_interp, omff_time, omemisff_l, &
                      trim(omff_emission_name(1)), is, js)
      case ('bond_2004') 
          call interpolator(omff_aerosol_interp, omff_time, omemisff_l, &
                      trim(omff_emission_name(2)), is, js)
    end select


!--------------------------------------------------------------------
!    define the time in the ombb data set from which data is to be 
!    taken. if ombb is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(ombb_time_serie_type .eq. 3) then
       if (ombb_negative_offset) then
         ombb_time = model_time - ombb_offset
       else
         ombb_time = model_time + ombb_offset
       endif
     else 
       if(ombb_time_serie_type .eq. 2 ) then
         call get_date (ombb_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(ombb_entry)
           if (dayspmn /= 29) then
             ombb_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             ombb_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           ombb_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         ombb_time = model_time
       endif
     endif
    select case (trim(ombb_source))
      case ('cooke_and_wilson_1996')
        call interpolator(ombb_aerosol_interp, ombb_time, omemisbb_l, &
                           trim(ombb_emission_name(1)), is, js)
        omemisbb_l(:,:) = omemisbb_l(:,:)*7.0
        omemisbb_h(:,:) = 0.0
      case ('bond_2004') 
        call interpolator(ombb_aerosol_interp, ombb_time, omemisbb_l, &
                       trim(ombb_emission_name(2)), is, js)
      case ('reddy_and_boucher_2004') 
          call interpolator(ombb_aerosol_interp, ombb_time, omemisbb_l, &
                       trim(ombb_emission_name(3)), is, js)
          call interpolator(ombb_aerosol_interp, ombb_time, omemisbb_h, &
                       trim(ombb_emission_name(4)), is, js)
    end select
!--------------------------------------------------------------------
!    define the time in the omob data set from which data is to be 
!    taken. if omob is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(omob_time_serie_type .eq. 3) then
       if (omob_negative_offset) then
         omob_time = model_time - omob_offset
       else
         omob_time = model_time + omob_offset
       endif
     else 
       if(omob_time_serie_type .eq. 2 ) then
         call get_date (omob_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(omob_entry)
           if (dayspmn /= 29) then
             omob_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             omob_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           omob_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         omob_time = model_time
       endif
     endif
    select case (trim(omob_source))
      case ('bond_2004') 
        call interpolator(omob_aerosol_interp, omob_time, omemisob, &
                       trim(omob_emission_name(1)), is, js)
    end select
!--------------------------------------------------------------------
!    define the time in the omna data set from which data is to be 
!    taken. if omna is not time-varying, it is simply model_time.
!---------------------------------------------------------------------
     if(omna_time_serie_type .eq. 3) then
       if (omna_negative_offset) then
         omna_time = model_time - omna_offset
       else
         omna_time = model_time + omna_offset
       endif
     else 
       if(omna_time_serie_type .eq. 2 ) then
         call get_date (omna_entry, yr, dum,dum,dum,dum,dum)
         call get_date (model_time, mo_yr, mo, dy, hr, mn, sc)
         if (mo ==2 .and. dy == 29) then
           dayspmn = days_in_month(omna_entry)
           if (dayspmn /= 29) then
             omna_time = set_date (yr, mo, dy-1, hr, mn, sc)
           else
             omna_time = set_date (yr, mo, dy, hr, mn, sc)
           endif
         else
           omna_time = set_date (yr, mo, dy, hr, mn, sc)
         endif
       else
         omna_time = model_time
       endif
     endif
    select case (trim(omna_source))
      case ('gunther')
        call interpolator(omna_aerosol_interp, model_time, omemisnat, &
                       trim(omna_emission_name(1)), is, js)
    end select
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
          *(bcemisff_l(i,j)+ bcemisbb_l(i,j) + bcemisob(i,j) )/pwt(i,j,kd)

        bc_sourcephil(i,j,kd) =  frac_bc_philic  &
          *(bcemisff_l(i,j)+ bcemisbb_l(i,j) + bcemisob(i,j) )/pwt(i,j,kd)

        bc_sourcephob(i,j,:) =  bc_sourcephob(i,j,:) + frac_bc_phobic &
          *(fff(:)*bcemisff_h(i,j)+fbb(:)*bcemisbb_h(i,j)) /pwt(i,j,:)
        bc_sourcephil(i,j,:) =  bc_sourcephil(i,j,:) + frac_bc_philic &
          *(fff(:)*bcemisff_h(i,j)+fbb(:)*bcemisbb_h(i,j)) /pwt(i,j,:)
!
        om_sourcephob(i,j,kd) =  frac_om_phobic  &
          *(omemisff_l(i,j)+ omemisbb_l(i,j) + omemisob(i,j) &
            +omemisnat(i,j))/pwt(i,j,kd)

        om_sourcephil(i,j,kd) =  frac_om_philic  &
          *(omemisff_l(i,j)+ omemisbb_l(i,j) + omemisob(i,j) &
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
real, dimension(:,:),    intent(in) :: lonb, latb
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

!----------------------------------------------------------------------
!    initialize namelist entries
!----------------------------------------------------------------------
        bcff_offset = set_time (0,0)
        bcbb_offset = set_time (0,0)
        bcob_offset = set_time (0,0)
        omff_offset = set_time (0,0)
        ombb_offset = set_time (0,0)
        omob_offset = set_time (0,0)
        omna_offset = set_time (0,0)

        bcff_entry = set_time (0,0)
        bcbb_entry = set_time (0,0)
        bcob_entry = set_time (0,0)
        omff_entry = set_time (0,0)
        ombb_entry = set_time (0,0)
        omob_entry = set_time (0,0)
        omna_entry = set_time (0,0)

        bcff_negative_offset = .false.
        bcbb_negative_offset = .false.
        bcob_negative_offset = .false.
        omff_negative_offset = .false.
        ombb_negative_offset = .false.
        omob_negative_offset = .false.
        omna_negative_offset = .false.

        bcff_time_serie_type = 1
        bcbb_time_serie_type = 1
        bcob_time_serie_type = 1
        omff_time_serie_type = 1
        ombb_time_serie_type = 1
        omob_time_serie_type = 1
        omna_time_serie_type = 1
!----------------------------------------------------------------------
!    define the model base time  (defined in diag_table)
!----------------------------------------------------------------------
        model_init_time = get_base_time()
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(bcff_time_dependency_type) == 'constant' ) then
        bcff_time_serie_type = 1
        bcff_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcff are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for bcff is selected.
!---------------------------------------------------------------------
      else if (trim(bcff_time_dependency_type) == 'time_varying') then
        bcff_time_serie_type = 3
        if (bcff_dataset_entry(1) == 1 .and. &
            bcff_dataset_entry(2) == 1 .and. &
            bcff_dataset_entry(3) == 1 .and. &
            bcff_dataset_entry(4) == 0 .and. &
            bcff_dataset_entry(5) == 0 .and. &
            bcff_dataset_entry(6) == 0 ) then
          bcff_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcff_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          bcff_entry  = set_date (bcff_dataset_entry(1), &
                                  bcff_dataset_entry(2), &
                                  bcff_dataset_entry(3), &
                                  bcff_dataset_entry(4), &
                                  bcff_dataset_entry(5), &
                                  bcff_dataset_entry(6))
        endif
        call print_date (bcff_entry , str= &
          'Data from bcff timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        bcff_offset = bcff_entry - model_init_time
        if (model_init_time > bcff_entry) then
          bcff_negative_offset = .true.
        else
          bcff_negative_offset = .false.
        endif
      else if (trim(bcff_time_dependency_type) == 'fixed_year') then
        bcff_time_serie_type = 2
        if (bcff_dataset_entry(1) == 1 .and. &
            bcff_dataset_entry(2) == 1 .and. &
            bcff_dataset_entry(3) == 1 .and. &
            bcff_dataset_entry(4) == 0 .and. &
            bcff_dataset_entry(5) == 0 .and. &
            bcff_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set bcff_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcff_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        bcff_entry  = set_date (bcff_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'bcff is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcff correspond to year :', &
                    bcff_dataset_entry(1)
        endif
     endif
     call interpolator_init (bcff_aerosol_interp,             &
                             trim(bcff_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = bcff_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(bcbb_time_dependency_type) == 'constant' ) then
        bcbb_time_serie_type = 1
        bcbb_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcbb are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for bcbb is selected.
!---------------------------------------------------------------------
      else if (trim(bcbb_time_dependency_type) == 'time_varying') then
        bcbb_time_serie_type = 3
        if (bcbb_dataset_entry(1) == 1 .and. &
            bcbb_dataset_entry(2) == 1 .and. &
            bcbb_dataset_entry(3) == 1 .and. &
            bcbb_dataset_entry(4) == 0 .and. &
            bcbb_dataset_entry(5) == 0 .and. &
            bcbb_dataset_entry(6) == 0 ) then
          bcbb_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcbb_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          bcbb_entry  = set_date (bcbb_dataset_entry(1), &
                                  bcbb_dataset_entry(2), &
                                  bcbb_dataset_entry(3), &
                                  bcbb_dataset_entry(4), &
                                  bcbb_dataset_entry(5), &
                                  bcbb_dataset_entry(6))
        endif
        call print_date (bcbb_entry , str= &
          'Data from bcbb timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        bcbb_offset = bcbb_entry - model_init_time
        if (model_init_time > bcbb_entry) then
          bcbb_negative_offset = .true.
        else
          bcbb_negative_offset = .false.
        endif
      else if (trim(bcbb_time_dependency_type) == 'fixed_year') then
        bcbb_time_serie_type = 2
        if (bcbb_dataset_entry(1) == 1 .and. &
            bcbb_dataset_entry(2) == 1 .and. &
            bcbb_dataset_entry(3) == 1 .and. &
            bcbb_dataset_entry(4) == 0 .and. &
            bcbb_dataset_entry(5) == 0 .and. &
            bcbb_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set bcbb_dataset_entry when using fixed_year source', FATAL)
        endif
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcbb_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        bcbb_entry  = set_date (bcbb_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'bcbb is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcff correspond to year :', &
                    bcbb_dataset_entry(1)
        endif
     endif
     call interpolator_init (bcbb_aerosol_interp,             &
                             trim(bcbb_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = bcbb_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(bcob_time_dependency_type) == 'constant' ) then
        bcob_time_serie_type = 1
        bcob_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcob are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for bcob is selected.
!---------------------------------------------------------------------
      else if (trim(bcob_time_dependency_type) == 'time_varying') then
        bcob_time_serie_type = 3
        if (bcob_dataset_entry(1) == 1 .and. &
            bcob_dataset_entry(2) == 1 .and. &
            bcob_dataset_entry(3) == 1 .and. &
            bcob_dataset_entry(4) == 0 .and. &
            bcob_dataset_entry(5) == 0 .and. &
            bcob_dataset_entry(6) == 0 ) then
          bcob_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcob_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          bcob_entry  = set_date (bcob_dataset_entry(1), &
                                  bcob_dataset_entry(2), &
                                  bcob_dataset_entry(3), &
                                  bcob_dataset_entry(4), &
                                  bcob_dataset_entry(5), &
                                  bcob_dataset_entry(6))
        endif
        call print_date (bcob_entry , str= &
          'Data from bcob timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        bcob_offset = bcob_entry - model_init_time
        if (model_init_time > bcob_entry) then
          bcob_negative_offset = .true.
        else
          bcob_negative_offset = .false.
        endif
      else if (trim(bcob_time_dependency_type) == 'fixed_year') then
        bcob_time_serie_type = 2
        if (bcob_dataset_entry(1) == 1 .and. &
            bcob_dataset_entry(2) == 1 .and. &
            bcob_dataset_entry(3) == 1 .and. &
            bcob_dataset_entry(4) == 0 .and. &
            bcob_dataset_entry(5) == 0 .and. &
            bcob_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set bcob_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to bcob_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        bcob_entry  = set_date (bcob_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'bcob is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'bcob correspond to year :', &
                    bcob_dataset_entry(1)
        endif
     endif
     call interpolator_init (bcob_aerosol_interp,             &
                             trim(bcob_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = bcob_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(omff_time_dependency_type) == 'constant' ) then
        omff_time_serie_type = 1
        omff_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omff are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for omff is selected.
!---------------------------------------------------------------------
      else if (trim(omff_time_dependency_type) == 'time_varying') then
        omff_time_serie_type = 3
        if (omff_dataset_entry(1) == 1 .and. &
            omff_dataset_entry(2) == 1 .and. &
            omff_dataset_entry(3) == 1 .and. &
            omff_dataset_entry(4) == 0 .and. &
            omff_dataset_entry(5) == 0 .and. &
            omff_dataset_entry(6) == 0 ) then
          omff_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omff_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          omff_entry  = set_date (omff_dataset_entry(1), &
                                  omff_dataset_entry(2), &
                                  omff_dataset_entry(3), &
                                  omff_dataset_entry(4), &
                                  omff_dataset_entry(5), &
                                  omff_dataset_entry(6))
        endif
        call print_date (omff_entry , str= &
          'Data from omff timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        omff_offset = omff_entry - model_init_time
        if (model_init_time > omff_entry) then
          omff_negative_offset = .true.
        else
          omff_negative_offset = .false.
        endif
      else if (trim(omff_time_dependency_type) == 'fixed_year') then
        omff_time_serie_type = 2
        if (omff_dataset_entry(1) == 1 .and. &
            omff_dataset_entry(2) == 1 .and. &
            omff_dataset_entry(3) == 1 .and. &
            omff_dataset_entry(4) == 0 .and. &
            omff_dataset_entry(5) == 0 .and. &
            omff_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set omff_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omff_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        omff_entry  = set_date (omff_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'omff is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omff correspond to year :', &
                    omff_dataset_entry(1)
        endif
     endif
     call interpolator_init (omff_aerosol_interp,             &
                             trim(omff_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = omff_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(bcob_time_dependency_type) == 'constant' ) then
        ombb_time_serie_type = 1
        ombb_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'ombb are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for ombb is selected.
!---------------------------------------------------------------------
      else if (trim(ombb_time_dependency_type) == 'time_varying') then
        ombb_time_serie_type = 3
        if (ombb_dataset_entry(1) == 1 .and. &
            ombb_dataset_entry(2) == 1 .and. &
            ombb_dataset_entry(3) == 1 .and. &
            ombb_dataset_entry(4) == 0 .and. &
            ombb_dataset_entry(5) == 0 .and. &
            ombb_dataset_entry(6) == 0 ) then
          ombb_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to ombb_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          ombb_entry  = set_date (ombb_dataset_entry(1), &
                                  ombb_dataset_entry(2), &
                                  ombb_dataset_entry(3), &
                                  ombb_dataset_entry(4), &
                                  ombb_dataset_entry(5), &
                                  ombb_dataset_entry(6))
        endif
        call print_date (ombb_entry , str= &
          'Data from ombb timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        ombb_offset = ombb_entry - model_init_time
        if (model_init_time > ombb_entry) then
          ombb_negative_offset = .true.
        else
          ombb_negative_offset = .false.
        endif
      else if (trim(ombb_time_dependency_type) == 'fixed_year') then
        ombb_time_serie_type = 2
        if (ombb_dataset_entry(1) == 1 .and. &
            ombb_dataset_entry(2) == 1 .and. &
            ombb_dataset_entry(3) == 1 .and. &
            ombb_dataset_entry(4) == 0 .and. &
            ombb_dataset_entry(5) == 0 .and. &
            ombb_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set ombb_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to ombb_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        ombb_entry  = set_date (ombb_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'ombb is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'ombb correspond to year :', &
                    ombb_dataset_entry(1)
        endif
     endif
     call interpolator_init (ombb_aerosol_interp,             &
                             trim(ombb_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = ombb_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(omob_time_dependency_type) == 'constant' ) then
        omob_time_serie_type = 1
        omob_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omob are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for omob is selected.
!---------------------------------------------------------------------
      else if (trim(omob_time_dependency_type) == 'time_varying') then
        omob_time_serie_type = 3
        if (omob_dataset_entry(1) == 1 .and. &
            omob_dataset_entry(2) == 1 .and. &
            omob_dataset_entry(3) == 1 .and. &
            omob_dataset_entry(4) == 0 .and. &
            omob_dataset_entry(5) == 0 .and. &
            omob_dataset_entry(6) == 0 ) then
          omob_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omob_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          omob_entry  = set_date (omob_dataset_entry(1), &
                                  omob_dataset_entry(2), &
                                  omob_dataset_entry(3), &
                                  omob_dataset_entry(4), &
                                  omob_dataset_entry(5), &
                                  omob_dataset_entry(6))
        endif
        call print_date (omob_entry , str= &
          'Data from omob timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        omob_offset = omob_entry - model_init_time
        if (model_init_time > omob_entry) then
          omob_negative_offset = .true.
        else
          omob_negative_offset = .false.
        endif
      else if (trim(omob_time_dependency_type) == 'fixed_year') then
        omob_time_serie_type = 2
        if (omob_dataset_entry(1) == 1 .and. &
            omob_dataset_entry(2) == 1 .and. &
            omob_dataset_entry(3) == 1 .and. &
            omob_dataset_entry(4) == 0 .and. &
            omob_dataset_entry(5) == 0 .and. &
            omob_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set omob_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omob_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        omob_entry  = set_date (omob_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'omob is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omob correspond to year :', &
                    omob_dataset_entry(1)
        endif
     endif
     call interpolator_init (omob_aerosol_interp,             &
                             trim(omob_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = omob_emission_name,        &
                             vert_interp=(/INTERP_WEIGHTED_P/)  )
!---------------------------------------------------------------------
!    Set time for input file base on selected time dependency.
!---------------------------------------------------------------------
      if (trim(omna_time_dependency_type) == 'constant' ) then
        omna_time_serie_type = 1
        omna_offset = set_time(0, 0)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omna are constant in sulfate module'
        endif
!---------------------------------------------------------------------
!    a dataset entry point must be supplied when the time dependency
!    for omna is selected.
!---------------------------------------------------------------------
      else if (trim(omna_time_dependency_type) == 'time_varying') then
        omna_time_serie_type = 3
        if (omna_dataset_entry(1) == 1 .and. &
            omna_dataset_entry(2) == 1 .and. &
            omna_dataset_entry(3) == 1 .and. &
            omna_dataset_entry(4) == 0 .and. &
            omna_dataset_entry(5) == 0 .and. &
            omna_dataset_entry(6) == 0 ) then
          omna_entry = model_init_time
        else
!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omna_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
          omna_entry  = set_date (omna_dataset_entry(1), &
                                  omna_dataset_entry(2), &
                                  omna_dataset_entry(3), &
                                  omna_dataset_entry(4), &
                                  omna_dataset_entry(5), &
                                  omna_dataset_entry(6))
        endif
        call print_date (omna_entry , str= &
          'Data from omna timeseries at time:')
        call print_date (model_init_time , str= &
          'This data is mapped to model time:')
        omna_offset = omna_entry - model_init_time
        if (model_init_time > omna_entry) then
          omna_negative_offset = .true.
        else
          omna_negative_offset = .false.
        endif
      else if (trim(omna_time_dependency_type) == 'fixed_year') then
        omna_time_serie_type = 2
        if (omna_dataset_entry(1) == 1 .and. &
            omna_dataset_entry(2) == 1 .and. &
            omna_dataset_entry(3) == 1 .and. &
            omna_dataset_entry(4) == 0 .and. &
            omna_dataset_entry(5) == 0 .and. &
            omna_dataset_entry(6) == 0 ) then
           call error_mesg ('atmos_carbon_aerosol_mod', &
            'must set omna_dataset_entry when using fixed_year source', FATAL)
        endif

!----------------------------------------------------------------------
!    define the offset from model base time (obtained from diag_table)
!    to omna_dataset_entry as a time_type variable.
!----------------------------------------------------------------------
        omna_entry  = set_date (omna_dataset_entry(1), &
                                  2,1,0,0,0)
        call error_mesg ('atmos_carbon_aerosol_mod', &
           'omna is defined from a single annual cycle &
                &- no interannual variation', NOTE)
        if (mpp_pe() == mpp_root_pe() ) then
          print *, 'omna correspond to year :', &
                    omna_dataset_entry(1)
        endif
     endif
     call interpolator_init (omna_aerosol_interp,             &
                             trim(omna_filename),           &
                             lonb, latb,                        &
                             data_out_of_bounds=  (/CONSTANT/), &
                             data_names = omna_emission_name,        &
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
      call interpolator_end ( bcff_aerosol_interp)
      call interpolator_end ( bcbb_aerosol_interp)
      call interpolator_end ( bcob_aerosol_interp)
      call interpolator_end ( omff_aerosol_interp)
      call interpolator_end ( ombb_aerosol_interp)
      call interpolator_end ( omob_aerosol_interp)
      call interpolator_end ( omna_aerosol_interp)
      module_is_initialized = .FALSE.

 end subroutine atmos_carbon_aerosol_end
!</SUBROUTINE>
!#######################################################################
end module atmos_carbon_aerosol_mod
