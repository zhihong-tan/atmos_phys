module tropchem_driver_mod
!
! <CONTACT EMAIL="Larry.Horowitz@noaa.gov">
!   Larry W. Horowitz
! </CONTACT>

! <OVERVIEW>
!     This code calculates tracer tendencies due to tropospheric chemistry
! </OVERVIEW>

! <DESCRIPTION>
!
! This code calculates chemical production and loss of tracers due
! to tropospheric chemistry. It also includes dry deposition, upper
! boundary conditions, emissions. Off-line sulfate concentrations are
! read in for use in calculating heterogeneous reaction rates (if SO4
! is not included as a tracer).
!
! This module is only activated is do_tropchem=T in tropchem_driver_nml
!
! </DESCRIPTION>


!-----------------------------------------------------------------------

use                    fms_mod, only : file_exist,   &
                                       write_version_number, &
                                       mpp_pe,  &
                                       mpp_root_pe, &
                                       lowercase,   &
                                       open_namelist_file, &
                                       close_file,   &
                                       stdlog, &
                                       check_nml_error, &
                                       error_mesg, &
                                       FATAL, &
                                       WARNING
use           time_manager_mod, only : time_type, &
                                       get_date_julian
use           diag_manager_mod, only : send_data,            &
                                       register_diag_field,  &
                                       register_static_field
use         tracer_manager_mod, only : get_tracer_index,     &
                                       get_tracer_names,     &
                                       query_method,         &
                                       check_if_prognostic
use          field_manager_mod, only : MODEL_ATMOS,          &
                                       parse
use atmos_tracer_utilities_mod, only : dry_deposition
use              constants_mod, only : grav, rdgas
use                    mpp_mod, only : mpp_clock_id,         &
                                       mpp_clock_begin,      &
                                       mpp_clock_end
use           interpolator_mod, only : interpolate_type,     &
                                       interpolator_init,    &
                                       interpolator_end,     &
                                       interpolator,         &
                                       query_interpolator,   &
                                       init_clim_diag,       &
                                       CONSTANT,             &
                                       INTERP_WEIGHTED_P  
use              mo_chemdr_mod, only : chemdr
use             mo_chemini_mod, only : chemini
use             M_TRACNAME_MOD, only : tracnam         
use                MO_GRID_MOD, only : pcnstm1 
use               MOZ_HOOK_MOD, only : moz_hook_init

implicit none

private

!-----------------------------------------------------------------------
!     ... interfaces
!-----------------------------------------------------------------------
public  tropchem_driver, tropchem_driver_init

!-----------------------------------------------------------------------
!     ...  declare type that will store the field infomation for the 
!          emission file
!-----------------------------------------------------------------------
type,public :: field_init_type
   character(len=64), pointer :: field_names(:)
end type field_init_type


!-----------------------------------------------------------------------
!     ... namelist
!-----------------------------------------------------------------------
integer, parameter :: maxinv = 100
real               :: relaxed_dt = 86400.*10.,        & ! relaxation timescale (sec) for the upper boundary values
                      ub_pres = 100.e2                  ! pressure (Pa) above which to apply chemical upper boundary conditions
character(len=64)  :: file_sulfate = 'sulfate.nc',    & ! NetCDF file for sulfate concentrations
                      file_conc = 'conc_all.nc',      & ! NetCDF file for tracer concentrations (initial and fixed)
                      file_emis_1 = 'emissions.',     & ! NetCDF file name (beginning) for emissions
                      file_emis_2 = '.nc',            & ! NetCDF file name (end) for emissions
                      file_ub = 'ub_vals.nc'            ! NetCDF file for chemical upper boundary conditions
character(len=64)  :: file_dry = 'depvel.nc',         & ! NetCDF file for dry deposition velocities
                      file_aircraft = 'aircraft.nc'     ! NetCDF file for aircraft emissions
character(len=10), dimension(maxinv) :: inv_list =''    ! list of invariant (fixed) tracers
real               :: lght_no_prd_factor = 1.           ! lightning NOx scale factor
logical            :: do_tropchem = .false.             ! Do tropospheric chemistry?
 
namelist /tropchem_driver_nml/    &
                               relaxed_dt, &
                               ub_pres, &
                               file_sulfate, &
                               file_conc, &
                               file_emis_1, &
                               file_emis_2, &
                               file_ub, &
                               file_dry, &
                               inv_list, & 
                               file_aircraft,&
                               lght_no_prd_factor, &
                               do_tropchem

character(len=7), parameter :: module_name = 'tracers'
real, parameter :: MW_air     = 28.9644,  & !molecular weightof air(g/mole)
                   Navo       = 6.023e23, & !Avogadro's constant(molec/mole)
                   g_to_kg    = 1.e-3,    & !conversion factor (kg/g)
                   m2_to_cm2  = 1.e4        !conversion factor (cm2/m2)
real, parameter :: emis_cons = MW_air * g_to_kg * m2_to_cm2 / Navo        
logical, dimension(pcnstm1) :: is_emis = .false.
type(interpolate_type),dimension(pcnstm1), save :: inter_emis
type(interpolate_type),dimension(pcnstm1), save :: inter_aircraft_emis
type(interpolate_type), save :: airc_default
type(field_init_type),dimension(pcnstm1) :: init_type
logical, dimension(pcnstm1) :: is_ub = .false.
logical, dimension(pcnstm1) :: is_airc = .false.
character(len=64),dimension(pcnstm1) :: ub_names,airc_names
real, parameter :: small = 1.e-50

!-----------------------------------------------------------------------
!     ... identification numbers for diagnostic fields
!-----------------------------------------------------------------------
integer :: id_sul, id_temp
integer :: inqa, inql, inqi !index of the three speices(nqa, nql, nqi)
logical :: module_is_initialized=.false.

integer, dimension(pcnstm1) :: indices, id_prod, id_loss, id_chem_tend, id_emiss, id_ub, id_airc

type(interpolate_type), save :: conc !to be used to read in the concentration of OH and CH4
type(interpolate_type), save :: sulfate !to be used to read in the data for sulate
type(interpolate_type), save :: ub_default !to be used for the upper bound data
type(interpolate_type),dimension(pcnstm1), save :: ub

type(interpolate_type), save :: drydep_data_default
integer :: clock_id,ndiag

!---- version number -----
character(len=128), parameter :: version     = '$Id: tropchem_driver.F90,v 13.0 2006/03/28 21:16:41 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: memphis_2006_08 $'
!-----------------------------------------------------------------------

contains


!#######################################################################

! <SUBROUTINE NAME="tropchem_driver">
!   <OVERVIEW>
!     Tropospheric chemistry driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calculates the sources and sinks of tracers
!     due to tropospheric chemistry. It is called from atmos_tracer_driver.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tropchem_driver (lon, lat, land, pwt, r, chem_dt,          &
!                           Time, phalf, pfull, t, is, js, dt,    &
!                           z_half, z_full, q, tsurf, albedo, coszen, &
!                           Time_next, rdiag, kbot)
!   </TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="land" TYPE="real" DIM="(:,:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     Pressure weighting (air mass) for each layer (kg/m2)
!   </IN>
!   <IN NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer mixing ratios (tropchem tracers in VMR)
!   </IN>
!   <IN NAME="Time, Time_next" TYPE="time_type">
!     Model time
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model half levels (Pa)
!   </IN>
!   <IN NAME="pfull" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model full levels (Pa)
!   </IN>
!   <IN NAME="t" TYPE="real" DIM="(:,:,:)">
!     Temperature.
!   </IN>
!   <IN NAME="is, js" TYPE="integer">
!     Local domain start indices
!   </IN>
!   <IN NAME="dt" TYPE="real">
!     Model physics timestep (s)
!   </IN>
!   <IN NAME="z_half" TYPE="real" DIM="(:,:,:)">
!     Height at model half levels (m)
!   </IN>
!   <IN NAME="z_full" TYPE="real" DIM="(:,:,:)">
!     Height at model full levels (m)
!   </IN>
!   <IN NAME="q" TYPE="real" DIM="(:,:,:)">
!     Specific humidity (kg/kg)
!   </IN>
!   <IN NAME="tsurf" TYPE="real" DIM="(:,:)">
!     Surface temperature (K)
!   </IN>
!   <IN NAME="albedo" TYPE="real" DIM="(:,:)">
!     Surface albedo
!   </IN>
!   <IN NAME="coszen" TYPE="real" DIM="(:,:)">
!     Cosine of the solar zenith angle
!   </IN>
!   <OUT NAME="chem_dt" TYPE="real" DIM="(:,:,:,:)">
!     Tracer tendencies from tropospheric chemistry (VMR/s)
!   </OUT>
!   <INOUT NAME="rdiag" TYPE="real" DIM="(:,:,:,:)">
!     Diagnostic tracer mixing ratios (tropchem tracers in VMR),
!     updated on output
!   </INOUT>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

subroutine tropchem_driver( lon, lat, land, pwt, r, chem_dt,          &
                            Time, phalf, pfull, t, is, js, dt,    &
                            z_half, z_full, q, tsurf, albedo, coszen, &
                            Time_next, rdiag, kbot)

!-----------------------------------------------------------------------
   real, intent(in),    dimension(:,:)            :: lon, lat
   real, intent(in),    dimension(:,:)            :: land
   real, intent(in),    dimension(:,:,:)          :: pwt
   real, intent(in),    dimension(:,:,:,:)        :: r
   real, intent(out),   dimension(:,:,:,:)        :: chem_dt
   type(time_type), intent(in)                    :: Time, Time_next     
   integer, intent(in)                            :: is,js
   real, intent(in),    dimension(:,:,:)          :: phalf,pfull,t
   real, intent(in)                               :: dt !to be passed into chemdr
   real, intent(in),    dimension(:,:,:)          :: z_half !height in meters at half levels
   real, intent(in),    dimension(:,:,:)          :: z_full !height in meters at full levels
   real, intent(in),    dimension(:,:,:)          :: q !specific humidity at current time step in kg/kg
   real, intent(in),    dimension(:,:)            :: tsurf !surface temperature
   real, intent(in),    dimension(:,:)            :: albedo
   real, intent(in),    dimension(:,:)            :: coszen
   real, intent(inout), dimension(:,:,:,:)        :: rdiag
   integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(r,1),size(r,2),size(r,3)) :: sulfate_data,o3_data
   real, dimension(size(r,1),size(r,2),size(r,3)) :: ub_temp,rno
   real, dimension(size(r,1),size(r,2),size(r,3),maxinv) :: inv_data
   real, dimension(size(r,1),size(r,2)) :: emis, temp_data
   integer :: i,j,k,n,kb,id,jd,kd,ninv,ntp
!  integer :: nno,nno2
   integer :: inv_index
   integer :: plonl
   logical :: used
   real,  dimension(size(r,1),size(r,3)) :: pdel
   real, dimension(size(r,1),size(r,2),size(r,3),pcnstm1)  ::r_temp,emis_source,r_ub,airc_emis
   real, dimension(size(land,1), size(land,2)) :: oro ! 0 and 1 rep. of land
   character(len=64) :: name
   real, dimension(size(r,1),size(r,2),size(r,3),pcnstm1) :: prod, loss
!-----------------------------------------------------------------------

!<ERROR MSG="tropchem_driver_init must be called first." STATUS="FATAL">
!   Tropchem_driver_init needs to be called before tracer_driver.
!</ERROR>
   if (.not. module_is_initialized)  &
      call error_mesg ('Tropchem_driver','tropchem_driver_init must be called first.', FATAL)

   ntp = size(r,4)
   plonl = size(r,1)

   where(land(:,:) >= 0.5)
      oro(:,:) = 1.
   elsewhere
      oro(:,:) = 0.
   endwhere

   id=size(r,1); jd=size(r,2); kd=size(r,3)
   
   ninv=0
   do n = 1, size(inv_list)
      if(inv_list(n) /= '') then
         ninv = ninv + 1
      else
         exit
      end if
   end do
 
   emis_source(:,:,:,:) = 0.0
   airc_emis(:,:,:,:) = 0.0

   do n = 1, pcnstm1
!-----------------------------------------------------------------------
!     ... read in the surface emissions, using interpolator
!-----------------------------------------------------------------------
      if (is_emis(n)) then
         call read_2D_emis_data( inter_emis(n), emis, Time, &
                                 init_type(n)%field_names, &
                                 is, js, id_emiss(n) )
      
         if (present(kbot)) then
            do j=1,jd
               do i=1,id
                  kb=kbot(i,j)
                  emis_source(i,j,kb,n) = emis(i,j)/pwt(i,j,kb) &
                                        * emis_cons
                  
               end do
            end do
         else
            emis_source(:,:,kd,n) = emis(:,:)/pwt(:,:,kd) &
                                  * emis_cons
         end if
      end if

!-----------------------------------------------------------------------
!     ... read in the aircraft emissions
!-----------------------------------------------------------------------
      if(is_airc(n)) then
         call interpolator( inter_aircraft_emis(n), Time, phalf, &
                            airc_emis(:,:,:,n), trim(airc_names(n)),is,js)
         if(id_airc(n) > 0)&
              used = send_data(id_airc(n),airc_emis(:,:,:,n),Time, is_in=is, js_in=js)
    
         airc_emis(:,:,:,n) = airc_emis(:,:,:,n)/pwt(:,:,:)*emis_cons
      end if
   end do

!-----------------------------------------------------------------------
!     ... read in the concentrations of "invariant" (i.e., prescribed)
!         species
!-----------------------------------------------------------------------
   do n = 1,ninv
      call interpolator( conc, Time, phalf, inv_data(:,:,:,n), &
                         trim(inv_list(n)), is, js)
      inv_index = get_tracer_index( MODEL_ATMOS, trim(inv_list(n)) ) - ntp
      rdiag(:,:,:,inv_index) = inv_data(:,:,:,n)
   end do
      
!-----------------------------------------------------------------------
!     ... read in the sulfate aerosol concentrations
!-----------------------------------------------------------------------
   call interpolator(sulfate, Time, phalf, sulfate_data, 'sulfate', is,js)
   used = send_data(id_sul, sulfate_data, Time, is_in=is, js_in=js)

!  call mpp_clock_begin(clock_id)

   chem_dt(:,:,:,:) =0.
    
!-----------------------------------------------------------------------
!     ... assign concentrations of prognostic (r) and diagnostic (rdiag)
!         species to r_temp
!-----------------------------------------------------------------------
   do n = 1,pcnstm1
      if(indices(n) <= ntp) then
         r_temp(:,:,:,n) = r(:,:,:,indices(n))
      else
         r_temp(:,:,:,n) = rdiag(:,:,:,indices(n)-ntp)
      end if
   end do
    
   r_temp(:,:,:,:) = MAX(r_temp(:,:,:,:),small)
  
   do j = 1,jd
      do k = 1,kd
         pdel(:,k) = phalf(:,j,k+1) - phalf(:,j,k)
      end do
      
!-----------------------------------------------------------------------
!     ... call chemistry driver
!-----------------------------------------------------------------------
      call chemdr(r_temp(:,j,:,:),             & ! species volume mixing ratios (VMR)
!                 0,                           & ! time step index
                  Time_next,                   & ! time
                  lat(:,j),                    & ! latitude
                  lon(:,j),                    & ! longitude
                  dt,                          & ! timestep in seconds
                  maxval(phalf(:,j,:),dim=2),  & ! surface press ( pascals )  
                  pfull(:,j,:),                & ! midpoint press ( pascals )
                  pdel,                        & ! delta press across midpoints
!                 oro(:,j),                    & ! surface orography flag
!                 tsurf(:,j),                  & ! surface temperature
                  z_full(:,j,:),               & ! height at midpoints ( m )
                  z_half(:,j,:),               & ! height at interfaces ( m )
                  MAX(r(:,j,:,inqa),0.),       & ! cloud fraction
                  MAX(r(:,j,:,inql)+r(:,j,:,inqi),0.), & ! total cloud water (kg/kg)
                  t(:,j,:),                    & ! temperature
                  inv_data(:,j,:,:),           & ! invariant species
                  q(:,j,:),                    & ! specific humidity ( kg/kg )
                  albedo(:,j),                 & ! surface albedo
                  coszen(:,j),                 & ! cosine of solar zenith angle
                  prod(:,j,:,:),               & ! chemical production rate
                  loss(:,j,:,:),               & ! chemical loss rate
                  sulfate_data(:,j,:),         & ! sulfate aerosol
                  plonl )                        ! number of longitudes
      
   end do

   r_temp(:,:,:,:) = MAX( r_temp(:,:,:,:), small )

   
!-----------------------------------------------------------------------
!     ... output diagnostics
!-----------------------------------------------------------------------
   do n = 1,pcnstm1
      if(id_prod(n)>0) then
         used = send_data(id_prod(n),prod(:,:,:,n),Time,is_in=is,js_in=js)
      end if
      if(id_loss(n)>0) then
         used = send_data(id_loss(n),loss(:,:,:,n),Time,is_in=is,js_in=js)
      end if
      
      if(indices(n) <= ntp) then
!-----------------------------------------------------------------------
!     ... prognostic species
!-----------------------------------------------------------------------
         if(id_chem_tend(n)>0) then
            used = send_data( id_chem_tend(n), &
                              ( r_temp(:,:,:,n) - MAX(r(:,:,:,indices(n)),0.) )/dt,&
                              Time, is_in=is,js_in=js)
         end if
!-----------------------------------------------------------------------
!     ... compute tendency
!-----------------------------------------------------------------------
         chem_dt(:,:,:,indices(n)) = airc_emis(:,:,:,n) + emis_source(:,:,:,n) + &
                                     ( r_temp(:,:,:,n) - MAX(r(:,:,:,indices(n)),0.) ) / dt

      else
!-----------------------------------------------------------------------
!     ... diagnostic species
!-----------------------------------------------------------------------
         if(id_chem_tend(n)>0) then
            used = send_data( id_chem_tend(n), &
                              ( r_temp(:,:,:,n) - MAX(rdiag(:,:,:,indices(n)-ntp),0.) )/dt,&
                              Time, is_in=is,js_in=js)
         end if
         rdiag(:,:,:,indices(n)-ntp) = r_temp(:,:,:,n)
      end if

     
!-----------------------------------------------------------------------
!     ... apply upper boundary condition
!-----------------------------------------------------------------------
      if(is_ub(n)) then
         call interpolator(ub(n), Time, phalf, r_ub(:,:,:,n), trim(ub_names(n)), is, js)
         if(id_ub(n)>0) then
            used = send_data(id_ub(n), r_ub(:,:,:,n), Time, is_in=is, js_in=js)
         end if
         where (pfull(:,:,:) < ub_pres)            
            chem_dt(:,:,:,indices(n)) = (r_ub(:,:,:,n) - r(:,:,:,indices(n))) / relaxed_dt
         endwhere
      end if
      
   end do
   
!-----------------------------------------------------------------------
!     ... special case(nox = no + no2)
!-----------------------------------------------------------------------
!  nno = get_tracer_index(MODEL_ATMOS,'no')
!  nno2 = get_tracer_index(MODEL_ATMOS,'no2')
!  if((nno /= 0) .and. (nno2 /= 0)) then
!     rno(:,:,:) = r(:,:,:,nno)/ MAX((r(:,:,:,nno) + r(:,:,:,nno2)),1.e-30)
     
!     call interpolator(ub, Time,phalf,ub_temp,'nox',is,js)

!     where(pfull(:,:,:) < ub_pres)
!        chem_dt(:,:,:,nno) =((rno(:,:,:)*ub_temp(:,:,:))-r(:,:,:,nno)) / relaxed_dt
!        chem_dt(:,:,:,nno2) = (((1.-rno(:,:,:))*ub_temp(:,:,:))-r(:,:,:,nno2)) / &
!             relaxed_dt
!     endwhere
!  end if

!  call mpp_clock_end(clock_id)
   
!-----------------------------------------------------------------------
    
end subroutine tropchem_driver
!</SUBROUTINE>

!#######################################################################

! <FUNCTION NAME="tropchem_driver_init">
!   <OVERVIEW>
!     Initializes the tropospheric chemistry driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the tropospheric chemistry module.
!     It is called from atmos_tracer_driver_init.
!     Data sets are read in for dry deposition, upper boundary conditions,
!     and emissions. Off-line sulfate concentrations are also read in for
!     use in calculating heterogeneous reaction rates (if SO4 is not
!     included as a tracer).
!   </DESCRIPTION>
!   <TEMPLATE>
!     Ltropchem = tropchem_driver_init( r, mask, axes, Time, &
!                                       lonb_mod, latb_mod, phalf, &
!                                       drydep_data )
!   </TEMPLATE>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask that designates which grid points
!      are above (1) or below (0) the ground
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="lonb_mod" TYPE="real" DIM="(:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="latb_mod" TYPE="real" DIM="(:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model half levels (Pa)
!   </IN>
!   <OUT NAME="drydep_data" TYPE="interpolate_type" DIM="(:)">
!     Tracer dry deposition velocities
!   </OUT>
!   <OUT NAME="Ltropchem" TYPE="logical">
!     Do tropospheric chemistry? (Output as function value)
!   </OUT>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer mixing ratios (tropchem tracers in VMR)
!   </INOUT>

function tropchem_driver_init( r, mask, axes, Time, &
                               lonb_mod, latb_mod, phalf, &
                               drydep_data ) result(Ltropchem)

!-----------------------------------------------------------------------
!
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
   real, intent(inout), dimension(:,:,:,:) :: r
   real, intent(in),    dimension(:,:,:), optional :: mask
   type(time_type), intent(in) :: Time
   integer        , intent(in) :: axes(4)
   real, intent(in), dimension(:) :: lonb_mod
   real, intent(in), dimension(:) :: latb_mod
   real, intent(in),dimension(:,:,:) :: phalf
   type(interpolate_type), intent(out) :: drydep_data(:)

   logical :: Ltropchem
   integer :: flag_file,flag_spec
   integer :: n,i,j,nfields,nt
   integer :: ierr, io
   character(len=64) :: nc_file,filename,specname
   character(len=64) :: units
   character(len=64) :: control=''
   character(len=64) :: name=''
   type(interpolate_type) :: init_conc
   character(len=64),dimension(pcnstm1) :: emis_files=''
   character(len=64),dimension(pcnstm1) :: conc_files=''
   character(len=64),dimension(pcnstm1) :: ub_files=''
   character(len=64),dimension(pcnstm1) :: dry_files
   character(len=64),dimension(pcnstm1) :: wet_ind, conc_names, dry_names, airc_files

   integer :: unit
   character(len=16) ::  fld

!
!-----------------------------------------------------------------------
!
   if (module_is_initialized) return

!-----------------------------------------------------------------------
!     ... write version number
!-----------------------------------------------------------------------
   call write_version_number(version, tagname)
    
!-----------------------------------------------------------------------
!     ... read namelist
!-----------------------------------------------------------------------
   if(file_exist('input.nml')) then
      unit = open_namelist_file('input.nml')
      ierr=1; do while (ierr /= 0)
      read(unit, nml = tropchem_driver_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'tropchem_driver_nml')
      end do
10    call close_file(unit)
   end if
  
   if(mpp_pe() == mpp_root_pe()) then       
      write(stdlog(), nml=tropchem_driver_nml)
   end if

   Ltropchem = do_tropchem
   if (.not. Ltropchem) then
      return
   end if
     
!-----------------------------------------------------------------------
!     ... Setup sulfate input/interpolation
!-----------------------------------------------------------------------
   call interpolator_init( sulfate, trim(file_sulfate), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/),      &
                           vert_interp=(/INTERP_WEIGHTED_P/) )

!-----------------------------------------------------------------------
!     ... Initialize chemistry driver
!-----------------------------------------------------------------------
   call CHEMINI(0.)
   
!-----------------------------------------------------------------------
!     ... set initial value of indices
!-----------------------------------------------------------------------
   indices(:) = 0
   do i=1,pcnstm1
      n = get_tracer_index(MODEL_ATMOS, tracnam(i))
      if (n >0) then
         indices(i) = n
         if (indices(i) > 0 .and. mpp_pe() == mpp_root_pe()) then 
            write(*,30) tracnam(i),indices(i)
            write(stdlog(),30) trim(tracnam(i)),indices(i)
         end if
      else
!<ERROR MSG="Tropospheric chemistry tracer not found in field table" STATUS="WARNING">
!   A tropospheric chemistry tracer was not included in the field table
!</ERROR>
         call error_mesg ('tropchem_driver_init', trim(tracnam(i)) // ' is not found', WARNING)
      end if
   end do
30 format (A,' was initialized as tracer number ',i2)

!-----------------------------------------------------------------------
!     ... Setup dry deposition
!-----------------------------------------------------------------------
   call tropchem_drydep_init( dry_files, dry_names, &
                              lonb_mod, latb_mod, &
                              drydep_data )

!-----------------------------------------------------------------------
!     ... Setup upper boundary condition data
!-----------------------------------------------------------------------
   call interpolator_init( ub_default, trim(file_ub), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/),          &
                           vert_interp=(/INTERP_WEIGHTED_P/) )

!-----------------------------------------------------------------------
!     ... Set up concentration input/interpolation
!-----------------------------------------------------------------------
   call interpolator_init( conc, trim(file_conc), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/),&
                           vert_interp=(/INTERP_WEIGHTED_P/) )

!-----------------------------------------------------------------------
!     ... Set up aircraft emissions interpolation
!-----------------------------------------------------------------------
   call interpolator_init( airc_default, trim(file_aircraft), lonb_mod, latb_mod,&
                           data_out_of_bounds=(/CONSTANT/), &
                           vert_interp=(/INTERP_WEIGHTED_P/))

!-----------------------------------------------------------------------
!     ... Setup emissions input/interpolation
!-----------------------------------------------------------------------
   do i = 1,pcnstm1
      nc_file = trim(file_emis_1)//lowercase(trim(tracnam(i)))//trim(file_emis_2)
      call init_2D_emis_data(inter_emis(i),MODEL_ATMOS,'emissions',indices(i),nc_file,&
              lonb_mod,latb_mod,init_type(i),is_emis(i))
      if( is_emis(i) ) emis_files(i) = trim(nc_file)
        
!-----------------------------------------------------------------------
!     ... Upper boundary condition
!-----------------------------------------------------------------------
      if( query_method('upper_bound', MODEL_ATMOS,indices(i),name,control) ) then
         if( trim(name)=='file' ) then
            flag_file = parse(control, 'file',filename)
            flag_spec = parse(control, 'name',specname)

            if( flag_file > 0 .and. trim(filename) /= trim(file_ub) ) then
               ub_files(i) = trim(filename)
               call interpolator_init(ub(i), trim(filename), lonb_mod, latb_mod, &
                       data_out_of_bounds=(/CONSTANT/),          &
                       vert_interp=(/INTERP_WEIGHTED_P/))
            else
               ub_files(i) = trim(file_ub)
               ub(i) = ub_default
            end if
            if(flag_spec > 0) then
               ub_names(i) = trim(specname)
            else
               ub_names(i) = trim(lowercase(tracnam(i)))
            end if

            is_ub(i) = .true.
              
         end if
      end if

!-----------------------------------------------------------------------
!     ... Initial conditions
!-----------------------------------------------------------------------
      if(.not. file_exist('INPUT/atmos_tracers.res.nc'))then
        if( .not. file_exist('INPUT/tracer_'//trim(lowercase(tracnam(i)))//'.res') ) then
           if( query_method('init_conc',MODEL_ATMOS,indices(i),name,control) ) then
             if( trim(name) == 'file' ) then
               flag_file = parse(control, 'file',filename)
               flag_spec = parse(control, 'name',specname)
                  
               if( flag_file>0 .and. trim(filename) /= trim(file_conc) ) then
                     conc_files(i) = trim(filename)
                     call interpolator_init( init_conc,trim(filename),lonb_mod,latb_mod,&
                                             data_out_of_bounds=(/CONSTANT/), &
                                             vert_interp=(/INTERP_WEIGHTED_P/) )
               else
                     conc_files(i) = trim(file_conc)
                     init_conc = conc
               end if
                  
               if( flag_spec > 0 ) then
                  conc_names(i) = trim(lowercase(specname))
                  specname = lowercase(specname)
               else
                  conc_names(i) = trim(lowercase(tracnam(i)))
                  specname = trim(lowercase(tracnam(i)))
               end if
                  
               call interpolator(init_conc, Time, phalf,r(:,:,:,indices(i)),trim(specname))                  
             end if
           end if
        end if
      end if
             
!-----------------------------------------------------------------------
!     ... Aircraft emissions
!-----------------------------------------------------------------------
      if( query_method('aircraft_emis',MODEL_ATMOS,indices(i),name,control) ) then
         is_airc(i) = .true.
         if( trim(name) == 'file' ) then
            flag_file = parse(control,'file',filename)
            flag_spec = parse(control,'name',specname)

            if( flag_file >0 .and. trim(filename) /= trim(lowercase(file_aircraft)) ) then
               airc_files(i) = trim(filename)
               call interpolator_init( inter_aircraft_emis(i),trim(filename), lonb_mod, latb_mod, &
                                       data_out_of_bounds=(/CONSTANT/), &
                                       vert_interp=(/INTERP_WEIGHTED_P/) )
            else
               airc_files(i) = trim(file_aircraft)
               inter_aircraft_emis(i) = airc_default
            end if
               
            if( flag_spec >0 ) then
               airc_names(i) = trim(specname)
            else
               airc_names(i) = trim(lowercase(tracnam(i)))
            end if
         end if
      end if

             
!-----------------------------------------------------------------------
!     ... Wet deposition
!-----------------------------------------------------------------------
      if( query_method('wet_deposition',MODEL_ATMOS,indices(i),name,control) ) then
         wet_ind(i) = 'This species has wet deposition'
      else
         wet_ind(i) = ''
      end if
         
   end do
   
!-----------------------------------------------------------------------
!     ... Print out settings for tracer
!-----------------------------------------------------------------------
   if( mpp_pe() == mpp_root_pe() ) then
      write(stdlog(),*) '---------------------------------------------------------------------------------------'
      do i = 1,pcnstm1
         write(stdlog(),*) 'The tracname index is ',i
         write(stdlog(),*) 'The tracname is ',tracnam(i)
         if(check_if_prognostic(MODEL_ATMOS,indices(i))) then
            write(stdlog(),*) 'This is a prognostic tracer.'
         else
            write(stdlog(),*) 'This is a diagnostic tracer.'
         end if
         if(emis_files(i) /= '') then
            write(stdlog(),*)'Emissions from file: ',trim(emis_files(i))
         end if
         if(ub_files(i) /= '') then
            write(stdlog(),*)'Upper BC from file: ',trim(ub_files(i)), &
                             ', with the name of ',trim(ub_names(i))
         end if
         if(conc_files(i) /= '') then
            write(stdlog(),*)'Concentration from file: ',trim(conc_files(i)), &
                             ', with the name of ',trim(conc_names(i))
         end if
         if(dry_files(i) /= '') then
            write(stdlog(),*)'Dry deposition velocity from file: ',trim(dry_files(i)), &
                             ' with the name of '//trim(dry_names(i))
         end if
         if(wet_ind(i) /= '') then
            write(stdlog(),*) wet_ind(i)
         end if
         if(is_airc(i)) then
            write(stdlog(),*)'Aircraft emissions from file: ',trim(airc_files(i)), &
                             ' with the name of '//trim(airc_names(i))
         end if
         write(stdlog(),*) '---------------------------------------------------------------------------------------'
      end do
   end if


!-----------------------------------------------------------------------
!     ... Get the index number for the cloud variables
!-----------------------------------------------------------------------
   inqa = get_tracer_index(MODEL_ATMOS,'cld_amt')!cloud fraction
   inql = get_tracer_index(MODEL_ATMOS,'liq_wat') !cloud liquid specific humidity
   inqi = get_tracer_index(MODEL_ATMOS,'ice_wat') !cloud ice water specific humidity
      

!-----------------------------------------------------------------------
!     ... Call the chemistry hook init routine
!-----------------------------------------------------------------------
   call moz_hook_init(lght_no_prd_factor,Time,axes)

!-----------------------------------------------------------------------
!     ... Register diagnostic fields for species tendencies
!-----------------------------------------------------------------------
     
   id_sul =register_diag_field( 'tracers', 'sulfate', axes(1:3), Time, 'sulfate', 'VMR' )
   do i=1,pcnstm1
      id_chem_tend(i) = register_diag_field( 'tracers', trim(tracnam(i))//'_chem_dt', axes(1:3), &
                                             Time, trim(tracnam(i))//'_chem_dt','VMR/s' )
      id_prod(i) = register_diag_field( 'tracers', trim(tracnam(i))//'_prod', axes(1:3), &
                                        Time, trim(tracnam(i))//'_prod','VMR/s')
      id_loss(i) = register_diag_field( 'tracers', trim(tracnam(i))//'_loss', axes(1:3), &
                                        Time, trim(tracnam(i))//'_loss','VMR/s')
      if( is_emis(i) ) then
         id_emiss(i) = register_diag_field( 'tracers', trim(tracnam(i))//'_emis', axes(1:2), &
                                            Time, trim(tracnam(i))//'_emis', 'molec/cm2/s')
      else
         id_emiss(i) = 0
      end if
      if( is_ub(i) ) then
         id_ub(i) = register_diag_field( 'tracers', trim(tracnam(i))//'_up', axes(1:3), &
                                         Time, trim(tracnam(i))//'_up','VMR' )
      else
         id_ub(i) = 0
      end if
      if( is_airc(i) ) then
         id_airc(i) = register_diag_field( 'tracers', trim(tracnam(i))//'_airc_emis', axes(1:3), &
                                           Time, trim(tracnam(i))//'_airc_emis','molec/cm2/s' )
      else
         id_airc(i) = 0
      end if
   end do
    
!-----------------------------------------------------------------------
!     ... initialize mpp clock id
!-----------------------------------------------------------------------
!  clock_id = mpp_clock_id('Chemistry')
      
   module_is_initialized = .true.
      
      
!-----------------------------------------------------------------------
      
end function tropchem_driver_init
!</FUNCTION>
 
      
subroutine tropchem_driver_end

!-----------------------------------------------------------------------
!     ... initialize mpp clock id
!-----------------------------------------------------------------------
      
   module_is_initialized = .false.
      
      
!-----------------------------------------------------------------------
      
end subroutine tropchem_driver_end

!#######################################################################

! <SUBROUTINE NAME="read_2D_emis_data">
!   <OVERVIEW>
!     Read emissions file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Reads tracer surface emissions from a NetCDF file
!   </DESCRIPTION>
!   <TEMPLATE>
!     call read_2D_emis_data( emis_type, emis, Time, field_names, is, js, id_emis ) 
!   </TEMPLATE>

subroutine read_2D_emis_data( emis_type, emis, Time, field_names, is, js, id_emis )
    
   type(interpolate_type),intent(inout) :: emis_type
   real, dimension(:,:),intent(out) :: emis
   type(time_type),intent(in) :: Time
   character(len=*),dimension(:), intent(in) :: field_names
   integer, intent(in) :: is, js
   integer, intent(in),optional :: id_emis ! id for diagnostic


   integer :: k
   logical :: used
   real, dimension(size(emis,1),size(emis,2)) :: temp_data

   emis(:,:) = 0.
   temp_data(:,:) = 0.
   do k = 1,size(field_names)
      call interpolator(emis_type,Time,temp_data,field_names(k),is,js)
      emis(:,:) = emis(:,:) + temp_data(:,:)
   end do

   if(present(id_emis) .and. (id_emis > 0)) then
      used = send_data(id_emis,emis,Time,is_in=is,js_in=js)
   end if
end subroutine read_2D_emis_data
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="init_2D_emis_data">
!   <OVERVIEW>
!     Open emissions file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Opens NetCDF file of tracer surface emissions for reading,
!     and set up interpolation to model grid/time
!   </DESCRIPTION>
!   <TEMPLATE>
!     call init_2D_emis_data( emis_type, model, method_type, index, file_name, &
!                             lonb_mod, latb_mod, field_type, flag )
!   </TEMPLATE>

subroutine init_2D_emis_data( emis_type, model, method_type, index, file_name, &
                              lonb_mod, latb_mod, field_type, flag )
    
   type(interpolate_type),intent(inout) :: emis_type
   integer, intent(in) :: model,index
   character(len=*),intent(in) :: method_type
   character(len=*),intent(inout) ::file_name
   real,intent(in),dimension(:) :: lonb_mod,latb_mod
   type(field_init_type),intent(out) :: field_type
   logical, intent(out) :: flag
    
   character(len=64) :: name, control
   integer :: nfields
   integer :: flag_name, flag_file
   character(len=64) :: emis_name, emis_file

   flag = .false.
   control = ''
   if( query_method(trim(method_type),model,index,name,control) ) then
!++lwh
!     if( trim(name) == 'true' ) then
      if( trim(name) == 'file' ) then
         flag = .true.
         flag_file = parse(control, 'file', emis_file)
         flag_name = parse(control, 'name', emis_name)
         if(flag_file > 0) then
            file_name = emis_file
         else if (flag_name > 0) then
            file_name  = trim(file_emis_1)//trim(emis_name)//trim(file_emis_2)
         end if

!        if( trim(control) /='' ) then
!           if( control(len(control)-3:len(control)) /='.nc' ) then
!              file_name  = trim(file_emis_1)//trim(control)//trim(file_emis_2)
!           else
!              file_name = trim(control)
!           end if
!        end if
         call interpolator_init( emis_type, trim(file_name), &
                                 lonb_mod, latb_mod,  &
                                 data_out_of_bounds=(/CONSTANT/), &
                                 vert_interp=(/INTERP_WEIGHTED_P/) )
         call query_interpolator(emis_type,nfields=nfields)
         allocate(field_type%field_names(nfields))
         call query_interpolator(emis_type,field_names=field_type%field_names)
!--lwh
      end if
   end if
end subroutine init_2D_emis_data
!</SUBROUTINE>

!############################################################################

! <SUBROUTINE NAME="tropchem_drydep_init">
!   <OVERVIEW>
!     Open dry deposition file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Opens NetCDF file of tracer dry deposition velocities for reading,
!     and set up interpolation to model grid/time
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tropchem_drydep_init( dry_files, dry_names, &
!                                lonb_mod, latb_mod, &
!                                drydep_data )
!   </TEMPLATE>

subroutine tropchem_drydep_init( dry_files, dry_names, &
                                 lonb_mod, latb_mod, &
                                 drydep_data )

!-----------------------------------------------------------------------

   real,                   intent(in),  dimension(:) :: lonb_mod, latb_mod
   character(len=64),      intent(out), dimension(:) :: dry_files, dry_names
   type(interpolate_type), intent(out)               :: drydep_data(:)

!-----------------------------------------------------------------------

   integer :: i
   integer :: flag_file,flag_spec
   character(len=64) :: filename,specname
   character(len=64) :: name='', control=''

!-----------------------------------------------------------------------

   do i = 1,pcnstm1
      dry_files(i) = ''
      dry_names(i) = ''
      if( query_method('dry_deposition',MODEL_ATMOS,indices(i),name,control) )then
         if( trim(name) == 'file' ) then
            flag_file = parse(control, 'file',filename)
            flag_spec = parse(control, 'name',specname)
            if(flag_file > 0 .and. trim(filename) /= trim(file_dry)) then
               dry_files(i) = trim(filename)
            else
               dry_files(i) = trim(file_dry)
            end if
            if(flag_spec >0) then
               dry_names(i) = trim(specname)
            else
               dry_names(i) = trim(lowercase(tracnam(i)))
            end if
         end if
      end if
   end do

!---------- Set interpolator type for dry deposition
   call interpolator_init( drydep_data_default, trim(file_dry), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/), &
                           vert_interp=(/INTERP_WEIGHTED_P/))

   do i = 1,pcnstm1
      if( query_method('dry_deposition',MODEL_ATMOS,i,name,control) ) then
         if( trim(name) =='file' ) then
            flag_file = parse(control, 'file',filename)
            if( flag_file>0 ) then
               control = trim(filename)
            end if
            if( (trim(control)==trim(file_dry)) .or. (trim(control)=='') ) then
               drydep_data(i) = drydep_data_default
            else
               call interpolator_init( drydep_data(i), trim(control), lonb_mod, latb_mod,&
                                       data_out_of_bounds=(/CONSTANT/), &
                                       vert_interp=(/INTERP_WEIGHTED_P/))
            end if
         end if
      end if
   end do
   
end subroutine tropchem_drydep_init
!</SUBROUTINE>

!############################################################################
end module tropchem_driver_mod
