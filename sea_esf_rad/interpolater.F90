!#define FEZ
module interpolater_mod
!
! Purpose: Module to interpolate climatology data to model grid.
!
! author: William Cooke wfc@gfdl.noaa.gov
!

use mpp_mod,           only : mpp_error, &
                              FATAL,     &
                              mpp_pe,    &
                              mpp_init,  &
                              mpp_exit,  &
                              mpp_npes,  &
                              WARNING,   &
                              NOTE
use mpp_io_mod,        only : mpp_open,          &
                              mpp_close,         & 
                              mpp_get_times,     &
                              mpp_get_atts,      &
                              mpp_get_info,      &
                              mpp_read,          &
                              mpp_get_axes,      &
                              mpp_get_axis_data, &
                              mpp_get_fields,    &
                              fieldtype,         &
                              atttype,           &
                              axistype,          &
                              MPP_RDONLY,        &
                              MPP_NETCDF,        &
                              MPP_MULTI,         &
                              MPP_APPEND,        &
                              MPP_SINGLE
use mpp_domains_mod,   only : mpp_domains_init,      &
                              mpp_update_domains,    &
                              mpp_define_domains,    &
                              mpp_global_field,      &
                              domain2d,              &
                              mpp_define_layout,     &
                              mpp_get_compute_domain
use diag_manager_mod
#ifdef FEZ
use utilities_mod,     only : lowercase, write_version_number, &
                              file_exist
#else
use fms_mod,           only : lowercase, write_version_number, &
                              file_exist, mpp_root_pe, stdlog
#endif
use horiz_interp_mod,  only : horiz_interp_type, &
                              horiz_interp_init, &
                              horiz_interp,      &
                              horiz_interp_end
use time_manager_mod,  only : time_type,   &
                              set_time,    &
                              set_date,    &
                              operator(+)
#ifdef FEZ
use time_interp_mod,   only : time_interp
#else
use time_interp_mod,   only : time_interp, YEAR
#endif
use constants_mod,     only : grav

implicit none
private 

public interpolater_init, &
       interpolater,      &
       interpolater_end,  &
       init_clim_diag,    &
       read_data

character(len=128) :: version = '$Id: interpolater.F90,v 1.2 2002/07/16 22:35:17 fms Exp $'
character(len=128) :: tagname = '$Name: havana $'
logical            :: module_is_initialized = .false.
logical            :: clim_diag_initialized = .false.

type, public  :: interpolate_type
private
!Redundant data between fields
!All climatology data
real, pointer            :: lat(:)
real, pointer            :: lon(:)
real, pointer            :: latb(:)
real, pointer            :: lonb(:)
real, pointer            :: levs(:)
real, pointer            :: halflevs(:) 
type(horiz_interp_type)  :: interph
type(time_type), pointer :: time_slice(:) ! An array of the times within the climatology.
integer                  :: unit          ! Unit number on which file is being read.
character(len=64)        :: file_name     ! Climatology filename
integer                  :: TIME_FLAG     ! Linear or seaonal interpolation?
integer                  :: level_type    ! Pressure or Sigma level
integer                  :: is,ie,js,je

!Field specific data  for nfields
type(fieldtype),   pointer :: field_type(:)   ! NetCDF field type
character(len=64), pointer :: field_name(:)   ! name of this field
integer,           pointer :: time_init(:,:)  ! second index is the number of time_slices being kept. 2 or ntime.
integer,           pointer :: mr(:)           ! Flag for conversion of climatology to mixing ratio. 
integer,           pointer :: out_of_bounds(:)! Flag for when surface pressure is out of bounds.
!++lwh
integer,           pointer :: vert_interp(:)  ! Flag for type of vertical interpolation.
!--lwh
real,              pointer :: data(:,:,:,:,:) ! (nlatmod,nlonmod,nlevclim,size(time_init,2),nfields)
end type interpolate_type


integer :: ndim, nvar,natt,ntime
integer :: nlat,nlatb,nlon,nlonb,nlev,nlevh
integer :: i, j, k, len, ntime_in, n, num_fields
type(axistype), allocatable :: axes(:)
type(axistype) :: time_axis
type(fieldtype), allocatable :: varfields(:)

real, allocatable :: time_in(:)
real, allocatable :: climdata(:,:,:)

character(len=32) :: name, units

integer, parameter :: max_diag_fields = 30

!++lwh
! Flags to indicate whether the time interpolation should be linear or some other scheme for seasonal data.
integer, parameter :: LINEAR = 1, SEASONAL = 2

! Flags to indicate where climatology pressure levels are pressure or sigma levels
integer, parameter :: PRESSURE = 1, SIGMA = 2 

! Flags to indicate whether the climatology units are mixing ratio (kg/kg) or column integral (kg/m2).
! Vertical interpolation scheme requires mixing ratio at this time.
integer, parameter :: NO_CONV = 1, KG_M2 = 2 

! Flags to indicate what to do when the model surface pressure exceeds the  climatology surface pressure level.
integer, parameter, public :: CONSTANT = 1, ZERO = 2 

! Flags to indicate the type of vertical interpolation
integer, parameter, public :: INTERP_WEIGHTED_P = 10, INTERP_LINEAR_P = 20, INTERP_LOG_P = 30
!--lwh

integer :: num_clim_diag = 0
character(len=64) :: climo_diag_name(max_diag_fields)
integer :: climo_diag_id(max_diag_fields), hinterp_id(max_diag_fields)
real ::  missing_value = -1.e10
integer :: itaum, itaup

contains
!
!#######################################################################
!
subroutine interpolater_init( clim_type, file_name, lonb_mod, latb_mod, &
                              data_names, data_out_of_bounds,           &
                              vert_interp, clim_units )
type(interpolate_type), intent(out) :: clim_type
character(len=*), intent(in)            :: file_name
real            , intent(in)            :: lonb_mod(:), latb_mod(:)
character(len=*), intent(in) , optional :: data_names(:)
!++lwh
integer         , intent(in)            :: data_out_of_bounds(:) 
integer         , intent(in), optional  :: vert_interp(:) 
!--lwh
character(len=*), intent(out), optional :: clim_units(:)
!
! INTENT IN
!  file_name  :: Climatology filename
!  lonb_mod   :: The boundaries of the model grid-box longitudes.
!  latb_mod   :: The boundaries of the model grid_box latitudes.
!  data_names :: A list of the names of components within the climatology file which you wish to read.
!  data_out_of_bounds :: A list of the flags that are to be used in determining what to do if the pressure levels in the model
!                        go out of bounds from those of the climatology.
!  vert_interp:: Flag to determine type of vertical interpolation
!
! INTENT OUT
!  clim_type  :: An interpolate type containing the necessary file and field data to be passed to the interpolater routine.
!  clim_units :: A list of the units for the components listed in data_names.
!

integer                      :: unit, log_unit
character(len=64)            :: src_file, cart, namelev
!++lwh
integer                      :: num_files
real                         :: dlat, dlon
!--lwh
type(axistype), allocatable  :: axes_field(:)
type(horiz_interp_type)      :: interph
type(time_type), allocatable :: time_slice(:)
type(time_type)              :: base_time
logical                      :: NAME_PRESENT
real                         :: dtr,tpi

tpi = 4.*acos(0.)
dtr = tpi/360.

num_fields = 0

!--------------------------------------------------------------------
! open source file containing fields to be interpolated
!--------------------------------------------------------------------
src_file = 'INPUT/'//trim(file_name)
#ifdef FEZ
! Note: For fez code mpp_open adds a .nc extension to your file. Therefore you must pass
! mpp_open a bare name.
if(file_exist(trim(src_file)//'.nc')) then
#else
! For Galway code the .nc can be included in file_name.
if(file_exist(trim(src_file))) then
#endif
   call mpp_open( unit, trim(src_file), action=MPP_RDONLY, &
                  form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE )
else
!Climatology file doesn't exist, so exit
   call mpp_error(FATAL,'Interpolater_init : Data file '//trim(src_file)//' does not exist')
endif

!Find the number of variables (nvar) in this file
call mpp_get_info(unit, ndim, nvar, natt, ntime)
clim_type%unit      = unit
clim_type%file_name = trim(file_name)

num_fields = nvar
if(present(data_names)) num_fields= size(data_names)

! -------------------------------------------------------------------
! Allocate space for the number of axes in the data file.
! -------------------------------------------------------------------
allocate(axes(ndim))
call mpp_get_axes(unit, axes, time_axis)

nlon=0 ! Number of longitudes (center-points) in the climatology.
nlat=0 ! Number of latitudes (center-points) in the climatology.
nlev=0 ! Number of levels (center-points) in the climatology.
nlatb=0 ! Number of longitudes (boundaries) in the climatology.
nlonb=0 ! Number of latitudes (boundaries) in the climatology.
nlevh=0 ! Number of levels (boundaries) in the climatology.

clim_type%level_type = 0 ! Default value
do i = 1, ndim
 call mpp_get_atts(axes(i), name=name,len=len,units=units)
select case(name)
  case('lat')
    nlat=len
    allocate(clim_type%lat(nlat))
    call mpp_get_axis_data(axes(i),clim_type%lat)
    if(lowercase(units(1:6)) /= "degree") &
       call mpp_error(FATAL, "interpolater_init : Units for lat are not degrees")
  case('lon')
    nlon=len
    allocate(clim_type%lon(nlon))
    call mpp_get_axis_data(axes(i),clim_type%lon)
    if(lowercase(units(1:6)) /= "degree") &
       call mpp_error(FATAL, "interpolater_init : Units for lon are not degrees")
  case('latb')
    nlatb=len
    allocate(clim_type%latb(nlatb))
    call mpp_get_axis_data(axes(i),clim_type%latb)
    if(lowercase(units(1:6)) /= "degree") &
       call mpp_error(FATAL, "interpolater_init : Units for latb are not degrees")
  case('lonb')
    nlonb=len
    allocate(clim_type%lonb(nlonb))
    call mpp_get_axis_data(axes(i),clim_type%lonb)
    if(lowercase(units(1:6)) /= "degree") &
       call mpp_error(FATAL, "interpolater_init : Units for lonb are not degrees")
  case('pfull')
    nlev=len
    allocate(clim_type%levs(nlev))
    call mpp_get_axis_data(axes(i),clim_type%levs)
    clim_type%level_type = PRESSURE
! Convert to Pa
    if( chomp(units) == "mb" .or. chomp(units) == "hPa") then
       clim_type%levs = clim_type%levs * 100.
    end if
  case('phalf')
    nlevh=len
    allocate(clim_type%halflevs(nlevh))
    call mpp_get_axis_data(axes(i),clim_type%halflevs)
    clim_type%level_type = PRESSURE
! Convert to Pa
    if( chomp(units) == "mb" .or. chomp(units) == "hPa") then
       clim_type%halflevs = clim_type%halflevs * 100.
    end if
  case('sigma_full')
    nlev=len
    allocate(clim_type%levs(nlev))
    call mpp_get_axis_data(axes(i),clim_type%levs)
    clim_type%level_type = SIGMA
  case('sigma_half')
    nlevh=len
    allocate(clim_type%halflevs(nlevh))
    call mpp_get_axis_data(axes(i),clim_type%halflevs)
    clim_type%level_type = SIGMA
 end select

enddo

deallocate(axes)


! In the case where only the midpoints of the longitudes are defined we force the definition
! of the boundaries to be half-way between the midpoints.
if (.not. associated(clim_type%lon) .and. .not. associated(clim_type%lonb)) &
   call mpp_error(FATAL,'Interpolater_init : There appears to be no longitude axis in file '//file_name)
if (.not. associated(clim_type%lonb) ) then
allocate(clim_type%lonb(size(clim_type%lon)+1))
dlon = (clim_type%lon(2)-clim_type%lon(1))/2.0
clim_type%lonb(1) = clim_type%lon(1) - dlon
clim_type%lonb(2:) = clim_type%lon(1:) + dlon
endif

clim_type%lonb=clim_type%lonb*dtr 
! This assumes the lonb are in degrees in the NetCDF file!

if (.not. associated(clim_type%lat) .and. .not. associated(clim_type%latb)) &
   call mpp_error(FATAL,'Interpolater_init : There appears to be no latitude axis in file '//file_name)
! In the case where only the grid midpoints of the latitudes are defined we force the 
! definition of the boundaries to be half-way between the midpoints.
if (.not. associated(clim_type%latb) ) then
!++lwh
   allocate(clim_type%latb(nlat+1))
   dlat = (clim_type%lat(2)-clim_type%lat(1)) * 0.5
   clim_type%latb(1) = min( 90., max(-90., clim_type%lat(1) - dlat) )
   clim_type%latb(2:nlat) = ( clim_type%lat(1:nlat-1) + clim_type%lat(2:nlat) ) * 0.5
   dlat = ( clim_type%lat(nlat) - clim_type%lat(nlat-1) ) * 0.5
   clim_type%latb(nlat+1) = min( 90., max(-90., clim_type%lat(nlat) - dlat) )
!--lwh
endif
clim_type%latb=clim_type%latb*dtr

!Assume that the horizontal interpolation within a file is the same for each variable.

 call horiz_interp_init(clim_type%interph, &
                        clim_type%lonb, clim_type%latb, &
                        lonb_mod, latb_mod)

! Now get the time data for the file. 
! Need to figure this out from the time axis attributes. For now make it 1/1/1980
base_time = set_date(1980,1,1,0,0,0)
ntime_in = 1
if (ntime > 0) then
   allocate(time_in(ntime), clim_type%time_slice(ntime))
   call mpp_get_times(clim_type%unit, time_in)
   ntime_in = ntime
   do i = 1, ntime
     !Assume that the times in the data file correspond to days only.
     clim_type%time_slice(i) = set_time(0,INT(time_in(i)))! + base_time
   enddo
else
   allocate(time_in(1), clim_type%time_slice(1))
   time_in = 0.0
   clim_type%time_slice = set_time(0,0) + base_time
endif
deallocate(time_in)

select case(ntime)
 case (5:)
! We have more than 4 timelevels 
! Assume we have monthly or higher time resolution datasets (climatology or time series)
! So we only need to read 2 datasets and apply linear temporal interpolation.
   allocate(clim_type%data(size(lonb_mod)-1, size(latb_mod)-1, nlev, 2, num_fields))
! This will be the climatology data horizontally interpolated, so it will be on the model 
! horizontal grid, but it will still be on the climatology vertical grid.
   clim_type%TIME_FLAG = LINEAR
 case (1:4) 
! Assume we have seasonal data and read in all the data.
! We can apply sine curves to these data.
 
   allocate(clim_type%data(size(lonb_mod)-1, size(latb_mod)-1, nlev, ntime, num_fields))
   clim_type%TIME_FLAG = SEASONAL
! case (default)
   
end select


!-----------------------------------------------------------------------------
!Allocate space for the single time level of the climatology on its grid size.
!-----------------------------------------------------------------------------

   if(clim_type%TIME_FLAG .eq. LINEAR ) then
   allocate(clim_type%time_init(num_fields,2))
   else
   allocate(clim_type%time_init(num_fields,ntime))
   endif
   clim_type%time_init(:,:) = 0

allocate(clim_type%field_name(num_fields))
allocate(clim_type%field_type(num_fields))
allocate(clim_type%mr(num_fields))
allocate(clim_type%out_of_bounds(num_fields))
clim_type%out_of_bounds(:)=0
allocate(clim_type%vert_interp(num_fields))
clim_type%vert_interp(:)=0
!--------------------------------------------------------------------
!Allocate the space for the fields within the climatology data file.
allocate(varfields(nvar))
!--------------------------------------------------------------------
! Get the variable names out of the file.
call mpp_get_fields(clim_type%unit, varfields)

if(present(data_names)) then

!++lwh
   if ( size(data_out_of_bounds) /= size(data_names) .and. size(data_out_of_bounds) /= 1 ) &
      call mpp_error(FATAL,'interpolater_init : The size of the data_out_of_bounds array must be 1&
                            & or size(data_names)')
   if (present(vert_interp) .and. &
       size(vert_interp) /= size(data_names) .and. size(vert_interp) /= 1 ) &
      call mpp_error(FATAL,'interpolater_init : The size of the vert_interp array must be 1&
                            & or size(data_names)')
! Only read the fields named in data_names
   do j=1,size(data_names)
      NAME_PRESENT = .FALSE.
      do i=1,nvar
         call mpp_get_atts(varfields(i),name=name,ndim=ndim,units=units)
         if( name == data_names(j) ) then
            units=chomp(units)
            if (mpp_pe() == 0 ) write(*,*) 'Initializing src field : ',trim(name)
            clim_type%field_name(j) = name
            clim_type%field_type(j) = varfields(i)
            clim_type%mr(j)         = check_climo_units(units)
            NAME_PRESENT = .TRUE.
            if (present(clim_units)) clim_units(j) = units
            clim_type%out_of_bounds(j) = data_out_of_bounds( MIN(j,SIZE(data_out_of_bounds)) )
            if( clim_type%out_of_bounds(j) /= CONSTANT .and. &
                clim_type%out_of_bounds(j) /= ZERO ) &
               call mpp_error(FATAL,"Interpolater_init: data_out_of_bounds must be&
                                    & set to ZERO or CONSTANT")               
            if( present(vert_interp) ) then
               clim_type%vert_interp(j) = vert_interp( MIN(j,SIZE(vert_interp)) )
               if( clim_type%vert_interp(j) /= INTERP_WEIGHTED_P .and. &
                   clim_type%vert_interp(j) /= INTERP_LINEAR_P ) &
                  call mpp_error(FATAL,"Interpolater_init: vert_interp must be&
                                       & set to INTERP_WEIGHTED_P or INTERP_LINEAR_P")
            else
               clim_type%vert_interp(j) = INTERP_WEIGHTED_P
            end if
         endif
      enddo
      if(.not. NAME_PRESENT) &
         call mpp_error(FATAL,'interpolater_init : Check names of fields being passed. ' &
                              //trim(data_names(j))//' does not exist.')
   enddo
else

   if ( size(data_out_of_bounds) /= nvar .and. size(data_out_of_bounds) /= 1 ) &
      call mpp_error(FATAL,'interpolater_init : The size of the out of bounds array must be 1&
                           & or the number of fields in the climatology dataset')
   if ( present(vert_interp) .and. size(vert_interp) /= nvar .and. size(vert_interp) /= 1 ) &
      call mpp_error(FATAL,'interpolater_init : The size of the vert_interp array must be 1&
                           & or the number of fields in the climatology dataset')
! Read all the fields within the climatology data file.
   do i=1,nvar
      call mpp_get_atts(varfields(i),name=name,ndim=ndim,units=units)
         if (mpp_pe() ==0 ) write(*,*) 'Initializing src field : ',trim(name)
         clim_type%field_name(i) = lowercase(trim(name))
         clim_type%field_type(i) = varfields(i)
         clim_type%mr(i)         = check_climo_units(units)
         if (present(clim_units)) clim_units(i) = units
         clim_type%out_of_bounds(i) = data_out_of_bounds( MIN(i,SIZE(data_out_of_bounds)) )
         if( clim_type%out_of_bounds(i) /= CONSTANT .and. &
             clim_type%out_of_bounds(i) /= ZERO ) &
            call mpp_error(FATAL,"Interpolater_init: data_out_of_bounds must be&
                                 & set to ZERO or CONSTANT")
         if( present(vert_interp) ) then
            clim_type%vert_interp(i) = vert_interp( MIN(i,SIZE(vert_interp)) )
            if( clim_type%vert_interp(i) /= INTERP_WEIGHTED_P .and. &
                clim_type%vert_interp(i) /= INTERP_LINEAR_P ) &
               call mpp_error(FATAL,"Interpolater_init: vert_interp must be&
                                    & set to INTERP_WEIGHTED_P or INTERP_LINEAR_P")
         else
            clim_type%vert_interp(i) = INTERP_WEIGHTED_P
         end if
   end do
!--lwh
endif

deallocate(varfields)


if( clim_type%TIME_FLAG .eq. SEASONAL ) then
! Read all the data at this point.
   do i=1,num_fields
      do n = 1, ntime
         call read_data( clim_type, clim_type%field_type(i), &
                         clim_type%data(:,:,:,n,i), n, i, base_time )
      enddo
   enddo
endif


module_is_initialized = .true.

#ifdef FEZ
call mpp_open(log_unit,file='logfile.out', action=MPP_APPEND)
if ( mpp_pe() == 0 ) then
   write (log_unit,'(/,80("="),/(a))') trim(version), trim(tagname)
endif
call mpp_close(log_unit)
#else
call write_version_number (version, tagname)
#endif

end subroutine interpolater_init
!
!#######################################################################
!
function check_climo_units(units)
! Function to check the units that the climatology data is using. 
! This is needed to allow for conversion of datasets to mixing ratios which is what the 
! vertical interpolation scheme requires
! The default is to assume no conversion is needed.
! If the units are those of a column burden (kg/m2) then conversion to mixing ratio is flagged.
!
character(len=*), intent(in) :: units

integer :: check_climo_units

check_climo_units = NO_CONV
select case(chomp(units))
  case('kg/m2')
     check_climo_units = KG_M2
  case('kg/m^2')
     check_climo_units = KG_M2
  case('kg/m**2')
     check_climo_units = KG_M2
  case('kg m^-2')
     check_climo_units = KG_M2
  case('kg m**-2')
     check_climo_units = KG_M2  
end select

end function check_climo_units
!
!#######################################################################
!
subroutine init_clim_diag(clim_type, mod_axes, init_time)
!
! Routine to register diagnostic fields for the climatology file. 
! This routine calculates the domain decompostion of the climatology fields 
! for later export through send_data.
! The ids created here are for column burdens that will diagnose the vertical interpolation routine.
! climo_diag_id : 'module_name = climo' is intended for use with the model vertical resolution.
! hinterp_id    : 'module_name = 'hinterp' is intended for use with the climatology vertical resolution.

! INTENT INOUT :
!    clim_type : The interpolate type containing the names of the fields in the climatology file.
!
! INTENT IN    :
!   mod_axes   : The axes of the model.
!   init_time  : The model initialization time.
!
type(interpolate_type), intent(inout)  :: clim_type
integer               , intent(in)     :: mod_axes(:)
type(time_type)       , intent(in)     :: init_time

integer :: axes(2),nxd,nyd,ndivs,i
type(domain2d) :: domain
integer :: domain_layout(2), iscomp, iecomp,jscomp,jecomp


if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
   call mpp_error(FATAL, "init_clim_diag : You must call interpolater_init before calling init_clim_diag")


ndivs = mpp_npes()
nxd = size(clim_type%lon)
nyd = size(clim_type%lat)

! Define the domain decomposition of the climatology file. This may be (probably is) different from the model domain.
call mpp_define_layout ((/1,nxd,1,nyd/), ndivs, domain_layout)
call mpp_define_domains((/1,nxd,1,nyd/),domain_layout, domain,xhalo=0,yhalo=0)  
call mpp_get_compute_domain (domain, iscomp, iecomp, jscomp, jecomp)
   axes(1) = diag_axis_init(clim_type%file_name(1:5)//'x',clim_type%lon,units='degrees',cart_name='x',domain2=domain)
   axes(2) = diag_axis_init(clim_type%file_name(1:5)//'y',clim_type%lat,units='degrees',cart_name='y',domain2=domain)
clim_type%is = iscomp
clim_type%ie = iecomp
clim_type%js = jscomp
clim_type%je = jecomp

!init_time = set_date(1980,1,1,0,0,0)

if ((num_clim_diag + size(clim_type%field_name)) .gt. max_diag_fields )  &
   call mpp_error(FATAL, "init_clim_diag : Trying to set up too many diagnostic fields for the climatology data")
do i=1,size(clim_type%field_name)
climo_diag_name(i+num_clim_diag) = clim_type%field_name(i)
climo_diag_id(i+num_clim_diag) =  register_diag_field('climo',clim_type%field_name(i),axes(1:2),init_time,&
                                'climo_'//clim_type%field_name(i), 'kg/kg', missing_value)
hinterp_id(i+num_clim_diag) =  register_diag_field('hinterp',clim_type%field_name(i),mod_axes(1:2),init_time,&
                                'interp_'//clim_type%field_name(i),'kg/kg' , missing_value)
enddo
! Total number of climatology diagnostics (num_clim_diag). This can be from multiple climatology fields with different spatial axes. 
! It is simply a holder for the diagnostic indices.
num_clim_diag = num_clim_diag+size(clim_type%field_name)

clim_diag_initialized = .true.

end subroutine init_clim_diag
!
!#######################################################################
!
subroutine interpolater(clim_type, Time, phalf, interp_data,field_name, is,js, clim_units)
!
! INTENT INOUT
!   clim_type   : The interpolate type previously defined by a call to interpolater_init
!
! INTENT IN
!   field_name  : The name of the field that you wish to interpolate.
!   Time        : The model time that you wish to interpolate to.
!   phalf       : The half level model pressure field.
!   is, js      : The indices of the physics window.
!
! INTENT OUT
!   interp_data : The model field with the interpolated climatology data.
!   clim_units  : The units of field_name
!
type(interpolate_type), intent(inout)  :: clim_type
character(len=*)      , intent(in)  :: field_name
type(time_type)       , intent(in)  :: Time
real, dimension(:,:,:), intent(in)  :: phalf
real, dimension(:,:,:), intent(out) :: interp_data
integer               , intent(in) , optional :: is,js
character(len=*)      , intent(out), optional :: clim_units
real :: tweight
integer :: taum, taup, ilon,year1,year2
real :: hinterp_data(size(interp_data,1),size(interp_data,2),size(clim_type%levs))
real :: p_fact(size(interp_data,1),size(interp_data,2))
real :: col_data(size(interp_data,1),size(interp_data,2))
real :: pclim(size(clim_type%halflevs))
integer :: istart,iend,jstart,jend
logical :: result, found


if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
   call mpp_error(FATAL, "interpolater : You must call interpolater_init before calling interpolater")

istart = 1
if (present(is)) istart = is
iend = istart - 1 + size(interp_data,1)

jstart = 1
if (present(js)) jstart = js
jend = jstart - 1 + size(interp_data,2)

do i= 1,size(clim_type%field_name)
!++lwh
  if ( field_name == clim_type%field_name(i) ) then
!--lwh

    if(present(clim_units)) then
      call mpp_get_atts(clim_type%field_type(i),units=clim_units)
      clim_units = chomp(clim_units)
    endif
    if(size(clim_type%time_slice).le. 12 ) then
#ifdef FEZ
      call time_interp(Time, tweight, year1, year2, taum, taup )
#else
       call time_interp(Time, clim_type%time_slice, tweight, taum, taup, modtime=YEAR )
#endif
    else
#ifndef FEZ
       call time_interp(Time, clim_type%time_slice, tweight, taum, taup )
#endif
    endif

    if(clim_type%TIME_FLAG .ne. LINEAR ) then
      itaum=taum
      itaup=taup
    endif


    if(clim_type%TIME_FLAG .eq. LINEAR ) then
! We need 2 time levels. Check we have the correct data.
      itaum=0
      itaup=0
      do n=1,size(clim_type%time_init,2)
        if (clim_type%time_init(i,n) .eq. taum ) itaum = n
        if (clim_type%time_init(i,n) .eq. taup ) itaup = n
      enddo

      if (itaum.eq.0 .and. itaup.eq.0) then
!Neither time is set so we need to read 2 time slices.
!Set up 
! field(:,:,:,1) as the previous time slice.
! field(:,:,:,2) as the next time slice.
    call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,1,i), taum,i,Time)
          clim_type%time_init(i,1) = taum
          itaum = 1
    call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,2,i), taup,i,Time)
          clim_type%time_init(i,2) = taup
          itaup = 2
      endif ! itaum.eq.itaup.eq.0
      if (itaum.eq.0 .and. itaup.ne.0) then
! Can't think of a situation where we would have the next time level but not the previous.
 call mpp_error(FATAL,'interpolater : No data from the previous climatology time but we have the next time. How did this happen?')
      endif
      if (itaum.ne.0 .and. itaup.eq.0) then
!We have the previous time step but not the next time step data
        itaup = 1
        if (itaum .eq. 1 ) itaup = 2
        call read_data(clim_type,clim_type%field_type(i), clim_type%data(:,:,:,itaup,i), taup,i, Time)
        clim_type%time_init(i,itaup)=taup
      endif



    endif! TIME_FLAG

!select case(clim_type%TIME_FLAG)
!  case (LINEAR)
    hinterp_data = (1-tweight)*clim_type%data(istart:iend,jstart:jend,:,itaum,i) &
    + tweight*clim_type%data(istart:iend,jstart:jend,:,itaup,i)
!  case (SEASONAL)
! Do sine fit to data at this point

!end select

select case(clim_type%level_type)
  case(PRESSURE)
    p_fact = 1.0
  case(SIGMA)
    p_fact = maxval(phalf,3)! max pressure in the column !(:,:,size(phalf,3))
end select

col_data(:,:)=0.0
select case(clim_type%mr(i))
  case(NO_CONV)
    do k = 1,size(hinterp_data,3)
   col_data(:,:) = col_data(:,:) + hinterp_data(:,:,k)* &
      (clim_type%halflevs(k+1)-clim_type%halflevs(k))/grav
    enddo
    
  case(KG_M2)
    do k = 1,size(hinterp_data,3)
       col_data(:,:) = col_data(:,:) + hinterp_data(:,:,k)
       hinterp_data(:,:,k) = hinterp_data(:,:,k)/ &
         ((clim_type%halflevs(k+1)-clim_type%halflevs(k))*p_fact)
    enddo
end select

found = .false.
do j = 1,size(climo_diag_name)
  if (climo_diag_name(j) .eq. clim_type%field_name(i)) then
    found = .true.
    exit
  endif
enddo

if (found) then
  if (hinterp_id(j) > 0 ) then
       result = send_data(hinterp_id(j),col_data,Time)
  endif
endif


!++lwh
do j = 1, size(phalf,2)
   do ilon=1,size(phalf,1)
      pclim = p_fact(ilon,j)*clim_type%halflevs
      if ( maxval(phalf(ilon,j,:)) > maxval(pclim) ) then
         call mpp_error(NOTE,"Interpolater: model surface pressure&
                             & is greater than climatology surface pressure for "&
                             // trim(clim_type%file_name))
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( maxloc(pclim) ) = maxval( phalf(ilon,j,:) )
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
         end select
      endif
      if ( minval(phalf(ilon,j,:)) < minval(pclim) ) then
         call mpp_error(NOTE,"Interpolater: model top pressure&
                             & is less than climatology top pressure for "&
                             // trim(clim_type%file_name))
         select case(clim_type%out_of_bounds(i))
            case(CONSTANT)
               pclim( minloc(pclim) ) = minval( phalf(ilon,j,:) )
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
         end select
      endif
      select case(clim_type%vert_interp(i))
         case(INTERP_WEIGHTED_P)
            call interp_weighted_scalar(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
         case(INTERP_LINEAR_P)
            call interp_linear(pclim, phalf(ilon,j,:),hinterp_data(ilon,j,:),interp_data(ilon,j,:))
!        case(INTERP_LOG)
      end select
   enddo
enddo
!--lwh

select case(clim_type%mr(i))
  case(KG_M2)
    do k = 1,size(interp_data,3)
       interp_data(:,:,k) = interp_data(:,:,k)*(phalf(:,:,k+1)-phalf(:,:,k))
    enddo
end select

  endif !field_name
enddo !End of i loop

end subroutine interpolater
!
!#######################################################################
!
subroutine interpolater_end(clim_type)
! Subroutine to deallocate the interpolate type clim_type.
!
! INTENT INOUT
!  clim_type : allocate type whose components will be deallocated.
!
type(interpolate_type), intent(inout) :: clim_type
integer :: log_unit

#ifdef FEZ
call mpp_open(log_unit,file='logfile.out', action=MPP_APPEND)
if ( mpp_pe() == 0 ) then
   write (log_unit,'(/,(a))') 'Exiting interpolater, have a nice day ...'
end if
call mpp_close(log_unit)
#else
if ( mpp_pe() == mpp_root_pe() ) then
   write (stdlog(),'(/,(a))') 'Exiting interpolater, have a nice day ...'
end if
#endif

deallocate(clim_type%lat)
deallocate(clim_type%lon)
deallocate(clim_type%latb)
deallocate(clim_type%lonb)
deallocate(clim_type%levs)
deallocate(clim_type%halflevs) 
call horiz_interp_end(clim_type%interph)
deallocate(clim_type%time_slice)
deallocate(clim_type%field_type)
deallocate(clim_type%field_name)
deallocate(clim_type%time_init)
deallocate(clim_type%mr)
deallocate(clim_type%data)

call mpp_close(clim_type%unit)

module_is_initialized = .false.

end subroutine interpolater_end
!
!#######################################################################
!
subroutine read_data(clim_type,src_field, hdata, nt,i, Time)
!
!  INTENT IN
!    clim_type : The interpolate type which contains the data 
!    src_field : The field type 
!    nt        : The index of the time slice of the climatology that you wish to read.
!    i         : The index of the field name that you are trying to read. (optional)
!    Time      : The model time. Used for diagnostic purposes only. (optional)
!
!  INTENT OUT
!
!    hdata     : The horizontally interpolated climatology field. This 
!                field will still be on the climatology vertical grid.
!
type(interpolate_type)   , intent(in)  :: clim_type
type(fieldtype)          , intent(in)  :: src_field
integer                  , intent(in)  :: nt
real                     , intent(out) :: hdata(:,:,:)
integer        , optional, intent(in)  :: i
type(time_type), optional, intent(in)  :: Time

allocate(climdata(size(clim_type%lon),size(clim_type%lat),size(clim_type%levs)))
      call mpp_read(clim_type%unit,src_field, climdata,nt)
      call horiz_interp(clim_type%interph, climdata, hdata)
if(clim_diag_initialized) &
  call diag_read_data(clim_type,climdata,i, Time)
deallocate(climdata)


end subroutine read_data
!
!#######################################################################
!
subroutine diag_read_data(clim_type,model_data, i, Time)
!
! A routine to diagnose the data read in by read_data
!
!  INTENT IN
!    clim_type  : The interpolate type.
!    model_data : The data read in from file that is being diagnosed.
!    i          : The index of the field name that you are diagnosing.
!    Time       : The model time
!
type(interpolate_type), intent(in) :: clim_type
real                  , intent(in) :: model_data(:,:,:)
integer               , intent(in) :: i
type(time_type)       , intent(in) :: Time

integer :: j,k
real :: col_data(size(model_data,1),size(model_data,2))
logical :: result, found


found = .false.
do j = 1,size(climo_diag_name)
  if (climo_diag_name(j) .eq. clim_type%field_name(i)) then
      found = .true.
      exit
  endif
enddo

if(found) then
  if(climo_diag_id(j)>0) then
  col_data(:,:)=0.0
    do k=1,size(model_data,3)
      col_data(:,:) = col_data(:,:) + &
        model_data(:,:,k)* &
        (clim_type%halflevs(k+1)-clim_type%halflevs(k))/grav
    enddo
    result = send_data(climo_diag_id(j),col_data(clim_type%is:clim_type%ie,clim_type%js:clim_type%je),Time)
  endif
endif

end subroutine diag_read_data
!
!#######################################################################
!
function chomp(string)
!
! A function to remove CHAR(0) from the end of strings read from NetCDF files.
!
character(len=*), intent(in) :: string
character(len=64) :: chomp

integer :: len

len = len_trim(string)
if (string(len:len) == CHAR(0)) len = len -1

chomp = string(:len)

end function chomp
!
!#################################################################
!
subroutine interp_weighted_scalar ( grdin, grdout, datin, datout )
real, intent(in),  dimension(:) :: grdin, grdout, datin
real, intent(out), dimension(:) :: datout

integer :: j, k

if (size(grdin).ne. (size(datin)+1)) &
 call mpp_error(FATAL,'interp_weighted_scalar : input data and pressure do not have the same number of levels')
if (size(grdout).ne. (size(datout)+1)) &
 call mpp_error(FATAL,'interp_weighted_scalar : output data and pressure do not have the same number of levels')

  do k = 1, size(datout)

     datout(k) = 0.0

     do j = 1, size(datin)

        if ( grdin(j)   <= grdout(k) .and. &
             grdin(j+1) >= grdout(k) .and. &
             grdin(j+1) <= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdin(j+1)-grdout(k))

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j)   <= grdout(k+1) .and. &
                  grdin(j+1) >= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdout(k+1)-grdin(j))

        else if ( grdin(j)   >= grdout(k)   .and. &
                  grdin(j+1) <= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdin(j+1)-grdin(j))

        else if ( grdin(j)   <= grdout(k)   .and. &
                  grdin(j+1) >= grdout(k+1) ) then

           datout(k) = datout(k) + datin(j)*(grdout(k+1)-grdout(k))

        endif

     enddo

     datout(k) = datout(k)/(grdout(k+1)-grdout(k))

  enddo

end subroutine interp_weighted_scalar
!
!#################################################################
!
subroutine interp_linear ( grdin, grdout, datin, datout )
real, intent(in),  dimension(:) :: grdin, grdout, datin
real, intent(out), dimension(:) :: datout

integer :: j, k, n
real    :: wt


if (size(grdin).ne. (size(datin)+1)) &
 call mpp_error(FATAL,'interp_linear : input data and pressure do not have the same number of levels')
if (size(grdout).ne. (size(datout)+1)) &
 call mpp_error(FATAL,'interp_linear : output data and pressure do not have the same number of levels')


  n = size(grdin)

  do k= 1, size(datout)

   ! ascending grid values
     if (grdin(1) < grdin(n)) then
         do j = 2, size(grdin)-1
           if (grdout(k) <= grdin(j)) exit
         enddo
   ! descending grid values
     else
         do j = size(grdin), 3, -1
           if (grdout(k) <= grdin(j-1)) exit
         enddo
     endif

   ! linear interpolation
     wt = (grdout(k)-grdin(j-1)) / (grdin(j)-grdin(j-1))
!print '(a,2i3,4f6.1)', 'k,j=',k,j,grdout(k),grdin(j-1),grdin(j),wt
   ! constant value extrapolation
   ! wt = min(max(wt,0.),1.)

     datout(k) = (1.-wt)*datin(j-1) + wt*datin(j)
     
  enddo

end subroutine interp_linear
!
!########################################################################
!
end module interpolater_mod
!
!#######################################################################
!
#ifdef test_interp
program test

use mpp_mod
use mpp_io_mod
use mpp_domains_mod
use time_manager_mod
use diag_manager_mod!, only : diag_axis_init, file_exist, MPP_NPES, &
                    !  MPP_PE, REGISTER_DIAG_FIELD, SEND_DATA, SET_DATE,&
                    !  SET_TIME

use interpolater_mod
!use sulfate_mod
!use ozone_mod
use constants_mod, only : grav

implicit none
integer, parameter :: nxd = 144, nyd = 90, ntsteps = 30, delt = 1, two_delt = 2*delt
integer, parameter :: max_fields = 20 ! maximum number of fields to be interpolated

integer :: i,k,n,level
integer :: unit, io_status
integer :: ndivs
integer :: jscomp, jecomp, iscomp, iecomp, isd,ied,jsd,jed
integer :: numfields, domain_layout(2)
integer :: num_nbrs, nbins,axes(3), interp_diagnostic_id
integer :: column_diagnostic_id1, column_diagnostic_id(max_fields)

real ::  missing_value = -1.e10

character(len=1) :: dest_grid
character(len=128) :: src_file, file_out, title, units, colaer
logical :: vector_field=.false., result

type(axistype), allocatable, dimension(:)  :: axes_out, axes_src
type(axistype) :: time_axis
type(fieldtype), allocatable, dimension(:) :: fields
type(fieldtype) :: dest_field(max_fields), src_field(max_fields), field_geolon_t, &
     field_geolat_t, field_geolon_c, field_geolat_c
type(atttype), allocatable, dimension(:) :: global_atts
type(domain2d) :: domain
type(time_type) :: model_time

type(interpolate_type) :: o3, aerosol

real, dimension(:,:), allocatable :: col_data
real, dimension(:,:,:), allocatable :: model_data, p_half, p_full
real, dimension(:), allocatable :: latb_mod(:),lonb_mod(:),lon_mod(:),lat_mod(:)
real :: dx,dy
real :: dtr,tpi
real :: p_bot,p_top,lambda
character(len=64) :: names(13)
data names(:) /"so4_anthro","so4_natural","organic_carbon","black_carbon","sea_salt",&
"anthro_dust_0.2","anthro_dust_0.8","anthro_dust_2.0","anthro_dust_8.0",&
"natural_dust_0.2","natural_dust_0.8","natural_dust_2.0","natural_dust_8.0"/

integer :: out_of_bounds(2)
data out_of_bounds / CONSTANT, CONSTANT/!, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, &
!ZERO, ZERO, ZERO, ZERO /

namelist /interpolater_nml/ src_file

level = 18
tpi = 4.*acos(0.)
dtr = tpi/360.

src_file = 'src_file'  ! input file containing fields to be interpolated

! initialize communication modules

call mpp_init
call mpp_io_init
call mpp_domains_init
call set_calendar_type(JULIAN)
call diag_manager_init

model_time = set_date(1980,1,1,0,0,0)

!if (numfields.ne.2.and.vector_field) call mpp_error(FATAL,'2 components of vector field not specified')
!if (numfields.gt.1.and..not.vector_field) call mpp_error(FATAL,'only 1 scalar at a time')
!if (numfields .gt. max_fields) call mpp_error(FATAL,'max num fields exceeded')

!--------------------------------------------------------------------
! namelist input
!--------------------------------------------------------------------

call mpp_open(unit, 'input.nml',  action=MPP_RDONLY, form=MPP_ASCII)
read  (unit, interpolater_nml,iostat=io_status)
#ifndef FEZ
write (stdlog(), nml=interpolater_nml)
#endif
if (io_status .gt. 0) then
    call mpp_error(FATAL,'=>Error reading interpolater_nml')
endif
call mpp_close(unit)


! decompose model grid points
! mapping can get expensive so we distribute the task at this level

ndivs = mpp_npes()

call mpp_define_layout ((/1,nxd,1,nyd/), ndivs, domain_layout)
call mpp_define_domains((/1,nxd,1,nyd/),domain_layout, domain,xhalo=0,yhalo=0)  
call mpp_get_data_domain(domain,isd,ied,jsd,jed)
call mpp_get_compute_domain (domain, iscomp, iecomp, jscomp, jecomp)

allocate(lonb_mod(nxd+1),lon_mod(nxd))
allocate(latb_mod(nyd+1),lat_mod(nyd))
allocate(col_data(isd:ied,jsd:jed))
allocate(p_half(isd:ied,jsd:jed,level+1),p_full(isd:ied,jsd:jed,level))
p_top = 1.0
p_bot = 101325.0 !Model level in Pa
lambda = -1.0*log(p_top/p_bot)/(level+1)

p_half(:,:,level+1) = p_bot
do i=level,1,-1
  p_half(:,:,i)=p_half(:,:,i+1)*exp(-1.0*lambda)
enddo
do i=1,level
  p_full(:,:,i)=(p_half(:,:,i+1)+p_half(:,:,i))/2.0
enddo

allocate(model_data(isd:ied,jsd:jed,level))

dx = 360./nxd
dy = 180./nyd
do i = 1,nxd+1
  lonb_mod(i) = (i-1)*dx 
enddo
do i = 1,nyd+1
  latb_mod(i) = -90. + (i-1)*dy 
enddo
do i=1,nxd
  lon_mod(i)=(lonb_mod(i+1)+lonb_mod(i))/2.0
enddo
do i=1,nyd
  lat_mod(i)=(latb_mod(i+1)+latb_mod(i))/2.0
enddo

lonb_mod = lonb_mod * dtr
latb_mod = latb_mod * dtr

   axes(1) = diag_axis_init('x',lon_mod,units='degrees',cart_name='x',domain2=domain)
   axes(2) = diag_axis_init('y',lat_mod,units='degrees',cart_name='y',domain2=domain)
   axes(3) = diag_axis_init('z',p_full(isd,jsd,:),units='mb',cart_name='z')

interp_diagnostic_id =  register_diag_field('interp','ozone',axes(1:3),model_time,&
                                'interpolated_ozone_clim', 'kg/kg', missing_value)      
column_diagnostic_id1 =  register_diag_field('interp','colozone',axes(1:2),model_time,&
                                'column_ozone_clim', 'kg/m2', missing_value)      

do i=1,size(names)
colaer = 'col'//trim(names(i))
column_diagnostic_id(i) =  register_diag_field('interp',colaer,axes(1:2),model_time,&
                                'column_aerosol_clim', 'kg/m2', missing_value)      
enddo


call ozone_init(o3,lonb_mod(isd:ied+1), latb_mod(jsd:jed+1), axes, model_time, data_out_of_bounds=out_of_bounds)
call sulfate_init(aerosol,lonb_mod(isd:ied+1), latb_mod(jsd:jed+1), names, data_out_of_bounds=(/CONSTANT/) )
call init_clim_diag(aerosol, axes, model_time)

do n=1,ntsteps
#ifdef FEZ
if( mpp_pe().eq.0) write(*,*) n
#else
if( mpp_pe() == mpp_root_pe() ) write(*,*) n
#endif

call get_ozone(o3,model_time,p_half,model_data)

if(interp_diagnostic_id>0) &
     result = send_data(interp_diagnostic_id,&
          model_data(iscomp:iecomp,jscomp:jecomp,:),model_time)

if(column_diagnostic_id1>0) then

col_data(iscomp:iecomp,jscomp:jecomp)=0.0
do k=1,level
   col_data(iscomp:iecomp,jscomp:jecomp)= col_data(iscomp:iecomp,jscomp:jecomp)+ &
      model_data(iscomp:iecomp,jscomp:jecomp,k)* &
      (p_half(iscomp:iecomp,jscomp:jecomp,k+1)-p_half(iscomp:iecomp,jscomp:jecomp,k))/grav
enddo
     result = send_data(column_diagnostic_id1,col_data(:,:),model_time)
endif



do i=1,size(names)

call get_anthro_sulfate(aerosol,model_time,p_half,names(i),model_data,clim_units=units)

if(column_diagnostic_id(i)>0) then

col_data(iscomp:iecomp,jscomp:jecomp)=0.0
do k=1,level
if (trim(units) .eq. 'kg/m^2') then
   col_data(iscomp:iecomp,jscomp:jecomp)= col_data(iscomp:iecomp,jscomp:jecomp)+ &
      model_data(iscomp:iecomp,jscomp:jecomp,k)
else
   col_data(iscomp:iecomp,jscomp:jecomp)= col_data(iscomp:iecomp,jscomp:jecomp)+ &
      model_data(iscomp:iecomp,jscomp:jecomp,k)* &
      (p_half(iscomp:iecomp,jscomp:jecomp,k+1)-p_half(iscomp:iecomp,jscomp:jecomp,k))/grav
endif
enddo
     result = send_data(column_diagnostic_id(i),&
col_data(iscomp:iecomp,jscomp:jecomp),model_time)
endif

enddo

   model_time = model_time + set_time(0,delt)      

   if (n.eq. ntsteps) call diag_manager_end(model_time)

enddo

call interpolater_end(aerosol)
call interpolater_end(o3)

deallocate(lonb_mod, lon_mod, latb_mod,lat_mod, col_data, p_half, p_full, model_data)

call mpp_exit

contains
!
!#######################################################################
!
subroutine sulfate_init(aerosol,lonb, latb, names, data_out_of_bounds, vert_interp, units)
type(interpolate_type), intent(inout)         :: aerosol
real,                   intent(in)            :: lonb(:),latb(:)
character(len=64),      intent(in)            :: names(:)
integer,                intent(in)            :: data_out_of_bounds(:) 
integer,                intent(in), optional  :: vert_interp(:)
character(len=*),       intent(out),optional  :: units(:)

#ifdef FEZ
call interpolater_init( aerosol, "aerosol.climatology", lonb, latb, &
#else
call interpolater_init( aerosol, "aerosol.climatology.nc", lonb, latb, &
#endif
                        data_names=names, data_out_of_bounds=data_out_of_bounds, &
                        vert_interp=vert_interp, clim_units=units )

end subroutine sulfate_init
!
!#######################################################################
!
subroutine get_anthro_sulfate( sulfate, model_time, p_half, name, model_data, is, js, clim_units )
type(interpolate_type), intent(out) :: sulfate
type(time_type), intent(in) :: model_time
real, intent(in)           :: p_half(:,:,:)
character(len=*), intent(in) :: name
character(len=*), intent(out), optional :: clim_units
real, intent(out) :: model_data(:,:,:)
integer, intent(in), optional :: is,js

call interpolater( sulfate, model_time, p_half, model_data, name, is, js, clim_units)

end subroutine get_anthro_sulfate
!
!#######################################################################
!
subroutine ozone_init( o3, lonb, latb, axes, model_time, data_out_of_bounds, vert_interp )
real,                  intent(in)           :: lonb(:),latb(:)
integer,               intent(in)           :: axes(:)
type(time_type),       intent(in)           :: model_time
type(interpolate_type),intent(out)          :: o3
integer,               intent(in)           :: data_out_of_bounds(:)
integer,               intent(in), optional :: vert_interp(:)

#ifdef FEZ
call interpolater_init( o3, "o3.climatology", lonb, latb, &
#else
call interpolater_init( o3, "o3.climatology.nc", lonb, latb, &
#endif
                        data_out_of_bounds=data_out_of_bounds, vert_interp=vert_interp )

end subroutine ozone_init
!
!#######################################################################
!
subroutine get_ozone( o3, model_time, p_half, model_data, is, js )
type(interpolate_type),intent(inout) :: o3
type(time_type), intent(in) :: model_time
real, intent(in)           :: p_half(:,:,:)
real, intent(out) :: model_data(:,:,:)
integer, intent(in), optional :: is,js

call interpolater( o3, model_time, p_half, model_data, "ozone", is, js)

end subroutine get_ozone

end program test

#endif
