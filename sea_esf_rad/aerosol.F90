                    module aerosol_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!  smf
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Code to initialize/allocate aerosol climatology
! </OVERVIEW>
! <DESCRIPTION>
!  This code initializes prescribed aerosol climatology from input file,
!  allocates necessary memory space and interpolate the aerosol climatology
!  to the model specification. Afterwards the memory space is deallocated,
!  the aerosol climatology information is freed.
! </DESCRIPTION>
!  shared modules:

use time_manager_mod,  only: time_type, time_manager_init
use fms_mod,           only: open_namelist_file, fms_init, &
                             mpp_pe, mpp_root_pe, stdlog, &
                             file_exist, write_version_number, &
                             check_nml_error, error_mesg, &
                             FATAL, NOTE, WARNING, close_file
use interpolator_mod,  only: interpolate_type, interpolator_init, &
                             interpolator, interpolator_end, &
                             CONSTANT, INTERP_WEIGHTED_P
use mpp_io_mod,        only: mpp_open, mpp_close, MPP_RDONLY,   &
                             MPP_ASCII, MPP_SEQUENTIAL, MPP_MULTI,  &
                             MPP_SINGLE, mpp_io_init

!  shared radiation package modules:

use rad_utilities_mod, only: aerosol_type, rad_utilities_init   

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    aersosol_mod defines the aerosol distribution that is to be seen
!    by the radiation package.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128) :: version = '$Id: aerosol.F90,v 10.0 2003/10/24 22:00:39 fms Exp $'
character(len=128) :: tagname = '$Name: jakarta $'


!-----------------------------------------------------------------------
!------  interfaces -------

public          &
        aerosol_init, aerosol_driver, aerosol_end, &
        aerosol_dealloc

!private         &


!---------------------------------------------------------------------
!------  namelist ------

integer, parameter     ::      &
        MAX_DATA_FIELDS = 100     ! maximum number of aerosol species
character(len=64)      ::      &
        data_names(MAX_DATA_FIELDS) = '  ' 
                                  ! names of active aerosol species 
character(len=64)      ::      &
        filename = '  '           ! name of netcdf file containing 
                                  ! aerosol species to be activated


namelist /aerosol_nml/                            &
                           data_names, filename               

!---------------------------------------------------------------------
!---- public data ----


!---------------------------------------------------------------------
!---- private data ----


!---------------------------------------------------------------------
!   the following is an interpolate_type variable containing the
!   information about the aerosol species.
!---------------------------------------------------------------------
type(interpolate_type), save  :: Aerosol_interp

!--------------------------------------------------------------------
!    miscellaneous variables
!--------------------------------------------------------------------
integer  :: nfields=0                        ! number of active aerosol 
                                             ! species
logical  :: module_is_initialized = .false.  ! module has been 
                                             ! initialized  ?



                           contains


 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------------
! <SUBROUTINE NAME="aerosol_init">
!  <OVERVIEW>
!   Subroutine to initialize/interpolate prescribed aerosol climatology
!  </OVERVIEW>
!  <DESCRIPTION>
!   Subroutine to initialize/interpolate prescribed aerosol climatology
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_init(lonb, latb, aerosol_names)
!  </TEMPLATE>
!  <IN NAME="lonb" TYPE="real">
!   Array of model longitudes on cell boundaries in [radians]
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   Array of model latitudes on cell boundaries in [radians]
!  </IN>
!  <IN NAME="aerosol_names" TYPE="character">
!   names of the activated aerosol species
!  </IN>
! </SUBROUTINE>
!
subroutine aerosol_init (lonb, latb, aerosol_names)

!-----------------------------------------------------------------------
!    aerosol_init is the constructor for aerosol_mod.
!-----------------------------------------------------------------------

real, dimension(:),              intent(in)  :: lonb,latb
character(len=64), dimension(:), pointer     :: aerosol_names

!----------------------------------------------------------------------
!  intent(in) variables:
!
!       lonb           array of model longitudes on cell boundaries 
!                      [ radians ]
!       latb           array of model latitudes at cell boundaries 
!                      [ radians ]
!
!   pointer variable:
!
!       aerosol_names  names of the activated aerosol species
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!  local variables:
      
      integer   ::   unit, ierr, io       
      integer   ::   n

!---------------------------------------------------------------------
!    local variables:
!
!         unit       io unit number used to read namelist file
!         ierr       error code
!         io         error status returned from io operation
!         n          do-loop index
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
      call time_manager_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosol_nml, iostat=io,  &
               end=10)
        ierr = check_nml_error(io,'aerosol_nml')
        end do
10      call close_file (unit)   
      endif                      
                                  
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
                        write (stdlog(), nml=aerosol_nml)

!-----------------------------------------------------------------------
!    count number of activated aerosol species. allocate a pointer 
!    array to return the names of the activated species to the calling 
!    routine.
!-----------------------------------------------------------------------
      do n=1,MAX_DATA_FIELDS
        if (data_names(n) /= ' '  ) then
          nfields = n
        else
          exit
        endif
      end do
      allocate (aerosol_names(nfields))
      aerosol_names (:) = data_names(1:nfields)

!-----------------------------------------------------------------------
!    initialize the interpolation module if any aerosol species have
!    been activated.
!-----------------------------------------------------------------------
      if (nfields > 0) then
        call interpolator_init (Aerosol_interp, filename, lonb, latb, &
                                data_names(:nfields),   &
                                data_out_of_bounds=(/CONSTANT/), &
                                vert_interp=(/INTERP_WEIGHTED_P/) )
      endif

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------



end subroutine aerosol_init



!######################################################################
! <SUBROUTINE NAME="aerosol_driver">
!  <OVERVIEW>
!   Interpolate aerosol verical profile based on prescribed aerosol
!   climatology input and model set up.
!  </OVERVIEW>
!  <TEMPLATE>
!   call aerosol_driver (is, js, model_time, p_half, Aerosol)
!  </TEMPLATE>
!  <INOUT NAME="Aerosol" TYPE="aerosol_type">
!   Aerosol climatology input
!  </INOUT>
!  <IN NAME="model_time" TYPE="time_type">
!   The internal model simulation time, i.e. Jan. 1 1982
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   The array of model layer pressure values
!  </IN>
!  <IN NAME="is" TYPE="integer">
!   The longitude index of model physics window domain
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   The latitude index of model physics window domain
!  </IN>
! </SUBROUTINE>
!
subroutine aerosol_driver (is, js, model_time, p_half, Aerosol)

!-----------------------------------------------------------------------
!    aerosol_driver returns the names and concentrations of activated 
!    aerosol species at model grid points at the model_time to the 
!    calling routine in aerosol_type variable Aerosol. 
!-----------------------------------------------------------------------

integer,                  intent(in)     :: is,js
type(time_type),          intent(in)     :: model_time
real, dimension(:,:,:),   intent(in)     :: p_half
type(aerosol_type),       intent(inout)  :: Aerosol

!--------------------------------------------------------------------
!   intent(in) variables:
!
!       is, js           starting subdomain i,j indices of data in 
!                        the physics_window being integrated
!       model_time       time for which aerosol data is desired
!                        [ time_type ]
!       p_half           model pressure at interface levels 
!                        [ Pa ]
!      
!   intent(inout) variables:
!
!       Aerosol    aerosol_type variable. the following components will
!                  be returned from this routine:
!                   aerosol      concentration of each active aerosol 
!                                species at each model grid point
!                                [ kg / m**2 ]
!                   aerosol_names 
!                                names assigned to each active species
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer :: n    ! do-loop index

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosol_mod',   &
                         'module has not been initialized',FATAL )
      endif

!--------------------------------------------------------------------
!    allocate and define an array to hold the activated aerosol names.
!---------------------------------------------------------------------
      allocate (Aerosol%aerosol_names (nfields))
      Aerosol%aerosol_names = data_names(:nfields) 

!---------------------------------------------------------------------
!    allocate an array to hold the aerosol amounts for each species at
!    each grid point. call interpolator to obtain data valid at 
!    model_time for each activated aerosol species. store the aerosol
!    amounts in Aerosol%aerosol.
!----------------------------------------------------------------------
      allocate (Aerosol%aerosol(size(p_half,1), size(p_half,2), &
                                size(p_half,3) - 1, nfields))
      do n=1,nfields
        call interpolator (Aerosol_interp, model_time, p_half, &
                           Aerosol%aerosol(:,:,:,n),    &
                           Aerosol%aerosol_names(n), is, js)
      end do 

!---------------------------------------------------------------------


end subroutine aerosol_driver



!#####################################################################
! <SUBROUTINE NAME="aerosol_end">
!  <OVERVIEW>
!   aerosol_end is the destructor for aerosol_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!   aerosol_end is the destructor for aerosol_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine aerosol_end

!----------------------------------------------------------------------
!    aerosol_end is the destructor for aerosol_mod.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosol_mod',   &
                         'module has not been initialized',FATAL )
      endif

!---------------------------------------------------------------------
!    call interpolator_end to release the interpolate_type variable 
!    used in this module.
!---------------------------------------------------------------------
      if (nfields > 0) then
        call interpolator_end (Aerosol_interp)
      endif

!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------



end subroutine aerosol_end



!#####################################################################
! <SUBROUTINE NAME="aerosol_dealloc">
!  <OVERVIEW>
!    aerosol_dealloc deallocates the array components of an 
!    aersol_type derived type variable.
!  </OVERVIEW>
!  <DESCRIPTION>
!    aerosol_dealloc deallocates the array components of an 
!    aersol_type derived type variable.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call aerosol_dealloc
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine aerosol_dealloc (Aerosol)

!---------------------------------------------------------------------
!    aerosol_dealloc deallocates the array components of an 
!    aersol_type derived type variable.
!---------------------------------------------------------------------

type(aerosol_type), intent(inout) :: Aerosol

!---------------------------------------------------------------------
!  intent(inout) variables:
! 
!      Aerosol       aerosol_type variable containing information on
!                    the activated aerosol species
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('aerosol_mod',   &
                         'module has not been initialized',FATAL )
      endif

!---------------------------------------------------------------------
!    deallocate the components of the aerosol_type variable.
!---------------------------------------------------------------------
      deallocate (Aerosol%aerosol)
      deallocate (Aerosol%aerosol_names)
 
!----------------------------------------------------------------------


end subroutine aerosol_dealloc


!#####################################################################



                    end module aerosol_mod


!=======================================================================



#ifdef test_aerosol

program main

use aerosol_mod
use mpp_mod
use mpp_io_mod
use mpp_domains_mod
use time_manager_mod
use diag_manager_mod
use rad_utilities_mod




implicit none

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
integer, parameter :: NLON=20, NLAT=10,NLEV=8
integer, parameter :: MAX_AERSOL_NAMES = 100
real :: latb(NLAT+1),lonb(NLON+1),pi,phalf(NLON,NLAT,NLEV+1)
integer :: i,nspecies
type(time_type) :: model_time
character(len=64), dimension(MAX_AEROSOL_NAMES) :: names
type(aerosol_type)  :: Aerosol

pi = 4.*atan(1.)

call mpp_init
call mpp_io_init
call mpp_domains_init
call diag_manager_init
call set_calendar_type(JULIAN)

do i = 1,NLAT+1
   latb(i) = -90. + 180.*REAL(i-1)/REAL(NLAT)
end do
do i = 1,NLON+1
   lonb(i) = -180. + 360.*REAL(i-1)/REAL(NLON)
end do

latb(:) = latb(:) * pi/180.
lonb(:) = lonb(:) * pi/180.

do i = 1,NLEV+1
   phalf(:,:,i) = 101325. * REAL(i-1) / REAL(NLEV)
end do

model_time = set_date(1980,1,1,0,0,0)

call aerosol_init (lonb, latb, names)

call aerosol_driver (1,1,model_time, phalf, Aerosol)

call aerosol_dealloc (Aerosol)

call aerosol_end

call mpp_exit

end program main

#endif
