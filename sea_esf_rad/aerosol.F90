
                    module aerosol_mod

!=======================================================================

use time_manager_mod,    only: time_type
use utilities_mod,       only: error_mesg, FATAL, file_exist,   &
                               check_nml_error, open_file,      &
                               get_my_pe, close_file
use interpolator_mod,    only: interpolate_type, interpolator_init, &
                               interpolator, interpolator_end, &
                               CONSTANT, INTERP_WEIGHTED_P
use mpp_io_mod,          only: mpp_open, mpp_close, MPP_RDONLY,   &
                               MPP_ASCII, MPP_SEQUENTIAL, MPP_MULTI,  &
			       MPP_SINGLE
use rad_utilities_mod, only : aerosol_type, radiation_control_type, &
			      aerosol_properties_type, Aerosol_props, &
                              Rad_control

implicit none
private

public aerosol_init, &
       aerosol_driver, aerosol_end, &
       get_aerosol_info, &
       get_aerosol_optical_info, get_aerosol_optical_index, &
!      interp_aerosol_end, &
       aerosol_alloc, aerosol_dealloc
!=======================================================================

character(len=128) :: version = '$Id: aerosol.F90,v 1.2 2003/04/09 20:58:36 fms Exp $'
character(len=128) :: tag = '$Name: inchon $'

!-----------------------------------------------------------------------



type(interpolate_type) :: aerosol_interp_type
integer, parameter     :: MAX_DATA_FIELDS = 100
character(len=64)      :: data_names(MAX_DATA_FIELDS), filename, &
                          optical_filename
integer                :: nfields=0
logical                :: aerosol_initialized = .false.

!-----------------------------------------------------------------------
! ... Aerosol optical property data
!-----------------------------------------------------------------------
integer, parameter                 :: MAX_OPTICAL_FIELDS = 100
integer                            :: NAERMODELS, num_wavenumbers
real, dimension(:,:), allocatable  :: aeroextivl, &
                                      aerossalbivl, &
                                      aeroasymmivl
integer, dimension(:), allocatable :: endaerwvnsf
integer :: sulfate_index(0:100), optical_index(MAX_DATA_FIELDS)

character(len=64)                  :: aerosol_optical_names(MAX_OPTICAL_FIELDS)


namelist /aerosol_nml/ data_names,            & 
                                   ! names of species to read from file
                       filename,              & 
                                 ! NetCDF filename (omit ".nc" for Fez)
                       aerosol_optical_names, & 
		       ! Names of aerosol optical property categories
                       num_wavenumbers, &   
             ! Number of wavenumber bins for aerosol optical properties
		       optical_filename   
	       	                    ! Optical property filename (ASCII)

CONTAINS
 
!=======================================================================

subroutine aerosol_init(lonb, latb, nlev)
!-----------------------------------------------------------------------
! ... Initialize aerosol interpolation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------
real,    intent(in) :: lonb(:),latb(:)
integer, intent(in) :: nlev

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
integer :: unit, ierr, io_status, n

if( aerosol_initialized ) return

!--------------------------------------------------------------------
! ... Initialize namelist variables
!---------------------------------------------------------------------

data_names(:)            = ' '
aerosol_optical_names(:) = ' '
filename                 = ' '
optical_filename         = ' '
naermodels = 0
num_wavenumbers = 0
sulfate_index = 0
optical_index =0

!-----------------------------------------------------------------------
! ... Read namelist
!-----------------------------------------------------------------------
if ( file_exist('input.nml')) then
   call mpp_open(unit, 'input.nml',  action=MPP_RDONLY, form=MPP_ASCII)
   read  (unit, aerosol_nml,iostat=io_status)
   ierr = check_nml_error(io_status,'aerosol_nml')
   call mpp_close(unit)
end if

  unit = open_file ('logfile.out', action='append')
      if (get_my_pe() == 0) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=aerosol_nml)
      endif
     call close_file (unit)

!-----------------------------------------------------------------------
! ... Count number of species
!-----------------------------------------------------------------------
do n=1,MAX_DATA_FIELDS
   if( data_names(n) /= ' '  ) then
      nfields = n
   else
      exit
   end if
end do
   if (Rad_control%do_aerosol .and. nfields == 0) then
     call error_mesg ('aerosol_mod', &
     ' aerosols desired but no aerosol data_names supplied', FATAL)
   endif

!-----------------------------------------------------------------------
! ... Initialize interpolation
!-----------------------------------------------------------------------

     if (Rad_control%do_aerosol) then

 call interpolator_init( aerosol_interp_type, filename, lonb, latb, &
                         data_names(:nfields), data_out_of_bounds=(/CONSTANT/), &
                         vert_interp=(/INTERP_WEIGHTED_P/) )

!-----------------------------------------------------------------------
! ... Initialize aerosol optical properties
!-----------------------------------------------------------------------
 call aerosol_optical_init

endif

Aerosol_props%nfields = nfields
allocate (Aerosol_props%aerosol_names(nfields))
Aerosol_props%aerosol_names(:) = data_names(:nfields)

aerosol_initialized = .true.

end subroutine aerosol_init

!=======================================================================

subroutine aerosol_optical_init
!-----------------------------------------------------------------------
! ... Initialize aerosol optical properties
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
integer           :: n, noptical
integer           :: unit, num_input_categories
character(len=64) :: name_in, target_name
real, allocatable, dimension(:)    :: aeroext_in, aerossalb_in,   &
                                      aeroasymm_in
logical, allocatable, dimension(:) :: found

!----------------------------------------------------------------------
! ... Count number of aerosol optical property categories           
!----------------------------------------------------------------------
do n=1,MAX_OPTICAL_FIELDS
   if( aerosol_optical_names(n) /= ' '  ) then
      NAERMODELS = n
   else
      exit
   end if
end do

!----------------------------------------------------------------------
! ... Open ASCII input file
!----------------------------------------------------------------------

call mpp_open( unit, 'INPUT/'//optical_filename, MPP_RDONLY, MPP_ASCII, &
               MPP_SEQUENTIAL, MPP_MULTI, MPP_SINGLE )


!----------------------------------------------------------------------
! ... Read dimension info from input file
!----------------------------------------------------------------------

read( unit,* ) num_wavenumbers
read( unit,* ) num_input_categories

!----------------------------------------------------------------------
! ... Read wavenumber limits from input file
!     (for shettle & fenn aerosol intervals)
!----------------------------------------------------------------------
 allocate( endaerwvnsf(num_wavenumbers) )
 read( unit,* )
 read( unit,* ) endaerwvnsf
 
!----------------------------------------------------------------------
! ... Allocate module data arrays
!----------------------------------------------------------------------
allocate( aeroextivl(  num_wavenumbers, NAERMODELS), &
          aerossalbivl(num_wavenumbers, NAERMODELS), &
          aeroasymmivl(num_wavenumbers, NAERMODELS) )

!----------------------------------------------------------------------
! ... Allocate local working arrays
!----------------------------------------------------------------------
 allocate( aeroext_in( num_wavenumbers ),             &
           aerossalb_in( num_wavenumbers ),           &
           aeroasymm_in( num_wavenumbers ),           &
           found( NAERMODELS ) )

!----------------------------------------------------------------------
! ... Match names of optical property categories from input file with
!     those specified in namelist
!----------------------------------------------------------------------
 found(:) = .false.

do n = 1,num_input_categories
   read( unit,* ) name_in
    read( unit,* )
   read( unit,* ) aeroext_in
    read( unit,* )
   read( unit,* ) aerossalb_in
   read( unit,* )
   read( unit,* ) aeroasymm_in
   do noptical = 1,NAERMODELS
     if( aerosol_optical_names(noptical) == name_in ) then
          aeroextivl(:,noptical)   = aeroext_in
          aerossalbivl(:,noptical) = aerossalb_in
          aeroasymmivl(:,noptical) = aeroasymm_in
          found( noptical ) = .true.
          exit
       end if
    end do
end do

!----------------------------------------------------------------------
! ... Close ASCII input file
!----------------------------------------------------------------------
 call mpp_close( unit )

!----------------------------------------------------------------------
! ... Check to make sure data for all aerosol optical property
!     categories specified in namelist were contained in ASCII
!     input file
!----------------------------------------------------------------------
 do noptical = 1,NAERMODELS
    if( .not. found( noptical ) ) &
       call error_mesg( 'aerosol_optical_init', &
                   'Cannot find aerosol optical properties for ' // &
                    TRIM(aerosol_optical_names(noptical)), &
                     FATAL )
 end do

    allocate ( Aerosol_props%endaerwvnsf(num_wavenumbers))
allocate( Aerosol_props%aeroextivl(  num_wavenumbers, NAERMODELS), &
          Aerosol_props%aerossalbivl(num_wavenumbers, NAERMODELS), &
          Aerosol_props%aeroasymmivl(num_wavenumbers, NAERMODELS) )
    Aerosol_props%naermodels = naermodels
    Aerosol_props%num_wavenumbers = num_wavenumbers
    Aerosol_props%endaerwvnsf     = endaerwvnsf       
   Aerosol_props%aeroextivl  = aeroextivl
   Aerosol_props%aerossalbivl = aerossalbivl                    
   Aerosol_props%aeroasymmivl = aeroasymmivl                    

!----------------------------------------------------------------------
! ... Deallocate local working arrays
!----------------------------------------------------------------------
 deallocate( aeroext_in, aerossalb_in, aeroasymm_in, found )


 


!----------------------------------------------------------------------
! ... Match aerosol optical property indices with aerosol indices
!    Sulfate aerosols are handled separately (below) with RH dependence.
!----------------------------------------------------------------------
do n = 1,nfields
   name_in = data_names( n )
   optical_index(n) = 0
   if( name_in /= "so4_anthro" .and. name_in /= "so4_natural" ) then
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
      do noptical = 1,NAERMODELS
         if( aerosol_optical_names(noptical) == target_name ) then
            optical_index(n) = noptical
            exit
         end if
      end do
      if( optical_index(n) == 0 ) &
         call error_mesg( 'aerosol_optical_init', &
                          'Cannot find aerosol optical model = ' // &
                          TRIM( target_name ), &
                          FATAL )
   end if
end do
!----------------------------------------------------------------------
! ... Set up RH-dependent sulfate aerosol optical property indices
!----------------------------------------------------------------------
do n = 0,100
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
   sulfate_index(n) = 0
   do noptical = 1,NAERMODELS
      if( aerosol_optical_names(noptical) == target_name ) then
         sulfate_index(n) = noptical
         exit
      end if
   end do
   if( sulfate_index(n) == 0 ) &
      call error_mesg( 'aerosol_optical_init', &
                       'Cannot find aerosol optical model = ' // &
                       TRIM( target_name), &
                       FATAL )
end do


end subroutine aerosol_optical_init

!=======================================================================

subroutine interp_aerosol(Aerosol, model_time, p_half, is, js)
!-----------------------------------------------------------------------
! ... Interpolate aerosol fields to current model time and pressures
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------
type(aerosol_type), intent(inout)   :: Aerosol
type(time_type), intent(in)           :: model_time
real,            intent(in)           :: p_half(:,:,:)
integer,         intent(in), optional :: is,js

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
integer :: n

if( .not. aerosol_initialized ) &
   call error_mesg( 'interp_aerosol','must first call aerosol_init',FATAL )

!-----------------------------------------------------------------------
! ... Call interpolator for each species
!-----------------------------------------------------------------------
 do n = 1,nfields
    call interpolator( aerosol_interp_type, model_time, p_half, &
                       Aerosol%aerosol(:,:,:,n),    &
		       Aerosol_props%aerosol_names(n), is, js)
 end do 

end subroutine interp_aerosol

!=======================================================================

subroutine aerosol_driver( Aerosol, model_time, p_half, is, js)
!-----------------------------------------------------------------------
! ... Interpolate aerosol fields to current model time and pressures
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------
type(time_type), intent(in)           :: model_time
real,            intent(in)           :: p_half(:,:,:)
type(aerosol_type), intent(inout) :: Aerosol
integer,         intent(in), optional :: is,js

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
integer :: n
integer  :: ix, jx, kx

if( .not. aerosol_initialized ) &
   call error_mesg( 'aerosol_driver','must first call aerosol_init',FATAL )

  ix = size (p_half,1)
  jx = size (p_half,2)
  kx = size (p_half,3) - 1

  call aerosol_alloc ( Aerosol, ix, jx, kx)

!-----------------------------------------------------------------------
! ... Call aerosol interpolator
!-----------------------------------------------------------------------
call interp_aerosol( Aerosol, model_time, p_half, is, js )

!----- send aerosol field for diagnostic output -----



end subroutine aerosol_driver


!-----------------------------------------------------------------------

subroutine get_aerosol_all(aerosol_out)
!-----------------------------------------------------------------------
! ... Return aerosol fields (for all aerosol species)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------
real, dimension(:,:,:,:), intent(out)    :: aerosol_out

if( .not. aerosol_initialized ) &
   call error_mesg( 'get_aerosol_all','must first call aerosol_init',FATAL )

!aerosol_out(:,:,:,:) = aerosol3(:,:,:,:)

end subroutine get_aerosol_all


subroutine get_aerosol_single(aerosol_out, index)
!-----------------------------------------------------------------------
! ... Return aerosol fields (for single aerosol species, by index number)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------
real, dimension(:,:,:), intent(out)      :: aerosol_out
integer,                intent(in)       :: index

if( .not. aerosol_initialized ) &
   call error_mesg( 'get_aerosol_single','must first call aerosol_init',FATAL )

!aerosol_out(:,:,:) = aerosol3(:,:,:,index)

end subroutine get_aerosol_single


subroutine get_aerosol_single_name(aerosol_out, name)
!-----------------------------------------------------------------------
! ... Return aerosol fields (for single aerosol species, by name)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------
real, dimension(:,:,:), intent(out)      :: aerosol_out
character(len=*),       intent(in)       :: name

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
integer :: index
logical :: found

if( .not. aerosol_initialized ) &
   call error_mesg( 'get_aerosol_single_name','must first call aerosol_init',FATAL )

found = .false.
do index = 1,nfields
   if( data_names(index) == name ) then
!      aerosol_out(:,:,:) = aerosol3(:,:,:,index)
      found = .true.
   end if
end do

if( .not. found ) then
   call error_mesg( 'get_aerosol_single_name','field '//TRIM(name)//' not found',FATAL )
end if

end subroutine get_aerosol_single_name


!=======================================================================

subroutine get_aerosol_info( names, num_fields )
!-----------------------------------------------------------------------
! ... Return aerosol field information
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------
character(len=*), dimension(:), intent(out), optional :: names
integer,                        intent(out), optional :: num_fields

if( .not. aerosol_initialized ) &
   call error_mesg( 'get_aerosol_info','must first call aerosol_init',FATAL )

if( present(names) )      names(:nfields) = data_names(:nfields)
if( present(num_fields) ) num_fields      = nfields

end subroutine get_aerosol_info

!=======================================================================

subroutine get_aerosol_optical_info( num_categories, nwavenumbers, &
                                     names, wavenumbers, &
                                     aer_ext, aer_ss_alb, aer_asymm)
!-----------------------------------------------------------------------
! ... Return aerosol optical property information
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------
integer,                        intent(out), optional :: num_categories,  &
                                                         nwavenumbers
character(len=*), dimension(:), intent(out), optional :: names
integer, dimension(:),          intent(out), optional :: wavenumbers
real, dimension(:,:),           intent(out), optional :: aer_ext, &
                                                         aer_ss_alb, &
                                                         aer_asymm

if( .not. aerosol_initialized ) &
   call error_mesg( 'get_aerosol_info','must first call aerosol_init',FATAL )

if( present(num_categories) ) num_categories            = NAERMODELS
if( present(nwavenumbers))    nwavenumbers              = num_wavenumbers
if( present(names) )          names(:NAERMODELS)        = aerosol_optical_names(:NAERMODELS)
if( present(wavenumbers) )    wavenumbers(:num_wavenumbers) = endaerwvnsf(:num_wavenumbers)
if( present(aer_ext) )        aer_ext(:,:)              = aeroextivl(:,:)
if( present(aer_ss_alb) )     aer_ss_alb(:,:)           = aerossalbivl(:,:)
if( present(aer_asymm) )      aer_asymm(:,:)            = aeroasymmivl(:,:)

end subroutine get_aerosol_optical_info

!=======================================================================

function get_aerosol_optical_index( naerosol, rh ) result(index)
!-----------------------------------------------------------------------
! ... Return aerosol optical property index for given aerosol number and RH
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------
integer,                        intent(in) :: naerosol
real,                           intent(in) :: rh

!-----------------------------------------------------------------------
! ... Function value
!-----------------------------------------------------------------------
integer                                    :: index

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
integer             :: irh

if( .not. aerosol_initialized ) &
   call error_mesg( 'get_aerosol_optical_index','must first call aerosol_init',FATAL )

if( naerosol > nfields ) then
   call error_mesg( 'get_aerosol_optical_index','aerosol index exceeds number of aerosol fields',FATAL )
end if

if( optical_index(naerosol) == 0 ) then
! ... Sulfate aerosol (RH-dependent)
   irh = MIN( 100, MAX( 0, NINT( 100. * rh ) ) )
   index = sulfate_index( irh )
else
! ... Other aerosols
   index = optical_index(naerosol)
end if

if( index == 0 ) &
   call error_mesg( 'get_aerosol_optical_index', &
                    'Cannot find aerosol optical properties for species = ' // &
                    TRIM( data_names(naerosol) ), &
                    FATAL )

end function get_aerosol_optical_index

!=======================================================================

subroutine aerosol_end
!-----------------------------------------------------------------------
! ... End aerosol interpolation
!-----------------------------------------------------------------------

if( .not. aerosol_initialized ) &
   call error_mesg( 'aerosol_end','must first call aerosol_init',FATAL )

call interpolator_end( aerosol_interp_type )

deallocate (aeroextivl, &
            aerossalbivl, &
            aeroasymmivl, &
            endaerwvnsf )
aerosol_initialized = .false.

!end subroutine interp_aerosol_end
end subroutine aerosol_end

!=======================================================================

subroutine aerosol_alloc (Aerosol, ix, jx, kx)

type(aerosol_type), intent(inout) :: Aerosol
integer, intent(in), optional :: ix, jx, kx

!-----------------------------------------------------------------------
! ... Allocate array for aerosol data
!-----------------------------------------------------------------------

if( .not. aerosol_initialized ) &
   call error_mesg( 'aerosol_alloc','must first call aerosol_init',FATAL )

   if (present (ix) .and. present(jx) .and. present(kx) ) then
     allocate( Aerosol%aerosol( ix, jx, kx,            nfields ) )
   endif 

allocate (Aerosol%aerosol_optical_names (MAX_OPTICAL_FIELDS))
allocate (Aerosol%sulfate_index (0:100          ))
allocate (Aerosol%optical_index(MAX_DATA_FIELDS))

   Aerosol%MAX_DATA_FIELDS = MAX_DATA_FIELDS
   Aerosol%nfields         = nfields         
   Aerosol%aerosol_optical_names = aerosol_optical_names
   Aerosol%sulfate_index   = sulfate_index
   Aerosol%optical_index   = optical_index



end subroutine aerosol_alloc

subroutine aerosol_dealloc
!-----------------------------------------------------------------------
! ... Allocate array for aerosol data
!-----------------------------------------------------------------------

if( .not. aerosol_initialized ) &
   call error_mesg( 'aerosol_dealloc','must first call aerosol_init',FATAL )


end subroutine aerosol_dealloc


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




implicit none

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
integer, parameter :: nlon=20, nlat=10,nlev=8
real :: latb(nlat+1),lonb(nlon+1),pi,phalf(nlon,nlat,nlev+1)
integer :: nlev,i,nspecies
type(time_type) :: model_time
real, allocatable :: aerosol(:,:,:,:)

pi = 4.*atan(1.)

call mpp_init
call mpp_io_init
call mpp_domains_init
call diag_manager_init
call set_calendar_type(JULIAN)

do i = 1,nlat+1
   latb(i) = -90. + 180.*REAL(i-1)/REAL(nlat)
end do
do i = 1,nlon+1
   lonb(i) = -180. + 360.*REAL(i-1)/REAL(nlon)
end do

latb(:) = latb(:) * pi/180.
lonb(:) = lonb(:) * pi/180.

do i = 1,nlev+1
   phalf(:,:,i) = 101325. * REAL(i-1) / REAL(nlev)
end do

!call interp_aerosol_init( lonb, latb, nlev )
call aerosol_init( lonb, latb, nlev )
call aerosol_alloc

call get_aerosol_info( num_fields=nspecies )
ALLOCATE( aerosol(nlon,nlat,nlev,nspecies) )

model_time = set_date(1980,1,1,0,0,0)
!call interp_aerosol(model_time, phalf)
call aerosol_driver(model_time, phalf)

call get_aerosol( aerosol )

call aerosol_dealloc
!call interp_aerosol_end
call aerosol_end
call mpp_exit

end program main

#endif
