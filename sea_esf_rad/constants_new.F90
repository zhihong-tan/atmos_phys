
                   module constants_new_mod


use utilities_mod,      only: open_file, file_exist,    &
                              check_nml_error, error_mesg, &
                              print_version_number, FATAL, NOTE, &
			      WARNING, get_my_pe, close_file, &
			      read_data, write_data
use rad_utilities_mod,  only: Environment, environment_type

!--------------------------------------------------------------------
 
implicit none
private

!------------------------------------------------------------------
!                   module with constants
!
!-------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!     character(len=5), parameter  ::  version_number = 'v0.08'
      character(len=128)  :: version =  '$Id: constants_new.F90,v 1.3 2001/10/25 17:48:25 fms Exp $'
      character(len=128)  :: tag     =  '$Name: fez $'



!---------------------------------------------------------------------
!-------  interfaces --------

public  constants_new_init, define_constants

!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=5)   :: dummy    = '     '


namelist /constants_new_nml/   &
			      dummy





!--------------------------------------------------------------------
!----- public ---------------


!--------------------------------------------------------------------
!the following constants are currently used in the longwave code:
!--------------------------------------------------------------------

real, public, parameter :: wtmair = 2.896440E+01
real, public, parameter :: wtmh2o = 1.801534E+01
real, public, parameter :: diffac = 1.660000E+00
real, public, parameter :: rh2oair  = wtmh2o/wtmair


real, public, parameter :: secday  = 8.640000E+04
real, public, parameter :: cpair   = 1.004840E+07
real, public, parameter :: grav    = 9.806650E+02


!!!! USE THIS VALUE FOR CONSISTENCY WITH VALUE IN Constants_mod, 
!!!! beginning with fez release:
real, public, parameter :: sigma    = 5.6734E-05 
!real, public, parameter :: sigma   = 5.670E-05

real, public, parameter :: avogno = 6.023000E+23
!--------------------------------------------------------------------
!  the following value may be a newer value to use  (ifdef newavog)?? : 
!real, public, parameter :: avogno = 6.022045E+23
!--------------------------------------------------------------------

real, public, parameter :: rgas = 8.314320E+07
!--------------------------------------------------------------------
!  the following value may be a newer value to use  (ifdef newrgas)?? : 
!real, public, parameter :: rgas = 8.314410E+07
!--------------------------------------------------------------------

real, public, parameter :: pstd    = 1.013250E+06
real, public, parameter :: rearth  = 6.356766E+08  ! earth radius in cm

real, public            :: radcon

!------------------------------------------------------------------
!the following constants are currently used in the esf shortwave code:
!---------------------------------------------------------------------

real, public, parameter :: o2mixrat    = 2.0953E-01
real, public, parameter :: gasconst    = 287.053    
real, public, parameter :: rhoair      = 1.292269  
real, public, parameter :: grav_mks    = 9.80665   
real, public, parameter :: cpair_mks   = 1.004840E+03
real, public, parameter :: radcon_mks  = (grav_mks/cpair_mks)*secday
real, public, parameter :: pstd_mks    = 101325.0


real, public, parameter :: alogmin     = -50.0      

real, public, parameter :: frezdk      = 273.16

real, public            :: pie, radians_to_degrees

!--------------------------------------------------------------------
!the following constants were defined previously in packconst, 
!though may not currently be used in the radiation code:
!---------------------------------------------------------------------

real, public, parameter :: calyear     = 365.2500

!real, public, parameter :: sigmasb     = 5.670320E-05 ! (ifdef newsigma)
!!!! USE THIS VALUE FOR CONSISTENCY WITH VALUE IN Constants_mod, 
!!!! beginning with fez release:
real, public, parameter :: sigmasb     = 5.6734E-05 
!--------------------------------------------------------------------
!real, public, parameter :: sigmasb     = 5.673000E-05 ! (older value) 
!--------------------------------------------------------------------

real, public, parameter :: boltz       = 1.380662E-23    ! [Joules/K]
real, public, parameter :: amu         = 1.66E-27        ! [kg]
real, public, parameter :: clight      = 2.99792458E+10
real, public, parameter :: planck      = 6.626176E-27  
real, public, parameter :: secmin      = 60.0000  
real, public, parameter :: sechour     = 3600.00       
real, public, parameter :: wtmo3       = 47.99820E+01

integer, public,parameter :: bytes_per_word = 8

!--------------------------------------------------------------------
! ---  private ---

real,         parameter :: radcon_mgroup  = (grav/cpair)*secday
!!!  for FMS MIMIC :
real,         parameter :: radcon_fms  = 8.427 ! instead of 8.4321 
 

!--------------------------------------------------------------------
!--------------------------------------------------------------------


                         contains



subroutine constants_new_init

!------------------------------------------------------------------
     integer    ::  unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=constants_new_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'constants_new_nml')
       enddo
10     call close_file (unit)
     endif

     unit = open_file ('logfile.out', action = 'append')
!    call print_version_number(unit, 'constants_new', version_number)
     if (get_my_pe() == 0)  then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=constants_new_nml)
     endif
     call close_file (unit)


     pie   = 4.0E+00*ATAN(1.0E+00)
     radians_to_degrees = 180.0/pie


!------------------------------------------------------------------



end subroutine constants_new_init


!####################################################################

subroutine define_constants

    if (Environment%using_fms_periphs) then
      radcon = radcon_fms
    else if (Environment%using_sky_periphs) then
      radcon = radcon_mgroup
    endif

end subroutine define_constants



                 end module constants_new_mod
