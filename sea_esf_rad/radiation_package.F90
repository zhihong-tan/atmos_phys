
		  module radiation_package_mod


use optical_path_mod,         only: optical_path_setup, optical_dealloc
use  utilities_mod,           only: open_file, file_exist,     &
                                    check_nml_error, error_mesg, &
                                    print_version_number, FATAL, NOTE, &
				    WARNING, get_my_pe, close_file, &
				    get_num_pes
use rad_utilities_mod,        only: Environment, environment_type
use radiation_diag_mod,       only: radiag_driver
use longwave_driver_mod,      only: lwrad
use shortwave_driver_mod,     only: shortwave
use longwave_setup_mod,       only: longwave_setup_driver,   &
                                    longwave_setup_dealloc
use cool_to_space_mod,        only: cool_to_space_alloc,  &
                                    cool_to_space_dealloc  

!------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!                    radiation package driver module
!
!--------------------------------------------------------------------



!--------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.09'
    character(len=128)  :: version =  '$Id: radiation_package.F90,v 1.2 2001/08/30 15:15:08 fms Exp $'
    character(len=128)  :: tag     =  '$Name: fez $'



!---------------------------------------------------------------------
!-------  interfaces --------

public  radpack_init, radpack   



!-------------------------------------------------------------------
!-------- namelist ----------


real     :: dummy = 0.0

namelist /radiation_package_nml/    &
                                   dummy


!--------------------------------------------------------------------
!------ public data -------





!--------------------------------------------------------------------
!------ private data -------
  




!--------------------------------------------------------------------
!--------------------------------------------------------------------

        contains


subroutine radpack_init


      integer :: unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=radiation_package_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'radiation_package_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'radiation_package',   &
!					       version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=radiation_package_nml)
      endif
      call close_file (unit)


end subroutine radpack_init


!####################################################################

subroutine radpack 

!-----------------------------------------------------------------------
!     radpack calls the routines which calculate the longwave and short-
!     wave radiational heating terms and fluxes and the routine radiag
!     which provides radiation package diagnostics.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------

!---------------------------------------------------------------------- 
      call longwave_setup_driver 
      call cool_to_space_alloc 
      call optical_path_setup 

!----------------------------------------------------------------------
!     compute longwave radiation.
!----------------------------------------------------------------------
      call lwrad 

!----------------------------------------------------------------------
!     compute shortwave radiation.
!----------------------------------------------------------------------
      call shortwave 

!--------------------------------------------------------------------
!   call radiag to compute radiation diagnostics at desired points
!--------------------------------------------------------------------
      call radiag_driver 

!--------------------------------------------------------------------
!   deallocate various arrays
!--------------------------------------------------------------------
      call optical_dealloc
      call cool_to_space_dealloc
      call longwave_setup_dealloc

!---------------------------------------------------------------------


end subroutine radpack


!###################################################################

                      end module radiation_package_mod

