
                    module radiation_step_mod


use time_manager_mod,          only:  time_type
use  utilities_mod,            only:  open_file, file_exist,   & 
                                      check_nml_error, error_mesg, &
                                      print_version_number, FATAL, &
				      NOTE, WARNING, get_my_pe, &
				      close_file
use rad_utilities_mod,         only:  Environment, environment_type
use radiation_package_mod,     only:  radpack
use ozone_mod,                 only:  ozone_drvr, ozone_dealloc, &
				      ozone_alloc, ozone_time_vary
use longwave_driver_mod,       only:  longwave_driver_alloc, &
                                      longwave_driver_dealloc
use shortwave_driver_mod,      only:  shortwave_driver_alloc, &
				      shortwave_driver_dealloc, &
				      solar_constant
use radiative_gases_mod,       only:  radiative_gases_time_vary
use cloudrad_package_mod,      only:  clouddrvr, cldmarch,    &
				      cloudrad_package_alloc, &
				      cloudrad_package_dealloc
use surface_albedo_mod,        only:  sfcalbedo, sfcalbedo_alloc, &
				      sfcalbedo_dealloc
use donner_ice_mod,            only:  inquire_donner_ice, &
				      donner_iceclouds_rad_output

!-------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!       module which executes physics that is executed only on a 
!                        radiation time step
!
!--------------------------------------------------------------------
  



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!    character(len=5), parameter  ::  version_number = 'v0.09'
     character(len=128)  :: version =  '$Id: radiation_step.F90,v 1.2 2001/08/30 15:12:57 fms Exp $'
     character(len=128)  :: tag     =  '$Name: galway $'



!---------------------------------------------------------------------
!-------  interfaces --------

public        &
        radiation_step_init, radiation_step_time_vary, &
	radiation_step_dr,   radiation_step_alloc,    &
	radiation_step_dealloc


!--------------------------------------------------------------------
!----    namelist -----

logical :: do_donner_cloud_forcing=.false.


namelist / radiation_step_nml /         &
                                 do_donner_cloud_forcing

!--------------------------------------------------------------------
!-------- public data  -----



!--------------------------------------------------------------------
!------ private data ------




!-------------------------------------------------------------------
!-------------------------------------------------------------------




contains




subroutine radiation_step_init

!------------------------------------------------------------------
      integer  :: unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=radiation_step_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'radiation_step_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'radiation_step', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=radiation_step_nml)
      endif
      call close_file (unit)

!--------------------------------------------------------------------

end subroutine radiation_step_init

!####################################################################

subroutine radiation_step_time_vary (Time)

type(time_type), intent(in), optional  :: Time
! Time will be present when running in FMS

      call solar_constant
      call radiative_gases_time_vary
      call cldmarch
      if (present (Time) ) then
        call ozone_time_vary (Time)
      else
        call ozone_time_vary
      endif


end subroutine radiation_step_time_vary




!####################################################################

subroutine radiation_step_dr (is, ie, js, je, Time_next, kbot, mask) 

integer,                 intent(in)           :: is, ie, js, je
type(time_type),         intent(in), optional :: Time_next
integer, dimension(:,:), intent(in), optional :: kbot
real, dimension(:,:,:),  intent(in), optional :: mask

!---------------------------------------------------------------------
!     radiation_step_dr calls the routines needed on 
!     radiation time steps.
!----------------------------------------------------------------------

!--------------------------------------------------------------------
      logical :: do_donner_ice

!--------------------------------------------------------------------
!  call sfcalbedo to obtain surface albedo information for the rad-
!  iation calculation
!--------------------------------------------------------------------
      call sfcalbedo

!--------------------------------------------------------------------
!  call clouddrvr to obtain cloud information for the radiation cal-
!  culation
!--------------------------------------------------------------------
      if (present(kbot) ) then
    call clouddrvr ( is, ie, js, je, Time_next=Time_next, kbot=kbot,   &
			 mask=mask)
      else
        call clouddrvr (is, ie, js, je, Time_next=Time_next)
      endif

!--------------------------------------------------------------------
!  call ozonedrvr to define the (i,j,k) ozone field to be used for the
!  radiation calculation
!--------------------------------------------------------------------
      call ozone_drvr

!--------------------------------------------------------------------
!  call radpack to perform the radiation calculation
!--------------------------------------------------------------------
      call radpack  

!--------------------------------------------------------------------
!  if donner_ice is active, save some radiation package outputs.
!--------------------------------------------------------------------
      call inquire_donner_ice (do_donner_ice)
      if (do_donner_ice) then
        call donner_iceclouds_rad_output
      endif

!--------------------------------------------------------------------
   

end subroutine radiation_step_dr


!####################################################################

subroutine radiation_step_alloc

!---------------------------------------------------------------------

     call sfcalbedo_alloc
     call cloudrad_package_alloc
     call ozone_alloc
     call shortwave_driver_alloc
     call longwave_driver_alloc

!---------------------------------------------------------------------

end subroutine radiation_step_alloc




!####################################################################

subroutine radiation_step_dealloc

!---------------------------------------------------------------------

     call longwave_driver_dealloc
     call shortwave_driver_dealloc
     call ozone_dealloc
     call cloudrad_package_dealloc
     call sfcalbedo_dealloc

!---------------------------------------------------------------------

end subroutine radiation_step_dealloc



!#####################################################################



	        end module radiation_step_mod

