
		 module longwave_params_mod


use utilities_mod,      only:  open_file, file_exist,    &
                               check_nml_error, error_mesg, &
                               print_version_number, FATAL, NOTE, &
			       WARNING, get_my_pe, close_file


!--------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!          module containing parameters for longwave code
!
!------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.08'
    character(len=128)  :: version =  '$Id: longwave_params.F90,v 1.3 2002/07/16 22:35:43 fms Exp $'
    character(len=128)  :: tag     =  '$Name: havana $'

!--------------------------------------------------------------------
!----- interfaces ------

public longwave_params_init

!---------------------------------------------------------------------
!-------- namelist  ---------

!!3 character(len=5)        :: dummy  = '     '
character(len=8)        :: dummy  = '     '


namelist /longwave_params_nml/    &
                               dummy                 



!-------------------------------------------------------------------
!----- public data --------

integer, parameter, public   :: NBCO215     = 3
integer, parameter, public   :: NBLY_ORIG   = 16
integer, parameter, public   :: NBLY_CKD2P1 = 48
integer, parameter, public   :: NBLW        = 300
integer, parameter, public   :: NBLX        = 48




!-------------------------------------------------------------------
!----- private data --------




!---------------------------------------------------------------------
!---------------------------------------------------------------------



contains



subroutine longwave_params_init

!------------------------------------------------------------------
     integer    ::  unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=longwave_params_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'longwave_params_nml')
       enddo
10     call close_file (unit)
     endif

     unit = open_file ('logfile.out', action='append')
!    call print_version_number (unit, 'longwave_params', version_number)
     if (get_my_pe() == 0) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=longwave_params_nml)
       write (unit,9000) NBCO215, NBLY_ORIG, NBLY_CKD2P1, NBLW, NBLX 
     endif
     call close_file (unit)


9000 format ( 'NBCO215=', i3,'  NBLY_ORIG=', i4,   &
              '  NBLY_CKD2P1=', i4, '  NBLW= ', i4, '  NBLX=', i4 )
!------------------------------------------------------------------

end  subroutine longwave_params_init



!####################################################################



	      end module longwave_params_mod
