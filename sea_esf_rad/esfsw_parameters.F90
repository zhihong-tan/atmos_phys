
		 module esfsw_parameters_mod


use utilities_mod,      only:  open_file, file_exist,    &
                               check_nml_error, error_mesg, &
                               print_version_number, FATAL, NOTE, &
			       WARNING, get_my_pe, close_file


!--------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!          module containing parameters for esf shortwave code,
!               including solar flux band information
!
!------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.08'
    character(len=128)  :: version =  '$Id: esfsw_parameters.F90,v 1.3 2002/07/16 22:35:06 fms Exp $'
    character(len=128)  :: tag     =  '$Name: inchon $'

!--------------------------------------------------------------------
!----- interfaces ------

public esfsw_parameters_init,   &
       put_solarfluxes, get_solarfluxes

!---------------------------------------------------------------------
!-------- namelist  ---------

!!3 character(len=12)         :: sw_resolution = ' '
character(len=16)         :: sw_resolution = ' '
integer                   :: sw_diff_streams = 0


namelist /esfsw_parameters_nml/    &
                                 sw_resolution, sw_diff_streams



!-------------------------------------------------------------------
!----- public data --------

integer, public   :: nbands, nfrqpts, nstreams, nh2obands

integer, parameter, public   :: TOT_WVNUMS  = 57600

!-------------------------------------------------------------------
!----- private data --------

real,    dimension(TOT_WVNUMS)        :: solarfluxtoa

real,    dimension(:), allocatable     :: solflxband
integer, dimension(:), allocatable     :: endwvnbands


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                      contains



subroutine esfsw_parameters_init

!------------------------------------------------------------------
     integer    ::  unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=esfsw_parameters_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'esfsw_parameters_nml')
       enddo
10     call close_file (unit)
     endif

     if (trim(sw_resolution) == 'high') then
       nbands = 25
       nfrqpts = 72
       nh2obands = 14
     elseif (trim(sw_resolution) == 'low') then
       nbands = 18
       nfrqpts = 38
       nh2obands = 9
     else
       call error_mesg ( 'esfsw_parameters_init',   &
	    ' sw_resolution must be specified as "high" or "low".', &
	    FATAL)
     endif

     if (sw_diff_streams == 4) then
       nstreams = 4
     else if (sw_diff_streams == 1) then
       nstreams = 1
     else
       call error_mesg ( 'esfsw_parameters_init',   &
          ' sw_diff_streams must be specified as either 1 or 4.', FATAL)
     endif

     unit = open_file ('logfile.out', action='append')
!    call print_version_number (unit, 'esfsw_parameters',  &
!							version_number)
     if (get_my_pe() == 0) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,9000) NBANDS, NFRQPTS, NSTREAMS, NH2OBANDS 
       write (unit,nml=esfsw_parameters_nml)
     endif
     call close_file (unit)

9000   format ( 'NBANDS=  ', i4, '  NFRQPTS=', i4, &
		'NSTREAMS= ', i4, 'NH2OBANDS= ', i4 )

     allocate (solflxband (nbands) )
     allocate (endwvnbands (0:nbands) )

!------------------------------------------------------------------

end subroutine esfsw_parameters_init



!####################################################################

subroutine put_solarfluxes (solarfluxtoa_in, solflxband_in,     &
			    endwvnbands_in)

!----------------------------------------------------------------------
real,    dimension(:), intent(in)  :: solarfluxtoa_in
real,    dimension(:),     intent(in)  :: solflxband_in
integer, dimension(:),   intent(in)  :: endwvnbands_in

!----------------------------------------------------------------------

    solarfluxtoa = solarfluxtoa_in
    solflxband = solflxband_in
    endwvnbands = endwvnbands_in


end subroutine put_solarfluxes 


!###################################################################

subroutine get_solarfluxes (solarfluxtoa_out, solflxband_out,     &
			    endwvnbands_out)

!----------------------------------------------------------------------
real,    dimension(:), intent(out)  :: solarfluxtoa_out
real,    dimension(:),     intent(out)  :: solflxband_out
integer, dimension(:),   intent(out)  :: endwvnbands_out
!----------------------------------------------------------------------

    solarfluxtoa_out = solarfluxtoa
    solflxband_out = solflxband
    endwvnbands_out = endwvnbands


end subroutine get_solarfluxes 


!###################################################################

	      end module esfsw_parameters_mod
