
              module   co2_mod


use  utilities_mod,     only:  open_file, file_exist,    &
			       check_nml_error, error_mesg, &
			       print_version_number, FATAL, NOTE, &
			       WARNING, get_my_pe, close_file, &
			       read_data, write_data
use  constants_new_mod, only:  wtmair
use lw_gases_stdtf_mod, only:  co2_lblinterp
use rad_utilities_mod,  only:  Environment, environment_type

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                      co2 gas module
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!     character(len=5), parameter  ::  version_number = 'v0.09'
      character(len=128)  :: version =  '$Id: co2.F90,v 1.2 2001/07/05 17:28:41 fms Exp $'
      character(len=128)  :: tag     =  '$Name: galway $'


!---------------------------------------------------------------------
!-------  interfaces --------

public      co2_init, co2_time_vary

!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=5)      :: co2_amt  = 'FIXED'


namelist /co2_nml/                      &
                         co2_amt



!---------------------------------------------------------------------
!------- public data ------


real, public        ::    rrco2


!---------------------------------------------------------------------
!------- private data ------

logical    ::  do_co2_init = .true.
real       ::  rco2air, rco2



!---------------------------------------------------------------------
!---------------------------------------------------------------------




                        contains


subroutine co2_init (data_source, rrvco2_in) 

!---------------------------------------------------------------------
real, intent(in), optional :: rrvco2_in
character(len=*)    ::  data_source

!---------------------------------------------------------------------
!    molecular weight of co2
!---------------------------------------------------------------------
     real                 ::  wtmco2  = 44.00995

!---------------------------------------------------------------------
! optional supplied initial trace gas volume mixing ratios in (no./no.)
!---------------------------------------------------------------------
     real                 ::  rco2_icrccm   = 3.00000E-04

     real                 ::  rco2_ipcc_92  = 3.56000E-04
     real                 ::  rco2_ipcc_80  = 3.37320E-04
     real                 ::  rco2_ipcc_98  = 3.69400E-04

     real                 ::  rco2_330ppm   = 3.30000E-04
     real                 ::  rco2_660ppm   = 6.60000E-04

!---------------------------------------------------------------------
     integer              :: unit, ierr, io, inrad

!---------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=co2_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'co2_nml')
       enddo
10     call close_file (unit)
     endif

     unit = open_file ('logfile.out', action='append')
!    call print_version_number (unit, 'co2', version_number)
     if (get_my_pe() == 0) then 
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit, nml=co2_nml)
     endif
     call close_file (unit)

!---------------------------------------------------------------------
!   define co2 mixing ratio conversion factor.
!   (an additional fms mimicking change, if desired)
!---------------------------------------------------------------------
!   if (Environment%using_fms_periphs) then
!     rco2air = 1.519449738
!   else if (Environment%using_sky_periphs) then
      rco2air = wtmco2/wtmair
!   endif

!---------------------------------------------------------------------
!   define initial co2 volume mixing ratio to be used.
!---------------------------------------------------------------------

    if (trim(data_source) == 'icrccm') then
      rco2   = rco2_icrccm
    else if (trim(data_source) == 'ipcc_92') then
      rco2   = rco2_ipcc_92
    else if (trim(data_source) == 'ipcc_98') then
      rco2   = rco2_ipcc_98
    else if (trim(data_source) == 'ipcc_80') then
      rco2   = rco2_ipcc_80
    else if (trim(data_source) == '330ppm') then
      rco2   = rco2_330ppm
    else if (trim(data_source) == '660ppm') then
      rco2   = rco2_660ppm
    else if (trim(data_source) == 'input') then
      if (file_exist ('INPUT/id1co2') ) then
        inrad = open_file ('INPUT/id1co2', form= 'formatted', &
                            action= 'read')
        read (inrad, FMT = '(5e18.10)')  rco2
        call close_file (inrad)
      else
	call error_mesg ( 'co2_init', &
	          'desired co2 input file is not present', FATAL)
      endif
    else if (trim(data_source) == 'restart') then
       rco2 = rrvco2_in
    else 
      call error_mesg ('co2_init', &
              'no valid data source is specified for co2 input', FATAL)
    endif
	


end subroutine co2_init



!####################################################################

subroutine co2_time_vary ( rrvco2)

!---------------------------------------------------------------------
real, intent(inout)    ::  rrvco2

!---------------------------------------------------------------------
        real    ::    co2_vmr

!---------------------------------------------------------------------
!   the co2 volume mixing ratio is set to the initial value (rco2) and
!   the mass mixing ratio is defined on the first access of this 
!   routine. after first access, this routine does nothing.
!---------------------------------------------------------------------
        if (trim(co2_amt) == 'FIXED' ) then
	  if (do_co2_init) then
            rrvco2 = rco2
            rrco2  = rrvco2*rco2air
            co2_vmr = rrvco2*1.0E+06
            call co2_lblinterp  (co2_vmr, co2_amt)
            do_co2_init = .false.
          endif

!---------------------------------------------------------------------
!  NOTE: TIME VARIATION OF RADIATIVE GASES NOT CURRENTLY AVAILABLE. 
!---------------------------------------------------------------------
        else if (trim(co2_amt) == 'VARY') then
!---------------------------------------------------------------------
!    this is where the time variation of co2 will be added 
!    define rrvco2, the volume mixing ratio, a function of rco2 and
!    time, and then convert it into  rrco2, the mass mixing ratio.
!         rrvco2 = ?????
!         rrco2  = rrvco2*rco2air
!    also then calculate new co2 lw transmission functions
!         co2_vmr = rrvco2*1.0E+06
!         call co2_lblinterp  (co2_vmr, co2_amt)
!---------------------------------------------------------------------
          call error_mesg ( 'co2_time_vary',   &
	      'time-varying co2 not yet implemented', FATAL)
	else
          call error_mesg ( 'co2_time_vary',   &
	      'co2_amt has unacceptable value', FATAL)
        endif


end subroutine co2_time_vary



!###################################################################


		  end module co2_mod
