

              module   co2_mod


use  utilities_mod,     only:  open_file, file_exist,    &
			       check_nml_error, error_mesg, &
			       print_version_number, FATAL, NOTE, &
			       WARNING, get_my_pe, close_file, &
			       read_data, write_data
!use  constants_mod, only:  wtmair
use lw_gases_stdtf_mod, only:  co2_lblinterp

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                  interface module to calculate co2 
!                      transmission functions
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

      character(len=128)  :: version =  '$Id: co2.F90,v 1.4 2003/04/09 20:59:02 fms Exp $'
      character(len=128)  :: tag     =  '$Name: inchon $'


!---------------------------------------------------------------------
!-------  interfaces --------

public      co2_init, co2_time_vary

!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=8)      :: dummy    = '     '



namelist /co2_nml/                      &
                         dummy  



!---------------------------------------------------------------------
!------- public data ------




!---------------------------------------------------------------------
!------- private data ------

real       ::  rco2air
real                ::    rrco2



!---------------------------------------------------------------------
!---------------------------------------------------------------------




                        contains


subroutine co2_init 

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!    molecular weight of co2
!---------------------------------------------------------------------
     real                 ::  wtmco2  = 44.00995


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
     if (get_my_pe() == 0) then 
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit, nml=co2_nml)
     endif
     call close_file (unit)

!---------------------------------------------------------------------
!   define co2 mixing ratio conversion factor.
!   (an additional fms mimicking change, if desired)
!---------------------------------------------------------------------
!      rco2air = wtmco2/wtmair

end subroutine co2_init



!####################################################################

subroutine co2_time_vary ( rrvco2            )

!---------------------------------------------------------------------
real, intent(in   )    ::  rrvco2

!---------------------------------------------------------------------
        real    ::    co2_vmr

        rrco2 = rrvco2*rco2air
        co2_vmr = rrvco2*1.0E+06

        call co2_lblinterp  (co2_vmr            )



end subroutine co2_time_vary



!###################################################################


		  end module co2_mod
