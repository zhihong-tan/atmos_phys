
                  module radiative_gases_mod


use rad_utilities_mod,       only: longwave_control_type, Lw_control
use utilities_mod,           only: open_file, file_exist,    &
			           check_nml_error, error_mesg, &
			           print_version_number, FATAL, NOTE, &
				   WARNING, get_my_pe, close_file 
use  ch4_n2o_mod,            only: ch4_n2o_init, ch4_n2o_time_vary
use  cfc_mod,                only: cfc_init, cfc_time_vary
use co2_mod,                 only: co2_init, co2_time_vary
use radiation_diag_mod,      only: radiag_from_radgases
use rad_output_file_mod,     only: hold_gases
use gas_tf_mod,              only: get_control_gas_tf
use lw_gases_stdtf_mod,      only: lw_gases_stdtf_time_vary, &
				   lw_gases_stdtf_dealloc
use longwave_setup_mod,      only: longwave_parameter_type,    &
				   Lw_parameters
use constants_new_mod,       only: bytes_per_word

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                  radiative gases module
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128)  :: version =  '$Id: radiative_gases.F90,v 1.2 2001/07/05 17:32:30 fms Exp $'
  character(len=128)  :: tag     =  '$Name: fez $'

!---------------------------------------------------------------------
!-------  interfaces --------

public   radiative_gases_init,      &
	 radiative_gases_time_vary, &
	 radiative_gases_end


private  read_restart_radiative_gases, &
         write_restart_radiative_gases


!---------------------------------------------------------------------
!-------- namelist  ---------

logical              :: direct_access_read  = .false.
logical              :: direct_access_write = .false.
logical              :: do_ch4_n2o          = .true.
logical              :: do_f11              = .true.
logical              :: do_f12              = .true.
logical              :: do_f113             = .true.
logical              :: do_f22              = .true.
logical              :: do_co2              = .true.
integer              :: initial_record_number = 1
character(len=16)    :: ch4n2o_data_source = '   '
character(len=16)    :: cfc_data_source = '   '
character(len=16)    :: co2_data_source = '   '


namelist /radiative_gases_nml/                      &
				do_ch4_n2o,      &
				do_f11, do_f12, do_f113, do_f22, &
				do_co2,            &
                                direct_access_read, &
                                direct_access_write,   &
				initial_record_number, &
				ch4n2o_data_source, &
				cfc_data_source, &
				co2_data_source


!---------------------------------------------------------------------
!------- public data ------

real, public            :: rrvco2, rrvf11, rrvf12, rrvf113,   &
                           rrvf22, rrvch4, rrvn2o     


!---------------------------------------------------------------------
!------- private data ------

integer,parameter        ::  ngasses=7                  
real                     ::  rrvco2_rst, rrvf11_rst, rrvf12_rst, &
			     rrvf113_rst, rrvf22_rst, rrvch4_rst, &
			     rrvn2o_rst

integer, dimension(1)    ::  restart_versions = (/ 1 /)

!---------------------------------------------------------------------
!---------------------------------------------------------------------




contains


!####################################################################

subroutine radiative_gases_init

!---------------------------------------------------------------------
     integer                       :: unit, ierr, io
     logical                       :: do_cfc

!---------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=radiative_gases_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'radiative_gases_nml')
       enddo
10     call close_file (unit)
     endif

     unit = open_file ('logfile.out', action='append')
!    call print_version_number (unit, 'radiative_gases', version_number)
     if (get_my_pe() == 0) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=radiative_gases_nml)
     endif
     call close_file (unit)

     if (do_f11 .or. do_f12 .or. do_f113 .or. do_f22) then
       do_cfc=.true.
     else
       do_cfc=.false.
     endif

     Lw_control%do_cfc = do_cfc
     Lw_control%do_ch4_n2o = do_ch4_n2o
     Lw_control%do_co2 = do_co2

!---------------------------------------------------------------------
!   read the radiative gases restart file, if present. if restart
!   files are not present values for gas mixing ratios will be set 
!   negative, and the model will abort below.
!---------------------------------------------------------------------
     call read_restart_radiative_gases

!---------------------------------------------------------------------
!   initialize constituent gases,  if they are to be activated.
!---------------------------------------------------------------------
     if (do_ch4_n2o) then
       if (trim(ch4n2o_data_source) /= 'restart' ) then
         call ch4_n2o_init (ch4n2o_data_source)
       else
	 if (rrvch4_rst >= 0.0) then
	   rrvch4 = rrvch4_rst
	   rrvn2o = rrvn2o_rst
	 else
	   call error_mesg ('radiative_gases_init', &
           'cannot use restart ch4n2o values without a restart file', &
							     FATAL)
	 endif
         call ch4_n2o_init (ch4n2o_data_source, rrvch4, rrvn2o)
       endif
       Lw_parameters%NBTRG  = 1
       Lw_parameters%NBTRGE = 1
     else
       rrvch4 = 0.0
       rrvn2o = 0.0
       Lw_parameters%NBTRG  = 0
       Lw_parameters%NBTRGE = 0
     endif 
 
     if (Lw_control%do_cfc) then
       if (trim(cfc_data_source) /= 'restart' ) then
         call cfc_init (cfc_data_source, do_f11, do_f12, do_f113,   &
			do_f22)
       else
	 if (rrvf11_rst >= 0.0) then
	   rrvf11  = rrvf11_rst
	   rrvf12  = rrvf12_rst
	   rrvf113 = rrvf113_rst
	   rrvf22  = rrvf22_rst
	 else
	   call error_mesg ('radiative_gases_init', &
             'cannot use restart cfc values without a restart file', &
							    FATAL)
	 endif
         call cfc_init (cfc_data_source, do_f11, do_f12, do_f113,   &
			do_f22, rrvf11, rrvf12, rrvf113, rrvf22)
       endif
     else
       rrvf11  = 0.0
       rrvf12  = 0.0
       rrvf113 = 0.0
       rrvf22  = 0.0
     endif

     if (do_co2) then
       if (trim(co2_data_source) /= 'restart' ) then
         call co2_init (co2_data_source)
       else
	 if (rrvco2_rst >= 0.0) then
	   rrvco2 = rrvco2_rst
	 else
	   call error_mesg ('radiative_gases_init', &
             'cannot use restart co2 values without a restart file', &  
							     FATAL)
	 endif
         call co2_init (co2_data_source, rrvco2)
       endif
     else 
       rrvco2 = 0.0
     endif

!---------------------------------------------------------------------


end subroutine radiative_gases_init



!####################################################################

subroutine radiative_gases_time_vary 

!---------------------------------------------------------------------
      logical             ::  calc_co2, calc_n2o, calc_ch4
!---------------------------------------------------------------------

      call get_control_gas_tf (calc_co2, calc_ch4, calc_n2o)

      call lw_gases_stdtf_time_vary

      if (do_ch4_n2o) then
	call ch4_n2o_time_vary (rrvch4, rrvn2o)
      endif

      if (Lw_control%do_cfc) then
	call cfc_time_vary (rrvf11, rrvf12, rrvf113, rrvf22, &
			    do_f11, do_f12, do_f113, do_f22)
      endif

      if (do_co2) then
	call co2_time_vary (rrvco2)
      endif

      call radiag_from_radgases (rrvf11, rrvf12, rrvf113, rrvf22, &
	  		         rrvch4, rrvn2o, rrvco2)
      
      call hold_gases (rrvco2, rrvch4, rrvn2o, rrvf11, rrvf12,   &
		       rrvf113, rrvf22)
		       
      if (calc_co2 .or. calc_ch4 .or. calc_n2o) then
	call lw_gases_stdtf_dealloc
      endif


end subroutine radiative_gases_time_vary


!####################################################################


subroutine radiative_gases_end


    call write_restart_radiative_gases



end subroutine radiative_gases_end


!####################################################################

subroutine write_restart_radiative_gases

!---------------------------------------------------------------------
       integer                 ::  ierr, io, unit
       character(len=64)       ::  name
       integer                 ::  recnum

!---------------------------------------------------------------------
!    open unit and write radiative gas restart parameters file
!---------------------------------------------------------------------
       if (get_my_pe() == 0) then
       name = 'RESTART/radiative_gases.parameters.res'
       unit = open_file (name, form= 'unformatted',     &
			 action= 'write')
       write (unit) restart_versions(size(restart_versions)), ngasses
       call close_file (unit)
       endif

!---------------------------------------------------------------------
!    open unit and write radiative gas restart file
!---------------------------------------------------------------------
       if (get_my_pe() == 0) then
         name = 'RESTART/radiative_gases.res'
         if (direct_access_write) then
           unit = open_file (name, form= 'unformatted',  &
			     action= 'write', access= 'direct',   &
                             recl= ngasses*bytes_per_word)
           write (unit, rec=1) rrvco2, rrvf11, rrvf12, rrvf113,   &
			       rrvf22, rrvch4, rrvn2o
         else
           unit = open_file (name, form= 'unformatted', action= 'write')
           write (unit) rrvco2
           write (unit) rrvf11, rrvf12, rrvf113, rrvf22
           write (unit) rrvch4, rrvn2o
         endif
         call close_file (unit)
       endif
!----------------------------------------------------------------------
      


end subroutine write_restart_radiative_gases


!####################################################################

subroutine read_restart_radiative_gases 

!---------------------------------------------------------------------

      character(len=64)        ::  name
      integer                 ::  vers 
      integer                  ::  radiative_gas_reclen=-1  ! in words
      integer                  ::  unit
      character(len=4)        ::  chvers

!---------------------------------------------------------------------
!    open unit to read radiative gas restart parameters file
!---------------------------------------------------------------------
     name = 'INPUT/radiative_gases.parameters.res'

     if (file_exist(name) ) then
       unit = open_file (name, form= 'unformatted',     &
		         action= 'read')
!----------------------------------------------------------------------
! read and check restart version number and data record length.
!----------------------------------------------------------------------
       read (unit) vers, radiative_gas_reclen
       if ( .not. any(vers == restart_versions) ) then
         write (chvers,'(i4)') vers
         call error_mesg ('sea_esf_rad_init', &
                      'restart version '//chvers//' cannot be read &
                       &by this module version', FATAL)
       endif
       if (get_my_pe() == 0) write (*,9000)     vers
9000 format (/, 'restart_version for radiative_gas_restart_file =  ',  &
	         i5)
       call close_file (unit)
     endif

!---------------------------------------------------------------------
!    open unit to read radiative gas restart file
!    read radiative gas restart file
!---------------------------------------------------------------------
     name = 'INPUT/radiative_gases.res'
     if (file_exist(name) ) then
       if (direct_access_read) then
	 if (radiative_gas_reclen > 0) then
           unit = open_file (name, form= 'unformatted',     &
	                     action= 'read',         &
	                     access= 'direct',      &
                             recl= radiative_gas_reclen*bytes_per_word)
           read (unit, rec=     1) rrvco2_rst, rrvf11_rst, rrvf12_rst, &
		  		   rrvf113_rst, rrvf22_rst,  &
	               	           rrvch4_rst, rrvn2o_rst
	 else
	   call error_mesg ('read_restart_radiative_gases',  &
	    	          'negative value for radiative_gas_reclen', &
							       FATAL)
	 endif
       else
         unit = open_file (name, form= 'unformatted',     &
 			   action= 'read')
         read (unit) rrvco2_rst
         read (unit) rrvf11_rst, rrvf12_rst, rrvf113_rst, rrvf22_rst
         read (unit) rrvch4_rst, rrvn2o_rst
       endif
       call close_file (unit)
     else
       rrvco2_rst = -1.0
       rrvf11_rst = -1.0
       rrvf12_rst = -1.0
       rrvf113_rst = -1.0
       rrvf22_rst = -1.0
       rrvch4_rst = -1.0
       rrvn2o_rst = -1.0
     endif

!---------------------------------------------------------------------


end subroutine read_restart_radiative_gases


!###################################################################




		  end module radiative_gases_mod
