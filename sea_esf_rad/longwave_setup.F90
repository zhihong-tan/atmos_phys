
		 module longwave_setup_mod


use rad_step_setup_mod,      only: temp, rh2o, press, pflux, tflux,  &
				   ISRAD, IERAD, JSRAD, JERAD,  &
				   KSRAD, KERAD
use rad_utilities_mod,       only: locate_in_table, table_axis_type, &
				   longwave_control_type, Lw_control, &
				   radiation_control_type, Rad_control
use utilities_mod,           only: open_file, file_exist,    &
                                   check_nml_error, error_mesg,  &
			           print_version_number, FATAL, NOTE, &
				   WARNING, get_my_pe, close_file
use radiation_diag_mod,      only: radiag_from_setup,   &
	                           radiag_setup_dealloc

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                 longwave setup module
!  this module calculates and stores often-used quantities for the
!  longwave radiation modules.
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.08'
    character(len=128)  :: version =  '$Id: longwave_setup.F90,v 1.2 2001/08/30 15:13:37 fms Exp $'
    character(len=128)  :: tag     =  '$Name: galway $'



!---------------------------------------------------------------------
!------    interfaces   ------

public         &
	  longwave_setup_init, &
	  longwave_setup_driver,  &
	  longwave_setup_dealloc

public longwave_parameter_type


type longwave_parameter_type
      integer             :: offset  
      integer             :: NBTRG
      integer             :: NBTRGE
      integer             :: NLWCLDB
end type longwave_parameter_type


private longwave_tables_fill


!---------------------------------------------------------------------
!------     namelist  -----


logical                 :: do_ckd2p1 = .true.
logical                 :: do_ch4n2olbltmpint   = .false.


namelist / longwave_setup_nml /  &
				    do_ckd2p1,   &
				    do_ch4n2olbltmpint


!----------------------------------------------------------------------
!----  public data -------



type (table_axis_type),        public   ::    &
     temp_1 = table_axis_type(1, 100.0, 370.0, 10.0), &
     mass_1 = table_axis_type(1, -16.0,   1.9,  0.1)

type (longwave_parameter_type), public  :: Lw_parameters

real,    dimension(:,:,:), allocatable, public  ::  pdflux, pdfinv, &
						    tpl1, tpl2
real,    dimension(:,:,:), allocatable, public  ::  dte1, dte2
integer, dimension(:,:,:), allocatable, public  ::  ixoe1, ixoe2


!----------------------------------------------------------------------
!----  private data -------





!---------------------------------------------------------------------
!---------------------------------------------------------------------


contains



subroutine longwave_setup_init

!---------------------------------------------------------------------
   integer                    :: unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ---------

   if (file_exist('input.nml')) then
      unit =  open_file ('input.nml', action='read')
      ierr=1; do while (ierr /= 0)
      read (unit, nml=longwave_setup_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'longwave_setup_nml')
      enddo
10    call close_file (unit)
   endif

   unit = open_file ('logfile.out', action='append')
!  call print_version_number (unit, 'longwave_setup', version_number)
   if (get_my_pe() == 0)  then
     write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
     write (unit,nml=longwave_setup_nml)
   endif
   call close_file (unit)


   Lw_control%do_ckd2p1          = do_ckd2p1
   Lw_control%do_ch4n2olbltmpint = do_ch4n2olbltmpint    

   if (do_ckd2p1) then
     Lw_parameters%offset = 32
   else
     Lw_parameters%offset = 0
   endif

end subroutine longwave_setup_init



!#####################################################################

subroutine longwave_setup_driver                                       

!---------------------------------------------------------------------
     integer       :: k, j

!---------------------------------------------------------------------
     allocate ( pdflux (ISRAD:IERAD, JSRAD:JERAD,   KSRAD:KERAD   ) )
     allocate ( pdfinv (ISRAD:IERAD, JSRAD:JERAD,   KSRAD:KERAD   ) )
     allocate ( tpl1   (ISRAD:IERAD, JSRAD:JERAD,   KSRAD:KERAD+1 ) )
     allocate ( tpl2   (ISRAD:IERAD, JSRAD:JERAD,   KSRAD:KERAD+1 ) )

!-------------------------------------------------------------------- 
!     compute difference between the flux level pressures (pdflux)
!     which amounts to the layer size.
!-------------------------------------------------------------------- 
     do k=KSRAD,KERAD
       pdflux(:,:,k) = pflux(:,:,k+1) - pflux(:,:,k)
       pdfinv(:,:,k) = 1.0E+00/pdflux(:,:,k)
     enddo

!-------------------------------------------------------------------- 
!     compute mean temperature in the "nearby layer" between a flux
!     level and the first data level below the flux level (tpl1) or the
!     first data level above the flux level (tpl2)
!---------------------------------------------------------------------
     tpl1(:,:,KSRAD)         = temp(:,:,KERAD)
     tpl1(:,:,KSRAD+1:KERAD) = tflux(:,:,KSRAD+1:KERAD) 
     tpl1(:,:,KERAD+1)       = 0.5E+00*(tflux(:,:,KERAD+1) +   &
	   			        temp(:,:,KERAD))
     tpl2(:,:,KSRAD+1:KERAD) = tflux(:,:,KSRAD+1:KERAD) 
     tpl2(:,:,KERAD+1)       = 0.5E+00*(tflux(:,:,KERAD) +    &
				        temp(:,:,KERAD))

!--------------------------------------------------------------------
!    send temps to routine which will use them to locate indices in 
!    longwave tables.
!--------------------------------------------------------------------
     call longwave_tables_fill (temp, tflux)

!--------------------------------------------------------------------
!send  necessary data to radiation_diag_mod
!--------------------------------------------------------------------
     if (Rad_control%do_diagnostics) then
       do j=JSRAD, JERAD
         call radiag_from_setup (j, pdflux, pdfinv, pflux, temp,   &
                                 press, rh2o) 
       end do
     endif
!--------------------------------------------------------------------


end subroutine longwave_setup_driver




!####################################################################

subroutine longwave_setup_dealloc


     if (Rad_control%do_diagnostics) then
       call radiag_setup_dealloc
     endif

     deallocate (ixoe2       )
     deallocate (dte2        )
     deallocate (ixoe1       )
     deallocate (dte1        )
     deallocate (tpl2        )
     deallocate (tpl1        )
     deallocate (pdfinv      )
     deallocate (pdflux      )
!--------------------------------------------------------------------


end subroutine longwave_setup_dealloc



!####################################################################

subroutine longwave_tables_fill (temp, tflux)

!--------------------------------------------------------------------
real, dimension (:,:,:), intent(in)  :: temp, tflux
!--------------------------------------------------------------------

  allocate ( dte1  (ISRAD:IERAD, JSRAD:JERAD,   KSRAD:KERAD+1) )
  allocate ( ixoe1 (ISRAD:IERAD, JSRAD:JERAD,   KSRAD:KERAD+1) )

  call locate_in_table(temp_1, temp, dte1, ixoe1, KSRAD, KERAD+1)

  allocate ( dte2  (ISRAD:IERAD, JSRAD:JERAD,   KSRAD:KERAD+1) )
  allocate ( ixoe2 (ISRAD:IERAD, JSRAD:JERAD,   KSRAD:KERAD+1) )

  call locate_in_table(temp_1, tflux, dte2, ixoe2, KSRAD, KERAD+1)

  ixoe2(:,:,KSRAD:KERAD) = ixoe2(:,:,KSRAD+1:KERAD+1)
  dte2 (:,:,KSRAD:KERAD) = dte2 (:,:,KSRAD+1:KERAD+1)
  ixoe2(:,:,KERAD+1)     = ixoe1(:,:,KERAD)
  dte2 (:,:,KERAD+1)     = dte1 (:,:,KERAD)



end subroutine longwave_tables_fill




!####################################################################



	      end module longwave_setup_mod


