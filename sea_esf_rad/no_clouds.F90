
                 module no_clouds_mod


use utilities_mod,         only:  open_file, file_exist,   &
                                  check_nml_error, error_mesg,   &
                                  print_version_number, FATAL, NOTE, &
				  WARNING, get_my_pe, close_file
!use rad_step_setup_mod,    only:  ISRAD, IERAD, JSRAD, JERAD, & 
!                                  KSRAD, KERAD

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!	   module to properly define cloud radiative properties 
!          when model is run without clouds using eta coordinates
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.08'
   character(len=128)  :: version =  '$Id: no_clouds.F90,v 1.3 2002/07/16 22:36:06 fms Exp $'
   character(len=128)  :: tag     =  '$Name: inchon $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
	  no_clouds_init, no_clouds_calc


private         &
	  step_mtn_clouds


!---------------------------------------------------------------------
!-------- namelist  ---------

real                         :: dummy = 0.0




namelist /no_clouds_nml /     &
			       dummy


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------


       integer :: israd, ierad, jsrad, jerad, ksrad, kerad

!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





subroutine no_clouds_init 

      integer      :: unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------
  
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=no_clouds_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'no_clouds_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'no_clouds', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=no_clouds_nml)
      endif
      call close_file (unit)


end subroutine no_clouds_init


!#################################################################

subroutine no_clouds_calc(ncldsw, nrndlw, camtsw, crndlw, &
			  emrndlw, kbot)

!--------------------------------------------------------------------
integer, dimension(:,:), intent(inout)        ::  ncldsw, nrndlw
real   , dimension(:,:,:), intent(inout)      ::  camtsw, crndlw
real   , dimension(:,:,:,:), intent(inout)    ::  emrndlw   
integer, dimension(:,:), intent(in), optional ::  kbot
!--------------------------------------------------------------------

       israd = 1
       ierad = size (camtsw, 1)
       jsrad = 1
       jerad = size (camtsw, 2)
       ksrad = 1
       kerad = size (camtsw, 3)

!--------------------------------------------------------------------
!!  IS THIS (OR SOMETHING LIKE IT) NEEDED IN THIS FORMULATION ?

     if (present (kbot)) call step_mtn_clouds (kbot, &
		               ncldsw, nrndlw, camtsw, crndlw, emrndlw)

end subroutine no_clouds_calc



!####################################################################

subroutine step_mtn_clouds (kbot, ncldsw, nrndlw, camtsw, &
			    crndlw, emrndlw)

!--------------------------------------------------------------------
!!! ???? IS THIS ROUTINE NEEDED AND DOING THE RIGHT THINGS ???

integer, intent(in),    dimension(:,:)        :: kbot
integer, dimension(:,:), intent(inout)        ::  ncldsw, nrndlw
real   , dimension(:,:,:), intent(inout)      ::  camtsw, crndlw
real   , dimension(:,:,:,:), intent(inout)    ::  emrndlw   

!--------------------------------------------------------------------
       integer  ::   i, j, k

       do j=JSRAD,JERAD
         do i=ISRAD,IERAD
	   do k=kbot(i,j)+1, KERAD
             ncldsw(i,j)=ncldsw(i,j)+1
	     nrndlw(i,j) = nrndlw(i,j) + 1
	     camtsw(i,j,k) = 1.0
	     crndlw(i,j,k) = 1.0
	     emrndlw(i,j,k,:) = 1.0
           enddo
         enddo
       enddo



end subroutine step_mtn_clouds 



!###################################################################




	       end module no_clouds_mod
