	     module rad_output_file_mod


use  utilities_mod,       only:  open_file, file_exist,    &
                                 check_nml_error, error_mesg, &
                                 print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file, &
				 get_num_pes, get_domain_decomp
use time_manager_mod,     only:  time_type
use diag_manager_mod,     only:  register_diag_field, send_data
use rad_utilities_mod,    only:  Environment, environment_type, &
				 Rad_control, radiation_control_type

!------------------------------------------------------------------

implicit none
private

!------------------------------------------------------------------
!                  collects data for, writes and reads 
!                         a radiation output file
!
!------------------------------------------------------------------



!-------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: rad_output_file.F90,v 1.2 2001/07/05 17:32:37 fms Exp $'
   character(len=128)  :: tag     =  '$Name: fez $'

!---------------------------------------------------------------------
!-------  interfaces --------

public   &
            rad_output_file_init, write_rad_output_file,   &
	    rad_output_file_dealloc, rad_output_file_end,  &
            hold_qo3, hold_sfcalbedo, hold_cloud, hold_sfctmp, &
            hold_lw, hold_sw, hold_gases, hold_clouds, hold_tendencies,&
	    put_cloudpackage_type
     

private   &
	    read_history10, register_fields


!-------------------------------------------------------------------
!-------- namelist  ---------

logical :: write_hist10=.false.
logical :: test_read=.false.
integer :: i_prt=0, j_prt=0, jd_prt=0
logical :: write_netcdf_file=.false.


namelist  / rad_output_file_nml /  &
				    test_read, write_hist10, &
				    i_prt, j_prt, jd_prt, &
				    write_netcdf_file
				      

!------------------------------------------------------------------
!----public data------
!  (these are used in SKYHI by archive_package.F)

real, dimension(:,:,:), public, allocatable ::    &
				       qo3, cmxolw, crndlw, flxnet, &
				       heatra, flxnetcf, heatracf, &
				       fsw, ufsw, fswcf, ufswcf, &
				       temp, rh2o, radp, radswp, &
				       radpcf, radswpcf

real, dimension(:,:), public, allocatable ::    &
				       cirrfgd, cvisrfgd, tmpsfc, psj,&
				       tot_clds, cld_isccp_hi, &
				       cld_isccp_mid, cld_isccp_low

real, public                              :: rrvco2, rrvch4, rrvn2o,  &
					     rrvf11, rrvf12, rrvf113, &
					     rrvf22



!------------------------------------------------------------------
!----private data------

integer                     :: file10
integer                     :: idim, jdim, kdim
integer                     :: entry_counter=0
logical                     :: do_cloudforcing=.false.
logical                     :: do_isccp=.false.
integer                     :: x(4), y(4)

!--------------------------------------------------------------------
! netcdf diagnostics field variables
!--------------------------------------------------------------------
character(len=9), parameter :: mod_name='radiation'
real                        :: missing_value = -999.
integer                     :: id_radswp, id_radp, id_temp, id_rh2o, &
			       id_qo3, id_cmxolw, id_crndlw, id_flxnet,&
			       id_fsw, id_ufsw, id_psj, id_tmpsfc, &
			       id_cvisrfgd, id_cirrfgd, id_rrvco2, &
			       id_rrvch4, id_rrvn2o, id_rrvf11,  &
			       id_rrvf12, id_rrvf113, id_rrvf22, &
			       id_tot_clds, id_cld_isccp_hi,  &
			       id_cld_isccp_mid, id_cld_isccp_low, &
			       id_radswpcf, id_radpcf, id_flxnetcf, &
			       id_fswcf, id_ufswcf



!---------------------------------------------------------------------
!---------------------------------------------------------------------


contains



subroutine rad_output_file_init (axes, Time)

integer, dimension(4), intent(in), optional    :: axes
type(time_type),       intent(in), optional    :: Time

integer   :: unit, io, ierr

!--------------------------------------------------------------------
!-----  read namelist  ------

    if (file_exist('input.nml')) then
      unit =  open_file ('input.nml', action='read')
      ierr=1; do while (ierr /= 0)
      read (unit, nml=rad_output_file_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'rad_output_file_nml')
      enddo
10    call close_file (unit)
    endif

    unit = open_file ('logfile.out', action='append')
!   call print_version_number (unit, 'rad_output_file', version_number)
    if (get_my_pe() == 0) then
      write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
      write (unit,nml=rad_output_file_nml)
    endif
    call close_file (unit)

    if (write_hist10 .and. Environment%running_standalone) then
      call error_mesg ('rad_output_file_init', &
	 'cannot write history file 10 when running standalone', FATAL)
    endif

    if (write_hist10 .and. write_netcdf_file) then
      call error_mesg ('rad_output_file_init', &
	 'writing both history10 AND netcdf file in rad_output_mod', &
								 NOTE)
    endif

    if (write_hist10) then
      if (get_num_pes() /= 1) then
	call error_mesg ('rad_output_file_init', &
          'history file 10 not currently writable when npes > 1', FATAL)
      else
        file10 = open_file ('history10', action='write', form='ieee32')
      endif
    endif

    if (test_read .and. get_num_pes() /= 1) then
      call error_mesg ( 'rad_output_file_init',  &
                " can't test-read a file written by more than 1 pe ",&  
								FATAL)
    endif

    if (test_read .and. .not. write_hist10) then
      call error_mesg ( 'rad_output_file_init',  &
                " can't test-read a file which is not being written",&  
								FATAL)
    endif

    call get_domain_decomp (x, y)

    if (write_hist10 .and. test_read) then
      if (i_prt <= 0 .or. j_prt <= 0 .or. jd_prt <= 0   .or. &
	  i_prt > x(2) .or. jd_prt > y(2)) then
        call error_mesg ( 'rad_output_file_init',  &
                ' desired column for test reading is invalid', FATAL)
      endif
    endif

    if (write_netcdf_file) then
      if (Environment%running_gcm .and. Environment%running_fms) then
        call register_fields (Time, axes)
      else
	call error_mesg ('rad_output_file_init', &
      ' can only write netcdf file in radiation when running in fms', &
							       FATAL)
      endif
    endif

end subroutine rad_output_file_init



!################################################################

subroutine write_rad_output_file (safe_to_dealloc, Time, is, js)

!----------------------------------------------------------------
type(time_type), intent(in), optional :: Time
integer,         intent(in), optional :: is, js
logical,         intent(in)           :: safe_to_dealloc
!------------------------------------------------------------------

    logical             :: used


    if (write_hist10) then
      write (file10) radswp, radp, temp, rh2o, qo3, cmxolw, crndlw   
      write (file10) flxnet, fsw, ufsw   
      write (file10) psj, tmpsfc, cvisrfgd, cirrfgd
      write (file10) rrvco2, rrvch4, rrvn2o, rrvf11, rrvf12, rrvf113, &
 	             rrvf22
      if (do_isccp) then
        write (file10) tot_clds, cld_isccp_hi, cld_isccp_mid, &
   	               cld_isccp_low
      endif
      if (do_cloudforcing) then
        write (file10) radswpcf, radpcf, flxnetcf, fswcf, ufswcf
      endif

    endif

  if (write_netcdf_file) then
      
      if (id_radswp > 0 ) then
	used = send_data (id_radswp, radswp, Time, is, js, 1)
      endif

      if (id_radp > 0 ) then
	used = send_data (id_radp, radp, Time, is, js, 1)
      endif

      if (id_temp > 0 ) then
	used = send_data (id_temp, temp, Time, is, js, 1)
      endif

      if (id_rh2o > 0 ) then
	used = send_data (id_rh2o, rh2o, Time, is, js, 1)
      endif

      if (id_qo3  > 0 ) then
	used = send_data (id_qo3, qo3, Time, is, js, 1)
      endif

      if (id_cmxolw > 0 ) then
	used = send_data (id_cmxolw, cmxolw, Time, is, js, 1)
      endif

      if (id_crndlw > 0 ) then
	used = send_data (id_crndlw, crndlw, Time, is, js, 1)
      endif

      if (id_flxnet > 0 ) then
	used = send_data (id_flxnet, flxnet, Time, is, js, 1)
      endif

      if (id_fsw    > 0 ) then
	used = send_data (id_fsw   , fsw   , Time, is, js, 1)
      endif

      if (id_ufsw   > 0 ) then
	used = send_data (id_ufsw  , ufsw  , Time, is, js, 1)
      endif

      if (id_psj    > 0 ) then
	used = send_data (id_psj   , psj   , Time, is, js)
      endif

      if (id_tmpsfc > 0 ) then
	used = send_data (id_tmpsfc, tmpsfc, Time, is, js)
      endif

      if (id_cvisrfgd > 0 ) then
	used = send_data (id_cvisrfgd, cvisrfgd, Time, is, js)
      endif

      if (id_cirrfgd > 0 ) then
	used = send_data (id_cirrfgd , cirrfgd , Time, is, js)
      endif

    if (do_isccp) then
      if (id_tot_clds > 0 ) then
	used = send_data (id_tot_clds , tot_clds, Time, is, js)
      endif

      if (id_cld_isccp_hi > 0 ) then
	used = send_data (id_cld_isccp_hi, cld_isccp_hi, Time, is, js)
      endif

      if (id_cld_isccp_mid > 0 ) then
	used = send_data (id_cld_isccp_mid, cld_isccp_mid, Time, is, js)
      endif

      if (id_cld_isccp_low > 0 ) then
	used = send_data (id_cld_isccp_low, cld_isccp_low, Time, is, js)
      endif
    endif

    if (Rad_control%do_totcld_forcing) then
      if (id_radswpcf > 0 ) then
	used = send_data (id_radswpcf , radswpcf, Time, is, js, 1)
      endif

      if (id_radpcf > 0 ) then
	used = send_data (id_radpcf, radpcf, Time, is, js, 1)
      endif

      if (id_flxnetcf > 0 ) then
	used = send_data (id_flxnetcf, flxnetcf, Time, is, js, 1)
      endif

      if (id_fswcf  > 0 ) then
	used = send_data (id_fswcf , fswcf , Time, is, js, 1)
      endif

      if (id_ufswcf  > 0 ) then
	used = send_data (id_ufswcf , ufswcf , Time, is, js, 1)
      endif
    endif
 endif

!--------------------------------------------------------------------
!   in skyhi do not deallocate since file 10 will be written by the
!   archive_package and will need these variables. in fms, the only
!   file using these variables will have just been written, and so 
!   they may be deallocated.
!-------------------------------------------------------------------
    if (safe_to_dealloc) then
      call rad_output_file_dealloc
    endif


!------------------------------------------------------------------


end subroutine write_rad_output_file




!############################################################


subroutine rad_output_file_dealloc

!-------------------------------------------------------------------
    if (write_hist10 .or. write_netcdf_file) then
      deallocate ( qo3      )
      if (allocated (cmxolw) ) then
        deallocate ( cmxolw   )
        deallocate ( crndlw   )
      endif
      deallocate ( flxnet   )
      deallocate ( heatra   )
      if (allocated (fsw) ) then
        deallocate ( fsw      )
        deallocate ( ufsw     )
        deallocate ( cvisrfgd )
        deallocate ( cirrfgd  )
      endif
      deallocate ( tmpsfc   )
      deallocate ( temp     )
      deallocate ( rh2o     )
      deallocate ( psj      )
      if (allocated (radpcf) ) then
        deallocate ( flxnetcf  )
        deallocate ( heatracf  )
	if (allocated (fswcf) ) then
          deallocate ( fswcf     )
          deallocate ( ufswcf    )
	endif
        deallocate ( radpcf    )
        deallocate ( radswpcf    )
      endif
      deallocate ( radp    )
      deallocate ( radswp    )
      if (allocated (tot_clds) ) then
        deallocate ( tot_clds    )
        deallocate ( cld_isccp_hi )
        deallocate ( cld_isccp_mid)
        deallocate ( cld_isccp_low)
      endif
    endif
!-----------------------------------------------------------------



end subroutine rad_output_file_dealloc


!#####################################################################

subroutine rad_output_file_end


   if (write_hist10) then
     call close_file (file10)
     if (test_read) then
       file10 = open_file ('history10', action='read', form='ieee32')
       call read_history10
       call close_file (file10)
     endif
   endif

end subroutine rad_output_file_end



!###################################################################

subroutine hold_qo3 (qo3_in)

real, dimension(:,:,:), intent(in) :: qo3_in

 if (write_hist10 .or. write_netcdf_file) then
   allocate ( qo3( size(qo3_in,1), size(qo3_in,2), size(qo3_in,3)) )
   qo3(:,:,:) = qo3_in(:,:,:)
 endif

end subroutine hold_qo3



!#####################################################################

subroutine hold_sfcalbedo (cirrfgd_in, cvisrfgd_in)

real, dimension(:,:), intent(in)  :: cirrfgd_in, cvisrfgd_in


 if (write_hist10 .or. write_netcdf_file) then
   allocate ( cirrfgd( size(cirrfgd_in,1), size(cirrfgd_in,2)) )
   allocate ( cvisrfgd( size(cvisrfgd_in,1), size(cvisrfgd_in,2)) )

   cirrfgd (:,:) = cirrfgd_in (:,:)
   cvisrfgd(:,:) = cvisrfgd_in(:,:)
 endif

end subroutine hold_sfcalbedo



!##################################################################

subroutine hold_cloud (cmxolw_in, crndlw_in)


real, dimension(:,:,:), intent(in) :: cmxolw_in, crndlw_in


 if (write_hist10 .or. write_netcdf_file) then
   allocate ( cmxolw( size(cmxolw_in,1), size(cmxolw_in,2),    &
                      size(cmxolw_in,3)) )
   allocate ( crndlw( size(crndlw_in,1), size(crndlw_in,2),    &
                      size(crndlw_in,3)) )

   cmxolw(:,:,:) = cmxolw_in(:,:,:)
   crndlw(:,:,:) = crndlw_in(:,:,:)
 endif

end subroutine hold_cloud



!##################################################################

subroutine hold_sfctmp (temp_in, rh2o_in, press_in)

real,intent(in),  dimension(:,:,:)  ::  temp_in, rh2o_in, press_in 


 integer :: KERAD

 if (write_hist10 .or. write_netcdf_file) then
   allocate ( tmpsfc ( size(temp_in,1), size(temp_in,2)) )
   allocate ( psj    ( size(press_in,1), size(press_in,2)) )
   allocate ( temp   ( size(temp_in,1), size(temp_in,2), &
                       size(temp_in,3)-1) )
   allocate ( rh2o   ( size(rh2o_in,1), size(rh2o_in,2), &
                       size(rh2o_in,3)  ) )


   entry_counter = entry_counter + 1
   KERAD = ubound(temp_in,3) - 1
   idim = size(temp_in,1)
   jdim = size(temp_in,2)
   kdim = size(temp_in,3) -1 

   tmpsfc(:,:) = temp_in(:,:,KERAD+1)
   psj   (:,:) = press_in(:,:,KERAD+1)
   temp(:,:,:) = temp_in(:,:,1:KERAD)
   rh2o(:,:,:) = rh2o_in(:,:,:)
 endif

end subroutine hold_sfctmp



!##################################################################

subroutine hold_lw (flxnet_in, heatra_in, flxnetcf_in, &
                    heatracf_in)

real, dimension(:,:,:), intent(in)           :: flxnet_in, heatra_in
real, dimension(:,:,:), intent(in), optional :: flxnetcf_in, heatracf_in



 if (write_hist10 .or. write_netcdf_file) then
   allocate ( flxnet( size(flxnet_in,1), size(flxnet_in,2),    &
                      size(flxnet_in,3)) )
   allocate ( heatra( size(heatra_in,1), size(heatra_in,2),    &
                      size(heatra_in,3)) )
   flxnet(:,:,:) = flxnet_in(:,:,:)
   heatra(:,:,:) = heatra_in(:,:,:)

   if (present(flxnetcf_in) ) then
     allocate ( flxnetcf( size(flxnetcf_in,1), size(flxnetcf_in,2),   & 
                          size(flxnetcf_in,3)) )
     allocate ( heatracf( size(heatracf_in,1), size(heatracf_in,2),   &
                          size(heatracf_in,3)) )
     flxnetcf(:,:,:) = flxnetcf_in(:,:,:)
     heatracf(:,:,:) = heatracf_in(:,:,:)
   endif
 endif

end subroutine hold_lw



!##################################################################

subroutine hold_sw   (fsw_in, ufsw_in, fswcf_in, ufswcf_in)

real, dimension(:,:,:), intent(in)           :: fsw_in, ufsw_in
real, dimension(:,:,:), intent(in), optional :: fswcf_in, ufswcf_in


   if (write_hist10 .or. write_netcdf_file) then
     allocate ( fsw   ( size(fsw_in,1), size(fsw_in,2),    &
                        size(fsw_in,3)) )
     allocate ( ufsw  ( size(ufsw_in,1), size(ufsw_in,2),    &
                        size(ufsw_in,3)) )
     fsw(:,:,:) = fsw_in(:,:,:)
     ufsw(:,:,:) = ufsw_in(:,:,:)

     if (present(fswcf_in) ) then
       allocate ( fswcf( size(fswcf_in,1), size(fswcf_in,2),    &
                         size(fswcf_in,3)) )
       allocate ( ufswcf( size(ufswcf_in,1), size(ufswcf_in,2),    &
                          size(ufswcf_in,3)) )
       fswcf(:,:,:) = fswcf_in(:,:,:)
       ufswcf(:,:,:) = ufswcf_in(:,:,:)
     endif
   endif

end subroutine hold_sw
       


!##################################################################

subroutine hold_gases (rrvco2_in, rrvch4_in, rrvn2o_in, rrvf11_in, &
		       rrvf12_in, rrvf113_in, rrvf22_in)

real, intent(in)    :: rrvco2_in, rrvch4_in, rrvn2o_in, rrvf11_in, &
		       rrvf12_in, rrvf113_in, rrvf22_in

   rrvco2 = rrvco2_in
   rrvch4 = rrvch4_in
   rrvn2o = rrvn2o_in
   rrvf11 = rrvf11_in
   rrvf12 = rrvf12_in
   rrvf113 = rrvf113_in
   rrvf22 = rrvf22_in


end subroutine hold_gases
       


!##################################################################

subroutine hold_clouds (tot_clds_in, cld_isccp_hi_in, cld_isccp_mid_in,&
			cld_isccp_low_in)

real, dimension(:,:), intent(in)  :: tot_clds_in, cld_isccp_hi_in, &
				     cld_isccp_mid_in, cld_isccp_low_in

   if (write_hist10 .or. write_netcdf_file) then
     allocate ( tot_clds  ( size(tot_clds_in,1), size(tot_clds_in,2)))
     allocate ( cld_isccp_hi (size(cld_isccp_hi_in,1),  &
                              size(cld_isccp_hi_in,2)))
     allocate ( cld_isccp_mid(size(cld_isccp_mid_in,1),  &
                              size(cld_isccp_mid_in,2)))
     allocate ( cld_isccp_low(size(cld_isccp_low_in,1),  &
                              size(cld_isccp_low_in,2)))

     tot_clds(:,:) = tot_clds_in(:,:)
     cld_isccp_hi (:,:) = cld_isccp_hi_in (:,:)
     cld_isccp_mid(:,:) = cld_isccp_mid_in(:,:)
     cld_isccp_low(:,:) = cld_isccp_low_in(:,:)
   endif

end subroutine hold_clouds
       

!##################################################################

subroutine hold_tendencies (radp_in, radswp_in, radpcf_in, radswpcf_in)

real, dimension(:,:,:), intent(in)           :: radp_in, radswp_in
real, dimension(:,:,:), intent(in), optional :: radpcf_in, radswpcf_in

   if (write_hist10 .or. write_netcdf_file) then
     allocate ( radp  ( size(radp_in,1), size(radp_in,2),    &
                        size(radp_in,3)) )
     allocate ( radswp( size(radswp_in,1), size(radswp_in,2),    &
                        size(radswp_in,3)) )
     radp(:,:,:) = radp_in(:,:,:)
     radswp(:,:,:) = radswp_in(:,:,:)

     if (present(radpcf_in) ) then
       allocate (radpcf(size(radpcf_in,1), size(radpcf_in,2),    &
                        size(radpcf_in,3)) )
       allocate (radswpcf(size(radswpcf_in,1), size(radswpcf_in,2), &
                          size(radswpcf_in,3)) )

       radpcf(:,:,:) = radpcf_in(:,:,:)
       radswpcf(:,:,:) = radswpcf_in(:,:,:)
       do_cloudforcing = .true.
     endif
   endif

end subroutine hold_tendencies
       

!###################################################################


subroutine put_cloudpackage_type (form)

character(len=*), intent(in)    :: form

     if (trim(form) == 'rh') then
       do_isccp = .true.
     else
       do_isccp = .false.
     endif
     
end subroutine put_cloudpackage_type


!####################################################################

subroutine read_history10         

!--------------------------------------------------------------------
   integer ::  jj, kt, k
   integer ::  time_levels_in_file, nrecs, jind, jrem, rows_per_rec
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   allocate arrays to hold file variables
!---------------------------------------------------------------------
   allocate ( radswp      (idim, jdim, kdim  ) )
   allocate ( radp        (idim, jdim, kdim  ) )
   allocate ( temp        (idim, jdim, kdim  ) )
   allocate ( rh2o        (idim, jdim, kdim  ) )
   allocate ( qo3         (idim, jdim, kdim  ) )
   allocate ( cmxolw      (idim, jdim, kdim  ) )
   allocate ( crndlw      (idim, jdim, kdim  ) )
   allocate ( flxnet      (idim, jdim, kdim+1) )
   allocate ( fsw         (idim, jdim, kdim+1) )
   allocate ( ufsw        (idim, jdim, kdim+1) )
   allocate ( psj         (idim, jdim        ) )
   allocate ( tmpsfc      (idim, jdim        ) )
   allocate ( cvisrfgd    (idim, jdim        ) )
   allocate ( cirrfgd     (idim, jdim        ) )
   if (do_isccp) then
     allocate ( tot_clds    (idim, jdim        ) )
     allocate ( cld_isccp_hi(idim, jdim        ) )
     allocate ( cld_isccp_mid(idim, jdim        ) )
     allocate ( cld_isccp_low(idim, jdim        ) )
   endif
   if (do_cloudforcing) then
     allocate ( radswpcf    (idim, jdim, kdim  ) )
     allocate ( radpcf      (idim, jdim, kdim  ) )
     allocate ( flxnetcf    (idim, jdim, kdim+1) )
     allocate ( fswcf       (idim, jdim, kdim+1) )
     allocate ( ufswcf      (idim, jdim, kdim+1) )
   end if

!--------------------------------------------------------------------
!  rewind file, determine how many time levels are in it
!--------------------------------------------------------------------

   time_levels_in_file = entry_counter/y(2)

!--------------------------------------------------------------------
!  determine number of records in file and which one will have 
!  specified grid column to be printed
!  read file, print out data at specified grid column
!--------------------------------------------------------------------
   nrecs = y(2)/jdim
   rows_per_rec = jdim
   jind = jd_prt/rows_per_rec
   jrem = jd_prt -jind*rows_per_rec
   if (jrem > 0) jind = jind + 1
   j_prt = jd_prt - (jind-1)*rows_per_rec 
   do kt = 1,time_levels_in_file
     write (*, FMT = "( 'time level # ',  i4 )" )  kt
     do jj=1,nrecs
       read  (file10) radswp, radp, temp, rh2o,  qo3, cmxolw, crndlw 
       read  (file10) flxnet, fsw, ufsw   
       read  (file10) psj, tmpsfc, cvisrfgd, cirrfgd
       read  (file10) rrvco2, rrvch4, rrvn2o, rrvf11, rrvf12,  &
                      rrvf113, rrvf22
       if (do_isccp) then
         read  (file10) tot_clds, cld_isccp_hi, cld_isccp_mid, &
                        cld_isccp_low
       endif
       if (do_cloudforcing) then
         read  (file10) radswpcf, radpcf, flxnetcf, fswcf, ufswcf
       endif

       if (jj == jind) then
         write (*,FMT = '(a10)' ) 'qo3' 
         write (*,FMT = '(5e18.10)' ) (qo3     (i_prt,j_prt,k), &
	   				             k=1,kdim)
         write (*,FMT = '(a10)' ) 'cmxolw'
         write (*,FMT = '(5e18.10)' ) (cmxolw  (i_prt,j_prt,k), &
						       k=1,kdim)
         write (*,FMT = '(a10)' ) 'crndlw'
         write (*,FMT = '(5e18.10)' ) (crndlw  (i_prt,j_prt,k), &
						       k=1,kdim)
         write (*,FMT = '(a10)' ) 'flxnet'
         write (*,FMT = '(5e18.10)' ) (flxnet  (i_prt,j_prt,k), &
						       k=1,kdim+1)
         write (*,FMT = '(a10)' ) 'fsw'
         write (*,FMT = '(5e18.10)' ) (fsw     (i_prt,j_prt,k), &
						       k=1,kdim+1)
         write (*,FMT = '(a10)' ) 'ufsw'
         write (*,FMT = '(5e18.10)' ) (ufsw    (i_prt,j_prt,k), &
						       k=1,kdim+1)
         write (*,FMT = '(a10)' ) 'temp'
         write (*,FMT = '(5e18.10)' ) (temp    (i_prt,j_prt,k), &
						       k=1,kdim)
         write (*,FMT = '(a10)' ) 'rh2o'
         write (*,FMT = '(5e18.10)' ) (rh2o    (i_prt,j_prt,k), &
					       k=1,kdim)
         write (*,FMT = '(a10)' ) 'tmpsfc'
         write (*,FMT = '(5e18.10)' )  tmpsfc  (i_prt,j_prt)
         write (*,FMT = '(a10)' ) 'psj'
         write (*,FMT = '(5e18.10)' )  psj     (i_prt,j_prt)
         write (*,FMT = '(a10)' ) 'cirrfgd'
         write (*,FMT = '(5e18.10)' )  cirrfgd (i_prt,j_prt)
         write (*,FMT = '(a10)' ) 'cvisrfgd'
         write (*,FMT = '(5e18.10)' ) cvisrfgd (i_prt,j_prt)
         if (allocated(radpcf)) then
	   write (*,FMT = '(a10)' ) 'radswpcf'
	   write (*,FMT = '(5e18.10)' ) (radswpcf(i_prt,j_prt,k), &
						       k=1,kdim)
	   write (*,FMT = '(a10)' ) 'radpcf'
	   write (*,FMT = '(5e18.10)' ) (radpcf  (i_prt,j_prt,k), &
						       k=1,kdim)
	   write (*,FMT = '(a10)' ) 'flxnetcf'
	   write (*,FMT = '(5e18.10)' ) (flxnetcf(i_prt,j_prt,k), &
						       k=1,kdim+1)
	   write (*,FMT = '(a10)' ) 'fswcf'
	   write (*,FMT = '(5e18.10)' ) (fswcf   (i_prt,j_prt,k), &
						       k=1,kdim+1)
	   write (*,FMT = '(a10)' ) 'ufswcf'
	   write (*,FMT = '(5e18.10)' ) (ufswcf  (i_prt,j_prt,k), &
						       k=1,kdim+1)
         endif
         if (allocated (tot_clds)) then
	   write (*,FMT = '(a10)' ) 'tot_clds'
	   write (*,FMT = '(5e18.10)' ) tot_clds (i_prt,j_prt)
	   write (*,FMT = '(a10)' ) 'isccp_hi'
	   write (*,FMT = '(5e18.10)' ) cld_isccp_hi(i_prt,j_prt)
	   write (*,FMT = '(a10)' ) 'isccp_mid'
	   write (*,FMT = '(5e18.10)' ) cld_isccp_mid(i_prt,j_prt)
	   write (*,FMT = '(a10)' ) 'isccp_low'
	   write (*,FMT = '(5e18.10)' ) cld_isccp_low(i_prt,j_prt)
         endif
         write (*,FMT = '(a16)' ) 'radiative_gases'
         write (*,FMT = '(5e18.10)' ) rrvco2, rrvch4, rrvn2o, &
	      		                   rrvf11, rrvf12, rrvf113, &
					   rrvf22
       endif
     end do
   end do

!--------------------------------------------------------------------


end subroutine read_history10 



!#################################################################

subroutine register_fields (Time, axes)

integer, dimension(4), intent(in)          :: axes
type(time_type),       intent(in)          :: Time
!---------------------------------------------------------------------

     integer, dimension(4)    :: bxes
!---------------------------------------------------------------------
 
     bxes(1:2) = axes(1:2)
     bxes(3) = axes(4)
     bxes(4) = axes(4)

     id_radswp = &
	 register_diag_field (mod_name, 'radswp', axes(1:3), Time, &
			      'temperature tendency for SW radiation', &
			      'deg_K/sec', missing_value=missing_value)

     id_radp = &
	 register_diag_field (mod_name, 'radp', axes(1:3), Time, &
			      'temperature tendency for radiation', &
			      'deg_K/sec', missing_value=missing_value)

     id_temp   = &
	 register_diag_field (mod_name, 'temp', axes(1:3), Time, &
			      'temperature field', &
			      'deg_K', missing_value=missing_value)

     id_rh2o   = &
	 register_diag_field (mod_name, 'rh2o', axes(1:3), Time, &
			      'water vapor mixing ratio', &
			      'kg/kg', missing_value=missing_value)

     id_qo3    = &
	 register_diag_field (mod_name, 'qo3', axes(1:3), Time, &
			      'ozone mixing ratio', &
			      'kg/kg', missing_value=missing_value)

     id_cmxolw = &
	 register_diag_field (mod_name, 'cmxolw', axes(1:3), Time, &
			      'maximum overlap cloud amount', &
			      'unknown', missing_value=missing_value)

     id_crndlw = &
	 register_diag_field (mod_name, 'crndlw', axes(1:3), Time, &
			      'random overlap cloud amount', &
			      'unknown', missing_value=missing_value)

     id_flxnet = &
	 register_diag_field (mod_name, 'flxnet', bxes(1:3), Time, &
			      'net longwave radiative flux', &
			      'unknown', missing_value=missing_value)

     id_fsw    = &
	 register_diag_field (mod_name, 'fsw', bxes(1:3), Time, &
			      'net shortwave radiative flux', &
			      'unknown', missing_value=missing_value)

     id_ufsw   = &
	 register_diag_field (mod_name, 'ufsw', bxes(1:3), Time, &
			      'upward shortwave radiative flux ', &
			      'unknown', missing_value=missing_value)

     id_psj    = &
	 register_diag_field (mod_name, 'psj', axes(1:2), Time, &
			      'surface pressure', &
			      'pressure', missing_value=missing_value)

     id_tmpsfc = &
	 register_diag_field (mod_name, 'tmpsfc', axes(1:2), Time, &
			      'surface temperature', &
			      'deg_K', missing_value=missing_value)

     id_cvisrfgd = &
	 register_diag_field (mod_name, 'cvisrfgd', axes(1:2), Time, &
			      'visible surface albedo', &
			      'unknown', missing_value=missing_value)

     id_cirrfgd = &
	 register_diag_field (mod_name, 'cirrfgd', axes(1:2), Time, &
			      'visible infra-red albedo', &
			      'unknown', missing_value=missing_value)

   if (do_isccp) then
     id_tot_clds = &
	 register_diag_field (mod_name, 'tot_clds', axes(1:2), Time, &
			      'total isccp cloud cover', &
			      'unknown', missing_value=missing_value)

     id_cld_isccp_hi = &
       register_diag_field (mod_name, 'cld_isccp_hi', axes(1:2), Time, &
			      'isccp high clouds', &
			      'unknown', missing_value=missing_value)

     id_cld_isccp_mid = &
      register_diag_field (mod_name, 'cld_isccp_mid', axes(1:2), Time, &
			      'isccp middle clouds', &
			      'unknown', missing_value=missing_value)

     id_cld_isccp_low = &
      register_diag_field (mod_name, 'cld_isccp_low', axes(1:2), Time, &
			      'isccp low clouds', &
			      'unknown', missing_value=missing_value)

   endif 

   if (Rad_control%do_totcld_forcing) then
     id_radswpcf = &
	 register_diag_field (mod_name, 'radswpcf', axes(1:3), Time, &
			    'temperature forcing from sw w/o clouds', &
			      'deg_K/sec', missing_value=missing_value)

     id_radpcf = &
	 register_diag_field (mod_name, 'radpcf', axes(1:3), Time, &
			      'temperature forcing w/o clouds', &
			      'deg_K/sec', missing_value=missing_value)

     id_flxnetcf = &
	 register_diag_field (mod_name, 'flxnetcf', bxes(1:3), Time, &
			      'net longwave flux w/o clouds', &
			      'unknown', missing_value=missing_value)

     id_fswcf = &
	 register_diag_field (mod_name, 'fswcf', bxes(1:3), Time, &
			      'net shortwave flux w/o clouds', &
			      'unknown', missing_value=missing_value)

     id_ufswcf   = &
	 register_diag_field (mod_name, 'ufswcf', bxes(1:3), Time, &
			      'upward shortwave flux w/o clouds', &
			      'unknown', missing_value=missing_value)

   endif

end subroutine register_fields


!####################################################################





	     end module rad_output_file_mod
