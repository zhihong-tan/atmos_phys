
                 module donner_deep_clouds_W_mod

use time_manager_mod,       only: time_type
use donner_deep_mod,        only: donner_deep_avg
use utilities_mod,          only: open_file, file_exist,   &
                                  check_nml_error, error_mesg,   &
                                  print_version_number, FATAL, NOTE, &
				  WARNING, get_my_pe, close_file
!use rad_step_setup_mod,     only: temp, rh2o, press, pflux, jabs,   &
!				  iabs, ISRAD, IERAD, JSRAD, JERAD, & 
!                                  KSRAD, KERAD, land, &
!				  cloud_ice, cloud_water
use rad_utilities_mod,      only: Environment, environment_type, &
                                  longwave_control_type, Lw_control, &
				  shortwave_control_type, Sw_control,&
				  cld_diagnostics_type, &
                                  cloudrad_control_type, Cldrad_control
use microphys_rad_mod,      only: microphys_rad_driver

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!	          donner deep cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: donner_deep_clouds_W.F90,v 1.3 2003/04/09 20:59:10 fms Exp $'
   character(len=128)  :: tag     =  '$Name: inchon $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          donner_deep_clouds_init, donner_deep_clouds_calc

!---------------------------------------------------------------------
!-------- namelist  ---------

logical   :: using_dge_lw = .true.
logical   :: using_dge_sw = .true.



namelist /donner_deep_clouds_W_nml /     &
			       using_dge_sw, using_dge_lw


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!character(len=10)     :: swform
!character(len=16)     :: swform
logical               :: do_lwcldemiss
logical               :: do_esfsw         
logical               :: do_lhsw        


!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





subroutine donner_deep_clouds_init 

      integer            :: unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------
  
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=donner_deep_clouds_W_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'donner_deep_clouds_W_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!      call print_version_number (unit, 'donner_deep_clouds_W', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=donner_deep_clouds_W_nml)
      endif
      call close_file (unit)

!     swform = Sw_control%sw_form
      do_esfsw = Sw_control%do_esfsw
      do_lhsw = Sw_control%do_lhsw 
      do_lwcldemiss = Lw_control%do_lwcldemiss


end subroutine donner_deep_clouds_init



!#################################################################

subroutine donner_deep_clouds_calc (             &
                  is,ie,js,je,deltaz,press,temp, Cld_diagnostics, &
                                    cld_cell,               &
			 cldext_cell, cldsct_cell, cldasymm_cell,  &
				    abscoeff_cell,          &
                                    cld_meso,               &
			 cldext_meso, cldsct_meso, cldasymm_meso,  &
				    abscoeff_meso)


integer, intent(in) :: is,ie,js,je
real, dimension(:,:,:), intent(in) :: deltaz, press, temp
type(cld_diagnostics_type), intent(inout) :: Cld_diagnostics
real, dimension(:,:,:), intent(inout) :: cld_cell, cld_meso
real, dimension(:,:,:,:), intent(out) :: cldext_cell, cldsct_cell,  &
                                         cldasymm_cell, abscoeff_cell
real, dimension(:,:,:,:), intent(out) :: cldext_meso, cldsct_meso,  &
                                         cldasymm_meso, abscoeff_meso


!-------------------------------------------------------------------
!   local variables
!-------------------------------------------------------------------

integer   :: idim, jdim, kdim
real, dimension(size(cld_cell,1),size(cld_cell,2),size(cld_cell,3)) :: &
      cell_liquid_amt, cell_ice_amt, cell_liquid_size, cell_ice_size,  &
      meso_liquid_amt, meso_ice_amt, meso_liquid_size, meso_ice_size
integer :: unit

     idim = size(cld_cell,1)
     jdim = size(cld_cell,2)
     kdim = size(cld_cell,3)

if (Environment%running_gcm) then
!--------------------------------------------------------------------
! if (Environment%running_fms) then


!---------------------------------------------------------------------
!     obtain the deep cloud areas and properties
!---------------------------------------------------------------------


!     call donner_deep_avg(  iabs(ISRAD), jabs(JSRAD),               &
!     call donner_deep_avg(  is, js,               &
      call donner_deep_avg(  is, ie, js, je,           &
                             cell_cloud_frac_out  = cld_cell,        &
                             cell_liquid_amt_out  = cell_liquid_amt, &
                             cell_liquid_size_out = cell_liquid_size,&
                             cell_ice_amt_out     = cell_ice_amt,    &
                             cell_ice_size_out    = cell_ice_size,   &
                             meso_cloud_frac_out  = cld_meso,        &
                             meso_liquid_amt_out  = meso_liquid_amt, &
                             meso_liquid_size_out = meso_liquid_size,&
                             meso_ice_amt_out     = meso_ice_amt,    &
                             meso_ice_size_out    = meso_ice_size    )
	 
!      unit = open_file ('fort.149', action='append',threading='multi')
!      call print_version_number (unit, 'microphys_rad', version_number)
!        write (unit,*) ' donner_deep_clouds_w'
! write (unit,*) idim,jdim,kdim
! write (unit,*) ' jabs = ', jabs
!        write (unit,*) ' cld_cell'
! write (unit,*) cld_cell
!        write (unit,*) ' cell_liquid_amt'
! write (unit,*) cell_liquid_amt
!        write (unit,*) ' cell_liquid_size'
! write (unit,*) cell_liquid_size
!        write (unit,*) ' cell_ice_amt'
! write (unit,*) cell_ice_amt
!        write (unit,*) ' cell_ice_size'
! write (unit,*) cell_ice_size
!      call close_file (unit)
  
!
!--------------------------------------------------------------------
!  if microphysics is being used for either sw or lw calculation, call
!  microphys_rad_driver to obtain cloud radiative properties.
!--------------------------------------------------------------------

          call microphys_rad_driver(                         &
                  is,ie,js,je,deltaz,press,temp, Cld_diagnostics, &
	                   conc_drop_in=cell_liquid_amt,        &
			   conc_ice_in=cell_ice_amt,            &
			   size_drop_in=cell_liquid_size,       &
			   size_ice_in=cell_ice_size,           &
			   cldext=cldext_cell,               &
			   cldsct=cldsct_cell,               &
			   cldasymm=cldasymm_cell,           &
			   abscoeff=abscoeff_cell,           &
			   using_dge_sw=using_dge_sw,        &
			   using_dge_lw=using_dge_lw)

          call microphys_rad_driver(                         &
                  is,ie,js,je,deltaz,press,temp, Cld_diagnostics, &
	                   conc_drop_in=meso_liquid_amt,        &
			   conc_ice_in=meso_ice_amt,            &
			   size_drop_in=meso_liquid_size,       &
			   size_ice_in=meso_ice_size,           &
			   cldext=cldext_meso,               &
			   cldsct=cldsct_meso,               &
			   cldasymm=cldasymm_meso,           &
			   abscoeff=abscoeff_meso,           &
			   using_dge_sw=using_dge_sw,        &
			   using_dge_lw=using_dge_lw)


!     else  ! (gcm, not running fms)

!       call error_mesg ('donner_deep_clouds_W',   &
!         'donner_deep_cloud only available in FMS', FATAL)

!   endif ! (running_gcm)

  else if (Environment%running_standalone) then ! (running_standalone)
 
     if (Cldrad_control%do_pred_cld_microphys) then

!
!--------------------------------------------------------------------
!  if microphysics is being used for either sw or lw calculation, call
!  microphys_rad_driver to obtain cloud radiative properties.
!--------------------------------------------------------------------
          call microphys_rad_driver(                         &
                  is,ie,js,je,deltaz,press,temp, Cld_diagnostics, &
	                   conc_drop_in=cell_liquid_amt,        &
			   conc_ice_in=cell_ice_amt,            &
			   size_drop_in=cell_liquid_size,       &
			   size_ice_in=cell_ice_size,           &
			   cldext=cldext_cell,               &
			   cldsct=cldsct_cell,               &
			   cldasymm=cldasymm_cell,           &
			   abscoeff=abscoeff_cell,           &
			   using_dge_sw=using_dge_sw,        &
			   using_dge_lw=using_dge_lw)

          call microphys_rad_driver(                         &
                  is,ie,js,je,deltaz,press,temp, Cld_diagnostics, &
	                   conc_drop_in=meso_liquid_amt,        &
			   conc_ice_in=meso_ice_amt,            &
			   size_drop_in=meso_liquid_size,       &
			   size_ice_in=meso_ice_size,           &
			   cldext=cldext_meso,               &
			   cldsct=cldsct_meso,               &
			   cldasymm=cldasymm_meso,           &
			   abscoeff=abscoeff_meso,           &
			   using_dge_sw=using_dge_sw,        &
			   using_dge_lw=using_dge_lw)

  


  else !  (standalone, not with pred_cld_microphys)

       call error_mesg ('donner_deep_clouds_W',   &
    ' standalone donner_deep_cloud only available with pred microphys', FATAL)
 
  
  endif  ! (standalone, microphys)

 endif  ! (gcm or standalone)


end subroutine donner_deep_clouds_calc


!####################################################################


	       end module donner_deep_clouds_W_mod



