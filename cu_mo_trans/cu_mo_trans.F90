
module cu_mo_trans_mod

!=======================================================================
!
!                 CUNVECTIVE MOMENTUM TRANSPORT MODULE
!
!=======================================================================

  use   constants_mod, only:  GRAV, RDGAS, RVGAS
 

  use         fms_mod, only: file_exist, check_nml_error,    &
                             open_namelist_file, close_file, &
                             write_version_number,           &
                             mpp_pe, mpp_root_pe, stdlog,    &
                             error_mesg, FATAL, NOTE

  use  Diag_Manager_Mod, ONLY: register_diag_field, send_data
  use  Time_Manager_Mod, ONLY: time_type

implicit none
private


! public interfaces
!=======================================================================
public :: cu_mo_trans_init, &
          cu_mo_trans,      &
          cu_mo_trans_end

!=======================================================================

! form of interfaces
!=======================================================================

      
logical :: initialized = .false.


!---------------diagnostics fields------------------------------------- 

integer :: id_diff_cmt, id_udt_cmt, id_vdt_cmt, id_tdt_cmt

character(len=11) :: mod_name = 'cu_mo_trans'

real :: missing_value = -999.

!--------------------- namelist variables with defaults -------------

real    :: c         =   0.7
real    :: diff_norm =   0.1
logical :: diffuse   = .true.
logical :: conserve  = .true.
logical :: add_on    = .true.

namelist/cu_mo_trans_nml/ c, diff_norm, diffuse, conserve, add_on


!--------------------- version number ---------------------------------

character(len=128) :: version = &
'$Id: cu_mo_trans.F90,v 1.2 2003/04/09 20:54:17 fms Exp $'
character(len=128) :: tag = &
'$Name: inchon $'

contains

!#######################################################################

subroutine cu_mo_trans_init( axes, Time )

 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time

integer :: unit, ierr, io

!------ read namelist ------

   if ( file_exist('input.nml')) then
         unit = open_namelist_file ( )
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=cu_mo_trans_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'cu_mo_trans_nml')
      enddo
 10   call close_file (unit)
   endif

!--------- write version number and namelist ------------------

      call write_version_number ( version, tag )
      if ( mpp_pe() == mpp_root_pe() ) &
      write ( stdlog(), nml=cu_mo_trans_nml )

! --- initialize quantities for diagnostics output -------------

   id_tdt_cmt = &
   register_diag_field ( mod_name, 'tdt_cmt', axes(1:3), Time,   &
                        'Temperature tendency from cu_mo_trans', &
                        'deg_K/s', missing_value=missing_value   )
   id_udt_cmt = &
   register_diag_field ( mod_name, 'udt_cmt', axes(1:3), Time,   &
                        'Zonal wind tendency from cu_mo_trans',  &
                        'm/s2', missing_value=missing_value      )

   id_vdt_cmt = &
   register_diag_field ( mod_name, 'vdt_cmt', axes(1:3), Time,       &
                        'Meridional wind tendency from cu_mo_trans', &
                        'm/s2', missing_value=missing_value          )

   id_diff_cmt = &
   register_diag_field ( mod_name, 'diff_cmt', axes(1:3), Time,    &
                        'cu_mo_trans coeff for momentum',  'm2/s', &
                         missing_value=missing_value               )

!--------------------------------------------------------------

  initialized = .true.


end subroutine cu_mo_trans_init

!#######################################################################

subroutine cu_mo_trans_end()

end subroutine cu_mo_trans_end

!#######################################################################

subroutine cu_mo_trans (is, js, Time, dt, mass_flux, u, v, t, q, &
                        p_half, p_full, z_half, z_full,          &
                        dt_u_in, dt_v_in, dt_t_in, diff)

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: is, js

real, intent(in)                      :: dt
real, intent(in)   , dimension(:,:,:) :: mass_flux, u, v, t, q, &
                                         p_half, p_full, z_half, z_full
real, intent(inout), dimension(:,:,:) :: dt_u_in, dt_v_in, dt_t_in
real, intent(out), dimension(:,:,:) :: diff                      

real, dimension(size(u,1),size(u,2),size(u,3)) ::  dt_u, dt_v, dt_t, &
                                                  del_u, del_v, rho
real, dimension(size(u,1),size(u,2))           :: ubot, vbot, zbot, ztop, &
                                                  zeros
real    :: rho, x, xx
integer :: k, nlev
logical :: used

!-----------------------------------------------------------------------

 if (.not.initialized) call error_mesg ('cu_mo_trans',  &
                      'cu_mo_trans_init has not been called.', FATAL)

!-----------------------------------------------------------------------

nlev = size(u,3)

zeros = 0.0

ubot = 0.0
vbot = 0.0
zbot = z_half(:,:,nlev+1)
ztop = z_half(:,:,nlev+1)
  
do k = 2, nlev
  where(mass_flux(:,:,k) .ne. 0.0 .and. mass_flux(:,:,k+1) == 0.0) 
    ubot =      u(:,:,k)        
    vbot =      v(:,:,k)
    zbot = z_half(:,:,k)
  endwhere
  where(mass_flux(:,:,k-1) == 0.0 .and. mass_flux(:,:,k) .ne. 0.0) 
    ztop = z_half(:,:,k)
  endwhere
end do

dt_u = 0.0
dt_v = 0.0

rho  = p_full/(RDGAS*t)

if(diffuse) then
  diff(:,:,1) = 0.0
  do k = 2, nlev
    diff(:,:,k) = diff_norm*mass_flux(:,:,k)*(ztop-zbot)/(rho(:,:,k)*GRAV)
  end do

endif


! --- Extra diagnostics
     if ( id_diff_cmt > 0 ) then
        used = send_data ( id_diff_cmt, diff, Time, is, js, 1 )
     endif

end subroutine cu_mo_trans


!#######################################################################


end module cu_mo_trans_mod

