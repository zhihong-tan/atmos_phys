
module damping_driver_mod

!-----------------------------------------------------------------------
!
!       This module controls two functions:
!
!   (1) rayleigh friction applied to momentum fields at levels
!       1 to kbot (i.e., momentum is damped toward zero).
!
!   (2) mountain gravity wave drag module may be called
!
!-----------------------------------------------------------------------

 use      mg_drag_mod, only:  mg_drag, mg_drag_init, mg_drag_end
 use    utilities_mod, only:  file_exist, open_file, error_mesg,     &
                              check_nml_error, print_version_number, &
                              get_my_pe, FATAL, close_file
 use diag_manager_mod, only:  register_diag_field, send_data
 use time_manager_mod, only:  time_type

 implicit none
 private

 public   damping_driver, damping_driver_init, damping_driver_end

!-----------------------------------------------------------------------
!---------------------- namelist ---------------------------------------

   real     :: trayfric = 0.
   integer  :: nlev_rayfric = 1
   logical  :: do_mg_drag = .false.

   namelist /damping_driver_nml/  trayfric,   nlev_rayfric,  &
                                  do_mg_drag

!
!   trayfric = damping time in seconds for rayleigh damping momentum
!              in the top nlev_rayfric layers (if trayfric < 0 then time
!              in days)
!                 [real, default: trayfric=0.]
!
!   nlev_rayfric = number of levels at the top of the model where
!                  rayleigh friction of momentum is performed, if
!                  trayfric=0. then nlev_rayfric has no effect
!                    [integer, default: nlev_rayfric=1]
!
!-----------------------------------------------------------------------
!----- id numbers for diagnostic fields -----

integer :: id_udt_rdamp,  id_vdt_rdamp,                 &
           id_udt_gwd,    id_vdt_gwd,    id_taub

!----- missing value for all fields ------

real :: missing_value = -999.

character(len=7) :: mod_name = 'damping'

!-----------------------------------------------------------------------

 logical :: do_rayleigh

 real, parameter ::  daypsec=1./86400.
 logical :: do_init=.true.

 real :: rfactr

!   note:  
!     rfactr = coeff. for damping momentum at the top level

 character(len=4) :: vers_num = 'v2.1'

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine damping_driver (is, js, Time, pfull, phalf, zfull, zhalf, &
                            u, v, t, q, r,  udt, vdt, tdt, qdt, rdt,  &
                            mask, kbot)
 
!-----------------------------------------------------------------------
 integer,         intent(in)                :: is, js
 type(time_type), intent(in)                :: Time
 real,    intent(in),    dimension(:,:,:)   :: pfull, phalf, &
                                               zfull, zhalf, &
                                               u, v, t, q
 real,    intent(in),    dimension(:,:,:,:) :: r
 real,    intent(inout), dimension(:,:,:)   :: udt,vdt,tdt,qdt
 real,    intent(inout), dimension(:,:,:,:) :: rdt
 real,    intent(in),    dimension(:,:,:), optional :: mask
 integer, intent(in),    dimension(:,:),   optional :: kbot
!-----------------------------------------------------------------------
 real, dimension(size(udt,1),size(udt,2))             :: taub
 real, dimension(size(udt,1),size(udt,2),size(udt,3)) :: utnd, vtnd, p2
 logical :: used
!-----------------------------------------------------------------------

   if (do_init) call error_mesg ('damping_driver',  &
                     'damping_driver_init must be called first', FATAL)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!----------------- r a y l e i g h   d a m p i n g ---------------------
!-----------------------------------------------------------------------
   if (do_rayleigh) then

       p2 = pfull * pfull
       call rayleigh (p2, u, v, utnd, vtnd)
       udt = udt + utnd
       vdt = vdt + vtnd

!----- diagnostics -----

       if ( id_udt_rdamp > 0 ) then
            used = send_data ( id_udt_rdamp, utnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_vdt_rdamp > 0 ) then
            used = send_data ( id_vdt_rdamp, vtnd, Time, is, js, 1, &
                               rmask=mask )
       endif

   endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--------- m t n   g r a v i t y   w a v e   d r a g -------------------
!-----------------------------------------------------------------------
   if (do_mg_drag) then

       call mg_drag (is, js, u, v, t, pfull, phalf, zfull, zhalf,  &
                     utnd, vtnd, taub, kbot)
       udt = udt + utnd
       vdt = vdt + vtnd

!----- diagnostics -----

       if ( id_udt_gwd > 0 ) then
            used = send_data ( id_udt_gwd, utnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_vdt_gwd > 0 ) then
            used = send_data ( id_vdt_gwd, vtnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_taub > 0 ) then
            used = send_data ( id_taub, taub, Time, is, js )
       endif

   endif
!-----------------------------------------------------------------------

 end subroutine damping_driver

!#######################################################################

 subroutine damping_driver_init ( nlon, nlat, axes, Time )

 integer,         intent(in) :: nlon, nlat, axes(4)
 type(time_type), intent(in) :: Time
!-----------------------------------------------------------------------
 integer :: unit, ierr, io

!-----------------------------------------------------------------------
!----------------- namelist (read & write) -----------------------------

   if (file_exist('input.nml')) then
      unit = open_file ('input.nml', action='read')
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=damping_driver_nml, iostat=io, end=10)
         ierr = check_nml_error (io, 'damping_driver_nml')
      enddo
 10   call close_file (unit)
   endif

   unit = open_file ('logfile.out', action='append')
   call print_version_number (unit, 'damping_driver', vers_num)
   if ( get_my_pe() == 0 ) write (unit,nml=damping_driver_nml)
   call close_file (unit)

!-----------------------------------------------------------------------
!--------- rayleigh friction ----------

   do_rayleigh=.false.

   if (abs(trayfric) > 0.0001 .and. nlev_rayfric > 0) then
      if (trayfric > 0.0) then
         rfactr=(1./trayfric)
      else
         rfactr=(1./abs(trayfric))*daypsec
      endif
         do_rayleigh=.true.
   else
         rfactr=0.0
   endif

!-----------------------------------------------------------------------
!----- mountain gravity wave drag -----

   if (do_mg_drag) call mg_drag_init (nlon, nlat, ierr)

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----

if (do_rayleigh) then

   id_udt_rdamp = &
   register_diag_field ( mod_name, 'udt_rdamp', axes(1:3), Time,       &
                       'u wind tendency for Rayleigh damping', 'm/s2', &
                        missing_value=missing_value               )

   id_vdt_rdamp = &
   register_diag_field ( mod_name, 'vdt_rdamp', axes(1:3), Time,       &
                       'v wind tendency for Rayleigh damping', 'm/s2', &
                        missing_value=missing_value               )
endif

if (do_mg_drag) then

   id_udt_gwd = &
   register_diag_field ( mod_name, 'udt_gwd', axes(1:3), Time,        &
                     'u wind tendency for gravity wave drag', 'm/s2', &
                        missing_value=missing_value               )

   id_vdt_gwd = &
   register_diag_field ( mod_name, 'vdt_gwd', axes(1:3), Time,        &
                     'v wind tendency for gravity wave drag', 'm/s2', &
                        missing_value=missing_value               )

   id_taub = &
   register_diag_field ( mod_name, 'taub', axes(1:2), Time,        &
                     'base flux for gravity wave drag', 'kg/m/s2', &
                        missing_value=missing_value               )
endif

!-----------------------------------------------------------------------

   do_init=.false.

!******************** end of initialization ****************************
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 end subroutine damping_driver_init

!#######################################################################

 subroutine damping_driver_end

     if (do_mg_drag) call mg_drag_end

 end subroutine damping_driver_end

!#######################################################################

 subroutine rayleigh (p2, u, v, udt, vdt)

  real,    intent(in),  dimension(:,:,:)   :: p2, u, v
  real,    intent(out), dimension(:,:,:)   :: udt, vdt

  real, dimension(size(u,1),size(u,2)) :: fact
  integer :: k
!-----------------------------------------------------------------------
!--------------rayleigh damping of momentum (to zero)-------------------

   do k = 1, nlev_rayfric
     fact(:,:) = rfactr*(1.+(p2(:,:,1)-p2(:,:,k))/(p2(:,:,1)+p2(:,:,k)))
     udt(:,:,k) = -u(:,:,k)*fact(:,:)
     vdt(:,:,k) = -v(:,:,k)*fact(:,:)
   enddo

   do k = nlev_rayfric+1, size(u,3)
     udt(:,:,k) = 0.0
     vdt(:,:,k) = 0.0
   enddo

!-----------------------------------------------------------------------

 end subroutine rayleigh

!#######################################################################

end module damping_driver_mod

