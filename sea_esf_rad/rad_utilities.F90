 
		     module rad_utilities_mod

use utilities_mod,      only:  open_file, file_exist,    &
                               check_nml_error, error_mesg, &
                               print_version_number, FATAL, NOTE, &
			       WARNING, get_my_pe, close_file


!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!         module containing radiation table search 
!      routines and derived types used in radiation code
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: rad_utilities.F90,v 1.2 2001/07/05 17:32:54 fms Exp $'
   character(len=128)  :: tag     =  '$Name: eugene $'

!---------------------------------------------------------------------
!-------  interfaces --------


public      &
	  rad_utilities_init, define_environment, &
	  locate_in_table, looktab, table_alloc

interface looktab
    module procedure  looktab_type1, looktab_type2, looktab_type3
end interface

interface table_alloc
   module procedure    table1_alloc, table2_alloc, table3_alloc
end interface



public longwave_tables1_type, longwave_tables2_type, &
       longwave_tables3_type, table_axis_type


type longwave_tables1_type
    real, dimension(:,:), pointer      ::  vae, td, md, cd
end type longwave_tables1_type


type longwave_tables2_type
    real, dimension(:,:,:), pointer    ::  vae, td, md, cd
end type longwave_tables2_type


type longwave_tables3_type
     real,  dimension(:,:), pointer    ::  vae, td          
end type longwave_tables3_type


type table_axis_type
  integer :: first_col
  real    :: min_val
  real    :: max_val
  real    :: tab_inc
end type table_axis_type


public longwave_control_type 

type longwave_control_type
    character(len=10)         :: lw_form
    logical                   :: do_cfc
    logical                   :: do_lwaerosol
    logical                   :: do_ckd2p1
    logical                   :: do_ch4_n2o
    logical                   :: do_ch4n2olbltmpint
    logical                   :: do_co2
    logical                   :: do_lwcldemiss
end type longwave_control_type


public shortwave_control_type

type shortwave_control_type
    character(len=10)         :: sw_form
    character(len=12)         :: swaerosol_form
    logical                   :: do_swaerosol
end type shortwave_control_type

public radiation_control_type

type radiation_control_type
    logical                        :: do_diagnostics
    logical                        :: do_totcld_forcing
    logical, dimension(:), pointer :: do_raddg
end type radiation_control_type

public cloudrad_control_type

type cloudrad_control_type
    logical                        :: do_pred_cld_microphys
    logical                        :: do_presc_cld_microphys
    logical                        :: do_no_cld_microphys
end type cloudrad_control_type

public environment_type

type environment_type
    logical                         ::  running_fms
    logical                         ::  running_skyhi
    logical                         ::  using_sky_periphs
    logical                         ::  using_fms_periphs
    logical                         ::  running_gcm     
    logical                         ::  running_standalone
    character(len=11)               ::  column_type
end type environment_type

!---------------------------------------------------------------------
!-------- namelist  ---------


real   ::  dummy = 1.0

namelist / rad_utilities_nml /   &
                                       dummy


!---------------------------------------------------------------------
!------- public data ------


type (longwave_control_type),  public   ::    &
   Lw_control = longwave_control_type( '          ', &
				      .false., .false., .false.,  &
				      .false., .false., .true., .false.)

type (shortwave_control_type), public   ::  &
   Sw_control = shortwave_control_type( '          ', '            ', &
					.false.  )

type (radiation_control_type), public   ::  Rad_control

type (cloudrad_control_type), public   ::   &
   Cldrad_control = cloudrad_control_type(.false., .false., .false.)

type (environment_type), public   ::   &
   Environment = environment_type(.false., .false., .false., &
				  .false., .false., .false., &
				  '      ')

!---------------------------------------------------------------------
!------- private data ------



!---------------------------------------------------------------------
!---------------------------------------------------------------------




contains





subroutine rad_utilities_init

!------------------------------------------------------------------
      integer    ::  unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=rad_utilities_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'rad_utilities_nml')
       enddo
10     call close_file (unit)
     endif

     unit = open_file ('logfile.out', action='append')
!    call print_version_number (unit, 'rad_utilities', version_number)
     if (get_my_pe() == 0) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=rad_utilities_nml)
     endif
     call close_file (unit)

!------------------------------------------------------------------

end  subroutine rad_utilities_init



!####################################################################

subroutine define_environment (model_name, peripherals_source,  &
                               application_type, column_type_in)

!--------------------------------------------------------------------
character(len=*),  intent(in)    :: model_name, peripherals_source, &
				    application_type, &
				    column_type_in
!---------------------------------------------------------------------

   if (trim(model_name) == 'skyhi') then
     Environment%running_skyhi = .true.
     Environment%running_fms   = .false.
   else if (trim(model_name) == 'fms') then
     Environment%running_skyhi = .false.
     Environment%running_fms   = .true.
   else
     call error_mesg ('define_environment', &
       ' model_name not recognized', FATAL)
   endif

   if (trim(peripherals_source) == 'skyhi')  then
     Environment%using_sky_periphs = .true.
     Environment%using_fms_periphs = .false.
   else if (trim(peripherals_source) == 'fms') then
     Environment%using_sky_periphs = .false.
     Environment%using_fms_periphs = .true.
   else
     call error_mesg ('define_environment', &
       ' peripherals_source not recognized', FATAL)
   endif

   if (trim(application_type) == 'gcm') then
     Environment%running_gcm = .true.
     Environment%running_standalone = .false.
   else if (trim(application_type) == 'standalone') then
     Environment%running_gcm = .false.
     Environment%running_standalone = .true.
     if (Environment%using_fms_periphs)   call error_mesg &
       ( 'define_environment', &
         'currently must use skyhi peripherals with standalone', FATAL)
   else
     call error_mesg ('define_environment', &
       ' application_type not acceptable', FATAL)
   endif

   Environment%column_type = column_type_in

!-------------------------------------------------------------------

end subroutine define_environment



!####################################################################

subroutine locate_in_table (table_axis, x, dx, ix, k_min, k_max) 

!-----------------------------------------------------------------------
!
!     given array x and an arithmetic sequence of table column headings
!     tabxmin, tabxmin+tabdeltax, ..., corresponding to column ixlow, 
!     ixlow+1, ..., ixupp, Locate returns the array ix is column 
!     indices and the array dx of residuals.
!
!     author: c. h. goldberg
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!
!-----------------------------------------------------------------------
type(table_axis_type),     intent(in)  :: table_axis
real,    dimension(:,:,:), intent(in)  :: x
integer,                   intent(in)  :: k_min, k_max
real,    dimension(:,:,:), intent(out) :: dx
integer, dimension(:,:,:), intent(out) :: ix

!-----------------------------------------------------------------------
real, dimension (size(x,1), size(x,2), size(x,3))  ::  fx
real                                               ::  table_min,  &
						       table_inc
integer                                            ::  k, table_col


      do k=k_min,k_max
        fx (:,:,k) = AINT((x(:,:,k) - table_axis%min_val )/  &
				      table_axis%tab_inc)
        ix (:,:,k) = INT(fx(:,:,k)) + table_axis%first_col
        dx (:,:,k) = x(:,:,k) - fx(:,:,k)*table_axis%tab_inc - &
				      table_axis%min_val
      end do
      
!---------------------------------------------------------------------


end subroutine locate_in_table



!####################################################################

subroutine looktab_type1 (tab, ix, iy, dx, dy, answer, k_min, k_max)   

!----------------------------------------------------------------------
!
!     given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!     arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!     y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!     from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!
!     author: c. h. goldberg
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!     
!--------------------------------------------------------------------
type(longwave_tables1_type), intent(in)  :: tab
integer,dimension(:,:,:),    intent(in)  :: ix, iy
real,   dimension(:,:,:),    intent(in)  :: dx, dy   
real,   dimension(:,:,:),    intent(out) :: answer
integer,                     intent(in)  :: k_min, k_max

!--------------------------------------------------------------------
      integer    ::  i_min, i_max, j_min, j_max, i,j, k
  
      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
 	do j=j_min, j_max
          do i=i_min, i_max
	    answer(i,j,k) =                                          &
                                    tab%vae (ix(i,j,k), iy(i,j,k)) + &
                          dx(i,j,k)*tab%td  (ix(i,j,k), iy(i,j,k)) + &
                          dy(i,j,k)*tab%md  (ix(i,j,k), iy(i,j,k)) + &
                dx(i,j,k)*dy(i,j,k)*tab%cd  (ix(i,j,k), iy(i,j,k))
	  end do
        end do
      end do
!---------------------------------------------------------------------


end subroutine looktab_type1



!#####################################################################

subroutine looktab_type2 (tab, ix, iy, dx, dy, answer, k_min, k_max, m)

!-------------------------------------------------------------------
!
!     given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!     arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!     y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!     from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!
!     author: c. h. goldberg
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!     
!--------------------------------------------------------------------
type(longwave_tables2_type), intent(in)   :: tab
integer, dimension (:,:,:),  intent(in)   :: ix, iy
integer,                     intent(in)   :: m
real, dimension (:,:,:),     intent(in)   :: dx, dy
real, dimension (:,:,:),     intent(out)  :: answer
integer,                     intent(in)   :: k_min, k_max

!--------------------------------------------------------------------
       integer    ::    i, j, k, i_min, i_max, j_min, j_max
  
!-------------------------------------------------------------------

      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
	do j=j_min, j_max
          do i=i_min, i_max
	    answer(i,j,k) =                                           &
                                   tab%vae (ix(i,j,k), iy(i,j,k),m) + &
                         dx(i,j,k)*tab%td (ix(i,j,k), iy(i,j,k),m) + &
                         dy(i,j,k)*tab%md (ix(i,j,k), iy(i,j,k),m) + &
               dx(i,j,k)*dy(i,j,k)*tab%cd   (ix(i,j,k), iy(i,j,k),m)
	  end do
        end do
      end do

!--------------------------------------------------------------------

end subroutine looktab_type2



!###################################################################

subroutine looktab_type3 (tab, ix, dx,  answer, k_min, k_max, n)

!----------------------------------------------------------------------
!
!     given arrays ix(:,:,:) and dx(:,:,:) of integer subscripts and!
!     differences from x(:,:,:) and constant column subscript iyconst, 
!     calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:)) from four tables
!     of values f, df/dx, df/dy, and d2f/dxdy.
!
!     author: c. h. goldberg
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!     
!-----------------------------------------------------------------------
type(longwave_tables3_type), intent(in)  :: tab
integer, dimension (:,:,:),  intent(in)  :: ix
integer,                     intent(in)  :: n
real,    dimension(:,:,:),   intent(in)  :: dx
real,    dimension(:,:,:),   intent(out) :: answer 
integer,                     intent(in)  :: k_min, k_max

!----------------------------------------------------------------------
      integer    :: i, j, k, i_min, i_max, j_min, j_max
   
!-----------------------------------------------------------------
      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
	do j=j_min, j_max
          do i=i_min, i_max
	    answer(i,j,k) =                                 &
                                      tab%vae (ix(i,j,k),n) +   &
                            dx(i,j,k)*tab%td(ix(i,j,k),n)
	  end do
        end do
      end do

!------------------------------------------------------------------

end subroutine  looktab_type3


!#####################################################################

subroutine table1_alloc (tab, dim1, dim2)

type(longwave_tables1_type), intent (inout) :: tab
integer,                     intent(in)  :: dim1, dim2

      allocate (tab%vae(dim1, dim2))
      allocate (tab%td (dim1, dim2))
      allocate (tab%md (dim1, dim2))
      allocate (tab%cd (dim1, dim2))

end subroutine table1_alloc

!####################################################################

subroutine table2_alloc (tab, dim1, dim2, dim3)

type(longwave_tables2_type), intent (inout) :: tab
integer,                     intent(in)  :: dim1, dim2, dim3

      allocate (tab%vae(dim1, dim2, dim3))
      allocate (tab%td (dim1, dim2, dim3))
      allocate (tab%md (dim1, dim2, dim3))
      allocate (tab%cd (dim1, dim2, dim3))

end subroutine table2_alloc

!#####################################################################

subroutine table3_alloc (tab, dim1, dim2)

type(longwave_tables3_type), intent (inout) :: tab
integer,                     intent(in)  :: dim1, dim2

      allocate (tab%vae(dim1, dim2))
      allocate (tab%td (dim1, dim2))

end subroutine table3_alloc


!##################################################################




		    end module rad_utilities_mod


