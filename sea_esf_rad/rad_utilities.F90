 
		     module rad_utilities_mod

use utilities_mod,      only:  open_file, file_exist,    &
                               check_nml_error, error_mesg, &
                               print_version_number, FATAL, NOTE, &
			       get_num_pes, &
			       WARNING, get_my_pe, close_file
!use constants_mod,      only : radian


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
   character(len=128)  :: version =  '$Id: rad_utilities.F90,v 1.4 2003/04/09 21:01:26 fms Exp $'
   character(len=128)  :: tag     =  '$Name: inchon $'

!---------------------------------------------------------------------
!-------  interfaces --------


public      &
	  rad_utilities_init, define_environment, &
	  map_global_indices,   &
          locate_in_table, looktab, table_alloc
!                   looktab, table_alloc

interface looktab
    module procedure  looktab_type1, looktab_type2, looktab_type3
end interface

interface table_alloc
   module procedure    table1_alloc, table2_alloc, table3_alloc
end interface



public longwave_tables1_type, longwave_tables2_type, &
       longwave_tables3_type, table_axis_type


type longwave_tables1_type
    real, dimension(:,:), pointer      ::  vae=>NULL(), td=>NULL(), &
                                           md=>NULL(), cd=>NULL()
end type longwave_tables1_type


type longwave_tables2_type
    real, dimension(:,:,:), pointer    ::  vae=>NULL(), td=>NULL(),  &
                                           md=>NULL(), cd=>NULL()
end type longwave_tables2_type


type longwave_tables3_type
     real,  dimension(:,:), pointer    ::  vae=>NULL(), td=>NULL()          
end type longwave_tables3_type


type table_axis_type
  integer :: first_col
  real    :: min_val
  real    :: max_val
  real    :: tab_inc
end type table_axis_type


public longwave_control_type 

type longwave_control_type
!    character(len=10)         :: lw_form
    character(len=16)         :: lw_form
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
!   character(len=10)         :: sw_form
!   character(len=16)         :: sw_form
    logical                   :: do_lhsw
    logical                   :: do_esfsw
!   character(len=12)         :: swaerosol_form
    character(len=16)         :: swaerosol_form
    logical                   :: do_swaerosol
end type shortwave_control_type

public radiation_control_type

type radiation_control_type
!   logical                        :: do_diagnostics
    logical                        :: do_totcld_forcing
    logical                        :: do_aerosol
!   logical, dimension(:), pointer :: do_raddg
end type radiation_control_type

!------------------------------------------------------------------

public cloudrad_control_type

type cloudrad_control_type
    logical                        :: do_pred_cld_microphys
    logical                        :: do_presc_cld_microphys
    logical                        :: do_no_cld_microphys
    logical                        :: do_rh_clouds        
    logical                        :: do_strat_clouds        
    logical                        :: do_zonal_clouds        
    logical                        :: do_mgroup_prescribed
    logical                        :: do_obs_clouds        
    logical                        :: do_no_clouds        
    logical                        :: do_diag_clouds        
    logical                        :: do_specified_clouds        
    logical                        :: do_donner_deep_clouds
    integer                        :: nlwcldb
end type cloudrad_control_type

!------------------------------------------------------------------

public environment_type

type environment_type
!   logical                         ::  running_fms
!   logical                         ::  running_skyhi
    logical                         ::  using_sky_periphs
    logical                         ::  using_fms_periphs
    logical                         ::  running_gcm     
    logical                         ::  running_standalone
    logical                         ::  running_sa_model
!   character(len=11)               ::  column_type
    character(len=16)               ::  column_type
end type environment_type

!------------------------------------------------------------------

public   radiative_gases_type
 
type radiative_gases_type
     real     :: rrvch4, rrvn2o, rrvco2,    &
                 rrvf11, rrvf12, rrvf113, rrvf22, &
		 rf11air, rf12air, rf113air, rf22air
     real, dimension(:,:,:), pointer :: qo3=>NULL()
!    logical  :: do_co2, do_ch4_n2o, do_cfc
     logical  :: time_varying_co2, time_varying_f11, &
                 time_varying_f12, time_varying_f113, &
                 time_varying_f22,  &
                 time_varying_ch4, time_varying_n2o
end type radiative_gases_type

!------------------------------------------------------------------

public   aerosol_type              
  
type aerosol_type
     real, dimension(:,:,:,:), pointer :: aerosol=>NULL()
     integer                  :: max_data_fields 
!     character(len=64), dimension(:), pointer :: data_names
     integer    :: nfields
     character(len=64), dimension(:), pointer ::   &
                                          aerosol_optical_names=>NULL()
     integer, dimension(:), pointer :: sulfate_index=>NULL()
     integer, dimension(:), pointer :: optical_index=>NULL()
end type aerosol_type              
!------------------------------------------------------------------
 
public   aerosol_properties_type

type aerosol_properties_type
     integer    :: NAERMODELS
     integer    :: num_wavenumbers
     integer    :: nfields
     integer, dimension(:), pointer :: endaerwvnsf=>NULL()
     character(len=64), dimension(:), pointer :: aerosol_names=>NULL()
     real, dimension(:,:), pointer :: aeroextivl=>NULL(),  &
                                      aerossalbivl=>NULL(), &
                                      aeroasymmivl=>NULL()
     real, dimension(:,:), pointer :: aerextband=>NULL(),   &
                                      aerssalbband=>NULL(), &
                                      aerasymmband=>NULL()
end type aerosol_properties_type

!------------------------------------------------------------------

public optical_path_type

type optical_path_type
     real, dimension (:,:,:,:), pointer :: empl1f=>NULL(),  &
                                           empl2f=>NULL(),  &
                                           vrpfh2o=>NULL(), &
!                                          totch2o, totch2obd, &
                                           xch2obd=>NULL(),  &
                                           tphfh2o=>NULL(), &
                                           avephif=>NULL(), &
                                           totaerooptdep=>NULL()
     real, dimension (:,:,:), pointer   :: empl1=>NULL(), &
                                           empl2=>NULL(),  &
					   var1=>NULL(), &
					   var2=>NULL(), &
                                          emx1f=>NULL(),   &
					  emx2f=>NULL(),   &
					  totvo2=>NULL(),  &
					  avephi=>NULL(),&
                                         totch2obdwd=>NULL(), &
					 xch2obdwd=>NULL(), &
                                           totphi=>NULL(),   &
					   cntval=>NULL(), &
					   toto3=>NULL(),   &
!                                        tphio3, var3, var4, sh2o,  &
!                                          tmpexp, rvh2o, wk, rhoave, &
                                         tphio3=>NULL(),  &
					 var3=>NULL(),  &
					 var4=>NULL(),        &
                                         wk=>NULL(),         &
                                         rh2os=>NULL(),  &
					 rfrgn=>NULL(),  &
					 tfac=>NULL(), &
                                         totaerooptdep_15=>NULL(), &
                                         totf11=>NULL(),   &
					 totf12=>NULL(),  &
					 totf113=>NULL(),   &
					 totf22=>NULL()
      real, dimension (:,:), pointer     :: emx1=>NULL(),  &
                                            emx2=>NULL(),  &
					    csfah2o=>NULL(), &
                                            aerooptdep_KE_15=>NULL()
end type optical_path_type

!------------------------------------------------------------------

public gas_tf_type

type gas_tf_type
     real, dimension(:,:,:), pointer :: tdav=>NULL(),   &
                                        tlsqu=>NULL(),   &
                                        tmpdiff=>NULL(),   &
                                        tstdav=>NULL(),  &
                                        co2nbl=>NULL(),   &
                                        n2o9c=>NULL(),   &
                                        tn2o17=>NULL()
     real, dimension(:,:,:,:), pointer :: co2spnb=>NULL()
     real, dimension(:,:), pointer :: a1=>NULL(), a2=>NULL()
end type gas_tf_type

!------------------------------------------------------------------

public lw_output_type

type lw_output_type
     real, dimension(:,:,:), pointer :: heatra=>NULL(), &
                                        flxnet=>NULL(),  &
                                        heatracf=>NULL(), &
                                        flxnetcf=>NULL()
end type lw_output_type

!------------------------------------------------------------------

public lw_clouds_type

type lw_clouds_type
     real, dimension(:,:,:,:),   pointer :: taucld_rndlw=>NULL(), &
                                            taucld_mxolw=>NULL(), &
                                            taunbl_mxolw=>NULL()
end type lw_clouds_type

!------------------------------------------------------------------

public lw_diagnostics_type

type lw_diagnostics_type
     real, dimension(:,:),   pointer :: flx1e1=>NULL(), gxcts=>NULL()
     real, dimension(:,:,:), pointer :: flx1e1f=>NULL(), excts=>NULL(),&
                                        fctsg=>NULL()
     real, dimension(:,:,:,:), pointer :: fluxn=>NULL(),   &
                                          fluxncf=>NULL(),   &
					  exctsn=>NULL(),  &
                                          cts_out=>NULL(), &
					  cts_outcf=>NULL()
end type lw_diagnostics_type

!------------------------------------------------------------------

public lw_table_type

type lw_table_type
     real, dimension(:),   pointer :: bdlocm=>NULL(),   &
                                      bdhicm=>NULL(),  &
				      bandlo=>NULL(),  &
				      bandhi=>NULL()
     integer, dimension(:), pointer :: iband=>NULL()
end type lw_table_type

!------------------------------------------------------------------

public cld_space_properties_type

type cld_space_properties_type
     real, dimension(:,:,:),   pointer :: camtswkc=>NULL()        
!    real, dimension(:,:,:,:),   pointer :: cirabswkc, cirrfswkc, &
     real, dimension(:,:,:),   pointer :: cirabswkc=>NULL(),  &
                                          cirrfswkc=>NULL(), &
                                            cvisrfswkc=>NULL()
     integer, dimension(:,:,:), pointer :: ktopswkc=>NULL(),   &
                                            kbtmswkc=>NULL()
end type cld_space_properties_type

!------------------------------------------------------------------

public cld_diagnostics_type

type cld_diagnostics_type
     real, dimension(:,:,:),   pointer :: lwpath=>NULL(),   &
                                          iwpath=>NULL(),   &
					  size_drop=>NULL(), &
                                          size_ice=>NULL()
     real, dimension(:,:),     pointer :: cld_isccp_hi=>NULL(),   &
                                          cld_isccp_mid=>NULL(), &
                                          cld_isccp_low=>NULL(),   &
                                          tot_clds=>NULL()
end type cld_diagnostics_type

!------------------------------------------------------------------

public sw_output_type

type sw_output_type
     real, dimension(:,:,:), pointer :: dfsw=>NULL(), ufsw=>NULL(),  &
                                        fsw=>NULL(), hsw=>NULL(), &
                                        dfswcf=>NULL(), ufswcf=>NULL(),&
                                        fswcf=>NULL(), hswcf=>NULL()
end type sw_output_type

!-------------------------------------------------------------------

public atmos_input_type

type atmos_input_type
     real, dimension(:,:,:), pointer :: press=>NULL(),   &
                                        temp=>NULL(), &
					rh2o=>NULL(),  &
					pflux=>NULL(), &
                                        tflux=>NULL(),  &
					deltaz=>NULL(),  &
					cloud_ice=>NULL(), &
                                        cloud_water=>NULL(), &
					phalf=>NULL(),   &
					rel_hum=>NULL(), &
                                        cloudtemp=>NULL(),   &
					clouddeltaz=>NULL(), &
                                        cloudvapor=>NULL()
     real, dimension(:,:),   pointer :: psfc=>NULL(), tsfc=>NULL(),   &
                                        asfc=>NULL(), land=>NULL()
endtype atmos_input_type

!-------------------------------------------------------------------

public fsrad_output_type

type fsrad_output_type
     real, dimension(:,:,:), pointer :: tdtsw=>NULL(), &
                                        tdtlw=>NULL(),  &
					tdtsw_clr=>NULL(),  &
                                       tdtlw_clr=>NULL()
     real, dimension(:,:),   pointer :: swdns=>NULL(),   &
                                        swups=>NULL(),  &
					lwups=>NULL(), &
					lwdns=>NULL(), &
                                       swin=>NULL(), &
				       swout=>NULL(), &
				       olr=>NULL(), &
                                       swdns_clr=>NULL(),  &
				       swups_clr=>NULL(),  &
				       lwups_clr=>NULL(),&
                                       lwdns_clr=>NULL(),   &
				       swin_clr=>NULL(),  &
				       swout_clr=>NULL(), &
                                       olr_clr=>NULL()
     integer      :: npass
end type fsrad_output_type

!-------------------------------------------------------------------

public rad_output_type

type rad_output_type
     real, dimension(:,:,:), pointer :: tdt_rad=>NULL(),  &
                                        tdt_rad_clr=>NULL(), &
                                        tdtsw=>NULL(),   &
                                        tdtsw_clr=>NULL(),  &
                                        tdtlw=>NULL()
     real, dimension(:,:), pointer :: flux_sw_surf=>NULL(), &
                                      flux_lw_surf=>NULL(), &
                                      coszen_angle=>NULL()
end type rad_output_type

!-------------------------------------------------------------------

public astronomy_type
    
type astronomy_type
     logical :: do_diurnal, do_annual, do_daily_mean
     real    :: rrsun, solar_constant
     real, dimension(:,:), pointer  :: solar=>NULL(),   &
                                       cosz=>NULL(),  &
				       fracday=>NULL()
end type astronomy_type

!--------------------------------------------------------------------

public cldrad_properties_type

type cldrad_properties_type
     real, dimension(:,:,:,:), pointer :: cldext=>NULL(),   &
                                          cldasymm=>NULL(), &
					  cldsct=>NULL(), &
                                          emmxolw=>NULL(), &
					  emrndlw=>NULL(),  &
					  cirabsw=>NULL(), &
                                          abscoeff=>NULL(),  &
					  cldemiss=>NULL(), &
                                          cirrfsw=>NULL(),   &
					  cvisrfsw=>NULL()
     real, dimension(:,:,:), pointer :: camtsw=>NULL(),  &
                                        cmxolw=>NULL(),  &
					crndlw=>NULL()      
     Integer, dimension(:,:), pointer :: ncldsw=>NULL(),  &
                                         nmxolw=>NULL(),  &
					 nrndlw=>NULL()     
!     integer                          :: NLWCLDB
end type cldrad_properties_type


public longwave_parameter_type


type longwave_parameter_type
     integer             :: offset
     integer             :: NBTRG
     integer             :: NBTRGE
!     integer             :: NLWCLDB
     integer             :: NBLY
     integer             :: n_lwaerosol_bands
end type longwave_parameter_type

!---------------------------------------------------------------------
!-------- namelist  ---------


!real   ::  dummy = 1.0
!character(len=12)              ::                    &
character(len=16)              ::                    &
                                   peripherals_source=' ', &
                                   application_type = '  ', &
                                   column_type = '    '


namelist / rad_utilities_nml /   &
!                                      dummy
                                   peripherals_source, &
                                   application_type, &
                                   column_type


!---------------------------------------------------------------------
!------- public data ------


type (longwave_control_type),  public   ::    &
   Lw_control = longwave_control_type( '          ', &
				      .false., .false., .false.,  &
				      .false., .false., .true., .false.)

type (shortwave_control_type), public   ::  &
   Sw_control = shortwave_control_type( .false., .false.,    &
                                        '            ',  .false.  )

type (radiation_control_type), public   ::  Rad_control

type (aerosol_properties_type), public   ::  Aerosol_props

type (cloudrad_control_type), public   ::   &
   Cldrad_control = cloudrad_control_type(.false., .false., .false., &
                                          .false., .false., .false., &
					  .false., .false., .false., &
					  .false., .false., .false., &
					  0 )

type (environment_type), public   ::   &
!  Environment = environment_type(.false., .false., .false., &
   Environment = environment_type(                  .false., &
				  .false., .false., .false., &
                                  .false.,  '      ')

type (longwave_parameter_type), public  ::   &
!   Lw_parameters = longwave_parameter_type(0, 0, 0, 0, 0)
   Lw_parameters = longwave_parameter_type(0, 0, 0, 0, 0)

type (table_axis_type),        public   ::    &
     temp_1 = table_axis_type(1, 100.0, 370.0, 10.0), &
     mass_1 = table_axis_type(1, -16.0,   1.9,  0.1)

!---------------------------------------------------------------------
!------- private data ------


logical :: rad_utilities_initialized=.false.

!---------------------------------------------------------------------
!---------------------------------------------------------------------




contains





subroutine rad_utilities_init

!------------------------------------------------------------------
      integer    ::  unit, ierr, io


      if (rad_utilities_initialized) return
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

     call define_environment


     rad_utilities_initialized = .true.

!------------------------------------------------------------------

end  subroutine rad_utilities_init



!####################################################################

!subroutine define_environment (model_name, peripherals_source,  &
subroutine define_environment            
!subroutine define_environment (            peripherals_source,  &
!                               application_type, column_type_in)

!--------------------------------------------------------------------
!character(len=*),  intent(in)    :: model_name, peripherals_source, &
!haracter(len=*),  intent(in)    ::             peripherals_source, &
!			    application_type, &
!			    column_type_in
!---------------------------------------------------------------------

!  if (trim(model_name) == 'skyhi') then
!    Environment%running_skyhi = .true.
!    Environment%running_fms   = .false.
!  else if (trim(model_name) == 'fms') then
!    Environment%running_skyhi = .false.
!    Environment%running_fms   = .true.
!  else
!    call error_mesg ('define_environment', &
!      ' model_name not recognized', FATAL)
!  endif

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
   else
     call error_mesg ('define_environment', &
       ' application_type not acceptable', FATAL)
   endif

!  Environment%column_type = column_type_in
   Environment%column_type = column_type

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


!subroutine map_global_indices (jmax_aerfile, latb, jindx2)
subroutine map_global_indices (global_rows , latb,    &
                               global_index_array )

integer, intent(in)   :: global_rows 
real, dimension(:), intent(in) :: latb
integer, dimension(:), intent(out) :: global_index_array

       real   :: lat_start
       real, dimension(global_rows) :: data_lat
       integer   :: j, jst, jj

       if (get_num_pes()*(size(latb,1) -1) == global_rows) then
!       lat_start = -90.0 + 180./(2.*float(global_rows))
        lat_start = -0.5*acos(-1.0) + acos(-1.0)/(2.*float(global_rows))
       do j=1,global_rows
!        data_lat(j) = lat_start + real(j-1)*180.0/real(global_rows)
         data_lat(j) = lat_start + real(j-1)*acos(-1.0)/real(global_rows)
       end do
       jst = 1
       do jj=1, size(latb,1) - 1
         do j = jst,global_rows
!   if (data_lat(j) >= latb(jj)*radian ) then
	   if (data_lat(j) >= latb(jj)                    ) then
	     global_index_array(jj) = j
	     jst = j
             exit
           endif
         end do
       end do
       else
         call error_mesg ('esfsw_scattering_mod', &
            'resolution of data input file doesn''t match model size --&
	    &  must provide inter(extra)polation not yet implemented', &
            FATAL)
       endif

end subroutine map_global_indices 


!####################################################################


		    end module rad_utilities_mod




