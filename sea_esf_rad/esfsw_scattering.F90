
		   module esfsw_scattering_mod


!use rad_step_setup_mod,    only: ISRAD, IERAD, JSRAD,JERAD,  &
!			         KSRAD, KERAD, iabs, jabs
!use rad_step_setup_mod,    only:  iabs, jabs
use  utilities_mod,        only: open_file, file_exist,    &
                                 check_nml_error, error_mesg, &
                                 get_num_pes,  &
                                 print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file
!			 WARNING, get_my_pe, close_file, &
!			 get_domain_decomp
use esfsw_parameters_mod,  only: nbands, get_solarfluxes, &
			         TOT_WVNUMS
use rad_utilities_mod,     only: shortwave_control_type, Sw_control, &
			         Environment, environment_type, &
				 map_global_indices
use microphys_rad_mod,     only: thickavg, thinavg
use constants_new_mod, only: radians_to_degrees

!-------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!         shortwave aerosol scattering parameterization module
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.09'
    character(len=128)  :: version =  '$Id: esfsw_scattering.F90,v 1.3 2002/07/16 22:35:10 fms Exp $'
    character(len=128)  :: tag     =  '$Name: havana $'



!---------------------------------------------------------------------
!-------  interfaces --------

public           &
       esfsw_scattering_init, scatpar, get_swaerprops_from_scattering
     

private          &
       aeropar, aerosf


!---------------------------------------------------------------------
!-------- namelist  ---------

!character(len=19)            :: swaer_model_type = '                  '
!character(len=12)            :: swaerosol_form = '            '
character(len=48)            :: swaer_model_type = '                  '
character(len=16)            :: swaerosol_form = '            '
integer                      :: swaer_model_nintvls = 0
logical                      :: do_swaerosol = .false.
integer                      :: imax_aerfile = 0
integer                      :: jmax_aerfile = 0
integer                      :: kmax_aerfile = 0
integer                      :: ktop_aer = 0
integer                      :: kbot_aer = 0

!--------------------------------------------------------------------
! current data sets:

! shettle_fenn_mar_0 :  swaerosol_model_type = 'prescribed'
!                       swaerosol_form = 'shettle_fenn_mar_0'
!                       swaer_model_nintvls=40
!                       ktop_aer=36, kbot_aer=40
!                       imax= 1  , jmax=   1, kmax =   40
! shettle_fenn_mar_95 :  swaerosol_model_type = 'prescribed'
!                       swaerosol_form = 'shettle_fenn_mar_95'
!                       swaer_model_nintvls=40
!                       ktop_aer=36, kbot_aer=40
!                       imax=   1 , jmax=   1, kmax =   40
! robock              :  swaerosol_model_type = 'prescribed'
!                       swaerosol_form = 'robock'
!                       swaer_model_nintvls=242
!                       ktop_aer=17, kbot_aer=26
!                       imax=   1 , jmax=  60, kmax =   40
! ramachandran        :  swaerosol_model_type = 'prescribed'
!                       swaerosol_form = 'ramachandran'
!                       swaer_model_nintvls=26
!                       ktop_aer=1 , kbot_aer=40
!                       imax=   1 , jmax=  60, kmax =   40

!---------------------------------------------------------------------

namelist / esfsw_scattering_nml /    &
                                 do_swaerosol, &
				 swaerosol_form, &
				 swaer_model_type, &
				 swaer_model_nintvls, &
				 imax_aerfile, &
				 jmax_aerfile, &
				 kmax_aerfile, &
				 ktop_aer, &
				 kbot_aer




!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

!---------------------------------------------------------------------
! aerextband      = the parameterization band values of the extinction c
!                   coefficient for aerosols                           c
!                                                                      c
! aerssalbband    = the parameterization band values of the single-    c
!                   scattering albedo for aerosols                     c
!                                                                      c
! aerasymmband    = the parameterization band values of the asymmetry  c
!                   factor for aerosols                                c
!                                                                      c
! aerextrefstr    = reference extinction coefficient. used for scaling c
!                   
! aeramtsc    = scaling factor (non-dimensional) to vary the aerosol   c
!               amount relative to its reference value for a given     c
!               size distribution                                     
!---------------------------------------------------------------------

real, dimension(:,:,:,:), allocatable    :: aerextband, &
					    aerssalbband, &
					    aerasymmband
real, dimension(:,:,:),   allocatable    :: aeramtsc

!--------------------------------------------------------------------
integer, parameter                       :: NICECLDIVLS=25
integer, parameter                       :: NLIQCLDIVLS=24
integer, parameter                       :: NAERIVLS_SF=40
integer, parameter                       :: NAERMODELS=2

real, dimension(:,:), allocatable        :: solivlaero  
integer, dimension(:), allocatable       :: nivl1aero, nivl2aero
integer                                  :: NAERIVLS

integer, dimension(:), allocatable       :: endaerwvn, endaerwvnstr
real,    dimension(:), allocatable       :: solflxband

!--------------------------------------------------------------------
!    reference extinction coefficient is from 0.55 micron, which is 
!    band (xx) in the robock model
!--------------------------------------------------------------------
real     :: aerextrefstr = 1.58614E-02



!-------------------------------------------------------------------
!nteger  :: kmin, kmax
!integer  :: x(4), y(4)
      integer, dimension(:), allocatable :: jindx2



!---------------------------------------------------------------------
!---------------------------------------------------------------------





                   contains


!subroutine esfsw_scattering_init (kmin_in, kmax_in, latb)
subroutine esfsw_scattering_init (                  latb)
 
!integer, intent(in)   :: kmin_in, kmax_in
real, dimension(:), intent(in) :: latb
 

!----------------------------------------------------------------------c
!  define the spectral limits for aerosol single   
!  scattering properties.                                              c
!  note: the last wavenumber value must be the same as the             c
!        parameterization's last band limit.                           c
!-----------------------------------------------------------------------
 
!----------------------------------------------------------------------c
! wavenumber limits for shettle & fenn aerosol intervals               c
!----------------------------------------------------------------------c
 
     integer, dimension(NAERIVLS_SF)     :: endaerwvnsf

     data endaerwvnsf / 250,   333,   400,   469,   541,   581,   610, &
                       667,   676,   800,   870,   909,   944,  1000, &
                      1087,  1111,  1149,  1220,  1266,  1389,  1538, &
                      1613,  1667,  1818,  2000,  2222,  2667,  2948, &
                      3333,  3704,  4000,  4444,  5000,  6510,  9434,&
                     14409, 18182, 29674, 33333, 57600 /
 
 !----------------------------------------------------------------------
     integer                 :: nband,nw, nivl3
     real                    :: sumsol3
     integer                 :: unit, io, ierr
     integer                 :: ni, k, j, i
     character(len=64)       :: name

     real,    dimension(:), allocatable :: solarfluxtoa
     real,    dimension(:), allocatable :: data_lat
!     integer, dimension(:), allocatable :: jindx2
     integer, dimension(:), allocatable :: endwvnbands
     real :: lat_start
     integer :: j, jst, jj

!--------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=esfsw_scattering_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'esfsw_scattering_nml')
       enddo
10     call close_file (unit)
     endif

     unit = open_file ('logfile.out', action='append')
!    call print_version_number (unit, 'esfsw_scattering',    &
!						      version_number)
     if (get_my_pe() == 0)  then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=esfsw_scattering_nml)
     endif
     call close_file (unit)

!------------------------------------------------------------------
!  retrieve model dimensions
!---------------------------------------------------------------------

!    kmin = kmin_in
!    kmax = kmax_in
!     call get_domain_decomp (x, y)


     if (Environment%running_gcm) then
     if (trim (swaerosol_form) /= ' ' .and. jmax_aerfile /= 0) then
        allocate (jindx2  (size(latb,1)-1))
        call map_global_indices (jmax_aerfile, latb, jindx2)



!     if (jmax_aerfile /= 0) then
!      if (get_num_pes()*(size(latb,1) -1) == jmax_aerfile) then
!       lat_start = -90.0 + 180./(2.*float(jmax_aerfile))
!       allocate (data_lat(jmax_aerfile))
!      do j=1,jmax_aerfile
!        data_lat(j) = lat_start + real(j-1)*180.0/real(jmax_aerfile)
!      end do
!      jst = 1
!      do jj=1, size(latb,1) - 1
!        do j = jst,jmax_aerfile
!   if (data_lat(j) >= latb(jj)*radians_to_degrees ) then
!     jindx2(jj) = j
!     jst = j
!            exit
!          endif
!        end do
!      end do
!      else
!        call error_mesg ('esfsw_scattering_mod', &
!           'resolution of data input file doesn''t match model size --&
!    &  must provide inter(extra)polation not yet implemented', &
!           FATAL)
!      endif
     endif
     endif

!---------------------------------------------------------------------
!   define existence and type of shortwave aerosols in model in the
!   derived data type Sw_control
!----------------------------------------------------------------------

     Sw_control%do_swaerosol = do_swaerosol
     if (.not. do_swaerosol) swaerosol_form = '     '
     Sw_control%swaerosol_form = swaerosol_form

!----------------------------------------------------------------------c
!   allocate module variables
!-------------------------------------------------------------------
       allocate ( nivl1aero(nbands) )
       allocate ( nivl2aero(nbands) )
     if (trim(swaerosol_form) /= '     '  ) then
       if (imax_aerfile <= 0 .or. jmax_aerfile <= 0 .or.  &
	   kmax_aerfile <= 0) then
         call error_mesg ('esfsw_scattering_init',  &
            'must specify i,j,k extents of sw aerosol type', FATAL)
       endif
       if (Environment%running_gcm) then
         allocate (aerextband   (imax_aerfile, jmax_aerfile,     &
	               	         kmax_aerfile, nbands) )
         allocate (aerssalbband (imax_aerfile, jmax_aerfile,   &
	               	         kmax_aerfile, nbands) )
         allocate (aerasymmband (imax_aerfile, jmax_aerfile, & 
	               	         kmax_aerfile, nbands) )
         allocate (aeramtsc     (imax_aerfile, jmax_aerfile,   &
	               	         kmax_aerfile) )
       else if (Environment%running_standalone) then
         allocate (aerextband   (imax_aerfile, 1,     &
	               	         kmax_aerfile, nbands) )
         allocate (aerssalbband (imax_aerfile, 1,   &
	               	         kmax_aerfile, nbands) )
         allocate (aerasymmband (imax_aerfile, 1, & 
	               	         kmax_aerfile, nbands) )
         allocate (aeramtsc     (imax_aerfile, 1,   &
	               	         kmax_aerfile) )
       endif
     else
       imax_aerfile = 0
       jmax_aerfile = 0
       kmax_aerfile = 0
     endif

!--------------------------------------------------------------------
!    retrieve sw aerosol data. 
!--------------------------------------------------------------------

     if (do_swaerosol .and. trim(swaerosol_form) == 'predicted' ) then
!--------------------------------------------------------------------
!    at the moment, predicted shortwave aerosols are not available, 
!    and their request leads to model shutdown.
!-------------------------------------------------------------------
       call error_mesg ('esfsw_scattering_init',  &
                'predicted sw aerosols not currently implemented', &  
							      FATAL)

     else if (do_swaerosol .and. trim(swaerosol_form)  == 'observed' &
	            .or. trim(swaerosol_form) == 'prescribed' )   then

!---------------------------------------------------------------------
!    prescribed / observed aerosol distributions
!---------------------------------------------------------------------
       if ( trim(swaer_model_type) == 'shettle_fenn_mar_0')  then
!--------------------------------------------------------------------
!    shettle-fenn aerosol type (dry atmosphere). this is a global, 
!    tropospheric, maritime, 0 percent relative humidity aerosol 
!    profile defined on NAERIVLS_SF frequency intervals
!-------------------------------------------------------------------
         NAERIVLS = swaer_model_nintvls
         allocate (endaerwvn(NAERIVLS))
         endaerwvn(:) = endaerwvnsf(:)
 
       else if (trim(swaer_model_type) == 'shettle_fenn_mar_95') then
!------------------------------------------------------------------
!    shettle-fenn aerosol type (humid atmosphere). this is a global, 
!    tropospheric, maritime, 95 percent relative humidity aerosol 
!    profile defined on NAERIVLS_SF frequency intervals
!------------------------------------------------------------------
         NAERIVLS = swaer_model_nintvls
         allocate (endaerwvn(NAERIVLS))
         endaerwvn(:) = endaerwvnsf(:)
 
       else if (trim(swaer_model_type) == 'robock') then
!----------------------------------------------------------------------
!   robock type aerosols. this is a global, stratospheric aerosol 
!   profile on about 242 frequency intervals
!-------------------------------------------------------------------
         NAERIVLS = swaer_model_nintvls
         allocate (endaerwvn(NAERIVLS))
 
         name = 'INPUT/swstratendaerdata'
	 if (file_exist(name) ) then
           unit = open_file (name, form= 'unformatted',     &
                             action= 'read')
	   read (unit) endaerwvn
	   call close_file (unit)
	 else
	   call error_mesg ('esfsw_scattering_init',  &
	             'no robock swstratendaerdata file present', FATAL)
	 endif

       else if (trim(swaer_model_type)  == 'ramachandran') then
!------------------------------------------------------------------
!    ramachandran type aerosol. this is a zonal, stratospheric  
!    aerosol profile on about 26 frequency intervals
!------------------------------------------------------------------
         NAERIVLS = swaer_model_nintvls
         allocate (endaerwvn(NAERIVLS))
 
         name = 'INPUT/swstratendramadata'
	 if (file_exist(name) ) then
           unit = open_file (name, form= 'formatted',     &
                             action= 'read')
	   read (unit, FMT = '(20i6)' ) endaerwvn
	   call close_file (unit)
	 else
	   call error_mesg ('esfsw_scattering_init',  &
	           'no rama swstratendaerdata file present', FATAL)
	 endif
 
!-------------------------------------------------------------------
!  error condition -- prescribed type is not acceptable
!-------------------------------------------------------------------
       else
         call error_mesg( 'esfsw_scattering_init',  &
	   ' shortwave aerosol type is not an acceptable value.', FATAL)
       end if
 
!----------------------------------------------------------------------c
! define the solar weights and interval counters that are used to      c
! determine the single-scattering properties for the parameterization  c
! band spectral intervals, from the specified spectral intervals for   c
! aerosols.  
!----------------------------------------------------------------------c
       allocate ( solivlaero(nbands,NAERIVLS))
       allocate ( solflxband(nbands) )
       allocate ( solarfluxtoa(TOT_WVNUMS) )
       allocate ( endwvnbands(0:nbands) )
 
       call get_solarfluxes (solarfluxtoa, solflxband, endwvnbands)

       nivl3 = 1
       sumsol3 = 0.0
       nband = 1
       solivlaero(:,:) = 0.0
       nivl1aero(1) = 1
 
       do nw = 1,endwvnbands(nbands)
 
	 sumsol3 = sumsol3 + solarfluxtoa(nw)
 
         if ( nw.eq.endaerwvn(nivl3) ) then
           solivlaero(nband,nivl3) = sumsol3
           sumsol3 = 0.0
         end if
 
         if ( nw.eq.endwvnbands(nband) ) then
 
           if ( nw.ne.endaerwvn(nivl3) ) then
             solivlaero(nband,nivl3) = sumsol3 
             sumsol3 = 0.0
           end if
 
           nivl2aero(nband) = nivl3
 
           nband = nband + 1
 
           if ( nband.le.nbands ) then
 
             if ( nw.eq.endaerwvn(nivl3) ) then
               nivl1aero(nband) = nivl3 + 1
             else
               nivl1aero(nband) = nivl3
             end if
 
           end if
 
         end if
 
         if ( nw.eq.endaerwvn(nivl3) ) nivl3 = nivl3 + 1
 
       end do
 
       deallocate (solarfluxtoa)
       deallocate (endwvnbands)

!------------------------------------------------------------------
! obtain observed or prescribed aerosol data
!------------------------------------------------------------------
       call aeropar
     endif

!---------------------------------------------------------------------


end subroutine esfsw_scattering_init


!####################################################################

!subroutine map_global_indices (global_rows , latb,    &
!                              global_index_array )

!integer, intent(in)   :: global_rows 
!real, dimension(:), intent(in) :: latb
!integer, dimension(:), intent(out) :: global_index_array

!      real   :: lat_start
!      real, dimension(global_rows) :: data_lat
!      integer   :: j, jst, jj

!      if (get_num_pes()*(size(latb,1) -1) == global_rows) then
!       lat_start = -90.0 + 180./(2.*float(global_rows))
!      do j=1,global_rows
!        data_lat(j) = lat_start + real(j-1)*180.0/real(global_rows)
!      end do
!      jst = 1
!      do jj=1, size(latb,1) - 1
!        do j = jst,global_rows
!   if (data_lat(j) >= latb(jj)*radians_to_degrees ) then
!     global_index_array(jj) = j
!     jst = j
!            exit
!          endif
!        end do
!      end do
!      else
!        call error_mesg ('esfsw_scattering_mod', &
!           'resolution of data input file doesn''t match model size --&
!    &  must provide inter(extra)polation not yet implemented', &
!           FATAL)
!      endif

!end subroutine map_global_indices 



!####################################################################

subroutine scatpar 

!------------------------------------------------------------------- 
! this routine is provided as the mechanism ot obtain the parameteriz-
! ation band values of the single scattering parameters (extinction 
! coefficient, single-scattering albedo and asymmetry factor) for drops,
! ice particles and aerosols, from the drop size, ice particle size and 
! the aerosol type, respectively, WHEN PREDICTED SW AEROSOLS ARE USED.
! CURRENTLY IT IS NOT ACCESSIBLE.
!--------------------------------------------------------------------
                                                                        

!real, dimension(:,:,:,:), intent(out) :: aerext_out, aerssalb_out, &
!				  aerasymm_out
!real, dimension(:,:,:),   intent(out) :: aeramtsc_out
!-------------------------------------------------------------------


!     aerext_out  (:,:,:,:) = 0.
!     aerssalb_out(:,:,:,:) = 0.
!     aerasymm_out(:,:,:,:) = 0.
!     aeramtsc_out(:,:,:)   = 0.

 
end  subroutine scatpar




!#####################################################################

subroutine get_swaerprops_from_scattering ( js, ix, jx, n,   &
           aerextband_out, aerssalbband_out, aerasymmband_out, &
	   aeramtsc_out)
 
!---------------------------------------------------------------------
integer,                intent(in)  ::   js, ix, jx, n
real, dimension(:,:,:), intent(out) ::   aerextband_out,    &
					 aerssalbband_out,   &
					 aerasymmband_out,  &
                                         aeramtsc_out
!---------------------------------------------------------------------

     integer    :: i, j, jindx
     integer    :: israd, ierad, jsrad, jerad


     israd = 1
     jsrad = 1
     ierad = ix
     jerad = jx

!--------------------------------------------------------------------
!   initialize output fields
!--------------------------------------------------------------------
     aerextband_out   (:,:,:) = 0.0
     aerssalbband_out (:,:,:) = 0.0  
     aerasymmband_out (:,:,:) = 0.0  
     aeramtsc_out     (:,:,:) = 1.0 

!-------------------------------------------------------------------
!   no specification of swaerosols (default case, return with above
!   values)
!-------------------------------------------------------------------
     if (imax_aerfile == 0 .and. jmax_aerfile == 0 .and.   &
         kmax_aerfile == 0 ) then  
!-------------------------------------------------------------------
!   specification of swaerosols as function of k only 
!   assumption made that aerosols are specified at model levels
!   if not, interpolation required
!-------------------------------------------------------------------
     else if (imax_aerfile == 1 .and. jmax_aerfile == 1 ) then
       do j=JSRAD,JERAD
         do i=ISRAD, IERAD
           aerextband_out(i,j,ktop_aer:kbot_aer) =     &
               aerextband(1,1,ktop_aer:kbot_aer,n)

           aerssalbband_out(i,j,ktop_aer:kbot_aer) =     &
               aerssalbband(1,1,ktop_aer:kbot_aer,n)

           aerasymmband_out(i,j,ktop_aer:kbot_aer) =     &
               aerasymmband(1,1,ktop_aer:kbot_aer,n)

           aeramtsc_out (i,j,ktop_aer:kbot_aer) = &
	        aeramtsc(1,1,ktop_aer:kbot_aer)
         enddo
       enddo

!-------------------------------------------------------------------
!  specification of swaerosols as function of k and j only 
!  assumption made that aerosols are specified at model levels
!  if not, interpolation required
!  assume a 1 to 1 transformation in the latitude data, with
!  JSRAD data being the data at latitude jabs(j)+y(3)-1. a
!  transformation of data to model grid is needed if this
!  assumption is not correct.
!-------------------------------------------------------------------
     else if (imax_aerfile == 1 ) then
       do j=JSRAD,JERAD
         if (Environment%running_standalone) then
	   jindx = 1
         else if (Environment%running_gcm) then
!	   jindx = jabs(j) + y(3) - 1
!	   jindx = js+j-1  + y(3) - 1
           jindx = jindx2(js+j-1)
         endif
         do i=ISRAD, IERAD
           aerextband_out(i,j,ktop_aer:kbot_aer) =     &
               aerextband(1, jindx, ktop_aer:kbot_aer,n)

           aerssalbband_out(i,j,ktop_aer:kbot_aer) =     &
               aerssalbband(1, jindx, ktop_aer:kbot_aer,n)

           aerasymmband_out(i,j,ktop_aer:kbot_aer) =     &
               aerasymmband(1, jindx, ktop_aer:kbot_aer,n)

           aeramtsc_out (i,j,ktop_aer:kbot_aer) = &
	        aeramtsc(1,jindx,ktop_aer:kbot_aer)
         enddo
       enddo

!---------------------------------------------------------------
!  assume that the (i,j,k) data points in the"... band" arrays map to 
!  the model grid points.  otherwise a transformation of the data 
!  from its grid layout to the model grid must be made here before 
!  returning to the calling routine.
!---------------------------------------------------------------
     else 
       do j=JSRAD,JERAD
         do i=ISRAD, IERAD
           aerextband_out(i,j,ktop_aer:kbot_aer) =     &
               aerextband(i,j, ktop_aer:kbot_aer,n)

           aerssalbband_out(i,j,ktop_aer:kbot_aer) =     &
                aerssalbband(i,j, ktop_aer:kbot_aer,n)

           aerasymmband_out(i,j,ktop_aer:kbot_aer) =     &
               aerasymmband(i,j, ktop_aer:kbot_aer,n)

           aeramtsc_out (i,j,ktop_aer:kbot_aer) = &
	        aeramtsc(i,j,ktop_aer:kbot_aer)
         enddo
       enddo
     endif
!------------------------------------------------------------------


end subroutine get_swaerprops_from_scattering





!####################################################################

subroutine aeropar
 
!----------------------------------------------------------------------c
! determine the parameterization band values of the single scattering  c
! parameters (extinction coefficient, single-scattering albedo and     c
! asymmetry factor) for aerosols, from the aerosol type, respectively. c
!----------------------------------------------------------------------c

!----------------------------------------------------------------------
! local variables:                                                     c
!----------------------------------------------------------------------c

!--------------------------------------------------------------------
! aermodel    = model number for the specified single-scattering      
!               properties                                           
! aerextivl   = the specified spectral values of the extinction        c
!               coefficient for aerosols                               c
! aerssalbivl = the specified spectral values of the single-           c
!               scattering albedo ifor aerosols                        c
! aerasymmivl = the specified spectral values of the asymmetry         c
!               factor for aerosols                                    c
!--------------------------------------------------------------------
                                                                       
     real, dimension(:,:,:,:), allocatable  :: aerextivl,    &
					       aerssalbivl, &
	            			       aerasymmivl
     real, dimension(:),       allocatable  :: aerextivlstr,     &
				               aerssalbivlstr, &
                                               aerasymmivlstr
     integer, dimension(:,:,:), allocatable :: aermodel

     integer                 :: naerivl
     integer                 :: unit, io, ierr
     integer                 :: j, k, lwb
     character(len=64)       :: name
	
!---------------------------------------------------------------------- 
!   allocate local variables
!---------------------------------------------------------------------
     if (Environment%running_gcm) then
       allocate (aermodel     (imax_aerfile, jmax_aerfile,     &
			       kmax_aerfile))
       allocate (aerextivl    (imax_aerfile, jmax_aerfile,    &
			       kmax_aerfile, NAERIVLS) )
       allocate (aerssalbivl  (imax_aerfile, jmax_aerfile,    &
			       kmax_aerfile, NAERIVLS) )
       allocate (aerasymmivl  (imax_aerfile, jmax_aerfile,   &
			       kmax_aerfile, NAERIVLS) )
     else if (Environment%running_standalone) then
       allocate (aermodel     (imax_aerfile, 1,     &
			       kmax_aerfile))
       allocate (aerextivl    (imax_aerfile, 1,    &
			       kmax_aerfile, NAERIVLS) )
       allocate (aerssalbivl  (imax_aerfile, 1,    &
			       kmax_aerfile, NAERIVLS) )
       allocate (aerasymmivl  (imax_aerfile, 1,   &
			       kmax_aerfile, NAERIVLS) )
     endif

!--------------------------------------------------------------------
!  define model type for shettle & fenn distributions.
!--------------------------------------------------------------------
     if (trim(swaer_model_type) == 'shettle_fenn_mar_0')   then
       aermodel(:,:,:) = 1
     else if (trim(swaer_model_type) == 'shettle_fenn_mar_95')  then
       aermodel(:,:,:) = 2
     endif

!---------------------------------------------------------------------- 
! define the single scattering parameters for aerosols based on the     
! model type using the shettle & fenn distributions.                    
! references:                                                           
! shettle, e.p. and r.w. fenn, models for the aerosols of the lower     
!      atmosphere and the effects of humidity variations on their       
!      optical properties,afgl-tr-79-0214,1979,94pp.                    
!---------------------------------------------------------------------- 
     if (trim(swaer_model_type) == 'shettle_fenn_mar_0' .or. & 
         trim(swaer_model_type) == 'shettle_fenn_mar_95')  then
       call aerosf (aermodel, aerextivl, aerssalbivl, aerasymmivl)

!---------------------------------------------------------------------
!   define default aerosol scaling values for the shettle_fenn
!   maritime profile(s). the assumption is that the aerosol is
!   in layers 36-40 of the l40 SKYHI level structure.
!---------------------------------------------------------------------
       aeramtsc(:,:,1:35) = 0.0
       aeramtsc(:,:,36) = 0.68180196
       aeramtsc(:,:,37) = 0.90226578
       aeramtsc(:,:,38) = 1.35693120
       aeramtsc(:,:,39) = 2.59581635
       aeramtsc(:,:,40) = 6.57827625

!----------------------------------------------------------------------
!   define single-scattering parameters for robock type aerosols
!----------------------------------------------------------------------
     else if (trim(swaer_model_type)  == 'robock') then
       allocate (aerextivlstr(NAERIVLS),       &
                 aerssalbivlstr(NAERIVLS),       &
                 aerasymmivlstr(NAERIVLS))

       name = 'INPUT/swstrataerosoldata'
       if (file_exist(name) ) then
         unit = open_file (name, form= 'unformatted',     &
                           action= 'read')
         do naerivl = 1,NAERIVLS
           read (unit)                                    &
               aerextivlstr(naerivl),                   &
               aerssalbivlstr(naerivl),                 &
               aerasymmivlstr(naerivl)
         enddo
         call close_file (unit)
       else
         call error_mesg ('esfsw_scattering_init',  &
	           'no robock aerosol data file present', FATAL)
       endif
       do naerivl = 1,NAERIVLS
         aerextivl  (:,:,:,naerivl) = aerextivlstr(naerivl)
         aerssalbivl(:,:,:,naerivl) = aerssalbivlstr(naerivl)
	 aerasymmivl(:,:,:,naerivl) = aerasymmivlstr(naerivl)
       enddo
!----------------------------------------------------------------------
!   define default aerosol scaling values for the robock profile. the 
!   assumption is that the aerosol is in layers 17-26 of the l40 SKYHI 
!   level structure. NO NORMALIZATION TO A FIXED OPTICAL DEPTH IS DONE.
!   (to do so, use the 1.0E-1*0.1/aerextrefstr here, and in swresf
!    do not multiply by deltaz). (this assumes a tot aero opt dep of
!    1.0E-1).
!--------------------------------------------------------------------
       aeramtsc(:,:,1:16) =  0.0
       aeramtsc(:,:,17:26) = 0.1
       aeramtsc(:,:,27:40) = 0.0

       deallocate (aerextivlstr)
       deallocate (aerssalbivlstr)
       deallocate (aerasymmivlstr)

!----------------------------------------------------------------------
!    ranmachandran type aerosols
!----------------------------------------------------------------------
     else if (trim(swaer_model_type) == 'ramachandran') then

!--------------------------------------------------------------------
!     read in the time-independent shortwave aerosol extinction and
!     albedo coefficients for rama's code. these are given for the
!     frequency band structure appropriate to the present shortwave
!     radiation code, so there is no need to do a frequency inter-
!     polation, as in the shettle-fenn and robock cases.
!    the sw aerosol data are for (jmax_aerfile) latitudes.input data
!    are arranged as if there are (freq,alt,lat) arrays. the number of
!    frequency bands in the data is 1 less than the number of frequency
!    bands in the sw program (bands 2-nbands have data, band 1 data is
!    set to zero). the alt bands are SKYHI-model dependent, at present
!    40. the number of latitude bands is also SKYHI-model dependent, at
!    present 60 (ie, the current value of JD).
!
!    read and process the extinction data
!---------------------------------------------------------------------

       name = 'INPUT/swaerosolextdata'
       if (file_exist(name) ) then
         allocate (aerextivlstr(NAERIVLS) )
         unit = open_file (name, form= 'ieee32',     &
                           action= 'read')
!---------------------------------------------------------------------
!   use latitude 40 (~30N) as sample for 1-column code
!---------------------------------------------------------------------
         if (Environment%running_standalone) then
           do j = 1,jmax_aerfile
             do k = 1, kmax_aerfile
               read (unit)                                   &
                 (aerextivlstr(lwb), lwb=2,NAERIVLS)
	       if (j == 40) then
 	         do lwb = 2,NAERIVLS
 	           aerextivl(:,1,k,lwb) = aerextivlstr(lwb)
   	         enddo
		 aerextivl(:,1,k,1) = 0.0
	       endif
             end do
           end do

!---------------------------------------------------------------------
!   use all latitudes for other codes
!---------------------------------------------------------------------
         else if (Environment%running_gcm) then
           do j = 1,jmax_aerfile
             do k = 1, kmax_aerfile
               read (unit)                                  &
                   (aerextivlstr(lwb), lwb=2,NAERIVLS)
   	       do lwb = 2,NAERIVLS
 	         aerextivl(:,j,k,lwb) = aerextivlstr(lwb)
 	       enddo
	       aerextivl(:,j,k,1) = 0.0
             end do
           end do
         endif
	 call close_file (unit)
         deallocate (aerextivlstr)
       else
	 call error_mesg ('esfsw_scattering_init',  &
	         'no rama extinction aerosol data file present', FATAL)
       endif
 
!---------------------------------------------------------------------
!    read and process the ssalb data
!---------------------------------------------------------------------
       name = 'INPUT/swaerosolssalbdata'
       if (file_exist(name) ) then
         allocate ( aerssalbivlstr(NAERIVLS)  )
         unit = open_file (name, form= 'ieee32',     &
                           action= 'read')
!---------------------------------------------------------------------
!   use latitude 40 (~30N) as sample for 1-column code
!---------------------------------------------------------------------
         if (Environment%running_standalone) then
           do j = 1,jmax_aerfile
             do k = 1, kmax_aerfile
               read (unit)                                   &
                   (aerssalbivlstr(lwb), lwb=2,NAERIVLS)
	       if (j == 40) then
 	         do lwb = 2,NAERIVLS
 	           aerssalbivl(:,1,k,lwb) = aerssalbivlstr(lwb)
 	         enddo
		 aerssalbivl(:,1,k,1) = 0.0
	       endif
             end do
           end do

!---------------------------------------------------------------------
!   use all latitudes for other codes
!---------------------------------------------------------------------
         else if (Environment%running_gcm) then
           do j = 1,jmax_aerfile
             do k = 1, kmax_aerfile
               read (unit)                                  &
                  (aerssalbivlstr(lwb), lwb=2,NAERIVLS)
 	       do lwb = 2,NAERIVLS
 	         aerssalbivl(:,j,k,lwb) = aerssalbivlstr(lwb)
   	       enddo
 	       aerssalbivl(:,j,k,1) = 0.0
             end do
           end do
         endif
	 call close_file (unit)
         deallocate (aerssalbivlstr)
       else
         call error_mesg ('esfsw_scattering_init',  &
	           'no rama ssalb aerosol data file present', FATAL)
       endif

!---------------------------------------------------------------------
!    read and process the asymm data
!---------------------------------------------------------------------
       name = 'INPUT/swaerosolasymmdata'
       if (file_exist(name) ) then
         allocate ( aerasymmivlstr(NAERIVLS))
         unit = open_file (name, form= 'ieee32',     &
                           action= 'read')
!---------------------------------------------------------------------
!   use latitude 40 (~30N) as sample for 1-column code
!---------------------------------------------------------------------
         if (Environment%running_standalone) then
           do j = 1,jmax_aerfile
             do k = 1, kmax_aerfile
               read (unit)                                   &
                  (aerasymmivlstr(lwb), lwb=2,NAERIVLS)
	       if (j == 40) then
 	         do lwb = 2,NAERIVLS
 	           aerasymmivl(:,1,k,lwb) = aerasymmivlstr(lwb)
 	         enddo
 	         aerasymmivl(:,1,k,1) = 0.0
	       endif
             end do
           end do

!---------------------------------------------------------------------
!   use all latitudes for other codes
!---------------------------------------------------------------------
         else if (Environment%running_gcm) then
           do j = 1,jmax_aerfile
             do k = 1, kmax_aerfile
               read (unit)                                  &
                  (aerasymmivlstr(lwb), lwb=2,NAERIVLS)
 	       do lwb = 2,NAERIVLS
 	         aerasymmivl(:,j,k,lwb) = aerasymmivlstr(lwb)
 	       enddo
 	       aerasymmivl(:,j,k,1) = 0.0
             end do
           end do
         endif
	 call close_file (unit)
         deallocate (aerasymmivlstr)
       else
         call error_mesg ('esfsw_scattering_init',  &
	             'no rama asymm aerosol data file present', FATAL)
       endif

!---------------------------------------------------------------------
!   the default aerosol scaling value for the ramachandran profile(s)
!   is unity.
!---------------------------------------------------------------------
       aeramtsc = 1.0

     endif
 
!----------------------------------------------------------------------c
! use the thick-averaging technique to define the single-scattering    c
! properties of the parameterization band spectral intervals from the  c
! specified spectral intervals for aerosols. 
!                                                                      c
! references:                                                          c
!                                                                      c
! edwards,j.m. and a. slingo, studies with a flexible new radiation    c
!      code I: choosing a configuration for a large-scale model.,      c
!      q.j.r. meteorological society, 122, 689-719, 1996.              c
!                                                                      c
! note: a thin-averaging technique (subroutine thinavg in the same 
!       module) is also available.                                    
!----------------------------------------------------------------------c
 
     call thickavg (nivl1aero , nivl2aero   , NAERIVLS   ,   &
                     aerextivl , aerssalbivl , aerasymmivl,  &
                     solivlaero, solflxband,       &
                     aerextband , aerssalbband, aerasymmband)

!---------------------------------------------------------------------
!  NOTE: aerextband, aerssalband and aerasymmband are left on the y
!        grid on which they were input. any interpolation will be done
!        later when used. All pes retain the entire latitude row set
!        at this time.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   deallocate local variables
!----------------------------------------------------------------------
     deallocate  (aermodel)
     deallocate  (aerextivl)
     deallocate  (aerssalbivl)
     deallocate  (aerasymmivl)
 
 !--------------------------------------------------------------------


end subroutine aeropar




!###################################################################

subroutine aerosf (aermodel, aerextivl_out, aerssalbivl_out,  &
		   aerasymmivl_out)      
 
!----------------------------------------------------------------------c
! define the single scattering parameters for aerosols based on the    c
! model type using the shettle & fenn distributions.                   c
!                                                                      c
! references:                                                          c
!                                                                      c
! shettle, e.p. and r.w. fenn, models for the aerosols of the lower    c
!      atmosphere and the effects of humidity variations on their      c
!      optical properties,afgl-tr-79-0214,1979,94pp.                   c
!----------------------------------------------------------------------c
 
integer, dimension(:,:,:), intent(in) :: aermodel
real, dimension(:,:,:,:), intent(out) :: aerextivl_out,     &
					 aerssalbivl_out, &
				         aerasymmivl_out
 
!----------------------------------------------------------------------c
! local variables:
!----------------------------------------------------------------------c
      integer     :: na
      real, dimension(NAERIVLS_SF,NAERMODELS)  :: aeroasymmivl,    &
						  aeroextivl,  &
				       	          aerossalbivl
 
!----------------------------------------------------------------------c
! model 1: shettle & fenn: maritime model 0% relative humidity.        c
!----------------------------------------------------------------------c
 
      data ( aeroextivl(na,1), na=1,NAERIVLS_SF )     &
                   / 3.486E-3,5.609E-3,4.310E-3,4.698E-3,5.397E-3,  &
                     5.957E-3,5.529E-3,4.435E-3,3.724E-3,3.758E-3, &
                     4.735E-3,5.663E-3,6.278E-3,7.415E-3,9.918E-3, &
                     1.209E-2,1.225E-2,1.028E-2,8.271E-3,8.761E-3, &
                     1.077E-2,1.460E-2,1.448E-2,1.257E-2,1.512E-2,  &
                     1.838E-2,2.143E-2,2.436E-2,2.973E-2,3.084E-2,  &
                     2.927E-2,3.205E-2,3.511E-2,4.087E-2,5.056E-2,  &
                     6.318E-2,7.503E-2,9.018E-2,1.032E-1,1.087E-1 /
      data ( aerossalbivl(na,1), na=1,NAERIVLS_SF )  &
                   / 0.0484,0.1534,0.2619,0.3552,0.4398,0.4941,0.4954, &
                     0.5350,0.6132,0.7150,0.8147,0.8428,0.8437,0.8428,&
                     0.7983,0.7684,0.7758,0.7659,0.7981,0.8602,0.9105, &
                     0.8987,0.8777,0.9236,0.9633,0.9758,0.9823,0.9812, &
                     0.9505,0.9143,0.9308,0.9644,0.9775,0.9771,0.9734, &
                     0.9791,0.9826,0.9796,0.9722,0.8786 /
      data ( aeroasymmivl(na,1), na=1,NAERIVLS_SF )   &
                   / 0.1172,0.3144,0.4158,0.4504,0.4728,0.4888,0.5151, &
                     0.5556,0.5800,0.5962,0.6033,0.6001,0.6007,0.5947, &
                     0.5831,0.5702,0.5747,0.6223,0.6645,0.6681,0.6611, &
                     0.6326,0.6491,0.6873,0.6759,0.6629,0.6691,0.6781, &
                     0.6533,0.6827,0.7272,0.7143,0.7082,0.6995,0.6890, &
                     0.6806,0.6761,0.6838,0.6940,0.7238 /
 
!----------------------------------------------------------------------c
! model 2: shettle & fenn: maritime model 95% relative humidity.       c
!----------------------------------------------------------------------c
 
      data ( aeroextivl(na,2), na=1,NAERIVLS_SF )    &
                   / 4.626E-2,1.018E-1,1.226E-1,1.443E-1,1.645E-1, &
                     1.772E-1,1.801E-1,1.788E-1,1.764E-1,1.572E-1, &
                     1.217E-1,1.006E-1,9.670E-2,1.047E-1,1.259E-1,&
                     1.425E-1,1.494E-1,1.582E-1,1.674E-1,1.834E-1, &
                     2.104E-1,2.370E-1,2.302E-1,2.210E-1,2.454E-1, &
                     2.706E-1,3.011E-1,3.341E-1,3.432E-1,3.038E-1, &
                     2.958E-1,3.355E-1,3.537E-1,3.695E-1,3.876E-1, &
                     4.083E-1,4.287E-1,4.596E-1,4.869E-1,5.047E-1 /
      data ( aerossalbivl(na,2), na=1,NAERIVLS_SF )   &
                   / 0.1121,0.2646,0.3080,0.3096,0.3033,0.2922,0.2814, &
                     0.2697,0.2618,0.2467,0.2481,0.3079,0.4126,0.5364, &
                     0.6388,0.6862,0.7034,0.7244,0.7473,0.7636,0.7578,&
                     0.6747,0.5564,0.7006,0.8895,0.8826,0.9182,0.8903, &
                     0.6345,0.6133,0.8741,0.9795,0.9839,0.9857,0.9950, &
                     0.9966,0.9970,0.9961,0.9943,0.9643 /
      data ( aeroasymmivl(na,2), na=1,NAERIVLS_SF )   &
                   / 0.2479,0.5266,0.5698,0.5947,0.6165,0.6336,0.6480, &
                     0.6684,0.6846,0.7250,0.7887,0.8213,0.8307,0.8295, &
                     0.8213,0.8151,0.8133,0.8138,0.8142,0.8115,0.8069, &
                     0.8069,0.8345,0.8392,0.8121,0.8056,0.7896,0.7737, &
                     0.8181,0.8847,0.8756,0.8318,0.8155,0.8055,0.7960, &
                     0.7941,0.7970,0.8010,0.8034,0.8057 /
 
      integer      :: j,k,i, nmodel, ni, jst, jend

!---------------------------------------------------------------------
!   define j bounds of arrays - will differ between model and standalone
!---------------------------------------------------------------------
      jst  = lbound(aerextivl_out,2)
      jend = ubound(aerextivl_out,2)

!-------------------------------------------------------------------
!   fill arrays with data from above
!-------------------------------------------------------------------
      do ni = 1,NAERIVLS_SF
        do k = 1,kmax_aerfile
	  do j=jst,jend
	    do i=1,imax_aerfile
              nmodel = aermodel(i,j,k)
              if (nmodel == 0 ) then
                aerextivl_out(i,j,k,ni) = 0.0
                aerssalbivl_out(i,j,k,ni) = 0.0
                aerasymmivl_out(i,j,k,ni) = 0.0
              else
                aerextivl_out(i,j,k,ni) = aeroextivl(ni,nmodel) 
                aerssalbivl_out(i,j,k,ni) = aerossalbivl(ni,nmodel) 
                aerasymmivl_out(i,j,k,ni) = aeroasymmivl(ni,nmodel) 
              end if
            end do
          end do
        end do
      end do
!-------------------------------------------------------------------
 

end subroutine aerosf 


!####################################################################



               end module esfsw_scattering_mod
