
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
				 aerosol_type,  &
			 aerosol_properties_type, Aerosol_props,&
				 map_global_indices
use microphys_rad_mod,     only: thickavg, thinavg
!use aerosol_mod,           only:                                &
!                                  aerosol_alloc

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
    character(len=128)  :: version =  '$Id: esfsw_scattering.F90,v 1.4 2003/04/09 20:59:23 fms Exp $'
    character(len=128)  :: tag     =  '$Name: inchon $'



!---------------------------------------------------------------------
!-------  interfaces --------

public           &
       esfsw_scattering_init, scatpar, get_swaerprops_from_scattering
     

private          &
       aeropar


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


!--------------------------------------------------------------------
integer, parameter                       :: NICECLDIVLS=25
integer, parameter                       :: NLIQCLDIVLS=24
integer, parameter                       :: NAERIVLS_SF=40
integer                                  :: NAERMODELS

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


!---------------------------------------------------------------------
!---------------------------------------------------------------------





                   contains


subroutine esfsw_scattering_init (                  latb)
 
real, dimension(:), intent(in) :: latb
 

!----------------------------------------------------------------------c
!  define the spectral limits for aerosol single   
!  scattering properties.                                              c
!  note: the last wavenumber value must be the same as the             c
!        parameterization's last band limit.                           c
!-----------------------------------------------------------------------
 
 
 !----------------------------------------------------------------------
     integer                 :: nband,nw, nivl3
     real                    :: sumsol3
     integer                 :: unit, io, ierr, unit2
     integer                 :: ni, k, j, i
     character(len=64)       :: name, name2

     real,    dimension(:), allocatable :: solarfluxtoa
     real,    dimension(:), allocatable :: data_lat
     integer, dimension(:), allocatable :: endwvnbands
     integer, dimension(:), allocatable :: endwvnb, endaerwvn8
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

         naermodels = Aerosol_props%naermodels
       allocate (Aerosol_props%aerextband   (nbands, NAERMODELS) )
        allocate (Aerosol_props%aerssalbband (nbands, NAERMODELS) )
        allocate (Aerosol_props%aerasymmband (nbands, NAERMODELS) )

       Aerosol_props%aerextband   = 0.                      
       Aerosol_props%aerssalbband = 0.                    
       Aerosol_props%aerasymmband = 0.                    
     else
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
           naerivls = Aerosol_props%num_wavenumbers
          allocate( endaerwvn(NAERIVLS) )
           endaerwvn = Aerosol_props%endaerwvnsf

       else if (trim(swaer_model_type) == 'shettle_fenn_mar_95') then
!------------------------------------------------------------------
!    shettle-fenn aerosol type (humid atmosphere). this is a global, 
!    tropospheric, maritime, 95 percent relative humidity aerosol 
!    profile defined on NAERIVLS_SF frequency intervals
!------------------------------------------------------------------
           naerivls = Aerosol_props%num_wavenumbers
        allocate( endaerwvn(NAERIVLS) )
           endaerwvn = Aerosol_props%endaerwvnsf
 
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
         call error_mesg ('esfsw_scattering_mod', &
           ' new ramachandran input file is needed -- code now &
	   &separates aersol amount and aerosol properties', FATAL)
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
           aerextband_out, aerssalbband_out, aerasymmband_out)
 
!---------------------------------------------------------------------
integer,                intent(in)  ::   js, ix, jx, n
real, dimension(:), intent(out) ::   aerextband_out,    &
					 aerssalbband_out,   &
					 aerasymmband_out
!---------------------------------------------------------------------


!--------------------------------------------------------------------
!   ... Assuming that optical properties are a function only of
!       aerosol type and wavenumber (not spatially varying)
!--------------------------------------------------------------------

    if (trim(swaerosol_form) /= '  ' ) then
      aerextband_out(:)   = Aerosol_props%aerextband(n,:)
     aerssalbband_out(:) = Aerosol_props%aerssalbband(n,:)
     aerasymmband_out(:) = Aerosol_props%aerasymmband(n,:)
    else
      aerextband_out(:)   = 0.                 
     aerssalbband_out(:) = 0.                  
     aerasymmband_out(:) = 0.                    
     endif



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
! aerextivl   = the specified spectral values of the extinction        c
!               coefficient for aerosols                               c
! aerssalbivl = the specified spectral values of the single-           c
!               scattering albedo ifor aerosols                        c
! aerasymmivl = the specified spectral values of the asymmetry         c
!               factor for aerosols                                    c
!--------------------------------------------------------------------
                                                                       
     real, dimension(:,:), allocatable  :: aerextivl,    &
                                           aerssalbivl, &
                                           aerasymmivl
  
      real, dimension(:), allocatable  :: aerextivl_in,    &
                                          aerssalbivl_in, &
                                          aerasymmivl_in, &
                                          aerextband_out1,&
                                          aerssalbband_out1,&
                                          aerasymmband_out1

     real, dimension(:),       allocatable  :: aerextivlstr,     &
				               aerssalbivlstr, &
                                               aerasymmivlstr

     integer                 :: naerivl
     integer                 :: unit, io, ierr, unit2
     integer                 :: j, k, lwb
     integer                 :: nmodel
     character(len=64)       :: name, name2
	
!---------------------------------------------------------------------- 
!   allocate local variables
!---------------------------------------------------------------------
     allocate (aerextivl    (NAERIVLS, NAERMODELS) )
     allocate (aerssalbivl  (NAERIVLS, NAERMODELS) )
     allocate (aerasymmivl  (NAERIVLS, NAERMODELS) )
     allocate (aerextivl_in  (NAERIVLS) )
     allocate (aerssalbivl_in (NAERIVLS) )
     allocate (aerasymmivl_in  (NAERIVLS) )
     allocate (aerextband_out1  (nbands) )
     allocate (aerssalbband_out1 (nbands) )
     allocate (aerasymmband_out1  (nbands) )



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
       aerextivl = Aerosol_props%aeroextivl
       aerssalbivl = Aerosol_props%aerossalbivl
       aerasymmivl = Aerosol_props%aeroasymmivl


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
         aerextivl  (naerivl,1) = aerextivlstr(naerivl)
         aerssalbivl(naerivl,1) = aerssalbivlstr(naerivl)
         aerasymmivl(naerivl,1) = aerasymmivlstr(naerivl)
       enddo
!----------------------------------------------------------------------
!   define default aerosol scaling values for the robock profile. the 
!   assumption is that the aerosol is in layers 17-26 of the l40 SKYHI 
!   level structure. NO NORMALIZATION TO A FIXED OPTICAL DEPTH IS DONE.
!   (to do so, use the 1.0E-1*0.1/aerextrefstr here, and in swresf
!    do not multiply by deltaz). (this assumes a tot aero opt dep of
!    1.0E-1).
!--------------------------------------------------------------------

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
!	           aerextivl(:,1,k,lwb) = aerextivlstr(lwb)
                   aerextivl(lwb,1) = aerextivlstr(lwb)
   	         enddo
!	 aerextivl(:,1,k,1) = 0.0
		 aerextivl(1,1) = 0.0
	       endif
             end do
           end do

!---------------------------------------------------------------------
!   use all latitudes for other codes
!---------------------------------------------------------------------
         else if (Environment%running_gcm) then
!          do j = 1,jmax_aerfile
             do k = 1, kmax_aerfile
               read (unit)                                  &
                   (aerextivlstr(lwb), lwb=2,NAERIVLS)
   	       do lwb = 2,NAERIVLS
!	         aerextivl(:,j,k,lwb) = aerextivlstr(lwb)
		 aerextivl(lwb,1) = aerextivlstr(lwb)
 	       enddo
!       aerextivl(:,j,k,1) = 0.0
	       aerextivl(1,1) = 0.0
             end do
!          end do
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
!	           aerssalbivl(:,1,k,lwb) = aerssalbivlstr(lwb)
 	           aerssalbivl(lwb,1) = aerssalbivlstr(lwb)
 	         enddo
!	 aerssalbivl(:,1,k,1) = 0.0
		 aerssalbivl(1,1) = 0.0
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
!	         aerssalbivl(:,j,k,lwb) = aerssalbivlstr(lwb)
 	           aerssalbivl(lwb,1) = aerssalbivlstr(lwb)
   	       enddo
!	       aerssalbivl(:,j,k,1) = 0.0
 	       aerssalbivl(1,1) = 0.0
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
!	           aerasymmivl(:,1,k,lwb) = aerasymmivlstr(lwb)
 	           aerasymmivl(lwb,1) = aerasymmivlstr(lwb)
 	         enddo
!	         aerasymmivl(:,1,k,1) = 0.0
 	         aerasymmivl(1,1) = 0.0
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
!	         aerasymmivl(:,j,k,lwb) = aerasymmivlstr(lwb)
 	           aerasymmivl(lwb,1) = aerasymmivlstr(lwb)
 	       enddo
!	       aerasymmivl(:,j,k,1) = 0.0
 	       aerasymmivl(1,1) = 0.0
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
!      aeramtsc = 1.0

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
 
   do nmodel = 1,NAERMODELS

       aerextivl_in(:) = aerextivl(:,nmodel)
       aerssalbivl_in(:) = aerssalbivl(:,nmodel)
       aerasymmivl_in(:) = aerasymmivl(:,nmodel)
!
       call thickavg (nivl1aero , nivl2aero   , NAERIVLS   ,   &
                      aerextivl_in , aerssalbivl_in , aerasymmivl_in,  &
                      solivlaero, solflxband,       &
                aerextband_out1 , aerssalbband_out1, aerasymmband_out1)
         Aerosol_props%aerextband(:,nmodel) = aerextband_out1(:)
       Aerosol_props%aerssalbband(:,nmodel) = aerssalbband_out1(:)
      Aerosol_props%aerasymmband(:,nmodel) = aerasymmband_out1(:)
   end do


!---------------------------------------------------------------------
!  NOTE: aerextband, aerssalband and aerasymmband are left on the y
!        grid on which they were input. any interpolation will be done
!        later when used. All pes retain the entire latitude row set
!        at this time.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   deallocate local variables
!----------------------------------------------------------------------
     deallocate  (aerextivl)
     deallocate  (aerssalbivl)
     deallocate  (aerasymmivl)
     deallocate (aerextivl_in  )
     deallocate (aerssalbivl_in  )
     deallocate (aerasymmivl_in   )
     deallocate (aerextband_out1  )
     deallocate (aerssalbband_out1  )
     deallocate (aerasymmband_out1   )
 
 !--------------------------------------------------------------------


end subroutine aeropar




!####################################################################



               end module esfsw_scattering_mod
