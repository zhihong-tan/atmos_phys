
                    module rad_step_setup_mod


use rad_output_file_mod,    only:  hold_sfctmp
use rad_utilities_mod,      only:  Rad_control, radiation_control_type,&
				   Environment, environment_type, &
				   cloudrad_control_type, Cldrad_control
use  utilities_mod,         only:  open_file, file_exist,   & 
                                   check_nml_error, error_mesg, &
                                   print_version_number, FATAL, NOTE, &
				   WARNING, get_my_pe, close_file, &
				   get_num_pes, get_domain_decomp
use time_manager_mod,       only:  time_type
use constants_new_mod,      only:  wtmair, grav, rgas

!-------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!       module to retrieve model variables needed within the modules
!       that are activated on a radiation time step
!
!--------------------------------------------------------------------
  



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!    character(len=5), parameter  ::  version_number = 'v0.09'
     character(len=128)  :: version =  '$Id: rad_step_setup.F90,v 1.2 2001/07/05 17:32:45 fms Exp $'
     character(len=128)  :: tag     =  '$Name: galway $'



!---------------------------------------------------------------------
!-------  interfaces --------

public     &
        rad_step_setup_init, &
        rad_step_setup_dr,  &
	rad_step_setup_dealloc


private    &
	rad_layers


!--------------------------------------------------------------------
!----    namelist -----

logical         :: do_bounds_chk = .false.
logical         :: all_column_radiation = .true.
logical         :: all_level_radiation = .true.
integer         :: topmost_radiation_level=-99


namelist / rad_step_setup_nml /         &
                                          do_bounds_chk,    &
					  all_level_radiation, &
					  topmost_radiation_level, &
					  all_column_radiation

!--------------------------------------------------------------------
!-------- public data  -----

!---------------------------------------------------------------------
!     temp    = temperature at data levels of model where radiation is
!               to be calculated.
!     press   =  pressure at data levels of model where radiation is to
!                be calculated.
!     rh2o    =  mass mixing ratio of h2o at model data levels where 
!                radiation is to be calculated.
!     jabs    = global j indices for each local j value.
!     iabs    = global i indices for each local i value.
!---------------------------------------------------------------------

real, dimension(:,:,:), allocatable, public    :: temp, press, rh2o, &
						  pflux, tflux, deltaz,&
						  cloud_ice, &
						  cloud_water
integer, dimension(:), allocatable,  public    :: jabs, iabs
integer,                             public    :: IMINP, IMAXP,   &
						  JMINP, JMAXP
integer,                             public    :: ISRAD, IERAD, &
						  JSRAD, JERAD, &
						  KSRAD, KERAD
real, dimension(:,:), allocatable,   public    :: sealp, snow, land
type(time_type),                     public    :: Rad_time_sv
real, dimension(:,:), allocatable,   public    :: lat_sv, albedo_sv

!--------------------------------------------------------------------
!------ private data ------

integer                ::   x(4), y(4)


!-------------------------------------------------------------------
!-------------------------------------------------------------------




	   contains






!####################################################################

subroutine rad_step_setup_init (kmax_in) 

integer, intent(in)  :: kmax_in

!------------------------------------------------------------------
      integer  :: unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=rad_step_setup_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'rad_step_setup_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'rad_step_setup', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=rad_step_setup_nml)
      endif
      call close_file (unit)

      if (all_level_radiation) then
	KSRAD = 1
	KERAD = kmax_in
      else
	KSRAD = topmost_radiation_level
	if (KSRAD <= 0) then
	  call error_mesg ('rad_step_setup_init', &
	    ' topmost_radiation_level MUST be specified > 0', FATAL)
	endif
	KERAD = kmax_in
      endif

!---------------------------------------------------------------------
!   define global and subdomain dimensions.
!--------------------------------------------------------------------
      call get_domain_decomp (x, y)

!--------------------------------------------------------------------

end subroutine rad_step_setup_init

!####################################################################

subroutine rad_step_setup_dr (is, ie, js, je, ttt_l, rrr_l, ts_l,  &
			      jabsp_l, iabsp_l, press3d_l, snow_l,  &
			      sealp_l, Rad_time, lat, albedo, land_in, &
	  		      cloud_ice_in, cloud_water_in)

!--------------------------------------------------------------------
integer,                   intent(in)           :: is, ie, js, je
real,    dimension(:,:,:), intent(in)           :: ttt_l, rrr_l,   &
						   press3d_l
real,    dimension(:,:),   intent(in)           :: ts_l
integer, dimension(:),     intent(in)           :: iabsp_l, jabsp_l
real,    dimension(:,:),   intent(in), optional :: snow_l, sealp_l
type(time_type),           intent(in), optional :: Rad_time
real,    dimension(:,:),   intent(in), optional :: lat, albedo, land_in
real,    dimension(:,:,:), intent(in), optional :: cloud_ice_in, &
                                                   cloud_water_in
!---------------------------------------------------------------------

      integer :: k, ip, i, j, iomsgs, np

!---------------------------------------------------------------------
!   define local indices  and dimensions
!--------------------------------------------------------------------
      IMINP = 1
      IMAXP = ie - is + 1
      JMINP = 1
      JMAXP = je - js + 1

!---------------------------------------------------------------------
!   define (i,j) extents over which radiation is to be calculated.
!---------------------------------------------------------------------

      if (all_column_radiation) then
	ISRAD = IMINP
	IERAD = IMAXP
	JSRAD = JMINP
	JERAD = JMAXP
      else   !  if other choice is desired, specify it here
        call error_mesg ('rad_step_setup_dr', &
               ' set of (i,j) points to do radiation not specified', &
							    FATAL)
      endif

      if (Environment%running_gcm) then

!---------------------------------------------------------------------
!   save time, latitude and surface albedo input fields
!---------------------------------------------------------------------
        if (Environment%running_fms) then
          Rad_time_sv = Rad_time
          allocate ( lat_sv (lbound(lat,1):ubound(lat,1),  &
	  		     lbound(lat,2):ubound(lat,2) ) )
          allocate (albedo_sv (lbound(albedo,1):ubound(albedo,1),  &
		               lbound(albedo,2):ubound(albedo,2) ) )
          allocate (land      (lbound(land_in,1):ubound(land_in,1),  &
		               lbound(land_in,2):ubound(land_in,2) ) )
          lat_sv = lat
          albedo_sv = albedo
	  land= land_in
        endif
      else if (Environment%running_standalone  .and.     &
               Cldrad_control%do_pred_cld_microphys ) then

!-------------------------------------------------------------------
!   save land, cloud ice, cloud water fields for use in
!   obtaining cloud microphysics parameters in standalone version
!-------------------------------------------------------------------
        allocate (land      (lbound(land_in,1):ubound(land_in,1),  &
                             lbound(land_in,2):ubound(land_in,2) ) )
        allocate (cloud_ice                                        &
                 (lbound(cloud_ice_in,1):ubound(cloud_ice_in,1), &
                  lbound(cloud_ice_in,2):ubound(cloud_ice_in,2), &
                  lbound(cloud_ice_in,3):ubound(cloud_ice_in,3)  ) )
        allocate (cloud_water                                      &
                 (lbound(cloud_water_in,1):ubound(cloud_water_in,1), &
                  lbound(cloud_water_in,2):ubound(cloud_water_in,2), &
                  lbound(cloud_water_in,3):ubound(cloud_water_in,3)  ) )

        land = land_in
        cloud_ice = cloud_ice_in
        cloud_water = cloud_water_in
      endif

!--------------------------------------------------------------------
!   allocate needed arrays
!--------------------------------------------------------------------
      allocate (temp (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )
      allocate (rh2o (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD  ) )
      allocate (press(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1) )
      allocate (jabs (JSRAD:JERAD) )
      allocate (iabs (ISRAD:IERAD) )
      if (Environment%running_gcm) then
        if (Environment%running_skyhi ) then
          allocate (sealp(ISRAD:IERAD, JSRAD:JERAD) )
          allocate (snow (ISRAD:IERAD, JSRAD:JERAD) )
        endif
      endif

!--------------------------------------------------------------------
!   define global index arrays for only the computational points
!--------------------------------------------------------------------
      jabs(JSRAD:JERAD) = jabsp_l(:)
      iabs(ISRAD:IERAD) = iabsp_l(:)

!--------------------------------------------------------------------
!   define the temperature field (kelvin), including the surface, the
!   mixing ratio field, and the pressure field, including the surface,
!   the snow cover field and the land-sea mask (if present).
!--------------------------------------------------------------------
      temp(:,:,KSRAD:KERAD)    = ttt_l(:,:,KSRAD:KERAD)
      temp (:,:,KERAD+1)       = ts_l(:,:)
      rh2o(:,:,KSRAD:KERAD)    = rrr_l(:,:,KSRAD:KERAD)
      press(:,:,KSRAD:KERAD+1) = press3d_l(:,:,KSRAD:KERAD+1)
      if (Environment%running_gcm) then
        if (Environment%running_skyhi)  then
          snow(:,:) = snow_l(:,:)
          sealp(:,:) = sealp_l(:,:)
        endif
      endif

!-----------------------------------------------------------------------
!     create a record of those values of mixing ratio and temperature 
!     which were out of bounds of the radiation package tables.
!     this code will run as is on both t90 and t3e.
!-----------------------------------------------------------------------
      if (do_bounds_chk) then
        if (Environment%running_skyhi) then
	  np = get_num_pes()
          do ip=0,np-1
	    if (get_my_pe()  == ip) then
	      iomsgs = open_file ('msg_file', form='formatted',     &
		  		   action='append')
              do k=KSRAD,KERAD
                do j=JSRAD,JERAD
                  do i=ISRAD,IERAD
                    if (rh2o(i,j,k) .LT. 2.0E-07) then
                      write (iomsgs,9010)  rh2o(i,j,k),    &
					   x(3)+iabs(i)-1, &
					   y(3)+jabs(j)-1, k
                    endif
                    if (temp(i,j,k) .LT. 100.0E+00 .OR.   &
                      temp(i,j,k) .GT. 370.0E+00) then
                      write (iomsgs,9020)  temp(i,j,k), &
					   x(3)+iabs(i)-1,  &
					   y(3)+jabs(j)-1, k
                    endif
                  end do
                end do
              end do
	      call close_file (iomsgs)
            else
	      call barrier()
            endif
          end do
        else if (Environment%running_fms) then
	  if (get_num_pes() > 1) then
            iomsgs = open_file ('msg_file', form='formatted',     &
			       threading='multi', action='append')
	  else
            iomsgs = open_file ('msg_file', form='formatted',     &
			       action='append')
	  endif
          do k=KSRAD,KERAD
            do j=JSRAD,JERAD
              do i=ISRAD,IERAD
                if (rh2o(i,j,k) .LT. 2.0E-07) then
                  write (iomsgs,9010)  rh2o(i,j,k),    &
				       x(3)+iabs(i)-1,  &
				       y(3)+jabs(j)-1, k
                endif
                if (temp(i,j,k) .LT. 100.0E+00 .OR.   &
                    temp(i,j,k) .GT. 370.0E+00) then
                  write (iomsgs,9020)  temp(i,j,k),    &
				       x(3)+iabs(i)-1,  &
				       y(3)+jabs(j)-1, k
                endif
              end do
            end do
          end do
          call close_file (iomsgs)
        endif
      endif

!-----------------------------------------------------------------------
!     be sure that the mixing ratio is no less than 2.0E-07 in Skyhi  
!     (3.0E-06 in FMS) and that the temperature lies between 100K and 
!     370K, the limits of the tables referenced within the radiation 
!     package.
!-----------------------------------------------------------------------
      if (Environment%using_fms_periphs) then
        rh2o(:,:,KSRAD:KERAD) = MAX(rh2o(:,:,KSRAD:KERAD), 3.0E-06)
      else if (Environment%using_sky_periphs) then
        rh2o(:,:,KSRAD:KERAD) = MAX(rh2o(:,:,KSRAD:KERAD), 2.0E-07)
      endif

      temp(:,:,KSRAD:KERAD) = MAX(temp(:,:,KSRAD:KERAD), 100.0E+00)
      temp(:,:,KSRAD:KERAD) = MIN(temp(:,:,KSRAD:KERAD), 370.0E+00)

!---------------------------------------------------------------------
!     define variable to indicate if radiation diagnostics are to be
!     done at any points on the latitude rows in this chunk
!--------------------------------------------------------------------
      Rad_control%do_diagnostics = .false.
      do j=JMINP,JMAXP
        if (Rad_control%do_raddg(jabs(j)) ) then
          Rad_control%do_diagnostics = .true.
        endif
      end do

!--------------------------------------------------------------------
!     compute pressure, temperature, altitude arrays for layer
!     boundaries that are used in radiation routines and in interface
!     routines between radiation and other physics (clouds, aerosols)
!--------------------------------------------------------------------
      call rad_layers

!--------------------------------------------------------------------
!     pass the surface temp field used by the radiation code to 
!     the archive package
!--------------------------------------------------------------------
      call hold_sfctmp (temp, rh2o, press)

!------------------------------------------------------------------

9010 format ( 'value of mixing ratio too small for radiation tables:',&
               / , e12.5,' at i= ', i4, ' , j= ', i4, ' and k= ', i4)
9020 format ( 'value of temperature out of range of radiation tables:',&
               / , e12.5,' at i= ', i4, ' , j= ', i4, ' and k= ', i4)
!------------------------------------------------------------------

 
end subroutine rad_step_setup_dr


!#####################################################################

subroutine rad_step_setup_dealloc

!-------------------------------------------------------------------
!  deallocate arrays
!-------------------------------------------------------------------

      deallocate (pflux)
      deallocate (tflux)
      deallocate (deltaz)
      if (Environment%running_gcm) then
        if (Environment%running_skyhi) then
          deallocate (snow)
          deallocate (sealp)
        endif
      endif
      deallocate (iabs)
      deallocate (jabs)
      deallocate (press )
      deallocate (rh2o )
      deallocate (temp )

      if (Environment%running_gcm) then
        if (Environment%running_fms) then
          deallocate ( lat_sv)
          deallocate ( albedo_sv)
	  deallocate (land)
        endif
      endif


end subroutine rad_step_setup_dealloc


!####################################################################

subroutine rad_layers


!---------------------------------------------------------------------
   integer       :: k,  j

!---------------------------------------------------------------------
   allocate ( pflux  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1 ) )
   allocate ( tflux  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1 ) )
   allocate ( deltaz (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD   ) )

!--------------------------------------------------------------------
!     define flux level pressures (pflux) as midway between data level
!     (layer-mean) pressures. specify temperatures at flux levels
!     (tflux).
!--------------------------------------------------------------------
   if (KSRAD == 1) then
     do k=KSRAD+1,KERAD
       pflux(:,:,k) = 0.5E+00*(press(:,:,k-1) + press(:,:,k))
       tflux(:,:,k) = 0.5E+00*(temp (:,:,k-1) + temp (:,:,k))
     enddo
     pflux(:,:,KSRAD  ) = 0.0E+00
     pflux(:,:,KERAD+1) = press(:,:,KERAD+1)
     tflux(:,:,KSRAD  ) = temp (:,:,KSRAD)
     tflux(:,:,KERAD+1) = temp (:,:,KERAD+1)
   else
     do k=KSRAD,KERAD
       pflux(:,:,k) = 0.5E+00*(press(:,:,k-1) + press(:,:,k))
       tflux(:,:,k) = 0.5E+00*(temp (:,:,k-1) + temp (:,:,k))
     enddo
     pflux(:,:,KERAD+1) = press(:,:,KERAD+1)
     tflux(:,:,KERAD+1) = temp (:,:,KERAD+1)
   endif
  
!-------------------------------------------------------------------
!   define deltaz in meters
!-------------------------------------------------------------------
   if (KSRAD == 1) then
     deltaz(:,:,KSRAD) = 2.0E+00*                                &
            1.0E-02*rgas*temp(:,:,KSRAD) / (grav*wtmair)
 
     do k = KSRAD+1,KERAD
       deltaz(:,:,k) = alog(pflux(:,:,k+1)/pflux(:,:,k))*        &
                       1.0E-02*rgas*temp(:,:,k) / (grav*wtmair)
     end do
   else
     do k = KSRAD,KERAD
       deltaz(:,:,k) = alog(pflux(:,:,k+1)/pflux(:,:,k))*        &
                       1.0E-02*rgas*temp(:,:,k) / (grav*wtmair)
     end do
   endif
!----------------------------------------------------------------------


end subroutine rad_layers



!####################################################################



	        end module rad_step_setup_mod

