module entrain_mod

!=======================================================================
!
!
!
!      K-PROFILE BOUNDARY LAYER SCHEME WITH CLOUD TOP ENTRAINMENT
!
!
!      January 2003
!      Contact person: Steve Klein
!
!      This routine calculates diffusivity coefficients for vertical
!      diffusion using a K-profile approach.  This scheme is modelled
!      after:
!
!      Lock, A.P., A.R. Brown, M.R. Bush, G.M. Martin, and R.N.B. Smith, 
!          2000: A new boundary layer mixing scheme. Part I: Scheme 
!          description and single-column modeling tests. Mon. Wea. Rev.,
!          128, 3187-3199.
!
!      The key part is the parameterization of entrainment at the top
!      convective layers.  This entrainment rate is parameterized as:
!
!      we = entrainment velocity 
!
!      we * delta_slv * grav / slv  = entrainment buoyancy consumption 
!
!           =  surface forcing + cloud top radiative cooling forcing
!
!           =  beta_surf*(u_star*b_star+(Ashear*(u_star**3)/zsml)) + 
!
!              beta_rad*grav*delta-F / rho/slv
! 
!-----------------------------------------------------------------------
!
! outside modules used
!

use      constants_mod, only: grav,vonkarm,cp_air,rdgas,rvgas,hlv,hls, &
                              tfreeze,radian 

use      utilities_mod, only: file_exist, open_file, error_mesg, FATAL,&
                              get_my_pe, close_file, read_data, &
			      write_data

use   diag_manager_mod, only: register_diag_field, send_data
        
use   time_manager_mod, only: time_type, get_date, month_name
 
use sat_vapor_pres_mod, only: lookup_es, lookup_des

use  monin_obukhov_mod, only: mo_diff

implicit none
private

!-----------------------------------------------------------------------
!
!      public interfaces

public entrain, entrain_init, entrain_end, entrain_tend, entrain_on

!-----------------------------------------------------------------------
!
!      global storage variable
!

real, allocatable, dimension(:,:,:) :: tdtlw  ! LW heating rate (K/s)

!-----------------------------------------------------------------------
!
!      set default values to namelist parameters       
!
real :: akmax       =  1.e4 ! maximum value for a diffusion coefficient 
                            ! (m2/s)
real :: wentrmax    =  0.05 ! maximum entrainment rate (m/s)
real :: parcel_buoy =  2.0  ! scaling factor for surface parcel buoyancy
real :: frac_inner  =  0.1  ! surface layer height divided by pbl height
real :: beta_surf   =  0.2  ! scaling of surface buoyancy flux for 
                            ! convective pbl entrainment
real :: Ashear      = 25.0  ! scaling of surface shear contribution to
                            ! entrainment
real :: beta_rad    =  0.23 ! entrainment scaling factor for radiative
                            ! cooling (from Lock et al. 2000)
real :: radfmin     = 30.   ! minimum radiative forcing for entrainment 
                            ! to be effective (W/m2)
real :: qdotmin     = 10.   ! minimum longwave cooling rate (K/day) for
                            ! entrainment to be effective, used only if
			    ! no layers were found using radfmin			    
real :: radperturb  =  0.3  ! parcel perturbation looking for depth of
                            ! radiatively driven layer (K)	
real :: critjump    =  0.3  ! critical jump for finding stable inter-
                            ! faces to bound convective layers, or to 
			    ! identify ambigous layers (K)
real :: zcldtopmax  =  3.e3 ! maximum altitude for cloud top of 
                            ! radiatively driven convection (m)			  		  
real :: pr          =  0.75 ! prandtl # (k_m/k_t) for radiatively driven 
                            ! convection and simple mixing scheme
real :: qamin       =  0.3  ! minimum cloud fraction for cloud top 
                            ! radiative cooling entrainment and kprofile
			    ! from radiative cooling to occur
real :: asympt_len  =  1.5e2! asymptotic mixing length for simple mixing
                            ! scheme
real :: rich_crit   =  0.25 ! critical richardson number for simple
                            ! mixing scheme
logical :: do_jump_exit = .true.
                            ! should an internal stable layer limit
			    ! the depth of the radiatively driven
			    ! convection?
logical :: prof_recon_on = .false.
                            ! should profile reconstruction be applied?			    
logical :: apply_entrain = .true. 
                            ! logical controlling whether results of
                            ! of entrainment module are applied:
	                    ! if F, then no diffusion coefficients 
			    ! from entrain_mod will be applied to
			    ! the actual diffusion coefficients;
		            ! i.e. the module is purely diagnostic

!-----------------------------------------------------------------------   
!
!  Stuff needed to write out extradiagnostics from a single point
!
			    
integer, dimension(2) :: ent_pts = 0 ! the global indices for i,j
                                     ! at which diagnostics will 
                                     ! print out
logical   :: do_print = .false.      ! should selected variables 
                                     ! be sent to logfile
logical   :: column_match = .false.  ! should this column be printed 
                                     ! out?
integer   :: dpu = 0                 ! unit # for do_print output
integer   :: n_print_levels = 14     ! how many of the lowest levels 
                                     ! should be printed out
				 
	
integer, parameter                 :: MAX_PTS = 20
integer, dimension (MAX_PTS)       :: i_entprt_gl=0, j_entprt_gl=0
real, dimension(MAX_PTS)           :: lat_entprt=999., lon_entprt=999.
integer                            :: num_pts_ij = 0
integer                            :: num_pts_latlon = 0

			  				   
namelist /entrain_nml/ wentrmax, parcel_buoy, frac_inner, beta_surf,   &
                       Ashear, beta_rad, radfmin, qdotmin, radperturb, &
		       critjump, zcldtopmax, pr, qamin, asympt_len,    &
		       do_jump_exit, prof_recon_on, apply_entrain,     &
		       ent_pts,  i_entprt_gl, j_entprt_gl, num_pts_ij, &
		       num_pts_latlon, lat_entprt, lon_entprt

integer     :: num_pts           !  total number of columns in which
                                 !  diagnostics are desired
		       
!-----------------------------------------------------------------------
!    deglon1 and deglat1 are the longitude and latitude of the columns
!    at which diagnostics will be calculated (degrees).
!-----------------------------------------------------------------------
real,    dimension(:), allocatable  :: deglon1, deglat1
 
!-----------------------------------------------------------------------
!    iradprt and jradprt are the processor-based i and j coordinates 
!    of the desired diagnostics columns.
!-----------------------------------------------------------------------
integer, dimension(:), allocatable  :: j_entprt, i_entprt

!-----------------------------------------------------------------------
!    do_raddg is an array of logicals indicating which latitude rows
!    belonging to the processor contain diagnostics columns.
!-----------------------------------------------------------------------
logical, dimension(:), allocatable  :: do_ent_dg

!-----------------------------------------------------------------------
!
!      diagnostic fields       
!

character(len=10) :: mod_name = 'entrain'
real              :: missing_value = 0.
integer           :: id_wentr_rad, id_wentr_pbl, id_radf,id_parcelkick,&
                     id_k_t_entr,  id_k_m_entr,  id_k_rad,  id_zsml,   &             
		     id_k_t_troen, id_k_m_troen, id_radfq,  id_pblfq,  &
		     id_zradbase,  id_zradtop,   id_vrad,   id_zradml, &
		     id_convpbl,   id_radpbl,    id_svpcp,  id_k_t_sim,&
		     id_k_m_sim,   id_simfq
       
      
!-----------------------------------------------------------------------
!
!      set default values to parameters       
!

logical         :: entrain_on = .false.
real, parameter :: small  = 1.e-4			      
real, parameter :: d608 = (rvgas-rdgas)/rdgas

!-----------------------------------------------------------------------
!
! declare version number 
!

character(len=128) :: Version = '$Id: entrain.F90,v 1.2 2003/04/09 20:55:21 fms Exp $'
character(len=128) :: Tag = '$Name: inchon $'
      
!-----------------------------------------------------------------------
!
! Subroutines include:
!
!      entrain         main driver program of the module
!
!      entrain_init    initialization routine       
!
!      entrain_tend    adds in the longwave heating rate to the 
!                      global storage variable
!
!      entrain_end     ending routine
!
!      pbl_depth       routine to calculate the depth of surface driven
!                      mixed layer
!
!      prof_recon      subroutine to perform profile reconstruction
!
!      radml_depth     subroutine to calculate the depth of the cloud
!                      topped radiatively driven mixed layer
!     
!      diffusivity_pbl subroutine to calculate diffusivity coefficients
!                      for surface driven mixed layer
!
!      diffusivity_sim subroutine to calculate simple diffusion 
!                      coefficients, used for all levels WHERE 
!                      surface driven unstable PBL and radiative driven
!                      convection do not occur

contains



!======================================================================= 
!
!      subroutine entrain_init 
!        
!
!      this subroutine reads the namelist file and restart data
!      and initializes some constants.
!        

subroutine entrain_init(lonb, latb, axes,time,idim,jdim,kdim)

!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
! 
!      idim,jdim,kdim    size of the first 3 dimensions 
!      axes, time        variables needed for netcdf diagnostics
!      latb, lonb        latitudes and longitudes at grid box boundaries
!
!
!      --------
!      internal
!      --------
! 
!      unit              unit number for namelist and restart file
!      io                internal variable for reading of namelist file
!      full              indices for full level axes coordinates
!      half              indices for half level axes coordinates
!
!-----------------------------------------------------------------------

integer,            intent(in) :: idim,jdim,kdim,axes(4)
type(time_type),    intent(in) :: time
real, dimension(:), intent(in) :: lonb, latb

integer                        :: unit,io
integer, dimension(3)          :: full = (/1,2,3/), half = (/1,2,4/)
integer                        :: nn, i, j, iloc, jloc
real                           :: dellat, dellon

!-----------------------------------------------------------------------
!
!      namelist functions

       If (File_Exist('input.nml')) Then
            unit = Open_File ('input.nml', action='read')
            io=1
            Do While (io .ne. 0)
                 Read  (unit, nml=entrain_nml, iostat=io, End=10)
            EndDo
  10        Call Close_File (unit)
       EndIf

       unit = Open_File ('logfile.out', action='append')
       if ( get_my_pe() == 0 ) then
            Write (unit,'(/,80("="),/(a))') trim(Version), trim(Tag)
            Write (unit,nml=entrain_nml)
       endif
       Call Close_File (unit)

!-----------------------------------------------------------------------
!
!      Stuff to do single point diagnostics


!      if (ent_pts(1) > 0 .and.ent_pts (2) >0) do_print = .true.

!      if (do_print) then
!!RSH  unit = Open_File ('entrain.out', action='write')
!      dpu  = Open_File ('entrain.out', threading='multi', action='write')
!      if ( get_my_pe() == 0 ) then
!           Write (unit,'(/,80("="),/(a))') trim(Version), trim(Tag)
!           Write (unit,nml=ent_nml)
!      endif
!!     Call Close_File (unit)
!      end if
       
!-----------------------------------------------------------------------
!    allocate and initialize a flag array which indicates the latitudes
!    containing columns where radiation diagnostics are desired.
!-----------------------------------------------------------------------
      allocate (do_ent_dg (size(latb)-1) )
      do_ent_dg(:) = .false.

!-----------------------------------------------------------------------
!    define the total number of points at which diagnostics are desired.
!    points may be specified either by lat-lon pairs or by global index
!    pairs. 
!-----------------------------------------------------------------------
      num_pts = num_pts_latlon + num_pts_ij

!-----------------------------------------------------------------------
!    continue on only if diagnostics are desired in at least one column.
!-----------------------------------------------------------------------
      if (num_pts > 0) then

!-----------------------------------------------------------------------
!    if more points are desired than space has been reserved for, print 
!    a message.
!-----------------------------------------------------------------------
        if (num_pts > MAX_PTS) then
          call error_mesg ( 'entrain_mod', &
         'must reset MAX_PTS or reduce number of diagnostics points',  &
                                                             FATAL)
        endif

!-----------------------------------------------------------------------
!    allocate space for arrays which will contain the lat and lon and
!    processor-local i and j indices.
!-----------------------------------------------------------------------
        allocate ( deglon1 (num_pts))
        allocate ( deglat1 (num_pts))
        allocate ( j_entprt (num_pts))
        allocate ( i_entprt (num_pts))

!-----------------------------------------------------------------------
!    if any points for diagnostics are specified by (i,j) global 
!    indices, determine their lat-lon coordinates. assumption is made 
!    that the deltas of latitude and longitude are uniform over 
!    the globe.
!-----------------------------------------------------------------------
        do nn=1,num_pts_ij
          dellat = latb(2) - latb(1)
          dellon = lonb(2) - lonb(1)
          lat_entprt(nn + num_pts_latlon) =     &
                      (-0.5*acos(-1.0) + (j_entprt_gl(nn) - 0.5)*  &
                                           dellat) * radian
          lon_entprt(nn + num_pts_latlon) =                & 
                       (i_entprt_gl(nn) - 0.5)*dellon*radian
        end do

!-----------------------------------------------------------------------
!    determine if the lat/lon values are within the global grid,
!    latitude between -90 and 90 degrees and longitude between 0 and
!    360 degrees.
!-----------------------------------------------------------------------
        do nn=1,num_pts
          j_entprt(nn) = 0
          i_entprt(nn) = 0
          deglat1(nn) = 0.0
          deglon1(nn) = 0.0
          if (lat_entprt(nn) .ge. -90. .and. &
              lat_entprt(nn) .le.  90.) then
          else
            call error_mesg ('entrain_mod', &
                ' invalid latitude for entrain diagnostics ', FATAL)
          endif

          if (lon_entprt(nn) .ge. 0. .and. &
              lon_entprt(nn) .le. 360.) then
          else
            call error_mesg ('entrain_mod', &
                ' invalid longitude for entrain diagnostics ', FATAL)
          endif

!-----------------------------------------------------------------------
!    determine if the diagnostics column is within the current 
!    processor's domain. if so, set a logical flag indicating the
!    presence of a diagnostic column on the particular row, define the 
!    i and j processor-coordinates and the latitude and longitude of 
!    the diagnostics column.
!-----------------------------------------------------------------------
          do j=1,size(latb) - 1
            if (lat_entprt(nn) .ge. latb(j)*radian .and.   &
                lat_entprt(nn) .lt. latb(j+1)*radian) then
              do i=1,size(lonb) - 1
                if (lon_entprt(nn) .ge. lonb(i)*radian     &
                                  .and.&
                    lon_entprt(nn) .lt. lonb(i+1)*radian)  &
                                   then
                  do_ent_dg(j) = .true.
                  j_entprt(nn) = j
                  i_entprt(nn) = i
                  deglon1(nn) = 0.5*(lonb(i) + lonb(i+1))*radian
                  deglat1(nn) = 0.5*(latb(j) + latb(j+1))*radian
                  exit
                endif
              end do
              exit
            endif
          end do
        end do

!-----------------------------------------------------------------------
!    open a unit for the entrain diagnostics output.
!-----------------------------------------------------------------------
        dpu = open_file ('entrain.out', action='write', &
                                 threading='multi', form='formatted')
       do_print = .true.
       if ( get_my_pe() == 0 ) then
            Write (dpu ,'(/,80("="),/(a))') trim(Version), trim(Tag)
            Write (dpu ,nml=entrain_nml)
       endif
      endif     ! (num_pts > 0)

!-----------------------------------------------------------------------
!
!      initialize entrain_on

       entrain_on = .TRUE.

       
!-----------------------------------------------------------------------
!
!      handle global storage
        
       if (allocated(tdtlw)) deallocate (tdtlw)
       allocate(tdtlw(idim,jdim,kdim))
                
       if (File_Exist('INPUT/entrain.res')) then
            unit = Open_File (FILE='INPUT/entrain.res', &
                 FORM='native', ACTION='read')
            call read_data (unit, tdtlw)
            call Close_File (unit)
       else
            tdtlw   (:,:,:) = 0.
       endif

!-----------------------------------------------------------------------
!
! register diagnostic fields       

       id_zsml = register_diag_field (mod_name, 'zsml', axes(1:2),     &
	    time, 'depth of surface well-mixed layer', 'm',            &
	    missing_value=missing_value )

       id_parcelkick = register_diag_field (mod_name, 'parcelkick',    &
            axes(1:2),time, 'surface parcel excess', 'K',              &
	    missing_value=missing_value )

       id_wentr_pbl = register_diag_field (mod_name,'wentr_pbl',       &
            axes(1:2), time,                                           &
	    'Entrainment velocity from surface buoyancy flux',         &
	    'meters per second', missing_value=missing_value )
       
       id_wentr_rad = register_diag_field (mod_name, 'wentr_rad',      &
            axes(1:2), time,                                           &
	    'Entrainment velocity from cloud top radiative cooling',   &
	    'meters per second', missing_value=missing_value )
       
       id_convpbl = register_diag_field (mod_name, 'convpbl_fq',       &
            axes(1:2), time, 'Frequency of convective boundary layer', &
	    'none')
       
       id_pblfq = register_diag_field (mod_name,'entr_pbl_fq',         &
            axes(half), time,                                          &
	    'Frequency of surface driven entrainment turbulent layer', &
	    'none')
       
       id_radfq = register_diag_field (mod_name, 'entr_rad_fq',        &
            axes(half), time,                                          &
	   'Frequency of radiative driven entrainment turbulent layer',&
	    'none')
       
       id_k_t_troen = register_diag_field (mod_name, 'k_t_troen',      &
            axes(half), time, 'Heat entrainment diffusivity from PBL', &
	    'meters squared per second', missing_value=missing_value )
       
       id_k_m_troen = register_diag_field (mod_name, 'k_m_troen',      &
            axes(half), time,                                          &
	    'Momentum entrainment diffusivity from PBL',               &
	    'meters squared per second', missing_value=missing_value )
       
       id_svpcp = register_diag_field (mod_name, 'svpcp',              &
            axes(1:2), time,                                           &
	    'Liquid water virtual potential temperature at cloud top', &
	    'K',missing_value=missing_value )

       id_zradbase = register_diag_field (mod_name, 'zradbase',        &
            axes(1:2), time,                                           &
	    'base of radiatively driven well-mixed layer', 'm',        &
	    missing_value=missing_value )

       id_zradtop = register_diag_field (mod_name, 'zradtop',          &
            axes(1:2), time,                                           &
	    'top of radiatively driven well-mixed layer', 'm',         &
	    missing_value=missing_value )

       id_zradml = register_diag_field (mod_name, 'zradml',            &
            axes(1:2), time,                                           &
	    'depth of radiatively driven well-mixed layer', 'm',       &
	    missing_value=missing_value )

       id_radpbl = register_diag_field (mod_name, 'radpbl_fq',         &
            axes(1:2), time,                                           &
	    'Frequency of radiatively driven turbulent layer',         &
	    'none' )
       
       id_vrad = register_diag_field (mod_name, 'vrad',                &
            axes(1:2), time,                                           &
	    'radiatively driven layer velocity scale', 'm',            &
	    missing_value=missing_value )

       id_k_rad = register_diag_field (mod_name, 'k_rad',              &
            axes(half), time,                                          &
	    'diffusivity from cloud top radiative cooling',            &
	    'meters squared per second', missing_value=missing_value )
       
       id_radf  = register_diag_field (mod_name, 'radf', axes(1:2),    &
	    time, 'Longwave jump in cloud top radiation',              &
	    'Watts/m**2', missing_value=missing_value )
 
       id_k_m_sim = register_diag_field (mod_name, 'k_m_sim',          &
            axes(half), time,                                          &
	    'Momentum entrainment diffusivity from simple mixing',     &
	    'meters squared per second', missing_value=missing_value )
       
       id_k_t_sim = register_diag_field (mod_name, 'k_t_sim',          &
            axes(half), time,                                          &
	    'Heat entrainment diffusivity from simple mixing',         &
	    'meters squared per second', missing_value=missing_value )
       
       id_simfq = register_diag_field (mod_name, 'entr_sim_fq',        &
            axes(half), time,                                          &
	   'Frequency of simple stable mixing scheme layer',           &
	    'none')
       
       id_k_t_entr = register_diag_field (mod_name, 'k_t_entr',        &
            axes(half), time,                                          &
	    'Heat diffusivity from entrainment module',                &
	    'meters squared per second')
       
       id_k_m_entr = register_diag_field (mod_name, 'k_m_entr',        &
            axes(half), time,                                          &
	    'Momentum diffusivity from entrainment module',            &
	    'meters squared per second' )
       
                        
!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine entrain_init

!
!======================================================================= 




!======================================================================= 
!
!      subroutine entrain
!        

subroutine entrain(is,ie,js,je,time,u_star,b_star,t,qv,ql,qi,qa,u,v,   &
                   zfull,pfull,zhalf,phalf,diff_m,diff_t,k_m_entr,     &
		   k_t_entr,kbot)

!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
!
!      is,ie,js,je  i,j indices marking the slab of model working on
!      time      variable needed for netcdf diagnostics
!
!      u_star    friction velocity (m/s)
!      b_star    buoyancy scale (m/s**2)
!
!      three dimensional fields on model full levels, reals dimensioned
!      (:,:,pressure), third index running from top of atmosphere to 
!      bottom
!          
!      t         temperature (K)
!      qv        water vapor specific humidity (kg vapor/kg air)
!      ql        liquid water specific humidity (kg cond/kg air)
!      qi        ice water specific humidity (kg cond/kg air)
!      qa        cloud fraction 
!      zfull     height of full levels (m)
!      pfull     pressure (Pa)
!      u         zonal wind (m/s)
!      v         meridional wind (m/s)
!
!      the following two fields are on the model half levels, with
!      size(zhalf,3) = size(t,3) +1, zhalf(:,:,size(zhalf,3)) 
!      must be height of surface (if you are not using eta-model)
!
!      zhalf     height at half levels (m)
!      phalf     pressure at half levels (Pa)
!
!      ------------
!      input/output
!      ------------
!
!      the following variables are defined at half levels and are
!      dimensions 1:nlev
!
!      diff_t   input and output heat diffusivity (m2/sec)
!      diff_m   input and output momentum diffusivity (m2/sec)
!
!      The diffusivity coefficient output from the routine includes
!      the modifications to use the internally calculated diffusivity
!      coefficients should apply_entrain = .true.  Otherwise, the
!      input and output values of the diffusivity coefficients will
!      be the same.
!
!      ------
!      output
!      ------
!
!      The following variables are defined at half levels and are
!      dimensions 1:nlev.
!
!      k_t_entr  heat diffusivity coefficient (m**2/s)
!      k_m_entr  momentum diffusivity coefficient (m**2/s)
!     
!      These are the diffusivity coefficients calculated by this 
!      routine.  This is sent out for diagnostics purposes only.
!      The diffusivity coefficients actually used by the model are
!      the output values of diff_t and diff_m.
!
!      --------------
!      optional input
!      --------------
!
!      kbot      integer indicating the lowest true layer of atmosphere
!                this is used only for eta coordinate model
!
!      --------
!      internal
!      --------
!
!
!      General variables
!      -----------------
!
!      zsurf       height of surface (m)
!      zfull_ag    height of full model levels above the surface (m)
!      slv         virtual static energy (J/kg)      
!      density     air density (kg/m3)
!      hleff       effective latent heat of vaporization/sublimation 
!                  (J/kg)
!      mask        real array indicating the point is above the surface
!                  if equal to 1.0 and indicating the point is below 
!                  the surface if equal to 0. (used for eta coordinate 
!                  model)
!      zhalf_ag    height of half model levels above the surface (m)
!
!
!
!      Variables related to surface driven convective layers
!      -----------------------------------------------------
!
!      zsml        height of surface driven mixed layer (m)
!      parcelkick  buoyancy kick to surface parcel (K)
!      wentr_pbl   surface driven entrainment rate (m/s)
!      convpbl     1 is surface driven convective layer present
!                  0 otherwise
!      pblfq       1 if the half level is part of a surface driven
!                  layer, 0 otherwise
!      k_m_troen   momentum diffusion coefficient (m2/s)
!      k_t_troen   heat diffusion coefficient (m2/s)
!
!
!      Variables related to cloud top driven radiatively driven layers
!      ---------------------------------------------------------------
!
!      zradbase    height of base of radiatively driven mixed layer (m)
!      zradtop     height of top of radiatively driven mixed layer (m)
!      zradml      depth of radiatively driven mixed layer (m)
!      vrad        radiatively driven velocity scale (m/s)
!      radf        longwave jump at cloud top (W/m2) -- the radiative 
!                  forcing for cloud top driven mixing.
!      wentr_rad   cloud top driven entrainment (m/s)
!      svpcp       cloud top value of liquid water virtual static energy 
!                  divided by cp (K)
!      radpbl      1 if cloud top radiatively driven layer is present
!                  0 otherwise
!      radfq       1 if the half level is part of a radiatively driven
!                  layer, 0 otherwise
!      k_rad       radiatively driven diffusion coefficient (m2/s)
!
!
!      Variables related to simple mixing scheme
!      ------------------------------------------
!
!      k_m_sim     momentum diffusion coefficient (m2/s)
!      k_t_sim     heat diffusion coefficient (m2/s)
!      simfq       frequency of stable simple mixing scheme
!
!-----------------------------------------------------------------------

integer,         intent(in)                      :: is,ie,js,je
type(time_type), intent(in)                      :: time
real,            intent(in),    dimension(:,:)   :: u_star,b_star
real,            intent(in),    dimension(:,:,:) :: t,qv,ql,qi,qa
real,            intent(in),    dimension(:,:,:) :: u,v,zfull,pfull
real,            intent(in),    dimension(:,:,:) :: zhalf, phalf
real,            intent(inout), dimension(:,:,:) :: diff_m,diff_t
real,            intent(out),   dimension(:,:,:) :: k_m_entr,k_t_entr
integer,  intent(in),   dimension(:,:), optional :: kbot


integer                                         :: i,j,k,ibot,itmp
integer                                         :: nlev,nlat,nlon,ipbl
integer                                         :: kmax,kcldtop
logical                                         :: used
real                                            :: maxradf, tmpradf
real                                            :: maxqdot, tmpqdot
real                                            :: wentr_tmp
real                                            :: k_entr_tmp,tmpjump
real                                            :: dslvcptmp,ztmp
real, dimension(size(t,1),size(t,2))            :: zsurf,zsml,parcelkick
real, dimension(size(t,1),size(t,2))            :: zradbase,zradtop
real, dimension(size(t,1),size(t,2))            :: vrad,radf,svpcp
real, dimension(size(t,1),size(t,2))            :: zradml
real, dimension(size(t,1),size(t,2))            :: wentr_rad,wentr_pbl
real, dimension(size(t,1),size(t,2))            :: convpbl, radpbl
real, dimension(size(t,1),size(t,2),size(t,3))  :: slv, density
real, dimension(size(t,1),size(t,2),size(t,3))  :: mask,hleff
real, dimension(size(t,1),size(t,2),size(t,3))  :: zfull_ag
real, dimension(size(t,1),size(t,2),size(t,3)+1):: zhalf_ag
real, dimension(size(t,1),size(t,2),size(t,3))  :: radfq,pblfq,simfq
real, dimension(size(t,1),size(t,2),size(t,3)+1):: mask3,rtmp
real, dimension(size(t,1),size(t,2),size(t,3))  :: k_m_troen,k_t_troen
real, dimension(size(t,1),size(t,2),size(t,3))  :: k_rad,k_t_sim,k_m_sim
real, dimension(size(t,3))                      :: use_entr

integer                                         :: ipt,jpt
integer, dimension(MAX_PTS) :: nsave
integer :: iloc(MAX_PTS), jloc(MAX_PTS), nn, kk, npts, nnsave
integer :: year, month, day, hour, minute, second
character(len=16) :: mon

!-----------------------------------------------------------------------
!
!      open ascii entrain.out file if do_print

!      if (do_print) then
!           if ( ent_pts(1) >= is .and. ent_pts(1) <= ie  .and.        &
!                ent_pts(2) >= js .and. ent_pts(2) <= je) then
!                ipt=ent_pts(1)-is+1
!	 jpt=ent_pts(2)-js+1
!!RSH            dpu = open_file ('entrain.out', action='append')
!           else
!         ipt = 0
!	 jpt = 0
!           endif
!      endif
                    
!-----------------------------------------------------------------------
!
!      initialize variables

       convpbl    = 0.0
       wentr_pbl  = missing_value
       pblfq      = 0.0
       ipbl       = 0
       zsml       = missing_value
       parcelkick = missing_value
       k_t_troen  = missing_value
       k_m_troen  = missing_value
       radpbl     = 0.0
       svpcp      = missing_value
       zradbase   = missing_value
       zradtop    = missing_value
       zradml     = missing_value
       vrad       = missing_value
       radf       = missing_value
       radfq      = 0.0
       wentr_rad  = missing_value
       k_rad      = missing_value
       simfq      = 0.0
       k_t_sim    = missing_value
       k_m_sim    = missing_value
       k_t_entr   = 0.0
       k_m_entr   = 0.0
                  
!-----------------------------------------------------------------------
!
!      compute height above surface
!

       nlev = size(t,3)
       nlat = size(t,2)
       nlon = size(t,1)
       
       mask = 1.0
                   
       if (present(kbot)) then
            do j=1,nlat
            do i=1,nlon
                 zsurf(i,j) = zhalf(i,j,kbot(i,j)+1)
            enddo
            enddo
       else
            zsurf(:,:) = zhalf(:,:,nlev+1)
       end if

       do k = 1, nlev
            zfull_ag(:,:,k) = zfull(:,:,k) - zsurf(:,:)
            zhalf_ag(:,:,k) = zhalf(:,:,k) - zsurf(:,:)
       end do
       zhalf_ag(:,:,nlev+1) = zhalf(:,:,nlev+1) - zsurf(:,:)
       
     	
!-----------------------------------------------------------------------
!
!      set up specific humidities and static energies  	
!      compute airdensity
!

       hleff   = (min(1.,max(0.,0.05*(t       -tfreeze+20.)))*hlv + &
                  min(1.,max(0.,0.05*(tfreeze -t          )))*hls)
     
       slv     = cp_air*t + grav*zfull_ag - hleff*(ql + qi)
       slv     = slv*(1+d608*(qv+ql+qi))
       density = pfull/rdgas/(t *(1.+d608*qv-ql-qi))              


!-----------------------------------------------------------------------
!      
!      compute simple mixing scheme
!

       call diffusivity_sim(slv/cp_air, u, v, zfull_ag, zhalf_ag,             &
                            k_m_sim, k_t_sim)
			    
!-----------------------------------------------------------------------
! 
!      big loop over points
!
       
       ibot = nlev
              
       do j=1,nlat	
         npts = 0
       if (do_ent_dg(j+js-1) ) then       
	 do nn=1,num_pts
         if (                          &
	       js == j_entprt(nn) .and.  &
                 i_entprt(nn) >= is .and. i_entprt(nn) <= ie) then
	      iloc(npts+1) = i_entprt(nn) - is + 1
	      jloc(npts+1) = j_entprt(nn) - js + 1
	      nsave(npts+1) = nn
	      npts = npts + 1
         endif
        end do    ! (num_points)
       else
          ipt = 0
           jpt = 0
	   column_match = .false.
       endif
       do i=1,nlon
          if (npts > 0) then
		   do nn=1,npts

		   ipt = iloc(nn)
		   jpt = jloc(nn)
!                  if (i == ipt .and. j == jpt) then
                   if (i == ipt ) then
!	 if (i == ipt .and. j == jpt .and.  &
!	     js == j_entprt(nn) ) then
	         column_match = .true.
!	 nnsave = nn
		 nnsave = nsave(nn)
		 exit
		 else
	         column_match = .false.
		 endif
               end do
	       nn = nnsave
	       else 
	         column_match = .false.
		 nn = 0
	       endif
            !-----------------------------------------------------------
	    !
            ! should diagnostics be printed out for this column
            !
	    
!    if (i .eq. ipt .and. j .eq. jpt) then
!         column_match = .true.
!           else
!         column_match = .false.		 
!           end if

            if (column_match) then
	    call get_date(Time, year, month, day, hour, minute, second)
	    mon = month_name (month)
            write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' '
            write (dpu,'(a)')  '===================================='//&
	         '=================='
            write (dpu,'(a)')  ' '
            write (dpu,'(a)')  '               ENTERING ENTRAIN    '
            write (dpu,'(a)')  ' '
            write (dpu,'(a, i6,a,i4,i4,i4,i4)')  ' time stamp:',   &
	                                       year, trim(mon), day, &
	                                       hour, minute, second
            write (dpu,'(a)')  '  DIAGNOSTIC POINT COORDINATES :'
            write (dpu,'(a)')  ' '
	    write (dpu,'(a,f8.3,a,f8.3)') ' longitude = ', deglon1(nn),&
	                                  ' latitude  = ', deglat1(nn)
	    write (dpu,'(a,i6,a,i6)')    &
                                   ' global i =', i_entprt_gl(nn), &
                                   ' global j = ', j_entprt_gl(nn)
	    write (dpu,'(a,i6,a,i6)')    &
                                   ' processor i =', i_entprt(nn),     &
                                   ' processor j = ',j_entprt(nn)
	    write (dpu,'(a,i6,a,i6)')     &
                                   ' window    i =', ipt,          &
                                   ' window    j = ',jpt
            write (dpu,'(a)')  ' '
	    write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' k      T         u         v       '//&
	         '  qv        qt ' 
            write (dpu,'(a)')  '       (K)      (m/s)     (m/s)     '//&
	         '(g/kg)    (g/kg)'
            write (dpu,'(a)')  '------------------------------------'//&
	         '-----------------'
            write (dpu,'(a)')  ' '
            do kk = nlev-n_print_levels,nlev
                 write(dpu,18) kk,t(i,j,kk),u(i,j,kk),v(i,j,kk),       &
		      1000.*qv(i,j,kk), 1000.*(qv(i,j,kk)+ql(i,j,kk)+  &
		      qi(i,j,kk))
            end do
18          format(1X,i2,1X,5(f9.4,1X))
            write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' k      qa        qc      sliv/cp_air'//&
	         'density    tdtlw'
            write (dpu,'(a)')  '                (g/kg)     (K)      '//&
	         ' (kg/m3)   (K/day)'
            write (dpu,'(a)')  '------------------------------------'//&
	         '-------------------'
            write (dpu,'(a)')  ' '
            do kk = nlev-n_print_levels,nlev
                 write(dpu,19) kk,qa(i,j,kk),1000.*                    &
		   (ql(i,j,kk)+qi(i,j,kk)),slv(i,j,kk)/cp_air,         &
		   density(i,j,kk),tdtlw(is-1+i,js-1+j,kk)*86400.
            enddo	    
19          format(1X,i2,1X,5(f9.4,1X))
            write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' k   z_full    z_half    p_full    p'//&
	         '_half  '
            write (dpu,'(a)')  '      (m)      (m)        (mb)      '//&
	         '(mb)'
            write (dpu,'(a)')  '------------------------------------'//&
	         '--------'
            write (dpu,'(a)')  ' '
            do kk = nlev-n_print_levels,nlev
                 write(dpu,619) kk,zfull_ag(i,j,kk),zhalf_ag(i,j,kk+1),&
	              pfull(i,j,kk)/100.,phalf(i,j,kk+1)/100.
            enddo
619         format(1X,i2,1X,4(f9.4,1X))
            write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' '
            end if
	    
            !-----------------------------------------------------------
	    !  BEGIN BASIC CODE
	    
	    use_entr = 0.
	    if (present(kbot)) ibot = kbot(i,j)
	    
	    !-----------------------------------------------------------
	    !
	    ! SURFACE DRIVEN CONVECTIVE LAYERS
	    !
	    ! Note this part is done only if b_star > 0., that is,
	    ! upward surface buoyancy flux
	 
	    
	    if (b_star(i,j) .gt. 0.) then
	    
	         call pbl_depth(slv(i,j,1:ibot)/cp_air, &
		      zfull_ag(i,j,1:ibot),u_star(i,j),b_star(i,j),    &
		      ipbl,zsml(i,j),parcelkick(i,j))
            
                 
		 !------------------------------------------------------
		 ! Following Lock et al. 2000, limit height of surface
		 ! well mixed layer if interior stable interface is
		 ! found.  An interior stable interface is diagnosed if
		 ! the slope between 2 full levels is greater than
		 ! the namelist parameter critjump
		  
		 if (ipbl .lt. ibot) then 
                      do k = ibot, ipbl+1, -1
		           tmpjump =(slv(i,j,k-1)-slv(i,j,k))/cp_air 
		           if (tmpjump .gt. critjump) then
			        ipbl = k
				zsml(i,j) = zhalf_ag(i,j,ipbl)
                                exit
                           end if
                      enddo
                 end if		      			   				

                 !-------------------------------------
		 ! compute entrainment rate
		 !
		 
                 wentr_tmp= min(wentrmax,max(0.,beta_surf*             &
		      (u_star(i,j)*b_star(i,j)  +                      &
		      (Ashear*(u_star(i,j)**3.)/zsml(i,j)))*           &
		      slv(i,j,ipbl)/grav/                              &
		      max((slv(i,j,ipbl-1)-slv(i,j,ipbl)),0.1) ))  
                                           
                 k_entr_tmp = wentr_tmp*(zfull_ag(i,j,ipbl-1)-         &
			                 zfull_ag(i,j,ipbl))  
                 k_entr_tmp = min ( k_entr_tmp, akmax )
			
                 pblfq(i,j,ipbl:ibot) = 1.
		 convpbl(i,j)         = 1.
		 use_entr(ipbl:ibot)  = 1.
	         wentr_pbl(i,j)       = wentr_tmp
                 k_t_troen(i,j,ipbl)  = k_entr_tmp
		 k_m_troen(i,j,ipbl)  = k_entr_tmp
		 k_t_entr (i,j,ipbl)  = k_t_entr(i,j,ipbl) + k_entr_tmp
	         k_m_entr (i,j,ipbl)  = k_m_entr(i,j,ipbl) + k_entr_tmp
	         
	         if (ipbl .lt. ibot) then
		      
		      call diffusivity_pbl(zsml(i,j),u_star(i,j),      &
			   b_star(i,j), slv(i,j,(ipbl+1):ibot)/cp_air, &
		                  zhalf_ag(i,j,(ipbl+1):ibot),        &
			          k_m_troen(i,j,(ipbl+1):ibot),        &
			          k_t_troen(i,j,(ipbl+1):ibot))
                		      
		      k_t_entr(i,j,(ipbl+1):ibot) =                    & 
		           k_t_entr(i,j,(ipbl+1):ibot) +               &
			   k_t_troen(i,j,(ipbl+1):ibot)
		     
		      k_m_entr(i,j,(ipbl+1):ibot) =                    & 
		           k_m_entr(i,j,(ipbl+1):ibot) +               &
			   k_m_troen(i,j,(ipbl+1):ibot)
		         
                 end if
	         
            end if		      
	    
	    
	    !-----------------------------------------------------------
	    !
	    ! LW RADIATIVELY DRIVEN CONVECTIVE LAYERS
	    !
	    ! This part is done only if a level can be found with
	    ! greater than radfmin (typically 30 W/m2) longwave
	    ! divergence and if the level is at a lower altitude
	    ! than zcldtopmax (typically 3000 m).
	    !
	    ! Note that if no layer is found with radiative divergence
	    ! greater than radfmin, a check is made on the heating rates
	    ! themselves and compared to qdotmin (10 K/day)

            !--------------------------
	    ! find level of zcldtopmax
            kmax = ibot+1
	    do k = 1, ibot
                if( zhalf_ag(i,j,k) < zcldtopmax) then
                     kmax = k
                     exit
                end if
            end do
	   	    
	    !-----------------------------------------------------------
	    ! compute radiative driving
	    !
	    ! look at heating rate itself if no levels with radiative 
	    ! forcing greater than radfmin are found.
	    
	    kcldtop = ibot+1
	    maxradf = radfmin
            do k = kmax, ibot
	         tmpradf = -1.*tdtlw(is-1+i,js-1+j,k)*cp_air*              &
		           ((phalf(i,j,k+1)-phalf(i,j,k))/grav)
                 if (tmpradf .gt. maxradf) then
		      kcldtop = k
		      maxradf = tmpradf
		 end if       			   
            enddo              
	    
	    !-----------------------------------------------------------
	    ! Second try to find a radiatively driven level
	    !
	    ! exit if no levels with longwave cooling rate greater than
	    ! qdotmin are found.
	    
	    if (kcldtop .eq. ibot+1) then
	    
	         kcldtop = ibot+1
	         maxradf = radfmin
		 maxqdot = qdotmin
                 
		 do k = kmax, ibot
		      tmpqdot = -1.*86400*tdtlw(is-1+i,js-1+j,k)
	              tmpradf = -1.*tdtlw(is-1+i,js-1+j,k)*cp_air*         &
		           ((phalf(i,j,k+1)-phalf(i,j,k))/grav)
                      if (tmpqdot .gt. maxqdot) then
		           kcldtop = k
		           maxradf = tmpradf
			   maxqdot = tmpqdot
		      end if       			   
                 enddo              
                 
            if (kcldtop .eq. ibot+1) go to 55
	   
	    end if
	   
	    !-----------------------------------------------------------
	    ! following Lock for stable layer one level down;
	    ! move cld top there if it exists
	    	    
	    !if (kcldtop .lt. ibot) then
	    !    tmpjump =(slv(i,j,kcldtop)-slv(i,j,kcldtop+1))/cp_air
            !    if (tmpjump .gt. critjump) kcldtop = kcldtop+1		           
            !end if	    
	    
	    !-----------------------------------------------------------
	    ! if layer is unstable move up a layer
	    ! if that layer is also unstable exit
	    
	    if (slv(i,j,kcldtop-1) .lt. slv(i,j,kcldtop)) then  
	         kcldtop = kcldtop - 1
		 if (slv(i,j,kcldtop-1) .lt. slv(i,j,kcldtop)) then
		      go to 55
		 end if     
            end if
	    		 
            !-----------------------------------------------------------
            ! exit if no cloud is present

            if ( qa(i,j,kcldtop-1)           .lt. qamin .and.          &
	         qa(i,j,kcldtop  )           .lt. qamin .and.          &
		 qa(i,j,min(kcldtop+1,ibot)) .lt. qamin ) go to 55
		  	    		 
            !-----------------------------------------------------------
	    ! determine if kcldtop is in an ambiguous layer
            !
            ! if not, assume inversion is at zhalf_ag(kcldtop) and
	    ! proceed normally.
	    !
	    ! if so, then do profile reconstruction. 
	    !
	    ! NOTE THAT	THE ENTRAINMENT RATE IS NOT APPLIED TO THE 
	    ! DIFFUSION COEFFICIENTS IF KCLDTOP IS IN AN AMBIGUOUS 
	    ! LAYER.
	    
	    radf(i,j) = maxradf
	    
	    if ( prof_recon_on .and.  kcldtop .lt. ibot .and.          &
	         ((slv(i,j,kcldtop)-slv(i,j,kcldtop+1))/cp_air) .gt.   &
		 critjump) then
                 
		 !---------------------------
	         ! do profile reconstruction
	    
                 call prof_recon(density(i,j,kcldtop),                 &
		                 slv(i,j,(kcldtop-2):(kcldtop+1))/cp_air, &
		              pfull(i,j,(kcldtop-2):(kcldtop+1)),      &
		              phalf(i,j, kcldtop   :(kcldtop+1)),      &
		              zradtop(i,j),dslvcptmp)		            
		 
		 zradtop(i,j)   = zradtop(i,j)+zhalf_ag(i,j,kcldtop+1)
		 wentr_rad(i,j) = min(wentrmax,beta_rad*maxradf/cp_air/    &
		                max(dslvcptmp,0.1)/density(i,j,kcldtop))
		     
                 svpcp(i,j)         = slv(i,j,kcldtop+1)/cp_air
		 radpbl(i,j)        = 1.
                 use_entr(kcldtop)  = 1.
	    
            else
	    	    
	         !---------------------------
		 ! compute entrainment rate
	         !
	         ! compute cloud top temperature, svpcp. Ensure that
		 ! that the cloud top parcel temperature, equal to
		 ! svpcp minus radperturb, is no greater than
		 ! the temperature of level kcldtop+1. This is
		 ! done so that if kcldtop is in an ambiguous layer 
		 ! then there will not be a sudden jump in cloud top
		 ! properties
		 
		 svpcp(i,j) =min((slv(i,j,kcldtop  )/cp_air),              &
		                 (slv(i,j,kcldtop+1)/cp_air) ) 
	    
	         dslvcptmp = (slv(i,j,kcldtop-1)/cp_air)-svpcp(i,j)
	    
	         wentr_rad(i,j) = min(wentrmax,beta_rad*maxradf/cp_air/    &
		                max(dslvcptmp,0.1)/density(i,j,kcldtop))
		 
                 k_entr_tmp = wentr_rad(i,j)*                          &
		       (zfull_ag(i,j,kcldtop-1)-zfull_ag(i,j,kcldtop))
                 k_entr_tmp = min ( k_entr_tmp, akmax )

                 zradtop(i,j) = zhalf_ag(i,j,kcldtop)
                 
	         radfq(i,j,kcldtop) = 1.
	         radpbl(i,j)        = 1.
                 use_entr(kcldtop)  = 1.
	         k_rad(i,j,kcldtop) = k_entr_tmp
	         k_t_entr (i,j,kcldtop) = k_t_entr(i,j,kcldtop)        &
		                        + k_entr_tmp
	         k_m_entr (i,j,kcldtop) = k_m_entr(i,j,kcldtop)        &
		                        + pr*k_entr_tmp
	    
	    end if

            		                  		 
            !-----------------------------------------------------------
	    ! find depth of radiatively driven convection 
 
            if (kcldtop .lt. ibot) then 
                 call radml_depth(svpcp(i,j),zradtop(i,j),             &
		      slv(i,j,kcldtop:ibot)/cp_air,                        &
		      zfull_ag(i,j,kcldtop:ibot),                      &
		      zhalf_ag(i,j,kcldtop:ibot),zradbase(i,j),        &
		      zradml(i,j))	      
            else
	         zradbase(i,j) = 0.0
                 zradml(i,j)   = zradtop(i,j)  
            end if                   
                 
            !-----------------------------------------------------------
	    ! compute radiation driven scale
	    !
	    ! Vrad**3 = g*zradml*radf/density/slv
	    
	    vrad(i,j) = grav*zradml(i,j)*maxradf/density(i,j,kcldtop)/ &
	           slv(i,j,kcldtop)
		   
            vrad(i,j) = vrad(i,j) ** (1./3.)		   
	    		 
            !-----------------------------------------------------------
	    ! if there are any interior layers to calculate diffusivity
	      
	    if ( kcldtop .lt. ibot ) then 		 						 

                 do k = kcldtop+1,ibot
		 
		     ztmp = max(0.,(zhalf_ag(i,j,k)-zradbase(i,j))/    &
		                    zradml(i,j) )
	             
		     if (ztmp.gt.0.) then
		     
		          radfq(i,j,k) = 1.
                          use_entr(k)  = 1.
			  k_entr_tmp = 0.85*vonkarm*vrad(i,j)*ztmp*    &
		                  zradml(i,j)*ztmp*((1.-ztmp)**0.5)
	                  k_entr_tmp = min ( k_entr_tmp, akmax )
	                  k_rad(i,j,k) = k_entr_tmp
	                  k_t_entr (i,j,k) = k_t_entr(i,j,k)           &
			                   + k_entr_tmp
	                  k_m_entr (i,j,k) = k_m_entr(i,j,k)           &
			                   + pr*k_entr_tmp
	             
		     end if		                  		 
                enddo  
		 
            end if
	
55          continue
	    
	    !-----------------------------------------------------------
	    !
	    ! Patch in simple mixing coefficients where surface convec-
	    ! tive pbl and radiatively driven pbl are not giving you
	    ! a coefficient
	    
	    do k = 2, ibot
	         if (use_entr(k).eq.0. .and.(k_t_sim(i,j,k).gt.0. .or. &
		                             k_m_sim(i,j,k).gt.0.)) then
		      simfq(i,j,k)    = 1.
		      k_t_entr(i,j,k) = k_t_sim(i,j,k)
		      k_m_entr(i,j,k) = k_m_sim(i,j,k)
	         else
		      simfq(i,j,k)    = 0.
		      k_t_sim(i,j,k)  = 0.
		      k_m_sim(i,j,k)  = 0.
		 end if
            enddo		 

	    !-----------------------------------------------------------
	    !
	    ! Modify diffusivity coefficients if apply_entrain = T
	    	    
            if (apply_entrain) then	   
	         diff_t(i,j,:) = k_t_entr(i,j,:)
	         diff_m(i,j,:) = k_m_entr(i,j,:)
            end if
	    
            !-----------------------------------------------------------
	    !
	    ! Selected points printout
	    

            if (column_match) then
	    
	    write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' '
	    write (dpu,'(a,f10.4)')  ' u_star = ', u_star(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' b_star = ', b_star(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' '
	    write (dpu,'(a,f10.4)')  ' convpbl = ', convpbl(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a,i3)')  ' ipbl = ', ipbl
            write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' parcelkick = ', parcelkick(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' zsml = ', zsml(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' wentr_pbl = ', wentr_pbl(i,j)
            write (dpu,'(a)')  ' '
	    write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' radpbl = ', radpbl(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a,i3)')  ' kcldtop = ', kcldtop
            write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' svpcp = ', svpcp(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' zradtop = ', zradtop(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' zradbase = ', zradbase(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' wentr_rad = ', wentr_rad(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' vrad = ', vrad(i,j)
            write (dpu,'(a)')  ' '
            write (dpu,'(a,f10.4)')  ' radf = ', radf(i,j)
            write (dpu,'(a)')  ' '
       
       
            write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' k  use_entr    diff_t    diff_m  '//  &
	         '  k_t_entr  k_m_entr' 
            write (dpu,'(a)')  '                (m2/s)    (m2/s)    '//&
	         '(m2/s)    (m2/s)'
            write (dpu,'(a)')  '------------------------------------'//&
	         '------------------'
            write (dpu,'(a)')  ' '
            do kk = nlev-n_print_levels,nlev
                 write(dpu,947) kk,use_entr(kk), diff_t(i,j,kk),       &
		      diff_m(i,j,kk), k_t_entr(i,j,kk), k_m_entr(i,j,kk)
            end do
947         format(1X,i2,3X,f3.1,4X,4(f9.4,1X))
            write (dpu,'(a)')  ' '
	    write (dpu,'(a)')  ' '
            write (dpu,'(a)')  ' k   pblfq  radfq  simfq     '//&
	         'k_t_troen k_m_troen   k_rad    k_t_sim   k_m_sim' 
            write (dpu,'(a)')  '                              (m2/s)'//   &
	         '    (m2/s)     (m2/s)    (m2/s)    (m2/s) '
            write (dpu,'(a)')  '------------------------------------'//&
	         '-----------------------------------------------'
            write (dpu,'(a)')  ' '
            do kk = nlev-n_print_levels,nlev
                 write(dpu,949) kk,pblfq(i,j,kk),radfq(i,j,kk),        &
		     simfq(i,j,kk),k_t_troen(i,j,kk),k_m_troen(i,j,kk),&
		     k_rad(i,j,kk),k_t_sim(i,j,kk),k_m_sim(i,j,kk)
            end do
949         format(1X,i2,3X,3(f3.1,4X),5(f9.4,1X))
            end if

      enddo
      enddo
	    
!-----------------------------------------------------------------------
! 
!      Diagnostics
!

       if ( id_convpbl    > 0 .or. id_wentr_pbl > 0 .or.               &
            id_k_m_entr   > 0 .or. id_wentr_rad > 0 .or.               &
            id_zsml       > 0 .or. id_k_t_entr  > 0 .or.               &
	    id_pblfq      > 0 .or. id_radfq     > 0 .or.               &
	    id_parcelkick > 0 .or. id_k_t_troen > 0 .or.               &
	    id_k_rad      > 0 .or. id_k_m_troen > 0 .or.               &
	    id_radf       > 0 .or. id_vrad      > 0 .or.               &
	    id_zradbase   > 0 .or. id_radpbl    > 0 .or.               &
	    id_zradtop    > 0 .or. id_zradml    > 0 .or.               &
	    id_svpcp      > 0 .or. id_k_t_sim   > 0 .or.               &
	    id_k_m_sim    > 0 .or. id_simfq     > 0 ) then 

            mask3(:,:,1:(nlev+1)) = 1.
            if (present(kbot)) then
                 where (zhalf_ag < 0.)
                      mask3(:,:,:) = 0.
                 end where
            endif
            
	    !----------------------------------------------
	    !
	    ! Convective PBL diagnostics
	    !
	    
	    if ( id_convpbl > 0 ) then
                 used = send_data ( id_convpbl, convpbl, time, is, js )
            end if
	    
	    if ( id_zsml > 0 ) then
                 used = send_data ( id_zsml, zsml, time, is, js )
            end if
	    
	    if ( id_parcelkick > 0 ) then
                 used = send_data ( id_parcelkick, parcelkick, time,   &
		                    is, js )
            end if
	    
	    if ( id_wentr_pbl > 0 ) then
                 used = send_data ( id_wentr_pbl, wentr_pbl, time,     &
		                    is, js)
            end if
            
	    if ( id_pblfq > 0 ) then
                 rtmp = 0.
	         rtmp(:,:,1:nlev) = pblfq
	         used = send_data ( id_pblfq, rtmp, time, is, js,      &
		                    1, rmask=mask3 )
            end if
            
	    if ( id_k_t_troen > 0 ) then
                 rtmp = 0.
	         rtmp(:,:,1:nlev) = k_t_troen
	         used = send_data ( id_k_t_troen, rtmp, time, is, js,  &
		                    1, rmask=mask3 )
            end if
            
	    if ( id_k_m_troen > 0 ) then
                 rtmp = 0.
	         rtmp(:,:,1:nlev) = k_m_troen
	         used = send_data ( id_k_m_troen, rtmp, time, is, js,  &
		                    1, rmask=mask3 )
            end if
            
	    !----------------------------------------------
	    !
	    ! Cloud top radiative cooling diagnostics
	    !
	    
	    if ( id_radpbl > 0 ) then
                 used = send_data ( id_radpbl, radpbl, time, is, js )
            end if
	    
	    if ( id_vrad > 0 ) then
                 used = send_data ( id_vrad, vrad, time, is, js )
            end if
	    
	    if ( id_zradml > 0 ) then
                 used = send_data ( id_zradml, zradml, time, is, js )
            end if
	    
	    if ( id_svpcp > 0 ) then
                 used = send_data ( id_svpcp, svpcp, time, is, js )
            end if
	    
	    if ( id_zradtop > 0 ) then
                 used = send_data ( id_zradtop, zradtop, time, is, js )
            end if
	    
	    if ( id_zradbase > 0 ) then
                 used = send_data ( id_zradbase, zradbase, time, is, js)
            end if
	    
	    if ( id_radf > 0 ) then
                 used = send_data ( id_radf, radf, time, is, js )            
            end if
            
	    if ( id_wentr_rad > 0 ) then
                 used = send_data ( id_wentr_rad, wentr_rad, time, is, &
		                    js)
            end if
            
	    if ( id_radfq > 0 ) then
                 rtmp = 0.
	         rtmp(:,:,1:nlev) = radfq
	         used = send_data ( id_radfq, rtmp, time, is, js,      &
		                    1, rmask=mask3 )
            end if
            
	    if ( id_k_rad > 0 ) then
                 rtmp = 0.
	         rtmp(:,:,1:nlev) = k_rad
	         used = send_data ( id_k_rad, rtmp, time, is, js, 1,   &
		                    rmask=mask3 )
            end if
            
	    !----------------------------------------------
	    !
	    ! Simple mixing diffusion coefficients
	    !
	    
	    if ( id_k_t_sim > 0 ) then
                 rtmp = 0.
	         rtmp(:,:,1:nlev) = k_t_sim
	         used = send_data ( id_k_t_sim, rtmp, time, is, js,  &
		                    1, rmask=mask3 )
            end if
            
	    if ( id_k_m_sim > 0 ) then
                 rtmp = 0.
	         rtmp(:,:,1:nlev) = k_m_sim
	         used = send_data ( id_k_m_sim, rtmp, time, is, js,  &
		                    1, rmask=mask3 )
            end if
            
	    if ( id_simfq > 0 ) then
                 rtmp = 0.
	         rtmp(:,:,1:nlev) = simfq
	         used = send_data ( id_simfq, rtmp, time, is, js,      &
		                    1, rmask=mask3 )
            end if
            
	    !----------------------------------------------
	    !
	    ! Total diffusivity coefficients
	    !
	    
	    if ( id_k_t_entr > 0 ) then
                 rtmp = 0.
	         rtmp(:,:,1:nlev) = k_t_entr
	         used = send_data ( id_k_t_entr, rtmp, time, is, js,   &
		                    1, rmask=mask3 )
            end if
            
	    if ( id_k_m_entr > 0 ) then
                 rtmp = 0.
	         rtmp(:,:,1:nlev) = k_m_entr
	         used = send_data ( id_k_m_entr, rtmp, time, is, js,   &
		                    1, rmask=mask3 )
            end if
            
       end if  ! do diagnostics if
       
!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine entrain

!
!======================================================================= 

!======================================================================= 
!
!  Subroutine to calculate pbl depth
!

subroutine pbl_depth(t, z, u_star, b_star, ipbl, h, parcelkick)

!
!  -----
!  INPUT
!  -----
!
!  t (= slv/cp)  liquid water virtual static energy divided by cp (K)
!  u_star        friction velocity (m/s)
!  b_star        buoyancy scale (m/s**2)
!       
!  ------
!  OUTPUT
!  ------
!
!  ipbl          half level containing pbl height
!  h             pbl height (m)
!  parcelkick    surface parcel excess (K)

real,    intent(in) ,  dimension(:) :: t, z
real,    intent(in)                 :: u_star, b_star
integer, intent(out)                :: ipbl
real,    intent(out)                :: h,parcelkick

real     :: svp,h1,h2,svp,t1,t2
real     :: ws,k_t_ref
integer  :: k,nlev

!initialize zsml
h = 0.

!compute # of levels
nlev = size(t,1)

!calculate surface parcel properties
svp  = t(nlev)
h1   = z(nlev)
call mo_diff(h1, u_star, b_star, ws, k_t_ref)
ws = max(small,ws/vonkarm/h1)
svp  = svp*(1.+(parcel_buoy*u_star*b_star/grav/ws) )
parcelkick = svp*parcel_buoy*u_star*b_star/grav/ws

!search for level where this is exceeded              
h    = h1
t1   = t(nlev)
do k = nlev-1 , 2, -1
     h2 = z(k)
     t2 = t(k)
     if (t2.gt.svp) then
          h = h2 + (h1 - h2)*(t2 - svp)/(t2 - t1 )
	  ipbl = k+1
          return
     end if
     h1 = h2
     t1 = t2
enddo

!one shouldn't end up here but nonetheless for safety this is put here
h = h2
ipbl = k+1

return

end subroutine pbl_depth

!=======================================================================

!======================================================================= 
!
!  Subroutine to do profile reconstuction
!

subroutine prof_recon(rho,t,pf,ph,zt,dt)

!
!  -----
!  INPUT
!  -----
!
!  rho    air density (kg/m3)
!  t      liquid water virtual static energy divided by cp (K)
!  pf     full level pressure (Pa)
!  ph     half level pressure (Pa)
!       
!  ------
!  OUTPUT
!  ------
!
!  zt     top of radiatively driven layer in distance relative to
!         boundary between cloud top layer and the level below (m)
!  dt     cloud top jump in liquid water virtual static energy divided 
!         by cp (K)
!
		 
real,   intent(in)                    :: rho
real,   intent(in) ,  dimension(-2:1) :: t, pf
real,   intent(in) ,  dimension( 0:1) :: ph
real,   intent(out)                   :: zt, dt

real, dimension(-2:1) :: pfp
real, dimension( 0:1) :: php

real                         :: svpar,slope,textrap
real                         :: a,b,c,det,pinv,ttop

!-----------------------------------------
! calculate all pressure relative to ph(0)
!
! pfp = full level relative pressures
! php = half level relative pressures
!
!  pfp(-2) < pfp(-1) < php(0) < pfp(0) < php(1) < pfp(1)
! 
! Note the following coordinate system
!
!            - - - - - - - - - - - - - -    
!
!                       *                   pfp(-2)
!
!            - - - - - - - - - - - - - - 
!
!                       *                   pfp(-1)
!
!            - - - - - - - - - - - - - -    php(0)
!
!  ambiguous layer ---> *                   pfp(0)
!
!            - - - - - - - - - - - - - -    php(1)
!
!                       *                   pfp(1)
!
!            - - - - - - - - - - - - - -
!

pfp = pf - ph(0)
php = ph - ph(0)

!----------------
! determine slope

slope = min ( (t(-2)-t(-1))/(pfp(-2)-pfp(-1)) , 0. )

! if this slope is such that the mean temperature of level 0 would
! exceed its actual temperature then assume a minimum protusion of
! mixed layer into ambiguous layer
!
! otherwise compute height of inversion using normal method


textrap = t(-1)+slope*(pfp(0)-pfp(-1))

if (textrap .lt. t(0)) then

     zt = 0.1*php(1)/rho/grav
     dt = t(0)-t(1)

else

     a = 0.5*slope
     b = t(-1)-t(1)-slope*pfp(-1)
     c = - php(1)*(t(0)-t(1))

     det = b*b - 4*a*c

     if (a.lt.0.) then
          pinv = (-b+sqrt(det))/(2*a)
     else
          pinv = c/b
     end if
     
     zt = (php(1)-pinv)/rho/grav
     ttop = t(-1) + slope*(pinv-pfp(-1))
     
     dt = ttop - t(1)
          		  
end if
 
return

end subroutine prof_recon

!=======================================================================

!======================================================================= 
!
!  Subroutine to calculate bottom and depth of radiatively driven mixed
!  layer
!
subroutine radml_depth(svp, zt, t, zf, zh, zb, zml)

!
!  -----
!  INPUT
!  -----
!
!  svp    cloud top liquid water virtual static energy divided by cp (K)
!  zt     top of radiatively driven layer (m)
!  t      liquid water virtual static energy divided by cp (K)
!  zf     full level height above ground (m)
!  zh     half level height above ground (m)
!       
!  ------
!  OUTPUT
!  ------
!
!  zb      base height of radiatively driven mixed layer (m)
!  zml     depth of radiatively driven mixed layer (m)


real,   intent(in)                 :: svp, zt
real,   intent(in) ,  dimension(:) :: t, zf, zh
real,   intent(out)                :: zb, zml

real    :: svpar,h1,h2,t1,t2
integer :: k,nlev

!initialize zml
zml = 0.

!compute # of levels
nlev = size(t,1)

!calculate cloud top parcel properties
svpar  = svp - radperturb
h1   = zf(1)
t1   = t(1)

!search for level where this is exceeded              
do k = 2,nlev
     h2 = zf(k)
     t2 = t(k)
     
     if (t2.lt.svpar) then
          zb = h2 + (h1 - h2)*(svpar - t2)/(t1 - t2)
	  zml = zt - zb
	  return
     end if
     
     if (do_jump_exit .and. (t1-t2) .gt. critjump .and. k .gt. 2) then
          zb = zh(k)
	  zml = zt - zb
	  return
     end if
     
     h1 = h2
     t1 = t2
enddo

zb = 0.
zml = zt
	  
return
end subroutine radml_depth

!=======================================================================


!=======================================================================

subroutine diffusivity_pbl(h, u_star, b_star, t, zm, k_m, k_t)
 	
real,    intent(in)                :: h, u_star, b_star
real,    intent(in),  dimension(:) :: t,zm
real,    intent(out), dimension(:) :: k_m, k_t

real    :: k_m_ref, k_t_ref, factor, hinner
integer :: k, kk, nlev

nlev = size(t,1)

k_m = 0.0
k_t = 0.0

hinner = frac_inner*h
kk = nlev+1
do k = 1, nlev
  if( zm(k) < hinner) then
      kk = k
      exit
  end if
end do

call mo_diff(hinner, u_star, b_star, k_m_ref, k_t_ref)

if (kk .lt. nlev+1) then 
     call mo_diff(zm(kk:nlev), u_star, b_star, k_m(kk:nlev),           &
                                               k_t(kk:nlev))
end if

if (kk .gt. 1) then
     do k = 1,kk-1
        factor = (zm(k)/hinner)* (1.0 -(zm(k)-hinner)/(h-hinner))**2
        k_m(k) = min( k_m_ref*factor, akmax )
        k_t(k) = min( k_t_ref*factor, akmax )
     end do
end if

return
end subroutine diffusivity_pbl

!
!======================================================================= 




!=======================================================================

subroutine diffusivity_sim(t, u, v, z, zz, k_m, k_t)

!-------------------------------
!
!  INPUT
!
!  defined on full levels 1:nlev
!
!  t      static energy of some type            (J/kg) 
!  u,v    horizontal wind speeds                (m/s)
!  z      full level height relative to surface (m)
!  
!  define on half levels 1:nlev+1
!
!  zz     half level height relative to surface (m)
!
!  OUTPUT 
!
!  defined on half levels 1:nlev
!
!  k_m    momentum diffusion coefficient        (m)
!  k_t    heat diffusion coefficient            (m)
!
!  k_m/k_t are defined so that k_m(k) is the diffusion coefficient
!  between full levels k-1 and k.
!

real, intent(in)    , dimension(:,:,:) :: t, u, v, z, zz
real, intent(out)   , dimension(:,:,:) :: k_m, k_t

real, dimension(size(t,1),size(t,2))   :: dz, b, speed2, rich
real, dimension(size(t,1),size(t,2))   :: fri, mix_len
integer                                :: k

do k = 2, size(t,3)

!----------------------------------------------------------------------
!  define the richardson number. set it to zero if it is negative. save
!  a copy of it for later use (rich2).
!----------------------------------------------------------------------
  
  dz     = z(:,:,k-1) - z(:,:,k)
  b      = grav*(t(:,:,k-1)-t(:,:,k))/t(:,:,k)
  speed2 = (u(:,:,k-1) - u(:,:,k))**2 + (v(:,:,k-1) - v(:,:,k))**2 
  rich= b*dz/(speed2+small)
  rich = max(rich, 0.0)

!-----------------------------------------------------------------------
! compute mixing length
!-----------------------------------------------------------------------

  mix_len(:,:) = (1./max(0.1,vonkarm*zz(:,:,k))) + (1./asympt_len)
  mix_len(:,:) = 1. / mix_len(:,:)
  
!-----------------------------------------------------------------------
!   compute the richardson number factor to be used in the eddy 
!   mixing coefficient. 
!-----------------------------------------------------------------------

  fri(:,:)   = 0.0
  where (rich .lt. rich_crit) 
       fri(:,:)   = (1.0 - rich/rich_crit)**2
  end where

!-----------------------------------------------------------------------
!   compute the eddy mixing coefficients. Values are obtained only when 
!   the richardson number is sub-critical
!---------------------------------------------------------------------

  where (rich < rich_crit)
     k_m(:,:,k) = mix_len*mix_len*sqrt(speed2)*fri(:,:)/dz
     k_t(:,:,k) = k_m(:,:,k)/pr
  end where
  
!---------------------------------------------------------------------
!   end loop over levels
!---------------------------------------------------------------------
  
end do

!---------------------------------------------------------------------
! bound diffusion coefficients
!---------------------------------------------------------------------

k_m = min ( max( 0., k_m ), akmax )
k_t = min ( max( 0., k_t ), akmax )

return
end subroutine diffusivity_sim

!
!======================================================================= 


!======================================================================= 
!
!      subroutine entrain_tend
!        
!
!      this subroutine takes the longwave heating rate and assigns it
!      to tdtlw
!        

subroutine entrain_tend(is,ie,js,je,tend)

!-----------------------------------------------------------------------
!
!      variables
!
!      -----
!      input
!      -----
!
!      is,ie,js,je       i,j indices marking the slab of model 
!      tend              longwave heating rate (deg K/sec)
!
!-----------------------------------------------------------------------

integer, intent(in)                   :: is,ie,js,je
real,    intent(in), dimension(:,:,:) :: tend

!-----------------------------------------------------------------------
!
!      assign tendency
!

       if (.not. entrain_on) return
       tdtlw(is:ie,js:je,:)=tend(:,:,:)

!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine entrain_tend

!
!======================================================================= 




!======================================================================= 
!
!      subroutine entrain_end
!        
!
!      this subroutine writes out the restart field
!        

subroutine entrain_end()

!-----------------------------------------------------------------------
!
!      variables
!
!      --------
!      internal
!      --------
!
!      unit              unit number for namelist and restart file
!
!-----------------------------------------------------------------------

integer :: unit

!-----------------------------------------------------------------------
!
!      write out restart file
!
       unit = Open_File ('RESTART/entrain.res', &
            FORM='native', ACTION='write')
       call write_data (unit, tdtlw)
       Call Close_File (unit)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 
!      subroutine end
!

end subroutine entrain_end

!
!=======================================================================

end module entrain_mod
