
                    module longwave_aerosol_mod

 
use  rad_utilities_mod,    only: longwave_control_type, Lw_control, &
                                  optical_path_type, &
                                 map_global_indices, &
	       aerosol_properties_type, Aerosol_props,  &
                                  longwave_parameter_type,   &
                                 aerosol_type, &
				  atmos_input_type,   &
                                 Lw_parameters, &
				 Environment, environment_type
use  utilities_mod,        only: open_file, file_exist,    &
			         check_nml_error, error_mesg, &
			         print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file
use constants_mod,         only: diffac
!use aerosol_mod,           only:                   &
!				aerosol_alloc
use longwave_params_mod,   only: NBLW

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                    longwave aerosol module
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

   character(len=128)  :: version =  '$Id: longwave_aerosol.F90,v 1.4 2003/04/09 20:59:42 fms Exp $'
   character(len=128)  :: tag     =  '$Name: inchon $'


!---------------------------------------------------------------------
!----- interfaces  -----
           
public         longwave_aerosol_init, &
               optical_depth_aerosol


!---------------------------------------------------------------------
!---- namelist   -----



logical           :: do_lwaerosol    = .false.
			  


namelist / longwave_aerosol_nml /    &
			          do_lwaerosol


!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


real, dimension (:,:), allocatable       :: aerextbandlw, aerssalbbandlw


integer, parameter           :: N_AEROSOL_BANDS_FR = 8
integer, parameter           :: N_AEROSOL_BANDS_CO = 1
integer, parameter           :: N_AEROSOL_BANDS = N_AEROSOL_BANDS_FR + N_AEROSOL_BANDS_CO
!---------------------------------------------------------------------
!---------------------------------------------------------------------


logical                      :: do_init = .false.


                           contains




subroutine longwave_aerosol_init


!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
     real, dimension (:,:,:), allocatable :: aerextivlstrlw,   &
		                             aerssalbivlstrlw,    &
			  	             aerasymmivlstrlw
    real, dimension (:,:), allocatable :: aerextivl, aerssalbivl
 real, dimension (:,:), allocatable       :: aerextbandlw_fr, aerssalbbandlw_fr
real, dimension (:,:), allocatable       :: aerextbandlw_co, aerssalbbandlw_co
    real, dimension (:,:), allocatable :: sflwwts, planckivlaer
    real, dimension (:,:), allocatable :: planckivlaer_fr, planckivlaer_co
    real, dimension (:,:), allocatable :: sflwwts_fr, sflwwts_co


     integer         :: unit, ierr, io

     integer         :: NAERIVLS, NAERMODELS
     integer, dimension(:), allocatable :: endsfbands, iendsfbands
    real            :: del, xtemv, sumplanck
     real, dimension(NBLW)           :: c1, centnb, sc, src1nb,  &
                                       x, x1
      integer         :: ib, nw, nivl, nband, n, ni
    logical         :: do_band1

!----------------------------------------------------------------------
!     parameters for longwave aerosol parameterizations
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   N_AEROSOL_BANDS = number of  infrared frequency bands over which 
!                   frequency-dependent aerosol optical parameters
!                   using appropriate parameterizations.
!----------------------------------------------------------------------


 integer, dimension (N_AEROSOL_BANDS) :: nivl1aer, nivl2aer
integer, dimension (N_AEROSOL_BANDS_FR) :: nivl1aer_fr, nivl2aer_fr
integer, dimension (N_AEROSOL_BANDS_CO) :: nivl1aer_co, nivl2aer_co
 real,    dimension (N_AEROSOL_BANDS) :: planckaerband


!----------------------------------------------------------------------
!   wavenumber ranges with separate aerosol optical properties in the
!   infrared parameterization. these may be changed only by the
!   keeper of the radiation code.
!----------------------------------------------------------------------

!    the order of the frequency bands corresponds to the order used
!    in the lw radiation code
 real, dimension (N_AEROSOL_BANDS_FR)        :: aerbandlo_fr =  &
      (/ 560.0, 630.0, 700.0, 800.0, 900.0,  990.0, 1070.0, 1200.0 /)
 real, dimension (N_AEROSOL_BANDS_FR)        :: aerbandhi_fr =  &
      (/ 630.0, 700.0, 800.0, 900.0, 990.0, 1070.0, 1200.0, 1400.0 /)
 integer, dimension (N_AEROSOL_BANDS_FR)        :: istartaerband_fr =  &
      (/ 57,  64,  71,  81,  91, 100, 108, 121 /)
 integer, dimension (N_AEROSOL_BANDS_FR)        :: iendaerband_fr =  &
      (/ 63,  70,  80,  90,  99, 107, 120, 140 /)
 real, dimension (N_AEROSOL_BANDS_CO)        :: aerbandlo_co =  &
      (/ 560.0 /)
 real, dimension (N_AEROSOL_BANDS_CO)        :: aerbandhi_co =  &
      (/ 800.0 /)
 integer, dimension (N_AEROSOL_BANDS_CO)     :: istartaerband_co =  &
     (/ 57  /)
 integer, dimension (N_AEROSOL_BANDS_CO)     :: iendaerband_co =  &
      (/ 80  /)
 real, dimension(N_AEROSOL_BANDS)       :: aerbandlo, aerbandhi
 integer, dimension(N_AEROSOL_BANDS)    :: istartaerband, iendaerband


          do n=1,N_AEROSOL_BANDS_FR
         aerbandlo(n) = aerbandlo_fr(n)
        aerbandhi(n) = aerbandhi_fr(n)
         istartaerband(n) = istartaerband_fr(n)
         iendaerband(n) = iendaerband_fr(n)
      enddo
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
         aerbandlo(n) = aerbandlo_co(n-N_AEROSOL_BANDS_FR)
        aerbandhi(n) = aerbandhi_co(n-N_AEROSOL_BANDS_FR)
         istartaerband(n) = istartaerband_co(n-N_AEROSOL_BANDS_FR)
         iendaerband(n) = iendaerband_co(n-N_AEROSOL_BANDS_FR)
      enddo

!---------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=longwave_aerosol_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'longwave_aerosol_nml')
       enddo
10     call close_file (unit)
     endif

     unit = open_file ('logfile.out', action='append')
!    call print_version_number (unit, 'longwave_aerosol',   &
!						 version_number)
     if (get_my_pe() == 0) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=longwave_aerosol_nml)
     endif
     call close_file (unit)

     Lw_control%do_lwaerosol = do_lwaerosol

     Lw_parameters%n_lwaerosol_bands = N_AEROSOL_BANDS

!-------------------------------------------------------------------
!   ****** the remainder of this subroutine need be executed only if 
!   ****** do_lwaerosol is true.
!-------------------------------------------------------------------
 

!-------------------------------------------------------------------
!   if longwave aerosols are activated, use lwaerosol_form to define 
!   whether or not predicted lwaerosols are to be used
!-------------------------------------------------------------------
     if (do_lwaerosol) then


!--------------------------------------------------------------------
!     compute band-averaged coefficients for aerosols in infrared
!     frequency ranges.
!     the actual extinction coefficients (to become emissivities)
!      (and other coefficients) are calculated as time-dependent
!     quantities.
!
!--------------------------------------------------------------------
 
     NAERIVLS = Aerosol_props%num_wavenumbers
     allocate (endsfbands(NAERIVLS) )
     allocate (iendsfbands(NAERIVLS) )
     endsfbands = Aerosol_props%endaerwvnsf
     NAERMODELS = Aerosol_props%NAERMODELS

     allocate (aerextivl    (NAERIVLS, NAERMODELS) )
     allocate (aerssalbivl  (NAERIVLS, NAERMODELS) )

     allocate (aerextbandlw   (N_AEROSOL_BANDS, NAERMODELS) )
     allocate (aerextbandlw_fr   (N_AEROSOL_BANDS_FR, NAERMODELS) )
     allocate (aerextbandlw_co   (N_AEROSOL_BANDS_CO, NAERMODELS) )
     allocate (aerssalbbandlw (N_AEROSOL_BANDS, NAERMODELS) )
     allocate (aerssalbbandlw_fr (N_AEROSOL_BANDS_FR, NAERMODELS) )
     allocate (aerssalbbandlw_co (N_AEROSOL_BANDS_CO, NAERMODELS) )
 
     allocate (planckivlaer (N_AEROSOL_BANDS, NAERIVLS) )
     allocate (planckivlaer_fr (N_AEROSOL_BANDS_FR, NAERIVLS) )
     allocate (planckivlaer_co (N_AEROSOL_BANDS_CO, NAERIVLS) )
     allocate (sflwwts (N_AEROSOL_BANDS, NAERIVLS) )
     allocate (sflwwts_fr (N_AEROSOL_BANDS_FR, NAERIVLS) )
     allocate (sflwwts_co (N_AEROSOL_BANDS_CO, NAERIVLS) )
!-------------------------------------------------------------------
!   get aerosol optical property data from aerosol module
!-------------------------------------------------------------------
    aerextivl = Aerosol_props%aeroextivl
    aerssalbivl = Aerosol_props%aerossalbivl
 

!--------------------------------------------------------------------
!    calculation for aerosols (Shettle-Fenn wavelength intervals)
!--------------------------------------------------------------------
 
        iendsfbands(:) = INT((endsfbands(:)+0.01)/10.0)

!--------------------------------------------------------------------
!      compute weighting function.  For aerosols, I am using
!      the Planck function at 10C.
!--------------------------------------------------------------------
        do n=1,NBLW 
          del  = 10.0E+00
	  xtemv = 283.15
	  centnb(n) = 5.0 + (n - 1)*del
          c1(n)     = (3.7412E-05)*centnb(n)**3
          x(n)      = 1.4387E+00*centnb(n)/xtemv
          x1(n)     = EXP(x(n))
          sc(n)     = c1(n)/(x1(n) - 1.0E+00)
          src1nb(n) = del*sc(n)
        enddo
 
!--------------------------------------------------------------------
!      compute summed weighting function over the (N_AEROSOL_BANDS) 
!      aerosol bands
!--------------------------------------------------------------------

        planckaerband(:) = 0.0E+00
        do n = 1,N_AEROSOL_BANDS
  	  do ib = istartaerband(n),iendaerband(n)
	    planckaerband(n) = planckaerband(n) + src1nb(ib)
	  enddo
        enddo

 
        nivl = 1
        sumplanck = 0.0
        nband = 1
        planckivlaer_fr(:,:) = 0.0
        nivl1aer_fr(1) = 1
	do_band1 = .true.
 
do nw = 1,NBLW
  	  sumplanck = sumplanck + src1nb(nw)
          if ( nw.eq.iendsfbands(nivl) ) then
            planckivlaer_fr(nband,nivl) = sumplanck
            sumplanck = 0.0
          end if
 if ( nw.eq.iendaerband_fr(nband) ) then
            if ( nw.ne.iendsfbands(nivl) ) then
              planckivlaer_fr(nband,nivl) = sumplanck 
              sumplanck = 0.0
            end if
            nivl2aer_fr(nband) = nivl
            nband = nband + 1
            if ( nband.le.N_AEROSOL_BANDS_FR ) then
              if ( nw.eq.iendsfbands(nivl) ) then
                nivl1aer_fr(nband) = nivl + 1
              else
                nivl1aer_fr(nband) = nivl
              end if
            end if
 end if
 if ( nw .eq. iendsfbands(nivl) ) then
	    nivl = nivl + 1
	    if (do_band1 .and. nband .eq. 1 .and. iendsfbands(nivl-1) .ge. istartaerband_fr(1) .and.  &
	        iendsfbands(nivl-1) .lt. iendaerband_fr(1)) then
	       nivl1aer_fr(nband) = nivl-1
	       do_band1 = .false.
            endif
 endif
 if ( nw .ge. iendaerband_fr(N_AEROSOL_BANDS_FR) ) then
	    exit
 endif
end do
!--------------------------------------------------------------------
!     compute planck-weighted band weights for Shettle-Fenn aerosol
!     calculations
!--------------------------------------------------------------------

        nivl = 1
        sumplanck = 0.0
        nband = 1
        planckivlaer_co(:,:) = 0.0
        nivl1aer_co(1) = 1
	do_band1 = .true.
 
        do nw = 1,NBLW
  	  sumplanck = sumplanck + src1nb(nw)
          if ( nw.eq.iendsfbands(nivl) ) then
            planckivlaer_co(nband,nivl) = sumplanck
            sumplanck = 0.0
          end if
          if ( nw.eq.iendaerband_co(nband) ) then
            if ( nw.ne.iendsfbands(nivl) ) then
              planckivlaer_co(nband,nivl) = sumplanck 
              sumplanck = 0.0
            end if
            nivl2aer_co(nband) = nivl
            nband = nband + 1
            if ( nband.le.N_AEROSOL_BANDS_CO ) then
              if ( nw.eq.iendsfbands(nivl) ) then
                nivl1aer_co(nband) = nivl + 1
              else
                nivl1aer_co(nband) = nivl
              end if
            end if
          end if
          if ( nw .eq. iendsfbands(nivl) ) then
	    nivl = nivl + 1
	    if (do_band1 .and. nband .eq. 1 .and. iendsfbands(nivl-1) .ge. istartaerband_co(1) .and.  &
	        iendsfbands(nivl-1) .lt. iendaerband_co(1)) then
	       nivl1aer_co(nband) = nivl-1
	       do_band1 = .false.
            endif
	  endif
          if ( nw .ge. iendaerband_co(N_AEROSOL_BANDS_CO) ) then
	    exit
	  endif
        end do
!--------------------------------------------------------------------
!     compute planck-weighted band weights for Shettle-Fenn aerosol
!     calculations
!--------------------------------------------------------------------

        sflwwts_fr(:,:) = 0.0E+00
        do n=1,N_AEROSOL_BANDS_FR
          do ni=nivl1aer_fr(n),nivl2aer_fr(n)
	    sflwwts_fr(n,ni) = planckivlaer_fr(n,ni)/planckaerband(n)
	  enddo
        enddo
        sflwwts_co(:,:) = 0.0E+00
        do n=1,N_AEROSOL_BANDS_CO
          do ni=nivl1aer_co(n),nivl2aer_co(n)
	    sflwwts_co(n,ni) = planckivlaer_co(n,ni)/planckaerband(N_AEROSOL_BANDS_FR+n)
	  enddo
        enddo

      do n=1,N_AEROSOL_BANDS_FR
        do ni = 1,NAERIVLS
	  sflwwts(n,ni) = sflwwts_fr(n,ni)
        enddo
      enddo
      do n=N_AEROSOL_BANDS_FR+1,N_AEROSOL_BANDS
        do ni = 1,NAERIVLS
	  sflwwts(n,ni) = sflwwts_co(n-N_AEROSOL_BANDS_FR,ni)
        enddo
      enddo
!-----------------------------------------------------------------------
!    use band weighting factors 
!    to derive band quantities for aerosol optical properties
!-----------------------------------------------------------------------
      aerextbandlw = 0.0E+00
      aerssalbbandlw = 0.0E+00
    do nw = 1, NAERMODELS     !   loop on aerosol types
      do n = 1,N_AEROSOL_BANDS  ! loop on radiation code aerosol freq bands
	do ni = 1,NAERIVLS ! loop on aerosol data freq bands (should be a function of aerosol type)
	  aerextbandlw(n,nw) = aerextbandlw(n,nw) +      &
               aerextivl(ni,nw)*sflwwts(n,ni)
 	  aerssalbbandlw(n,nw) = aerssalbbandlw(n,nw) +  &
               aerssalbivl(ni,nw)*sflwwts(n,ni)
	enddo
      enddo
    enddo
!---------------------------------------------------------------------


      endif
!-----------------------------------------------------------------
!  set flag to indicate that the routine has been completed.
!----------------------------------------------------------------

      do_init = .true.



end subroutine longwave_aerosol_init


!#####################################################################

subroutine optical_depth_aerosol (Atmos_input, n, naerosol,   &
                                  Aerosol,                   Optical)

type(atmos_input_type), intent(in) :: Atmos_input
integer,                  intent(in)  :: n, naerosol
type(aerosol_type), intent(in)  :: Aerosol
type(optical_path_type), intent(inout) :: Optical


integer                               :: i,j,k
integer   :: ix, jx, kx
integer                               :: nsc, opt_index
real, dimension( size(Aerosol%aerosol,1),  &
                 size(Aerosol%aerosol,2),  &
                 size(Aerosol%aerosol,3), &
                       Aerosol%nfields)   :: aerooptdepspec
real, dimension( size(Aerosol%aerosol,1),  &
                 size(Aerosol%aerosol,2),  &
                 size(Aerosol%aerosol,3))  :: aerooptdep
integer, dimension( size(Aerosol%aerosol,1),  &
                 size(Aerosol%aerosol,2),  &
                 size(Aerosol%aerosol,3))  :: irh, opt_index_v
 real :: asum
 real, dimension(size (Aerosol%aerosol,3)+1) :: bsum

    ix = size (Aerosol%aerosol,1)
    jx = size (Aerosol%aerosol,2)
    kx = size (Aerosol%aerosol,3)

   aerooptdep(:,:,:) = 0.0
   Optical%totaerooptdep(:,:,:,n) = 0.0
!
   do k = 1,kx         
   do j = 1,jx         
   do i = 1,ix           
    irh(i,j,k) = MIN(100, MAX(0, NINT(100.*Atmos_input%rel_hum(i,j,k))))
   opt_index_v(i,j,k) = Aerosol%sulfate_index( irh(i,j,k) )
   enddo
   enddo
   enddo



   do k = 1,kx         
   do j = 1,jx         
   do i = 1,ix           

     opt_index = opt_index_v(i,j,k)
      do nsc = 1,Aerosol%nfields  ! loop on aerosol species

!  using relative humidity criterion (where necessary) determine the
! aerosol category (as an index) appropriate for the aerosol species

 if( Aerosol%optical_index(nsc) == 0 ) then
! ... Sulfate aerosol (RH-dependent)



        aerooptdepspec(i,j,k,nsc) = diffac*Aerosol%aerosol(i,j,k,nsc)*   &
	  (1.0 - aerssalbbandlw(n,opt_index))*aerextbandlw(n,opt_index)

    endif
      enddo
   enddo
   enddo
   enddo

   do k = 1,kx         
   do j = 1,jx         
   do i = 1,ix           

      do nsc = 1,Aerosol%nfields  ! loop on aerosol species

!  using relative humidity criterion (where necessary) determine the
! aerosol category (as an index) appropriate for the aerosol species

  if( Aerosol%optical_index(nsc) == 0 ) then
  else
! ... Other aerosols
    opt_index = Aerosol%optical_index(nsc)

 if( opt_index == 0 ) &
   call error_mesg( 'get_aerosol_optical_index', &
                  'Cannot find aerosol optical properties for species = ' // &
!                   TRIM( Aerosol%data_names(nsc) ), &
                    TRIM( Aerosol_props%aerosol_names(nsc) ), &
               FATAL )


        aerooptdepspec(i,j,k,nsc) = diffac*Aerosol%aerosol(i,j,k,nsc)*   &
	  (1.0 - aerssalbbandlw(n,opt_index))*aerextbandlw(n,opt_index)

    end if
      enddo
   enddo
   enddo
   enddo

!   sum optical depths over all species and obtain column optical depth
   do k=1,kx
   do j=1,jx
   do i=1,ix
     asum = 0.0
   do nsc = 1,Aerosol%nfields
   asum = asum + aerooptdepspec(i,j,k,nsc)
   enddo
   aerooptdep(i,j,k) = asum                                         
   enddo
   enddo
   enddo

   do j=1,jx
   do i=1,ix
     bsum(1) = 0.0
   do k = 2,kx+1         
      bsum(k) = bsum(k-1) + aerooptdep(i,j,k-1)
   enddo
   do k = 2,kx+1         
          Optical%totaerooptdep(i,j,k,n) = bsum(k)
   enddo

   enddo
   enddo

! continuum band is the last indx:
    if ( n == N_AEROSOL_BANDS) then
    Optical%aerooptdep_KE_15(:,:) = aerooptdep(:,:,kx)
    endif
    

end subroutine optical_depth_aerosol


!#####################################################################


	end module longwave_aerosol_mod


