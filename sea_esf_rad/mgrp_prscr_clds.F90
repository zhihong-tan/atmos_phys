
                 module mgrp_prscr_clds_mod

use utilities_mod,          only: open_file, file_exist,   &
                                  check_nml_error, error_mesg,   &
                                  print_version_number, FATAL, NOTE, &
				  WARNING, get_my_pe, close_file, &
				  get_domain_decomp
!use rad_step_setup_mod,     only: jabs, ISRAD, IERAD, JSRAD, JERAD, & 
!                                  KSRAD, KERAD, get_std_pressures
!use rad_step_setup_mod,     only: jabs,                             & 
!use rad_step_setup_mod,     only:   get_std_pressures
!use longwave_setup_mod,     only: longwave_parameter_type, &    
!                                  Lw_parameters
!use std_pressures_mod,      only: get_std_pressures
use constants_mod,          only: pstd, radian
use rad_utilities_mod,      only: shortwave_control_type, Sw_control, &
                                  map_global_indices, &
				  cld_diagnostics_type, &
                                  longwave_parameter_type, &    
                                  Lw_parameters, &
				  longwave_control_type, Lw_control 
use microphys_rad_mod,      only: microphys_rad_init, &
                                  microphys_rad_driver,&
                                  microphys_presc_conc,   &
                                  lwemiss_calc


!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!	       mgroup prescribed cloud properties module
!               (this module runnable in SKYHI and FMS; 
!                zonal_clouds_mod is FMS native equivalent)
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

! character(len=5), parameter  ::  version_number = 'v0.09'
  character(len=128)  :: version =  '$Id: mgrp_prscr_clds.F90,v 1.4 2003/04/09 21:00:29 fms Exp $'
  character(len=128)  :: tag     =  '$Name: inchon $'



!---------------------------------------------------------------------
!-------  interfaces --------

public     mgrp_prscr_init, mgrp_prscr_calc


private    cldht, cldint  


!---------------------------------------------------------------------
!-------- namelist  ---------



integer    ::     dummy=0



namelist /mgrp_prscr_clds_nml /     &
                                           dummy


!----------------------------------------------------------------------
!----  public data -------

!----------------------------------------------------------------------
!----  private data -------


!--------------------------------------------------------------------
!     these variables define cloud amounts and radiative properties  
!     on a global (j,k) model grid.
!--------------------------------------------------------------------
  
real, dimension(:,:),   allocatable :: zcamtmxo, zcamtrnd, zcrfvis,  &
			  	       zcrfir, zcabir
real, dimension(:,:,:), allocatable :: zemmxo, zemrnd
real, dimension(:),     allocatable :: cldhm_abs, cldml_abs

 
!--------------------------------------------------------------------
!     these variables define cloud tops and bottoms, amounts and rad-
!     iative properties on an input (LATOBS,NOFCLDS_SP) grid. the input 
!     values are the original Skyhi values. 
!--------------------------------------------------------------------
 
integer, parameter                    :: NOFCLDS_SP=3  
integer, parameter                    :: NOFMXOLW=1  
integer, parameter                    :: NOFRNDLW=2  
integer, parameter                    :: LATOBS=19

integer, dimension(:,:), allocatable  :: kkbh, kkth
real,    dimension(NOFCLDS_SP)        :: crfvis, crfir, cabir, cldem

!-------------------------------------------------------------------
!   default low, middle, high cloud properties (index = 1, 2, 3)
!   cldem    : infrared emissivity
!   crfvis   : visible band reflectivity
!   crfir    : near-ir band reflectivity
!   cabir    : near-ir band absorptivity
!-------------------------------------------------------------------

data cldem  / 1.00E+00, 1.00E+00, 1.00E+00 / 
data crfvis / 0.69E+00, 0.48E+00, 0.21E+00 / 
data crfir  / 0.69E+00, 0.48E+00, 0.21E+00 / 
data cabir  / 0.35E-01, 0.02E+00, 0.05E-01 /

!-------------------------------------------------------------------
!   prescribed high, mid, low cloud amounts on 19 latitudes
!-------------------------------------------------------------------

integer                                 ::  jj, kkc
real, dimension(LATOBS,NOFCLDS_SP)      ::  ccd           
real, dimension(LATOBS)                 ::  cloud_lats

data cloud_lats / -90., -80., -70., -60., -50., -40., -30., -20., &
                  -10., 0.0, 10., 20., 30., 40., 50., 60., 70., 80., &
                   90. /

data ((ccd(jj,kkc),jj=1,19),kkc=1,3)/    &
     &  0.360E+00, 0.401E+00, 0.439E+00, 0.447E+00, 0.417E+00, &
     &  0.343E+00, 0.269E+00, 0.249E+00, 0.290E+00, 0.330E+00, &
     &  0.290E+00, 0.249E+00, 0.269E+00, 0.343E+00, 0.417E+00, &
     &  0.447E+00, 0.439E+00, 0.401E+00, 0.360E+00,            &
     &  0.090E+00, 0.102E+00, 0.117E+00, 0.128E+00, 0.122E+00, &
     &  0.095E+00, 0.070E+00, 0.060E+00, 0.068E+00, 0.080E+00, &
     &  0.068E+00, 0.060E+00, 0.070E+00, 0.095E+00, 0.122E+00, &
     &  0.128E+00, 0.117E+00, 0.102E+00, 0.090E+00,            &
     &  0.198E+00, 0.231E+00, 0.254E+00, 0.250E+00, 0.227E+00, &
     &  0.192E+00, 0.159E+00, 0.168E+00, 0.205E+00, 0.241E+00, &
     &  0.205E+00, 0.168E+00, 0.159E+00, 0.192E+00, 0.227E+00, &
     &  0.250E+00, 0.254E+00, 0.231E+00, 0.198E+00/ 
!----------------------------------------------------------------------


real, dimension(:), allocatable :: qlevel
  
!----------------------------------------------------------------------
!    NLWCLDB is the number of frequency bands for which lw
!    emissitivies are defined.
!----------------------------------------------------------------------
!integer            :: NLWCLDB 

real               :: psbar  
!integer            :: x(4), y(4), jdf, jd
!character(len=10)  :: swform
!character(len=16)  :: swform
logical    :: do_lhsw, do_esfsw
logical            :: do_lwcldemiss
integer            :: ksrad, kerad


!----------------------------------------------------------------------
!----------------------------------------------------------------------




	 contains 


!subroutine mgrp_prscr_init (kx, pd, latb, qlyr)
subroutine mgrp_prscr_init (kx, pref, latb, qlyr)

!------------------------------------------------------------------
integer, intent(in) :: kx
!real, dimension(:), intent(in)             :: pd, latb      
real, dimension(:), intent(in)             ::  latb      
real, dimension(:,:), intent(in)             :: pref          
real, dimension(:), intent(in), optional   :: qlyr
!--------------------------------------------------------------------

      integer            :: unit, ierr, io
      integer            :: j, k, kc, ktop, kbot, li
      integer            :: jdf
      real               :: fl
         

!     real, dimension(:),     allocatable :: pd, cldhm, cldml, &
      real, dimension(:),     allocatable ::     cldhm, cldml, &
				             cldhm_abs_gl, cldml_abs_gl
      real, dimension(:,:),   allocatable :: emmxo19, emrnd19,    &
					     crfvis19, crfir19,   &
					     cabir19, camtmxo19,&
                                             camtrnd19, zcamtmxo_g,  &
					     zcamtrnd_g, zcrfvis_g,  &
	  				     zcrfir_g, zcabir_g
      real, dimension(:,:,:), allocatable :: zemmxo_g, zemrnd_g
      integer, dimension(:), allocatable :: jindx2

!---------------------------------------------------------------------
!-----  read namelist  ------
  
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=mgrp_prscr_clds_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'mgrp_prscr_clds_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'mgrp_prscr_clds',  &
!						version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=mgrp_prscr_clds_nml)
      endif
      call close_file (unit)


!--------------------------------------------------------------------
!  retrieve module variables that come from other modules
!--------------------------------------------------------------------

      ksrad = 1
      kerad = kx
!     call get_domain_decomp (x, y)
!     jd = y(2)
!     jdf = y(4) - y(3) + 1
      jdf = size(latb,1) - 1

!      if (Environment%running_gcm) then
!     if (trim (swaerosol_form) /= ' ' .and. jmax_aerfile /= 0) then
         allocate (jindx2  (size(latb,1)-1))
         call find_nearest_index (latb, jindx2)
!        print *, 'jindx2', get_my_pe(), jindx2 

!---------------------------------------------------------------------
!    save the sigma levels that have been input for later use within 
!    this module.
!---------------------------------------------------------------------
      allocate (qlevel (KSRAD:KERAD) )
      if (present(qlyr)) then
        qlevel(:) = qlyr(KSRAD:KERAD)
      else
!       allocate ( pd(KSRAD:KERAD) )
!       call get_std_pressures (pd_out=pd)
        psbar = pstd*1.0E-03   !  convert from cgs to mb
!qlevel(KSRAD:KERAD) = pd(KSRAD:KERAD)/psbar
!qlevel(KSRAD:KERAD) = (1.0E-02*pd(KSRAD:KERAD))/psbar
	qlevel(KSRAD:KERAD) = (1.0E-02*pref(KSRAD:KERAD,1))/psbar
!deallocate (pd)
      endif

!---------------------------------------------------------------------
!    define the number of cloud emissivity bands for use in this module.
!---------------------------------------------------------------------

!     NLWCLDB = Lw_parameters%NLWCLDB
!      swform = Sw_control%sw_form
      do_lhsw = Sw_control%do_lhsw
      do_esfsw = Sw_control%do_esfsw
      do_lwcldemiss = Lw_control%do_lwcldemiss

!---------------------------------------------------------------------
!   allocate space to hold the cloud radiative properties on the model
!   grid (j,k)
!---------------------------------------------------------------------

      allocate( zcamtmxo    (jdf , KSRAD:KERAD))
      allocate( zcamtrnd    (jdf , KSRAD:KERAD))
!     allocate( zemmxo      (jdf , KSRAD:KERAD, NLWCLDB))
!     allocate( zemrnd      (jdf , KSRAD:KERAD, NLWCLDB))
      allocate( zemmxo      (jdf , KSRAD:KERAD, 1      ))
      allocate( zemrnd      (jdf , KSRAD:KERAD, 1      ))
      allocate( zcrfvis     (jdf , KSRAD:KERAD))
      allocate( zcrfir      (jdf , KSRAD:KERAD))
      allocate( zcabir      (jdf , KSRAD:KERAD))
!     allocate( zcamtmxo_g  (jd , KSRAD:KERAD))
!     allocate( zcamtrnd_g  (jd , KSRAD:KERAD))
!     allocate( zemmxo_g    (jd , KSRAD:KERAD, NLWCLDB))
!     allocate( zemrnd_g    (jd , KSRAD:KERAD, NLWCLDB))
!     allocate( zcrfvis_g   (jd , KSRAD:KERAD))
!     allocate( zcrfir_g    (jd , KSRAD:KERAD))
!     allocate( zcabir_g    (jd , KSRAD:KERAD))

!--------------------------------------------------------------------
!   allocate arrays to hold the cloud radiative properties at the 
!   input latitudes and model  radiation levels.
!--------------------------------------------------------------------

      allocate ( emmxo19   (LATOBS, KSRAD:KERAD))
      allocate ( emrnd19   (LATOBS, KSRAD:KERAD) ) 
      allocate ( crfvis19  (LATOBS, KSRAD:KERAD) )
      allocate ( crfir19   (LATOBS, KSRAD:KERAD)) 
      allocate ( cabir19   (LATOBS, KSRAD:KERAD))
      allocate ( camtmxo19 (LATOBS, KSRAD:KERAD))
      allocate ( camtrnd19 (LATOBS, KSRAD:KERAD) )

!---------------------------------------------------------------------
!     define index arrays for cloud tops and cloud bottoms(kkth, kkbh). 
!     a program courtesy of Ron Stouffer is used. cldht specifies cloud 
!     height data in mb for each cloud type, at 10 deg intervals from 
!     90S to 90N.(ie, 19 lats).
!---------------------------------------------------------------------
 
      allocate ( kkbh(LATOBS, NOFCLDS_SP) )
      allocate ( kkth(LATOBS, NOFCLDS_SP) )
      call cldht 

!-------------------------------------------------------------------
!    compute cloud quantities at the 19 canonical latitudes. set values
!    to 0.0 at points where clouds are not prescribed.
!-------------------------------------------------------------------

      do j=1,19
	do k=KSRAD,KERAD
	  camtrnd19(j,k) = 0.0E+00
	  camtmxo19(j,k) = 0.0E+00
          emrnd19(j,k)   = 0.0E+00
	  emmxo19(j,k)   = 0.0E+00
  	  crfvis19(j,k)  = 0.0E+00
	  crfir19(j,k)   = 0.0E+00
	  cabir19(j,k)   = 0.0E+00
	enddo
        do kc = 1,NOFCLDS_SP
          ktop = kkth(j,kc)
	  kbot = kkbh(j,kc)
	  do k=ktop,kbot
	    if (ktop .NE. kbot) then
	      camtmxo19(j,k) = ccd(j,kc)
	      emmxo19(j,k)   = cldem(kc)
            else
	      camtrnd19(j,k) = ccd(j,kc)
	      emrnd19(j,k)   = cldem(kc)
            endif
	    crfvis19(j,k) = crfvis(kc)
	    crfir19(j,k) = crfir(kc)
	    cabir19(j,k) = cabir(kc)
	  enddo
	enddo
      enddo

!---------------------------------------------------------------------
!    perform latitude "interpolation" to the (JD) model latitudes
!    in default case, no interpolation is actually done; the nearest
!    latitude available (using NINT function) is used.
!---------------------------------------------------------------------

!     do j=1,jd
      do j=1,jdf
!       fl = 9.0E+00 - 9.0E+00*(FLOAT(jd+1 -j - jd/2) -    &
!            0.5E+00)/FLOAT(jd/2)
!       li = NINT(fl) + 1
        li = jindx2(j)
	do k=KSRAD,KERAD
!  zcamtmxo_g(j,k) = camtmxo19(li,k)
!  zcamtrnd_g(j,k) = camtrnd19(li,k)
!  zemmxo_g(j,k,:) = emmxo19(li,k)
!  zemrnd_g(j,k,:) = emrnd19(li,k)
!  zcrfvis_g(j,k)  = crfvis19(li,k)
!  zcrfir_g(j,k)   = crfir19(li,k)
!  zcabir_g(j,k)   = cabir19(li,k)
	  zcamtmxo(j,k) = camtmxo19(li,k)
	  zcamtrnd(j,k) = camtrnd19(li,k)
!  zemmxo(j,k,:) = emmxo19(li,k)
!  zemrnd(j,k,:) = emrnd19(li,k)
	  zemmxo(j,k,1) = emmxo19(li,k)
	  zemrnd(j,k,1) = emrnd19(li,k)
	  zcrfvis(j,k)  = crfvis19(li,k)
	  zcrfir(j,k)   = crfir19(li,k)
	  zcabir(j,k)   = cabir19(li,k)
        end do
      end do
!     do j=1,jdf
!do k=KSRAD,KERAD
!         zcamtmxo(j,k) = zcamtmxo_g(j+y(3)-1,k)	
!         zcamtrnd(j,k) = zcamtrnd_g(j+y(3)-1,k)	
!         zemmxo(j,k,:) = zemmxo_g(j+y(3)-1,k,:)	
!         zemrnd(j,k,:) = zemrnd_g(j+y(3)-1,k,:)	
!         zcrfvis (j,k) = zcrfvis_g(j+y(3)-1,k)	
!         zcrfir  (j,k) = zcrfir_g (j+y(3)-1,k)	
!         zcabir  (j,k) = zcabir_g (j+y(3)-1,k)	
!       end do
!     end do
!     deallocate (zcamtmxo_g)
!     deallocate (zcamtrnd_g)
!     deallocate (zemmxo_g )
!     deallocate (zemrnd_g )
!     deallocate (zcrfvis_g)
!     deallocate (zcrfir_g )
!     deallocate (zcabir_g )

      deallocate ( kkbh)
      deallocate ( kkth)
      deallocate ( emmxo19 )
      deallocate ( emrnd19  ) 
      deallocate ( crfvis19 )
      deallocate ( crfir19 ) 
      deallocate ( cabir19 )
      deallocate ( camtmxo19 )
      deallocate ( camtrnd19  )

!--------------------------------------------------------------------
!   allocate arrays for the sigmas of the low-mid and mid-high cloud
!   interfaces
!---------------------------------------------------------------------

      allocate ( cldhm_abs(jdf) )
      allocate ( cldml_abs(jdf) )
      allocate ( cldhm    (LATOBS) )
      allocate ( cldml    (LATOBS) )
      
!---------------------------------------------------------------------
!      this block supplies the standard skyhi values for these param-
!      eters, to be used for checkout when needed.
!----------------------------------------------------------------------
      do jj=1,LATOBS
        cldhm(jj) = 0.52
        cldml(jj) = 0.70
      end do
      cldml(1) = 0.78
      cldml(2) = 0.78
      cldml(18) = 0.78
      cldml(19) = 0.78

!---------------------------------------------------------------------
!    perform latitude "interpolation" to the (JD) model latitudes
!    in default case, no interpolation is actually done; the nearest
!    latitude available (using NINT function) is used.
!---------------------------------------------------------------------
  
!     allocate ( cldhm_abs_gl (jd) )
!     allocate ( cldml_abs_gl (jd) )
!     do j=1,jd
      do j=1,jdf
!         fl = 9.0E+00 - 9.0E+00*(FLOAT(JD+1 -j - JD/2) -    &
!              0.5E+00)/FLOAT(JD/2)
!         li = NINT(fl) + 1
          li = jindx2(j)
!         cldhm_abs_gl(j) = cldhm(li)
!         cldml_abs_gl(j) = cldml(li)
          cldhm_abs   (j) = cldhm(li)
          cldml_abs   (j) = cldml(li)
      end do
!     do j=1,jdf
!         cldhm_abs(j) =  cldhm_abs_gl(j+y(3)-1)
!         cldml_abs(j) =  cldml_abs_gl(j+y(3)-1)
!     end do
!     deallocate (cldhm_abs_gl)
!     deallocate (cldml_abs_gl)

!---------------------------------------------------------------------
!    if a cloud microphysics scheme is to be employed with the cloud
!    scheme, initialize the microphysics_rad module.
!--------------------------------------------------------------------

!     if (trim(swform) == 'esfsw99' .or. do_lwcldemiss) then
      if (do_esfsw                  .or. do_lwcldemiss) then
        call microphys_rad_init (cldhm_abs, cldml_abs)
      endif

!---------------------------------------------------------------------
!    deallocate local variables
!---------------------------------------------------------------------
      deallocate (cldhm)
      deallocate (cldml)


end subroutine mgrp_prscr_init


!#####################################################################

subroutine find_nearest_index (latb, jindx2)

real, dimension(:), intent(in) :: latb
integer, dimension(:), intent(out)  :: jindx2


      integer :: jd, j, jj
      real   :: diff_low, diff_high
      real, dimension(size(latb,1)-1) :: lat


      jd = size(latb,1) - 1

      do j = 1,jd
        lat(j) = 0.5*(latb(j) + latb(j+1))
      do jj=1, LATOBS         
        if (lat(j)*radian >= cloud_lats(jj)) then
          diff_low = lat(j)*radian - cloud_lats(jj)    
          diff_high = cloud_lats(jj+1) - lat(j)*radian
          if (diff_high <= diff_low) then
            jindx2(j) = jj+1
          else
            jindx2(j) = jj
          endif
        endif
      end do
      end do






end subroutine find_nearest_index 





!######################################################################

subroutine mgrp_prscr_calc (is, ie, js, je, Cld_diagnostics, deltaz, press, temp, &
                            camtsw, cmxolw,  &
                             crndlw, ncldsw, &
				   nmxolw, nrndlw,                    &
				            emmxolw, emrndlw, &
				   cirabsw, cirrfsw, cvisrfsw, cldext, &
				   cldsct, cldasymm)

integer, intent(in) :: is, ie, js, je
real,    dimension(:,:,:),   intent(in) :: deltaz, press, temp     
real,    dimension(:,:,:),   intent(inout) :: camtsw, cmxolw, crndlw
real,    dimension(:,:,:,:), intent(inout) ::   &                  
					      emmxolw, emrndlw
type(cld_diagnostics_type), intent(inout) :: Cld_diagnostics
real,    dimension(:,:,:,:), intent(inout), optional  ::    &
					      cirabsw, cvisrfsw, &
					      cirrfsw, cldext, cldsct, &
					      cldasymm
integer, dimension(:,:),     intent(inout) :: ncldsw, nmxolw, nrndlw

!-------------------------------------------------------------------
!   these arrays define the basic radiative properties of the clouds.

!     cmxolw  =  amounts of maximally overlapped longwave clouds in
!                layers from KSRAD to KERAD.
!     crndlw  =  amounts of randomly overlapped longwave clouds in
!                layers from KSRAD to KERAD.
!     nmxolw  =  number of maximally overlapped longwave clouds
!                in each grid column.
!     nrndlw  =  number of maximally overlapped longwave clouds
!                in each grid column.
!     emmxolw =  longwave cloud emissivity for maximally overlapped
!                clouds through layers from KSRAD to KERAD.
!     emrndlw =  longwave cloud emissivity for randomly overlapped
!                clouds through layers from KSRAD to KERAD.
!     camtsw  =  shortwave cloud amounts. the sum of the maximally
!                overlapped and randomly overlapped longwave
!                cloud amounts.
!     ncldsw  =  number of shortwave clouds in each grid column.
!     cirabsw =  absorptivity of clouds in the infrared frequency band.
!                may be zenith angle dependent.
!     cirrfsw =  reflectivity of clouds in the infrared frequency band.
!                may be zenith angle dependent.
!    cvisrfsw =  reflectivity of clouds in the visible frequency band.
!                may be zenith angle dependent.
!---------------------------------------------------------------------
 
       integer :: i,j,k
       real, dimension(:,:,:,:), allocatable       :: abscoeff, cldemiss
       real, dimension(:,:,:), allocatable    :: conc_drop, conc_ice, &
                                                 size_drop, size_ice

       integer :: israd, ierad, jsrad, jerad


       israd = 1
       ierad = size(camtsw, 1)
       jsrad = 1
       jerad = size(camtsw,2)

!---------------------------------------------------------------------
!     define the number of clouds in each column. the assumption is made
!     that all grid columns have NOFCLDS_SP clouds, with NOFMXOLW of 
!     these being maximally overlapped and the remainder (NOFRNDLW) 
!     randomly overlapped. The default case, corresponding to all 
!     previous model simulations, is 3 clouds (1 maximally overlapped, 
!     2 randomly overlapped).
!----------------------------------------------------------------------

       nmxolw(:,:) = NOFMXOLW
       nrndlw(:,:) = NOFRNDLW
       ncldsw(:,:) = NOFCLDS_SP 

       do k=KSRAD,KERAD
         do j=JSRAD,JERAD
           do i=ISRAD,IERAD
!	     cmxolw(i,j,k)     = zcamtmxo(jabs (j),k)
!	     crndlw(i,j,k)     = zcamtrnd(jabs (j),k)
            cmxolw(i,j,k)     = zcamtmxo(j+js-1  ,k)
            crndlw(i,j,k)     = zcamtrnd(j+js-1  ,k)
             camtsw(i,j,k)     = cmxolw(i,j,k) + crndlw(i,j,k)
           end do
         end do
       end do

       if (.not. do_lwcldemiss) then
         do k=KSRAD,KERAD
           do j=JSRAD,JERAD
             do i=ISRAD,IERAD
! 	       emmxolw(i,j,k,:)  = zemmxo  (jabs (j),k,:)
!              emrndlw(i,j,k,:)  = zemrnd  (jabs (j),k,:)
!              emmxolw(i,j,k,:)  = zemmxo  (j+js-1  ,k,:)
!              emrndlw(i,j,k,:)  = zemrnd  (j+js-1  ,k,:)
               emmxolw(i,j,k,:)  = zemmxo  (j+js-1  ,k,1)
               emrndlw(i,j,k,:)  = zemrnd  (j+js-1  ,k,1)
             end do
           end do
         end do
       endif

!      if (trim(swform) /= 'esfsw99' ) then
       if (.not. do_esfsw            ) then
         do k=KSRAD,KERAD
           do j=JSRAD,JERAD
             do i=ISRAD,IERAD
! 	       cirabsw(i,j,k,:)  = zcabir  (jabs (j),k)
! 	       cirrfsw(i,j,k,:)  = zcrfir  (jabs (j),k)
!       cvisrfsw(i,j,k,:) = zcrfvis (jabs (j),k)
                cirabsw(i,j,k,:)  = zcabir  (j+js-1  ,k)
                cirrfsw(i,j,k,:)  = zcrfir  (j+js-1  ,k)
                cvisrfsw(i,j,k,:) = zcrfvis (j+js-1  ,k)
             end do
           end do
         end do
       endif

!--------------------------------------------------------------------
!  if microphysics is being used for either sw or lw calculation, call
!  microphys_rad_driver to obtain cloud radiative properties.
!--------------------------------------------------------------------
!      if (do_esfsw                  ) then
!        call microphys_rad_driver (is, ie, js, je, deltaz, &
!       press, temp,         camtsw, emmxolw, &
!                                   emrndlw, Cld_diagnostics, &
!                                   cldext=cldext, cldsct=cldsct,    &
!                                   cldasymm=cldasymm)
!      else if (do_lwcldemiss) then
!        call microphys_rad_driver (is, ie, js, je, deltaz, &
!             press, temp,    camtsw, emmxolw,  &
!                                    emrndlw, Cld_diagnostics)
!      endif


 if (do_esfsw                  .or. do_lwcldemiss) then



  allocate (abscoeff ( SIZE(emmxolw,1), SIZE(emmxolw,2),SIZE(emmxolw,3), &
                       SIZE(emmxolw,4)) )
  allocate (cldemiss ( SIZE(emmxolw,1), SIZE(emmxolw,2),SIZE(emmxolw,3), &
                     SIZE(emmxolw,4)) )
allocate (conc_drop (SIZE(emmxolw,1), SIZE(emmxolw,2),SIZE(emmxolw,3)) )
  allocate (conc_ice (SIZE(emmxolw,1), SIZE(emmxolw,2),SIZE(emmxolw,3)) )
  allocate (size_drop (SIZE(emmxolw,1), SIZE(emmxolw,2),SIZE(emmxolw,3)) )
  allocate (size_ice (SIZE(emmxolw,1), SIZE(emmxolw,2),SIZE(emmxolw,3)) )

       call microphys_presc_conc (  is,ie,js,je,             &
                             camtsw,  deltaz, press, temp,       &
                           conc_drop, size_drop, conc_ice, size_ice)


      if (do_esfsw                  ) then


         call microphys_rad_driver (is,ie,js,je,deltaz, press, temp, &
                           Cld_diagnostics, &
                           conc_drop_in=conc_drop, conc_ice_in=conc_ice, &
                      size_drop_in=size_drop, size_ice_in=size_ice, &
                               cldext=cldext, cldsct=cldsct,    &
                              cldasymm=cldasymm,               &
                              abscoeff=abscoeff)

    else if (do_lwcldemiss) then
!       call microphys_rad_driver (                           &
        call microphys_rad_driver (is,ie,js,je,deltaz,press,temp, &
                           Cld_diagnostics, &
                   conc_drop_in=conc_drop, conc_ice_in=conc_ice, &
                  size_drop_in=size_drop, size_ice_in=size_ice, &
                      abscoeff=abscoeff)
     endif  ! do_esfsw

     if (do_lwcldemiss) then
   call lwemiss_calc( deltaz,                             &
                         abscoeff,                         &
                         cldemiss)
      endif  
      endif  ! do_esfsw or lwcldemiss


   if (do_lwcldemiss) then
!   as of 3 august 2001, we are assuming randomly overlapped clouds
!  only. the cloud and emissivity settings follow:

     emmxolw = cldemiss
     emrndlw = cldemiss
 endif  





      if (ALLOCATED (abscoeff)) then
    deallocate (abscoeff)
    deallocate (cldemiss)
   deallocate (conc_drop)

   deallocate (size_drop)
  deallocate (conc_ice)
   deallocate (size_ice)
 endif


end subroutine mgrp_prscr_calc


!##################################################################


subroutine cldht 
 
!------------------------------------------------------------------
!  This subroutine computes the heights of the cloud tops
!  and bottoms for the fixed cloud model.  The observed data
!  are from London (1954, 1957).  This data is a function of 10 deg.
!  latitude bands (0-10, 10-20 and etc.), season and height 
!  in the orginal paper and only for
!  the Northern Hemisphere for various cloud types.
!  Dick and Suki averaged the four seasons together to get annual 
!  mean cloud heights for three type of clouds (hi, middle and low).
!  Somebody also interpolated the data from the 10 deg latitude
!  bands to 5 deg bands.  At the equator, this interpolation
!  was more like an extrapolation.
!  These heights were then put in pressure coordinates using
!  a Skew-T diagram which assumes a "standard atmosphere".
!
!  The Dick and Suki (Ron's pressures) data follow:
!
! TYPE:      hi      |   middle    |            low
!     |              |             |     ctop    |   cbase
! Lat |   ht   press |  ht   press |  ht   press |  ht   press
! 0   |  9.80  272   | 4.35  590   | 3.00  700   | 1.40  855
! 5   |  9.82  271   | 4.40  586   | 3.04  698   | 1.47  848
! 10  | 10.13  259   | 4.45  581   | 3.08  693   | 1.61  836
! 15  | 10.35  250   | 4.50  578   | 3.08  693   | 1.70  828 
! 20  | 10.50  244   | 4.50  578   | 3.01  699   | 1.72  825
! 25  | 10.50  244   | 4.41  584   | 2.91  710   | 1.71  826
! 30  | 10.38  248   | 4.26  595   | 2.80  719   | 1.70  828
! 35  | 10.03  263   | 4.10  614   | 2.70  729   | 1.65  830
! 40  |  9.44  285   | 3.92  621   | 2.60  735   | 1.58  839
! 45  |  8.65  322   | 3.79  633   | 2.47  750   | 1.50  846
! 50  |  7.97  357   | 3.67  647   | 2.35  760   | 1.40  859
! 55  |  7.55  379   | 3.56  651   | 2.24  770   | 1.31  867
! 60  |  7.29  392   | 3.51  657   | 2.17  780   | 1.25  871
! 65  |  7.13  401   | 3.50  658   | 2.10  783   | 1.20  875
! 70  |  7.03  406   | 3.48  659   | 2.03  788   | 1.12  881
! 75  |  7.01  409   | 3.44  660   | 1.98  795   | 1.05  890
! 80  |  6.99  410   | 3.43  661   | 1.91  800   | 1.02  891
! 85  |  6.98  411   | 3.43  661   | 1.88  803   | 1.00  896
! 90  |  6.98  411   | 3.43  661   | 1.87  804   | 1.00  896
!
!  Note that the heights are in kilometers and the pressures
!  in millibars.
!--------------------------------------------------------------------

      real, dimension(10)   :: ciht, asht, cltop, clbase
      integer               :: n, nl, j, k
 
!------------------------------------------------------------------
!  The data in the following arrays run pole to equator,
!  starting at 90, 80, 70....10, 0.
!------------------------------------------------------------------
 
!  Observed high cloud heights running pole to equator (mb)
      data ciht / 411, 410, 406, 392, 357, 285, 248, 244, 259, 272 /
 
!  Observed middle cloud heights running pole to equator (mb)
      data asht / 661, 661, 659, 657, 647, 621, 595, 578, 581, 590 /
 
!  Observed low cloud top heights running pole to equator (mb)
!     data cltop / 804, 800, 788, 780, 760, 735, 719, 699, 693, 700 /
      data cltop / 804, 800, 788, 780, 766, 735, 719, 699, 693, 700 /
!  The above data statement was changed so that this code would
!  reproduce exactly the indexes used in the 9 level model.  The 6 mb
!  error should not affect the results very much....if at all.
 
!  Observed low cloud bottom heights running pole to equator (mb)
      data clbase / 896, 891, 881, 871, 859, 839, 828, 825, 836, 855 /
 
 
!--------------------------------------------------------------------
!   for (at present) unexplained reasons, the middle cloud 
!   specification does not agree with the current gcm specification.
!   a glance at the telegadas and london paper suggests differences
!   caused by 1) use of cldmid rather than clotop and cldbase heights
!   2) seeming errors in polar cloud height locations. for now,
!   the foregoing kluge gives cloud positions in agreement with
!   SKYHI, PROVIDED:
!    it is realized that these indexes are the true values, which
!   differ by 1 from the skyhiindices. in the gcm, a subtraction
!   is done to get the indices to be correct (see radmn.F).
!--------------------------------------------------------------------

      do n=8,10
	asht(n) = asht(n) - 25.
      enddo
 
 
!-------------------------------------------------------------------
!  First compute the cloud top indexes for the high clouds
!  nl is the index for cloud type nl=3 => high
!--------------------------------------------------------------------  

      nl = 3
      call Cldint(ciht ,kkth, nl)
 
!---------------------------------------------------------------------
!  Set cloud bottom height equal to cloud top height.
!  This assumes cirrus clouds are one level thick.
!---------------------------------------------------------------------

      do j=1,latobs
       kkbh(j,nl) = kkth(j,nl)
      end do
 
!-------------------------------------------------------------------
!  Second compute the cloud top indexes for the middle clouds
!  nl is the index for cloud type nl=2 => middle
!-------------------------------------------------------------------
  
      nl = 2
      call Cldint(asht,kkth,  nl)
 
!-------------------------------------------------------------------
!  Set cloud bottom height equal to cloud top height.
!  This assumes middle clouds are one level thick.
!-------------------------------------------------------------------

      do j=1,latobs
       kkbh(j,nl) = kkth(j,nl)
      end do
 
!-------------------------------------------------------------------
!  Third compute the cloud top indexes for the low clouds
!  nl is the index for cloud type nl=1 => low
!-------------------------------------------------------------------
  
      nl = 1
      call Cldint(cltop,kkth,  nl)
 
!-------------------------------------------------------------------
!  Lastly compute the cloud bottom indexes for the low clouds
!  This assumes that low clouds can be thicker than one level.
!  nl is the index for cloud type nl=1 => low
!-------------------------------------------------------------------
  
      nl = 1
      call Cldint(clbase, kkbh, nl)
 
end subroutine cldht


!##################################################################

subroutine cldint(cldobs, kindex, nl)
 
!-------------------------------------------------------------------
!  This subroutine computes the indexes for the heights of the cloud 
!  tops and bottoms for the fixed cloud model. 
!-------------------------------------------------------------------
 
real, dimension(10), intent(in)                    :: cldobs
integer, dimension(LATOBS,NOFCLDS_SP), intent(out) :: kindex
integer,                               intent(in)  :: nl
!------------------------------------------------------------------
 
      real, dimension(LATOBS)     :: cldlat
      real                        :: alevel, prsmid
      integer                     :: j, k
 
!---------------------------------------------------------------------
!  Fill in Southern hemisphere cloud heights
!---------------------------------------------------------------------

      do j=1,10
       cldlat(j) = cldobs(j)
      end do  
      do j=1,9
       cldlat(j+10) = cldobs(10-j)
      end do  
  
!-------------------------------------------------------------------
!  Start latitude loop to compute index at each latitude.
!-------------------------------------------------------------------
      do j=1,LATOBS
 
!-------------------------------------------------------------------
!  Find first place where the pressure on the model level
!  is greater than the pressure of the cloud height.  Starting
!  from the top of the atm and going down the column.
!-------------------------------------------------------------------
        if (qlevel(KERAD)*psbar .lt. cldlat(j)) then       
        call error_mesg ( 'cldint', &
         'no level found with pressure greater than cloud pressure', & 
								FATAL)
	endif
        if (qlevel(1)*psbar .gt. cldlat(j)) then       
	  call error_mesg ( 'cldint',  &
	                ' cloud is above highest model level', FATAL)
        endif
	do k=KSRAD,KERAD
	  alevel = qlevel(k)*psbar
          if (alevel .gt. cldlat(j)) then       
!-------------------------------------------------------------------
!  k is the index of the first model level below the cloud height.
!  compute the pressure half way between the model levels
!-------------------------------------------------------------------
            prsmid = (qlevel(k)+qlevel(k-1))*psbar*0.5   
 
!-------------------------------------------------------------------
!  If prsmid is greater than cldlat (cloud height) then the
!  level above is closer to cloud height, otherwise it is the
!  level below.
!-------------------------------------------------------------------
            if (prsmid .gt. cldlat(j)) then
              kindex(j,nl) = k-1
            else
              kindex(j,nl) = k
            endif
            exit
	  endif
        end do
      end do
!--------------------------------------------------------------------

 

end subroutine cldint




!####################################################################

	       end module mgrp_prscr_clds_mod




