
		 module radiation_diag_mod


use rad_utilities_mod,       only: longwave_control_type, Lw_control, &
				   shortwave_control_type, Sw_control,&
				   Environment, environment_type, &
				   radiation_control_type, Rad_control
use  utilities_mod,          only: open_file, file_exist,    &
                                   check_nml_error, error_mesg, &
				   print_version_number, FATAL, NOTE, &
				   WARNING, get_my_pe, close_file, &
				   get_domain_decomp, get_num_pes
use constants_new_mod,       only: radcon, radians_to_degrees
use esfsw_parameters_mod,    only: nbands
use rad_step_setup_mod,      only: jabs, iabs, ISRAD, IERAD, JSRAD,   &
				   JERAD, KS=>KSRAD, KE=>KERAD

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                 radiation diagnostics module
!          this module is the interface between the radiation code 
!                 and the model diagnostics handler
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.09'
    character(len=128)  :: version =  '$Id: radiation_diag.F90,v 1.2 2001/07/05 17:32:22 fms Exp $'
    character(len=128)  :: tag     =  '$Name: eugene $'

!---------------------------------------------------------------------
!------    interfaces   ------

public           &
	  radiation_diag_init, radiation_diag_end, &
	  radiag_from_lwtables, radiag_driver, &
	  radiag_from_setup, radiag_setup_dealloc,   &
	  radiag_from_driver, radiag_from_cts, &
	  radiag_from_radgases,  radiation_diag_dealloc, &
	  radiag_from_fluxes,   &
	  radiag_from_ctsappx, radiag_from_sw_driver_init, &
	  radiag_from_sw_driver, radiag_from_clouds_lh, &
	  radiag_from_clouds_esf, &
          radiag_from_astronomy, radiag_from_sfcalbedo, &
	  radiag_from_ozone, radiag_from_cloudrad


private          &
	  radiag

!---------------------------------------------------------------------
!------     namelist  -----

integer, parameter                 :: max_pts = 20
logical                            :: do_totcld_forcing = .false.
integer, dimension (max_pts)       :: iradprt=0, jradprt=0
integer                            :: num_pts = 0
character(len=12)                  :: radiag_output_file=' '
logical                            :: write_radiag_file=.false.


namelist / radiation_diag_nml /  &
                                  do_totcld_forcing, &
				  iradprt, jradprt, &
				  num_pts, write_radiag_file, &
				  radiag_output_file


!----------------------------------------------------------------------
!---  public variables ---


	       
!----------------------------------------------------------------------
!---  private variables ---


!---------------------------------------------------------------------
!   the following variables are passed into this module via subroutine
!   calls from the modules where they reside, mimicking the anticip-
!   ated diagnostics handler.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! from longwave_driver_mod:
!     flxnet  =  net longwave flux at model flux levels (including the 
!                ground and the top of the atmosphere).
!     heatra  =  heating rate at data levels.
!     flx1e1  =  flux at top of atmosphere for 0-160, 1200-2200
!                cm-1 range.
!     flx1e1f =  flux at top of atmosphere for NBTRGE bands in 1200-
!                1400 cm-1 range.
!     cts     =  approximate cool-to-space heating rates for 160-560
!                and 800-990, 1070-1200 cm-1 ranges.
!     ctsco2  =  approximate cool-to-space heating rates for 560-800
!                cm-1 range.
!     ctso3   =  approximate cool-to-space heating rates for 990-1070
!                cm-1 range.
!     excts    =  exact cool-to-space heating rates for 160-1200 cm-1
!                 range.
!     gxcts    =  flux at top of atmosphere for 160-1200 cm-1 range. 
!---------------------------------------------------------------------
real, dimension(:,:),      allocatable  :: cts, ctsco2, ctso3
real, dimension(:,:),      allocatable  :: excts 
real, dimension(:),        allocatable  :: gxcts                
real, dimension(:),        allocatable  :: flx1e1
real, dimension(:,:),      allocatable  :: flx1e1f
real, dimension(:,:),      allocatable  :: heatra, heatracf, flxnet,&
			                   flxnetcf


!--------------------------------------------------------------------
! from longwave_fluxes_mod:
!if(.not. do_ch4n2o), then
!        flx1    =  net emissivity flux for 0-160, 1200-2200 cm-1 band.
!        flx7    =  flux for 1200-1400 cm-1 band (NBTRGE bands).
!else if (do_ch4_n2o) then
!        flx1    =  net emissivity flux for 0-160, 1400-2200 cm-1 band.
!endif
!        flx2    =  flux for 560-800 cm-1 band (as one band).
!        flx3    =  flux for 800-900 cm-1 band.
!        flx4    =  flux for 900-990 cm-1 band.
!        flx5    =  flux for 990-1070 cm-1 band.
!        flx6    =  flux for 1070-1200 cm-1 band.
!---------------------------------------------------------------------
!if(do_totcld_forcing) then
!     flx1cf, flx2cf, flx3cf, flx4cf, flx5cf, flx6cf = cloud-free
!     emissivity flux for bands corresponding to flx1, flx2, flx3, flx4,
!     flx5, flx6.
!     (if do_ch4n2o option is true)
!     flx7cf = cloud-free emissivity flux for flx7 band.
!     endif
!--------------------------------------------------------------------
real, dimension(:,:),      allocatable  :: flx1, flx2, flx3, flx4, &
   			                   flx5, flx6
real, dimension(:,:,:),    allocatable  :: flx7         
real, dimension(:,:),      allocatable  :: flx1cf, flx2cf, flx3cf, &
				           flx4cf, flx5cf, flx6cf
real, dimension(:,:,:),    allocatable  :: flx7cf       


!--------------------------------------------------------------------
! from longwave_setup_mod:
!     pflux   =  pressure at flux levels of model.
!     press   =  pressure at data levels of model.
!     rh2o    =  mass mixing ratio of h2o at model data levels.
!     temp    =  temperature at data levels of model.
!     pdfinv  =  inverse of pressure difference between flux levels.
!     pdflux  =  pressure difference between flux levels.
!--------------------------------------------------------------------
real, dimension(:,:),      allocatable  :: pdflux, pdfinv, pflux, &
				           temp, press, rh2o


!--------------------------------------------------------------------
! from cool_to_space_mod:
!     exctsn   =  exact cool-to-space heating rates for each band.
!     fctsg    =  cool-to-space flux at the ground for each band.
!------------------------------------------------------------------
real, dimension(:,:),      allocatable  :: fctsg       
real, dimension(:,:,:),    allocatable  :: exctsn       


!--------------------------------------------------------------------
! from radiative_gases_mod:
!--------------------------------------------------------------------
real                                    :: rrvf11, rrvf12, rrvf113, &
				           rrvf22, rrvch4, rrvn2o,  &
					   rrvco2


!--------------------------------------------------------------------
! from longwave_tables_mod:
!--------------------------------------------------------------------
real, dimension(:),        allocatable  :: bandlo, bandhi, bdlocm, &
				           bdhicm
integer, dimension(:),     allocatable  :: iband


!---------------------------------------------------------------------
! from cloudrad_mod:
!     deltaz  =  altitude difference between flux levels
!--------------------------------------------------------------------
real, dimension (:,:,:),   allocatable  :: cldext, cldssalb, cldasymm
real, dimension (:,:),     allocatable  :: deltaz
real, dimension (:,:),     allocatable  :: lwpath, iwpath,       &
                                           size_drop, size_ice


!--------------------------------------------------------------------
! from shortwave_driver_mod:
!     dfsw    =  downward short-wave radiation at pressure levels
!     fsw     =  net radiation (up-down) at all pressure levels.
!     gwt     = gaussian weights
!     hsw     =  radiation heating rates at all pressure layers.
!     nsolwg  =  number of gaussian points(or 1 if not gaussian).
!     ssolar  =  solar constant (may vary over one year). units:W/m2.
!     ufsw    =  upward radiation at all pressure levels.
!--------------------------------------------------------------------
real                                    ::  ssolar
logical                                 ::  lswg, ldiurn, do_lhsw, &
					    do_esfsw, do_annual
integer                                 ::  nsolwg
real, dimension(:), allocatable         ::  gwt
real, dimension(:,:), allocatable       ::  dfsw, fsw, hsw, ufsw, &
					    dfswcf, fswcf, hswcf, &
					    ufswcf


!--------------------------------------------------------------------
! from cloud_package:
!     cmxolw  =  amounts of maximally overlapped longwave clouds in
!                layers from KS to KE.
!     crndlw  =  amounts of randomly overlapped longwave clouds in
!                layers from KS to KE.
!     camtsw  =  shortwave cloud amounts. their locations are specified
!                in the ktopsw/kbtmsw indices. when no clouds are 
!                present camtsw should be given a value of zero.
!     cirabsw =  absorptivity of clouds in the infrared frequency band.
!                when no clouds are present cirabsw should be given a
!                value of zero.
!     cirrfsw =  reflectivity of clouds in the infrared frequency band.
!                when no clouds are present cirrfsw should be given a
!                value of zero.
!     cvisrfsw=  reflectivity of clouds in the visible frequency band.
!                when no clouds are present cvisrfsw should be given a
!                value of zero.
!     emmxolw =  longwave cloud emissivity for maximally overlapped
!                clouds through layers from KS to KE. (default is one).
!     emrndlw =  longwave cloud emissivity for randomly overlapped
!                clouds through layers from KS to KE. (default is one).
!     kbtmsw  =  index of flux level pressure of cloud bottom.  when no
!                clouds are present kbtmsw should be given a value of
!                KS.
!     ktopsw  =  index of flux level pressure of cloud top.  when no
!                clouds are present ktopsw should be given a value of
!                KS.
!     nmxolw  =  number of maximally overlapped longwave clouds 
!                at each grid point.
!     nrndlw  =  number of maximally overlapped longwave clouds 
!                at each grid point.
!     ncldsw  =  number of clouds at each grid point.
!--------------------------------------------------------------------
real, dimension(:,:,:), allocatable     :: cirabsw, cirrfsw, cvisrfsw,&
					   emmxolw, emrndlw
real, dimension(:,:),   allocatable     :: camtsw, cmxolw, crndlw
integer, dimension(:),  allocatable     :: ncldsw, nrndlw, nmxolw
integer, dimension(:,:),allocatable     :: kbtmsw, ktopsw


!--------------------------------------------------------------------
! from astronomy_package:
!     fracday = fraction of averaging period that has daylight.
!     zenith  = mean cosine of zenith angle for all longitudes.
!--------------------------------------------------------------------
real, dimension(:,:), allocatable     :: cosang                      
real, dimension(:),   allocatable     :: fracday


!--------------------------------------------------------------------
! from sfcalbedo_package:
!--------------------------------------------------------------------
real, dimension(:),    allocatable    :: cvisrfgd, cirrfgd


!--------------------------------------------------------------------
! from ozone_package:
!     qo3     =  mass mixing ratio of o3 at model data levels.
!--------------------------------------------------------------------
real, dimension(:,:),  allocatable    :: qo3


!-------------------------------------------------------------------

real,    dimension(:),      allocatable ::  deglon1, deglat1
integer, dimension (max_pts)            :: jradprt_gl, iradprt_gl


logical             :: cld_flg=.false.
integer             :: NBLY, NBTRGE, NBLW
integer             :: radiag_unit
integer             :: x(4), y(4)
integer             :: stdout1=6  


                       contains





subroutine radiation_diag_init (rlat, rlong)

!---------------------------------------------------------------------
real, dimension(:),   intent(in)    :: rlat
real, dimension(:,:), intent(in)    :: rlong
!---------------------------------------------------------------------

   integer       :: unit, ierr, io, nn, j, jdf

!---------------------------------------------------------------------
!----  read namelist

   if (file_exist('input.nml')) then
     unit =  open_file ('input.nml', action='read')
     ierr=1; do while (ierr /= 0)
     read (unit, nml=radiation_diag_nml, iostat=io, end=10)
     ierr = check_nml_error (io, 'radiation_diag_nml')
     enddo
10   call close_file (unit)
   endif
   unit = open_file ('logfile.out', action='append')
!  call print_version_number (unit, 'radiation_diag', version_number)
   if (get_my_pe() == 0) then
     write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
     write (unit,nml=radiation_diag_nml)
   endif
   call close_file (unit)

   Rad_control%do_totcld_forcing = do_totcld_forcing

   call get_domain_decomp (x, y)
   jdf = y(4) - y(3) + 1

!-------------------------------------------------------------------
!     define latitude rows where diagnostics are desired
!-------------------------------------------------------------------
   if (num_pts > 0) then
     allocate ( deglon1 (num_pts))
     allocate ( deglat1 (num_pts))
   endif

   allocate (Rad_control%do_raddg (jdf) )

   Rad_control%do_raddg = .false.
   do nn=1, num_pts
     jradprt_gl(nn) = jradprt(nn)
     iradprt_gl(nn) = iradprt(nn)
     iradprt(nn)    = 0
     jradprt(nn)    = 0
     do j=y(3), y(4)
       if (jradprt_gl(nn) == j .and.   &
	     iradprt_gl(nn) >= x(3) .and. iradprt_gl(nn) <= x(4) ) then
         Rad_control%do_raddg(j-y(3)+1) = .true.
	 jradprt(nn) = j-y(3)+1
	 iradprt(nn) = iradprt_gl(nn)-x(3)+1
     deglon1(nn) = rlong(iradprt(nn), jradprt(nn))*radians_to_degrees
     deglat1(nn) = rlat(jradprt(nn))*radians_to_degrees
         exit
       endif
     end do
   end do

   if (num_pts > 0) then
     if (write_radiag_file) then
       radiag_unit = open_file (radiag_output_file, action='write', &
				threading='multi', form='formatted')
     else
       radiag_unit = stdout1
     endif
   endif    
!--------------------------------------------------------------------


end subroutine radiation_diag_init





!####################################################################

subroutine radiation_diag_end

  if (num_pts > 0 ) then
    if (write_radiag_file) then
      call close_file (radiag_unit)
      if (get_my_pe() == 0) write (*, 90006)
    endif
  endif


90006 format (//, ' closed radiag_output_file')

end subroutine radiation_diag_end




!#####################################################################

subroutine radiag_from_lwtables (bdlocm_in, bdhicm_in, iband_in,  &
				 bandlo_in, bandhi_in)

!--------------------------------------------------------------------
real, dimension(:), intent(in)    :: bdlocm_in, bdhicm_in, bandlo_in, &
				     bandhi_in
integer, dimension(:), intent(in) :: iband_in
!--------------------------------------------------------------------

   NBLY = size(bdlocm_in)
   NBLW = size(bandlo_in)

   allocate ( bdlocm (NBLY) )
   allocate ( bdhicm (NBLY) )
   allocate ( iband  (40  ) )
   allocate ( bandlo (NBLW) )
   allocate ( bandhi (NBLW) )

!--------------------------------------------------------------------
   bdlocm(:) = bdlocm_in(:)
   bdhicm(:) = bdhicm_in(:)
   iband (:) = iband_in(:)
   bandlo(:) = bandlo_in(:)
   bandhi(:) = bandhi_in(:)
   
end subroutine radiag_from_lwtables



!#####################################################################

subroutine radiag_driver

      integer              :: ip, np
      character(len=13)    :: acc
!--------------------------------------------------------------------

      if (Environment%running_skyhi) then
        call barrier()
        if (get_my_pe() ==  0) then
          inquire (stdout1, access = acc)
          if (acc /= 'UNDEFINED    ') then
            call flush (stdout1)
          endif
        end if
        np = get_num_pes()
        do ip = 0,np-1
          if (get_my_pe() == ip) then
	    call radiag
            inquire (stdout1, access = acc)
            if (acc /= 'UNDEFINED    ') then
              call flush (stdout1)
            endif
          endif
          call barrier()
        end do
      else if (Environment%running_fms) then
        inquire (stdout1, access = acc)
         if (acc /= 'UNDEFINED    ') then
           call flush (stdout1)
         endif

!--------------------------------------------------------------------
!   call radiag to compute radiation diagnostics at desired points
!--------------------------------------------------------------------
         call radiag

         inquire (stdout1, access = acc)
         if (acc /= 'UNDEFINED    ') then
           call flush (stdout1)
         endif
      endif

end subroutine radiag_driver



!##################################################################

subroutine radiag_from_setup (jloc, pdflux_in, pdfinv_in, pflux_in, &
    		              temp_in, press_in, rh2o_in)

!--------------------------------------------------------------------
real, dimension(:,:,:), intent(in)     :: pdflux_in, pdfinv_in, pflux_in
real, dimension(:,:,:), intent(in)     :: press_in, temp_in, rh2o_in
integer,                intent(in)     :: jloc  
!--------------------------------------------------------------------

  integer    ::  k, nn, iloc

!--------------------------------------------------------------------
!allocate arrays to receive desired data from longwave_setup_mod
!--------------------------------------------------------------------
  if (.not. allocated (temp)) then
  allocate ( temp   (KS:KE+1, num_pts) )
  allocate ( press  (KS:KE+1, num_pts) )
  allocate ( rh2o   (KS:KE,   num_pts) )
  allocate ( pdfinv (KS:KE,   num_pts) )
  allocate ( pdflux (KS:KE,   num_pts) )
  allocate ( pflux  (KS:KE+1, num_pts) )
  endif

!---------------------------------------------------------------------
!load desired data
!---------------------------------------------------------------------
  do nn=1, num_pts
    if (jabs(jloc) == jradprt(nn)) then
      if ( (iradprt(nn) >= iabs(ISRAD)) .and.   &
	   (iradprt(nn) <= iabs(IERAD)) ) then
	iloc = iradprt(nn) - iabs(ISRAD) + 1
        do k=KS,KE
          pdflux(k,nn) = pdflux_in(iloc,jloc,k)
          pdfinv(k,nn) = pdfinv_in(iloc,jloc,k)
          rh2o(k,nn)   = rh2o_in(iloc,jloc,k)
        end do
        do k=KS,KE+1
          pflux(k,nn)  = 0.10*pflux_in(iloc,jloc,k)
          press(k,nn)  = 0.10*press_in(iloc,jloc,k)
          temp (k,nn)  = temp_in(iloc,jloc,k)
        end do
      endif
    endif
  end do

!--------------------------------------------------------------------
   


end subroutine radiag_from_setup 

 



!##################################################################

subroutine radiag_setup_dealloc

!---------------------------------------------------------------------
!deallocate longwave_setup_mod arrays
!---------------------------------------------------------------------

  deallocate (pflux )
  deallocate (pdflux)
  deallocate (pdfinv)
  deallocate (rh2o    )
  deallocate (press   )
  deallocate (temp    )


end subroutine radiag_setup_dealloc


!#####################################################################

subroutine radiag_from_driver (jloc, flx1e1_in, flx1e1f_in,   &
			       heatra_in, heatracf_in, flxnet_in,  &
			       flxnetcf_in)

!--------------------------------------------------------------------
real, dimension(:,:),     intent(in) :: flx1e1_in
real, dimension(:,:,:),   intent(in) :: flx1e1f_in, heatra_in, &
			                heatracf_in, flxnet_in,   &
					flxnetcf_in
integer,                  intent(in) :: jloc                
!--------------------------------------------------------------------

     integer    ::  k,         nn, m, iloc

!---------------------------------------------------------------------
!allocate arrays to store needed data from longwave_driver_mod
!---------------------------------------------------------------------
     if (.not. allocated (heatra) ) then
       allocate ( flx1e1(        num_pts) )
       if (Lw_control%do_ch4_n2o) then
         allocate ( flx1e1f (  NBTRGE, num_pts) )
       endif

       allocate (heatra  (KS:KE, num_pts) )
       allocate (flxnet  (KS:KE+1, num_pts) )
       if (do_totcld_forcing) then
         allocate (heatracf  (KS:KE, num_pts) )
         allocate (flxnetcf  (KS:KE+1, num_pts) )
       endif
     endif

!---------------------------------------------------------------------
!move desired data into the prepared arrays
!---------------------------------------------------------------------
     do nn=1, num_pts
       if (jabs(jloc) == jradprt(nn)) then
	 if (iradprt(nn) >= iabs(ISRAD) .and.    &
	     iradprt(nn) <= iabs(IERAD)) then
	   iloc = iradprt(nn) - iabs(ISRAD) + 1
           do k=KS,KE+1
             flxnet(k,nn) = flxnet_in(iloc,jloc, k)
             if (do_totcld_forcing) then
               flxnetcf(k,nn) = flxnetcf_in(iloc, jloc, k)
             endif
           end do
           do k=KS,KE
             heatra(k,nn) = heatra_in(iloc, jloc, k)
             if (do_totcld_forcing) then
               heatracf(k,nn) = heatracf_in(iloc, jloc, k)
             endif
           end do
           flx1e1(nn) = flx1e1_in(iloc, jloc)
           if (Lw_control%do_ch4_n2o) then
             do m=1,NBTRGE
               flx1e1f(m,nn) = flx1e1f_in(iloc, jloc,m)
             end do
           endif
         endif
       endif
     end do
!---------------------------------------------------------------------


end subroutine radiag_from_driver




!#####################################################################

subroutine radiag_from_cts (jloc, exctsn_in, fctsg_in, excts_in, &
			    gxcts_in)

!--------------------------------------------------------------------
real, dimension(:,:,:,:), intent(in) :: exctsn_in
real, dimension(:,:,:  ), intent(in) :: fctsg_in, excts_in
real, dimension(:,:    ), intent(in) :: gxcts_in
integer,                  intent(in) :: jloc
!--------------------------------------------------------------------

   integer    ::  k,         nn, n, iloc

!--------------------------------------------------------------------
!allocate space for data coming from cool_to_space_exact 
!--------------------------------------------------------------------
   if (.not. allocated (fctsg) ) then
     allocate ( exctsn (KS:KE, NBLY, num_pts) )
     allocate ( fctsg  (       NBLY, num_pts) )
     allocate ( excts ( KS:KE, num_pts) )
     allocate ( gxcts (        num_pts) )
   endif


!--------------------------------------------------------------------
!move data into prepared arrays
!--------------------------------------------------------------------
   do nn = 1, num_pts
     if (jabs(jloc) == jradprt(nn)) then
       if (iradprt(nn) >= iabs(ISRAD) .and.   &
	   iradprt(nn) <= iabs(IERAD)) then
         iloc = iradprt(nn) - iabs(ISRAD) + 1
         do n=1,NBLY
           do k=KS,KE
             exctsn(k,n,nn) = exctsn_in(iloc,jloc,k,n)
           end do
           fctsg(n,nn)  = fctsg_in(iloc,jloc,n)
         end do
         do k=KS,KE
           excts (k,nn) = excts_in(iloc, jloc, k)
         end do
         gxcts (nn) = gxcts_in(iloc, jloc)
       endif
     endif
   end do
!--------------------------------------------------------------------
   


end subroutine radiag_from_cts 



!###################################################################

subroutine radiag_from_radgases (rrvf11_in, rrvf12_in, rrvf113_in,   &
	                 	 rrvf22_in, rrvch4_in, rrvn2o_in ,  &
				 rrvco2_in)

!--------------------------------------------------------------------
real, intent(in) :: rrvf11_in, rrvf12_in, rrvf113_in, rrvf22_in,   &
		    rrvch4_in, rrvn2o_in, rrvco2_in
!--------------------------------------------------------------------

     rrvf11 = rrvf11_in
     rrvf12 = rrvf12_in
     rrvf113 = rrvf113_in
     rrvf22 = rrvf22_in
     rrvch4 = rrvch4_in
     rrvn2o = rrvn2o_in
     rrvco2 = rrvco2_in

!--------------------------------------------------------------------

   

end subroutine radiag_from_radgases



!###################################################################

subroutine radiation_diag_dealloc 

!-------------------------------------------------------------------

     if (do_totcld_forcing) then
       deallocate ( heatracf )
       deallocate ( flxnetcf )
     endif

     deallocate ( flxnet )
     deallocate ( heatra )
   
     if (Lw_control%do_ch4_n2o) then
       deallocate ( flx1e1f  )
     endif

     deallocate ( flx1e1 )
     deallocate ( gxcts  )
     deallocate ( excts  )

     if (do_totcld_forcing) then
       if (Lw_control%do_ch4_n2o) then
         deallocate ( flx7cf )
       endif
       deallocate ( flx6cf )
       deallocate ( flx5cf )
       deallocate ( flx4cf )
       deallocate ( flx3cf )
       deallocate ( flx2cf )
       deallocate ( flx1cf )
     endif

     if (Lw_control%do_ch4_n2o) then
       deallocate ( flx7 )
     endif

     deallocate ( flx6 )
     deallocate ( flx5 )
     deallocate ( flx4 )
     deallocate ( flx3 )
     deallocate ( flx2 )
     deallocate ( flx1 )
     
     deallocate ( ctso3  )
     deallocate ( ctsco2 )
     deallocate ( cts    )

     deallocate ( fctsg  )
     deallocate ( exctsn )

     if (allocated (dfsw) ) then
       deallocate ( dfsw )
       deallocate (  fsw )
       deallocate (  hsw )
       deallocate ( ufsw )

       if (allocated (dfswcf)) then
         deallocate ( dfswcf )
         deallocate (  fswcf )
         deallocate (  hswcf )
         deallocate ( ufswcf )
       endif
     endif

     if (allocated (cldext) ) then
       deallocate (cldext    )
       deallocate (cldssalb  )
       deallocate (cldasymm  )
       deallocate (deltaz    )
     endif

     if (allocated (lwpath) ) then
       deallocate (lwpath)
       deallocate (iwpath)
       deallocate (size_drop)
       deallocate (size_ice)
     endif

     if (allocated (nrndlw) ) then
       if (do_lhsw) then
         deallocate (ktopsw  )
         deallocate (kbtmsw  )
         deallocate ( cvisrfsw   )
         deallocate ( cirrfsw    )
         deallocate ( cirabsw    )
       endif
       deallocate ( nrndlw)
       deallocate ( nmxolw)
       deallocate ( ncldsw)
       deallocate ( emrndlw    )
       deallocate ( emmxolw    )
       deallocate ( crndlw     )
       deallocate ( cmxolw     )
       deallocate ( camtsw     )
     endif

     if (allocated (fracday)) then
       deallocate (fracday  )
       deallocate (cosang )
     endif

     deallocate (cvisrfgd)
     deallocate (cirrfgd )

     deallocate (qo3  )
!---------------------------------------------------------------------


end subroutine radiation_diag_dealloc




!#################################################################

subroutine radiag_from_fluxes (jloc, flx_in, flxcf_in, NBTRGE_pas)

!--------------------------------------------------------------------
integer,                  intent(in) :: NBTRGE_pas
real, dimension(:,:,:,:), intent(in) :: flx_in, flxcf_in
integer,                  intent(in) :: jloc                

!--------------------------------------------------------------------
     integer    ::  k, nn, m, iloc

!--------------------------------------------------------------------
     NBTRGE = NBTRGE_pas

!---------------------------------------------------------------------
!allocate arrays to store needed data from longwave_driver_mod
!---------------------------------------------------------------------
     if (.not. allocated (flx2) ) then
     allocate ( flx1 ( KS:KE+1,num_pts) )
     allocate ( flx2 ( KS:KE+1,num_pts) )
     allocate ( flx3 ( KS:KE+1,num_pts) )
     allocate ( flx4 ( KS:KE+1,num_pts) )
     allocate ( flx5 ( KS:KE+1,num_pts) )
     allocate ( flx6 ( KS:KE+1,num_pts) )
     if (Lw_control%do_ch4_n2o) then
       allocate ( flx7 ( KS:KE+1,NBTRGE,num_pts))
     endif
     if (do_totcld_forcing) then
       allocate ( flx1cf ( KS:KE+1,num_pts) )
       allocate ( flx2cf ( KS:KE+1,num_pts) )
       allocate ( flx3cf ( KS:KE+1,num_pts) )
       allocate ( flx4cf ( KS:KE+1,num_pts) )
       allocate ( flx5cf ( KS:KE+1,num_pts) )
       allocate ( flx6cf ( KS:KE+1,num_pts) )
       if (Lw_control%do_ch4_n2o) then
         allocate ( flx7cf ( KS:KE+1,NBTRGE,num_pts) )
       endif
     endif
     endif

!---------------------------------------------------------------------
!move desired data into the prepared arrays
!---------------------------------------------------------------------
     do nn=1, num_pts
       if (jabs(jloc) == jradprt(nn)) then
	 if (iradprt(nn) >= iabs(ISRAD) .and.   &
	     iradprt(nn) <= iabs(IERAD)) then
	   iloc = iradprt(nn) - iabs(ISRAD) + 1
           do k=KS,KE+1
             flx1(k,nn) = flx_in(iloc,jloc,k,1)
             flx2(k,nn) = flx_in(iloc,jloc,k,2)
             flx3(k,nn) = flx_in(iloc,jloc,k,3)
             flx4(k,nn) = flx_in(iloc,jloc,k,4)
             flx5(k,nn) = flx_in(iloc,jloc,k,5)
             flx6(k,nn) = flx_in(iloc,jloc,k,6)
             if (do_totcld_forcing) then
               flx1cf(k,nn) = flxcf_in(iloc, jloc, k,1)
               flx2cf(k,nn) = flxcf_in(iloc, jloc, k,2)
               flx3cf(k,nn) = flxcf_in(iloc, jloc, k,3)
               flx4cf(k,nn) = flxcf_in(iloc, jloc, k,4)
               flx5cf(k,nn) = flxcf_in(iloc, jloc, k,5)
               flx6cf(k,nn) = flxcf_in(iloc, jloc, k,6)
             endif
             if (Lw_control%do_ch4_n2o) then
               do m=1,NBTRGE
                 if (do_totcld_forcing) then
                   flx7cf(k,m,nn) = flxcf_in(iloc, jloc, k,6+m)
                 endif
                 flx7(k,m,nn) = flx_in(iloc, jloc, k,6+m)
               end do
             endif
           end do
         endif
       endif
     end do
!---------------------------------------------------------------------


end subroutine radiag_from_fluxes



!####################################################################

subroutine radiag_from_ctsappx (jloc, cts_out_in,  index)

!--------------------------------------------------------------------
integer,                       intent(in) :: index
real, dimension(:,:,:),        intent(in) :: cts_out_in
integer,                       intent(in) :: jloc                
!--------------------------------------------------------------------

   integer    ::  k, nn, n, iloc

!--------------------------------------------------------------------
!allocate space for data coming from cool_to_space_approx
!--------------------------------------------------------------------
   if (.not. allocated (cts) ) then
   if (index == 1) then
     allocate ( cts    (KS:KE, num_pts) )
     allocate ( ctsco2 (KS:KE, num_pts) )
     allocate ( ctso3  (KS:KE, num_pts) )
   endif
   endif

!--------------------------------------------------------------------
!move data into prepared arrays
!--------------------------------------------------------------------
   do nn = 1, num_pts
     if (jabs(jloc) == jradprt(nn)) then
       if (iradprt(nn) >= iabs(ISRAD) .and.    &
	   iradprt(nn) <= iabs(IERAD)) then
         iloc = iradprt(nn) - iabs(ISRAD) + 1
         do k=KS,KE
	   if (index == 1) then
             ctsco2(k,nn) = cts_out_in(iloc,jloc,k)
           else if (index == 2) then
	     ctso3(k,nn) = cts_out_in (iloc,jloc,k)
	   else if (index == 3) then
	     cts (k,nn) = cts_out_in (iloc,jloc,k)
           else
	     cts (k,nn) = cts(k,nn) + cts_out_in (iloc,jloc,k)
	   endif
         end do
       endif
     endif
   end do

!--------------------------------------------------------------------

end subroutine radiag_from_ctsappx 



!###################################################################

subroutine radiag_from_sw_driver_init (ldiurn_in, lswg_in, nsolwg_in,  &
				       gwt_in, do_lhsw_in, do_esfsw_in,&
				       do_annual_in)

!--------------------------------------------------------------------
integer,                     intent(in)   :: nsolwg_in
logical,                     intent(in)   :: lswg_in, ldiurn_in, &
					     do_lhsw_in, do_esfsw_in, &
					     do_annual_in
real, dimension(:),          intent(in)   :: gwt_in
!---------------------------------------------------------------------

    allocate ( gwt(nsolwg_in) )

    ldiurn  = ldiurn_in
    lswg    = lswg_in
    nsolwg    = nsolwg_in
    do_lhsw = do_lhsw_in
    do_esfsw = do_esfsw_in
    do_annual = do_annual_in
    gwt(:)  = gwt_in(:)

!---------------------------------------------------------------------

end subroutine radiag_from_sw_driver_init 




!###################################################################

subroutine radiag_from_sw_driver (jloc, &
				   ssolar_in, dfsw_in,  &
				  fsw_in, hsw_in, ufsw_in, dfswcf_in, &
				  fswcf_in, hswcf_in, ufswcf_in)

!--------------------------------------------------------------------
real,                     intent(in)   :: ssolar_in
real, dimension(:,:,:),   intent(in)   :: dfsw_in, fsw_in, hsw_in, &
					  ufsw_in
real, dimension(:,:,:),intent(in),optional :: dfswcf_in, fswcf_in,  &
					      hswcf_in, ufswcf_in
integer,                  intent(in) ::   jloc
!---------------------------------------------------------------------

    integer     ::  nn, k, iloc

    if (.not. allocated (dfsw) ) then
    allocate ( dfsw( KS:KE+1, num_pts) )
    allocate (  fsw( KS:KE+1, num_pts) )
    allocate (  hsw( KS:KE  , num_pts) )
    allocate ( ufsw( KS:KE+1, num_pts) )

    if (present(dfswcf_in)) then
      allocate ( dfswcf( KS:KE+1, num_pts) )
      allocate (  fswcf( KS:KE+1, num_pts) )
      allocate (  hswcf( KS:KE  , num_pts) )
      allocate ( ufswcf( KS:KE+1, num_pts) )
    endif
    endif

    ssolar = ssolar_in
    do nn = 1, num_pts
      if (jabs(jloc) == jradprt(nn)) then
        if (iradprt(nn) >= iabs(ISRAD) .and.  &
	    iradprt(nn) <= iabs(IERAD)) then
          iloc = iradprt(nn) - iabs(ISRAD) + 1
          do k=KS,KE+1
            dfsw(k,nn) = dfsw_in(iloc, jloc, k)
            fsw(k,nn) =  fsw_in(iloc, jloc, k)
            ufsw(k,nn) = ufsw_in(iloc, jloc, k)
          end do
          do k=KS,KE
            hsw(k,nn) =  hsw_in(iloc, jloc, k)
          end do
          if (present(dfswcf_in)) then
            do k=KS,KE+1
              dfswcf(k,nn) = dfswcf_in(iloc, jloc, k)
              fswcf(k,nn) =  fswcf_in(iloc, jloc, k)
              ufswcf(k,nn) = ufswcf_in(iloc, jloc, k)
            end do
            do k=KS,KE
              hswcf(k,nn) =  hswcf_in(iloc, jloc, k)
            end do
	  endif
        endif
      endif
    end do

!---------------------------------------------------------------------

end subroutine radiag_from_sw_driver 




!###################################################################


subroutine radiag_from_clouds_lh(jloc, camtsw_in, cirabsw_in,    &
				 cirrfsw_in, &
	    		         cmxolw_in, crndlw_in, cvisrfsw_in, &
			         emmxolw_in, emrndlw_in, ncldsw_in, &
			         nmxolw_in, nrndlw_in, ktopsw_in, &
			         kbtmsw_in) 

!--------------------------------------------------------------------
real, dimension(:,:,:,:), intent(in) :: cirabsw_in, cirrfsw_in, &
					cvisrfsw_in  
real, dimension(:,:,:,:), intent(in) :: emmxolw_in, emrndlw_in
real, dimension(:,:,:), intent(in)   :: cmxolw_in, crndlw_in
real, dimension(:,:,:), intent(in)   :: camtsw_in
integer, dimension(:,:),intent(in)   :: ncldsw_in, nmxolw_in, nrndlw_in
integer,                intent(in)   :: jloc
integer, dimension(:,:,:),intent(in) :: ktopsw_in, kbtmsw_in           
!--------------------------------------------------------------------

  integer  ::  k, nn, iloc, nlwcldb, KCAMTL, KCAMTU, nsolwg_loc, KCSW

!--------------------------------------------------------------------
  cld_flg = .true.
  nlwcldb = size(emmxolw_in,4)
  KCAMTL = lbound(camtsw_in,3)
  KCAMTU = ubound(camtsw_in,3)
  nsolwg_loc = size(cvisrfsw_in,4)
  KCSW = ubound(ktopsw_in,3)

!--------------------------------------------------------------------
!allocate arrays to receive desired data from cloud_package     
!--------------------------------------------------------------------
  if (.not. allocated (camtsw) ) then
  allocate ( camtsw (KCAMTL:KCAMTU, num_pts) )
  allocate ( cirabsw(1:KCSW , nsolwg_loc, num_pts) )
  allocate ( cirrfsw(1:KCSW, nsolwg_loc,   num_pts) )
  allocate ( cmxolw (KS:KE,   num_pts) )
  allocate ( crndlw (KS:KE,   num_pts) )
  allocate ( cvisrfsw(1:KCSW ,nsolwg_loc,   num_pts) )
  allocate ( emmxolw(KS:KE, nlwcldb,  num_pts) )
  allocate ( emrndlw(KS:KE, nlwcldb,num_pts) )
  allocate ( ncldsw (         num_pts) )
  allocate ( nmxolw (         num_pts) )
  allocate ( nrndlw (         num_pts) )
  allocate (ktopsw  ( 1:KCSW, num_pts) )
  allocate (kbtmsw  ( 1:KCSW+1, num_pts) )
  endif

!---------------------------------------------------------------------
!load desired data
!---------------------------------------------------------------------
  do nn=1, num_pts
    if (jabs(jloc) == jradprt(nn)) then
      if ( (iradprt(nn) >= iabs(ISRAD)) .and.   &
	   (iradprt(nn) <= iabs(IERAD)) ) then
	iloc = iradprt(nn) - iabs(ISRAD) + 1
        do k=KCAMTL,KCAMTU
          camtsw(k,nn) = camtsw_in(iloc,jloc,k)
        end do
        do k=1,KCSW         
	  cirabsw(k,:,nn) = cirabsw_in(iloc, jloc,k,:)
	  cirrfsw(k,:,nn) = cirrfsw_in(iloc, jloc,k,:)
	  cvisrfsw(k,:,nn) = cvisrfsw_in(iloc, jloc,k,:)
	  ktopsw(k,nn) = ktopsw_in(iloc, jloc, k)
	  kbtmsw(k,nn) = kbtmsw_in(iloc, jloc, k)
        end do
	do k=KS,KE
	  cmxolw(k,nn) = cmxolw_in(iloc, jloc, k)
	  crndlw(k,nn) = crndlw_in(iloc, jloc, k)
        end do
	do k= KS,KE
	  emmxolw(k,:,nn) = emmxolw_in(iloc, jloc, k, :)
	  emrndlw(k,:,nn) = emrndlw_in(iloc, jloc, k, :)
	end do
	ncldsw(nn) = ncldsw_in(iloc, jloc)
	nmxolw(nn) = nmxolw_in(iloc, jloc)
	nrndlw(nn) = nrndlw_in(iloc, jloc)
      endif
    endif
  end do

!--------------------------------------------------------------------
   

end subroutine radiag_from_clouds_lh




!###################################################################

subroutine radiag_from_clouds_esf(jloc, camtsw_in,                &
	    		         cmxolw_in, crndlw_in,              &
			         emmxolw_in, emrndlw_in, ncldsw_in, &
			         nmxolw_in, nrndlw_in)          

!--------------------------------------------------------------------
real, dimension(:,:,:,:), intent(in) :: emmxolw_in,  emrndlw_in
real, dimension(:,:,:), intent(in)   :: cmxolw_in, crndlw_in
real, dimension(:,:,:), intent(in)   :: camtsw_in                 
integer, dimension(:,:),intent(in)   :: ncldsw_in, nmxolw_in, nrndlw_in
integer,                intent(in)   :: jloc
!--------------------------------------------------------------------

  integer    ::  k, nn, iloc, nlwcldb, KCAMTL, KCAMTU

!--------------------------------------------------------------------
  cld_flg=.true.
  nlwcldb = size(emmxolw_in,4)
  KCAMTL = lbound(camtsw_in,3)
  KCAMTU = ubound(camtsw_in,3)

!--------------------------------------------------------------------
!allocate arrays to receive desired data from cloud_package     
!--------------------------------------------------------------------
  if (.not. allocated (camtsw) ) then
  allocate ( camtsw (KCAMTL:KCAMTU, num_pts) )
  allocate ( cmxolw (KS:KE,   num_pts) )
  allocate ( crndlw (KS:KE,   num_pts) )
  allocate ( emmxolw(KS:KE, nlwcldb,  num_pts) )
  allocate ( emrndlw(KS:KE, nlwcldb,num_pts) )
  allocate ( ncldsw (         num_pts) )
  allocate ( nmxolw (         num_pts) )
  allocate ( nrndlw (         num_pts) )
  endif

!---------------------------------------------------------------------
!load desired data
!---------------------------------------------------------------------
  do nn=1, num_pts
    if (jabs(jloc) == jradprt(nn)) then
      if ( (iradprt(nn) >= iabs(ISRAD)) .and.   &
	   (iradprt(nn) <= iabs(IERAD)) ) then
	iloc = iradprt(nn) - iabs(ISRAD) + 1
        do k=KCAMTL,KCAMTU
          camtsw(k,nn) = camtsw_in(iloc,jloc,k)
        end do
	do k=KS,KE
	  cmxolw(k,nn) = cmxolw_in(iloc, jloc, k)
	  crndlw(k,nn) = crndlw_in(iloc, jloc, k)
        end do
	do k= KS,KE
	  emmxolw(k,:,nn) = emmxolw_in(iloc, jloc, k, :)
	  emrndlw(k,:,nn) = emrndlw_in(iloc, jloc, k, :)
	end do
	ncldsw(nn) = ncldsw_in(iloc, jloc)
	nmxolw(nn) = nmxolw_in(iloc, jloc)
	nrndlw(nn) = nrndlw_in(iloc, jloc)
      endif
    endif
  end do

!--------------------------------------------------------------------
   

end subroutine radiag_from_clouds_esf

  



!###################################################################


subroutine radiag_from_astronomy(cosang_in, fracday_in, nsolwg_in, &
 			         jrow, jloc, icolst, icolend)

!--------------------------------------------------------------------
real, dimension(:,:,:), intent(in) :: cosang_in
real, dimension(:,:),   intent(in) :: fracday_in
integer,                intent(in) :: jrow, jloc, icolst, icolend, &
				      nsolwg_in
!--------------------------------------------------------------------

   integer    ::  nn, iloc

!--------------------------------------------------------------------
!allocate arrays to receive desired data from cloud_package     
!--------------------------------------------------------------------
   if (.not. allocated (cosang) ) then
   allocate ( cosang (nsolwg_in , num_pts) )
   allocate (fracday (num_pts) )
   endif

!---------------------------------------------------------------------
!load desired data
!---------------------------------------------------------------------
   do nn=1, num_pts
     if (jrow == jradprt(nn)) then
       if ( (iradprt(nn) >= icolst) .and.    &
	    (iradprt(nn) <= icolend) ) then
         iloc = iradprt(nn) - icolst + 1
         cosang(:, nn) = cosang_in(iloc, jloc, :)
         fracday(nn) = fracday_in (iloc, jloc)
       endif
     endif
   end do


!--------------------------------------------------------------------
   

end subroutine radiag_from_astronomy


!###################################################################


subroutine radiag_from_sfcalbedo(jloc, cvisrfgd_in, cirrfgd_in)

!--------------------------------------------------------------------
real, dimension(:,:),   intent(in) :: cvisrfgd_in, cirrfgd_in
integer,                intent(in) :: jloc
!--------------------------------------------------------------------

    integer    ::  nn, iloc

!--------------------------------------------------------------------
! allocate arrays to receive desired data from cloud_package     
!--------------------------------------------------------------------

    if (.not. allocated (cvisrfgd) ) then
      allocate (cvisrfgd(num_pts) )
      allocate (cirrfgd (num_pts) )
    endif

!---------------------------------------------------------------------
! load desired data
!---------------------------------------------------------------------

    do nn=1, num_pts
      if (jabs(jloc) == jradprt(nn)) then
        if ( (iradprt(nn) >= iabs(ISRAD)) .and.   &
	     (iradprt(nn) <= iabs(IERAD)) ) then
	  iloc = iradprt(nn) - iabs(ISRAD) + 1
	  cvisrfgd(nn) = cvisrfgd_in (iloc, jloc)
	  cirrfgd(nn) = cirrfgd_in (iloc, jloc)
        endif
      endif
    end do

!--------------------------------------------------------------------
   

end subroutine radiag_from_sfcalbedo


!###################################################################

subroutine radiag_from_ozone (jloc, qo3_in)

!--------------------------------------------------------------------
real, dimension(:,:,:),   intent(in) :: qo3_in
integer,                  intent(in) :: jloc
!--------------------------------------------------------------------

    integer    ::  nn, iloc

!--------------------------------------------------------------------
!allocate arrays to receive desired data from ozone_mod     
!--------------------------------------------------------------------
    if (.not. allocated (qo3)) then
      allocate (qo3     ( KS:KE, num_pts) )
    endif

!---------------------------------------------------------------------
!load desired data
!--------------------------------------------------------------------
    do nn=1, num_pts
      if (jabs(jloc) == jradprt(nn)) then
        if ( (iradprt(nn) >= iabs(ISRAD)) .and.   &
	     (iradprt(nn) <= iabs(IERAD))) then
	  iloc = iradprt(nn) - iabs(ISRAD) + 1
	  qo3(:,nn) = qo3_in (iloc, jloc,:)
        endif
      endif
    end do

!--------------------------------------------------------------------
   

end subroutine radiag_from_ozone


!#####################################################################

subroutine radiag_from_cloudrad (jloc, cldext_in, cldsct_in,    &
				 cldasymm_in, deltaz_in,     &
				 lwpath_in, iwpath_in,        &
                                 size_drop_in, size_ice_in)

!--------------------------------------------------------------------
real,    dimension(:,:,:,:), intent(in) ::  cldext_in, cldsct_in,  &
					    cldasymm_in
real,    dimension(:,:,:)  , intent(in) ::  deltaz_in
integer,                     intent(in) ::  jloc
real,  dimension(:,:,:), optional, intent(in) ::                   &
                                              lwpath_in, iwpath_in,  &
                                              size_drop_in, size_ice_in
!--------------------------------------------------------------------

  integer    ::  k, nn, iloc, n
  character(len=10)  :: swform

!--------------------------------------------------------------------
! allocate arrays to receive desired data from cloudrad_package_mod
!--------------------------------------------------------------------
  if (.not. allocated (cldext) ) then
    allocate ( cldext   (KS:KE, nbands, num_pts) )
    allocate ( cldssalb (KS:KE, nbands, num_pts) )
    allocate ( cldasymm (KS:KE, nbands, num_pts) )
    allocate ( deltaz   (KS:KE,         num_pts) )
  endif

  if ( present(lwpath_in) .and. (.not. allocated (lwpath)) ) then
    allocate (lwpath    (KS:KE,          num_pts) )
    allocate (iwpath    (KS:KE,          num_pts) )
    allocate (size_drop (KS:KE,          num_pts) )
    allocate (size_ice  (KS:KE,          num_pts) )
  endif

!---------------------------------------------------------------------
! load desired data : deltaz in meters, cloud extinction optical depth, 
! single scattering albedo, asymmetry parameter for each shortwave band 
!---------------------------------------------------------------------
  swform = Sw_control%sw_form
  if (trim(swform) == 'esfsw99') then
    do nn=1, num_pts
      if (jabs(jloc) == jradprt(nn)) then
        if (iradprt(nn) >= iabs(ISRAD) .and.   &
	    iradprt(nn) <= iabs(IERAD)) then
          iloc = iradprt(nn) - iabs(ISRAD) + 1
          do k=KS,KE
	    deltaz(k,nn) = deltaz_in(iloc,jloc,k)
	  enddo
	  if (allocated (lwpath)) then
            do k=KS,KE
              lwpath(k,nn) = lwpath_in(iloc,jloc,k)
              iwpath(k,nn) = iwpath_in(iloc,jloc,k)
              size_drop(k,nn) = size_drop_in(iloc,jloc,k)
              size_ice(k,nn) = size_ice_in(iloc,jloc,k)
            enddo
          endif
	  do n = 1,nbands
            do k=KS,KE
              cldext(k,n,nn) = cldext_in(iloc,jloc,k,n)*            &
			       deltaz(k,nn) / 1.0E+03
	      if (cldext(k,n,nn) .GT. 0.0) then                      
                cldssalb(k,n,nn) = cldsct_in(iloc,jloc,k,n) /       &
		  		   cldext_in(iloc,jloc,k,n)
              else
	        cldssalb(k,n,nn) = 0.0
              endif
              cldasymm(k,n,nn) = cldasymm_in(iloc,jloc,k,n)
	    enddo
	  enddo
        endif
      endif
    end do
  endif

!--------------------------------------------------------------------
   


end subroutine radiag_from_cloudrad 

 



!##################################################################

subroutine radiag

!-----------------------------------------------------------------------
!     radiag writes out all quantities required for diagnosis of short
!     and longwave radiation quantities.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------

real, dimension(:), allocatable    :: ctst, dfswd, flwsw, flxem, &
                                      flxemch4n2o, flxnetd, fswd,&
 			              ftopac, ftopn, hlwsw, htem,&
			              htem1, htem2, htem3, htem4,&
			              htem5, htem6, ftopef,  &
				      htem7t, ufswd, vsumac,   &
				      dfswdcf, flwswcf,   &
				      flxnetdcf, fswdcf,   &
				      hlwswcf,  ufswdcf
integer, dimension(:), allocatable :: indx_kbtmsw, indx_ktopsw
real, dimension(:,:), allocatable  :: htem7 

integer :: ioffset, k, m, n, nw, no_writes, nbi, nbf, index_k,  &
	   kc, nprt, ngp, ny, nx, NLWCLDB, nn, nd, stdout, jindx,&
	   KCSW
real    :: ftopc, ftope, ftop, fgrd, fdiff, qsum, ftopeft, zenith2
!--------------------------------------------------------------------

do jindx=JSRAD,JERAD
  if ( .not. Rad_control%do_raddg(jabs(jindx))) then
    cycle 
  else
    stdout = radiag_unit
    if (do_lhsw) then
      KCSW = ubound(ktopsw,1)
    endif
    nd = 0
    do nn = 1,num_pts
      if (jradprt(nn) == jabs(jindx) .and. iradprt(nn) /= 0 .and. &
	   iradprt(nn) .ge. iabs(israd) .and.     &
	   iradprt(nn) .le. iabs(ierad) ) then
        nd = nn
        write(stdout,99000) iradprt_gl(nn), jradprt_gl(nn),   &
			    deglon1(nn), deglat1(nn)
        if (lswg) then
          write(stdout,99010) nsolwg
        else if (ldiurn) then
          write(stdout,99020)
        else if (do_annual) then
          write(stdout,99025)
        else
          write(stdout,99030)
        endif

!-------------------------------------------------------------------
!allocate some local arrays
!-------------------------------------------------------------------

        allocate ( ctst        (KS:KE  ))
        allocate ( dfswd       (KS:KE+1))
        allocate ( flwsw       (KS:KE+1))
        allocate ( flxem       (KS:KE+1))
        allocate ( flxemch4n2o (KS:KE+1))
        allocate ( flxnetd     (KS:KE+1))
        allocate ( fswd        (KS:KE+1))
        allocate ( ftopac      (NBLY       ))
        allocate ( ftopn       (NBLY       ))
        allocate ( hlwsw       (KS:KE  ))
        allocate ( htem        (KS:KE  ))
        allocate ( htem1       (KS:KE  ))
        allocate ( htem2       (KS:KE  ))
        allocate ( htem3       (KS:KE  ))
        allocate ( htem4       (KS:KE  ))
        allocate ( htem5       (KS:KE  ))
        allocate ( htem6       (KS:KE  ))
        if (Lw_control%do_ch4_n2o) then
          allocate ( ftopef      (NBTRGE     ))
          allocate ( htem7       (KS:KE, NBTRGE  ))
          allocate ( htem7t      (KS:KE  ))
        endif
        if (do_lhsw) then
          allocate ( indx_kbtmsw (KCSW       ))
          allocate ( indx_ktopsw (KCSW       ))
        endif
        allocate ( ufswd       (KS:KE+1))
        allocate ( vsumac      (NBLY       ))
        if (do_totcld_forcing) then
          allocate ( dfswdcf     (KS:KE+1))
          allocate ( flwswcf     (KS:KE+1))
          allocate ( flxnetdcf   (KS:KE+1))
          allocate ( fswdcf      (KS:KE+1))
          allocate ( hlwswcf     (KS:KE  ))
          allocate ( ufswdcf     (KS:KE+1))
        endif
      
!----------------------------------------------------------------------
!obtain emissivity heating rates and fluxes.
!
!if(not do_ch4_n2o) then
!     htem1 = emissivity heating rate for 0-160, 1200-2200 cm-1 band, 
!else 
!     htem1 = emissivity heating rate for 0-160, 1400-2200 cm-1 band.
!   
!     htem2 = emissivity heating rate for 560-800 cm-1 band.
!     htem3 = emissivity heating rate for 800-900 cm-1 band
!     htem4 = emissivity heating rate for 900-990 cm-1 band
!     htem5 = emissivity heating rate for 990-1070 cm-1 band.
!     htem6 = emissivity heating rate for 1070-1200 cm-1 band
!----------------------------------------------------------------------
        do k=KS,KE
          htem1(k) = radcon*(flx1(k+1,nd) - flx1(k,nd))*  &
                     pdfinv(k,nd)
          htem2(k) = radcon*(flx2(k+1, nd) - flx2(k,nd))*  &
                     pdfinv(k,nd)
          htem3(k) = radcon*(flx3(k+1,nd) - flx3(k,nd))*  &
                     pdfinv(k,nd)
          htem4(k) = radcon*(flx4(k+1,nd) - flx4(k,nd))*   &
                     pdfinv(k,nd)
          htem5(k) = radcon*(flx5(k+1,nd) - flx5(k,nd))*   &
                     pdfinv(k,nd)
          htem6(k) = radcon*(flx6(k+1,nd) - flx6(k,nd))*  &
                     pdfinv(k,nd)
        enddo

!----------------------------------------------------------------------
!     htem7 = emissivity heating rate for 1200-1400 cm-1 range
!             (NBTRGE bands)
!     htem7t = sum of emissivity heating rate over 1200-1400 cm-1 range.
!----------------------------------------------------------------------
        if (Lw_control%do_ch4_n2o) then
          do m=1,NBTRGE
            do k=KS,KE
              htem7(k,m) = radcon*(flx7(k+1,m,nd) - flx7(k,m,nd))* &
                           pdfinv(k,nd)
            enddo
          enddo
          do k=KS,KE
	    htem7t(k) = 0.0E+00
	    do m=1,NBTRGE
	      htem7t(k) = htem7t(k) + htem7(k,m)
	    enddo
          enddo
        endif

!----------------------------------------------------------------------
!     htem = approximate heating rate summed over all frequency ranges.
!----------------------------------------------------------------------
        do k=KS,KE
          htem(k) = htem1(k) + htem2(k) + htem3(k) + htem4(k) +   &
                    htem5(k) + htem6(k)
        enddo
        if (Lw_control%do_ch4_n2o) then
          do m=1,NBTRGE
            do k=KS,KE
  	      htem(k) = htem(k) + htem7(k,m)
            enddo
          enddo
        endif

!----------------------------------------------------------------------
!     compute fluxes at top, ground and net flux at all levels.
!----------------------------------------------------------------------
        ftopc     = gxcts (nd)*1.0E-03
        ftope     = flx1e1(nd)*1.0E-03
        do k=KS,KE+1
          flxem(k) = flx1(k,nd) + flx2(k,nd) +   &
                     flx3(k,nd) + flx4(k,nd) +  &
                     flx5(k,nd) + flx6(k,nd)
        enddo

        if (Lw_control%do_ch4_n2o) then
          ftopeft   = 0.0E+00
          flxemch4n2o(:) = 0.0E+00
          do m=1,NBTRGE
            ftopef(m) = flx1e1f(m,nd)*1.0E-03
            ftopeft   = ftopeft + ftopef(m)
	    do k=KS,KE+1 
	      flxemch4n2o(k) = flxemch4n2o(k) + flx7(k,m,nd)
	    enddo
          enddo
          do k=KS,KE+1
	    flxem(k) = flxem(k) + flxemch4n2o(k)
          enddo
        endif

        ftop      = ftopc + ftope
        fgrd      = flxnet(        KE+1,nd)*1.0E-03
        fdiff     = ftop - fgrd
        do k=KS,KE+1
	  flxnetd(k) = flxnet(k,nd)*1.0E-03
          fswd (k) = fsw (k,nd)*1.0E-03
          dfswd(k) = dfsw(k,nd)*1.0E-03
          ufswd(k) = ufsw(k,nd)*1.0E-03
        enddo
        if (do_totcld_forcing) then
          do k=KS,KE+1
	    flxnetdcf(k) = flxnetcf(k,nd)*1.0E-03
            fswdcf (k) = fswcf (k,nd)*1.0E-03
            dfswdcf(k) = dfswcf(k, nd)*1.0E-03
            ufswdcf(k) = ufswcf(k, nd)*1.0E-03
          enddo
        endif

!----------------------------------------------------------------------
!     compute net radiative heating and net flux (up-down).
!----------------------------------------------------------------------
        do k=KS,KE+1
          flwsw(k) = flxnetd(k) + fswd(k)
        enddo
        do k=KS,KE
          hlwsw(k) = hsw(k,nd) + heatra(k,nd)
        enddo

        if (do_totcld_forcing) then
          do k=KS,KE+1
            flwswcf(k) = flxnetdcf(k) + fswdcf(k)
          enddo
          do k=KS,KE
            hlwswcf(k) = hswcf(k,nd) + heatracf(k,nd)
          enddo
        endif

!----------------------------------------------------------------------
!     compute diagnosis of cool-to-space quantities.     
!----------------------------------------------------------------------
        do n=1,NBLY-1
          qsum = 0.0E+00
          do k=KS,KE
            qsum = qsum + exctsn(k,n,nd)/(pdfinv(k,nd)*radcon)
          enddo
          ftopn(n) = fctsg(n, nd) - qsum
        enddo

!----------------------------------------------------------------------
!     compute 4.3 um band contribution.
!---------------------------------------------------------------------- 
        ftopn(NBLY) = 0.0E+00
        do n=1,NBLY
          ftopn(n) = ftopn(n)*1.0E-03
        enddo
        ftopac(1) = ftopn(1)
        do n=2,NBLY
          ftopac(n) = ftopac(n-1) + ftopn(n)
        enddo

!----------------------------------------------------------------------
!     compute 4.3 um band contribution.
!----------------------------------------------------------------------
        fctsg(NBLY,nd) = 0.0E+00

!----------------------------------------------------------------------
!     set excts contribution of band (NBLY) (4.3um) to zero.
!---------------------------------------------------------------------- 
        do k=KS,KE
          exctsn(k,NBLY,nd) = 0.0E+00
        enddo
        do n=1,NBLY
          fctsg(n,nd) = fctsg(n,nd)*1.0E-03
        enddo
        vsumac(1) = fctsg(1,nd)
        do n=2,NBLY
          vsumac(n) = vsumac(n-1) + fctsg(n,nd)
        enddo

!----------------------------------------------------------------------
!     compute approximate cool-to-space heating rates.     
!----------------------------------------------------------------------
        do k=KS,KE
          ctst(k) = ctso3(k,nd) + ctsco2(k,nd) + cts(k,nd)
        enddo

!----------------------------------------------------------------------
!     write point, cloud data (if any), and surface albedo.
!----------------------------------------------------------------------
        if (cld_flg)  then
          write(stdout,9010) nmxolw(nd), nrndlw(nd)
          if (nmxolw(nd) .GT. 0 .OR. nrndlw(nd) .GT. 0) then
            NLWCLDB = size (emmxolw,2)
            if (NLWCLDB == 1) then
              write (stdout,9020)
              do k = KS,KE
                if (cmxolw(k,nd) .GT. 0.0   .or.   &
                    crndlw(k,nd) .GT. 0.0        )  then
	          index_k = k - KS + 1
	          write (stdout,9030)   &
                        index_k, cmxolw(k,nd), emmxolw(k,1,nd),  &
                        crndlw(k,nd), emrndlw(k,1,nd)
	        endif
              enddo
            else
              no_writes = (NLWCLDB - 1)/8 + 1
              do nw = 1,no_writes
                if (nw .lt. no_writes) then
  	          nbi = (nw-1)*8 + 1
	          nbf = (nw-1)*8 + 8
                else
  	          nbi = (nw-1)*8 + 1
	          nbf = NLWCLDB
                endif
                write (stdout,9022) (n,n=nbi,nbf)
                do k = KS,KE
	          if (cmxolw(k,nd) .GT. 0.0 )     then
                    index_k = k - KS + 1
	            write (stdout,9032)    &
                        index_k, cmxolw(k,nd), (emmxolw(k,n,nd),  &
                                                 n=nbi,nbf)
	          endif
                enddo
              enddo
              do nw = 1,no_writes
                if (nw .lt. no_writes) then
  	          nbi = (nw-1)*8 + 1
	          nbf = (nw-1)*8 + 8
                else
  	          nbi = (nw-1)*8 + 1
	          nbf = NLWCLDB
                endif
                write (stdout,9024) (n,n=nbi,nbf)
                do k = KS,KE
	          if (crndlw(k,nd) .GT. 0.0 )     then
     	            index_k = k - KS + 1
	            write (stdout,9032)    &
                           index_k, crndlw(k,nd), (emrndlw(k,n,nd),   &
                                                    n=nbi,nbf)
	          endif
                enddo
              enddo
            endif 
          endif 

!--------------------------------------------------------------------
          if (do_lhsw) then
            if (ncldsw(nd) .NE. 0) then
              do kc = ncldsw(nd),1,-1
                indx_ktopsw(kc) = ktopsw(kc, nd) - (KS - 1)
                indx_kbtmsw(kc) = kbtmsw(kc,nd) - (KS - 1)
              enddo
              write(stdout,9035) 
              write(stdout,9036)    &
                          (kc, camtsw (kc, nd), indx_ktopsw(kc),   &
                           indx_kbtmsw(kc)   , cvisrfsw(kc,1,nd),&
                           cirrfsw(kc,1, nd), cirabsw (kc,1,nd),&
                                             kc=ncldsw(nd),1,-1)
            endif
          else if (do_esfsw) then
            if (ncldsw(nd) .NE. 0) then
	      if (allocated (lwpath)) then
                write (stdout, 9510)
                do k=KS,KE
                  if (camtsw(k,nd) .GT. 0.0) then
                    write (stdout, 9520)                              &
                     k,lwpath(k,nd), iwpath(k,nd), size_drop(k,nd),  &
                     size_ice(k,nd)
                  endif
                enddo
              endif
              do n=1,nbands
                write(stdout,9040) n
                do k=KS,KE
                  if (camtsw(k,nd) .GT. 0.0) then
                    write (stdout,9050)                            &
                            k, camtsw (k,nd), cldext(k,n,nd),       &
                            cldssalb(k,n,nd),                       &
                            cldasymm(k,n,nd)
                  endif
                enddo
              enddo
            endif
          else
	    call error_mesg ('radiation_diag', &
                           'no shortwave clouds are activated', FATAL)
          endif

!----------------------------------------------------------------------
! if there are no clouds in the column, write a message
!----------------------------------------------------------------------
        else   ! (cld_flg = .false., i.e., radiag is being called w/o 
	       !  clouds, or no clouds are present in this column or 
	       !  chunk)
          write (stdout, 9052)
        endif

        write(stdout,9060) cvisrfgd(nd), cirrfgd(nd)

!----------------------------------------------------------------------
!     write co2 amount.     
!--------------------------------------------------------------------
        write(stdout,9070) rrvco2

!--------------------------------------------------------------------
!     write cfc amounts
!--------------------------------------------------------------------
        write (stdout,9071) rrvf11
        write (stdout,9072) rrvf12
        write (stdout,9075) rrvf113
        write (stdout,9076) rrvf22

!---------------------------------------------------------------------
!     write ch4 and n2o amounts
!--------------------------------------------------------------------
        write (stdout,9073) rrvch4
        write (stdout,9074) rrvn2o
 
!----------------------------------------------------------------------
!     write solar input, zenith angle, and day fraction.      
!----------------------------------------------------------------------
        zenith2 = 0.0E+00 
        do ngp=1,nsolwg
          zenith2 = zenith2 + gwt(ngp)*cosang(ngp,nd)
        end do
        write(stdout,9080) ssolar, zenith2, fracday(nd)

!----------------------------------------------------------------------
!     write input data and overall heating rates and fluxes.
!----------------------------------------------------------------------
        write(stdout,9090)
        write(stdout,9100)
        write(stdout,9110) (k-KS+1, press (    k,nd), temp  (k,nd),&
                                    rh2o  (k,nd), qo3   (k,nd),&
                                    heatra(k,nd), flxnetd(k), &
                                    pflux (        k,nd), k=KS,KE)

        write(stdout,9120) press (KE+1,nd), &
                           temp  (KE+1,nd), &
                           flxnetd(KE+1), pflux(        KE+1,nd)
        write(stdout,9130)
        write(stdout,9140)
        write(stdout,9150) (k-KS+1, press(k,nd), hsw(k,nd), &
                                    fswd(k), dfswd(k), ufswd(k),&
                                    pflux(        k,nd), k=KS,KE)

        write(stdout,6556) press(KE+1,nd), fswd(KE+1), dfswd(KE+1), &
                           ufswd(KE+1), pflux(        KE+1,nd)
        write(stdout,9160)
        write(stdout,9170)
        write(stdout,9190) (k-KS+1, press(k,nd), hlwsw(k), flwsw(k),&
                                    pflux(        k,nd), k=KS,KE)

        write(stdout,9180) press(KE+1,nd), flwsw(KE+1),   &
                           pflux(        KE+1,nd)

        if (do_totcld_forcing) then

          write(stdout,9400)
          write(stdout,9100)
          write(stdout,9110) (k-KS+1, press (k,nd), temp  (k,nd),&
                                      rh2o  (k,nd), qo3   (k,nd),&
                                      heatracf(k,nd), flxnetdcf(k), &
                                      pflux (        k,nd), k=KS,KE)

          write(stdout,9120) press (KE+1,nd),  &
                             temp  (KE+1,nd), &
                             flxnetdcf(KE+1), pflux(        KE+1,nd)
          write(stdout,9410)
          write(stdout,9140)
          write(stdout,9150) (k-KS+1, press(k,nd), hswcf(k,nd), &
                                      fswdcf(k), dfswdcf(k),   &
				      ufswdcf(k), pflux(k,nd), k=KS,KE)

          write(stdout,6556) press(KE+1,nd), fswdcf(KE+1),  &
                             dfswdcf(KE+1),ufswdcf(KE+1), pflux(KE+1,nd)
          write(stdout,9420)
          write(stdout,9170)
          write(stdout,9190) (k-KS+1, press(k,nd), hlwswcf(k),  &
                              flwswcf(k), pflux(        k,nd), k=KS,KE)

          write(stdout,9180) press(KE+1,nd), flwswcf(KE+1),  &
                             pflux(        KE+1,nd)
        endif

!----------------------------------------------------------------------
!     write approximate heating rates.
!----------------------------------------------------------------------
        if (.not. Lw_control%do_ch4_n2o) then
          write(stdout,9200)
          write(stdout,9210) (k-KS+1, press(k,nd), htem1(k), htem2(k),&
                                      htem3(k), htem4(k), htem5(k),  &
		 		      htem6(k), htem(k), k=KS,KE)
        else
          if (NBTRGE .EQ. 1) then
            write(stdout,9201)
            write(stdout,9211) (k-KS+1, press(k,nd), htem1(k),  &
                                        htem2(k), htem3(k), htem4(k), &
                                        htem5(k), htem6(k), &
                                        (htem7(k,n),n=1,NBTRGE), &
                                        htem(k),k=KS,KE)
          elseif (NBTRGE .EQ. 2) then
            write(stdout,9202)
            write(stdout,9212) (k-KS+1, press(k,nd), htem1(k),  &
                                        htem2(k), htem3(k), htem4(k), &
                                        htem5(k), htem6(k),   &
		  		        (htem7(k,n),n=1,NBTRGE),  &
				        htem7t(k), htem(k),k=KS,KE)
          elseif (NBTRGE .EQ. 4) then
            write(stdout,9203)
            write(stdout,9213) (k-KS+1, press(k,nd), htem1(k),   &
                                        htem2(k), htem3(k), htem4(k),&
                                        htem5(k), htem6(k),   &
                                        (htem7(k,n),n=1,NBTRGE),   &
		  		        htem7t(k), htem(k),k=KS,KE)
          endif
        endif

        write(stdout,9220)
        write(stdout,9230) (k-KS+1, press(k,nd), cts(k,nd),  &
                                    ctsco2(k,nd), ctso3(k,nd), &
                                    ctst(k), k=KS,KE)
        if (NBLY == 16) then
          write(stdout,9240)
          write(stdout,9250) (k-KS+1, press(k,nd), excts(k,nd),  &
                              exctsn(k,1,nd), exctsn(k,2,nd), &
                              exctsn(k,3,nd), exctsn(k,4,nd), &
                              exctsn(k,5,nd), exctsn(k,6,nd), &
                              exctsn(k,7,nd), k=KS,KE)
          write(stdout,9260)
          write(stdout,9250) (k-KS+1, press(k,nd),   &
                              exctsn(k,8,nd), exctsn(     k,9,nd) , &
                              exctsn(k,10,nd), exctsn(k,11,nd), &
                              exctsn(k,12,nd), exctsn(k,13,nd), &
                              exctsn(k,14,nd), exctsn(k,15,nd), &
                              k=KS,KE)
        else if (NBLY == 48) then
          write(stdout,9240)
          write(stdout,9250) (k-KS+1, press(k,nd), excts(k, nd), &
                              exctsn(k,1,nd), exctsn(k,2,nd), &
                              exctsn(k,3,nd), exctsn(k,4,nd), &
                              exctsn(k,5,nd), exctsn(k,6,nd), &
                              exctsn(k,7,nd), k=KS,KE)
          write(stdout,9260)
          write(stdout,9250) (k-KS+1, press(k,nd),  &
                              exctsn(k,8,nd), exctsn(        k,9,nd) , &
                              exctsn(k,10,nd), exctsn(k,11,nd), &
                              exctsn(k,12,nd), exctsn(k,13,nd), &
                              exctsn(k,14,nd), exctsn(k,15,nd), &
                              k=KS,KE)
          write(stdout,9261)
          write(stdout,9250) (k-KS+1, press(k,nd),   &
                              exctsn(k,16,nd), exctsn(k,17,nd) , &
                              exctsn(k,18,nd), exctsn(k,19,nd), &
                              exctsn(k,20,nd), exctsn(k,21,nd), &
                              exctsn(k,22,nd), exctsn(k,23,nd), &
                              k=KS,KE)
          write(stdout,9262)
          write(stdout,9250) (k-KS+1, press(k,nd),   &
                              exctsn(k,24,nd), exctsn(k,25,nd) , &
                              exctsn(k,26,nd), exctsn(k,27,nd), &
                              exctsn(k,28,nd), exctsn(k,29,nd), &
                              exctsn(k,30,nd), exctsn(k,31,nd), &
                              k=KS,KE)
          write(stdout,9263)
          write(stdout,9250) (k-KS+1, press(k,nd),   &
                              exctsn(k,32,nd), exctsn(k,33,nd) , &
                              exctsn(k,34,nd), exctsn(k,35,nd), &
                              exctsn(k,36,nd), exctsn(k,37,nd), &
                              exctsn(k,38,nd), exctsn(k,39,nd), &
                              k=KS,KE)
          write(stdout,9264)
          write(stdout,9250) (k-KS+1, press(k,nd),   &
                              exctsn(k,40,nd), exctsn(k,41,nd) , &
                              exctsn(k,42,nd), exctsn(k,43,nd), &
                              exctsn(k,44,nd), exctsn(k,45,nd), &
                              exctsn(k,46,nd), exctsn(k,47,nd), &
                              k=KS,KE)
        else
          call error_mesg ( 'radiag',  &
                  'printout is set up only for NBLY=16 or NBLY=48', &  
								 FATAL)
        endif

!----------------------------------------------------------------------
!     write fluxes.
!----------------------------------------------------------------------
        write(stdout,9270) ftopc, ftope, ftop, fgrd, fdiff
        if (Lw_control%do_ch4_n2o) then
          do m=1,NBTRGE
            write (stdout,9271) m,ftopef(m)
          enddo
          write (stdout,9272) ftopeft
        endif
        write(stdout,9280)

!----------------------------------------------------------------------
!#ifndef ckd2p1
!     write out for 8 combined bands 160-560 cm-1.
!#else   ckd2p1
!     write out for 40 combined bands 160-560 cm-1.
!#endif  ckd2p1
!----------------------------------------------------------------------
        if (NBLY == 48) then
          ioffset = 32
        else if (NBLY == 16) then
          ioffset = 0
        else
          call error_mesg ('radiag', &
                    ' NBLY must be either 48 or 16', FATAL)
        endif
        
        do ny = 1, 8+ioffset
          nprt = 1
          do nx=1,40
            if (iband(nx) .EQ. ny) then
              if (nprt .EQ. 1) then
                write(stdout,9290) ny, bandlo(nx+16),bandhi(nx+16),  &
                                   ftopn(ny),ftopac(ny),fctsg(ny,nd),&
                                   vsumac(ny)
                nprt = 0
              else
                write(stdout,9300) bandlo(nx+16), bandhi(nx+16)
              endif
            endif
          enddo
        enddo

!----------------------------------------------------------------------
!     write out for remaining bands.
!----------------------------------------------------------------------
        do ny =9+ioffset, NBLY
          write(stdout,9290) ny, bdlocm(ny), bdhicm(ny), ftopn(ny),  &
                             ftopac(ny), fctsg(ny,nd), vsumac(ny)
        enddo

!----------------------------------------------------------------------
!     write out emissivity fluxes
!----------------------------------------------------------------------
        write (stdout,9310)
        if (.not. Lw_control%do_ch4_n2o) then
          write (stdout,9320) 
          write (stdout,9330) (k-KS+1, flx1(k,nd),flx2(k,nd),&
                                       flx3(k,nd), flx4(k,nd), &
                                       flx5(k,nd), flx6(k,nd), &
                                       flxem(k),  k=KS,KE+1)
        else
          if (NBTRGE .EQ. 1) then
            write (stdout,9321)
            write (stdout,9331) (k-KS+1, flx1(k,nd),flx2(k,nd), &
                                         flx3(k,nd),flx4(k,nd), &
                                         flx5(k,nd), flx6(k,nd), &
                                         flxemch4n2o(k), flxem(k), &
                                         k=KS,KE+1)
          elseif (NBTRGE .EQ. 2) then
            write (stdout,9322)
            write (stdout,9332) (k-KS+1, flx1(k,nd),flx2(k,nd), &
                                         flx3(k,nd),flx4(k,nd), &
                                         flx5(k,nd), flx6(k,nd), &
                                         (flx7(k,m,nd),m=1,NBTRGE), &
                                         flxemch4n2o(k), flxem(k), &
                                         k=KS,KE+1)
          elseif (NBTRGE .EQ. 4) then
            write (stdout,9323)
            write (stdout,9333) (k-KS+1, flx1(k,nd),flx2(k,nd), &
                                         flx3(k,nd),flx4(k,nd), &
                                         flx5(k,nd), flx6(k,nd),&
                                        (flx7(k,m,nd),m=1,NBTRGE),&
                                         flxemch4n2o(k), flxem(k), &
                                         k=KS,KE+1)
          endif
        endif

!-----------------------------------------------------------------
!deallocate local arrays
!-----------------------------------------------------------------

        if (Rad_control%do_totcld_forcing) then
          deallocate ( ufswdcf     )
          deallocate ( hlwswcf     )
          deallocate ( fswdcf      )
          deallocate ( flxnetdcf   )
          deallocate ( flwswcf     )
          deallocate ( dfswdcf     )
        endif
        deallocate ( vsumac      )
        deallocate ( ufswd       )
        if (do_lhsw) then
          deallocate ( indx_ktopsw )
          deallocate ( indx_kbtmsw )
        endif
        if (Lw_control%do_ch4_n2o) then
          deallocate ( htem7t      )
          deallocate ( htem7       )
          deallocate ( ftopef      )
        endif
        deallocate ( htem6       )
        deallocate ( htem5       )
        deallocate ( htem4       )
        deallocate ( htem3       )
        deallocate ( htem2       )
        deallocate ( htem1       )
        deallocate ( htem        )
        deallocate ( hlwsw       )
        deallocate ( ftopn       )
        deallocate ( ftopac      )
        deallocate ( fswd        )
        deallocate ( flxnetd     )
        deallocate ( flxemch4n2o )
        deallocate ( flxem       )
        deallocate ( flwsw       )
        deallocate (dfswd)
        deallocate (ctst)
      endif
    end do
    call radiation_diag_dealloc 
  endif
end do

cld_flg=.false.

!----------------------------------------------------------------------
!     format statements.
!----------------------------------------------------------------------

99000  format (' COLUMN = ', I3, 2X, ' ROW = ', I3, &
       ' LON = ', F10.5, 2X, ' LAT = ', F10.5)
99010  format (' SHORTWAVE RESULTS ', I3, ' GAUSSIAN ZENITH ANGLES ')
99020  format (' SHORTWAVE RESULTS, DIURNALLY VARYING ZENITH ANGLES')
99025  format (' SHORTWAVE RESULTS, ANNUAL MEAN ZENITH ANGLES')
99030  format (' SHORTWAVE RESULTS, DIURNALLY AVERAGED ZENITH ANGLES ')
9000  format (//,'  RADIATION RESULTS FOR IP= ',I3)
9010  format (/,' NO. MAX OVERLAP CLOUDS= ',I2, ' NO. RANDOM OVERLAP',&
       ' CLOUDS= ',I2)
9020  format (37X,' LW CLOUD DATA'/,' LEVEL',8X,  &
       'MAX OVERLAP CLD. AMT.',2X,'MAX OVERLAP CLD. EMISS.',2X, &
       'RND OVERLAP CLD. AMT.',2X,'RND OVERLAP CLD. EMISS.')
9021  format (37X,' SW CLOUD DATA'/,' LEVEL',8X,  &
       '            CLD. AMT.')
9022  format (37X,' LW MAX OVERLAP CLOUD DATA'/,' LEVEL',8X, &
       ' CLOUD AMOUNT',20X,'CLOUD EMISSIVITY (BANDS)',/29x,8(i2,6x))
9024  format (37X,' LW RND OVERLAP CLOUD DATA'/,' LEVEL',8X, &
       ' CLOUD AMOUNT',20X,'CLOUD EMISSIVITY (BANDS)',/29x,8(i2,6x))
9030  format (I4,7X,4F20.6)
9032  format (i4,7x,f12.6,2x,8f8.4)
9035  format (37X,' SW CLOUD DATA'/,' CLD. NO',8X,'CLD. AMT.',2X, &
       'CLD TOP INDEX',2X,'CLD BOT INDEX',2X,'VIS. REFL',3X, &
       ' IR REFL',4X,' IR ABS.')  
9036  format (I5,7X,F12.6,I8,I15,6X,3F12.6)
9040  format (27X,' SW CLOUD DATA, BAND = ',i2,/,' MDL. LVL',7X, &
      'CLD. AMT.',7X,'EXT OP DEP.',3X, &
      ' SSALB.',2X,' ASYMM. PAR.')
9050  format (I5,7X,F12.6,6X,3F12.6)
9052  format (/, ' THIS RADIATION CALL SEES NO CLOUDS.')
9060  format (/,10X,'VIS. SFC. ALBEDO=',F12.6,' IR SFC. ALBEDO=',  &
              F12.6)
9070  format (/,' CO2 VOL. MIXING RATIO= ',F14.6,/)
9071  format (/,' F11 VOL. MIXING RATIO= ',12PF14.2,' pptv',/)
9072  format (/,' F12 VOL. MIXING RATIO= ',12PF14.2,' pptv',/)
9075  format (/,' F113 VOL. MIXING RATIO= ',12PF14.2,' pptv',/)
9076  format (/,' F22 VOL. MIXING RATIO= ',12PF14.2,' pptv',/)
9073  format (/,' CH4 VOL. MIXING RATIO= ',9PF14.2,' ppbv',/)
9074  format (/,' N2O VOL. MIXING RATIO= ',9PF14.2,' ppbv',/)
9080  format (/,' INCOMING SOLAR FLUX =',F12.6,' W/M**2',/,  &
       ' COS(AZIMUTH)=',F12.6,10X,' FRACTION SUNUP=',F12.6)
9090  format (/,20X,' LW HEATING RATES AND FLUXES',/)
9100  format ('  LVL',' PRESSURE   ',4X,' TEMP.     ','H2O MMR',5X,&
       'O3 MMR',7X,'HEAT RATE',2X,'NET FLUX',3X,'FLUX PRESS.')
9110  format (I4,E13.6,F12.4,2E12.5,2F12.6,E13.6)
9120  format (4X,E13.6,F12.4,36X,F12.6,E13.6)
9130  format (/,20X,' SW HEATING RATES AND FLUXES',/)
9140  format ('  LVL',' PRESSURE    ',3X,'HEAT RATE',2X,'NET FLUX',  &
       4X,'DN FLUX',6X,'UP FLUX',3X,'FLUX PRESS.') 
6556  format (4X,E13.6,12X,3F12.6,E13.6)
9150  format (I4,E13.6,4F12.6,E13.6)
9160  format (/,20X,' COMBINED HEATING RATES AND FLUXES',/)
9170  format ('  LVL',' PRESSURE    ',4X,'HEAT RATE',2X,'NET FLUX',  &
       3X,'FLUX PRESS.') 
9180  format (4X,E13.6,12X,F12.6,E13.6)
9190  format (I4,E13.6,2F12.6,E13.6)
9200  format (/,37X,'APPROXIMATE HEATING RATES        (Q(APPROX))'/  &
       '  LVL',' PRESSURE   ',5X,'  0-160,1200-2200 ', &
      '    560-800       ','     800-900      ','      900-990     ', &
      '    990-1070      ','     1070-1200    ',  &
       '       TOTAL')
9210  format (I4,E13.6,7F18.6)
9201  format (/,37X,'APPROXIMATE HEATING RATES        (Q(APPROX))'/  &
       '  LVL',' PRESSURE   ',5X,'  0-160,1400-2200 ',  &
      '    560-800       ','     800-900      ','      900-990     ',&
      '    990-1070      ','     1070-1200    ','     1200-1400    ', &
       '       TOTAL')
9211  format (I4,E13.6,8F18.6)
9202  format (/,37X,'APPROXIMATE HEATING RATES        (Q(APPROX))'/  &
       '  LVL',' PRESSURE   ',5X,'  0-160,1400-2200 ',   &
      '    560-800       ','     800-900      ','      900-990     ', &
      '    990-1070      ','     1070-1200    ','  1200-1300 ', &
       '  1300-1400 ','  1200-1400 ','       TOTAL')
9212  format (I4,E13.6,6F18.6,3F12.6,F18.6)
9203  format (/,37X,'APPROXIMATE HEATING RATES        (Q(APPROX))'/   &
       '  LVL',' PRESSURE   ',5X,'  0-160,1400-2200 ',  &
      '    560-800       ','     800-900      ','      900-990     ', &
      '    990-1070      ','     1070-1200    ','  1200-1250 ',  &
       '  1250-1300 ','  1300-1350 ','  1350-1400 ', &
       '  1200-1400 ','       TOTAL')
9213  format (I4,E13.6,6F18.6,5F12.6,F18.6)
9220  format(/,37X,'APPROXIMATE CTS HEATING RATES'/'  LVL',' PRESSURE',&
       7X,' H2O BANDS    ',' 15 UM BAND   ',' 9.6 UM BAND  ',' TOTAL')
9230  format (I4,E13.6,4F14.6)
9240  format (/,37X,'EXACT CTS HEATING RATES, BY BAND'/   &
       '  LVL',' PRESSURE   ','    TOTAL    ',5X,'1',11X,'2',11X,'3',  &
         11X,'4',11X,'5',11X,'6',11X,'7',/)
9250  format (I4,E13.6,8F12.6)
9260  format ('  LVL',' PRESSURE   ',7X,'8',11X,'9',10X,'10',10X,'11', &
           10X,'12',10X,'13',10X,'14',10X,'15')
9261  format('  LVL',' PRESSURE   ',6X,'16',10X,'17',10X,'18',10X,'19',&
           10X,'20',10X,'21',10X,'22',10X,'23')
9262  format('  LVL',' PRESSURE   ',6X,'24',10X,'25',10X,'26',10X,'27',&
           10X,'28',10X,'29',10X,'30',10X,'31')
9263  format('  LVL',' PRESSURE   ',6X,'32',10X,'33',10X,'34',10X,'35',&
           10X,'36',10X,'37',10X,'38',10X,'39')
9264  format('  LVL',' PRESSURE   ',6X,'40',10X,'41',10X,'42',10X,'43',&
           10X,'44',10X,'45',10X,'46',10X,'47')
9270  format ( 40X,'   FLUXES'/   &
               ' FLUX AT TOP,160-1200 CM-1       =',F14.6,' W/M**2'/ &
               ' FLUX AT TOP,0-160,1200-2200 CM-1=',F14.6,' W/M**2'/&
               ' FLUX AT TOP,0-2200 CM-1         =',F14.6,' W/M**2'/&
               ' NET FLUX AT GROUND,0-2200 CM-1  =',F14.6,' W/M**2'/  &
               ' NET FLUX DIFFERENCE,0-2200 CM-1 =',F14.6,' W/M**2')
9271  format ( /,  &
             ' FLUX AT TOP, BAND ',I2,' 1200-1400 CM-1 RANGE =',F14.6, &
                 ' W/M**2')
9272  format ( /,   &
               ' FLUX AT TOP, 1200-1400 CM-1 BAND  =',F14.6, &
                 ' W/M**2')
9280  format (40X,'CTS FLUXES'/   &
            1X,'BAND NO',8X,'LOFREQ',9X,'HIFREQ',9X,'F(1)',11X,'ACCUM. &
    &   F(1)',4X,'CTS F(GRD)',5X,'ACCUM. CTS F(GRD)')
9290  format (I11,6F15.6)
9300  format (11X,2F15.6)
9310  format (40x, 'EMISSIVITY FLUXES')
9320  format (/,2x,' lvl ',2x,   &
       ' h2o emiss ',' 560-800   ',' 800-900   ',  &
       ' 900-990   ',' 990-1070  ',' 1070-1200 ','  total    ') 
9321  format (/,2x,' lvl ',2x,   &
       ' h2o emiss ',' 560-800   ',' 800-900   ',' 900-990   ',   &
       ' 990-1070  ',' 1070-1200 ',' 1200-1400 ',  &
       '  total    ')
9322  format (/,2x,' lvl ',2x,  &
       ' h2o emiss ',' 560-800   ',' 800-900   ',' 900-990   ', &
       ' 990-1070  ',' 1070-1200 ',' 1200-1300 ', &
       ' 1300-1400 ',' 1200-1400 ',' total     ')
9323  format (/,2x,' lvl ',2x,   &
       ' h2o emiss ',' 560-800   ',' 800-900   ',' 900-990   ', &
       ' 990-1070  ',' 1070-1200 ',' 1200-1250 ',' 1250-1300 ', &
       ' 1300-1350 ',' 1350-1400 ',' 1200-1400 ',  &
       '  total    ')
9330  format (i5,-3p7f11.5)
9331  format (i5,-3p8f11.5)
9332  format (i5,-3p10f11.5)
9333  format (i5,-3p12f11.5)
9400  format (/,20X,' CLEAR-SKY LW HEATING RATES AND FLUXES',/)
9410  format (/,20X,' CLEAR-SKY SW HEATING RATES AND FLUXES',/)
9420  format (/,20X,' COMBINED CLEAR-SKY HEATING RATES AND FLUXES',/)
9510  format (27X,'CLOUD MICROPHYSICAL PARAMETERS ',/, 2x,'lyr',8x, &
      'liq water path', 3x, 'ice water path',3x, 'eff diam water', &
      ' eff diam ice')
9520  format (I5,7X,F14.6,3x,F14.6,3x,F14.6,3x,F14.6)


end subroutine radiag


!##################################################################



	      end module radiation_diag_mod

