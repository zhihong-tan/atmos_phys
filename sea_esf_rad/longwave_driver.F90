 
                           module longwave_driver_mod


use rad_output_file_mod,   only: hold_lw
use rad_step_setup_mod,    only: temp, rh2o, press, jabs, iabs, &
			         ISRAD, IERAD, JSRAD, JERAD, &
			         KS => KSRAD, KE => KERAD, lat_sv, &
			         Rad_time_sv
use  utilities_mod,        only: open_file, file_exist,    &
                                 check_nml_error, error_mesg, &
			         print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file
use constants_new_mod,     only: radcon
use longwave_setup_mod,    only: longwave_setup_dealloc, pdflux,  &
			         pdfinv, longwave_parameter_type, &
			         Lw_parameters
use longwave_clouds_mod,   only: longwave_clouds_dealloc,          &
                                 cldtau, cloud, thickcld 
use radiation_diag_mod,    only: radiation_diag_dealloc,   & 
			         radiag_from_driver
use cfc_mod,               only: cfc_indx8, cfc_indx8_part        
use longwave_fluxes_mod,   only: longwave_fluxes_ks, & 
			         longwave_fluxes_k_down,  &
                                 longwave_fluxes_KE_KEp1,  &
                                 longwave_fluxes_diag,   &
			         longwave_fluxes_dealloc,&
	               	         longwave_fluxes_alloc,  &
			         longwave_fluxes_sum
use cool_to_space_mod,     only: cool_to_space_approx,    &  
                                 cool_to_space_exact,  &
			         cool_to_space_get_hra, &
			         cool_to_space_dealloc
use longwave_tables_mod,   only: e1e290, e290, enear
use optical_path_mod,      only: optical_dealloc,  &
			         optical_trans_funct_from_KS, &
                                 optical_trans_funct_k_down, &
			         optical_trans_funct_KE, &
                                 optical_trans_funct_diag
use co2_source_mod,        only: co2_source_calc,     &
			         co2_sorc_dealloc
use gas_tf_mod,            only: co2coef, transcolrow, transfn, &
			         gas_dealloc, trans_nearby
use longwave_aerosol_mod,  only: get_totaerooptdep
use rad_utilities_mod,     only: longwave_control_type, Lw_control, &
			         radiation_control_type, Rad_control
use ozone_mod,             only: get_ozone
use cloudrad_package_mod,  only: get_clouds_for_lwrad

!------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!                        longwave radiation driver
!
!
!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.09'
    character(len=128)  :: version =  '$Id: longwave_driver.F90,v 1.2 2001/07/05 17:31:55 fms Exp $'
    character(len=128)  :: tag     =  '$Name: eugene $'

!---------------------------------------------------------------------
!-------  interfaces --------

public       &
	 longwave_driver_init, longwave_driver_alloc,  & 
	 longwave_driver_dealloc, lwrad, &
	 get_lw_output, get_lwflux_toa



!---------------------------------------------------------------------
!-------- namelist  ---------

logical     :: do_thick = .false.

namelist / longwave_driver_nml /    &
                                  do_thick


!---------------------------------------------------------------------
!------- public data ------




!---------------------------------------------------------------------
!------- private data ------

integer, parameter                  :: max_cld_bands = 7
integer, dimension (max_cld_bands)  :: cld_band_no, cld_indx
integer                             :: NLWCLDB, NBTRGE

!---------------------------------------------------------------------
!     flxnet  =  net longwave flux at model flux levels (including the
!                ground and the top of the atmosphere).
!     heatra  =  longwave heating rates in model layers.
!     flxnetcf, heatracf = values for corresponding variables  
!                          (without ...cf) computed for cloud-free case
!---------------------------------------------------------------------

real, dimension(:,:,:), allocatable           :: flxnet, flxnetcf, &
                                                 heatra, heatracf



!---------------------------------------------------------------------
!---------------------------------------------------------------------




contains





subroutine longwave_driver_init
 
!---------------------------------------------------------------------
    integer unit, ierr, io, nn

!---------------------------------------------------------------------
!----  read namelist ----

    if (file_exist('input.nml')) then
      unit = open_file ('input.nml', action='read')
      ierr=1; do while (ierr /= 0)
      read (unit, nml=longwave_driver_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'longwave_driver_nml')
      enddo
10    call close_file (unit)
    endif

    unit = open_file ('logfile.out', action='append')
!   call print_version_number(unit, 'longwave_driver', version_number)
    if (get_my_pe() == 0) then
      write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
      write (unit,nml=longwave_driver_nml)
    endif
    call close_file (unit)

!--------------------------------------------------------------------
!    define the cloud band index to be used in the longwave transmission
!    calculations for each cloud band. 
!--------------------------------------------------------------------
    NLWCLDB = Lw_parameters%NLWCLDB
    do nn=1,max_cld_bands
      cld_band_no(nn) = nn
      if (nn > NLWCLDB) then
        cld_indx(nn) = NLWCLDB
      else
        cld_indx(nn) = nn
      endif
    end do

    NBTRGE = Lw_parameters%NBTRGE

!--------------------------------------------------------------------
   

end  subroutine longwave_driver_init




!#####################################################################

subroutine longwave_driver_alloc
 
!-------------------------------------------------------------------
    
   allocate (flxnet ( ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
   allocate (heatra ( ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )

!--------------------------------------------------------------
!  initialization is needed so that avg history files are valid
!--------------------------------------------------------------
   flxnet(:,:,:) = 0.0
   heatra(:,:,:) = 0.0

   if (Rad_control%do_totcld_forcing) then 		   
     allocate (flxnetcf ( ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
     allocate (heatracf ( ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )

!--------------------------------------------------------------
!  initialization is needed so that avg history files are valid
!--------------------------------------------------------------
     flxnetcf(:,:,:) = 0.0
     heatracf(:,:,:) = 0.0
   endif


end  subroutine longwave_driver_alloc



!#####################################################################

subroutine longwave_driver_dealloc

   if (Rad_control%do_totcld_forcing) then 		   
     deallocate ( heatracf)
     deallocate ( flxnetcf)
   endif

   deallocate (heatra)
   deallocate (flxnet)

end subroutine longwave_driver_dealloc


!#####################################################################

subroutine lwrad

!--------------------------------------------------------------------
!
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1991), 9075-9096.
!
!     (2) schwarzkopf, m. d., and s. b. fels, "improvements to the
!         algorithm for computing co2 transmissivities and cooling
!         rates," journal geophysical research, 90 (1985) 10541-10550.
!
!     (3) fels, s.b., "simple strategies for inclusion of voigt
!         effects in infrared cooling calculations," application
!         optics, 18 (1979), 2634-2637.
!
!     (4) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: m. d. schwarzkopf
!
!     revised: 10/7/93
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------
!     intent in:
!
!     cmxolw  =  amounts of maximally overlapped longwave clouds in
!                layers from KS to KE.
!
!     crndlw  =  amounts of randomly overlapped longwave clouds in
!                layers from KS to KE.
!
!     emmxolw =  longwave cloud emissivity for maximally overlapped
!                clouds through layers from KS to KE. (default is one).
!
!     emrndlw =  longwave cloud emissivity for randomly overlapped
!                clouds through layers from KS to KE. (default is one).
!
!     kmxolw =  maximum number of maximally overlapped longwave clouds. 
!
!     krndlw =  maximum number of randomly overlapped longwave clouds. 
!
!     nmxolw  =  number of maximally overlapped longwave clouds 
!                at each grid point.
!
!     nrndlw  =  number of maximally overlapped longwave clouds 
!                at each grid point.
!
!     press  =  pressure at data levels of model.
!
!     qo3    =  mass mixing ratio of o3 at model data levels.
!
!     rh2o   =  mass mixing ratio of h2o at model data levels.
!
!     temp   =  temperature at data levels of model. 
!-----------------------------------------------------------------------
!     intent out:
!
!     flxnet  =  net longwave flux at model flux levels (including the 
!                ground and the top of the atmpsphere).
!
!     heatra  =  heating rate at data levels.
!
!     pflux   =  pressure at flux levels of model.
!
!     tflux   =  temperature assigned to model flux levels. 
!-----------------------------------------------------------------------
!     intent local:
!
!     atmden   =  atmospheric density, in gm/cm**2, for each of the
!                 KMAX layers.
!
!     co21c    =  transmission function for the 560-800 cm-1 band 
!                 from levels k through KE+1 to level k. used for flux
!                 at level k arising from exchange with levels k 
!                 through KE+1. includes co2 (from tables), h2o (after
!                 multiplication with over).
!                
!     co21diag =  transmission function for co2 only in the 560-800
!                 cm-1 band, from levels k to k, where k ranges from
!                 KS to KE+1. 
!
!     co21r    =  transmission function for the 560-800 cm-1 band 
!                 from level k to levels k+1 through KE+1. used for 
!                 flux at levels k+1 through KE+1 arising from
!                 exchange with level k. includes co2 (from tables),
!                 h2o (after multiplication with over).
!                
!     co2spnb  =  transmission functions from level KS to levels KS+1
!                 through KE+1, for the (NBCO215) narrow bands 
!                 comprising the 15um range (560-800 cm-1). used for
!                 exact CTS computations. includes co2 and possibly
!                 n2o (in the first such band).
!
!       to3cnt =  transmission function for the 990-1070 cm-1 band
!                 including o3(9.6 um) + h2o continuum (no lines) 
!                 and possibly cfcs.
!
!       sorc15 =  planck function for 560-800 cm-1 bands (sum over
!                 bands 9 and 10).
! 
!       sorc   =  planck function, at model temperatures, for all
!                 bands;  used in cool-to-space calculations.
!
!       tmp1   =  temporary array, used for computational purposes
!                 in various places. should have no consequences
!                 beyond the immediate module wherein it is defined.
!    tmp2,tmp3 =  temporary arrays, similar to tmp1
!
!         vv   =  layer-mean pressure in atmospheres. due to quadra-
!                 ture considerations, this does not equal the pressure
!                 at the data level (press).
!---------------------------------------------------------------------

   real, dimension(:,:,:,:), allocatable :: clddiag, cldtf
   real, dimension(:,:,:  ), allocatable :: cnttaub1, cnttaub2,  &
	 				    cnttaub3, co21c, co21diag, &
					    co21r, dsorc15, dsorc93, &
					    dsorcb1, dsorcb2, dsorcb3, &
					    dt4, emiss, heatem, overod,&
					    sorc15, t4, tmp1, tmp2, &
					    to3cnt, ch41c, n2o1c, n2o17c
   real, dimension(:,:,:),   allocatable :: emspec
   real, dimension(:,:),     allocatable :: flx1e1, s1a, flxge1, gxcts
   real, dimension(:,:,:,:), allocatable :: contdg
   real, dimension(:,:,:),   allocatable :: cfc_tf
   real, dimension (:,: ),   allocatable :: flx1e1cf, flxge1cf, &
					    gxctscf
   integer, dimension (:,: ),allocatable :: nmxolw, nrndlw
   real, dimension(:,:,:),   allocatable :: contodb1, contodb2,     &
				            contodb3, e1cts1, e1cts2,&
					    e1ctw1, e1ctw2, e1flx,   &
					    emisdg, emissb, flxcf, &
	    			            heatemcf, flx, flx1e1f, &
					    flxge1f, flx1e1fcf,   &
					    flxge1fcf, cts_sum, &
					    cts_sumcf, soe2, soe3,  &
					    soe4, soe5, to3dg, &
					    taero8, taero8kp,   &
					    totaer_tmp, tcfc8, &
					    tch4n2o_diag
   real, dimension(:,:,:),   allocatable :: cmxolw, crndlw
   real, dimension(:,:,:),   allocatable :: qo3
   real, dimension(:,:,:,:), allocatable :: emmxolw, emrndlw     
   real, dimension(:,:,:,:), allocatable :: e1cts1f, e1cts2f, e1flxf,&
					    emissf, emissbf,   &
					    emisdgf, emspecf, &
					    tch4n2oe

   integer   ::          n, k, kp, m, j, kmxolw, krndlw

!---------------------------------------------------------------------
!  allocate local arrays
!---------------------------------------------------------------------
   allocate ( clddiag   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NLWCLDB))
   allocate ( cldtf     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NLWCLDB))
   allocate ( cnttaub1  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( cnttaub2  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( cnttaub3  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( co21c     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( co21diag  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( co21r     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( dsorc15   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( dsorc93   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( dsorcb1   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( dsorcb2   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( dsorcb3   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( dt4       (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( emiss     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( heatem    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( overod    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( sorc15    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( t4        (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( tmp1      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( tmp2      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( to3cnt    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( ch41c     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( n2o1c     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( n2o17c    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1         ))
   allocate ( emspec    (ISRAD:IERAD, JSRAD:JERAD, 2               ))
   allocate ( flx1e1    (ISRAD:IERAD, JSRAD:JERAD                  ))
   allocate ( s1a       (ISRAD:IERAD, JSRAD:JERAD                  ))
   allocate ( flxge1    (ISRAD:IERAD, JSRAD:JERAD                  ))
   allocate ( gxcts     (ISRAD:IERAD, JSRAD:JERAD                  ))
   allocate ( contdg    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, 3      ))
   allocate ( cfc_tf    (ISRAD:IERAD, JSRAD:JERAD, KS:KE           ))

!---------------------------------------------------------------------
!     allocate space for and then retrieve ozne field from ozone_mod
!---------------------------------------------------------------------
   allocate (qo3 (ISRAD:IERAD, JSRAD:JERAD, KS:KE) )
   call get_ozone (qo3)

!-------------------------------------------------------------------
!     get the cloud properties needed by the longwave radiation code
!-------------------------------------------------------------------
   allocate ( cmxolw (ISRAD:IERAD, JSRAD:JERAD, KS:KE ) )
   allocate ( crndlw (ISRAD:IERAD, JSRAD:JERAD, KS:KE ) )
   allocate (emmxolw (ISRAD:IERAD, JSRAD:JERAD, KS:KE, NLWCLDB) )
   allocate (emrndlw (ISRAD:IERAD, JSRAD:JERAD, KS:KE, NLWCLDB) )
   allocate (nmxolw  (ISRAD:IERAD, JSRAD:JERAD) )
   allocate (nrndlw  (ISRAD:IERAD, JSRAD:JERAD) )
   call get_clouds_for_lwrad (cmxolw, crndlw, emmxolw, emrndlw, &
			      nmxolw, nrndlw)

!--------------------------------------------------------------------
!     find the maximum number of maximally and randomly overlapped
!     clouds for longwave radiation.
!-------------------------------------------------------------------
   kmxolw = MAXVAL(nmxolw)
   krndlw = MAXVAL(nrndlw) ! this value never used

!---------------------------------------------------------------------
!     allocate space in longwave_fluxes_mod
!-------------------------------------------------------------------
   call longwave_fluxes_alloc (6+NBTRGE)
       
!--------------------------------------------------------------------
!     allocate some arrays for cloud forcing and ch4n2o options
!-------------------------------------------------------------------
   if (Rad_control%do_totcld_forcing) then
     allocate ( flx1e1cf (ISRAD:IERAD, JSRAD:JERAD ) )
     allocate ( flxge1cf (ISRAD:IERAD, JSRAD:JERAD ) )
     allocate ( gxctscf  (ISRAD:IERAD, JSRAD:JERAD ) )
     if (Lw_control%do_ch4_n2o) then
       allocate ( flx1e1fcf (ISRAD:IERAD, JSRAD:JERAD, NBTRGE) )
       allocate ( flxge1fcf (ISRAD:IERAD, JSRAD:JERAD, NBTRGE) )
     endif
   endif

   if (Lw_control%do_ch4_n2o) then
     allocate( flx1e1f (ISRAD:IERAD, JSRAD:JERAD, NBTRGE) )
     allocate( flxge1f (ISRAD:IERAD, JSRAD:JERAD, NBTRGE) )
   end if

!--------------------------------------------------------------------
!     call co2coef to compute some co2 temperature and pressure   
!     interpolation quantities.
!-------------------------------------------------------------------
   call co2coef (press, temp)

!----------------------------------------------------------------------
!     call cldtau to compute cloud layer transmission functions for 
!     all layers.
!-------------------------------------------------------------------
   call cldtau (cmxolw, crndlw, emmxolw, emrndlw)

!---------------------------------------------------------------------
!     allocate some aerosol arrays and some ch4n2o arrays, if needed.
!---------------------------------------------------------------------
   if (Lw_control%do_lwaerosol) then
     allocate (totaer_tmp  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1 ) )
     allocate (taero8      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1 ) )
     allocate (taero8kp    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1 ) )
   endif
   if (Lw_control%do_ch4_n2o) then
     allocate (emissf      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
     allocate (emspecf     (ISRAD:IERAD, JSRAD:JERAD, 2, NBTRGE ) )
     allocate (tch4n2oe    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
     allocate (tch4n2o_diag(ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
     if (Lw_control%do_cfc) then
       allocate (tcfc8     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
     endif
   endif

!----------------------------------------------------------------------
!     call transfn to compute temperature-corrected co2 transmission 
!     functions (co2spnb and co2nbl). these fields remain in gas_tf_mod,
!     from where they are accessed as needed.
!---------------------------------------------------------------------
   call transfn 
  
!-----------------------------------------------------------------
!     compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!     mission functions and n2o 560-670 cm-1 transmission functions
!     appropriate for level KS. if active, save ch4n2o diagonal 
!     elements for use in nearby layer transmission function calcul-
!     ations. 
!---------------------------------------------------------------------
   if (Lw_control%do_ch4_n2o) then
     call transcolrow (KS, KS, KS, KE+1, KS+1, KE+1, co21c, co21r,  &
 	   	       tch4n2oe)
     tch4n2o_diag(:,:,KS:KE+1) = tch4n2oe(:,:,KS:KE+1,1)
   else
     call transcolrow (KS, KS, KS, KE+1, KS+1, KE+1, co21c, co21r)
   endif

!---------------------------------------------------------------------
!     save co2 diagonal elements for use in nearby layer co2
!     transmission function calculations. 
!---------------------------------------------------------------------
   co21diag    (:,:,KS)      = co21c   (:,:,KS)

!---------------------------------------------------------------------
!     go into optical_path_mod to obtain the optical path functions 
!     needed for use from level KS.
!     the 15um band transmission functions between levels KS and KS+1
!     are stored in overod and co2nbl; they will not be overwritten,
!     as long as calculations are made for pressure levels increasing
!     from KS.
!---------------------------------------------------------------------
   call optical_trans_funct_from_KS (to3cnt, overod, cnttaub1,    &
				     cnttaub2, cnttaub3)

!----------------------------------------------------------------------
!    obtain combined transmission functions for 560-800 cm-1 band,
!    pertaining to level KS.
!---------------------------------------------------------------------
   co21r (:,:,KS+1:KE+1) = overod(:,:,KS:KE)*co21r(:,:,KS+1:KE+1)
   co21c (:,:,KS+1:KE+1) = overod(:,:,KS:KE)*co21c(:,:,KS+1:KE+1)
                          
!-----------------------------------------------------------------------
!    compute cloud transmission functions between level KS and all
!    other levels.
!-----------------------------------------------------------------------
   call cloud (KS, cmxolw, crndlw, cldtf)

!----------------------------------------------------------------------
!    save "nearby layer" cloud transmission function for later use.
!-----------------------------------------------------------------------
   do n = 1,NLWCLDB
     clddiag(:,:,KS,n) = cldtf(:,:,KS,n)
   enddo

!----------------------------------------------------------------------
!    allocate some arrays needed for the source function. call 
!    co2_source_calc (co2_source_mod) to calculate the source function.
!    only those components of sorc that are needed in lwrad are returned
!-----------------------------------------------------------------------
   allocate (soe2 (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
   allocate (soe3 (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
   allocate (soe4 (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
   allocate (soe5 (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )

   call co2_source_calc (press, sorc15, soe2, soe3, soe4, soe5)
 
!----------------------------------------------------------------------
!    define "source function" appropriate for emissivity calculations
!    (temp**4), source functions for selected ranges including more 
!    than 1 frequency band (sorc15 for 15 um co2 band)
!    and differences in source functions (deltab) over
!    pressure layers.
!  
!    note: the values of sorc, sorc15, sorcwin, and derivatives 
!    depend on the no. of freq. bands!
!
!-----------------------------------------------------------------------
   t4(:,:,KS  :KE+1) = temp (:,:,KS:KE+1)**4

   dsorc93 (:,:,KS+1:KE+1) = soe2  (:,:,KS+1:KE+1) - soe2  (:,:,KS:KE)
   dsorcb1 (:,:,KS+1:KE+1) = soe3  (:,:,KS+1:KE+1) - soe3  (:,:,KS:KE)
   dsorcb2 (:,:,KS+1:KE+1) = soe4  (:,:,KS+1:KE+1) - soe4  (:,:,KS:KE)
   dsorcb3 (:,:,KS+1:KE+1) = soe5  (:,:,KS+1:KE+1) - soe5  (:,:,KS:KE)
   dsorc15 (:,:,KS+1:KE+1) = sorc15(:,:,KS+1:KE+1) - sorc15(:,:,KS:KE)
   dt4     (:,:,KS+1:KE+1) = t4    (:,:,KS+1:KE+1) -  t4   (:,:,KS:KE)

!-----------------------------------------------------------------------
!     obtain exact cool-to-space for water co2, and approximate cool-to-
!     space for co2 and o3.
!     REDO THIS COMMENT TO REFLECT CHANGES************************
!     obtain cool-to-space flux at the top by integration of heating
!     rates and using cool-to-space flux at the bottom (current value 
!     of gxcts).  note that the pressure quantities and conversion
!     factors have not been included either in excts or in gxcts.
!     these cancel out, thus reducing computations.
!-----------------------------------------------------------------------
   call cool_to_space_exact (cldtf, press, temp, to3cnt, gxcts,   &
			     gxctscf)

!---------------------------------------------------------------------
!     free up some no longer needed space
!----------------------------------------------------------------------
   call co2_sorc_dealloc
 
!----------------------------------------------------------------------
!     allocate space for emissivity flux arrays. compute the emissivity
!     fluxes for k=KS.
!----------------------------------------------------------------------
   allocate (e1cts1  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
   allocate (e1cts2  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
   allocate (e1ctw1  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
   allocate (e1ctw2  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
   allocate (e1flx   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
   if (Lw_control%do_ch4_n2o) then
     allocate (e1cts1f (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE))
     allocate (e1cts2f (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE))
     allocate (e1flxf  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE))
     call e1e290 (e1cts1, e1cts2, e1ctw1, e1ctw2, e1flx, emiss,  &
                  e1cts1f, e1cts2f, e1flxf, emissf)
   else
     call e1e290 (e1cts1, e1cts2, e1ctw1, e1ctw2, e1flx, emiss)
   endif

!---------------------------------------------------------------------
!    add the effects of other radiative gases on these flux arrays.
!    the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!    thus, tch4e and tn2oe are obtained directly from the transmission
!    functions. 
!----------------------------------------------------------------------
   if (Lw_control%do_ch4_n2o) then
     e1flxf (:,:,KS+1:KE+1,1) = e1flxf(:,:,KS+1:KE+1,1)*  &
                                tch4n2oe(:,:,KS+1:KE+1,1)
     e1cts1f(:,:,KS  :KE+1,1) = e1cts1f(:,:,KS:KE+1,1)*  &
                                tch4n2oe(:,:,KS:KE+1,1)
     e1cts2f(:,:,KS  :KE,  1) = e1cts2f(:,:,KS:KE,1)*   &
                                tch4n2oe(:,:,KS+1:KE+1,1)
     emissf (:,:,KS  :KE,  1) = emissf(:,:,KS:KE,1)*   &
                                tch4n2oe(:,:,KS+1:KE+1,1)
 
!----------------------------------------------------------------------
!    add cfc transmissivities if species which absorb in this fre-
!    quency range ( presently index 8) are present.
!----------------------------------------------------------------------
     if (Lw_control%do_cfc) then
       call cfc_indx8 (8, tcfc8)
       e1flxf (:,:,KS+1:KE+1,1) = e1flxf(:,:,KS+1:KE+1,1)*  &
                                  tcfc8(:,:,KS+1:KE+1)
       e1cts1f(:,:,KS  :KE+1,1) = e1cts1f(:,:,KS:KE+1,1)*  &
                                  tcfc8(:,:,KS:KE+1)
       e1cts2f(:,:,KS  :KE,  1) = e1cts2f(:,:,KS:KE,1)*  &
                                  tcfc8(:,:,KS+1:KE+1)
       emissf (:,:,KS  :KE,  1) = emissf(:,:,KS:KE,1)*   &
                                  tcfc8(:,:,KS+1:KE+1)
     endif 

!----------------------------------------------------------------------
!    compute aerosol transmission function for 1200-1400 cm-1 region
!----------------------------------------------------------------------
     if (Lw_control%do_lwaerosol) then
       call get_totaerooptdep(8, totaer_tmp)
       taero8(:,:,KS:KE+1) = EXP(-1.0E+00*totaer_tmp(:,:,KS:KE+1))
       e1flxf (:,:,KS+1:KE+1,1) = e1flxf(:,:,KS+1:KE+1,1)*   &
                                  taero8(:,:,KS+1:KE+1)
       e1cts1f(:,:,KS  :KE+1,1) = e1cts1f(:,:,KS:KE+1,1)*  &
                                  taero8(:,:,KS:KE+1)
       e1cts2f(:,:,KS  :KE,  1) = e1cts2f(:,:,KS:KE,1)*   &
                                  taero8(:,:,KS+1:KE+1)
       emissf (:,:,KS  :KE,  1) = emissf(:,:,KS:KE,1)*   &
                                  taero8(:,:,KS+1:KE+1)
     endif
   endif

!----------------------------------------------------------------------
!     the following is a rewrite of the original code largely to
!     eliminate three-dimensional arrays.  the code works on the
!     following principles.  let k be a fixed flux level and kp be
!     a varying flux level, then
! 
!     flux(k) = sum(deltab(kp)*tau(kp,k)) for kp=KS,KE+1.
!
!     if we assume a symmetrical array tau(k,kp)=tau(kp,k), we can
!     break down the calculations for k=KS,KE+1 as follows:
! 
!     flux(k) = sum(deltab(kp)*tau(kp,k)) for kp=k+1,KE+1            (1)
!
!     flux(kp) =   (deltab(k )*tau(kp,k)) for kp=k+1,KE+1.           (2)
!
!     plus deltab(k)*tau(k,k) for all k.
!
!     if we compute a one-dimensional array tauod(kp) for 
!     kp=k+1,KE+1, equations (1) and (2) become:
!
!     tauod(kp) = tau(kp,k)                                          (3)
!
!     flux (k ) = sum(deltab(kp)*tauod(kp)) for kp=k+1,KE+1          (4)
!
!     flux (kp) =    (deltab(k )*tauod(kp)) for kp=k+1,KE+1          (5)
!
!     where tau(k,k) and nearby layer terms are handled separately.
!
!     compute fluxes at level k = KS
!     compute the terms for flux at levels KS+1 to KE+1 from level KS.
!     compute terms for flux at level KS from level KS.
!     compute the terms for flux at level KS due to levels KP from KS+1 
!     to KE+1.
!
!-----------------------------------------------------------------------
    call longwave_fluxes_ks (t4, e1flx, 1, dt4, emiss, 0, &
			     cldtf, cld_indx(1) , cld_band_no(1))  

    call longwave_fluxes_ks (sorc15, co21r, 1, dsorc15, co21c, 1,   &
			     cldtf, cld_indx(2), cld_band_no(2))

    call longwave_fluxes_ks (soe3, cnttaub1, 0, dsorcb1, cnttaub1, 0,  &
			     cldtf, cld_indx(3) , cld_band_no(3))

    call longwave_fluxes_ks (soe4, cnttaub2, 0, dsorcb2, cnttaub2, 0,  &
			     cldtf, cld_indx(4), cld_band_no(4))

    call longwave_fluxes_ks (soe2, to3cnt, 0, dsorc93, to3cnt, 0,  &
			     cldtf, cld_indx(5), cld_band_no(5))

    call longwave_fluxes_ks (soe5, cnttaub3, 0, dsorcb3, cnttaub3, 0,  &
			     cldtf, cld_indx(6), cld_band_no(6))

    if (Lw_control%do_ch4_n2o) then
      do m=1,NBTRGE
        call longwave_fluxes_ks (t4, e1flxf(:,:,:,m), 1, dt4,    &
				 emissf(:,:,:,m), 0,cldtf,   &
				 cld_indx(7), cld_band_no(7) + m - 1)
      end do
    endif
 
!-----------------------------------------------------------------------
!     compute approximate cool-to-space heating rates for 1 wide band
!     in the 15um  range (560-800 cm-1) (ctsco2) and for 1 band in 
!     the 9.6 um band (ctso3).
!----------------------------------------------------------------------
    call cool_to_space_approx (1, sorc15, co21r, 0,    &
                               cldtf(:,:,:,cld_indx(2)))
    call cool_to_space_approx (2, soe2(:,:,:), to3cnt, 1,    &
	            	       cldtf(:,:,:,cld_indx(5)))

!-----------------------------------------------------------------------
!     compute the emissivity cool-to-space heating rates for the 
!     frequency ranges: 160-560, 800-990, and 1070-1200 cm-1. (the
!     latter 2 are combined in calculations). 
!----------------------------------------------------------------------
    call cool_to_space_approx (3, t4, e1ctw2, 1,    &
			       cldtf(:,:,:,cld_indx(1)),  &
			       e1ctw1, 0)
    call cool_to_space_approx (4, soe3, cnttaub1, 1,    &
			       cldtf(:,:,:,cld_indx(3)))
    call cool_to_space_approx (5, soe4, cnttaub2, 1,   &
				 cldtf(:,:,:,cld_indx(4)))
    call cool_to_space_approx (6, soe5, cnttaub3, 1,    &
			       cldtf(:,:,:,cld_indx(6)))
      
!-----------------------------------------------------------------------
!     obtain the flux at the top of the atmosphere in the 0-160, 
!     1200-2200 cm-1 frequency ranges, where heating rates and fluxes
!     are derived from h2o emissivity calculations (flx1e1) by:
!     1) obtaining the surface flux (flxge1); 2) summing the
!     emissivity flux divergence for these ranges (tmp1) over all 
!     pressure layers.
!#ifdef ch4n2o
!     if the 1200-1400 cm-1 range is computed separately, flux calcu-
!     lations are done separately in this range, then combined with
!     those from the other frequency range.
!#endif ch4n2o
!----------------------------------------------------------------------
    s1a   (:,:) = t4(:,:,KE+1)*(e1cts1(:,:,KE+1) - e1ctw1(:,:,KE+1))
    flxge1(:,:) = s1a(:,:)*cldtf(:,:,KE+1,1)
    tmp1(:,:,KS:KE) = t4(:,:,KS:KE)*    &
                      (e1cts1(:,:,KS:KE) - e1ctw1(:,:,KS:KE)) 
    tmp2(:,:,KS:KE) = t4(:,:,KS:KE)*   &
                      (e1cts2(:,:,KS:KE) - e1ctw2(:,:,KS:KE))
    flx1e1(:,:) = flxge1(:,:)
    do k=KS,KE
      flx1e1(:,:) = flx1e1(:,:) + tmp1(:,:,k)*cldtf(:,:,k,1) -   &
                    tmp2(:,:,k)*cldtf(:,:,k+1,1)
    enddo
    if (Rad_control%do_totcld_forcing) then
      flxge1cf(:,:) = s1a(:,:)
      flx1e1cf(:,:) = flxge1cf(:,:)
      do k=KS,KE
        flx1e1cf(:,:) = flx1e1cf(:,:) + tmp1(:,:,k) - tmp2(:,:,k)
      enddo
    endif
    if (Lw_control%do_ch4_n2o) then
      do m=1,NBTRGE
        s1a(:,:) = t4(:,:,KE+1)*e1cts1f(:,:,KE+1,m)
        flxge1f(:,:,m) = s1a(:,:)*cldtf(:,:,KE+1,cld_indx(7))
        flx1e1f(:,:,m) = flxge1f(:,:,m)
        do k=KS,KE
          tmp1(:,:,k) = t4(:,:,k)*e1cts1f(:,:,k,m)
          tmp2(:,:,k) = t4(:,:,k)*e1cts2f(:,:,k,m)
          flx1e1f(:,:,m) = flx1e1f(:,:,m) + tmp1(:,:,k)*   &
                           cldtf(:,:,k,cld_indx(7)) - tmp2(:,:,k)*  &
                           cldtf(:,:,k+1,cld_indx(7))
        end do
      end do
      do m=1,NBTRGE
        flx1e1(:,:) = flx1e1(:,:) + flx1e1f(:,:,m)
      enddo
      if (Rad_control%do_totcld_forcing) then
        do m=1,NBTRGE
	  flxge1fcf(:,:,m) = s1a(:,:)
	  flx1e1fcf(:,:,m) = s1a(:,:)
	  do k=KS,KE
            flx1e1fcf(:,:,m) = flx1e1fcf(:,:,m) + tmp1(:,:,k) -   &
		 			 	  tmp2(:,:,k)
	  end do
        end do
        do m=1,NBTRGE
	  flx1e1cf(:,:) = flx1e1cf(:,:) + flx1e1fcf(:,:,m)
        enddo
      endif
    endif
  
!-----------------------------------------------------------------------
!     deallocate some no longer needed emissitivy and source arrays 
!     before allocating some new ones 
!----------------------------------------------------------------------
    if (Lw_control%do_ch4_n2o) then
      deallocate (e1flxf )
      deallocate (e1cts2f)
      deallocate (e1cts1f)
    endif

    deallocate (e1flx)
    deallocate (e1ctw2)
    deallocate (e1ctw1)
    deallocate (e1cts2)
    deallocate (e1cts1)
    deallocate (soe5)
    deallocate (soe4)
    deallocate (soe3)
    deallocate (soe2)

    allocate (contodb1(ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
    allocate (contodb2(ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
    allocate (contodb3(ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
    allocate (emissb  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
    if (Lw_control%do_ch4_n2o) then
      allocate (emissbf (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE))
    endif

!-----------------------------------------------------------------------
!     perform flux calculations for the flux levels KS+1 to KE-1. calcu-
!     lations for flux levels KE and KE+1 are done separately, as all
!     calculations are special cases or nearby layers.
!----------------------------------------------------------------------
    do k=KS+1,KE-1

!--------------------------------------------------------------------
!     compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!     mission functions and n2o 560-670 cm-1 transmission functions
!     appropriate for level k. if activated, save ch4n2o tf term for 
!     use in nearby layer transmission function calculations.
!--------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        call transcolrow (k, k, k, KE+1, k+1, KE+1, co21c, co21r,    &
			  tch4n2oe)
        tch4n2o_diag(:,:,k+1:KE+1) = tch4n2oe(:,:,k+1:KE+1,1)
      else
        call transcolrow (k, k, k, KE+1, k+1, KE+1, co21c, co21r)
      endif
 
!---------------------------------------------------------------------
!     save co2 diagonal element for use in nearby layer co2 trans-
!     mission function calculations.
!---------------------------------------------------------------------
      co21diag(:,:,k) = co21c(:,:,k)

!-------------------------------------------------------------------
!     the 15 um band transmission functions between levels k and k+1
!     are stored in overod and co2nbl; they will not be overwritten,
!     as long as calculations are made for pressure levels increasing
!     from k.
!---------------------------------------------------------------------
      call optical_trans_funct_k_down (k, cnttaub1, cnttaub2, &
                                       cnttaub3, to3cnt, overod,  &
	 		               contodb1, contodb2, contodb3)

!---------------------------------------------------------------------
!     add the overod factor to then co21 transmission functions.
!---------------------------------------------------------------------
      co21c(:,:,k+1:KE+1) = overod(:,:,k:KE)*co21c(:,:,k+1:KE+1)
      co21r(:,:,k+1:KE+1) = overod(:,:,k:KE)*co21r(:,:,k+1:KE+1)

!-----------------------------------------------------------------------
!     compute cloud transmission functions between level k and all
!     other levels greater or equal to k.
!---------------------------------------------------------------------
      call cloud (k,cmxolw, crndlw, cldtf)

!-----------------------------------------------------------------------
!     save "nearby layer" cloud transmission function for later use.
!---------------------------------------------------------------------
      do n = 1,NLWCLDB
        clddiag(:,:,k,n) = cldtf(:,:,k,n)
      enddo

!-----------------------------------------------------------------------
!     compute the exchange terms in the flux equation (except the 
!     nearby layer (k,k) terms, done later).
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        call e290 (k, emiss, emissb, emissf, emissbf)
      else
        call e290 (k, emiss, emissb)
      endif

!----------------------------------------------------------------------
!    add the effects of other radiative gases on these flux arrays.
!    the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!    thus, tch4e and tn2oe are obtained directly from the transmission
!    functions. 
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        emissbf(:,:,k:KE,1) = emissbf(:,:,k:KE,1)*  &
	  		      tch4n2oe(:,:,k+1:KE+1,1)
        emissf(:,:,k:KE,1) = emissf(:,:,k:KE,1)*tch4n2oe(:,:,k+1:KE+1,1)

!--------------------------------------------------------------------
!    add cfc transmissivities if species which absorb in this fre-
!    quency range are present.
!--------------------------------------------------------------------
        if (Lw_control%do_cfc) then
          call cfc_indx8_part (8, tcfc8, k)
          emissbf(:,:,k:KE,1) = emissbf(:,:,k:KE,1)*tcfc8(:,:,k+1:KE+1)
          emissf(:,:,k:KE,1) = emissf(:,:,k:KE,1)*tcfc8(:,:,k+1:KE+1)
        endif

!--------------------------------------------------------------------
!     compute aerosol transmission function for 1200-1400 cm-1 region
!    (as quotient of 2 exponentials)
!     taero8kp(k) contains the (k+1,k) transmissivities for all k
!     in the 1200-1400 cm-1 frequency range.
!---------------------------------------------------------------------
        if (Lw_control%do_lwaerosol) then
          do kp = k+1,KE+1
            taero8kp(:,:,kp) = taero8(:,:,kp)/taero8(:,:,k)
          enddo
          emissbf(:,:,k:KE,1) = emissbf(:,:,k:KE,1)*  &
		                taero8kp(:,:,k+1:KE+1)
          emissf(:,:,k:KE,1) = emissf(:,:,k:KE,1)*  &
	     	               taero8kp(:,:,k+1:KE+1)
        endif
      endif

!-----------------------------------------------------------------------
!     compute the terms for flux at levels k+1 to KE+1 from level k.
!     compute the terms for flux at level k due to levels 
!     kp from k+1 to KE+1.
!----------------------------------------------------------------------
      call longwave_fluxes_k_down  (k, dt4, emissb, emiss, 0,  &
                                    cldtf(:,:,:,cld_indx(1)),  &
                                    cld_band_no(1)          )  
      call longwave_fluxes_k_down  (k, dsorc15, co21r, co21c, 1,    &
                                    cldtf(:,:,:,cld_indx(2)),  &
                                    cld_band_no(2)          )  
      call longwave_fluxes_k_down  (k, dsorcb1, contodb1, contodb1, 0, &
                                    cldtf(:,:,:,cld_indx(3)),  &
                                    cld_band_no(3)          )  
      call longwave_fluxes_k_down  (k, dsorcb2, contodb2, contodb2, 0, &
                                    cldtf(:,:,:,cld_indx(4)),  &
                                    cld_band_no(4)          )  
      call longwave_fluxes_k_down  (k, dsorc93, to3cnt, to3cnt, 0,  &
                                    cldtf(:,:,:,cld_indx(5)),  &
                                    cld_band_no(5)          )  
      call longwave_fluxes_k_down  (k, dsorcb3, contodb3, contodb3, 0, &
                                    cldtf(:,:,:,cld_indx(6)),  &
                                    cld_band_no(6)          )  
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          call longwave_fluxes_k_down (k, dt4, emissbf(:,:,:,m),   &
				       emissf(:,:,:,m), 0,  &
                                       cldtf(:,:,:,cld_indx(7)),  &
                                       cld_band_no(7) + m - 1)  
        end do
      endif
   end do   ! (end of k=KS+1,KE-1 loop)

!-----------------------------------------------------------------------
!     compute remaining flux terms. these include:
!       1) the (k,k) terms, for pressure levels k from KS+1 to KE-1 
!          (the KS,KS term was handled earlier);
!       2) terms for pressure level KE. these include the (KE,KE) term,
!          computed as in (1), and the (KE,KE+1) and (KE+1,KE) terms,
!          computed somewhat differently from the similar terms at
!          higher levels, owing to the proximity to the surface layer
!          KE+1;
!       3) the term for pressure level KE+1 (the (KE+1,KE+1 term).
!
!     compute k=KE case.  since the kp loop is length one, many 
!     simplifications occur.  the co2 quantities and the emissivity
!     quantities) are computed in the nbl section. therefore, we want
!     to compute over, to3cnt, and contod; according to our notation    
!     over(:,:,KE), to3cnt(:,:,KE), and contod(:,:,KE).  the boundary
!     layer and nearby layer corrections to the transmission functions 
!     are obtained above.  the following ratios are used in various nbl
!     nbl calculations.  the remaining calculations are for:
!
!       1) the (k,k) terms, k=KS+1,KE-1;
!       2) the (KE,KE    ) term;
!       3) the (KE,KE+1  ) term;
!       4) the (KE+1,KE  ) term;
!       5) the (KE+1,KE+1) term.
!
!     each is uniquely handled.  different flux terms are computed
!     differently the fourth section obtains water transmission 
!     functions used in q(approximate) calculations and also makes nbl 
!     corrections:
!  
!       1) emiss (:,:) is the transmission function matrix obtained 
!          using E2spec;
! 
!       2) "nearby layer" corrections (emiss(i,i)) are obtained
!          using E3v88;
! 
!       3) special values at the surface (emiss(KE,KE+1),
!          emiss(KE+1,KE), emiss(KE+1,KE+1)) are calculated.
!
!
!     compute temperature and/or scaled amount indices and residuals 
!     for nearby layer and special layer lookup tables.
!
!          calculation for special cases (KE,KE+1) and (KE+1,KE)
!
!     compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!     mission functions and n2o 560-670 cm-1 transmission functions
!     appropriate for level KE. if activated, save ch4n2o tf term.
!--------------------------------------------------------------------
    if (Lw_control%do_ch4_n2o) then
      call transcolrow (KE, KE, KE, KE+1, KE+1, KE+1, co21c, co21r,   &
		        tch4n2oe)
      tch4n2o_diag(:,:,KE+1) = tch4n2oe(:,:,KE+1,1)
    else
      call transcolrow (KE, KE, KE, KE+1, KE+1, KE+1, co21c, co21r)
    endif

!----------------------------------------------------------------------
!     get optical path terms for KE
!----------------------------------------------------------------------
    call optical_trans_funct_KE (cnttaub1, cnttaub2, cnttaub3, to3cnt, &
				overod, contodb1, contodb2, contodb3)

!-----------------------------------------------------------------------
!     compute cloud transmission functions between level KE and KE and
!     KE+1
!----------------------------------------------------------------------
    call cloud (KE, cmxolw, crndlw, cldtf)

!-----------------------------------------------------------------------
!     save "nearby layer" cloud transmission function for later use.
!----------------------------------------------------------------------
    do n = 1,NLWCLDB
      clddiag(:,:,KE,n) = cldtf(:,:,KE,n)
    enddo
   
!---------------------------------------------------------------------
!     allocate space for emissivity arrays for level KE
!----------------------------------------------------------------------
    allocate (emisdg  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
    if (Lw_control%do_ch4_n2o) then
      allocate (emisdgf (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1,NBTRGE) )
    endif

!---------------------------------------------------------------------
!     call enear to calculate emissivity arrays
!----------------------------------------------------------------------
    if (Lw_control%do_ch4_n2o) then
      call enear (emisdg, emspec, emisdgf, emspecf)
    else
      call enear (emisdg, emspec)
    endif

!--------------------------------------------------------------------
!    add the effects of other radiative gases on these flux arrays.
!    the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!    thus, tch4e and tn2oe are obtained directly from the transmission
!    functions. 
!----------------------------------------------------------------------
    if (Lw_control%do_ch4_n2o) then
      do k=KS,KS+1
	emspecf(:,:,K,1) = emspecf(:,:,K,1)*tch4n2oe(:,:,KE+1,1)
      end do
      emisdgf(:,:,KS+1:KE+1,1) = emisdgf(:,:,KS+1:KE+1,1) *   &
                                 tch4n2o_diag(:,:,KS+1:KE+1)

!--------------------------------------------------------------------
!     add cfc transmissivities if species which absorb in this fre-
!    quency range are present.
!----------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_indx8_part (8, tcfc8, KE)
	do k=KS,KS+1
	  emspecf(:,:,K,1) = emspecf(:,:,K,1)*tcfc8(:,:,KE+1)
	end do
	emisdgf(:,:,KS+1:KE+1,1) = emisdgf(:,:,KS+1:KE+1,1) *   &
                                   tcfc8(:,:,KS+1:KE+1)
      endif
    endif

!----------------------------------------------------------------------
!     compute nearby layer transmission functions for 15 um band, cont-
!     inuum bands, and 9.3 um band in subroutine Nearbylyrtf. trans-
!     mission functions for the special cases (KE,KE+1) and (KE+1,KE)
!     are also computed for the 15 um band.
!----------------------------------------------------------------------
    call trans_nearby (press, overod, co21c, co21diag, co21r)

!-----------------------------------------------------------------------
!     obtain fluxes for the two terms (KE,KE+1) and (KE+1,KE), both 
!     using the same cloud transmission functions (from layer KE)
!----------------------------------------------------------------------
    call longwave_fluxes_KE_KEp1 (dt4, emspec(:,:,KS+1),    &
				  emspec(:,:,KS),   &
				  cldtf(:,:,:,cld_indx(1)),  &
				  cld_band_no(1) )
    call longwave_fluxes_KE_KEp1 (dsorc15, co21c(:,:,KE+1), &
				  co21r(:,:,KE+1),  &
                                  cldtf(:,:,:,cld_indx(2)), &
				  cld_band_no(2) )
    call longwave_fluxes_KE_KEp1 (dsorcb1, contodb1(:,:,KE),  &
				  contodb1(:,:,KE),&
                                  cldtf(:,:,:,cld_indx(3)),  &
				  cld_band_no(3) )
    call longwave_fluxes_KE_KEp1 (dsorcb2, contodb2(:,:,KE),  &
				  contodb2(:,:,KE),&
                                  cldtf(:,:,:,cld_indx(4)),  &
				  cld_band_no(4) )
    call longwave_fluxes_KE_KEp1 (dsorc93, to3cnt(:,:,KE),   &
				  to3cnt(:,:,KE), &
                                  cldtf(:,:,:,cld_indx(5)),  &
				  cld_band_no(5) )
    call longwave_fluxes_KE_KEp1 (dsorcb3, contodb3(:,:,KE),   &
				  contodb3(:,:,KE),&
                                  cldtf(:,:,:,cld_indx(6)),  &
				  cld_band_no(6) )
    if (Lw_control%do_ch4_n2o) then
      do m=1,NBTRGE
        call longwave_fluxes_KE_KEp1 (dt4, emspecf(:,:,KS+1,m),   &
                                      emspecf(:,:,KS,m),   &
                                      cldtf(:,:,:,cld_indx(7)),    &
	  			      cld_band_no(7) + m - 1)
      end do
    endif

!------------------------------------------------------------------
!      allocate array needed for diagonal terms and for heating rate
!      if cloud forcing is activated
!------------------------------------------------------------------
    allocate (to3dg   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1))
    if (Rad_control%do_totcld_forcing) then
      allocate (heatemcf(ISRAD:IERAD, JSRAD:JERAD, KS:KE+1))
    endif

!-------------------------------------------------------------------
!     obtain optical path transmission functions for diagonal terms
!------------------------------------------------------------------
    call optical_trans_funct_diag (press, contdg, to3dg)
 
!-----------------------------------------------------------------------
!     compute cloud transmission functions between level KE+1 and KE+1
!------------------------------------------------------------------
    call cloud (KE+1, cmxolw, crndlw, cldtf)
 
!-----------------------------------------------------------------------
!     save "nearby layer" cloud transmission function for later use.
!------------------------------------------------------------------
    do n = 1,NLWCLDB
      clddiag(:,:,KE+1,n) = cldtf(:,:,KE+1,n)
    enddo

!-----------------------------------------------------------------------
!     obtain fluxes for the diagonal terms at all levels.
!-----------------------------------------------------------------------
    call longwave_fluxes_diag (dt4, emisdg,   &
			       clddiag(:,:,:,cld_indx(1)),   &
			       cld_band_no(1) )
    call longwave_fluxes_diag (dsorc15, co21diag,   &
                               clddiag(:,:,:,cld_indx(2)),  &
			       cld_band_no(2) )
    call longwave_fluxes_diag (dsorcb1, contdg(:,:,:,1), &
                               clddiag(:,:,:,cld_indx(3)),  &
			       cld_band_no(3) )
    call longwave_fluxes_diag (dsorcb2, contdg(:,:,:,2),  &
                               clddiag(:,:,:,cld_indx(4)),   &
			       cld_band_no(4) )
    call longwave_fluxes_diag (dsorc93, to3dg ,   &
                               clddiag(:,:,:,cld_indx(5)),   &
			       cld_band_no(5) )
    call longwave_fluxes_diag (dsorcb3, contdg(:,:,:,3),  &
                               clddiag(:,:,:,cld_indx(6)),  &
			       cld_band_no(6) )
    if (Lw_control%do_ch4_n2o) then
      do m=1, NBTRGE
        call longwave_fluxes_diag (dt4, emisdgf(:,:,:,m),    &
                                   clddiag(:,:,:,cld_indx(7)),  &
				   cld_band_no(7) + m - 1)
      end do
    endif

!--------------------------------------------------------------------
!      allocate arrays and sum up fluxes over bands
!-----------------------------------------------------------------------
    allocate (flx (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
    if (Rad_control%do_totcld_forcing) then
      allocate (flxcf (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
    endif 

    if (Rad_control%do_totcld_forcing) then
      call longwave_fluxes_sum (flx, NBTRGE, flxcf)
    else
      call longwave_fluxes_sum (flx, NBTRGE)
    endif

!-----------------------------------------------------------------------
!     compute emissivity heating rates.
!-----------------------------------------------------------------------
    heatem(:,:,KS:KE) = radcon*(flx(:,:,KS+1:KE+1) -    &
                        flx(:,:,KS:KE))*pdfinv(:,:,KS:KE)
    if (Rad_control%do_totcld_forcing) then 		   
      heatemcf(:,:,KS:KE) = radcon*(flxcf(:,:,KS+1:KE+1) -    &
                            flxcf(:,:,KS:KE))*pdfinv(:,:,KS:KE)
    endif

!--------------------------------------------------------------------
!     deallocate some arrays
!-----------------------------------------------------------------------
    if (Rad_control%do_totcld_forcing) then
      deallocate (flxcf)
    endif
    deallocate (flx)

!-----------------------------------------------------------------------
!     compute total heating rates.
!-----------------------------------------------------------------------
    allocate ( cts_sum (ISRAD:IERAD, JSRAD:JERAD, KS:KE) )
    if (Rad_control%do_totcld_forcing) then
      allocate ( cts_sumcf (ISRAD:IERAD, JSRAD:JERAD, KS:KE) )
    endif 

    if (Rad_control%do_totcld_forcing) then 		   
      call cool_to_space_get_hra (cts_sum, cts_sumcf)
    else
      call cool_to_space_get_hra (cts_sum)
    endif

    heatra(:,:,KS:KE) = heatem(:,:,KS:KE) + cts_sum  (:,:,KS:KE) 
    if (Rad_control%do_totcld_forcing) then 		   
      heatracf(:,:,KS:KE) = heatemcf(:,:,KS:KE) + cts_sumcf(:,:,KS:KE) 
      deallocate (cts_sumcf)
    endif
    deallocate (cts_sum)

!-----------------------------------------------------------------------
!     compute the flux at each flux level using the flux at the
!     top (flx1e1 + gxcts) and the integral of the heating rates 
!-----------------------------------------------------------------------
    flxnet(:,:,KS   ) = flx1e1(:,:) + gxcts(:,:)
    tmp1 (:,:,KS:KE) = heatra(:,:,KS:KE)*pdflux(:,:,KS:KE)/radcon 
    do k=KS+1,KE+1
      flxnet(:,:,k) = flxnet(:,:,k-1) + tmp1(:,:,k-1)
    enddo

   if (Rad_control%do_totcld_forcing) then 		   
     flxnetcf(:,:,KS   ) = flx1e1cf(:,:) + gxctscf(:,:)
     tmp1 (:,:,KS:KE) = heatracf(:,:,KS:KE)*pdflux(:,:,KS:KE)/radcon 
     do k=KS+1,KE+1
       flxnetcf(:,:,k) = flxnetcf(:,:,k-1) + tmp1(:,:,k-1)
     enddo
   endif

!-----------------------------------------------------------------------
!    call thickcld to perform "pseudo-convective adjustment" for
!    maximally overlapped clouds, if desired.
!-----------------------------------------------------------------------
   if (do_thick) then
     call thickcld (cmxolw, emmxolw, kmxolw, flxnet, heatra)
   endif

!---------------------------------------------------------------------
!    if diagnostics are desired for this chunk, send the data to the
!     radiation_diag_mod.
!-----------------------------------------------------------------------
   if (Rad_control%do_diagnostics) then
     do j=JSRAD,JERAD
       call radiag_from_driver (j, flx1e1, flx1e1f, heatra, heatracf,  &
				flxnet, flxnetcf)
     end do
   endif

!-------------------------------------------------------------------
!    deallocate the remaining longwave arrays
!-----------------------------------------------------------------------
   if (Rad_control%do_totcld_forcing) then
     deallocate (heatemcf)
   endif
   deallocate (to3dg)
   if (Lw_control%do_ch4_n2o) then
     deallocate (emisdgf)
   endif
   deallocate (emisdg)
   if (Lw_control%do_ch4_n2o) then
     deallocate (emissbf)
   endif
   deallocate (emissb)
   deallocate (contodb3)
   deallocate (contodb2)
   deallocate (contodb1)
   if (Lw_control%do_ch4_n2o) then
     if (Lw_control%do_cfc) then
       deallocate (tcfc8)
     endif
     deallocate (tch4n2oe)
     deallocate (tch4n2o_diag)
     deallocate (emspecf)
     deallocate (emissf )
   endif
   if (Lw_control%do_lwaerosol) then 
     deallocate (taero8kp)
     deallocate (taero8  )
     deallocate (totaer_tmp )
   endif
   if (Lw_control%do_ch4_n2o) then
     deallocate (flxge1f )
     deallocate (flx1e1f )
   end if
   if (Rad_control%do_totcld_forcing) then
     if (Lw_control%do_ch4_n2o) then
       deallocate (flxge1fcf)
       deallocate (flx1e1fcf)
     endif
     deallocate (gxctscf )
     deallocate (flxge1cf)
     deallocate (flx1e1cf)
   endif

!------------------------------------------------------------------
!    go into longwave_fluxes and deallocate the arrays allocated there
!-----------------------------------------------------------------------
   call longwave_fluxes_dealloc

   deallocate (cmxolw  )
   deallocate (crndlw  )
   deallocate (emmxolw )
   deallocate (emrndlw )
   deallocate (nmxolw  )
   deallocate (nrndlw  )
   deallocate (qo3)
   deallocate (cfc_tf      )
   deallocate (contdg      )
   deallocate (gxcts       )
   deallocate (flxge1      )
   deallocate (s1a         )
   deallocate (flx1e1      )
   deallocate (emspec      )
   deallocate (n2o17c      )
   deallocate (n2o1c       )
   deallocate (ch41c       )
   deallocate (to3cnt      )
   deallocate (tmp2        )
   deallocate (tmp1        )
   deallocate (t4          )
   deallocate (sorc15      )
   deallocate (overod      )
   deallocate (heatem      )
   deallocate (emiss       )
   deallocate (dt4         )
   deallocate (dsorcb3     )
   deallocate (dsorcb2     )
   deallocate (dsorcb1     )
   deallocate (dsorc93     )
   deallocate (dsorc15     )
   deallocate (co21r       )
   deallocate (co21diag    )
   deallocate (co21c       )
   deallocate (cnttaub3    )
   deallocate (cnttaub2    )
   deallocate (cnttaub1    )
   deallocate (cldtf       )
   deallocate (clddiag     )

!------------------------------------------------------------------
!    now deallocate the remaining arrays allocated in the longwave    
!    routines
!------------------------------------------------------------------
   call longwave_clouds_dealloc
   call gas_dealloc

!--------------------------------------------------------------------
!   pass the output fluxes and heating rates to the archive_package
!   for later use
!--------------------------------------------------------------------
   if (Rad_control%do_totcld_forcing) then
     call hold_lw (flxnet, heatra, flxnetcf, heatracf)
   else
     call hold_lw (flxnet, heatra)
   endif
!--------------------------------------------------------------------


end subroutine lwrad




!####################################################################

subroutine get_lw_output (heatra_out, flxnet_out, heatracf_out, &
			  flxnetcf_out)

real, dimension(:,:,:), intent(out)              :: heatra_out,   &
						    flxnet_out
real, dimension(:,:,:), intent(out), optional    :: heatracf_out,   &
						    flxnetcf_out

!----------------------------------------------------------------------
!  here the "..._out" arrays have model dimensions while the right-side
!  arrays are dimensioned by the radiation limits (ISRAD, etc.). if 
!  these differ, the proper mapping must be made here.
!----------------------------------------------------------------------

    heatra_out(:,:,:) = heatra(:,:,:)
    flxnet_out(:,:,:) = flxnet(:,:,:)

    if (Rad_control%do_totcld_forcing) then
      heatracf_out(:,:,:) = heatracf(:,:,:)
      flxnetcf_out(:,:,:) = flxnetcf(:,:,:)
    endif

end subroutine get_lw_output

!---------------------------------------------------------------------

!####################################################################

subroutine get_lwflux_toa (flxnet_out )

real, dimension(:,:), intent(out)              :: flxnet_out

!----------------------------------------------------------------------
!  here the "..._out" arrays have model dimensions while the right-side
!  arrays are dimensioned by the radiation limits (ISRAD, etc.). if 
!  these differ, the proper mapping must be made here.
!----------------------------------------------------------------------

   flxnet_out(:,:) = flxnet(:,:,KS)


end subroutine get_lwflux_toa

!---------------------------------------------------------------------


		  end module longwave_driver_mod
