 
                           module sealw99_mod


use  utilities_mod,        only: open_file, file_exist,    &
                                 check_nml_error, error_mesg, &
			         print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file, &
				read_data, write_data
use constants_mod,         only: diffac, radcon_mks, seconds_per_day, radcon
use longwave_clouds_mod,   only: longwave_clouds_init, &
                                 cldtau, cloud, thickcld 
!use cfc_mod,               only: cfc_indx8, cfc_indx8_part, &
!                                cfc_exact
use longwave_fluxes_mod,   only: longwave_fluxes_ks, & 
				 longwave_fluxes_init, &
			         longwave_fluxes_k_down,  &
                                 longwave_fluxes_KE_KEp1,  &
                                 longwave_fluxes_diag,   &
			         longwave_fluxes_sum
use longwave_tables_mod,   only:                      &
                                 longwave_tables_init
use optical_path_mod,      only:                   &
				 optical_path_setup, &
				 optical_path_init,  &
			         optical_trans_funct_from_KS, &
                                 optical_trans_funct_k_down, &
			         optical_trans_funct_KE, &
                                 optical_trans_funct_diag, &
                                              get_totvo2,   &
				get_totch2o, get_totch2obd
use gas_tf_mod,            only: co2coef, transcolrow,          &
				  get_control_gas_tf, &
                                  gas_tf_init, &
                                 transcol, &
                                 gas_tf_end,  &
					      trans_sfc, &
			                      trans_nearby
use rad_utilities_mod,     only: longwave_control_type, Lw_control, &
                                 cldrad_properties_type, &
				 lw_output_type, &
		                 longwave_tables1_type,  &
				 longwave_tables2_type,  &
				 longwave_tables3_type, &
				 atmos_input_type, &
				 radiative_gases_type, &
                                 aerosol_type, &
				 optical_path_type, &
                                 gas_tf_type, &
				 lw_table_type, &
			                 longwave_parameter_type, &
			         Lw_parameters, &
				 lw_diagnostics_type, &
				 lw_clouds_type, &
				locate_in_table, &
				 environment_type, Environment, &
				 looktab,  &
				mass_1, temp_1, &
				table_axis_type, &
			         radiation_control_type, Rad_control
use lw_gases_stdtf_mod, only : lw_gases_stdtf_init, &
                                 cfc_indx8, cfc_indx8_part, &
                                cfc_exact , &
                                lw_gases_stdtf_time_vary, &
                                lw_gases_stdtf_dealloc, &
                                 co2_lblinterp, &
                                ch4_lblinterp, n2o_lblinterp

use longwave_params_mod, only : longwave_params_init, &
                                NBCO215, &
                                NBLY_CKD2P1, NBLY_ORIG

!------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!                        sealw99 radiation code  
!
!
!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

    character(len=128)  :: version =  '$Id: sealw99.F90,v 1.3 2003/04/09 21:01:44 fms Exp $'
    character(len=128)  :: tag     =  '$Name: inchon $'

!---------------------------------------------------------------------
!-------  interfaces --------

public       &
	 sealw99_init,   sealw99_end,  & 
	                          sealw99

private   &
!         longwave_driver_alloc, &
	                      cool_to_space_approx,   &
	  cool_to_space_exact, &
          e1e290, e290, enear, esfc, &
                           co2_source_calc, &
          nlte, co2curt


!---------------------------------------------------------------------
!-------- namelist  ---------

logical     :: do_thick = .false.
logical     :: do_ckd2p1 = .true.
logical     :: do_ch4n2olbltmpint   = .false.
logical                 :: do_nlte = .false.
 

namelist / sealw99_nml /    &
                                  do_thick, &
                                     do_ckd2p1,   &
                                    do_nlte, &
                                 do_ch4n2olbltmpint


!---------------------------------------------------------------------
!------- public data ------




!---------------------------------------------------------------------
!------- private data ------

!integer, parameter                  :: max_cld_bands = 7
!integer, dimension (max_cld_bands)  :: cld_band_no, cld_indx
!integer, dimension (max_cld_bands)  ::              cld_indx
integer, dimension (:), allocatable  ::              cld_indx
integer                             :: NLWCLDB, NBTRGE
!integer, parameter :: NTTABH2O = 28

!---------------------------------------------------------------------
!     flxnet  =  net longwave flux at model flux levels (including the
!                ground and the top of the atmosphere).
!     heatra  =  longwave heating rates in model layers.
!     flxnetcf, heatracf = values for corresponding variables  
!                          (without ...cf) computed for cloud-free case
!---------------------------------------------------------------------


       integer ::                             ks, ke

type(longwave_tables3_type) :: tabsr
type (longwave_tables1_type)             :: tab1, tab2, tab3, tab1w
type (longwave_tables2_type)             :: tab1a, tab2a, tab3a

!---------------------------------------------------------------------
!     apcm, bpcm   =   capphi coefficients for NBLY bands.
!     atpcm, btpcm =   cappsi coefficients for NBLY bands.
!     acomb        =   random "a" parameter for NBLY bands.
!     bcomb        =   random "b" parameter for NBLY bands.
!---------------------------------------------------------------------
      real, dimension (:), allocatable    ::  apcm, bpcm, atpcm, btpcm,&
					      acomb, bcomb


!-------------------------------------------------------------------
      integer, parameter                  ::  no_combined_bands = 8
      real, dimension(no_combined_bands)  ::  band_no_start, band_no_end
      integer, dimension(NBLY_CKD2P1-1)   ::  cld_indx_table

!---------------------------------------------------------------------
!     NBLY    =  number of frequency bands for exact cool-to-
!                space computations.
!---------------------------------------------------------------------
      integer                             ::  NBLY
!---------------------------------------------------------------------
real, dimension (:),      allocatable       :: c1b7, c2b7


integer            :: ixprnlte
!logical :: do_diagnostics
!---------------------------------------------------------------------

logical :: do_ch4_n2o_init = .true.
logical    ::  do_co2_init = .true.



contains




!#####################################################################

subroutine sealw99_init (latb, lonb, pref, Lw_tables)
 
real, dimension(:), intent(in) :: latb, lonb
real, dimension(:,:), intent(in) :: pref
type(lw_table_type), intent(inout) :: Lw_tables

!---------------------------------------------------------------------
!    integer                 nn

   real, dimension(size(pref,1) ) :: plm
   real, dimension (NBCO215) :: cent, del

   integer         :: unit, ierr, io, k, n, m, nn
   integer         :: ioffset
   real            :: prnlte
   integer         ::     kmax, kmin

       real, dimension(NBLY_ORIG)   :: apcm_orig, bpcm_orig,     &
				       atpcm_orig, btpcm_orig,   &
				       acomb_orig, bcomb_orig
       real, dimension(NBLY_CKD2P1) :: apcm_ckd, bpcm_ckd, atpcm_ckd,  &
			               btpcm_ckd, acomb_ckd, bcomb_ckd 

!      integer                      ::   m

!----------------------------------------------------------------------
!    2) 160-560 (as 8 bands using combined bands). program gasbnd is
!    used as 40 bands (160-560,10 cm-1 bandwidth) with ifdef icomb
!    on. combined bands defined according to index iband.
!----------------------------------------------------------------------
 
       data (acomb_orig(k),k=1,8)    /     &
         0.151896E+05,  0.331743E+04,  0.526549E+03,  0.162813E+03,  &
         0.268531E+03,  0.534040E+02,  0.267777E+02,  0.123003E+02/
       data (bcomb_orig(k),k=1,8)    /    &
         0.153898E+00,  0.116352E+00,  0.102814E+00,  0.966280E-01,  &
         0.128586E+00,  0.119393E+00,  0.898891E-01,  0.645340E-01/
       data (apcm_orig(k),k=1,8)    /      &
        -0.531605E-03,  0.679588E-02,  0.145133E-01,  0.928230E-02,  &
         0.120727E-01,  0.159210E-01,  0.175268E-01,  0.150281E-01/
       data (bpcm_orig(k),k=1,8)    /     &
        -0.122002E-04, -0.335081E-04, -0.449999E-04, -0.231290E-04,  &  
        -0.385298E-04, -0.157170E-04,  0.183434E-04, -0.501718E-04/
       data (atpcm_orig(k),k=1,8)    /   &
        -0.156433E-02,  0.611348E-02,  0.134469E-01,  0.871421E-02,  &
         0.133191E-01,  0.164697E-01,  0.199640E-01,  0.163219E-01/
       data (btpcm_orig(k),k=1,8)    /    &
        -0.698856E-05, -0.295269E-04, -0.503674E-04, -0.268392E-04,  &
        -0.496663E-04, -0.333122E-04, -0.337346E-04, -0.258706E-04/

!----------------------------------------------------------------------
!    3) 560-1200 and 4.3 um band (8 bands, frequency range given
!    by bdlocm-bdhicm). program gasbnd is used with 8 specified
!    bandwidths.
!----------------------------------------------------------------------
 
       data (acomb_orig(k),k=9,NBLY_orig)    /    &
         0.790655E+01,  0.208292E+01,  0.435231E+00,  0.501124E-01,   &
         0.167058E-01,  0.213731E-01,  0.211154E+00,  0.458804E-02/
       data (bcomb_orig(k),k=9,NBLY_orig)    /     &
         0.676982E-01,  0.691090E-01,  0.660834E-01,  0.730815E-01,   &
         0.654709E-01,  0.848528E-01,  0.852302E-01,  0.242787E+00/
       data (apcm_orig(k),k=9,NBLY_orig)    /      &
         0.168815E-01,  0.211690E-01,  0.224373E-01,  0.257869E-01,    &
         0.282310E-01,  0.293583E-01,  0.204657E-01,  0.383510E-01/
       data (bpcm_orig(k),k=9,NBLY_orig)    /     &
        -0.560024E-04, -0.645285E-04, -0.559339E-04, -0.456444E-04,    &
        -0.788833E-04, -0.105692E-03, -0.811350E-04, -0.103490E-03/
       data (atpcm_orig(k),k=9,NBLY_orig)    /    &
         0.167628E-01,  0.195779E-01,  0.228872E-01,  0.274006E-01,   &
         0.305678E-01,  0.292862E-01,  0.201527E-01,  0.373711E-01/
       data (btpcm_orig(k),k=9,NBLY_orig)    /     &
        -0.510036E-04, -0.643516E-04, -0.722669E-04, -0.742132E-04,   &
        -0.903238E-04, -0.102388E-03, -0.702017E-04, -0.120756E-03/

!----------------------------------------------------------------------
!    2) 160-560 (as 40 bands). program gasbnd is used with 10 cm-1
!    bandwidth. iband is straightforward mapping.
!----------------------------------------------------------------------
 
       data (acomb_ckd(k),k=1,40)    /      &
         0.849130E+03,  0.135587E+05,  0.286836E+04,  0.169580E+04,  &
         0.208642E+05,  0.126034E+04,  0.109494E+05,  0.335111E+03, &
         0.488806E+04,  0.860045E+04,  0.537333E+03,  0.437769E+04, &
         0.345836E+04,  0.129538E+03,  0.463562E+04,  0.251489E+03,&
         0.256264E+04,  0.485091E+03,  0.889881E+03,  0.116598E+04,&
         0.125244E+03,  0.457264E+03,  0.142197E+03,  0.444551E+03, &
         0.301446E+02,  0.392750E+03,  0.436426E+02,  0.347449E+02, &
         0.612509E+02,  0.142506E+03,  0.103643E+02,  0.721701E+02 , &
         0.315040E+02,  0.941759E+01,  0.540473E+02,  0.350084E+02, &
         0.300816E+02,  0.379020E+01,  0.125727E+02,  0.545869E+01/
       data (bcomb_ckd(k),k=1,40)    /  &
         0.175174E+00,  0.176667E+00,  0.109512E+00,  0.111893E+00,  &
         0.145289E+00,  0.203190E+00,  0.151547E+00,  0.911103E-01,&
         0.151444E+00,  0.850764E-01,  0.756520E-01,  0.100377E+00,&
         0.171557E+00,  0.125429E+00,  0.105915E+00,  0.816331E-01, &
         0.149472E+00,  0.857054E-01,  0.107092E+00,  0.185458E+00, &
         0.753818E-01,  0.108639E+00,  0.123934E+00,  0.178712E+00,  &
         0.833855E-01,  0.119886E+00,  0.133082E+00,  0.935851E-01,  &
         0.156848E+00,  0.166457E+00,  0.162215E+00,  0.114845E+00, &
         0.724304E-01,  0.740525E-01,  0.734090E-01,  0.141319E+00, &
         0.359408E-01,  0.833614E-01,  0.128919E+00,  0.996329E-01/
       data (apcm_ckd(k),k=1,40)    /    &
         0.549325E-02, -0.150653E-02,  0.268788E-02,  0.138495E-01, &
        -0.714528E-03,  0.112319E-01,  0.113418E-02,  0.215116E-01,&
         0.388898E-02,  0.398385E-02,  0.931768E-02,  0.655185E-02,&
         0.735642E-02,  0.190346E-01,  0.104483E-01,  0.917671E-02,&
         0.108668E-01,  0.305797E-02,  0.163975E-01,  0.147718E-01,&
         0.485502E-02,  0.223258E-01,  0.567357E-02,  0.197808E-01, &
         0.245634E-01,  0.116045E-01,  0.269989E-01,  0.176298E-01,&
         0.128961E-01,  0.134788E-01,  0.391238E-01,  0.117165E-01,  &
         0.691808E-02,  0.202443E-01,  0.137798E-01,  0.215153E-01, &
         0.154358E-01,  0.850256E-02,  0.111306E-01,  0.185757E-01/
       data (bpcm_ckd(k),k=1,40)    /   &
        -0.305151E-04, -0.244741E-05, -0.203093E-04, -0.736015E-04,  &
        -0.158662E-04, -0.381826E-04, -0.197166E-04, -0.984160E-04,&
        -0.222455E-04, -0.346880E-04, -0.395593E-04, -0.426165E-04,&
        -0.410312E-04, -0.848479E-04, -0.597304E-04, -0.318474E-04,&
        -0.450295E-04, -0.284497E-04, -0.772035E-04, -0.545821E-04,&
        -0.242272E-04, -0.105653E-03, -0.854473E-05, -0.672510E-04,&
        -0.109627E-03, -0.330446E-04, -0.682675E-04,  0.479154E-04, &
         0.411211E-04, -0.554504E-04, -0.145967E-03, -0.425913E-04, &
         0.413272E-05, -0.531586E-04, -0.429814E-04, -0.847248E-04, &
        -0.733456E-04,  0.403362E-05, -0.389712E-04, -0.531450E-04/
       data (atpcm_ckd(k),k=1,40)    /   &
         0.541966E-02, -0.153876E-02,  0.158818E-02,  0.133698E-01, &
        -0.270746E-02,  0.691660E-02,  0.485749E-05,  0.199036E-01,  &
         0.319826E-02,  0.220802E-02,  0.985921E-02,  0.462151E-02,  &
         0.554947E-02,  0.149315E-01,  0.911982E-02,  0.696417E-02, &
         0.892579E-02,  0.222579E-02,  0.123105E-01,  0.129295E-01, &
         0.745423E-02,  0.203878E-01,  0.597427E-02,  0.170838E-01, &
         0.231443E-01,  0.127692E-01,  0.222678E-01,  0.165331E-01, &
         0.141333E-01,  0.141307E-01,  0.312569E-01,  0.114137E-01, &
         0.126050E-01,  0.242966E-01,  0.145196E-01,  0.224802E-01, &
         0.150399E-01,  0.173815E-01,  0.103438E-01,  0.188690E-01/
       data (btpcm_ckd(k),k=1,40)    /     &
        -0.202513E-04,  0.948663E-06, -0.130243E-04, -0.678688E-04, &
        -0.986181E-05, -0.258818E-04, -0.139292E-04, -0.916890E-04, &
        -0.138148E-04, -0.275165E-04, -0.298451E-04, -0.310005E-04, &
        -0.314745E-04, -0.695971E-04, -0.486158E-04, -0.198383E-04, &
        -0.421494E-04, -0.102396E-04, -0.591516E-04, -0.575729E-04, &
        -0.238464E-04, -0.938688E-04, -0.885666E-05, -0.728295E-04, &
        -0.938897E-04, -0.448622E-04, -0.642904E-04, -0.102699E-04, &
        -0.348743E-05, -0.568427E-04, -0.104122E-03, -0.313943E-04, &
         0.939109E-05, -0.631332E-04, -0.325944E-04, -0.757699E-04,&
        -0.398066E-04,  0.285026E-04, -0.355222E-04, -0.266443E-04/
 
!----------------------------------------------------------------------
!    3) 560-1200 and 4.3 um band (8 bands, frequency range given
!    by bdlocm-bdhicm). program gasbnd is used with 8 specified
!    bandwidths.
!----------------------------------------------------------------------
 
       data (acomb_ckd(k),k=41,NBLY_ckd2p1)    /     &
     &   0.788160E+01,  0.207769E+01,  0.433980E+00,  0.499718E-01, &
     &   0.166671E-01,  0.213073E-01,  0.210228E+00,  0.458804E-02/
       data (bcomb_ckd(k),k=41,NBLY_ckd2p1)    /   & 
     &   0.676888E-01,  0.690856E-01,  0.660847E-01,  0.730849E-01, &
     &   0.654438E-01,  0.848294E-01,  0.852442E-01,  0.242787E+00/
       data (apcm_ckd(k),k=41,NBLY_ckd2p1)    /   &
     &   0.168938E-01,  0.211805E-01,  0.224487E-01,  0.257983E-01, &
     &   0.282393E-01,  0.293678E-01,  0.204784E-01,  0.383510E-01/
       data (bpcm_ckd(k),k=41,NBLY_ckd2p1)    /    &
     &  -0.560515E-04, -0.645812E-04, -0.559626E-04, -0.456728E-04, &
     &  -0.789113E-04, -0.105722E-03, -0.811760E-04, -0.103490E-03/
       data (atpcm_ckd(k),k=41,NBLY_ckd2p1)    /  &
     &   0.167743E-01,  0.195884E-01,  0.228968E-01,  0.274106E-01, &
     &   0.305775E-01,  0.292968E-01,  0.201658E-01,  0.373711E-01/
       data (btpcm_ckd(k),k=41,NBLY_ckd2p1)    /   &
     &  -0.510390E-04, -0.643921E-04, -0.722986E-04, -0.742483E-04,  &
     &  -0.903545E-04, -0.102420E-03, -0.702446E-04, -0.120756E-03/
!---------------------------------------------------------------------

       real, dimension (no_combined_bands)  :: band_no_start_orig,   &
					   band_no_start_ckd2p1, &
	  		                   band_no_end_orig,    &
					   band_no_end_ckd2p1


       data band_no_start_ckd2p1 / 1, 25, 41, 42, 43, 44, 46, 47 /
       data band_no_end_ckd2p1   / 24,40, 41, 42, 43, 45, 46, 47 /

       data band_no_start_orig / 1, 5, 9, 10, 11, 12, 14, 15 /
       data band_no_end_orig   / 4, 8, 9, 10, 11, 13, 14, 15 /

       real, dimension (NBLY_CKD2P1-1) :: cld_indx_table_lwclde, &
			              cld_indx_table_orig

       data cld_indx_table_lwclde /40*1, 2, 2, 2, 3, 4, 5, 6 /

       data cld_indx_table_orig   / 47*1 /



!---------------------------------------------------------------------
!----  read namelist ----

    if (file_exist('input.nml')) then
      unit = open_file ('input.nml', action='read')
      ierr=1; do while (ierr /= 0)
      read (unit, nml=sealw99_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'sealw99_nml')
      enddo
10    call close_file (unit)
    endif

    unit = open_file ('logfile.out', action='append')
    if (get_my_pe() == 0) then
      write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
      write (unit,nml=sealw99_nml)
    endif
    call close_file (unit)

    kmax = size(pref,1) - 1   ! radiation grid size
     ks = 1
     ke = kmax
    call longwave_params_init
    Lw_control%do_ckd2p1          = do_ckd2p1
    Lw_control%do_ch4n2olbltmpint = do_ch4n2olbltmpint
 
    if (do_ckd2p1) then
      Lw_parameters%offset = 32
      NBLY = NBLY_CKD2P1
   else
     Lw_parameters%offset = 0
      NBLY = NBLY_ORIG  
   endif
     call lw_gases_stdtf_init   (pref)
!    call optical_path_init (kmax, latb)
    call optical_path_init 

    call longwave_tables_init (Lw_tables, tabsr, &
                              tab1, tab2, tab3, tab1w, &
                              tab1a, tab2a, tab3a)
 

   if (do_nlte) then
!---------------------------------------------------------------------
!    define pressure-dependent index values used in the infrared
!    radiation code. by this manner, the coefficients are defined
!    at execution time (not dependent on the choice of vertical
!    layers)
!     prnlte : pressure (mb) below which non-LTE code (Nlte.F) affects
!              CO2 transmissivities
!---------------------------------------------------------------------
     prnlte = 0.1
!--------------------------------------------------------------------
!    abort execution if trying to run with modified radiation grid
!    code must be added to properly map plm (on the model grid)
!    to the radiation grid. if simply dropping upper levels, then is
!    fairly easy.  abort here so that problem may be addressed and to
!    prevent uncertified results.
!    solution is likely to be to pass in plm on radiation grid, whatever
!    that is.
!--------------------------------------------------------------------
!    if (ks /= 1) then
!      call error_mesg ( 'co2_source_mod', &
!          'cannot currently execute code with all_level_radiation&
!     & = .false -- must handle nlte calculation properly', &
!             FATAL)
!    endif
      kmin = 1
!     kmax = size(pref,1) - 1

!    pd (:) = pref(:,1)
!    pd8(:) = pref(:,2)
     plm (kmin) = 0.
!   plm8(kmin) = 0.
     do k=kmin+1,kmax
!      plm (k) = 0.5*(pd (k-1) + pd (k))
       plm (k) = 0.5*(pref (k-1,1) + pref (k,1))
!       plm8(k) = 0.5*(pd8(k-1) + pd8(k))
     enddo
!    plm (kmax+1) = pd (kmax+1)
     plm (kmax+1) = pref (kmax+1,1)
!   plm8(kmax+1) = pd8(kmax+1)
!--------------------------------------------------------------------
!   convert pressure specification for bottom (flux) pressure level
!   for nlte calculation into an index (ixprnlte)
!-------------------------------------------------------------------
!    allocate ( plm(kmin:kmax) )
!    call get_std_pressures (plm_out=plm)

!!! CAN THE MODEL TOP BE AT A PRESSURE THAT IS NOT ZERO ??
!!! if not, then must define plm(ks) always to be 0.0
!!! implications for lw_gas tf calcs ??
!      kmax = size(plm,1) - 1

!     ks = 1



!    ixprnlte = KS 
     ixprnlte = 1 
!! ixprnlte in radiation grid coordinates
!    do k=KS+1, KE
!    do k=KS+1, kmax
     do k=ks+1, kmax
!      if ((plm(k) - prnlte) .LT. 0.0) then
!! plm*1.0E-02 = mb
!      if ((plm(k)*1.0E-02 - prnlte) .LT. 0.0) then
       if (plm(k)*1.0E-02  .LT. prnlte) then
!        ixprnlte = k 
         ixprnlte = k-ks+1 
       else
         exit
       endif
     enddo







!---------------------------------------------------------------------
!   allocate and obtain elements of the source function for bands in 
!   the 15 um range (used in nlte)
!---------------------------------------------------------------------
     ioffset = Lw_parameters%offset
     allocate ( c1b7    (NBCO215) )
     allocate ( c2b7    (NBCO215) )
     do n=1,NBCO215 
       cent(n) = 0.5E+00*(Lw_tables%bdlocm(n+8+ioffset) + Lw_tables%bdhicm(n+8+ioffset))
       del (n) = Lw_tables%bdhicm(n+8+ioffset) - Lw_tables%bdlocm(n+8+ioffset)
       c1b7(n) = (3.7412E-05)*cent(n)*cent(n)*cent(n)*del(n) 
       c2b7(n) = (1.4387E+00)*cent(n)
     end do
   endif

!--------------------------------------------------------------------



   call gas_tf_init           (pref     )
!    call lw_gases_stdtf_init   (pref)
     call longwave_clouds_init

!---------------------------------------------------------------------



!---------------------------------------------------------------------

      if (Lw_control%do_ckd2p1) then
        band_no_start = band_no_start_ckd2p1
        band_no_end   = band_no_end_ckd2p1
        NBLY = NBLY_CKD2P1
      else
        band_no_start = band_no_start_orig  
        band_no_end   = band_no_end_orig  
        NBLY = NBLY_ORIG 
      endif

      Lw_parameters%NBLY = NBLY

!!!!!!&&&&&&&&&&&&&&&^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!! THIS IS TEMPORARY UNTIL PROPER ORDER MAY BE MADE!!!!!!
!   NLWCLDB = Lw_parameters%NLWCLDB
    NLWCLDB = 7
!   NLWCLDB = 1

!     if (Lw_parameters%NLWCLDB == 7) then
      if (              NLWCLDB == 7) then
        cld_indx_table = cld_indx_table_lwclde
      else
        cld_indx_table = cld_indx_table_orig
      endif
!!!!!!&&&&&&&&&&&&&&&^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      allocate (apcm(NBLY))
      allocate (bpcm(NBLY))
      allocate (atpcm(NBLY))
      allocate (btpcm(NBLY))
      allocate (acomb(NBLY))
      allocate (bcomb(NBLY))

      if (Lw_control%do_ckd2p1) then
	apcm = apcm_ckd
	atpcm = atpcm_ckd
	bpcm = bpcm_ckd
	btpcm = btpcm_ckd
	acomb = acomb_ckd
	bcomb = bcomb_ckd
      else
	apcm = apcm_orig
	atpcm = atpcm_orig
	bpcm = bpcm_orig
	btpcm = btpcm_orig
	acomb = acomb_orig
	bcomb = bcomb_orig
      endif

!------------------------------------------------------------------

!--------------------------------------------------------------------
!    define the cloud band index to be used in the longwave transmission
!    calculations for each cloud band. 
!--------------------------------------------------------------------
    NBTRGE = Lw_parameters%NBTRGE




!!!!!!&&&&&&&&&&&&&&&^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!! THIS IS TEMPORARY UNTIL PROPER ORDER MAY BE MADE!!!!!!
!   NLWCLDB = Lw_parameters%NLWCLDB
    NLWCLDB = 7
!   NLWCLDB = 1
!!!!!!&&&&&&&&&&&&&&&^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




!   if (6+NBTRGE > max_cld_bands) then
!     call error_mesg ( 'longwave_driver_mod', &
!      ' max_cld_bands < 6 + NBTRGE; reset one of these', FATAL)
!   endif
    allocate (cld_indx(6+NBTRGE) )
    do nn=1,6+NBTRGE       
!      cld_band_no(nn) = nn
      if (nn > NLWCLDB) then
        cld_indx(nn) = NLWCLDB
      else
        cld_indx(nn) = nn
      endif
    end do


!--------------------------------------------------------------------
   
    call longwave_fluxes_init

end subroutine sealw99_init


!#####################################################################

subroutine sealw99 (is, ie, js, je, Atmos_input, Rad_gases, &
                    Aerosol, Cldrad_props, Lw_output, Lw_diagnostics)

integer, intent(in) :: is, ie, js, je
type(atmos_input_type), intent(in) :: Atmos_input  
type(radiative_gases_type), intent(in) :: Rad_gases   
type(aerosol_type), intent(in) :: Aerosol      
type(cldrad_properties_type), intent(in) :: Cldrad_props
type(lw_output_type), intent(inout) :: Lw_output   
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics

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

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3), NLWCLDB) ::    &
                                                                cldtf

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3)) ::    &
                                            cnttaub1, cnttaub2,  &
	 				    cnttaub3, co21c, co21diag, &
					    co21r, dsorc15, dsorc93, &
					    dsorcb1, dsorcb2, dsorcb3, &
 					    dt4, emiss, heatem, overod,&
					    sorc15, t4, tmp1, tmp2, &
					    to3cnt, ch41c, n2o1c, n2o17c

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), 2)  ::  emspec
      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2))   ::  &
            s1a, flxge1, co21c_KEp1, co21r_KEp1   
      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3), 3) ::    &
                                contdg
      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3) -1) ::    &
                               cfc_tf, pdflux

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2))   ::  &
                                         flx1e1cf, flxge1cf, gxctscf

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3)) ::    &
                                                              emissb, &
                            e1cts1, e1cts2, e1ctw1, e1ctw2,        &
                            soe2, soe3, soe4, soe5, emisdg, flxcf, &
                            heatemcf, flx, to3dg, taero8, taero8kp,   &
                            totaer_tmp, tcfc8

                       
      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3) - 1) ::    &
                         cts_sum, cts_sumcf

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2),    &
                       size(Atmos_input%press,3), NBTRGE) ::    &
                                emissbf, e1cts1f, e1cts2f,         &
                                 emisdgf, emissf, tch4n2oe

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), NBTRGE) ::    &
                               flx1e1fcf,  flxge1f, flxge1fcf        


      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
		       size(Atmos_input%press,3),   &
                           8):: sorc

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
		       size(Atmos_input%press,3),   &
		          6+NBTRGE     ) :: source_band, dsrcdp_band

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
		       size(Atmos_input%press,3) ,  &
		          6+NBTRGE     ) ::  &
			      trans_band1, trans_band2

      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), &
		          6+NBTRGE     ) ::  &
			      trans_b2d1, trans_b2d2


      real, dimension (size(Atmos_input%press,1),    &
                       size(Atmos_input%press,2), 2, NBTRGE) ::    &
                                      emspecf

   integer   ::          n, k, kp, m, j
     integer  :: kx, ix, jx
   real, dimension (size(Atmos_input%press,1),    &
                 size(Atmos_input%press,2),  &     
                                size(Atmos_input%press,3)-1 )  :: &
                                                    pdfinv

     type(lw_clouds_type) :: Lw_clouds
     type(optical_path_type) :: Optical
     type(gas_tf_type) :: Gas_tf   
      logical             ::  calc_co2, calc_n2o, calc_ch4
      real                ::  ch4_vmr, n2o_vmr, co2_vmr

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      call get_control_gas_tf (calc_co2, calc_ch4, calc_n2o)

      call lw_gases_stdtf_time_vary

!       if (do_ch4_n2o .and. ( time_varying_ch4        .or. &
!                                  time_varying_n2o        ) ) then
!       if (time_varying_ch4 .or. time_varying_n2o) then
        if (Rad_gases%time_varying_ch4 .or. Rad_gases%time_varying_n2o) then
           if (Rad_gases%time_varying_ch4        .and. .not. calc_ch4) then
          call error_mesg ('radiative_gases_mod', &
          ' if ch4 amount is to vary in time, ch4 tfs must be &
	    &recalculated', FATAL)
           endif
           if (Rad_gases%time_varying_n2o        .and. .not. calc_n2o) then
          call error_mesg ('radiative_gases_mod', &
          ' if n2o amount is to vary in time, n2o tfs must be &
	    &recalculated', FATAL)
           endif

        endif



!       if (do_ch4_n2o .and. do_ch4_n2o_init) then
        if (Lw_control%do_ch4_n2o .and. (do_ch4_n2o_init .or.    &
	                      Rad_gases%time_varying_ch4 .or. &
			      Rad_gases%time_varying_n2o) ) then
!          call ch4_n2o_time_vary (rrvch4, rrvn2o)
           call ch4_n2o_time_vary (Rad_gases%rrvch4, Rad_gases%rrvn2o)
           do_ch4_n2o_init = .false.
       endif

!       if (do_co2 .and. ( time_varying_co2       )) then 
        if (Rad_gases%time_varying_co2) then 



!       if (time_varying_co2  .and. .not. calc_co2) then
        if (.not. calc_co2) then
          call error_mesg ('radiative_gases_mod', &
          ' if co2 amount is to vary in time, co2 tfs must be &
	    &recalculated', FATAL)
        endif



        endif


        if (Lw_control%do_co2 .and. (do_co2_init  .or.  &
	                  Rad_gases%time_varying_co2) ) then
        call co2_time_vary (Rad_gases%rrvco2)
           do_co2_init = .false.
       endif





		       
      if (calc_co2 .or. calc_ch4 .or. calc_n2o) then
	call lw_gases_stdtf_dealloc
      endif

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


     kx = size(Atmos_input%press,3) - 1
     ix = ie-is+1
     jx = je-js+1

!--------------------------------------------------------------------
!   call sealw99_alloc to allocate component arrays of 
!   lw_diagnostics_type variables.
!----------------------------------------------------------------------
     call sealw99_alloc (ix, jx, kx,             &
                                Lw_diagnostics)

!      do_diagnostics = .false.
!   do j=js,je
!     if (Rad_control%do_raddg(j)) then
!          do_diagnostics = .true.
!          exit
!      endif
!  end do



     call optical_path_setup (is, ie, js, je,     Atmos_input  , &
                               Rad_gases, Aerosol, Optical)


!----------------------------------------------------------------------
!     call cldtau to compute cloud layer transmission functions for 
!     all layers.
!-------------------------------------------------------------------
      call cldtau (Cldrad_props, Lw_clouds)

!--------------------------------------------------------------------
!     call co2coef to compute some co2 temperature and pressure   
!     interpolation quantities and to compute temperature-corrected co2
!     transmission functions (co2spnb and co2nbl). 
!-------------------------------------------------------------------
     call co2coef (Atmos_input, Gas_tf)

!----------------------------------------------------------------------
!    call co2_source_calc to calculate the source function.
!-----------------------------------------------------------------------
     call co2_source_calc (Atmos_input, Rad_gases, sorc, Gas_tf, &
                           source_band, dsrcdp_band)

!---------------------------------------------------------------------
!   BEGIN LEVEL KS CALCULATIONS
!---------------------------------------------------------------------


!-----------------------------------------------------------------
!     compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!     mission functions and n2o 560-670 cm-1 transmission functions
!     appropriate for level KS. 
!---------------------------------------------------------------------
     call transcolrow (Gas_tf, KS, KS, KS, KE+1, KS+1, KE+1,   &
                       co21c, co21r, tch4n2oe)

!---------------------------------------------------------------------
!     go into optical_path_mod to obtain the optical path functions 
!     needed for use from level KS.
!     to3cnt contains values in the 990 - 1070 cm-1 range (ozone, water
!     vapor continuum, aerosol, cfc).
!     overod contains values in the 560 - 800 cm-1 range (water vapor 
!     lines, continuum, aerosol, 17 um n2o band, cfc).
!     cnttaub1 is continuum band 4 (water vapor continuum, aerosol, cfc)
!     cnttaub2 is continuum band 5 (water vapor continuum, aerosol, cfc)
!     cnttaub3 is continuum band 7 (water vapor continuum, aerosol, cfc)
!     the 15um band transmission functions between levels KS and KS+1

!     are stored in overod and co2nbl; they will not be overwritten,
!     as long as calculations are made for pressure levels increasing
!     from KS.
!---------------------------------------------------------------------
   call optical_trans_funct_from_KS (Gas_tf, to3cnt, overod,  &
                                     Optical, cnttaub1,    &
				     cnttaub2, cnttaub3)


!-----------------------------------------------------------------------
!    compute cloud transmission functions between level KS and all
!    other levels.
!-----------------------------------------------------------------------
   call cloud (KS, Cldrad_props, Lw_clouds, cldtf)


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


     call cool_to_space_exact (                cldtf, Atmos_input,   &
                          Optical,  Gas_tf,  sorc,        &
                            to3cnt, Lw_diagnostics,   &
			     cts_sum, cts_sumcf, gxctscf)

 
!----------------------------------------------------------------------
!     compute the emissivity fluxes for k=KS.
!----------------------------------------------------------------------


     call e1e290 (Atmos_input,  e1ctw1, e1ctw2,  &
                   trans_band1, trans_band2,  Optical, tch4n2oe, &
                   source_band(:,:,:,1), Lw_diagnostics, cldtf, &
                   cld_indx, flx1e1cf,  tcfc8)

 
!! trans_band1:
!   index 1 = e1flx
!   index 2 = co21r*overod
!   index 3 = cnttaub1
!   index 4 = cnttaub2
!   index 5 = to3cnt
!   index 6 = cnttaub3
!   index 7 = e1flxf




     trans_band1(:,:,KS:KE+1,2) = co21r   (:,:,KS:KE+1)*  &
                                               overod(:,:,KS:KE+1)
     trans_band1(:,:,KS:KE+1,3) = cnttaub1(:,:,KS:KE+1)
     trans_band1(:,:,KS:KE+1,4) = cnttaub2(:,:,KS:KE+1)
     trans_band1(:,:,KS:KE+1,5) = to3cnt  (:,:,KS:KE+1)
     trans_band1(:,:,KS:KE+1,6) = cnttaub3(:,:,KS:KE+1)

!! trans_band2:
!   index 1 = emiss
!   index 2 = co21c*overod
!   index 3 = cnttaub1
!   index 4 = cnttaub2
!   index 5 = to3cnt
!   index 6 = cnttaub3
!   index 7 = emissf


     trans_band2(:,:,KS+1:KE+1,2) = co21c(:,:,KS+1:KE+1)*   &
                                       overod(:,:,KS+1:KE+1)
!    trans_band2(:,:,KS+1:KE+1,3) = trans_band1(:,:,KS+1:KE+1,3)
!    trans_band2(:,:,KS+1:KE+1,4) = trans_band1(:,:,KS+1:KE+1,4)
!    trans_band2(:,:,KS+1:KE+1,5) = trans_band1(:,:,KS+1:KE+1,5)
!    trans_band2(:,:,KS+1:KE+1,6) = trans_band1(:,:,KS+1:KE+1,6)
     trans_band2(:,:,KS+1:KE+1,3:6) = trans_band1(:,:,KS+1:KE+1,3:6)

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
    call longwave_fluxes_ks (source_band, trans_band1, dsrcdp_band,  &
                             trans_band2, cldtf, cld_indx , &
			                  Lw_diagnostics)  

 
     trans_band1(:,:,KS:KE+1,1) = e1ctw2(:,:,KS:KE+1)

!-----------------------------------------------------------------------
!     compute approximate cool-to-space heating rates for 1 wide band
!     in the 15um  range (560-800 cm-1) (ctsco2) and for 1 band in 
!     the 9.6 um band (ctso3).
!----------------------------------------------------------------------

    call cool_to_space_approx ( Atmos_input%pflux,    &
    source_band, trans_band1,       &
                               cldtf, cld_indx ,  &
                                          Lw_diagnostics, e1ctw1)

!-----------------------------------------------------------------------
!     perform flux calculations for the flux levels KS+1 to KE-1. calcu-
!     lations for flux levels KE and KE+1 are done separately, as all
!     calculations are special cases or nearby layers.
!----------------------------------------------------------------------
    do k=KS+1,KE-1

!--------------------------------------------------------------------
!     compute co2 560-800 cm-1, ch4 and n2o 1200-1400 cm-1 trans-
!     mission functions and n2o 560-670 cm-1 transmission functions
!     appropriate for level k. 
!--------------------------------------------------------------------
        call transcolrow (Gas_tf, k, k, k, KE+1, k+1, KE+1,  &
                          co21c, co21r, tch4n2oe)
 

!-------------------------------------------------------------------
!     the 15 um band transmission functions between levels k and k+1
!     are stored in overod and co2nbl; they will not be overwritten,
!     as long as calculations are made for pressure levels increasing
!     from k.
!---------------------------------------------------------------------
      call optical_trans_funct_k_down (Gas_tf, k,                     &
                                             to3cnt, overod, Optical)



!-----------------------------------------------------------------------
!     compute cloud transmission functions between level k and all
!     other levels greater or equal to k.
!---------------------------------------------------------------------
      call cloud (k,Cldrad_props, Lw_clouds, cldtf)


!-----------------------------------------------------------------------
!     compute the exchange terms in the flux equation (except the 
!     nearby layer (k,k) terms, done later).
!---------------------------------------------------------------------
        call e290 (Atmos_input, k, trans_band2, trans_band1,   &
                   Optical,  tch4n2oe,  tcfc8)

!! trans_band1:
!   index 1 = emissb
!   index 2 = co21r*overod
!   index 3 = contodb1
!   index 4 = contodb2
!   index 5 = to3cnt
!   index 6 = contodb3
!   index 7 = emissbf

       do kp=k,KE
      trans_band1(:,:,kp+1,3) = cnttaub1(:,:,kp+1  )/cnttaub1(:,:,k  )
      trans_band1(:,:,kp+1,4) = cnttaub2(:,:,kp+1  )/cnttaub2(:,:,k  )
      trans_band1(:,:,kp+1,6) = cnttaub3(:,:,kp+1  )/cnttaub3(:,:,k  )
       end do


     trans_band1(:,:,k+1:KE+1,2) = co21r(:,:,k+1:KE+1)*   &
                                       overod(:,:,k+1:KE+1)
     trans_band1(:,:,k+1:KE+1,5) = to3cnt(:,:,k+1:KE+1)


!! trans_band2:
!   index 1 = emiss
!   index 2 = co21c*overod
!   index 3 = contodb1
!   index 4 = contodb2
!   index 5 = to3cnt
!   index 6 = contodb3
!   index 7 = emissf

     trans_band2(:,:,k+1:KE+1,2) = co21c(:,:,k+1:KE+1)*   &
                                            overod(:,:,k+1:KE+1)
!    trans_band2(:,:,k+1:KE+1,3) = trans_band1(:,:,k+1:KE+1,3)
!    trans_band2(:,:,k+1:KE+1,3) = trans_band1(:,:,k+1:KE+1,3)
!    trans_band2(:,:,k+1:KE+1,4) = trans_band1(:,:,k+1:KE+1,4)
!    trans_band2(:,:,k+1:KE+1,5) = trans_band1(:,:,k+1:KE+1,5)
     trans_band2(:,:,k+1:KE+1,3:6) = trans_band1(:,:,k+1:KE+1,3:6)

!-----------------------------------------------------------------------
!     compute the terms for flux at levels k+1 to KE+1 from level k.
!     compute the terms for flux at level k due to levels 
!     kp from k+1 to KE+1.
!----------------------------------------------------------------------
      call longwave_fluxes_k_down  (k, dsrcdp_band,   &
               trans_band1, trans_band2,    &
                                    cldtf, cld_indx,  &
                                     Lw_diagnostics          )  


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





      call transcolrow (Gas_tf, KE, KE, KE, KE+1, KE+1, KE+1,  &
                        co21c, co21r, tch4n2oe)

!----------------------------------------------------------------------
!     get optical path terms for KE
!----------------------------------------------------------------------
    call optical_trans_funct_KE (Gas_tf, to3cnt, Optical, overod)

   trans_b2d1(:,:,3  ) = cnttaub1(:,:,KE+1)/cnttaub1(:,:,KE  )
   trans_b2d1(:,:,4  ) = cnttaub2(:,:,KE+1)/cnttaub2(:,:,KE  )
   trans_b2d1(:,:,6  ) = cnttaub3(:,:,KE+1)/cnttaub3(:,:,KE  )

!-----------------------------------------------------------------------
!     compute cloud transmission functions between level KE and KE and
!     KE+1
!----------------------------------------------------------------------
    call cloud (KE, Cldrad_props, Lw_clouds, cldtf)

   
!-------------------------------------------------------------------- 
!     compute mean temperature in the "nearby layer" between a flux
!     level and the first data level below the flux level (tpl1) or the
!     first data level above the flux level (tpl2)
!---------------------------------------------------------------------


      call esfc  (Atmos_input, emspec, Optical, &
                           emspecf, tch4n2oe, tcfc8)

!----------------------------------------------------------------------
!     compute nearby layer transmission functions for 15 um band, cont-
!     inuum bands, and 9.3 um band in subroutine Nearbylyrtf. trans-
!     mission functions for the special cases (KE,KE+1) and (KE+1,KE)
!     are also computed for the 15 um band.
!----------------------------------------------------------------------

    call trans_sfc    (Gas_tf, Atmos_input, overod, co21c_KEp1, &
                       co21r_KEp1)



     trans_b2d1(:,:,1) = emspec(:,:,KS+1)
     trans_b2d1(:,:,2) = co21c_KEp1(:,:)
     trans_b2d1(:,:,5) = to3cnt(:,:,KE+1)

     do m=1,NBTRGE
     trans_b2d1(:,:,6+m) = emspecf(:,:,KS+1,m)
     end do

     trans_b2d2(:,:,1) = emspec(:,:,KS)
     trans_b2d2(:,:,2) = co21r_KEP1(:,:)
!    trans_b2d2(:,:,3) = trans_b2d1(:,:,3)
!    trans_b2d2(:,:,4) = trans_b2d1(:,:,4)
!    trans_b2d2(:,:,5) = to3cnt(:,:,KE+1)
!    trans_b2d2(:,:,6) = trans_b2d1(:,:,6)
     trans_b2d2(:,:,3:6) = trans_b2d1(:,:,3:6)
     do m=1,NBTRGE
     trans_b2d2(:,:,6+m) = emspecf(:,:,KS,m)
     end do

!! trans_band1:
!   index 1 = emspec(KS+1)
!   index 2 = co21c_KEp1  
!   index 3 = contodb1
!   index 4 = contodb2
!   index 5 = to3cnt
!   index 6 = contodb3
!   index 7 = emspecf(KS+1)

!! trans_band2:
!   index 1 = emspec(KS)
!   index 2 = co21r_KEp1   
!   index 3 = contodb1
!   index 4 = contodb2
!   index 5 = to3cnt
!   index 6 = contodb3
!   index 7 = emspecf(KS)

!-----------------------------------------------------------------------
!     obtain fluxes for the two terms (KE,KE+1) and (KE+1,KE), both 
!     using the same cloud transmission functions (from layer KE)
!----------------------------------------------------------------------
    call longwave_fluxes_KE_KEp1 (dsrcdp_band, &
            trans_b2d1, trans_b2d2,    &
				  cldtf, cld_indx,  &
				                  Lw_diagnostics )

!---------------------------------------------------------------------
!     call enear to calculate emissivity arrays
!----------------------------------------------------------------------
      call enear (Atmos_input, emisdg,                     Optical, &
                  emisdgf , tch4n2oe, tcfc8        )

!-------------------------------------------------------------------
!     obtain optical path transmission functions for diagonal terms
!------------------------------------------------------------------
    call optical_trans_funct_diag (Atmos_input, contdg, to3dg, Optical)
 
!-----------------------------------------------------------------------
!     compute cloud transmission functions between level KE+1 and KE+1
!------------------------------------------------------------------
    call cloud (KE+1, Cldrad_props, Lw_clouds, cldtf)
 
!----------------------------------------------------------------------
!     compute nearby layer transmission functions for 15 um band, cont-
!     inuum bands, and 9.3 um band in subroutine Nearbylyrtf. trans-
!     mission functions for the special cases (KE,KE+1) and (KE+1,KE)
!     are also computed for the 15 um band.
!----------------------------------------------------------------------
     call trans_nearby (Gas_tf, Atmos_input, overod,  co21c)


     trans_band1(:,:,ks+1:KE+1,1) = emisdg(:,:,ks+1:KE+1)
     trans_band1(:,:,ks+1:KE+1,2) = co21c(:,:,ks+1:KE+1)
     trans_band1(:,:,ks+1:KE+1,3) = contdg(:,:,ks+1:KE+1,1)
     trans_band1(:,:,ks+1:KE+1,4) = contdg(:,:,ks+1:KE+1,2)
     trans_band1(:,:,ks+1:KE+1,5) = to3dg(:,:,ks+1:KE+1)
     trans_band1(:,:,ks+1:KE+1,6) = contdg(:,:,ks+1:KE+1,3)


     do m=1,NBTRGE
     trans_band1(:,:,ks+1:KE+1,6+m) = emisdgf(:,:,ks+1:KE+1,m)
     end do

!-----------------------------------------------------------------------
!     obtain fluxes for the diagonal terms at all levels.
!-----------------------------------------------------------------------
    call longwave_fluxes_diag (dsrcdp_band,           &
                          trans_band1,   &
			       cldtf , cld_indx,   &
			                       Lw_diagnostics )


!--------------------------------------------------------------------
!      sum up fluxes over bands
!-----------------------------------------------------------------------

    if (Rad_control%do_totcld_forcing) then
      call longwave_fluxes_sum (is, ie, js, je, flx, NBTRGE, Lw_diagnostics, flxcf)
    else
      call longwave_fluxes_sum (is, ie, js, je, flx, NBTRGE,Lw_diagnostics)
    endif

!-----------------------------------------------------------------------
!     compute emissivity heating rates.
!-----------------------------------------------------------------------

     pdfinv(:,:,ks:ke) = 1.0/(Atmos_input%pflux(:,:,ks+1:ke+1) -  &
                            Atmos_input%pflux(:,:,ks:ke))



!6  heatem(:,:,KS:KE) = radcon_mks*(flx(:,:,KS+1:KE+1) -    &
!6                      flx(:,:,KS:KE))*pdfinv(:,:,KS:KE)
    heatem(:,:,KS:KE) = (radcon_mks*(flx(:,:,KS+1:KE+1) -    &
                        flx(:,:,KS:KE))*pdfinv(:,:,KS:KE))*1.0e-03
    if (Rad_control%do_totcld_forcing) then 		   
!6    heatemcf(:,:,KS:KE) = radcon_mks*(flxcf(:,:,KS+1:KE+1) -    &
!6                          flxcf(:,:,KS:KE))*pdfinv(:,:,KS:KE)
      heatemcf(:,:,KS:KE) = (radcon_mks*(flxcf(:,:,KS+1:KE+1) -    &
                            flxcf(:,:,KS:KE))*pdfinv(:,:,KS:KE)*1.0e-03)
    endif


!-----------------------------------------------------------------------
!     compute total heating rates.
!-----------------------------------------------------------------------
!--------------------------------------------------------------------
!     cts_sum is the sum of the values from cool_to_space_exact and
!     the values defined here in cool_to_space_approx. it will be used
!     by longwave_driver_mod.
!--------------------------------------------------------------------
    cts_sum(:,:,:) = ((((((cts_sum(:,:,:) -   &
               Lw_diagnostics%cts_out(:,:,:,2)) -  & 
               Lw_diagnostics%cts_out(:,:,:,5) )-  &
               Lw_diagnostics%cts_out(:,:,:,1))  - &
               Lw_diagnostics%cts_out(:,:,:,3)) - &
               Lw_diagnostics%cts_out(:,:,:,4))- &
               Lw_diagnostics%cts_out(:,:,:,6))

    if (Rad_control%do_totcld_forcing) then
    cts_sumcf(:,:,:) = ((((((cts_sumcf(:,:,:) -  &
                            Lw_diagnostics%cts_outcf(:,:,:,2)) - &
                            Lw_diagnostics%cts_outcf(:,:,:,5)) - &
                            Lw_diagnostics%cts_outcf(:,:,:,1)) - &
                            Lw_diagnostics%cts_outcf(:,:,:,3)) - &
                            Lw_diagnostics%cts_outcf(:,:,:,4) ) - &
                            Lw_diagnostics%cts_outcf(:,:,:,6))
    endif

!6    Lw_output%heatra(:,:,KS:KE) = 1.0e-03*(heatem(:,:,KS:KE) +   &
!6             cts_sum  (:,:,KS:KE) )
      Lw_output%heatra(:,:,KS:KE) =          heatem(:,:,KS:KE) +   &
               cts_sum  (:,:,KS:KE)  

    if (Rad_control%do_totcld_forcing) then 		   
!6    Lw_output%heatracf(:,:,KS:KE) = 1.0E-03*(heatemcf(:,:,KS:KE) +   &
!6                         cts_sumcf(:,:,KS:KE)) 
      Lw_output%heatracf(:,:,KS:KE) =          heatemcf(:,:,KS:KE) +   &
                           cts_sumcf(:,:,KS:KE) 
    endif

!-----------------------------------------------------------------------
!     compute the flux at each flux level using the flux at the
!     top (flx1e1 + gxcts) and the integral of the heating rates 
!-----------------------------------------------------------------------
      pdflux(:,:,KS:KE) = Atmos_input%pflux(:,:,KS+1:KE+1) -   &
                              Atmos_input%pflux(:,:,KS:KE)

      Lw_output%flxnet(:,:,KS   ) = Lw_diagnostics%flx1e1(:,:) + Lw_diagnostics%gxcts(:,:)



! convert values to mks (1.0e-03 factor) 
      Lw_diagnostics%gxcts(:,:) = 1.0e-03*Lw_diagnostics%gxcts(:,:)
      Lw_diagnostics%flx1e1(:,:) = 1.0e-03*Lw_diagnostics%flx1e1(:,:)
      Lw_diagnostics%flx1e1f(:,:,:) =    &
                    1.0e-03*Lw_diagnostics%flx1e1f(:,:,:)


! convert mks values to cgs (1.0e03 factor) so can be summed with cgs value

      tmp1 (:,:,KS:KE) = 1.0e03*Lw_output%heatra(:,:,KS:KE)*pdflux(:,:,KS:KE)/radcon_mks
    do k=KS+1,KE+1
      Lw_output%flxnet(:,:,k) = Lw_output%flxnet(:,:,k-1) + tmp1(:,:,k-1)
    enddo

   if (Rad_control%do_totcld_forcing) then 		   
     Lw_output%flxnetcf(:,:,KS   ) = flx1e1cf(:,:) + gxctscf(:,:)
! convert mks values to cgs (1.0e03 factor) so can be summed with cgs value
     tmp1 (:,:,KS:KE) = 1.0e03*Lw_output%heatracf(:,:,KS:KE)*pdflux(:,:,KS:KE)/radcon_mks 
     do k=KS+1,KE+1
       Lw_output%flxnetcf(:,:,k) = Lw_output%flxnetcf(:,:,k-1) + tmp1(:,:,k-1)
     enddo
   endif

!-----------------------------------------------------------------------
!    call thickcld to perform "pseudo-convective adjustment" for
!    maximally overlapped clouds, if desired.
!-----------------------------------------------------------------------
   if (do_thick) then

       call thickcld (  Atmos_input%pflux, Cldrad_props, Lw_output)
   endif  ! (do_thick)


!--------------------------------------------------------------------
!   convert lw fluxes to mks units.
!---------------------------------------------------------------------
     Lw_output%flxnet(:,:,:) = 1.0E-03*Lw_output%flxnet(:,:,:)
     if (Rad_control%do_totcld_forcing) then
        Lw_output%flxnetcf(:,:,:) = 1.0E-03*Lw_output%flxnetcf(:,:,:)
      endif

!---------------------------------------------------------------------
!    call subroutine to deallocate components of stack resident derived
!    type variables.
!---------------------------------------------------------------------
      call deallocate_arrays (Lw_clouds, Optical, Gas_tf)
      
end subroutine sealw99 


!#####################################################################

subroutine deallocate_arrays (Lw_clouds, Optical, Gas_tf)

type(lw_clouds_type),    intent(in)       :: Lw_clouds
type(optical_path_type), intent(in)       :: Optical
type(gas_tf_type),       intent(in)       :: Gas_tf   

!--------------------------------------------------------------------
!    deallocate component arrays of Lw_clouds.
!--------------------------------------------------------------------
      deallocate (Lw_clouds%taucld_rndlw)
      deallocate (Lw_clouds%taucld_mxolw)
      deallocate (Lw_clouds%taunbl_mxolw)

!--------------------------------------------------------------------
!    deallocate component arrays of Gas_tf.
!--------------------------------------------------------------------
      deallocate (Gas_tf%tdav)
      deallocate (Gas_tf%tlsqu         )
      deallocate (Gas_tf%tmpdiff       )
      deallocate (Gas_tf%tstdav        )
      deallocate (Gas_tf%co2nbl        )
      deallocate (Gas_tf%n2o9c         )
      deallocate (Gas_tf%tn2o17        )
      deallocate (Gas_tf%co2spnb       )
      deallocate (Gas_tf%a1            )
      deallocate (Gas_tf%a2            )

!--------------------------------------------------------------------
!    deallocate component arrays of Optical.
!--------------------------------------------------------------------
      deallocate (Optical%empl1          )
      deallocate (Optical%empl2          )
      deallocate (Optical%var1           )
      deallocate (Optical%var2           )
      deallocate (Optical%avephi         )
      deallocate (Optical%totphi         )
      deallocate (Optical%emx1           )
      deallocate (Optical%emx2           )

      if (Lw_control%do_ch4_n2o) then
        deallocate (Optical%avephif        )
        deallocate (Optical%emx1f          )
        deallocate (Optical%emx2f          )
        deallocate (Optical%empl1f         )
        deallocate (Optical%empl2f         )
        deallocate (Optical%vrpfh2o        )
        deallocate (Optical%tphfh2o         )
      endif

      if (Lw_control%do_ckd2p1) then
        deallocate (Optical%xch2obd        )
        deallocate (Optical%totch2obdwd    )
        deallocate (Optical%xch2obdwd      )
      else
        deallocate (Optical%cntval         )
        deallocate (Optical%totvo2         )
      endif

      deallocate (Optical%toto3          )
      deallocate (Optical%tphio3         )
      deallocate (Optical%var3           )
      deallocate (Optical%var4           )
      deallocate (Optical%wk             )
      deallocate (Optical%rh2os          )
      deallocate (Optical%rfrgn          )
      deallocate (Optical%tfac           )

      if (Lw_control%do_cfc) then
        deallocate (Optical%totf11         )
        deallocate (Optical%totf12         )
        deallocate (Optical%totf113         )
        deallocate (Optical%totf22         )
      endif

      if(Lw_control%do_lwaerosol) then
        deallocate (Optical%totaerooptdep  )
        deallocate (Optical%aerooptdep_KE_15  )
      endif

end subroutine deallocate_arrays 



!#####################################################################
subroutine sealw99_end                  

      call gas_tf_end

end subroutine sealw99_end                  


!####################################################################

subroutine cool_to_space_approx (     pflux_in,        source,  &
                                 trans,      cld_trans, cld_ind, &
                                 Lw_diagnostics, &
                                 trans2      )

!---------------------------------------------------------------------
integer, dimension (:),  intent(in)           ::  cld_ind         
real, dimension (:,:,:),  intent(in)           ::  pflux_in        
real, dimension (:,:,:,:),  intent(in)           ::  source, trans, &
					           cld_trans
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
real, dimension (:,:,:),  intent(in), optional ::  trans2

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
    integer                               ::  j
    integer  :: index, nbands
    real, dimension(size(pflux_in,1), size(pflux_in,2), &
                    size(pflux_in,3)-1) :: pdfinv

     nbands = size(source,4) - NBTRGE


     pdfinv(:,:,KS:KE) = 1.0/(pflux_in(:,:,KS+1:KE+1) - pflux_in(:,:,KS:KE))
!---------------------------------------------------------------------

     do index=1,nbands

    if (index == 1     ) then
!6    Lw_diagnostics%cts_out(:,:,KS  :KE,index) = radcon_mks*pdfinv(:,:,KS  :KE)*     &
      Lw_diagnostics%cts_out(:,:,KS  :KE,index) = (radcon_mks*pdfinv(:,:,KS  :KE)*     &
			     source(:,:,KS  :KE,index)*  &
                             (trans(:,:,KS      :KE ,index)*   &
	   		     cld_trans(:,:,KS+1:KE+1, cld_ind(index)) -   &
                             trans2(:,:,KS        :KE)*     &
!6	     cld_trans(:,:,KS  :KE, cld_ind(index)))
			     cld_trans(:,:,KS  :KE, cld_ind(index)))*1.0e-03)
    else
!6    Lw_diagnostics%cts_out(:,:,KS  :KE,index) = radcon_mks*pdfinv(:,:,KS  :KE)*     &
      Lw_diagnostics%cts_out(:,:,KS  :KE,index) = (radcon_mks*pdfinv(:,:,KS  :KE)*     &
			     source(:,:,KS  :KE, index)*  &
                             (trans(:,:,KS+1    :KE+1,index    )*   &
			     cld_trans(:,:,KS+1:KE+1, cld_ind(index)) -   &
                             trans(:,:,KS       :KE,index    )*     &
!6	     cld_trans(:,:,KS  :KE, cld_ind(index)))
			     cld_trans(:,:,KS  :KE, cld_ind(index)))*1.0e-03)
    endif

    if (Rad_control%do_totcld_forcing) then

    if (index == 1     ) then
!6      Lw_diagnostics%cts_outcf(:,:,KS  :KE,index) = radcon_mks*pdfinv(:,:,KS  :KE)*     &
        Lw_diagnostics%cts_outcf(:,:,KS  :KE,index) = (radcon_mks*pdfinv(:,:,KS  :KE)*     &
			         source(:,:,KS  :KE,index)* &
                                 (trans(:,:,KS      :KE,index      ) -     &
!6	          trans2(:,:,KS       :KE     ))
			          trans2(:,:,KS       :KE     ))*1.0e-03)
      else
!6      Lw_diagnostics%cts_outcf(:,:,KS  :KE,index) = radcon_mks*pdfinv(:,:,KS  :KE)*     &
        Lw_diagnostics%cts_outcf(:,:,KS  :KE,index) = (radcon_mks*pdfinv(:,:,KS  :KE)*     &
			         source(:,:,KS  :KE,index)* &
                                 (trans(:,:,KS+1    :KE+1, index    ) -     &
!6	          trans(:,:,KS      :KE, index    ))
			          trans(:,:,KS      :KE, index    ))*1.0e-03)
      endif
    endif

   end do  ! (index loop)


!--------------------------------------------------------------------

end subroutine cool_to_space_approx



!####################################################################

!subroutine cool_to_space_exact (is, ie, js, je, cldtf,          &
subroutine cool_to_space_exact (                cldtf,          &
                                Atmos_input, Optical, Gas_tf,  &
                              sorc,        to3cnt, Lw_diagnostics, &
                                cts_sum, cts_sumcf, &
                                gxctscf) 

!-----------------------------------------------------------------------
!      cool_to_space calculates the cool-to-space cooling rate for 
!      a band n.
!
!     author: m. d. schwarzkopf
!
!     revised: 7/21/94
!
!     certified:  radiation version 2.0
! 
!-----------------------------------------------------------------------
!     intent in:
!
!     cldtf    = cloud transmission function between levels k and 
!                  level KS.
!
!     co2spnb  =  co2 transmission functions between levels k and 
!                  level KS for the freq bands in the 15 um range. 
!                 the first band may include n2o 17um transmissivities.
!
!     pdfinv   =  inverse of pressure difference between flux levels.
!
!     pflux    =  pressure at flux levels of model.
!
!     press    =  pressure at data levels of model.
!
!     sorc     =   band-integrated Planck function, for each combined
!                band in the 160-1200 cm-1 region.
!
!     temp     =  temperature at data levels of model.
!
!     to3cnt   =   transmission functions between levels k and
!                   level KS for the 990-1070 cm-1 range.
!
!     totvo2   =  summed h2o continuum path from top of atmosphere to
!                   flux level k.
!
!     var1     =  h2o optical path in model layers.
!
!     var2     =  pressure-weighted h2o optical path in model layers.
!-----------------------------------------------------------------------

!integer, intent(in) :: is, ie, js, je
real, dimension (:,:,:,:), intent(in)     :: cldtf, sorc
real, dimension (:,:,:),   intent(in)     ::                     to3cnt
real, dimension(:,:,:),      intent(inout)  :: cts_sum, cts_sumcf
real, dimension(:,:),      intent(inout)  ::  gxctscf
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
type(optical_path_type), intent(inout) :: Optical
type(gas_tf_type), intent(inout) :: Gas_tf 
type(atmos_input_type), intent(in) :: Atmos_input

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
      integer        :: n, k, j, ioffset
    real, dimension(size(Atmos_input%pflux,1), size(Atmos_input%pflux,2), &
                    size(Atmos_input%pflux,3)-1) :: pdfinv, pdfinv2
    real, dimension(size(Atmos_input%pflux,1), size(Atmos_input%pflux,2), &
                    size(Atmos_input%pflux,3)  ) :: dte1, &
                      press, temp, pflux
    integer, dimension(size(Atmos_input%pflux,1), size(Atmos_input%pflux,2), &
                    size(Atmos_input%pflux,3)  ) :: ixoe1 

    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2))   ::   &
!                     gxctsbd, pfac1, pfac2, gxctsbdcf
                               pfac1, pfac2
    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2),        &
                    size(Atmos_input%pflux,3)) :: &
                       sorc_tmp,        ctmp, totch2o_tmp, totaer_tmp
    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2),        &
                    2:size(Atmos_input%pflux,3)) :: &
                        totvo2_tmp

    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2),        &
                    size(Atmos_input%pflux,3)-1) :: &
		                  exctscf, tt, x, y, topm, &
                      topphi, phitmp, psitmp, ag, agg, f, ff, &
                      tmp1, tmp2, fac1, fac2, cfc_tf
    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2),        &
                    size(Atmos_input%pflux,3)-1, NBLY) :: &
                            exctsncf

    real, dimension(size(Atmos_input%pflux,1),   &
                    size(Atmos_input%pflux,2),        &
                                                 NBLY) :: &
		                fctsgcf




!  convert press and pflux to cgs.
        press(:,:,:) = 10.0*Atmos_input%press(:,:,:)
       pflux(:,:,:) = 10.0*Atmos_input%pflux(:,:,:)
      temp(:,:,:) = Atmos_input%temp(:,:,:)




      pdfinv2(:,:,KS:KE) = 1.0/(pflux(:,:,KS+1:KE+1) - pflux(:,:,KS:KE))
     pdfinv(:,:,KS:KE) = 1.0/(Atmos_input%pflux(:,:,KS+1:KE+1) -   &
                               Atmos_input%pflux(:,:,KS:KE))
!----------------------------------------------------------------------
      ioffset = Lw_parameters%offset

!-----------------------------------------------------------------------
!     initialize quantities.
!-----------------------------------------------------------------------
      Lw_diagnostics%excts(:,:,KS:KE) = 0.0E+00
      Lw_diagnostics%gxcts(:,:)       = 0.0E+00

      if (Rad_control%do_totcld_forcing) then
        exctscf(:,:,KS:KE) = 0.0E+00
        gxctscf(:,:)       = 0.0E+00
      endif


!-----------------------------------------------------------------------
!     compute temperature quantities.
!-----------------------------------------------------------------------
      x(:,:,KS:KE) = temp(:,:,KS:KE) - 2.5E+02
      y(:,:,KS:KE) = x(:,:,KS:KE)*x(:,:,KS:KE)
      ctmp(:,:,KS) = 1.0E+00



       call locate_in_table(temp_1, temp, dte1, ixoe1, KS, KE+1)




      Lw_diagnostics%fctsg(:,:,NBLY) = 0.0
      do n=1,NBLY-1
!-----------------------------------------------------------------------
!     obtain temperature correction capphi, cappsi, then multiply
!     by optical path var1, var2 to compute temperature-corrected
!     optical path and mean pressure for a layer: phitmp, psitmp.
!-----------------------------------------------------------------------
        f     (:,:,KS:KE) = 0.44194E-01*(apcm (n)*x(:,:,KS:KE) +  &
                            bpcm (n)*y(:,:,KS:KE)) 
        ff    (:,:,KS:KE) = 0.44194E-01*(atpcm(n)*x(:,:,KS:KE) +   &
                            btpcm(n)*y(:,:,KS:KE))
        ag    (:,:,KS:KE) = (1.418191E+00 + f (:,:,KS:KE))*   &
                            f (:,:,KS:KE) + 1.0E+00
        agg   (:,:,KS:KE) = (1.418191E+00 + ff(:,:,KS:KE))*   &
                            ff(:,:,KS:KE) + 1.0E+00 
        phitmp(:,:,KS:KE) = Optical%var1(:,:,KS:KE)*     &
                            ((((ag (:,:,KS:KE)*        &
                            ag (:,:,KS:KE))**2)**2)**2)
        psitmp(:,:,KS:KE) = Optical%var2(:,:,KS:KE)*     &
                            ((((agg(:,:,KS:KE)*       &
                            agg(:,:,KS:KE))**2)**2)**2)

!-----------------------------------------------------------------------
!     obtain optical path and mean pressure from the top of the 
!     atmosphere to the level k.
!-----------------------------------------------------------------------
        topm  (:,:,KS) = phitmp(:,:,KS) 
        topphi(:,:,KS) = psitmp(:,:,KS) 
        do k=KS+1,KE
          topm  (:,:,k) = topm  (:,:,k-1) + phitmp(:,:,k) 
          topphi(:,:,k) = topphi(:,:,k-1) + psitmp(:,:,k) 
        enddo
!-----------------------------------------------------------------------
!     tt is the cloud-free cool-to-space transmission function.
!-----------------------------------------------------------------------

        fac1(:,:,KS:KE    ) = acomb(n)*topm(:,:,KS:KE)
        fac2(:,:,KS:KE    ) = fac1(:,:,KS:KE)*topm(:,:,KS:KE)/   &
                              (bcomb(n)*topphi(:,:,KS:KE))
        tmp1(:,:,KS:KE) = fac1(:,:,KS:KE)/SQRT(1.0E+00 +     &
                          fac2(:,:,KS:KE))

        


!-----------------------------------------------------------------------
        if (n .GE. band_no_start(1) .and. n .le. band_no_end(1)) then
          if (Lw_control%do_ckd2p1) then         !  combined bands 1-24.
	    call get_totch2o (n, Optical, totch2o_tmp, dte1, ixoe1)
            tt(:,:,KS:KE) = EXP(-1.0*(tmp1(:,:,KS:KE) + diffac*   &
                            totch2o_tmp(:,:,KS+1:KE+1)))
          else                                   !  combined bands 1-4.
            tt(:,:,KS:KE) = EXP(-1.0*tmp1(:,:,KS:KE)) 
          endif

!-----------------------------------------------------------------------
	else &
        if (n .GE. band_no_start(2) .and. n .le. band_no_end(2)) then
          if (Lw_control%do_ckd2p1) then    !   combined bands 25-40.
	    call get_totch2o (n, Optical, totch2o_tmp, dte1, ixoe1)
            tt(:,:,KS:KE) = EXP(-1.0*(tmp1(:,:,KS:KE) + diffac*   &
                            totch2o_tmp(:,:,KS+1:KE+1)))
	  else                              !  combined bands 5-8.
	    call get_totvo2 (n, Optical, totvo2_tmp)
            tt(:,:,KS:KE) = EXP(-1.0*(tmp1(:,:,KS:KE) +  &
                                totvo2_tmp(:,:,KS+1:KE+1))) 
          endif

!-----------------------------------------------------------------------
	else &                            !  first co2 band
        if (n .GE. band_no_start(3) .and. n .le. band_no_end(3)) then
          if (Lw_control%do_ckd2p1) then           !  combined band 41 
	    call get_totch2obd (n-40, Optical, totch2o_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*     &
                            totch2o_tmp(:,:,KS+1:KE+1     )
          else                                     !   combined band 9.
	    call get_totvo2 (n, Optical, totvo2_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) +               &
                                totvo2_tmp(:,:,KS+1:KE+1)
          endif
	  if (Lw_control%do_lwaerosol) then
!    call get_totaerooptdep(1, Optical, totaer_tmp)
            totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,1)
	    tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +        &
                            totaer_tmp  (:,:,KS+1:KE+1  )
          endif
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))*    &
                                (Gas_tf%co2spnb(:,:,KS+1:KE+1,1)* &
                                Gas_tf%tn2o17(:,:,KS+1:KE+1))

!-----------------------------------------------------------------------
	else  &                      ! second co2 band
        if (n .GE. band_no_start(4) .and. n .le. band_no_end(4)) then
          if (Lw_control%do_ckd2p1) then         !   combined band 42.
	    call get_totch2obd (n-40, Optical, totch2o_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*     &
                            totch2o_tmp(:,:,KS+1:KE+1     )
          else                                   ! combined band 10. 
	    call get_totvo2 (n, Optical, totvo2_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) +               &
                                totvo2_tmp(:,:,KS+1:KE+1)
          endif
	  if (Lw_control%do_lwaerosol) then
!    call get_totaerooptdep(2, Optical, totaer_tmp)
            totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,2)
	    tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +          &
                            totaer_tmp   (:,:,KS+1:KE+1  )
          endif
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))*     &
                                Gas_tf%co2spnb(:,:,KS+1:KE+1,2)

!-----------------------------------------------------------------------
	else &                                   !  third co2 band
        if (n .GE. band_no_start(5) .and. n .le. band_no_end(5)) then
          if (Lw_control%do_ckd2p1) then      !   combined band 43.
	    call get_totch2obd (n-40, Optical, totch2o_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*      &
                            totch2o_tmp(:,:,KS+1:KE+1     )
          else                                !   combined band 11. 
	    call get_totvo2 (n, Optical, totvo2_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) +                &
                                totvo2_tmp(:,:,KS+1:KE+1)
          endif
	  if (Lw_control%do_lwaerosol) then
!	    call get_totaerooptdep(3, Optical, totaer_tmp)
            totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,3)
	    tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +       &
                            totaer_tmp   (:,:,KS+1:KE+1  )
          endif
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))*     &
                                Gas_tf%co2spnb(:,:,KS+1:KE+1,3)

!-----------------------------------------------------------------------
	else &
        if (n .GE. band_no_start(6) .and. n .le. band_no_end(6)) then
          if (Lw_control%do_ckd2p1) then     ! combined bands 44-45.
	    call get_totch2obd (n-40, Optical, totch2o_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*    &
                            totch2o_tmp(:,:,KS+1:KE+1     )
          else                               ! combined bands 12-13.
	    call get_totvo2 (n, Optical, totvo2_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) +               &
                                totvo2_tmp(:,:,KS+1:KE+1)
          endif
	  if (Lw_control%do_lwaerosol) then
!	    call get_totaerooptdep(n-8-ioffset, Optical, totaer_tmp)
            totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,n-8-ioffset)
	    tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +         &
                            totaer_tmp   (:,:,KS+1:KE+1            )
          endif
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))

!-----------------------------------------------------------------------
	else   &                             !  combined band 14 or 46.
        if (n .GE. band_no_start(7) .and. n .le. band_no_end(7)) then
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp1(:,:,KS:KE))*      &
!                              to3cnt (:,:,KS:KE)
                               to3cnt (:,:,KS+1:KE+1)

!-----------------------------------------------------------------------
	else   &
        if (n .GE. band_no_start(8) .and. n .le. band_no_end(8)) then
          if (Lw_control%do_ckd2p1) then        !  combined band 47.
	    call get_totch2obd (n-40, Optical, totch2o_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*    &
                            totch2o_tmp(:,:,KS+1:KE+1     )
          else                                  !  combined band 15.
	    call get_totvo2 (n, Optical, totvo2_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) +              &
                                totvo2_tmp(:,:,KS+1:KE+1)
          endif
	  if (Lw_control%do_lwaerosol) then
            totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,7)
	    tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE)    +                &
                            totaer_tmp   (:,:,KS+1:KE+1  )
          endif
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))*   &
                                Gas_tf%n2o9c  (:,:,KS+1:KE+1)
        endif
!--------------------------------------------------------------------
!     calculate or retrieve the source function for the current band.
!--------------------------------------------------------------------
        if (n <= 8 + ioffset) then
	  call looktab (tabsr, ixoe1, dte1, sorc_tmp, KS, KE+1, n)
        else
         sorc_tmp(:,:,:) = sorc(:,:,:,n-8-ioffset)
        endif

!---------------------------------------------------------------------
!     retrieve the cfc effect if cfcs are activated.
!---------------------------------------------------------------------
        if (Lw_control%do_cfc .and. n >= 9+ioffset) then
          call cfc_exact(n-8-ioffset, Optical, cfc_tf)
	  tt(:,:,KS:KE) = tt(:,:,KS:KE)*cfc_tf(:,:,KS:KE)
        endif

!------------------------------------------------------------------
!     define some near-surface pressure functions that are needed
!------------------------------------------------------------------
        pfac1(:,:) = 0.5*pdfinv2(:,:,KE)*(pflux(:,:,KE+1) - &
                     press(:,:,KE))*tt(:,:,KE-1)
        pfac2(:,:) = 0.5*pdfinv2(:,:,KE)*(pflux(:,:,KE+1) +    &
	             press(:,:,KE) - 2.0*pflux(:,:,KE))*tt(:,:,KE)

!--------------------------------------------------------------------
!     calculate the ground fluxes (?)
!--------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
!         gxctsbdcf(:,:) = tt(:,:,KE)*sorc_tmp(:,:,KE) +   &
          fctsgcf(:,:,n) = tt(:,:,KE)*sorc_tmp(:,:,KE) +   &
	                   (pfac1(:,:) + pfac2(:,:))*   &
		           (sorc_tmp(:,:,KE+1) - sorc_tmp(:,:,KE))
!         gxctsbd(:,:) = cldtf(:,:,KE+1,cld_indx_table(n+32-ioffset))* &
    Lw_diagnostics%fctsg(:,:,n) = cldtf(:,:,KE+1,cld_indx_table(n+32-ioffset))* &
!                 gxctsbdcf(:,:)
	                 fctsgcf(:,:,n)
!         gxctscf(:,:) = gxctscf(:,:) + gxctsbdcf(:,:)
          gxctscf(:,:) = gxctscf(:,:) + fctsgcf(:,:,n)
!         Lw_diagnostics%gxcts(:,:) = Lw_diagnostics%gxcts(:,:) + gxctsbd(:,:)
          Lw_diagnostics%gxcts(:,:) = Lw_diagnostics%gxcts(:,:) +   &
	              Lw_diagnostics%fctsg(:,:,n)
        else
!         gxctsbd(:,:) = cldtf(:,:,KE+1,cld_indx_table(n+32-ioffset))* &
Lw_diagnostics%fctsg(:,:,n) = cldtf(:,:,KE+1,cld_indx_table(n+32-ioffset))* &
                         (tt(:,:,KE)*sorc_tmp(:,:,KE) +  &
	                 (pfac1(:,:) + pfac2(:,:))*   &
		         (sorc_tmp(:,:,KE+1) - sorc_tmp(:,:,KE)))
!         Lw_diagnostics%gxcts(:,:) = Lw_diagnostics%gxcts(:,:) + gxctsbd(:,:)
          Lw_diagnostics%gxcts(:,:) = Lw_diagnostics%gxcts(:,:) +   &
	  Lw_diagnostics%fctsg(:,:,n)
        endif
!7 convert to mks
        Lw_diagnostics%fctsg(:,:,n) = 1.0e-03*   &
	                              Lw_diagnostics%fctsg(:,:,n)

!! the following array only needed when diagnostics on
!       if (            do_diagnostics) then
!         Lw_diagnostics%fctsg(:,:,n) = gxctsbd(:,:)
!  if (Rad_control%do_totcld_forcing) then
!           fctsgcf(:,:,n) = gxctsbdcf(:,:)
!         endif
!       endif

!--------------------------------------------------------------------
!    include the effect of the cloud transmission function
!--------------------------------------------------------------------
        ctmp(:,:,KS+1:KE+1) = tt(:,:,KS:KE)*        &
                      cldtf(:,:,KS+1:KE+1,cld_indx_table(n+32-ioffset)) 

!---------------------------------------------------------------------
!    if diagnostics is on, save each band's contribution separately.
!    exctsn is the cool-to-space heating rate for each frequency
!    band. fctsg is the "exact" surface flux for each frequency
!    band in the 160-1200 cm-1 range.
!---------------------------------------------------------------------
!! the following array only needed when diagnostics on
!       if (            do_diagnostics) then
          Lw_diagnostics%exctsn(:,:,KS:KE,n) = sorc_tmp(:,:,KS:KE)*     &
                                (ctmp(:,:,KS+1:KE+1) - ctmp(:,:,KS:KE))
          if (Rad_control%do_totcld_forcing) then
	    exctsncf(:,:,KS,n) =  sorc_tmp(:,:,KS)*         &
                                  (tt(:,:,KS) - 1.0E+00)
            exctsncf(:,:,KS+1:KE,n) = sorc_tmp(:,:,KS+1:KE)*     &
                                     (tt(:,:,KS+1:KE) - tt(:,:,KS:KE-1))
          endif
!       endif

!-----------------------------------------------------------------------
!     excts is the cool-to-space cooling rate accumulated over
!     frequency bands.
!-----------------------------------------------------------------------
!       Lw_diagnostics%excts(:,:,KS:KE) = Lw_diagnostics%excts(:,:,KS:KE) + sorc_tmp(:,:,KS:KE)*    &
!                          (ctmp(:,:,KS+1:KE+1) - ctmp(:,:,KS:KE))
        Lw_diagnostics%excts(:,:,KS:KE) =    &
	       Lw_diagnostics%excts(:,:,KS:KE) +   &
	       Lw_diagnostics%exctsn(:,:,KS:KE,n)  
        if (Rad_control%do_totcld_forcing) then
!         exctscf(:,:,KS) = exctscf(:,:,KS) + sorc_tmp(:,:,KS)*     &
!                           (tt(:,:,KS) - 1.0E+00)
!         exctscf(:,:,KS+1:KE) = exctscf(:,:,KS+1:KE) +     &
!                                sorc_tmp(:,:,KS+1:KE)*      &
!                                (tt(:,:,KS+1:KE) - tt(:,:,KS:KE-1))
          exctscf(:,:,KS) = exctscf(:,:,KS) +     &
	                     exctsncf(:,:,KS,n)
          exctscf(:,:,KS+1:KE) = exctscf(:,:,KS+1:KE) +     &
	                     exctsncf(:,:,KS+1:KE,n)
        endif
      end do

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     gxcts is the "exact" surface flux accumulated over
!     the frequency bands in the 160-1200 cm-1 range.
!     obtain cool-to-space flux at the top by integration of heating
!     rates and using cool-to-space flux at the bottom (current value 
!     of gxcts).  note that the pressure quantities and conversion
!     factors have not been included either in excts or in gxcts.
!     these cancel out, thus reducing computations.
!-----------------------------------------------------------------------
      do k=KS,KE
        Lw_diagnostics%gxcts(:,:) = Lw_diagnostics%gxcts(:,:) - Lw_diagnostics%excts(:,:,k)
      enddo
      if (Rad_control%do_totcld_forcing) then
        do k=KS,KE
          gxctscf(:,:) = gxctscf(:,:) - exctscf(:,:,k)
        enddo  
      endif

!-----------------------------------------------------------------------
!     now scale the cooling rate excts by including the pressure 
!     factor pdfinv and the conversion factor radcon.
!-----------------------------------------------------------------------
!! the following array only needed when diagnostics on
!     if (            do_diagnostics) then
        do n=1,NBLY-1
	  Lw_diagnostics%exctsn(:,:,KS:KE,n) = 1.0e-03*(Lw_diagnostics%exctsn(:,:,KS:KE,n)*radcon_mks*     &
                                pdfinv(:,:,KS:KE))
        enddo
        if (Rad_control%do_totcld_forcing) then
          do n=1,NBLY-1
	    exctsncf(:,:,KS:KE,n) = 1.0e-03*(exctsncf(:,:,KS:KE,n)*radcon_mks*    &
                                    pdfinv(:,:,KS:KE) )
          enddo
        endif
!     endif
!6    Lw_diagnostics%excts(:,:,KS:KE) = Lw_diagnostics%excts(:,:,KS:KE)*radcon_mks*pdfinv(:,:,KS:KE)

      Lw_diagnostics%excts(:,:,KS:KE) = (Lw_diagnostics%excts(:,:,KS:KE)*radcon_mks*pdfinv(:,:,KS:KE))*1.0e-03

      if (Rad_control%do_totcld_forcing) then
!6        exctscf(:,:,KS:KE) = exctscf(:,:,KS:KE)*radcon_mks*pdfinv(:,:,KS:KE)
        exctscf(:,:,KS:KE) = (exctscf(:,:,KS:KE)*radcon_mks*pdfinv(:,:,KS:KE))*1.0e-03
      endif



!--------------------------------------------------------------------
!    save the heating rates to be later sent to longwave_driver_mod.
!--------------------------------------------------------------------
      cts_sum(:,:,:) = Lw_diagnostics%excts(:,:,:)
      cts_sumcf(:,:,:) = exctscf(:,:,:)


!--------------------------------------------------------------------


end  subroutine cool_to_space_exact




!####################################################################

subroutine e1e290 (Atmos_input,                 e1ctw1, e1ctw2,   &
                   trans_band1, trans_band2, Optical, tch4n2oe, &
                   t4, Lw_diagnostics, cldtf, cld_indx, flx1e1cf, &
                                                     tcfc8)

!-----------------------------------------------------------------------
real, dimension (:,:,:,:), intent(in)   ::  tch4n2oe                  
real, dimension (:,:,:,:), intent(out)   ::  trans_band1, trans_band2  
real, dimension (:,:,:), intent(out)   ::                  e1ctw1, &
					   e1ctw2
real, dimension(:,:,:), intent(in) :: t4
real, dimension(:,:), intent(out) :: flx1e1cf
real, dimension(:,:,:,:), intent(in) :: cldtf
integer, dimension(:), intent(in) :: cld_indx
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
type(optical_path_type), intent(inout) :: Optical
type(atmos_input_type), intent(in) :: Atmos_input
real, dimension (:,:,:),           intent(inout) ::  tcfc8           
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     E1e290 computes two different quantities.
!     
!     1) emissivities used to compute the exchange terms for flux at the
!     top of the atmosphere (level KS). (the top layer, isothermal by
!     assumption, does not contribute to photon exchanges with other
!     layers). these terms are obtained using precomputed e2 functions
!     (see ref. (2)).
!
!     2) emissivities used to obtain the cool-to-space heating rates
!     for all pressure layers. these are obtained using precomputed
!     e1 functions (see ref. (2)).
!
!     the frequency ranges for the e2 calculations are 0-560 and 1200-
!     2200 cm-1. the CTS calculations also require calculations in the
!     160-560 cm-1 range. (see refs. (1) and (2)).
!ifdef ch4n2o
!
!     if ch4 and n2o are included, the frequency range for emissivities
!     is 1400-2200 cm-1, with separate emissivity quantities for the
!     1200-1400 cm-1 range.
!endif ch4n2o
!
!     the reason for combining these calculations is that both use
!     the same scaled h2o amount (avephi) as input, thus reducing
!     some calculation time for obtaining index quantities.
!   
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1981), 9075-9096.
!
!     (2) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------

      real, dimension (size(trans_band2,1), size(trans_band2,2), &
                         size(trans_band2,3)) :: dte1, dte2

      integer, dimension (size(trans_band2,1), size(trans_band2,2), &
                         size(trans_band2,3)) :: ixoe1, ixoe2
      real, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)) ::   &
                          temp, tflux,        totaer_tmp, taero8, &
                           tmp1, tmp2, &
                         e1cts1, e1cts2, &
                                            avphilog, dt1, du, dup

      integer, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)) ::   &
                                                  ixo1, iyo, iyop 

      real, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2)) :: &
                        s1a, flxge1, flxge1cf
      real, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2), NBTRGE) :: &
                        flx1e1fcf, flxge1f, flxge1fcf
      real, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3), NBTRGE) ::   &
                         e1cts1f, e1cts2f

      integer  :: k, m

     tflux(:,:,:) = Atmos_input%tflux(:,:,:)
     temp(:,:,:) = Atmos_input%temp(:,:,:)

!---------------------------------------------------------------------
!     obtain the "exchange" emissivities as a function of temperature 
!     (fxo) and water amount (fyo). the temperature indices have
!     been obtained in longwave_setup_mod.
!-------------------------------------------------------------------
  call locate_in_table(temp_1, temp, dte1, ixoe1, KS, KE+1)
  call locate_in_table(temp_1, tflux, dte2, ixoe2, KS, KE+1)

  ixoe2(:,:,KS:KE) = ixoe2(:,:,KS+1:KE+1)
  dte2 (:,:,KS:KE) = dte2 (:,:,KS+1:KE+1)
  ixoe2(:,:,KE+1)     = ixoe1(:,:,KE)
  dte2 (:,:,KE+1)     = dte1 (:,:,KE)



      avphilog(:,:,KS:KE+1) = LOG10(Optical%avephi(:,:,KS:KE+1))
      call locate_in_table (mass_1, avphilog, du, iyo, KS, KE+1)
      call looktab (tab2, ixoe2, iyo, dte2, du, trans_band2(:,:,:,1), KS, KE+1)

!-----------------------------------------------------------------------
!     the special case emiss(:,:,KE+1) for layer KE+1 is obtained by 
!     averaging the values for KE and KE+1.
!---------------------------------------------------------------------
      trans_band2(:,:,KE+1,1) = 0.5E+00*(trans_band2(:,:,KE,1) +  &
                          trans_band2(:,:,KE+1, 1))
      trans_band2(:,:,KS+1:KE,1) = trans_band2(:,:,KS:KE-1,1)
 
!---------------------------------------------------------------------
!     perform calculations for the e1 function. the terms involving top
!     layer du are not known.  we use index two to represent index one
!     in previous calculations.
!--------------------------------------------------------------------
      iyop(:,:,KS)        = 1
      iyop(:,:,KS+1:KE+1) = iyo(:,:,KS:KE)
      dup (:,:,KS)        = 0.0E+00
      dup (:,:,KS+1:KE+1) = du (:,:,KS:KE)
      do k=KS,KE+1
        ixo1(:,:,k) = ixoe1(:,:,KS)
        dt1 (:,:,k) = dte1 (:,:,KS)
      enddo

!-----------------------------------------------------------------------
!     e1flx(:,:,KS) equals e1cts1(:,:,KS).
!-----------------------------------------------------------------------
      call looktab (tab1, ixoe1, iyop, dte1, dup, e1cts1, KS, KE+1)
      call looktab (tab1, ixoe1, iyo, dte1, du, e1cts2, KS, KE)
      call looktab (tab1, ixo1, iyop, dt1, dup, trans_band1(:,:,:,1), KS, KE+1)
      call looktab (tab1w, ixoe1, iyop, dte1, dup, e1ctw1, KS, KE+1)
      call looktab (tab1w, ixoe1, iyo, dte1, du, e1ctw2, KS, KE)

!--------------------------------------------------------------------
!     calculations with ch4 and n2o require NBTRGE separate emissivity
!     bands for h2o.
!--------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          avphilog(:,:,KS:KE+1) = LOG10(Optical%avephif(:,:,KS:KE+1,m))
          call locate_in_table (mass_1, avphilog, du, iyo, KS , KE+1)
          iyop(:,:,KS)        = 1
          iyop(:,:,KS+1:KE+1) = iyo(:,:,KS:KE)
          dup (:,:,KS)        = 0.0E+00
          dup (:,:,KS+1:KE+1) = du (:,:,KS:KE)
          call looktab (tab2a, ixoe2, iyo, dte2, du, trans_band2(:,:,:,6+m), &
		        KS, KE+1, m)
          call looktab (tab1a, ixoe1, iyop, dte1, dup,    &
		        e1cts1f(:,:,:,m), KS, KE+1, m)
          call looktab (tab1a, ixoe1, iyo, dte1, du, e1cts2f(:,:,:,m),&
		        KS, KE, m)
          call looktab (tab1a, ixo1, iyop, dt1, dup,   &
                        trans_band1(:,:,:,6+m), KS, KE+1, m)
        enddo

!--------------------------------------------------------------------
!     the special case emissf(:,:,KE+1,m) for layer KE+1 is obtained by 
!     averaging the values for KE and KE+1.
!--------------------------------------------------------------------
        do m=1,NBTRGE
          trans_band2(:,:,KE+1,6+m) = 0.5E+00*    &
                             (trans_band2(:,:,KE,6+m) +  &
			           trans_band2(:,:,KE+1,6+m))
      trans_band2(:,:,KS+1:KE,6+m) = trans_band2(:,:,KS:KE-1,6+m)
        enddo
      endif

!---------------------------------------------------------------------
!    add the effects of other radiative gases on these flux arrays.
!    the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!    thus, tch4e and tn2oe are obtained directly from the transmission
!    functions. 
!----------------------------------------------------------------------
   if (Lw_control%do_ch4_n2o) then
     trans_band1 (:,:,KS+1:KE+1,6+1) = trans_band1(:,:,KS+1:KE+1,6+1)*  &
                                tch4n2oe(:,:,KS+1:KE+1,1)
     e1cts1f(:,:,KS  :KE+1,1) = e1cts1f(:,:,KS:KE+1,1)*  &
                                tch4n2oe(:,:,KS:KE+1,1)
     e1cts2f(:,:,KS  :KE,  1) = e1cts2f(:,:,KS:KE,1)*   &
                                tch4n2oe(:,:,KS+1:KE+1,1)
     trans_band2 (:,:,KS+1:KE+1,  6+1) = trans_band2(:,:,KS+1:KE+1,6+1)*   &
                                tch4n2oe(:,:,KS+1:KE+1,1)
 
!----------------------------------------------------------------------
!    add cfc transmissivities if species which absorb in this fre-
!    quency range ( presently index 8) are present.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
     if (Lw_control%do_cfc) then
       call cfc_indx8 (8, Optical, tcfc8)
       trans_band1 (:,:,KS+1:KE+1,6+1) = trans_band1(:,:,KS+1:KE+1,6+1)*  &
                                  tcfc8(:,:,KS+1:KE+1)
       e1cts1f(:,:,KS  :KE+1,1) = e1cts1f(:,:,KS:KE+1,1)*  &
                                  tcfc8(:,:,KS:KE+1)
       e1cts2f(:,:,KS  :KE,  1) = e1cts2f(:,:,KS:KE,1)*  &
                                  tcfc8(:,:,KS+1:KE+1)
       trans_band2 (:,:,KS+1:KE+1,6+1) = trans_band2(:,:,KS+1:KE+1,6+1)*   &
                                  tcfc8(:,:,KS+1:KE+1)
     endif 

!----------------------------------------------------------------------
!    compute aerosol transmission function for 1200-1400 cm-1 region
!----------------------------------------------------------------------
     if (Lw_control%do_lwaerosol) then
        totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,8)
       taero8(:,:,KS:KE+1) = EXP(-1.0E+00*totaer_tmp(:,:,KS:KE+1))
       trans_band1(:,:,KS+1:KE+1,6+1) = trans_band1(:,:,KS+1:KE+1,6+1)*   &
                                  taero8(:,:,KS+1:KE+1)
       e1cts1f(:,:,KS  :KE+1,1) = e1cts1f(:,:,KS:KE+1,1)*  &
                                  taero8(:,:,KS:KE+1)
       e1cts2f(:,:,KS  :KE,  1) = e1cts2f(:,:,KS:KE,1)*   &
                                  taero8(:,:,KS+1:KE+1)
       trans_band2 (:,:,KS+1:KE+1,6+1) = trans_band2(:,:,KS+1:KE+1,6+1)*   &
                                  taero8(:,:,KS+1:KE+1)
     endif
   endif

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
!   t4(:,:,:) = source_band(:,:,:,1)
    s1a   (:,:) = t4(:,:,KE+1)*(e1cts1(:,:,KE+1) - e1ctw1(:,:,KE+1))
    flxge1(:,:) = s1a(:,:)*cldtf(:,:,KE+1,1)
    tmp1(:,:,KS:KE) = t4(:,:,KS:KE)*    &
                      (e1cts1(:,:,KS:KE) - e1ctw1(:,:,KS:KE)) 
    tmp2(:,:,KS:KE) = t4(:,:,KS:KE)*   &
                       (e1cts2(:,:,KS:KE) - e1ctw2(:,:,KS:KE))
    Lw_diagnostics%flx1e1(:,:) = flxge1(:,:)
    do k=KS,KE
      Lw_diagnostics%flx1e1(:,:) = Lw_diagnostics%flx1e1(:,:) + tmp1(:,:,k)*cldtf(:,:,k,1) -   &
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
        Lw_diagnostics%flx1e1f(:,:,m) = flxge1f(:,:,m)
        do k=KS,KE
          tmp1(:,:,k) = t4(:,:,k)*e1cts1f(:,:,k,m)
          tmp2(:,:,k) = t4(:,:,k)*e1cts2f(:,:,k,m)
          Lw_diagnostics%flx1e1f(:,:,m) = Lw_diagnostics%flx1e1f(:,:,m) + tmp1(:,:,k)*   &
                           cldtf(:,:,k,cld_indx(7)) - tmp2(:,:,k)*  &
                           cldtf(:,:,k+1,cld_indx(7))
        end do
      end do
      do m=1,NBTRGE
        Lw_diagnostics%flx1e1(:,:) = Lw_diagnostics%flx1e1(:,:) + Lw_diagnostics%flx1e1f(:,:,m)
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

end  subroutine e1e290




!###################################################################

subroutine e290 (Atmos_input, k, trans_band2, trans_band1, Optical,  tch4n2oe, &
                                  tcfc8)

!-----------------------------------------------------------------------
integer,                  intent(in)            ::  k
real, dimension(:,:,:,:),   intent(in)           :: tch4n2oe       
real, dimension(:,:,:,:),   intent(inout)           :: trans_band2, &
                                                        trans_band1
!real, dimension(:,:,:),   intent(out)           :: emiss, emissb  
!real, dimension(:,:,:),   intent(out)           ::        emissb  
real, dimension(:,:,:),   intent(inout)           :: tcfc8          
type(optical_path_type), intent(inout) ::  Optical
type(atmos_input_type), intent(in) ::  Atmos_input
!real, dimension(:,:,:,:), intent(out) :: emissf, emissbf
!real, dimension(:,:,:,:), intent(out) ::         emissbf
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     e290 computes the exchange terms in the flux equation for longwave
!     radiation for all terms except the exchange with the top of the
!     atmosphere.  the method is a table lookup on a pre-computed e2
!     function (defined in reference (2)).  calculation are done in the
!     frequency range: 0-560, 1200-2200 cm-1 for q(approximate).
!     motivation for these calculations is in references (1) and (2).
!
!     references:
!
!     (1) schwarzkopf, m. d., and s. b. fels, "the simplified
!         exchange method revisited: an accurate, rapid method for
!         computation of infrared cooling rates and fluxes," journal
!         of geophysical research, 96 (1981), 9075-9096.
!
!     (2) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: c. h. goldberg
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!
!-----------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables
!-------------------------------------------------------------------
      integer      :: kp, m

      real, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)) :: temp, tflux, &
                                  totaer_tmp, taero8, taero8kp, &
                                                  avphilog, dtk, du

      integer, dimension (size(Atmos_input%temp,1),    &
                       size(Atmos_input%temp,2), &
                       size(Atmos_input%temp,3)) ::              &
                                                  ixok, iyo          

      real, dimension (size(trans_band2,1), size(trans_band2,2), &
                         size(trans_band2,3)) :: dte1, dte2

      integer, dimension (size(trans_band2,1), size(trans_band2,2), &
                         size(trans_band2,3)) :: ixoe1, ixoe2

      temp = Atmos_input%temp
      tflux = Atmos_input%tflux

!-----------------------------------------------------------------------
!     obtain the "exchange" emissivities as a function of temperature 
!     (fxo) and water amount (avephi). the temperature indices have
!     been obtained in Lwrad. calculations are for flux level k, with
!     kp = k+1 to KE+1. the case k=KS is excluded (done in E1e290).
!     calculations are also made for flux levels k to KE, for
!     contributions from flux level k. in this case, the temperature
!     index (ixok) represents tflux(:,:,k-1); the water index (iyo)
!     has the same values as in the case with varying kp.
!---------------------------------------------------------------------

!----------------------------------------------------------------------

  call locate_in_table(temp_1, temp, dte1, ixoe1, KS, KE+1)
  call locate_in_table(temp_1, tflux, dte2, ixoe2, KS, KE+1)

  ixoe2(:,:,KS:KE) = ixoe2(:,:,KS+1:KE+1)
  dte2 (:,:,KS:KE) = dte2 (:,:,KS+1:KE+1)
  ixoe2(:,:,KE+1)     = ixoe1(:,:,KE)
  dte2 (:,:,KE+1)     = dte1 (:,:,KE)


      do kp=k,KE
        ixok(:,:,kp) = ixoe2(:,:,k-1)
        dtk (:,:,kp) = dte2 (:,:,k-1)
      end do

      avphilog(:,:,k:KE+1) = LOG10(Optical%avephi(:,:,k:KE+1))
      call locate_in_table (mass_1, avphilog, du, iyo,k, KE+1)
!     call looktab (tab2, ixoe2, iyo, dte2, du, emiss, k, KE+1)
      call looktab (tab2, ixoe2, iyo, dte2, du, trans_band2(:,:,:,1), k, KE+1)
!     call looktab (tab2, ixok, iyo, dtk, du, emissb, k, KE)
      call looktab (tab2, ixok, iyo, dtk, du, trans_band1(:,:,:,1), k, KE)
       trans_band1(:,:,k+1:KE+1,1) = trans_band1(:,:,k:KE,1)

!--------------------------------------------------------------------
!     the special case emiss(:,:,KE) for layer KE is obtained by 
!     averaging the values for KE and KE+1. note that emiss(:,:,KE+1) 
!     is not useful after this point.
!-------------------------------------------------------------------
!      emiss(:,:,KE) = 0.5E+00*(emiss(:,:,KE) +emiss(:,:,KE+1))
!     emiss(:,:,KE+1) = 0.5E+00*(emiss(:,:,KE) +emiss(:,:,KE+1))
!     emiss(:,:,k+1:KE) = emiss(:,:,k:KE-1)
      trans_band2(:,:,KE+1,1) = 0.5E+00*(trans_band2(:,:,KE,1) +  &
                         trans_band2(:,:,KE+1,1))
      trans_band2(:,:,k+1:KE,1) = trans_band2(:,:,k:KE-1,1)
 
!--------------------------------------------------------------------
!     calculations with ch4 and n2o require NBTRGE separate emissivity
!     bands for h2o. reqults are in emissf (flux level k) and
!     emissbf (other levels).
!-------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          avphilog(:,:,k:KE+1) = LOG10(Optical%avephif(:,:,k:KE+1,m))
          call locate_in_table (mass_1, avphilog, du, iyo, k, KE+1)
!         call looktab (tab2a, ixoe2, iyo, dte2, du, emissf(:,:,:,m), &
          call looktab (tab2a, ixoe2, iyo, dte2, du, trans_band2(:,:,:,6+m), &
			k, KE+1, m)
!          call looktab (tab2a, ixok, iyo, dtk, du, emissbf(:,:,:,m), &
          call looktab (tab2a, ixok, iyo, dtk, du, trans_band1(:,:,:,6+m), &
			k, KE, m)
       trans_band1(:,:,k+1:KE+1,6+m) = trans_band1(:,:,k:KE,6+m)
        enddo

!----------------------------------------------------------------------
!     the special case emissf(:,:,KE) for layer KE is obtained by 
!     averaging the values for KE and KE+1. note that emissf(:,:,KE+1,m)
!     is not useful after this point.
!----------------------------------------------------------------------
        do m=1,NBTRGE
!         emissf(:,:,KE,m) = 0.5E+00*  &
!                            (emissf(:,:,KE,m) +emissf(:,:,KE+1,m))
!         emissf(:,:,KE+1,m) = 0.5E+00*  &
!                            (emissf(:,:,KE,m) +emissf(:,:,KE+1,m))
!         emissf(:,:,k+1:KE,m) = emissf(:,:,k:KE-1,m)
          trans_band2(:,:,KE+1,6+m) = 0.5E+00*  &
                             (trans_band2(:,:,KE,6+m) +   &
			       trans_band2(:,:,KE+1,6+m))
          trans_band2(:,:,k+1:KE,6+m) = trans_band2(:,:,k:KE-1,6+m)
        enddo
      endif

!----------------------------------------------------------------------
!    add the effects of other radiative gases on these flux arrays.
!    the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!    thus, tch4e and tn2oe are obtained directly from the transmission
!    functions. 
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
!       emissbf(:,:,k:KE,1) = emissbf(:,:,k:KE,1)*  &
        trans_band1(:,:,k+1:KE+1,6+1) = trans_band1(:,:,k+1:KE+1,6+1)*  &
   		      tch4n2oe(:,:,k+1:KE+1,1)
!       emissf(:,:,k:KE,1) = emissf(:,:,k:KE,1)*tch4n2oe(:,:,k+1:KE+1,1)
!        emissf(:,:,k+1:KE+1,1) = emissf(:,:,k+1:KE+1,1)*tch4n2oe(:,:,k+1:KE+1,1)
        trans_band2(:,:,k+1:KE+1,6+1) = trans_band2(:,:,k+1:KE+1,6+1)*tch4n2oe(:,:,k+1:KE+1,1)

!--------------------------------------------------------------------
!    add cfc transmissivities if species which absorb in this fre-
!    quency range are present.
!--------------------------------------------------------------------
        if (Lw_control%do_cfc) then
          call cfc_indx8_part (8, Optical, tcfc8, k)
!         emissbf(:,:,k:KE,1) = emissbf(:,:,k:KE,1)*tcfc8(:,:,k+1:KE+1)
          trans_band1(:,:,k+1:KE+1,6+1) = trans_band1(:,:,k+1:KE+1,6+1)*tcfc8(:,:,k+1:KE+1)
!         emissf(:,:,k:KE,1) = emissf(:,:,k:KE,1)*tcfc8(:,:,k+1:KE+1)
!         emissf(:,:,k+1:KE+1,1) = emissf(:,:,k+1:KE+1,1)*tcfc8(:,:,k+1:KE+1)
          trans_band2(:,:,k+1:KE+1,6+1) = trans_band2(:,:,k+1:KE+1,6+1)*tcfc8(:,:,k+1:KE+1)
        endif

!--------------------------------------------------------------------
!     compute aerosol transmission function for 1200-1400 cm-1 region
!    (as quotient of 2 exponentials)
!     taero8kp(k) contains the (k+1,k) transmissivities for all k
!     in the 1200-1400 cm-1 frequency range.
!---------------------------------------------------------------------
        if (Lw_control%do_lwaerosol) then
        totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,8)
       taero8(:,:,KS:KE+1) = EXP(-1.0E+00*totaer_tmp(:,:,KS:KE+1))
          do kp = k+1,KE+1
            taero8kp(:,:,kp) = taero8(:,:,kp)/taero8(:,:,k)
          enddo
!         emissbf(:,:,k:KE,1) = emissbf(:,:,k:KE,1)*  &
          trans_band1(:,:,k+1:KE+1,6+1) = trans_band1(:,:,k+1:KE+1,6+1)*  &
                  taero8kp(:,:,k+1:KE+1)
!         emissf(:,:,k:KE,1) = emissf(:,:,k:KE,1)*  &
!         emissf(:,:,k+1:KE+1,1) = emissf(:,:,k+1:KE+1,1)*  &
          trans_band2(:,:,k+1:KE+1,6+1) = trans_band2(:,:,k+1:KE+1,6+1)*  &
                  taero8kp(:,:,k+1:KE+1)
        endif
      endif


end subroutine e290



!####################################################################

subroutine esfc  (Atmos_input,         emspec,             Optical, &
                           emspecf, tch4n2oe, tcfc8 ) 
   
!--------------------------------------------------------------------
type(atmos_input_type), intent(in) :: Atmos_input
real, dimension (:,:,:,:),   intent(in)           ::   tch4n2oe
real, dimension (:,:,:),   intent(inout)           ::   tcfc8   
real, dimension (:,:,:),   intent(out)           ::    emspec
type(optical_path_type), intent(inout) :: Optical
real, dimension (:,:,:,:), intent(out)           ::    emspecf
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------

      integer   :: m
      integer :: k



      real, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                         size(Atmos_input%temp,3)) :: temp, tflux, &
                                            tpl1, tpl2, &
!		              dte1, dte2, tcfc8, &
		              dte1, dte2,        &
                  dxsp, ylog, dysp, emiss, emd1, emd2

      integer, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                         size(Atmos_input%temp,3)) :: ixsp, iysp, & 
			      ixoe1, ixoe2

      real, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                         size(Atmos_input%temp,3), NBTRGE) ::    &

                            emissf, emd2f,  emd1f
!     israd = 1
!     jsrad = 1

!     ks    = 1


!     ierad = size (Atmos_input%temp, 1)
!     jerad = size (Atmos_input%temp, 2)

!     ke    = size (Atmos_input%temp, 3) - 1

     tflux(:,:,:) = Atmos_input%tflux(:,:,:)
      temp(:,:,:) = Atmos_input%temp(:,:,:)
 


      tpl1(:,:,KS)         = temp(:,:,KE)
     tpl1(:,:,KS+1:KE) = tflux(:,:,KS+1:KE)
     tpl1(:,:,KE+1)       = 0.5E+00*(tflux(:,:,KE+1) +   &
                                     temp(:,:,KE))
  tpl2(:,:,KS+1:KE) = tflux(:,:,KS+1:KE)
    tpl2(:,:,KE+1)       = 0.5E+00*(tflux(:,:,KE) +    &
                                temp(:,:,KE))


!--------------------------------------------------------------------

  call locate_in_table(temp_1, temp, dte1, ixoe1, KS, KE+1)
  call locate_in_table(temp_1, tflux, dte2, ixoe2, KS, KE+1)

  ixoe2(:,:,KS:KE) = ixoe2(:,:,KS+1:KE+1)
  dte2 (:,:,KS:KE) = dte2 (:,:,KS+1:KE+1)
  ixoe2(:,:,KE+1)     = ixoe1(:,:,KE)
  dte2 (:,:,KE+1)     = dte1 (:,:,KE)




      ixsp(:,:,KE)   = ixoe2(:,:,KE-1)
      ixsp(:,:,KE+1) = ixoe1(:,:,KE-1)
      dxsp(:,:,KE)   = dte2(:,:,KE-1)
      dxsp(:,:,KE+1) = dte1(:,:,KE-1)

      ylog(:,:,KE  ) = ALOG10(Optical%var2(:,:,KE))
      ylog(:,:,KE+1) = ALOG10(Optical%var2(:,:,KE) + Optical%empl1(:,:,KE))

      call locate_in_table (mass_1, ylog, dysp, iysp, KE, KE+1)

!--------------------------------------------------------------------
!     compute exchange terms in the flux equation for two terms used
!     for nearby layer computations.
!--------------------------------------------------------------------
      call looktab (tab2, ixsp, iysp, dxsp, dysp, emiss, KE, KE+1)

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          ylog(:,:,KE  ) = ALOG10(Optical%vrpfh2o(:,:,KE,m))
          ylog(:,:,KE+1) = ALOG10(Optical%vrpfh2o(:,:,KE,m) + Optical%empl1f(:,:,KE,m))

          call locate_in_table (mass_1, ylog, dysp, iysp, KE, KE+1)

!-----------------------------------------------------------------------
!     compute exchange terms in the flux equation for two terms used
!     for nearby layer computations.
!---------------------------------------------------------------------
          call looktab (tab2a, ixsp, iysp, dxsp, dysp, emissf(:,:,:,m),&
		        KE, KE+1, m)
        enddo
      endif
!-----------------------------------------------------------------------
!     compute nearby layer transmissivities for h2o.
!--------------------------------------------------------------------
      call locate_in_table (temp_1, tpl1, dxsp, ixsp, KS, KE+1)
      ylog(:,:,KS:KE+1) = ALOG10(Optical%empl1(:,:,KS:KE+1))
      call locate_in_table (mass_1, ylog, dysp, iysp, KS, KE+1)
      call looktab (tab3, ixsp, iysp, dxsp, dysp, emd1, KS, KE+1)

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          ylog(:,:,KS:KE+1) = ALOG10(Optical%empl1f(:,:,KS:KE+1,m))
          call locate_in_table (mass_1, ylog, dysp, iysp, KS, KE+1)
          call looktab (tab3a, ixsp, iysp, dxsp, dysp, emd1f(:,:,:,m), &
			KS, KE+1, m)
        enddo
      endif

!---------------------------------------------------------------------
      call locate_in_table (temp_1, tpl2, dxsp, ixsp, KS+1, KE+1)
      ylog(:,:,KS+1:KE+1) = ALOG10(Optical%empl2(:,:,KS+1:KE+1))
      call locate_in_table (mass_1, ylog, dysp, iysp, KS+1, KE+1)
      call looktab (tab3, ixsp, iysp, dxsp, dysp, emd2, KS+1, KE+1)

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          ylog(:,:,KS+1:KE+1) = ALOG10(Optical%empl2f(:,:,KS+1:KE+1,m))
          call locate_in_table (mass_1, ylog, dysp, iysp, KS+1, KE+1)
          call looktab (tab3a, ixsp, iysp, dxsp, dysp, emd2f(:,:,:,m), &
			KS+1, KE+1, m)
        enddo
      endif

!---------------------------------------------------------------------- 
!     compute nearby layer and special-case transmissivities for
!     emissivity using methods for h2o given in reference (4).
!-------------------------------------------------------------------- 
      emspec(:,:,KS     ) = (emd1(:,:,KS)*Optical%empl1(:,:,KS) -    &
                             emd1(:,:,KE+1)*Optical%empl1(:,:,KE+1))/  &
                             Optical%emx1(:,:) + 0.25E+00*(emiss(:,:,KE) +   &
                             emiss(:,:,KE+1))
      emspec(:,:,KS+1) = 2.0E+00*(emd1(:,:,KS)*Optical%empl1(:,:,KS) -    &
                         emd2(:,:,KE+1)*Optical%empl2(:,:,KE+1))/  &
                           Optical%emx2(:,:)

     if (Lw_control%do_ch4_n2o) then
       do m=1,NBTRGE
         emspecf(:,:,KS,m   ) = (emd1f(:,:,KS,m)*Optical%empl1f(:,:,KS,m) -   &
                              emd1f(:,:,KE+1,m)*Optical%empl1f(:,:,KE+1,m))/   &
                   Optical%emx1f(:,:,m) + 0.25E+00*(emissf(:,:,KE,m) +  &
                             emissf(:,:,KE+1,m))
         emspecf(:,:,KS+1,m) = 2.0E+00*    &
                                 (emd1f(:,:,KS,m)*Optical%empl1f(:,:,KS,m) -  &
                              emd2f(:,:,KE+1,m)*Optical%empl2f(:,:,KE+1,m)) / &
                              Optical%emx2f(:,:,m)

       enddo
     endif
!--------------------------------------------------------------------

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

!--------------------------------------------------------------------
!     add cfc transmissivities if species which absorb in this fre-
!    quency range are present.
!----------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_indx8_part (8, Optical, tcfc8, KE)
	do k=KS,KS+1
	  emspecf(:,:,K,1) = emspecf(:,:,K,1)*tcfc8(:,:,KE+1)
	end do
      endif
    endif

!------------------------------------------------------------------



end subroutine esfc 


!######################################################################

subroutine enear (Atmos_input, emisdg,                     Optical, &
                  emisdgf, tch4n2oe, tcfc8          ) 
   
!--------------------------------------------------------------------
type(atmos_input_type), intent(in) :: Atmos_input
real, dimension (:,:,:,:),   intent(in)           ::   tch4n2oe
real, dimension (:,:,:),   intent(inout)           ::   tcfc8       
real, dimension (:,:,:),   intent(out)           ::   emisdg 
type(optical_path_type), intent(inout) :: Optical
real, dimension (:,:,:,:), intent(out)           ::   emisdgf
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------

      integer   :: m
      integer   :: k
!      integer    :: israd, ierad, jsrad, jerad, ks, ke
!     integer    ::                             ks, ke
      real, dimension (size(emisdg,1), size(emisdg,2), &
                         size(emisdg,3)) :: dte1, dte2

      integer, dimension (size(emisdg,1), size(emisdg,2), &
                         size(emisdg,3)) :: ixoe1, ixoe2

      real, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                         size(Atmos_input%temp,3)) :: temp, tflux, &
                                            tpl1, tpl2, &
                  dxsp, ylog, dysp, emiss, emd1, emd2

      integer, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                         size(Atmos_input%temp,3)) :: ixsp, iysp    

      real, dimension (size(Atmos_input%temp,1),   &
                       size(Atmos_input%temp,2), &
                         size(Atmos_input%temp,3), NBTRGE) ::    &

                            emissf, emd2f,  emd1f
!     israd = 1
!     jsrad = 1

!     !ks    = 1

!     ierad = size (emisdg, 1)
!     jerad = size (emisdg, 2)

!      ke    = size (emisdg, 3) - 1

     tflux(:,:,:) = Atmos_input%tflux(:,:,:)
      temp(:,:,:) = Atmos_input%temp(:,:,:)
 


      tpl1(:,:,KS)         = temp(:,:,KE)
     tpl1(:,:,KS+1:KE) = tflux(:,:,KS+1:KE)
     tpl1(:,:,KE+1)       = 0.5E+00*(tflux(:,:,KE+1) +   &
                                     temp(:,:,KE))
  tpl2(:,:,KS+1:KE) = tflux(:,:,KS+1:KE)
    tpl2(:,:,KE+1)       = 0.5E+00*(tflux(:,:,KE) +    &
                                temp(:,:,KE))


!--------------------------------------------------------------------

  call locate_in_table(temp_1, temp, dte1, ixoe1, KS, KE+1)
  call locate_in_table(temp_1, tflux, dte2, ixoe2, KS, KE+1)

  ixoe2(:,:,KS:KE) = ixoe2(:,:,KS+1:KE+1)
  dte2 (:,:,KS:KE) = dte2 (:,:,KS+1:KE+1)
  ixoe2(:,:,KE+1)     = ixoe1(:,:,KE)
  dte2 (:,:,KE+1)     = dte1 (:,:,KE)




      ixsp(:,:,KE)   = ixoe2(:,:,KE-1)
      ixsp(:,:,KE+1) = ixoe1(:,:,KE-1)
      dxsp(:,:,KE)   = dte2(:,:,KE-1)
      dxsp(:,:,KE+1) = dte1(:,:,KE-1)

!--------------------------------------------------------------------
!     compute exchange terms in the flux equation for two terms used
!     for nearby layer computations.
!--------------------------------------------------------------------

!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     compute nearby layer transmissivities for h2o.
!--------------------------------------------------------------------
      call locate_in_table (temp_1, tpl1, dxsp, ixsp, KS, KE+1)
      ylog(:,:,KS:KE+1) = ALOG10(Optical%empl1(:,:,KS:KE+1))
      call locate_in_table (mass_1, ylog, dysp, iysp, KS, KE+1)
      call looktab (tab3, ixsp, iysp, dxsp, dysp, emd1, KS, KE+1)

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          ylog(:,:,KS:KE+1) = ALOG10(Optical%empl1f(:,:,KS:KE+1,m))
          call locate_in_table (mass_1, ylog, dysp, iysp, KS, KE+1)
          call looktab (tab3a, ixsp, iysp, dxsp, dysp, emd1f(:,:,:,m), &
			KS, KE+1, m)
        enddo
      endif

!---------------------------------------------------------------------
      call locate_in_table (temp_1, tpl2, dxsp, ixsp, KS+1, KE+1)
      ylog(:,:,KS+1:KE+1) = ALOG10(Optical%empl2(:,:,KS+1:KE+1))
      call locate_in_table (mass_1, ylog, dysp, iysp, KS+1, KE+1)
      call looktab (tab3, ixsp, iysp, dxsp, dysp, emd2, KS+1, KE+1)

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          ylog(:,:,KS+1:KE+1) = ALOG10(Optical%empl2f(:,:,KS+1:KE+1,m))
          call locate_in_table (mass_1, ylog, dysp, iysp, KS+1, KE+1)
          call looktab (tab3a, ixsp, iysp, dxsp, dysp, emd2f(:,:,:,m), &
			KS+1, KE+1, m)
        enddo
      endif

!---------------------------------------------------------------------- 
!     compute nearby layer and special-case transmissivities for
!     emissivity using methods for h2o given in reference (4).
!-------------------------------------------------------------------- 
      emisdg(:,:,KS+1:KE) = emd2(:,:,KS+1:KE) + emd1(:,:,KS+1:KE)
      emisdg(:,:,KE+1) = 2.0E+00*emd1(:,:,KE+1)

     if (Lw_control%do_ch4_n2o) then
       do m=1,NBTRGE
         emisdgf(:,:,KS+1:KE,m) =    &
                            emd2f(:,:,KS+1:KE,m) + emd1f(:,:,KS+1:KE,m)
         emisdgf(:,:,KE+1,m) = 2.0E+00*emd1f(:,:,KE+1,m)
       enddo
     endif
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    add the effects of other radiative gases on these flux arrays.
!    the lbl transmissivities imply (at present) NBTRG = NBTRGE = 1).
!    thus, tch4e and tn2oe are obtained directly from the transmission
!    functions. 
!----------------------------------------------------------------------
    if (Lw_control%do_ch4_n2o) then
      emisdgf(:,:,KS+1:KE+1,1) = emisdgf(:,:,KS+1:KE+1,1) *   &
                                 tch4n2oe(:,:,KS+1:KE+1,1)

!--------------------------------------------------------------------
!     add cfc transmissivities if species which absorb in this fre-
!    quency range are present.
!----------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_indx8_part (8, Optical, tcfc8, KE)
	emisdgf(:,:,KS+1:KE+1,1) = emisdgf(:,:,KS+1:KE+1,1) *   &
                                   tcfc8(:,:,KS+1:KE+1)
      endif
    endif

!------------------------------------------------------------------



end subroutine enear



!####################################################################

subroutine co2_source_calc (Atmos_input, Rad_gases, sorc,  Gas_tf, &
                            source_band, dsrcdp_band)

!--------------------------------------------------------------------
real, dimension(:,:,:,:), intent(out) ::  sorc                        
real, dimension(:,:,:,:), intent(out) ::  source_band, dsrcdp_band      
type(gas_tf_type), intent(in) :: Gas_tf
type(atmos_input_type), intent(in) :: Atmos_input
type(radiative_gases_type), intent(in) :: Rad_gases
!----------------------------------------------------------------------

      integer            ::   n, ioffset, m
      integer            ::   NBLY                    
      real :: rrvco2
      real, dimension (size(Atmos_input%temp,1), size(Atmos_input%temp,2), &
                           size(Atmos_input%temp,3)) :: dte1, &
                         press, pflux, temp
      integer, dimension (size(Atmos_input%temp,1), size(Atmos_input%temp,2), &
                           size(Atmos_input%temp,3)) :: ixoe1

!--------------------------------------------------------------------
      ioffset = Lw_parameters%offset
      NBLY = 16+ioffset


 !  convert press and pflux to cgs.
        press(:,:,:) = 10.0*Atmos_input%press(:,:,:)
        pflux(:,:,:) = 10.0*Atmos_input%pflux(:,:,:)
     temp(:,:,:) = Atmos_input%temp(:,:,:)
      rrvco2 = Rad_gases%rrvco2


       
!----------------------------------------------------------------------
!     compute source function for frequency bands (9+ioffset to NBLY-1) 
!     at layer temperatures using table lookup.
!----------------------------------------------------------------------
      call locate_in_table(temp_1, temp, dte1, ixoe1, KS, Ke+1)
      do n=9+ioffset,NBLY-1
        call looktab (tabsr, ixoe1, dte1,   & 
                      sorc(:,:,:,n-ioffset-8), KS, ke+1, n)
      enddo

!-----------------------------------------------------------------------
!     compute the nlte source function for co2.
!-----------------------------------------------------------------------
      if (do_nlte) then
        call nlte (pflux, press, rrvco2, sorc, Gas_tf)
      endif

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
      source_band(:,:,:,1) =  Atmos_input%temp (:,:,KS:KE+1)**4
      source_band(:,:,:,2) = sorc(:,:,:, 1) + &
                             sorc(:,:,:, 2) + &
                             sorc(:,:,:, 3 )
      source_band(:,:,:,3)  = sorc(:,:,:,4 )
      source_band(:,:,:,4)=  sorc(:,:,:,5 )
      source_band(:,:,:,5) =  sorc(:,:,:, 6 )
      source_band(:,:,:,6) =  sorc(:,:,:,7 )
      do m=1,NBTRGE
      source_band(:,:,:,6+m) =source_band(:,:,:,1)
      end do

      do n=1, 6+NBTRGE       
       dsrcdp_band(:,:,KS+1:KE+1,n) =  source_band(:,:,KS+1:KE+1,n) - &
                               source_band(:,:,KS:KE,n)
      end do

!-------------------------------------------------------------------


end subroutine co2_source_calc




!#####################################################################

subroutine nlte (pflux, press, rrvco2, sorc, Gas_tf)

!-----------------------------------------------------------------------
!     nlte is the present formulation of an nlte calculation of the 
!     source function in the 15 um region (two bands).
!
!     the essential theory is:
!
!           phi = C*j
!             j = b + E*phi
!
!     where
!             C = Curtis matrix
!	      E = NLTE contribution (diagonal matrix)
!           phi = heating rate vector
!             b = LTE source function vector
!             j = NLTE source function vector
!
!             j = b (by assumption) for pressure layers > ixnltr
!             j = b (by assumption) for pressure layers > ixprnlte
!      E is obtained using a formulation devised by Fels (denoted
!      Ri in his notes).
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------
 real, dimension (:,:,:), intent(in)  ::  press, pflux               
 real, dimension (:,:,:,:), intent(inout)  ::  sorc               
 real,                    intent(in)  ::  rrvco2
type(gas_tf_type), intent(in) :: Gas_tf

!-----------------------------------------------------------------------
!     intent local:
!
!     degeneracy factor = 0.5
!
!                fnlte  = NLTE contribution: (E in above notes)
!
!                phifx  = fixed portion of PHI (contributions from
!                         layers > ixnltr, where j(k) = b(k))
!                         layers > ixprnlte, where j(k) = b(k))
!
!                phivar = varying portion of PHI (contributions
!                         from layers <= ixprnlte).
!                         from layers <= ixnltr).
!-----------------------------------------------------------------------
      real                                   :: degen = 0.5
      integer                                :: n, k, inb, kp, ioffset
      real, dimension (size(press,1), size(press,2), &
                       ixprnlte) ::  &
		                ag, az, bdenom, cdiag, &
	                                        tcoll, phifx, phivar
      real, dimension (size(press,1), size(press,2), &
                       ixprnlte, NBCO215) ::  &
                                    fnlte
      real, dimension (size(press,1), size(press,2), &
                       size(press,3)-1, ixprnlte ) ::  &
		                     cmtrx

!---------------------------------------------------------------------
      ioffset =  Lw_parameters%offset

!--------------------------------------------------------------------

!-----------------------------------------------------------------------
!     compute curtis matrix for both frequency bands.
!-----------------------------------------------------------------------
      call co2curt (pflux, cmtrx, Gas_tf)

      do k=KS,ixprnlte
        cdiag(:,:,k) = cmtrx(:,:,k,k)
      end do

!-----------------------------------------------------------------------
!   collisional relaxation time (see fels notes for "tcoll")
!-----------------------------------------------------------------------
      do k=KS,ixprnlte
        tcoll(:,:,k) = degen*1.5E-05*press(:,:,KE+1)/   &
                       (seconds_per_day*press(:,:,k)) 
      end do

!-----------------------------------------------------------------------
!   compute NLTE contribution for eack band at each pressure level
!   <= ixprnlte. fnlte = zero by assumption at other levels.
!-----------------------------------------------------------------------
      do n=1,NBCO215
        fnlte (:,:,KS:ixprnlte,n) = 3.5E+00*tcoll(:,:,KS:ixprnlte)*  &
				    c1b7(n)/(rrvco2*c2b7(n)) 
      enddo

!-----------------------------------------------------------------------
!     begin computations for (NBCO215) bands in 15um range.
!-----------------------------------------------------------------------
      do inb = 1,NBCO215
        bdenom(:,:,KS:ixprnlte) = 1.0E+00/   &
              (1.0E+00 - fnlte(:,:,KS:ixprnlte,inb)*   &
			cdiag(:,:,KS:ixprnlte))
        phifx(:,:,KS:ixprnlte) = 0.0E+00
        do k=KS,ixprnlte
          do kp=ixprnlte+1,KE
            phifx(:,:,k) = phifx(:,:,k) +   &
                           cmtrx(:,:,kp,k)*sorc(:,:,kp,inb           )
          end do
        end do
        az(:,:,KS:ixprnlte) = sorc (:,:,KS:ixprnlte,inb           ) +  &
                     fnlte(:,:,KS:ixprnlte,inb)*phifx(:,:,KS:ixprnlte)

!----------------------------------------------------------------------
!     first iteration. (J(k) = B(k)) as initial guess)
!-----------------------------------------------------------------------
        phivar(:,:,KS:ixprnlte) = 0.0E+00
        do k=KS,ixprnlte
          do kp=KS,ixprnlte
            phivar(:,:,k) = phivar(:,:,k) +   &
                            cmtrx(:,:,kp,k)*sorc(:,:,kp,inb           )
          end do
        end do
        ag  (:,:,KS:ixprnlte) = fnlte(:,:,KS:ixprnlte,inb)*   &
                                (phivar(:,:,KS:ixprnlte) -   &
                                 cdiag(:,:,KS:ixprnlte)*  &
				 sorc(:,:,KS:ixprnlte,inb           ))

        sorc(:,:,KS:ixprnlte,inb           ) = bdenom(:,:,KS:ixprnlte)*&
                                               (az(:,:,KS:ixprnlte) + &
						ag(:,:,KS:ixprnlte)) 

!-----------------------------------------------------------------------
!     second iteration.  (J(k) = result of first iteration as guess)
!-----------------------------------------------------------------------
        phivar(:,:,KS:ixprnlte) = 0.0E+00
        do k=KS,ixprnlte
          do kp=KS,ixprnlte
            phivar(:,:,k) = phivar(:,:,k) +    &
                            cmtrx(:,:,kp,k)*sorc(:,:,kp,inb           )
          end do
        end do
        ag  (:,:,KS:ixprnlte) = fnlte(:,:,KS:ixprnlte,inb)*   &
                        (phivar(:,:,KS:ixprnlte) -   &
            cdiag(:,:,KS:ixprnlte)*sorc(:,:,KS:ixprnlte,inb           ))

        sorc(:,:,KS:ixprnlte,inb           ) = bdenom(:,:,KS:ixprnlte)*&
                                               (az(:,:,KS:ixprnlte) +  &
						ag(:,:,KS:ixprnlte)) 
      enddo

!-----------------------------------------------------------------------


end subroutine nlte



!#####################################################################

subroutine co2curt (pflux, cmtrx, Gas_tf)

!----------------------------------------------------------------------
!     co2curt computes Curtis matrix elements derived from co2
!     transmission functions.
!     functions.
!
!     author: m. d. schwarzkopf
!
!     revised: 8/18/94
!
!     certified:  radiation version 1.0
!
!----------------------------------------------------------------------
real, dimension(:,:, :), intent(in)  ::  pflux                  
type(gas_tf_type), intent(in) :: Gas_tf
real, dimension(:,:, :,:), intent(out)  ::  cmtrx                  

!-----------------------------------------------------------------------
!     intent out:
!
!       co21c  = column of transmission functions.
!
!       co21r  = row of transmission functions. 
!
!       cmtrx  = cutris matrix.
!-----------------------------------------------------------------------
     integer                               :: k, krow, kp
     real, dimension(size(pflux,1), size(pflux,2), &
                   size(pflux,3)-1) :: pdfinv
     real, dimension(size(pflux,1), size(pflux,2), &
                   size(pflux,3)) :: co2row, co2rowp

!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!     compute co2 transmission functions.
!-----------------------------------------------------------------------
      co2row(:,:,KS:KE+1)  = 1.0E+00
      co2rowp(:,:,KS:KE+1) = 1.0E+00

!-----------------------------------------------------------------------
!    compute curtis matrix for rows from KS to ixprnlte
!-----------------------------------------------------------------------
      do k = KS,ixprnlte
	krow = k
        pdfinv(:,:,k) = 1.0/(pflux(:,:,k+1) - pflux(:,:,k))

        call transcol ( KS, krow, KS, KE+1, co2row, Gas_tf)        

        call transcol ( KS, krow+1, KS, KE+1, co2rowp, Gas_tf)        

        do kp=KS,KE-1 
          cmtrx(:,:,kp,k) = radcon*pdfinv(:,:,k)*   &
                            (co2rowp(:,:,kp) - co2rowp(:,:,kp+1) -  &
                             co2row(:,:,kp) + co2row(:,:,kp+1)) 
        end do

        cmtrx(:,:,KE,k) = radcon*pdfinv(:,:,k)*   &
                          (co2rowp(:,:,KE) - co2row(:,:,KE)) 
      enddo

!--------------------------------------------------------------------


end subroutine co2curt




!#####################################################################

subroutine sealw99_alloc (ix, jx, kx, Lw_diagnostics)

!--------------------------------------------------------------------
!    longwave_driver_alloc allocates and initializes the components
!    of the lw_output_type variable Lw_output which holds the longwave
!    output needed by radiation_driver_mod.
!--------------------------------------------------------------------

integer,                   intent(in)  :: ix, jx, kx
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics

!--------------------------------------------------------------------
!   intent(in) variables:
!
!      ix,jx,kx     (i,j,k) lengths of radiation arrays to be allocated
!
!
!   intent(out) variables:
!
!      Lw_diagnostics
!                   lw_diagnostics_type variable containing diagnostic
!                   longwave output used by the radiation diagnostics
!                   module
!  
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      integer ::  NBTRGE, NBLY
!     integer ::  NBLY

    
!---------------------------------------------------------------------
!    allocate (and initialize where necessary) lw_diagnostics_type 
!    component arrays.
!---------------------------------------------------------------------
      NBTRGE = Lw_parameters%NBTRGE
      NBLY   = Lw_parameters%NBLY

      allocate ( Lw_diagnostics%flx1e1   (ix, jx                ) )
      allocate ( Lw_diagnostics%fluxn    (ix, jx, kx+1, 6+NBTRGE) )
      allocate (Lw_diagnostics%cts_out   (ix, jx, kx,   6       ) )
      allocate (Lw_diagnostics%cts_outcf (ix, jx, kx,   6       ) )
      allocate (Lw_diagnostics%gxcts     (ix, jx                ) )
      allocate (Lw_diagnostics%excts     (ix, jx, kx            ) )
      allocate (Lw_diagnostics%exctsn    (ix, jx, kx,   NBLY    ) )
      allocate (Lw_diagnostics%fctsg     (ix, jx,       NBLY    ) )

      Lw_diagnostics%flx1e1   = 0.
      Lw_diagnostics%cts_out    = 0.
      Lw_diagnostics%cts_outcf = 0.
      Lw_diagnostics%gxcts    = 0.
      Lw_diagnostics%excts  = 0.
      Lw_diagnostics%exctsn   = 0.
      Lw_diagnostics%fctsg   = 0.

      Lw_diagnostics%fluxn  (:,:,:,:) = 0.0

      if (Rad_control%do_totcld_forcing) then
        allocate ( Lw_diagnostics%fluxncf (ix, jx, kx+1, 6+NBTRGE) )
        Lw_diagnostics%fluxncf(:,:,:,:) = 0.0
      endif

      if (Lw_control%do_ch4_n2o) then
        allocate( Lw_diagnostics%flx1e1f  (ix, jx,       NBTRGE  ) )
         Lw_diagnostics%flx1e1f  = 0.
      end if

!--------------------------------------------------------------------

end subroutine sealw99_alloc


!####################################################################

subroutine co2_time_vary ( rrvco2            )

!---------------------------------------------------------------------
real, intent(in   )    ::  rrvco2

!---------------------------------------------------------------------
        real    ::    co2_vmr

        co2_vmr = rrvco2*1.0E+06

        call co2_lblinterp  (co2_vmr            )



end subroutine co2_time_vary



!####################################################################

subroutine ch4_n2o_time_vary (rrvch4, rrvn2o)

!---------------------------------------------------------------------
real, intent(in) :: rrvch4, rrvn2o               

!---------------------------------------------------------------------
     real             ::  ch4_vmr, n2o_vmr

!---------------------------------------------------------------------
!  the ch4 volume mixing ratio is set to the initial value (rch4) and 
!  the mass mixing ratio is defined on the first access of this routine.
!  then the lbl transmission function is calculated. after first access,
!  this routine does nothing. 
!--------------------------------------------------------------------
         ch4_vmr = rrvch4*1.0E+09
         call Ch4_lblinterp  (ch4_vmr)

!---------------------------------------------------------------------
!  the n2o volume mixing ratio is set to initial value (rn2o) and the 
!  mass mixing ratio is defined on the first access of this routine. 
!  routines are called to calculate the lbl transmission functions for 
!  n2o. after first access, this routine does nothing. 
!--------------------------------------------------------------------
         n2o_vmr = rrvn2o*1.0E+09
         call N2o_lblinterp (n2o_vmr)


end subroutine ch4_n2o_time_vary



!####################################################################




		  end module sealw99_mod


