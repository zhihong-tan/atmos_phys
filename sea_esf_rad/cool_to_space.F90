
                 module cool_to_space_mod


use rad_step_setup_mod,   only: KS=>KSRAD, KE=>KERAD, ISRAD, &
			        IERAD, JSRAD, JERAD, pflux
use longwave_params_mod,  only: NBLY_ORIG, NBLY_CKD2P1
use  rad_utilities_mod,   only: longwave_control_type, Lw_control, &
				looktab, longwave_tables3_type, &
				radiation_control_type, Rad_control
use  utilities_mod,       only: open_file, file_exist,       &
                                check_nml_error, error_mesg, &
			        print_version_number, FATAL, NOTE, &
				WARNING, get_my_pe, close_file, &
				read_data, write_data
use constants_new_mod,    only: diffac, radcon
use optical_path_mod,     only: get_var1var2, get_totvo2,   &
				get_totch2o, get_totch2obd
use longwave_aerosol_mod, only: get_totaerooptdep
use longwave_setup_mod,   only: pdfinv, ixoe1, dte1, &
				longwave_parameter_type,  &
				Lw_parameters
use co2_source_mod,       only: get_sorc
use cfc_mod,              only: cfc_exact
use longwave_tables_mod,  only: tabsr
use radiation_diag_mod,   only: radiag_from_cts, radiag_from_ctsappx
use gas_tf_mod,           only: get_gastf_for_cts

!-------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                  cool_to_space heating rates module
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.08'
    character(len=128)  :: version =  '$Id: cool_to_space.F90,v 1.2 2001/08/30 15:14:30 fms Exp $'
    character(len=128)  :: tag     =  '$Name: eugene $'

!---------------------------------------------------------------------
!-------  interfaces --------

public      &
	  cool_to_space_init, cool_to_space_approx,   &
	  cool_to_space_exact, &
	  cool_to_space_alloc, cool_to_space_dealloc, &
	  cool_to_space_get_hra

!---------------------------------------------------------------------
!-------- namelist  ---------

real                :: dummy = 1.0

namelist / cool_to_space_nml /        &
                               dummy


!---------------------------------------------------------------------
!------- public data ------



!---------------------------------------------------------------------
!------- private data ------


!---------------------------------------------------------------------
!     apcm, bpcm   =   capphi coefficients for NBLY bands.
!     atpcm, btpcm =   cappsi coefficients for NBLY bands.
!     acomb        =   random "a" parameter for NBLY bands.
!     bcomb        =   random "b" parameter for NBLY bands.
!---------------------------------------------------------------------
      real, dimension (:), allocatable    ::  apcm, bpcm, atpcm, btpcm,&
					      acomb, bcomb

      real, dimension(:,:,:), allocatable :: cts_sum, cts_sumcf

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
!---------------------------------------------------------------------



                           contains


subroutine cool_to_space_init

!---------------------------------------------------------------------
       real, dimension(NBLY_ORIG)   :: apcm_orig, bpcm_orig,     &
				       atpcm_orig, btpcm_orig,   &
				       acomb_orig, bcomb_orig
       real, dimension(NBLY_CKD2P1) :: apcm_ckd, bpcm_ckd, atpcm_ckd,  &
			               btpcm_ckd, acomb_ckd, bcomb_ckd 

       integer                      ::  k, m, n

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

       integer unit, ierr, io
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
         read (unit, nml=cool_to_space_nml, iostat=io, end=10)
         ierr = check_nml_error (io, 'cool_to_space_nml')
         enddo
 10      call close_file (unit)
       endif

       unit = open_file ('logfile.out', action='append')
!      call print_version_number (unit, 'cool_to_space', version_number)
       if (get_my_pe() == 0) then
	 write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	 write (unit,nml=cool_to_space_nml)
       endif
       call close_file (unit)

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

      if (Lw_parameters%NLWCLDB == 7) then
        cld_indx_table = cld_indx_table_lwclde
      else
        cld_indx_table = cld_indx_table_orig
      endif

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

end subroutine cool_to_space_init



!####################################################################

subroutine cool_to_space_approx (index, source, trans, iof, cld_trans, &
                                 trans2, iof2)

!---------------------------------------------------------------------
integer,                  intent(in)           ::  iof, index
real, dimension (:,:,:),  intent(in)           ::  source, trans, &
					           cld_trans
integer,                  intent(in), optional ::  iof2
real, dimension (:,:,:),  intent(in), optional ::  trans2

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
    real, dimension(:,:,:), allocatable   ::  cts_out, cts_outcf
    integer                               ::  j

!---------------------------------------------------------------------
    allocate ( cts_out  (ISRAD:IERAD, JSRAD:JERAD,     KS:KE) )
    allocate ( cts_outcf(ISRAD:IERAD, JSRAD:JERAD,     KS:KE) )

    if (present(trans2)) then
      cts_out(:,:,KS) = radcon*pdfinv(:,:,KS)*source(:,:,KS)*  &
                        (trans(:,:,KS+1-iof)*cld_trans(:,:,KS+1) -   &
                        trans2(:,:,KS       )*cld_trans(:,:,KS) )
    else
      cts_out(:,:,KS) = radcon*pdfinv(:,:,KS)*source(:,:,KS)*  &
                        (trans(:,:,KS+1-iof)*cld_trans(:,:,KS+1) - 1.0)
    endif

    if (present(trans2)) then
      cts_out(:,:,KS+1:KE) = radcon*pdfinv(:,:,KS+1:KE)*     &
			     source(:,:,KS+1:KE)*  &
                             (trans(:,:,KS+2-iof:KE+1-iof)*   &
	   		     cld_trans(:,:,KS+2:KE+1) -   &
                             trans2(:,:,KS+1-iof2 :KE-iof2)*     &
			     cld_trans(:,:,KS+1:KE))
    else
      cts_out(:,:,KS+1:KE) = radcon*pdfinv(:,:,KS+1:KE)*     &
			     source(:,:,KS+1:KE)*  &
                             (trans(:,:,KS+2-iof:KE+1-iof)*   &
			     cld_trans(:,:,KS+2:KE+1) -   &
                             trans(:,:,KS+1-iof :KE-iof)*     &
			     cld_trans(:,:,KS+1:KE))
    endif

    if (Rad_control%do_totcld_forcing) then
      if (present(trans2)) then
        cts_outcf(:,:,KS) = radcon*pdfinv(:,:,KS)*source(:,:,KS)* &
                            (trans(:,:,KS+1-iof) -       &
		  	     trans2(:,:,KS  -iof2))
      else
        cts_outcf(:,:,KS) = radcon*pdfinv(:,:,KS)*source(:,:,KS)* &
                            (trans(:,:,KS+1-iof) -  1.0)
      endif

      if (present(trans2)) then
        cts_outcf(:,:,KS+1:KE) = radcon*pdfinv(:,:,KS+1:KE)*     &
			         source(:,:,KS+1:KE)* &
                                 (trans(:,:,KS+2-iof:KE+1-iof) -     &
			          trans2(:,:,KS+1-iof2:KE-iof2))
      else
        cts_outcf(:,:,KS+1:KE) = radcon*pdfinv(:,:,KS+1:KE)*     &
			         source(:,:,KS+1:KE)* &
                                 (trans(:,:,KS+2-iof:KE+1-iof) -     &
			          trans(:,:,KS+1-iof:KE-iof))
      endif
    endif

!--------------------------------------------------------------------
!     send the contribution from the cool_to_space_approx calculation 
!     to radiation_diag_mod.
!--------------------------------------------------------------------
    if (Rad_control%do_diagnostics) then
      do j=JSRAD,JERAD
        call radiag_from_ctsappx (j, cts_out, index)
      end do
    endif

!--------------------------------------------------------------------
!     cts_sum is the sum of the values from cool_to_space_exact and
!     the values defined here in cool_to_space_approx. it will be used
!     by longwave_driver_mod.
!--------------------------------------------------------------------
    cts_sum(:,:,:) = cts_sum(:,:,:) - cts_out(:,:,:)
    cts_sumcf(:,:,:) = cts_sumcf(:,:,:) - cts_outcf(:,:,:)

!----------------------------------------------------------------
    deallocate (cts_out  )
    deallocate (cts_outcf)

!--------------------------------------------------------------------

end subroutine cool_to_space_approx



!####################################################################

subroutine cool_to_space_exact (cldtf, press, temp , to3cnt, gxcts,  &
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

real, dimension (:,:,:,:), intent(in)     :: cldtf
real, dimension (:,:,:),   intent(in)     :: temp, press, to3cnt
real, dimension(:,:),      intent(inout)  :: gxcts, gxctscf

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
      real, dimension(:,:),     allocatable   :: gxctsbd, pfac1,   &
						 pfac2, gxctsbdcf
      real, dimension(:,:,:),   allocatable   ::    &
	                 ag, agg, f, fac1, fac2, ff, phitmp, psitmp,  &
			 tmp1, tmp2, topm, topphi, tt, x, y, ctmp, &
			 cfc_tf, totch2o_tmp, totaer_tmp, &
			 excts, exctscf, fctsgcf, fctsg, n2o9c, &
			 totvo2_tmp, var1, var2, sorc_tmp     
      real, dimension(:,:,:,:), allocatable   :: co2spnb, exctsncf,  &
						 exctsn
      integer        :: n, k, j, ioffset

!----------------------------------------------------------------------
      ioffset = Lw_parameters%offset

!-----------------------------------------------------------------------
!    allocate needed arrays 
!-----------------------------------------------------------------------
      allocate ( gxctsbd (ISRAD:IERAD, JSRAD:JERAD            ))
      allocate ( pfac1   (ISRAD:IERAD, JSRAD:JERAD            ))
      allocate ( pfac2   (ISRAD:IERAD, JSRAD:JERAD            ))
      allocate ( sorc_tmp(ISRAD:IERAD, JSRAD:JERAD, KS:KE+1   ))
      allocate ( co2spnb (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, 3))
      allocate (n2o9c    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1   )) 
      allocate (excts    (ISRAD:IERAD, JSRAD:JERAD, KS:KE     ))
      allocate (exctscf  (ISRAD:IERAD, JSRAD:JERAD, KS:KE     ))
      if (Rad_control%do_totcld_forcing) then
        allocate (gxctsbdcf(ISRAD:IERAD, JSRAD:JERAD ))
      endif
      if (Rad_control%do_diagnostics) then
        allocate (exctsn  (ISRAD:IERAD, JSRAD:JERAD, KS:KE, NBLY))
        allocate (fctsg   (ISRAD:IERAD, JSRAD:JERAD,        NBLY)) 
	if (Rad_control%do_totcld_forcing) then
          allocate (exctsncf(ISRAD:IERAD, JSRAD:JERAD, KS:KE, NBLY))
          allocate (fctsgcf (ISRAD:IERAD, JSRAD:JERAD,        NBLY))
        endif
      endif

      allocate (tt          (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (x           (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))       
      allocate (y           (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))       
      allocate (ctmp        (ISRAD:IERAD, JSRAD:JERAD,  KS:KE+1))       
      allocate (topm        (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))      
      allocate (topphi      (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (phitmp      (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (psitmp      (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (ag          (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (agg         (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (f           (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (ff          (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (tmp1        (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (tmp2        (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (fac1        (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (fac2        (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (totch2o_tmp (ISRAD:IERAD, JSRAD:JERAD,  KS:KE+1))
      allocate (totvo2_tmp  (ISRAD:IERAD, JSRAD:JERAD,  KS+1:KE+1))
      allocate (totaer_tmp  (ISRAD:IERAD, JSRAD:JERAD,  KS:KE+1))
      allocate (cfc_tf      (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (var1        (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))
      allocate (var2        (ISRAD:IERAD, JSRAD:JERAD,  KS:KE  ))

!-----------------------------------------------------------------------
!     initialize quantities.
!-----------------------------------------------------------------------
      excts(:,:,KS:KE) = 0.0E+00
      gxcts(:,:)       = 0.0E+00

      if (Rad_control%do_totcld_forcing) then
        exctscf(:,:,KS:KE) = 0.0E+00
        gxctscf(:,:)       = 0.0E+00
      endif

!---------------------------------------------------------------------
!   retrieve co2spnb, tn2o17 and n2o9c from gas_tf_mod.
!   retrieve var1, var2 from optical_path_mod.
!-----------------------------------------------------------------------
      call get_gastf_for_cts (co2spnb, n2o9c)
      call get_var1var2 (var1, var2)

!-----------------------------------------------------------------------
!     compute temperature quantities.
!-----------------------------------------------------------------------
      x(:,:,KS:KE) = temp(:,:,KS:KE) - 2.5E+02
      y(:,:,KS:KE) = x(:,:,KS:KE)*x(:,:,KS:KE)
      ctmp(:,:,KS) = 1.0E+00

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
        phitmp(:,:,KS:KE) = var1(:,:,KS:KE)*     &
                            ((((ag (:,:,KS:KE)*        &
                            ag (:,:,KS:KE))**2)**2)**2)
        psitmp(:,:,KS:KE) = var2(:,:,KS:KE)*     &
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
	    call get_totch2o (n, totch2o_tmp)
            tt(:,:,KS:KE) = EXP(-1.0*(tmp1(:,:,KS:KE) + diffac*   &
                            totch2o_tmp(:,:,KS+1:KE+1)))
          else                                   !  combined bands 1-4.
            tt(:,:,KS:KE) = EXP(-1.0*tmp1(:,:,KS:KE)) 
          endif

!-----------------------------------------------------------------------
	else &
        if (n .GE. band_no_start(2) .and. n .le. band_no_end(2)) then
          if (Lw_control%do_ckd2p1) then    !   combined bands 25-40.
	    call get_totch2o (n, totch2o_tmp)
            tt(:,:,KS:KE) = EXP(-1.0*(tmp1(:,:,KS:KE) + diffac*   &
                            totch2o_tmp(:,:,KS+1:KE+1)))
	  else                              !  combined bands 5-8.
	    call get_totvo2 (n, totvo2_tmp)
            tt(:,:,KS:KE) = EXP(-1.0*(tmp1(:,:,KS:KE) +  &
                                totvo2_tmp(:,:,KS+1:KE+1))) 
          endif

!-----------------------------------------------------------------------
	else &                            !  first co2 band
        if (n .GE. band_no_start(3) .and. n .le. band_no_end(3)) then
          if (Lw_control%do_ckd2p1) then           !  combined band 41 
	    call get_totch2obd (n-40, totch2o_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*     &
                            totch2o_tmp(:,:,KS+1:KE+1     )
          else                                     !   combined band 9.
	    call get_totvo2 (n, totvo2_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) +               &
                                totvo2_tmp(:,:,KS+1:KE+1)
          endif
	  if (Lw_control%do_lwaerosol) then
	    call get_totaerooptdep(1, totaer_tmp)
	    tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +        &
                            totaer_tmp  (:,:,KS+1:KE+1  )
          endif
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))*    &
                                co2spnb(:,:,KS+1:KE+1,1)

!-----------------------------------------------------------------------
	else  &                      ! second co2 band
        if (n .GE. band_no_start(4) .and. n .le. band_no_end(4)) then
          if (Lw_control%do_ckd2p1) then         !   combined band 42.
	    call get_totch2obd (n-40, totch2o_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*     &
                            totch2o_tmp(:,:,KS+1:KE+1     )
          else                                   ! combined band 10. 
	    call get_totvo2 (n, totvo2_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) +               &
                                totvo2_tmp(:,:,KS+1:KE+1)
          endif
	  if (Lw_control%do_lwaerosol) then
	    call get_totaerooptdep(2, totaer_tmp)
	    tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +          &
                            totaer_tmp   (:,:,KS+1:KE+1  )
          endif
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))*     &
                                co2spnb(:,:,KS+1:KE+1,2)

!-----------------------------------------------------------------------
	else &                                   !  third co2 band
        if (n .GE. band_no_start(5) .and. n .le. band_no_end(5)) then
          if (Lw_control%do_ckd2p1) then      !   combined band 43.
	    call get_totch2obd (n-40, totch2o_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*      &
                            totch2o_tmp(:,:,KS+1:KE+1     )
          else                                !   combined band 11. 
	    call get_totvo2 (n, totvo2_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) +                &
                                totvo2_tmp(:,:,KS+1:KE+1)
          endif
	  if (Lw_control%do_lwaerosol) then
	    call get_totaerooptdep(3, totaer_tmp)
	    tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +       &
                            totaer_tmp   (:,:,KS+1:KE+1  )
          endif
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))*     &
                                co2spnb(:,:,KS+1:KE+1,3)

!-----------------------------------------------------------------------
	else &
        if (n .GE. band_no_start(6) .and. n .le. band_no_end(6)) then
          if (Lw_control%do_ckd2p1) then     ! combined bands 44-45.
	    call get_totch2obd (n-40, totch2o_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*    &
                            totch2o_tmp(:,:,KS+1:KE+1     )
          else                               ! combined bands 12-13.
	    call get_totvo2 (n, totvo2_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) +               &
                                totvo2_tmp(:,:,KS+1:KE+1)
          endif
	  if (Lw_control%do_lwaerosol) then
	    call get_totaerooptdep(n-8-ioffset, totaer_tmp)
	    tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +         &
                            totaer_tmp   (:,:,KS+1:KE+1            )
          endif
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))

!-----------------------------------------------------------------------
	else   &                             !  combined band 14 or 46.
        if (n .GE. band_no_start(7) .and. n .le. band_no_end(7)) then
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp1(:,:,KS:KE))*      &
                               to3cnt (:,:,KS:KE)

!-----------------------------------------------------------------------
	else   &
        if (n .GE. band_no_start(8) .and. n .le. band_no_end(8)) then
          if (Lw_control%do_ckd2p1) then        !  combined band 47.
	    call get_totch2obd (n-40, totch2o_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*    &
                            totch2o_tmp(:,:,KS+1:KE+1     )
          else                                  !  combined band 15.
	    call get_totvo2 (n, totvo2_tmp)
            tmp2(:,:,KS:KE) = tmp1(:,:,KS:KE) +              &
                                totvo2_tmp(:,:,KS+1:KE+1)
          endif
	  if (Lw_control%do_lwaerosol) then
	    call get_totaerooptdep(7, totaer_tmp)
	    tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE)    +                &
                            totaer_tmp   (:,:,KS+1:KE+1  )
          endif
	  tt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))*   &
                                n2o9c  (:,:,KS+1:KE+1)
        endif
!--------------------------------------------------------------------
!     calculate or retrieve the source function for the current band.
!--------------------------------------------------------------------
        if (n <= 8 + ioffset) then
	  call looktab (tabsr, ixoe1, dte1, sorc_tmp, KS, KE+1, n)
        else
	  call get_sorc (sorc_tmp, n)
        endif

!---------------------------------------------------------------------
!     retrieve the cfc effect if cfcs are activated.
!---------------------------------------------------------------------
        if (Lw_control%do_cfc .and. n >= 9+ioffset) then
          call cfc_exact(n-8-ioffset, cfc_tf)
	  tt(:,:,KS:KE) = tt(:,:,KS:KE)*cfc_tf(:,:,KS:KE)
        endif

!------------------------------------------------------------------
!     define some near-surface pressure functions that are needed
!------------------------------------------------------------------
        pfac1(:,:) = 0.5*pdfinv(:,:,KE)*(pflux(:,:,KE+1) - &
                     press(:,:,KE))*tt(:,:,KE-1)
        pfac2(:,:) = 0.5*pdfinv(:,:,KE)*(pflux(:,:,KE+1) +    &
	             press(:,:,KE) - 2.0*pflux(:,:,KE))*tt(:,:,KE)

!--------------------------------------------------------------------
!     calculate the ground fluxes (?)
!--------------------------------------------------------------------
        if (Rad_control%do_totcld_forcing) then
          gxctsbdcf(:,:) = tt(:,:,KE)*sorc_tmp(:,:,KE) +   &
	                   (pfac1(:,:) + pfac2(:,:))*   &
		           (sorc_tmp(:,:,KE+1) - sorc_tmp(:,:,KE))
          gxctsbd(:,:) = cldtf(:,:,KE+1,cld_indx_table(n+32-ioffset))* &
	                 gxctsbdcf(:,:)
          gxctscf(:,:) = gxctscf(:,:) + gxctsbdcf(:,:)
          gxcts(:,:) = gxcts(:,:) + gxctsbd(:,:)
        else
          gxctsbd(:,:) = cldtf(:,:,KE+1,cld_indx_table(n+32-ioffset))* &
                         (tt(:,:,KE)*sorc_tmp(:,:,KE) +  &
	                 (pfac1(:,:) + pfac2(:,:))*   &
		         (sorc_tmp(:,:,KE+1) - sorc_tmp(:,:,KE)))
          gxcts(:,:) = gxcts(:,:) + gxctsbd(:,:)
        endif

        if (Rad_control%do_diagnostics) then
          fctsg(:,:,n) = gxctsbd(:,:)
	  if (Rad_control%do_totcld_forcing) then
            fctsgcf(:,:,n) = gxctsbdcf(:,:)
          endif
        endif

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
        if (Rad_control%do_diagnostics) then
          exctsn(:,:,KS:KE,n) = sorc_tmp(:,:,KS:KE)*     &
                                (ctmp(:,:,KS+1:KE+1) - ctmp(:,:,KS:KE))
          if (Rad_control%do_totcld_forcing) then
	    exctsncf(:,:,KS,n) =  sorc_tmp(:,:,KS)*         &
                                  (tt(:,:,KS) - 1.0E+00)
            exctsncf(:,:,KS+1:KE,n) = sorc_tmp(:,:,KS+1:KE)*     &
                                     (tt(:,:,KS+1:KE) - tt(:,:,KS:KE-1))
          endif
        endif

!-----------------------------------------------------------------------
!     excts is the cool-to-space cooling rate accumulated over
!     frequency bands.
!-----------------------------------------------------------------------
        excts(:,:,KS:KE) = excts(:,:,KS:KE) + sorc_tmp(:,:,KS:KE)*    &
                           (ctmp(:,:,KS+1:KE+1) - ctmp(:,:,KS:KE))
        if (Rad_control%do_totcld_forcing) then
          exctscf(:,:,KS) = exctscf(:,:,KS) + sorc_tmp(:,:,KS)*     &
                            (tt(:,:,KS) - 1.0E+00)
          exctscf(:,:,KS+1:KE) = exctscf(:,:,KS+1:KE) +     &
                                 sorc_tmp(:,:,KS+1:KE)*      &
                                 (tt(:,:,KS+1:KE) - tt(:,:,KS:KE-1))
        endif
      end do

!-----------------------------------------------------------------------
      deallocate (cfc_tf    )
      deallocate ( totaer_tmp )
      deallocate ( totch2o_tmp )
      deallocate ( totvo2_tmp )
      deallocate (fac2)
      deallocate (fac1)
      deallocate (tmp2)
      deallocate (tmp1)
      deallocate  (ff)
      deallocate (f)
      deallocate (agg)
      deallocate (ag)
      deallocate (psitmp    )
      deallocate (phitmp    )
      deallocate (topphi    )
      deallocate (topm      )
      deallocate (var1      )
      deallocate (var2      )

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
        gxcts(:,:) = gxcts(:,:) - excts(:,:,k)
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
      if (Rad_control%do_diagnostics) then
        do n=1,NBLY-1
	  exctsn(:,:,KS:KE,n) = exctsn(:,:,KS:KE,n)*radcon*     &
                                pdfinv(:,:,KS:KE)
        enddo
        if (Rad_control%do_totcld_forcing) then
          do n=1,NBLY-1
	    exctsncf(:,:,KS:KE,n) = exctsncf(:,:,KS:KE,n)*radcon*    &
                                    pdfinv(:,:,KS:KE)
          enddo
        endif
      endif
      excts(:,:,KS:KE) = excts(:,:,KS:KE)*radcon*pdfinv(:,:,KS:KE)

      if (Rad_control%do_totcld_forcing) then
        exctscf(:,:,KS:KE) = exctscf(:,:,KS:KE)*radcon*pdfinv(:,:,KS:KE)
      endif

!-----------------------------------------------------------------------
!     this is the end of the exact cool-to-space computations; at this
!     point excts has its appropriate value.
!-----------------------------------------------------------------------
      deallocate (ctmp      )
      deallocate (y         )
      deallocate (x         )
      deallocate (tt        )

!--------------------------------------------------------------------
      if (Rad_control%do_diagnostics) then
        if (Rad_control%do_totcld_forcing) then
	  deallocate (fctsgcf   )
	  deallocate (exctsncf  )
        endif
      endif

!--------------------------------------------------------------------
!    send data to radiation_diag_mod.
!--------------------------------------------------------------------
      if (Rad_control%do_diagnostics) then
        do j=JSRAD,JERAD
          call radiag_from_cts (j, exctsn, fctsg, excts, gxcts)
        end do
        deallocate (fctsg     )
        deallocate (exctsn    )
      endif

!--------------------------------------------------------------------
!    save the heating rates to be later sent to longwave_driver_mod.
!--------------------------------------------------------------------
      cts_sum(:,:,:) = excts(:,:,:)
      cts_sumcf(:,:,:) = exctscf(:,:,:)

!--------------------------------------------------------------------
      deallocate (exctscf)
      deallocate (excts  )
      deallocate (n2o9c  )

      if (Rad_control%do_totcld_forcing) then
        deallocate (gxctsbdcf)
      endif

      deallocate (co2spnb  )
      deallocate (sorc_tmp)
      deallocate (pfac2   )
      deallocate (pfac1   )
      deallocate (gxctsbd )

!--------------------------------------------------------------------


end  subroutine cool_to_space_exact




!####################################################################

subroutine cool_to_space_alloc 

!-------------------------------------------------------------------
        allocate (cts_sum (ISRAD:IERAD, JSRAD:JERAD, KS:KE ))

        if (Rad_control%do_totcld_forcing) then
          allocate (cts_sumcf (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ))
        endif
!--------------------------------------------------------------------


end subroutine cool_to_space_alloc



!####################################################################

subroutine cool_to_space_dealloc 

!--------------------------------------------------------------------
   deallocate (cts_sum)

   if (Rad_control%do_totcld_forcing) then
     deallocate (cts_sumcf)
   endif
!-------------------------------------------------------------------

end subroutine cool_to_space_dealloc 


!####################################################################

subroutine cool_to_space_get_hra (cts_sumd, cts_sumcfd)

!---------------------------------------------------------------------
real, dimension(:,:,:),           intent(out)  :: cts_sumd
real, dimension(:,:,:), optional, intent(out)  :: cts_sumcfd

!------------------------------------------------------------------
!     send the summed cool_to_space heating rates from _exact and 
!     _approx to the longwave_driver_mod.
!------------------------------------------------------------------
    cts_sumd  (:,:,:) = cts_sum  (:,:,:)

    if ( present(cts_sumcfd) ) then
      cts_sumcfd(:,:,:) = cts_sumcf(:,:,:)
    endif

end subroutine cool_to_space_get_hra 



!####################################################################



	            end module cool_to_space_mod
