
		 module longwave_tables_mod


use rad_step_setup_mod,    only: KS=>KSRAD, KE=>KERAD, ISRAD, IERAD, &
				 JSRAD, JERAD
use longwave_params_mod,   only: NBLW, NBLX, NBLY_ORIG, NBLY_CKD2P1, &
				 NBCO215
use  utilities_mod,        only: open_file, file_exist,    &
                                 check_nml_error, error_mesg, &
                                 print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file
use radiation_diag_mod,    only: radiag_from_lwtables
use ch4_n2o_mod,           only: ch4_n2o_input_type, CN_basic
use optical_path_mod,      only: get_path_for_enear, get_path_for_e90 
use longwave_setup_mod,    only: tpl1, tpl2, dte1, dte2, ixoe1, &
		                 ixoe2, temp_1, mass_1,    &
				 Lw_parameters,longwave_parameter_type
use rad_utilities_mod,     only: locate_in_table, looktab,   &
		                 longwave_tables1_type,  &
				 longwave_tables2_type,  &
	                  	 longwave_tables3_type,  &
				 table_alloc, table_axis_type, &
		                 longwave_control_type, Lw_control


!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                 longwave tables module
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: longwave_tables.F90,v 1.2 2001/08/30 15:13:04 fms Exp $'
   character(len=128)  :: tag     =  '$Name: galway $'

!---------------------------------------------------------------------
!------    interfaces   ------

public   &
       	   longwave_tables_init, &
	   e1e290, e290, enear, &
	   get_bds

private  &
	   idrbtsh2o, id2h2o, table

!---------------------------------------------------------------------
!------     namelist  -----

real                     :: dummy = 1.0

namelist / longwave_tables_nml /  &
                                   dummy


!---------------------------------------------------------------------
!---- public data -------


type (longwave_tables3_type), public     :: tabsr


!---------------------------------------------------------------------
!---- private data -------

 
!--------------------------------------------------------------------
!  define continuum coefficients over special bands, the choices 
!  depend on model architecture. the program gasbnd is used.
!--------------------------------------------------------------------
real, dimension(:), allocatable :: afah2o
real, dimension(:), allocatable :: afach4, afan2o
              
real, dimension(:), allocatable :: bdlah2o, bdhah2o
real, dimension(:), allocatable :: bdlahcn, bdhahcn

real, dimension(:), allocatable :: bfah2o             
real, dimension(:), allocatable :: bfach4, bfan2o             

real, dimension(:), allocatable :: dch4, dn2o, ech4, en2o
real                            :: d171n2o, e171n2o



real, dimension (:), allocatable :: acomb, bcomb, apcm, bpcm, atpcm,  &
				    btpcm, bdlocm, bdhicm

type (longwave_tables1_type)     :: tab1, tab2, tab3, tab1w
type (longwave_tables2_type)     :: tab1a, tab2a, tab3a

integer, parameter               :: NTTABH2O   = 28
integer, parameter               :: NUTABH2O   = 180

real, dimension (NBLW)           :: bandlo, bandhi, arndm, brndm, betad
integer, dimension(40)           :: iband
real, dimension(3)               :: ao3cm, bo3cm
real, dimension(2)               :: ab15cm

integer                          :: NBTRG, NBTRGE, NBLY
real                             :: apwd, bpwd, atpwd, btpwd, bdlowd, &
				    bdhiwd 


!---------------------------------------------------------------------
!---------------------------------------------------------------------




contains



subroutine longwave_tables_init

!-----------------------------------------------------------------------
    integer              :: k4, n4

!---------------------------------------------------------------------
!  define continuum coefficients over special bands, the choices 
!  depend on model architecture. the program gasbnd is used.
!---------------------------------------------------------------------
    real                          :: apwd_c, bpwd_c, atpwd_c,    &
				     btpwd_c, bdlowd_c, bdhiwd_c
    integer, dimension(40)        :: iband_c
    real, dimension (NBLY_CKD2P1) :: acomb_c, bcomb_c, apcm_c,  &
				     bpcm_c, atpcm_c,   &
	                	     btpcm_c, bdlocm_c,  bdhicm_c

 
!    1)  560-800 as 1 band

      data apwd_c            /   &
         0.180524E-01/
      data bpwd_c            /  &
        -0.556925E-04/
      data atpwd_c            /  &
         0.186362E-01/
      data btpwd_c            /   &
        -0.572090E-04/
      data bdlowd_c /    &
         0.560000E+03/
      data bdhiwd_c /    &
         0.800000E+03/
 
!    2) 160-560 (as 40 bands). program gasbnd is used with 10 cm-1
!    bandwidth. iband is straightforward mapping.
 
      data iband_c /    &
          1,   2,   3,   4,   5,   6,   7,   8,   9,  10,   &
         11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  &
         21,  22,  23,  24,  25,  26,  27,  28,  29,  30,   &
         31,  32,  33,  34,  35,  36,  37,  38,  39,  40/ 
!
      data (acomb_c(k4),k4=1,40)    /   &
         0.849130E+03,  0.135587E+05,  0.286836E+04,  0.169580E+04, &
         0.208642E+05,  0.126034E+04,  0.109494E+05,  0.335111E+03,  &
         0.488806E+04,  0.860045E+04,  0.537333E+03,  0.437769E+04, &
         0.345836E+04,  0.129538E+03,  0.463562E+04,  0.251489E+03, &
         0.256264E+04,  0.485091E+03,  0.889881E+03,  0.116598E+04,  &
         0.125244E+03,  0.457264E+03,  0.142197E+03,  0.444551E+03, &
         0.301446E+02,  0.392750E+03,  0.436426E+02,  0.347449E+02,  &
         0.612509E+02,  0.142506E+03,  0.103643E+02,  0.721701E+02,  &
         0.315040E+02,  0.941759E+01,  0.540473E+02,  0.350084E+02,  &
         0.300816E+02,  0.379020E+01,  0.125727E+02,  0.545869E+01/

      data (bcomb_c(k4),k4=1,40)    /   &
         0.175174E+00,  0.176667E+00,  0.109512E+00,  0.111893E+00,  &
         0.145289E+00,  0.203190E+00,  0.151547E+00,  0.911103E-01,  &
         0.151444E+00,  0.850764E-01,  0.756520E-01,  0.100377E+00,  &
         0.171557E+00,  0.125429E+00,  0.105915E+00,  0.816331E-01,  &
         0.149472E+00,  0.857054E-01,  0.107092E+00,  0.185458E+00,  &
         0.753818E-01,  0.108639E+00,  0.123934E+00,  0.178712E+00,  &
         0.833855E-01,  0.119886E+00,  0.133082E+00,  0.935851E-01,  &
         0.156848E+00,  0.166457E+00,  0.162215E+00,  0.114845E+00, &
         0.724304E-01,  0.740525E-01,  0.734090E-01,  0.141319E+00,  &
         0.359408E-01,  0.833614E-01,  0.128919E+00,  0.996329E-01/

      data (apcm_c(k4),k4=1,40)    /   &
         0.549325E-02, -0.150653E-02,  0.268788E-02,  0.138495E-01,  &
        -0.714528E-03,  0.112319E-01,  0.113418E-02,  0.215116E-01, &
         0.388898E-02,  0.398385E-02,  0.931768E-02,  0.655185E-02, &
         0.735642E-02,  0.190346E-01,  0.104483E-01,  0.917671E-02, &
         0.108668E-01,  0.305797E-02,  0.163975E-01,  0.147718E-01, &
         0.485502E-02,  0.223258E-01,  0.567357E-02,  0.197808E-01,  &
         0.245634E-01,  0.116045E-01,  0.269989E-01,  0.176298E-01, &
         0.128961E-01,  0.134788E-01,  0.391238E-01,  0.117165E-01,  &
         0.691808E-02,  0.202443E-01,  0.137798E-01,  0.215153E-01, &
         0.154358E-01,  0.850256E-02,  0.111306E-01,  0.185757E-01/

      data (bpcm_c(k4),k4=1,40)    /   &
        -0.305151E-04, -0.244741E-05, -0.203093E-04, -0.736015E-04,  &
        -0.158662E-04, -0.381826E-04, -0.197166E-04, -0.984160E-04,  &
        -0.222455E-04, -0.346880E-04, -0.395593E-04, -0.426165E-04,  &
        -0.410312E-04, -0.848479E-04, -0.597304E-04, -0.318474E-04, &
        -0.450295E-04, -0.284497E-04, -0.772035E-04, -0.545821E-04, &
        -0.242272E-04, -0.105653E-03, -0.854473E-05, -0.672510E-04, &
        -0.109627E-03, -0.330446E-04, -0.682675E-04,  0.479154E-04,  &
         0.411211E-04, -0.554504E-04, -0.145967E-03, -0.425913E-04,  &
         0.413272E-05, -0.531586E-04, -0.429814E-04, -0.847248E-04, &
        -0.733456E-04,  0.403362E-05, -0.389712E-04, -0.531450E-04/

      data (atpcm_c(k4),k4=1,40)    /   &
         0.541966E-02, -0.153876E-02,  0.158818E-02,  0.133698E-01,  &
        -0.270746E-02,  0.691660E-02,  0.485749E-05,  0.199036E-01, &
         0.319826E-02,  0.220802E-02,  0.985921E-02,  0.462151E-02, &
         0.554947E-02,  0.149315E-01,  0.911982E-02,  0.696417E-02, &
         0.892579E-02,  0.222579E-02,  0.123105E-01,  0.129295E-01, &
         0.745423E-02,  0.203878E-01,  0.597427E-02,  0.170838E-01, &
         0.231443E-01,  0.127692E-01,  0.222678E-01,  0.165331E-01, &
         0.141333E-01,  0.141307E-01,  0.312569E-01,  0.114137E-01, &
         0.126050E-01,  0.242966E-01,  0.145196E-01,  0.224802E-01, &
         0.150399E-01,  0.173815E-01,  0.103438E-01,  0.188690E-01/

      data (btpcm_c(k4),k4=1,40)    /   &
        -0.202513E-04,  0.948663E-06, -0.130243E-04, -0.678688E-04, &
        -0.986181E-05, -0.258818E-04, -0.139292E-04, -0.916890E-04, &
        -0.138148E-04, -0.275165E-04, -0.298451E-04, -0.310005E-04, &
        -0.314745E-04, -0.695971E-04, -0.486158E-04, -0.198383E-04, &
        -0.421494E-04, -0.102396E-04, -0.591516E-04, -0.575729E-04, &
        -0.238464E-04, -0.938688E-04, -0.885666E-05, -0.728295E-04, &
        -0.938897E-04, -0.448622E-04, -0.642904E-04, -0.102699E-04, &
        -0.348743E-05, -0.568427E-04, -0.104122E-03, -0.313943E-04, &
         0.939109E-05, -0.631332E-04, -0.325944E-04, -0.757699E-04, &
        -0.398066E-04,  0.285026E-04, -0.355222E-04, -0.266443E-04/

      data (bdlocm_c(k4),k4=1,40)  /   &
         0.160000E+03,  0.170000E+03,  0.180000E+03,  0.190000E+03, &
         0.200000E+03,  0.210000E+03,  0.220000E+03,  0.230000E+03, &
         0.240000E+03,  0.250000E+03,  0.260000E+03,  0.270000E+03, &
         0.280000E+03,  0.290000E+03,  0.300000E+03,  0.310000E+03, &
         0.320000E+03,  0.330000E+03,  0.340000E+03,  0.350000E+03,& 
         0.360000E+03,  0.370000E+03,  0.380000E+03,  0.390000E+03,&
         0.400000E+03,  0.410000E+03,  0.420000E+03,  0.430000E+03, &
         0.440000E+03,  0.450000E+03,  0.460000E+03,  0.470000E+03, &
         0.480000E+03,  0.490000E+03,  0.500000E+03,  0.510000E+03,&
         0.520000E+03,  0.530000E+03,  0.540000E+03,  0.550000E+03/

      data (bdhicm_c(k4),k4=1,40)  /   &
         0.170000E+03,  0.180000E+03,  0.190000E+03,  0.200000E+03,  &
         0.210000E+03,  0.220000E+03,  0.230000E+03,  0.240000E+03, &
         0.250000E+03,  0.260000E+03,  0.270000E+03,  0.280000E+03, &
         0.290000E+03,  0.300000E+03,  0.310000E+03,  0.320000E+03, &
         0.330000E+03,  0.340000E+03,  0.350000E+03,  0.360000E+03, &
         0.370000E+03,  0.380000E+03,  0.390000E+03,  0.400000E+03, &
         0.410000E+03,  0.420000E+03,  0.430000E+03,  0.440000E+03, &
         0.450000E+03,  0.460000E+03,  0.470000E+03,  0.480000E+03, &
         0.490000E+03,  0.500000E+03,  0.510000E+03,  0.520000E+03, &
         0.530000E+03,  0.540000E+03,  0.550000E+03,  0.560000E+03/
 
!    3) 560-1200 and 4.3 um band (8 bands, frequency range given
!    by bdlocm-bdhicm). program gasbnd is used with 8 specified
!    bandwidths.
 
      data (acomb_c(k4),k4=41,NBLY_CKD2P1)    /   &
         0.788160E+01,  0.207769E+01,  0.433980E+00,  0.499718E-01,  &
         0.166671E-01,  0.213073E-01,  0.210228E+00,  0.458804E-02/

      data (bcomb_c(k4),k4=41,NBLY_CKD2P1)    / &
         0.676888E-01,  0.690856E-01,  0.660847E-01,  0.730849E-01,   &
         0.654438E-01,  0.848294E-01,  0.852442E-01,  0.242787E+00/

      data (apcm_c(k4),k4=41,NBLY_CKD2P1)    /  &
         0.168938E-01,  0.211805E-01,  0.224487E-01,  0.257983E-01,  &
         0.282393E-01,  0.293678E-01,  0.204784E-01,  0.383510E-01/

      data (bpcm_c(k4),k4=41,NBLY_CKD2P1)    /  &
        -0.560515E-04, -0.645812E-04, -0.559626E-04, -0.456728E-04,  &
        -0.789113E-04, -0.105722E-03, -0.811760E-04, -0.103490E-03/

      data (atpcm_c(k4),k4=41,NBLY_CKD2P1)    /  &
         0.167743E-01,  0.195884E-01,  0.228968E-01,  0.274106E-01,  &
         0.305775E-01,  0.292968E-01,  0.201658E-01,  0.373711E-01/

      data (btpcm_c(k4),k4=41,NBLY_CKD2P1)    /   &
        -0.510390E-04, -0.643921E-04, -0.722986E-04, -0.742483E-04, &
        -0.903545E-04, -0.102420E-03, -0.702446E-04, -0.120756E-03/

      data (bdlocm_c(k4),k4=41,NBLY_CKD2P1)  /   &
         0.560000E+03,  0.630000E+03,  0.700000E+03,  0.800000E+03,  &
         0.900000E+03,  0.990000E+03,  0.107000E+04,  0.227000E+04/

      data (bdhicm_c(k4),k4=41,NBLY_CKD2P1)  /  &
         0.630000E+03,  0.700000E+03,  0.800000E+03,  0.900000E+03,  &
         0.990000E+03,  0.107000E+04,  0.120000E+04,  0.238000E+04/

!----------------------------------------------------------------------

      real                        ::   apwd_n, bpwd_n, atpwd_n,   &
				       btpwd_n, bdlowd_n, bdhiwd_n
      integer, dimension(40)      ::   iband_n
      real, dimension (NBLY_ORIG) ::   acomb_n, bcomb_n, apcm_n,  &
				       bpcm_n, atpcm_n, btpcm_n,  &
				       bdlocm_n,  bdhicm_n
 
!    define random band parameters for special bands. the choices 
!    depend on model architecture. the program gasbnd is used.
 
!    1)  560-800 as 1 band
 
      data apwd_n            /  &
        0.180398E-01/
      data bpwd_n            /   &
       -0.556428E-04/
      data atpwd_n            /   &
        0.186251E-01/
      data btpwd_n            /   &
       -0.571712E-04/
      data bdlowd_n /   &
        0.560000E+03/
      data bdhiwd_n /   &
        0.800000E+03/
 
!    2) 160-560 (as 8 bands using combined bands). program gasbnd is
!    used as 40 bands (160-560,10 cm-1 bandwidth) with ifdef icomb
!    on. combined bands defined according to index iband.
 
      data iband_n /   &
          2,   1,   2,   2,   1,   2,   1,   3,   2,   2,   &
          3,   2,   2,   4,   2,   4,   2,   3,   3,   2,  &
          4,   3,   4,   3,   7,   5,   6,   7,   6,   5,  &
          7,   6,   7,   8,   6,   6,   8,   8,   8,   8/
 
      data (acomb_n(k4),k4=1,8)    /    &
         0.151896E+05,  0.331743E+04,  0.526549E+03,  0.162813E+03,   &
         0.268531E+03,  0.534040E+02,  0.267777E+02,  0.123003E+02/

      data (bcomb_n(k4),k4=1,8)    /   &
         0.153898E+00,  0.116352E+00,  0.102814E+00,  0.966280E-01,  &
         0.128586E+00,  0.119393E+00,  0.898891E-01,  0.645340E-01/

      data (apcm_n(k4),k4=1,8)    /   &
        -0.531605E-03,  0.679588E-02,  0.145133E-01,  0.928230E-02,   &
         0.120727E-01,  0.159210E-01,  0.175268E-01,  0.150281E-01/

      data (bpcm_n(k4),k4=1,8)    /   &
        -0.122002E-04, -0.335081E-04, -0.449999E-04, -0.231290E-04,  &
        -0.385298E-04, -0.157170E-04,  0.183434E-04, -0.501718E-04/

      data (atpcm_n(k4),k4=1,8)    /   &
        -0.156433E-02,  0.611348E-02,  0.134469E-01,  0.871421E-02,  &
         0.133191E-01,  0.164697E-01,  0.199640E-01,  0.163219E-01/

      data (btpcm_n(k4),k4=1,8)    / &
        -0.698856E-05, -0.295269E-04, -0.503674E-04, -0.268392E-04,  &
        -0.496663E-04, -0.333122E-04, -0.337346E-04, -0.258706E-04/

      data (bdlocm_n(k4),k4=1,8)  /  &
         0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  &
         0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/ 

      data (bdhicm_n(k4),k4=1,8)  /  &
         0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  &
         0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/ 
 
!    3) 560-1200 and 4.3 um band (8 bands, frequency range given
!    by bdlocm-bdhicm). program gasbnd is used with 8 specified
!    bandwidths.
 
      data (acomb_n(k4),k4=9,NBLY_ORIG)    /   &
         0.790655E+01,  0.208292E+01,  0.435231E+00,  0.501124E-01,  &
         0.167058E-01,  0.213731E-01,  0.211154E+00,  0.458804E-02/

      data (bcomb_n(k4),k4=9,NBLY_ORIG)    /  &
         0.676982E-01,  0.691090E-01,  0.660834E-01,  0.730815E-01,   &
         0.654709E-01,  0.848528E-01,  0.852302E-01,  0.242787E+00/

      data (apcm_n(k4),k4=9,NBLY_ORIG)    /  &
         0.168815E-01,  0.211690E-01,  0.224373E-01,  0.257869E-01,  &
         0.282310E-01,  0.293583E-01,  0.204657E-01,  0.383510E-01/

      data (bpcm_n(k4),k4=9,NBLY_ORIG)    /   &
        -0.560024E-04, -0.645285E-04, -0.559339E-04, -0.456444E-04,  &
        -0.788833E-04, -0.105692E-03, -0.811350E-04, -0.103490E-03/

      data (atpcm_n(k4),k4=9,NBLY_ORIG)    /   &
         0.167628E-01,  0.195779E-01,  0.228872E-01,  0.274006E-01,   &
         0.305678E-01,  0.292862E-01,  0.201527E-01,  0.373711E-01/

      data (btpcm_n(k4),k4=9,NBLY_ORIG)    /   &
        -0.510036E-04, -0.643516E-04, -0.722669E-04, -0.742132E-04,   &
        -0.903238E-04, -0.102388E-03, -0.702017E-04, -0.120756E-03/

      data (bdlocm_n(k4),k4=9,NBLY_ORIG)  /   &
         0.560000E+03,  0.630000E+03,  0.700000E+03,  0.800000E+03,  &
         0.900000E+03,  0.990000E+03,  0.107000E+04,  0.227000E+04/

      data (bdhicm_n(k4),k4=9,NBLY_ORIG)  /   &
         0.630000E+03,  0.700000E+03,  0.800000E+03,  0.900000E+03,   &
         0.990000E+03,  0.107000E+04,  0.120000E+04,  0.238000E+04/

      integer    ::  unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=longwave_tables_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'longwave_tables_nml')
       enddo
10     call close_file (unit)
     endif

     unit = open_file ( 'logfile.out', action='append')
!    call print_version_number (unit, 'longwave_tables', version_number)
     if (get_my_pe() == 0) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=longwave_tables_nml)
     endif
     call close_file (unit)

!----------------------------------------------------------------------
     NBTRG  = Lw_parameters%NBTRG
     NBTRGE = Lw_parameters%NBTRGE

     if (Lw_control%do_ch4_n2o) then
       allocate ( afah2o (NBTRGE) )
       allocate ( afach4 (NBTRG ) )
       allocate ( afan2o (NBTRG ) )
       allocate ( bdlah2o(NBTRGE) )
       allocate ( bdlahcn(NBTRG ) )
       allocate ( bdhah2o(NBTRGE) )
       allocate ( bdhahcn(NBTRG ) )
       allocate ( bfah2o (NBTRGE) )
       allocate ( bfach4 (NBTRG ) )
       allocate ( bfan2o (NBTRG ) )
       allocate ( dch4   (NBTRG ) )
       allocate ( dn2o   (NBTRG ) )
       allocate ( ech4   (NBTRG ) )
       allocate ( en2o   (NBTRG ) )
     endif

!----------------------------------------------------------------------
     if (Lw_control%do_ckd2p1) then
       NBLY = NBLY_CKD2P1
       call id2h2o ('INPUT/id2h2obdckd2p1')
     else
       NBLY = NBLY_ORIG
       call id2h2o ('INPUT/id2h2obdfull')
!------------------------------------------------------------------
!  read roberts continuum data for self-broadened h2o continuum
!-------------------------------------------------------------------
       call idrbtsh2o
     endif

     allocate  ( acomb(NBLY))
     allocate  ( bcomb(NBLY))
     allocate  ( apcm (NBLY))
     allocate  ( bpcm (NBLY))
     allocate  ( atpcm(NBLY))
     allocate  ( btpcm(NBLY))
     allocate  (bdlocm(NBLY))
     allocate  (bdhicm(NBLY))

     if (Lw_control%do_ckd2p1) then
       apwd = apwd_c
       bpwd = bpwd_c
       atpwd = atpwd_c
       btpwd = btpwd_c
       bdlowd = bdlowd_c
       bdhiwd = bdhiwd_c
       iband = iband_c
       acomb = acomb_c
       bcomb = bcomb_c
       apcm = apcm_c
       bpcm = bpcm_c
       atpcm = atpcm_c
       btpcm = btpcm_c
       bdlocm = bdlocm_c
       bdhicm = bdhicm_c
     else
       apwd = apwd_n
       bpwd = bpwd_n
       atpwd = atpwd_n
       btpwd = btpwd_n
       bdlowd = bdlowd_n
       bdhiwd = bdhiwd_n
       iband = iband_n
       acomb = acomb_n
       bcomb = bcomb_n
       apcm = apcm_n
       bpcm = bpcm_n
       atpcm = atpcm_n
       btpcm = btpcm_n
       bdlocm = bdlocm_n
       bdhicm = bdhicm_n
     endif

     call table_alloc (tab1 , NTTABH2O, NUTABH2O)
     call table_alloc (tab2 , NTTABH2O, NUTABH2O)
     call table_alloc (tab3 , NTTABH2O, NUTABH2O)
     call table_alloc (tab1w, NTTABH2O, NUTABH2O)
     if (Lw_control%do_ch4_n2o) then
       call table_alloc (tab1a, NTTABH2O, NUTABH2O, NBTRGE)
       call table_alloc (tab2a, NTTABH2O, NUTABH2O, NBTRGE)
       call table_alloc (tab3a, NTTABH2O, NUTABH2O, NBTRGE)
     endif
     call table_alloc (tabsr, NTTABH2O, NBLY    )

     call table

     call radiag_from_lwtables (bdlocm, bdhicm, iband, bandlo, bandhi)
!---------------------------------------------------------------------


end subroutine longwave_tables_init




!####################################################################

subroutine e1e290 (e1cts1, e1cts2, e1ctw1, e1ctw2, e1flx, emiss,  &
                   e1cts1f, e1cts2f, e1flxf, emissf)

!-----------------------------------------------------------------------
real, dimension (:,:,:), intent(out)   ::  e1cts1, e1cts2, e1ctw1, &
					   e1ctw2, e1flx, emiss
real, dimension (:,:,:, :), optional, intent(out) :: e1cts1f, e1cts2f, &
                     			             e1flxf, emissf
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
      real,    dimension(:,:,:),   allocatable :: avphilog, dt1, du, dup
      integer, dimension(:,:,:),   allocatable :: ixo1, iyo, iyop 
      real,    dimension(:,:,:),   allocatable :: avephi
      real,    dimension(:,:,:,:), allocatable :: avephif
      integer  :: k, m

!---------------------------------------------------------------------
!  allocate and retrieve variables from optical_path_mod
!---------------------------------------------------------------------
      allocate (avephi    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      if (Lw_control%do_ch4_n2o) then
        allocate (avephif (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
        call get_path_for_e90 (avephi, avephif)
      else
        call get_path_for_e90 (avephi)
      endif

!---------------------------------------------------------------------
!  allocatate local variables
!---------------------------------------------------------------------
      allocate (avphilog (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (dt1      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (du       (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (dup      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (ixo1     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (iyo      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (iyop     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )

!---------------------------------------------------------------------
!     obtain the "exchange" emissivities as a function of temperature 
!     (fxo) and water amount (fyo). the temperature indices have
!     been obtained in longwave_setup_mod.
!-------------------------------------------------------------------
      avphilog(:,:,KS:KE+1) = LOG10(avephi(:,:,KS:KE+1))
      call locate_in_table (mass_1, avphilog, du, iyo, KS, KE+1)
      call looktab (tab2, ixoe2, iyo, dte2, du, emiss, KS, KE+1)

!-----------------------------------------------------------------------
!     the special case emiss(:,:,KE) for layer KE is obtained by 
!     averaging the values for KE and KE+1.
!---------------------------------------------------------------------
      emiss(:,:,KE) = 0.5E+00*(emiss(:,:,KE) +emiss(:,:,KE+1))
 
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
      call looktab (tab1, ixo1, iyop, dt1, dup, e1flx, KS, KE+1)
      call looktab (tab1w, ixoe1, iyop, dte1, dup, e1ctw1, KS, KE+1)
      call looktab (tab1w, ixoe1, iyo, dte1, du, e1ctw2, KS, KE)

!--------------------------------------------------------------------
!     calculations with ch4 and n2o require NBTRGE separate emissivity
!     bands for h2o.
!--------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          avphilog(:,:,KS:KE+1) = LOG10(avephif(:,:,KS:KE+1,m))
          call locate_in_table (mass_1, avphilog, du, iyo, KS , KE+1)
          iyop(:,:,KS)        = 1
          iyop(:,:,KS+1:KE+1) = iyo(:,:,KS:KE)
          dup (:,:,KS)        = 0.0E+00
          dup (:,:,KS+1:KE+1) = du (:,:,KS:KE)
          call looktab (tab2a, ixoe2, iyo, dte2, du, emissf(:,:,:,m), &
		        KS, KE+1, m)
          call looktab (tab1a, ixoe1, iyop, dte1, dup,    &
		        e1cts1f(:,:,:,m), KS, KE+1, m)
          call looktab (tab1a, ixoe1, iyo, dte1, du, e1cts2f(:,:,:,m),&
		        KS, KE, m)
          call looktab (tab1a, ixo1, iyop, dt1, dup,   &
                        e1flxf(:,:,:,m), KS, KE+1, m)
        enddo

!--------------------------------------------------------------------
!     the special case emissf(:,:,KE,m) for layer KE is obtained by 
!     averaging the values for KE and KE+1.
!--------------------------------------------------------------------
        do m=1,NBTRGE
          emissf(:,:,KE,m) = 0.5E+00*    &
                             (emissf(:,:,KE,m) +emissf(:,:,KE+1,m))
        enddo
      endif

!--------------------------------------------------------------------
      deallocate (avephi)
      deallocate ( iyop       )
      deallocate ( iyo        )
      deallocate ( ixo1       )
      deallocate ( dup        )
      deallocate ( du         )
      deallocate ( dt1        )
      deallocate ( avphilog   )
      if (Lw_control%do_ch4_n2o) then
	deallocate (avephif)
      endif


end  subroutine e1e290




!###################################################################

subroutine e290 (k, emiss, emissb, emissf, emissbf)

!-----------------------------------------------------------------------
integer,                  intent(in)            ::  k
real, dimension(:,:,:),   intent(out)           :: emiss, emissb  
real, dimension(:,:,:,:), intent(out), optional :: emissf, emissbf
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
      real,    dimension(:,:,:),   allocatable :: avphilog, dtk, du
      real,    dimension(:,:,:),   allocatable :: avephi
      real,    dimension(:,:,:,:), allocatable :: avephif
      integer, dimension(:,:,:),   allocatable :: ixok, iyo          
      integer      :: kp, m

!--------------------------------------------------------------------
!   allocate and retrieve needed variables from optical_path_mod
!--------------------------------------------------------------------
      allocate (avephi    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      if (Lw_control%do_ch4_n2o) then
        allocate (avephif (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
        call get_path_for_e90 (avephi, avephif)
      else
        call get_path_for_e90 (avephi)
      endif

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
      allocate ( avphilog (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate ( dtk      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate ( du       (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate ( ixok     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate ( iyo      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )

!----------------------------------------------------------------------
      do kp=k,KE
        ixok(:,:,kp) = ixoe2(:,:,k-1)
        dtk (:,:,kp) = dte2 (:,:,k-1)
      end do

      avphilog(:,:,k:KE+1) = LOG10(avephi(:,:,k:KE+1))
      call locate_in_table (mass_1, avphilog, du, iyo,k, KE+1)
      call looktab (tab2, ixoe2, iyo, dte2, du, emiss, k, KE+1)
      call looktab (tab2, ixok, iyo, dtk, du, emissb, k, KE)

!--------------------------------------------------------------------
!     the special case emiss(:,:,KE) for layer KE is obtained by 
!     averaging the values for KE and KE+1. note that emiss(:,:,KE+1) 
!     is not useful after this point.
!-------------------------------------------------------------------
      emiss(:,:,KE) = 0.5E+00*(emiss(:,:,KE) +emiss(:,:,KE+1))
 
!--------------------------------------------------------------------
!     calculations with ch4 and n2o require NBTRGE separate emissivity
!     bands for h2o. reqults are in emissf (flux level k) and
!     emissbf (other levels).
!-------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          avphilog(:,:,k:KE+1) = LOG10(avephif(:,:,k:KE+1,m))
          call locate_in_table (mass_1, avphilog, du, iyo, k, KE+1)
          call looktab (tab2a, ixoe2, iyo, dte2, du, emissf(:,:,:,m), &
			k, KE+1, m)
          call looktab (tab2a, ixok, iyo, dtk, du, emissbf(:,:,:,m), &
			k, KE, m)
        enddo

!----------------------------------------------------------------------
!     the special case emissf(:,:,KE) for layer KE is obtained by 
!     averaging the values for KE and KE+1. note that emissf(:,:,KE+1,m)
!     is not useful after this point.
!----------------------------------------------------------------------
        do m=1,NBTRGE
          emissf(:,:,KE,m) = 0.5E+00*  &
                             (emissf(:,:,KE,m) +emissf(:,:,KE+1,m))
        enddo
      endif

!-------------------------------------------------------------------
      deallocate (avephi)
      deallocate (iyo        )
      deallocate (ixok       )
      deallocate (du         )
      deallocate (dtk        )
      deallocate (avphilog   )
      if (Lw_control%do_ch4_n2o) then
	deallocate (avephif)
      endif
!--------------------------------------------------------------------


end subroutine e290



!####################################################################

subroutine enear (emisdg, emspec, emisdgf, emspecf ) 
   
!--------------------------------------------------------------------
real, dimension (:,:,:),   intent(out)           ::   emisdg , emspec
real, dimension (:,:,:,:), intent(out), optional ::   emisdgf, emspecf
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
      real,    dimension(:,:,:),   allocatable :: dxsp, ylog, dysp, &
						  emiss,  emd1, emd2
      real,    dimension(:,:,:,:), allocatable :: emissf, emd2f, emd1f
      real,    dimension(:,:,:,:), allocatable :: empl1f, empl2f,  &
						  vrpfh2o
      real,    dimension(:,:,:  ), allocatable :: empl1, empl2, var2, &
						  emx1f, emx2f
      real,    dimension(:,:    ), allocatable :: emx1, emx2
      integer, dimension(:,:,:),   allocatable :: ixsp, iysp           

      integer   :: m

!----------------------------------------------------------------------
!  allocate and retrieve needed variables from optical_path_mod
!----------------------------------------------------------------------
      allocate (  empl1    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (  empl2    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (  var2     (ISRAD:IERAD, JSRAD:JERAD, KS:KE  )  )
      allocate (  emx2    (ISRAD:IERAD, JSRAD:JERAD                 ) )
      allocate (  emx1    (ISRAD:IERAD, JSRAD:JERAD                 ) )
      if (Lw_control%do_ch4_n2o) then
        allocate (  empl2f (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
        allocate (  empl1f (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
        allocate (  vrpfh2o(ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
        allocate (  emx2f  (ISRAD:IERAD, JSRAD:JERAD,          NBTRGE) )
        allocate (  emx1f  (ISRAD:IERAD, JSRAD:JERAD,          NBTRGE) )
        call get_path_for_enear (var2, emx1, emx2, empl1, empl2, &
			         vrpfh2o_out=vrpfh2o, emx1f_out=emx1f, &
			         emx2f_out=emx2f, empl1f_out=empl1f,  &
			         empl2f_out=empl2f )
      else
        call get_path_for_enear (var2, emx1, emx2, empl1, empl2)
      endif

!-------------------------------------------------------------------
!  allocate local variables
!-------------------------------------------------------------------
      allocate (  ixsp     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (  iysp     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (  dxsp     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (  ylog     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (  dysp     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (  emiss    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (  emd1     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (  emd2     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
      allocate (  emissf   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
      allocate (  emd2f    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
      allocate (  emd1f    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )

!--------------------------------------------------------------------
      ixsp(:,:,KE)   = ixoe2(:,:,KE-1)
      ixsp(:,:,KE+1) = ixoe1(:,:,KE-1)
      dxsp(:,:,KE)   = dte2(:,:,KE-1)
      dxsp(:,:,KE+1) = dte1(:,:,KE-1)

      ylog(:,:,KE  ) = ALOG10(var2(:,:,KE))
      ylog(:,:,KE+1) = ALOG10(var2(:,:,KE) + empl1(:,:,KE))

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
          ylog(:,:,KE  ) = ALOG10(vrpfh2o(:,:,KE,m))
          ylog(:,:,KE+1) = ALOG10(vrpfh2o(:,:,KE,m) + empl1f(:,:,KE,m))

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
      ylog(:,:,KS:KE+1) = ALOG10(empl1(:,:,KS:KE+1))
      call locate_in_table (mass_1, ylog, dysp, iysp, KS, KE+1)
      call looktab (tab3, ixsp, iysp, dxsp, dysp, emd1, KS, KE+1)

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          ylog(:,:,KS:KE+1) = ALOG10(empl1f(:,:,KS:KE+1,m))
          call locate_in_table (mass_1, ylog, dysp, iysp, KS, KE+1)
          call looktab (tab3a, ixsp, iysp, dxsp, dysp, emd1f(:,:,:,m), &
			KS, KE+1, m)
        enddo
      endif

!---------------------------------------------------------------------
      call locate_in_table (temp_1, tpl2, dxsp, ixsp, KS+1, KE+1)
      ylog(:,:,KS+1:KE+1) = ALOG10(empl2(:,:,KS+1:KE+1))
      call locate_in_table (mass_1, ylog, dysp, iysp, KS+1, KE+1)
      call looktab (tab3, ixsp, iysp, dxsp, dysp, emd2, KS+1, KE+1)

!----------------------------------------------------------------------
!     obtain index values of h2o pressure-scaled mass for each band
!     in the 1200-1400 range.
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
          ylog(:,:,KS+1:KE+1) = ALOG10(empl2f(:,:,KS+1:KE+1,m))
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
      emspec(:,:,KS     ) = (emd1(:,:,KS)*empl1(:,:,KS) -    &
                             emd1(:,:,KE+1)*empl1(:,:,KE+1))/  &
                             emx1(:,:) + 0.25E+00*(emiss(:,:,KE) +   &
                             emiss(:,:,KE+1))
      emspec(:,:,KS+1) = 2.0E+00*(emd1(:,:,KS)*empl1(:,:,KS) -    &
                         emd2(:,:,KE+1)*empl2(:,:,KE+1))/emx2(:,:)

     if (Lw_control%do_ch4_n2o) then
       do m=1,NBTRGE
         emisdgf(:,:,KS+1:KE,m) =    &
                            emd2f(:,:,KS+1:KE,m) + emd1f(:,:,KS+1:KE,m)
         emisdgf(:,:,KE+1,m) = 2.0E+00*emd1f(:,:,KE+1,m)
         emspecf(:,:,KS,m   ) = (emd1f(:,:,KS,m)*empl1f(:,:,KS,m) -   &
                              emd1f(:,:,KE+1,m)*empl1f(:,:,KE+1,m))/   &
                          emx1f(:,:,m) + 0.25E+00*(emissf(:,:,KE,m) +  &
                             emissf(:,:,KE+1,m))
         emspecf(:,:,KS+1,m) = 2.0E+00*    &
                                 (emd1f(:,:,KS,m)*empl1f(:,:,KS,m) -  &
                              emd2f(:,:,KE+1,m)*empl2f(:,:,KE+1,m)) / &
                              emx2f(:,:,m)
       enddo
     endif
!--------------------------------------------------------------------

     deallocate (var2    )
     deallocate (emx1    )
     deallocate (emx2    )
     deallocate (empl1   )
     deallocate (empl2   )

     if (Lw_control%do_ch4_n2o) then
       deallocate (vrpfh2o )
       deallocate (emx1f   )
       deallocate (emx2f   )
       deallocate (empl1f  )
       deallocate (empl2f  )
     endif

     deallocate (emd1f  )
     deallocate (emd2f  )
     deallocate (emissf )
     deallocate (emd2   )
     deallocate (emd1   )
     deallocate (emiss  )
     deallocate (dysp   )
     deallocate (ylog   )
     deallocate (dxsp   )
     deallocate (iysp   )
     deallocate (ixsp   )
!------------------------------------------------------------------



end subroutine enear
 




!#####################################################################

subroutine get_bds (bdlo_out, bdhi_out)

real, dimension(:), intent(out) :: bdlo_out, bdhi_out

      integer   :: ioffset

      ioffset = Lw_parameters%offset
      bdlo_out(:) = bdlocm(9+ioffset:8+NBCO215+ioffset) 
      bdhi_out(:) = bdhicm(9+ioffset:8+NBCO215+ioffset) 

end subroutine get_bds



!####################################################################

subroutine idrbtsh2o

!-----------------------------------------------------------------------
!     idrbtsh2o reads h2o roberts continuum quantities used in
!     longwave radiation
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/96
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------
      integer               :: inrad, k

!-----------------------------------------------------------------------
!  the following roberts continuum coefficients are computed using the
!  program (gasbnd) over the 0-3000 cm-1 range with 10 cm-1 bandwidth.
!-----------------------------------------------------------------------
      inrad = open_file ('INPUT/id2h2orbts', form='formatted', &
			 action='read')
      read (inrad, FMT = '(5e14.6)') (betad(k),k=1,NBLW)
      call close_file (inrad)
!---------------------------------------------------------------------
 
end subroutine idrbtsh2o


!#####################################################################

subroutine id2h2o (filename)

!---------------------------------------------------------------------
!     id2h2o reads h2o random model band parameters used for 
!     longwave radiation
!     references:
!
!     (1) fels, s. b., and m. d. schwarzkopf, "the simplified exchange
!         approximation: a new method for radiative transfer
!         calculations," journal atmospheric science, 32 (1975),
!         1475-1488.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/96
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------
character(len=*), intent(in)   :: filename

!----------------------------------------------------------------------
      integer                :: inrad, k
      real, dimension (NBLW) :: dummy

!-----------------------------------------------------------------------
!     the following h2o random band parameters are obtained from the
!     afgl 1992 HITRAN tape. parameters are obtained using an auxi-
!     liary program (gasbnd). values depend on assumptions as to
!     line shape, line strength and width. The inputted values span
!     the 0-3000 cm-1 range, with 10 cm-1 bandwidth. other parameter
!     values used in the program are obtained separately.
!-----------------------------------------------------------------------
      inrad = open_file (filename, form='formatted', action='read')

      read (inrad,9000) (arndm(k),k=1,NBLW)
      read (inrad,9000) (brndm(k),k=1,NBLW)
      read (inrad,9000) (dummy(k),k=1,NBLW)
      read (inrad,9000) (dummy(k),k=1,NBLW)
      read (inrad,9000) (dummy(k),k=1,NBLW)
      read (inrad,9000) (dummy(k),k=1,NBLW)
      read (inrad,9000) (bandlo(k),k=1,NBLW)
      read (inrad,9000) (bandhi(k),k=1,NBLW)

      call close_file (inrad)

9000  format(5e14.6)


end subroutine id2h2o




!#####################################################################

subroutine table 

!---------------------------------------------------------------------
!     table computes table entries used in longwave radiation.  
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!---------------------------------------------------------------------
 
      real, dimension (NBLW)               :: alfanb, anb, arotnb,   &
					      betanb, bnb, centnb, delnb
      real, dimension (30)                 :: cnusb, dnusb
      real, dimension (NTTABH2O,NBLW)      :: dbdtnb, src1nb
      real, dimension (NTTABH2O,NBLX)      :: srcwd        
      real, dimension (NTTABH2O, NUTABH2O) :: sumdbe, sum, sum3, sumwde
      real, dimension (NTTABH2O)           :: ddsc, fortcu, r1, r1wd, &
					      r2, s2, sc, srcs, sum4, &
					      sum4wd, sum6, sum7, sum8,&
					      t3, tfour, x, x1, xtemv
      real, dimension (NUTABH2O)           :: expo, fac, x2, zmass, &
					      zroot
      integer                              :: n, m, ioffset, itab,   &
					      jtab, nsubds, nsb,  iter
      real                                 :: zmassincr, cent, del,&
					      bdlo, bdhi, anu, c1,   &
					      freq_cutoff

      real, dimension (:,:), allocatable   :: r1a, r2a, s2a, t3a,   &
					      sum4a, sum6a, sum7a, sum8a
      real, dimension(:,:,:),allocatable   :: suma, sumdbea, sum3a 

!----------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
	allocate ( r1a     (NTTABH2O,NBTRGE) )
	allocate (r2a      (NTTABH2O,NBTRGE) )
	allocate (s2a      (NTTABH2O,NBTRGE) )
	allocate ( t3a     (NTTABH2O,NBTRGE) )
	allocate ( suma    (NTTABH2O,NUTABH2O,NBTRGE) )
	allocate ( sumdbea (NTTABH2O,NUTABH2O,NBTRGE) )
	allocate ( sum3a   (NTTABH2O,NUTABH2O,NBTRGE) )
	allocate ( sum4a   (NTTABH2O,NBTRGE) )
	allocate ( sum6a   (NTTABH2O,NBTRGE) )
	allocate ( sum7a   (NTTABH2O,NBTRGE) )
	allocate (sum8a    (NTTABH2O,NBTRGE) )
      endif

      ioffset = Lw_parameters%offset

!--------------------------------------------------------------------- 
!     compute local quantities and ao3, bo3, ab15 for narrow bands.
!---------------------------------------------------------------------
      do n=1,NBLW
        anb   (n) = arndm(n) 
        bnb   (n) = brndm(n) 
        centnb(n) = 0.5E+00*(bandlo(n) + bandhi(n)) 
        delnb (n) = bandhi(n) - bandlo(n)
        betanb(n) = betad(n)
      enddo

!---------------------------------------------------------------------
!    compute a*b and sqrt(a*b) for all 10 cm-1 frequency bands.
!---------------------------------------------------------------------
      do n=1,NBLW
        alfanb(n) = bnb(n)*anb(n) 
        arotnb(n) = SQRT(alfanb(n))
      enddo

   if (Lw_control%do_ch4_n2o) then
 
!---------------------------------------------------------------------
!     select, from random band coefficients in 1200-1400 cm-1 range,
!     those appropriate for NBTRGE h2o bands.
!---------------------------------------------------------------------
      if (NBTRGE .EQ. 1) then
        do m=1,NBTRGE
          afah2o(m) = CN_basic%af1h2o
          bfah2o(m) = CN_basic%bf1h2o
          bdlah2o(m) = CN_basic%bdl1h2o
          bdhah2o(m) = CN_basic%bdh1h2o
        end do
      elseif (NBTRGE .EQ. 2) then
        do  m=1,NBTRGE
          afah2o(m) = CN_basic%af2h2o(m)
          bfah2o(m) = CN_basic%bf2h2o(m)
          bdlah2o(m) = CN_basic%bdl2h2o(m)
          bdhah2o(m) = CN_basic%bdh2h2o(m)
        end do
      elseif (NBTRGE .EQ. 4) then
        do  m=1,NBTRGE
          afah2o(m) = CN_basic%af4h2o(m)
          bfah2o(m) = CN_basic%bf4h2o(m)
          bdlah2o(m) = CN_basic%bdl4h2o(m)
          bdhah2o(m) = CN_basic%bdh4h2o(m)
        end do
      elseif (NBTRGE .EQ. 10) then
        do m=1,NBTRGE
          afah2o(m) = CN_basic%af10h2o(m)
          bfah2o(m) = CN_basic%bf10h2o(m)
          bdlah2o(m) = CN_basic%bdl10h2o(m)
          bdhah2o(m) = CN_basic%bdh10h2o(m)
        enddo
      elseif (NBTRGE .EQ. 20) then
        do m=1,NBTRGE
          afah2o(m) = arndm(m+120)
          bfah2o(m) = brndm(m+120)
          bdlah2o(m) = bandlo(m+120)
          bdhah2o(m) = bandhi(m+120)
        enddo
      else
        call error_mesg ( 'longwave_tables',   &
	   'NBTRGE is inconsistent with available data', FATAL)
      endif
!-------------------------------------------------------------------
!   define critical frequency (cutoff for wide band ?? )
!------------------------------------------------------------------
      freq_cutoff = 1400.

  else

      freq_cutoff = 1200.

  endif

!---------------------------------------------------------------------
!     begin table computations here.  compute temperatures and masses
!     for table entries.
! 
!     note: the dimensioning and initialization of xtemv and other
!     arrays with dimension of NTTABH2O=28 imply a restriction of model 
!     temperatures from 100k to 370k.
! 
!     the dimensioning of zmass, zroot and other arrays with 
!     dimension of NUTABH2O=180 imply a restriction of model h2o amounts
!     such that optical paths are between 10**-16 and 10**2, in cgs 
!     units.
!---------------------------------------------------------------------
      zmass(1) = 10.0E+00**mass_1%min_val
      zmassincr = 10.0E+00**mass_1%tab_inc
 
!---------------------------------------------------------------------
!   the definition of zmassincr as 10**0.1 is slightly different from
!   all previous versions, in which it is 1.258925411E+00. This
!   produces slightly different answers (fluxes differ by 1.0e-6 W/m2).
!---------------------------------------------------------------------
      do jtab=2,NUTABH2O
        zmass(jtab) = zmass(jtab-1)*zmassincr
      enddo
      do jtab=1,NUTABH2O
        zroot(jtab) = SQRT(zmass(jtab))
      enddo 
      do itab=1,NTTABH2O
        xtemv (itab) = temp_1%min_val + temp_1%tab_inc*(itab-1)
        tfour (itab) = xtemv(itab)**4
        fortcu(itab) = 4.0E+00*xtemv(itab)**3
      enddo
      
!---------------------------------------------------------------------
!     the computation of source, dsrce is needed only for the combined 
!     wide band case.  to obtain them,  the source must be computed 
!     for each of the NBLX wide bands srcwd then combined using iband
!     into source.
!---------------------------------------------------------------------
      do n=1,NBLY
        do itab=1,NTTABH2O
          tabsr%vae  (itab,n) = 0.0E+00
        enddo
      enddo
      do n=1,NBLX
        do itab=1,NTTABH2O
          srcwd(itab,n) = 0.0E+00
        enddo
      enddo

!---------------------------------------------------------------------
!     begin frequency loop.
!---------------------------------------------------------------------
      do n=1,NBLX 
  
!---------------------------------------------------------------------
!     the 160-560 cm-1 region
!---------------------------------------------------------------------
        if (n .LE. 40) then
          cent = centnb(n+16) 
          del  = delnb (n+16) 
          bdlo = bandlo(n+16) 
          bdhi = bandhi(n+16) 
 
!---------------------------------------------------------------------
!      the 560-1200 cm-1 region, and the 2270-2380 cm-1 region
!---------------------------------------------------------------------
        else
          cent = 0.5E+00*(bdlocm(n-32+ioffset) + bdhicm(n-32+ioffset))
          del  = bdhicm(n-32+ioffset) - bdlocm(n-32+ioffset)
          bdlo = bdlocm(n-32+ioffset)
          bdhi = bdhicm(n-32+ioffset)
        endif

!---------------------------------------------------------------------
!     for purposes of accuracy, all evaluations of planck functions
!     are made on 10 cm-1 intervals, then summed into the NBLX wide 
!     bands.  the last subband may be narrower than 10 cm-1.
!---------------------------------------------------------------------
        nsubds = (del - 1.0E-03)/10 + 1
        do nsb=1,nsubds 
          if(nsb .NE. nsubds) then 
            cnusb(nsb) = 10.0E+00*(nsb - 1) + bdlo + 5.0E+00
            dnusb(nsb) = 10.0E+00
          else
            cnusb(nsb) = 0.5E+00*(10.0E+00*(nsb - 1) + bdlo + bdhi)
            dnusb(nsb) = bdhi -  (10.0E+00*(nsb - 1) + bdlo)
          endif 
          c1 = 3.7412E-05*cnusb(nsb)**3

!---------------------------------------------------------------------
!     begin temperature loop.
!---------------------------------------------------------------------
          do itab=1,NTTABH2O
            x    (itab)   = 1.4387E+00*cnusb(nsb)/xtemv(itab)
            x1   (itab)   = EXP(x(itab)) 
            srcs (itab)   = c1/(x1(itab) - 1.0E+00)
            srcwd(itab,n) = srcwd(itab,n) + srcs(itab)*dnusb(nsb)
          enddo
        enddo
      enddo

!---------------------------------------------------------------------
!     the following loops create the combined wide band quantities 
!#ifndef ckd2p1
!     source and dsrce.  the first 40 bands map to bands 1 to 8 in
!#else   ckd2p1
!     source and dsrce.  the first 40 bands map to bands 1 to 40 in
!#endif ckd2p1
!     source and dsrce.
!---------------------------------------------------------------------
      do n=1,40
        do itab=1,NTTABH2O 
          tabsr%vae  (itab,iband(n)) = tabsr%vae(itab,iband(n)) +      &
                                  srcwd(itab,n)
        enddo
      enddo
      do n=9+ioffset,NBLY  
        do itab=1,NTTABH2O
          tabsr%vae  (itab,n) = srcwd(itab,n+32-ioffset)
        enddo
      enddo
      do n=1,NBLY
        do itab=1,NTTABH2O-1 
          tabsr%td(itab,n) = (tabsr%vae(itab+1,n) -   &
	                      tabsr%vae(itab,n))*0.1E+00
        enddo
      enddo

!---------------------------------------------------------------------
!     first compute planck functions src1nb and derivatives dbdtnb 
!     for use in table evaluations.  these are different from source,
!     dsrce because different frequency points are used in evaluation,
!     the frequency ranges are different, and the derivative algorithm
!     is different.
!---------------------------------------------------------------------
      do n=1,NBLW 
        cent = centnb(n)
        del  = delnb (n)

!---------------------------------------------------------------------
!     note: at present, the iter loop is only used for iter=2.  the 
!     loop structure is kept so that in the future, we may use a
!     quadrature scheme for the planck function evaluation, rather
!     than use the mid-band frequency.
!---------------------------------------------------------------------
        do iter=2,2
          anu = cent + 0.5E+00*(iter - 2)*del 
          c1  = (3.7412E-05)*anu**3
!---------------------------------------------------------------------
!     temperature loop.
!---------------------------------------------------------------------
          do itab=1,NTTABH2O
            x  (itab)      = 1.4387E+00*anu/xtemv(itab)
            x1 (itab)      = EXP(x(itab))
            sc (itab)      = c1/((x1(itab) - 1.0E+00) + 1.0E-20) 
            sc (itab)      = c1/(x1(itab) - 1.0E+00)
            ddsc(itab)     = sc(itab)/(x1(itab)-1.0E+00)*x1(itab)*  &
                             x(itab)/xtemv(itab)
            src1nb(itab,n) = del*sc (itab)
            dbdtnb(itab,n) = del*ddsc(itab)
          enddo
        enddo
      enddo

!---------------------------------------------------------------------
!     next compute r1, r2, s2, and t3 coefficients used for e3 
!     function when the optical path is less than 10**-4.  in this 
!     case, we assume a different dependence on zmass.  also obtain 
!     r1wd, which is r1 summed over the 160-560 cm-1 range.
!---------------------------------------------------------------------
      do itab=1,NTTABH2O
        sum4  (itab) = 0.0E+00
        sum6  (itab) = 0.0E+00
        sum7  (itab) = 0.0E+00
        sum8  (itab) = 0.0E+00
        sum4wd(itab) = 0.0E+00
      enddo

      if (Lw_control%do_ch4_n2o) then
        sum4a (:,:)    = 0.0E+00
        sum6a (:,:)    = 0.0E+00
        sum7a (:,:)    = 0.0E+00
        sum8a (:,:)    = 0.0E+00
      endif

      do n=1,NBLW 
        cent = centnb(n)
!---------------------------------------------------------------------
!#ifndef ch4n2o
!     perform summations for frequency ranges of 0-560, 1200-2200 cm-1 
!#else   ch4n2o
!     perform summations for frequency ranges of 0-560, 1400-2200 cm-1 
!#endif ch4n2o
!---------------------------------------------------------------------
        if (cent .LT. 5.6E+02 .OR. cent .GT. freq_cutoff .AND.   &
            cent .LE. 2.2E+03) then
          do itab=1,NTTABH2O 
            sum4(itab) = sum4(itab) + src1nb(itab,n)
            sum6(itab) = sum6(itab) + dbdtnb(itab,n)
            sum7(itab) = sum7(itab) + dbdtnb(itab,n)*arotnb(n)
            sum8(itab) = sum8(itab) + dbdtnb(itab,n)*alfanb(n)
          enddo
        endif

       if (Lw_control%do_ch4_n2o) then

!---------------------------------------------------------------------
!     perform summations for frequency range of 1200-1400 cm-1
!     for sum4a, sum6a, sum7a, and sum8a. the computation depends
!     on the value of NBTRGE.
!---------------------------------------------------------------------
        if (cent .GT. 1.2E+03 .AND. cent .LE. 1.4E+03) then
          do m=1,NBTRGE
            if (cent .GT. bdlah2o(m) .AND. cent .LE. bdhah2o(m)) then
              sum4a(:,m) = sum4a(:,m) + src1nb(:,n)
              sum6a(:,m) = sum6a(:,m) + dbdtnb(:,n)
              sum7a(:,m) = sum7a(:,m) + dbdtnb(:,n)*arotnb(n)
              sum8a(:,m) = sum8a(:,m) + dbdtnb(:,n)*alfanb(n)
            endif
          enddo
        endif
       endif

!---------------------------------------------------------------------
!     perform summations over 160-560 cm-1 frequency range for e1 
!     calculations sum4wd.
!---------------------------------------------------------------------
        if(cent .GT. 1.6E+02 .AND. cent .LT. 5.6E+02) then
          do itab=1,NTTABH2O 
            sum4wd(itab) = sum4wd(itab) + src1nb(itab,n)
          enddo
        endif 
      enddo
      do itab=1,NTTABH2O
        r1(itab)   = sum4(itab)/tfour (itab)
        r2(itab)   = sum6(itab)/fortcu(itab) 
        s2(itab)   = sum7(itab)/fortcu(itab) 
        t3(itab)   = sum8(itab)/fortcu(itab) 
        r1wd(itab) = sum4wd(itab)/tfour(itab)
      enddo
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          sum   (itab,jtab) = 0.0E+00 
          sumdbe(itab,jtab) = 0.0E+00
          sum3  (itab,jtab) = 0.0E+00
          sumwde(itab,jtab) = 0.0E+00
        enddo
      enddo
      if (Lw_control%do_ch4_n2o) then
      do m=1,NBTRGE
        r1a(:,m)   = sum4a(:,m)/tfour(:)
        r2a(:,m)   = sum6a(:,m)/fortcu(:)
        s2a(:,m)   = sum7a(:,m)/fortcu(:)
        t3a(:,m)   = sum8a(:,m)/fortcu(:)
      enddo
        suma   (:,:,:) = 0.0E+00 
        sumdbea(:,:,:) = 0.0E+00
        sum3a  (:,:,:) = 0.0E+00
      endif

!---------------------------------------------------------------------
!     frequency loop begins.
!---------------------------------------------------------------------
      do n=1,NBLW 
        cent = centnb(n)

!---------------------------------------------------------------------
!     perform calculations for frequency ranges of 0-560, 
!#ifndef ch4n2o
!     1200-2200 cm-1.
!#else   ch4n2o
!     1400-2200 cm-1.
!#endif  ch4n2o
!---------------------------------------------------------------------
        if(cent .LT. 5.6E+02 .OR. cent .GT. freq_cutoff .AND.    &
           cent .LE. 2.2E+03) then
          do jtab=1,NUTABH2O 
            x2  (jtab) = arotnb(n)*zroot(jtab) 
            expo(jtab) = EXP( - x2(jtab))
          enddo
	  do jtab=121,NUTABH2O
            fac(jtab) = (1.0E+00 - (1.0E+00 + x2(jtab))*expo(jtab))/ &
                        (alfanb(n)*zmass(jtab))
          enddo
          do jtab=1,NUTABH2O 
            do itab=1,NTTABH2O
              sum   (itab,jtab) = sum   (itab,jtab) +   &
                                  src1nb(itab,n)*expo(jtab)
              sumdbe(itab,jtab) = sumdbe(itab,jtab) +    &
                                  dbdtnb(itab,n)*expo(jtab)
            enddo
          enddo 
          do jtab=121,NUTABH2O
            do itab=1,NTTABH2O 
              sum3(itab,jtab) = sum3(itab,jtab) +    &
                                dbdtnb(itab,n)*fac(jtab)
            enddo 
          enddo 
        endif

!-------------------------------------------------------------------
!   perform calculations over the frequency range 1200-1400 cm-1. 
!   the calculations depend on the value of NBTRGE.
!-------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        if (cent .GT. 1.2E+03 .AND. cent .LE. 1.4E+03) then
          do m=1,NBTRGE
            if (cent .GT. bdlah2o(m) .AND. cent .LE. bdhah2o(m)) then
              x2  (:) = arotnb(n)*zroot(:) 
              expo(:) = EXP( - x2(:))
              do jtab=121,NUTABH2O
          fac(jtab) = (1.0E+00 - (1.0E+00 + x2(jtab))*expo(jtab))/   &
                          (alfanb(n)*zmass(jtab))
	      enddo
              do jtab=1,NUTABH2O 
                suma(:,jtab,m)    = suma(:,jtab,m) +  &
                                    src1nb(:,n)*expo(jtab)
                sumdbea(:,jtab,m) = sumdbea(:,jtab,m) +   &
                                    dbdtnb(:,n)*expo(jtab)
	      enddo
              do jtab=121,NUTABH2O 
                sum3a(:,jtab,m)   = sum3a(:,jtab,m) +   &
                                    dbdtnb(:,n)*fac(jtab)
	      enddo
            endif
          enddo
        endif
      endif

!---------------------------------------------------------------------
!     compute sum over 160-560 cm-1 range for use in e1 calculations
!     sumwde.
!---------------------------------------------------------------------
        if(cent .GT. 1.6E+02 .AND. cent .LT. 5.6E+02) then 
          do jtab=1,NUTABH2O
            do itab=1,NTTABH2O 
              sumwde(itab,jtab) = sumwde(itab,jtab) +     &
                                  src1nb(itab,n)*expo(jtab)
            enddo
          enddo 
        endif 
      enddo

!--------------------------------------------------------------------
!     frequency loop ends
!--------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          tab1%vae      (itab,jtab) = sum(itab,jtab)/tfour(itab)
          tab2%vae(itab,jtab) = sumdbe(itab,jtab)/fortcu(itab)
        enddo 
      enddo
      do jtab=121,NUTABH2O
        do itab=1,NTTABH2O
          tab3%vae(itab,jtab) = sum3(itab,jtab)/fortcu(itab)
        enddo
      enddo
      do jtab=1,2
        do itab=1,NTTABH2O
          tab1%vae      (itab,jtab) = r1(itab)
        enddo
      enddo
      do jtab=1,120
        do itab=1,NTTABH2O
          tab3%vae(itab,jtab) = r2(itab)/2.0E+00 -    &
                                s2(itab)*zroot(jtab)/3.0E+00 +   &
                                t3(itab)*zmass(jtab)/8.0E+00
        enddo
      enddo

!---------------------------------------------------------------------
!     compute e1 tables for 160-560 cm-1 bands.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          tab1w%vae      (itab,jtab) = sumwde(itab,jtab)/tfour(itab)
        enddo
      enddo
      do jtab=1,2
        do itab=1,NTTABH2O
          tab1w%vae      (itab,jtab) = r1wd(itab)
        enddo 
      enddo

!---------------------------------------------------------------------
!     initialize all derivative table entries.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O
          tab1%td (itab,jtab) = 0.0E+00
          tab1w%td(itab,jtab) = 0.0E+00
          tab2%td  (itab,jtab) = 0.0E+00
          tab3%td (itab,jtab) = 0.0E+00
          tab1%md (itab,jtab) = 0.0E+00
          tab1w%md(itab,jtab) = 0.0E+00
          tab2%md  (itab,jtab) = 0.0E+00
          tab3%md (itab,jtab) = 0.0E+00
          tab1%cd (itab,jtab) = 0.0E+00
          tab1w%cd(itab,jtab) = 0.0E+00
          tab2%cd  (itab,jtab) = 0.0E+00
          tab3%cd (itab,jtab) = 0.0E+00
        enddo
      enddo

!---------------------------------------------------------------------
!     compute table entries for temperature derivatives.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O
        do itab=1,NTTABH2O-1

           tab1%td  (itab,jtab) =    &
          (tab1%vae(itab+1,jtab) - tab1%vae (itab,jtab))/temp_1%tab_inc

          tab1w%td(itab,jtab) =     &
         (tab1w%vae(itab+1,jtab) - tab1w%vae(itab,jtab))/temp_1%tab_inc

          tab2%td  (itab,jtab) =    &
       (tab2%vae  (itab+1,jtab) - tab2%vae  (itab,jtab))/temp_1%tab_inc

          tab3%td (itab,jtab) =     &
         (tab3%vae (itab+1,jtab) - tab3%vae (itab,jtab))/temp_1%tab_inc

        enddo
      enddo

!---------------------------------------------------------------------
!     compute table entries for mass derivatives.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O-1
        do itab=1,NTTABH2O

          tab1%md (itab,jtab) =   &
     (tab1%vae (itab,jtab+1) - tab1%vae (itab,jtab))/mass_1%tab_inc

          tab1w%md(itab,jtab) =   &
    (tab1w%vae(itab,jtab+1) - tab1w%vae(itab,jtab))/mass_1%tab_inc

          tab2%md  (itab,jtab) =    &
   (tab2%vae  (itab,jtab+1) - tab2%vae  (itab,jtab))/mass_1%tab_inc

          tab3%md (itab,jtab) =    &
   (tab3%vae (itab,jtab+1) - tab3%vae (itab,jtab))/mass_1%tab_inc

        enddo
      enddo

!---------------------------------------------------------------------
!     compute table entries for cross derivatives.
!---------------------------------------------------------------------
      do jtab=1,NUTABH2O-1
        do itab=1,NTTABH2O-1

      tab1%cd (itab,jtab) =    &
             (tab1%vae (itab+1,jtab+1) - tab1%vae (itab+1,jtab) -   &
              tab1%vae (itab  ,jtab+1) + tab1%vae (itab  ,jtab))/   &
             (temp_1%tab_inc*mass_1%tab_inc)

      tab1w%cd(itab,jtab) =    &
             (tab1w%vae(itab+1,jtab+1) - tab1w%vae(itab+1,jtab) -   &
              tab1w%vae(itab  ,jtab+1) + tab1w%vae(itab  ,jtab))/   &
             (temp_1%tab_inc*mass_1%tab_inc)


!!  THIS NEVER USED :
!     tab2%cd  (itab,jtab) =     &
!            (tab2%vae  (itab+1,jtab+1) - tab2%vae  (itab+1,jtab) -   &
!             tab2%vae  (itab  ,jtab+1) + tab2%vae  (itab  ,jtab))/   &
!            (DTTABH2O*DUTABH2O)
!            (temp_1%tab_inc*mass_1%tab_inc)


      tab3%cd (itab,jtab) =     &
             (tab3%vae (itab+1,jtab+1) - tab3%vae (itab+1,jtab) -   &
              tab3%vae (itab  ,jtab+1) + tab3%vae (itab  ,jtab))/   &
             (temp_1%tab_inc*mass_1%tab_inc)

        enddo
      enddo
      if (Lw_control%do_ch4_n2o) then
      do m=1,NBTRGE
        do jtab=1,NUTABH2O
          tab1a%vae      (:,jtab,m) = suma(:,jtab,m)/tfour(:)
          tab2a%vae (:,jtab,m) = sumdbea(:,jtab,m)/fortcu(:) 
        enddo
        do jtab=121,NUTABH2O
          tab3a%vae(:,jtab,m) = sum3a(:,jtab,m)/fortcu(:)
        enddo
        do jtab=1,2
          tab1a%vae      (:,jtab,m) = r1a(:,m)
        enddo
        do jtab=1,120
          tab3a%vae(:,jtab,m) = r2a(:,m)/2.0E+00 -     &
                                s2a(:,m)*zroot(jtab)/3.0E+00 +   &
                                t3a(:,m)*zmass(jtab)/8.0E+00
        enddo
      enddo

      tab1a%td (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
      tab2a%td  (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
      tab3a%td (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
      tab1a%md (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
      tab2a%md  (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
      tab3a%md (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
      tab1a%cd (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
      tab2a%cd (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00
      tab3a%cd (1:NTTABH2O,1:NUTABH2O,:) = 0.0E+00

      tab1a%td(1:NTTABH2O-1,1:NUTABH2O,:) =    &
                 (tab1a%vae(2:NTTABH2O,1:NUTABH2O,:) -   &
                  tab1a%vae(1:NTTABH2O-1,1:NUTABH2O,:))/temp_1%tab_inc

      tab2a%td (1:NTTABH2O-1,1:NUTABH2O,:) =   &
                (tab2a%vae (2:NTTABH2O,1:NUTABH2O,:) -    &
                 tab2a%vae (1:NTTABH2O-1,1:NUTABH2O,:))/temp_1%tab_inc

      tab3a%td(1:NTTABH2O-1,1:NUTABH2O,:) =    &
                  (tab3a%vae(2:NTTABH2O,1:NUTABH2O,:) -  &
                   tab3a%vae(1:NTTABH2O-1,1:NUTABH2O,:))/temp_1%tab_inc

      tab1a%md(1:NTTABH2O,1:NUTABH2O-1,:) =     &
                  (tab1a%vae(1:NTTABH2O,2:NUTABH2O,:) -   &
                   tab1a%vae(1:NTTABH2O,1:NUTABH2O-1,:))/mass_1%tab_inc

      tab2a%md (1:NTTABH2O,1:NUTABH2O-1,:) =   &
                 (tab2a%vae (1:NTTABH2O,2:NUTABH2O,:) -   &
                  tab2a%vae (1:NTTABH2O,1:NUTABH2O-1,:))/mass_1%tab_inc

      tab3a%md(1:NTTABH2O,1:NUTABH2O-1,:) =     &
                  (tab3a%vae(1:NTTABH2O,2:NUTABH2O,:) -   &
                   tab3a%vae(1:NTTABH2O,1:NUTABH2O-1,:))/mass_1%tab_inc

      tab1a%cd(1:NTTABH2O-1,1:NUTABH2O-1,:) =     &
                        (tab1a%vae(2:NTTABH2O,2:NUTABH2O,:) -    &
                         tab1a%vae(2:NTTABH2O,1:NUTABH2O-1,:)   -  &
                         tab1a%vae(1:NTTABH2O-1,2:NUTABH2O,:)   +  &
                         tab1a%vae(1:NTTABH2O-1,1:NUTABH2O-1,:))/  &
                                         (temp_1%tab_inc*mass_1%tab_inc)

      tab3a%cd(1:NTTABH2O-1,1:NUTABH2O-1,:) =    &
                        (tab3a%vae(2:NTTABH2O,2:NUTABH2O,:) -    &
                         tab3a%vae(2:NTTABH2O,1:NUTABH2O-1,:)   -  &
                         tab3a%vae(1:NTTABH2O-1,2:NUTABH2O,:)   +  &
                         tab3a%vae(1:NTTABH2O-1,1:NUTABH2O-1,:))/  &
                                        (temp_1%tab_inc*mass_1%tab_inc)

	deallocate ( r1a    )
	deallocate (r2a     )
	deallocate (s2a     )
	deallocate ( t3a     )
	deallocate ( suma    )
	deallocate ( sumdbea )
	deallocate ( sum3a   )
	deallocate ( sum4a   )
	deallocate ( sum6a   )
	deallocate ( sum7a   )
	deallocate (sum8a  )
      endif



end subroutine table

 
	      end module longwave_tables_mod
		
