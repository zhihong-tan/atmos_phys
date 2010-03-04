 
!---------------------------------------------------------------------
!------------ FMS version number and tagname for this file -----------
        
! $Id: optics_lib.f90,v 18.0 2010/03/02 23:29:38 fms Exp $
! $Name: riga $
 
! OPTICS_LIB: Optical proecures for for F90
! Compiled/Modified:
!   07/01/06  John Haynes (haynes@atmos.colostate.edu)
!
! m_wat (subroutine)
! m_ice (subroutine)
! mie_int (subroutine)
  
  module optics_lib
  implicit none

  contains

! ----------------------------------------------------------------------------
! subroutine M_WAT
! ----------------------------------------------------------------------------
  subroutine m_wat(freq, t, n_r, n_i)
  implicit none
!  
! Purpose:
!   compute complex index of refraction of liquid water
!
! Inputs:
!   [freq]    frequency (GHz)
!   [t]       temperature (C)
!
! Outputs:
!   [n_r]     real part index of refraction
!   [n_i]     imaginary part index of refraction
!
! Reference:
!   Based on the work of Ray (1972)
!
! Coded:
!   03/22/05  John Haynes (haynes@atmos.colostate.edu)
  
! ----- INPUTS -----
  real*8, intent(in) :: freq,t
  
! ----- OUTPUTS -----
  real*8, intent(out) :: n_r, n_i

! ----- INTERNAL -----    
  real*8 ld,es,ei,a,ls,sg,tm1,cos1,sin1
  real*8 e_r,e_i
  real*8 pi
  complex*16 e_comp, sq

  ld = 100.*2.99792458E8/(freq*1E9)
  es = 78.54*(1-(4.579E-3*(t-25.)+1.19E-5*(t-25.)**2 &
       -2.8E-8*(t-25.)**3))
  ei = 5.27137+0.021647*t-0.00131198*t**2
  a = -(16.8129/(t+273.))+0.0609265
  ls = 0.00033836*exp(2513.98/(t+273.))
  sg = 12.5664E8

  tm1 = (ls/ld)**(1-a)
  pi = acos(-1.D0)
  cos1 = cos(0.5*a*pi)
  sin1 = sin(0.5*a*pi)

  e_r = ei + (((es-ei)*(1.+tm1*sin1))/(1.+2*tm1*sin1+tm1**2))
  e_i = (((es-ei)*tm1*cos1)/(1.+2*tm1*sin1+tm1**2)) &
        +((sg*ld)/1.885E11)

  e_comp = dcmplx(e_r,e_i)
  sq = sqrt(e_comp)
  
  n_r = real(sq)
  n_i = aimag(sq)      
  
  return
  end subroutine m_wat

! ----------------------------------------------------------------------------
! subroutine M_ICE
! ----------------------------------------------------------------------------
  subroutine m_ice(freq,t,n_r,n_i)
  implicit none
!
! Purpose:
!   compute complex index of refraction of ice
!
! Inputs:
!   [freq]    frequency (GHz)
!   [t]       temperature (C)
!
! Outputs:
!   [n_r]     real part index of refraction
!   [n_i]     imaginary part index of refraction
!
! Reference:
!    Fortran 90 port from IDL of REFICE by Stephen G. Warren
!
! Modified:
!   05/31/05  John Haynes (haynes@atmos.colostate.edu)

! ----- INPUTS -----
  real*8, intent(in) :: freq, t
  
! ----- OUTPUTS -----  
  real*8, intent(out) :: n_r,n_i

! Parameters:
  integer*2 :: i,lt1,lt2,nwl,nwlt
  parameter(nwl=468,nwlt=62)

  real*8 :: alam,cutice,pi,t1,t2,tk,wlmax,wlmin, &
            x,x1,x2,y,y1,y2,ylo,yhi

  real*8 :: &
       tabim(nwl),tabimt(nwlt,4),tabre(nwl),tabret(nwlt,4),temref(4), &
       wl(nwl),wlt(nwlt)

! Defines wavelength dependent complex index of refraction for ice.
! Allowable wavelength range extends from 0.045 microns to 8.6 meter
! temperature dependence only considered beyond 167 microns.
! 
! interpolation is done     n_r  vs. log(xlam)
!                           n_r  vs.        t
!                       log(n_i) vs. log(xlam)
!                       log(n_i) vs.        t
!
! Stephen G. Warren - 1983
! Dept. of Atmospheric Sciences
! University of Washington
! Seattle, Wa  98195
!
! Based on
!
!    Warren,S.G.,1984.
!    Optical constants of ice from the ultraviolet to the microwave.
!    Applied Optics,23,1206-1225
!
! Reference temperatures are -1.0,-5.0,-20.0, and -60.0 deg C
 
      data temref/272.16,268.16,253.16,213.16/
 
      data wlmin,wlmax/0.045,8.6e6/
      data cutice/167.0/
 
      data (wl(i),i=1,114)/ &
      0.4430e-01,0.4510e-01,0.4590e-01,0.4680e-01,0.4770e-01,0.4860e-01, &
      0.4960e-01,0.5060e-01,0.5170e-01,0.5280e-01,0.5390e-01,0.5510e-01, &
      0.5640e-01,0.5770e-01,0.5900e-01,0.6050e-01,0.6200e-01,0.6360e-01, &
      0.6530e-01,0.6700e-01,0.6890e-01,0.7080e-01,0.7290e-01,0.7380e-01, &
      0.7510e-01,0.7750e-01,0.8000e-01,0.8270e-01,0.8550e-01,0.8860e-01, &
      0.9180e-01,0.9300e-01,0.9540e-01,0.9920e-01,0.1033e+00,0.1078e+00, &
      0.1100e+00,0.1127e+00,0.1140e+00,0.1181e+00,0.1210e+00,0.1240e+00, &
      0.1272e+00,0.1295e+00,0.1305e+00,0.1319e+00,0.1333e+00,0.1348e+00, &
      0.1362e+00,0.1370e+00,0.1378e+00,0.1387e+00,0.1393e+00,0.1409e+00, &
      0.1425e+00,0.1435e+00,0.1442e+00,0.1450e+00,0.1459e+00,0.1468e+00, &
      0.1476e+00,0.1480e+00,0.1485e+00,0.1494e+00,0.1512e+00,0.1531e+00, &
      0.1540e+00,0.1550e+00,0.1569e+00,0.1580e+00,0.1589e+00,0.1610e+00, &
      0.1625e+00,0.1648e+00,0.1669e+00,0.1692e+00,0.1713e+00,0.1737e+00, &
      0.1757e+00,0.1779e+00,0.1802e+00,0.1809e+00,0.1821e+00,0.1833e+00, &
      0.1843e+00,0.1850e+00,0.1860e+00,0.1870e+00,0.1880e+00,0.1890e+00, &
      0.1900e+00,0.1910e+00,0.1930e+00,0.1950e+00,0.2100e+00,0.2500e+00, &
      0.3000e+00,0.3500e+00,0.4000e+00,0.4100e+00,0.4200e+00,0.4300e+00, &
      0.4400e+00,0.4500e+00,0.4600e+00,0.4700e+00,0.4800e+00,0.4900e+00, &
      0.5000e+00,0.5100e+00,0.5200e+00,0.5300e+00,0.5400e+00,0.5500e+00/
      data (wl(i),i=115,228)/ &
      0.5600e+00,0.5700e+00,0.5800e+00,0.5900e+00,0.6000e+00,0.6100e+00, &
      0.6200e+00,0.6300e+00,0.6400e+00,0.6500e+00,0.6600e+00,0.6700e+00, &
      0.6800e+00,0.6900e+00,0.7000e+00,0.7100e+00,0.7200e+00,0.7300e+00, &
      0.7400e+00,0.7500e+00,0.7600e+00,0.7700e+00,0.7800e+00,0.7900e+00, &
      0.8000e+00,0.8100e+00,0.8200e+00,0.8300e+00,0.8400e+00,0.8500e+00, &
      0.8600e+00,0.8700e+00,0.8800e+00,0.8900e+00,0.9000e+00,0.9100e+00, &
      0.9200e+00,0.9300e+00,0.9400e+00,0.9500e+00,0.9600e+00,0.9700e+00, &
      0.9800e+00,0.9900e+00,0.1000e+01,0.1010e+01,0.1020e+01,0.1030e+01, &
      0.1040e+01,0.1050e+01,0.1060e+01,0.1070e+01,0.1080e+01,0.1090e+01, &
      0.1100e+01,0.1110e+01,0.1120e+01,0.1130e+01,0.1140e+01,0.1150e+01, &
      0.1160e+01,0.1170e+01,0.1180e+01,0.1190e+01,0.1200e+01,0.1210e+01, &
      0.1220e+01,0.1230e+01,0.1240e+01,0.1250e+01,0.1260e+01,0.1270e+01, &
      0.1280e+01,0.1290e+01,0.1300e+01,0.1310e+01,0.1320e+01,0.1330e+01, &
      0.1340e+01,0.1350e+01,0.1360e+01,0.1370e+01,0.1380e+01,0.1390e+01, &
      0.1400e+01,0.1410e+01,0.1420e+01,0.1430e+01,0.1440e+01,0.1449e+01, &
      0.1460e+01,0.1471e+01,0.1481e+01,0.1493e+01,0.1504e+01,0.1515e+01, &
      0.1527e+01,0.1538e+01,0.1563e+01,0.1587e+01,0.1613e+01,0.1650e+01, &
      0.1680e+01,0.1700e+01,0.1730e+01,0.1760e+01,0.1800e+01,0.1830e+01, &
      0.1840e+01,0.1850e+01,0.1855e+01,0.1860e+01,0.1870e+01,0.1890e+01/
      data (wl(i),i=229,342)/ &
      0.1905e+01,0.1923e+01,0.1942e+01,0.1961e+01,0.1980e+01,0.2000e+01, &
      0.2020e+01,0.2041e+01,0.2062e+01,0.2083e+01,0.2105e+01,0.2130e+01, &
      0.2150e+01,0.2170e+01,0.2190e+01,0.2220e+01,0.2240e+01,0.2245e+01, &
      0.2250e+01,0.2260e+01,0.2270e+01,0.2290e+01,0.2310e+01,0.2330e+01, &
      0.2350e+01,0.2370e+01,0.2390e+01,0.2410e+01,0.2430e+01,0.2460e+01, &
      0.2500e+01,0.2520e+01,0.2550e+01,0.2565e+01,0.2580e+01,0.2590e+01, &
      0.2600e+01,0.2620e+01,0.2675e+01,0.2725e+01,0.2778e+01,0.2817e+01, &
      0.2833e+01,0.2849e+01,0.2865e+01,0.2882e+01,0.2899e+01,0.2915e+01, &
      0.2933e+01,0.2950e+01,0.2967e+01,0.2985e+01,0.3003e+01,0.3021e+01, &
      0.3040e+01,0.3058e+01,0.3077e+01,0.3096e+01,0.3115e+01,0.3135e+01, &
      0.3155e+01,0.3175e+01,0.3195e+01,0.3215e+01,0.3236e+01,0.3257e+01, &
      0.3279e+01,0.3300e+01,0.3322e+01,0.3345e+01,0.3367e+01,0.3390e+01, &
      0.3413e+01,0.3436e+01,0.3460e+01,0.3484e+01,0.3509e+01,0.3534e+01, &
      0.3559e+01,0.3624e+01,0.3732e+01,0.3775e+01,0.3847e+01,0.3969e+01, &
      0.4099e+01,0.4239e+01,0.4348e+01,0.4387e+01,0.4444e+01,0.4505e+01, &
      0.4547e+01,0.4560e+01,0.4580e+01,0.4719e+01,0.4904e+01,0.5000e+01, &
      0.5100e+01,0.5200e+01,0.5263e+01,0.5400e+01,0.5556e+01,0.5714e+01, &
      0.5747e+01,0.5780e+01,0.5814e+01,0.5848e+01,0.5882e+01,0.6061e+01, &
      0.6135e+01,0.6250e+01,0.6289e+01,0.6329e+01,0.6369e+01,0.6410e+01/
      data (wl(i),i=343,456)/ &
      0.6452e+01,0.6494e+01,0.6579e+01,0.6667e+01,0.6757e+01,0.6897e+01, &
      0.7042e+01,0.7143e+01,0.7246e+01,0.7353e+01,0.7463e+01,0.7576e+01, &
      0.7692e+01,0.7812e+01,0.7937e+01,0.8065e+01,0.8197e+01,0.8333e+01, &
      0.8475e+01,0.8696e+01,0.8929e+01,0.9091e+01,0.9259e+01,0.9524e+01, &
      0.9804e+01,0.1000e+02,0.1020e+02,0.1031e+02,0.1042e+02,0.1053e+02, &
      0.1064e+02,0.1075e+02,0.1087e+02,0.1100e+02,0.1111e+02,0.1136e+02, &
      0.1163e+02,0.1190e+02,0.1220e+02,0.1250e+02,0.1282e+02,0.1299e+02, &
      0.1316e+02,0.1333e+02,0.1351e+02,0.1370e+02,0.1389e+02,0.1408e+02, &
      0.1429e+02,0.1471e+02,0.1515e+02,0.1538e+02,0.1563e+02,0.1613e+02, &
      0.1639e+02,0.1667e+02,0.1695e+02,0.1724e+02,0.1818e+02,0.1887e+02, &
      0.1923e+02,0.1961e+02,0.2000e+02,0.2041e+02,0.2083e+02,0.2222e+02, &
      0.2260e+02,0.2305e+02,0.2360e+02,0.2460e+02,0.2500e+02,0.2600e+02, &
      0.2857e+02,0.3100e+02,0.3333e+02,0.3448e+02,0.3564e+02,0.3700e+02, &
      0.3824e+02,0.3960e+02,0.4114e+02,0.4276e+02,0.4358e+02,0.4458e+02, &
      0.4550e+02,0.4615e+02,0.4671e+02,0.4736e+02,0.4800e+02,0.4878e+02, &
      0.5003e+02,0.5128e+02,0.5275e+02,0.5350e+02,0.5424e+02,0.5500e+02, &
      0.5574e+02,0.5640e+02,0.5700e+02,0.5746e+02,0.5840e+02,0.5929e+02, &
      0.6000e+02,0.6100e+02,0.6125e+02,0.6250e+02,0.6378e+02,0.6467e+02, &
      0.6558e+02,0.6655e+02,0.6760e+02,0.6900e+02,0.7053e+02,0.7300e+02/
      data (wl(i),i=457,468)/ &
      0.7500e+02,0.7629e+02,0.8000e+02,0.8297e+02,0.8500e+02,0.8680e+02, &
      0.9080e+02,0.9517e+02,0.1000e+03,0.1200e+03,0.1500e+03,0.1670e+03/
      data  wlt/ &
                                       0.1670e+03,0.1778e+03,0.1884e+03, &
      0.1995e+03,0.2113e+03,0.2239e+03,0.2371e+03,0.2512e+03,0.2661e+03, &
      0.2818e+03,0.2985e+03,0.3162e+03,0.3548e+03,0.3981e+03,0.4467e+03, &
      0.5012e+03,0.5623e+03,0.6310e+03,0.7943e+03,0.1000e+04,0.1259e+04, &
      0.2500e+04,0.5000e+04,0.1000e+05,0.2000e+05,0.3200e+05,0.3500e+05, &
      0.4000e+05,0.4500e+05,0.5000e+05,0.6000e+05,0.7000e+05,0.9000e+05, &
      0.1110e+06,0.1200e+06,0.1300e+06,0.1400e+06,0.1500e+06,0.1600e+06, &
      0.1700e+06,0.1800e+06,0.2000e+06,0.2500e+06,0.2900e+06,0.3200e+06, &
      0.3500e+06,0.3800e+06,0.4000e+06,0.4500e+06,0.5000e+06,0.6000e+06, &
      0.6400e+06,0.6800e+06,0.7200e+06,0.7600e+06,0.8000e+06,0.8400e+06, &
      0.9000e+06,0.1000e+07,0.2000e+07,0.5000e+07,0.8600e+07/
      data (tabre(i),i=1,114)/ &
         0.83441,   0.83676,   0.83729,   0.83771,   0.83827,   0.84038, &
         0.84719,   0.85522,   0.86047,   0.86248,   0.86157,   0.86093, &
         0.86419,   0.86916,   0.87764,   0.89296,   0.91041,   0.93089, &
         0.95373,   0.98188,   1.02334,   1.06735,   1.11197,   1.13134, &
         1.15747,   1.20045,   1.23840,   1.27325,   1.32157,   1.38958, &
         1.41644,   1.40906,   1.40063,   1.40169,   1.40934,   1.40221, &
         1.39240,   1.38424,   1.38075,   1.38186,   1.39634,   1.40918, &
         1.40256,   1.38013,   1.36303,   1.34144,   1.32377,   1.30605, &
         1.29054,   1.28890,   1.28931,   1.30190,   1.32025,   1.36302, &
         1.41872,   1.45834,   1.49028,   1.52128,   1.55376,   1.57782, &
         1.59636,   1.60652,   1.61172,   1.61919,   1.62522,   1.63404, &
         1.63689,   1.63833,   1.63720,   1.63233,   1.62222,   1.58269, &
         1.55635,   1.52453,   1.50320,   1.48498,   1.47226,   1.45991, &
         1.45115,   1.44272,   1.43498,   1.43280,   1.42924,   1.42602, &
         1.42323,   1.42143,   1.41897,   1.41660,   1.41434,   1.41216, &
         1.41006,   1.40805,   1.40423,   1.40067,   1.38004,   1.35085, &
         1.33394,   1.32492,   1.31940,   1.31854,   1.31775,   1.31702, &
         1.31633,   1.31569,   1.31509,   1.31452,   1.31399,   1.31349, &
         1.31302,   1.31257,   1.31215,   1.31175,   1.31136,   1.31099/
      data (tabre(i),i=115,228)/ &
         1.31064,   1.31031,   1.30999,   1.30968,   1.30938,   1.30909, &
         1.30882,   1.30855,   1.30829,   1.30804,   1.30780,   1.30756, &
         1.30733,   1.30710,   1.30688,   1.30667,   1.30646,   1.30625, &
         1.30605,   1.30585,   1.30566,   1.30547,   1.30528,   1.30509, &
         1.30491,   1.30473,   1.30455,   1.30437,   1.30419,   1.30402, &
         1.30385,   1.30367,   1.30350,   1.30333,   1.30316,   1.30299, &
         1.30283,   1.30266,   1.30249,   1.30232,   1.30216,   1.30199, &
         1.30182,   1.30166,   1.30149,   1.30132,   1.30116,   1.30099, &
         1.30082,   1.30065,   1.30048,   1.30031,   1.30014,   1.29997, &
         1.29979,   1.29962,   1.29945,   1.29927,   1.29909,   1.29891, &
         1.29873,   1.29855,   1.29837,   1.29818,   1.29800,   1.29781, &
         1.29762,   1.29743,   1.29724,   1.29705,   1.29686,   1.29666, &
         1.29646,   1.29626,   1.29605,   1.29584,   1.29563,   1.29542, &
         1.29521,   1.29499,   1.29476,   1.29453,   1.29430,   1.29406, &
         1.29381,   1.29355,   1.29327,   1.29299,   1.29272,   1.29252, &
         1.29228,   1.29205,   1.29186,   1.29167,   1.29150,   1.29130, &
         1.29106,   1.29083,   1.29025,   1.28962,   1.28891,   1.28784, &
         1.28689,   1.28623,   1.28521,   1.28413,   1.28261,   1.28137, &
         1.28093,   1.28047,   1.28022,   1.27998,   1.27948,   1.27849/
      data (tabre(i),i=229,342)/ &
         1.27774,   1.27691,   1.27610,   1.27535,   1.27471,   1.27404, &
         1.27329,   1.27240,   1.27139,   1.27029,   1.26901,   1.26736, &
         1.26591,   1.26441,   1.26284,   1.26036,   1.25860,   1.25815, &
         1.25768,   1.25675,   1.25579,   1.25383,   1.25179,   1.24967, &
         1.24745,   1.24512,   1.24266,   1.24004,   1.23725,   1.23270, &
         1.22583,   1.22198,   1.21548,   1.21184,   1.20790,   1.20507, &
         1.20209,   1.19566,   1.17411,   1.14734,   1.10766,   1.06739, &
         1.04762,   1.02650,   1.00357,   0.98197,   0.96503,   0.95962, &
         0.97269,   0.99172,   1.00668,   1.02186,   1.04270,   1.07597, &
         1.12954,   1.21267,   1.32509,   1.42599,   1.49656,   1.55095, &
         1.59988,   1.63631,   1.65024,   1.64278,   1.62691,   1.61284, &
         1.59245,   1.57329,   1.55770,   1.54129,   1.52654,   1.51139, &
         1.49725,   1.48453,   1.47209,   1.46125,   1.45132,   1.44215, &
         1.43366,   1.41553,   1.39417,   1.38732,   1.37735,   1.36448, &
         1.35414,   1.34456,   1.33882,   1.33807,   1.33847,   1.34053, &
         1.34287,   1.34418,   1.34634,   1.34422,   1.33453,   1.32897, &
         1.32333,   1.31800,   1.31432,   1.30623,   1.29722,   1.28898, &
         1.28730,   1.28603,   1.28509,   1.28535,   1.28813,   1.30156, &
         1.30901,   1.31720,   1.31893,   1.32039,   1.32201,   1.32239/
      data (tabre(i),i=343,456)/ &
         1.32149,   1.32036,   1.31814,   1.31705,   1.31807,   1.31953, &
         1.31933,   1.31896,   1.31909,   1.31796,   1.31631,   1.31542, &
         1.31540,   1.31552,   1.31455,   1.31193,   1.30677,   1.29934, &
         1.29253,   1.28389,   1.27401,   1.26724,   1.25990,   1.24510, &
         1.22241,   1.19913,   1.17150,   1.15528,   1.13700,   1.11808, &
         1.10134,   1.09083,   1.08734,   1.09254,   1.10654,   1.14779, &
         1.20202,   1.25825,   1.32305,   1.38574,   1.44478,   1.47170, &
         1.49619,   1.51652,   1.53328,   1.54900,   1.56276,   1.57317, &
         1.58028,   1.57918,   1.56672,   1.55869,   1.55081,   1.53807, &
         1.53296,   1.53220,   1.53340,   1.53289,   1.51705,   1.50097, &
         1.49681,   1.49928,   1.50153,   1.49856,   1.49053,   1.46070, &
         1.45182,   1.44223,   1.43158,   1.41385,   1.40676,   1.38955, &
         1.34894,   1.31039,   1.26420,   1.23656,   1.21663,   1.20233, &
         1.19640,   1.19969,   1.20860,   1.22173,   1.24166,   1.28175, &
         1.32784,   1.38657,   1.46486,   1.55323,   1.60379,   1.61877, &
         1.62963,   1.65712,   1.69810,   1.72065,   1.74865,   1.76736, &
         1.76476,   1.75011,   1.72327,   1.68490,   1.62398,   1.59596, &
         1.58514,   1.59917,   1.61405,   1.66625,   1.70663,   1.73713, &
         1.76860,   1.80343,   1.83296,   1.85682,   1.87411,   1.89110/
      data (tabre(i),i=457,468)/ &
         1.89918,   1.90432,   1.90329,   1.88744,   1.87499,   1.86702, &
         1.85361,   1.84250,   1.83225,   1.81914,   1.82268,   1.82961/
      data (tabret(i,1),i=1,nwlt)/ &
                                          1.82961,   1.83258,   1.83149, &
         1.82748,   1.82224,   1.81718,   1.81204,   1.80704,   1.80250, &
         1.79834,   1.79482,   1.79214,   1.78843,   1.78601,   1.78434, &
         1.78322,   1.78248,   1.78201,   1.78170,   1.78160,   1.78190, &
         1.78300,   1.78430,   1.78520,   1.78620,   1.78660,   1.78680, &
         1.78690,   1.78700,   1.78700,   1.78710,   1.78710,   1.78720, &
         1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720, &
         1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720, &
         1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720, &
         1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720, &
         1.78720,   1.78720,   1.78720,   1.78720,   1.78800/
      data (tabret(i,2),i=1,nwlt)/ &
                               1.82961,   1.83258,   1.83149,   1.82748, &
         1.82224,   1.81718,   1.81204,   1.80704,   1.80250,   1.79834, &
         1.79482,   1.79214,   1.78843,   1.78601,   1.78434,   1.78322, &
         1.78248,   1.78201,   1.78170,   1.78160,   1.78190,   1.78300, &
         1.78430,   1.78520,   1.78610,   1.78630,   1.78640,   1.78650, &
         1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650, &
         1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650, &
         1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650, &
         1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650, &
         1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650, &
         1.78650,   1.78650,   1.78650,   1.78720/
      data(tabret(i,3),i=1,nwlt)/ &
                    1.82961,   1.83258,   1.83149,   1.82748,   1.82224, &
         1.81718,   1.81204,   1.80704,   1.80250,   1.79834,   1.79482, &
         1.79214,   1.78843,   1.78601,   1.78434,   1.78322,   1.78248, &
         1.78201,   1.78160,   1.78140,   1.78160,   1.78220,   1.78310, &
         1.78380,   1.78390,   1.78400,   1.78400,   1.78400,   1.78400, &
         1.78400,   1.78390,   1.78380,   1.78370,   1.78370,   1.78370, &
         1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370, &
         1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370, &
         1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370, &
         1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370, &
         1.78370,   1.78400,   1.78450/
      data (tabret(i,4),i=1,nwlt)/ &
         1.82961,   1.83258,   1.83149,   1.82748,   1.82224,   1.81718, &
         1.81204,   1.80704,   1.80250,   1.79834,   1.79482,   1.79214, &
         1.78843,   1.78601,   1.78434,   1.78322,   1.78248,   1.78201, &
         1.78150,   1.78070,   1.78010,   1.77890,   1.77790,   1.77730, &
         1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720, &
         1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720, &
         1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720, &
         1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720, &
         1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720, &
         1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720, &
         1.77720,   1.77800/
      data(tabim(i),i=1,114)/ &
      0.1640e+00,0.1730e+00,0.1830e+00,0.1950e+00,0.2080e+00,0.2230e+00, &
      0.2400e+00,0.2500e+00,0.2590e+00,0.2680e+00,0.2790e+00,0.2970e+00, &
      0.3190e+00,0.3400e+00,0.3660e+00,0.3920e+00,0.4160e+00,0.4400e+00, &
      0.4640e+00,0.4920e+00,0.5170e+00,0.5280e+00,0.5330e+00,0.5340e+00, &
      0.5310e+00,0.5240e+00,0.5100e+00,0.5000e+00,0.4990e+00,0.4680e+00, &
      0.3800e+00,0.3600e+00,0.3390e+00,0.3180e+00,0.2910e+00,0.2510e+00, &
      0.2440e+00,0.2390e+00,0.2390e+00,0.2440e+00,0.2470e+00,0.2240e+00, &
      0.1950e+00,0.1740e+00,0.1720e+00,0.1800e+00,0.1940e+00,0.2130e+00, &
      0.2430e+00,0.2710e+00,0.2890e+00,0.3340e+00,0.3440e+00,0.3820e+00, &
      0.4010e+00,0.4065e+00,0.4050e+00,0.3890e+00,0.3770e+00,0.3450e+00, &
      0.3320e+00,0.3150e+00,0.2980e+00,0.2740e+00,0.2280e+00,0.1980e+00, &
      0.1720e+00,0.1560e+00,0.1100e+00,0.8300e-01,0.5800e-01,0.2200e-01, &
      0.1000e-01,0.3000e-02,0.1000e-02,0.3000e-03,0.1000e-03,0.3000e-04, &
      0.1000e-04,0.3000e-05,0.1000e-05,0.7000e-06,0.4000e-06,0.2000e-06, &
      0.1000e-06,0.6377e-07,0.3750e-07,0.2800e-07,0.2400e-07,0.2200e-07, &
      0.1900e-07,0.1750e-07,0.1640e-07,0.1590e-07,0.1325e-07,0.8623e-08, &
      0.5504e-08,0.3765e-08,0.2710e-08,0.2510e-08,0.2260e-08,0.2080e-08, &
      0.1910e-08,0.1540e-08,0.1530e-08,0.1550e-08,0.1640e-08,0.1780e-08, &
      0.1910e-08,0.2140e-08,0.2260e-08,0.2540e-08,0.2930e-08,0.3110e-08/
      data(tabim(i),i=115,228)/ &
      0.3290e-08,0.3520e-08,0.4040e-08,0.4880e-08,0.5730e-08,0.6890e-08, &
      0.8580e-08,0.1040e-07,0.1220e-07,0.1430e-07,0.1660e-07,0.1890e-07, &
      0.2090e-07,0.2400e-07,0.2900e-07,0.3440e-07,0.4030e-07,0.4300e-07, &
      0.4920e-07,0.5870e-07,0.7080e-07,0.8580e-07,0.1020e-06,0.1180e-06, &
      0.1340e-06,0.1400e-06,0.1430e-06,0.1450e-06,0.1510e-06,0.1830e-06, &
      0.2150e-06,0.2650e-06,0.3350e-06,0.3920e-06,0.4200e-06,0.4440e-06, &
      0.4740e-06,0.5110e-06,0.5530e-06,0.6020e-06,0.7550e-06,0.9260e-06, &
      0.1120e-05,0.1330e-05,0.1620e-05,0.2000e-05,0.2250e-05,0.2330e-05, &
      0.2330e-05,0.2170e-05,0.1960e-05,0.1810e-05,0.1740e-05,0.1730e-05, &
      0.1700e-05,0.1760e-05,0.1820e-05,0.2040e-05,0.2250e-05,0.2290e-05, &
      0.3040e-05,0.3840e-05,0.4770e-05,0.5760e-05,0.6710e-05,0.8660e-05, &
      0.1020e-04,0.1130e-04,0.1220e-04,0.1290e-04,0.1320e-04,0.1350e-04, &
      0.1330e-04,0.1320e-04,0.1320e-04,0.1310e-04,0.1320e-04,0.1320e-04, &
      0.1340e-04,0.1390e-04,0.1420e-04,0.1480e-04,0.1580e-04,0.1740e-04, &
      0.1980e-04,0.2500e-04,0.5400e-04,0.1040e-03,0.2030e-03,0.2708e-03, &
      0.3511e-03,0.4299e-03,0.5181e-03,0.5855e-03,0.5899e-03,0.5635e-03, &
      0.5480e-03,0.5266e-03,0.4394e-03,0.3701e-03,0.3372e-03,0.2410e-03, &
      0.1890e-03,0.1660e-03,0.1450e-03,0.1280e-03,0.1030e-03,0.8600e-04, &
      0.8220e-04,0.8030e-04,0.8500e-04,0.9900e-04,0.1500e-03,0.2950e-03/
      data(tabim(i),i=229,342)/ &
      0.4687e-03,0.7615e-03,0.1010e-02,0.1313e-02,0.1539e-02,0.1588e-02, &
      0.1540e-02,0.1412e-02,0.1244e-02,0.1068e-02,0.8414e-03,0.5650e-03, &
      0.4320e-03,0.3500e-03,0.2870e-03,0.2210e-03,0.2030e-03,0.2010e-03, &
      0.2030e-03,0.2140e-03,0.2320e-03,0.2890e-03,0.3810e-03,0.4620e-03, &
      0.5480e-03,0.6180e-03,0.6800e-03,0.7300e-03,0.7820e-03,0.8480e-03, &
      0.9250e-03,0.9200e-03,0.8920e-03,0.8700e-03,0.8900e-03,0.9300e-03, &
      0.1010e-02,0.1350e-02,0.3420e-02,0.7920e-02,0.2000e-01,0.3800e-01, &
      0.5200e-01,0.6800e-01,0.9230e-01,0.1270e+00,0.1690e+00,0.2210e+00, &
      0.2760e+00,0.3120e+00,0.3470e+00,0.3880e+00,0.4380e+00,0.4930e+00, &
      0.5540e+00,0.6120e+00,0.6250e+00,0.5930e+00,0.5390e+00,0.4910e+00, &
      0.4380e+00,0.3720e+00,0.3000e+00,0.2380e+00,0.1930e+00,0.1580e+00, &
      0.1210e+00,0.1030e+00,0.8360e-01,0.6680e-01,0.5400e-01,0.4220e-01, &
      0.3420e-01,0.2740e-01,0.2200e-01,0.1860e-01,0.1520e-01,0.1260e-01, &
      0.1060e-01,0.8020e-02,0.6850e-02,0.6600e-02,0.6960e-02,0.9160e-02, &
      0.1110e-01,0.1450e-01,0.2000e-01,0.2300e-01,0.2600e-01,0.2900e-01, &
      0.2930e-01,0.3000e-01,0.2850e-01,0.1730e-01,0.1290e-01,0.1200e-01, &
      0.1250e-01,0.1340e-01,0.1400e-01,0.1750e-01,0.2400e-01,0.3500e-01, &
      0.3800e-01,0.4200e-01,0.4600e-01,0.5200e-01,0.5700e-01,0.6900e-01, &
      0.7000e-01,0.6700e-01,0.6500e-01,0.6400e-01,0.6200e-01,0.5900e-01/
      data(tabim(i),i=343,456)/ &
      0.5700e-01,0.5600e-01,0.5500e-01,0.5700e-01,0.5800e-01,0.5700e-01, &
      0.5500e-01,0.5500e-01,0.5400e-01,0.5200e-01,0.5200e-01,0.5200e-01, &
      0.5200e-01,0.5000e-01,0.4700e-01,0.4300e-01,0.3900e-01,0.3700e-01, &
      0.3900e-01,0.4000e-01,0.4200e-01,0.4400e-01,0.4500e-01,0.4600e-01, &
      0.4700e-01,0.5100e-01,0.6500e-01,0.7500e-01,0.8800e-01,0.1080e+00, &
      0.1340e+00,0.1680e+00,0.2040e+00,0.2480e+00,0.2800e+00,0.3410e+00, &
      0.3790e+00,0.4090e+00,0.4220e+00,0.4220e+00,0.4030e+00,0.3890e+00, &
      0.3740e+00,0.3540e+00,0.3350e+00,0.3150e+00,0.2940e+00,0.2710e+00, &
      0.2460e+00,0.1980e+00,0.1640e+00,0.1520e+00,0.1420e+00,0.1280e+00, &
      0.1250e+00,0.1230e+00,0.1160e+00,0.1070e+00,0.7900e-01,0.7200e-01, &
      0.7600e-01,0.7500e-01,0.6700e-01,0.5500e-01,0.4500e-01,0.2900e-01, &
      0.2750e-01,0.2700e-01,0.2730e-01,0.2890e-01,0.3000e-01,0.3400e-01, &
      0.5300e-01,0.7550e-01,0.1060e+00,0.1350e+00,0.1761e+00,0.2229e+00, &
      0.2746e+00,0.3280e+00,0.3906e+00,0.4642e+00,0.5247e+00,0.5731e+00, &
      0.6362e+00,0.6839e+00,0.7091e+00,0.6790e+00,0.6250e+00,0.5654e+00, &
      0.5433e+00,0.5292e+00,0.5070e+00,0.4883e+00,0.4707e+00,0.4203e+00, &
      0.3771e+00,0.3376e+00,0.3056e+00,0.2835e+00,0.3170e+00,0.3517e+00, &
      0.3902e+00,0.4509e+00,0.4671e+00,0.4779e+00,0.4890e+00,0.4899e+00, &
      0.4873e+00,0.4766e+00,0.4508e+00,0.4193e+00,0.3880e+00,0.3433e+00/
      data(tabim(i),i=457,468)/ &
      0.3118e+00,0.2935e+00,0.2350e+00,0.1981e+00,0.1865e+00,0.1771e+00, &
      0.1620e+00,0.1490e+00,0.1390e+00,0.1200e+00,0.9620e-01,0.8300e-01/
      data(tabimt(i,1),i=1,nwlt)/ &
                                       0.8300e-01,0.6900e-01,0.5700e-01, &
      0.4560e-01,0.3790e-01,0.3140e-01,0.2620e-01,0.2240e-01,0.1960e-01, &
      0.1760e-01,0.1665e-01,0.1620e-01,0.1550e-01,0.1470e-01,0.1390e-01, &
      0.1320e-01,0.1250e-01,0.1180e-01,0.1060e-01,0.9540e-02,0.8560e-02, &
      0.6210e-02,0.4490e-02,0.3240e-02,0.2340e-02,0.1880e-02,0.1740e-02, &
      0.1500e-02,0.1320e-02,0.1160e-02,0.8800e-03,0.6950e-03,0.4640e-03, &
      0.3400e-03,0.3110e-03,0.2940e-03,0.2790e-03,0.2700e-03,0.2640e-03, &
      0.2580e-03,0.2520e-03,0.2490e-03,0.2540e-03,0.2640e-03,0.2740e-03, &
      0.2890e-03,0.3050e-03,0.3150e-03,0.3460e-03,0.3820e-03,0.4620e-03, &
      0.5000e-03,0.5500e-03,0.5950e-03,0.6470e-03,0.6920e-03,0.7420e-03, &
      0.8200e-03,0.9700e-03,0.1950e-02,0.5780e-02,0.9700e-02/
      data(tabimt(i,2),i=1,nwlt)/ &
                            0.8300e-01,0.6900e-01,0.5700e-01,0.4560e-01, &
      0.3790e-01,0.3140e-01,0.2620e-01,0.2240e-01,0.1960e-01,0.1760e-01, &
      0.1665e-01,0.1600e-01,0.1500e-01,0.1400e-01,0.1310e-01,0.1230e-01, &
      0.1150e-01,0.1080e-01,0.9460e-02,0.8290e-02,0.7270e-02,0.4910e-02, &
      0.3300e-02,0.2220e-02,0.1490e-02,0.1140e-02,0.1060e-02,0.9480e-03, &
      0.8500e-03,0.7660e-03,0.6300e-03,0.5200e-03,0.3840e-03,0.2960e-03, &
      0.2700e-03,0.2520e-03,0.2440e-03,0.2360e-03,0.2300e-03,0.2280e-03, &
      0.2250e-03,0.2200e-03,0.2160e-03,0.2170e-03,0.2200e-03,0.2250e-03, &
      0.2320e-03,0.2390e-03,0.2600e-03,0.2860e-03,0.3560e-03,0.3830e-03, &
      0.4150e-03,0.4450e-03,0.4760e-03,0.5080e-03,0.5400e-03,0.5860e-03, &
      0.6780e-03,0.1280e-02,0.3550e-02,0.5600e-02/
      data(tabimt(i,3),i=1,nwlt)/ &
                 0.8300e-01,0.6900e-01,0.5700e-01,0.4560e-01,0.3790e-01, &
      0.3140e-01,0.2620e-01,0.2190e-01,0.1880e-01,0.1660e-01,0.1540e-01, &
      0.1470e-01,0.1350e-01,0.1250e-01,0.1150e-01,0.1060e-01,0.9770e-02, &
      0.9010e-02,0.7660e-02,0.6520e-02,0.5540e-02,0.3420e-02,0.2100e-02, &
      0.1290e-02,0.7930e-03,0.5700e-03,0.5350e-03,0.4820e-03,0.4380e-03, &
      0.4080e-03,0.3500e-03,0.3200e-03,0.2550e-03,0.2120e-03,0.2000e-03, &
      0.1860e-03,0.1750e-03,0.1660e-03,0.1560e-03,0.1490e-03,0.1440e-03, &
      0.1350e-03,0.1210e-03,0.1160e-03,0.1160e-03,0.1170e-03,0.1200e-03, &
      0.1230e-03,0.1320e-03,0.1440e-03,0.1680e-03,0.1800e-03,0.1900e-03, &
      0.2090e-03,0.2160e-03,0.2290e-03,0.2400e-03,0.2600e-03,0.2920e-03, &
      0.6100e-03,0.1020e-02,0.1810e-02/
      data(tabimt(i,4),i=1,nwlt)/ &
      0.8300e-01,0.6900e-01,0.5700e-01,0.4450e-01,0.3550e-01,0.2910e-01, &
      0.2440e-01,0.1970e-01,0.1670e-01,0.1400e-01,0.1235e-01,0.1080e-01, &
      0.8900e-02,0.7340e-02,0.6400e-02,0.5600e-02,0.5000e-02,0.4520e-02, &
      0.3680e-02,0.2990e-02,0.2490e-02,0.1550e-02,0.9610e-03,0.5950e-03, &
      0.3690e-03,0.2670e-03,0.2510e-03,0.2290e-03,0.2110e-03,0.1960e-03, &
      0.1730e-03,0.1550e-03,0.1310e-03,0.1130e-03,0.1060e-03,0.9900e-04, &
      0.9300e-04,0.8730e-04,0.8300e-04,0.7870e-04,0.7500e-04,0.6830e-04, &
      0.5600e-04,0.4960e-04,0.4550e-04,0.4210e-04,0.3910e-04,0.3760e-04, &
      0.3400e-04,0.3100e-04,0.2640e-04,0.2510e-04,0.2430e-04,0.2390e-04, &
      0.2370e-04,0.2380e-04,0.2400e-04,0.2460e-04,0.2660e-04,0.4450e-04, &
      0.8700e-04,0.1320e-03/
 
  pi = acos(-1.0)
  n_r=0.0
  n_i=0.0

! // convert frequency to wavelength (um)
  alam=3E5/freq
  if((alam < wlmin) .or. (alam > wlmax)) then
    print *, 'm_ice: wavelength out of bounds'
    stop
  endif

! // convert temperature to K
  tk = t + 273.16

  if (alam < cutice) then

!   // region from 0.045 microns to 167.0 microns - no temperature depend
    do i=2,nwl
      if(alam < wl(i)) continue
    enddo
    x1=log(wl(i-1))
    x2=log(wl(i))
    y1=tabre(i-1)
    y2=tabre(i)
    x=log(alam)
    y=((x-x1)*(y2-y1)/(x2-x1))+y1
    n_r=y
    y1=log(abs(tabim(i-1)))
    y2=log(abs(tabim(i)))
    y=((x-x1)*(y2-y1)/(x2-x1))+y1
    n_i=exp(y)

  else

!   // region from 167.0 microns to 8.6 meters - temperature dependence
    if(tk > temref(1)) tk=temref(1)
    if(tk < temref(4)) tk=temref(4)
    do 11 i=2,4
      if(tk.ge.temref(i)) go to 12
    11 continue
    12 lt1=i
    lt2=i-1
    do 13 i=2,nwlt
      if(alam.le.wlt(i)) go to 14
    13 continue
    14 x1=log(wlt(i-1))
    x2=log(wlt(i))
    y1=tabret(i-1,lt1)
    y2=tabret(i,lt1)
    x=log(alam)
    ylo=((x-x1)*(y2-y1)/(x2-x1))+y1
    y1=tabret(i-1,lt2)
    y2=tabret(i,lt2)
    yhi=((x-x1)*(y2-y1)/(x2-x1))+y1
    t1=temref(lt1)
    t2=temref(lt2)
    y=((tk-t1)*(yhi-ylo)/(t2-t1))+ylo
    n_r=y
    y1=log(abs(tabimt(i-1,lt1)))
    y2=log(abs(tabimt(i,lt1)))
    ylo=((x-x1)*(y2-y1)/(x2-x1))+y1
    y1=log(abs(tabimt(i-1,lt2)))
    y2=log(abs(tabimt(i,lt2)))
    yhi=((x-x1)*(y2-y1)/(x2-x1))+y1
    y=((tk-t1)*(yhi-ylo)/(t2-t1))+ylo
    n_i=exp(y)

  endif

  end subroutine m_ice

! ----------------------------------------------------------------------------
! subroutine MIEINT
! ----------------------------------------------------------------------------
!
!     General purpose Mie scattering routine for single particles
!     Author: R Grainger 1990
!     History:
!     G Thomas, March 2005: Added calculation of Phase function and
!     code to ensure correct calculation of backscatter coeficient
!     Options/Extend_Source
!
      Subroutine MieInt(Dx, SCm, Inp, Dqv, Dqxt, Dqsc, Dbsc, Dg, Xs1, Xs2, DPh, Error)

      Integer * 2  Imaxx
      Parameter (Imaxx = 12000)
      Real * 4     RIMax          ! largest real part of refractive index
      Parameter (RIMax = 2.5)
      Real * 4     IRIMax         ! largest imaginary part of refractive index
      Parameter (IRIMax = -2)
      Integer * 2  Itermax
      Parameter (Itermax = 12000 * 2.5)
                                ! must be large enough to cope with the
                                ! largest possible nmx = x * abs(scm) + 15
                                ! or nmx =  Dx + 4.05*Dx**(1./3.) + 2.0
      Integer * 2  Imaxnp
      Parameter (Imaxnp = 10000)  ! Change this as required
!     INPUT
      Real * 8     Dx
      Complex * 16  SCm
      Integer * 4  Inp
      Real * 8     Dqv(Inp)
!     OUTPUT
      Complex * 16  Xs1(InP)
      Complex * 16  Xs2(InP)
      Real * 8     Dqxt
      Real * 8     Dqsc
      Real * 8     Dg
      Real * 8     Dbsc
      Real * 8     DPh(InP)
      Integer * 4  Error
!     LOCAL
      Integer * 2  I
      Integer * 2  NStop
      Integer * 2  NmX
      Integer * 4  N    ! N*N > 32767 ie N > 181
      Integer * 4  Inp2
      Real * 8     Chi,Chi0,Chi1
      Real * 8     APsi,APsi0,APsi1
      Real * 8     Pi0(Imaxnp)
      Real * 8     Pi1(Imaxnp)
      Real * 8     Taun(Imaxnp)
      Real * 8     Psi,Psi0,Psi1
      Complex * 8  Ir
      Complex * 16 Cm
      Complex * 16 A,ANM1,APB
      Complex * 16 B,BNM1,AMB
      Complex * 16 D(Itermax)
      Complex * 16 Sp(Imaxnp)
      Complex * 16 Sm(Imaxnp)
      Complex * 16 Xi,Xi0,Xi1
      Complex * 16 Y
!     ACCELERATOR VARIABLES
      Integer * 2  Tnp1
      Integer * 2  Tnm1
      Real * 8     Dn
      Real * 8     Rnx
      Real * 8     S(Imaxnp)
      Real * 8     T(Imaxnp)
      Real * 8     Turbo
      Real * 8     A2
      Complex * 16 A1
      
      If ((Dx.Gt.Imaxx) .Or. (InP.Gt.ImaxNP)) Then
        Error = 1
        Return
      EndIf
      Cm = SCm
      Ir = 1 / Cm
      Y =  Dx * Cm
      If (Dx.Lt.0.02) Then
         NStop = 2
      Else
         If (Dx.Le.8.0) Then
            NStop = Dx + 4.00*Dx**(1./3.) + 2.0
         Else
            If (Dx.Lt. 4200.0) Then
               NStop = Dx + 4.05*Dx**(1./3.) + 2.0
            Else
               NStop = Dx + 4.00*Dx**(1./3.) + 2.0
            End If
         End If
      End If
      NmX = Max(Real(NStop),Real(Abs(Y))) + 15.
      If (Nmx .gt. Itermax) then
          Error = 1
          Return
      End If
      Inp2 = Inp+1
      D(NmX) = Dcmplx(0,0)
      Do N = Nmx-1,1,-1
         A1 = (N+1) / Y
         D(N) = A1 - 1/(A1+D(N+1))
      End Do
      Do I =1,Inp2
         Sm(I) = Dcmplx(0,0)
         Sp(I) = Dcmplx(0,0)
         Pi0(I) = 0
         Pi1(I) = 1
      End Do
      Psi0 = Cos(Dx)
      Psi1 = Sin(Dx)
      Chi0 =-Sin(Dx)
      Chi1 = Cos(Dx)
      APsi0 = Psi0
      APsi1 = Psi1
      Xi0 = Dcmplx(APsi0,Chi0)
      Xi1 = Dcmplx(APsi1,Chi1)
      Dg = 0
      Dqsc = 0
      Dqxt = 0
      Tnp1 = 1
      Do N = 1,Nstop
         DN = N
         Tnp1 = Tnp1 + 2
         Tnm1 = Tnp1 - 2
         A2 = Tnp1 / (DN*(DN+1D0))
         Turbo = (DN+1D0) / DN
         Rnx = DN/Dx
         Psi = Dble(Tnm1)*Psi1/Dx - Psi0
         APsi = Psi
         Chi = Tnm1*Chi1/Dx       - Chi0
         Xi = Dcmplx(APsi,Chi)
         A = ((D(N)*Ir+Rnx)*APsi-APsi1) / ((D(N)*Ir+Rnx)*  Xi-  Xi1)
         B = ((D(N)*Cm+Rnx)*APsi-APsi1) / ((D(N)*Cm+Rnx)*  Xi-  Xi1)
         Dqxt = Tnp1 *      Dble(A + B)          + Dqxt
         Dqsc = Tnp1 * (A*Conjg(A) + B*Conjg(B)) + Dqsc
         If (N.Gt.1) then
	    Dg = Dg + (dN*dN - 1) * Dble(ANM1*Conjg(A) + BNM1 * Conjg(B)) / dN + TNM1 * Dble(ANM1*Conjg(BNM1)) / (dN*dN - dN)
         End If
         Anm1 = A
         Bnm1 = B
         APB = A2 * (A + B)
         AMB = A2 * (A - B)
         Do I = 1,Inp2
            If (I.GT.Inp) Then
               S(I) = -Pi1(I)
            Else
               S(I) = Dqv(I) * Pi1(I)
            End If
            T(I) = S(I) - Pi0(I)
            Taun(I) = N*T(I) - Pi0(I)
            Sp(I) = APB * (Pi1(I) + Taun(I)) + Sp(I)
            Sm(I) = AMB * (Pi1(I) - Taun(I)) + Sm(I)
            Pi0(I) = Pi1(I)
            Pi1(I) = S(I) + T(I)*Turbo
         End Do
         Psi0 = Psi1
         Psi1 = Psi
         Apsi1 = Psi1
         Chi0 = Chi1
         Chi1 = Chi
         Xi1 = Dcmplx(APsi1,Chi1)
      End Do
      If (Dg .GT.0) Dg = 2 * Dg / Dqsc
      Dqsc =  2 * Dqsc / Dx**2
      Dqxt =  2 * Dqxt / Dx**2
      Do I = 1,Inp
         Xs1(I) = (Sp(I)+Sm(I)) / 2
         Xs2(I) = (Sp(I)-Sm(I)) / 2
         Dph(I) = 2 * Dble(Xs1(I)*Conjg(Xs1(I)) + Xs2(I)*Conjg(Xs2(I))) / (Dx**2 * Dqsc)
      End Do
      Dbsc = 4 * Abs(( (Sp(Inp2)+Sm(Inp2))/2 )**2) / Dx**2
      Error = 0
      Return
      End subroutine MieInt

  end module optics_lib
