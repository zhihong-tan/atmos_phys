        module aer_ccn_act_mod

implicit none
private
    private Loading, CalcG, aer_ccn_act_init, erff, CalcAlphaGamma, &
      CalcBeta
      
    public aer_ccn_act, aer_ccn_act2

!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id: aer_ccn_act.F90,v 13.0.2.1 2006/05/26 02:46:18 wfc Exp $'
character(len=128) :: tagname = '$Name: memphis_2006_08 $'

!---------------- private data -------------------

  integer, parameter :: TY = 3  !  Number of aerosol types
  integer, parameter :: MD = 2  !  Number of lognormal modes

  real :: T = 283.15 !  Temperature (K)
  real :: P = 0.800e5 !  Pressure (Pa)
  real, parameter :: R = 8.314  !  Gas constant (J mol-1 K-1)
  real, parameter :: ZERO = 273.15 !  Zero degree C (K)
  real, parameter :: ATM = 1.01325e5 !  Standard atmosphere pressure (Pa)
  real, parameter :: PI = 3.1415926

  
!NO 1 Ammonium Sulfate  NO 2 Sea Salt NO 3 Organics
  
  real :: B_term (TY) =(/0.7822,0.6342,1.3764/) ! 2 * 0.3492/(Bprim)**(1/3)
  real :: Mass_scal (TY) =(/0.15896,0.198,0.1241/) ! scaling mass (ug m-3)

  real :: N(MD) =(/340., 60./) ! Total Number Concen (cm-3)
  real :: Dm(MD) =(/0.01, 0.07/) ! Geometric mean diameter (micron)
  real :: LNSIGMA(MD) =(/0.47, 0.6931/) ! ln( Sigma (St. Dev.) )

!Parameters for look-up tables

  integer, parameter :: res = 20 !
  real ::  lowup=0.3 !m/s
  real ::  highup=10.
  real, dimension(res,res,res,res) :: droplets2

  integer, parameter :: res2 = 20 !
  real ::  lowup2=0.0001 !m/s
  real ::  highup2=0.3
  real ::  lowmass2=0.01 !ug m-3
  real ::  highmass2=1000.
  real ::  lowmass3=0.01 !ug m-3
  real ::  highmass3=1000.
  real :: lowT2=243.15 !K
  real :: highT2=308.15
  real, dimension(res2,res2,res2,res2) :: droplets


!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

logical :: module_is_initialized  = .false.
 
contains

subroutine aer_ccn_act (T1, P1, Updraft1, TotalMass, Drop)
real, dimension(TY), intent(inout) :: TotalMass
real, intent(in) :: T1, P1, Updraft1
real, intent(inout) :: Drop
    
real number, tmass, tmass2, updr, temp, sum
integer nomass, nomass2, noup, noT, i
real, dimension(3) :: Drop1
        
  if(.not. module_is_initialized) call aer_ccn_act_init()

  tmass=(TotalMass(1)+TotalMass(2)+TotalMass(3))*1.e12
    
  if (Updraft1>lowup2 .and. tmass>lowmass2) then

    tmass=(TotalMass(1)+TotalMass(2))*1.e12
    tmass2=TotalMass(3)*1.e12
    
    if (Updraft1>highup2) then
    
!      updr=max(min(Updraft1,highup),lowup)
!      temp=max(min(T1,highT2),lowT2)
      updr=max(min(Updraft1,highup-1.e-5),lowup)
      temp=max(min(T1,highT2-1.e-5),lowT2)
    
      noup= log(updr/lowup)/log(highup/lowup)*(res2-1.)
      noT= (temp-lowT2)*(res2-1)/(highT2-lowT2)

!      tmass=max(min(tmass,highmass2),lowmass2)
      tmass=max(min(tmass,highmass2-1.e-5),lowmass2)
      nomass= log(tmass/lowmass2)/log(highmass2/lowmass2)*(res2-1.)

!      tmass2=max(min(tmass2,highmass3),lowmass3)
      tmass2=max(min(tmass2,highmass3-1.e-5),lowmass3)
      nomass2= log(tmass2/lowmass3)/log(highmass3/lowmass3)*(res2-1.)
    
      Drop = 0.2*(droplets2(nomass+1,nomass2+1,noup+1,noT+1)+&
                  droplets2(nomass+2,nomass2+1,noup+1,noT+1)+ &
                  droplets2(nomass+1,nomass2+2,noup+1,noT+1)+ &
                  droplets2(nomass+1,nomass2+1,noup+2,noT+1)+ &
                  droplets2(nomass+1,nomass2+1,noup+1,noT+2))
  
    else

      updr=max(min(Updraft1,highup2),lowup2)
      temp=max(min(T1,highT2),lowT2)
    
      noup= log(updr/lowup2)/log(highup2/lowup2)*(res2-1.)
      noT= (temp-lowT2)*(res2-1)/(highT2-lowT2)

      tmass=max(min(tmass,highmass2),lowmass2)
      nomass= log(tmass/lowmass2)/log(highmass2/lowmass2)*(res2-1.)

      tmass2=max(min(tmass2,highmass3),lowmass3)
      nomass2= log(tmass2/lowmass3)/log(highmass3/lowmass3)*(res2-1.)
    
      Drop = 0.2*(droplets(nomass+1,nomass2+1,noup+1,noT+1)+droplets(nomass+2,nomass2+1,noup+1,noT+1)+ &
                  droplets(nomass+1,nomass2+2,noup+1,noT+1)+droplets(nomass+1,nomass2+1,noup+2,noT+1)+ &
                  droplets(nomass+1,nomass2+1,noup+1,noT+2))
  
    endif
  
  
  else
  
    Drop=0.
  
  endif
  
end subroutine aer_ccn_act

subroutine aer_ccn_act2 (T1, P1, Updraft1, TotalMass, mu,airdens,Nc,qc,qt,qe,tc,te,Drop)

!T1 temperature (K)
!P1 pressure (Pa)
!Updraft1 updraft velocity (m/s)
!TotalMass aerosol mass ()
!mu entrainment coef. (/s)
!airdens air density (kg/m3 air)
!Nc droplet mixing ratio (#/kg air)
!qc in-cloud vapor mixing ratio (kg water/kg air)
!qt in-cloud total water mixing ratio qc + ql (kg water/kg air)
!qe environment vapor mixing ratio (kg water/kg air)
!tc in-cloud temperature (K)
!te environment temperature (K)
!Drop droplet number concentration (#/cc)

real, dimension(TY), intent(in) :: TotalMass
real, intent(in) :: T1, P1, Updraft1, mu,airdens, Nc, qc, qt, qe, tc, te
real, intent(inout) :: Drop

real :: G, alpha, gamma, Smax
real :: Diam, beta, Le_cpa, Dcut
integer :: i, j
        
  if(.not. module_is_initialized) call aer_ccn_act_init()

  Drop=0.
  
  if (Nc > 0.) then    
    T = T1
    P = P1

    call CalcAlphaGamma(alpha, gamma)
    call CalcBeta(beta, Le_cpa)
              
! Diam  average diameter of droplets (micron)
    if (qt > qc) then
      Diam= ((qt-qc)/Nc*1.91e15)**(1./3.)    
    else
      Diam= 20.
    endif

!set the upper and lower limits of Diam
    if (Diam < 10.) Diam=10.
  
    call calcG(Diam, G)
    
    Smax=(alpha-gamma*mu*(qt-qe)*airdens+beta*mu*(Le_cpa*(qc-qe)+(tc-te)))*Updraft1/ &
         (gamma)/(0.5*3.1415*1.e3*G*Diam*1.e-6*Nc*airdens)
  
    if (Smax>0.) then
      do i=1,TY
        Dcut=B_term(i)/T/(Smax**(2./3.))
        do j=1, MD
          Drop=Drop+TotalMass(i)/Mass_scal(i)*N(j)*0.5* &
               (1.-erff(log(Dcut/Dm(j))/LNSIGMA(j)*0.707107))      
        end do
      end do
    endif
  endif

end subroutine aer_ccn_act2

subroutine aer_ccn_act_init ()
    call Loading()
    module_is_initialized  = .true.
end subroutine


subroutine Loading()

real xx
integer i, j, k, l

  open(11, FILE='INPUT/droplets.dat')
  do k=1,res2
    do i=1,res2
      do j=1, res2
        do l=1, res2
          read(11,*) xx
          droplets(k,i,j,l)=xx
        end do
      end do
    end do
  end do
  close(11)

  open(11, FILE='INPUT/droplets2.dat')
  do k=1,res2
    do i=1,res2
      do j=1, res2
        do l=1, res2
          read(11,*) xx
          droplets2(k,i,j,l)=xx
        end do
      end do
    end do
  end do
  close(11)

end subroutine Loading

subroutine CalcG(Dp, G)

real, intent(inout) :: Dp
real, intent(inout) :: G
real :: rhow = 1.0e3  ! density of water (Kg m-3)
real :: Mw = 0.018  ! molecular weight of water (Kg mol-1)
real :: alpc = 1.0  ! mass accomodation coef. 
real :: alpt = 0.97  ! thermal accomodation coef.
!real :: alpc = 0.042  ! mass accomodation coef. 
!real :: alpt = 1.  ! thermal accomodation coef.
!real :: alpc = 0.2  ! mass accomodation coef. 
!real :: alpt = 1.  ! thermal accomodation coef.
real :: delt = 2.16e-1 !thermal jump (micron)
real :: delv = 1.096e-1 !vapor jump (micron)
real vpres, Dv, ka, Le, mass, heat, TC
      
      
      Dv = 0.211/(P/ATM)*(T/ZERO)**1.94*1e-4  ! diffusivity of water vapor (m2 s-1)
      Dv = Dv/(Dp/(Dp+delv*2.)+2*Dv/(alpc*(Dp*1e-6))*(2.*PI*Mw/R/T)**0.5)
      ka = 1e-3*(4.39+0.071*T)  ! thermal conductivity (J m-1 s-1 K-1)
      ka = ka/(Dp/(Dp+delt*2.)+2*ka/(alpt*(Dp*1e-6)*1.007e3*P)*(2.*PI*R*T/0.028965)**0.5)
      TC = T-ZERO
      vpres = (6.107799961+TC*(4.436518521e-1+TC*(1.428945805e-2+TC*(2.650648471e-4 &
              +TC*(3.031240396e-6+TC*(2.034080948e-8+6.136820929e-11*TC))))))*1e2  ! saturated water vapor pressure(Pa)
      Le = 597.3*(ZERO/T)**(0.167+3.67e-4*T)*4.182*1e3  ! latent heat of water (J Kg -1) 
      
      mass = rhow*R*T/(vpres*Dv*Mw)
      heat = Le*rhow/(ka*T)*(Le*Mw/T/R-1)
!      print *, Dv, vpres, Mw, rhow, mass, heat
      G = 4./(mass+heat) ! (m2 s-1)
end subroutine CalcG

recursive function erff(x) RESULT(y)

! Error function from Numerical Recipes.
! erf(x) = 1 - erfc(x)

real dumerfc, x
real t, z, y


  z = abs(x)
  t = 1.0 / ( 1.0 + 0.5 * z )

  dumerfc =     t * exp(-z * z - 1.26551223 + t *      &
            ( 1.00002368 + t * ( 0.37409196 + t *    &
            ( 0.09678418 + t * (-0.18628806 + t *    &
            ( 0.27886807 + t * (-1.13520398 + t *    &
            ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))

  if ( x.lt.0.0 ) dumerfc = 2.0 - dumerfc
 
  y = 1.0 - dumerfc

end function erff


subroutine CalcAlphaGamma(alpha, gamma)

  real, intent(inout) :: alpha, gamma
  real :: rhow = 1.0e3  ! density of water (Kg m-3)
  real rhoa ! density of air (Kg m-3)
  real :: Cpa = 1.007e3 ! specific heat of air (J Kg-1 K-1)
  real :: Mw = 0.018  ! molecular weight of water (Kg mol-1)
  real :: Ma = 0.028965  ! molecular weight of air (Kg mol-1)
  real :: alpc = 1.  ! mass accomodation coef. 
  real :: g = 9.815 ! gravitational acceleration (m s-2) 
  real vpres, Dv, ka, Le, TC
  
  rhoa = P*Ma/R/T  ! (Kg m-3)
  Dv = 0.211/(P/ATM)*(T/ZERO)**1.94*1e-4  ! diffusivity of water vapor (m2 s-1)
!  Dv = Dv/(1+2*Dv/(alpc*Dp)*(2.*PI*Mw/R/T)**0.5)
  ka = 1e-3*(4.39+0.071*T)  ! thermal conductivity (J m-1 s-1 K-1)
  TC = T-ZERO
  vpres = (6.107799961+TC*(4.436518521e-1+TC*(1.428945805e-2+TC*(2.650648471e-4 &
          +TC*(3.031240396e-6+TC*(2.034080948e-8+6.136820929e-11*TC))))))*1e2  ! saturated water vapor pressure (Pa)
  Le = 597.3*(ZERO/T)**(0.167+3.67e-4*T)*4.182*1e3  ! latent heat of water (J Kg -1)
  alpha = g*Mw*Le/(Cpa*R*T**2.)-g*Ma/(R*T) ! (m-1)
  gamma = R*T/(vpres*Mw)+Mw*Le**2./(Cpa*P*Ma*T) ! (m3 Kg-1)
end subroutine CalcAlphaGamma

subroutine CalcBeta(beta, Le_cpa)

  real, intent(inout) :: beta, Le_cpa 
  real :: rhow = 1.0e3  ! density of water (Kg m-3)
  real rhoa ! density of air (Kg m-3)
  real :: Cpa = 1.007e3 ! specific heat of air (J Kg-1 K-1)
  real :: Mw = 0.018  ! molecular weight of water (Kg mol-1)
  real :: Ma = 0.028965  ! molecular weight of air (Kg mol-1)
  real :: alpc = 1.  ! mass accomodation coef. 
  real :: g = 9.815 ! gravitational acceleration (m s-2) 
  real vpres, Dv, ka, Le, TC
  
  rhoa = P*Ma/R/T  ! (Kg m-3)
  Dv = 0.211/(P/ATM)*(T/ZERO)**1.94*1e-4  ! diffusivity of water vapor (m2 s-1)
!      Dv = Dv/(1+2*Dv/(alpc*Dp)*(2.*PI*Mw/R/T)**0.5)
  ka = 1e-3*(4.39+0.071*T)  ! thermal conductivity (J m-1 s-1 K-1)
  TC = T-ZERO
  vpres = (6.107799961+TC*(4.436518521e-1+TC*(1.428945805e-2+TC*(2.650648471e-4 &
          +TC*(3.031240396e-6+TC*(2.034080948e-8+6.136820929e-11*TC))))))*1e2  ! saturated water vapor pressure (Pa)
  Le = 597.3*(ZERO/T)**(0.167+3.67e-4*T)*4.182*1e3  ! latent heat of water (J Kg -1)
  Le_cpa = Le/Cpa
  beta = Mw*Le/(R*T**2.) ! (K-1)
end subroutine CalcBeta


end module aer_ccn_act_mod
