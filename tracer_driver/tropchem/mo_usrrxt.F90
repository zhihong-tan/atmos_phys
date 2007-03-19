      module mo_usrrxt_mod

      use sat_vapor_pres_mod, only : escomp
      use constants_mod, only : rdgas, rvgas
      
implicit none
      private
      public :: usrrxt_init, usrrxt

!     save

      integer :: usr1_ndx, usr2_ndx, usr3_ndx, usr5_ndx, usr6_ndx, usr7_ndx, &
                 usr8_ndx, usr9_ndx, usr11_ndx, usr12_ndx, usr14_ndx, usr15_ndx, &
                 usr16_ndx, usr17_ndx, usr21_ndx, usr22_ndx, &
                 usr24_ndx, usr25_ndx, &
                 so4_ndx, h2o_ndx, &
                 strat37_ndx, strat38_ndx, strat72_ndx, strat73_ndx, strat74_ndx, &
                 strat75_ndx, strat76_ndx, strat77_ndx, strat78_ndx, strat79_ndx, &
                 strat80_ndx

      real, parameter :: d622 = rdgas/rvgas
      real, parameter :: d378 = 1. - d622     

character(len=128), parameter :: version     = '$Id: mo_usrrxt.F90,v 14.0 2007/03/15 22:11:18 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: nalanda $'
logical                       :: module_is_initialized = .false.

      contains

      subroutine usrrxt_init
!-----------------------------------------------------------------
!        ... Intialize the user reaction constants module
!-----------------------------------------------------------------

      use mo_chem_utls_mod, only : get_rxt_ndx, get_spc_ndx

      implicit none

      usr1_ndx = get_rxt_ndx( 'usr1' )
      usr2_ndx = get_rxt_ndx( 'usr2' )
      usr3_ndx = get_rxt_ndx( 'usr3' )
      usr5_ndx = get_rxt_ndx( 'usr5' )
      usr6_ndx = get_rxt_ndx( 'usr6' )
      usr7_ndx = get_rxt_ndx( 'usr7' )
      usr8_ndx = get_rxt_ndx( 'usr8' )
      usr9_ndx = get_rxt_ndx( 'usr9' )
      usr11_ndx = get_rxt_ndx( 'usr11' )
      usr12_ndx = get_rxt_ndx( 'usr12' )
      usr14_ndx = get_rxt_ndx( 'usr14' )
      usr15_ndx = get_rxt_ndx( 'usr15' )
      usr16_ndx = get_rxt_ndx( 'usr16' )
      usr17_ndx = get_rxt_ndx( 'usr17' )
      usr21_ndx = get_rxt_ndx( 'usr21' )
      usr22_ndx = get_rxt_ndx( 'usr22' )
      so4_ndx = get_spc_ndx( 'SO4' )
      h2o_ndx = get_spc_ndx( 'H2O' )
      usr24_ndx = get_rxt_ndx( 'usr24' )
      usr25_ndx = get_rxt_ndx( 'usr25' )
      strat37_ndx = get_rxt_ndx( 'strat37' )
      strat38_ndx = get_rxt_ndx( 'strat38' )
      strat72_ndx = get_rxt_ndx( 'strat72' )
      strat73_ndx = get_rxt_ndx( 'strat73' )
      strat74_ndx = get_rxt_ndx( 'strat74' )
      strat75_ndx = get_rxt_ndx( 'strat75' )
      strat76_ndx = get_rxt_ndx( 'strat76' )
      strat77_ndx = get_rxt_ndx( 'strat77' )
      strat78_ndx = get_rxt_ndx( 'strat78' )
      strat79_ndx = get_rxt_ndx( 'strat79' )
      strat80_ndx = get_rxt_ndx( 'strat80' )

      write(*,*) ' '
      write(*,*) 'usrrxt_init: diagnostics '
      write(*,'(10i5)') usr1_ndx, usr2_ndx, usr3_ndx, usr5_ndx, usr6_ndx, usr7_ndx, &
                 usr8_ndx, usr9_ndx, usr11_ndx, usr12_ndx, usr14_ndx, usr15_ndx, &
                 usr16_ndx, usr17_ndx, usr21_ndx, usr22_ndx, &
                 usr24_ndx, usr25_ndx, &
                 strat37_ndx, strat38_ndx, strat72_ndx, strat73_ndx, strat74_ndx, &
                 strat75_ndx, strat76_ndx, strat77_ndx, strat78_ndx, strat79_ndx, &
                 strat80_ndx

      end subroutine usrrxt_init

      subroutine usrrxt( rxt, temp, invariants, h2ovmr,  &
                         pmid, m, sulfate, qin, &
                         sh, &
                         plonl )
!-----------------------------------------------------------------
!        ... set the user specified reaction rates
!-----------------------------------------------------------------

      use chem_mods_mod, only : nfs, rxntot, indexh2o, indexm
      use constants_mod, only : PI
      use m_rxt_id_mod

      implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in) ::   qin(:,:,:)           ! transported species ( vmr )
      real, intent(in) ::   temp(:,:), &         ! temperature
                            m(:,:), &            ! total atm density
                            sulfate(:,:), &      ! sulfate aerosol vmr
                            h2ovmr(:,:), &       ! water vapor vmr
                            pmid(:,:), &         ! midpoint pressure in pa
                            sh(:,:), &           ! specific humidity
                            invariants(:,:,:)    ! invariants density
      real, intent(inout) ::  rxt(:,:,:)         ! gas phase rates
      
!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------
      real, parameter :: boltz = 1.38044e-16            ! erg / k
      real, parameter :: avo   = 6.023e23               ! molecules/mole
!-----------------------------------------------------------------
!        ... density of sulfate aerosol
!-----------------------------------------------------------------
!     real, parameter :: gam1 = 0.04                    ! n2o5+sul ->2hno3
      real, parameter :: gam1 = 0.10                    ! n2o5+sul ->2hno3
      real, parameter :: gam4 = 0.05                    ! NH3 +SUL ->NH4SO4 (Dentener 1994)
      real, parameter :: wso4 = 98.
      real, parameter :: den  = 1.15                    ! each molecule of so4(aer) density g/cm3
!-------------------------------------------------
!         ... volume of sulfate particles
!           assuming mean rm 
!           continient 0.05um  0.07um  0.09um
!           ocean      0.09um  0.25um  0.37um
!                      0.16um                  blake jgr,7195, 1995
!-------------------------------------------------
      real, parameter :: rm1  = 0.16*1.e-4                   ! mean radii in cm
      real, parameter :: fare = 4.*3.14*rm1*rm1              ! each mean particle(r=0.1u) area   cm2/cm3
      real, parameter :: dg   = 0.1                          ! mole diffusion =0.1 cm2 (Dentener, 1993)

      integer  ::  i, k
      real     ::  amas
      real, dimension( SIZE(temp,1) ) :: &
                   tp, &                    ! 300/t
                   tinv, &                  ! 1/t
                   ko, &
                   kinf, &
                   fc, &
!                  relhum, &                ! relative humidity
                   satq, &                  ! saturation specific humidity
                   satv, &                  ! saturation vapor pressure
                   xr, &                    ! factor to increase particle radii depending on rel hum
                   sur, &                   ! sulfate particle surface area (cm^2/cm^3)
                   exp_fac                  ! vector exponential
      integer :: plev
      real, dimension(SIZE(temp,1),SIZE(temp,2)) :: &
                   relhum                   ! relative humidity


      plev = SIZE(temp,2)
      
      amas = 4.*PI*rm1**3*den/3.            ! each mean particle(r=0.1u) mass (g)
!-----------------------------------------------------------------
!        ... o + o2 + m --> o3 + m
!-----------------------------------------------------------------
      do k = 1,plev
         tinv(:)           = 1. / temp(:,k)
         tp(:)             = 300. * tinv(:)
         if( usr1_ndx > 0 ) then
            rxt(:,k,usr1_ndx) = 6.e-34 * tp(:)**2.4
         end if
#ifdef IBM
!-----------------------------------------------------------------
!        ... n2o5 + m --> no2 + no3 + m
!-----------------------------------------------------------------
         if( usr3_ndx > 0 ) then
            if( usr2_ndx > 0 ) then
               call vexp( exp_fac, -10991.*tinv, plonl )
               rxt(:,k,usr3_ndx) = rxt(:,k,usr2_ndx) * 3.333e26 * exp_fac(:)
            else
               rxt(:,k,usr3_ndx) = 0.
            end if
         end if

!-----------------------------------------------------------------
!        set rates for:
!         ... hno3 + oh --> no3 + h2o
!           ho2no2 + m --> ho2 + no2 + m
!           co + oh --> co2 + ho2
!-----------------------------------------------------------------
         if( usr5_ndx > 0 ) then
            call vexp( exp_fac, 1335.*tinv, plonl )
            ko(:) = m(:,k) * 6.5e-34 * exp_fac(:)
            call vexp( exp_fac, 2199.*tinv, plonl )
            ko(:) = ko(:) / (1. + ko(:)/(2.7e-17*exp_fac(:)))
            call vexp( exp_fac, 460.*tinv, plonl )
            rxt(:,k,usr5_ndx) = ko(:) + 2.4e-14*exp_fac(:)
         end if
         if( usr7_ndx > 0 ) then
            if( usr6_ndx > 0 ) then
               call vexp( exp_fac, -10900.*tinv, plonl )
               rxt(:,k,usr7_ndx) = rxt(:,k,usr6_ndx) * exp_fac(:) / 2.1e-27
            else
               rxt(:,k,usr7_ndx) = 0.
            end if
         end if
         if( usr8_ndx > 0 ) then
            rxt(:,k,usr8_ndx) = 1.5e-13 * (1. + 6.e-7*boltz*m(:,k)*temp(:,k))
         end if

!-----------------------------------------------------------------
!        ... ho2 + ho2 --> h2o2
!        note: this rate involves the water vapor number density
!-----------------------------------------------------------------
         if( usr9_ndx > 0 ) then
            if( indexh2o > 0 ) then
               call vexp( exp_fac, 2200.*tinv, plonl )
               fc(:)   = 1. + 1.4e-21 * invariants(:,k,indexh2o) * exp_fac(:)
            else if( h2o_ndx > 0 ) then
               call vexp( exp_fac, 2200.*tinv, plonl )
               fc(:)   = 1. + 1.4e-21 * qin(:,k,h2o_ndx) * invariants(:,k,indexm) * exp_fac(:)
            else
               fc(:) = 1.
            end if
            call vexp( exp_fac, 600.*tinv, plonl )
            ko(:)   = 2.3e-13 * exp_fac(:)
            call vexp( exp_fac, 1000.*tinv, plonl )
            kinf(:) = 1.7e-33 * m(:,k) * exp_fac(:)
            rxt(:,k,usr9_ndx) = (ko(:) + kinf(:)) * fc(:)
         end if

!-----------------------------------------------------------------
!            ... mco3 + no2 -> mpan
!-----------------------------------------------------------------
         if( usr14_ndx > 0 ) then
            rxt(:,k,usr14_ndx) = 1.1e-11 * tp(:) / m(:,k)
         end if

!-----------------------------------------------------------------
!        ... pan + m --> ch3co3 + no2 + m
!-----------------------------------------------------------------
         call vexp( exp_fac, -14000.*tinv, plonl )
         if( usr12_ndx > 0 ) then
            if( usr11_ndx > 0 ) then
               rxt(:,k,usr12_ndx) = rxt(:,k,usr11_ndx) * 1.111e28 * exp_fac(:)
            else
               rxt(:,k,usr12_ndx) = 0.
            end if
         end if

!-----------------------------------------------------------------
!        ... mpan + m --> mco3 + no2 + m
!-----------------------------------------------------------------
         if( usr15_ndx > 0 ) then
            if( usr14_ndx > 0 ) then
               rxt(:,k,usr15_ndx) = rxt(:,k,usr14_ndx) * 1.111e28 * exp_fac(:)
            else
               rxt(:,k,usr15_ndx) = 0.
            end if
         end if

!-----------------------------------------------------------------
!       ... xooh + oh -> h2o + oh
!-----------------------------------------------------------------
         if( usr21_ndx > 0 ) then
            call vexp( exp_fac, 253.*tinv, plonl )
            rxt(:,k,usr21_ndx) = temp(:,k)**2 * 7.69e-17 * exp_fac(:)
         end if

!-----------------------------------------------------------------
!       ... ch3coch3 + oh -> ro2 + h2o
!-----------------------------------------------------------------
         if( usr22_ndx > 0 ) then
            call vexp( exp_fac, -1320.*tinv, plonl )
            call vexp( xr, 423.*tinv, plonl )
            rxt(:,k,usr22_ndx) = 8.8e-12 * exp_fac(:) + 1.7e-14 * xr(:)
         end if
!-----------------------------------------------------------------
!       ... DMS + OH -> .75 * SO2
!-----------------------------------------------------------------
         if( usr24_ndx > 0 ) then
            call vexp( exp_fac, 7810.*tinv, plonl )
            call vexp( xr, 7460.*tinv, plonl )
            ko(:) = 1. + 5.5e-31 * xr * invariants(:,k,indexm) * 0.21
            ko(:) = 1.7e-42 * exp_fac * invariants(:,k,indexm) * 0.21 / ko(:)
            call vexp( exp_fac, -234.*tinv, plonl )
            rxt(:,k,usr24_ndx) = 0.75*ko(:) + 9.6e-12 * exp_fac
         end if

!-----------------------------------------------------------------
!        ... Cl2O2 + M -> 2*ClO + M
!-----------------------------------------------------------------
         if( strat38_ndx > 0 ) then
            if( strat37_ndx > 0 ) then
               call vexp( exp_fac, -8744.*tinv, plonl )
               rxt(:,k,strat38_ndx) = rxt(:,k,strat37_ndx) * 7.874e26 * exp_fac(:)
            else
               rxt(:,k,strat38_ndx) = 0.
            end if
         end if
#else
!-----------------------------------------------------------------
!        ... n2o5 + m --> no2 + no3 + m
!-----------------------------------------------------------------
         if( usr3_ndx > 0 ) then
            if( usr2_ndx > 0 ) then
               rxt(:,k,usr3_ndx) = rxt(:,k,usr2_ndx) * 3.333e26 * exp( -10991.*tinv(:) )
            else
               rxt(:,k,usr3_ndx) = 0.
            end if
         end if

!-----------------------------------------------------------------
!        set rates for:
!         ... hno3 + oh --> no3 + h2o
!           ho2no2 + m --> ho2 + no2 + m
!           co + oh --> co2 + ho2
!-----------------------------------------------------------------
         if( usr5_ndx > 0 ) then
            ko(:) = m(:,k) * 6.5e-34 * exp( 1335.*tinv(:) )
            ko(:) = ko(:) / (1. + ko(:)/(2.7e-17*exp( 2199.*tinv(:) )))
            rxt(:,k,usr5_ndx) = ko(:) + 2.4e-14*exp( 460.*tinv(:) )
         end if
         if( usr7_ndx > 0 ) then
            if( usr6_ndx > 0 ) then
               rxt(:,k,usr7_ndx) = rxt(:,k,usr6_ndx) * exp( -10900.*tinv(:) )/ 2.1e-27
            else
               rxt(:,k,usr7_ndx) = 0.
            end if
         end if
         if( usr8_ndx > 0 ) then
            rxt(:,k,usr8_ndx) = 1.5e-13 * (1. + 6.e-7*boltz*m(:,k)*temp(:,k))
         end if

!-----------------------------------------------------------------
!        ... ho2 + ho2 --> h2o2
!        note: this rate involves the water vapor number density
!-----------------------------------------------------------------
         if( usr9_ndx > 0 ) then
            if( indexh2o > 0 ) then
               fc(:)   = 1. + 1.4e-21 * invariants(:,k,indexh2o) * exp( 2200.*tinv(:) )
            else if( h2o_ndx > 0 ) then
               fc(:)   = 1. + 1.4e-21 * qin(:,k,h2o_ndx) * invariants(:,k,indexm) * exp( 2200.*tinv(:) )
            else
               fc(:) = 1.
            end if
            ko(:)   = 2.3e-13 * exp( 600.*tinv(:) )
            kinf(:) = 1.7e-33 * m(:,k) * exp( 1000.*tinv(:) )
            rxt(:,k,usr9_ndx) = (ko(:) + kinf(:)) * fc(:)
         end if

!-----------------------------------------------------------------
!            ... mco3 + no2 -> mpan
!-----------------------------------------------------------------
         if( usr14_ndx > 0 ) then
            rxt(:,k,usr14_ndx) = 1.1e-11 * tp(:) / m(:,k)
         end if

!-----------------------------------------------------------------
!        ... pan + m --> ch3co3 + no2 + m
!-----------------------------------------------------------------
         exp_fac(:) = exp( -14000.*tinv(:) )
         if( usr12_ndx > 0 ) then
            if( usr11_ndx > 0 ) then
               rxt(:,k,usr12_ndx) = rxt(:,k,usr11_ndx) * 1.111e28 * exp_fac(:)
            else
               rxt(:,k,usr12_ndx) = 0.
            end if
         end if

!-----------------------------------------------------------------
!        ... mpan + m --> mco3 + no2 + m
!-----------------------------------------------------------------
         if( usr15_ndx > 0 ) then
            if( usr14_ndx > 0 ) then
               rxt(:,k,usr15_ndx) = rxt(:,k,usr14_ndx) * 1.111e28 * exp_fac(:)
            else
               rxt(:,k,usr15_ndx) = 0.
            end if
         end if

!-----------------------------------------------------------------
!       ... xooh + oh -> h2o + oh
!-----------------------------------------------------------------
         if( usr21_ndx > 0 ) then
            rxt(:,k,usr21_ndx) = temp(:,k)**2 * 7.69e-17 * exp( 253.*tinv(:) )
         end if

!-----------------------------------------------------------------
!       ... ch3coch3 + oh -> ro2 + h2o
!-----------------------------------------------------------------
         if( usr22_ndx > 0 ) then
            rxt(:,k,usr22_ndx) = 8.8e-12 * exp( -1320.*tinv(:) ) &
                                 + 1.7e-14 * exp( 423.*tinv(:) )
         end if
!-----------------------------------------------------------------
!       ... DMS + OH -> .75 * SO2
!-----------------------------------------------------------------
         if( usr24_ndx > 0 ) then
            ko(:) = 1. + 5.5e-31 * exp( 7460.*tinv(:) ) * invariants(:,k,indexm) * 0.21
            ko(:) = 1.7e-42 * exp( 7810.*tinv(:) ) * invariants(:,k,indexm) * 0.21 / ko(:)
            rxt(:,k,usr24_ndx) = 0.75*ko(:) + 9.6e-12*exp( -234.*tinv(:) )
         end if

!-----------------------------------------------------------------
!        ... Cl2O2 + M -> 2*ClO + M
!-----------------------------------------------------------------
         if( strat38_ndx > 0 ) then
            if( strat37_ndx > 0 ) then
               rxt(:,k,strat38_ndx) = rxt(:,k,strat37_ndx) * 7.874e26 * exp( -8744.*tinv(:) )
            else
               rxt(:,k,strat38_ndx) = 0.
            end if
         end if
#endif
         if( usr16_ndx > 0 .or. usr17_ndx > 0 .or. usr25_ndx > 0 ) then
!-----------------------------------------------------------------
!         ... n2o5 --> 2*hno3
!             no3 --> hno3
!-----------------------------------------------------------------
!        ... first compute the relative humidity
!-----------------------------------------------------------------
!           call aqsat( temp(1,k), pmid(1,k), satv, satq, plonl, &
!                       plonl, 1, 1, 1 )
!           relhum(:) = .622 * h2ovmr(:,k) / satq(:)
!           relhum(:) = max( 0.,min( 1.,relhum(:) ) )
            call rh_calc( pmid(:,k), temp(:,k), sh(:,k), relhum(:,k) )
!-------------------------------------------------------------------------
!         ... estimate humidity effect on aerosols (from shettle and fenn, 1979)
!           xr is a factor of the increase aerosol radii with hum (hum=0., factor=1)
!-------------------------------------------------------------------------
            xr(:)     = .999151 + relhum(:,k)*(1.90445 + relhum(:,k)*(-6.35204 + relhum(:,k)*5.32061))
!-------------------------------------------------------------------------
!         ... estimate sulfate particles surface area (cm2/cm3) in each grid
!-------------------------------------------------------------------------
            if( so4_ndx > 0 ) then
               sur(:)    = qin(:,k,so4_ndx)
            else
               sur(:)    = sulfate(:,k)
            end if
            sur(:)    = sur(:)*m(:,k)/avo*wso4 &              ! xform mixing ratio to g/cm3
                        / amas &                                    ! xform g/cm3 to num particels/cm3
                        * fare &                                    ! xform num particels/cm3 to cm2/cm3
                        * xr(:)*xr(:)                               ! humidity factor
!-----------------------------------------------------------------
!        ... compute the "aerosol" reaction rates
!-----------------------------------------------------------------
!             k = gam * a * velo/4
!
!       where velo = sqrt[ 8*bk*t/pi/(w/av) ]
!             bk = 1.381e-16
!             av = 6.02e23
!             w  = 108 (n2o5)  ho2(33)  ch2o (30)  nh3(15)  
!
!       so that velo = 1.40e3*sqrt(t)  (n2o5)   gama=0.1
!       so that velo = 2.53e3*sqrt(t)  (ho2)    gama>0.2
!       so that velo = 2.65e3*sqrt(t)  (ch2o)   gama>0.022
!       so that velo = 3.75e3*sqrt(t)  (nh3)    gama=0.4
!--------------------------------------------------------
!           xr(:) = .25 * gam1 * sur(:) * 1.40e3 * sqrt( temp(:,k) )
            xr(:) = 1./(rm1/dg + 4./(gam1+1.e-30)/(1.40e3 * sqrt( temp(:,k))))*sur(:)
            if( usr16_ndx > 0 ) then
               rxt(:,k,usr16_ndx) = xr(:)
            end if
            if( usr17_ndx > 0 ) then
               rxt(:,k,usr17_ndx) = xr(:)
            end if
            if( usr25_ndx > 0 ) then
               rxt(:,k,usr25_ndx) = &
                  1./(rm1/dg + 4./(gam4+1.e-30)/(3.75e3 * sqrt( temp(:,k))))*sur(:)
            end if
         end if

         if( strat72_ndx > 0 .or. strat73_ndx > 0 .or. strat74_ndx > 0 .or. &
             strat75_ndx > 0 .or. strat76_ndx > 0 .or. strat77_ndx > 0 .or. &
             strat78_ndx > 0 .or. strat79_ndx > 0 .or. strat80_ndx > 0 ) then

            if( strat72_ndx > 0 ) then
               rxt(:,k,strat72_ndx) = 0.
            end if
            if( strat73_ndx > 0 ) then
               rxt(:,k,strat73_ndx) = 0.
            end if
            if( strat74_ndx > 0 ) then
               rxt(:,k,strat74_ndx) = 0.
            end if
            if( strat75_ndx > 0 ) then
               rxt(:,k,strat75_ndx) = 0.
            end if
            if( strat76_ndx > 0 ) then
               rxt(:,k,strat76_ndx) = 0.
            end if
            if( strat77_ndx > 0 ) then
               rxt(:,k,strat77_ndx) = 0.
            end if
            if( strat78_ndx > 0 ) then
               rxt(:,k,strat78_ndx) = 0.
            end if
            if( strat79_ndx > 0 ) then
               rxt(:,k,strat79_ndx) = 0.
            end if
            if( strat80_ndx > 0 ) then
               rxt(:,k,strat80_ndx) = 0.
            end if

         end if

      end do

      end subroutine usrrxt

      subroutine rh_calc(pmid, temp, sh, rh)
              
        implicit none
        
        real, intent(in), dimension(:) :: pmid, temp, sh
        real, intent(out), dimension(:) :: rh
        
        real, dimension(size(temp)) :: esat
!-----------------------------------------------------------------------
!       Calculate RELATIVE humidity.
!       This is calculated according to the formula:
!
!       RH   = qv / (epsilon*esat/ [pfull  -  (1.-epsilon)*esat])
!
!       Where epsilon = Rdgas/RVgas = d622
!
!       and where 1- epsilon = d378
!
!       Note that rh does not have its proper value
!       until all of the following code has been executed.  That
!       is, rh is used to store intermediary results
!       in forming the full solution.
!-----------------------------------------------------------------------
        
!-----------------------------------------------------------------------
!calculate water saturated vapor pressure
!-----------------------------------------------------------------------
        call escomp(temp, esat)
        
!-----------------------------------------------------------------------
!calulate denominator in qsat formula
!-----------------------------------------------------------------------
        rh(:) = pmid(:) - d378 * esat(:)
        
!-----------------------------------------------------------------------
!limit denominator to esat, and thus qs to epsilon
!this is done to avoid blow up in the upper stratosphere
!where pfull ~ esat
!-----------------------------------------------------------------------
        rh(:) = MAX(RH(:),esat(:))
        
!-----------------------------------------------------------------------
!calculate rh
!-----------------------------------------------------------------------
        rh(:) = sh(:) / (d622 * esat(:) / rh(:))
        
      end subroutine rh_calc
        
      end module mo_usrrxt_mod
