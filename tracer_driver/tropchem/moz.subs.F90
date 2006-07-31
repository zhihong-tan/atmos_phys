      module mo_setrxt_mod

implicit none
      private
      public :: setrxt

character(len=128), parameter :: version     = '$Id: moz.subs.F90,v 13.0 2006/03/28 21:16:37 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: memphis_2006_07 $'
logical                       :: module_is_initialized = .false.

      contains

!++lwh
      subroutine setrxt( rate, temp, m, plonl, plev, plnplv )
!--lwh

      use chem_mods_mod, only : rxntot
      use mo_jpl_mod,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... Dummy arguments
!-------------------------------------------------------
!++lwh
      integer, intent(in) :: plonl, plev, plnplv
!--lwh
      real, intent(in)    :: temp(plonl,plev), m(plonl,plev)
      real, intent(inout) :: rate(plonl,plev,rxntot)

!-------------------------------------------------------
!       ... Local variables
!-------------------------------------------------------
      real  ::  itemp(plonl,plev), exp_fac(plonl,plev)
      real, dimension(plonl,plev) :: ko, kinf

      rate(:,:,25) = 2.2e-10
      rate(:,:,33) = 0.
      rate(:,:,41) = 1.5e-10
      rate(:,:,48) = 1e-11
      rate(:,:,50) = 1.1e-10
      rate(:,:,66) = 1e-12
      rate(:,:,73) = 2.e-13
      rate(:,:,74) = 6.8e-14
      rate(:,:,93) = 4.e-14
      rate(:,:,94) = 7.1e-6
      rate(:,:,95) = 7.1e-6
      itemp(:,:) = 1. / temp(:,:)
      rate(:,:,22) = 8e-12 * exp( -2060. * itemp(:,:) )
      rate(:,:,23) = 1.8e-11 * exp( 110. * itemp(:,:) )
      rate(:,:,24) = 3.2e-11 * exp( 70. * itemp(:,:) )
      exp_fac(:,:) = exp( 250. * itemp(:,:) )
      rate(:,:,26) = 3.5e-12 * exp_fac(:,:)
      rate(:,:,57) = 4.8e-11 * exp_fac(:,:)
      rate(:,:,27) = 3e-12 * exp( -1500. * itemp(:,:) )
      exp_fac(:,:) = exp( 180. * itemp(:,:) )
      rate(:,:,28) = 5.6e-12 * exp_fac(:,:)
      rate(:,:,77) = 2.2e-12 * exp_fac(:,:)
      rate(:,:,82) = 4.2e-12 * exp_fac(:,:)
      rate(:,:,87) = 4.2e-12 * exp_fac(:,:)
      rate(:,:,29) = 1.2e-13 * exp( -2450. * itemp(:,:) )
      exp_fac(:,:) = exp( 170. * itemp(:,:) )
      rate(:,:,30) = 2.3e-12 * exp_fac(:,:)
      rate(:,:,36) = 1.5e-11 * exp_fac(:,:)
      rate(:,:,38) = 1.3e-12 * exp( 380. * itemp(:,:) )
      rate(:,:,40) = 2.45e-12 * exp( -1775. * itemp(:,:) )
      rate(:,:,42) = 3.e-12 * exp( 280. * itemp(:,:) )
      rate(:,:,43) = 5.e-13 * exp( -424. * itemp(:,:) )
      rate(:,:,44) = 1.9e-14 * exp( 706. * itemp(:,:) )
      rate(:,:,45) = 3.8e-13 * exp( 800. * itemp(:,:) )
      exp_fac(:,:) = exp( 200. * itemp(:,:) )
      rate(:,:,46) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,52) = 3e-11 * exp_fac(:,:)
      rate(:,:,75) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,85) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,89) = 3.8e-12 * exp_fac(:,:)
      rate(:,:,47) = 6.0e-13 * exp( -2058. * itemp(:,:) )
      rate(:,:,51) = 2.2e-11 * exp( 120. * itemp(:,:) )
      rate(:,:,53) = 1.5e-12 * exp( -880. * itemp(:,:) )
      rate(:,:,54) = 2e-14 * exp( -680. * itemp(:,:) )
      rate(:,:,56) = 2.9e-12 * exp( -160. * itemp(:,:) )
      rate(:,:,58) = 4.2e-12 * exp( -240. * itemp(:,:) )
      exp_fac(:,:) = exp( -2000. * itemp(:,:) )
      rate(:,:,59) = 5.5e-12 * exp_fac(:,:)
      rate(:,:,69) = 1.05e-14 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270. * itemp(:,:) )
      rate(:,:,60) = 5.6e-12 * exp_fac(:,:)
      rate(:,:,62) = 8.1e-12 * exp_fac(:,:)
      rate(:,:,61) = 1.4e-12 * exp( -1900. * itemp(:,:) )
      rate(:,:,64) = 4.3e-13 * exp( 1040. * itemp(:,:) )
      rate(:,:,65) = 1.3e-12 * exp( 640. * itemp(:,:) )
      rate(:,:,68) = 2.5e-12 * exp( 500. * itemp(:,:) )
      rate(:,:,70) = 8.7e-12 * exp( -1070. * itemp(:,:) )
      rate(:,:,71) = 2.6e-12 * exp( 365. * itemp(:,:) )
      exp_fac(:,:) = exp( 700. * itemp(:,:) )
      rate(:,:,72) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,78) = 8.e-13 * exp_fac(:,:)
      rate(:,:,83) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,88) = 7.5e-13 * exp_fac(:,:)
      rate(:,:,76) = 2.54e-11 * exp( 410. * itemp(:,:) )
      rate(:,:,81) = 1.0e-11 * exp( -660. * itemp(:,:) )
      rate(:,:,84) = 3.75e-13 * exp( -40. * itemp(:,:) )
      rate(:,:,90) = 3.03e-12 * exp( -446. * itemp(:,:) )
      rate(:,:,91) = 8.4e-13 * exp( 830. * itemp(:,:) )
      rate(:,:,92) = 1.4e-12 * exp( -1860. * itemp(:,:) )
      rate(:,:,98) = 1.9e-13 * exp( 520. * itemp(:,:) )
      rate(:,:,100) = 1.7e-12 * exp( -710. * itemp(:,:) )

      itemp(:,:) = 300. * itemp(:,:)

      ko(:,:) = 2.e-30 * itemp(:,:)**4.4
      kinf(:,:) = 1.4e-12 * itemp(:,:)**.7
      call jpl( rate(1,1,31), m, .6, ko, kinf, plnplv )

      ko(:,:) = 2.4e-30 * itemp(:,:)**3.1
      kinf(:,:) = 1.7e-11 * itemp(:,:)**2.1
      call jpl( rate(1,1,34), m, .6, ko, kinf, plnplv )

      ko(:,:) = 1.8e-31 * itemp(:,:)**3.2
      kinf(:,:) = 4.7e-12 * itemp(:,:)**1.4
      call jpl( rate(1,1,37), m, .6, ko, kinf, plnplv )

      ko(:,:) = 8.5e-29 * itemp(:,:)**6.5
      kinf(:,:) = 1.1e-11 * itemp(:,:)
      call jpl( rate(1,1,63), m, .6, ko, kinf, plnplv )

      ko(:,:) = 3.e-31 * itemp(:,:)**3.3
      kinf(:,:) = 1.5e-12
      call jpl( rate(1,1,96), m, 0.6, ko, kinf, plnplv )

      end subroutine setrxt

      end module mo_setrxt_mod

      module mo_adjrxt_mod

      private
      public :: adjrxt

      contains

      subroutine adjrxt( rate, inv, m, plnplv )

      use chem_mods_mod, only : nfs, rxntot

      implicit none

!--------------------------------------------------------------------
!       ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real, intent(in)    :: inv(plnplv,nfs)
      real, intent(in)    :: m(plnplv)
      real, intent(inout) :: rate(plnplv,rxntot)

!--------------------------------------------------------------------
!       ... Local variables
!--------------------------------------------------------------------
      real    ::  im(plnplv)

      rate(:, 23) = rate(:, 23) * inv(:, 2)
      rate(:, 24) = rate(:, 24) * inv(:, 3)
      rate(:, 25) = rate(:, 25) * inv(:, 4)
      rate(:, 31) = rate(:, 31) * inv(:, 1)
      rate(:, 32) = rate(:, 32) * inv(:, 1)
      rate(:, 33) = rate(:, 33) * inv(:, 4)
      rate(:, 34) = rate(:, 34) * inv(:, 1)
      rate(:, 37) = rate(:, 37) * inv(:, 1)
      rate(:, 39) = rate(:, 39) * inv(:, 1)
      rate(:, 40) = rate(:, 40) * inv(:, 5)
      rate(:, 41) = rate(:, 41) * inv(:, 5)
      rate(:, 50) = rate(:, 50) * inv(:, 8)
      rate(:, 59) = rate(:, 59) * inv(:, 8)
      rate(:, 63) = rate(:, 63) * inv(:, 1)
      rate(:, 67) = rate(:, 67) * inv(:, 1)
      rate(:, 70) = rate(:, 70) * inv(:, 6)
      rate(:, 81) = rate(:, 81) * inv(:, 7)
      rate(:, 96) = rate(:, 96) * inv(:, 1)
      rate(:, 21) = rate(:, 21) * inv(:, 3) * inv(:, 1)
      rate(:, 22) = rate(:, 22) * m(:)
      rate(:, 26) = rate(:, 26) * m(:)
      rate(:, 27) = rate(:, 27) * m(:)
      rate(:, 28) = rate(:, 28) * m(:)
      rate(:, 29) = rate(:, 29) * m(:)
      rate(:, 30) = rate(:, 30) * m(:)
      rate(:, 31) = rate(:, 31) * m(:)
      rate(:, 34) = rate(:, 34) * m(:)
      rate(:, 35) = rate(:, 35) * m(:)
      rate(:, 36) = rate(:, 36) * m(:)
      rate(:, 37) = rate(:, 37) * m(:)
      rate(:, 38) = rate(:, 38) * m(:)
      rate(:, 42) = rate(:, 42) * m(:)
      rate(:, 43) = rate(:, 43) * m(:)
      rate(:, 44) = rate(:, 44) * m(:)
      rate(:, 45) = rate(:, 45) * m(:)
      rate(:, 46) = rate(:, 46) * m(:)
      rate(:, 47) = rate(:, 47) * m(:)
      rate(:, 48) = rate(:, 48) * m(:)
      rate(:, 49) = rate(:, 49) * m(:)
      rate(:, 51) = rate(:, 51) * m(:)
      rate(:, 52) = rate(:, 52) * m(:)
      rate(:, 53) = rate(:, 53) * m(:)
      rate(:, 54) = rate(:, 54) * m(:)
      rate(:, 55) = rate(:, 55) * m(:)
      rate(:, 56) = rate(:, 56) * m(:)
      rate(:, 57) = rate(:, 57) * m(:)
      rate(:, 58) = rate(:, 58) * m(:)
      rate(:, 60) = rate(:, 60) * m(:)
      rate(:, 61) = rate(:, 61) * m(:)
      rate(:, 62) = rate(:, 62) * m(:)
      rate(:, 63) = rate(:, 63) * m(:)
      rate(:, 64) = rate(:, 64) * m(:)
      rate(:, 65) = rate(:, 65) * m(:)
      rate(:, 66) = rate(:, 66) * m(:)
      rate(:, 68) = rate(:, 68) * m(:)
      rate(:, 69) = rate(:, 69) * m(:)
      rate(:, 71) = rate(:, 71) * m(:)
      rate(:, 72) = rate(:, 72) * m(:)
      rate(:, 73) = rate(:, 73) * m(:)
      rate(:, 74) = rate(:, 74) * m(:)
      rate(:, 75) = rate(:, 75) * m(:)
      rate(:, 76) = rate(:, 76) * m(:)
      rate(:, 77) = rate(:, 77) * m(:)
      rate(:, 78) = rate(:, 78) * m(:)
      rate(:, 82) = rate(:, 82) * m(:)
      rate(:, 83) = rate(:, 83) * m(:)
      rate(:, 84) = rate(:, 84) * m(:)
      rate(:, 85) = rate(:, 85) * m(:)
      rate(:, 86) = rate(:, 86) * m(:)
      rate(:, 87) = rate(:, 87) * m(:)
      rate(:, 88) = rate(:, 88) * m(:)
      rate(:, 89) = rate(:, 89) * m(:)
      rate(:, 90) = rate(:, 90) * m(:)
      rate(:, 91) = rate(:, 91) * m(:)
      rate(:, 92) = rate(:, 92) * m(:)
      rate(:, 93) = rate(:, 93) * m(:)
      rate(:, 96) = rate(:, 96) * m(:)
      rate(:, 97) = rate(:, 97) * m(:)
      rate(:, 98) = rate(:, 98) * m(:)
      rate(:,100) = rate(:,100) * m(:)

      end subroutine adjrxt

      end module mo_adjrxt_mod

      module mo_phtadj_mod

      private
      public :: phtadj

      contains

      subroutine phtadj( p_rate, inv, m, plnplv )

      use chem_mods_mod, only : nfs, phtcnt

      implicit none

!--------------------------------------------------------------------
!       ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real, intent(in)    :: inv(plnplv,nfs)
      real, intent(in)    :: m(plnplv)
      real, intent(inout) :: p_rate(plnplv,phtcnt)

!--------------------------------------------------------------------
!       ... Local variables
!--------------------------------------------------------------------
      real    ::  im(plnplv)

      im(:) = 1. / m(:)
      p_rate(:,  1) = p_rate(:,  1)  * inv(:, 3) * im(:)

      end subroutine phtadj

      end module mo_phtadj_mod

      module mo_rxt_mod

      private
      public :: rxt_mod

      contains

      subroutine rxt_mod( rate, het_rates, grp_ratios, plnplv )

      use chem_mods_mod, only : rxntot, hetcnt, grpcnt

      implicit none

!---------------------------------------------------------------------------
!       ... Dummy arguments
!---------------------------------------------------------------------------
      integer, intent(in) ::  plnplv
      real, intent(inout) ::  rate(plnplv,rxntot)
      real, intent(inout) ::  het_rates(plnplv,hetcnt)
      real, intent(in)    ::  grp_ratios(plnplv,grpcnt)


      end subroutine rxt_mod

      end module mo_rxt_mod

      module mo_make_grp_vmr_mod

      private
      public :: mak_grp_vmr

      contains

      subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs, plonl )

      use mo_grid_mod,   only : plev, pcnstm1
      use chem_mods_mod, only : grpcnt

      implicit none

!----------------------------------------------------------------------------
!        ... Dummy arguments
!----------------------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in)    :: vmr(plonl,plev,pcnstm1)
      real, intent(in)    :: group_ratios(plonl,plev,grpcnt)
      real, intent(out)   :: group_vmrs(plonl,plev,grpcnt)

!----------------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------------
      integer ::  k

      end subroutine mak_grp_vmr

      end module mo_make_grp_vmr_mod
