      module mo_exp_prod_loss_mod

implicit none
character(len=128), parameter :: version     = '$Id: moz.mat.F90,v 13.0 2006/03/28 21:16:31 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: memphis_2006_07 $'
logical                       :: module_is_initialized = .false.

      contains

      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )

      use chem_mods_mod, only : clscnt1, rxntot, hetcnt
      use mo_grid_mod,   only : pcnstm1

      implicit none

!--------------------------------------------------------------------
!     ... Dummy args                                                                      
!--------------------------------------------------------------------
      real, dimension(:,:), intent(out) :: &
            prod, &
            loss
      real, intent(in)    ::  y(:,:)
      real, intent(in)    ::  rxt(:,:)
      real, intent(in)    ::  het_rates(:,:)


!--------------------------------------------------------------------
!       ... Loss and production for Explicit method
!--------------------------------------------------------------------

      loss(:,1) = (rxt(:,49)* y(:,14))* y(:,13)
      prod(:,1) = 0.
      loss(:,2) = ( + rxt(:,94))* y(:,32)
      prod(:,2) = 0.
      loss(:,3) = 0.
      prod(:,3) =rxt(:,94)*y(:,32)
      loss(:,4) = ( + rxt(:,95))* y(:,34)
      prod(:,4) = 0.
      loss(:,5) = 0.
      prod(:,5) =rxt(:,95)*y(:,34)

      end subroutine exp_prod_loss

      end module mo_exp_prod_loss_mod

      module mo_imp_prod_loss_mod

      contains

      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )

      use chem_mods_mod, only : clscnt4, rxntot, hetcnt, clsze
      use mo_grid_mod,   only : pcnstm1

      implicit none

!--------------------------------------------------------------------
!     ... Dummy args                                                                      
!--------------------------------------------------------------------
      real, dimension(:,:), intent(out) :: &
            prod, &
            loss
      real, intent(in)    ::  y(:,:)
      real, intent(in)    ::  rxt(:,:)
      real, intent(in)    ::  het_rates(:,:)


!--------------------------------------------------------------------
!     ... Local variables                                                                 
!--------------------------------------------------------------------
      integer :: k

!--------------------------------------------------------------------
!       ... Loss and production for Implicit method
!--------------------------------------------------------------------

      do k = 1,clsze
         loss(k,28) = (rxt(k,22)* y(k,2) +rxt(k,27)* y(k,4) +rxt(k,29)* y(k,5) &
                  +rxt(k,53)* y(k,14) +rxt(k,54)* y(k,15) +rxt(k,69)* y(k,18) &
                  + rxt(k,2) + rxt(k,3))* y(k,1)
         prod(k,28) =rxt(k,21)*y(k,2) +.890*rxt(k,7)*y(k,6) +.300*rxt(k,64)*y(k,20) &
                 *y(k,15)
         loss(k,7) = ( + rxt(k,23) + rxt(k,24) + rxt(k,25) + rxt(k,41) + rxt(k,50)) &
                 * y(k,3)
         prod(k,7) =rxt(k,2)*y(k,1)
         loss(k,21) = (rxt(k,22)* y(k,1) +rxt(k,28)* y(k,5) +rxt(k,51)* y(k,14) &
                  +rxt(k,52)* y(k,15) + rxt(k,21))* y(k,2)
         prod(k,21) = (rxt(k,23) +rxt(k,24))*y(k,3) +rxt(k,3)*y(k,1) +rxt(k,4)*y(k,5) &
                  +rxt(k,58)*y(k,14)*y(k,14)
         loss(k,33) = (rxt(k,27)* y(k,1) +rxt(k,36)* y(k,6) +rxt(k,42)* y(k,10) &
                  +rxt(k,26)* y(k,15) +rxt(k,62)* y(k,20) +rxt(k,77)* y(k,23) &
                  +rxt(k,71)* y(k,24) +rxt(k,82)* y(k,26) +rxt(k,87)* y(k,30))* y(k,4)
         prod(k,33) = (rxt(k,4) +rxt(k,28)*y(k,2))*y(k,5) +.110*rxt(k,7)*y(k,6)
         loss(k,30) = (rxt(k,29)* y(k,1) +rxt(k,28)* y(k,2) +rxt(k,31)* y(k,6) &
                  +rxt(k,34)* y(k,14) +rxt(k,37)* y(k,15) +rxt(k,63)* y(k,20) &
                  + rxt(k,4))* y(k,5)
         prod(k,30) = (rxt(k,26)*y(k,15) +rxt(k,27)*y(k,1) +2.000*rxt(k,36)*y(k,6) + &
                 rxt(k,42)*y(k,10) +rxt(k,62)*y(k,20) +rxt(k,71)*y(k,24) + &
                 rxt(k,77)*y(k,23) +rxt(k,82)*y(k,26) +rxt(k,87)*y(k,30))*y(k,4) &
                  + (.890*rxt(k,7) +rxt(k,30)*y(k,15) +rxt(k,90)*y(k,18) + &
                 rxt(k,98)*y(k,38))*y(k,6) + (rxt(k,8) +rxt(k,39) +rxt(k,38)*y(k,14)) &
                 *y(k,8) + (rxt(k,5) +rxt(k,32))*y(k,9) + (.600*rxt(k,15) +rxt(k,67)) &
                 *y(k,22) +rxt(k,6)*y(k,7)
         loss(k,34) = (rxt(k,36)* y(k,4) +rxt(k,31)* y(k,5) +rxt(k,47)* y(k,12) &
                  +rxt(k,30)* y(k,15) +rxt(k,90)* y(k,18) +rxt(k,61)* y(k,19) &
                  +rxt(k,92)* y(k,31) +rxt(k,98)* y(k,38) + rxt(k,7) + rxt(k,80)) &
                 * y(k,6)
         prod(k,34) = (rxt(k,5) +rxt(k,32))*y(k,9) + (rxt(k,35)*y(k,7) + &
                 rxt(k,93)*y(k,22))*y(k,14) +rxt(k,29)*y(k,5)*y(k,1) +.400*rxt(k,15) &
                 *y(k,22)
         loss(k,17) = (rxt(k,35)* y(k,14) + rxt(k,6))* y(k,7)
         prod(k,17) = (rxt(k,80) +rxt(k,47)*y(k,12) +rxt(k,61)*y(k,19) + &
                 rxt(k,92)*y(k,31))*y(k,6) + (2.000*rxt(k,33) +2.000*rxt(k,79))*y(k,9) &
                  +rxt(k,34)*y(k,14)*y(k,5)
         loss(k,10) = (rxt(k,38)* y(k,14) + rxt(k,8) + rxt(k,39))* y(k,8)
         prod(k,10) =rxt(k,37)*y(k,15)*y(k,5)
         loss(k,8) = ( + rxt(k,5) + rxt(k,32) + rxt(k,33) + rxt(k,79))* y(k,9)
         prod(k,8) =rxt(k,31)*y(k,6)*y(k,5)
         loss(k,32) = (rxt(k,42)* y(k,4) + 2.*(rxt(k,43) +rxt(k,44))* y(k,10) &
                  +rxt(k,45)* y(k,15) +rxt(k,65)* y(k,20) +rxt(k,73)* y(k,24) &
                  +rxt(k,84)* y(k,26))* y(k,10)
         prod(k,32) = (rxt(k,62)*y(k,4) +.900*rxt(k,65)*y(k,10) + &
                 2.000*rxt(k,68)*y(k,20))*y(k,20) + (rxt(k,40) + &
                 .700*rxt(k,46)*y(k,11))*y(k,14) +.750*rxt(k,41)*y(k,3) +rxt(k,13) &
                 *y(k,19) +rxt(k,14)*y(k,21) +.400*rxt(k,15)*y(k,22) +rxt(k,19) &
                 *y(k,28)
         loss(k,16) = (rxt(k,46)* y(k,14) + rxt(k,9))* y(k,11)
         prod(k,16) =rxt(k,45)*y(k,15)*y(k,10)
         loss(k,27) = (rxt(k,47)* y(k,6) +rxt(k,48)* y(k,14) + rxt(k,10) + rxt(k,11)) &
                 * y(k,12)
         prod(k,27) = (rxt(k,42)*y(k,4) +2.000*rxt(k,43)*y(k,10) +rxt(k,44)*y(k,10) + &
                 rxt(k,65)*y(k,20) +.700*rxt(k,73)*y(k,24) +rxt(k,84)*y(k,26))*y(k,10) &
                  + (.300*rxt(k,46)*y(k,11) +.500*rxt(k,66)*y(k,21) + &
                 rxt(k,93)*y(k,22))*y(k,14) +.250*rxt(k,41)*y(k,3) +rxt(k,87)*y(k,30) &
                 *y(k,4) +rxt(k,9)*y(k,11) +rxt(k,18)*y(k,29)
         loss(k,29) = (rxt(k,53)* y(k,1) +rxt(k,51)* y(k,2) +rxt(k,34)* y(k,5) &
                  +rxt(k,35)* y(k,7) +rxt(k,38)* y(k,8) +rxt(k,46)* y(k,11) +rxt(k,48) &
                 * y(k,12) +rxt(k,49)* y(k,13) + 2.*rxt(k,58)* y(k,14) +rxt(k,57) &
                 * y(k,15) +rxt(k,56)* y(k,16) +rxt(k,76)* y(k,18) +rxt(k,60)* y(k,19) &
                  +rxt(k,66)* y(k,21) +rxt(k,93)* y(k,22) +rxt(k,75)* y(k,25) &
                  +rxt(k,85)* y(k,27) +rxt(k,86)* y(k,28) +rxt(k,89)* y(k,29) &
                  +rxt(k,91)* y(k,31) +rxt(k,96)* y(k,36) +rxt(k,97)* y(k,38) &
                  +rxt(k,100)* y(k,39) + rxt(k,40) + rxt(k,59) + rxt(k,70) &
                  + rxt(k,81))* y(k,14)
         prod(k,29) = (rxt(k,26)*y(k,4) +rxt(k,30)*y(k,6) +rxt(k,52)*y(k,2) + &
                 rxt(k,54)*y(k,1))*y(k,15) + (2.000*rxt(k,25) +.750*rxt(k,41) + &
                 rxt(k,50))*y(k,3) + (rxt(k,9) +.300*rxt(k,46)*y(k,14))*y(k,11) &
                  + (rxt(k,16) +.500*rxt(k,75)*y(k,14))*y(k,25) +rxt(k,6)*y(k,7) &
                  +2.000*rxt(k,12)*y(k,16) +rxt(k,14)*y(k,21) +rxt(k,17)*y(k,27) &
                  +rxt(k,18)*y(k,29)
         loss(k,31) = (rxt(k,54)* y(k,1) +rxt(k,52)* y(k,2) +rxt(k,26)* y(k,4) &
                  +rxt(k,37)* y(k,5) +rxt(k,30)* y(k,6) +rxt(k,45)* y(k,10) +rxt(k,57) &
                 * y(k,14) + 2.*rxt(k,55)* y(k,15) +rxt(k,64)* y(k,20) +rxt(k,78) &
                 * y(k,23) +rxt(k,72)* y(k,24) +rxt(k,83)* y(k,26) +rxt(k,88)* y(k,30) &
                 )* y(k,15)
         prod(k,31) = (rxt(k,59) +rxt(k,49)*y(k,13) +rxt(k,48)*y(k,12) + &
                 rxt(k,51)*y(k,2) +rxt(k,53)*y(k,1) +rxt(k,56)*y(k,16))*y(k,14) &
                  + (rxt(k,42)*y(k,4) +2.000*rxt(k,43)*y(k,10) + &
                 .900*rxt(k,65)*y(k,20) +rxt(k,73)*y(k,24) +rxt(k,84)*y(k,26))*y(k,10) &
                  + (rxt(k,71)*y(k,24) +rxt(k,77)*y(k,23) +rxt(k,82)*y(k,26))*y(k,4) &
                  + (.400*rxt(k,41) +rxt(k,50))*y(k,3) + (rxt(k,8) +rxt(k,39))*y(k,8) &
                  + (2.000*rxt(k,10) +rxt(k,47)*y(k,6))*y(k,12) +rxt(k,69)*y(k,18) &
                 *y(k,1) +rxt(k,9)*y(k,11) +rxt(k,13)*y(k,19) +1.200*rxt(k,74)*y(k,24) &
                 *y(k,24) +rxt(k,16)*y(k,25) +rxt(k,17)*y(k,27) +rxt(k,20)*y(k,31)
         loss(k,6) = (rxt(k,56)* y(k,14) + rxt(k,12))* y(k,16)
         prod(k,6) =rxt(k,55)*y(k,15)*y(k,15)
         loss(k,22) = (rxt(k,69)* y(k,1) +rxt(k,90)* y(k,6) +rxt(k,76)* y(k,14)) &
                 * y(k,18)
         prod(k,22) = 0.
         loss(k,26) = (rxt(k,61)* y(k,6) +rxt(k,60)* y(k,14) + rxt(k,13))* y(k,19)
         prod(k,26) = (rxt(k,71)*y(k,4) +.800*rxt(k,73)*y(k,10) + &
                 1.600*rxt(k,74)*y(k,24))*y(k,24) + (rxt(k,16) + &
                 .500*rxt(k,75)*y(k,14))*y(k,25) +.270*rxt(k,82)*y(k,26)*y(k,4)
         loss(k,35) = (rxt(k,62)* y(k,4) +rxt(k,63)* y(k,5) +rxt(k,65)* y(k,10) &
                  +rxt(k,64)* y(k,15) + 2.*rxt(k,68)* y(k,20))* y(k,20)
         prod(k,35) = (rxt(k,61)*y(k,19) +rxt(k,90)*y(k,18) +rxt(k,92)*y(k,31))*y(k,6) &
                  + (rxt(k,60)*y(k,19) +.500*rxt(k,66)*y(k,21) +rxt(k,91)*y(k,31)) &
                 *y(k,14) + (rxt(k,77)*y(k,23) +rxt(k,87)*y(k,30))*y(k,4) &
                  + (.600*rxt(k,15) +rxt(k,67))*y(k,22) +rxt(k,69)*y(k,18)*y(k,1) &
                  +rxt(k,19)*y(k,28) +rxt(k,18)*y(k,29) +rxt(k,20)*y(k,31)
         loss(k,12) = (rxt(k,66)* y(k,14) + rxt(k,14))* y(k,21)
         prod(k,12) =.700*rxt(k,64)*y(k,20)*y(k,15)
         loss(k,19) = (rxt(k,93)* y(k,14) + rxt(k,15) + rxt(k,67))* y(k,22)
         prod(k,19) =rxt(k,63)*y(k,20)*y(k,5)
         loss(k,18) = (rxt(k,77)* y(k,4) +rxt(k,78)* y(k,15))* y(k,23)
         prod(k,18) =rxt(k,76)*y(k,18)*y(k,14)
         loss(k,23) = (rxt(k,71)* y(k,4) +rxt(k,73)* y(k,10) +rxt(k,72)* y(k,15) &
                  + 2.*rxt(k,74)* y(k,24))* y(k,24)
         prod(k,23) = (rxt(k,70) +.500*rxt(k,75)*y(k,25))*y(k,14)
         loss(k,13) = (rxt(k,75)* y(k,14) + rxt(k,16))* y(k,25)
         prod(k,13) =rxt(k,72)*y(k,24)*y(k,15)
         loss(k,25) = (rxt(k,82)* y(k,4) +rxt(k,84)* y(k,10) +rxt(k,83)* y(k,15)) &
                 * y(k,26)
         prod(k,25) = (rxt(k,81) +rxt(k,85)*y(k,27))*y(k,14)
         loss(k,14) = (rxt(k,85)* y(k,14) + rxt(k,17))* y(k,27)
         prod(k,14) =rxt(k,83)*y(k,26)*y(k,15)
         loss(k,20) = (rxt(k,86)* y(k,14) + rxt(k,19))* y(k,28)
         prod(k,20) = (.820*rxt(k,82)*y(k,4) +.820*rxt(k,84)*y(k,10))*y(k,26) &
                  +.820*rxt(k,17)*y(k,27)
         loss(k,15) = (rxt(k,89)* y(k,14) + rxt(k,18))* y(k,29)
         prod(k,15) =rxt(k,88)*y(k,30)*y(k,15)
         loss(k,24) = (rxt(k,87)* y(k,4) +rxt(k,88)* y(k,15))* y(k,30)
         prod(k,24) = (rxt(k,86)*y(k,28) +rxt(k,89)*y(k,29))*y(k,14)
         loss(k,11) = (rxt(k,92)* y(k,6) +rxt(k,91)* y(k,14) + rxt(k,20))* y(k,31)
         prod(k,11) = 0.
         loss(k,5) = (rxt(k,96)* y(k,14))* y(k,36)
         prod(k,5) = (rxt(k,97)*y(k,14) +rxt(k,98)*y(k,6))*y(k,38)
         loss(k,1) = 0.
         prod(k,1) =rxt(k,96)*y(k,36)*y(k,14)
         loss(k,9) = (rxt(k,98)* y(k,6) +rxt(k,97)* y(k,14))* y(k,38)
         prod(k,9) = 0.
         loss(k,4) = (rxt(k,100)* y(k,14) + rxt(k,99))* y(k,39)
         prod(k,4) = 0.
         loss(k,2) = 0.
         prod(k,2) = 0.
         loss(k,3) = 0.
         prod(k,3) =rxt(k,99)*y(k,39)
      end do

      end subroutine imp_prod_loss

      end module mo_imp_prod_loss_mod

      module mo_rodas_prod_loss_mod

      contains

      subroutine rodas_prod_loss( prod, loss, y, rxt, het_rates )

      use chem_mods_mod, only : clscnt5, rxntot, hetcnt, clsze
      use mo_grid_mod,   only : pcnstm1

      implicit none

!--------------------------------------------------------------------
!     ... Dummy args                                                                      
!--------------------------------------------------------------------
      real, dimension(:,:), intent(out) :: &
            prod, &
            loss
      real, intent(in)    ::  y(:,:)
      real, intent(in)    ::  rxt(:,:)
      real, intent(in)    ::  het_rates(:,:)


      end subroutine rodas_prod_loss

      end module mo_rodas_prod_loss_mod

      module mo_indprd_mod

      private
      public :: indprd

      contains

      subroutine indprd( class, prod, y, extfrc, rxt )

      implicit none

!--------------------------------------------------------------------
!       ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      real, intent(in)    :: y(:,:)
      real, intent(in)    :: rxt(:,:)
      real, intent(in)    :: extfrc(:,:)
      real, intent(inout) :: prod(:,:)

!--------------------------------------------------------------------
!       ... "Independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then
         prod(:,1) = (rxt(:,10) +rxt(:,11) +rxt(:,47)*y(:,6) +rxt(:,48)*y(:,14)) &
                 *y(:,12) + (rxt(:,20) +rxt(:,91)*y(:,14) +rxt(:,92)*y(:,6))*y(:,31) &
                  + (2.000*rxt(:,77)*y(:,4) +rxt(:,78)*y(:,15))*y(:,23) +rxt(:,13) &
                 *y(:,19)
                                                                                          
         prod(:,2) = 0.
                                                                                          
         prod(:,3) = 0.
                                                                                          
         prod(:,4) = 0.
                                                                                          
         prod(:,5) = 0.
                                                                                          
!--------------------------------------------------------------------
!       ... "Independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,28) = 0.
                                                                                          
         prod(:,7) = 0.
                                                                                          
         prod(:,21) =2.000*rxt(:,1)
                                                                                          
         prod(:,33) = 0.
                                                                                          
         prod(:,30) = 0.
                                                                                          
         prod(:,34) = 0.
                                                                                          
         prod(:,17) = 0.
                                                                                          
         prod(:,10) = 0.
                                                                                          
         prod(:,8) = 0.
                                                                                          
         prod(:,32) = 0.
                                                                                          
         prod(:,16) = 0.
                                                                                          
         prod(:,27) = 0.
                                                                                          
         prod(:,29) = 0.
                                                                                          
         prod(:,31) = 0.
                                                                                          
         prod(:,6) = 0.
                                                                                          
         prod(:,22) = 0.
                                                                                          
         prod(:,26) = 0.
                                                                                          
         prod(:,35) = 0.
                                                                                          
         prod(:,12) = 0.
                                                                                          
         prod(:,19) = 0.
                                                                                          
         prod(:,18) = 0.
                                                                                          
         prod(:,23) = 0.
                                                                                          
         prod(:,13) = 0.
                                                                                          
         prod(:,25) = 0.
                                                                                          
         prod(:,14) = 0.
                                                                                          
         prod(:,20) = 0.
                                                                                          
         prod(:,15) = 0.
                                                                                          
         prod(:,24) = 0.
                                                                                          
         prod(:,11) = 0.
                                                                                          
         prod(:,5) = 0.
                                                                                          
         prod(:,1) = 0.
                                                                                          
         prod(:,9) = 0.
                                                                                          
         prod(:,4) = 0.
                                                                                          
         prod(:,2) = 0.
                                                                                          
         prod(:,3) = 0.
                                                                                          
      end if                                                                              
                                                                                          
      end subroutine INDPRD                                                               
                                                                                          
      end module MO_INDPRD_MOD                                                                

      module MO_IMP_LIN_MATRIX_MOD

      CONTAINS

      subroutine IMP_LINMAT01( mat, y, rxt, het_rates )
!----------------------------------------------
!       ... Linear Matrix entries for Implicit species
!----------------------------------------------

      use MO_GRID_MOD,   only : pcnstm1
      use CHEM_MODS_MOD, only : rxntot, hetcnt, imp_nzcnt, clsze

      implicit none

!----------------------------------------------
!       ... Dummy args
!----------------------------------------------
      real, intent(in)    ::  y(clsze,pcnstm1)
      real, intent(in)    ::  rxt(clsze,rxntot)
      real, intent(in)    ::  het_rates(clsze,hetcnt)
      real, intent(inout) ::  mat(clsze,imp_nzcnt)
!----------------------------------------------
!       ... Local variables
!----------------------------------------------
      integer :: k


      do k = 1,clsze

         mat(k,146) = -( rxt(k,2) + rxt(k,3) )
         mat(k,259) = .890*rxt(k,7)
         mat(k,85) = rxt(k,21)

         mat(k,13) = -( rxt(k,23) + rxt(k,24) + rxt(k,25) + rxt(k,41) + rxt(k,50) )
         mat(k,142) = rxt(k,2)

         mat(k,84) = -( rxt(k,21) )
         mat(k,143) = rxt(k,3)
         mat(k,189) = rxt(k,4)
         mat(k,14) = rxt(k,23) + rxt(k,24)

         mat(k,196) = rxt(k,4)
         mat(k,264) = .110*rxt(k,7)

         mat(k,193) = -( rxt(k,4) )
         mat(k,21) = rxt(k,5) + rxt(k,32)
         mat(k,65) = rxt(k,6)
         mat(k,261) = .890*rxt(k,7)
         mat(k,30) = rxt(k,8) + rxt(k,39)
         mat(k,75) = .600*rxt(k,15) + rxt(k,67)

         mat(k,265) = -( rxt(k,7) + rxt(k,80) )
         mat(k,22) = rxt(k,5) + rxt(k,32)
         mat(k,77) = .400*rxt(k,15)

         mat(k,63) = -( rxt(k,6) )
         mat(k,20) = 2.000*rxt(k,33) + 2.000*rxt(k,79)
         mat(k,255) = rxt(k,80)

         mat(k,28) = -( rxt(k,8) + rxt(k,39) )

         mat(k,19) = -( rxt(k,5) + rxt(k,32) + rxt(k,33) + rxt(k,79) )

         mat(k,133) = rxt(k,13)
         mat(k,41) = rxt(k,14)
         mat(k,76) = .400*rxt(k,15)
         mat(k,82) = rxt(k,19)
         mat(k,181) = rxt(k,40)
         mat(k,18) = .750*rxt(k,41)

         mat(k,58) = -( rxt(k,9) )

         mat(k,137) = -( rxt(k,10) + rxt(k,11) )
         mat(k,59) = rxt(k,9)
         mat(k,55) = rxt(k,18)
         mat(k,15) = .250*rxt(k,41)

         mat(k,178) = -( rxt(k,40) + rxt(k,59) + rxt(k,70) + rxt(k,81) + rxt(k,49)*y(k,13))
         mat(k,64) = rxt(k,6)
         mat(k,60) = rxt(k,9)
         mat(k,11) = 2.000*rxt(k,12)
         mat(k,40) = rxt(k,14)
         mat(k,46) = rxt(k,16)
         mat(k,51) = rxt(k,17)
         mat(k,56) = rxt(k,18)
         mat(k,16) = 2.000*rxt(k,25) + .750*rxt(k,41) + rxt(k,50)

         mat(k,31) = rxt(k,8) + rxt(k,39)
         mat(k,61) = rxt(k,9)
         mat(k,140) = 2.000*rxt(k,10)
         mat(k,132) = rxt(k,13)
         mat(k,47) = rxt(k,16)
         mat(k,52) = rxt(k,17)
         mat(k,35) = rxt(k,20)
         mat(k,17) = .400*rxt(k,41) + rxt(k,50)
         mat(k,180) = rxt(k,59) + rxt(k,49)*y(k,13)

         mat(k,10) = -( rxt(k,12) )


         mat(k,129) = -( rxt(k,13) )
         mat(k,45) = rxt(k,16)

         mat(k,78) = .600*rxt(k,15) + rxt(k,67)
         mat(k,57) = rxt(k,18)
         mat(k,83) = rxt(k,19)
         mat(k,37) = rxt(k,20)

         mat(k,38) = -( rxt(k,14) )

         mat(k,72) = -( rxt(k,15) + rxt(k,67) )


         mat(k,172) = rxt(k,70)

         mat(k,43) = -( rxt(k,16) )

         mat(k,174) = rxt(k,81)

         mat(k,48) = -( rxt(k,17) )

         mat(k,79) = -( rxt(k,19) )
         mat(k,49) = .820*rxt(k,17)

         mat(k,53) = -( rxt(k,18) )


         mat(k,32) = -( rxt(k,20) )




         mat(k,5) = -( rxt(k,99) )


         mat(k,4) = rxt(k,99)

      end do

      end subroutine IMP_LINMAT01

      subroutine IMP_LINMAT( mat, y, rxt, het_rates )
!----------------------------------------------
!       ... Linear Matrix entries for Implicit species
!----------------------------------------------

      use MO_GRID_MOD,   only : pcnstm1
      use CHEM_MODS_MOD, only : rxntot, hetcnt, imp_nzcnt, clsze

      implicit none

!----------------------------------------------
!       ... Dummy args
!----------------------------------------------
      real, intent(in)    ::  y(clsze,pcnstm1)
      real, intent(in)    ::  rxt(clsze,rxntot)
      real, intent(in)    ::  het_rates(clsze,hetcnt)
      real, intent(inout) ::  mat(clsze,imp_nzcnt)

      call IMP_LINMAT01( mat, y, rxt, het_rates )

      end subroutine IMP_LINMAT

      end module MO_IMP_LIN_MATRIX_MOD

      module MO_ROD_LIN_MATRIX_MOD

      CONTAINS

      subroutine ROD_LINMAT( mat, y, rxt, het_rates )
!----------------------------------------------
!       ... Linear Matrix entries for Implicit species
!----------------------------------------------

      use MO_GRID_MOD,   only : pcnstm1
      use CHEM_MODS_MOD, only : rxntot, hetcnt, rod_nzcnt, clsze

      implicit none

!----------------------------------------------
!       ... Dummy args
!----------------------------------------------
      real, intent(in)    ::  y(clsze,pcnstm1)
      real, intent(in)    ::  rxt(clsze,rxntot)
      real, intent(in)    ::  het_rates(clsze,hetcnt)
      real, intent(inout) ::  mat(clsze,rod_nzcnt)


      end subroutine ROD_LINMAT

      end module MO_ROD_LIN_MATRIX_MOD

      module MO_IMP_NLN_MATRIX_MOD

      CONTAINS

      subroutine IMP_NLNMAT01( mat, y, rxt )

      use MO_GRID_MOD,   only : pcnstm1
      use CHEM_MODS_MOD, only : rxntot, imp_nzcnt, clsze

      implicit none

!----------------------------------------------
!       ... Dummy args
!----------------------------------------------
      real, intent(in)    ::  y(clsze,pcnstm1)
      real, intent(in)    ::  rxt(clsze,rxntot)
      real, intent(inout) ::  mat(clsze,imp_nzcnt)


!----------------------------------------------
!       ... Local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... Complete matrix entries Implicit species
!----------------------------------------------

      do k = 1,clsze
         mat(k,146) = -(rxt(k,22)*y(k,2) + rxt(k,27)*y(k,4) + rxt(k,29)*y(k,5) &
                      + rxt(k,53)*y(k,14) + rxt(k,54)*y(k,15) + rxt(k,69)*y(k,18))
         mat(k,85) = -rxt(k,22)*y(k,1)
         mat(k,243) = -rxt(k,27)*y(k,1)
         mat(k,191) = -rxt(k,29)*y(k,1)
         mat(k,177) = -rxt(k,53)*y(k,1)
         mat(k,214) = -rxt(k,54)*y(k,1)
         mat(k,92) = -rxt(k,69)*y(k,1)

         mat(k,214) = mat(k,214) + .300*rxt(k,64)*y(k,20)
         mat(k,270) = .300*rxt(k,64)*y(k,15)


         mat(k,84) = -(rxt(k,22)*y(k,1) + rxt(k,28)*y(k,5) + rxt(k,51)*y(k,14) &
                      + rxt(k,52)*y(k,15))
         mat(k,143) = -rxt(k,22)*y(k,2)
         mat(k,189) = -rxt(k,28)*y(k,2)
         mat(k,170) = -rxt(k,51)*y(k,2)
         mat(k,208) = -rxt(k,52)*y(k,2)

         mat(k,170) = mat(k,170) + 2.000*rxt(k,58)*y(k,14)

         mat(k,248) = -(rxt(k,26)*y(k,15) + rxt(k,27)*y(k,1) + rxt(k,36)*y(k,6) &
                      + rxt(k,42)*y(k,10) + rxt(k,62)*y(k,20) + rxt(k,71)*y(k,24) &
                      + rxt(k,77)*y(k,23) + rxt(k,82)*y(k,26) + rxt(k,87)*y(k,30))
         mat(k,219) = -rxt(k,26)*y(k,4)
         mat(k,151) = -rxt(k,27)*y(k,4)
         mat(k,264) = -rxt(k,36)*y(k,4)
         mat(k,233) = -rxt(k,42)*y(k,4)
         mat(k,275) = -rxt(k,62)*y(k,4)
         mat(k,107) = -rxt(k,71)*y(k,4)
         mat(k,70) = -rxt(k,77)*y(k,4)
         mat(k,126) = -rxt(k,82)*y(k,4)
         mat(k,114) = -rxt(k,87)*y(k,4)

         mat(k,89) = rxt(k,28)*y(k,5)
         mat(k,196) = rxt(k,28)*y(k,2)

         mat(k,193) = -(rxt(k,28)*y(k,2) + rxt(k,29)*y(k,1) + rxt(k,31)*y(k,6) &
                      + rxt(k,34)*y(k,14) + rxt(k,37)*y(k,15) + rxt(k,63)*y(k,20))
         mat(k,87) = -rxt(k,28)*y(k,5)
         mat(k,148) = -rxt(k,29)*y(k,5)
         mat(k,261) = -rxt(k,31)*y(k,5)
         mat(k,179) = -rxt(k,34)*y(k,5)
         mat(k,216) = -rxt(k,37)*y(k,5)
         mat(k,272) = -rxt(k,63)*y(k,5)

         mat(k,148) = mat(k,148) + rxt(k,27)*y(k,4)
         mat(k,245) = rxt(k,27)*y(k,1) + 2.000*rxt(k,36)*y(k,6) + rxt(k,42)*y(k,10)  &
                      + rxt(k,26)*y(k,15) + rxt(k,62)*y(k,20) + rxt(k,77)*y(k,23)  &
                      + rxt(k,71)*y(k,24) + rxt(k,82)*y(k,26) + rxt(k,87)*y(k,30)
         mat(k,261) = mat(k,261) + 2.000*rxt(k,36)*y(k,4) + rxt(k,30)*y(k,15)  &
                      + rxt(k,90)*y(k,18) + rxt(k,98)*y(k,38)
         mat(k,30) = rxt(k,38)*y(k,14)
         mat(k,230) = rxt(k,42)*y(k,4)
         mat(k,179) = mat(k,179) + rxt(k,38)*y(k,8)
         mat(k,216) = mat(k,216) + rxt(k,26)*y(k,4) + rxt(k,30)*y(k,6)
         mat(k,94) = rxt(k,90)*y(k,6)
         mat(k,272) = mat(k,272) + rxt(k,62)*y(k,4)
         mat(k,68) = rxt(k,77)*y(k,4)
         mat(k,104) = rxt(k,71)*y(k,4)
         mat(k,123) = rxt(k,82)*y(k,4)
         mat(k,112) = rxt(k,87)*y(k,4)
         mat(k,26) = rxt(k,98)*y(k,6)

         mat(k,265) = -(rxt(k,30)*y(k,15) + rxt(k,31)*y(k,5) + rxt(k,36)*y(k,4) &
                      + rxt(k,47)*y(k,12) + rxt(k,61)*y(k,19) + rxt(k,90)*y(k,18) &
                      + rxt(k,92)*y(k,31) + rxt(k,98)*y(k,38))
         mat(k,220) = -rxt(k,30)*y(k,6)
         mat(k,197) = -rxt(k,31)*y(k,6)
         mat(k,249) = -rxt(k,36)*y(k,6)
         mat(k,141) = -rxt(k,47)*y(k,6)
         mat(k,134) = -rxt(k,61)*y(k,6)
         mat(k,97) = -rxt(k,90)*y(k,6)
         mat(k,36) = -rxt(k,92)*y(k,6)
         mat(k,27) = -rxt(k,98)*y(k,6)

         mat(k,152) = rxt(k,29)*y(k,5)
         mat(k,197) = mat(k,197) + rxt(k,29)*y(k,1)
         mat(k,66) = rxt(k,35)*y(k,14)
         mat(k,183) = rxt(k,35)*y(k,7) + rxt(k,93)*y(k,22)
         mat(k,77) = rxt(k,93)*y(k,14)

         mat(k,63) = -(rxt(k,35)*y(k,14))
         mat(k,166) = -rxt(k,35)*y(k,7)

         mat(k,187) = rxt(k,34)*y(k,14)
         mat(k,255) = rxt(k,47)*y(k,12) + rxt(k,61)*y(k,19) + rxt(k,92)*y(k,31)
         mat(k,136) = rxt(k,47)*y(k,6)
         mat(k,166) = mat(k,166) + rxt(k,34)*y(k,5)
         mat(k,128) = rxt(k,61)*y(k,6)
         mat(k,33) = rxt(k,92)*y(k,6)

         mat(k,28) = -(rxt(k,38)*y(k,14))
         mat(k,159) = -rxt(k,38)*y(k,8)

         mat(k,186) = rxt(k,37)*y(k,15)
         mat(k,200) = rxt(k,37)*y(k,5)


         mat(k,185) = rxt(k,31)*y(k,6)
         mat(k,252) = rxt(k,31)*y(k,5)

         mat(k,232) = -(rxt(k,42)*y(k,4) + (4.*rxt(k,43) + 4.*rxt(k,44)) * y(k,10) &
                      + rxt(k,45)*y(k,15) + rxt(k,65)*y(k,20) + rxt(k,73)*y(k,24) &
                      + rxt(k,84)*y(k,26))
         mat(k,247) = -rxt(k,42)*y(k,10)
         mat(k,218) = -rxt(k,45)*y(k,10)
         mat(k,274) = -rxt(k,65)*y(k,10)
         mat(k,106) = -rxt(k,73)*y(k,10)
         mat(k,125) = -rxt(k,84)*y(k,10)

         mat(k,247) = mat(k,247) + rxt(k,62)*y(k,20)
         mat(k,232) = mat(k,232) + .900*rxt(k,65)*y(k,20)
         mat(k,62) = .700*rxt(k,46)*y(k,14)
         mat(k,181) = .700*rxt(k,46)*y(k,11)
         mat(k,274) = mat(k,274) + rxt(k,62)*y(k,4) + .900*rxt(k,65)*y(k,10)  &
                      + 4.000*rxt(k,68)*y(k,20)

         mat(k,58) = -(rxt(k,46)*y(k,14))
         mat(k,165) = -rxt(k,46)*y(k,11)

         mat(k,222) = rxt(k,45)*y(k,15)
         mat(k,205) = rxt(k,45)*y(k,10)

         mat(k,137) = -(rxt(k,47)*y(k,6) + rxt(k,48)*y(k,14))
         mat(k,258) = -rxt(k,47)*y(k,12)
         mat(k,176) = -rxt(k,48)*y(k,12)

         mat(k,242) = rxt(k,42)*y(k,10) + rxt(k,87)*y(k,30)
         mat(k,228) = rxt(k,42)*y(k,4) + (4.000*rxt(k,43)+2.000*rxt(k,44))*y(k,10)  &
                      + rxt(k,65)*y(k,20) + .700*rxt(k,73)*y(k,24) + rxt(k,84)*y(k,26)
         mat(k,59) = .300*rxt(k,46)*y(k,14)
         mat(k,176) = mat(k,176) + .300*rxt(k,46)*y(k,11) + .500*rxt(k,66)*y(k,21)  &
                      + rxt(k,93)*y(k,22)
         mat(k,269) = rxt(k,65)*y(k,10)
         mat(k,39) = .500*rxt(k,66)*y(k,14)
         mat(k,73) = rxt(k,93)*y(k,14)
         mat(k,102) = .700*rxt(k,73)*y(k,10)
         mat(k,121) = rxt(k,84)*y(k,10)
         mat(k,110) = rxt(k,87)*y(k,4)

         mat(k,178) = -(rxt(k,34)*y(k,5) + rxt(k,35)*y(k,7) + rxt(k,38)*y(k,8) &
                      + rxt(k,46)*y(k,11) + rxt(k,48)*y(k,12) + rxt(k,51)*y(k,2) &
                      + rxt(k,53)*y(k,1) + rxt(k,56)*y(k,16) + rxt(k,57)*y(k,15) &
                      + 4.*rxt(k,58)*y(k,14) + rxt(k,60)*y(k,19) + rxt(k,66)*y(k,21) &
                      + rxt(k,75)*y(k,25) + rxt(k,76)*y(k,18) + rxt(k,85)*y(k,27) &
                      + rxt(k,86)*y(k,28) + rxt(k,89)*y(k,29) + rxt(k,91)*y(k,31) &
                      + rxt(k,93)*y(k,22) + rxt(k,96)*y(k,36) + rxt(k,97)*y(k,38) &
                      + rxt(k,100)*y(k,39))
         mat(k,192) = -rxt(k,34)*y(k,14)
         mat(k,64) = -rxt(k,35)*y(k,14)
         mat(k,29) = -rxt(k,38)*y(k,14)
         mat(k,60) = -rxt(k,46)*y(k,14)
         mat(k,138) = -rxt(k,48)*y(k,14)
         mat(k,86) = -rxt(k,51)*y(k,14)
         mat(k,147) = -rxt(k,53)*y(k,14)
         mat(k,11) = -rxt(k,56)*y(k,14)
         mat(k,215) = -rxt(k,57)*y(k,14)
         mat(k,130) = -rxt(k,60)*y(k,14)
         mat(k,40) = -rxt(k,66)*y(k,14)
         mat(k,46) = -rxt(k,75)*y(k,14)
         mat(k,93) = -rxt(k,76)*y(k,14)
         mat(k,51) = -rxt(k,85)*y(k,14)
         mat(k,81) = -rxt(k,86)*y(k,14)
         mat(k,56) = -rxt(k,89)*y(k,14)
         mat(k,34) = -rxt(k,91)*y(k,14)
         mat(k,74) = -rxt(k,93)*y(k,14)
         mat(k,9) = -rxt(k,96)*y(k,14)
         mat(k,25) = -rxt(k,97)*y(k,14)
         mat(k,6) = -rxt(k,100)*y(k,14)

         mat(k,147) = mat(k,147) + rxt(k,54)*y(k,15)
         mat(k,86) = mat(k,86) + rxt(k,52)*y(k,15)
         mat(k,244) = rxt(k,26)*y(k,15)
         mat(k,260) = rxt(k,30)*y(k,15)
         mat(k,60) = mat(k,60) + .300*rxt(k,46)*y(k,14)
         mat(k,178) = mat(k,178) + .300*rxt(k,46)*y(k,11) + .500*rxt(k,75)*y(k,25)
         mat(k,215) = mat(k,215) + rxt(k,54)*y(k,1) + rxt(k,52)*y(k,2) + rxt(k,26) &
                      *y(k,4) + rxt(k,30)*y(k,6)
         mat(k,46) = mat(k,46) + .500*rxt(k,75)*y(k,14)

         mat(k,217) = -(rxt(k,26)*y(k,4) + rxt(k,30)*y(k,6) + rxt(k,37)*y(k,5) &
                      + rxt(k,45)*y(k,10) + rxt(k,52)*y(k,2) + rxt(k,54)*y(k,1) &
                      + 4.*rxt(k,55)*y(k,15) + rxt(k,57)*y(k,14) + rxt(k,64)*y(k,20) &
                      + rxt(k,72)*y(k,24) + rxt(k,78)*y(k,23) + rxt(k,83)*y(k,26) &
                      + rxt(k,88)*y(k,30))
         mat(k,246) = -rxt(k,26)*y(k,15)
         mat(k,262) = -rxt(k,30)*y(k,15)
         mat(k,194) = -rxt(k,37)*y(k,15)
         mat(k,231) = -rxt(k,45)*y(k,15)
         mat(k,88) = -rxt(k,52)*y(k,15)
         mat(k,149) = -rxt(k,54)*y(k,15)
         mat(k,180) = -rxt(k,57)*y(k,15)
         mat(k,273) = -rxt(k,64)*y(k,15)
         mat(k,105) = -rxt(k,72)*y(k,15)
         mat(k,69) = -rxt(k,78)*y(k,15)
         mat(k,124) = -rxt(k,83)*y(k,15)
         mat(k,113) = -rxt(k,88)*y(k,15)

         mat(k,149) = mat(k,149) + rxt(k,53)*y(k,14) + rxt(k,69)*y(k,18)
         mat(k,88) = mat(k,88) + rxt(k,51)*y(k,14)
         mat(k,246) = mat(k,246) + rxt(k,42)*y(k,10) + rxt(k,77)*y(k,23) + rxt(k,71) &
                      *y(k,24) + rxt(k,82)*y(k,26)
         mat(k,262) = mat(k,262) + rxt(k,47)*y(k,12)
         mat(k,231) = mat(k,231) + rxt(k,42)*y(k,4) + 4.000*rxt(k,43)*y(k,10)  &
                      + .900*rxt(k,65)*y(k,20) + rxt(k,73)*y(k,24) + rxt(k,84)*y(k,26)
         mat(k,140) = rxt(k,47)*y(k,6) + rxt(k,48)*y(k,14)
         mat(k,180) = mat(k,180) + rxt(k,53)*y(k,1) + rxt(k,51)*y(k,2) + rxt(k,48) &
                      *y(k,12) + rxt(k,56)*y(k,16)
         mat(k,12) = rxt(k,56)*y(k,14)
         mat(k,95) = rxt(k,69)*y(k,1)
         mat(k,273) = mat(k,273) + .900*rxt(k,65)*y(k,10)
         mat(k,69) = mat(k,69) + rxt(k,77)*y(k,4)
         mat(k,105) = mat(k,105) + rxt(k,71)*y(k,4) + rxt(k,73)*y(k,10)  &
                      + 2.400*rxt(k,74)*y(k,24)
         mat(k,124) = mat(k,124) + rxt(k,82)*y(k,4) + rxt(k,84)*y(k,10)

         mat(k,10) = -(rxt(k,56)*y(k,14))
         mat(k,157) = -rxt(k,56)*y(k,16)

         mat(k,199) = 2.000*rxt(k,55)*y(k,15)

         mat(k,91) = -(rxt(k,69)*y(k,1) + rxt(k,76)*y(k,14) + rxt(k,90)*y(k,6))
         mat(k,144) = -rxt(k,69)*y(k,18)
         mat(k,171) = -rxt(k,76)*y(k,18)
         mat(k,256) = -rxt(k,90)*y(k,18)

         mat(k,129) = -(rxt(k,60)*y(k,14) + rxt(k,61)*y(k,6))
         mat(k,175) = -rxt(k,60)*y(k,19)
         mat(k,257) = -rxt(k,61)*y(k,19)

         mat(k,241) = rxt(k,71)*y(k,24) + .270*rxt(k,82)*y(k,26)
         mat(k,227) = .800*rxt(k,73)*y(k,24)
         mat(k,175) = mat(k,175) + .500*rxt(k,75)*y(k,25)
         mat(k,101) = rxt(k,71)*y(k,4) + .800*rxt(k,73)*y(k,10) + 3.200*rxt(k,74) &
                      *y(k,24)
         mat(k,45) = .500*rxt(k,75)*y(k,14)
         mat(k,120) = .270*rxt(k,82)*y(k,4)

      end do

      end subroutine IMP_NLNMAT01

      subroutine IMP_NLNMAT02( mat, y, rxt )

      use MO_GRID_MOD,   only : pcnstm1
      use CHEM_MODS_MOD, only : rxntot, imp_nzcnt, clsze

      implicit none

!----------------------------------------------
!       ... Dummy args
!----------------------------------------------
      real, intent(in)    ::  y(clsze,pcnstm1)
      real, intent(in)    ::  rxt(clsze,rxntot)
      real, intent(inout) ::  mat(clsze,imp_nzcnt)


!----------------------------------------------
!       ... Local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... Complete matrix entries Implicit species
!----------------------------------------------

      do k = 1,clsze
         mat(k,277) = -(rxt(k,62)*y(k,4) + rxt(k,63)*y(k,5) + rxt(k,64)*y(k,15) &
                      + rxt(k,65)*y(k,10) + 4.*rxt(k,68)*y(k,20))
         mat(k,250) = -rxt(k,62)*y(k,20)
         mat(k,198) = -rxt(k,63)*y(k,20)
         mat(k,221) = -rxt(k,64)*y(k,20)
         mat(k,235) = -rxt(k,65)*y(k,20)

         mat(k,153) = rxt(k,69)*y(k,18)
         mat(k,250) = mat(k,250) + rxt(k,77)*y(k,23) + rxt(k,87)*y(k,30)
         mat(k,266) = rxt(k,90)*y(k,18) + rxt(k,61)*y(k,19) + rxt(k,92)*y(k,31)
         mat(k,184) = rxt(k,60)*y(k,19) + .500*rxt(k,66)*y(k,21) + rxt(k,91)*y(k,31)
         mat(k,98) = rxt(k,69)*y(k,1) + rxt(k,90)*y(k,6)
         mat(k,135) = rxt(k,61)*y(k,6) + rxt(k,60)*y(k,14)
         mat(k,42) = .500*rxt(k,66)*y(k,14)
         mat(k,71) = rxt(k,77)*y(k,4)
         mat(k,115) = rxt(k,87)*y(k,4)
         mat(k,37) = rxt(k,92)*y(k,6) + rxt(k,91)*y(k,14)

         mat(k,38) = -(rxt(k,66)*y(k,14))
         mat(k,161) = -rxt(k,66)*y(k,21)

         mat(k,201) = .700*rxt(k,64)*y(k,20)
         mat(k,267) = .700*rxt(k,64)*y(k,15)

         mat(k,72) = -(rxt(k,93)*y(k,14))
         mat(k,168) = -rxt(k,93)*y(k,22)

         mat(k,188) = rxt(k,63)*y(k,20)
         mat(k,268) = rxt(k,63)*y(k,5)

         mat(k,67) = -(rxt(k,77)*y(k,4) + rxt(k,78)*y(k,15))
         mat(k,236) = -rxt(k,77)*y(k,23)
         mat(k,206) = -rxt(k,78)*y(k,23)

         mat(k,167) = rxt(k,76)*y(k,18)
         mat(k,90) = rxt(k,76)*y(k,14)

         mat(k,100) = -(rxt(k,71)*y(k,4) + rxt(k,72)*y(k,15) + rxt(k,73)*y(k,10) &
                      + 4.*rxt(k,74)*y(k,24))
         mat(k,238) = -rxt(k,71)*y(k,24)
         mat(k,209) = -rxt(k,72)*y(k,24)
         mat(k,224) = -rxt(k,73)*y(k,24)

         mat(k,172) = .500*rxt(k,75)*y(k,25)
         mat(k,44) = .500*rxt(k,75)*y(k,14)

         mat(k,43) = -(rxt(k,75)*y(k,14))
         mat(k,162) = -rxt(k,75)*y(k,25)

         mat(k,202) = rxt(k,72)*y(k,24)
         mat(k,99) = rxt(k,72)*y(k,15)

         mat(k,119) = -(rxt(k,82)*y(k,4) + rxt(k,83)*y(k,15) + rxt(k,84)*y(k,10))
         mat(k,240) = -rxt(k,82)*y(k,26)
         mat(k,211) = -rxt(k,83)*y(k,26)
         mat(k,226) = -rxt(k,84)*y(k,26)

         mat(k,174) = rxt(k,85)*y(k,27)
         mat(k,50) = rxt(k,85)*y(k,14)

         mat(k,48) = -(rxt(k,85)*y(k,14))
         mat(k,163) = -rxt(k,85)*y(k,27)

         mat(k,203) = rxt(k,83)*y(k,26)
         mat(k,116) = rxt(k,83)*y(k,15)

         mat(k,79) = -(rxt(k,86)*y(k,14))
         mat(k,169) = -rxt(k,86)*y(k,28)

         mat(k,237) = .820*rxt(k,82)*y(k,26)
         mat(k,223) = .820*rxt(k,84)*y(k,26)
         mat(k,117) = .820*rxt(k,82)*y(k,4) + .820*rxt(k,84)*y(k,10)

         mat(k,53) = -(rxt(k,89)*y(k,14))
         mat(k,164) = -rxt(k,89)*y(k,29)

         mat(k,204) = rxt(k,88)*y(k,30)
         mat(k,108) = rxt(k,88)*y(k,15)

         mat(k,109) = -(rxt(k,87)*y(k,4) + rxt(k,88)*y(k,15))
         mat(k,239) = -rxt(k,87)*y(k,30)
         mat(k,210) = -rxt(k,88)*y(k,30)

         mat(k,173) = rxt(k,86)*y(k,28) + rxt(k,89)*y(k,29)
         mat(k,80) = rxt(k,86)*y(k,14)
         mat(k,54) = rxt(k,89)*y(k,14)

         mat(k,32) = -(rxt(k,91)*y(k,14) + rxt(k,92)*y(k,6))
         mat(k,160) = -rxt(k,91)*y(k,31)
         mat(k,254) = -rxt(k,92)*y(k,31)

         mat(k,8) = -(rxt(k,96)*y(k,14))
         mat(k,156) = -rxt(k,96)*y(k,36)

         mat(k,251) = rxt(k,98)*y(k,38)
         mat(k,156) = mat(k,156) + rxt(k,97)*y(k,38)
         mat(k,23) = rxt(k,98)*y(k,6) + rxt(k,97)*y(k,14)


         mat(k,154) = rxt(k,96)*y(k,36)
         mat(k,7) = rxt(k,96)*y(k,14)

         mat(k,24) = -(rxt(k,97)*y(k,14) + rxt(k,98)*y(k,6))
         mat(k,158) = -rxt(k,97)*y(k,38)
         mat(k,253) = -rxt(k,98)*y(k,38)

         mat(k,5) = -(rxt(k,100)*y(k,14))
         mat(k,155) = -rxt(k,100)*y(k,39)



      end do

      end subroutine IMP_NLNMAT02

      subroutine IMP_NLNMAT_FINIT( mat, lmat, dti )

      use MO_GRID_MOD,   only : pcnstm1
      use CHEM_MODS_MOD, only : rxntot, imp_nzcnt, clsze

      implicit none

!----------------------------------------------
!       ... Dummy args
!----------------------------------------------
      real, intent(in)    ::  dti
      real, intent(in)    ::  lmat(clsze,imp_nzcnt)
      real, intent(inout) ::  mat(clsze,imp_nzcnt)


!----------------------------------------------
!       ... Local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... Complete matrix entries Implicit species
!----------------------------------------------

      do k = 1,clsze
         mat(k,  4) = lmat(k,  4)
         mat(k,  5) = mat(k,  5) + lmat(k,  5)
         mat(k, 10) = mat(k, 10) + lmat(k, 10)
         mat(k, 11) = mat(k, 11) + lmat(k, 11)
         mat(k, 13) = lmat(k, 13)
         mat(k, 14) = lmat(k, 14)
         mat(k, 15) = lmat(k, 15)
         mat(k, 16) = lmat(k, 16)
         mat(k, 17) = lmat(k, 17)
         mat(k, 18) = lmat(k, 18)
         mat(k, 19) = lmat(k, 19)
         mat(k, 20) = lmat(k, 20)
         mat(k, 21) = lmat(k, 21)
         mat(k, 22) = lmat(k, 22)
         mat(k, 28) = mat(k, 28) + lmat(k, 28)
         mat(k, 30) = mat(k, 30) + lmat(k, 30)
         mat(k, 31) = lmat(k, 31)
         mat(k, 32) = mat(k, 32) + lmat(k, 32)
         mat(k, 35) = lmat(k, 35)
         mat(k, 37) = mat(k, 37) + lmat(k, 37)
         mat(k, 38) = mat(k, 38) + lmat(k, 38)
         mat(k, 40) = mat(k, 40) + lmat(k, 40)
         mat(k, 41) = lmat(k, 41)
         mat(k, 43) = mat(k, 43) + lmat(k, 43)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 46) = mat(k, 46) + lmat(k, 46)
         mat(k, 47) = lmat(k, 47)
         mat(k, 48) = mat(k, 48) + lmat(k, 48)
         mat(k, 49) = lmat(k, 49)
         mat(k, 51) = mat(k, 51) + lmat(k, 51)
         mat(k, 52) = lmat(k, 52)
         mat(k, 53) = mat(k, 53) + lmat(k, 53)
         mat(k, 55) = lmat(k, 55)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 57) = lmat(k, 57)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 60) = mat(k, 60) + lmat(k, 60)
         mat(k, 61) = lmat(k, 61)
         mat(k, 63) = mat(k, 63) + lmat(k, 63)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 72) = mat(k, 72) + lmat(k, 72)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = lmat(k, 76)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 82) = lmat(k, 82)
         mat(k, 83) = lmat(k, 83)
         mat(k, 84) = mat(k, 84) + lmat(k, 84)
         mat(k, 85) = mat(k, 85) + lmat(k, 85)
         mat(k,129) = mat(k,129) + lmat(k,129)
         mat(k,132) = lmat(k,132)
         mat(k,133) = lmat(k,133)
         mat(k,137) = mat(k,137) + lmat(k,137)
         mat(k,140) = mat(k,140) + lmat(k,140)
         mat(k,142) = lmat(k,142)
         mat(k,143) = mat(k,143) + lmat(k,143)
         mat(k,146) = mat(k,146) + lmat(k,146)
         mat(k,172) = mat(k,172) + lmat(k,172)
         mat(k,174) = mat(k,174) + lmat(k,174)
         mat(k,178) = mat(k,178) + lmat(k,178)
         mat(k,180) = mat(k,180) + lmat(k,180)
         mat(k,181) = mat(k,181) + lmat(k,181)
         mat(k,189) = mat(k,189) + lmat(k,189)
         mat(k,193) = mat(k,193) + lmat(k,193)
         mat(k,196) = mat(k,196) + lmat(k,196)
         mat(k,255) = mat(k,255) + lmat(k,255)
         mat(k,259) = lmat(k,259)
         mat(k,261) = mat(k,261) + lmat(k,261)
         mat(k,264) = mat(k,264) + lmat(k,264)
         mat(k,265) = mat(k,265) + lmat(k,265)
         mat(k,  1) = 0.
         mat(k,  2) = 0.
         mat(k,  3) = 0.
         mat(k, 96) = 0.
         mat(k,103) = 0.
         mat(k,111) = 0.
         mat(k,118) = 0.
         mat(k,122) = 0.
         mat(k,127) = 0.
         mat(k,131) = 0.
         mat(k,139) = 0.
         mat(k,145) = 0.
         mat(k,150) = 0.
         mat(k,182) = 0.
         mat(k,190) = 0.
         mat(k,195) = 0.
         mat(k,207) = 0.
         mat(k,212) = 0.
         mat(k,213) = 0.
         mat(k,225) = 0.
         mat(k,229) = 0.
         mat(k,234) = 0.
         mat(k,263) = 0.
         mat(k,271) = 0.
         mat(k,276) = 0.
         mat(k,  1) = -dti
         mat(k,  2) = -dti
         mat(k,  3) = -dti
         mat(k,  5) = mat(k,  5) - dti
         mat(k,  8) = mat(k,  8) - dti
         mat(k, 10) = mat(k, 10) - dti
         mat(k, 13) = mat(k, 13) - dti
         mat(k, 19) = mat(k, 19) - dti
         mat(k, 24) = mat(k, 24) - dti
         mat(k, 28) = mat(k, 28) - dti
         mat(k, 32) = mat(k, 32) - dti
         mat(k, 38) = mat(k, 38) - dti
         mat(k, 43) = mat(k, 43) - dti
         mat(k, 48) = mat(k, 48) - dti
         mat(k, 53) = mat(k, 53) - dti
         mat(k, 58) = mat(k, 58) - dti
         mat(k, 63) = mat(k, 63) - dti
         mat(k, 67) = mat(k, 67) - dti
         mat(k, 72) = mat(k, 72) - dti
         mat(k, 79) = mat(k, 79) - dti
         mat(k, 84) = mat(k, 84) - dti
         mat(k, 91) = mat(k, 91) - dti
         mat(k,100) = mat(k,100) - dti
         mat(k,109) = mat(k,109) - dti
         mat(k,119) = mat(k,119) - dti
         mat(k,129) = mat(k,129) - dti
         mat(k,137) = mat(k,137) - dti
         mat(k,146) = mat(k,146) - dti
         mat(k,178) = mat(k,178) - dti
         mat(k,193) = mat(k,193) - dti
         mat(k,217) = mat(k,217) - dti
         mat(k,232) = mat(k,232) - dti
         mat(k,248) = mat(k,248) - dti
         mat(k,265) = mat(k,265) - dti
         mat(k,277) = mat(k,277) - dti
      end do

      end subroutine IMP_NLNMAT_FINIT

      subroutine IMP_NLNMAT( mat, y, rxt, lmat, dti )

      use MO_GRID_MOD,   only : pcnstm1
      use CHEM_MODS_MOD, only : rxntot, imp_nzcnt, clsze

      implicit none

!----------------------------------------------
!       ... Dummy args
!----------------------------------------------
      real, intent(in)    ::  dti
      real, intent(in)    ::  lmat(clsze,imp_nzcnt)
      real, intent(in)    ::  y(clsze,pcnstm1)
      real, intent(in)    ::  rxt(clsze,rxntot)
      real, intent(inout) ::  mat(clsze,imp_nzcnt)

      call IMP_NLNMAT01( mat, y, rxt )
      call IMP_NLNMAT02( mat, y, rxt )
      call IMP_NLNMAT_FINIT( mat, lmat, dti )

      end subroutine IMP_NLNMAT

      end module MO_IMP_NLN_MATRIX_MOD

      module MO_ROD_NLN_MATRIX_MOD

      CONTAINS

      subroutine ROD_NLNMAT( mat, y, rxt, lmat, dti )

      use MO_GRID_MOD,   only : pcnstm1
      use CHEM_MODS_MOD, only : rxntot, rod_nzcnt, clsze

      implicit none

!----------------------------------------------
!       ... Dummy args
!----------------------------------------------
      real, intent(in)    ::  dti
      real, intent(in)    ::  lmat(clsze,rod_nzcnt)
      real, intent(in)    ::  y(clsze,pcnstm1)
      real, intent(in)    ::  rxt(clsze,rxntot)
      real, intent(inout) ::  mat(clsze,rod_nzcnt)


      end subroutine ROD_NLNMAT

      end module MO_ROD_NLN_MATRIX_MOD

      module MO_IMP_FACTOR_MOD

      CONTAINS
                                                                        
      subroutine IMP_LU_FAC01( lu )
                                                                        
      use CHEM_MODS_MOD, only : imp_nzcnt, clsze
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real, intent(inout) ::   lu(clsze,imp_nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,clsze
         lu(k,1) = 1. / lu(k,1)
                                                                        
         lu(k,2) = 1. / lu(k,2)
                                                                        
         lu(k,3) = 1. / lu(k,3)
                                                                        
         lu(k,5) = 1. / lu(k,5)
         lu(k,6) = lu(k,6) * lu(k,5)
         lu(k,178) = lu(k,178) - lu(k,6) * lu(k,155)
                                                                        
         lu(k,8) = 1. / lu(k,8)
         lu(k,9) = lu(k,9) * lu(k,8)
         lu(k,25) = lu(k,25) - lu(k,9) * lu(k,23)
         lu(k,178) = lu(k,178) - lu(k,9) * lu(k,156)
         lu(k,260) = lu(k,260) - lu(k,9) * lu(k,251)
                                                                        
         lu(k,10) = 1. / lu(k,10)
         lu(k,11) = lu(k,11) * lu(k,10)
         lu(k,12) = lu(k,12) * lu(k,10)
         lu(k,178) = lu(k,178) - lu(k,11) * lu(k,157)
         lu(k,180) = lu(k,180) - lu(k,12) * lu(k,157)
         lu(k,215) = lu(k,215) - lu(k,11) * lu(k,199)
         lu(k,217) = lu(k,217) - lu(k,12) * lu(k,199)
                                                                        
         lu(k,13) = 1. / lu(k,13)
         lu(k,14) = lu(k,14) * lu(k,13)
         lu(k,15) = lu(k,15) * lu(k,13)
         lu(k,16) = lu(k,16) * lu(k,13)
         lu(k,17) = lu(k,17) * lu(k,13)
         lu(k,18) = lu(k,18) * lu(k,13)
         lu(k,143) = lu(k,143) - lu(k,14) * lu(k,142)
         lu(k,145) = - lu(k,15) * lu(k,142)
         lu(k,147) = lu(k,147) - lu(k,16) * lu(k,142)
         lu(k,149) = lu(k,149) - lu(k,17) * lu(k,142)
         lu(k,150) = - lu(k,18) * lu(k,142)
                                                                        
         lu(k,19) = 1. / lu(k,19)
         lu(k,20) = lu(k,20) * lu(k,19)
         lu(k,21) = lu(k,21) * lu(k,19)
         lu(k,22) = lu(k,22) * lu(k,19)
         lu(k,187) = lu(k,187) - lu(k,20) * lu(k,185)
         lu(k,193) = lu(k,193) - lu(k,21) * lu(k,185)
         lu(k,197) = lu(k,197) - lu(k,22) * lu(k,185)
         lu(k,255) = lu(k,255) - lu(k,20) * lu(k,252)
         lu(k,261) = lu(k,261) - lu(k,21) * lu(k,252)
         lu(k,265) = lu(k,265) - lu(k,22) * lu(k,252)
                                                                        
         lu(k,24) = 1. / lu(k,24)
         lu(k,25) = lu(k,25) * lu(k,24)
         lu(k,26) = lu(k,26) * lu(k,24)
         lu(k,27) = lu(k,27) * lu(k,24)
         lu(k,178) = lu(k,178) - lu(k,25) * lu(k,158)
         lu(k,179) = lu(k,179) - lu(k,26) * lu(k,158)
         lu(k,183) = lu(k,183) - lu(k,27) * lu(k,158)
         lu(k,260) = lu(k,260) - lu(k,25) * lu(k,253)
         lu(k,261) = lu(k,261) - lu(k,26) * lu(k,253)
         lu(k,265) = lu(k,265) - lu(k,27) * lu(k,253)
                                                                        
         lu(k,28) = 1. / lu(k,28)
         lu(k,29) = lu(k,29) * lu(k,28)
         lu(k,30) = lu(k,30) * lu(k,28)
         lu(k,31) = lu(k,31) * lu(k,28)
         lu(k,178) = lu(k,178) - lu(k,29) * lu(k,159)
         lu(k,179) = lu(k,179) - lu(k,30) * lu(k,159)
         lu(k,180) = lu(k,180) - lu(k,31) * lu(k,159)
         lu(k,192) = lu(k,192) - lu(k,29) * lu(k,186)
         lu(k,193) = lu(k,193) - lu(k,30) * lu(k,186)
         lu(k,194) = lu(k,194) - lu(k,31) * lu(k,186)
         lu(k,215) = lu(k,215) - lu(k,29) * lu(k,200)
         lu(k,216) = lu(k,216) - lu(k,30) * lu(k,200)
         lu(k,217) = lu(k,217) - lu(k,31) * lu(k,200)
                                                                        
         lu(k,32) = 1. / lu(k,32)
         lu(k,33) = lu(k,33) * lu(k,32)
         lu(k,34) = lu(k,34) * lu(k,32)
         lu(k,35) = lu(k,35) * lu(k,32)
         lu(k,36) = lu(k,36) * lu(k,32)
         lu(k,37) = lu(k,37) * lu(k,32)
         lu(k,166) = lu(k,166) - lu(k,33) * lu(k,160)
         lu(k,178) = lu(k,178) - lu(k,34) * lu(k,160)
         lu(k,180) = lu(k,180) - lu(k,35) * lu(k,160)
         lu(k,183) = lu(k,183) - lu(k,36) * lu(k,160)
         lu(k,184) = lu(k,184) - lu(k,37) * lu(k,160)
         lu(k,255) = lu(k,255) - lu(k,33) * lu(k,254)
         lu(k,260) = lu(k,260) - lu(k,34) * lu(k,254)
         lu(k,262) = lu(k,262) - lu(k,35) * lu(k,254)
         lu(k,265) = lu(k,265) - lu(k,36) * lu(k,254)
         lu(k,266) = lu(k,266) - lu(k,37) * lu(k,254)
                                                                        
         lu(k,38) = 1. / lu(k,38)
         lu(k,39) = lu(k,39) * lu(k,38)
         lu(k,40) = lu(k,40) * lu(k,38)
         lu(k,41) = lu(k,41) * lu(k,38)
         lu(k,42) = lu(k,42) * lu(k,38)
         lu(k,176) = lu(k,176) - lu(k,39) * lu(k,161)
         lu(k,178) = lu(k,178) - lu(k,40) * lu(k,161)
         lu(k,181) = lu(k,181) - lu(k,41) * lu(k,161)
         lu(k,184) = lu(k,184) - lu(k,42) * lu(k,161)
         lu(k,213) = - lu(k,39) * lu(k,201)
         lu(k,215) = lu(k,215) - lu(k,40) * lu(k,201)
         lu(k,218) = lu(k,218) - lu(k,41) * lu(k,201)
         lu(k,221) = lu(k,221) - lu(k,42) * lu(k,201)
         lu(k,269) = lu(k,269) - lu(k,39) * lu(k,267)
         lu(k,271) = - lu(k,40) * lu(k,267)
         lu(k,274) = lu(k,274) - lu(k,41) * lu(k,267)
         lu(k,277) = lu(k,277) - lu(k,42) * lu(k,267)
                                                                        
         lu(k,43) = 1. / lu(k,43)
         lu(k,44) = lu(k,44) * lu(k,43)
         lu(k,45) = lu(k,45) * lu(k,43)
         lu(k,46) = lu(k,46) * lu(k,43)
         lu(k,47) = lu(k,47) * lu(k,43)
         lu(k,100) = lu(k,100) - lu(k,44) * lu(k,99)
         lu(k,101) = lu(k,101) - lu(k,45) * lu(k,99)
         lu(k,103) = - lu(k,46) * lu(k,99)
         lu(k,105) = lu(k,105) - lu(k,47) * lu(k,99)
         lu(k,172) = lu(k,172) - lu(k,44) * lu(k,162)
         lu(k,175) = lu(k,175) - lu(k,45) * lu(k,162)
         lu(k,178) = lu(k,178) - lu(k,46) * lu(k,162)
         lu(k,180) = lu(k,180) - lu(k,47) * lu(k,162)
         lu(k,209) = lu(k,209) - lu(k,44) * lu(k,202)
         lu(k,212) = - lu(k,45) * lu(k,202)
         lu(k,215) = lu(k,215) - lu(k,46) * lu(k,202)
         lu(k,217) = lu(k,217) - lu(k,47) * lu(k,202)
                                                                        
         lu(k,48) = 1. / lu(k,48)
         lu(k,49) = lu(k,49) * lu(k,48)
         lu(k,50) = lu(k,50) * lu(k,48)
         lu(k,51) = lu(k,51) * lu(k,48)
         lu(k,52) = lu(k,52) * lu(k,48)
         lu(k,117) = lu(k,117) - lu(k,49) * lu(k,116)
         lu(k,119) = lu(k,119) - lu(k,50) * lu(k,116)
         lu(k,122) = - lu(k,51) * lu(k,116)
         lu(k,124) = lu(k,124) - lu(k,52) * lu(k,116)
         lu(k,169) = lu(k,169) - lu(k,49) * lu(k,163)
         lu(k,174) = lu(k,174) - lu(k,50) * lu(k,163)
         lu(k,178) = lu(k,178) - lu(k,51) * lu(k,163)
         lu(k,180) = lu(k,180) - lu(k,52) * lu(k,163)
         lu(k,207) = - lu(k,49) * lu(k,203)
         lu(k,211) = lu(k,211) - lu(k,50) * lu(k,203)
         lu(k,215) = lu(k,215) - lu(k,51) * lu(k,203)
         lu(k,217) = lu(k,217) - lu(k,52) * lu(k,203)
                                                                        
      end do
                                                                        
      end subroutine IMP_LU_FAC01
                                                                        
      subroutine IMP_LU_FAC02( lu )
                                                                        
      use CHEM_MODS_MOD, only : imp_nzcnt, clsze
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real, intent(inout) ::   lu(clsze,imp_nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,clsze
         lu(k,53) = 1. / lu(k,53)
         lu(k,54) = lu(k,54) * lu(k,53)
         lu(k,55) = lu(k,55) * lu(k,53)
         lu(k,56) = lu(k,56) * lu(k,53)
         lu(k,57) = lu(k,57) * lu(k,53)
         lu(k,109) = lu(k,109) - lu(k,54) * lu(k,108)
         lu(k,110) = lu(k,110) - lu(k,55) * lu(k,108)
         lu(k,111) = - lu(k,56) * lu(k,108)
         lu(k,115) = lu(k,115) - lu(k,57) * lu(k,108)
         lu(k,173) = lu(k,173) - lu(k,54) * lu(k,164)
         lu(k,176) = lu(k,176) - lu(k,55) * lu(k,164)
         lu(k,178) = lu(k,178) - lu(k,56) * lu(k,164)
         lu(k,184) = lu(k,184) - lu(k,57) * lu(k,164)
         lu(k,210) = lu(k,210) - lu(k,54) * lu(k,204)
         lu(k,213) = lu(k,213) - lu(k,55) * lu(k,204)
         lu(k,215) = lu(k,215) - lu(k,56) * lu(k,204)
         lu(k,221) = lu(k,221) - lu(k,57) * lu(k,204)
                                                                        
         lu(k,58) = 1. / lu(k,58)
         lu(k,59) = lu(k,59) * lu(k,58)
         lu(k,60) = lu(k,60) * lu(k,58)
         lu(k,61) = lu(k,61) * lu(k,58)
         lu(k,62) = lu(k,62) * lu(k,58)
         lu(k,176) = lu(k,176) - lu(k,59) * lu(k,165)
         lu(k,178) = lu(k,178) - lu(k,60) * lu(k,165)
         lu(k,180) = lu(k,180) - lu(k,61) * lu(k,165)
         lu(k,181) = lu(k,181) - lu(k,62) * lu(k,165)
         lu(k,213) = lu(k,213) - lu(k,59) * lu(k,205)
         lu(k,215) = lu(k,215) - lu(k,60) * lu(k,205)
         lu(k,217) = lu(k,217) - lu(k,61) * lu(k,205)
         lu(k,218) = lu(k,218) - lu(k,62) * lu(k,205)
         lu(k,228) = lu(k,228) - lu(k,59) * lu(k,222)
         lu(k,229) = - lu(k,60) * lu(k,222)
         lu(k,231) = lu(k,231) - lu(k,61) * lu(k,222)
         lu(k,232) = lu(k,232) - lu(k,62) * lu(k,222)
                                                                        
         lu(k,63) = 1. / lu(k,63)
         lu(k,64) = lu(k,64) * lu(k,63)
         lu(k,65) = lu(k,65) * lu(k,63)
         lu(k,66) = lu(k,66) * lu(k,63)
         lu(k,130) = lu(k,130) - lu(k,64) * lu(k,128)
         lu(k,131) = - lu(k,65) * lu(k,128)
         lu(k,134) = lu(k,134) - lu(k,66) * lu(k,128)
         lu(k,138) = lu(k,138) - lu(k,64) * lu(k,136)
         lu(k,139) = - lu(k,65) * lu(k,136)
         lu(k,141) = lu(k,141) - lu(k,66) * lu(k,136)
         lu(k,178) = lu(k,178) - lu(k,64) * lu(k,166)
         lu(k,179) = lu(k,179) - lu(k,65) * lu(k,166)
         lu(k,183) = lu(k,183) - lu(k,66) * lu(k,166)
         lu(k,192) = lu(k,192) - lu(k,64) * lu(k,187)
         lu(k,193) = lu(k,193) - lu(k,65) * lu(k,187)
         lu(k,197) = lu(k,197) - lu(k,66) * lu(k,187)
         lu(k,260) = lu(k,260) - lu(k,64) * lu(k,255)
         lu(k,261) = lu(k,261) - lu(k,65) * lu(k,255)
         lu(k,265) = lu(k,265) - lu(k,66) * lu(k,255)
                                                                        
         lu(k,67) = 1. / lu(k,67)
         lu(k,68) = lu(k,68) * lu(k,67)
         lu(k,69) = lu(k,69) * lu(k,67)
         lu(k,70) = lu(k,70) * lu(k,67)
         lu(k,71) = lu(k,71) * lu(k,67)
         lu(k,94) = lu(k,94) - lu(k,68) * lu(k,90)
         lu(k,95) = lu(k,95) - lu(k,69) * lu(k,90)
         lu(k,96) = - lu(k,70) * lu(k,90)
         lu(k,98) = lu(k,98) - lu(k,71) * lu(k,90)
         lu(k,179) = lu(k,179) - lu(k,68) * lu(k,167)
         lu(k,180) = lu(k,180) - lu(k,69) * lu(k,167)
         lu(k,182) = - lu(k,70) * lu(k,167)
         lu(k,184) = lu(k,184) - lu(k,71) * lu(k,167)
         lu(k,216) = lu(k,216) - lu(k,68) * lu(k,206)
         lu(k,217) = lu(k,217) - lu(k,69) * lu(k,206)
         lu(k,219) = lu(k,219) - lu(k,70) * lu(k,206)
         lu(k,221) = lu(k,221) - lu(k,71) * lu(k,206)
         lu(k,245) = lu(k,245) - lu(k,68) * lu(k,236)
         lu(k,246) = lu(k,246) - lu(k,69) * lu(k,236)
         lu(k,248) = lu(k,248) - lu(k,70) * lu(k,236)
         lu(k,250) = lu(k,250) - lu(k,71) * lu(k,236)
                                                                        
         lu(k,72) = 1. / lu(k,72)
         lu(k,73) = lu(k,73) * lu(k,72)
         lu(k,74) = lu(k,74) * lu(k,72)
         lu(k,75) = lu(k,75) * lu(k,72)
         lu(k,76) = lu(k,76) * lu(k,72)
         lu(k,77) = lu(k,77) * lu(k,72)
         lu(k,78) = lu(k,78) * lu(k,72)
         lu(k,176) = lu(k,176) - lu(k,73) * lu(k,168)
         lu(k,178) = lu(k,178) - lu(k,74) * lu(k,168)
         lu(k,179) = lu(k,179) - lu(k,75) * lu(k,168)
         lu(k,181) = lu(k,181) - lu(k,76) * lu(k,168)
         lu(k,183) = lu(k,183) - lu(k,77) * lu(k,168)
         lu(k,184) = lu(k,184) - lu(k,78) * lu(k,168)
         lu(k,190) = - lu(k,73) * lu(k,188)
         lu(k,192) = lu(k,192) - lu(k,74) * lu(k,188)
         lu(k,193) = lu(k,193) - lu(k,75) * lu(k,188)
         lu(k,195) = - lu(k,76) * lu(k,188)
         lu(k,197) = lu(k,197) - lu(k,77) * lu(k,188)
         lu(k,198) = lu(k,198) - lu(k,78) * lu(k,188)
         lu(k,269) = lu(k,269) - lu(k,73) * lu(k,268)
         lu(k,271) = lu(k,271) - lu(k,74) * lu(k,268)
         lu(k,272) = lu(k,272) - lu(k,75) * lu(k,268)
         lu(k,274) = lu(k,274) - lu(k,76) * lu(k,268)
         lu(k,276) = - lu(k,77) * lu(k,268)
         lu(k,277) = lu(k,277) - lu(k,78) * lu(k,268)
                                                                        
         lu(k,79) = 1. / lu(k,79)
         lu(k,80) = lu(k,80) * lu(k,79)
         lu(k,81) = lu(k,81) * lu(k,79)
         lu(k,82) = lu(k,82) * lu(k,79)
         lu(k,83) = lu(k,83) * lu(k,79)
         lu(k,118) = - lu(k,80) * lu(k,117)
         lu(k,122) = lu(k,122) - lu(k,81) * lu(k,117)
         lu(k,125) = lu(k,125) - lu(k,82) * lu(k,117)
         lu(k,127) = - lu(k,83) * lu(k,117)
         lu(k,173) = lu(k,173) - lu(k,80) * lu(k,169)
         lu(k,178) = lu(k,178) - lu(k,81) * lu(k,169)
         lu(k,181) = lu(k,181) - lu(k,82) * lu(k,169)
         lu(k,184) = lu(k,184) - lu(k,83) * lu(k,169)
         lu(k,210) = lu(k,210) - lu(k,80) * lu(k,207)
         lu(k,215) = lu(k,215) - lu(k,81) * lu(k,207)
         lu(k,218) = lu(k,218) - lu(k,82) * lu(k,207)
         lu(k,221) = lu(k,221) - lu(k,83) * lu(k,207)
         lu(k,225) = - lu(k,80) * lu(k,223)
         lu(k,229) = lu(k,229) - lu(k,81) * lu(k,223)
         lu(k,232) = lu(k,232) - lu(k,82) * lu(k,223)
         lu(k,235) = lu(k,235) - lu(k,83) * lu(k,223)
         lu(k,239) = lu(k,239) - lu(k,80) * lu(k,237)
         lu(k,244) = lu(k,244) - lu(k,81) * lu(k,237)
         lu(k,247) = lu(k,247) - lu(k,82) * lu(k,237)
         lu(k,250) = lu(k,250) - lu(k,83) * lu(k,237)
                                                                        
         lu(k,84) = 1. / lu(k,84)
         lu(k,85) = lu(k,85) * lu(k,84)
         lu(k,86) = lu(k,86) * lu(k,84)
         lu(k,87) = lu(k,87) * lu(k,84)
         lu(k,88) = lu(k,88) * lu(k,84)
         lu(k,89) = lu(k,89) * lu(k,84)
         lu(k,146) = lu(k,146) - lu(k,85) * lu(k,143)
         lu(k,147) = lu(k,147) - lu(k,86) * lu(k,143)
         lu(k,148) = lu(k,148) - lu(k,87) * lu(k,143)
         lu(k,149) = lu(k,149) - lu(k,88) * lu(k,143)
         lu(k,151) = lu(k,151) - lu(k,89) * lu(k,143)
         lu(k,177) = lu(k,177) - lu(k,85) * lu(k,170)
         lu(k,178) = lu(k,178) - lu(k,86) * lu(k,170)
         lu(k,179) = lu(k,179) - lu(k,87) * lu(k,170)
         lu(k,180) = lu(k,180) - lu(k,88) * lu(k,170)
         lu(k,182) = lu(k,182) - lu(k,89) * lu(k,170)
         lu(k,191) = lu(k,191) - lu(k,85) * lu(k,189)
         lu(k,192) = lu(k,192) - lu(k,86) * lu(k,189)
         lu(k,193) = lu(k,193) - lu(k,87) * lu(k,189)
         lu(k,194) = lu(k,194) - lu(k,88) * lu(k,189)
         lu(k,196) = lu(k,196) - lu(k,89) * lu(k,189)
         lu(k,214) = lu(k,214) - lu(k,85) * lu(k,208)
         lu(k,215) = lu(k,215) - lu(k,86) * lu(k,208)
         lu(k,216) = lu(k,216) - lu(k,87) * lu(k,208)
         lu(k,217) = lu(k,217) - lu(k,88) * lu(k,208)
         lu(k,219) = lu(k,219) - lu(k,89) * lu(k,208)
                                                                        
         lu(k,91) = 1. / lu(k,91)
         lu(k,92) = lu(k,92) * lu(k,91)
         lu(k,93) = lu(k,93) * lu(k,91)
         lu(k,94) = lu(k,94) * lu(k,91)
         lu(k,95) = lu(k,95) * lu(k,91)
         lu(k,96) = lu(k,96) * lu(k,91)
         lu(k,97) = lu(k,97) * lu(k,91)
         lu(k,98) = lu(k,98) * lu(k,91)
         lu(k,146) = lu(k,146) - lu(k,92) * lu(k,144)
         lu(k,147) = lu(k,147) - lu(k,93) * lu(k,144)
         lu(k,148) = lu(k,148) - lu(k,94) * lu(k,144)
         lu(k,149) = lu(k,149) - lu(k,95) * lu(k,144)
         lu(k,151) = lu(k,151) - lu(k,96) * lu(k,144)
         lu(k,152) = lu(k,152) - lu(k,97) * lu(k,144)
         lu(k,153) = lu(k,153) - lu(k,98) * lu(k,144)
         lu(k,177) = lu(k,177) - lu(k,92) * lu(k,171)
         lu(k,178) = lu(k,178) - lu(k,93) * lu(k,171)
         lu(k,179) = lu(k,179) - lu(k,94) * lu(k,171)
         lu(k,180) = lu(k,180) - lu(k,95) * lu(k,171)
         lu(k,182) = lu(k,182) - lu(k,96) * lu(k,171)
         lu(k,183) = lu(k,183) - lu(k,97) * lu(k,171)
         lu(k,184) = lu(k,184) - lu(k,98) * lu(k,171)
         lu(k,259) = lu(k,259) - lu(k,92) * lu(k,256)
         lu(k,260) = lu(k,260) - lu(k,93) * lu(k,256)
         lu(k,261) = lu(k,261) - lu(k,94) * lu(k,256)
         lu(k,262) = lu(k,262) - lu(k,95) * lu(k,256)
         lu(k,264) = lu(k,264) - lu(k,96) * lu(k,256)
         lu(k,265) = lu(k,265) - lu(k,97) * lu(k,256)
         lu(k,266) = lu(k,266) - lu(k,98) * lu(k,256)
                                                                        
      end do
                                                                        
      end subroutine IMP_LU_FAC02
                                                                        
      subroutine IMP_LU_FAC03( lu )
                                                                        
      use CHEM_MODS_MOD, only : imp_nzcnt, clsze
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real, intent(inout) ::   lu(clsze,imp_nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,clsze
         lu(k,100) = 1. / lu(k,100)
         lu(k,101) = lu(k,101) * lu(k,100)
         lu(k,102) = lu(k,102) * lu(k,100)
         lu(k,103) = lu(k,103) * lu(k,100)
         lu(k,104) = lu(k,104) * lu(k,100)
         lu(k,105) = lu(k,105) * lu(k,100)
         lu(k,106) = lu(k,106) * lu(k,100)
         lu(k,107) = lu(k,107) * lu(k,100)
         lu(k,175) = lu(k,175) - lu(k,101) * lu(k,172)
         lu(k,176) = lu(k,176) - lu(k,102) * lu(k,172)
         lu(k,178) = lu(k,178) - lu(k,103) * lu(k,172)
         lu(k,179) = lu(k,179) - lu(k,104) * lu(k,172)
         lu(k,180) = lu(k,180) - lu(k,105) * lu(k,172)
         lu(k,181) = lu(k,181) - lu(k,106) * lu(k,172)
         lu(k,182) = lu(k,182) - lu(k,107) * lu(k,172)
         lu(k,212) = lu(k,212) - lu(k,101) * lu(k,209)
         lu(k,213) = lu(k,213) - lu(k,102) * lu(k,209)
         lu(k,215) = lu(k,215) - lu(k,103) * lu(k,209)
         lu(k,216) = lu(k,216) - lu(k,104) * lu(k,209)
         lu(k,217) = lu(k,217) - lu(k,105) * lu(k,209)
         lu(k,218) = lu(k,218) - lu(k,106) * lu(k,209)
         lu(k,219) = lu(k,219) - lu(k,107) * lu(k,209)
         lu(k,227) = lu(k,227) - lu(k,101) * lu(k,224)
         lu(k,228) = lu(k,228) - lu(k,102) * lu(k,224)
         lu(k,229) = lu(k,229) - lu(k,103) * lu(k,224)
         lu(k,230) = lu(k,230) - lu(k,104) * lu(k,224)
         lu(k,231) = lu(k,231) - lu(k,105) * lu(k,224)
         lu(k,232) = lu(k,232) - lu(k,106) * lu(k,224)
         lu(k,233) = lu(k,233) - lu(k,107) * lu(k,224)
         lu(k,241) = lu(k,241) - lu(k,101) * lu(k,238)
         lu(k,242) = lu(k,242) - lu(k,102) * lu(k,238)
         lu(k,244) = lu(k,244) - lu(k,103) * lu(k,238)
         lu(k,245) = lu(k,245) - lu(k,104) * lu(k,238)
         lu(k,246) = lu(k,246) - lu(k,105) * lu(k,238)
         lu(k,247) = lu(k,247) - lu(k,106) * lu(k,238)
         lu(k,248) = lu(k,248) - lu(k,107) * lu(k,238)
                                                                        
         lu(k,109) = 1. / lu(k,109)
         lu(k,110) = lu(k,110) * lu(k,109)
         lu(k,111) = lu(k,111) * lu(k,109)
         lu(k,112) = lu(k,112) * lu(k,109)
         lu(k,113) = lu(k,113) * lu(k,109)
         lu(k,114) = lu(k,114) * lu(k,109)
         lu(k,115) = lu(k,115) * lu(k,109)
         lu(k,121) = lu(k,121) - lu(k,110) * lu(k,118)
         lu(k,122) = lu(k,122) - lu(k,111) * lu(k,118)
         lu(k,123) = lu(k,123) - lu(k,112) * lu(k,118)
         lu(k,124) = lu(k,124) - lu(k,113) * lu(k,118)
         lu(k,126) = lu(k,126) - lu(k,114) * lu(k,118)
         lu(k,127) = lu(k,127) - lu(k,115) * lu(k,118)
         lu(k,176) = lu(k,176) - lu(k,110) * lu(k,173)
         lu(k,178) = lu(k,178) - lu(k,111) * lu(k,173)
         lu(k,179) = lu(k,179) - lu(k,112) * lu(k,173)
         lu(k,180) = lu(k,180) - lu(k,113) * lu(k,173)
         lu(k,182) = lu(k,182) - lu(k,114) * lu(k,173)
         lu(k,184) = lu(k,184) - lu(k,115) * lu(k,173)
         lu(k,213) = lu(k,213) - lu(k,110) * lu(k,210)
         lu(k,215) = lu(k,215) - lu(k,111) * lu(k,210)
         lu(k,216) = lu(k,216) - lu(k,112) * lu(k,210)
         lu(k,217) = lu(k,217) - lu(k,113) * lu(k,210)
         lu(k,219) = lu(k,219) - lu(k,114) * lu(k,210)
         lu(k,221) = lu(k,221) - lu(k,115) * lu(k,210)
         lu(k,228) = lu(k,228) - lu(k,110) * lu(k,225)
         lu(k,229) = lu(k,229) - lu(k,111) * lu(k,225)
         lu(k,230) = lu(k,230) - lu(k,112) * lu(k,225)
         lu(k,231) = lu(k,231) - lu(k,113) * lu(k,225)
         lu(k,233) = lu(k,233) - lu(k,114) * lu(k,225)
         lu(k,235) = lu(k,235) - lu(k,115) * lu(k,225)
         lu(k,242) = lu(k,242) - lu(k,110) * lu(k,239)
         lu(k,244) = lu(k,244) - lu(k,111) * lu(k,239)
         lu(k,245) = lu(k,245) - lu(k,112) * lu(k,239)
         lu(k,246) = lu(k,246) - lu(k,113) * lu(k,239)
         lu(k,248) = lu(k,248) - lu(k,114) * lu(k,239)
         lu(k,250) = lu(k,250) - lu(k,115) * lu(k,239)
                                                                        
         lu(k,119) = 1. / lu(k,119)
         lu(k,120) = lu(k,120) * lu(k,119)
         lu(k,121) = lu(k,121) * lu(k,119)
         lu(k,122) = lu(k,122) * lu(k,119)
         lu(k,123) = lu(k,123) * lu(k,119)
         lu(k,124) = lu(k,124) * lu(k,119)
         lu(k,125) = lu(k,125) * lu(k,119)
         lu(k,126) = lu(k,126) * lu(k,119)
         lu(k,127) = lu(k,127) * lu(k,119)
         lu(k,175) = lu(k,175) - lu(k,120) * lu(k,174)
         lu(k,176) = lu(k,176) - lu(k,121) * lu(k,174)
         lu(k,178) = lu(k,178) - lu(k,122) * lu(k,174)
         lu(k,179) = lu(k,179) - lu(k,123) * lu(k,174)
         lu(k,180) = lu(k,180) - lu(k,124) * lu(k,174)
         lu(k,181) = lu(k,181) - lu(k,125) * lu(k,174)
         lu(k,182) = lu(k,182) - lu(k,126) * lu(k,174)
         lu(k,184) = lu(k,184) - lu(k,127) * lu(k,174)
         lu(k,212) = lu(k,212) - lu(k,120) * lu(k,211)
         lu(k,213) = lu(k,213) - lu(k,121) * lu(k,211)
         lu(k,215) = lu(k,215) - lu(k,122) * lu(k,211)
         lu(k,216) = lu(k,216) - lu(k,123) * lu(k,211)
         lu(k,217) = lu(k,217) - lu(k,124) * lu(k,211)
         lu(k,218) = lu(k,218) - lu(k,125) * lu(k,211)
         lu(k,219) = lu(k,219) - lu(k,126) * lu(k,211)
         lu(k,221) = lu(k,221) - lu(k,127) * lu(k,211)
         lu(k,227) = lu(k,227) - lu(k,120) * lu(k,226)
         lu(k,228) = lu(k,228) - lu(k,121) * lu(k,226)
         lu(k,229) = lu(k,229) - lu(k,122) * lu(k,226)
         lu(k,230) = lu(k,230) - lu(k,123) * lu(k,226)
         lu(k,231) = lu(k,231) - lu(k,124) * lu(k,226)
         lu(k,232) = lu(k,232) - lu(k,125) * lu(k,226)
         lu(k,233) = lu(k,233) - lu(k,126) * lu(k,226)
         lu(k,235) = lu(k,235) - lu(k,127) * lu(k,226)
         lu(k,241) = lu(k,241) - lu(k,120) * lu(k,240)
         lu(k,242) = lu(k,242) - lu(k,121) * lu(k,240)
         lu(k,244) = lu(k,244) - lu(k,122) * lu(k,240)
         lu(k,245) = lu(k,245) - lu(k,123) * lu(k,240)
         lu(k,246) = lu(k,246) - lu(k,124) * lu(k,240)
         lu(k,247) = lu(k,247) - lu(k,125) * lu(k,240)
         lu(k,248) = lu(k,248) - lu(k,126) * lu(k,240)
         lu(k,250) = lu(k,250) - lu(k,127) * lu(k,240)
                                                                        
         lu(k,129) = 1. / lu(k,129)
         lu(k,130) = lu(k,130) * lu(k,129)
         lu(k,131) = lu(k,131) * lu(k,129)
         lu(k,132) = lu(k,132) * lu(k,129)
         lu(k,133) = lu(k,133) * lu(k,129)
         lu(k,134) = lu(k,134) * lu(k,129)
         lu(k,135) = lu(k,135) * lu(k,129)
         lu(k,178) = lu(k,178) - lu(k,130) * lu(k,175)
         lu(k,179) = lu(k,179) - lu(k,131) * lu(k,175)
         lu(k,180) = lu(k,180) - lu(k,132) * lu(k,175)
         lu(k,181) = lu(k,181) - lu(k,133) * lu(k,175)
         lu(k,183) = lu(k,183) - lu(k,134) * lu(k,175)
         lu(k,184) = lu(k,184) - lu(k,135) * lu(k,175)
         lu(k,215) = lu(k,215) - lu(k,130) * lu(k,212)
         lu(k,216) = lu(k,216) - lu(k,131) * lu(k,212)
         lu(k,217) = lu(k,217) - lu(k,132) * lu(k,212)
         lu(k,218) = lu(k,218) - lu(k,133) * lu(k,212)
         lu(k,220) = lu(k,220) - lu(k,134) * lu(k,212)
         lu(k,221) = lu(k,221) - lu(k,135) * lu(k,212)
         lu(k,229) = lu(k,229) - lu(k,130) * lu(k,227)
         lu(k,230) = lu(k,230) - lu(k,131) * lu(k,227)
         lu(k,231) = lu(k,231) - lu(k,132) * lu(k,227)
         lu(k,232) = lu(k,232) - lu(k,133) * lu(k,227)
         lu(k,234) = - lu(k,134) * lu(k,227)
         lu(k,235) = lu(k,235) - lu(k,135) * lu(k,227)
         lu(k,244) = lu(k,244) - lu(k,130) * lu(k,241)
         lu(k,245) = lu(k,245) - lu(k,131) * lu(k,241)
         lu(k,246) = lu(k,246) - lu(k,132) * lu(k,241)
         lu(k,247) = lu(k,247) - lu(k,133) * lu(k,241)
         lu(k,249) = lu(k,249) - lu(k,134) * lu(k,241)
         lu(k,250) = lu(k,250) - lu(k,135) * lu(k,241)
         lu(k,260) = lu(k,260) - lu(k,130) * lu(k,257)
         lu(k,261) = lu(k,261) - lu(k,131) * lu(k,257)
         lu(k,262) = lu(k,262) - lu(k,132) * lu(k,257)
         lu(k,263) = - lu(k,133) * lu(k,257)
         lu(k,265) = lu(k,265) - lu(k,134) * lu(k,257)
         lu(k,266) = lu(k,266) - lu(k,135) * lu(k,257)
                                                                        
         lu(k,137) = 1. / lu(k,137)
         lu(k,138) = lu(k,138) * lu(k,137)
         lu(k,139) = lu(k,139) * lu(k,137)
         lu(k,140) = lu(k,140) * lu(k,137)
         lu(k,141) = lu(k,141) * lu(k,137)
         lu(k,147) = lu(k,147) - lu(k,138) * lu(k,145)
         lu(k,148) = lu(k,148) - lu(k,139) * lu(k,145)
         lu(k,149) = lu(k,149) - lu(k,140) * lu(k,145)
         lu(k,152) = lu(k,152) - lu(k,141) * lu(k,145)
         lu(k,178) = lu(k,178) - lu(k,138) * lu(k,176)
         lu(k,179) = lu(k,179) - lu(k,139) * lu(k,176)
         lu(k,180) = lu(k,180) - lu(k,140) * lu(k,176)
         lu(k,183) = lu(k,183) - lu(k,141) * lu(k,176)
         lu(k,192) = lu(k,192) - lu(k,138) * lu(k,190)
         lu(k,193) = lu(k,193) - lu(k,139) * lu(k,190)
         lu(k,194) = lu(k,194) - lu(k,140) * lu(k,190)
         lu(k,197) = lu(k,197) - lu(k,141) * lu(k,190)
         lu(k,215) = lu(k,215) - lu(k,138) * lu(k,213)
         lu(k,216) = lu(k,216) - lu(k,139) * lu(k,213)
         lu(k,217) = lu(k,217) - lu(k,140) * lu(k,213)
         lu(k,220) = lu(k,220) - lu(k,141) * lu(k,213)
         lu(k,229) = lu(k,229) - lu(k,138) * lu(k,228)
         lu(k,230) = lu(k,230) - lu(k,139) * lu(k,228)
         lu(k,231) = lu(k,231) - lu(k,140) * lu(k,228)
         lu(k,234) = lu(k,234) - lu(k,141) * lu(k,228)
         lu(k,244) = lu(k,244) - lu(k,138) * lu(k,242)
         lu(k,245) = lu(k,245) - lu(k,139) * lu(k,242)
         lu(k,246) = lu(k,246) - lu(k,140) * lu(k,242)
         lu(k,249) = lu(k,249) - lu(k,141) * lu(k,242)
         lu(k,260) = lu(k,260) - lu(k,138) * lu(k,258)
         lu(k,261) = lu(k,261) - lu(k,139) * lu(k,258)
         lu(k,262) = lu(k,262) - lu(k,140) * lu(k,258)
         lu(k,265) = lu(k,265) - lu(k,141) * lu(k,258)
         lu(k,271) = lu(k,271) - lu(k,138) * lu(k,269)
         lu(k,272) = lu(k,272) - lu(k,139) * lu(k,269)
         lu(k,273) = lu(k,273) - lu(k,140) * lu(k,269)
         lu(k,276) = lu(k,276) - lu(k,141) * lu(k,269)
                                                                        
         lu(k,146) = 1. / lu(k,146)
         lu(k,147) = lu(k,147) * lu(k,146)
         lu(k,148) = lu(k,148) * lu(k,146)
         lu(k,149) = lu(k,149) * lu(k,146)
         lu(k,150) = lu(k,150) * lu(k,146)
         lu(k,151) = lu(k,151) * lu(k,146)
         lu(k,152) = lu(k,152) * lu(k,146)
         lu(k,153) = lu(k,153) * lu(k,146)
         lu(k,178) = lu(k,178) - lu(k,147) * lu(k,177)
         lu(k,179) = lu(k,179) - lu(k,148) * lu(k,177)
         lu(k,180) = lu(k,180) - lu(k,149) * lu(k,177)
         lu(k,181) = lu(k,181) - lu(k,150) * lu(k,177)
         lu(k,182) = lu(k,182) - lu(k,151) * lu(k,177)
         lu(k,183) = lu(k,183) - lu(k,152) * lu(k,177)
         lu(k,184) = lu(k,184) - lu(k,153) * lu(k,177)
         lu(k,192) = lu(k,192) - lu(k,147) * lu(k,191)
         lu(k,193) = lu(k,193) - lu(k,148) * lu(k,191)
         lu(k,194) = lu(k,194) - lu(k,149) * lu(k,191)
         lu(k,195) = lu(k,195) - lu(k,150) * lu(k,191)
         lu(k,196) = lu(k,196) - lu(k,151) * lu(k,191)
         lu(k,197) = lu(k,197) - lu(k,152) * lu(k,191)
         lu(k,198) = lu(k,198) - lu(k,153) * lu(k,191)
         lu(k,215) = lu(k,215) - lu(k,147) * lu(k,214)
         lu(k,216) = lu(k,216) - lu(k,148) * lu(k,214)
         lu(k,217) = lu(k,217) - lu(k,149) * lu(k,214)
         lu(k,218) = lu(k,218) - lu(k,150) * lu(k,214)
         lu(k,219) = lu(k,219) - lu(k,151) * lu(k,214)
         lu(k,220) = lu(k,220) - lu(k,152) * lu(k,214)
         lu(k,221) = lu(k,221) - lu(k,153) * lu(k,214)
         lu(k,244) = lu(k,244) - lu(k,147) * lu(k,243)
         lu(k,245) = lu(k,245) - lu(k,148) * lu(k,243)
         lu(k,246) = lu(k,246) - lu(k,149) * lu(k,243)
         lu(k,247) = lu(k,247) - lu(k,150) * lu(k,243)
         lu(k,248) = lu(k,248) - lu(k,151) * lu(k,243)
         lu(k,249) = lu(k,249) - lu(k,152) * lu(k,243)
         lu(k,250) = lu(k,250) - lu(k,153) * lu(k,243)
         lu(k,260) = lu(k,260) - lu(k,147) * lu(k,259)
         lu(k,261) = lu(k,261) - lu(k,148) * lu(k,259)
         lu(k,262) = lu(k,262) - lu(k,149) * lu(k,259)
         lu(k,263) = lu(k,263) - lu(k,150) * lu(k,259)
         lu(k,264) = lu(k,264) - lu(k,151) * lu(k,259)
         lu(k,265) = lu(k,265) - lu(k,152) * lu(k,259)
         lu(k,266) = lu(k,266) - lu(k,153) * lu(k,259)
         lu(k,271) = lu(k,271) - lu(k,147) * lu(k,270)
         lu(k,272) = lu(k,272) - lu(k,148) * lu(k,270)
         lu(k,273) = lu(k,273) - lu(k,149) * lu(k,270)
         lu(k,274) = lu(k,274) - lu(k,150) * lu(k,270)
         lu(k,275) = lu(k,275) - lu(k,151) * lu(k,270)
         lu(k,276) = lu(k,276) - lu(k,152) * lu(k,270)
         lu(k,277) = lu(k,277) - lu(k,153) * lu(k,270)
                                                                        
         lu(k,178) = 1. / lu(k,178)
         lu(k,179) = lu(k,179) * lu(k,178)
         lu(k,180) = lu(k,180) * lu(k,178)
         lu(k,181) = lu(k,181) * lu(k,178)
         lu(k,182) = lu(k,182) * lu(k,178)
         lu(k,183) = lu(k,183) * lu(k,178)
         lu(k,184) = lu(k,184) * lu(k,178)
         lu(k,193) = lu(k,193) - lu(k,179) * lu(k,192)
         lu(k,194) = lu(k,194) - lu(k,180) * lu(k,192)
         lu(k,195) = lu(k,195) - lu(k,181) * lu(k,192)
         lu(k,196) = lu(k,196) - lu(k,182) * lu(k,192)
         lu(k,197) = lu(k,197) - lu(k,183) * lu(k,192)
         lu(k,198) = lu(k,198) - lu(k,184) * lu(k,192)
         lu(k,216) = lu(k,216) - lu(k,179) * lu(k,215)
         lu(k,217) = lu(k,217) - lu(k,180) * lu(k,215)
         lu(k,218) = lu(k,218) - lu(k,181) * lu(k,215)
         lu(k,219) = lu(k,219) - lu(k,182) * lu(k,215)
         lu(k,220) = lu(k,220) - lu(k,183) * lu(k,215)
         lu(k,221) = lu(k,221) - lu(k,184) * lu(k,215)
         lu(k,230) = lu(k,230) - lu(k,179) * lu(k,229)
         lu(k,231) = lu(k,231) - lu(k,180) * lu(k,229)
         lu(k,232) = lu(k,232) - lu(k,181) * lu(k,229)
         lu(k,233) = lu(k,233) - lu(k,182) * lu(k,229)
         lu(k,234) = lu(k,234) - lu(k,183) * lu(k,229)
         lu(k,235) = lu(k,235) - lu(k,184) * lu(k,229)
         lu(k,245) = lu(k,245) - lu(k,179) * lu(k,244)
         lu(k,246) = lu(k,246) - lu(k,180) * lu(k,244)
         lu(k,247) = lu(k,247) - lu(k,181) * lu(k,244)
         lu(k,248) = lu(k,248) - lu(k,182) * lu(k,244)
         lu(k,249) = lu(k,249) - lu(k,183) * lu(k,244)
         lu(k,250) = lu(k,250) - lu(k,184) * lu(k,244)
         lu(k,261) = lu(k,261) - lu(k,179) * lu(k,260)
         lu(k,262) = lu(k,262) - lu(k,180) * lu(k,260)
         lu(k,263) = lu(k,263) - lu(k,181) * lu(k,260)
         lu(k,264) = lu(k,264) - lu(k,182) * lu(k,260)
         lu(k,265) = lu(k,265) - lu(k,183) * lu(k,260)
         lu(k,266) = lu(k,266) - lu(k,184) * lu(k,260)
         lu(k,272) = lu(k,272) - lu(k,179) * lu(k,271)
         lu(k,273) = lu(k,273) - lu(k,180) * lu(k,271)
         lu(k,274) = lu(k,274) - lu(k,181) * lu(k,271)
         lu(k,275) = lu(k,275) - lu(k,182) * lu(k,271)
         lu(k,276) = lu(k,276) - lu(k,183) * lu(k,271)
         lu(k,277) = lu(k,277) - lu(k,184) * lu(k,271)
                                                                        
      end do
                                                                        
      end subroutine IMP_LU_FAC03
                                                                        
      subroutine IMP_LU_FAC04( lu )
                                                                        
      use CHEM_MODS_MOD, only : imp_nzcnt, clsze
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real, intent(inout) ::   lu(clsze,imp_nzcnt)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
      do k = 1,clsze
         lu(k,193) = 1. / lu(k,193)
         lu(k,194) = lu(k,194) * lu(k,193)
         lu(k,195) = lu(k,195) * lu(k,193)
         lu(k,196) = lu(k,196) * lu(k,193)
         lu(k,197) = lu(k,197) * lu(k,193)
         lu(k,198) = lu(k,198) * lu(k,193)
         lu(k,217) = lu(k,217) - lu(k,194) * lu(k,216)
         lu(k,218) = lu(k,218) - lu(k,195) * lu(k,216)
         lu(k,219) = lu(k,219) - lu(k,196) * lu(k,216)
         lu(k,220) = lu(k,220) - lu(k,197) * lu(k,216)
         lu(k,221) = lu(k,221) - lu(k,198) * lu(k,216)
         lu(k,231) = lu(k,231) - lu(k,194) * lu(k,230)
         lu(k,232) = lu(k,232) - lu(k,195) * lu(k,230)
         lu(k,233) = lu(k,233) - lu(k,196) * lu(k,230)
         lu(k,234) = lu(k,234) - lu(k,197) * lu(k,230)
         lu(k,235) = lu(k,235) - lu(k,198) * lu(k,230)
         lu(k,246) = lu(k,246) - lu(k,194) * lu(k,245)
         lu(k,247) = lu(k,247) - lu(k,195) * lu(k,245)
         lu(k,248) = lu(k,248) - lu(k,196) * lu(k,245)
         lu(k,249) = lu(k,249) - lu(k,197) * lu(k,245)
         lu(k,250) = lu(k,250) - lu(k,198) * lu(k,245)
         lu(k,262) = lu(k,262) - lu(k,194) * lu(k,261)
         lu(k,263) = lu(k,263) - lu(k,195) * lu(k,261)
         lu(k,264) = lu(k,264) - lu(k,196) * lu(k,261)
         lu(k,265) = lu(k,265) - lu(k,197) * lu(k,261)
         lu(k,266) = lu(k,266) - lu(k,198) * lu(k,261)
         lu(k,273) = lu(k,273) - lu(k,194) * lu(k,272)
         lu(k,274) = lu(k,274) - lu(k,195) * lu(k,272)
         lu(k,275) = lu(k,275) - lu(k,196) * lu(k,272)
         lu(k,276) = lu(k,276) - lu(k,197) * lu(k,272)
         lu(k,277) = lu(k,277) - lu(k,198) * lu(k,272)
                                                                        
         lu(k,217) = 1. / lu(k,217)
         lu(k,218) = lu(k,218) * lu(k,217)
         lu(k,219) = lu(k,219) * lu(k,217)
         lu(k,220) = lu(k,220) * lu(k,217)
         lu(k,221) = lu(k,221) * lu(k,217)
         lu(k,232) = lu(k,232) - lu(k,218) * lu(k,231)
         lu(k,233) = lu(k,233) - lu(k,219) * lu(k,231)
         lu(k,234) = lu(k,234) - lu(k,220) * lu(k,231)
         lu(k,235) = lu(k,235) - lu(k,221) * lu(k,231)
         lu(k,247) = lu(k,247) - lu(k,218) * lu(k,246)
         lu(k,248) = lu(k,248) - lu(k,219) * lu(k,246)
         lu(k,249) = lu(k,249) - lu(k,220) * lu(k,246)
         lu(k,250) = lu(k,250) - lu(k,221) * lu(k,246)
         lu(k,263) = lu(k,263) - lu(k,218) * lu(k,262)
         lu(k,264) = lu(k,264) - lu(k,219) * lu(k,262)
         lu(k,265) = lu(k,265) - lu(k,220) * lu(k,262)
         lu(k,266) = lu(k,266) - lu(k,221) * lu(k,262)
         lu(k,274) = lu(k,274) - lu(k,218) * lu(k,273)
         lu(k,275) = lu(k,275) - lu(k,219) * lu(k,273)
         lu(k,276) = lu(k,276) - lu(k,220) * lu(k,273)
         lu(k,277) = lu(k,277) - lu(k,221) * lu(k,273)
                                                                        
         lu(k,232) = 1. / lu(k,232)
         lu(k,233) = lu(k,233) * lu(k,232)
         lu(k,234) = lu(k,234) * lu(k,232)
         lu(k,235) = lu(k,235) * lu(k,232)
         lu(k,248) = lu(k,248) - lu(k,233) * lu(k,247)
         lu(k,249) = lu(k,249) - lu(k,234) * lu(k,247)
         lu(k,250) = lu(k,250) - lu(k,235) * lu(k,247)
         lu(k,264) = lu(k,264) - lu(k,233) * lu(k,263)
         lu(k,265) = lu(k,265) - lu(k,234) * lu(k,263)
         lu(k,266) = lu(k,266) - lu(k,235) * lu(k,263)
         lu(k,275) = lu(k,275) - lu(k,233) * lu(k,274)
         lu(k,276) = lu(k,276) - lu(k,234) * lu(k,274)
         lu(k,277) = lu(k,277) - lu(k,235) * lu(k,274)
                                                                        
         lu(k,248) = 1. / lu(k,248)
         lu(k,249) = lu(k,249) * lu(k,248)
         lu(k,250) = lu(k,250) * lu(k,248)
         lu(k,265) = lu(k,265) - lu(k,249) * lu(k,264)
         lu(k,266) = lu(k,266) - lu(k,250) * lu(k,264)
         lu(k,276) = lu(k,276) - lu(k,249) * lu(k,275)
         lu(k,277) = lu(k,277) - lu(k,250) * lu(k,275)
                                                                        
         lu(k,265) = 1. / lu(k,265)
         lu(k,266) = lu(k,266) * lu(k,265)
         lu(k,277) = lu(k,277) - lu(k,266) * lu(k,276)
                                                                        
         lu(k,277) = 1. / lu(k,277)
                                                                        
      end do
                                                                        
      end subroutine IMP_LU_FAC04
                                                                        
      subroutine IMP_LU_FAC( lu )
                                                                        
      use CHEM_MODS_MOD, only : imp_nzcnt, clsze
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real, intent(inout) ::   lu(clsze,imp_nzcnt)
                                                                        
      call IMP_LU_FAC01( lu )
      call IMP_LU_FAC02( lu )
      call IMP_LU_FAC03( lu )
      call IMP_LU_FAC04( lu )
                                                                        
      end subroutine IMP_LU_FAC
                                                                        
      end module MO_IMP_FACTOR_MOD

      module MO_ROD_FACTOR_MOD

      CONTAINS
                                                                        
      subroutine ROD_LU_FAC( lu )
                                                                        
      use CHEM_MODS_MOD, only : rod_nzcnt, clsze
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real, intent(inout) ::   lu(clsze,rod_nzcnt)
                                                                        
                                                                        
      end subroutine ROD_LU_FAC
                                                                        
      end module MO_ROD_FACTOR_MOD

      module MO_IMP_SOLVE_MOD

      CONTAINS
                                                                        
      subroutine IMP_LU_SLV01( lu, b )
                                                                        
      use CHEM_MODS_MOD, only : imp_nzcnt, clsze, clscnt4
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real, intent(in)    ::   lu(clsze,imp_nzcnt)
      real, intent(inout) ::   b(clsze,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... Solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,clsze
                                                                        
                                                                        
                                                                        
         b(k,29) = b(k,29) - lu(k,6) * b(k,4)
                                                                        
         b(k,29) = b(k,29) - lu(k,9) * b(k,5)
                                                                        
         b(k,29) = b(k,29) - lu(k,11) * b(k,6)
         b(k,31) = b(k,31) - lu(k,12) * b(k,6)
                                                                        
         b(k,21) = b(k,21) - lu(k,14) * b(k,7)
         b(k,27) = b(k,27) - lu(k,15) * b(k,7)
         b(k,29) = b(k,29) - lu(k,16) * b(k,7)
         b(k,31) = b(k,31) - lu(k,17) * b(k,7)
         b(k,32) = b(k,32) - lu(k,18) * b(k,7)
                                                                        
         b(k,17) = b(k,17) - lu(k,20) * b(k,8)
         b(k,30) = b(k,30) - lu(k,21) * b(k,8)
         b(k,34) = b(k,34) - lu(k,22) * b(k,8)
                                                                        
         b(k,29) = b(k,29) - lu(k,25) * b(k,9)
         b(k,30) = b(k,30) - lu(k,26) * b(k,9)
         b(k,34) = b(k,34) - lu(k,27) * b(k,9)
                                                                        
         b(k,29) = b(k,29) - lu(k,29) * b(k,10)
         b(k,30) = b(k,30) - lu(k,30) * b(k,10)
         b(k,31) = b(k,31) - lu(k,31) * b(k,10)
                                                                        
         b(k,17) = b(k,17) - lu(k,33) * b(k,11)
         b(k,29) = b(k,29) - lu(k,34) * b(k,11)
         b(k,31) = b(k,31) - lu(k,35) * b(k,11)
         b(k,34) = b(k,34) - lu(k,36) * b(k,11)
         b(k,35) = b(k,35) - lu(k,37) * b(k,11)
                                                                        
         b(k,27) = b(k,27) - lu(k,39) * b(k,12)
         b(k,29) = b(k,29) - lu(k,40) * b(k,12)
         b(k,32) = b(k,32) - lu(k,41) * b(k,12)
         b(k,35) = b(k,35) - lu(k,42) * b(k,12)
                                                                        
         b(k,23) = b(k,23) - lu(k,44) * b(k,13)
         b(k,26) = b(k,26) - lu(k,45) * b(k,13)
         b(k,29) = b(k,29) - lu(k,46) * b(k,13)
         b(k,31) = b(k,31) - lu(k,47) * b(k,13)
                                                                        
         b(k,20) = b(k,20) - lu(k,49) * b(k,14)
         b(k,25) = b(k,25) - lu(k,50) * b(k,14)
         b(k,29) = b(k,29) - lu(k,51) * b(k,14)
         b(k,31) = b(k,31) - lu(k,52) * b(k,14)
                                                                        
         b(k,24) = b(k,24) - lu(k,54) * b(k,15)
         b(k,27) = b(k,27) - lu(k,55) * b(k,15)
         b(k,29) = b(k,29) - lu(k,56) * b(k,15)
         b(k,35) = b(k,35) - lu(k,57) * b(k,15)
                                                                        
         b(k,27) = b(k,27) - lu(k,59) * b(k,16)
         b(k,29) = b(k,29) - lu(k,60) * b(k,16)
         b(k,31) = b(k,31) - lu(k,61) * b(k,16)
         b(k,32) = b(k,32) - lu(k,62) * b(k,16)
                                                                        
         b(k,29) = b(k,29) - lu(k,64) * b(k,17)
         b(k,30) = b(k,30) - lu(k,65) * b(k,17)
         b(k,34) = b(k,34) - lu(k,66) * b(k,17)
                                                                        
         b(k,30) = b(k,30) - lu(k,68) * b(k,18)
         b(k,31) = b(k,31) - lu(k,69) * b(k,18)
         b(k,33) = b(k,33) - lu(k,70) * b(k,18)
         b(k,35) = b(k,35) - lu(k,71) * b(k,18)
                                                                        
         b(k,27) = b(k,27) - lu(k,73) * b(k,19)
         b(k,29) = b(k,29) - lu(k,74) * b(k,19)
         b(k,30) = b(k,30) - lu(k,75) * b(k,19)
         b(k,32) = b(k,32) - lu(k,76) * b(k,19)
         b(k,34) = b(k,34) - lu(k,77) * b(k,19)
         b(k,35) = b(k,35) - lu(k,78) * b(k,19)
                                                                        
         b(k,24) = b(k,24) - lu(k,80) * b(k,20)
         b(k,29) = b(k,29) - lu(k,81) * b(k,20)
         b(k,32) = b(k,32) - lu(k,82) * b(k,20)
         b(k,35) = b(k,35) - lu(k,83) * b(k,20)
                                                                        
         b(k,28) = b(k,28) - lu(k,85) * b(k,21)
         b(k,29) = b(k,29) - lu(k,86) * b(k,21)
         b(k,30) = b(k,30) - lu(k,87) * b(k,21)
         b(k,31) = b(k,31) - lu(k,88) * b(k,21)
         b(k,33) = b(k,33) - lu(k,89) * b(k,21)
                                                                        
         b(k,28) = b(k,28) - lu(k,92) * b(k,22)
         b(k,29) = b(k,29) - lu(k,93) * b(k,22)
         b(k,30) = b(k,30) - lu(k,94) * b(k,22)
         b(k,31) = b(k,31) - lu(k,95) * b(k,22)
         b(k,33) = b(k,33) - lu(k,96) * b(k,22)
         b(k,34) = b(k,34) - lu(k,97) * b(k,22)
         b(k,35) = b(k,35) - lu(k,98) * b(k,22)
                                                                        
         b(k,26) = b(k,26) - lu(k,101) * b(k,23)
         b(k,27) = b(k,27) - lu(k,102) * b(k,23)
         b(k,29) = b(k,29) - lu(k,103) * b(k,23)
         b(k,30) = b(k,30) - lu(k,104) * b(k,23)
         b(k,31) = b(k,31) - lu(k,105) * b(k,23)
         b(k,32) = b(k,32) - lu(k,106) * b(k,23)
         b(k,33) = b(k,33) - lu(k,107) * b(k,23)
                                                                        
         b(k,27) = b(k,27) - lu(k,110) * b(k,24)
         b(k,29) = b(k,29) - lu(k,111) * b(k,24)
         b(k,30) = b(k,30) - lu(k,112) * b(k,24)
         b(k,31) = b(k,31) - lu(k,113) * b(k,24)
         b(k,33) = b(k,33) - lu(k,114) * b(k,24)
         b(k,35) = b(k,35) - lu(k,115) * b(k,24)
                                                                        
         b(k,26) = b(k,26) - lu(k,120) * b(k,25)
         b(k,27) = b(k,27) - lu(k,121) * b(k,25)
         b(k,29) = b(k,29) - lu(k,122) * b(k,25)
         b(k,30) = b(k,30) - lu(k,123) * b(k,25)
         b(k,31) = b(k,31) - lu(k,124) * b(k,25)
         b(k,32) = b(k,32) - lu(k,125) * b(k,25)
         b(k,33) = b(k,33) - lu(k,126) * b(k,25)
         b(k,35) = b(k,35) - lu(k,127) * b(k,25)
                                                                        
         b(k,29) = b(k,29) - lu(k,130) * b(k,26)
         b(k,30) = b(k,30) - lu(k,131) * b(k,26)
         b(k,31) = b(k,31) - lu(k,132) * b(k,26)
         b(k,32) = b(k,32) - lu(k,133) * b(k,26)
         b(k,34) = b(k,34) - lu(k,134) * b(k,26)
         b(k,35) = b(k,35) - lu(k,135) * b(k,26)
                                                                        
         b(k,29) = b(k,29) - lu(k,138) * b(k,27)
         b(k,30) = b(k,30) - lu(k,139) * b(k,27)
         b(k,31) = b(k,31) - lu(k,140) * b(k,27)
         b(k,34) = b(k,34) - lu(k,141) * b(k,27)
                                                                        
         b(k,29) = b(k,29) - lu(k,147) * b(k,28)
         b(k,30) = b(k,30) - lu(k,148) * b(k,28)
         b(k,31) = b(k,31) - lu(k,149) * b(k,28)
         b(k,32) = b(k,32) - lu(k,150) * b(k,28)
         b(k,33) = b(k,33) - lu(k,151) * b(k,28)
         b(k,34) = b(k,34) - lu(k,152) * b(k,28)
         b(k,35) = b(k,35) - lu(k,153) * b(k,28)
                                                                        
         b(k,30) = b(k,30) - lu(k,179) * b(k,29)
         b(k,31) = b(k,31) - lu(k,180) * b(k,29)
         b(k,32) = b(k,32) - lu(k,181) * b(k,29)
         b(k,33) = b(k,33) - lu(k,182) * b(k,29)
         b(k,34) = b(k,34) - lu(k,183) * b(k,29)
         b(k,35) = b(k,35) - lu(k,184) * b(k,29)
                                                                        
         b(k,31) = b(k,31) - lu(k,194) * b(k,30)
         b(k,32) = b(k,32) - lu(k,195) * b(k,30)
         b(k,33) = b(k,33) - lu(k,196) * b(k,30)
         b(k,34) = b(k,34) - lu(k,197) * b(k,30)
         b(k,35) = b(k,35) - lu(k,198) * b(k,30)
                                                                        
         b(k,32) = b(k,32) - lu(k,218) * b(k,31)
         b(k,33) = b(k,33) - lu(k,219) * b(k,31)
         b(k,34) = b(k,34) - lu(k,220) * b(k,31)
         b(k,35) = b(k,35) - lu(k,221) * b(k,31)
                                                                        
         b(k,33) = b(k,33) - lu(k,233) * b(k,32)
         b(k,34) = b(k,34) - lu(k,234) * b(k,32)
         b(k,35) = b(k,35) - lu(k,235) * b(k,32)
                                                                        
         b(k,34) = b(k,34) - lu(k,249) * b(k,33)
         b(k,35) = b(k,35) - lu(k,250) * b(k,33)
                                                                        
         b(k,35) = b(k,35) - lu(k,266) * b(k,34)
                                                                        
      end do
                                                                        
      end subroutine IMP_LU_SLV01
                                                                        
      subroutine IMP_LU_SLV02( lu, b )
                                                                        
      use CHEM_MODS_MOD, only : imp_nzcnt, clsze, clscnt4
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real, intent(in)    ::   lu(clsze,imp_nzcnt)
      real, intent(inout) ::   b(clsze,clscnt4)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: k
                                                                        
!-----------------------------------------------------------------------
!       ... Solve L * y = b
!-----------------------------------------------------------------------
      do k = 1,clsze
                                                                        
!-----------------------------------------------------------------------
!       ... Solve U * x = y
!-----------------------------------------------------------------------
         b(k,35) = b(k,35) * lu(k,277)
         b(k,34) = b(k,34) - lu(k,276) * b(k,35)
         b(k,33) = b(k,33) - lu(k,275) * b(k,35)
         b(k,32) = b(k,32) - lu(k,274) * b(k,35)
         b(k,31) = b(k,31) - lu(k,273) * b(k,35)
         b(k,30) = b(k,30) - lu(k,272) * b(k,35)
         b(k,29) = b(k,29) - lu(k,271) * b(k,35)
         b(k,28) = b(k,28) - lu(k,270) * b(k,35)
         b(k,27) = b(k,27) - lu(k,269) * b(k,35)
         b(k,19) = b(k,19) - lu(k,268) * b(k,35)
         b(k,12) = b(k,12) - lu(k,267) * b(k,35)
                                                                        
         b(k,34) = b(k,34) * lu(k,265)
         b(k,33) = b(k,33) - lu(k,264) * b(k,34)
         b(k,32) = b(k,32) - lu(k,263) * b(k,34)
         b(k,31) = b(k,31) - lu(k,262) * b(k,34)
         b(k,30) = b(k,30) - lu(k,261) * b(k,34)
         b(k,29) = b(k,29) - lu(k,260) * b(k,34)
         b(k,28) = b(k,28) - lu(k,259) * b(k,34)
         b(k,27) = b(k,27) - lu(k,258) * b(k,34)
         b(k,26) = b(k,26) - lu(k,257) * b(k,34)
         b(k,22) = b(k,22) - lu(k,256) * b(k,34)
         b(k,17) = b(k,17) - lu(k,255) * b(k,34)
         b(k,11) = b(k,11) - lu(k,254) * b(k,34)
         b(k,9) = b(k,9) - lu(k,253) * b(k,34)
         b(k,8) = b(k,8) - lu(k,252) * b(k,34)
         b(k,5) = b(k,5) - lu(k,251) * b(k,34)
                                                                        
         b(k,33) = b(k,33) * lu(k,248)
         b(k,32) = b(k,32) - lu(k,247) * b(k,33)
         b(k,31) = b(k,31) - lu(k,246) * b(k,33)
         b(k,30) = b(k,30) - lu(k,245) * b(k,33)
         b(k,29) = b(k,29) - lu(k,244) * b(k,33)
         b(k,28) = b(k,28) - lu(k,243) * b(k,33)
         b(k,27) = b(k,27) - lu(k,242) * b(k,33)
         b(k,26) = b(k,26) - lu(k,241) * b(k,33)
         b(k,25) = b(k,25) - lu(k,240) * b(k,33)
         b(k,24) = b(k,24) - lu(k,239) * b(k,33)
         b(k,23) = b(k,23) - lu(k,238) * b(k,33)
         b(k,20) = b(k,20) - lu(k,237) * b(k,33)
         b(k,18) = b(k,18) - lu(k,236) * b(k,33)
                                                                        
         b(k,32) = b(k,32) * lu(k,232)
         b(k,31) = b(k,31) - lu(k,231) * b(k,32)
         b(k,30) = b(k,30) - lu(k,230) * b(k,32)
         b(k,29) = b(k,29) - lu(k,229) * b(k,32)
         b(k,27) = b(k,27) - lu(k,228) * b(k,32)
         b(k,26) = b(k,26) - lu(k,227) * b(k,32)
         b(k,25) = b(k,25) - lu(k,226) * b(k,32)
         b(k,24) = b(k,24) - lu(k,225) * b(k,32)
         b(k,23) = b(k,23) - lu(k,224) * b(k,32)
         b(k,20) = b(k,20) - lu(k,223) * b(k,32)
         b(k,16) = b(k,16) - lu(k,222) * b(k,32)
                                                                        
         b(k,31) = b(k,31) * lu(k,217)
         b(k,30) = b(k,30) - lu(k,216) * b(k,31)
         b(k,29) = b(k,29) - lu(k,215) * b(k,31)
         b(k,28) = b(k,28) - lu(k,214) * b(k,31)
         b(k,27) = b(k,27) - lu(k,213) * b(k,31)
         b(k,26) = b(k,26) - lu(k,212) * b(k,31)
         b(k,25) = b(k,25) - lu(k,211) * b(k,31)
         b(k,24) = b(k,24) - lu(k,210) * b(k,31)
         b(k,23) = b(k,23) - lu(k,209) * b(k,31)
         b(k,21) = b(k,21) - lu(k,208) * b(k,31)
         b(k,20) = b(k,20) - lu(k,207) * b(k,31)
         b(k,18) = b(k,18) - lu(k,206) * b(k,31)
         b(k,16) = b(k,16) - lu(k,205) * b(k,31)
         b(k,15) = b(k,15) - lu(k,204) * b(k,31)
         b(k,14) = b(k,14) - lu(k,203) * b(k,31)
         b(k,13) = b(k,13) - lu(k,202) * b(k,31)
         b(k,12) = b(k,12) - lu(k,201) * b(k,31)
         b(k,10) = b(k,10) - lu(k,200) * b(k,31)
         b(k,6) = b(k,6) - lu(k,199) * b(k,31)
                                                                        
         b(k,30) = b(k,30) * lu(k,193)
         b(k,29) = b(k,29) - lu(k,192) * b(k,30)
         b(k,28) = b(k,28) - lu(k,191) * b(k,30)
         b(k,27) = b(k,27) - lu(k,190) * b(k,30)
         b(k,21) = b(k,21) - lu(k,189) * b(k,30)
         b(k,19) = b(k,19) - lu(k,188) * b(k,30)
         b(k,17) = b(k,17) - lu(k,187) * b(k,30)
         b(k,10) = b(k,10) - lu(k,186) * b(k,30)
         b(k,8) = b(k,8) - lu(k,185) * b(k,30)
                                                                        
         b(k,29) = b(k,29) * lu(k,178)
         b(k,28) = b(k,28) - lu(k,177) * b(k,29)
         b(k,27) = b(k,27) - lu(k,176) * b(k,29)
         b(k,26) = b(k,26) - lu(k,175) * b(k,29)
         b(k,25) = b(k,25) - lu(k,174) * b(k,29)
         b(k,24) = b(k,24) - lu(k,173) * b(k,29)
         b(k,23) = b(k,23) - lu(k,172) * b(k,29)
         b(k,22) = b(k,22) - lu(k,171) * b(k,29)
         b(k,21) = b(k,21) - lu(k,170) * b(k,29)
         b(k,20) = b(k,20) - lu(k,169) * b(k,29)
         b(k,19) = b(k,19) - lu(k,168) * b(k,29)
         b(k,18) = b(k,18) - lu(k,167) * b(k,29)
         b(k,17) = b(k,17) - lu(k,166) * b(k,29)
         b(k,16) = b(k,16) - lu(k,165) * b(k,29)
         b(k,15) = b(k,15) - lu(k,164) * b(k,29)
         b(k,14) = b(k,14) - lu(k,163) * b(k,29)
         b(k,13) = b(k,13) - lu(k,162) * b(k,29)
         b(k,12) = b(k,12) - lu(k,161) * b(k,29)
         b(k,11) = b(k,11) - lu(k,160) * b(k,29)
         b(k,10) = b(k,10) - lu(k,159) * b(k,29)
         b(k,9) = b(k,9) - lu(k,158) * b(k,29)
         b(k,6) = b(k,6) - lu(k,157) * b(k,29)
         b(k,5) = b(k,5) - lu(k,156) * b(k,29)
         b(k,4) = b(k,4) - lu(k,155) * b(k,29)
         b(k,1) = b(k,1) - lu(k,154) * b(k,29)
                                                                        
         b(k,28) = b(k,28) * lu(k,146)
         b(k,27) = b(k,27) - lu(k,145) * b(k,28)
         b(k,22) = b(k,22) - lu(k,144) * b(k,28)
         b(k,21) = b(k,21) - lu(k,143) * b(k,28)
         b(k,7) = b(k,7) - lu(k,142) * b(k,28)
                                                                        
         b(k,27) = b(k,27) * lu(k,137)
         b(k,17) = b(k,17) - lu(k,136) * b(k,27)
                                                                        
         b(k,26) = b(k,26) * lu(k,129)
         b(k,17) = b(k,17) - lu(k,128) * b(k,26)
                                                                        
         b(k,25) = b(k,25) * lu(k,119)
         b(k,24) = b(k,24) - lu(k,118) * b(k,25)
         b(k,20) = b(k,20) - lu(k,117) * b(k,25)
         b(k,14) = b(k,14) - lu(k,116) * b(k,25)
                                                                        
         b(k,24) = b(k,24) * lu(k,109)
         b(k,15) = b(k,15) - lu(k,108) * b(k,24)
                                                                        
         b(k,23) = b(k,23) * lu(k,100)
         b(k,13) = b(k,13) - lu(k,99) * b(k,23)
                                                                        
         b(k,22) = b(k,22) * lu(k,91)
         b(k,18) = b(k,18) - lu(k,90) * b(k,22)
                                                                        
         b(k,21) = b(k,21) * lu(k,84)
                                                                        
         b(k,20) = b(k,20) * lu(k,79)
                                                                        
         b(k,19) = b(k,19) * lu(k,72)
                                                                        
         b(k,18) = b(k,18) * lu(k,67)
                                                                        
         b(k,17) = b(k,17) * lu(k,63)
                                                                        
         b(k,16) = b(k,16) * lu(k,58)
                                                                        
         b(k,15) = b(k,15) * lu(k,53)
                                                                        
         b(k,14) = b(k,14) * lu(k,48)
                                                                        
         b(k,13) = b(k,13) * lu(k,43)
                                                                        
         b(k,12) = b(k,12) * lu(k,38)
                                                                        
         b(k,11) = b(k,11) * lu(k,32)
                                                                        
         b(k,10) = b(k,10) * lu(k,28)
                                                                        
         b(k,9) = b(k,9) * lu(k,24)
         b(k,5) = b(k,5) - lu(k,23) * b(k,9)
                                                                        
         b(k,8) = b(k,8) * lu(k,19)
                                                                        
         b(k,7) = b(k,7) * lu(k,13)
                                                                        
         b(k,6) = b(k,6) * lu(k,10)
                                                                        
         b(k,5) = b(k,5) * lu(k,8)
         b(k,1) = b(k,1) - lu(k,7) * b(k,5)
                                                                        
         b(k,4) = b(k,4) * lu(k,5)
         b(k,3) = b(k,3) - lu(k,4) * b(k,4)
                                                                        
         b(k,3) = b(k,3) * lu(k,3)
                                                                        
         b(k,2) = b(k,2) * lu(k,2)
                                                                        
         b(k,1) = b(k,1) * lu(k,1)
                                                                        
      end do
                                                                        
      end subroutine IMP_LU_SLV02
                                                                        
      subroutine IMP_LU_SLV( lu, b )
                                                                        
      use CHEM_MODS_MOD, only : imp_nzcnt, clsze, clscnt4
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real, intent(in)    ::   lu(clsze,imp_nzcnt)
      real, intent(inout) ::   b(clsze,clscnt4)
                                                                        
      call IMP_LU_SLV01( lu, b )
      call IMP_LU_SLV02( lu, b )
                                                                        
      end subroutine IMP_LU_SLV
                                                                        
      end module MO_IMP_SOLVE_MOD

      module MO_ROD_SOLVE_MOD

      CONTAINS
                                                                        
      subroutine ROD_LU_SLV( lu, b )
                                                                        
      use CHEM_MODS_MOD, only : rod_nzcnt, clsze, clscnt5
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real, intent(in)    ::   lu(clsze,rod_nzcnt)
      real, intent(inout) ::   b(clsze,clscnt5)
                                                                        
                                                                        
      end subroutine ROD_LU_SLV
                                                                        
      end module MO_ROD_SOLVE_MOD
