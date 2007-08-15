      module MO_SETINV_MOD

implicit none
character(len=128), parameter :: version     = '$Id: mo_setinv.F90,v 14.0 2007/03/15 22:11:13 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: omsk $'
logical                       :: module_is_initialized = .false.

      CONTAINS

      subroutine setinv( invariants, tfld, h2ovmr, pmid, inv_data, plonl )
!-----------------------------------------------------------------
!        ... Set the invariant densities (molecules/cm**3)
!-----------------------------------------------------------------
      
      implicit none

!-----------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------
      integer, intent(in) ::    plonl
      real, intent(in)  ::      tfld(:,:)           ! temperature
      real, intent(in)  ::      h2ovmr(:,:)         ! water vapor vmr
      real, intent(in)  ::      pmid(:,:)           ! pressure
      real, intent(in)  ::      inv_data(:,:,:)     ! invariants
      real, intent(out) ::      invariants(:,:,:)   ! invariant array
!-----------------------------------------------------------------
!        .. Local variables
!-----------------------------------------------------------------
      real, parameter ::  boltz = 1.38044e-16      ! erg/K

      integer :: k,j
      integer :: plev
      
      plev = SIZE(tfld,2)

!-----------------------------------------------------------------
!        NOTE: Invariants are in cgs density units.
!              The pmid array is in pascals and must be
!	       mutiplied by 10. to yield dynes/cm**2.
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!	... Set M, N2, O2, and H2O densities
!-----------------------------------------------------------------
      do k = 1,plev
         invariants(:,k,1) = 10. * pmid(:,k) / (boltz*tfld(:,k))
         invariants(:,k,2) = .79 * invariants(:,k,1)
         invariants(:,k,3) = .21 * invariants(:,k,1)
!        invariants(:,k,4) = h2ovmr(:,k) * invariants(:,k,1)
!        do j = 5, size(invariants,3)
         do j = 4, size(invariants,3)
            invariants(:,k,j) = inv_data(:,k,j-3) * invariants(:,k,1)
         enddo
      end do

      end subroutine SETINV

      end module MO_SETINV_MOD
