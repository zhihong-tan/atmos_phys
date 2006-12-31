module strat_chem_utilities_mod

use       mpp_io_mod, only : mpp_open, mpp_close, MPP_RDONLY
use          mpp_mod, only : mpp_pe, mpp_root_pe, stdout
use    constants_mod, only : DEG_TO_RAD
use time_manager_mod, only : time_type, get_date, days_in_month

implicit none
private

integer, parameter :: nlon_input=144, nlat_input=90, nlev_input=48, &
                      nspecies_age=8, nspecies_lbc=15, &
                      ntime_tropc=151, year_start_tropc=1950, nspecies_tropc=9
real, parameter :: agefact1 = 1.5, &
                   agefact2 = 1.25, &
                   clweight(7) = (/ 3., 2., 3., 4., 1., 3., 1. /)

real :: dfdage(nlat_input,nlev_input,nspecies_age), &
        tropc(ntime_tropc,nspecies_tropc)
real :: lat_input(nlat_input)
real, parameter :: tfact = 1./(365.25*86400.)
integer :: jstart

!-----------------------------------------------------------------------
!     ... interfaces
!-----------------------------------------------------------------------
public strat_chem_utilities_init, strat_chem_dcly_dt


!---- version number -----
character(len=128), parameter :: version     = ''
character(len=128), parameter :: tagname     = ''

logical :: module_is_initialized=.false.

CONTAINS

subroutine strat_chem_utilities_init(latb)

   implicit none
! dummy arguments
   real, intent(in) :: latb(:)
   
! local variables
   real :: age_dummy(nlat_input,nlev_input,nspecies_age), &
           chlb_dummy(nlat_input,nspecies_lbc), &
           ozb_dummy(nlon_input, nlat_input, 12)
   integer :: unit, nc, j
   
   if (module_is_initialized) return

!  read in chemical lower boundary 
   call mpp_open( unit, 'INPUT/chemlbf',action=MPP_RDONLY )
   if (mpp_pe() == mpp_root_pe()) WRITE(stdout(),*) 'INPUT/chemlbf'
   do nc = 1,15                                           
     read(unit,'(6E13.6)') chlb_dummy(:,nc)
   end do
   read(unit,'(6E13.6)') ozb_dummy
   read(unit,'(6e13.6)') tropc
   call mpp_close(unit)

!  read in data for Cly and Bry computation
   call mpp_open( unit, 'INPUT/ageair_fms_90.dat', action=MPP_RDONLY )
   if (mpp_pe() == mpp_root_pe()) &
      write(stdout(),*) 'strat_chem_utilities_init: Reading from INPUT/ageair_fms_90.dat'
   read(unit,'(6e13.6)') age_dummy
   read(unit,'(6e13.6)') dfdage
   call mpp_close(unit)

   jstart = 0   
   do j = 1,nlat_input
      lat_input(j) = ( -90. + (180./(nlat_input-1))*(j-1) ) * DEG_TO_RAD
      if (lat_input(j) >= latb(1) .and. lat_input(j) <= latb(2)) jstart = j
      if (mpp_pe() == mpp_root_pe()) &
         write(stdout(),*) 'strat_chem_utilities_init: jstart=',j,' on PE ',mpp_pe()
   end do
!
   module_is_initialized = .true.
   
end subroutine strat_chem_utilities_init

subroutine strat_chem_dcly_dt(Time, js, age, cly, bry, dclydt, dbrydt)

implicit none

type(time_type),        intent(in)  :: Time
integer,                intent(in)  :: js
real, dimension(:,:,:), intent(in)  :: age, cly, bry
real, dimension(:,:,:), intent(out) :: dclydt, dbrydt

! local variables

integer :: iyear, imon, iday, ihour, imin, isec
real :: time0, monthfrac, dt1, factor
integer :: it1,it2, imon2
integer :: ic, i, j, k, il, jl, kl
real :: clytot, brytot
real, dimension(size(age,1),size(age,2),nspecies_age) :: dfdtau
real, dimension(size(age,1),size(age,2),nspecies_tropc) :: cfc



call get_date( Time, iyear, imon, iday, ihour, imin, isec )

il = size(age,1)
jl = size(age,2)
kl = size(age,3)

!
!  Compute multiplying factor for missing CFCs, and include factor for 
!  conversion of rates to a per second rate. 
!
time0 = iyear + REAL(imon-1)/12.
it1 = INT(time0-year_start_tropc) + 1
it1 = min(max(it1,1),ntime_tropc-1)
it2 = it1+1
dt1 = time0 - (year_start_tropc-1) - it1

! sum1 = 0.
! sum2 = 0.
! do ic = 1,7
!    sum1 = sum1 + clweight(ic) * tropc(it1,ic)
!    sum2 = sum2 + clweight(ic) * tropc(it2,ic)
! end do
! factor = ((1-dt1)*tropc(it1,8) + dt1*tropc(it2,8))*tfact / (sum1*(1-dt1) + sum2*dt1)

monthfrac = REAL(iday)/REAL(days_in_month(Time))
imon2 = MOD(imon,12) + 1

level_loop: &
do k = 1,kl


! Copy dfdage to locally indexed variable
   do j = 1,jl
   do ic = 1,8
      dfdtau(:,j,ic) = dfdage(j+jstart-1+js-1,k,ic)
   end do
   end do

!
! Compute CFCs at time t - age
!

   do j = 1,jl
   do i = 1,il
      time0 = iyear + (imon+monthfrac-1)/12. - age(i,j,k)*agefact1
      it1 = INT(time0-year_start_tropc) + 1
      it1 = min(max(it1,1),ntime_tropc-1)
      it2 = it1 + 1
      dt1 = time0 - (year_start_tropc-1) - it1
      cfc(i,j,:) = tropc(it1,:)*(1-dt1) + tropc(it2,:)*dt1
      factor = cfc(i,j,8) / SUM(cfc(i,j,1:7)*clweight(1:7))
      dclydt(i,j,k) = 0.
      do ic = 1,7
         dclydt(i,j,k) = dclydt(i,j,k) &
                       + factor * 1.e-12 * tfact * agefact2 &
                       * dfdtau(i,j,ic) * clweight(ic) * cfc(i,j,ic)
      end do
      clytot = 1.e-12*cfc(i,j,8)
      if (cly(i,j,k) >= clytot) dclydt(i,j,k) = 0.
      dbrydt(i,j,k) = 1.e-12 * tfact * dfdtau(i,j,8)*cfc(i,j,9)
      brytot = 1.e-12*cfc(i,j,9)
      if (bry(i,j,k) >= brytot) dbrydt(i,j,k) = 0.
   end do
   end do
   
end do level_loop


end subroutine strat_chem_dcly_dt

end module strat_chem_utilities_mod
