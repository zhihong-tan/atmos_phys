        module aer_ccn_act_mod

use fms_mod,             only: error_mesg, FATAL, open_namelist_file, &
                               mpp_pe, mpp_root_pe, stdlog, &
                               file_exist, write_version_number, &
                               check_nml_error, close_file
use aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act2_k, &
                               aer_ccn_act_wpdf_k

implicit none
private
    private Loading, aer_ccn_act_init
      
    public aer_ccn_act, aer_ccn_act2, aer_ccn_act_wpdf ! cjg

!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id: aer_ccn_act.F90,v 15.0 2007/08/14 03:56:38 fms Exp $'
character(len=128) :: tagname = '$Name: omsk_2007_12 $'

!---------------- private data -------------------


!Parameters for look-up tables

  integer, parameter :: res = 20 !
  real, dimension(res,res,res,res) :: droplets

  integer, parameter :: res2 = 20 !
  real, dimension(res2,res2,res2,res2) :: droplets2


!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

logical  :: nooc = .false.   ! include organic aerosols as ccns ?

namelist /aer_ccn_act_nml/ nooc


logical :: module_is_initialized  = .false.
 
contains

subroutine aer_ccn_act (T1, P1, Updraft1, TotalMass, Drop)
real, dimension(:), intent(inout) :: TotalMass
real, intent(in) :: T1, P1, Updraft1
real, intent(inout) :: Drop
    
  integer :: tym, ier
  character(len=256) :: ermesg

  if(.not. module_is_initialized) call aer_ccn_act_init()

  tym = size (totalmass,1)

  call aer_ccn_act_k (T1, P1, Updraft1, TotalMass, tym, droplets,  &
                      droplets2, res, res2, nooc, Drop, ier, ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act', ermesg, FATAL)

  
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

real, dimension(:), intent(in) :: TotalMass
real, intent(in) :: T1, P1, Updraft1, mu,airdens, Nc, qc, qt, qe, tc, te
real, intent(inout) :: Drop

  integer :: tym, ier
  character(len=256) :: ermesg
        
  if(.not. module_is_initialized) call aer_ccn_act_init()
  tym = size (totalmass,1)

  call aer_ccn_act2_k (T1, P1, Updraft1, TotalMass, tym, mu,  &
                       airdens,Nc,qc,qt,qe,tc,te,Drop, ier, ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act2', ermesg, FATAL)

end subroutine aer_ccn_act2

!-->cjg: addition
!
! Additional subroutines to compute CCN activation by integrating
! over an assumed subgrid-scale PDF of w

subroutine aer_ccn_act_wpdf(T, p, wm, wp2, totalmass, drop)

! Compute CCN activation assuming a normal distribution of w
! given by its mean (wm) and second moment (wp2)

real :: T, p, wm, wp2, totalmass(3)
real :: drop

  integer :: tym, ier
  character(len=256) :: ermesg

  if(.not. module_is_initialized) call aer_ccn_act_init()
  tym = size (totalmass,1)

   call aer_ccn_act_wpdf_k(T, p, wm, wp2, totalmass, tym, droplets, &
            droplets2, res, res2, nooc, drop, ier, ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act_wpdf', ermesg, FATAL)

end subroutine aer_ccn_act_wpdf


subroutine aer_ccn_act_init ()

!--------------------------------------------------------------------  
!  local variables:
      
      integer   ::   unit, ierr, io       
      integer   ::   n
!--------------------------------------------------------------------- 
!    read namelist.
!--------------------------------------------------------------------
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=aer_ccn_act_nml, iostat=io,  &
               end=10)
        ierr = check_nml_error(io,'aer_ccn_act_nml')
        end do
10      call close_file (unit)   
      endif                      
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!--------------------------------------------------------------------
       call write_version_number (version, tagname)
       if (mpp_pe() == mpp_root_pe() ) &
                        write (stdlog(), nml=aer_ccn_act_nml)

       call Loading()

       module_is_initialized  = .true.

end subroutine aer_ccn_act_init


subroutine Loading()

real xx
integer i, j, k, l

  open(11, FILE='INPUT/droplets.dat')
  do k=1,res
    do i=1,res
      do j=1, res
        do l=1, res
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



end module aer_ccn_act_mod
