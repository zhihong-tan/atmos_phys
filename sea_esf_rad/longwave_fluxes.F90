
                  module longwave_fluxes_mod


!use rad_step_setup_mod,   only:  ISRAD, IERAD, JSRAD, JERAD, &
!			         KS=>KSRAD, KE=>KERAD
use  utilities_mod,       only:  open_file, file_exist,       &
                                 check_nml_error, error_mesg, &
			         print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file
use rad_utilities_mod,    only:  Rad_control, radiation_control_type, &
                                 lw_diagnostics_type
!use radiation_diag_mod,   only:  radiag_from_fluxes

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                  longwave fluxes module
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.08'
   character(len=128)  :: version =  '$Id: longwave_fluxes.F90,v 1.4 2002/07/16 22:35:39 fms Exp $'
   character(len=128)  :: tag     =  '$Name: inchon $'

!---------------------------------------------------------------------
!-------  interfaces --------

public    &
       longwave_fluxes_init, &
!      longwave_fluxes_alloc,  &
       longwave_fluxes_ks, longwave_fluxes_k_down, &
       longwave_fluxes_KE_KEp1, longwave_fluxes_diag, &
!      longwave_fluxes_dealloc, &
       longwave_fluxes_sum



!---------------------------------------------------------------------
!-------- namelist  ---------

real      ::  dummy = 1.0                    

namelist / longwave_fluxes_nml /        &
                                    dummy




!---------------------------------------------------------------------
!------- public data ------




!---------------------------------------------------------------------
!------- private data ------


!real, dimension(:,:,:,:), allocatable :: fluxn, fluxncf

!integer   ::   flux_index4

!        integer :: israd, ierad, jsrad, jerad, ks, ke

!---------------------------------------------------------------------
!---------------------------------------------------------------------



                         contains


subroutine longwave_fluxes_init 

!---------------------------------------------------------------------
     integer    ::  unit, ierr, io

!---------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=longwave_fluxes_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'longwave_fluxes_nml')
       enddo
10     call close_file (unit)
     endif

     unit = open_file ('logfile.out', action='append')
!    call print_version_number (unit, 'longwave_fluxes', version_number)
     if (get_my_pe() == 0) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=longwave_fluxes_nml)
     endif
     call close_file (unit)

!---------------------------------------------------------------------


end subroutine longwave_fluxes_init



!#####################################################################

!subroutine longwave_fluxes_alloc (ix, jx, kx, index4, Lw_diagnostics)

!---------------------------------------------------------------------
!integer,                intent(in) :: ix, jx, kx, index4
!type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics

!---------------------------------------------------------------------
!         integer :: israd, ierad, jsrad, jerad, ks, ke

!     israd = 1
!     ierad = ix
!     jsrad = 1
!     jerad = jx
!     ks = 1
!     ke = kx

! allocate ( Lw_diagnostics%fluxn   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, index4) )
! Lw_diagnostics%fluxn  (:,:,:,:) = 0.0
!
! if (Rad_control%do_totcld_forcing) then
!   allocate ( Lw_diagnostics%fluxncf (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, index4) )
!   Lw_diagnostics%fluxncf(:,:,:,:) = 0.0
! endif

! flux_index4 = index4

!---------------------------------------------------------------------


!end subroutine longwave_fluxes_alloc





!####################################################################

subroutine longwave_fluxes_ks ( source, trans, source2, trans2,  &
			        cld_trans, cld_ind, Lw_diagnostics)

!---------------------------------------------------------------------
integer, dimension(:),     intent(in)   ::  cld_ind
real, dimension (:,:,:,:),   intent(in)   :: source
real, dimension (:,:,:,:),   intent(in)   :: source2
real, dimension (:,:,:,:),   intent(in)   :: trans2, trans
real, dimension (:,:,:,:), intent(in)   :: cld_trans
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
!---------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables
!------------------------------------------------------------------
   real, dimension (size(source2,1), size(source2,2), &
          size(source2,3) ) :: flux_tmp, flux_tmp2
        integer   ::   k, ks, ke, nbands, m

!---------------------------------------------------------------------
        ks =1
	ke = size(source2,3)-1
	nbands = size(source,4)

   do m = 1, nbands

	do k=KS+1, KE+1
           flux_tmp(:,:,k) = source(:,:,KS,m)*trans(:,:,k,m    )
           flux_tmp2(:,:,k) = source2(:,:,k,m)*trans2(:,:,k ,m    )
        end do

        if ((m  == 1) .or. (m  >= 7)) then
          Lw_diagnostics%fluxn(:,:,KS,m) = Lw_diagnostics%fluxn(:,:,KS,m) + source(:,:,KS,m)*  &
				   trans(:,:,KS,m)
        else
          Lw_diagnostics%fluxn(:,:,KS,m) = Lw_diagnostics%fluxn(:,:,KS,m) + source(:,:,KS,m)
        endif

        do k=KS+1,KE+1
          Lw_diagnostics%fluxn(:,:,k,m) = Lw_diagnostics%fluxn(:,:,k,m) + flux_tmp(:,:,k)* &
	                   cld_trans(:,:,k,cld_ind(m))
          Lw_diagnostics%fluxn(:,:,KS,m) = Lw_diagnostics%fluxn(:,:,KS,m) + flux_tmp2(:,:,k)*   &
		            cld_trans(:,:,k,cld_ind(m))
        end do

        if (Rad_control%do_totcld_forcing) then
          if ((m  == 1) .or. (m  >= 7)) then
            Lw_diagnostics%fluxncf(:,:,KS,m) =  source(:,:,KS,m)*trans(:,:,KS,m)
          else
            Lw_diagnostics%fluxncf(:,:,KS,m) =  source(:,:,KS,m)
          endif
          do k=KS+1,KE+1
	    Lw_diagnostics%fluxncf(:,:,k,m) = Lw_diagnostics%fluxncf(:,:,k,m) + flux_tmp(:,:,k)
	    Lw_diagnostics%fluxncf(:,:,KS,m) = Lw_diagnostics%fluxncf(:,:,KS,m) + flux_tmp2(:,:,k)
          end do
        endif

     end do   ! (m loop)
!---------------------------------------------------------------------


end subroutine longwave_fluxes_ks



!####################################################################

subroutine longwave_fluxes_k_down (klevel, source, trans, trans2,   &
!				        cld_trans, m, Lw_diagnostics)
				        cld_trans, cld_ind,   Lw_diagnostics)

!---------------------------------------------------------------------
!integer,                    intent(in)     ::  klevel, m
integer,                    intent(in)     ::  klevel
!real,    dimension (:,:,:), intent(in)     :: cld_trans, source
real,    dimension (:,:,:,:), intent(in)     :: cld_trans, source
!real,    dimension (:,:,:), intent(in)     :: trans2, trans
real,    dimension (:,:,:,:), intent(in)     :: trans2, trans
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
integer, dimension(:), intent(in) :: cld_ind
!---------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables
!------------------------------------------------------------------
   real, dimension (size(source,1), size(source,2)) :: flux4, flux4a
       real :: flux3a, flux_tmp, flux_tmp2

       integer     ::  kp, i, j, k, nn, ntot, israd, ierad, jsrad, jerad
       integer :: ke
       integer :: m, nbands

       ierad = size(source,1)
       jerad = size(source,2)
       israd=1
       jsrad = 1
       ke = size(source,3)-1
       nbands = size(trans,4)
!---------------------------------------------------------------------

       do m=1, nbands

       do kp=klevel+1,KE+1
        do j=jsrad,jerad
          do i=israd,ierad   
         flux_tmp         = source(i,j,klevel,m)*trans(i,j,kp,m    )
         Lw_diagnostics%fluxn(i,j,kp,m) = Lw_diagnostics%fluxn(i,j,kp,m) +flux_tmp        *   &
                           cld_trans(i,j,kp, cld_ind(m))
       if (Rad_control%do_totcld_forcing) then
           Lw_diagnostics%fluxncf(i,j,kp,m) = Lw_diagnostics%fluxncf(i,j,kp,m) + flux_tmp         

       endif
       end do
       end do
         end do

          flux4(:,:       ) = 0.0
          flux4a(:,:       ) = 0.0
       do kp=klevel+1,KE+1
        do j=jsrad,jerad
          do i=israd,ierad   

         flux_tmp2         = source(i,j,kp,m)*trans2(i,j,kp,m      )
          flux4(i,j         ) = flux4(i,j         ) +    &
                        flux_tmp2        *cld_trans(i,j,kp, cld_ind(m))
       if (Rad_control%do_totcld_forcing) then
           flux4a (i,j         ) = flux4a (i,j         ) +  &
                                   flux_tmp2        
          endif
            end do
            end do
            end do
        do j=jsrad,jerad
          do i=israd,ierad   
           Lw_diagnostics%fluxn  (i,j,klevel,m) = Lw_diagnostics%fluxn  (i,j,klevel,m) +  &
                                   flux4    (i,j       )
       if (Rad_control%do_totcld_forcing) then
           Lw_diagnostics%fluxncf(i,j,klevel,m) = Lw_diagnostics%fluxncf(i,j,klevel,m) +  &
                                   flux4a   (i,j       )
       endif
       end do
       end do


   end do  ! (nbands loop)
!---------------------------------------------------------------------


end subroutine longwave_fluxes_k_down



!####################################################################

subroutine longwave_fluxes_KE_KEp1 (source, trans, trans2, cld_trans,&
!                                      m, Lw_diagnostics)
                                       cld_ind, Lw_diagnostics)

!---------------------------------------------------------------------
!integer,                      intent(in)    :: m
!real,    dimension (:,:),     intent(in)    :: trans2, trans
real,    dimension (:,:,:),     intent(in)    :: trans2, trans
!real,    dimension (:,:,:),   intent(in)    :: source, cld_trans
real,    dimension (:,:,:,:),   intent(in)    :: source, cld_trans
integer, dimension(:), intent(in) :: cld_ind
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables
!--------------------------------------------------------------------
   real, dimension (size(trans,1), size(trans,2)) :: flux_tmp, flux_tmp2
      integer :: ke
      integer :: m, nbands

!-------------------------------------------------------------------
        ke = size(source,3) - 1
	nbands = size(trans,3) 

do m=1,nbands
       flux_tmp(:,:) = source(:,:,KE+1,m)*trans(:,:,m)
       flux_tmp2(:,:) = source(:,:,KE,m)*trans2(:,:,m)

       Lw_diagnostics%fluxn(:,:,KE,m) = Lw_diagnostics%fluxn(:,:,KE,m) +flux_tmp(:,:)*  &
                         cld_trans(:,:,KE+1,cld_ind(m))
       Lw_diagnostics%fluxn(:,:,KE+1,m) = Lw_diagnostics%fluxn(:,:,KE+1,m) +   &
                           flux_tmp2(:,:)*cld_trans(:,:,KE+1,cld_ind(m))

       if (Rad_control%do_totcld_forcing) then
         Lw_diagnostics%fluxncf(:,:,KE,m) = Lw_diagnostics%fluxncf(:,:,KE,m) + flux_tmp(:,:)
         Lw_diagnostics%fluxncf(:,:,KE+1,m) = Lw_diagnostics%fluxncf(:,:,KE+1,m) + flux_tmp2(:,:)
       endif

end do  ! (nbands loop)
!---------------------------------------------------------------------

end subroutine longwave_fluxes_KE_KEp1



!####################################################################

!subroutine longwave_fluxes_diag (source, trans, cld_trans, m, Lw_diagnostics)
subroutine longwave_fluxes_diag (source, trans, cld_trans, cld_ind, &
                                 Lw_diagnostics)

!---------------------------------------------------------------------
!real, dimension (:,:,:),   intent(in)    :: cld_trans, source, trans
real, dimension (:,:,:,:),   intent(in)    :: cld_trans, source, trans
!integer,                   intent(in)    :: m
integer, dimension(:),     intent(in)    :: cld_ind
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
   real, dimension (size(trans,1), size(trans,2), &
            size(trans,3)) :: flux_tmp
        integer      ::   k, ks, ke
	integer :: m, nbands

!---------------------------------------------------------------------
         ks = 1
         ke = size(trans,3) - 1
	 nbands = size(trans,4)

do m=1,nbands


        do k=KS+1, KE+1
          flux_tmp(:,:,k) = source(:,:,k,m)*trans(:,:,k,m)
        end do

        do k=KS+1,KE+1
          Lw_diagnostics%fluxn(:,:,k,m) = Lw_diagnostics%fluxn(:,:,k,m) + flux_tmp(:,:,k)*   &
                           cld_trans(:,:,k,cld_ind(m))
        end do

        if (Rad_control%do_totcld_forcing) then
          do k=KS+1,KE+1
	    Lw_diagnostics%fluxncf(:,:,k,m) = Lw_diagnostics%fluxncf(:,:,k,m) +   &
                               flux_tmp(:,:,k)
          end do
        endif
!---------------------------------------------------------------------

end do ! (m loop)


end subroutine longwave_fluxes_diag



!####################################################################

!subroutine longwave_fluxes_dealloc

!    deallocate (fluxn)
!   if (Rad_control%do_totcld_forcing) then
!     deallocate (fluxncf)
!   endif


!end subroutine longwave_fluxes_dealloc



!###################################################################

subroutine longwave_fluxes_sum (is, ie, js, je, flux, NBTRGE,         &
                                 Lw_diagnostics, fluxcf)

!--------------------------------------------------------------------
integer,                          intent(in)    :: is, ie, js, &
                                                   je, NBTRGE
real, dimension(:,:,:),           intent(out)   :: flux
real, dimension(:,:,:), optional, intent(out)   :: fluxcf
type(lw_diagnostics_type), intent(in) :: Lw_diagnostics
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
    integer       ::  m, j

!-------------------------------------------------------------------
    flux = 0.
!   do m= 1, flux_index4            
    do m= 1, 6+NBTRGE               
      flux(:,:,:) = flux(:,:,:) + Lw_diagnostics%fluxn(:,:,:,m)
    end do

    if (Rad_control%do_totcld_forcing) then 
      fluxcf = 0.
!     do m= 1, flux_index4            
      do m= 1, 6+NBTRGE               
	fluxcf(:,:,:) = fluxcf(:,:,:) + Lw_diagnostics%fluxncf(:,:,:,m)
      end do
    endif

!     do j=js,je         
!      if (Rad_control%do_raddg(j)) then
!       call radiag_from_fluxes (Lw_diagnostics%fluxn,   &
!                 Lw_diagnostics%fluxncf, NBTRGE, j, j-js+1, &
!                                 is, ie)
!      endif
!     end do
!--------------------------------------------------------------------

end subroutine longwave_fluxes_sum


!#####################################################################


                end module longwave_fluxes_mod
