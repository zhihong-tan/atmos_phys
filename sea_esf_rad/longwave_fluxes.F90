
                  module longwave_fluxes_mod


use rad_step_setup_mod,   only:  ISRAD, IERAD, JSRAD, JERAD, &
			         KS=>KSRAD, KE=>KERAD
use  utilities_mod,       only:  open_file, file_exist,       &
                                 check_nml_error, error_mesg, &
			         print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file
use rad_utilities_mod,    only:  Rad_control, radiation_control_type
use radiation_diag_mod,   only:  radiag_from_fluxes

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
   character(len=128)  :: version =  '$Id: longwave_fluxes.F90,v 1.3 2001/10/25 17:48:33 fms Exp $'
   character(len=128)  :: tag     =  '$Name: galway $'

!---------------------------------------------------------------------
!-------  interfaces --------

public    &
       longwave_fluxes_init, &
       longwave_fluxes_alloc,  &
       longwave_fluxes_ks, longwave_fluxes_k_down, &
       longwave_fluxes_KE_KEp1, longwave_fluxes_diag, &
       longwave_fluxes_dealloc, &
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


real, dimension(:,:,:,:), allocatable :: fluxn, fluxncf

integer   ::   flux_index4


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

subroutine longwave_fluxes_alloc (index4)

!---------------------------------------------------------------------
integer,                intent(in) :: index4

!---------------------------------------------------------------------

 allocate ( fluxn   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, index4) )
 fluxn  (:,:,:,:) = 0.0

 if (Rad_control%do_totcld_forcing) then
   allocate ( fluxncf (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, index4) )
   fluxncf(:,:,:,:) = 0.0
 endif

 flux_index4 = index4

!---------------------------------------------------------------------


end subroutine longwave_fluxes_alloc





!####################################################################

subroutine longwave_fluxes_ks (source, trans, iof, source2, trans2,  &
			       iof2, cld_trans, cld_ind, m)

!---------------------------------------------------------------------
integer,                   intent(in)   :: iof, iof2, m, &
	   		   	           cld_ind
real, dimension (:,:,:),   intent(in)   :: source
real, dimension (:,:,:),   intent(in)   :: source2
real, dimension (:,:,:),   intent(in)   :: trans2, trans
real, dimension (:,:,:,:), intent(in)   :: cld_trans
!---------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables
!------------------------------------------------------------------
        real, dimension(:,:,:), allocatable :: flux_tmp, flux_tmp2
        integer   ::   k

!---------------------------------------------------------------------
        allocate ( flux_tmp  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
        allocate ( flux_tmp2 (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )

	do k=KS+1, KE+1
          flux_tmp(:,:,k) = source(:,:,KS)*trans(:,:,k-1+iof)
          flux_tmp2(:,:,k) = source2(:,:,k)*trans2(:,:,k-1+iof2)
        end do

        if ((m  == 1) .or. (m  >= 7)) then
          fluxn(:,:,KS,m) = fluxn(:,:,KS,m) + source(:,:,KS)*  &
							   trans(:,:,KS)
        else
          fluxn(:,:,KS,m) = fluxn(:,:,KS,m) + source(:,:,KS)
        endif

        do k=KS+1,KE+1
          fluxn(:,:,k,m) = fluxn(:,:,k,m) + flux_tmp(:,:,k)* &
	                   cld_trans(:,:,k,cld_ind)
          fluxn(:,:,KS,m) = fluxn(:,:,KS,m) + flux_tmp2(:,:,k)*   &
		            cld_trans(:,:,k,cld_ind)
        end do

        if (Rad_control%do_totcld_forcing) then
          if ((m  == 1) .or. (m  >= 7)) then
            fluxncf(:,:,KS,m) =  source(:,:,KS)*trans(:,:,KS)
          else
            fluxncf(:,:,KS,m) =  source(:,:,KS)
          endif
          do k=KS+1,KE+1
	    fluxncf(:,:,k,m) = fluxncf(:,:,k,m) + flux_tmp(:,:,k)
	    fluxncf(:,:,KS,m) = fluxncf(:,:,KS,m) + flux_tmp2(:,:,k)
          end do
        endif

        deallocate (flux_tmp  )
        deallocate (flux_tmp2 )
!---------------------------------------------------------------------


end subroutine longwave_fluxes_ks



!####################################################################

subroutine longwave_fluxes_k_down (klevel, source, trans, trans2,   &
				   iof, cld_trans, m)

!---------------------------------------------------------------------
integer,                    intent(in)     :: iof, klevel, m
real,    dimension (:,:,:), intent(in)     :: cld_trans, source
real,    dimension (:,:,:), intent(in)     :: trans2, trans
!---------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables
!------------------------------------------------------------------
       real, dimension(:,:  ), allocatable  ::  flux4, flux4a
       real :: flux3a, flux_tmp, flux_tmp2

       integer     ::  kp, i, j, k, nn, ntot

!---------------------------------------------------------------------
       allocate ( flux4      (ISRAD:IERAD, JSRAD:JERAD         ) )
       allocate ( flux4a     (ISRAD:IERAD, JSRAD:JERAD         ) )


       do kp=klevel+1,KE+1
        do j=jsrad,jerad
          do i=israd,ierad   
         flux_tmp         = source(i,j,klevel)*trans(i,j,kp-1+iof)
         fluxn(i,j,kp,m) = fluxn(i,j,kp,m) +flux_tmp        *   &
                           cld_trans(i,j,kp)
       if (Rad_control%do_totcld_forcing) then
           fluxncf(i,j,kp,m) = fluxncf(i,j,kp,m) + flux_tmp         

       endif
       end do
       end do
         end do

          flux4(:,:       ) = 0.0
          flux4a(:,:       ) = 0.0
       do kp=klevel+1,KE+1
        do j=jsrad,jerad
          do i=israd,ierad   

         flux_tmp2         = source(i,j,kp)*trans2(i,j,kp-1+iof)
          flux4(i,j         ) = flux4(i,j         ) +    &
                               flux_tmp2        *cld_trans(i,j,kp)
       if (Rad_control%do_totcld_forcing) then
           flux4a (i,j         ) = flux4a (i,j         ) +  &
                                   flux_tmp2        
          endif
            end do
            end do
            end do
        do j=jsrad,jerad
          do i=israd,ierad   
           fluxn  (i,j,klevel,m) = fluxn  (i,j,klevel,m) +  &
                                   flux4    (i,j       )
       if (Rad_control%do_totcld_forcing) then
           fluxncf(i,j,klevel,m) = fluxncf(i,j,klevel,m) +  &
                                   flux4a   (i,j       )
       endif
       end do
       end do

       deallocate (flux4     )
       deallocate (flux4a    )
!---------------------------------------------------------------------


end subroutine longwave_fluxes_k_down



!####################################################################

subroutine longwave_fluxes_KE_KEp1 (source, trans, trans2, cld_trans, m)

!---------------------------------------------------------------------
integer,                      intent(in)    :: m
real,    dimension (:,:),     intent(in)    :: trans2, trans
real,    dimension (:,:,:),   intent(in)    :: source, cld_trans
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables
!--------------------------------------------------------------------
       real, dimension(:,:), allocatable :: flux_tmp, flux_tmp2

!-------------------------------------------------------------------
       allocate ( flux_tmp   (ISRAD:IERAD, JSRAD:JERAD)  )
       allocate ( flux_tmp2  (ISRAD:IERAD, JSRAD:JERAD)  )

       flux_tmp(:,:) = source(:,:,KE+1)*trans(:,:)
       flux_tmp2(:,:) = source(:,:,KE)*trans2(:,:)

       fluxn(:,:,KE,m) = fluxn(:,:,KE,m) +flux_tmp(:,:)*  &
                         cld_trans(:,:,KE+1)
       fluxn(:,:,KE+1,m) = fluxn(:,:,KE+1,m) +   &
                           flux_tmp2(:,:)*cld_trans(:,:,KE+1)

       if (Rad_control%do_totcld_forcing) then
         fluxncf(:,:,KE,m) = fluxncf(:,:,KE,m) + flux_tmp(:,:)
         fluxncf(:,:,KE+1,m) = fluxncf(:,:,KE+1,m) + flux_tmp2(:,:)
       endif

       deallocate (flux_tmp   )
       deallocate (flux_tmp2  )
!---------------------------------------------------------------------

end subroutine longwave_fluxes_KE_KEp1



!####################################################################

subroutine longwave_fluxes_diag (source, trans, cld_trans, m)

!---------------------------------------------------------------------
real, dimension (:,:,:),   intent(in)    :: cld_trans, source, trans
integer,                   intent(in)    :: m
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
        real, dimension(:,:,:), allocatable :: flux_tmp
        integer      ::   k

!---------------------------------------------------------------------
        allocate (flux_tmp (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
        do k=KS+1, KE+1
          flux_tmp(:,:,k) = source(:,:,k)*trans(:,:,k)
        end do

        do k=KS+1,KE+1
          fluxn(:,:,k,m) = fluxn(:,:,k,m) + flux_tmp(:,:,k)*   &
                           cld_trans(:,:,k)
        end do

        if (Rad_control%do_totcld_forcing) then
          do k=KS+1,KE+1
	    fluxncf(:,:,k,m) = fluxncf(:,:,k,m) +   &
                               flux_tmp(:,:,k)
          end do
        endif
        deallocate (flux_tmp  )
!---------------------------------------------------------------------


end subroutine longwave_fluxes_diag



!####################################################################

subroutine longwave_fluxes_dealloc

    deallocate (fluxn)
    if (Rad_control%do_totcld_forcing) then
      deallocate (fluxncf)
    endif


end subroutine longwave_fluxes_dealloc



!###################################################################

subroutine longwave_fluxes_sum (flux, NBTRGE, fluxcf)

!--------------------------------------------------------------------
integer,                          intent(in)    :: NBTRGE
real, dimension(:,:,:),           intent(out)   :: flux
real, dimension(:,:,:), optional, intent(out)   :: fluxcf
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
    integer       ::  m, j

!-------------------------------------------------------------------
    flux = 0.
    do m= 1, flux_index4            
      flux(:,:,:) = flux(:,:,:) + fluxn(:,:,:,m)
    end do

    if (Rad_control%do_totcld_forcing) then 
      fluxcf = 0.
      do m= 1, flux_index4            
	fluxcf(:,:,:) = fluxcf(:,:,:) + fluxncf(:,:,:,m)
      end do
    endif

    if (Rad_control%do_diagnostics) then
      do j=JSRAD, JERAD
        call radiag_from_fluxes (j, fluxn, fluxncf, NBTRGE)
      end do
    endif
!--------------------------------------------------------------------

end subroutine longwave_fluxes_sum


!#####################################################################


                end module longwave_fluxes_mod
