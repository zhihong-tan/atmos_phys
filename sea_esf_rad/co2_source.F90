
		 module co2_source_mod


use rad_step_setup_mod,   only:  KS=>KSRAD, KE=>KERAD, ISRAD, IERAD,&
			         JSRAD, JERAD
use longwave_params_mod,  only:  NBCO215
use rad_utilities_mod,    only:  looktab, longwave_tables3_type, &
				 Environment, environment_type
use  utilities_mod,       only:  open_file, file_exist,    &
                                 check_nml_error, error_mesg, &
			         print_version_number, FATAL, NOTE, &
				 WARNING, close_file, get_my_pe, &
				 read_data, write_data
use constants_new_mod,    only : secday, radcon
use longwave_setup_mod,   only : pdfinv, dte1, ixoe1, &
			         longwave_parameter_type, Lw_parameters
use radiative_gases_mod,  only : rrvco2      
use std_pressures_mod,    only : get_std_pressures
use longwave_tables_mod,  only : get_bds, tabsr
use gas_tf_mod,           only : transcol

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                 co2 longwave source module
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!     character(len=5), parameter  ::  version_number = 'v0.09'
      character(len=128)  :: version =  '$Id: co2_source.F90,v 1.2 2001/07/05 17:28:53 fms Exp $'
      character(len=128)  :: tag     =  '$Name: eugene $'

!---------------------------------------------------------------------
!------    interfaces   ------

public    co2_source_init, co2_source_calc,    &
	  co2_sorc_dealloc, get_sorc

private   nlte, co2curt


!---------------------------------------------------------------------
!------     namelist  -----

logical                 :: do_nlte = .false.


namelist / co2_source_nml /  &
				  do_nlte               



!---------------------------------------------------------------------
!---- public data -------



!---------------------------------------------------------------------
!---- private data -------


real, dimension(:,:,:,:), allocatable       :: sorc
real, dimension (:),      allocatable       :: c1b7, c2b7


integer            :: ixprnlte


!---------------------------------------------------------------------
!---------------------------------------------------------------------




                         contains




subroutine co2_source_init (kmin, kmax)

integer, intent(in)     :: kmin, kmax

!--------------------------------------------------------------------
   real, dimension(:), allocatable :: plm, bdlocm, bdhicm, cent, del

   integer         :: unit, ierr, io, k, n
   integer         :: ioffset
   real            :: prnlte
 
!---------------------------------------------------------------------
!-----  read namelist  ------
 
   if (file_exist('input.nml')) then
     unit =  open_file ('input.nml', action='read')
     ierr=1; do while (ierr /= 0)
     read (unit, nml=co2_source_nml, iostat=io, end=10)
     ierr = check_nml_error (io, 'co2_source_nml')
     enddo
10   call close_file (unit)
   endif

   unit = open_file ('logfile.out', action='append')
!  call print_version_number(unit, 'co2_source', version_number)
   if (get_my_pe() == 0) then
     write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
     write (unit,nml=co2_source_nml)
   endif
   call close_file (unit)

   if (do_nlte) then
!---------------------------------------------------------------------
!    define pressure-dependent index values used in the infrared
!    radiation code. by this manner, the coefficients are defined
!    at execution time (not dependent on the choice of vertical
!    layers)
!     prnlte : pressure (mb) below which non-LTE code (Nlte.F) affects
!              CO2 transmissivities
!---------------------------------------------------------------------
     prnlte = 0.1

!--------------------------------------------------------------------
!   convert pressure specification for bottom (flux) pressure level
!   for nlte calculation into an index (ixprnlte)
!-------------------------------------------------------------------
     allocate ( plm(kmin:kmax) )
     call get_std_pressures (plm_out=plm)

     ixprnlte = KS 
     do k=KS+1, KE
       if ((plm(k) - prnlte) .LT. 0.0) then
         ixprnlte = k 
       else
         exit
       endif
     enddo

     deallocate (plm)
   endif

!---------------------------------------------------------------------
!   allocate and obtain elements of the source function for bands in 
!   the 15 um range (used in nlte)
!---------------------------------------------------------------------
   if (do_nlte) then
     ioffset = Lw_parameters%offset
     allocate ( bdlocm  (9+ioffset:8+NBCO215+ioffset) )
     allocate ( bdhicm  (9+ioffset:8+NBCO215+ioffset) )
     allocate ( c1b7    (NBCO215) )
     allocate ( c2b7    (NBCO215) )
     allocate ( cent    (NBCO215) )
     allocate ( del     (NBCO215) )
     call get_bds (bdlocm, bdhicm)
     do n=1,NBCO215 
       cent(n) = 0.5E+00*(bdlocm(n+8+ioffset) + bdhicm(n+8+ioffset))
       del (n) = bdhicm(n+8+ioffset) - bdlocm(n+8+ioffset)
       c1b7(n) = (3.7412E-05)*cent(n)*cent(n)*cent(n)*del(n) 
       c2b7(n) = (1.4387E+00)*cent(n)
     end do
     deallocate (bdlocm)
     deallocate (bdhicm)
     deallocate (del  )
     deallocate (cent )
   endif

!--------------------------------------------------------------------



end subroutine co2_source_init



!####################################################################

subroutine co2_source_calc (press, soe1, soe2, soe3, soe4, soe5)

!--------------------------------------------------------------------
real, dimension(:,:,:), intent(in)  ::  press         
real, dimension(:,:,:), intent(out) ::  soe1, soe2, soe3, soe4, soe5
!----------------------------------------------------------------------

      integer            ::   n, ioffset
      integer            ::   NBLY                    

!--------------------------------------------------------------------
      ioffset = Lw_parameters%offset
      NBLY = 16+ioffset

      allocate (sorc(ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, 9+ioffset:NBLY))
       
!----------------------------------------------------------------------
!     compute source function for frequency bands (9+ioffset to NBLY-1) 
!     at layer temperatures using table lookup.
!----------------------------------------------------------------------
      do n=9+ioffset,NBLY-1
        call looktab (tabsr, ixoe1, dte1,   & 
                      sorc(:,:,:,n), KS, KE+1, n)
      enddo

!-----------------------------------------------------------------------
!     compute the nlte source function for co2.
!-----------------------------------------------------------------------
      if (do_nlte) then
        call nlte (press)
      endif

!--------------------------------------------------------------------
!     pass the values needed in longwave_driver back there. sorc will 
!     also be used in cool_to_space_exact.
!--------------------------------------------------------------------
      soe1(:,:,:)  = sorc(:,:,:, 9+ioffset ) + &
                     sorc(:,:,:,10+ioffset ) + &
                     sorc(:,:,:,11+ioffset ) 
      soe2(:,:,:)  = sorc(:,:,:,14+ioffset ) 
      soe3(:,:,:)  = sorc(:,:,:,12+ioffset ) 
      soe4(:,:,:)  = sorc(:,:,:,13+ioffset ) 
      soe5(:,:,:)  = sorc(:,:,:,15+ioffset ) 

!-------------------------------------------------------------------


end subroutine co2_source_calc




!#####################################################################

subroutine co2_sorc_dealloc

     deallocate (sorc)

end subroutine co2_sorc_dealloc



!#####################################################################

subroutine get_sorc (sorc_tmp_out, n)

!---------------------------------------------------------------------
real, dimension(:,:,:), intent(out)   :: sorc_tmp_out
integer,                intent(in)    :: n
!--------------------------------------------------------------------

     sorc_tmp_out(:,:,:) = sorc(:,:,:,n)
 
!---------------------------------------------------------------

end subroutine get_sorc



!#####################################################################

subroutine nlte (press)

!-----------------------------------------------------------------------
!     nlte is the present formulation of an nlte calculation of the 
!     source function in the 15 um region (two bands).
!
!     the essential theory is:
!
!           phi = C*j
!             j = b + E*phi
!
!     where
!             C = Curtis matrix
!	      E = NLTE contribution (diagonal matrix)
!           phi = heating rate vector
!             b = LTE source function vector
!             j = NLTE source function vector
!
!             j = b (by assumption) for pressure layers > ixnltr
!             j = b (by assumption) for pressure layers > ixprnlte
!      E is obtained using a formulation devised by Fels (denoted
!      Ri in his notes).
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!-----------------------------------------------------------------------
 real, dimension (:,:,:), intent(in)  ::  press               

!-----------------------------------------------------------------------
!     intent local:
!
!     degeneracy factor = 0.5
!
!                fnlte  = NLTE contribution: (E in above notes)
!
!                phifx  = fixed portion of PHI (contributions from
!                         layers > ixnltr, where j(k) = b(k))
!                         layers > ixprnlte, where j(k) = b(k))
!
!                phivar = varying portion of PHI (contributions
!                         from layers <= ixprnlte).
!                         from layers <= ixnltr).
!-----------------------------------------------------------------------
      real, dimension(:,:,:), allocatable    :: ag, az, bdenom, cdiag, &
	                                        tcoll, phifx, phivar
      real, dimension (:,:,:,:), allocatable :: cmtrx, fnlte
      real                                   :: degen = 0.5
      integer                                :: n, k, inb, kp, ioffset

!---------------------------------------------------------------------
      ioffset =  Lw_parameters%offset

!--------------------------------------------------------------------
      allocate ( ag      (ISRAD:IERAD, JSRAD:JERAD,     KS:ixprnlte)) 
      allocate ( az      (ISRAD:IERAD, JSRAD:JERAD,     KS:ixprnlte)) 
      allocate ( bdenom  (ISRAD:IERAD, JSRAD:JERAD,     KS:ixprnlte)) 
      allocate ( cdiag   (ISRAD:IERAD, JSRAD:JERAD,     KS:ixprnlte)) 
      allocate ( tcoll   (ISRAD:IERAD, JSRAD:JERAD,     KS:ixprnlte)) 
      allocate ( phifx   (ISRAD:IERAD, JSRAD:JERAD,     KS:ixprnlte)) 
      allocate ( phivar  (ISRAD:IERAD, JSRAD:JERAD,     KS:ixprnlte)) 
      allocate ( fnlte   (ISRAD:IERAD, JSRAD:JERAD,     KS:ixprnlte  , &
						            NBCO215  )) 
      allocate ( cmtrx   (ISRAD:IERAD, JSRAD:JERAD, KS:KE, KS:ixprnlte))

!-----------------------------------------------------------------------
!     compute curtis matrix for both frequency bands.
!-----------------------------------------------------------------------
      call co2curt (cmtrx)

      do k=KS,ixprnlte
        cdiag(:,:,k) = cmtrx(:,:,k,k)
      end do

!-----------------------------------------------------------------------
!   collisional relaxation time (see fels notes for "tcoll")
!-----------------------------------------------------------------------
      do k=KS,ixprnlte
        tcoll(:,:,k) = degen*1.5E-05*press(:,:,KE+1)/   &
                       (secday*press(:,:,k)) 
      end do

!-----------------------------------------------------------------------
!   compute NLTE contribution for eack band at each pressure level
!   <= ixprnlte. fnlte = zero by assumption at other levels.
!-----------------------------------------------------------------------
      do n=1,NBCO215
        fnlte (:,:,KS:ixprnlte,n) = 3.5E+00*tcoll(:,:,KS:ixprnlte)*  &
				    c1b7(n)/(rrvco2*c2b7(n)) 
      enddo

!-----------------------------------------------------------------------
!     begin computations for (NBCO215) bands in 15um range.
!-----------------------------------------------------------------------
      do inb = 1,NBCO215
        bdenom(:,:,KS:ixprnlte) = 1.0E+00/   &
              (1.0E+00 - fnlte(:,:,KS:ixprnlte,inb)*   &
			cdiag(:,:,KS:ixprnlte))
        phifx(:,:,KS:ixprnlte) = 0.0E+00
        do k=KS,ixprnlte
          do kp=ixprnlte+1,KE
            phifx(:,:,k) = phifx(:,:,k) +   &
                           cmtrx(:,:,kp,k)*sorc(:,:,kp,inb+8+ioffset )
          end do
        end do
        az(:,:,KS:ixprnlte) = sorc (:,:,KS:ixprnlte,inb+8+ioffset ) +  &
                     fnlte(:,:,KS:ixprnlte,inb)*phifx(:,:,KS:ixprnlte)

!----------------------------------------------------------------------
!     first iteration. (J(k) = B(k)) as initial guess)
!-----------------------------------------------------------------------
        phivar(:,:,KS:ixprnlte) = 0.0E+00
        do k=KS,ixprnlte
          do kp=KS,ixprnlte
            phivar(:,:,k) = phivar(:,:,k) +   &
                            cmtrx(:,:,kp,k)*sorc(:,:,kp,inb+8+ioffset )
          end do
        end do
        ag  (:,:,KS:ixprnlte) = fnlte(:,:,KS:ixprnlte,inb)*   &
                                (phivar(:,:,KS:ixprnlte) -   &
                                 cdiag(:,:,KS:ixprnlte)*  &
				 sorc(:,:,KS:ixprnlte,inb+8+ioffset ))

        sorc(:,:,KS:ixprnlte,inb+8+ioffset ) = bdenom(:,:,KS:ixprnlte)*&
                                               (az(:,:,KS:ixprnlte) + &
						ag(:,:,KS:ixprnlte)) 

!-----------------------------------------------------------------------
!     second iteration.  (J(k) = result of first iteration as guess)
!-----------------------------------------------------------------------
        phivar(:,:,KS:ixprnlte) = 0.0E+00
        do k=KS,ixprnlte
          do kp=KS,ixprnlte
            phivar(:,:,k) = phivar(:,:,k) +    &
                            cmtrx(:,:,kp,k)*sorc(:,:,kp,inb+8+ioffset )
          end do
        end do
        ag  (:,:,KS:ixprnlte) = fnlte(:,:,KS:ixprnlte,inb)*   &
                        (phivar(:,:,KS:ixprnlte) -   &
            cdiag(:,:,KS:ixprnlte)*sorc(:,:,KS:ixprnlte,inb+8+ioffset ))

        sorc(:,:,KS:ixprnlte,inb+8+ioffset ) = bdenom(:,:,KS:ixprnlte)*&
                                               (az(:,:,KS:ixprnlte) +  &
						ag(:,:,KS:ixprnlte)) 
      enddo

!-----------------------------------------------------------------------
      deallocate ( cmtrx)
      deallocate ( fnlte) 
      deallocate ( phivar) 
      deallocate ( phifx) 
      deallocate ( tcoll) 
      deallocate ( cdiag) 
      deallocate ( bdenom)
      deallocate ( az   ) 
      deallocate ( ag   ) 

!---------------------------------------------------------------------


end subroutine nlte



!#####################################################################

subroutine co2curt (cmtrx)

!----------------------------------------------------------------------
!     co2curt computes Curtis matrix elements derived from co2
!     transmission functions.
!     functions.
!
!     author: m. d. schwarzkopf
!
!     revised: 8/18/94
!
!     certified:  radiation version 1.0
!
!----------------------------------------------------------------------
real, dimension(:,:, :,:), intent(out)  ::  cmtrx                  

!-----------------------------------------------------------------------
!     intent out:
!
!       co21c  = column of transmission functions.
!
!       co21r  = row of transmission functions. 
!
!       cmtrx  = cutris matrix.
!-----------------------------------------------------------------------
     real, dimension(:,:,:), allocatable   ::  co2row, co2rowp
     integer                               :: k, krow, kp

!---------------------------------------------------------------------
      allocate ( co2row ( ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1) )
      allocate ( co2rowp( ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1) )

!-----------------------------------------------------------------------
!     compute co2 transmission functions.
!-----------------------------------------------------------------------
      co2row(:,:,KS:KE+1)  = 1.0E+00
      co2rowp(:,:,KS:KE+1) = 1.0E+00

!-----------------------------------------------------------------------
!    compute curtis matrix for rows from KS to ixprnlte
!-----------------------------------------------------------------------
      do k = KS,ixprnlte
	krow = k

        call transcol ( KS, krow, KS, KE+1, co2row)        

        call transcol ( KS, krow+1, KS, KE+1, co2rowp)        

        do kp=KS,KE-1 
          cmtrx(:,:,kp,k) = radcon*pdfinv(:,:,k)*   &
                            (co2rowp(:,:,kp) - co2rowp(:,:,kp+1) -  &
                             co2row(:,:,kp) + co2row(:,:,kp+1)) 
        end do

        cmtrx(:,:,KE,k) = radcon*pdfinv(:,:,k)*   &
                          (co2rowp(:,:,KE) - co2row(:,:,KE)) 
      enddo

!-----------------------------------------------------------------
      deallocate (co2rowp)
      deallocate (co2row )

!--------------------------------------------------------------------


end subroutine co2curt




!#####################################################################


	      end module co2_source_mod

