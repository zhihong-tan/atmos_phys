
                    module longwave_aerosol_mod

 
use rad_step_setup_mod,    only: KSRAD, KERAD, ISRAD, IERAD, JSRAD, &
			         JERAD, jabs, deltaz
use  rad_utilities_mod,    only: longwave_control_type, Lw_control, &
				 Environment, environment_type
use  utilities_mod,        only: open_file, file_exist,    &
			         check_nml_error, error_mesg, &
			         print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file, &
				 get_domain_decomp
use constants_new_mod,     only: diffac

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                    longwave aerosol module
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: longwave_aerosol.F90,v 1.2 2001/07/05 17:31:50 fms Exp $'
   character(len=128)  :: tag     =  '$Name: fez $'


!---------------------------------------------------------------------
!----- interfaces  -----
           
public         longwave_aerosol_init, &
	       aertau, longwave_aerosol_dealloc, &
               get_totaerooptdep, get_totaerooptdep_15, &
	       get_aerooptdep_KE_15


!---------------------------------------------------------------------
!---- namelist   -----


!---------------------------------------------------------------------
!    strat_lwmodel_type : identifier giving specific lw aerosol 
!                         type (max length = 12)
!    strat_lwmodel_nintvls : number of band intervals in the aerosol
!                            data (not necessarily the same as the  
!                            number of band intervals in the radiation 
!                            code) 
!    ktop_aer  : top model layer where aerosol extinction is to
!                be included (for all latitudes).
!    kbot_aer  : bottom model layer where aerosol extinction is to
!                be included (for all latitudes).
!    jmax_aerfile : number of latitudes in the aerosol input file
!    kmax_aerfile : number of vertical layers in the aerosol input file
!---------------------------------------------------------------------


logical           :: do_lwaerosol    = .false.
character(len=12) :: strat_lwmodel_type  =  '            '
character(len=12) :: lwaerosol_form      =  '          '
integer           :: strat_lwmodel_nintvls = 1
integer           :: imax_aerfile = -1
integer           :: jmax_aerfile = -1
integer           :: kmax_aerfile = -1 
integer           :: ktop_aer = -1
integer           :: kbot_aer = -10
			  


namelist / longwave_aerosol_nml /    &
			          do_lwaerosol,  &
				  lwaerosol_form, &
                                  strat_lwmodel_type, &
                                  strat_lwmodel_nintvls, &
                                  imax_aerfile, &
                                  jmax_aerfile, &
                                  kmax_aerfile, &
                                  ktop_aer, &
                                  kbot_aer


!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------


real, dimension (:,:,:), allocatable       :: totaerooptdep_15, &
                                              strlwext, strlwalb 
real, dimension (:,:,:,:), allocatable     :: aeralb, aerext,   &
					      aerconc, totaerooptdep
real, dimension (:,:), allocatable         :: strlwext_15, strlwalb_15,&
                                              aerooptdep_KE_15

integer                      :: NLWAERB
integer, parameter           :: NLWAERMODELSSTR = 1
logical                      :: do_prdlwaerosol
integer                      :: x(4), y(4)

!---------------------------------------------------------------------
!---------------------------------------------------------------------




                           contains




subroutine longwave_aerosol_init (kmin, kmax)

integer, intent(in)    :: kmin, kmax

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
     real, dimension (:,:,:), allocatable :: aerextivlstrlw,   &
		                             aerssalbivlstrlw,    &
			  	             aerasymmivlstrlw

     integer         :: unit, ierr, io
     integer         :: iounit, inaer1, inaer2, k, j, lwb, jd, jdf

!---------------------------------------------------------------------
!-----  read namelist  ------

     if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=longwave_aerosol_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'longwave_aerosol_nml')
       enddo
10     call close_file (unit)
     endif

     unit = open_file ('logfile.out', action='append')
!    call print_version_number (unit, 'longwave_aerosol',   &
!						 version_number)
     if (get_my_pe() == 0) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=longwave_aerosol_nml)
     endif
     call close_file (unit)

     Lw_control%do_lwaerosol = do_lwaerosol

!-------------------------------------------------------------------
!   ****** the remainder of this subroutine need be executed only if 
!   ****** do_lwaerosol is true.
!-------------------------------------------------------------------
 

!-------------------------------------------------------------------
!   if longwave aerosols are activated, use lwaerosol_form to define 
!   whether or not predicted lwaerosols are to be used
!-------------------------------------------------------------------
     if (do_lwaerosol) then

!-------------------------------------------------------------------
!  lw aerosol amounts and properties are inputted from measurements
!  or are prescribed
!-------------------------------------------------------------------
       if (trim(lwaerosol_form) == 'prescribed' .or.    &
           trim(lwaerosol_form) == 'observed')       then
         do_prdlwaerosol = .false.

!-------------------------------------------------------------------
!  lw aerosol amounts and properties are predicted by the model
!  currently this option is not implemented.
!-------------------------------------------------------------------
       else if (trim(lwaerosol_form) == 'predicted') then
         call error_mesg ('longwave_aerosol_init',   &
		'predicted longwave aerosols not yet implemented.', &
		                                             FATAL)
       else
!-------------------------------------------------------------------
!  error condition
!-------------------------------------------------------------------
         call error_mesg( 'longwave_aerosol_init',  &
                     ' lwaerosol_form is not an acceptable value.', & 
							     FATAL)
       endif

!--------------------------------------------------------------------
!   retrieve needed variables from other modules.
!--------------------------------------------------------------------
       call get_domain_decomp (x, y)
       jd = y(2)
       jdf = y(4) -y(3) + 1

!--------------------------------------------------------------------
!   check to see if the model configuration matches the form of the
!   data which is available (n30 SKYHI only, until interpolation 
!   becomes available).
!--------------------------------------------------------------------
       if (Environment%running_gcm) THEN
         if (kmax /= 40) call error_mesg ('longwave_aerosol_init', &
           ' cannot currently activate lw aerosols with this kerad', & 
							  FATAL)
         if (jd  /= 60) call error_mesg ('longwave_aerosol_init', &
            ' cannot currently activate lw aerosols with this jd', & 
							    FATAL)
       endif

!---------------------------------------------------------------------
!   define the number of band intervals in the aerosol data. allocate
!   arrays that are needed.
!---------------------------------------------------------------------
       NLWAERB = strat_lwmodel_nintvls

       allocate ( aerextivlstrlw  (NLWAERB,kmax_aerfile, jmax_aerfile) )
       allocate ( aerssalbivlstrlw(NLWAERB,kmax_aerfile, jmax_aerfile) )
       allocate ( aerasymmivlstrlw(NLWAERB,kmax_aerfile, jmax_aerfile) )

       allocate ( strlwext    (1:jdf, ktop_aer:kbot_aer, NLWAERB ) )
       allocate ( strlwalb    (1:jdf, ktop_aer:kbot_aer, NLWAERB ) )
       allocate ( strlwext_15 (1:jdf, ktop_aer:kbot_aer          ) )
       allocate ( strlwalb_15 (1:jdf, ktop_aer:kbot_aer          ) )

       if (do_prdlwaerosol) then
         allocate ( aeralb (1:jdf, ktop_aer:kbot_aer, NLWAERB,  &
						      NLWAERMODELSSTR))
         allocate ( aerext (1:jdf, ktop_aer:kbot_aer, NLWAERB,    &
						      NLWAERMODELSSTR))
         allocate ( aerconc(1:jdf, ktop_aer:kbot_aer, NLWAERB,     &
						      NLWAERMODELSSTR))
       endif
 
!--------------------------------------------------------------------
!     read in the time-independent longwave aerosol extinction and
!     albedo coefficients for rama's code. these are given for the
!     frequency band structure appropriate to the present longwave
!     radiation code, so there is no need (as there is in the 
!     shortwave code) to do a frequency interpolation.
!--------------------------------------------------------------------
       inaer1 = open_file ('INPUT/lwaerosolextdata',     &
			   form = 'ieee32', action='read')
       inaer2 = open_file ('INPUT/lwaerosolssalbdata',    &
			   form = 'ieee32', action='read')
 
!--------------------------------------------------------------------
!    the lw aerosol data are for jmax_aerfile latitudes. data are 
!    stored as (freq.(NLWAERB), lat.(jmax_aerfile), alt.(kmax_aerfile)) 
!    arrays.
!--------------------------------------------------------------------
       do j = 1,jmax_aerfile
         do k = 1,kmax_aerfile
           read (inaer1)  (aerextivlstrlw  (lwb,k,j), lwb=1,NLWAERB)
	   read (inaer2)  (aerssalbivlstrlw(lwb,k,j), lwb=1,NLWAERB)
         end do
       end do
       call close_file (inaer1)
       call close_file (inaer2)

!--------------------------------------------------------------------
!    (if needed): interpolate the aerosol data into model latitudes and
!    vertical layers
!----note----interpolation code does not presently exist, thus
!    requiring the restriction that jd=60 and KMAX=40 which is
!    implemented above-----------------------------------------
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   if there is a need to integrate over the frequency bands given
!   for the aerosol data (numbering strat_model_nintvls(model)) to the
!   model frequency ranges (numbering NLWAERB) do so here. (the method
!   should resemble code in Initialswr95.F; another array then is used
!   in the do loop that follows).
!--------------------------------------------------------------------
 
!-------------------------------------------------------------------
!   define arrays to hold the input data stored in the preferred 
!   index order.
!--------------------------------------------------------------------
       if (Environment%running_gcm) THEN
         do lwb = 1, NLWAERB
           do k = ktop_aer, kbot_aer
             do j = 1,jdf
	       strlwext (j,k,lwb) =  aerextivlstrlw (lwb,k,y(3)+j-1)
    	       strlwalb (j,k,lwb) =  aerssalbivlstrlw (lwb,k,y(3)+j-1)
             end do
           end do
         end do

       else if (Environment%running_standalone) then
!--------------------------------------------------------------------
!   use data for j=40 (~30n) to fill up the jd rows (usually 1)
!   used in the single column. this latitude could change!
!--------------------------------------------------------------------
       	 do lwb = 1, NLWAERB
           do j = 1,jdf
             do k = ktop_aer, kbot_aer
	       strlwext (j,k,lwb) = aerextivlstrlw (lwb,k,40)
      	       strlwalb (j,k,lwb) = aerssalbivlstrlw (lwb,k,40)
             end do
           end do
         end do
       endif
 
!--------------------------------------------------------------------
!     compute extinction, albedo for 560-800 cm-1 band (corresponding
!     to aerosol bands 1-3. this is done without use of a weighting
!     function.
!--------------------------------------------------------------------
       do j = 1,jdf
	 do k = ktop_aer, kbot_aer
	   strlwext_15(j,k) =    &
      	                    (70.*strlwext(j,k,1) +   &
                             70.*strlwext(j,k,2) + &
                          100.*strlwext(j,k,3)) / 240.
	   strlwalb_15(j,k) =    &
                            (70.*strlwalb(j,k,1) + &
                             70.*strlwalb(j,k,2) + &
                          100.*strlwalb(j,k,3)) / 240.
	 enddo
       enddo
!--------------------------------------------------------------------
!   deallocate local arrays
!--------------------------------------------------------------------

       deallocate (aerextivlstrlw)
       deallocate (aerssalbivlstrlw)
       deallocate (aerasymmivlstrlw)

     endif


end subroutine longwave_aerosol_init




!###################################################################

subroutine aertau (aeralb, aerext, aerconc)
 
!---------------------------------------------------------------------
real, dimension(:,:,:,:),  intent(in), optional   :: aeralb, aerext, &
						     aerconc

!---------------------------------------------------------------------
      integer                                :: j, k, n, nm
      real, dimension (:,:,:,:), allocatable :: aerooptdep
      real, dimension (:,:,:  ), allocatable :: aerooptdep_15

!---------------------------------------------------------------------
!allocate needed arrays
!---------------------------------------------------------------------
      allocate (totaerooptdep (ISRAD:IERAD, JSRAD:JERAD,         &
                                           KSRAD:KERAD+1,  NLWAERB) )
      allocate (aerooptdep_KE_15 (ISRAD:IERAD, JSRAD:JERAD           ) )
      allocate (totaerooptdep_15 (ISRAD:IERAD, JSRAD:JERAD,   &
					               KSRAD:KERAD+1))
      allocate (aerooptdep (ISRAD:IERAD, JSRAD:JERAD,    &
						KSRAD:KERAD, NLWAERB) )
      allocate (aerooptdep_15 (ISRAD:IERAD, JSRAD:JERAD,  KSRAD:KERAD) )

!---------------------------------------------------------------------
      aerooptdep = 0.0
      aerooptdep_15 = 0.0

      do j=JSRAD,JERAD
        do n = 1,NLWAERB
          do k=ktop_aer,kbot_aer
            aerooptdep(:,j,k,n) = aerooptdep(:,j,k,n) + diffac*  &
                                  (1.0 - strlwalb(jabs(j),k,n))* &
                                  strlwext(jabs(j),k,n)*  &
		                  deltaz(:,j,k)
          enddo
        enddo
        do k=ktop_aer,kbot_aer
          aerooptdep_15(:,j,k) = aerooptdep_15(:,j,k) + diffac*  &
                                 (1.0 - strlwalb_15(jabs(j),k))* &
                                 strlwext_15(jabs(j),k)*   &
			         deltaz(:,j,k)
        enddo
      enddo
 
      do n = 1,NLWAERB
        totaerooptdep (:,:,KSRAD,n) = 0.0E+00
        do k = KSRAD+1,KERAD+1
          totaerooptdep(:,:,k,n) = totaerooptdep(:,:,k-1,n) +  &
                                   aerooptdep(:,:,k-1,n)
        enddo
      enddo
      totaerooptdep_15 (:,:,KSRAD) = 0.0E+00
      do k = KSRAD+1,KERAD+1
	totaerooptdep_15(:,:,k) = totaerooptdep_15(:,:,k-1) +    &
                                  aerooptdep_15(:,:,k-1)
      enddo

      aerooptdep_KE_15(:,:) = aerooptdep_15(:,:,KERAD)

!---------------------------------------------------------------------

      deallocate ( aerooptdep_15 )
      deallocate ( aerooptdep )
 

end  subroutine aertau





!####################################################################

subroutine longwave_aerosol_dealloc

      deallocate ( totaerooptdep_15 )
      deallocate ( aerooptdep_KE_15 )
      deallocate ( totaerooptdep    )

end subroutine longwave_aerosol_dealloc


!####################################################################

subroutine get_totaerooptdep (n, totaer)

!----------------------------------------------------------------------
integer,                    intent(in)    :: n
real, dimension(:,:,:),     intent(out)   :: totaer

!----------------------------------------------------------------------

      totaer(:,:,:) = totaerooptdep(:,:,:,n)

!-------------------------------------------------------------------


end subroutine get_totaerooptdep




!#####################################################################

subroutine get_totaerooptdep_15 (totaer)

!----------------------------------------------------------------------
real, dimension(:,:,:),     intent(out)   :: totaer

!----------------------------------------------------------------------

      totaer(:,:,:) = totaerooptdep_15(:,:,:)

!-------------------------------------------------------------------


end subroutine get_totaerooptdep_15




!#####################################################################

subroutine get_aerooptdep_KE_15 (totaer)

!----------------------------------------------------------------------
real, dimension(:,:),       intent(out)   :: totaer

!----------------------------------------------------------------------

      totaer(:,:) = aerooptdep_KE_15(:,:)

!-------------------------------------------------------------------


end subroutine get_aerooptdep_KE_15




!#####################################################################


	end module longwave_aerosol_mod

