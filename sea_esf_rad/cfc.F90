
                      module   cfc_mod


use rad_step_setup_mod, only:  KSRAD, KERAD, ISRAD, IERAD, JSRAD, JERAD
use utilities_mod,      only:  open_file, file_exist,    &
			       check_nml_error, error_mesg, &
			       print_version_number, FATAL, NOTE, &
			       get_my_pe, read_data, write_data, &
			       close_file
use constants_new_mod,  only:  wtmair

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                       cfc gases module
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!        character(len=5), parameter  ::  version_number = 'v0.09'
         character(len=128)  :: version =  '$Id: cfc.F90,v 1.2 2001/07/05 17:27:24 fms Exp $'
         character(len=128)  :: tag     =  '$Name: galway $'



!---------------------------------------------------------------------
!-------  interfaces --------

public                                    &
	 cfc_init, cfc_time_vary, cfc_optical_depth,   &
	 cfc_exact, cfc_exact_part, cfc_indx8, cfc_indx8_part, &
	 cfc_overod, cfc_overod_part, cfc_dealloc



!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=5)           :: f11_amt  = 'FIXED'
character(len=5)           :: f12_amt  = 'FIXED'
character(len=5)           :: f113_amt = 'FIXED'
character(len=5)           :: f22_amt  = 'FIXED'


namelist /cfc_nml/    &
                       f11_amt, f12_amt, f113_amt, f22_amt



!---------------------------------------------------------------------
!------- public data ------




!---------------------------------------------------------------------
!------- private data ------

!--------------------------------------------------------------------
!       NBLWCFC =  number of frequency bands with cfc band strengths
!                  included. The bands have the same frequency ranges
!                  as those used for h2o calculations
!--------------------------------------------------------------------
integer, parameter :: NBLWCFC = 8


!--------------------------------------------------------------------
!   data for averaged f11 band strength
!--------------------------------------------------------------------
real strf11(NBLWCFC) 

data  strf11 /       &
         0.000000E+00,  0.000000E+00,  0.527655E+02,  0.297523E+04,  &
         0.134488E+03,  0.247279E+03,  0.710717E+03,  0.000000E+00/

!--------------------------------------------------------------------
!   data for averaged f12 band strength
!--------------------------------------------------------------------
real strf12(NBLWCFC) 

data strf12 /       &
         0.552499E+01,  0.136436E+03,  0.243867E+02,  0.612532E+03, &
         0.252378E+04,  0.438226E+02,  0.274950E+04,  0.000000E+00/

!--------------------------------------------------------------------
!   data for averaged f113 band strength
!--------------------------------------------------------------------
real strf113(NBLWCFC)

data strf113 /     &
         0.627223E+01,  0.690936E+02,  0.506764E+02,  0.122039E+04,  &
         0.808762E+03,  0.742843E+03,  0.109485E+04,  0.194768E+03/

!--------------------------------------------------------------------
!   data for averaged f22 band strength
!--------------------------------------------------------------------
real strf22(NBLWCFC) 

data strf22 /    &
         0.301881E+02,  0.550826E+01,  0.397496E+03,  0.124802E+04,  &
         0.190285E+02,  0.460065E+02,  0.367359E+04,  0.508838E+03/

!--------------------------------------------------------------------
!   data for averaged f11 560-800 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf1115=0.219856E+02

!--------------------------------------------------------------------
!   data for averaged f12 560-800 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf1215=0.515665E+02

!--------------------------------------------------------------------
!   data for averaged f113 560-800 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf11315=0.430969E+02

!--------------------------------------------------------------------
!   data for averaged f22 560-800 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf2215=0.176035E+03

!--------------------------------------------------------------------
!   data for averaged f11 800-990, 1070-1200 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf11ct=0.125631E+04

!--------------------------------------------------------------------
!   data for averaged f12 800-990, 1070-1200 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf12ct=0.201821E+04

!--------------------------------------------------------------------
!   data for averaged f113 800-990, 1070-1200 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf113ct=0.105362E+04

!--------------------------------------------------------------------
!   data for averaged f22 800-990, 1070-1200 cm-1 band strength
!--------------------------------------------------------------------
real  :: sf22ct=0.188775E+04

!--------------------------------------------------------------------

real, allocatable, dimension (:,:,:)   ::  totf11, totf12,   &
					   totf113, totf22

logical      ::  do_f11_init = .true.
logical      ::  do_f12_init = .true.
logical      ::  do_f113_init = .true.
logical      ::  do_f22_init = .true.
real         ::  rf11air, rf12air,  rf113air, rf22air, &
		 rrf11, rrf12, rrf113, rrf22,  &
		 rf11, rf12, rf113, rf22



!---------------------------------------------------------------------
!---------------------------------------------------------------------




                       contains


subroutine cfc_init (data_source, do_f11, do_f12, do_f113, do_f22, &
		     rrvf11_in, rrvf12_in, rrvf113_in, rrvf22_in) 

!---------------------------------------------------------------------
character(len=*), intent(in)    ::  data_source
logical,          intent(in)    ::  do_f11, do_f12, do_f113, do_f22
real, intent(in),optional :: rrvf11_in, rrvf12_in, rrvf113_in, rrvf22_in

!---------------------------------------------------------------------
!   the values of the molecular weights of f11 and f12 are derived
!   from elemental atomic weights adopted by the International Union of 
!   Pure and Applied Chemistry in 1961. These values are also used in 
!   the US Standard Atmosphere, 1976.
!   some previous radiative calculations at gfdl have used the
!   values  137.5, 121.0 for the molecular weights of f11 and f12.
!---------------------------------------------------------------------
      real       ::  wtmf11  = 137.36855
      real       ::  wtmf12  = 120.91395
      real       ::  wtmf113 = 187.3765
      real       ::  wtmf22  =  86.46892

!---------------------------------------------------------------------
!  optional supplied initial trace gas volume mixing ratios in (no./no.)
!---------------------------------------------------------------------
      real       ::  rf11_icrccm   = 1.00000E-09
      real       ::  rf12_icrccm   = 1.00000E-09
      real       ::  rf113_icrccm  = 1.00000E-09
      real       ::  rf22_icrccm   = 1.00000E-09

      real       ::  rf11_ipcc_92  = 2.68000E-10
      real       ::  rf12_ipcc_92  = 5.03000E-10
      real       ::  rf113_ipcc_92 = 8.20000E-11
      real       ::  rf22_ipcc_92  = 1.05000E-10

      real       ::  rf11_ipcc_98  = 2.68960E-10
      real       ::  rf12_ipcc_98  = 5.31510E-10
      real       ::  rf113_ipcc_98 = 8.58100E-11
      real       ::  rf22_ipcc_98  = 1.26520E-10

      real       ::  rf11_ipcc_80  = 1.57500E-10
      real       ::  rf12_ipcc_80  = 2.72500E-10
      real       ::  rf113_ipcc_80 = 2.31400E-11
      real       ::  rf22_ipcc_80  = 6.20200E-11

!---------------------------------------------------------------------
     integer    ::  unit, ierr, io, inrad


!---------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=cfc_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'cfc_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'cfc', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=cfc_nml)
      endif
      call close_file (unit)

!--------------------------------------------------------------------
!   define cfc mixing ratio conversion factors.
!--------------------------------------------------------------------
      rf11air  = wtmf11/wtmair
      rf12air  = wtmf12/wtmair
      rf113air = wtmf113/wtmair
      rf22air  = wtmf22/wtmair

!--------------------------------------------------------------------
!   define initial cfc volume mixing ratio to be used.
!--------------------------------------------------------------------
      if (trim(data_source) == 'icrccm') then
        rf11   = rf11_icrccm
        rf12   = rf12_icrccm
        rf113  = rf113_icrccm
        rf22   = rf22_icrccm
      else if (trim(data_source) == 'ipcc_92') then
        rf11   = rf11_ipcc_92  
        rf12   = rf12_ipcc_92
        rf113  = rf113_ipcc_92
        rf22   = rf22_ipcc_92
      else if (trim(data_source) == 'ipcc_98') then
        rf11   = rf11_ipcc_98
        rf12   = rf12_ipcc_98
        rf113  = rf113_ipcc_98
        rf22   = rf22_ipcc_98
      else if (trim(data_source) == 'ipcc_80') then
        rf11   = rf11_ipcc_80
        rf12   = rf12_ipcc_80
        rf113  = rf113_ipcc_80
        rf22   = rf22_ipcc_80
      else if (trim(data_source) == 'input') then
	if (file_exist ('INPUT/id1cfc') ) then
          inrad = open_file ('INPUT/id1cfc', form= 'formatted', &
                             action= 'read')
          read (inrad, FMT = '(5e18.10)')  rf11
          read (inrad, FMT = '(5e18.10)')  rf12
          read (inrad, FMT = '(5e18.10)')  rf113
          read (inrad, FMT = '(5e18.10)')  rf22
          call close_file (inrad)
	else
	  call error_mesg ( 'cfc_init', &
                      'desired cfc input file is not present ', FATAL)
	endif
      else if (trim(data_source) == 'restart') then
	 rf11 = rrvf11_in
	 rf12 = rrvf12_in
	 rf113 = rrvf113_in
	 rf22 = rrvf22_in
      else
	call error_mesg ('cfc_init', &
           'no valid data source is specified for cfc input ', FATAL)
      endif
	
!---------------------------------------------------------------------
!   if an individual cfc gas is not to be used in this realization, 
!   set its initial volume mixing ratio to zero.
!---------------------------------------------------------------------
      if (.not. do_f11 )   rf11  = 0.0
      if (.not. do_f12 )   rf12  = 0.0
      if (.not. do_f113)   rf113 = 0.0
      if (.not. do_f22 )   rf22  = 0.0

!---------------------------------------------------------------------

end subroutine cfc_init




!####################################################################

subroutine cfc_time_vary ( rrvf11, rrvf12, rrvf113, rrvf22, &
			   do_f11, do_f12, do_f113, do_f22)

!---------------------------------------------------------------------
real, intent(inout)    ::  rrvf11, rrvf12, rrvf113, rrvf22
logical, intent(in)    ::  do_f11, do_f12, do_f113, do_f22


!---------------------------------------------------------------------
!  define current volume and mass mixing ratios of f11
!---------------------------------------------------------------------
    if (do_f11) then

!---------------------------------------------------------------------
!   the f11 volume mixing ratio is set to the initial value (rf11) and
!   the mass mixing ratio is defined on the first access of this 
!   routine. after first access, this routine does nothing.
!---------------------------------------------------------------------
      if (trim(f11_amt) == 'FIXED' ) then
	if (do_f11_init) then
          rrvf11 = rf11
          rrf11  = rrvf11*rf11air
          do_f11_init = .false.
        endif

!---------------------------------------------------------------------
!   NOTE: TIME VARIATION OF RADIATIVE GASES NOT CURRENTLY AVAILABLE
!---------------------------------------------------------------------
      else if (trim(f11_amt) == 'VARY') then
!---------------------------------------------------------------------
!   this is where the time-variation of f11 will be added
!   define rrvf11, the volume mixing ratio = function of (rf11,time ) 
!   and then convert it into rrf11, the mass mixing ratio.
!        rrvf11 =   ?????
!        rrf11  = rrvf11*rf11air
!---------------------------------------------------------------------
	call error_mesg ('cfc_time_vary',  &
	       'time-varying cfcs not yet implemented', FATAL)
      else 
	call error_mesg ('cfc_time_vary',  &
	       'f11_amt has unacceptable value', FATAL)
      endif

!--------------------------------------------------------------------
!   if f11 turned off, set mixing ratios to zero.
!--------------------------------------------------------------------
    else
      rrf11  = 0.0
      rrvf11 = 0.0
    endif

!---------------------------------------------------------------------
!  define current volume and mass mixing ratios of f12
!---------------------------------------------------------------------
    if (do_f12) then

!---------------------------------------------------------------------
!   the f12 volume mixing ratio is set to the initial value (rf12) and
!   the mass mixing ratio is defined on the first access of this 
!   routine. after first access, this routine does nothing.
!---------------------------------------------------------------------
      if (trim(f12_amt) == 'FIXED' ) then
	if (do_f12_init) then
          rrvf12 = rf12
          rrf12  = rrvf12*rf12air
          do_f12_init = .false.
        endif

!---------------------------------------------------------------------
!   NOTE: TIME VARIATION OF RADIATIVE GASES NOT CURRENTLY AVAILABLE
!---------------------------------------------------------------------
      else if (trim(f12_amt) == 'VARY') then
!---------------------------------------------------------------------
!   this is where the time-variation of f12 will be added
!   define rrvf12, the volume mixing ratio = function of (rf12,time ) 
!   and then convert it into rrf12, the mass mixing ratio.
!        rrvf12 =   ?????
!        rrf12  = rrvf12*rf12air
!---------------------------------------------------------------------
	call error_mesg ('cfc_time_vary',  &
	       'time-varying cfcs not yet implemented', FATAL)
      else 
	call error_mesg ('cfc_time_vary',  &
	       'f12_amt has unacceptable value', FATAL)
      endif

!--------------------------------------------------------------------
!   if f12 turned off, set mixing ratios to zero.
!--------------------------------------------------------------------
    else
      rrf12  = 0.0
      rrvf12 = 0.0
    endif

!---------------------------------------------------------------------
!  define current volume and mass mixing ratios of f113
!---------------------------------------------------------------------
    if (do_f113) then

!---------------------------------------------------------------------
!   the f113 volume mixing ratio is set to the initial value (rf113) and
!   the mass mixing ratio is defined on the first access of this 
!   routine. after first access, this routine does nothing.
!---------------------------------------------------------------------
      if (trim(f113_amt) == 'FIXED' ) then
	if (do_f113_init) then
          rrvf113 = rf113
          rrf113  = rrvf113*rf113air
          do_f113_init = .false.
        endif

!---------------------------------------------------------------------
!   NOTE: TIME VARIATION OF RADIATIVE GASES NOT CURRENTLY AVAILABLE
!---------------------------------------------------------------------
      else if (trim(f113_amt) == 'VARY') then
!---------------------------------------------------------------------
!   this is where the time-variation of f113 will be added
!   define rrvf113, the volume mixing ratio = function of (rf113,time ) 
!   and then convert it into rrf113, the mass mixing ratio.
!        rrvf113 =   ?????
!        rrf113  = rrvf113*rf113air
!---------------------------------------------------------------------
	call error_mesg ('cfc_time_vary',  &
	       'time-varying cfcs not yet implemented', FATAL)
      else 
	call error_mesg ('cfc_time_vary',  &
	       'f113_amt has unacceptable value', FATAL)
      endif

!--------------------------------------------------------------------
!   if f113 turned off, set mixing ratios to zero.
!--------------------------------------------------------------------
    else
      rrf113  = 0.0
      rrvf113 = 0.0
    endif

!---------------------------------------------------------------------
!  define current volume and mass mixing ratios of f22
!---------------------------------------------------------------------
    if (do_f22) then

!---------------------------------------------------------------------
!   the f22 volume mixing ratio is set to the initial value (rf22) and
!   the mass mixing ratio is defined on the first access of this 
!   routine. after first access, this routine does nothing.
!---------------------------------------------------------------------
      if (trim(f22_amt) == 'FIXED' ) then
	if (do_f22_init) then
          rrvf22 = rf22
          rrf22  = rrvf22*rf22air
          do_f22_init = .false.
        endif

!---------------------------------------------------------------------
!   NOTE: TIME VARIATION OF RADIATIVE GASES NOT CURRENTLY AVAILABLE
!---------------------------------------------------------------------
      else if (trim(f22_amt) == 'VARY') then
!---------------------------------------------------------------------
!   this is where the time-variation of f22 will be added
!   define rrvf22, the volume mixing ratio = function of (rf22,time ) 
!   and then convert it into rrf22, the mass mixing ratio.
!        rrvf22 =   ?????
!        rrf22  = rrvf22*rf22air
!---------------------------------------------------------------------
	call error_mesg ('cfc_time_vary',  &
	       'time-varying cfcs not yet implemented', FATAL)
      else 
	call error_mesg ('cfc_time_vary',  &
	       'f22_amt has unacceptable value', FATAL)
      endif

!--------------------------------------------------------------------
!   if f22 turned off, set mixing ratios to zero.
!--------------------------------------------------------------------
    else
      rrf22  = 0.0
      rrvf22 = 0.0
    endif
!---------------------------------------------------------------------


end subroutine cfc_time_vary



!####################################################################

subroutine cfc_optical_depth (density)

!----------------------------------------------------------------------
!     cfc_optical_depth computes optical paths for cfc. The code assumes
!     a constant mixing ratio throughout the atmosphere.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!----------------------------------------------------------------------

real, dimension (:,:,:), intent(in)    :: density 

!----------------------------------------------------------------------
      integer          ::      k

!----------------------------------------------------------------------
      allocate ( totf11 (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1     ) )
      allocate ( totf12 (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1     ) )
      allocate ( totf113(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1     ) )
      allocate ( totf22 (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1     ) )
 
!----------------------------------------------------------------------
!   compute summed optical paths for f11,f12, f113 and f22  with the 
!   diffusivity factor of 2 (appropriate for weak-line absorption 
!   limit).
!----------------------------------------------------------------------
      totf11(:,:,KSRAD) = 0.0E+00
      totf12(:,:,KSRAD) = 0.0E+00
      totf113(:,:,KSRAD) = 0.0E+00
      totf22 (:,:,KSRAD) = 0.0E+00
      do k=KSRAD+1,KERAD+1
        totf11(:,:,k) = totf11(:,:,k-1) + density(:,:,k-1)*rrf11*2.0E+00
	totf12(:,:,k) = totf12(:,:,k-1) + density(:,:,k-1)*rrf12*2.0E+00
        totf113(:,:,k) = totf113(:,:,k-1) + density(:,:,k-1)*rrf113* &  
                                                      	      2.0E+00
        totf22(:,:,k) = totf22(:,:,k-1) + density(:,:,k-1)*rrf22*2.0E+00
      enddo
       
!--------------------------------------------------------------------


end subroutine cfc_optical_depth




!####################################################################

subroutine cfc_exact (index, cfc_tf)

!----------------------------------------------------------------------
!     cfc_exact computes exact cool-to-space transmission function 
!     for cfc for the desired band (given by index). 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index
real, dimension (:,:,:), intent(out)   :: cfc_tf

!----------------------------------------------------------------------
      cfc_tf(:,:,KSRAD:KERAD) = 1.0E+00 -    &
                      strf113(index)*totf113(:,:,KSRAD+1:KERAD+1) -   &
                      strf22 (index)*totf22 (:,:,KSRAD+1:KERAD+1) -   &
                      strf11 (index)*totf11 (:,:,KSRAD+1:KERAD+1) -   &
                      strf12 (index)*totf12 (:,:,KSRAD+1:KERAD+1)    


end subroutine cfc_exact




!####################################################################

subroutine cfc_exact_part (index, cfc_tf, klevel)

!----------------------------------------------------------------------
!     cfc_exact computes exact cool-to-space transmission function 
!     at levels below klevel for cfc for the band given by index. 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index, klevel
real, dimension (:,:,:), intent(out)   :: cfc_tf

!----------------------------------------------------------------------
      integer     ::  k

!----------------------------------------------------------------------
      do k=klevel,KERAD
        cfc_tf(:,:,k) = 1.0E+00 -    &
                         strf113(index)*    &
		          (totf113(:,:,k+1) - totf113(:,:,klevel)) -   &
                         strf22 (index)*   &
		          (totf22(:,:,k+1) - totf22(:,:,klevel)) -   &
                         strf11 (index)*                          &
		          (totf11(:,:,k+1) - totf11(:,:,klevel)) -   &
                         strf12 (index)*      &   
		          (totf12(:,:,k+1) - totf12(:,:,klevel)) 
      end do


end subroutine cfc_exact_part



!####################################################################

subroutine cfc_indx8 (index, tcfc8)

!----------------------------------------------------------------------
!     cfc_indx8 computes transmission function for cfc for the band 8. 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index
real, dimension (:,:,:), intent(out)   :: tcfc8

!----------------------------------------------------------------------
      tcfc8 (:,:,KSRAD:KERAD+1) = 1.0E+00 -    &
                       strf113(index)*totf113(:,:,KSRAD:KERAD+1) -   &
                       strf22 (index)*totf22 (:,:,KSRAD:KERAD+1) 


end subroutine cfc_indx8



!####################################################################

subroutine cfc_indx8_part (index, tcfc8, klevel)

!----------------------------------------------------------------------
!     cfc_indx8_part computes transmission function for cfc for 
!     the band 8. 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index, klevel
real, dimension (:,:,:), intent(out)   :: tcfc8

!----------------------------------------------------------------------
    integer     ::      k

!----------------------------------------------------------------------
      do k=klevel,KERAD
        tcfc8 (:,:,k+1) = 1.0E+00 -  strf113(index)*    &
	                  (totf113(:,:,k+1) - totf113(:,:,klevel)) -   &
                          strf22 (index)*  &
		  	  (totf22(:,:,k+1) - totf22(:,:,klevel)) 
      end do


end subroutine cfc_indx8_part




!####################################################################

subroutine cfc_overod (cfc_tf)

!----------------------------------------------------------------------
!     cfc_overod computes transmission function for cfc that is used   
!     with overod variable.
!----------------------------------------------------------------------

real, dimension (:,:,:), intent(out)   :: cfc_tf

!----------------------------------------------------------------------
      cfc_tf(:,:,KSRAD:KERAD) = 1.0E+00 -    &
                           sf11315*totf113(:,:,KSRAD+1:KERAD+1) -   &
                           sf2215 *totf22 (:,:,KSRAD+1:KERAD+1) -  &
		           sf1115*totf11  (:,:,KSRAD+1:KERAD+1) -  &
		           sf1215*totf12  (:,:,KSRAD+1:KERAD+1)


end subroutine cfc_overod




!####################################################################

subroutine cfc_overod_part (cfc_tf, klevel)

!----------------------------------------------------------------------
!     cfc_overod_part computes transmission function for cfc that is 
!     used with overod variable from klevel down.
!----------------------------------------------------------------------

real, dimension (:,:,:), intent(out)   :: cfc_tf
integer,                 intent(in)    :: klevel

!----------------------------------------------------------------------
      integer     ::      k

!----------------------------------------------------------------------
      do k=klevel,KERAD
        cfc_tf(:,:,k) = 1.0E+00 - sf11315*      &
	   	        (totf113(:,:,k+1) - totf113(:,:,klevel)) - &
                           sf2215 *  &
		        (totf22 (:,:,k+1) - totf22 (:,:,klevel)) -   & 
		           sf1115*   &
		        (totf11 (:,:,k+1) - totf11 (:,:,klevel)) -   & 
		           sf1215*  &
		        (totf12 (:,:,k+1) - totf12 (:,:,klevel))   
      end do



end subroutine cfc_overod_part



!####################################################################

subroutine cfc_dealloc

    deallocate (totf22 )
    deallocate (totf113)
    deallocate (totf12 )
    deallocate (totf11 )

end subroutine cfc_dealloc



!####################################################################


		  end module cfc_mod
