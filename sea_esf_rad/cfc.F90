
                      module   cfc_mod


use utilities_mod,      only:  open_file, file_exist,    &
			       check_nml_error, error_mesg, &
			       print_version_number, FATAL, NOTE, &
			       get_my_pe, read_data, write_data, &
			       close_file
use rad_utilities_mod,  only:  radiative_gases_type, &
                               optical_path_type

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!            module which supplies cfc-dependent contributions
!                     to longwave radiative fluxes      
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

         character(len=128)  :: version =  '$Id: cfc.F90,v 1.4 2003/04/09 20:58:44 fms Exp $'
         character(len=128)  :: tag     =  '$Name: inchon $'



!---------------------------------------------------------------------
!-------  interfaces --------

public                                    &
	 cfc_init,                  &
	 cfc_exact, cfc_exact_part, cfc_indx8, cfc_indx8_part, &
	 cfc_overod, cfc_overod_part



!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=8)           :: dummy    = '     '



namelist /cfc_nml/    &
                       dummy                                



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

!real, allocatable, dimension (:,:,:)   ::  totf11, totf12,   &
!					   totf113, totf22

integer      ::  kx   !  kx = # of half-levels
real         ::  rf11air, rf12air,  rf113air, rf22air




!---------------------------------------------------------------------
!---------------------------------------------------------------------




                       contains


subroutine cfc_init  

!---------------------------------------------------------------------

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
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=cfc_nml)
      endif
      call close_file (unit)

!---------------------------------------------------------------------

end subroutine cfc_init


!####################################################################

subroutine cfc_optical_depth (density, Rad_gases, Optical)

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
type(radiative_gases_type), intent(in) :: Rad_gases
type(optical_path_type), intent(inout) :: Optical 

!----------------------------------------------------------------------
      integer          ::      k
      real             :: rrf11, rrf12, rrf113, rrf22

!----------------------------------------------------------------------
      allocate ( Optical%totf11 (size(density,1), size(density,2),    &
                         size(density,3) ) )
      allocate ( Optical%totf12 (size(density,1), size(density,2),    &
                         size(density,3) ) )
      allocate ( Optical%totf113(size(density,1), size(density,2),    &
                         size(density,3) ) )
      allocate ( Optical%totf22 (size(density,1), size(density,2),    &
                         size(density,3) ) )
 
      kx = size(density,3)

      rrf11 = Rad_gases%rrvf11*rf11air
      rrf12 = Rad_gases%rrvf12*rf12air
      rrf113 = Rad_gases%rrvf113*rf113air
      rrf22 = Rad_gases%rrvf22*rf22air

!----------------------------------------------------------------------
!   compute summed optical paths for f11,f12, f113 and f22  with the 
!   diffusivity factor of 2 (appropriate for weak-line absorption 
!   limit).
!----------------------------------------------------------------------
      Optical%totf11(:,:,1) = 0.0E+00
      Optical%totf12(:,:,1) = 0.0E+00
      Optical%totf113(:,:,1) = 0.0E+00
      Optical%totf22 (:,:,1) = 0.0E+00
      do k=2,kx           
        Optical%totf11(:,:,k) = Optical%totf11(:,:,k-1) + density(:,:,k-1)*rrf11*2.0E+00
        Optical%totf12(:,:,k) = Optical%totf12(:,:,k-1) + density(:,:,k-1)*rrf12*2.0E+00
        Optical%totf113(:,:,k) = Optical%totf113(:,:,k-1) + density(:,:,k-1)*rrf113* &  
                                                      	      2.0E+00
        Optical%totf22(:,:,k) = Optical%totf22(:,:,k-1) + density(:,:,k-1)*rrf22*2.0E+00

      enddo
       
!--------------------------------------------------------------------


end subroutine cfc_optical_depth




!####################################################################

subroutine cfc_exact (index, Optical, cfc_tf)

!----------------------------------------------------------------------
!     cfc_exact computes exact cool-to-space transmission function 
!     for cfc for the desired band (given by index). 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index
type(optical_path_type), intent(in) :: Optical
real, dimension (:,:,:), intent(out)   :: cfc_tf

!----------------------------------------------------------------------
     kx = size (Optical%totf11,3) 
      cfc_tf(:,:,:          ) = 1.0E+00 -    &
                      strf113(index)*Optical%totf113(:,:,2:kx           ) -   &
                      strf22 (index)*Optical%totf22 (:,:,2:kx           ) -   &
                      strf11 (index)*Optical%totf11 (:,:,2:kx           ) -   &
                      strf12 (index)*Optical%totf12 (:,:,2:kx           )    


end subroutine cfc_exact




!####################################################################

subroutine cfc_exact_part (index, Optical, cfc_tf, klevel)

!----------------------------------------------------------------------
!     cfc_exact computes exact cool-to-space transmission function 
!     at levels below klevel for cfc for the band given by index. 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index, klevel
type(optical_path_type), intent(in) :: Optical
real, dimension (:,:,:), intent(out)   :: cfc_tf

!----------------------------------------------------------------------
      integer     ::  k

!----------------------------------------------------------------------
     kx = size (Optical%totf11,3) 
      do k=klevel,kx-1 
        cfc_tf(:,:,k) = 1.0E+00 -    &
                         strf113(index)*    &
		          (Optical%totf113(:,:,k+1) -  &
			   Optical%totf113(:,:,klevel)) -   &
                         strf22 (index)*   &
		          (Optical%totf22(:,:,k+1) -   &
			   Optical%totf22(:,:,klevel)) -   &
                         strf11 (index)*                          &
		          (Optical%totf11(:,:,k+1) -   &
			    Optical%totf11(:,:,klevel)) -   &
                         strf12 (index)*      &   
		          (Optical%totf12(:,:,k+1) -   &
			   Optical%totf12(:,:,klevel)) 
      end do


end subroutine cfc_exact_part



!####################################################################

subroutine cfc_indx8 (index, Optical, tcfc8)

!----------------------------------------------------------------------
!     cfc_indx8 computes transmission function for cfc for the band 8. 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index
type(optical_path_type), intent(in) :: Optical
real, dimension (:,:,:), intent(out)   :: tcfc8

!----------------------------------------------------------------------
      tcfc8 (:,:, :           ) = 1.0E+00 -    &
                       strf113(index)*Optical%totf113(:,:, :           ) -   &
                       strf22 (index)*Optical%totf22 (:,:, :           ) 


end subroutine cfc_indx8



!####################################################################

subroutine cfc_indx8_part (index, Optical, tcfc8, klevel)

!----------------------------------------------------------------------
!     cfc_indx8_part computes transmission function for cfc for 
!     the band 8. 
!----------------------------------------------------------------------

integer,                 intent(in)    :: index, klevel
type(optical_path_type), intent(in) :: Optical
real, dimension (:,:,:), intent(out)   :: tcfc8

!----------------------------------------------------------------------
    integer     ::      k

!----------------------------------------------------------------------
     kx = size (Optical%totf11,3) 
      do k=klevel,kx-1 
        tcfc8 (:,:,k+1) = 1.0E+00 -  strf113(index)*    &
	                  (Optical%totf113(:,:,k+1) -  &
			   Optical%totf113(:,:,klevel)) -   &
                          strf22 (index)*  &
		  	  (Optical%totf22(:,:,k+1) -   &
			    Optical%totf22(:,:,klevel)) 
      end do


end subroutine cfc_indx8_part




!####################################################################

subroutine cfc_overod (Optical, cfc_tf)

!----------------------------------------------------------------------
!     cfc_overod computes transmission function for cfc that is used   
!     with overod variable.
!----------------------------------------------------------------------

type(optical_path_type), intent(in) :: Optical
real, dimension (:,:,:), intent(out)   :: cfc_tf

!----------------------------------------------------------------------
     kx = size (Optical%totf11,3) 
      cfc_tf(:,:, :         ) = 1.0E+00 -    &
                           sf11315*Optical%totf113(:,:,2:kx           ) -   &
                           sf2215 *Optical%totf22 (:,:,2:kx           ) -  &
		           sf1115*Optical%totf11  (:,:,2:kx           ) -  &
		           sf1215*Optical%totf12  (:,:,2:kx           )


end subroutine cfc_overod




!####################################################################

subroutine cfc_overod_part (Optical, cfc_tf, klevel)

!----------------------------------------------------------------------
!     cfc_overod_part computes transmission function for cfc that is 
!     used with overod variable from klevel down.
!----------------------------------------------------------------------

type(optical_path_type), intent(in) :: Optical
real, dimension (:,:,:), intent(out)   :: cfc_tf
integer,                 intent(in)    :: klevel

!----------------------------------------------------------------------
      integer     ::      k

!----------------------------------------------------------------------
     kx = size (Optical%totf11,3) 
      do k=klevel,kx-1 
        cfc_tf(:,:,k) = 1.0E+00 - sf11315*      &
	   	        (Optical%totf113(:,:,k+1) - Optical%totf113(:,:,klevel)) - &
                           sf2215 *  &
		        (Optical%totf22 (:,:,k+1) - Optical%totf22 (:,:,klevel)) -   & 
		           sf1115*   &
		        (Optical%totf11 (:,:,k+1) - Optical%totf11 (:,:,klevel)) -   & 
		           sf1215*  &
		        (Optical%totf12 (:,:,k+1) - Optical%totf12 (:,:,klevel))   
      end do



end subroutine cfc_overod_part



!####################################################################



		  end module cfc_mod
