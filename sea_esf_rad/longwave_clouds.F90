
                     module longwave_clouds_mod
 

use utilities_mod,         only:  open_file, file_exist,    &
                                  check_nml_error, error_mesg, &
				  print_version_number, FATAL, NOTE, &
				  WARNING, get_my_pe, close_file
use constants_mod,         only:  radcon
!use longwave_setup_mod,    only:  pdfinv,  Lw_parameters, &
!use longwave_setup_mod,    only:           Lw_parameters, &
!				  longwave_parameter_type
use rad_utilities_mod,     only:  Lw_control, longwave_control_type, &
				  lw_output_type, &
		  cloudrad_control_type, Cldrad_control, &
                                  lw_clouds_type, &
                                  Lw_parameters, &
				  longwave_parameter_type, &
                                  cldrad_properties_type


!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                  longwave cloud module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128)  :: version =  '$Id: longwave_clouds.F90,v 1.4 2003/04/09 20:59:50 fms Exp $'
  character(len=128)  :: tag     =  '$Name: inchon $'


!---------------------------------------------------------------------
!----- interfaces  -----
           
public       &
		longwave_clouds_init, &
		cldtau, &
		cloud, &
		thickcld

!---------------------------------------------------------------------
!---- namelist   -----

!logical    :: do_lwcldemiss = .false.
logical    :: dummy         = .false.


namelist / longwave_clouds_nml /   &
                                     dummy
!                                    do_lwcldemiss

!----------------------------------------------------------------------
!--- public data ---------




!----------------------------------------------------------------------
!---   private ---------

integer        :: NLWCLDB


!---------------------------------------------------------------------
!---------------------------------------------------------------------





      contains





subroutine longwave_clouds_init

!--------------------------------------------------------------------
        integer               :: unit, ierr, io

!--------------------------------------------------------------------
!-----  read namelist  ------
  
    if (file_exist('input.nml')) then
       unit =  open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
       read (unit, nml=longwave_clouds_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'longwave_clouds_nml')
       enddo
10     call close_file (unit)
    endif

    unit = open_file ('logfile.out', action='append')
!   call print_version_number (unit, 'longwave_clouds', version_number)
    if (get_my_pe() == 0) then
      write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
      write (unit,nml=longwave_clouds_nml)
    endif
    call close_file (unit)

!   Lw_control%do_lwcldemiss = do_lwcldemiss

!--------------------------------------------------------------------
!     NLWCLDB =  number of infrared bands with varying cloud
!                properties (emissivity). at present either unity
!                or a model-dependent amount (ifdef lwcldemiss)
!--------------------------------------------------------------------
!   if (Lw_control%do_lwcldemiss) then
!     NLWCLDB = 7
!   else
!     NLWCLDB = 1
!   endif 

!   Lw_parameters%NLWCLDB = NLWCLDB

!--------------------------------------------------------------------


end subroutine longwave_clouds_init




!####################################################################

subroutine cldtau (Cldrad_props, Lw_clouds)
 
!--------------------------------------------------------------------
type(cldrad_properties_type), intent(in) :: Cldrad_props
type(lw_clouds_type), intent(inout) :: Lw_clouds

!---------------------------------------------------------------------
      integer   :: n, k, i, j
      integer  :: is, ie, js, je, ks, ke


      is = 1
      ie = size(Cldrad_props%cmxolw, 1)
      js = 1
      je = size(Cldrad_props%cmxolw, 2)
      ks = 1
!      NLWCLDB = Lw_parameters%NLWCLDB
      NLWCLDB = Cldrad_control%NLWCLDB
      ke = size(Cldrad_props%cmxolw, 3)

!--------------------------------------------------------------------
      allocate (Lw_clouds%taucld_rndlw (IS:IE, JS:JE,KS:KE, NLWCLDB) )
      allocate (Lw_clouds%taunbl_mxolw (IS:IE, JS:JE,KS:KE, NLWCLDB) )
      allocate (Lw_clouds%taucld_mxolw (IS:IE, JS:JE,KS:KE, NLWCLDB) )

!----------------------------------------------------------------------
!    define max overlap layer transmission function over layers KS,KE
!----------------------------------------------------------------------
      do n = 1,NLWCLDB
        do k = KS,KE
          Lw_clouds%taucld_mxolw(:,:,k,n) = 1.0E+00 - Cldrad_props%emmxolw(:,:,k,n)
        enddo
      enddo
 
!----------------------------------------------------------------------
!    define "weighted random cloud" layer transmission function
!    over layers KS,KE
!----------------------------------------------------------------------
      do n = 1,NLWCLDB
        do k = KS,KE
          do j = JS,JE
	    do i = IS,IE
	      if (Cldrad_props%crndlw(i,j,k) > 0.0E+00) then
	        Lw_clouds%taucld_rndlw(i,j,k,n) =   &
                  (Cldrad_props%crndlw(i,j,k)/(1.0E+00 - &
		   Cldrad_props%cmxolw(i,j,k)))*  &
                   (1.0E+00 - Cldrad_props%emrndlw(i,j,k,n)) + &
                   1.0E+00 - Cldrad_props%crndlw(i,j,k)/   &
		   (1.0E+00 - Cldrad_props%cmxolw(i,j,k))
	      else
	        Lw_clouds%taucld_rndlw(i,j,k,n) = 1.0E+00
	      endif
	    enddo
	  enddo
        enddo
      enddo
 
!--------------------------------------------------------------------
!     define "nearby layer" cloud transmission function for max
!     overlapped clouds (if emissivity not equal to one)
!--------------------------------------------------------------------
      do n = 1,NLWCLDB
        do k=KS,KE
          Lw_clouds%taunbl_mxolw(:,:,k,n) = 0.0E+00
        enddo
      enddo
 
!------------------------------------------------------------------

  
end subroutine cldtau




!######################################################################

subroutine cloud (kl, Cldrad_props, Lw_clouds, cldtf)

!---------------------------------------------------------------------
type(cldrad_properties_type), intent(in) :: Cldrad_props
type(lw_clouds_type), intent(in) :: Lw_clouds
real, dimension(:,:,:,:),   intent(out) :: cldtf        
integer,                    intent(in)  :: kl
!--------------------------------------------------------------------

    integer   ::   n, i, j, kp
    integer   :: is, ie, js, je, ks, ke

    real, dimension(size(Cldrad_props%cmxolw,1), &
                    size(Cldrad_props%cmxolw,2), & 
		    size(Cldrad_props%cmxolw,3) + 1) ::   &
		                         cldtfmo, cldtfrd


    is = 1
    ie = size (Cldrad_props%cmxolw, 1)
    js = 1
    je = size(Cldrad_props%cmxolw,2)
    ks = 1
    ke = size(Cldrad_props%cmxolw, 3)

!---------------------------------------------------------------------
!    the definition of "within a max overlapped cloud" is:
!    at pressure level k (separating layers k and (k-1)), the max
!    overlap cloud amounts for layers k and (k-1) must be 1) nonzero
!    (else no such cloud) and 2) equal (else one such cloud ends at
!    level k and another begins). Another way to define this is: if
!    considering the transmission across layer kp (between levels
!    kp and (kp+1)) the max overlap cloud amounts for layers kp and
!    (kp-1) must be nonzero and equal.
!---------------------------------------------------------------------
    do n = 1,NLWCLDB

!--------------------------------------------------------------------
!   cloud "nearby layer" transmission functions
!--------------------------------------------------------------------
      cldtfmo(:,:,kl) = 0.0
      cldtfrd(:,:,kl) = 1.0

!--------------------------------------------------------------------
!   if level kl is within a maximum overlapped cloud, the cloud
!   "nearby layer" transmission function may be non-unity. Exception:
!   at levels KS,KE+1  the function must be unity.
!--------------------------------------------------------------------
      if (kl > KS .AND. kl < KE+1) then
        do j=JS,JE
          do i=IS,IE
            if ( Cldrad_props%cmxolw(i,j,kl-1) /= 0.0 .and.     &
                 Cldrad_props%cmxolw(i,j,kl) == Cldrad_props%cmxolw(i,j,kl-1)) then
              cldtfmo(i,j,kl) = Cldrad_props%cmxolw(i,j,kl)*Lw_clouds%taunbl_mxolw(i,j,kl,n)
              cldtfrd(i,j,kl) = 1.0 - Cldrad_props%cmxolw(i,j,kl)
            endif
          enddo
        enddo
      endif

      cldtf(:,:,kl,n) = cldtfmo(:,:,kl) + cldtfrd(:,:,kl)

!--------------------------------------------------------------------
!     cloud transmission functions between level kl and higher
!     levels ( when kl le KE+1)
!--------------------------------------------------------------------
      if (kl .LT. KE+1) then
        cldtfmo(:,:,kl) = 0.0
        cldtfrd(:,:,kl) = 1.0

!--------------------------------------------------------------------
!    for first layer below  level kl, assume flux at level kl
!   is unity and is apportioned between (cmxolw) max. overlap cld,
!   (crndlw) rnd overlap cld, and remainder as clear sky.
!--------------------------------------------------------------------
        cldtfmo(:,:,kl+1) = Cldrad_props%cmxolw(:,:,kl)*Lw_clouds%taucld_mxolw(:,:,kl,n)
        cldtfrd(:,:,kl+1) = (1.0 - Cldrad_props%cmxolw(:,:,kl))*  &
                            Lw_clouds%taucld_rndlw(:,:,kl,n)
        cldtf(:,:,kl+1,n) = cldtfmo(:,:,kl+1) + cldtfrd(:,:,kl+1)

!--------------------------------------------------------------------
!     if layers above and below level (kp-1) have no max overlap cloud,
!     or their amounts differ (ie, either top of new max overlap cld or
!     no max overlap cld at all), then apportion total "flux" (or,
!     cloud tf (cldtf)) between any max overlap cloud in layer(kp-1),
!     any rnd overlap cloud and clear sky.
!--------------------------------------------------------------------
        do kp = kl+2, KE+1
          do j=JS,JE
	    do i=IS,IE
	      if (Cldrad_props%cmxolw(i,j,kp-2) .eq. 0. .or.    &
                  Cldrad_props%cmxolw(i,j,kp-2) .ne. Cldrad_props%cmxolw(i,j,kp-1)) then
	        cldtfmo(i,j,kp) = cldtf(i,j,kp-1,n)*   &
                             Cldrad_props%cmxolw(i,j,kp-1)*Lw_clouds%taucld_mxolw(i,j,kp-1,n)
	        cldtfrd(i,j,kp) = cldtf(i,j,kp-1,n)*   &
                    (1.0 - Cldrad_props%cmxolw(i,j,kp-1))*Lw_clouds%taucld_rndlw(i,j,kp-1,n)
	        cldtf(i,j,kp,n) = cldtfmo(i,j,kp) + cldtfrd(i,j,kp)

!--------------------------------------------------------------------
!    if layer above level (kp-1) has a max overlap cloud, and layer
!    layer below level (kp-1) also does (ie, within max overlap cld)
!    obtain separate cloud tfs for max overlap cloud and for 
!    remainder (which may contain a random overlap cloud).
!--------------------------------------------------------------------
              else 
	        cldtfmo(i,j,kp) = cldtfmo(i,j,kp-1)*   &
                                  Lw_clouds%taucld_mxolw(i,j,kp-1,n)
                cldtfrd(i,j,kp) = cldtfrd(i,j,kp-1)*   &
                                  Lw_clouds%taucld_rndlw(i,j,kp-1,n)
                cldtf(i,j,kp,n) = cldtfmo(i,j,kp) + cldtfrd(i,j,kp)
              endif
	    enddo
	  enddo
        enddo
      endif
    enddo

end  subroutine cloud



!####################################################################

subroutine thickcld (pflux_in, Cldrad_props, Lw_output)

!------------------------------------------------------------------
real,   dimension (:,:,:),   intent(in)    ::  pflux_in
type(cldrad_properties_type), intent(in) :: Cldrad_props
type(lw_output_type), intent(inout) :: Lw_output      

!----------------------------------------------------------------------
!     input variables:
!
!     cmxolw  =  amounts of maximally overlapped longwave clouds in
!                layers from KS to KE.
!
!     emmxolw =  longwave cloud emissivity for maximally overlapped
!                clouds through layers from KS to KE. (default is one).
!
!     kmxolw  =  maximum number of maximally overlapped longwave clouds.
!
!     pdfinv  =  inverse of pressure difference between flux levels.
!
!     pflux   =  pressure at flux levels of model.
!
!----------------------------------------------------------------------
!     output variables:
!
!     flxnet  =  net flux at flux levels of model (also an input
!                variable)
!
!     heatra  =  heating rate at data levels. (also an input variable)
!----------------------------------------------------------------------

      real, dimension (size(pflux_in,1), size(pflux_in,2))  :: &
                                delptc, fbtm, ftop, pbtm, ptop
      integer, dimension (size(pflux_in,1), size(pflux_in,2))  :: &
                            itopmxo, ibtmmxo
      integer, dimension (size(pflux_in,1), size(pflux_in,2),  &
                          size(pflux_in,3)-1)  :: &
                               ktopmxo, kbtmmxo
      real, dimension (size(pflux_in,1), size(pflux_in,2),  &
                          size(pflux_in,3)-1)  :: &
!                                tmp1
                                 tmp1, pdfinv
      real, dimension (size(pflux_in,1), size(pflux_in,2),  &
                          size(pflux_in,3)  )  :: pflux


      integer :: i,j, k, kc, kc1, kc2
    integer   :: is, ie, js, je, ks, ke
 integer                                    ::  kmxolw

!--------------------------------------------------------------------

       kmxolw = MAXVAL(Cldrad_props%nmxolw)

    is = 1
    ie = size (pflux_in, 1)
    js = 1
    je = size(pflux_in,2)
    ks = 1
    ke = size(pflux_in, 3) - 1

! convert pflux to cgs
    pflux = 10.0*pflux_in


    pdfinv(:,:,ks:ke) = 1.0 / (pflux(:,:,ks+1:ke+1) - pflux(:,:,ks:ke))

!----------------------------------------------------------------------
!     this module recomputes cloud fluxes in "thick" clouds assuming
!     that df/dp is constant. the effect is to reduce top-of-cloud
!     cooling rates, thus performing a "pseudo-convective adjustment"
!     by heating (in a relative sense) the cloud top.
!
!     NOTE: this module cannot handle a frequency-dependent emissivity.
!     Therefore, it assumes that emissivity quantities (emmxolw) are
!     from frequency band 1 (normally unity).
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   determine levels at which max overlap clouds start and stop
!--------------------------------------------------------------------
      itopmxo(:,:) = 0
      ibtmmxo(:,:) = 0
      ktopmxo(:,:,:) = 0
      kbtmmxo(:,:,:) = 0

!--------------------------------------------------------------------
!   max overlap cloud in first layer (not likely)
!--------------------------------------------------------------------
      do j = JS,JE
        do i = IS,IE
          if ( Cldrad_props%cmxolw(i,j,KS) .GT. 0.0E+00) then
                 itopmxo(i,j) = itopmxo(i,j) + 1
		 ktopmxo(i,j,itopmxo(i,j)) = KS
	  endif
	enddo
      enddo

!--------------------------------------------------------------------
!   k-level for which top of max overlap cloud is defined
!--------------------------------------------------------------------
      do k = KS+1,KE
        do j = JS,JE
          do i = IS,IE
            if (Cldrad_props%cmxolw(i,j,k) .GT. 0.0E+00 .AND.   &
                Cldrad_props%cmxolw(i,j,k-1) .NE. Cldrad_props%cmxolw(i,j,k)) then
                  itopmxo(i,j) = itopmxo(i,j) + 1
	          ktopmxo(i,j,itopmxo(i,j)) = k
            endif
	  enddo
        enddo
      enddo

!--------------------------------------------------------------------
!   k-level for which bottom of max overlap cloud is defined
!--------------------------------------------------------------------
      do k = KS,KE-1
        do j = JS,JE
          do i = IS,IE
            if (CLdrad_props%cmxolw(i,j,k) .GT. 0.0E+00 .AND.    &
                Cldrad_props%cmxolw(i,j,k+1) .NE. Cldrad_props%cmxolw(i,j,k)) then
	      ibtmmxo(i,j) = ibtmmxo(i,j) + 1
	      kbtmmxo(i,j,ibtmmxo(i,j)) = k+1
            endif
          enddo
        enddo
      enddo

!--------------------------------------------------------------------
!   bottom of max overlap cloud in KE'th level
!--------------------------------------------------------------------
     do j = JS,JE
       do i = IS,IE
         if (Cldrad_props%cmxolw(i,j,KE) .GT. 0.0E+00) then
           ibtmmxo(i,j) = ibtmmxo(i,j) + 1
	   kbtmmxo(i,j,ibtmmxo(i,j)) = KE+1
         endif
       enddo
     enddo
 
!---------------------------------------------------------------------- 
!    obtain the pressures and fluxes of the top and bottom of the cloud.
!---------------------------------------------------------------------- 
      if (kmxolw .NE. 0) then
        do kc=1,kmxolw 
          do j=JS,JE
            do i=IS,IE
	      if (kbtmmxo(i,j,kc) > ktopmxo(i,j,kc)) then
                kc1 = ktopmxo(i,j,kc)
                kc2 = kbtmmxo(i,j,kc)
                ptop(i,j) = pflux (i,j,kc1) 
                pbtm(i,j) = pflux (i,j,kc2)
                ftop(i,j) = Lw_output%flxnet(i,j,kc1)
                fbtm(i,j) = Lw_output%flxnet(i,j,kc2)

!-----------------------------------------------------------------------
!      compute the "flux derivative" df/dp delptc.
!-----------------------------------------------------------------------
                delptc(i,j) = (ftop(i,j) - fbtm(i,j))/   &
                              (ptop(i,j) - pbtm(i,j))
!-----------------------------------------------------------------------
!      compute the total flux change from the top of the cloud.
!-----------------------------------------------------------------------
                do k=kc1+1,kc2-1
                  tmp1(i,j,k) = ftop(i,j) + (pflux(i,j,k) - ptop(i,j))*&
				delptc(i,j) 
                  Lw_output%flxnet(i,j,k) = Lw_output%flxnet(i,j,k)*(1.0E+00 -    &
				  Cldrad_props%cmxolw(i,j,k)*Cldrad_props%emmxolw(i,j,k,1)) +  &
				  tmp1(i,j,k)*Cldrad_props%cmxolw(i,j,k)*   &
				  Cldrad_props%emmxolw(i,j,k,1)
                end do
              endif
            end do
          end do
        end do
      endif

!-----------------------------------------------------------------------
!     recompute the heating rates based on the revised fluxes.
!-----------------------------------------------------------------------
      Lw_output%heatra(:,:,KS:KE) = radcon*(Lw_output%flxnet(:,:,KS+1:KE+1) -   &
                          Lw_output%flxnet(:,:,KS:KE))*pdfinv(:,:,KS:KE)




end subroutine thickcld
  
   

!#####################################################################


                   end module longwave_clouds_mod

