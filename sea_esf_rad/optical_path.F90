
                    module optical_path_mod
 

use  longwave_params_mod,   only: NBLW, NBCO215, NBLY_ORIG
use  utilities_mod,         only: open_file, file_exist,    &
                                  check_nml_error, error_mesg, &
                                  print_version_number, FATAL, NOTE, &
				  WARNING, get_my_pe, close_file
use rad_utilities_mod,      only: looktab, longwave_tables3_type, &
				  radiative_gases_type, &
				  atmos_input_type, &
				  Environment_type, Environment, &
			                       Lw_parameters, &
			          longwave_parameter_type, &
				  optical_path_type, &
				  gas_tf_type, &
				  table_alloc, longwave_control_type,&
				  Lw_control
!use ch4_n2o_mod,            only: ch4_n2o_input_type, CN_basic
!use cfc_mod,                only: cfc_optical_depth, cfc_exact,&
!use cfc_mod,                only:                    cfc_exact,&
use lw_gases_stdtf_mod,     only:                    cfc_exact,&
				  cfc_overod, cfc_overod_part,   &
!			  cfc_exact_part, cfc_dealloc
				  cfc_exact_part
use constants_new_mod,      only: wtmair, avogno, wtmh2o, grav,    &
				  rgas, pstd, diffac, rh2oair
!use longwave_aerosol_mod,   only: aertau, longwave_aerosol_dealloc, &
use longwave_aerosol_mod,   only: aertau, &
                                  longwave_aerosol_init
!                                 get_totaerooptdep,   &
!		          get_totaerooptdep_15,   &
!			  get_aerooptdep_KE_15
!use gas_tf_mod,             only: get_gastf_for_optpath

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                   optical depth calculation module
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

   character(len=128)  :: version =  '$Id: optical_path.F90,v 1.3 2002/07/16 22:36:15 fms Exp $'
   character(len=128)  :: tag     =  '$Name: havana $'


!---------------------------------------------------------------------
!----- interfaces  -----
           
public     &
!	      optical_path_init, optical_dealloc,   &
	      optical_path_init,                    &
	      optical_path_setup,     &
	      optical_trans_funct_from_KS,    &
	      optical_trans_funct_k_down, &
	      optical_trans_funct_KE,  &
	      optical_trans_funct_diag, &
	      get_totch2o, get_totch2obd, &
	      get_totvo2         


private    &
	      optical_ckd2p1_init, optical_path_ckd2p1, &
	      optical_o3, optical_rbts, optical_h2o, cfc_optical_depth
!---------------------------------------------------------------------
!---- namelist   -----

real                :: dummy = 0.0


 namelist / optical_path_nml /    &
				dummy



!-------------------------------------------------------------------
!-----  public data ----




!---------------------------------------------------------------------
!----- private  data  ----

!--------------------------------------------------------------------
!    data from former block data bs296 for self-broadened continuum
!    at 296K, band-integrated, in 5 - 19995 cm-1 range.
!               06/28/82
!               units of (cm**3/mol) * 1.E-20
!-------------------------------------------------------------------
real       :: v1sh2o_296, v2sh2o_296, dvsh2o_296,   &
              ssh2o_296(2000)
integer    :: nptsh2o_296

!--------------------------------------------------------------------
!  data from former block data bfh2o for foreign-broadened continuum
!    band-integrated, in 5 - 19995 cm-1 range.
!               06/28/82
!               units of (cm**3/mol) * 1.E-20
!--------------------------------------------------------------------
real        ::  v1fh2o, v2fh2o, dvfh2o, sfh2o(2000)
integer     ::  nptfh2o

!--------------------------------------------------------------------
!    array sfac is the frequency-dependent multiplicative factor used
!    to change the original self-broadened continuum coefficients
!    to those used in ckd2.1 (including intermediate changes).
!
!    array fscal is the frequency-dependent multiplicative factor used
!    to change the original foreign-broadened continuum coefficients
!    to those used in ckd2.1 (including intermediate changes).
!
!    array tmpfctrs is the logarithmic temperature dependence (per K)
!    of the self-broadened continuum coefficient, as a function of
!    frequency, used in all AFGL continuum models, including ckd2.1.
!    the frequency ranges and intervals are as in sh2o_296.
!-----------------------------------------------------------------------
real         :: sfac(2000), fscal(2000), tmpfctrs(2000)

!----------------------------------------------------------------------
!         the radfunc function (1 - exp(-h*nu/kt))/(1 + exp(-h*nu/kt))
!    is tabulated from 5 to 2995 cm-1 at intervals of 10 cm-1,
!    and from 100K to 490K at 10K intervals. note that the
!    radfn function used in ckd2.1 equals the radfunc function
!    defined above, multiplied by nu (in cm-1).
!        the temperature derivative (at 105K to 485K, with the final
!    array value set to zero) is obtained from radfunc, and stored in
!    radfuncderiv.
!        tktab and vjtab are the respective temperature and frequency
!    points at which tabulations occurred.
!----------------------------------------------------------------------
type (longwave_tables3_type)  :: radfunc
integer                       :: ioffh2o, nptch2o 
real                          :: vvj(2000)

!---------------------------------------------------------------------
!        fvj = foreign-broadened ckd 2.1 coefficient (including
!              all corrections), averaged over 7 specified wide
!              frequency bands in the 560-1200 cm-1 range. The average
!              is weighted by the frequency of the individual 10 cm-1
!              bands used in the averaging process.
!     fvjinw = band-averaged foreign coefficient (as in fvj) over
!              the 900-990,1070-1200 cm-1 range.
!      fvjwd = band-averaged foreign coefficient (as in fvj) over
!              the 560-800 cm-1 range.
!        svj = self-broadened ckd 2.1 coefficient (including
!              all corrections), averaged over 7 specified wide
!              frequency bands in the 560-1200 cm-1 range. The average
!              is weighted by the frequency of the individual 10 cm-1
!              bands used in the averaging process.
!     fvjinw = band-averaged self coefficient (as in svj) over
!              the 900-990,1070-1200 cm-1 range.
!      svjwd = band-averaged self coefficient (as in svj) over
!              the 560-800 cm-1 range.
!    radfnbd = the radiation function (radfn) averaged over each of
!              the 7 frequency bands: assumed to be altitude-independent
! radfnbdinw = same as radfnbd, but for the 560-800 cm-1 range.
!  radfnbdwd = same as radfnbd, but for the 900-990,1070-1200 cm-1 range
!----------------------------------------------------------------------
real      ::    fvj(7), fvjinw, fvjwd, svj(7), svjinw, svjwd,    &
                radfnbd(7), radfnbdinw, radfnbdwd

real      ::    ao3rnd(3), bo3rnd(3)
data ao3rnd    / 0.516411E+02,  0.237768E+04,  0.331449E+02/
data bo3rnd    /  0.670035E+01,  0.955445E+01,  0.856424E+01/

real      :: betawd, ab15wd 

!! NOTE: this variable not used     real    :: betinw
!----------------------------------------------------------------------
!   2) 800-990, 1070-1200 as 1 band. program gasbnd is used as 2
!    bands (800-990,1070-1200) with ifdef comb2 on.
!----------------------------------------------------------------------
!     data betinw  /0.766811E+01/

!---------------------------------------------------------------------
!  define continuum coefficients over special bands, the choices
!  depend on model architecture. the program gasbnd is used.
!
!    1) 560-800 as 1 band
!--------------------------------------------------------------------
data betawd   /   0.347839E+02/

!----------------------------------------------------------------------
!    3) 160-560 (as 8 bands using combined bands). program gasbnd is
!    used as 40 bands (160-560,10 cm-1 bandwidth) with ifdef icomb on.
!--------------------------------------------------------------------
integer                               :: n
real, dimension (NBLY_ORIG)           :: betacm
data (betacm(n),n=1,8)  /      &
      0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,   &
      0.188625E+03,  0.144293E+03,  0.174098E+03,  0.909366E+02/

!---------------------------------------------------------------------
!    4) 560-1200 and 4.3 um band (8 bands, frequency range given
!    by bdlocm-bdhicm). program gasbnd is used with 8 specified
!    bandwidths.
!---------------------------------------------------------------------
data (betacm(n),n=9,16) /    &
      0.565430E+02,  0.343645E+02,  0.198461E+02,  0.113124E+02,   &
      0.754174E+01,  0.589554E+01,  0.495227E+01,  0.000000E+00/

!---------------------------------------------------------------------

  real, allocatable, dimension (:,:)     ::             csfah2o

integer   :: NBTRG, NBTRGE

      integer :: israd, ierad, jsrad, jerad, ks, ke

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
!---------------------------------------------------------------------
        real, dimension(2,10)  :: cpf10h2o, csf10h2o
        real, dimension(2, 4)  :: cpf4h2o, csf4h2o
        real, dimension(2, 2)  :: cpf2h2o, csf2h2o
        real, dimension(2   )  :: cpf1h2o, csf1h2o

!---------------------------------------------------------------------
!   data for h2o for 10 20 cm-1 wide bands in 1200-1400 cm-1 range
!---------------------------------------------------------------------

      data          cpf10h2o            / &
     &   0.189824E-01, -0.660707E-04, &
     &   0.175350E-01, -0.275910E-04, &
     &   0.142808E-01, -0.274332E-04, &
     &   0.142165E-01, -0.415313E-04, &
     &   0.260793E-01, -0.524774E-04, &
     &   0.182821E-01, -0.573604E-04, &
     &   0.199166E-01, -0.721160E-04, &
     &   0.147737E-01, -0.477693E-04, &
     &   0.124860E-01, -0.483700E-04, &
     &   0.951027E-02, -0.441301E-04/

      data          csf10h2o            / &
     &   0.179076E-01, -0.493954E-04, &
     &   0.188315E-01, -0.282847E-04, &
     &   0.165058E-01, -0.155592E-04, &
     &   0.163390E-01, -0.389195E-04, &
     &   0.198982E-01, -0.380731E-04, &
     &   0.156372E-01, -0.462165E-04, &
     &   0.143535E-01, -0.410551E-04, &
     &   0.116195E-01, -0.260344E-04, &
     &   0.100036E-01, -0.367164E-04, &
     &   0.846707E-02, -0.385911E-04/


!---------------------------------------------------------------------
!   data for h2o for 4 50 cm-1 wide bands in 1200-1400 cm-1 range
!---------------------------------------------------------------------
      data          cpf4h2o /    &
     &   0.166392E-01, -0.428283E-04,  &
     &   0.175692E-01, -0.295058E-04,  &
     &   0.175389E-01, -0.566292E-04, &
     &   0.105481E-01, -0.443731E-04/

      data          csf4h2o /  &
     &   0.184822E-01, -0.337057E-04,  &
     &   0.177179E-01, -0.351393E-04,  &
     &   0.140172E-01, -0.386980E-04,  &
     &   0.938797E-02, -0.366436E-04/

!---------------------------------------------------------------------
!   data for h2o for 2 100 cm-1 wide bands in 1200-1400 cm-1 range
!---------------------------------------------------------------------
      data          cpf2h2o /  &
     &   0.174097E-01, -0.317609E-04,  &
     &   0.114313E-01, -0.432460E-04/

      data          csf2h2o /   &
     &   0.179228E-01, -0.347266E-04,  &
     &   0.109031E-01, -0.361374E-04/

!---------------------------------------------------------------------
!   data for h2o for 1 200 cm-1 wide band in 1200-1400 cm-1 range
!---------------------------------------------------------------------
      data          cpf1h2o /   &
     &   0.115749E-01, -0.425616E-04/

      data          csf1h2o /    &
     &   0.119546E-01, -0.343610E-04/



contains





!subroutine optical_path_init
subroutine optical_path_init (kmax, latb)

integer, intent(in) :: kmax
real, dimension(:), intent(in) :: latb

!--------------------------------------------------------------------
       real, dimension (NBLW)  :: dummy
       real, dimension (NBLW)  :: ap, bp, atp, btp
       real                    :: awide_c, bwide_c, awide_n, bwide_n, &
				  awide, bwide
       integer                 ::  unit, ierr, io
       integer                 :: inrad, k, m

!---------------------------------------------------------------------
!    define random band parameters for special bands. the choices 
!    depend on model architecture. the program gasbnd is used.
!
!    1)  560-800 as 1 band
!---------------------------------------------------------------------
      data awide_c              /  &
         0.308562E+01/
      data bwide_c              /  &
         0.503352E-01/

      data awide_n              /  &
        0.309494E+01/
      data bwide_n              /  &
        0.503391E-01/ 

!---------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=optical_path_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'optical_path_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'optical_path', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=optical_path_nml)
      endif
      call close_file (unit)

!---------------------------------------------------------------------
      NBTRG  = Lw_parameters%NBTRG
      NBTRGE = Lw_parameters%NBTRGE

      if (Lw_control%do_ckd2p1 ) then
	inrad = open_file ('INPUT/id2h2obdckd2p1', form= 'formatted', &
			    action='read')
      else
	inrad = open_file ('INPUT/id2h2obdfull', form= 'formatted', &
			    action='read')
      endif

      read (inrad,9000) (dummy(k),k=1,NBLW)
      read (inrad,9000) (dummy(k),k=1,NBLW)
      read (inrad,9000) (ap(k),k=1,NBLW)
      read (inrad,9000) (bp(k),k=1,NBLW)
      read (inrad,9000) (atp(k),k=1,NBLW)
      read (inrad,9000) (btp(k),k=1,NBLW)
9000  format(5e14.6)

      call close_file (inrad)

!---------------------------------------------------------------------
      call longwave_aerosol_init (kmax, latb)

      if (Lw_control%do_ckd2p1 ) then
        awide = awide_c
        bwide = bwide_c
      else
        awide = awide_n
        bwide = bwide_n
      endif

!---------------------------------------------------------------------
!    compute a*b for computational frequency bands for the 15 um
!    region, as 1 band (ab15wd)
!---------------------------------------------------------------------
      ab15wd = awide*bwide


    if (Lw_control%do_ch4_n2o) then
      allocate ( csfah2o(2, NBTRGE) )

!----------------------------------------------------------------------
!     select, from random band coefficients in 1200-1400 cm-1 range,
!     those appropriate for NBTRGE h2o bands.
!---------------------------------------------------------------------
      if (NBTRGE .EQ. 1) then
        do m=1,NBTRGE
          csfah2o(1,m) =          csf1h2o(1)
          csfah2o(2,m) =          csf1h2o(2)
	end do
      elseif (NBTRGE .EQ. 2) then
        do m=1,NBTRGE
          csfah2o(1,m) =          csf2h2o(1,m)
          csfah2o(2,m) =          csf2h2o(2,m)
        end do
      elseif (NBTRGE .EQ. 4) then
        do  m=1,NBTRGE
          csfah2o(1,m) =          csf4h2o(1,m)
          csfah2o(2,m) =          csf4h2o(2,m)
        end do
      elseif (NBTRGE .EQ. 10) then
        do m=1,NBTRGE
          csfah2o(1,m) =          csf10h2o(1,m)
          csfah2o(2,m) =          csf10h2o(2,m)
        enddo
      elseif (NBTRGE .EQ. 20) then
        do m=1,NBTRGE
          csfah2o(1,m) = atp(m+120)
          csfah2o(2,m) = btp(m+120)
        enddo
      else
	call error_mesg ('optical_path_init', &
	   'NBTRGE is inconsistent with available data', FATAL)
      endif
    endif

!------------------------------------------------------------------
    if (Lw_control%do_ckd2p1) then
      call optical_ckd2p1_init
    endif

!------------------------------------------------------------------


end subroutine optical_path_init





!###################################################################

!subroutine optical_path_setup (is, ie, js, je, kx, Atmos_input, &
subroutine optical_path_setup (is, ie, js, je,  Atmos_input, &
                                Rad_gases, Optical) 

!integer, intent(in) ::  is, ie, js, je, kx
integer, intent(in) ::  is, ie, js, je
type(radiative_gases_type), intent(in) :: Rad_gases
type(atmos_input_type), intent(in) :: Atmos_input
type(optical_path_type), intent(inout) :: Optical     

!---------------------------------------------------------------------
     integer      :: k, i
     integer      :: ix, jx, kx
real, dimension (size(Atmos_input%press,1), size(Atmos_input%press,2), &
                        size(Atmos_input%press,3) )  :: press, pflux, &
                                                        temp, tflux, &
                                                    atmden, vv

real, dimension (size(Atmos_input%press,1), size(Atmos_input%press,2), &
                  size(Atmos_input%press,3)-1 )  ::  rh2o, deltaz

     ix = ie -is + 1
     jx = je -js +1
     israd = 1
     ierad = ix
     jsrad = 1
     jerad = jx
     ks    = 1
     kx = size(Atmos_input%press,3) - 1
     ke    = kx
!  convert press and pflux to cgs.
        press(:,:,:) = 10.0*Atmos_input%press(:,:,:)
        pflux(:,:,:) = 10.0*Atmos_input%pflux(:,:,:)
      deltaz = Atmos_input%deltaz
      temp = Atmos_input%temp
      rh2o = Atmos_input%rh2o
      tflux = Atmos_input%tflux
!--------------------------------------------------------------------
!     atmden   =  atmospheric density, in gm/cm**2, for each of the
!                 KMAX layers.
!-------------------------------------------------------------------
!     allocate (Optical%sh2o        (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
!     allocate (Optical%tmpexp      (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
!     allocate (Optical%rvh2o       (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (Optical%wk          (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
!     allocate (Optical%rhoave      (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (Optical%rh2os       (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (Optical%rfrgn       (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (Optical%tfac        (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (Optical%avephi      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )

     if (Lw_control%do_ch4_n2o) then
       allocate (Optical%avephif   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
     endif
!---------------------------------------------------------------------
 
!----------------------------------------------------------------------
!     define the layer-mean pressure in atmospheres (vv) and the layer 
!     density (atmden). 
!----------------------------------------------------------------------
     do k=KS,KE
       atmden(:,:,k) = (pflux(:,:,k+1) - pflux(:,:,k))/grav
       vv(:,:,k)     = 0.5E+00*(pflux(:,:,k+1) + pflux(:,:,k)  )/pstd
     enddo

!----------------------------------------------------------------------
!     compute optical paths.
!----------------------------------------------------------------------
     call optical_h2o (pflux, atmden, vv, press, temp, rh2o, tflux, Optical)

!---------------------------------------------------------------------
     if (.not. Lw_control%do_ckd2p1) then
       call optical_rbts  (temp, rh2o, Optical)
     else
!---------------------------------------------------------------------
!    call optical_ckd2p1 to determine self- and foreign-broadened h2o
!    continuum paths, for the given temperature, pressure and mixing
!    ratio, over the predetermined frequency range.
!---------------------------------------------------------------------
       call optical_path_ckd2p1  (atmden, press, temp, rh2o, Optical)
     endif

!---------------------------------------------------------------------
     call optical_o3 (atmden, Rad_gases%qo3, vv, Optical)

!--------------------------------------------------------------------
     if (Lw_control%do_cfc) then
       call cfc_optical_depth (atmden, Rad_gases, Optical) 
     endif

!---------------------------------------------------------------------
!     compute aerosol layer transmission functions for all layers.
!     option predlwaer is planned, but not yet available. when  it 
!     becomes available,  aeralb, aerext and aerconc will be additional 
!     arguments going to aertau.
!---------------------------------------------------------------------
     if (Lw_control%do_lwaerosol) then
       call aertau (ix, jx, kx, js, deltaz, Optical) 
     endif
!---------------------------------------------------------------------
       

end subroutine  optical_path_setup



!####################################################################

subroutine optical_trans_funct_from_KS (Gas_tf, to3cnt, overod, Optical, &
				        cnttaub1, cnttaub2, cnttaub3)

!---------------------------------------------------------------------
real, dimension (:,:,:), intent(out) ::  to3cnt, overod, &
                                         cnttaub1, cnttaub2, cnttaub3
type(optical_path_type), intent(inout) :: Optical
type(gas_tf_type), intent(inout) :: Gas_tf 

!---------------------------------------------------------------------
   real, dimension (size(to3cnt,1), size(to3cnt,2), &
                 size(to3cnt,3)) ::   &
                                              tmp1, tmp2, tmp3,   &
					              totch2o_tmp,  &
					      totaer_tmp, tn2o17,  &
					      totaerooptdep_15

   real, dimension (size(to3cnt,1), size(to3cnt,2), &
                 size(to3cnt,3)-1) :: cfc_tf

      integer    :: m


!-----------------------------------------------------------------------
!   compute transmission functions in 990-1070 cm-1 range, including
!   ozone and h2o continuum, from level KS to all other levels. 
!------------------------------------------------------------------
      tmp1  (:,:,KS:KE) = bo3rnd(2)*Optical%tphio3(:,:,KS+1:KE+1)/  &
                           Optical%toto3(:,:,KS+1:KE+1)
      tmp2(:,:,KS:KE) = 0.5*(tmp1(:,:,KS:KE)*(SQRT(1.0E+00 + (4.0E+00*&
                        ao3rnd(2)*Optical%toto3(:,:,KS+1:KE+1))/  &
			tmp1(:,:,KS:KE))  - 1.0E+00))
      if (Lw_control%do_ckd2p1) then
	call get_totch2obd(6, Optical, totch2o_tmp)
	tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) + diffac*   &
		          totch2o_tmp(:,:,KS+1:KE+1)
      else 
        tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) + betacm(14)*   &
                          Optical%totvo2(:,:,KS+1:KE+1)
      endif

      if (Lw_control%do_lwaerosol) then
!	call get_totaerooptdep(6, Optical, totaer_tmp)
	totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,6)
	tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +    &
                          totaer_tmp   (:,:,KS+1:KE+1)
      endif
!     to3cnt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))
      to3cnt(:,:,KS) = 1.0
      to3cnt(:,:,KS+1:KE+1) = EXP(-1.0E+00*tmp2(:,:,KS:KE))
 
!--------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in to3cnt.
!---------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_exact (6, Optical, cfc_tf)
!       to3cnt(:,:,KS:KE) = to3cnt(:,:,KS:KE)* cfc_tf(:,:,KS:KE)
        to3cnt(:,:,KS+1:KE+1) = to3cnt(:,:,KS+1:KE+1)* cfc_tf(:,:,KS:KE)
      endif

!---------------------------------------------------------------------
!   compute transmission function in the 560-800 cm-1 range
!   evaluate  optical depth contributions 
!   add contributions from h2o(lines) and h2o(continuum).
!   h2o(continuum) contributions are either Roberts or CKD2.1
!---------------------------------------------------------------------
      tmp1(:,:,KS:KE) = SQRT(ab15wd*Optical%totphi(:,:,KS+1:KE+1)) 
      if (Lw_control%do_ckd2p1) then
	tmp1(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*   &
                          Optical%totch2obdwd(:,:,KS+1:KE+1)
      else
	tmp1(:,:,KS:KE) = tmp1(:,:,KS:KE) + betawd*   &
                          Optical%totvo2     (:,:,KS+1:KE+1)
      endif

!--------------------------------------------------------------------
!   add contribution from longwave aerosols (if desired).
!--------------------------------------------------------------------
      if (Lw_control%do_lwaerosol) then
!	call get_totaerooptdep_15 (Optical, totaerooptdep_15)
	totaerooptdep_15(:,:,:) = Optical%totaerooptdep_15(:,:,:)
	tmp1(:,:,KS:KE) = tmp1(:,:,KS:KE) +    &
                          totaerooptdep_15(:,:,KS+1:KE+1)
      endif
 
!----------------------------------------------------------------------
!    compute transmission function due to these contributions. the
!    effects of co2, n2o  and  cfc's (not exponentials) are added
!    later.
!---------------------------------------------------------------------
!     overod(:,:,KS:KE) = EXP(-1.0E+00*tmp1     (:,:,KS:KE))
      overod(:,:,KS) = 1.0
      overod(:,:,KS+1:KE+1) = EXP(-1.0E+00*tmp1     (:,:,KS:KE))

!---------------------------------------------------------------------
!   add contribution from the 17 um n2o band (if desired).
!   the expression with tn2o17 retains the 560-630 cm-1 equi-
!   valent widths in evaluating 560-800 cm-1 transmissivities.
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
	  tn2o17(:,:,ks+1:ke+1) = Gas_tf%tn2o17(:,:,ks+1:ke+1)

        if (NBCO215 .EQ. 2) then
!         overod(:,:,KS:KE) = overod(:,:,KS:KE) *    &
          overod(:,:,KS+1:KE+1) = overod(:,:,KS+1:KE+1) *    &
                              (130./240. + 110./240.*   &
			       tn2o17(:,:,KS+1:KE+1))
        elseif (NBCO215 .EQ. 3) then
!         overod(:,:,KS:KE) = overod(:,:,KS:KE) * (170./240. +    &
          overod(:,:,KS+1:KE+1) = overod(:,:,KS+1:KE+1) * (170./240. +    &
			      70./240.*tn2o17(:,:,KS+1:KE+1))
        endif
      endif

!--------------------------------------------------------------------- 
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod .
!--------------------------------------------------------------------
      if (Lw_control%do_cfc) then
	call cfc_overod (Optical, cfc_tf)
!overod(:,:,KS:KE) = overod(:,:,KS:KE)*cfc_tf(:,:,KS:KE)
	overod(:,:,KS+1:KE+1) = overod(:,:,KS+1:KE+1)*cfc_tf(:,:,KS:KE)
      endif 

!----------------------------------------------------------------------
!     compute continuum band transmission functions from level KS to
!     other levels (cnttau). the continuum transmission function from
!     level k to kp (contod) equals cnttau for k=KS, so is not
!     evaluated here. for all other levels k, contod is obtained by
!     division of relevant values of cnttau.
!---------------------------------------------------------------------
      if (Lw_control%do_ckd2p1) then 
	call get_totch2obd(4, Optical, totch2o_tmp)
	tmp1(:,:,KS:KE) = diffac*totch2o_tmp(:,:,KS+1:KE+1)
	call get_totch2obd(5, Optical, totch2o_tmp)
	tmp2(:,:,KS:KE) = diffac*totch2o_tmp(:,:,KS+1:KE+1)
	call get_totch2obd(7, Optical, totch2o_tmp)
	tmp3(:,:,KS:KE) = diffac*totch2o_tmp(:,:,KS+1:KE+1)
      else
	tmp1(:,:,KS:KE) = betacm(12)*Optical%totvo2(:,:,KS+1:KE+1)
	tmp2(:,:,KS:KE) = betacm(13)*Optical%totvo2(:,:,KS+1:KE+1)
	tmp3(:,:,KS:KE) = betacm(15)*Optical%totvo2(:,:,KS+1:KE+1)
      endif

      if (Lw_control%do_lwaerosol) then
!	call get_totaerooptdep(4, Optical, totaer_tmp)
	totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,4)
	tmp1(:,:,KS:KE) =  tmp1(:,:,KS:KE) +    &
                           totaer_tmp   (:,:,KS+1:KE+1  )
!	call get_totaerooptdep(5, Optical, totaer_tmp)
	totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,5)
	tmp2(:,:,KS:KE) =  tmp2(:,:,KS:KE) +    &
                           totaer_tmp   (:,:,KS+1:KE+1)
!	call get_totaerooptdep(7, Optical, totaer_tmp)
	totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,7)
	tmp3(:,:,KS:KE) =  tmp3(:,:,KS:KE) +    &
                           totaer_tmp   (:,:,KS+1:KE+1)
      endif

!     cnttaub1(:,:,KS:KE) = EXP(-1.0*tmp1(:,:,KS:KE))
!     cnttaub2(:,:,KS:KE) = EXP(-1.0*tmp2(:,:,KS:KE))
!     cnttaub3(:,:,KS:KE) = EXP(-1.0*tmp3(:,:,KS:KE))
      cnttaub1(:,:,KS) = 1.0                       
      cnttaub2(:,:,KS) = 1.0                       
      cnttaub3(:,:,KS) = 1.0                       
      cnttaub1(:,:,KS+1:KE+1) = EXP(-1.0*tmp1(:,:,KS:KE))
      cnttaub2(:,:,KS+1:KE+1) = EXP(-1.0*tmp2(:,:,KS:KE))
      cnttaub3(:,:,KS+1:KE+1) = EXP(-1.0*tmp3(:,:,KS:KE))

!---------------------------------------------------------------------
!     if cfcs are included, add transmission functions for f11, f12,
!     f113, and f22.
!---------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_exact (4, Optical, cfc_tf)
!       cnttaub1(:,:,KS:KE) = cnttaub1(:,:,KS:KE)*cfc_tf(:,:,KS:KE)
        cnttaub1(:,:,KS+1:KE+1) = cnttaub1(:,:,KS+1:KE+1)*cfc_tf(:,:,KS:KE)
        call cfc_exact (5, Optical, cfc_tf)
!       cnttaub2(:,:,KS:KE) = cnttaub2(:,:,KS:KE)*cfc_tf(:,:,KS:KE)
        cnttaub2(:,:,KS+1:KE+1) = cnttaub2(:,:,KS+1:KE+1)*cfc_tf(:,:,KS:KE)
        call cfc_exact (7, Optical, cfc_tf)
!       cnttaub3(:,:,KS:KE) = cnttaub3(:,:,KS:KE)*cfc_tf(:,:,KS:KE)
        cnttaub3(:,:,KS+1:KE+1) = cnttaub3(:,:,KS+1:KE+1)*cfc_tf(:,:,KS:KE)
      endif 
 
!----------------------------------------------------------------------
!     evaluate h2o (mbar*phibar) between level KS and other levels.
!----------------------------------------------------------------------
      Optical%avephi(:,:,KS:KE) = Optical%totphi(:,:,KS+1:KE+1)
 
!----------------------------------------------------------------------
!     the evaluation of emiss over the layer between data level (KS)
!     and flux level (KE+1) is done by averaging E2 functions referring
!     to the top and bottom of the layer. a special value of (mbar*
!     phibar) is required; it is stored in the (otherwise vacant)
!     KE+1'th position of avephi.
!----------------------------------------------------------------------
      Optical%avephi(:,:,KE+1) = Optical%avephi(:,:,KE-1) + Optical%emx1(:,:)

!----------------------------------------------------------------------
!     if h2o lines in the 1200-1400 range are assumed to have a temp-
!     erature dependent intensity, similar evaluation for (mbar*phibar)
!     is performed, with a special value for the lowest layer
!----------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
         Optical%avephif(:,:,KS:KE,m) = Optical%tphfh2o(:,:,KS+1:KE+1,m)
        enddo
        do m=1,NBTRGE
         Optical%avephif(:,:,KE+1,m) = Optical%avephif(:,:,KE-1,m) +  &
	                Optical%emx1f(:,:,m)
        enddo
      endif
!----------------------------------------------------------------------


end subroutine optical_trans_funct_from_KS




!####################################################################

!subroutine optical_trans_funct_k_down (Gas_tf, k, cnttaub1, cnttaub2, &
!			       cnttaub3, to3cnt, overod, Optical, &
!			               contodb1, contodb2, contodb3)
subroutine optical_trans_funct_k_down (Gas_tf, k,                     &
			                 to3cnt, overod, Optical)   

!----------------------------------------------------------------------
integer,                 intent (in) :: k
!real, dimension(:,:,:),  intent(in)  :: cnttaub1, cnttaub2, cnttaub3
!real, dimension (:,:,:), intent(out) :: to3cnt, overod, contodb1,   &
real, dimension (:,:,:), intent(out) :: to3cnt, overod
!				contodb2, contodb3
type(optical_path_type), intent(inout) :: Optical
type(gas_tf_type), intent(inout) :: Gas_tf 
!---------------------------------------------------------------------

   real, dimension (size(to3cnt,1), size(to3cnt,2), &
                 size(to3cnt,3)) ::    &
                                                avmo3, avpho3, tmp1, &
!				        tmp2, cfc_tf, avvo2,  &
					        tmp2,         avvo2,  &
					        avckdwd, avckdo3, &
	                                        avaero3, totch2o_tmp,  &
					        totaer_tmp, tn2o17, &
	                      		        totaerooptdep_15

   real, dimension (size(to3cnt,1), size(to3cnt,2), &
                 size(to3cnt,3)-1) :: cfc_tf

       integer       :: kp, m
!--------------------------------------------------------------------



       if (Lw_control%do_ckd2p1) then
         call get_totch2obd(6, Optical, totch2o_tmp)
       endif
       if (Lw_control%do_lwaerosol) then
!         call get_totaerooptdep(6, Optical, totaer_tmp)
	totaer_tmp(:,:,:) = Optical%totaerooptdep(:,:,:,6)
       endif

       do kp=1,KE+1-k
         avmo3 (:,:,kp+k-1) = Optical%toto3 (:,:,kp+k) - Optical%toto3 (:,:,k)
         avpho3(:,:,kp+k-1) = Optical%tphio3(:,:,kp+k) - Optical%tphio3(:,:,k) 
         if (.not. Lw_control%do_ckd2p1) then
           avvo2 (:,:,kp+k-1) = Optical%totvo2(:,:,kp+k) - Optical%totvo2(:,:,k)
         else
           avckdwd(:,:,kp+k-1) = Optical%totch2obdwd(:,:,kp+k) -   &
                                 Optical%totch2obdwd(:,:,k)
           avckdo3(:,:,kp+k-1) = totch2o_tmp(:,:,kp+k) -  &
                                 totch2o_tmp(:,:,k)
	 endif			 
         if (Lw_control%do_lwaerosol) then
           avaero3(:,:,kp+k-1) =  &
                         totaer_tmp   (:,:,kp+k) - totaer_tmp   (:,:,k)
	 endif
       enddo

       do kp=1,KE+1-k
!        contodb1(:,:,kp+k-1) = cnttaub1(:,:,kp+k-1)/cnttaub1(:,:,k-1)
!        contodb2(:,:,kp+k-1) = cnttaub2(:,:,kp+k-1)/cnttaub2(:,:,k-1)
!        contodb3(:,:,kp+k-1) = cnttaub3(:,:,kp+k-1)/cnttaub3(:,:,k-1)
         Optical%avephi  (:,:,kp+k-1) = Optical%totphi(:,:,kp+k) -  &
	                                 Optical%totphi(:,:,k)
       end do
       Optical%avephi (:,:,KE+1) = Optical%avephi(:,:,KE-1) + Optical%emx1(:, :)


!---------------------------------------------------------------------
!     if h2o lines in the 1200-1400 range are assumed to have a temp-
!     erature dependent intensity, similar evaluation for (mbar*phibar)
!     is performed, with a special value for the lowest layer
!---------------------------------------------------------------------
       if (Lw_control%do_ch4_n2o) then
         do m=1,NBTRGE
           do kp=1,KE+1-k
             Optical%avephif(:,:,kp+k-1,m) =  Optical%tphfh2o(:,:,kp+k,m) -  &
                                  Optical%tphfh2o(:,:,k,   m)
           enddo
           Optical%avephif(:,:,KE+1,m) = Optical%avephif(:,:,KE-1,m) + &
	           Optical%emx1f(:,:,m)
         enddo
       endif

!----------------------------------------------------------------------
!   compute transmission function in the 560-800 cm-1 range
!   evaluate  optical depth contributions 
!
!   add contributions from h2o(lines) and h2o(continuum).
!   h2o(continuum) contributions are either Roberts or CKD2.1
!----------------------------------------------------------------------
       tmp1     (:,:,k:KE) = SQRT(ab15wd*Optical%avephi(:,:,k:KE)) 
       if (Lw_control%do_ckd2p1) then
	 tmp1(:,:,k:KE) = tmp1(:,:,k:KE) + diffac*   &
                          avckdwd    (:,:,k:KE)
       else
	 tmp1(:,:,k:KE) = tmp1(:,:,k:KE) + betawd*   &
                          avvo2      (:,:,k:KE)
       endif

!-------------------------------------------------------------------
!      add contribution from longwave aerosols (if desired).
!-------------------------------------------------------------------
       if (Lw_control%do_lwaerosol) then
!	 call get_totaerooptdep_15 (Optical, totaerooptdep_15)
	totaerooptdep_15(:,:,:) = Optical%totaerooptdep_15(:,:,:)
	 do kp=k,KE
	   tmp1(:,:,kp) = tmp1(:,:,kp) +    &
                          (totaerooptdep_15(:,:,kp+1) -   &
                           totaerooptdep_15(:,:,k) )
	 end do
       endif

!---------------------------------------------------------------------
!    compute transmission function due to these contributions. the
!    effects of co2, n2o  and  cfc's (not exponentials) are added
!    later.
!--------------------------------------------------------------------
!      overod(:,:,k:KE) = EXP(-1.0E+00*tmp1     (:,:,k:KE))
       overod(:,:,k+1:KE+1) = EXP(-1.0E+00*tmp1     (:,:,k:KE))

!----------------------------------------------------------------------
!       add contribution from the 17 um n2o band (if desired).
!   the expression with tn2o17 retains the 560-630 cm-1 equi-
!   valent widths in evaluating 560-800 cm-1 transmissivities.
!---------------------------------------------------------------------
       if (Lw_control%do_ch4_n2o) then
!         call get_gastf_for_optpath (Gas_tf, k+1, KE+1, tn2o17)
!	   tn2o17_d(:,:,kst:kend) = Gas_tf%tn2o17(:,:,kst:kend)
	   tn2o17(:,:,k+1:ke+1) = Gas_tf%tn2o17(:,:,k+1:ke+1)

         if (NBCO215 .EQ. 2) then
!          overod(:,:,k:KE) = overod(:,:,k:KE) *(130./240. +  &
           overod(:,:,k+1:KE+1) = overod(:,:,k+1:KE+1) *(130./240. +  &
		 	      110./240.*tn2o17(:,:,k+1:KE+1))
         elseif (NBCO215 .EQ. 3) then
!          overod(:,:,k:KE) = overod(:,:,k:KE)*(170./240. + &
           overod(:,:,k+1:KE+1) = overod(:,:,k+1:KE+1)*(170./240. + &
			      70./240.*tn2o17(:,:,k+1:KE+1))
         endif
       endif

!----------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod .
!----------------------------------------------------------------------
       if (Lw_control%do_cfc) then
         call cfc_overod_part ( Optical, cfc_tf, k)
!        overod(:,:,k:KE) = overod(:,:,k:KE)*cfc_tf(:,:,k:KE)
         overod(:,:,k+1:KE+1) = overod(:,:,k+1:KE+1)*cfc_tf(:,:,k:KE)
       endif

!--------------------------------------------------------------------
!   compute transmission functions in 990-1070 cm-1 range, including
!   ozone and h2o continuum, from level k to all other levels. 
!---------------------------------------------------------------------
       tmp1  (:,:,k:KE) = bo3rnd(2)*avpho3(:,:,k:KE)/avmo3(:,:,k:KE)
       tmp2(:,:,k:KE) = 0.5*(tmp1(:,:,k:KE)*(SQRT(1.0E+00 + (4.0E+00* &
                        ao3rnd(2)*avmo3(:,:,k:KE))/tmp1(:,:,k:KE))  &
                        - 1.0E+00))
       if (Lw_control%do_ckd2p1) then
	 tmp2(:,:,k:KE) = tmp2(:,:,k:KE) + diffac*   &
                          avckdo3  (:,:,k:KE) 
       else
	 tmp2(:,:,k:KE) = tmp2(:,:,k:KE) + betacm(14)*   &
                          avvo2 (:,:,k:KE)
       endif
       if (Lw_control%do_lwaerosol) then
	 tmp2(:,:,k:KE) = tmp2(:,:,k:KE) +   &
                          avaero3      (:,:,k:KE)
       endif
!      to3cnt(:,:,k:KE) = EXP(-1.0E+00*tmp2(:,:,k:KE))
       to3cnt(:,:,k+1:KE+1) = EXP(-1.0E+00*tmp2(:,:,k:KE))

!---------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in to3cnt.
!---------------------------------------------------------------------
       if (Lw_control%do_cfc) then
         call cfc_exact_part (6, Optical, cfc_tf, k)
!        to3cnt(:,:,k:KE) = to3cnt(:,:,k:KE)*cfc_tf(:,:,k:KE)
         to3cnt(:,:,k+1:KE+1) = to3cnt(:,:,k+1:KE+1)*cfc_tf(:,:,k:KE)
       endif 
!---------------------------------------------------------------------


end subroutine optical_trans_funct_k_down



!#################################################################

!subroutine optical_trans_funct_KE (Gas_tf, cnttaub1, cnttaub2,    &
!				   cnttaub3, to3cnt, Optical,   &
!				   overod, contodb1, &
!				   contodb2, contodb3)
subroutine optical_trans_funct_KE (Gas_tf,                        &
				             to3cnt, Optical,   &
				   overod)               

!---------------------------------------------------------------------
!real, dimension (:,:,:), intent(in) ::  cnttaub1, cnttaub2,  cnttaub3
real, dimension (:,:,:), intent(out) :: to3cnt, overod
!real, dimension (:,:,:), intent(out) :: to3cnt, overod, contodb1, &
!					  contodb2, contodb3
type(optical_path_type), intent(inout) :: Optical
type(gas_tf_type), intent(inout) :: Gas_tf 

!---------------------------------------------------------------------

   real, dimension (size(to3cnt,1), size(to3cnt,2), &
                 size(to3cnt,3)) ::    &
                                             tmp1, tmp2, tn2o17

   real, dimension (size(to3cnt,1), size(to3cnt,2), &
                 size(to3cnt,3)-1) ::    &
                                          cfc_tf

   real, dimension (size(to3cnt,1), size(to3cnt,2)) :: &
                                             aerooptdep_KE_15

!---------------------------------------------------------------------


!-----------------------------------------------------------------------
!   compute transmission function in the 560-800 cm-1 range
!
!   evaluate  optical depth contributions 
!
!       add contributions from h2o(lines) and h2o(continuum).
!       h2o(continuum) contributions are either Roberts or CKD2.1
!----------------------------------------------------------------------
      tmp1     (:,:,KE) = SQRT(ab15wd*Optical%var2  (:,:,KE)) 
      if (Lw_control%do_ckd2p1) then
	tmp1(:,:,KE) = tmp1(:,:,KE) + diffac*   &
                         Optical%xch2obdwd   (:,:,KE)
      else
	tmp1(:,:,KE) = tmp1(:,:,KE) + betawd*  &
                          Optical%cntval     (:,:,KE)
      endif

!---------------------------------------------------------------------
!      add contribution from longwave aerosols (if desired).
!---------------------------------------------------------------------
      if (Lw_control%do_lwaerosol) then
!       call get_aerooptdep_KE_15 (Optical, aerooptdep_KE_15)
 	aerooptdep_KE_15(:,:) = Optical%aerooptdep_KE_15(:,:)
	tmp1(:,:,KE) = tmp1(:,:,KE) +     &
                          aerooptdep_KE_15(:,:)  
      endif
 
!---------------------------------------------------------------------
!    compute transmission function due to these contributions. the
!    effects of co2, n2o  and  cfc's (not exponentials) are added
!    later.
!---------------------------------------------------------------------
!     overod(:,:,KE) = EXP(-1.0E+00*tmp1     (:,:,KE))
      overod(:,:,KE+1) = EXP(-1.0E+00*tmp1     (:,:,KE))
 
!---------------------------------------------------------------------
!       add contribution from the 17 um n2o band (if desired).
!   the expression with tn2o17 retains the 560-630 cm-1 equi-
!   valent widths in evaluating 560-800 cm-1 transmissivities.
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
!       call get_gastf_for_optpath (Gas_tf, KE+1, KE+1, tn2o17)
!	  tn2o17_d(:,:,kst:kend) = Gas_tf%tn2o17(:,:,kst:kend)
	  tn2o17(:,:,ke+1    ) = Gas_tf%tn2o17(:,:,ke+1)
        if (NBCO215 .EQ. 2) then
!         overod(:,:,KE) = overod(:,:,KE) *  &
          overod(:,:,KE+1) = overod(:,:,KE+1) *  &
                           (130./240. + 110./240.*tn2o17(:,:,KE+1))
        elseif (NBCO215 .EQ. 3) then
!         overod(:,:,KE) = overod(:,:,KE) *   &
          overod(:,:,KE+1) = overod(:,:,KE+1) *   &
                           (170./240. + 70./240.*tn2o17(:,:,KE+1))
        endif
      endif

!---------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod .
!---------------------------------------------------------------------
     if (Lw_control%do_cfc) then
       call cfc_overod_part (Optical, cfc_tf, KE)
!      overod(:,:,KE) = overod(:,:,KE)*cfc_tf(:,:,KE)
       overod(:,:,KE+1) = overod(:,:,KE+1)*cfc_tf(:,:,KE)
     endif 

!-----------------------------------------------------------------------
!   compute transmission functions in 990-1070 cm-1 range, including
!   ozone and h2o continuum, from level KS to all other levels. 
!---------------------------------------------------------------------
     tmp1  (:,:,KE) = bo3rnd(2)*Optical%var4(:,:,KE)/Optical%var3(:,:,KE)
     tmp2(:,:,KE) = 0.5*(tmp1(:,:,KE)*(SQRT(1.0E+00 + (4.0E+00*  &
                    ao3rnd(2)*Optical%var3 (:,:,KE))/tmp1(:,:,KE))  - 1.0E+00))

     if (Lw_control%do_ckd2p1) then
       tmp2(:,:,KE) = tmp2(:,:,KE) + diffac*Optical%xch2obd  (:,:,KE,6) 
     else
       tmp2(:,:,KE) = tmp2(:,:,KE) + betacm(14)*Optical%cntval (:,:,KE)
     endif

!    to3cnt(:,:,KE) = EXP(-1.0E+00*tmp2(:,:,KE))
     to3cnt(:,:,KE+1) = EXP(-1.0E+00*tmp2(:,:,KE))

!---------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod and to3cnt.
!---------------------------------------------------------------------
     if (Lw_control%do_cfc) then
       call cfc_exact_part (6, Optical, cfc_tf, KE)
!      to3cnt(:,:,KE) = to3cnt(:,:,KE)*cfc_tf(:,:,KE)
       to3cnt(:,:,KE+1) = to3cnt(:,:,KE+1)*cfc_tf(:,:,KE)
     endif

!    contodb1(:,:,KE  ) = cnttaub1(:,:,KE)/cnttaub1(:,:,KE-1)
!    contodb2(:,:,KE  ) = cnttaub2(:,:,KE)/cnttaub2(:,:,KE-1)
!    contodb3(:,:,KE  ) = cnttaub3(:,:,KE)/cnttaub3(:,:,KE-1)

!-------------------------------------------------------------------


end subroutine optical_trans_funct_KE




!####################################################################

!subroutine optical_trans_funct_diag ( pflux, press, contdg, to3dg, &
subroutine optical_trans_funct_diag ( Atmos_input, contdg, to3dg, &
                                      Optical)

!---------------------------------------------------------------------
!real, dimension (:,:,:), intent(in)     :: pflux, press
real, dimension (:,:,:), intent(out)    :: to3dg                
real, dimension (:,:,:,:), intent(out)  :: contdg               
type(optical_path_type), intent(inout) :: Optical
type(atmos_input_type), intent(in) :: Atmos_input
					
!---------------------------------------------------------------------
!   local variables:
!---------------------------------------------------------------------
     real, dimension (size(Atmos_input%pflux,1), size(Atmos_input%pflux,2), &
                         size(Atmos_input%pflux,3)-1) :: pdfinv

     real, dimension (size(Atmos_input%pflux,1), size(Atmos_input%pflux,2), &
                         size(Atmos_input%pflux,3)) ::  &
                          press, pflux, &
                   ca, cb, csuba, csubb, ctmp2, ctmp3, delpr1, delpr2

!  convert press and pflux to cgs.
      press(:,:,:) = 10.0*Atmos_input%press(:,:,:)
       pflux(:,:,:) = 10.0*Atmos_input%pflux(:,:,:)



     pdfinv(:,:,ks:ke) = 1.0/(pflux(:,:,ks+1:ke+1) - pflux(:,:,ks:ke))

!----------------------------------------------------------------------


!---------------------------------------------------------------------
    delpr1(:,:,KS+1:KE)   = pdfinv (:,:,KS+1:KE)*(press(:,:,KS+1:KE) -&
                              pflux(:,:,KS+1:KE)) 
    delpr2(:,:,KS+1:KE+1) = pdfinv(:,:,KS:KE)*   &
                              (pflux(:,:,KS+1:KE+1) - press(:,:,KS:KE)) 

!-----------------------------------------------------------------------
!     compute nearby-layer transmissivities for the o3 band and for the
!     one-band continuum band.  the sf function is used.
!     the method is the same as described for co2 in reference(4).
!-----------------------------------------------------------------------
    if ( .not. Lw_control%do_ckd2p1) then
      ctmp2(:,:,KS+1:KE)  = Optical%cntval(:,:,KS+1:KE)*delpr1(:,:,KS+1:KE) 
      ctmp3(:,:,KS+1:KE)  = Optical%cntval(:,:,KS:KE-1)*delpr2(:,:,KS+1:KE) 
    endif
    
!-----------------------------------------------------------------------
!    compute sf2.
!    continuum band 1
!-----------------------------------------------------------------------
    if ( .not. Lw_control%do_ckd2p1) then
      csuba(:,:,KS+1:KE)  = betacm(12)*ctmp2(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = betacm(12)*ctmp3(:,:,KS+1:KE)
	else
      csuba(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS+1:KE,4)*  &
                            delpr1(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS:KE-1,4)*  &
                            delpr2(:,:,KS+1:KE)
    endif
    ca    (:,:,KS+1:KE) = csuba(:,:,KS+1:KE)*(-0.5E+00 +    &
                          csuba(:,:,KS+1:KE)*(0.166666E+00 -  &
                          csuba(:,:,KS+1:KE)*0.416666E-01))   
    cb    (:,:,KS+1:KE) = csubb(:,:,KS+1:KE)*(-0.5E+00 +  &
                          csubb(:,:,KS+1:KE)*(0.166666E+00 - &
                          csubb(:,:,KS+1:KE)*0.416666E-01)) 
    contdg(:,:,KE+1,1)    = 1.0E+00 + cb (:,:,KE)
    contdg(:,:,KS+1:KE,1) = 1.0E+00 + 0.5E+00*(ca (:,:,KS+1:KE) +  &
                            cb (:,:,KS+1:KE))
!--------------------------------------------------------------------
!    continuum band 2
!---------------------------------------------------------------------
    if ( .not. Lw_control%do_ckd2p1) then
      csuba(:,:,KS+1:KE)  = betacm(13)*ctmp2(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = betacm(13)*ctmp3(:,:,KS+1:KE)
    else
      csuba(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS+1:KE,5)*   &
                            delpr1(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS:KE-1,5)*  &
                            delpr2(:,:,KS+1:KE)
    endif
    ca    (:,:,KS+1:KE) = csuba(:,:,KS+1:KE)*(-0.5E+00 +  &
                          csuba(:,:,KS+1:KE)*(0.166666E+00 -   &
                          csuba(:,:,KS+1:KE)*0.416666E-01)) 
    cb    (:,:,KS+1:KE) = csubb(:,:,KS+1:KE)*(-0.5E+00 +   &
                          csubb(:,:,KS+1:KE)*(0.166666E+00 -   &
                          csubb(:,:,KS+1:KE)*0.416666E-01)) 
    contdg(:,:,KE+1,2)    = 1.0E+00 + cb (:,:,KE)
    contdg(:,:,KS+1:KE,2) = 1.0E+00 + 0.5E+00*(ca (:,:,KS+1:KE) +  &
                            cb (:,:,KS+1:KE))
!--------------------------------------------------------------------
!    continuum band 3
!--------------------------------------------------------------------
    if ( .not. Lw_control%do_ckd2p1) then
      csuba(:,:,KS+1:KE)  = betacm(15)*ctmp2(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = betacm(15)*ctmp3(:,:,KS+1:KE)
    else
      csuba(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS+1:KE,7)*   &
                            delpr1(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS:KE-1,7)*  &
                            delpr2(:,:,KS+1:KE)
    endif
    ca    (:,:,KS+1:KE) = csuba(:,:,KS+1:KE)*(-0.5E+00 +    &
                          csuba(:,:,KS+1:KE)*(0.166666E+00 -  &
                          csuba(:,:,KS+1:KE)*0.416666E-01)) 
    cb    (:,:,KS+1:KE) = csubb(:,:,KS+1:KE)*(-0.5E+00 +   &
                          csubb(:,:,KS+1:KE)*(0.166666E+00 -  &
                          csubb(:,:,KS+1:KE)*0.416666E-01)) 
    contdg(:,:,KE+1,3)    = 1.0E+00 + cb (:,:,KE)
    contdg(:,:,KS+1:KE,3) = 1.0E+00 + 0.5E+00*(ca (:,:,KS+1:KE) +   &
                            cb (:,:,KS+1:KE))
!--------------------------------------------------------------------
!    ozone band
!--------------------------------------------------------------------
    if ( .not. Lw_control%do_ckd2p1) then
      csuba(:,:,KS+1:KE)  = betacm(14)*ctmp2(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = betacm(14)*ctmp3(:,:,KS+1:KE)
    else
      csuba(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS+1:KE,6)*   &
                            delpr1(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = diffac*Optical%xch2obd(:,:,KS:KE-1,6)*  &
                            delpr2(:,:,KS+1:KE)
    endif
    ca   (:,:,KS+1:KE)  = csuba(:,:,KS+1:KE)*(-0.5E+00 +   &
                          csuba(:,:,KS+1:KE)*   &
                          (0.166666E+00 - csuba(:,:,KS+1:KE)*  &
                          0.416666E-01)) 
    cb   (:,:,KS+1:KE)  = csubb(:,:,KS+1:KE)*(-0.5E+00 +  &
                          csubb(:,:,KS+1:KE)*   &
                          (0.166666E+00 - csubb(:,:,KS+1:KE)*   &
                          0.416666E-01)) 
    to3dg (:,:,KE+1)    = 1.0E+00 + cb(:,:,KE)
    to3dg (:,:,KS+1:KE) = 1.0E+00 + 0.5E+00*(ca(:,:,KS+1:KE) +   &
                           cb(:,:,KS+1:KE))
!-------------------------------------------------------------------



end subroutine optical_trans_funct_diag


!###################################################################

subroutine get_totch2o (n, Optical, totch2o, dte1, ixoe1)

!------------------------------------------------------------------
real, dimension(:,:,:), intent(in)     :: dte1    
type(optical_path_type), intent(inout) :: Optical
integer, dimension(:,:,:), intent(in)     :: ixoe1   
real, dimension(:,:,:), intent(out)     :: totch2o
integer,                intent(in)      :: n

!-----------------------------------------------------------------

         real, dimension(size(Optical%tfac,1), size(Optical%tfac,2), &
	             size(Optical%tfac,3)) ::  radf, sh2o , tmpexp

         real               ::  t0 = 296.0
         integer            ::  k, nu
         real               ::  fh2o0, sh2o0

!--------------------------------------------------------------------
!     compute self-broadened temperature-dependent continuum coefficient
!     using the single coefficient -.013 for all frequencies in
!     the 160-560 cm-1 range. experiments with the mid-latitude
!     summer profile show errors of < .01 W/m**2 (in the net broadband
!     flux, 0-2200 cm-1) using this value. this value is used instead
!     of tmpfctrs at each frequency band.
!--------------------------------------------------------------------
!         Optical%tmpexp(:,:,KS:KE) = EXP(-.013*Optical%tfac(:,:,KS:KE))
                 tmpexp(:,:,KS:KE) = EXP(-.013*Optical%tfac(:,:,KS:KE))

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!     compute source function for frequency bands (ioffh2o+1 to ioffh2o
!     +nptch2o) at layer temperatures using table lookup.
!     note that ixoe1 can be used for temp index, and dte1 for deltat,
!     as the table extent for radf is the same as for the e1 tables
!     of the model.
!--------------------------------------------------------------------
         nu = n
         call looktab (radfunc, ixoe1, dte1, radf(:,:,:), KS, KE,   &
		       nu+ioffh2o)
         sh2o0 = ssh2o_296(nu+ioffh2o)*sfac(nu+ioffh2o)

         do k=KS,KE 
! 	   Optical%sh2o(:,:,k) = sh2o0*Optical%tmpexp(:,:,k)
 	           sh2o(:,:,k) = sh2o0*        tmpexp(:,:,k)
         enddo
 
!--------------------------------------------------------------------
!     compute h2o self- and foreign- broadened continuum optical path,
!     summed from the top of the atmosphere through layer k.
!--------------------------------------------------------------------
         fh2o0 = sfh2o(nu+ioffh2o)*fscal(nu+ioffh2o)
         totch2o(:,:,1) = 0.0E+00
	 do k = KS+1,KE+1
	   totch2o(:,:,k) = Optical%wk(:,:,k-1)*1.0e-20*   &
!                            (Optical%sh2o(:,:,k-1)*Optical%rh2os(:,:,k-1) +    &
                            (        sh2o(:,:,k-1)*Optical%rh2os(:,:,k-1) +    &
	                    fh2o0*Optical%rfrgn(:,:,k-1))* &
                            vvj(nu)*radf(:,:,k-1   )    +   &
                            totch2o(:,:,k-1)
         enddo

!------------------------------------------------------------------

!------------------------------------------------------------------

end subroutine get_totch2o



!#####################################################################


subroutine get_totch2obd (n, Optical, totch2obd)

!------------------------------------------------------------------
real, dimension(:,:,:), intent(out)     :: totch2obd
integer,                intent(in)      :: n
type(optical_path_type), intent(inout) :: Optical

!-----------------------------------------------------------------
        real               ::  t0 = 296.0
        integer            ::  k, nu
        real               ::  fh2o0, sh2o0

!---------------------------------------------------------------------
!     compute h2o self- and foreign- broadened continuum optical path 
!     for each layer k (xch2obd, xch2obdinw, xch2obdwd) and summed from
!     the top of the atmosphere through layer k (totch2obd,
!     totch2obdinw, totch2obdwd).
!---------------------------------------------------------------------
        nu = n     
	totch2obd(:,:,1) = 0.0E+00
	do k = KS+1,KE+1
	  totch2obd(:,:,k) = totch2obd(:,:,k-1) +   &
                             Optical%xch2obd(:,:,k-1,nu)
	enddo
!--------------------------------------------------------------------
 
end subroutine get_totch2obd




!#####################################################################

subroutine get_totvo2 (n, Optical, totvo2_out) 

!------------------------------------------------------------------
integer,                intent(in)         :: n
type(optical_path_type), intent(inout) :: Optical
real, dimension(:,:,:), intent(out)        :: totvo2_out

!-----------------------------------------------------------------

        totvo2_out(:,:,:) = betacm(n)*Optical%totvo2(:,:,KS+1:KE+1)

end subroutine get_totvo2 


!###################################################################

subroutine optical_ckd2p1_init

!-----------------------------------------------------------------------
!     Idckdh2o reads ckd2.1 self and foreign-broadened h2o continuum
!     coefficients, corrections, and coefficients for temperature
!     dependence of the self-continuum. these are tabulated at 10
!     cm-1 intervals from 0 to 20000 cm-1 (this is subject to change;
!     the above information is as of 2/12/96).
!
!     references:
!
!     (1) clough, s. a.  et al. "line shape and the water vapor
!         continuum," atmospheric research, 23 (1989) 229-241.
!
!
!     author: m. d. schwarzkopf
!-----------------------------------------------------------------------
        integer  :: inrad, k, j, ihih2o, iloh2o, nu
        real     :: fh2o0, sh2o0

!-------------------------------------------------------------------
!   data from former block data bs260 for self-broadened continuum
!    at 260K, band-integrated, in 5 - 19995 cm-1 range.
!               06/28/82
!               units of (cm**3/mol) * 1.E-20
!---------------------------------------------------------------------
        real    ::  v1sh2o_260, v2sh2o_260, dvsh2o_260,    &
                    ssh2o_260(2000)
        integer ::  nptsh2o_260

!--------------------------------------------------------------------
!        tktab and vjtab are the respective temperature and frequency
!    points at which tabulations occurred.
!---------------------------------------------------------------------
        real   ::   tktab(40),  vjtab(300)
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!    call routine to allocate radfunc table
!---------------------------------------------------------------------
        call table_alloc (radfunc, 40, 300)

!--------------------------------------------------------------------
!    read h2o (original) data
!    data are at frequencies 5 - 19995 cm-1, at 10 cm-1 intervals
!-------------------------------------------------------------------
        inrad = open_file ('INPUT/h2ockd2.1_data', form='formatted', &
			    action='read')
        read (inrad,9001) v1sh2o_296, v2sh2o_296, dvsh2o_296,  &
                          nptsh2o_296
	read (inrad,9002) (ssh2o_296(k),k=1,2000)
        read (inrad,9001) v1sh2o_260, v2sh2o_260, dvsh2o_260,   &
                          nptsh2o_260
	read (inrad,9002) (ssh2o_260(k),k=1,2000)
        read (inrad,9001) v1fh2o, v2fh2o, dvfh2o, nptfh2o
	read (inrad,9002) (sfh2o(k),k=1,2000)
9001    format (3f12.1,i8)
9002    format (5e14.5)
 
        call close_file (inrad)

!--------------------------------------------------------------------
!    read h2o corrected data
!--------------------------------------------------------------------
        inrad = open_file ('INPUT/h2ockd2.1_corrdata',   &
			    form='formatted', action='read')
        read (inrad,9007) (sfac(k),k=1,2000)
        read (inrad,9007) (fscal(k),k=1,2000)
        read (inrad,9007) (tmpfctrs(k),k=1,2000)
9007    format (5e13.6)
 
	call close_file (inrad)

!--------------------------------------------------------------------
!    read radfn data
!--------------------------------------------------------------------
        inrad = open_file ('INPUT/radfn_5-2995_100-490k',   &
			    form='formatted', action='read')
        read (inrad,9000) ((radfunc%vae(k,j),radfunc%td(k,j),k=1,40), &
                                                           j=1,300)
9000    format (8f14.6)
        call close_file (inrad)

!---------------------------------------------------------------------
        do k=1,40
  	  tktab(k) = 100. + 10.*(k-1)
        enddo
        do j=1,300
	  vjtab(j) = 5. + 10.*(j-1)
        enddo
 
!--------------------------------------------------------------------
!    compute range to use in datasets for actual frequency intervals
!    used in model.
!
!     freqlo = 160.
!     freqhi = 560.
!
!     define initial offset and number of data points to use
!     for the 3 h2o continua over the frequency range of the
!     calculations (freqlo,freqhi). note: we assume no interpolation
!     is needed. if interp. was required, these limits would be
!     expanded. values are put into commons in include file tab.h
!     for transmission into Optical_ckd2.1.F.
!
!     ioff is the offset from the absorption tables (starting at 5)
!     needed for proper freq computations. first index used then
!     is (ioff+1). for calculations with the first band beginning
!     at 160 cm-1, this number is 16, and the index number of the
!     band ending at 560 cm-1 is 56.
!-----------------------------------------------------------------------
        ioffh2o = 16

!---------------------------------------------------------------------
!   the final index number used in the calculation is (ihi)
!--------------------------------------------------------------------
        ihih2o  = 56

!--------------------------------------------------------------------
!   nptc is the number of frequency points used in the calculation.
!    ( = ihi - (ioff+1) + 1)
!---------------------------------------------------------------------
        nptch2o = ihih2o - ioffh2o

!---------------------------------------------------------------------
!   vvj are the frequencies for calculation of h2o coefficients. by
!   assumption, no other frequencies are used.
!----------------------------------------------------------------------
        do j=1,nptch2o
	  vvj(j) = v1sh2o_296 + dvsh2o_296*float(j+ioffh2o-1)
	enddo

!---------------------------------------------------------------------
!     compute h2o coefficients averaged over the broad bands used
!     in the 560 -1200 cm-1 range. until the frequency bands are read
!     in, I will re-list them here, rather than use rnddta.H variables
!     (where they are stored).
!
!     the required wide bands are:
!        560-630 cm-1
!	 630-700   (assuming 3 bands in 15um complex)
!	 700-800
!	 560-800   (1 band for entire complex)
!	 800-900
!	 900-990
!	 990-1070
!	 1070-1200
!	 800-900,1070-1200   (until this band is broken into 2)
!        we assume that, for best accuracy:
!     the quantity required is <svj> and <fvj) where angle brackets are
!     averages over frequency, s and f are self- and foreign coeff-
!     icients, including corrections, and vj is frequency (from vjtab).
!     notations for special bands attempt similarity with that
!     previously used in the radiation code.
!        we also assume that one value may be used (at all altitudes)
!     for the radiation correction term radfn, in each frequency band.
!     the values used below result from experimentation.
!---------------------------------------------------------------------
        svj = 0.0
        fvj = 0.0
        svjwd = 0.0
        fvjwd = 0.0
        svjinw = 0.0
        fvjinw = 0.0
!--------------------------------------------------------------------
!        560-630 band:
!--------------------------------------------------------------------
        do j=57,63
	  svj(1) = svj(1) + vjtab(j)*ssh2o_296(j)*sfac(j)/7.
	  fvj(1) = fvj(1) + vjtab(j)*sfh2o(j)*fscal(j)/7.
        enddo
        radfnbd(1) = 0.90
!--------------------------------------------------------------------
!        630-700 band:
!--------------------------------------------------------------------
        do j=64,70
	  svj(2) = svj(2) + vjtab(j)*ssh2o_296(j)*sfac(j)/7.
	  fvj(2) = fvj(2) + vjtab(j)*sfh2o(j)*fscal(j)/7.
        enddo
        radfnbd(2) = 0.92
!--------------------------------------------------------------------
!        700-800 band:
!--------------------------------------------------------------------
        do j=71,80
	  svj(3) = svj(3) + vjtab(j)*ssh2o_296(j)*sfac(j)/10.
	  fvj(3) = fvj(3) + vjtab(j)*sfh2o(j)*fscal(j)/10.
        enddo
        radfnbd(3) = 0.95
!--------------------------------------------------------------------
!        800-900 band:
!--------------------------------------------------------------------
        do j=81,90
	  svj(4) = svj(4) + vjtab(j)*ssh2o_296(j)*sfac(j)/10.
	  fvj(4) = fvj(4) + vjtab(j)*sfh2o(j)*fscal(j)/10.
        enddo
        radfnbd(4) = 0.97
!--------------------------------------------------------------------
!        900-990 band:
!--------------------------------------------------------------------
        do j=91,99
	  svj(5) = svj(5) + vjtab(j)*ssh2o_296(j)*sfac(j)/9.
	  fvj(5) = fvj(5) + vjtab(j)*sfh2o(j)*fscal(j)/9.
        enddo
        radfnbd(5) = 0.98
!--------------------------------------------------------------------
!        990-1070 band:
!--------------------------------------------------------------------
        do j=100,107
	  svj(6) = svj(6) + vjtab(j)*ssh2o_296(j)*sfac(j)/8.
	  fvj(6) = fvj(6) + vjtab(j)*sfh2o(j)*fscal(j)/8.
        enddo
        radfnbd(6) = 0.99
!--------------------------------------------------------------------
!        1070-1200 band:
!--------------------------------------------------------------------
        do j=108,120
	  svj(7) = svj(7) + vjtab(j)*ssh2o_296(j)*sfac(j)/13.
	  fvj(7) = fvj(7) + vjtab(j)*sfh2o(j)*fscal(j)/13.
        enddo
        radfnbd(7) = 0.992
!--------------------------------------------------------------------
!        560-800 combined band:
!--------------------------------------------------------------------
        do j=57,80
	  svjwd = svjwd + vjtab(j)*ssh2o_296(j)*sfac(j)/24.
	  fvjwd = fvjwd + vjtab(j)*sfh2o(j)*fscal(j)/24.
        enddo
        radfnbdwd = 0.92
!--------------------------------------------------------------------
!        800-990,1070-1200 combined band:
!--------------------------------------------------------------------
        do j=81,99
	  svjinw = svjinw + vjtab(j)*ssh2o_296(j)*sfac(j)/22.
	  fvjinw = fvjinw + vjtab(j)*sfh2o(j)*fscal(j)/32.
        enddo
        do j=108,120
	  svjinw = svjinw + vjtab(j)*ssh2o_296(j)*sfac(j)/32.
	  fvjinw = fvjinw + vjtab(j)*sfh2o(j)*fscal(j)/32.
        enddo
        radfnbdinw = 0.98
!--------------------------------------------------------------------


end subroutine optical_ckd2p1_init




!###################################################################

subroutine optical_path_ckd2p1 (atmden, press, temp, rh2o, Optical) 

!---------------------------------------------------------------------
!    subroutine optical_ckd2p1 computes h2o continuum optical paths
!    (self + foreign) over the frequency range specified by
!    ioffh2o and nptch2o using the ckd2.1 algorithm, modified for 
!    the gcm parameterization.
!    (this routine is previously called contnm.F)
!---------------------------------------------------------------------
real, dimension (:,:,:), intent(in)       ::  atmden, press, temp, rh2o
type(optical_path_type), intent(inout) :: Optical

!---------------------------------------------------------------------
      real, dimension(size(press,1), size(press,2), &
              size(press,3)) :: totch2obdinw

      real, dimension(size(press,1), size(press,2), &
              size(press,3)-1) :: xch2obdinw, tmpexp, rvh2o, rhoave

      real                                ::  t0 = 296.0
      integer                             ::  n, k, nu
      real                                ::  fh2o0, sh2o0

!----------------------------------------------------------------------

      allocate (Optical%xch2obd    (ISRAD:IERAD, JSRAD:JERAD,   KS:KE  , 7) )
      allocate (Optical%totch2obdwd(ISRAD:IERAD, JSRAD:JERAD,   KS:KE+1   ) )
      allocate (Optical%xch2obdwd  (ISRAD:IERAD, JSRAD:JERAD,   KS:KE     ) )

!--------------------------------------------------------------------
!    define the volume mixing ratio of h2o
!---------------------------------------------------------------------
!      Optical%rvh2o(:,:,KS:KE) = rh2o(:,:,KS:KE)*wtmair/wtmh2o
              rvh2o(:,:,KS:KE) = rh2o(:,:,KS:KE)*wtmair/wtmh2o

!---------------------------------------------------------------------
!    define input arguments to optical_ckd2p1
!-------------------------------------------------------------------
!      Optical%wk(:,:,KS:KE) = Optical%rvh2o(:,:,KS:KE)*avogno/wtmair*   &
      Optical%wk(:,:,KS:KE) =         rvh2o(:,:,KS:KE)*avogno/wtmair*   &
                          atmden(:,:,KS:KE)
!      Optical%rhoave(:,:,KS:KE) = (press(:,:,KS:KE)/pstd)*(t0/temp(:,:,KS:KE))
              rhoave(:,:,KS:KE) = (press(:,:,KS:KE)/pstd)*(t0/temp(:,:,KS:KE))
!      Optical%rh2os(:,:,KS:KE) = Optical%rvh2o(:,:,KS:KE)*Optical%rhoave(:,:,KS:KE)
!      Optical%rh2os(:,:,KS:KE) =         rvh2o(:,:,KS:KE)*Optical%rhoave(:,:,KS:KE)
      Optical%rh2os(:,:,KS:KE) =         rvh2o(:,:,KS:KE)*        rhoave(:,:,KS:KE)
!      Optical%rfrgn(:,:,KS:KE) = Optical%rhoave(:,:,KS:KE) - Optical%rh2os(:,:,KS:KE)
      Optical%rfrgn(:,:,KS:KE) =         rhoave(:,:,KS:KE) - Optical%rh2os(:,:,KS:KE)
      Optical%tfac(:,:,KS:KE) = temp(:,:,KS:KE) - t0

!--------------------------------------------------------------------
!     compute self-broadened temperature-dependent continuum coefficient
!     using the single coefficient -.020 for all frequencies in
!     the 560-1200 cm-1 range. experiments with the mid-latitude
!     summer profile show errors of < .01 W/m**2 (in the net broadband
!     flux, 0-2200 cm-1) using this value. this value is used instead
!     of tmpfctrs at each frequency band.
!-------------------------------------------------------------------
!      Optical%tmpexp(:,:,KS:KE) = EXP(-.020*Optical%tfac(:,:,KS:KE))
              tmpexp(:,:,KS:KE) = EXP(-.020*Optical%tfac(:,:,KS:KE))
 
!-------------------------------------------------------------------
!     compute h2o self- and foreign- broadened continuum optical path 
!     for each layer k (xch2obd, xch2obdinw, xch2obdwd) and summed from
!     the top of the atmosphere through layer k (totch2obd,
!     totch2obdinw, totch2obdwd).
!--------------------------------------------------------------------
      do nu = 1,7
        do k = KS,KE 
      Optical%xch2obd(:,:,k,nu) = Optical%wk(:,:,k)*1.0e-20*   &
            (svj(nu)*Optical%rh2os(:,:,k)*&
!                Optical%tmpexp(:,:,k) + fvj(nu)*Optical%rfrgn(:,:,k))* &
                        tmpexp(:,:,k) + fvj(nu)*Optical%rfrgn(:,:,k))* &
                              radfnbd(nu)
        enddo
      enddo
 
      do k = KS,KE 
        xch2obdinw(:,:,k) = Optical%wk(:,:,k)*1.0e-20*(svjinw*  &
	       Optical%rh2os(:,:,k)* &
!	    Optical%tmpexp(:,:,k) + fvjinw*Optical%rfrgn(:,:,k))* &
	            tmpexp(:,:,k) + fvjinw*Optical%rfrgn(:,:,k))* &
                            radfnbdinw
        Optical%xch2obdwd(:,:,k) = Optical%wk(:,:,k)*1.0e-20*(svjwd* &
	                   Optical%rh2os(:,:,k)* &
!	   Optical%tmpexp(:,:,k) + fvjwd*Optical%rfrgn(:,:,k))*  &
	           tmpexp(:,:,k) + fvjwd*Optical%rfrgn(:,:,k))*  &
                           radfnbdwd
      enddo
 
      totch2obdinw(:,:,1) = 0.0E+00
      Optical%totch2obdwd(:,:,1) = 0.0E+00
      do k = KS+1,KE+1
	totch2obdinw(:,:,k) = totch2obdinw(:,:,k-1) +    &
			      xch2obdinw(:,:,k-1)
	Optical%totch2obdwd(:,:,k) = Optical%totch2obdwd(:,:,k-1) + &
	                             Optical%xch2obdwd(:,:,k-1)
      enddo

!----------------------------------------------------------------------

 
end subroutine optical_path_ckd2p1
 


!################################################################## 

subroutine optical_o3 (atmden, qo3, vv, Optical)

!----------------------------------------------------------------------
!     optical_o3 computes optical paths for o3.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!----------------------------------------------------------------------
!
!     intent in:
!
!     qo3     =  mass mixing ratio of o3 at model data levels.
!
!-------------------------------------------------------------------    
real, dimension(:,:,:), intent(in) ::  atmden, qo3, vv
type(optical_path_type), intent(inout) :: Optical

!---------------------------------------------------------------------
!
!     intent out:
!
!       toto3  =  summed o3 optical path from the top of atmosphere to 
!                 flux level k.
!
!       tphio3 =  summed pressure weighted o3 optical path from top of
!                 atmosphere to flux level k.
! 
!       var3   =  o3 optical path in model layers.
! 
!       var4   =  pressure-weighted o3 optical path in model layers.
!  
!---------------------------------------------------------------------
      integer                  ::    k

!--------------------------------------------------------------------

!--------------------------------------------------------------------
      allocate (Optical%toto3 (ISRAD:IERAD, JSRAD:JERAD, KS:KE      +1))
      allocate (Optical%tphio3(ISRAD:IERAD, JSRAD:JERAD, KS:KE      +1))
      allocate (Optical%var3  (ISRAD:IERAD, JSRAD:JERAD, KS:KE        ))
      allocate (Optical%var4  (ISRAD:IERAD, JSRAD:JERAD, KS:KE        ))

!-----------------------------------------------------------------------
!     compute optical paths for o3, using the diffusivity 
!     approximation 1.66 for the angular integration.  obtain 
!     unweighted values var3 and weighted values  var4.
!     the quantities  0.003 (.003) appearing in the
!     var4 expression are the approximate voigt corrections
!     for o3.
!---------------------------------------------------------------------  
      Optical%var3(:,:,KS:KE) = atmden(:,:,KS:KE)*qo3(:,:,KS:KE)*diffac
      Optical%var4(:,:,KS:KE) = Optical%var3(:,:,KS:KE)*(vv(:,:,KS:KE) + 3.0E-03)

!----------------------------------------------------------------------
!     compute summed optical paths for o3.
!----------------------------------------------------------------------
      Optical%toto3 (:,:,KS) = 0.0E+00
      Optical%tphio3(:,:,KS) = 0.0E+00
      do k=KS+1,KE+1
        Optical%toto3 (:,:,k) = Optical%toto3 (:,:,k-1) + Optical%var3  (:,:,k-1) 
        Optical%tphio3(:,:,k) = Optical%tphio3(:,:,k-1) + Optical%var4  (:,:,k-1) 
      enddo
!----------------------------------------------------------------------


end  subroutine optical_o3




!#####################################################################

subroutine optical_rbts (temp, rh2o, Optical) 

real, dimension(:,:,:), intent(in) :: temp, rh2o
type(optical_path_type), intent(inout) :: Optical

!----------------------------------------------------------------------
!     optical_rbts computes optical paths for h2o rbts comtinuum.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!
!-------------------------------------------------------------------    
!
!       cntval =  h2o continuum path in model layers for the 800-990
!                 and 1070-1200 cm-1 combined band.
!
!       totvo2 =  summed h2o continuum path from top of atmosphere to
!                 flux level k.
!-----------------------------------------------------------------------


      real, dimension(size(temp,1), size(temp,2), &
                                     size(temp,3)) :: texpsl

      integer     :: k, i

!---------------------------------------------------------------------

      allocate (Optical%cntval(ISRAD:IERAD, JSRAD:JERAD,   KS:KE+1     ))
      allocate (Optical%totvo2(ISRAD:IERAD, JSRAD:JERAD,   KS:KE+1     ))

!----------------------------------------------------------------------
!     compute argument for constant temperature coefficient (this is 
!     1.800E+03/(1.0E+00/temp - 1.0E+00/2.960E+02)).
!---------------------------------------------------------------------- 
      texpsl(:,:,KS:KE+1)=EXP(1.800E+03/temp(:,:,KS:KE+1) -   &
                              6.081081081E+00) 

!----------------------------------------------------------------------
!     compute optical path for the h2o continuum, using roberts 
!     coefficients betinw, and temperature correction texpsl. 
!     the diffusivity approximation (which cancels out in this
!     expression) is assumed to be 1.66.  the use of the diffusivity
!     factor has been shown to be a significant source of error in the
!     continuum calculations, however, the time penalty of an angular
!     integration is severe.
!---------------------------------------------------------------------  
      Optical%cntval(:,:,KS:KE) = texpsl(:,:,KS:KE)*rh2o(:,:,KS:KE)*   &
                          Optical%var2(:,:,KS:KE)/(rh2o(:,:,KS:KE) +  &
                          rh2oair)

!----------------------------------------------------------------------
!     compute summed optical paths for h2o roberts continuum.
!----------------------------------------------------------------------
      Optical%totvo2(:,:,KS) = 0.0E+00
      do k=KS+1,KE+1
        Optical%totvo2(:,:,k) = Optical%totvo2(:,:,k-1) +   &
	                         Optical%cntval(:,:,k-1) 
      enddo
!----------------------------------------------------------------------



end  subroutine optical_rbts



!####################################################################

subroutine optical_h2o (pflux, atmden, vv, press, temp, rh2o, tflux, &
                        Optical) 

!----------------------------------------------------------------------
!     optical_h2o computes optical paths for h2o.
!
!     author: m. d. schwarzkopf
!
!     revised: 1/1/93
!
!     certified:  radiation version 1.0
!----------------------------------------------------------------------
!     intent in:
!
!     press   =  pressure at data levels of model.
!               
!     pflux   =  pressure at flux levels of model.
! 
!     rh2o    =  mass mixing ratio of h2o at model data levels 
!
!     temp    =  temperature at data levels of model. 
!
!     tpl1    =  temperature at "upper" nearby layers of model.
!
!     tpl2    =  temperature at "lower" nearby layers of model.
!-------------------------------------------------------------------    
real, dimension (:,:,:), intent(in) ::  pflux, atmden, vv, press, &
                      temp, rh2o, tflux
type(optical_path_type), intent(inout) :: Optical
!-----------------------------------------------------------------------
!     intent out:
!
!       empl1  =  h2o pressure scaled optical path between flux 
!                 level k and nearest data level k?
!
!       empl2  =  h2o pressure scaled optical path between flux 
!                 level k and nearest data level k+1?
!
!       emx1   =  h2o pressure scaled optical path between flux level
!                 KE and data level KE.
!
!       emx2   =  h2o pressure scaled optical path between flux level
!                 KE+1 and data level KE.
!
!       totm   =  summed h2o optical path from the top of atmosphere to 
!                 flux level k.
!
!       totphi =  summed pressure weighted h2o optical path from top
!                 of atmosphere to flux level k.
!
!       var1   =  h2o optical path in model layers.
!
!       var2   =  pressure-weighted h2o optical path in model layers.
!-----------------------------------------------------------------------
!     intent local:
!
!       qh2o    =  h2o mass mixing ratio, multiplied by the diffusivity
!                  factor diffac.
!-----------------------------------------------------------------------
      real, dimension (size(pflux,1), size(pflux,2), &
                     size(pflux,3)) :: tpl1, tpl2, &
                                             qh2o, tdif, tdif2

!!!  DAN NOTE ::
!!!   totfh2o, vrfh2o not used
!!  Are they still needed ??

!     real, dimension (size(temp,1), size(temp,2), size(temp,3), &
!		  NBTRGE)::  totfh2o

!     real, dimension (size(temp,1), size(temp,2), size(temp,3)-1,  &
!	       NBTRGE)::   vrfh2o

      integer    ::  m, k

!-------------------------------------------------------------------- 
!     compute mean temperature in the "nearby layer" between a flux
!     level and the first data level below the flux level (tpl1) or the
!     first data level above the flux level (tpl2)
!---------------------------------------------------------------------
      tpl1(:,:,KS   )         = temp(:,:,KE   )
      tpl1(:,:,KS   +1:KE   ) = tflux(:,:,KS   +1:KE   )
      tpl1(:,:,KE   +1)       = 0.5E+00*(tflux(:,:,KE   +1) +   &
                                       temp(:,:,KE   ))
   tpl2(:,:,KS   +1:KE   ) = tflux(:,:,KS   +1:KE   )
      tpl2(:,:,KE   +1)       = 0.5E+00*(tflux(:,:,KE   ) +    &
                                       temp(:,:,KE   ))

!---------------------------------------------------------------------
      allocate (Optical%empl1  (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1       ))
      allocate (Optical%empl2  (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1       ))
      allocate (Optical%totphi (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1       ))
      allocate (Optical%var1   (ISRAD:IERAD, JSRAD:JERAD  , KS:KE         ))
      allocate (Optical%var2   (ISRAD:IERAD, JSRAD:JERAD  , KS:KE         ))
      allocate (Optical%emx1   (ISRAD:IERAD, JSRAD:JERAD                  ))
      allocate (Optical%emx2   (ISRAD:IERAD, JSRAD:JERAD                  ))
!----------------------------------------------------------------------
!     compute optical paths for h2o, using the diffusivity 
!     approximation 1.66 for the angular integration.  obtain 
!     unweighted values var1, and weighted values var2.
!     the quantities 0.0003 (.0003) appearing in the
!     var2 expressions are the approximate voigt corrections
!     for h2o.  vv is the layer-mean pressure (in 
!     atmosphere), which is not the same as the level pressure press.
!---------------------------------------------------------------------  
      qh2o(:,:,KS:KE) = rh2o(:,:,KS:KE)*diffac
      Optical%var1(:,:,KS:KE) = atmden(:,:,KS:KE)*qh2o(:,:,KS:KE)
      Optical%var2(:,:,KS:KE) = Optical%var1(:,:,KS:KE)*(vv(:,:,KS:KE) + 3.0E-04)

!----------------------------------------------------------------------
!     compute summed optical paths for h2o.
!----------------------------------------------------------------------
      Optical%totphi(:,:,KS) = 0.0E+00
      do k=KS+1,KE+1
        Optical%totphi(:,:,k) = Optical%totphi(:,:,k-1) + Optical%var2  (:,:,k-1) 
      enddo

!----------------------------------------------------------------------
!     emx1 is the additional pressure-scaled mass from press(KE) to 
!     pflux(KE).  it is used in nearby layer and emiss calculations.
!     emx2 is the additional pressure-scaled mass from press(KE) to 
!     pflux(KE+1).  it is used in calculations between flux levels k
!     and KE+1.
!----------------------------------------------------------------------
      Optical%emx1(:,:) = qh2o(:,:,KE)*press(:,:,KE)*(press(:,:,KE) -   &
                  pflux(:,:,KE))/(grav*pstd)
      Optical%emx2(:,:) = qh2o(:,:,KE)*press(:,:,KE)*(pflux(:,:,KE+1) -  &
                  press(:,:,KE))/(grav*pstd)

!----------------------------------------------------------------------
!     empl is the pressure scaled mass from pflux(k) to press(k) or to 
!     press(k+1).
!----------------------------------------------------------------------
      Optical%empl1(:,:,KS) = Optical%var2(:,:,KE)
      Optical%empl1(:,:,KS+1:KE+1) =   &
                 qh2o(:,:,KS:KE)*pflux(:,:,KS+1:KE+1)*   &
                 (pflux(:,:,KS+1:KE+1) - press(:,:,KS:KE))/(grav*pstd)
      Optical%empl2(:,:,KS+1:KE) =    &
                 qh2o(:,:,KS+1:KE)*pflux(:,:,KS+1:KE)*   &
                 (press(:,:,KS+1:KE) - pflux(:,:,KS+1:KE))/(grav*pstd)
      Optical%empl2(:,:,KE+1) = Optical%empl2(:,:,KE) 

!---------------------------------------------------------------------
  if (Lw_control%do_ch4_n2o) then
    allocate ( Optical%empl1f (ISRAD:IERAD , JSRAD:JERAD , KS:KE+1,  NBTRGE ) )
    allocate ( Optical%empl2f (ISRAD:IERAD , JSRAD:JERAD , KS:KE+1,  NBTRGE ) )
    allocate ( Optical%tphfh2o(ISRAD:IERAD , JSRAD:JERAD , KS:KE+1,  NBTRGE ) ) 
    allocate ( Optical%vrpfh2o(ISRAD:IERAD , JSRAD:JERAD , KS:KE+1,  NBTRGE ) )
    allocate ( Optical%emx1f  (ISRAD:IERAD , JSRAD:JERAD ,           NBTRGE ) )
    allocate ( Optical%emx2f  (ISRAD:IERAD , JSRAD:JERAD ,           NBTRGE ) )

!----------------------------------------------------------------------
!   compute h2o optical paths for use in the 1200-1400 cm-1 range if
!   temperature dependence of line intensities is accounted for.
!----------------------------------------------------------------------
    tdif(:,:,KS:KE) = temp(:,:,KS:KE)-2.5E+02

    do m=1,NBTRGE
      Optical%vrpfh2o(:,:,KS:KE,m) = Optical%var2(:,:,KS:KE) *    &
                             EXP(csfah2o(1,m)*(tdif(:,:,KS:KE)) +   &
                                 csfah2o(2,m)*(tdif(:,:,KS:KE))**2 )
    enddo

    do m=1,NBTRGE
      Optical%tphfh2o(:,:,KS,m) = 0.0E+00
      do k=KS+1,KE+1
        Optical%tphfh2o(:,:,k,m) = Optical%tphfh2o(:,:,k-1,m) + Optical%vrpfh2o(:,:,k-1,m)
      enddo
    enddo

    tdif2(:,:,KS+1:KE+1) = tpl2(:,:,KS+1:KE+1)-2.5E+02
    tdif (:,:,KS+1:KE+1) = tpl1(:,:,KS+1:KE+1)-2.5E+02

!---------------------------------------------------------------------
!   compute this additional mass, for use in the 1200-1400 cm-1 range,
!   if temperature dependence of line intensities is accounted for.
!--------------------------------------------------------------------
    do m=1,NBTRGE
      Optical%emx1f(:,:,m) = Optical%emx1(:,:) *    &
                     EXP(csfah2o(1,m)*(tdif2(:,:,KE+1)        ) +   &
                         csfah2o(2,m)*(tdif2(:,:,KE+1)        )**2 )
      Optical%emx2f(:,:,m) = Optical%emx2(:,:) *    &
                     EXP(csfah2o(1,m)*(tdif (:,:,KE+1)        ) +  &
                         csfah2o(2,m)*(tdif (:,:,KE+1)        )**2 )
    enddo

!----------------------------------------------------------------------
!   compute this additional mass, for use in the 1200-1400 cm-1 range,
!   if temperature dependence of line intensities is accounted for.
!----------------------------------------------------------------------
    do m=1,NBTRGE
      Optical%empl1f(:,:,KS+1:KE+1,m) = Optical%empl1(:,:,KS+1:KE+1) *   &
                          EXP(csfah2o(1,m)*(tdif(:,:,KS+1:KE+1)) +    &
                              csfah2o(2,m)*(tdif(:,:,KS+1:KE+1))**2 )
      Optical%empl2f(:,:,KS+1:KE,m) = Optical%empl2(:,:,KS+1:KE) *   &
                          EXP(csfah2o(1,m)*(tdif2(:,:,KS+1:KE)) +   &
                              csfah2o(2,m)*(tdif2(:,:,KS+1:KE))**2 )
      Optical%empl1f(:,:,KS ,m) = Optical%vrpfh2o(:,:,KE,m)
      Optical%empl2f(:,:,KE+1,m) = Optical%empl2f(:,:,KE,m)
    enddo
  endif

!---------------------------------------------------------------------



end subroutine optical_h2o



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
      integer          :: kx
      real             :: rrf11, rrf12, rrf113, rrf22
     real :: rf11air, rf12air, rf113air, rf22air

!----------------------------------------------------------------------
!     allocate ( totf11 (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1     ) )
!     allocate ( totf12 (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1     ) )
!     allocate ( totf113(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1     ) )
!     allocate ( totf22 (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD+1     ) )
      allocate ( Optical%totf11 (size(density,1), size(density,2),    &
                         size(density,3) ) )
      allocate ( Optical%totf12 (size(density,1), size(density,2),    &
                         size(density,3) ) )
      allocate ( Optical%totf113(size(density,1), size(density,2),    &
                         size(density,3) ) )
      allocate ( Optical%totf22 (size(density,1), size(density,2),    &
                         size(density,3) ) )
 
      kx = size(density,3)

!--------------------------------------------------------------------
!   define cfc mixing ratio conversion factors.
!--------------------------------------------------------------------
      rf11air  = wtmf11/wtmair
       rf12air  = wtmf12/wtmair
    rf113air = wtmf113/wtmair
      rf22air  = wtmf22/wtmair

      rrf11 = Rad_gases%rrvf11*rf11air
      rrf12 = Rad_gases%rrvf12*rf12air
      rrf113 = Rad_gases%rrvf113*rf113air
      rrf22 = Rad_gases%rrvf22*rf22air

!----------------------------------------------------------------------
!   compute summed optical paths for f11,f12, f113 and f22  with the 
!   diffusivity factor of 2 (appropriate for weak-line absorption 
!   limit).
!----------------------------------------------------------------------
!     totf11(:,:,KSRAD) = 0.0E+00
!     totf12(:,:,KSRAD) = 0.0E+00
!     totf113(:,:,KSRAD) = 0.0E+00
!     totf22 (:,:,KSRAD) = 0.0E+00
      Optical%totf11(:,:,1) = 0.0E+00
      Optical%totf12(:,:,1) = 0.0E+00
      Optical%totf113(:,:,1) = 0.0E+00
      Optical%totf22 (:,:,1) = 0.0E+00
!     do k=KSRAD+1,KERAD+1
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


                   end module optical_path_mod

