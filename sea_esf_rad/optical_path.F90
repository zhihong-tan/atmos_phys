
                    module optical_path_mod
 

use rad_step_setup_mod,     only: temp, press, rh2o, KS=>KSRAD, &
				  KE=>KERAD, ISRAD, IERAD, JSRAD, &
				  JERAD, Rad_time_sv, lat_sv, pflux
use  longwave_params_mod,   only: NBLW, NBCO215, NBLY_ORIG
use  utilities_mod,         only: open_file, file_exist,    &
                                  check_nml_error, error_mesg, &
                                  print_version_number, FATAL, NOTE, &
				  WARNING, get_my_pe, close_file
use rad_utilities_mod,      only: looktab, longwave_tables3_type, &
				  table_alloc, longwave_control_type,&
				  Lw_control
use ch4_n2o_mod,            only: ch4_n2o_input_type, CN_basic
use cfc_mod,                only: cfc_optical_depth, cfc_exact,&
				  cfc_overod, cfc_overod_part,   &
				  cfc_exact_part, cfc_dealloc
use constants_new_mod,      only: wtmair, avogno, wtmh2o, grav,    &
				  rgas, pstd, diffac, rh2oair
use longwave_setup_mod,     only: pdfinv, pdflux, tpl1, tpl2,  &
			          dte1, ixoe1, Lw_parameters, &
			          longwave_parameter_type
use longwave_aerosol_mod,   only: aertau, longwave_aerosol_dealloc, &
                                  get_totaerooptdep,   &
			          get_totaerooptdep_15,   &
				  get_aerooptdep_KE_15
use gas_tf_mod,             only: get_gastf_for_optpath
use ozone_mod,              only: get_ozone

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!                   optical depth calculation module
!
!---------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: optical_path.F90,v 1.2 2001/08/30 15:13:20 fms Exp $'
   character(len=128)  :: tag     =  '$Name: galway $'


!---------------------------------------------------------------------
!----- interfaces  -----
           
public     &
	      optical_path_init, optical_dealloc,   &
	      optical_path_setup,     &
	      optical_trans_funct_from_KS,    &
	      optical_trans_funct_k_down, &
	      optical_trans_funct_KE,  &
	      optical_trans_funct_diag, &
	      get_totch2o, get_totch2obd, &
	      get_totvo2, get_var1var2, &
	      get_path_for_enear, get_path_for_e90


private    &
	      optical_ckd2p1_init, optical_path_ckd2p1, &
	      optical_o3, optical_rbts, optical_h2o
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

real, allocatable, dimension (:,:,:,:) :: empl1f, empl2f, vrpfh2o, &
                                          totch2o, totch2obd, &
                                          xch2obd, tphfh2o, avephif 
real, allocatable, dimension (:,:,:)   :: empl1, empl2,  var1, var2, &
                                          emx1f, emx2f, totvo2, avephi,&
                                          totch2obdwd, xch2obdwd, & 
                                          totphi, cntval,toto3,   &
					  tphio3, var3, var4, sh2o,  &
					  tmpexp, rvh2o, wk, rhoave, &
					  rh2os,  rfrgn, tfac
real, allocatable, dimension (:,:)     :: emx1, emx2, csfah2o

integer   :: NBTRG, NBTRGE



!---------------------------------------------------------------------
!---------------------------------------------------------------------



contains





subroutine optical_path_init

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
          csfah2o(1,m) = CN_basic%csf1h2o(1)
          csfah2o(2,m) = CN_basic%csf1h2o(2)
	end do
      elseif (NBTRGE .EQ. 2) then
        do m=1,NBTRGE
          csfah2o(1,m) = CN_basic%csf2h2o(1,m)
          csfah2o(2,m) = CN_basic%csf2h2o(2,m)
        end do
      elseif (NBTRGE .EQ. 4) then
        do  m=1,NBTRGE
          csfah2o(1,m) = CN_basic%csf4h2o(1,m)
          csfah2o(2,m) = CN_basic%csf4h2o(2,m)
        end do
      elseif (NBTRGE .EQ. 10) then
        do m=1,NBTRGE
          csfah2o(1,m) = CN_basic%csf10h2o(1,m)
          csfah2o(2,m) = CN_basic%csf10h2o(2,m)
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

subroutine optical_dealloc
      
!-------------------------------------------------------------------


      if (Lw_control%do_lwaerosol) then
        call longwave_aerosol_dealloc
      endif 

      if (Lw_control%do_cfc) then
        call cfc_dealloc
      endif

      deallocate  (var4 )
      deallocate  (var3 )
      deallocate  (tphio3)
      deallocate  (toto3)

      if (Lw_control%do_ckd2p1) then
        deallocate ( xch2obdwd)
        deallocate ( totch2obdwd)
        deallocate ( xch2obd)
      else
        deallocate (totvo2)
        deallocate (cntval)
      endif

      if (Lw_control%do_ch4_n2o) then
        deallocate ( emx2f  )
        deallocate ( emx1f  )
        deallocate ( vrpfh2o)
        deallocate ( tphfh2o)
        deallocate ( empl2f )
        deallocate ( empl1f )
      endif

      deallocate ( emx2   )
      deallocate ( emx1  )
      deallocate ( var2   )
      deallocate ( var1   )
      deallocate ( totphi )
      deallocate ( empl2  )
      deallocate ( empl1  )

      if (Lw_control%do_ch4_n2o) then
        deallocate (avephif)
      endif 

      deallocate (avephi)
      deallocate (tfac)
      deallocate (rfrgn)
      deallocate (rh2os)
      deallocate (rhoave)
      deallocate (wk  )
      deallocate (rvh2o)
      deallocate (tmpexp)
      deallocate (sh2o)

!--------------------------------------------------------------------


end subroutine optical_dealloc




!###################################################################

subroutine optical_path_setup 

!---------------------------------------------------------------------
     real, dimension(:,:,:), allocatable :: atmden, vv, qo3
     integer      :: k, i

!--------------------------------------------------------------------
!     atmden   =  atmospheric density, in gm/cm**2, for each of the
!                 KMAX layers.
!-------------------------------------------------------------------
     allocate (sh2o        (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (tmpexp      (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (rvh2o       (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (wk          (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (rhoave      (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (rh2os       (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (rfrgn       (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (tfac        (ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )
     allocate (avephi      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )

     if (Lw_control%do_ch4_n2o) then
       allocate (avephif   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1, NBTRGE) )
     endif
!---------------------------------------------------------------------
 
     allocate (atmden  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
     allocate (vv      (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )

     allocate (qo3 (ISRAD:IERAD, JSRAD:JERAD, KS:KE) )
     call get_ozone (qo3)

!----------------------------------------------------------------------
!     define the layer-mean pressure in atmospheres (vv) and the layer 
!     density (atmden). 
!----------------------------------------------------------------------
     do k=KS,KE
       atmden(:,:,k) = pdflux(:,:,k)/grav
       vv(:,:,k)     = 0.5E+00*(pflux(:,:,k+1) + pflux(:,:,k)  )/pstd
     enddo

!----------------------------------------------------------------------
!     compute optical paths.
!----------------------------------------------------------------------
     call optical_h2o (atmden, vv)

!---------------------------------------------------------------------
     if (.not. Lw_control%do_ckd2p1) then
       call optical_rbts 
     else
!---------------------------------------------------------------------
!    call optical_ckd2p1 to determine self- and foreign-broadened h2o
!    continuum paths, for the given temperature, pressure and mixing
!    ratio, over the predetermined frequency range.
!---------------------------------------------------------------------
       call optical_path_ckd2p1  (atmden)
     endif

!---------------------------------------------------------------------
     call optical_o3 (atmden, qo3, vv)

!--------------------------------------------------------------------
     if (Lw_control%do_cfc) then
       call cfc_optical_depth (atmden) 
     endif

!---------------------------------------------------------------------
!     compute aerosol layer transmission functions for all layers.
!     option predlwaer is planned, but not yet available. when  it 
!     becomes available,  aeralb, aerext and aerconc will be additional 
!     arguments going to aertau.
!---------------------------------------------------------------------
     if (Lw_control%do_lwaerosol) then
       call aertau 
     endif
!---------------------------------------------------------------------
       
     deallocate (qo3)
     deallocate (vv      )
     deallocate (atmden  )
 

end subroutine  optical_path_setup



!####################################################################

subroutine optical_trans_funct_from_KS (to3cnt, overod, &
				        cnttaub1, cnttaub2, cnttaub3)

!---------------------------------------------------------------------
real, dimension (:,:,:), intent(out) ::  to3cnt, overod, &
                                         cnttaub1, cnttaub2, cnttaub3

!---------------------------------------------------------------------
      real, dimension (:,:,:), allocatable :: tmp1, tmp2, tmp3,   &
					      cfc_tf, totch2o_tmp,  &
					      totaer_tmp, tn2o17,  &
					      totaerooptdep_15
      integer    :: m

!--------------------------------------------------------------------
!  allocate localk arrays
!--------------------------------------------------------------------
      allocate ( tmp1(ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1  ) )
      allocate ( tmp2(ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1  ) )
      allocate ( tmp3(ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1  ) )
      allocate ( cfc_tf(ISRAD:IERAD, JSRAD:JERAD,       KS:KE  ) )
      allocate ( totch2o_tmp(ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1) )
      allocate ( totaer_tmp(ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1) )
      allocate ( tn2o17(ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1  ) )

!-----------------------------------------------------------------------
!   compute transmission functions in 990-1070 cm-1 range, including
!   ozone and h2o continuum, from level KS to all other levels. 
!------------------------------------------------------------------
      tmp1  (:,:,KS:KE) = bo3rnd(2)*tphio3(:,:,KS+1:KE+1)/  &
                           toto3(:,:,KS+1:KE+1)
      tmp2(:,:,KS:KE) = 0.5*(tmp1(:,:,KS:KE)*(SQRT(1.0E+00 + (4.0E+00*&
                        ao3rnd(2)*toto3(:,:,KS+1:KE+1))/  &
			tmp1(:,:,KS:KE))  - 1.0E+00))
      if (Lw_control%do_ckd2p1) then
	call get_totch2obd(6, totch2o_tmp)
	tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) + diffac*   &
		          totch2o_tmp(:,:,KS+1:KE+1)
      else 
        tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) + betacm(14)*   &
                          totvo2(:,:,KS+1:KE+1)
      endif

      if (Lw_control%do_lwaerosol) then
	call get_totaerooptdep(6, totaer_tmp)
	tmp2(:,:,KS:KE) = tmp2(:,:,KS:KE) +    &
                          totaer_tmp   (:,:,KS+1:KE+1)
      endif
      to3cnt(:,:,KS:KE) = EXP(-1.0E+00*tmp2(:,:,KS:KE))
 
!--------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in to3cnt.
!---------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_exact (6, cfc_tf)
        to3cnt(:,:,KS:KE) = to3cnt(:,:,KS:KE)* cfc_tf(:,:,KS:KE)
      endif

!---------------------------------------------------------------------
!   compute transmission function in the 560-800 cm-1 range
!   evaluate  optical depth contributions 
!   add contributions from h2o(lines) and h2o(continuum).
!   h2o(continuum) contributions are either Roberts or CKD2.1
!---------------------------------------------------------------------
      tmp1(:,:,KS:KE) = SQRT(ab15wd*totphi(:,:,KS+1:KE+1)) 
      if (Lw_control%do_ckd2p1) then
	tmp1(:,:,KS:KE) = tmp1(:,:,KS:KE) + diffac*   &
                          totch2obdwd(:,:,KS+1:KE+1)
      else
	tmp1(:,:,KS:KE) = tmp1(:,:,KS:KE) + betawd*   &
                          totvo2     (:,:,KS+1:KE+1)
      endif

!--------------------------------------------------------------------
!   add contribution from longwave aerosols (if desired).
!--------------------------------------------------------------------
      if (Lw_control%do_lwaerosol) then
	allocate ( totaerooptdep_15(ISRAD:IERAD, JSRAD:JERAD,   &
						        KS:KE+1))
	call get_totaerooptdep_15 (totaerooptdep_15)
	tmp1(:,:,KS:KE) = tmp1(:,:,KS:KE) +    &
                          totaerooptdep_15(:,:,KS+1:KE+1)
	deallocate (totaerooptdep_15)
      endif
 
!----------------------------------------------------------------------
!    compute transmission function due to these contributions. the
!    effects of co2, n2o  and  cfc's (not exponentials) are added
!    later.
!---------------------------------------------------------------------
      overod(:,:,KS:KE) = EXP(-1.0E+00*tmp1     (:,:,KS:KE))

!---------------------------------------------------------------------
!   add contribution from the 17 um n2o band (if desired).
!   the expression with tn2o17 retains the 560-630 cm-1 equi-
!   valent widths in evaluating 560-800 cm-1 transmissivities.
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        call get_gastf_for_optpath (KS+1, KE+1, tn2o17)

        if (NBCO215 .EQ. 2) then
          overod(:,:,KS:KE) = overod(:,:,KS:KE) *    &
                              (130./240. + 110./240.*   &
			       tn2o17(:,:,KS+1:KE+1))
        elseif (NBCO215 .EQ. 3) then
          overod(:,:,KS:KE) = overod(:,:,KS:KE) * (170./240. +    &
			      70./240.*tn2o17(:,:,KS+1:KE+1))
        endif
      endif

!--------------------------------------------------------------------- 
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod .
!--------------------------------------------------------------------
      if (Lw_control%do_cfc) then
	call cfc_overod (cfc_tf)
	overod(:,:,KS:KE) = overod(:,:,KS:KE)*cfc_tf(:,:,KS:KE)
      endif 

!----------------------------------------------------------------------
!     compute continuum band transmission functions from level KS to
!     other levels (cnttau). the continuum transmission function from
!     level k to kp (contod) equals cnttau for k=KS, so is not
!     evaluated here. for all other levels k, contod is obtained by
!     division of relevant values of cnttau.
!---------------------------------------------------------------------
      if (Lw_control%do_ckd2p1) then 
	call get_totch2obd(4, totch2o_tmp)
	tmp1(:,:,KS:KE) = diffac*totch2o_tmp(:,:,KS+1:KE+1)
	call get_totch2obd(5, totch2o_tmp)
	tmp2(:,:,KS:KE) = diffac*totch2o_tmp(:,:,KS+1:KE+1)
	call get_totch2obd(7, totch2o_tmp)
	tmp3(:,:,KS:KE) = diffac*totch2o_tmp(:,:,KS+1:KE+1)
      else
	tmp1(:,:,KS:KE) = betacm(12)*totvo2(:,:,KS+1:KE+1)
	tmp2(:,:,KS:KE) = betacm(13)*totvo2(:,:,KS+1:KE+1)
	tmp3(:,:,KS:KE) = betacm(15)*totvo2(:,:,KS+1:KE+1)
      endif

      if (Lw_control%do_lwaerosol) then
	call get_totaerooptdep(4, totaer_tmp)
	tmp1(:,:,KS:KE) =  tmp1(:,:,KS:KE) +    &
                           totaer_tmp   (:,:,KS+1:KE+1  )
	call get_totaerooptdep(5, totaer_tmp)
	tmp2(:,:,KS:KE) =  tmp2(:,:,KS:KE) +    &
                           totaer_tmp   (:,:,KS+1:KE+1)
	call get_totaerooptdep(7, totaer_tmp)
	tmp3(:,:,KS:KE) =  tmp3(:,:,KS:KE) +    &
                           totaer_tmp   (:,:,KS+1:KE+1)
      endif

      cnttaub1(:,:,KS:KE) = EXP(-1.0*tmp1(:,:,KS:KE))
      cnttaub2(:,:,KS:KE) = EXP(-1.0*tmp2(:,:,KS:KE))
      cnttaub3(:,:,KS:KE) = EXP(-1.0*tmp3(:,:,KS:KE))

!---------------------------------------------------------------------
!     if cfcs are included, add transmission functions for f11, f12,
!     f113, and f22.
!---------------------------------------------------------------------
      if (Lw_control%do_cfc) then
        call cfc_exact (4, cfc_tf)
        cnttaub1(:,:,KS:KE) = cnttaub1(:,:,KS:KE)*cfc_tf(:,:,KS:KE)
        call cfc_exact (5, cfc_tf)
        cnttaub2(:,:,KS:KE) = cnttaub2(:,:,KS:KE)*cfc_tf(:,:,KS:KE)
        call cfc_exact (7, cfc_tf)
        cnttaub3(:,:,KS:KE) = cnttaub3(:,:,KS:KE)*cfc_tf(:,:,KS:KE)
      endif 
 
!----------------------------------------------------------------------
!     evaluate h2o (mbar*phibar) between level KS and other levels.
!----------------------------------------------------------------------
      avephi(:,:,KS:KE) = totphi(:,:,KS+1:KE+1)
 
!----------------------------------------------------------------------
!     the evaluation of emiss over the layer between data level (KS)
!     and flux level (KE+1) is done by averaging E2 functions referring
!     to the top and bottom of the layer. a special value of (mbar*
!     phibar) is required; it is stored in the (otherwise vacant)
!     KE+1'th position of avephi.
!----------------------------------------------------------------------
      avephi(:,:,KE+1) = avephi(:,:,KE-1) + emx1(:,:)

!----------------------------------------------------------------------
!     if h2o lines in the 1200-1400 range are assumed to have a temp-
!     erature dependent intensity, similar evaluation for (mbar*phibar)
!     is performed, with a special value for the lowest layer
!----------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        do m=1,NBTRGE
         avephif(:,:,KS:KE,m) = tphfh2o(:,:,KS+1:KE+1,m)
        enddo
        do m=1,NBTRGE
         avephif(:,:,KE+1,m) = avephif(:,:,KE-1,m) + emx1f(:,:,m)
        enddo
      endif
!----------------------------------------------------------------------

      deallocate (tn2o17)
      deallocate (totaer_tmp)
      deallocate (totch2o_tmp)
      deallocate (cfc_tf)
      deallocate (tmp3 )
      deallocate (tmp2  )
      deallocate (tmp1  )


end subroutine optical_trans_funct_from_KS




!####################################################################

subroutine optical_trans_funct_k_down (k, cnttaub1, cnttaub2, &
				       cnttaub3, to3cnt, overod,  &
			               contodb1, contodb2, contodb3)

!----------------------------------------------------------------------
integer,                 intent (in) :: k
real, dimension(:,:,:),  intent(in)  :: cnttaub1, cnttaub2, cnttaub3
real, dimension (:,:,:), intent(out) :: to3cnt, overod, contodb1,   &
					contodb2, contodb3
!---------------------------------------------------------------------

       real, dimension(:,:,:), allocatable  ::  avmo3, avpho3, tmp1, &
					        tmp2, cfc_tf, avvo2,  &
					        avckdwd, avckdo3, &
	                                        avaero3, totch2o_tmp,  &
					        totaer_tmp, tn2o17, &
	                      		        totaerooptdep_15

       integer       :: kp, m
!--------------------------------------------------------------------

       allocate (avmo3   (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
       allocate (avpho3  (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
       allocate (tmp1    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
       allocate (tmp2    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
       allocate (cfc_tf  (ISRAD:IERAD, JSRAD:JERAD, KS:KE  )  )

!---------------------------------------------------------------------
       if (Lw_control%do_ckd2p1) then
	 allocate (avckdwd (ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1  ))
	 allocate (avckdo3 (ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1  ))
	 allocate (totch2o_tmp (ISRAD:IERAD, JSRAD:JERAD,      KS:KE+1))
       else
	 allocate (avvo2   (ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1  ))
       endif

       if (Lw_control%do_lwaerosol) then
	 allocate (totaer_tmp (ISRAD:IERAD, JSRAD:JERAD,      KS:KE+1))
	 allocate (avaero3 (ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1  ))
       endif

       allocate (tn2o17  (ISRAD:IERAD, JSRAD:JERAD,       KS:KE+1  ))

       if (Lw_control%do_ckd2p1) then
         call get_totch2obd(6, totch2o_tmp)
       endif
       if (Lw_control%do_lwaerosol) then
         call get_totaerooptdep(6, totaer_tmp)
       endif

       do kp=1,KE+1-k
         avmo3 (:,:,kp+k-1) = toto3 (:,:,kp+k) - toto3 (:,:,k)
         avpho3(:,:,kp+k-1) = tphio3(:,:,kp+k) - tphio3(:,:,k) 
         if (.not. Lw_control%do_ckd2p1) then
           avvo2 (:,:,kp+k-1) = totvo2(:,:,kp+k) - totvo2(:,:,k)
         else
           avckdwd(:,:,kp+k-1) = totch2obdwd(:,:,kp+k) -   &
                                 totch2obdwd(:,:,k)
           avckdo3(:,:,kp+k-1) = totch2o_tmp(:,:,kp+k) -  &
                                 totch2o_tmp(:,:,k)
	 endif			 
         if (Lw_control%do_lwaerosol) then
           avaero3(:,:,kp+k-1) =  &
                         totaer_tmp   (:,:,kp+k) - totaer_tmp   (:,:,k)
	 endif
       enddo

       do kp=1,KE+1-k
         contodb1(:,:,kp+k-1) = cnttaub1(:,:,kp+k-1)/cnttaub1(:,:,k-1)
         contodb2(:,:,kp+k-1) = cnttaub2(:,:,kp+k-1)/cnttaub2(:,:,k-1)
         contodb3(:,:,kp+k-1) = cnttaub3(:,:,kp+k-1)/cnttaub3(:,:,k-1)
         avephi  (:,:,kp+k-1) = totphi(:,:,kp+k) - totphi(:,:,k)
       end do
       avephi (:,:,KE+1) = avephi(:,:,KE-1) + emx1(:, :)


!---------------------------------------------------------------------
!     if h2o lines in the 1200-1400 range are assumed to have a temp-
!     erature dependent intensity, similar evaluation for (mbar*phibar)
!     is performed, with a special value for the lowest layer
!---------------------------------------------------------------------
       if (Lw_control%do_ch4_n2o) then
         do m=1,NBTRGE
           do kp=1,KE+1-k
             avephif(:,:,kp+k-1,m) =  tphfh2o(:,:,kp+k,m) -  &
                                  tphfh2o(:,:,k,   m)
           enddo
           avephif(:,:,KE+1,m) = avephif(:,:,KE-1,m) + emx1f(:,:,m)
         enddo
       endif

!----------------------------------------------------------------------
!   compute transmission function in the 560-800 cm-1 range
!   evaluate  optical depth contributions 
!
!   add contributions from h2o(lines) and h2o(continuum).
!   h2o(continuum) contributions are either Roberts or CKD2.1
!----------------------------------------------------------------------
       tmp1     (:,:,k:KE) = SQRT(ab15wd*avephi(:,:,k:KE)) 
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
	 allocate ( totaerooptdep_15(ISRAD:IERAD, JSRAD:JERAD,   &
						        KS:KE+1))
	 call get_totaerooptdep_15 (totaerooptdep_15)
	 do kp=k,KE
	   tmp1(:,:,kp) = tmp1(:,:,kp) +    &
                          (totaerooptdep_15(:,:,kp+1) -   &
                           totaerooptdep_15(:,:,k) )
	 end do
	 deallocate (totaerooptdep_15)
       endif

!---------------------------------------------------------------------
!    compute transmission function due to these contributions. the
!    effects of co2, n2o  and  cfc's (not exponentials) are added
!    later.
!--------------------------------------------------------------------
       overod(:,:,k:KE) = EXP(-1.0E+00*tmp1     (:,:,k:KE))

!----------------------------------------------------------------------
!       add contribution from the 17 um n2o band (if desired).
!   the expression with tn2o17 retains the 560-630 cm-1 equi-
!   valent widths in evaluating 560-800 cm-1 transmissivities.
!---------------------------------------------------------------------
       if (Lw_control%do_ch4_n2o) then
         call get_gastf_for_optpath (k+1, KE+1, tn2o17)

         if (NBCO215 .EQ. 2) then
           overod(:,:,k:KE) = overod(:,:,k:KE) *(130./240. +  &
		 	      110./240.*tn2o17(:,:,k+1:KE+1))
         elseif (NBCO215 .EQ. 3) then
           overod(:,:,k:KE) = overod(:,:,k:KE)*(170./240. + &
			      70./240.*tn2o17(:,:,k+1:KE+1))
         endif
       endif

!----------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod .
!----------------------------------------------------------------------
       if (Lw_control%do_cfc) then
         call cfc_overod_part ( cfc_tf, k)
         overod(:,:,k:KE) = overod(:,:,k:KE)*cfc_tf(:,:,k:KE)
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
       to3cnt(:,:,k:KE) = EXP(-1.0E+00*tmp2(:,:,k:KE))

!---------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in to3cnt.
!---------------------------------------------------------------------
       if (Lw_control%do_cfc) then
         call cfc_exact_part (6, cfc_tf, k)
         to3cnt(:,:,k:KE) = to3cnt(:,:,k:KE)*cfc_tf(:,:,k:KE)
       endif 
!---------------------------------------------------------------------

       deallocate (tn2o17)
       if (allocated (avvo2)) then 
	 deallocate (avvo2)
       endif
       if (allocated (avckdwd)) then
	 deallocate (avckdwd)
	 deallocate (avckdo3)
	 deallocate (totch2o_tmp)
       endif
       if (allocated (avaero3)) then
	 deallocate (avaero3)
	 deallocate (totaer_tmp )
       endif

       deallocate (cfc_tf   )
       deallocate (tmp2     )
       deallocate (tmp1     )
       deallocate (avpho3   )
       deallocate (avmo3    )


end subroutine optical_trans_funct_k_down



!#################################################################

subroutine optical_trans_funct_KE (cnttaub1, cnttaub2,    &
				   cnttaub3, to3cnt, overod, contodb1, &
				   contodb2, contodb3)

!---------------------------------------------------------------------
real, dimension (:,:,:), intent(in) ::  cnttaub1, cnttaub2,  cnttaub3
real, dimension (:,:,:), intent(out) :: to3cnt, overod, contodb1, &
					  contodb2, contodb3

!---------------------------------------------------------------------
      real, dimension(:,:,:), allocatable :: tmp1, tmp2, tn2o17, cfc_tf
      real,  dimension(:,:) , allocatable :: aerooptdep_KE_15

!---------------------------------------------------------------------

      allocate (tmp1    ( ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate (tmp2    ( ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate (tn2o17  ( ISRAD:IERAD, JSRAD:JERAD, KS:KE+1) )
      allocate (cfc_tf  ( ISRAD:IERAD, JSRAD:JERAD, KS:KE  ) )

!-----------------------------------------------------------------------
!   compute transmission function in the 560-800 cm-1 range
!
!   evaluate  optical depth contributions 
!
!       add contributions from h2o(lines) and h2o(continuum).
!       h2o(continuum) contributions are either Roberts or CKD2.1
!----------------------------------------------------------------------
      tmp1     (:,:,KE) = SQRT(ab15wd*var2  (:,:,KE)) 
      if (Lw_control%do_ckd2p1) then
	tmp1(:,:,KE) = tmp1(:,:,KE) + diffac*   &
                         xch2obdwd   (:,:,KE)
      else
	tmp1(:,:,KE) = tmp1(:,:,KE) + betawd*  &
                          cntval     (:,:,KE)
      endif

!---------------------------------------------------------------------
!      add contribution from longwave aerosols (if desired).
!---------------------------------------------------------------------
      if (Lw_control%do_lwaerosol) then
	allocate (aerooptdep_KE_15(ISRAD:IERAD, JSRAD:JERAD) )
	call get_aerooptdep_KE_15 (aerooptdep_KE_15)
	tmp1(:,:,KE) = tmp1(:,:,KE) +     &
                          aerooptdep_KE_15(:,:)  
	deallocate (aerooptdep_KE_15)
      endif
 
!---------------------------------------------------------------------
!    compute transmission function due to these contributions. the
!    effects of co2, n2o  and  cfc's (not exponentials) are added
!    later.
!---------------------------------------------------------------------
      overod(:,:,KE) = EXP(-1.0E+00*tmp1     (:,:,KE))
 
!---------------------------------------------------------------------
!       add contribution from the 17 um n2o band (if desired).
!   the expression with tn2o17 retains the 560-630 cm-1 equi-
!   valent widths in evaluating 560-800 cm-1 transmissivities.
!---------------------------------------------------------------------
      if (Lw_control%do_ch4_n2o) then
        call get_gastf_for_optpath (KE+1, KE+1, tn2o17)
        if (NBCO215 .EQ. 2) then
          overod(:,:,KE) = overod(:,:,KE) *  &
                           (130./240. + 110./240.*tn2o17(:,:,KE+1))
        elseif (NBCO215 .EQ. 3) then
          overod(:,:,KE) = overod(:,:,KE) *   &
                           (170./240. + 70./240.*tn2o17(:,:,KE+1))
        endif
      endif

!---------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod .
!---------------------------------------------------------------------
     if (Lw_control%do_cfc) then
       call cfc_overod_part (cfc_tf, KE)
       overod(:,:,KE) = overod(:,:,KE)*cfc_tf(:,:,KE)
     endif 

!-----------------------------------------------------------------------
!   compute transmission functions in 990-1070 cm-1 range, including
!   ozone and h2o continuum, from level KS to all other levels. 
!---------------------------------------------------------------------
     tmp1  (:,:,KE) = bo3rnd(2)*var4(:,:,KE)/var3(:,:,KE)
     tmp2(:,:,KE) = 0.5*(tmp1(:,:,KE)*(SQRT(1.0E+00 + (4.0E+00*  &
                    ao3rnd(2)*var3 (:,:,KE))/tmp1(:,:,KE))  - 1.0E+00))

     if (Lw_control%do_ckd2p1) then
       tmp2(:,:,KE) = tmp2(:,:,KE) + diffac*xch2obd  (:,:,KE,6) 
     else
       tmp2(:,:,KE) = tmp2(:,:,KE) + betacm(14)*cntval (:,:,KE)
     endif

     to3cnt(:,:,KE) = EXP(-1.0E+00*tmp2(:,:,KE))

!---------------------------------------------------------------------
!    if cfcs are included, also include the transmission functions for
!    f11, f12, f113, and f22 in overod and to3cnt.
!---------------------------------------------------------------------
     if (Lw_control%do_cfc) then
       call cfc_exact_part (6, cfc_tf, KE)
       to3cnt(:,:,KE) = to3cnt(:,:,KE)*cfc_tf(:,:,KE)
     endif

     contodb1(:,:,KE  ) = cnttaub1(:,:,KE)/cnttaub1(:,:,KE-1)
     contodb2(:,:,KE  ) = cnttaub2(:,:,KE)/cnttaub2(:,:,KE-1)
     contodb3(:,:,KE  ) = cnttaub3(:,:,KE)/cnttaub3(:,:,KE-1)

!-------------------------------------------------------------------
     deallocate (cfc_tf  )
     deallocate (tn2o17  )
     deallocate (tmp2    )
     deallocate (tmp1    )


end subroutine optical_trans_funct_KE




!####################################################################

subroutine optical_trans_funct_diag ( press, contdg, to3dg)

!---------------------------------------------------------------------
real, dimension (:,:,:), intent(in)     :: press
real, dimension (:,:,:), intent(out)    :: to3dg                
real, dimension (:,:,:,:), intent(out)  :: contdg               
					
!---------------------------------------------------------------------
!   local variables:
!---------------------------------------------------------------------
    real, dimension(:,:,:), allocatable ::    &
                   ca, cb, csuba, csubb, ctmp2, ctmp3, delpr1, delpr2

!----------------------------------------------------------------------

    allocate (  ca        (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
    allocate (  cb        (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
    allocate (  csuba     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
    allocate (  csubb     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
    allocate (  ctmp2     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
    allocate (  ctmp3     (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
    allocate (  delpr1    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )
    allocate (  delpr2    (ISRAD:IERAD, JSRAD:JERAD, KS:KE+1)  )

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
      ctmp2(:,:,KS+1:KE)  = cntval(:,:,KS+1:KE)*delpr1(:,:,KS+1:KE) 
      ctmp3(:,:,KS+1:KE)  = cntval(:,:,KS:KE-1)*delpr2(:,:,KS+1:KE) 
    endif
    
!-----------------------------------------------------------------------
!    compute sf2.
!    continuum band 1
!-----------------------------------------------------------------------
    if ( .not. Lw_control%do_ckd2p1) then
      csuba(:,:,KS+1:KE)  = betacm(12)*ctmp2(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = betacm(12)*ctmp3(:,:,KS+1:KE)
	else
      csuba(:,:,KS+1:KE)  = diffac*xch2obd(:,:,KS+1:KE,4)*  &
                            delpr1(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = diffac*xch2obd(:,:,KS:KE-1,4)*  &
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
      csuba(:,:,KS+1:KE)  = diffac*xch2obd(:,:,KS+1:KE,5)*   &
                            delpr1(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = diffac*xch2obd(:,:,KS:KE-1,5)*  &
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
      csuba(:,:,KS+1:KE)  = diffac*xch2obd(:,:,KS+1:KE,7)*   &
                            delpr1(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = diffac*xch2obd(:,:,KS:KE-1,7)*  &
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
      csuba(:,:,KS+1:KE)  = diffac*xch2obd(:,:,KS+1:KE,6)*   &
                            delpr1(:,:,KS+1:KE)
      csubb(:,:,KS+1:KE)  = diffac*xch2obd(:,:,KS:KE-1,6)*  &
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

     deallocate ( ca       )
     deallocate ( cb       )
     deallocate ( csuba    )
     deallocate ( csubb    )
     deallocate ( ctmp2    )
     deallocate ( ctmp3    )
     deallocate ( delpr1   )
     deallocate ( delpr2   )


end subroutine optical_trans_funct_diag


!###################################################################

subroutine get_totch2o (n, totch2o)

!------------------------------------------------------------------
real, dimension(:,:,:), intent(out)     :: totch2o
integer,                intent(in)      :: n

!-----------------------------------------------------------------
         real, allocatable  ::  radf (:,:,:)
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
         tmpexp(:,:,KS:KE) = EXP(-.013*tfac(:,:,KS:KE))

!--------------------------------------------------------------------
!     compute source function for frequency bands (ioffh2o+1 to ioffh2o
!     +nptch2o) at layer temperatures using table lookup.
!     note that ixoe1 can be used for temp index, and dte1 for deltat,
!     as the table extent for radf is the same as for the e1 tables
!     of the model.
!--------------------------------------------------------------------
         allocate (radf (ISRAD:IERAD, JSRAD:JERAD,         KS:KE) )
         nu = n
         call looktab (radfunc, ixoe1, dte1, radf(:,:,:), KS, KE,   &
		       nu+ioffh2o)
         sh2o0 = ssh2o_296(nu+ioffh2o)*sfac(nu+ioffh2o)

         do k=KS,KE 
 	   sh2o(:,:,k) = sh2o0*tmpexp(:,:,k)
         enddo
 
!--------------------------------------------------------------------
!     compute h2o self- and foreign- broadened continuum optical path,
!     summed from the top of the atmosphere through layer k.
!--------------------------------------------------------------------
         fh2o0 = sfh2o(nu+ioffh2o)*fscal(nu+ioffh2o)
         totch2o(:,:,1) = 0.0E+00
	 do k = KS+1,KE+1
	   totch2o(:,:,k) = wk(:,:,k-1)*1.0e-20*   &
                            (sh2o(:,:,k-1)*rh2os(:,:,k-1) +    &
	                    fh2o0*rfrgn(:,:,k-1))* &
                            vvj(nu)*radf(:,:,k-1   )    +   &
                            totch2o(:,:,k-1)
         enddo

!------------------------------------------------------------------
         deallocate (radf)

!------------------------------------------------------------------

end subroutine get_totch2o



!#####################################################################


subroutine get_totch2obd (n, totch2obd)

!------------------------------------------------------------------
real, dimension(:,:,:), intent(out)     :: totch2obd
integer,                intent(in)      :: n

!-----------------------------------------------------------------
        real, allocatable  ::  radf (:,:,:)
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
                             xch2obd(:,:,k-1,nu)
	enddo
!--------------------------------------------------------------------
 
end subroutine get_totch2obd




!#####################################################################

subroutine get_totvo2 (n, totvo2_out) 

!------------------------------------------------------------------
integer,                intent(in)         :: n
real, dimension(:,:,:), intent(out)        :: totvo2_out

!-----------------------------------------------------------------

        totvo2_out(:,:,:) = betacm(n)*totvo2(:,:,KS+1:KE+1)

end subroutine get_totvo2 


!###################################################################

subroutine get_var1var2 (var1_out, var2_out)

!------------------------------------------------------------------
real, dimension(:,:,:), intent(out)        :: var1_out, var2_out

!-----------------------------------------------------------------

	var1_out(:,:,:) = var1(:,:,:)
	var2_out(:,:,:) = var2(:,:,:)

end subroutine get_var1var2              


!###################################################################

subroutine get_path_for_enear ( var2_out, emx1_out, emx2_out,  &
			       empl1_out,  empl2_out, &
			       emx1f_out, emx2f_out, empl1f_out, &
			       empl2f_out, vrpfh2o_out)

!---------------------------------------------------------------------
real, dimension(:,:,:),   intent(out)            :: empl1_out,  &
						    empl2_out, &
                                                    var2_out 
real, dimension(:,:),     intent(out)            :: emx1_out, emx2_out
real, dimension(:,:,:,:), intent(out), optional  :: empl1f_out,   &
					            empl2f_out, &
					            vrpfh2o_out
real, dimension(:,:,:),   intent(out), optional  :: emx1f_out, &
                                                    emx2f_out
!---------------------------------------------------------------------

      if (Lw_control%do_ch4_n2o) then
         empl1f_out    = empl1f
         empl2f_out    = empl2f
	 vrpfh2o_out   = vrpfh2o
	 emx1f_out     = emx1f
	 emx2f_out     = emx2f
      endif

      empl1_out     = empl1
      empl2_out     = empl2
      var2_out      = var2
      emx1_out      = emx1
      emx2_out      = emx2

!--------------------------------------------------------------------


end subroutine get_path_for_enear



!#####################################################################

subroutine get_path_for_e90 (avephi_out, avephif_out)

real, dimension(:,:,:),   intent(out)             :: avephi_out
real, dimension(:,:,:,:), intent(out), optional   :: avephif_out


      avephi_out(:,:,:)    = avephi(:,:,:)

      if (Lw_control%do_ch4_n2o) then
        avephif_out(:,:,:,:) = avephif(:,:,:,:)
      endif

end subroutine get_path_for_e90



!#####################################################################

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

subroutine optical_path_ckd2p1 (atmden) 

!---------------------------------------------------------------------
!    subroutine optical_ckd2p1 computes h2o continuum optical paths
!    (self + foreign) over the frequency range specified by
!    ioffh2o and nptch2o using the ckd2.1 algorithm, modified for 
!    the gcm parameterization.
!    (this routine is previously called contnm.F)
!---------------------------------------------------------------------
real, dimension (:,:,:), intent(in)       ::  atmden

!---------------------------------------------------------------------
      real, dimension(:,:,:), allocatable :: xch2obdinw, totch2obdinw, &
					     radf
      real                                ::  t0 = 296.0
      integer                             ::  n, k, nu
      real                                ::  fh2o0, sh2o0

!----------------------------------------------------------------------

      allocate (xch2obd    (ISRAD:IERAD, JSRAD:JERAD,   KS:KE  , 7) )
      allocate (totch2obdwd(ISRAD:IERAD, JSRAD:JERAD,   KS:KE+1   ) )
      allocate (xch2obdwd  (ISRAD:IERAD, JSRAD:JERAD,   KS:KE     ) )
      allocate (xch2obdinw  (ISRAD:IERAD, JSRAD:JERAD,   KS:KE  ) )
      allocate (totch2obdinw(ISRAD:IERAD, JSRAD:JERAD,   KS:KE+1) )

!--------------------------------------------------------------------
!    define the volume mixing ratio of h2o
!---------------------------------------------------------------------
      rvh2o(:,:,KS:KE) = rh2o(:,:,KS:KE)*wtmair/wtmh2o

!---------------------------------------------------------------------
!    define input arguments to optical_ckd2p1
!-------------------------------------------------------------------
      wk(:,:,KS:KE) = rvh2o(:,:,KS:KE)*avogno/wtmair*   &
                          atmden(:,:,KS:KE)
      rhoave(:,:,KS:KE) = (press(:,:,KS:KE)/pstd)*(t0/temp(:,:,KS:KE))
      rh2os(:,:,KS:KE) = rvh2o(:,:,KS:KE)*rhoave(:,:,KS:KE)
      rfrgn(:,:,KS:KE) = rhoave(:,:,KS:KE) - rh2os(:,:,KS:KE)
      tfac(:,:,KS:KE) = temp(:,:,KS:KE) - t0

!--------------------------------------------------------------------
!     compute self-broadened temperature-dependent continuum coefficient
!     using the single coefficient -.020 for all frequencies in
!     the 560-1200 cm-1 range. experiments with the mid-latitude
!     summer profile show errors of < .01 W/m**2 (in the net broadband
!     flux, 0-2200 cm-1) using this value. this value is used instead
!     of tmpfctrs at each frequency band.
!-------------------------------------------------------------------
      tmpexp(:,:,KS:KE) = EXP(-.020*tfac(:,:,KS:KE))
 
!-------------------------------------------------------------------
!     compute h2o self- and foreign- broadened continuum optical path 
!     for each layer k (xch2obd, xch2obdinw, xch2obdwd) and summed from
!     the top of the atmosphere through layer k (totch2obd,
!     totch2obdinw, totch2obdwd).
!--------------------------------------------------------------------
      do nu = 1,7
        do k = KS,KE 
          xch2obd(:,:,k,nu) = wk(:,:,k)*1.0e-20* (svj(nu)*rh2os(:,:,k)*&
			      tmpexp(:,:,k) + fvj(nu)*rfrgn(:,:,k))* &
                              radfnbd(nu)
        enddo
      enddo
 
      do k = KS,KE 
        xch2obdinw(:,:,k) = wk(:,:,k)*1.0e-20*(svjinw*rh2os(:,:,k)* &
			    tmpexp(:,:,k) + fvjinw*rfrgn(:,:,k))* &
                            radfnbdinw
        xch2obdwd(:,:,k) = wk(:,:,k)*1.0e-20*(svjwd*rh2os(:,:,k)* &
			   tmpexp(:,:,k) + fvjwd*rfrgn(:,:,k))*  &
                           radfnbdwd
      enddo
 
      totch2obdinw(:,:,1) = 0.0E+00
      totch2obdwd(:,:,1) = 0.0E+00
      do k = KS+1,KE+1
	totch2obdinw(:,:,k) = totch2obdinw(:,:,k-1) +    &
			      xch2obdinw(:,:,k-1)
	totch2obdwd(:,:,k) = totch2obdwd(:,:,k-1) + xch2obdwd(:,:,k-1)
      enddo

!----------------------------------------------------------------------
      deallocate  (totch2obdinw )
      deallocate  (xch2obdinw   )

 
end subroutine optical_path_ckd2p1
 


!################################################################## 

subroutine optical_o3 (atmden, qo3, vv)

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
      allocate (toto3 (ISRAD:IERAD, JSRAD:JERAD, KS:KE      +1))
      allocate (tphio3(ISRAD:IERAD, JSRAD:JERAD, KS:KE      +1))
      allocate (var3  (ISRAD:IERAD, JSRAD:JERAD, KS:KE        ))
      allocate (var4  (ISRAD:IERAD, JSRAD:JERAD, KS:KE        ))

!-----------------------------------------------------------------------
!     compute optical paths for o3, using the diffusivity 
!     approximation 1.66 for the angular integration.  obtain 
!     unweighted values var3 and weighted values  var4.
!     the quantities  0.003 (.003) appearing in the
!     var4 expression are the approximate voigt corrections
!     for o3.
!---------------------------------------------------------------------  
      var3(:,:,KS:KE) = atmden(:,:,KS:KE)*qo3(:,:,KS:KE)*diffac
      var4(:,:,KS:KE) = var3(:,:,KS:KE)*(vv(:,:,KS:KE) + 3.0E-03)

!----------------------------------------------------------------------
!     compute summed optical paths for o3.
!----------------------------------------------------------------------
      toto3 (:,:,KS) = 0.0E+00
      tphio3(:,:,KS) = 0.0E+00
      do k=KS+1,KE+1
        toto3 (:,:,k) = toto3 (:,:,k-1) + var3  (:,:,k-1) 
        tphio3(:,:,k) = tphio3(:,:,k-1) + var4  (:,:,k-1) 
      enddo
!----------------------------------------------------------------------


end  subroutine optical_o3




!#####################################################################

subroutine optical_rbts 

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

      real, dimension(:,:,:), allocatable :: texpsl
      integer     :: k, i

!---------------------------------------------------------------------

      allocate (cntval(ISRAD:IERAD, JSRAD:JERAD,   KS:KE+1     ))
      allocate (totvo2(ISRAD:IERAD, JSRAD:JERAD,   KS:KE+1     ))
      allocate (texpsl(ISRAD:IERAD, JSRAD:JERAD,   KS:KE+1     ))

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
      cntval(:,:,KS:KE) = texpsl(:,:,KS:KE)*rh2o(:,:,KS:KE)*   &
                          var2(:,:,KS:KE)/(rh2o(:,:,KS:KE) +  &
                          rh2oair)

!----------------------------------------------------------------------
!     compute summed optical paths for h2o roberts continuum.
!----------------------------------------------------------------------
      totvo2(:,:,KS) = 0.0E+00
      do k=KS+1,KE+1
        totvo2(:,:,k) = totvo2(:,:,k-1) + cntval(:,:,k-1) 
      enddo
!----------------------------------------------------------------------

      deallocate (texpsl  )


end  subroutine optical_rbts



!####################################################################

subroutine optical_h2o (atmden, vv) 

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
real, dimension (:,:,:), intent(in) ::  atmden, vv 
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
      real, dimension(:,:,:), allocatable :: qh2o, tdif, tdif2
      real, dimension(:,:,:,:), allocatable :: totfh2o, vrfh2o   

!!!  DAN NOTE ::
!!!   totfh2o, vrfh2o not used
!!  Are they still needed ??

!     real, dimension (size(temp,1), size(temp,2), size(temp,3), &
!		  NBTRGE)::  totfh2o

!     real, dimension (size(temp,1), size(temp,2), size(temp,3)-1,  &
!	       NBTRGE)::   vrfh2o

      integer    ::  m, k

!---------------------------------------------------------------------
      allocate (empl1  (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1       ))
      allocate (empl2  (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1       ))
      allocate (totphi (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1       ))
      allocate (var1   (ISRAD:IERAD, JSRAD:JERAD  , KS:KE         ))
      allocate (var2   (ISRAD:IERAD, JSRAD:JERAD  , KS:KE         ))
      allocate (emx1   (ISRAD:IERAD, JSRAD:JERAD                  ))
      allocate (emx2   (ISRAD:IERAD, JSRAD:JERAD                  ))
      allocate (qh2o   (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1       ))
      allocate (tdif   (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1       ))
      allocate (tdif2  (ISRAD:IERAD, JSRAD:JERAD  , KS:KE+1       ))

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
      var1(:,:,KS:KE) = atmden(:,:,KS:KE)*qh2o(:,:,KS:KE)
      var2(:,:,KS:KE) = var1(:,:,KS:KE)*(vv(:,:,KS:KE) + 3.0E-04)

!----------------------------------------------------------------------
!     compute summed optical paths for h2o.
!----------------------------------------------------------------------
      totphi(:,:,KS) = 0.0E+00
      do k=KS+1,KE+1
        totphi(:,:,k) = totphi(:,:,k-1) + var2  (:,:,k-1) 
      enddo

!----------------------------------------------------------------------
!     emx1 is the additional pressure-scaled mass from press(KE) to 
!     pflux(KE).  it is used in nearby layer and emiss calculations.
!     emx2 is the additional pressure-scaled mass from press(KE) to 
!     pflux(KE+1).  it is used in calculations between flux levels k
!     and KE+1.
!----------------------------------------------------------------------
      emx1(:,:) = qh2o(:,:,KE)*press(:,:,KE)*(press(:,:,KE) -   &
                  pflux(:,:,KE))/(grav*pstd)
      emx2(:,:) = qh2o(:,:,KE)*press(:,:,KE)*(pflux(:,:,KE+1) -  &
                  press(:,:,KE))/(grav*pstd)

!----------------------------------------------------------------------
!     empl is the pressure scaled mass from pflux(k) to press(k) or to 
!     press(k+1).
!----------------------------------------------------------------------
      empl1(:,:,KS) = var2(:,:,KE)
      empl1(:,:,KS+1:KE+1) =   &
                 qh2o(:,:,KS:KE)*pflux(:,:,KS+1:KE+1)*   &
                 (pflux(:,:,KS+1:KE+1) - press(:,:,KS:KE))/(grav*pstd)
      empl2(:,:,KS+1:KE) =    &
                 qh2o(:,:,KS+1:KE)*pflux(:,:,KS+1:KE)*   &
                 (press(:,:,KS+1:KE) - pflux(:,:,KS+1:KE))/(grav*pstd)
      empl2(:,:,KE+1) = empl2(:,:,KE) 

!---------------------------------------------------------------------
  if (Lw_control%do_ch4_n2o) then
    allocate ( empl1f (ISRAD:IERAD , JSRAD:JERAD , KS:KE+1,  NBTRGE ) )
    allocate ( empl2f (ISRAD:IERAD , JSRAD:JERAD , KS:KE+1,  NBTRGE ) )
    allocate ( tphfh2o(ISRAD:IERAD , JSRAD:JERAD , KS:KE+1,  NBTRGE ) ) 
    allocate ( vrpfh2o(ISRAD:IERAD , JSRAD:JERAD , KS:KE+1,  NBTRGE ) )
    allocate ( emx1f  (ISRAD:IERAD , JSRAD:JERAD ,           NBTRGE ) )
    allocate ( emx2f  (ISRAD:IERAD , JSRAD:JERAD ,           NBTRGE ) )

!----------------------------------------------------------------------
!   compute h2o optical paths for use in the 1200-1400 cm-1 range if
!   temperature dependence of line intensities is accounted for.
!----------------------------------------------------------------------
    tdif(:,:,KS:KE) = temp(:,:,KS:KE)-2.5E+02

    do m=1,NBTRGE
      vrpfh2o(:,:,KS:KE,m) = var2(:,:,KS:KE) *    &
                             EXP(csfah2o(1,m)*(tdif(:,:,KS:KE)) +   &
                                 csfah2o(2,m)*(tdif(:,:,KS:KE))**2 )
    enddo

    do m=1,NBTRGE
      tphfh2o(:,:,KS,m) = 0.0E+00
      do k=KS+1,KE+1
        tphfh2o(:,:,k,m) = tphfh2o(:,:,k-1,m) + vrpfh2o(:,:,k-1,m)
      enddo
    enddo

    tdif2(:,:,KS+1:KE+1) = tpl2(:,:,KS+1:KE+1)-2.5E+02
    tdif (:,:,KS+1:KE+1) = tpl1(:,:,KS+1:KE+1)-2.5E+02

!---------------------------------------------------------------------
!   compute this additional mass, for use in the 1200-1400 cm-1 range,
!   if temperature dependence of line intensities is accounted for.
!--------------------------------------------------------------------
    do m=1,NBTRGE
      emx1f(:,:,m) = emx1(:,:) *    &
                     EXP(csfah2o(1,m)*(tdif2(:,:,KE+1)        ) +   &
                         csfah2o(2,m)*(tdif2(:,:,KE+1)        )**2 )
      emx2f(:,:,m) = emx2(:,:) *    &
                     EXP(csfah2o(1,m)*(tdif (:,:,KE+1)        ) +  &
                         csfah2o(2,m)*(tdif (:,:,KE+1)        )**2 )
    enddo

!----------------------------------------------------------------------
!   compute this additional mass, for use in the 1200-1400 cm-1 range,
!   if temperature dependence of line intensities is accounted for.
!----------------------------------------------------------------------
    do m=1,NBTRGE
      empl1f(:,:,KS+1:KE+1,m) = empl1(:,:,KS+1:KE+1) *   &
                          EXP(csfah2o(1,m)*(tdif(:,:,KS+1:KE+1)) +    &
                              csfah2o(2,m)*(tdif(:,:,KS+1:KE+1))**2 )
      empl2f(:,:,KS+1:KE,m) = empl2(:,:,KS+1:KE) *   &
                          EXP(csfah2o(1,m)*(tdif2(:,:,KS+1:KE)) +   &
                              csfah2o(2,m)*(tdif2(:,:,KS+1:KE))**2 )
      empl1f(:,:,KS ,m) = vrpfh2o(:,:,KE,m)
      empl2f(:,:,KE+1,m) = empl2f(:,:,KE,m)
    enddo
  endif

!---------------------------------------------------------------------
      deallocate (tdif2 )
      deallocate (tdif  )
      deallocate (qh2o  )



end subroutine optical_h2o



!####################################################################



                   end module optical_path_mod


