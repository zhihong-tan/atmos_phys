
                 module microphys_rad_mod

use utilities_mod,       only: open_file, file_exist,   &
                               check_nml_error, error_mesg,   &
                               print_version_number, FATAL, NOTE, &
			       WARNING, get_my_pe, close_file, &
			       get_domain_decomp
use rad_step_setup_mod,  only: temp, press, jabs, iabs, deltaz,&
			       IMINP, IMAXP, JMINP, JMAXP, &
                               ISRAD, IERAD, JSRAD, JERAD, & 
                               KSRAD, KERAD
use rad_utilities_mod,   only: longwave_control_type, Lw_control, &
                               shortwave_control_type, Sw_control,&
			       radiation_control_type, Rad_control
use longwave_setup_mod,  only: longwave_parameter_type, &    
                               Lw_parameters
use radiation_diag_mod,  only: radiag_from_cloudrad
use longwave_params_mod, only: NBLW 
use esfsw_parameters_mod,only: nbands, get_solarfluxes, &
			       TOT_WVNUMS
use constants_new_mod,   only: diffac

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!	       module to produce radiative cloud properties 
!                 based upon microphysical properties
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* -----------------------

!   character(len=5), parameter  ::  version_number = 'v0.09'
    character(len=128)  :: version =  '$Id: microphys_rad.F90,v 1.2 2001/07/05 17:32:05 fms Exp $'
    character(len=128)  :: tag     =  '$Name: eugene $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
	  microphys_rad_init, microphys_rad_driver,   &
	  thickavg, thinavg


private         &
	  cloud_lwpar, furainlw, fusnowlw, cliqlw, el, cloudpar, &
	  slingo, savijarvi, fu, snowsw, snowlw, icesolar, &
	  cloud_lwem_oneband


!---------------------------------------------------------------------
!-------- namelist  ---------


integer                      :: n_prsc_clds=3
real                         :: wtr_cld_reff=10.
real                         :: ice_cld_reff=50.
real                         :: rain_reff=250.
logical                      :: using_fu=.false. ! MUST be .false.
character(len=10)            :: lwem_form=' '



namelist /microphys_rad_nml /     &
	        		  n_prsc_clds, &
				  wtr_cld_reff, &
				  ice_cld_reff, &
				  using_fu,   &
				  rain_reff, &
				  lwem_form


!----------------------------------------------------------------------
!----  public data -------

!----------------------------------------------------------------------
!----  private data -------

!----------------------------------------------------------------------
!     parameters for longwave cloud-radiation parameterizations
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!     N_EMISS_BDS = number of  infrared frequency bands over which 
!                   frequency-dependent emissivities are computed 
!                   using appropriate parameterizations. whether
!                   or not these emissivities are actually employed
!                   by the model depends on the value of NLWCLDB
!                   which is allowed to have the value of N_EMISS_BDS
!                   or of unity (when these emissivities are not
!                   used).
!    NLWCLDB      = the actual number of frequency bands for which lw
!                   emissitivies are defined. 
!----------------------------------------------------------------------
 
integer, parameter       :: N_EMISS_BDS = 7
integer                  :: NLWCLDB 

real,    dimension (N_EMISS_BDS)     :: cldbandlo, cldbandhi
integer, dimension (N_EMISS_BDS + 1) :: istartcldband, iendcldband
integer, dimension (N_EMISS_BDS + 1) :: nivl1lwicecld, nivl2lwicecld
real,    dimension (N_EMISS_BDS)     :: planckcldband

!----------------------------------------------------------------------
!   wavenumber ranges with separate cloud emissivity values in the
!   infrared parameterization. these may be changed only by the
!   keeper of the radiation code.
!----------------------------------------------------------------------

data cldbandlo /                    &
             0.0, 560.0, 800.0, 900.0, 990.0, 1070.0, 1200.0/
data cldbandhi /                    &
           560.0, 800.0, 900.0, 990.0, 1070.0, 1200.0, 1400.0/

!----------------------------------------------------------------------
!   note: the cloud properties for wavelengths beyond 1400 wavenumbers
!   are included in the results for the first band, ie, that band
!   actually is 0-560, 1400-2200 cm-1. thus the indices include an
!   extra band.
!----------------------------------------------------------------------

data istartcldband /                &
               1,   57,   81,    91,   100,    108,   121,   141/
data iendcldband /                  &
              56,   80,   90,    99,   107,    120,   140,   220/

!----------------------------------------------------------------------
!     add parameters for Fu and Liou lw snow water parameterization
!
!     NBFL = number of frequency bands in the lw snow water para-
!     meterization.corresponds to bands 7-18 in Fu and Liou (Table 2).
!     NBA, NBB, NBC = number of terms in parameterization for series
!    expansions (not counting the 0th power term) for ai, bi, ci
!----------------------------------------------------------------------
 
integer, parameter  :: NBFL= 12
integer, parameter  :: NBA = 3
integer, parameter  :: NBB = 3
integer, parameter  :: NBC = 3
integer, parameter  :: NBD = 3

integer, dimension (NBFL)        :: iendfubands
real,    dimension (NBFL)        :: endfubands
 
!----------------------------------------------------------------------
!     wavenumber ranges  for Fu-Liou ice crystal parameterizations
!       these apply to the ice crystal (El), cloud rain (Furainlw)
!       and cloud snow (Fusnowlw) parameterizations.
!     note: the cloud liquid drop parameterization (Cliqlw) is 
!       frequency-independent,
!----------------------------------------------------------------------
 
!----------------------------------------------------------------------
!     endfubands = high wavenumber boundary of wavenumber bands used
!             in Fu-Liou parameterization. the order of increasing
!             wavenumber, thus bands 1-12 correspond to Fu's
!             bands 18-7.
!     iendfubands = index of model 10 cm-1 bands corresponding to
!              the value of endfubands. computed in the code.
!----------------------------------------------------------------------
 
data endfubands/                  &
                280,  400,  540,  670,  800,  980,  1100,  1250,    &
                1400,   1700,   1900,   2200/

!----------------------------------------------------------------------
!    weighting factors relating fu lw bands to model lw frequency
!    bands
!----------------------------------------------------------------------

real,    dimension (N_EMISS_BDS, NBFL)     :: fulwwts
real,    dimension (N_EMISS_BDS + 1, NBFL) :: planckivlicecld
 
!----------------------------------------------------------------------
!     parameters for shortwave cloud-radiation parameterizations
!----------------------------------------------------------------------
 
!----------------------------------------------------------------------
! NICECLDIVLS  = the number of scattering spectral intervals for the    
!                ice crystals                                           
!                                                                       
! NLIQCLDIVLS  = the number of scattering spectral intervals for the    
!                cloud drops                                            
!                                                                       
! NRAINCLDIVLS = the number of scattering spectral intervals for the    
!                rain drops                                             
!                                                                       
! NSNOWCLDIVLS = the number of scattering spectral intervals for snow   
!----------------------------------------------------------------------
                                                                        
integer, parameter         ::   NICECLDIVLS = 25 
integer, parameter         ::   NICESOLARCLDIVLS = 6
integer, parameter         ::   NLIQCLDIVLS = 24 
integer, parameter         ::   NRAINCLDIVLS = 4 
integer, parameter         ::   NSNOWCLDIVLS = 6 
 
integer, dimension (NICECLDIVLS)       :: endicecldwvn
integer, dimension (NICESOLARCLDIVLS)  :: endicesolcldwvn
integer, dimension (NLIQCLDIVLS)       :: endliqcldwvn
integer, dimension (NRAINCLDIVLS)      :: endraincldwvn
integer, dimension (NSNOWCLDIVLS)      :: endsnowcldwvn

!----------------------------------------------------------------------c
!  define the spectral limits for drop, rain, ice and snow single      c
!  scattering properties in shortwave frequency ranges.                c
!                                                                      c
!  note: the last wavenumber value must be the same as the             c
!        parameterization's last band limit.                           
!----------------------------------------------------------------------
 
!----------------------------------------------------------------------c
! wavenumber limits for slingo cloud drop intervals                    c
!----------------------------------------------------------------------c
 
data endliqcldwvn / 2924, 3437, 4202, 4695, 6098, 6536, 7813,   &
                          8404, 9091,10000,11494,12821,13333,14493,   &
                         15625,17544,19231,20833,22727,25000,27778,   &
                         30303,33333,57600 /
 
!----------------------------------------------------------------------c
! wavenumber limits for Savijarvi rain drop intervals                  c
!----------------------------------------------------------------------c
 
data endraincldwvn / 4202,8403,14493,57600 /
 
!---------------------------------------------------------------------- 
! wavenumber limits for icesolar ice crystal intervals                
!---------------------------------------------------------------------- 
 
data endicesolcldwvn / 2857, 4000, 5263, 7692, 14493, 57600 /
 
!----------------------------------------------------------------------c
! wavenumber limits for fu ice crystal intervals                       c
!----------------------------------------------------------------------c
 
data endicecldwvn / 2000, 2924, 3437, 4202, 4695, 6098, 6536,   &
                          7092, 8404, 9091,10000,11494,12821,13333,   &
                         14493,15625,17544,19231,20833,22727,25000,   &
                         27778,30303,33333,57600 /
 
!----------------------------------------------------------------------c
! wavenumber limits for the Fu snow intervals                          c
!----------------------------------------------------------------------c
 
data endsnowcldwvn / 2857,4000,5263,7692,14493,57600 /
 
!----------------------------------------------------------------------
integer, dimension(:), allocatable      :: nivl1liqcld,   &
				           nivl1icecld,   &
				           nivl1icesolcld,   &
				           nivl1raincld,  &
				           nivl1snowcld,  &
				           nivl2liqcld,   &
				           nivl2icecld,   &
				           nivl2icesolcld,   &
				           nivl2raincld,  &
				           nivl2snowcld 

real, dimension(:,:), allocatable       :: solivlicecld,   &
                                           solivlicesolcld, &
                                           solivlliqcld, &
                                           solivlraincld , &
                                           solivlsnowcld


real,    dimension (:), allocatable     :: solflxband


!--------------------------------------------------------------------
!     these variables define the boundaries (in sigma coordinates) 
!     between high and middle and middle and low clouds. 
!-----------------------------------------------------------------

real, dimension(:), allocatable   ::  cldhm_abs, cldml_abs

!----------------------------------------------------------------------
! define the liquid water path for high, middle and low clouds in units
! of grams per meter**2.  
!----------------------------------------------------------------------
real   ::  lwpath_hi  = 6.313929
real   ::  lwpath_mid = 18.94179
real   ::  lwpath_low = 75.76714

!-------------------------------------------------------------------
character(len=10)    :: swform
logical              :: do_init=.false.
integer              :: xdom(4), y(4), jdf 


!--------------------------------------------------------------------
!--------------------------------------------------------------------




                 contains 




subroutine microphys_rad_init (cldhm_abs_in, cldml_abs_in)

!------------------------------------------------------------------
real, dimension(:), intent(in),optional :: cldhm_abs_in, cldml_abs_in
!--------------------------------------------------------------------

!----------------------------------------------------------------------
! local variables:                                                  
!---------------------------------------------------------------------

      integer           :: unit, ierr, io
      integer           :: k, j, i, n
      integer           :: ni
      integer           :: iounit, nf, nw1, nw2, nintsolar, &
  	                   nivl1, nivl2, nivl3
      integer           :: ib, nw, nivl, nband  
      integer           :: nivl1, nivl2, nivl3, nivl4, nivl5
      real              :: sumsol1, sumsol2, sumsol3, sumsol
      real              :: sumsol4, sumsol5
      real              :: del, xtemv, sumplanck 

      real,    dimension (NBLW)            :: c1, centnb, sc, src1nb,  &
                                              x, x1
      integer, dimension(:), allocatable   :: endwvnbands
      real   , dimension(:), allocatable   :: solflxband_in
      real   , dimension(:), allocatable   :: solarfluxtoa        

!---------------------------------------------------------------------
!-----  read namelist  ------
  
      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=microphys_rad_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'microphys_rad_nml')
        enddo
10      call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'microphys_rad', version_number)
      if (get_my_pe() == 0)  then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=microphys_rad_nml)
      endif
      call close_file (unit)

 !--------------------------------------------------------------------
 !   perform consistency checks between lwem_form and desired lwcld 
 !   emissivity.
 !--------------------------------------------------------------------
      if (Lw_control%do_lwcldemiss)    then
        if (trim(lwem_form) == 'fuliou') then
        else if (trim(lwem_form) == 'ebertcurry') then
          call error_mesg('microphys_rad_driver',  &
              'ebert-curry not yet implemented for multi-bands', FATAL)
        else
          call error_mesg('microphys_rad_driver',  &
           'incorrect specification of lwem_form for multi-band', FATAL)
        endif
      else
        if (trim(lwem_form) == 'fuliou') then
          call error_mesg('microphys_rad_driver',  &
          'fu parameterization implemented only for multi-band', FATAL)
        else if (trim(lwem_form) == 'ebertcurry') then
        else
           call error_mesg('microphys_rad_driver',  &
         'incorrect specification of lwem_form for single band', FATAL)
        endif
      endif

!--------------------------------------------------------------------
!  retrieve module variables that come from other modules. save input 
!  arguments as module variables. define the number of cloud emissivity 
!  bands actually being used in this experiment.
!--------------------------------------------------------------------

     if (present (cldhm_abs_in)) then
      allocate (cldhm_abs (1:jdf) )
      cldhm_abs = cldhm_abs_in
     endif
     if (present (cldml_abs_in)) then
      allocate (cldml_abs (1:jdf) )
      cldml_abs = cldml_abs_in
     endif

!--------------------------------------------------------------------
!   allocate module variables.
!--------------------------------------------------------------------
      allocate ( nivl1liqcld (nbands),    &
                 nivl1icecld (nbands),   &
                 nivl1icesolcld(nbands),   &
                 nivl1raincld(nbands),  &
                 nivl1snowcld(nbands),  &
 	         nivl2liqcld(nbands),   &
 	         nivl2icecld(nbands),   &
 	         nivl2icesolcld(nbands),   &
 	         nivl2raincld(nbands),   & 
                 nivl2snowcld(nbands) ) 

      allocate ( solivlicecld(nbands, nicecldivls),  &
                 solivlicesolcld(nbands, nicesolarcldivls), &
                 solivlliqcld(nbands, nliqcldivls), &
                 solivlraincld(nbands, nraincldivls), &
                 solivlsnowcld(nbands, nsnowcldivls) )

      allocate ( solflxband(nbands) )

      NLWCLDB = Lw_parameters%NLWCLDB
!------------------------------------------------------------------

      if (Lw_control%do_lwcldemiss) then
 
!--------------------------------------------------------------------
!     compute band-averaged coefficients for microphysics-to-radiation
!     interface for cloud species in infrared frequency ranges.
!     the actual extinction coefficients (to become emissivities)
!      (and other coefficients) are calculated as time-dependent
!     quantities.
!
!     at present, the species included are:
!     1)   cloud drops
!     2)   ice crystals
!     3)   cloud rain
!     4)   cloud snow
!--------------------------------------------------------------------
 
 
!--------------------------------------------------------------------
!    calculation for ice crystals (Fu-Liou)
!--------------------------------------------------------------------
 
        iendfubands(:) = INT((endfubands(:)+0.01)/10.0)

!--------------------------------------------------------------------
!      compute weighting function. according to Fu and Liou, this
!      should be the Planck function at -40C.
!--------------------------------------------------------------------
        do n=1,NBLW 
          del  = 10.0E+00
	  xtemv = 233.15
	  centnb(n) = 5.0 + (n - 1)*del
          c1(n)     = (3.7412E-05)*centnb(n)**3
          x(n)      = 1.4387E+00*centnb(n)/xtemv
          x1(n)     = EXP(x(n))
          sc(n)     = c1(n)/(x1(n) - 1.0E+00)
          src1nb(n) = del*sc(n)
        enddo
 
!--------------------------------------------------------------------
!      compute summed weighting function over the (N_EMISS_BDS) cloud
!      bands
!--------------------------------------------------------------------

        planckcldband(:) = 0.0E+00
        do n = 1,N_EMISS_BDS
  	  do ib = istartcldband(n),iendcldband(n)
	    planckcldband(n) = planckcldband(n) + src1nb(ib)
	  enddo
        enddo

!--------------------------------------------------------------------
!     add contribution of 1400-2200 cm-1 region to first band
!--------------------------------------------------------------------
        do ib = istartcldband(N_EMISS_BDS+1),iendcldband(N_EMISS_BDS+1)
          planckcldband(1) = planckcldband(1) + src1nb(ib)
        enddo
 
        nivl = 1
        sumplanck = 0.0
        nband = 1
        planckivlicecld(:,:) = 0.0
        nivl1lwicecld(1) = 1
 
        do nw = 1,NBLW
  	  sumplanck = sumplanck + src1nb(nw)
          if ( nw.eq.iendfubands(nivl) ) then
            planckivlicecld(nband,nivl) = sumplanck
            sumplanck = 0.0
          end if
          if ( nw.eq.iendcldband(nband) ) then
            if ( nw.ne.iendfubands(nivl) ) then
              planckivlicecld(nband,nivl) = sumplanck 
              sumplanck = 0.0
            end if
            nivl2lwicecld(nband) = nivl
            nband = nband + 1
            if ( nband.le.N_EMISS_BDS+1 ) then
              if ( nw.eq.iendfubands(nivl) ) then
                nivl1lwicecld(nband) = nivl + 1
              else
                nivl1lwicecld(nband) = nivl
              end if
            end if
          end if
          if ( nw .eq. iendfubands(nivl) ) nivl = nivl + 1
          if ( nw .ge. iendcldband(N_EMISS_BDS+1) ) then
	    exit
	  endif
        end do
!--------------------------------------------------------------------
!     compute planck-weighted band weights for Fu lw microphysics
!     calculations
!--------------------------------------------------------------------
        fulwwts(:,:) = 0.0E+00
        do n=1,N_EMISS_BDS
          do ni=nivl1lwicecld(n),nivl2lwicecld(n)
	    fulwwts(n,ni) = planckivlicecld(n,ni)/planckcldband(n)
	  enddo
        enddo
!--------------------------------------------------------------------
!     add band (N_EMISS_BDS+1) to band 1 weights
!--------------------------------------------------------------------
        do ni=nivl1lwicecld(N_EMISS_BDS+1),nivl2lwicecld(N_EMISS_BDS+1)
          fulwwts(1,ni) = planckivlicecld(N_EMISS_BDS+1,ni)/   &
                          planckcldband(1)
        enddo
      endif

!--------------------------------------------------------------------
!    compute shortwave microphysics bands and weights
!--------------------------------------------------------------------

      swform = Sw_control%sw_form
      if (trim(swform) == 'esfsw99') then

	allocate (solarfluxtoa(TOT_WVNUMS) )
	allocate (solflxband_in(nbands) )
	allocate (endwvnbands(0:nbands) )

        call get_solarfluxes (solarfluxtoa, solflxband_in, endwvnbands)
	solflxband = solflxband_in

 
!----------------------------------------------------------------------c
! define the solar weights and interval counters that are used to      c
! determine the single-scattering properties for the parameterization  c
! band spectral intervals, from the specified spectral intervals for   c
! drops, ice particles and aerosols.                                   c
!----------------------------------------------------------------------c
 
        nivl1 = 1
        nivl2 = 1
        nivl3 = 1
        nivl4 = 1
        nivl5 = 1
        sumsol1 = 0.0
        sumsol2 = 0.0
        sumsol3 = 0.0
        sumsol4 = 0.0
        sumsol5 = 0.0
        nband = 1
        solivlliqcld(:,:) = 0.0
        solivlicecld(:,:) = 0.0
        solivlicesolcld(:,:) = 0.0
        solivlraincld(:,:) = 0.0
        solivlsnowcld(:,:) = 0.0
        nivl1liqcld(1) = 1
        nivl1icecld(1) = 1
        nivl1icesolcld(1) = 1
        nivl1raincld(1) = 1
        nivl1snowcld(1) = 1
 
        do nw = 1,endwvnbands(nbands)
	  sumsol1 = sumsol1 + solarfluxtoa(nw) 
	  sumsol2 = sumsol2 + solarfluxtoa(nw) 
	  sumsol3 = sumsol3 + solarfluxtoa(nw) 
	  sumsol4 = sumsol4 + solarfluxtoa(nw) 
	  sumsol5 = sumsol5 + solarfluxtoa(nw) 
 
          if ( nw.eq.endliqcldwvn(nivl1) ) then
            solivlliqcld(nband,nivl1) = sumsol1
            sumsol1 = 0.0
          end if
          if ( nw.eq.endicecldwvn(nivl2) ) then
            solivlicecld(nband,nivl2) = sumsol2
            sumsol2 = 0.0
          end if
          if ( nw.eq.endraincldwvn(nivl3) ) then
            solivlraincld(nband,nivl3) = sumsol3
            sumsol3 = 0.0
          end if
          if ( nw.eq.endsnowcldwvn(nivl4) ) then
            solivlsnowcld(nband,nivl4) = sumsol4
            sumsol4 = 0.0
          end if
          if ( nw.eq.endicesolcldwvn(nivl5) ) then
            solivlicesolcld(nband,nivl5) = sumsol5
            sumsol5 = 0.0
          end if
 

          if ( nw.eq.endwvnbands(nband) ) then
 
            if ( nw.ne.endliqcldwvn(nivl1) ) then
              solivlliqcld(nband,nivl1) = sumsol1 
              sumsol1 = 0.0
            end if
            if ( nw.ne.endicecldwvn(nivl2) ) then
              solivlicecld(nband,nivl2) = sumsol2 
              sumsol2 = 0.0
            end if
            if ( nw.ne.endraincldwvn(nivl3) ) then
              solivlraincld(nband,nivl3) = sumsol3 
              sumsol3 = 0.0
            end if
            if ( nw.ne.endsnowcldwvn(nivl4) ) then
              solivlsnowcld(nband,nivl4) = sumsol4 
              sumsol4 = 0.0
            end if
            if ( nw.ne.endicesolcldwvn(nivl5) ) then
              solivlicesolcld(nband,nivl5) = sumsol5 
              sumsol5 = 0.0
            end if
 
            nivl2liqcld(nband) = nivl1
            nivl2icecld(nband) = nivl2
            nivl2raincld(nband) = nivl3
            nivl2snowcld(nband) = nivl4
            nivl2icesolcld(nband) = nivl5
 
            nband = nband + 1
 
            if ( nband.le.nbands ) then
 
              if ( nw.eq.endliqcldwvn(nivl1) ) then
                nivl1liqcld(nband) = nivl1 + 1
              else
                nivl1liqcld(nband) = nivl1
              end if
              if ( nw.eq.endicecldwvn(nivl2) ) then
                nivl1icecld(nband) = nivl2 + 1
              else
                nivl1icecld(nband) = nivl2
              end if
              if ( nw.eq.endraincldwvn(nivl3) ) then
                nivl1raincld(nband) = nivl3 + 1
              else
                nivl1raincld(nband) = nivl3
              end if
              if ( nw.eq.endsnowcldwvn(nivl4) ) then
                nivl1snowcld(nband) = nivl4 + 1
              else
                nivl1snowcld(nband) = nivl4
              end if
              if ( nw.eq.endicesolcldwvn(nivl5) ) then
                nivl1icesolcld(nband) = nivl5 + 1
              else
                nivl1icesolcld(nband) = nivl5
              end if
 
            end if
 
          end if
 
          if ( nw.eq.endliqcldwvn(nivl1) ) nivl1 = nivl1 + 1
          if ( nw.eq.endicecldwvn(nivl2) ) nivl2 = nivl2 + 1
          if ( nw.eq.endraincldwvn(nivl3) ) nivl3 = nivl3 + 1
          if ( nw.eq.endsnowcldwvn(nivl4) ) nivl4 = nivl4 + 1
          if ( nw.eq.endicesolcldwvn(nivl5) ) nivl5 = nivl5 + 1
 
        end do
 
        deallocate  (solarfluxtoa)
        deallocate  (solflxband_in)
        deallocate  (endwvnbands)
 

      endif

!-----------------------------------------------------------------
!  obtain global and subdomain dimensions for later use.
!----------------------------------------------------------------
      call get_domain_decomp (xdom, y)
      jdf = y(4) - y(3) + 1

!-----------------------------------------------------------------
!  set flag to indicate that the routine has been completed.
!----------------------------------------------------------------

      do_init = .true.

!-----------------------------------------------------------------

 101  format( 12f10.4 )
 102  format( 32i4 )
 103  format( 20i6 )
 104  format( 12f10.2 )
 105  format( 1p,16e8.1 )
 106  format( 1p,3e16.6,l16 )
 107  format( i5,1p,e14.5 )


 !-----------------------------------------------------------------


end subroutine microphys_rad_init




!#################################################################

subroutine microphys_rad_driver (camtsw, emmxolw, emrndlw, &
			         cldext, cldsct, cldasymm, &
                                 conc_drop_in, conc_ice_in,   &
				 conc_rain_in, conc_snow_in,    &
	                         size_drop_in, size_ice_in,   &
				 size_rain_in, size_snow_in)

!---------------------------------------------------------------------
real, dimension(:,:,:),   intent(inout) :: camtsw
real, dimension(:,:,:,:), intent(inout) :: emmxolw, emrndlw
real, dimension(:,:,:,:), intent(inout), optional  ::  &
					   cldext, cldsct, cldasymm
real, dimension (:,:,:), intent(inout), optional  ::  &
                                           size_drop_in, size_ice_in,  &
                                           size_rain_in, size_snow_in
real, dimension (:,:,:), intent(inout), optional  ::  &
                                           conc_drop_in, conc_ice_in,  &
                                           conc_rain_in, conc_snow_in

!----------------------------------------------------------------------
! local variables:                                                  
!---------------------------------------------------------------------

     real, dimension (:,:,:), allocatable  :: lwpath, iwpath
     real, dimension (:,:,:), allocatable  :: conc_drop, conc_ice,  &
                                              conc_rain, conc_snow
     real, dimension (:,:,:), allocatable  :: size_drop, size_ice,  &
                                              size_rain, size_snow
     integer, dimension (:,:), allocatable :: nhi_clouds,    &
					      nmid_clouds, &
                                              nlow_clouds
     real, dimension (:,:), allocatable    :: press_hm, press_ml
     real, dimension (:,:,:), allocatable  :: conc


     integer       :: k, j, i

!-------------------------------------------------------------------
!  be sure module has been initialized.
!--------------------------------------------------------------------
     if (.not. do_init) then
       call error_mesg ('microphys_rad_driver',  &
            'initialization routine of this module was never called', &
							    FATAL)
     endif

!-------------------------------------------------------------------- 
!  allocate arrays for size and concentration of microphys species.
!-------------------------------------------------------------------- 
     allocate (conc_drop   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
     allocate (conc_ice    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
     allocate (conc_rain   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
     allocate (conc_snow   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
     allocate (size_drop   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
     allocate (size_ice    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
     allocate (size_rain   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
     allocate (size_snow   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
     allocate ( lwpath     (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
     allocate ( iwpath     (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )

!----------------------------------------------------------------------
!   define sizes and concentrations of microphysical species. use 
!   supplied default value if no input argument is provided. note 
!   default drop and rain values are radii, multiply by 2 to produce 
!   diameter as desired in this modules microphysics routines. assume
!   precip particles are non-existent, if values are not supplied via 
!   input.
!----------------------------------------------------------------------

     if ( present(size_drop_in) )  then
       size_drop(:,:,:) = size_drop_in(:,:,:)
     else
       size_drop(:,:,:) = 2.0*wtr_cld_reff
     endif

     if ( present(size_ice_in) )  then
       size_ice(:,:,:) = size_ice_in(:,:,:)
     else
       size_ice(:,:,:) = ice_cld_reff
     endif

     if ( present(size_rain_in) )  then
       size_rain(:,:,:) = size_rain_in(:,:,:)
     else
       size_rain(:,:,:) = 2.0*rain_reff
     endif

     if ( present(conc_rain_in) )  then
       conc_rain(:,:,:) = conc_rain_in(:,:,:)
     else
       conc_rain = 0.0
     endif

     if ( present(conc_snow_in) )  then
       conc_snow(:,:,:) = conc_snow_in(:,:,:)
     else
       conc_snow = 0.0
     endif

     if ( present(conc_ice_in) .and. present(conc_drop_in) )  then
       conc_ice(:,:,:) = conc_ice_in(:,:,:)
       conc_drop(:,:,:) = conc_drop_in(:,:,:)
     else  if (present(conc_ice_in) ) then
       call error_mesg ('microphys_rad_driver',  &
	 'ice concentration present but no drop concentration', FATAL)
     else  if (present(conc_drop_in) ) then
       call error_mesg ('microphys_rad_driver',  &
	 'drop concentration present but no ice concentration', FATAL)
     else
!---------------------------------------------------------------------
!  if concentrations of drops and ice are not input from some
!  microphysical model (external source) then they will be calculated
!  here, based on some prescribed values. allocate arrays which will
!  be used in defining these concentrations.
!---------------------------------------------------------------------
     
       allocate (nhi_clouds  (ISRAD:IERAD, JSRAD:JERAD) )
       allocate (nmid_clouds (ISRAD:IERAD, JSRAD:JERAD) )
       allocate (nlow_clouds (ISRAD:IERAD, JSRAD:JERAD) )
       allocate (press_hm    (ISRAD:IERAD, JSRAD:JERAD) )
       allocate (press_ml    (ISRAD:IERAD, JSRAD:JERAD) )
       allocate (conc        (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )

!----------------------------------------------------------------------
! assume that the water path ( = conc.* path(km) ) is preset at fixed 
! values (lwpath_hi, _mid, _low) for "high", "mid", "low" clouds. the 
! lwpath in each cloud layer within "hi", "mid" "low" pressure specif-
! ication is that lwpath_... divided by the number of clouds actually 
! in that pressure interval. the unit of liquid water path is grams per 
! meter**2. 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! define the number of high, middle, low clouds according to
! Wetherald's criterion
!----------------------------------------------------------------------
    if (.not. allocated (cldhm_abs) .or.   &
        .not. allocated(cldml_abs)) then
           call error_mesg ('microphys_rad_driver', &
	 'must specify h-m-l cloud boundaries if not supplying micro&
	 &physical properties', FATAL)
     endif
       do j=JSRAD,JERAD
         do i=ISRAD,IERAD
           nhi_clouds(i,j) = 0
           nmid_clouds(i,j) = 0
           nlow_clouds(i,j) = 0
           press_hm(i,j) = cldhm_abs(jabs(j))*press(i,j,KERAD+1)
           press_ml(i,j) = cldml_abs(jabs(j))*press(i,j,KERAD+1)
           do k=KSRAD,KERAD
             if (camtsw(i,j,k) > 0.0) then
               if (press(i,j,k) .LE. press_hm(i,j)) then     
                 nhi_clouds(i,j) = nhi_clouds(i,j) + 1
	       endif
	       if (press(i,j,k) .GT. press_hm(i,j) .AND.      &
                   press(i,j,k) .LE. press_ml(i,j)) then
	         nmid_clouds(i,j) = nmid_clouds(i,j) + 1
	       endif
	       if (press(i,j,k) .GT. press_ml(i,j) ) then       
	         nlow_clouds(i,j) = nlow_clouds(i,j) + 1
	       endif
             endif
           enddo
         enddo
       enddo
      
!----------------------------------------------------------------------
!   compute conc in each layer as (water path / layer geometric path)
!----------------------------------------------------------------------
       conc(:,:,:) = 0.0E+00
       do j=JSRAD,JERAD
 	 do i=ISRAD,IERAD
           do k=KSRAD,KERAD
             if (camtsw(i,j,k) > 0.0) then
               if (press(i,j,k) .LE. press_hm(i,j) ) then
	 	 conc(i,j,k) = lwpath_hi     /                 &
                               (nhi_clouds(i,j)*deltaz(i,j,k))
	       endif
	       if (press(i,j,k) .GT. press_hm(i,j) .AND.          &
                   press(i,j,k) .LE. press_ml(i,j) )  then          
		 conc(i,j,k) = lwpath_mid    /                   &
                               (nmid_clouds(i,j)*deltaz(i,j,k))
	       endif
	       if (press(i,j,k) .GT. press_ml(i,j)) then
		 conc(i,j,k) = lwpath_low    /                   &
                               (nlow_clouds(i,j)*deltaz(i,j,k))
	       endif
             endif
	   enddo
	 enddo
       enddo

!----------------------------------------------------------------------
!   split conc into conc_ice and conc_drop, depending on temperature
!   criterion (T < 273.16); add default size_ice, size_drop as
!   needed.
!----------------------------------------------------------------------
       do k = KSRAD,KERAD
         do j = JSRAD,JERAD
	   do i = ISRAD,IERAD
	     if (temp(i,j,k) .LT. 273.16) then
	       conc_ice(i,j,k) = conc(i,j,k)
	       conc_drop(i,j,k) = 0.0E+00
	     else
	       conc_ice(i,j,k) = 0.0E+00
	       conc_drop(i,j,k) = conc(i,j,k)
	     endif
	   enddo
	 enddo
       enddo
     endif


!---------------------------------------------------------------------
!   define liquid water path (lwpath) and ice water path (iwpath) from
!   conc_drop(ice) and deltaz. use units kg/m**2 accounting for
!   factor of 1/1000. however for diagnostic output we will desire
!   lwpath, iwpath in units of g/m**2
!---------------------------------------------------------------------
 
     lwpath(:,:,:) = 1.0E-03*conc_drop(:,:,:)*deltaz(:,:,:)
    iwpath(:,:,:) = 1.0E-03*conc_ice(:,:,:)*deltaz(:,:,:)

!---------------------------------------------------------------------
!   define microphysically-based lw cloud emissivities. use
!   cloud_lwpar to compute multi-band emissivities based on fu
!   parameterizations for drop, ice, snow and rain. use
!   cloud_lwem_oneband to compute a single value for the lw emissivity
!   (including effects of drops and ice) based on ebert and curry.
!---------------------------------------------------------------------
    if (Lw_control%do_lwcldemiss)    then
      if (trim(lwem_form) == 'fuliou') then
        call cloud_lwpar                                           &
                    (size_drop, size_ice, size_rain,              &
                     conc_drop, conc_ice, conc_rain, conc_snow,    &
                     deltaz,                                        &
                     emmxolw,   emrndlw)             
      endif
    else
      if (trim(lwem_form) == 'ebertcurry') then
        call cloud_lwem_oneband                                 &
                   (lwpath, iwpath, size_drop, size_ice,        &
                    emmxolw, emrndlw)
      endif
    endif
 
!---------------------------------------------------------------------
!   call cloudpar to define microphysically-based sw cloud 
!   properties.
!---------------------------------------------------------------------
     if (trim(swform) == 'esfsw99')    then
       call cloudpar                                                   &
                    (size_drop, size_ice, size_rain,                  &
                     conc_drop, conc_ice, conc_rain, conc_snow,       &
                     cldext, cldsct,  cldasymm)
     endif
 
!---------------------------------------------------------------------
!   deallocate arrays
!---------------------------------------------------------------------
    if (allocated (nhi_clouds) ) then
      deallocate (nhi_clouds   )
      deallocate (nmid_clouds  )
      deallocate (nlow_clouds  )
      deallocate (press_hm     )
      deallocate (press_ml     )
      deallocate (conc         )
    endif

    deallocate (conc_drop    )
    deallocate (conc_ice     )
    deallocate (conc_rain    )
    deallocate (conc_snow    )
    deallocate (size_rain    )
    deallocate (size_snow    )


!-------------------------------------------------------------------- 
!send  necessary data to radiation_diag_mod
!----------------------------------------------------------------------
    if (Rad_control%do_diagnostics) then
      do j=JSRAD, JERAD
        call radiag_from_cloudrad (j, cldext, cldsct, cldasymm,   &
	                         deltaz, &
	          lwpath_in=1000.*lwpath, iwpath_in=1000.*iwpath,   &
                     size_drop_in=size_drop, size_ice_in=size_ice)
      end do
    endif

!----------------------------------------------------------------------
! deallocate the lwpath, iwpath, size_drop, size_ice arrays after
!  their values have been shipped to the diagnostic routine
!--------------------------------------------------------------------

      deallocate (lwpath     )
     deallocate (iwpath     )
      deallocate (size_drop  )
     deallocate (size_ice   )


end subroutine microphys_rad_driver  


!####################################################################

subroutine thickavg (nivl1    , nivl2     , nivls   ,   &
                     extivl   , ssalbivl  , asymmivl, solflxivl, &
                     solflxband, extband  , ssalbband , asymmband)
 
!----------------------------------------------------------------------c
! use the thick-averaging technique to define the single-scattering    c
! properties of the parameterization band spectral intervals from the  c
! specified spectral intervals of the particular scatterer.            c
!                                                                      c
! references:                                                          c
!                                                                      c
! edwards,j.m. and a. slingo, studies with a flexible new radiation    c
!      code I: choosing a configuration for a large-scale model.,      c
!      q.j.r. meteorological society, 122, 689-719, 1996.              c
!                                                                      c
! note: the 1.0E-100 factor to calculate asymmband is to prevent        
!       division by zero.                                              c
!----------------------------------------------------------------------c
!                                                                      c
! intent in:                                                           c
!                                                                      c
! nivl1      = interval number for the specified single-scattering     c
!              properties corresponding to the first psuedo-           c
!              monochromatic frequency in a given parameterization     c
!              band                                                    c
!                                                                      c
! nivl2      = interval number for the specified single-scattering     c
!              properties corresponding to the last psuedo-            c
!              monochromatic frequency in a given parameterization     c
!              band                                                    c
!                                                                      c
! nivls      = number of specified scattering spectral intervals       c
!                                                                      c
! extivl     = the specified spectral values of the extinction         c
!              coefficient                                             c
!                                                                      c
! ssalbivl   = the specified spectral values of the single-            c
!              scattering albedo                                       c
!                                                                      c
! asymmivl   = the specified spectral values of the asymmetry          c
!              factor                                                  c
!                                                                      c
! solflxivl  = the solar flux in each specified scattering spectral    c
!              interval                                                c
!                                                                      c
! solflxband = the solar flux in each parameterization band            c
!                                                                      c
!----------------------------------------------------------------------c
 
integer, dimension(:),    intent(in)       :: nivl1, nivl2
integer,                  intent(in)       :: nivls
real, dimension(:,:,:,:), intent(in)       :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout)    :: ssalbivl
real, dimension(:,:),     intent(in)       :: solflxivl             
real, dimension(:),       intent(in)       :: solflxband            
 
!----------------------------------------------------------------------c
!                                                                      c
! intent out:                                                          c
!                                                                      c
! extband   =  the parameterization band values of the extinction      c
!              coefficient                                             c
!                                                                      c
! ssalbband =  the parameterization band values of the single-         c
!              scattering albedo                                       c
!                                                                      c
! asymmband =  the parameterization band values of the asymmetry       c
!              factor                                                  c
!                                                                      c
!----------------------------------------------------------------------c
 
real, dimension(:,:,:,:), intent(out)  :: extband, ssalbband, asymmband
 
!----------------------------------------------------------------------c
! local variables:                                                     c
!----------------------------------------------------------------------c
 
   real, dimension(:,:,:,:), allocatable ::   refband
   real, dimension(:,:,:  ), allocatable ::   refthick, sp, sumk,   &
	   				      sumomegak, sumomegakg, &
					      sumrefthick
   integer  :: nband, ni, j, k, i
   integer  :: ilb, iub, jlb, jub, klb, kub
 
!--------------------------------------------------------------------
!   allocate local arrays
!----------------------------------------------------------------------

      ilb = lbound(ssalbivl,1)
      iub = ubound(ssalbivl,1)
      jlb = lbound(ssalbivl,2)
      jub = ubound(ssalbivl,2)
      klb = lbound(ssalbivl,3)
      kub = ubound(ssalbivl,3)

      allocate( refband     (ilb:iub, jlb:jub, klb:kub, nbands))
      allocate( refthick    (ilb:iub, jlb:jub, klb:kub) )
      allocate( sp          (ilb:iub, jlb:jub, klb:kub) )
      allocate( sumk        (ilb:iub, jlb:jub, klb:kub) )
      allocate( sumomegak   (ilb:iub, jlb:jub, klb:kub) )
      allocate( sumomegakg  (ilb:iub, jlb:jub, klb:kub) )
      allocate( sumrefthick (ilb:iub, jlb:jub, klb:kub) )
	  
      do nband = 1,nbands
  
        sumk(:,:,:) = 0.0
	sumomegak(:,:,:) = 0.0
	sumomegakg(:,:,:) = 0.0
	sumrefthick(:,:,:) = 0.0
 
        do ni = nivl1(nband),nivl2(nband)
	  do k=klb,kub
	    do j=jlb,jub    
	      do i=ilb,iub    
	        if ((ssalbivl(i,j,k,ni) +    &
		     asymmivl(i,j,k,ni)) .ne. 0.0) then
	          ssalbivl(i,j,k,ni) = MIN(ssalbivl(i,j,k,ni), 1.0)
                  sp(i,j,k) = sqrt( ( 1.0 - ssalbivl(i,j,k,ni) ) /    &
                                    ( 1.0 - ssalbivl(i,j,k,ni) *      &
                                      asymmivl(i,j,k,ni) ) )
                  refthick(i,j,k) = (1.0 - sp(i,j,k))/(1.0 + sp(i,j,k))
	          sumrefthick(i,j,k) = sumrefthick(i,j,k) +    &
				       refthick(i,j,k)*  &
				       solflxivl(nband,ni)
	          sumk(i,j,k) = sumk(i,j,k) + extivl(i,j,k,ni) *   &
                                solflxivl(nband,ni)
	          sumomegak(i,j,k) = sumomegak(i,j,k) +     &
				     ssalbivl(i,j,k,ni)*   &
				     extivl(i,j,k,ni) *   &
				     solflxivl(nband,ni)
                  sumomegakg(i,j,k) = sumomegakg(i,j,k) +    &
				      ssalbivl(i,j,k,ni)*&
                                      extivl(i,j,k,ni)*  &
				      asymmivl(i,j,k,ni) * &
                                      solflxivl(nband,ni)
	        endif
              end do
            end do
          end do
        end do
   
        do k=klb,kub
	  do j=jlb,jub    
	    do i=ilb,iub 
              extband(i,j,k,nband) = sumk(i,j,k) / solflxband(nband)
              asymmband(i,j,k,nband) = sumomegakg(i,j,k) /         &
                                       ( sumomegak(i,j,k) + 1.0E-100)
	      refband(i,j,k,nband) = sumrefthick(i,j,k)/  &
				     solflxband(nband)
              ssalbband(i,j,k,nband) = 4.0 * refband(i,j,k,nband) / &
                                       ((1.0 +    &
	            		       refband(i,j,k,nband)) ** 2 -&
                                       asymmband(i,j,k,nband) *     &
                                       (1.0 - refband(i,j,k,nband))**2 )
            end do
          end do
        end do

      end do

!-------------------------------------------------------------------
!  deallocate local arrays
!------------------------------------------------------------------
      deallocate (sumrefthick )
      deallocate (sumomegakg  )
      deallocate (sumomegak   )
      deallocate (sumk        )
      deallocate (sp          )
      deallocate (refthick    )
      deallocate (refband     )



!---------------------------------------------------------------------
  
end subroutine thickavg



!####################################################################

subroutine thinavg (nivl1    , nivl2     , nivls   ,   &
                    extivl   , ssalbivl  , asymmivl,  solflxivl, &
                    solflxband, extband  , ssalbband , asymmband)
 
!----------------------------------------------------------------------c
! use the thin-averaging technique to define the single-scattering     c
! properties of the parameterization band spectral intervals from the  c
! specified spectral intervals of the particular scatterer.            c
!                                                                      c
! references:                                                          c
!                                                                      c
! edwards,j.m. and a. slingo, studies with a flexible new radiation    c
!      code I: choosing a configuration for a large-scale model.,      c
!      q.j.r. meteorological society, 122, 689-719, 1996.              c
!                                                                      c
! note: the 1.0E-100 factor to calculate ssalbband and asymmband is to  
!       prevent division by zero.                                      c
!----------------------------------------------------------------------c
!                                                                      c
! intent in:                                                           c
!                                                                      c
! nivl1      = interval number for the specified single-scattering     c
!              properties corresponding to the first psuedo-           c
!              monochromatic frequency in a given parameterization     c
!              band                                                    c
!                                                                      c
! nivl2      = interval number for the specified single-scattering     c
!              properties corresponding to the last psuedo-            c
!              monochromatic frequency in a given parameterization     c
!              band                                                    c
!                                                                      c
! nivls      = number of specified scattering spectral intervals       c
!                                                                      c
! extivl     = the specified spectral values of the extinction         c
!              coefficient                                             c
!                                                                      c
! ssalbivl   = the specified spectral values of the single-            c
!              scattering albedo                                       c
!                                                                      c
! asymmivl   = the specified spectral values of the asymmetry          c
!              factor                                                  c
!                                                                      c
! solflxivl  = the solar flux in each specified scattering spectral    c
!              interval                                                c
!                                                                      c
! solflxband = the solar flux in each parameterization band            c
!                                                                      c
!----------------------------------------------------------------------c
 
integer, dimension(:),    intent(in)    :: nivl1, nivl2
integer,                  intent(in)    :: nivls
real, dimension(:,:,:,:), intent(in)    :: extivl, asymmivl
real, dimension(:,:,:,:), intent(inout) :: ssalbivl
real, dimension(:,:),     intent(in)    :: solflxivl                  
real, dimension(:),       intent(in)    :: solflxband                 

!----------------------------------------------------------------------c
!                                                                      c
! intent out:                                                          c
!                                                                      c
! extband   =  the parameterization band values of the extinction      c
!              coefficient                                             c
!                                                                      c
! ssalbband =  the parameterization band values of the single-         c
!              scattering albedo                                       c
!                                                                      c
! asymmband =  the parameterization band values of the asymmetry       c
!              factor                                                  c
!                                                                      c
!----------------------------------------------------------------------c
 
real, dimension(:,:,:,:), intent(out)  :: extband, ssalbband, asymmband
 
!----------------------------------------------------------------------c
! local variables:                                                     c
!----------------------------------------------------------------------c

   real, dimension(:,:,:,:), allocatable ::   refband
   real, dimension(:,:,:  ), allocatable ::   sumk,   &
	   				      sumomegak, sumomegakg
 
   integer   ::   nband, ni, i, j, k
   integer   ::   ilb, iub, jlb, jub, klb, kub

!--------------------------------------------------------------------
!   allocate local arrays
!----------------------------------------------------------------------

      ilb = lbound(ssalbivl,1)
      iub = ubound(ssalbivl,1)
      jlb = lbound(ssalbivl,2)
      jub = ubound(ssalbivl,2)
      klb = lbound(ssalbivl,3)
      kub = ubound(ssalbivl,3)

      allocate( sumk        (ilb:iub, jlb:jub, klb:kub) )
      allocate( sumomegak   (ilb:iub, jlb:jub, klb:kub) )
      allocate( sumomegakg  (ilb:iub, jlb:jub, klb:kub) )

!---------------------------------------------------------------------
      do nband = 1,nbands
  
        sumk(:,:,:) = 0.0
	sumomegak(:,:,:) = 0.0
	sumomegakg(:,:,:) = 0.0
 
        do ni = nivl1(nband),nivl2(nband)
	  do k=klb,kub
	    do j=jlb,jub    
	      do i=ilb,iub    
	        if ((ssalbivl(i,j,k,ni) +    &
		     asymmivl(i,j,k,ni)) .ne. 0.0) then
	          ssalbivl(i,j,k,ni) = MIN(ssalbivl(i,j,k,ni), 1.0)
	          sumk(i,j,k) = sumk(i,j,k) + extivl(i,j,k,ni) *   &
                                solflxivl(nband,ni)
	          sumomegak(i,j,k) = sumomegak(i,j,k) +    &
				     ssalbivl(i,j,k,ni) *  &
                                     extivl(i,j,k,ni) *   &
				     solflxivl(nband,ni)
	          sumomegakg(i,j,k) = sumomegakg(i,j,k) +    &
				      ssalbivl(i,j,k,ni) * & 
				      extivl(i,j,k,ni) *   &
				      asymmivl(i,j,k,ni) *  &
                                      solflxivl(nband,ni)
                endif
              end do
            end do
          end do
        end do
   
	do k=klb,kub
	  do j=jlb,jub    
	    do i=ilb,iub    
              extband(i,j,k,nband) = sumk(i,j,k) / solflxband(nband)
              asymmband(i,j,k,nband) = sumomegakg(i,j,k) /    &
                                       ( sumomegak(i,j,k) + 1.0E-100 )
              ssalbband(i,j,k,nband) = sumomegak(i,j,k) /   &
                                       ( sumk(i,j,k) + 1.0E-100 )
            end do
          end do
        end do

      end do

!-------------------------------------------------------------------
!  deallocate local arrays
!------------------------------------------------------------------
      deallocate (sumomegakg  )
      deallocate (sumomegak   )
      deallocate (sumk        )

!-------------------------------------------------------------------
  
end subroutine thinavg 


!####################################################################

subroutine cloud_lwpar                                    &
                    (size_drop, size_ice, size_rain,            &
                     conc_drop, conc_ice, conc_rain, conc_snow, &
                     deltaz,                                    &
                     emmxolw,   emrndlw )
 
!----------------------------------------------------------------------
! determine the infrared cloud emissivities for specified wavenumber    
! bands from parameterizations for absorption coefficients due to       
! cloud drops, cloud ice crystals, rain and snow. conceptually one      
! could have separate concentrations and sizes for "thin" or randomly   
! overlapped and for maximally overlapped clouds. for now, there is     
! one concentration and size, therefore the two emissivities are set    
! equal.                                                                
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!                                                                   
! intent in:                                                        
!                                                                   
! size_drop = the cloud drop effective diameter in microns          
!                                                                   
! size_ice  = the ice crystal effective size in microns             
!                                                                   
! size_rain = the rain drop effective diameter in microns           
!                                                                   
! conc_drop = the cloud drop liquid water concentration in grams /  
!             meter**3                                              
!                                                                   
! conc_ice = the ice water concentation in grams / meter**3         
!                                                                   
! conc_rain = the rain drop water concentration in grams / meter**3 
!                                                                   
! conc_snow = the snow concentration in grams / meter**3            
!                                                                   
! deltaz    = the thickness of pressure layers (in meters)
!----------------------------------------------------------------------
 
real, dimension (:,:,:), intent(in)     ::   size_drop, size_ice,    &
                                             size_rain,              &
					     conc_drop, conc_ice,    &
					     conc_rain, conc_snow,   &
					     deltaz
 
!----------------------------------------------------------------------
!                                                                   
! intent out:                                                       
!                                                                   
! emmxolw     = the infrared cloud emissivity for the (NLWCLDB)     
!               bands. used for maximally overlapped clouds.        
!                                                                   
! emrndlw     = the infrared cloud emissivity for the (NLWCLDB)     
!               bands. used for randomly overlapped clouds.         
!                                                                   
!----------------------------------------------------------------------
 
real, dimension (:,:,:,:), intent(out)  ::  emmxolw, emrndlw
 
!----------------------------------------------------------------------
! local variables:                                                  
!---------------------------------------------------------------------

        integer       :: k, j, i, n

!---------------------------------------------------------------------
!
! cldextbndrainlw = the specified values of the extinction          
!                 coefficient for rain water in kilometers**(-1)    
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldssalbbndrainlw = the specified values of the single-           
!                  scattering albedo for rain water                 
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldasymmbndrainlw = the specified values of the asymmetry         
!                  factor for rain water                            
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldextbndsnowlw = the specified values of the extinction          
!                 coefficient for snow water in kilometers**(-1)    
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldssalbbndsnowlw = the specified values of the single-           
!                  scattering albedo for snow water                 
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldasymmbndsnowlw = the specified values of the asymmetry         
!                  factor for snow water                            
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldextbndicelw = the specified values of the extinction           
!                 coefficient for ice particles in kilometers**(-1) 
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldssalbbndicelw = the specified values of the single-            
!                  scattering albedo for ice particles              
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldasymmbndicelw = the specified values of the asymmetry          
!                  factor for ice particles                         
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldextbnddroplw = the specified values of the extinction          
!                 coefficient for cloud drops in kilometers**(-1)   
!                 over wavenumber bands used by the radiation code  
!                                                                   
!----------------------------------------------------------------------
      real, dimension (:,:,:,:), allocatable   ::                   &
             cldextbndrainlw, cldssalbbndrainlw, cldasymmbndrainlw, &
             cldextbndsnowlw, cldssalbbndsnowlw, cldasymmbndsnowlw, &
             cldextbndicelw, cldssalbbndicelw, cldasymmbndicelw,    &
	     cldextbnddroplw,                                       &
	     abscoeff
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!     allocate allocatable arrays
!----------------------------------------------------------------------
      allocate ( cldextbndrainlw                                      &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )
      allocate ( cldssalbbndrainlw                                    &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )
      allocate ( cldasymmbndrainlw                                    &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )
      allocate ( cldextbndsnowlw                                      &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )
      allocate ( cldssalbbndsnowlw                                    &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )
      allocate ( cldasymmbndsnowlw                                    &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )
      allocate ( cldextbndicelw                                       &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )
      allocate ( cldssalbbndicelw                                     &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )
      allocate ( cldasymmbndicelw                                    &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )
      allocate ( cldextbnddroplw                                      &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )
      allocate ( abscoeff                                             &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, N_EMISS_BDS) )

!--------------------------------------------------------------------
!     compute extinction coefficient, single scattering coefficient
!     asymmetry parameter for rain
!-------------------------------------------------------------------
      call furainlw                                                  &
                (conc_rain    ,                                      &
                 cldextbndrainlw, cldssalbbndrainlw, cldasymmbndrainlw)

!----------------------------------------------------------------------
!     compute extinction coefficient, single scattering coefficient
!     asymmetry parameter for snow
!----------------------------------------------------------------------
      call fusnowlw                                                  &
                (conc_snow    ,                                      &
                 cldextbndsnowlw, cldssalbbndsnowlw, cldasymmbndsnowlw)
 
!----------------------------------------------------------------------
!     compute extinction coefficient for cloud drops
!----------------------------------------------------------------------
      call cliqlw                                                    &
                (conc_drop    ,                                      &
                 cldextbnddroplw)
 
!----------------------------------------------------------------------
!     compute extinction coefficient, single scattering coefficient
!     asymmetry parameter for cloud ice crystals
!----------------------------------------------------------------------
      call el                                                        &
                   (conc_ice    , size_ice      ,                    &
                    cldextbndicelw, cldssalbbndicelw, cldasymmbndicelw)
 
!----------------------------------------------------------------------
!       compute absorption coefficient (in km-1) for each species
!     as the product of the extinction coefficient and (1 - 
!     single scattering albedo).  the total absorption coefficient is
!     the sum of the species absorption coefficients. 
!        over a single frequency band (see goody and yung, eq. 6.72)
!     the emissivity in a layer is defined as:
!         1 - T(f)    where
!     T(f) is the flux transmissivity, which may be computed as
!         EXP(-(1.66)*(abs. coeff)*(layer thickness))
!     the factor 1.66 is the diffusivity factor (diffac).
!----------------------------------------------------------------------
      abscoeff = cldextbndicelw*(1.0E+00 - cldssalbbndicelw)   +    &
                 cldextbnddroplw                            +       &
                 cldextbndsnowlw*(1.0E+00 - cldssalbbndsnowlw) +    &
                 cldextbndrainlw*(1.0E+00 - cldssalbbndrainlw)
!----------------------------------------------------------------------
!     1.0E-3 is conversion factor from (m) to (km).
!----------------------------------------------------------------------
      do n=1,NLWCLDB
        emmxolw(:,:,:,n)  = 1.0E+00 -                               &
                 EXP(-diffac*abscoeff(:,:,:,n)*deltaz(:,:,:)*1.0E-03)
      enddo
      emrndlw  = emmxolw
 
!----------------------------------------------------------------------
!      deallocate deallocatable arrays
!----------------------------------------------------------------------
      deallocate (cldextbndrainlw)
      deallocate (cldssalbbndrainlw)
      deallocate (cldasymmbndrainlw)
      deallocate (cldextbndsnowlw)
      deallocate (cldssalbbndsnowlw)
      deallocate (cldasymmbndsnowlw)
      deallocate (cldextbndicelw)
      deallocate (cldssalbbndicelw)
      deallocate (cldasymmbndicelw)
      deallocate (cldextbnddroplw)
      deallocate (abscoeff)
 

end subroutine cloud_lwpar


!######################################################################


 subroutine cloud_lwem_oneband                                    &
                     (lwpath, iwpath, size_drop, size_ice,        &
                      emmxolw, emrndlw)

!----------------------------------------------------------------------
! determine the infrared cloud emissivities for a single broadband
! from parameterizations for absorption coefficients due to
! cloud drops and cloud ice crystals. the parameterization comes from
! from the 5-band formulation given by Ebert and Curry (1992,
! J. Geophys. Res., vol. 97, pp. 3831-3836). S. Klein  derived the
! coefficients for the 1-band version used here.
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
! intent in:
!
! size_drop = the cloud drop effective diameter in microns
!
! size_ice  = the ice crystal effective size in microns
!
! lwpath    = the liquid water path, in kg / meters**2
!
! iwpath    = the ice water path, in kg / meters**2
!
!----------------------------------------------------------------------

real, dimension (:,:,:), intent(in)     ::   size_drop, size_ice,    &
                                             lwpath, iwpath

!----------------------------------------------------------------------
!
!
! intent out:  
!
! emmxolw     = the one-band infrared cloud emissivity.
!               used for maximally overlapped clouds.
!
! emrndlw     = the one-band infrared cloud emissivity.
!               used for randomly overlapped clouds.
!
!----------------------------------------------------------------------
 
real, dimension (:,:,:,:), intent(out)  ::  emmxolw, emrndlw

!----------------------------------------------------------------------
! local variables:
!---------------------------------------------------------------------

  integer       ::  n

  real, dimension(:,:,:), allocatable :: Reff_ice,  k_liq, k_ice,  &
                                        em_lw
 
!---------------   COMPUTE LONGWAVE EMISSIVITY ----------------------!

  allocate  (  Reff_ice (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD),  &
               k_liq    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD),  &
              k_ice    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD),  &
              em_lw    (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD)   )

!   Reff_ice is the effective diameter obtained using an expression
!    provided by S. Klein

!   Reff and Deff are in microns
       Reff_ice =       &
    ((0.1033741*size_ice*size_ice + 0.2115169*(size_ice**2.272))**0.5)
  
      k_liq(:,:,:) = 140.
      k_ice(:,:,:) = 4.83591 + 1758.511/Reff_ice(:,:,:)

      
! compute combined emissivity
      em_lw(:,:,:) = 1. - exp(-1.*(k_liq(:,:,:)*lwpath(:,:,:) +   &
                                   k_ice(:,:,:)*iwpath(:,:,:) ) )

!----------------------------------------------------------------------
     do n=1,NLWCLDB
       emmxolw(:,:,:,n)  = em_lw(:,:,:)
     enddo
     emrndlw  = emmxolw

!----------------------------------------------------------------------
!      deallocate deallocatable arrays
!----------------------------------------------------------------------

     deallocate (em_lw)
     deallocate (k_ice)
     deallocate (k_liq)
     deallocate (Reff_ice)

end subroutine cloud_lwem_oneband


!######################################################################


subroutine furainlw                                         &
                (conc_rain    ,                                   &
                 cldextbndrainlw, cldssalbbndrainlw, cldasymmbndrainlw)
 
!----------------------------------------------------------------------
!      Calculates absorption coefficient for cloud rain water for
!      longwave radiation (Fu et al., 1995, J. Atmos. Sci.)
!      To be used for rain water with radii between 60 um and 1.8 mm.
!      See also notes from Q. Fu (4 Sept 98)
!      note: the size of the rain water from the microphysical
!      model (if existing) does not enter into the calculations.
!
!      Leo Donner, GFDL, 20 Mar 99
!
!----------------------------------------------------------------------
!
! intent in:                                                        
!
! conc_rain = the rain drop water concentration in grams / meter**3 
!
!----------------------------------------------------------------------

real, dimension (:,:,:), intent(in)     ::   conc_rain
 
!----------------------------------------------------------------------
! intent out:                                                       
!
! cldextbndrainlw = the specified values of the extinction          
!                 coefficient for rain water in kilometers**(-1)    
!                 over wavenumber bands used by the radiation code  
! cldssalbbndrainlw = the specified values of the single-           
!                  scattering albedo for rain water                 
!                 over wavenumber bands used by the radiation code  
! cldasymmbndrainlw = the specified values of the asymmetry         
!                  factor for rain water                            
!                 over wavenumber bands used by the radiation code  
!----------------------------------------------------------------------

real, dimension (:,:,:,:), intent(out)  :: cldextbndrainlw,    &
                                           cldssalbbndrainlw,   &
					   cldasymmbndrainlw

!----------------------------------------------------------------------
! local variables:                                                  
!---------------------------------------------------------------------
      integer       :: k, j, i, n
      integer       :: ni
      real          :: rwc0

      real, dimension (NBFL)                   :: brn, wrnf, grn
      real, dimension (:,:,:,:), allocatable   :: cldextivlrain,    &
						  cldssalbivlrain,    &
						  cldasymmivlrain
 
!---------------------------------------------------------------------
!     Boundaries for spectral bands(cm**(-1)) (7-18):
!     2200,1900,1700,1400,1250,1100,980,800,670,540,400,280,0
!
!     brn  = empirical coefficients for extinction coefficient
!             parameterization (km**-1)
!     wrnf  = empirical coefficients for single scattering albedo
!             parameterization
!     grn   = empirical coefficients for asymmetry parameter
!             parameterization
!     rwc0  = rain water content (g/m**3) used to obtain above
!             empirical coefficients. 
!     note: since model wavenumber bands are in order of increasing
!     wavenumber, the Fu coefficients have been reversed (thus are
!     in order  bands 18-7)
!
! cldextivlrain   = the specified spectral values of the extinction 
!                 coefficient for rain water in kilometers**(-1)    
!                                                                   
! cldssalbivlrain   = the specified spectral values of the single-  
!                  scattering albedo for rain water                 
!                                                                   
! cldasymmivlrain   = the specified spectral values of the asymmetry   
!                  factor for rain water                            
!---------------------------------------------------------------------
 
      data brn /                                                     &
             1.6765,  1.6149,  1.5993,  1.5862,  1.5741,  1.5647,    &
             1.5642,  1.5600,  1.5559,  1.5512,  1.5478,  1.5454/
      data wrnf /                                                    &
              .55218,  .55334,  .55488,  .55169,  .53859,  .51904,   &
              .52321,  .52716,  .52969,  .53192,  .52884,  .53233/
      data grn /                                                     &
              .90729,  .92990,  .93266,  .94218,  .96374,  .98584,   &
              .98156,  .97745,  .97467,  .97216,  .97663,  .97226/
      data rwc0 /                                                    &
                 0.5/
                                                                    
!---------------------------------------------------------------------
!     allocate allocatable arrays
!---------------------------------------------------------------------
      allocate ( cldextivlrain                                         &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NBFL) )
      allocate ( cldssalbivlrain                                       &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NBFL) )
      allocate ( cldasymmivlrain                                       &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NBFL) )
 
!-----------------------------------------------------------------------
!
!     Calculate extinction coefficient (km**(-1)) over wavenumber
!     bands of the Fu-Liou parameterization (not the radiation
!     code wavenumber bands)
!     as of 4/8/99, the asymmetry parameter is not used in the
!     infrared code. therefore the calculations are commented out.
!
!-----------------------------------------------------------------------
      do n=1,NBFL
	cldextivlrain(:,:,:,n) = brn(n)*conc_rain(:,:,:)/rwc0
	cldssalbivlrain(:,:,:,n) = wrnf(n)
!	cldasymmivlrain(:,:,:,n) = grn(n)
      enddo
 
!-----------------------------------------------------------------------
!    use band weighting factors (computed in initialization routine)
!    to derive band quantities for these quantities
!-----------------------------------------------------------------------
      cldextbndrainlw = 0.0E+00
      cldssalbbndrainlw = 0.0E+00
!     cldasymmbndrainlw = 0.0E+00
      do n = 1,NLWCLDB
	do ni = 1,NBFL
	  cldextbndrainlw(:,:,:,n) = cldextbndrainlw(:,:,:,n) +      &
               cldextivlrain(:,:,:,ni)*fulwwts(n,ni)
 	  cldssalbbndrainlw(:,:,:,n) = cldssalbbndrainlw(:,:,:,n) +  &
               cldssalbivlrain(:,:,:,ni)*fulwwts(n,ni)
!	  cldasymmbndrainlw(:,:,:,n) = cldasymmbndrainlw(:,:,:,n) +  &
!              cldasymmivlrain(:,:,:,ni)*fulwwts(n,ni)
	enddo
      enddo
!---------------------------------------------------------------------
 
      deallocate (cldextivlrain)
      deallocate (cldssalbivlrain)
      deallocate (cldasymmivlrain)

!-------------------------------------------------------------------
 
end subroutine furainlw


!#####################################################################

subroutine fusnowlw                                         &
                (conc_snow    ,                                   &
                 cldextbndsnowlw, cldssalbbndsnowlw, cldasymmbndsnowlw)
 
!-----------------------------------------------------------------------
!
!      Calculates absorption coefficient for cloud snow water for
!      longwave radiation (Fu et al., 1995, J. Atmos. Sci.)
!      To be used for snow water with radii between 60 um and 5.0 mm.
!      See also notes from Q. Fu (4 Sept 98)
!      note: the size of the snow water from the microphysical
!      model (if existing) does not enter into the calculations.
!
!      Leo Donner, GFDL, 20 Mar 99
!
!-----------------------------------------------------------------------
                                                                    
!--------------------------------------------------------------------
! intent in:                                                        
!
! conc_snow = the snow drop water concentration in grams / meter**3 
!---------------------------------------------------------------------

real, dimension (:,:,:), intent(in)     ::   conc_snow
 
!----------------------------------------------------------------------
!                                                                   
! intent out:                                                       
!                                                                   
! cldextbndsnowlw = the specified values of the extinction          
!                 coefficient for snow water in kilometers**(-1)    
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldssalbbndsnowlw = the specified values of the single-           
!                  scattering albedo for snow water                 
!                 over wavenumber bands used by the radiation code  
!                                                                   
! cldasymmbndsnowlw = the specified values of the asymmetry         
!                  factor for snow water                            
!                 over wavenumber bands used by the radiation code  
!----------------------------------------------------------------------
 
real, dimension (:,:,:,:), intent(out)  :: cldextbndsnowlw,    &
                             cldssalbbndsnowlw, cldasymmbndsnowlw
 
!----------------------------------------------------------------------
! local variables:                                                  
!---------------------------------------------------------------------
        integer                :: k, j, i,  n
        integer                :: ni
        real                   :: swc0

        real, dimension (NBFL)                   :: brn, wrnf, grn
        real, dimension (:,:,:,:), allocatable   :: cldextivlsnow,   &
						    cldssalbivlsnow,   &
						    cldasymmivlsnow

!---------------------------------------------------------------------
!     Boundaries for spectral bands(cm**(-1)) (7-18):
!     2200,1900,1700,1400,1250,1100,980,800,670,540,400,280,0
!
!
!     brn  = empirical coefficients for extinction coefficient
!             parameterization (km**-1)
!     wrnf  = empirical coefficients for single scattering albedo
!             parameterization
!     grn   = empirical coefficients for asymmetry parameter
!             parameterization
!     swc0  = snow water content (g/m**3) used to obtain above
!             empirical coefficients. 
!     note: since model wavenumber bands are in order of increasing
!     wavenumber, the Fu coefficients have been reversed (thus are
!     in order  bands 18-7)
!
! cldextivlsnow   = the specified spectral values of the extinction 
!                 coefficient for snow water in kilometers**(-1)    
!                                                                   
! cldssalbivlsnow   = the specified spectral values of the single-  
!                  scattering albedo for snow water                 
!                                                                   
! cldasymmivlsnow   = the specified spectral values of the asymmetry   
!                  factor for snow water                            
!---------------------------------------------------------------------

      data brn /                                                     &
              .87477,  .85421,  .84825,  .84418,  .84286,  .84143,   &
              .84097,  .84058,  .84029,  .83995,  .83979,  .83967/
      data wrnf /                                                    &
              .55474,  .53160,  .54307,  .55258,  .54914,  .52342,   &
              .52446,  .52959,  .53180,  .53182,  .53017,  .53296/
      data grn /                                                     &
              .93183,  .97097,  .95539,  .94213,  .94673,  .98396,   &
              .98274,  .97626,  .97327,  .97330,  .97559,  .97173/
      data swc0 /                                                    &
                 0.5/
                                                                    
!---------------------------------------------------------------------
!     allocate allocatable arrays
!---------------------------------------------------------------------
       allocate ( cldextivlsnow                                   & 
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NBFL) )
       allocate ( cldssalbivlsnow                                 &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NBFL) )
       allocate ( cldasymmivlsnow                                 &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NBFL) )
 
!-----------------------------------------------------------------------
!     Calculate extinction coefficient (km**(-1)) over wavenumber
!     bands of the Fu-Liou parameterization (not the radiation
!     code wavenumber bands)
!     as of 4/8/99, the asymmetry parameter is not used in the
!     infrared code. therefore the calculations are commented out.
!-----------------------------------------------------------------------
      do n=1,NBFL
	cldextivlsnow(:,:,:,n) = brn(n)*conc_snow(:,:,:)/swc0
	cldssalbivlsnow(:,:,:,n) = wrnf(n)
!	cldasymmivlsnow(:,:,:,n) = grn(n)
      enddo
 
!-----------------------------------------------------------------------
!    use band weighting factors (computed in initialization routine)
!    to derive band quantities for these quantities
!-----------------------------------------------------------------------
      cldextbndsnowlw = 0.0E+00
      cldssalbbndsnowlw = 0.0E+00
!     cldasymmbndsnowlw = 0.0E+00
      do n = 1,NLWCLDB
	do ni = 1,NBFL
	  cldextbndsnowlw(:,:,:,n) = cldextbndsnowlw(:,:,:,n) +      &
               cldextivlsnow(:,:,:,ni)*fulwwts(n,ni)
 	  cldssalbbndsnowlw(:,:,:,n) = cldssalbbndsnowlw(:,:,:,n) +  &
               cldssalbivlsnow(:,:,:,ni)*fulwwts(n,ni)
!	  cldasymmbndsnowlw(:,:,:,n) = cldasymmbndsnowlw(:,:,:,n) +  &
!              cldasymmivlsnow(:,:,:,ni)*fulwwts(n,ni)
	enddo
      enddo
!----------------------------------------------------------------------
 
      deallocate (cldextivlsnow)
      deallocate (cldssalbivlsnow)
      deallocate (cldasymmivlsnow)

!---------------------------------------------------------------------
 
end subroutine fusnowlw



!####################################################################

subroutine cliqlw (conc_drop, cldextbnddroplw)
 
!-----------------------------------------------------------------------
!
!     Calculates longwave absorption optical depth for liquid.
!     Follows Held et al. (J. Atmos. Sci., 1993).
!     
!     Leo Donner, GFDL, 1 Feb 1999
!
!-----------------------------------------------------------------------
!
! intent in:                                                       
!                                                                  
!     conc_drop = the cloud drop concentration in grams / meter**3  
!                                                                  
!-----------------------------------------------------------------------
 
real, dimension (:,:,:), intent(in)     ::   conc_drop
                                                                   
!----------------------------------------------------------------------
! intent out:                                                      
!                                                                  
! cldextbnddroplw = the specified values of the extinction         
!                 coefficient for cloud drops in kilometers**(-1)  
!                 over wavenumber bands used by the radiation code 
!                                                                  
!-----------------------------------------------------------------------

real, dimension (:,:,:,:), intent(out)  ::   cldextbnddroplw
 
!---------------------------------------------------------------------
! local variables:                                                   
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!     alpha     = frequency-independent parameter (in m**2/g) for
!                 absorption due to cloud drops in the infrared.
!                 this value is given in held et al, JAS, 1993.
!--------------------------------------------------------------------

      integer     ::   n
      real, save  ::   alpha=0.1           
!---------------------------------------------------------------------

      do n=1,NLWCLDB
	cldextbnddroplw(:,:,:,n) = 1.0E+03*alpha*conc_drop(:,:,:)
      enddo
 
end subroutine cliqlw

!####################################################################

subroutine el  (conc_ice    , size_ice      ,                   &
                cldextbndicelw, cldssalbbndicelw, cldasymmbndicelw)
 
!-----------------------------------------------------------------------
!
!     calculates total optical depth and scattering optical depth
!     for infrared radiation using Fu and Liou (1993,
!     JAS). To be used for crystal effective sizes from 20 to 130 um.
!     limits changed to 18.6 to 130.2 um on 2 April 1999 to
!     match shortwave limits.
!
!
!        Leo Donner, GFDL, 29 Aug 98
!
!-----------------------------------------------------------------------
!
! intent in:                                                       
!                                                                  
! conc_ice = the ice crystal concentation in grams / meter**3      
!                                                                  
! size_ice = the ice crystal effective size in microns             
!---------------------------------------------------------------------
 
real, dimension (:,:,:), intent(in)  ::  conc_ice, size_ice
                                                                   
!----------------------------------------------------------------------
!                                                                  
! intent out:                                                      
!                                                                  
! cldextbndicelw = the specified values of the extinction          
!                 coefficient for ice particles in kilometers**(-1)
!                 over wavenumber bands used by the radiation code 
!                                                                  
! cldssalbbndicelw = the specified values of the single-           
!                  scattering albedo for ice particles             
!                 over wavenumber bands used by the radiation code 
!                                                                  
! cldasymmbndicelw = the specified values of the asymmetry         
!                  factor for ice particles                        
!                 over wavenumber bands used by the radiation code 
!--------------------------------------------------------------------- 
 
real, dimension (:,:,:,:), intent(out)   ::  cldextbndicelw,   &
                                             cldssalbbndicelw,    &
					     cldasymmbndicelw
                                                                   
!---------------------------------------------------------------------
! local variables:                                                   
!--------------------------------------------------------------------
 
!----------------------------------------------------------------------
! cldextivlice = the specified spectral values of the extinction   
!                 coefficient for ice particles in kilometers**(-1)
!                                                                  
! cldssalbivlice = the specified spectral values of the single-    
!                  scattering albedo for ice particles             
!                                                                  
! cldasymmivlice = the specified spectral values of the asymmetry  
!                  factor for ice particles                        
!----------------------------------------------------------------------
        real, dimension (:,:,:,:), allocatable   ::   cldextivlice, &
						      cldssalbivlice, &
						      cldasymmivlice

!----------------------------------------------------------------------
!     a0,a1,a2 = empirical coefficients for extinction coefficient
!             parameterization
!     b     = empirical coefficients for single scattering albedo
!             parameterization
!     cpr   = empirical coefficients for asymmetry parameter
!             parameterization
!     note: since model wavenumber bands are in order of increasing
!     wavenumber, the Fu coefficients have been reversed (thus are
!     in order  bands 18-7)
!----------------------------------------------------------------------
        real, dimension (1:NBFL,0:NBB)     ::   b
        real, dimension (1:NBFL,0:NBC)     ::   cpr
        real, dimension (NBFL)             ::   a0, a1, a2
 
        data a0 /                                                    &
          -7.752E-03,  -1.741E-02,  -1.704E-02,  -1.151E-02,         &
          -1.026E-02,  -8.294E-03,  -1.153E-02,  -9.609E-03,         &
          -9.061E-03,  -8.441E-03,  -8.088E-03,  -7.770E-03/
        data a1 /                                                    &
           4.624,   5.541,   4.830,   4.182,   4.105,   3.925,       &
           4.109,   3.768,   3.741,   3.715,   3.717,   3.734/
        data a2 /                                                    &
         -42.010, -58.420,  16.270,  31.130,  16.360,   1.315,       &
          17.320,  34.110,  26.480,  19.480,  17.170,  11.850/
        data (b(n,0),n=1,NBFL) /                                     &
          0.8079,   0.3964,   0.1028,   0.3254,   0.5207,   0.5631,  &
          0.2307,   0.2037,   0.3105,   0.3908,   0.3014,   0.1996/
        data (b(n,1),n=1,NBFL) /                                     &
         -0.7004E-02, -0.3155E-02,  0.5019E-02,  0.3434E-02,         &
         -0.9778E-03, -0.1434E-02,  0.3830E-02,  0.4247E-02,         &
          0.2603E-02,  0.1272E-02,  0.2639E-02,  0.3780E-02/
        data (b(n,2),n=1,NBFL) /                                     &
          0.5209E-04,  0.6417E-04, -0.2024E-04, -0.3081E-04,         &
          0.3725E-05,  0.6298E-05, -0.1616E-04, -0.1810E-04,         &
          -0.1139E-04, -0.5564E-05, -0.1116E-04, -0.1491E-04/
        data (b(n,3),n=1,NBFL) /                                     &
         -0.1425E-06, -0.2979E-06,  0.0000E+00,  0.9143E-07,         &
          0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00,         &
          0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00/
        data (cpr(n,0),n=1,NBFL) /                                   &
          0.2292,   0.7024,   0.7290,   0.7678,   0.8454,   0.9092,  &
          0.9167,   0.8815,   0.8765,   0.8915,   0.8601,   0.7955/
        data (cpr(n,1),n=1,NBFL) /                                   &
           1.724E-02,   4.581E-03,   2.132E-03,   2.571E-03,         &
           1.429E-03,   9.295E-04,   5.499E-04,   9.858E-04,         &
           1.198E-03,   1.060E-03,   1.599E-03,   2.524E-03/
        data (cpr(n,2),n=1,NBFL) /                                   &
          -1.573E-04,  -3.054E-05,  -5.584E-06,  -1.041E-05,         &
          -5.859E-06,  -3.877E-06,  -1.507E-06,  -3.116E-06,         &
          -4.485E-06,  -4.171E-06,  -6.465E-06,  -1.022E-05/
        data (cpr(n,3),n=1,NBFL) /                                   &
           4.995E-07,   6.684E-08,   0.000E+00,   0.000E+00,         &
           0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00,         &
           0.000E+00,   0.000E+00,   0.000E+00,   0.000E+00/
!---------------------------------------------------------------------

        integer     :: k, j, i, n, ni


!---------------------------------------------------------------------
!     allocate allocatable arrays
!---------------------------------------------------------------------
      allocate ( cldextivlice                                      &
                 (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NBFL) )
      allocate ( cldssalbivlice                                    &
                 (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NBFL) )
      allocate ( cldasymmivlice                                    &
                  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NBFL) )

!-----------------------------------------------------------------------
!     Calculate extinction coefficient (km**(-1)) over wavenumber
!     bands of the Fu-Liou parameterization (not the radiation
!     code wavenumber bands)
!-----------------------------------------------------------------------
      do n=1,NBFL                                                      
        cldextivlice(:,:,:,n) = 1.0E+03*conc_ice(:,:,:)*              &
            (a0(n) + a1(n)/size_ice(:,:,:) + a2(n)/size_ice(:,:,:)**2)
      enddo

!-----------------------------------------------------------------------
!     Calculate single-scattering albedo and asymmetry parameter.
!     these are dimensionless. as of 4/8/99, the asymmetry parameter is
!     not used in the infrared code. therefore the calculations are
!     commented out.
!-----------------------------------------------------------------------
      do n=1,NBFL                                                      
 	cldssalbivlice(:,:,:,n) = 1.0E+00 -                           &
        (b(n,0) + b(n,1)*size_ice(:,:,:) + b(n,2)*size_ice(:,:,:)**2 + &
          b(n,3)*size_ice(:,:,:)**3)
      enddo
 
!     do n=1,NBFL                                                      
! 	cldasymmivlice(:,:,:,n) = cpr(n,0) +                          &
!        cpr(n,1)*size_ice(:,:,:) + cpr(n,2)*size_ice(:,:,:)**2 +     &
!        cpr(n,3) * size_ice(:,:,:)**3
!     enddo
 
!-----------------------------------------------------------------------
!    use band weighting factors (computed in initialization routine)
!    to derive band quantities for these quantities
!-----------------------------------------------------------------------
      cldextbndicelw = 0.0E+00
!     cldasymmbndicelw = 0.0E+00
      cldssalbbndicelw = 0.0E+00
      do n = 1,NLWCLDB                                                 
	do ni = 1,NBFL                                                 
	  cldextbndicelw(:,:,:,n) = cldextbndicelw(:,:,:,n) +         &
               cldextivlice(:,:,:,ni)*fulwwts(n,ni)
 	  cldssalbbndicelw(:,:,:,n) = cldssalbbndicelw(:,:,:,n) +     &
               cldssalbivlice(:,:,:,ni)*fulwwts(n,ni)
!	  cldasymmbndicelw(:,:,:,n) = cldasymmbndicelw(:,:,:,n) +     &
!              cldasymmivlice(:,:,:,ni)*fulwwts(n,ni)
	enddo
      enddo
!--------------------------------------------------------------------
 

end subroutine el


!######################################################################

subroutine cloudpar                                            &
                    (size_drop, size_ice, size_rain,                 &
                     conc_drop, conc_ice, conc_rain, conc_snow,      &
                     cldext, cldsct, cldasymm)
 
!----------------------------------------------------------------------
! determine the parameterization band values of the single scattering   
! parameters (extinction coefficient, scattering coefficient and   
! asymmetry factor) for clouds from the size and/or concentration of    
! each constituent (cloud drops, rain drops, ice crystals and snow)     
! present.                                                              
!----------------------------------------------------------------------
!                                                                  
! intent in:                                                       
!                                                                  
! size_drop = the cloud drop effective diameter in microns         
!                                                                  
! size_ice  = the ice crystal effective size in microns            
!                                                                  
! size_rain = the rain drop effective diameter in microns          
!                                                                  
! conc_drop = the cloud drop liquid water concentration in grams / 
!             meter**3                                             
!                                                                  
! conc_ice = the ice water concentation in grams / meter**3        
!                                                                  
! conc_rain = the rain drop water concentration in grams / meter**3
!                                                                  
! conc_snow = the snow concentration in grams / meter**3           
!----------------------------------------------------------------------

real, dimension (:,:,:), intent(in)        ::   size_drop, size_ice,  &
                                                size_rain,             &
					        conc_drop, conc_ice,   &
					        conc_rain, conc_snow
 
!----------------------------------------------------------------------
! intent inout:                                                      
!                                                                  
! cldext      = the parameterization band values of the cloud      
!               extinction coefficient in kilometer**(-1)          
!                                                                  
! cldsct      = the parameterization band values of the cloud      
!               scattering coefficient in kilometer**(-1)          
!                                                                  
! cldasymm    = the parameterization band values of the asymmetry  
!               factor                                             
!----------------------------------------------------------------------
 real, dimension (:,:,:,:), intent(inout)   ::  cldasymm, cldext,  &
                                                cldsct
 
!---------------------------------------------------------------------
! local variables:                                                   
!--------------------------------------------------------------------
       real, dimension (:,:,:,:), allocatable ::  cldextivlliq,     &
                                  cldssalbivlliq, cldasymmivlliq,   &
				  cldextivlice,            &
                                  cldssalbivlice, cldasymmivlice,     &
				  cldextivlrain,           &
                                  cldssalbivlrain, cldasymmivlrain,   &
				  cldextivlsnow,           &
                                  cldssalbivlsnow, cldasymmivlsnow
       real, dimension (:,:,:,:), allocatable ::  cldextbandliq,     &
                                  cldssalbbandliq, cldasymmbandliq,   &
				  cldextbandice,           &
                                  cldssalbbandice, cldasymmbandice,   &
				  cldextbandrain,          &
                                  cldssalbbandrain, cldasymmbandrain, &
				  cldextbandsnow,          &
                                  cldssalbbandsnow, cldasymmbandsnow
       real, dimension (:,:,:), allocatable   ::  cldscatice,     &
				  cldscatliq,  &
                                  cldscatrain, cldscatsnow
 
       integer       :: k, j, i, n, nband

!----------------------------------------------------------------------
! define the single scattering parameters for cloud drops.         
!----------------------------------------------------------------------
      allocate (cldextbandliq  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                nbands) )
      allocate (cldssalbbandliq(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                nbands) )
      allocate (cldasymmbandliq(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                nbands) )
      allocate (cldextivlliq   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                NLIQCLDIVLS) )
      allocate (cldssalbivlliq (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                NLIQCLDIVLS) )
      allocate (cldasymmivlliq (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                              NLIQCLDIVLS) )
      call slingo                                                     &
                    (conc_drop   , size_drop,                         &
                     cldextivlliq, cldssalbivlliq, cldasymmivlliq)      
 
!----------------------------------------------------------------------
! use the thick-averaging technique to define the single-scattering
! properties of the parameterization band spectral intervals from the   
! specified spectral intervals for cloud drops.
!----------------------------------------------------------------------
      call thickavg                                                  &
                  (nivl1liqcld , nivl2liqcld    , NLIQCLDIVLS   ,     &
                   cldextivlliq, cldssalbivlliq , cldasymmivlliq,     &
                   solivlliqcld, solflxband,       &
                   cldextbandliq, cldssalbbandliq, cldasymmbandliq)
      deallocate (      cldextivlliq    )
      deallocate (      cldssalbivlliq  )
      deallocate (      cldasymmivlliq  )

!----------------------------------------------------------------------
! define the single scattering parameters for rain drops.          
!----------------------------------------------------------------------
      allocate (cldextbandrain  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 nbands) )
      allocate (cldssalbbandrain(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 nbands) )
      allocate (cldasymmbandrain(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 nbands) )
      allocate (cldextivlrain   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 NRAINCLDIVLS) )
      allocate (cldssalbivlrain (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 NRAINCLDIVLS) )
      allocate (cldasymmivlrain (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 NRAINCLDIVLS) )
      call savijarvi                                                  &
                    (conc_rain    , size_rain      ,                  &
                     cldextivlrain, cldssalbivlrain, cldasymmivlrain) 
 
!----------------------------------------------------------------------
! use the thick-averaging technique to define the single-scattering
! properties of the parameterization band spectral intervals from the   
! specified spectral intervals for rain drops.
!----------------------------------------------------------------------
      call thickavg                                                  &
                  (nivl1raincld , nivl2raincld    , NRAINCLDIVLS   ,  &
                   cldextivlrain, cldssalbivlrain , cldasymmivlrain,  &
                   solivlraincld, solflxband,         &
                  cldextbandrain, cldssalbbandrain, cldasymmbandrain)
      deallocate (      cldextivlrain   )
      deallocate (      cldssalbivlrain )
      deallocate (      cldasymmivlrain )

!----------------------------------------------------------------------
! define the single scattering parameters for ice crystals.        
!----------------------------------------------------------------------
      allocate (cldextbandice  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                nbands) )
      allocate (cldssalbbandice(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                nbands) )
      allocate (cldasymmbandice(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                nbands) )
    if (using_fu) then
      allocate (cldextivlice   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                              NICECLDIVLS) )
      allocate (cldssalbivlice (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                NICECLDIVLS) )
      allocate (cldasymmivlice (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                NICECLDIVLS) )
      call fu                                                         &
                    (conc_ice    , size_ice      ,                    &
                     cldextivlice, cldssalbivlice, cldasymmivlice)      
 
!----------------------------------------------------------------------
! use the thick-averaging technique to define the single-scattering
! properties of the parameterization band spectral intervals from the   
! specified spectral intervals for ice crystals.                                      
!----------------------------------------------------------------------
      call thickavg                                                  &
                  (nivl1icecld , nivl2icecld    , NICECLDIVLS   ,     &
                   cldextivlice, cldssalbivlice , cldasymmivlice,     &
                   solivlicecld, solflxband,         &
                  cldextbandice, cldssalbbandice, cldasymmbandice)
    else
      allocate (cldextivlice   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                              NICESOLARCLDIVLS) )
      allocate (cldssalbivlice (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                NICESOLARCLDIVLS) )
      allocate (cldasymmivlice (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, &
                                NICESOLARCLDIVLS) )
      call icesolar                                                   &
                    (conc_ice    , size_ice      ,                    &
                     cldextivlice, cldssalbivlice, cldasymmivlice)      
 
!----------------------------------------------------------------------
! use the thick-averaging technique to define the single-scattering
! properties of the parameterization band spectral intervals from the   
! specified spectral intervals for ice crystals.                                      
!----------------------------------------------------------------------
      call thickavg                                                  &
                  (nivl1icesolcld, nivl2icesolcld, NICESOLARCLDIVLS, &
                   cldextivlice, cldssalbivlice , cldasymmivlice,     &
                   solivlicesolcld, solflxband,         &
                  cldextbandice, cldssalbbandice, cldasymmbandice)
    endif

      deallocate (      cldextivlice    )
      deallocate (      cldssalbivlice  )
      deallocate (      cldasymmivlice  )

!----------------------------------------------------------------------
! define the single scattering parameters for snow.                
!----------------------------------------------------------------------
      allocate (cldextbandsnow  (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 nbands) )
      allocate (cldssalbbandsnow(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 nbands) )
      allocate (cldasymmbandsnow(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 nbands) )
      allocate (cldextivlsnow   (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 NSNOWCLDIVLS) )
      allocate (cldssalbivlsnow (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 NSNOWCLDIVLS) )
      allocate (cldasymmivlsnow (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD,&
                                 NSNOWCLDIVLS) )
      call snowsw                                                     &
                    (conc_snow,                                       &
                     cldextivlsnow, cldssalbivlsnow, cldasymmivlsnow)
 
!----------------------------------------------------------------------
! use the thick-averaging technique to define the single-scattering
! properties of the parameterization band spectral intervals from the   
! specified spectral intervals for snow.                    
!----------------------------------------------------------------------
      call thickavg                                                  &
                  (nivl1snowcld , nivl2snowcld    , NSNOWCLDIVLS   ,  &
                   cldextivlsnow, cldssalbivlsnow , cldasymmivlsnow,  &
                   solivlsnowcld, solflxband,        &
                  cldextbandsnow, cldssalbbandsnow, cldasymmbandsnow)
 
      deallocate (      cldextivlsnow   )
      deallocate (      cldssalbivlsnow )
      deallocate (      cldasymmivlsnow )

!----------------------------------------------------------------------
! combine the single-scattering properties for all the constituents     
! to define the corresponding values for the clouds.                    
!----------------------------------------------------------------------
      allocate (   cldscatice      (ISRAD:IERAD, JSRAD:JERAD,    &
				    KSRAD:KERAD) )
      allocate (   cldscatliq      (ISRAD:IERAD, JSRAD:JERAD,    &
				    KSRAD:KERAD) )
      allocate (   cldscatrain     (ISRAD:IERAD, JSRAD:JERAD, &
				    KSRAD:KERAD) )
      allocate (   cldscatsnow     (ISRAD:IERAD, JSRAD:JERAD,   &
				    KSRAD:KERAD) )
      do nband = 1,nbands
	do k=KSRAD,KERAD
	  do j=JSRAD,JERAD
	    do i=ISRAD,IERAD
	      cldscatliq(i,j,k) = cldssalbbandliq(i,j,k,nband) *   &
                                  cldextbandliq(i,j,k,nband) 
              cldscatrain(i,j,k) = cldssalbbandrain(i,j,k,nband) *   &
                                   cldextbandrain(i,j,k,nband) 
              cldscatice(i,j,k) = cldssalbbandice(i,j,k,nband) *      &
                                  cldextbandice(i,j,k,nband) 
              cldscatsnow(i,j,k) = cldssalbbandsnow(i,j,k,nband) *    &
                                   cldextbandsnow(i,j,k,nband) 
              cldext(i,j,k,nband) =    &
	      			   cldextbandliq(i,j,k,nband) + &
                                   cldextbandrain(i,j,k,nband) +  &
                                   cldextbandice(i,j,k,nband) +   &
                                   cldextbandsnow(i,j,k,nband)
	      cldsct(i,j,k,nband) =   &
				   cldscatliq(i,j,k) + &
                                   cldscatrain(i,j,k) +           &
                                   cldscatice(i,j,k) +          &
                                   cldscatsnow(i,j,k) 
	      cldasymm(i,j,k,nband) =    &
			           (cldasymmbandliq(i,j,k,nband) *    &
                                    cldscatliq(i,j,k) +           &
                                    cldasymmbandrain(i,j,k,nband) *  &
                                    cldscatrain(i,j,k) +           &
                                    cldasymmbandice(i,j,k,nband) *  &
                                    cldscatice(i,j,k) +             &
                                    cldasymmbandsnow(i,j,k,nband) *  &
                                    cldscatsnow(i,j,k) ) /          &
                                    ( cldsct(i,j,k,nband) + 1.0E-100 )
            end do
          end do
        end do
      end do
!---------------------------------------------------------------------
 
      deallocate (      cldscatice      )
      deallocate (      cldscatliq      )
      deallocate (      cldscatrain     )
      deallocate (      cldscatsnow     )

      deallocate (      cldextbandliq   )
      deallocate (      cldssalbbandliq )
      deallocate (      cldasymmbandliq )
      deallocate (      cldextbandice   )
      deallocate (      cldssalbbandice )
      deallocate (      cldasymmbandice )
      deallocate (      cldextbandrain  )
      deallocate (      cldssalbbandrain)
      deallocate (      cldasymmbandrain)
      deallocate (      cldextbandsnow  )
      deallocate (      cldssalbbandsnow)
      deallocate (      cldasymmbandsnow)
!--------------------------------------------------------------------


end subroutine cloudpar



!######################################################################

subroutine slingo                                               &
                    (conc_drop   , size_drop     ,                    &
                     cldextivlliq, cldssalbivlliq, cldasymmivlliq)      
 
!----------------------------------------------------------------------
! define the single scattering parameters for cloud drops using the     
! Slingo parameterization for his spectral intervals.                   
!                                                                       
! references:                                                           
!                                                                       
! slingo, a., a gcm parameterization of the shortwave properties of     
!      water clouds., j. atmos. sci.,46, 1419-1427, 1989.               
!                                                                       
! notes: size_drop is inputted as a diameter value; it is divided in    
!        half to get an effective radius, the quantity that is used.    
!                                                                       
!        the cloud drop effective radius can only be 4.2 <= re <=       
!        16.6 microns.                                                  
!                                                                       
!        the single scattering properties for wavenumbers < 2500 cm-1   
!        are assigned the values in the first interval, since the       
!        formulation is not defined for those wavenumbers.              
!                                                                       
!        the extinction coefficient is converted to kilometer**(-1)     
!        the unit utilized by the shortwave routine Swresf.             
!                                                                       
!        a value of 1.0E-100 is added to the size so that no division   
!        by zero occurs when the size is zero, in defining the          
!        extinction coefficient.                                        
!----------------------------------------------------------------------
!                                                                       
! intent in:                                                            
!                                                                       
! conc_drop = the cloud drop liquid water concentration in grams /      
!             meter**3                                                  
!                                                                       
! size_drop = the cloud drop effective diameter in microns              
!----------------------------------------------------------------------
 
real, dimension (:,:,:), intent(in)     ::   conc_drop, size_drop
 
!----------------------------------------------------------------------
! intent out:                                                           
!                                                                       
! cldextivlliq   = the specified spectral values of the extinction      
!                  coefficient in kilometer**(-1) for drops             
!                                                                       
! cldssalbivlliq = the specified spectral values of the single-         
!                  scattering albedo for drops                          
!                                                                       
! cldasymmivlliq = the specified spectral values of the asymmetry       
!                  factor for drops                                     
!----------------------------------------------------------------------
 
real, dimension (:,:,:,:), intent(out)  ::  cldextivlliq,             &
                                       cldssalbivlliq, cldasymmivlliq
 
!---------------------------------------------------------------------
! local variables:                                                   
!--------------------------------------------------------------------

      integer                                 :: k, j, i, n, ni
      real, dimension (NLIQCLDIVLS)           ::  a, b, c, d, e, f
      real, dimension (:,:,:), allocatable    ::  size
 
      data a /-1.023E+00, 1.950E+00, 1.579E+00, 1.850E+00, 1.970E+00, &
               2.237E+00, 2.463E+00, 2.551E+00, 2.589E+00, 2.632E+00, &
               2.497E+00, 2.622E+00, 2.650E+00, 3.115E+00, 2.895E+00, &
               2.831E+00, 2.838E+00, 2.672E+00, 2.698E+00, 2.668E+00, &
               2.801E+00, 3.308E+00, 2.944E+00, 3.094E+00 /
      data b / 1.933E+00, 1.540E+00, 1.611E+00, 1.556E+00, 1.501E+00, &
               1.452E+00, 1.420E+00, 1.401E+00, 1.385E+00, 1.365E+00, &
               1.376E+00, 1.362E+00, 1.349E+00, 1.244E+00, 1.315E+00, &
               1.317E+00, 1.300E+00, 1.320E+00, 1.315E+00, 1.307E+00, &
               1.293E+00, 1.246E+00, 1.270E+00, 1.252E+00 /
      data c / 2.500E-02, 4.490E-01, 1.230E-01, 1.900E-04, 1.200E-03, &
               1.200E-04, 2.400E-04, 6.200E-05,-2.800E-05,-4.600E-05, &
               9.800E-06, 3.300E-06, 2.300E-06,-2.700E-07,-1.200E-07, &
              -1.200E-06, 0.000E+00, 0.000E+00, 1.000E-06, 0.000E+00, &
               1.000E-06,-3.000E-07,-6.500E-07, 7.900E-07 /
      data d / 1.220E-02, 1.540E-03, 9.350E-03, 2.540E-03, 2.160E-03, &
               6.670E-04, 8.560E-04, 2.600E-04, 8.000E-05, 5.000E-05, &
               2.100E-05, 2.800E-06, 1.700E-06, 1.400E-06, 4.400E-07, &
               4.000E-07, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
               0.000E+00, 2.360E-07, 4.330E-07, 3.690E-07 /
      data e / 7.260E-01, 8.310E-01, 8.510E-01, 7.690E-01, 7.400E-01, &
               7.490E-01, 7.540E-01, 7.730E-01, 7.800E-01, 7.840E-01, &
               7.830E-01, 8.060E-01, 8.090E-01, 8.040E-01, 8.180E-01, &
               8.280E-01, 8.250E-01, 8.280E-01, 8.200E-01, 8.400E-01, &
               8.360E-01, 8.390E-01, 8.410E-01, 8.440E-01 /
      data f / 6.652E+00, 6.102E+00, 2.814E+00, 5.171E+00, 7.469E+00, &
               6.931E+00, 6.555E+00, 5.405E+00, 4.989E+00, 4.745E+00, &
               5.035E+00, 3.355E+00, 3.387E+00, 3.520E+00, 2.989E+00, &
               2.492E+00, 2.776E+00, 2.467E+00, 3.004E+00, 1.881E+00, &
               2.153E+00, 1.946E+00, 1.680E+00, 1.558E+00 /

!---------------------------------------------------------------------
!  allocate local array.
!---------------------------------------------------------------------
        allocate     ( size(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
 
 !-----------------------------------------------------------------
 !  compute scattering parameters for cloud drops.  convert input size 
 !  from diameter to radius.  bypass calculations if no drops present.
 !-----------------------------------------------------------------
        do k = KSRAD,KERAD
	  do j = JSRAD,JERAD
	    do i = ISRAD,IERAD
	      if (conc_drop(i,j,k) == 0.0) then
	        cldextivlliq(i,j,k,:) = 0.0
	        cldasymmivlliq(i,j,k, :) = 0.0
	        cldssalbivlliq(i,j,k, :) = 0.0
              else
                size(i,j,k) = 0.5*size_drop(i,j,k)
	        if ((size(i,j,k) >= 4.2E+00 .and.  &
		     size(i,j,k) <= 1.66E+01 )) then                  
                  do ni = 1,NLIQCLDIVLS
                    cldextivlliq(i,j,k,ni) = 1.0E+03*conc_drop(i,j,k)* &
                                             (1.0E-02*a(ni) + (b(ni)/  &
                                             size(i,j,k) + 1.0E-100 ) )
                    cldssalbivlliq(i,j,k,ni) = 1.0 - ( c(ni) + d(ni)* &
					       size(i,j,k) )
                    cldasymmivlliq(i,j,k,ni) = e(ni) + 1.0E-03*f(ni)*  &
                                               size(i,j,k)
                  end do
	        else
                  call error_mesg('slingo',  &
	       	           'cloud drop size out of range', FATAL)
	        endif			     
              endif			     
            end do
          end do
        end do

        deallocate (size)

!-------------------------------------------------------------------


end subroutine slingo




!#####################################################################

subroutine savijarvi                                          &
                    (conc_rain    , size_rain      ,                &
                     cldextivlrain, cldssalbivlrain, cldasymmivlrain)      
 
!----------------------------------------------------------------------
! define the single scattering parameters for rain drops using the      
! Savijarvi parameterization for his spectral intervals.                
!                                                                       
! references:                                                           
!                                                                       
! savijarvi, h., shortwave optical properties of rain., tellus, 49a,    
!      177-181, 1997.                                                   
!                                                                       
! notes: size_rain is inputted as a diameter value; it is divided in    
!        half to get an effective radius, the quantity that is used.    
!                                                                       
!        the rain drop effective radius can only be 16.6 < re < 5000    
!        microns.
!                                                                       
!        the single scattering properties for wavenumbers < 2500 cm-1   
!        are assigned the values in the first interval, since the       
!        formulation is not defined for those wavenumbers.              
!                                                                       
!        the extinction coefficient is converted to kilometer**(-1)     
!        the unit utilized by the shortwave routine Swresf.             
!                                                                       
!        a value of 1.0E-100 is added to the size so that no division   
!        by zero occurs when the size is zero, in defining the          
!        extinction coefficient.                                        
!----------------------------------------------------------------------
!                                                                      
! intent in:                                                           
!                                                                      
! conc_rain = the rain drop water concentration in grams / meter**3    
!                                                                      
! size_rain = the rain drop effective diameter in microns              
!----------------------------------------------------------------------
 
real, dimension (:,:,:), intent(in)     ::   conc_rain, size_rain
 
!----------------------------------------------------------------------
! intent out:                                                           
!                                                                       
! cldextivlrain   = the specified spectral values of the extinction     
!                   coefficient for rain in kilometers**(-1)            
!                                                                       
! cldssalbivlrain = the specified spectral values of the single-        
!                   scattering albedo for rain                          
!                                                                       
! cldasymmivlrain = the specified spectral values of the asymmetry      
!                   factor for rain                                     
!----------------------------------------------------------------------
 
real, dimension (:,:,:,:), intent(out)  ::  cldextivlrain,            &
                                            cldssalbivlrain,   &
					    cldasymmivlrain

!----------------------------------------------------------------------c
! local variables:                                                     c
!----------------------------------------------------------------------c
 
      integer                                 ::  k, j, i, n, ni
      real, dimension (NRAINCLDIVLS)          ::  a, asymm, b
      real, dimension (:,:,:), allocatable    ::  rcap, size
 
      data a     / 4.65E-01,2.64E-01,1.05E-02,8.00E-05 /
      data b     / 1.00E-03,9.00E-02,2.20E-01,2.30E-01 /
      data asymm / 9.70E-01,9.40E-01,8.90E-01,8.80E-01 /
 !--------------------------------------------------------------------
 
      allocate ( size(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
      allocate ( rcap(ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD) )
 
 !-----------------------------------------------------------------
 !  compute scattering parameters for rain drops.  convert input size 
 !  from diameter to radius.  bypass calculations if no drops present.
 !-----------------------------------------------------------------
      do k = KSRAD,KERAD
        do j = JSRAD,JERAD
          do i = ISRAD,IERAD
            if (conc_rain(i,j,k) == 0.0) then
              cldextivlrain(i,j,k,:) = 0.0
              cldasymmivlrain(i,j,k,:) = 0.0
              cldssalbivlrain(i,j,k,:) = 0.0
            else
              size(i,j,k) = 5.0E-01 * size_rain(i,j,k) 
              if (size(i,j,k) .gt. 1.66E+01 .and.              &
                  size(i,j,k) .le. 5.0E+03 ) then                       
                rcap(i,j,k) = (size(i,j,k)/5.0E+02) ** 4.348E+00
                do ni = 1,NRAINCLDIVLS
                  cldextivlrain(i,j,k,ni) = 1.00E+03*1.505E+00*     &
                                            conc_rain(i,j,k) /   &
		  			    (size(i,j,k) + 1.0E-100 )
                  cldssalbivlrain(i,j,k,ni) = 1.0E+00 - ( a(ni) *    &
                                              (rcap(i,j,k)**b(ni) ) )
                  cldasymmivlrain(i,j,k,ni) = asymm(ni)
	        end do
              else
		call error_mesg ('savijarvi', &
 	                         'rain drop size out of range', FATAL)
              endif
            endif
          end do
        end do
      end do
!---------------------------------------------------------------------
 
      deallocate ( rcap )
      deallocate ( size )


end subroutine savijarvi



!####################################################################

subroutine fu                                                 &
                    (conc_ice    , size_ice      ,                  &
                     cldextivlice, cldssalbivlice, cldasymmivlice)      
 
!----------------------------------------------------------------------
! define the single scattering parameters for ice crystals using the    
! Fu parameterization for his spectral intervals.                       
!                                                                       
! references:                                                           
!                                                                       
! fu, q., an accurate parameterization of the solar radiative           
!      properties of cirrus clouds for climate models., j. climate,     
!      9, 2058-2082, 1996.                                              
!                                                                       
! notes: the ice crystal effective size (D^sub^ge in his paper) can     
!        only be 18.6 <= D^sub^ge <= 130.2 microns.                     
!                                                                       
!        the single scattering properties for wavenumbers < 2000 cm-1   
!        are assigned the values in the first interval, since the       
!        formulation is not defined for those wavenumbers.              
!                                                                       
!        the extinction coefficient is converted to kilometer**(-1)     
!        the unit utilized by the shortwave routine Swresf.             
!                                                                       
!        a value of 1.0E-100 is added to the size so that no division   
!        by zero occurs when the size is zero, in defining the          
!        extinction coefficient.                                        
!----------------------------------------------------------------------
! intent in:                                                            
!                                                                       
! conc_ice = the ice water concentation in grams / meter**3             
!                                                                       
! size_ice = the ice crystal effective size in microns                  
!----------------------------------------------------------------------
 
real, dimension (:,:,:), intent(in)     ::   conc_ice, size_ice
 
!----------------------------------------------------------------------
! intent out:                                                           
!                                                                       
! cldextivlice   = the specified spectral values of the extinction      
!                  coefficient for ice particles in kilometers**(-1)    
!                                                                       
! cldssalbivlice = the specified spectral values of the single-         
!                  scattering albedo for ice particles                  
!                                                                       
! cldasymmivlice = the specified spectral values of the asymmetry       
!                  factor for ice particles                             
!----------------------------------------------------------------------
 
real, dimension (:,:,:,:), intent(out)  ::  cldextivlice,             &
                                            cldssalbivlice,    &
					    cldasymmivlice
 
!----------------------------------------------------------------------c
! local variables:                                                     c
!----------------------------------------------------------------------c
      integer                       ::  k, j, i, n, ni
      real, dimension (NICECLDIVLS) ::  a0fu, a1fu,             &
                                        b0fu, b1fu, b2fu, b3fu,       &
                                        c0fu, c1fu, c2fu, c3fu
 
      data a0fu / -2.54823E-04, 1.87598E-04, 2.97295E-04, 2.34245E-04, &
                   4.89477E-04,-8.37325E-05, 6.44675E-04,-8.05155E-04, &
                   6.51659E-05, 4.13595E-04,-6.14288E-04, 7.31638E-05, &
                   8.10443E-05, 2.26539E-04,-3.04991E-04, 1.61983E-04, &
                   9.82244E-05,-3.03108E-05,-9.45458E-05, 1.29121E-04, &
                  -1.06451E-04,-2.58858E-04,-2.93599E-04,-2.66955E-04, &
                  -2.36447E-04 /
      data a1fu /  2.52909E+00, 2.51396E+00, 2.48895E+00, 2.48573E+00, &
                   2.48776E+00, 2.52504E+00, 2.47060E+00, 2.57600E+00, &
                   2.51660E+00, 2.48783E+00, 2.56520E+00, 2.51051E+00, &
                   2.51619E+00, 2.49909E+00, 2.54412E+00, 2.50746E+00, &
                   2.50875E+00, 2.51805E+00, 2.52061E+00, 2.50410E+00, &
                   2.52684E+00, 2.53815E+00, 2.54540E+00, 2.54179E+00, &
                   2.53817E+00 /
      data b0fu /  2.60155E-01, 1.96793E-01, 4.64416E-01, 9.05631E-02, &
                   5.83469E-04, 2.53234E-03, 2.01931E-03,-2.85518E-05, &
                  -1.48012E-07, 6.47675E-06,-9.38455E-06,-2.32733E-07, &
                  -1.57963E-07,-2.75031E-07, 3.12168E-07,-7.78001E-08, &
                  -8.93276E-08, 9.89368E-08, 5.08447E-07, 7.10418E-07, &
                   3.25057E-08,-1.98529E-07, 1.82299E-07,-1.00570E-07, &
                  -2.69916E-07 /
      data b1fu/   5.45547E-03, 5.75235E-03, 2.04716E-05, 2.93035E-03, &
                   1.18127E-03, 1.75078E-03, 1.83364E-03, 1.71993E-03, &
                   9.02355E-05, 2.18111E-05, 1.77414E-05, 6.41602E-06, &
                   1.72475E-06, 9.72285E-07, 4.93304E-07, 2.53360E-07, &
                   1.14916E-07, 5.44286E-08, 2.73206E-08, 1.42205E-08, &
                   5.43665E-08, 9.39480E-08, 1.12454E-07, 1.60441E-07, &
                   2.12909E-07 /
      data b2fu / -5.58760E-05,-5.29220E-05,-4.60375E-07,-1.89176E-05, &
                  -3.40011E-06,-8.00994E-06,-7.00232E-06,-7.43697E-06, &
                  -1.98190E-08, 1.83054E-09,-1.13004E-09, 1.97733E-10, &
                   9.02156E-11,-2.23685E-10, 1.79019E-10,-1.15489E-10, &
                  -1.62990E-10,-1.00877E-10, 4.96553E-11, 1.99874E-10, &
                  -9.24925E-11,-2.54540E-10,-1.08031E-10,-2.05663E-10, &
                  -2.65397E-10 /
      data b3fu /  1.97086E-07, 1.76618E-07, 2.03198E-09, 5.93361E-08, &
                   8.78549E-09, 2.31309E-08, 1.84287E-08, 2.09647E-08, &
                   4.01914E-11,-8.28710E-12, 2.37196E-12,-6.96836E-13, &
                  -3.79423E-13, 5.75512E-13,-7.31058E-13, 4.65084E-13, &
                   6.53291E-13, 4.56410E-13,-1.86001E-13,-7.81101E-13, &
                   4.53386E-13, 1.10876E-12, 4.99801E-13, 8.88595E-13, &
                   1.12983E-12 /
      data c0fu /  7.99084E-01, 7.59183E-01, 9.19599E-01, 8.29283E-01, &
                   7.75916E-01, 7.58748E-01, 7.51497E-01, 7.52528E-01, &
                   7.51277E-01, 7.52292E-01, 7.52048E-01, 7.51715E-01, &
                   7.52318E-01, 7.51779E-01, 7.53393E-01, 7.49693E-01, &
                   7.52131E-01, 7.51135E-01, 7.49856E-01, 7.48613E-01, &
                   7.47054E-01, 7.43546E-01, 7.40926E-01, 7.37809E-01, &
                   7.33260E-01 /
      data c1fu /  4.81706E-03, 4.93765E-03, 5.03025E-04, 2.06865E-03, &
                   1.74517E-03, 2.02709E-03, 2.05963E-03, 1.95748E-03, &
                   1.29824E-03, 1.14395E-03, 1.12044E-03, 1.10166E-03, &
                   1.04224E-03, 1.03341E-03, 9.61630E-04, 1.05446E-03, &
                   9.37763E-04, 9.09208E-04, 8.89161E-04, 8.90545E-04, &
                   8.86508E-04, 9.08674E-04, 8.90216E-04, 8.97515E-04, &
                   9.18317E-04 /
      data c2fu / -5.13220E-05,-4.84059E-05,-5.74771E-06,-1.59247E-05, &
                  -9.21314E-06,-1.17029E-05,-1.12135E-05,-1.02495E-05, &
                  -4.99075E-06,-3.27944E-06,-3.11826E-06,-2.91300E-06, &
                  -2.26618E-06,-2.13121E-06,-1.32519E-06,-2.32576E-06, &
                  -9.72292E-07,-6.34939E-07,-3.49578E-07,-3.44038E-07, &
                  -2.59305E-07,-4.65326E-07,-1.87919E-07,-2.17099E-07, &
                  -4.22974E-07 /
      data c3fu /  1.84420E-07, 1.65801E-07, 2.01731E-08, 5.01791E-08, &
                   2.15003E-08, 2.95195E-08, 2.73998E-08, 2.35479E-08, &
                   6.33757E-09,-2.42583E-10,-5.70868E-10,-1.37242E-09, &
                  -3.68283E-09,-4.24308E-09,-7.17071E-09,-3.58307E-09, &
                  -8.62063E-09,-9.84390E-09,-1.09913E-08,-1.10117E-08, &
                  -1.13305E-08,-1.05786E-08,-1.16760E-08,-1.16090E-08, &
                  -1.07976E-08 /
 
 !-----------------------------------------------------------------
 !  compute scattering parameters for ice crystals. bypass calculations
 !  if no crystals are present.
 !-----------------------------------------------------------------
      do k = KSRAD,KERAD
        do j = JSRAD,JERAD
          do i = ISRAD,IERAD
            if (conc_ice (i,j,k) == 0.0) then
              cldextivlice (i,j,k,:) = 0.0
              cldasymmivlice (i,j,k,:) = 0.0
              cldssalbivlice (i,j,k,:) = 0.0
            else
	      if (size_ice(i,j,k).ge.1.86E+01 .and.                    &
                  size_ice(i,j,k).le.1.302E+02 ) then                   
                do ni = 1,NICECLDIVLS
                  cldextivlice(i,j,k,ni) = 1.0E+03*conc_ice(i,j,k) *  &
                                           (a0fu(ni) + ( a1fu(ni) /    &
                                            size_ice(i,j,k) + 1.0E-100))
                  cldssalbivlice(i,j,k,ni) = 1.0 -    &
					     (b0fu(ni) + b1fu(ni)* &
                                             size_ice(i,j,k) +    &
					     b2fu(ni)*  &
	   				     size_ice(i,j,k) ** 2 +    &
					     b3fu(ni)*   &
					     size_ice(i,j,k)**3)
                  cldasymmivlice(i,j,k,ni) = c0fu(ni) + c1fu(ni) * &
                                             size_ice(i,j,k) +    &
					     c2fu(ni)* &
                                             size_ice(i,j,k) ** 2 +   &
	  				     c3fu(ni)*   &
					     size_ice(i,j,k) ** 3
	        end do
	      else
	        call error_mesg ('fu', &
		   'ice crystal size out of range', FATAL)
              endif
            endif
          end do
        end do
      end do
!---------------------------------------------------------------------
 

end subroutine fu


!#####################################################################

subroutine snowsw                                            &
                    (conc_snow,                                    &
                     cldextivlsnow, cldssalbivlsnow, cldasymmivlsnow)      
 
!----------------------------------------------------------------------
! define the single scattering parameters for snow using the Fu         
! parameterization for his spectral intervals.                          
!                                                                       
! author: leo donner, gfdl, 11 Sept 98                                  
!                                                                       
! references:                                                           
!                                                                       
! fu, q., et al., (See notes from Kuo-Nan Liou, 1 Sept 98). (SNOW)      
!                                                                       
! notes: the single scattering properties for wavenumbers < 2500 cm-1   
!        are assigned the values in the first interval, since the       
!        formulation is not defined for those wavenumbers.              
!                                                                       
!        the extinction coefficient is in units of kilometer**(-1)      
!----------------------------------------------------------------------
! intent in:                                                            
!                                                                       
! conc_snow = the snow concentration in grams / meter**3                
!----------------------------------------------------------------------
 
real, dimension (:,:,:), intent(in)     ::   conc_snow
 
!----------------------------------------------------------------------
! intent out:                                                           
!                                                                       
! cldextivlsnow   = the specified spectral values of the extinction     
!                   coefficient for snow in kilometers**(-1)            
!                                                                       
! cldssalbivlsnow = the specified spectral values of the single-        
!                   scattering albedo for snow                          
!                                                                       
! cldasymmivlsnow = the specified spectral values of the asymmetry      
!                   factor for snow                                     
!----------------------------------------------------------------------
 
real, dimension (:,:,:,:), intent(out)  :: cldextivlsnow,            &
                                           cldssalbivlsnow,     &
					   cldasymmivlsnow
 
!----------------------------------------------------------------------c
! local variables:                                                     c
!----------------------------------------------------------------------c
      integer                         ::  k, j, i, n, ni
      real, dimension (NSNOWCLDIVLS)  ::  asymm, ext, ssalb
      real                            ::  conc_ref=0.5
 
      data asymm / 9.6373E-01,9.8141E-01,9.7816E-01,9.6820E-01,      &
                   8.9940E-01,8.9218E-01 /
      data ext   / 8.3951E-01,8.3946E-01,8.3941E-01,8.3940E-01,      &
                   8.3940E-01,8.3939E-01 /
      data ssalb / 5.3846E-01,5.2579E-01,5.3156E-01,5.6192E-01,      &
                   9.7115E-01,9.99911E-01 /

!-----------------------------------------------------------------
!  compute scattering parameters for snow flakes. bypass calculations
!  if no snow flakes are present.
!-----------------------------------------------------------------
      do k=KSRAD,KERAD
        do j=JSRAD,JERAD
          do i=ISRAD,IERAD
            if (conc_snow(i,j,k) == 0.0) then
              cldextivlsnow(i,j,k,:) = 0.0
              cldasymmivlsnow(i,j,k,:) = 0.0
              cldssalbivlsnow(i,j,k,:) = 0.0
            else
              do ni = 1,NSNOWCLDIVLS
                cldextivlsnow(i,j,k,ni) = ext(ni)*conc_snow(i,j,k)/   &
					  conc_ref
                cldssalbivlsnow(i,j,k,ni) = ssalb(ni)
                cldasymmivlsnow(i,j,k,ni) = asymm(ni)
              end do
            endif
          end do
        end do
      end do
!--------------------------------------------------------------------
 

end subroutine snowsw


!######################################################################

subroutine snowlw
!subroutine snowlw(riwp,
!                        tau)
!
!-----------------------------------------------------------------------
!
!     Calculates emissivity
!     for longwave radiation using Fu et al. (1995,
!     JAS). (See notes from Kuo-Nan Liou, 1 Sept 98).
!
!-----------------------------------------------------------------------
!
!     On Input:
!
!        riwp   snow path  (g/(m**2))
!
!     On Output:
!
!        tau   absorption optical depth
!
!        Leo Donner, GFDL, 11 Sept 98
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!     Calculate extinction coefficient (m**(-1)) and optical depth
!     for absorption. Extinction coefficient is taken as product of
!     total extinction coefficient (.84 km**(-1)) and single-scattering
!     albedo. See Kuo-Nan Liou notes 
!     (1 Sept 98).
!
!-----------------------------------------------------------------------
!
!     riwc0=0.5
!     ext=.4
!     if (riwp .eq. 0.) then
!       tau=0.
!     else
!       tau=ext*.001*riwp/riwc0
!     end if
!     emis=1.-exp(-tau)
!-----------------------------------------------------------------------
end subroutine snowlw

!####################################################################

subroutine icesolar                                           &
                    (conc_ice    , size_ice      ,                  &
                     cldextivlice, cldssalbivlice, cldasymmivlice)      
 
!---------------------------------------------------------------------- 
! define the single scattering parameters for ice crystals using the    
! Fu parameterization for his spectral intervals.                       
!                                                                       
! references:                                                           
!                                                                       
! fu and Liou (1993, JAS) 
!                                                                       
! notes: the ice crystal effective size (D^sub^e in paper) can          
!        only be 18.6 <= D^sub^e <= 130.2 microns.                     
!                                                                       
!        the single scattering properties for wavenumbers < 2000 cm-1   
!        are assigned the values in the first interval, since the       
!        formulation is not defined for those wavenumbers.              
!                                                                       
!        the extinction coefficient is converted to kilometer**(-1)     
!        the unit utilized by the shortwave routine Swresf.             
!                                                                       
!        a value of 1.0E-100 is added to the size so that no division  
!        by zero occurs when the size is zero, in defining the          
!        extinction coefficient.                                        
!---------------------------------------------------------------------- 
! intent in:                                                            
!                                                                       
! conc_ice = the ice water concentation in grams / meter**3             
!                                                                       
! size_ice = the ice crystal effective size in microns                  
! Corresponds to minimum dimension of hexagonal crystal.
!                                                                       
!----------------------------------------------------------------------
 
real, dimension (:,:,:), intent(in)     ::   conc_ice, size_ice
 
!----------------------------------------------------------------------
! intent out:                                                           
!                                                                       
! cldextivlice   = the specified spectral values of the extinction      
!                  coefficient for ice particles in kilometers**(-1)    
!                                                                       
! cldssalbivlice = the specified spectral values of the single-         
!                  scattering albedo for ice particles                  
!                                                                       
! cldasymmivlice = the specified spectral values of the asymmetry       
!                  factor for ice particles                             
!                                                                       
!----------------------------------------------------------------------
 
real, dimension (:,:,:,:), intent(out)  ::  cldextivlice,             &
                                            cldssalbivlice, &
					    cldasymmivlice
 
!---------------------------------------------------------------------- 
! local variables:                                                      
!---------------------------------------------------------------------- 
      real, dimension (:,:,:,:), allocatable      :: fdel, fgam
      real, dimension (1:NICESOLARCLDIVLS, 0:NBB) :: b
      real, dimension (1:NICESOLARCLDIVLS, 0:NBC) :: c
      real, dimension (1:NICESOLARCLDIVLS, 0:NBD) :: d

      integer :: k, j, i, n, ni
      real    :: a0 = -6.656e-03
      real    :: a1 = 3.686

      data b/.10998e-05,.20208e-04,.1359e-03,-.16598e-02,.4618,      &
             .42362e-1, -.26101e-7,.96483e-05,.73453e-3,.20933e-2,   &
             .24471e-3,.86425e-2,.10896e-8,.83009e-7,.28281e-5,      &
            -.13977e-5,-.27839e-5,-.75519e-4,-.47387e-11,-.32217e-9, &
            -.18272e-7,-.18703e-7,.10379e-7,.24056e-6 /
      data c/2.211,2.2151,2.2376,2.3012,2.7975,1.9655,               &
         -.10398e-2,-.77982e-3,.10293e-2,.33854e-2,.29741e-2,        &
         .20094e-1,.65199e-4,.6375e-4,.50842e-4,.23528e-4,-.32344e-4,&
         -.17067e-3,-.34498e-6,-.34466e-6,-.30135e-6,-.20068e-6,     &
         .11636e-6,.50806e-6 /
      data d/.12495,.12363,.12117,.11581,-.15968e-3,.1383,           &
        -.43582e-3,-.44419e-3,-.48474e-3,-.55031e-3,.10115e-4,       &
        -.18921e-2,.14092e-4,.14038e-4,.12495e-4,.98776e-5,          &
        -.12472e-6,.1203e-4,-.69565e-7,-.68851e-7,-.62411e-7,        &
        -.50193e-7,.48667e-9,-.31698e-7/

!----------------------------------------------------------------------
!     allocate allocatable arrays
!----------------------------------------------------------------------
      allocate ( fgam                  &
             (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NICESOLARCLDIVLS))
      allocate ( fdel        &                        
             (ISRAD:IERAD, JSRAD:JERAD, KSRAD:KERAD, NICESOLARCLDIVLS))
 
 !-----------------------------------------------------------------
 !  compute scattering parameters for ice crystals. bypass calculations
 !  if no crystals are present.
 !-----------------------------------------------------------------
      do k = KSRAD,KERAD
        do j = JSRAD,JERAD
          do i = ISRAD,IERAD
            if (conc_ice (i,j,k) == 0.0) then
              cldextivlice (i,j,k,:) = 0.0
              cldasymmivlice (i,j,k,:) = 0.0
              cldssalbivlice (i,j,k,:) = 0.0
            else
	      if (size_ice(i,j,k).ge.1.86E+01 .and.                    &
                  size_ice(i,j,k).le.1.302E+02 ) then                   
                do ni = 1,NICESOLARCLDIVLS
                  cldextivlice(i,j,k,ni) = 1.0E+03 * conc_ice(i,j,k) *&
                                           (a0 + (a1/size_ice(i,j,k) + &
					   1.0E-100 ) )
                  cldssalbivlice(i,j,k,ni) = 1.0 - (b(7-ni,0) +   &
					     b(7-ni,1) *    &
                                             size_ice(i,j,k) +    &
					     b(7-ni,2) *      &
                                             size_ice(i,j,k) ** 2 + &
					     b(7-ni,3) * &
                                             size_ice(i,j,k) ** 3 )
                  fgam(i,j,k,ni)           = c(7-ni,0) + c(7-ni,1) *  &
                                             size_ice(i,j,k) +   &
					     c(7-ni,2) *      &
                                             size_ice(i,j,k) ** 2 +  &
					     c(7-ni,3) * &
                                             size_ice(i,j,k) ** 3
                  fdel(i,j,k,ni)           = d(7-ni,0) + d(7-ni,1)*  &
                                             size_ice(i,j,k) +   &
					     d(7-ni,2) *      &
                                             size_ice(i,j,k) ** 2 +  &
					     d(7-ni,3) * &
                                             size_ice(i,j,k) ** 3
	          cldasymmivlice(i,j,k,ni) =                 &
	                                     ((1.-fdel(i,j,k,ni))*  &
					     fgam(i,j,k,ni) +    &
                                             3.*fdel(i,j,k,ni))/3.
                end do
              else
	        call error_mesg ('icesolar',   &
		   'ice crystal size out of range', FATAL)
              endif
            endif
          end do
        end do
      end do

!-------------------------------------------------------------------
!  deallocate  arrays
!------------------------------------------------------------------
      deallocate (fgam )
      deallocate (fdel )

!------------------------------------------------------------------



end subroutine icesolar


!###################################################################


	       end module microphys_rad_mod

