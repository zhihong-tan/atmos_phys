module cloud_generator_mod

!   shared modules:
  use sat_vapor_pres_mod, only: lookup_es, lookup_des
  use constants_mod,      only: hlv, hls, cp_air, tfreeze, &
                                rvgas, rdgas
  use fms_mod,            only: open_namelist_file, mpp_pe,       &
                                mpp_root_pe, stdlog,              &
                                write_version_number, file_exist, &
                                check_nml_error, error_mesg,      &
                                FATAL, close_file
  use random_numbers_mod, only: randomNumberStream,        &
                                getRandomNumbers
  use beta_dist_mod,      only: beta_dist_init, beta_dist_end, &
                                incomplete_beta, beta_deviate
!--------------------------------------------------------------------
  !
  ! Given a profile of cloud fraction, produce a set of columns indicating 
  !   the presence or absence of cloud consistent with on of four overlap  
  !   assumptions: random, maximum, maximum-random, and one allowing 
  !   the rank correlation of the variable to be specified between layers.
  ! The module uses a random number generation module which can be used to wrap
  !   an arbitrary random number generator. The module defines a type that keeps 
  !   the state (the seed, a generator indicator, or whatever).  
  ! Each function takes cloud fraction as a function of height. You can either supply 
  !   a vector or a 3D array of dimensions (nX, ny, nLevels); in either case the 
  !   first element in the level dimension is the highest model layer. 
  ! Each function returns an array nSamples by nLevels long, where nLevels
  !   is determined by the size of the cloud fraction array. 
  ! The neighbor-to-neighbor correlation routine takes an additional 
  !   parameters described below. 
  !
!--------------------------------------------------------------------
  implicit none
  private

!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id: cloud_generator.F90,v 11.0 2004/09/28 19:14:22 fms Exp $'
character(len=128)  :: tagname =  '$Name: khartoum $'

!---------------------------------------------------------------------
!-------  interfaces --------

  interface genRandomOverlapSamples
    module procedure genRandomOverlapSamples_1D, genRandomOverlapSamples_3D
  end interface ! genRandomOverlapSamples
  
  interface genMaximumOverlapSamples
    module procedure genMaximumOverlapSamples_1D, genMaximumOverlapSamples_3D
  end interface ! genMaximumOverlapSamples
  
  interface genMaxRanOverlapSamples
    module procedure genMaxRanOverlapSamples_1D, genMaxRanOverlapSamples_3D
  end interface ! genMaxRanOverlapSamples
  
  interface genWeightedOverlapSamples
    module procedure genWeightedOverlapSamples_1D, genWeightedOverlapSamples_3D
  end interface ! genWeightedOverlapSamples
  
  public :: cloud_generator_init, &
            cloud_generator_end,  &
            generate_stochastic_clouds, &
            do_cloud_generator,   &
            compute_overlap_weighting
  
!---------------------------------------------------------------------
!-------- namelist  ---------

  ! Minimum values for cloud fraction, water, ice contents
  !   Taken from cloud_rad. Perhaps these should be namelist parameters? 
  real, parameter :: qmin = 1.E-10, qamin = 1.E-2
  
  ! Pressure scale height - for converting pressure differentials to height
  real, parameter :: pressureScaleHeight = 7.3 ! km
  
  ! Overlap parameter: 1 - Maximum, 2 - Random, 3 - Maximum/Random
  integer         :: defaultOverlap = 2 
  real            :: overlapLengthScale = 2.0  ! km
  ! These control the option to pull cloud condensate from a symmetric  
  !    beta distribution with p = q = betaP, adjusting 
  !    the upper and lower bounds of the ditribution to match 
  !    cloud fraction and cloud condensate. 
  logical         :: do_inhomogeneous_clouds = .false. 
  integer         :: betaP = 5
  
  namelist /cloud_generator_nml/  defaultOverlap, overlapLengthScale, &
                                  do_inhomogeneous_clouds, betaP 

!----------------------------------------------------------------------
!----  private data -------

logical :: module_is_initialized = .false.  ! module is initialized ?
logical :: cloud_generator_on    = .false.  ! is module being operated?

!----------------------------------------------------------------------

                              contains 
                              
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              
!######################################################################
subroutine cloud_generator_init


!---------------------------------------------------------------------
!    cloud_generator_init is the constructor for 
!    cloud_generator_mod.

!----------------------------------------------------------------------
!   local variables:
      integer   ::   unit, ierr, io

!--------------------------------------------------------------------
!   local variables:
!
!      unit     io unit for reading nml file and writing logfile
!      ierr     error code
!      io       error status returned from io operation  
!
!--------------------------------------------------------------------

      if (.not. module_is_initialized) then
!---------------------------------------------------------------------
!    read namelist.         
!---------------------------------------------------------------------
        if (file_exist('input.nml')) then
          unit =  open_namelist_file ( )
          ierr=1; do while (ierr /= 0)
          read (unit, nml=cloud_generator_nml, iostat=io, end=10) 
          ierr = check_nml_error (io, 'cloud_generator_nml')
          enddo                       
10        call close_file (unit)      
        endif                         
!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
        call write_version_number (version, tagname)
        if (mpp_pe() == mpp_root_pe() ) &
                   write (stdlog(), nml=cloud_generator_nml)
                   
!---------------------------------------------------------------------
!    Initialize the beta distribution module if we're going to need it. 
!---------------------------------------------------------------------
        if(do_inhomogeneous_clouds) call beta_dist_init

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
        module_is_initialized = .true.
        cloud_generator_on    = .true.
     end if

end subroutine cloud_generator_init
!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine generate_stochastic_clouds(streams, ql, qi, qa,         &
                                      overlap, pFull, temperature, &
                                      cld_thickness, ql_stoch, qi_stoch, qa_stoch)
!--------------------------------------------------------------------
!   intent(in) variables:
!
  type(randomNumberStream), &
           dimension(:, :),     intent(inout) :: streams
  ! Dimension nx, ny, nz
  real,    dimension(:, :, :),    intent( in) :: ql, qi, qa
  integer,                     optional, &
                                  intent( in) :: overlap
  real,    dimension(:, :, :), optional, &
                                 intent( in)  :: pFull, temperature
  ! Dimension nx, ny, nz, nCol = nBands
  integer, dimension(:, :, :, :), intent(out) :: cld_thickness 
  real,    dimension(:, :, :, :), intent(out) :: ql_stoch, &
                                                 qi_stoch, qa_stoch
  ! ---------------------------------------------------------
  ! Local variables
  real,    dimension(size(ql_stoch, 1), &
                     size(ql_stoch, 2), &
                     size(ql_stoch, 3)) :: qa_local, ql_local, qi_local
  real,    dimension(size(ql_stoch, 1), &
                     size(ql_stoch, 2), &
                     size(ql_stoch, 3)) :: heightDifference, &
                                           overlapWeighting  ! 1 for max, 0 for random
                         
  real,    dimension(size(ql_stoch, 1), &
                     size(ql_stoch, 2), &
                     size(ql_stoch, 3), &
                     size(ql_stoch, 4)) :: pdfRank ! continuously from 0 to 1
                                  
  logical, dimension(size(ql_stoch, 1), &
                     size(ql_stoch, 2), &
                     size(ql_stoch, 3), &
                     size(ql_stoch, 4)) :: isCloudy

  ! These arrays could be declared allocatable and used only when 
  !    do_inhomogeneous_clouds is true. 
  real,    dimension(size(ql_stoch, 1), &  ! Quantities for estimating condensate variability
                     size(ql_stoch, 2), &  !   from a beta distribution
                     size(ql_stoch, 3)) :: aThermo, qlqcRatio, qs_norm, deltaQ 

  real,    dimension(size(ql_stoch, 1), &
                     size(ql_stoch, 2), &
                     size(ql_stoch, 3), &
                     size(ql_stoch, 4)) :: qc_stoch
                     
  integer :: nLev, nCol, lev
  integer :: overlapToUse
  integer :: i, j, k, n
 
  ! ---------------------------------------------------------

  nLev = size(ql_stoch, 3)
  nCol = size(ql_stoch, 4)
  
  !
  ! Normally, we use the overlap specified in the namelist for this module, but
  !   this can be overridden with an optional argument. 
  !
  if(present(overlap)) then
    overlapToUse = overlap
  else
    overlapToUse = defaultOverlap
  end if
  
  !
  ! Ensure that cloud fraction, water, and ice contents are in bounds
  !   After similar code in cloud_summary3
  !
  qa_local(:,:,:) = 0.
  do k=1, nLev
       do j=1, size(qa_stoch,2)
       do i=1, size(qa_stoch,1)
         if (qa(i,j,k) >= qamin .and. ql(i,j,k) > qmin) then
           qa_local(i,j,k) = qa(i,j,k)
           ql_local(i,j,k) = ql(i,j,k)
         else
           ql_local(i,j,k) = 0.0           
         endif
         if (qa(i,j,k) >= qamin .and. qi(i,j,k) >= qmin) then
           qa_local(i,j,k) = qa(i,j,k)
           qi_local(i,j,k) = qi(i,j,k)
         else
           qi_local(i,j,k) = 0.0           
         endif
       end do
       end do
       end do

  
  !
  ! Apply overlap assumption
  !
   if (overlapToUse == 2) then
      pdfRank(:, :, :, :) = genRandomOverlapSamples( qa_local(:, :, :), nCol, streams(:, :))
   else if (overlapToUse == 1) then
      pdfRank(:, :, :, :) = genMaximumOverlapSamples(qa_local(:, :, :), nCol, streams(:, :))
   else if (overlapToUse == 3) then
      pdfRank(:, :, :, :) = genMaxRanOverlapSamples( qa_local(:, :, :), nCol, streams(:, :))
   else if (overlapToUse == 4) then
      if(.not. present(pFull)) call error_mesg("cloud_generator_mod", &
                                               "Need to provide pFull when using overlap = 4", FATAL)
      !
      ! Height difference from hydrostatic equation with fixed scale height for T = 250 K
      !
      heightDifference(:, :, :nLev-1) = (log(pFull(:, :, 2:nLev)) - log(pFull(:, :, 1:nLev-1))) * &
                                         pressureScaleHeight 
      heightDifference(:, :, nLev) = heightDifference(:, :, nLev-1)
      !
      ! Overlap is weighted between max and random with parameter overlapWeighting (0 = random, 
      !    1 = max), which decreases exponentially with the separation distance. 
      !
      overlapWeighting(:, :, :) = exp( - heightDifference(:, :, :) / overlapLengthScale )
      pdfRank(:, :, :, :) = genWeightedOverlapSamples(qa_local(:, :, :),         &
                                                      overlapWeighting(:, :, :), & 
                             nCol, streams(:, :))
    else
      call error_mesg("cloud_generator_mod", "unknown overlap parameter", FATAL)
     
     endif



  
  if(.not. do_inhomogeneous_clouds) then
    !
    ! The clouds are uniform, so every cloudy cell at a given level 
    !   gets the same ice and/or liquid concentration (subject to a minumum). 
    ! We're looking for in-cloud ice and water contents, so we 
    !   divide by cloud fraction
    !  
    do n=1,nCol
      do k=1,nLev
       do j=1, size(qa_stoch,2)
       do i=1, size(qa_stoch,1)
!  a "true" for the following test indicates the presence of cloudiness
         if (pdfRank(i,j,k,n) > (1.-qa_local(i,j,k))) then
           cld_thickness(i,j,k,n   ) = 1.
           qa_stoch(i,j,k,n   ) = 1. 
           ql_stoch(i,j,k,n   ) = ql_local(i,j,k)/qa_local(i,j,k)
           qi_stoch(i,j,k,n   ) = qi_local(i,j,k)/qa_local(i,j,k)
         else
           cld_thickness(i,j,k,n   ) = 0.
           qa_stoch(i,j,k,n   ) = 0. 
           ql_stoch(i,j,k,n   ) = 0.                              
           qi_stoch(i,j,k,n   ) = 0.                              
         endif
        end do
        end do
        end do
        end do
           

  else
    if(.not. present(pFull) .or. .not. present(temperature)) &
      call error_mesg("cloud_generator_mod",                 &
      "Need to provide pFull and temperature when using inhomogenous clouds", FATAL)
    !
    ! Assume that total water in each grid cell follows a symmetric beta distribution with  
    !   exponents p=q set in the namelist. Determine the normalized amount of condensate 
    !   (qc - qmin)/(qmax - qmin) from the incomplete beta distribution at the pdf rank. 
    !   Convert to physical units based on three equations
    !   qs_norm = (qs - qmin)/(qmax - qmin)
    !   qa = 1. - betaIncomplete(qs_norm, p, q)
    !   qc_mean/(qmax - qmin) = aThermo( p/(p+q) (1 - betaIncomplete(qs_norm, p+1, q)) - 
    !                                    qs_norm * (1 - cf) )
    !   The latter equation comes from integrating aThermo * (qtot - qsat) times a beta distribution 
    !   from qsat to qmax; see, for example, from Tompkins 2002 (Eq. 14), but including a thermodynamic 
    !   term aThermo = 1/(1 + L/cp dqs/dt) evaluated at the "frozen temperature", as per SAK. 
    !   
    where (qa_local(:, :, :) > qamin) 
      qlqcRatio(:, :, :) = ql_local(:, :, :) / (ql_local(:, :, :) + qi_local(:, :, :))
    elsewhere
      qlqcRatio(:, :, :) = 1 ! Shouldn't be used. 
    end where
    call computeThermodynamics(temperature, pFull, ql, qi, aThermo)

    qs_norm(:, :, :) = 1. ! This assignment is made so the values of qs_norm 
                          ! are always valid; in practice these should be masked out 
  ! in the regions below. 
    where(qa_local(:, :, :) < qamin)
      !
      ! Not cloudy, so deltaQ is irrelevant
      !
      qs_norm(:, :, :) = 1. 
      deltaQ(:, :, :) = 0.
    elsewhere
      !
      ! Hey, is this the right test for fully cloudy conditions? Is cloud fraction ever 1.? 
      !
      where (qa_local(:, :, :) >= 1.) 
        !
        ! In fully cloudy conditions we don't have any information about the bounds
        !   of the total water distribution. We arbitrarily set the lower bound to qsat.
        !   For a symmetric distribution the upper bound to qsat plus twice the amount of 
        !   condensate. 
        !
          qs_norm(:, :, :) = 0.
        deltaQ(:, :, :) = (2/aThermo(:, :, :)) * &
                          (ql_local(:, :, :) + qi_local(:, :, :)) / qa_local(:, :, :)
        ! qMin(:, :, :) = qsat(:, :, :)
      elsewhere 
        !
        ! Partially cloudy conditions - diagnose width of distribution from mean 
        !   condensate amount and cloud fraction. 
        !   The factor 1/2 = p/(p+q)
        !
        qs_norm(:, :, :) = beta_deviate(1. - qa_local(:, :, :), p = betaP, q = betaP)
        deltaQ(:, :, :) =                                                               &
          (ql_local(:, :, :) + qi_local(:, :, :)) /                                     &
          (aThermo(:, :, :) * ((1./2. * (1. - incomplete_beta(qs_norm(:, :, :),         &
                                                          p = betaP + 1, q = betaP))) - &
                               qs_norm(:, :, :) * qa_local(:, :, :) ))
        ! qMin(:, :, :) = qsat(:, :, :) - qs_norm(:, :, :) * deltaQ(:, :, :)
      end where 
    end where
  
    do n=1,nCol
      do k=1,nLev
       do j=1, size(qa_stoch,2)
       do i=1, size(qa_stoch,1)
!  a "true" for the following test indicates the presence of cloudiness
         if (pdfRank(i,j,k,n) > (1.-qa_local(i,j,k))) then
           cld_thickness(i,j,k,n   ) = 1.
           qa_stoch(i,j,k,n   ) = 1. 
      !
      ! Hey, do we need to account for cloud fraction here, as we do in the homogeneous case? 
      !   Also, does this seem like the right way to go from the mean to individual samples? 
      !
      qc_stoch(i, j, k, n) = aThermo(i,j,k  ) * deltaQ(i,j,k  ) * &
            (beta_deviate(pdfRank(i,j,k,n   ), p = betaP, q = betaP) - &
                   qs_norm(i,j,k) ) 
      ! 
      ! The proportion of ice and water in each sample is the same as the mean proportion.  
      !
           ql_stoch(i,j,k,n   ) = qc_stoch(i,j,k,n)* qlqcRatio(i,j,k)
           qi_stoch(i,j,k,n   ) = qc_stoch(i,j,k,n)*   &
                                     (1.-qlqcRatio(i,j,k))
         else
           cld_thickness(i,j,k,n   ) = 0.
           qa_stoch(i,j,k,n   ) = 0. 
           ql_stoch(i,j,k,n   ) = 0.                              
           qi_stoch(i,j,k,n   ) = 0.                              
         endif
        end do
        end do
        end do
        end do
           
  end if 
end subroutine generate_stochastic_clouds

  ! ---------------------------------------------------------
  !  Function to return the weighting between maximum and random overlap 
  !    given the pressure difference 
  ! 
  !  Note pPlus is the pressure at a higher altitude
  !  (i.e. pPlus < pMinus)
  ! ---------------------------------------------------------
function compute_overlap_weighting(qaPlus, qaMinus, pPlus, pMinus) result(weighting)
  real, dimension(:, :), intent( in) :: qaPlus, qaMinus, pPlus, pMinus
  real, dimension(size(pPlus, 1), &
                  size(pPlus, 2))    :: weighting
        
  select case(defaultOverlap)
    case(1) ! Maximum overlap
      weighting(:, :) = 1.
    case(2) ! Random overlap
      weighting(:, :) = 0.
    case(3) ! Maximum-random
      where(qaPlus(:, :) > qamin) 
        weighting(:, :) = 1.
      elsewhere
        weighting(:, :) = 0.
      end where
    case(4)
      !
      ! Overlap is weighted between max and random with parameter overlapWeighting (0 = random, 
      !    1 = max), which decreases exponentially with the separation distance. 
      !
      weighting(:, :) = exp(-abs(log(pMinus(:, :)) - log(pPlus(:, :))) * &
                              pressureScaleHeight / overlapLengthScale)
    case default 
      call error_mesg("cloud_generator_mod", "unknown overlap parameter", FATAL)
    end select
      
end function compute_overlap_weighting

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!  These generate cloud samples according to specific overlap rules.  
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine computeThermodynamics(temperature, pressure, ql, qi, aThermo)
    real, dimension(:, :, :), intent( in) :: temperature, pressure, ql, qi
    real, dimension(:, :, :), intent(out) :: aThermo
    
    integer :: lev, nLev
    real, parameter :: d608 = (rvgas - rdgas)/rdgas, &
                       d622 = rdgas/rvgas,           &
                       d378 = 1. - d622
    
    ! Compute saturation mixing ratio and gamma = 1/(1 + L/cp dqs/dT) evaluated
    !   at the ice water temperature
    ! Taken from strat_cloud_mod
    
    ! Local variables
    real, dimension(size(temperature, 1), &
                    size(temperature, 2), &
                    size(temperature, 3)) :: Tl, L, esat, desdT, dqsdT
    
    !
    ! Ice water temperature - ql and qi are grid cell means
    !
    Tl(:, :, :) =  temperature(:, :, :) -       &
                   (hlv/cp_air) * ql(:, :, :) - &
                   (hls/cp_air) * qi(:, :, :)
  
    !calculate water saturated vapor pressure and its derivative from table
    call lookup_es( Tl(:, :, :),  esat(:, :, :))
    call lookup_des(Tl(:, :, :), desdT(:, :, :))
 
    !calculate dqsdT
    !limit denominator to esat, and thus qs to d622
    esat(:, :, :) = max(pressure(:, :, :) - d378 * esat(:, :, :), esat(:, :, :))
    !this is done to avoid blow up in the upper stratosphere
    dqsdT(:, :, :) = d622 * pressure(:, :, :) * desdT(:, :, :) / esat(:, :, :)**2
 
    ! Latent heat of phase change, varying from that of water to ice with temperature.
    ! 
    L(:, :, :) = (min(1., max(0., 0.05*(temperature(:, :, :) - tfreeze + 20.)))*hlv + &
                  min(1., max(0., 0.05*(tfreeze - temperature(:, :, :)      )))*hls)
    aThermo(:, :, :) = 1./ (1. + L(:, :, :)/cp_air * dqsdT(:, :, :)) 
  end subroutine computeThermodynamics
  ! ---------------------------------------------------------
  ! Random overlap - the value is chosen randomly from the distribution at 
  !   every height. 
  ! ---------------------------------------------------------
  function genRandomOverlapSamples_1D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:),    intent(in   ) :: cloudFraction
    integer,                  intent(in   ) :: nSamples
    type(randomNumberStream), intent(inout) :: stream
    real,    dimension(size(cloudFraction), nSamples) &
                                            :: randomValues
    ! -------------------
    call getRandomNumbers(stream, randomValues)
  end  function genRandomOverlapSamples_1D
  ! ---------------------------------------------------------
  function genRandomOverlapSamples_3D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:, :, :),    intent(in   ) :: cloudFraction
    integer,                        intent(in   ) :: nSamples
    type(randomNumberStream), &
             dimension(:, :),       intent(inout) :: stream
    real,    dimension(size(cloudFraction, 1), &
                       size(cloudFraction, 2), &
                       size(cloudFraction, 3), &
                       nSamples)                  :: randomValues
    ! Local variables
    integer :: i, j
    ! -------------------
    do j = 1, size(cloudFraction, 2)
      do i = 1, size(cloudFraction, 1)
        call getRandomNumbers(stream(i, j), randomValues(i, j, :, :))
      end do
    end do
  end  function genRandomOverlapSamples_3D
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  ! Maximum overlap - the position in the PDF is the same 
  !   at every height in a given column (though it varies from 
  !   column to column). 
  ! ---------------------------------------------------------

  function genMaximumOverlapSamples_1D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:),    intent(in   ) :: cloudFraction
    integer,                  intent(in   ) :: nSamples
    type(randomNumberStream), intent(inout) :: stream
    real,     dimension(size(cloudFraction), nSamples) &
                                            :: randomValues
    
    
    ! -------------------
    call getRandomNumbers(stream, randomValues(1, :))
    randomValues(:, :)  = spread(randomValues(1, :), &
                                 dim = 1, nCopies = size(cloudFraction))
  end  function genMaximumOverlapSamples_1D
  ! ---------------------------------------------------------
  function genMaximumOverlapSamples_3D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:, :, :),    intent(in   ) :: cloudFraction
    integer,                        intent(in   ) :: nSamples
    type(randomNumberStream), &
             dimension(:, :),       intent(inout) :: stream
    real,    dimension(size(cloudFraction, 1), &
                       size(cloudFraction, 2), &
                       size(cloudFraction, 3), &
                       nSamples)                  :: randomValues
    ! -------------------
    ! Local variables
    integer :: i, j, nX, nY, nLev, nCol
    ! -------------------
    nX   = size(cloudFraction, 1); nY  = size(cloudFraction, 2)
    nLev = size(cloudFraction, 3)
    
    do j = 1, nY
      do i = 1, nX
        call getRandomNumbers(stream(i, j), randomValues(i, j, 1, :))
      end do
    end do 
    randomValues(:, :, :, :)  = spread(randomValues(:, :, 1, :), dim = 3, nCopies = nLev)

  end  function genMaximumOverlapSamples_3D
  ! ---------------------------------------------------------
  
  ! ---------------------------------------------------------
  ! Meximum-random overlap. 
  ! Within each column, the value in the top layer is chosen 
  !   at random. We then walk down one layer at a time. If the layer above is cloudy
  !   we use the same random deviate in this layer; otherwise
  !   we choose a new value. 
  ! ---------------------------------------------------------
  
  function genMaxRanOverlapSamples_1D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:),    intent(in   ) :: cloudFraction
    integer,                  intent(in   ) :: nSamples
    type(randomNumberStream), intent(inout) :: stream
    real,    dimension(size(cloudFraction), nSamples) &
                                            :: randomValues
    ! Local variables
    integer                                 :: level
    
    ! -------------------
    call getRandomNumbers(stream, randomValues)
    do level = 2, size(cloudFraction)
      where(randomValues(:, level - 1) > 1. - cloudFraction(level - 1))
        randomValues(level, :) = randomValues(level - 1, :)
      elsewhere
        randomValues(level, :) = randomValues(level, :) * (1. - cloudFraction(level - 1))
      end where
    end do
  end  function genMaxRanOverlapSamples_1D
  ! ---------------------------------------------------------
  function genMaxRanOverlapSamples_3D(cloudFraction, nSamples, stream) &
           result(randomValues)
    real,    dimension(:, :, :),    intent(in   ) :: cloudFraction
    integer,                        intent(in   ) :: nSamples
    type(randomNumberStream), &
             dimension(:, :),       intent(inout) :: stream
    real,    dimension(size(cloudFraction, 1), &
                       size(cloudFraction, 2), &
                       size(cloudFraction, 3), &
                       nSamples)                  :: randomValues
    ! -------------------
    ! Local variables
    integer :: i, j, level, nX, nY, nLev, nCol
    ! -------------------
    nX   = size(cloudFraction, 1); nY  = size(cloudFraction, 2)
    nLev = size(cloudFraction, 3)
    
    do j = 1, nY
      do i = 1, nX
        call getRandomNumbers(stream(i, j), randomValues(i, j, :, :))
      end do
    end do 
              
    do level = 2, nLev
      where(randomValues(:, :, level - 1, :) > &
            spread(1. - cloudFraction(:, :, level - 1), dim = 3, nCopies = nSamples))
        randomValues(:, :, level, :) = randomValues(:, :, level - 1, :)
      elsewhere
        randomValues(:, :, level, :) = randomValues(:, :, level, :) * &
                                       spread(1. - cloudFraction(:, :, level - 1), &
                                              dim = 3, nCopies = nSamples)
      end where
    end do
  end  function genMaxRanOverlapSamples_3D
  ! ---------------------------------------------------------
  
  ! Neighbor to neighbor rank correlation, which gives exponential dependence 
  !   of rank correlation if the correlation is fixed. The correlation coefficient 
  !   is the array alpha. 
  ! Two streams of random numbers are generated. The first corresponds to the postion in 
  !   the PDF, and the second is used to enforce the correlation . If the value of 
  !   the second stream at one level in one column is less than alpha at that level, 
  !   the same relative position in the PDF is chosen in the lower layer as the upper.  
  ! ---------------------------------------------------------
  function genWeightedOverlapSamples_1D(cloudFraction, alpha, nSamples, stream) &
           result(randomValues)
    real,    dimension(:),    intent(in   ) :: cloudFraction, alpha
    integer,                  intent(in   ) :: nSamples
    type(randomNumberStream), intent(inout) :: stream
    real,    dimension(size(cloudFraction), nSamples) &
                                            :: randomValues
    
    ! Local variables
    real, dimension(size(cloudFraction), nSamples) :: randomValues2
    integer                                        :: level
    
    ! -------------------
    call getRandomNumbers(stream, randomValues)
    call getRandomNumbers(stream, randomValues2)
    
    do level = 1, size(cloudFraction) - 1 
      where(randomValues2(level + 1, :) < alpha(level)) &
        randomValues(level + 1, :) = randomValues(level, :)
    end do   
  end  function genWeightedOverlapSamples_1D
  ! ---------------------------------------------------------
  function genWeightedOverlapSamples_3D(cloudFraction, alpha, nSamples, stream) &
           result(randomValues)
    real,    dimension(:, :, :), intent(in   ) :: cloudFraction, alpha
    integer,                     intent(in   ) :: nSamples
    type(randomNumberStream), &
             dimension(:, :),    intent(inout) :: stream
    real,    dimension(size(cloudFraction, 1), &
                       size(cloudFraction, 2), &
                       size(cloudFraction, 3), &
                       nSamples)               :: randomValues
     
    ! Local variables
    real, dimension(size(cloudFraction, 1), &
                    size(cloudFraction, 2), &
                    size(cloudFraction, 3), &
                                  nSamples) :: randomValues2
    integer                                 :: i, j, nX, nY, nLev, level
    
    ! -------------------
    nX   = size(cloudFraction, 1); nY  = size(cloudFraction, 2)
    nLev = size(cloudFraction, 3)
    
    do j = 1, nY
      do i = 1, nX
        call getRandomNumbers(stream(i, j), randomValues (i, j, :, :))
      end do
    end do 
    do j = 1, nY
      do i = 1, nX
        call getRandomNumbers(stream(i, j), randomValues2(i, j, :, :))
      end do
    end do 
     
    do level = 1, nLev - 1
      where(randomValues2(:, :, level + 1, :) < spread(alpha(:, :, level),           &
                                                       dim = 3, nCopies = nSamples)) &
        randomValues(:, :, level + 1, :) = randomValues(:, :, level, :)
    end do
  end  function genWeightedOverlapSamples_3D
  ! ---------------------------------------------------------

  subroutine cloud_generator_end       
  !----------------------------------------------------------------------
  !    cloud_generator_end is the destructor for cloud_generator_mod.
  !----------------------------------------------------------------------
          
  !---------------------------------------------------------------------
  !    be sure module has been initialized.
  !---------------------------------------------------------------------
        if (.not. module_is_initialized ) then
          call error_mesg ('cloud_generator_mod',   &
               'module has not been initialized', FATAL )
        endif
        
        if(do_inhomogeneous_clouds) call beta_dist_end
  !---------------------------------------------------------------------
  !    mark the module as not initialized.
  !---------------------------------------------------------------------
        module_is_initialized = .false.
  end subroutine cloud_generator_end
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !
  !  Function to report if the cloud generator is being used. 
  !
  function do_cloud_generator()
    logical :: do_cloud_generator
    
    do_cloud_generator = cloud_generator_on
  end function do_cloud_generator
  !--------------------------------------------------------------------

end module cloud_generator_mod
