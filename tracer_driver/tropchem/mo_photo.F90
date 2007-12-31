      module MO_PHOTO_MOD
!----------------------------------------------------------------------
!        ... Photolysis interp table and related arrays
!----------------------------------------------------------------------
      use mpp_mod,          only : mpp_error, FATAL
      use mpp_io_mod,       only : mpp_open, MPP_RDONLY, MPP_ASCII,MPP_MULTI, &
                                   MPP_SINGLE, mpp_close
      use time_manager_mod, only : time_type, get_date
      use constants_mod,    only : PI

      implicit none

      private
      public :: prate_init, photo, set_ub_col, setcol, sundis

      save

      integer, parameter :: jdim     = 19
      integer, parameter :: altdim   = 18
      integer, parameter :: zangdim  = 8
      integer, parameter :: o3ratdim = 7
      integer, parameter :: albdim   = 4
      integer, parameter :: t500dim  = 3
      integer, parameter :: t200dim  = 2
      integer, parameter :: tabdim   = jdim*altdim*zangdim*o3ratdim*albdim*t500dim*t200dim

      integer ::  offset(7)
      integer ::  indexer(jdim)
      integer ::  jno_ndx, jpooh_ndx, jc2h5ooh_ndx, jc3h7ooh_ndx, jrooh_ndx, &
                  jch3co3h_ndx, jmpan_ndx, jmacr_a_ndx, jmacr_b_ndx, jonitr_ndx, &
                  jxooh_ndx, jisopooh_ndx, jglyald_ndx, jhyac_ndx, jch3ooh_ndx, &
                  jh2o2_ndx, jpan_ndx, jch3cho_ndx, &
                  jn2o5_ndx, jo3p_ndx, jno2_ndx, jno3_ndx, &
                  jclono2_ndx, jhocl_ndx, jcl2o2_ndx, jbrono2_ndx, jhobr_ndx, &
                  jbrcl_ndx, jbro_ndx, jcl2_ndx, jh2o_ndx
      integer ::  ox_ndx, o3_ndx
      real    ::  ajl(jdim,altdim,zangdim,o3ratdim,albdim,t500dim,t200dim) = 0.
      real    ::  vo3(0:50)
      real    ::  delvo3(0:49)
      real    ::  delz(altdim-1)
      real    ::  delang(zangdim-1)
      real    ::  delv(o3ratdim-1)
      real    ::  delalb(albdim-1)
      real    ::  delt500(t500dim-1)
      real    ::  delt200(t200dim-1)

      real :: zz(altdim) = &
             (/ 0., 1., 2., 3., 5., 7., 9., 12., 15., 18., 21., 24., &
                27., 30., 35., 40., 45., 50. /)
      real :: vsec(zangdim) = (/ 1., 1.3, 1.6, 2., 3., 6., 12., 50. /)
      real :: xv3(o3ratdim) = (/ .5, .75, 1., 1.25, 1.5, 2., 5. /)
      real :: albev(albdim) = (/ .05, .2, .5, 1. /)
      real :: t500(t500dim) = (/ 228., 248., 268. /)
      real :: t200(t200dim) = (/ 205., 225. /)

      integer, parameter :: &
         TAB_NDX_JO2     = 1, &
         TAB_NDX_JO1D    = 2, &
         TAB_NDX_JO3P    = 3, &
         TAB_NDX_JN2O    = 4, &
         TAB_NDX_JNO2    = 5, &
         TAB_NDX_JN2O5   = 6, &
         TAB_NDX_JHNO3   = 7, &
         TAB_NDX_JNO3    = 8, &
         TAB_NDX_JHO2NO2 = 9, &
         TAB_NDX_JCH3OOH = 10, &
         TAB_NDX_JCH2Oa  = 11, &
         TAB_NDX_JCH2Ob  = 12, &
         TAB_NDX_JH2O2   = 13, &
         TAB_NDX_JCH3CHO = 14, &
         TAB_NDX_JPAN    = 15, &
         TAB_NDX_JMACRa  = 16, &
         TAB_NDX_JMVK    = 17, &
         TAB_NDX_JACET   = 18, &
         TAB_NDX_JMGLY   = 19

character(len=128), parameter :: version     = '$Id: mo_photo.F90,v 14.0 2007/03/15 22:11:05 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: omsk_2007_12 $'
logical                       :: module_is_initialized = .false.

      CONTAINS
      
      subroutine PRATE_INIT( filename, lpath, mspath )
!----------------------------------------------------------------------
!     ... Read in the photorate tables and arrays
!         Results are "returned" via the common block photo_tables
!          This is for the new, expanded chemistry (11/21/94)
!----------------------------------------------------------------------
        
      use mo_chem_utls_mod,  only : get_spc_ndx, get_rxt_ndx
      implicit none

!----------------------------------------------------------------------
!        ... Dummy args
!----------------------------------------------------------------------
      character(len=*), intent(in) :: filename, lpath, mspath

!----------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------
      integer    :: it500, it200, izen, ialb, idob
      integer    :: ios, ofl, k
      integer    :: retval, unit
      logical    :: cosb
      real       :: temp(tabdim)
      character(len=128) :: msg

!     unit = navu()
!----------------------------------------------------------------------
!        ... Only masternode gets photorate table
!----------------------------------------------------------------------
!     if( masternode ) then
!        retval = ATTACH( TRIM(lpath) // TRIM( filename ), &
!                         TRIM( mspath ) // TRIM( filename ), &
!                         unit, &
!                         .false., &                 ! non binary dataset
!                         cosb )
!        if( retval /= 0 ) then
!           write(*,*) 'PRATE_INIT: Failure opening file ',TRIM(lpath) // TRIM( filename )
!           write(*,*) 'Error code = ',retval
!           call ENDRUN
!        end if
!        CLOSE( unit )
!     end if
#ifdef USE_MPI
!----------------------------------------------------------------------
!        ... All compute nodes wait for masternode to acquire file
!----------------------------------------------------------------------
!     call MPI_BARRIER( mpi_comm_comp, ios )
!     if( ios /= MPI_SUCCESS ) then
!        write(*,*) 'PRATE_INIT: Mpi barrier failed; error = ',ios
!        call ENDRUN
!     end if
#endif
!----------------------------------------------------------------------
!        ... all compute nodes open file
!----------------------------------------------------------------------
!     OPEN( unit   = unit, &
!           file   = TRIM(lpath) // TRIM(filename), &
!           status = 'old', &
!           form   = 'formatted', &
!           recl   = 4500, &
!          iostat = ios )
!     if( ios /= 0 ) then
!----------------------------------------------------------------------
!        ... Open error exit
!----------------------------------------------------------------------
!        write(6,'('' PRATE_INIT : Error ('',i5,'') opening file '',a)') &
!           ios, TRIM(lpath) // TRIM(filename)
!        call ENDRUN
!     end if

!----------------------------------------------------------------------
!            ... open file using mpp_open
!----------------------------------------------------------------------

      call mpp_open( unit, trim(filename), MPP_RDONLY, MPP_ASCII, &
                     threading = MPP_MULTI, fileset = MPP_SINGLE, &
                 recl = 4500)

!----------------------------------------------------------------------
!        ... Readin the reference o3 column and photorate table
!----------------------------------------------------------------------
      read(unit,*,iostat=ios) vo3
      if( ios /= 0 ) then
         msg = ' PRATE_INIT: Failed to read o3 column'
         call ENDRUN(msg)
      end if

      do it500 = 1,t500dim
         do it200 = 1,t200dim
            do izen = 1,zangdim
               do ialb = 1,albdim
                  do idob = 1,o3ratdim
                     read(unit,*,iostat=ios) ajl(:,:,izen,idob,ialb,it500,it200)
                     if( ios /= 0 ) then
                        msg = 'PRATE_INIT: Failed to read photo table; ' // &
                              'error = '//char(ios)
                        call ENDRUN(msg)
                     end if
                  end do
               end do
            end do
         end do
      end do

!----------------------------------------------------------------------
!        ... Set module variables
!----------------------------------------------------------------------
      delz(:altdim-1) = 1. / (zz(2:altdim) - zz(:altdim-1))
      delvo3(0:49) = vo3(1:50) - vo3(0:49)
      delang(:zangdim-1)  = 1. / (vsec(2:zangdim) - vsec(:zangdim-1))
      delv(:o3ratdim-1)   = 1. / (xv3(2:o3ratdim) - xv3(:o3ratdim-1))
      delalb(:albdim-1)   = 1. / (albev(2:albdim) - albev(:albdim-1))
      delt500(:t500dim-1) = 1. / (t500(2:t500dim) - t500(:t500dim-1))
      delt200(:t200dim-1) = 1. / (t200(2:t200dim) - t200(:t200dim-1))

      offset(1) = jdim
      offset(2) = offset(1)*altdim
      offset(3) = offset(2)*zangdim
      offset(4) = offset(3)*o3ratdim
      offset(5) = offset(4)*albdim
      offset(6) = offset(5)*t500dim
      offset(7) = SUM( offset(1:6) )

!     close( unit )
      call mpp_close( unit )

!-----------------------------------------------------------------
!           ... setup mapping array, indexer, from table to model
!-----------------------------------------------------------------
      indexer(TAB_NDX_JO2)     = get_rxt_ndx( 'jo2' )
      indexer(TAB_NDX_JO1D)    = get_rxt_ndx( 'jo1d' )
      indexer(TAB_NDX_JO3P)    = get_rxt_ndx( 'jo3p' )
      indexer(TAB_NDX_JN2O)    = get_rxt_ndx( 'jn2o' )
      indexer(TAB_NDX_JNO2)    = get_rxt_ndx( 'jno2' )
      indexer(TAB_NDX_JN2O5)   = get_rxt_ndx( 'jn2o5' )
      indexer(TAB_NDX_JHNO3)   = get_rxt_ndx( 'jhno3' )
      indexer(TAB_NDX_JNO3)    = get_rxt_ndx( 'jno3' )
      indexer(TAB_NDX_JHO2NO2) = get_rxt_ndx( 'jho2no2' )
      indexer(TAB_NDX_JCH3OOH) = get_rxt_ndx( 'jch3ooh' )
      indexer(TAB_NDX_JCH2Oa)  = get_rxt_ndx( 'jch2o_a' )
      indexer(TAB_NDX_JCH2Ob)  = get_rxt_ndx( 'jch2o_b' )
      indexer(TAB_NDX_JH2O2)   = get_rxt_ndx( 'jh2o2' )
      indexer(TAB_NDX_JCH3CHO) = get_rxt_ndx( 'jch3cho' )
      indexer(TAB_NDX_JPAN)    = get_rxt_ndx( 'jpan' )
      indexer(TAB_NDX_JMACRa)  = get_rxt_ndx( 'jmacr_a' )
      indexer(TAB_NDX_JMVK)    = get_rxt_ndx( 'jmvk' )
      indexer(TAB_NDX_JACET)   = get_rxt_ndx( 'jacet' )
      indexer(TAB_NDX_JMGLY)   = get_rxt_ndx( 'jmgly' )

      jno_ndx = get_rxt_ndx( 'jno' )
      jpooh_ndx = get_rxt_ndx( 'jpooh' )
      jc2h5ooh_ndx = get_rxt_ndx( 'jc2h5ooh' )
      jc3h7ooh_ndx = get_rxt_ndx( 'jc3h7ooh' )
      jrooh_ndx = get_rxt_ndx( 'jrooh' )
      jch3co3h_ndx = get_rxt_ndx( 'jch3co3h' )
      jmpan_ndx = get_rxt_ndx( 'jmpan' )
      jmacr_a_ndx = get_rxt_ndx( 'jmacr_a' )
      jmacr_b_ndx = get_rxt_ndx( 'jmacr_b' )
      jonitr_ndx = get_rxt_ndx( 'jonitr' )
      jxooh_ndx = get_rxt_ndx( 'jxooh' )
      jisopooh_ndx = get_rxt_ndx( 'jisopooh' )
      jglyald_ndx = get_rxt_ndx( 'jglyald' )
      jhyac_ndx = get_rxt_ndx( 'jhyac' )
      jch3ooh_ndx = get_rxt_ndx( 'jch3ooh' )
      jh2o2_ndx = get_rxt_ndx( 'jh2o2' )
      jpan_ndx = get_rxt_ndx( 'jpan' )
      jch3cho_ndx = get_rxt_ndx( 'jch3cho' )
      jn2o5_ndx = get_rxt_ndx( 'jn2o5' )
      jo3p_ndx = get_rxt_ndx( 'jo3p' )
      jno2_ndx = get_rxt_ndx( 'jno2' )
      jno3_ndx = get_rxt_ndx( 'jno3' )
      jclono2_ndx = get_rxt_ndx( 'jclono2' )
      jhocl_ndx = get_rxt_ndx( 'jhocl' )
      jcl2o2_ndx = get_rxt_ndx( 'jcl2o2' )
      jbrono2_ndx = get_rxt_ndx( 'jbrono2' )
      jhobr_ndx = get_rxt_ndx( 'jhobr' )
      jbrcl_ndx = get_rxt_ndx( 'jbrcl' )
      jbro_ndx = get_rxt_ndx( 'jbro' )
      jcl2_ndx = get_rxt_ndx( 'jcl2' )
      jh2o_ndx = get_rxt_ndx( 'jh2o' )

      ox_ndx = get_spc_ndx( 'OX' )
      o3_ndx = get_spc_ndx( 'O3' )

      end subroutine PRATE_INIT

      subroutine PHOTO( photos, pmid, pdel, temper, zmid, &
                        col_dens, &
                        coszen,  & 
                         srf_alb, lwc, clouds, &
                        esfact, &
                        plonl &
                        )

      use CHEM_MODS_MOD, only : ncol_abs, phtcnt
!      use M_RXT_ID_MOD

      implicit none

!-----------------------------------------------------------------
!           ... Dummy arguments
!-----------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in) ::   esfact                           ! earth sun distance factor
      real, intent(in) ::   col_dens(:,:,:), & ! column densities
                            coszen(:), &              ! solar zenith angle
                            srf_alb(:), &                ! surface albedo
                            lwc(:,:), &               ! liquid water content (mass mr)
                            clouds(:,:), &            ! cloud fraction
                            pmid(:,:), &              ! midpoint pressure in pascals
                            pdel(:,:), &              ! del pressure about midpoint in pascals
                            zmid(:,:), &              ! midpoint height
                            temper(:,:)               ! midpoint temperature
      real, intent(out) ::  photos(:,:,:)        ! photodissociation rates

!-----------------------------------------------------------------
!            ... Local variables
!-----------------------------------------------------------------
      integer  ::  i, k, m                 ! indicies
      integer  ::  plev
      logical  ::  zagtz(size(coszen))             ! zenith angle > 0 flag array
      real    ::   secant
      real    ::   t500, t200             ! 500 & 200 mb temperatures
      real, dimension(size(zmid,2)) :: &
                   fac1, &                ! work space for J(no) calc
                   fac2, &                ! work space for J(no) calc
                   colo3, &               ! vertical o3 column density
                   zarg, &                ! vertical height array
                   pline, &               ! vertical pressure array
                   tline, &               ! vertical temperature array
                   cld_line, &            ! vertical cloud array
                   lwc_line, &            ! vertical lwc array
                   eff_alb, &             ! effective albedo from cloud modifications
                   cld_mult               ! clould multiplier
      real, dimension(plonl,size(zmid,2)) :: &
                   tmp, &                        ! wrk array
                   tmp_jch3ooh, &                ! wrk array
                   tmp_jpan, &                   ! wrk array
                   tmp_jh2o2, &                  ! wrk array
                   tmp_jch3cho, &                ! wrk array
                   tmp_jmacr_a, &                ! wrk array
                   tmp_jn2o5, tmp_jo3p, tmp_jno2, tmp_jno3, &
                   tmp_jno
      real    ::   prates(jdim,size(zmid,2))        ! photorates matrix

      plev = SIZE(zmid,2)
!-----------------------------------------------------------------
!        ... Zero all photorates
!-----------------------------------------------------------------
      do m = 1,max(1,phtcnt)
         do k = 1,plev
            photos(:,k,m) = 0.
            tmp_jch3ooh(:,k) = 0.
            tmp_jpan(:,k)    = 0.
            tmp_jh2o2(:,k)   = 0.
            tmp_jch3cho(:,k) = 0.
            tmp_jmacr_a(:,k) = 0.
            tmp_jn2o5(:,k)   = 0.
            tmp_jo3p(:,k)    = 0.
            tmp_jno2(:,k)    = 0.
            tmp_jno3(:,k)    = 0.
            tmp_jno(:,k)     = 0.
         end do
      end do
      zagtz(:) = coszen(:) > 0. 

      do i = 1,plonl
         if( zagtz(i) ) then
            secant = 1. / coszen(i)
            if( secant <= 50. ) then
               zarg(:)     = zmid(i,:)
               colo3(:)    = col_dens(i,:,1)
               pline(:)    = pmid(i,:)
               fac1(:)     = pdel(i,:)
               tline(:)    = temper(i,:)
               lwc_line(:) = lwc(i,:)
               cld_line(:) = clouds(i,:)
               call cloud_mod( coszen(i), cld_line, lwc_line, fac1, srf_alb(i), &
                               eff_alb, cld_mult )
               call T_INT( pline, tline, t500, t200 )
               call PHOTO_INTERP( zarg, secant, colo3, eff_alb, t500, &
                                  t200, prates )
               do m = 1,jdim
                  if( indexer(m) > 0 ) then
                  photos(i,:,indexer(m)) = esfact *prates(m,:) * cld_mult(:)
                  else
                     select case( m )
                        case( TAB_NDX_JCH3OOH )
                           tmp_jch3ooh(i,:) = esfact *prates(m,:) * cld_mult(:)
                        case( TAB_NDX_JH2O2 )
                           tmp_jh2o2(i,:) = esfact *prates(m,:) * cld_mult(:)
                        case( TAB_NDX_JCH3CHO )
                           tmp_jch3cho(i,:) = esfact *prates(m,:) * cld_mult(:)
                        case( TAB_NDX_JPAN )
                           tmp_jpan(i,:) = esfact *prates(m,:) * cld_mult(:)
                        case( TAB_NDX_JMACRa )
                           tmp_jmacr_a(i,:) = esfact *prates(m,:) * cld_mult(:)
                        case( TAB_NDX_JN2O5 )
                           tmp_jn2o5(i,:) = esfact *prates(m,:) * cld_mult(:)
                        case( TAB_NDX_JO3P )
                           tmp_jo3p(i,:) = esfact *prates(m,:) * cld_mult(:)
                        case( TAB_NDX_JNO2 )
                           tmp_jno2(i,:) = esfact *prates(m,:) * cld_mult(:)
                        case( TAB_NDX_JNO3 )
                           tmp_jno3(i,:) = esfact *prates(m,:) * cld_mult(:)
                     end select
                  end if
               end do
               fac1(:) = 1.e-8  * (col_dens(i,:,2)/coszen(i))**.38
               fac2(:) = 5.e-19 * col_dens(i,:,1) / coszen(i)
               if( jno_ndx > 0 ) then
!-----------------------------------------------------------------
!        ... Calculate J(no) from formula
!-----------------------------------------------------------------
                  photos(i,:,jno_ndx) = 4.5e-6 * esfact * exp( -(fac1(:) + fac2(:)) ) * cld_mult(:)
               else
                  tmp_jno(i,:) = 4.5e-6 * esfact * exp( -(fac1(:) + fac2(:)) ) * cld_mult(:)
               end if
            end if
         end if
      end do
        
!-----------------------------------------------------------------
!        ... Set J(pooh) from J(ch3ooh)
!                J(c2h5ooh) from J(ch3ooh)
!                J(c3h7ooh) from J(ch3ooh)
!                J(rooh) from J(ch3ooh)
!                J(ch3coooh) = .28 * J(h2o2)
!                J(mpan) from J(pan)
!                J(macr_a) and J(macr_b) = 1/2 * J(macr_total)
!               J(onitr) from j(ch3cho)
!               J(xooh) from J(ch3ooh)
!               J(isopooh) from J(ch3ooh)
!                J(glyald) = 3 * J(ch3cho)
!               J(hyac) from J(ch3cho)
!-----------------------------------------------------------------
      if( jch3ooh_ndx > 0 ) then
         tmp(:,:) = photos(:,:,jch3ooh_ndx)
      else
         tmp(:,:) = tmp_jch3ooh(:,:)
      end if
      if( jpooh_ndx > 0 ) then
         photos(:,:,jpooh_ndx)    = tmp(:,:)
      end if
      if( jc2h5ooh_ndx > 0 ) then
         photos(:,:,jc2h5ooh_ndx) = tmp(:,:)
      end if
      if( jc3h7ooh_ndx > 0 ) then
         photos(:,:,jc3h7ooh_ndx) = tmp(:,:)
      end if
      if( jrooh_ndx > 0 ) then
         photos(:,:,jrooh_ndx)    = tmp(:,:)
      end if
      if( jxooh_ndx > 0 ) then
        photos(:,:,jxooh_ndx)    = tmp(:,:)
      end if
      if( jisopooh_ndx > 0 ) then
         photos(:,:,jisopooh_ndx) = tmp(:,:)
      end if
      if( jch3co3h_ndx > 0 ) then
         if( jh2o2_ndx > 0 ) then
            photos(:,:,jch3co3h_ndx) = .28 * photos(:,:,jh2o2_ndx)
         else
            photos(:,:,jch3co3h_ndx) = .28 * tmp_jh2o2(:,:)
         end if
      end if
      if( jmpan_ndx > 0 ) then
         if( jpan_ndx > 0 ) then
            photos(:,:,jmpan_ndx)    = photos(:,:,jpan_ndx)
         else
            photos(:,:,jmpan_ndx)    = tmp_jpan(:,:)
         end if
      end if
      if( jmacr_a_ndx > 0 ) then
         photos(:,:,jmacr_a_ndx)  = .5 * photos(:,:,jmacr_a_ndx)
      end if
      if( jmacr_b_ndx > 0 ) then
         if( jmacr_a_ndx > 0 ) then
            photos(:,:,jmacr_b_ndx)  = photos(:,:,jmacr_a_ndx)
         else
            photos(:,:,jmacr_b_ndx)  = .5 * tmp_jmacr_a(:,:)
         end if
      end if
      if( jonitr_ndx > 0 ) then
         if( jch3cho_ndx > 0 ) then
            photos(:,:,jonitr_ndx)   = photos(:,:,jch3cho_ndx)
         else
            photos(:,:,jonitr_ndx)   = tmp_jch3cho(:,:)
         end if
      end if
      if( jglyald_ndx > 0 ) then
         if( jch3cho_ndx > 0 ) then
            photos(:,:,jglyald_ndx)  = 3. * photos(:,:,jch3cho_ndx)
         else
            photos(:,:,jglyald_ndx)   = 3. *tmp_jch3cho(:,:)
         end if
      end if
      if( jhyac_ndx > 0 ) then
         if( jch3cho_ndx > 0 ) then
            photos(:,:,jhyac_ndx)    = photos(:,:,jch3cho_ndx)
         else
            photos(:,:,jhyac_ndx)    = tmp_jch3cho(:,:)
         end if
      end if
      if( jclono2_ndx > 0 ) then
         if( jn2o5_ndx > 0 ) then
            photos(:,:,jclono2_ndx) = photos(:,:,jn2o5_ndx)
         else
            photos(:,:,jclono2_ndx) = tmp_jn2o5(:,:)
         end if
      end if
      if( jhocl_ndx > 0 ) then
         if( jo3p_ndx > 0 ) then
            photos(:,:,jhocl_ndx) = photos(:,:,jo3p_ndx)
         else
            photos(:,:,jhocl_ndx) = tmp_jo3p(:,:)
         end if
      end if
      if( jcl2o2_ndx > 0 ) then
         if( jno2_ndx > 0 ) then
            photos(:,:,jcl2o2_ndx) = photos(:,:,jno2_ndx)
         else
            photos(:,:,jcl2o2_ndx) = tmp_jno2(:,:)
         end if
      end if
      if( jbrono2_ndx > 0 ) then
         if( jo3p_ndx > 0 ) then
            photos(:,:,jbrono2_ndx) = 2*photos(:,:,jo3p_ndx)
         else
            photos(:,:,jbrono2_ndx) = 2*tmp_jo3p(:,:)
         end if
      end if
      if( jhobr_ndx > 0 ) then
         if( jo3p_ndx > 0 ) then
            photos(:,:,jhobr_ndx) = photos(:,:,jo3p_ndx)
         else
            photos(:,:,jhobr_ndx) = tmp_jo3p(:,:)
         end if
      end if
      if( jbrcl_ndx > 0 ) then
         if( jno2_ndx > 0 ) then
            photos(:,:,jbrcl_ndx) = photos(:,:,jno2_ndx)
         else
            photos(:,:,jbrcl_ndx) = tmp_jno2(:,:)
         end if
      end if
      if( jbro_ndx > 0 ) then
         if( jno3_ndx > 0 ) then
            photos(:,:,jbro_ndx) = photos(:,:,jno3_ndx)
         else
            photos(:,:,jbro_ndx) = tmp_jno3(:,:)
         end if
      end if
      if( jcl2_ndx > 0 ) then
         if( jno2_ndx > 0 ) then
            photos(:,:,jcl2_ndx) = photos(:,:,jno2_ndx)
         else
            photos(:,:,jcl2_ndx) = tmp_jno2(:,:)
         end if
      end if
      if( jh2o_ndx > 0 ) then
         if( jno_ndx > 0 ) then
            photos(:,:,jh2o_ndx) = 0.1*photos(:,:,jno_ndx)
         else
            photos(:,:,jh2o_ndx) = 0.1*tmp_jno(:,:)
         end if
      end if

      end subroutine PHOTO

      subroutine cloud_mod( coszen, clouds, lwc, delp, srf_alb, &
                            eff_alb, cld_mult )
!-----------------------------------------------------------------------
!         ... Cloud alteration factors for photorates and albedo
!-----------------------------------------------------------------------


      implicit none

      real, parameter :: gi = 1./9.80616

!-----------------------------------------------------------------------
!         ... Dummy arguments
!-----------------------------------------------------------------------
      real, intent(in)    ::  coszen             ! cosine of zenith angle
      real, intent(in)    ::  srf_alb            ! surface albedo
      real, intent(in)    ::  clouds(:)          ! cloud fraction
      real, intent(in)    ::  lwc(:)             ! liquid water content (mass mr)
      real, intent(in)    ::  delp(:)            ! del press about midpoint in pascals
      real, intent(out)   ::  eff_alb(:)         ! effective albedo
      real, intent(out)   ::  cld_mult(:)        ! photolysis mult factor

!-----------------------------------------------------------------------
!         ... Local variables
!-----------------------------------------------------------------------
      integer :: k
      integer :: plev, plevm
      real    :: coschi
      real    :: del_lwp(SIZE(clouds,1))
      real    :: del_tau(SIZE(clouds,1))
      real    :: above_tau(SIZE(clouds,1))
      real    :: below_tau(SIZE(clouds,1))
      real    :: above_cld(SIZE(clouds,1))
      real    :: below_cld(SIZE(clouds,1))
      real    :: above_tra(SIZE(clouds,1))
      real    :: below_tra(SIZE(clouds,1))
      real    :: fac1(SIZE(clouds,1))
      real    :: fac2(SIZE(clouds,1))
      
      plev = SIZE(clouds,1)
      plevm = plev-1
!---------------------------------------------------------
!        ... Modify lwc for cloud fraction and form
!            liquid water path for each layer
!---------------------------------------------------------
      where( clouds(:) /= 0. )
         del_lwp(:) = gi * lwc(:) * delp(:) * 1.e3 / clouds(:)
      elsewhere
         del_lwp(:) = 0.
      endwhere
!---------------------------------------------------------
!            ... Form tau for each model layer
!---------------------------------------------------------
      where( clouds(:) /= 0. )
         del_tau(:) = del_lwp(:) *.155 * clouds(:)**1.5
      elsewhere
         del_tau(:) = 0.
      end where
!---------------------------------------------------------
!            ... Form integrated tau from top down
!---------------------------------------------------------
      above_tau(1) = 0.
      do k = 1,plevm
         above_tau(k+1) = del_tau(k) + above_tau(k)
      end do
!---------------------------------------------------------
!            ... Form integrated tau from bottom up
!---------------------------------------------------------
      below_tau(plev) = 0.
      do k = plevm,1,-1
         below_tau(k) = del_tau(k+1) + below_tau(k+1)
      end do
!---------------------------------------------------------
!        ... Form vertically averaged cloud cover above and below
!---------------------------------------------------------
      above_cld(1) = 0.
      do k = 1,plevm
         above_cld(k+1) = clouds(k) * del_tau(k) + above_cld(k)
      end do
      do k = 2,plev
         if( above_tau(k) /= 0. ) then
            above_cld(k) = above_cld(k) / above_tau(k)
         else
            above_cld(k) = above_cld(k-1)
         end if
      end do
      below_cld(plev) = 0.
      do k = plevm,1,-1
         below_cld(k) = clouds(k+1) * del_tau(k+1) + below_cld(k+1)
      end do
      do k = plevm,1,-1
         if( below_tau(k) /= 0. ) then
            below_cld(k) = below_cld(k) / below_tau(k)
         else
            below_cld(k) = below_cld(k+1)
         end if
      end do
!---------------------------------------------------------
!        ... Modify above_tau and below_tau via jfm
!---------------------------------------------------------
      where( above_cld(2:plev) /= 0. )
         above_tau(2:plev) = above_tau(2:plev) / above_cld(2:plev)
      end where
      where( below_cld(:plevm) /= 0. )
         below_tau(:plevm) = below_tau(:plevm) / below_cld(:plevm)
      end where
      where( above_tau(2:plev) < 5. )
            above_cld(2:plev) = 0.
      end where
      where( below_tau(:plevm) < 5. )
         below_cld(:plevm) = 0.
      end where
!---------------------------------------------------------
!        ... Form transmission factors
!---------------------------------------------------------
      above_tra(:) = 11.905 / (9.524 + above_tau(:))
      below_tra(:) = 11.905 / (9.524 + below_tau(:))
!---------------------------------------------------------
!        ... Form effective albedo
!---------------------------------------------------------
      where( below_cld(:) /= 0. )
         eff_alb(:) = srf_alb + below_cld(:) * (1. - below_tra(:)) &
                                             * (1. - srf_alb)
      elsewhere
         eff_alb(:) = srf_alb
      end where
      coschi = max( coszen,.5 )
      where( del_lwp(:)*.155 < 5. )
         fac1(:) = 0.
      elsewhere
         fac1(:) = 1.4 * coschi - 1.
      end where
      fac2(:)     = MIN( 0.,1.6*coschi*above_tra(:) - 1. )
      cld_mult(:) = 1. + fac1(:) * clouds(:) + fac2(:) * above_cld(:)
      cld_mult(:) = MAX( .05,cld_mult(:) )

      end subroutine cloud_mod

      subroutine T_INT( p, t, t500, t200 )
!----------------------------------------------------------------
!        ... Interpolate for temperature on 500 and 200mb surfaces
!----------------------------------------------------------------


      implicit none

!----------------------------------------------------------------
!        ... Dummy args
!----------------------------------------------------------------
      real, intent(in)  ::  p(:)              ! pressure in pascals
      real, intent(in)  ::  t(:)              ! temperature on grid
      real, intent(out) ::  t500, t200        ! temp at 500 and 200mb

!----------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------
      integer :: k, k1
      real    :: delp
      integer :: plev, plevp, plevm
      
      plev = SIZE(p)
      plevp = plev+1
      plevm = plev-1

      if( p(plev) < 500.e2 ) then
         t500 = t(plev)
         k1 = plevp
      else
         do k = plevm,1,-1
            if( p(k) < 500.e2 ) then
               k1 = k
               exit
            end if
         end do
         delp = LOG( 500.e2/p(k) ) / LOG( p(k+1)/p(k) )
         t500 = t(k) + delp * (t(k+1) - t(k))
      end if
      do k = k1-1,1,-1
         if( p(k) < 200.e2 ) then
            exit
         end if
      end do
      delp = LOG( 200.e2/p(k) ) / LOG( p(k+1)/p(k) )
      t200 = t(k) + delp * (t(k+1) - t(k))

      end subroutine T_INT

      subroutine PHOTO_INTERP( zin, sin, vin, albin, t500in, &
                               t200in, ajout )
!----------------------------------------------------------------------
!           ... Loglinear interpolation for the photodissociation rates
!            Note: this subroutine computes photorates for a vertical
!                  column at a given longitude and latitude
!           This routine uses a six parameter table via a Taylor
!           series expansion. 
!           Stacy Walters, Sep 30, 1996.  Changed code to strictly limit
!           the 200mb and 500mb temperature interpolation to the table
!           endpoints; i.e. no extrapolation beyond the table is allowed.
!----------------------------------------------------------------------


      implicit none

!----------------------------------------------------------------------
!        ... Dummy arguments
!----------------------------------------------------------------------
      real, intent(in)  ::   zin(:), &              ! geo height of midpoints
                             sin, &                 ! secant solar zenith angle
                             vin(:), &              ! o3 column density
                             albin(:), &            ! surface albedo
                             t500in, &              ! temp on 500mb surface
                             t200in                 ! temp on 200mb surface
      real, intent(out) ::   ajout(:,:)             ! photodissociation rates

!----------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------
      integer, parameter ::  it200 = 1, it200p1 = 2
      integer  ::  plev
      integer  ::  iz, is, iv, ial, nn, it500
      integer  ::  izp1, isp1, ivp1, ialp1, it500p1
      integer  ::  i, k
      integer  ::  izl
      integer, dimension(SIZE(zin)) :: altind, ratind, albind
      real     ::  wght0
      real     ::  v3std
      real     ::  dels(6)
      real, dimension(SIZE(zin)) :: v3rat
      
      plev = SIZE(zin)

!----------------------------------------------------------------------
!        ... Find the zenith angle index ( same for all levels )
!----------------------------------------------------------------------
      do is = 1,zangdim
         if( vsec(is) > sin ) then
            exit
         end if
      end do
      is       = MAX( MIN( is,zangdim ) - 1,1 )
      isp1     = is + 1
      dels(2)  = MIN( 1.,MAX( 0.,(sin - vsec(is)) * delang(is) ) )

!----------------------------------------------------------------------
!        ... Find the 500mb temp index ( same for all levels )
!----------------------------------------------------------------------
      do it500 = 1,t500dim
         if( t500(it500) > t500in ) then
            exit
         end if
      end do
      it500    = MAX( MIN( it500,t500dim ) - 1,1 )
      it500p1  = it500 + 1
      dels(5)  = MIN( 1.,MAX( 0.,(t500in - t500(it500)) * delt500(it500) ) )

!----------------------------------------------------------------------
!        ... Find the 200mb temp index ( same for all levels )
!----------------------------------------------------------------------
      dels(6)  = MIN( 1.,MAX( 0.,(t200in - t200(it200)) * delt200(it200) )) 

      izl = 1
      do k = plev,1,-1
!----------------------------------------------------------------------
!        ... Find albedo indicies
!----------------------------------------------------------------------
         do ial = 1,albdim
            if( albev(ial) > albin(k) ) then
               exit
            end if
         end do
         albind(k) = MAX( MIN( ial,albdim ) - 1,1 )
!----------------------------------------------------------------------
!        ... Find level indicies
!----------------------------------------------------------------------
         do iz = izl,altdim
            if( zz(iz) > zin(k) ) then
               izl = iz
               exit
            end if
         end do
         altind(k) = MAX( MIN( iz,altdim ) - 1,1 )
!----------------------------------------------------------------------
!        ... Find "o3 ratio" indicies
!----------------------------------------------------------------------
         i        = MAX( MIN( 49,INT( zin(k) ) ),0 )
         v3std    = vo3(i) + (zin(k) - REAL(i)) * delvo3(i)
         v3rat(k) = vin(k) / v3std
         do iv = 1,o3ratdim
            if( xv3(iv) > v3rat(k) ) then
               exit
            end if
         end do
         ratind(k) = MAX( MIN( iv,o3ratdim ) - 1,1 )
      end do
Vert_loop : &
      do k = 1,plev
         iz    = altind(k)
         izp1  = iz + 1
         iv    = ratind(k)
         ivp1  = iv + 1
         ial   = albind(k)
         ialp1 = ial + 1
!----------------------------------------------------------------------
!        ... Interval deltas and primary weight
!----------------------------------------------------------------------
         dels(1)  = MIN( 1.,MAX( 0.,(zin(k) - zz(iz)) * delz(iz) ) )
         dels(3)  = MIN( 1.,MAX( 0.,(v3rat(k) - xv3(iv)) * delv(iv) ) )
         dels(4)  = MIN( 1.,MAX( 0.,(albin(k) - albev(ial)) * delalb(ial) ) )
         wght0    = 1. - SUM( dels )
Rate_loop : &
         do nn = 1,jdim
            ajout(nn,k) = EXP( wght0     * ajl(nn,iz,is,iv,ial,it500,it200) &
                               + dels(1) * ajl(nn,izp1,is,iv,ial,it500,it200) &
                               + dels(2) * ajl(nn,iz,isp1,iv,ial,it500,it200) &
                               + dels(3) * ajl(nn,iz,is,ivp1,ial,it500,it200) &
                               + dels(4) * ajl(nn,iz,is,iv,ialp1,it500,it200) &
                               + dels(5) * ajl(nn,iz,is,iv,ial,it500p1,it200) &
                               + dels(6) * ajl(nn,iz,is,iv,ial,it500,it200p1) )
         end do Rate_loop
      end do Vert_loop

      end subroutine PHOTO_INTERP

      subroutine set_ub_col( col_delta, vmr, invariants, pdel, plonl )
!---------------------------------------------------------------
!        ... Set the column densities at the upper boundary
!---------------------------------------------------------------

      use CHEM_MODS_MOD, only : ncol_abs

      implicit none

!---------------------------------------------------------------
!        ... Dummy args
!---------------------------------------------------------------
      integer, intent(in) ::  plonl
      real, intent(out)   ::  col_delta(:,0:,:)  ! /cm**2
      real, intent(in)    ::  vmr(:,:,:), &               ! xported species vmr
                              invariants(:,:,:), &        ! invariant species
                              pdel(:,:)

!---------------------------------------------------------------
!        NOTE: xfactor = 10.*R/(K*g) in cgs units.
!              The factor 10. is to convert pdel
!              from pascals to dyne/cm**2.
!---------------------------------------------------------------
      real, parameter :: pa_to_dyncm2 = 10. !unit = (dyn/cm2)/pa
      real, parameter :: mw_air = 28.9644   !g/mole
      real, parameter :: grav = 981.         !cm/s2
      real, parameter :: navo = 6.023e23   ! molec/mole
      real, parameter :: xfactor = pa_to_dyncm2 * navo /(grav * mw_air)
      integer :: k, spc_ndx
      integer :: plev
      
      plev = SIZE(invariants,2)

!---------------------------------------------------------------
!        ... Assign column density at the upper boundary
!            The first column is O3 and the second is O2.
!            Add 10 DU O3 column above top of model.
!---------------------------------------------------------------
      spc_ndx = ox_ndx
      if( spc_ndx < 1 ) then
         spc_ndx = o3_ndx
      end if
      if( spc_ndx > 0 ) then
!        col_delta(:,0,1) = 2.687e16*10.
         col_delta(:,0,1) = 2.687e16*0.1
         do k = 1,plev
            col_delta(:,k,1) = xfactor * pdel(:,k) * vmr(:,k,spc_ndx) ! O3
         end do
      end if
      col_delta(:,0,2) = 2.8e22
      do k = 1,plev
         col_delta(:,k,2) = xfactor * pdel(:,k) * invariants(:,k,3)/invariants(:,k,1) ! O2
      end do

      end subroutine set_ub_col

      subroutine setcol( col_delta, col_dens, pdel, plonl )
!---------------------------------------------------------------
!             ... Set the column densities
!---------------------------------------------------------------

      use CHEM_MODS_MOD, only : ncol_abs

      implicit none

!---------------------------------------------------------------
!             ... Dummy arguments
!---------------------------------------------------------------
      integer, intent(in) :: plonl
!     real, intent(in)  ::   vmr(plonl,plev,pcnstm1)           ! xported species vmr
      real, intent(in)  ::   pdel(:,:)                         ! delta about midpoints
      real, intent(in)  ::   col_delta(:,0:,:)                 ! layer column densities (molecules/cm^2)
      real, intent(out) ::   col_dens(:,:,:)                   ! column densities ( /cm**2 )

!---------------------------------------------------------------
!        The local variables
!---------------------------------------------------------------
      integer  ::   i, k, km1      ! long, alt indicies
      integer  ::   spc_ndx
      integer  ::   plev
      
!---------------------------------------------------------------
!        NOTE: xfactor = 10.*R/(K*g) in cgs units.
!              The factor 10. is to convert pdel
!              from pascals to dyne/cm**2.
!---------------------------------------------------------------
!     real, parameter :: xfactor = 2.8704e21/(9.80616*1.38044)

      plev = SIZE(pdel,2)

!---------------------------------------------------------------
!           ... Compute column densities down to the
!           current eta index in the calling routine.
!           The first column is O3 and the second is O2.
!---------------------------------------------------------------
      spc_ndx = ox_ndx
      if( spc_ndx < 1 ) then
         spc_ndx = o3_ndx
      end if
      if( spc_ndx > 0 ) then
         col_dens(:,1,1) = col_delta(:,0,1) + .5 * col_delta(:,1,1)
         do k = 2,plev
            km1 = k - 1
            col_dens(:,k,1) = col_dens(:,km1,1) + .5 * (col_delta(:,km1,1) + col_delta(:,k,1))
         end do
      end if
      col_dens(:,1,2) = col_delta(:,0,2) + .5 * col_delta(:,1,2)
      do k = 2,plev
         km1 = k - 1
         col_dens(:,k,2) = col_dens(:,km1,2) + .5 * (col_delta(:,km1,2) + col_delta(:,k,2))
      end do

      end subroutine SETCOL

!     subroutine diurnal_geom( ip, lat, time_of_year, polar_night, polar_day, &
!                              sunon, sunoff, loc_angle, zen_angle, plonl )
!------------------------------------------------------------------
!            ... Diurnal geometry factors
!------------------------------------------------------------------


!     implicit none

!------------------------------------------------------------------
!            ... Dummy arguments
!------------------------------------------------------------------
!     integer, intent(in)  ::     ip                 ! longitude index
!     integer, intent(in)  ::     lat                ! latitude index
!     integer, intent(in)  ::     plonl
!     real, intent(in)     ::     time_of_year       ! time of year
!     real, intent(out)    ::     sunon           ! sunrise angle in radians
!     real, intent(out)    ::     sunoff          ! sunset angle in radians
!     real, intent(out)    ::     zen_angle(plonl) ! solar zenith angle
!     real, intent(out)    ::     loc_angle(plonl) ! "local" time angle
!     logical, intent(out) ::     polar_day       ! continuous daylight flag
!     logical, intent(out) ::     polar_night     ! continuous night flag

!------------------------------------------------------------------
!        ... Local variables
!------------------------------------------------------------------
!     integer ::  i
!     real    ::  dec_max
!     real    ::  declination
!     real    ::  latitude
!     real    ::  doy_loc            ! day of year
!     real    ::  tod                ! time of day
!     real    ::  sin_dec, cos_dec   ! sin, cos declination
!     real    ::  cosphi             ! cos latitude
!     real    ::  sinphi             ! sin latitude

!     dec_max     = 23.45 * d2r
!     latitude    = phi(base_lat + lat)
!     sinphi      = sin( latitude )
!     cosphi      = cos( latitude )
!     polar_day   = .false.
!     polar_night = .false.
!------------------------------------------------------------------
!        Note: this formula assumes a 365 day year !
!------------------------------------------------------------------
!     doy_loc     = aint( time_of_year )
!     declination = dec_max * cos((doy_loc - 172.)*twopi/dayspy)
!------------------------------------------------------------------
!        Determine if in polar day or night
!        If NOT in polar day or night then
!        calculate terminator longitudes
!------------------------------------------------------------------
!     if( abs(latitude) >= (pid2 - abs(declination)) ) then
!         if( sign(1.,declination) == sign(1.,latitude) ) then
!            polar_day = .true.
!            sunoff    = 2.*twopi
!            sunon     = -twopi
!        else
!            polar_night  = .true.
!           zen_angle(:) = -1.0
!            return
!        end if
!     else
!        sunoff = acos( -tan(declination)*tan(latitude) )
!        sunon  = twopi - sunoff
!     end if

!     sin_dec = sin( declination )
!     cos_dec = cos( declination )
!------------------------------------------------------------------
!        ... Compute base for zenith angle
!------------------------------------------------------------------
!     tod = (time_of_year - doy_loc) + .5
!-------------------------------------------------------------------
!        Note: Longitude 0 (Greenwich) at 0:00 hrs
!              maps to local angle = pi
!-------------------------------------------------------------------
!     loc_angle(:) = (/ ((tod + real(i+(ip-1)*plonl-1)/real(plong))*twopi,i = 1,plonl) /)
!     loc_angle(:) = mod( loc_angle(:),twopi )

!     if( polar_day ) then
!         zen_angle(:) = acos( sinphi*sin_dec + cosphi*cos_dec*cos(loc_angle(:)) )
!     else
!         where( loc_angle(:) <= sunoff .or. loc_angle(:) >= sunon )
!            zen_angle(:) = acos( sinphi*sin_dec + cosphi*cos_dec*cos(loc_angle(:)) )
!         elsewhere
!           zen_angle(:) = -1.
!         endwhere
!     end if

!     end subroutine diurnal_geom

      real function SUNDIS( Time )
!-----------------------------------------------------------------------------
!=  PURPOSE:                                                                 =*
!=  Calculate Earth-Sun distance variation for a given date.  Based on       =*
!=  Fourier coefficients originally from:  Spencer, J.W., 1971, Fourier      =*
!=  series representation of the position of the sun, Search, 2:172          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  IDATE  - INTEGER, specification of the date, from YYMMDD              (I)=*
!=  ESRM2  - REAL, variation of the Earth-sun distance                    (O)=*
!=           ESRM2 = (average e/s dist)^2 / (e/s dist on day IDATE)^2        =*
!-----------------------------------------------------------------------------*
!=  EDIT HISTORY:                                                            =*
!=  01/95  Changed computation of trig function values                       =*
!-----------------------------------------------------------------------------*
!= This program is free software;  you can redistribute it and/or modify     =*
!= it under the terms of the GNU General Public License as published by the  =*
!= Free Software Foundation;  either version 2 of the license, or (at your   =*
!= option) any later version.                                                =*
!= The TUV package is distributed in the hope that it will be useful, but    =*
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!= License for more details.                                                 =*
!= To obtain a copy of the GNU General Public License, write to:             =*
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!-----------------------------------------------------------------------------*
!= To contact the authors, please mail to:                                   =*
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!= send email to:  sasha@ucar.edu                                            =*
!-----------------------------------------------------------------------------*
!= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!-----------------------------------------------------------------------------


      implicit none

!-----------------------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------------------
      type(time_type), intent(in) :: Time             ! time

!-----------------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------------
      integer :: iyear, imonth, iday, ihour, iminute, isecond
      integer :: mday, month, jday
      integer, save :: imn(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
      real    :: dayn, thet0
      real    :: sinth, costh, sin2th, cos2th
      character(len=128) :: msg
      
!-----------------------------------------------------------------------------
!         ... Parse date to find day number (Julian day)
!-----------------------------------------------------------------------------
      call get_date( Time, iyear, imonth, iday, ihour, iminute, isecond )
      if( imonth > 12 ) then
         write(msg,*) 'Month in date exceeds 12, month = ',imonth
         call endrun(msg)
      end if

      if( MOD(iyear,4) == 0 ) then
         imn(2) = 29
      else
         imn(2) = 28
      end if

      if( iday > imn(imonth) ) then
         write(msg,*) 'Day in date exceeds days in month, day = ',iday,', month = ',imonth
         call endrun(msg)
      end if

      mday = 0
      do month = 1,imonth-1
         mday = mday + imn(month)                     
      end do
      jday = mday + iday
      dayn = REAL(jday - 1) + .5

!-----------------------------------------------------------------------------
!         ... Define angular day number and compute esrm2:
!-----------------------------------------------------------------------------
      thet0 = 2.*PI*dayn/365.

!-----------------------------------------------------------------------------
!         ... Calculate SIN(2*thet0), COS(2*thet0) 
!-----------------------------------------------------------------------------
      sinth   = SIN( thet0 )
      costh   = COS( thet0 )
      sin2th  = 2.*sinth*costh
      cos2th  = costh*costh - sinth*sinth
      SUNDIS  = 1.000110 + .034221*costh  +  .001280*sinth + .000719*cos2th +  .000077*sin2th

      end function SUNDIS

      subroutine endrun(msg)

      character(len=128), intent(in) :: msg
      call mpp_error(FATAL, msg)
      
      end subroutine endrun        

      end module MO_PHOTO_MOD
