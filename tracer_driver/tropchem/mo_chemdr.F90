
      module mo_chemdr_mod

      implicit none

      private
      public :: chemdr

!     save

      real :: esfact = 1.           ! earth sun distance factor

character(len=128), parameter :: version     = '$Id: mo_chemdr.F90,v 13.0 2006/03/28 21:16:03 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: memphis_2006_07 $'
logical                       :: module_is_initialized = .false.

      contains

      subroutine chemdr( vmr, &
                         Time, &
                         lat, lon, &
                         delt, &
                         ps, pmid, pdel, &
                         zma, zi, &
                         cldfr, cwat, tfld, inv_data, sh, &
                         albedo, coszen, prod_out, loss_out, sulfate, &
                         plonl )
!-----------------------------------------------------------------------
!     ... Chem_solver advances the volumetric mixing ratio
!         forward one time step via a combination of explicit,
!         ebi, hov, fully implicit, and/or rodas algorithms.
!-----------------------------------------------------------------------

      use chem_mods_mod,    only : indexm, nadv_mass, phtcnt, gascnt, rxntot, clscnt1, clscnt4, clscnt5, &
                                   ncol_abs, grpcnt, nfs, extcnt, hetcnt
      use mo_photo_mod,     only : set_ub_col, setcol, photo, &
                                   sundis
      use mo_exp_sol_mod,   only : exp_sol
      use mo_imp_sol_mod,   only : imp_sol
      use mo_rodas_sol_mod, only : rodas_sol
      use mo_usrrxt_mod,    only : usrrxt
      use mo_setinv_mod,    only : setinv
      use mo_setrxt_mod,    only : setrxt
      use mo_adjrxt_mod,    only : adjrxt
      use mo_phtadj_mod,    only : phtadj
      use mo_setsox_mod,    only : setsox
      use mo_chem_utls_mod, only : inti_mr_xform, adjh2o, negtrc, mmr2vmr, vmr2mmr, &
                                   get_spc_ndx, get_grp_mem_ndx
      use time_manager_mod, only : time_type

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------------
!     integer, intent(in) ::  nstep                   ! time index
      type(time_type), intent(in) :: Time             ! time
      real,    intent(in) ::  lat(:), lon(:)          ! latitude, longitude
      integer, intent(in) ::  plonl
      real,    intent(in) ::  delt                    ! timestep in seconds
      real, intent(inout) ::  vmr(:,:,:)              ! transported species ( vmr )
      real, dimension(:), intent(in) :: &
                              ps, &                   ! surface press ( pascals )
!                             oro, &                  ! surface orography flag
                              albedo, &               ! surface albedo
                              coszen                  ! cosine of solar zenith angle
!                             tsurf, &                ! surface temperature
!                             phis, &                 ! surf geopot
!                             cldtop                  ! cloud top level ( 1 ... plev )
      real, dimension(:,:), intent(in) :: &
                              pmid, &                 ! midpoint press ( pascals )
                              pdel, &                 ! delta press across midpoints
                              zma, &                  ! abs geopot height at midpoints ( m )
                              cldfr, &                ! cloud fraction
!                             cmfdqr, &               ! dq/dt for convective rainout
!                             nrain, &                ! release of strt precip ( 1/s )
!                             nevapr, &               ! evap precip ( 1/s )
                              cwat, &                 ! total cloud water (kg/kg)
                              tfld, &                 ! midpoint temperature
                              sh, &                   ! specific humidity ( kg/kg )
                              sulfate                 ! sulfate aerosol
      real, dimension(:,:), intent(in) :: &
                              zi                      ! abs geopot height at interfaces ( m )
      real, dimension(:,:,:), intent(out) :: &
                              prod_out, &             ! chemical production rate
                              loss_out                ! chemical loss rate
      real, dimension(:,:,:), intent(in) :: &
                              inv_data                ! invariant species

!-----------------------------------------------------------------------
!             ... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: inst = 1, avrg = 2
      integer  :: idate
      integer  ::  i, k, m, n, hndx, file
      integer  ::  ox_ndx, o3_ndx
      integer  ::  so2_ndx, so4_ndx
!     real     ::  sunon, &
!                  sunoff, &
!                  caldayn               ! day of year at end of time step
      real     ::  invariants(plonl,SIZE(vmr,2),max(1,nfs))
      real     ::  group_ratios(plonl,SIZE(vmr,2),max(1,grpcnt))
      real     ::  group_vmr(plonl,SIZE(vmr,2),max(1,grpcnt))
      real     ::  col_dens(plonl,SIZE(vmr,2),max(1,ncol_abs))                  ! column densities (molecules/cm^2)
      real     ::  col_delta(plonl,0:SIZE(vmr,2),max(1,ncol_abs))               ! layer column densities (molecules/cm^2)
      real     ::  het_rates(plonl,SIZE(vmr,2),max(1,hetcnt))
      real     ::  extfrc(plonl,SIZE(vmr,2),max(1,extcnt))
!     real     ::  vmr(plonl,SIZE(vmr,2),pcnstm1)                               ! xported species ( vmr )
      real     ::  reaction_rates(plonl,SIZE(vmr,2),rxntot)
      real     ::  nas(plonl,SIZE(vmr,2),max(1,grpcnt))                         ! non-advected species( mmr )
      real, dimension(plonl,SIZE(vmr,2)) :: &
                   h2ovmr, &             ! water vapor volume mixing ratio
                   mbar, &               ! mean wet atmospheric mass ( amu )
                   zmid                  ! midpoint geopotential in km
!                  sulfate               ! sulfate aerosols
      real, dimension(plonl,SIZE(zi,2)) :: &
                   zint                  ! interface geopotential in km
      real, dimension(plonl)  :: &
                   zen_angle, &
                   loc_angle, &
                   albs
      logical  ::  polar_night, &
                   polar_day
!     logical  ::  group_write(moz_file_cnt)
      character(len=32) :: fldname
      integer :: plev, plevp, plnplv, num_invar
      integer :: nstep

      plev = SIZE(vmr,2)
      plevp = SIZE(zi,2)
      plnplv = plonl*plev
      num_invar = SIZE(invariants,3)
      nstep = 0
      
      
!     caldayn = caldayr( ncdate, ncsec )
      if( phtcnt /= 0 ) then
!-----------------------------------------------------------------------      
!        ... Calculate parameters for diurnal geometry
!-----------------------------------------------------------------------      
!        call diurnal_geom( ip, lat, caldayn, polar_night, polar_day, &
!                           sunon, sunoff, loc_angle,  zen_angle, plonl )
      end if
!-----------------------------------------------------------------------      
!        ... Initialize xform between mass and volume mixing ratios
!-----------------------------------------------------------------------      
      call inti_mr_xform( sh, mbar, plonl )
!-----------------------------------------------------------------------      
!        ... Xform from mmr to vmr
!-----------------------------------------------------------------------      
!     call mmr2vmr( vmr, mmr, mbar, plonl )
!-----------------------------------------------------------------------      
!        ... Xform water vapor from mmr to vmr and adjust in stratosphere
!-----------------------------------------------------------------------      
      call adjh2o( h2ovmr, sh, mbar, vmr, plonl )
!-----------------------------------------------------------------------      
!        ... Xform geopotential height from m to km 
!            and pressure from hPa to mb
!-----------------------------------------------------------------------      
      do k = 1,plev
         zmid(:,k) = 1.e-3 * zma(:,k)
         zint(:,k) = 1.e-3 * zi(:,k)
      end do
      zint(:,plevp) = 1.e-3 * zi(:,plevp)

      if( nfs > 0 ) then
!-----------------------------------------------------------------------      
!        ... Set the "invariants"
!-----------------------------------------------------------------------      
         call setinv( invariants, tfld, h2ovmr, pmid, inv_data, plonl )
      end if
      if( ncol_abs > 0 .and. phtcnt > 0 ) then
!-----------------------------------------------------------------------      
!        ... Xform family ox assuming that all ox is o3
!-----------------------------------------------------------------------      
         ox_ndx = get_spc_ndx( 'OX' )
         if( ox_ndx > 0 ) then
            o3_ndx = get_grp_mem_ndx( 'O3' )
            if( o3_ndx > 0 ) then
!              vmr(:,:,ox_ndx) = mbar(:,:) * mmr(:,:,ox_ndx) / nadv_mass(o3_ndx)
            end if
         end if
!-----------------------------------------------------------------------      
!        ... Set the column densities at the upper boundary
!-----------------------------------------------------------------------      
         call set_ub_col( col_delta, vmr, invariants, pdel, plonl )
      end if
      if( gascnt > 0 ) then
!-----------------------------------------------------------------------      
!       ...  Set rates for "tabular" and user specified reactions
!-----------------------------------------------------------------------      
         call setrxt( reaction_rates, tfld, invariants(:,:,indexm), plonl, plev, plnplv )
!        call sulf_interp( lat, ip, pmid, caldayn, sulfate, plonl )
         call usrrxt( reaction_rates, tfld, invariants, h2ovmr, &
                      pmid, invariants(:,:,indexm), sulfate, vmr, sh, plonl )
!-----------------------------------------------------------------------      
!       ...  History output for instantaneous reaction rates
!-----------------------------------------------------------------------      
!        do file = 1,moz_file_cnt
!           if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(9,inst) > 0 ) then
!              do m = 1,hfile(file)%histout_cnt(9,inst)
!                  fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(9,inst)+m-1)
!                  hndx    = hfile(file)%inst_map(hfile(file)%histout_ind(9,inst)+m-1)
!                 call outfld( fldname, reaction_rates(1,1,hndx+phtcnt), plonl, ip, lat, file )
!              end do
!           end if
!-----------------------------------------------------------------------      
!       ...  History output for time averaged reaction rates
!-----------------------------------------------------------------------      
!           do m = 1,hfile(file)%histout_cnt(9,avrg)
!               fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(9,avrg)+m-1)
!               hndx = hfile(file)%timav_map(hfile(file)%histout_ind(9,avrg)+m-1)
!              call outfld( fldname, reaction_rates(1,1,hndx+phtcnt), plonl, ip, lat, file )
!           end do
!        end do
         call adjrxt( reaction_rates, invariants, invariants(:,:,indexm), &
                      plnplv )
      end if
      
      if( phtcnt > 0 ) then
!-----------------------------------------------------------------------
!        ... Compute the photolysis rates at time = t(n+1)
!-----------------------------------------------------------------------      
!        if( polar_night ) then
!           reaction_rates(:,:,1:phtcnt) = 0.
!        else
            esfact = sundis( Time )
            if( ncol_abs > 0 ) then
!-----------------------------------------------------------------------      
!             ... Set the column densities
!-----------------------------------------------------------------------      
               call setcol( col_delta, col_dens, pdel, plonl )
            end if
!-----------------------------------------------------------------------      
!             ... Calculate the surface albedo
!-----------------------------------------------------------------------      
!            call srfalb( lat, ip, albs, caldayn, tsurf, plonl )
!-----------------------------------------------------------------------      
!             ... Calculate the photodissociation rates
!-----------------------------------------------------------------------      
            call photo( reaction_rates(:,:,:phtcnt), pmid, pdel, tfld, zmid, &
                        col_dens, &
!                       zen_angle, albs, &
                        coszen, albedo, &
                        cwat, cldfr, &
!                       sunon, sunoff, &
                        esfact, plonl )
!        end if
!-----------------------------------------------------------------------      
!       ...  History output for instantaneous photo rates
!-----------------------------------------------------------------------      
!        do file = 1,moz_file_cnt
!           if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(8,inst) > 0 ) then
!              do m = 1,hfile(file)%histout_cnt(8,inst)
!                 fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(8,inst)+m-1)
!                  hndx    = hfile(file)%inst_map(hfile(file)%histout_ind(8,inst)+m-1)
!                 call outfld( fldname, reaction_rates(1,1,hndx), plonl, ip, lat, file )
!              end do
!           end if
!-----------------------------------------------------------------------      
!       ...  History output for time averaged photo rates
!-----------------------------------------------------------------------      
!           do m = 1,hfile(file)%histout_cnt(8,avrg)
!               fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(8,avrg)+m-1)
!               hndx    = hfile(file)%timav_map(hfile(file)%histout_ind(8,avrg)+m-1)
!              call outfld( fldname, reaction_rates(1,1,hndx), plonl, ip, lat, file )
!           end do
!        end do
!-----------------------------------------------------------------------      
!             ... Adjust the photodissociation rates
!-----------------------------------------------------------------------      
         call phtadj( reaction_rates, invariants, invariants(:,:,indexm), &
                      plnplv )
      end if
      if( hetcnt > 0 ) then
!-----------------------------------------------------------------------
!        ... Compute the heterogeneous rates at time = t(n+1)
!-----------------------------------------------------------------------      
!        call sethet( het_rates, pmid, lat, zmid, phis, &
!                     tfld, cmfdqr, nrain, nevapr, delt, &
!                     invariants(1,1,indexm), vmr, plonl )
!-----------------------------------------------------------------------      
!       ...  History output for instantaneous wet removal rates
!-----------------------------------------------------------------------      
!        do file = 1,moz_file_cnt
!           if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(10,inst) > 0 ) then
!              do m = 1,hfile(file)%histout_cnt(10,inst)
!                 fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(10,inst)+m-1)
!                     hndx    = hfile(file)%inst_map(hfile(file)%histout_ind(10,inst)+m-1)
!                    call outfld( fldname, het_rates(1,1,hndx), plonl, ip, lat, file )
!              end do
!           end if
!-----------------------------------------------------------------------      
!       ...  History output for time averaged wet removal rates
!-----------------------------------------------------------------------      
!           do m = 1,hfile(file)%histout_cnt(10,avrg)
!               fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(10,avrg)+m-1)
!               hndx    = hfile(file)%timav_map(hfile(file)%histout_ind(10,avrg)+m-1)
!              call outfld( fldname, het_rates(1,1,hndx), plonl, ip, lat, file )
!           end do
!        end do
      end if
      if( extcnt > 0 ) then
!-----------------------------------------------------------------------
!        ... Compute the extraneous frcing at time = t(n+1)
!-----------------------------------------------------------------------      
!        call setext( extfrc, lat, ip, zint, cldtop, plonl )
!-----------------------------------------------------------------------      
!       ...  History output for instantaneous external forcing rates
!-----------------------------------------------------------------------      
!        do file = 1,moz_file_cnt
!           if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(11,inst) > 0 ) then
!              do m = 1,hfile(file)%histout_cnt(11,inst)
!                  fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(11,inst)+m-1)
!                  hndx = hfile(file)%inst_map(hfile(file)%histout_ind(11,inst)+m-1)
!                 call outfld( fldname, extfrc(1,1,hndx), plonl, ip, lat, file )
!              end do
!           end if
!-----------------------------------------------------------------------      
!       ...  History output for time averaged external forcing rates
!-----------------------------------------------------------------------      
!           do m = 1,hfile(file)%histout_cnt(11,avrg)
!               fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(11,avrg)+m-1)
!               hndx    = hfile(file)%timav_map(hfile(file)%histout_ind(11,avrg)+m-1)
!              call outfld( fldname, extfrc(1,1,hndx), plonl, ip, lat, file )
!           end do
!        end do
!        do m = 1,max(1,extcnt)
!           do k = 1,SIZE(vmr,2)
!               extfrc(:,k,m) = extfrc(:,k,m) / invariants(:,k,indexm)
!           end do
!        end do
      end if
      if( grpcnt > 0 ) then
!-----------------------------------------------------------------------
!        ... Set the group ratios
!-----------------------------------------------------------------------      
!        call set_grp_ratios( group_ratios, reaction_rates, vmr, mmr, nas, &
!                              mbar, invariants, plonl )
!-----------------------------------------------------------------------
!             ... Modify the reaction rate of any reaction
!           with group member or proportional reactant(s)
!-----------------------------------------------------------------------
!        call rxt_mod( reaction_rates, het_rates, group_ratios, plnplv )
      end if

!=======================================================================
!        ... Call the class solution algorithms
!=======================================================================
      if( clscnt1 > 0 .and. rxntot > 0 ) then
!-----------------------------------------------------------------------
!        ... Solve for "explicit" species
!-----------------------------------------------------------------------
         call exp_sol( vmr, reaction_rates, &
                       het_rates, extfrc, &
                       nstep, delt, &
!                      invariants(1,1,indexm), &
                       prod_out, loss_out, &
                       plonl, plnplv )
      end if
      if( clscnt4 > 0 .and. rxntot > 0 ) then
!-----------------------------------------------------------------------
!        ... Solve for "Implicit" species
!-----------------------------------------------------------------------
         call imp_sol( vmr, reaction_rates, &
                       het_rates, extfrc, &
                       nstep, delt, &
!                      invariants(1,1,indexm), &
                       lat, lon, &
                       prod_out, loss_out, &
                       plonl, plnplv )
      end if
      if( clscnt5 > 0 .and. rxntot > 0 ) then
!-----------------------------------------------------------------------
!        ... Solve for "Rodas" species
!-----------------------------------------------------------------------
         call rodas_sol( vmr, reaction_rates, &
                         het_rates, extfrc, &
                         nstep, delt, &
!                        invariants(1,1,indexm), &
                         plonl, plnplv )
      end if
!-----------------------------------------------------------------------
!       ... Heterogeneous chemistry
!-----------------------------------------------------------------------
      so2_ndx = get_spc_ndx( 'SO2' )
      so4_ndx = get_spc_ndx( 'SO4' )
      if( so2_ndx > 0 .and. so4_ndx > 0 ) then
         call setsox( pmid, plonl, delt, tfld, sh, &
!                     nrain, nevapr, cmfdqr, &
                      cwat, invariants(:,:,indexm), &
                      vmr )
      end if
!-----------------------------------------------------------------------      
!         ... Check for negative values and reset to zero
!-----------------------------------------------------------------------      
!     call negtrc( lat, 'After chemistry ', vmr, plonl )
!-----------------------------------------------------------------------      
!         ... Set values near upper boundary
!-----------------------------------------------------------------------      
!     call set_ub_vals( lat, ip, vmr, pmid, zmid, &
!                        tfld, caldayn, plonl )
!-----------------------------------------------------------------------      
!         ... Output instantaneous "wet" advected volume mixing
!-----------------------------------------------------------------------      
!     do file = 1,moz_file_cnt
!        if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(1,inst) > 0 ) then
!           do m = 1,hfile(file)%histout_cnt(1,inst)
!               fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(1,inst)+m-1)
!               hndx    = hfile(file)%inst_map(hfile(file)%histout_ind(1,inst)+m-1)
!              call outfld( fldname, vmr(1,1,hndx), plonl, ip, lat, file )
!           end do
!        end if
!-----------------------------------------------------------------------      
!         ... Output time averaged "wet" advected volume mixing ratios
!-----------------------------------------------------------------------      
!        do m = 1,hfile(file)%histout_cnt(1,avrg)
!            fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(1,avrg)+m-1)
!            hndx    = hfile(file)%timav_map(hfile(file)%histout_ind(1,avrg)+m-1)
!           call outfld( fldname, vmr(1,1,hndx), plonl, ip, lat, file )
!        end do
!     end do
!-----------------------------------------------------------------------      
!         ... Output instantaneous "wet" non-advected volume mixing
!-----------------------------------------------------------------------      
!     group_write(:moz_file_cnt) = hfile(:moz_file_cnt)%wrhstts .and. &
!                                  hfile(:moz_file_cnt)%histout_cnt(2,inst) > 0
!     if( ANY( group_write(:moz_file_cnt) ) .or. &
!         ANY( hfile(:moz_file_cnt)%histout_cnt(2,avrg) > 0 ) ) then
!        call mak_grp_vmr( vmr, group_ratios(1,1,1), group_vmr(1,1,1), plonl )
!     end if
!     do file = 1,moz_file_cnt
!        if( hfile(file)%wrhstts .and. hfile(file)%histout_cnt(2,inst) > 0 ) then
!           do m = 1,hfile(file)%histout_cnt(2,inst)
!               fldname = hfile(file)%hist_inst(hfile(file)%histout_ind(2,inst)+m-1)
!               hndx    = hfile(file)%inst_map(hfile(file)%histout_ind(2,inst)+m-1)
!              call outfld( fldname, group_vmr(1,1,hndx), plonl, ip, lat, file )
!           end do
!        end if
!-----------------------------------------------------------------------      
!         ... Output time averaged "wet" non-advected volume mixing ratios
!-----------------------------------------------------------------------      
!        do m = 1,hfile(file)%histout_cnt(2,avrg)
!            fldname = hfile(file)%hist_timav(hfile(file)%histout_ind(2,avrg)+m-1)
!            hndx    = hfile(file)%timav_map(hfile(file)%histout_ind(2,avrg)+m-1)
!           call outfld( fldname, group_vmr(1,1,hndx), plonl, ip, lat, file )
!        end do
!     end do
!-----------------------------------------------------------------------      
!         ... Xform from vmr to mmr
!-----------------------------------------------------------------------      
!     call vmr2mmr( vmr, mmr, nas, group_ratios, mbar, plonl )

      end subroutine chemdr

      end module mo_chemdr_mod
