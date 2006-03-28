      module mo_rodas_sol_mod

      use chem_mods_mod, only : clscnt5

      implicit none

      private
      public :: rodas_slv_init, rodas_sol

!     save

      integer :: ox_ndx
      integer :: oh_ndx, ho2_ndx, ch3o2_ndx, po2_ndx, ch3co3_ndx, &
                 c2h5o2_ndx, isopo2_ndx, macro2_ndx, mco3_ndx, c3h7o2_ndx, &
                 ro2_ndx, xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx, &
                 c2h4_ndx, c3h6_ndx, isop_ndx, mvk_ndx, c10h16_ndx, n_ndx
      real :: epsilon(max(1,clscnt5))
      real :: err_wghts(max(1,clscnt5))

character(len=128), parameter :: version     = '$Id: mo_rodas_slv.F90,v 13.0 2006/03/28 21:16:23 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: memphis $'
logical                       :: module_is_initialized = .false.

      contains

      subroutine rodas_slv_init
!-----------------------------------------------------------------------      
!        ... initialize the implict solver
!-----------------------------------------------------------------------      

      use chem_mods_mod,  only : rodas
      use mo_grid_mod,    only : pcnstm1
      use mo_chem_utls_mod, only : get_spc_ndx, get_rxt_ndx

      implicit none

!-----------------------------------------------------------------------      
!        ... local variables
!-----------------------------------------------------------------------      
      real, parameter :: rel_err      = 1.e-2
      real, parameter :: high_rel_err = 1.e-3
      integer :: m
      real    :: eps(pcnstm1)
      real    :: wghts(pcnstm1)

      eps(:) = rel_err
      ox_ndx = get_spc_ndx( 'OX' )
      if( ox_ndx > 0 ) then
         eps(ox_ndx) = high_rel_err
      else
         m = get_spc_ndx( 'O3' )
         if( m > 0 ) then
            eps(m) = high_rel_err
         end if
      end if
      m = get_spc_ndx( 'NO' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'NO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'NO3' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'HNO3' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'HO2NO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'N2O5' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'OH' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if
      m = get_spc_ndx( 'HO2' )
      if( m > 0 ) then
         eps(m) = high_rel_err
      end if

      wghts(:) = 1.
      n_ndx = get_spc_ndx( 'N' )
      if( n_ndx > 0 ) then
         wghts(n_ndx) = 0.
      end if
      do m = 1,rodas%clscnt
         epsilon(m)   = eps(rodas%clsmap(m))
         err_wghts(m) = wghts(rodas%clsmap(m))
      end do

      end subroutine rodas_slv_init

      subroutine rodas_sol( base_sol, reaction_rates, &
                            het_rates, extfrc, &
                            nstep, delt, &
                            plonl, plnplv )
!-----------------------------------------------------------------------
!              ... rodas_sol advances the volumetric mixing ratio
!           forward one time step via the implicit runge-kutta rosenbrock scheme
!-----------------------------------------------------------------------

      use chem_mods_mod,           only : rod_nzcnt, clscnt5, clsze, &
                                      rxntot, hetcnt, extcnt, rodas
      use mo_grid_mod,             only : pcnstm1
      use mo_indprd_mod,           only : indprd
      use mo_rodas_prod_loss_mod,  only : rodas_prod_loss
      use mo_rod_lin_matrix_mod,   only : rod_linmat
      use mo_rod_nln_matrix_mod,   only : rod_nlnmat
      use mo_rod_factor_mod,       only : rod_lu_fac
      use mo_rod_solve_mod,        only : rod_lu_slv

      implicit none

!-----------------------------------------------------------------------
!             ... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::   nstep                     ! time step index (zero based)
      integer, intent(in) ::   plonl                     ! longitude tile dimension
      integer, intent(in) ::   plnplv                    ! plonl*plev
      real, intent(in)    ::   delt                      ! time step (s)
      real, intent(in)    ::   reaction_rates(plnplv,rxntot)
      real, intent(in)    ::   het_rates(plnplv,max(1,hetcnt)), &
                               extfrc(plnplv,max(1,extcnt))
      real, intent(inout) ::   base_sol(plnplv,pcnstm1)

!-----------------------------------------------------------------------
!             ... local variables
!-----------------------------------------------------------------------
      integer, parameter :: att_limit = 5
      real, parameter    :: hmin      = 1.
      real, parameter    :: min_val   = 1.e-30
      real, parameter    :: con3      = 8./3.

      integer ::   i, isec, j, k, m
      integer ::   lev, ofl, ofu
      integer ::   attempts, failures, tsteps, step_fail_cnt
      real    ::   con1, con2
      real, dimension(clsze,max(1,rod_nzcnt)) :: sys_jac, lin_jac
      real, dimension(clsze,max(1,clscnt5))   :: yn, prod, loss, &
                                                 u1, u2, u3, u4, &
                                                 ind_prd
      real, dimension(plnplv,max(1,clscnt5))  :: gl_ind_prd
      real, dimension(clsze,max(1,rxntot))    :: lrxt
      real, dimension(clsze,max(1,hetcnt))    :: lhet
      real, dimension(clsze,max(1,pcnstm1))   :: lsol, y_temp
      real, dimension(max(1,clscnt5))         :: spc_err
      real, dimension(clsze)                  :: err, h_pred
      real    ::   timer 
      real    ::   hfull, hinv, interval
      real    ::   h
      real    ::   hused(2)
      logical ::   interval_done

!-----------------------------------------------------------------------      
!        ... if there is "independent" production put it in the forcing
!        ... set the iteration invariant part of the function f(y)
!-----------------------------------------------------------------------      
      if( rodas%indprd_cnt /= 0 .or. extcnt > 0 ) then
         call indprd( 5, gl_ind_prd, base_sol, extfrc, reaction_rates )
      else
         do m = 1,max(1,clscnt5)
            gl_ind_prd(:,m) = 0.
         end do
      end if
!level_loop : &
!     do lev = 1,plev
lon_tile_loop : &
         do isec = 1,plnplv/clsze
!        do isec = 1,plonl/clsze
!            ofl  = (lev - 1)*plonl + (isec - 1)*clsze + 1
            ofl  = (isec - 1)*clsze + 1
            ofu  = ofl + clsze - 1
            h    = delt
            hinv = 1./h
            interval      = 0.
            step_fail_cnt = 0
            tsteps        = 0
            hused(1)      = 1.e36
            hused(2)      = -1.e36
            do m = 1,rxntot
               lrxt(:,m) = reaction_rates(ofl:ofu,m) 
            end do
            if( hetcnt > 0 ) then
               do m = 1,max(1,hetcnt)
                   lhet(:,m) = het_rates(ofl:ofu,m) 
                end do
            end if
            if( rodas%indprd_cnt /= 0 .or. extcnt > 0 ) then
               do m = 1,clscnt5
                  ind_prd(:,m) = gl_ind_prd(ofl:ofu,m) 
               end do
            end if
!-----------------------------------------------------------------------      
!        ... full timestep loop
!-----------------------------------------------------------------------      
full_time_step_loop : &
            do
               interval_done = .false.
               failures      = 0
               tsteps        = tsteps + 1
               hused(1)      = min( hused(1),h )
               hused(2)      = max( hused(2),h )
!-----------------------------------------------------------------------      
!        ... transfer from base to local work arrays
!-----------------------------------------------------------------------      
               do m = 1,pcnstm1
                  lsol(:,m)   = base_sol(ofl:ofu,m) 
                  y_temp(:,m) = lsol(:,m)
               end do
!----------------------------------------------------------------------      
!        ... store values at t(n)
!-----------------------------------------------------------------------      
               do k = 1,clscnt5
                  j       = rodas%clsmap(k)
                  m       = rodas%permute(k)
                  yn(:,m) = lsol(:,j)
               end do
!-----------------------------------------------------------------------      
!        ... attemp step size
!-----------------------------------------------------------------------      
sub_step_loop : &
               do attempts = 1,att_limit
                  con1 = 2.*hinv
                  con2 = 4.*hinv
!-----------------------------------------------------------------------      
!        ... the linear component
!-----------------------------------------------------------------------      
                  do m = 1,rod_nzcnt
                     lin_jac(:,m) = 0.
                  end do
                  if( rodas%lin_rxt_cnt > 0 ) then
                     call rod_linmat( lin_jac, lsol, lrxt, lhet )
                  else
                     do j = 1,clscnt5
                        m = rodas%diag_map(j)
                        lin_jac(:,m) = -con1
                     end do
                  end if
!-----------------------------------------------------------------------      
!        ... the non-linear component
!-----------------------------------------------------------------------      
                  do m = 1,rod_nzcnt
                     sys_jac(:,m) = 0.
                  end do
                  if( rodas%nln_rxt_cnt > 0 ) then
                     call rod_nlnmat( sys_jac, lsol, lrxt, lin_jac, con1 )
                  else
                     do m = 1,rod_nzcnt
                        sys_jac(:,m) = lin_jac(:,m)
                     end do
                  end if
!-----------------------------------------------------------------------      
!         ... factor the "system" matrix
!-----------------------------------------------------------------------      
                  call rod_lu_fac( sys_jac )
!-----------------------------------------------------------------------      
!           ... form dy/dt = prod - loss
!-----------------------------------------------------------------------      
                  call rodas_prod_loss( prod, loss, lsol, lrxt, lhet )
                   if( rodas%indprd_cnt > 0 .or. extcnt > 0 ) then
                     do m = 1,clscnt5
                        u1(:,m) = loss(:,m) - (prod(:,m) + ind_prd(:,m))
                     end do
                  else
                     do m = 1,clscnt5
                        u1(:,m) = loss(:,m) - prod(:,m)
                     end do
                  end if
                  do m = 1,clscnt5
                     u2(:,m) = u1(:,m)
                  end do
!-----------------------------------------------------------------------      
!           ... solve for the first intermediate
!-----------------------------------------------------------------------      
                  call rod_lu_slv( sys_jac, u1 )
!-----------------------------------------------------------------------      
!           ... solve for the second intermediate
!-----------------------------------------------------------------------      
                  do m = 1,clscnt5
                     u2(:,m) = u2(:,m) - con2*u1(:,m)
                  end do
                  call rod_lu_slv( sys_jac, u2 )
!-----------------------------------------------------------------------      
!           ... solve for the third intermediate
!-----------------------------------------------------------------------      
                  do k = 1,clscnt5
                     j = rodas%clsmap(k)
                     m = rodas%permute(k)
                     lsol(:,j) = yn(:,m) + 2.*u1(:,m)
                  end do
                  call rodas_prod_loss( prod, loss, lsol, lrxt, lhet )
                   if( rodas%indprd_cnt > 0 .or. extcnt > 0 ) then
                     do m = 1,clscnt5
                        u3(:,m) = loss(:,m) - (prod(:,m) + ind_prd(:,m) + hinv*(u1(:,m) - u2(:,m)))
                     end do
                  else
                     do m = 1,clscnt5
                        u3(:,m) = loss(:,m) - (prod(:,m) + hinv*(u1(:,m) - u2(:,m)))
                     end do
                  end if
                  call rod_lu_slv( sys_jac, u3 )
!-----------------------------------------------------------------------      
!           ... solve for the fourth intermediate
!-----------------------------------------------------------------------      
                  do k = 1,clscnt5
                     j = rodas%clsmap(k)
                     m = rodas%permute(k)
                     lsol(:,j) = yn(:,m) + 2.*u1(:,m) + u3(:,m)
                  end do
                  call rodas_prod_loss( prod, loss, lsol, lrxt, lhet )
                   if( rodas%indprd_cnt > 0 .or. extcnt > 0 ) then
                     do m = 1,clscnt5
                        u4(:,m) = loss(:,m) - (prod(:,m) + ind_prd(:,m) + hinv*(u1(:,m) - u2(:,m) - con3*u3(:,m)))
                     end do
                  else
                     do m = 1,clscnt5
                        u4(:,m) = loss(:,m) - (prod(:,m) + hinv*(u1(:,m) - u2(:,m) - con3*u3(:,m)))
                     end do
                  end if
                  call rod_lu_slv( sys_jac, u4 )
!-----------------------------------------------------------------------      
!           ... form y(n+1) from intermediates
!-----------------------------------------------------------------------      
                  do k = 1,clscnt5
                     j = rodas%clsmap(k)
                     m = rodas%permute(k)
                     lsol(:,j) = yn(:,m) + 2.*u1(:,m) + u3(:,m) + u4(:,m)
                  end do
!-----------------------------------------------------------------------      
!           ... form estimated trunc error
!-----------------------------------------------------------------------      
                  err(:) = 0.
                  do k = 1,clscnt5
                     j = rodas%clsmap(k)
                     m = rodas%permute(k)
                     do i = 1,clsze
                        if( lsol(i,j) > min_val ) then
                           spc_err(k) = err_wghts(k) * u4(i,m) / (epsilon(k)*lsol(i,j))
                           err(i)     = err(i) + spc_err(k)*spc_err(k)
                        end if
                     end do
                  end do
                  do i = 1,clsze
                     err(i) = sqrt( err(i)/real(clscnt5) )
                  end do
                  if( all( err(:) < 1. ) ) then
                     if( h == delt ) then
                        interval_done = .true.
                        hfull         = h
                        exit
                     end if
                     interval  = interval + h
                     h_pred(:) = h * min( 10.,max( .1,1./(err(:)**.33) ) )
                     h         = minval( h_pred(:) )
                     h         = max( hmin,h )
                     hfull     = h
                     if( abs( interval - delt ) > 1.e-6*delt ) then
                        h    = min( delt-interval,h )
                        hinv = 1. / h
                     end if
                     exit
                  else
                     if( h == hmin ) then
                        interval      = interval + h
                        hfull         = h
                        step_fail_cnt = step_fail_cnt + 1
                        exit
                     end if
                     failures = failures + 1
                     if( attempts == att_limit ) then
                        interval = interval + h
                     end if
                     if( failures >= 2 ) then
                        h = .1 * h
                     else
                        h_pred(:) = h * min( 10.,max( .1,.5/(err(:)**.33) ) )
                        h         = minval( h_pred(:) )
                     end if
                     h = max( hmin,h )
                     h = min( delt-interval,h )
                     hinv = 1. / h
                     if( attempts == att_limit ) then
                        hfull         = h
                        step_fail_cnt = step_fail_cnt + 1
                        exit
                     end if
                     lsol(:,:) = y_temp(:,:)
                  end if
               end do sub_step_loop
               do m = 1,pcnstm1
                  base_sol(ofl:ofu,m) = lsol(:,m)
               end do
               if( interval_done .or. abs( interval - delt ) <= 1.e-6*delt ) then
                  h = min( hfull,delt )
                  exit
               end if
            end do full_time_step_loop
         end do lon_tile_loop
!     end do level_loop

      end subroutine rodas_sol

      end module mo_rodas_sol_mod
