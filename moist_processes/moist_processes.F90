
                    module moist_processes_mod

!-----------------------------------------------------------------------
!
!         interface module for moisture processes
!         ---------------------------------------
!             moist convective adjustment
!             relaxed arakawa-schubert
!             large-scale condensation
!             stratiform prognostic cloud scheme (coming soon)
!
!-----------------------------------------------------------------------

use      moist_conv_mod, only: moist_conv, moist_conv_init
use     lscale_cond_mod, only: lscale_cond, lscale_cond_init
use  sat_vapor_pres_mod, only: tcheck,escomp

use    time_manager_mod, only: time_type

use    diag_manager_mod, only: register_diag_field, send_data

use       utilities_mod, only: file_exist, error_mesg,      &
                               check_nml_error, open_file,  &
                               get_my_pe, FATAL, NOTE,      &
                               close_file

use             ras_mod, only: ras, ras_init

use         dry_adj_mod, only: dry_adj, dry_adj_init

use     strat_cloud_mod, only: strat_init, strat_driv, strat_end, &
                               strat_cloud_sum

use       rh_clouds_mod, only: rh_clouds_init, rh_clouds_end, &
                               rh_clouds_sum

use diag_integral_mod, only:     diag_integral_field_init, &
                             sum_diag_integral_field

use       constants_mod, only: grav, rdgas, rvgas


implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

   public   moist_processes, moist_processes_init, moist_processes_end

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

   integer           :: num_calls, len_pat, tot_pts, num_pts
      real,parameter :: epst=200.
   logical           :: do_init=.true.

   real, parameter :: d622 = rdgas/rvgas
   real, parameter :: d378 = 1.-d622

!--------------------- version number ----------------------------------
   character(len=128) :: version = '$Id: moist_processes.F90,v 1.3 2000/08/04 19:18:46 fms Exp $'
   character(len=128) :: tag = '$Name: bombay $'
!-----------------------------------------------------------------------
!-------------------- namelist data (private) --------------------------


   logical :: do_mca=.true., do_lsc=.true., do_ras=.false.,  &
              do_strat=.false., do_dryadj=.false., &
              do_rh_clouds=.false.,  &
              use_tau=.false.
   character(len=8) :: call_pat = 'x       '

!------ tracer mapping for stratiform cloud scheme ------

   integer :: nql=0, nqi=0, nqa=0

   real :: pdepth = 150.e2
   real :: tfreeze = 273.16

!---------------- namelist variable definitions ------------------------
!
!   do_mca   = switch to turn on/off moist convective adjustment;
!                [logical, default: do_mca=true ]
!   do_lsc   = switch to turn on/off large scale condensation
!                [logical, default: do_lsc=true ]
!   do_ras   = switch to turn on/off relaxed arakawa shubert
!                [logical, default: do_ras=false ]
!   do_strat = switch to turn on/off stratiform cloud scheme
!                [logical, default: do_strat=false ]
!  do_dryadj = switch to turn on/off dry adjustment scheme
!                [logical, default: do_dryadj=false ]
!   use_tau  = switch to determine whether current time level (tau)
!                will be used or else future time level (tau+1).
!                if use_tau = true then the input values for t,q, and r
!                are used; if use_tau = false then input values
!                tm+tdt*dt, etc. are used.
!                [logical, default: use_tau=false ]
!
!   pdepth   = boundary layer depth in pascals for determining mean
!                temperature tfreeze (used for snowfall determination)
!   tfreeze  = mean temperature used for snowfall determination (deg k)
!                [real, default: tfreeze=273.16]
!   nql,nqi,nqa  = tracer number (index) for cloud liquid water,ice,and
!                fraction,respectively.  [integer, default: nql=nqi=nqa=0]
!   notes: 1) do_mca and do_ras   cannot both be true
!          2) do_lsc and do_strat cannot both be true
!          3) pdepth and tfreeze are used to determine liquid vs. solid
!             precipitation for mca, lsc, and ras schemes, the 
!             stratiform scheme determines it's own precipitation type.
!          4) if do_strat=true then the tracer numbers (nql,nqi,nqa) must
!             be set to different values > 0.
!
!-----------------------------------------------------------------------

namelist /moist_processes_nml/ do_mca, do_lsc, do_ras, do_strat,  &
                               do_dryadj, pdepth, tfreeze,        &
                               nql, nqi, nqa, use_tau, call_pat,  &
                               do_rh_clouds

!-----------------------------------------------------------------------
!-------------------- diagnostics fields -------------------------------

integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_tdt_ls  , id_qdt_ls  , id_prec_ls  , id_snow_ls  , &
           id_precip  , id_WVP, id_LWP, id_IWP, &
           id_tdt_dadj, id_rh,  id_mc

character(len=5) :: mod_name = 'moist'

real :: missing_value = -999.

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine moist_processes (is, ie, js, je, Time, dt, land,      &
                            phalf, pfull, omega,                 &
                            t, q, r, u, v, tm, qm, rm, um, vm,   &
                            tdt, qdt, rdt, udt, vdt,             &
                            lprec, fprec, mask, kbot)

!-----------------------------------------------------------------------
!
!    in:  is,ie      starting and ending i indices for window
!
!         js,je      starting and ending j indices for window
!
!         dt         time step (from t(n-1) to t(n+1) if leapfrog)
!                    in seconds   [real]
!
!         land        fraction of surface covered by land
!                     [real, dimension(nlon,nlat)]
!
!         phalf      pressure at half levels in pascals
!                      [real, dimension(nlon,nlat,nlev+1)]
!
!         pfull      pressure at full levels in pascals
!                      [real, dimension(nlon,nlat,nlev)]
!
!         omega      omega (vertical velocity) at full levels
!                    in pascals per second
!                      [real, dimension(nlon,nlat,nlev)]
!
!         t, q       temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         r          tracer fields at full model levels,
!                    at the current time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         u, v,      zonal and meridional wind [m/s] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
! 
!         tm, qm     temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the previous time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rm         tracer fields at full model levels,
!                    at the previous time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         um, vm     zonal and meridional wind [m/s]at full model levels,
!                    at the previous time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!  
! inout:  tdt, qdt   temperature (tdt) [deg k/sec] and specific
!                    humidity of water vapor (qdt) tendency [1/sec]
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rdt        tracer tendencies 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         udt, vdt   zonal and meridional wind tendencies [m/s/s]
! 
!   out:  lprec      liquid precipitiaton rate (rain) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
!
!         fprec      frozen precipitation rate (snow) in kg/m2/s
!                      [real, dimension(nlon,nlat)]
! 
!       optional
!  -----------------
! 
!    in:  mask       mask (1. or 0.) for grid boxes above or below
!                    the ground   [real, dimension(nlon,nlat,nlev)]
!
!         kbot       index of the lowest model level
!                      [integer, dimension(nlon,nlat)]
!
!
!-----------------------------------------------------------------------
integer,         intent(in)              :: is,ie,js,je
type(time_type), intent(in)              :: Time
   real, intent(in)                      :: dt
   real, intent(in) , dimension(:,:,:)   :: phalf, pfull, omega,  &
                                            t, q, u, v, tm, qm, um, vm
   real, intent(in) , dimension(:,:,:,:) :: r, rm
   real, intent(inout),dimension(:,:,:)  :: tdt, qdt, udt, vdt
   real, intent(inout),dimension(:,:,:,:) :: rdt
   real, intent(out), dimension(:,:)   :: lprec, fprec

   real, intent(in)   , dimension(:,:), optional :: land

   real, intent(in) , dimension(:,:,:), optional :: mask
integer, intent(in) , dimension(:,:),   optional :: kbot

!-----------------------------------------------------------------------
real, dimension(size(t,1),size(t,2),size(t,3)) :: tin,qin,ttnd,qtnd
real, dimension(size(t,1),size(t,2))           :: tsnow,snow
logical,dimension(size(t,1),size(t,2))         :: coldT
real, dimension(size(t,1),size(t,2),size(t,3)) :: utnd,vtnd,uin,vin

real, dimension(size(t,1),size(t,2),size(t,3)) :: qltnd,qitnd,qatnd, &
                                                  qlin, qiin, qain
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: mc,mask3
real, dimension(size(t,1),size(t,2),size(t,3)) :: RH, pmass
real, dimension(size(t,1),size(t,2))           :: rain, precip
real, dimension(size(t,1),size(t,2))           :: wvp,lwp,iwp
integer unit,n

integer :: i, j, k, ix, jx, kx, nt, ipat, ip
real    :: dtinv, fdt
logical :: use_mask, do_adjust, used
!-----------------------------------------------------------------------

      if (do_init) call error_mesg ('moist_processes',  &
                     'moist_processes_init has not been called.', FATAL)

!-----------------------------------------------------------------------
!-------- input array size and position in global storage --------------

      ix=size(t,1); jx=size(t,2); kx=size(t,3); nt=size(r,4)

!------- should moisture adjustments be made for this time ??? ---------

      ipat = mod(num_calls,len_pat)+1
      if ( call_pat(ipat:ipat) == 'x' .or. &
           call_pat(ipat:ipat) == 'x' ) then
           do_adjust = .true.
      else
           do_adjust = .false.
           ttnd=0.0; qtnd=0.0; rain=0.0; snow=0.0
      endif

!------- determine time step factor (for strat scheme) ----------
!!    if (do_strat) then
!!       if (do_adjust .and. len_pat > 1) then
!!          fdt=0.
!!          do i=1,len_pat
!!             ip=ipat-i; if (ip < 1) ip=ip+len_pat
!!             fdt=fdt+1.0
!!             if (call_pat(ip:ip) == 'x' .or.  &
!!                 call_pat(ip:ip) == 'x') exit
!!          enddo
!!       else
            fdt=1.
!!       endif
!!    endif
               


!    ----- increment number of calls ----
      num_pts = num_pts + ix*jx
      if (num_pts == tot_pts) then
         num_calls = num_calls +1
         num_pts = 0
      endif

!-----------------------------------------------------------------------

      use_mask=.false.
      if (present(mask) .and. present(kbot))  use_mask=.true.

      dtinv=1./dt
      lprec=0.0; fprec=0.0; precip=0.0
!-----------------------------------------------------------------------
!del if (do_adjust) then
if (do_adjust .or. do_strat) then
!------------------ setup input data -----------------------------------

   if (use_tau) then
      tin (:,:,:)=t(:,:,:)
      qin (:,:,:)=q(:,:,:)
      uin (:,:,:)=u(:,:,:)
      vin (:,:,:)=v(:,:,:)
   else
      tin (:,:,:)=tm(:,:,:)+tdt(:,:,:)*dt
      qin (:,:,:)=qm(:,:,:)+qdt(:,:,:)*dt
      uin (:,:,:)=um(:,:,:)+udt(:,:,:)*dt
      vin (:,:,:)=vm(:,:,:)+vdt(:,:,:)*dt
   endif

   if (use_mask) then
      tin (:,:,:)=mask(:,:,:)*tin(:,:,:)+(1.0-mask(:,:,:))*epst
      qin (:,:,:)=mask(:,:,:)*qin(:,:,:)
      uin (:,:,:)=mask(:,:,:)*uin(:,:,:)
      vin (:,:,:)=mask(:,:,:)*vin(:,:,:)
   endif
   
!------ check for then setup ql & qi & qa, also reset mc -------

   if (do_strat) then
      if (nt == 0 .or. nt < max(nql,nqi,nqa))  call error_mesg ( &
                 'moist_processes', &
                 'number of tracers less than nql or nqi or nqa', FATAL)
      if (nql == nqi .or. nqa == nqi .or. nql == nqa) call error_mesg (&
       'moist_processes',  &
       'tracers indices cannot be the same (i.e., nql=nqi=nqa).', FATAL)

      if (use_tau) then
         qlin (:,:,:)=r(:,:,:,nql)
         qiin (:,:,:)=r(:,:,:,nqi)
         qain (:,:,:)=r(:,:,:,nqa)
      else
         qlin (:,:,:)=rm(:,:,:,nql)+rdt(:,:,:,nql)*dt
         qiin (:,:,:)=rm(:,:,:,nqi)+rdt(:,:,:,nqi)*dt
         qain (:,:,:)=rm(:,:,:,nqa)+rdt(:,:,:,nqa)*dt
      endif

      if (use_mask) then
         qlin (:,:,:)=mask(:,:,:)*qlin(:,:,:)
         qiin (:,:,:)=mask(:,:,:)*qiin(:,:,:)
         qain (:,:,:)=mask(:,:,:)*qain(:,:,:)
      endif
     
      mc(:,:,:) = 0.0

   endif

!------ check temperatures (must be with range of es lookup table) -----

      call tempcheck (is,ie,js,je,tin)

!----------------- mean temp in lower atmosphere -----------------------
!----------------- used for determination of rain vs. snow -------------
!----------------- use input temp ? ------------------------------------

   call tempavg (pdepth, phalf, t, tsnow, mask)

   where (tsnow <= tfreeze)
          coldT=.TRUE.
   elsewhere
          coldT=.FALSE.
   endwhere
      
endif
!-----------------------------------------------------------------------
!***********************************************************************
!----------------- dry adjustment scheme -------------------------------

if (do_dryadj) then

     if (do_adjust) then
         call dry_adj (tin, pfull, phalf, ttnd, mask)
         tin=tin+ttnd
         ttnd=ttnd*dtinv
         tdt = tdt + ttnd
     endif
!------------- diagnostics for dt/dt_dry_adj ---------------------------
     if ( id_tdt_dadj > 0 ) then
        used = send_data ( id_tdt_dadj, ttnd, Time, is, js, 1, &
                           rmask=mask )
     endif
! ----------------------------------------------------------------------
end if

!-----------------------------------------------------------------------
!***********************************************************************
!----------------- moist convective adjustment -------------------------

                        if (do_mca) then
!-----------------------------------------------------------------------

if (do_adjust) then
   
   if (do_strat) then
         call moist_conv (tin,qin,pfull,phalf,coldT,&
                          ttnd,qtnd,rain,snow,Lbot=kbot,&
                          cf=qain,cfdel=qatnd,qldel=qltnd,qidel=qitnd)
   else
         call moist_conv (tin,qin,pfull,phalf,coldT,&
                          ttnd,qtnd,rain,snow,Lbot=kbot)
   end if

!------- update input values and compute tendency -------
                    
      tin=tin+ttnd;    qin=qin+qtnd
      
      ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
      rain=rain*dtinv; snow=snow*dtinv
      
!------- add on tendency ----------
     tdt=tdt+ttnd; qdt=qdt+qtnd

!------- update input values , compute and add on tendency -----------
!-------              in the case of strat                 -----------

      if (do_strat) then
         qlin=qlin+qltnd; qiin=qiin+qitnd; qain=qain+qatnd

         qltnd=qltnd*dtinv; qitnd=qitnd*dtinv; qatnd=qatnd*dtinv
         
         rdt(:,:,:,nql)=rdt(:,:,:,nql)+qltnd(:,:,:)
         rdt(:,:,:,nqi)=rdt(:,:,:,nqi)+qitnd(:,:,:)
         rdt(:,:,:,nqa)=rdt(:,:,:,nqa)+qatnd(:,:,:)

      end if

!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

endif
!-----------------------------------------------------------------------
                           endif

!-----------------------------------------------------------------------
!***********************************************************************
!----------- relaxed arakawa/schubert cumulus param scheme -------------

                        if (do_ras) then
!-----------------------------------------------------------------------
if (do_adjust) then

   if (do_strat) then
      call ras (tin,qin,uin,vin,pfull,phalf,coldT,dt,ttnd,qtnd,&
         utnd,vtnd,rain,snow,kbot,qlin,qiin,qain,mc,qltnd,qitnd,qatnd)
   else if ( id_mc > 0 ) then
      call ras (tin,qin,uin,vin,pfull,phalf,coldT,dt,ttnd,qtnd,&
         utnd,vtnd,rain,snow,kbot=kbot,mc0=mc)
   else
      call ras (tin,qin,uin,vin,pfull,phalf,coldT,dt,ttnd,qtnd,&
         utnd,vtnd,rain,snow,kbot=kbot)
   end if
!------- update input values and compute tendency -----------

      tin=tin+ttnd;    qin=qin+qtnd
      ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
      uin=uin+utnd;    vin=vin+vtnd
      utnd=utnd*dtinv; vtnd=vtnd*dtinv
      rain=rain*dtinv; snow=snow*dtinv

!------- add on tendency ----------
      tdt=tdt+ttnd; qdt=qdt+qtnd
      udt=udt+utnd; vdt=vdt+vtnd

!------- update input values , compute and add on tendency -----------
!-------              in the case of strat                 -----------

      if (do_strat) then
         qlin=qlin+qltnd; qiin=qiin+qitnd; qain=qain+qatnd

         qltnd=qltnd*dtinv; qitnd=qitnd*dtinv; qatnd=qatnd*dtinv
         
         rdt(:,:,:,nql)=rdt(:,:,:,nql)+qltnd(:,:,:)
         rdt(:,:,:,nqi)=rdt(:,:,:,nqi)+qitnd(:,:,:)
         rdt(:,:,:,nqa)=rdt(:,:,:,nqa)+qatnd(:,:,:)

      end if

!------- save total precip and snow ---------
      
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

endif
!-----------------------------------------------------------------------
                           endif
!-----------------------------------------------------------------------
!***********************************************************************
!--------------- DIAGNOSTICS FOR CONVECTIVE SCHEME ---------------------
!-----------------------------------------------------------------------
 if ( do_mca .or. do_ras ) then
!------- diagnostics for dt/dt_ras -------
      if ( id_tdt_conv > 0 ) then
        used = send_data ( id_tdt_conv, ttnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for dq/dt_ras -------
      if ( id_qdt_conv > 0 ) then
        used = send_data ( id_qdt_conv, qtnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for precip_ras -------
      if ( id_prec_conv > 0 ) then
        used = send_data ( id_prec_conv, rain+snow, Time, is, js )
      endif
!------- diagnostics for snow_ras -------
      if ( id_snow_conv > 0 ) then
        used = send_data ( id_snow_conv, snow, Time, is, js )
      endif

!------- diagnostics for cumulus mass flux from ras -------
      if ( id_mc > 0 .and. do_ras ) then
           !------- set up mask --------
           mask3(:,:,1:(kx+1)) = 1.
           if (present(mask)) then
              WHERE (mask(:,:,1:kx) <= 0.5)
                     mask3(:,:,1:kx) = 0.
              END WHERE
           endif
           used = send_data ( id_mc, mc, Time, is, js, 1, rmask=mask3 )
      endif

 endif
!-----------------------------------------------------------------------
!***********************************************************************
!----------------- large-scale condensation ----------------------------

                      if (do_lsc) then
!-----------------------------------------------------------------------
if (do_adjust) then

   call lscale_cond (tin,qin,pfull,phalf,coldT,rain,snow,ttnd,qtnd,&
                     mask=mask)

!------- (update input values and) compute tendency -----
      tin=tin+ttnd;    qin=qin+qtnd
      
      ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
      rain=rain*dtinv; snow=snow*dtinv
      
!------- add on tendency ----------
     tdt=tdt+ttnd; qdt=qdt+qtnd

!------- compute rh clouds if desired ------
     if (do_rh_clouds) then
           
           !calculate relative humidity
           call rh_calc(pfull,tin,qin,RH,mask)
           
           !pass RH to rh_clouds_sum
           call rh_clouds_sum (is, js, RH)
           
     end if

!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

endif

!-----------------------------------------------------------------------
                           endif
!-----------------------------------------------------------------------
!***********************************************************************
!----------------- stratiform precip/cloud scheme ----------------------
                     if (do_strat) then
!-----------------------------------------------------------------------
!del if (do_adjust) then
         call strat_driv (is,ie,js,je,dt*fdt,pfull,phalf,&
                          tin,qin,qlin,qiin,qain,&
                          omega,mc,land,ttnd,qtnd,qltnd,qitnd,qatnd,&
                          rain,snow,mask=mask)
      
!------- update fields before passing to the clouds
      qlin=qlin+qltnd; qiin=qiin+qitnd; qain=qain+qatnd
      tin=tin+ttnd; qin=qin+qtnd
      call strat_cloud_sum (is,js,qlin,qiin,qain)
!     call average_clouds (is,js,qlin,qiin,qain)
!     call average_clouds (is,ie,js,je,qlin,qiin,qain)

!------- compute and add on tendencies ----------

      !note mask multiplication is not necessary since the
      !strat scheme already does this
      ttnd=ttnd*dtinv; qtnd=qtnd*dtinv
      qltnd=qltnd*dtinv; qitnd=qitnd*dtinv; qatnd=qatnd*dtinv
      rain=rain*dtinv; snow=snow*dtinv
    
      tdt=tdt+ttnd; qdt=qdt+qtnd
      rdt(:,:,:,nql)=rdt(:,:,:,nql)+qltnd(:,:,:)
      rdt(:,:,:,nqi)=rdt(:,:,:,nqi)+qitnd(:,:,:)
      rdt(:,:,:,nqa)=rdt(:,:,:,nqa)+qatnd(:,:,:)

!------- save total precip and snow ---------
      lprec=lprec+rain
      fprec=fprec+snow
      precip=precip+rain+snow

!del endif
!-----------------------------------------------------------------------
                           endif
!-----------------------------------------------------------------------
!***********************************************************************
!--------------- DIAGNOSTICS FOR LARGE-SCALE SCHEME --------------------
!-----------------------------------------------------------------------
 if ( do_lsc .or. do_strat ) then
!------- diagnostics for dt/dt_strat -------
      if ( id_tdt_ls > 0 ) then
        used = send_data ( id_tdt_ls, ttnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for dq/dt_strat -------
      if ( id_qdt_ls > 0 ) then
        used = send_data ( id_qdt_ls, qtnd, Time, is, js, 1, &
                           rmask=mask )
      endif
!------- diagnostics for precip_strat -------
      if ( id_prec_ls > 0 ) then
        used = send_data ( id_prec_ls, rain+snow, Time, is, js )
      endif
!------- diagnostics for snow_strat -------
      if ( id_snow_ls > 0 ) then
        used = send_data ( id_snow_ls, snow, Time, is, js )
      endif

 endif
!-----------------------------------------------------------------------
!***********************************************************************
!--------------------- GENERAL DIAGNOSTICS -----------------------------
!-----------------------------------------------------------------------

!------- diagnostics for total precip -------
   if ( id_precip > 0 ) then
        used = send_data ( id_precip, precip, Time, is, js )
   endif


!-----------------------------------------------------------------------
!------- diagnostics for column water vapor, liquid water path and
!------- ice water path

   if ( id_WVP > 0 .or. id_LWP > 0 .or. id_IWP > 0 ) then
        do k=1,kx
          pmass(:,:,k) = (phalf(:,:,k+1)-phalf(:,:,k))/grav
        end do
   endif

!-- compute and write out water vapor path --
   if ( id_WVP > 0 ) then
        wvp(:,:)=0.
        do k=1,kx
          wvp(:,:) = wvp(:,:) + qin(:,:,k)*pmass(:,:,k)
        end do

        used = send_data ( id_WVP, wvp, Time, is, js )
   end if

!-- compute and write out liquid and ice water path --
   if ( id_LWP > 0 .and. do_strat ) then
        lwp(:,:)=0.
        do k=1,kx
          lwp(:,:) = lwp(:,:) + qlin(:,:,k)*pmass(:,:,k)
        end do

        used = send_data ( id_LWP, lwp, Time, is, js )
   end if

!-- compute and write out liquid and ice water path --
   if ( id_IWP > 0 .and. do_strat ) then
        iwp(:,:)=0.
        do k=1,kx
          iwp(:,:) = iwp(:,:) + qiin(:,:,k)*pmass(:,:,k)
        end do

        used = send_data ( id_IWP, iwp, Time, is, js )
   end if

!-----------------------------------------------------------------------
!------- diagnostics for relative humidity -------

   if ( id_rh > 0 ) then
      if (.not.do_rh_clouds) call rh_calc (pfull,tin,qin,RH,mask)
      used = send_data ( id_rh, RH*100., Time, is, js, 1, rmask=mask )
   endif

!-----------------------------------------------------------------------
!---- accumulate global integral of precipiation (mm/day) -----

call sum_diag_integral_field ('prec', precip*86400., is, js)

!-----------------------------------------------------------------------

end subroutine moist_processes

!#######################################################################

subroutine moist_processes_init ( id, jd, kd, axes, Time )

!-----------------------------------------------------------------------
integer,         intent(in) :: id, jd, kd, axes(4)
type(time_type), intent(in) :: Time
!-----------------------------------------------------------------------
!
!      input
!     --------
!
!      id, jd        number of horizontal grid points in the global
!                    fields along the x and y axis, repectively.
!                      [integer]
!
!      kd            number of vertical points in a column of atmosphere
!-----------------------------------------------------------------------

integer  unit,io,ierr,nt
!-----------------------------------------------------------------------

       if ( file_exist('input.nml')) then

         unit = open_file ('input.nml', action='read')
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=moist_processes_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'moist_processes_nml')
         enddo
  10     call close_file (unit)

!------------------- dummy checks --------------------------------------

         if ( do_mca .and. do_ras ) call error_mesg   &
                   ('moist_processes_init',  &
                    'both do_mca and do_ras cannot be specified', FATAL)

         if ( do_lsc .and. do_strat ) call error_mesg   &
                 ('moist_processes_init',  &
                  'both do_lsc and do_strat cannot be specified', FATAL)
         if ( do_rh_clouds .and. do_strat .and. get_my_pe() == 0) &
              call error_mesg ('moist_processes_init',            &
         'both do_rh_clouds and do_strat should not be specified', NOTE)

         if (do_strat) then
         endif

      endif

!----------- determine length of calling pattern ------------

      len_pat=len_trim(call_pat)
      tot_pts=id*jd; num_pts=0

!--------- write namelist ------------------

      unit = open_file ('logfile.out', action='append')
      if ( get_my_pe() == 0 ) then
           write (unit, '(/,80("="),/(a))') trim(version),trim(tag)
           write (unit, nml=moist_processes_nml)
      endif
      call close_file (unit)

!-------- read restart data --------

      if (file_exist('INPUT/moist_processes.res')) then
        unit = open_file ('INPUT/moist_processes.res', &
                          action='read', form='native')
        read  (unit) num_calls
        call close_file (unit)
      else
        num_calls = 0
        if (get_my_pe() == 0) call error_mesg ('moist_processes_init', &
                         'restart data (num_calls) set to zero.',NOTE)
      endif

!------------ initialize various schemes ----------
      if (do_mca)    call  moist_conv_init ()
      if (do_lsc) then
                     call lscale_cond_init ()
                     if (do_rh_clouds) call rh_clouds_init (id,jd,kd)
      endif
      if (do_ras)    call         ras_init ()
      if (do_strat)  call       strat_init (id,jd,kd)
      if (do_dryadj) call     dry_adj_init ()

!----- initialize quantities for global integral package -----

   call diag_integral_field_init ('prec', 'f6.3')

!----- initialize quantities for diagnostics output -----

   call diag_field_init ( axes, Time )


   do_init = .false.

!-----------------------------------------------------------------------

end subroutine moist_processes_init

!#######################################################################

subroutine moist_processes_end

integer  unit
!-----------------------------------------------------------------------

      unit = open_file ('RESTART/moist_processes.res', &
                        action='write', form='native')
      if ( get_my_pe() == 0 ) write (unit)  num_calls
      call close_file (unit)

!----------------close various schemes-----------------

      if (do_strat)     call     strat_end
      if (do_rh_clouds) call rh_clouds_end

!-----------------------------------------------------------------------

end subroutine moist_processes_end

!#######################################################################

      subroutine tempavg (pdepth,phalf,temp,tsnow,mask)

!-----------------------------------------------------------------------
!
!    computes a mean atmospheric temperature for the bottom
!    "pdepth" pascals of the atmosphere.
!
!   input:  pdepth     atmospheric layer in pa.
!           phalf      pressure at model layer interfaces
!           temp       temperature at model layers
!           mask       data mask at model layers (0.0 or 1.0)
!
!   output:  tsnow     mean model temperature in the lowest
!                      "pdepth" pascals of the atmosphere
!
!-----------------------------------------------------------------------
      real, intent(in)  :: pdepth
      real, intent(in) , dimension(:,:,:) :: phalf,temp
      real, intent(out), dimension(:,:)   :: tsnow
      real, intent(in) , dimension(:,:,:), optional :: mask
!-----------------------------------------------------------------------
 real, dimension(size(temp,1),size(temp,2)) :: prsum, done, pdel, pdep
 real  sumdone
 integer  k
!-----------------------------------------------------------------------

      tsnow=0.0; prsum=0.0; done=1.0; pdep=pdepth

      do k=size(temp,3),1,-1

         if (present(mask)) then
           pdel(:,:)=(phalf(:,:,k+1)-phalf(:,:,k))*mask(:,:,k)*done(:,:)
         else
           pdel(:,:)=(phalf(:,:,k+1)-phalf(:,:,k))*done(:,:)
         endif

         where ((prsum(:,:)+pdel(:,:))  >  pdep(:,:))
            pdel(:,:)=pdepth-prsum(:,:)
            done(:,:)=0.0
            pdep(:,:)=0.0
         endwhere

         tsnow(:,:)=tsnow(:,:)+pdel(:,:)*temp(:,:,k)
         prsum(:,:)=prsum(:,:)+pdel(:,:)

         sumdone=sum(done(:,:))
         if (sumdone < 1.e-4) exit

      enddo

         tsnow(:,:)=tsnow(:,:)/prsum(:,:)

!-----------------------------------------------------------------------

      end subroutine tempavg

!#######################################################################

      subroutine tempcheck (is,ie,js,je,temp)

!-----------------------------------------------------------------------
      integer, intent(in)                   :: is,ie,js,je
         real, intent(in), dimension(:,:,:) :: temp
!-----------------------------------------------------------------------
      integer numbad,i,j,k

!------ check temperatures (must be with range of es lookup table) -----

              call tcheck (temp,numbad)

              if (numbad > 0) then

!del               print *, 'pe,numbad,is,js=',get_my_pe(),numbad,is,js

                   do k=1,size(temp,3)
                     call tcheck (temp(:,:,k),numbad)
!del                 print *, 'pe,k,numbad=',get_my_pe(),k,numbad
                   enddo

                 if (size(temp,2) > 1) then
                   do j=1,size(temp,2)
                     call tcheck (temp(:,j,:),numbad)
!del                 print *, 'pe,j,numbad=',get_my_pe(),j,numbad
                   enddo
                 endif

                 call error_mesg ('moist_processes',  &
                   'temperatures not within range of es lookup table', &
                    FATAL)

              endif

!-----------------------------------------------------------------------

      end subroutine tempcheck

!#######################################################################

      subroutine rh_calc(pfull,T,qv,RH,MASK)

        IMPLICIT NONE


        REAL, INTENT (IN),    DIMENSION(:,:,:) :: pfull,T,qv
        REAL, INTENT (OUT),   DIMENSION(:,:,:) :: RH
        REAL, INTENT (IN), OPTIONAL, DIMENSION(:,:,:) :: MASK

        REAL, DIMENSION(SIZE(T,1),SIZE(T,2),SIZE(T,3)) :: esat
        
!-----------------------------------------------------------------------
!       Calculate RELATIVE humidity.
!       This is calculated according to the formula:
!
!       RH   = qv / (epsilon*esat/ [pfull  -  (1.-epsilon)*esat])
!
!       Where epsilon = Rdgas/RVgas = d622
!
!       and where 1- epsilon = d378
!
!       Note that RH does not have its proper value
!       until all of the following code has been executed.  That
!       is, RH is used to store intermediary results
!       in forming the full solution.

        !calculate water saturated vapor pressure from table
        !and store temporarily in the variable esat
        CALL ESCOMP(T(:,:,:),esat(:,:,:))
        
        !calculate denominator in qsat formula
        RH(:,:,:) = pfull(:,:,:)-d378*esat(:,:,:)
     
        !limit denominator to esat, and thus qs to epsilon
        !this is done to avoid blow up in the upper stratosphere
        !where pfull ~ esat
        RH(:,:,:) = MAX(RH(:,:,:),esat(:,:,:)) 
        
        !calculate RH
        RH(:,:,:)=qv(:,:,:)/(d622*esat(:,:,:)/RH(:,:,:))
      
        !IF MASK is present set RH to zero
        IF (present(MASK)) RH(:,:,:)=MASK(:,:,:)*RH(:,:,:)


END SUBROUTINE rh_calc



!#######################################################################

subroutine diag_field_init ( axes, Time )

  integer,         intent(in) :: axes(4)
  type(time_type), intent(in) :: Time

  integer, dimension(3) :: half = (/1,2,4/)

!------------ initializes diagnostic fields in this module -------------

if ( do_mca ) then

   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency from moist conv adj',     'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Spec humidity tendency from moist conv adj',   'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate from moist conv adj',       'kg/m2/s' )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate from moist conv adj',       'kg/m2/s' )

endif

if ( do_ras ) then

   id_tdt_conv = register_diag_field ( mod_name, &
     'tdt_conv', axes(1:3), Time, &
     'Temperature tendency from RAS',                'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_conv = register_diag_field ( mod_name, &
     'qdt_conv', axes(1:3), Time, &
     'Spec humidity tendency from RAS',              'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_conv = register_diag_field ( mod_name, &
     'prec_conv', axes(1:2), Time, &
    'Precipitation rate from RAS',                  'kg/m2/s' )

   id_snow_conv = register_diag_field ( mod_name, &
     'snow_conv', axes(1:2), Time, &
    'Frozen precip rate from RAS',                  'kg/m2/s' )

endif

if ( do_lsc ) then

   id_tdt_ls = register_diag_field ( mod_name, &
     'tdt_ls', axes(1:3), Time, &
       'Temperature tendency from large-scale cond',   'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_ls = register_diag_field ( mod_name, &
     'qdt_ls', axes(1:3), Time, &
     'Spec humidity tendency from large-scale cond', 'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_ls = register_diag_field ( mod_name, &
     'prec_ls', axes(1:2), Time, &
    'Precipitation rate from large-scale cond',     'kg/m2/s' )

   id_snow_ls = register_diag_field ( mod_name, &
     'snow_ls', axes(1:2), Time, &
    'Frozen precip rate from large-scale cond',     'kg/m2/s' )

endif

if ( do_strat ) then

   id_tdt_ls = register_diag_field ( mod_name, &
     'tdt_ls', axes(1:3), Time, &
     'Temperature tendency from strat cloud',        'deg_K/s',  &
                        missing_value=missing_value               )

   id_qdt_ls = register_diag_field ( mod_name, &
     'qdt_ls', axes(1:3), Time, &
     'Spec humidity tendency from strat cloud',      'kg/kg/s',  &
                        missing_value=missing_value               )

   id_prec_ls = register_diag_field ( mod_name, &
     'prec_ls', axes(1:2), Time, &
    'Precipitation rate from strat cloud',          'kg/m2/s' )

   id_snow_ls = register_diag_field ( mod_name, &
     'snow_ls', axes(1:2), Time, &
    'Frozen precip rate from strat cloud',          'kg/m2/s' )

endif

   id_precip = register_diag_field ( mod_name, &
     'precip', axes(1:2), Time, &
     'Total precipitation rate',                     'kg/m2/s' )

   id_WVP = register_diag_field ( mod_name, &
     'WVP', axes(1:2), Time, &
        'Column integrated water vapor',                'kg/m2'  )

if ( do_strat ) then

   id_LWP = register_diag_field ( mod_name, &
     'LWP', axes(1:2), Time, &
        'Liquid water path',                            'kg/m2'   )

   id_IWP = register_diag_field ( mod_name, &
     'IWP', axes(1:2), Time, &
        'Ice water path',                               'kg/m2'   )

endif

   id_tdt_dadj = register_diag_field ( mod_name, &
     'tdt_dadj', axes(1:3), Time, &
   'Temperature tendency from dry conv adj',       'deg_K/s',  &
                        missing_value=missing_value               )

   id_rh = register_diag_field ( mod_name, &
     'rh', axes(1:3), Time, &
         'relative humidity',                            'percent',  & 
                        missing_value=missing_value               )

if ( do_ras ) then

   id_mc = register_diag_field ( mod_name, &
     'mc', axes(half), Time, &
         'Cumulus Mass Flux from RAS',                   'kg/m2/s', &
                        missing_value=missing_value               )

endif
   
!-----------------------------------------------------------------------

end subroutine diag_field_init

!#######################################################################

                 end module moist_processes_mod

