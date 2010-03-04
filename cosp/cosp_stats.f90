
!---------------------------------------------------------------------
!------------ FMS version number and tagname for this file -----------

! $Id: cosp_stats.f90,v 18.0 2010/03/02 23:29:04 fms Exp $
! $Name: riga $

! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Jul 2008 - A. Bodas-Salcedo - Added capability of producing outputs in standard grid
! Oct 2008 - J.-L. Dufresne   - Bug fixed. Assignment of Npoints,Nlevels,Nhydro,Ncolumns in COSP_STATS
! Oct 2008 - H. Chepfer       - Added PARASOL reflectance arguments
!
! 
MODULE MOD_COSP_STATS
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  USE MOD_LLNL_STATS
  USE MOD_LMD_IPSL_STATS
  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_STATS ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_STATS(gbx,sgx,cfg,sgradar,sglidar,vgrid,stradar,stlidar)
  
   ! Input arguments
   type(cosp_gridbox),intent(in) :: gbx
   type(cosp_subgrid),intent(in) :: sgx
   type(cosp_config),intent(in)  :: cfg
   type(cosp_sgradar),intent(in) :: sgradar
   type(cosp_sglidar),intent(in) :: sglidar
   type(cosp_vgrid),intent(in)   :: vgrid
   ! Output arguments
   type(cosp_radarstats),intent(inout) :: stradar ! Summary statistics for radar
   type(cosp_lidarstats),intent(inout) :: stlidar ! Summary statistics for lidar 
   
   ! Local variables 
   integer :: Npoints  !# of grid points
   integer :: Nlevels  !# of levels
   integer :: Nhydro   !# of hydrometeors
   integer :: Ncolumns !# of columns
   integer :: Nlr
   logical :: ok_lidar_cfad = .false.
   real,dimension(:,:,:),allocatable :: Ze_out,betatot_out,betamol_in,betamol_out,ph_in,ph_out
   real,dimension(:,:),allocatable :: ph_c,betamol_c
 
   Npoints  = gbx%Npoints
   Nlevels  = gbx%Nlevels
   Nhydro   = gbx%Nhydro
   Ncolumns = gbx%Ncolumns
   Nlr      = vgrid%Nlvgrid
  
   if (cfg%Lcfad_lidarsr532) ok_lidar_cfad=.true.

   if (vgrid%use_vgrid) then ! Statistics in a different vertical grid
        allocate(Ze_out(Npoints,Ncolumns,Nlr),betatot_out(Npoints,Ncolumns,Nlr), &
                 betamol_in(Npoints,1,Nlevels),betamol_out(Npoints,1,Nlr),betamol_c(Npoints,Nlr), &
                 ph_in(Npoints,1,Nlevels),ph_out(Npoints,1,Nlr),ph_c(Npoints,Nlr))
        Ze_out = 0.0
        betatot_out  = 0.0
        betamol_in(:,1,:) = sglidar%beta_mol(:,:)
        betamol_out= 0.0
        betamol_c  = 0.0
        ph_in(:,1,:)  = gbx%ph(:,:)
        ph_out  = 0.0
        ph_c    = 0.0
        !++++++++++++ Radar CFAD ++++++++++++++++
        if (cfg%Lradar_sim) then
            call cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,gbx%zlev,gbx%zlev_half,sgradar%Ze_tot, &
                                           Nlr,vgrid%zl,vgrid%zu,Ze_out,log_units=.true.)
            stradar%cfad_ze = cosp_cfad(Npoints,Ncolumns,Nlr,DBZE_BINS,Ze_out, &
                                        DBZE_MIN,DBZE_MAX,CFAD_ZE_MIN,CFAD_ZE_WIDTH)
        endif
        !++++++++++++ Lidar CFAD ++++++++++++++++
        if (cfg%Llidar_sim) then
            call cosp_change_vertical_grid(Npoints,1,Nlevels,gbx%zlev,gbx%zlev_half,betamol_in, &
                                           Nlr,vgrid%zl,vgrid%zu,betamol_out)
            call cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,gbx%zlev,gbx%zlev_half,sglidar%beta_tot, &
                                           Nlr,vgrid%zl,vgrid%zu,betatot_out)
            call cosp_change_vertical_grid(Npoints,1,Nlevels,gbx%zlev,gbx%zlev_half,ph_in, &
                                           Nlr,vgrid%zl,vgrid%zu,ph_out)
            ph_c(:,:) = ph_out(:,1,:)
            betamol_c(:,:) = betamol_out(:,1,:)
            ! Stats from lidar_stat_summary
            call diag_lidar(Npoints,Ncolumns,Nlr,SR_BINS,PARASOL_NREFL &
                            ,betatot_out,betamol_c,sglidar%refl,gbx%land,ph_c &
                            ,LIDAR_UNDEF,ok_lidar_cfad &
                            ,stlidar%cfad_sr,stlidar%srbval &
                            ,LIDAR_NCAT,stlidar%lidarcld,stlidar%cldlayer,stlidar%parasolrefl)
        endif
        !++++++++++++ Lidar-only cloud amount and lidar&radar total cloud mount ++++++++++++++++
        if (cfg%Lradar_sim.and.cfg%Llidar_sim) call cosp_lidar_only_cloud(Npoints,Ncolumns,Nlr, &
                                    betatot_out,betamol_c,Ze_out, &
                                    stradar%lidar_only_freq_cloud,stradar%radar_lidar_tcc)   
        ! Deallocate arrays at coarse resolution
        deallocate(Ze_out,betatot_out,betamol_in,betamol_out,betamol_c,ph_in,ph_out,ph_c)
   else ! Statistics in model levels
        !++++++++++++ Radar CFAD ++++++++++++++++
        if (cfg%Lradar_sim) stradar%cfad_ze = cosp_cfad(Npoints,Ncolumns,Nlr,DBZE_BINS,sgradar%Ze_tot, &
                                        DBZE_MIN,DBZE_MAX,CFAD_ZE_MIN,CFAD_ZE_WIDTH)
        !++++++++++++ Lidar CFAD ++++++++++++++++
        ! Stats from lidar_stat_summary
        if (cfg%Llidar_sim) call diag_lidar(Npoints,Ncolumns,Nlr,SR_BINS,PARASOL_NREFL &
                        ,sglidar%beta_tot,sglidar%beta_mol,sglidar%refl,gbx%land,gbx%ph &
                        ,LIDAR_UNDEF,ok_lidar_cfad &
                        ,stlidar%cfad_sr,stlidar%srbval &
                        ,LIDAR_NCAT,stlidar%lidarcld,stlidar%cldlayer,stlidar%parasolrefl)
        !++++++++++++ Lidar-only cloud amount and lidar&radar total cloud mount ++++++++++++++++
        if (cfg%Lradar_sim.and.cfg%Llidar_sim) call cosp_lidar_only_cloud(Npoints,Ncolumns,Nlr, &
                                    sglidar%beta_tot,sglidar%beta_mol,sgradar%Ze_tot, &
                                    stradar%lidar_only_freq_cloud,stradar%radar_lidar_tcc)   
   endif
   ! Replace undef
   where (stlidar%cfad_sr   == LIDAR_UNDEF) stlidar%cfad_sr   = R_UNDEF 
   where (stlidar%lidarcld  == LIDAR_UNDEF) stlidar%lidarcld  = R_UNDEF 
   where (stlidar%cldlayer  == LIDAR_UNDEF) stlidar%cldlayer  = R_UNDEF 
   where (stlidar%parasolrefl == LIDAR_UNDEF) stlidar%parasolrefl = R_UNDEF 

END SUBROUTINE COSP_STATS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------- SUBROUTINE COSP_CHANGE_VERTICAL_GRID ----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_CHANGE_VERTICAL_GRID(Npoints,Ncolumns,Nlevels,zfull,zhalf,y,M,zl,zu,r,log_units)
   implicit none
   ! Input arguments
   integer,intent(in) :: Npoints  !# of grid points
   integer,intent(in) :: Nlevels  !# of levels
   integer,intent(in) :: Ncolumns !# of columns
   real,dimension(Npoints,Nlevels),intent(in) :: zfull ! Height at model levels [m] (Bottom of model layer)
   real,dimension(Npoints,Nlevels),intent(in) :: zhalf ! Height at half model levels [m] (Bottom of model layer)
   real,dimension(Npoints,Ncolumns,Nlevels),intent(in) :: y     ! Variable to be changed to a different grid
   integer,intent(in) :: M  !# levels in the new grid
   real,dimension(M),intent(in) :: zl ! Lower boundary of new levels  [m]
   real,dimension(M),intent(in) :: zu ! Upper boundary of new levels  [m]
   logical,optional,intent(in) :: log_units ! log units, need to convert to linear units
   ! Output
   real,dimension(Npoints,Ncolumns,M),intent(out) :: r ! Variable on new grid

   ! Local variables
   integer :: i,j,k
   logical :: lunits
   real :: ws
   real,dimension(Nlevels) :: xl,xu ! Lower and upper boundaries of model grid
   real,dimension(M) :: dz          ! Layer depth
   real,dimension(Nlevels,M) :: w   ! Weights to do the mean at each point
   real,dimension(Ncolumns,Nlevels) :: yp  ! Variable to be changed to a different grid.
                                           ! Local copy at a particular point.
                                           ! This allows for change of units.
   
   lunits=.false.
   if (present(log_units)) lunits=log_units
   
   r = R_UNDEF
   do i=1,Npoints
     ! Vertical grid at that point
     xl = zhalf(i,:)
     xu(1:Nlevels-1) = xl(2:Nlevels)
     xu(Nlevels) = zfull(i,Nlevels) +  zfull(i,Nlevels) - zhalf(i,Nlevels) ! Top level symmetric
     dz = zu - zl
     yp = y(i,:,:) ! Temporary variable to regrid
     ! Find weights
     w = 0.0
     do k=1,M
       do j=1,Nlevels
         if ((xl(j) < zl(k)).and.(xu(j) > zl(k)).and.(xu(j) <= zu(k))) then
           !xl(j)-----------------xu(j)
           !      zl(k)------------------------------zu(k)
           w(j,k) = xu(j) - zl(k)
         else if ((xl(j) >= zl(k)).and.(xu(j) <= zu(k))) then
           !           xl(j)-----------------xu(j)
           !      zl(k)------------------------------zu(k)
           w(j,k) = xu(j) - xl(j)
         else if ((xl(j) >= zl(k)).and.(xl(j) < zu(k)).and.(xu(j) >= zu(k))) then
           !                           xl(j)-----------------xu(j)
           !      zl(k)------------------------------zu(k)
           w(j,k) = zu(k) - xl(j)
         else if ((xl(j) <= zl(k)).and.(xu(j) >= zu(k))) then
           !  xl(j)---------------------------xu(j)
           !        zl(k)--------------zu(k)
           w(j,k) = dz(j)
         endif
       enddo
     enddo
     ! Check for dBZ and change if necessary
     if (lunits) then
        where (yp /= R_UNDEF)
          yp = 10.0**(yp/10.0)
        elsewhere
          yp = 0.0
        end where
     endif
     ! Do the weighted mean
     do j=1,Ncolumns
       do k=1,M
          ws = sum(w(:,k))
          if (ws > 0.0) r(i,j,k) = sum(w(:,k)*yp(j,:))/ws
       enddo
     enddo
     ! Check for dBZ and change if necessary
     if (lunits) then
        where (r(i,:,:) <= 0.0)
          r(i,:,:) = R_UNDEF
        elsewhere
          r(i,:,:) = 10.0*log10(r(i,:,:))
        end where
     endif
   enddo
 
 
   
END SUBROUTINE COSP_CHANGE_VERTICAL_GRID 

END MODULE MOD_COSP_STATS
