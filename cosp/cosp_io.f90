
!---------------------------------------------------------------------
!------------ FMS version number and tagname for this file -----------
 
! $Id: cosp_io.f90,v 18.0 2010/03/02 23:28:59 fms Exp $
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
! Jul 2008 - A. Bodas-Salcedo - Initial version
! Oct 2008 - S. Bony - In nc_write_cosp_1d and nc_write_cosp_2d :
!                      the label of layered cloud fractions was wrong -> corrected
!                      (before: low was actually mid, mid was high, high was total,
!                      total was low)
!
 
MODULE MOD_COSP_IO
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
! USE cmor_users_functions
  USE netcdf

  use fms_mod, only: open_namelist_file, open_file, close_file,   &
                     file_exist, mpp_pe, mpp_root_pe,   &
                     check_nml_error, write_version_number, stdlog
  
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  
!---------------------------------------------------------------------
!----------- version number for this module --------------------------     
character(len=128)  :: versiona =  '$Id: cosp_io.f90,v 18.0 2010/03/02 23:28:59 fms Exp $'
character(len=128)  :: tagnamea =  '$Name: riga $'

  ! Types to be used as arrays of pointers
  TYPE var1d
     character(len=16) :: name
     character(len=16) :: units
     integer :: dimsid(3)
     integer :: dimssz(2)
     real,pointer,dimension(:) :: pntr
  END TYPE
  TYPE var2d
     character(len=16) :: name
     character(len=16) :: units
     integer :: dimsid(4)
     integer :: dimssz(3)
     real,pointer,dimension(:,:) :: pntr
  END TYPE
  TYPE var3d
     character(len=16) :: name
     character(len=16) :: units
     integer :: dimsid(5)
     integer :: dimssz(4)
     real,pointer,dimension(:,:,:) :: pntr
  END TYPE
CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE CONSTRUCT_VAR1D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_VAR1D(name,dimsid,dimssz,pntr,y,units)
     ! Input arguments
     character(len=*),intent(in) :: name
     integer,intent(in) :: dimsid(3)
     integer,intent(in) :: dimssz(2)
     real,dimension(:),target,intent(in) :: pntr
     type(var1d),intent(out) :: y
     character(len=*),optional,intent(in) :: units
     
     y%name =  name
     if (present(units)) y%units   =  units
     y%dimsid =  dimsid
     y%dimssz =  dimssz
     y%pntr => pntr
  
  END SUBROUTINE CONSTRUCT_VAR1D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE CONSTRUCT_VAR2D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_VAR2D(name,dimsid,dimssz,pntr,y,units)
     ! Input arguments
     character(len=*),intent(in) :: name
     integer,intent(in) :: dimsid(4)
     integer,intent(in) :: dimssz(3)
     real,dimension(:,:),target,intent(in) :: pntr
     type(var2d),intent(out) :: y
     character(len=*),optional,intent(in) :: units
     
     y%name =  name
     if (present(units)) y%units   =  units
     y%dimsid =  dimsid
     y%dimssz =  dimssz
     y%pntr => pntr
  
  END SUBROUTINE CONSTRUCT_VAR2D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE CONSTRUCT_VAR3D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_VAR3D(name,dimsid,dimssz,pntr,y,units)
     ! Input arguments
     character(len=*),intent(in) :: name
     integer,intent(in) :: dimsid(5)
     integer,intent(in) :: dimssz(4)
     real,dimension(:,:,:),target,intent(in) :: pntr
     type(var3d),intent(out) :: y
     character(len=*),optional,intent(in) :: units
     
     y%name =  name
     if (present(units)) y%units   =  units
     y%dimsid =  dimsid
     y%dimssz =  dimssz
     y%pntr => pntr
  
  END SUBROUTINE CONSTRUCT_VAR3D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE MAP_POINT_TO_LL---------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE MAP_POINT_TO_LL(Nx,Ny,geomode,x1,x2,x3,x4,y2,y3,y4,y5)
     ! Input arguments
     integer,intent(in) :: Nx,Ny,geomode
     real,intent(in),optional :: x1(:),x2(:,:),x3(:,:,:), &
                                 x4(:,:,:,:)
     real,intent(out),optional :: y2(:,:),y3(:,:,:), &
                                  y4(:,:,:,:),y5(:,:,:,:,:)
     ! Local variables
     integer :: Npoints
     integer :: px(Nx*Ny),py(Nx*Ny)
     integer :: i,j,k,l,m
     integer :: Ni,Nj,Nk,Nl
     integer :: Mi,Mj,Mk,Ml,Mm
     character(len=128) :: proname='MAP_POINT_TO_LL'

     Npoints = Nx*Ny
     
     px=0
     py=0
     ! Obtain pointers to do the mapping
     if (geomode == 2) then ! (lon,lat) mode
      do j=1,Ny
        do i=1,Nx
            k = (j-1)*Nx+i
            px(k) = i  
            py(k) = j  
        enddo
      enddo
     else if (geomode == 3) then ! (lon,lat) mode
      do j=1,Nx
        do i=1,Ny
            k = (j-1)*Ny+i
            px(k) = j
            py(k) = i  
        enddo
      enddo
     else
       print *, ' -- '//trim(proname)//': geomode not supported, ',geomode
       stop
     endif

     if (present(x1).and.present(y2)) then
        Ni = size(x1,1)
        Mi = size(y2,1)
        Mj = size(y2,2)
        if (Mi*Mj /= Ni) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 1)'
          stop
        endif
        do i=1,Npoints
          y2(px(i),py(i)) = x1(i)
        enddo
     else if (present(x2).and.present(y3)) then
        Ni = size(x2,1)
        Nj = size(x2,2)
        Mi = size(y3,1)
        Mj = size(y3,2)
        Mk = size(y3,3)
        if (Mi*Mj /= Ni) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 2)'
          stop
        endif
        if (Nj /= Mk) then
          print *, ' -- '//trim(proname)//': Nj /= Mk (opt 2)'
          stop
        endif
        do k=1,Mk
         do i=1,Npoints
            y3(px(i),py(i),k) = x2(i,k)
         enddo
        enddo
     else if (present(x3).and.present(y4)) then
        Ni = size(x3,1)
        Nj = size(x3,2)
        Nk = size(x3,3)
        Mi = size(y4,1)
        Mj = size(y4,2)
        Mk = size(y4,3)
        Ml = size(y4,4)
        if (Mi*Mj /= Ni) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 3)'
          stop
        endif
        if (Nj /= Mk) then
          print *, ' -- '//trim(proname)//': Nj /= Mk (opt 3)'
          stop
        endif
        if (Nk /= Ml) then
          print *, ' -- '//trim(proname)//': Nk /= Ml (opt 3)'
          stop
        endif
        do l=1,Ml
         do k=1,Mk
          do i=1,Npoints
            y4(px(i),py(i),k,l) = x3(i,k,l)
          enddo
         enddo
        enddo
     else if (present(x4).and.present(y5)) then
        Ni = size(x4,1)
        Nj = size(x4,2)
        Nk = size(x4,3)
        Nl = size(x4,4)
        Mi = size(y5,1)
        Mj = size(y5,2)
        Mk = size(y5,3)
        Ml = size(y5,4)
        Mm = size(y5,5)
        if (Mi*Mj /= Ni) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 4)'
          stop
        endif
        if (Nj /= Mk) then
          print *, ' -- '//trim(proname)//': Nj /= Mk (opt 4)'
          stop
        endif
        if (Nk /= Ml) then
          print *, ' -- '//trim(proname)//': Nk /= Ml (opt 4)'
          stop
        endif
        if (Nl /= Mm) then
          print *, ' -- '//trim(proname)//': Nl /= Mm (opt 4)'
          stop
        endif
        do m=1,Mm
         do l=1,Ml
          do k=1,Mk
            do i=1,Npoints
                y5(px(i),py(i),k,l,m) = x4(i,k,l,m)
            enddo
          enddo
         enddo
        enddo
     else
        print *, ' -- '//trim(proname)//': wrong option'
        stop
     endif

     
  END SUBROUTINE MAP_POINT_TO_LL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE MAP_LL_TO_POINT---------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE MAP_LL_TO_POINT(Nx,Ny,Np,x2,x3,x4,x5,y1,y2,y3,y4)
     ! Input arguments
     integer,intent(in) :: Nx,Ny,Np
     real,intent(in),optional :: x2(:,:),x3(:,:,:), &
                                 x4(:,:,:,:),x5(:,:,:,:,:)
     real,intent(out),optional :: y1(:),y2(:,:),y3(:,:,:), &
                                 y4(:,:,:,:)
     ! Local variables
     integer :: px(Nx*Ny),py(Nx*Ny)
     integer :: i,j,k,l,m
     integer :: Ni,Nj,Nk,Nl,Nm
     integer :: Mi,Mj,Mk,Ml
     character(len=128) :: proname='MAP_LL_TO_POINT'
     
     px=0
     py=0
     if (Nx*Ny < Np) then
       print *, ' -- '//trim(proname)//': Nx*Ny < Np'
       stop
     endif
     do j=1,Ny
       do i=1,Nx
          k = (j-1)*Nx+i
          px(k) = i  
          py(k) = j  
       enddo
     enddo
     
     if (present(x2).and.present(y1)) then
        Ni = size(x2,1)
        Nj = size(x2,2)
        Mi = size(y1,1)
        if (Ni*Nj < Mi) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 1)'
          stop
        endif
        do j=1,Np
          y1(j) = x2(px(j),py(j))
        enddo
     else if (present(x3).and.present(y2)) then
        Ni = size(x3,1)
        Nj = size(x3,2)
        Nk = size(x3,3)
        Mi = size(y2,1)
        Mj = size(y2,2)
        if (Ni*Nj < Mi) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 2)'
          stop
        endif
        if (Nk /= Mj) then
          print *, ' -- '//trim(proname)//': Nk /= Mj (opt 2)'
          stop
        endif
        do k=1,Nk
          do j=1,Np
            y2(j,k) = x3(px(j),py(j),k)
          enddo
        enddo
     else if (present(x4).and.present(y3)) then
        Ni = size(x4,1)
        Nj = size(x4,2)
        Nk = size(x4,3)
        Nl = size(x4,4)
        Mi = size(y3,1)
        Mj = size(y3,2)
        Mk = size(y3,3)
        if (Ni*Nj < Mi) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 3)'
          stop
        endif
        if (Nk /= Mj) then
          print *, ' -- '//trim(proname)//': Nk /= Mj (opt 3)'
          stop
        endif
        if (Nl /= Mk) then
          print *, ' -- '//trim(proname)//': Nl /= Mk (opt 3)'
          stop
        endif
        do l=1,Nl
         do k=1,Nk
          do j=1,Np
            y3(j,k,l) = x4(px(j),py(j),k,l)
          enddo
         enddo
        enddo
     else if (present(x5).and.present(y4)) then
        Ni = size(x5,1)
        Nj = size(x5,2)
        Nk = size(x5,3)
        Nl = size(x5,4)
        Nm = size(x5,5)
        Mi = size(y4,1)
        Mj = size(y4,2)
        Mk = size(y4,3)
        Ml = size(y4,4)
        if (Ni*Nj < Mi) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 4)'
          stop
        endif
        if (Nk /= Mj) then
          print *, ' -- '//trim(proname)//': Nk /= Mj (opt 4)'
          stop
        endif
        if (Nl /= Mk) then
          print *, ' -- '//trim(proname)//': Nl /= Mk (opt 4)'
          stop
        endif
        if (Nm /= Ml) then
          print *, ' -- '//trim(proname)//': Nm /= Ml (opt 4)'
          stop
        endif
        do m=1,Nm
         do l=1,Nl
          do k=1,Nk
           do j=1,Np
            y4(j,k,l,m) = x5(px(j),py(j),k,l,m)
           enddo
          enddo
         enddo
        enddo
     else
        print *, ' -- '//trim(proname)//': wrong option'
        stop
     endif
  
  END SUBROUTINE MAP_LL_TO_POINT
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINE NC_READ_INPUT_FILE -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE NC_READ_INPUT_FILE(fname,Npnts,Nl,Nhydro,lon,lat,p,ph,z,zh,T,qv,rh,tca,cca, &
            mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl, &
            fl_ccrain,fl_ccsnow,Reff,dtau_s,dtau_c,dem_s,dem_c,skt,landmask,sfc_height, &
            mr_ozone,u_wind,v_wind,emsfc_lw,mode,Nlon,Nlat)
    
    !Arguments
    character(len=512),intent(in) :: fname ! File name
    integer,intent(in) :: Npnts,Nl,Nhydro
    real,dimension(Npnts),intent(out) :: lon,lat
    real,dimension(Npnts,Nl),target,intent(out) :: p,ph,z,zh,T,qv,rh,tca,cca, &
                  mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl, &
                  fl_ccrain,fl_ccsnow,dtau_s,dtau_c,dem_s,dem_c,mr_ozone
    real,dimension(Npnts,Nl,Nhydro),intent(out) :: Reff
    real,dimension(Npnts),intent(out) :: skt,landmask,sfc_height,u_wind,v_wind
    real,intent(out) :: emsfc_lw
    integer,intent(out) :: mode,Nlon,Nlat
    
        
    !Local variables
    integer :: Npoints,Nlevels,i,j,k
    character(len=128) :: vname
    integer,parameter :: NMAX_DIM=5
    integer :: vrank,vdimid(NMAX_DIM)
    character(len=256) :: dimname(NMAX_DIM) ! 256 hardcoded, instead of MAXNCNAM. This works for NetCDF 3 and 4.
    integer :: ncid,vid,ndims,nvars,ngatts,recdim,dimsize(NMAX_DIM)
    integer :: errst
    logical :: Llat,Llon,Lpoint
    integer :: Na,Nb,Nc,Nd,Ne
    real,dimension(Npnts) :: ll
    integer,dimension(:),allocatable :: plon,plat
    real,allocatable :: x1(:),x2(:,:),x3(:,:,:),x4(:,:,:,:),x5(:,:,:,:,:) ! Temporary arrays
    
    mode = 0
    Nlon = 0
    Nlat = 0
    
    Npoints = Npnts
    Nlevels = Nl
    
    ! Open file
    errst = nf90_open(fname, nf90_nowrite, ncid)
    
    ! Get information about dimensions. Curtain mode or lat/lon mode?
    Llat  =.false.
    Llon  =.false.
    Lpoint=.false.
    errst = nf90_inquire(ncid, ndims, nvars, ngatts, recdim)
    do i = 1,ndims
       errst = nf90_Inquire_Dimension(ncid,i,name=dimname(i),len=dimsize(i))
       if ((trim(dimname(i)).eq.'level').and.(Nlevels > dimsize(i))) then
         print *,' --- NC_READ_INPUT_FILE: number of levels selected is greater than in input file '//trim(fname)
         stop
       endif
       if (trim(dimname(i)).eq.'point') then
         Lpoint = .true.
         if (Npnts > dimsize(i)) then
           print *,' --- NC_READ_INPUT_FILE: number of points selected is greater than in input file '//trim(fname)
           stop
         endif
       endif
       if (trim(dimname(i)).eq.'lon') then
         Llon = .true.
         Nlon = dimsize(i)
       endif
       if (trim(dimname(i)).eq.'lat') then
         Llat = .true.
         Nlat = dimsize(i)
       endif
    enddo

    ! Get lon and lat
    if (Llon.and.Llat) then ! 2D mode
        if ((Npnts) > Nlon*Nlat) Npoints=Nlon*Nlat
        lon = R_UNDEF
        lat = R_UNDEF
        mode = 2 ! Don't know yet if (lon,lat) or (lat,lon) at this point
    else if (Lpoint) then ! 1D mode
        Nlon = Npoints
        Nlat = Npoints
        mode = 1
    else
        print *, ' -- NC_READ_INPUT_FILE: '//trim(fname)//' file contains wrong dimensions'
        stop
    endif
    errst = nf90_inq_varid(ncid, 'lon', vid)
    errst = nf90_get_var(ncid, vid, lon, start = (/1/), count = (/Nlon/))
    errst = nf90_inq_varid(ncid, 'lat', vid)
    errst = nf90_get_var(ncid, vid, lat, start = (/1/), count = (/Nlat/))
    
    ! Get all variables
    do vid = 1,nvars
       vdimid=0
       errst = nf90_Inquire_Variable(ncid, vid, name=vname, ndims=vrank, dimids=vdimid)
       ! Read in into temporary array of correct shape
       if (vrank == 1) then
          Na = dimsize(vdimid(1))
          allocate(x1(Na))
          errst = nf90_get_var(ncid, vid, x1, start=(/1/), count=(/Na/))
       endif
       if (vrank == 2) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          allocate(x2(Na,Nb))
          errst = nf90_get_var(ncid, vid, x2, start=(/1,1/), count=(/Na,Nb/))
       endif
       if (vrank == 3) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          Nc = dimsize(vdimid(3))
          allocate(x3(Na,Nb,Nc))
          errst = nf90_get_var(ncid, vid, x3, start=(/1,1,1/), count=(/Na,Nb,Nc/))
          if ((mode == 2).or.(mode == 3)) then
            if ((Na == Nlon).and.(Nb == Nlat)) then
              mode = 2
            else if ((Na == Nlat).and.(Nb == Nlon)) then
              mode = 3
            else
              print *, '  -- NC_READ_INPUT_FILE: wrong mode for variable '//trim(vname)
              stop
            endif
          endif
       endif
       if (vrank == 4) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          Nc = dimsize(vdimid(3))
          Nd = dimsize(vdimid(4))
          allocate(x4(Na,Nb,Nc,Nd))
          errst = nf90_get_var(ncid, vid, x4, start=(/1,1,1,1/), count=(/Na,Nb,Nc,Nd/))
       endif
       if (vrank == 5) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          Nc = dimsize(vdimid(3))
          Nd = dimsize(vdimid(4))
          Ne = dimsize(vdimid(5))
          allocate(x5(Na,Nb,Nc,Nd,Ne))
          errst = nf90_get_var(ncid, vid, x5, start=(/1,1,1,1,1/), count=(/Na,Nb,Nc,Nd,Ne/))
       endif
       ! Map to the right input argument
       print *, 'Reading '//trim(vname)//' ...'
       select case (trim(vname))
       case ('pfull')
         if (Lpoint) then
           p(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=p)
         endif
       case ('phalf')
         if (Lpoint) then
           ph(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=ph)
         endif
       case ('height')
         if (Lpoint) then
           z(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=z)
         endif
       case ('height_half')
         if (Lpoint) then
           zh(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=zh)
         endif
       case ('T_abs')
         if (Lpoint) then
           T(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=T)
         endif
       case ('qv')
         if (Lpoint) then
           qv(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=qv)
         endif
       case ('rh')
         if (Lpoint) then
           rh(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=rh)
         endif
       case ('tca')
         if (Lpoint) then
           tca(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=tca)
         endif
       case ('cca')
         if (Lpoint) then
           cca(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=cca)
         endif
       case ('mr_lsliq')
         if (Lpoint) then
           mr_lsliq(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_lsliq)
         endif
       case ('mr_lsice')
         if (Lpoint) then
           mr_lsice(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_lsice)
         endif
       case ('mr_ccliq')
         if (Lpoint) then
           mr_ccliq(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ccliq)
         endif
       case ('mr_ccice')
         if (Lpoint) then
           mr_ccice(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ccice)
         endif
       case ('fl_lsrain')
         if (Lpoint) then
           fl_lsrain(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lsrain)
         endif
       case ('fl_lssnow')
         if (Lpoint) then
           fl_lssnow(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lssnow)
         endif
       case ('fl_lsgrpl')
         if (Lpoint) then
           fl_lsgrpl(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lsgrpl)
         endif
       case ('fl_ccrain')
         if (Lpoint) then
           fl_ccrain(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_ccrain)
         endif
       case ('fl_ccsnow')
         if (Lpoint) then
           fl_ccsnow(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_ccsnow)
         endif
       case ('dtau_s')
         if (Lpoint) then
           dtau_s(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dtau_s)
         endif
       case ('dtau_c')
         if (Lpoint) then
           dtau_c(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dtau_c)
         endif
       case ('dem_s')
         if (Lpoint) then
           dem_s(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dem_s)
         endif
       case ('dem_c')
         if (Lpoint) then
           dem_c(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dem_c)
         endif
       case ('Reff')
         if (Lpoint) then
           Reff(1:Npoints,:,:) = x3(1:Npoints,1:Nlevels,:)
         else
           call map_ll_to_point(Na,Nb,Npoints,x4=x4,y3=Reff)
         endif
       case ('skt')
         if (Lpoint) then
           skt(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=skt)
         endif
       case ('landmask')
         if (Lpoint) then
           landmask(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=landmask)
         endif
       case ('orography')
         if (Lpoint) then
           sfc_height(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=sfc_height)
         endif
       case ('mr_ozone')
         if (Lpoint) then
           mr_ozone(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
         else
           call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ozone)
         endif
       case ('u_wind')
         if (Lpoint) then
           u_wind(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=u_wind)
         endif
       case ('v_wind')
         if (Lpoint) then
           v_wind(1:Npoints) = x1(1:Npoints)
         else
           call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=v_wind)
         endif
       end select
       ! Free memory
       if (vrank == 1) deallocate(x1)
       if (vrank == 2) deallocate(x2)
       if (vrank == 3) deallocate(x3)
       if (vrank == 4) deallocate(x4)
       if (vrank == 5) deallocate(x5)
    enddo
       
    ! SFC emissivity
    errst = nf90_inq_varid(ncid, 'emsfc_lw', vid)
    errst = nf90_get_var(ncid, vid, emsfc_lw)
    
    ! Fill in the lat/lon vectors with the right values for 2D modes
    ! This might be helpful if the inputs are 2D (gridded) and 
    ! you want outputs in 1D mode
    allocate(plon(Npoints),plat(Npoints))
    if (mode == 2) then !(lon,lat)
      ll = lat
      do j=1,Nb
        do i=1,Na
          k = (j-1)*Na + i
          plon(k) = i  
          plat(k) = j
        enddo
      enddo
      lon(1:Npoints) = lon(plon(1:Npoints))
      lat(1:Npoints) = ll(plat(1:Npoints))
    else if (mode == 3) then !(lat,lon)
      ll = lon
      do j=1,Nb
        do i=1,Na
          k = (j-1)*Na + i
          lon(k) = ll(j)
          lat(k) = lat(i)
        enddo
      enddo
      lon(1:Npoints) = ll(plon(1:Npoints))
      lat(1:Npoints) = lat(plat(1:Npoints))
    endif
    deallocate(plon,plat)
    
    ! Close file
!     call ncclos(ncid,errst)
    errst = nf90_close(ncid)

  END SUBROUTINE NC_READ_INPUT_FILE


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE READ_COSP_OUTPUT_NL -------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SUBROUTINE READ_COSP_OUTPUT_NL(cosp_nl,cfg)
  character(len=*),intent(in) :: cosp_nl
  type(cosp_config),intent(out) :: cfg
  ! Local variables
  integer :: i
  logical ::   Lradar_sim,Llidar_sim,Lisccp_sim,Lmisr_sim,Lrttov_sim, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,Lcfad_dbze94, &
             Lcfad_lidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp2,Lcllcalipso, &
             Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lctpisccp,Ldbze94,Ltauisccp,Ltclisccp, &
             Llongitude,Llatitude,Lparasol_refl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfrac_out,Lbeta_mol532,Ltbrttov
  namelist/COSP_OUTPUT/Lradar_sim,Llidar_sim,Lisccp_sim,Lmisr_sim,Lrttov_sim, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,Lcfad_dbze94, &
             Lcfad_lidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp2, &
             Lcllcalipso,Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lctpisccp,Ldbze94,Ltauisccp, &
             Ltclisccp,Llongitude,Llatitude,Lparasol_refl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfrac_out,Lbeta_mol532,Ltbrttov

  integer :: unit, io, ierr
  
  do i=1,N_OUT_LIST
    cfg%out_list(i)=''
  enddo
! open(10,file=cosp_nl,status='old')
! read(10,nml=cosp_output)
! close(10)
!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
    if ( file_exist('input.nml')) then
      unit =  open_namelist_file ()
      ierr=1; do while (ierr /= 0)
      read  (unit, nml=cosp_output, iostat=io, end=10)
      ierr = check_nml_error(io,'cosp_output')
      enddo
10      call close_file (unit)
    endif
 
!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
     call write_version_number (versiona, tagnamea)
     if (mpp_pe() == mpp_root_pe() )    &
                        write (stdlog(), nml=cosp_output)

  
  ! Deal with dependencies
  if (.not.Lradar_sim) then
    Lcfad_dbze94   = .false.
    Lclcalipso2    = .false.
    Lcltlidarradar = .false.
    Ldbze94        = .false.
  endif
  if (.not.Llidar_sim) then
    Latb532 = .false.
    Lcfad_lidarsr532 = .false.
    Lclcalipso2      = .false.
    Lclcalipso       = .false.
    Lclhcalipso      = .false.
    Lcllcalipso      = .false.
    Lclmcalipso      = .false.
    Lcltcalipso      = .false.
    Lcltlidarradar   = .false.
    Lparasol_refl    = .false.
    Lbeta_mol532     = .false.
  endif
  if (.not.Lisccp_sim) then
    Lalbisccp       = .false.
    Lboxptopisccp   = .false.
    Lboxtauisccp    = .false.
    Lclisccp2       = .false.
    Lctpisccp       = .false.
    Ltauisccp       = .false.
    Ltclisccp       = .false.
    Lmeantbisccp    = .false.
    Lmeantbclrisccp = .false.
  endif
  if (.not.Lmisr_sim) then
    LclMISR = .false.
  endif
  if (.not.Lrttov_sim) then
    Ltbrttov = .false.
  endif
  if ((.not.Lradar_sim).and.(.not.Llidar_sim).and. &
      (.not.Lisccp_sim).and.(.not.Lmisr_sim)) then
    Lfrac_out = .false.
  endif
  
  ! Diagnostics that use Radar and Lidar
  if (((Lclcalipso2).or.(Lcltlidarradar)).and.((Lradar_sim).or.(Llidar_sim))) then
    Lclcalipso2    = .true.
    Lcltlidarradar = .true.
    Llidar_sim     = .true.
    Lradar_sim     = .true.
  endif
  
  cfg%Lstats = .false.
  if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) cfg%Lstats = .true.
  
  ! Copy instrument flags to cfg structure
  cfg%Lradar_sim = Lradar_sim
  cfg%Llidar_sim = Llidar_sim
  cfg%Lisccp_sim = Lisccp_sim
  cfg%Lmisr_sim  = Lmisr_sim
  cfg%Lrttov_sim = Lrttov_sim
  
  ! Flag to control output to file
  cfg%Lwrite_output = .false.
  if (cfg%Lstats.or.cfg%Lmisr_sim.or.cfg%Lrttov_sim) then
    cfg%Lwrite_output = .true.
  endif
  
  ! Output diagnostics
  i = 1
  if (Lalbisccp)        cfg%out_list(i) = 'albisccp'
  i = i+1
  if (Latb532)          cfg%out_list(i) = 'atb532'
  i = i+1
  if (Lboxptopisccp)    cfg%out_list(i) = 'boxptopisccp'
  i = i+1
  if (Lboxtauisccp)     cfg%out_list(i) = 'boxtauisccp'
  i = i+1
  if (Lcfad_dbze94)     cfg%out_list(i) = 'cfad_dbze94'
  i = i+1
  if (Lcfad_lidarsr532) cfg%out_list(i) = 'cfad_lidarsr532'
  i = i+1
  if (Lclcalipso2)      cfg%out_list(i) = 'clcalipso2'
  i = i+1
  if (Lclcalipso)       cfg%out_list(i) = 'clcalipso'
  i = i+1
  if (Lclhcalipso)      cfg%out_list(i) = 'clhcalipso'
  i = i+1
  if (Lclisccp2)        cfg%out_list(i) = 'clisccp2'
  i = i+1
  if (Lcllcalipso)      cfg%out_list(i) = 'cllcalipso'
  i = i+1
  if (Lclmcalipso)      cfg%out_list(i) = 'clmcalipso'
  i = i+1
  if (Lcltcalipso)      cfg%out_list(i) = 'cltcalipso'
  i = i+1
  if (Lcltlidarradar)   cfg%out_list(i) = 'cltlidarradar'
  i = i+1
  if (Lctpisccp)        cfg%out_list(i) = 'ctpisccp'
  i = i+1
  if (Ldbze94)          cfg%out_list(i) = 'dbze94'
  i = i+1
  if (Ltauisccp)        cfg%out_list(i) = 'tauisccp'
  i = i+1
  if (Ltclisccp)        cfg%out_list(i) = 'tclisccp'
  i = i+1
  if (Llongitude)       cfg%out_list(i) = 'lon'
  i = i+1
  if (Llatitude)        cfg%out_list(i) = 'lat'
  i = i+1
  if (Lparasol_refl)    cfg%out_list(i) = 'parasol_refl'
  i = i+1
  if (LclMISR)          cfg%out_list(i) = 'clMISR'
  i = i+1
  if (Lmeantbisccp)     cfg%out_list(i) = 'meantbisccp'
  i = i+1
  if (Lmeantbclrisccp)  cfg%out_list(i) = 'meantbclrisccp'
  i = i+1
  if (Lfrac_out)        cfg%out_list(i) = 'frac_out'
  i = i+1
  if (Lbeta_mol532)     cfg%out_list(i) = 'beta_mol532'
  i = i+1
  if (Ltbrttov)         cfg%out_list(i) = 'tbrttov'

  if (i /= N_OUT_LIST) then
     print *, 'COSP_IO: wrong number of output diagnostics'
     stop
  endif

  ! Copy diagnostic flags to cfg structure
  cfg%Lalbisccp = Lalbisccp
  cfg%Latb532 = Latb532
  cfg%Lboxptopisccp = Lboxptopisccp
  cfg%Lboxtauisccp = Lboxtauisccp
  cfg%Lcfad_dbze94 = Lcfad_dbze94
  cfg%Lcfad_lidarsr532 = Lcfad_lidarsr532
  cfg%Lclcalipso2 = Lclcalipso2
  cfg%Lclcalipso = Lclcalipso
  cfg%Lclhcalipso = Lclhcalipso
  cfg%Lclisccp2 = Lclisccp2
  cfg%Lcllcalipso = Lcllcalipso
  cfg%Lclmcalipso = Lclmcalipso
  cfg%Lcltcalipso = Lcltcalipso
  cfg%Lcltlidarradar = Lcltlidarradar
  cfg%Lctpisccp = Lctpisccp
  cfg%Ldbze94 = Ldbze94
  cfg%Ltauisccp = Ltauisccp
  cfg%Ltclisccp = Ltclisccp
  cfg%Llongitude = Llongitude
  cfg%Llatitude = Llatitude
  cfg%Lparasol_refl = Lparasol_refl
  cfg%LclMISR = LclMISR
  cfg%Lmeantbisccp = Lmeantbisccp
  cfg%Lmeantbclrisccp = Lmeantbclrisccp
  cfg%Lfrac_out = Lfrac_out
  cfg%Lbeta_mol532 = Lbeta_mol532
  cfg%Ltbrttov = Ltbrttov
  
 END SUBROUTINE READ_COSP_OUTPUT_NL
   
END MODULE MOD_COSP_IO
