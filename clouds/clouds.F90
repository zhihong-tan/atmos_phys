
                      module clouds_mod

!=======================================================================
!
!            determines cloud properties necessary for 
!                    fels-schwartzkopf radiation
!
!=======================================================================

use    cloud_rad_mod, only:  cloud_summary
use  cloud_zonal_mod, only:  cloud_zonal
use    cloud_obs_mod, only:  cloud_obs, cloud_obs_init
use time_manager_mod, only:  time_type
use    utilities_mod, only:  error_mesg, FATAL, file_exist,   &
                             check_nml_error, open_file,      &
                             print_version_number, get_my_pe, &
                             close_file
use    rh_clouds_mod, only:  do_rh_clouds, rh_clouds, rh_clouds_avg
use  strat_cloud_mod, only:  do_strat_cloud, strat_cloud_avg

implicit none
private

!------------------- public interfaces ---------------------------------

public   clouds, clouds_init, clouds_end

!-----------------------------------------------------------------------
!--------------------- version number ----------------------------------
            character(len=4), parameter :: vers_num = 'v2.0'
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!   note:  the fels-schwarzkopf radiation code permits bi-spectral
!          cloud reflectivity associated with cloud cdwtr droplets:
!            -->  visible band - (cao3sw and cuvrf);
!            -->  near infra-red band - (cah2sw and cirrf).
!          the f-s code assumes that all gaseous absorption by
!          cdwtr vapor occurs in the near infra-red band.
!          thus, original code contains cbsw and cirab.
!          we shall include cbo3sw and cuvab and let cbsw = cbh2sw.
!          however, these spectral absorptivities will be set to zero.

      real, dimension(3) :: cao3sw = (/ 0.210, 0.450, 0.590 /)
      real, dimension(3) :: cah2sw = (/ 0.210, 0.450, 0.590 /)
      real, dimension(3) :: cbsw   = (/ 0.005, 0.020, 0.035 /)


!-----------------------------------------------------------------------

      logical :: do_init=.true.

!-----------------------------------------------------------------------
!------------------------- namelist ------------------------------------

      logical :: do_zonal_clouds = .false.
      logical :: do_obs_clouds   = .false.
      logical :: do_no_clouds    = .false.
      
      namelist /clouds_nml/ do_zonal_clouds,  &
                            do_obs_clouds,    &
                            do_no_clouds

!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine clouds  (is, js, clear_sky, time, lat, land,       &
                    pfull, phalf, t,  cosz,                   &
                    nclds, ktopsw, kbtmsw, ktoplw, kbtmlw,    &
                    cldamt, cuvrf, cirrf, cirab, emcld, kbot, tca)

!-----------------------------------------------------------------------
        integer, intent(in)                    :: is, js
        logical, intent(in)                    :: clear_sky
type(time_type), intent(in)                    :: time

   real, intent(in), dimension(:,:)    :: lat
   real, intent(in), dimension(:,:)    :: land
   real, intent(in), dimension(:,:,:)  :: pfull,phalf,t
   real, intent(in), dimension(:,:)    :: cosz
integer, intent(out), dimension(:,:)   :: nclds
integer, intent(out), dimension(:,:,:) :: ktopsw,kbtmsw,ktoplw,kbtmlw
   real, intent(out), dimension(:,:,:) :: cldamt,cuvrf,cirrf,cirab,emcld
integer, intent(in),  dimension(:,:),optional :: kbot
   real, intent(inout), dimension(:,:),optional :: tca
!-----------------------------------------------------------------------
   real,dimension(size(cirab,1),size(cirab,2),size(cirab,3)) :: cuvab
integer,dimension(size(ktoplw,1),size(ktoplw,2),size(ktoplw,3)) ::  &
                       ktop, kbtm
   real,dimension(size(pfull,1),size(pfull,2),size(pfull,3)) :: ql,qi,cf,rh
   real,dimension(size(phalf,1),size(phalf,2),size(phalf,3)) :: phaf
integer  i,j,k,kb,kx,kp1,n,ierr
!-----------------------------------------------------------------------
  ierr = 1

  kx  = size(ktopsw,3)-1
  kp1 = kx+1
    
  if (kx /= size(pfull,3)) call error_mesg ('clouds in clouds_mod', &
                       'input arrays have the incorrect size.',FATAL)

!-----------------------------------------------------------------------
!----------- reset value of tca ----------
  
  if (PRESENT(tca)) tca=0.

!-----------------------------------------------------------------------
!----------- default clouds values ----------

      call default_clouds (nclds,ktopsw,kbtmsw,ktoplw,kbtmlw,  &
                           cldamt,cuvrf,cirrf,cuvab,cirab,emcld)

!-------------- no clouds ----------------------------------------------

      if (clear_sky .or. do_no_clouds) then
          if (present(kbot)) call step_mtn_clouds (kx,kbot,          &
                                 nclds,ktopsw,kbtmsw,ktoplw,kbtmlw,  &
                                 cldamt,cuvrf,cirrf,cirab,emcld)
          return
      endif

!-----------------------------------------------------------------------
!--------------- determine rh cloud properties -----------------
 if ( do_rh_clouds() ) then
!-----------------------------------------------------------------------
!---- compute rh_clouds -----

     call rh_clouds_avg (is, js, rh, ierr)

     if (ierr == 0) then
         call rh_clouds(rh,pfull,phalf(:,:,kx+1),cosz,lat*180./3.14159,&
                        nclds,ktop(:,:,2:kp1),kbtm(:,:,2:kp1),  &
                        cldamt(:,:,2:kp1),cuvrf(:,:,2:kp1),  &
                        cirrf(:,:,2:kp1),cuvab(:,:,2:kp1),  &
                        cirab(:,:,2:kp1),emcld(:,:,2:kp1))
     endif     

!-----------------------------------------------------------------------
 endif
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!--------------- determine prognostic cloud properties -----------------
 if ( do_strat_cloud() ) then
!-----------------------------------------------------------------------
      
     call strat_cloud_avg (is, js, ql, qi, cf, ierr)

     if (ierr == 0) then
         call cloud_summary (land,ql,qi,cf,pfull,phalf,t,&
                             cosz,&
                             nclds,ktop(:,:,2:kp1),kbtm(:,:,2:kp1),  &
                             cldamt(:,:,2:kp1),cuvrf(:,:,2:kp1),  &
                             cirrf(:,:,2:kp1),cuvab(:,:,2:kp1),  &
                             cirab(:,:,2:kp1),emcld(:,:,2:kp1))
     endif     

!-----------------------------------------------------------------------
 endif     
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ----------- zonal or observed clouds ??? --------------
!  (also do if avg cloud properties could not be returned)
!-----------------------------------------------------------------------
 if ( (do_zonal_clouds .or. do_obs_clouds) .and. ierr /= 0  ) then
!-----------------------------------------------------------------------

!    ---- constrain phalf to:  0. <= phalf <= 101325. ----
      if (present(kbot)) then
         do k=1,kx+1; do j=1,size(phalf,2); do i=1,size(phalf,1)
            kb=kbot(i,j)
            phaf(i,j,k)=101325.*phalf(i,j,k)/phalf(i,j,kb+1)
         enddo; enddo; enddo
      else
         do k=1,kx+1
            phaf(:,:,k)=101325.*phalf(:,:,k)/phalf(:,:,kx+1)
         enddo
      endif

!   ---- get three cloud levels (store in k=2,4) -----

      call cloud_zonal (time, lat, phaf, nclds,  &
                        ktopsw(:,:,2:4), kbtmsw(:,:,2:4), &
                        ktoplw(:,:,2:4), kbtmlw(:,:,2:4), &
                        cldamt(:,:,2:4), cuvrf(:,:,2:4),  &
                        cirrf(:,:,2:4), cirab(:,:,2:4),  &
                        emcld(:,:,2:4))

      if (do_obs_clouds) call cloud_obs (is, js, time, cldamt(:,:,2:4))

!-----------------------------------------------------------------------
 endif
!-----------------------------------------------------------------------


!----- store longwave and shortwave indices -----

   if ( (do_rh_clouds() .or. do_strat_cloud()) .and. ierr == 0 ) then

         do j=1,size(ktoplw,2)
         do i=1,size(ktoplw,1)
            if (nclds(i,j) > 0) then
               n=nclds(i,j)
               ktoplw(i,j,2:n+1)=ktop(i,j,2:n+1)
               kbtmlw(i,j,2:n+1)=kbtm(i,j,2:n+1)
               ktopsw(i,j,2:n+1)=ktop(i,j,2:n+1)
               kbtmsw(i,j,2:n+1)=kbtm(i,j,2:n+1)+1
            endif
         enddo
         enddo
   endif

!-----------------------------------------------------------------------
!---- total cloud diagnostic ----

     if (present(tca)) call compute_tca_random (nclds,              &
                                                cldamt(:,:,2:kp1),  &
                                                tca                 )

!---- step mountain clouds in underground levels ----

     if (present(kbot)) call step_mtn_clouds (kx,kbot,              &
                                nclds,ktopsw,kbtmsw,ktoplw,kbtmlw,  &
                                cldamt,cuvrf,cirrf,cirab,emcld)

!-----------------------------------------------------------------------

      end subroutine clouds

!#######################################################################

 subroutine compute_tca_random ( nclds, cldamt, tca )

   integer, intent(in)  :: nclds (:,:)
   real,    intent(in)  :: cldamt(:,:,:)
   real,    intent(out) :: tca   (:,:)

   integer :: k, max_cld

!---- find maximum number of clouds -----

    max_cld = maxval(nclds)

!---- compute total cloud amount assuming that -----
!       independent clouds overlap randomly

    tca = 1.0

    do k = 1, max_cld
       tca(:,:) = tca(:,:) * (1. - cldamt(:,:,k))
    enddo

    tca = 1. - tca


 end subroutine compute_tca_random

!#######################################################################

      subroutine default_clouds (nclds,ktopsw,kbtmsw,ktop,kbtm,  &
                                 cldamt,cuvrf,cirrf,cuvab,cirab,emcld)

!----------------------------------------------------------------------
   integer, intent(inout), dimension(:,:)   :: nclds
   integer, intent(inout), dimension(:,:,:) :: ktopsw,kbtmsw,ktop,kbtm
      real, intent(inout), dimension(:,:,:) :: cldamt,cuvrf,cirrf,  &
                                               cuvab,cirab,emcld
!----------------------------------------------------------------------
      integer  kp1

      kp1=size(ktopsw,3)

      nclds(:,:)=0

      cldamt=0.0; emcld =1.0
      cuvrf =0.0; cirrf =0.0; cuvab =0.0; cirab =0.0
      ktopsw=kp1; kbtmsw=kp1
!     ktop  =kp1; kbtm  =kp1
      ktop  =kp1-1; kbtm  =kp1-1
!     ---- reset top properties ----
      ktopsw(:,:,1)=1
      kbtmsw(:,:,1)=0
      ktop  (:,:,1)=1
      kbtm  (:,:,1)=0

!----------------------------------------------------------------------

      end subroutine default_clouds

!#######################################################################

      subroutine step_mtn_clouds (kmax,kbot,nclds,  &
                                  ktopsw,kbtmsw,ktop,kbtm,  &
                                  cldamt,cuvrf,cirrf,cirab,emcld)

!-----------------------------------------------------------------------
   integer, intent(in)                      :: kmax
   integer, intent(in),    dimension(:,:)   :: kbot
   integer, intent(inout), dimension(:,:)   :: nclds
   integer, intent(inout), dimension(:,:,:) :: ktopsw,kbtmsw,ktop,kbtm
      real, intent(inout), dimension(:,:,:) :: cldamt,cuvrf,cirrf,  &
                                                      cirab,emcld
!-----------------------------------------------------------------------
   integer  i,j,n,kp1
integer,save :: jrow=0

   kp1=kmax+1

         do j=1,size(kbot,2)
         do i=1,size(kbot,1)
            if (kbot(i,j) < kmax) then
               nclds(i,j)=nclds(i,j)+1
               n=nclds(i,j)+1
               ktopsw(i,j,n)=kbot(i,j)+1
               kbtmsw(i,j,n)=kp1
               ktop  (i,j,n)=kbot(i,j)+1
               kbtm  (i,j,n)=kmax
               cldamt(i,j,n)=1.0
               cuvrf (i,j,n)=0.0
               cirrf (i,j,n)=0.0
               cirab (i,j,n)=0.0
               emcld (i,j,n)=1.0
            endif
         enddo
         enddo

!-----------------------------------------------------------------------

      end subroutine step_mtn_clouds

!#######################################################################

      subroutine remove_cloud_overlap (nclds,ktop,kbtm)

!-----------------------------------------------------------------------
   integer, intent(in),    dimension(:,:)   :: nclds
   integer, intent(inout), dimension(:,:,:) :: ktop,kbtm
!-----------------------------------------------------------------------
!                   removes sw cloud overlap
!
!    nclds = number of clouds
!    ktop  = sw indice for cloud top
!    kbtm  = sw indice for cloud bottom
!
!-----------------------------------------------------------------------
   integer  i,j,kc

         do j=1,size(nclds,2)
         do i=1,size(nclds,1)
            if (nclds(i,j) <= 1) cycle

            do kc=2,nclds(i,j)
               if (ktop(i,j,kc+1) >= kbtm(i,j,kc)) cycle
!    case 1: thin or thick upper cloud, thick lower cloud
                  if (ktop(i,j,kc+1) <  kbtm(i,j,kc+1)) then
                     ktop(i,j,kc+1)=ktop(i,j,kc+1)+1
                  else
!    case 2: thick upper cloud, thin lower cloud
                     kbtm(i,j,kc)=kbtm(i,j,kc)-1
                  endif
            enddo

         enddo
         enddo

!-----------------------------------------------------------------------

      end subroutine remove_cloud_overlap

!#######################################################################

      subroutine clouds_init (lonb, latb)

!-----------------------------------------------------------------------
         real, intent(in), dimension(:) :: lonb, latb
!-----------------------------------------------------------------------
      integer  unit,io,ierr,id,jd

!-------------- read namelist --------------

      if ( file_exist('input.nml')) then
         unit = open_file (file='input.nml', action='read')
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=clouds_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'clouds_nml')
         enddo
  10     call close_file (unit)
      endif

!      ----- write namelist -----

      unit = open_file ('logfile.out', action='append')
      call print_version_number (unit, 'clouds', vers_num)
      if ( get_my_pe() == 0 ) write (unit, nml=clouds_nml)
      call close_file (unit)


      if (do_obs_clouds) call cloud_obs_init (lonb, latb)

      do_init=.false.

!-----------------------------------------------------------------------

      end subroutine clouds_init

!#######################################################################

      subroutine clouds_end

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

      end subroutine clouds_end

!#######################################################################

end module clouds_mod

