
                module surface_albedo_mod


use rad_output_file_mod,   only: hold_sfcalbedo
use utilities_mod,         only: open_file, file_exist,    &
                                 check_nml_error, error_mesg, &
                                 print_version_number, FATAL, NOTE, &
				 WARNING, get_my_pe, close_file, &
				 get_domain_decomp
use radiation_diag_mod,    only: radiag_from_sfcalbedo
use rad_utilities_mod,     only: Rad_control, radiation_control_type, &
				 Environment, environment_type
use astronomy_package_mod, only: get_astronomy_for_sfcalbedo
use rad_step_setup_mod,    only: jabs, iabs, snow, sealp,  &
                                 IMINP,IMAXP, JMINP, JMAXP, &
				 albedo_sv
use store_rad_output_mod,  only: store_sfcalbedo_wrap

!-------------------------------------------------------------------

implicit none
private

!------------------------------------------------------------------
!                    surface albedo module
!
!-------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!  character(len=5), parameter  ::  version_number = 'v0.09'
   character(len=128)  :: version =  '$Id: surface_albedo.F90,v 1.2 2001/07/05 17:33:39 fms Exp $'
   character(len=128)  :: tag     =  '$Name: galway $'



!---------------------------------------------------------------------
!-------  interfaces --------

public     &
        surface_albedo_init, get_surfacealbedo_for_swrad,  &
        sfcalbedo, sfcalbedo_alloc, sfcalbedo_dealloc


!--------------------------------------------------------------------
!----- namelist -------

real               ::   single_column_albedo = 0.1  
character(len=10)  ::   albedo_type='climap'

namelist /surface_albedo_nml/    &
                                     albedo_type, &
				     single_column_albedo



!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

real, dimension(:,:), allocatable :: cirabgd, cirrfgd, cvisrfgd, &
				     zoalnew, zoalnew_gl
real, dimension(:),   allocatable :: zoal, zoas, zoal_gl, zoas_gl


real, dimension(19)               :: oalm, oasm
integer                           :: jj

data (oalm(jj),jj=1,19)/    &
        0.160E+00, 0.160E+00, 0.160E+00, 0.160E+00, 0.160E+00,  &
        0.149E+00, 0.168E+00, 0.167E+00, 0.100E+00, 0.080E+00, &
        0.100E+00, 0.167E+00, 0.168E+00, 0.149E+00, 0.160E+00, &
        0.160E+00, 0.160E+00, 0.160E+00, 0.160E+00 /

data (oasm(jj),jj=1,19)/   &
        0.110E+00, 0.110E+00, 0.110E+00, 0.097E+00, 0.089E+00, &
        0.076E+00, 0.068E+00, 0.063E+00, 0.060E+00, 0.060E+00, &
        0.060E+00, 0.063E+00, 0.068E+00, 0.076E+00, 0.089E+00, &
        0.097E+00, 0.110E+00, 0.110E+00, 0.110E+00 /

real      :: albgla=0.7  ! albedo of glacier points without snowcover 
logical   :: skyhi_res=.false.
integer   :: x(4), y(4), id, jd, jdf, idf
logical   :: do_climap_albedo, do_zonal_albedo 



!---------------------------------------------------------------------
!---------------------------------------------------------------------




                         contains





subroutine surface_albedo_init
 
!---------------------------------------------------------------------
      integer    ::  unit, ierr, io, iounit
      integer    ::  j, ll, i
      real       ::  fl, than
!---------------------------------------------------------------------

      real, dimension(:,:), allocatable  :: alb_n18, alb_n30, alb_n45, &
					    alb_n90
 
!---------------------------------------------------------------------
!-----  read namelist  ------

      if (file_exist('input.nml')) then
        unit =  open_file ('input.nml', action='read')
        ierr=1; do while (ierr /= 0)
        read (unit, nml=surface_albedo_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'surface_albedo_nml')
        enddo
 10     call close_file (unit)
      endif

      unit = open_file ('logfile.out', action='append')
!     call print_version_number (unit, 'surface_albedo', version_number)
      if (get_my_pe() == 0) then
	write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
	write (unit,nml=surface_albedo_nml)
      endif
      call close_file (unit)

      if (trim(albedo_type) == 'climap') then
	do_climap_albedo = .true.
	do_zonal_albedo  = .false.
      else if (trim(albedo_type) == 'zonal') then
	do_climap_albedo = .false.
	do_zonal_albedo  = .true.
      else
	call error_mesg ('surface_albedo_input', &
	  ' improper specification for surface albedo type', FATAL)
      endif

      call get_domain_decomp (x, y)
      idf = x(4) -x(3) + 1
      jdf = y(4) -y(3) + 1
      id  = x(2) -x(1) + 1
      jd  = y(2) -y(1) + 1

!---------------------------------------------------------------------
!    determine if skyhi albedo formulation is usable in this model 
!    configuration.
!--------------------------------------------------------------------
      if (Environment%running_gcm) then
        if ((JD == 36  .and. ID == 60)   .or.  &
	    (JD == 60  .and. ID == 100)  .or.  &
	    (JD == 90  .and. ID == 150)  .or.  &
	    (JD == 180 .and. ID == 300)) then
          skyhi_res = .true.
        else
          skyhi_res = .false.
        endif

!---------------------------------------------------------------------
!     interpolate the surface infra-red reflectivities to model latit-
!     udes. arrays zoal and zoas are defined at 90, 80, 70,... degrees
!     of latitude, from south pole to north pole.
!----------------------------------------------------------------------
        if (skyhi_res) then
          if (do_zonal_albedo) then
	    allocate (zoal (jdf) )
	    allocate (zoas (jdf) )
	    allocate (zoal_gl (jd) )
	    allocate (zoas_gl (jd) )
            do j=1,JD/2
	      fl     = 18.0/float(JD)*float(j) - 9.0/float(jd)
              ll     = fl
              than   = fl - FLOAT(ll)
              ll     = ll + 1
              zoal_gl(j) = oalm(ll) + than*(oalm(ll+1) - oalm(ll))
              zoal_gl(JD+1-j) = zoal_gl(j)
              zoas_gl(j) = oasm(ll) + than*(oasm(ll+1) - oasm(ll))
              zoas_gl(JD+1-j) = zoas_gl(j)
            end do
	    do j=1,jdf
	      zoal(j) = zoal_gl(j+y(3)-1)
	      zoas(j) = zoas_gl(j+y(3)-1)
            end do
	    deallocate (zoal_gl)
	    deallocate (zoas_gl)

!----------------------------------------------------------------------c
! initialize the surface albedo with the CLIMAP land albedo data       c
! interpolated to the SKYHI grid.                                      c
! note: the surface albedo for non-land (ocean or glacial) grid points c
!       are assigned the values for adjacent land points. this is done c
!       to avoid any mismatch between what the CLIMAP data and the     c
!       model consider a land point near the continental boundaries.   c
!----------------------------------------------------------------------c
          else if (do_climap_albedo) then
	    iounit = open_file ('INPUT/sfcalb', action='read', &
		   	        form='unformatted')
!    if (get_my_pe() == 0) then
!    iounit2 = open_file ('INPUT/sfcalb_o3k', action='append', &
!	   	        form='unformatted')
!    endif
	    allocate ( alb_n18 (60,  36 ) )
	    allocate ( alb_n30 (100, 60 ) )
	    allocate ( alb_n45 (150, 90 ) )
	    allocate ( alb_n90 (300, 180) )
            read (iounit) alb_n18
            read (iounit) alb_n30
            read (iounit) alb_n45
            read (iounit) alb_n90
!    if (get_my_pe() == 0) then
!           write(iounit2) alb_n18
!           write(iounit2) alb_n30
!           write(iounit2) alb_n45
!           write(iounit2) alb_n90
!    endif
            call close_file (iounit)
!    if (get_my_pe() == 0) then
!           call close_file (iounit2)
!    endif
 
!---------------------------------------------------------------------
!   allocate the proper size array and store the proper data in it.
!--------------------------------------------------------------------
            if (JD == 36) then
              allocate (zoalnew_gl(1:60, 1:36) )
              zoalnew_gl(:,:) = alb_n18(:,:)
            else if (JD == 60) then
              allocate (zoalnew_gl(1:100, 1:60) )
              zoalnew_gl(:,:) = alb_n30(:,:)
            else if (JD == 90) then
              allocate (zoalnew_gl(1:150, 1:90) )
              zoalnew_gl(:,:) = alb_n45(:,:)
            else if (JD == 180) then
              allocate (zoalnew_gl(1:300, 1:180) )
              zoalnew_gl(:,:) = alb_n90(:,:)
            else
	      call error_mesg ('surface_albedo_init', &
               'no climap albedo data for the specified resolution', &
							   FATAL)
            endif
            allocate (zoalnew(idf, jdf) )
	    do j=1,jdf
              do i=1,idf
		zoalnew(i,j) = zoalnew_gl(i+x(3)-1, j+y(3)-1)
              end do
            end do
	    deallocate (zoalnew_gl)
!---------------------------------------------------------------------
!   to avoid the (rare) event of a land point (re SKYHI) having a
!   zero land albedo (which may happen in the Arctic) and then
!   a zero surface albedo (when snowmelt occurs) all grid points are
!   assigned a land albedo of 0.15. This value will be overridden for
!   all sea, sea ice or land glacier points and all land points 
!   specified according to CLIMAP.
!---------------------------------------------------------------------
            do j=1,jdf
	      do i=1,idf
	        if (zoalnew(i,j) .LT. 1.0E-04) then
	          zoalnew(i,j) = 0.15
	        endif
	      enddo
            enddo

	    deallocate (alb_n18)
	    deallocate (alb_n30)
	    deallocate (alb_n45)
	    deallocate (alb_n90)
          endif ! (do_zonal_albedo)
        endif   ! (skyhi_res)
      endif
!---------------------------------------------------------------------


end  subroutine surface_albedo_init



!###################################################################

subroutine get_surfacealbedo_for_swrad (cirabgd_out, cirrfgd_out,  &
                                        cvisrfgd_out)

!-------------------------------------------------------------------
real, dimension(:,:),intent(out) :: cirabgd_out, cirrfgd_out,  &
				    cvisrfgd_out
!--------------------------------------------------------------------

      integer  :: is, ie, js, je

      is = lbound(cirabgd_out,1)
      ie = ubound(cirabgd_out,1)
      js = lbound(cirabgd_out,2)
      je = ubound(cirabgd_out,2)
      
      cirabgd_out (:,:) = cirabgd (is:ie, js:je)
      cirrfgd_out (:,:) = cirrfgd (is:ie, js:je)
      cvisrfgd_out(:,:) = cvisrfgd(is:ie, js:je)

end  subroutine get_surfacealbedo_for_swrad


!##################################################################

subroutine sfcalbedo

!-------------------------------------------------------------------
      real, dimension(:,:), allocatable :: snowsq, oal, oas
      real, dimension(:,:), allocatable :: zenith
      integer                           :: j, i

!-----------------------------------------------------------------------
!     these are the assumed reflectivities for seaice, seaice plus snow,
!     glacier, and glacier plus snow. 
!-----------------------------------------------------------------------
      real, dimension(4)                :: cafx   
      data cafx  /0.70E+00, 0.35E+00, 0.70E+00, 0.60E+00/

      if (Environment%running_gcm) then
        if (Environment%running_skyhi) then

!---------------------------------------------------------------------
!     allocate local arrays.
!---------------------------------------------------------------------
          allocate (zenith    (IMINP:IMAXP, JMINP:JMAXP)) 
          allocate (snowsq    (IMINP:IMAXP, JMINP:JMAXP)) 
          allocate (oal       (IMINP:IMAXP, JMINP:JMAXP)) 
          allocate (oas       (IMINP:IMAXP, JMINP:JMAXP)) 

!---------------------------------------------------------------------
!     define 2d surface albedo arrays, either from zonal values or 
!     climap data.
!---------------------------------------------------------------------
          if (do_zonal_albedo) then
            do j=JMINP,JMAXP
	      do i=IMINP,IMAXP
                oal(i,j) = zoal(jabs (j))
                oas(i,j) = zoas(jabs (j))
              end do
            end do

!---------------------------------------------------------------------
!    obtain the zenith angle from the astronomy module.
!--------------------------------------------------------------------
          else if (do_climap_albedo) then
            call get_astronomy_for_sfcalbedo (zenith)

!-----------------------------------------------------------------------
!     define new land and ocean albedos, ocean values use the Edwards-
!     Slingo parameterization.
!-----------------------------------------------------------------------
            do j=JMINP,JMAXP
              do i=IMINP,IMAXP
                oal(i,j) = zoalnew(iabs (i),jabs (j))
                oas(i,j) = 3.7E-02/(1.1E+00*zenith(i,j)**1.4E+00 +  &
			   1.5E-01)
              end do
            end do
          endif

!-----------------------------------------------------------------------
!     define the reflectivities and absorptivities at the ground.
!-----------------------------------------------------------------------
          do j=JMINP,JMAXP
            do i=IMINP,IMAXP
              if (snow(i,j     ) .GT. 1.0E-15) then
                snowsq(i,j)  = SQRT(snow(i,j     ))
              else
                snowsq(i,j)  = 0.0E+00
              endif

              if (sealp(i,j) .EQ. 1.0E+00) then
                cirrfgd(i,j) = oas(i,j)
              else if(sealp(i,j) .EQ. 0.0E+00) then
                if (snow(i,j     ) .LT. 1.0E+00) then
                  cirrfgd(i,j) = oal(i,j) + snowsq(i,j)*  &
				 (cafx(3) - oal(i,j))
                else
                  cirrfgd(i,j) = cafx(3)
                endif
              else if (sealp(i,j) .EQ. 2.0E+00) then
                if (snow(i,j     ) .LT. 1.0E+00) then
                  cirrfgd(i,j) = albgla + snowsq(i,j)*(cafx(3) - albgla)
                else
                  cirrfgd(i,j) = cafx(3)
                endif
              else
                cirrfgd(i,j) = cafx(1)
              endif
              cirabgd (i,j) = 0.0E+00
              cvisrfgd(i,j) = cirrfgd(i,j)
            end do
          end do
        else ! (not running skyhi)
	  cirabgd(:,:) = 0.0
	  cirrfgd(:,:) = albedo_sv(:,:)
	  cvisrfgd(:,:) = albedo_sv(:,:)
        endif
      else if (Environment%running_standalone) then 
        do j=JMINP,JMAXP
          do i=IMINP,IMAXP
            cirrfgd(i,j) = single_column_albedo
          enddo
        enddo
        cirabgd(:,:) = 0.0E+00
        cvisrfgd(:,:) = cirrfgd(:,:)
      endif

!------------------------------------------------------------------
!    send the surface albedo components to the radiation_diag_mod for
!    processing.
!------------------------------------------------------------------
      if (Environment%running_gcm) then
        do j=JMINP,JMAXP
          if (Rad_control%do_raddg(jabs (j))) then
            if (Environment%running_skyhi) then
              call radiag_from_sfcalbedo (j, cvisrfgd, cirrfgd)       
	    else 
              call radiag_from_sfcalbedo (j, albedo_sv, albedo_sv)   
	    endif
          endif
        end do
      else if (Environment%running_standalone) then
        do j=JMINP,JMAXP
          if (Rad_control%do_raddg(jabs (j))) then
            call radiag_from_sfcalbedo (j, cvisrfgd, cirrfgd)           
          endif
        end do
      endif

!-------------------------------------------------------------------
!    send the surface albedo components to a routine where the value 
!    of albp may be placed in a model common block for use elsewhere 
!    in the model.
!-------------------------------------------------------------------
      if (Environment%running_gcm) then
        if (Environment%running_skyhi) then
          call store_sfcalbedo_wrap (cirrfgd)
        endif
      endif

!-------------------------------------------------------------------
!    deallocate the local arrays.
!------------------------------------------------------------------
      if (Environment%running_skyhi .and. Environment%running_gcm) then
        deallocate (zenith)
        deallocate (snowsq)
        deallocate (oal   )
        deallocate (oas   )
      endif

!-------------------------------------------------------------------


end  subroutine sfcalbedo


!#####################################################################

subroutine sfcalbedo_alloc

      allocate (cirabgd   (IMINP:IMAXP, JMINP:JMAXP)) 
      allocate (cirrfgd   (IMINP:IMAXP, JMINP:JMAXP)) 
      allocate (cvisrfgd  (IMINP:IMAXP, JMINP:JMAXP)) 

end subroutine sfcalbedo_alloc


!#####################################################################

subroutine sfcalbedo_dealloc

      deallocate (cvisrfgd ) 
      deallocate (cirrfgd ) 
      deallocate (cirabgd ) 

end subroutine sfcalbedo_dealloc


!###################################################################


               end module surface_albedo_mod
