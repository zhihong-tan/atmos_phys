
       module donner_ice_mod

!----------------------------------------------------------------
implicit none
private

!---------------------------------------------------------------
!         this is a stub module for use with SKYHI -- it will
!         be filled in when this module code is made available
!
!---------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!    character(len=5), parameter  ::  version_number = 'v0.08'
     character(len=128)  :: version =  '$Id: donner_ice.F90,v 1.2 2001/08/30 15:13:44 fms Exp $'
     character(len=128)  :: tag     =  '$Name: eugene $'

!--------------------------------------------------------------------
!-------  interfaces --------

public  inquire_donner_ice, sw_albedo_zen_angle, &
	donner_iceclouds_rad_output


!---------------------------------------------------------------------
!------- private data ------

!  this will ultimately be a namelist variable
logical      ::   do_donner_ice=.false.

!-------------------------------------------------------------------
!-------------------------------------------------------------------



       contains



subroutine inquire_donner_ice (do_donner_ice_out)

logical, intent(out) :: do_donner_ice_out

       do_donner_ice_out = do_donner_ice

end subroutine inquire_donner_ice



!################################################################

subroutine sw_albedo_zen_angle (ncldsw, cvis, cirr, nsolwg)

integer, dimension(:,:), intent(in)   :: ncldsw
real, dimension(:,:,:,:),intent(in)   :: cvis, cirr
integer,                 intent(in)   :: nsolwg

    return


end subroutine sw_albedo_zen_angle 
  


!##################################################################

subroutine donner_iceclouds_rad_output

     return


end subroutine donner_iceclouds_rad_output




!##################################################################



          end module donner_ice_mod
