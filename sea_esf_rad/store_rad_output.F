
             module store_rad_output_mod

!----------------------------------------------------------------
       implicit none
       private

!---------------------------------------------------------------
!      this is a stub module for use when running FMS or the stand-
!      alone code -- it provides dummy interfaces in these applications
!      to take the place of actual interfaces used in SKYHI
!
!---------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!       character(len=5), parameter  ::  version_number = 'v0.08'
character(len=128)  :: version =  '$Id: store_rad_output.F,v 1.2 2001/08/30 15:16:03 fms Exp $'
character(len=128)  :: tag     =  '$Name: galway $'

!--------------------------------------------------------------------
!-------  interfaces --------

        public  store_sfcalbedo_wrap,  store_radiation_output 


!-------------------------------------------------------------------
!-------------------------------------------------------------------



                          contains


       subroutine store_sfcalbedo_wrap (cirrfgd)

          real, dimension(:,:), intent(in) :: cirrfgd


       end subroutine store_sfcalbedo_wrap



!################################################################

       subroutine store_radiation_output


       end subroutine store_radiation_output



!####################################################################




             end module store_rad_output_mod

