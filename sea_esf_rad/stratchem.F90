
                        module stratchem_mod

implicit none
private

!------------------------------------------------------------------
!            temporary module used to supply needed hooks
!                  until actual module is created   
!
!------------------------------------------------------------------



!-------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

!   character(len=5), parameter  ::  version_number = 'v0.08'
    character(len=128)  :: version =  '$Id: stratchem.F90,v 1.2 2001/08/30 15:12:30 fms Exp $'
    character(len=128)  :: tag     =  '$Name: eugene $'

!---------------------------------------------------------------------
!-------  interfaces --------

public  inquire_stratchem

!---------------------------------------------------------------------
!------- public data -----

logical :: do_stratchem=.false.  !  will be namelist variable once  
				 !  mod  is developed


!---------------------------------------------------------------------
!---------------------------------------------------------------------


contains



subroutine inquire_stratchem (do_stratchem_out)

logical, intent(out)     ::   do_stratchem_out

      do_stratchem_out = do_stratchem

end subroutine inquire_stratchem



!#####################################################################


                  end module stratchem_mod
