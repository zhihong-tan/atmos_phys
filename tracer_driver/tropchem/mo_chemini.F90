
      module mo_chemini_mod

implicit none
      private
      public :: chemini

character(len=128), parameter :: version     = '$Id: mo_chemini.F90,v 13.0 2006/03/28 21:16:05 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: memphis_2006_12 $'
logical                       :: module_is_initialized = .false.

      contains

      subroutine chemini( calday )
!-----------------------------------------------------------------------
! 	... Chemistry module intialization
!-----------------------------------------------------------------------

      use MO_PHOTO_MOD,      only : PRATE_INIT
      use mo_chem_utls_mod,  only : chem_utls_init
      use mo_usrrxt_mod,     only : usrrxt_init
      use CHEM_MODS_mod,     only : grpcnt, clscnt1, clscnt4, clscnt5, CHEM_MODS_INIT
      use MO_EXP_SOL_mod,    only : EXP_SLV_INIT
      use MO_IMP_SOL_mod,    only : IMP_SLV_INIT
      use MO_RODAS_SOL_mod,  only : RODAS_SLV_INIT

      use MO_READ_SIM_CHM_mod, only : read_sim_chm

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      real, intent(in)    ::  calday

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      character(len=80) ::   lpath
      character(len=80) ::   mspath
      character(len=32) ::   filename
      character(len=30) ::   emires, surfres
      
      character(len=128) ::  sim
      integer :: sim_file_cnt

!-----------------------------------------------------------------------
! 	... Allocate variables
!-----------------------------------------------------------------------
      call chem_mods_init

!-----------------------------------------------------------------------
! 	... Read sim.dat
!-----------------------------------------------------------------------
      sim = 'INPUT/sim.dat'
      call read_sim_chm( sim, sim_file_cnt )

!-----------------------------------------------------------------------
! 	... Diagnostics initialization
!-----------------------------------------------------------------------
!     call diags_init( tracnam, plonl, platl, pplon )

!-----------------------------------------------------------------------
! 	... Initialize photorate module
!-----------------------------------------------------------------------
!     filename = photo_flsp%nl_filename
!     lpath    = photo_flsp%local_path
!     mspath   = photo_flsp%remote_path
      filename = 'INPUT/jvals.v5'
      lpath = ''
      mspath = ''
      call prate_init( filename, lpath, mspath )

!-----------------------------------------------------------------------
! 	... Read time-independent airplane emissions
!-----------------------------------------------------------------------
!     emires = emis_flsp%hor_res
!     if( emires(1:1) /= '.' ) then
!        emires = '.' // emires
!     end if
!     lpath    = emis_flsp%local_path
!     mspath   = emis_flsp%remote_path
!     filename = 'emissions.aircraft' // TRIM(emires) // '.nc'
!     call airpl_src( filename, lpath, mspath, plonl, platl, pplon )

!-----------------------------------------------------------------------
! 	... Initialize the chem utils module
!-----------------------------------------------------------------------
      call chem_utls_init

!-----------------------------------------------------------------------
! 	... Read time-dependent surface flux dataset
!-----------------------------------------------------------------------
!     call srf_emis_init( plonl, platl, pplon )

!-----------------------------------------------------------------------
! 	... Intialize the het rates module
!-----------------------------------------------------------------------
!     call sethet_init

!-----------------------------------------------------------------------
! 	... Intialize the ext frcing module
!-----------------------------------------------------------------------
!     call setext_init

!-----------------------------------------------------------------------
! 	... Intialize the rxt rate constant module
!-----------------------------------------------------------------------
      call usrrxt_init

!-----------------------------------------------------------------------
! 	... Intialize the grp ratios module
!-----------------------------------------------------------------------
!     call set_grp_ratios_init

!-----------------------------------------------------------------------
! 	... Read time-dependent surface variables dataset
!-----------------------------------------------------------------------
!     surfres = surf_flsp%hor_res
!     if( surfres(1:1) /= '.' ) then
!        surfres = '.' // surfres
!     end if
!     filename = 'surfdata' // TRIM(surfres) // '.nc'
!     lpath    = surf_flsp%local_path
!     mspath   = surf_flsp%remote_path
!     call surf_init( filename, lpath, mspath, plonl, platl, pplon )

!-----------------------------------------------------------------------
! 	... Read time-dependent upper boundary values
!-----------------------------------------------------------------------
!     filename = ubc_flsp%nl_filename
!     lpath    = ubc_flsp%local_path
!     mspath   = ubc_flsp%remote_path
!     call ub_init( platl, filename, lpath, mspath )

!-----------------------------------------------------------------------
! 	... Read time-dependent sulfate dataset
!	    NOTE : This is now a netcdf dataset
!-----------------------------------------------------------------------
!     filename = 'sulfate.M1.nc'
!     lpath    = sulf_flsp%local_path
!     mspath   = sulf_flsp%remote_path
!     call sulf_init( plonl, platl, pplon, filename, lpath, mspath )

      if( clscnt1 > 0 ) then
!-----------------------------------------------------------------------
!	... Initialize the explicit solver
!-----------------------------------------------------------------------
         call exp_slv_init
      end if
      if( clscnt4 > 0 ) then
!-----------------------------------------------------------------------
!	... Initialize the implicit solver
!-----------------------------------------------------------------------
         call imp_slv_init
      end if
      if( clscnt5 > 0 ) then
!-----------------------------------------------------------------------
!	... Initialize the implicit solver
!-----------------------------------------------------------------------
         call rodas_slv_init
      end if

      end subroutine chemini

      end module mo_chemini_mod
