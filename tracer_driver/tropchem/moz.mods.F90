      module mo_grid_mod
!---------------------------------------------------------------------
! 	... Basic grid point resolution parameters
!---------------------------------------------------------------------
      implicit none

      save

character(len=128), parameter :: version     = '$Id: moz.mods.F90,v 13.0 2006/03/28 21:16:35 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: memphis_2006_08 $'
logical                       :: module_is_initialized = .false.

      integer, parameter :: &
                pcnst    =    41+1, &     ! number of advected constituents including cloud water
                pcnstm1  =    41, &     ! number of advected constituents excluding cloud water
                plev     =   1, &         ! number of vertical levels
                plevp    = plev+1, &      ! plev plus 1
                plevm    = plev-1, &      ! plev minus 1
                plon     =   1, &         ! number of longitudes
                plat     =   1            ! number of latitudes

      integer, parameter :: &
                pnats    =     0    ! number of non-advected trace species



      integer :: nodes                ! mpi task count
      integer :: plonl                ! longitude tile dimension
      integer :: pplon                ! longitude tile count
      integer :: plnplv               ! plonl * plev

      end module mo_grid_mod

      module chem_mods_mod
!--------------------------------------------------------------
!     	... basic chemistry array parameters
!--------------------------------------------------------------

      use mo_grid_mod, only : pcnstm1
!++lwh
      use mpp_mod,     only : mpp_error, FATAL
!--lwh

      implicit none

      save

      integer, parameter :: hetcnt     =     0, &    ! number of heterogeneous processes
                            phtcnt     =    20, &    ! number of photo processes
                            rxntot     =   100, &    ! number of total reactions
                            gascnt     =    80, &    ! number of gas phase reactions
                            nfs        =     8, &       ! number of "fixed" species
                            relcnt     =     0, &    ! number of relationship species
                            grpcnt     =     0, &    ! number of group members
                            imp_nzcnt  =   277, &     ! number of non-zero implicit matrix entries
                            rod_nzcnt  =     0, &     ! number of non-zero rodas matrix entries
                            extcnt     =     0, &    ! number of species with external forcing
                            clscnt1    =     5, &  ! number of species in explicit class
                            clscnt2    =     0, &  ! number of species in hov class
                            clscnt3    =     0, &  ! number of species in ebi class
                            clscnt4    =    35, &  ! number of species in implicit class
                            clscnt5    =     0, &  ! number of species in rodas class
                            indexm     =     1, &    ! index of total atm density in invariant array
                            ncol_abs   =     2, &    ! number of column densities
                            indexh2o   =     4, &    ! index of water vapor density
                            clsze      =   1       ! loop length for implicit chemistry

      integer ::            ngrp       = 0
      integer ::            drydep_cnt = 0
      integer ::            srfems_cnt = 0
      integer ::            rxt_alias_cnt = 0
      integer, allocatable :: grp_mem_cnt(:)
      integer, allocatable :: rxt_alias_map(:)
      real      :: adv_mass(pcnstm1)
      real      :: nadv_mass(grpcnt)
      character(len=16), allocatable :: rxt_alias_lst(:)
      character(len=8), allocatable  :: drydep_lst(:)
      character(len=8), allocatable  :: srfems_lst(:)
      character(len=8), allocatable  :: grp_lst(:)
      character(len=8)               :: het_lst(max(1,hetcnt))
      character(len=8)               :: extfrc_lst(max(1,extcnt))

      type solver_class
	 integer :: clscnt
	 integer :: lin_rxt_cnt
	 integer :: nln_rxt_cnt
	 integer :: indprd_cnt
	 integer :: iter_max
         integer :: cls_rxt_cnt(4)
         integer, pointer :: permute(:)
         integer, pointer :: diag_map(:)
         integer, pointer :: clsmap(:)
      end type solver_class

      type(solver_class) :: explicit, implicit, rodas

      contains

      subroutine endrun(msg)

      implicit none

      character(len=128), intent(in), optional  :: msg
      call mpp_error(FATAL, msg)

      end subroutine endrun 

      subroutine chem_mods_init
!--------------------------------------------------------------
!     	... intialize the class derived type
!--------------------------------------------------------------

      implicit none

      integer :: astat

      explicit%clscnt       =     5
      explicit%indprd_cnt   =    10

      implicit%clscnt       =    35
      implicit%lin_rxt_cnt  =    37
      implicit%nln_rxt_cnt  =    60
      implicit%indprd_cnt   =     1
      implicit%iter_max     =    11

      rodas%clscnt          =     0
      rodas%lin_rxt_cnt     =     0
      rodas%nln_rxt_cnt     =     0
      rodas%indprd_cnt      =     0

      if( explicit%clscnt > 0 ) then
	 allocate( explicit%clsmap(explicit%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_init: failed to allocate explicit%clsmap ; error = ',astat
	    call endrun
	 end if
         explicit%clsmap(:)  = 0
      end if
      if( implicit%clscnt > 0 ) then
	 allocate( implicit%permute(implicit%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_init: failed to allocate implicit%permute ; error = ',astat
	    call endrun
	 end if
         implicit%permute(:)  = 0
	 allocate( implicit%diag_map(implicit%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_init: failed to allocate implicit%diag_map ; error = ',astat
	    call endrun
	 end if
         implicit%diag_map(:)  = 0
	 allocate( implicit%clsmap(implicit%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_init: failed to allocate implicit%clsmap ; error = ',astat
	    call endrun
	 end if
         implicit%clsmap(:)  = 0
      end if
      if( rodas%clscnt > 0 ) then
	 allocate( rodas%permute(rodas%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_init: failed to allocate rodas%permute ; error = ',astat
	    call endrun
	 end if
         rodas%permute(:)  = 0
	 allocate( rodas%diag_map(rodas%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_init: failed to allocate rodas%diag_map ; error = ',astat
	    call endrun
	 end if
         rodas%diag_map(:)  = 0
	 allocate( rodas%clsmap(rodas%clscnt),stat=astat )
	 if( astat /= 0 ) then
	    write(*,*) 'chem_mods_init: failed to allocate rodas%clsmap ; error = ',astat
	    call endrun
	 end if
         rodas%clsmap(:)  = 0
      end if

      end subroutine chem_mods_init

      end module chem_mods_mod
  
      module M_SPC_ID_MOD
  
      implicit none                                                             
  
      integer, parameter :: id_O3 =   1
      integer, parameter :: id_O =   2
      integer, parameter :: id_O1D =   3
      integer, parameter :: id_NO =   4
      integer, parameter :: id_NO2 =   5
      integer, parameter :: id_NO3 =   6
      integer, parameter :: id_HNO3 =   7
      integer, parameter :: id_HO2NO2 =   8
      integer, parameter :: id_N2O5 =   9
      integer, parameter :: id_CH3O2 =  10
      integer, parameter :: id_CH3OOH =  11
      integer, parameter :: id_CH2O =  12
      integer, parameter :: id_CO =  13
      integer, parameter :: id_OH =  14
      integer, parameter :: id_HO2 =  15
      integer, parameter :: id_H2O2 =  16
      integer, parameter :: id_C3H6 =  17
      integer, parameter :: id_ISOP =  18
      integer, parameter :: id_CH3CHO =  19
      integer, parameter :: id_CH3CO3 =  20
      integer, parameter :: id_CH3COOOH =  21
      integer, parameter :: id_PAN =  22
      integer, parameter :: id_ISOPO2 =  23
      integer, parameter :: id_C2H5O2 =  24
      integer, parameter :: id_C2H5OOH =  25
      integer, parameter :: id_C3H7O2 =  26
      integer, parameter :: id_C3H7OOH =  27
      integer, parameter :: id_CH3COCH3 =  28
      integer, parameter :: id_ROOH =  29
      integer, parameter :: id_RO2 =  30
      integer, parameter :: id_CH3COCHO =  31
      integer, parameter :: id_BC1 =  32
      integer, parameter :: id_BC2 =  33
      integer, parameter :: id_OC1 =  34
      integer, parameter :: id_OC2 =  35
      integer, parameter :: id_SO2 =  36
      integer, parameter :: id_SO4 =  37
      integer, parameter :: id_DMS =  38
      integer, parameter :: id_NH3 =  39
      integer, parameter :: id_NH4NO3 =  40
      integer, parameter :: id_NH4 =  41
  
  
      end module M_SPC_ID_MOD
                                                                                
      module M_RXT_ID_MOD
                                                                                
      implicit none                                                             
                                                                                
      integer, parameter :: rid_jo2 =    1                                      
      integer, parameter :: rid_jo1d =    2                                     
      integer, parameter :: rid_jo3p =    3                                     
      integer, parameter :: rid_jno2 =    4                                     
      integer, parameter :: rid_jn2o5 =    5                                    
      integer, parameter :: rid_jhno3 =    6                                    
      integer, parameter :: rid_jno3 =    7                                     
      integer, parameter :: rid_jho2no2 =    8                                  
      integer, parameter :: rid_jch3ooh =    9                                  
      integer, parameter :: rid_jch2o_a =   10                                  
      integer, parameter :: rid_jch2o_b =   11                                  
      integer, parameter :: rid_jh2o2 =   12                                    
      integer, parameter :: rid_jch3cho =   13                                  
      integer, parameter :: rid_jch3co3h =   14                                 
      integer, parameter :: rid_jpan =   15                                     
      integer, parameter :: rid_jc2h5ooh =   16                                 
      integer, parameter :: rid_jc3h7ooh =   17                                 
      integer, parameter :: rid_jrooh =   18                                    
      integer, parameter :: rid_jacet =   19                                    
      integer, parameter :: rid_jmgly =   20                                    
      integer, parameter :: rid_usr1 =   21                                     
      integer, parameter :: rid_o1d_n2 =   23                                   
      integer, parameter :: rid_o1d_o2 =   24                                   
      integer, parameter :: rid_ox_l1 =   25                                    
      integer, parameter :: rid_ox_p1 =   26                                    
      integer, parameter :: rid_usr2 =   31                                     
      integer, parameter :: rid_usr3 =   32                                     
      integer, parameter :: rid_usr4 =   34                                     
      integer, parameter :: rid_usr5 =   35                                     
      integer, parameter :: rid_usr6 =   37                                     
      integer, parameter :: rid_usr7 =   39                                     
      integer, parameter :: rid_ox_p2 =   42                                    
      integer, parameter :: rid_usr8 =   49                                     
      integer, parameter :: rid_ox_l2 =   53                                    
      integer, parameter :: rid_ox_l3 =   54                                    
      integer, parameter :: rid_usr9 =   55                                     
      integer, parameter :: rid_ox_p4 =   62                                    
      integer, parameter :: rid_usr11 =   63                                    
      integer, parameter :: rid_usr12 =   67                                    
      integer, parameter :: rid_ox_l5 =   69                                    
      integer, parameter :: rid_ox_p5 =   71                                    
      integer, parameter :: rid_usr16 =   79                                    
      integer, parameter :: rid_usr17 =   80                                    
      integer, parameter :: rid_ox_p9 =   82                                    
      integer, parameter :: rid_usr22 =   86                                    
      integer, parameter :: rid_ox_p10 =   87                                   
      integer, parameter :: rid_usr24 =   97                                    
      integer, parameter :: rid_usr25 =   99                                    
                                                                                
      integer, parameter :: rid_r0022 =   22                                    
      integer, parameter :: rid_r0027 =   27                                    
      integer, parameter :: rid_r0028 =   28                                    
      integer, parameter :: rid_r0029 =   29                                    
      integer, parameter :: rid_r0030 =   30                                    
      integer, parameter :: rid_r0033 =   33                                    
      integer, parameter :: rid_r0036 =   36                                    
      integer, parameter :: rid_r0038 =   38                                    
      integer, parameter :: rid_r0040 =   40                                    
      integer, parameter :: rid_r0041 =   41                                    
      integer, parameter :: rid_r0043 =   43                                    
      integer, parameter :: rid_r0044 =   44                                    
      integer, parameter :: rid_r0045 =   45                                    
      integer, parameter :: rid_r0046 =   46                                    
      integer, parameter :: rid_r0047 =   47                                    
      integer, parameter :: rid_r0048 =   48                                    
      integer, parameter :: rid_r0050 =   50                                    
      integer, parameter :: rid_r0051 =   51                                    
      integer, parameter :: rid_r0052 =   52                                    
      integer, parameter :: rid_r0056 =   56                                    
      integer, parameter :: rid_r0057 =   57                                    
      integer, parameter :: rid_r0058 =   58                                    
      integer, parameter :: rid_r0059 =   59                                    
      integer, parameter :: rid_r0060 =   60                                    
      integer, parameter :: rid_r0061 =   61                                    
      integer, parameter :: rid_r0064 =   64                                    
      integer, parameter :: rid_r0065 =   65                                    
      integer, parameter :: rid_r0066 =   66                                    
      integer, parameter :: rid_r0068 =   68                                    
      integer, parameter :: rid_r0070 =   70                                    
      integer, parameter :: rid_r0072 =   72                                    
      integer, parameter :: rid_r0073 =   73                                    
      integer, parameter :: rid_r0074 =   74                                    
      integer, parameter :: rid_r0075 =   75                                    
      integer, parameter :: rid_r0076 =   76                                    
      integer, parameter :: rid_r0077 =   77                                    
      integer, parameter :: rid_r0078 =   78                                    
      integer, parameter :: rid_r0081 =   81                                    
      integer, parameter :: rid_r0083 =   83                                    
      integer, parameter :: rid_r0084 =   84                                    
      integer, parameter :: rid_r0085 =   85                                    
      integer, parameter :: rid_r0088 =   88                                    
      integer, parameter :: rid_r0089 =   89                                    
      integer, parameter :: rid_r0090 =   90                                    
      integer, parameter :: rid_r0091 =   91                                    
      integer, parameter :: rid_r0092 =   92                                    
      integer, parameter :: rid_r0093 =   93                                    
      integer, parameter :: rid_r0094 =   94                                    
      integer, parameter :: rid_r0095 =   95                                    
      integer, parameter :: rid_r0096 =   96                                    
      integer, parameter :: rid_r0098 =   98                                    
      integer, parameter :: rid_r0100 =  100                                    
                                                                                
      end module M_RXT_ID_MOD
                                                                                
      module M_HET_ID_MOD
                                                                                
      implicit none                                                             
                                                                                
                                                                                
      end module M_HET_ID_MOD
