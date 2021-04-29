!$Id: rtm_common_mod.f90 4052 2020-11-19 13:55:22Z awalther $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: rtm_common.f90 (src)
!       RTM_COMMON_MOD (program)
!
! PURPOSE: RADIATIVE TRANSFER MODEL COMMON MEMORY MODULE
!
! DESCRIPTION: RTM_COMMON_MOD stores the RTM data structure makes available to other
!              modules
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! Dependencies:  (Following are names of other CLAVR-x modules)
!  CONSTANTS_MOD
! 
!   Members of RTM structure
!   Flag = flag to indicate where this RTM cell is populated (0=no,1=yes)
!   Sfc_Level = level in the rtm profiles the closest to but just above the surface
!   Tropo_Level = level in the rtm pressure prof that is closest to tropopause
!   Inversion_Level = level in rtm pressure that is highest inversion
!   Level440 =  level in rtm pressure profile closest to 440 hPa
!   Level850 =  level in rtm pressure profile closest to 850 hPa
!   type (Chan_Rtm_Prof), dimension(:), ALLOCATABLE :: d
!   T_Prof = Temperature Profile (101 levels) (K)
!   Z_Prof = Geopotential height Profile (km)
!   Wvmr_Prof = Water Vapor Mixing Ratio Profile (101 levels) (g/kg)
!   Tpw_Prof = Total Precipitable Water Profile (101 levels) (g/cm^2 or cm)
!   Ozmr_Prof = Ozone Mixing Ratio Profile (101 levels) (g/kg)

!   Members of RTM Struture that computed for each angle bin (within d sub-structure)
!   Trans_Atm_Profile = Profile of atmospheric transmission from a level to space
!                       along viewing zenith path
!   Trans_Atm_Solar_Profile = Profile of atmospheric transmission from a level to
!                            space along the solar zenith path
!   Trans_Atm_Total_Profile = Profile of atmospheric tranmission from a level to
!                             space for both solar and view zenith paths (aka two-way trans)
!   Rad_Atm_Profile = Profile of atmospheric emitted radiance from level to space
!                     along the viewing zenith path
!   Rad_Atm_Dwn_Profile = Profile of atmospheric emitted radiance from level to surface
!                     along the viewing zenith path
!   Rad_BB_Cloud_Profile = Profile ot atmospheric emitted radiance and black-body cloud
!                          emitted radaince from level to space along the viewing zenith path
!--------------------------------------------------------------------------------------
module RTM_COMMON_MOD
 
   use CONSTANTS_MOD

   implicit none
   private
   
   public ::  ALLOCATE_RTM
   public ::  DEALLOCATE_RTM

   integer, parameter, public:: NLevels_Rtm = 101

   real, dimension(NLevels_Rtm), parameter, public:: P_Std_Rtm = (/ &
          0.0050,    0.0161,    0.0384,    0.0769,    0.1370,  &  
          0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    1.2972,  &
          1.6872,    2.1526,    2.7009,    3.3398,    4.0770,    4.9204,  &
          5.8776,    6.9567,    8.1655,    9.5119,   11.0038,   12.6492,  &
         14.4559,   16.4318,   18.5847,   20.9224,   23.4526,   26.1829,  &
         29.1210,   32.2744,   35.6505,   39.2566,   43.1001,   47.1882,  &
         51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   77.2396,  &
         83.2310,   89.5204,   96.1138,  103.0172,  110.2366,  117.7775,  &
        125.6456,  133.8462,  142.3848,  151.2664,  160.4959,  170.0784,  &
        180.0183,  190.3203,  200.9887,  212.0277,  223.4415,  235.2338,  &
        247.4085,  259.9691,  272.9191,  286.2617,  300.0000,  314.1369,  &
        328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  407.4738,  &
        424.4698,  441.8819,  459.7118,  477.9607,  496.6298,  515.7200,  &
        535.2322,  555.1669,  575.5248,  596.3062,  617.5112,  639.1398,  &
        661.1920,  683.6673,  706.5654,  729.8857,  753.6275,  777.7897,  &
        802.3714,  827.3713,  852.7880,  878.6201,  904.8659,  931.5236,  &
        958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, 1100.0000/)
                                                                                                                                                    
    real, dimension(NLevels_Rtm), parameter, public:: T_Std_Rtm = (/    &
                                190.19, 203.65, 215.30, 226.87, 237.83, &
        247.50, 256.03, 263.48, 267.09, 270.37, 266.42, 261.56, 256.40, & 
        251.69, 247.32, 243.27, 239.56, 236.07, 232.76, 230.67, 228.71, &
        227.35, 226.29, 225.28, 224.41, 223.61, 222.85, 222.12, 221.42, &
        220.73, 220.07, 219.44, 218.82, 218.23, 217.65, 217.18, 216.91, &
        216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, &
        216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.71, &
        216.71, 216.72, 216.81, 217.80, 218.77, 219.72, 220.66, 222.51, &
        224.57, 226.59, 228.58, 230.61, 232.61, 234.57, 236.53, 238.48, &
        240.40, 242.31, 244.21, 246.09, 247.94, 249.78, 251.62, 253.45, &
        255.26, 257.04, 258.80, 260.55, 262.28, 264.02, 265.73, 267.42, &
        269.09, 270.77, 272.43, 274.06, 275.70, 277.32, 278.92, 280.51, &
        282.08, 283.64, 285.20, 286.74, 288.25, 289.75, 291.22, 292.68/) 

    real, dimension(NLevels_Rtm), parameter, public:: Wvmr_Std_Rtm = (/ &
                                 0.001,  0.001,  0.002,  0.003,  0.003, &
         0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, &
         0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, &
         0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, &
         0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, &
         0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, &
         0.003,  0.003,  0.004,  0.004,  0.005,  0.005,  0.007,  0.009, &
         0.011,  0.012,  0.014,  0.020,  0.025,  0.030,  0.035,  0.047, &
         0.061,  0.075,  0.089,  0.126,  0.162,  0.197,  0.235,  0.273, &
         0.310,  0.356,  0.410,  0.471,  0.535,  0.601,  0.684,  0.784, &
         0.886,  0.987,  1.094,  1.225,  1.353,  1.519,  1.686,  1.852, &
         2.036,  2.267,  2.496,  2.721,  2.947,  3.170,  3.391,  3.621, &
         3.848,  4.084,  4.333,  4.579,  4.822,  5.061,  5.296,  5.528/)
                                                                                                                                                    
    real, dimension(NLevels_Rtm), parameter, public:: Ozmr_Std_Rtm = (/ &
                               0.47330,0.27695,0.28678,0.51816,0.83229, &
       1.18466,1.69647,2.16633,3.00338,3.76287,4.75054,5.61330,6.33914, &
       7.03675,7.50525,7.75612,7.81607,7.69626,7.56605,7.28440,7.01002, &
       6.72722,6.44629,6.17714,5.92914,5.69481,5.47387,5.26813,5.01252, &
       4.68941,4.35141,4.01425,3.68771,3.37116,3.06407,2.77294,2.50321, &
       2.24098,1.98592,1.74840,1.54451,1.34582,1.17824,1.02513,0.89358, &
       0.78844,0.69683,0.62654,0.55781,0.50380,0.45515,0.42037,0.38632, &
       0.35297,0.32029,0.28832,0.25756,0.22739,0.19780,0.16877,0.14901, &
       0.13190,0.11511,0.09861,0.08818,0.07793,0.06786,0.06146,0.05768, &
       0.05396,0.05071,0.04803,0.04548,0.04301,0.04081,0.03983,0.03883, &
       0.03783,0.03685,0.03588,0.03491,0.03395,0.03368,0.03349,0.03331, &
       0.03313,0.03292,0.03271,0.03251,0.03190,0.03126,0.03062,0.02990, &
       0.02918,0.02850,0.02785,0.02721,0.02658,0.02596,0.02579,0.02579/)

!---------------------------------------------------------------------
! RTM structure definition
!---------------------------------------------------------------------
  type, public :: Rtm_Prof
    logical :: is_set
    logical :: is_allocated
    real, dimension(:), allocatable :: Trans_Atm_Profile
    real, dimension(:), allocatable :: Trans_Atm_Solar_Profile
    real, dimension(:), allocatable :: Trans_Atm_Total_Profile
    real, dimension(:), allocatable :: Rad_Atm_Profile
    real, dimension(:), allocatable :: Rad_Atm_Dwn_Profile
    real, dimension(:), allocatable :: Rad_BB_Cloud_Profile
  contains
    procedure :: alloc => rtm_prof__alloc
    procedure :: dealloc => rtm_prof__dealloc    
    
  end type Rtm_Prof

  type, public :: Chan_Rtm_Prof
    logical :: is_set = .false.
    logical :: is_allocated = .false.
    type(Rtm_Prof):: ch(45)
  contains
    procedure :: alloc => Chan_Rtm_Prof__alloc
  end type Chan_Rtm_Prof

  type, public :: Rtm_Params
    logical :: is_set = .false.
    logical :: is_allocated = .false.
    integer (kind=int1) :: Sfc_Level = 0
    integer (kind=int1) :: Tropo_Level = 0
    integer (kind=int1) :: Inversion_Level = 0
    integer (kind=int1) :: Level440 = 0
    integer (kind=int1) :: Level850 = 0
    integer (kind=int1) :: Vza_Idx_Min = Missing_Value_Int1
    integer (kind=int1) :: Vza_Idx_Max = Missing_Value_Int1
    type (Chan_Rtm_Prof), dimension(:), ALLOCATABLE :: d
    real, dimension(:), allocatable :: T_Prof
    real, dimension(:), allocatable :: Z_Prof
    real, dimension(:), allocatable :: Wvmr_Prof
    real, dimension(:), allocatable :: Tpw_Prof
    real, dimension(:), allocatable :: Ozmr_Prof
    real :: satzen_mid_bin
    integer (kind=int1), dimension(:), allocatable :: Inver_Prof
    
  contains
    procedure :: alloc => rtm_params__alloc
    procedure :: dealloc => rtm_params__dealloc  
  end type Rtm_Params

  type (Rtm_Params), public, target, dimension(:, :), ALLOCATABLE :: Rtm
  
  
contains  

   !====================================================================
   !   subroutine Name: ALLOCATE_RTM
   !
   ! Function:
   ! Subroutine that allocates memory for nlon x nlat Rtm structure.
   !
   !====================================================================
   subroutine ALLOCATE_RTM(nx,ny)
      integer (kind=int4), intent(in) :: nx, ny

      integer :: Alloc_Status

      allocate(Rtm(nx,ny),stat=Alloc_Status)
      if (Alloc_Status /= 0) then
         print "(a,'Not enough memory to allocate Rtm_Params structure.')",EXE_PROMPT
         stop
      end if

   end subroutine ALLOCATE_RTM

   !====================================================================
   ! subroutine Name: DEALLOCATE_RTM
   !
   ! Function:
   !  Subroutine that deallocates memory used for nlon x nlat Rtm structure.
   !
   !====================================================================
   subroutine DEALLOCATE_RTM()

      integer :: Alloc_Status

      deallocate(Rtm,stat=Alloc_Status)
      if (Alloc_Status /= 0) then
         print "(a,'Error deallocating Rtm_Params structure.')",EXE_PROMPT
         stop
      end if

   end subroutine DEALLOCATE_RTM


  subroutine rtm_params__alloc ( this, n_levels, n_angles)
    class ( Rtm_params ) :: this
    integer , intent(in) :: n_levels
    integer , intent(in) :: n_angles
    allocate (this % d(n_angles))
    allocate (this % T_Prof(n_levels))
    allocate (this % Z_Prof(n_levels))
    allocate (this % Wvmr_Prof(n_levels))
    allocate (this % Tpw_Prof(n_levels))
    allocate (this % Ozmr_Prof(n_levels))
    allocate (this % Inver_Prof(n_levels))
    this % is_set = .false.
    this % is_allocated = .true.
  end  subroutine rtm_params__alloc
  
   subroutine rtm_params__dealloc ( this )
    class ( Rtm_Params ) :: this
    if ( .not. this % is_allocated) return
    if ( allocated (this % d ) ) deallocate (this % d)
    if ( allocated (this % T_Prof) ) deallocate (this % T_Prof)
    if ( allocated (this % Z_Prof) ) deallocate (this % Z_Prof)
    if ( allocated (this % Wvmr_Prof) ) deallocate (this % Wvmr_Prof)
    if ( allocated (this % Tpw_Prof) ) deallocate (this % Tpw_Prof)
    if ( allocated (this %  Ozmr_Prof) ) deallocate (this % Ozmr_Prof)
    if ( allocated (this %  Inver_Prof) ) deallocate (this %Inver_Prof)
    
    this % Sfc_Level = 0
    this % Inversion_Level = 0
      this % Tropo_Level = 0
      this % Level440 = 0
      this % Level850 = 0
    
    
    this % is_set = .false.
    this % is_allocated = .false.
  end  subroutine rtm_params__dealloc 
  
  
  subroutine Chan_Rtm_Prof__alloc (this, chn_on,nlevels)
      class ( Chan_Rtm_Prof ) :: this
      integer(int1), intent(in) :: chn_on(:)
      integer, intent(in) :: nlevels
      integer :: i
      
      do i = 1, 45 
        if ( chn_on(i) == 1 ) then
            call this%ch(i)%alloc(nlevels)
        end if
      
      end do
      this % is_set = .false.
       this % is_allocated = .true.
  end subroutine Chan_Rtm_Prof__alloc

  
  subroutine rtm_prof__alloc ( this, n_levels)
    class ( Rtm_Prof ) :: this
    integer , intent(in) :: n_levels
    allocate (this % Trans_Atm_Profile(n_levels))
    allocate (this % Trans_Atm_Solar_Profile(n_levels))
    allocate (this % Trans_Atm_Total_Profile(n_levels))
    allocate (this % Rad_Atm_Profile(n_levels))
    allocate (this % Rad_Atm_Dwn_Profile(n_levels))
    allocate (this % Rad_BB_Cloud_Profile(n_levels))
    
    this % is_set = .false.
    this % is_allocated = .true.
  end  subroutine rtm_prof__alloc
  
   subroutine rtm_prof__dealloc ( this )
    class ( Rtm_Prof ) :: this
    if ( allocated (this % Trans_Atm_Profile) ) deallocate (this % Trans_Atm_Profile)
    if ( allocated (this % Trans_Atm_Solar_Profile) ) deallocate (this % Trans_Atm_Solar_Profile)
    if ( allocated (this % Trans_Atm_Total_Profile) ) deallocate (this % Trans_Atm_Total_Profile)
    if ( allocated (this % Rad_Atm_Profile) ) deallocate (this % Rad_Atm_Profile)
    if ( allocated (this % Rad_Atm_Dwn_Profile) ) deallocate (this % Rad_Atm_Dwn_Profile)
    if ( allocated (this %  Rad_BB_Cloud_Profile) ) deallocate (this % Rad_BB_Cloud_Profile)
    
    this % is_set = .false.
    this % is_allocated = .false.
  end  subroutine rtm_prof__dealloc 
  

end module RTM_COMMON_MOD
