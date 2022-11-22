!$Id: acha_clavrx_services_module.f90 4092 2021-03-02 22:05:14Z heidinger $
!------------------------------------------------------------------------------
!this module holds all the dependencies for ACHA for the various frameworks
!------------------------------------------------------------------------------
module ACHA_SERVICES_MOD

 use PLANCK_MOD
 use CONSTANTS_MOD
 use NWP_COMMON_MOD
 use RTM_COMMON_MOD
 use PIXEL_COMMON_MOD, only: &
       Ch, &
       Nav, &
       Geo, &
       Sensor, &
       Image, &
       ACHA, &
       CLDMASK, &
       NWP_PIX, &
       Sfc, &
       Bad_Pixel_Mask, &
       I_Lrc, &
       J_Lrc, &
       Zen_Idx_Rtm, &
       Process_Undetected_Cloud_Flag, &
       Cld_Type, &
       Tc_Cirrus_Background, &
       Diag_Pix_Array_1, &
       Diag_Pix_Array_2, &
       Diag_Pix_Array_3, &
       Elem_Abs_Idx_Acha_Dump, &
       Line_Abs_Idx_Acha_Dump, &
       Ancil_Data_Dir, &
       Cld_Phase_Uncertainty, &
       Tc_Opaque_Cloud, Pc_Opaque_Cloud, Zc_Opaque_Cloud
 
 use NUMERICAL_ROUTINES_MOD, only: INVERT_MATRIX, LOCATE

 use CX_STRING_TOOLS_MOD, only: COUNTSUBSTRING

 implicit none

 public:: ACHA_FETCH_PIXEL_NWP_RTM 
 private:: INTERPOLATE_PROFILE_ACHA

 integer(KIND=INT4), PRIVATE, PARAMETER :: Num_Levels_Rtm_Prof = 101

!ACHA dump structure
 type, public :: acha_dump_struct
  integer:: Segment_Number
  integer:: Number_Lines_Per_Segment
  integer:: Elem_Abs_Idx
  integer:: Line_Abs_Idx
 end type acha_dump_struct

!ACHA diagnostic structure
 type, public :: acha_diag_struct
  real (kind=real4), dimension(:,:), pointer:: Array_1
  real (kind=real4), dimension(:,:), pointer:: Array_2
  real (kind=real4), dimension(:,:), pointer:: Array_3
 end type acha_diag_struct

!ACHA input structure
 type, public :: acha_input_struct
 character (len=50) :: ACHA_Mode_Flag_In
 logical :: Use_10_4
 integer (kind=int4):: Number_of_Elements
 integer (kind=int4):: Number_Of_Lines
 integer (kind=int4):: Smooth_Nwp_Fields_Flag
 integer (kind=int4):: Process_Undetected_Cloud_Flag
 real (kind=real4):: Sensor_Resolution_KM
 integer (kind=int4):: WMO_Id

 !-- local channel indices
 integer:: Chan_Idx_038um
 integer:: Chan_Idx_062um
 integer:: Chan_Idx_067um
 integer:: Chan_Idx_073um
 integer:: Chan_Idx_085um
 integer:: Chan_Idx_097um
 integer:: Chan_Idx_104um
 integer:: Chan_Idx_110um
 integer:: Chan_Idx_120um
 integer:: Chan_Idx_133um 
 integer:: Chan_Idx_136um 
 integer:: Chan_Idx_139um 
 integer:: Chan_Idx_142um 
 integer:: Chan_On_038um
 integer:: Chan_On_062um
 integer:: Chan_On_067um
 integer:: Chan_On_073um
 integer:: Chan_On_085um
 integer:: Chan_On_097um
 integer:: Chan_On_104um
 integer:: Chan_On_110um
 integer:: Chan_On_120um
 integer:: Chan_On_133um
 integer:: Chan_On_136um
 integer:: Chan_On_139um
 integer:: Chan_On_142um

 !-- local pointers that point to global variables
 integer (kind=int1), dimension(:,:), pointer:: Invalid_Data_Mask
 integer (kind=int1), dimension(:,:), pointer:: DQF_038um
 integer (kind=int1), dimension(:,:), pointer:: DQF_062um
 integer (kind=int1), dimension(:,:), pointer:: DQF_067um
 integer (kind=int1), dimension(:,:), pointer:: DQF_073um
 integer (kind=int1), dimension(:,:), pointer:: DQF_085um
 integer (kind=int1), dimension(:,:), pointer:: DQF_097um
 integer (kind=int1), dimension(:,:), pointer:: DQF_104um
 integer (kind=int1), dimension(:,:), pointer:: DQF_110um
 integer (kind=int1), dimension(:,:), pointer:: DQF_120um
 integer (kind=int1), dimension(:,:), pointer:: DQF_133um
 integer (kind=int1), dimension(:,:), pointer:: DQF_136um
 integer (kind=int1), dimension(:,:), pointer:: DQF_139um
 integer (kind=int1), dimension(:,:), pointer:: DQF_142um
 real, dimension(:,:), pointer:: Bt_038um
 real, dimension(:,:), pointer:: Bt_062um
 real, dimension(:,:), pointer:: Bt_067um
 real, dimension(:,:), pointer:: Bt_073um
 real, dimension(:,:), pointer:: Bt_085um
 real, dimension(:,:), pointer:: Bt_097um
 real, dimension(:,:), pointer:: Bt_104um
 real, dimension(:,:), pointer:: Bt_110um
 real, dimension(:,:), pointer:: Bt_120um
 real, dimension(:,:), pointer:: Bt_133um
 real, dimension(:,:), pointer:: Bt_136um
 real, dimension(:,:), pointer:: Bt_139um
 real, dimension(:,:), pointer:: Bt_142um
 real, dimension(:,:), pointer:: Rad_038um
 real, dimension(:,:), pointer:: Rad_062um
 real, dimension(:,:), pointer:: Rad_067um
 real, dimension(:,:), pointer:: Rad_073um
 real, dimension(:,:), pointer:: Rad_085um
 real, dimension(:,:), pointer:: Rad_097um
 real, dimension(:,:), pointer:: Rad_104um
 real, dimension(:,:), pointer:: Rad_110um
 real, dimension(:,:), pointer:: Rad_120um
 real, dimension(:,:), pointer:: Rad_133um
 real, dimension(:,:), pointer:: Rad_136um
 real, dimension(:,:), pointer:: Rad_139um
 real, dimension(:,:), pointer:: Rad_142um
 real, dimension(:,:), pointer:: Cosine_Zenith_Angle
 real, dimension(:,:), pointer:: Sensor_Zenith_Angle
 real, dimension(:,:), pointer:: Sensor_Azimuth_Angle
 real, dimension(:,:), pointer:: Surface_Temperature
 real, dimension(:,:), pointer:: Surface_Air_Temperature
 real, dimension(:,:), pointer:: Tropopause_Temperature
 real, dimension(:,:), pointer:: Tropopause_Height
 real, dimension(:,:), pointer:: Tropopause_Pressure
 real, dimension(:,:), pointer:: Surface_Pressure
 real, dimension(:,:), pointer:: Surface_Elevation
 real, dimension(:,:), pointer:: Latitude
 real, dimension(:,:), pointer:: Longitude
 real, dimension(:,:), pointer:: Rad_Clear_038um
 real, dimension(:,:), pointer:: Rad_Clear_062um
 real, dimension(:,:), pointer:: Rad_Clear_067um
 real, dimension(:,:), pointer:: Rad_Clear_073um
 real, dimension(:,:), pointer:: Rad_Clear_085um
 real, dimension(:,:), pointer:: Rad_Clear_097um
 real, dimension(:,:), pointer:: Rad_Clear_104um
 real, dimension(:,:), pointer:: Rad_Clear_110um
 real, dimension(:,:), pointer:: Rad_Clear_120um
 real, dimension(:,:), pointer:: Rad_Clear_133um
 real, dimension(:,:), pointer:: Rad_Clear_136um
 real, dimension(:,:), pointer:: Rad_Clear_139um
 real, dimension(:,:), pointer:: Rad_Clear_142um
 real, dimension(:,:), pointer:: Surface_Emissivity_038um
 real, dimension(:,:), pointer:: Surface_Emissivity_062um 
 real, dimension(:,:), pointer:: Surface_Emissivity_067um 
 real, dimension(:,:), pointer:: Surface_Emissivity_073um 
 real, dimension(:,:), pointer:: Surface_Emissivity_085um
 real, dimension(:,:), pointer:: Surface_Emissivity_097um
 real, dimension(:,:), pointer:: Surface_Emissivity_104um
 real, dimension(:,:), pointer:: Surface_Emissivity_110um
 real, dimension(:,:), pointer:: Surface_Emissivity_120um
 real, dimension(:,:), pointer:: Surface_Emissivity_133um
 real, dimension(:,:), pointer:: Surface_Emissivity_136um
 real, dimension(:,:), pointer:: Surface_Emissivity_139um
 real, dimension(:,:), pointer:: Surface_Emissivity_142um
 integer (kind=int1),dimension(:,:), pointer:: Snow_Class
 integer (kind=int1),dimension(:,:), pointer:: Surface_Type
 integer (kind=int1),dimension(:,:), pointer:: Cloud_Mask
 real, dimension(:,:), pointer:: Cloud_Probability
 real, dimension(:,:), pointer:: Ice_Cloud_Probability
 real, dimension(:,:), pointer:: Cloud_Phase_Uncertainty
 integer (kind=int1),dimension(:,:), pointer:: Cloud_Type
 integer (kind=int4), dimension(:,:), pointer:: Elem_Idx_NWP 
 integer (kind=int4), dimension(:,:), pointer:: Line_Idx_NWP 
 integer (kind=int4), dimension(:,:), pointer:: Elem_Idx_Opposite_Corner_NWP 
 integer (kind=int4), dimension(:,:), pointer:: Line_Idx_Opposite_Corner_NWP 
 integer (kind=int4), dimension(:,:), pointer:: Viewing_Zenith_Angle_Idx_Rtm
 real (kind=real4), dimension(:,:), pointer:: Latitude_Interp_Weight_NWP
 real (kind=real4), dimension(:,:), pointer:: Longitude_Interp_Weight_NWP
 integer(kind=int4), dimension(:,:), pointer :: Elem_Idx_LRC_Input
 integer(kind=int4), dimension(:,:), pointer :: Line_Idx_LRC_Input
 real (kind=real4), dimension(:,:), pointer:: Tc_Cirrus_Sounder
 real (kind=real4), dimension(:,:), pointer:: Tc_Nwp
 real (kind=real4), dimension(:,:), pointer:: Tc_Opaque
 real (kind=real4), dimension(:,:), pointer:: Zc_Opaque
 real (kind=real4), dimension(:,:), pointer:: Pc_Opaque
 real (kind=real4), dimension(:,:), pointer:: Tc_Ap
 real (kind=real4), dimension(:,:), pointer:: Ec_Ap
 real (kind=real4), dimension(:,:), pointer:: Beta_Ap
 real (kind=real4), dimension(:,:), pointer:: Ice_Prob_Ap
 real (kind=real4), dimension(:,:), pointer:: Tc_Lower_Ap
 
 end type acha_input_struct

 !---RTM and NWP pixel level structure
 type, public :: acha_rtm_nwp_struct

   integer:: Smooth_Nwp_Fields_Flag

   !-- NWP Levels
   integer:: Sfc_Level
   integer:: Tropo_Level

   !-- RTM profiles
   real, dimension(:), pointer :: Atm_Rad_Prof_038um
   real, dimension(:), pointer :: Atm_Rad_Prof_062um
   real, dimension(:), pointer :: Atm_Rad_Prof_067um
   real, dimension(:), pointer :: Atm_Rad_Prof_073um
   real, dimension(:), pointer :: Atm_Rad_Prof_085um
   real, dimension(:), pointer :: Atm_Rad_Prof_097um
   real, dimension(:), pointer :: Atm_Rad_Prof_104um
   real, dimension(:), pointer :: Atm_Rad_Prof_110um
   real, dimension(:), pointer :: Atm_Rad_Prof_120um
   real, dimension(:), pointer :: Atm_Rad_Prof_133um
   real, dimension(:), pointer :: Atm_Rad_Prof_136um
   real, dimension(:), pointer :: Atm_Rad_Prof_139um
   real, dimension(:), pointer :: Atm_Rad_Prof_142um
   real, dimension(:), pointer :: Atm_Trans_Prof_038um
   real, dimension(:), pointer :: Atm_Trans_Prof_062um
   real, dimension(:), pointer :: Atm_Trans_Prof_067um
   real, dimension(:), pointer :: Atm_Trans_Prof_073um
   real, dimension(:), pointer :: Atm_Trans_Prof_085um
   real, dimension(:), pointer :: Atm_Trans_Prof_097um
   real, dimension(:), pointer :: Atm_Trans_Prof_104um
   real, dimension(:), pointer :: Atm_Trans_Prof_110um
   real, dimension(:), pointer :: Atm_Trans_Prof_120um
   real, dimension(:), pointer :: Atm_Trans_Prof_133um
   real, dimension(:), pointer :: Atm_Trans_Prof_136um
   real, dimension(:), pointer :: Atm_Trans_Prof_139um
   real, dimension(:), pointer :: Atm_Trans_Prof_142um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_038um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_062um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_067um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_073um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_085um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_097um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_104um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_110um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_120um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_133um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_136um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_139um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_142um

   !-- NWP profiles
   real, dimension(:), pointer :: T_Prof
   real, dimension(Num_Levels_Rtm_Prof) :: P_Prof
   real, dimension(:), pointer :: Z_Prof

end type acha_rtm_nwp_struct

!output structure
 type, public :: acha_output_struct
   real, dimension(:,:), pointer:: Latitude_Pc
   real, dimension(:,:), pointer:: Longitude_Pc
   real, dimension(:,:), pointer:: Tc_Ap_Uncer
   real, dimension(:,:), pointer:: Tc_Ap
   real, dimension(:,:), pointer:: Tc
   real, dimension(:,:), pointer:: Tfm
   real, dimension(:,:), pointer:: Es
   real, dimension(:,:), pointer:: Zc_rtm
   real, dimension(:,:), pointer:: Zs
   real, dimension(:,:), pointer:: Ts
   real, dimension(:,:), pointer:: Ec
   real, dimension(:,:), pointer:: Beta
   real, dimension(:,:), pointer:: Pc
   real, dimension(:,:), pointer:: Zc
   real, dimension(:,:), pointer:: Zc_Base
   real, dimension(:,:), pointer:: Tau
   real, dimension(:,:), pointer:: Tau_Uncertainty
   real, dimension(:,:), pointer:: Reff
   real, dimension(:,:), pointer:: Pc_Eff
   real, dimension(:,:), pointer:: Tc_Eff
   real, dimension(:,:), pointer:: Zc_Eff
   real, dimension(:,:), pointer:: Ice_Probability
   real, dimension(:,:), pointer:: Tc_Uncertainty
   real, dimension(:,:), pointer:: Ec_Uncertainty
   real, dimension(:,:), pointer:: Beta_Uncertainty
   real, dimension(:,:), pointer:: Pc_Uncertainty
   real, dimension(:,:), pointer:: Zc_Uncertainty
   real, dimension(:,:), pointer:: Lower_Tc_Uncertainty
   real, dimension(:,:), pointer:: Lower_Zc_Uncertainty
   real, dimension(:,:), pointer:: Lower_Pc_Uncertainty
   real, dimension(:,:), pointer:: Ice_Probability_Uncertainty
   real, dimension(:,:), pointer:: Lower_Pc
   real, dimension(:,:), pointer:: Lower_Tc
   real, dimension(:,:), pointer:: Lower_Zc
   real, dimension(:,:), pointer:: Lower_Tc_Ap
   real, dimension(:,:), pointer:: Conv_Test
   real, dimension(:,:), pointer:: Cost
   real, dimension(:,:), pointer:: Goodness
   real, dimension(:,:), pointer:: Ec_038um
   real, dimension(:,:), pointer:: Ec_062um
   real, dimension(:,:), pointer:: Ec_067um
   real, dimension(:,:), pointer:: Ec_073um
   real, dimension(:,:), pointer:: Ec_085um
   real, dimension(:,:), pointer:: Ec_097um
   real, dimension(:,:), pointer:: Ec_104um
   real, dimension(:,:), pointer:: Ec_110um
   real, dimension(:,:), pointer:: Ec_120um
   real, dimension(:,:), pointer:: Ec_133um
   real, dimension(:,:), pointer:: Ec_136um
   real, dimension(:,:), pointer:: Ec_139um
   real, dimension(:,:), pointer:: Ec_142um
   integer (kind=int1), dimension(:,:), pointer:: Qf
   integer (kind=int1), dimension(:,:,:), pointer:: OE_Qf
   integer (kind=int1), dimension(:,:), pointer :: Packed_Qf
   integer (kind=int1), dimension(:,:), pointer :: Packed_Meta_Data
   integer(kind=int1), dimension(:,:), pointer :: Processing_Order   
   integer(kind=int1), dimension(:,:), pointer :: Inversion_Flag
   integer(kind=int1), dimension(:,:), pointer :: Cloud_Type
  end type acha_output_struct
  
!Symbol stucture

 type, public :: acha_symbol_struct
    integer(kind=int1) :: CLOUDY
    integer(kind=int1) :: PROB_CLOUDY
    integer(kind=int1) :: PROB_CLEAR
    integer(kind=int1) :: CLEAR

    integer(kind=int1) :: NO
    integer(kind=int1) :: YES

    integer(kind=int1) :: WATER_SFC
    integer(kind=int1) :: EVERGREEN_NEEDLE_SFC
    integer(kind=int1) :: EVERGREEN_BROAD_SFC
    integer(kind=int1) :: DECIDUOUS_NEEDLE_SFC
    integer(kind=int1) :: DECIDUOUS_BROAD_SFC
    integer(kind=int1) :: MIXED_FORESTS_SFC
    integer(kind=int1) :: WOODLANDS_SFC
    integer(kind=int1) :: WOODED_GRASS_SFC
    integer(kind=int1) :: CLOSED_SHRUBS_SFC
    integer(kind=int1) :: OPEN_SHRUBS_SFC
    integer(kind=int1) :: GRASSES_SFC
    integer(kind=int1) :: CROPLANDS_SFC
    integer(kind=int1) :: BARE_SFC
    integer(kind=int1) :: URBAN_SFC

    integer(kind=int1) :: SHALLOW_OCEAN
    integer(kind=int1) :: LAND
    integer(kind=int1) :: COASTLINE
    integer(kind=int1) :: SHALLOW_INLAND_WATER
    integer(kind=int1) :: EPHEMERAL_WATER
    integer(kind=int1) :: DEEP_INLAND_WATER
    integer(kind=int1) :: MODERATE_OCEAN
    integer(kind=int1) :: DEEP_OCEAN

    integer(kind=int1) :: NO_SNOW
    integer(kind=int1) :: SEA_ICE
    integer(kind=int1) :: SNOW

    integer(kind=int1) :: CLEAR_TYPE
    integer(kind=int1) :: PROB_CLEAR_TYPE
    integer(kind=int1) :: FOG_TYPE
    integer(kind=int1) :: WATER_TYPE
    integer(kind=int1) :: SUPERCOOLED_TYPE
    integer(kind=int1) :: MIXED_TYPE
    integer(kind=int1) :: OPAQUE_ICE_TYPE
    integer(kind=int1) :: TICE_TYPE
    integer(kind=int1) :: CIRRUS_TYPE
    integer(kind=int1) :: OVERLAP_TYPE
    integer(kind=int1) :: OVERSHOOTING_TYPE
    integer(kind=int1) :: UNKNOWN_TYPE
    integer(kind=int1) :: DUST_TYPE
    integer(kind=int1) :: SMOKE_TYPE
    integer(kind=int1) :: FIRE_TYPE

    integer(kind=int1) :: CLEAR_PHASE
    integer(kind=int1) :: WATER_PHASE
    integer(kind=int1) :: SUPERCOOLED_PHASE
    integer(kind=int1) :: MIXED_PHASE
    integer(kind=int1) :: ICE_PHASE
    integer(kind=int1) :: UNKNOWN_PHASE
 end type acha_symbol_struct

 logical, public :: ABI_Use_104um_Flag
 
 contains

!-------------------------------------------------------------------------------
! This subroutine gathers the necessary NWP and RTM profiles used for a given
! pixel for ACHA. 
!--------------------------------------------------------------------------------
 subroutine  ACHA_FETCH_PIXEL_NWP_RTM(Acha_Input, Symbol, &
                                      Elem_Idx, Line_Idx, Acha_RTM_NWP,Error_Flag)
                                      
   type(acha_input_struct), intent(inout) :: Acha_Input
   type(acha_rtm_nwp_struct), intent(inout) :: Acha_RTM_NWP
   type(acha_symbol_struct), intent(inout) :: Symbol
   integer, intent(in) :: Elem_Idx
   integer, intent(in) :: Line_Idx
   integer, intent(out):: Error_Flag
   integer:: Ivza
   integer:: Inwp
   integer:: Jnwp
   integer:: Inwp_x
   integer:: Jnwp_x
   real:: Inwp_Weight
   real:: Jnwp_Weight

   Error_Flag = 1   !initialize as bad

   Inwp = Acha_Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx)
   Jnwp = Acha_Input%Line_Idx_Nwp(Elem_Idx,Line_Idx)
   
   Inwp_x = Acha_Input%Elem_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
   Jnwp_x = Acha_Input%Line_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
   
   Inwp_Weight = Acha_Input%Longitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
   Jnwp_Weight = Acha_Input%Latitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)

   Ivza =  Acha_Input%Viewing_Zenith_Angle_Idx_Rtm(Elem_Idx,Line_Idx)

   !--- populate height and temperature profiles
   if (Inwp <= 0 .or. Jnwp <= 0) then
     print *, "ACHA_FETCH_PIXEL_NWP_RTM: Bad Nwp Indices"
     return
   endif
   if (Allocated(Rtm(Inwp,Jnwp)%T_Prof) .eqv. .false.) then
     print *, "ACHA_FETCH_PIXEL_NWP_RTM: Error, T_Prof not allocated"
     return
   endif

   !initialize smooth NWP flag 
   Acha_RTM_NWP%Sfc_Level = Rtm(Inwp,Jnwp)%Sfc_Level
   Acha_RTM_NWP%Tropo_Level = Rtm(Inwp,Jnwp)%Tropo_Level
   
   !--- do various 101 level NWP Profiles
   Acha_RTM_NWP%P_Prof = P_Std_Rtm
   Acha_RTM_NWP%T_Prof => Rtm(Inwp,Jnwp)%T_Prof 
   Acha_RTM_NWP%Z_Prof => Rtm(Inwp,Jnwp)%Z_Prof 

   ACHA_RTM_NWP%Smooth_Nwp_Fields_Flag = Acha_Input%Smooth_Nwp_Fields_Flag
   !------------------------------------------------------
   ! Before smoothing profiles, ensure that all required
   ! rtm profiles are populated, if not, skip smoothing
   !------------------------------------------------------
   
    !do smoothing routines here 
    if ((ACHA_RTM_NWP%Smooth_Nwp_Fields_Flag == Symbol%YES) .and. &
        Rtm(Inwp,Jnwp)%is_set .and. &
        Rtm(Inwp_x,Jnwp)%is_set .and. &
        Rtm(Inwp,Jnwp_x)%is_set .and. &
        Rtm(Inwp_x,Jnwp_x)%is_set) then

       ACHA_RTM_NWP%Z_Prof = INTERPOLATE_PROFILE_ACHA(Rtm(Inwp,Jnwp)%Z_Prof, &
                                                 Rtm(Inwp_x,Jnwp)%Z_Prof, & 
                                                 Rtm(Inwp,Jnwp_x)%Z_Prof, & 
                                                 Rtm(Inwp_x,Jnwp_x)%Z_Prof, &
                                                 Inwp_Weight,Jnwp_Weight)

       ACHA_RTM_NWP%T_Prof = INTERPOLATE_PROFILE_ACHA(Rtm(Inwp,Jnwp)%T_Prof, &
                                                 Rtm(Inwp_x,Jnwp)%T_Prof, & 
                                                 Rtm(Inwp,Jnwp_x)%T_Prof, & 
                                                 Rtm(Inwp_x,Jnwp_x)%T_Prof, &
                                                 Inwp_Weight,Jnwp_Weight)
    endif

   !---- RTM profiles
 
   !--- populate radiance and transmission profiles -- CAN'T WE USE CHAN_IDX here?
   if (Acha_Input%Chan_On_038um == sym%YES) then
     Acha_RTM_NWP%Atm_Rad_Prof_038um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(20)%Rad_Atm_Profile
     Acha_RTM_NWP%Atm_Trans_Prof_038um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(20)%Trans_Atm_Profile
     Acha_RTM_NWP%Black_Body_Rad_Prof_038um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(20)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_062um == sym%YES) then
     Acha_RTM_NWP%Atm_Rad_Prof_062um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(37)%Rad_Atm_Profile
     Acha_RTM_NWP%Atm_Trans_Prof_062um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(37)%Trans_Atm_Profile
     Acha_RTM_NWP%Black_Body_Rad_Prof_062um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(37)%Rad_BB_Cloud_Profile
   endif
 
   if (Acha_Input%Chan_On_067um == sym%YES) then
     Acha_RTM_NWP%Atm_Rad_Prof_067um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(27)%Rad_Atm_Profile
     Acha_RTM_NWP%Atm_Trans_Prof_067um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(27)%Trans_Atm_Profile
     Acha_RTM_NWP%Black_Body_Rad_Prof_067um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(27)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_073um == sym%YES) then
     Acha_RTM_NWP%Atm_Rad_Prof_073um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(28)%Rad_Atm_Profile
     Acha_RTM_NWP%Atm_Trans_Prof_073um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(28)%Trans_Atm_Profile
     Acha_RTM_NWP%Black_Body_Rad_Prof_073um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(28)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_085um == sym%YES) then
     Acha_RTM_NWP%Atm_Rad_Prof_085um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(29)%Rad_Atm_Profile
     Acha_RTM_NWP%Atm_Trans_Prof_085um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(29)%Trans_Atm_Profile
     Acha_RTM_NWP%Black_Body_Rad_Prof_085um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(29)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_097um == sym%YES) then
     Acha_RTM_NWP%Atm_Rad_Prof_097um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(30)%Rad_Atm_Profile
     Acha_RTM_NWP%Atm_Trans_Prof_097um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(30)%Trans_Atm_Profile
     Acha_RTM_NWP%Black_Body_Rad_Prof_097um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(30)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_104um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_104um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(38)%Rad_Atm_Profile
      Acha_RTM_NWP%Atm_Trans_Prof_104um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(38)%Trans_Atm_Profile
      Acha_RTM_NWP%Black_Body_Rad_Prof_104um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(38)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_110um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_110um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Rad_Atm_Profile
      Acha_RTM_NWP%Atm_Trans_Prof_110um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Trans_Atm_Profile
      Acha_RTM_NWP%Black_Body_Rad_Prof_110um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Rad_BB_Cloud_Profile
   endif
   
   if (Acha_Input%Chan_On_120um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_120um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(32)%Rad_Atm_Profile
      Acha_RTM_NWP%Atm_Trans_Prof_120um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(32)%Trans_Atm_Profile
      Acha_RTM_NWP%Black_Body_Rad_Prof_120um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(32)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_133um == sym%YES) then
       Acha_RTM_NWP%Atm_Rad_Prof_133um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(33)%Rad_Atm_Profile
       Acha_RTM_NWP%Atm_Trans_Prof_133um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(33)%Trans_Atm_Profile
       Acha_RTM_NWP%Black_Body_Rad_Prof_133um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(33)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_136um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_136um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(34)%Rad_Atm_Profile
      Acha_RTM_NWP%Atm_Trans_Prof_136um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(34)%Trans_Atm_Profile
      Acha_RTM_NWP%Black_Body_Rad_Prof_136um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(34)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_139um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_139um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(35)%Rad_Atm_Profile
      Acha_RTM_NWP%Atm_Trans_Prof_139um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(35)%Trans_Atm_Profile
      Acha_RTM_NWP%Black_Body_Rad_Prof_139um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(35)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_142um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_142um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(36)%Rad_Atm_Profile
      Acha_RTM_NWP%Atm_Trans_Prof_142um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(36)%Trans_Atm_Profile
      Acha_RTM_NWP%Black_Body_Rad_Prof_142um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(36)%Rad_BB_Cloud_Profile
   endif

   !--- set error flag to good
   Error_Flag = 0

 end subroutine ACHA_FETCH_PIXEL_NWP_RTM

!----------------------------------------------------------------------------
! Function INTERPOLATE_PROFILE_ACHA
!
! general interpoLation routine for profiles
!
! input:
! lonx - longitude weighting factor
! Latx = Latitude weighting factor
! z1 = data(ilon, iLat)
! z2 = data(ilonx,iLat)
! z3 = data(ilon,iLatx)
! z4 = data(ilonx,iLatx)
!
! output:
! z = interpoLated profile
!
!
!---------------------------------------------------------------------------
 function INTERPOLATE_PROFILE_ACHA(z1,z2,z3,z4,lonx,Latx) result(z)

  real, dimension(:), intent(in):: z1
  real, dimension(:), intent(in):: z2
  real, dimension(:), intent(in):: z3
  real, dimension(:), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: Latx
  real, dimension(size(z1)):: z

  !--- linear inteprpoLation scheme
  z =  (1.0-lonx) * ((1.0-Latx) * z1 + (Latx)* z3) + &
           (lonx) * ((1.0-Latx) * z2 + (Latx)* z4)

 end function INTERPOLATE_PROFILE_ACHA


end module ACHA_SERVICES_MOD
