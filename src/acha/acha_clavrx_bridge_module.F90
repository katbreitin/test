!$Id: acha_clavrx_bridge_module.f90 4092 2021-03-02 22:05:14Z heidinger $
!------------------------------------------------------------------------------
!  NOAA AWG Cloud Height Algorithm (ACHA) Bridge Code
!
!  This module houses the routines that serve as a bridge between
!  the CLAVR-x processing system and the ACHA code.
!
!------------------------------------------------------------------------------
module ACHA_CLAVRX_BRIDGE

 use AWG_CLOUD_HEIGHT
 use ACHA_COMP
 use ACHA_SHADOW
 use ACHA_SERVICES_MOD
 use CLAVRX_MESSAGE_MOD, only: MESG
 use LEVEL2_STRUCTURES_MOD, only: &
     Clavrx_Global_Attr

 use CX_ABI_LHP_MOD_ALG_FUNCTION

 use CALIBRATION_CONSTANTS_MOD, only: &
     ABI_FPT_Thresh_038um, ABI_FPT_Thresh_062um, ABI_FPT_Thresh_067um &
   , ABI_FPT_Thresh_073um, ABI_FPT_Thresh_085um, ABI_FPT_Thresh_097um &
   , ABI_FPT_Thresh_104um, ABI_FPT_Thresh_110um, ABI_FPT_Thresh_120um &
   , ABI_FPT_Thresh_133um
   
 implicit none

 public :: AWG_CLOUD_HEIGHT_BRIDGE
 private:: SET_SYMBOL, SET_INPUT, SET_OUTPUT, NULL_INPUT, NULL_OUTPUT
 private:: SET_DIAG, NULL_DIAG, SET_DUMP
 private :: LHP_CHN_CHECK   !WCS3

 !--------------------------------------------------------------------
 ! define structures that will be arguments to ACHA
 !--------------------------------------------------------------------
 type(acha_symbol_struct), private :: Symbol
 type(acha_input_struct), private :: Input
 type(acha_output_struct), private :: Output
 type(acha_diag_struct), private :: Diag
 type(acha_dump_struct), private :: Dump

 contains

!----------------------------------------------------------------------
! ACHA BRIDGE subroutine
!---------------------------------------------------------------------- 
 subroutine AWG_CLOUD_HEIGHT_BRIDGE()
 
   implicit none

   logical, save:: First_Call = .true.

   logical:: DQF_Set_Flag

   if (First_Call .eqv. .true.) then
       call MESG('ACHA starts ', color = 46)
   endif

   !---null pointers before filling them
   call NULL_INPUT()
   call NULL_OUTPUT()

   !-------------------------------------------
   !--- initialize structures
   !-------------------------------------------

   !--- store integer values
   call SET_INPUT()

   !---- initalize Output structure
   call SET_OUTPUT()
  
   !----set symbols to local values
   call SET_SYMBOL()

   !---- initialize diagnostic structure
   call SET_DIAG()

   !---- initialize dump structure
   call SET_DUMP()

   !--- Determine if there is any DQF information
   DQF_Set_Flag = .false.
   if (Input%Chan_On_104um == 1_int4) then
     if (maxval(Input%DQF_104um) > Missing_Value_Int1) then
      DQF_Set_Flag = .true.
     endif
   endif
   if (Input%Chan_On_110um == 1_int4) then
     if (maxval(Input%DQF_110um) > Missing_Value_Int1) then
      DQF_Set_Flag = .true.
     endif
   endif

   !--- Determine 10.4 vs 11 micron channel and DQF impacts on Mode (GOES-RU only)
   Input%Use_10_4 = .false.
   ABI_Use_104um_Flag = .false.

   !--- if DQFs not set, skip using them for acha mode 
   if (DQF_Set_Flag) then

      if (Sensor%WMO_id == 271) then

        !--- Turn off GOES-17 channels if focal plane temperatures are too warm.
        CALL LHP_CHN_CHECK()

        call MODifY_MODE_USING_LHP_THRESHOLDS(Input%Acha_Mode_Flag_In, &
                                              Input%Chan_On_038um, &
                                              Input%Chan_On_067um, &
                                              Input%Chan_On_085um, &
                                              Input%Chan_On_104um, &
                                              Input%Chan_On_110um, &
                                              Input%Chan_On_120um, &
                                              Input%Chan_On_133um)

        !call SET_10_4_Flag(Input%Chan_On_104um, Input%Chan_On_110um, Input%Chan_On_120um, Input%Chan_On_133um, &
        !                   Input%DQF_104um,Input%DQF_110um,Input%DQF_120um,Input%DQF_133um,Input%Use_10_4)

        !call MODifY_MODE_USING_DQF(Input%Invalid_Data_Mask, &
        !                           Input%Chan_On_038um, Input%Chan_On_067um, Input%Chan_On_085um, &
        !                           Input%Chan_On_120um, Input%Chan_On_133um, &
        !                           Input%DQF_038um, Input%DQF_067um,  &
        !                           Input%DQF_085um, Input%DQF_120um, Input%DQF_133um, &
        !                           Input%Acha_Mode_Flag_In)
      endif
   endif

   !-----------------------------------------------------------------------
   !--- Call to AWG Cloud Height Algorithm (ACHA)
   !-----------------------------------------------------------------------
#ifdef ACHADIAG
#ifdef ACHADUMP
   call AWG_CLOUD_HEIGHT_ALGORITHM(Input, Symbol, Output, Dump=Dump, Diag=Diag)
#else
   call AWG_CLOUD_HEIGHT_ALGORITHM(Input, Symbol, Output, Diag=Diag)
#endif
#else
#ifdef ACHADUMP
   call AWG_CLOUD_HEIGHT_ALGORITHM(Input, Symbol, Output, Dump=Dump)
#else
   call AWG_CLOUD_HEIGHT_ALGORITHM(Input, Symbol, Output)
#endif
#endif
   

   !-----------------------------------------------------------------------
   !--- Call algorithm to make ACHA optical and microphysical properties
   !-----------------------------------------------------------------------
   call ACHA_COMP_ALGORITHM(Input, Symbol, Output)


#ifdef SHADOWON
   !-----------------------------------------------------------------------
   !--- Call to Geometrical Shadow Algorithm
   !-----------------------------------------------------------------------
  
   call CLOUD_SHADOW_RETR (  &
           ACHA%Zc &
         , Geo%Solaz &
         , Geo%Solzen &
         , Nav%Lat &
         , Nav%Lon &
         , Nav%Lat_Pc &
         , Nav%Lon_Pc &
         , CLDMASK%Shadow_Mask ) 
 
   !!---- copy shadow result into cloud mask test bits
   where (CLDMASK%Shadow_Mask == 1 .and. CLDMASK%Cld_Mask == 0 )  
           CLDMASK%Cld_Test_Vector_Packed ( 2 , :, : )  = ibset (CLDMASK%Cld_Test_Vector_Packed ( 2 , :, : )  , 6 )
   end where
#endif

   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_INPUT()
   call NULL_OUTPUT()

   !-----------------------------------------------------------------------
   !---  read CVS Tag from ACHA and store in global variable for output
   !-----------------------------------------------------------------------
   call SET_ACHA_VERSION(Acha_Version)
  
   First_Call = .false.

 end subroutine AWG_CLOUD_HEIGHT_BRIDGE

 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding input 
 !-----------------------------------------------------------------------------
 subroutine NULL_INPUT()
     Input%Invalid_Data_Mask =>  null()
     Input%DQF_038um =>  null()
     Input%DQF_062um =>  null()
     Input%DQF_067um =>  null()
     Input%DQF_073um =>  null()
     Input%DQF_085um =>  null()
     Input%DQF_097um =>  null()
     Input%DQF_104um =>  null()
     Input%DQF_110um =>  null()
     Input%DQF_120um =>  null()
     Input%DQF_133um =>  null()
     Input%DQF_136um =>  null()
     Input%DQF_139um =>  null()
     Input%DQF_142um =>  null()
     Input%Bt_038um =>  null()
     Input%Bt_067um =>  null()
     Input%Bt_085um =>  null()
     Input%Bt_097um =>  null()
     Input%Bt_104um =>  null()
     Input%Bt_110um =>  null()
     Input%Bt_120um =>  null()
     Input%Bt_133um =>  null()
     Input%Rad_038um =>  null()
     Input%Rad_067um =>  null()
     Input%Rad_085um =>  null()
     Input%Rad_097um =>  null()
     Input%Rad_104um =>  null()
     Input%Rad_110um =>  null()
     Input%Rad_120um =>  null()
     Input%Rad_133um =>  null()
     Input%Rad_136um =>  null()
     Input%Rad_139um =>  null()
     Input%Rad_142um =>  null()
     Input%Cosine_Zenith_Angle =>  null()
     Input%Sensor_Zenith_Angle =>  null()
     Input%Sensor_Azimuth_Angle =>  null()
     Input%Surface_Temperature => null()
     Input%Surface_Air_Temperature =>  null()
     Input%Tropopause_Temperature =>  null()
     Input%Tropopause_Height =>  null()
     Input%Tropopause_Pressure =>  null()
     Input%Surface_Pressure =>  null()
     Input%Surface_Elevation =>  null()
     Input%Latitude =>  null()
     Input%Longitude =>  null()
     Input%Rad_Clear_038um =>  null()
     Input%Rad_Clear_062um =>  null()
     Input%Rad_Clear_067um =>  null()
     Input%Rad_Clear_073um =>  null()
     Input%Rad_Clear_085um =>  null()
     Input%Rad_Clear_097um =>  null()
     Input%Rad_Clear_104um =>  null()
     Input%Rad_Clear_110um =>  null()
     Input%Rad_Clear_120um =>  null()
     Input%Rad_Clear_133um =>  null()
     Input%Rad_Clear_136um =>  null()
     Input%Rad_Clear_139um =>  null()
     Input%Rad_Clear_142um =>  null()
     Input%Surface_Emissivity_038um =>  null()
     Input%Surface_Emissivity_062um =>  null()
     Input%Surface_Emissivity_067um =>  null()
     Input%Surface_Emissivity_073um =>  null()
     Input%Surface_Emissivity_085um =>  null()
     Input%Surface_Emissivity_097um =>  null()
     Input%Surface_Emissivity_104um =>  null()
     Input%Surface_Emissivity_110um =>  null()
     Input%Surface_Emissivity_120um =>  null()
     Input%Surface_Emissivity_133um =>  null()
     Input%Surface_Emissivity_136um =>  null()
     Input%Surface_Emissivity_139um =>  null()
     Input%Surface_Emissivity_142um =>  null()
     Input%Snow_Class =>  null()
     Input%Surface_Type =>  null()
     Input%Cloud_Mask =>  null()
     Input%Cloud_Probability => null()
     Input%Ice_Cloud_Probability => null()
     Input%Cloud_Phase_Uncertainty => null()
     Input%Cloud_Type =>  null()
     Input%Elem_Idx_Nwp =>   null()
     Input%Line_Idx_Nwp =>  null()
     Input%Elem_Idx_Opposite_Corner_NWP =>  null()
     Input%Line_Idx_Opposite_Corner_NWP =>  null()
     Input%Viewing_Zenith_Angle_Idx_Rtm =>  null()
     Input%Latitude_Interp_Weight_NWP =>  null()
     Input%Longitude_Interp_Weight_NWP =>  null()
     Input%Elem_Idx_LRC_Input =>  null()
     Input%Line_Idx_LRC_Input =>   null()
     Input%Tc_Cirrus_Sounder =>   null()
     Input%Tc_NWP =>   null()
 end subroutine NULL_INPUT
 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding input 
 !-----------------------------------------------------------------------------
 subroutine NULL_DIAG()
     Diag%Array_1 =>  null()
     Diag%Array_2 =>  null()
     Diag%Array_3 =>  null()
 end subroutine NULL_DIAG
 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding output to ACHA
 !-----------------------------------------------------------------------------
 subroutine NULL_OUTPUT()
     Output%Latitude_Pc =>  null()
     Output%Longitude_Pc =>  null()
     Output%Tc_Ap =>  null()
     Output%Tc =>  null()
     Output%Ec =>  null()
     Output%Beta =>  null()
     Output%Pc =>  null()
     Output%Zc =>  null()
     Output%Zc_Base =>  null()
     Output%Pc_Eff =>  null()
     Output%Tc_Eff =>  null()
     Output%Zc_Eff =>  null()
     Output%Tau =>  null()
     Output%Tau_Uncertainty =>  null()
     Output%Reff =>  null()
     Output%Ice_Probability =>  null()
     Output%Tc_Uncertainty =>  null()
     Output%Ec_Uncertainty =>  null()
     Output%Beta_Uncertainty =>  null()
     Output%Pc_Uncertainty =>  null()
     Output%Zc_Uncertainty =>  null()
     Output%Ice_Probability_Uncertainty =>  null()
     Output%Lower_Pc =>  null()
     Output%Lower_Tc =>  null()
     Output%Lower_Zc =>  null()
     Output%Lower_Tc_Ap =>  null()
     Output%Conv_Test =>  null()
     Output%Cost =>  null()
     Output%Goodness =>  null()
     Output%Qf =>  null()
     Output%OE_Qf =>  null()
     Output%Packed_Qf =>  null()
     Output%Packed_Meta_Data =>  null()
     Output%Processing_Order  =>  null()
     Output%Inversion_Flag  =>  null()
     Output%Ec_038um =>  null()
     Output%Ec_062um =>  null()
     Output%Ec_067um =>  null()
     Output%Ec_073um =>  null()
     Output%Ec_085um =>  null()
     Output%Ec_097um =>  null()
     Output%Ec_104um =>  null()
     Output%Ec_110um =>  null()
     Output%Ec_120um =>  null()
     Output%Ec_133um =>  null()
     Output%Ec_136um =>  null()
     Output%Ec_139um =>  null()
     Output%Ec_142um =>  null()
     Output%Cloud_Type =>  null()
 end subroutine NULL_OUTPUT
 !-----------------------------------------------------------------------------
 ! Copy needed Symbol elements
 !-----------------------------------------------------------------------------
 subroutine SET_SYMBOL()
   Symbol%CLOUDY = sym%CLOUDY
   Symbol%PROB_CLOUDY = sym%PROB_CLOUDY
   Symbol%PROB_CLEAR = sym%PROB_CLEAR
   Symbol%CLEAR = sym%CLEAR

   Symbol%NO = sym%NO
   Symbol%YES = sym%YES

   Symbol%WATER_SFC = sym%WATER_SFC
   Symbol%EVERGREEN_NEEDLE_SFC = sym%EVERGREEN_NEEDLE_SFC
   Symbol%EVERGREEN_BROAD_SFC = sym%EVERGREEN_BROAD_SFC
   Symbol%DECIDUOUS_NEEDLE_SFC = sym%DECIDUOUS_NEEDLE_SFC
   Symbol%DECIDUOUS_BROAD_SFC = sym%DECIDUOUS_BROAD_SFC
   Symbol%MIXED_FORESTS_SFC = sym%MIXED_FORESTS_SFC
   Symbol%WOODLandS_SFC = sym%WOODLandS_SFC
   Symbol%WOODED_GRASS_SFC = sym%WOODED_GRASS_SFC
   Symbol%CLOSED_SHRUBS_SFC = sym%CLOSED_SHRUBS_SFC
   Symbol%OPEN_SHRUBS_SFC = sym%OPEN_SHRUBS_SFC
   Symbol%GRASSES_SFC = sym%GRASSES_SFC
   Symbol%CROPLandS_SFC = sym%CROPLandS_SFC
   Symbol%BARE_SFC = sym%BARE_SFC
   Symbol%URBAN_SFC = sym%URBAN_SFC

   Symbol%SHALLOW_OCEAN = sym%SHALLOW_OCEAN
   Symbol%Land = sym%Land
   Symbol%COASTLINE = sym%COASTLINE
   Symbol%SHALLOW_INLand_WATER = sym%SHALLOW_INLand_WATER
   Symbol%EPHEMERAL_WATER = sym%EPHEMERAL_WATER
   Symbol%DEEP_INLand_WATER = sym%DEEP_INLand_WATER
   Symbol%MODERATE_OCEAN = sym%MODERATE_OCEAN
   Symbol%DEEP_OCEAN = sym%DEEP_OCEAN

   Symbol%NO_SNOW = sym%NO_SNOW
   Symbol%SEA_ICE = sym%SEA_ICE
   Symbol%SNOW = sym%SNOW

   Symbol%CLEAR_TYPE = sym%CLEAR_TYPE
   Symbol%PROB_CLEAR_TYPE = sym%PROB_CLEAR_TYPE
   Symbol%FOG_TYPE = sym%FOG_TYPE
   Symbol%WATER_TYPE = sym%WATER_TYPE
   Symbol%SUPERCOOLED_TYPE = sym%SUPERCOOLED_TYPE
   Symbol%MIXED_TYPE = sym%MIXED_TYPE
   Symbol%OPAQUE_ICE_TYPE = sym%OPAQUE_ICE_TYPE
   Symbol%TICE_TYPE = sym%TICE_TYPE
   Symbol%CIRRUS_TYPE = sym%CIRRUS_TYPE
   Symbol%OVERLAP_TYPE = sym%OVERLAP_TYPE
   Symbol%OVERSHOOTING_TYPE = sym%OVERSHOOTING_TYPE
   Symbol%UNKNOWN_TYPE = sym%UNKNOWN_TYPE
   Symbol%DUST_TYPE = sym%DUST_TYPE
   Symbol%SMOKE_TYPE = sym%SMOKE_TYPE
   Symbol%FIRE_TYPE = sym%FIRE_TYPE

   Symbol%CLEAR_PHASE = sym%CLEAR_PHASE
   Symbol%WATER_PHASE = sym%WATER_PHASE
   Symbol%SUPERCOOLED_PHASE = sym%SUPERCOOLED_PHASE
   Symbol%MIXED_PHASE = sym%MIXED_PHASE
   Symbol%ICE_PHASE = sym%ICE_PHASE
   Symbol%UNKNOWN_PHASE = sym%UNKNOWN_PHASE
 end subroutine SET_SYMBOL

 subroutine SET_OUTPUT()
   Output%Latitude_Pc => Nav%Lat_Pc
   Output%Longitude_Pc => Nav%Lon_Pc
   Output%Tc => ACHA%Tc
   Output%Tc_Ap => ACHA%Tc_Ap
   Output%Ec => ACHA%Ec
   Output%Beta => ACHA%Beta
   Output%Pc => ACHA%Pc
   Output%Zc => ACHA%Zc
   Output%Pc_Eff => ACHA%Pc_Eff
   Output%Tc_Eff => ACHA%Tc_Eff
   Output%Zc_Eff => ACHA%Zc_Eff
   Output%Zc_Base => ACHA%Zc_Base
   Output%Tau => ACHA%Tau
   Output%Tau_Uncertainty => ACHA%Tau_Uncer
   Output%Reff => ACHA%Reff
   Output%Ice_Probability => ACHA%Ice_Probability
   Output%Ice_Probability_Uncertainty => ACHA%Ice_Probability_Uncertainty
   Output%Tc_Uncertainty => ACHA%Tc_Uncertainty
   Output%Ec_Uncertainty => ACHA%Ec_Uncertainty
   Output%Beta_Uncertainty => ACHA%Beta_Uncertainty
   Output%Pc_Uncertainty => ACHA%Pc_Uncertainty
   Output%Zc_Uncertainty => ACHA%Zc_Uncertainty
   Output%Lower_Tc_Uncertainty => ACHA%Lower_Tc_Uncertainty
   Output%Lower_Zc_Uncertainty => ACHA%Lower_Zc_Uncertainty
   Output%Lower_Pc_Uncertainty => ACHA%Lower_Pc_Uncertainty
   Output%Lower_Pc => ACHA%Lower_Pc
   Output%Lower_Tc => ACHA%Lower_Tc
   Output%Lower_Zc => ACHA%Lower_Zc
   Output%Lower_Tc_Ap => ACHA%Lower_Tc_Ap
   Output%Conv_Test  => ACHA%Conv_Test
   Output%Cost  => ACHA%Cost
   Output%Goodness  => ACHA%Goodness
   Output%Qf => ACHA%Quality_Flag
   Output%OE_Qf => ACHA%OE_Quality_Flags
   Output%Packed_Qf => ACHA%Packed_Quality_Flags
   Output%Packed_Meta_Data => ACHA%Packed_Meta_Data_Flags
   Output%Processing_Order  => ACHA%Processing_Order
   Output%Inversion_Flag  => ACHA%Inversion_Flag
   Output%Ec_038um => ACHA%Ec_375um
   Output%Ec_062um => ACHA%Ec_62um
   Output%Ec_067um => ACHA%Ec_67um
   Output%Ec_073um => ACHA%Ec_73um
   Output%Ec_085um => ACHA%Ec_85um
   Output%Ec_097um => ACHA%Ec_97um
   Output%Ec_104um => ACHA%Ec_104um
   Output%Ec_110um => ACHA%Ec_11um
   Output%Ec_120um => ACHA%Ec_12um
   Output%Ec_133um => ACHA%Ec_133um
   Output%Ec_136um => ACHA%Ec_136um
   Output%Ec_139um => ACHA%Ec_139um
   Output%Ec_142um => ACHA%Ec_142um
   Output%Cloud_Type => ACHA%Cld_Type
 end subroutine SET_OUTPUT
!--------------------------------------------------------
 subroutine SET_INPUT()

   Input%ACHA_Mode_Flag_In = ACHA%Mode
   Input%Number_of_Elements = Image%Number_Of_Elements
   Input%Number_of_Lines = Image%Number_Of_Lines_Per_Segment
   Input%Smooth_Nwp_Fields_Flag = NWP_PIX%Smooth_Nwp_Flag
   Input%Process_Undetected_Cloud_Flag = Process_Undetected_Cloud_Flag
   Input%Sensor_Resolution_KM = Sensor%Spatial_Resolution_Meters/1000.0
   Input%WMO_Id = Sensor%WMO_Id
   Input%Chan_Idx_038um = 20      !channel number for 3.75
   Input%Chan_Idx_062um = 37      !channel number for 6.2
   Input%Chan_Idx_067um = 27      !channel number for 6.7
   Input%Chan_Idx_073um = 28      !channel number for 7.3
   Input%Chan_Idx_085um = 29      !channel number for 8.5
   Input%Chan_Idx_097um = 30      !channel number for 9.7
   Input%Chan_Idx_104um = 38      !channel number for 10.4
   Input%Chan_Idx_110um = 31      !channel number for 11
   Input%Chan_Idx_120um = 32      !channel number for 12
   Input%Chan_Idx_133um = 33      !channel number for 13.3
   Input%Chan_Idx_136um = 34      !channel number for 13.6
   Input%Chan_Idx_139um = 35      !channel number for 13.9
   Input%Chan_Idx_142um = 36      !channel number for 14.2
   Input%Chan_On_038um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_038um)
   Input%Chan_On_062um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_062um)
   Input%Chan_On_067um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_067um)
   Input%Chan_On_073um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_073um)
   Input%Chan_On_085um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_085um)
   Input%Chan_On_097um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_097um)
   Input%Chan_On_104um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_104um)
   Input%Chan_On_110um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_110um)
   Input%Chan_On_120um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_120um)
   Input%Chan_On_133um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_133um)
   Input%Chan_On_136um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_136um)
   Input%Chan_On_139um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_139um)
   Input%Chan_On_142um = Sensor%Chan_On_Flag_Default(Input%Chan_Idx_142um)
   Input%Invalid_Data_Mask => Bad_Pixel_Mask
   Input%DQF_038um => ch(Input%Chan_Idx_038um)%DQF
   Input%DQF_062um => ch(Input%Chan_Idx_062um)%DQF
   Input%DQF_067um => ch(Input%Chan_Idx_067um)%DQF
   Input%DQF_073um => ch(Input%Chan_Idx_073um)%DQF
   Input%DQF_085um => ch(Input%Chan_Idx_085um)%DQF
   Input%DQF_097um => ch(Input%Chan_Idx_097um)%DQF
   Input%DQF_104um => ch(Input%Chan_Idx_104um)%DQF
   Input%DQF_110um => ch(Input%Chan_Idx_110um)%DQF
   Input%DQF_120um => ch(Input%Chan_Idx_120um)%DQF
   Input%DQF_133um => ch(Input%Chan_Idx_133um)%DQF
   Input%DQF_136um => ch(Input%Chan_Idx_136um)%DQF
   Input%DQF_139um => ch(Input%Chan_Idx_139um)%DQF
   Input%DQF_142um => ch(Input%Chan_Idx_142um)%DQF
   Input%Bt_038um => ch(Input%Chan_Idx_038um)%Bt_Toa
   Input%Bt_062um => ch(Input%Chan_Idx_062um)%Bt_Toa
   Input%Bt_067um => ch(Input%Chan_Idx_067um)%Bt_Toa
   Input%Bt_073um => ch(Input%Chan_Idx_073um)%Bt_Toa
   Input%Bt_085um => ch(Input%Chan_Idx_085um)%Bt_Toa
   Input%Bt_097um => ch(Input%Chan_Idx_097um)%Bt_Toa
   Input%Bt_104um => ch(Input%Chan_Idx_104um)%Bt_Toa
   Input%Bt_110um => ch(Input%Chan_Idx_110um)%Bt_Toa
   Input%Bt_120um => ch(Input%Chan_Idx_120um)%Bt_Toa
   Input%Bt_133um => ch(Input%Chan_Idx_133um)%Bt_Toa
   Input%Bt_136um => ch(Input%Chan_Idx_136um)%Bt_Toa
   Input%Bt_139um => ch(Input%Chan_Idx_139um)%Bt_Toa
   Input%Bt_142um => ch(Input%Chan_Idx_142um)%Bt_Toa
   Input%Rad_038um => ch(Input%Chan_Idx_038um)%Rad_Toa
   Input%Rad_062um => ch(Input%Chan_Idx_062um)%Rad_Toa
   Input%Rad_067um => ch(Input%Chan_Idx_067um)%Rad_Toa
   Input%Rad_073um => ch(Input%Chan_Idx_073um)%Rad_Toa
   Input%Rad_085um => ch(Input%Chan_Idx_085um)%Rad_Toa
   Input%Rad_097um => ch(Input%Chan_Idx_097um)%Rad_Toa
   Input%Rad_104um => ch(Input%Chan_Idx_104um)%Rad_Toa
   Input%Rad_110um => ch(Input%Chan_Idx_110um)%Rad_Toa
   Input%Rad_120um => ch(Input%Chan_Idx_120um)%Rad_Toa
   Input%Rad_133um => ch(Input%Chan_Idx_133um)%Rad_Toa
   Input%Rad_136um => ch(Input%Chan_Idx_136um)%Rad_Toa
   Input%Rad_139um => ch(Input%Chan_Idx_139um)%Rad_Toa
   Input%Rad_142um => ch(Input%Chan_Idx_142um)%Rad_Toa
   Input%Cosine_Zenith_Angle => Geo%Coszen
   Input%Sensor_Zenith_Angle => Geo%Satzen
   Input%Sensor_Azimuth_Angle => Geo%Sataz
   Input%Surface_Temperature => NWP_PIX%Tsfc
   Input%Surface_Air_Temperature =>  NWP_PIX%Tair
   Input%Tropopause_Temperature =>  NWP_PIX%Ttropo
   Input%Tropopause_Height =>  NWP_PIX%Ztropo
   Input%Tropopause_Pressure =>  NWP_PIX%Ptropo
   Input%Surface_Pressure =>  NWP_PIX%Psfc
   Input%Surface_Elevation => Sfc%Zsfc
   Input%Latitude => Nav%Lat
   Input%Longitude => Nav%Lon
   Input%Rad_Clear_038um => ch(Input%Chan_Idx_038um)%Rad_Toa_Clear
   Input%Rad_Clear_062um => ch(Input%Chan_Idx_062um)%Rad_Toa_Clear
   Input%Rad_Clear_067um => ch(Input%Chan_Idx_067um)%Rad_Toa_Clear
   Input%Rad_Clear_073um => ch(Input%Chan_Idx_073um)%Rad_Toa_Clear
   Input%Rad_Clear_085um => ch(Input%Chan_Idx_085um)%Rad_Toa_Clear
   Input%Rad_Clear_097um => ch(Input%Chan_Idx_097um)%Rad_Toa_Clear
   Input%Rad_Clear_104um => ch(Input%Chan_Idx_104um)%Rad_Toa_Clear
   Input%Rad_Clear_110um => ch(Input%Chan_Idx_110um)%Rad_Toa_Clear
   Input%Rad_Clear_120um => ch(Input%Chan_Idx_120um)%Rad_Toa_Clear
   Input%Rad_Clear_133um => ch(Input%Chan_Idx_133um)%Rad_Toa_Clear
   Input%Rad_Clear_136um => ch(Input%Chan_Idx_136um)%Rad_Toa_Clear
   Input%Rad_Clear_139um => ch(Input%Chan_Idx_139um)%Rad_Toa_Clear
   Input%Rad_Clear_142um => ch(Input%Chan_Idx_142um)%Rad_Toa_Clear
   Input%Surface_Emissivity_038um => ch(Input%Chan_Idx_038um)%Sfc_Emiss
   Input%Surface_Emissivity_062um => ch(Input%Chan_Idx_062um)%Sfc_Emiss
   Input%Surface_Emissivity_067um => ch(Input%Chan_Idx_067um)%Sfc_Emiss
   Input%Surface_Emissivity_073um => ch(Input%Chan_Idx_073um)%Sfc_Emiss
   Input%Surface_Emissivity_085um => ch(Input%Chan_Idx_085um)%Sfc_Emiss
   Input%Surface_Emissivity_097um => ch(Input%Chan_Idx_097um)%Sfc_Emiss
   Input%Surface_Emissivity_104um => ch(Input%Chan_Idx_104um)%Sfc_Emiss
   Input%Surface_Emissivity_110um => ch(Input%Chan_Idx_110um)%Sfc_Emiss
   Input%Surface_Emissivity_120um => ch(Input%Chan_Idx_120um)%Sfc_Emiss
   Input%Surface_Emissivity_133um => ch(Input%Chan_Idx_133um)%Sfc_Emiss
   Input%Surface_Emissivity_136um => ch(Input%Chan_Idx_136um)%Sfc_Emiss
   Input%Surface_Emissivity_139um => ch(Input%Chan_Idx_139um)%Sfc_Emiss
   Input%Surface_Emissivity_142um => ch(Input%Chan_Idx_142um)%Sfc_Emiss
   Input%Snow_Class => Sfc%Snow
   Input%Surface_Type => Sfc%Sfc_Type
   Input%Cloud_Mask => CLDMASK%Cld_Mask
   Input%Cloud_Probability => CLDMASK%Posterior_Cld_Probability
   Input%Ice_Cloud_Probability =>CLDMASK%Posterior_Ice_Probability 
   Input%Cloud_Phase_Uncertainty => Cld_Phase_Uncertainty
   Input%Cloud_Type => Cld_Type
   Input%Elem_Idx_Nwp =>  NWP_PIX%I_Nwp
   Input%Line_Idx_Nwp => NWP_PIX%J_Nwp
   Input%Elem_Idx_Opposite_Corner_NWP => NWP_PIX%I_Nwp_x
   Input%Line_Idx_Opposite_Corner_NWP => NWP_PIX%J_Nwp_x
   Input%Viewing_Zenith_Angle_Idx_Rtm => Zen_Idx_Rtm
   Input%Latitude_Interp_Weight_NWP => NWP_PIX%Lat_Nwp_Fac
   Input%Longitude_Interp_Weight_NWP => NWP_PIX%Lon_Nwp_Fac
   Input%Elem_Idx_LRC_Input => I_LRC
   Input%Line_Idx_LRC_Input =>  J_LRC
   Input%Tc_Cirrus_Sounder =>  Tc_Cirrus_Background
   Input%Tc_Nwp =>  NWP_PIX%Tc
   Input%Tc_Opaque =>  Tc_Opaque_Cloud
   Input%Pc_Opaque =>  Pc_Opaque_Cloud
   Input%Zc_Opaque =>  Zc_Opaque_Cloud

 end subroutine SET_INPUT
!----------------------------------------------------------------------
! Set up the ACHA Diagnostic Structure
!----------------------------------------------------------------------
 subroutine SET_DIAG
     Diag%Array_1 => Diag_Pix_Array_1 
     Diag%Array_2 => Diag_Pix_Array_2 
     Diag%Array_3 => Diag_Pix_Array_3 
 end subroutine SET_DIAG
 !----------------------------------------------------------------------
 ! Set up the ACHA Diagnostic (one pixel) Dump Structure
 !----------------------------------------------------------------------
 subroutine SET_DUMP
     Dump%Segment_Number = Image%Segment_Number
     Dump%Number_Lines_Per_Segment =  Image%Number_Of_Lines_Per_Segment
     Dump%Elem_Abs_Idx = Elem_Abs_Idx_Acha_Dump
     Dump%Line_Abs_Idx = Line_Abs_Idx_Acha_Dump
 end subroutine SET_DUMP
 !----------------------------------------------------------------------
 ! 
 !----------------------------------------------------------------------
 subroutine LHP_CHN_CHECK ()
   REAL, DIMENSION(Nchan_Clavrx) :: LHP_THRESH
   REAL, DIMENSION(Nchan_Clavrx) :: LHP_THRESH_INIT
   REAL :: MFPT_104, MFPT_110
   INTEGER :: N_chan
   CHARACTER(3000) :: Algo_Thresh_File
   REAL, PARAMETER:: SOUNDER_LHP_THRESH = 999.0  !  K

   !--- Initalize LHP THRESHOLD
   LHP_THRESH (:) = MISSING_VALUE_REAL4
   LHP_THRESH_INIT (:) = MISSING_VALUE_REAL4

   !--- Initialize 10.3 and 11um FPT local variables
   MFPT_104 = MISSING_VALUE_REAL4
   MFPT_110 = MISSING_VALUE_REAL4

   !--- Initialize Use_104um_Flag
   ABI_Use_104um_Flag = .FALSE.

   !--- Set local variables for easier reading.
   MFPT_104 = Ch(38)%Max_Focal_Plane_Temp
   MFPT_110 = Ch(31)%Max_Focal_Plane_Temp

   !--- Initialize LHP to default thresholds
   LHP_THRESH_INIT(20) = ABI_FPT_Thresh_038um
   LHP_THRESH_INIT(37) = ABI_FPT_Thresh_062um
   LHP_THRESH_INIT(27) = ABI_FPT_Thresh_067um
   LHP_THRESH_INIT(28) = ABI_FPT_Thresh_073um
   LHP_THRESH_INIT(29) = ABI_FPT_Thresh_085um
   LHP_THRESH_INIT(30) = ABI_FPT_Thresh_097um
   LHP_THRESH_INIT(38) = ABI_FPT_Thresh_104um
   LHP_THRESH_INIT(31) = ABI_FPT_Thresh_110um
   LHP_THRESH_INIT(32) = ABI_FPT_Thresh_120um
   LHP_THRESH_INIT(33) = ABI_FPT_Thresh_133um

   Algo_Thresh_File =trim(Ancil_Data_Dir)// &
                      "static/algo_lhp_thresh/"//&
                      "goes17_thermal_limits_acha.txt"

   !--- Read in LHP Thresholds
   if (Image%Segment_Number == 1 ) then
     print*, "Opening ", trim(Algo_Thresh_File)
   endif
   CALL READ_LHP_THRESH_FILE(TRIM(Algo_Thresh_File), LHP_THRESH)

   !--- Check if LHP_THRESH is missing
   do N_Chan = 1, Nchan_Clavrx
     if (LHP_THRESH(N_Chan) == MISSING_VALUE_REAL4) &
        LHP_THRESH(N_Chan) = LHP_THRESH_INIT(N_Chan)
   end do

   !--- Nullify LHP_THRESH is sounder (aka fusion) data is used
   do N_Chan = 1, Nchan_Clavrx
      if (Sensor%Chan_On_Flag_Default(N_Chan) == sym%NO) cycle
      if (maxval(Ch(N_Chan)%Source) == 1) then
         LHP_THRESH(N_Chan) = SOUNDER_LHP_THRESH
      endif
   enddo

   !--- For now limit to ABI channels
   Input%Chan_On_038um = LHP_LOCAL_CHAN_ON(20, LHP_THRESH(20))
   Input%Chan_On_062um = LHP_LOCAL_CHAN_ON(37, LHP_THRESH(37))
   Input%Chan_On_067um = LHP_LOCAL_CHAN_ON(27, LHP_THRESH(27))
   Input%Chan_On_073um = LHP_LOCAL_CHAN_ON(28, LHP_THRESH(28))
   Input%Chan_On_085um = LHP_LOCAL_CHAN_ON(29, LHP_THRESH(29))
   Input%Chan_On_097um = LHP_LOCAL_CHAN_ON(30, LHP_THRESH(30))
   Input%Chan_On_104um = LHP_LOCAL_CHAN_ON(38, LHP_THRESH(38))
   Input%Chan_On_110um = LHP_LOCAL_CHAN_ON(31, LHP_THRESH(31))
   Input%Chan_On_120um = LHP_LOCAL_CHAN_ON(32, LHP_THRESH(32))
   Input%Chan_On_133um = LHP_LOCAL_CHAN_ON(33, LHP_THRESH(33))

   !--- 10.4um check
   if ((MFPT_104 < LHP_THRESH(38)) .and. &
      (MFPT_110 > LHP_THRESH(31))) then
      ABI_Use_104um_Flag = .true.
   endif

 END subroutine LHP_CHN_CHECK

end module ACHA_CLAVRX_BRIDGE
