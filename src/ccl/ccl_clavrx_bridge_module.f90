!$Id: acha_clavrx_bridge_module.f90 1789 2016-09-28 22:20:51Z heidinger $
!------------------------------------------------------------------------------
!  NOAA AWG Cloud Height Algorithm (ACHA) Bridge Code
!
!  This module houses the routines that serve as a bridge between
!  the CLAVR-x processing system and the ACHA code.
!
!------------------------------------------------------------------------------
module CCL_CLAVRX_BRIDGE

 use CCL_SERVICES_MOD
 use CCL_MODULE
 use CLAVRX_MESSAGE_MOD, only: MESG
   
 implicit none

 public :: CCL_BRIDGE
 private:: SET_SYMBOL, SET_INPUT, SET_OUTPUT, NULL_INPUT, NULL_OUTPUT
 private:: SET_DIAG, NULL_DIAG

 !--------------------------------------------------------------------
 ! define structures that will be arguments to ACHA
 !--------------------------------------------------------------------
 type(ccl_symbol_struct), private :: Symbol
 type(ccl_input_struct), private :: Input
 type(ccl_output_struct), private :: Output
 type(ccl_diag_struct), private :: Diag

 contains

!----------------------------------------------------------------------
! CCL BRIDGE SUBROUTINE
!---------------------------------------------------------------------- 
 subroutine CCL_BRIDGE()
 
   implicit none

   logical, save:: First_Call = .true.

   if (CCL%MODE == 0) return

   if (First_Call .eqv. .true.) then
       call MESG('CCL starts ', color = 46)
   endif

   !---null pointers before filling them
   call NULL_INPUT()
   call NULL_OUTPUT()

   !---------------------------------------------------------------
   ! allocate layer fractions here since they depend on CCL_Type
   !---------------------------------------------------------------

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

   if (Image%Segment_Number ==1) then 
     call SETUP_CCL(Input%CCL_Type)
   endif

   !--- cloud cover layers
   call COMPUTE_CLOUD_COVER_LAYERS(Input, Symbol, Output) !, Diag)

   if (Image%Segment_Number == Image%Number_Of_Segments) then 
      call DESTROY_CCL()
   endif

   !--- copy output into CLAVR-x variables

!NEEDED?
!  CCL%Cloud_Fraction = Output%Total_Cloud_Fraction
!  CCL%Cloud_Fraction_Uncer = Output%Total_Cloud_Fraction_Uncer
!  CCL%High_Cloud_Fraction = Output%High_Cloud_Fraction
!  CCL%Mid_Cloud_Fraction = Output%Mid_Cloud_Fraction
!  CCL%Low_Cloud_Fraction = Output%Low_Cloud_Fraction
!  CCL%Cld_Layer = Output%Cloud_Layer

   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_INPUT()
   call NULL_OUTPUT()

   First_Call = .false.

 end subroutine CCL_BRIDGE

 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding input 
 !-----------------------------------------------------------------------------
 subroutine NULL_INPUT()
     Input%Invalid_Data_Mask =>  null()
     Input%Latitude =>  null()
     Input%Longitude =>  null()
     Input%Surface_Type =>  null()
     Input%Cloud_Mask =>  null()
     Input%Cloud_Probability => null()
     Input%Cloud_Type =>  null()
     Input%Pc =>  null()
     Input%Pc_Base =>  null()
!ynoh (cira/csu) for ccl mode 3
     Input%Pc_Lower_Base =>  null()
     Input%Pc_Lower =>  null()
     Input%Alt =>  null()
     Input%Alt_Base =>  null()
     Input%Alt_Lower =>  null()
     Input%Alt_Lower_Base =>  null()
     Input%Elem_Idx_Nwp =>   null()
     Input%Line_Idx_Nwp =>  null()
     Input%Elem_Idx_Opposite_Corner_NWP =>  null()
     Input%Line_Idx_Opposite_Corner_NWP =>  null()
     Input%Viewing_Zenith_Angle_Idx_Rtm =>  null()
     Input%Latitude_Interp_Weight_NWP =>  null()
     Input%Longitude_Interp_Weight_NWP =>  null()
     Input%ACHA_QF => null()
     Input%Freezing_Level_Pressure => null()
     Input%Tc => null()
     Input%Tc_Base => null()
     Input%Tc_Lower => null()
     Input%Zc => null()
     Input%Free_Convection_Height => null()
     Input%Bt_110um => null()
     Input%Bt_067um => null()
     Input%Tsfc => null()
     Input%Emiss_Tropo_11 => null()

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
     Output%Total_Cloud_Fraction =>  null()
     Output%Total_Cloud_Fraction_Uncer =>  null()
     Output%Fraction_Layer1 => null()
     Output%Fraction_Layer2 => null()
     Output%Fraction_Layer3 => null()
     Output%Fraction_Layer4 => null()
     Output%Fraction_Layer5 => null()
     Output%Cloud_Layer =>  null()
     Output%Supercooled_Cloud_Layer => null()
     Output%Conv_Cloud_Layer => null()
     Output%Supercooled_Total_Fraction => null()
     Output%Supercooled_Fraction_Layer1 => null()
     Output%Supercooled_Fraction_Layer2 => null()
     Output%Supercooled_Fraction_Layer3 => null()
     Output%Supercooled_Fraction_Layer4 => null()
     Output%Supercooled_Fraction_Layer5 => null()
     Output%Conv_Total_Fraction => null()
     Output%Conv_Fraction_Layer1 => null()
     Output%Conv_Fraction_Layer2 => null()
     Output%Conv_Fraction_Layer3 => null()
     Output%Conv_Fraction_Layer4 => null()
     Output%Conv_Fraction_Layer5 => null()
     Output%QF => null()
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
   Symbol%WOODLANDS_SFC = sym%WOODLANDS_SFC
   Symbol%WOODED_GRASS_SFC = sym%WOODED_GRASS_SFC
   Symbol%CLOSED_SHRUBS_SFC = sym%CLOSED_SHRUBS_SFC
   Symbol%OPEN_SHRUBS_SFC = sym%OPEN_SHRUBS_SFC
   Symbol%GRASSES_SFC = sym%GRASSES_SFC
   Symbol%CROPLANDS_SFC = sym%CROPLANDS_SFC
   Symbol%BARE_SFC = sym%BARE_SFC
   Symbol%URBAN_SFC = sym%URBAN_SFC

   Symbol%SHALLOW_OCEAN = sym%SHALLOW_OCEAN
   Symbol%LAND = sym%LAND
   Symbol%COASTLINE = sym%COASTLINE
   Symbol%SHALLOW_INLAND_WATER = sym%SHALLOW_INLAND_WATER
   Symbol%EPHEMERAL_WATER = sym%EPHEMERAL_WATER
   Symbol%DEEP_INLAND_WATER = sym%DEEP_INLAND_WATER
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
   Output%Total_Cloud_Fraction => CCL%Cloud_Fraction
   Output%Total_Cloud_Fraction_Uncer => CCL%Cloud_Fraction_Uncer
   Output%Cloud_Layer => CCL%Cloud_Layer
   Output%Fraction_Layer1 => CCL%Cloud_Fraction_Layer1
   Output%Fraction_Layer2 => CCL%Cloud_Fraction_Layer2
   Output%Fraction_Layer3 => CCL%Cloud_Fraction_Layer3
   Output%Fraction_Layer4 => CCL%Cloud_Fraction_Layer4
   Output%Fraction_Layer5 => CCL%Cloud_Fraction_Layer5
   Output%Supercooled_Total_Fraction => CCL%Supercooled_Cloud_Fraction
   Output%Supercooled_Cloud_Layer => CCL%Supercooled_Cloud_Layer
   Output%Supercooled_Fraction_Layer1 => CCL%Supercooled_Cloud_Fraction_Layer1
   Output%Supercooled_Fraction_Layer2 => CCL%Supercooled_Cloud_Fraction_Layer2
   Output%Supercooled_Fraction_Layer3 => CCL%Supercooled_Cloud_Fraction_Layer3
   Output%Supercooled_Fraction_Layer4 => CCL%Supercooled_Cloud_Fraction_Layer4
   Output%Supercooled_Fraction_Layer5 => CCL%Supercooled_Cloud_Fraction_Layer5
   Output%Conv_Total_Fraction => CCL%Conv_Cloud_Fraction
   Output%Conv_Cloud_Layer => CCL%Conv_Cloud_Layer
   Output%Conv_Fraction_Layer1 => CCL%Conv_Cloud_Fraction_Layer1
   Output%Conv_Fraction_Layer2 => CCL%Conv_Cloud_Fraction_Layer2
   Output%Conv_Fraction_Layer3 => CCL%Conv_Cloud_Fraction_Layer3
   Output%Conv_Fraction_Layer4 => CCL%Conv_Cloud_Fraction_Layer4
   Output%Conv_Fraction_Layer5 => CCL%Conv_Cloud_Fraction_Layer5
   Output%QF => CCL%QF

 end subroutine SET_OUTPUT
!--------------------------------------------------------
 subroutine SET_INPUT()

   Input%CCL_Mode = CCL%Mode
   Input%CCL_Type = CCL%Type

   Input%Number_of_Elements = Image%Number_Of_Elements
   Input%Number_of_Lines = Image%Number_Of_Lines_Per_Segment
   Input%Sensor_Resolution_KM = Sensor%Spatial_Resolution_Meters/1000.0
   Input%Smooth_Nwp_Fields_Flag = NWP_PIX%Smooth_Nwp_Flag

   Input%Invalid_Data_Mask =>  Bad_Pixel_Mask
   Input%Latitude => Nav%Lat
   Input%Longitude => Nav%Lon
   Input%Surface_Type => Sfc%Sfc_Type
   Input%Cloud_Mask => CLDMASK%Cld_Mask
   Input%Cloud_Probability => CLDMASK%Posterior_Cld_Probability
   Input%Cloud_Type => Cld_Type
   Input%Pc => ACHA%Pc
   Input%Pc_Base => BASE%Pc_Base
!ynoh (cira/csu) for ccl mode 3
   Input%Pc_Lower_Base => BASE%Lower_Base_Pc
   Input%Pc_Lower => ACHA%Lower_Pc
   Input%Alt => ACHA%Alt
   Input%Alt_Base => BASE%Base_Alt
   Input%Alt_Lower => ACHA%Lower_Alt
   Input%Alt_Lower_Base => BASE%Lower_Base_Alt
   Input%Elem_Idx_Nwp =>  NWP_PIX%I_Nwp
   Input%Line_Idx_Nwp => NWP_PIX%J_Nwp
   Input%Elem_Idx_Opposite_Corner_NWP => NWP_PIX%I_Nwp_x
   Input%Line_Idx_Opposite_Corner_NWP => NWP_PIX%J_Nwp_x
   Input%Viewing_Zenith_Angle_Idx_Rtm => Zen_Idx_Rtm
   Input%Latitude_Interp_Weight_NWP => NWP_PIX%Lat_Nwp_Fac
   Input%Longitude_Interp_Weight_NWP => NWP_PIX%Lon_Nwp_Fac
   Input%ACHA_QF => ACHA%Quality_Flag
   Input%Freezing_Level_Pressure => NWP_PIX%FrzPre
   Input%Tc => ACHA%Tc
   Input%Tc_Base => BASE%Tc_Base
   Input%Tc_Lower => ACHA%Lower_Tc
   Input%Zc => ACHA%Zc
   Input%Free_Convection_Height => NWP_PIX%LFC_Height
   Input%Bt_110um => ch(31)%Bt_TOA
   Input%Bt_067um => ch(27)%Bt_TOA
   Input%Tsfc => NWP_PIX%Tsfc
   Input%Emiss_Tropo_11 => Ch(31)%Emiss_Tropo ! ACHA%Ec !Ch(31)%Emiss_Tropo
   Input%Chan_On_067um = Sensor%Chan_On_Flag_Default(27)
   Input%Chan_On_110um = Sensor%Chan_On_Flag_Default(31)

 end subroutine SET_INPUT
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
 subroutine SET_DIAG
     Diag%Array_1 => Diag_Pix_Array_1 
     Diag%Array_2 => Diag_Pix_Array_2 
     Diag%Array_3 => Diag_Pix_Array_3 
 end subroutine SET_DIAG

end module CCL_CLAVRX_BRIDGE
