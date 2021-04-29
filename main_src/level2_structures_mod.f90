! $Id: level2_structures_mod.f90 3974 2020-09-13 19:39:15Z awalther $
!--------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------
module LEVEL2_STRUCTURES_MOD

  use CONSTANTS_MOD , only: &
      real4, int2, int4, int1, real8

  implicit none

  private
  public SET_CLAVRX_GLOBAL_ATTRIBUTES

  type, public :: L2_Glob_Attr_Definition
    integer:: Hdf_File_Id
    character(len=20):: Data_Type
    character(len=1020):: File_Name
    character(len=1020):: File_1b
    integer(kind=int2):: Start_Year
    integer(kind=int2):: End_Year
    integer(kind=int2):: Start_Day
    integer(kind=int2):: End_Day
    real(kind=real4):: Start_Time
    real(kind=real4):: End_Time
    integer(kind=int4):: Num_Cells
    integer(kind=int4):: Num_Cells_With_Data
    character(len=20):: Grid_Format
    real(kind=real4):: Grid_Resolution
    integer:: Therm_Cal_1b
    integer:: Ref_Cal_1b
    integer:: Nav_Opt
    integer:: Use_Sst_Anal
    integer:: Modis_Clr_Alb_Flag
    integer:: Nwp_Opt

    real(kind=real4):: resolution_km
    real(kind=real4):: dlat
    real(kind=real4):: ch1_gain_low
    real(kind=real4):: ch1_gain_high
    real(kind=real4):: ch2_gain_low
    real(kind=real4):: ch2_gain_high
    real(kind=real4):: ch3a_gain_low
    real(kind=real4):: ch3a_gain_high
    real(kind=real4):: ch1_switch_count
    real(kind=real4):: ch1_dark_count

    real(kind=real4):: ch2_switch_count
    real(kind=real4):: ch2_dark_count
    real(kind=real4):: ch3a_switch_count
    real(kind=real4):: ch3a_dark_count
    real(kind=real4):: sun_earth_distance
    real(kind=real4):: c1
    real(kind=real4):: c2
    real(kind=real4):: a_20
    real(kind=real4):: b_20
    real(kind=real4):: nu_20
    real(kind=real4):: a_31
    real(kind=real4):: b_31
    real(kind=real4):: nu_31
    real(kind=real4):: a_32
    real(kind=real4):: b_32
    real(kind=real4):: nu_32
    real(kind=real4):: solar_ch20_nu
    real(kind=real4):: geo_sub_lon
    real(kind=real4):: geo_sub_lat
    real(kind=real4):: timerr_seconds
    character(len=10):: mask_mode  
    character(len=50):: acha_mode  
    character(len=50):: acha_mode_user
    integer(kind=int4):: dcomp_mode
    integer(kind=int4):: wmo_sc_code
    character(len=20):: platform_name
    character(len=20):: sensor_name
    character(len=120):: dark_name
    character(len=120):: mask_name
    integer(kind=int4):: subset_flag
    real(kind=real4):: lat_south_subset
    real(kind=real4):: lat_north_subset
    real(kind=real4):: lon_west_subset
    real(kind=real4):: lon_east_subset
    character(len=20):: subset_name
    character(len=36):: creator
    character(len=36):: plang
    character(len=36):: timestamp
    character(len=80) :: hdf_ver
    character(len=36) :: machine
    real(kind=real4)::  LWIR_Focal_Plane_Temperature
  end type L2_Glob_Attr_Definition

!level2_mod
  type, public :: Sds_Struct
    character(len=100):: Sds_Name
    logical:: On_Flag
    integer(kind=int4):: Sds_Idx
    integer(kind=int4):: Input_Data_Type_HDF
    integer(kind=int4):: Level2_Data_Type_HDF
    integer(kind=int4):: Level2_Data_Type_NETCDF
    integer(kind=int4):: Rank
    integer(kind=int4):: Data_Type
    integer(kind=int1):: Scaling_Type
    character(len=100):: Units
    character(len=300):: Long_Name
    character(len=100):: Standard_Name
    character(len=1020):: Flags_String
    integer:: Number_Of_Flags
    real(kind=real4):: Actual_Missing
    real(kind=real4):: Scale_Factor
    real(kind=real4):: Add_Offset
    real(kind=real4), dimension(2):: Actual_Range
    real(kind=real4):: Fill_Value
    real(kind=real4), dimension(2):: Valid_Range
    integer(kind=int4):: Sds_Id
    !-- pointers to hold all possible data types
    integer(kind=int1), dimension(:), pointer:: Sds_Data_1d_I1
    integer(kind=int4), dimension(:), pointer:: Sds_Data_1d_I4
    integer(kind=int1), dimension(:,:), pointer:: Sds_Data_2d_I1
    integer(kind=int2), dimension(:,:), pointer:: Sds_Data_2d_I2
    integer(kind=int4), dimension(:,:), pointer:: Sds_Data_2d_I4
    integer(kind=int1), dimension(:,:,:), pointer:: Sds_Data_3d_I1
    real, dimension(:), pointer:: Sds_Data_1d_R4
    real, dimension(:,:), pointer:: Sds_Data_2d_R4
  end type Sds_Struct

  type (L2_Glob_Attr_Definition), public:: Clavrx_Global_Attr

contains
!-------------------------------------------------------------------------
! Set the Global Attribute Structure
!-------------------------------------------------------------------------
 subroutine SET_CLAVRX_GLOBAL_ATTRIBUTES(data_type,file_name,file_1b, &
                           resolution_km, &
                           start_year,end_year,start_day,end_day,start_time,end_time,&
                           num_cells,num_cells_with_data,grid_format,dlat, &
                           therm_cal_1b,Ref_cal_1b,nav_opt,use_sst_anal, &
                           modis_clr_alb_flag, nwp_opt, ch1_gain_low, ch1_gain_high, &
                           ch1_switch_count, ch1_dark_count, &
                           ch2_gain_low, ch2_gain_high, &
                           ch2_switch_count, ch2_dark_count, &
                           ch3a_gain_low, ch3a_gain_high, &
                           ch3a_switch_count, ch3a_dark_count, &
                           sun_earth_distance, geo_sub_lon, geo_sub_lat, &
                           c1, c2, a1_20, a2_20, nu_20, &
                           a1_31, a2_31, nu_31, a1_32, a2_32, nu_32, &
                           solar_Ch20_nu,timerr_seconds, mask_mode, &
                           acha_mode, acha_mode_user, dcomp_mode, wmo_sc_code, &
                           platform_name,sensor_name,dark_name,mask_name, &
                           subset_flag, lat_south_subset, lat_north_subset, &
                           lon_west_subset, lon_east_subset, subset_name, &
                           lwir_focal_plane_temperature)

    integer(kind=int4), intent(in):: num_cells,num_cells_with_data
    integer, intent(in):: therm_cal_1b,Ref_cal_1b,nav_opt,use_sst_anal, &
                       modis_clr_alb_flag, nwp_opt

    integer(kind=int4), intent(in):: start_time
    integer(kind=int4), intent(in):: end_time
    character(len=*), intent(in):: mask_mode
    character(len=*), intent(in):: acha_mode, acha_mode_user
    integer(kind=int4), intent(in):: dcomp_mode
    integer(kind=int2), intent(in)::  start_year,end_year,start_day,end_day
    real(kind=real4), intent(in)::  resolution_km, dlat, ch1_gain_low, ch1_gain_high, &
                           ch2_gain_low, ch2_gain_high, &
                           ch3a_gain_low, ch3a_gain_high, &
                           ch1_switch_count, ch1_dark_count, &
                           ch2_switch_count, ch2_dark_count, &
                           ch3a_switch_count, ch3a_dark_count, &
                           sun_earth_distance, &
                           c1, c2, a1_20, a2_20, nu_20, &
                           a1_31, a2_31, nu_31, a1_32, a2_32, nu_32, &
                           solar_Ch20_nu,timerr_seconds, geo_sub_lon, geo_sub_lat, &
                           lwir_focal_plane_temperature
    integer(kind=int4):: wmo_sc_code

    character(len=*), intent(in):: file_name
    character(len=*), intent(in):: file_1b
    character(len=*), intent(in):: data_type
    character(len=*), intent(in):: grid_format
    character(len=*), intent(in):: platform_name
    character(len=*), intent(in):: sensor_name
    character(len=*), intent(in):: dark_name
    character(len=*), intent(in):: mask_name
    integer(kind=int4), intent(in):: subset_flag
    real(kind=real4), intent(in):: lat_south_subset, lat_north_subset, lon_west_subset, lon_east_subset
    character(len=*), intent(in):: subset_name

    Clavrx_Global_Attr%Data_Type = Data_Type
    Clavrx_Global_Attr%File_Name = File_Name
    Clavrx_Global_Attr%File_1b = File_1b
    Clavrx_Global_Attr%Resolution_KM = Resolution_KM
    Clavrx_Global_Attr%Start_Year = Start_Year
    Clavrx_Global_Attr%End_Year = End_Year
    Clavrx_Global_Attr%Start_Day = Start_Day
    Clavrx_Global_Attr%End_Day = End_Day
    Clavrx_Global_Attr%Start_Time = Start_Time/3600000.0
    Clavrx_Global_Attr%End_Time = End_Time/3600000.0
    Clavrx_Global_Attr%Num_Cells = Num_Cells
    Clavrx_Global_Attr%Num_Cells_With_Data = Num_Cells_With_Data
    Clavrx_Global_Attr%Grid_Format = Grid_Format
    Clavrx_Global_Attr%dLat = dLat
    Clavrx_Global_Attr%Therm_Cal_1b = therm_cal_1b
    Clavrx_Global_Attr%Ref_Cal_1b = ref_cal_1b
    Clavrx_Global_Attr%Nav_Opt = nav_opt
    Clavrx_Global_Attr%Use_SST_Anal = use_sst_anal
    Clavrx_Global_Attr%Modis_Clr_Alb_Flag = Modis_Clr_alb_Flag
    Clavrx_Global_Attr%Nwp_Opt = Nwp_Opt
    Clavrx_Global_Attr%Ch1_Gain_Low = Ch1_Gain_Low
    Clavrx_Global_Attr%Ch1_Gain_High = Ch1_Gain_High
    Clavrx_Global_Attr%Ch1_Switch_Count = Ch1_Switch_Count
    Clavrx_Global_Attr%Ch1_Dark_Count = Ch1_Dark_Count
    Clavrx_Global_Attr%Ch2_Gain_Low = Ch2_Gain_Low
    Clavrx_Global_Attr%Ch2_Gain_High = Ch2_Gain_High
    Clavrx_Global_Attr%Ch2_Switch_Count = Ch2_Switch_Count
    Clavrx_Global_Attr%Ch2_Dark_Count = Ch2_Dark_Count
    Clavrx_Global_Attr%Ch3a_Gain_Low = Ch3a_Gain_Low
    Clavrx_Global_Attr%Ch3a_Gain_High = Ch3a_Gain_High
    Clavrx_Global_Attr%Ch3a_Switch_Count = Ch3a_Switch_Count
    Clavrx_Global_Attr%Ch3a_Dark_Count = Ch3a_Dark_Count
    Clavrx_Global_Attr%Sun_Earth_Distance = Sun_Earth_Distance
    Clavrx_Global_Attr%C1 = C1
    Clavrx_Global_Attr%C2 = C2
    Clavrx_Global_Attr%a_20 = A1_20
    Clavrx_Global_Attr%b_20 = A2_20
    Clavrx_Global_Attr%Nu_20 = Nu_20
    Clavrx_Global_Attr%a_31 = A1_31
    Clavrx_Global_Attr%b_31 = A2_31
    Clavrx_Global_Attr%Nu_31 = Nu_31
    Clavrx_Global_Attr%a_32 = A1_32
    Clavrx_Global_Attr%b_32 = A2_32
    Clavrx_Global_Attr%Nu_32 = Nu_32
    Clavrx_Global_Attr%Solar_Ch20_Nu = Solar_Ch20_Nu
    Clavrx_Global_Attr%Timerr_Seconds = Timerr_Seconds
    Clavrx_Global_Attr%Geo_Sub_Lon = Geo_Sub_Lon
    Clavrx_Global_Attr%Geo_Sub_Lat = Geo_Sub_Lat
    Clavrx_Global_Attr%Mask_Mode = Mask_Mode
    Clavrx_Global_Attr%Acha_Mode = Acha_Mode
    Clavrx_Global_Attr%Acha_Mode_User = Acha_Mode_User
    Clavrx_Global_Attr%Dcomp_Mode = Dcomp_Mode
    Clavrx_Global_Attr%WMO_Sc_Code = WMO_Sc_Code
    Clavrx_Global_Attr%Platform_Name = Platform_Name
    Clavrx_Global_Attr%Sensor_Name = Sensor_Name
    Clavrx_Global_Attr%Dark_Name = Dark_Name
    Clavrx_Global_Attr%Mask_Name = Mask_Name
    Clavrx_Global_Attr%Subset_Flag = Subset_Flag
    Clavrx_Global_Attr%Lat_South_Subset = Lat_South_Subset
    Clavrx_Global_Attr%Lat_North_Subset = Lat_North_Subset
    Clavrx_Global_Attr%Lon_West_Subset = Lon_West_Subset
    Clavrx_Global_Attr%Lon_East_Subset = Lon_East_Subset
    Clavrx_Global_Attr%Subset_Name = Subset_Name
    Clavrx_Global_Attr%LWIR_Focal_Plane_Temperature = LWIR_Focal_Plane_Temperature

 end subroutine SET_CLAVRX_GLOBAL_ATTRIBUTES

end module LEVEL2_STRUCTURES_MOD
