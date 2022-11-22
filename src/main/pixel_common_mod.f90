! $Id: pixel_common_mod.f90 4128 2021-04-19 02:03:04Z heidinger $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: pixel_common_mod.f90 (src)
!       PIXEL_COMMON_MOD (program)
!
! PURPOSE: This module houses routines defining the global set of variables and
!          arrays used in CLAVR-x
!
! DESCRIPTION:
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
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
! Public routines in this module:
!
! CREATE_PIXEL_ARRAYS - allocate memory for pixel level arrays
! DESTROY_PIXEL_ARRAYS - deallocate memory for pixel level arrays
! RESET_PIXEL_ARRAYS_TO_MISSING - set pixel arrays to missing
!
! File I/O: None
!
!  CLAVR-x uses MODIS channel numbers for all sensors
!
!   CLAVR-x  Modis    Avhrr   ABI    AHI    VIIRS   EPS-SG  Wavelength  Obs_Type
!     01       1       1       2      3      M5       3       0.659      solar
!     02       2       2       3      4      M7       6       0.865      solar
!     03       3       -       1      1      M3       1       0.470      solar
!     04       4       -       -      2      M4       2       0.555      solar
!     05       5       -       -      -      M8       8       1.240      solar
!     06       6       3a      5      5     M10       10      1.640      solar
!     07       7       -       6      6     M11       11      2.130      solar
!     08       8       -       -      -      M1       -       0.415      solar
!     09       9       -       -      -      M2       -       0.443      solar
!     10       10      -       -      -       -       -       0.490      solar
!     11       11      -       -      -       -       -       0.531      solar
!     12       12      -       -      -       -       -       0.565      solar
!     13       13      -       -      -       -       -       0.653      solar
!     14       14      -       -      -       -       -       0.681      solar
!     15       15      -       -      -      M6       4       0.750      solar
!     16       16      -       -      -       -       -       0.865      solar
!     17       17      -       -      -       -       7       0.905      solar
!     18       18      -       -      -       -       -       0.936      solar
!     19       19      -       -      -       -       -       0.940      solar
!     20       20      3b      7      7     M12       12      3.750      mixed
!     21       21      -       -      -       -       13      3.959      mixed
!     22       22      -       -      -     M13       -       3.959      mixed
!     23       23      -       -      -       -       14      4.050      therm
!     24       24      -       -      -       -       -       4.465      therm
!     25       25      -       -      -       -       -       4.515      therm
!     26       26      -       4      -      M9       9       1.375      solar
!     27       27      -       9      9       -       15      6.715      therm
!     28       28      -      10     10       -       16      7.325      therm
!     29       29      -      11     11     M14       29      8.550      therm
!     30       30      -      12     12       -       -       9.730      therm
!     31       31      4      14     14     M15       18     11.030      therm
!     32       32      5      15     15     M16       19     12.020      therm
!     33       33      -      16     16       -       20     13.335      therm
!     34       34      -       -      -       -       -      13.635      therm
!     35       35      -       -      -       -       -      13.935      therm
!     36       36      -       -      -       -       -      14.235      therm
!     37       -       -       8      8       -       -       6.200      therm
!     38       -       -      13     13       -       -      10.400      therm
!     39       -       -       -      -      I1       -       0.640      solar
!     40       -       -       -      -      I2       -       0.865      solar
!     41       -       -       -      -      I3       -       1.610      solar
!     42       -       -       -      -      I4       -       3.740      mixed
!     43       -       -       -      -      I5       -      11.450      therm
!     44       -       -       -      -     DNB       -       0.700      lunar
!     45       -       -       -      -       -       5       0.763      solar
!
!     * = a pseudo 13.3 channel only AVHRR/HIRS and VIIRS/CRIS IFF
!
!  Description of variables in "ch" structure:
!
!  Rad_Toa = Observed Top of Atmosphere Radiance
!  Bt_Toa = Observed Top of Atmosphere Brightness Temperature
!  Bt_Toa_Min_Sub = Top of Atmosphere BT Minimum within Pixel
!  Bt_Toa_Max_Sub = Top of Atmosphere BT Maximum within Pixel
!  Bt_Toa_Mean_Sub = Top of Atmosphere BT Maximum within Pixel Mean of sub-pixels
!  Bt_Toa_Std_Sub = Top of Atmosphere BT Stddev within Pixel
!  Rad_Toa_Min = Top of Atmosphere Radiance Minimum within Pixel
!  Rad_Toa_Max = Top of Atmosphere Radiance Maximum within Pixel
!  Rad_Toa_Clear = Simulated Top of Atmosphere Radiance under clear-sky
!  Rad_Atm = Simulated Radiance at Toa from atmospheric cloud-free emission
!  Rad_Atm_Dwn_Sfc = Simulated Downward Radiance at Sfc from atmospheric cloud-free emission
!  Trans_Atm = Simulated Transmission from Surface to Toa for cloud-free atmosphere
!              along viewing zenith angle path
!  Trans_Atm_Total = Simulated Transmission from Toa to Surface to Toa for cloud-free atmosphere
!                      along solar and viewing zenith angle path
!  Bt_Toa_Clear = Simulated Toa Brightness Temperature under clear-sky
!  Ref_Toa = Top of Atmosphere Reflectance
!  Ref_Toa_Min = Top of Atmosphere Reflectance Minimum within Pixel
!  Ref_Toa_Max = Top of Atmosphere Reflectance Maximum within Pixel
!  Ref_Toa_Std = Top of Atmosphere Reflectance Standard Deviation within Pixel
!  Ref_Sfc = Observed Reflectance adjusted as if measured at surface level
!  Sfc_Ref_White_Sky - surface reflectance under diffuse illumination
!  Emiss_Tropo - emissity of cloud placed at Tropopause needed to match toa  radiance
!  Emiss_Rel_11um - emissity of cloud relative to 11 microns
!  Emiss_Rel_10_4um - emissity of cloud relative to 10_4 microns
!  Emiss_Rel_11um_Clear - clear sky emissity of cloud relative to 11 microns
!  Emiss_Rel_10_4um_Clear - clear sky emissity of cloud relative to 10_4 microns
!  Bt_Toa_Clear = Simulated Toa Reflectance under clear-sky
!  DQF = Data Quaity Flag
!  Source = Data Source (Imager = 0, Sounder = 1)
!  Opd = Optical Depth
!  CSBT_Mask = Clear Sky Brighntess Temperature Mask
!  Opaque_Height = Maximum Height at which trans to space is 0.
!  Obs_Type = type of observation (solar,lunar,mixed,thermal).  Controls what is allocated
!
!  ISCCP-NG Vars
!  WMO_Idx_L1g =
!  Layer_Idx_L1g =
!  Sample_Mod_L1g =
!--------------------------------------------------------------------------------------
module PIXEL_COMMON_MOD

  use CONSTANTS_MOD, only: &
    real4 &
  , real8 &
  , int4 &
  , int2 &
  , int1 &
  , THERMAL_OBS_TYPE &
  , LUNAR_OBS_TYPE &
  , SOLAR_OBS_TYPE &
  , MIXED_OBS_TYPE &
  , Missing_Value_Int2 &
  , Missing_Value_Int4 &
  , Missing_Value_Int1 &
  , Missing_Value_Real4 &
  , sym &
  , nchan_clavrx &
  , Max_Num_Cld_Test_Bytes

  use CLASS_TIME_DATE, only: &
      date_type

  implicit none
  private
  public:: CREATE_PIXEL_ARRAYS, &
           DESTROY_PIXEL_ARRAYS, &
           RESET_PIXEL_ARRAYS_TO_MISSING

  private:: CREATE_NAV_ARRAYS, RESET_NAV_ARRAYS, DESTROY_NAV_ARRAYS
  private:: CREATE_GEO_ARRAYS, RESET_GEO_ARRAYS, DESTROY_GEO_ARRAYS
  private:: CREATE_SENSOR_ARRAYS, RESET_SENSOR_ARRAYS, DESTROY_SENSOR_ARRAYS
  private:: CREATE_AVHRR_ANCHOR_ARRAYS, RESET_AVHRR_ANCHOR_ARRAYS, DESTROY_AVHRR_ANCHOR_ARRAYS
  private:: CREATE_NWP_PIX_ARRAYS, RESET_NWP_PIX_ARRAYS, DESTROY_NWP_PIX_ARRAYS
  private:: CREATE_REF_CHANNEL_ARRAYS, RESET_REF_CHANNEL_ARRAYS, DESTROY_REF_CHANNEL_ARRAYS
  private:: CREATE_THERM_CHANNEL_ARRAYS, RESET_THERM_CHANNEL_ARRAYS, DESTROY_THERM_CHANNEL_ARRAYS
  private:: CREATE_EXTRA_CHANNEL_ARRAYS, RESET_EXTRA_CHANNEL_ARRAYS, DESTROY_EXTRA_CHANNEL_ARRAYS
  private:: CREATE_BTD_ARRAYS, RESET_BTD_ARRAYS, DESTROY_BTD_ARRAYS
  private:: CREATE_SURFACE_ARRAYS, RESET_SURFACE_ARRAYS, DESTROY_SURFACE_ARRAYS
  private:: CREATE_ACHA_ARRAYS, RESET_ACHA_ARRAYS, DESTROY_ACHA_ARRAYS
  private:: CREATE_BASE_ARRAYS, RESET_BASE_ARRAYS, DESTROY_BASE_ARRAYS
  private:: CREATE_DCOMP_ARRAYS, RESET_DCOMP_ARRAYS, DESTROY_DCOMP_ARRAYS
  private:: CREATE_NLCOMP_ARRAYS, RESET_NLCOMP_ARRAYS, DESTROY_NLCOMP_ARRAYS
  private:: CREATE_SASRAB_ARRAYS, RESET_SASRAB_ARRAYS, DESTROY_SASRAB_ARRAYS
  private:: CREATE_OLR_ARRAYS, RESET_OLR_ARRAYS, DESTROY_OLR_ARRAYS
  private:: CREATE_AEROSOL_ARRAYS, RESET_AEROSOL_ARRAYS, DESTROY_AEROSOL_ARRAYS
  private:: CREATE_CLOUD_MASK_ARRAYS, RESET_CLOUD_MASK_ARRAYS, DESTROY_CLOUD_MASK_ARRAYS
  private:: CREATE_CLOUD_TYPE_ARRAYS, RESET_CLOUD_TYPE_ARRAYS, DESTROY_CLOUD_TYPE_ARRAYS
  private:: CREATE_DIAGNOSTIC_ARRAYS, RESET_DIAGNOSTIC_ARRAYS, DESTROY_DIAGNOSTIC_ARRAYS
  private:: CREATE_SFC_PROD_ARRAYS, RESET_SFC_PROD_ARRAYS, DESTROY_SFC_PROD_ARRAYS
  private:: CREATE_CLOUD_PROD_ARRAYS, RESET_CLOUD_PROD_ARRAYS, DESTROY_CLOUD_PROD_ARRAYS
  private:: CREATE_NUCAPS_ARRAYS, RESET_NUCAPS_ARRAYS, DESTROY_NUCAPS_ARRAYS
  private:: CREATE_CALIOP_ARRAYS, RESET_CALIOP_ARRAYS, DESTROY_CALIOP_ARRAYS
  private:: CREATE_L1G_ARRAYS,RESET_L1G_ARRAYS,DESTROY_L1G_ARRAYS

  !--- arrays to keep track of files written to temporary directory, max 100 assumed
  integer, public, save:: Number_Of_Temporary_Files
  character(len=1020),dimension(100), public, save:: Temporary_File_Name

  !--- variables to control one-pixel diagnostic dump for ACHA (if both > 0,
  !--- dump is activated
  integer, public, save:: Elem_Abs_Idx_Acha_Dump
  integer, public, save:: Line_Abs_Idx_Acha_Dump

  integer, public :: Use_Land_IR_Emiss
  integer, public :: WMO_Id_ISCCPNG
  character(len=20), public :: Sensor_Name_ISCCPNG

  !---------------------------------------------------------------------------------
  ! CLAVR-x file list variables
  !---------------------------------------------------------------------------------
  type :: observations
    integer (kind=int1), dimension(:,:), allocatable:: DQF
    real, dimension(:,:), allocatable:: Bt_Toa
    real, dimension(:,:), allocatable:: Bt_Toa_Min_3x3
    real, dimension(:,:), allocatable:: Bt_Toa_Max_3x3
    real, dimension(:,:), allocatable:: Bt_Toa_Mean_3x3
    real, dimension(:,:), allocatable:: Bt_Toa_Std_3x3
    real, dimension(:,:), allocatable:: Bt_Toa_Min_Sub
    real, dimension(:,:), allocatable:: Bt_Toa_Max_Sub
    real, dimension(:,:), allocatable:: Bt_Toa_Mean_Sub
    real, dimension(:,:), allocatable:: Bt_Toa_Std_Sub
    real, dimension(:,:), allocatable:: Rad_Toa
    real, dimension(:,:), allocatable:: Rad_Toa_Min_3x3
    real, dimension(:,:), allocatable:: Rad_Toa_Max_3x3
    real, dimension(:,:), allocatable:: Rad_Toa_Min_Sub
    real, dimension(:,:), allocatable:: Rad_Toa_Max_Sub
    real, dimension(:,:), allocatable:: Rad_Toa_Clear
    real, dimension(:,:), allocatable:: Rad_Atm
    real, dimension(:,:), allocatable:: Rad_Atm_Dwn_Sfc
    real, dimension(:,:), allocatable:: Trans_Atm
    real, dimension(:,:), allocatable:: Trans_Atm_Total
    real, dimension(:,:), allocatable:: Bt_Toa_Clear
    real, dimension(:,:), allocatable:: Ref_Toa_Min_3x3
    real, dimension(:,:), allocatable:: Ref_Toa_Max_3x3
    real, dimension(:,:), allocatable:: Ref_Toa_Mean_3x3
    real, dimension(:,:), allocatable:: Ref_Toa_Std_3x3
    real, dimension(:,:), allocatable:: Ref_Toa_Min_Sub
    real, dimension(:,:), allocatable:: Ref_Toa_Max_Sub
    real, dimension(:,:), allocatable:: Ref_Toa_Std_Sub
    real, dimension(:,:), allocatable:: Ref_Toa
    real, dimension(:,:), allocatable:: Ref_Toa_Clear
    real, dimension(:,:), allocatable:: Ref_Toa_Clear_Min_3x3
    real, dimension(:,:), allocatable:: Ref_Toa_Clear_Max_3x3
    real, dimension(:,:), allocatable:: Ref_Toa_Clear_Std_3x3
    real, dimension(:,:), allocatable:: Ref_Sfc
    real, dimension(:,:), allocatable:: Ref_Lunar_Toa
    real, dimension(:,:), allocatable:: Ref_Lunar_Toa_Clear
    real, dimension(:,:), allocatable:: Ref_Lunar_Sfc
    real, dimension(:,:), allocatable:: Ref_Lunar_Min_3x3
    real, dimension(:,:), allocatable:: Ref_Lunar_Mean_3x3
    real, dimension(:,:), allocatable:: Ref_Lunar_Max_3x3
    real, dimension(:,:), allocatable:: Ref_Lunar_Std_3x3
    real, dimension(:,:), allocatable:: Sfc_Emiss
    real, dimension(:,:), allocatable:: Emiss_Tropo
    real, dimension(:,:), allocatable:: Emiss_Rel_11um
    real, dimension(:,:), allocatable:: Emiss_Rel_10_4um
    real, dimension(:,:), allocatable:: Emiss_Rel_11um_Clear
    real, dimension(:,:), allocatable:: Emiss_Rel_10_4um_Clear
    real, dimension(:,:), allocatable:: Sfc_Ref_White_Sky
    real, dimension(:,:), allocatable:: Sfc_Ref_White_Sky_Mean_3x3
    real, dimension(:,:), allocatable:: Opd
    integer (kind=int1), dimension(:,:), allocatable:: CSBT_Mask
    integer (kind=int1), dimension(:,:), allocatable:: Source
    real (kind=real4), dimension(:,:), allocatable:: Opaque_Height
    character(len=5) :: Obs_Type
    integer (kind=int4):: WMO_Id
    character(len=10) :: Satellite_Name
    character(len=10) :: Sensor_Name
    logical :: Fusion_Flag
    logical :: Sub_Pixel_On_Flag
    integer (kind=int1):: Native_Chan_Idx
    real (kind=real4):: Max_Focal_Plane_Temp
  end type observations

  type :: geometry_definition
     real (kind=real4), dimension(:,:), allocatable:: Satzen
     real (kind=real4), dimension(:,:), allocatable:: Solzen
     real (kind=real4), dimension(:,:), allocatable:: Solaz
     real (kind=real4), dimension(:,:), allocatable:: Sataz
     real (kind=real4), dimension(:,:), allocatable:: Relaz
     real (kind=real4), dimension(:,:), allocatable:: Glintzen
     real (kind=real4), dimension(:,:), allocatable:: Seczen
     real (kind=real4), dimension(:,:), allocatable:: Coszen
     real (kind=real4), dimension(:,:), allocatable:: CosSolzen
     real (kind=real4), dimension(:,:), allocatable:: Scatangle
     real (kind=real4), dimension(:,:), allocatable:: Airmass
     real (kind=real4), dimension(:,:), allocatable:: Glintzen_Lunar
     real (kind=real4), dimension(:,:), allocatable:: Scatangle_Lunar
     real (kind=real4), dimension(:,:), allocatable:: Lunzen
     real (kind=real4), dimension(:,:), allocatable:: Lunaz
     real (kind=real4), dimension(:,:), allocatable:: LunRelaz
     real (kind=real4), dimension(:,:), allocatable:: LunFrac
      logical, dimension(:,:), allocatable:: Space_Mask
     double precision:: Moon_Phase_Angle
     real (kind=real4):: Moon_Illum_Frac
     real (kind=real4):: Solzen_Min_Limit
     real (kind=real4):: Solzen_Max_Limit
     real (kind=real4):: Satzen_Min_Limit
     real (kind=real4):: Satzen_Max_Limit
  end type geometry_definition

  type :: navigation_definition
   integer (kind=int1), dimension(:), allocatable:: Ascend
   integer (kind=int1), dimension(:,:), allocatable:: Sounder_Fov
   integer (kind=int1), dimension(:,:), allocatable:: Sounder_Fov_Mask
   integer (kind=int2), dimension(:,:), allocatable:: Sounder_X
   integer (kind=int2), dimension(:,:), allocatable:: Sounder_Y
   integer (kind=int2), dimension(:,:), allocatable:: Sounder_Fov_Segment_Idx
   integer (kind=int4), dimension(:,:), allocatable:: X
   integer (kind=int4), dimension(:,:), allocatable:: Y
   real (kind=real4), dimension(:,:), allocatable:: Lat
   real (kind=real4), dimension(:,:), allocatable:: Lon
   real (kind=real4), dimension(:,:), allocatable:: Lat_Pc
   real (kind=real4), dimension(:,:), allocatable:: Lon_Pc
   real (kind=real4), dimension(:,:), allocatable:: Lat_1b
   real (kind=real4), dimension(:,:), allocatable:: Lon_1b
   integer (kind=int4):: Limit_Flag
   real(kind=real4):: Lat_South_Limit
   real(kind=real4):: Lat_North_Limit
   real(kind=real4):: Lon_West_Limit
   real(kind=real4):: Lon_East_Limit
   real(kind=real4):: Timerr_Seconds
   logical :: lon_lat_limits_set
   character(len=20):: Domain_Name
  end type navigation_definition

  type :: surface_definition
     integer(kind=int1), dimension(:,:), allocatable:: Land
     integer(kind=int1), dimension(:,:), allocatable:: Land_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Coast
     integer(kind=int1), dimension(:,:), allocatable:: Coast_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Coast_Mask_Nwp
     integer(kind=int1), dimension(:,:), allocatable:: Glint_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Glint_Mask_Lunar
     integer(kind=int1), dimension(:,:), allocatable:: Forward_Scatter_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Forward_Scatter_Mask_Lunar
     integer(kind=int1), dimension(:,:), allocatable:: Desert_Mask
     integer(kind=int1), dimension(:,:), allocatable:: City_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Volcano_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Snow_OISST
     integer(kind=int1), dimension(:,:), allocatable:: Snow_NWP
     integer(kind=int1), dimension(:,:), allocatable:: Snow_IMS
     integer(kind=int1), dimension(:,:), allocatable:: Snow_GLOB
     integer(kind=int1), dimension(:,:), allocatable:: Snow
     integer(kind=int1), dimension(:,:), allocatable:: Sfc_Type
     real (kind=real4), dimension(:,:), allocatable:: Zsfc
     real (kind=real4), dimension(:,:), allocatable:: Zsfc_Max
     real (kind=real4), dimension(:,:), allocatable:: Zsfc_Std
     real (kind=real4), dimension(:,:), allocatable:: Zsfc_Hires
  end type surface_definition

  type :: sensor_definition
    character(len=32):: Sensor_Name
    integer(kind=int4):: Spatial_Resolution_Meters
    character(len=32):: Platform_Name
    integer(kind=int4):: WMO_Id
    integer(kind=int4):: WMO_Id_Previous
    character(len=1020):: Instr_Const_File
    real(kind=real8):: Geo_Sub_Satellite_Longitude
    real(kind=real8):: Geo_Sub_Satellite_Latitude
    integer(kind=int1), dimension(Nchan_Clavrx):: Chan_On_Flag_Default
    integer(kind=int1), dimension(:,:), allocatable:: Chan_On_Flag_Per_Line
    integer(kind=int4):: Num_Chan_Sensor
    integer, dimension(:), allocatable:: CLAVRx_Chan_Map
    integer,dimension(:),  allocatable:: Chan_Stride
    logical :: Fusion_Flag
    integer(kind=int1), dimension(Nchan_Clavrx):: Chan_On_Flag_Initial_G17_Only
    integer(kind=int1), dimension(Nchan_Clavrx):: Chan_On_Flag_Initial_G17_Only_Previous
  end type sensor_definition

  type :: image_definition
    character(len=1020):: Level1b_Name
    character(len=1020):: Level1b_Full_Name
    character(len=1020):: Level1b_Path
    character(len=1020):: Level1b_Fusion_Name
    integer(kind=int4):: Orbit_Number
    integer(kind=int4):: Number_Of_Elements
    integer(kind=int4):: Number_Of_Lines
    integer(kind=int4):: Number_Of_Lines_Per_Segment
    integer(kind=int4):: Number_Of_Lines_Read_This_Segment
    integer(kind=int4):: Number_Of_Segments
    integer(kind=int4):: Segment_Number
    integer(kind=int2):: Start_Year
    integer(kind=int2):: Start_Doy
    integer(kind=int4):: Start_Time
    integer(kind=int2):: End_Year
    integer(kind=int2):: End_Doy
    integer(kind=int4):: End_Time
    type(date_type) :: time_start
    type(date_type) :: time_end
    integer(kind=int4):: X_Stride
    integer(kind=int4):: Y_Stride
    integer(kind=int4):: Chan_Average_Flag
    real(kind=real4):: Start_Time_Hours
    real(kind=real4):: End_Time_Hours
    real(kind=real4):: Mean_Time_Hours
    character(len=1020) :: Auxiliary_Cloud_Mask_File_Name
    character(len=1020) :: Auxiliary_Cloud_Product_File_Name
    character(len=1020) :: Auxiliary_Cloud_Type_File_Name
    character(len=1020) :: Auxiliary_Cloud_Height_File_Name
    character(len=1020) :: Auxiliary_Geolocation_File_Name
    integer (kind=int4), dimension(:), allocatable:: Scan_Time_Ms
    integer (kind=int4), dimension(:), allocatable:: Scan_Number
    real (kind=real4), dimension(:), allocatable:: Utc_Scan_Time_Hours
    logical:: Nc_Format_Flag !needed for HSD
    logical:: Area_Format_Flag
    logical:: Mixed_Resolution_Flag
    logical:: Static_Nav_Flag
    logical:: DB_Flag
    !--->type(date_type):: Start_Date, End_Date !causes compile error in lhp
  end type image_definition

  type :: cloud_mask_definition
     character(len=2020):: Classifiers_Names_Attr
     integer (kind=int1):: N_Classifiers
     integer (kind=int1),dimension(:,:),allocatable:: Cld_Mask
     integer (kind=int1),dimension(:,:),allocatable:: Cld_Mask_Binary
     integer (kind=int1),dimension(:,:),allocatable:: Cld_Mask_IR
     integer (kind=int1),dimension(:,:),allocatable:: Cld_Mask_Binary_IR
     integer (kind=int1),dimension(:,:),allocatable:: Cld_Mask_Aux
     integer (kind=int1),dimension(:,:),allocatable:: Adj_Pix_Cld_Mask
     integer (kind=int1),dimension(:,:),allocatable:: Cld_Mask_Qf
     integer (kind=int1),dimension(:,:),allocatable:: TUT
     integer (kind=int1),dimension(:,:),allocatable:: RUT
     real (kind=real4),dimension(:,:),allocatable:: Prior_Cld_Probability
     real (kind=real4),dimension(:,:),allocatable:: Posterior_Cld_Probability
     real (kind=real4),dimension(:,:),allocatable:: Posterior_Cld_Probability_Uncer
     real (kind=real4),dimension(:,:),allocatable:: Posterior_Cld_Probability_IR
     real (kind=real4),dimension(:,:),allocatable:: Posterior_Cld_Probability_Aux
     real (kind=real4),dimension(:,:),allocatable:: Posterior_Ice_Probability
     real (kind=real4),dimension(:,:),allocatable:: Posterior_Water_Probability
     integer(kind=int1), dimension(:,:), allocatable:: Bayes_Mask_Sfc_Type
     integer (kind=int1), dimension(:,:,:), allocatable:: Cld_Test_Vector_Packed
     integer(kind=int1), dimension(:,:), allocatable:: Shadow_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Dust_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Smoke_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Fire_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Thin_Cirr_Mask
     real(kind=real4), dimension(:,:), allocatable:: Dust_Prob
  end type cloud_mask_definition

  type :: base_definition
    real (kind=real4), dimension(:,:), allocatable:: Zc_Base
    real (kind=real4), dimension(:,:), allocatable:: Pc_Base
    real (kind=real4), dimension(:,:), allocatable:: Tc_Base
    real (kind=real4), dimension(:,:), allocatable:: Geo_Thickness
    real (kind=real4), dimension(:,:), allocatable:: Base_Alt
    real (kind=real4), dimension(:,:), allocatable:: Lower_Base_Alt
    real (kind=real4), dimension(:,:), allocatable:: Lower_Base_Pc
    integer (kind=int1), dimension(:,:), allocatable:: Base_Quality_Flag
    real (kind=real4), dimension(:,:), allocatable:: Lower_Base_Pc_Uncertainty
  end type base_definition

  type :: acha_definition
    character(len=50):: Mode
    character(len=50):: Mode_User
    real (kind=real4), dimension(:,:), allocatable:: Tc_Ap
    real (kind=real4), dimension(:,:), allocatable:: Ec_Ap
    real (kind=real4), dimension(:,:), allocatable:: Beta_Ap
    real (kind=real4), dimension(:,:), allocatable:: Ice_Prob_Ap
    real (kind=real4), dimension(:,:), allocatable:: Lower_Tc_Ap
    real (kind=real4), dimension(:,:), allocatable:: Tc_Ap_Uncer
    real (kind=real4), dimension(:,:), allocatable:: Ec_Ap_Uncer
    real (kind=real4), dimension(:,:), allocatable:: Beta_Ap_Uncer
    real (kind=real4), dimension(:,:), allocatable:: Ice_Prob_Ap_Uncer
    real (kind=real4), dimension(:,:), allocatable:: Lower_Tc_Ap_Uncer
    real (kind=real4), dimension(:,:), allocatable:: Tc
    real (kind=real4), dimension(:,:), allocatable:: Tfm
    real (kind=real4), dimension(:,:), allocatable:: Es
    real (kind=real4), dimension(:,:), allocatable:: Zc_rtm
    real (kind=real4), dimension(:,:), allocatable:: Zs
    real (kind=real4), dimension(:,:), allocatable:: Ts
    real (kind=real4), dimension(:,:), allocatable:: Ec
    real (kind=real4), dimension(:,:), allocatable:: Pc
    real (kind=real4), dimension(:,:), allocatable:: Pc_Median
    real (kind=real4), dimension(:,:), allocatable:: Zc
    real (kind=real4), dimension(:,:), allocatable:: Zc_Base
    real (kind=real4), dimension(:,:), allocatable:: Pc_Base
    real (kind=real4), dimension(:,:), allocatable:: Beta
    real (kind=real4), dimension(:,:), allocatable:: Tau
    real (kind=real4), dimension(:,:), allocatable:: Tau_Uncer
    real (kind=real4), dimension(:,:), allocatable:: Reff
    real (kind=real4), dimension(:,:), allocatable:: Ice_Probability
    real (kind=real4), dimension(:,:), allocatable:: Tc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Ec_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Beta_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Zc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Pc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Lower_Tc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Lower_Zc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Lower_Pc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Ice_Probability_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Alt
    real (kind=real4), dimension(:,:), allocatable:: Cost
    real (kind=real4), dimension(:,:), allocatable:: Goodness
    real (kind=real4), dimension(:,:), allocatable:: Conv_Test
    real (kind=real4), dimension(:,:), allocatable:: Lower_Pc
    real (kind=real4), dimension(:,:), allocatable:: Lower_Zc
    real (kind=real4), dimension(:,:), allocatable:: Lower_Tc
    real (kind=real4), dimension(:,:), allocatable:: Lower_Alt
    integer(kind=int1), dimension(:,:), allocatable:: Processing_Order
    integer(kind=int1), dimension(:,:), allocatable:: Inversion_Flag
    integer (kind=int1), dimension(:,:), allocatable:: Quality_Flag
    integer (kind=int1), dimension(:,:), allocatable:: Meta_Data
    integer (kind=int1), dimension(:,:,:), allocatable:: OE_Quality_Flags
    integer (kind=int1), dimension(:,:), allocatable:: Packed_Quality_Flags
    integer (kind=int1), dimension(:,:), allocatable:: Packed_Meta_Data_Flags
    real (kind=real4), dimension(:,:), allocatable:: Conv_Cld_Prob
    real (kind=real4), dimension(:,:), allocatable:: Supercooled_Cld_Prob
    real(kind=real4):: Success_Fraction
    real(kind=real4):: Processed_Count
    real(kind=real4):: Valid_Count
    real (kind=real4), dimension(:,:), allocatable:: Ec_375um
    real (kind=real4), dimension(:,:), allocatable:: Ec_62um
    real (kind=real4), dimension(:,:), allocatable:: Ec_67um
    real (kind=real4), dimension(:,:), allocatable:: Ec_73um
    real (kind=real4), dimension(:,:), allocatable:: Ec_85um
    real (kind=real4), dimension(:,:), allocatable:: Ec_97um
    real (kind=real4), dimension(:,:), allocatable:: Ec_104um
    real (kind=real4), dimension(:,:), allocatable:: Ec_11um
    real (kind=real4), dimension(:,:), allocatable:: Ec_12um
    real (kind=real4), dimension(:,:), allocatable:: Ec_133um
    real (kind=real4), dimension(:,:), allocatable:: Ec_136um
    real (kind=real4), dimension(:,:), allocatable:: Ec_139um
    real (kind=real4), dimension(:,:), allocatable:: Ec_142um
    integer (kind=int1), dimension(:,:), allocatable:: Cld_Type
    real (kind=real4), dimension(:,:), allocatable:: Pc_Eff
    real (kind=real4), dimension(:,:), allocatable:: Tc_Eff
    real (kind=real4), dimension(:,:), allocatable:: Zc_Eff
  end type acha_definition

  type :: ccl_definition
    integer:: Mode
    integer:: Type
    integer (kind=int1), dimension(:,:), allocatable:: Cloud_Layer
    real(kind=real4), dimension(:,:), allocatable, public :: Cloud_Fraction
    real(kind=real4), dimension(:,:), allocatable, public :: Cloud_Fraction_Uncer
    integer(kind=int1), dimension(:,:), allocatable, public :: Cloud_Fraction_Layer1
    integer(kind=int1), dimension(:,:), allocatable, public :: Cloud_Fraction_Layer2
    integer(kind=int1), dimension(:,:), allocatable, public :: Cloud_Fraction_Layer3
    integer(kind=int1), dimension(:,:), allocatable, public :: Cloud_Fraction_Layer4
    integer(kind=int1), dimension(:,:), allocatable, public :: Cloud_Fraction_Layer5
    integer(kind=int1), dimension(:,:), allocatable, public :: QF
    integer (kind=int1), dimension(:,:), allocatable:: Supercooled_Cloud_Layer
    real(kind=real4), dimension(:,:), allocatable, public :: Supercooled_Cloud_Fraction
    integer(kind=int1), dimension(:,:), allocatable, public :: Supercooled_Cloud_Fraction_Layer1
    integer(kind=int1), dimension(:,:), allocatable, public :: Supercooled_Cloud_Fraction_Layer2
    integer(kind=int1), dimension(:,:), allocatable, public :: Supercooled_Cloud_Fraction_Layer3
    integer(kind=int1), dimension(:,:), allocatable, public :: Supercooled_Cloud_Fraction_Layer4
    integer(kind=int1), dimension(:,:), allocatable, public :: Supercooled_Cloud_Fraction_Layer5
    integer (kind=int1), dimension(:,:), allocatable:: Conv_Cloud_Layer
    real(kind=real4), dimension(:,:), allocatable, public :: Conv_Cloud_Fraction
    integer(kind=int1), dimension(:,:), allocatable, public :: Conv_Cloud_Fraction_Layer1
    integer(kind=int1), dimension(:,:), allocatable, public :: Conv_Cloud_Fraction_Layer2
    integer(kind=int1), dimension(:,:), allocatable, public :: Conv_Cloud_Fraction_Layer3
    integer(kind=int1), dimension(:,:), allocatable, public :: Conv_Cloud_Fraction_Layer4
    integer(kind=int1), dimension(:,:), allocatable, public :: Conv_Cloud_Fraction_Layer5
  end type ccl_definition

  type :: asos_definition
    integer (kind=int1):: Mode
    integer (kind=int1), dimension(:,:), allocatable, public:: Code
    real(kind=real4), dimension(:,:), allocatable, public:: ECA
    real(kind=real4), dimension(:,:), allocatable, public:: Zmin
    real(kind=real4), dimension(:,:), allocatable, public:: Zmax
  end type asos_definition

  type :: nucaps_definition
    real (kind=real4), dimension(:,:), allocatable, public:: Cld_Press_Layer1
    real (kind=real4), dimension(:,:), allocatable, public:: Cld_Press_Layer2
    real (kind=real4), dimension(:,:), allocatable, public:: Cld_Fraction_Layer1
    real (kind=real4), dimension(:,:), allocatable, public:: Cld_Fraction_Layer2
    real (kind=real4), dimension(:,:), allocatable, public:: Cld_Fraction
    real (kind=real4), dimension(:,:), allocatable, public:: Cld_Press
    real (kind=real4), dimension(:,:), allocatable, public:: Cld_Temp
    real (kind=real4), dimension(:,:), allocatable, public:: Cld_Height
    real (kind=real4), dimension(:,:), allocatable, public:: Cld_Temp_Smoothed
  end type nucaps_definition

  type :: nwp_pix_definition
    integer (kind=int1):: Mode
    integer:: Nwp_Opt
    integer:: Smooth_Nwp_Flag
    integer*4, allocatable, dimension(:,:):: I_Nwp
    integer*4, allocatable, dimension(:,:):: J_Nwp
    integer*4, allocatable, dimension(:,:):: I_Nwp_x
    integer*4, allocatable, dimension(:,:):: J_Nwp_x
    real, allocatable, dimension(:,:):: Lon_Nwp_Fac
    real, allocatable, dimension(:,:):: Lat_Nwp_Fac
    real (kind=real4), dimension(:,:), allocatable:: Tsfc
    real (kind=real4), dimension(:,:), allocatable:: Tair
    real (kind=real4), dimension(:,:), allocatable:: Uth
    real (kind=real4), dimension(:,:), allocatable:: Rhsfc
    real (kind=real4), dimension(:,:), allocatable:: Rh300
    real (kind=real4), dimension(:,:), allocatable:: Pmsl
    real (kind=real4), dimension(:,:), allocatable:: Psfc ! changed to target
    real (kind=real4), dimension(:,:), allocatable:: Weasd
    real (kind=real4), dimension(:,:), allocatable:: Sea_Ice_Frac
    real (kind=real4), dimension(:,:), allocatable:: Tpw
    real (kind=real4), dimension(:,:), allocatable:: Ozone
    real (kind=real4), dimension(:,:), allocatable:: Wnd_Spd_10m
    real (kind=real4), dimension(:,:), allocatable:: Wnd_Dir_10m
    real (kind=real4), dimension(:,:), allocatable:: Wnd_Spd_Cld_Top
    real (kind=real4), dimension(:,:), allocatable:: Wnd_Dir_Cld_Top
    real (kind=real4), dimension(:,:), allocatable:: Inversion_Base
    real (kind=real4), dimension(:,:), allocatable:: Inversion_Top
    real (kind=real4), dimension(:,:), allocatable:: Inversion_Strength
    real (kind=real4), dimension(:,:), allocatable:: Ttropo
    real (kind=real4), dimension(:,:), allocatable:: Ztropo
    real (kind=real4), dimension(:,:), allocatable:: Ptropo
    real (kind=real4), dimension(:,:), allocatable:: FrzPre
    real (kind=real4), dimension(:,:), allocatable:: FrzAlt
    real (kind=real4), dimension(:,:), allocatable:: HomoFrzPre
    real (kind=real4), dimension(:,:), allocatable:: HomoFrzAlt
    real (kind=real4), dimension(:,:), allocatable:: Div_Sfc
    real (kind=real4), dimension(:,:), allocatable:: Div_200

    !-- nwp parameters computed if NWP_Mode > 0
    real (kind=real4), dimension(:,:), allocatable:: K_Index
    real (kind=real4), dimension(:,:), allocatable:: Sc_Lwp
    real (kind=real4), dimension(:,:), allocatable:: Lwp
    real (kind=real4), dimension(:,:), allocatable:: Iwp
    real (kind=real4), dimension(:,:), allocatable:: Cwp
    real (kind=real4), dimension(:,:), allocatable:: Pc
    real (kind=real4), dimension(:,:), allocatable:: Tc
    real (kind=real4), dimension(:,:), allocatable:: LCL_Height
    real (kind=real4), dimension(:,:), allocatable:: CCL_Height
    real (kind=real4), dimension(:,:), allocatable:: LFC_Height
    real (kind=real4), dimension(:,:), allocatable:: EL_Height
    real (kind=real4), dimension(:,:), allocatable:: CAPE
    real (kind=real4), dimension(:,:), allocatable:: Cfrac
    integer (kind=int1), dimension(:,:), allocatable:: Ncld_Layers
    integer (kind=int1), dimension(:,:), allocatable:: Cld_Type
    real (kind=real4), dimension(:,:), allocatable:: Tpw_Above_Cloud
  end type nwp_pix_definition

  type l1g_definition
     integer(kind=2), dimension(:,:), allocatable:: WMO_Id
     integer(kind=1), dimension(:,:), allocatable:: Layer_Idx
     integer(kind=1), dimension(:,:), allocatable:: Sample_Mode
  end type l1g_definition


  !---- declare structures using above types
  type(observations), dimension(Nchan_Clavrx), public, save, target :: Ch
  type(sensor_definition), public, save, target :: Sensor
  type(image_definition), public, save, target :: Image
  type(geometry_definition), public, save, target :: Geo
  type(navigation_definition), public, save, target :: Nav
  type(surface_definition), public, save, target :: Sfc
  type(acha_definition), public, save, target :: ACHA
  type(base_definition), public, save, target :: BASE
  type(ccl_definition), public, save, target :: CCL
  type(asos_definition), public, save, target :: ASOS
  type(cloud_mask_definition), public, save, target :: CLDMASK
  type(nucaps_definition), public, save, target :: NUCAPS
  type(nwp_pix_definition), public, save, target :: NWP_PIX
  type(l1g_definition), public, save, target :: L1g

  !---- declare other global variables
  integer,public, save:: Use_Aux_Flag
  integer,public, save:: Verbose_Level_Flag
  integer,public, save:: Cloud_Mask_Aux_Read_Flag
  integer,public, save:: Cloud_Prob_Aux_Read_Flag
  integer,public, save:: Cloud_Type_Aux_Read_Flag
  integer,public, save:: Cloud_Height_Aux_Read_Flag
  integer,public, save:: Cloud_Mask_Bayesian_Flag
  integer,public, save:: Ref_cal_1b
  integer,public, save:: Therm_cal_1b
  integer,public, save:: Nav_Opt       !0=level1b,1=clevernav,2=reposnx
  integer,public, save:: Output_Format_Flag
  integer,public, save:: Level2_File_Flag
  integer,public, save:: Use_Sst_Anal
  logical,public, save:: Use_IR_Cloud_Type_Flag
  integer,public, save:: L1b_Gzip
  integer,public, save:: L1b_Bzip2
  integer,public, save:: X_Sample_Offset
  integer,public, save:: Y_Sample_Offset

  logical,public, save:: Use_Iband
  integer,public, save:: sfc_emiss_option
  integer,public, save:: Use_Sea_IR_Emiss
  integer,public, save:: Use_ABI_Dust
  integer,public, save:: Use_Sc_Prob
  integer,public, save:: Read_Volcano_Mask
  integer,public, save:: Read_Land_Mask
  integer,public, save:: Read_Coast_Mask
  integer,public, save:: Read_Surface_Elevation
  integer,public, save:: Read_Hires_Sfc_Type
  integer,public, save:: Read_Snow_Mask
  integer,public, save:: Read_GLOBSnow_Mask
  integer,public, save:: Read_Dark_Comp
  integer,public, save:: Machine_Byte_Ordering
  integer,public, save:: LRC_Flag  !local radiative center flag
  integer,public, save:: Process_Undetected_Cloud_Flag
  integer,public, save:: DCOMP_Mode
  integer,public, save:: NLCOMP_Mode
  integer,public, save:: Aerosol_Mode
  integer,public, save:: Mask_Mode
  integer,public, save:: NWP_Mode
  integer,public, save:: Cld_Flag
  integer,public, save:: Tracer_Flag
  integer,public, save:: Skip_Output
  integer,public, save:: Blank_Flag
  integer,public, save:: Sasrab_Flag
  integer,public, save:: Modis_Clr_Alb_Flag
  integer,public, save:: Rtm_Opt
  integer,public, save:: Compress_Flag
  character(len=10),public, save:: Cloud_Mask_Mode

  !---------------------------------------------------------------------------------
  ! Flags Computed within CLAVR-x that describe the sensor data
  !---------------------------------------------------------------------------------
  integer,public, save:: AVHRR_GAC_Flag
  integer,public, save:: AVHRR_KLM_Flag
  integer,public, save:: AVHRR_AAPP_Flag
  integer,public, save:: AVHRR_1_Flag
  integer,public, save:: AVHRR_IFF_Flag
  integer,public, save:: Goes_Scan_Line_Flag
  logical,public, save:: AVHRR_Fusion_Flag
  logical,public, save:: ABI_Use_104um_Flag

  !---------------------------------------------------------------------------------
  ! Internal Flags to communicate ancillary data information
  !---------------------------------------------------------------------------------
  integer,public, save:: Failed_IMS_Snow_Mask_Flag
  integer,public, save:: Failed_GLOB_Snow_Mask_Flag
  integer,public, save:: Output_Scaled_Reflectances
  integer,public, save:: Ncdc_Level2_Flag


  !---------------------------------------------------------------------------------
  ! variables that are computed to serve as attributes in the output files
  !---------------------------------------------------------------------------------
  real(kind=real4), public, save:: Orbital_Processing_Time_Minutes
  real(kind=real4), public, save:: DCOMP_Success_Fraction
  real(kind=real4), public, save:: DCOMP_Processed_Count
  real(kind=real4), public, save:: DCOMP_Valid_Count
  real(kind=real4), public, save:: Nonconfident_Cloud_Mask_Fraction
  real(kind=real4), public, save:: Nonconfident_Cloud_Mask_Count
  real(kind=real4), public, save:: Cloud_Mask_Count

  integer(kind=int4), public, save:: Byte_Swap_1b
  integer(kind=int4), public, save:: Lun_Level1b

  !---------------------------------------------------------------------------------
  ! CLAVR-x file list variables
  !---------------------------------------------------------------------------------
  character(len=1020),public,save:: Ancil_Data_Dir
  character(len=1020),public,save:: Gfs_Data_Dir
  character(len=1020),public,save:: Ncep_Data_Dir
  character(len=1020),public,save:: Cfsr_Data_Dir
  character(len=1020),public,save:: Merra_Data_Dir
  character(len=1020),public,save:: Gdas_Data_Dir
  character(len=1020),public,save:: Erai_Data_Dir
  character(len=1020),public,save:: Oisst_Data_Dir
  character(len=1020),public,save:: Snow_Data_Dir
  character(len=1020),public,save:: GLOBSnow_Data_Dir
  character(len=1020),public,save:: Dark_Comp_Data_Dir
  character(len=1020),public,save:: QRNN_CTP_Data_Dir
  character(len=1020),public,save:: Temporary_Data_Dir
  character(len=1020),public,save:: File_Nav

  character(len=1020),public,save:: Dir_Rtm
  character(len=1020),public,save:: Dir_Level2
  character(len=1020),public,save:: Bayesian_Cloud_Mask_Name
  character(len=1020),public,save:: Dark_Composite_Name
  character(len=1020),public,save:: QRNN_CTP_Name

  !----- IFF data files
  character(len=1020),public,save:: IFF_File

  !----- Static Navigation files
  character(len=1020),public,save:: Static_Nav_File

  !----- CALIOP collocation
  character(len=1020),public,save:: Caliop_Dir
  logical, public,save:: Caliop_Flag
  logical, public,save:: Skip_L1b_File_Flag
  real (kind=real4), dimension(:,:), allocatable,save,public,target:: Caliop_Num_Cld_Layers
  real (kind=real4), dimension(:,:), allocatable,save,public,target:: Caliop_Cod
  real (kind=real4), dimension(:,:), allocatable,save,public,target:: Caliop_Cld_Height

  !----- Command Line Options
  logical, public,save:: Nucaps_Flag
  logical, public,save:: Static_Dark_Sky_Flag
  logical, public,save:: QRNN_CTP_Flag

  !------------------------------------------------------------------
  !--- variables pertaining to scanline size
  !------------------------------------------------------------------
  integer, public, save:: Line_Idx_Min_Segment
  integer, public, save:: Line_Idx_Max_Segment

  integer(kind=int4), public, save:: L1b_Rec_Length
  integer(kind=int4), public, save:: Num_Anchors
  integer , public, save :: Goes_Stride

  real(kind=real4), public, save:: dLat_hist2d

  !------- pixel array declarations
  integer (kind=int1), dimension(:), allocatable, public, target:: Bad_Scan_Flag
  integer (kind=int1), dimension(:), allocatable, public:: Ch3a_On_AVHRR

  integer (kind=int2), dimension(:), allocatable, public:: Scan_Day
  integer (kind=int2), dimension(:), allocatable, public:: Scan_Year
  real (kind=real4), dimension(:,:), allocatable, public:: Pixel_Local_Time_Hours
  real (kind=real4), dimension(:,:), allocatable, public:: Pixel_Time

  real (kind=real4), dimension(:,:), allocatable,save,public:: Solzen_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Lat_Anchor_1b
  real (kind=real4), dimension(:,:), allocatable,save,public:: Lon_Anchor_1b
  real (kind=real4), dimension(:,:), allocatable,save,public:: Satzen_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Relaz_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Solaz_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Sataz_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Glintzen_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Scatangle_Anchor

  !---------------------------------------------------------------------------
  !--- parameters from avhrr Header
  !---------------------------------------------------------------------------
  integer(kind=int4), public, save:: Num_Scans_Level2_Hdf
  integer(kind=int2), public, save:: Sc_Id_Avhrr,AVHRR_Ver_1b, &
                                     AVHRR_Data_Type, Num_Loc, Tip_Parity, Aux_Sync,  &
                                     Ramp_Auto_Cal
  character(len=6), public, save:: Sc_Id_Char
  character(len=7),public,save:: Proc_Block_Id

  !--- instrument counts
  integer (kind=int2), dimension(:,:), allocatable, public,save, target:: Ch1_Counts
  integer (kind=int2), dimension(:,:), allocatable, public,save, target:: Ch2_Counts
  integer (kind=int2), dimension(:,:), allocatable, public,save, target:: Ch6_Counts
  real (kind=real4), dimension(:,:), allocatable, public,save:: Ch20_Counts_Filtered

  !--- sounder brightness temperatures
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_375um_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_11um_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_12um_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_145um_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_147um_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_149um_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_445um_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_457um_Sounder

  !--- MJH HIRS/AVHRR aux fields
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Cld_Temp_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Cld_Press_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Cld_Height_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Cld_Emiss_Sounder

  !--- calibrated observations
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Mean_3x3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Max_3x3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Min_3x3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Std_3x3

  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Sst_Anal
  real (kind=real4), dimension(:,:), allocatable, public, save:: Sst_Anal_Err
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Sst_Anal_Uni
  real (kind=real4), dimension(:,:), allocatable, public, save:: Sst_Anal_Cice
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tsfc_Retrieved
  real (kind=real4), dimension(:,:), allocatable, public, save:: Trad_Retrieved

  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Temp_Pix_Array_1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Temp_Pix_Array_2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Temp_Pix_Array_3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Diag_Pix_Array_1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Diag_Pix_Array_2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Diag_Pix_Array_3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Missing_Pixel_Array_Real4

  !uniformity metrics
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch43_Max_Sub_3x3   !max of ch43(iband) over 3x3 Mband
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Ems_Ch20_Median_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch20_Median_5x5
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch32_Bt_Ch31_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch38_Ch32_Bt_Ch38_Max_3x3

  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Covar_Ch27_Ch31_5x5
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Covar_Ch37_Ch31_5x5
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Covar_Ch27_Ch38_5x5
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Covar_Ch37_Ch38_5x5
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Elem_Idx_Max_Bt_Ch31_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Line_Idx_Max_Bt_Ch31_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Elem_Idx_Min_Bt_Ch31_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Line_Idx_Min_Bt_Ch31_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Elem_Idx_Max_Bt_Ch38_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Line_Idx_Max_Bt_Ch38_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Elem_Idx_Min_Bt_Ch38_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Line_Idx_Min_Bt_Ch38_3x3


  real (kind=real4), dimension(:,:), allocatable, public, target:: Sst_Retrieved   !sst used in cld Mask
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ndvi_Toa
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ndsi_Toa
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ndsi_Sfc
  real (kind=real4), dimension(:,:), allocatable, public, target:: Nddi_Toa
  real (kind=real4), dimension(:,:), allocatable, public, target:: Btd_Ch31_Ch32
  real (kind=real4), dimension(:,:), allocatable, public:: Btd_Ch20_Ch31
  real (kind=real4), dimension(:,:), allocatable, public:: Btd_Ch20_Ch32
  real (kind=real4), dimension(:,:), allocatable, public, target:: Btd_Ch38_Ch32
  real (kind=real4), dimension(:,:), allocatable, public:: Btd_Ch20_Ch38

  integer(kind=int1), dimension(:,:), allocatable, public, target:: Solar_Contamination_Mask
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Bad_Pixel_Mask

  real, public:: Segment_Valid_Fraction

  !--- viirs arrays
  integer(kind=int4), dimension(:,:), allocatable, public:: Gap_Pixel_Mask_Pattern
  integer(kind=int4), dimension(:,:), allocatable, public:: Gap_Line_Idx_Pattern
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Gap_Pixel_Mask
  integer(kind=int4), dimension(:,:), allocatable, public:: Gap_Line_Idx
  integer(kind=int1), dimension(:,:), allocatable, public:: IFF_Gap_Mask

  !--- cloud phase arrays
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Type
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Phase
  real (kind=real4),dimension(:,:),allocatable, public, save, target:: Cld_Phase_Uncertainty
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Type_IR
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Phase_IR
 !integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Ctp_Multilayer_Flag

  !--- Auxilliary variables
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Metadata_Aux
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Type_Aux
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Phase_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tc_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ec_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_Top1_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_Top2_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_Uncertainty1_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_Uncertainty2_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Cost_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tau_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Reff_Aux

  !--- pixel level cloud props

     !--- direct heights
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_H2O
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tc_H2O
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_H2O
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ec_H2O
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zclr_H2O_Peak
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_CO2IRW
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_CO2IRW
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tc_CO2IRW
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ec_CO2IRW
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_SplitWin
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_SplitWin
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tc_SplitWin
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ec_SplitWin

     !-- DCOMP cloud algorithm results
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_Ap
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: vis_Ref_fm
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Iwp_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Iwp_Tau_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Lwp_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Refl_Asym_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_Ice_Layer_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_Water_Layer_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_Scwater_Layer_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Iwc_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Lwc_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Rain_Rate_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Hcld_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cdnc_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_Cost
     real (kind=real4), dimension(:,:), allocatable, public, target, save:: Reff_DCOMP_Cost
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_Qf
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP_Qf
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: DCOMP_Quality_Flag
     integer (kind=int2), dimension(:,:), allocatable, public,target, save:: DCOMP_Info_Flag
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_Fit
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP_Fit
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_1
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_2
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_3
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP_1
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP_2
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP_3

     !-- Nlcomp cloud algorithm results
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_Nlcomp
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_Nlcomp
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_Nlcomp_Cost
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_Nlcomp_Cost
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: Nlcomp_Quality_Flag
     integer (kind=int2), dimension(:,:), allocatable, public,target, save:: Nlcomp_Info_Flag

     !--- pixel level aerosol props
     real (kind=real4), dimension(:,:), allocatable, public, target, save:: Aot1
     real (kind=real4), dimension(:,:), allocatable, public, target, save:: Aot2
     real (kind=real4), dimension(:,:), allocatable, public, target, save:: Aot3a
     integer (kind=int1), dimension(:,:), allocatable, public, target, save:: Aot_Qf

     !--- non-cloud properties
     integer (kind=int1), dimension(:,:), allocatable, public, save:: Ndvi_Qf
     integer (kind=int1), dimension(:,:), allocatable, public, save:: Tsfc_Qf

     !--- pixel level radiative flux props from DCOMP
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Olr
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_VIS_Albedo
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_VIS_Spherical_Albedo
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Above_Cloud_VIS_Transmission
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_VIS_Transmission_View
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_VIS_Transmission_Solar
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_VIS_Transmission_AC
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_NIR_Spherical_Albedo
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Above_Cloud_NIR_Transmission
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_NIR_Transmission_View
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_NIR_Transmission_Solar
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_NIR_Transmission_AC
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Insolation_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Insolation_Diffuse_DCOMP

     real (kind=real4), dimension(:,:), allocatable, public, save,target:: refl_vis_fm_dcomp
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: refl_nir_fm_dcomp

     !--- SASRAB output
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_All_Sky
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_All_Sky_Diffuse
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_Clear_Sky
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_Cld_Opd
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_Aer_Opd

     !---- DCOMP params
     real (kind=real4), dimension(:,:), allocatable, public, save:: Cost_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public, save:: Error_Cov_Matrix_Cod
     real (kind=real4), dimension(:,:), allocatable, public, save:: Error_Cov_Matrix_Ref
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_1
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_2
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_3
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_4
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Wv1
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Wv2
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Toc_Rfl1
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Toc_Rfl2
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Virt_Alb1
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Virt_Alb2
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Toc_Rfl_Unc1
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Toc_Rfl_Unc2

     INTEGER, public, save::DCOMP_VIS_CHN
     INTEGER, public, save::DCOMP_IR_CHN
     character(200), public, save:: DCOMP_version_string

     !--- pixel level radiative flux props
     real (kind=real4), dimension(:,:), allocatable, public, save:: Rsr
     integer (kind=int1), dimension(:,:), allocatable, public, save:: Rsr_Qf

     !-- pixel level atmospheric reflectance components
     real (kind=real4), dimension(:,:), allocatable, public, target:: Ref_Ch1_Dark_Composite
     real (kind=real4), dimension(:,:), allocatable, public, target:: Static_Ref_065um_Dark_Composite
     real (kind=real4), dimension(:,:), allocatable, public, target:: Static_Ref_065um_Dark_Composite_Stddev
     real (kind=real4), dimension(:,:), allocatable, public, target:: Subpixel_Cloud_Fraction
     real (kind=real4), dimension(:,:), allocatable, public, target:: Ndvi_Sfc

  !--- scratch arrays
  integer(kind=int1), dimension(:,:), public,save,allocatable, target:: One_Byte_Temp
  integer(kind=int2), dimension(:,:), public,save,allocatable, target:: Two_Byte_Temp
  integer(kind=int1),dimension(:,:),allocatable, public, save, target:: Temp_Mask

  !--- nwp parameters
  integer*4, allocatable, dimension(:,:), public, save, target :: Zen_Idx_Rtm

  !--- local radiative center
  integer (kind=int1), allocatable, dimension(:,:), public, save, target :: Mask_LRC
  integer, allocatable, dimension(:,:), public, save, target :: i_LRC
  integer, allocatable, dimension(:,:), public, save, target :: j_LRC

  real (kind=real4), dimension(:,:), allocatable, public:: Trans_Atm_Ch20_Solar_Total_Rtm
  real (kind=real4), dimension(:,:), allocatable, public:: Trans_Atm_Ch20_Solar_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Emiss_11um_Tropo_Nadir_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_12um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_104um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_85um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_67um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_133um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Emiss_104um_Tropo_Nadir_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_104um_12um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_104um_11um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_104um_85um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_104um_67um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_104um_133um_Tropo_Rtm


  real (kind=real4), dimension(:,:), allocatable, public, target:: Pc_Opaque_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Zc_Opaque_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Tc_Opaque_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Tc_Co2
  real (kind=real4), dimension(:,:), allocatable, public, target:: Pc_Co2
  real (kind=real4), dimension(:,:), allocatable, public, target:: Zc_Co2
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ec_Co2
  real (kind=real4), dimension(:,:), allocatable, public, target:: Tc_Cirrus_Background
  real (kind=real4), dimension(:,:), allocatable, public, target:: Zc_Cirrus_Background
  real (kind=real4), dimension(:,:), allocatable, public, target:: Cloud_Fraction_Background

  !--- modis white sky albedo maps
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ndvi_Sfc_White_Sky

  !--- BCM Variables
  real (kind=real4), dimension(:,:), allocatable, public, target:: BTD_11_12um_NWC
  real (kind=real4), dimension(:,:), allocatable, public, target:: Emiss_39um_NWC
  real (kind=real4), dimension(:,:), allocatable, public, target:: Emiss_Tropo_11um_LRC
  real (kind=real4), dimension(:,:), allocatable, public, target:: Emiss_Tropo_104um_LRC

  !---  QRNN CTP
  real (kind=real4), dimension(:,:), allocatable, public, target:: QRNN_CTP
  real (kind=real4), dimension(:,:), allocatable, public, target:: QRNN_CTP_UNCER


  !------------------------------------------------------------------------------------
  !---- other static arrays carried by this module
  !------------------------------------------------------------------------------------


  !--- Solar RTM Terms
  type, public :: solar_rtm_struct
      real, dimension(Nchan_Clavrx,3):: Tau_H2O_Coef
      real, dimension(Nchan_Clavrx):: Tau_Ray
      real, dimension(Nchan_Clavrx):: Tau_O2
      real, dimension(Nchan_Clavrx):: Tau_O3
      real, dimension(Nchan_Clavrx):: Tau_CH4
      real, dimension(Nchan_Clavrx):: Tau_CO2
      real, dimension(Nchan_Clavrx):: Tau_Aer
      real, dimension(Nchan_Clavrx):: Wo_Aer
      real, dimension(Nchan_Clavrx):: G_Aer
  end type solar_rtm_struct

  type (solar_rtm_struct), public, save:: Solar_Rtm

  !--- flags for using clavrxorb_Default_file
  integer ,public, save :: Use_Default

  !--- clavrxorb_File_List filename
  character(len=1020), public, save:: File_List
  character(len=1020), public, save:: Level2_List

 contains

!----------------------------------------------------------------------------
! This routine allocate the memory for the pixel arrays
!----------------------------------------------------------------------------
subroutine CREATE_PIXEL_ARRAYS()
  integer:: dim1
  integer:: dim2
  integer:: idx

  !--- allocate pixel arrays
  dim1 = Image%Number_Of_Elements
  dim2 = Image%Number_Of_Lines_Per_Segment

  !---- new
!----------------------------------------------------------------------
!  begin allocation of ch structure
!----------------------------------------------------------------------

  !--- set obs type for each channel
  Ch(1:19)%Obs_Type = SOLAR_OBS_TYPE
  Ch(20)%Obs_Type = MIXED_OBS_TYPE
  Ch(21)%Obs_Type = MIXED_OBS_TYPE
  Ch(22:25)%Obs_Type = THERMAL_OBS_TYPE
  Ch(26)%Obs_Type = SOLAR_OBS_TYPE
  Ch(27:38)%Obs_Type = THERMAL_OBS_TYPE
  Ch(39:41)%Obs_Type = SOLAR_OBS_TYPE
  Ch(42)%Obs_Type = MIXED_OBS_TYPE
  Ch(43)%Obs_Type = THERMAL_OBS_TYPE
  Ch(44)%Obs_Type = LUNAR_OBS_TYPE
  Ch(45)%Obs_Type = SOLAR_OBS_TYPE

  !--- loop through each that is on, allocate fields based on obs type
  do idx = 1,Nchan_Clavrx
      if (Sensor%Chan_On_Flag_Default(idx) == sym%YES) then

        allocate(Ch(idx)%DQF(dim1,dim2))
        allocate(Ch(idx)%Source(dim1,dim2))

        select case (ch(idx)%Obs_Type)

        case(SOLAR_OBS_TYPE)
            allocate(Ch(idx)%Ref_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Sfc(dim1,dim2))
            allocate(Ch(idx)%Ref_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm_Total(dim1,dim2))
            if (idx == 1 .or. idx == 6 .or. idx == 26) allocate(Ch(idx)%Opd(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Ref_Toa_Min_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Ref_Toa_Max_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Ref_Toa_Std_Sub(dim1,dim2))

        case(LUNAR_OBS_TYPE)
            allocate(Ch(idx)%Rad_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Lunar_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Lunar_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Ref_Lunar_Sfc(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm_Total(dim1,dim2))
            allocate(Ch(idx)%Opd(dim1,dim2))

        case(THERMAL_OBS_TYPE)
            allocate(Ch(idx)%Rad_Toa(dim1,dim2))
            allocate(Ch(idx)%Bt_Toa(dim1,dim2))
            allocate(Ch(idx)%Rad_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Bt_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Rad_Atm(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm(dim1,dim2))
            allocate(Ch(idx)%Sfc_Emiss(dim1,dim2))
            if (idx == 31) allocate(Ch(idx)%Rad_Atm_Dwn_Sfc(dim1,dim2))
            if (idx == 38) allocate(Ch(idx)%Rad_Atm_Dwn_Sfc(dim1,dim2))
            !-- compute emiss tropo for these specific ir channels
            !-- make sure these channels have emiss_tropo allocated in pixel_common
            if ( any ( idx ==  [27,29,31,32,33,37,38] ) ) then
                allocate(Ch(idx)%Emiss_Tropo(dim1,dim2))
            endif
            if (idx >= 27 .and. idx <= 38)  allocate(Ch(idx)%CSBT_Mask(dim1,dim2))
            if (idx >= 27 .and. idx <= 38)  allocate(Ch(idx)%Opaque_Height(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Bt_Toa_Min_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Bt_Toa_Max_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Bt_Toa_Mean_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Bt_Toa_Std_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Rad_Toa_Min_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Rad_Toa_Max_Sub(dim1,dim2))

        case(MIXED_OBS_TYPE)
            allocate(Ch(idx)%Ref_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Sfc(dim1,dim2))
            allocate(Ch(idx)%Ref_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm_Total(dim1,dim2))
            allocate(Ch(idx)%Rad_Toa(dim1,dim2))
            allocate(Ch(idx)%Bt_Toa(dim1,dim2))
            allocate(Ch(idx)%Rad_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Bt_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Rad_Atm(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm(dim1,dim2))
            allocate(Ch(idx)%Sfc_Emiss(dim1,dim2))
            allocate(Ch(idx)%Emiss_Rel_11um(dim1,dim2))
            allocate(Ch(idx)%Emiss_Rel_10_4um(dim1,dim2))
            allocate(Ch(idx)%Emiss_Rel_11um_Clear(dim1,dim2))
            allocate(Ch(idx)%Emiss_Rel_10_4um_Clear(dim1,dim2))
            if (idx == 20) allocate(Ch(idx)%Opd(dim1,dim2))
            if (idx == 20) allocate(Ch(idx)%Emiss_Tropo(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Bt_Toa_Min_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Bt_Toa_Max_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Bt_Toa_Mean_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Bt_Toa_Std_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Rad_Toa_Min_Sub(dim1,dim2))
            if (Ch(idx)%Sub_Pixel_On_Flag) allocate(Ch(idx)%Rad_Toa_Max_Sub(dim1,dim2))

        case default

         print *, "CLAVR-x: Arrays not allocated for unknown channel = ",idx

        end select

      endif

  enddo

  !--- force allocation of channel 20 surface emissivity for use in mask
  if (.not. allocated(Ch(20)%Sfc_Emiss)) allocate(Ch(20)%Sfc_Emiss(dim1,dim2))

!----------------------------------------------------------------------
!  end allocation of ch structure
!----------------------------------------------------------------------

   if ((Sensor%Chan_On_Flag_Default(27) == sym%YES) .and.   &
       (Sensor%Chan_On_Flag_Default(31) == sym%YES)) then
           allocate(Covar_Ch27_Ch31_5x5(dim1,dim2))
   endif

   if ((Sensor%Chan_On_Flag_Default(37) == sym%YES) .and.   &
       (Sensor%Chan_On_Flag_Default(31) == sym%YES)) then
           allocate(Covar_Ch37_Ch31_5x5(dim1,dim2))
   endif

   if ((Sensor%Chan_On_Flag_Default(27) == sym%YES) .and.   &
       (Sensor%Chan_On_Flag_Default(38) == sym%YES)) then
           allocate(Covar_Ch27_Ch38_5x5(dim1,dim2))
   endif

   if ((Sensor%Chan_On_Flag_Default(37) == sym%YES) .and.   &
       (Sensor%Chan_On_Flag_Default(38) == sym%YES)) then
           allocate(Covar_Ch37_Ch38_5x5(dim1,dim2))
   endif

   allocate(Bad_Pixel_Mask(dim1,dim2))


   allocate(Solar_Contamination_Mask(dim1,dim2))

   !--- VIIRS Arrays
   allocate(IFF_Gap_Mask(dim1,dim2))
   allocate(Gap_Pixel_Mask(dim1,dim2))
   allocate(Gap_Line_Idx(dim1,dim2))

   allocate(   &
          Sst_Anal(dim1,dim2), &
          Sst_Anal_Err(dim1,dim2), &
          Sst_Anal_Cice(dim1,dim2), &
          Sst_Anal_Uni(dim1,dim2))

  if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
      (Sensor%Chan_On_Flag_Default(27) == sym%YES)) then
      allocate(Beta_11um_67um_Tropo_Rtm(dim1,dim2))
  endif
  if ((Sensor%Chan_On_Flag_Default(38) == sym%YES) .and. &
      (Sensor%Chan_On_Flag_Default(27) == sym%YES)) then
      allocate(Beta_104um_67um_Tropo_Rtm(dim1,dim2))
  endif
  if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
      (Sensor%Chan_On_Flag_Default(29) == sym%YES)) then
      allocate(Beta_11um_85um_Tropo_Rtm(dim1,dim2))
  endif
  if ((Sensor%Chan_On_Flag_Default(38) == sym%YES) .and. &
      (Sensor%Chan_On_Flag_Default(29) == sym%YES)) then
      allocate(Beta_104um_85um_Tropo_Rtm(dim1,dim2))
  endif
  if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
      (Sensor%Chan_On_Flag_Default(38) == sym%YES)) then
      allocate(Beta_11um_104um_Tropo_Rtm(dim1,dim2))
  endif
  if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
      (Sensor%Chan_On_Flag_Default(38) == sym%YES)) then
      allocate(Beta_104um_11um_Tropo_Rtm(dim1,dim2))
  endif
  if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
      (Sensor%Chan_On_Flag_Default(32) == sym%YES)) then
      allocate(Beta_11um_12um_Tropo_Rtm(dim1,dim2))
  endif
  if ((Sensor%Chan_On_Flag_Default(38) == sym%YES) .and. &
      (Sensor%Chan_On_Flag_Default(32) == sym%YES)) then
      allocate(Beta_104um_12um_Tropo_Rtm(dim1,dim2))
  endif
  if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
      (Sensor%Chan_On_Flag_Default(33) == sym%YES)) then
      allocate(Beta_11um_133um_Tropo_Rtm(dim1,dim2))
  endif
  if ((Sensor%Chan_On_Flag_Default(38) == sym%YES) .and. &
      (Sensor%Chan_On_Flag_Default(33) == sym%YES)) then
      allocate(Beta_104um_133um_Tropo_Rtm(dim1,dim2))
  endif

  allocate(Temp_Mask(dim1,dim2))

  !--- sensor arrays
  call  CREATE_SENSOR_ARRAYS(Nchan_Clavrx,dim2)
  !--- navigation arrays
  call  CREATE_NAV_ARRAYS(dim1, dim2)
  !--- geometry arrays
  call  CREATE_GEO_ARRAYS(dim1, dim2)
  !--- anchor point arrays
  call  CREATE_AVHRR_ANCHOR_ARRAYS(Num_Anchors, dim2)
  !--- nwp fields interpolated to the pixel level
  call  CREATE_NWP_PIX_ARRAYS(dim1, dim2)
  call  CREATE_SURFACE_ARRAYS(dim1, dim2)
  !--- ch1 arrays
  call  CREATE_REF_CHANNEL_ARRAYS(dim1, dim2)
  call  CREATE_THERM_CHANNEL_ARRAYS(dim1, dim2)
  call  CREATE_EXTRA_CHANNEL_ARRAYS(dim1, dim2)
  call  CREATE_BTD_ARRAYS(dim1, dim2)
  call  CREATE_ACHA_ARRAYS(dim1, dim2)
  call  CREATE_BASE_ARRAYS(dim1, dim2)
  call  CREATE_CCL_ARRAYS(dim1, dim2)
  if (ASOS%MODE > 0) call  CREATE_ASOS_ARRAYS(dim1, dim2)
  call  CREATE_DCOMP_ARRAYS(dim1, dim2)
  call  CREATE_NLCOMP_ARRAYS(dim1, dim2)
  call  CREATE_SASRAB_ARRAYS(dim1, dim2)
  call  CREATE_OLR_ARRAYS(dim1, dim2)
  call  CREATE_AEROSOL_ARRAYS(dim1, dim2)
  call  CREATE_CLOUD_MASK_ARRAYS(dim1, dim2, Max_Num_Cld_Test_Bytes)
  call  CREATE_CLOUD_TYPE_ARRAYS(dim1, dim2)
  call  CREATE_DIAGNOSTIC_ARRAYS(dim1, dim2)
  call  CREATE_SFC_PROD_ARRAYS(dim1, dim2)
  call  CREATE_CLOUD_PROD_ARRAYS(dim1, dim2)
  call  CREATE_NUCAPS_ARRAYS(dim1, dim2)
  call  CREATE_CALIOP_ARRAYS(dim1, dim2)
  call  CREATE_L1G_ARRAYS(dim1, dim2)

  !--- pixel level parameters
   allocate(Zen_Idx_Rtm(dim1,dim2))
   allocate(Mask_LRC(dim1,dim2))
   allocate(i_LRC(dim1,dim2))
   allocate(j_LRC(dim1,dim2))

  allocate(Ch3a_On_AVHRR(dim2), &
           Bad_Scan_Flag(dim2), &
           Image%Scan_Number(dim2), &
           Image%Scan_Time_Ms(dim2), &
           Scan_Day(dim2), &
           Scan_Year(dim2), &
           Image%Utc_Scan_Time_Hours(dim2), &
           Pixel_Local_Time_Hours(dim1,dim2), &
           Pixel_Time(dim1,dim2))

  ! Other BCM output arrays
  allocate (BTD_11_12um_NWC(dim1,dim2), &
            Emiss_39um_NWC(dim1,dim2), &
            Emiss_Tropo_11um_LRC(dim1,dim2), &
            Emiss_Tropo_104um_LRC(dim1,dim2))

  !--- QRNN CTP
  allocate (QRNN_CTP(dim1,dim2), &
            QRNN_CTP_UNCER(dim1,dim2))

  !--------------------------------------------------------------------------------
  ! Initialize variables that are not reset for each segment
  !--------------------------------------------------------------------------------

  !--- metrics - needed initialize counts to be zero but they accumulate through orbit
  DCOMP_Processed_Count = 0
  DCOMP_Valid_Count = 0
  DCOMP_Success_Fraction = Missing_Value_Real4
  Nonconfident_Cloud_Mask_Fraction = Missing_Value_Real4
  Nonconfident_Cloud_Mask_Count = 0
  Cloud_Mask_Count = 0

  !--- other
  DCOMP_version_string = "not DCOMP source id  available"


end subroutine CREATE_PIXEL_ARRAYS
!------------------------------------------------------------------------------
subroutine DESTROY_PIXEL_ARRAYS()

  integer:: idx

  !DESTROY CH Array
  do idx = 1,Nchan_Clavrx
      if (allocated(Ch(idx)%DQF)) deallocate(Ch(idx)%DQF)
      if (allocated(Ch(idx)%Source)) deallocate(Ch(idx)%Source)
      if (allocated(Ch(idx)%Bt_Toa)) deallocate(Ch(idx)%Bt_Toa)
      if (allocated(Ch(idx)%Bt_Toa_Min_3x3)) deallocate(Ch(idx)%Bt_Toa_Min_3x3)
      if (allocated(Ch(idx)%Bt_Toa_Max_3x3)) deallocate(Ch(idx)%Bt_Toa_Max_3x3)
      if (allocated(Ch(idx)%Bt_Toa_Mean_3x3)) deallocate(Ch(idx)%Bt_Toa_Mean_3x3)
      if (allocated(Ch(idx)%Bt_Toa_Std_3x3)) deallocate(Ch(idx)%Bt_Toa_Std_3x3)
      if (allocated(Ch(idx)%Bt_Toa_Min_Sub)) deallocate(Ch(idx)%Bt_Toa_Min_Sub)
      if (allocated(Ch(idx)%Bt_Toa_Max_Sub)) deallocate(Ch(idx)%Bt_Toa_Max_Sub)
      if (allocated(Ch(idx)%Bt_Toa_Mean_Sub)) deallocate(Ch(idx)%Bt_Toa_Mean_Sub)
      if (allocated(Ch(idx)%Bt_Toa_Std_Sub)) deallocate(Ch(idx)%Bt_Toa_Std_Sub)
      if (allocated(Ch(idx)%Rad_Toa)) deallocate(Ch(idx)%Rad_Toa)
      if (allocated(Ch(idx)%Rad_Toa_Min_3x3)) deallocate(Ch(idx)%Rad_Toa_Min_3x3)
      if (allocated(Ch(idx)%Rad_Toa_Max_3x3)) deallocate(Ch(idx)%Rad_Toa_Max_3x3)
      if (allocated(Ch(idx)%Rad_Toa_Min_Sub)) deallocate(Ch(idx)%Rad_Toa_Min_Sub)
      if (allocated(Ch(idx)%Rad_Toa_Max_Sub)) deallocate(Ch(idx)%Rad_Toa_Max_Sub)
      if (allocated(Ch(idx)%Rad_Toa_Clear)) deallocate(Ch(idx)%Rad_Toa_Clear)
      if (allocated(Ch(idx)%Rad_Atm)) deallocate(Ch(idx)%Rad_Atm)
      if (allocated(Ch(idx)%Rad_Atm_Dwn_Sfc)) deallocate(Ch(idx)%Rad_Atm_Dwn_Sfc)
      if (allocated(Ch(idx)%Trans_Atm)) deallocate(Ch(idx)%Trans_Atm)
      if (allocated(Ch(idx)%Trans_Atm_Total)) deallocate(Ch(idx)%Trans_Atm_Total)
      if (allocated(Ch(idx)%Bt_Toa_Clear)) deallocate(Ch(idx)%Bt_Toa_Clear)
      if (allocated(Ch(idx)%Ref_Toa_Min_3x3)) deallocate(Ch(idx)%Ref_Toa_Min_3x3)
      if (allocated(Ch(idx)%Ref_Toa_Max_3x3)) deallocate(Ch(idx)%Ref_Toa_Max_3x3)
      if (allocated(Ch(idx)%Ref_Toa_Mean_3x3)) deallocate(Ch(idx)%Ref_Toa_Mean_3x3)
      if (allocated(Ch(idx)%Ref_Toa_Std_3x3)) deallocate(Ch(idx)%Ref_Toa_Std_3x3)
      if (allocated(Ch(idx)%Ref_Toa_Min_Sub)) deallocate(Ch(idx)%Ref_Toa_Min_Sub)
      if (allocated(Ch(idx)%Ref_Toa_Max_Sub)) deallocate(Ch(idx)%Ref_Toa_Max_Sub)
      if (allocated(Ch(idx)%Ref_Toa_Std_Sub)) deallocate(Ch(idx)%Ref_Toa_Std_Sub)
      if (allocated(Ch(idx)%Ref_Toa_Clear)) deallocate(Ch(idx)%Ref_Toa_Clear)
      if (allocated(Ch(idx)%Ref_Toa_Clear_Min_3x3)) deallocate(Ch(idx)%Ref_Toa_Clear_Min_3x3)
      if (allocated(Ch(idx)%Ref_Toa_Clear_Max_3x3)) deallocate(Ch(idx)%Ref_Toa_Clear_Max_3x3)
      if (allocated(Ch(idx)%Ref_Toa_Clear_Std_3x3)) deallocate(Ch(idx)%Ref_Toa_Clear_Std_3x3)
      if (allocated(Ch(idx)%Ref_Toa)) deallocate(Ch(idx)%Ref_Toa)
      if (allocated(Ch(idx)%Ref_Sfc)) deallocate(Ch(idx)%Ref_Sfc)
      if (allocated(Ch(idx)%Ref_Lunar_Toa)) deallocate(Ch(idx)%Ref_Lunar_Toa)
      if (allocated(Ch(idx)%Ref_Lunar_Toa_Clear)) deallocate(Ch(idx)%Ref_Lunar_Toa_Clear)
      if (allocated(Ch(idx)%Ref_Lunar_Sfc)) deallocate(Ch(idx)%Ref_Lunar_Sfc)
      if (allocated(Ch(idx)%Ref_Lunar_Min_3x3)) deallocate(Ch(idx)%Ref_Lunar_Min_3x3)
      if (allocated(Ch(idx)%Ref_Lunar_Mean_3x3)) deallocate(Ch(idx)%Ref_Lunar_Mean_3x3)
      if (allocated(Ch(idx)%Ref_Lunar_Max_3x3)) deallocate(Ch(idx)%Ref_Lunar_Max_3x3)
      if (allocated(Ch(idx)%Ref_Lunar_Std_3x3)) deallocate(Ch(idx)%Ref_Lunar_Std_3x3)
      if (allocated(Ch(idx)%Sfc_Emiss)) deallocate(Ch(idx)%Sfc_Emiss)
      if (allocated(Ch(idx)%Emiss_Tropo)) deallocate(Ch(idx)%Emiss_Tropo)
      if (allocated(Ch(idx)%Emiss_Rel_11um)) deallocate(Ch(idx)%Emiss_Rel_11um)
      if (allocated(Ch(idx)%Emiss_Rel_10_4um)) deallocate(Ch(idx)%Emiss_Rel_10_4um)
      if (allocated(Ch(idx)%Emiss_Rel_11um_Clear)) deallocate(Ch(idx)%Emiss_Rel_11um_Clear)
      if (allocated(Ch(idx)%Emiss_Rel_10_4um_Clear)) deallocate(Ch(idx)%Emiss_Rel_10_4um_Clear)
      if (allocated(Ch(idx)%Sfc_Ref_White_Sky)) deallocate(Ch(idx)%Sfc_Ref_White_Sky)
      if (allocated(Ch(idx)%Sfc_Ref_White_Sky_Mean_3x3)) deallocate(Ch(idx)%Sfc_Ref_White_Sky_Mean_3x3)
      if (allocated(Ch(idx)%Opd)) deallocate(Ch(idx)%Opd)
      if (allocated(Ch(idx)%CSBT_Mask)) deallocate(Ch(idx)%CSBT_Mask)
      if (allocated(Ch(idx)%Opaque_Height)) deallocate(Ch(idx)%Opaque_Height)
  enddo

  if (allocated(Temp_Mask)) deallocate(Temp_Mask)

  deallocate(Ch3a_On_AVHRR, &
             Bad_Scan_Flag,  &
             Image%Scan_Number, &
             Image%Scan_Time_Ms,Scan_Day,Scan_Year, &
             Image%Utc_Scan_Time_Hours, &
             Pixel_Local_Time_Hours,Pixel_Time)

  if (allocated(Covar_Ch27_Ch31_5x5)) deallocate(Covar_Ch27_Ch31_5x5)
  if (allocated(Covar_Ch37_Ch31_5x5)) deallocate(Covar_Ch37_Ch31_5x5)
  if (allocated(Covar_Ch27_Ch38_5x5)) deallocate(Covar_Ch27_Ch38_5x5)
  if (allocated(Covar_Ch37_Ch38_5x5)) deallocate(Covar_Ch37_Ch38_5x5)

  deallocate(Bad_Pixel_Mask)

  call DESTROY_SENSOR_ARRAYS()
  call DESTROY_NAV_ARRAYS()
  call DESTROY_GEO_ARRAYS()
  call DESTROY_AVHRR_ANCHOR_ARRAYS()
  call DESTROY_NWP_PIX_ARRAYS()
  call DESTROY_SURFACE_ARRAYS()
  call DESTROY_REF_CHANNEL_ARRAYS()
  call DESTROY_THERM_CHANNEL_ARRAYS()
  call DESTROY_EXTRA_CHANNEL_ARRAYS()
  call DESTROY_BTD_ARRAYS()
  call DESTROY_ACHA_ARRAYS()
  call DESTROY_BASE_ARRAYS()
  call DESTROY_CCL_ARRAYS()
  call DESTROY_ASOS_ARRAYS()
  call DESTROY_DCOMP_ARRAYS()
  call DESTROY_NLCOMP_ARRAYS()
  call DESTROY_SASRAB_ARRAYS()
  call DESTROY_OLR_ARRAYS()
  call DESTROY_AEROSOL_ARRAYS()
  call DESTROY_CLOUD_MASK_ARRAYS()
  call DESTROY_CLOUD_TYPE_ARRAYS()
  call DESTROY_DIAGNOSTIC_ARRAYS()
  call DESTROY_SFC_PROD_ARRAYS()
  call DESTROY_CLOUD_PROD_ARRAYS()
  call DESTROY_NUCAPS_ARRAYS()
  call DESTROY_CALIOP_ARRAYS()
  call DESTROY_L1G_ARRAYS()

  deallocate(Sst_Anal)
  deallocate(Sst_Anal_Err)
  deallocate(Sst_Anal_Cice)
  deallocate(Sst_Anal_Uni)

  deallocate(Solar_Contamination_Mask)

  !--- nwp and rtm indices
  if (allocated(Zen_Idx_Rtm)) deallocate(Zen_Idx_Rtm)

  !--- local radiative center indices
  if (allocated(Mask_LRC)) deallocate(Mask_LRC)
  if (allocated(i_LRC)) deallocate(i_LRC)
  if (allocated(j_LRC)) deallocate(j_LRC)

  !--- ir cloud layer
  if (allocated(Beta_11um_12um_Tropo_Rtm)) deallocate(Beta_11um_12um_Tropo_Rtm)
  if (allocated(Beta_11um_67um_Tropo_Rtm)) deallocate(Beta_11um_67um_Tropo_Rtm)
  if (allocated(Beta_11um_85um_Tropo_Rtm)) deallocate(Beta_11um_85um_Tropo_Rtm)
  if (allocated(Beta_11um_104um_Tropo_Rtm)) deallocate(Beta_11um_104um_Tropo_Rtm)
  if (allocated(Beta_11um_133um_Tropo_Rtm)) deallocate(Beta_11um_133um_Tropo_Rtm)
  if (allocated(Beta_104um_12um_Tropo_Rtm)) deallocate(Beta_104um_12um_Tropo_Rtm)
  if (allocated(Beta_104um_67um_Tropo_Rtm)) deallocate(Beta_104um_67um_Tropo_Rtm)
  if (allocated(Beta_104um_85um_Tropo_Rtm)) deallocate(Beta_104um_85um_Tropo_Rtm)
  if (allocated(Beta_104um_11um_Tropo_Rtm)) deallocate(Beta_104um_11um_Tropo_Rtm)
  if (allocated(Beta_104um_133um_Tropo_Rtm)) deallocate(Beta_104um_133um_Tropo_Rtm)

  !--- VIIRS
  if (allocated(Gap_Pixel_Mask)) deallocate(Gap_Pixel_Mask)
  if (allocated(Gap_Line_Idx)) deallocate(Gap_Line_Idx)
  if (allocated(Gap_Pixel_Mask_Pattern)) deallocate(Gap_Pixel_Mask_Pattern)
  if (allocated(Gap_Line_Idx_Pattern)) deallocate(Gap_Line_Idx_Pattern)
  if (allocated(IFF_Gap_Mask)) deallocate(IFF_Gap_Mask)

  !--- BCM arrays
  deallocate (BTD_11_12um_NWC, &
              Emiss_39um_NWC, &
              Emiss_Tropo_11um_LRC, &
              Emiss_Tropo_104um_LRC)

  deallocate (QRNN_CTP, QRNN_CTP_UNCER)

end subroutine DESTROY_PIXEL_ARRAYS

!-----------------------------------------------------------------------------
!  reset pixel arrays to missing
!-----------------------------------------------------------------------------
subroutine RESET_PIXEL_ARRAYS_TO_MISSING()

      Bad_Scan_Flag = sym%NO        !not initialized to missing
      Bad_Pixel_Mask = sym%NO      !not initialized to missing
      Ch3a_On_AVHRR = Missing_Value_Int1

      Cloud_Mask_Aux_Read_Flag = sym%NO
      Cloud_Type_Aux_Read_Flag = sym%NO
      Cloud_Height_Aux_Read_Flag = sym%NO

      call RESET_SENSOR_ARRAYS()
      call RESET_NAV_ARRAYS()
      call RESET_GEO_ARRAYS()
      call RESET_AVHRR_ANCHOR_ARRAYS()
      call RESET_NWP_PIX_ARRAYS()
      call RESET_SURFACE_ARRAYS()
      call RESET_REF_CHANNEL_ARRAYS()
      call RESET_THERM_CHANNEL_ARRAYS()
      call RESET_EXTRA_CHANNEL_ARRAYS()
      call RESET_BTD_ARRAYS()
      call RESET_CLOUD_MASK_ARRAYS()
      call RESET_ACHA_ARRAYS()
      call RESET_BASE_ARRAYS()
      call RESET_CCL_ARRAYS()
      if (ASOS%MODE > 0) call RESET_ASOS_ARRAYS()
      call RESET_DCOMP_ARRAYS()
      call RESET_NLCOMP_ARRAYS()
      call RESET_SASRAB_ARRAYS()
      call RESET_OLR_ARRAYS()
      call RESET_AEROSOL_ARRAYS()
      call RESET_CLOUD_MASK_ARRAYS()
      call RESET_CLOUD_TYPE_ARRAYS()
      call RESET_DIAGNOSTIC_ARRAYS()
      call RESET_SFC_PROD_ARRAYS()
      call RESET_CLOUD_PROD_ARRAYS()
      call RESET_NUCAPS_ARRAYS()
      call RESET_CALIOP_ARRAYS()
      call RESET_L1G_ARRAYS()

      Sst_Anal = Missing_Value_Real4
      Sst_Anal_Err = Missing_Value_Real4
      Sst_Anal_Cice = Missing_Value_Real4
      Sst_Anal_Uni = Missing_Value_Real4

      Zen_Idx_Rtm = Missing_Value_Int1

      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(32) == sym%YES)) then
        Beta_11um_12um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(38) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(32) == sym%YES)) then
        Beta_104um_12um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(27) == sym%YES)) then
        Beta_11um_67um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(38) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(27) == sym%YES)) then
        Beta_104um_67um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(29) == sym%YES)) then
        Beta_11um_85um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(38) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(29) == sym%YES)) then
        Beta_104um_85um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(38) == sym%YES)) then
        Beta_11um_104um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(38) == sym%YES)) then
        Beta_104um_11um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(33) == sym%YES)) then
        Beta_11um_133um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(38) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(33) == sym%YES)) then
        Beta_104um_133um_Tropo_Rtm = Missing_Value_Real4
      endif

      Image%Scan_Number = Missing_Value_Int4
      Image%Scan_Time_Ms = Missing_Value_Int4
      Image%Utc_Scan_Time_Hours = Missing_Value_Real4
      Pixel_Local_Time_Hours = Missing_Value_Real4
      Pixel_Time = Missing_Value_Real4

      Mask_LRC = Missing_Value_Int1
      i_LRC = Missing_Value_Int4
      j_LRC = Missing_Value_Int4

      !--- note do not reset the patterns
      Gap_Pixel_Mask = sym%NO
      Gap_Line_Idx = Missing_Value_Int4
      IFF_Gap_Mask = sym%NO

      !--- BCM arrays
      BTD_11_12um_NWC = Missing_Value_Real4
      Emiss_39um_NWC = Missing_Value_Real4
      Emiss_Tropo_11um_LRC = Missing_Value_Real4
      Emiss_Tropo_104um_LRC = Missing_Value_Real4

      !--- QRNN
      QRNN_CTP = Missing_Value_Real4
      QRNN_CTP_UNCER = Missing_Value_Real4

end subroutine RESET_PIXEL_ARRAYS_TO_MISSING
!-------------------------------------------------------------
! The following routines provide an internal structure to
! to this module
!
! What follows are a pattern of 3 routines per data group
! create - allocates arrays
! reset - sets arrays to missing
! destroy - deallocate arrays
!-------------------------------------------------------------
subroutine CREATE_SENSOR_ARRAYS(Nchan,dim2)
   integer, intent(in):: Nchan, dim2
   allocate(Sensor%Chan_On_Flag_Per_Line(Nchan,dim2))
end subroutine CREATE_SENSOR_ARRAYS
subroutine RESET_SENSOR_ARRAYS()
   if (allocated(Sensor%Chan_On_Flag_Per_Line)) Sensor%Chan_On_Flag_Per_Line = Missing_Value_Int1
end subroutine RESET_SENSOR_ARRAYS
subroutine DESTROY_SENSOR_ARRAYS()
   if (allocated(Sensor%Chan_On_Flag_Per_Line)) deallocate(Sensor%Chan_On_Flag_Per_Line)
end subroutine DESTROY_SENSOR_ARRAYS
!-------------------------------------------------------------
!
!-------------------------------------------------------------
subroutine CREATE_NAV_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(Nav%Ascend(dim2))
   allocate(Nav%X(dim1,dim2))
   allocate(Nav%Y(dim1,dim2))
   allocate(Nav%Lat(dim1,dim2))
   allocate(Nav%Lon(dim1,dim2))
   allocate(Nav%Lat_1b(dim1,dim2))
   allocate(Nav%Lon_1b(dim1,dim2))
   allocate(Nav%Lat_Pc(dim1,dim2))
   allocate(Nav%Lon_Pc(dim1,dim2))
   if (index(Sensor%Sensor_Name,'IFF') > 0) then
     allocate(Nav%Sounder_Fov(dim1,dim2))
     allocate(Nav%Sounder_Fov_Mask(dim1,dim2))
     allocate(Nav%Sounder_Fov_Segment_Idx(dim1,dim2))
     allocate(Nav%Sounder_X(dim1,dim2))
     allocate(Nav%Sounder_Y(dim1,dim2))
   endif
end subroutine CREATE_NAV_ARRAYS
subroutine DESTROY_NAV_ARRAYS()
  if (allocated(Nav%Ascend)) deallocate(Nav%Ascend)
  if (allocated(Nav%Sounder_Fov)) deallocate(Nav%Sounder_Fov)
  if (allocated(Nav%Sounder_Fov_Mask)) deallocate(Nav%Sounder_Fov_Mask)
  if (allocated(Nav%Sounder_Fov_Segment_Idx)) deallocate(Nav%Sounder_Fov_Segment_Idx)
  if (allocated(Nav%Sounder_X)) deallocate(Nav%Sounder_X)
  if (allocated(Nav%Sounder_Y)) deallocate(Nav%Sounder_Y)
  if (allocated(Nav%X)) deallocate(Nav%X)
  if (allocated(Nav%Y)) deallocate(Nav%Y)
  if (allocated(Nav%Lat)) deallocate(Nav%Lat)
  if (allocated(Nav%Lon)) deallocate(Nav%Lon)
  if (allocated(Nav%Lat_1b)) deallocate(Nav%Lat_1b)
  if (allocated(Nav%Lon_1b)) deallocate(Nav%Lon_1b)
  if (allocated(Nav%Lat_Pc)) deallocate(Nav%Lat_Pc)
  if (allocated(Nav%Lon_Pc)) deallocate(Nav%Lon_Pc)
end subroutine
subroutine RESET_NAV_ARRAYS()
  if (allocated(Nav%Ascend)) Nav%Ascend = Missing_Value_Int1
  if (allocated(Nav%Sounder_Fov)) Nav%Sounder_Fov = Missing_Value_Int1
  if (allocated(Nav%Sounder_Fov_Mask)) Nav%Sounder_Fov_Mask = Missing_Value_Int1
  if (allocated(Nav%Sounder_Fov_Segment_Idx)) Nav%Sounder_Fov_Segment_Idx = Missing_Value_Int2
  if (allocated(Nav%Sounder_X)) Nav%Sounder_X = Missing_Value_Int2
  if (allocated(Nav%Sounder_Y)) Nav%Sounder_Y = Missing_Value_Int2
  if (allocated(Nav%X)) Nav%X = Missing_Value_Int4
  if (allocated(Nav%Y)) Nav%Y = Missing_Value_Int4
  if (allocated(Nav%Lat)) Nav%Lat = Missing_Value_Real4
  if (allocated(Nav%Lon)) Nav%Lon = Missing_Value_Real4
  if (allocated(Nav%Lat_1b)) Nav%Lat_1b = Missing_Value_Real4
  if (allocated(Nav%Lon_1b)) Nav%Lon_1b = Missing_Value_Real4
  if (allocated(Nav%Lat_Pc)) Nav%Lat_Pc = Missing_Value_Real4
  if (allocated(Nav%Lon_Pc)) Nav%Lon_Pc = Missing_Value_Real4
end subroutine RESET_NAV_ARRAYS
!------------------------------------------------------------------------------
!  routines to create, destroy and reset geo structure
!------------------------------------------------------------------------------
subroutine CREATE_GEO_ARRAYS(dim1,dim2)
    integer, intent(in):: dim1, dim2
    allocate (Geo%Satzen(dim1,dim2))
    allocate (Geo%Solzen(dim1,dim2))
    allocate (Geo%Sataz(dim1,dim2))
    allocate (Geo%Solaz(dim1,dim2))
    allocate (Geo%Relaz(dim1,dim2))
    allocate (Geo%Glintzen(dim1,dim2))
    allocate (Geo%Seczen(dim1,dim2))
    allocate (Geo%Coszen(dim1,dim2))
    allocate (Geo%Cossolzen(dim1,dim2))
    allocate (Geo%Scatangle(dim1,dim2))
    allocate (Geo%Airmass(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then
        allocate(Geo%Lunzen(dim1,dim2))
        allocate(Geo%Lunaz(dim1,dim2))
        allocate(Geo%LunRelaz(dim1,dim2))
        allocate(Geo%LunFrac(dim1,dim2))
        allocate(Geo%Scatangle_Lunar(dim1,dim2))
        allocate(Geo%Glintzen_Lunar(dim1,dim2))
    endif
    allocate (Geo%Space_Mask(dim1,dim2))
end subroutine CREATE_GEO_ARRAYS
subroutine DESTROY_GEO_ARRAYS()
  if (allocated(Geo%Satzen)) deallocate(Geo%Satzen)
  if (allocated(Geo%Solzen)) deallocate(Geo%Solzen)
  if (allocated(Geo%Sataz)) deallocate(Geo%Sataz)
  if (allocated(Geo%Solaz)) deallocate(Geo%Solaz)
  if (allocated(Geo%Relaz)) deallocate(Geo%Relaz)
  if (allocated(Geo%Glintzen)) deallocate(Geo%Glintzen)
  if (allocated(Geo%Seczen)) deallocate(Geo%Seczen)
  if (allocated(Geo%Coszen)) deallocate(Geo%Coszen)
  if (allocated(Geo%Cossolzen)) deallocate(Geo%Cossolzen)
  if (allocated(Geo%Scatangle)) deallocate(Geo%Scatangle)
  if (allocated(Geo%Airmass)) deallocate(Geo%Airmass)
  if (allocated(Geo%Lunzen))deallocate(Geo%Lunzen)
  if (allocated(Geo%Lunaz)) deallocate(Geo%Lunaz)
  if (allocated(Geo%LunRelaz)) deallocate(Geo%LunRelaz)
  if (allocated(Geo%LunFrac)) deallocate(Geo%LunFrac)
  if (allocated(Geo%Scatangle_Lunar)) deallocate(Geo%Scatangle_Lunar)
  if (allocated(Geo%Glintzen_Lunar)) deallocate(Geo%Glintzen_Lunar)
  if (allocated(Geo%Space_Mask)) deallocate(Geo%Space_Mask)
end subroutine DESTROY_GEO_ARRAYS
subroutine RESET_GEO_ARRAYS()
  if (allocated(Geo%Satzen)) Geo%Satzen = Missing_Value_Real4
  if (allocated(Geo%Solzen)) Geo%Solzen = Missing_Value_Real4
  if (allocated(Geo%Sataz)) Geo%Sataz = Missing_Value_Real4
  if (allocated(Geo%Solaz)) Geo%Solaz = Missing_Value_Real4
  if (allocated(Geo%Relaz)) Geo%Relaz = Missing_Value_Real4
  if (allocated(Geo%Glintzen)) Geo%Glintzen = Missing_Value_Real4
  if (allocated(Geo%Seczen)) Geo%Seczen = Missing_Value_Real4
  if (allocated(Geo%Coszen)) Geo%Coszen = Missing_Value_Real4
  if (allocated(Geo%Cossolzen)) Geo%Cossolzen = Missing_Value_Real4
  if (allocated(Geo%Scatangle)) Geo%Scatangle = Missing_Value_Real4
  if (allocated(Geo%Airmass)) Geo%Airmass = Missing_Value_Real4
  if (allocated(Geo%Lunzen)) Geo%Lunzen = Missing_Value_Real4
  if (allocated(Geo%Lunaz)) Geo%Lunaz = Missing_Value_Real4
  if (allocated(Geo%LunRelaz)) Geo%LunRelaz = Missing_Value_Real4
  if (allocated(Geo%LunFrac)) Geo%LunFrac = Missing_Value_Real4
  if (allocated(Geo%Scatangle_Lunar)) Geo%Scatangle_Lunar = Missing_Value_Real4
  if (allocated(Geo%Glintzen_Lunar)) Geo%Glintzen_Lunar = Missing_Value_Real4
  if (allocated(Geo%Space_Mask)) Geo%Space_Mask = .false.
  Geo%Moon_Phase_Angle = Missing_Value_Real4
  Geo%Moon_Illum_Frac = Missing_Value_Real4
end subroutine RESET_GEO_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_AVHRR_ANCHOR_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(Lat_Anchor_1b(dim1,dim2))
   allocate(Lon_Anchor_1b(dim1,dim2))
   allocate(Solzen_Anchor(dim1,dim2))
   allocate(Satzen_Anchor(dim1,dim2))
   allocate(Scatangle_Anchor(dim1,dim2))
   allocate(Glintzen_Anchor(dim1,dim2))
   allocate(Relaz_Anchor(dim1,dim2))
   allocate(Solaz_Anchor(dim1,dim2))
   allocate(Sataz_Anchor(dim1,dim2))
end subroutine CREATE_AVHRR_ANCHOR_ARRAYS
subroutine RESET_AVHRR_ANCHOR_ARRAYS()
  Lat_Anchor_1b = Missing_Value_Real4
  Lon_Anchor_1b = Missing_Value_Real4
  Solzen_Anchor = Missing_Value_Real4
  Satzen_Anchor = Missing_Value_Real4
  Scatangle_Anchor = Missing_Value_Real4
  Glintzen_Anchor = Missing_Value_Real4
  Relaz_Anchor = Missing_Value_Real4
  Solaz_Anchor = Missing_Value_Real4
  Sataz_Anchor = Missing_Value_Real4
end subroutine RESET_AVHRR_ANCHOR_ARRAYS
subroutine DESTROY_AVHRR_ANCHOR_ARRAYS
  deallocate(Lat_Anchor_1b)
  deallocate(Lon_Anchor_1b)
  deallocate(Solzen_Anchor)
  deallocate(Satzen_Anchor)
  deallocate(Scatangle_Anchor)
  deallocate(Glintzen_Anchor)
  deallocate(Relaz_Anchor)
  deallocate(Solaz_Anchor)
  deallocate(Sataz_Anchor)
end subroutine DESTROY_AVHRR_ANCHOR_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_NWP_PIX_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(NWP_PIX%Tair(dim1,dim2))
   allocate(NWP_PIX%Rhsfc(dim1,dim2))
   allocate(NWP_PIX%Rh300(dim1,dim2))
   allocate(NWP_PIX%Uth(dim1,dim2))
   allocate(NWP_PIX%LCL_Height(dim1,dim2))
   allocate(NWP_PIX%CCL_Height(dim1,dim2))
   allocate(NWP_PIX%LFC_Height(dim1,dim2))
   allocate(NWP_PIX%EL_Height(dim1,dim2))
   allocate(NWP_PIX%CAPE(dim1,dim2))
   allocate(NWP_PIX%Tsfc(dim1,dim2))
   allocate(NWP_PIX%Pmsl(dim1,dim2))
   allocate(NWP_PIX%Psfc(dim1,dim2))
   allocate(NWP_PIX%K_Index(dim1,dim2))
   allocate(NWP_PIX%Sc_Lwp(dim1,dim2))
   allocate(NWP_PIX%Lwp(dim1,dim2))
   allocate(NWP_PIX%Iwp(dim1,dim2))
   allocate(NWP_PIX%Cwp(dim1,dim2))
   allocate(NWP_PIX%Pc(dim1,dim2))
   allocate(NWP_PIX%Tc(dim1,dim2))
   allocate(NWP_PIX%Cfrac(dim1,dim2))
   allocate(NWP_PIX%Ncld_Layers(dim1,dim2))
   allocate(NWP_PIX%Cld_Type(dim1,dim2))
   allocate(NWP_PIX%Weasd(dim1,dim2))
   allocate(NWP_PIX%Sea_Ice_Frac(dim1,dim2))
   allocate(NWP_PIX%Tpw(dim1,dim2))
   allocate(NWP_PIX%Ozone(dim1,dim2))
   allocate(NWP_PIX%Ttropo(dim1,dim2))
   allocate(NWP_PIX%Ztropo(dim1,dim2))
   allocate(NWP_PIX%Ptropo(dim1,dim2))
   allocate(NWP_PIX%FrzPre(dim1,dim2))
   allocate(NWP_PIX%FrzAlt(dim1,dim2))
   allocate(NWP_PIX%HomoFrzPre(dim1,dim2))
   allocate(NWP_PIX%HomoFrzAlt(dim1,dim2))
   allocate(NWP_PIX%Div_Sfc(dim1,dim2))
   allocate(NWP_PIX%Div_200(dim1,dim2))
   allocate(NWP_PIX%Wnd_Spd_10m(dim1,dim2))
   allocate(NWP_PIX%Wnd_Dir_10m(dim1,dim2))
   allocate(NWP_PIX%Wnd_Spd_Cld_Top(dim1,dim2))
   allocate(NWP_PIX%Wnd_Dir_Cld_Top(dim1,dim2))
   allocate(NWP_PIX%Inversion_Base(dim1,dim2))
   allocate(NWP_PIX%Inversion_Top(dim1,dim2))
   allocate(NWP_PIX%Inversion_Strength(dim1,dim2))
   allocate(NWP_PIX%Tpw_Above_Cloud(dim1,dim2))
   allocate(NWP_PIX%I_Nwp(dim1,dim2))
   allocate(NWP_PIX%J_Nwp(dim1,dim2))
   allocate(NWP_PIX%I_Nwp_x(dim1,dim2))
   allocate(NWP_PIX%J_Nwp_x(dim1,dim2))
   allocate(NWP_PIX%Lon_Nwp_Fac(dim1,dim2))
   allocate(NWP_PIX%Lat_Nwp_Fac(dim1,dim2))
end subroutine CREATE_NWP_PIX_ARRAYS
subroutine RESET_NWP_PIX_ARRAYS()
   NWP_PIX%Tair = Missing_Value_Real4
   NWP_PIX%Rhsfc = Missing_Value_Real4
   NWP_PIX%Rh300 = Missing_Value_Real4
   NWP_PIX%Uth = Missing_Value_Real4
   NWP_PIX%LCL_Height = Missing_Value_Real4
   NWP_PIX%CCL_Height = Missing_Value_Real4
   NWP_PIX%LFC_Height = Missing_Value_Real4
   NWP_PIX%EL_Height = Missing_Value_Real4
   NWP_PIX%CAPE = Missing_Value_Real4
   NWP_PIX%Tsfc = Missing_Value_Real4
   NWP_PIX%Tsfc = Missing_Value_Real4
   NWP_PIX%Pmsl = Missing_Value_Real4
   NWP_PIX%Psfc = Missing_Value_Real4
   NWP_PIX%K_Index = Missing_Value_Real4
   NWP_PIX%Sc_Lwp = Missing_Value_Real4
   NWP_PIX%Lwp = Missing_Value_Real4
   NWP_PIX%Iwp = Missing_Value_Real4
   NWP_PIX%Cwp = Missing_Value_Real4
   NWP_PIX%Pc = Missing_Value_Real4
   NWP_PIX%Tc = Missing_Value_Real4
   NWP_PIX%Cfrac = Missing_Value_Real4
   NWP_PIX%Ncld_Layers = Missing_Value_Int1
   NWP_PIX%Cld_Type = Missing_Value_Int1
   NWP_PIX%Sea_Ice_Frac = Missing_Value_Real4
   NWP_PIX%Weasd = Missing_Value_Real4
   NWP_PIX%Tpw = Missing_Value_Real4
   NWP_PIX%Ozone = Missing_Value_Real4
   NWP_PIX%Ttropo = Missing_Value_Real4
   NWP_PIX%Ztropo = Missing_Value_Real4
   NWP_PIX%Ptropo = Missing_Value_Real4
   NWP_PIX%FrzPre = Missing_Value_Real4
   NWP_PIX%FrzAlt = Missing_Value_Real4
   NWP_PIX%HomoFrzPre = Missing_Value_Real4
   NWP_PIX%HomoFrzAlt = Missing_Value_Real4
   NWP_PIX%Div_Sfc = Missing_Value_Real4
   NWP_PIX%Div_200 = Missing_Value_Real4
   NWP_PIX%Wnd_Spd_10m = Missing_Value_Real4
   NWP_PIX%Wnd_Dir_10m = Missing_Value_Real4
   NWP_PIX%Wnd_Spd_Cld_Top = Missing_Value_Real4
   NWP_PIX%Wnd_Dir_Cld_Top = Missing_Value_Real4
   NWP_PIX%Inversion_Base = Missing_Value_Real4
   NWP_PIX%Inversion_Top = Missing_Value_Real4
   NWP_PIX%Inversion_Strength = Missing_Value_Real4
   NWP_PIX%Tpw_Above_Cloud = Missing_Value_Real4
   NWP_PIX%I_Nwp = Missing_Value_Int4
   NWP_PIX%J_Nwp = Missing_Value_Int4
   NWP_PIX%I_Nwp_x = Missing_Value_Int4
   NWP_PIX%J_Nwp_x = Missing_Value_Int4
   NWP_PIX%Lon_Nwp_Fac = Missing_Value_Real4
   NWP_PIX%Lat_Nwp_Fac = Missing_Value_Real4
end subroutine RESET_NWP_PIX_ARRAYS
subroutine DESTROY_NWP_PIX_ARRAYS()
   deallocate(NWP_PIX%Tair)
   deallocate(NWP_PIX%Rhsfc)
   deallocate(NWP_PIX%Rh300)
   deallocate(NWP_PIX%Uth)
   deallocate(NWP_PIX%LCL_Height)
   deallocate(NWP_PIX%CCL_Height)
   deallocate(NWP_PIX%LFC_Height)
   deallocate(NWP_PIX%EL_Height)
   deallocate(NWP_PIX%CAPE)
   deallocate(NWP_PIX%Tsfc)
   deallocate(NWP_PIX%Pmsl)
   deallocate(NWP_PIX%Psfc)
   deallocate(NWP_PIX%K_Index)
   deallocate(NWP_PIX%Sc_Lwp)
   deallocate(NWP_PIX%Lwp)
   deallocate(NWP_PIX%Iwp)
   deallocate(NWP_PIX%Cwp)
   deallocate(NWP_PIX%Pc)
   deallocate(NWP_PIX%Tc)
   deallocate(NWP_PIX%Cfrac)
   deallocate(NWP_PIX%Ncld_Layers)
   deallocate(NWP_PIX%Cld_Type)
   deallocate(NWP_PIX%Sea_Ice_Frac)
   deallocate(NWP_PIX%Weasd)
   deallocate(NWP_PIX%Tpw)
   deallocate(NWP_PIX%Ozone)
   deallocate(NWP_PIX%Ttropo)
   deallocate(NWP_PIX%Ztropo)
   deallocate(NWP_PIX%Ptropo)
   deallocate(NWP_PIX%FrzPre)
   deallocate(NWP_PIX%FrzAlt)
   deallocate(NWP_PIX%HomoFrzPre)
   deallocate(NWP_PIX%HomoFrzAlt)
   deallocate(NWP_PIX%Div_Sfc)
   deallocate(NWP_PIX%Div_200)
   deallocate(NWP_PIX%Wnd_Spd_10m)
   deallocate(NWP_PIX%Wnd_Dir_10m)
   deallocate(NWP_PIX%Wnd_Spd_Cld_Top)
   deallocate(NWP_PIX%Wnd_Dir_Cld_Top)
   deallocate(NWP_PIX%Inversion_Base)
   deallocate(NWP_PIX%Inversion_Top)
   deallocate(NWP_PIX%Inversion_Strength)
   deallocate(NWP_PIX%Tpw_Above_Cloud)
   deallocate(NWP_PIX%I_Nwp)
   deallocate(NWP_PIX%J_Nwp)
   deallocate(NWP_PIX%I_Nwp_x)
   deallocate(NWP_PIX%J_Nwp_x)
   deallocate(NWP_PIX%Lon_Nwp_Fac)
   deallocate(NWP_PIX%Lat_Nwp_Fac)
end subroutine DESTROY_NWP_PIX_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_L1G_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(L1g%WMO_Id(dim1,dim2))
   allocate(L1g%Layer_Idx(dim1,dim2))
   allocate(L1g%Sample_Mode(dim1,dim2))
end subroutine CREATE_L1G_ARRAYS

subroutine RESET_L1G_ARRAYS()
   L1g%WMO_Id = MISSING_VALUE_INT2
   L1g%Layer_Idx = MISSING_VALUE_INT1
   L1g%Sample_Mode = MISSING_VALUE_INT1
end subroutine RESET_L1G_ARRAYS

subroutine DESTROY_L1G_ARRAYS()
   deallocate(L1g%WMO_Id)
   deallocate(L1g%Layer_Idx)
   deallocate(L1g%Sample_Mode)
end subroutine DESTROY_L1G_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_REF_CHANNEL_ARRAYS(dim1,dim2)

   integer, intent(in):: dim1, dim2

   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
           allocate(Ch1_Counts(dim1,dim2))
           allocate(Ref_Ch1_Dark_Composite(dim1,dim2))
           allocate(Static_Ref_065um_Dark_Composite(dim1,dim2))
           allocate(Static_Ref_065um_Dark_Composite_Stddev(dim1,dim2))
           allocate(Subpixel_Cloud_Fraction(dim1,dim2))
   endif

   if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
      allocate(Ch2_Counts(dim1,dim2))
   endif

   if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
      allocate(Ch6_Counts(dim1,dim2))
   endif

end subroutine CREATE_REF_CHANNEL_ARRAYS

subroutine RESET_REF_CHANNEL_ARRAYS
   integer:: idx

   do idx = 1,Nchan_Clavrx
      if (allocated(Ch(idx)%Rad_Toa)) Ch(idx)%Rad_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Toa)) Ch(idx)%Bt_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Toa_Min_Sub)) Ch(idx)%Bt_Toa_Min_Sub = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Toa_Max_Sub)) Ch(idx)%Bt_Toa_Max_Sub = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Toa_Mean_Sub)) Ch(idx)%Bt_Toa_Mean_Sub = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Toa_Std_Sub)) Ch(idx)%Bt_Toa_Std_Sub = Missing_Value_Real4
      if (allocated(Ch(idx)%Rad_Toa_Clear)) Ch(idx)%Rad_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Toa_Clear)) Ch(idx)%Bt_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Rad_Atm)) Ch(idx)%Rad_Atm = Missing_Value_Real4
      if (allocated(Ch(idx)%Rad_Atm_Dwn_Sfc)) Ch(idx)%Rad_Atm_Dwn_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Trans_Atm)) Ch(idx)%Trans_Atm = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa)) Ch(idx)%Ref_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa_Min_Sub)) Ch(idx)%Ref_Toa_Min_Sub = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa_Max_Sub)) Ch(idx)%Ref_Toa_Max_Sub = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa_Std_Sub)) Ch(idx)%Ref_Toa_Std_Sub = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa_Clear)) Ch(idx)%Ref_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Sfc)) Ch(idx)%Ref_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Toa)) Ch(idx)%Ref_Lunar_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Toa_Clear)) Ch(idx)%Ref_Lunar_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Sfc)) Ch(idx)%Ref_Lunar_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Sfc)) Ch(idx)%Ref_Lunar_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Mean_3x3)) Ch(idx)%Ref_Lunar_Mean_3x3 = Missing_Value_Real4
      if (allocated(Ch(idx)%Opd)) Ch(idx)%Opd = Missing_Value_Real4
      if (allocated(Ch(idx)%Sfc_Emiss)) Ch(idx)%Sfc_Emiss = Missing_Value_Real4
      if (allocated(Ch(idx)%Sfc_Ref_White_Sky)) Ch(idx)%Sfc_Ref_White_Sky = Missing_Value_Real4
      if (allocated(Ch(idx)%Sfc_Ref_White_Sky_Mean_3x3)) Ch(idx)%Sfc_Ref_White_Sky_Mean_3x3 = Missing_Value_Real4
      if (allocated(Ch(idx)%Emiss_Tropo)) Ch(idx)%Emiss_Tropo = Missing_Value_Real4
      if (allocated(Ch(idx)%Emiss_Rel_11um)) Ch(idx)%Emiss_Rel_11um = Missing_Value_Real4
      if (allocated(Ch(idx)%Emiss_Rel_10_4um)) Ch(idx)%Emiss_Rel_10_4um = Missing_Value_Real4
      if (allocated(Ch(idx)%Emiss_Rel_11um_Clear)) Ch(idx)%Emiss_Rel_11um_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Emiss_Rel_10_4um_Clear)) Ch(idx)%Emiss_Rel_10_4um_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%DQF)) Ch(idx)%DQF = Missing_Value_Int1
      if (allocated(Ch(idx)%Source)) Ch(idx)%Source = Missing_Value_Int1
      if (allocated(Ch(idx)%CSBT_Mask)) Ch(idx)%CSBT_Mask = Missing_Value_Int1
      if (allocated(Ch(idx)%Opaque_Height)) Ch(idx)%Opaque_Height = Missing_Value_Real4
   enddo

   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
      Ch1_Counts = Missing_Value_Int2
      Ref_Ch1_Dark_Composite = Missing_Value_Real4
      Static_Ref_065um_Dark_Composite = Missing_Value_Real4
      Static_Ref_065um_Dark_Composite_Stddev = Missing_Value_Real4
      Subpixel_Cloud_Fraction = Missing_Value_Real4
   endif

   if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
      Ch2_Counts = Missing_Value_Int2
   endif

   if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
      Ch6_Counts = Missing_Value_Int2
   endif

end subroutine RESET_REF_CHANNEL_ARRAYS
subroutine DESTROY_REF_CHANNEL_ARRAYS

  if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
   if (allocated(Ch1_Counts)) deallocate (Ch1_Counts)
   if (allocated(Ref_Ch1_Dark_Composite)) deallocate (Ref_Ch1_Dark_Composite)
   if (allocated(Static_Ref_065um_Dark_Composite)) deallocate (Static_Ref_065um_Dark_Composite)
   if (allocated(Static_Ref_065um_Dark_Composite_Stddev)) deallocate (Static_Ref_065um_Dark_Composite_Stddev)
   if (allocated(Subpixel_Cloud_Fraction)) deallocate (Subpixel_Cloud_Fraction)
  endif

   if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
      deallocate(Ch2_Counts)
   endif

   if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
      deallocate(Ch6_Counts)
   endif

end subroutine DESTROY_REF_CHANNEL_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_THERM_CHANNEL_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2

   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       allocate(Bt_Ch20_Median_5x5(dim1,dim2))
       allocate(Ems_Ch20_Median_3x3(dim1,dim2))
       allocate(Ch20_Counts_Filtered(dim1,dim2))
       allocate(Trans_Atm_Ch20_Solar_Rtm(dim1,dim2))
       allocate(Trans_Atm_Ch20_Solar_Total_Rtm(dim1,dim2))
   endif

   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      allocate(Emiss_11um_Tropo_Nadir_Rtm(dim1,dim2))
   endif

   if (Sensor%Chan_On_Flag_Default(38) == sym%YES) then
      allocate(Emiss_104um_Tropo_Nadir_Rtm(dim1,dim2))
   endif
   if (.not. allocated(Elem_Idx_Max_Bt_Ch31_3x3))   allocate(Elem_Idx_Max_Bt_Ch31_3x3(dim1,dim2))
   if (.not. allocated(Line_Idx_Max_Bt_Ch31_3x3))   allocate(Line_Idx_Max_Bt_Ch31_3x3(dim1,dim2))
   if (.not. allocated(Elem_Idx_Min_Bt_Ch31_3x3))   allocate(Elem_Idx_Min_Bt_Ch31_3x3(dim1,dim2))
   if (.not. allocated(Line_Idx_Min_Bt_Ch31_3x3))   allocate(Line_Idx_Min_Bt_Ch31_3x3(dim1,dim2))
   if (.not. allocated(Elem_Idx_Max_Bt_Ch38_3x3))   allocate(Elem_Idx_Max_Bt_Ch38_3x3(dim1,dim2))
   if (.not. allocated(Line_Idx_Max_Bt_Ch38_3x3))   allocate(Line_Idx_Max_Bt_Ch38_3x3(dim1,dim2))
   if (.not. allocated(Elem_Idx_Min_Bt_Ch38_3x3))   allocate(Elem_Idx_Min_Bt_Ch38_3x3(dim1,dim2))
   if (.not. allocated(Line_Idx_Min_Bt_Ch38_3x3))   allocate(Line_Idx_Min_Bt_Ch38_3x3(dim1,dim2))

   if (Sensor%Chan_On_Flag_Default(43) == sym%YES) then
       allocate(Bt_Ch43_Max_Sub_3x3(dim1,dim2))
   endif

end subroutine CREATE_THERM_CHANNEL_ARRAYS

subroutine RESET_THERM_CHANNEL_ARRAYS()

   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       Bt_Ch20_Median_5x5 = Missing_Value_Real4
       Ems_Ch20_Median_3x3 = Missing_Value_Real4
       Ch20_Counts_Filtered = Missing_Value_Real4
       Trans_Atm_Ch20_Solar_Rtm = Missing_Value_Real4
       Trans_Atm_Ch20_Solar_Total_Rtm = Missing_Value_Real4
   endif

   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      Emiss_11um_Tropo_Nadir_Rtm = Missing_Value_Real4
   endif

   if (Sensor%Chan_On_Flag_Default(38) == sym%YES) then
      Emiss_104um_Tropo_Nadir_Rtm = Missing_Value_Real4
   endif
   if (allocated(Elem_Idx_Max_Bt_Ch31_3x3))   Elem_Idx_Max_Bt_Ch31_3x3 = Missing_Value_Real4
   if (allocated(Line_Idx_Max_Bt_Ch31_3x3))   Line_Idx_Max_Bt_Ch31_3x3 = Missing_Value_Real4
   if (allocated(Elem_Idx_Min_Bt_Ch31_3x3))   Elem_Idx_Min_Bt_Ch31_3x3 = Missing_Value_Real4
   if (allocated(Line_Idx_Min_Bt_Ch31_3x3))   Line_Idx_Min_Bt_Ch31_3x3 = Missing_Value_Real4

   if (Sensor%Chan_On_Flag_Default(43) == sym%YES) then
       Bt_Ch43_Max_Sub_3x3 = Missing_Value_Real4
   endif

end subroutine RESET_THERM_CHANNEL_ARRAYS

subroutine DESTROY_THERM_CHANNEL_ARRAYS()

   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
!      deallocate(Rad_Ch20_Ems)
       deallocate(Bt_Ch20_Median_5x5)
       deallocate(Ems_Ch20_Median_3x3)
       deallocate(Ch20_Counts_Filtered)
       deallocate(Trans_Atm_Ch20_Solar_Rtm)
       deallocate(Trans_Atm_Ch20_Solar_Total_Rtm)
   endif

   if (allocated(Elem_Idx_Max_Bt_Ch31_3X3)) deallocate(Elem_Idx_Max_Bt_Ch31_3x3)
   if (allocated(Line_Idx_Max_Bt_Ch31_3x3)) deallocate(Line_Idx_Max_Bt_Ch31_3x3)
   if (allocated(Elem_Idx_Min_Bt_Ch31_3X3)) deallocate(Elem_Idx_Min_Bt_Ch31_3x3)
   if (allocated(Line_Idx_Min_Bt_Ch31_3x3)) deallocate(Line_Idx_Min_Bt_Ch31_3x3)
   if (allocated(Emiss_11um_Tropo_Nadir_Rtm)) deallocate(Emiss_11um_Tropo_Nadir_Rtm)

   if (allocated(Elem_Idx_Max_Bt_Ch38_3X3)) deallocate(Elem_Idx_Max_Bt_Ch38_3x3)
   if (allocated(Line_Idx_Max_Bt_Ch38_3x3)) deallocate(Line_Idx_Max_Bt_Ch38_3x3)
   if (allocated(Elem_Idx_Min_Bt_Ch38_3X3)) deallocate(Elem_Idx_Min_Bt_Ch38_3x3)
   if (allocated(Line_Idx_Min_Bt_Ch38_3x3)) deallocate(Line_Idx_Min_Bt_Ch38_3x3)
   if (allocated(Emiss_104um_Tropo_Nadir_Rtm)) deallocate(Emiss_104um_Tropo_Nadir_Rtm)
   if (allocated(Bt_Ch43_Max_Sub_3x3)) deallocate(Bt_Ch43_Max_Sub_3x3)

end subroutine DESTROY_THERM_CHANNEL_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_EXTRA_CHANNEL_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (index(Sensor%Sensor_Name,'IFF') > 0 .or. Sensor%Fusion_Flag) then
           allocate(Bt_375um_Sounder(dim1,dim2))
           allocate(Bt_11um_Sounder(dim1,dim2))
           allocate(Bt_12um_Sounder(dim1,dim2))
   endif
   if (index(Sensor%Sensor_Name,'FUSION') > 0) then
           allocate(Bt_375um_Sounder(dim1,dim2))
           allocate(Bt_11um_Sounder(dim1,dim2))
           allocate(Bt_12um_Sounder(dim1,dim2))
           allocate(Bt_145um_Sounder(dim1,dim2))
           allocate(Bt_147um_Sounder(dim1,dim2))
           allocate(Bt_149um_Sounder(dim1,dim2))
           allocate(Bt_445um_Sounder(dim1,dim2))
           allocate(Bt_457um_Sounder(dim1,dim2))
   endif
   if (index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
           ! MJH HIRS/AVHRR aux fields.
           ! Does this belong in CREATE_EXTRA_CHANNEL_ARRAYS??
           allocate(Cld_Temp_Sounder(dim1,dim2))
           allocate(Cld_Press_Sounder(dim1,dim2))
           allocate(Cld_Height_Sounder(dim1,dim2))
           allocate(Cld_Emiss_Sounder(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then
           allocate(Ref_ChDNB_Lunar_Mean_3x3(dim1,dim2))
           allocate(Ref_ChDNB_Lunar_Max_3x3(dim1,dim2))
           allocate(Ref_ChDNB_Lunar_Min_3x3(dim1,dim2))
           allocate(Ref_ChDNB_Lunar_Std_3x3(dim1,dim2))
   endif
end subroutine CREATE_EXTRA_CHANNEL_ARRAYS

subroutine RESET_EXTRA_CHANNEL_ARRAYS()
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) Ref_ChDNB_Lunar_Mean_3x3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) Ref_ChDNB_Lunar_Max_3x3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) Ref_ChDNB_Lunar_Min_3x3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) Ref_ChDNB_Lunar_Std_3x3 = Missing_Value_Real4
      if (index(Sensor%Sensor_Name,'IFF') > 0 .or. Sensor%Fusion_Flag) then
          Bt_375um_Sounder = Missing_Value_Real4
          Bt_11um_Sounder = Missing_Value_Real4
          Bt_12um_Sounder = Missing_Value_Real4
      endif
      if (index(Sensor%Sensor_Name,'FUSION') > 0) then
          Bt_375um_Sounder = Missing_Value_Real4
          Bt_11um_Sounder = Missing_Value_Real4
          Bt_12um_Sounder = Missing_Value_Real4
          Bt_145um_Sounder = Missing_Value_Real4
          Bt_147um_Sounder = Missing_Value_Real4
          Bt_149um_Sounder = Missing_Value_Real4
          Bt_445um_Sounder = Missing_Value_Real4
          Bt_457um_Sounder = Missing_Value_Real4
      endif
      if (index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
          Cld_Temp_Sounder = Missing_Value_Real4 ! MJH
          Cld_Press_Sounder = Missing_Value_Real4
          Cld_Height_Sounder = Missing_Value_Real4
          Cld_Emiss_Sounder = Missing_Value_Real4
      endif
end subroutine RESET_EXTRA_CHANNEL_ARRAYS

subroutine DESTROY_EXTRA_CHANNEL_ARRAYS
  if (allocated(Ref_ChDNB_Lunar_Mean_3x3)) deallocate(Ref_ChDNB_Lunar_Mean_3x3)
  if (allocated(Ref_ChDNB_Lunar_Min_3x3)) deallocate(Ref_ChDNB_Lunar_Min_3x3)
  if (allocated(Ref_ChDNB_Lunar_Max_3x3)) deallocate(Ref_ChDNB_Lunar_Max_3x3)
  if (allocated(Ref_ChDNB_Lunar_Std_3x3)) deallocate(Ref_ChDNB_Lunar_Std_3x3)
  if (allocated(Bt_375um_Sounder)) deallocate(Bt_375um_Sounder)
  if (allocated(Bt_11um_Sounder)) deallocate(Bt_11um_Sounder)
  if (allocated(Bt_12um_Sounder)) deallocate(Bt_12um_Sounder)
  if (allocated(Bt_145um_Sounder)) deallocate(Bt_145um_Sounder)
  if (allocated(Bt_147um_Sounder)) deallocate(Bt_147um_Sounder)
  if (allocated(Bt_149um_Sounder)) deallocate(Bt_149um_Sounder)
  if (allocated(Bt_445um_Sounder)) deallocate(Bt_445um_Sounder)
  if (allocated(Bt_457um_Sounder)) deallocate(Bt_457um_Sounder)
  if (allocated(Cld_Temp_Sounder)) deallocate(Cld_Temp_Sounder) ! MJH
  if (allocated(Cld_Press_Sounder)) deallocate(Cld_Press_Sounder)
  if (allocated(Cld_Height_Sounder)) deallocate(Cld_Height_Sounder)
  if (allocated(Cld_Emiss_Sounder)) deallocate(Cld_Emiss_Sounder)
end subroutine DESTROY_EXTRA_CHANNEL_ARRAYS

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_BTD_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      allocate(Btd_Ch20_Ch31(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
      allocate(Btd_Ch20_Ch38(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      allocate(Btd_Ch20_Ch32(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      allocate(Btd_Ch31_Ch32(dim1,dim2))
      allocate(Btd_Ch31_Ch32_Bt_Ch31_Max_3x3(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(38) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      allocate(Btd_Ch38_Ch32(dim1,dim2))
      allocate(Btd_Ch38_Ch32_Bt_Ch38_Max_3x3(dim1,dim2))
   endif
end subroutine CREATE_BTD_ARRAYS
subroutine RESET_BTD_ARRAYS()
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      Btd_Ch20_Ch31 = Missing_Value_Real4
   endif
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
      Btd_Ch20_Ch38 = Missing_Value_Real4
   endif
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      Btd_Ch20_Ch32 = Missing_Value_Real4
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      Btd_Ch31_Ch32 = Missing_Value_Real4
      Btd_Ch31_Ch32_Bt_Ch31_Max_3x3 = Missing_Value_Real4
   endif
   if (Sensor%Chan_On_Flag_Default(38) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      Btd_Ch38_Ch32 = Missing_Value_Real4
      Btd_Ch38_Ch32_Bt_Ch38_Max_3x3 = Missing_Value_Real4
   endif
end subroutine RESET_BTD_ARRAYS
subroutine DESTROY_BTD_ARRAYS()
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      deallocate(Btd_Ch20_Ch31)
   endif
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
      deallocate(Btd_Ch20_Ch38)
   endif
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      deallocate(Btd_Ch20_Ch32)
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      deallocate(Btd_Ch31_Ch32)
      deallocate(Btd_Ch31_Ch32_Bt_Ch31_Max_3x3)
   endif
   if (Sensor%Chan_On_Flag_Default(38) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      deallocate(Btd_Ch38_Ch32)
      deallocate(Btd_Ch38_Ch32_Bt_Ch38_Max_3x3)
   endif
end subroutine DESTROY_BTD_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_SURFACE_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(Sfc%Land(dim1,dim2))
   allocate(Sfc%Land_Mask(dim1,dim2))
   allocate(Sfc%Coast(dim1,dim2))
   allocate(Sfc%Coast_Mask(dim1,dim2))
   allocate(Sfc%Coast_Mask_Nwp(dim1,dim2))
   allocate(Sfc%Glint_Mask(dim1,dim2))
   allocate(Sfc%Glint_Mask_Lunar(dim1,dim2))
   allocate(Sfc%Forward_Scatter_Mask(dim1,dim2))
   allocate(Sfc%Forward_Scatter_Mask_Lunar(dim1,dim2))
   allocate(Sfc%Desert_Mask(dim1,dim2))
   allocate(Sfc%City_Mask(dim1,dim2))
   allocate(Sfc%Volcano_Mask(dim1,dim2))
   allocate(Sfc%Snow_NWP(dim1,dim2))
   allocate(Sfc%Snow_OISST(dim1,dim2))
   allocate(Sfc%Snow_IMS(dim1,dim2))
   allocate(Sfc%Snow_GLOB(dim1,dim2))
   allocate(Sfc%Snow(dim1,dim2))
   allocate(Sfc%Sfc_Type(dim1,dim2))
   allocate(Sfc%Zsfc(dim1,dim2))
   allocate(Sfc%Zsfc_Max(dim1,dim2))
   allocate(Sfc%Zsfc_Std(dim1,dim2))
   allocate(Sfc%Zsfc_Hires(dim1,dim2))
end subroutine CREATE_SURFACE_ARRAYS
subroutine RESET_SURFACE_ARRAYS
   Sfc%Land = Missing_Value_Int1
   Sfc%Land_Mask = Missing_Value_Int1
   Sfc%Coast = Missing_Value_Int1
   Sfc%Coast_Mask = Missing_Value_Int1
   Sfc%Coast_Mask_Nwp = Missing_Value_Int1
   Sfc%Glint_Mask = Missing_Value_Int1
   Sfc%Glint_Mask_Lunar = Missing_Value_Int1
   Sfc%Forward_Scatter_Mask = Missing_Value_Int1
   Sfc%Forward_Scatter_Mask_Lunar = Missing_Value_Int1
   Sfc%Desert_Mask = Missing_Value_Int1
   Sfc%City_Mask = Missing_Value_Int1
   Sfc%Volcano_Mask = Missing_Value_Int1
   Sfc%Snow_NWP = Missing_Value_Int1
   Sfc%Snow_OISST = Missing_Value_Int1
   Sfc%Snow_IMS = Missing_Value_Int1
   Sfc%Snow_GLOB = Missing_Value_Int1
   Sfc%Snow = Missing_Value_Int1
   Sfc%Sfc_Type = Missing_Value_Int1
   Sfc%Zsfc = Missing_Value_Real4
   Sfc%Zsfc_Max = Missing_Value_Real4
   Sfc%Zsfc_Std = Missing_Value_Real4
   Sfc%Zsfc_Hires = Missing_Value_Real4
end subroutine RESET_SURFACE_ARRAYS
subroutine DESTROY_SURFACE_ARRAYS
   deallocate(Sfc%Land)
   deallocate(Sfc%Land_Mask)
   deallocate(Sfc%Coast)
   deallocate(Sfc%Coast_Mask)
   deallocate(Sfc%Coast_Mask_Nwp)
   deallocate(Sfc%Glint_Mask)
   deallocate(Sfc%Glint_Mask_Lunar)
   deallocate(Sfc%Forward_Scatter_Mask)
   deallocate(Sfc%Forward_Scatter_Mask_Lunar)
   deallocate(Sfc%Desert_Mask)
   deallocate(Sfc%City_Mask)
   deallocate(Sfc%Volcano_Mask)
   deallocate(Sfc%Snow_NWP)
   deallocate(Sfc%Snow_OISST)
   deallocate(Sfc%Snow_IMS)
   deallocate(Sfc%Snow_GLOB)
   deallocate(Sfc%Snow)
   deallocate(Sfc%Sfc_Type)
   deallocate(Sfc%Zsfc)
   deallocate(Sfc%Zsfc_Max)
   deallocate(Sfc%Zsfc_Std)
   deallocate(Sfc%Zsfc_Hires)
end subroutine DESTROY_SURFACE_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_ACHA_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then

    allocate(ACHA%Tc_Ap(dim1,dim2))
    allocate(ACHA%Ec_Ap(dim1,dim2))
    allocate(ACHA%Beta_Ap(dim1,dim2))
    allocate(ACHA%Ice_Prob_Ap(dim1,dim2))
    allocate(ACHA%Lower_Tc_Ap(dim1,dim2))
    allocate(ACHA%Tc_Ap_Uncer(dim1,dim2))
    allocate(ACHA%Ec_Ap_Uncer(dim1,dim2))
    allocate(ACHA%Beta_Ap_Uncer(dim1,dim2))
    allocate(ACHA%Ice_Prob_Ap_Uncer(dim1,dim2))
    allocate(ACHA%Lower_Tc_Ap_Uncer(dim1,dim2))
    allocate(ACHA%Tc(dim1,dim2))
    allocate(ACHA%Tfm(dim1,dim2))
    allocate(ACHA%Es(dim1,dim2))
    allocate(ACHA%Zc_rtm(dim1,dim2))
    allocate(ACHA%Zs(dim1,dim2))
    allocate(ACHA%Ts(dim1,dim2))
    allocate(ACHA%Ec(dim1,dim2))
    allocate(ACHA%Pc(dim1,dim2))
    allocate(ACHA%Pc_Median(dim1,dim2))
    allocate(ACHA%Zc(dim1,dim2))
    allocate(ACHA%Pc_Eff(dim1,dim2))
    allocate(ACHA%Tc_Eff(dim1,dim2))
    allocate(ACHA%Zc_Eff(dim1,dim2))
    allocate(Pc_Top1_Aux(dim1,dim2))
    allocate(Pc_Top2_Aux(dim1,dim2))
    allocate(ACHA%Zc_Base(dim1,dim2))
    allocate(ACHA%Pc_Base(dim1,dim2))
    allocate(ACHA%Beta(dim1,dim2))
    allocate(ACHA%Tau(dim1,dim2))
    allocate(ACHA%Tau_Uncer(dim1,dim2))
    allocate(ACHA%Reff(dim1,dim2))
    allocate(ACHA%Ice_Probability(dim1,dim2))
    allocate(ACHA%Tc_Uncertainty(dim1,dim2))
    allocate(ACHA%Ec_Uncertainty(dim1,dim2))
    allocate(ACHA%Beta_Uncertainty(dim1,dim2))
    allocate(ACHA%Zc_Uncertainty(dim1,dim2))
    allocate(ACHA%Pc_Uncertainty(dim1,dim2))
    allocate(ACHA%Lower_Tc_Uncertainty(dim1,dim2))
    allocate(ACHA%Lower_Pc_Uncertainty(dim1,dim2))
    allocate(ACHA%Lower_Zc_Uncertainty(dim1,dim2))
    allocate(ACHA%Ice_Probability_Uncertainty(dim1,dim2))
    allocate(Pc_Uncertainty1_Aux(dim1,dim2))
    allocate(Pc_Uncertainty2_Aux(dim1,dim2))
    allocate(ACHA%Alt(dim1,dim2))
    allocate(ACHA%Conv_Test(dim1,dim2))
    allocate(ACHA%Cost(dim1,dim2))
    allocate(ACHA%Goodness(dim1,dim2))
    allocate(Cost_Aux(dim1,dim2))
    allocate(ACHA%Lower_Pc(dim1,dim2))
    allocate(ACHA%Lower_Zc(dim1,dim2))
    allocate(ACHA%Lower_Tc(dim1,dim2))
    allocate(ACHA%Lower_Alt(dim1,dim2))
    allocate(ACHA%Processing_Order(dim1,dim2))
    allocate(ACHA%Inversion_Flag(dim1,dim2))
    allocate(ACHA%Quality_Flag(dim1,dim2))
    allocate(ACHA%Meta_Data(dim1,dim2))
    allocate(ACHA%OE_Quality_Flags(5,dim1,dim2))
    allocate(ACHA%Packed_Quality_Flags(dim1,dim2))
    allocate(ACHA%Packed_Meta_Data_Flags(dim1,dim2))
    allocate(ACHA%Conv_Cld_Prob(dim1,dim2))
    allocate(ACHA%Supercooled_Cld_Prob(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(20) == sym%YES) allocate(ACHA%Ec_375um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(37) == sym%YES) allocate(ACHA%Ec_62um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(27) == sym%YES) allocate(ACHA%Ec_67um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(28) == sym%YES) allocate(ACHA%Ec_73um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(29) == sym%YES) allocate(ACHA%Ec_85um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(30) == sym%YES) allocate(ACHA%Ec_97um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(38) == sym%YES) allocate(ACHA%Ec_104um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(31) == sym%YES) allocate(ACHA%Ec_11um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(32) == sym%YES) allocate(ACHA%Ec_12um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(33) == sym%YES) allocate(ACHA%Ec_133um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(34) == sym%YES) allocate(ACHA%Ec_136um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(35) == sym%YES) allocate(ACHA%Ec_139um(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(36) == sym%YES) allocate(ACHA%Ec_142um(dim1,dim2))
    allocate(ACHA%Cld_Type(dim1,dim2))
   endif

   !--- these accumulate through the whole image, do not reset with each segment
   ACHA%Processed_Count = 0
   ACHA%Valid_Count = 0
   ACHA%Success_Fraction = Missing_Value_Real4

end subroutine CREATE_ACHA_ARRAYS

subroutine RESET_ACHA_ARRAYS()

   if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
    ACHA%Tc_Ap = Missing_Value_Real4
    ACHA%Ec_Ap = Missing_Value_Real4
    ACHA%Beta_Ap = Missing_Value_Real4
    ACHA%Ice_Prob_Ap = Missing_Value_Real4
    ACHA%Lower_Tc_Ap = Missing_Value_Real4
    ACHA%Tc_Ap_Uncer = Missing_Value_Real4
    ACHA%Ec_Ap_Uncer = Missing_Value_Real4
    ACHA%Beta_Ap_Uncer = Missing_Value_Real4
    ACHA%Ice_Prob_Ap_Uncer = Missing_Value_Real4
    ACHA%Lower_Tc_Ap_Uncer = Missing_Value_Real4
    ACHA%Tc = Missing_Value_Real4
    ACHA%Tfm = Missing_Value_Real4
    ACHA%Es = Missing_Value_Real4
    ACHA%Zc_rtm = Missing_Value_Real4
    ACHA%Zs = Missing_Value_Real4
    ACHA%Ts = Missing_Value_Real4
    ACHA%Ec = Missing_Value_Real4
    ACHA%Pc = Missing_Value_Real4
    ACHA%Pc_Median = Missing_Value_Real4
    ACHA%Zc = Missing_Value_Real4
    ACHA%Pc_Eff = Missing_Value_Real4
    ACHA%Tc_Eff = Missing_Value_Real4
    ACHA%Zc_Eff = Missing_Value_Real4
    Pc_Top1_Aux = Missing_Value_Real4
    Pc_Top2_Aux = Missing_Value_Real4
    ACHA%Zc_Base  = Missing_Value_Real4
    ACHA%Pc_Base  = Missing_Value_Real4
    ACHA%Beta = Missing_Value_Real4
    ACHA%Tau = Missing_Value_Real4
    ACHA%Tau_Uncer = Missing_Value_Real4
    ACHA%Reff = Missing_Value_Real4
    ACHA%Ice_Probability = Missing_Value_Real4
    ACHA%Tc_Uncertainty = Missing_Value_Real4
    ACHA%Ec_Uncertainty = Missing_Value_Real4
    ACHA%Beta_Uncertainty = Missing_Value_Real4
    ACHA%Zc_Uncertainty = Missing_Value_Real4
    ACHA%Pc_Uncertainty = Missing_Value_Real4
    ACHA%Lower_Tc_Uncertainty = Missing_Value_Real4
    ACHA%Lower_Pc_Uncertainty = Missing_Value_Real4
    ACHA%Lower_Zc_Uncertainty = Missing_Value_Real4
    ACHA%Ice_Probability_Uncertainty = Missing_Value_Real4
    Pc_Uncertainty1_Aux = Missing_Value_Real4
    Pc_Uncertainty2_Aux = Missing_Value_Real4
    ACHA%Alt = Missing_Value_Real4
    ACHA%Cost = Missing_Value_Real4
    ACHA%Goodness = Missing_Value_Real4
    ACHA%Conv_Test = Missing_Value_Real4
    Cost_Aux = Missing_Value_Real4
    ACHA%Lower_Pc = Missing_Value_Real4
    ACHA%Lower_Zc = Missing_Value_Real4
    ACHA%Lower_Tc = Missing_Value_Real4
    ACHA%Lower_Alt = Missing_Value_Real4
    ACHA%Processing_Order = Missing_Value_Int1
    ACHA%Inversion_Flag = Missing_Value_Int1
    ACHA%Quality_Flag = Missing_Value_Int1
    ACHA%Meta_Data = 0
    ACHA%OE_Quality_Flags = 0
    ACHA%Packed_Quality_Flags = 0
    ACHA%Packed_Meta_Data_Flags = 0
    ACHA%Conv_Cld_Prob = Missing_Value_Real4
    ACHA%Supercooled_Cld_Prob = Missing_Value_Real4
    if (allocated(ACHA%Ec_375um)) ACHA%Ec_375um = Missing_Value_Real4
    if (allocated(ACHA%Ec_62um)) ACHA%Ec_62um = Missing_Value_Real4
    if (allocated(ACHA%Ec_67um)) ACHA%Ec_67um = Missing_Value_Real4
    if (allocated(ACHA%Ec_73um)) ACHA%Ec_73um = Missing_Value_Real4
    if (allocated(ACHA%Ec_85um)) ACHA%Ec_85um = Missing_Value_Real4
    if (allocated(ACHA%Ec_97um)) ACHA%Ec_97um = Missing_Value_Real4
    if (allocated(ACHA%Ec_104um)) ACHA%Ec_104um = Missing_Value_Real4
    if (allocated(ACHA%Ec_11um)) ACHA%Ec_11um = Missing_Value_Real4
    if (allocated(ACHA%Ec_12um)) ACHA%Ec_12um = Missing_Value_Real4
    if (allocated(ACHA%Ec_133um)) ACHA%Ec_133um = Missing_Value_Real4
    if (allocated(ACHA%Ec_136um)) ACHA%Ec_136um = Missing_Value_Real4
    if (allocated(ACHA%Ec_139um)) ACHA%Ec_139um = Missing_Value_Real4
    if (allocated(ACHA%Ec_142um)) ACHA%Ec_142um = Missing_Value_Real4
    ACHA%Cld_Type = Missing_Value_Int1
   endif

end subroutine RESET_ACHA_ARRAYS

subroutine DESTROY_ACHA_ARRAYS()

   if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
    deallocate(ACHA%Tc_Ap)
    deallocate(ACHA%Ec_Ap)
    deallocate(ACHA%Beta_Ap)
    deallocate(ACHA%Ice_Prob_Ap)
    deallocate(ACHA%Lower_Tc_Ap)
    deallocate(ACHA%Tc_Ap_Uncer)
    deallocate(ACHA%Ec_Ap_Uncer)
    deallocate(ACHA%Beta_Ap_Uncer)
    deallocate(ACHA%Ice_Prob_Ap_Uncer)
    deallocate(ACHA%Lower_Tc_Ap_Uncer)
    deallocate(ACHA%Tc)
    deallocate(ACHA%Tfm)
    deallocate(ACHA%Es)
    deallocate(ACHA%Zc_rtm)
    deallocate(ACHA%Zs)
    deallocate(ACHA%Ts)
    deallocate(ACHA%Ec)
    deallocate(ACHA%Pc)
    deallocate(ACHA%Pc_Median)
    deallocate(ACHA%Zc)
    deallocate(ACHA%Pc_Eff)
    deallocate(ACHA%Tc_Eff)
    deallocate(ACHA%Zc_Eff)
    deallocate(Pc_Top1_Aux)
    deallocate(Pc_Top2_Aux)
    deallocate(ACHA%Zc_Base)
    deallocate(ACHA%Pc_Base)
    deallocate(ACHA%Beta)
    deallocate(ACHA%Tau)
    deallocate(ACHA%Tau_Uncer)
    deallocate(ACHA%Reff)
    deallocate(ACHA%Ice_Probability)
    deallocate(ACHA%Tc_Uncertainty)
    deallocate(ACHA%Ec_Uncertainty)
    deallocate(ACHA%Beta_Uncertainty)
    deallocate(ACHA%Zc_Uncertainty)
    deallocate(ACHA%Pc_Uncertainty)
    deallocate(ACHA%Lower_Tc_Uncertainty)
    deallocate(ACHA%Lower_Zc_Uncertainty)
    deallocate(ACHA%Lower_Pc_Uncertainty)
    deallocate(ACHA%Ice_Probability_Uncertainty)
    deallocate(Pc_Uncertainty1_Aux)
    deallocate(Pc_Uncertainty2_Aux)
    deallocate(ACHA%Alt)
    deallocate(ACHA%Cost)
    deallocate(ACHA%Goodness)
    deallocate(ACHA%Conv_Test)
    deallocate(Cost_Aux)
    deallocate(ACHA%Lower_Pc)
    deallocate(ACHA%Lower_Zc)
    deallocate(ACHA%Lower_Tc)
    deallocate(ACHA%Lower_Alt)
    deallocate(ACHA%Processing_Order)
    deallocate(ACHA%Inversion_Flag)
    deallocate(ACHA%Quality_Flag)
    deallocate(ACHA%Meta_Data)
    deallocate(ACHA%OE_Quality_Flags)
    deallocate(ACHA%Packed_Quality_Flags)
    deallocate(ACHA%Packed_Meta_Data_Flags)
    deallocate(ACHA%Conv_Cld_Prob)
    deallocate(ACHA%Supercooled_Cld_Prob)
    if (allocated(ACHA%Ec_375um)) deallocate(ACHA%Ec_375um)
    if (allocated(ACHA%Ec_62um)) deallocate(ACHA%Ec_62um)
    if (allocated(ACHA%Ec_67um)) deallocate(ACHA%Ec_67um)
    if (allocated(ACHA%Ec_73um)) deallocate(ACHA%Ec_73um)
    if (allocated(ACHA%Ec_85um)) deallocate(ACHA%Ec_85um)
    if (allocated(ACHA%Ec_97um)) deallocate(ACHA%Ec_97um)
    if (allocated(ACHA%Ec_104um)) deallocate(ACHA%Ec_104um)
    if (allocated(ACHA%Ec_11um)) deallocate(ACHA%Ec_11um)
    if (allocated(ACHA%Ec_12um)) deallocate(ACHA%Ec_12um)
    if (allocated(ACHA%Ec_133um)) deallocate(ACHA%Ec_133um)
    if (allocated(ACHA%Ec_136um)) deallocate(ACHA%Ec_136um)
    if (allocated(ACHA%Ec_139um)) deallocate(ACHA%Ec_139um)
    if (allocated(ACHA%Ec_142um)) deallocate(ACHA%Ec_142um)
    if (allocated(ACHA%Cld_Type)) deallocate(ACHA%Cld_Type)
   endif
end subroutine DESTROY_ACHA_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_BASE_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
     allocate (BASE%Zc_Base(dim1,dim2))
     allocate (BASE%Pc_Base(dim1,dim2))
     allocate (BASE%Tc_Base(dim1,dim2))
     allocate (BASE%Geo_Thickness(dim1,dim2))
     allocate (BASE%Base_Alt(dim1,dim2))
     allocate (BASE%Lower_Base_Alt(dim1,dim2))
     allocate (BASE%Lower_Base_Pc(dim1,dim2))
     allocate (BASE%Base_Quality_Flag(dim1,dim2))
     allocate (BASE%Lower_Base_Pc_Uncertainty(dim1,dim2))
   endif
end subroutine CREATE_BASE_ARRAYS
subroutine RESET_BASE_ARRAYS()
    if (allocated(BASE%Zc_Base)) BASE%Zc_Base = Missing_Value_Real4
    if (allocated(BASE%Pc_Base)) BASE%Pc_Base = Missing_Value_Real4
    if (allocated(BASE%Tc_Base)) BASE%Tc_Base = Missing_Value_Real4
    if (allocated(BASE%Geo_Thickness)) BASE%Geo_Thickness = Missing_Value_Real4
    if (allocated(BASE%Base_Alt)) BASE%Base_Alt = Missing_Value_Real4
    if (allocated(BASE%Lower_Base_Alt)) BASE%Lower_Base_Alt = Missing_Value_Real4
    if (allocated(BASE%Lower_Base_Pc)) BASE%Lower_Base_Pc = Missing_Value_Real4
    if (allocated(BASE%Base_Quality_Flag)) BASE%Base_Quality_Flag = Missing_Value_Int1
    if (allocated(BASE%Lower_Base_Pc_Uncertainty)) BASE%Lower_Base_Pc_Uncertainty = Missing_Value_Real4
end subroutine RESET_BASE_ARRAYS
subroutine DESTROY_BASE_ARRAYS()
    if (allocated(BASE%Zc_Base)) deallocate(BASE%Zc_Base)
    if (allocated(BASE%Pc_Base)) deallocate(BASE%Pc_Base)
    if (allocated(BASE%Tc_Base)) deallocate(BASE%Tc_Base)
    if (allocated(BASE%Geo_Thickness)) deallocate(BASE%Geo_Thickness)
    if (allocated(BASE%Base_Alt)) deallocate(BASE%Base_Alt)
    if (allocated(BASE%Lower_Base_Alt)) deallocate(BASE%Lower_Base_Alt)
    if (allocated(BASE%Lower_Base_Pc)) deallocate(BASE%Lower_Base_Pc)
    if (allocated(BASE%Base_Quality_Flag)) deallocate(BASE%Base_Quality_Flag)
    if (allocated(BASE%Lower_Base_Pc_Uncertainty)) deallocate(BASE%Lower_Base_Pc_Uncertainty)
end subroutine DESTROY_BASE_ARRAYS
!------------------------------------------------------------------------------
! Cloud Cover Layers (CCL) data structure routines
!------------------------------------------------------------------------------
subroutine CREATE_CCL_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (CCL%Mode > 0) then
    allocate (CCL%Cloud_Layer(dim1,dim2))
    allocate (CCL%Cloud_Fraction(dim1,dim2))
    allocate (CCL%Cloud_Fraction_Uncer(dim1,dim2))
    allocate (CCL%Cloud_Fraction_Layer1(dim1,dim2))
    allocate (CCL%Cloud_Fraction_Layer2(dim1,dim2))
    allocate (CCL%Cloud_Fraction_Layer3(dim1,dim2))
    allocate (CCL%Supercooled_Cloud_Layer(dim1,dim2))
    allocate (CCL%Supercooled_Cloud_Fraction(dim1,dim2))
    allocate (CCL%Supercooled_Cloud_Fraction_Layer1(dim1,dim2))
    allocate (CCL%Supercooled_Cloud_Fraction_Layer2(dim1,dim2))
    allocate (CCL%Supercooled_Cloud_Fraction_Layer3(dim1,dim2))
    allocate (CCL%Conv_Cloud_Layer(dim1,dim2))
    allocate (CCL%Conv_Cloud_Fraction(dim1,dim2))
    allocate (CCL%Conv_Cloud_Fraction_Layer1(dim1,dim2))
    allocate (CCL%Conv_Cloud_Fraction_Layer2(dim1,dim2))
    allocate (CCL%Conv_Cloud_Fraction_Layer3(dim1,dim2))
    if (CCL%Type == 0) then
       allocate (CCL%Cloud_Fraction_Layer4(dim1,dim2))
       allocate (CCL%Cloud_Fraction_Layer5(dim1,dim2))
       allocate (CCL%Supercooled_Cloud_Fraction_Layer4(dim1,dim2))
       allocate (CCL%Supercooled_Cloud_Fraction_Layer5(dim1,dim2))
       allocate (CCL%Conv_Cloud_Fraction_Layer4(dim1,dim2))
       allocate (CCL%Conv_Cloud_Fraction_Layer5(dim1,dim2))
    endif
   endif
    allocate (CCL%QF(dim1,dim2))
end subroutine CREATE_CCL_ARRAYS
subroutine RESET_CCL_ARRAYS()
   if (CCL%Mode > 0) then
    if (allocated(CCL%Cloud_Layer)) CCL%Cloud_Layer = Missing_Value_Int1
    if (allocated(CCL%Cloud_Fraction)) CCL%Cloud_Fraction = Missing_Value_Real4
    if (allocated(CCL%Cloud_Fraction_Uncer)) CCL%Cloud_Fraction_Uncer = Missing_Value_Real4
    if (allocated(CCL%Cloud_Fraction_Layer1)) CCL%Cloud_Fraction_Layer1 = MISSING_VALUE_INT1
    if (allocated(CCL%Cloud_Fraction_Layer2)) CCL%Cloud_Fraction_Layer2 = MISSING_VALUE_INT1
    if (allocated(CCL%Cloud_Fraction_Layer3)) CCL%Cloud_Fraction_Layer3 = MISSING_VALUE_INT1
    if (allocated(CCL%Supercooled_Cloud_Layer)) CCL%Supercooled_Cloud_Layer = Missing_Value_Int1
    if (allocated(CCL%Supercooled_Cloud_Fraction)) CCL%Supercooled_Cloud_Fraction = Missing_Value_Real4
    if (allocated(CCL%Supercooled_Cloud_Fraction_Layer1)) CCL%Supercooled_Cloud_Fraction_Layer1 = Missing_Value_Real4
    if (allocated(CCL%Supercooled_Cloud_Fraction_Layer2)) CCL%Supercooled_Cloud_Fraction_Layer2 = Missing_Value_Real4
    if (allocated(CCL%Supercooled_Cloud_Fraction_Layer3)) CCL%Supercooled_Cloud_Fraction_Layer3 = Missing_Value_Real4
    if (allocated(CCL%Conv_Cloud_Layer)) CCL%Conv_Cloud_Layer = Missing_Value_Int1
    if (allocated(CCL%Conv_Cloud_Fraction)) CCL%Conv_Cloud_Fraction = Missing_Value_Real4
    if (allocated(CCL%Conv_Cloud_Fraction_Layer1)) CCL%Conv_Cloud_Fraction_Layer1 = Missing_Value_Real4
    if (allocated(CCL%Conv_Cloud_Fraction_Layer2)) CCL%Conv_Cloud_Fraction_Layer2 = Missing_Value_Real4
    if (allocated(CCL%Conv_Cloud_Fraction_Layer3)) CCL%Conv_Cloud_Fraction_Layer3 = Missing_Value_Real4
    if (CCL%Type == 0) then
      if (allocated(CCL%Cloud_Fraction_Layer4)) CCL%Cloud_Fraction_Layer4 = MISSING_VALUE_INT1
      if (allocated(CCL%Cloud_Fraction_Layer5)) CCL%Cloud_Fraction_Layer5 = MISSING_VALUE_INT1
      if (allocated(CCL%Supercooled_Cloud_Fraction_Layer4)) CCL%Supercooled_Cloud_Fraction_Layer4 = Missing_Value_Real4
      if (allocated(CCL%Supercooled_Cloud_Fraction_Layer5)) CCL%Supercooled_Cloud_Fraction_Layer5 = Missing_Value_Real4
      if (allocated(CCL%Conv_Cloud_Fraction_Layer4)) CCL%Conv_Cloud_Fraction_Layer4 = Missing_Value_Real4
      if (allocated(CCL%Conv_Cloud_Fraction_Layer5)) CCL%Conv_Cloud_Fraction_Layer5 = Missing_Value_Real4
    endif
   endif
    if (allocated(CCL%QF)) CCL%QF = MISSING_VALUE_INT1
end subroutine RESET_CCL_ARRAYS
subroutine DESTROY_CCL_ARRAYS()
    if (allocated(CCL%Cloud_Layer)) deallocate (CCL%Cloud_Layer)
    if (allocated(CCL%Cloud_Fraction)) deallocate (CCL%Cloud_Fraction)
    if (allocated(CCL%Cloud_Fraction_Uncer)) deallocate (CCL%Cloud_Fraction_Uncer)
    if (allocated(CCL%Cloud_Fraction_Layer1)) deallocate(CCL%Cloud_Fraction_Layer1)
    if (allocated(CCL%Cloud_Fraction_Layer2)) deallocate(CCL%Cloud_Fraction_Layer2)
    if (allocated(CCL%Cloud_Fraction_Layer3)) deallocate(CCL%Cloud_Fraction_Layer3)
    if (allocated(CCL%Cloud_Fraction_Layer4)) deallocate(CCL%Cloud_Fraction_Layer4)
    if (allocated(CCL%Cloud_Fraction_Layer5)) deallocate(CCL%Cloud_Fraction_Layer5)
    if (allocated(CCL%Supercooled_Cloud_Layer)) deallocate (CCL%Supercooled_Cloud_Layer)
    if (allocated(CCL%Supercooled_Cloud_Fraction)) deallocate (CCL%Supercooled_Cloud_Fraction)
    if (allocated(CCL%Supercooled_Cloud_Fraction_Layer1)) deallocate(CCL%Supercooled_Cloud_Fraction_Layer1)
    if (allocated(CCL%Supercooled_Cloud_Fraction_Layer2)) deallocate(CCL%Supercooled_Cloud_Fraction_Layer2)
    if (allocated(CCL%Supercooled_Cloud_Fraction_Layer3)) deallocate(CCL%Supercooled_Cloud_Fraction_Layer3)
    if (allocated(CCL%Supercooled_Cloud_Fraction_Layer4)) deallocate(CCL%Supercooled_Cloud_Fraction_Layer4)
    if (allocated(CCL%Supercooled_Cloud_Fraction_Layer5)) deallocate(CCL%Supercooled_Cloud_Fraction_Layer5)
    if (allocated(CCL%Conv_Cloud_Layer)) deallocate (CCL%Conv_Cloud_Layer)
    if (allocated(CCL%Conv_Cloud_Fraction)) deallocate (CCL%Conv_Cloud_Fraction)
    if (allocated(CCL%Conv_Cloud_Fraction_Layer1)) deallocate(CCL%Conv_Cloud_Fraction_Layer1)
    if (allocated(CCL%Conv_Cloud_Fraction_Layer2)) deallocate(CCL%Conv_Cloud_Fraction_Layer2)
    if (allocated(CCL%Conv_Cloud_Fraction_Layer3)) deallocate(CCL%Conv_Cloud_Fraction_Layer3)
    if (allocated(CCL%Conv_Cloud_Fraction_Layer4)) deallocate(CCL%Conv_Cloud_Fraction_Layer4)
    if (allocated(CCL%Conv_Cloud_Fraction_Layer5)) deallocate(CCL%Conv_Cloud_Fraction_Layer5)
    if (allocated(CCL%QF)) deallocate(CCL%QF)
end subroutine DESTROY_CCL_ARRAYS
!------------------------------------------------------------------------------
! Automatic Surface Observing System (ASOS) data structure routines
!------------------------------------------------------------------------------
subroutine CREATE_ASOS_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
    allocate(ASOS%Code(dim1,dim2))
    allocate(ASOS%ECA(dim1,dim2))
    allocate(ASOS%Zmax(dim1,dim2))
    allocate(ASOS%Zmin(dim1,dim2))
   endif
end subroutine CREATE_ASOS_ARRAYS
subroutine RESET_ASOS_ARRAYS()
   if (allocated(ASOS%Code)) ASOS%Code = Missing_Value_Int1
   if (allocated(ASOS%ECA)) ASOS%ECA = Missing_Value_Real4
   if (allocated(ASOS%Zmax)) ASOS%Zmax = Missing_Value_Real4
   if (allocated(ASOS%Zmin)) ASOS%Zmin = Missing_Value_Real4
end subroutine RESET_ASOS_ARRAYS
subroutine DESTROY_ASOS_ARRAYS()
   if (allocated(ASOS%Code)) deallocate(ASOS%Code)
   if (allocated(ASOS%ECA)) deallocate(ASOS%ECA)
   if (allocated(ASOS%Zmax)) deallocate(ASOS%Zmax)
   if (allocated(ASOS%Zmin)) deallocate(ASOS%Zmin)
end subroutine DESTROY_ASOS_ARRAYS

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_DCOMP_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
      allocate(Tau_DCOMP_1(dim1,dim2))
      allocate(Tau_DCOMP_2(dim1,dim2))
      allocate(Tau_DCOMP_3(dim1,dim2))
      allocate(Tau_DCOMP(dim1,dim2))
      allocate(Tau_Aux(dim1,dim2))
      allocate(Reff_Aux(dim1,dim2))
      allocate(Tau_DCOMP_Ap(dim1,dim2))
      allocate(Vis_Ref_Fm(dim1,dim2))
      allocate(Reff_DCOMP(dim1,dim2))
      allocate(Reff_DCOMP_1(dim1,dim2))
      allocate(Reff_DCOMP_2(dim1,dim2))
      allocate(Reff_DCOMP_3(dim1,dim2))
      allocate(Lwp_DCOMP(dim1,dim2))
      allocate (refl_asym_dcomp(dim1,dim2))
      allocate(Iwp_DCOMP(dim1,dim2))
      allocate(Iwp_Tau_DCOMP(dim1,dim2))
      allocate(Cwp_DCOMP(dim1,dim2))
      allocate(Cwp_Fit(dim1,dim2))
      allocate(Reff_DCOMP_Fit(dim1,dim2))
      allocate(Cwp_Ice_Layer_DCOMP(dim1,dim2))
      allocate(Cwp_Water_Layer_DCOMP(dim1,dim2))
      allocate(Cwp_Scwater_Layer_DCOMP(dim1,dim2))
      allocate(Iwc_DCOMP(dim1,dim2))
      allocate(Lwc_DCOMP(dim1,dim2))
      allocate(Rain_Rate_DCOMP(dim1,dim2))
      allocate(Hcld_DCOMP(dim1,dim2))
      allocate(Cdnc_DCOMP(dim1,dim2))
      allocate(Tau_DCOMP_Cost(dim1,dim2))
      allocate(Reff_DCOMP_Cost(dim1,dim2))
      allocate(Tau_DCOMP_Qf(dim1,dim2))
      allocate(Reff_DCOMP_Qf(dim1,dim2))
      allocate(DCOMP_Quality_Flag(dim1,dim2))
      allocate(DCOMP_Info_Flag(dim1,dim2))
      allocate(Cloud_VIS_Albedo(dim1,dim2))
      allocate(refl_vis_fm_dcomp(dim1,dim2))
      allocate(refl_nir_fm_dcomp(dim1,dim2))
      allocate(Cloud_VIS_Spherical_Albedo(dim1,dim2))
      allocate(Above_Cloud_VIS_Transmission(dim1,dim2))
      allocate(Cloud_VIS_Transmission_View(dim1,dim2))
      allocate(Cloud_VIS_Transmission_Solar(dim1,dim2))
      allocate(Cloud_VIS_Transmission_AC(dim1,dim2))
      allocate(Cloud_NIR_Spherical_Albedo(dim1,dim2))
      allocate(Above_Cloud_NIR_Transmission(dim1,dim2))
      allocate(Cloud_NIR_Transmission_View(dim1,dim2))
      allocate(Cloud_NIR_Transmission_Solar(dim1,dim2))
      allocate(Cloud_NIR_Transmission_AC(dim1,dim2))
      allocate(Insolation_DCOMP(dim1,dim2))
      allocate(Insolation_Diffuse_DCOMP(dim1,dim2))
      allocate(Cost_DCOMP(dim1,dim2))
      allocate(Error_Cov_Matrix_Cod(dim1,dim2))
      allocate(Error_Cov_Matrix_Ref(dim1,dim2))
      allocate(DCOMP_Diag_1(dim1,dim2))
      allocate(DCOMP_Diag_2(dim1,dim2))
      allocate(DCOMP_Diag_3(dim1,dim2))
      allocate(DCOMP_Diag_4(dim1,dim2))
      allocate(DCOMP_Diag_Wv1(dim1,dim2))
      allocate(DCOMP_Diag_Wv2(dim1,dim2))
      allocate(DCOMP_Diag_Toc_Rfl1(dim1,dim2))
      allocate(DCOMP_Diag_Toc_Rfl2(dim1,dim2))
      allocate(DCOMP_Diag_Virt_Alb1(dim1,dim2))
      allocate(DCOMP_Diag_Virt_Alb2(dim1,dim2))
      allocate(DCOMP_Diag_Toc_Rfl_Unc1(dim1,dim2))
      allocate(DCOMP_Diag_Toc_Rfl_Unc2(dim1,dim2))
   endif
end subroutine CREATE_DCOMP_ARRAYS
subroutine RESET_DCOMP_ARRAYS()
   if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
      Tau_DCOMP = Missing_Value_Real4
      Tau_DCOMP_1 = Missing_Value_Real4
      Tau_DCOMP_2 = Missing_Value_Real4
      Tau_DCOMP_3 = Missing_Value_Real4
      Tau_Aux = Missing_Value_Real4
      Reff_Aux = Missing_Value_Real4
      Tau_DCOMP_Ap = Missing_Value_Real4
      Vis_Ref_Fm = Missing_Value_Real4
      Reff_DCOMP = Missing_Value_Real4
      Reff_DCOMP_1 = Missing_Value_Real4
      Reff_DCOMP_2 = Missing_Value_Real4
      Reff_DCOMP_3 = Missing_Value_Real4
      Lwp_DCOMP = Missing_Value_Real4
      refl_asym_dcomp = Missing_Value_Real4
      Iwp_DCOMP = Missing_Value_Real4
      Iwp_Tau_DCOMP = Missing_Value_Real4
      Cwp_DCOMP = Missing_Value_Real4
      Cwp_Fit = Missing_Value_Real4
      Reff_DCOMP_Fit = Missing_Value_Real4
      Cwp_Ice_Layer_DCOMP = Missing_Value_Real4
      Cwp_Water_Layer_DCOMP = Missing_Value_Real4
      Cwp_Scwater_Layer_DCOMP = Missing_Value_Real4
      Iwc_DCOMP = Missing_Value_Real4
      Lwc_DCOMP = Missing_Value_Real4
      Rain_Rate_DCOMP = Missing_Value_Real4
      Hcld_DCOMP = Missing_Value_Real4
      Cdnc_DCOMP = Missing_Value_Real4
      Tau_DCOMP_Cost = Missing_Value_Real4
      Reff_DCOMP_Cost = Missing_Value_Real4
      Tau_DCOMP_Qf = Missing_Value_Int1
      Reff_DCOMP_Qf = Missing_Value_Int1
      DCOMP_Quality_Flag = 0
      DCOMP_Info_Flag = 0
      Cloud_VIS_Albedo = Missing_Value_Real4
      refl_vis_fm_dcomp = Missing_Value_Real4
      refl_nir_fm_dcomp = Missing_Value_Real4
      Cloud_VIS_Spherical_Albedo = Missing_Value_Real4
      Above_Cloud_VIS_Transmission = Missing_Value_Real4
      Cloud_VIS_Transmission_View = Missing_Value_Real4
      Cloud_VIS_Transmission_Solar = Missing_Value_Real4
      Cloud_VIS_Transmission_AC = Missing_Value_Real4
      Cloud_NIR_Spherical_Albedo = Missing_Value_Real4
      Above_Cloud_NIR_Transmission = Missing_Value_Real4
      Cloud_NIR_Transmission_View = Missing_Value_Real4
      Cloud_NIR_Transmission_Solar = Missing_Value_Real4
      Cloud_NIR_Transmission_AC = Missing_Value_Real4
      Insolation_DCOMP = Missing_Value_Real4
      Insolation_Diffuse_DCOMP = Missing_Value_Real4
      Cost_DCOMP = Missing_Value_Real4
      Error_Cov_Matrix_Cod = Missing_Value_Real4
      Error_Cov_Matrix_Ref = Missing_Value_Real4
      DCOMP_Diag_1 = Missing_Value_Real4
      DCOMP_Diag_2 = Missing_Value_Real4
      DCOMP_Diag_3 = Missing_Value_Real4
      DCOMP_Diag_4 = Missing_Value_Real4
      DCOMP_Diag_Wv1 = Missing_Value_Real4
      DCOMP_Diag_Wv2 = Missing_Value_Real4
      DCOMP_Diag_Toc_Rfl1 = Missing_Value_Real4
      DCOMP_Diag_Toc_Rfl2 = Missing_Value_Real4
      DCOMP_Diag_Virt_Alb1 = Missing_Value_Real4
      DCOMP_Diag_Virt_Alb2 = Missing_Value_Real4
      DCOMP_Diag_Toc_Rfl_Unc1 = Missing_Value_Real4
      DCOMP_Diag_Toc_Rfl_Unc2 = Missing_Value_Real4
   endif
end subroutine RESET_DCOMP_ARRAYS
subroutine DESTROY_DCOMP_ARRAYS()
   if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
      deallocate(Tau_DCOMP)
      deallocate(Tau_DCOMP_1)
      deallocate(Tau_DCOMP_2)
      deallocate(Tau_DCOMP_3)
      deallocate(Tau_Aux)
      deallocate(Reff_Aux)
      deallocate(Tau_DCOMP_Ap)
      deallocate(Vis_Ref_Fm)
      deallocate(Reff_DCOMP)
      deallocate(Reff_DCOMP_1)
      deallocate(Reff_DCOMP_2)
      deallocate(Reff_DCOMP_3)
      deallocate(Lwp_DCOMP)
      deallocate ( refl_asym_dcomp)
      deallocate(Iwp_DCOMP)
      deallocate(Iwp_Tau_DCOMP)
      deallocate(Cwp_DCOMP)
      deallocate(Cwp_Fit)
      deallocate(Reff_DCOMP_Fit)
      deallocate(Cwp_Ice_Layer_DCOMP)
      deallocate(Cwp_Water_Layer_DCOMP)
      deallocate(Cwp_Scwater_Layer_DCOMP)
      deallocate(Iwc_DCOMP)
      deallocate(Lwc_DCOMP)
      deallocate(Rain_Rate_DCOMP)
      deallocate(Hcld_DCOMP)
      deallocate(Cdnc_DCOMP)
      deallocate(Tau_DCOMP_Cost)
      deallocate(Reff_DCOMP_Cost)
      deallocate(Tau_DCOMP_Qf)
      deallocate(Reff_DCOMP_Qf)
      deallocate(DCOMP_Quality_Flag)
      deallocate(DCOMP_Info_Flag)
      deallocate(Cloud_VIS_Albedo)
      deallocate(refl_vis_fm_dcomp)
      deallocate(refl_nir_fm_dcomp)
      deallocate(Cloud_VIS_Spherical_Albedo)
      deallocate(Above_Cloud_VIS_Transmission)
      deallocate(Cloud_VIS_Transmission_View)
      deallocate(Cloud_VIS_Transmission_Solar)
      deallocate(Cloud_VIS_Transmission_AC)
      deallocate(Cloud_NIR_Spherical_Albedo)
      deallocate(Above_Cloud_NIR_Transmission)
      deallocate(Cloud_NIR_Transmission_View)
      deallocate(Cloud_NIR_Transmission_Solar)
      deallocate(Cloud_NIR_Transmission_AC)
      deallocate(Insolation_DCOMP)
      deallocate(Insolation_Diffuse_DCOMP)
      deallocate(Cost_DCOMP)
      deallocate(Error_Cov_Matrix_Cod)
      deallocate(Error_Cov_Matrix_Ref)
      deallocate(DCOMP_Diag_1)
      deallocate(DCOMP_Diag_2)
      deallocate(DCOMP_Diag_3)
      deallocate(DCOMP_Diag_4)
      deallocate(DCOMP_Diag_Wv1)
      deallocate(DCOMP_Diag_Wv2)
      deallocate(DCOMP_Diag_Toc_Rfl1)
      deallocate(DCOMP_Diag_Toc_Rfl2)
      deallocate(DCOMP_Diag_Virt_Alb1)
      deallocate(DCOMP_Diag_Virt_Alb2)
      deallocate(DCOMP_Diag_Toc_Rfl_Unc1)
      deallocate(DCOMP_Diag_Toc_Rfl_Unc2)
   endif
end subroutine DESTROY_DCOMP_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_SASRAB_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Sasrab_Flag == sym%YES) then
      allocate(Insolation_All_Sky(dim1,dim2))
      allocate(Insolation_All_Sky_Diffuse(dim1,dim2))
      allocate(Insolation_Clear_Sky(dim1,dim2))
      allocate(Insolation_Cld_Opd(dim1,dim2))
      allocate(Insolation_Aer_Opd(dim1,dim2))
   endif
end subroutine CREATE_SASRAB_ARRAYS
subroutine RESET_SASRAB_ARRAYS()
   if (Sasrab_Flag == sym%YES) then
      Insolation_All_Sky = Missing_Value_Real4
      Insolation_All_Sky_Diffuse = Missing_Value_Real4
      Insolation_Clear_Sky = Missing_Value_Real4
      Insolation_Cld_Opd = Missing_Value_Real4
      Insolation_Aer_Opd = Missing_Value_Real4
   endif
end subroutine RESET_SASRAB_ARRAYS
subroutine DESTROY_SASRAB_ARRAYS()
   if (sasrab_Flag == sym%YES) then
      deallocate(Insolation_All_Sky)
      deallocate(Insolation_All_Sky_Diffuse)
      deallocate(Insolation_Clear_Sky)
      deallocate(Insolation_Cld_Opd)
      deallocate(Insolation_Aer_Opd)
   endif
end subroutine DESTROY_SASRAB_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_OLR_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2

      allocate(Olr(dim1,dim2))

end subroutine CREATE_OLR_ARRAYS
subroutine RESET_OLR_ARRAYS

      Olr = Missing_Value_Real4

end subroutine RESET_OLR_ARRAYS
subroutine DESTROY_OLR_ARRAYS

      deallocate(Olr)

end subroutine DESTROY_OLR_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_AEROSOL_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  allocate(Aot_Qf(dim1,dim2))
  if (Aerosol_Mode > 0) then
     allocate(Aot1(dim1,dim2))
     allocate(Aot2(dim1,dim2))
     allocate(Aot3a(dim1,dim2))
  endif
end subroutine CREATE_AEROSOL_ARRAYS
subroutine RESET_AEROSOL_ARRAYS()
   Aot_Qf = 0
   if (Aerosol_Mode > 0) then
    Aot1 = Missing_Value_Real4
    Aot2 = Missing_Value_Real4
    Aot3a = Missing_Value_Real4
   endif
end subroutine RESET_AEROSOL_ARRAYS
subroutine DESTROY_AEROSOL_ARRAYS()
  deallocate(Aot_Qf)
  if (Aerosol_Mode > 0) then
     deallocate(Aot1)
     deallocate(Aot2)
     deallocate(Aot3a)
  endif
end subroutine DESTROY_AEROSOL_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_CLOUD_MASK_ARRAYS(dim1,dim2,dim3)
  integer, intent(in):: dim1, dim2, dim3
  allocate(CLDMASK%Cld_Mask_Qf(dim1,dim2))
  if (Cld_Flag == sym%YES .OR. Tracer_Flag==1) then
     allocate(CLDMASK%Cld_Mask(dim1,dim2))
     allocate(CLDMASK%Cld_Mask_Binary(dim1,dim2))
     allocate(CLDMASK%Cld_Mask_IR(dim1,dim2))
     allocate(CLDMASK%Cld_Mask_Binary_IR(dim1,dim2))
     allocate(CLDMASK%Cld_Mask_Aux(dim1,dim2))
     allocate(CLDMASK%Adj_Pix_Cld_Mask(dim1,dim2))
     allocate(CLDMASK%TUT(dim1,dim2))
     allocate(CLDMASK%RUT(dim1,dim2))
     allocate(CLDMASK%Posterior_Cld_Probability(dim1,dim2))
     allocate(CLDMASK%Posterior_Cld_Probability_Uncer(dim1,dim2))
     allocate(CLDMASK%Posterior_Cld_Probability_IR(dim1,dim2))
     allocate(CLDMASK%Posterior_Cld_Probability_Aux(dim1,dim2))
     allocate(CLDMASK%Posterior_Water_Probability(dim1,dim2))
     allocate(CLDMASK%Posterior_Ice_Probability(dim1,dim2))
     allocate(CLDMASK%Prior_Cld_Probability(dim1,dim2))
     allocate(CLDMASK%Bayes_Mask_Sfc_Type(dim1,dim2))
     allocate(CLDMASK%Cld_Test_Vector_Packed(dim3,dim1,dim2))
     allocate(CLDMASK%Shadow_Mask(dim1,dim2))
     allocate(CLDMASK%Dust_Mask(dim1,dim2))
     allocate(CLDMASK%Smoke_Mask(dim1,dim2))
     allocate(CLDMASK%Fire_Mask(dim1,dim2))
     allocate(CLDMASK%Thin_Cirr_Mask(dim1,dim2))
     allocate(CLDMASK%Dust_Prob(dim1,dim2))
  endif
end subroutine CREATE_CLOUD_MASK_ARRAYS
subroutine RESET_CLOUD_MASK_ARRAYS()
  if (allocated(CLDMASK%Cld_Mask_Qf)) CLDMASK%Cld_Mask_Qf = Missing_Value_Int1
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
     CLDMASK%Cld_Mask = Missing_Value_Int1
     CLDMASK%Cld_Mask_Binary = Missing_Value_Int1
     CLDMASK%Cld_Mask_IR = Missing_Value_Int1
     CLDMASK%Cld_Mask_Binary_IR = Missing_Value_Int1
     CLDMASK%Cld_Mask_Aux = Missing_Value_Int1
     CLDMASK%Adj_Pix_Cld_Mask = Missing_Value_Int1
     CLDMASK%TUT = Missing_Value_Int1
     CLDMASK%RUT = Missing_Value_Int1
     CLDMASK%Prior_Cld_Probability = Missing_Value_Real4
     CLDMASK%Posterior_Cld_Probability = Missing_Value_Real4
     CLDMASK%Posterior_Cld_Probability_Uncer = Missing_Value_Real4
     CLDMASK%Posterior_Cld_Probability_IR = Missing_Value_Real4
     CLDMASK%Posterior_Cld_Probability_Aux = Missing_Value_Real4
     CLDMASK%Posterior_Ice_Probability = Missing_Value_Real4
     CLDMASK%Posterior_Water_Probability = Missing_Value_Real4
     CLDMASK%Cld_Test_Vector_Packed = 0
     CLDMASK%Bayes_Mask_Sfc_Type = Missing_Value_Int1
     CLDMASK%Shadow_Mask = Missing_Value_Int1
     CLDMASK%Dust_Mask = Missing_Value_Int1
     CLDMASK%Smoke_Mask = Missing_Value_Int1
     CLDMASK%Fire_Mask = Missing_Value_Int1
     CLDMASK%Thin_Cirr_Mask = Missing_Value_Int1
     CLDMASK%Dust_Prob = Missing_Value_Real4
  endif
end subroutine RESET_CLOUD_MASK_ARRAYS
subroutine DESTROY_CLOUD_MASK_ARRAYS()
  deallocate(CLDMASK%Cld_Mask_Qf)
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
     deallocate(CLDMASK%Cld_Mask)
     deallocate(CLDMASK%Cld_Mask_Binary)
     deallocate(CLDMASK%Cld_Mask_IR)
     deallocate(CLDMASK%Cld_Mask_Binary_IR)
     deallocate(CLDMASK%Cld_Mask_Aux)
     deallocate(CLDMASK%Adj_Pix_Cld_Mask)
     deallocate(CLDMASK%TUT)
     deallocate(CLDMASK%RUT)
     deallocate(CLDMASK%Posterior_Cld_Probability)
     deallocate(CLDMASK%Posterior_Cld_Probability_Uncer)
     deallocate(CLDMASK%Posterior_Cld_Probability_IR)
     deallocate(CLDMASK%Posterior_Cld_Probability_Aux)
     deallocate(CLDMASK%Posterior_Water_Probability)
     deallocate(CLDMASK%Posterior_Ice_Probability)
     deallocate(CLDMASK%Prior_Cld_Probability)
     deallocate(CLDMASK%Cld_Test_Vector_Packed)
     deallocate(CLDMASK%Bayes_Mask_Sfc_Type)
     deallocate(CLDMASK%Shadow_Mask)
     deallocate(CLDMASK%Dust_Mask)
     deallocate(CLDMASK%Smoke_Mask)
     deallocate(CLDMASK%Fire_Mask)
     deallocate(CLDMASK%Thin_Cirr_Mask)
     deallocate(CLDMASK%Dust_Prob)
  endif
end subroutine DESTROY_CLOUD_MASK_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_CLOUD_TYPE_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
     allocate(Metadata_Aux(dim1,dim2))
     allocate(Cld_Type_Aux(dim1,dim2))
     allocate(Cld_Phase_Aux(dim1,dim2))
     allocate(Cld_Phase(dim1,dim2))
     allocate(Cld_Phase_Uncertainty(dim1,dim2))
     allocate(Cld_Type(dim1,dim2))
     allocate(Cld_Phase_IR(dim1,dim2))
     allocate(Cld_Type_IR(dim1,dim2))
!    allocate(Ctp_Multilayer_Flag(dim1,dim2))
     allocate(Zc_Aux(dim1,dim2))
     allocate(Ec_Aux(dim1,dim2))
     allocate(Tc_Aux(dim1,dim2))
  endif
end subroutine CREATE_CLOUD_TYPE_ARRAYS
subroutine RESET_CLOUD_TYPE_ARRAYS()
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
      Cld_Phase = Missing_Value_Int1
      Cld_Phase_Uncertainty = Missing_Value_Real4
      Cld_Type = Missing_Value_Int1
      Cld_Phase_IR = Missing_Value_Int1
      Cld_Type_IR = Missing_Value_Int1
      Cld_Phase_Aux = Missing_Value_Int1
      Metadata_Aux = Missing_Value_Int1
      Cld_Type_Aux = Missing_Value_Int1
!     Ctp_Multilayer_Flag = Missing_Value_Int1
      Zc_Aux = Missing_Value_Real4
      Ec_Aux = Missing_Value_Real4
      Tc_Aux = Missing_Value_Real4
  endif
end subroutine RESET_CLOUD_TYPE_ARRAYS
subroutine DESTROY_CLOUD_TYPE_ARRAYS
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
     deallocate(Metadata_Aux)
     deallocate(Cld_Type_Aux)
     deallocate(Cld_Phase_Aux)
     deallocate(Cld_Phase)
     deallocate(Cld_Phase_Uncertainty)
     deallocate(Cld_Type)
     deallocate(Cld_Phase_IR)
     deallocate(Cld_Type_IR)
!    deallocate(Ctp_Multilayer_Flag)
     deallocate(Zc_Aux)
     deallocate(Ec_Aux)
     deallocate(Tc_Aux)
  endif
end subroutine DESTROY_CLOUD_TYPE_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_DIAGNOSTIC_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  allocate(Temp_Pix_Array_1(dim1,dim2))
  allocate(Temp_Pix_Array_2(dim1,dim2))
  allocate(Temp_Pix_Array_3(dim1,dim2))
  allocate(Diag_Pix_Array_1(dim1,dim2))
  allocate(Diag_Pix_Array_2(dim1,dim2))
  allocate(Diag_Pix_Array_3(dim1,dim2))
  allocate(Missing_Pixel_Array_Real4(dim1,dim2))
  allocate(One_Byte_Temp(dim1,dim2))
  allocate(Two_Byte_Temp(dim1,dim2))
end subroutine CREATE_DIAGNOSTIC_ARRAYS
subroutine RESET_DIAGNOSTIC_ARRAYS()
  Temp_Pix_Array_1 = Missing_Value_Real4
  Temp_Pix_Array_2 = Missing_Value_Real4
  Temp_Pix_Array_3 = Missing_Value_Real4
  Diag_Pix_Array_1 = Missing_Value_Real4
  Diag_Pix_Array_2 = Missing_Value_Real4
  Diag_Pix_Array_3 = Missing_Value_Real4
  Missing_Pixel_Array_Real4 = Missing_Value_Real4
  One_Byte_Temp = Missing_Value_Int1
  Two_Byte_Temp = Missing_Value_Int2
end subroutine RESET_DIAGNOSTIC_ARRAYS
subroutine DESTROY_DIAGNOSTIC_ARRAYS()
  deallocate(Temp_Pix_Array_1)
  deallocate(Temp_Pix_Array_2)
  deallocate(Temp_Pix_Array_3)
  deallocate(Diag_Pix_Array_1)
  deallocate(Diag_Pix_Array_2)
  deallocate(Diag_Pix_Array_3)
  deallocate(Missing_Pixel_Array_Real4)
  deallocate(One_Byte_Temp)
  deallocate(Two_Byte_Temp)
end subroutine DESTROY_DIAGNOSTIC_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_SFC_PROD_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  allocate(Tsfc_Retrieved(dim1,dim2))
  allocate(Trad_Retrieved(dim1,dim2))
  allocate(Tsfc_Qf(dim1,dim2))
  allocate(Ndsi_Toa(dim1,dim2))
  allocate(Ndsi_Sfc(dim1,dim2))
  allocate(Ndvi_Toa(dim1,dim2))
  allocate(Ndvi_Qf(dim1,dim2))
  allocate(Ndvi_Sfc(dim1,dim2))
  allocate(Ndvi_Sfc_White_Sky(dim1,dim2))
  allocate(Nddi_Toa(dim1,dim2))
  allocate(Rsr(dim1,dim2))
  allocate(Rsr_Qf(dim1,dim2))
  allocate(Sst_Retrieved(dim1,dim2))
end subroutine CREATE_SFC_PROD_ARRAYS
subroutine RESET_SFC_PROD_ARRAYS()
  Tsfc_Retrieved = Missing_Value_Real4
  Trad_Retrieved = Missing_Value_Real4
  Tsfc_Qf = Missing_Value_Int1
  Ndsi_Toa = Missing_Value_Real4
  Ndsi_Sfc = Missing_Value_Real4
  Ndvi_Toa = Missing_Value_Real4
  Ndvi_Qf = Missing_Value_Int1
  Ndvi_Sfc = Missing_Value_Real4
  Ndvi_Sfc_White_Sky = Missing_Value_Real4
  Nddi_Toa = Missing_Value_Real4
  Rsr = Missing_Value_Real4
  Rsr_Qf = Missing_Value_Int1
  Sst_Retrieved = Missing_Value_Real4
end subroutine RESET_SFC_PROD_ARRAYS
subroutine DESTROY_SFC_PROD_ARRAYS()
  if (allocated(Tsfc_Retrieved))  deallocate(Tsfc_Retrieved)
  if (allocated(Trad_Retrieved))  deallocate(Trad_Retrieved)
  if (allocated(Tsfc_Qf))  deallocate(Tsfc_Qf)
  if (allocated(Ndsi_Toa)) deallocate(Ndsi_Toa)
  if (allocated(Ndsi_Sfc)) deallocate(Ndsi_Sfc)
  if (allocated(Ndvi_Toa)) deallocate(Ndvi_Toa)
  if (allocated(Ndvi_Qf)) deallocate(Ndvi_Qf)
  if (allocated(Ndvi_Sfc)) deallocate(Ndvi_Sfc)
  if (allocated(Ndvi_Sfc_White_Sky)) deallocate(Ndvi_Sfc_White_Sky)
  if (allocated(Nddi_Toa)) deallocate(Nddi_Toa)
  if (allocated(Rsr)) deallocate(Rsr)
  if (allocated(Rsr_Qf)) deallocate(Rsr_Qf)
  if (allocated(Sst_Retrieved)) deallocate(Sst_Retrieved)
end subroutine DESTROY_SFC_PROD_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_CLOUD_PROD_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  allocate(Pc_Opaque_Cloud(dim1,dim2))
  allocate(Zc_Opaque_Cloud(dim1,dim2))
  allocate(Tc_Opaque_Cloud(dim1,dim2))
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
    allocate(Pc_H2O(dim1,dim2))
    allocate(Tc_H2O(dim1,dim2))
    allocate(Zc_H2O(dim1,dim2))
    allocate(Ec_H2O(dim1,dim2))
    allocate(Zclr_H2O_Peak(dim1,dim2))
    allocate(Zc_CO2IRW(dim1,dim2))
    allocate(Pc_CO2IRW(dim1,dim2))
    allocate(Tc_CO2IRW(dim1,dim2))
    allocate(Ec_CO2IRW(dim1,dim2))
    allocate(Zc_SplitWin(dim1,dim2))
    allocate(Pc_SplitWin(dim1,dim2))
    allocate(Tc_SplitWin(dim1,dim2))
    allocate(Ec_SplitWin(dim1,dim2))
    allocate(Tc_Co2(dim1,dim2))
    allocate(Pc_Co2(dim1,dim2))
    allocate(Zc_Co2(dim1,dim2))
    allocate(Ec_Co2(dim1,dim2))
    allocate(Tc_Cirrus_Background(dim1,dim2))
    allocate(Zc_Cirrus_Background(dim1,dim2))
    allocate(Cloud_Fraction_Background(dim1,dim2))
  endif
end subroutine CREATE_CLOUD_PROD_ARRAYS
subroutine RESET_CLOUD_PROD_ARRAYS()
  Pc_Opaque_Cloud = Missing_Value_Real4
  Zc_Opaque_Cloud = Missing_Value_Real4
  Tc_Opaque_Cloud = Missing_Value_Real4
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
     Pc_H2O = Missing_Value_Real4
     Tc_H2O = Missing_Value_Real4
     Zc_H2O = Missing_Value_Real4
     Ec_H2O = Missing_Value_Real4
     Zclr_H2O_Peak = Missing_Value_Real4
     Pc_CO2IRW = Missing_Value_Real4
     Tc_CO2IRW = Missing_Value_Real4
     Zc_CO2IRW = Missing_Value_Real4
     Ec_CO2IRW = Missing_Value_Real4
     Pc_SplitWin = Missing_Value_Real4
     Zc_SplitWin = Missing_Value_Real4
     Tc_SplitWin = Missing_Value_Real4
     Ec_SplitWin = Missing_Value_Real4
     Tc_Co2 = Missing_Value_Real4
     Pc_Co2 = Missing_Value_Real4
     Zc_Co2 = Missing_Value_Real4
     Ec_Co2 = Missing_Value_Real4
     Tc_Cirrus_Background = Missing_Value_Real4
     Zc_Cirrus_Background = Missing_Value_Real4
     Cloud_Fraction_Background = Missing_Value_Real4
  endif
end subroutine RESET_CLOUD_PROD_ARRAYS
subroutine DESTROY_CLOUD_PROD_ARRAYS()
  deallocate(Pc_Opaque_Cloud)
  deallocate(Zc_Opaque_Cloud)
  deallocate(Tc_Opaque_Cloud)
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
     deallocate(Pc_H2O)
     deallocate(Tc_H2O)
     deallocate(Zc_H2O)
     deallocate(Ec_H2O)
     deallocate(Zclr_H2O_Peak)
     deallocate(Zc_CO2IRW)
     deallocate(Tc_CO2IRW)
     deallocate(Pc_CO2IRW)
     deallocate(Ec_CO2IRW)
     deallocate(Zc_SplitWin)
     deallocate(Tc_SplitWin)
     deallocate(Pc_SplitWin)
     deallocate(Ec_SplitWin)
     deallocate(Tc_Co2)
     deallocate(Pc_Co2)
     deallocate(Zc_Co2)
     deallocate(Ec_Co2)
     deallocate(Tc_Cirrus_Background)
     deallocate(Zc_Cirrus_Background)
     deallocate(Cloud_Fraction_Background)
  endif
end subroutine DESTROY_CLOUD_PROD_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_CALIOP_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  if (Caliop_Flag) then
    allocate(Caliop_Num_Cld_Layers(dim1,dim2))
    allocate(Caliop_Cod(dim1,dim2))
    allocate(Caliop_Cld_Height(dim1,dim2))
  endif
end subroutine CREATE_CALIOP_ARRAYS
subroutine RESET_CALIOP_ARRAYS()
  if (Caliop_Flag) then
    Caliop_Num_Cld_Layers = Missing_Value_Real4
    Caliop_Cod = Missing_Value_Real4
    Caliop_Cld_Height = Missing_Value_Real4
  endif
end subroutine RESET_CALIOP_ARRAYS
subroutine DESTROY_CALIOP_ARRAYS()
  if (Caliop_Flag) then
     deallocate(Caliop_Num_Cld_Layers)
     deallocate(Caliop_Cod)
     deallocate(Caliop_Cld_Height)
  endif
end subroutine DESTROY_CALIOP_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_NLCOMP_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
     allocate(Tau_Nlcomp(dim1,dim2))
     allocate(Reff_Nlcomp(dim1,dim2))
     allocate(Nlcomp_Info_Flag(dim1,dim2))
     allocate(Nlcomp_Quality_Flag(dim1,dim2))
     allocate(Tau_Nlcomp_Cost(dim1,dim2))
     allocate(Reff_Nlcomp_Cost(dim1,dim2))
  endif
end subroutine CREATE_NLCOMP_ARRAYS
subroutine RESET_NLCOMP_ARRAYS()
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
      Tau_Nlcomp = Missing_Value_Real4
      Reff_Nlcomp = Missing_Value_Real4
      Tau_Nlcomp_Cost = Missing_Value_Real4
      Reff_Nlcomp_Cost = Missing_Value_Real4
      Nlcomp_Quality_Flag = 0
      Nlcomp_Info_Flag = 0
  endif
end subroutine RESET_NLCOMP_ARRAYS
subroutine DESTROY_NLCOMP_ARRAYS()
  if (Cld_Flag == sym%YES .OR. Tracer_Flag == 1) then
     deallocate(Tau_Nlcomp)
     deallocate(Reff_Nlcomp)
     deallocate(Nlcomp_Info_Flag)
     deallocate(Nlcomp_Quality_Flag)
     deallocate(Tau_Nlcomp_Cost)
     deallocate(Reff_Nlcomp_Cost)
  endif
end subroutine DESTROY_NLCOMP_ARRAYS
!-----------------------------------------------------------
! NUCAPS
!-----------------------------------------------------------
subroutine CREATE_NUCAPS_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  if (Nucaps_Flag) then
    allocate(NUCAPS%Cld_Press_Layer1(dim1,dim2))
    allocate(NUCAPS%Cld_Press_Layer2(dim1,dim2))
    allocate(NUCAPS%Cld_Press(dim1,dim2))
    allocate(NUCAPS%Cld_Fraction_Layer1(dim1,dim2))
    allocate(NUCAPS%Cld_Fraction_Layer2(dim1,dim2))
    allocate(NUCAPS%Cld_Fraction(dim1,dim2))
    allocate(NUCAPS%Cld_Temp(dim1,dim2))
    allocate(NUCAPS%Cld_Height(dim1,dim2))
    allocate(NUCAPS%Cld_Temp_Smoothed(dim1,dim2))
  endif
end subroutine CREATE_NUCAPS_ARRAYS
subroutine RESET_NUCAPS_ARRAYS()
  if (Nucaps_Flag) then
    if (allocated(NUCAPS%Cld_Press_Layer1)) NUCAPS%Cld_Press_Layer1 = Missing_Value_Real4
    if (allocated(NUCAPS%Cld_Press_Layer2)) NUCAPS%Cld_Press_Layer2 = Missing_Value_Real4
    if (allocated(NUCAPS%Cld_Fraction_Layer1)) NUCAPS%Cld_Fraction_Layer1 = Missing_Value_Real4
    if (allocated(NUCAPS%Cld_Fraction_Layer2)) NUCAPS%Cld_Fraction_Layer2 = Missing_Value_Real4
    if (allocated(NUCAPS%Cld_Fraction)) NUCAPS%Cld_Fraction = Missing_Value_Real4
    if (allocated(NUCAPS%Cld_Press)) NUCAPS%Cld_Press = Missing_Value_Real4
    if (allocated(NUCAPS%Cld_Temp)) NUCAPS%Cld_Temp = Missing_Value_Real4
    if (allocated(NUCAPS%Cld_Height)) NUCAPS%Cld_Height = Missing_Value_Real4
    if (allocated(NUCAPS%Cld_Temp_Smoothed)) NUCAPS%Cld_Temp_Smoothed = Missing_Value_Real4
  endif
end subroutine RESET_NUCAPS_ARRAYS
subroutine DESTROY_NUCAPS_ARRAYS()
  if (Nucaps_Flag) then
    if (allocated(NUCAPS%Cld_Press_Layer1)) deallocate(NUCAPS%Cld_Press_Layer1)
    if (allocated(NUCAPS%Cld_Press_Layer2)) deallocate(NUCAPS%Cld_Press_Layer2)
    if (allocated(NUCAPS%Cld_Fraction_Layer1)) deallocate(NUCAPS%Cld_Fraction_Layer1)
    if (allocated(NUCAPS%Cld_Fraction_Layer2)) deallocate(NUCAPS%Cld_Fraction_Layer2)
    if (allocated(NUCAPS%Cld_Fraction)) deallocate(NUCAPS%Cld_Fraction)
    if (allocated(NUCAPS%Cld_Press)) deallocate(NUCAPS%Cld_Press)
    if (allocated(NUCAPS%Cld_Temp)) deallocate(NUCAPS%Cld_Temp)
    if (allocated(NUCAPS%Cld_Height)) deallocate(NUCAPS%Cld_Height)
    if (allocated(NUCAPS%Cld_Temp_Smoothed)) deallocate(NUCAPS%Cld_Temp_Smoothed)
  endif
end subroutine DESTROY_NUCAPS_ARRAYS
!-----------------------------------------------------------
! end of module
!-----------------------------------------------------------
end module PIXEL_COMMON_MOD
