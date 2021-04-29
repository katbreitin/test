! $Id: constants_mod.f90 3825 2020-05-05 14:44:19Z dbotambekov $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: constant.f90 (src)
!       CONSTANTS (program)
!
! PURPOSE: store and serve various constants for use in the CLAVR-x system
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE public
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED public USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! (c) This code is copyrighted by the author and all NOAA restrictions apply
! 
! Reference: CLAVR-x system description document
!
! Dependencies:  None
!
! Calling Sequece:
!   use CONSTANTS_MOD
!
! Public Routines within in this Module: None
!--------------------------------------------------------------------------------------
module CONSTANTS_MOD
  implicit none
  private
  integer, parameter, public:: int1 = selected_int_kind(1)
  integer, parameter, public:: int2 = selected_int_kind(3)
  integer, parameter, public:: int4 = selected_int_kind(8)
  integer, parameter, public:: int8 = selected_int_kind(10)
  integer, parameter, public:: real4 = selected_real_kind(6,37)
  integer, parameter, public:: real8 = selected_real_kind(15,307)
  integer, parameter, public:: ipre = real4

  !--- Common Numerical Constants
  real (kind=real4), parameter, public:: g = 9.8
  real (kind=real4), parameter, public:: pi = 3.14159265
  real (kind=real4), parameter, public:: dtor = pi/180.0

  !--- Missing Values
  real (kind=real4), parameter, public:: Missing_Value_Real4 = -999.0
  real (kind=real8), parameter, public:: Missing_Value_Real8 = -999.0
  integer(kind=int1), parameter, public:: Missing_Value_Int1 = -128_int1
  integer(kind=int2), parameter, public:: Missing_Value_Int2 = -32768_int2
  integer(kind=int4), parameter, public:: Missing_Value_Int4 = -999_int4
  real(kind=real4), parameter, public:: No_Attribute_Missing_Value = -888.0

  !-- constants used for scaling
  integer(kind=int4), parameter, public:: one_byte_max = 127, &   !(2**8)/2 -1
                                         one_byte_min = -127     !-(2**8)/2
  integer(kind=int4), parameter, public:: two_byte_max = 32767, & !(2**15)/2 - 1
                                         two_byte_min = -32767   !-(2**15)/2

  !--- other often used constants
  real(kind=real4), parameter, public:: Day_Solzen_Thresh_Mask = 85.0
  real (kind=real4), parameter, public:: TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH = 60.0
  integer(kind=int4), parameter, public:: Msec_Per_Day = 24*60*60*1000

  !--- parameters used in retrievals
  real, parameter, public::  Cossolzen_Min_Solar =  0.1    !solar angle limit for daytime

  !--- channel numbers
  integer, parameter, public:: Nchan_Avhrr = 6
  integer, parameter, public:: Nchan_Clavrx = 45

  !--- maximum number of cloud mask tests - used to dimension arrays
  integer, parameter, public:: Max_Num_Cld_Tests = 37
  integer, parameter, public:: Max_Num_Cld_Test_Bytes = 7

  !--- error prompt id
  character(*), parameter, public :: EXE_PROMPT = "CLAVR-x>> "

  !--- useful constants
  integer(kind=int4), public, parameter :: MAX_STR_LEN = 256
  integer(kind=int4), public, parameter :: LARGE_HDF_NUMBER = 999999999

  !--- cvs strings to be written as attributes
  character(120), public :: ACHA_Version
  character(120), public :: DCOMP_Version
  character(120), public :: Cloud_Mask_Version
  character(120), public :: Cloud_Mask_Lut_Version
  character(120), public :: Cloud_Type_Version
  character(120), public :: Cloud_Type_IR_Version
  
  !--- define sds names in hdf files of relevant static ancillary data
  character(*), parameter, public :: SFC_TYPE_SDS_NAME = "surface_type"
  character(*), parameter, public :: COAST_MASK_SDS_NAME = "coast_mask"
  character(*), parameter, public :: VOLCANO_MASK_SDS_NAME = "volcano_mask"
  character(*), parameter, public :: LAND_MASK_SDS_NAME = "land_sea_mask"
  character(*), parameter, public :: SURFACE_ELEV_SDS_NAME = "surface_elevation"
  character(*), parameter, public ::  &
                     MODIS_ALB_0_66_SDS_NAME="Albedo_Map_0.659"
  character(*), parameter, public ::  &
                     MODIS_ALB_0_86_SDS_NAME="Albedo_Map_0.858"
  character(*), parameter, public ::  &
                     MODIS_ALB_1_24_SDS_NAME="Albedo_Map_1.24"
  character(*), parameter, public ::  &
                     MODIS_ALB_1_64_SDS_NAME="Albedo_Map_1.64"
  character(*), parameter, public ::  &
                     MODIS_ALB_2_13_SDS_NAME="Albedo_Map_2.13"
  character(*), parameter, public :: SNOW_MASK_SDS_NAME = "snow_ice_cover"

  !--- define a structure of symbols of clarity
  TYPE, public:: symbol_struct
    integer(kind=int1) :: CLOUDY = 3_int1
    integer(kind=int1) :: PROB_CLOUDY = 2_int1
    integer(kind=int1) :: PROB_CLEAR = 1_int1
    integer(kind=int1) :: CLEAR = 0_int1
    
    integer(kind=int1) :: CLEAR_BINARY = 0_int1
    integer(kind=int1) :: CLOUDY_BINARY = 1_int1

    integer(kind=int1) :: CLEAR_TYPE = 0_int1
    integer(kind=int1) :: PROB_CLEAR_TYPE = 1_int1
    integer(kind=int1) :: FOG_TYPE = 2_int1
    integer(kind=int1) :: WATER_TYPE = 3_int1
    integer(kind=int1) :: SUPERCOOLED_TYPE = 4_int1
    integer(kind=int1) :: MIXED_TYPE = 5_int1
    integer(kind=int1) :: OPAQUE_ICE_TYPE = 6_int1
    integer(kind=int1) :: TICE_TYPE = 6_int1
    integer(kind=int1) :: CIRRUS_TYPE = 7_int1
    integer(kind=int1) :: OVERLAP_TYPE = 8_int1
    integer(kind=int1) :: OVERSHOOTING_TYPE = 9_int1
    integer(kind=int1) :: UNKNOWN_TYPE = 10_int1
    integer(kind=int1) :: DUST_TYPE = 11_int1
    integer(kind=int1) :: SMOKE_TYPE = 12_int1
    integer(kind=int1) :: FIRE_TYPE = 13_int1
                                                  
    integer(kind=int1) :: CLEAR_PHASE = 0_int1
    integer(kind=int1) :: WATER_PHASE = 1_int1
    integer(kind=int1) :: SUPERCOOLED_PHASE = 2_int1
    integer(kind=int1) :: MIXED_PHASE = 3_int1
    integer(kind=int1) :: ICE_PHASE = 4_int1
    integer(kind=int1) :: UNKNOWN_PHASE = 5_int1
                                                                                
    integer(kind=int1) :: NO_SPACE = 0_int1
    integer(kind=int1) :: SPACE = 1_int1
                                                                                
    integer(kind=int1) :: NO = 0_int1
    integer(kind=int1) :: YES = 1_int1
    integer(kind=int1) :: ECM1 = 1_int1 !needed for secondary ECM
    integer(kind=int1) :: ECM2 = 2_int1 !needed for secondary ECM

    integer(kind=int1) :: NO_AUX = 0_int1
    integer(kind=int1) :: USE_AUX = 1_int1
    integer(kind=int1) :: READ_BUT_DO_NOT_USE_AUX = 2_int1
    integer(kind=int1) :: USE_AUX_MODAWG = 3_int1
                                                                                
    !--- this apply to the sfc_type array
    integer(kind=int1) :: WATER_SFC = 0_int1
    integer(kind=int1) :: EVERGREEN_NEEDLE_SFC = 1_int1
    integer(kind=int1) :: EVERGREEN_BROAD_SFC = 2_int1
    integer(kind=int1) :: DECIDUOUS_NEEDLE_SFC = 3_int1
    integer(kind=int1) :: DECIDUOUS_BROAD_SFC = 4_int1
    integer(kind=int1) :: MIXED_FORESTS_SFC = 5_int1
    integer(kind=int1) :: WOODLANDS_SFC = 6_int1
    integer(kind=int1) :: WOODED_GRASS_SFC = 7_int1
    integer(kind=int1) :: CLOSED_SHRUBS_SFC = 8_int1
    integer(kind=int1) :: OPEN_SHRUBS_SFC = 9_int1
    integer(kind=int1) :: GRASSES_SFC = 10_int1
    integer(kind=int1) :: CROPLANDS_SFC = 11_int1
    integer(kind=int1) :: BARE_SFC = 12_int1
    integer(kind=int1) :: URBAN_SFC = 13_int1
                                                                                
    integer(kind=int1) :: NO_DESERT = 0_int1
    integer(kind=int1) :: NIR_DESERT = 1_int1
    integer(kind=int1) :: BRIGHT_DESERT = 2_int1

    !--- this apply to the land flags (land array)
    integer(kind=int1) :: SHALLOW_OCEAN = 0_int1
    integer(kind=int1) :: LAND = 1_int1
    integer(kind=int1) :: COASTLINE = 2_int1
    integer(kind=int1) :: SHALLOW_INLAND_WATER = 3_int1
    integer(kind=int1) :: EPHEMERAL_WATER = 4_int1
    integer(kind=int1) :: DEEP_INLAND_WATER = 5_int1
    integer(kind=int1) :: MODERATE_OCEAN = 6_int1
    integer(kind=int1) :: DEEP_OCEAN = 7_int1

    integer(kind=int1) :: NO_VOLCANO = 0_int1
    integer(kind=int1) :: CLOSE_VOLCANO = 1_int1
    integer(kind=int1) :: VERY_CLOSE_VOLCANO = 2_int1

    integer(kind=int1) :: NO_COAST = 0_int1
    integer(kind=int1) :: COAST_1KM = 1_int1
    integer(kind=int1) :: COAST_2KM = 2_int1
    integer(kind=int1) :: COAST_3KM = 3_int1
    integer(kind=int1) :: COAST_4KM = 4_int1
    integer(kind=int1) :: COAST_5KM = 5_int1
    integer(kind=int1) :: COAST_6KM = 6_int1
    integer(kind=int1) :: COAST_7KM = 7_int1
    integer(kind=int1) :: COAST_8KM = 8_int1
    integer(kind=int1) :: COAST_9KM = 9_int1
    integer(kind=int1) :: COAST_10KM = 10_int1 

    integer(kind=int1) :: WATER_GEN = 0_int1
    integer(kind=int1) :: COAST_GEN = 1_int1
    integer(kind=int1) :: LAND_GEN  = 2_int1
    integer(kind=int1) :: DESERT_GEN = 3_int1
    integer(kind=int1) :: SNOW_GEN   = 4_int1

    integer(kind=int1) :: NO_SNOW = 1_int1
    integer(kind=int1) :: SEA_ICE = 2_int1
    integer(kind=int1) :: SNOW = 3_int1

    integer(kind=int1) :: READ_SNOW_NWP = 0_int1
    integer(kind=int1) :: READ_SNOW_HIRES = 1_int1
    integer(kind=int1) :: READ_SNOW_GLOB = 2_int1

    integer(kind=int1) :: SUCCESS = 0_int1
    integer(kind=int1) :: FAILURE = 1_int1
    integer(kind=int1) :: INFORMATION = 2_int1
    integer(kind=int1) :: WARNING = 3_int1
    integer(kind=int1) :: EOF = 4_int1
    integer(kind=int1) :: UNDEFINED = 5_int1
    integer(kind=int1) :: EXISTS =  6_int1
!   integer(kind=int1) :: EXIT = 7_int1

    integer(kind=int1) :: SNOW_NOT_AVAILABLE = 1_int1
    integer(kind=int1) :: NWP_SNOW = 2_int1
    integer(kind=int1) :: IMS_SNOW = 3_int1

    integer(kind=int1) :: CONSTANT_SFC_EMISS = 1_int1
    integer(kind=int1) :: TABLE_SFC_EMISS    = 2_int1
    integer(kind=int1) :: SEEBOR_SFC_EMISS   =3_int1

    integer(kind=int1) :: CONSTANT_SFC_ALB   = 1_int1
    integer(kind=int1) :: TABLE_SFC_ALB      = 2_int1
    integer(kind=int1) :: MODIS_SFC_ALB      = 3_int1

    integer(kind=int4) :: LITTLE_ENDIAN      = 0_int1
    integer(kind=int4) :: BIG_ENDIAN         = 1_int1
    integer(kind=int4) :: SIGNED             = 1_int1
    integer(kind=int4) :: UNSIGNED           = 0_int1
    integer(kind=int4) :: SWAP               = 1_int1
    integer(kind=int4) :: NOSWAP             = 0_int1

    !--- scaling options
    integer(kind=int1) :: NO_SCALING = 0_int1
    integer(kind=int1) :: LINEAR_SCALING = 1_int1
    integer(kind=int1) :: LOG10_SCALING = 2_int1
    integer(kind=int1) :: SQUARE_ROOT_SCALING = 3_int1

    !--- ECM DQF - Don't know exact enumeration for SAPF. Will need to fix
    integer(kind=int1) :: ECM_GOOD_QF = 0_int1
    integer(kind=int1) :: ECM_BAD_QF = 1_int1
    integer(kind=int1) :: ECM_DEGRADED_QF = 2_int1
    integer(kind=int1) :: ECM_SPACE_QF = 3_int1
    integer(kind=int1) :: ECM_FILL_QF = -127_int1

    !General product quality flags - VOLCAT
    integer(kind=int4) :: INVALID_PRODUCT = 0_int1
    integer(kind=int4) :: UNKNOWN_PRODUCT_QUALITY = 1_int1
    integer(kind=int4) :: LOW_PRODUCT_QUALITY = 2_int1
    integer(kind=int4) :: HIGH_PRODUCT_QUALITY_NOintERP = 3_int1
    integer(kind=int4) :: HIGH_PRODUCT_QUALITY_intERP = 4_int1
    integer(kind=int4) :: INVALID_FILLED = 5_int1
    integer(kind=int4) :: INVALID_REPLACED = 6_int1

  END TYPE symbol_struct

!ccm  character(len=4), parameter:: SOLAR_OBS_TYPE = "SOLAR"
  character(len=5), public, parameter:: SOLAR_OBS_TYPE = "SOLAR"
  character(len=5),  public, parameter:: LUNAR_OBS_TYPE = "LUNAR"
  character(len=5),  public, parameter:: MIXED_OBS_TYPE = "MIXED"
  character(len=5),  public, parameter:: THERMAL_OBS_TYPE = "THERM"
     

  TYPE(symbol_struct), public, save :: sym
  
end module CONSTANTS_MOD
