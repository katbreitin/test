! $Id: cx_hdf4_mod.f90 3840 2020-05-13 17:39:08Z heidinger $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: hdf_params.f90 (src)
!       HDF_PARAMS (program)
!
! PURPOSE: This module contains routines used to read and write to the hdf
!          output files from CLAVR-X
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!
! COPYRIGHT
! (c) This code is copyrighted by the author and all NOAA restrictions apply
!
! Dependencies:  (The following are names of other CLAVR-x modules)
!  CONSTANTS_MOD
!  HDF_MOD
!  SCALING_PARAMETERS_MOD
!
! Calling Sequence:
!  use HDF_PARAMS
!
! Public Routines within this module
!  SCALE_VECTOR_I1_RANK1
!  SCALE_VECTOR_I1_RANK2
!  SCALE_VECTOR_I1_RANK3
!  SCALE_VECTOR_I2_RANK1
!  SCALE_VECTOR_I2_RANK2
!  SCALE_VECTOR_I2_RANK3
!  UNSCALE_VECTOR_I1_RANK1
!  WRITE_CLAVRX_HDF4_SDS
!  HDF_TSTAMP
!  WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES
!  READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES
!  READ_CLAVRX_HDF4_SDS_RANK1
!
!--------------------------------------------------------------------------------------
module CX_HDF4_MOD

use CONSTANTS_MOD, only: &
   int1 &
   , int2 &
   , int4 &
   , real4 &
   , real8 &
   , MISSING_VALUE_INT1 &
   , MISSING_VALUE_INT2 &
   , MISSING_VALUE_INT4 &
   , MISSING_VALUE_REAL4 &
   , MISSING_VALUE_REAL8 &
   , ONE_BYTE_MIN &
   , ONE_BYTE_MAX &
   , TWO_BYTE_MIN &
   , TWO_BYTE_MAX

use univ_fp_comparison_mod, only: operator(.EQfp.)
   
use HDF, only: &
   DFNT_INT16 &
   , DFNT_CHAR8 &
   , DFNT_FLOAT32 &
   , DFNT_FLOAT64 &
   , DFNT_INT32 &
   , DFNT_INT8 &
   , DFACC_READ &
   , FAIL &
   , SUCCEED &
   , MAX_RANK_HDF &
   , DFNT_CHAR &
   , DFNT_FLOAT64 &
   , DFNT_UINT8 &
   , DFNT_UINT16 &
   , DFNT_UINT32

use LEVEL2_STRUCTURES_MOD, only: Sds_Struct, L2_Glob_Attr_Definition, Clavrx_Global_Attr

implicit none

private

 private:: CHECK_EDGE_VS_OUTPUT_SHAPE

 public:: SCALE_VECTOR_I1_RANK1, &
          SCALE_VECTOR_I1_RANK2, &
          SCALE_VECTOR_I1_RANK3, &
          SCALE_VECTOR_I2_RANK1, &
          SCALE_VECTOR_I2_RANK2, &
          SCALE_VECTOR_I2_RANK3, &
          UNSCALE_VECTOR_I1_RANK1, &
          WRITE_CLAVRX_HDF4_SDS, &
          WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES,  &
          READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES, &
          READ_CLAVRX_HDF4_SDS_RANK1

 public:: OPEN_FILE_HDF_READ
 public:: HDF_SDS_DIMENSIONS_READER
 public:: HDF_SDS_READER
 public:: HDF_SDS_ATTRIBUTE_READER
 public:: CLOSE_FILE_HDF_READ

 public:: WRITE_CLAVRX_HDF_SDS_1D
 public:: WRITE_CLAVRX_HDF_SDS_2D
 public:: WRITE_CLAVRX_HDF_SDS_3D
 public:: READ_CLAVRX_HDF_SDS_3D
 public:: READ_CLAVRX_HDF_SDS_2D
 public:: READ_CLAVRX_HDF_SDS_1D
 public:: READ_HDF_GLOBAL_ATTRIBUTE_NUM
 public:: READ_HDF_GLOBAL_ATTRIBUTE_STR

  interface HDF_SDS_READER
    module procedure READ_HDF_SDS_INT8_1D, &
                     READ_HDF_SDS_INT16_1D, &
                     READ_HDF_SDS_INT32_1D, &
                     READ_HDF_SDS_FLOAT32_1D, &
                     READ_HDF_SDS_FLOAT64_1D, &
                     READ_HDF_SDS_INT8_2D, &
                     READ_HDF_SDS_INT16_2D, &
                     READ_HDF_SDS_INT32_2D, &
                     READ_HDF_SDS_FLOAT32_2D, &
                     READ_HDF_SDS_FLOAT64_2D, &
                     READ_HDF_SDS_INT8_3D, &
                     READ_HDF_SDS_INT16_3D, &
                     READ_HDF_SDS_INT32_3D, &
                     READ_HDF_SDS_FLOAT32_3D, &
                     READ_HDF_SDS_FLOAT64_3D
  end interface

  interface HDF_SDS_ATTRIBUTE_READER
    module procedure READ_HDF_ATTRIBUTE_CHAR8_SCALAR,   &
                     READ_HDF_ATTRIBUTE_INT8_SCALAR,    &
                     READ_HDF_ATTRIBUTE_INT16_SCALAR,   &
                     READ_HDF_ATTRIBUTE_INT32_SCALAR,   &
                     READ_HDF_ATTRIBUTE_FLOAT32_SCALAR, &
                     READ_HDF_ATTRIBUTE_FLOAT64_SCALAR, &
                     READ_HDF_ATTRIBUTE_INT8_VECTOR,    &
                     READ_HDF_ATTRIBUTE_INT16_VECTOR,   &
                     READ_HDF_ATTRIBUTE_INT32_VECTOR,   &
                     READ_HDF_ATTRIBUTE_FLOAT32_VECTOR, &
                     READ_HDF_ATTRIBUTE_FLOAT64_VECTOR
  end interface


 interface WRITE_CLAVRX_HDF4_SDS 
     module procedure  &
         WRITE_CLAVRX_HDF4_SDS_RANK1,  &
         WRITE_CLAVRX_HDF4_SDS_RANK2,  &
         WRITE_CLAVRX_HDF4_SDS_RANK3
 end interface

 interface READ_CLAVRX_HDF_SDS
     module procedure  &
         READ_CLAVRX_HDF_SDS_1D,  &
         READ_CLAVRX_HDF_SDS_2D,  &
         READ_CLAVRX_HDF_SDS_3D
 end interface

 interface WRITE_CLAVRX_HDF_SDS
      module procedure  &
         WRITE_CLAVRX_HDF_SDS_1D,  &
         WRITE_CLAVRX_HDF_SDS_2D,  &
         WRITE_CLAVRX_HDF_SDS_3D
 end interface

 !---- begin module variable definition
 character(len=256), save, public:: renav_data_from
 character(len=11), parameter :: EXE_PROMPT = 'CX_HDF4_MOD'


 !--- scaling options
 integer(kind=int1), parameter, public:: NO_SCALING = 0_int1
 integer(kind=int1), parameter, public:: LINEAR_SCALING = 1_int1
 integer(kind=int1), parameter, public:: LOG10_SCALING = 2_int1
 integer(kind=int1), parameter, public:: SQUARE_ROOT_SCALING = 3_int1

 real, parameter, private:: DEFAULT_MISSING_VALUE= MISSING_VALUE_real4

 ! HDF function declarations
 integer(kind=int4), external:: sfcreate, sfdimid, sfsdmname, sfwdata, sfendacc, &
                                sfscatt, sfsnatt, sfschnk, sfselect, sfn2index,   &
                                sfrnatt,sfrcatt,sffattr,sfrdata,sfginfo, sfstart, &
                                sfend

contains

!-------------------------------------------------------------------------
! routine to write global attributes to clavrx files
!
!-------------------------------------------------------------------------
 subroutine WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES(hdf_file_id)

 integer(kind=int4), intent(in):: hdf_file_id

 integer:: sfscatt, sfsnatt, hglibver
 integer:: major_v, minor_v, release
 !integer(kind=int2):: Istatus = 0
 integer(kind=int4):: Istatus = 0
 !character(len=80) :: hdf_ver
 !character(len=36) :: machine


 character(len=6), parameter:: Conventions_String = "CF-1.6"
 character(len=38), parameter:: Metadata_Conventions_String = "CF-1.6, Unidata Dataset Discovery v1.0"
 character(len=42), parameter:: Standard_Name_Vocabulary_String = "CF Standard Name Table (v25, 05 July 2013)"
 character(len=13), parameter:: Naming_Authority_String = "gov.noaa.ncdc"
 character(len=37), parameter:: License_String = "No constraints on data access or use."


! Definition of strings used as HDF attributes.
!
! version history 4.1 - delivered to OSDPD in November 2006
! version history 4.2 - demonstrated on METOP
! version history 4.3 - included surface emissivity fields
! version history 4.4 - included lrc
! version history 5.0 - included modis white sky and ash protoype 
! version history 5.1 - rtm structures now 101 levels and 
!                       reorganized level-1b ingest to 
!                       read segment all at once prior to processing
!                       first version with working volcanic ash
! version history 5.2    bayesian cloud mask and DCOMP
! version history 6.0  - MODIS capability begin
! version history 6.5  - VIIRS capability begin

 character(len=36), parameter :: creator0 = "CLAVR-x + PATMOS-x "

 character(len=100):: Title_String
 character(len=100):: Calibration_String
 character(len=100):: Product_Version_String
 character(len=100):: Status_String
 character(len=100):: Institution_String
 character(len=100):: Program_String
 character(len=500):: Summary_String 
 character(len=200):: Variable_String
 character(len=500):: Keywords_String
 character(len=200):: Keywords_Vocabulary_String
 character(len=100):: Time_Coverage_Resolution_String
 character(len=100):: Metadata_Link_String
 character(len=100):: Spatial_Resolution_String

 include 'default_xDF_version_info.inc'

 !complete the creator string with the version number
 Clavrx_Global_Attr%creator = trim(creator0)//trim(Product_Version_String)


 call getenv ("HOST",Clavrx_Global_Attr%machine)
 if (len_trim(Clavrx_Global_Attr%machine) == 0) call getenv ("HOSTNAME",Clavrx_Global_Attr%machine)
 if (len_trim(Clavrx_Global_Attr%machine) == 0) Clavrx_Global_Attr%machine = "unknown"

 Clavrx_Global_Attr%plang = "F90"
 Clavrx_Global_Attr%timestamp = trim(hdf_timestamp())

!--- determine HDF library version
if (hglibver(major_v, minor_v, release, Clavrx_Global_Attr%hdf_ver) /= SUCCEED) then
   print *, "could not determine HDF library version"
   stop 961
end if

!---- describe CLAVR-x
Istatus = sfscatt(hdf_file_id, "HDF_LIB_VERSION", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%hdf_ver), trim(Clavrx_Global_Attr%hdf_ver))+Istatus
Istatus = sfscatt(hdf_file_id, "MACHINE", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%machine), trim(Clavrx_Global_Attr%machine))+Istatus
Istatus = sfscatt(hdf_file_id, "PROGLANG", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%plang), trim(Clavrx_Global_Attr%plang))+Istatus

!--- CF compliant global attributes required for NCDC delivery
Istatus = sfscatt(hdf_file_id, "date_created", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%timestamp), trim(Clavrx_Global_Attr%timestamp))+Istatus
Istatus = sfscatt(hdf_file_id, "product_version", DFNT_CHAR8,  &
                                len_trim(Product_Version_String), trim(Product_Version_String))+Istatus
Istatus = sfscatt(hdf_file_id,"summary", DFNT_CHAR8,len_trim(Summary_String),trim(Summary_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"cdr_variable", DFNT_CHAR8,len_trim(Variable_String),trim(Variable_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"institution", DFNT_CHAR8,len_trim(Institution_String),trim(Institution_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"cdr_program", DFNT_CHAR8, len_trim(Program_String),trim(Program_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"title", DFNT_CHAR8,len_trim(Title_String),trim(Title_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"calibration_version", DFNT_CHAR8,len_trim(Calibration_String),trim(Calibration_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"keywords", DFNT_CHAR8,len_trim(Keywords_String),trim(Keywords_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"keywords_vocabulary", DFNT_CHAR8,len_trim(Keywords_Vocabulary_String),trim(Keywords_Vocabulary_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"time_coverage_resolution", DFNT_CHAR8,len_trim(Time_Coverage_Resolution_String),trim(Time_Coverage_Resolution_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"metadata_link", DFNT_CHAR8,len_trim(Metadata_Link_String),trim(Metadata_Link_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"spatial_resolution", DFNT_CHAR8,len_trim(Spatial_Resolution_String),trim(Spatial_Resolution_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"Conventions", DFNT_CHAR8,len_trim(Conventions_String),trim(Conventions_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"title", DFNT_CHAR8,len_trim(Title_String),Title_String) + Istatus
Istatus = sfscatt(hdf_file_id,"Metadata_Conventions", DFNT_CHAR8,  &
                   len_trim(Metadata_Conventions_String),trim(Metadata_Conventions_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"standard_name_vocabulary", DFNT_CHAR8,  &
                   len_trim(Standard_Name_Vocabulary_String),trim(Standard_Name_Vocabulary_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"naming_authority", DFNT_CHAR8, &
                   len_trim(Naming_Authority_String),trim(Naming_Authority_String)) + Istatus
Istatus = sfscatt(hdf_file_id,"license", DFNT_CHAR8, len_trim(License_String),trim(License_String)) + Istatus

!--- describe the data
Istatus = sfscatt(hdf_file_id, "sensor", DFNT_CHAR8,len_trim(Clavrx_Global_Attr%sensor_name),trim(Clavrx_Global_Attr%sensor_name))+Istatus
Istatus = sfscatt(hdf_file_id, "platform", DFNT_CHAR8,len_trim(Clavrx_Global_Attr%platform_name),trim(Clavrx_Global_Attr%platform_name))+Istatus

Istatus = sfscatt(hdf_file_id, "FILENAME", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%file_name), trim(Clavrx_Global_Attr%file_name))+Istatus
Istatus = sfscatt(hdf_file_id, "L1B", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%file_1b), trim(Clavrx_Global_Attr%file_1b))+Istatus
Istatus = sfsnatt(hdf_file_id, "RESOLUTION_KM", DFNT_FLOAT32, 1, Clavrx_Global_Attr%resolution_km)+Istatus
Istatus = sfsnatt(hdf_file_id, "START_YEAR", DFNT_INT16, 1, Clavrx_Global_Attr%start_year)+Istatus
Istatus = sfsnatt(hdf_file_id, "START_DAY", DFNT_INT16, 1, Clavrx_Global_Attr%start_day)+Istatus
Istatus = sfsnatt(hdf_file_id, "START_TIME", DFNT_FLOAT32, 1, Clavrx_Global_Attr%start_time)+Istatus
Istatus = sfsnatt(hdf_file_id, "END_YEAR", DFNT_INT16, 1, Clavrx_Global_Attr%end_year)+Istatus
Istatus = sfsnatt(hdf_file_id, "END_DAY", DFNT_INT16, 1, Clavrx_Global_Attr%end_day)+Istatus
Istatus = sfsnatt(hdf_file_id, "END_TIME", DFNT_FLOAT32, 1, Clavrx_Global_Attr%end_time)+Istatus
Istatus = sfscatt(hdf_file_id, "CLOUD_MASK_MODE", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%Mask_mode), trim(Clavrx_Global_Attr%Mask_mode))+Istatus
Istatus = sfscatt(hdf_file_id, "ACHA_MODE_FINAL", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%acha_mode), trim(Clavrx_Global_Attr%acha_mode))+Istatus
Istatus = sfscatt(hdf_file_id, "ACHA_MODE_USER", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%acha_mode_user), trim(Clavrx_Global_Attr%acha_mode_user))+Istatus
Istatus = sfsnatt(hdf_file_id, "DCOMP_MODE", DFNT_INT32, 1, Clavrx_Global_Attr%dcomp_mode)+Istatus
Istatus = sfsnatt(hdf_file_id, "WMO_SATELLITE_CODE", DFNT_INT32, 1, Clavrx_Global_Attr%wmo_sc_code)+Istatus
Istatus = sfscatt(hdf_file_id, "REFL_0_65UM_NOM_DARK_COMPOSITE_NAME", DFNT_CHAR8, &
                               len_trim(Clavrx_Global_Attr%dark_name),trim(Clavrx_Global_Attr%dark_name))+Istatus
Istatus = sfscatt(hdf_file_id, "NAIVE_BAYESIAN_CLOUD_MASK_NAME", DFNT_CHAR8, &
                               len_trim(Clavrx_Global_Attr%mask_name),trim(Clavrx_Global_Attr%mask_name))+Istatus

Istatus = sfsnatt(hdf_file_id, "GEO_SUB_LON", DFNT_FLOAT32, 1, Clavrx_Global_Attr%geo_sub_lon)+Istatus
Istatus = sfsnatt(hdf_file_id, "GEO_SUB_LAT", DFNT_FLOAT32, 1, Clavrx_Global_Attr%geo_sub_lat)+Istatus
!--- describe subsetting
Istatus = sfsnatt(hdf_file_id, "SUBSET_FLAG", DFNT_INT32, 1, Clavrx_Global_Attr%subset_flag)+Istatus
Istatus = sfsnatt(hdf_file_id, "LAT_SOUTH_SUBSET", DFNT_FLOAT32, 1, Clavrx_Global_Attr%lat_south_subset)+Istatus
Istatus = sfsnatt(hdf_file_id, "LAT_NORTH_SUBSET", DFNT_FLOAT32, 1, Clavrx_Global_Attr%lat_north_subset)+Istatus
Istatus = sfsnatt(hdf_file_id, "LON_WEST_SUBSET", DFNT_FLOAT32, 1, Clavrx_Global_Attr%lon_west_subset)+Istatus
Istatus = sfsnatt(hdf_file_id, "LON_EAST_SUBSET", DFNT_FLOAT32, 1, Clavrx_Global_Attr%lon_east_subset)+Istatus
Istatus = sfscatt(hdf_file_id, "SUBSET_NAME", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%subset_name), trim(Clavrx_Global_Attr%subset_name))+Istatus

!--- data type
Istatus = sfscatt(hdf_file_id, "DATA_TYPE", DFNT_CHAR8,len_trim(Clavrx_Global_Attr%data_type),trim(Clavrx_Global_Attr%data_type))+Istatus


!--- NCDC Attributes

!-- other global attributes for the level3 file.
     if (Clavrx_Global_Attr%data_type(1:4) == "GRID") then
      Istatus = sfsnatt(hdf_file_id, "NUM_CELLS_TOTAL", DFNT_INT32,1,Clavrx_Global_Attr%num_cells)+Istatus
      Istatus = sfsnatt(hdf_file_id, "NUM_CELLS_WITH_DATA", DFNT_INT32,1,Clavrx_Global_Attr%num_cells_with_data)+Istatus
      Istatus = sfsnatt(hdf_file_id, "GRIDCELL_RESOLUTION", DFNT_FLOAT32,1,Clavrx_Global_Attr%dlat)+Istatus
      Istatus = sfscatt(hdf_file_id, "GRIDCELL_RESOLUTION_UNIT", DFNT_CHAR8,6,"degree")+Istatus
      if (Clavrx_Global_Attr%grid_format(1:10) == "EQUAL_AREA") then
       Istatus = sfscatt(hdf_file_id, "GRIDCELL_FORMAT", DFNT_CHAR8,10,"EQUAL_AREA")+Istatus
      else
       Istatus = sfscatt(hdf_file_id, "GRIDCELL_FORMAT", DFNT_CHAR8,11,"EQUAL_ANGLE")+Istatus
      endif
    endif

!---- processing flags
 Istatus = sfscatt(hdf_file_id, "USE_1B_THERMAL_CALIBRATION_FLAG", DFNT_INT32,1,Clavrx_Global_Attr%therm_cal_1b)+Istatus
 Istatus = sfscatt(hdf_file_id, "USE_1B_REFLECTANCE_CALIBRATION_FLAG", DFNT_INT32,1,Clavrx_Global_Attr%Ref_cal_1b)+Istatus
! Istatus = sfscatt(hdf_file_id, "RENAVIGATION_DATA_FROM", DFNT_CHAR8, &
!     len_trim(renav_data_from), trim(renav_data_from))+Istatus
 Istatus = sfscatt(hdf_file_id, "RENAVIGATION_FLAG", DFNT_INT32,1,Clavrx_Global_Attr%nav_opt)+Istatus
 Istatus = sfscatt(hdf_file_id, "USE_SST_ANALYSIS_FLAG", DFNT_INT32,1,Clavrx_Global_Attr%use_sst_anal)+Istatus
 Istatus = sfsnatt(hdf_file_id, "NWP_OPT", DFNT_INT32,1,Clavrx_Global_Attr%nwp_opt)+Istatus
 Istatus = sfsnatt(hdf_file_id, "MODIS_CLEAR_SKY_REFLECTANCE_FLAG", DFNT_INT32,1,Clavrx_Global_Attr%modis_clr_alb_flag)+Istatus

!-- reflectance channel calibration
Istatus = sfsnatt(hdf_file_id, "CH1_GAIN_LOW", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch1_gain_low)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH1_GAIN_HIGH", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch1_gain_high)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH1_SWITCH_COUNT", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch1_switch_count)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH1_DARK_COUNT", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch1_dark_count)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH2_GAIN_LOW", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch2_gain_low)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH2_GAIN_HIGH", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch2_gain_high)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH2_SWITCH_COUNT", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch2_switch_count)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH2_DARK_COUNT", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch2_dark_count)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH3A_GAIN_LOW", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch3a_gain_low)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH3A_GAIN_HIGH", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch3a_gain_high)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH3A_SWITCH_COUNT", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch3a_switch_count)+Istatus
Istatus = sfsnatt(hdf_file_id, "CH3A_DARK_COUNT", DFNT_FLOAT32,1,Clavrx_Global_Attr%ch3a_dark_count)+Istatus
Istatus = sfsnatt(hdf_file_id, "SUN_EARTH_DISTANCE", DFNT_FLOAT32,1,Clavrx_Global_Attr%sun_earth_distance)+Istatus

!--- thermal calibration constants
Istatus = sfsnatt(hdf_file_id, "C1", DFNT_FLOAT32,1,Clavrx_Global_Attr%c1)+Istatus
Istatus = sfsnatt(hdf_file_id, "C2", DFNT_FLOAT32,1,Clavrx_Global_Attr%c2)+Istatus
Istatus = sfsnatt(hdf_file_id, "A_20", DFNT_FLOAT32,1,Clavrx_Global_Attr%a_20)+Istatus
Istatus = sfsnatt(hdf_file_id, "B_20", DFNT_FLOAT32,1,Clavrx_Global_Attr%b_20)+Istatus
Istatus = sfsnatt(hdf_file_id, "NU_20", DFNT_FLOAT32,1,Clavrx_Global_Attr%nu_20)+Istatus
Istatus = sfsnatt(hdf_file_id, "A_31", DFNT_FLOAT32,1,Clavrx_Global_Attr%a_31)+Istatus
Istatus = sfsnatt(hdf_file_id, "B_31", DFNT_FLOAT32,1,Clavrx_Global_Attr%b_31)+Istatus
Istatus = sfsnatt(hdf_file_id, "NU_31", DFNT_FLOAT32,1,Clavrx_Global_Attr%nu_31)+Istatus
Istatus = sfsnatt(hdf_file_id, "A_32", DFNT_FLOAT32,1,Clavrx_Global_Attr%a_32)+Istatus
Istatus = sfsnatt(hdf_file_id, "B_32", DFNT_FLOAT32,1,Clavrx_Global_Attr%b_32)+Istatus
Istatus = sfsnatt(hdf_file_id, "NU_32", DFNT_FLOAT32,1,Clavrx_Global_Attr%nu_32)+Istatus
Istatus = sfsnatt(hdf_file_id, "SOLAR_20_NU", DFNT_FLOAT32,1,Clavrx_Global_Attr%solar_Ch20_nu)+Istatus
Istatus = sfsnatt(hdf_file_id, "LWIR_FOCAL_PLANE_TEMPERATURE", DFNT_FLOAT32,1,Clavrx_Global_Attr%LWIR_Focal_Plane_Temperature)+Istatus

Istatus = sfsnatt(hdf_file_id, "TIME_ERROR_SECONDS", DFNT_FLOAT32,1,Clavrx_Global_Attr%timerr_seconds)+Istatus

if (Istatus /= 0) then
  print *, "error writing run-control flags as global HDF attributes"
  stop 963
endif

contains

  function hdf_timestamp() result(string)
     character (len = 36) ::string
     character(len=10), dimension(3):: ctime

     call date_and_time(ctime(1), ctime(2), ctime(3))

     ! Timestamp string format in accordance with the ISO-8601 standard.
     string = ctime(1)(1:4)//"-"//ctime(1)(5:6)//"-"//ctime(1)(7:8)&
                //"T"//ctime(2)(1:2)//":"//ctime(2)(3:4)//":"//ctime(2)(5:6)&
                //ctime(3)(1:3)//":"//ctime(3)(4:5)
     return
   end function hdf_timestamp

 end subroutine WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES

!------------------------------------------------------------------------------------------
! READ GLOBAL ATTRIBUTES FROM A CLAVR-x HDF FILE
!------------------------------------------------------------------------------------------
 subroutine READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES(hdf_file_id)

 integer(kind=int4), intent(in):: hdf_file_id

 !integer(kind=int2):: Istatus = 0
 integer(kind=int4):: Istatus = 0
 integer(kind=int4):: blank_int4
 real(kind=real4):: blank_real4

 !--- hdf  calls
 integer:: sfrcatt, sfrnatt, sffattr

 blank_int4 = 0
 blank_real4 = 0.0

 !--- attributes about the code used to make this data
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"PROCESSOR"), Clavrx_Global_Attr%creator) + Istatus
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"CREATED"), Clavrx_Global_Attr%timestamp) + Istatus
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"HDF_LIB_VERSION"), Clavrx_Global_Attr%hdf_ver(1:48)) + Istatus
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"PROGLANG"), Clavrx_Global_Attr%plang) + Istatus

 !--- file names
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"FILENAME"), Clavrx_Global_Attr%file_name) + Istatus
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"L1B"), Clavrx_Global_Attr%file_1b) + Istatus

 !--- temporal attributes
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "START_YEAR"), Clavrx_Global_Attr%start_year) + Istatus
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "START_DAY"), Clavrx_Global_Attr%start_day) + Istatus
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "START_TIME"), Clavrx_Global_Attr%start_time) + Istatus
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "END_YEAR"), Clavrx_Global_Attr%end_year) + Istatus
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "END_DAY"), Clavrx_Global_Attr%end_day) + Istatus
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "END_TIME"), Clavrx_Global_Attr%end_time) + Istatus
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "WMO_SATELLITE_CODE"),Clavrx_Global_Attr%wmo_sc_code) + Istatus
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id, "platform"),Clavrx_Global_Attr%platform_name) + Istatus           
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id, "sensor"),Clavrx_Global_Attr%sensor_name) + Istatus          

 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id, "REFL_0_65UM_NOM_DARK_COMPOSITE_NAME"), &
                                                     Clavrx_Global_Attr%dark_name) + Istatus
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id, "NAIVE_BAYESIAN_CLOUD_MASK_NAME"), &
                                                     Clavrx_Global_Attr%mask_name) + Istatus

!--- describe subsetting
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id,"SUBSET_FLAG"), Clavrx_Global_Attr%subset_flag)+Istatus
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id,"LAT_SOUTH_SUBSET"), Clavrx_Global_Attr%lat_south_subset)+Istatus
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id,"LAT_NORTH_SUBSET"), Clavrx_Global_Attr%lat_north_subset)+Istatus
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id,"LON_WEST_SUBSET"), Clavrx_Global_Attr%lon_west_subset)+Istatus
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id,"LON_EAST_SUBSET"), Clavrx_Global_Attr%lon_east_subset)+Istatus
Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"SUBSET_NAME"), Clavrx_Global_Attr%subset_name)+Istatus

 !--- algorithm options
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "ACHA_MODE"), Clavrx_Global_Attr%acha_mode) + Istatus
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "DCOMP_MODE"), Clavrx_Global_Attr%dcomp_mode) + Istatus

!--- data type
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"DATA_TYPE"), Clavrx_Global_Attr%data_type) + Istatus

!--- grid variables
if (Clavrx_Global_Attr%data_type(1:4) == "GRID") then
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NUM_CELLS_WITH_DATA"), Clavrx_Global_Attr%num_cells_with_data) + Istatus
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NUM_CELLS_TOTAL"), Clavrx_Global_Attr%num_cells) + Istatus
 Istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id, "GRIDCELL_FORMAT"), Clavrx_Global_Attr%grid_format) + Istatus
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "GRIDCELL_RESOLUTION"), Clavrx_Global_Attr%grid_resolution) + Istatus
 Clavrx_Global_Attr%resolution_km = -999.0
else
 Clavrx_Global_Attr%num_cells_with_data = blank_int4
 Clavrx_Global_Attr%num_cells = blank_int4
 Clavrx_Global_Attr%grid_format = "     " 
 Clavrx_Global_Attr%grid_resolution = blank_real4
 Clavrx_Global_Attr%resolution_km = -999.0 !if RESOLUTION_KM global attribute is not present, missing value
 Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "RESOLUTION_KM"), Clavrx_Global_Attr%resolution_km) + Istatus
endif

!--- flags
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "USE_1B_THERMAL_CALIBRATION_FLAG"), Clavrx_Global_Attr%therm_cal_1b)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "USE_1B_REFLECTANCE_CALIBRATION_FLAG"), Clavrx_Global_Attr%Ref_cal_1b)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "USE_SST_ANALYSIS_FLAG"), Clavrx_Global_Attr%use_sst_anal)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NWP_FLAG"), Clavrx_Global_Attr%nwp_opt)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "MODIS_CLEAR_SKY_REFLECTANCE_FLAG"), Clavrx_Global_Attr%modis_clr_alb_flag)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "RENAVIGATION_FLAG"), Clavrx_Global_Attr%nav_opt)

!--- calibration attributes
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "C1"), Clavrx_Global_Attr%c1)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "C2"), Clavrx_Global_Attr%c2)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "A_20"), Clavrx_Global_Attr%a_20)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "B_20"), Clavrx_Global_Attr%b_20)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NU_20"), Clavrx_Global_Attr%nu_20)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "A_31"), Clavrx_Global_Attr%a_31)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "B_31"), Clavrx_Global_Attr%b_31)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NU_31"), Clavrx_Global_Attr%nu_31)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "A_32"), Clavrx_Global_Attr%a_32)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "B_32"), Clavrx_Global_Attr%b_32)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NU_32"), Clavrx_Global_Attr%nu_32)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "SOLAR_20_NU"), Clavrx_Global_Attr%solar_Ch20_nu)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "TIME_ERROR_SECONDS"), Clavrx_Global_Attr%timerr_seconds)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "SUN_EARTH_DISTANCE"), Clavrx_Global_Attr%sun_earth_distance)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH1_GAIN_LOW"), Clavrx_Global_Attr%ch1_gain_low)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH1_GAIN_HIGH"), Clavrx_Global_Attr%ch1_gain_high)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH1_SWITCH_COUNT"), Clavrx_Global_Attr%ch1_switch_count)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH1_DARK_COUNT"), Clavrx_Global_Attr%ch1_dark_count)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH2_GAIN_LOW"), Clavrx_Global_Attr%ch2_gain_low)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH2_GAIN_HIGH"), Clavrx_Global_Attr%ch2_gain_high)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH2_SWITCH_COUNT"), Clavrx_Global_Attr%ch2_switch_count)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH2_DARK_COUNT"), Clavrx_Global_Attr%ch2_dark_count)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH6_GAIN_LOW"), Clavrx_Global_Attr%ch3a_gain_low)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH6_GAIN_HIGH"), Clavrx_Global_Attr%ch3a_gain_high)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH6_SWITCH_COUNT"), Clavrx_Global_Attr%ch3a_switch_count)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH6_DARK_COUNT"), Clavrx_Global_Attr%ch3a_dark_count)
Istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "LWIR_FOCAL_PLANE_TEMPERATURE"), Clavrx_Global_Attr%LWIR_Focal_Plane_Temperature)

end subroutine READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES

!-----------------------------------------------------------------------------------------------------
! VECTOR SCALING ROUTINES
!
! iscaled = 0 = no scaling
!           1 = linear
!           2 = log10
!           3 = sqrt
!-----------------------------------------------------------------------------------------------------
 subroutine SCALE_VECTOR_I1_RANK1(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i1)
   real, dimension(:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int1), dimension(:),  intent(out):: temp_i1
   real, dimension(size(temp_r4,1)):: scratch_r4

   scratch_r4 = 0.0

   !---- linear
   if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = int(real(ONE_BYTE_MIN) + scratch_r4 * real(ONE_BYTE_MAX - ONE_BYTE_MIN),kind=int1)
   endif
                                                                                                                                                          
   !---- log10
   if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i1 = int(real(ONE_BYTE_MIN) + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * real(ONE_BYTE_MAX - ONE_BYTE_MIN),kind=int1)
    endif
    !---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = int(real(ONE_BYTE_MIN) + sqrt(scratch_r4) * real(ONE_BYTE_MAX - ONE_BYTE_MIN),kind=int1)
    endif
                                                                                                                                                          
    !--- set scaled missing values
    where (temp_r4 .EQfp. unscaled_missing)
         temp_i1 = MISSING_VALUE_INT1
    endwhere

 end subroutine SCALE_VECTOR_I1_RANK1

 subroutine SCALE_VECTOR_I1_RANK2(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i1)
   real, dimension(:,:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int1), dimension(:,:),  intent(out):: temp_i1
   real, dimension(size(temp_r4,1),size(temp_r4,2)):: scratch_r4

   scratch_r4 = 0.0
   !---- linear
   if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = int(real(ONE_BYTE_MIN) + scratch_r4 * real(ONE_BYTE_MAX - ONE_BYTE_MIN),kind=int1)
   endif

   !---- log10
   if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i1 = int(real(ONE_BYTE_MIN) + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * real(ONE_BYTE_MAX - ONE_BYTE_MIN),kind=int1)
   endif

   !---- square root
   if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = int(real(ONE_BYTE_MIN) + sqrt(scratch_r4) * real(ONE_BYTE_MAX - ONE_BYTE_MIN),kind=int1)
   endif

   !--- set scaled missing values
   where (temp_r4 .EQfp. unscaled_missing)
         temp_i1 = MISSING_VALUE_INT1
   endwhere

 end subroutine SCALE_VECTOR_I1_RANK2

 subroutine SCALE_VECTOR_I1_RANK3(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i1)
   real, dimension(:,:,:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int1), dimension(:,:,:),  intent(out):: temp_i1
   real, dimension(size(temp_r4,1),size(temp_r4,2),size(temp_r4,3)):: scratch_r4

   scratch_r4 = 0.0
   !---- linear
   if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = int(real(ONE_BYTE_MIN) + scratch_r4 * real(ONE_BYTE_MAX - ONE_BYTE_MIN),kind=int1)
   endif
   !---- log10
   if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i1 = int(real(ONE_BYTE_MIN) + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * real(ONE_BYTE_MAX - ONE_BYTE_MIN),kind=int1)
    endif
    !---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = int(real(ONE_BYTE_MIN) + sqrt(scratch_r4) * real(ONE_BYTE_MAX - ONE_BYTE_MIN),kind=int1)
    endif
    !--- set scaled missing values
    where (temp_r4 .EQfp. unscaled_missing)
         temp_i1 = MISSING_VALUE_INT1
    endwhere
 end subroutine SCALE_VECTOR_I1_RANK3


 subroutine SCALE_VECTOR_I2_RANK1(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i2)
   real, dimension(:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int2), dimension(:), intent(out):: temp_i2
   real, dimension(size(temp_r4,1)):: scratch_r4

   !---- linear
   if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = int(real(TWO_BYTE_MIN) + scratch_r4 * real(TWO_BYTE_MAX - TWO_BYTE_MIN),kind=int2)
   endif
   !---- log10
   if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i2 = int(real(TWO_BYTE_MIN) + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * real(TWO_BYTE_MAX - TWO_BYTE_MIN),kind=int2)
   endif
   !---- square root
   if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = int(real(TWO_BYTE_MIN) + sqrt(scratch_r4) * real(TWO_BYTE_MAX - TWO_BYTE_MIN),kind=int2)
   endif
   !--- set scaled missing values
   where (temp_r4 .EQfp. unscaled_missing)
         temp_i2 = MISSING_VALUE_INT2
   endwhere
 end subroutine SCALE_VECTOR_I2_RANK1

 subroutine SCALE_VECTOR_I2_RANK2(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i2)
   real, dimension(:,:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int2), dimension(:,:), intent(out):: temp_i2
   real, dimension(size(temp_r4,1),size(temp_r4,2)):: scratch_r4

    !---- linear
    if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = int(real(TWO_BYTE_MIN) + scratch_r4 * real(TWO_BYTE_MAX - TWO_BYTE_MIN),kind=int2)
    endif

    !---- log10
    if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      !temp_i2 = TWO_BYTE_MIN + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
      !              (unscaled_max - unscaled_min)) * (TWO_BYTE_MAX - TWO_BYTE_MIN)
      temp_i2 = int(real(TWO_BYTE_MIN) + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * real(TWO_BYTE_MAX - TWO_BYTE_MIN),kind=int2)
    endif

    !---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = int(real(TWO_BYTE_MIN) + sqrt(scratch_r4) * real(TWO_BYTE_MAX - TWO_BYTE_MIN),kind=int2)
    endif

    !--- set scaled missing values
    where (temp_r4 .EQfp. unscaled_missing)
         temp_i2 = MISSING_VALUE_INT2
    endwhere

 end subroutine SCALE_VECTOR_I2_RANK2

 subroutine SCALE_VECTOR_I2_RANK3(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i2)
   real, dimension(:,:,:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int2), dimension(:,:,:), intent(out):: temp_i2
   real, dimension(size(temp_r4,1),size(temp_r4,2),size(temp_r4,3)):: scratch_r4

    !---- linear
    if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = int(real(TWO_BYTE_MIN) + scratch_r4 * real(TWO_BYTE_MAX - TWO_BYTE_MIN),kind=int2)
    endif

    !---- log10
    if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i2 = int(real(TWO_BYTE_MIN) + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * real(TWO_BYTE_MAX - TWO_BYTE_MIN),kind=int2)
    endif

    !---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = int(real(TWO_BYTE_MIN) + sqrt(scratch_r4) * real(TWO_BYTE_MAX - TWO_BYTE_MIN),kind=int2)
    endif

    !--- set scaled missing values
    where (temp_r4 .EQfp. unscaled_missing)
         temp_i2 = MISSING_VALUE_INT2
    endwhere

 end subroutine SCALE_VECTOR_I2_RANK3

!-----------------------------------------------------------------------------------------------------
! VECTOR UNSCALING ROUTINES
!-----------------------------------------------------------------------------------------------------
 subroutine UNSCALE_VECTOR_I1_RANK1(temp_i1,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_r4)
   integer(kind=int1), dimension(:), intent(in):: temp_i1
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   real, dimension(:),  intent(out):: temp_r4
   real, dimension(size(temp_r4,1)):: scratch_r4
   integer (kind=int1):: scaled_min,scaled_max,scaled_missing

   scaled_min = int(ONE_BYTE_MIN,kind=int1)
   scaled_max = int(ONE_BYTE_MAX,kind=int1)
   scaled_missing = MISSING_VALUE_INT1

   scratch_r4 = 0.0
   !---- linear
   if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,real(temp_i1 - scaled_min)/real(scaled_max - scaled_min)))
       temp_r4 = unscaled_min + scratch_r4 * (unscaled_max - unscaled_min)
   endif
                                                                                                                                                         
   !---- log10
   if (iscaled == 2) then
       scratch_r4 = min(1.0,max(0.0,real(temp_i1 - scaled_min)/real(scaled_max - scaled_min)))
       temp_r4 = unscaled_min + scratch_r4 * (unscaled_max - unscaled_min)
       temp_r4 = 10**(temp_r4)
   endif

   !---- square root
   if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,real(temp_i1 - scaled_min)/real(scaled_max - scaled_min)))
       temp_r4 = unscaled_min + (scratch_r4**2) * (unscaled_max - unscaled_min)
   endif
                                                                                                                                                         
   !--- set scaled missing values
   where (temp_i1 == scaled_missing)
         temp_r4 = unscaled_missing
   endwhere
                                                                                                                                                         
 end subroutine UNSCALE_VECTOR_I1_RANK1


!----------------------------------------------------------------------------------------------------
! HDF4 WRITE ROUTINES
!-----------------------------------------------------------------------------------------------------
 subroutine WRITE_CLAVRX_HDF4_SDS_RANK1(Sd_Id,sds_data,Sds_Name,Sds_Type,scaled,sds_min,sds_max,&
                                        sds_units,sds_missing,sds_dim1,compress_flag,Istatus)
  real, intent(in), dimension(:):: sds_data
  real, intent(in):: sds_min, sds_max,sds_missing
  integer, intent(in):: Sd_Id,Sds_Type,compress_flag
  integer(kind=int1), intent(in):: scaled
  character(len=*), intent(in):: Sds_Name, sds_units,sds_dim1

  integer, intent(out):: Istatus
  integer:: scaled_min, scaled_max, scaled_missing
  integer:: sds_rank,units_len
  integer, dimension(1):: sds_dims,sds_edge,sds_chunk_size
  integer:: Sds_Id
  integer(kind=int4), dimension(2):: comp_prm
  integer(kind=int4):: comp_type
  integer(kind=int1),dimension(size(sds_data,1)):: temp_i1
  integer(kind=int2),dimension(size(sds_data,1)):: temp_i2
  integer(kind=int4),dimension(size(sds_data,1)):: temp_i4
  integer(kind=int4):: i,n

! HDF function declarations
 integer:: sfcreate, sfdimid, sfsdmname, sfwdata, sfendacc, &
           sfscatt, sfsnatt, sfschnk

  n = size(sds_data,1)
  sds_rank = 1
  sds_dims(1) = n
  sds_edge(1) = n
  units_len = len_trim(sds_units)
  Istatus = 0

!-------------------------------------------------------------
! define compression here
!-------------------------------------------------------------
 sds_chunk_size(1) = size(sds_data,1)     
 comp_type = 0                  !no compression
 comp_prm(1) = 0
 comp_prm(2) = 0

 if (compress_flag == 1) then  !gzip compression
   comp_type = 4
   comp_prm(1) = 6
   comp_prm(2) = 0
 endif

 if (compress_flag == 2) then  !szip compression
   comp_type = 5
   comp_prm(1) = 32
   comp_prm(2) = 2
 endif

!----------------------------------------------------------------------------
! create sds
!----------------------------------------------------------------------------

!-- write initial sds attributes
     Sds_Id = sfcreate(Sd_Id,Sds_Name,Sds_Type,sds_rank,sds_dims)
     Istatus = sfscatt(Sds_Id, "UNITS", DFNT_CHAR8, units_len, sds_units)
     Istatus = sfsnatt(Sds_Id, "SCALED", DFNT_INT8, 1, scaled) + Istatus
     Istatus = sfsnatt(Sds_Id, "RANGE_MISSING", DFNT_FLOAT32, 1, sds_missing) + Istatus
     Istatus = sfsdmname(sfdimid(Sds_Id,0),sds_dim1)

!--- compression and chunking
     Istatus = sfschnk(Sds_Id,sds_chunk_size,comp_type,comp_prm)+Istatus

!--- determined if a scaled sds, if so write needed attributes for scaling
     if (scaled > 0) then

!--- determine scaled ranges based on Sds_Type
       if (Sds_Type == DFNT_INT8) then
           scaled_min = ONE_BYTE_MIN
           scaled_max = ONE_BYTE_MAX
           scaled_missing = MISSING_VALUE_INT1
       elseif (Sds_Type == DFNT_INT16) then
           scaled_min = TWO_BYTE_MIN
           scaled_max = TWO_BYTE_MAX
           scaled_missing = MISSING_VALUE_INT2
       endif

!--- write remaining attributes
      Istatus = sfsnatt(Sds_Id, "RANGE_MIN", DFNT_FLOAT32, 1, sds_min) + Istatus
      Istatus = sfsnatt(Sds_Id, "RANGE_MAX", DFNT_FLOAT32, 1, sds_max) + Istatus
      Istatus = sfsnatt(Sds_Id, "SCALED_MIN", DFNT_INT32, 1, scaled_min) + Istatus
      Istatus = sfsnatt(Sds_Id, "SCALED_MAX", DFNT_INT32, 1, scaled_max) + Istatus
      Istatus = sfsnatt(Sds_Id, "SCALED_MISSING", DFNT_INT32, 1, scaled_missing) + Istatus

   endif

!-----------------------------------------------------------------------------------
! write data
!------------------------------------------------------------------------------------

!--- write unscaled arrays
 if (scaled == 0) then
   if (Sds_Type == DFNT_FLOAT32) then
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, sds_data) + Istatus
   endif

   if (Sds_Type == DFNT_INT8) then
    temp_i1 = int(sds_data,kind=int1)
    do i = 1,n
      if (sds_data(i) .EQfp. sds_missing) then
        temp_i1(i) = MISSING_VALUE_INT1
      endif
    end do
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i1) + Istatus
   endif

   if (Sds_Type == DFNT_INT16) then
    temp_i2 = int(sds_data,kind=int2)
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i2) + Istatus
   endif
   if (Sds_Type == DFNT_INT32) then
    temp_i4 = int(sds_data,kind=int4)
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i4) + Istatus
   endif
 endif

!--- write scaled arrays
 if (scaled > 0) then
   if (Sds_Type == DFNT_INT8) then
     call SCALE_VECTOR_I1_RANK1(sds_data,scaled,sds_min,sds_max,sds_missing,temp_i1)
     Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i1) + Istatus
   elseif (Sds_Type == DFNT_INT16) then
     call SCALE_VECTOR_I2_RANK1(sds_data,scaled,sds_min,sds_max,sds_missing,temp_i2) 
     Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i2) + Istatus
   else
    print *, "unsupported scaled / Sds_Type combination, stopping"
    stop
   endif
 endif

!--------------------------------------------------------------------
! close this record
!--------------------------------------------------------------------
   Istatus = sfendacc(Sds_Id) + Istatus

end subroutine WRITE_CLAVRX_HDF4_SDS_RANK1

!------------------------------------------------------------------------------------------------
!--- Write Routine for Rank=2
!-----------------------------------------------------------------------------------------------
subroutine WRITE_CLAVRX_HDF4_SDS_RANK2(Sd_Id,sds_data,Sds_Name,Sds_Type,scaled,sds_min,sds_max,&
                                    sds_units,sds_missing,sds_dim1,sds_dim2,compress_flag,Istatus)
  real, intent(in), dimension(:,:):: sds_data
  real, intent(in):: sds_min, sds_max,sds_missing
  integer, intent(in):: Sd_Id,Sds_Type,compress_flag
  integer(kind=int1), intent(in):: scaled
  character(len=*), intent(in):: Sds_Name, sds_units,sds_dim1,sds_dim2
  integer, intent(out):: Istatus

  real:: scaled_min, scaled_max, scaled_missing
  integer:: sds_rank,units_len
  integer, dimension(2):: sds_dims,sds_edge,sds_chunk_size
  integer:: Sds_Id
  integer(kind=int4), dimension(2):: comp_prm
  integer(kind=int4):: comp_type
  integer(kind=int1),dimension(size(sds_data,1),size(sds_data,2)):: temp_i1
  integer(kind=int2),dimension(size(sds_data,1),size(sds_data,2)):: temp_i2
  integer(kind=int4):: i1,i2,n1,n2

! HDF function declarations
 integer:: sfcreate, sfdimid, sfsdmname, sfwdata, sfendacc, &
           sfscatt, sfsnatt

  sds_rank = 1
  n1 = size(sds_data,1)
  n2 = size(sds_data,2)
  sds_dims(1) = n1
  sds_edge(1) = n1
  sds_dims(2) = n2
  sds_edge(2) = n2
  units_len = len_trim(sds_units)
  Istatus = 0

!-------------------------------------------------------------
! define compression here
!-------------------------------------------------------------
 sds_chunk_size(1) = size(sds_data,1)      
 sds_chunk_size(2) = size(sds_data,2)      
 comp_type = 0                  !no compression
 comp_prm(1) = 0
 comp_prm(2) = 0

 if (compress_flag == 1) then  !gzip compression
   comp_type = 4
   comp_prm(1) = 6
   comp_prm(2) = 0
 endif

 if (compress_flag == 2) then  !szip compression
   comp_type = 5
   comp_prm(1) = 32
   comp_prm(2) = 2
 endif
!----------------------------------------------------------------------------
! create sds
!----------------------------------------------------------------------------

!-- write initial sds attributes
     Sds_Id = sfcreate(Sd_Id,Sds_Name,Sds_Type,sds_rank,sds_dims)
     Istatus = sfsdmname(sfdimid(Sds_Id,0),sds_dim1)
     Istatus = sfsdmname(sfdimid(Sds_Id,1),sds_dim2)
     Istatus = sfscatt(Sds_Id, "UNITS", DFNT_CHAR8, units_len, sds_units)
     Istatus = sfsnatt(Sds_Id, "SCALED", DFNT_INT8, 1, scaled) + Istatus
     Istatus = sfsnatt(Sds_Id, "RANGE_MISSING", DFNT_FLOAT32, 1, sds_missing) + Istatus

!--- determined if a scaled sds, if so write needed attributes for scaling
     if (scaled > 0) then

!--- determine scaled ranges based on Sds_Type
       if (Sds_Type == DFNT_INT8) then
           scaled_min = ONE_BYTE_MIN
           scaled_max = ONE_BYTE_MAX
           scaled_missing = MISSING_VALUE_INT1
       elseif (Sds_Type == DFNT_INT16) then
           scaled_min = TWO_BYTE_MIN
           scaled_max = TWO_BYTE_MAX
           scaled_missing = MISSING_VALUE_INT2
       endif

!--- write remaining attributes
      Istatus = sfsnatt(Sds_Id, "RANGE_MIN", DFNT_FLOAT32, 1, sds_min) + Istatus
      Istatus = sfsnatt(Sds_Id, "RANGE_MAX", DFNT_FLOAT32, 1, sds_max) + Istatus
      Istatus = sfsnatt(Sds_Id, "SCALED_MIN", DFNT_INT32, 1, scaled_min) + Istatus
      Istatus = sfsnatt(Sds_Id, "SCALED_MAX", DFNT_INT32, 1, scaled_max) + Istatus
      Istatus = sfsnatt(Sds_Id, "SCALED_MISSING", DFNT_INT32, 1, scaled_missing) + Istatus

   endif

!-----------------------------------------------------------------------------------
! write data
!------------------------------------------------------------------------------------

!--- write unscaled arrays
 if (scaled == 0) then
   if (Sds_Type == DFNT_FLOAT32) then
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, sds_data) + Istatus
   endif
   if (Sds_Type == DFNT_INT8) then
    temp_i1 = int(sds_data,kind=int1)
    do i1 = 1, n1
      do i2 = 1, n2
        if (sds_data(i1,i2) .EQfp. sds_missing) then
            temp_i1(i1,i2) = MISSING_VALUE_INT1
        endif
      enddo
    enddo
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i1) + Istatus
   endif
   if (Sds_Type == DFNT_INT16) then
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, int(sds_data,kind=int2)) + Istatus
   endif
   if (Sds_Type == DFNT_INT32) then
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, int(sds_data,kind=int4)) + Istatus
   endif
 endif

!--- write scaled arrays
 if (scaled > 0) then
   if (Sds_Type == DFNT_INT8) then
     call SCALE_VECTOR_I1_RANK2(sds_data,scaled,sds_min,sds_max,sds_missing,temp_i1)
     Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i1) + Istatus
   elseif (Sds_Type == DFNT_INT16) then
     call SCALE_VECTOR_I2_RANK2(sds_data,scaled,sds_min,sds_max,sds_missing,temp_i2) 
     Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i2) + Istatus
   else
    print *, "unsupported scaled / Sds_Type combination, stopping"
    stop
   endif
 endif

!--------------------------------------------------------------------
! close this record
!--------------------------------------------------------------------
   Istatus = sfendacc(Sds_Id) + Istatus

end subroutine WRITE_CLAVRX_HDF4_SDS_RANK2

!--------------------------------------------------------------------------------
!---- Write Routine for Rank=3
!--------------------------------------------------------------------------------
subroutine WRITE_CLAVRX_HDF4_SDS_RANK3(Sd_Id,sds_data,Sds_Name,Sds_Type,scaled,sds_min,sds_max,&
                                    sds_units,sds_missing,sds_dim1,sds_dim2,sds_dim3,compress_flag,Istatus)
  real, intent(in), dimension(:,:,:):: sds_data
  real, intent(in):: sds_min, sds_max, sds_missing
  integer, intent(in):: Sd_Id,Sds_Type,compress_flag
  integer(kind=int1), intent(in):: scaled
  character(len=*), intent(in):: Sds_Name, sds_units,sds_dim1,sds_dim2,sds_dim3
  integer, intent(out):: Istatus

  real:: scaled_min, scaled_max, scaled_missing
  integer:: sds_rank,units_len
  integer, dimension(3):: sds_dims,sds_edge,sds_chunk_size
  integer:: Sds_Id
  integer(kind=int4), dimension(2):: comp_prm
  integer(kind=int4):: comp_type
  integer(kind=int1),dimension(size(sds_data,1),size(sds_data,2),size(sds_data,3)):: temp_i1
  integer(kind=int2),dimension(size(sds_data,1),size(sds_data,2),size(sds_data,3)):: temp_i2
  integer(kind=int4):: i1,i2,i3,n1,n2,n3

! HDF function declarations
 integer:: sfcreate, sfdimid, sfsdmname, sfwdata, sfendacc,  &
           sfscatt, sfsnatt

  sds_rank = 1

  n1 = size(sds_data,1)
  n2 = size(sds_data,2)
  n3 = size(sds_data,3)

  sds_dims(1) = n1
  sds_edge(1) = n1
  sds_dims(2) = n2
  sds_edge(2) = n2
  sds_dims(3) = n3
  sds_edge(3) = n3
  units_len = len_trim(sds_units)
  Istatus = 0

!-------------------------------------------------------------
! define compression here
!-------------------------------------------------------------
 sds_chunk_size(1) = size(sds_data,1)
 sds_chunk_size(2) = size(sds_data,2)
 sds_chunk_size(3) = size(sds_data,3)
 comp_type = 0                  !no compression
 comp_prm(1) = 0
 comp_prm(2) = 0

 if (compress_flag == 1) then  !gzip compression
   comp_type = 4
   comp_prm(1) = 6
   comp_prm(2) = 0
 endif

 if (compress_flag == 2) then  !szip compression
   comp_type = 5
   comp_prm(1) = 32
   comp_prm(2) = 2
 endif

!----------------------------------------------------------------------------
! create sds
!----------------------------------------------------------------------------

!-- write initial sds attributes
     Sds_Id = sfcreate(Sd_Id,Sds_Name,Sds_Type,sds_rank,sds_dims)
     Istatus = sfsdmname(sfdimid(Sds_Id,0),sds_dim1)
     Istatus = sfsdmname(sfdimid(Sds_Id,1),sds_dim2)
     Istatus = sfsdmname(sfdimid(Sds_Id,2),sds_dim3)
     Istatus = sfscatt(Sds_Id, "UNITS", DFNT_CHAR8, units_len, sds_units)
     Istatus = sfsnatt(Sds_Id, "SCALED", DFNT_INT8, 1, scaled) + Istatus
     Istatus = sfsnatt(Sds_Id, "RANGE_MISSING", DFNT_FLOAT32, 1, sds_missing) + Istatus

!--- determined if a scaled sds, if so write needed attributes for scaling
     if (scaled > 0) then

!--- determine scaled ranges based on Sds_Type
       if (Sds_Type == DFNT_INT8) then
           scaled_min = ONE_BYTE_MIN
           scaled_max = ONE_BYTE_MAX
           scaled_missing = MISSING_VALUE_INT1
       elseif (Sds_Type == DFNT_INT16) then
           scaled_min = TWO_BYTE_MIN
           scaled_max = TWO_BYTE_MAX
           scaled_missing = MISSING_VALUE_INT2
       endif

!--- write remaining attributes
      Istatus = sfsnatt(Sds_Id, "RANGE_MIN", DFNT_FLOAT32, 1, sds_min) + Istatus
      Istatus = sfsnatt(Sds_Id, "RANGE_MAX", DFNT_FLOAT32, 1, sds_max) + Istatus
      Istatus = sfsnatt(Sds_Id, "SCALED_MIN", DFNT_INT32, 1, scaled_min) + Istatus
      Istatus = sfsnatt(Sds_Id, "SCALED_MAX", DFNT_INT32, 1, scaled_max) + Istatus
      Istatus = sfsnatt(Sds_Id, "SCALED_MISSING", DFNT_INT32, 1, scaled_missing) + Istatus

   endif

!-----------------------------------------------------------------------------------
! write data
!------------------------------------------------------------------------------------

!--- write unscaled arrays
 if (scaled == 0) then
   if (Sds_Type == DFNT_FLOAT32) then
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, sds_data) + Istatus
   endif
   if (Sds_Type == DFNT_INT8) then
    temp_i1 = int(sds_data,kind=int1)
    do i1 = 1, n1
     do i2 = 1, n2
       do i3 = 1, n3
        if (sds_data(i1,i2,i3) .EQfp. sds_missing) then
            temp_i1(i2,i2,i3) = MISSING_VALUE_INT1
        endif
       enddo
      enddo
    enddo
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i1) + Istatus
 
   endif
   if (Sds_Type == DFNT_INT16) then
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, int(sds_data,kind=int2)) + Istatus
   endif
   if (Sds_Type == DFNT_INT32) then
    Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, int(sds_data,kind=int4)) + Istatus
   endif
 endif

!--- write scaled arrays
 if (scaled > 0) then
   if (Sds_Type == DFNT_INT8) then
     call SCALE_VECTOR_I1_RANK3(sds_data,scaled,sds_min,sds_max,sds_missing,temp_i1)
     Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i1) + Istatus
   elseif (Sds_Type == DFNT_INT16) then
     call SCALE_VECTOR_I2_RANK3(sds_data,scaled,sds_min,sds_max,sds_missing,temp_i2) 
     Istatus = sfwdata(Sds_Id, 0, 1, sds_edge, temp_i2) + Istatus
   else
    print *, "unsupported scaled / Sds_Type combination, stopping"
    stop
   endif
 endif

!--------------------------------------------------------------------
! close this record
!--------------------------------------------------------------------
   Istatus = sfendacc(Sds_Id) + Istatus

end subroutine WRITE_CLAVRX_HDF4_SDS_RANK3
!----------------------------------------------------------------------------------------------------
! HDF4 READ ROUTINES
! this routine reads in level3 sds's, it assumes
! 1. - you know the size of the array (num_cells_with_data)
!-----------------------------------------------------------------------------------------------------
 subroutine READ_CLAVRX_HDF4_SDS_RANK1(Sd_Id,sds_dim_input,Sds_Name,sds_data, &
                                       sds_data_type,scaled,sds_units, &
                                       unscaled_min,unscaled_max,unscaled_missing,&
                                       scaled_min,scaled_max,scaled_missing,&
                                       Istatus)
  integer(kind=int4), intent(in):: Sd_Id,sds_dim_input
  character(len=*), intent(in):: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  real, intent(out):: unscaled_min, unscaled_max, unscaled_missing
  integer(kind=int4), intent(out):: scaled_min, scaled_max, scaled_missing
  real, intent(out), dimension(:):: sds_data
  integer(kind=int1), intent(out):: scaled
  character(len=*), intent(out):: sds_units
  integer(kind=int4), intent(out):: sds_data_type,Istatus
  integer(kind=int4):: Sds_Id,sds_dim1
  integer(kind=int4):: num_attrs, sds_rank
  integer(kind=int4), dimension(1):: dimsizes

  real, dimension(size(sds_data)):: sds_data_temp

  integer(kind=4), dimension(1):: sds_dims_1d, sds_start_1d, sds_stride_1d, sds_edges_1d


  integer(kind=int1), dimension(:), allocatable:: temp_i1
  integer(kind=int2), dimension(:), allocatable:: temp_i2
  integer(kind=int4), dimension(:), allocatable:: temp_i4
  real(kind=real4), dimension(:), allocatable:: temp_r4

  Istatus = 0

!----------------------------------------------------------------------------
! open sds for reading
!----------------------------------------------------------------------------
Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,Sds_Name))

Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, dimsizes, sds_data_type, num_attrs) + Istatus
sds_dim1 = dimsizes(1)
sds_dims_1d = (/sds_dim1/)
sds_start_1d = (/ 0 /)      
sds_stride_1d = (/ 1 /)
sds_edges_1d = (/ sds_dim1 /)   
sds_units = " "

!--- read sds attributes
 Istatus = sfrnatt(Sds_Id, sffattr(Sds_Id,"SCALED"), scaled)
 Istatus = sfrcatt(Sds_Id, sffattr(Sds_Id,"UNITS"), sds_units)
 Istatus = sfrnatt(Sds_Id, sffattr(Sds_Id,"RANGE_MISSING"), unscaled_missing)

!-- if scaled, read attributes that allow unscaling
 if (scaled > 0) then
   Istatus = sfrnatt(Sds_Id, sffattr(Sds_Id,"SCALED_MISSING"), scaled_missing)
   Istatus = sfrnatt(Sds_Id, sffattr(Sds_Id,"SCALED_MIN"), scaled_min)
   Istatus = sfrnatt(Sds_Id, sffattr(Sds_Id,"SCALED_MAX"), scaled_max)
   Istatus = sfrnatt(Sds_Id, sffattr(Sds_Id,"RANGE_MIN"), unscaled_min)
   Istatus = sfrnatt(Sds_Id, sffattr(Sds_Id,"RANGE_MAX"), unscaled_max)
 endif

!-- check dimension against expectations
    if (sds_dim_input /= sds_dim1) then
      print *, "error, sds dimension differs from expectations, stopping", sds_dim_input, sds_dim1
      stop
     endif

!--- allocate arrays for holding data, read data and store in output array
if (sds_data_type == DFNT_INT8) then
    allocate(temp_i1(sds_dim1))
    Istatus = sfrdata(Sds_Id, sds_start_1d, sds_stride_1d, sds_edges_1d, temp_i1) + Istatus   
    sds_data = real(temp_i1)
elseif (sds_data_type == DFNT_INT16) then
    allocate(temp_i2(sds_dim1))
    Istatus = sfrdata(Sds_Id, sds_start_1d, sds_stride_1d, sds_edges_1d, temp_i2) + Istatus   
    sds_data = real(temp_i2)
elseif (sds_data_type == DFNT_INT32) then
    allocate(temp_i4(sds_dim1))
    Istatus = sfrdata(Sds_Id, sds_start_1d, sds_stride_1d, sds_edges_1d, temp_i4) + Istatus   
    sds_data = real(temp_i4)
elseif (sds_data_type == DFNT_FLOAT32) then
    allocate(temp_r4(sds_dim1))
    Istatus = sfrdata(Sds_Id, sds_start_1d, sds_stride_1d, sds_edges_1d, temp_r4) + Istatus   
    sds_data = temp_r4
else
    print *, "attempt to read unsupported data type, stopping"
    stop
endif

!---deallocate temp arrays
   if (allocated(temp_i1)) deallocate(temp_i1)
   if (allocated(temp_i2)) deallocate(temp_i2)
   if (allocated(temp_i4)) deallocate(temp_i4)
   if (allocated(temp_r4)) deallocate(temp_r4)

!--- close sds
   Istatus = sfendacc(Sds_Id) + Istatus

!--- unscale sds

if (scaled > 0) then 

    sds_data_temp = sds_data
!---- linear
    if (scaled == 1) then
       sds_data_temp = min(1.0,max(0.0,real(sds_data_temp - scaled_min)/real(scaled_max - scaled_min)))
       sds_data_temp = unscaled_min + sds_data_temp * (unscaled_max - unscaled_min)
    endif

!---- log10
    if (scaled == 2) then
       sds_data_temp = min(1.0,max(0.0,real(sds_data_temp - scaled_min)/real(scaled_max - scaled_min)))
       sds_data_temp = unscaled_min + sds_data_temp * (unscaled_max - unscaled_min)
       sds_data_temp = 10**(sds_data_temp)
     endif

!---- square root
    if (scaled == 3) then
       sds_data_temp = min(1.0,max(0.0,real(sds_data_temp - scaled_min)/real(scaled_max - scaled_min)))
       sds_data_temp = unscaled_min + (sds_data_temp**2) * (unscaled_max - unscaled_min)
    endif

!--- set scaled missing values
    !where (sds_data == scaled_missing)
    where (sds_data .EQfp. real(scaled_missing))
         sds_data_temp = unscaled_missing
    endwhere

    sds_data = sds_data_temp

 endif
 
end subroutine READ_CLAVRX_HDF4_SDS_RANK1

!---------------------------------------------------------------------------------------------------------
! Begin of New Routines to this Module
!---------------------------------------------------------------------------------------------------------

!=======================================================================
! read a global numerical attribute
!=======================================================================
function READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id,Attribute_Name) result(Attr_Value)

integer, intent(in):: Sd_Id
character(len=*), intent(in):: Attribute_Name
real(kind=real4):: Attr_Value
integer:: Istatus
integer:: Attr_Id
character(len=100):: Attr_Name
integer:: Attr_Type
integer:: Attr_Count
integer(kind=int1):: Attr_Value_I1
integer(kind=int2):: Attr_Value_I2
integer(kind=int4):: Attr_Value_I4
real(kind=real4):: Attr_Value_R4
real(kind=real8):: Attr_Value_R8

!--- hdf commands
integer:: sffattr
integer:: sfrnatt
integer:: sfgainfo

!------------------------------------------------------------------------------------
! begin executable code
!------------------------------------------------------------------------------------

    Istatus = 0

    !--- open attribute
    Attr_Id = sffattr(Sd_Id,trim(Attribute_Name))
    if (Attr_Id < 0) then
         Istatus = 1
         Attr_Value = Missing_Value_Real4
         return
    endif

    Istatus = sfgainfo(Sd_Id,Attr_Id,Attr_Name,Attr_Type,Attr_Count) + Istatus

    if (Attr_Type == DFNT_INT8) then
       Istatus = sfrnatt(Sd_Id, Attr_Id, Attr_Value_I1) + Istatus
       Attr_Value = real(Attr_Value_I1)
    elseif (Attr_Type == DFNT_INT16) then
       Istatus = sfrnatt(Sd_Id, Attr_Id, Attr_Value_I2) + Istatus
       Attr_Value = real(Attr_Value_I2)
    elseif (Attr_Type == DFNT_INT32) then
       Istatus = sfrnatt(Sd_Id, Attr_Id, Attr_Value_I4) + Istatus
       Attr_Value = real(Attr_Value_I4)
    elseif (Attr_Type == DFNT_FLOAT32) then
       Istatus = sfrnatt(Sd_Id, Attr_Id, Attr_Value_R4) + Istatus
       Attr_Value = Attr_Value_R4
    elseif (Attr_Type == DFNT_FLOAT64) then
       Istatus = sfrnatt(Sd_Id, Attr_Id, Attr_Value_R8) + Istatus
       Attr_Value = real(Attr_Value_R8)
    else
       print *, "Unknown Attribute Type, stopping"
       stop
    endif

    !--- check for final error status
    if (Istatus /= 0) then
       print *, "Errors on read of ", trim(Attribute_Name)
    endif

end function READ_HDF_GLOBAL_ATTRIBUTE_NUM

!=======================================================================
! read a global string attribute
!=======================================================================
function READ_HDF_GLOBAL_ATTRIBUTE_STR(Sd_Id,Attribute_Name) result(Attr_Value)

integer, intent(in):: Sd_Id
character(len=*), intent(in):: Attribute_Name
character(len=500):: Attr_Value
integer:: Istatus
integer:: Attr_Id
!character(len=100):: Attr_Name

!--- hdf commands
integer:: sffattr
integer:: sfrnatt

!------------------------------------------------------------------------------------
! begin executable code
!------------------------------------------------------------------------------------

    Istatus = 0

    !--- open attribute
    Attr_Id = sffattr(Sd_Id,trim(Attribute_Name))
    if (Attr_Id < 0) then
         Istatus = 1
         Attr_Value = ''
         return
    endif

    !--- read attribute
    Istatus = sfrnatt(Sd_Id, Attr_Id, Attr_Value) + Istatus

    !--- check for final error status
    if (Istatus /= 0) then
       print *, "Errors on read of ", trim(Attribute_Name)
    endif

end function READ_HDF_GLOBAL_ATTRIBUTE_STR

!=======================================================================
! read a one dimensional sds
!=======================================================================
subroutine READ_CLAVRX_HDF_SDS_1D(file_name,Sds_Name,Unscaled_Sds_Data,Istatus,quiet, &
                                  Sds_Start,Sds_Stride,Sds_Edges)

character(len=*), intent(in):: file_name
character(len=*), intent(in):: Sds_Name
real(kind=real4), intent(out),dimension(:), allocatable:: Unscaled_Sds_Data
character(len=*), intent(in), optional:: quiet

integer:: Sd_Id
integer, dimension(1):: Sds_Dims
integer, optional, dimension(1):: Sds_Start
integer, optional, dimension(1):: Sds_Stride
integer, optional, dimension(1):: Sds_Edges
integer(kind=int1), dimension(:), allocatable:: Temp_I1
integer(kind=int2), dimension(:), allocatable:: Temp_I2
integer(kind=int4), dimension(:), allocatable:: Temp_I4
real(kind=real4), dimension(:), allocatable:: Temp_R4
real(kind=real4), dimension(:), allocatable:: Scaled_Sds_Data
integer:: Num_Attrs
integer, intent(out):: Istatus
integer:: sds_index
integer:: Attr_Index
character(72):: Sds_Name_Temp
integer(kind=int1):: dummy_i1
integer(kind=int2):: dummy_i2

type(Sds_Struct):: Sds

!------------------------------------------------------------------------------------
! begin executable code
!------------------------------------------------------------------------------------
Sds%Sds_Name = Sds_Name

Istatus = 0

!--- open the file
Sd_Id = sfstart(trim(file_name), DFACC_READ)

!--- handle a failure of the opening
if (Sd_Id < 0) then
    Istatus = 1
    return
endif

!--- find sds in the file
sds_index = sfn2index(Sd_Id,trim(Sds%Sds_Name))

!--- handle failure of finding sds
if (sds_index < 0) then
    Istatus = 1
    if (.not. present(quiet)) print *, "Error could not find sds named ", trim(Sds%Sds_Name)
    return
endif

!--- open the sds
Sds%Sds_Idx = sfselect(Sd_Id, sds_index)

!--- handle failure of opening the sds
if (Sds%Sds_Idx < 0) then
    Istatus = 1
    if (.not. present(quiet)) print *, "Error opening sds named ", trim(Sds%Sds_Name)
    return
endif


!--- get information on this file
Istatus = sfginfo(Sds%Sds_Idx, Sds_Name_Temp, Sds%Rank, Sds_Dims, Sds%Data_Type, Num_Attrs) + Istatus


!---
if (.not. present(Sds_Start)) Sds_Start = (/0/)
if (.not. present(Sds_Stride)) Sds_Stride = (/ 1 /)
if (.not. present(Sds_Edges)) Sds_Edges = Sds_Dims

!--- set missing value to default
Sds%Actual_Missing = DEFAULT_MISSING_VALUE

!--- read scaling attribute (0=no scaling, 1 = linear)
Attr_Index = sffattr(Sds%Sds_Idx,"SCALED")
if (Attr_Index == -1) then
    Sds%Scaling_Type = 0_int1
else
    Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, Sds%Scaling_Type) + Istatus
endif

!--- scale factor
Sds%Scale_Factor = 1.0
if (Sds%Scaling_Type == 1) then 
    Attr_Index = sffattr(Sds%Sds_Idx,"scale_factor")
    if (Attr_Index == -1) then
       if (.not. present(quiet)) print *, "scale factor missing, assumed 1.0 for sds named ", trim(Sds%Sds_Name)
    else
       Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, Sds%Scale_Factor) + Istatus
    endif
endif

!--- add_offset
Sds%Add_Offset = 0.0
if (Sds%Scaling_Type == 1) then 
    Attr_Index = sffattr(Sds%Sds_Idx,"add_offset")
    if (Attr_Index == -1) then
       if (.not. present(quiet)) print *, "add offset missing, assumed 0.0 for sds named ", trim(Sds%Sds_Name)
    else
       Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, Sds%Add_Offset) + Istatus
    endif
endif

!--- read fill value
Sds%Fill_Value = 0
if (Sds%Scaling_Type == 1) then 
    Attr_Index = sffattr(Sds%Sds_Idx,'_FillValue')
    if (Attr_Index == -1) then
       if (.not. present(quiet)) print *, "fill value missing, assumed for sds named ", trim(Sds%Sds_Name)
       if (Sds%Data_Type == DFNT_INT8) Sds%Fill_Value = -128.0
       if (Sds%Data_Type == DFNT_INT16) Sds%Fill_Value = -32768.0
    else
       if (Sds%Data_Type == DFNT_INT8) then
           Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, dummy_i1) + Istatus
           Sds%Fill_Value = real(dummy_i1)
       endif
       if (Sds%Data_Type == DFNT_INT16) then
           Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, dummy_i2) + Istatus
           Sds%Fill_Value = real(dummy_i2)
       endif
    endif
endif


!--- allocate arrays for holding data, read data and store in output array
allocate(Scaled_Sds_Data(Sds_Dims(1)))
allocate(Unscaled_Sds_Data(Sds_Dims(1)))

if (Sds%Data_Type == DFNT_INT8) then
       allocate(Temp_I1(Sds_Dims(1)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
       Scaled_Sds_Data = real(Temp_I1)
elseif (Sds%Data_Type == DFNT_INT16) then
       allocate(Temp_I2(Sds_Dims(1)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
       Scaled_Sds_Data = real(Temp_I2)
elseif (Sds%Data_Type == DFNT_INT32) then
       allocate(Temp_I4(Sds_Dims(1)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
       Scaled_Sds_Data = real(Temp_I4)
elseif (Sds%Data_Type == DFNT_FLOAT32) then
       allocate(Temp_R4(Sds_Dims(1)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
       Scaled_Sds_Data = Temp_R4
else
       if (.not. present(quiet)) print *, "Possibly fatal error at location 2:"
       if (.not. present(quiet)) print *, "data type was ", Sds%Data_Type
       if (.not. present(quiet)) print *, "sds complete value is "
!      if (.not. present(quiet)) print *, Sds
       print *, "attempt to read unsupported data type, stopping"
       Scaled_Sds_Data = MISSING_VALUE_REAL4
       Sds%Data_Type = -999
       stop
endif

!---deallocate temp arrays
if (allocated(Temp_I1)) deallocate(Temp_I1)
if (allocated(Temp_I2)) deallocate(Temp_I2)
if (allocated(Temp_I4)) deallocate(Temp_I4)
if (allocated(Temp_R4)) deallocate(Temp_R4)

!--- close Sds
Istatus = sfendacc(Sds%Sds_Idx) + Istatus

!----- close input hdf file
Istatus = sfend(Sd_Id) + Istatus

!--- unscale Sds
Unscaled_Sds_Data = Scaled_Sds_Data * Sds%Scale_Factor + Sds%Add_Offset

!--- set scaled missing values (unless its a packed data set)
if (Sds%Scaling_Type /= 0 .or. index(Sds%Sds_Name, 'packed') > 0) then
 where (Scaled_Sds_Data .EQfp. Sds%Fill_Value)
    Unscaled_Sds_Data = Sds%Actual_Missing
 endwhere
endif

if (allocated(Scaled_Sds_Data)) deallocate(Scaled_Sds_Data)

end subroutine READ_CLAVRX_HDF_SDS_1D


!=======================================================================
! read a two dimensional sds
!=======================================================================
subroutine READ_CLAVRX_HDF_SDS_2D(file_name,Sds_Name,Unscaled_Sds_Data,Istatus,quiet, &
                                  Sds_Start,Sds_Stride,Sds_Edges)

character(len=*), intent(in):: file_name
character(len=*), intent(in):: Sds_Name
real(kind=real4), intent(out),dimension(:,:), allocatable:: Unscaled_Sds_Data
character(len=*), intent(in), optional:: quiet

integer:: Sd_Id
integer, dimension(2):: Sds_Dims
integer, optional, dimension(2):: Sds_Start
integer, optional, dimension(2):: Sds_Stride
integer, optional, dimension(2):: Sds_Edges
integer(kind=int1), dimension(:,:), allocatable:: Temp_I1
integer(kind=int2), dimension(:,:), allocatable:: Temp_I2
integer(kind=int4), dimension(:,:), allocatable:: Temp_I4
real(kind=real4), dimension(:,:), allocatable:: Temp_R4
real(kind=real4), dimension(:,:), allocatable:: Scaled_Sds_Data
integer:: Num_Attrs
integer, intent(out):: Istatus
integer:: sds_index
integer:: Attr_Index
character(72):: Sds_Name_Temp
integer(kind=int1):: dummy_i1
integer(kind=int2):: dummy_i2

type(Sds_Struct):: Sds

!------------------------------------------------------------------------------------
! begin executable code
!------------------------------------------------------------------------------------
Sds%Sds_Name = Sds_Name

Istatus = 0

!--- open the file
Sd_Id = sfstart(trim(file_name), DFACC_READ)

!--- handle a failure of the opening
if (Sd_Id < 0) then
    Istatus = 1
    return
endif

!--- find sds in the file
sds_index = sfn2index(Sd_Id,trim(Sds%Sds_Name))

!--- handle failure of finding sds
if (sds_index < 0) then
    Istatus = 1
    if (.not. present(quiet)) print *, "Error could not find sds named ", trim(Sds%Sds_Name)
    return
endif

!--- open the sds
Sds%Sds_Idx = sfselect(Sd_Id, sds_index)

!--- handle failure of opening the sds
if (Sds%Sds_Idx < 0) then
    Istatus = 1
    if (.not. present(quiet)) print *, "Error opening sds named ", trim(Sds%Sds_Name)
    return
endif

!--- get information on this file
Istatus = sfginfo(Sds%Sds_Idx, Sds_Name_Temp, Sds%Rank, Sds_Dims, Sds%Data_Type, Num_Attrs) + Istatus

!---
if (.not. present(Sds_Start)) Sds_Start  = (/ 0, 0 /)
if (.not. present(Sds_Stride)) Sds_Stride = (/ 1, 1 /)
if (.not. present(Sds_Edges)) Sds_Edges = Sds_Dims

!--- set missing value to default
Sds%Actual_Missing = DEFAULT_MISSING_VALUE

!--- read scaling attribute (0=no scaling, 1 = linear)
Attr_Index = sffattr(Sds%Sds_Idx,"SCALED")
if (Attr_Index == -1) then
    Sds%Scaling_Type = 0_int1
else
    Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, Sds%Scaling_Type) + Istatus
endif

!--- scale factor
Sds%Scale_Factor = 1.0
if (Sds%Scaling_Type == 1) then 
    Attr_Index = sffattr(Sds%Sds_Idx,"scale_factor")
    if (Attr_Index == -1) then
       if (.not. present(quiet)) print *, "scale factor missing, assumed 1.0 for sds named ", trim(Sds%Sds_Name)
    else
       Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, Sds%Scale_Factor) + Istatus
    endif
endif

!--- add_offset
Sds%Add_Offset = 0.0
if (Sds%Scaling_Type == 1) then 
    Attr_Index = sffattr(Sds%Sds_Idx,"add_offset")
    if (Attr_Index == -1) then
       if (.not. present(quiet)) print *, "add offset missing, assumed 0.0 for sds named ", trim(Sds%Sds_Name)
    else
       Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, Sds%Add_Offset) + Istatus
    endif
endif

!--- read fill value
Sds%Fill_Value = 0
if (Sds%Scaling_Type == 1) then 
    Attr_Index = sffattr(Sds%Sds_Idx,'_FillValue')
    if (Attr_Index == -1) then
       if (.not. present(quiet)) print *, "fill value missing, assumed for sds named ", trim(Sds%Sds_Name)
       if (Sds%Data_Type == DFNT_INT8) Sds%Fill_Value = -128
       if (Sds%Data_Type == DFNT_INT16) Sds%Fill_Value = -32768
    else
       if (Sds%Data_Type == DFNT_INT8) then
           Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, dummy_i1) + Istatus
           Sds%Fill_Value = real(dummy_i1)
       endif
       if (Sds%Data_Type == DFNT_INT16) then
           Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, dummy_i2) + Istatus
           Sds%Fill_Value = real(dummy_i2)
       endif
    endif
endif


!--- allocate arrays for holding data, read data and store in output array
allocate(Scaled_Sds_Data(Sds_Dims(1),Sds_Dims(2)))
allocate(Unscaled_Sds_Data(Sds_Dims(1),Sds_Dims(2)))

if (Sds%Data_Type == DFNT_INT8) then
       allocate(Temp_I1(Sds_Dims(1),Sds_Dims(2)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
       Scaled_Sds_Data = real(Temp_I1)
elseif (Sds%Data_Type == DFNT_INT16) then
       allocate(Temp_I2(Sds_Dims(1), Sds_Dims(2)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
       Scaled_Sds_Data = real(Temp_I2)
elseif (Sds%Data_Type == DFNT_INT32) then
       allocate(Temp_I4(Sds_Dims(1),Sds_Dims(2)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
       Scaled_Sds_Data = real(Temp_I4)
elseif (Sds%Data_Type == DFNT_FLOAT32) then
       allocate(Temp_R4(Sds_Dims(1),Sds_Dims(2)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
       Scaled_Sds_Data = Temp_R4
else
       if (.not. present(quiet)) print *, "Possibly fatal error at location 2:"
       if (.not. present(quiet)) print *, "data type was ", Sds%Data_Type
       if (.not. present(quiet)) print *, "sds complete value is "
!      if (.not. present(quiet)) print *, Sds
       print *, "attempt to read unsupported data type, stopping"
       Scaled_Sds_Data = MISSING_VALUE_REAL4
       Sds%Data_Type = -999
       stop
endif

!---deallocate temp arrays
if (allocated(Temp_I1)) deallocate(Temp_I1)
if (allocated(Temp_I2)) deallocate(Temp_I2)
if (allocated(Temp_I4)) deallocate(Temp_I4)
if (allocated(Temp_R4)) deallocate(Temp_R4)

!--- close Sds
Istatus = sfendacc(Sds%Sds_Idx) + Istatus

!----- close input hdf file
Istatus = sfend(Sd_Id) + Istatus

!--- unscale Sds
Unscaled_Sds_Data = Scaled_Sds_Data * Sds%Scale_Factor + Sds%Add_Offset

!--- set scaled missing values (unless its a packed data set)
if (Sds%Scaling_Type /= 0 .or. index(Sds%Sds_Name, 'packed') > 0) then
 where (Scaled_Sds_Data .EQfp. Sds%Fill_Value)
    Unscaled_Sds_Data = Sds%Actual_Missing
 endwhere
endif

if (allocated(Scaled_Sds_Data)) deallocate(Scaled_Sds_Data)

end subroutine READ_CLAVRX_HDF_SDS_2D

!=======================================================================
! read a three dimensional sds
!=======================================================================
subroutine READ_CLAVRX_HDF_SDS_3D(file_name,Sds_Name,Unscaled_Sds_Data,Istatus,quiet, &
                                  Sds_Start,Sds_Stride,Sds_Edges)

character(len=*), intent(in):: file_name
character(len=*), intent(in):: Sds_Name
real(kind=real4), intent(out),dimension(:,:,:), allocatable:: Unscaled_Sds_Data
character(len=*), intent(in), optional:: quiet

integer:: Sd_Id
integer, dimension(3):: Sds_Dims
integer, optional, dimension(3):: Sds_Start
integer, optional, dimension(3):: Sds_Stride
integer, optional, dimension(3):: Sds_Edges
integer(kind=int1), dimension(:,:,:), allocatable:: Temp_I1
integer(kind=int2), dimension(:,:,:), allocatable:: Temp_I2
integer(kind=int4), dimension(:,:,:), allocatable:: Temp_I4
real(kind=real4), dimension(:,:,:), allocatable:: Temp_R4
real(kind=real4), dimension(:,:,:), allocatable:: Scaled_Sds_Data
integer:: Num_Attrs
integer, intent(out):: Istatus
integer:: sds_index
integer:: Attr_Index
character(72):: Sds_Name_Temp
integer(kind=int1):: dummy_i1
integer(kind=int2):: dummy_i2

type(Sds_Struct):: Sds

!------------------------------------------------------------------------------------
! begin executable code
!------------------------------------------------------------------------------------
Sds%Sds_Name = Sds_Name

Istatus = 0

!--- open the file
Sd_Id = sfstart(trim(file_name), DFACC_READ)

!--- handle a failure of the opening
if (Sd_Id < 0) then
    Istatus = 1
    return
endif

!--- find sds in the file
sds_index = sfn2index(Sd_Id,trim(Sds%Sds_Name))

!--- handle failure of finding sds
if (sds_index < 0) then
    Istatus = 1
    if (.not. present(quiet)) print *, "Error could not find sds named ", trim(Sds%Sds_Name)
    return
endif

!--- open the sds
Sds%Sds_Idx = sfselect(Sd_Id, sds_index)

!--- handle failure of opening the sds
if (Sds%Sds_Idx < 0) then
    Istatus = 1
    if (.not. present(quiet)) print *, "Error opening sds named ", trim(Sds%Sds_Name)
    return
endif

!--- get information on this file
Istatus = sfginfo(Sds%Sds_Idx, Sds_Name_Temp, Sds%Rank, Sds_Dims, Sds%Data_Type, Num_Attrs) + Istatus

!---
if (.not. present(Sds_Start)) Sds_Start  = (/ 0, 0, 0 /)
if (.not. present(Sds_Stride)) Sds_Stride = (/ 1, 1, 1 /)
if (.not. present(Sds_Edges)) Sds_Edges = Sds_Dims

!--- set missing value to default
Sds%Actual_Missing = DEFAULT_MISSING_VALUE

!--- read scaling attribute (0=no scaling, 1 = linear)
Attr_Index = sffattr(Sds%Sds_Idx,"SCALED")
if (Attr_Index == -1) then
    Sds%Scaling_Type = 0_int1
else
    Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, Sds%Scaling_Type) + Istatus
endif

!--- scale factor
Sds%Scale_Factor = 1.0
if (Sds%Scaling_Type == 1) then 
    Attr_Index = sffattr(Sds%Sds_Idx,"scale_factor")
    if (Attr_Index == -1) then
       if (.not. present(quiet)) print *, "scale factor missing, assumed 1.0 for sds named ", trim(Sds%Sds_Name)
    else
       Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, Sds%Scale_Factor) + Istatus
    endif
endif

!--- add_offset
Sds%Add_Offset = 0.0
if (Sds%Scaling_Type == 1) then 
    Attr_Index = sffattr(Sds%Sds_Idx,"add_offset")
    if (Attr_Index == -1) then
       if (.not. present(quiet)) print *, "add offset missing, assumed 0.0 for sds named ", trim(Sds%Sds_Name)
    else
       Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, Sds%Add_Offset) + Istatus
    endif
endif

!--- read fill value
Sds%Fill_Value = 0
if (Sds%Scaling_Type == 1) then 
    Attr_Index = sffattr(Sds%Sds_Idx,'_FillValue')
    if (Attr_Index == -1) then
       if (.not. present(quiet)) print *, "fill value missing, assumed for sds named ", trim(Sds%Sds_Name)
       if (Sds%Data_Type == DFNT_INT8) Sds%Fill_Value = -128
       if (Sds%Data_Type == DFNT_INT16) Sds%Fill_Value = -32768
    else
       if (Sds%Data_Type == DFNT_INT8) then
           Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, dummy_i1) + Istatus
           Sds%Fill_Value = real(dummy_i1)
       endif
       if (Sds%Data_Type == DFNT_INT16) then
           Istatus = sfrnatt(Sds%Sds_Idx, Attr_Index, dummy_i2) + Istatus
           Sds%Fill_Value = real(dummy_i2)
       endif
    endif
endif


!--- allocate arrays for holding data, read data and store in output array
allocate(Scaled_Sds_Data(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
allocate(Unscaled_Sds_Data(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))

if (Sds%Data_Type == DFNT_INT8) then
       allocate(Temp_I1(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_I1) + Istatus
       Scaled_Sds_Data = real(Temp_I1)
elseif (Sds%Data_Type == DFNT_INT16) then
       allocate(Temp_I2(Sds_Dims(1), Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_I2) + Istatus
       Scaled_Sds_Data = real(Temp_I2)
elseif (Sds%Data_Type == DFNT_INT32) then
       allocate(Temp_I4(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_I4) + Istatus
       Scaled_Sds_Data = real(Temp_I4)
elseif (Sds%Data_Type == DFNT_FLOAT32) then
       allocate(Temp_R4(Sds_Dims(1),Sds_Dims(2),Sds_Dims(3)))
       Istatus = sfrdata(Sds%Sds_Idx, Sds_Start, Sds_Stride, Sds_Edges, Temp_R4) + Istatus
       Scaled_Sds_Data = Temp_R4
else
       if (.not. present(quiet)) print *, "Possibly fatal error at location 2:"
       if (.not. present(quiet)) print *, "data type was ", Sds%Data_Type
       if (.not. present(quiet)) print *, "sds complete value is "
!      if (.not. present(quiet)) print *, Sds
       print *, "attempt to read unsupported data type, stopping"
       Scaled_Sds_Data = MISSING_VALUE_REAL4
       Sds%Data_Type = -999
       stop
endif

!---deallocate temp arrays
if (allocated(Temp_I1)) deallocate(Temp_I1)
if (allocated(Temp_I2)) deallocate(Temp_I2)
if (allocated(Temp_I4)) deallocate(Temp_I4)
if (allocated(Temp_R4)) deallocate(Temp_R4)

!--- close Sds
Istatus = sfendacc(Sds%Sds_Idx) + Istatus

!----- close input hdf file
Istatus = sfend(Sd_Id) + Istatus

!--- unscale Sds
Unscaled_Sds_Data = Scaled_Sds_Data * Sds%Scale_Factor + Sds%Add_Offset

!--- set scaled missing values (unless its a packed data set)
if (Sds%Scaling_Type /= 0 .or. index(Sds%Sds_Name, 'packed') > 0) then
 where (Scaled_Sds_Data .EQfp. Sds%Fill_Value)
    Unscaled_Sds_Data = Sds%Actual_Missing
 endwhere
endif

if (allocated(Scaled_Sds_Data)) deallocate(Scaled_Sds_Data)

end subroutine READ_CLAVRX_HDF_SDS_3D

!======================================================================
! routine to write unscaled 1d sds 
!======================================================================
subroutine  WRITE_CLAVRX_HDF_SDS_1D(Sd_Id,Sds_Name,Sds_Data,Istatus)

integer, intent(in):: Sd_Id
character (len=*), intent(in):: Sds_Name
real, dimension(:),intent(in):: Sds_Data
integer, intent(out):: Istatus
integer, dimension(1):: Sds_Dims
integer, dimension(1):: Sds_Start
integer, dimension(1):: Sds_Edges
integer, dimension(1):: Sds_Stride
integer:: Nx

integer:: Sds_Rank
integer:: Sds_Id
integer:: sfcreate
integer:: sfsnatt
integer:: sfwdata
integer:: sfendacc

Istatus = 0
Sds_Rank = 1

Nx = size(Sds_Data)

Sds_Dims = (/Nx/)
Sds_Start = (/0/)
Sds_Stride = (/1/)
Sds_Edges = (/Nx/)

Sds_Id = sfcreate(Sd_Id,Sds_Name,DFNT_FLOAT32,Sds_Rank,Sds_Dims)

if (Sds_Id < 0) then
   Istatus = 1
   return
endif
Istatus = sfsnatt(Sds_Id, "SCALED", DFNT_INT8, 1, 0) + Istatus
Istatus = sfsnatt(Sds_Id, "_Fill_Value", DFNT_FLOAT32, 1, Missing_Value_Real4) + Istatus
Istatus = sfwdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, Sds_Data) + Istatus
Istatus = sfendacc(Sds_Id) + Istatus

end subroutine WRITE_CLAVRX_HDF_SDS_1D

!======================================================================
! routine to write unscaled 2d sds 
!======================================================================
subroutine  WRITE_CLAVRX_HDF_SDS_2D(Sd_Id,Sds_Name,Sds_Data,Istatus)

integer, intent(in):: Sd_Id
character (len=*), intent(in):: Sds_Name
real, dimension(:,:),intent(in):: Sds_Data
integer, intent(out):: Istatus
integer, dimension(2):: Sds_Dims
integer, dimension(2):: Sds_Start
integer, dimension(2):: Sds_Edges
integer, dimension(2):: Sds_Stride
integer:: Nx
integer:: Ny
integer:: Sds_Rank
integer:: Sds_Id
integer:: sfcreate
integer:: sfsnatt
integer:: sfwdata
integer:: sfendacc
integer:: sfdimid
integer:: sfsdmname

Istatus = 0
Sds_Rank = 2

Nx = size(Sds_Data(:,1))
Ny = size(Sds_Data(1,:))

Sds_Dims = (/Nx, Ny/)
Sds_Start = (/0,0/)
Sds_Stride = (/1,1/)
Sds_Edges = (/Nx, Ny/)

Sds_Id = sfcreate(Sd_Id,Sds_Name,DFNT_FLOAT32,Sds_Rank,Sds_Dims)
if (Sds_Id < 0) then
   Istatus = 1
   return
endif
Istatus = sfsnatt(Sds_Id, "SCALED", DFNT_INT8, 1, 0) + Istatus
Istatus = sfsnatt(Sds_Id, "_Fill_Value", DFNT_FLOAT32, 1, Missing_Value_Real4) + Istatus
Istatus = sfsdmname(sfdimid(Sds_Id, 0),"longitude index") + Istatus
Istatus = sfsdmname(sfdimid(Sds_Id, 1),"latitude index") + Istatus
Istatus = sfwdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, Sds_Data) + Istatus
Istatus = sfendacc(Sds_Id) + Istatus

end subroutine WRITE_CLAVRX_HDF_SDS_2D

!======================================================================
! routine to write unscaled 3d sds 
!======================================================================
subroutine  WRITE_CLAVRX_HDF_SDS_3D(Sd_Id,Sds_Name,Sds_Data,Istatus)

integer, intent(in):: Sd_Id
character (len=*), intent(in):: Sds_Name
real, dimension(:,:,:),intent(in):: Sds_Data
integer, intent(out):: Istatus
integer, dimension(3):: Sds_Dims
integer, dimension(3):: Sds_Start
integer, dimension(3):: Sds_Edges
integer, dimension(3):: Sds_Stride
integer:: Nx
integer:: Ny
integer:: Nz
integer:: Sds_Rank
integer:: Sds_Id
integer:: sfcreate
integer:: sfsnatt
integer:: sfwdata
integer:: sfendacc
integer:: sfdimid
integer:: sfsdmname

Istatus = 0
Sds_Rank = 3

Nx = size(Sds_Data(:,1,1))
Ny = size(Sds_Data(1,:,1))
Nz = size(Sds_Data(1,1,:))

Sds_Dims = (/Nx, Ny, Nz/)
Sds_Start = (/0,0,0/)
Sds_Stride = (/1,1,1/)
Sds_Edges = (/Nx, Ny, Nz/)

Sds_Id = sfcreate(Sd_Id,Sds_Name,DFNT_FLOAT32,Sds_Rank,Sds_Dims)
if (Sds_Id < 0) then
   Istatus = 1
   return
endif
Istatus = sfsnatt(Sds_Id, "SCALED", DFNT_INT8, 1, 0) + Istatus
Istatus = sfsnatt(Sds_Id, "_Fill_Value", DFNT_FLOAT32, 1, Missing_Value_Real4) + Istatus
Istatus = sfsdmname(sfdimid(Sds_Id, 0),"longitude index") + Istatus
Istatus = sfsdmname(sfdimid(Sds_Id, 1),"latitude index") + Istatus
Istatus = sfwdata(Sds_Id, Sds_Start, Sds_Stride, Sds_Edges, Sds_Data) + Istatus
Istatus = sfendacc(Sds_Id) + Istatus

end subroutine WRITE_CLAVRX_HDF_SDS_3D


!======================================================================
! opens file
!======================================================================
function OPEN_FILE_HDF_READ(filename,file_id) result(Error_Status)
  character(*), intent(in) :: filename
  integer(kind=int4), intent(out) :: file_id

  integer(kind=int4) :: Error_Status
  integer :: sfstart

  Error_Status = SUCCEED

  ! --- open file
  file_id = sfstart(trim(filename),DFACC_READ)

  if (file_id == FAIL) then
    print "(a,'Cannot Open HDF file, ',a)", trim(filename)
    Error_Status = FAIL
  endif

  return

end function OPEN_FILE_HDF_READ


!======================================================================
! closes file
!======================================================================
function CLOSE_FILE_HDF_READ(file_id,filename) result(Error_Status)
  integer(kind=int4), intent(in) :: file_id
  character(*), intent(in) :: filename

  integer :: Error_Status
  integer :: sfend

  Error_Status = SUCCEED

  ! --- close file
  Error_Status = sfend(file_id)
  if (Error_Status == FAIL) then
    print "(a,'Cannot close HDF file, ',a,' - aborting')", &
           trim(filename)
    stop
  endif

  return

end function CLOSE_FILE_HDF_READ


!======================================================================
! reads dimensions
!======================================================================
function HDF_SDS_DIMENSIONS_READER(id, Sds_Name, Rank, Dims) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), intent(out) :: Rank
  integer(kind=int4), dimension(MAX_RANK_HDF), intent(out) :: Dims

  integer(kind=int4) :: Error_Status
  integer(kind=int4) :: Sds_Id, Istatus, Sds_Type, Sds_Nattr

  Error_Status = SUCCEED

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, Rank, Dims, Sds_Type, Sds_Nattr)

  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",trim(Sds_Name),id
    Error_Status = FAIL
  endif

  Istatus = sfendacc(Sds_Id)
  return

end function HDF_SDS_DIMENSIONS_READER


!-------------------------------------------------------------------
! Subroutine to read 1D int8 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_INT8_1D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(1), intent(inout) :: istart, istride
  integer(kind=int4), dimension(1), intent(inout) :: iedge
  integer(kind=int1), dimension(:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status
  integer(kind=int4), dimension(1) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_INT8 .and. &
      Sds_Type /= DFNT_UINT8 .and. &
      Sds_Type /= DFNT_CHAR8 .and. &
      Sds_Type /= DFNT_CHAR) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 1) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (istart(1) < 0) istart(1) = 0

  if (istride(1) < 1) istride(1) = 1

  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))

  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  iedge(1) = min(max_iedge(1),iedge(1))

  if (iedge(1) .eq. 0) iedge(1) = 1

  if (istart(1) .eq. sds_dims(1)) istart(1) = istart(1) - 1

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_INT8_1D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_INT8_1D


!-------------------------------------------------------------------
! Subroutine to read 1D int16 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_INT16_1D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id                                                                                                                               
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(1), intent(inout) :: istart, istride
  integer(kind=int4), dimension(1), intent(inout) :: iedge
  integer(kind=int2), dimension(:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status
  integer(kind=int4), dimension(1) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_INT16) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 1) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (istart(1) < 0) istart(1) = 0

  if (istride(1) < 1) istride(1) = 1

  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
                                                                                                                                                                        
  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  iedge(1) = min(max_iedge(1),iedge(1))

  if (iedge(1) .eq. 0) iedge(1) = 1

  if (istart(1) .eq. sds_dims(1)) istart(1) = istart(1) - 1

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_INT16_1D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_INT16_1D


!-------------------------------------------------------------------
! Subroutine to read 1D int32 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_INT32_1D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(1), intent(inout) :: istart, istride
  integer(kind=int4), dimension(1), intent(inout) :: iedge
  integer(kind=int4), dimension(:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status
  integer(kind=int4), dimension(1) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_INT32) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 1) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (istart(1) < 0) istart(1) = 0

  if (istride(1) < 1) istride(1) = 1

  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))

  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  iedge(1) = min(max_iedge(1),iedge(1))

  if (iedge(1) .eq. 0) iedge(1) = 1

  if (istart(1) .eq. sds_dims(1)) istart(1) = istart(1) - 1

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then                                                                                                                                            
        print "(a,'Error deallocating 1d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_INT32_1D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_INT32_1D


!-------------------------------------------------------------------
! Subroutine to read 1D float32 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_FLOAT32_1D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(1), intent(inout) :: istart, istride                                                                                                    
  integer(kind=int4), dimension(1), intent(inout) :: iedge
  real(kind=real4), dimension(:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status
  integer(kind=int4), dimension(1) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)

  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_FLOAT32) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 1) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (istart(1) < 0) istart(1) = 0

  if (istride(1) < 1) istride(1) = 1

  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))

  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  iedge(1) = min(max_iedge(1),iedge(1))                                                                                                                                 

  if (iedge(1) .eq. 0) iedge(1) = 1

  if (istart(1) .eq. sds_dims(1)) istart(1) = istart(1) - 1

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_FLOAT32_2D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif
 
  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_FLOAT32_1D


!-------------------------------------------------------------------
! Subroutine to read 1D float64 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_FLOAT64_1D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(1), intent(inout) :: istart, istride
  integer(kind=int4), dimension(1), intent(inout) :: iedge
  real(kind=real8), dimension(:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status
  integer(kind=int4), dimension(1) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_FLOAT64) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 1) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (istart(1) < 0) istart(1) = 0

  if (istride(1) < 1) istride(1) = 1

  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))

  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  iedge(1) = min(max_iedge(1),iedge(1))

  if (iedge(1) .eq. 0) iedge(1) = 1

  if (istart(1) .eq. sds_dims(1)) istart(1) = istart(1) - 1

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_FLOAT64_1D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_FLOAT64_1D


!-------------------------------------------------------------------
! Subroutine to read 2D int8 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_INT8_2D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(2), intent(inout) :: istart
  integer(kind=int4), dimension(2), intent(inout) :: istride
  integer(kind=int4), dimension(2), intent(inout) :: iedge
  integer(kind=int1), dimension(:,:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status
  integer(kind=int4), dimension(2) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_INT8 .and. &
      Sds_Type /= DFNT_UINT8 .and. &
      Sds_Type /= DFNT_CHAR8 .and. &
      Sds_Type /= DFNT_CHAR) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 2) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (istart(1) < 0) istart(1) = 0
  if (istart(2) < 0) istart(2) = 0

  if (istride(1) < 1) istride(1) = 1
  if (istride(2) < 1) istride(2) = 1

  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  if (iedge(2) < 0) iedge(2) = int(ceiling(real(sds_dims(2))/real(istride(2))))

  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  max_iedge(2) = int(ceiling(real(sds_dims(2) - istart(2))/real(istride(2))))
  iedge(1) = min(max_iedge(1),iedge(1))
  iedge(2) = min(max_iedge(2),iedge(2))

  if (iedge(1) .eq. 0) iedge(1) = 1
  if (iedge(2) .eq. 0) iedge(2) = 1

  if (istart(1) .eq. sds_dims(1)) istart(1) = istart(1) - 1
  if (istart(2) .eq. sds_dims(2)) istart(2) = istart(2) - 1

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 2d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_INT8_2D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_INT8_2D


!-------------------------------------------------------------------
! Subroutine to read 2D int16 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_INT16_2D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(2), intent(inout) :: istart, istride
  integer(kind=int4), dimension(2), intent(inout) :: iedge
  integer(kind=int2), dimension(:,:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status
  integer(kind=int4), dimension(2) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_INT16) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 2) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (istart(1) < 0) istart(1) = 0
  if (istart(2) < 0) istart(2) = 0

  if (istride(1) < 1) istride(1) = 1
  if (istride(2) < 1) istride(2) = 1

  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  if (iedge(2) < 0) iedge(2) = int(ceiling(real(sds_dims(2))/real(istride(2))))

  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  max_iedge(2) = int(ceiling(real(sds_dims(2) - istart(2))/real(istride(2))))
  iedge(1) = min(max_iedge(1),iedge(1))
  iedge(2) = min(max_iedge(2),iedge(2))

  if (iedge(1) .eq. 0) iedge(1) = 1
  if (iedge(2) .eq. 0) iedge(2) = 1

  if (istart(1) .eq. sds_dims(1)) istart(1) = istart(1) - 1
  if (istart(2) .eq. sds_dims(2)) istart(2) = istart(2) - 1

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 2d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_INT16_2D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_INT16_2D

!-------------------------------------------------------------------
! Subroutine to read 2D int32 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_INT32_2D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(2), intent(inout) :: istart, istride
  integer(kind=int4), dimension(2), intent(inout) :: iedge
  integer(kind=int4), dimension(:,:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status
  integer(kind=int4), dimension(2) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_INT32) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 2) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (istart(1) < 0) istart(1) = 0
  if (istart(2) < 0) istart(2) = 0

  if (istride(1) < 1) istride(1) = 1
  if (istride(2) < 1) istride(2) = 1

  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  if (iedge(2) < 0) iedge(2) = int(ceiling(real(sds_dims(2))/real(istride(2))))

  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  max_iedge(2) = int(ceiling(real(sds_dims(2) - istart(2))/real(istride(2))))
  iedge(1) = min(max_iedge(1),iedge(1))
  iedge(2) = min(max_iedge(2),iedge(2))

  if (iedge(1) .eq. 0) iedge(1) = 1
  if (iedge(2) .eq. 0) iedge(2) = 1

  if (istart(1) .eq. sds_dims(1)) istart(1) = istart(1) - 1
  if (istart(2) .eq. sds_dims(2)) istart(2) = istart(2) - 1

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 2d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_INT32_2D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_INT32_2D


!-------------------------------------------------------------------
! Subroutine to read 2D float32 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_FLOAT32_2D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(2), intent(inout) :: istart, istride
  integer(kind=int4), dimension(2), intent(inout) :: iedge
  real(kind=real4), dimension(:,:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status
  integer(kind=int4), dimension(2) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_FLOAT32) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 2) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (istart(1) < 0) istart(1) = 0
  if (istart(2) < 0) istart(2) = 0

  if (istride(1) < 1) istride(1) = 1
  if (istride(2) < 1) istride(2) = 1

  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  if (iedge(2) < 0) iedge(2) = int(ceiling(real(sds_dims(2))/real(istride(2))))

  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  max_iedge(2) = int(ceiling(real(sds_dims(2) - istart(2))/real(istride(2))))
  iedge(1) = min(max_iedge(1),iedge(1))
  iedge(2) = min(max_iedge(2),iedge(2))

  if (iedge(1) .eq. 0) iedge(1) = 1
  if (iedge(2) .eq. 0) iedge(2) = 1

  if (istart(1) .eq. sds_dims(1)) istart(1) = istart(1) - 1
  if (istart(2) .eq. sds_dims(2)) istart(2) = istart(2) - 1

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 2d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_FLOAT32_2D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_FLOAT32_2D


!-------------------------------------------------------------------
! Subroutine to read 2D float64 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_FLOAT64_2D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(2), intent(inout) :: istart, istride
  integer(kind=int4), dimension(2), intent(inout) :: iedge
  real(kind=real8), dimension(:,:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status
  integer(kind=int4), dimension(2) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_FLOAT64) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 2) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (istart(1) < 0) istart(1) = 0
  if (istart(2) < 0) istart(2) = 0

  if (istride(1) < 1) istride(1) = 1
  if (istride(2) < 1) istride(2) = 1

  if (iedge(1) < 0) iedge(1) = int(ceiling(real(sds_dims(1))/real(istride(1))))
  if (iedge(2) < 0) iedge(2) = int(ceiling(real(sds_dims(2))/real(istride(2))))

  max_iedge(1) = int(ceiling(real(sds_dims(1) - istart(1))/real(istride(1))))
  max_iedge(2) = int(ceiling(real(sds_dims(2) - istart(2))/real(istride(2))))
  iedge(1) = min(max_iedge(1),iedge(1))
  iedge(2) = min(max_iedge(2),iedge(2))

  if (iedge(1) .eq. 0) iedge(1) = 1
  if (iedge(2) .eq. 0) iedge(2) = 1

  if (istart(1) .eq. sds_dims(1)) istart(1) = istart(1) - 1
  if (istart(2) .eq. sds_dims(2)) istart(2) = istart(2) - 1

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 2d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 2d buffer.')",EXE_PROMPT
      stop
    endif
  endif
 
  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_FLOAT64_2D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_FLOAT64_2D


!-------------------------------------------------------------------
! Subroutine to read 3D int8 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_INT8_3D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(3), intent(inout) :: istart, istride
  integer(kind=int4), dimension(3), intent(inout) :: iedge
  integer(kind=int1), dimension(:,:,:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status, idim
  integer(kind=int4), dimension(3) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_INT8 .and. &
      Sds_Type /= DFNT_UINT8 .and. &
      Sds_Type /= DFNT_CHAR8 .and. &
      Sds_Type /= DFNT_CHAR) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 3) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  do idim=1, sds_rank

    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1

    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))

  end do

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. & 
        size(buffer,3) < iedge(3)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 3d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 3d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_INT8_3D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_INT8_3D


!-------------------------------------------------------------------
! Subroutine to read 3D int16 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_INT16_3D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(3), intent(inout) :: istart, istride
  integer(kind=int4), dimension(3), intent(inout) :: iedge
  integer(kind=int2), dimension(:,:,:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status, idim
  integer(kind=int4), dimension(3) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_INT16) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 3) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  do idim=1, sds_rank

    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1

    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))

  end do

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
       size(buffer,3) < iedge(3)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 3d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 3d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_INT16_3D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_INT16_3D


!-------------------------------------------------------------------
! Subroutine to read 3D int32 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_INT32_3D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(3), intent(inout) :: istart, istride
  integer(kind=int4), dimension(3), intent(inout) :: iedge
  integer(kind=int4), dimension(:,:,:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status, idim
  integer(kind=int4), dimension(3) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_INT32) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 3) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  do idim=1, sds_rank

    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1

    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))

  end do

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
        size(buffer,3) < iedge(3)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 3d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 3d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_INT32_2D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_INT32_3D


!-------------------------------------------------------------------
! Subroutine to read 3D float32 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_FLOAT32_3D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(3), intent(inout) :: istart, istride
  integer(kind=int4), dimension(3), intent(inout) :: iedge
  real(kind=real4), dimension(:,:,:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status, idim
  integer(kind=int4), dimension(3) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_FLOAT32) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 3) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  do idim=1, sds_rank

    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1

    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))

  end do

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. & 
        size(buffer,3) < iedge(3)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 3d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 3d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_FLOAT32_3D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_FLOAT32_3D


!-------------------------------------------------------------------
! Subroutine to read 3D float64 hdf data.
!-------------------------------------------------------------------

function READ_HDF_SDS_FLOAT64_3D(Sd_Id, Sds_Name, istart, istride, iedge, buffer, type) result(Error_Status)
  integer(kind=int4), intent(in) :: Sd_Id
  character(*), intent(in) :: Sds_Name
  character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int4), dimension(3), intent(inout) :: istart, istride
  integer(kind=int4), dimension(3), intent(inout) :: iedge
  real(kind=real8), dimension(:,:,:), allocatable, intent(inout) :: buffer
  integer(kind=int4), intent(out), optional :: type
  integer(kind=int4) :: astatus, Istatus, Sds_Id, sds_rank, Sds_Type, Sds_Nattr, Error_Status, idim
  integer(kind=int4), dimension(3) :: sds_dims, max_iedge
  logical:: size_check

  Error_Status = SUCCEED

  Sds_Id = sfselect(Sd_Id, sfn2index(Sd_Id,trim(Sds_Name)))
  if (Sds_Id == FAIL) then
    print "(a,'Error selecting ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfginfo(Sds_Id, Sds_Name_Temp, sds_rank, sds_dims, Sds_Type, Sds_Nattr)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name_Temp),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_FLOAT64) then
    print "(a,'Error reading (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (sds_rank /= 3) then
    print "(a,'Error reading (rank mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  do idim=1, sds_rank

    if (istart(idim) < 0) istart(idim) = 0
    if (istride(idim) < 1) istride(idim) = 1

    if (iedge(idim) < 0) iedge(idim) = int(ceiling(real(sds_dims(idim))/real(istride(idim))))
    max_iedge(idim) = int(ceiling(real(sds_dims(idim) - istart(idim))/real(istride(idim))))
    iedge(idim) = min(max_iedge(idim),iedge(idim))

  end do

  if (allocated(buffer)) then
    if (size(buffer,1) < iedge(1) .or. &
        size(buffer,2) < iedge(2) .or. &
        size(buffer,3) < iedge(3)) then
      deallocate(buffer,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 3d HDF buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(buffer))) then
    allocate(buffer(iedge(1),iedge(2),iedge(3)),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 3d buffer.')",EXE_PROMPT
      stop
    endif
  endif

  !--- check for conistency in size of input and output
  size_check = CHECK_EDGE_VS_OUTPUT_SHAPE(iedge,shape(buffer))

  !-- if not consistent size, return and don't attempt to read
  if (.not. size_check) then
    print *, "Size inconsistency in READ_HDF_SDS_FLOAT64_3D, skipping read"
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return   
  endif

  Istatus = sfrdata(Sds_Id, istart, istride, iedge, buffer)
  if (Istatus /= 0) then
    print "(a,'Error reading ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Sds_Name),Sd_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_SDS_FLOAT64_3D


!-------------------------------------------------------------------
! This routine is used to read char8 HDF SDS attributes.
!-------------------------------------------------------------------

function READ_HDF_ATTRIBUTE_CHAR8_SCALAR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  !character(len=len(Sds_Name)):: Sds_Name_Temp
  character(len=*), intent(inout) :: attr
  character(len=1000), dimension(1) :: buffer
  integer(kind=int4), intent(out), optional :: type

  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status
  integer :: sffattr, sfrcatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id

  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif

  if (Sds_Type /= DFNT_INT8 .and. &
      Sds_Type /= DFNT_CHAR8 .and. &
      Sds_Type /= DFNT_CHAR) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  Istatus = sfrcatt(Sds_Id, Attr_Index, buffer)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    attr = " "
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  attr = TRIM(buffer(1)(1:count))

  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_ATTRIBUTE_CHAR8_SCALAR


!-------------------------------------------------------------------
! This routine is used to read int8 HDF SDS attributes.
!-------------------------------------------------------------------

function READ_HDF_ATTRIBUTE_INT8_SCALAR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  !character(len=len(Sds_Name)):: Sds_Name_Temp
  integer(kind=int1), intent(inout) :: attr
  integer(kind=int4), intent(out), optional :: type
 
  integer(kind=int1), dimension(1) :: buffer
  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status
  integer :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED
  attr = MISSING_VALUE_INT1

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id
 
  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Error_Status = FAIL                                                                                                                                          
    return
  endif
 
  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif
 
  if (Sds_Type /= DFNT_INT8 .and. &
      Sds_Type /= DFNT_CHAR8 .and. &
      Sds_Type /= DFNT_CHAR) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
 
  if (count /= 1) then
    print "(a,'Error reading attribute (count mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfrnatt(Sds_Id, Attr_Index, buffer)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  Istatus = sfendacc(Sds_Id)

  attr = buffer(1)

  return

end function READ_HDF_ATTRIBUTE_INT8_SCALAR


!-------------------------------------------------------------------
! This routine is used to read int16 HDF SDS attributes.
!-------------------------------------------------------------------
 
function READ_HDF_ATTRIBUTE_INT16_SCALAR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  integer(kind=int2), intent(inout) :: attr
  integer(kind=int4), intent(out), optional :: type
 
  integer(kind=int2), dimension(1) :: buffer
  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status
  integer :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED
  attr = MISSING_VALUE_INT2

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id
 
  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif
 
  if (Sds_Type /= DFNT_INT16) then
    print*,trim(Sds_Name),type,DFNT_INT16
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
 
  if (count /= 1) then
    print "(a,'Error reading attribute (count mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfrnatt(Sds_Id, Attr_Index, buffer)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  Istatus = sfendacc(Sds_Id)

  attr = buffer(1)

  return

end function READ_HDF_ATTRIBUTE_INT16_SCALAR


!-------------------------------------------------------------------
! This routine is used to read int32 HDF SDS attributes.
!-------------------------------------------------------------------

function READ_HDF_ATTRIBUTE_INT32_SCALAR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  integer(kind=int4), intent(inout) :: attr
  integer(kind=int4), intent(out), optional :: type
 
  integer(kind=int4), dimension(1) :: buffer
  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status
  integer :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED
  attr = MISSING_VALUE_INT4

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id
 
  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL                                                                                                                                          
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif
 
  if (Sds_Type /= DFNT_INT32) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
 
  if (count /= 1) then
    print "(a,'Error reading attribute (count mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfrnatt(Sds_Id, Attr_Index, buffer)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  Istatus = sfendacc(Sds_Id)

  attr = buffer(1)

  return

end function READ_HDF_ATTRIBUTE_INT32_SCALAR


!-------------------------------------------------------------------
! This routine is used to read float32 HDF SDS attributes.
!-------------------------------------------------------------------
 
function READ_HDF_ATTRIBUTE_FLOAT32_SCALAR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  real(kind=real4), intent(inout) :: attr
  integer(kind=int4), intent(out), optional :: type
 
  real(kind=real4), dimension(1) :: buffer
  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status
  integer :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED
  attr = MISSING_VALUE_REAL4

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id
 
  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif
 
  if (Sds_Type /= DFNT_FLOAT32) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
  
  if (count /= 1) then
    print "(a,'Error reading attribute (count mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfrnatt(Sds_Id, Attr_Index, buffer)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  Istatus = sfendacc(Sds_Id)

  attr = buffer(1)

  return

end function READ_HDF_ATTRIBUTE_FLOAT32_SCALAR


!-------------------------------------------------------------------
! This routine is used to read float64 HDF SDS attributes.
!-------------------------------------------------------------------

function READ_HDF_ATTRIBUTE_FLOAT64_SCALAR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  real(kind=real8), intent(inout) :: attr
  integer(kind=int4), intent(out), optional :: type
 
  real(kind=real8), dimension(1) :: buffer
  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status
  integer :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED
  attr = MISSING_VALUE_REAL8

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id
 
  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Error_Status = FAIL
    return
  endif

  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif                                                                                                                                                                 

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif
 
  if (Sds_Type /= DFNT_FLOAT64) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
 
  if (count /= 1) then
    print "(a,'Error reading attribute (count mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfrnatt(Sds_Id, Attr_Index, buffer)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  Istatus = sfendacc(Sds_Id)

  attr = buffer(1)

  return

end function READ_HDF_ATTRIBUTE_FLOAT64_SCALAR


!-------------------------------------------------------------------
! This routine is used to read int8 HDF SDS attributes.
!-------------------------------------------------------------------
 
function READ_HDF_ATTRIBUTE_INT8_VECTOR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  integer(kind=int1), dimension(:), allocatable, intent(inout) :: attr
  integer(kind=int4), intent(out), optional :: type

  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status, astatus
  integer :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id
 
  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif
 
  if (Sds_Type /= DFNT_INT8 .and. &
      Sds_Type /= DFNT_UINT8 .and. &
      Sds_Type /= DFNT_CHAR8 .and. &
      Sds_Type /= DFNT_CHAR) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (allocated(attr)) then
    if (size(attr,1) < count) then
      deallocate(attr,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF attr buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(attr))) then
    allocate(attr(count),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d attr buffer.')",EXE_PROMPT
      stop
    endif
  endif
 
  Istatus = sfrnatt(Sds_Id, Attr_Index, attr)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_ATTRIBUTE_INT8_VECTOR


!-------------------------------------------------------------------
! This routine is used to read int16 HDF SDS attributes.
!-------------------------------------------------------------------
 
function READ_HDF_ATTRIBUTE_INT16_VECTOR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  integer(kind=int2), dimension(:), allocatable, intent(inout) :: attr
  integer(kind=int4), intent(out), optional :: type 
 
  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status, astatus
  integer :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED

  ! attr = MISSING_VALUE_REAL4

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id 
 
  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id                                                                          
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif
 
  if (Sds_Type /= DFNT_INT16 .AND. &
      Sds_Type /= DFNT_UINT16) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (allocated(attr)) then
    if (size(attr,1) < count) then
      deallocate(attr,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF attr buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(attr))) then
    allocate(attr(count),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d attr buffer.')",EXE_PROMPT
      stop
    endif
  endif

  attr = MISSING_VALUE_INT2
 
  Istatus = sfrnatt(Sds_Id, Attr_Index, attr)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif                                                                                                                                                                 

  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_ATTRIBUTE_INT16_VECTOR


!-------------------------------------------------------------------
! This routine is used to read int32 HDF SDS attributes.
!-------------------------------------------------------------------
 
function READ_HDF_ATTRIBUTE_INT32_VECTOR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  integer(kind=int4), dimension(:), allocatable, intent(inout) :: attr
  integer(kind=int4), intent(out), optional :: type
 
  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status, astatus
  integer :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED

  ! attr = MISSING_VALUE_REAL4

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id
 
  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id                                                                          
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif 
 
  if (Sds_Type /= DFNT_INT32 .AND. &
      Sds_Type /= DFNT_UINT32) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (allocated(attr)) then
    if (size(attr,1) < count) then
      deallocate(attr,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF attr buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(attr))) then
    allocate(attr(count),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d attr buffer.')",EXE_PROMPT
      stop
    endif
  endif

  attr = MISSING_VALUE_INT4
 
  Istatus = sfrnatt(Sds_Id, Attr_Index, attr)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif                                                                                                                                                                 

  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_ATTRIBUTE_INT32_VECTOR


!-------------------------------------------------------------------
! This routine is used to read float32 HDF SDS attributes.
!-------------------------------------------------------------------
 
function READ_HDF_ATTRIBUTE_FLOAT32_VECTOR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  real(kind=real4), dimension(:), allocatable, intent(inout) :: attr
  integer(kind=int4), intent(out), optional :: type
 
  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status, astatus
  integer :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED
  !attr = MISSING_VALUE_REAL4

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id
 
  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Error_Status = FAIL
    return
  endif
 
  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif
 
  if (Sds_Type /= DFNT_FLOAT32) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (allocated(attr)) then
    if (size(attr,1) < count) then
      deallocate(attr,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF attr buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(attr))) then
    allocate(attr(count),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d attr buffer.')",EXE_PROMPT
      stop
    endif
  endif

  attr = MISSING_VALUE_REAL4
 
  Istatus = sfrnatt(Sds_Id, Attr_Index, attr)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  Istatus = sfendacc(Sds_Id)

  return

end function READ_HDF_ATTRIBUTE_FLOAT32_VECTOR


!-------------------------------------------------------------------
! This routine is used to read float32 HDF SDS attributes.
!-------------------------------------------------------------------

function READ_HDF_ATTRIBUTE_FLOAT64_VECTOR(id, Sds_Name, Attr_Name, attr, type) result(Error_Status)
  integer(kind=int4), intent(in) :: id
  character(len=*), intent(in) :: Sds_Name, Attr_Name
  real(kind=real8), dimension(:), allocatable, intent(inout) :: attr
  integer(kind=int4), intent(out), optional :: type
 
  character(len=1020) :: name
  integer(kind=int4) :: Istatus, Attr_Index, Sds_Type, count, Sds_Id, Error_Status, astatus
  integer :: sffattr, sfrnatt, sfgainfo, sfn2index, sfselect, sfendacc

  Error_Status = SUCCEED

  Sds_Id = sfselect(id, sfn2index(id,trim(Sds_Name)))
  if (Sds_Id == FAIL) Sds_Id = id

  Attr_Index = sffattr(Sds_Id, trim(Attr_Name))
  if (Attr_Index == FAIL) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Error_Status = FAIL
    return                                                                                                                                                              
  endif
 
  Istatus = sfgainfo(Sds_Id, Attr_Index, name, Sds_Type, count)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (present(type)) then
    type = Sds_Type
    Istatus = sfendacc(Sds_Id)
    Error_Status = SUCCEED
    return
  endif
 
  if (Sds_Type /= DFNT_FLOAT64) then
    print "(a,'Error reading attribute (type mismatch) ',a,' from Sd_Id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  if (allocated(attr)) then
    if (size(attr,1) < count) then
      deallocate(attr,stat=astatus)
      if (astatus /= 0) then
        print "(a,'Error deallocating 1d HDF attr buffer.')",EXE_PROMPT
        stop
      endif
    endif
  endif

  if ((.not. allocated(attr))) then
    allocate(attr(count),stat=astatus)
    if (astatus /= 0) then
      print "(a,'Not enough memory to allocate 1d attr buffer.')",EXE_PROMPT
      stop
    endif
  endif

  attr = MISSING_VALUE_REAL8
 
  Istatus = sfrnatt(Sds_Id, Attr_Index, attr)
  if (Istatus /= 0) then
    print "(a,'Attribute ',a,' reading error from id: ',i0)",EXE_PROMPT,trim(Attr_Name),Sds_Id
    Istatus = sfendacc(Sds_Id)
    Error_Status = FAIL
    return
  endif

  Istatus = sfendacc(Sds_Id)
                                                                                                                                                                        
  return

end function READ_HDF_ATTRIBUTE_FLOAT64_VECTOR
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
function CHECK_EDGE_VS_OUTPUT_SHAPE(size_in, size_out) result(size_check)
   integer(kind=int4), dimension(:), intent(in):: size_in
   integer(kind=int4), dimension(:), intent(in):: size_out
   logical:: size_check
   integer:: i, n

   n = size(size_in)

   size_check = .true.

   do i = 1,n
     if (size_in(i) /= size_out(i)) size_check = .false.
   enddo      

   return
end function CHECK_EDGE_VS_OUTPUT_SHAPE


!======================================================================

end module CX_HDF4_MOD

