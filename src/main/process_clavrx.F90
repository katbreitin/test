! $Id: process_clavrx.f90 4129 2021-04-19 19:30:41Z heidinger $
  program PROCESS_CLAVRX
!-----------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: process_clavrx.f90 (src)
!       PROCESS_CLAVRX (program)
!       CLAVRXORB (executable)
!
! PURPOSE:
!
! DESCRIPTION: This code serves as the NESDIS operational cloud
!      processing system (CLAVR-x). This code also serves as the cloud
!      climate data generation system (PATMOS-x)
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
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
! This copyright pertains to all routines in the CLAVR-x system unless stated
!
! REVISION HISTORY:
! Version 5.0 - 2010 GEWEX submission
! Version 5.2 - Code Delivered to NCDC
! Version 5.2.1 - MODIS Capability
! Version 5.2.2 - GOES imager Capability, GlobSnow and CFSR capability
! Version 5.2.3 - Istvan Laszlo's Insolation (SASRAB) added
! Version 5.2.4 - Volcanic Ash added
! Version 5.2.5 - VIIRS M-band Added
! Version 5.2.6 - VIIRS DNB Support
! Version 5.2.7 - MTSAT, COMS Support
! Version 5.2.8 - GOES Sounder Support
! Version 5.3.0 - NCDC Delivery
! Version 2015a - Stable Release for CSPP
! Version 6.0 - NCEI code, AVHRR+HIRS Support
!
!
! Basic Running Instruction
!  The input to this code is controlled through three mechanisms
!    1. command-line options  (type clavrxorb --help to see documentation)
!    2. a FILELIST - a list of level-1b files and directories (default name is
!      clavrxorb_File_list)
!    3. a OPTIONSLIST - a list of processing options (default is
!                     clavrxorb_Default_Options)
!
! Overview of capabilities.
!    CLAVRXORB can
!       - use AVHRR level-1b calibration or apply new calibration routines
!       - use AVHRR level-1b geolocation or apply new geolocation routines
!       - process NESDIS or AAPP AVHRR Level1b
!       - process MYD021KM or MYD02SSH MODIS Level1b files
!       - process band-separated AREA files from GOES Imager/Sounder,
!         SEVIRI, MTSAT-1R, MTSAT-2, COMS, FY-2
!       - generate pixel level cloud, aerosol and surface products
!       - write to a series of pixel-level hdf files
!       - write a level-3 file (gridded data for each orbit - AVHRR only)
!
! In general, CLAVRXORB uses global data arrays and structures to pass data
!
! Note, comments the begin with "Marker" refer to flowchart delivered to NCDC
!
! Web-page:  http://cimss.ssec.wisc.edu/clavr or
!            http://cimss.ssec.wisc.edu/patmosx
!
! Channels 1 - 36 refer to MODIS or their analogs on other sensors
! Channel 37 = ABI Channel 8 = 6.2 micron
! Channel 38 = ABI Channel 13 = 10.4 micron
! Channels 39-44 are defined only for VIIRS
! Channel 39 - VIIRS I1 - 0.64 micron
! Channel 40 - VIIRS I2 - 0.865 micron
! Channel 41 - VIIRS I3 - 1.61 micron
! Channel 42 - VIIRS I4 - 3.74 micron
! Channel 43 - VIIRS I5 - 11.45 micron
! Channel 44 - VIIRS DNB - 0.7 micron
!
!-------------------------------------------------------------------------

!*****************************************************************************
! Marker: ACCESS MODULES
!******************************************************************************
   use ACHA_CLAVRX_BRIDGE,only: &
      awg_cloud_height_bridge

   use AEROSOL_PROPERTIES, only: &
      pixel_aer_ret_ocean &
      , read_aer_ch123a_ref_luts

   use ASOS_CLAVRX_BRIDGE, only: ASOS_BRIDGE

   use AVHRR_REPOSITION_MOD,only: &
      interpolate_clock_error &
      , reposition_for_clock_error &
      , setup_clock_corrections

   use BASELINE_CLOUD_HEIGHT,only: &
      baseline_cloud_height_main &
      , populate_planck_tables &
      , populate_planck_tables

   use BASELINE_CLOUD_MASK, only: BASELINE_CLOUD_MASK_MAIN

   use CALIBRATION_CONSTANTS_MOD, only: &
      Sun_Earth_Distance &
      , launch_date &
      , c1, c2, planck_a1 , planck_a2, planck_nu &
      , Ch1_Gain_Low,Ch1_gain_High &
      , Ch2_Gain_Low,Ch2_gain_High &
      , Ch1_Switch_Count_Cal,Ch1_Dark_Count_Cal &
      , Ch2_Switch_Count_Cal,Ch2_Dark_Count_Cal &
      , Ch3a_Gain_low,Ch3a_Gain_High &
      , Ch3a_Switch_Count_Cal,Ch3a_Dark_Count_Cal &
      , Solar_Ch20_Nu &
      , MITIGATE_CH20_NOISE

   use CCL_CLAVRX_BRIDGE, only: CCL_BRIDGE

   use CLAVRX_MESSAGE_MOD, only: MESG, VERB_LEV

   use CLAVRX_OLR_MOD, only: COMPUTE_OLR, SETUP_OLR

   use CLAVRX_SST_MOD, only: SETUP_SST, COMPUTE_SST

   use CLOUD_BASE_CLAVRX_BRIDGE,only: CLOUD_BASE_BRIDGE

   use CLOUD_HEIGHT_ROUTINES, only: &
      co2_slicing_cloud_height &
      , co2irw_cloud_height &
      , compute_altitude_from_pressure &
      , compute_cloud_top_level_nwp_wind_and_tpw &
      , compute_csbt_cloud_masks &
      , convective_cloud_probability &
      , H2O_INTERCEPT_CLOUD_HEIGHT &
      , MAKE_CIRRUS_PRIOR_TEMPERATURE_FROM_CO2 &
      , mode_zero_cloud_height &
      , modify_cloud_type_with_sounder &
      , opaque_cloud_height &
      , opaque_transmission_height &
      , sounder_emissivity &
      , splitwin_cloud_height &
      , supercooled_cloud_probability

   use CLOUD_TYPE_BRIDGE_MODULE,only: &
    cloud_type_bridge &
    , set_cloud_type_version

   use CX_SPATIAL_METRICS_MOD, only: &
        COMPUTE_MEDIAN_METRICS_L1B &
      , COMPUTE_MEDIAN_METRICS_L2 &
      , COMPUTE_MIN_MAX_MEAN_STD_METRICS &
      , COMPUTE_SPATIAL_CORRELATION_ARRAYS &
      , COMPUTE_RADIATIVE_CENTER_ARRAYS

   use univ_fp_comparison_mod, only: operator(.EQfp.), operator(.NEfp.)

   use CONSTANTS_MOD

   use CX_DATE_TIME_TOOLS_MOD, only: &
      LEAP_YEAR_FCT &
      , COMPUTE_MONTH &
      , COMPUTE_DAY &
      , COMPUTE_TIME_HOURS

   use CX_MURI_CLAVRX_BRIDGE_MOD, only: &
        MURI &
      , CX_MURI_ALGORITHM

   use CX_SEA_IR_EMISS_MOD, only: &
    forget_sea_ir_emiss

   use CX_DUST_MOD, only: &
        FORGET_ABI_DUST &
      , GET_SEGMENT_ABI_DUST_PROB &
      , READ_ABI_DUST_LUT

   use CX_NUCAPS_MOD, only: &
        VIIRS_NUCAPS &
      , CONVERT_SMOOTH_NUCAPS_TEMP

   use DCOMP_DERIVED_PRODUCTS_MOD,only: &
        ADJUST_DCOMP_LWP &
      , COMPUTE_ADIABATIC_CLOUD_PROPS &
      , COMPUTE_CLOUD_WATER_PATH &
      , COMPUTE_DCOMP_INSOLATION &
      , COMPUTE_PRECIPITATION &
      , COMPUTE_PRECIPITATION_AHI &
      , COMPUTE_MASS_CONCENTRATION &
      , COMPUTE_SUBPIXEL_MAX_MIN_COD

   use DNB_RETRIEVALS_MOD, only: &
        COMPUTE_LUNAR_REFLECTANCE &
      , lunar_reflectance_nasa_data_adjustment

   use DNCOMP_CLAVRX_BRIDGE_MOD, only: &
      AWG_CLOUD_DNCOMP_ALGORITHM


   use FILE_UTILS, only: FILE_TEST, GET_LUN

   use GFS_HDF_MOD, only: READ_GFS_DATA

   use GLOBSNOW_READ_ROUTINES,only: &
      GET_GLOBSNOW_FILENAME &
     , GET_PIXEL_GLOBSNOW_ANALYSIS &
     , READ_GLOBSNOW_ANALYSIS_MAP

   use GOES_MOD,only: &
      Area_Struct &
      , GVAR_NAV &
      , DARK_COMPOSITE_CLOUD_MASK &
      , DETERMINE_DARK_COMPOSITE_NAME &
      , POST_PROCESS_GOES_DARK_COMPOSITE

   use HIRS_FUSION_MOD, only: &
      SET_REPLACED_AVHRR_TO_MISSING &
    , HIRS_AVHRR_FUSION_PREPERATION

   use IR_CLOUD_TYPE_BAUM_MODULE, only: &
       IR_CLOUD_TYPE_BAUM &
       , SET_IR_CLOUD_TYPE_VERSION

   use LAND_SFC_PROPERTIES_MOD, only: &
       Land_grid_Description &
    , OPEN_LAND_SFC_HDF &
    ,  GET_SNOW_MAP_FILENAME &
    , close_land_sfc_hdf &
    , read_land_sfc_hdf

    use LASZLO_INSOLATION_MOD,only: &
      insolation

   use LEVEL2_MOD,only: &
        CLOSE_PIXEL_HDF_FILES &
      , DEFINE_HDF_FILE_STRUCTURES &
      , WRITE_ALGORITHM_ATTRIBUTES &
      , WRITE_SEGMENT_LEVEL2 &
      , READ_LEVEL2_VAR_LIST &
      , SETUP_LEVEL2_SDS_INFO

   use NB_CLOUD_MASK_CLAVRX_BRIDGE, only: NB_CLOUD_MASK_BRIDGE

   use ECM2_CLOUD_MASK_CLAVRX_BRIDGE, only: ECM2_CLOUD_MASK_BRIDGE

   use NCEP_REANALYSIS,only: READ_NCEP_REANALYSIS_DATA

   use NWP_COMMON_MOD, only: &
       NWP &
       , compute_coast_mask_nwp &
       , compute_nwp_levels_segment &
       , compute_pixel_nwp_parameters &
       , compute_segment_nwp_cloud_parameters &
       , destroy_nwp_arrays &
       , map_pixel_nwp &
       , modify_tsfc_nwp_pix &
       , qc_nwp &
       , temporal_interp_tmpsfc_nwp

   use OCA_MOD, only: &
       Read_OCA

   use OISST_ANALYSIS,only: &
        Get_OISST_Map_Filename &
      , Get_pixel_sst_analysis &
      , Read_oisst_analysis_map

   use PIXEL_COMMON_MOD,only: &
      ! - major global structures
      Ch &
    , Geo &
    , Nav &
    , Sfc &
    , Image &
    , Sensor &
    , ACHA &
    , BASE &
    , ASOS &
    , NUCAPS &
    , CLDMASK &
    , NWP_PIX &
      ! - routines
    , Destroy_Pixel_Arrays &
    , Create_Pixel_Arrays &
    , Reset_Pixel_Arrays_To_Missing  &
    , initial_pixel_common_alloc &
    , nav_opt &
    , file_list &
    , dir_level2 &
    , Temporary_Data_Dir &
    , Use_Sst_Anal &
    , Use_Sea_IR_Emiss &
    , Use_Land_IR_Emiss &
    , Use_ABI_Dust &
    , OiSst_Data_Dir &
    , Gdas_Data_Dir &
    , Cfsr_Data_Dir &
    , Ncep_Data_Dir &
    , Gfs_Data_Dir &
    , L1b_Gzip,L1b_Bzip2 &
    , Line_Idx_Min_Segment &
    , Line_Idx_Max_Segment  &
    , dark_composite_name &
    , read_dark_comp &
    , Ancil_data_dir &
    , Modis_Clr_Alb_Flag &
    , number_of_temporary_files &
    , Merra_Data_Dir &
    , Ref_Ch1_Dark_Composite &
    , Bad_Pixel_Mask &
    , Bt_Ch20_Median_5x5 &
    , lrc_flag &
    , SOLAR_CONTAMINATION_MASK &
    , Ref_Cal_1b &
    , Level2_File_Flag &
    , AVHRR_Fusion_Flag &
    , Pc_Co2,Tc_Co2,Zc_Co2, Ec_Co2 &
    , Pc_H2O,Tc_H2O,Zc_H2O  &
    , Nonconfident_Cloud_Mask_Count &
    , cloud_mask_count &
    , Cloud_Mask_Bayesian_Flag &
    , Cloud_Mask_Mode &
    , Cld_Type_Aux &
    , cld_type &
    , Cloud_Type_Aux_Read_Flag &
    , Cloud_Mask_Aux_Read_Flag &
    , cld_flag &
    , Skip_Output &
    , Cld_Phase_IR &
    , Use_IR_Cloud_Type_Flag &
    , Cld_Phase &
    , Cld_Temp_Sounder &
    , Cld_Emiss_Sounder &
    , Tc_Cirrus_Background &
    , Zc_Cirrus_Background &
    , Cloud_Fraction_Background &
    , cld_press_sounder &
    , DCOMP_Mode, NLCOMP_Mode &
    , Aerosol_Mode &
    , DCOMP_Processed_Count &
    , DCOMP_Valid_Count &
    , sasrab_flag &
    , Temporary_File_Name &
    , zen_idx_rtm &
    , cld_type_IR &
    , Cld_Type_Aux &
    , Cld_Phase_Aux &
    , Pc_Top1_Aux &
    , Pc_Top2_Aux &
    , Pc_Uncertainty1_Aux &
    , Pc_Uncertainty2_Aux &
    , Cost_Aux &
    , Tau_Aux &
    , sst_anal_cice &
    , Read_Volcano_Mask &
    , Read_Land_Mask &
    , Read_Coast_Mask &
    , Read_Surface_Elevation &
    , Read_Hires_Sfc_Type &
    , Read_Snow_Mask &
    , Read_Dark_Comp &
    , two_byte_temp &
    , GlobSnow_Data_Dir &
    , Failed_Glob_Snow_Mask_Flag &
    , snow_data_dir &
    , FAILED_IMS_SNOW_MASK_FLAG &
    , Orbital_Processing_Time_Minutes &
    , Therm_Cal_1b &
    , Num_Scans_Level2_Hdf &
    , Use_Aux_Flag &
    , Goes_Scan_Line_Flag &
    , Erai_Data_Dir &
    , Skip_L1b_File_Flag &
    , Nucaps_Flag &
    , Use_Sst_Anal &
    , Sst_Anal &
    , ABI_Use_104um_Flag &
    , Static_Ref_065um_Dark_Composite &
    , Subpixel_Cloud_Fraction &
    , Static_Dark_Sky_Flag &
    , Verbose_Level_Flag

   use SNOW_ROUTINES_MOD, only: &
      COMPUTE_SNOW_CLASS &
      , COMPUTE_SNOW_CLASS_NWP &
      , COMPUTE_SNOW_CLASS_OISST

   use PIXEL_ROUTINES_MOD,only: &
      SET_CHAN_ON_FLAG &
      , DESERT_MASK_FOR_CLOUD_DETECTION &
      , CITY_MASK_FOR_CLOUD_DETECTION  &
      , ADJACENT_PIXEL_CLOUD_MASK &
      , ASSIGN_CLEAR_SKY_QUALITY_FLAGS &
      , CH20_PSEUDO_REFLECTANCE &
      , Compute_cloud_mask_performance_metrics &
      , Compute_acha_Performance_metrics &
      , Compute_DCOMP_Performance_metrics &
      , Compute_Glint &
      , Compute_Pixel_Arrays &
      , Determine_Level1b_Compression &
      , Expand_space_mask_for_user_limits &
      , Merge_nwp_hires_zsfc &
      , Modify_land_class_with_ndvi &
      , Normalize_reflectances &
      , Quality_control_ancillary_data &
      , Set_bad_pixel_mask &
      , Set_Chan_On_Flag &
      , Set_Solar_Contamination_Mask &
      , Surface_Remote_Sensing &
      , Read_MODIS_White_Sky_Albedo &
      , MODIFY_AUX_CLOUD_TYPE &
      , DAYTIME_REFL_BALANCE_CLOUD_FRACTION

   use CITY_MASK_MOD, only: &
      READ_CITY_MASK

   use PLANCK_MOD, only: &
    populate_planck_tables_sounder

   use RTM_COMMON_MOD, only: &
      NLEVELS_RTM &
      , p_std_rtm  &
      , rtm &
      , ALLOCATE_RTM &
      , DEALLOCATE_RTM

   use CX_ATMOSPHERIC_CORRECTION_VIS_MOD , only: &
       ATMOS_CORR &
      , SETUP_SOLAR_RTM

   use CX_NWP_RTM_MOD, only: &
      MAP_NWP_RTM &
      , CREATE_TEMP_NWP_VECTORS &
      , DESTROY_TEMP_NWP_VECTORS

   use RT_UTILITIES_MOD, only: &
        RTM_NVZEN &
      , GET_PIXEL_NWP_RTM

   use SENSOR_MOD, only: &
      detect_sensor_from_file &
    , output_image_to_screen &
    , output_processing_limits_to_screen &
    , output_sensor_to_screen &
    , read_instr_constants &
    , read_level1b_data &
    , set_data_date_and_time &
    , set_file_dimensions, &
      SET_SENSOR_CHANNEL_MAPPING

   use SFC_EMISS, only: close_seebor_emiss
#ifdef LIBRTTOV
    use CX_RTTOV_SFC_EMISS, only: &
        destroy_rttov_emiss
#endif
   use SIMPLE_COD_065um_MOD,only: &
    compute_simple_solar_cod_065um

   use SIMPLE_COD_LUNAR_MOD, only: &
      compute_simple_lunar_cod_065um

   use SIMPLE_COD_138um_MOD,only: &
    compute_simple_solar_cod_138um

   use SIMPLE_COD_160um_MOD,only: &
    compute_simple_solar_cod_160um

   use SURFACE_PROPERTIES_MOD,only: &
      COMPUTE_BINARY_LAND_COAST_MASKS &
      , setup_umd_props &
      , GET_PIXEL_SFC_REFL_FROM_SFC_TYPE

    use USER_OPTIONS,only: &
      setup_user_defined_options &
    , update_configuration

   use CX_SFC_EMISSIVITY_MOD, only: &
    cx_sfc_emiss_populate_ch &
    , cx_sfc_emiss_correct_for_sfctype

   use TIMER_MOD

   use CX_TIMER_MOD, only: timer_set_up, timer_set_up_all &
      , chronos_rttov, chronos_acha

   use UNIVERSAL_CLOUD_TYPE_MODULE, only: UNIVERSAL_CLOUD_TYPE

   use CLAVRX_STATIC_NAV_MODULE, only: &
       DESTROY_READ_LEVEL1B_FIXED_GRID_STATIC_NAV

   use CX_ABI_LHP_MOD, only: &
       SET_ABI_USE_104um_FLAG

   use cleanup, only: cleanup_tempdir, cleanup_tempdir__exit

   use tracer, only: waitpoint, Set_Tracer_Flag, maybe_clone, update_skip_processing

   implicit none

   !***********************************************************************
   ! Marker: DECLARE VARIABLES
   !***********************************************************************
   integer(kind=int4):: Nrec_Avhrr_Header
   integer(kind=int4):: Ifile
   integer(kind=int4):: File_List_Lun                 !logical unit number for File_list
   integer(kind=int4):: Segment_Number
   integer(kind=int4):: Skip_Processing_Flag
   integer(kind=int4):: Lat_Idx
   integer(kind=int4):: Lon_Idx
   integer(kind=int4):: Zen_Idx
   integer(kind=int4):: Elem_Idx  !generic pixel (along scan) index
   integer(kind=int4):: Line_Idx  !generic line (across scan) index
   logical:: Level1b_Exists

   character(len=1020):: File_1b_Temp
   integer(kind=int4):: erstat
   real(kind=real4):: Time_Since_Launch
   integer(kind=int4):: err_reposnx_Flag


   integer(kind=int4):: ios, nc
   integer(kind=int4):: File_Number
   integer(kind=int4):: Ierror_Level1b
   integer(kind=int4):: ierror_Nwp
   integer(kind=int4):: iperiod16
   integer(kind=int4) :: ierror
   character(len=3):: Day_String
   character(len=1020):: Modis_White_Sky_0_66_Name
   character(len=1020):: Modis_White_Sky_0_86_Name
   character(len=1020):: Modis_White_Sky_1_24_Name
   character(len=1020):: Modis_White_Sky_1_64_Name
   character(len=1020):: Modis_White_Sky_2_13_Name
   character(len=1020):: Snow_Mask_File_Name
   character(len=1020):: oiSst_File_Name
   character(len=4096) :: cmd


   integer(kind=int4):: Emiss_File_Id = missing_value_int4
   integer(kind=int4):: Coast_Mask_Id = missing_value_int4
   integer(kind=int4):: Land_Mask_Id = missing_value_int4
   integer(kind=int4):: Sfc_Type_Id = missing_value_int4
   integer(kind=int4):: Volcano_Mask_Id = missing_value_int4
   integer(kind=int4):: Surface_Elev_Id = missing_value_int4
   integer(kind=int4):: Modis_Alb_0_66_Id = missing_value_int4
   integer(kind=int4):: Modis_Alb_0_86_Id = missing_value_int4
   integer(kind=int4):: Modis_Alb_1_24_Id = missing_value_int4
   integer(kind=int4):: Modis_Alb_1_64_Id = missing_value_int4
   integer(kind=int4):: Modis_Alb_2_13_Id = missing_value_int4
   integer(kind=int4):: Snow_Mask_Id = missing_value_int4
   type(Land_grid_Description) :: Coast_Mask_Str
   type(Land_grid_Description) :: Sfc_Type_Str
   type(Land_grid_Description) :: Land_Mask_Str
   type(Land_grid_Description) :: Volcano_Mask_Str
   type(Land_grid_Description) :: Surface_Elev_Str
   type(Land_grid_Description) :: Modis_Alb_0_66_Str
   type(Land_grid_Description) :: Modis_Alb_0_86_Str
   type(Land_grid_Description) :: Modis_Alb_1_24_Str
   type(Land_grid_Description) :: Modis_Alb_1_64_Str
   type(Land_grid_Description) :: Modis_Alb_2_13_Str
   type(Land_grid_Description) :: Snow_Mask_Str

   ! GOES header structures
   TYPE (AREA_STRUCT) :: AREAstr
   TYPE (GVAR_NAV)    :: NAVstr

   logical :: NLCOMP_Run
   logical :: DCOMP_Run

   character (len = 30) :: string_30
   character (len = 100) :: string_100

   integer:: Chan_Idx
   integer, parameter :: Chan_Idx_Min = 1
   integer, parameter :: Chan_Idx_Max = 44


   character(len = 1024) ::string
   integer :: narg,cptArg
   character(len=30) :: arg_name

   type(timer) :: chrono
   type(timer) :: chrono_all

   !***********************************************************************
   ! Begin Executable Code
   !***********************************************************************


   call timer_set_up_all (chrono_all)
   call timer_set_up (chrono)

   narg=command_argument_count()

   if (narg>0) then
     !loop across options
      do cptArg=1,narg
         call get_command_argument(cptArg,arg_name)

         select case(adjustl(arg_name))
         case("--compile_info","-ci")
             write(*,*)"This is CLAVR-x"
             write(*,*)"Binary compiled on the ",__DATE__," at ",__TIME__
             stop
         case default

         end select
      end do
   end if

   call MESG( '<----------  Start of CLAVRXORB ---------->  $' &
      , level = verb_lev % DEFAULT , color = 4 )
   write(string,*)"Compiled on ",__DATE__//' '//__TIME__
   call MESG( trim(string),level=verb_lev % DEFAULT, color=2)

   !----------------------------------------------------------------------------
   ! Initialize some flags
   !----------------------------------------------------------------------------
   Number_Of_Temporary_Files = 0
   Skip_Processing_Flag = sym%NO
   CALL Set_Tracer_Flag()

   !----------------------------------------------------------------------------
   ! Determine time of the start of all processing
   !----------------------------------------------------------------------------


   call chrono_all%tic(1)

   !------------------------------------------------------------------------------
   ! initialize previous date variables
   !-------------------------------------------------------------------------------

   Sensor%WMO_Id_Previous = 0


   !*************************************************************************
   ! Marker: Read and Quality Check User Defined Options
   !*************************************************************************

   call SETUP_USER_DEFINED_OPTIONS()


  call chronos_rttov % init(['clear_sky.','sfc_emis .','pfaast   .' &
          ,'seebor   .' ] &
     , off=Verbose_Level_Flag .lt. 9)

  call chronos_acha % init([&
          'ACHA SETUP .' &
         ,'ACHA HEIGT .' &
         ,'ACHA COMP  .'&
         ,'Shadow     .'&
         ,'Alloc Main .'&
         ,'LLR Center .'&
         ,'Pass loop  .'&
         ,'Post Proc  .'&
         ,'ACHA total .'&
         ,'FULL RETRI .'&
         ,'NON-FULL   .'&
         ] , off = Verbose_Level_Flag .lt. 9)



   !--- make directory for temporary files created during this run
   nc = len_trim(Temporary_Data_Dir)
   call univ_mkdir_p_f(nc, trim(Temporary_Data_Dir), ierror)

   ! SIGTERM
   call univ_reg_sigterm_handler(cleanup_tempdir__exit)
   ! SIGINT
   call univ_reg_sigint_handler(cleanup_tempdir__exit)

   call initial_pixel_common_alloc
   
   !*************************************************************************
   ! Marker: Open high spatial resolution ancillary data files
   !*************************************************************************
   call OPEN_STATIC_ANCIL_FILES()

   !-----------------------------------------------------------------------
   !--- set up surface radiative properties
   !-----------------------------------------------------------------------
   call SETUP_UMD_PROPS()


   !--- NO LOGIC FOR CH 38.
   call READ_ABI_DUST_LUT()

   !--------------------------------------------------------------------
   !--- setup clock corrections in memory
   !--------------------------------------------------------------------
   if (nav_opt== 2) then
     call SETUP_CLOCK_CORRECTIONS()
   endif

   !**********************************************************************
   ! Marker: Read level2 variables in LEVEL2_LIST
   !**********************************************************************
   call READ_LEVEL2_VAR_LIST()

   !**********************************************************************
   ! Marker: Read file directories in FILE_LIST
   !**********************************************************************

   !--- print to screen which file list is used
   call MESG( "CLAVR-x FILE LIST FILE USED: "//trim(File_list) , level = verb_lev % VERBOSE )

   !--- open file containing list of level1b data to process
   File_List_Lun = GET_LUN()
   open(unit=File_List_Lun, file = trim(File_List),status="old",action="read",iostat=ios)
   if (ios /= 0) then
      write ( string_30, '(i8)') ios
      call MESG("ERROR: Opening clavrxorb_file_list, iostat = "//trim(string_30) , level = verb_lev % ERROR)
      write(*,*) "Stopped in ",__FILE__," at line ",__LINE__
      stop
   endif

   !--- read directories from clavrxorb_input_Files
   read(unit=File_List_Lun,fmt="(a)") Image%Level1b_Path
   read(unit=File_List_Lun,fmt="(a)") Dir_Level2

   !--- reset file counter
   File_Number = 1

   !----------------------------------------------------------------------
   ! Marker: BEGIN LOOP OVER FILES
   !----------------------------------------------------------------------
   File_Loop: do
     !----------------------------------------------------------------------------
     ! Determine time of the start of the processing of this orbit
     !----------------------------------------------------------------------------
      call chrono%tic(16)

      call chrono%tic(17)



      !----------------------------------------------------------------------
      ! Marker: READ IN CLAVRXORB_FILE_LIST AND SET FLAGS
      !----------------------------------------------------------------------
      read(unit=File_List_Lun,fmt="(a)",iostat=Ios) File_1b_Temp
      if ( File_1b_Temp == "") exit
      if (Ios /= 0) then
         if (Ios /= -1) then
            !-- non eof error
            Erstat = 8
            call MESG( "ERROR: Problem reading orbit names from control file" &
               , level = verb_lev % QUIET , color = 1 )
            call cleanup_tempdir
            stop 8
         else
            !-- end of orbits
            if (File_Number == 1) then
               call MESG( "ERROR: No orbits to process, stopping" &
                , level = verb_lev % QUIET , color = 1 )
               call cleanup_tempdir
               stop 404
            endif
            exit
         endif
      endif

      ! --- Initialize Skip_L1b_File_Flag, if true will skipp file
      Skip_L1b_File_Flag = .false.

      !--------------------------------------------------------------
      ! Determine if this level-1b file can be opended, if not skip
      !--------------------------------------------------------------

      !**********************************************************************
      ! Marker: Prepare to read Level-1b file
      !**********************************************************************

      !-- see if level-1b file exists
      Level1b_Exists = file_test(trim(Image%Level1b_Path)//trim(File_1b_Temp))
      if (.not. Level1b_Exists ) then
         call MESG( "ERROR: Level-1b file not found, skipping this file", level=5, color=1)
         print*,trim(Image%Level1b_Path)//trim(File_1b_Temp)
         cycle file_loop
      endif

      ! ---
      !  check if this file is a AVHRR-HIRS File
      !  if so Level1b file is changed to AVHRR file, HIRS File is added to
      !   Image structure to Image%Level1b_Fusion_Name
      if (index(trim(File_1b_Temp),'fusion') > 0) then
         call HIRS_AVHRR_FUSION_PREPERATION(File_1b_Temp)    !mispelled
      end if

      !--------------------------------------------------------------
      !  Determine based on the L1b file's extension if the input file
      !  is gzipped. Announce this fact thru a boolean var.
      !--------------------------------------------------------------
      call DETERMINE_LEVEL1B_COMPRESSION(File_1b_Temp,L1b_Gzip,L1b_Bzip2)

      !------------------------------------------------------------------------
      ! Determine from which sensor this file comes from (MODIS,AVHRR or VIIRS)
      ! and populate sensor structure
      !------------------------------------------------------------------------
      call DETECT_SENSOR_FROM_FILE(AREAstr,NAVstr,Ierror)

      if (Ierror == sym%YES) then
         call MESG ("ERROR: Sensor could not be detected, skipping file " &
            , level = verb_lev % ERROR)
         cycle file_loop
      endif

      !------------------------------------------------------------------------
      ! Having determined the SENSOR, populate information about the
      ! Sensor and Satellite in the sensor and ch structures.  This is needed for each
      ! channel since CLAVR-x supports multi-instrument files (i.e. SSEC Fusion)
      !------------(------------------------------------------------------------
      call SET_SENSOR_CHANNEL_MAPPING()

      !--- print to screen the file name
      call MESG (" " )
      call MESG ("<------------- Next Orbit ---------------> "  &
            ,level = verb_lev % DEFAULT , color = 4)

      !-------------------------------------------------------
      ! reset record counters
      !-------------------------------------------------------
      File_Number = File_Number + 1

      !-----------------------------------------------------------------------
      ! knowing the file type, determine the expected number of elements and
      ! lines in the file
      ! for AVHRR, determine file type and some record lengths
      ! AVHRR Header is read here
      !-----------------------------------------------------------------------

      !--- read in Instrument Constants from appropriate file
      !--- Need this here for static navigation.
      call READ_INSTR_CONSTANTS()

      !------------------------------------------------------------------
      ! update settings according sensor ( algo mode and channel settings
      ! including turn-on and off)
      !------------------------------------------------------------------

      call UPDATE_CONFIGURATION (Sensor%Sensor_Name)

      call SET_FILE_DIMENSIONS(Image%Level1b_Full_Name &
               ,AREAstr &
               ,Nrec_Avhrr_Header &
               ,Ierror)

      if (Ierror == sym%YES) then
         call MESG ("ERROR:  Could not set file dimensions, skipping file " , &
                     level = verb_lev % ERROR)
         cycle file_loop
      endif

      if (Image%Number_Of_Lines <= 0) then
         call MESG ("File dimensions were not set correctly for this sensor "// &
                     sensor%sensor_name , level = verb_lev % ERROR)
         cycle file_loop
      endif

      !-----------------------------------------------------------------------
      !--- set up pixel level arrays (size depends on sensor)
      !-----------------------------------------------------------------------
      !--- determine segment size here
      Line_Idx_Min_Segment = 1
      Line_Idx_Max_Segment = Image % Number_Of_Lines_Per_Segment

      !*************************************************************************
      ! Marker:  READ IN HEADER AND DETERMINE SOME CONSTANTS
      !*************************************************************************

      !----------------------------------------------------------------------
      ! Knowing the sensor, interogate files to start, end date and time
      !----------------------------------------------------------------------
      call SET_DATA_DATE_AND_TIME ( AREAstr)

      !----------------------------------------------------------------------
      ! Output sensor and image parameters to screen
      !----------------------------------------------------------------------
      call OUTPUT_SENSOR_TO_SCREEN()
      call OUTPUT_IMAGE_TO_SCREEN()
      call OUTPUT_PROCESSING_LIMITS_TO_SCREEN()

      !------------------------------------------------------------------
      ! Setup Solar-channel RTM terms for this particular sensor
      !------------------------------------------------------------------
      call SETUP_SOLAR_RTM(Sensor%WMO_Id)
      call SETUP_OLR()
      call SETUP_SST()

      !------------------------------------------------------------------
      ! Create pixel arrays which data for this segment
      !------------------------------------------------------------------
      call CREATE_PIXEL_ARRAYS()

      !------------------------------------------------------------------
      ! Read in Dark Sky Composite
      !------------------------------------------------------------------
      Dark_Composite_Name = "no_file"
      if (Read_Dark_Comp == sym%YES) then
         call DETERMINE_DARK_COMPOSITE_NAME(AREAstr)
      endif


      !*************************************************************************
      ! Marker:  READ IN SENSOR-SPECIFIC CONSTANTS
      !*************************************************************************

      !*************************************************************************
      ! Marker:  Open non-static high spatial resolution ancillary data
      !*************************************************************************

      !--------------------------------------------------------------------
      ! read in sst analysis file
      !--------------------------------------------------------------------

      !Set to use default options. This will be turned off if there is NO daily OISST data
      Oisst_File_Name= GET_OISST_MAP_FILENAME(Image%Start_Year,Image%Start_Doy, &
                                                 trim(OiSst_Data_Dir) )
      if (trim(Oisst_File_Name) == "no_file") then
          Use_Sst_Anal = sym%NO
          call MESG ("WARNING: Could not find daily OISST file", level = verb_lev % WARNING )
      else
          call READ_OISST_ANALYSIS_MAP(OISST_File_Name)
          Use_Sst_Anal = sym%YES
      endif

      !----------------------------------------------------------------------
      ! Open Modis White Sky Albedo Map appropriate for this day
      !----------------------------------------------------------------------
      if (Modis_Clr_Alb_Flag == sym%YES) then
         call OPEN_MODIS_WHITE_SKY_SFC_REFLECTANCE_FILES()
         call MESG("Modis clear albedo map opened successfully")
      endif
      !------------------------------------------------------------------
      ! Open Snow mask file
      !------------------------------------------------------------------
      call GET_SNOW_MASK()

      !*************************************************************************
      ! Marker:  READ IN NWP DATA
      !*************************************************************************

      !--- GFS
      if (NWP_PIX%Nwp_Opt == 1) then
         call READ_GFS_DATA(NWP_PIX%Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time, &
                           Image%End_Year, Image%End_Doy, Image%End_Time, Gfs_Data_Dir, ierror_Nwp)
      endif

      !--- NCEP Reanalysis
      if (NWP_PIX%Nwp_Opt == 2) then
         call READ_NCEP_REANALYSIS_DATA(Image%Start_Year, Image%Start_Doy, Image%Start_Time,  &
                                        Image%End_Doy, Image%End_Time, Ncep_Data_Dir)
      endif

      !--- CFSR
      if (NWP_PIX%Nwp_Opt == 3) then
         call READ_GFS_DATA(NWP_PIX%Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time, &
                            Image%End_Year, Image%End_Doy, Image%End_Time, Cfsr_Data_Dir, ierror_Nwp)
      endif

      !--- GDAS
      if (NWP_PIX%Nwp_Opt == 4) then
         call READ_GFS_DATA(NWP_PIX%Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time,  &
                            Image%End_Year, Image%End_Doy, Image%End_Time, Gdas_Data_Dir, ierror_Nwp)
      endif

      !--- MERRA
      if (NWP_PIX%Nwp_Opt == 5) then
         call READ_GFS_DATA(NWP_PIX%Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time,  &
                            Image%End_Year, Image%End_Doy, Image%End_Time, Merra_Data_Dir, ierror_Nwp)
      endif

      !--- ERA INTERIM ANALYSIS
      if (NWP_PIX%Nwp_Opt == 6) then
         call READ_GFS_DATA(NWP_PIX%Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time,  &
                            Image%End_Year, Image%End_Doy, Image%End_Time, Erai_Data_Dir, ierror_Nwp)
      endif

      !--- GFS AIT Style NWP interpolation
      if (NWP_PIX%Nwp_Opt == 7) then
         call READ_GFS_DATA(NWP_PIX%Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time, &
                           Image%End_Year, Image%End_Doy, Image%End_Time, Gfs_Data_Dir, ierror_Nwp)
      endif

      !--- GFS FV3
      if (NWP_PIX%Nwp_Opt == 8) then
         call READ_GFS_DATA(NWP_PIX%Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time, &
                           Image%End_Year, Image%End_Doy, Image%End_Time, Gfs_Data_Dir, ierror_Nwp)
      endif

      !---- if NWP is being read in, then proceeed in allocating RTM, NWP arrays
      if (NWP_PIX%Nwp_Opt /= 0) then

         !--- Quality control NWP fields
         call QC_NWP()

         !--- create temporary NWP vectors needed for RTM
         call CREATE_TEMP_NWP_VECTORS()

         !--- Compute mappings for NWP and RTM vertical coordinates
         call MAP_NWP_RTM(NWP%NLevels, &
                          NWP%P_Std, &
                          NLevels_Rtm, &
                          P_Std_Rtm)

         !--- allocate RTM structures
         call ALLOCATE_RTM(NWP%Nlon,NWP%Nlat)

      endif


      !*************************************************************************
      ! Marker: Populate or read in other lookup tables
      !*************************************************************************

      !--- planck and aerosol tables
      if (Sensor%WMO_Id /= Sensor%WMO_Id_Previous) then

         call POPULATE_PLANCK_TABLES()

         !--- planck for 11 and 12um sounder ch
         if (trim(Sensor%Sensor_Name) == 'AVHRR-IFF' .or. &
             trim(Sensor%Sensor_Name) == 'AVHRR-FUSION' .or. &
             trim(Sensor%Sensor_Name) == 'VIIRS-IFF' .or. &
             Sensor%Fusion_Flag) then
            call POPULATE_PLANCK_TABLES_SOUNDER()
         endif

         if (Aerosol_Mode > 0 .and. index(Sensor%Sensor_Name,'AVHRR') > 0) then
            call READ_AER_CH123A_REF_LUTS(Ancil_Data_Dir,Sensor%WMO_Id)
         endif

         Sensor%WMO_Id_Previous = Sensor%WMO_Id

      endif

      !--- compute Sun-Earth distance
      Sun_Earth_Distance = 1.0 - 0.016729*cos(0.9856*(Image%Start_Doy-4.0)*dtor)

      !--- compute time since launch needed for reflectance calibration
      Time_Since_Launch = Image%Start_Year + Image%Start_Doy / 365.25  - launch_Date

      !--------------------------------------------------------------
      !-- Marker: Interpolate a clock correction for this orbit
      !--------------------------------------------------------------
      if ((trim(Sensor%Sensor_Name) == 'AVHRR-1') .or. &
          (trim(Sensor%Sensor_Name) == 'AVHRR-2') .or. &
          (trim(Sensor%Sensor_Name) == 'AVHRR-3') .or. &
          (trim(Sensor%Sensor_Name) == 'AVHRR-FUSION') ) then

         if(Nav_Opt == 2) then
            call INTERPOLATE_CLOCK_ERROR(Image%Start_Year, Image%Start_Time,  &
                                   Image%End_Year, Image%End_Time, &
                                   Sensor%WMO_Id,Nav%Timerr_Seconds)
           print *, EXE_PROMPT, "Clock Error = ", Nav%Timerr_Seconds
         endif
      endif

      !*************************************************************************
      ! READ IN DATA RECORDS AND PERFORM QUALITY CHECKING
      !*************************************************************************

      !-------------------------------------------------------------------------------------------
      ! Marker: Begin loop over orbit segments
      !-------------------------------------------------------------------------------------------
      call MESG ( " ")
      call MESG ("Started Processing All Orbital Segments")

      !--- compute number of segments in this orbit
      if (mod(Image%Number_Of_Lines,Image%Number_Of_Lines_Per_Segment) == 0) then
         Image%Number_Of_Segments = Image%Number_Of_Lines / Image%Number_Of_Lines_Per_Segment
      else
         Image%Number_Of_Segments = Image%Number_Of_Lines / Image%Number_Of_Lines_Per_Segment + 1
      endif

      if (trim(Sensor%Sensor_Name) .eq. 'FCI') Image%Number_Of_Segments = 40

      call chrono%tac(17)


      Segment_loop: do Segment_Number = 1,Image%Number_Of_Segments

         !--- save the segment number in a global structure
         Image%Segment_Number  = Segment_Number

         !--- reset skip processing flag
         Skip_Processing_Flag = sym%NO

         !--- reset Goes_Scan_Line_Flag
         Goes_Scan_Line_Flag = sym%NO

         !--- reset pixel arrays to missing for this segment
         call RESET_PIXEL_ARRAYS_TO_MISSING()

         !-----------------------------------------------------------------
         !---- Marker: Read level-1b data
         !-----------------------------------------------------------------

         call chrono % tic(1)

         call READ_LEVEL1B_DATA(Image%Level1b_Full_Name,Segment_Number, &
                                Time_Since_Launch,AREAstr,NAVstr,Nrec_Avhrr_Header &
                                ,Ierror_Level1b)
         if (Ierror_Level1b /= 0) then
            call MESG ("ERROR:  Error reading level1b, skipping this file ",level = verb_lev% ERROR)
            print*,Image%Level1b_Full_Name
            exit
         endif

         if (trim(Sensor%Sensor_Name) == 'SEVIRI' .and. Use_Aux_Flag /= sym%NO_AUX) then
            call READ_OCA(trim(Image%Level1b_Path),trim(Image%Level1b_Name), &
                 Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment, &
                 Segment_Number,Image%Number_Of_Lines,Image%Number_Of_Lines_Read_This_Segment, &
                 Cost_Aux,Pc_Top1_Aux,Pc_Top2_Aux,Pc_Uncertainty1_Aux,Pc_Uncertainty2_Aux, &
                 Tau_Aux,Cld_Phase_Aux)
         endif

         !------------------------------------------------------------------
         ! Apply spatial limits
         !------------------------------------------------------------------
         call EXPAND_SPACE_MASK_FOR_USER_LIMITS(Segment_Number, Geo%Space_Mask)

         !------------------------------------------------------------------
         ! --- if CALIOP file not found skip level1b file
         !------------------------------------------------------------------
         if (Skip_L1b_File_Flag) then
             exit
         endif

         !-------------------------------------------------------------------
         ! Modify Chan_On flags to account for channels read in
         !-------------------------------------------------------------------
         call SET_CHAN_ON_FLAG(Sensor%Chan_On_Flag_Default, Sensor%Chan_On_Flag_Per_Line &
            , is_last_segment = segment_number .EQ. Image%Number_Of_Segments )

         !-------------------------------------------------------------------
         ! Compute Lunar Reflectance
         !-------------------------------------------------------------------
         if ((trim(Sensor%Sensor_Name) == 'VIIRS' .or. trim(Sensor%Sensor_Name) == 'VIIRS-NASA') &
             .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then

           ! - check the angles if this is a good lunar scene
           ! - lun and solar zenith angle
            call COMPUTE_LUNAR_REFLECTANCE (Ch (44) % Rad_Toa &
                      , Geo % Solzen, Geo%Lunzen &
                      , Image % time_start % Year, Image % time_start % Month &
                      , Image % time_start % Day, Image % time_start % msec_of_day &
                      , Geo % Moon_Phase_Angle  &
                      , Ancil_Data_Dir &
                      , Ch(44)%Ref_Lunar_Toa)


           !TEST - EMPIRICAL FIT TO NASA Reflectances to match NOAA - AKH
            if (trim(Sensor%Sensor_Name) == 'VIIRS-NASA') then
               call lunar_reflectance_nasa_data_adjustment(Ch(44)%Ref_Lunar_Toa)
            endif

         end if



         call maybe_clone()
         call waitpoint(1)
         call update_skip_processing(Skip_Processing_Flag)

         !--- update time summation for level-1b processing
         call chrono % tac(1)
         !--- check to see that some data was read in
         if (Image%Number_Of_Lines_Read_This_Segment == 0) then
            call MESG ("WARNING: no scans read in, exiting segment processing loop ",level = verb_lev% WARNING)
            Skip_Processing_Flag = sym%YES
            !exit
         endif


         call chrono % tic(2)
         !---- go no further if no data is read in
         if (Skip_Processing_Flag == sym%NO) then    !skip_Processing_Flag

            !*******************************************************************
            ! Marker: Recompute geolocation
            !*******************************************************************
            !--- use level 1b navigation
            Nav%Lat = Nav%Lat_1b
            Nav%Lon = Nav%Lon_1b

            !---  AVHRR Repositioning
            if (trim(Sensor%Sensor_Name) == 'AVHRR-1' .or. &
                trim(Sensor%Sensor_Name) == 'AVHRR-2' .or. &
                trim(Sensor%Sensor_Name) == 'AVHRR-3' ) then

               if (Nav_Opt == 2) then

                  if ((Nav%Timerr_Seconds .NEfp. Missing_Value_Real4) .and. &
                   (Nav%Timerr_Seconds .NEfp. 0.0)) then

                     call REPOSITION_FOR_CLOCK_ERROR(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment, &
                                                  Nav%Timerr_Seconds,Err_Reposnx_Flag)
                  else
                     Nav%Lat = Nav%Lat_1b
                     Nav%Lon = Nav%Lon_1b
                  endif

               endif

            endif

            if (NWP_PIX%Nwp_Opt /= 0) then

               !--- map each each into correct NWP cell
               call MAP_PIXEL_NWP(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)

               !--- compute needed NWP levels (sfc, tropo, inversion, ...)
               call COMPUTE_NWP_LEVELS_SEGMENT(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)

               !--- compute desired nwp parameters
               call COMPUTE_SEGMENT_NWP_CLOUD_PARAMETERS()

            endif

            !--- determine a pixel-level mask to exclude data solar contaminiation
            call SET_SOLAR_CONTAMINATION_MASK(Solar_Contamination_Mask)


            !--- determine a pixel-level mask to exclude bad or unprocessible data
            !--- call SET_BAD_PIXEL_MASK(Bad_Pixel_Mask)
            call SET_BAD_PIXEL_MASK(Bad_Pixel_Mask,ABI_Use_104um_Flag)

            !*******************************************************************
            ! Marker: Interpolate ancillary data to each pixel in this segment
            !*******************************************************************

            !--- map nwp parameters to pixel projection
            call COMPUTE_PIXEL_NWP_PARAMETERS(NWP_PIX%Smooth_Nwp_Flag)


            !--- map non-nwp ancillary data to the pixel projection
            call MAP_ANCIL_DATA_TO_PIXEL_GRID()

             ! - set surface emissivity
            call  CX_SFC_EMISS_POPULATE_CH

            !--- post process dark composite if one read in
            if (Read_Dark_Comp == sym%YES .and. Dark_Composite_Name /= "no_file") then
               call POST_PROCESS_GOES_DARK_COMPOSITE(Ref_Ch1_Dark_Composite)
            endif

            !--- check ancillary data and modify Bad_Pixel_Mask accordingly
            call QUALITY_CONTROL_ANCILLARY_DATA(Bad_Pixel_Mask)

            !-----------------------------------------------------------------------
            !--- Marker: Compute some fundamental pixel-level arrays
            !-----------------------------------------------------------------------
            !--- compute some common used pixel arrays
            call COMPUTE_PIXEL_ARRAYS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            !--- normalize reflectances by the solar zenith angle and sun-earth distance
            call NORMALIZE_REFLECTANCES(Sun_Earth_Distance)

            !--- do we want to change reflectances only if 10_4 flag is used?
            !--- compute the channel 20 pseudo reflectance
            if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and.  Sensor%Chan_On_Flag_Default(31) == sym%YES) then
              call CH20_PSEUDO_REFLECTANCE(Solar_Ch20_Nu,Geo%CosSolzen,ch(20)%Rad_Toa,31,ch(31)%Bt_Toa, &
                                           Sun_Earth_Distance,ch(20)%Ref_Toa,Ch(20)%Emiss_Rel_11um)
            endif
            if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and.  Sensor%Chan_On_Flag_Default(38) == sym%YES) then
              call CH20_PSEUDO_REFLECTANCE(Solar_Ch20_Nu,Geo%CosSolzen,ch(20)%Rad_Toa,38,ch(38)%Bt_Toa, &
                                           Sun_Earth_Distance,ch(20)%Ref_Toa,Ch(20)%Emiss_Rel_10_4um)
            endif

            !--- COMPUTE SPATIAL METRICS FROM LEVEL-1B Terms

            !--- populate basic spatial uniformity arrays(min, max, mean and std)
            call COMPUTE_MIN_MAX_MEAN_STD_METRICS()

            !--- call spatial correlation routines
            call COMPUTE_SPATIAL_CORRELATION_ARRAYS()

            !--- apply median filters
            call COMPUTE_MEDIAN_METRICS_L1B()

            !--- correct for noise in 3.75 micron (important for older AVHRRs)
            if (Sensor%Chan_On_Flag_Default(20)==sym%YES .and. &
                Sensor%Chan_On_Flag_Default(31) == sym%YES) then
                call MITIGATE_CH20_NOISE(Sensor%WMO_Id,ch(20)%Bt_Toa,Bt_Ch20_Median_5x5,ch(31)%Bt_Toa)
            endif

            !--- compute pixel level Snow map based on all ancillary data

            !----- IS ALL DATA AVAILABLE HERE?
            if (NWP_PIX%Nwp_Opt /= 0) then
               call COMPUTE_SNOW_CLASS(Sfc%Snow_NWP,Sfc%Snow_OISST, &
                                       Sfc%Snow_IMS,Sfc%Snow_GLOB, &
                                       Sfc%Land,Sfc%Snow,ABI_Use_104um_Flag)
            endif




            !--- SST if possible for this sensor (needs snow info for masking)
            !--- This subroutine uses bt85, bt104, bt11, bt12).  Will fail in
            !--- LHP events.

            call COMPUTE_SST()

            !--- fill in holes in surface emiss and surface refl with default
            !--- values based on surface type
            call CX_SFC_EMISS_CORRECT_FOR_SFCTYPE()

            !--- compute desert mask cloud detection
            Sfc%Desert_Mask =  DESERT_MASK_FOR_CLOUD_DETECTION(ch(20)%Sfc_Emiss, Nav%Lat, Sfc%Snow, Sfc%Sfc_Type)

            !--- compute city mask cloud detection
            if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then
               !Sfc%City_Mask =  CITY_MASK_FOR_CLOUD_DETECTION(ch(44)%Rad_Toa, Sfc%Sfc_Type)
               call READ_CITY_MASK()
            endif

            !--- update time summation for level-1b processing
            call chrono % tac(2)
            !*******************************************************************
            ! Marker: Compute nwp mapping and rtm values for each pixel in segment
            !*******************************************************************

            !--- needs an NWP to run.
            if (NWP_PIX%Nwp_Opt /= 0) then

               call chrono % tic(3)
               !-- temporally interp skin temp for each segment (only ncep reanalysis)
               if (NWP_PIX%Nwp_Opt == 2) then
                  call TEMPORAL_INTERP_TMPSFC_NWP(Image%Scan_Time_Ms(Line_Idx_Min_Segment),  &
                                Image%Scan_Time_Ms(Line_Idx_Min_Segment+Image%Number_Of_Lines_Read_This_Segment-1))
               endif

               !--- compute a surface temperature from the NWP
               call MODIFY_TSFC_NWP_PIX(1,Image%Number_Of_Elements,Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)


                !--- if an sst analysis is available, use that
                ! TODO find a better place
               if   (Use_Sst_Anal == sym%YES) then
                  where ((Sfc%Land_Mask == sym%NO) .and.  &
                    (Sfc%Snow == sym%NO_SNOW) .and.  &
                    (Sst_Anal > 270.0 ))

                     NWP_PIX%Tsfc = Sst_Anal

                  end where
               end if

               !--- compute pixel-level rtm parameters
               call GET_PIXEL_NWP_RTM(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
               call COMPUTE_MEDIAN_METRICS_L2()
               if (ABI_Use_104um_Flag) then
                   call COMPUTE_RADIATIVE_CENTER_ARRAYS(LRC_Flag,38)
               else
                   call COMPUTE_RADIATIVE_CENTER_ARRAYS(LRC_Flag,31)
               endif

               call OPAQUE_CLOUD_HEIGHT(ABI_Use_104um_Flag)

               !--- apply atmospheric correction - needs rtm results
               call ATMOS_CORR(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

               !--- compute the channel 20 pseudo reflectance for clear-skies
               if (Sensor%Chan_On_Flag_Default(20) == sym%YES &
                  .and.  Sensor%Chan_On_Flag_Default(31) == sym%YES) then

                  call CH20_PSEUDO_REFLECTANCE(Solar_Ch20_Nu &
                                             , Geo%CosSolzen &
                                             , ch(20)%Rad_Toa_Clear &
                                             , 31 &
                                             , ch(31)%Bt_Toa_Clear &
                                             , Sun_Earth_Distance &
                                             , ch(20)%Ref_Toa_Clear &
                                             , ch(20)%Emiss_Rel_11um_Clear)
               endif

               if (Sensor%Chan_On_Flag_Default(20) == sym%YES &
                   .and.  Sensor%Chan_On_Flag_Default(38) == sym%YES) then

                  call CH20_PSEUDO_REFLECTANCE(Solar_Ch20_Nu &
                                              ,Geo%CosSolzen &
                                              ,ch(20)%Rad_Toa_Clear &
                                              ,38 &
                                              ,ch(38)%Bt_Toa_Clear &
                                              ,Sun_Earth_Distance &
                                              ,ch(20)%Ref_Toa_Clear &
                                              ,ch(20)%Emiss_Rel_10_4um_Clear)
               endif

               !--- compute surface products (Tsfc,Ndvi,Rsr ...)
               call SURFACE_REMOTE_SENSING(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment,ABI_Use_104um_Flag)

               call chrono % tac(3)
            endif

            !*******************************************************************
            ! Marker: Spatial metrics processing
            !*******************************************************************
            call chrono % tic(4)
            !------------------------------------------------------------------------------
            ! compute glint masks for use in cloud detection
            !------------------------------------------------------------------------------

            if (ABI_Use_104um_Flag) then
               if (Sensor%Chan_On_Flag_Default(38) == sym%YES) then
                  !--- solar glint mask
                  if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
                     call COMPUTE_GLINT(Geo%Glintzen, Ch(1)%Ref_Toa, ch(1)%Ref_Toa_Std_3x3, &
                                     Sfc%Glint_Mask,ABI_Use_104um_Flag)
                  endif

                  !--- lunar glint mask
                  if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then
                     call COMPUTE_GLINT(Geo%Glintzen_Lunar,ch(44)%Ref_Lunar_Toa,  &
                                     ch(44)%Ref_Lunar_Std_3x3, Sfc%Glint_Mask_Lunar, &
                                     ABI_Use_104um_Flag)
                  endif

               endif
            else
               if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
                  !--- solar glint mask
                  if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
                     call COMPUTE_GLINT(Geo%Glintzen, Ch(1)%Ref_Toa, ch(1)%Ref_Toa_Std_3x3, &
                                     Sfc%Glint_Mask,ABI_Use_104um_Flag)
                  endif

                  !--- lunar glint mask
                  if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then
                     call COMPUTE_GLINT(Geo%Glintzen_Lunar,ch(44)%Ref_Lunar_Toa,  &
                                     ch(44)%Ref_Lunar_Std_3x3, Sfc%Glint_Mask_Lunar, &
                                     ABI_Use_104um_Flag)
                  endif

               endif
            endif

            !*******************************************************************
            ! Marker: Definition of pixel hdf files
            !  if first segment, define pixel hdf files
            ! this has to be here because cal coefficients are only known
            ! after data has been read in and calibrated
            !*******************************************************************
            if ((Segment_Number == 1) .and. (Skip_Output == 0)) then

               !--- place algorithm cvs tags into global strings for output
               call SET_CLOUD_TYPE_VERSION()
               call SET_IR_CLOUD_TYPE_VERSION()

               Num_Scans_Level2_Hdf = 0

               !----------------------------------------------------------------------
               ! Define Level2 Variables - needs to be called after CREATE_PIXEL_ARRAYS
               !----------------------------------------------------------------------
               ! -----------
               !  All level2 output variables have to allocated before Setup_level2_sds_info
               !  Also Muri
               ! -

               if ( (trim(Sensor%Sensor_Name) == 'AHI' .or. trim(Sensor%Sensor_Name) == 'AHI9') .and. Aerosol_Mode == 1) then
                  call muri % allocate (Image%Number_of_elements , Image%Number_Of_Lines_Read_This_Segment)
               end if

               !----------------------------------------------------------------------
               ! set up level2 file
               !----------------------------------------------------------------------
               call SETUP_LEVEL2_SDS_INFO()

               call DEFINE_HDF_FILE_STRUCTURES(Image%Number_Of_Lines, &
                              Dir_Level2, &
                              Image%Level1b_Name, &
                              Level2_File_Flag, &
                              c1,c2,planck_a1(20),planck_a2(20),planck_nu(20), &
                              planck_a1(31),planck_a2(31),planck_nu(31), &
                              planck_a1(32),planck_a2(32),planck_nu(32),Solar_Ch20_Nu,&
                              Sun_Earth_Distance,Therm_Cal_1b, &
                              Ref_Cal_1b,Nav_Opt,Use_Sst_Anal, &
                              Modis_Clr_Alb_Flag,NWP_PIX%Nwp_Opt, &
                              Ch1_Gain_Low,Ch1_gain_High, &
                              Ch1_Switch_Count_Cal,Ch1_Dark_Count_Cal, &
                              Ch2_Gain_low,Ch2_Gain_High, &
                              Ch2_Switch_Count_Cal,Ch2_Dark_Count_Cal, &
                              Ch3a_Gain_low,Ch3a_Gain_High, &
                              Ch3a_Switch_Count_Cal,Ch3a_Dark_Count_Cal, &
                              Image%Start_Year,Image%End_Year,Image%Start_Doy,Image%End_Doy,&
                              Image%Start_Time,Image%End_Time)
            endif


            call chrono % tac(4)
            !*******************************************************************
            ! Marker: Generate pixel-level products
            !*******************************************************************

            !---- pixel level aerosol
            if (index(Sensor%Sensor_Name,'AVHRR') > 0 .and. Aerosol_Mode > 0) then


               call chrono % tic(5)
               call PIXEL_AER_RET_OCEAN(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
               call chrono % tac(5)
            endif

            !--- only apply cloud mask and type routines if nwp/rtm information available
            if (Cld_Flag == sym%YES .and. NWP_PIX%Nwp_Opt > 0) then


               call chrono % tic(6)
               !--- simple cloud optical depths (no effective radius or atmos corr). Used for masking
               if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
                  call COMPUTE_SIMPLE_SOLAR_COD_065um(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)
               endif
               if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
                  call COMPUTE_SIMPLE_SOLAR_COD_160um(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)
               endif
               if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then
                  call COMPUTE_SIMPLE_LUNAR_COD_065um(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)
               endif
               if (Sensor%Chan_On_Flag_Default(26) == sym%YES) then
                  call COMPUTE_SIMPLE_SOLAR_COD_138um(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)
               endif

               !--- direct cloud height retrievals used as inputs to typing and acha
               call CO2IRW_CLOUD_HEIGHT()
               call SPLITWIN_CLOUD_HEIGHT()
               call H2O_INTERCEPT_CLOUD_HEIGHT(Image%Number_Of_Elements,Line_Idx_Min_Segment, &
                                               Image%Number_Of_Lines_Read_This_Segment, &
                                               P_Std_Rtm, &
                                               Pc_H2O,Tc_H2O,Zc_H2O)
               call CO2_SLICING_CLOUD_HEIGHT(Image%Number_Of_Elements,Line_Idx_Min_Segment, &
                                             Image%Number_Of_Lines_Read_This_Segment, &
                                             P_Std_Rtm, &
                                             Pc_Co2,Tc_Co2,Zc_Co2)

               if (Nucaps_Flag) then   !Sets Tc_Cirrus_Background
                   call VIIRS_NUCAPS(Segment_Number)
                   Cloud_Fraction_Background = NUCAPS%Cld_Fraction
               endif

               !--- cloud mask
               if (Use_Aux_Flag == sym%USE_AUX .and. Cloud_Mask_Aux_Read_Flag == sym%YES) then

                  CLDMASK%Cld_Mask = CLDMASK%Cld_Mask_Aux

               else

                  if (Cloud_Mask_Bayesian_Flag == sym%ECM1) then
                     Cloud_Mask_Mode = 'ecm1'
                     call NB_CLOUD_MASK_BRIDGE(Segment_Number)
                  elseif (Cloud_Mask_Bayesian_Flag == sym%NO) then
                     Cloud_Mask_Mode = 'bcm'
                     call BASELINE_CLOUD_MASK_MAIN(Segment_Number)
                  elseif (Cloud_Mask_Bayesian_Flag == sym%ECM2) then
                     Cloud_Mask_Mode = 'ecm2'
                     call ECM2_CLOUD_MASK_BRIDGE(Segment_Number)

                     ! if use MODAWG set (3) and it is read in save it as cloud mask (Denis B. 2020-08-03)
                     if (Use_Aux_Flag == sym%USE_AUX_MODAWG .and. Cloud_Mask_Aux_Read_Flag == sym%YES) then
                        CLDMASK%Cld_Mask = CLDMASK%Cld_Mask_Aux
                     endif

                  else
                     call MESG("Only the Enterprise and Baseline Cloud Masks are available, check selection",level=verb_lev % ERROR)
                     stop
                  endif

                  if (Use_Aux_Flag == sym%USE_AUX .and. Cloud_Mask_Aux_Read_Flag == sym%NO) then
                     call MESG("Auxilliary Cloud Mask Selected but Unavailable ", level=verb_lev % WARNING)
                  endif

               endif

               !--- Compute Adjacent pixels Cloud Mask
               call ADJACENT_PIXEL_CLOUD_MASK(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

               !--- compute datyime reflectance balance subpixel cloud fraction

               if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. &
                   Ch(1)%Sub_Pixel_On_Flag .and. Static_Dark_Sky_Flag ) then

                   call DAYTIME_REFL_BALANCE_CLOUD_FRACTION(Ch(1)%Ref_Toa, &
                                                            Ch(1)%Ref_Toa_Max_Sub, &
                                                            Static_Ref_065um_Dark_Composite, &
                                                            Geo%SolZen, &
                                                            Subpixel_Cloud_Fraction)
               endif

               !--- Dust Mask
               if (Use_ABI_Dust == sym%YES) then
                  if (Sensor%WMO_Id == 270 .or. Sensor%WMO_Id == 271 .or. Sensor%WMO_Id == 272) then
                    CLDMASK%Dust_Mask = sym%NO
                    call GET_SEGMENT_ABI_DUST_PROB()
                  endif
               endif

               !--- if dark composite available for GOES-1km data, do this
               if (trim(Sensor%Sensor_Name) == 'GOES-IL-IMAGER' .or. trim(Sensor%Sensor_Name) == 'GOES-MP-IMAGER') then
                 if (Read_Dark_Comp == sym%YES .and. Dark_Composite_Name /= "no_file") then
                   call DARK_COMPOSITE_CLOUD_MASK(CLDMASK%Cld_Mask)
                 endif
               endif

               !--- accumulate cloud mask metrics
               call COMPUTE_CLOUD_MASK_PERFORMANCE_METRICS(Cloud_Mask_Count,Nonconfident_Cloud_Mask_Count)

               call chrono % tac(6)

               !--- cloud type

               call chrono % tic(7)
               !--- make sure MODIS aux phase/type conforms to CLAVR-x expectations
               if (Cloud_Type_Aux_Read_Flag == sym%YES .and. Sensor%Sensor_Name == 'MODIS') then
                  call MODIFY_AUX_CLOUD_TYPE()
               endif

               if (Use_Aux_Flag == sym%USE_AUX .and. Cloud_Type_Aux_Read_Flag == sym%YES) then

                  Cld_Type = Cld_Type_Aux
                  Cld_Phase = Cld_Phase_Aux

               else

                  !--- if ECM2 is not called, compute type and phase here
                  if (Cloud_Mask_Bayesian_Flag /= sym%ECM2) then
                     call CLOUD_TYPE_BRIDGE()
                     !Cld_Type_Aux = Cld_Type
                     !call UNIVERSAL_CLOUD_TYPE()

                     !--- WAITPOINT Cloud Phase Complete (if not ECM2)
                     call waitpoint(4)
                  endif

                  !--- Bryan Baum MODIS C6 Cloud Type
                  call IR_CLOUD_TYPE_BAUM()

                  if (Use_IR_Cloud_Type_Flag) then
                     Cld_Type = Cld_Type_IR
                     Cld_Phase = Cld_Phase_IR
                  endif

              endif


              call chrono % tac(7)
            endif   !end of Cld_Flag check

            !--------------------------------------------------------------------
            !   Compute Cloud Properties (Height, Optical Depth, ...)
            !--------------------------------------------------------------------
            if (Cld_Flag == sym%YES .and. NWP_PIX%Nwp_Opt > 0) then

               !---------------------------------------------------------------------
               ! ACHA Section
               !---------------------------------------------------------------------

               call chrono % tic(8)
               !--->call CTP_MULTILAYER()  !WHAT IS THIS???

               !-------------------------------------------------------------------
               ! make co2 slicing height from sounder with using sounder/imager IFF
               !-------------------------------------------------------------------
               if (index(Sensor%Sensor_Name,'IFF') > 0) then

                  if (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') then

                      call SOUNDER_EMISSIVITY()
                      call MODIFY_CLOUD_TYPE_WITH_SOUNDER (Cld_Temp_Sounder, Cld_Emiss_Sounder, Cld_Type)
                      call MAKE_CIRRUS_PRIOR_TEMPERATURE_FROM_CO2(Cld_Temp_Sounder, Cld_Press_Sounder, Cld_Emiss_Sounder,  &
                                                                  Tc_Cirrus_Background, Zc_Cirrus_Background)

                  else

                    call MODIFY_CLOUD_TYPE_WITH_SOUNDER (Tc_CO2, Ec_CO2, Cld_Type)
                    call MAKE_CIRRUS_PRIOR_TEMPERATURE_FROM_CO2(Tc_Co2, Pc_Co2, Ec_Co2,  &
                                                                Tc_Cirrus_Background, Zc_Cirrus_Background)

                  endif

               endif

               if (Nucaps_Flag) then   !Sets Tc_Cirrus_Background
                   call CONVERT_SMOOTH_NUCAPS_TEMP
                   Tc_Cirrus_Background = NUCAPS%Cld_Temp_Smoothed
                   Zc_Cirrus_Background = NUCAPS%Cld_Height
               endif

               !--- Single Channel (11 um) opaque solution
               if (trim(ACHA%Mode) == "off") then
                  call MODE_ZERO_CLOUD_HEIGHT(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
               endif

               !--- GOES-R Baseline (11,12,13.3 um)
               if (trim(ACHA%Mode) == "baseline") then
                  call BASELINE_CLOUD_HEIGHT_MAIN()
               endif

               !--- General ACHA Solution
               if (trim(ACHA%Mode) /= "off" .and. trim(ACHA%MODE) /= "baseline") then

                  !--- AWG CLoud Height Algorithm (ACHA) and associated products
                  call AWG_CLOUD_HEIGHT_BRIDGE()

                  if (ASOS%Mode > 0) call ASOS_BRIDGE()

                  !--- interpolate NWP wind and tpw profiles at cloud-top level
                  call COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND_AND_TPW(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

                  !--- interpolate ACHA cloud heights to flight level altitude.
                  call COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment,ACHA%Pc,ACHA%Alt)

                  !--- interpolate freezing level to flight level altitude.
                  call COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment, &
                                                      NWP_PIX%FrzPre,NWP_PIX%FrzAlt)
                  call COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment, &
                                                      NWP_PIX%HomoFrzPre,NWP_PIX%HomoFrzAlt)

                  !--accumulate performance metrics
                  call COMPUTE_ACHA_PERFORMANCE_METRICS(ACHA%Processed_Count,ACHA%Valid_Count,ACHA%Success_Fraction)

                  !-- make CSBT masks (Clear Sky Brightness Temperature)
                  call OPAQUE_TRANSMISSION_HEIGHT()
                  call COMPUTE_CSBT_CLOUD_MASKS()

               endif

               !--- Convective Cloud Probability
               if (Sensor%Chan_On_Flag_Default(27) == sym%YES) then
                  call CONVECTIVE_CLOUD_PROBABILITY(Bad_Pixel_Mask,  &
                       ch(31)%Bt_TOA, ch(27)%Bt_TOA, Ch(31)%Emiss_Tropo,  &
                       NWP_PIX%Tsfc, ACHA%Conv_Cld_Prob)
               else  ! ch(27)%Bt_TOA is not allocated/used, pass ch(31) instead
                  call CONVECTIVE_CLOUD_PROBABILITY(Bad_Pixel_Mask,  &
                       ch(31)%Bt_TOA, ch(31)%Bt_TOA, Ch(31)%Emiss_Tropo,  &
                       NWP_PIX%Tsfc, ACHA%Conv_Cld_Prob)
               end if
               call SUPERCOOLED_CLOUD_PROBABILITY(Bad_Pixel_Mask,Cld_Type,ACHA%Tc,ACHA%Supercooled_Cld_Prob)

               call chrono % tac(8)
               !-----------------------------------------------------------------------------------
               ! NLCOMP Section
               !-----------------------------------------------------------------------------------


               NLCOMP_Run = .false.
               if ((trim(Sensor%Sensor_Name) == 'VIIRS' &
                    .or. trim(Sensor%Sensor_Name) == 'VIIRS-NASA') &
                    .and. Sensor%Chan_On_Flag_Default(44) == sym % yes &
                    .and. NLCOMP_Mode > 0) then
                  if ( count (ch(44)%Ref_Lunar_Toa > 0) > 0 ) then
                     call chrono % tic(9)
                     call AWG_CLOUD_DNCOMP_ALGORITHM( Iseg_In = Segment_Number , NLCOMP_Mode = .true. &
                        , algorithm_started = NLCOMP_Run)
                     call chrono % tac(9)
                  endif
               endif


               !-----------------------------------------------------------------------------------
               ! DCOMP Section
               !-----------------------------------------------------------------------------------

               DCOMP_Run = .false.

               if (DCOMP_Mode > 0) then
                 call chrono % tic(10)
                 call AWG_CLOUD_DNCOMP_ALGORITHM( Iseg_In = Segment_Number , algorithm_started = dcomp_run, version = dcomp_version )
                 call chrono % tac(10)
               endif

               !-----------------------------------------------------------------------------------
               ! Derived Product Section
               !-----------------------------------------------------------------------------------


               if ( DCOMP_Run .or. NLCOMP_Run) then
                  call chrono % tic(11)
                     call COMPUTE_CLOUD_WATER_PATH(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
                     call COMPUTE_DCOMP_INSOLATION(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment,Sun_Earth_Distance)
                     call COMPUTE_ADIABATIC_CLOUD_PROPS(Line_Idx_Min_segment,Image%Number_Of_Lines_Read_This_Segment)
                     call COMPUTE_DCOMP_PERFORMANCE_METRICS(DCOMP_Processed_Count,DCOMP_Valid_Count)
                     call ADJUST_DCOMP_LWP()
                     call COMPUTE_PRECIPITATION(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

                  if (trim(Sensor%Sensor_Name) == 'AHI' .or. trim(Sensor%Sensor_Name) == 'AHI9') then
                     call COMPUTE_PRECIPITATION_AHI(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
                  end if
                  call chrono % tac(11)
               endif

               !-----------------------------------------------------------------------------------
               ! Cloud Altitude, Base and CCL Section (if ACHA was executed)
               !-----------------------------------------------------------------------------------

               call chrono % tic(12)
               if (trim(ACHA%Mode) /= "off") then

                 call CLOUD_BASE_BRIDGE()

                 !---Calculate flight level altitude of the base for AWC
                 call COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment,BASE%Pc_Base,BASE%Base_Alt)

                 !---Calculate flight level altitude of the lower top and base for AWC
                 call COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment,ACHA%Lower_Pc,ACHA%Lower_Alt)
                 call COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment,BASE%Lower_Base_Pc,BASE%Lower_Base_Alt)

                 call CCL_BRIDGE()

                 !--- compute ice and liquid water contents
                 call COMPUTE_MASS_CONCENTRATION()

               endif

               call chrono % tac(12)
            endif

            !-----------------------------------------------------------------------------------
            ! - MURI = Multidisciplinary Reseach Initiative  Aerosol
            !-----------------------------------------------------------------------------------
            if ( (trim(Sensor%Sensor_Name) == 'AHI' .or. trim(Sensor%Sensor_Name) == 'AHI9') .and. Aerosol_Mode == 1) then

               call chrono % tic(13)
               call CX_MURI_ALGORITHM (Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)


                   call chrono % tac(13)
            endif


            !-----------------------------------------------------------------------------------
            !--- radiative flux parameters
            !-----------------------------------------------------------------------------------

            call chrono % tic(14)
            !---  OLR
            call COMPUTE_OLR()

            !---  SASRAB
            if ( Sasrab_Flag == sym%YES) then
              call INSOLATION(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
            end if

            call chrono % tac(14)
            !--- assign quality flags for things like sst and aerosol (clear-sky)
            if (NWP_PIX%Nwp_Opt > 0 .and. Cld_Flag == sym%YES) then
               !--- assign clear sky quality flags
               call ASSIGN_CLEAR_SKY_QUALITY_FLAGS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
            endif

            !-----------------------------------------------------------------------------------
            !--- for fusion files reset replaced 3.75 bt to missing so output isn't confusing
            !-----------------------------------------------------------------------------------
            if (AVHRR_Fusion_Flag) then
               call SET_REPLACED_AVHRR_TO_MISSING()
            endif

            !-- WAITPOINT before write
            call waitpoint(3)

            !*******************************************************************
            ! Marker: Write to output files (pixel-level)
            !*******************************************************************

            call chrono % tic(15)
            if(Skip_Output == 0) then
            call WRITE_SEGMENT_LEVEL2(Level2_File_Flag)
            endif
            call chrono % tac(15)
            !*************************************************************************
            ! Marker: RTM Structure Memory Deallocation
            !*************************************************************************

            !--- deallocate rtm profile arrays (only do if NWP is used)
            if (NWP_PIX%Nwp_Opt /= 0) then

               do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
                  do Elem_Idx = 1, Image%Number_Of_Elements
                     Lon_Idx = NWP_PIX%I_Nwp(Elem_Idx,Line_Idx)
                     Lat_Idx = NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)
                     Zen_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)
                     !--- check to see if Lon_Idx and Lat_Idx are valid, if deallocate
                     if ((Lon_Idx >= 0) .and. (Lon_Idx <= NWP%Nlon) .and. &
                         (Lat_Idx >= 0) .and. (Lat_Idx <= NWP%Nlat) .and. &
                         (Zen_Idx >= 0) .and. (Zen_Idx <= Rtm_Nvzen)) then

                        if (Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx_Rtm(Elem_Idx,Line_Idx))%is_set ) then
                          do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max
                              call Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%dealloc
                          end do
                          Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%is_set = .false.
                        endif

                     endif
                  end do
               end do

               !--- deallocate rtm cells
               do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
                  do Elem_Idx = 1, Image%Number_Of_Elements

                     Lon_Idx = NWP_PIX%I_Nwp(Elem_Idx,Line_Idx)
                     Lat_Idx = NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)

                     !--- check to see if Lon_Idx and Lat_Idx are valid, if deallocate
                     if ((Lon_Idx >= 0) .and. (Lon_Idx <= NWP%Nlon) .and. &
                           & (Lat_Idx >= 0) .and. (Lat_Idx <= NWP%Nlat)) then
                          call Rtm(Lon_Idx,Lat_Idx)% DEALLOC()

                     endif
                  end do
               end do

            endif

            !--- screen output to mark progress through orbit
           write ( string_100, '(A22, I3, A16, I5 , 4X , I5)')  &
                        "processed segment #",Segment_Number," scanlines = ",  &
                        Image%Scan_Number(Line_Idx_Min_Segment), &
                        Image%Scan_Number(Image%Number_Of_Lines_Read_This_Segment)
           call MESG  ( string_100 )

        !--- WAITPOINT saving complete
        call waitpoint(6)

         endif   !end Skip_Processing_Flag condition
        !*************************************************************************
        ! Marker: End of loop over orbital segments
        !*************************************************************************
       call chronos_acha % summary(title = 'ACHA', sort=.true.)
       call chronos_rttov % summary(title = 'RTTOV', sort=.true.)
      end do Segment_loop

      call MESG ( "Finished Processing All Orbital Segments")
      call MESG ( " ")

      !--- WAITPOINT finished all segments
      call waitpoint(7)

      !*************************************************************************
      !   Marker: Close output pixel-level files
      !*************************************************************************

      !*************************************************************************
      ! Marker: Close non-static ancillary data files
      !*************************************************************************

      if (Emiss_File_Id > 0) call CLOSE_SEEBOR_EMISS(Emiss_File_Id)
      if (Modis_Alb_0_66_Id > 0) call CLOSE_LAND_SFC_HDF(Modis_Alb_0_66_Id)
      if (Modis_Alb_0_86_Id > 0) call CLOSE_LAND_SFC_HDF(Modis_Alb_0_86_Id)
      if (Modis_Alb_1_24_Id > 0) call CLOSE_LAND_SFC_HDF(Modis_Alb_1_24_Id)
      if (Modis_Alb_1_64_Id > 0) call CLOSE_LAND_SFC_HDF(Modis_Alb_1_64_Id)
      if (Modis_Alb_2_13_Id > 0) call CLOSE_LAND_SFC_HDF(Modis_Alb_2_13_Id)
#ifdef LIBRTTOV
      if (Use_Land_IR_Emiss > 1) then
            call DESTROY_RTTOV_EMISS
            call MESG("DESTROY RTTOV_EMIS ")
       endif
#endif

      !*************************************************************************
      !Marker: Deallocate remaining arrays
      !*************************************************************************
      !--- main RTM structures and NWP arrays
      if (NWP_PIX%Nwp_Opt > 0) then
         call DESTROY_NWP_ARRAYS()
         call DESTROY_TEMP_NWP_VECTORS()
         call DEALLOCATE_RTM()
      endif

      !--- Deallocate memory from pixels arrays
      call DESTROY_PIXEL_ARRAYS()

      !--- Deallocate memory from static navigation
      if (Image%Static_Nav_Flag) then
         call DESTROY_READ_LEVEL1B_FIXED_GRID_STATIC_NAV
      endif

      !--- remove files in temporary directory and then the directory
      if (Number_Of_Temporary_Files > 0) then
         do Ifile = 1, Number_Of_Temporary_Files
            call MESG("Removing Temporary File: "//trim(Temporary_File_Name(Ifile)))
            cmd = trim(Temporary_Data_Dir)//  &
                 trim(Temporary_File_Name(Ifile)); nc = len_trim(cmd)
            call univ_remove_f(nc, trim(cmd), ierror)
         end do
      endif
      Number_Of_Temporary_Files = 0   !reset for next file


      !*************************************************************************
      ! Marker: End loop over files
      !*************************************************************************

      call chrono % tac(16)
      Orbital_Processing_Time_Minutes = chrono % get(16,minute=.true.)
      call chrono % summary(sort = .true.,title = 'GENERAL')
      call chrono % reset()

      call chronos_rttov % summary(title = 'RTTOV')
      call chronos_rttov % reset()

      call chronos_acha % summary(title = 'ACHA',sort=.true.)
      call chronos_acha % reset()

      !--- write algorithm attributes to level2
      call WRITE_ALGORITHM_ATTRIBUTES()

      !--- close pixel level hdf files
      call CLOSE_PIXEL_HDF_FILES(Level2_File_Flag)

   end do File_Loop

   !*************************************************************************
   ! Marker: Close FILELIST
   !*************************************************************************
   close(unit=File_List_Lun)

   !----------------------------------------------------------------------
   ! DEALLOCATE MEMORY AND END PROGRAM
   !----------------------------------------------------------------------

   !*************************************************************************
   !-- Marker: Close static ancillary data files
   !*************************************************************************

   !--- high resolution hdf ancillary data
   if (Sfc_Type_Id > 0) call CLOSE_LAND_SFC_HDF(Sfc_Type_Id)
   if (Coast_Mask_Id > 0)  call CLOSE_LAND_SFC_HDF(Coast_Mask_Id)
   if (Land_Mask_Id > 0)  call CLOSE_LAND_SFC_HDF(Land_Mask_Id)
   if (Volcano_Mask_Id > 0)  call CLOSE_LAND_SFC_HDF(Volcano_Mask_Id)
   if (Surface_Elev_Id > 0)  call CLOSE_LAND_SFC_HDF(Surface_Elev_Id)
   if (Snow_Mask_Id > 0)  call CLOSE_LAND_SFC_HDF(Snow_Mask_Id)
   if (Use_Sea_Ir_Emiss == sym%YES) call FORGET_SEA_IR_EMISS()
   if (Use_ABI_Dust == sym%YES) call FORGET_ABI_DUST()

   !--- remove directory for temporary files
   cmd = 'rmdir '//trim(Temporary_Data_Dir); nc = len_trim(cmd)
   call univ_system_cmd_f(nc, trim(cmd), ierror)

   !*************************************************************************
   ! Marker: Final remaining memory deallocation
   !*************************************************************************

   !--- Determine time of the start of all processing

   call chrono_all % tac(1)
   call chrono_all % summary(minute=.true.,title = 'ALL FILES')

   !---- print to screen that processing is done
   call MESG (  "<--------- End of CLAVRXORB ---------->",level=verb_lev % MINIMAL)

   stop

!*************************************************************************
! Marker: End of main code
!*************************************************************************

contains

!*************************************************************************
! Marker: Begin of Subroutines
!*************************************************************************

!*************************************************************************
! open modis white sky reflectance files
!*************************************************************************
subroutine OPEN_MODIS_WHITE_SKY_SFC_REFLECTANCE_FILES()

         !--- determine 16 day period and its string value
         iperiod16 = 16 * ((Image%Start_Doy-1) / 16) + 1
         write(Day_String,fmt="(i3.3)") iperiod16

         !--- Open Modis white sky 0.66 um albedo
         if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
            Modis_White_Sky_0_66_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".0.659_x4.hdf"
            Modis_Alb_0_66_Str%sds_Name = MODIS_ALB_0_66_SDS_NAME
            Modis_Alb_0_66_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        trim(Modis_White_Sky_0_66_Name), &
                                        grid_Str=Modis_Alb_0_66_Str)
         endif

         !--- Open Modis white sky 0.86 um albedo
         if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
            Modis_White_Sky_0_86_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".0.858_x4.hdf"
            Modis_Alb_0_86_Str%sds_Name = MODIS_ALB_0_86_SDS_NAME
            Modis_Alb_0_86_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        trim(Modis_White_Sky_0_86_Name), &
                                        grid_Str=Modis_Alb_0_86_Str)
         endif

         !--- Open Modis white sky 1.24 um albedo
         if (Sensor%Chan_On_Flag_Default(5) == sym%YES) then
            Modis_White_Sky_1_24_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".1.24_x4.hdf"
            Modis_Alb_1_24_Str%sds_Name = MODIS_ALB_1_24_SDS_NAME
            Modis_Alb_1_24_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        trim(Modis_White_Sky_1_24_Name), &
                                        grid_Str=Modis_Alb_1_24_Str)
         endif

         !--- Open Modis white sky 1.66 um albedo
         if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
            Modis_White_Sky_1_64_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".1.64_x4.hdf"
            Modis_Alb_1_64_Str%sds_Name = MODIS_ALB_1_64_SDS_NAME
            Modis_Alb_1_64_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        trim(Modis_White_Sky_1_64_Name), &
                                        grid_Str=Modis_Alb_1_64_Str)
         endif

         !--- Open Modis white sky 2.13 um albedo
         if (Sensor%Chan_On_Flag_Default(7) == sym%YES) then
            Modis_White_Sky_2_13_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".2.13_x4.hdf"
         Modis_Alb_2_13_Str%sds_Name = MODIS_ALB_2_13_SDS_NAME
         Modis_Alb_2_13_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        trim(Modis_White_Sky_2_13_Name), &
                                        grid_Str=Modis_Alb_2_13_Str)
         endif

  end subroutine OPEN_MODIS_WHITE_SKY_SFC_REFLECTANCE_FILES
  !****************************************************************************
  ! acquire the information for the snow classication from multiple sources
  !****************************************************************************
  subroutine GET_SNOW_MASK()

      if (Read_Snow_Mask == sym%READ_SNOW_HIRES) then
         Failed_IMS_Snow_Mask_Flag = sym%NO

         Snow_Mask_Str%sds_Name = SNOW_MASK_SDS_NAME
         Snow_Mask_File_Name = GET_SNOW_MAP_FILENAME(Image%Start_Year,Image%Start_Doy, &
                                                 trim(Snow_Data_Dir))
         if (trim(Snow_Mask_File_Name) == "no_file") then
            Failed_IMS_Snow_Mask_Flag = sym%YES
            call MESG ( "WARNING: Could not find Snow mask file ==> "// &
              Snow_Mask_File_Name, level = verb_lev % WARNING)
         else
            Snow_Mask_Id = OPEN_LAND_SFC_HDF(trim(Snow_Data_Dir), &
                                      Snow_Mask_File_Name, &
                                      grid_Str=Snow_Mask_Str)
            if (Snow_Mask_Id > 0) then
               Failed_IMS_Snow_Mask_Flag = sym%NO
               call MESG("Snow mask file opened successfully ")
            else
               Failed_IMS_Snow_Mask_Flag = sym%YES
               call MESG("WARNING: Snow mask file open failed ", level = verb_lev % WARNING)
            endif
         endif
      endif


      if (Read_Snow_Mask == sym%READ_SNOW_GLOB) THEN
         Failed_Glob_Snow_Mask_Flag = sym%NO

         Snow_Mask_File_Name = GET_GLOBSNOW_FILENAME(Image%Start_Year,Image%Start_Doy, &
                                                 trim(GlobSnow_Data_Dir))
         if (trim(Snow_Mask_File_Name) == "no_file") THEN
            Failed_Glob_Snow_Mask_Flag = sym%YES
            call MESG ( "WARNING:  Could not find GlobSnow mask file, using NWP ", level = verb_lev % WARNING )
         else
            CALL READ_GLOBSNOW_ANALYSIS_MAP(trim(GlobSnow_Data_Dir)//Snow_Mask_File_Name)
         endif
      endif

  end subroutine GET_SNOW_MASK
  !****************************************************************************
  ! open the hdf files that hold ancillary data
  !****************************************************************************
  subroutine OPEN_STATIC_ANCIL_FILES()
   !------------------------------------------------------------------------------
   !--- Read elevation data
   !------------------------------------------------------------------------------
   if (Read_Surface_Elevation /= 0) then
      call MESG  ( "Opening surface elevation file", level = verb_lev % VERBOSE)
      Surface_Elev_Str%sds_Name = SURFACE_ELEV_SDS_NAME

      !--- read in which elevation type that is specified.
      if (Read_Surface_Elevation == 1) then
         Surface_Elev_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"static/sfc_data/", &
                                      "GLOBE_1km_digelev.hdf", &
                                      grid_Str=Surface_Elev_Str)
      else ! low resolution, Read_Surface_Elevation = 2
         Surface_Elev_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"static/sfc_data/", &
                                     "GLOBE_8km_digelev.hdf", &
                                      grid_Str=Surface_Elev_Str)
      endif

   endif

   !------------------------------------------------------------------
   ! Open coast mask file
   !------------------------------------------------------------------
   if (Read_Coast_Mask == sym%YES) then
      call MESG ( "Opening coast file", level = verb_lev % VERBOSE)
      Coast_Mask_Str%sds_Name = COAST_MASK_SDS_NAME
      Coast_Mask_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"static/sfc_data/", &
                                      "coast_mask_1km.hdf", &
                                      grid_Str=Coast_Mask_Str)
   endif

   !------------------------------------------------------------------
   ! Open land surface type file
   !------------------------------------------------------------------
   call MESG ( "Opening land surface type file", level = verb_lev % VERBOSE)
   Sfc_Type_Str%sds_Name = SFC_TYPE_SDS_NAME

   if (Read_Hires_sfc_type == sym%YES) then
      Sfc_Type_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                     "gl-latlong-1km-landcover.hdf", &
                                      grid_Str=Sfc_Type_Str)
   else
      Sfc_Type_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                     "gl-latlong-8km-landcover.hdf", &
                                      grid_Str=Sfc_Type_Str)
   endif

   !------------------------------------------------------------------
   ! Open land mask file
   !------------------------------------------------------------------
   if (Read_Land_Mask == sym%YES) then
      call MESG  ( "Opening land mask file" ,level = verb_lev % VERBOSE )
      Land_Mask_Str%sds_Name = LAND_MASK_SDS_NAME
      Land_Mask_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                    "lw_geo_2001001_v03m.hdf", &
                                     grid_Str=Land_Mask_Str)
   endif

   !------------------------------------------------------------------
   ! Open volcano mask file
   !------------------------------------------------------------------
   if (Read_Volcano_Mask == sym%YES) then
      call MESG( "Opening volcano mask file",level = verb_lev % VERBOSE )
      Volcano_Mask_Str%sds_Name = VOLCANO_MASK_SDS_NAME
      Volcano_Mask_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        "volcano_mask_1km.hdf", &
                                        grid_Str=Volcano_Mask_Str)
   endif

  end subroutine OPEN_STATIC_ANCIL_FILES
  !*******************************************************************************************************
  ! Read in various ancillary data and map it to the pixel grid for each segment
  !*******************************************************************************************************
  subroutine  MAP_ANCIL_DATA_TO_PIXEL_GRID()

        ! - emiss surface is gone to cx_sfc_emissivity_mod.f90

        !--- mandatory fields - check for substitution of Bad_Pixel for space

        !--- surface type
        call READ_LAND_SFC_HDF(Sfc_Type_Id, Sfc_Type_Str &
            , Nav%Lat, Nav%Lon, Geo%Space_Mask, Sfc%Sfc_Type)

        !--- surface elevation
        if (Read_Surface_Elevation /= 0) then

               !--- read the high res data
               call READ_LAND_SFC_HDF(Surface_Elev_Id, Surface_Elev_Str &
                 , Nav%Lat, Nav%Lon, Geo%Space_Mask, Two_Byte_Temp)

               !---  convert to a real number
               Sfc%Zsfc_Hires = real(two_byte_temp,kind=real4)
               !--- values over water are missing, set to zero
               where(Sfc%Zsfc_Hires .EQfp. Missing_Value_Real4)
                  Sfc%Zsfc_Hires = 0.0
               end where

        endif

        !--- merge with nwp surface elevation
        if (NWP_PIX%Nwp_Opt /= 0) then
                call MERGE_NWP_HIRES_ZSFC(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
        endif

        !--- read coast mask
        if (Read_Coast_Mask == sym%YES) then
             call READ_LAND_SFC_HDF(Coast_Mask_Id, Coast_Mask_Str, Nav%Lat, Nav%Lon, Geo%Space_Mask, Sfc%Coast)
        endif

        !--- read land mask
        if (Read_Land_Mask == sym%YES) then
            call READ_LAND_SFC_HDF(Land_Mask_Id, Land_Mask_Str, Nav%Lat, Nav%Lon, Geo%Space_Mask, Sfc%Land)
        endif

        !--- modify land class with ndvi if available (helps mitigate navigation errors)
        call MODIFY_LAND_CLASS_WITH_NDVI(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

        !--- read volcano mask
        if (Read_Volcano_Mask == sym%YES) then
            call READ_LAND_SFC_HDF(Volcano_Mask_Id, Volcano_Mask_Str, Nav%Lat, Nav%Lon, Geo%Space_Mask, Sfc%Volcano_Mask)
        endif

        !--- read Snow mask
        if (Read_Snow_Mask == sym%READ_SNOW_HIRES .and. Failed_IMS_Snow_Mask_Flag == sym%NO) then
             call READ_LAND_SFC_HDF(Snow_Mask_Id, Snow_Mask_Str, Nav%Lat, Nav%Lon, Geo%Space_Mask, Sfc%Snow_IMS)
        endif

        if (Read_Snow_Mask == sym%READ_SNOW_GLOB .and. Failed_Glob_Snow_Mask_Flag == sym%NO ) then
            call GET_PIXEL_GLOBSNOW_ANALYSIS(Nav%Lat,Nav%Lon,Sfc%Land,Bad_Pixel_Mask,Sfc%Snow_Glob)
        endif

        !--- define binary land and coast masks (yes/no) from land and coast flags
        call COMPUTE_BINARY_LAND_COAST_MASKS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

        !--- interpolate sst analyses to each pixel
        if (Use_Sst_Anal == 1) then
               call GET_PIXEL_SST_ANALYSIS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
               call COMPUTE_SNOW_CLASS_OISST(SST_Anal_Cice,Sfc%Snow_OISST)
        else
               Sst_Anal = NWP_PIX%Tsfc
        endif

        !--- compute a coast mask relative to nwp data
        if (NWP_PIX%Nwp_Opt /= 0) then
            call COMPUTE_COAST_MASK_NWP(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
        endif

        !--- compute a snow classification from NWP
        if (NWP_PIX%Nwp_Opt /= 0) then
            call COMPUTE_SNOW_CLASS_NWP(NWP_PIX%Weasd, NWP_PIX%Sea_Ice_Frac,Sfc%Snow_NWP)
        endif

        !---- ensure missing values for space scenes
        where (Geo%Space_Mask )
               Sfc%Zsfc_Hires = Missing_Value_Real4
               Sfc%Coast = Missing_Value_Int1
               Sfc%Land = Missing_Value_Int1
               Sfc%Snow_IMS = Missing_Value_Int1
               Sfc%Snow_GLOB = Missing_Value_Int1
               Sfc%Snow_NWP = Missing_Value_Int1
               Sfc%Snow_OISST = Missing_Value_Int1
               Sfc%Volcano_Mask = Missing_Value_Int1
               Sfc%Sfc_Type = Missing_Value_Int1
        end where


        !--- interpolate white sky albedoes to each pixel in segment
        if (Modis_Clr_Alb_Flag == sym%YES) then

           if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
                if (.not. allocated(ch(1)%Sfc_Ref_White_Sky)) then
                     allocate(ch(1)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                     ch(1)%Sfc_Ref_White_Sky = Missing_Value_Real4
                endif
                call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_0_66_Id, Modis_Alb_0_66_Str, &
                                        ch(1)%Sfc_Ref_White_Sky)
           endif
           if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
                  if (.not. allocated(ch(2)%Sfc_Ref_White_Sky)) then
                     allocate(ch(2)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                     ch(2)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_0_86_Id, Modis_Alb_0_86_Str, &
                                          ch(2)%Sfc_Ref_White_Sky)
           endif
           if (Sensor%Chan_On_Flag_Default(5) == sym%YES) then
                  if (.not. allocated(ch(5)%Sfc_Ref_White_Sky)) then
                     allocate(ch(5)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                     ch(5)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_1_24_Id, Modis_Alb_1_24_Str, &
                                          ch(5)%Sfc_Ref_White_Sky)
           endif
           if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
                  if (.not. allocated(ch(6)%Sfc_Ref_White_Sky)) then
                     allocate(ch(6)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                     ch(6)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_1_64_Id, Modis_Alb_1_64_Str, &
                                          ch(6)%Sfc_Ref_White_Sky)
           endif
           if (Sensor%Chan_On_Flag_Default(7) == sym%YES) then
                  if (.not. allocated(ch(7)%Sfc_Ref_White_Sky)) then
                     allocate(ch(7)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                     ch(7)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_2_13_Id, Modis_Alb_2_13_Str, &
                                          ch(7)%Sfc_Ref_White_Sky)
           endif

        endif

    end subroutine  MAP_ANCIL_DATA_TO_PIXEL_GRID

!******************************************************************************
! End of Code
!******************************************************************************

 end program PROCESS_CLAVRX
