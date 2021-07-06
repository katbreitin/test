!$Id: sensor_mod.f90 3928 2020-07-26 15:41:26Z heidinger $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: sensor_mod.f90 (src)
!       SENSOR_MOD (program)
!
! PURPOSE: This module houses routines that apply to multiple sensors
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
! NOTES:
! 
! Routines in this module and their purpose:
!
!  the routines are called from process_clavrx.f90 in this order
!   file_loop
!      DETECT_SENSOR_FROM_FILE()
!      SET_FILE_DIMENSIONS()
!      SET_DATA_DATE_AND_TIME()
!      READ_INSTR_CONSTANTS()
!      READ_LEVEL1B_DATA()
!       ...
!       ... 
! end loop
!--------------------------------------------------------------------------------------
module SENSOR_MOD
   use PIXEL_COMMON_MOD, only: &
      AVHRR_Data_Type &
      , AVHRR_Ver_1b &
      , AVHRR_GAC_FLAG &
      , AVHRR_KLM_FLAG &
      , AVHRR_AAPP_FLAG &
      , AVHRR_IFF_Flag &
      , AVHRR_Fusion_Flag &
      , ABI_Use_104um_Flag &
      , Sc_Id_AVHRR &
      , Byte_Swap_1b &
      , AVHRR_1_Flag &
      , Cloud_Mask_Bayesian_Flag &
      , L1b_Rec_Length &
      , Dark_Composite_Name &
      , Two_Byte_Temp &
      , Ref_Ch1_Dark_Composite &
      , Therm_Cal_1b &
      , Static_Nav_File &
      , ch &
      , sensor &
      , image &
      , nav &
      , geo &
      , ancil_data_dir &
      , use_aux_flag &
      , cloud_mask_aux_read_flag
      
   use CALIBRATION_CONSTANTS_MOD,only: &  
      planck_a1,planck_a2 &
      , goes_input_time &
      , goes_epoch_time &
      , sat_name &
      , solar_ch20 &
      , ew_ch20 &
      , solar_ch20_nu &
      , planck_nu
   
   use ALGORITHM_CONSTANTS_MOD,only:
   
   use CONSTANTS_MOD,only:  &
      real4 &
      , int1 &
      , int4 &
      , int1 &
      , sym &
      , MISSING_VALUE_INT4  &
      , MISSING_VALUE_REAL4 &
      , NCHAN_CLAVRX &
      , DTOR &
      , EXE_PROMPT
   
   use FILE_TOOLS,only: &
      get_lun
   
   use AVHRR_MOD,only: &
      assign_avhrr_sat_id_num_internal &
      , define_1b_data &
      , determine_avhrr_file_type &
      , read_avhrr_instr_constants &
      , read_avhrr_level1b_data &
      , read_avhrr_level1b_header
   
   use GOES_MOD, only: &
         area_struct &
         , gvar_nav &
         , goes_xstride &
         , goes_sndr_xstride &
         , calibrate_goes_dark_composite &
         , get_goes_headers &
         , lmodel &
         , read_dark_composite_counts &
         , read_goes &
         , read_goes_instr_constants &
         , read_goes_sndr &
         , read_goes_sndr_instr_constants
   
   use MODIS_MOD, only : &
       DETERMINE_MODIS_CLOUD_MASK_FILE &
     , READ_MODIS_INSTR_CONSTANTS &
     , READ_MODIS &
     , READ_MODIS_SIZE_ATTR &
     , DETERMINE_MODIS_GEOLOCATION_FILE &
     , READ_MODIS_TIME_ATTR

   use FY2_MOD, only: &
      READ_FY &
    , READ_FY_INSTR_CONSTANTS

   use FY4_MOD, only: &
      READ_FY4_INSTR_CONSTANTS &
    , READ_FY4_DATE_TIME & 
    , READ_FY4_LEVEL1B_DATA

   use ABI_MOD, only: &
      READ_ABI &
    , READ_ABI_INSTR_CONSTANTS &
    , READ_NAVIGATION_BLOCK_ABI

   use COMS_MOD,only: &
    read_coms &
    , read_coms_instr_constants &
    , read_navigation_block_coms
   
   use HIRS_FUSION_MOD,only: READ_FUSION_HIRS_INSTR_CONSTANTS, &
                             READ_HIRS_DATA, &
                             REPLACE_AVHRR_WITH_HIRS
   
   use IFF_CLAVRX_BRIDGE , only : &
      READ_IFF_DATA &
      , READ_IFF_VIIRS_INSTR_CONSTANTS &
      , READ_IFF_AVHRR_INSTR_CONSTANTS &
      , READ_IFF_DATE_TIME &
      , GET_IFF_DIMS_BRIDGE
      
   use MTSAT_MOD, only: &
      mtsat_xstride &
      , calibrate_mtsat_dark_composite &
      , read_mtsat &
      , read_mtsat_instr_constants &
      , read_navigation_block_mtsat_fy
   
   use SEVIRI_MOD, only: &
      seviri_xstride &
      , calibrate_seviri_dark_composite &
      , read_msg_instr_constants &
      , read_navigation_block_seviri &
      , read_seviri
   
#ifdef HDF5LIBS
   use VIIRS_CLAVRX_BRIDGE , only : &
       READ_VIIRS_DATE_TIME &
       , READ_VIIRS_DATA &
       , GET_NUMBER_OF_SCANS_FROM_VIIRS_BRIDGE &
       , READ_VIIRS_INSTR_CONSTANTS
       
   use AHI_CLAVRX_BRIDGE, only: &
      read_ahi_data
   
   use VIIRS_NASA_READ_MODULE, only : &
       READ_VIIRS_NASA_DATE_TIME &
       , READ_VIIRS_NASA_DATA &
       , READ_NUMBER_OF_SCANS_VIIRS_NASA &
       , CHECK_IF_FUSION &
       , READ_FUSION_INSTR_CONSTANTS
  
  
   use VIIRS_NASA_HRES_READ_MOD, only : &
      READ_VIIRS_NASA_HRES_DATA &
      , viirs_nasa_hres_config_type
       
  
   use FY3D_READ_MODULE, only : &
       READ_FY3D_DATA &
       , READ_FY3D_DATE_TIME &
       , READ_NUMBER_OF_SCANS_FY3D &
       , READ_FY3D_INSTR_CONSTANTS

#endif

   use CLAVRX_MESSAGE_MOD,only: &
      mesg &
      , verb_lev
   
   use MVCM_READ_MOD, only: &
      determine_mvcm_name &
      ,read_mvcm_data
   
   use SAPF_READ_MOD, only: DETERMINE_SAPF_NAME, READ_SAPF_DATA 

   use AHI_MOD

   use CLAVRX_STATIC_NAV_MODULE, only: &
      SETUP_READ_LEVEL1B_FIXED_GRID_STATIC_NAV, &
      READ_HSD_FIXED_GRID_STATIC_NAV, &
      READ_LEVEL1B_FIXED_GRID_STATIC_NAV

   use CX_VGAC_MOD, only: &
      READ_NUMBER_OF_SCANS_VGAC, &
      READ_VGAC_DATE_TIME, &
      READ_VGAC_DATA

   use CX_EPS_SG_MOD, only: &
      READ_NUMBER_OF_SCANS_EPS_SG, &
      READ_EPS_SG_DATE_TIME, &
      READ_EPS_SG_INSTR_CONSTANTS, &
      READ_EPS_SG_DATA

   use CX_NETCDF4_MOD, only: &
      READ_NETCDF_GLOBAL_ATTRIBUTE &
      , READ_NETCDF_DIMENSION
   
    use iso_fortran_env
   
   implicit none

   private
   public :: SET_DATA_DATE_AND_TIME   
   public :: READ_INSTR_CONSTANTS
   public :: DETECT_SENSOR_FROM_FILE
   public :: SET_FILE_DIMENSIONS
   public :: READ_LEVEL1B_DATA 
   public :: OUTPUT_SENSOR_TO_SCREEN
   public :: OUTPUT_IMAGE_TO_SCREEN
   public :: OUTPUT_PROCESSING_LIMITS_TO_SCREEN
   public :: SET_SENSOR_CHANNEL_MAPPING
   private :: COMPUTE_GOES_RU_STATIC_NAV_FILE_NAME

   

   character(24), parameter, private :: MODULE_PROMPT = " SENSOR_MODULE: "
   character(38) :: Orbit_Identifier
  
   character (len = 3), private :: string_3
   
   character (len = 6), private :: string_6
   character (len = 7), private :: string_7
   character (len = 30), private :: string_30

   contains

   !==============================================================================
   !   Determine date and time of the data and store in image structure
   !==============================================================================
   subroutine SET_DATA_DATE_AND_TIME(AREAstr)

      use class_time_date, only: date_type
      
      use CX_READ_AHI_MOD, only: &
         ahi_time_from_filename, &
         ahi_hcast_time_from_filename
      
      type ( date_type ) :: time0_obj, time1_obj
      
      ! - this is only needed/used for AVHRR
      type (AREA_STRUCT), intent(in) :: AREAstr

      integer(kind=int32):: Start_Year_Tmp
      integer(kind=int4):: Start_Day_Tmp
      integer(kind=int4):: End_Year_Tmp
      integer(kind=int4):: End_Day_Tmp
      integer(kind=int4):: Start_Time_Tmp
      integer(kind=int4):: End_Time_Tmp
      integer(kind=int4):: Orbit_Number_Tmp
      integer(kind=int4):: Hour
      integer(kind=int4):: Minute
      integer(kind=int4):: Second
      integer :: year, doy

      !--- Static Navigation additions
      type ( date_type ) :: abi_time_start
      type ( date_type ) :: abi_time_end
      character(len=1020):: L1b_Full_File_Name
      character(len=1020):: Start_Year_Mon_Day_HH_MM_Tmp
      character(len=1020):: End_Year_Mon_Day_HH_MM_Tmp
      character(len=4):: Start_Year_Tmp_Str
      character(len=2):: Start_Mon_Tmp_Str
      character(len=2):: Start_DD_Tmp_Str
      character(len=2):: Start_HH_Tmp_Str
      character(len=2):: Start_Min_Tmp_Str
      character(len=2):: Start_Sec_Tmp_Str
      character(len=1):: Start_Msec_Tmp_Str
      integer(kind=int4):: Start_Mon_Tmp
      integer(kind=int4):: Start_DD_Tmp
      integer(kind=int4):: Start_HH_Tmp
      integer(kind=int4):: Start_Min_Tmp
      integer(kind=int4):: Start_Sec_Tmp
      integer(kind=int4):: Start_Msec_Tmp
      character(len=4):: End_Year_Tmp_Str
      character(len=2):: End_Mon_Tmp_Str
      character(len=2):: End_DD_Tmp_Str
      character(len=2):: End_HH_Tmp_Str
      character(len=2):: End_Min_Tmp_Str
      character(len=2):: End_Sec_Tmp_Str
      character(len=1):: End_Msec_Tmp_Str
      integer(kind=int4):: End_Mon_Tmp
      integer(kind=int4):: End_DD_Tmp
      integer(kind=int4):: End_HH_Tmp
      integer(kind=int4):: End_Min_Tmp
      integer(kind=int4):: End_Sec_Tmp
      integer(kind=int4):: End_Msec_Tmp

      !----------------------------------------------
      ! for AVHRR, this is read in with level-1b data
      !----------------------------------------------

      !----------------------------------------------
      ! for Modis, take time from file name
      !----------------------------------------------
      if (index(Sensor%Sensor_Name,'MODIS') > 0) then

         call READ_MODIS_TIME_ATTR(trim(Image%Level1b_Path), trim(Image%Level1b_Name), &
                            Image%Start_Year, Image%Start_Doy, Image%Start_Time, &
                            Image%End_Year, Image%End_Doy, Image%End_Time)
  
      end if

      !----------------------------------------------
      ! read VIIRS-NOAA time from GMTCO
      !----------------------------------------------
      if (trim(Sensor%Sensor_Name) == 'VIIRS') then
#ifdef HDF5LIBS
         call READ_VIIRS_DATE_TIME(trim(Image%Level1b_Path),trim(Image%Level1b_Name), &
                             Start_Year_Tmp,Start_Day_Tmp,Start_Time_Tmp, &
                             End_Time_Tmp,Orbit_Number_Tmp,Orbit_Identifier, &
                             End_Year_Tmp , End_Day_Tmp)
         Image%Start_Year = Start_Year_Tmp
         Image%End_Year = End_Year_Tmp
         Image%Start_Doy = Start_Day_Tmp
         Image%End_Doy = End_Day_Tmp
         Image%Orbit_Number = Orbit_Number_Tmp
 
         Image%Start_Time = Start_Time_Tmp
         Image%End_Time = End_Time_Tmp
#else
         PRINT *, "No HDF5 libraries installed, stopping"
         stop
#endif
      end if 
      
      if (trim(Sensor%Sensor_Name) == 'VIIRS-NASA-HRES') then
          print*
          print*,'  +++++++++++++++++++++++++++++++++++  +++++++++++++++++++++++++++++++++++++='
          print*,'read date time has to be written is fake... READ_VIIRS_NASA_DATE_TIME..  File: ', __FILE__,' Line: ',__LINE__
         Image%Start_Year = 2020
         Image%End_Year = 2020
         Image%Start_Doy = 118
         Image%End_Doy = 118
         Image%Orbit_Number = 7889887

         Image%Start_Time = 0.5
         Image%End_Time = 0.9
          print*,'  +++++++++++++++++++++++++++++++++++  ++++++++++++++++++++++++++++++++++++++++='
          
          
         
      end if 
       
        
      !----------------------------------------------
      ! read VIIRS-NASA time from VGEOM
      !----------------------------------------------
      if (trim(Sensor%Sensor_Name) == 'VIIRS-NASA') then
#ifdef HDF5LIBS
         call READ_VIIRS_NASA_DATE_TIME(trim(Image%Level1b_Path),trim(Image%Level1b_Name), &
                             Start_Year_Tmp,Start_Day_Tmp,Start_Time_Tmp, &
                             End_Time_Tmp,Orbit_Number_Tmp, End_Year_Tmp , End_Day_Tmp)
         Image%Start_Year = Start_Year_Tmp
         Image%End_Year = End_Year_Tmp
         Image%Start_Doy = Start_Day_Tmp
         Image%End_Doy = End_Day_Tmp
         Image%Orbit_Number = Orbit_Number_Tmp

         Image%Start_Time = Start_Time_Tmp
         Image%End_Time = End_Time_Tmp

#else
         PRINT *, "No HDF5 libraries installed, stopping"
         stop
#endif
      endif

      !----------------------------------------------
      ! 
      !----------------------------------------------
      !if (index(Sensor%Sensor_Name,'AHI') > 0 .OR. index(Sensor%Sensor_Name,'AHI9') > 0) then
      if (index(Sensor%Sensor_Name,'AHI') > 0) then
         
         
         ! Needed since HCAST native is not named the same - WCS3
         IF ((Image%DB_Flag) .AND. (Image%Nc_Format_Flag .EQV. .false.)) THEN 
            call ahi_hcast_time_from_filename( trim(Image%Level1b_Name) , &
                                                time0_obj, time1_obj )
         else
         
            call ahi_time_from_filename ( trim(Image%Level1b_Name) , &
                                          time0_obj, time1_obj )
         endif
        
         call time0_obj % get_date ( year =  year &
                               , doy = doy  &
                               , msec_of_day = Image%Start_Time  )
         
         call time1_obj % get_date ( msec_of_day = Image%End_Time  )                                                

         Image%Start_Year  = year
         Image%Start_Doy   = doy   
         Image%End_Year  = year
         Image%End_Doy   = doy  

      endif
      
      !----------------------------------------------
      ! --- VIIRS Subsampled files
      !----------------------------------------------
      if (index(Sensor%Sensor_Name,'VGAC') > 0) then
         call READ_VGAC_DATE_TIME(Start_Year_Tmp,Start_Day_Tmp,Start_Time_Tmp,&
                                  End_Year_Tmp,End_Day_Tmp,End_Time_Tmp)
         Image%Start_Year = Start_Year_Tmp
         Image%End_Year = End_Year_Tmp
         Image%Start_Doy = Start_Day_Tmp
         Image%End_Doy = End_Day_Tmp
         Image%Start_Time = Start_Time_Tmp
         Image%End_Time = End_Time_Tmp
      endif

      !----------------------------------------------
      ! --- EUMETSAT Polar System - Secong Generation (EPS SG)
      !----------------------------------------------
      if (index(Sensor%Sensor_Name,'METIMAGE') > 0) then
         call READ_EPS_SG_DATE_TIME(trim(Image%Level1b_Path)//trim(Image%Level1b_Name),&
                                  Start_Year_Tmp,Start_Day_Tmp,Start_Time_Tmp,&
                                  End_Year_Tmp,End_Day_Tmp,End_Time_Tmp)
         Image%Start_Year = Start_Year_Tmp
         Image%End_Year = End_Year_Tmp
         Image%Start_Doy = Start_Day_Tmp
         Image%End_Doy = End_Day_Tmp
         Image%Start_Time = Start_Time_Tmp
         Image%End_Time = End_Time_Tmp
      endif

      !----------------------------------------------
      ! for IFF take time and set some constants
      ! could be VIIRS, MODIS AVHRR sensor
      !----------------------------------------------
      if (index(Sensor%Sensor_Name,'IFF') > 0) then
         call READ_IFF_DATE_TIME(trim(Image%Level1b_Path), trim(Image%Level1b_Name),Start_Year_Tmp, &
                      Start_Day_Tmp,Start_Time_Tmp, End_Year_Tmp,End_Day_Tmp,End_Time_Tmp)
         Image%Start_Year = Start_Year_Tmp
         Image%End_Year = End_Year_Tmp
         Image%Start_Doy = Start_Day_Tmp
         Image%End_Doy = End_Day_Tmp
         Image%Start_Time = Start_Time_Tmp
         Image%End_Time = End_Time_Tmp
      endif

       !----------------------------------------------
       ! read time from FY3D
       !----------------------------------------------
       if (trim(Sensor%Sensor_Name) == 'MERSI-2') then
#ifdef HDF5LIBS
          call READ_FY3D_DATE_TIME(trim(Image%Level1b_Path),trim(Image%Level1b_Name), &
                              Start_Year_Tmp,Start_Day_Tmp,Start_Time_Tmp, &
                              End_Time_Tmp,Orbit_Number_Tmp, End_Year_Tmp , End_Day_Tmp)
          Image%Start_Year = Start_Year_Tmp
          Image%End_Year = End_Year_Tmp
          Image%Start_Doy = Start_Day_Tmp
          Image%End_Doy = End_Day_Tmp
          Image%Orbit_Number = Orbit_Number_Tmp

          Image%Start_Time = Start_Time_Tmp
          Image%End_Time = End_Time_Tmp

#else
          PRINT *, "No HDF5 libraries installed, stopping"
          stop
#endif
       endif


      !----------------------------------------------
      ! for FY-4A sensor AGRI
      !----------------------------------------------
      if (index(Sensor%Sensor_Name,'AGRI') > 0) then
         call READ_FY4_DATE_TIME(trim(Image%Level1b_Path),trim(Image%Level1b_Name), &
                             Start_Year_Tmp,Start_Day_Tmp,Start_Time_Tmp, &
                             End_Year_Tmp, End_Day_Tmp, End_Time_Tmp)
         Image%Start_Year = Start_Year_Tmp
         Image%End_Year = End_Year_Tmp
         Image%Start_Doy = Start_Day_Tmp
         Image%End_Doy = End_Day_Tmp
         Image%Start_Time = Start_Time_Tmp
         Image%End_Time = End_Time_Tmp
      endif


      !----------------------------------------------
      ! for GOES, MTSAT and MSG, take time from AREAstr
      !----------------------------------------------

      !-- initialize this to false, set true below
      if (index(Sensor%Sensor_Name,'GOES') > 0 .or.  &
          index(Sensor%Sensor_Name,'COMS') > 0 .or.  &
          index(Sensor%Sensor_Name,'MTSAT') > 0 .or.  &
          index(Sensor%Sensor_Name,'SEVIRI') > 0 .or.  &
          index(Sensor%Sensor_Name,'FY2') > 0) then

        !--- check if area file
        if (AREAstr%Version_Num == 4) then

         Image%Start_Year = 1900 + int(AREAstr%img_Date / 1000)
         Image%End_Year = Image%Start_Year
         Image%Start_Doy = AREAstr%img_Date - (Image%Start_Year - 1900) * 1000
         Image%End_Doy = Image%Start_Doy
         hour = AREAstr%img_Time / 10000 
         minute = (AREAstr%img_Time - hour * 10000) / 100
         second = (AREAstr%img_Time - hour * 10000 - minute * 100) / 100
         Image%Start_Time = ((hour * 60 + minute) * 60 + second) * 1000 !millisec
         Image%End_Time = Image%Start_Time

        endif

        if (.not. Image%Area_Format_Flag) then

          Image%Mixed_Resolution_Flag = .true.
          Image%Static_Nav_Flag = .true. 

          !--- check if GOES ABI static navigation.
          !--- 2018-05-14T18:30:40.5Z
          L1b_Full_File_Name = trim(Image%Level1b_Path)//trim(Image%Level1b_Name)

          !--- Image Start
          call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(L1b_Full_File_Name), 'time_coverage_start', Start_Year_Mon_Day_HH_MM_Tmp)
          Start_Year_Tmp_Str = Start_Year_Mon_Day_HH_MM_Tmp(1:4)
          read(Start_Year_Tmp_Str,*) Start_Year_Tmp                 ! convert to integer
          Start_Mon_Tmp_Str = Start_Year_Mon_Day_HH_MM_Tmp(6:7)
          read(Start_Mon_Tmp_Str,*) Start_Mon_Tmp
          Start_DD_Tmp_Str = Start_Year_Mon_Day_HH_MM_Tmp(9:10)
          read(Start_DD_Tmp_Str,*) Start_DD_Tmp
          Start_HH_Tmp_Str = Start_Year_Mon_Day_HH_MM_Tmp(12:13)
          read(Start_HH_Tmp_Str,*) Start_HH_Tmp
          Start_Min_Tmp_Str = Start_Year_Mon_Day_HH_MM_Tmp(15:16)
          read(Start_Min_Tmp_Str,*) Start_Min_Tmp
          Start_Sec_Tmp_Str = Start_Year_Mon_Day_HH_MM_Tmp(18:19)
          read(Start_Sec_Tmp_Str,*) Start_Sec_Tmp
          Start_Msec_Tmp_Str = Start_Year_Mon_Day_HH_MM_Tmp(21:21)
          read(Start_Msec_Tmp_Str,*) Start_Msec_Tmp
          call abi_time_start % set_date(Start_Year_Tmp, Start_Mon_Tmp, Start_DD_Tmp, Start_HH_Tmp, Start_Min_Tmp, Start_Sec_Tmp)
          Image%Start_Year = abi_time_start%year
          Image%Start_Doy = abi_time_start%dayOfYear
          Image%Start_Time = abi_time_start%msec_of_day !millisec

          !--- Image End
          call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(L1b_Full_File_Name), 'time_coverage_end', End_Year_Mon_Day_HH_MM_Tmp)
          End_Year_Tmp_Str = End_Year_Mon_Day_HH_MM_Tmp(1:4)
          read(End_Year_Tmp_Str,*) End_Year_Tmp                 ! convert to integer
          End_Mon_Tmp_Str = End_Year_Mon_Day_HH_MM_Tmp(6:7)
          read(End_Mon_Tmp_Str,*) End_Mon_Tmp
          End_DD_Tmp_Str = End_Year_Mon_Day_HH_MM_Tmp(9:10)
          read(End_DD_Tmp_Str,*) End_DD_Tmp
          End_HH_Tmp_Str = End_Year_Mon_Day_HH_MM_Tmp(12:13)
          read(End_HH_Tmp_Str,*) End_HH_Tmp
          End_Min_Tmp_Str = End_Year_Mon_Day_HH_MM_Tmp(15:16)
          read(End_Min_Tmp_Str,*) End_Min_Tmp
          End_Sec_Tmp_Str = End_Year_Mon_Day_HH_MM_Tmp(18:19)
          read(End_Sec_Tmp_Str,*) End_Sec_Tmp
          End_Msec_Tmp_Str = End_Year_Mon_Day_HH_MM_Tmp(21:21)
          read(End_Msec_Tmp_Str,*) End_Msec_Tmp
          call abi_time_end % set_date(End_Year_Tmp, End_Mon_Tmp, End_DD_Tmp, End_HH_Tmp, End_Min_Tmp, End_Sec_Tmp)
          Image%End_Year = abi_time_end%year
          Image%End_Doy = abi_time_end%dayOfYear
          Image%End_Time = abi_time_end%msec_of_day !millisec

        endif
 
      end if

      !-------------------------------------------------------------------------
      ! compute start, end and mean time in units of hours
      !-------------------------------------------------------------------------
      Image%Start_Time_Hours = Image%Start_Time / 60.0 / 60.0 / 1000.0
      Image%End_Time_Hours = Image%End_Time / 60.0 / 60.0 / 1000.0
      Image%Mean_Time_Hours = (Image%Start_Time_Hours + Image%End_Time_Hours) / 2.0
      if (Image%End_Time_Hours < Image%Start_Time_Hours) then
        Image%Mean_Time_Hours = (Image%Start_Time_Hours + Image%End_Time_Hours + 24.0) / 2.0
      endif
      if (Image%Mean_Time_Hours > 24.0) then
        Image%Mean_Time_Hours = Image%Mean_Time_Hours - 24.0
      endif


      !----------------------------------------------------------------------------
      ! Replace with a date_type variables from class_time_date.f90
      !----------------------------------------------------------------------------

   end subroutine SET_DATA_DATE_AND_TIME
   !--------------------------------------------------------------------------------------
   !   screen output of sensor structure
   !--------------------------------------------------------------------------------------
   subroutine OUTPUT_SENSOR_TO_SCREEN()

      call mesg ( " ",level = verb_lev % DEFAULT)
      call mesg ( "SENSOR DEFINITION",level = verb_lev % DEFAULT)

      call mesg ( "Satellite = "//trim(Sensor%Platform_Name), level = verb_lev % DEFAULT)
      call mesg ( "Sensor = "//trim(Sensor%Sensor_Name), level = verb_lev % DEFAULT)

      write(string_3,'(i3)' ) Sensor%WMO_ID 
      call mesg ( "Spacecraft WMO number = "//trim(string_3) , level = verb_lev % DEFAULT)

      write ( string_6,'(i6)')   Sensor%Spatial_Resolution_Meters
      call mesg ( "Pixel Resolution (m) = "//string_6, level = verb_lev % DEFAULT)

      !--- some avhrr specific output
      if (index(Sensor%Sensor_Name,'AVHRR-1') > 0 .or. &
          index(Sensor%Sensor_Name,'AVHRR-2') > 0 .or. &
          index(Sensor%Sensor_Name,'AVHRR-3') > 0 .or. &
          index(Sensor%Sensor_Name,'AVHRR-FUSION') > 0) then

         write ( string_3,'(i3)')   AVHRR_Data_Type
         call mesg ( "AVHRR data type = "//string_3, level = verb_lev % DEFAULT)

         write ( string_3,'(i3)')   AVHRR_Ver_1b
         call mesg ( "AVHRR Level1b Version = "//string_3, level = verb_lev % DEFAULT)

         write ( string_3,'(i3)')   AVHRR_GAC_FLAG
         call mesg ( "AVHRR GAC Flag = "//string_3, level = verb_lev % DEFAULT)

         write ( string_3,'(i3)')   AVHRR_KLM_FLAG
         call mesg ( "AVHRR KLM Flag = "//string_3, level = verb_lev % DEFAULT)

         write ( string_3,'(i3)')   AVHRR_AAPP_FLAG
         call mesg ( "AVHRR AAPP Flag = "//string_3, level = verb_lev % DEFAULT)
    
      end if

   end subroutine OUTPUT_SENSOR_TO_SCREEN
   !--------------------------------------------------------------------------------------
   !   screen output of some members of the image structure
   !--------------------------------------------------------------------------------------
   subroutine OUTPUT_IMAGE_TO_SCREEN()

      call mesg ( " ",level = verb_lev % DEFAULT)
      call mesg ( "IMAGE DEFINITION",level = verb_lev % DEFAULT)

      call mesg ("Level1b Name = "//trim(Image%Level1b_Name) , level = verb_lev % MINIMAL )

      write(string_6,'(i6)' ) Image%Number_Of_Elements
      call mesg ( "Number of Elements Per Line = "//string_6,level = verb_lev % DEFAULT)

      write(string_6,'(i6)' ) Image%Number_Of_Lines
      call mesg ( "Number of Lines in File = "//string_6,level = verb_lev % DEFAULT)

      write(string_6,'(i6)' ) Image%Number_Of_Lines_Per_Segment
      call mesg ( "Number of Lines in each Segment = "//string_6,level = verb_lev % DEFAULT)

      write ( string_30, '(I4,1X,I3,1X,F9.5)') Image%Start_Year,Image%Start_Doy, &
             Image%Start_Time/60.0/60.0/1000.0
      call mesg ("Start Year, Doy, Time = "//string_30,level = verb_lev % DEFAULT)

      write ( string_30, '(I4,1X,I3,1X,F9.5)') Image%End_Year,Image%End_Doy, &
             Image%End_Time/60.0/60.0/1000.0
      call mesg ("End Year, Doy, Time = "//string_30,level = verb_lev % DEFAULT)

      if (.not. Image%Area_Format_Flag) call mesg("Area Format Flag = F",level = verb_lev % DEFAULT)
      if (Image%Area_Format_Flag) call mesg("Area Format Flag = T",level = verb_lev % DEFAULT)

      if (.not. Image%Static_Nav_Flag) call mesg("Static Nav Flag = F",level = verb_lev % DEFAULT)
      if (Image%Static_Nav_Flag) call mesg("Static Nav Flag = T",level = verb_lev % DEFAULT)

      if (.not. Image%Mixed_Resolution_Flag) call mesg("Mixed Resolution Flag = F",level = verb_lev % DEFAULT)
      if (Image%Mixed_Resolution_Flag) call mesg("Mixed Resolution Flag = T",level = verb_lev % DEFAULT)

      if (.not. Image%DB_Flag) call mesg("Direct Broadcast Flag = F",level = verb_lev % DEFAULT)
      if (Image%DB_Flag) call mesg("Direct Broadcast Flag = T",level = verb_lev % DEFAULT)

      call mesg ( " ",level = verb_lev % DEFAULT)

   end subroutine OUTPUT_IMAGE_TO_SCREEN

   !--------------------------------------------------------------------------------------
   !   screen output of the spatial or viewing limits imposed on this processing
   !--------------------------------------------------------------------------------------
   subroutine OUTPUT_PROCESSING_LIMITS_TO_SCREEN()
      call mesg ( "PROCESSING LIMITS",level = verb_lev % DEFAULT)

      write(string_7,'(f7.2)' ) Nav%Lon_East_Limit
      call mesg ( "Eastern Longitude for Processing = "//string_7, level = verb_lev % DEFAULT)

      write(string_7,'(f7.2)' ) Nav%Lon_West_Limit
      call mesg ( "Western Longitude for Processing = "//string_7, level = verb_lev % DEFAULT)

      write(string_7,'(f7.1)' ) Nav%Lat_North_Limit
      call mesg ( "Northern Latitude for Processing = "//string_7, level = verb_lev % DEFAULT)

      write(string_7,'(f7.1)' ) Nav%Lat_South_Limit
      call mesg ( "Southern Latitude for Processing = "//string_7, level = verb_lev % DEFAULT)

      write(string_7,'(f7.1)' ) Geo%Satzen_Max_Limit
      call mesg ( "Maximum Sensor Zenith Angle for Processing = "//string_7, level = verb_lev % DEFAULT)

      write(string_7,'(f7.1)' ) Geo%Satzen_Min_Limit
      call mesg ( "Minimum Sensor Zenith Angle for Processing = "//string_7, level = verb_lev % DEFAULT)

      write(string_7,'(f7.1)' ) Geo%Solzen_Max_Limit
      call mesg ( "Maximum Solar Zenith Angle for Processing = "//string_7, level = verb_lev % DEFAULT)

      write(string_7,'(f7.1)' ) Geo%Solzen_Min_Limit
      call mesg ( "Minimum Solar Zenith Angle for Processing = "//string_7, level = verb_lev % DEFAULT)

      call mesg ( " ",level = verb_lev % DEFAULT)

   end subroutine OUTPUT_PROCESSING_LIMITS_TO_SCREEN

   !--------------------------------------------------------------------------------------------------
   !  Read the values from instrument constant files
   !--------------------------------------------------------------------------------------------------
   subroutine READ_INSTR_CONSTANTS()

      AVHRR_IFF_Flag = 0

      select case(trim(Sensor%Sensor_Name))

         case('AVHRR-1','AVHRR-2','AVHRR-3')
              call READ_AVHRR_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('MODIS','MODIS-MAC','MODIS-CSPP','AQUA-IFF')
              call READ_MODIS_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('GOES-IL-IMAGER','GOES-MP-IMAGER')
              call READ_GOES_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('GOES-IP-SOUNDER')
              call READ_GOES_SNDR_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('GOES-RU-IMAGER')
           call READ_ABI_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('SEVIRI')
              call READ_MSG_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('MTSAT-IMAGER')
              call READ_MTSAT_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('FY2-IMAGER')
              call READ_FY_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('AGRI')
              call READ_FY4_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('MERSI-2')
                call READ_FY3D_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('COMS-IMAGER')
              call READ_COMS_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('METIMAGE')
              call READ_EPS_SG_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('AHI')
              call READ_AHI_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
!        case('AHI9')
!             call READ_AHI_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('VIIRS','VIIRS-NASA','VGAC','VIIRS-NASA-HRES')
              if (.not. Sensor%Fusion_Flag) then
                call READ_VIIRS_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
              else
                call READ_FUSION_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
              endif
         case('VIIRS-IFF')
            call READ_IFF_VIIRS_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         case('AVHRR-IFF')
            call READ_IFF_AVHRR_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
            AVHRR_IFF_Flag = 1
         case('AVHRR-FUSION')
            call READ_AVHRR_INSTR_CONSTANTS(trim(Sensor%Instr_Const_File))
         if (index(trim(Sensor%Instr_Const_File),'avhrr_5_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_5_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_6_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_6_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_7_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_7_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_8_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_8_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_9_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_9_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_10_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_10_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_11_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_11_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_12_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_12_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_14_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_14_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_15_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_15_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_16_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_16_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_17_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_17_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_18_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_18_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_19_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_19_instr.dat'))            
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_2_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_2_instr.dat'))
         endif
         if (index(trim(Sensor%Instr_Const_File),'avhrr_1_instr.dat') > 0) then
            call READ_FUSION_HIRS_INSTR_CONSTANTS(trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim('hirs_1_instr.dat'))
         endif
      end select

   end subroutine READ_INSTR_CONSTANTS

   !--------------------------------------------------------------------------------------------------
   !  Apply various tests to determine from which sensor this data comes
   !
   !  output is sesnorname and platform name
   !
   !
   !   AVHRR 
   !        sensors:  AVHRR-1, AVHRR-2, AVHRR-3
   !        platforms: NOAA-5 - NOAA-19, METOP-A - METOP-C
   !        spatial_resolution:  1.1, 4
   !  
   !   GOES
   !        sensors:   GOES-IL-IMAGER, GOES-MP-IMAGER, GOES-IP-SOUNDER
   !        platforms: GOES-8 - GOES-15
   !        spatial_resolution:  1, 4
   !
   !   METEOSAT
   !        sensors:  SEVIRI
   !        platform:  Meteosat-8 - Meteosat-11
   !        spatial_resolution:  3
   ! 
   !   MTSAT
   !        sensors:  MTSAT-IMAGER
   !        platform:  MTSAT-1R, MTSAT-2
   !        spatial_resolution:  4
   !
   !   COMS
   !        sensors:  COMS-IMAGER
   !        platform:  COMS-1
   !        spatial_resolution:  4
   !
   !   MODIS
   !        sensors:  MODIS, AQUA-IFF, MODIS-MAC, MODIS-CSPP
   !        platform:  AQUA, TERRA
   !        spatial_resolution:  1, 5
   !
   !   VIIRS
   !        sensors:  VIIRS, VIIRS-NASA, VIIRS-IFF
   !        platform:  SNPP
   !        spatial_resolution:  0.75
   !        
   !
   !  Output: The sensor structure which is global does not appear as an argument
   !--------------------------------------------------------------------------------------------------
   subroutine DETECT_SENSOR_FROM_FILE( &
           AREAstr &
         , NAVstr &
         , Ierror)

      !use CX_HDF4_MOD, only &
      !   OPEN_FILE_HDF_READ, &
      !   READ_HDF_GLOBAL_ATTRIBUTE_NUM, &
      !   CLOSE_FILE_HDF_READ

      implicit none

      TYPE (AREA_STRUCT), intent(out) :: AREAstr
      TYPE (GVAR_NAV), intent(out)    :: NAVstr
      integer(kind=int4) :: Ierror
      integer(kind=int4) :: Ifound

      !-------------------------------------------------------------------------
      !-- Initialize Sensor Structure
      !-------------------------------------------------------------------------
      Sensor%Sensor_Name = ''
      Sensor%Spatial_Resolution_Meters = Missing_Value_Int4
      Sensor%Platform_Name = ''
      Sensor%WMO_Id = Missing_Value_Int4
      Sensor%Instr_Const_File = 'no_file'
      Sensor%Geo_Sub_Satellite_Longitude = Missing_Value_Real4
      Sensor%Geo_Sub_Satellite_Latitude = Missing_Value_Real4
      Image%Area_Format_Flag = .false.
      Image%Mixed_Resolution_Flag = .false.
      Image%Static_Nav_Flag = .false.
      Image%DB_Flag = .false.

      !-------------------------------------------------------------------------
      !-- Loop through tests for each file type, if not found assume avhrr
      !-------------------------------------------------------------------------
      Ierror = sym%NO
      ifound = sym%NO

      test_loop: do while (ifound == sym%NO)
      
      !--- HIMAWARI-8 AHI Test
      if (index(Image%Level1b_Name, 'HS_H08') > 0) then

        Sensor%Sensor_Name = 'AHI'
        Sensor%Platform_Name = 'HIM8'
        Sensor%Spatial_Resolution_Meters = 2000
        Sensor%WMO_Id = 173
        Sensor%Instr_Const_File = 'him8_instr.dat'
        Sensor%Geo_Sub_Satellite_Longitude = 140.0   ! 140.66
        Sensor%Geo_Sub_Satellite_Latitude = 0.0
        Static_Nav_File = 'himawari_central_ahi_fulldisk_static_nav.nc'
        Image%Area_Format_Flag = .false.
        Image%Nc_Format_Flag = .false.
              
        !Done for netCDF files - WCS3
        if (index(Image%Level1b_Name, '.nc') > 1) then
            Image%Nc_Format_Flag = .true.
            call READ_NETCDF_DIMENSION(trim(Image%Level1b_Full_Name), 'x', Image%Number_Of_Elements)
            call READ_NETCDF_DIMENSION(trim(Image%Level1b_Full_Name), 'y', Image%Number_Of_Lines)
        endif
        
        !Done for HSD files - WCS3
        if (index(Image%Level1b_Name, '.nc') < 1) then
#ifdef LIBHIM
            Image%Nc_Format_Flag = .false.
            !For HSD, we will always do static nav. Because of that, we don't need 
            ! Image%Number_Of_Elements and Image%Number_Of_Lines. These are determined 
            ! in DETERMINE_BOUNDS_STATIC_NAV for a given segment of data 
#else
            call MESG( "LibHimawari not installed. Cannot process HSD. Stopping", level = verb_lev % ERROR , color = 4 )
            stop
#endif
        endif

        !END - WCS3

        Image%Mixed_Resolution_Flag = .true.
        Image%Static_Nav_Flag = .true.
        
        exit test_loop

      endif

      !--- HIMAWARI-9 AHI Test
      if (index(Image%Level1b_Name, 'HS_H09') > 0) then
        Sensor%Sensor_Name = 'AHI'
        Sensor%Platform_Name = 'HIM9'
        Sensor%Spatial_Resolution_Meters = 2000
        Sensor%WMO_Id = 174
        Sensor%Instr_Const_File = 'him9_instr.dat'
        Sensor%Geo_Sub_Satellite_Longitude = 140.0   ! 140.72
        Sensor%Geo_Sub_Satellite_Latitude = 0.0
        Static_Nav_File = 'himawari_central_ahi_fulldisk_static_nav.nc'
        Image%Area_Format_Flag = .false.
        Image%Nc_Format_Flag = .false.
              
        !Done for netCDF files - WCS3
        if (index(Image%Level1b_Name, '.nc') > 1) then
            Image%Nc_Format_Flag = .true.
            call READ_NETCDF_DIMENSION(trim(Image%Level1b_Full_Name), 'x', Image%Number_Of_Elements)
            call READ_NETCDF_DIMENSION(trim(Image%Level1b_Full_Name), 'y', Image%Number_Of_Lines)
        endif
        
        !Done for HSD files - WCS3
        if (index(Image%Level1b_Name, '.nc') < 1) then
#ifdef LIBHIM
            Image%Nc_Format_Flag = .false.
            !For HSD, we will always do static nav. Because of that, we don't need 
            ! Image%Number_Of_Elements and Image%Number_Of_Lines. These are determined 
            ! in DETERMINE_BOUNDS_STATIC_NAV for a given segment of data 
#else
            call MESG( "LibHimawari not installed. Cannot process HSD. Stopping", level = verb_lev % ERROR , color = 4 )
            stop
#endif

        endif
        !END - WCS3
                
        Image%Mixed_Resolution_Flag = .true.
        Image%Static_Nav_Flag = .true.
        
        exit test_loop
      endif

      !--- HIMAWARI-8 HCAST AHI netCDF Test
      if (index(Image%Level1b_Name, 'HC_H08') > 0) then
        Sensor%Sensor_Name = 'AHI'
        Sensor%Platform_Name = 'HIM8'
        Sensor%Spatial_Resolution_Meters = 4000
        Sensor%WMO_Id = 173
        Sensor%Instr_Const_File = 'him8_instr.dat'
        Sensor%Geo_Sub_Satellite_Longitude = 140.0   ! 140.66
        Sensor%Geo_Sub_Satellite_Latitude = 0.0
        Static_Nav_File = 'him8_hcast_static_nav_fulldisk.nc'
        Image%Area_Format_Flag = .false.
        Image%DB_Flag = .true.
        Image%Nc_Format_Flag = .true.
        call READ_NETCDF_DIMENSION(trim(Image%Level1b_Full_Name), 'x', Image%Number_Of_Elements)
        call READ_NETCDF_DIMENSION(trim(Image%Level1b_Full_Name), 'y', Image%Number_Of_Lines)
        Image%Mixed_Resolution_Flag = .true.
        Image%Static_Nav_Flag = .true.
        if (Image%Number_Of_Elements ==  2750) then
            Image%Mixed_Resolution_Flag = .false.
            Image%Static_Nav_Flag = .false.
        endif
        exit test_loop
      endif

      !--- HIMAWARI-8 HCAST AHI Raw data Test
      ! WCS 12/21/19
      if ((index(Image%Level1b_Name, 'IMG_DK') > 0) .AND. &
           (index(Image%Level1b_Name, '.nc') < 1) ) then !ensure no netCDF files
#ifdef LIBHIM
        Sensor%Sensor_Name = 'AHI'
        Sensor%Platform_Name = 'HIM8'
        Sensor%Spatial_Resolution_Meters = 4000
        Sensor%WMO_Id = 173
        Sensor%Instr_Const_File = 'him8_instr.dat'
        Sensor%Geo_Sub_Satellite_Longitude = 140.0   ! 140.66
        Sensor%Geo_Sub_Satellite_Latitude = 0.0
        Static_Nav_File = 'him8_hcast_static_nav_fulldisk.nc'
        Image%Area_Format_Flag = .false.
        Image%DB_Flag = .true.
        Image%Mixed_Resolution_Flag = .true.
        Image%Static_Nav_Flag = .true.
        Image%Nc_Format_Flag = .false.
        call MESG( "HCAST DATA", level = verb_lev % ERROR , color = 4 )
            !For HSD, we will always do static nav. Because of that, we don't need 
            ! Image%Number_Of_Elements and Image%Number_Of_Lines. These are determined 
            ! in DETERMINE_BOUNDS_STATIC_NAV for a given segment of data 
#else
        call MESG( "LibHimawari not installed. Cannot process HSD. Stopping", level = verb_lev % ERROR , color = 4 )
            stop


#endif

        exit test_loop
      endif


      !--- GOES-16 ABI Static Nav Test
      if (index(Image%Level1b_Name, 'OR_ABI') > 0 .AND. INDEX(trim(Image%Level1b_Name), 'G16') > 0 &
          .AND. .NOT. INDEX(trim(Image%Level1b_Name) , '_HIGH_RES_') > 0  ) then
        Sensor%Sensor_Name = 'GOES-RU-IMAGER'
        Sensor%Platform_Name = 'GOES-16'
        Sensor%Spatial_Resolution_Meters = 2000
        Sensor%WMO_Id = 270
        Sensor%Instr_Const_File = 'goes_16_instr.dat'
        Sensor%Geo_Sub_Satellite_Longitude = -75.2   ! -75.2- not projection sub-point
        Sensor%Geo_Sub_Satellite_Latitude = 0.0
        call COMPUTE_GOES_RU_STATIC_NAV_FILE_NAME()
        Image%Static_Nav_Flag = .true.
        Image%Mixed_Resolution_Flag = .true.
        Image%Area_Format_Flag = .false.
        Image%Nc_Format_Flag = .true. !added in to make sure that static nav will run due to HSD test - WCS3
        ABI_Use_104um_Flag = .false.
        exit test_loop
      endif
      
      
      
            !--- GOES-16 ABI Static Nav Test FD HIGH RES processing
      if (index(Image%Level1b_Name, 'OR_ABI-L1b-RadC') > 0 .AND. INDEX(trim(Image%Level1b_Name), 'G16') > 0 &
          .AND. INDEX(trim(Image%Level1b_Name) , '_HIGH_RES_') > 0 ) then
        Sensor%Sensor_Name = 'GOES-RU-IMAGER'
        Sensor%Platform_Name = 'GOES-16'
        Sensor%Spatial_Resolution_Meters = 500
        Sensor%WMO_Id = 270
        Sensor%Instr_Const_File = 'goes_16_instr.dat'
        Sensor%Geo_Sub_Satellite_Longitude = -75.2   ! -75.2- not projection sub-point
        Sensor%Geo_Sub_Satellite_Latitude = 0.0
        Static_Nav_File = "goes_east_abi_conus_static_nav_05km.nc"
        Image%Static_Nav_Flag = .true.
        Image%Mixed_Resolution_Flag = .true.
        Image%Area_Format_Flag = .false.
        Image%Nc_Format_Flag = .true. 
        exit test_loop
      endif
      
      
                  !--- GOES-16 ABI Static Nav Test FD HIGH RES processing
      if (index(Image%Level1b_Name, 'OR_ABI-L1b-RadF') > 0 .AND. INDEX(trim(Image%Level1b_Name), 'G16') > 0 &
          .AND. INDEX(trim(Image%Level1b_Name) , '_HIGH_RES_') > 0 ) then
        Sensor%Sensor_Name = 'GOES-RU-IMAGER'
        Sensor%Platform_Name = 'GOES-16'
        Sensor%Spatial_Resolution_Meters = 500
        Sensor%WMO_Id = 270
        Sensor%Instr_Const_File = 'goes_16_instr.dat'
        Sensor%Geo_Sub_Satellite_Longitude = -75.2   ! -75.2- not projection sub-point
        Sensor%Geo_Sub_Satellite_Latitude = 0.0
        Static_Nav_File = "goes_east_abi_fulldisk_static_nav_05km.nc"
        Image%Static_Nav_Flag = .true.
        Image%Mixed_Resolution_Flag = .true.
        Image%Area_Format_Flag = .false.
        Image%Nc_Format_Flag = .true. 
        Image%Number_Of_Lines = 21696
        Image%Number_Of_Elements = 21696
        exit test_loop
      endif
      

      !--- GOES-17 ABI Static Nav Test
      if (index(Image%Level1b_Name, 'OR_ABI') > 0 .and. INDEX(trim(Image%Level1b_Name), 'G17') > 0 &
          .AND. .NOT. INDEX(trim(Image%Level1b_Name) , '_HIGH_RES_') > 0  ) then 
        Sensor%Sensor_Name = 'GOES-RU-IMAGER'
        Sensor%Platform_Name = 'GOES-17'
        Sensor%Spatial_Resolution_Meters = 2000
        Sensor%WMO_Id = 271
        Sensor%Instr_Const_File = 'goes_17_instr_lhp.dat'
        Sensor%Geo_Sub_Satellite_Longitude = -137.2   ! -137.2 - not projection sub-point
        Sensor%Geo_Sub_Satellite_Latitude = 0.0
        call COMPUTE_GOES_RU_STATIC_NAV_FILE_NAME()
        Image%Static_Nav_Flag = .true.
        Image%Mixed_Resolution_Flag = .true.
        Image%Area_Format_Flag = .false.
        Image%Nc_Format_Flag = .true. !added in to make sure that static nav will run due to HSD test - WCS3
        ABI_Use_104um_Flag = .false.
        exit test_loop
      endif
      
         !--- GOES-17 ABI Static Nav Test
      if (index(Image%Level1b_Name, 'OR_ABI') > 0 .and. INDEX(trim(Image%Level1b_Name), 'G17') > 0  &
          .AND. INDEX(trim(Image%Level1b_Name) , '_HIGH_RES_') > 0 ) then
        Sensor%Sensor_Name = 'GOES-RU-IMAGER'
        Sensor%Platform_Name = 'GOES-17'
        Sensor%Spatial_Resolution_Meters = 500
        Sensor%WMO_Id = 271
        Sensor%Instr_Const_File = 'goes_17_instr_lhp.dat'
        Sensor%Geo_Sub_Satellite_Longitude = -137.2   ! -137.2 - not projection sub-point
        Sensor%Geo_Sub_Satellite_Latitude = 0.0
        Static_Nav_File = "goes_west_abi_conus_static_nav_05km.nc"
        Image%Static_Nav_Flag = .true.
        Image%Mixed_Resolution_Flag = .true.
        Image%Area_Format_Flag = .false.
        Image%Nc_Format_Flag = .true. !added in to make sure that static nav will run due to HSD test - WCS3
        ABI_Use_104um_Flag = .false.
        exit test_loop
      endif


      !--- FY-4a AGRI Test
      if (index(Image%Level1b_Name, 'FY4A') > 0) then
        Sensor%Sensor_Name = 'AGRI'
        Sensor%Platform_Name = 'FY4A'
        Sensor%Spatial_Resolution_Meters = 4000
        Sensor%WMO_Id = 530
        Sensor%Instr_Const_File = 'fy4a_instr.dat'
        Sensor%Geo_Sub_Satellite_Longitude = 105.0
        Sensor%Geo_Sub_Satellite_Latitude = 0.0
        Static_Nav_File = 'fy4_static_nav_fulldisk.nc'
        call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(Image%Level1b_Full_Name),'RegWidth', Image%Number_Of_Elements) 
        call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(Image%Level1b_Full_Name),'RegLength', Image%Number_Of_Lines) 
        Image%Area_Format_Flag = .false.
        Image%Mixed_Resolution_Flag = .false.
        Image%Static_Nav_Flag = .true.
        exit test_loop
      endif
      
      !--- MODIS Test
      if (index(Image%Level1b_Name, 'MYD021KM') > 0) then
        Sensor%Sensor_Name = 'MODIS'
        Sensor%Platform_Name = 'AQUA'
        Sensor%Spatial_Resolution_Meters = 1000
        Sensor%WMO_Id = 784
        Sensor%Instr_Const_File = 'modis_aqua_instr.dat'
        exit test_loop
      endif

      if (index(Image%Level1b_Name, 'MYD02SSH') > 0) then
        Sensor%Sensor_Name = 'MODIS'
        Sensor%Platform_Name = 'AQUA'
        Sensor%Spatial_Resolution_Meters = 5000
        Sensor%Instr_Const_File = 'modis_aqua_instr.dat'
        Sensor%WMO_Id = 784
        exit test_loop
      endif

      if (index(Image%Level1b_Name, 'MYDATML') > 0) then
        Sensor%Sensor_Name = 'MODIS'
        Sensor%Platform_Name = 'AQUA'
        Sensor%Spatial_Resolution_Meters = 5000
        Sensor%Instr_Const_File = 'modis_aqua_instr.dat'
        Sensor%WMO_Id = 784
        exit test_loop
      endif

      if (index(Image%Level1b_Name, 'MOD021KM') > 0) then
        Sensor%Sensor_Name = 'MODIS'
        Sensor%Platform_Name = 'TERRA'
        Sensor%Spatial_Resolution_Meters = 1000
        Sensor%WMO_Id = 783
        Sensor%Instr_Const_File = 'modis_terra_instr.dat'
        exit test_loop
      endif

      if (index(Image%Level1b_Name, 'MOD02SSH') > 0) then
        Sensor%Sensor_Name = 'MODIS'
        Sensor%Platform_Name = 'TERRA'
        Sensor%Spatial_Resolution_Meters = 5000
        Sensor%WMO_Id = 783
        Sensor%Instr_Const_File = 'modis_terra_instr.dat'
        exit test_loop
      endif

      if (index(Image%Level1b_Name, 'MAC02') > 0) then
        Sensor%Sensor_Name = 'MODIS-MAC'
        Sensor%Platform_Name = 'AQUA'
        Sensor%Spatial_Resolution_Meters = 1000
        Sensor%WMO_Id = 784
        Sensor%Instr_Const_File = 'modis_aqua_instr.dat'
        exit test_loop
      endif

      if (index(Image%Level1b_Name, 'a1.') > 0) then
        Sensor%Sensor_Name = 'MODIS-CSPP'
        Sensor%Platform_Name = 'AQUA'
        Sensor%Spatial_Resolution_Meters = 1000
        Sensor%WMO_Id = 784
        Sensor%Instr_Const_File = 'modis_aqua_instr.dat'
        exit test_loop
      endif

      if (index(Image%Level1b_Name, 't1.') > 0) then
        Sensor%Sensor_Name = 'MODIS-CSPP'
        Sensor%Platform_Name = 'TERRA'
        Sensor%Spatial_Resolution_Meters = 1000
        Sensor%WMO_Id = 783
        Sensor%Instr_Const_File = 'modis_terra_instr.dat'
        exit test_loop
      endif

      !---  Check Geostationary (assumed to be areafiles)
      call GET_GOES_HEADERS(Image%Level1b_Full_Name, AREAstr, NAVstr)

      Image%Area_Format_Flag = .false.

      if (AREAstr%Version_Num == 4) then                          !begin valid Areafile test

         Image%Area_Format_Flag = .true.

         !--- set spatial resolution  !(AKH??? - Is this valid for all sensors)
         Sensor%Spatial_Resolution_Meters = int(AREAstr%Elem_Res,kind=int4)

         !--- based on McIdas Id, set sensor struc parameters
         select case(AREAstr%Sat_Id_Num)
 
            !test for SEVIRI
            case (51:53,354)
               Sensor%Sensor_Name = 'SEVIRI'
               Sensor%Spatial_Resolution_Meters = 3000
              
               if (AREAstr%Sat_Id_Num == 51) then
                  Sensor%Platform_Name = 'Meteosat-8'
                  Sensor%WMO_Id = 55
                  Sensor%Instr_Const_File = 'met8_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 52)  then
                  Sensor%Platform_Name = 'Meteosat-9'
                  Sensor%WMO_Id = 56
                  Sensor%Instr_Const_File = 'met9_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 53) then
                  Sensor%Platform_Name = 'Meteosat-10'
                  Sensor%WMO_Id = 57
                  Sensor%Instr_Const_File = 'met10_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 354) then
                  Sensor%Platform_Name = 'Meteosat-11'
                  Sensor%WMO_Id = 70
                  Sensor%Instr_Const_File = 'met11_instr.dat'
                  exit test_loop
               endif

            case (84,85)
               Sensor%Sensor_Name = 'MTSAT-IMAGER'
               Sensor%Spatial_Resolution_Meters = 4000
               if (AREAstr%Sat_Id_Num == 84) then
                  Sensor%Platform_Name = 'MTSAT-1R'
                  Sensor%WMO_Id = 171
                  Sensor%Instr_Const_File = 'mtsat1r_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 85) then
                  Sensor%Platform_Name = 'MTSAT-2'
                  Sensor%WMO_Id = 172
                  Sensor%Instr_Const_File = 'mtsat2_instr.dat'
                  exit test_loop
               endif

            case (36, 37)
               Sensor%Sensor_Name = 'FY2-IMAGER'
               Sensor%Spatial_Resolution_Meters = 4000
               if (AREAstr%Sat_Id_Num == 36) then
                  Sensor%Platform_Name = 'FY-2D'
                  Sensor%WMO_Id = 514
                  Sensor%Instr_Const_File = 'fy2d_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 37) then
                  Sensor%Platform_Name = 'FY-2E'
                  Sensor%WMO_Id = 515
                  Sensor%Instr_Const_File = 'fy2e_instr.dat'
                  exit test_loop
               endif

            case (250)
               Sensor%Sensor_Name = 'COMS-IMAGER'
               Sensor%Spatial_Resolution_Meters = 4000
               Sensor%Platform_Name = 'COMS-1'
               Sensor%WMO_Id = 810
               Sensor%Instr_Const_File = 'coms1_instr.dat'
               exit test_loop

            case (70,72,74,76,78,180,182,184,186)

               if (AREAstr%Sat_Id_Num == 70) then
                  Sensor%Sensor_Name = 'GOES-IL-IMAGER'
                  Sensor%Platform_Name = 'GOES-8'
                  Sensor%Spatial_Resolution_Meters = 4000
                  Sensor%WMO_Id = 252
                  Sensor%Instr_Const_File = 'goes_08_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 72) then
                  Sensor%Sensor_Name = 'GOES-IL-IMAGER'
                  Sensor%Platform_Name = 'GOES-9'
                  Sensor%Spatial_Resolution_Meters = 4000
                  Sensor%WMO_Id = 253
                  Sensor%Instr_Const_File = 'goes_09_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 74) then
                  Sensor%Sensor_Name = 'GOES-IL-IMAGER'
                  Sensor%Platform_Name = 'GOES-10'
                  Sensor%Spatial_Resolution_Meters = 4000
                  Sensor%WMO_Id = 254
                  Sensor%Instr_Const_File = 'goes_10_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 76) then
                  Sensor%Sensor_Name = 'GOES-IL-IMAGER'
                  Sensor%Platform_Name = 'GOES-11'
                  Sensor%Spatial_Resolution_Meters = 4000
                  Sensor%WMO_Id = 255
                  Sensor%Instr_Const_File = 'goes_11_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 78) then
                  Sensor%Sensor_Name = 'GOES-MP-IMAGER'
                  Sensor%Platform_Name = 'GOES-12'
                  Sensor%Spatial_Resolution_Meters = 4000
                  Sensor%WMO_Id = 256
                  Sensor%Instr_Const_File = 'goes_12_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 180) then
                  Sensor%Sensor_Name = 'GOES-MP-IMAGER'
                  Sensor%Platform_Name = 'GOES-13'
                  Sensor%Spatial_Resolution_Meters = 4000
                  Sensor%WMO_Id = 257
                  Sensor%Instr_Const_File = 'goes_13_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 182) then
                  Sensor%Sensor_Name = 'GOES-MP-IMAGER'
                  Sensor%Platform_Name = 'GOES-14'
                  Sensor%Spatial_Resolution_Meters = 4000
                  Sensor%WMO_Id = 258
                  Sensor%Instr_Const_File = 'goes_14_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 184) then
                  Sensor%Sensor_Name = 'GOES-MP-IMAGER'
                  Sensor%Platform_Name = 'GOES-15'
                  Sensor%Spatial_Resolution_Meters = 4000
                  Sensor%WMO_Id = 259
                  Sensor%Instr_Const_File = 'goes_15_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 186) then
                  Sensor%Sensor_Name = 'GOES-RU-IMAGER'
                  Sensor%Spatial_Resolution_Meters = 2000
                  Sensor%Platform_Name = 'GOES-16'
                  Sensor%WMO_Id = 270
                  Sensor%Instr_Const_File = 'goes_16_instr.dat'
                  exit test_loop
               endif

            case (71,73,75,77,79,181,183,185)
               if (AREAstr%Sat_Id_Num == 71) then
                  Sensor%Sensor_Name = 'GOES-IL-SOUNDER'
                  Sensor%Spatial_Resolution_Meters = 10000
                  Sensor%Platform_Name = 'GOES-8'
                  Sensor%WMO_Id = 252
                  Sensor%Instr_Const_File = 'goes_08_sndr_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 73) then
                  Sensor%Sensor_Name = 'GOES-IL-SOUNDER'
                  Sensor%Spatial_Resolution_Meters = 10000
                  Sensor%Platform_Name = 'GOES-9'
                  Sensor%WMO_Id = 253
                  Sensor%Instr_Const_File = 'goes_09_sndr_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 75) then
                  Sensor%Sensor_Name = 'GOES-IL-SOUNDER'
                  Sensor%Spatial_Resolution_Meters = 10000
                  Sensor%Platform_Name = 'GOES-10'
                  Sensor%WMO_Id = 254
                  Sensor%Instr_Const_File = 'goes_10_sndr_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 77) then
                  Sensor%Sensor_Name = 'GOES-IL-SOUNDER'
                  Sensor%Spatial_Resolution_Meters = 10000
                  Sensor%Platform_Name = 'GOES-11'
                  Sensor%WMO_Id = 255
                  Sensor%Instr_Const_File = 'goes_11_sndr_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 79) then
                  Sensor%Sensor_Name = 'GOES-IL-SOUNDER'
                  Sensor%Spatial_Resolution_Meters = 10000
                  Sensor%Platform_Name = 'GOES-12'
                  Sensor%WMO_Id = 256
                  Sensor%Instr_Const_File = 'goes_12_sndr_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 181) then
                  Sensor%Sensor_Name = 'GOES-IL-SOUNDER'
                  Sensor%Spatial_Resolution_Meters = 10000
                  Sensor%Platform_Name = 'GOES-13'
                  Sensor%WMO_Id = 257
                  Sensor%Instr_Const_File = 'goes_13_sndr_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 183) then
                  Sensor%Sensor_Name = 'GOES-IL-SOUNDER'
                  Sensor%Spatial_Resolution_Meters = 10000
                  Sensor%Platform_Name = 'GOES-14'
                  Sensor%WMO_Id = 258
                  Sensor%Instr_Const_File = 'goes_14_sndr_instr.dat'
                  exit test_loop
               endif
               if (AREAstr%Sat_Id_Num == 185) then
                  Sensor%Sensor_Name = 'GOES-IL-SOUNDER'
                  Sensor%Spatial_Resolution_Meters = 100000
                  Sensor%Platform_Name = 'GOES-15'
                  Sensor%WMO_Id = 259
                  Sensor%Instr_Const_File = 'goes_15_sndr_instr.dat'
                  exit test_loop
               endif

            end select
      endif

      !---  VIIRS SNPP
      if (index(Image%Level1b_Name, 'GMTCO_npp') > 0) then 
         Sensor%Sensor_Name = 'VIIRS'
         Sensor%Spatial_Resolution_Meters = 750
         Sensor%Platform_Name = 'SNPP'
         Sensor%WMO_Id = 224
         Sensor%Instr_Const_File = 'viirs_npp_instr.dat'
         exit test_loop
      endif

      !--- NPP/JPSS IFF 
      if (index(Image%Level1b_Name, 'IFFSDR_npp') > 0 .or. &
          index(Image%Level1b_Name, 'IFFSVM_npp') > 0 .or. &
          index(Image%Level1b_Name, 'IFF_npp') > 0) then
         Sensor%Sensor_Name = 'VIIRS-IFF'
         Sensor%Spatial_Resolution_Meters = 750
         Sensor%Platform_Name = 'SNPP'
         Sensor%WMO_Id = 224
         Sensor%Instr_Const_File = 'iff_viirs_npp_instr.dat'
         exit test_loop
      endif

      !---  VIIRS-NASA SNPP
      if (index(Image%Level1b_Name, 'VGEOM_snpp') > 0 .or. &
          index(Image%Level1b_Name, 'VNP03MOD') > 0) then
         Sensor%Sensor_Name = 'VIIRS-NASA'
         Sensor%Spatial_Resolution_Meters = 750
         Sensor%Platform_Name = 'SNPP'
         Sensor%WMO_Id = 224
         ! - check if it is FUSION
         Sensor%Fusion_Flag = .false.
         call CHECK_IF_FUSION(Sensor%Fusion_Flag)
         if (Sensor%Fusion_Flag) then
           Sensor%Instr_Const_File = 'fusion_npp_instr.dat'
           exit test_loop
         endif
         Sensor%Instr_Const_File = 'viirs_npp_instr.dat'
         exit test_loop
      endif
      
      
            !---  VIIRS-NASA SNPP HIGH RES
      if (index(Image%Level1b_Name, 'VNP02MOD') > 0 .or. &
          index(Image%Level1b_Name, 'highres') > 0) then
         Sensor%Sensor_Name = 'VIIRS-NASA-HRES'
         Sensor%Spatial_Resolution_Meters = 375
         Sensor%Platform_Name = 'SNPP'
         Sensor%WMO_Id = 224
         ! - check if it is FUSION
         Sensor%Fusion_Flag = .false.
        
         Sensor%Instr_Const_File = 'viirs_npp_instr.dat'
         exit test_loop
      endif
           
      !---  VIIRS NOAA-20 (JPSS-1)
      if (index(Image%Level1b_Name, 'GMTCO_j01') > 0) then
         Sensor%Sensor_Name = 'VIIRS'
         Sensor%Spatial_Resolution_Meters = 750
         Sensor%Platform_Name = 'NOAA-20'
         Sensor%WMO_Id = 225
         Sensor%Instr_Const_File = 'viirs_20_instr.dat'
         exit test_loop
      endif

      !---  VIIRS-NASA NOAA-20 (JPSS-1)
      if (index(Image%Level1b_Name, 'VJ103MOD') > 0) then
         Sensor%Sensor_Name = 'VIIRS-NASA'
         Sensor%Spatial_Resolution_Meters = 750
         Sensor%Platform_Name = 'NOAA-20'
         Sensor%WMO_Id = 225
         ! - check if it is FUSION
         Sensor%Fusion_Flag = .false.
         call CHECK_IF_FUSION(Sensor%Fusion_Flag)
         if (Sensor%Fusion_Flag) then
           Sensor%Instr_Const_File = 'fusion_20_instr.dat'
           exit test_loop
         endif
         Sensor%Instr_Const_File = 'viirs_20_instr.dat'
         exit test_loop
      endif

      !---  VGAC NOAA-20
      if (index(Image%Level1b_Name, 'VGAC_VJ1') > 0) then 
         Sensor%Sensor_Name = 'VGAC'
         Sensor%Spatial_Resolution_Meters = 3900
         Sensor%Platform_Name = 'NOAA-20'
         Sensor%WMO_Id = 225
         Sensor%Instr_Const_File = 'viirs_20_instr.dat'
         exit test_loop
      endif

      !---  EPS-SG / MetImage
      if (index(Image%Level1b_Name, 'SGA1-VII') > 0) then
         Sensor%Sensor_Name = 'METIMAGE'
         Sensor%Spatial_Resolution_Meters = 500
         Sensor%Platform_Name = 'EPS-SG'
         Sensor%WMO_Id = 384 ! TODO This needs to be change to real #, Fake for now
         !Sensor%Instr_Const_File = 'metimage_1_instr.dat'
         Sensor%Instr_Const_File = 'metimage_1_instr_terra.dat'
         exit test_loop
      endif

      !---  MERSI FY3D
      if (index(Image%Level1b_Name, 'FY3D') > 0) then
         Sensor%Sensor_Name = 'MERSI-2'
         Sensor%Spatial_Resolution_Meters = 1000
         Sensor%Platform_Name = 'FY3D'
         Sensor%WMO_Id = 523
         Sensor%Instr_Const_File = 'fy3d_instr.dat'
         exit test_loop
      endif

      !--- AQUA IFF
      if (index(Image%Level1b_Name, 'IFFSDR_aqua') > 0) then
         Sensor%Sensor_Name = 'AQUA-IFF'
         Sensor%Spatial_Resolution_Meters = 1000
         Sensor%Platform_Name = 'AQUA'
         Sensor%WMO_Id = 784
         Sensor%Instr_Const_File = 'modis_aqua_instr.dat'
         exit test_loop
      endif

      !--- AVHRR IFF
      !--- MJH: is it ok to use avhrr algo files here? or do we need new
      !         ones for AVHRR/HIRS?
      if (index(Image%Level1b_Name, 'IFF_noaa') > 0 .or. &
          index(Image%Level1b_Name, 'IFF_metop') > 0) then
            Sensor%Spatial_Resolution_Meters = 4000
            if (index(Image%Level1b_Name, 'IFF_noaa06') == 1) then
               AVHRR_KLM_Flag = sym%NO
               Sc_Id_AVHRR = 2_int1
               Sensor%Platform_Name = 'NOAA-6'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 706
               Sensor%Instr_Const_File = "iff_avhrr_6_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa07') == 1) then
               AVHRR_KLM_Flag = sym%NO
               Sc_Id_AVHRR = 4_int1
               Sensor%Platform_Name = 'NOAA-7'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 707
               Sensor%Instr_Const_File = "iff_avhrr_7_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa08') == 1) then
               AVHRR_KLM_Flag = sym%NO
               Sc_Id_AVHRR = 6_int1
               Sensor%Platform_Name = 'NOAA-8'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 200
               Sensor%Instr_Const_File = "iff_avhrr_8_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa09') == 1) then
               AVHRR_KLM_Flag = sym%NO
               Sc_Id_AVHRR = 7_int1
               Sensor%Platform_Name = 'NOAA-9'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 201
               Sensor%Instr_Const_File = "iff_avhrr_9_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa10') == 1) then
               AVHRR_KLM_Flag = sym%NO
               Sc_Id_AVHRR = 8_int1
               Sensor%Platform_Name = 'NOAA-10'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 202
               Sensor%Instr_Const_File = "iff_avhrr_10_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa11') == 1) then
               AVHRR_KLM_Flag = sym%NO
               Sc_Id_AVHRR = 1_int1
               Sensor%Platform_Name = 'NOAA-11'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 203
               Sensor%Instr_Const_File = "iff_avhrr_11_instr.dat"
               exit test_loop
            endif
            if (index(Image%Level1b_Name, 'IFF_noaa12') == 1) then
               AVHRR_KLM_Flag = sym%NO
               Sc_Id_AVHRR = 5_int1
               Sensor%Platform_Name = 'NOAA-12'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 204
               Sensor%Instr_Const_File = "iff_avhrr_12_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa14') == 1) then
               AVHRR_KLM_Flag = sym%NO
               Sc_Id_AVHRR = 3_int1
               Sensor%Platform_Name = 'NOAA-14'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 205
               Sensor%Instr_Const_File = "iff_avhrr_14_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa15') == 1) then
               AVHRR_KLM_Flag = sym%YES
               Sc_Id_AVHRR = 4_int1
               Sensor%Platform_Name = 'NOAA-15'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 206
               Sensor%Instr_Const_File = "iff_avhrr_15_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa16') == 1) then
               AVHRR_KLM_Flag = sym%YES
               Sc_Id_AVHRR = 2_int1
               Sensor%Platform_Name = 'NOAA-16'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 207
               Sensor%Instr_Const_File = "iff_avhrr_16_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa17') == 1) then
               AVHRR_KLM_Flag = sym%YES
               Sc_Id_AVHRR = 6_int1
               Sensor%Platform_Name = 'NOAA-17'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 208
               Sensor%Instr_Const_File = "iff_avhrr_17_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa18') == 1) then
               AVHRR_KLM_Flag = sym%YES
               Sc_Id_AVHRR = 7_int1
               Sensor%Platform_Name = 'NOAA-18'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 209
               Sensor%Instr_Const_File = "iff_avhrr_18_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_noaa19') == 1) then
               AVHRR_KLM_Flag = sym%YES
               Sc_Id_AVHRR = 8_int1
               Sensor%Platform_Name = 'NOAA-19'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 223
               Sensor%Instr_Const_File = "iff_avhrr_19_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_metop02') == 1) then
               AVHRR_KLM_Flag = sym%YES
               Sc_Id_AVHRR = 12_int1
               Sensor%Platform_Name = 'METOP-A'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 4
               Sensor%Instr_Const_File = "iff_avhrr_2_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_metop01') == 1) then
               AVHRR_KLM_Flag = sym%YES
               Sc_Id_AVHRR = 11_int1
               Sensor%Platform_Name = 'METOP-B'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 3
               Sensor%Instr_Const_File = "iff_avhrr_1_instr.dat"
               exit test_loop
            end if
            if (index(Image%Level1b_Name, 'IFF_metop03') == 1) then
               AVHRR_KLM_Flag = sym%YES
               Sc_Id_AVHRR = 13_int1 ! Metop-C Sc_Id numbers are not known at this time
               Sensor%Platform_Name = 'METOP-C'
               Sensor%Sensor_Name = 'AVHRR-IFF'
               Sensor%WMO_Id = 5
               Sensor%Instr_Const_File = "iff_avhrr_3_instr.dat"
               exit test_loop
            end if

      endif
      
      !-------------------------------------------------------------------------------
      !---if sensor not detected, assume AVHRR
      !-------------------------------------------------------------------------------

      !--- from a preliminary header read, set some global AVHRR flags and parameters
      call DETERMINE_AVHRR_FILE_TYPE(trim(Image%Level1b_Full_Name), &
                                     AVHRR_GAC_FLAG,AVHRR_KLM_Flag,AVHRR_AAPP_Flag, &
                                     AVHRR_Ver_1b,AVHRR_Data_Type,Byte_Swap_1b,AVHRR_1_Flag)

      !--- knowing Sc_Id_AVHRR and the above flags, populate sensor structure for AVHRR
      call ASSIGN_AVHRR_SAT_ID_NUM_INTERNAL()

      if (AVHRR_GAC_Flag == sym%YES) then
         Sensor%Spatial_Resolution_Meters = 4000
      else
         Sensor%Spatial_Resolution_Meters = 1000
      endif
      
      if ( AVHRR_Fusion_Flag) then
         sensor % sensor_name = 'AVHRR-FUSION'
      end if 
     
      ifound = sym%YES   ! force exit need to develop logic for setting Ierror

      enddo test_loop

      
      !---------------------------------------------------------------------------------
      ! Set sub-satellite point for geostationary satellites that are Areafiles
      !---------------------------------------------------------------------------------
      call DETERMINE_GEO_SUB_SATELLITE_POSITION(trim(Image%Level1b_Full_Name),AREAstr,NAVstr)

      !---------------------------------------------------------------------------------
      ! append full path on constant files
      !---------------------------------------------------------------------------------
      Sensor%Instr_Const_File = trim(Ancil_Data_Dir)//"static/clavrx_constant_files/"//trim(Sensor%Instr_Const_File)

      !---------------------------------------------------------------------------------
      ! For MODIS, determine names of additional files for level-1b processing
      !---------------------------------------------------------------------------------

      !-- for 1 km MODIS, determine name of separate geolocation file
      if ((trim(Sensor%Sensor_Name) == 'MODIS' .or. trim(Sensor%Sensor_Name) == 'MODIS-CSPP') &
          .and.  Sensor%Spatial_Resolution_Meters == 1000) then
         call DETERMINE_MODIS_GEOLOCATION_FILE(Image%Level1b_Name,Image%Level1b_Path,Image%Auxiliary_Geolocation_File_Name)
         if (trim(Image%Auxiliary_Geolocation_File_Name) == "no_file") then
            Ierror = sym%YES
         endif
      endif


      !-- determine modis cloud mask name
      if ((trim(Sensor%Sensor_Name) == 'MODIS' .or. trim(Sensor%Sensor_Name) == 'MODIS-CSPP') &
          .and. Use_Aux_Flag /= sym%NO_AUX) then

         call DETERMINE_MODIS_CLOUD_MASK_FILE(Image%Level1b_Name,Image%Level1b_Path,Image%Auxiliary_Cloud_Mask_File_Name )
         if (trim(Image%Auxiliary_Cloud_Mask_File_Name) == "no_file" .and. &
                  Cloud_Mask_Bayesian_Flag == sym%NO) then
            Ierror = sym%YES
         endif

      endif


   end subroutine DETECT_SENSOR_FROM_FILE

   !=====================================================================
   !
   !=====================================================================
   subroutine DETERMINE_GEO_SUB_SATELLITE_POSITION(Level1b_Full_Name,AREAstr,NAVstr)

    character(len=*), intent(in):: Level1b_Full_Name
    type (AREA_STRUCT), intent(inout) :: AREAstr
    type (GVAR_NAV), intent(inout)    :: NAVstr
    REAL (KIND=REAL4)                 :: Lat_temp, Lon_temp
    INTEGER (KIND=INT4)               :: Year_temp

    Sensor%Geo_Sub_Satellite_Longitude = Missing_Value_Real4
    Sensor%Geo_Sub_Satellite_Latitude = Missing_Value_Real4

    ! Calculate year of the image from the McIDAS AREA file.
    Year_temp = 1900 + int(AREAstr%img_Date / 1000)

    if (AREAstr%Version_Num == 4) then                          !begin valid Areafile test

      select case(AREAstr%Sat_Id_Num)
 
            !test for SEVIRI
            case (51:53,354)
               ! Read the satellite sub longitude point from the AREA file.
               call READ_NAVIGATION_BLOCK_SEVIRI(trim(Level1b_Full_Name), AREAstr, NAVstr)
               Sensor%Geo_Sub_Satellite_Latitude = NAVstr%sublat 
               Sensor%Geo_Sub_Satellite_Longitude = NAVstr%sublon
               ! Override the above for the operational sub satellite longitudes.
               ! Longitude of actual Sub-Satellite Point for Met-8 when it was operational.  For Met-8 Indian
               ! Ocean service, the subpoint from the AREA file is used.
               if (AREAstr%Sat_Id_Num == 51 .AND. Year_temp < 2016 ) Sensor%Geo_Sub_Satellite_Longitude = -3.477996 
               if (AREAstr%Sat_Id_Num == 52 ) Sensor%Geo_Sub_Satellite_Longitude = -0.159799     ! Longitude of actual Sub-Satellite Point for Met-9
               if (AREAstr%Sat_Id_Num == 53 ) Sensor%Geo_Sub_Satellite_Longitude = 0.06          ! Longitude of actual Sub-Satellite Point for Met-10
               if (AREAstr%Sat_Id_Num == 354 ) Sensor%Geo_Sub_Satellite_Longitude = 0.26         ! Longitude of actual Sub-Satellite Point for Met-11
               AREAstr%Cal_Offset = AREAstr%reserved(3)
                    
            !test for MTSAT
            case (84,85)
               !This is needed to determine type of navigation
               !as Nav coefficents specific to MTSAT (and FY2)
               call READ_NAVIGATION_BLOCK_MTSAT_FY(trim(Level1b_Full_Name), AREAstr, NAVstr)
               Sensor%Geo_Sub_Satellite_Latitude = NAVstr%sublat
               Sensor%Geo_Sub_Satellite_Longitude = NAVstr%sublon
          
            !WCS3 - FY2-D AREA files have the subsat lat/lon flipped. Fix here
            !        Fix with McIDAS-X is fixed.
            case (36, 37)
               !This is needed to determine type of navigation
               !as Nav coefficents specific to FY2D/E. They are stored in
               ! the same manner as MTSAT, hence using the same routine
               call READ_NAVIGATION_BLOCK_MTSAT_FY(trim(Level1b_Full_Name), AREAstr, NAVstr)
               
               !Some data from BOM has subpoints flipped, so need to fix that
               IF (NAVstr%sublat .GT. 10.0) THEN
                    Lat_temp = NAVstr%sublon                   
                    Lon_temp = NAVstr%sublat
                    
                    NAVstr%sublon = Lon_temp
                    NAVstr%sublat = Lat_temp
               
               endif
               
               
               Sensor%Geo_Sub_Satellite_Latitude = NAVstr%sublat
               Sensor%Geo_Sub_Satellite_Longitude = NAVstr%sublon
 
            !test for COMS
            case (250)
               !This is needed to determine type of navigation
               !as Nav coefficents specific to COMS
               call READ_NAVIGATION_BLOCK_COMS(trim(Level1b_Full_Name), AREAstr,NAVstr)
               Sensor%Geo_Sub_Satellite_Latitude = NAVstr%sublat
               Sensor%Geo_Sub_Satellite_Longitude = NAVstr%sublon

            !test for GOES-16
            case (186)
               ! This is needed to determine type of navigation
               ! as Nav coefficents specific to GOES-16. These
               ! coefficients are based on 1 km data, not 2 km.
               ! Navigation transformations in abi_mod.f90 will
               ! account for this difference.

               call READ_NAVIGATION_BLOCK_ABI(trim(Level1b_Full_Name), AREAstr,NAVstr)

               !--- There will need to be hard coded values for AREA file
               !--- sub_lon calculations. Currently, the AREA file contains
               !--- -75.0.  It should be -75.2. 

               if ((NAVstr%sublon <= -75.) .AND. (NAVstr%sublon > -80.)) then
                 Sensor%Geo_Sub_Satellite_Latitude = NAVstr%sublat
                 !Sensor%Geo_Sub_Satellite_Longitude = NAVstr%sublon
                 Sensor%Geo_Sub_Satellite_Longitude = -75.21
               else
                 Sensor%Geo_Sub_Satellite_Latitude = NAVstr%sublat
                 Sensor%Geo_Sub_Satellite_Longitude = NAVstr%sublon
               endif
   
            !test for GOES Imagers or Sounders
            case (70:79,180:185)

               call LMODEL(Goes_Input_Time,  &
                           Goes_Epoch_Time, &
                           NAVstr, &
                           Sensor%Geo_Sub_Satellite_Latitude, &
                           Sensor%Geo_Sub_Satellite_Longitude)

                      Sensor%Geo_Sub_Satellite_Longitude = Sensor%Geo_Sub_Satellite_Longitude / Dtor
                      Sensor%Geo_Sub_Satellite_Latitude = Sensor%Geo_Sub_Satellite_Latitude  / Dtor

            case default
                print *, "Could not determine geostationary sub-satellite point"
                stop

      end select

    endif

   end subroutine DETERMINE_GEO_SUB_SATELLITE_POSITION
   
   !---------------------------------------------------------------------------------------------
   ! Determine the number of elements (Image%Number_Of_Elements) and Number of Scans (Image%Number_Of_Lines)
   ! expected.  Also, 
   !
   !    the output will be written in global Image structure
   !---------------------------------------------------------------------------------------------
   subroutine SET_FILE_DIMENSIONS(Level1b_Full_Name,AREAstr,Nrec_Avhrr_Header, Ierror)

      use CX_READ_AHI_MOD, only : AHI_SEGMENT_INFORMATION_REGION, AHI_CONFIG_TYPE
                                                               
      character(len=*), intent(in) :: Level1b_Full_Name
      type (AREA_STRUCT), intent(in) :: AREAstr ! AVHRR only
      integer(kind=int4), intent(out) :: Nrec_Avhrr_Header ! AVHRR only
      integer(kind=int4), intent(out) :: Ierror

      integer(kind=int4) :: Nword_Clavr
      integer(kind=int4) :: Nword_Clavr_Start
      integer(kind=int4) :: Ierror_Nscans
      CHARACTER(len=1020) :: Dir_File
      
      type ( Ahi_Config_Type ) :: Ahi_Config
      integer :: Offset(2), count(2)

      Ierror = sym%NO

      if (index(Sensor%Sensor_Name,'MODIS') > 0) then
         call READ_MODIS_SIZE_ATTR(trim(Level1b_Full_Name),Image%Number_Of_Elements,Image%Number_Of_Lines)
      endif
   
      if (trim(Sensor%Sensor_Name) == 'MODIS-MAC') then
         Image%Number_Of_Elements =  11
         Image%Number_Of_Lines = 2030
      endif
      
      if ( trim(Sensor%Sensor_Name) == 'AHI') then
     
         if (Image%Static_Nav_Flag) then
            call SETUP_READ_LEVEL1B_FIXED_GRID_STATIC_NAV()
         else
            ahi_config % data_path = trim(Image%Level1b_Path)
            ahi_config % file_base = trim (Image%level1b_name)
            ahi_config % lon_range =[Nav%Lon_West_Limit,Nav%Lon_East_Limit]
            ahi_config % lat_range =[Nav%Lat_South_Limit,Nav%Lat_North_Limit]
            call ahi_segment_information_region ( ahi_config , offset, count )
        
            Image%Number_Of_Elements =  count(1)
            Image%Number_Of_Lines = count(2)
         endif


      end if
      
      if (trim(Sensor%Sensor_Name) == 'VIIRS') then
         Image%Number_Of_Elements = 3200
         Dir_File = trim(Image%Level1b_Path) // trim(Image%Level1b_Name)

         call GET_NUMBER_OF_SCANS_FROM_VIIRS_BRIDGE (trim(Dir_File),Image%Number_Of_Lines,Ierror_Nscans)

         ! If error reading, then go to next file
         if (Ierror_Nscans /= 0) then
            Ierror = sym%YES
            return      ! skips file
         endif

         ! Check if VIIRS Number of scans is regular (48) and calculate Number of y pixels
         if (Image%Number_Of_Lines .ge. 48) then
            Image%Number_Of_Lines = Image%Number_Of_Lines * 16      !16pix per Scan
         else if (Image%Number_Of_Lines == 47) then
            Image%Number_Of_Lines = (Image%Number_Of_Lines+1) * 16
         else
            Ierror = sym%YES
            return      !skips file
         end if

      end if
      
      if (trim(Sensor%Sensor_Name) == 'VIIRS-NASA-HRES') then
       
          Image%Number_Of_Elements = 6400
          Image%Number_Of_Lines    = 6464
      end if

      if (trim(Sensor%Sensor_Name) == 'VIIRS-NASA') then
       
          Image%Number_Of_Elements = 3200
          Dir_File = trim(Image%Level1b_Path) // trim(Image%Level1b_Name)

          call READ_NUMBER_OF_SCANS_VIIRS_NASA (trim(Dir_File),Image%Number_Of_Lines,Ierror_Nscans)

          ! If error reading, then go to next file
          if (Ierror_Nscans /= 0) then
            Ierror = sym%YES
            return      ! skips file
          endif
        
      endif

      if (trim(Sensor%Sensor_Name) == 'VGAC') then
         Dir_File = trim(Image%Level1b_Path) // trim(Image%Level1b_Name)

         call READ_NUMBER_OF_SCANS_VGAC (trim(Dir_File),Image%Number_Of_Lines, &
                                         Image%Number_Of_ELements,Ierror_Nscans)
         ! If error reading, then go to next file
         if (Ierror_Nscans /= 0) then
            Ierror = sym%YES
            return      ! skips file
         endif
      endif

      if (trim(Sensor%Sensor_Name) == 'METIMAGE') then
         Dir_File = trim(Image%Level1b_Path) // trim(Image%Level1b_Name)

         call READ_NUMBER_OF_SCANS_EPS_SG (trim(Dir_File),Image%Number_Of_Lines, &
                                         Image%Number_Of_ELements,Ierror_Nscans)
         ! If error reading, then go to next file
         if (Ierror_Nscans /= 0) then
            Ierror = sym%YES
            return      ! skips file
         endif
      endif

      if (trim(Sensor%Sensor_Name) == 'MERSI-2') then
         Image%Number_Of_Elements = 2048
         Dir_File = trim(Image%Level1b_Path) // trim(Image%Level1b_Name)
         call READ_NUMBER_OF_SCANS_FY3D(trim(Dir_File),Image%Number_Of_Lines,Ierror_Nscans)

         ! If error reading, then go to next file
         if (Ierror_Nscans /= 0) then
            Ierror = sym%YES
            return      ! skips file
         endif
      endif
      
      !--- if an IFF, call routine to determine dimensions from latitude sds
      if (index(Sensor%Sensor_Name,'IFF') > 0) then
         call GET_IFF_DIMS_BRIDGE(trim(Image%Level1b_Path)//trim(Image%Level1b_Name),Image%Number_Of_Elements,Image%Number_Of_Lines)
      end if
     
      !--- AVHRR
      if (trim(Sensor%Sensor_Name) == 'AVHRR-1' .or. &
          trim(Sensor%Sensor_Name) == 'AVHRR-2' .or. &
          trim(Sensor%Sensor_Name) == 'AVHRR-3' .or. &
          trim(Sensor%Sensor_Name) == 'AVHRR-FUSION') then
         
         !-------------------------------------------------------
         ! Determine the type of level 1b file
         !-------------------------------------------------------
         call DETERMINE_AVHRR_FILE_TYPE(trim(Level1b_Full_Name),AVHRR_GAC_FLAG,AVHRR_KLM_Flag,AVHRR_AAPP_Flag, &
                                        AVHRR_Ver_1b,AVHRR_Data_Type,Byte_Swap_1b,AVHRR_1_Flag)
   
         !-------------------------------------------------------------------
         !-- based on file type (AVHRR_KLM_Flag and Gac), determine parameters needed
         !-- to read in header and data records for this orbit
         !------------------------------------------------------------------- 
         call DEFINE_1B_DATA(AVHRR_GAC_Flag,AVHRR_KLM_Flag,AVHRR_AAPP_Flag,Nrec_Avhrr_Header, &
                             Nword_Clavr_Start,Nword_Clavr)
 
         !-------------------------------------------------------------------
         !-- read in header
         !-------------------------------------------------------------------
         call READ_AVHRR_LEVEL1B_HEADER(trim(Level1b_Full_Name))

      end if

      !------------------------------------------------------------------------
      !  if GOES, SEVIRI and MTSAT, use elements of AREAstr to determine filesize
      !------------------------------------------------------------------------
      if (trim(Sensor%Sensor_Name) == 'GOES-IL-IMAGER' .or. &
          trim(Sensor%Sensor_Name) == 'GOES-MP-IMAGER') then
         Image%Number_Of_Elements =  int(AREAstr%Num_Elem / Goes_Xstride)
         Image%Number_Of_Lines = AREAstr%Num_Line
         L1b_Rec_Length = AREAstr%Num_Byte_Ln_Prefix +  &
                     (AREAstr%Num_Elem*AREAstr%Bytes_Per_Pixel)
      end if

      !--- Does this break AREA files?
      if (trim(Sensor%Sensor_Name) == 'GOES-RU-IMAGER') then
        if (Image%Area_Format_Flag) then
         Image%Number_Of_Elements =  int(AREAstr%Num_Elem)
         Image%Number_Of_Lines = AREAstr%Num_Line
        else
         if (Image%Static_Nav_Flag) then
         
          call SETUP_READ_LEVEL1B_FIXED_GRID_STATIC_NAV()
          
         else
          print *, "ERROR, unknown GOES-RU data"
          stop
         endif
        endif
      end if

      if (trim(Sensor%Sensor_Name) == 'GOES_IP_SOUNDER') then
         Image%Number_Of_Elements =  int(AREAstr%Num_Elem / Goes_Sndr_Xstride)
         Image%Number_Of_Lines = AREAstr%Num_Line
         L1b_Rec_Length = AREAstr%Num_Byte_Ln_Prefix +  &
                     (AREAstr%Num_Elem*AREAstr%Bytes_Per_Pixel)
      end if

      if (trim(Sensor%Sensor_Name) == 'SEVIRI') then
         Image%Number_Of_Elements =  int(AREAstr%Num_Elem)
         Image%Number_Of_Lines = AREAstr%Num_Line
      end if

      if (trim(Sensor%Sensor_Name) == 'MTSAT-IMAGER' .or. &
          trim(Sensor%Sensor_Name) == 'FY2-IMAGER' .or. &
          trim(Sensor%Sensor_Name) == 'COMS-IMAGER') then
         Image%Number_Of_Elements =  int(AREAstr%Num_Elem)
         Image%Number_Of_Lines = AREAstr%Num_Line
      end if


   end subroutine SET_FILE_DIMENSIONS

   !--------------------------------------------------------------------------------------------------
   ! Main Level1b Read Routine for all Sensors
   !--------------------------------------------------------------------------------------------------
   subroutine READ_LEVEL1B_DATA(Level1b_Full_Name &
        ,Segment_Number &
        ,Time_Since_Launch &
        ,AREAstr &
        ,NAVstr &
        ,Nrec_Avhrr_Header &
        ,Ierror_Level1b)

      character(len=*), intent(in):: Level1b_Full_Name
      integer, intent(in):: Segment_Number
      integer, intent(in):: Nrec_Avhrr_Header
      TYPE (AREA_STRUCT), intent(in) :: AREAstr
      TYPE (GVAR_NAV), intent(in)    :: NAVstr
      real, intent(in):: Time_Since_Launch
      integer, intent(out):: Ierror_Level1b
      TYPE(viirs_nasa_hres_config_type) :: nasa_hres_config
      integer :: i_line
      
      
      Ierror_Level1b = 0
      Cloud_Mask_Aux_Read_Flag = sym%NO

      if (index(Sensor%Sensor_Name,'MODIS') > 0) then
         call READ_MODIS(Segment_Number,Ierror_Level1b)
         if (Ierror_Level1b /= 0) return
      end if

      select case (trim(Sensor%Sensor_Name))

      case('GOES-IL-IMAGER','GOES-MP-IMAGER')
         call READ_GOES(Segment_Number,Image%Level1b_Name, &
                     Image%Start_Doy, Image%Start_Time, &
                     Time_Since_Launch, &
                     AREAstr,NAVstr)

         if (Sensor%Chan_On_Flag_Default(1)==sym%YES) then
            call READ_DARK_COMPOSITE_COUNTS(Segment_Number, Goes_Xstride, &
                     Dark_Composite_Name,AREAstr,Two_Byte_Temp) 
            call CALIBRATE_GOES_DARK_COMPOSITE(Two_Byte_Temp,Time_Since_Launch,Ref_Ch1_Dark_Composite)
         end if

      case('GOES-IP-SOUNDER')
         call READ_GOES_SNDR(Segment_Number,Image%Level1b_Name, &
                     Image%Start_Doy, Image%Start_Time, &
                     
                     AREAstr,NAVstr)

      case('GOES-RU-IMAGER')
         if (Image%Static_Nav_Flag) then
            if (Image%Area_Format_Flag) then
               call READ_ABI(Segment_Number,Image%Level1b_Name, &
                           Image%Start_Doy, Image%Start_Time, &
                           AREAstr,NAVstr)
            else
                call READ_LEVEL1B_FIXED_GRID_STATIC_NAV()
            end if

            !--- read auxillary cloud mask and cloud type
            if (Use_Aux_Flag /= sym%NO_AUX) then

               call DETERMINE_SAPF_NAME(Segment_Number)
               call READ_SAPF_DATA(Segment_Number)

            end if
         else 
            print*,'read abi is to installed stopping'
            stop
         end if

      case('SEVIRI')
       !--------  MSG/SEVIRI
         call READ_SEVIRI(Segment_Number,Image%Level1b_Name, &
                     Image%Start_Doy, Image%Start_Time, &
                     AREAstr)
         call READ_DARK_COMPOSITE_COUNTS(Segment_Number,Seviri_Xstride, &
                     Dark_Composite_Name,AREAstr,Two_Byte_Temp) 
         call CALIBRATE_SEVIRI_DARK_COMPOSITE(Two_Byte_Temp,Ref_Ch1_Dark_Composite)

      case('MTSAT-IMAGER')
         call READ_MTSAT(Segment_Number,Image%Level1b_Name, &
                     Image%Start_Doy, Image%Start_Time, &
                     Time_Since_Launch, &
                     AREAstr,NAVstr)
         call READ_DARK_COMPOSITE_COUNTS(Segment_Number,Mtsat_Xstride, &
                     Dark_Composite_Name,AREAstr,Two_Byte_Temp) 
         call CALIBRATE_MTSAT_DARK_COMPOSITE(Two_Byte_Temp,Ref_Ch1_Dark_Composite)

      case('FY2-IMAGER')
         call READ_FY(Segment_Number,Image%Level1b_Name, &
                     Image%Start_Doy, Image%Start_Time, &
                     AREAstr,NAVstr)

      case('AGRI') ! FY4
         call READ_FY4_LEVEL1B_DATA(Segment_Number, trim(Image%Level1b_Name), Ierror_Level1b)
        
      case('COMS-IMAGER')
         call READ_COMS(Segment_Number,Image%Level1b_Name, &
                     Image%Start_Doy, Image%Start_Time, &
                     AREAstr,NAVstr)

      case('AVHRR-1','AVHRR-2','AVHRR-3')
         call READ_AVHRR_LEVEL1B_DATA(trim(Level1b_Full_Name), &
              AVHRR_KLM_Flag,AVHRR_AAPP_Flag,Therm_Cal_1b,&
              Time_Since_Launch,Nrec_Avhrr_Header,Segment_Number)

       case ( 'AVHRR-FUSION')
            call READ_AVHRR_LEVEL1B_DATA(trim(Level1b_Full_Name), &
              AVHRR_KLM_Flag,AVHRR_AAPP_Flag,Therm_Cal_1b,&
              Time_Since_Launch,Nrec_Avhrr_Header,Segment_Number)
            call READ_HIRS_DATA(Segment_Number)
            call REPLACE_AVHRR_WITH_HIRS()
        
       case('VIIRS')

         call READ_VIIRS_DATA (Segment_Number, trim(Image%Level1b_Name), Ierror_Level1b)
      
         ! If error reading, then go to next file
         if (Ierror_Level1b /= 0) return

         !--- read auxillary cloud mask and cloud type
         if (Use_Aux_Flag /= sym%NO_AUX) then
           call DETERMINE_SAPF_NAME(Segment_Number)
           call READ_SAPF_DATA(Segment_Number)
         endif
        
        
      case('VIIRS-NASA-HRES')
          print*,'read routine has to be finished '
          print*, 'File: ',__FILE__,' Line: ',__LINE__
          print*,' +++++++++++++++++++++++++++++++++++'
          
          nasa_hres_config % channel_on_modis(1:45) = Sensor%Chan_On_Flag_Default(1:45)  == sym%YES
          
          nasa_hres_config % sensor = 'npp'
          nasa_hres_config % filename = trim(Image%Level1b_Name)
          nasa_hres_config % path = trim(Image%Level1b_Path)
          nasa_hres_config % ny_start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
          nasa_hres_config % ny_end = min(Image%Number_Of_Lines, nasa_hres_config % ny_start + Image%Number_of_Lines_Per_Segment - 1)
          call READ_VIIRS_NASA_HRES_DATA(nasa_hres_config)
          
          print*,'read ready..'
        
           Image%Number_Of_Lines_Read_This_Segment = nasa_hres_config % ny_end - nasa_hres_config % ny_start + 1
           do i_line = 1, Image%Number_Of_Lines_Per_Segment
              Image%Scan_Number(i_line) =nasa_hres_config % ny_start + i_line - 1
           end do
       case('VIIRS-NASA')

         call READ_VIIRS_NASA_DATA (Segment_Number, trim(Image%Level1b_Name), Ierror_Level1b)

         !--- If error reading, then go to next file
         if (Ierror_Level1b /= 0) return

         !--- read auxillary cloud mask
         if (Use_Aux_Flag /= sym%NO_AUX) then 
          call DETERMINE_MVCM_NAME(Segment_Number)
          call READ_MVCM_DATA(Segment_Number)
         endif

       case('MERSI-2')
         call READ_FY3D_DATA (Segment_Number, trim(Image%Level1b_Name), Ierror_Level1b)

         !--- If error reading, then go to next file
         if (Ierror_Level1b /= 0) return

       case('AHI')

         if (Image%Static_Nav_Flag) then
            IF (Image%Nc_Format_Flag) then 
                call READ_LEVEL1B_FIXED_GRID_STATIC_NAV()
            else
#ifdef LIBHIM
                call READ_HSD_FIXED_GRID_STATIC_NAV()   
#else
            call MESG( "LibHimawari not installed. Cannot process HSD. Stopping", level = verb_lev % ERROR , color = 4 )
            stop
#endif
                 
            endif
            
         else
            call READ_AHI_DATA (Segment_Number, trim(Image%Level1b_Name), Ierror_Level1b)
         endif

      end select

      !--- IFF data (all sensors same format)
      if (index(Sensor%Sensor_Name,'IFF') > 0) then
         call READ_IFF_DATA (Segment_Number, trim(Level1b_Full_Name), Ierror_Level1b)

         ! If error reading, then go to next file
         if (Ierror_Level1b /= 0) return

         !---- determine auxilliary cloud mask name and read it
         if (Use_Aux_Flag /= sym%NO_AUX) then
           if (Segment_Number == 1) print *,'Searching and reading MVCM'
          call DETERMINE_MVCM_NAME(Segment_Number)
          call READ_MVCM_DATA(Segment_Number)
         endif

      end if

      !--- VIIRS GAC data 
      if (index(Sensor%Sensor_Name,'VGAC') > 0) then
         call READ_VGAC_DATA(Segment_Number, Ierror_Level1b)
         ! If error reading, then go to next file
         if (Ierror_Level1b /= 0) return
      endif

      !--- read EPS-SG data
      if (index(Sensor%Sensor_Name,'METIMAGE') > 0) then
         call READ_EPS_SG_DATA(Segment_Number, Ierror_Level1b)
         ! If error reading, then go to next file
         if (Ierror_Level1b /= 0) return
      endif
      
   end subroutine READ_LEVEL1B_DATA

   !-------------------------------------------------------------------
   ! determine the name of the static nav file for GOES RU data
   !-------------------------------------------------------------------
   subroutine COMPUTE_GOES_RU_STATIC_NAV_FILE_NAME()
      character(len=10):: Orbital_Slot, Scene_Id

      call READ_NETCDF_GLOBAL_ATTRIBUTE(Image%Level1b_Full_Name, 'orbital_slot', Orbital_Slot)
      call READ_NETCDF_GLOBAL_ATTRIBUTE(Image%Level1b_Full_Name, 'scene_id', Scene_Id)
          
      if (trim(Orbital_Slot)=="GOES-East") then
         if (trim(Scene_Id)=="Full Disk") then
             Static_Nav_File = "goes_east_abi_fulldisk_static_nav.nc"
         elseif (trim(Scene_Id)=="CONUS") then
             Static_Nav_File = "goes_east_abi_conus_static_nav.nc"
         elseif (trim(Scene_Id)=="Mesoscale") then
             Static_Nav_File = "goes_east_abi_meso_static_nav.nc"
         else
             Static_Nav_File = "no_file"
         endif
      elseif (trim(Orbital_Slot)=="GOES-Test") then
         if (trim(Scene_Id)=="Full Disk") then
             Static_Nav_File = "goes_central_abi_fulldisk_static_nav.nc"
         elseif (trim(Scene_Id)=="CONUS") then
             Static_Nav_File = "goes_central_abi_conus_static_nav.nc"
         else
             Static_Nav_File = "no_file"
         endif
      elseif (trim(Orbital_Slot)=="GOES-West") then
         if (trim(Scene_Id)=="Full Disk") then
             Static_Nav_File = "goes_west_abi_fulldisk_static_nav.nc"
         elseif (trim(Scene_Id)=="CONUS") then
             Static_Nav_File = "goes_west_abi_conus_static_nav.nc"
         else
             Static_Nav_File = "no_file"
         endif
      else
          Static_Nav_File = "no_file"
      endif

   end subroutine COMPUTE_GOES_RU_STATIC_NAV_FILE_NAME
   !--------------------------------------------------------------------------------------------------
   !
   !--------------------------------------------------------------------------------------------------
   subroutine SET_SENSOR_CHANNEL_MAPPING()

      select case(Sensor%WMO_Id)

      case(3:5,200:209,223,706:708) !AVHRR
         Sensor%Num_Chan_Sensor = 6
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = [1,2,6,20,31,32]

      case(55:57,70) !MSG-SEVIRI
         Sensor%Num_Chan_Sensor = 11
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = [1,2,6,20,27,28,29,30,31,32,33]
      
      case(171:172,810,514:515) !MTSAT,COMS,Fy2D,Fy2E
         Sensor%Num_Chan_Sensor = 5
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = [1,20,27,31,32]
      case(252:259) !GOES-I/P
         Sensor%Num_Chan_Sensor = 6
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = [1,20,27,31,32,33]
      case (173:174) !-- ahi
         Sensor%Num_Chan_Sensor = 16
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = (/3,4,1,2,6,7,20,37,27,28,29,30,38,31,32,33/)
      case(270:271) !-- abi
         Sensor%Num_Chan_Sensor = 16
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = (/3,1,2,26,6,7,20,37,27,28,29,30,38,31,32,33/)
      case(224:225) !VIIRS - SNPP
         Sensor%Num_Chan_Sensor = 16
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = (/8,9,3,4,1,15,2,5,26,6,7,20,22,29,31,32/)
      case(384) !ESP-SG - METIMAGE TODO Change to real WMO ID when available
         Sensor%Num_Chan_Sensor = 20
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = (/3,4,1,15,45,2,18,5,26,6,7,20,22,23,27,28,29,31,32,33/)
      case(530) !FY4A - AGRI
         Sensor%Num_Chan_Sensor = 14
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = (/3,1,2,26,6,7,21,20,27,28,29,31,32,33/)
      case(783,784) !MODIS 
         Sensor%Num_Chan_Sensor = 36
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, &
                                   20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36]
      case(523) !FY3D
         Sensor%Num_Chan_Sensor = 21
         if (.not. allocated(Sensor%CLAVRx_Chan_Map)) allocate(Sensor%Clavrx_Chan_Map(Sensor%Num_Chan_Sensor))
         Sensor%CLAVRx_Chan_Map = [3,4,1,2,26,6,7,8,9,10,15,17,18,19,5, &
                                   20,23,28,29,31,32]
      case default
         print*,'sensor for WMO number not found in SET_SENSOR_CHANNEL_MAPPING for ', Sensor%WMO_ID
         print*,'stopping ... Please fix this in sensor_mod.F90'
         print*,' better tell andi.walther@ssec.wisc.edu'
         stop    
      end select

   end subroutine SET_SENSOR_CHANNEL_MAPPING

end module SENSOR_MOD
