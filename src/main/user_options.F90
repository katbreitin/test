! $Id: user_options.f90 4076 2021-01-26 20:34:35Z dbotambekov $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: user_options.f90 (src)
!       USER_OPTIONS (program)
!
! PURPOSE: CLAVR-x Module to house routines dealing with user options from the
!          default_options file or the command line
!
! DESCRIPTION:
!             Public Routines in this module and their purpose:
!               SETUP_USER_DEFINED_OPTIONS:
!                    Reads user defined configuration
!           	      Called in process_clavrx.f90 outside any loop in the beginning
!
!                UPDATE_CONFIGURATION
!                     Updates configuarion for each file
!                     called in process_clavrx.f90 inside file loop
!                     Check algorithm modes and channel switches
!
!
! AUTHORS:
!  	Andrew Heidinger, Andrew.Heidinger@noaa.gov
!   	Andi Walther
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
!
!  HISTORY:
!        16 Dec 2014: Switch to new option file  (AW)
!                       code cleaning
!
!       30 Dec 2014: Submitting to trunk
!
!
!--------------------------------------------------------------------------------------
module USER_OPTIONS

   use CX_REAL_BOOLEAN_MOD

   use PIXEL_COMMON_MOD, only: &
        Sensor &
      , Geo &
      , Nav &
      , Ch &
      , Image  &
      , ACHA &
      , CCL &
      , ASOS &
      , NWP_PIX &
      , Cld_Flag &
      , Use_Aux_Flag &
      , Cloud_Mask_Bayesian_Flag &
      , Sasrab_Flag &
      , Output_Format_Flag &
      , Nav_Opt &
      , Ref_Cal_1b &
      , Rtm_Opt &
      , Temporary_Data_Dir &
      , Therm_Cal_1b &
      , Ancil_Data_Dir &
      , Cfsr_Data_Dir &
      , Merra_Data_Dir &
      , Gdas_Data_Dir &
      , Erai_Data_Dir &
      , File_List &
      , Level2_List &
      , Gfs_Data_Dir &
      , Level2_File_Flag &
      , Lrc_Flag &
      , Modis_Clr_Alb_Flag &
      , Ncep_Data_Dir &
      , Oisst_Data_Dir &
      , Output_Scaled_Reflectances &
      , Process_Undetected_Cloud_Flag &
      , Read_Coast_Mask &
      , Read_Hires_Sfc_Type &
      , Read_Land_Mask &
      , Read_Snow_Mask &
      , Read_Surface_Elevation &
      , Read_Volcano_Mask &
      , Snow_Data_Dir &
      , Use_Default &
      , Sfc_Emiss_Option &
      , Use_Sea_IR_Emiss &
      , Bayesian_Cloud_Mask_Name &
      , Compress_Flag &
      , Cloud_Mask_Mode &
      , Nlcomp_Mode &
      , Dcomp_Mode &
      , Aerosol_Mode &
      , Avhrr_1_flag &
      , Read_Dark_Comp &
      , Dark_Comp_Data_Dir &
      , Globsnow_Data_Dir &
      , Use_IR_Cloud_Type_Flag &
      , Caliop_Dir &
      , Caliop_Flag &
      , Verbose_Level_Flag &
      , Nucaps_Flag &
      , Static_Dark_Sky_Flag &
      , X_Sample_Offset &
      , Y_Sample_Offset &
      , Elem_Abs_Idx_ACHA_Dump &
      , Line_Abs_Idx_ACHA_Dump &
      , Use_Iband &
      , WMO_Id_ISCCPNG

   use CONSTANTS_MOD, only: &
      Sym &
      , Exe_Prompt &
      , Nchan_Clavrx &
      , MISSING_VALUE_INT4, int1

   use FILE_UTILS, only: &
      Get_Lun

   use LEVEL2B_MOD, only: &
      INIT_RANDOM_SEED

   use  CLAVRX_MESSAGE_MOD, only: &
      Mesg &
    , Mesg_1i &
    , Verb_Lev

   use AWG_CLOUD_HEIGHT, only: Num_ACHA_Modes, ACHA_Mode_Max_Length, ACHA_Mode_Values

   implicit none

   private

   public  :: SETUP_USER_DEFINED_OPTIONS
   public  :: UPDATE_CONFIGURATION

   private :: SET_DEFAULT_VALUES
   private :: READ_OPTION_FILE
   private :: DETERMINE_USER_CONFIG
   private :: QC_CLAVRXORB_OPTIONS
   private :: HELPER
   private :: DEFAULT_ACHA_MODE
   private :: DEFAULT_DCOMP_MODE
   private :: DEFAULT_NB_MASK_CLASSIFIER_FILE
   private :: CHECK_DCOMP_ALGORITHM_CHOICES
   private :: CHECK_ACHA_ALGORITHM_CHOICES
   private :: EXISTING_CHANNELS
   private :: CHANNEL_SWITCH_ON
   private :: SUB_PIXEL_CHANNEL_ON_SET
   private :: CHECK_USER_CHANNEL_CHOICES
   private :: EXPERT_MODE_CHANNEL_ALGORITHM_CHECK
   private :: DEFAULT_NBM_MASK_CLASSIFIER_FILE

   character(24), parameter, private :: MOD_PROMPT = " USER_OPTIONS_ROUTINES: "
   character ( len = 1020 ) :: Data_Base_Path
   integer :: DCOMP_Mode_User_Set
   character (len = ACHA_Mode_Max_Length)  :: ACHA_Mode_User_Set
   integer :: NLCOMP_Mode_User_Set
   integer :: Expert_Mode

   integer(kind=int1) :: Chan_On_Flag_Default_User_Set (Nchan_Clavrx)

   ! ---------------------------------------------------------------------------------
   ! Default Algorithm Modes -
   ! ---------------------------------------------------------------------------------
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_Avhrr = '110_120'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_Avhrr1 = '110'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_Avhrr_Fusion = '110_120_133'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_Avhrr1_Fusion = '110_133'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_Goes_IL = '067_110_120'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_Goes_MP = '067_110_133'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_Goes_SNDR = '067_110_133'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_Goes_RU = '110_120_133'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_COMS = '067_110_120'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_VIIRS = '085_110_120'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_MTSAT = '067_110_120'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_SEVIRI = '110_120_133'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_Modis = '110_120_133'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_FY2 = '067_110_120'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_FY4A_AGRI = '110_120_133'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_AHI = '110_120_133'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_Mersi2 = '085_110_120'
   character(len=ACHA_Mode_Max_Length),parameter:: ACHA_Mode_Default_METIMAGE = '110_120_133'

contains

   ! ---------------------------------------------------------------------------------
   !  wrapper for initial clavrx option read
   ! ---------------------------------------------------------------------------------
   subroutine SETUP_USER_DEFINED_OPTIONS()
      call SET_DEFAULT_VALUES
      call DETERMINE_USER_CONFIG()
      call QC_CLAVRXORB_OPTIONS()
   end subroutine SETUP_USER_DEFINED_OPTIONS

   ! ---------------------------------------------------------------------------------
   !
   ! ---------------------------------------------------------------------------------
   subroutine SET_DEFAULT_VALUES

      Modis_Clr_Alb_Flag = 1 ! use clear-sky MODIS albedo maps
      Output_Scaled_Reflectances = sym%NO !default is to output ref / cossolzen

      !--- default solar zenith limits
      Geo%Solzen_Min_Limit= 0
      Geo%Solzen_Max_Limit= 180.0
      Geo%Satzen_Min_Limit= 0
      Geo%Satzen_Max_Limit= 85.0

      !--- offset when sampling subpixels
      X_Sample_Offset = 0
      Y_Sample_Offset = 0
      Use_Iband = .false.

      !--- ACHA single pixel dump location
      Elem_Abs_Idx_ACHA_Dump = 0
      Line_Abs_Idx_ACHA_Dump = 0

      !--- default what can be changed for expert mode
      Verbose_Level_Flag = verb_lev % DEFAULT
      Cloud_Mask_Bayesian_Flag = 1
      Cloud_Mask_Mode = 'ecm2'
      Dcomp_Mode_User_Set = 3
      Acha_Mode_User_Set = 'default'
      CCL%Mode = 1
      CCL%Type = 0
      ASOS%Mode = 0_int1
      Nlcomp_Mode = 1
      Level2_File_Flag = 1
      Output_Format_Flag = 0
      Cld_Flag = 1
      Image%Number_Of_Lines_Per_Segment = 240
      Sasrab_Flag = 0
      NWP_PIX%Nwp_Opt = 1
      NWP_PIX%Mode = 0_int1
      Rtm_Opt = 1
      Compress_Flag = 1
      Use_Aux_Flag = 0
      bayesian_cloud_mask_name = 'default'
      sfc_emiss_option = 1
      Use_Sea_IR_Emiss = 1
      Read_Hires_Sfc_Type = 1
      Read_Land_Mask = 1
      Read_Coast_Mask = 1
      Read_Surface_Elevation = 1
      Read_Volcano_Mask = 0
      Read_Snow_Mask = 1
      Read_Dark_Comp = 0
      Ref_Cal_1b = 0
      Therm_Cal_1b = 0
      Nav_Opt = 0
      Lrc_Flag = 1
      NWP_PIX%Smooth_Nwp_Flag = 1
      Process_Undetected_Cloud_Flag = 0
      Chan_On_Flag_Default_User_Set(1:6) = int([1,1,1,1,1,1],kind=int1)
      Chan_On_Flag_Default_User_Set(7:12) = int([1,1,1,1,1,1],kind=int1)
      Chan_On_Flag_Default_User_Set(13:18) = int([1,1,1,1,1,1],kind=int1)
      Chan_On_Flag_Default_User_Set(19:24) = int([1,1,1,1,1,1],kind=int1)
      Chan_On_Flag_Default_User_Set(25:30) = int([1,1,1,1,1,1],kind=int1)
      Chan_On_Flag_Default_User_Set(31:36) = int([1,1,1,1,1,1],kind=int1)
      Chan_On_Flag_Default_User_Set(37:42) = int([1,1,0,0,0,0],kind=int1)
      Chan_On_Flag_Default_User_Set(43:45) = int([0,1,1],kind=int1)
      Nav%Lat_North_Limit = 90.0
      Nav%Lat_South_Limit = -90.0
      Nav%Lon_East_Limit = 180.0
      Nav%Lon_West_Limit = -180.0
      Geo%Satzen_Max_Limit = 85.0
      Geo%Satzen_Min_Limit = 0.0
      Geo%Solzen_Max_Limit = 180.0
      Geo%Solzen_Min_Limit = 0.0
      Image%Chan_Average_Flag = 0
      Image%X_Stride = 1
      Image%Y_Stride = 1
      Nav%Limit_Flag = 0
      Nav%Lon_Lat_Limits_Set = .false.
      Nav%Lat_North_Limit = 90.0
      Nav%Lat_South_Limit = -90.0
      Nav%Lon_East_Limit = 180.0
      Nav%Lon_West_Limit = -180.0
      Geo%Satzen_Max_Limit = 85.0
      Geo%Satzen_Min_Limit = 0.0
      Geo%Solzen_Max_Limit = 180.0
      Geo%Solzen_Min_Limit = 0.0
      Nav%Domain_Name = "UNKNOWN"

   end subroutine SET_DEFAULT_VALUES


   !---------------------------------------------------------------------------------
   ! Read Parameters from AVHRR INPUT files and check them for errors
   !
   ! parameters are pased in avhrr_pixel_common_mod public memory
   !---------------------------------------------------------------------------------
   subroutine READ_OPTION_FILE (File_Default)

#if defined(FC_INTEL)
      use ifport, only: hostnam, getpid
#endif

      character(len=*), intent(in):: File_Default
      integer::ios0
      integer::erstat
      integer:: Default_Lun
      integer:: PID
      character(len=7):: Pid_String
      character(len=1020):: Temporary_Data_Dir_Root
      character(len=1020):: hostname
      integer:: String_Length
      character(len=1):: Last_Char

      call MESG ("Option file to be read in: "//trim(File_Default),level = verb_lev % DEFAULT)

      Default_Lun = GET_LUN()

      open(unit=Default_Lun,file=trim(File_Default), &
          iostat = ios0, &
          action="read", &
          status="old", &
          position="rewind")

      erstat = 0
      if (ios0 /= 0) then    !avhrr_input_check
         erstat = 1
         call MESG(trim(EXE_PROMPT)//"error opening clavrx_options file",color=1,level=verb_lev % ERROR )
         stop 1
      end if

      read(unit=Default_Lun,fmt="(a)") Data_Base_Path
      read(unit=Default_Lun,fmt="(a)") Temporary_Data_Dir_Root
      read(unit=Default_Lun,fmt=*) Expert_Mode

      !--- add the PID to Temporary_Data_Dir
      string_length = len_trim(Temporary_Data_Dir_Root)
      last_char = Temporary_Data_Dir_Root(string_length:string_length)

      if (index(last_char, '/') > 0) then
        Temporary_Data_Dir_Root = Temporary_Data_Dir_Root(1:string_length-1)
      endif

#if defined(FC_GNU)
      call HOSTNM(hostname)
#elif defined(FC_INTEL)
      erstat = HOSTNAM(hostname)
#endif

      PID = getpid()
      write(Pid_String,'(I7.7)' ) pid
      Temporary_Data_Dir = trim(Temporary_Data_Dir_Root) // '_' // trim(hostname) // '_' // trim(Pid_String) // '/'

      !--- check expert mode, if 0 return
      if ( Expert_Mode  == 0 )  then
          close(unit=Default_Lun)
          return
      end if

      !--- contintue reading file if selected expert mode allows
      read(unit=Default_Lun,fmt=*) Verbose_Level_Flag
      read(unit=Default_Lun,fmt=*) Cloud_Mask_Bayesian_Flag
      read(unit=Default_Lun,fmt=*) Dcomp_Mode_User_Set
      read(unit=Default_Lun,fmt=*) ACHA_Mode_User_Set
      Acha_Mode_User_Set = trim(adjustl(Acha_Mode_User_Set))
      ACHA%Mode_User = ACHA_Mode_User_Set
      read(unit=Default_Lun,fmt=*) CCL%Mode
      read(unit=Default_Lun,fmt=*) CCL%Type
      read(unit=Default_Lun,fmt=*) ASOS%Mode

      read(unit=Default_Lun,fmt=*) Nlcomp_Mode
      read(unit=Default_Lun,fmt=*) Aerosol_Mode

      if ( Expert_Mode <= 1 )  then
          close(unit=Default_Lun)
          return
      end if

      read(unit=Default_Lun,fmt=*) Level2_File_Flag
      read(unit=Default_Lun,fmt=*) Output_Format_Flag
      read(unit=Default_Lun,fmt=*) Cld_Flag
      read(unit=Default_Lun,fmt=*) Image%Number_Of_Lines_Per_Segment
      read(unit=Default_Lun,fmt=*) Sasrab_Flag
      Dark_Comp_Data_Dir=trim(Data_Base_Path)//'/dynamic/goes_dark_sky_composites/'
      if ( Sasrab_Flag == 1 )  then
          Read_Dark_Comp=1
      end if

      read(unit=Default_Lun,fmt=*) NWP_PIX%Nwp_Opt
      read(unit=Default_Lun,fmt=*) NWP_PIX%Mode

      read(unit=Default_Lun,fmt=*) Rtm_Opt
      read(unit=Default_Lun,fmt=*) Nav_Opt

      read(unit=Default_Lun,fmt=*) Compress_Flag
      read(unit=Default_Lun,fmt=*) Use_Aux_Flag

      read(unit=Default_Lun,fmt="(a)") Bayesian_Cloud_Mask_Name

      if ( Expert_Mode <= 2 )  then
          close(unit=Default_Lun)
          return
      end if

      read(unit=Default_Lun,fmt=*) sfc_emiss_option
      read(unit=Default_Lun,fmt=*) Use_Sea_IR_Emiss
      read(unit=Default_Lun,fmt=*) Read_Hires_Sfc_Type
      read(unit=Default_Lun,fmt=*) Read_Land_Mask
      read(unit=Default_Lun,fmt=*) Read_Coast_Mask
      read(unit=Default_Lun,fmt=*) Read_Surface_Elevation
      read(unit=Default_Lun,fmt=*) Read_Volcano_Mask
      read(unit=Default_Lun,fmt=*) Read_Snow_Mask
      read(unit=Default_Lun,fmt=*) Read_Dark_Comp

      if ( Expert_Mode <= 3 ) then
          close(unit=Default_Lun)
          return
      end if

      read(unit=Default_Lun,fmt=*) Ref_Cal_1b
      read(unit=Default_Lun,fmt=*) Therm_Cal_1b

      if ( Expert_Mode <= 4 ) then
          close(unit=Default_Lun)
          return
      end if

      read(unit=Default_Lun,fmt=*) Lrc_Flag
      read(unit=Default_Lun,fmt=*) NWP_PIX%Smooth_Nwp_Flag
      read(unit=Default_Lun,fmt=*) Process_Undetected_Cloud_Flag

      if ( Expert_Mode <= 5 ) then
          close(unit=Default_Lun)
          return
      end if

      ! --- Read lat, lon and sun angle high - low limits
      read(unit=Default_Lun,fmt=*) Nav%Limit_Flag

      if (Nav%Limit_Flag == 0) then
         Nav%Lon_Lat_Limits_Set = .false.
         Nav%Lat_North_Limit = 90.0
         Nav%Lat_South_Limit = -90.0
         Nav%Lon_East_Limit = 180.0
         Nav%Lon_West_Limit = -180.0
         Geo%Satzen_Max_Limit = 85.0
         Geo%Satzen_Min_Limit = 0.0
         Geo%Solzen_Max_Limit = 180.0
         Geo%Solzen_Min_Limit = 0.0
         Nav%Domain_Name = "UNKNOWN"
      else
         backspace(unit=Default_Lun)
         Nav%Lon_Lat_Limits_Set = .true.
         read(unit=Default_Lun,fmt=*) Nav%Limit_Flag, &
                                      Nav%Lat_South_Limit, Nav%Lat_North_Limit, &
                                      Nav%Lon_West_Limit, Nav%Lon_East_Limit, &
                                      Geo%Satzen_Min_Limit, Geo%Satzen_Max_Limit, &
                                      Geo%Solzen_Min_Limit, Geo%Solzen_Max_Limit, &
                                      Nav%Domain_Name
      endif

      !--- constrain values
      Nav%Lat_North_Limit = min(Nav%Lat_North_Limit, 90.0)
      Nav%Lat_South_Limit = max(Nav%Lat_South_Limit, -90.0)
      Nav%Lon_East_Limit = min(Nav%Lon_East_Limit, 180.0)
      Nav%Lon_East_Limit = max(Nav%Lon_East_Limit, -180.0)
      Nav%Lon_West_Limit = min(Nav%Lon_West_Limit, 180.0)
      Nav%Lon_West_Limit = max(Nav%Lon_West_Limit, -180.0)
      Geo%Satzen_Max_Limit = min(Geo%Satzen_Max_Limit, 85.0)
      Geo%Satzen_Min_Limit = max(Geo%Satzen_Min_Limit, 0.0)
      Geo%Solzen_Max_Limit = min(Geo%Solzen_Max_Limit, 180.0)
      Geo%Solzen_Min_Limit = max(Geo%Solzen_Min_Limit, 0.0)


      !---- Read input on handling of native res processing and striding
      read(unit=Default_Lun,fmt=*) Image%Chan_Average_Flag, Image%X_Stride, Image%Y_Stride

      !--- Read channel on/off flags
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(1:6)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(7:12)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(13:18)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(19:24)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(25:30)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(31:36)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(37:42)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(43:45)

      !--- read in wmo id for selected ISSCP-NG L1g sensor (can be missing)
      read(unit=Default_Lun,fmt=*,iostat = erstat) WMO_Id_ISCCPNG
      if (erstat /= 0) then
         WMO_Id_ISCCPNG = -999
      endif

      !--- close the options file
      close(unit=Default_Lun)

   end subroutine READ_OPTION_FILE
   !------------------------------------------------------------------
   ! Command-line input option variables
   !------------------------------------------------------------------
   subroutine DETERMINE_USER_CONFIG()

      character(len=30) :: fargv
      character(len=30) :: junk
      character(len=1) :: temp_string
      logical:: back
      character(len=1020):: default_temp
      real :: int_temp
      integer :: fargc
      integer :: i
      integer :: Temp_Scans_Arg !--- temporary integer for number of scanlines
      integer :: iargc

      temp_string = '.'  !--- temporary string to search for in angle commandline

      !---- SET DEFAULT OPTIONS

      Use_Default = sym%YES
      Default_Temp="./clavrx_options"
      File_List = "./file_list"
      Level2_List = "./level2_list"
      Use_IR_Cloud_Type_Flag = .false.      !controls which type is used in algorithms

      Temp_Scans_Arg = 0
      fargc = iargc()

      !--- first we will check to see if the default file is used
      !--- also check to see if the help file is to be displayed

      do i=1, fargc
         call getarg(i,fargv)
         if (trim(fargv) == "-help" .or. &
               trim(fargv) == "-h") then
            call HELPER()
            stop
         else if  ( trim(fargv) == "-version" .or. &
             trim(fargv) == "-ver") then
            print*,&
            &'$Header: https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/trunk/main_src/user_options.f90 4076 2021-01-26 20:34:35Z dbotambekov $'
            stop

         else if (trim(fargv) == "-no_default") then
            Use_Default = sym%NO
            !Different default file used
         else if (trim(fargv) == "-default") then
            call getarg(i+1,default_temp)
            default_temp=trim(default_temp)
         end if
      end do


      !---- If the default file is used, read it in first, then
      !---- check for other command line changes to the options
      if(Use_Default == sym%YES)  then
         call READ_OPTION_FILE (Default_Temp)
      end if

      if(Use_Default == sym%NO)  then
          print *, EXE_PROMPT, "Using standard defaults and command line options"
      end if

      !--- Initialize caliop flag
      Caliop_Flag = .false.
      Nucaps_Flag = .false.
      Static_Dark_Sky_Flag = .false.

      do i=1, fargc

        call getarg(i,fargv)

        !Change Ref_Cal_1b flag
        if (trim(fargv) == "-Ref_Cal_1b") then
          Ref_Cal_1b = sym%YES
        elseif (trim(fargv) == "-no_Ref_Cal_1b") then
          Ref_Cal_1b = sym%NO

        !Change therm_Cal_1b flag
        elseif(trim(fargv) == "-therm_Cal_1b") then
          therm_Cal_1b = sym%YES
        elseif(trim(fargv) == "-no_Therm_Cal_1b") then
          therm_Cal_1b = sym%NO

        !Change Nav type
        elseif(trim(fargv) == "-l1bnav") then
          nav_opt = 0

        !Change level2 output flag
        elseif(trim(fargv) == "-level2_file") then
          level2_file_Flag = sym%YES
        elseif(trim(fargv) == "-no_level2_file") then
          level2_file_Flag = sym%NO

        !Change cloud mask
        elseif(trim(fargv) == "-use_aux") then
          Use_Aux_Flag = sym%YES
        elseif(trim(fargv) == "-no_aux") then
          Use_Aux_Flag = sym%NO

        !Change data compression flag
        elseif(trim(fargv) == "-no_output_comp") then
          Compress_Flag=0
        elseif(trim(fargv) == "-output_comp_gzip") then
          Compress_Flag=1
        elseif(trim(fargv) == "-output_comp_szip") then
          Compress_Flag=2

        ! - change clear-sky trans
        elseif(trim(fargv) == "-rttov") then
          rtm_opt = 2

        !Change lat max/min for processing
        elseif(trim(fargv) == "-lat_south_limit") then
           call getarg(i+1,junk)
           back = .true.
           int_temp = scan(junk,temp_string, back)
           if(int_temp .gtr. 0.0) read(junk,'(f6.3)') Nav%Lat_South_Limit
           if(int_temp .eqr. 0.0) read(junk,'(f6.0)') Nav%Lat_North_Limit
        elseif(trim(fargv) == "-lat_north_limit") then
          call getarg(i+1,junk)
           back = .true.
          int_temp = scan(junk,temp_string, back)
          if(int_temp .gtr. 0.0) read(junk,'(f6.3)') Nav%Lat_North_Limit
          if(int_temp .eqr. 0.0) read(junk,'(f6.0)') Nav%Lat_North_Limit

        !Change ancillary data directory
        elseif(trim(fargv) == "-ancil_data_dir") then
          call getarg(i+1,ancil_data_dir)
          ancil_data_dir=trim(ancil_data_dir)

        !Smooth/not smooth NWP data
        elseif(trim(fargv) == "-smooth_nwp") then
          NWP_PIX%Smooth_Nwp_Flag = sym%YES
        elseif(trim(fargv) == "-no_smooth_nwp") then
          NWP_PIX%Smooth_Nwp_Flag = sym%NO

        !Read/not read volcano mask
        elseif(trim(fargv) == "-read_volcano_mask") then
          read_volcano_mask = sym%YES
        elseif(trim(fargv) == "-no_volcano_mask") then
          read_volcano_mask = sym%NO

        !Output scaled reflectances or not
        elseif(trim(fargv) == "-output_scaled_ref") then
          output_scaled_reflectances = sym%YES

        elseif(trim(fargv) == "-no_output_scaled_ref") then
          output_scaled_reflectances = sym%NO

        ! change the number of scanlines per segment
        elseif(trim(fargv) == "-lines_per_seg") then
          call getarg(i+1,junk)
          read(junk,'(i4)') Temp_Scans_Arg
          if(Temp_Scans_Arg > 1) read(junk,'(i4)') Image%Number_of_Lines_Per_Segment

        !Change solar zenith angle limits
        elseif(trim(fargv) == "-solzen_min_limit") then
          call getarg(i+1,junk)
          int_temp = scan(junk,temp_string, back)
          back = .true.
          if(int_temp .gtr. 0.0) read(junk,'(f6.3)') Geo%Solzen_Min_Limit
          if(int_temp .eqr. 0.0) read(junk,'(f6.0)') Geo%Solzen_Min_Limit

        !Change dcomp mode
         elseif(trim(fargv) == "-DCOMP_Mode") then
          call getarg(i+1,junk)
          read(junk,'(i1)') Dcomp_Mode_User_Set

        elseif(trim(fargv) == "-solzen_max_limit") then
           call getarg(i+1,junk)
           back = .true.
           int_temp = scan(junk,temp_string, back)
           if(int_temp .gtr. 0.0) read(junk,'(f6.3)') Geo%Solzen_Max_Limit
           if(int_temp .eqr. 0.0) read(junk,'(f6.0)') Geo%Solzen_Max_Limit

        elseif (trim(fargv) == "-file_list") then
           call getarg(i+1,File_List)
           File_List=trim(File_List)

        elseif (trim(fargv) == "-level2_list") then
           call getarg(i+1,Level2_List)
           Level2_List=trim(Level2_List)

        elseif (trim(fargv) == "-use_ir_type") then
           Use_IR_Cloud_Type_Flag = .true.

        elseif (trim(fargv) == "-caliop") then
           call getarg(i+1,Caliop_Dir)
           Caliop_Flag = .true.

        elseif (trim(fargv) == "-nucaps") then
           Nucaps_Flag = .true.

        elseif (trim(fargv) == "-static_dark_sky") then
           Static_Dark_Sky_Flag = .true.

        elseif (trim(fargv) == "-x_off") then
           call getarg(i+1,junk)
           read(junk,'(I2)') X_Sample_Offset

        elseif (trim(fargv) == "-y_off") then
           call getarg(i+1,junk)
           read(junk,'(I2)') Y_Sample_Offset

        elseif (trim(fargv) == "-acha_x") then
           call getarg(i+1,junk)
           read(junk,'(I5)') Elem_Abs_Idx_ACHA_Dump
           print *, "acha_x ", Elem_Abs_Idx_ACHA_Dump

        elseif (trim(fargv) == "-acha_y") then
           call getarg(i+1,junk)
           read(junk,'(I5)') Line_Abs_Idx_ACHA_Dump
           print *, "acha_y ", Line_Abs_Idx_ACHA_Dump

        elseif (trim(fargv) == "-use_iband") then
           Use_Iband = .true.

        endif
      enddo

      !--- default ancillary data directory
      Ancil_Data_Dir = trim(Data_Base_Path)
      Gfs_Data_Dir = trim(Data_Base_Path)//'/dynamic/gfs/'
      if (NWP_PIX%Nwp_Opt == 8) Gfs_Data_Dir = trim(Data_Base_Path)//'/dynamic/gfs_fv3/'
      if (NWP_PIX%Nwp_Opt == 7) Gfs_Data_Dir = trim(Data_Base_Path)//'/dynamic/gfs_ait/'
      Ncep_Data_Dir = trim(Data_Base_Path)//'/dynamic/ncep-reanalysis/'
      Cfsr_Data_Dir = trim(Data_Base_Path)//'/dynamic/cfsr/'
      Merra_Data_Dir = trim(Data_Base_Path)//'/dynamic/merra/'
      Gdas_Data_Dir = trim(Data_Base_Path)//'/dynamic/gdas/'
      Erai_Data_Dir = trim(Data_Base_Path)//'/dynamic/erai/'
      Oisst_data_Dir = trim(Data_Base_Path)//'/dynamic/oisst_nc/'
      Snow_Data_Dir = trim(Data_Base_Path)//'/dynamic/snow/hires/'
      Globsnow_Data_Dir = trim(Data_Base_Path)//'/dynamic/snow/globsnow/'

      Sensor%Chan_On_Flag_Default = Chan_On_Flag_Default_User_Set

   end subroutine DETERMINE_USER_CONFIG

   !-------------------------------------------------------------------------------
   !--- QC options and modify as needed
   !-------------------------------------------------------------------------------
   subroutine QC_CLAVRXORB_OPTIONS()

      !---- Since the NWP controls everything, we first check if an NWP is being used
      !---- before anything else is checked.  If no nwp, we stop processing

      select case ( NWP_PIX%Nwp_Opt)
      case ( 0 )
         print *,  EXE_PROMPT, "No choice made for NWP data, will not run algoritms or orbital level3 files"
         Cld_Flag = sym%NO

         Sasrab_Flag = sym%NO
         Cloud_Mask_Bayesian_Flag = sym%NO
         Use_Aux_Flag = sym%NO ! this is to determine if the lut's are being read in
      case ( 1 )
         call MESG ("GFS data will be used",level = verb_lev % DEFAULT)
      case ( 2 )
         call MESG ( "NCEP Reanalysis data will be used",level = verb_lev % DEFAULT)
      case ( 3 )
         call MESG ( "NCEP Climate Forecast System Reanalysis data will be used",level = verb_lev % DEFAULT)
      case ( 4 )
         call MESG ( "GDAS Reanalysis data will be used",level = verb_lev % DEFAULT)
      case ( 5 )
         call MESG ( "MERRA Reanalysis data will be used",level = verb_lev % DEFAULT)
      case ( 6 )
         call MESG ( "ERA Interim Reanalysis data will be used",level = verb_lev % DEFAULT)
      case ( 7 )
         call MESG ( "GFS with AIT Interpolation will be used",level = verb_lev % DEFAULT)
      case ( 8 )
         call MESG ( "GFS FV3 will be used",level = verb_lev % DEFAULT)
      case default
         print *,  EXE_PROMPT, "unrecognized value for Nwp_Opt: ", NWP_PIX%Nwp_Opt
         stop "6-Nwp_Flag"

      end select

      if (Cloud_Mask_Bayesian_Flag == sym%YES) then
         call MESG  ("Bayesian cloud mask will be generated",level = 9)
      endif

      if (Ref_Cal_1b == sym%YES) then
         call MESG ("Reflectance Calibration within 1b will be used",level = 9)
      endif

      if (Therm_Cal_1b == sym%YES) then
         call MESG ("Thermal Calibration within 1b will be used",level = 9)
      endif

      if (Nav_Opt == 1) then
         call MESG ("CLEVERNAV geolocation no longer supported, using REPOSNX",level = 9)
         nav_Opt = 2
       endif

      if (Nav_Opt == 2) then
         call MESG( "REPOSNX geolocation adjustment done")
      endif

      if (Use_Aux_Flag == sym%YES) then
         call MESG("Cloud mask results will be read in from an aux file", level = 5)
      endif

      if (Rtm_Opt /= 1 .and. Rtm_Opt /= 2) then
         call MESG("Only PFAAST and RTTOV RTM are implemented, switching to rtm ==2 (RTTOV)", level = 5, color = 6)
         rtm_opt = 2
      endif

      if (Nav%Limit_Flag == sym%YES) then
        call MESG ("Static navigation will be used")
      endif


   end subroutine QC_CLAVRXORB_OPTIONS


   !--------------------------------------------------------------------------
   ! This subroutine outputs what each command line options
   ! are available to override the default file options or to set
   ! options different from the default options
   !--------------------------------------------------------------------------
   subroutine HELPER()

      print "(a,'help: option list')",EXE_PROMPT
      print *,"  This is a list of all of the command line options to override"
      print *,"     the file list options"
      print *," "

      print *,"  -default (default_temp)"
      print *, "  This option allows you to set which default file you want to use."
      print *, "  Initial default file is clavrx_options"
      print *," "

      print *,"  -file_list (file_list)"
      print *, "  This option allows you to set which file_list file you want to use."
      print *, "  Initial default file is namd file_list in the current directory"
      print *," "

      print *,"  -level2_list (level2_list)"
      print *, "  This option allows you to set which level2_list file you want to use."
      print *, "  Initial default file is namd level2_list in the current directory"
      print *," "

      print *,"  -lines_per_seg (Imager%Number_Of_Lines_Per_Segment)"
      print *, "  specify the number of lines per segment"
      print *," "

      print *,"  -ref_Cal_1b"
      print *,"  Use the reflectance cal in level 1b file."
      print *," "

      print *,"  -no_ref_Cal_1b"
      print *,"  AVHRR-ONLY: Do not us the reflectance cal in level 1b file."
      print *," "

      print *,"  -therm_Cal_1b"
      print *,"   AVHRR-ONLY: Use the thermal calibration in the level 1b file"
      print *," "

      print *,"  -no_therm_Cal_1b"
      print *,"   Do not use the thermal calibration in the level 1b file"
      print *," "

      print *,"  -rtm_file"
      print *, "   Output rtm data file. "
      print *," "
      print *,"  -no_rtm_file"
      print *, "   Do not output rtm data file. "
      print *," "

      print *,"  -level2_file"
      print *, "   Make Level-2 output."
      print *," "
      print *,"  -no_level2_file"
      print *, "   Don't make Level-2 output."
      print *," "

      print *,"  -cloud_mask_1b"
      print *, "   Read cloud mask from level 1b file."
      print *," "


      print *,"  -no_nwp"
      print *, "  Do not use NWP data. No algorithms will be processed. Also, no orbital data will be processed"
      print *," "

      print *,"  -gfs_nwp"
      print *, "  Use GFS data for NWP dataset. "
      print *," "

      print *,"  -ncep_nwp"
      print *, "  Use NCEP Reanalysis data"
      print *," "

      print *,"  -lat_south_limit (Lat_South_Limit)"
      print *, "  southern latitude for processing(degrees)"
      print *," "
      print *,"  -lat_north_limit (Lat_North_Limit)"
      print *, "  northern latitude for processing (degrees)"
      print *," "

      print *,"  -ancil_data_dir (ancil_data_dir)"
      print *, "  change the location of the ancillary data directory"
      print *," "

      print *,"  -smooth_nwp"
      print *, "  Smooth the NWP fields. "
      print *," "
      print *,"  -no_smooth_nwp"
      print *, "  Don't smooth the NWP fields. "
      print *," "

      print *,"  -use_seebor"
      print *, "  Use Seebor emissivity dataset. "
      print *," "
      print *,"  -no_seebor"
      print *, "  Don't use Seebor emissivity dataset. "
      print *," "

      print *,"  -no_surface_elevation"
      print *, "  Don't read surface elevation map in. "
      print *," "
      print *,"  -high_surface_elevation"
      print *, "  Use 1km res GLOBE global elevation map. "
      print *," "
      print *,"  -low_surface_elevation"
      print *, "  Use 8km res GLOBE global elevation map. "
      print *," "


      print *,"  -read_volcano_mask"
      print *, "  Read volcano map in. "
      print *," "
      print *,"  -no_volcano_mask"
      print *, "  Don't read volcano map in. "
      print *," "

      print *,"  -Solzen_min_limit (Solzen_min_limit)"
      print *, "  Solar zenith angle minimum limit"
      print *," "

      print *,"  -Solzen_max_limit (Solzen_max_limit)"
      print *, "  Solar zenith angle maximum limit."
      print *," "

      print *,"  -DCOMP_Mode"
      print *, "  dcomp mode 1,2 or 3 (1.6,2.2 or 3.8 micron)."
      print *," "

      print *,"  -caliop /PATH/TO/COLLOCATION/FILES/"
      print *, "  processes only +/- pixels along caliop collocated path"
      print *, " "

  end subroutine HELPER

   ! -----------------------------------------------------------------
   ! -e wrapper for all updating tools for a new file
   !     called from PROCESS_CLAVRX inside file loop
   ! ---------------------------------------------------
   subroutine UPDATE_CONFIGURATION (SensorName)

      character (len=*) , intent(in) :: SensorName

      if ( Expert_Mode == 0 ) then
         DCOMP_Mode_User_Set = Default_DCOMP_Mode ( )
      end if

      !--- turn off channels if not on sensor
      call CHECK_USER_CHANNEL_CHOICES(SensorName)

      !--- for expert modes < 6, turn on channels on sensor
      call CHANNEL_SWITCH_ON (SensorName)

      !--- sub pixel channel set
      call SUB_PIXEL_CHANNEL_ON_SET()

      !--- check if dcomp choice is supported by available channels
      call CHECK_DCOMP_ALGORITHM_CHOICES(SensorName)

      !--- check if acha choice is supported by available channels
      call CHECK_ACHA_ALGORITHM_CHOICES(SensorName)

      !--- select the cloud mask lut name if default is chosen
      if ( Expert_Mode < 2 .or. trim(Bayesian_Cloud_Mask_Name) == 'default' &
           .or. trim(Bayesian_Cloud_Mask_Name) == 'DEFAULT') then

         if (Cloud_Mask_Bayesian_Flag == sym%ECM1) then
            Bayesian_Cloud_Mask_Name = DEFAULT_NB_MASK_CLASSIFIER_FILE ( SensorName )
         elseif (Cloud_Mask_Bayesian_Flag == sym%ECM2) then
            Bayesian_Cloud_Mask_Name = DEFAULT_NBM_MASK_CLASSIFIER_FILE ( SensorName )
         endif

      end if

      !--- final check if algorithms can run
      call EXPERT_MODE_CHANNEL_ALGORITHM_CHECK ( )

      !--- force missing Himawari8/9 HCAST channels turned off
      if ((Image%DB_Flag) .and. (trim(SensorName) == 'AHI')) then
         Sensor%Chan_On_Flag_Default(3) = 0_int1
         Sensor%Chan_On_Flag_Default(4) = 0_int1
      endif

   end subroutine UPDATE_CONFIGURATION

   !----------------------------------------------------------------------
   !  returns default acha mode
   !----------------------------------------------------------------------
   character(len=ACHA_Mode_Max_Length) function DEFAULT_ACHA_MODE ( SensorName )

      character ( len =*) , intent(in) :: SensorName

      select case ( trim(SensorName))
      case ('FCI')
        DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Goes_RU
      case ( 'AVHRR-2')
         DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Avhrr
      case ( 'AVHRR-3')
         DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Avhrr
      case ( 'AVHRR-1')
         DEFAULT_ACHA_MODE = ACHA_Mode_Default_Avhrr1
      case ( 'GOES-MP-IMAGER')
         DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Goes_MP
      case ( 'GOES-IL-IMAGER')
         DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Goes_IL
      case ( 'GOES-IP-SOUNDER')
         DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Goes_SNDR
      case ( 'GOES-RU-IMAGER')
         DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Goes_RU
      case ( 'MTSAT-IMAGER')
          DEFAULT_ACHA_MODE  = ACHA_Mode_Default_MTSAT
      case ('SEVIRI')
         DEFAULT_ACHA_MODE  = ACHA_Mode_Default_SEVIRI
      case ('FY2-IMAGER')
         DEFAULT_ACHA_MODE  =  ACHA_Mode_Default_FY2
      case ('AGRI') ! - FY4A
         DEFAULT_ACHA_MODE  =  ACHA_Mode_Default_FY4A_AGRI
      case ('VIIRS','VIIRS-NASA','VIIRS-IFF','VIIRS-NASA-HRES','VGAC')
         DEFAULT_ACHA_MODE  =  ACHA_Mode_Default_VIIRS
      case ('AQUA-IFF')
          DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Modis
      case ('AVHRR-FUSION')
            if (Avhrr_1_flag == sym%YES) then
              DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Avhrr1_Fusion
            else
              DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Avhrr_Fusion
            endif
      case ('AVHRR-IFF')
         DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Avhrr
      case ('COMS-IMAGER')
         DEFAULT_ACHA_MODE  = ACHA_Mode_Default_COMS
      case ('MODIS')
          DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Modis
      case ('MODIS-MAC')
          DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Modis
      case ('MODIS-CSPP')
          DEFAULT_ACHA_MODE  = ACHA_Mode_Default_Modis
      case ('AHI')
          DEFAULT_ACHA_MODE  = ACHA_Mode_Default_AHI
      case ('MERSI-2')
          DEFAULT_ACHA_MODE = ACHA_Mode_Default_Mersi2
      case ('METIMAGE')
          DEFAULT_ACHA_MODE = ACHA_Mode_Default_Modis
      case default
          call MESG("sensor "//trim(SensorName)//" is not set in user_options.f90: check channels settings Inform andi.walther@ssec.wisc.edu")
      end select

   end function DEFAULT_ACHA_MODE

   !-----------------------------------------------------------------
   !   returns default dcomp mode
   ! -----------------------------------------------------------------
   integer function DEFAULT_DCOMP_MODE ( )

      implicit none

      default_DCOMP_Mode = 3

      if (Sensor%WMO_Id == 208) Default_Dcomp_Mode = 1 !NOAA-17
      if (Sensor%WMO_Id == 3) Default_Dcomp_Mode = 1   !METOP-A
      if (Sensor%WMO_Id == 4) Default_Dcomp_Mode = 1   !METOP-B
      if (Sensor%WMO_Id == 5) Default_Dcomp_Mode = 1   !METOP-C

!---- AKH - What about NOAA-16 pre 2002 or noaa-16 post 2017?

   end function DEFAULT_DCOMP_MODE

!-----------------------------------------------------------------
!   returns default classifier name
!-----------------------------------------------------------------
   function DEFAULT_NB_MASK_CLASSIFIER_FILE (SensorName) result (filename)
      character ( len = *) , intent(in) :: SensorName
      character ( len = 1020 ) :: filename

      select case ( trim(SensorName))
      case ('FCI')
        filename  = 'ahi_default_nb_cloud_mask_lut.nc'
      case ( 'AVHRR-1')
         filename  = 'avhrr_default_nb_cloud_mask_lut.nc'
      case ( 'AVHRR-2')
         filename  = 'avhrr_default_nb_cloud_mask_lut.nc'
      case ( 'AVHRR-3')
         filename  = 'avhrr_default_nb_cloud_mask_lut.nc'
      case ( 'GOES-MP-IMAGER')
         filename  = 'goesmp_default_nb_cloud_mask_lut.nc'
      case ( 'GOES-IL-IMAGER')
         filename  = 'goesil_default_nb_cloud_mask_lut.nc'
      case ( 'GOES-IP-SOUNDER')
         filename  = 'goesmp_default_nb_cloud_mask_lut.nc'
      case ( 'GOES-RU-IMAGER')
         filename  = 'ahi_default_nb_cloud_mask_lut.nc'
      case ( 'MTSAT-IMAGER')
          filename  = 'mtsat_default_nb_cloud_mask_lut.nc'
      case ('SEVIRI')
         filename  = 'seviri_default_nb_cloud_mask_lut.nc'
      case ('FY2-IMAGER')
         filename  = 'fy2_default_nb_cloud_mask_lut.nc'
      case ('AGRI') !FY4A
         filename  = 'ahi_default_nb_cloud_mask_lut.nc'
      case ('VIIRS','VIIRS-NASA','VIIRS-NASA-HRES','VIIRS-IFF','VGAC')
         filename  = 'viirs_default_nb_cloud_mask_lut_fw_10312018.nc'
      case ('AQUA-IFF')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'
      case ('AVHRR-IFF')
        filename  = 'avhrr_default_nb_cloud_mask_lut.nc'
      case ('AVHRR-FUSION')
         filename = 'avhrr_fusion_default_nb_cloud_mask_lut.nc'
      case ('COMS-IMAGER')
         filename  = 'coms_default_nb_cloud_mask_lut.nc'
      case ('MODIS')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'
      case ('MODIS-CSPP')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'
      case ('MODIS-MAC')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'
      case ('AHI')
          filename  = 'ahi_default_nb_cloud_mask_lut.nc'
      case ('MERSI-2')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'
      case ('METIMAGE')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'
      case default
         call MESG("sensor "//SensorName//" is not set in user_options.f90:  Inform andi.walther@ssec.wisc.edu")
         stop
      end select


     end function DEFAULT_NB_MASK_CLASSIFIER_FILE


!-----------------------------------------------------------------
!   returns default classifier name
!-----------------------------------------------------------------
   function DEFAULT_NBM_MASK_CLASSIFIER_FILE (SensorName) result (filename)
      character ( len = *) , intent(in) :: SensorName
      character ( len = 1020 ) :: filename

      select case ( trim(SensorName))

      case ('GOES-RU-IMAGER')
         filename  = 'ecm2_lut_abhi_default.nc'
      case ('MODIS','METIMAGE','MERSI-2')
         filename  = 'ecm2_lut_modis_default.nc'
      case ('VIIRS','VIIRS-NASA','VIIRS-NASA-HRES','VIIRS-IFF','VGAC')
         filename  = 'ecm2_lut_viirs_default.nc'
      case ('AHI')
          filename  = 'ecm2_lut_abhi_default.nc'
      case ('AVHRR-1','AVHRR-2','AVHRR-3')
         filename  = 'ecm2_lut_avhrr_default.nc'
      case ( 'MTSAT-IMAGER')
          filename  = 'ecm2_lut_mtsat2_default.nc'
      case ( 'SEVIRI')
          filename  = 'ecm2_lut_seviri_default.nc'
      case ( 'GOES-MP-IMAGER')
          filename  = 'ecm2_lut_goesmp_src-abhi_default.nc'
      case ( 'GOES-IL-IMAGER')
          filename  = 'ecm2_lut_goesil_src-abhi_default.nc'
      case default
         call MESG("sensor "//TRIM(SensorName)//" does not have ECM2 LUT:  Inform andrew.heidinger@noaa.gov")
         stop
      end select


     end function DEFAULT_NBM_MASK_CLASSIFIER_FILE

   !----------------------------------------------------------------------------
   !  check if algo mode set by user is possible
   !----------------------------------------------------------------------------
   subroutine CHECK_DCOMP_ALGORITHM_CHOICES(SensorName)
      character (len=*) , intent(in) :: SensorName

      integer :: Possible_DCOMP_Modes ( 4 )

      DCOMP_Mode = DCOMP_Mode_User_Set

      Possible_DCOMP_Modes = 0

      select case ( trim ( SensorName))
      case('FCI')
        Possible_DCOMP_Modes(1:4)    = [1,2,3,9]
      case ( 'AVHRR-3')
         Possible_DCOMP_Modes(1:3)    = [1,3,9]
      case ( 'AVHRR-2')
         Possible_DCOMP_Modes(1)    =  3
      case ( 'AVHRR-1')
         Possible_DCOMP_Modes(1)    =  3
      case ( 'GOES-MP-IMAGER')
         Possible_DCOMP_Modes(1)    =  3
      case ( 'GOES-IL-IMAGER')
         Possible_DCOMP_Modes(1)    =  3
      case ( 'GOES-IP-SOUNDER')
         Possible_DCOMP_Modes(1)    =  3
      case ( 'GOES-RU-IMAGER')
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9]
      case ( 'MTSAT-IMAGER')
         Possible_DCOMP_Modes(1)    =  3
      case ('SEVIRI')
         Possible_DCOMP_Modes(1:3)  =  [1, 3, 9]
      case ('FY2-IMAGER')
         Possible_DCOMP_Modes(1:1)  =  [3]
      case ('AGRI')
         Possible_DCOMP_Modes(1)    =  3
      case ('VIIRS','VIIRS-NASA','VIIRS-NASA-HRES')
         NLCOMP_Mode_User_Set       =  1
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9]
      case ('VIIRS-IFF')
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9]
      case ('VGAC')
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9]
      case ('AQUA-IFF')
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9]
      case ('AVHRR-IFF')
         Possible_DCOMP_Modes(1)    =  3
       case ('AVHRR-FUSION')
         Possible_DCOMP_Modes(1)    =  3
      case ('COMS-IMAGER')
         Possible_DCOMP_Modes(1:1)  =  [3]
      case ('MODIS')
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9 ]
      case ('MODIS-CSPP')
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9]
      case ('MODIS-MAC')
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9]
      case ('AHI')
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9]
      case ('MERSI-2')
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9]
      case ('METIMAGE')
         Possible_DCOMP_Modes(1:4)  =  [1, 2, 3, 9]
      case default
         call MESG("sensor "//SensorName//" is not set in check channels user_options settings Inform andi.walther@ssec.wisc.edu")
      end select

      if ( DCOMP_Mode_User_Set /= 0 .and. .not. ANY ( DCOMP_Mode_User_Set == Possible_DCOMP_Modes ) ) then
         DCOMP_Mode = Default_DCOMP_Mode ( )
         print*,EXE_PROMPT,'User set DCOMP mode ',DCOMP_Mode_user_set,' not possible for '// &
                trim(SensorName)//'  with WMO ID ',Sensor%WMO_Id, ' switched to default ', DCOMP_Mode
      end if

   end subroutine CHECK_DCOMP_ALGORITHM_CHOICES

   !----------------------------------------------------------------------------
   !  check if algo mode set by user is possible for ACHA
   !
   !  Allowed Modes
   !  "off" = do not run ACHA
   !  "default" use mode set by DEFAULT_ACHA_MODE()
   !  "maximum" use all relevant channels that are available
   !  "xxx_yyy_zzz" - direct acha mode choice xx < yy < zz
   !               examples, 067_110, 110_120_133, 085_110_120, ...
   !   default and maximum will choose an acha mode of the form above
   !----------------------------------------------------------------------------
   subroutine CHECK_ACHA_ALGORITHM_CHOICES(SensorName)
      character (len=*) , intent(in) :: SensorName

      character (len=ACHA_Mode_Max_Length):: ACHA_Mode_Temp
      integer:: strlen, Mode_Idx

      !--- start with user choice
      ACHA%Mode = ACHA_Mode_User_Set

      !--- if ACHA mode is set to off, there is nothing to do
      if (trim(ACHA%Mode) == "off") then
         return
      endif

      !--- if ACHA mode is set to default, then choose the default mode based on
      !--- sensor.  if a needed channel is missing, it will be handled below
      if (trim(ACHA%Mode) == "default") then
          ACHA%Mode =  DEFAULT_ACHA_MODE ( SensorName )
          call MESG("User selected ACHA default, ACHA Mode = "//trim(ACHA%Mode))
      endif

      !--- if ACHA mode is maximum, that use all channel available
      !--- note, only 375, 67, 85, 11, 12 and 133 are available
      !--- as more channel added to acha, these needs modification

      ACHA_Mode_Temp = ""

      if (trim(ACHA%Mode) == "maximum") then
         if (Sensor%Chan_On_Flag_Default(20) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_038"
         endif
         if (Sensor%Chan_On_Flag_Default(37) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_062"
         endif
         if (Sensor%Chan_On_Flag_Default(27) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_067"
         endif
         if (Sensor%Chan_On_Flag_Default(28) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_073"
         endif
         if (Sensor%Chan_On_Flag_Default(29) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_085"
         endif
         if (Sensor%Chan_On_Flag_Default(30) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_097"
         endif
         if (Sensor%Chan_On_Flag_Default(31) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_110"
         endif
         if (Sensor%Chan_On_Flag_Default(32) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_120"
         endif
         if (Sensor%Chan_On_Flag_Default(33) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_133"
         endif
         if (Sensor%Chan_On_Flag_Default(34) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_136"
         endif
         if (Sensor%Chan_On_Flag_Default(35) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_139"
         endif
         if (Sensor%Chan_On_Flag_Default(36) == 1) then
            Acha_Mode_Temp = trim(Acha_Mode_Temp) // "_142"
         endif
         !-- take off first "_"
         strlen  = len_trim(Acha_Mode_Temp)
         ACHA%Mode = Acha_Mode_Temp(2:strlen)
         print *, "User selected ACHA maximum, ACHA Mode = ",trim(ACHA%Mode)
      endif

      !---  see if ACHA mode calls for an ir channel that is turned off or unavailable
      !--- if this happens, first try ACHA default mode; if not working, then turn acha off
      if ((Sensor%Chan_On_Flag_Default(20) == 0 .and. index(ACHA%Mode,'038') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(37) == 0 .and. index(ACHA%Mode,'062') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(27) == 0 .and. index(ACHA%Mode,'067') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(28) == 0 .and. index(ACHA%Mode,'073') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(29) == 0 .and. index(ACHA%Mode,'085') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(30) == 0 .and. index(ACHA%Mode,'097') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(38) == 0 .and. index(ACHA%Mode,'104') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(31) == 0 .and. index(ACHA%Mode,'110') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(32) == 0 .and. index(ACHA%Mode,'120') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(33) == 0 .and. index(ACHA%Mode,'133') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(34) == 0 .and. index(ACHA%Mode,'136') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(35) == 0 .and. index(ACHA%Mode,'139') > 0) .or. &
          (Sensor%Chan_On_Flag_Default(36) == 0 .and. index(ACHA%Mode,'142') > 0) ) then

          if (ACHA_Mode_User_Set /= DEFAULT_ACHA_MODE(SensorName) ) then
             ACHA%Mode = DEFAULT_ACHA_MODE(SensorName)
             call MESG("user selected ACHA mode inconsistent with spectral information, try default ACHA Mode")
                   if ((Sensor%Chan_On_Flag_Default(20) == 0 .and. index(ACHA%Mode,'038') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(37) == 0 .and. index(ACHA%Mode,'062') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(27) == 0 .and. index(ACHA%Mode,'067') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(28) == 0 .and. index(ACHA%Mode,'073') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(29) == 0 .and. index(ACHA%Mode,'085') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(30) == 0 .and. index(ACHA%Mode,'097') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(38) == 0 .and. index(ACHA%Mode,'104') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(31) == 0 .and. index(ACHA%Mode,'110') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(32) == 0 .and. index(ACHA%Mode,'120') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(33) == 0 .and. index(ACHA%Mode,'133') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(34) == 0 .and. index(ACHA%Mode,'136') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(35) == 0 .and. index(ACHA%Mode,'139') > 0) .or. &
                       (Sensor%Chan_On_Flag_Default(36) == 0 .and. index(ACHA%Mode,'142') > 0) ) then
                          ACHA%Mode = "off"
                          call MESG("default ACHA mode inconsistent with spectral information, ACHA Mode = off")
                   endif
          else
             ACHA%Mode = "off"
             call MESG("user selected ACHA mode inconsistent with spectral information, ACHA Mode = off")
          endif
          return
      endif

      !--------------------------------------------------------------------------------
      ! check if mode is one that is allowed, if not, stop
      !--------------------------------------------------------------------------------
      do Mode_Idx = 1, Num_Acha_Modes
         if (trim(ACHA%Mode) == trim(ACHA_Mode_Values(Mode_Idx))) then
            exit
         endif
      enddo

      if (ACHA%Mode == 'baseline') then
         call MESG( "Running AWG Baseline Cloud Height ",level = verb_lev % DEFAULT)
         Mode_Idx = -1
      endif

      if (Mode_Idx == Num_Acha_Modes + 1) then
         call MESG("Error, Unsupported ACHA Mode = "//trim(ACHA%Mode))
         print *, "CLAVR-x:  Supported ACHA Modes = ",ACHA_Mode_Values
         stop
      endif

   end subroutine CHECK_ACHA_ALGORITHM_CHOICES

   ! ----------------------------------------------------------------------
   !    returns all available sensors for this sensors
   ! ----------------------------------------------------------------------
   function EXISTING_CHANNELS (SensorName)  result( Valid_Channels )
      character (len = *) , intent(in) :: SensorName

      integer , target :: Valid_Channels (Nchan_Clavrx)
      integer :: i

      Valid_Channels = MISSING_VALUE_INT4

      select case ( trim(SensorName))
      case ('FCI')
        Valid_Channels(1:16) = [1,6,7,9,11,16,19,20,26,28,29,30,31,32,33,37]
      case ( 'AVHRR-1')
         Valid_Channels (1:4) = [1,2,20,31]
      case ( 'AVHRR-2')
         Valid_Channels (1:5) = [1,2,20,31,32]
      case ( 'AVHRR-3')
         Valid_Channels (1:6) = [1,2,6,20,31,32]
      case ( 'GOES-IL-IMAGER')
         Valid_Channels (1:5) = [1,20,27,31,32]
      case ( 'GOES-MP-IMAGER')
         Valid_Channels (1:5) = [1,20,27,31,33]
      case ( 'GOES-IP-SOUNDER')
         Valid_Channels (1:18) = [1,20,21,23,24,25,30,31,32,33,34,35,36,37,38,39,40,41]
      case ('GOES-RU-IMAGER')
         Valid_Channels(1:16) = [1,2,3,6,7,20,26,27,28,29,30,31,32,33,37,38]
      case ( 'MTSAT-IMAGER')
         Valid_Channels (1:5) = [1,20,27,31,32]
      case ('SEVIRI')
         Valid_Channels (1:11) = [1,2,6,20,37,28,29,30,31,32,33]
      case ('FY2-IMAGER')
         Valid_Channels (1:5) = [1,20,27,31,32]
      case ('AGRI') ! FY4A
         Valid_Channels (1:14) = [1,2,3,6,7,20,21,26,27,28,29,31,32,33]
      case ('VIIRS')
         Valid_Channels (1:22) = [1,2,3,4,5,6,7,8,9,15,20,22,26,29,31,32,39,40,41,42,43,44]
      case ('VGAC')
         if (.not. Sensor%Fusion_Flag) then
           Valid_Channels (1:16) = [1,2,3,4,5,6,7,8,9,15,20,22,26,29,31,32]
         else
           Valid_Channels (1:25) = [1,2,3,4,5,6,7,8,9,15,20,22,24,25,26,27,28,29,30,31,32,33,34,35,36]
         endif
      case ('VIIRS-NASA','VIIRS-NASA-HRES')
         if (.not. Sensor%Fusion_Flag) then
           Valid_Channels (1:22) = [1,2,3,4,5,6,7,8,9,15,20,22,26,29,31,32,39,40,41,42,43,44]
         else
           Valid_Channels (1:32) = [1,2,3,4,5,6,7,8,9,15,20,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,39,40,41,42,43,44]
         endif
      case ('VIIRS-IFF')
         Valid_Channels (1:25) = [1,2,3,4,5,6,7,8,9,15,20,22,26,27,28,29,30,31,32,33,34,35,36,37,38]
      case ('AQUA-IFF')
         Valid_Channels (1:36) = [(i,i=1,36,1)]
      case ('AVHRR-IFF')
         Valid_Channels (1:21) = [1,2,6,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36,37,38]
      case ('AVHRR-FUSION')
         Valid_Channels (1:21) = [1,2,6,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36,37,38]
      case ('COMS-IMAGER')
         Valid_Channels (1:5) = [1,20,27,31,32]
      case ('MODIS')
         Valid_Channels(1:36) = [(i,i=1,36,1)]
      case ('MODIS-MAC')
         Valid_Channels(1:36) = [(i,i=1,36,1)]
      case ('MODIS-CSPP')
         Valid_Channels(1:36) = [(i,i=1,36,1)]
      case ('AHI')
         Valid_Channels(1:16) = [1,2,3,4,6,7,20,27,28,29,30,31,32,33,37,38]
      case ('MERSI-2')
         Valid_Channels(1:24) = [3,4,1,2,26,6,7,8,9,10,12,13,15,16,17,18,19,5,20,23,28,29,31,32]
         !Valid_Channels(1:21) = [3,4,1,2,26,6,7,8,9,10,15,17,18,19,5,20,23,28,29,31,32]
      case ('METIMAGE')
         Valid_Channels(1:20) = [1,2,3,4,5,6,7,15,18,20,22,23,26,27,28,29,31,32,33,45]
      case default
         call MESG("sensor "//SensorName//" is not set in check channels settings Inform andi.walther@ssec.wisc.edu")
      end select


   end function EXISTING_CHANNELS

   !----------------------------------------------------------------------
   !   Channel settings
   !     will not be done for full-experts  ( expert mode 7 and higher)
   !     this turns on all channels availble for this sensor
   !----------------------------------------------------------------------
   subroutine CHANNEL_SWITCH_ON (SensorName)
      character (len=*) , intent(in) :: SensorName
      integer :: Valid_Channels (Nchan_Clavrx)
      integer :: i

      ! expert can decide themselves
      if (Expert_Mode > 6 ) return

      Valid_Channels = Existing_Channels ( SensorName )

      Sensor%Chan_On_Flag_Default =  0_int1

      do i = 1, Nchan_Clavrx
         if (Valid_Channels (i) < 0 ) cycle
         Sensor%Chan_On_Flag_Default (Valid_Channels(i)) = 1_int1
      end do

   end subroutine CHANNEL_SWITCH_ON

   !----------------------------------------------------------------------
   ! Determine if sub-pixel variables are needed
   ! currently this is only possible for fixed-grid geo using the static nav
   ! files
   ! this is only setup for the ABI and AHI reflectance channels that are
   ! higher resolution than the thermal bands.
   ! Also, I-Bands on VIIRS
   !----------------------------------------------------------------------
   subroutine SUB_PIXEL_CHANNEL_ON_SET()

      Ch(1:NCHAN_CLAVRX)%Sub_Pixel_On_Flag = .false.

      !--- ISCCP-NG L1g
      if (index(Image%Level1b_Name,'ISCCP-NG_L1g') > 0) then
         if (Sensor%Chan_On_Flag_Default(1) == sym%YES) Ch(1)%Sub_Pixel_On_Flag = .true.  !0.65
         if (Sensor%Chan_On_Flag_Default(31) == sym%YES) Ch(31)%Sub_Pixel_On_Flag = .true.  !11
         return
      endif

      !--- also check that is not an area file
      if (Image%Mixed_Resolution_Flag) then
         if (Image%Chan_Average_Flag == 2) then
            if (Sensor%Chan_On_Flag_Default(1) == sym%YES) Ch(1)%Sub_Pixel_On_Flag = .true.  !0.65
            if (Sensor%Chan_On_Flag_Default(2) == sym%YES) Ch(2)%Sub_Pixel_On_Flag = .true.  !0.86
            if (Sensor%Chan_On_Flag_Default(3) == sym%YES) Ch(3)%Sub_Pixel_On_Flag = .true.  !0.47
            if (Sensor%Chan_On_Flag_Default(4) == sym%YES) Ch(4)%Sub_Pixel_On_Flag = .true.  !0.51
            if (Sensor%Chan_On_Flag_Default(6) == sym%YES) Ch(6)%Sub_Pixel_On_Flag = .true.  !1.60
         endif
      endif

      !--- Treat VIIRS I-Bands as higher resolution versions of analogous M-Bands
      if ((Sensor%WMO_Id == 224 .or. Sensor%WMO_Id == 225 .or. Sensor%WMO_Id == 226)) then

        !--- VGAC has no I-bands but it does report subpixel stddev of vis and irwin
        if (Sensor%Sensor_Name == "VGAC") then

            if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
              Ch(1)%Sub_Pixel_On_Flag = .true.  !0.65
            endif
            if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
              Ch(31)%Sub_Pixel_On_Flag = .true.  !11.4
            endif

        else   !--- normal viirs with I-bands

            !-- channel 1 is 0.65 M-band and channel 39 is 0.64 I-band
            if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. &
                Sensor%Chan_On_Flag_Default(39) == sym%YES) then
                Ch(1)%Sub_Pixel_On_Flag = .true.  !0.65
            endif

            !-- channel 2 is 0.86 M-band and channel 40 is 0.86 I-band
            if (Sensor%Chan_On_Flag_Default(2) == sym%YES .and. &
                Sensor%Chan_On_Flag_Default(40) == sym%YES) then
                Ch(2)%Sub_Pixel_On_Flag = .true.  !0.86
            endif

            !-- channel 6 is 1.6 M-band and channel 41 is 1.6 I-band
            if (Sensor%Chan_On_Flag_Default(6) == sym%YES .and. &
                Sensor%Chan_On_Flag_Default(41) == sym%YES) then
                Ch(6)%Sub_Pixel_On_Flag = .true.  !1.6
            endif

            !-- channel 20 is 3.75 M-band and channel 42 is 3.75 I-band
            if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. &
                Sensor%Chan_On_Flag_Default(42) == sym%YES) then
                Ch(20)%Sub_Pixel_On_Flag = .true.  !3.75
            endif

            !-- channel 31 is 11.0 M-band and channel 43 is 11.4 I-band
            if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
                Sensor%Chan_On_Flag_Default(43) == sym%YES) then
                Ch(31)%Sub_Pixel_On_Flag = .true.  !11.4
            endif

        endif

      endif

   end subroutine SUB_PIXEL_CHANNEL_ON_SET

   !----------------------------------------------------------------------
   ! turn off channels not support by sensor (even if selected by user)
   !----------------------------------------------------------------------
   subroutine CHECK_USER_CHANNEL_CHOICES (SensorName)
      character (len=*) , intent(in) :: SensorName
      integer :: Valid_Channels (Nchan_Clavrx)
      integer :: i

      Valid_Channels = EXISTING_CHANNELS(SensorName)

      do i = 1, Nchan_Clavrx
         if ( any ( i == Valid_Channels )) cycle
         Sensor%Chan_On_Flag_Default(i) = 0_int1
      end do

   end subroutine CHECK_USER_CHANNEL_CHOICES
   ! --------------------------------------------------------------------
   !  every incosistency between channel settings and algorithm mode
   ! --------------------------------------------------------------------
   subroutine  EXPERT_MODE_CHANNEL_ALGORITHM_CHECK ( )

      logical :: Not_Run_Flag

      if ( Expert_Mode < 6 ) return

      !--- check ACHA mode based on available channels
      Not_Run_Flag = .false.
      if (ACHA%Mode == "off") Not_Run_Flag = .true.

      if ( Not_Run_Flag ) then
         print *, EXE_PROMPT, 'ACHA Mode ', ACHA%Mode,' not possible with selected channels. ACHA and DCOMP  will not run.'
         Dcomp_Mode = 0
      end if

      !--- check based on available channels
      if (Dcomp_Mode == 1 .and. &
         (Sensor%Chan_On_Flag_Default(1) == sym%NO .or. Sensor%Chan_On_Flag_Default(6)==sym%NO)) then
         call MESG("DCOMP Mode 1 not possible with selected channels, DCOMP is now off")
      endif

      if (Dcomp_Mode == 2 .and. &
         (Sensor%Chan_On_Flag_Default(1) == sym%NO .or. Sensor%Chan_On_Flag_Default(7)==sym%NO)) then
         call MESG("DCOMP Mode 2 not possible with selected channels, DCOMP is now off")
      endif

      if (Dcomp_Mode == 3 .and. &
         (Sensor%Chan_On_Flag_Default(1) == sym%NO .or. Sensor%Chan_On_Flag_Default(20)==sym%NO)) then
         call MESG("DCOMP Mode 3 not possible with selected channels, DCOMP is now off")
      endif

   end subroutine EXPERT_MODE_CHANNEL_ALGORITHM_CHECK

end module USER_OPTIONS
