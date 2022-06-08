! $Id: level2_mod.f90 4128 2021-04-19 02:03:04Z heidinger $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: level2_mod.f90 (src)
!       LEVEL2_MOD (program)
!
! PURPOSE: Routines for creating, writing and closing pixel-level output files
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
! REVISION HISTORY:
!          5/09 - Created
!--------------------------------------------------------------------------------------
module LEVEL2_MOD
 use CONSTANTS_MOD
 use FILE_TOOLS,only: FILE_NR_LINES, GETLUN, FILE_TEST
 use PIXEL_COMMON_MOD
 use HDF
 use CX_STRING_TOOLS_MOD
 use CX_REAL_BOOLEAN_MOD
 use CX_HDF4_MOD
 use NETCDF
 use CX_NETCDF4_MOD
 use CLAVRX_MESSAGE_MOD
 use CX_MURI_CLAVRX_BRIDGE_MOD, only: muri
 use LEVEL2_STRUCTURES_MOD, only: Sds_Struct, Clavrx_Global_Attr, SET_CLAVRX_GLOBAL_ATTRIBUTES

 implicit none

 private

 public:: DEFINE_HDF_FILE_STRUCTURES, &
          WRITE_SEGMENT_LEVEL2, &
          CLOSE_PIXEL_HDF_FILES, &
          WRITE_ALGORITHM_ATTRIBUTES, &
          READ_LEVEL2_VAR_LIST, &
          SETUP_LEVEL2_SDS_INFO

 private:: WRITE_SDS_HDF, &
           DEFINE_LEVEL2_SDS_HDF, &
           WRITE_FLAG_VALUES_HDF, &
           WRITE_2D_REAL4_SCALED_SDS_HDF, &
           SET_CH_ATTRIBUTES, &
           WRITE_SDS_NETCDF, &
           DEFINE_LEVEL2_SDS_NETCDF
          
 !--- Level2 File Index
 integer, private, save:: Sd_Id_Level2

 !--- number of level2 variables in the level2_list-
 integer, private, save:: Num_Level2_Vars_List

 !--- Dimension Id
 !integer(kind=int4),private:: Dim_Id

 !--- compression variables
 integer(kind=int4), dimension(2),private, save::  Comp_Prm
 integer(kind=int4), private, save::  Comp_Type

 integer(kind=int4),parameter,private:: Sds_Rank_1d = 1
 integer(kind=int4),dimension(1),private:: Sds_Dims_1d

 integer(kind=int4),parameter,private:: Sds_Rank_2d = 2
 integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Dims_2d
 integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Start_2d
 integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Stride_2d
 integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Edge_2d
 integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Chunk_Size_2d

 integer(kind=int4),parameter,private:: Sds_Rank_3d = 3
 integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Dims_3d
 integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Start_3d
 integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Stride_3d
 integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Edge_3d
 integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Chunk_Size_3d

 character(len=11), private, parameter:: MOD_PROMPT = "LEVEL2:"
 character(len=18), private, parameter:: coordinates_string = "longitude latitude"

 type(Sds_Struct), private, save, dimension(:), allocatable:: Sds_Info 

!====================================================================
!
!====================================================================

 contains
!====================================================================
! SUBROUTINE Name: READ_LEVEL2_VAR_LIST
!
! Function:
!   Reads a text file of what is wanted in the level2
!
! Description:
!====================================================================
subroutine READ_LEVEL2_VAR_LIST()
   character(len=200):: Level2_Var_List_Header
   character(len=200) :: String_Dummy
   !integer:: Lun, Sds_Counter, Var_Idx, Num_Occurs
   integer:: Lun, Var_Idx, Num_Occurs
  
   if ( .not. file_test(trim(Level2_List))) then
      stop 'Please add a file level2_list'
   end if
   
   Num_Level2_Vars_List = file_nr_lines (trim(Level2_List)) - 1
   
   if ( Num_Level2_Vars_List < 1) then
      stop 'No variable names in level2_list'
   end if
   
   Lun = GETLUN()

   open(unit=Lun,file=trim(Level2_List),status='old')
   read(unit=lun,fmt=*) Level2_Var_List_Header
   read(unit=Lun,fmt=*) String_Dummy
   
   if (is_numeric(trim(string_dummy))) then
      print*
      print*,'        You still use the number of level2 products in your level2_list line 2'
      print*,'        This is not needed anymore. Please delete line 2 with the number from level2_list '
      print*,'        Tolerant exception  handling for this processing expired on Dec 31 2018!!'
      print*,'        We stop now... '
      print*
      stop ' PLEASE REMOVE THE LINE WITH THE NUMBER OF LEVEL2 PRODUCTS  in your level2_list first.. '
      Num_Level2_Vars_List = Num_Level2_Vars_List - 1
   else
      backspace(lun)
   end if
   
  
   allocate(Sds_Info(Num_Level2_Vars_List))

   do Var_Idx = 1, Num_Level2_Vars_List
      read(unit=Lun,fmt=*) Sds_Info(Var_Idx)%Sds_Name
   enddo

   !------------------------------------------------------------------
   ! Loop through variable names and look for and remove duplicates
   ! sds name will be changed to "duplicate" and should be ignored
   !------------------------------------------------------------------
   do Var_Idx = 2, Num_Level2_Vars_List
       Num_Occurs = count(Sds_Info(Var_Idx)%Sds_Name == Sds_Info(1:Var_Idx)%Sds_Name)
       if (Num_Occurs > 1) then 
          call MESG(MOD_PROMPT//"Duplicate level2 sds name in level2_list ==> "//trim(Sds_Info(Var_Idx)%Sds_Name))
          Sds_Info(Var_Idx)%Sds_Name = "duplicate"    !rename to duplicate
       endif
   enddo

   close(unit=Lun)

end subroutine READ_LEVEL2_VAR_LIST
!------------------------------------------------------------------------------
! Main routine to populate sds_info structure which defines each sds
!------------------------------------------------------------------------------
subroutine SETUP_LEVEL2_SDS_INFO()
   integer:: Var_Idx
   do Var_Idx = 1, Num_Level2_Vars_List

      !--- set default strcuture for a real4 2d scaled sds
      Sds_Info(Var_Idx)%On_Flag =  .true.
      Sds_Info(Var_Idx)%Sds_Idx = Var_Idx 
      Sds_Info(Var_Idx)%Rank = 2
      Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_FLOAT32
      Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT16
      Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_SHORT
      Sds_Info(Var_Idx)%Scaling_Type = 1_int1
      Sds_Info(Var_Idx)%Units =  "none"
      Sds_Info(Var_Idx)%Standard_Name = "none"
      Sds_Info(Var_Idx)%Long_Name = "none"
      Sds_Info(Var_Idx)%Flags_String = "none"
      Sds_Info(Var_Idx)%Number_Of_Flags = 0
      Sds_Info(Var_Idx)%Valid_Range = real([TWO_BYTE_MIN,TWO_BYTE_MAX])
      Sds_Info(Var_Idx)%Scale_Factor = 1.0
      Sds_Info(Var_Idx)%Add_Offset = 1.0
      Sds_Info(Var_Idx)%Actual_Range = [0,0]

      Sds_Info(Var_Idx)%Sds_Data_1d_I1 => null()
      Sds_Info(Var_Idx)%Sds_Data_1d_I4 => null()
      Sds_Info(Var_Idx)%Sds_Data_1d_R4 => null()
      Sds_Info(Var_Idx)%Sds_Data_2d_I1 => null()
      Sds_Info(Var_Idx)%Sds_Data_2d_I2 => null()
      Sds_Info(Var_Idx)%Sds_Data_2d_I4 => null()
      Sds_Info(Var_Idx)%Sds_Data_2d_R4 => null()
      Sds_Info(Var_Idx)%Sds_Data_3d_I1 => null()

      select case (trim(Sds_Info(Var_Idx)%Sds_Name))
        
         !----------------------------------------------------------------------------------------------------
         !  
         !----------------------------------------------------------------------------------------------------
         case("scan_line_number")
            Sds_Info(Var_Idx)%Standard_Name = "scan_line_number"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Rank = 1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT32
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT32
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF =  NF90_INT
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Image%Scan_Number)) Sds_Info(Var_Idx)%Sds_Data_1d_I4 => Image%Scan_Number
         case("scan_line_time")
            Sds_Info(Var_Idx)%Standard_Name = "scan_line_time"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Rank = 1
            Sds_Info(Var_Idx)%Long_Name="time for the scan line in fractional hours"
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_FLOAT32
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_FLOAT
            Sds_Info(Var_Idx)%Units = "hours"
            if (allocated(Image%Utc_Scan_Time_Hours)) Sds_Info(Var_Idx)%Sds_Data_1d_R4 => Image%Utc_Scan_Time_Hours
         case("bad_scan_line_flag")
            Sds_Info(Var_Idx)%Standard_Name = "bad_scan_line_flag"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Rank = 1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Bad_Scan_Flag)) Sds_Info(Var_Idx)%Sds_Data_1d_I1 => Bad_Scan_Flag
         !----------------------------------------------------------------------------------------------------
         !  packed arrays
         !----------------------------------------------------------------------------------------------------
         case("packed_land_cover")
            Sds_Info(Var_Idx)%Standard_Name = "packed_land_cover"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
         case("packed_pixel_meta_data")
            Sds_Info(Var_Idx)%Standard_Name = "pixel_quality_flags_packed_into_one_byte"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
         !----------------------------------------------------------------------------------------------------
         !  2d Int8 fields
         !----------------------------------------------------------------------------------------------------
         case("bad_pixel_mask")
            Sds_Info(Var_Idx)%Standard_Name = "bad_pixel_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Bad_Pixel_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Bad_Pixel_Mask
         case("gap_pixel_mask")
            Sds_Info(Var_Idx)%Standard_Name = "gap_pixel_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Gap_Pixel_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Gap_Pixel_Mask
         case("dust_mask")
            Sds_Info(Var_Idx)%Standard_Name = "binary_dust_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Long_Name = "dust mask (0=no,1=yes)"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CLDMASK%Dust_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Dust_Mask
         case("dust_probability")
            Sds_Info(Var_Idx)%Standard_Name = "dust_probability"
            Sds_Info(Var_Idx)%Long_Name = "probability of dust"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CLDMASK%Dust_Prob)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CLDMASK%Dust_Prob
         case("smoke_mask")
            Sds_Info(Var_Idx)%Standard_Name = "smoke_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CLDMASK%Smoke_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Smoke_Mask
         case("fire_mask")
            Sds_Info(Var_Idx)%Standard_Name = "fire_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CLDMASK%Fire_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Fire_Mask
         case("shadow_mask")
            Sds_Info(Var_Idx)%Standard_Name = "shadow_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CLDMASK%Shadow_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Shadow_Mask
         !----------------------------------------------------------------------------------------------------
         !  Unscaled 2D Reals
         !----------------------------------------------------------------------------------------------------
        
         case("radiance_dnb_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_radiance_dnb_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_FLOAT32
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_FLOAT
            Sds_Info(Var_Idx)%Units = "mW/m^2/cm^-1"
            if (allocated(ch(44)%Rad_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(44)%Rad_Toa   
         
         case("diagnostic_1")
            Sds_Info(Var_Idx)%Standard_Name = "diagnostic_1"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_FLOAT32
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_FLOAT
            Sds_Info(Var_Idx)%Units = "unknown"
            if (allocated(Diag_Pix_Array_1)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Diag_Pix_Array_1
         case("diagnostic_2")
            Sds_Info(Var_Idx)%Standard_Name = "diagnostic_2"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_FLOAT32
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_FLOAT
            Sds_Info(Var_Idx)%Units = "unknown"
            if (allocated(Diag_Pix_Array_2)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Diag_Pix_Array_2
         case("diagnostic_3")
            Sds_Info(Var_Idx)%Standard_Name = "diagnostic_3"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_FLOAT32
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_FLOAT
            Sds_Info(Var_Idx)%Units = "unknown"
            if (allocated(Diag_Pix_Array_3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Diag_Pix_Array_3
         !----------------------------------------------------------------------------------------------------
         ! Sfc Members
         !----------------------------------------------------------------------------------------------------
         case("glint_mask")
            Sds_Info(Var_Idx)%Standard_Name = "glint_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%Glint_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%Glint_Mask
         case("glint_mask_lunar")
            Sds_Info(Var_Idx)%Standard_Name = "glint_mask_lunar"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%Glint_Mask_Lunar)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%Glint_Mask_Lunar
         case("glint_zenith_angle")
            Sds_Info(Var_Idx)%Standard_Name = "glint_zenith_angle"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,90.0]
            Sds_Info(Var_Idx)%Units = "degrees"
            if (allocated(Geo%Glintzen)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Glintzen
         case("lunar_glint_zenith_angle")
            Sds_Info(Var_Idx)%Standard_Name = "lunar_glint_zenith_angle"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,90.0]
            Sds_Info(Var_Idx)%Units = "degrees"
            if (allocated(Geo%Glintzen_Lunar)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Glintzen_Lunar
         case("moon_illum_frac")
            Sds_Info(Var_Idx)%Standard_Name = "moon_illum_frac"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100.0]
            if (allocated(Geo%LunFrac)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%LunFrac
         case("scattering_angle")
            Sds_Info(Var_Idx)%Standard_Name = "scattering_angle"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,180.0]
            Sds_Info(Var_Idx)%Units = "degrees"
            if (allocated(Geo%Scatangle)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Scatangle
         case("lunar_scattering_angle")
            Sds_Info(Var_Idx)%Standard_Name = "lunar_scattering_angle"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,180.0]
            Sds_Info(Var_Idx)%Units = "degrees"
            if (allocated(Geo%Scatangle_Lunar)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Scatangle_Lunar
         case("city_mask")
            Sds_Info(Var_Idx)%Standard_Name = "city_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%City_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%City_Mask
         case("forward_scatter_mask")
            Sds_Info(Var_Idx)%Standard_Name = "forward_scatter_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%Forward_Scatter_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%Forward_Scatter_Mask
         case("forward_scatter_mask_lunar")
            Sds_Info(Var_Idx)%Standard_Name = "forward_scatter_mask_lunar"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%Forward_Scatter_Mask_Lunar)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%Forward_Scatter_Mask_Lunar
         case("surface_type")
            Sds_Info(Var_Idx)%Standard_Name = "area_type"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Number_Of_Flags = 14
            Sds_Info(Var_Idx)%Flags_String = "water "// &
                               " evergreen_needle "// &
                               " evergreen_broad "// &
                               " deciduous_needle "// &
                               " deciduous_broad "// &
                               " mixed_forest "// &
                               " woodlands "// &
                               " wooded_grass "// &
                               " closed_shrubs "// &
                               " open_shrubs "// &
                               " grasses "// &
                               " croplands "// &
                               " bare "// &
                               " urban "
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%Sfc_Type)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%Sfc_Type
         case("land_class")
            Sds_Info(Var_Idx)%Standard_Name = "land_class"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Number_Of_Flags = 8
            Sds_Info(Var_Idx)%Flags_String = "ocean "// &
                        " land "// &
                        " coastline "// &
                        " shallow_inland_water "// &
                        " ephemeral_water "// &
                        " deep_inland_water "// &
                        " moderate_ocean "// &
                        " deep_ocean "
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%Land)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%Land
         case("snow_class")
            Sds_Info(Var_Idx)%Standard_Name = "snow_class"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Number_Of_Flags = 3
            Sds_Info(Var_Idx)%Flags_String =  "no_snow_or_ice sea_ice snow "
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%Snow)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%Snow
         case("snow_class_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "snow_class_nwp"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Number_Of_Flags = 3
            Sds_Info(Var_Idx)%Flags_String =  "no_snow_or_ice sea_ice snow "
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%Snow_Nwp)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%Snow_Nwp
         case("coast_mask")
            Sds_Info(Var_Idx)%Standard_Name = "coast_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%Coast_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%Coast_Mask
         case("land_mask")
            Sds_Info(Var_Idx)%Standard_Name = "land_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Sfc%Land_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Sfc%Land_Mask
         case("surface_elevation")
            Sds_Info(Var_Idx)%Standard_Name = "surface_elevation"
            Sds_Info(Var_Idx)%Actual_Range = [-500.0,10000.0]
            Sds_Info(Var_Idx)%Units =  "meters"
            if (allocated(Sfc%Zsfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Sfc%Zsfc
         case("surface_elevation_max_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "surface_elevation_max_3x3"
            Sds_Info(Var_Idx)%Actual_Range = [-500.0,10000.0]
            Sds_Info(Var_Idx)%Units =  "meters"
            if (allocated(Sfc%Zsfc_Max)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Sfc%Zsfc_Max
         case("surface_elevation_stddev_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "surface_elevation_stddev_3x3"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,5000.0]
            Sds_Info(Var_Idx)%Units =  "meters"
            if (allocated(Sfc%Zsfc_Std)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Sfc%Zsfc_Std
         case("layer_l1g")
            Sds_Info(Var_Idx)%Standard_Name = "layer_l1g"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(L1g%Layer_Idx)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => L1g%Layer_Idx
         case("sample_mode_l1g")
            Sds_Info(Var_Idx)%Standard_Name = "sample_mode_l1g"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(L1g%Sample_Mode)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => L1g%Sample_Mode
         case("wmo_id_l1g")
            Sds_Info(Var_Idx)%Standard_Name = "wmo_id_l1g"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT16
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT16
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_SHORT
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(L1g%WMO_Id)) Sds_Info(Var_Idx)%Sds_Data_2d_I2 => L1g%WMO_Id

         !----------------------------------------------------------------------------------------------------
         ! Nav Members
         !----------------------------------------------------------------------------------------------------
         case("asc_des_flag") 
            Sds_Info(Var_Idx)%Standard_Name = "asc_des_flag"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Rank = 1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(Nav%Ascend)) Sds_Info(Var_Idx)%Sds_Data_1d_I1 => Nav%Ascend
         case("x")
            Sds_Info(Var_Idx)%Standard_Name = "x"
            Sds_Info(Var_Idx)%Long_Name = "x index of pixel"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Rank = 2
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT32
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT32
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_INT
            if (allocated(Nav%X)) Sds_Info(Var_Idx)%Sds_Data_2d_I4 => Nav%X
         case("y")
            Sds_Info(Var_Idx)%Standard_Name = "y"
            Sds_Info(Var_Idx)%Long_Name = "y index of pixel"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Rank = 2
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT32
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT32
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_INT
            if (allocated(Nav%Y)) Sds_Info(Var_Idx)%Sds_Data_2d_I4 => Nav%Y
         case("latitude")
            Sds_Info(Var_Idx)%Standard_Name = "latitude"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_FLOAT32
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_FLOAT
            Sds_Info(Var_Idx)%Actual_Range = [-90.0,90.0]
            Sds_Info(Var_Idx)%Valid_Range = [-90.0,90.0]
            Sds_Info(Var_Idx)%Units =  "degrees_north"
            if (allocated(Nav%Lat)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Nav%Lat
         case("longitude")
            Sds_Info(Var_Idx)%Standard_Name = "longitude"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_FLOAT32
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_FLOAT
            Sds_Info(Var_Idx)%Actual_Range = [-180.0,180.0]
            Sds_Info(Var_Idx)%Valid_Range = [-180.0,180.0]
            Sds_Info(Var_Idx)%Units =  "degrees_east"
            if (allocated(Nav%Lon)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Nav%Lon
         case("latitude_pc")
            Sds_Info(Var_Idx)%Standard_Name = "latitude_parallax_corrected"
            Sds_Info(Var_Idx)%Actual_Range = [-90.0,90.0]
            Sds_Info(Var_Idx)%Units =  "degrees_north"
            Sds_Info(Var_Idx)%Long_Name = "latitude_parallax_corrected_using_cloud_height"
            if (allocated(Nav%Lat_Pc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Nav%Lat_Pc
         case("longitude_pc")
            Sds_Info(Var_Idx)%Standard_Name = "longitude_parallax_corrected"
            Sds_Info(Var_Idx)%Actual_Range = [-180.0,180.0]
            Sds_Info(Var_Idx)%Units =  "degrees_east"
            Sds_Info(Var_Idx)%Long_Name = "longitude_parallax_corrected_using_cloud_height"
            if (allocated(Nav%Lon_Pc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Nav%Lon_Pc
         !----------------------------------------------------------------------------------------------------
         ! Geo Members
         !----------------------------------------------------------------------------------------------------
         case("sensor_zenith_angle")
            Sds_Info(Var_Idx)%Standard_Name = "sensor_zenith_angle"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,90.0]
            Sds_Info(Var_Idx)%Units =  "degrees"
            if (allocated(Geo%Satzen)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Satzen
         case("sensor_azimuth_angle")
            Sds_Info(Var_Idx)%Standard_Name = "sensor_azimuth_angle"
            Sds_Info(Var_Idx)%Actual_Range = [-180.0,180.0]
            Sds_Info(Var_Idx)%Units =  "degrees"
            if (allocated(Geo%Sataz)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Sataz
         case("solar_zenith_angle")
            Sds_Info(Var_Idx)%Standard_Name = "solar_zenith_angle"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,180.0]
            Sds_Info(Var_Idx)%Units =  "degrees"
            if (allocated(Geo%Solzen)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Solzen
         case("relative_azimuth_angle")
            Sds_Info(Var_Idx)%Standard_Name = "relative_sensor_azimuth_angle"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,180.0]
            Sds_Info(Var_Idx)%Units =  "degrees"
            if (allocated(Geo%Relaz)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Relaz
         case("solar_azimuth_angle") 
            Sds_Info(Var_Idx)%Standard_Name = "solar_azimuth_angle"
            Sds_Info(Var_Idx)%Actual_Range = [-180.0,180.0]
            Sds_Info(Var_Idx)%Units =  "degrees"
            if (allocated(Geo%Solaz)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Solaz
         case("lunar_zenith_angle")
            Sds_Info(Var_Idx)%Standard_Name = "lunar_zenith_angle"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,90.0]
            Sds_Info(Var_Idx)%Units =  "degrees"
            if (allocated(Geo%Lunzen)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Lunzen
         case("lunar_relative_azimuth_angle")
            Sds_Info(Var_Idx)%Standard_Name = "lunar_relative_azimuth_angle"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,180.0]
            Sds_Info(Var_Idx)%Units =  "degrees"
            if (allocated(Geo%LunRelaz)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%LunRelaz
         case("lunar_azimuth_angle")
            Sds_Info(Var_Idx)%Standard_Name = "lunar_azimuth_angle"
            Sds_Info(Var_Idx)%Actual_Range = [-180.0,180.0]
            Sds_Info(Var_Idx)%Units =  "degrees"
            if (allocated(Geo%Lunaz)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Geo%Lunaz

         !----------------------------------------------------------------------------------------------------
         ! Ch Members
         !----------------------------------------------------------------------------------------------------
         case("refl_0_65um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_65um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(1)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(1)%Ref_Toa
         case("refl_0_65um_nom_min_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(1)%Ref_Toa_Min_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(1)%Ref_Toa_Min_Sub
         case("refl_0_65um_nom_max_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(1)%Ref_Toa_Max_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(1)%Ref_Toa_Max_Sub
         case("refl_0_65um_nom_stddev_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(1)%Ref_Toa_Std_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(1)%Ref_Toa_Std_Sub
         case("refl_0_86um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_86um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(2)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(2)%Ref_Toa
         case("refl_0_86um_nom_min_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(2)%Ref_Toa_Min_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(2)%Ref_Toa_Min_Sub
         case("refl_0_86um_nom_max_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(2)%Ref_Toa_Max_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(2)%Ref_Toa_Max_Sub
         case("refl_0_47um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_47um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(3)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(3)%Ref_Toa
         case("refl_0_47um_nom_min_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(3)%Ref_Toa_Min_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(3)%Ref_Toa_Min_Sub
         case("refl_0_47um_nom_max_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(3)%Ref_Toa_Max_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(3)%Ref_Toa_Max_Sub
         case("refl_0_55um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_55um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(4)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(4)%Ref_Toa
         case("refl_0_55um_nom_min_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(4)%Ref_Toa_Min_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(4)%Ref_Toa_Min_Sub
         case("refl_0_55um_nom_max_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(4)%Ref_Toa_Max_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(4)%Ref_Toa_Max_Sub
         case("refl_1_2um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_1_2um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(5)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(5)%Ref_Toa
         case("refl_1_60um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_1_60um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(6)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(6)%Ref_Toa
         case("refl_1_60um_nom_min_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(6)%Ref_Toa_Min_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(6)%Ref_Toa_Min_Sub
         case("refl_1_60um_nom_max_sub")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(6)%Ref_Toa_Max_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(6)%Ref_Toa_Max_Sub
         case("refl_2_10um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_2_10um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(7)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(7)%Ref_Toa
         case("refl_0_41um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_41um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(8)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(8)%Ref_Toa
         case("refl_0_44um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_44um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(9)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(9)%Ref_Toa
         case("refl_0_75um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_75um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(15)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(15)%Ref_Toa
         case("refl_0_76um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_76um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(45)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(45)%Ref_Toa
         case("refl_0_95um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_95um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(17)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(17)%Ref_Toa
         case("refl_0_93um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_93um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(18)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(18)%Ref_Toa
         case("refl_0_94um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_0_94um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(19)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(19)%Ref_Toa
         case("refl_1_38um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_1_38um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(26)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(26)%Ref_Toa
         case("refl_sol_dnb_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_solar_dnb_nominal"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(44)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(44)%Ref_Toa
         case("refl_lunar_dnb_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_lunar_dnb_nominal"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(44)%Ref_Lunar_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(44)%Ref_Lunar_Toa
         case("refl_3_75um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_3_75um_nom"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(20)%Ref_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(20)%Ref_Toa
         case("temp_3_75um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_3_75um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(20)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(20)%Bt_Toa
         case("temp_3_75um_nom_sounder")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_3_75um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Bt_375um_Sounder)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Bt_375um_Sounder
         case("temp_3_75um_nom_min_sub")
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(20)%Bt_Toa_Min_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(20)%Bt_Toa_Min_Sub
         case("temp_3_75um_nom_max_sub")
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(20)%Bt_Toa_Max_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(20)%Bt_Toa_Max_Sub
         case("temp_3_9um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_3_9um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(21)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(21)%Bt_Toa
         case("temp_4_05um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_4_05um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(23)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(23)%Bt_Toa
         case("temp_4_45um_nom_sounder")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_4_45um_nom"
            Sds_Info(Var_Idx)%Units =  "K" 
            if (allocated(Bt_445um_Sounder)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Bt_445um_Sounder
         case("temp_4_46um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_4_46um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(24)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(24)%Bt_Toa
         case("temp_4_52um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_4_52um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(25)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(25)%Bt_Toa
         case("temp_4_57um_nom_sounder")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_4_57um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Bt_457um_Sounder)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Bt_457um_Sounder
         case("temp_6_7um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_6_7um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(27)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(27)%Bt_Toa
         case("temp_7_3um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_7_3um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(28)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(28)%Bt_Toa
         case("temp_8_5um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_8_5um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(29)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(29)%Bt_Toa
         case("temp_9_7um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_9_7um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(30)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(30)%Bt_Toa
         case("radiance_11_0um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_radiance_11_0um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,150.0]
            Sds_Info(Var_Idx)%Units = "mW/m^2/cm^-1"
            if (allocated(ch(31)%Rad_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Rad_Toa
         case("radiance_11_0um_nom_min_sub")
            Sds_Info(Var_Idx)%Standard_Name = "radiance_11_0um_nom_min_sub"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,150.0]
            Sds_Info(Var_Idx)%Units = "mW/m^2/cm^-1"
            if (allocated(ch(31)%Rad_Toa_Min_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Rad_Toa_Min_Sub
         case("radiance_11_0um_nom_max_sub")
            Sds_Info(Var_Idx)%Standard_Name = "radiance_11_0um_nom_max_sub"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,150.0]
            Sds_Info(Var_Idx)%Units = "mW/m^2/cm^-1"
            if (allocated(ch(31)%Rad_Toa_Max_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Rad_Toa_Max_Sub
         case("radiance_6_7um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_radiance_6_7um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,75.0]
            Sds_Info(Var_Idx)%Units = "mW/m^2/cm^-1"
            if (allocated(ch(27)%Rad_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(27)%Rad_Toa
       
         case("temp_11_0um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_11_0um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(31)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Bt_Toa
         case("temp_11_0um_nom_min_sub")
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(31)%Bt_Toa_Min_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Bt_Toa_Min_Sub
         case("temp_11_0um_nom_max_sub")
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(31)%Bt_Toa_Max_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Bt_Toa_Max_Sub
         case("temp_11_0um_nom_mean_sub")
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(31)%Bt_Toa_Mean_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Bt_Toa_Mean_Sub
         case("temp_11_0um_nom_stddev_sub")
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(31)%Bt_Toa_Std_Sub)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Bt_Toa_Std_Sub
         case("temp_11_0um_nom_sounder")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_11_0um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Bt_11um_Sounder)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Bt_11um_Sounder
         case("temp_11_0um_nom_max_3x3_sub")
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Bt_Ch43_Max_Sub_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Bt_Ch43_Max_Sub_3x3
         case("temp_12_0um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_12_0um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(32)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(32)%Bt_Toa
         case("temp_12_0um_nom_sounder")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_12_0um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Bt_12um_Sounder)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Bt_12um_Sounder
         case("temp_13_3um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_13_3um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(33)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(33)%Bt_Toa
         case("temp_13_6um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_13_6um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(34)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(34)%Bt_Toa
         case("temp_13_9um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_13_9um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(35)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(35)%Bt_Toa
         case("temp_14_2um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_14_2um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(36)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(36)%Bt_Toa
         case("temp_14_5um_nom_sounder")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_14_5um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Bt_145um_Sounder)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Bt_145um_Sounder
         case("temp_14_7um_nom_sounder")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_14_7um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Bt_147um_Sounder)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Bt_147um_Sounder
         case("temp_14_9um_nom_sounder")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_14_9um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Bt_149um_Sounder)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Bt_149um_Sounder
         case("temp_6_2um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_6_2um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(37)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(37)%Bt_Toa
         case("temp_10_4um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_10_4um_nom"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(38)%Bt_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(38)%Bt_Toa
         case("refl_lunar_dnb_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_assuming_clear_sky_lunar_dnb_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(44)%Ref_Lunar_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 =>  Ch(44)%Ref_Lunar_Toa_Clear
         case("refl_0_65um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_assuming_clear_sky_0_65_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(1)%Ref_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 =>  Ch(1)%Ref_Toa_Clear
         case("refl_3_75um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_assuming_clear_sky_3_75_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(20)%Ref_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 =>  Ch(20)%Ref_Toa_Clear
         case("temp_3_75um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_3_75_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(20)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(20)%Bt_Toa_Clear
         case("temp_6_7um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_6_7_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(27)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(27)%Bt_Toa_Clear
         case("temp_7_3um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_7_3_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(28)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(28)%Bt_Toa_Clear
         case("temp_8_5um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_8_5_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(29)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(29)%Bt_Toa_Clear
         case("temp_9_7um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_9_7_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(30)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(30)%Bt_Toa_Clear
         case("temp_11_0um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_11_0_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(31)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Bt_Toa_Clear
         case("temp_12_0um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_12_0_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(32)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(32)%Bt_Toa_Clear
         case("temp_13_3um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_13_3_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(33)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(33)%Bt_Toa_Clear
         case("temp_13_6um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_13_6_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(34)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(34)%Bt_Toa_Clear
         case("temp_13_9um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_13_9_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(35)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(35)%Bt_Toa_Clear
         case("temp_14_2um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_14_2_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(36)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(36)%Bt_Toa_Clear
         case("temp_6_2um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_6_2_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(37)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(37)%Bt_Toa_Clear
         case("temp_10_4um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_brightness_temperature_assuming_clear_sky_10_4_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ch(38)%Bt_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(38)%Bt_Toa_Clear
         case("cld_opd_mask_0_65um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_optical_thickness_due_to_cloud_assuming_water_phase"
            Sds_Info(Var_Idx)%Actual_Range = [-20.0,100.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(1)%Opd)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(1)%Opd
         case("cld_opd_mask_dnb_nom")
            Sds_Info(Var_Idx)%Standard_Name = "cld_opd_mask_dnb_nom"
            Sds_Info(Var_Idx)%Actual_Range = [-20.0,100.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(44)%Opd)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(44)%Opd
         case("cld_opd_mask_1_38um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_optical_thickness_due_to_cloud_assuming_water_phase"
            Sds_Info(Var_Idx)%Standard_Name = "cld_opd_mask_1_38um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [-20.0,100.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(26)%Opd)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(26)%Opd
         case("cld_opd_mask_1_60um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "cld_opd_mask_1_60um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [-20.0,100.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(6)%Opd)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(6)%Opd
         case("refl_sfc_white_sky_0_65um_nom")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(1)%Sfc_Ref_White_Sky)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(1)%Sfc_Ref_White_Sky
         case("refl_sfc_white_sky_0_86um_nom")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(2)%Sfc_Ref_White_Sky)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(2)%Sfc_Ref_White_Sky
         case("refl_sfc_white_sky_1_20um_nom")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(5)%Sfc_Ref_White_Sky)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(5)%Sfc_Ref_White_Sky
         case("refl_sfc_white_sky_1_60um_nom")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(6)%Sfc_Ref_White_Sky)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(6)%Sfc_Ref_White_Sky
         case("refl_sfc_white_sky_2_10um_nom")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(7)%Sfc_Ref_White_Sky)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(7)%Sfc_Ref_White_Sky
         case("emiss_tropo_3_75um_nom")
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(20)%Emiss_Tropo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(20)%Emiss_Tropo
         case("emiss_tropo_6_2um_nom")
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(37)%Emiss_Tropo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(37)%Emiss_Tropo
         case("emiss_tropo_6_7um_nom")
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(27)%Emiss_Tropo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(27)%Emiss_Tropo
         case("emiss_tropo_8_5um_nom")
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(29)%Emiss_Tropo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(29)%Emiss_Tropo
         case("emiss_tropo_10_4um_nom")
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(38)%Emiss_Tropo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(38)%Emiss_Tropo
         case("emiss_tropo_11_0um_nom")
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(31)%Emiss_Tropo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Emiss_Tropo
         case("emiss_tropo_12_0um_nom")
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(32)%Emiss_Tropo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(32)%Emiss_Tropo
         case("emiss_tropo_13_3um_nom")
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(33)%Emiss_Tropo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(33)%Emiss_Tropo
         case("trans_atm_11_0um_nom")
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(31)%Trans_Atm)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Trans_Atm

         case("csbt_mask_6_2um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "csbt_mask_6_2um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(37)%CSBT_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(37)%CSBT_Mask
         case("csbt_mask_6_7um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "csbt_mask_6_7um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(27)%CSBT_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(27)%CSBT_Mask
         case("csbt_mask_7_3um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "csbt_mask_7_3um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(28)%CSBT_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(28)%CSBT_Mask
         case("csbt_mask_8_5um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "csbt_mask_8_5um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(29)%CSBT_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(29)%CSBT_Mask
         case("csbt_mask_9_7um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "csbt_mask_9_7um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(30)%CSBT_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(30)%CSBT_Mask
         case("csbt_mask_10_4um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "csbt_mask_10_4um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(38)%CSBT_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(38)%CSBT_Mask
         case("csbt_mask_11_0um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "csbt_mask_11_0um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(31)%CSBT_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(31)%CSBT_Mask
         case("csbt_mask_12_0um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "csbt_mask_12_0um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(32)%CSBT_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(32)%CSBT_Mask
         case("csbt_mask_13_3um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "csbt_mask_13_3um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(33)%CSBT_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(33)%CSBT_Mask
         case("opaque_height_6_2um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_6_2um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(37)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(37)%Opaque_Height
         case("opaque_height_6_7um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_6_2um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(27)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(27)%Opaque_Height
         case("opaque_height_7_3um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_7_3um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(28)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(28)%Opaque_Height
         case("opaque_height_8_5um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_8_5um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(29)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(29)%Opaque_Height
         case("opaque_height_9_7um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_9_7um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(30)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(30)%Opaque_Height
         case("opaque_height_10_4um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_10_4um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(38)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(38)%Opaque_Height
         case("opaque_height_11_0um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_11_0um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(31)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(31)%Opaque_Height
         case("opaque_height_12_0um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_12_0um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(32)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(32)%Opaque_Height
         case("opaque_height_13_3um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_13_3um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(33)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(33)%Opaque_Height
         case("opaque_height_13_6um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_13_6um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(34)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(34)%Opaque_Height
         case("opaque_height_13_9um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_13_9um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(35)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(35)%Opaque_Height
         case("opaque_height_14_2um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "opaque_height_14_2um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units =  "m"
            if (allocated(ch(36)%Opaque_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(36)%Opaque_Height
         case("emiss_sfc_3_75um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "emiss_sfc_3_75um_nom"
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(20)%Sfc_Emiss)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(20)%Sfc_Emiss
         case("emiss_sfc_8_5um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "emiss_sfc_8_5um_nom"
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(29)%Sfc_Emiss)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(29)%Sfc_Emiss
         case("emiss_sfc_11_0um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "emiss_sfc_11_0um_nom"
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ch(31)%Sfc_Emiss)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(31)%Sfc_Emiss

         !--- data source flags
         case("3_75um_nom_source")
            Sds_Info(Var_Idx)%Standard_Name = "3_75um_nom_source"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(20)%Source)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(20)%Source

         !--- data quality flags 
         case("dqf_3_75um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "dqf_3_75um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(20)%DQF)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(20)%DQF
         case("dqf_6_2um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "dqf_6_2um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(37)%DQF)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(37)%DQF
         case("dqf_6_7um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "dqf_6_7um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(27)%DQF)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(27)%DQF
         case("dqf_7_3um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "dqf_7_3um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(28)%DQF)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(28)%DQF
         case("dqf_8_5um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "dqf_8_5um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(29)%DQF)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(29)%DQF
         case("dqf_9_7um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "dqf_9_7um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(30)%DQF)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(30)%DQF
         case("dqf_10_4um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "dqf_10_4um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(38)%DQF)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(38)%DQF
         case("dqf_11_0um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "dqf_11_0um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(31)%DQF)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(31)%DQF
         case("dqf_12_0um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "dqf_12_0um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(32)%DQF)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(32)%DQF
         case("dqf_13_3um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "dqf_13_3um_nom"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(ch(33)%DQF)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ch(33)%DQF
         !----------------------------------------------------------------------------------------------------
         ! Sfc Members
         !----------------------------------------------------------------------------------------------------
         case("refl_0_47um_nom_atmos_corr")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(3)%Ref_Sfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(3)%Ref_Sfc
         case("refl_0_55um_nom_atmos_corr")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(4)%Ref_Sfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(4)%Ref_Sfc
         case("refl_0_65um_nom_atmos_corr")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(1)%Ref_Sfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(1)%Ref_Sfc
         case("refl_0_86um_nom_atmos_corr")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(2)%Ref_Sfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(2)%Ref_Sfc
         case("refl_1_60um_nom_atmos_corr")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(6)%Ref_Sfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(6)%Ref_Sfc
         case("refl_2_10um_nom_atmos_corr")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(7)%Ref_Sfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(7)%Ref_Sfc
         case("refl_3_75um_nom_atmos_corr")
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(ch(20)%Ref_Sfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(20)%Ref_Sfc
         !----------------------------------------------------------------------------------------------------
         ! instrument counts
         !----------------------------------------------------------------------------------------------------
         case("refl_0_65um_nom_counts")
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1024.0]
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT16
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF = DFNT_INT16
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_SHORT
            Sds_Info(Var_Idx)%Units =  "none"
            Sds_Info(Var_Idx)%Standard_Name = "refl_0_65um_nom_counts"
            Sds_Info(Var_Idx)%Long_Name = "refl_0_65um_nom_counts"
            if (allocated(Ch1_Counts)) Sds_Info(Var_Idx)%Sds_Data_2d_I2 => Ch1_Counts
         case("refl_0_86um_nom_counts")
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1024.0]
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT16
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF = DFNT_INT16
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_SHORT
            Sds_Info(Var_Idx)%Units =  "none"
            Sds_Info(Var_Idx)%Standard_Name = "refl_0_86um_nom_counts"
            Sds_Info(Var_Idx)%Long_Name = "refl_0_86um_nom_counts"
            if (allocated(Ch2_Counts)) Sds_Info(Var_Idx)%Sds_Data_2d_I2 => Ch2_Counts
         case("refl_1_60um_nom_counts")
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1024.0]
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT16
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF = DFNT_INT16
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_SHORT
            Sds_Info(Var_Idx)%Units =  "none"
            Sds_Info(Var_Idx)%Standard_Name = "refl_1_60um_nom_counts"
            Sds_Info(Var_Idx)%Long_Name = "refl_1_60um_nom_counts"
            if (allocated(Ch6_Counts)) Sds_Info(Var_Idx)%Sds_Data_2d_I2 => Ch6_Counts
         !----------------------------------------------------------------------------------------------------
         ! 3x3 uniformity arrays
         !----------------------------------------------------------------------------------------------------
         case("refl_0_65um_nom_min_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "bidirectional_reflectance_0_65_micron_nom_min_3x3"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(1)%Ref_Toa_Min_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(1)%Ref_Toa_Min_3x3
         case("refl_0_65um_nom_max_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "bidirectional_reflectance_0_65_micron_nom_max_3x3"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(1)%Ref_Toa_Max_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(1)%Ref_Toa_Max_3x3
         case("refl_0_65um_nom_stddev_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "bidirectional_reflectance_0_65_micron_nom_stddev_3x3"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(1)%Ref_Toa_Std_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(1)%Ref_Toa_Std_3x3
         case("refl_1_38um_nom_stddev_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "bidirectional_reflectance_1_38_micron_nom_stddev_3x3"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(26)%Ref_Toa_Std_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(26)%Ref_Toa_Std_3x3
         case("temp_3_75um_nom_stddev_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "brightness_temperature_3_75_micron_nom_stddev_3x3"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Ch(20)%Bt_Toa_Std_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(20)%Bt_Toa_Std_3x3
         case("temp_11_0um_nom_stddev_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "brightness_temperature_11_0_micron_nom_stddev_3x3"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Ch(31)%Bt_Toa_Std_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(31)%Bt_Toa_Std_3x3
         case("temp_11_0um_nom_min_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "brightness_temperature_11_0_micron_nom_min_3x3"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Ch(31)%Bt_Toa_Min_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(31)%Bt_Toa_Min_3x3
         case("temp_11_0um_nom_max_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "brightness_temperature_11_0_micron_nom_max_3x3"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Ch(31)%Bt_Toa_Max_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(31)%Bt_Toa_Max_3x3
         case("refl_lunar_dnb_nom_mean_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "bidirectional_reflectance_lunar_dnb_nom_mean_3x3"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(44)%Ref_Lunar_Mean_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(44)%Ref_Lunar_Mean_3x3
         case("refl_lunar_dnb_nom_min_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "bidirectional_reflectance_lunar_dnb_nom_min_3x3"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(44)%Ref_Lunar_Min_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(44)%Ref_Lunar_Min_3x3
         case("refl_lunar_dnb_nom_max_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "bidirectional_reflectance_lunar_dnb_nom_max_3x3"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(44)%Ref_Lunar_Max_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(44)%Ref_Lunar_Max_3x3
         case("refl_lunar_dnb_nom_stddev_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "bidirectional_reflectance_lunar_dnb_nom_stddev_3x3"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(44)%Ref_Lunar_Std_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(44)%Ref_Lunar_Std_3x3
         case("temp_10_4um_nom_stddev_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "brightness_temperature_10_4_micron_nom_stddev_3x3"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Ch(38)%Bt_Toa_Std_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(38)%Bt_Toa_Std_3x3
         case("temp_10_4um_nom_min_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "brightness_temperature_10_4_micron_nom_min_3x3"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Ch(38)%Bt_Toa_Min_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(38)%Bt_Toa_Min_3x3
         case("temp_10_4um_nom_max_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "brightness_temperature_10_4_micron_nom_max_3x3"
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(Ch(38)%Bt_Toa_Max_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(38)%Bt_Toa_Max_3x3
         !----------------------------------------------------------------------------------------------------
         ! CLDMASK Members
         !----------------------------------------------------------------------------------------------------
         case("cloud_mask")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Long_Name = "integer classification of the cloud mask including clear=0, probably-clear=1, "// &
                                          "probably-cloudy=2, cloudy=3"
            Sds_Info(Var_Idx)%Number_Of_Flags = 4
            Sds_Info(Var_Idx)%Flags_String = "clear "// &
                                      " probably_clear "// &
                                      " probably_cloudy "// &
                                      " cloudy "
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Cld_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Cld_Mask
         case("cloud_mask_ir")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_mask_ir"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Long_Name = "integer classification of the cloud mask including clear=0, probably-clear=1, "// &
                                          "probably-cloudy=2, cloudy=3 from the infrared observations"
            Sds_Info(Var_Idx)%Number_Of_Flags = 4
            Sds_Info(Var_Idx)%Flags_String = "clear "// &
                                      " probably_clear "// &
                                      " probably_cloudy "// &
                                      " cloudy "
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Cld_Mask_IR)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Cld_Mask_IR
         case("cloud_mask_aux")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_mask_aux"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Cld_Mask_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Cld_Mask_Aux
         case("cloud_mask_binary")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_mask_binary"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Cld_Mask_Binary)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Cld_Mask_Binary
         case("cloud_mask_binary_ir")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_mask_binary_ir"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Cld_Mask_Binary_IR)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Cld_Mask_Binary_IR
         case("adj_pix_cld_mask")
            Sds_Info(Var_Idx)%Standard_Name = "adj_pix_cld_mask"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Adj_Pix_Cld_Mask)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Adj_Pix_Cld_Mask
         case("cloud_probability")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_probability"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Posterior_Cld_Probability)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CLDMASK%Posterior_Cld_Probability
         case("cloud_probability_uncertainty")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_probability_uncertainty"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,0.5]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Posterior_Cld_Probability_Uncer)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CLDMASK%Posterior_Cld_Probability_Uncer
         case("cloud_probability_ir")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_probability_ir"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Posterior_Cld_Probability_IR)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CLDMASK%Posterior_Cld_Probability_IR
         case("prior_cloud_probability")
            Sds_Info(Var_Idx)%Standard_Name = "prior_cloud_probability"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Prior_Cld_Probability)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CLDMASK%Prior_Cld_Probability
         case("cloud_probability_aux")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_probability_aux"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Posterior_Cld_Probability_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CLDMASK%Posterior_Cld_Probability_Aux
         case("bayes_mask_sfc_type")
            Sds_Info(Var_Idx)%Standard_Name = "bayes_mask_sfc_type"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Bayes_Mask_Sfc_Type)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Bayes_Mask_Sfc_Type
         case("cloud_mask_test_packed_results")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_mask_test_packed_results"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Rank = 3
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Cld_Test_Vector_Packed)) Sds_Info(Var_Idx)%Sds_Data_3d_I1 => CLDMASK%Cld_Test_Vector_Packed
         case("cloud_mask_tut")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_mask_tut"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%TUT)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%TUT
         case("cloud_mask_rut")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_mask_rut"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%RUT)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%RUT
         case("ice_cloud_probability")
            Sds_Info(Var_Idx)%Standard_Name = "ice_cloud_probability"
            Sds_Info(Var_Idx)%Long_Name = "ice cloud probability"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Posterior_Ice_Probability)) &
               Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CLDMASK%Posterior_Ice_Probability
         case("water_cloud_probability")
            Sds_Info(Var_Idx)%Standard_Name = "water_cloud_probability"
            Sds_Info(Var_Idx)%Long_Name = "water cloud probability"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(CLDMASK%Posterior_Water_Probability)) &
               Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CLDMASK%Posterior_Water_Probability
         case("cloud_mask_qf")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_mask_qf"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CLDMASK%Cld_Mask_Qf)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CLDMASK%Cld_Mask_Qf
         !----------------------------------------------------------------------------------------------------
         ! CLDPHASE Members
         !----------------------------------------------------------------------------------------------------
         case("cloud_phase")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_phase"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Long_Name = "integer classification of the cloud phase including clear "// &
                               "and aerosol type,0=clear,1=water,2=supercooled water,3=mixed,4=ice,5=unknown"
            Sds_Info(Var_Idx)%Number_Of_Flags = 6
            Sds_Info(Var_Idx)%Flags_String = "clear "// &
                                      " water "// &
                                      " supercooled water "// &
                                      " mixed "// &
                                      " ice "// &
                                      " unknown "
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(Cld_Phase)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Cld_Phase
         case("cloud_phase_uncertainty")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_phase_uncertainty"
            Sds_Info(Var_Idx)%Long_Name = "cloud phase uncertainty"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(Cld_Phase_Uncertainty)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Cld_Phase_Uncertainty
         case("cloud_type")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_type"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Long_Name = "integer classification of the cloud type including clear "// &
                               "and aerosol type,0=clear,1=probably clear,2=fog,3=water,4=supercooled water,"//&
                               "5=mixed,6=opaque_ice,7=cirrus,8=overlapping,9=overshooting,10=unknown,11=dust,12=smoke,13=fire"
            Sds_Info(Var_Idx)%Number_Of_Flags = 14
            Sds_Info(Var_Idx)%Flags_String = "clear "// &
                                      " probably clear "// &
                                      " fog "// &
                                      " water "// &
                                      " supercooled water "// &
                                      " mixed "// &
                                      " opaque_ice "// &
                                      " cirrus "// &
                                      " overlapping "// &
                                      " overshooting "// &
                                      " unknown "// &
                                      " dust "// &
                                      " smoke "// &
                                      " fire "
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(Cld_Type)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Cld_Type
         case("cloud_phase_ir")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_phase_ir"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(Cld_Phase_Ir)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Cld_Phase_Ir
         case("cloud_type_ir")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_type_ir"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(Cld_Type_Ir)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Cld_Type_Ir
         case("cloud_phase_aux")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_phase_aux"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(Cld_Phase_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Cld_Phase_Aux
         case("cloud_type_aux")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_type_aux"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(Cld_Type_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Cld_Type_Aux

         !----------------------------------------------------------------------------------------------------
         ! ACHA Members
         !----------------------------------------------------------------------------------------------------
         case("cld_temp_prior_acha")
            Sds_Info(Var_Idx)%Standard_Name = "air_temperature_at_cloud_top_prior_value"
            Sds_Info(Var_Idx)%Long_Name = "air_temperature_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ACHA%Tc_Ap)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Tc_Ap
         case("cld_temp_uncer_prior_acha")
            Sds_Info(Var_Idx)%Standard_Name = "air_temperature_uncertainty_at_cloud_top_prior_value"
            Sds_Info(Var_Idx)%Long_Name = "air_temperature_uncertainty_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100.0]
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ACHA%Tc_Ap_Uncer)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Tc_Ap_Uncer
         case("lower_cld_temp_prior_acha")
            Sds_Info(Var_Idx)%Standard_Name = "air_temperature_at_cloud_top_prior_value"
            Sds_Info(Var_Idx)%Long_Name = "air_temperature_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ACHA%Lower_Tc_Ap)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Lower_Tc_Ap
         case("cld_temp_acha")
            Sds_Info(Var_Idx)%Standard_Name = "air_temperature_at_cloud_top"
            Sds_Info(Var_Idx)%Long_Name = "air_temperature_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ACHA%Tc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Tc
         case("cld_temp_eff_acha")
            Sds_Info(Var_Idx)%Standard_Name = "air_temperature_at_eff_cloud_top"
            Sds_Info(Var_Idx)%Long_Name = "air_temperature_at_eff_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ACHA%Tc_Eff)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Tc_Eff
         case("cld_temp_lower_acha")
            Sds_Info(Var_Idx)%Standard_Name = "air_temperature_of_lower_cloud"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Long_Name = "air_temperature_at_cloud_top_of_lower_cloud"
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(ACHA%Lower_Tc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Lower_Tc
         case("cld_temp_uncer_acha")
            Sds_Info(Var_Idx)%Standard_Name = "cld_temp_uncer_acha"
            Sds_Info(Var_Idx)%Long_Name = "cld_temp_uncer_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100.0]
            Sds_Info(Var_Idx)%Units =  "K"
            if (allocated(ACHA%Tc_Uncertainty)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Tc_Uncertainty
         case("cld_emiss_acha")
            Sds_Info(Var_Idx)%Standard_Name = "convective_cloud_longwave_emissivity"
            Sds_Info(Var_Idx)%Long_Name = "convective_cloud_longwave_emissivity"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,1.2]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ACHA%Ec)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Ec
         case("cld_emiss_prior_acha")
            Sds_Info(Var_Idx)%Standard_Name = "convective_cloud_longwave_emissivity_prior"
            Sds_Info(Var_Idx)%Long_Name = "convective_cloud_longwave_emissivity_prior"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,1.2]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ACHA%Ec_Ap)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Ec_Ap
         case("cld_emiss_acha_uncer")
            Sds_Info(Var_Idx)%Standard_Name = "convective_cloud_longwave_emissivity_uncertainty"
            Sds_Info(Var_Idx)%Long_Name = "convective_cloud_longwave_emissivity_uncertainty"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ACHA%Ec)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Ec_Uncertainty
         case("cld_beta_acha")
            Sds_Info(Var_Idx)%Standard_Name = "cld_beta_acha"
            Sds_Info(Var_Idx)%Long_Name = "cld_beta_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.2,2.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ACHA%Beta)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Beta
         case("cld_beta_prior_acha")
            Sds_Info(Var_Idx)%Standard_Name = "cld_beta_prior_acha"
            Sds_Info(Var_Idx)%Long_Name = "cld_beta_prior_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.2,2.0]
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(ACHA%Beta_Ap)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Beta_Ap
         case("cld_press_acha")
            Sds_Info(Var_Idx)%Standard_Name = "air_pressure_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [50.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure from ACHA"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(ACHA%Pc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Pc
         case("cld_press_eff_acha")
            Sds_Info(Var_Idx)%Standard_Name = "air_pressure_at_effective_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [50.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "Effective Cloud Top Pressure from ACHA"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(ACHA%Pc_Eff)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Pc_Eff
         case("cld_press_lower_acha")
            Sds_Info(Var_Idx)%Standard_Name = "air_pressure_at_effective_cloud_top_lower_layer"
            Sds_Info(Var_Idx)%Actual_Range = [50.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "Effective Cloud Top Pressure for lower layer from ACHA"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(ACHA%Lower_Pc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Lower_Pc
         case("cld_press_uncer_acha")
            Sds_Info(Var_Idx)%Standard_Name = "cld_press_uncer_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure Uncertainty from ACHA"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(ACHA%Pc_Uncertainty)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Pc_Uncertainty
         case("cld_press_uncer1_aux")
            Sds_Info(Var_Idx)%Standard_Name = "cld_press_uncer1_aux"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure Uncertainty from AUX"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(Pc_Uncertainty1_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Pc_Uncertainty1_Aux
         case("cld_press_uncer2_aux")
            Sds_Info(Var_Idx)%Standard_Name = "cld_press_uncer2_aux"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure Uncertainty from layer 2 AUX"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(Pc_Uncertainty2_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Pc_Uncertainty2_Aux
         case("cld_height_acha")
            Sds_Info(Var_Idx)%Standard_Name = "height_at_cloud_top"
            !Sds_Info(Var_Idx)%Actual_Range = [75.0,20000.0]
            Sds_Info(Var_Idx)%Actual_Range = [0.0,40000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Height from ACHA"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(ACHA%Zc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Zc
         case("cld_height_eff_acha")
            Sds_Info(Var_Idx)%Standard_Name = "height_at_eff_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [75.0,20000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Eff Height from ACHA"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(ACHA%Zc_Eff)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Zc_Eff
         case("cld_height_lower_acha")
            Sds_Info(Var_Idx)%Standard_Name = "height_of_lower_cloud"
            Sds_Info(Var_Idx)%Actual_Range = [75.0,20000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Height of Lower Level Cloud from ACHA"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(ACHA%Lower_Zc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Lower_Zc
         case("cld_height_uncer_acha")
            Sds_Info(Var_Idx)%Standard_Name = "cld_height_uncer_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,10000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Height Uncertainty from ACHA"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(ACHA%Zc_Uncertainty)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Zc_Uncertainty
         case("cld_height_base_acha")
            Sds_Info(Var_Idx)%Standard_Name = "base_height_of_cloud"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Base Height from ACHA"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(ACHA%Zc_Base)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Zc_Base
         case("cld_altitude_acha")
            Sds_Info(Var_Idx)%Standard_Name = "altitude_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Altitude from ACHA"
            Sds_Info(Var_Idx)%Units = "feet"
            if (allocated(ACHA%Alt)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Alt
         case("lower_cld_altitude_acha")
            Sds_Info(Var_Idx)%Standard_Name = "altitude_at_lower_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100000.0]
            Sds_Info(Var_Idx)%Long_Name = "Lower Cloud Top Altitude from ACHA"
            Sds_Info(Var_Idx)%Units = "feet"
            if (allocated(ACHA%Lower_Alt)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Lower_Alt
         case("supercooled_prob_acha")
            Sds_Info(Var_Idx)%Standard_Name = "supercooled_cloud_probability_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "supercooled_cloud_probability_acha"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Supercooled_Cld_Prob)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Supercooled_Cld_Prob
         case("ice_prob_acha")
            Sds_Info(Var_Idx)%Standard_Name = "ice_probability_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "ice_probability_acha"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Ice_Probability)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Ice_Probability
         case("ice_prob_prior_acha")
            Sds_Info(Var_Idx)%Standard_Name = "prior_ice_probability_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "prior_ice_probability_acha"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Ice_Prob_Ap)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Ice_Prob_Ap
         case("ice_prob_uncer_acha")
            Sds_Info(Var_Idx)%Standard_Name = "ice_probability_uncertainty_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "ice_probability_uncertainty_acha"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Ice_Probability)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Ice_Probability
         case("acha_quality")
            Sds_Info(Var_Idx)%Standard_Name = "acha_quality"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Number_Of_Flags = 7
            Sds_Info(Var_Idx)%Flags_String =  "Processed "// &
                               " valid_Tc_retrieval "// &
                               " valid_ec_retrieval "// &
                               " valid_beta_retrieval "// &
                               " degraded_Tc_retrieval "// &
                               " degraded_ec_retrieval "// &
                               " degraded_beta_retrieval "
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Packed_Quality_Flags)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ACHA%Packed_Quality_Flags

         case("acha_info")
            Sds_Info(Var_Idx)%Standard_Name = "acha_info"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Number_Of_Flags = 7
            Sds_Info(Var_Idx)%Flags_String = "Cloud_Height_Attempted "// &
                               " Bias_Correction_Employed "// &
                               " Ice_Cloud_Retrieval "// &
                               " Local_Radiative_Center_Processing_Used "// &
                               " Multi_Layer_Retrieval "// &
                               " Lower_Cloud_Interpolation_Used "// &
                               " Boundary_Layer_Inversion_Assumed "
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Packed_Meta_Data_Flags)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ACHA%Packed_Meta_Data_Flags
         case("cost_acha")
            Sds_Info(Var_Idx)%Standard_Name = "oe_cost_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Cost)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Cost
         case("cost_aux")
            Sds_Info(Var_Idx)%Standard_Name = "oe_cost_aux"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,2000.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Cost_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Cost_Aux
         case("goodness_acha")
            Sds_Info(Var_Idx)%Standard_Name = "oe_goodness_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Goodness)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Goodness
         case("cloud_type_acha")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_type_acha"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units =  "none"
            if (allocated(Acha%Cld_Type)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Acha%Cld_Type
         case("metadata_aux")
            Sds_Info(Var_Idx)%Standard_Name = "metadata_aux"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF =  NF90_BYTE
!           Sds_Info(Var_Idx)%Long_Name = "metadata from auxilliary data (sensor specific)"
            Sds_Info(Var_Idx)%Long_Name = "modis ctp: " // &
                                          "1 -- CO2-slicing retrieval, bands 36/35 "// &
                                          "2 -- CO2-slicing retrieval, bands 35/34 3 -- CO2-slicing retrieval, bands 35/33 "// &
                                          "4 --CO2-slicing retrieval, bands 34/33 6 -- IR-window retrieval, band 31"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Metadata_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Metadata_Aux
         case("cld_height_aux")
            Sds_Info(Var_Idx)%Standard_Name = "height_at_cloud_top_aux"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(Zc_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Zc_Aux
         case("cld_emiss_aux")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_emissivity_aux"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,1.2]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Ec_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ec_Aux
         case("cld_temp_aux")
            Sds_Info(Var_Idx)%Standard_Name = "air_temperature_at_cloud_top_aux"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(Tc_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tc_Aux
         case("cld_press_aux")
            Sds_Info(Var_Idx)%Standard_Name = "air_pressure_at_cloud_top_aux"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(Pc_Top1_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Pc_Top1_Aux
         case("cld_press_l2_aux")
            Sds_Info(Var_Idx)%Standard_Name = "air_pressure_at_layer2_cloud_top_aux"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(Pc_Top2_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Pc_Top2_Aux
         case("acha_inversion_flag")
            Sds_Info(Var_Idx)%Standard_Name = "acha_inversion_flag"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Long_Name = "acha inversion flag (0=no,1=inversion)"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Inversion_Flag)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => ACHA%Inversion_Flag
         case("cld_opd_acha")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_optical_depth_from_acha"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,8.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Tau)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Tau
         case("cld_opd_uncer_acha")
            Sds_Info(Var_Idx)%Standard_Name = "uncertainty_of_cloud_optical_depth_from_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,8.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Tau_Uncer)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Tau_Uncer
         case("cld_reff_acha")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_effective_radius_from_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,160.0]
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(ACHA%Reff)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Reff
         case("conv_cloud_probability")
            Sds_Info(Var_Idx)%Standard_Name = "convective_cloud_probability_acha"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "convective_cloud_probability_acha"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(ACHA%Conv_Cld_Prob)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ACHA%Conv_Cld_Prob

         !----------------------------------------------------------------------------------------------------
         ! Other ACHA Related Variables
         !----------------------------------------------------------------------------------------------------
         case("cld_height_co2irw")
            Sds_Info(Var_Idx)%Standard_Name = "height_at_cloud_top_from_co2irw"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Height from CO2 and IRW method"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(Zc_CO2IRW)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Zc_CO2IRW
         case("cld_height_h2o")
            Sds_Info(Var_Idx)%Standard_Name = "height_at_cloud_top_from_h2o_intercept"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Height from H2O Intercept"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(Zc_H2O)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Zc_H2O
         case("cld_height_opaque")
            Sds_Info(Var_Idx)%Standard_Name = "height_at_cloud_top_assuming_opaque"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Height Assuming Opaque Solution"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(Zc_Opaque_Cloud)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Zc_Opaque_Cloud
         case("cld_press_opaque")
            Sds_Info(Var_Idx)%Standard_Name = "pressure_at_cloud_top_assuming_opaque"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure Assuming Opaque Solution"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(Pc_Opaque_Cloud)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Pc_Opaque_Cloud
         case("cld_press_h2o")
            Sds_Info(Var_Idx)%Standard_Name = "pressure_at_cloud_top_assuming_opaque"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure H2O Solution"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(Pc_H2O)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Pc_H2O
!        case("cld_press_em")
!           Sds_Info(Var_Idx)%Standard_Name = "pressure_at_cloud_top_assuming_eyre_menzel"
!           Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
!           Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure Eyre Menzel Solution"
!           Sds_Info(Var_Idx)%Units = "hPa"
!           if (allocated(Pc_EM)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Pc_EM
         case("cld_press_qrnn")
            Sds_Info(Var_Idx)%Standard_Name = "pressure_at_cloud_top_qrnn"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure QRNN Solution"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(QRNN_CTP)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => QRNN_CTP
         case("cld_press_uncer_qrnn")
            Sds_Info(Var_Idx)%Standard_Name = "pressure_uncertainty_at_cloud_top_qrnn"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,500.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure Uncertainty QRNN Solution"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(QRNN_CTP_Uncer)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => QRNN_CTP_Uncer
         case("cld_press_co2irw")
            Sds_Info(Var_Idx)%Standard_Name = "pressure_at_cloud_top_from_co2irw"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure CO2IRW Solution"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(Pc_CO2IRW)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Pc_CO2IRW
         case("cld_press_splitwin")
            Sds_Info(Var_Idx)%Standard_Name = "pressure_at_cloud_top_from_splitwin"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure Split Window Solution"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(Pc_Splitwin)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Pc_Splitwin
         case("cld_emiss_h2o")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_emissivity_from_h2o"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.2]
            if (allocated(Ec_H2O)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ec_H2O
         case("cld_emiss_co2irw")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_emissivity_from_co2irw"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.2]
            if (allocated(Ec_CO2IRW)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ec_CO2IRW
         case("cld_emiss_splitwin")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_emissivity_from_splitwin"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.2]
            if (allocated(Ec_Splitwin)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ec_Splitwin
!        case("cld_emiss_em")
!           Sds_Info(Var_Idx)%Standard_Name = "cloud_emissivity_from_eyre_menzel"
!           Sds_Info(Var_Idx)%Actual_Range = [0.0,1.2]
!           if (allocated(Ec_EM)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ec_EM
!        case("cld_frac_em")
!           Sds_Info(Var_Idx)%Standard_Name = "cloud_fraction_from_eyre_menzel"
!           Sds_Info(Var_Idx)%Actual_Range = [0.0,1.2]
!           if (allocated(N_EM)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => N_EM
!        case("cld_residual_em")
!           Sds_Info(Var_Idx)%Standard_Name = "residual_from_eyre_menzel"
!           Sds_Info(Var_Idx)%Actual_Range = [0.0,100.0]
!           if (allocated(Res_EM)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Res_EM
         case("cld_temp_cirrus_background")
            Sds_Info(Var_Idx)%Standard_Name = "temperature_at_cloud_top_for_cirrus_from_background"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(Tc_Cirrus_Background)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tc_Cirrus_Background
         case("cld_height_cirrus_background")
            Sds_Info(Var_Idx)%Standard_Name = "height_at_cloud_top_for_cirrus_from_background"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(Zc_Cirrus_Background)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Zc_Cirrus_Background
         case("cld_height_co2")
            Sds_Info(Var_Idx)%Standard_Name = "height_at_cloud_top_from_co2_slicing"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(Zc_Co2)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Zc_Co2
         case("cld_emiss_co2")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_emissivity_from_co2_slicing"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.2]
            if (allocated(Ec_Co2)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ec_Co2
         case("cld_height_lidar")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_height_from_lidar"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            if (allocated(Caliop_Cld_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Caliop_Cld_Height
         case("freezing_altitude")
            Sds_Info(Var_Idx)%Standard_Name = "freezing_altitude"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100000.0]
            Sds_Info(Var_Idx)%Long_Name = "altitude of freezing level"
            Sds_Info(Var_Idx)%Units = "feet"
            if (allocated(NWP_PIX%FrzAlt)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%FrzAlt
         case("homogenous_freezing_altitude")
            Sds_Info(Var_Idx)%Standard_Name = "homogenous_freezing_altitude"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100000.0]
            Sds_Info(Var_Idx)%Long_Name = "altitude of freezing level"
            Sds_Info(Var_Idx)%Units = "feet"
            if (allocated(NWP_PIX%HomoFrzAlt)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%HomoFrzAlt
         case("cape_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "cape_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,10000.0]
            Sds_Info(Var_Idx)%Long_Name = "Convective Available Potential Energy from NWP"
            Sds_Info(Var_Idx)%Units = "Jules/kg"
            if (allocated(NWP_PIX%CAPE)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%CAPE
         !----------------------------------------------------------------------------------------------------
         ! Base Members
         !----------------------------------------------------------------------------------------------------
         case("cld_height_base")
            Sds_Info(Var_Idx)%Standard_Name = "base_height_of_cloud"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Base Height from BASE"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(BASE%Zc_Base)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => BASE%Zc_Base
         case("cld_base_altitude")
            Sds_Info(Var_Idx)%Standard_Name = "altitude_at_cloud_base"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100000.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Base Altitude from BASE"
            Sds_Info(Var_Idx)%Units = "feet"
            if (allocated(BASE%Base_Alt)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => BASE%Base_Alt
         case("lower_cld_base_altitude")
            Sds_Info(Var_Idx)%Standard_Name = "altitude_at_lower_cloud_base"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100000.0]
            Sds_Info(Var_Idx)%Long_Name = "Lower Cloud Base Altitude from BASE"
            Sds_Info(Var_Idx)%Units = "feet"
            if (allocated(BASE%Lower_Base_Alt)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => BASE%Lower_Base_Alt
         case("cld_geo_thick")
            Sds_Info(Var_Idx)%Standard_Name = "cld_geo_thick"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,10000.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud geometrical thickness"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(BASE%Geo_Thickness)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => BASE%Geo_Thickness
         case("cld_temp_base")
            Sds_Info(Var_Idx)%Standard_Name = "base_temp_of_cloud"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Base Temperature from BASE"
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(BASE%Tc_Base)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => BASE%Tc_Base

         !----------------------------------------------------------------------------------------------------
         ! CCL Members
         !----------------------------------------------------------------------------------------------------
         case("cloud_fraction")
            Sds_Info(Var_Idx)%Standard_Name = "total_cloud_area_fraction"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "CCL total cloud fraction"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Cloud_Fraction)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CCL%Cloud_Fraction
         case("cloud_fraction_uncertainty")
            Sds_Info(Var_Idx)%Standard_Name = "total_cloud_fraction_uncertainty"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "CCL uncertainty in the total cloud fraction"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Cloud_Fraction_Uncer)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CCL%Cloud_Fraction_Uncer
         case("ccl_layer_flag")
            Sds_Info(Var_Idx)%Standard_Name = "ccl_layer_flag"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Cloud_Layer)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Cloud_Layer
         case("ccl_1")
            Sds_Info(Var_Idx)%Standard_Name = "ccl_1"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Cloud_Fraction_Layer1)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Cloud_Fraction_Layer1
         case("ccl_2")
            Sds_Info(Var_Idx)%Standard_Name = "ccl_2"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Cloud_Fraction_Layer2)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Cloud_Fraction_Layer2
         case("ccl_3")
            Sds_Info(Var_Idx)%Standard_Name = "ccl_3"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Cloud_Fraction_Layer3)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Cloud_Fraction_Layer3
         case("ccl_4")
            Sds_Info(Var_Idx)%Standard_Name = "ccl_4"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Cloud_Fraction_Layer4)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Cloud_Fraction_Layer4
         case("ccl_5")
            Sds_Info(Var_Idx)%Standard_Name = "ccl_5"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Cloud_Fraction_Layer5)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Cloud_Fraction_Layer5
         case("supercooled_cloud_fraction")
            Sds_Info(Var_Idx)%Standard_Name = "total_supercooled_cloud_area_fraction"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "CCL total supercooled cloud fraction"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Supercooled_Cloud_Fraction)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CCL%Supercooled_Cloud_Fraction
         case("supercooled_ccl_layer_flag")
            Sds_Info(Var_Idx)%Standard_Name = "supercooled_ccl_layer_flag"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Supercooled_Cloud_Layer)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Supercooled_Cloud_Layer
         case("supercooled_ccl_1")
            Sds_Info(Var_Idx)%Standard_Name = "supercooled_cloud_fraction_layer_1"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Supercooled_Cloud_Fraction_Layer1)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Supercooled_Cloud_Fraction_Layer1
         case("supercooled_ccl_2")
            Sds_Info(Var_Idx)%Standard_Name = "supercooled_cloud_fraction_layer_2"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Supercooled_Cloud_Fraction_Layer2)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Supercooled_Cloud_Fraction_Layer2
         case("supercooled_ccl_3")
            Sds_Info(Var_Idx)%Standard_Name = "supercooled_cloud_fraction_layer_3"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Supercooled_Cloud_Fraction_Layer3)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Supercooled_Cloud_Fraction_Layer3
         case("supercooled_ccl_4")
            Sds_Info(Var_Idx)%Standard_Name = "supercooled_cloud_fraction_layer_4"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Supercooled_Cloud_Fraction_Layer4)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Supercooled_Cloud_Fraction_Layer4
         case("supercooled_ccl_5")
            Sds_Info(Var_Idx)%Standard_Name = "supercooled_cloud_fraction_layer_5"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Supercooled_Cloud_Fraction_Layer5)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Supercooled_Cloud_Fraction_Layer5
         case("conv_cloud_fraction")
            Sds_Info(Var_Idx)%Standard_Name = "total_convective_cloud_area_fraction"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "CCL total convective cloud fraction"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Conv_Cloud_Fraction)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => CCL%Conv_Cloud_Fraction
         case("conv_ccl_layer_flag")
            Sds_Info(Var_Idx)%Standard_Name = "conv_ccl_layer_flag"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Conv_Cloud_Layer)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Conv_Cloud_Layer
         case("conv_ccl_1")
            Sds_Info(Var_Idx)%Standard_Name = "convective_cloud_fraction_layer_1"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Conv_Cloud_Fraction_Layer1)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Conv_Cloud_Fraction_Layer1
         case("conv_ccl_2")
            Sds_Info(Var_Idx)%Standard_Name = "convective_cloud_fraction_layer_2"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Conv_Cloud_Fraction_Layer2)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Conv_Cloud_Fraction_Layer2
         case("conv_ccl_3")
            Sds_Info(Var_Idx)%Standard_Name = "convective_cloud_fraction_layer_3"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Conv_Cloud_Fraction_Layer3)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Conv_Cloud_Fraction_Layer3
         case("conv_ccl_4")
            Sds_Info(Var_Idx)%Standard_Name = "convective_cloud_fraction_layer_4"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Conv_Cloud_Fraction_Layer4)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Conv_Cloud_Fraction_Layer4
         case("conv_ccl_5")
            Sds_Info(Var_Idx)%Standard_Name = "convective_cloud_fraction_layer_5"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(CCL%Conv_Cloud_Fraction_Layer5)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => CCL%Conv_Cloud_Fraction_Layer5

         case("iwc_dcomp")
            Sds_Info(Var_Idx)%Standard_Name = "iwc_dcomp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud ice content from DCOMP"
            Sds_Info(Var_Idx)%Units = "g m-3"
            if (allocated(Iwc_Dcomp)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Iwc_Dcomp

         case("lwc_dcomp")
            Sds_Info(Var_Idx)%Standard_Name = "lwc_dcomp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud water content from DCOMP"
            Sds_Info(Var_Idx)%Units = "g m-3"
            if (allocated(lwc_Dcomp)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Lwc_Dcomp

         !---------------------------------------------------------------------
         ! DCOMP
         !---------------------------------------------------------------------
         case("cld_opd_dcomp")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_optical_thickness_due_to_cloud"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,160.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud optical depth at the nominal wavelength of 0.65 microns, "//&
                                          "determined from DCOMP"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Tau_DCOMP)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tau_DCOMP
         case("cld_opd_dcomp_1")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_optical_thickness_due_to_cloud"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,160.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud optical depth at the nominal wavelength of 0.65 microns, "//&
                                          "determined from DCOMP Mode 1"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Tau_DCOMP_1)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tau_DCOMP_1
         case("cld_opd_dcomp_2")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_optical_thickness_due_to_cloud"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,160.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud optical depth at the nominal wavelength of 0.65 microns, "//&
                                          "determined from DCOMP Mode 2"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Tau_DCOMP_2)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tau_DCOMP_2
         case("cld_opd_dcomp_3")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_optical_thickness_due_to_cloud"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,160.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud optical depth at the nominal wavelength of 0.65 microns, "//&
                                          "determined from DCOMP Mode 3"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Tau_DCOMP_3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tau_DCOMP_3
         case("cld_opd_aux")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_optical_thickness_due_to_cloud_aux"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,160.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud optical depth at the nominal wavelength of 0.65 microns, "//&
                               "read from aux file"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Tau_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tau_Aux
         case("cld_cwp_dcomp")
            Sds_Info(Var_Idx)%Standard_Name = "cld_cwp_dcomp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1200.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud top water path from DCOMP"
            Sds_Info(Var_Idx)%Units = "g m-2"
            if (allocated(Cwp_Dcomp)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Cwp_Dcomp
         case("cld_reff_dcomp")
            Sds_Info(Var_Idx)%Standard_Name = "effective_radius_of_cloud_condensed_water_particles_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,160.0]
            Sds_Info(Var_Idx)%Long_Name = "effective radius of cloud particles determined from DCOMP; "//&
                                          "see attributes for channels used"
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(Reff_DCOMP)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Reff_DCOMP
         case("cld_reff_dcomp_1")
            Sds_Info(Var_Idx)%Standard_Name = "effective_radius_of_cloud_condensed_water_particles_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,160.0]
            Sds_Info(Var_Idx)%Long_Name = "effective radius of cloud particles determined from DCOMP Mode 1; "//&
                              "see attributes for channels used"
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(Reff_DCOMP_1)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Reff_DCOMP_1
         case("cld_reff_dcomp_2")
            Sds_Info(Var_Idx)%Standard_Name = "effective_radius_of_cloud_condensed_water_particles_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,160.0]
            Sds_Info(Var_Idx)%Long_Name = "effective radius of cloud particles determined from DCOMP Mode 2; "//&
                              "see attributes for channels used"
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(Reff_DCOMP_2)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Reff_DCOMP_2
         case("cld_reff_dcomp_3")
            Sds_Info(Var_Idx)%Standard_Name = "effective_radius_of_cloud_condensed_water_particles_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,160.0]
            Sds_Info(Var_Idx)%Long_Name = "effective radius of cloud particles determined from DCOMP Mode 3; "//&
                              "see attributes for channels used"
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(Reff_DCOMP_3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Reff_DCOMP_3
         case("cld_reff_dcomp_fit")
            Sds_Info(Var_Idx)%Standard_Name = "effective_radius_of_cloud_condensed_water_particles_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,160.0]
            Sds_Info(Var_Idx)%Long_Name = "effective radius of cloud particles determined from DCOMP after an amr2-derived"//&
                              " fit ; see attributes for channels used - NOAA CDR"
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(Reff_Dcomp_Fit)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Reff_Dcomp_Fit
         case("cld_reff_aux")
            Sds_Info(Var_Idx)%Standard_Name = "effective_radius_of_cloud_condensed_water_particles_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,160.0]
            Sds_Info(Var_Idx)%Long_Name = "effective radius of cloud particles read from aux file"
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(Reff_Aux)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Reff_Aux
         case("cld_opd_dcomp_unc")
            Sds_Info(Var_Idx)%Standard_Name = "cld_opd_dcomp_unc"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,160.0]
            Sds_Info(Var_Idx)%Long_Name = "uncertainty in the log10 cloud optical depth at the nominal wavelength "// &
                              "of 0.65 microns, determined from DCOMP; see attributes for channels used"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Tau_Dcomp_Cost)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tau_Dcomp_Cost
         case("cld_reff_dcomp_unc")
            Sds_Info(Var_Idx)%Standard_Name = "cld_reff_dcomp_unc"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,160.0]
            Sds_Info(Var_Idx)%Long_Name = "uncertainty in the log10 effective radius of cloud particle determined from DCOMP; "// &
                               "see attributes for channels used"
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(Reff_Dcomp_Cost)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Reff_Dcomp_Cost
         case("cld_opd_dcomp_qf")
            Sds_Info(Var_Idx)%Standard_Name = "cld_opd_dcomp_qf"
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Actual_Range = [0.0,0.0]
            Sds_Info(Var_Idx)%Long_Name = "quality flag for cloud optical depth from DCOMP "// &
                               "not attempted=0, failed=1, low quality=2, high quality=3"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Tau_Dcomp_Qf)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Tau_Dcomp_Qf
         case("cld_reff_dcomp_qf")
            Sds_Info(Var_Idx)%Standard_Name = "cld_reff_dcomp_qf"
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Actual_Range = [0.0,0.0]
            Sds_Info(Var_Idx)%Long_Name = "quality flag for cloud effective radius from DCOMP "// &
                               "not attempted=0, failed=1, low quality=2, high quality=3"
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(Reff_Dcomp_Qf)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Reff_Dcomp_Qf
         case("insolation_dcomp")
            Sds_Info(Var_Idx)%Standard_Name = "surface_downwelling_shortwave_flux_dcomp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1500.0]
            Sds_Info(Var_Idx)%Long_Name = "surface downwelling shortwave flux computed from the DCOMP cloud properties"
            Sds_Info(Var_Idx)%Units = "W m-2"
            if (allocated(Insolation_Dcomp)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Insolation_Dcomp
         case("insolation_diffuse_dcomp")
            Sds_Info(Var_Idx)%Standard_Name = "surface_downwelling_shortwave_flux_diffuse_dcomp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1500.0]
            Sds_Info(Var_Idx)%Long_Name = "diffuse component of the surface downwelling shortwave flux "// &
                              "computed from the DCOMP cloud properties"
            Sds_Info(Var_Idx)%Units = "W m-2"
            if (allocated(Insolation_Diffuse_Dcomp)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Insolation_Diffuse_Dcomp
         case("cdnc_dcomp")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_droplet_number_concentration"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1000.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud_droplet_number_concentration from DCOMP algorithm"
            Sds_Info(Var_Idx)%Units = "cm-3"
            if (allocated(Cdnc_Dcomp)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Cdnc_Dcomp
         case("hcld_dcomp")
            Sds_Info(Var_Idx)%Standard_Name = "geometrical_thickness_of_liquid_clouds"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,4000.0]
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(Hcld_Dcomp)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Hcld_Dcomp
         case("dcomp_quality")
            Sds_Info(Var_Idx)%Standard_Name =  "dcomp_quality_flags_packed"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Long_Name = "quality flags for DCOMP products "// &
                               "1:Processed (0=no,1=yes) "// &
                               "2:valid COD retrieval (0=yes,1=no) "// &
                               "3:valid REF retrieval (0=yes,1=no) "// &
                               "4:degraded COD retrieval (0=no,1=degraded) "// &
                               "5:degraded REF retrieval (0=no,1=degraded) "// &
                               "6:convergency (0=no,1=yes) "// &
                               "7:glint (0=no,1=yes)"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(DCOMP_Quality_Flag)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => DCOMP_Quality_Flag

         case("dcomp_info")
            Sds_Info(Var_Idx)%Standard_Name = "dcomp_information_flags_packed"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT16
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_SHORT
            Sds_Info(Var_Idx)%Long_Name = "processing flags for DCOMP "// &
                               "1:info flag set ? (0=no,1=yes) "// &
                               "2:land/sea mask (0=land,1=sea) "// &
                               "3:day/night mask (0=Day,1=Night) "// &
                               "4:twilight (65-82 solar zenith) (0=no,1=yes) "// &
                               "5:snow (0=no,1= snow) "// &
                               "6:sea-ice (0=no,1=sea-ice) "// &
                               "7:phase (0=liquid,1=ice) "// &
                               "8:thick_cloud (0=no,1=yes) " // &
                               "9:thin_cloud (0=no,1=yes) "
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(DCOMP_Info_Flag)) Sds_Info(Var_Idx)%Sds_Data_2d_I2 => DCOMP_Info_Flag

         case("cloud_transmission_0_65um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_transmission_0_65um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud transmission at 0.65 microns nominal from DCOMP"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Cloud_063um_Transmission_Solar)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Cloud_063um_Transmission_Solar

         case("cloud_albedo_0_65um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_albedo_0_65um_nom"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud albedo at 0.65 microns nominal from DCOMP"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Cloud_063um_Albedo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Cloud_063um_Albedo

         case("rain_rate")
            Sds_Info(Var_Idx)%Standard_Name = "rain_rate"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,32.0]
            Sds_Info(Var_Idx)%Long_Name = "derived rain rate from DCOMP"
            Sds_Info(Var_Idx)%Units = "mm h-1"
            if (allocated(Rain_Rate_DCOMP)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Rain_Rate_DCOMP

         case("refl_0_65um_nom_asym_dcomp")
            Sds_Info(Var_Idx)%Standard_Name = "refl_0_65um_nom_asym_dcomp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "asymptotic_reflectance_for_DCOMP_solution"
            Sds_Info(Var_Idx)%Units = "%"
            if (allocated(Refl_Asym_DCOMP)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Refl_Asym_DCOMP

         !------------------------------------------------------------------------------------------------
         ! NLCOMP
         !------------------------------------------------------------------------------------------------
         case("cld_opd_nlcomp")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_optical_thickness_due_to_cloud"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,160.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud optical depth at the nominal wavelength of 0.65 microns, "//&
                                          "determined from NLCOMP"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Tau_NLCOMP)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tau_NLCOMP
         case("cld_reff_nlcomp")
            Sds_Info(Var_Idx)%Standard_Name = "effective_radius_of_cloud_condensed_water_particles_at_cloud_top"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,160.0]
            Sds_Info(Var_Idx)%Long_Name = "effective radius of cloud particles determined from NLCOMP; "//&
                                          "see attributes for channels used"
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(Reff_NLCOMP)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Reff_NLCOMP
         case("cld_opd_nlcomp_unc")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_optical_thickness_due_to_cloud_uncertainty"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,160.0]
            Sds_Info(Var_Idx)%Long_Name = "uncertainty in cloud optical depth at the nominal wavelength "// &
                              "of 0.65 microns, determined from NLCOMP"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Tau_Nlcomp_Cost)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tau_Nlcomp_Cost
         case("cld_reff_nlcomp_unc")
            Sds_Info(Var_Idx)%Standard_Name = "effective_radius_of_cloud_particle_uncertainty" 
            Sds_Info(Var_Idx)%Actual_Range = [0.0,160.0]
            Sds_Info(Var_Idx)%Long_Name = "effective radius of cloud particle uncertainty " // &
                               "determined from NLCOMP; see attributes for channels used"
            Sds_Info(Var_Idx)%Units = "micron"
            if (allocated(Reff_Nlcomp_Cost)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Reff_Nlcomp_Cost
         case("nlcomp_quality")
            Sds_Info(Var_Idx)%Standard_Name =  "nlcomp_quality_flags_packed"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Long_Name =  "nlcomp_quality_flags_packed"
            if (allocated(NLCOMP_Quality_Flag)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => NLCOMP_Quality_Flag
            Sds_Info(Var_Idx)%Flags_String = "quality flags for NLCOMP products "// &
                               " see documentation http://cimss.ssec.wisc.edu/clavr/ "// &
                               "1:Processed (0=no,1=yes) "// &
                               "2:valid COD retrieval (0=yes,1=no) "// &
                               "3:valid REF retrieval (0=yes,1=no) "// &
                               "4:degraded COD retrieval (0=no,1=degraded) "// &
                               "5:degraded REF retrieval (0=no,1=degraded) "// &
                               "6:convergency (0=no,1=yes) "// &
                               "7:glint (0=no,1=yes) "
            Sds_Info(Var_Idx)%Units = "none"
         case("nlcomp_info")
            Sds_Info(Var_Idx)%Standard_Name =  "nlcomp_information_flags_packed"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT16
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_SHORT
            Sds_Info(Var_Idx)%Long_Name =  "nlcomp_information_flags_packed"
            Sds_Info(Var_Idx)%Flags_String =  "processing flags for NLCOMP "// &
                               "see http://cimss.ssec.wisc.edu/clavr/ "// &
                               "1: info flag set ? (0=no,1=yes) "// &
                               "2: land/sea mask (0=land,1=sea) "// &
                               "3: day/night mask (0=Day,1=Night) "// &
                               "4: twilight (65-82 solar zenith) (0=no,1=yes) "// &
                               "5: snow (0=no,1= snow) "// &
                               "6: sea-ice (0=no,1=sea-ice) "// &
                               "7: phase (0=liquid,1=ice) "// &
                               "8: thick_cloud (0=no,1=yes) " // &
                               "9: thin_cloud (0=no,1=yes) "
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(NLCOMP_Info_Flag)) Sds_Info(Var_Idx)%Sds_Data_2d_I2 => NLCOMP_Info_Flag

         !---------------------------------------------
         !  DARK SKY COMPOSITE
         !---------------------------------------------
         case ("refl_0_65um_nom_dark_static")
            Sds_info(Var_Idx)%Standard_Name = "dark_sky_composite_reflectance_0.65_micron"
            Sds_Info(Var_Idx)%Units =  "%"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,120.0]
            if (allocated(Static_Ref_065um_Dark_Composite)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Static_Ref_065um_Dark_Composite
         case ("refl_0_65um_nom_dark_static_stddev")
            Sds_info(Var_Idx)%Standard_Name = "dark_sky_composite_reflectance_0.65_micron_std"
            Sds_Info(Var_Idx)%Units =  "%"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,50.0]
            if (allocated(Static_Ref_065um_Dark_Composite_Stddev)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Static_Ref_065um_Dark_Composite_Stddev
         case("subpixel_cloud_fraction")
            Sds_Info(Var_Idx)%Standard_Name = "subpixel_total_cloud_area_fraction"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "subpixel total cloud fraction"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Subpixel_Cloud_Fraction)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Subpixel_Cloud_Fraction


         !------------------------------------------------------------------------------------------------
         ! NUCAPS
         !------------------------------------------------------------------------------------------------
         case("cld_press_nucaps_1")
            Sds_Info(Var_Idx)%Standard_Name = "cloud-top_pressure_from_NUCAPS_for_Layer_1"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud-top pressure from NUCAPS for Layer 1"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(NUCAPS%Cld_Press_Layer1)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NUCAPS%Cld_Press_Layer1
         case("cld_press_nucaps_2")
            Sds_Info(Var_Idx)%Standard_Name = "cloud-top_pressure_from_NUCAPS_for_Layer_2"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud-top pressure from NUCAPS for Layer 2"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(NUCAPS%Cld_Press_Layer2)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NUCAPS%Cld_Press_Layer2
         case("cld_frac_nucaps_1")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_fraction_from_NUCAPS_for_Layer_1"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud fraction from NUCAPS for Layer 1"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(NUCAPS%Cld_Fraction_Layer1)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NUCAPS%Cld_Fraction_Layer1
         case("cld_frac_nucaps_2")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_fraction_from_NUCAPS_for_Layer_2"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud fraction from NUCAPS for Layer 2"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(NUCAPS%Cld_Fraction_Layer2)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NUCAPS%Cld_Fraction_Layer2
         case("cld_frac_nucaps")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_fraction_from_NUCAPS"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud fraction from NUCAPS"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(NUCAPS%Cld_Fraction)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NUCAPS%Cld_Fraction
         case("cld_press_nucaps")
            Sds_Info(Var_Idx)%Standard_Name = "cloud-top_pressure_from_NUCAPS"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud-top pressure from NUCAPS"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(NUCAPS%Cld_Press)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NUCAPS%Cld_Press
         case("cld_temp_nucaps")
            Sds_Info(Var_Idx)%Standard_Name = "cloud-top_temperature_from_NUCAPS"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud-top temperature from NUCAPS"
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(NUCAPS%Cld_Temp)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NUCAPS%Cld_Temp
         case("cld_height_nucaps")
            Sds_Info(Var_Idx)%Standard_Name = "cloud-top_height_from_NUCAPS"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20000.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud-top height from NUCAPS"
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(NUCAPS%Cld_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NUCAPS%Cld_Height
         case("cld_temp_nucaps_smoothed")
            Sds_Info(Var_Idx)%Standard_Name = "cloud-top_temperature_from_NUCAPS_smoothed"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Long_Name = "cloud-top temperature from NUCAPS smoothed"
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(NUCAPS%Cld_Temp_Smoothed)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NUCAPS%Cld_Temp_Smoothed

         !------------------------------------------------------------------------------------------------
         ! Other
         !------------------------------------------------------------------------------------------------
         case("beta_11um_67um_tropopause")
            Sds_Info(Var_Idx)%Standard_Name = "beta_11um_67um_tropopause"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,2.0]
            Sds_Info(Var_Idx)%Long_Name = "beta ratio from 11um and 67um referenced to tropopause"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Beta_11um_67um_Tropo_Rtm)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Beta_11um_67um_Tropo_Rtm
         case("beta_11um_85um_tropopause")
            Sds_Info(Var_Idx)%Standard_Name = "beta_11um_85um_tropopause"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,2.0]
            Sds_Info(Var_Idx)%Long_Name = "beta ratio from 11um and 85um referenced to tropopause"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Beta_11um_85um_Tropo_Rtm)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Beta_11um_85um_Tropo_Rtm 
         case("beta_11um_104um_tropopause")
            Sds_Info(Var_Idx)%Standard_Name = "beta_11um_104um_tropopause"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,2.0]
            Sds_Info(Var_Idx)%Long_Name = "beta ratio from 11um and 10.4um referenced to tropopause"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Beta_11um_104um_Tropo_Rtm)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Beta_11um_104um_Tropo_Rtm 
         case("beta_11um_12um_tropopause")
            Sds_Info(Var_Idx)%Standard_Name = "beta_11um_12um_tropopause"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,2.0]
            Sds_Info(Var_Idx)%Long_Name = "beta ratio from 11um and 12um referenced to tropopause"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Beta_11um_12um_Tropo_Rtm)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Beta_11um_12um_Tropo_Rtm 
         case("beta_11um_133um_tropopause")
            Sds_Info(Var_Idx)%Standard_Name = "beta_11um_133um_tropopause"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,2.0]
            Sds_Info(Var_Idx)%Long_Name = "beta ratio from 11um and 133um referenced to tropopause"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Beta_11um_133um_Tropo_Rtm)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Beta_11um_133um_Tropo_Rtm

         !---------------------------------------------------------------------
         ! NWP
         !---------------------------------------------------------------------
         case("wind_direction_cloud_top_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "wind_direction_cloud_top_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,360.0]
            Sds_Info(Var_Idx)%Long_Name = "Wind Direction at Cloud-Top from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "degrees"
            if (allocated(NWP_PIX%Wnd_Dir_Cld_Top)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Wnd_Dir_Cld_Top
         case("wind_speed_cloud_top_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "wind_speed_cloud_top_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100.0]
            Sds_Info(Var_Idx)%Long_Name = "Wind Speed at Cloud-Top from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "m.s-1"
            if (allocated(NWP_PIX%Wnd_Spd_Cld_Top)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Wnd_Spd_Cld_Top
         case("surface_air_temperature_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "surface_air_temperature_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [220.0,340.0]
            Sds_Info(Var_Idx)%Long_Name = "Surface Air Temperature from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(NWP_PIX%Tair)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Tair
         case("surface_temperature_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "surface_temperature_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [220.0,340.0]
            Sds_Info(Var_Idx)%Long_Name = "Surface Temperature from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(NWP_PIX%Tsfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Tsfc
         case("cld_press_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "cld_press_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1100.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Pressure from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(NWP_PIX%Pc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Pc
         case("cld_temp_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "cld_temp_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Top Temperature from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(NWP_PIX%Tc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Tc
         case("cloud_fraction_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_fraction_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Fraction Temperature from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(NWP_PIX%Cfrac)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Cfrac
         case("cld_cwp_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "cld_cwp_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1200.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Water Path from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "g m-2"
         case("cld_iwp_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "cld_iwp_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1200.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Ice Path from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "g m-2"
            if (allocated(NWP_PIX%Iwp)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Iwp
         case("cld_lwp_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "cld_lwp_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1200.0]
            Sds_Info(Var_Idx)%Long_Name = "Cloud Liquid Water Path from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "g m-2"
            if (allocated(NWP_PIX%Lwp)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Lwp
         case("surface_relative_humidity_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "surface_relative_humidity_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,110.0]
            Sds_Info(Var_Idx)%Long_Name = "Surface Relative Humidity from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "%"
            if (allocated(NWP_PIX%Rhsfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Rhsfc
         case("rh300_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "300hpa_relative_humidity_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,110.0]
            Sds_Info(Var_Idx)%Long_Name = "300hpa Relative Humidity from NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "%"
            if (allocated(NWP_PIX%Rh300)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Rh300
         case("uth_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "upper_tropospheric_humidity_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,110.0]
            Sds_Info(Var_Idx)%Long_Name = "Upper Tropospheric Humidity NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "%"
            if (allocated(NWP_PIX%Uth)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Uth
         case("div_sfc_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "divergence_sfc_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [-50.0,50.0]
            Sds_Info(Var_Idx)%Long_Name = "Surface Level Divergence NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "1.0e06/s"
            if (allocated(NWP_PIX%Div_Sfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Div_Sfc
         case("div_200_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "divergence_200hPa_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [-50.0,50.0]
            Sds_Info(Var_Idx)%Long_Name = "200hPa Level Divergence NWP Ancillary Data"
            Sds_Info(Var_Idx)%Units = "1.0e06/s"
            if (allocated(NWP_PIX%Div_200)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Div_200
         case("surface_pressure_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "surface_pressure_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [700.0,1100.0]
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(NWP_PIX%Psfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Psfc
         case("ndvi_sfc_white_sky_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "ndvi_sfc_white_sky_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [-0.5,1.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Ndvi_Sfc_White_Sky)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ndvi_Sfc_White_Sky
         case("ndvi_sfc")
            Sds_Info(Var_Idx)%Standard_Name = "ndvi_sfc"
            Sds_Info(Var_Idx)%Actual_Range = [-0.5,1.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Ndvi_Sfc)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ndvi_Sfc
         case("ndvi_toa")
            Sds_Info(Var_Idx)%Standard_Name = "ndvi_toa"
            Sds_Info(Var_Idx)%Long_Name = "top_of_atmosphere_ndvi"
            Sds_Info(Var_Idx)%Actual_Range = [-0.5,1.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Ndvi_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ndvi_Toa
         case("ndsi_toa")
            Sds_Info(Var_Idx)%Standard_Name = "ndsi_toa"
            Sds_Info(Var_Idx)%Long_Name = "top_of_atmosphere_ndsi"
            Sds_Info(Var_Idx)%Actual_Range = [-0.5,1.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Ndsi_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ndsi_Toa
         case("nddi_toa")
            Sds_Info(Var_Idx)%Standard_Name = "nddi_toa"
            Sds_Info(Var_Idx)%Long_Name = "top_of_atmosphere_nddi"
            Sds_Info(Var_Idx)%Actual_Range = [-0.5,1.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Nddi_Toa)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Nddi_Toa
         case("total_precipitable_water_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_mass_content_of_water_vapor"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,10.0]
            Sds_Info(Var_Idx)%Units = "cm"
            if (allocated(NWP_PIX%Tpw)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Tpw
         case("total_column_ozone_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "atmosphere_mass_content_of_ozone"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1000.0]
            Sds_Info(Var_Idx)%Units = "dobson"
            if (allocated(NWP_PIX%Ozone)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Ozone
         case("tropopause_temperature_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "tropopause_temperature_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,260.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(NWP_PIX%Ttropo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Ttropo
         case("tropopause_pressure_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "tropopause_pressure_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [50.0,500.0]
            Sds_Info(Var_Idx)%Units = "hPa"
            if (allocated(NWP_PIX%Ptropo)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%Ptropo
         case("k_index_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "k_index_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,100.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(NWP_PIX%K_Index)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%K_Index
         case("lcl_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "lcl_nwp"
            Sds_Info(Var_Idx)%Long_Name = "lifting_condensation_level_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,10000.0]
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(NWP_PIX%LCL_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%LCL_Height
         case("ccl_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "ccl_nwp"
            Sds_Info(Var_Idx)%Long_Name = "convective_condensation_level_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,10000.0]
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(NWP_PIX%CCL_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%CCL_Height
         case("lfc_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "lfc_nwp"
            Sds_Info(Var_Idx)%Long_Name = "level_of_free_convection_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,10000.0]
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(NWP_PIX%LFC_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%LFC_Height
         case("el_nwp")
            Sds_Info(Var_Idx)%Standard_Name = "el_nwp"
            Sds_Info(Var_Idx)%Long_Name = "equilibrium_level_nwp"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,10000.0]
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(NWP_PIX%EL_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => NWP_PIX%EL_Height
         !---------------------------------------------------------------------
         ! Surface Temperature
         !---------------------------------------------------------------------
         case("sea_surface_temperature_retrieved")
            Sds_Info(Var_Idx)%Standard_Name = "sea_surface_temperature_retrieved"
            Sds_Info(Var_Idx)%Actual_Range = [265.0,315.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(Sst_Retrieved)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Sst_Retrieved
         case("surface_temperature_retrieved")
            Sds_Info(Var_Idx)%Standard_Name = "surface_temperature_retrieved"
            Sds_Info(Var_Idx)%Actual_Range = [220.0,340.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(Tsfc_Retrieved)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tsfc_Retrieved
         case("oisst")
            Sds_Info(Var_Idx)%Standard_Name = "oisst"
            Sds_Info(Var_Idx)%Actual_Range = [265.0,315.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(Sst_Anal)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Sst_Anal

         !---------------------------------------------------------------------
         ! Needed for Tuning
         !---------------------------------------------------------------------
         case("emiss_3_75um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "top_of_atmosphere_emissivity_3_75_micron_nominal"
            Sds_Info(Var_Idx)%Actual_Range = [0.5,3.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Ch(20)%Emiss_Rel_11um)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(20)%Emiss_Rel_11um
         case("emiss_3_75um_nom_rel_10_4um_nom")
            Sds_Info(Var_Idx)%Standard_Name = "top_of_atmosphere_emissivity_3_75_micron_nominal_relative_to_10_4_micron_nominal"
            Sds_Info(Var_Idx)%Actual_Range = [0.5,3.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Ch(20)%Emiss_Rel_10_4um)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(20)%Emiss_Rel_10_4um
         case("emiss_3_75um_nom_clear")
            Sds_Info(Var_Idx)%Standard_Name = "top_of_atmosphere_emissivity_3_75_micron_nominal_rel_11_micron_nominal_clear"
            Sds_Info(Var_Idx)%Actual_Range = [0.5,3.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Ch(20)%Emiss_Rel_11um_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(20)%Emiss_Rel_11um_Clear
         case("emiss_3_75um_nom_rel_10_4um_nom_clear")
            Sds_Info(Var_Idx)%Standard_Name = "top_of_atmosphere_emissivity_3_75_micron_nominal_rel_10_4_micron_nominal_clear"
            Sds_Info(Var_Idx)%Actual_Range = [0.5,3.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Ch(20)%Emiss_Rel_10_4um_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(20)%Emiss_Rel_10_4um_Clear
         case("emiss_3_75um_nom_median_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "emissivity_3_75_micron_nominal_median_3x3"
            Sds_Info(Var_Idx)%Actual_Range = [0.5,3.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Ems_Ch20_Median_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ems_Ch20_Median_3x3
         case("temp_11um_vs_67um_covar_5x5")
            Sds_Info(Var_Idx)%Standard_Name = "brightness_temperature_11_vs_67_micron_5x5_covariance"
            Sds_Info(Var_Idx)%Actual_Range = [-10.0,10.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(Covar_Ch27_Ch31_5x5)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Covar_Ch27_Ch31_5x5
         case("diff_ch31_ch32_bt_ch31_max_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "difference_11_minus_12_brightness_temperature_max_3x3"
            Sds_Info(Var_Idx)%Actual_Range = [-4.0,20.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(Btd_Ch31_Ch32_Bt_Ch31_Max_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Btd_Ch31_Ch32_Bt_Ch31_Max_3x3
         case("cld_temp_opaque")
            Sds_Info(Var_Idx)%Standard_Name = "cloud_top_temperature_opaque"
            Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(Tc_Opaque_Cloud)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tc_Opaque_Cloud
!        case("cld_temp_em")
!           Sds_Info(Var_Idx)%Standard_Name = "cloud_top_temperature_eyre_menzel"
!           Sds_Info(Var_Idx)%Actual_Range = [160.0,320.0]
!           Sds_Info(Var_Idx)%Units = "K"
!           if (allocated(Tc_EM)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Tc_EM
         case("temp_3_75um_nom_median_5x5")
            Sds_Info(Var_Idx)%Standard_Name = "brightness_temperature_3.7_micron_nominal_median_5x5"
            Sds_Info(Var_Idx)%Actual_Range = [180.0,340.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(Bt_Ch20_Median_5x5)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Bt_Ch20_Median_5x5
         case("sst_background_uni_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "sea_surface_skin_temperature_background_uni_3x3"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(Sst_Anal_Uni)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Sst_Anal_Uni
         case("refl_0_65um_nom_clear_sky_min_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_assuming_clear_sky_0_65_micron_nominal_min_3x3"
            Sds_Info(Var_Idx)%Actual_Range = [-2.0,120.0]
            Sds_Info(Var_Idx)%Units = "%"
            if (allocated(ch(1)%Ref_Toa_Clear_Min_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(1)%Ref_Toa_Clear_Min_3x3
         case("refl_0_65um_nom_clear_sky_max_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_assuming_clear_sky_0_65_micron_nominal_max_3x3"
            Sds_Info(Var_Idx)%Actual_Range = [-2.0,120.0]
            Sds_Info(Var_Idx)%Units = "%"
            if (allocated(ch(1)%Ref_Toa_Clear_Max_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => ch(1)%Ref_Toa_Clear_Max_3x3
         case("refl_0_65um_nom_clear_sky_std_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "refl_0_65um_nom_clear_sky_std_3x3"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,20.0]
            Sds_Info(Var_Idx)%Units = "%"
            if (allocated(Ch(1)%Ref_Toa_Clear_Std_3x3)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Ch(1)%Ref_Toa_Clear_Std_3x3
         case("refl_0_86um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_assuming_clear_sky_0_86_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(2)%Ref_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 =>  Ch(2)%Ref_Toa_Clear
         case("refl_1_38um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_assuming_clear_sky_1_38_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(26)%Ref_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 =>  Ch(26)%Ref_Toa_Clear
         case("refl_1_60um_nom_clear_sky")
            Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_assuming_clear_sky_1_60_micron_nominal"
            Sds_Info(Var_Idx)%Units =  "%"
            if (allocated(Ch(6)%Ref_Toa_Clear)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 =>  Ch(6)%Ref_Toa_Clear
         case("btd_11um_12um_nom_nwc")
            Sds_Info(Var_Idx)%Standard_Name = "brightness_temperature_difference_11um_12um_nom_neighboring_warm_center"
            Sds_Info(Var_Idx)%Actual_Range = [-100.0, 100.0]
            Sds_Info(Var_Idx)%Units = "K"
            if (allocated(BTD_11_12um_NWC)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => BTD_11_12um_NWC
         case("emiss_3_75um_nom_nwc")
            Sds_Info(Var_Idx)%Standard_Name = "emissivity_3_75um_nom_neighboring_warm_center"
            Sds_Info(Var_Idx)%Actual_Range = [0.75,1.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Emiss_39um_NWC)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Emiss_39um_NWC
         case("emissivity_11um_tropopause_lrc")
            Sds_Info(Var_Idx)%Standard_Name = "emissivity_11um_tropopause_local_radiative_center"
            Sds_Info(Var_Idx)%Actual_Range = [-0.5,1.2]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Emiss_Tropo_11um_LRC)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Emiss_Tropo_11um_LRC
         case("surface_elevation_std_3x3")
            Sds_Info(Var_Idx)%Standard_Name = "surface_elevation_std_3x3"
            Sds_Info(Var_Idx)%Actual_Range = [0.0,5000.0]
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(Sfc%Zsfc_Std)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Sfc%Zsfc_Std


         !---------------------------------------------------------------------
         ! CALIOP data
         !---------------------------------------------------------------------
         case("caliop_num_cld_layers")
            Sds_Info(Var_Idx)%Actual_Range = [0.0,10.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Caliop_Num_Cld_Layers)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Caliop_Num_Cld_Layers
         case("caliop_cod")
            Sds_Info(Var_Idx)%Actual_Range = [0.0,60.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Caliop_Cod)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Caliop_Cod
         case("caliop_cld_height")
            Sds_Info(Var_Idx)%Actual_Range = [-300.0,20000.0]
            Sds_Info(Var_Idx)%Units = "m"
            if (allocated(Caliop_Cld_Height)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Caliop_Cld_Height


         !---------------------------------------------
         !  LEGACY AEROSOL
         !---------------------------------------------
         case ("aot_0_65um_nom")
            Sds_info(Var_Idx)%Standard_Name = "optical_thickness_of_atmosphere_layer_due_to_aerosol_0.65_micron"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,5.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Aot1)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Aot1
         case ("aot_0_86um_nom")
            Sds_info(Var_Idx)%Standard_Name = "optical_thickness_of_atmosphere_layer_due_to_aerosol_0.86_micron"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,5.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Aot2)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Aot2
         case ("aot_1_60um_nom")
            Sds_info(Var_Idx)%Standard_Name = "optical_thickness_of_atmosphere_layer_due_to_aerosol_1.60_micron"
            Sds_Info(Var_Idx)%Actual_Range = [-0.2,5.0]
            if (allocated(Aot3a)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Aot3a
         case("aot_qf")
            Sds_Info(Var_Idx)%Standard_Name = "optical_thickness_of_atmosphere_layer_quality_flag"
            Sds_Info(Var_Idx)%Long_Name = "optical_thickness_of_atmosphere_layer_quality_flag"
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            if (allocated(Aot_Qf)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Aot_Qf

            
         ! --------------------------------------------
         !  MURI
         ! --------------------------------
         
         case ("aerosol_optical_depth_muri")
            Sds_info(Var_Idx)%Standard_Name = "aerosol_optical_depth_muri" 
            Sds_Info(Var_Idx)%Actual_Range = [0.0,5.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Muri%aod)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Muri%Aod
           
         case ("aerosol_coarse_mode_muri")
            Sds_info(Var_Idx)%Standard_Name = "aerosol_coarse_mode_muri" 
            
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Number_Of_Flags = 5
            Sds_Info(Var_Idx)%Flags_String = "Coarse Mode 1 "// &
                               " Coarse Mode 2 "// &
                               " Coarse Mode 3 "// &
                               " Coarse Mode 4 " // &
                               " Coarse Mode 5 "
                               
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Muri%cm_mode)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Muri%cm_mode
         
         case ("aerosol_fine_mode_muri")
            Sds_info(Var_Idx)%Standard_Name = "aerosol_fine_mode_muri" 
            
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Number_Of_Flags = 4
            Sds_Info(Var_Idx)%Flags_String = "Fine Mode 1 "// &
                               " Fine Mode 2 "// &
                               " Fine Mode 3 "// &
                               " Fine Mode 4 " 
                               
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Muri%fm_mode)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Muri%fm_mode
         
         
           case ("sediment_class_muri")
            Sds_info(Var_Idx)%Standard_Name = "sediment_class_muri" 
            
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Number_Of_Flags = 3
            Sds_Info(Var_Idx)%Flags_String = "Sediment Mode = -9 (fill in) "// &
                               " Sediment Mode = 0 (clear water)"// &
                               " Sediment Mode = 1 (turbid water)"
                               
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Muri%sediment_class)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Muri%sediment_class
         case ("aerosol_QA_muri")
            Sds_info(Var_Idx)%Standard_Name = "aerosol_QA_muri" 
            
            Sds_Info(Var_Idx)%Scaling_Type =  0_int1
            Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_HDF =  DFNT_INT8
            Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF = NF90_BYTE
            Sds_Info(Var_Idx)%Number_Of_Flags = 3
            Sds_Info(Var_Idx)%Flags_String = "Quality Assurance = -9 (fill in)"// &
                               " Quality Assurance = 0 (poor) "// &
                               " Quality Assurance = 1 (good) "
                               
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Muri%aerosol_QA)) Sds_Info(Var_Idx)%Sds_Data_2d_I1 => Muri%aerosol_QA
         
         case ("aerosol_fine_mode_ratio_muri")
            Sds_info(Var_Idx)%Standard_Name = "aerosol_fine_mode_ratio_muri" 
            Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
            Sds_Info(Var_Idx)%Units = "none"
            if (allocated(Muri%fmr)) Sds_Info(Var_Idx)%Sds_Data_2d_R4 => Muri%fmr
             
             
         case default
            call MESG ( "Unknown Level2 Variable ==> "//trim(Sds_Info(Var_Idx)%Sds_Name), level=verb_lev % ERROR )
               
      end select
   enddo

   !-------------------------------------------------------------------------
   ! set attributes for members of the Ch structure
   !-------------------------------------------------------------------------
   call SET_CH_ATTRIBUTES()

   !-------------------------------------------------------------------
   ! with the information set above, compute other needed attributes
   !-------------------------------------------------------------------
   do Var_Idx = 1, Num_Level2_Vars_List
      call STANDARD_SDS_INFO(Var_Idx)
   enddo

   !-------------------------------------------------------------------
   ! turn off sds that are irrelevant due to channel choices
   !-------------------------------------------------------------------
   do Var_Idx = 1, Num_Level2_Vars_List
     if (.not. associated(Sds_Info(Var_Idx)%Sds_Data_1d_I1) .and. &
         .not. associated(Sds_Info(Var_Idx)%Sds_Data_1d_I4) .and. &
         .not. associated(Sds_Info(Var_Idx)%Sds_Data_1d_R4) .and. &
         .not. associated(Sds_Info(Var_Idx)%Sds_Data_2d_I1) .and. &
         .not. associated(Sds_Info(Var_Idx)%Sds_Data_2d_I2) .and. &
         .not. associated(Sds_Info(Var_Idx)%Sds_Data_2d_I4) .and. &
         .not. associated(Sds_Info(Var_Idx)%Sds_Data_2d_R4) .and. &
         .not. associated(Sds_Info(Var_Idx)%Sds_Data_3d_I1)) then

         Sds_Info(Var_Idx)%On_Flag = .false.

     endif
   enddo

end subroutine SETUP_LEVEL2_SDS_INFO
!====================================================================
! SUBROUTINE Name: DEFINE_HDF_FILE_STRUCTURES
!
! Function:
!   Determines the structure of the level 2 files.
!
! Description:
!   This subroutine, given the inputs, determins and opens/creates the apppropriate
!       Level 2 (pixel level output) for various types of files. In addition
!       the HDF global attributes are put into the various HDF files 
!
!====================================================================
subroutine DEFINE_HDF_FILE_STRUCTURES(Num_Scans, &
   Dir_Level2, &
   File_1b, &
   Level2_File_Flag, &
   c1,c2,a1_20, &
   a2_20, &
   nu_20, &
   a1_31, &
   a2_31, &
   nu_31, &
   a1_32, &
   a2_32, &
   nu_32, &
   Solar_Ch20_nu, &
   Sun_Earth_Distance, &
   Therm_Cal_1b,  &
   Ref_Cal_1b, &
   nav_opt, &
   Use_Sst_anal, &
   Modis_clr_alb_flag, &
   Nwp_Opt, &
   Ch1_gain_low, &
   Ch1_gain_high, &
   Ch1_switch_count, &
   Ch1_dark_count,&
   Ch2_gain_low, &
   Ch2_gain_high, &
   Ch2_switch_count, &
   Ch2_dark_count, &
   Ch3a_Gain_low, &
   Ch3a_Gain_high, &
   Ch3a_Switch_Count, &
   Ch3a_Dark_Count,&
   Start_Year, &
   End_Year, &
   Start_Day, &
   End_Day, &
   Start_Time, &
   End_Time)


 character(len=*), intent(in):: Dir_Level2
 integer(kind=int4), intent(in):: Num_Scans
 character(len=*), intent(in):: File_1b
 integer, intent(in):: Level2_File_Flag
 
 integer(kind=int4), intent(in) :: Modis_Clr_Alb_Flag
 integer(kind=int4), intent(in) :: Nwp_Opt
 integer, intent(in):: therm_cal_1b
 integer, intent(in):: nav_opt
 integer, intent(in):: use_sst_anal
 integer, intent(in):: Ref_cal_1b
 real(kind=real4), intent(in):: c1,c2,a1_20,a2_20,nu_20,a1_31,a2_31,nu_31,a1_32, &
                                a2_32,nu_32,solar_Ch20_nu,sun_earth_distance, &
                                Ch1_gain_low,Ch1_gain_high,Ch1_switch_count, &
                                Ch1_dark_count, Ch2_gain_low,Ch2_gain_high, &
                                Ch2_switch_count,Ch2_dark_count, Ch3a_gain_low,&
                                Ch3a_gain_high,Ch3a_switch_count, &
                                Ch3a_dark_count
  integer(kind=int4), intent(in):: Start_Time
  integer(kind=int4), intent(in):: End_Time
  integer(kind=int2), intent(in):: Start_Year
  integer(kind=int2), intent(in):: End_Year
  integer(kind=int2), intent(in):: Start_Day
  integer(kind=int2), intent(in):: End_Day
  character(len=4):: l1b_ext
  character(len=1020):: File_1b_root
  character(len=1020):: File_Level2

  character(len=8):: Orbit_Number_String

  integer(kind=int4):: blank_int4
  character(len=1):: blank_char
  real:: blank_real4
  real(kind=real4):: Resolution_KM
  integer:: Istatus_Sum
  integer:: Istatus
  integer:: ipos
  integer:: ilen
  integer:: ilen_wmo_id
  integer:: Var_Idx
  integer:: StrLen

  character(len=3):: Wmo_String
  !integer:: n
 
  ! HDF function declarations
  integer:: sfstart
  !integer:: sfcreate
  !integer:: sfscatt
  !integer:: sfsnatt
  !integer:: erstat

  !--- begin executable code
  blank_int4 = 0
  blank_real4 = 0.0
  blank_char = " "

  !--- make a file name base for output
  File_1b_Root = trim(file_1b)

  !-----------------------------------------------------------------------------
  ! perform sensor specific modification of level-2 file name
  !-----------------------------------------------------------------------------

  sensor_search: do

    !--- ISCCP-NG
    if (index(File_1b_Root,'ISCCP-NG_L1g') > 0) then
      write(Wmo_String,fmt="(I0.3)") WMO_Id_ISCCPNG
      ipos = index(file_1b,"__")
      ilen = len(trim(file_1b))
      ilen_wmo_id = len('wmo_id')
      File_1b_Root = trim(file_1b(1:ipos)) // trim(file_1b(ipos+ilen_wmo_id+3:ilen-3)) // "_" // Wmo_String
      exit
    end if

    !--- FY4A - shorten name.
    if (File_1b_Root(1:4) == 'FY4A') then
      File_1b_Root = "FY4A_"//File_1b_Root(45:48)//File_1b_Root(49:52)//"_"//File_1b_Root(53:56)//"_"//File_1b_Root(75:79)//".hdf"
      exit
    end if

    !--- goes-r native format - shorten name
    if (File_1b_Root(1:6) == 'OR_ABI') then
      File_1b_Root = File_1b_Root(1:41)
      exit
    end if

    !--- special processing for modis - remove hdf suffix
    l1b_ext = File_1b_Root(len_trim(File_1b_Root)-3:len_trim(File_1b_Root))
    if (trim(l1b_ext) == ".hdf") then
      File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-4)
      exit
    end if

    !--- special processing for viirs - remove hdf suffix - this hard coded for
    if (trim(Sensor%Sensor_Name) == 'VIIRS') then
      StrLen = len_trim(File_1b_Root)-34 -7
      File_1b_Root(1:Strlen) = File_1b_Root(7:len_trim(File_1b_Root)-34)
      exit
    endif

    !--- special processing for viirs nasa - remove nc suffix - this hard coded for
    if (trim(Sensor%Sensor_Name) == 'VIIRS-NASA') then
      write (Orbit_Number_String, fmt='(I0.8)' )  Image%Orbit_Number
      if (File_1b_Root(1:5) == 'VNP03') then
        if (.not. Sensor%Fusion_Flag) then
          File_1b_Root = "snpp_viirs"//File_1b_Root(9:len_trim(File_1b_Root)-3)//"_B"//Orbit_Number_String
        else
          File_1b_Root = "snpp_fusion"//File_1b_Root(9:len_trim(File_1b_Root)-3)//"_B"//Orbit_Number_String
        end if
      else if (File_1b_Root(1:5) == 'VGEOM') then 
        File_1b_Root = "snpp_viirs"//File_1b_Root(11:len_trim(File_1b_Root)-3)//"_B"//Orbit_Number_String
      else if (File_1b_Root(1:5) == 'VJ103') then
        if (.not. Sensor%Fusion_Flag) then
          File_1b_Root = "j1_viirs"//File_1b_Root(9:len_trim(File_1b_Root)-3)//"_B"//Orbit_Number_String
        else
          File_1b_Root = "j1_fusion"//File_1b_Root(9:len_trim(File_1b_Root)-3)//"_B"//Orbit_Number_String
        end if
      end if
      exit
    end if

    if (trim(Sensor%Sensor_Name) == 'MERSI-2') then
      File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-4)
      exit
    end if

    !--- special processing for ahi - remove B01.nc suffix - this hard coded for
    if (trim(Sensor%Sensor_Name) == 'AHI' .OR. trim(Sensor%Sensor_Name) == 'AHI9') then
      StrLen = len_trim(File_1b_Root) - 12 - 4
      File_1b_Root(1:StrLen) = File_1b_Root(4:len_trim(File_1b_Root)-12)
      exit
    endif

    !--- do this for GOES names which are named goesxx_1_year_jday_hour.area
    if (trim(Sensor%Sensor_Name) == 'GOES-IL-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'GOES-MP-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'GOES-RU-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'COMS-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'MTSAT-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'FY2-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'SEVIRI') then

      !-- remove area suffix
      l1b_ext = File_1b_Root(len_trim(File_1b_Root)-3:len_trim(File_1b_Root))
      if (trim(l1b_ext) == "area") then
        File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-5)
      end if
      !-- remove channel number
      if (trim(l1b_ext) == "area") then
        ipos = index(File_1b_Root, "_1_")
        ilen = len(File_1b_Root)
        File_1b_Root = File_1b_Root(1:ipos-1) // "_"//File_1b_Root(ipos+3:ilen)
      end if
      exit
    end if

    !--- add 1km tag for 1km GOES files
    if (index(Sensor%Sensor_Name,'GOES') > 0 .and. Sensor%Spatial_Resolution_Meters == 1000) then
      File_1b_Root = trim(File_1b_Root)//".1km"
      exit
    endif

    !--- EPS-SG
    !SGA1-VII-1B-RAD_C_EUMT_20200223000200_G_D_20200103171000_20200103171500_T_N____.nc
    if (Index(Sensor%Sensor_Name,'METIMAGE') > 0) then
      File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root) - 11)
      exit
    endif
  
    exit

  enddo sensor_search

  !--- additional level2 name modifications

  !--- the Domain name for lat/lon subset is set, use it in level2 name
  if (Nav%Limit_Flag /= 0) then
    File_1b_Root = trim(File_1b_Root)//"_"//trim(Nav%Domain_Name)
  endif

  !--- CALIOP collocation files add sufix
  if (Caliop_Flag) File_1b_Root = trim(File_1b_Root)//".caliop"
  
  ! -- 
  if ( AVHRR_Fusion_Flag) File_1b_Root = trim(File_1b_Root)//".hirs_avhrr_fusion"

  !--- VGAC - remove .nc
  if (Index(Sensor%Sensor_Name,'VGAC') > 0) then
    File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-3)
  endif

  !--- add 'clavrx_' to the file name output
  File_1b_Root = 'clavrx_' // trim(File_1b_Root)

  !--- set Resolution_KM for global attribute
  Resolution_KM = Sensor%Spatial_Resolution_Meters / 1000.0

  !-------------------------------------------------------------
  ! define chunking here
  !-------------------------------------------------------------

  !--- 3d , first dim is sds dependent
  Sds_Dims_3d(1) = Max_Num_Cld_Test_Bytes
  Sds_Dims_3d(2) = Image%Number_Of_Elements
  Sds_Dims_3d(3) = Num_Scans
  Sds_Chunk_Size_3d(1) = Max_Num_Cld_Test_Bytes
  Sds_Chunk_Size_3d(2) = Image%Number_Of_Elements !200
  Sds_Chunk_Size_3d(3) = Num_Scans !100

  !-- dimension of 2d variables
  Sds_Dims_2d(1) = Image%Number_Of_Elements
  Sds_Dims_2d(2) = Image%Number_Of_Lines
  Sds_Chunk_Size_2d(1) = Image%Number_Of_Elements !200
  Sds_Chunk_Size_2d(2) = min(Image%Number_Of_Lines,Image%Number_Of_Lines_Per_Segment) !100

  !-- dimension of 1d variables
  Sds_Dims_1d(1) = Image%Number_Of_Lines

  !-------------------------------------------------------------
  ! define compression here
  !-------------------------------------------------------------
  Comp_Type = 0                  !no compression
  Comp_Prm(1) = 0
  Comp_Prm(2) = 0

  if (Compress_Flag == 1) then  !gzip compression
   Comp_Type = 4
   Comp_Prm(1) = 6
   Comp_Prm(2) = 0
  endif

  if (Compress_Flag == 2) then  !szip compression
   Comp_Type = 5
   Comp_Prm(1) = 32
   Comp_Prm(2) = 2
  endif


  !---------------------------------------------------------
  !-- define level2 file structure
  !---------------------------------------------------------
 if (Level2_File_Flag == sym%YES) then

     if (Output_Format_Flag == 0) File_Level2 = trim(File_1b_Root)//".level2.hdf"

     if (Output_Format_Flag == 1) File_Level2 = trim(File_1b_Root)//".level2.nc"

     call MESG(MOD_PROMPT//"creating level-2 file "//trim(File_Level2))

     if (Output_Format_Flag == 0) then 
        Sd_Id_Level2 = sfstart(trim(Dir_Level2)//trim(File_Level2),DFACC_CREATE)
     elseif (Output_Format_Flag == 1) then 
        Istatus_Sum = nf90_create(trim(Dir_Level2)//trim(File_Level2), &
                                  cmode=ior(NF90_CLOBBER,NF90_NETCDF4), &
                                  ncid = Sd_Id_Level2)
     endif

     if (Sd_Id_Level2 < 0) then
       call MESG(MOD_PROMPT//"Level-2 file creation failed. Exiting...",level = verb_lev % ERROR )
       stop 68
     endif

    !--- set up global attribute structure
    call SET_CLAVRX_GLOBAL_ATTRIBUTES("PIXEL",File_Level2,File_1b_Root, &
                           Resolution_KM, &
                           start_year,end_year,start_day,end_day,start_time,end_time,&
                           blank_int4,blank_int4,blank_char,blank_real4, &
                           therm_cal_1b,Ref_cal_1b,nav_opt,use_sst_anal, &
                           modis_clr_alb_flag, nwp_opt, Ch1_gain_low, Ch1_gain_high, &
                           Ch1_switch_count, Ch1_dark_count, &
                           Ch2_gain_low, Ch2_gain_high, &
                           Ch2_Switch_count, Ch2_Dark_count, &
                           Ch3a_gain_low, Ch3a_gain_high, &
                           Ch3a_switch_count, Ch3a_dark_count, &
                           sun_earth_distance, &
                           real(Sensor%Geo_Sub_Satellite_Longitude,kind=real4), &
                           real(Sensor%Geo_Sub_Satellite_Latitude,kind=real4), &
                           c1, c2, a1_20, a2_20, nu_20, &
                           a1_31, a2_31, nu_31, a1_32, a2_32, nu_32, &
                           solar_Ch20_nu,Nav%Timerr_seconds, Cloud_Mask_Mode, &
                           ACHA%Mode, ACHA%Mode_User, dcomp_mode,Sensor%WMO_Id, &
                           Sensor%Platform_Name, Sensor%Sensor_Name, &
                           Dark_Composite_Name,Bayesian_Cloud_Mask_Name, &
                           Nav%Limit_Flag,Nav%Lat_South_Limit, Nav%Lat_North_Limit, &
                           Nav%Lon_West_Limit, Nav%Lon_East_Limit,Nav%Domain_Name, &
                           Ch(31)%Max_Focal_Plane_Temp)

    !--- write attributes
    if (Output_Format_Flag == 0) call WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES(Sd_Id_Level2)
    if (Output_Format_Flag == 1) call WRITE_CLAVRX_NETCDF_GLOBAL_ATTRIBUTES(Sd_Id_Level2)

    !-- reset status flag for error checking
    Istatus_Sum = 0

    !------------------------------------------------------------------
    ! Loop through and define SDS using information in Sds_Info
    !------------------------------------------------------------------
    do Var_Idx = 1, Num_Level2_Vars_List
        if (Output_Format_Flag == 0) call DEFINE_LEVEL2_SDS_HDF(Var_Idx,Istatus)
        if (Output_Format_Flag == 1) call DEFINE_LEVEL2_SDS_NETCDF(Var_Idx,Istatus)
        if ( Istatus /= 0 ) then
          print *, EXE_PROMPT, MOD_PROMPT, "Error defining sds in level2 file for variable index "
          print *, Var_idx
          print *, EXE_PROMPT, MOD_PROMPT, Sds_Info(Var_Idx)%Sds_Name
          print *, 'Stop at line ', __LINE__ ,' in file ', __FILE__
          stop
        end if 
    enddo
    if (Output_Format_Flag == 1) then
        Istatus_Sum = nf90_enddef(Sd_Id_Level2) + Istatus_Sum
        if ( istatus_sum /= 0) then
          print *, EXE_PROMPT, MOD_PROMPT, "Error defining sds in level2 file for variable index "
          print *, EXE_PROMPT, MOD_PROMPT, Var_Idx
          print *, EXE_PROMPT, MOD_PROMPT, Sds_Info(Var_Idx)%Sds_Name 
          print *, 'Stop at line ', __LINE__ ,' in file ', __FILE__
          stop
        end if
    endif

    !--- check for and report errors
    if (Istatus_Sum /= 0) then
       print *, EXE_PROMPT, MOD_PROMPT, "Error defining sds in level2 file"
       stop
    endif

  endif

end subroutine DEFINE_HDF_FILE_STRUCTURES

!====================================================================
! SUBROUTINE Name: WRITE_SEGMENT_LEVEL2
!
! Function:
!   Writes output to appropriate Level 2 files
!
! Description:
!   This subroutine, given the flags that determine which files are 
!   being created, outputs the various pixel level outputs to the
!   appropriate level 2 files and appropriate SDSs for a given segment
!
!====================================================================
subroutine WRITE_SEGMENT_LEVEL2(Level2_File_Flag)

 integer, intent(in):: Level2_File_Flag
 integer:: Istatus
 integer:: Istatus_Sum
 !integer:: Line_Idx

 ! HDF function declarations
 !integer:: sfwdata
 integer:: Var_Idx

   !-----------------------------------------------------------------------
   ! Get time of each scan line and convert to scale
   !-----------------------------------------------------------------------
   !where(Image%Scan_Time_Ms == Missing_Value_Real4)
   where(Image%Scan_Time_Ms == Missing_Value_Int4)
      Image%Utc_Scan_Time_Hours = Missing_Value_Real4
   elsewhere
      Image%Utc_Scan_Time_Hours = Image%Scan_Time_Ms/60.0/60.0/1000.0
   endwhere

   !--------------------------------------------------------------------------------
   ! determine start and edges for writing this segment's data to pixel hdf file
   ! UNCLEAR IF HDF is zero based?  Modified for 1 based as NETCDF is one based.
   !-----------------------------------------------------------------------------
   if (Output_Format_Flag == 0) then  ! note that netcdf indices are 0-based
     Sds_Start_2d(1) = 0     !pixel dimension
     Sds_Start_2d(2) = Num_Scans_Level2_Hdf
     Sds_Stride_2d(1) = 1
     Sds_Stride_2d(2) = 1
     Sds_Edge_2d(1) = Image%Number_Of_Elements
     Sds_Edge_2d(2) = min(Image%Number_Of_Lines_Read_This_Segment,Image%Number_Of_Lines - Sds_Start_2d(2))
   endif

   if (Output_Format_Flag == 1) then  ! note that netcdf indices are 1-based
     Sds_Start_2d(1) = 1     !pixel dimension
     Sds_Start_2d(2) = Num_Scans_Level2_Hdf + 1
     Sds_Stride_2d(1) = 1
     Sds_Stride_2d(2) = 1
     Sds_Edge_2d(1) = Image%Number_Of_Elements
     Sds_Edge_2d(2) = min(Image%Number_Of_Lines_Read_This_Segment,Image%Number_Of_Lines)
   endif

   if (Sds_Edge_2d(2) <= 0) then
      return
   endif

    !--- update Num_Scans_Level2_Hdf
    Num_Scans_Level2_Hdf = min(Image%Number_Of_Lines,Num_Scans_Level2_Hdf +  &
                               Image%Number_Of_Lines_Read_This_Segment)
   !-------------------------------------------------------------------------
   ! write to level2 file
   !-------------------------------------------------------------------------
   if (Level2_File_Flag == sym%YES) then

      Istatus = 0
      Istatus_Sum = 0
      
      do Var_Idx = 1, Num_Level2_Vars_List

         if (.not. Sds_Info(Var_Idx)%On_Flag) cycle

         if (Output_Format_Flag == 0) then
            call WRITE_SDS_HDF(Var_Idx,Istatus)
         elseif (Output_Format_Flag == 1) then
            call WRITE_SDS_NETCDF(Var_Idx,Istatus)
         endif

      enddo

    !--- check for and report errors
    if (Istatus_Sum /= 0) then
       print *, EXE_PROMPT, MOD_PROMPT, "Error writing to level2 file: ", Istatus
       stop
    endif
    
   endif

end subroutine WRITE_SEGMENT_LEVEL2
!====================================================================
! SUBROUTINE Name: CLOSE_PIXEL_HDF_FILES
!
! Function:
!   Closes appropriate Level 2 files
!
! Description:
!   This subroutine, given the flags that determine which files have 
!   been created, closes the open output files.
!
!====================================================================
subroutine CLOSE_PIXEL_HDF_FILES(Level2_File_Flag)

 integer, intent(in):: Level2_File_Flag

 integer:: Var_Idx
 integer:: Istatus

! HDF function declarations
 integer:: sfsnatt
 integer:: sfendacc
 integer:: sfend
 integer:: sfscatt

!------------------------------------------------------------------------
!--- close level2 file
!------------------------------------------------------------------------
  Istatus = 0
  if (Output_Format_Flag == 0) then
     Istatus = sfsnatt(Sd_Id_Level2, "NUMBER_OF_ELEMENTS", DFNT_INT32,1,Image%Number_Of_Elements)+Istatus
     Istatus = sfsnatt(Sd_Id_Level2, "NUMBER_OF_SCANS_LEVEL1B", DFNT_INT32,1,Image%Number_Of_Lines)+Istatus
     Istatus = sfsnatt(Sd_Id_Level2, "NUMBER_OF_SCANS_LEVEL2", DFNT_INT32,1,Num_Scans_Level2_Hdf)+Istatus
     Istatus = sfsnatt(Sd_Id_Level2, "PROCESSING_TIME_MINUTES", DFNT_FLOAT32,1,Orbital_Processing_Time_Minutes)+Istatus
     Istatus = sfsnatt(Sd_Id_Level2, "NONCONFIDENT_CLOUD_MASK_FRACTION", DFNT_FLOAT32,1,NONCONFIDENT_CLOUD_MASK_Fraction)+Istatus
     Istatus = sfsnatt(Sd_Id_Level2, "ACHA_SUCCESS_FRACTION", DFNT_FLOAT32,1,ACHA%Success_Fraction)+Istatus
     Istatus = sfsnatt(Sd_Id_Level2, "DCOMP_SUCCESS_FRACTION", DFNT_FLOAT32,1,DCOMP_Success_Fraction)+Istatus
     Istatus = sfscatt(Sd_Id_Level2, "ACHA_MODE_FINAL", DFNT_CHAR8, len_trim(Clavrx_Global_Attr%acha_mode), trim(Clavrx_Global_Attr%acha_mode))+Istatus
  elseif (Output_Format_Flag == 1) then
     Istatus = nf90_put_att(Sd_Id_Level2, NF90_GLOBAL, "NUMBER_OF_ELEMENTS", Image%Number_Of_Elements)+Istatus
     Istatus = nf90_put_att(Sd_Id_Level2, NF90_GLOBAL, "NUMBER_OF_SCANS_LEVEL1B", Image%Number_Of_Lines)+Istatus
     Istatus = nf90_put_att(Sd_Id_Level2, NF90_GLOBAL, "NUMBER_OF_SCANS_LEVEL2", Num_Scans_Level2_Hdf)+Istatus
     Istatus = nf90_put_att(Sd_Id_Level2, NF90_GLOBAL, "PROCESSING_TIME_MINUTES", Orbital_Processing_Time_Minutes)+Istatus
     Istatus = nf90_put_att(Sd_Id_Level2, NF90_GLOBAL, "NONCONFIDENT_CLOUD_MASK_FRACTION", NONCONFIDENT_CLOUD_MASK_Fraction)+Istatus
     Istatus = nf90_put_att(Sd_Id_Level2, NF90_GLOBAL, "ACHA_SUCCESS_FRACTION", ACHA%Success_Fraction)+Istatus
     Istatus = nf90_put_att(Sd_Id_Level2, NF90_GLOBAL, "DCOMP_SUCCESS_FRACTION", DCOMP_Success_Fraction)+Istatus
     Istatus = nf90_put_att(Sd_Id_Level2, NF90_GLOBAL, "ACHA_MODE_FINAL", trim(Clavrx_Global_Attr%acha_mode))+Istatus
  endif

  if (Level2_File_Flag == sym%YES) then

   do Var_Idx = 1, Num_Level2_Vars_List


   ! --- Special case: add attributes to the ECM packed bits
   if (trim(Sds_Info(Var_Idx)%Sds_Name) == "cloud_mask_test_packed_results") then
     if (Output_Format_Flag == 0) then
       Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, 'number_classifiers', DFNT_INT32, 1, CLDMASK%N_Classifiers)+Istatus
       Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, 'classifier_names', DFNT_CHAR8, len_trim(CLDMASK%Classifiers_Names_Attr), &
                         (CLDMASK%Classifiers_Names_Attr))+Istatus
     elseif (Output_Format_Flag == 1) then
       Istatus = nf90_put_att(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, 'number_classifiers', CLDMASK%N_Classifiers)+Istatus
       Istatus = nf90_put_att(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, 'classifier_names', trim(CLDMASK%Classifiers_Names_Attr))+Istatus
     endif
   endif

      
      if (Output_Format_Flag == 0) then
         if (Sds_Info(Var_Idx)%Sds_Id /= 0 ) Istatus = sfendacc(Sds_Info(Var_Idx)%Sds_Id) + Istatus
      elseif (Output_Format_Flag == 1) then
         continue   !no netcdf equivalent ?
      endif
      Sds_Info(Var_Idx)%Sds_Id = 0

      !--- nullify sds_data points
      if (associated(Sds_Info(Var_Idx)%Sds_Data_1D_I1)) Sds_Info(Var_Idx)%Sds_Data_1D_I1 => null()
      if (associated(Sds_Info(Var_Idx)%Sds_Data_1D_I4)) Sds_Info(Var_Idx)%Sds_Data_1D_I4 => null()
      if (associated(Sds_Info(Var_Idx)%Sds_Data_2D_I1)) Sds_Info(Var_Idx)%Sds_Data_2D_I1 => null()
      if (associated(Sds_Info(Var_Idx)%Sds_Data_1D_R4)) Sds_Info(Var_Idx)%Sds_Data_1D_R4 => null()
      if (associated(Sds_Info(Var_Idx)%Sds_Data_2D_R4)) Sds_Info(Var_Idx)%Sds_Data_2D_R4 => null()
     
   enddo
   if (Output_Format_Flag == 0) then
      Istatus = sfend(Sd_Id_Level2) + Istatus
   elseif (Output_Format_Flag == 1) then
      !call CLOSE_NETCDF(Sd_Id_Level2)
      Istatus = nf90_close(Sd_Id_Level2)
   endif
  endif

end subroutine CLOSE_PIXEL_HDF_FILES

  !============================================================================
  !
  !============================================================================
  subroutine WRITE_ALGORITHM_ATTRIBUTES()

    integer:: sfscatt
    integer:: istatus
    
    istatus = 0
    istatus = sfscatt(Sd_Id_Level2, "CLOUD_MASK_VERSION", DFNT_CHAR8, len_trim(Cloud_Mask_Version), trim(Cloud_Mask_Version))+istatus
    istatus = sfscatt(Sd_Id_Level2, "CLOUD_MASK_THRESHOLDS_VERSION", DFNT_CHAR8,  &
                 len_trim(Cloud_Mask_Lut_Version), trim(Cloud_Mask_Lut_Version))+istatus 
    istatus = sfscatt(Sd_Id_Level2, "CLOUD_TYPE_VERSION", DFNT_CHAR8, len_trim(Cloud_Type_Version), trim(Cloud_Type_Version))+istatus 
    istatus = sfscatt(Sd_Id_Level2, "ACHA_VERSION", DFNT_CHAR8, len_trim(ACHA_Version), trim(ACHA_Version))+istatus 
    istatus = sfscatt(Sd_Id_Level2, "DCOMP_VERSION", DFNT_CHAR8, len_trim(dcomp_version), trim(dcomp_version))+istatus

  end subroutine WRITE_ALGORITHM_ATTRIBUTES

!============================================================================
! define variables and attributes for hdf level2 file
!============================================================================
subroutine DEFINE_LEVEL2_SDS_HDF(Var_Idx,Istatus_Sum)

   integer, intent(in):: Var_Idx

   integer, intent(out):: Istatus_Sum

   integer:: sfcreate, sfsnatt, sfscatt, sfdimid, sfsdmname, sfschnk
 
   integer:: Istatus

   integer:: Dim_Id

   Istatus = 0
   Istatus_Sum = 0

   select case(Sds_Info(Var_Idx)%Rank)

       case(1)
         Sds_Info(Var_Idx)%Sds_Id  = sfcreate(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Name,&
                     Sds_Info(Var_Idx)%Level2_Data_Type_HDF, &
                     Sds_Info(Var_Idx)%Rank,Sds_Dims_1d)
       case(2)
         Sds_Info(Var_Idx)%Sds_Id  = sfcreate(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Name,&
                     Sds_Info(Var_Idx)%Level2_Data_Type_HDF, &
                     Sds_Info(Var_Idx)%Rank,Sds_Dims_2d)
      case(3)
         Sds_Info(Var_Idx)%Sds_Id  = sfcreate(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Name,&
                     Sds_Info(Var_Idx)%Level2_Data_Type_HDF, &
                     Sds_Info(Var_Idx)%Rank,Sds_Dims_3d)
   end select

   !--- attributes
   Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "SCALED", DFNT_INT8, 1, Sds_Info(Var_Idx)%Scaling_Type) + Istatus

   Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "add_offset", DFNT_FLOAT32, 1, Sds_Info(Var_Idx)%Add_Offset) + Istatus

   Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "scale_factor", DFNT_FLOAT32, 1, Sds_Info(Var_Idx)%Scale_Factor) + Istatus

   Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "units", DFNT_CHAR8,  &
                     len_trim(Sds_Info(Var_Idx)%Units), trim(Sds_Info(Var_Idx)%Units)) + Istatus

   Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "standard_name", DFNT_CHAR8, len_trim(Sds_Info(Var_Idx)%Standard_Name),  &
                     trim(Sds_Info(Var_Idx)%Standard_Name)) + Istatus

   Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "long_name", DFNT_CHAR8, len_trim(Sds_Info(Var_Idx)%Long_Name),  &
                     trim(Sds_Info(Var_Idx)%Long_Name)) + Istatus

   if (trim(Sds_Info(Var_Idx)%Sds_Name) /= "latitude" .and. trim(Sds_Info(Var_Idx)%Sds_Name) /= "longitude") then
     Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "coordinates", DFNT_CHAR8, len_trim(coordinates_string),  &
                          trim(coordinates_string)) + Istatus

   endif


   !--- Fill Value
   select case(Sds_Info(Var_Idx)%Level2_Data_Type_HDF)

          case(DFNT_INT8)
            Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "_FillValue", DFNT_INT8, 1, int(Sds_Info(Var_Idx)%Fill_Value,kind=int1)) + Istatus

          case(DFNT_INT16)
            Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "_FillValue", DFNT_INT16, 1, int(Sds_Info(Var_Idx)%Fill_Value,kind=int2)) + Istatus

          case(DFNT_FLOAT32)
            Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "_FillValue", DFNT_FLOAT32, 1, real(Sds_Info(Var_Idx)%Fill_Value,kind=real4)) + Istatus

          case(DFNT_FLOAT64)
            Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "_FillValue", DFNT_FLOAT64, 1, real(Sds_Info(Var_Idx)%Fill_Value,kind=real8)) + Istatus
    end select

   !--- Valid Range
   select case(Sds_Info(Var_Idx)%Level2_Data_Type_HDF)

          case(DFNT_INT8)
            Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "valid_range", DFNT_INT8, 2, int(Sds_Info(Var_Idx)%Valid_Range,kind=int1)) + Istatus

          case(DFNT_INT16)
            Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "valid_range", DFNT_INT16, 2, int(Sds_Info(Var_Idx)%Valid_Range,kind=int2)) + Istatus

          case(DFNT_FLOAT32)
            Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "valid_range", DFNT_FLOAT32, 2, real(Sds_Info(Var_Idx)%Valid_Range,kind=real4)) + Istatus

          case(DFNT_FLOAT64)
            Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "valid_range", DFNT_FLOAT64, 2, real(Sds_Info(Var_Idx)%Valid_Range,kind=real8)) + Istatus
    end select

    !--- Actual Missing - NEEDED???
   if (Sds_Info(Var_Idx)%Scaling_Type > 0) then
    Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "actual_range", DFNT_FLOAT32, 2, Sds_Info(Var_Idx)%Actual_Range) + Istatus
    Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "actual_missing", DFNT_FLOAT32, 1, Sds_Info(Var_Idx)%Actual_Missing) + Istatus
   endif

   !--- dimension ids
   if (Sds_Info(Var_Idx)%Rank > 1) then
      Dim_Id = sfdimid(Sds_Info(Var_Idx)%Sds_Id, 0)
      Istatus = sfsdmname(Dim_Id,"pixel_elements_along_scan_direction") + Istatus
      Dim_Id = sfdimid(Sds_Info(Var_Idx)%Sds_Id, 1)
      Istatus = sfsdmname(Dim_Id,"scan_lines_along_track_direction") + Istatus
   endif

   if (Sds_Info(Var_Idx)%Rank == 1) then
      Dim_Id = sfdimid(Sds_Info(Var_Idx)%Sds_Id, 0)
      Istatus = sfsdmname(Dim_Id,"FakeDim1D") + Istatus
   endif

   !--- Compression (no compression of 1d output)
   if (Sds_Info(Var_Idx)%Rank == 2) then
       Istatus = sfschnk(Sds_Info(Var_Idx)%Sds_Id,Sds_Chunk_Size_2d,Comp_Type,Comp_Prm)+Istatus
   endif
   if (Sds_Info(Var_Idx)%Rank == 3) then
       Istatus = sfschnk(Sds_Info(Var_Idx)%Sds_Id,Sds_Chunk_Size_3d,Comp_Type,Comp_Prm)+Istatus
   endif

   !--- Flag Value REDO THIS!!!
   !call WRITE_FLAG_VALUES_HDF(Var_Idx,Istatus)

   if (Sds_Info(Var_Idx)%Number_Of_Flags > 0) then

      select case(Sds_Info(Var_Idx)%Number_Of_Flags)
         case(1)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 1, [0]) + Istatus
         case(2)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 2, [0,1]) + Istatus
         case(3)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 3, [0,1,2]) + Istatus
         case(4)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 4, [0,1,2,3]) + Istatus
         case(5)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 5, [0,1,2,3,4]) + Istatus
         case(6)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 6, [0,1,2,3,4,5]) + Istatus
         case(7)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 7, [0,1,2,3,4,5,6]) + Istatus
         case(8)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 8, [0,1,2,3,4,5,6,7]) + Istatus
         case(9)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 9, [0,1,2,3,4,5,6,7,8]) + Istatus
         case(10)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 10, [0,1,2,3,4,5,6,7,8,9]) + Istatus
         case(11)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 11, [0,1,2,3,4,5,6,7,8,9,10]) + Istatus
         case(12)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 12, [0,1,2,3,4,5,6,7,8,9,10,11]) + Istatus
         case(13)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 13, [0,1,2,3,4,5,6,7,8,9,10,11,12]) + Istatus
         case(14)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 14, [0,1,2,3,4,5,6,7,8,9,10,11,12,13]) + Istatus
         case(15)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 15, [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]) + Istatus
         case(16)
           Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 16, [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]) + Istatus
      end select


      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8,  &
                        len_trim(Sds_Info(Var_Idx)%Flags_String), Sds_Info(Var_Idx)%Flags_String) + Istatus

   endif

end subroutine DEFINE_LEVEL2_SDS_HDF

!============================================================================
! define variables and attributes for netcdf level2 file
!============================================================================
subroutine DEFINE_LEVEL2_SDS_NETCDF(Var_Idx,Istatus_Sum)

   integer, intent(in):: Var_Idx

   integer, intent(out):: Istatus_Sum

   integer:: Istatus, Deflate_Level_User

   integer, save:: Dim_Id_X, Dim_Id_Y, Dim_Id_Z

   integer:: Chunksize_X, Chunksize_Y, Chunksize_Z

   Istatus = 0
   Istatus_Sum = 0

   !--- define dimensions only needs to be done once
   if (Var_Idx == 1) then
      Istatus = nf90_def_dim(Sd_Id_Level2,"pixel_elements_along_scan_direction",Image%Number_Of_Elements, Dim_Id_X) + Istatus
      Istatus = nf90_def_dim(Sd_Id_Level2,"scan_lines_along_track_direction",Image%Number_Of_Lines, Dim_Id_Y) + Istatus
   endif

   Deflate_Level_User = 2   ! 1 to 9

   Chunksize_X = min(Image%Number_Of_Elements,200)
   Chunksize_Y = min(Image%Number_Of_Lines,100)
   Chunksize_Y = min(Image%Number_Of_Lines_Per_Segment,Chunksize_Y)
   Chunksize_Z = 1

   select case(Sds_Info(Var_Idx)%Rank)

       case(1)

         if (Compress_Flag == 0) then 
            Istatus = nf90_def_var(ncid = Sd_Id_Level2, &
                                name = trim(Sds_Info(Var_Idx)%Sds_Name), &
                                xtype = Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF, &
                                dimids = [Dim_Id_Y], &
                                varid = Sds_Info(Var_Idx)%Sds_Id) + Istatus
         else
            Istatus = nf90_def_var(ncid = Sd_Id_Level2, &
                                name = trim(Sds_Info(Var_Idx)%Sds_Name), &
                                xtype = Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF, &
                                dimids = [Dim_Id_Y], &
                                varid = Sds_Info(Var_Idx)%Sds_Id, &
                                chunksizes = [Chunksize_Y], &
                                deflate_level = Deflate_Level_User, &
                                shuffle = .true.) + Istatus
         endif

       case(2)

         if (Compress_Flag == 0) then 
            Istatus = nf90_def_var(Sd_Id_Level2,  &
                                name = Sds_Info(Var_Idx)%Sds_Name, &
                                xtype = Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF, &
                                dimids = [Dim_Id_X,Dim_Id_Y], &
                                varid = Sds_Info(Var_Idx)%Sds_Id) + Istatus
         else
            Istatus = nf90_def_var(Sd_Id_Level2,  &
                                name = Sds_Info(Var_Idx)%Sds_Name, &
                                xtype = Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF, &
                                dimids = [Dim_Id_X,Dim_Id_Y], &
                                varid = Sds_Info(Var_Idx)%Sds_Id, &
                                chunksizes = [Chunksize_X,Chunksize_Y], &
                                deflate_level = Deflate_Level_User, &
                                shuffle = .true.) + Istatus
         endif

      case(3)

         Istatus = nf90_def_dim(Sd_Id_Level2,"the_third_dimension",size(Sds_Info(Var_Idx)%Sds_Data_3d_I1,1), Dim_Id_Z) + Istatus
         if (Compress_Flag == 0) then 
            Istatus = nf90_def_var(Sd_Id_Level2,  &
                                name = Sds_Info(Var_Idx)%Sds_Name, &
                                xtype = Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF, &
                                dimids = [Dim_Id_Z,Dim_Id_X,Dim_Id_Y], &
                                varid = Sds_Info(Var_Idx)%Sds_Id) + Istatus
         else
            Istatus = nf90_def_var(Sd_Id_Level2,  &
                                name = Sds_Info(Var_Idx)%Sds_Name, &
                                xtype = Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF, &
                                dimids = [Dim_Id_Z,Dim_Id_X,Dim_Id_Y], &
                                varid = Sds_Info(Var_Idx)%Sds_Id, &
                                chunksizes = [Chunksize_Z,Chunksize_X,Chunksize_Y], &
                                deflate_level = Deflate_Level_User, &
                                shuffle = .true.) + Istatus
         endif

   end select

   !--- attributes
   Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "SCALED", Sds_Info(Var_Idx)%Scaling_Type) + Istatus
   if (Sds_Info(Var_Idx)%Scaling_Type == 1) then

   Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "add_offset", Sds_Info(Var_Idx)%Add_Offset) + Istatus

   Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "scale_factor", Sds_Info(Var_Idx)%Scale_Factor) + Istatus

   endif

   Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "units", trim(Sds_Info(Var_Idx)%Units)) + Istatus

   Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "standard_name", trim(Sds_Info(Var_Idx)%Standard_Name)) + Istatus

   Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "long_name", trim(Sds_Info(Var_Idx)%Long_Name)) + Istatus

   if (trim(Sds_Info(Var_Idx)%Sds_Name) /= "latitude" .and. trim(Sds_Info(Var_Idx)%Sds_Name) /= "longitude") then
     Istatus = nf90_put_att(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, "coordinates",trim(coordinates_string)) + Istatus
   endif


   !--- Fill Value
   select case(Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF)

          case(NF90_BYTE)
            Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "_FillValue", int(Sds_Info(Var_Idx)%Fill_Value,kind=int1)) + Istatus

          case(NF90_SHORT)
            Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "_FillValue", int(Sds_Info(Var_Idx)%Fill_Value,kind=int2)) + Istatus

          case(NF90_FLOAT)
            Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "_FillValue", real(Sds_Info(Var_Idx)%Fill_Value,kind=real4)) + Istatus

          case(NF90_DOUBLE)
            Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "_FillValue", real(Sds_Info(Var_Idx)%Fill_Value,kind=real8)) + Istatus
    end select

   !--- Valid Range
   if (Sds_Info(Var_Idx)%Scaling_Type == 1) then
   select case(Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF)

          case(NF90_BYTE)
            Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "valid_range",int(Sds_Info(Var_Idx)%Valid_Range,kind=int1)) + Istatus

          case(NF90_SHORT)
            Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "valid_range", int(Sds_Info(Var_Idx)%Valid_Range,kind=int2)) + Istatus

          case(NF90_FLOAT)
            Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "valid_range", int(Sds_Info(Var_Idx)%Valid_Range,kind=real4)) + Istatus

          case(NF90_DOUBLE)
            Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "valid_range", real(Sds_Info(Var_Idx)%Valid_Range,kind=real8)) + Istatus
    end select
    endif

    !--- Actual Missing - NEEDED???
   if (Sds_Info(Var_Idx)%Scaling_Type > 0) then
    Istatus = nf90_put_att(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, "actual_range", Sds_Info(Var_Idx)%Actual_Range) + Istatus
    Istatus = nf90_put_att(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, "actual_missing", Sds_Info(Var_Idx)%Actual_Missing) + Istatus
   endif

   !--- Compression (no compression of 1d output)
!  if (Sds_Info(Var_Idx)%Rank == 2) then
!      Istatus = sfschnk(Sds_Info(Var_Idx)%Sds_Id,Sds_Chunk_Size_2d,Comp_Type,Comp_Prm)+Istatus
!  endif
!  if (Sds_Info(Var_Idx)%Rank == 3) then
!      Istatus = sfschnk(Sds_Info(Var_Idx)%Sds_Id,Sds_Chunk_Size_3d,Comp_Type,Comp_Prm)+Istatus
!  endif

   !--- Flag Value REDO THIS!!!
   !call WRITE_FLAG_VALUES_HDF(Var_Idx,Istatus)

   if (Sds_Info(Var_Idx)%Number_Of_Flags > 0) then

      select case(Sds_Info(Var_Idx)%Number_Of_Flags)
         case(1)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0]) + Istatus
         case(2)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1]) + Istatus
         case(3)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2]) + Istatus
         case(4)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3]) + Istatus
         case(5)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4]) + Istatus
         case(6)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5]) + Istatus
         case(7)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5,6]) + Istatus
         case(8)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5,6,7]) + Istatus
         case(9)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5,6,7,8]) + Istatus
         case(10)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5,6,7,8,9]) + Istatus
         case(11)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5,6,7,8,9,10]) + Istatus
         case(12)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5,6,7,8,9,10,11]) + Istatus
         case(13)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5,6,7,8,9,10,11,12]) + Istatus
         case(14)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5,6,7,8,9,10,11,12,13]) + Istatus
         case(15)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]) + Istatus
         case(16)
           Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_values", [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]) + Istatus
      end select


      Istatus = nf90_put_att(Sd_Id_Level2,Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", Sds_Info(Var_Idx)%Flags_String) + Istatus

   endif
   
   Istatus_sum = Istatus

end subroutine DEFINE_LEVEL2_SDS_NETCDF

!---------------------------------------------------------------------------------
! compute some standard variable attributes
!---------------------------------------------------------------------------------
subroutine STANDARD_SDS_INFO(Var_Idx)

    integer, intent(in):: Var_Idx

    !--------------------------------------------------------------------
    ! Set Fill and Missing Values and Valid Range
    !--------------------------------------------------------------------
    if (Output_Format_Flag == 0)  then
       select case(Sds_Info(Var_Idx)%Level2_Data_Type_HDF)

          case(DFNT_INT8)
            Sds_Info(Var_Idx)%Fill_Value = real(ONE_BYTE_MIN)
            Sds_Info(Var_Idx)%Actual_Missing = MISSING_VALUE_REAL4

          case(DFNT_INT16)
            Sds_Info(Var_Idx)%Fill_Value = real(MISSING_VALUE_INT2)
            Sds_Info(Var_Idx)%Actual_Missing = MISSING_VALUE_REAL4

          case(DFNT_FLOAT32)
            Sds_Info(Var_Idx)%Fill_Value = MISSING_VALUE_REAL4
            Sds_Info(Var_Idx)%Actual_Missing = MISSING_VALUE_REAL4

          case(DFNT_FLOAT64)
            Sds_Info(Var_Idx)%Fill_Value = MISSING_VALUE_REAL4
            Sds_Info(Var_Idx)%Actual_Missing = MISSING_VALUE_REAL4

        end select

        if (Sds_Info(Var_Idx)%Scaling_Type /= 0) then
         select case(Sds_Info(Var_Idx)%Level2_Data_Type_HDF)

          case(DFNT_INT8)
            Sds_Info(Var_Idx)%Valid_Range = real([ONE_BYTE_MIN,ONE_BYTE_MAX])

          case(DFNT_INT16)
            Sds_Info(Var_Idx)%Valid_Range = real([TWO_BYTE_MIN,TWO_BYTE_MAX])

          case(DFNT_FLOAT32)
            Sds_Info(Var_Idx)%Valid_Range = real([-1.0*huge(real(1.0,kind=real4)),huge(real(1.0,kind=real4))])

          case(DFNT_FLOAT64)
            Sds_Info(Var_Idx)%Valid_Range = real([-1.0*huge(real(1.0,kind=real8)),huge(real(1.0,kind=real8))])
          end select

        endif
    endif

    if (Output_Format_Flag == 1)  then
       select case(Sds_Info(Var_Idx)%Level2_Data_Type_NETCDF)

          case(NF90_BYTE)
            Sds_Info(Var_Idx)%Fill_Value = real(ONE_BYTE_MIN)
            Sds_Info(Var_Idx)%Actual_Missing = MISSING_VALUE_REAL4

          case(NF90_SHORT)
            Sds_Info(Var_Idx)%Fill_Value = real(MISSING_VALUE_INT2)
            Sds_Info(Var_Idx)%Actual_Missing = MISSING_VALUE_REAL4

          case(NF90_FLOAT)
            Sds_Info(Var_Idx)%Fill_Value = MISSING_VALUE_REAL4
            Sds_Info(Var_Idx)%Actual_Missing = MISSING_VALUE_REAL4

          case(NF90_DOUBLE)
            Sds_Info(Var_Idx)%Fill_Value = MISSING_VALUE_REAL4
            Sds_Info(Var_Idx)%Actual_Missing = MISSING_VALUE_REAL4
        end select

        if (Sds_Info(Var_Idx)%Scaling_Type /= 0) then
         select case(Sds_Info(Var_Idx)%Level2_Data_Type_HDF)
          case(NF90_BYTE)
            Sds_Info(Var_Idx)%Valid_Range = real([ONE_BYTE_MIN,ONE_BYTE_MAX])

          case(NF90_SHORT)
            Sds_Info(Var_Idx)%Valid_Range = real([TWO_BYTE_MIN,TWO_BYTE_MAX])

          case(NF90_FLOAT)
            Sds_Info(Var_Idx)%Valid_Range = real([-1.0*huge(real(1.0,kind=real4)),huge(real(1.0,kind=real4))])

          case(NF90_DOUBLE)
            Sds_Info(Var_Idx)%Valid_Range = real([-1.0*huge(real(1.0,kind=real8)),huge(real(1.0,kind=real8))])
         end select
        endif

    endif
    !--------------------------------------------------------------------
    ! Set Add Offset and Scale Factor
    !--------------------------------------------------------------------
    if (Sds_Info(Var_Idx)%Scaling_Type == 1) then   !SCALED VARIABLES

       !if (Sds_Info(Var_Idx)%Actual_Range(1) /= Sds_Info(Var_Idx)%Actual_Range(2)) then
       if (Sds_Info(Var_Idx)%Actual_Range(1) .ner. Sds_Info(Var_Idx)%Actual_Range(2)) then
          Sds_Info(Var_Idx)%Scale_Factor =  &
                         (Sds_Info(Var_Idx)%Actual_Range(2) - Sds_Info(Var_Idx)%Actual_Range(1)) / &
                         (Sds_Info(Var_Idx)%Valid_Range(2) - Sds_Info(Var_Idx)%Valid_Range(1))
          Sds_Info(Var_Idx)%Add_Offset =   &
                         Sds_Info(Var_Idx)%Actual_Range(1) - &
                         Sds_Info(Var_Idx)%Scale_Factor*Sds_Info(Var_Idx)%Valid_Range(1)
       else
          Sds_Info(Var_Idx)%Scale_Factor = 1.0
          Sds_Info(Var_Idx)%Add_Offset = 0.0
       endif


    elseif (Sds_Info(Var_Idx)%Scaling_Type == 0) then   !UNSCALED VARIABLES

       Sds_Info(Var_Idx)%Scale_Factor = 1.0
       Sds_Info(Var_Idx)%Add_Offset = 0.0

    else

        print *, "Unknown Scaling Type"

    endif

end subroutine STANDARD_SDS_INFO

!----------------------------------------------------------------------------------------------------
! Routine for writing 2d scaled variables (int8 and int16)
!----------------------------------------------------------------------------------------------------
subroutine WRITE_2D_REAL4_SCALED_SDS_HDF(Z,Z_Idx,Istatus)

   real, dimension(:,:), intent(in):: Z
   integer, intent(in):: Z_Idx
   integer, intent(out):: Istatus

   integer:: sfwdata

   Istatus = 0

   if (Sds_Info(Z_Idx)%Level2_Data_Type_HDF == DFNT_INT16) then

      call SCALE_VECTOR_I2_RANK2(Z, &
                                sym%LINEAR_SCALING, &
                                Sds_Info(Z_Idx)%Actual_Range(1), &
                                Sds_Info(Z_Idx)%Actual_Range(2), &
                                Sds_Info(Z_Idx)%Actual_Missing, &
                                Two_Byte_Temp)

      Istatus = sfwdata(Sds_Info(Z_Idx)%Sds_Id,Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d, &
                Two_Byte_Temp(:,1:Sds_Edge_2d(2))) + Istatus

   elseif (Sds_Info(Z_Idx)%Level2_Data_Type_HDF == DFNT_INT8) then

      call SCALE_VECTOR_I1_RANK2(Z, &
                                sym%LINEAR_SCALING, &
                                Sds_Info(Z_Idx)%Actual_Range(1), &
                                Sds_Info(Z_Idx)%Actual_Range(2), &
                                Sds_Info(Z_Idx)%Actual_Missing, &
                                One_Byte_Temp)
      Istatus = sfwdata(Sds_Info(Z_Idx)%Sds_Id,Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d, &
                One_Byte_Temp(:,1:Sds_Edge_2d(2))) + Istatus

   else

     Istatus = 1

   endif
   
end subroutine WRITE_2D_REAL4_SCALED_SDS_HDF
!----------------------------------------------------------------------------------------------------
! Routine for writing 2d scaled variables (int8 and int16)
!----------------------------------------------------------------------------------------------------
subroutine WRITE_2D_REAL4_SCALED_SDS_NETCDF(Z,Z_Idx,Istatus)

   real, dimension(:,:), intent(in):: Z
   integer, intent(in):: Z_Idx
   integer, intent(out):: Istatus

   Istatus = 0

   if (Sds_Info(Z_Idx)%Level2_Data_Type_HDF == DFNT_INT16) then

      call SCALE_VECTOR_I2_RANK2(Z, &
                                sym%LINEAR_SCALING, &
                                Sds_Info(Z_Idx)%Actual_Range(1), &
                                Sds_Info(Z_Idx)%Actual_Range(2), &
                                Sds_Info(Z_Idx)%Actual_Missing, &
                                Two_Byte_Temp)
      Istatus = nf90_put_var(Sd_Id_Level2, Sds_Info(Z_Idx)%Sds_Id, &
                             Two_Byte_Temp(:,1:Sds_Edge_2d(2)), &
                             Sds_Start_2d, Sds_Edge_2d, Sds_Stride_2d) + Istatus

   elseif (Sds_Info(Z_Idx)%Level2_Data_Type_HDF == DFNT_INT8) then

      call SCALE_VECTOR_I1_RANK2(Z, &
                                sym%LINEAR_SCALING, &
                                Sds_Info(Z_Idx)%Actual_Range(1), &
                                Sds_Info(Z_Idx)%Actual_Range(2), &
                                Sds_Info(Z_Idx)%Actual_Missing, &
                                One_Byte_Temp)

      Istatus = nf90_put_var(Sd_Id_Level2, Sds_Info(Z_Idx)%Sds_Id, &
                             One_Byte_Temp(:,1:Sds_Edge_2d(2)), &
                             Sds_Start_2d, Sds_Edge_2d, Sds_Stride_2d) + Istatus

   else

     Istatus = 1

   endif
   
end subroutine WRITE_2D_REAL4_SCALED_SDS_NETCDF

!----------------------------------------------------------------------------------------------------
!Add CF flag_meanings and flag_values attributes to relevant sdss
!----------------------------------------------------------------------------------------------------
subroutine WRITE_FLAG_VALUES_HDF(Var_Idx,Istatus)

   integer, intent(in):: Var_Idx
   integer, intent(out):: Istatus
   !integer (kind = int1), parameter  :: temp_indgen_2_int1  (2)  = [0,1]
   !integer (kind = int1), parameter :: temp_indgen_4_int1(4) = [0,1,2,3]
   !integer (kind = int1), parameter :: temp_indgen_8_int1  (8) = [0,1,2,3,4,5,6,7]
   !integer (kind = int1), parameter  :: temp_indgen_13_int1 (13) = [0,1,2,3,4,5,6,7,8,9,10,11,12]
   !integer (kind = int1), parameter :: temp_indgen_14_int1 (14) = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
   integer (kind = int1), dimension(2), parameter  :: temp_indgen_2_int1  = int([0,1],kind=int1)
   integer (kind = int1), dimension(4), parameter :: temp_indgen_4_int1 = int([0,1,2,3],kind=int1)
   integer (kind = int1), dimension(8), parameter :: temp_indgen_8_int1 = int([0,1,2,3,4,5,6,7],kind=int1)
   integer (kind = int1), dimension(13), parameter  :: temp_indgen_13_int1 = int([0,1,2,3,4,5,6,7,8,9,10,11,12],kind=int1)
   integer (kind = int1), dimension(14), parameter :: temp_indgen_14_int1 = int([0,1,2,3,4,5,6,7,8,9,10,11,12,13],kind=int1)
   integer:: sfscatt, sfsnatt
   character(len=1020):: Flag_String



   !--- testing CF flag_meanings and flag_values attributes
   if (trim(Sds_Info(Var_Idx)%Sds_Name) == "acha_quality") then
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 2,  &
                        temp_indgen_2_int1) + Istatus
      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 143, &
                               "Processed "// &
                               " valid_Tc_retrieval "// &
                               " valid_ec_retrieval "// &
                               " valid_beta_retrieval "// &
                               " degraded_Tc_retrieval "// &
                               " degraded_ec_retrieval "// &
                               " degraded_beta_retrieval ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "acha_info") then
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 2, &
                        temp_indgen_2_int1) + Istatus
      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 196, &
                               "Cloud_Height_Attempted "// &
                               " Bias_Correction_Employed "// &
                               " Ice_Cloud_Retrieval "// &
                               " Local_Radiative_Center_Processing_Used "// &
                               " Multi_Layer_Retrieval "// &
                               " Lower_Cloud_Interpolation_Used "// &
                               " Boundary_Layer_Inversion_Assumed ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "dcomp_quality") then
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 2, &
                        temp_indgen_2_int1) + Istatus

      !--- which is correct?
      Flag_String = "Processed "// &
                    " valid_COD_retrieval "// &
                    " valid_REF_retrieval "// &
                    " degraded_COD_retrieval "// &
                    " degraded_REF_retrieval "// &
                    " convergency "// &
                    " glint "

      Flag_String = "quality flags for DCOMP products "// &
                    "1:Processed (0=no,1=yes) "// &
                    "2:valid COD retrieval (0=yes,1=no) "// &
                    "3:valid REF retrieval (0=yes,1=no) "// &
                    "4:degraded COD retrieval (0=no,1=degraded) "// &
                    "5:degraded REF retrieval (0=no,1=degraded) "// &
                    "6:convergency (0=no,1=yes) "// &
                    "7:glint (0=no,1=yes) "

      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8,  &
                        len_trim(Flag_String),trim(Flag_String)) + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "dcomp_info") then
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 2, &
                        temp_indgen_2_int1) + Istatus

      Flag_String = "processing flags for DCOMP "// &
                    "1:info flag set ? (0=no,1=yes) "// &
                    "2:land/sea mask (0=land,1=sea) "// &
                    "3:day/night mask (0=Day,1=Night) "// &
                    "4:twilight (65-82 solar zenith) (0=no,1=yes) "// &
                    "5:snow (0=no,1= snow) "// &
                    "6:sea-ice (0=no,1=sea-ice) "// &
                    "7:phase (0=liquid,1=ice) "// &
                    "8:thick_cloud (0=no,1=yes) " // &
                    "9:thin_cloud (0=no,1=yes) "

      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8,  &
                        len_trim(Flag_String),trim(Flag_String)) + Istatus

    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "cld_opd_dcomp_qf") then
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 4, &
                        temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 49, &
                               "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "cld_reff_dcomp_qf") then
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 4, &
                                 temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 49, &
                               "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "cld_temp_acha_qf") then
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 4, &
                        temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 49, &
                               "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "cld_emiss_acha_qf") then
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 4, &
                        temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 49, &
                               "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "cld_beta_acha_qf") then
      
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 4, &
                        temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 49, &
                               "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "cloud_mask") then
    
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 4,  &
                        temp_indgen_4_int1) + Istatus
            Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 47, &
                               "clear "// &
                               " probably_clear "// &
                               " probably_cloudy "// &
                               " cloudy ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "cloud_type") then
      
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 13, &
                        temp_indgen_13_int1) + Istatus
            Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 120, &
                               "clear "// &
                               " probably_clear "// &
                               " fog "// &
                               " water "// &
                               " supercooled_water "// &
                               " mixed "// &
                               " opaque_ice "// &
                               " cirrus "// &
                               " overlapping "// &
                               " overshooting "// &
                               " unknown "// &
                               " dust "// &
                               " smoke ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "surface_type") then
      
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 14, &
                        temp_indgen_14_int1) + Istatus
            Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 106, &
                               "water "// &
                               " evergreen_needle "// &
                               " evergreen_broad "// &
                               " deciduous_needle "// &
                               " deciduous_broad "// &
                               " mixed_forest "// &
                               " woodlands "// &
                               " wooded_grass "// &
                               " closed_shrubs "// &
                               " open_shrubs "// &
                               " grasses "// &
                               " croplands "// &
                               " bare "// &
                               " urban ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "land_class") then
      
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", DFNT_INT8, 8,  &
                         temp_indgen_8_int1) + Istatus

      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 109, &
                        "ocean "// &
                        " land "// &
                        " coastline "// &
                        " shallow_inland_water "// &
                        " ephemeral_water "// &
                        " deep_inland_water "// &
                        " moderate_ocean "// &
                        " deep_ocean ") + Istatus
    end if

    if (trim(Sds_Info(Var_Idx)%Sds_Name) == "snow_class") then
     
      Istatus = sfsnatt(Sds_Info(Var_Idx)%Sds_Id, "flag_values", &
                        DFNT_INT8, 3,  temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Info(Var_Idx)%Sds_Id, "flag_meanings", DFNT_CHAR8, 30, &
                       "no_snow_or_ice "// &
                      " sea_ice "// &
                      " snow ") + Istatus
    end if

 
end subroutine WRITE_FLAG_VALUES_HDF
!------------------------------------------------------------------------------
!  HDF4 Write Routine
!------------------------------------------------------------------------------
subroutine WRITE_SDS_HDF(Var_Idx,Istatus_Sum)

   integer,intent(in):: Var_Idx
   integer,intent(out):: Istatus_Sum
   integer:: Istatus
   integer:: Line_Idx
     
   ! HDF function declarations
   integer:: sfwdata


   Istatus = 0
   Istatus_Sum = 0

   !----------------------------------------------------------------------------------------
   ! special sds that involve packing
   !----------------------------------------------------------------------------------------
          if (trim(Sds_Info(Var_Idx)%Sds_Name) == "packed_pixel_meta_data") then
            One_Byte_Temp = 0_int1
            Temp_Mask = 0_int1
            do Line_Idx = 1, Image%Number_Of_Lines_Per_Segment
               Temp_Mask(:,Line_Idx) = Sensor%Chan_On_Flag_Per_Line(6,Line_Idx)
            enddo
            One_Byte_Temp = ishft(CLDMASK%Bayes_Mask_Sfc_Type,3) + ishft(Temp_Mask,2)+ &
                       ishft(Solar_Contamination_Mask,1) + bad_pixel_mask
            Istatus_Sum = sfwdata(Sds_Info(Var_Idx)%Sds_Id, Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                              One_Byte_Temp(:, 1:Sds_Edge_2d(2))) + Istatus_Sum
            return
         endif

         !--- packed land cover (land,snow,coast masks)
         if (trim(Sds_Info(Var_Idx)%Sds_Name) == "packed_land_cover") then
            One_Byte_Temp = 0_int1
            One_Byte_Temp = ishft(Sfc%Land,5) + ishft(Sfc%Snow,3) + Sfc%Coast_Mask
            Istatus_Sum = sfwdata(Sds_Info(Var_Idx)%Sds_Id, Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                              One_Byte_Temp(:, 1:Sds_Edge_2d(2))) + Istatus_Sum
            return
         endif

         !--- cld mask test vector (first byte - acm only)
         if (trim(Sds_Info(Var_Idx)%Sds_Name) == "cloud_mask_test_packed_results") then     
             if (Output_Format_Flag == 0) then
                Sds_Start_3d(1) = 0
                Sds_Start_3d(2) = 0
             endif 
             if (Output_Format_Flag == 1) then
                Sds_Start_3d(1) = 1
                Sds_Start_3d(2) = 1
             endif 

             Sds_Start_3d(3) = Sds_Start_2d(2)

             Sds_Stride_3d(1) = 1
             Sds_Stride_3d(2) = 1
             Sds_Stride_3d(3) = 1

             Sds_Edge_3d(1) = Max_Num_Cld_Test_Bytes
             Sds_Edge_3d(2) = Sds_Edge_2d(1)
             Sds_Edge_3d(3) = Sds_Edge_2d(2)

             Istatus_Sum = sfwdata(Sds_Info(Var_Idx)%Sds_Id,  &
                               Sds_Start_3d, &
                               Sds_Stride_3d, &
                               Sds_Edge_3d, &
                               CLDMASK%Cld_Test_Vector_Packed(:,:,1:Sds_Edge_2d(2))) + Istatus_Sum
            return
         endif

         !----------------------------------------------------------------------------------------
         ! Standard Sds's
         !----------------------------------------------------------------------------------------

         !--- 1d Real4 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 1 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_FLOAT32) then
            Istatus_Sum = sfwdata(Sds_Info(Var_Idx)%Sds_Id, Sds_Start_2d(2), Sds_Stride_2d(2), Sds_Edge_2d(2),  &
                              Sds_Info(Var_Idx)%Sds_Data_1d_R4(1:Sds_Edge_2d(2))) + Istatus_Sum
         endif

         !--- 1d Int8 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 1 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_INT8) then
            Istatus_Sum = sfwdata(Sds_Info(Var_Idx)%Sds_Id, Sds_Start_2d(2), Sds_Stride_2d(2), Sds_Edge_2d(2),  &
                              Sds_Info(Var_Idx)%Sds_Data_1d_I1(1:Sds_Edge_2d(2))) + Istatus_Sum
         endif

         !--- 1d Int4 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 1 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_INT32) then
            Istatus_Sum = sfwdata(Sds_Info(Var_Idx)%Sds_Id, Sds_Start_2d(2), Sds_Stride_2d(2), Sds_Edge_2d(2),  &
                              Sds_Info(Var_Idx)%Sds_Data_1d_I4(1:Sds_Edge_2d(2))) + Istatus_Sum
         endif

         !--- 2d Int1 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 2 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_INT8) then
            Istatus_Sum = sfwdata(Sds_Info(Var_Idx)%Sds_Id, Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                              Sds_Info(Var_Idx)%Sds_Data_2d_I1(:,1:Sds_Edge_2d(2))) + Istatus_Sum
         endif

         !--- 2d Int2 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 2 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_INT16) then
            Istatus_Sum = sfwdata(Sds_Info(Var_Idx)%Sds_Id, Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                              Sds_Info(Var_Idx)%Sds_Data_2d_I2(:,1:Sds_Edge_2d(2))) + Istatus_Sum
         endif

         !--- 2d Int4 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 2 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_INT32) then
            Istatus_Sum = sfwdata(Sds_Info(Var_Idx)%Sds_Id, Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                              Sds_Info(Var_Idx)%Sds_Data_2d_I4(:,1:Sds_Edge_2d(2))) + Istatus_Sum
         endif

         !--- 2d Real4 Unscaled
         if (Sds_Info(Var_Idx)%Rank == 2 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_FLOAT32) then
             Istatus_Sum = sfwdata(Sds_Info(Var_Idx)%Sds_Id, Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                              Sds_Info(Var_Idx)%Sds_Data_2d_R4(:,1:Sds_Edge_2d(2))) + Istatus_Sum
             Istatus_Sum = Istatus_Sum + Istatus
         endif

         !--- 2d Real4 Scaled
         if (Sds_Info(Var_Idx)%Rank == 2 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 1 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_FLOAT32) then

            call WRITE_2D_REAL4_SCALED_SDS_HDF(Sds_Info(Var_Idx)%Sds_Data_2d_R4,Sds_Info(Var_Idx)%Sds_Idx,Istatus)
            Istatus_Sum = Istatus_Sum + Istatus

         endif

end subroutine WRITE_SDS_HDF

!------------------------------------------------------------------------------
!  NETCDF Write Routine
!------------------------------------------------------------------------------
subroutine WRITE_SDS_NETCDF(Var_Idx,Istatus_Sum)

   integer,intent(in):: Var_Idx
   integer,intent(out):: Istatus_Sum
   integer:: Istatus
   integer:: Line_Idx
     
   Istatus = 0
   Istatus_Sum = 0


   !----------------------------------------------------------------------------------------
   ! special sds that involve packing
   !----------------------------------------------------------------------------------------
          if (trim(Sds_Info(Var_Idx)%Sds_Name) == "packed_pixel_meta_data") then
            One_Byte_Temp = 0_int1
            Temp_Mask = 0_int1
            do Line_Idx = 1, Image%Number_Of_Lines_Per_Segment
               Temp_Mask(:,Line_Idx) = Sensor%Chan_On_Flag_Per_Line(6,Line_Idx)
            enddo
            One_Byte_Temp = ishft(CLDMASK%Bayes_Mask_Sfc_Type,3) + ishft(Temp_Mask,2)+ &
                       ishft(Solar_Contamination_Mask,1) + bad_pixel_mask

            Istatus_Sum = nf90_put_var(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, &
                                       One_Byte_Temp(:, 1:Sds_Edge_2d(2)), &
                                       Sds_Start_2d, Sds_Edge_2d, Sds_Stride_2d) + Istatus_Sum

            return
         endif

         !--- packed land cover (land,snow,coast masks)
         if (trim(Sds_Info(Var_Idx)%Sds_Name) == "packed_land_cover") then
            One_Byte_Temp = 0_int1
            One_Byte_Temp = ishft(Sfc%Land,5) + ishft(Sfc%Snow,3) + Sfc%Coast_Mask
            Istatus_Sum = nf90_put_var(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, &
                                       One_Byte_Temp(:, 1:Sds_Edge_2d(2)), &
                                       Sds_Start_2d, Sds_Edge_2d, Sds_Stride_2d) + Istatus_Sum

            return
         endif

         !--- cld mask test vector (first byte - acm only)
         if (trim(Sds_Info(Var_Idx)%Sds_Name) == "cloud_mask_test_packed_results") then     
             if (Output_Format_Flag == 0) then
                Sds_Start_3d(1) = 0
                Sds_Start_3d(2) = 0
             endif 
             if (Output_Format_Flag == 1) then
                Sds_Start_3d(1) = 1
                Sds_Start_3d(2) = 1
             endif 
             Sds_Start_3d(3) = Sds_Start_2d(2)

             Sds_Stride_3d(1) = 1
             Sds_Stride_3d(2) = 1
             Sds_Stride_3d(3) = 1

             Sds_Edge_3d(1) = Max_Num_Cld_Test_Bytes
             Sds_Edge_3d(2) = Sds_Edge_2d(1)
             Sds_Edge_3d(3) = Sds_Edge_2d(2)

            Istatus_Sum = nf90_put_var(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, &
                                       CLDMASK%Cld_Test_Vector_Packed(:,:,1:Sds_Edge_2d(2)), &
                                       Sds_Start_3d, Sds_Edge_3d, Sds_Stride_3d) + Istatus_Sum
            return
         endif

         !----------------------------------------------------------------------------------------
         ! Standard Sds's
         !----------------------------------------------------------------------------------------

         !--- 1d Real4 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 1 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_FLOAT32) then
            Istatus_Sum = nf90_put_var(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, &
                                       Sds_Info(Var_Idx)%Sds_Data_1d_R4(1:Sds_Edge_2d(2)), &
                                       start = [Sds_Start_2d(2)], &
                                       count = [Sds_Edge_2d(2)],  &
                                       stride = [Sds_Stride_2d(2)]) + Istatus_Sum
         endif

         !--- 1d Int8 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 1 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_INT8) then
             Istatus_Sum = nf90_put_var(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, &
                                       Sds_Info(Var_Idx)%Sds_Data_1d_I1(1:Sds_Edge_2d(2)), &
                                       start = [Sds_Start_2d(2)], &
                                       count = [Sds_Edge_2d(2)],  &
                                       stride = [Sds_Stride_2d(2)]) + Istatus_Sum
         endif

         !--- 1d Int4 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 1 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_INT32) then
             Istatus_Sum = nf90_put_var(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, &
                                       Sds_Info(Var_Idx)%Sds_Data_1d_I4(1:Sds_Edge_2d(2)), &
                                       start = [Sds_Start_2d(2)], &
                                       count = [Sds_Edge_2d(2)], &
                                       stride = [Sds_Stride_2d(2)]) + Istatus_Sum
         endif

         !--- 2d Int1 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 2 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_INT8) then
             Istatus_Sum = nf90_put_var(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, &
                                       Sds_Info(Var_Idx)%Sds_Data_2d_I1(:,1:Sds_Edge_2d(2)), &
                                       start = Sds_Start_2d, &
                                       count = Sds_Edge_2d,  &
                                       stride = Sds_Stride_2d) + Istatus_Sum
         endif

         !--- 2d Int2 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 2 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_INT16) then
             Istatus_Sum = nf90_put_var(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, &
                                       Sds_Info(Var_Idx)%Sds_Data_2d_I2(:,1:Sds_Edge_2d(2)), &
                                       start = Sds_Start_2d, &
                                       count = Sds_Edge_2d,  &
                                       stride = Sds_Stride_2d) + Istatus_Sum
         endif

         !--- 2d Int4 UnScaled
         if (Sds_Info(Var_Idx)%Rank == 2 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_INT32) then
             Istatus_Sum = nf90_put_var(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, &
                                       Sds_Info(Var_Idx)%Sds_Data_2d_I4(:,1:Sds_Edge_2d(2)), &
                                       start = Sds_Start_2d, &
                                       count = Sds_Edge_2d,  &
                                       stride = Sds_Stride_2d) + Istatus_Sum
         endif

         !--- 2d Real4 Unscaled
         if (Sds_Info(Var_Idx)%Rank == 2 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 0 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_FLOAT32) then

             Istatus_Sum = nf90_put_var(Sd_Id_Level2, Sds_Info(Var_Idx)%Sds_Id, &
                                       Sds_Info(Var_Idx)%Sds_Data_2d_R4(:,1:Sds_Edge_2d(2)), &
                                       start = Sds_Start_2d, &
                                       count = Sds_Edge_2d,  &
                                       stride = Sds_Stride_2d) + Istatus_Sum
         endif

         !--- 2d Real4 Scaled
         if (Sds_Info(Var_Idx)%Rank == 2 .and. &
             Sds_Info(Var_Idx)%Scaling_Type == 1 .and. &
             Sds_Info(Var_Idx)%Input_Data_Type_HDF == DFNT_FLOAT32) then

            call WRITE_2D_REAL4_SCALED_SDS_NETCDF(Sds_Info(Var_Idx)%Sds_Data_2d_R4,Sds_Info(Var_Idx)%Sds_Idx,Istatus)

            Istatus_Sum = Istatus_Sum + Istatus

         endif

end subroutine WRITE_SDS_NETCDF
!-------------------------------------------------------------------
! set units and actual_ranges for reflectances and brightness temps
!-------------------------------------------------------------------
 subroutine SET_CH_ATTRIBUTES()

   integer:: Var_Idx
   character(len=20):: Temp_String_1
   !character(len=4):: Temp_String_2
   character(len=20):: Temp_String_2
   character(len=6):: Wvl_String
   integer:: ipos1, ipos2

   do Var_Idx = 1, Num_Level2_Vars_List

      !--- reflectances
      if (index(Sds_Info(Var_Idx)%Sds_Name, 'refl_') == 1 .and. &
          index(Sds_Info(Var_Idx)%Sds_Name, '_nom') >0) then 
         Sds_Info(Var_Idx)%Actual_Range = [-2.0,120.0]
         Sds_Info(Var_Idx)%Units =  "%"
         Sds_Info(Var_Idx)%Standard_Name =  "toa_bidirectional_reflectance"
         Sds_Info(Var_Idx)%Long_Name =  "top of atmosphere reflectance at the nominal wavelength of"

         if (index(Sds_Info(Var_Idx)%Sds_Name,'_std_') == 1) then
             Sds_Info(Var_Idx)%Actual_Range = [0.0,100.0]
             Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_reflectance_sub-pixel_stddev"
             Sds_Info(Var_Idx)%Long_Name =  "top of atmosphere sub-pixel reflectance sttdev at the nominal wavelength of"
         endif

         !--- add wavelength
         temp_string_1 ="refl_"
         temp_string_2 ="_nom"
         ipos1 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_1))  + len_trim(temp_string_1)
         ipos2 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_2)) - 1
         Wvl_String = trim(Sds_Info(Var_Idx)%Sds_Name(ipos1:ipos2))
         Sds_Info(Var_Idx)%Long_Name = trim(Sds_Info(Var_Idx)%Long_Name) // " " //trim(Wvl_String)
      endif

      !--- surface reflectances atmospherically corrected
      if (index(Sds_Info(Var_Idx)%Sds_Name, 'refl_') == 1 .and. &
          index(Sds_Info(Var_Idx)%Sds_Name, '_nom_atmos_corr') >0) then
         Sds_Info(Var_Idx)%Actual_Range = [-2.0,120.0]
         Sds_Info(Var_Idx)%Units =  "%"
         Sds_Info(Var_Idx)%Standard_Name = "toa_bidirectional_pseudo_reflectance_"
         Sds_Info(Var_Idx)%Long_Name =  "observed pseudo-reflectance, atmospherically corrected, " // &
                                         "at the nominal wavelength of"

         !--- add wavelength
         temp_string_1 ="refl_"
         temp_string_2 ="_nom_atmos_corr"
         ipos1 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_1))  + len_trim(temp_string_1)
         ipos2 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_2)) - 1
         Wvl_String = trim(Sds_Info(Var_Idx)%Sds_Name(ipos1:ipos2))
         Sds_Info(Var_Idx)%Long_Name = trim(Sds_Info(Var_Idx)%Long_Name) // " " //trim(Wvl_String)
         Sds_Info(Var_Idx)%Standard_Name = trim(Sds_Info(Var_Idx)%Standard_Name) // " " // &
                                           trim(Wvl_String) // "_atmos_corr"
      endif

      !--- 3.75 micron pseudo reflectance
      if (trim(Sds_Info(Var_Idx)%Sds_Name) == "refl_3_75um_nom") then
         Sds_Info(Var_Idx)%Actual_Range = [-20.0,80.0]
         Sds_Info(Var_Idx)%Units =  "%"
      endif

      !--- 3x3 stddev of reflectance
      if (index(Sds_Info(Var_Idx)%Sds_Name, 'refl_') == 1 .and. &
          index(Sds_Info(Var_Idx)%Sds_Name, '_nom') > 0 .and. & 
          index(Sds_Info(Var_Idx)%Sds_Name, '_stddev') > 0) then
          Sds_Info(Var_Idx)%Actual_Range = [0.0,20.0]
          Sds_Info(Var_Idx)%Units =  "%"
          Sds_Info(Var_Idx)%Standard_Name =  trim(Sds_Info(Var_Idx)%Standard_Name) // &
                                         " standard deviation over a 3x3 pixel array"
          Sds_Info(Var_Idx)%Long_Name =  trim(Sds_Info(Var_Idx)%Long_Name) // &
                                         " standard deviation over a 3x3 pixel array at the nominal wavelength of "
          !--- add wavelength
          temp_string_1 ="refl_"
          temp_string_2 ="_nom"
          ipos1 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_1))  + len_trim(temp_string_1)
          ipos2 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_2)) - 1
          Wvl_String = trim(Sds_Info(Var_Idx)%Sds_Name(ipos1:ipos2))
          Sds_Info(Var_Idx)%Long_Name = trim(Sds_Info(Var_Idx)%Long_Name) // " " //trim(Wvl_String)
      endif

      !--- brightness temperatures
      if (index(Sds_Info(Var_Idx)%Sds_Name, 'temp_') == 1 .and. &
          index(Sds_Info(Var_Idx)%Sds_Name, '_nom') >0) then 
          Sds_Info(Var_Idx)%Actual_Range = [160.0,340.0]
          Sds_Info(Var_Idx)%Units =  "K"
          Sds_Info(Var_Idx)%Standard_Name =  "toa_brightness_temperature"
          Sds_Info(Var_Idx)%Long_Name =  "top of atmosphere brightness temperature at the nominal wavelength of "

          !--- add wavelength
          temp_string_1 ="temp_"
          temp_string_2 ="_nom"
          ipos1 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_1))  + len_trim(temp_string_1)
          ipos2 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_2)) - 1
          Wvl_String = trim(Sds_Info(Var_Idx)%Sds_Name(ipos1:ipos2))
          Sds_Info(Var_Idx)%Long_Name = trim(Sds_Info(Var_Idx)%Long_Name) // " " //trim(Wvl_String)
      endif

      !--- 3x3 stddev of brightness temperature
      if (index(Sds_Info(Var_Idx)%Sds_Name, 'temp_') == 1 .and. &
          index(Sds_Info(Var_Idx)%Sds_Name, '_nom') > 0 .and. & 
          index(Sds_Info(Var_Idx)%Sds_Name, '_stddev') > 0) then
          Sds_Info(Var_Idx)%Actual_Range = [0.0,20.0]
          Sds_Info(Var_Idx)%Units =  "K"
          Sds_Info(Var_Idx)%Standard_Name =  trim(Sds_Info(Var_Idx)%Standard_Name) // &
                                         " standard deviation over a 3x3 pixel array"
          Sds_Info(Var_Idx)%Long_Name =  trim(Sds_Info(Var_Idx)%Long_Name) // &
                                         " standard deviation over a 3x3 pixel array "// &
                                         "at the nominal wavelength of "
          !--- add wavelength
          temp_string_1 ="temp_"
          temp_string_2 ="_nom"
          ipos1 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_1))  + len_trim(temp_string_1)
          ipos2 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_2)) - 1
          Wvl_String = trim(Sds_Info(Var_Idx)%Sds_Name(ipos1:ipos2))
          Sds_Info(Var_Idx)%Long_Name = trim(Sds_Info(Var_Idx)%Long_Name) // " " //trim(Wvl_String)
      endif

      !--- sfc emiss
      if (index(Sds_Info(Var_Idx)%Sds_Name, 'emiss_sfc_') == 1 .and. &
          index(Sds_Info(Var_Idx)%Sds_Name, '_nom') >0) then 
         Sds_Info(Var_Idx)%Actual_Range = [0.75,1.0]
         Sds_Info(Var_Idx)%Level2_Data_Type_HDF = DFNT_INT16
         Sds_Info(Var_Idx)%Standard_Name =  "surface_emissivity"
         Sds_Info(Var_Idx)%Long_Name =  "surface_emissivity at the nominal wavelength of"

         !--- add wavelength
         temp_string_1 ="emiss_sfc_"
         temp_string_2 ="_nom"
         ipos1 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_1))  + len_trim(temp_string_1)
         ipos2 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_2)) - 1
         Wvl_String = trim(Sds_Info(Var_Idx)%Sds_Name(ipos1:ipos2))
         Sds_Info(Var_Idx)%Long_Name = trim(Sds_Info(Var_Idx)%Long_Name) // " " //trim(Wvl_String)
      endif

      !--- tropopause emissivity 
      if (index(Sds_Info(Var_Idx)%Sds_Name, 'emiss_tropo_') == 1 .and. &
          index(Sds_Info(Var_Idx)%Sds_Name, '_nom') >0) then 
         Sds_Info(Var_Idx)%Actual_Range = [-0.5,1.2]
         Sds_Info(Var_Idx)%Level2_Data_Type_HDF = DFNT_INT16
         Sds_Info(Var_Idx)%Standard_Name =  "emissivity referenced to tropopause"
         Sds_Info(Var_Idx)%Long_Name =  "emissivity referenced to tropopause at nominal wavelength of"

         !--- add wavelength
         temp_string_1 ="emiss_tropo_"
         temp_string_2 ="_nom"
         ipos1 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_1))  + len_trim(temp_string_1)
         ipos2 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_2)) - 1
         Wvl_String = trim(Sds_Info(Var_Idx)%Sds_Name(ipos1:ipos2))
         Sds_Info(Var_Idx)%Long_Name = trim(Sds_Info(Var_Idx)%Long_Name) // " " //trim(Wvl_String)
      endif

      !---  white sky reflectance
      if (index(Sds_Info(Var_Idx)%Sds_Name, 'refl_sfc_white_sky_') == 1 .and. &
          index(Sds_Info(Var_Idx)%Sds_Name, '_nom') >0) then 
         Sds_Info(Var_Idx)%Actual_Range = [-20.0,100.0]
         Sds_Info(Var_Idx)%Level2_Data_Type_HDF = DFNT_INT8
         Sds_Info(Var_Idx)%Units =  "%"
         Sds_Info(Var_Idx)%Standard_Name =  "white_sky_surface_reflectance"
         Sds_Info(Var_Idx)%Long_Name =  "surface reflectance for white skies at nominal wavelength of"

         !--- add wavelength
         temp_string_1 ="refl_sfc_white_sky_"
         temp_string_2 ="_nom"
         ipos1 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_1))  + len_trim(temp_string_1)
         ipos2 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_2)) - 1
         Wvl_String = trim(Sds_Info(Var_Idx)%Sds_Name(ipos1:ipos2))
         Sds_Info(Var_Idx)%Long_Name = trim(Sds_Info(Var_Idx)%Long_Name) // " " //trim(Wvl_String)
      endif

      !---  atmospheric transmission
      if (index(Sds_Info(Var_Idx)%Sds_Name, 'trans_atm_') == 1 .and. &
          index(Sds_Info(Var_Idx)%Sds_Name, '_nom') >0) then 
         Sds_Info(Var_Idx)%Actual_Range = [0.0,1.0]
         Sds_Info(Var_Idx)%Level2_Data_Type_HDF = DFNT_INT8
         Sds_Info(Var_Idx)%Units =  "none"
         Sds_Info(Var_Idx)%Standard_Name = "path_transmission_from_clear_atmosphere"
         Sds_Info(Var_Idx)%Long_Name =  "path transmission at nominal wavelength of"

         !--- add wavelength
         temp_string_1 ="trans_atm_"
         temp_string_2 ="_nom"
         ipos1 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_1))  + len_trim(temp_string_1)
         ipos2 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_2)) - 1
         Wvl_String = trim(Sds_Info(Var_Idx)%Sds_Name(ipos1:ipos2))
         Sds_Info(Var_Idx)%Long_Name = trim(Sds_Info(Var_Idx)%Long_Name) // " " //trim(Wvl_String)
      endif

      !---  data quality flag
      if (index(Sds_Info(Var_Idx)%Sds_Name, 'dqf_') == 1 .and. &
          index(Sds_Info(Var_Idx)%Sds_Name, '_nom') >0) then 
         Sds_Info(Var_Idx)%Scaling_Type =  0_int1
         Sds_Info(Var_Idx)%Actual_Range = [0.0,4.0]
         Sds_Info(Var_Idx)%Input_Data_Type_HDF =  DFNT_INT8
         Sds_Info(Var_Idx)%Level2_Data_Type_HDF = DFNT_INT8
         Sds_Info(Var_Idx)%Units =  "none"
         Sds_Info(Var_Idx)%Standard_Name = "dqf"
         Sds_Info(Var_Idx)%Long_Name =  "data_quality_flag"

         !--- add wavelength
         temp_string_1 ="dqf_"
         temp_string_2 ="_nom"
         ipos1 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_1))  + len_trim(temp_string_1)
         ipos2 = index(trim(Sds_Info(Var_Idx)%Sds_Name), trim(temp_string_2)) - 1
         Wvl_String = trim(Sds_Info(Var_Idx)%Sds_Name(ipos1:ipos2))
         Sds_Info(Var_Idx)%Long_Name = trim(Sds_Info(Var_Idx)%Long_Name) // " " //trim(Wvl_String)
      endif

   enddo

 end subroutine SET_CH_ATTRIBUTES

!----------------------------------------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------------------------------------
end module LEVEL2_MOD

