! $Id: pixel_routines_mod.f90 3197 2019-03-23 01:10:50Z heidinger $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: snow_routines_mod.f90 (src)
!       SNOW_ROUTINES_MOD (program)
!
! PURPOSE: this module houses routines for computing some needed pixel-level arrays
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
! Public routines used in this MODULE:
!
!--------------------------------------------------------------------------------------
module SNOW_ROUTINES_MOD

 use CONSTANTS_MOD, only: Sym, Missing_Value_Int1, Missing_Value_Real4, &
                          int1, real4

 use PIXEL_COMMON_MOD, only: Ch, Sensor, Geo, Nwp_Pix, &
                             Diag_Pix_Array_1, &
                             Diag_Pix_Array_2, &
                             Diag_Pix_Array_3, &
                             Failed_IMS_Snow_Mask_Flag, &
                             Failed_Glob_Snow_Mask_Flag, &
                             Read_Snow_Mask

 implicit none

 private

 public:: COMPUTE_SNOW_CLASS, &
          COMPUTE_SNOW_CLASS_NWP, &
          COMPUTE_SNOW_CLASS_OISST

 private:: ADD_NDSI_160_SNOW, &
           ADD_NDSI_375_SNOW, &
           REMOVE_WARM_SNOW, &
           REMOVE_DARK_SNOW

 contains

!-------------------------------------------------------------------------------
!--- populate the snow_class array based on all available sources of Snow data
!--
!--- Input:
!---  NWP_Wat_Eqv_Snow_Depth - water equivalent snow depth from nwp
!---  NWP_Sea_Ice_Frac - sea ice fracion from nwp
!---  SST_Sea_Ice_Frac - sea ice fracion from sst data source
!---  Snow_Class_IMS - high resolution snow class field (highest priority)
!---  Snow_Class_Global - ESA GlobSnow products (lower priority)
!---
!--- Output:
!---  Snow_Class_Final - final classificiation
!---
!--- Symbology:
!---  1 = sym%NO_SNOW
!---  2 = sym%SEA_ICE
!---  3 = sym%SNOW
!-------------------------------------------------------------------------------
 subroutine COMPUTE_SNOW_CLASS(Snow_Class_NWP, Snow_Class_OISST, Snow_Class_IMS, &
                               Snow_Class_Glob,Land_Class,Snow_Class_Final, &
                               ABI_Use_104um_Flag)
 
   integer(kind=int1), intent(in), dimension(:,:):: Snow_Class_NWP
   integer(kind=int1), intent(in), dimension(:,:):: Snow_Class_OISST
   integer(kind=int1), intent(in), dimension(:,:):: Snow_Class_IMS
   integer(kind=int1), intent(in), dimension(:,:):: Snow_Class_Glob
   integer(kind=int1), intent(in), dimension(:,:):: Land_Class
   integer(kind=int1), intent(out), dimension(:,:):: Snow_Class_Final
   integer(kind=int1):: Hires_Snow_Flag
   logical, intent(in):: ABI_Use_104um_Flag

   Snow_Class_Final = Missing_Value_Int1
   Hires_Snow_Flag = 0

   !--- initialize to NWP
   Snow_Class_Final = Snow_Class_Nwp

   !--- High Res
   if (Read_Snow_Mask == sym%READ_SNOW_HIRES .and.     &
          Failed_IMS_Snow_Mask_flag == sym%NO) then
          Snow_Class_Final = Snow_Class_IMS
          Hires_Snow_Flag = 1
   endif

!--- maybe through away high res if missing in NWP?

   !-- GlobSnow - only use if selected, read in and HIRES not available
   if (Hires_Snow_Flag == 0) then
      if (Read_Snow_Mask == sym%READ_SNOW_GLOB .and.   &
          Failed_Glob_Snow_Mask_Flag == sym%NO) then
          Snow_Class_Final = Snow_Class_Glob
      endif
   endif

   !--- overwrite with oisst
   where(Snow_Class_OISST == sym%SEA_ICE)
       Snow_Class_Final = Snow_Class_OISST
   endwhere
       
   !-- check for consistency of land and snow masks
   where(Snow_Class_Final == sym%SNOW .and. Land_Class /= sym%LAND)
             Snow_Class_Final = sym%SEA_ICE
   endwhere
   where(Snow_Class_Final == sym%SEA_ICE .and. Land_Class == sym%LAND)
             Snow_Class_Final = sym%SNOW
   endwhere

   if (ABI_Use_104um_Flag) then
     !---- add snow under certain conditions
     if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(6) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(38) == sym%YES) then
           Snow_Class_Final = ADD_NDSI_160_SNOW(Snow_Class_Final,ch(6)%Ref_Toa,ch(1)%Ref_Toa,  &
                                            ch(38)%Emiss_Tropo,ch(38)%Bt_Toa_Std_3x3,NWP_PIX%Tsfc,Land_Class,Geo%Solzen)
     endif

     if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(20) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(38) == sym%YES) then
       Snow_Class_Final = ADD_NDSI_375_SNOW(Snow_Class_Final,ch(20)%Ref_Toa,ch(1)%Ref_Toa,  &
                                            ch(38)%Emiss_Tropo,ch(38)%Bt_Toa_Std_3x3,NWP_PIX%Tsfc,Land_Class,Geo%Solzen)
     endif

     !---- remove snow under certain conditions
     if (Sensor%Chan_On_Flag_Default(38) == sym%YES) then
       Snow_Class_Final = REMOVE_WARM_SNOW(Snow_Class_Final,ch(38)%Bt_Toa,NWP_PIX%Tsfc,Land_Class)
     endif
     if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       Snow_Class_Final = REMOVE_DARK_SNOW(Snow_Class_Final,ch(1)%Ref_Toa,Geo%Solzen,Land_Class)
     endif

   else
     !---- add snow under certain conditions
     if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(6) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(31) == sym%YES) then
         Snow_Class_Final = ADD_NDSI_160_SNOW(Snow_Class_Final,ch(6)%Ref_Toa,ch(1)%Ref_Toa,  &
                                            ch(31)%Emiss_Tropo,ch(31)%Bt_Toa_Std_3x3,NWP_PIX%Tsfc,Land_Class,Geo%Solzen)
     endif

     if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(20) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(31) == sym%YES) then
         Snow_Class_Final = ADD_NDSI_375_SNOW(Snow_Class_Final,ch(20)%Ref_Toa,ch(1)%Ref_Toa,  &
                                            ch(31)%Emiss_Tropo,ch(31)%Bt_Toa_Std_3x3,NWP_PIX%Tsfc,Land_Class,Geo%Solzen)
     endif

     !---- remove snow under certain conditions
     if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       Snow_Class_Final = REMOVE_WARM_SNOW(Snow_Class_Final,ch(31)%Bt_Toa,NWP_PIX%Tsfc,Land_Class)
     endif
     if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       Snow_Class_Final = REMOVE_DARK_SNOW(Snow_Class_Final,ch(1)%Ref_Toa,Geo%Solzen,Land_Class)
     endif

   endif

   ! --- make sure everything out of range is missing 
   where(Snow_Class_Final < sym%NO_SNOW .and. &
         Snow_Class_Final > sym%SNOW)
        Snow_Class_Final = Missing_Value_Int1
   endwhere

 end subroutine COMPUTE_SNOW_CLASS

 !---------------------------------------------------------------------------------------
 ! Compute a Snow Classification from the NWP
 ! 
 ! threshold are empirically derived by comparing NWP fields to IMS (A.  Heidinger)
 !---------------------------------------------------------------------------------------
 subroutine COMPUTE_SNOW_CLASS_NWP(NWP_Wat_Eqv_Snow_Depth,NWP_Sea_Ice_Frac, Snow_Class_Nwp)

   real(kind=real4), intent(in), dimension(:,:):: NWP_Wat_Eqv_Snow_Depth,  &
                                                  NWP_Sea_Ice_Frac
   integer(kind=int1), intent(out), dimension(:,:):: Snow_Class_NWP

   !--- initialize all pixels as no snow (GFS is now missing over ocean)
   Snow_Class_NWP = sym%NO_SNOW

   !--- detect snow  (note snow can cover sea-ice)
   where (NWP_Wat_Eqv_Snow_Depth > 0.1) 
          Snow_Class_NWP = sym%SNOW
   end where 

   !--- detect sea ice
   where (NWP_Sea_Ice_Frac > 0.5) 
          Snow_Class_NWP = sym%SEA_ICE
   end where

   !--- detect areas with missing information
   where (NWP_Wat_Eqv_Snow_Depth == Missing_Value_Real4 .and. &
          NWP_Sea_Ice_Frac == Missing_Value_Real4)
          Snow_Class_NWP = Missing_Value_Int1
   end where

 end subroutine COMPUTE_SNOW_CLASS_NWP
 !---------------------------------------------------------------------------------------
 ! Compute a Sea-Ice Classification from OISST Analysis
 ! Note =- OISST only provides Sea Ice, not Snow
 !---------------------------------------------------------------------------------------
 subroutine COMPUTE_SNOW_CLASS_OISST(SST_Sea_Ice_Frac, Snow_Class_OISST)

   real(kind=real4), intent(in), dimension(:,:):: SST_Sea_Ice_Frac
   integer(kind=int1), intent(out), dimension(:,:):: Snow_Class_OISST

   !--- initialize all to missing
   Snow_Class_OISST = Missing_Value_Int1

   !--- initialize valid values as no_snow
   where (SST_Sea_Ice_Frac >= 0.0) 
          Snow_Class_OISST = sym%NO_SNOW
   end where
   !--- detect sea ice
   where (SST_Sea_Ice_Frac > 0.5) 
          Snow_Class_OISST = sym%SEA_ICE
   end where

 end subroutine COMPUTE_SNOW_CLASS_OISST

!==============================================================================
!
! Add snow using NDSI from 1.6 and 0.65 um
!
!==============================================================================
 real elemental function ADD_NDSI_160_SNOW ( &
                                            snow_class, &
                                            toa_160_reflectance, &
                                            toa_063_reflectance, &
                                            tsfc, &
                                            emiss_tropo, &
                                            bt11_stddev, &
                                            land_class, &
                                            solar_zenith)
 

  real, intent(in):: toa_160_reflectance
  real, intent(in):: toa_063_reflectance
  real, intent(in):: tsfc
  real, intent(in):: emiss_tropo
  real, intent(in):: bt11_stddev
  real, intent(in):: solar_zenith
  integer(kind=int1), intent(in):: snow_class
  integer(kind=int1), intent(in):: land_class
  real, parameter:: solar_zenith_max_threshold = 85.0 ! orig 75.0
  real, parameter:: toa_063_reflectance_threshold = 10.0
  !real, parameter:: toa_ndsi_threshold = 0.35 ! orig 0.50 
  real, parameter:: toa_ndsi_threshold = 0.60
  real, parameter:: emiss_tropo_threshold = 0.5
  real, parameter:: bt11_stddev_threshold = 1.0  
  real, parameter:: tsfc_threshold = 277
  real:: toa_ndsi

  add_ndsi_160_snow = snow_class

  if (tsfc  > tsfc_threshold) return

  if (emiss_tropo > emiss_tropo_threshold) return

  if (bt11_stddev > bt11_stddev_threshold) return

  if (solar_zenith > solar_zenith_max_threshold) return

  if (toa_063_reflectance > toa_063_reflectance_threshold) return

  if (land_class /= sym%LAND) return

  toa_ndsi = (toa_063_reflectance - toa_160_reflectance) / (toa_063_reflectance + toa_160_reflectance)

  if (toa_ndsi > toa_ndsi_threshold) add_ndsi_160_snow = sym%SNOW

 end function ADD_NDSI_160_SNOW

!==============================================================================
 real elemental function ADD_NDSI_375_SNOW ( &
                                            snow_class, &
                                            toa_375_reflectance, &
                                            toa_063_reflectance, &
                                            tsfc, &
                                            emiss_tropo, &
                                            bt11_stddev, &
                                            land_class, &
                                            solar_zenith)
 

  real, intent(in):: toa_375_reflectance
  real, intent(in):: toa_063_reflectance
  real, intent(in):: tsfc
  real, intent(in):: emiss_tropo
  real, intent(in):: bt11_stddev
  integer(kind=int1), intent(in):: snow_class
  integer(kind=int1), intent(in):: land_class
  real, intent(in):: solar_zenith
  real, parameter:: solar_zenith_max_threshold = 85.0 ! orig 75.0
  real, parameter:: toa_063_reflectance_threshold = 10.0
  real, parameter:: toa_ndsi_threshold = 0.5  ! orig 0.75
  real, parameter:: emiss_tropo_threshold = 0.5
  real, parameter:: bt11_stddev_threshold = 1.0  
  real, parameter:: tsfc_threshold = 277
  real:: toa_ndsi

  add_ndsi_375_snow = snow_class

  if (tsfc  > tsfc_threshold) return

  if (emiss_tropo > emiss_tropo_threshold) return

  if (bt11_stddev > bt11_stddev_threshold) return

  if (solar_zenith > solar_zenith_max_threshold) return

  if (toa_063_reflectance > toa_063_reflectance_threshold) return

  if (land_class /= sym%LAND) return

  toa_ndsi = (toa_063_reflectance - toa_375_reflectance) / (toa_063_reflectance + toa_375_reflectance)

  if (toa_ndsi > toa_ndsi_threshold) add_ndsi_375_snow = sym%SNOW

 end function ADD_NDSI_375_SNOW

 !==============================================================================
 ! function to "melt" snow or ice if it is too warm
 !==============================================================================
 real elemental function REMOVE_WARM_SNOW ( &
                                            snow_class, &
                                            toa_11_bt, &
                                            surface_temperature, &
                                            land_class)
 

  integer(kind=int1), intent(in):: snow_class
  real, intent(in):: toa_11_bt
  real, intent(in):: surface_temperature
  integer(kind=int1), intent(in):: land_class
  real, parameter:: toa_11_bt_max_threshold = 290.0  ! orig 277.0
  real, parameter:: tsfc_max_threshold = 320.0  ! orig 277.0

  remove_warm_snow = snow_class

  if ((toa_11_bt > toa_11_bt_max_threshold .or. &
      surface_temperature > tsfc_max_threshold) .and. &
      land_class == sym%LAND) then
         remove_warm_snow = sym%NO_SNOW
  endif

 end function REMOVE_WARM_SNOW

 !==============================================================================
 ! function to "melt" snow or ice if it is too dark
 !==============================================================================
 real elemental function REMOVE_DARK_SNOW ( &
                                            snow_class, &
                                            toa_063_reflectance, &
                                            solar_zenith, &
                                            land_class)
  integer(kind=int1), intent(in):: snow_class
  real, intent(in):: toa_063_reflectance
  real, intent(in):: solar_zenith
  integer(kind=int1), intent(in):: land_class
  real, parameter:: solar_zenith_max_threshold = 85.0 !75.0
  real, parameter:: toa_063_reflectance_threshold = 10.0

  remove_dark_snow = snow_class

  if (solar_zenith < solar_zenith_max_threshold .and. &
      toa_063_reflectance < toa_063_reflectance_threshold .and. &
      land_class == sym%LAND) then
         remove_dark_snow = sym%NO_SNOW
  endif

 end function REMOVE_DARK_SNOW

!-----------------------------------------------------------
! end of MODULE
!-----------------------------------------------------------
end module SNOW_ROUTINES_MOD
