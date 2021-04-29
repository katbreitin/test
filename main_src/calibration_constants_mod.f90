!$Id: calibration_constants_mod.f90 3804 2020-04-21 15:14:30Z stevew $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: calibration_constants.f90 (src)
!       CALIBRATION_CONSTANTS (program)
!
! PURPOSE: This module serves as a common block for passing the 
!          instrument and calibration coefficients
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
!   File I/O - none
!
!   Public routines: none
!
!   Private routines: none
!--------------------------------------------------------------------------------------
module CALIBRATION_CONSTANTS_MOD

  use CONSTANTS_MOD, only: &
  real4 & 
  , real8 &
  , nchan_clavrx

 implicit none
 private

 public:: MITIGATE_CH20_NOISE

 !--- Planck Constants
 real (kind=real4), parameter, public:: c1 = 1.191062e-5, &
                                        c2 = 1.4387863

 real (kind=real4), dimension(20:Nchan_Clavrx), public, save:: Planck_A1 = 0.0
 real (kind=real4), dimension(20:Nchan_Clavrx), public, save:: Planck_A2 = 1.0

 real (kind=real4), dimension(1:Nchan_Clavrx), public, save:: Planck_Nu = 0.0
 real (kind=real4), public, save :: Planck_A1_375um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A2_375um_Sndr = 1.0
 real (kind=real4), public, save :: Planck_Nu_375um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A1_11um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A2_11um_Sndr = 1.0
 real (kind=real4), public, save :: Planck_Nu_11um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A1_12um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A2_12um_Sndr = 1.0
 real (kind=real4), public, save :: Planck_Nu_12um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A1_145um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A2_145um_Sndr = 1.0
 real (kind=real4), public, save :: Planck_Nu_145um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A1_147um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A2_147um_Sndr = 1.0
 real (kind=real4), public, save :: Planck_Nu_147um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A1_457um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A2_457um_Sndr = 1.0
 real (kind=real4), public, save :: Planck_Nu_457um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A1_445um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A2_445um_Sndr = 1.0
 real (kind=real4), public, save :: Planck_Nu_445um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A1_149um_Sndr = 0.0
 real (kind=real4), public, save :: Planck_A2_149um_Sndr = 1.0
 real (kind=real4), public, save :: Planck_Nu_149um_Sndr = 0.0


 !----- variables read in from instrument constant files (planck stored in constants module)
 real (kind=real4), save, public:: B0_3b = 0.0
 real (kind=real4), save, public:: B1_3b = 0.0
 real (kind=real4), save, public:: B2_3b = 0.0
 real (kind=real4), save, public:: B0_4 = 0.0
 real (kind=real4), save, public:: B1_4 = 0.0
 real (kind=real4), save, public:: B2_4 = 0.0
 real (kind=real4), save, public:: B0_5 = 0.0
 real (kind=real4), save, public:: B1_5 = 0.0
 real (kind=real4), save, public:: B2_5 = 0.0
 real (kind=real4), save, public:: Space_Rad_3b = 0.0
 real (kind=real4), save, public:: Space_Rad_4 = 0.0
 real (kind=real4), save, public:: Space_Rad_5 = 0.0

 real(kind=real4), dimension(0:4,4),public,save:: Prt_Coef = 0.0
 real(kind=real4), dimension(4),public,save:: Prt_Weight = 0.0

 real(kind=real4),save,public:: Ch1_Gain_Low = 0.0
 real(kind=real4),save,public:: Ch1_Gain_High = 0.0
 real(kind=real4),save,public:: Ch1_Gain_Low_0 = 0.0
 real(kind=real4),save,public:: Ch1_Gain_High_0 = 0.0
 real(kind=real4),save,public:: Ch1_Degrad_Low_1 = 0.0
 real(kind=real4),save,public:: Ch1_Degrad_High_1 = 0.0
 real(kind=real4),save,public:: Ch1_Degrad_Low_2 = 0.0
 real(kind=real4),save,public:: Ch1_Degrad_High_2 = 0.0
 real(kind=real4),save,public:: Ch1_Dark_Count = 0.0
 real(kind=real4),save,public:: Ch1_Dark_Count_Cal = 0.0
 real(kind=real4),save,public:: Ch1_Switch_Count = 0.0
 real(kind=real4),save,public:: Ch1_Switch_Count_Cal = 0.0
 real(kind=real4),save,public:: Ref_Ch1_Switch = 0.0

 real(kind=real4),save,public:: Ch2_Gain_Low = 0.0
 real(kind=real4),save,public:: Ch2_Gain_High = 0.0
 real(kind=real4),save,public:: Ch2_Gain_Low_0 = 0.0
 real(kind=real4),save,public:: Ch2_Gain_High_0 = 0.0
 real(kind=real4),save,public:: Ch2_Degrad_Low_1 = 0.0
 real(kind=real4),save,public:: Ch2_Degrad_High_1 = 0.0
 real(kind=real4),save,public:: Ch2_Degrad_Low_2 = 0.0
 real(kind=real4),save,public:: Ch2_Degrad_High_2 = 0.0
 real(kind=real4),save,public:: Ch2_Dark_Count = 0.0
 real(kind=real4),save,public:: Ch2_Dark_Count_Cal = 0.0
 real(kind=real4),save,public:: Ch2_Switch_Count = 0.0
 real(kind=real4),save,public:: Ch2_Switch_Count_Cal = 0.0
 real(kind=real4),save,public:: Ref_Ch2_Switch = 0.0

 real(kind=real4),save,public:: Ch3a_Gain_Low = 0.0
 real(kind=real4),save,public:: Ch3a_Gain_High = 0.0
 real(kind=real4),save,public:: Ch3a_Gain_Low_0 = 0.0
 real(kind=real4),save,public:: Ch3a_Gain_High_0 = 0.0
 real(kind=real4),save,public:: Ch3a_Degrad_Low_1 = 0.0
 real(kind=real4),save,public:: Ch3a_Degrad_High_1 = 0.0
 real(kind=real4),save,public:: Ch3a_Degrad_Low_2 = 0.0
 real(kind=real4),save,public:: Ch3a_Degrad_High_2 = 0.0
 real(kind=real4),save,public:: Ch3a_Dark_Count = 0.0
 real(kind=real4),save,public:: Ch3a_Dark_Count_Cal = 0.0
 real(kind=real4),save,public:: Ch3a_Switch_Count = 0.0
 real(kind=real4),save,public:: Ch3a_Switch_Count_Cal = 0.0
 real(kind=real4),save,public:: Ref_Ch6_Switch = 0.0

 real(kind=real4),save,public:: Sun_Earth_Distance = 0.0
 real(kind=real4),save,public:: Launch_Date = 0.0
 
 real(kind=real8),save,public:: Goes_Input_Time = 0
 real(kind=real8),save,public:: Goes_Epoch_Time = 0

 real(kind=real4),save,public:: Solar_Ch20 = 0.0
 real(kind=real4),save,public:: Solar_Ch20_Nu = 0.0
 real(kind=real4),save,public:: Ew_Ch20 = 0.0
 character(len=7),save,public:: Sat_Name = ' '

 !--- ABI specific constants
 real(kind=real4),save,public:: Band2_Correction_Start_Date    !GOES16/17 only
 real(kind=real4),save,public:: Band2_Correction_Factor        !GOES16/17 only

 !--- Initialize ABI Focal Plane Temperature objective thresholds from Schmit/Gunshor.
 !--- Values are now read in through the instrument file.
 real(kind=real4),save,public:: ABI_FPT_Thresh_038um = 0.0
 real(kind=real4),save,public:: ABI_FPT_Thresh_062um = 0.0
 real(kind=real4),save,public:: ABI_FPT_Thresh_067um = 0.0
 real(kind=real4),save,public:: ABI_FPT_Thresh_073um = 0.0
 real(kind=real4),save,public:: ABI_FPT_Thresh_085um = 0.0
 real(kind=real4),save,public:: ABI_FPT_Thresh_097um = 0.0
 real(kind=real4),save,public:: ABI_FPT_Thresh_104um = 0.0
 real(kind=real4),save,public:: ABI_FPT_Thresh_110um = 0.0
 real(kind=real4),save,public:: ABI_FPT_Thresh_120um = 0.0
 real(kind=real4),save,public:: ABI_FPT_Thresh_133um = 0.0

 !--- VIIRS specific constants
 real(kind=real4), dimension(11), save, public:: VIIRS_Correction_Factor

 !--- Variables needed for static navigation: conversions from Rad to Ref
 !--- AHI + ABI
 real (kind=real4), save, public:: Rad_to_Ref_Fac_0_47um = 0.0 ! ABI_1, AHI_1
 real (kind=real4), save, public:: Rad_to_Ref_Fac_0_55um = 0.0 ! AHI_2
 real (kind=real4), save, public:: Rad_to_Ref_Fac_0_65um = 0.0 ! ABI_2, AHI_3
 real (kind=real4), save, public:: Rad_to_Ref_Fac_0_86um = 0.0 ! ABI_3, AHI_4
 real (kind=real4), save, public:: Rad_to_Ref_Fac_1_60um = 0.0 ! ABI_5, AHI_5
 real (kind=real4), save, public:: Rad_to_Ref_Fac_2_10um = 0.0 ! ABI_6, AHI_6
 real (kind=real4), save, public:: Rad_to_Ref_Fac_1_38um = 0.0 ! ABI_4

 contains

!-------------------------------------------------------------------------
! NOAA-5 = 708, NOAA-6 = 706, NOAA-7 = 707
! NOAA-8 = 200, NOAA-9 = 201, NOAA-10 = 202, NOAA-11=203, NOAA-12=204
! threshold is assumed to 270K - should be refined to be sensor dependent
!-------------------------------------------------------------------------
 subroutine MITIGATE_CH20_NOISE(WMO_Id,Bt20,Bt20_Med,Bt11)
   integer, intent(in):: WMO_Id
   real, intent(inout), dimension(:,:):: Bt20
   real, intent(in), dimension(:,:):: Bt20_Med
   real, intent(in), dimension(:,:):: Bt11

   !--- handle noise in imager values
   select case (WMO_Id)
     case (200,201,202,706:708) 
       where(Bt11 < 300.0)
          Bt20 = Bt20_Med
       endwhere
     case default
   end select

 end subroutine MITIGATE_CH20_NOISE

end module CALIBRATION_CONSTANTS_MOD
