! $Id: clavrx_sst_mod.f90 3862 2020-06-03 18:50:40Z heidinger $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: clavrx_sst_module.f90 (src)
!       clavrx_sst_module.f90 (program)
!
! PURPOSE:
!       compute Sea Surface Temperature (SST)
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
! File I/O: None, note this a clavrx module makes use of the clavrx global
!           memory in PIXEL_COMMON_MOD and CONSTANTS
!
! Public Routines:
!  COMPUTE_SST
!
! Private Routines:
!  MCSST
!
!--------------------------------------------------------------------------------------
module CLAVRX_SST_MOD
use CONSTANTS_MOD
use PIXEL_COMMON_MOD, only: ch, sensor, geo, sfc, Bad_Pixel_Mask, Sst_Retrieved,  &
                        CLDMASK, Sst_Anal, NWP_Pix, Temp_Pix_Array_1, Use_Sst_Anal
use ALGORITHM_CONSTANTS_MOD

implicit none
private

public:: SETUP_SST, COMPUTE_SST
private::  MCSST

real, dimension(4), save, private:: sst_coef
real, dimension(8), save, private:: acspo_coef

contains

!-----------------------------------------------------------------
!  set sst_coef vector for the sensor
!
!-----------------------------------------------------------------
subroutine SETUP_SST()


!--- initialize
sst_coef = MISSING_VALUE_REAL4
acspo_coef = MISSING_VALUE_REAL4

!---
select case (Sensor%WMO_Id)

      case(4) !METOP-A
         sst_coef = (/ 1.0159, 2.2285, 1.2289, -277.2571/)

      case(3) !METOP-B
         sst_coef = (/ 1.0148, 1.7483, 1.2989, -276.8420/)

      case(5) !METOP-C

      case(55) !MSG-8

      case(56) !MSG-9

      case(57) !MSG-10

      case(70) !MSG-11

      case(171) !MTSAT-1R

      case(172) !MTSAT-2

      case(173) !AHI-8
!         acspo_coef = (/38.112282,0.876435,0.826084,0.511226,0.5617733,0.020591,0.034769,0.078092 /)
         acspo_coef = (/14.507927,0.963572,1.428572,0.641366,1.214019,0.076626,0.070244,0.098831 /)

      case(174) !AHI-9, but values are replicated from AHI-8 above.
         acspo_coef = (/38.112282,0.876435,0.826084,0.511226,0.5617733,0.020591,0.034769,0.078092 /)

      case(200) !NOAA-8

      case(201) !NOAA-9
         sst_coef = (/0.9639, 3.0260, 0.0539, -262.4647 /)

      case(202) !NOAA-10

      case(203) !NOAA-11
         sst_coef = (/0.9762,    2.1872,    1.1063, -265.5705/)

      case(204) !NOAA-12
         sst_coef = (/0.9972,    2.2811,    0.9336, -271.8081/)

      case(205) !NOAA-14
         sst_coef = (/0.9865,    2.2326,    0.8566, -269.5311/)

      case(206) !NOAA-15
         sst_coef = (/0.964243, 2.71296, 0.387491, -262.443/)

      case(207) !NOAA-16
         sst_coef = (/0.9687, 2.3471,  0.8959, -264.7613/)

      case(208) !NOAA-17
         sst_coef = (/0.992818, 2.49916, 0.915103, -271.206/)

      case(209) !NOAA-18
         sst_coef = (/ 1.012, 1.817,  1.171, -276.079/)

      case(223) !NOAA-19
         sst_coef = (/1.037, 1.574, 0.984, -283.552/)

      case(224,225,226) !VIIRS
         sst_coef = (/1.0563,    1.5878,    2.0775, -287.0832 /)

      case(252) !GOES-8

      case(253) !GOES-9

      case(254) !GOES-10

      case(255) !GOES-11

      case(256) !GOES-12

      case(257) !GOES-13

      case(258) !GOES-14

      case(259) !GOES-15

      case(270) !GOES-16
         acspo_coef = (/61.353560,0.795478,2.23466,0.4403643,1.1701886,0.070138,0.177908,0.185610/)

      case(271) !GOES-17
         acspo_coef = (/74.574503,0.751711,2.81401,0.2970021,1.3209088,0.014818,0.108972,0.153607/)

      case(272) !GOES-18 TODO FAKE NUMBER ADJUST WHEN KNOWN
         acspo_coef = (/74.574503,0.751711,2.81401,0.2970021,1.3209088,0.014818,0.108972,0.153607/)



      case(706) !NOAA-6

      case(707) !NOAA-7
         sst_coef = (/0.9822, 2.2534, 0.9260, -267.126/)

      case(708) !NOAA-5

      case(783) !MODIS
         sst_coef = (/1.0557,    1.3826,    3.8160, -287.1476/)

      case(784) !MODIS
         sst_coef = (/1.0252,    3.1114,    1.9334, -279.5697/)

      case(810) !COMS

      case(840) !EPS-SG 

      case(514) !FY2D

      case(515) !FY2E

      case(523) !FY3D

      case(530) !FY-4A

      case default
         print*,'sensor for WMO number not found in SST  ', Sensor%WMO_id
         print*,'SST will not be computed'
      end select

end subroutine SETUP_SST
!-----------------------------------------------------------------
!
!-----------------------------------------------------------------
subroutine COMPUTE_SST()

      Sst_Retrieved = MISSING_VALUE_REAL4

      !--- this the background sst (if needed) and initialize with nwp
      Temp_Pix_Array_1 = NWP_Pix%Tsfc

      !--- if an Sst_Anal is available, use it for background
      if (Use_Sst_Anal == sym%YES) then
         Temp_Pix_Array_1 = Sst_Anal
         !--- fill in holes with nwp
         where(Sst_Anal == MISSING_VALUE_REAL4)
           Temp_Pix_Array_1 = NWP_Pix%Tsfc
         endwhere
      endif

      if ((Sensor%Wmo_Id == 173 .or. Sensor%Wmo_Id == 174 .or. &
           Sensor%Wmo_Id == 270 .or. Sensor%Wmo_Id == 271) .and. &
          Sensor%Chan_On_Flag_Default(29)==sym%YES .and. &
          Sensor%Chan_On_Flag_Default(31)==sym%YES .and. &
          Sensor%Chan_On_Flag_Default(32)==sym%YES .and. &
          Sensor%Chan_On_Flag_Default(38)==sym%YES) then

          Sst_Retrieved = ACSPO_SST(ch(29)%Bt_Toa,ch(38)%Bt_Toa,ch(31)%Bt_Toa,ch(32)%Bt_Toa,Temp_Pix_Array_1,Geo%Seczen)

      elseif (Sensor%Chan_On_Flag_Default(31)==sym%YES .and. &
          Sensor%Chan_On_Flag_Default(32)==sym%YES .and. &
          maxval(sst_coef) /= MISSING_VALUE_REAL4) then

          Sst_Retrieved = MCSST(ch(31)%Bt_Toa,ch(32)%Bt_Toa,Geo%Seczen)

      endif



      !--- mask bad pixels, land and snow/ice
      where(Bad_Pixel_Mask == sym%YES .or.  &
            Sfc%Land == sym%Land .or. Sfc%Land == sym%COASTLINE .or. &
            Sfc%Snow /= sym%NO_SNOW)

              Sst_Retrieved = MISSING_VALUE_REAL4

      endwhere


end subroutine COMPUTE_SST

!----------------------------------------------------------------------
! function to compute sst from mcsst formulation
!----------------------------------------------------------------------
real elemental function mcsst(bt11,bt12,seczen)

  real, intent(in):: bt11
  real, intent(in):: bt12
  real, intent(in):: seczen

  !this makes the window eff flux temp (K) (win_Olr = sigma*win_temp^4)

  mcsst =  sst_coef(1) * bt11 +  &
           sst_coef(2) * (bt11-bt12) +  &
           sst_coef(3) * (bt11-bt12) * (seczen-1.0) +  &
           sst_coef(4) + 273.15

end function mcsst

!----------------------------------------------------------------------
! function to compute sst from acspo formulation
!----------------------------------------------------------------------
real elemental function acspo_sst(bt85,bt104,bt11,bt12,first_guess_sst,seczen)

  real, intent(in):: bt85
  real, intent(in):: bt104
  real, intent(in):: bt11
  real, intent(in):: bt12
  real, intent(in):: first_guess_sst
  real, intent(in):: seczen

  acspo_sst =  acspo_coef(1) +  &
               acspo_coef(2) * (bt104) +  &
               (-1.0)*acspo_coef(3) * (bt104 - bt12) +  &
               acspo_coef(4) * (bt104 - bt85)*seczen + &
               acspo_coef(5) * (bt104 - bt11)*seczen + &
               (-1.0)*acspo_coef(6) * (bt104 - bt85)*(first_guess_sst - 273.15) + &
               (-1.0)*acspo_coef(7) * (bt104 - bt11)*(first_guess_sst - 273.15) + &
               1.0*acspo_coef(8) * (bt104 - bt12)*(first_guess_sst - 273.15)

end function acspo_sst

end module CLAVRX_SST_MOD
