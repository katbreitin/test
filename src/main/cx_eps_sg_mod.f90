! $Id:$
!----------------------------------------------------------------------
! - CLAVR-x EUMETSAT Polar System - Secong Generation (EPS-SG)
!----------------------------------------------------------------------
!
! Ch #     Wavelength     MODIS Ch.  Obs_Type
!  1         0.443           3        solar
!  2         0.555           4        solar
!  3         0.668           1        solar
!  4         0.752          15        solar
!  5         0.763          45        solar
!  6         0.865           2        solar
!  7         0.914          17        solar
!  8         1.240           5        solar
!  9         1.375          26        solar
! 10         1.630           6        solar
! 11         2.250           7        solar
! 12         3.740          20        mixed
! 13         3.959          21        mixed
! 14         4.050          23        therm
! 15         6.725          27        therm
! 16         7.325          28        therm
! 17         8.540          29        therm
! 18        10.690          31        therm
! 19        12.020          32        therm
! 20        13.345          33        therm
!----------------------------------------------------------------------


module CX_EPS_SG_MOD

!--- use statements
use PIXEL_COMMON_MOD, only: Ch, Geo, Image, Nav, Sensor,  &
                            Temp_Pix_Array_1, Gap_Pixel_Mask

use CX_NETCDF4_MOD,only: OPEN_NETCDF, &
                         CLOSE_NETCDF, &
                         GET_GROUP_ID, &
                         READ_NETCDF_GLOBAL_ATTRIBUTE, &
                         READ_NETCDF_DIMENSION_2D, &
                         READ_NETCDF, &
                         READ_AND_UNSCALE_NETCDF_1D, &
                         READ_AND_UNSCALE_NETCDF_2D &
                         , read_abi_max_focal_plane_temp

use VIEWING_GEOMETRY_MOD, only: GLINT_ANGLE, SCATTERING_ANGLE, &
                                RELATIVE_AZIMUTH

use PLANCK_MOD,only: COMPUTE_BT_ARRAY  &
              , CONVERT_RADIANCE

use CONSTANTS_MOD, only: MSEC_PER_DAY, Sym, MISSING_VALUE_REAL4, &
                         MISSING_VALUE_REAL8,  &
                         SOLAR_OBS_TYPE, THERMAL_OBS_TYPE, &
                         MIXED_OBS_TYPE, LUNAR_OBS_TYPE, DTOR

use CALIBRATION_CONSTANTS_MOD,only: &
        Planck_A1 &
      , Planck_A2 &
      , Planck_Nu &
      , ew_ch20 &
      , sat_name &     ! overwritten here
      , solar_ch20 &  ! overwritten here
      , solar_ch20_nu ! overwritten here

use CLAVRX_MESSAGE_MOD,only: MESG, verb_lev

implicit none

private :: CONVERT_RAD_TO_REF
public :: READ_NUMBER_OF_SCANS_EPS_SG
public :: READ_EPS_SG_DATE_TIME
public :: READ_EPS_SG_INSTR_CONSTANTS
public:: READ_EPS_SG_DATA

!--- module wide variables and parameters
integer, private, save :: Ncid_Eps_Sg
real, dimension(11), private, save :: Solar_Corr_Coeff
character(len=14), parameter :: EXE_PROMPT="EPS_SG_MODULE:"

type eps_navigation
  logical :: is_calculated  = .FALSE.
  real, allocatable, dimension(:,:) :: longitude
  real, allocatable, dimension(:,:) :: latitude
  real, allocatable, dimension(:,:) :: observation_azimuth
  real, allocatable, dimension(:,:) :: observation_zenith
  real, allocatable, dimension(:,:) :: solar_azimuth
  real, allocatable, dimension(:,:) :: solar_zenith


  contains
  procedure :: set_nav => eps_navigation_set_nav
  procedure :: dealloc => eps_navigation_dealloc
end type eps_navigation

type ( eps_navigation) :: eps_nav

contains

!----------------------------------------------------------------------
! open and read dimensions
!----------------------------------------------------------------------
subroutine READ_NUMBER_OF_SCANS_EPS_SG(File_Name,Nscn,Npix,Ierror)
   character(len=*), intent(in):: File_Name
   integer, intent(out):: Nscn,Npix,Ierror

   integer:: Group_Id1,Group_Id2
   integer, dimension(2):: Dim_2d

   ! --- initialize
   Ierror = 0
   Dim_2d = -1

   ! --- open level1b file
   call OPEN_NETCDF(File_Name,Ncid_Eps_Sg)

   ! --- get group id
   call GET_GROUP_ID(Ncid_Eps_Sg, 'data', Group_Id1)
   call GET_GROUP_ID(Group_Id1, 'measurement_data', Group_Id2)

   ! --- read dimension
   call READ_NETCDF_DIMENSION_2D(Group_Id2,'delta_lat_N_dem', Dim_2d)
   Nscn = Dim_2d(2)
   Npix = Dim_2d(1)

end subroutine READ_NUMBER_OF_SCANS_EPS_SG

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine READ_EPS_SG_DATE_TIME(File_Name,Start_Year,Start_Doy,Start_Time,&
                               End_Year,End_Doy,End_Time)

   use CX_DATE_TIME_TOOLS_MOD, only: JULIAN


   character(len=*), intent(in):: File_Name

   integer, intent(out) :: Start_Year  !year
   integer, intent(out) :: Start_Doy   !day of year
   integer, intent(out) :: Start_Time  !millisec
   integer, intent(out) :: End_Year    !year
   integer, intent(out) :: End_Doy     !day of year
   integer, intent(out) :: End_Time    !millisec

   character(len=25) :: String_Tmp
   integer :: Start_Month, Start_Day, Start_Hour, Start_Minute, Start_Sec
   integer :: End_Month, End_Day, End_Hour, End_Minute, End_Sec
   integer :: Start_Msec, End_Msec
   !integer :: resloc, tempstrlen

   ! --- read start time
   call READ_NETCDF_GLOBAL_ATTRIBUTE(File_Name, 'sensing_start_time_utc', String_Tmp)
   read(String_Tmp(1:4), fmt="(I4)") Start_Year
   read(String_Tmp(6:7), fmt="(I2)") Start_Month
   read(String_Tmp(9:10), fmt="(I2)") Start_Day

!+++++++++++++++++++++++++ new
   !resloc = index(File_Name, 'W_')
   !tempstrlen = len('W_xx-eumetsat-darmstadt,SAT,SGA1-VII-1B-RAD_C_EUMT_20191001043852_G_D_')
   !read(File_Name(resloc+tempstrlen:resloc+tempstrlen+3), fmt="(I4)") Start_Year
   !read(File_Name(resloc+tempstrlen+4:resloc+tempstrlen+5), fmt="(I2)") Start_Month
   !read(File_Name(resloc+tempstrlen+6:resloc+tempstrlen+7), fmt="(I2)") Start_Day
!+++++++++++++++++++++++++

   !Start_Month
   !Start_Day
!print *, 'eps time ',File_Name, 'result ',resloc, String_Tmp, Start_Year,Start_Month,Start_Day
!stop
   !--- compute day of year
   call JULIAN (Start_Day, Start_Month, Start_Year, Start_Doy)

   read(String_Tmp(12:13), fmt="(I2)") Start_Hour
   read(String_Tmp(15:16), fmt="(I2)") Start_Minute
   read(String_Tmp(18:19), fmt="(I2)") Start_Sec
   read(String_Tmp(21:23), fmt="(I3)") Start_Msec

   !read(File_Name(resloc+tempstrlen+8:resloc+tempstrlen+9), fmt="(I2)") Start_Hour
   !read(File_Name(resloc+tempstrlen+10:resloc+tempstrlen+11), fmt="(I2)") Start_Minute
   !read(File_Name(resloc+tempstrlen+12:resloc+tempstrlen+13), fmt="(I2)") Start_Sec
   !Start_Msec = 0

   ! --- read end time
   call READ_NETCDF_GLOBAL_ATTRIBUTE(File_Name, 'sensing_end_time_utc', String_Tmp)

   read(String_Tmp(1:4), fmt="(I4)") End_Year
   read(String_Tmp(6:7), fmt="(I2)") End_Month
   read(String_Tmp(9:10), fmt="(I2)") End_Day

!+++++++++++++++++++++++++ new
   !read(File_Name(resloc+tempstrlen+15:resloc+tempstrlen+18), fmt="(I4)") End_Year
   !read(File_Name(resloc+tempstrlen+19:resloc+tempstrlen+20), fmt="(I2)") End_Month
   !read(File_Name(resloc+tempstrlen+21:resloc+tempstrlen+22), fmt="(I2)") End_Day
!+++++++++++++++++++++++++

   !Start_Month
   !Start_Day

   !--- compute day of year
   call JULIAN (End_Day, End_Month, End_Year, End_Doy)

   read(String_Tmp(12:13), fmt="(I2)") End_Hour
   read(String_Tmp(15:16), fmt="(I2)") End_Minute
   read(String_Tmp(18:19), fmt="(I2)") End_Sec
   read(String_Tmp(21:23), fmt="(I3)") End_Msec

   !read(File_Name(resloc+tempstrlen+23:resloc+tempstrlen+24), fmt="(I2)") End_Hour
   !read(File_Name(resloc+tempstrlen+25:resloc+tempstrlen+26), fmt="(I2)") End_Minute
   !read(File_Name(resloc+tempstrlen+27:resloc+tempstrlen+28), fmt="(I2)") End_Sec
   !End_Msec = 0

   ! --- Calculate start and end time
   Start_Time = ((Start_Hour * 60 + Start_Minute) * 60 + Start_Sec) * 1000 + Start_Msec
   End_Time = ((End_Hour * 60 + End_Minute) * 60 + End_Sec) * 1000 + End_Msec


end subroutine READ_EPS_SG_DATE_TIME

!----------------------------------------------------------------------
! read instrument file
!----------------------------------------------------------------------
subroutine READ_EPS_SG_INSTR_CONSTANTS(Instr_Const_File)

   use FILE_UTILS, only: GET_LUN

   character(len=*), intent(in):: Instr_Const_file

   integer:: Ios0, Erstat
   integer:: Instr_Const_Lun
   character(len=20):: header

   Instr_Const_Lun = GET_LUN()

   open(unit=Instr_Const_Lun,file=trim(Instr_Const_File),status="old",position="rewind",action="read",iostat=ios0)
   call MESG (EXE_PROMPT//" Opening "// trim(Instr_Const_File),level = verb_lev % DEFAULT)
   Erstat = 0
   if (Ios0 /= 0) then
      Erstat = 19
      call MESG (EXE_PROMPT//" Error Opening FUSION Constants File",level = verb_lev % DEFAULT)
      stop 19
   endif

   read(unit=Instr_Const_Lun,fmt="(a5)") sat_name
   read(unit=Instr_Const_Lun,fmt=*) Solar_Ch20
   read(unit=Instr_Const_Lun,fmt=*) Ew_Ch20
   read(unit=Instr_Const_lun,fmt=*) header
   read(unit=Instr_Const_lun,fmt=*) header
   read(unit=Instr_Const_Lun,fmt=*) Planck_A1(20), Planck_A2(20), Planck_Nu(20)
   read(unit=Instr_Const_Lun,fmt=*) Planck_A1(22), Planck_A2(22), Planck_Nu(22)
   read(unit=Instr_Const_Lun,fmt=*) Planck_A1(23), Planck_A2(23), Planck_Nu(23)
   read(unit=Instr_Const_Lun,fmt=*) Planck_A1(27), Planck_A2(27), Planck_Nu(27)
   read(unit=Instr_Const_Lun,fmt=*) Planck_A1(28), Planck_A2(28), Planck_Nu(28)
   read(unit=Instr_Const_Lun,fmt=*) Planck_A1(29), Planck_A2(29), Planck_Nu(29)
   read(unit=Instr_Const_Lun,fmt=*) Planck_A1(31), Planck_A2(31), Planck_Nu(31)
   read(unit=Instr_Const_Lun,fmt=*) Planck_A1(32), Planck_A2(32), Planck_Nu(32)
   read(unit=Instr_Const_Lun,fmt=*) Planck_A1(33), Planck_A2(33), Planck_Nu(33)
   read(unit=Instr_Const_lun,fmt=*) header
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(1)
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(2)
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(3)
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(4)
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(5)
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(6)
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(7)
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(8)
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(9)
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(10)
   read(unit=Instr_Const_lun,fmt=*) Solar_Corr_Coeff(11)
   close(unit=Instr_Const_Lun)

   !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
   Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

end subroutine READ_EPS_SG_INSTR_CONSTANTS

!----------------------------------------------------------------------
! read level1b
!----------------------------------------------------------------------
subroutine READ_EPS_SG_DATA(Segment_Number, Error_Status)
  use univ_kind_defs_mod, only: f8, i8, f4

  integer, intent(in):: Segment_Number
  integer, intent(out):: Error_Status

  integer, dimension(2) :: Sds_Start_2d
  integer, dimension(2) :: Sds_Stride_2d
  integer, dimension(2) :: Sds_Count_2d
  integer, dimension(1) :: Sds_Start_1d
  integer, dimension(1) :: Sds_Stride_1d
  integer, dimension(1) :: Sds_Count_1d
  integer :: Group_Id1,Group_Id2

  real(kind=f8), dimension(:), allocatable:: Sds_Data_1d
  real(kind=f4), dimension(:), allocatable:: Sds_Data_1d_real4
  real, dimension(:,:), allocatable :: Sds_Data_2d
  real, dimension(:,:), allocatable :: Sds_Data_2d_clavrx_grid
  real(kind=f8) :: Earth_Sun_Ratio
  integer, dimension (:), allocatable :: Time_Msec_Day
  integer(kind=i8), parameter :: Microsec_Per_Day =  86400000000_i8
  integer, parameter :: Num_Ch = 20
  integer, parameter :: Num_Sol_Ch = 11
  integer :: Num_Scan
  integer :: k
  integer :: x_count, y_count
  integer :: Nx_Start, Nx_End, Ny_Start, Ny_End
  integer :: Sds_Nx_Start, Sds_Nx_End, Sds_Ny_Start, Sds_Ny_End
  integer :: Chan_Idx, CLAVRx_Chan_Idx, Line_Idx, Elem_Idx
  integer :: Tmp_Start
  real, dimension(Num_Sol_Ch) :: Ch_Aver_Sol_Irrad
  character(len=5), dimension(20), parameter :: Ch_Name = &
   ["443  ","555  ","668  ","752  ","763  ","865  ","914  ", &
    "1240 ","1375 ","1630 ","2250 ","3740 ","3959 ","4050 ", &
    "6725 ","7325 ","8540 ","10690","12020","13345"]


  if (Segment_Number == 1) then
    call MESG (EXE_PROMPT//" Reading level1b",level = verb_lev % DEFAULT)
  end if

  ! --- determine location of segement data in the file
  ! -- dimensions in file are switched...
  ! -- file specific values start with sds
  Sds_Ny_Start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
  Sds_Ny_End = min(Image%Number_Of_Lines, sds_Ny_Start + Image%Number_of_Lines_Per_Segment - 1)
  Sds_Nx_Start = 1
  Sds_Nx_End = sds_Nx_Start + Image%Number_Of_Elements - 1

  Sds_Start_2d = [Sds_Nx_start,Sds_Ny_Start]
  Sds_Stride_2d = [1,1]
  Sds_Count_2d = [Sds_Nx_End - Sds_Nx_Start + 1, Sds_Ny_End - Sds_Ny_Start + 1]


  ! - data from anchor calculations

  ! - make calculations from anchor points for whole file
  call eps_nav % set_nav (sds_nx_start &
              , sds_ny_start &
              , sds_nx_end - sds_nx_start + 1 &
              , sds_ny_end - sds_ny_start + 1)


  Nav%Lon_1b(1:Sds_Count_2d(1),1:Sds_Count_2d(2))  = eps_nav%longitude
  Nav%Lat_1b(1:Sds_Count_2d(1),1:Sds_Count_2d(2))  = eps_nav%latitude
  geo % solaz(1:Sds_Count_2d(1),1:Sds_Count_2d(2)) = eps_nav%solar_azimuth
  Geo%Solzen(1:Sds_Count_2d(1),1:Sds_Count_2d(2))  = eps_nav%solar_zenith
  Geo%Satzen(1:Sds_Count_2d(1),1:Sds_Count_2d(2))   = eps_nav%observation_zenith
  where (Geo%Satzen(1:Sds_Count_2d(1),1:Sds_Count_2d(2)) < 0.0000001)
    Geo%Satzen(1:Sds_Count_2d(1),1:Sds_Count_2d(2)) = 0.0000001
  end where
  Geo%Sataz(1:Sds_Count_2d(1),1:Sds_Count_2d(2))  = eps_nav%observation_azimuth


  !--- make sure lat/lon are within the limits
!  where(Nav%Lon_1b <= -180.0 .or. Nav%Lon_1b >= 180.0)
  where(Nav%Lon_1b < -180.0 .or. Nav%Lon_1b > 180.0)
     Nav%Lon_1b = MISSING_VALUE_REAL4
  end where
!  where(Nav%Lat_1b <= -90.0 .or. Nav%Lat_1b >= 90.0)
  where(Nav%Lat_1b < -90.0 .or. Nav%Lat_1b > 90.0)
     Nav%Lat_1b = MISSING_VALUE_REAL4
  end where

  ! free  anchor values if last segment
    call eps_nav % dealloc

  !--- read earth-sun distance ratio
  allocate(Sds_Data_1d(1))
  call GET_GROUP_ID(Ncid_Eps_Sg, 'status', Group_Id1)
  call GET_GROUP_ID(Group_Id1, 'satellite', Group_Id2)
  call READ_NETCDF(Group_Id2, [1], [1], [1],"earth_sun_distance_ratio", Sds_Data_1d)
  Earth_Sun_Ratio = Sds_Data_1d(1)
  deallocate(Sds_Data_1d)


  !--- read solar ch averaged solar irradiance
  call GET_GROUP_ID(Ncid_Eps_Sg, 'data', Group_Id1)
  call GET_GROUP_ID(Group_Id1, 'calibration_data', Group_Id2)
  allocate(Sds_Data_1d_real4(Num_Sol_Ch))
  call READ_NETCDF(Group_Id2,[1], [1], [Num_Sol_Ch],"Band_averaged_solar_irradiance" &
      , Sds_Data_1d_real4)
    
  Ch_Aver_Sol_Irrad  = Sds_Data_1d_real4
  deallocate(Sds_Data_1d_real4)

  !--- allocate a temporary array to hold 2d data
  allocate(Sds_Data_2d(Sds_Count_2d(1),Sds_Count_2d(2)))

  ! --- get group id
  !call GET_GROUP_ID(Ncid_Eps_Sg, 'data', Group_Id1)
  call GET_GROUP_ID(Group_Id1, 'measurement_data', Group_Id2)


  ! --- read channels
  do Chan_Idx = 1, Sensor%Num_Chan_Sensor

     CLAVRx_Chan_Idx = Sensor%CLAVRx_Chan_Map(Chan_Idx)

     if (Sensor%Chan_On_Flag_Default(Clavrx_Chan_Idx) == sym%No) cycle

     call READ_AND_UNSCALE_NETCDF_2D(Group_Id2, Sds_Start_2d, Sds_Stride_2d, Sds_Count_2d, &
                                "vii_"//trim(Ch_Name(Chan_Idx)), Sds_Data_2d)
      ! call file_to_clavrx (Sds_Data_2d, Sds_Data_2d_clavrx_grid)

     !--- solar channels
     if (Ch(CLAVRx_Chan_Idx)%Obs_Type == SOLAR_OBS_TYPE) then

       !--- convert radiance to reflectance
       call CONVERT_RAD_TO_REF(Earth_Sun_Ratio, Ch_Aver_Sol_Irrad(Chan_Idx), Sds_Data_2d)
        !print *,Sds_Data_2d(1:5,1:5)


        x_count  = Sds_Count_2d(1)
        y_count = Sds_Count_2d(2)

       !--- store reflectance
       Ch(CLAVRx_Chan_Idx)%Ref_Toa(1:x_count,1:y_count) = &
                         Solar_Corr_Coeff(Chan_Idx) * Sds_Data_2d

      end if

      !--- thermal and mixed channels
     if (Ch(CLAVRx_Chan_Idx)%Obs_Type == THERMAL_OBS_TYPE .or. &
         Ch(CLAVRx_Chan_Idx)%Obs_Type == MIXED_OBS_TYPE) then
       !--- store radiance


       Ch(CLAVRx_Chan_Idx)%Rad_Toa(1:x_count,1:y_count) = Sds_Data_2d
       !--- convert radiance to noaa units
       call CONVERT_RADIANCE (Ch(CLAVRx_Chan_Idx)%Rad_Toa, Planck_Nu(CLAVRx_Chan_Idx),MISSING_VALUE_REAL4)
       !--- convert radiance to brightness temperature
       call COMPUTE_BT_ARRAY(Ch(CLAVRx_Chan_Idx)%Bt_Toa,ch(CLAVRx_Chan_Idx)%Rad_Toa,CLAVRx_Chan_Idx,MISSING_VALUE_REAL4)

     endif

  enddo ! read channels


  !--- adjust start and count for 1d
  Num_Scan = Image%Number_Of_Lines/24

  Sds_Start_1d = 1
  Sds_Count_1d = Num_Scan
  Sds_Stride_1d(1) = 1


  !--- allocate a temporary array to hold 2d data
  allocate(Sds_Data_1d(Sds_Count_1d(1)))
  allocate(Time_Msec_Day(Sds_Count_1d(1)))

  !--- read scanline time

  call READ_NETCDF(Group_Id2, [Sds_Start_1d], [Sds_Stride_1d], [Sds_Count_1d], "scan_start_time_utc", Sds_Data_1d)

  !--- make sure not exceed max msec per day
  where(Sds_Data_1d > MSEC_PER_DAY)
    Sds_Data_1d = Sds_Data_1d - MSEC_PER_DAY
  end where

  Sds_Data_1d = (/(Image%Start_Time+(Image%End_Time - Image%Start_Time)/Sds_Count_1d(1)*(k-1), k=1,Sds_Count_1d(1))/)

  Time_Msec_Day = real(mod(int(Sds_Data_1d, kind=i8),  &
       Microsec_Per_Day/1000_i8), kind=f8)
  Image%Scan_Time_Ms(1:Sds_Count_2d(2)) = (/(Time_Msec_Day((k - 1) /24 + 1), k=Sds_Ny_Start, Sds_Ny_End)/)

  !--- special case last segment
  if(Sds_Nx_End == Image%Number_Of_Lines) then
    Tmp_Start = Sds_Count_1d(1) * 24
    Image%Scan_Time_Ms(Tmp_Start:Sds_Count_2d(1)) = Time_Msec_Day(Sds_Count_1d(1))
  endif


  !--- deallocate temporary memory
  deallocate(Sds_Data_1d)
  deallocate(Sds_Data_2d)
  deallocate(Time_Msec_Day)


  !--- calculate other angles

  Geo%Relaz = RELATIVE_AZIMUTH (Geo%Solaz, Geo%Sataz)
  Geo%Glintzen = GLINT_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz )
  Geo%Scatangle = SCATTERING_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz)

  !--- There are no gaps in vgac
  Gap_Pixel_Mask = sym%NO

  !--- record number of scans read for this segment
  Image%Number_Of_Lines_Read_This_Segment = Sds_Count_2d(2)

  !--- populate scan line number
  do Line_Idx = 1, y_count
    Image%Scan_Number(Line_Idx) = sds_Ny_Start + Line_Idx - 1

  enddo

end subroutine READ_EPS_SG_DATA


!-------------------------------------------------------------------------------
! converting for solar channels radiance (Wm^-2 sr^-1 µm^-) to reflectance (%)
!-------------------------------------------------------------------------------
subroutine CONVERT_RAD_TO_REF(Earth_Sun_Rat, Aver_Sol_Irrad, Buffer_2d)
  use univ_kind_defs_mod, only: f8

  real(f8), intent(in) :: Earth_Sun_Rat
  real, intent(in) :: Aver_Sol_Irrad
  real, dimension(:,:), intent(inout) :: Buffer_2d

  real, parameter :: Pi = 3.14159265359
  integer :: size_arr(2)

  !--- convert only valid data
  size_arr = shape(Buffer_2d)

  where(Buffer_2d .ne. MISSING_VALUE_REAL4)
    Buffer_2d = (Pi * Buffer_2d) / (cos(Geo%Solzen(1:size_arr(1),1:size_arr(2)) * Pi / 180.0))  &
                * (1. / Aver_Sol_Irrad) * Earth_Sun_Rat**2 * 100.0
  end where


  !--- set to missing negative
  where(Buffer_2d .lt. 0.0)
    Buffer_2d = MISSING_VALUE_REAL4
  end where

  !--- set reflectance > 120. to 120.
  where(Buffer_2d .gt. 120.0)
    Buffer_2d = 120.0
  end where

end subroutine CONVERT_RAD_TO_REF


subroutine file_to_clavrx (arr_old, arr_new)
    real,dimension(:,:) :: arr_old
    real,dimension(:,:) :: arr_new
    integer :: xx, yy
    integer :: old_sh(2), new_sh(2)

    old_sh = shape(arr_old)

    do xx =1 , old_sh(1)
      do yy = 1, old_sh(2)
        arr_new (yy,xx) = arr_old (xx,yy)

      end do
   end do


end subroutine


!-------------------------------------------------------------------------------
!  This routine computes navigation
!   AW 2020 Mai 12
!-------------------------------------------------------------------------------
subroutine EPS_NAVIGATION_SET_NAV ( this , nx_start, ny_start, nx, ny )
  use univ_kind_defs_mod, only: f8, f4

  class ( eps_navigation ) :: this

  integer , intent(in) ::  nx_start, ny_start, nx, ny

  integer :: Group_Id1,Group_Id2
  integer :: n_x
  integer :: n_y

  integer, dimension(2) :: start2d, stride2d, count2d
  integer :: kk_alt, kk_act
  integer :: i_zone_alt, i_zone_act
  integer, parameter :: num_pix_alt = 24
!  integer, parameter :: num_scans = 35 !175
  integer :: num_scans
  !integer, parameter :: num_lines = num_pix_alt * num_scans
  integer, parameter :: num_elem = 3144

  integer, parameter :: zone_size_act = 8
  integer, parameter :: zone_size_alt = 8

  real(f8), parameter :: ECCENTRICITY =1/298.257223563
  real(f8), parameter :: SEMI_MAJOR = 6378137.
  real(f8), parameter :: R_EARTH = 6371.
  real(f8) :: SEMI_MINOR = SEMI_MAJOR*(1-ECCENTRICITY)
  character(len=120) :: var_name



  integer :: i_scan
  real :: kkk

  integer :: i_rel_alt, i_rel_act
  real(f4) :: val_for_interp

  real, allocatable, dimension(:,:) :: longitude_anchor
  real, allocatable, dimension(:,:) :: latitude_anchor
  real, allocatable, dimension(:,:) :: observation_azimuth_anchor
  real, allocatable, dimension(:,:) :: observation_zenith_anchor
  real, allocatable, dimension(:,:) :: solar_azimuth_anchor
  real, allocatable, dimension(:,:) :: solar_zenith_anchor
  real, allocatable, dimension(:,:) :: cart_x, cart_y, cart_z, cart_n

  real:: this_x, this_y,  this_z
  real :: phi, p, e_strich
  integer :: xx, yy

  integer, dimension(2):: Dim_2d


  if ( this % is_calculated ) return

  Dim_2d = -1
  call GET_GROUP_ID(Ncid_Eps_Sg, 'data', Group_Id1)
  call GET_GROUP_ID(Group_Id1, 'measurement_data', Group_Id2)

   call READ_NETCDF_DIMENSION_2D(Group_Id2,'delta_lat_N_dem', Dim_2d)

   num_scans = Dim_2d(2)/24


   ! status = nf90_inq_varid(group_id1, trim('num_scans'), nc_var_id)
   ! call read_netcdf_attribute_real(gropu_id1, nc_var_id, 'CFAC', attr)

  n_x = 394 !140 !num_scans * 4
  n_y = num_scans * 4 !394

  !n_x = num_scans * 4
  !n_y = 394

  start2d = [1,1]
  stride2d = [1,1]
  count2d = [n_x,n_y]

  allocate (longitude_anchor(n_x,n_y))
  allocate (latitude_anchor(n_x,n_y))
  allocate (observation_azimuth_anchor(n_x,n_y))
  allocate (observation_zenith_anchor(n_x,n_y))
  allocate (solar_azimuth_anchor(n_x,n_y))
  allocate (solar_zenith_anchor(n_x,n_y))
  allocate (cart_x(n_x,n_y))
  allocate (cart_y(n_x,n_y))
  allocate (cart_z(n_x,n_y))
  allocate (cart_n(n_x,n_y))


  allocate (this % longitude(nx,ny))
  allocate (this % latitude(nx,ny))
  allocate (this % observation_azimuth(nx,ny))
  allocate (this % observation_zenith(nx,ny))
  allocate (this % solar_azimuth(nx,ny))
  allocate (this % solar_zenith(nx,ny))


  ! - set to missing
  longitude_anchor = MISSING_VALUE_REAL4
  latitude_anchor = MISSING_VALUE_REAL4
  observation_azimuth_anchor = MISSING_VALUE_REAL4
  observation_zenith_anchor = MISSING_VALUE_REAL4
  solar_azimuth_anchor = MISSING_VALUE_REAL4
  solar_zenith_anchor = MISSING_VALUE_REAL4
  this % longitude = MISSING_VALUE_REAL4
  this % latitude = MISSING_VALUE_REAL4
  this % observation_azimuth = MISSING_VALUE_REAL4
  this % observation_zenith = MISSING_VALUE_REAL4
  this % solar_azimuth = MISSING_VALUE_REAL4
  this % solar_zenith = MISSING_VALUE_REAL4

  ! --- get group id
!  call GET_GROUP_ID(Ncid_Eps_Sg, 'data', Group_Id1)
  call GET_GROUP_ID(Group_Id1, 'measurement_data', Group_Id2)

  ! --- navigation

  call READ_AND_UNSCALE_NETCDF_2D(Group_Id2, Start2d, Stride2d, Count2d, &
                                "longitude", longitude_anchor)



 ! - MetImage is longitude form 0 to 360.
 ! - fix 2020/06/08 AW
 longitude_anchor = modulo((longitude_anchor+180.),360.)-180.

  call READ_AND_UNSCALE_NETCDF_2D(Group_Id2, Start2d, Stride2d, Count2d, &
                                "latitude",  latitude_anchor)

  call READ_AND_UNSCALE_NETCDF_2D(Group_Id2, Start2d, Stride2d, Count2d, &
                                "observation_azimuth",  observation_azimuth_anchor)

  call READ_AND_UNSCALE_NETCDF_2D(Group_Id2, Start2d, Stride2d, Count2d, &
                                "observation_zenith",  observation_zenith_anchor)

  call READ_AND_UNSCALE_NETCDF_2D(Group_Id2, Start2d, Stride2d, Count2d, &
                                "solar_azimuth",  solar_azimuth_anchor)

  call READ_AND_UNSCALE_NETCDF_2D(Group_Id2, Start2d, Stride2d, Count2d, &
                                "solar_zenith",  solar_zenith_anchor)



  ! - compute cartesian coordinates
  cart_n = SEMI_MAJOR/SQRT(1-ECCENTRICITY * (2. - ECCENTRICITY ) * sin ( latitude_anchor*DTOR) ** 2.)
  cart_x = cart_N * cos ( latitude_anchor*DTOR) * cos ( longitude_anchor*DTOR )
  cart_y = cart_N * cos ( latitude_anchor*DTOR) * sin ( longitude_anchor*DTOR )
  cart_z = (( 1- ECCENTRICITY) ** 2.) * cart_N * sin ( latitude_anchor*DTOR )

  where ( latitude_anchor .eq. -999.  .or. longitude_anchor .eq. -999. )
      cart_x = -999.
      cart_y = -999.
      cart_z = -999.
  end where




   do kk_alt = ny_start, ny_start + ny - 1  ! 200 (seg-size)

       i_zone_alt = (kk_alt -1) /  zone_size_alt + 1
       i_rel_alt = mod( kk_alt-1,zone_size_alt)
       i_scan = (kk_alt-1) / num_pix_alt

       do kk_act = nx_start ,nx_start + nx - 1 ! 3144

          i_zone_act = (kk_act-1) / zone_size_act + 1
          i_rel_act = mod( kk_act-1, zone_size_act)

          ! print*,kk_alt,kk_act,i_zone_alt,i_zone_act,i_rel_alt,i_rel_act

         !- transformation to Cartesian is needed to avoid dateline trouble

          this_x = eps_bilin_interp ( &
            cart_x &
              ,i_zone_alt  &
              , i_scan &
              , i_zone_act &
              , i_rel_act &
              , i_rel_alt &
              ,zone_size_act &
              ,zone_size_alt )

         ! stop


           this_y = eps_bilin_interp ( &
            cart_y &
              ,i_zone_alt  &
              , i_scan &
              , i_zone_act &
              , i_rel_act &
              , i_rel_alt &
              ,zone_size_act &
              ,zone_size_alt )


           this_z = eps_bilin_interp ( &
            cart_z &
              ,i_zone_alt  &
              , i_scan &
              , i_zone_act &
              , i_rel_act &
              , i_rel_alt &
              ,zone_size_act &
              ,zone_size_alt )

          p = SQRT(( this_x ** 2.) + (this_y ** 2.))
          e_strich = SQRT ( (SEMI_MAJOR ** 2. - SEMI_MINOR ** 2.) / SEMI_MINOR ** 2.)

          phi = atan(this_z * SEMI_major/ (p * SEMI_minor))


          xx = kk_act - nx_start + 1
          yy = kk_alt - ny_start + 1



          this % longitude(xx,yy) = (2 * ATAN(this_y/(this_x + p)))/DTOR

          this % latitude(xx,yy) =  ((ATAN((this_z  + e_strich ** 2. * SEMI_MINOR * sin(phi *DTOR ) ** 3.) &
                                                    /(p- ECCENTRICITY ** 2. * ACOS (phi *DTOR) **3.)))) / DTOR


       !  print*,xx,yy,this % longitude(xx,yy),this % longitude(xx,yy)
       !  if (kk_alt .gt. 28 .and. kk_act .gt. 28) stop

          if (this_x .eq. -999. .or. this_y .eq. -999. .or. this_z .eq. -999.) then
              this % longitude(xx,yy) = -999.
              this % latitude(xx,yy) = -999.
             ! print*,this_x,this_y, this_z
             ! print*,xx,yy
             ! print*,kk_act,kk_alt
             ! stop
          end if


         ! print*, xx,yy,this % longitude(xx,yy) , this % latitude(xx,yy)

          ! simpler version
          !this % longitude(xx,yy) =  ATAN2(this_y,this_x) / DTOR
          !this % latitude(xx,yy) =  ASIN ( this_z / SEMI_MAJOR)  / DTOR

         !print*, this % longitude(xx,yy) , this % latitude(xx,yy)

         !stop

           this % observation_azimuth(xx,yy) = eps_bilin_interp ( &
             observation_azimuth_anchor &
              ,i_zone_alt  &
              , i_scan &
              , i_zone_act &
              , i_rel_act &
              , i_rel_alt &
              ,zone_size_act &
              ,zone_size_alt )

           this % observation_zenith(xx,yy) = eps_bilin_interp ( &
             observation_zenith_anchor &
              ,i_zone_alt  &
              , i_scan &
              , i_zone_act &
              , i_rel_act &
              , i_rel_alt &
          ,zone_size_act &
          ,zone_size_alt )

       this % solar_azimuth(xx,yy) = eps_bilin_interp ( &
        solar_azimuth_anchor &
          ,i_zone_alt  &
          , i_scan &
          , i_zone_act &
          , i_rel_act &
          , i_rel_alt &
          ,zone_size_act &
          ,zone_size_alt )

         this % solar_zenith(xx,yy) = eps_bilin_interp ( &
       solar_zenith_anchor &
          ,i_zone_alt  &
          , i_scan &
              , i_zone_act &
              , i_rel_act &
              , i_rel_alt &
              ,zone_size_act &
              ,zone_size_alt )

       end do
       !stop
  end do

  !print*,minval(this % longitude), maxval(this % longitude)
  !print*,minval(this % latitude), maxval(this % latitude)





  this % is_calculated = .TRUE.

  deallocate (longitude_anchor)
  deallocate (latitude_anchor)
  deallocate (observation_azimuth_anchor)
  deallocate (observation_zenith_anchor)
  deallocate (solar_azimuth_anchor)
  deallocate (solar_zenith_anchor)
  deallocate (cart_x)
  deallocate (cart_y)
  deallocate (cart_z)
  deallocate (cart_n)


  contains

  function eps_bilin_interp ( val_array &
  , i_zone_alt  &
  , i_scan &
  , i_zone_act &
  , i_rel_act &
  , i_rel_alt &
  ,zone_size_act &
  ,zone_size_alt  )

   implicit none
  real, intent(in) :: val_array(:,:)
  integer, intent(in) :: i_zone_alt, i_scan, i_zone_act
  integer, intent(in) :: i_rel_act, i_rel_alt
  integer, intent(in) ::zone_size_act,zone_size_alt


  real :: eps_bilin_interp
  real :: r,s
  integer, dimension(2) :: box_a, box_b, box_c, box_d
  real,dimension(2,2) :: tbl
  integer :: size_arr (2)



   box_a = [i_zone_act,i_zone_alt + i_scan]
   box_b = [i_zone_act+1,i_zone_alt + i_scan]
   box_c = [i_zone_act+1,i_zone_alt + i_scan+1]
   box_d = [i_zone_act,i_zone_alt + i_scan+1]

 ! print*,box_a
 ! print*,box_b
 ! print*,box_c
 ! print*,box_d
   size_arr = shape(val_array)


   box_a(1) = minval([box_a(1),size_arr(1)])
   box_a(2) = minval([box_a(2),size_arr(2)])
   box_b(1) = minval([box_b(1),size_arr(1)])
   box_b(2) = minval([box_b(2),size_arr(2)])
   box_c(1) = minval([box_c(1),size_arr(1)])
   box_c(2) = minval([box_c(2),size_arr(2)])
   box_d(1) = minval([box_d(1),size_arr(1)])
   box_d(2) = minval([box_d(2),size_arr(2)])

   box_a(1) = maxval([box_a(1),1])
   box_a(2) = maxval([box_a(2),1])
   box_b(1) = maxval([box_b(1),1])
   box_b(2) = maxval([box_b(2),1])
   box_c(1) = maxval([box_c(1),1])
   box_c(2) = maxval([box_c(2),1])
   box_d(1) = maxval([box_d(1),1])
   box_d(2) = maxval([box_d(2),1])

   tbl(1,1) = val_array(box_a(1),box_a(2))
   tbl(2,1) = val_array(box_b(1),box_b(2))
   tbl(2,2) = val_array(box_c(1),box_c(2))
   tbl(1,2) = val_array(box_d(1),box_d(2))

   r  = i_rel_act/(1.*zone_size_act)
   s = i_rel_alt/(1.*zone_size_alt)


  ! print*,r,s
  ! print*,tbl
   eps_bilin_interp =  ( 1.0 - r ) * ( 1.0 - s ) * tbl( 1 , 1 ) &
         & + r * ( 1.0 - s )          * tbl( 2 , 1 ) &
         & + (1.0 - r ) * s           * tbl( 1 , 2 ) &
         & + r * s                    * tbl( 2 , 2 )

   if ( count ( tbl .eq. -999) .gt. 0) then
      print*,tbl
      print*,box_a
      print*,box_b
      print*,box_c
      print*,box_d
      print*,size_arr
      eps_bilin_interp =   -999.
      stop
   end if

   end function eps_bilin_interp

end subroutine EPS_NAVIGATION_SET_NAV


subroutine eps_navigation_dealloc (this)
  class(eps_navigation) :: this


  deallocate (this % longitude)
  deallocate (this % latitude)
  deallocate (this % observation_azimuth)
  deallocate (this % observation_zenith)
  deallocate (this % solar_azimuth)
  deallocate (this % solar_zenith)
  this % is_calculated = .FALSE.
end

!======================================================================

end module CX_EPS_SG_MOD
