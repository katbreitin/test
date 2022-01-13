!$Id:$
!----------------------------------------------------------------------
!CLAVR-x Next Generation ISCCP (ISCCPNG)
!
! Questions 
!   Do we keep logic of segments? aka  Will we fill memory if we don't?
!
!----------------------------------------------------------------------
module CX_ISCCPNG_MOD

!--- use statements
use PIXEL_COMMON_MOD, only: Ch, Geo, Image, Nav, Sensor,  L1g, &
                            Temp_Pix_Array_1, Gap_Pixel_Mask, WMO_Id_ISCCPNG, Pixel_Time

use NETCDF,only: nf90_inquire_dimension, nf90_inq_dimid, &
    nf90_noerr,nf90_inq_varid,nf90_get_var

Use CX_NETCDF4_MOD,only: OPEN_NETCDF, &
                         CLOSE_NETCDF, &
                         READ_NETCDF, &
                         READ_AND_UNSCALE_NETCDF_1D, &
                         READ_AND_UNSCALE_NETCDF_2D, &
                         READ_AND_UNSCALE_NETCDF_3D, &
                         READ_AND_UNSCALE_NETCDF_4D

use VIEWING_GEOMETRY_MOD, only: GLINT_ANGLE, SCATTERING_ANGLE, &
                                RELATIVE_AZIMUTH
use PLANCK_MOD
use CONSTANTS_MOD, only: MSEC_PER_DAY, Sym, MISSING_VALUE_REAL4, &
                         MISSING_VALUE_REAL8, &
                         SOLAR_OBS_TYPE, THERMAL_OBS_TYPE, &
                         MIXED_OBS_TYPE, LUNAR_OBS_TYPE

use CALIBRATION_CONSTANTS_MOD
use CX_REAL_BOOLEAN_MOD
use CLASS_TIME_DATE

implicit none

private

public:: READ_NUMBER_OF_SCANS_ISCCPNG
public:: READ_ISCCPNG_DATE_TIME
public:: READ_ISCCPNG_DATA

private:: READ_WMO_L1G
private:: READ_NAV_L1G
private:: READ_SINGLE_L1G
private:: READ_LAYERED_L1G
private:: CONVERT_WMO

integer,parameter, private:: Num_Layers_L1g = 3
integer,parameter, private:: Num_Time_L1g = 1

contains

!----------------------------------------------------------------------
! open and read dimensions
!----------------------------------------------------------------------
subroutine READ_NUMBER_OF_SCANS_ISCCPNG(Nlat,Nlon,Ierror)
   integer, intent(out):: Nlat,Nlon,Ierror
   integer :: status, dim_id, Ncid
   character(len=100):: dim_name_dummy

   Ierror = 0

   call OPEN_NETCDF(Image%Level1b_Full_Name,Ncid)

   if (Ncid < 0) then
     Ierror = 1
     return
   endif

   status = nf90_inq_dimid(Ncid, 'latitude', dim_id)
   if (status /= nf90_noerr) then
       print *, 'ERROR: nf90_inq_dimid(Ncid, latitude, dim_id)'
       Ierror = 1
       return
   endif
   status = nf90_inquire_dimension(Ncid,dim_id,dim_name_dummy,len=Nlat)
   if (status /= nf90_noerr) then
       print *, 'ERROR: nf90_inquire_dimension(...,len=Nlat)'
       Ierror = 1
       return
   endif

   status = nf90_inq_dimid(Ncid, 'longitude', dim_id)
   if (status /= nf90_noerr) then
       print *, 'ERROR: nf90_inq_dimid(Ncid, longitude, dim_id)'
       Ierror = 1
       return
   endif
   status = nf90_inquire_dimension(Ncid,dim_id,dim_name_dummy,len=Nlon)
   if (status /= nf90_noerr) then
       print *, 'ERROR: nf90_inquire_dimension(...,len=Nlon)'
       Ierror = 1
       return
   endif
   
   call CLOSE_NETCDF(Ncid)

end subroutine READ_NUMBER_OF_SCANS_ISCCPNG

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine READ_ISCCPNG_DATE_TIME(Start_Year,Start_Doy,Start_Time,&
                                  End_Year,End_Doy,End_Time)
   integer, intent(out) :: Start_Year  !year
   integer, intent(out) :: Start_Doy   !day of year
   integer, intent(out) :: Start_Time  !millisec
   integer, intent(out) :: End_Year    !year
   integer, intent(out) :: End_Doy     !day of year
   integer, intent(out) :: End_Time    !millisec
   integer:: Start_Hour, Start_Minute
   integer:: Ncid

   integer, dimension(1) :: var_start
   integer, dimension(1) :: var_stride
   integer, dimension(1) :: var_dim
   character(len=100) :: var_name
   real, dimension(1) :: var_output
   integer :: Start_Time_L1g, End_Time_L1g
   character(len=4):: Year_String
   character(len=3):: Doy_String
   character(len=2):: Hour_String
   character(len=2):: Minute_String
   integer:: Hour_Temp, Minute_Temp
   type(date_type) :: l1g_start_time
   type(date_type) :: l1g_end_time


   call OPEN_NETCDF(Image%Level1b_Full_Name,Ncid)

   !--- read dates and time from level1b
   var_start = 1
   var_stride = 1
   var_dim = 1

   !--- start time
   call READ_NETCDF(Ncid, var_start, var_stride, var_dim, "start_time", var_output)
   start_time_l1g = var_output(1)
   call l1g_start_time % SET_DATE(year=1970,month=1,hour=0,minute=0,second=0)
   call l1g_start_time % ADD_TIME(second = Start_Time_L1g) 
   Start_Year = l1g_start_time % year  
   Start_Doy = l1g_start_time % dayofyear
   Start_Time = l1g_start_time % msec_of_day

   !--- end time
   call READ_NETCDF(Ncid, var_start, var_stride, var_dim, "end_time", var_output)
   end_time_l1g = var_output(1)
   call l1g_end_time % SET_DATE(year=1970,month=1,hour=0,minute=0,second=0)
   call l1g_end_time % ADD_TIME(second = End_Time_L1g) 
   end_time_l1g = var_output(1)
   End_Year = l1g_end_time % year  
   End_Doy = l1g_end_time % dayofyear
   End_Time = l1g_end_time % msec_of_day

   call CLOSE_NETCDF(Ncid)

end subroutine READ_ISCCPNG_DATE_TIME
!--------------------------------------------------------------------------------------------------
! public read routine
!--------------------------------------------------------------------------------------------------
subroutine READ_ISCCPNG_DATA(Segment_Number, Error_Status)
  integer, intent(in):: Segment_Number
  integer, intent(out):: Error_Status
  integer, dimension(2) :: Sds_Start
  integer, dimension(2) :: Sds_Stride
  integer, dimension(2) :: Sds_Count
  real(kind=8), dimension(:), allocatable:: Sds_Data_1d
  real, dimension(:,:), allocatable:: Sds_Data_2d
  integer:: Nx_Start, Nx_End, Ny_Start, Ny_End
  integer:: Chan_Idx, CLAVRx_Chan_Idx, Line_Idx, Elem_Idx
  real(kind=8) :: proj_time0
  integer:: varid, status, Ncid
  integer(kind=2), dimension(:,:,:), allocatable:: WMO_L1g
  integer(kind=2), dimension(:,:), allocatable:: WMO_Temp
  integer(kind=1), dimension(:,:,:), allocatable:: Sample_Mode_L1g
  character(len=1020):: File_Prefix, File_Suffix
  integer:: ipos, ilen, ilen_wmo_id
  integer:: Lat_Idx, Layer_Idx

  !--- determine location of segement data in the file
  Nx_Start = 1
  Nx_End = Nx_Start + Image%Number_Of_Elements - 1
  Ny_Start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
  Ny_End = min(Image%Number_Of_Lines, Ny_Start + Image%Number_of_Lines_Per_Segment - 1)
  Image%Number_Of_Lines_Read_This_Segment = Ny_End - Ny_Start + 1

  Sds_Start = [1,Ny_Start]
  Sds_Stride = [1,1]
  Sds_Count = [Nx_End - Nx_Start + 1, Ny_End - Ny_Start + 1] 

  !--- read L1g data and make a field of WMO and Layer for future reads
  allocate(WMO_L1g(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment,Num_Layers_L1g))
  call READ_WMO_L1G(WMO_Id_ISCCPNG, Segment_Number, WMO_L1g) 

  !--- make 2d fields of WMO and Layer_Idx
  allocate(WMO_Temp(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
  do  Layer_Idx = 1,Num_Layers_L1g
     WMO_Temp = WMO_L1g(:,:,Layer_Idx)
     where (WMO_Temp == WMO_Id_ISCCPNG)
             L1g%WMO_Id = WMO_Temp
             L1g%Layer_Idx = Layer_Idx
     endwhere
  enddo

  !-- lat and lon from wmo_id file
  call READ_NAV_L1G(Segment_Number) 

  !---- manipulate file name into its two roots so we can insert sds_name and grab right file
  ipos = index(Image%Level1b_Full_Name,"__")
  ilen = len(trim(Image%Level1b_Full_Name))
  ilen_wmo_id = len('wmo_id')
  File_Prefix = Image%Level1b_Full_Name(1:ipos+1)
  File_Suffix = Image%Level1b_Full_Name(ipos+ilen_wmo_id+2:ilen)

  call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'sample_mode', Segment_Number, WMO_L1g, Temp_Pix_Array_1) 
  L1g%Sample_Mode = int(Temp_Pix_Array_1,kind=1)

  !--- generic call for constructng one 2d field for this WMO and all layers
  call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'satellite_zenith_angle', Segment_Number, WMO_L1g, Geo%Satzen) 
  call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'satellite_azimuth_angle', Segment_Number, WMO_L1g, Geo%Sataz) 

  !--- solar angles are a single layer
  call READ_SINGLE_L1G(File_Prefix, File_Suffix, 'solar_zenith_angle', Segment_Number, WMO_L1g, Geo%Solzen) 
  call READ_SINGLE_L1G(File_Prefix, File_Suffix, 'solar_azimuth_angle', Segment_Number, WMO_L1g, Geo%Solaz) 

  !--- based Space Mask on Satzen
  Geo%Space_Mask = sym%YES
  where(Geo%Satzen /= MISSING_VALUE_REAL4)
     Geo%Space_Mask = sym%NO
  endwhere

  !--- convert azimuths from L1g to CLAVR-x convention
  ! L1g is clockwise from West
  ! CLAVR-x is clockwise from North
  ! l1g = (clavrx + 90)%360.
  where(Geo%Sataz /= MISSING_VALUE_REAL4)
     Geo%Sataz = Geo%Sataz + 90.0
  endwhere
  where(Geo%Sataz > 180.0) 
     Geo%Sataz = Geo%Sataz - 360.0
  endwhere
  where(Geo%Solaz /= MISSING_VALUE_REAL4)
     Geo%Solaz = Geo%Solaz + 90.0
  endwhere
  where(Geo%Solaz > 180.0) 
     Geo%Solaz = Geo%Solaz - 360.0
  endwhere

  !--- make Relaz and Glint and Scattering angles here
  Geo%Relaz = RELATIVE_AZIMUTH (Geo%Solaz, Geo%Sataz)
  Geo%Glintzen = GLINT_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz )
  Geo%Scatangle = SCATTERING_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz)

  !--- read time and convert to clavrx standard
  call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'pixel_time', Segment_Number, WMO_L1g, Pixel_Time) 
  where(Pixel_Time /= MISSING_VALUE_REAL4)
       Pixel_Time = 1000.0*Pixel_Time
  endwhere

  !--- make scanline time
  do Lat_Idx = 1, Image%Number_Of_Lines_Per_Segment
     Image%Scan_Time_Ms(Lat_Idx) = maxval(Pixel_Time(:,Lat_Idx))
     Image%Scan_Number(Lat_Idx) = (Segment_Number - 1)*Image%Number_of_Lines_Per_Segment + Lat_Idx
  enddo

  !--- read radiometric data
  if (Sensor%Chan_On_Flag_Default(3) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'refl_00_47um', Segment_Number, WMO_L1g, Ch(3)%Ref_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(4) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'refl_00_51um', Segment_Number, WMO_L1g, Ch(4)%Ref_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(1) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'refl_00_65um', Segment_Number, WMO_L1g, Ch(1)%Ref_Toa) 
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'refl_00_65um_std', Segment_Number, WMO_L1g, Ch(1)%Ref_Toa_Std_Sub) 
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'refl_00_65um_min', Segment_Number, WMO_L1g, Ch(1)%Ref_Toa_Min_Sub) 
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'refl_00_65um_max', Segment_Number, WMO_L1g, Ch(1)%Ref_Toa_Max_Sub) 
  endif
  if (Sensor%Chan_On_Flag_Default(2) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'refl_00_86um', Segment_Number, WMO_L1g, Ch(2)%Ref_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(26) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'refl_01_38um', Segment_Number, WMO_L1g, Ch(26)%Ref_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(6) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'refl_01_60um', Segment_Number, WMO_L1g, Ch(6)%Ref_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(7) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'refl_02_20um', Segment_Number, WMO_L1g, Ch(7)%Ref_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(20) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_03_80um', Segment_Number, WMO_L1g, Ch(20)%Bt_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(37) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_06_20um', Segment_Number, WMO_L1g, Ch(37)%Bt_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(27) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_06_70um', Segment_Number, WMO_L1g, Ch(27)%Bt_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(28) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_07_30um', Segment_Number, WMO_L1g, Ch(28)%Bt_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(29) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_08_50um', Segment_Number, WMO_L1g, Ch(29)%Bt_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(30) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_09_70um', Segment_Number, WMO_L1g, Ch(30)%Bt_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(38) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_10_40um', Segment_Number, WMO_L1g, Ch(38)%Bt_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(31) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_11_00um', Segment_Number, WMO_L1g, Ch(31)%Bt_Toa) 
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_11_00um_std', Segment_Number, WMO_L1g, Ch(31)%Bt_Toa_Std_Sub) 
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_11_00um_min', Segment_Number, WMO_L1g, Ch(31)%Bt_Toa_Min_Sub) 
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_11_00um_max', Segment_Number, WMO_L1g, Ch(31)%Bt_Toa_Max_Sub) 
    call COMPUTE_RAD_ARRAY( Ch(31)%Bt_Toa,Ch(31)%Rad_Toa,31,MISSING_VALUE_REAL4)
  endif
  if (Sensor%Chan_On_Flag_Default(32) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_12_00um', Segment_Number, WMO_L1g, Ch(32)%Bt_Toa) 
  endif
  if (Sensor%Chan_On_Flag_Default(33) == sym%Yes) then
    call READ_LAYERED_L1G(File_Prefix, File_Suffix, 'temp_13_30um', Segment_Number, WMO_L1g, Ch(33)%Bt_Toa) 
  endif

  deallocate(WMO_L1g, WMO_Temp)

  !--- close file on last segment if still open
  if (Image%Segment_Number == Image%Number_Of_Segments .and. Ncid > 0) then
    call CLOSE_NETCDF(Ncid)
  endif

end subroutine READ_ISCCPNG_DATA
!--------------------------------------------------------------------------------------------------
! Read L1g variables that are layered (ie. brightness temperatures)
!--------------------------------------------------------------------------------------------------
subroutine READ_LAYERED_L1G(File_Name_Prefix, File_Name_Suffix, Sds_Name, Segment_Number, WMO_L1g_3d, Sds_Data)
  character(len=*), intent(in):: File_Name_Prefix, File_Name_Suffix, Sds_Name
  integer, intent(in):: Segment_Number
  real, intent(out), dimension(:,:) :: Sds_Data
  integer(kind=2), dimension(:,:,:), intent(in):: WMO_L1g_3d
  integer(kind=2), dimension(:,:), allocatable:: WMO_L1g_2d
  real, dimension(:,:), allocatable:: Sds_Data_2d
  real, dimension(:,:,:,:), allocatable:: Sds_Data_4d
  character(len=1020):: File_Name
  integer:: Nx_Start, Nx_End, Ny_Start, Ny_End, Layer_Idx
  integer, dimension(4):: Sds_Start, Sds_Count, Sds_Stride
  integer:: Ncid

  File_Name = trim(File_Name_Prefix)// trim(Sds_Name)//trim(File_Name_Suffix)

  call OPEN_NETCDF(File_Name,Ncid)

  Nx_Start = 1
  Nx_End = Nx_Start + Image%Number_Of_Elements - 1
  Ny_Start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
  Ny_End = min(Image%Number_Of_Lines, Ny_Start + Image%Number_of_Lines_Per_Segment - 1)

  Sds_Start = [1,Ny_Start,1,1]
  Sds_Stride = [1,1,1,1]
  Sds_Count = [Nx_End - Nx_Start + 1, Ny_End - Ny_Start + 1,Num_Layers_L1g,Num_Time_L1g]

  allocate(Sds_Data_4d(Image%Number_Of_Elements, Image%Number_of_Lines_Per_Segment, Num_Layers_L1g,Num_Time_L1g))
  allocate(Sds_Data_2d(Image%Number_Of_Elements, Image%Number_of_Lines_Per_Segment))
  allocate(WMO_L1g_2d(Image%Number_Of_Elements, Image%Number_of_Lines_Per_Segment))

  call READ_AND_UNSCALE_NETCDF_4D(Ncid, Sds_Start, Sds_Stride, Sds_Count, Sds_Name, Sds_Data_4d)

  Sds_Data = MISSING_VALUE_REAL4
  do  Layer_Idx = 1,Num_Layers_L1g
    WMO_L1g_2d = WMO_L1g_3d(:,:,Layer_Idx)
    Sds_Data_2d = Sds_Data_4d(:,:,Layer_Idx,1)
    where (WMO_L1g_2d == WMO_Id_ISCCPNG)
            Sds_Data =  Sds_Data_2d
    endwhere
  enddo

  call CLOSE_NETCDF(Ncid)

  deallocate(Sds_Data_4d, Sds_Data_2d, WMO_L1g_2d)

!  if (Segment_Number == 2 .and. trim(Sds_Name) == "temp_11_00um") stop

end subroutine READ_LAYERED_L1G
!--------------------------------------------------------------------------------------------------
! Read L1g variables that are not layered (ie. solar angles)
!--------------------------------------------------------------------------------------------------
subroutine READ_SINGLE_L1G(File_Name_Prefix, File_Name_Suffix, Sds_Name, Segment_Number, WMO_L1g_3d, Sds_Data)
  character(len=*), intent(in):: File_Name_Prefix, File_Name_Suffix, Sds_Name
  integer, intent(in):: Segment_Number
  real, intent(out), dimension(:,:) :: Sds_Data
  integer(kind=2), dimension(:,:,:), intent(in):: WMO_L1g_3d
  integer(kind=2), dimension(:,:), allocatable:: WMO_L1g_2d
  real, dimension(:,:), allocatable:: Sds_Data_2d
  real, dimension(:,:,:,:), allocatable:: Sds_Data_4d
  character(len=1020):: File_Name
  integer:: Nx_Start, Nx_End, Ny_Start, Ny_End, Layer_Idx
  integer, dimension(4):: Sds_Start, Sds_Count, Sds_Stride
  integer:: Ncid

  File_Name = trim(File_Name_Prefix)// trim(Sds_Name)//trim(File_Name_Suffix)

  call OPEN_NETCDF(File_Name,Ncid)

  Nx_Start = 1
  Nx_End = Nx_Start + Image%Number_Of_Elements - 1
  Ny_Start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
  Ny_End = min(Image%Number_Of_Lines, Ny_Start + Image%Number_of_Lines_Per_Segment - 1)

  Sds_Start = [1,Ny_Start,1,1]
  Sds_Stride = [1,1,1,1]
  Sds_Count = [Nx_End - Nx_Start + 1, Ny_End - Ny_Start + 1,1,Num_Time_L1g]

  allocate(Sds_Data_4d(Image%Number_Of_Elements, Image%Number_of_Lines_Per_Segment, 1,Num_Time_L1g))
  allocate(Sds_Data_2d(Image%Number_Of_Elements, Image%Number_of_Lines_Per_Segment))
  allocate(WMO_L1g_2d(Image%Number_Of_Elements, Image%Number_of_Lines_Per_Segment))

  call READ_AND_UNSCALE_NETCDF_4D(Ncid, Sds_Start, Sds_Stride, Sds_Count, Sds_Name, Sds_Data_4d)

  Sds_Data = MISSING_VALUE_REAL4
  Sds_Data_2d = Sds_Data_4d(:,:,1,1)
  do  Layer_Idx = 1,Num_Layers_L1g
    WMO_L1g_2d = WMO_L1g_3d(:,:,Layer_Idx)
    where (WMO_L1g_2d == WMO_Id_ISCCPNG)
            Sds_Data =  Sds_Data_2d
    endwhere
  enddo

  call CLOSE_NETCDF(Ncid)

  deallocate(Sds_Data_4d, Sds_Data_2d, WMO_L1g_2d)

end subroutine READ_SINGLE_L1G
!--------------------------------------------------------------------------------------------------
! Construct from L1g a field of WMO and Layer where the L1g WMO equals the chosen WMO
!--------------------------------------------------------------------------------------------------
subroutine READ_WMO_L1G(WMO_Idx,Segment_Number, WMO_L1g)

  integer, intent(in) :: WMO_Idx, Segment_Number
  integer(kind=2), intent(out), dimension(:,:,:):: WMO_L1g

  integer:: Nx_Start, Nx_End, Ny_Start, Ny_End, Layer_Idx
  integer:: Ncid
  integer, dimension(4):: Sds_Start, Sds_Count, Sds_Stride
  integer(kind=2), dimension(:,:,:,:), allocatable:: Sds_Data_Temp

  call OPEN_NETCDF(Image%Level1b_Full_Name,Ncid)

  Nx_Start = 1
  Nx_End = Nx_Start + Image%Number_Of_Elements - 1
  Ny_Start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
  Ny_End = min(Image%Number_Of_Lines, Ny_Start + Image%Number_of_Lines_Per_Segment - 1)

  Sds_Start = [1,Ny_Start,1,1]
  Sds_Stride = [1,1,1,1]
  Sds_Count = [Nx_End - Nx_Start + 1, Ny_End - Ny_Start + 1,Num_Layers_L1g,Num_Time_L1g]

  allocate(Sds_Data_Temp(Image%Number_Of_Elements, Image%Number_of_Lines_Per_Segment, Num_Layers_L1g, Num_Time_L1g))

  Sds_Data_Temp = MISSING_VALUE_REAL4

  call READ_NETCDF(Ncid, Sds_Start, Sds_Stride, Sds_Count, "wmo_id", Sds_Data_Temp)

  WMO_L1g = Sds_Data_Temp(:,:,:,1)

  WMO_L1g = CONVERT_WMO(WMO_L1g)

  call CLOSE_NETCDF(Ncid)

  deallocate(Sds_Data_Temp)

end subroutine READ_WMO_L1G
!--------------------------------------------------------------------------------------------------
! Construct from L1g a field of WMO and Layer where the L1g WMO equals the chosen WMO
!--------------------------------------------------------------------------------------------------
subroutine READ_NAV_L1G(Segment_Number)

  integer, intent(in) :: Segment_Number

  integer:: Nx_Start, Nx_End, Ny_Start, Ny_End
  integer:: Ncid, Lon_Idx, Lat_Idx
  real, dimension(:), allocatable:: Lon_1d, Lat_1d

  Nx_Start = 1
  Nx_End = Nx_Start + Image%Number_Of_Elements - 1
  Ny_Start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
  Ny_End = min(Image%Number_Of_Lines, Ny_Start + Image%Number_of_Lines_Per_Segment - 1)

  call OPEN_NETCDF(Image%Level1b_Full_Name,Ncid)

  allocate(Lon_1d(Image%Number_of_Elements))
  call READ_AND_UNSCALE_NETCDF_1D(Ncid, [Nx_Start], [1],[Nx_End-Nx_Start+1] , "longitude", Lon_1d)

  allocate(Lat_1d(Image%Number_of_Lines_Per_Segment))
  call READ_AND_UNSCALE_NETCDF_1D(Ncid, [Ny_Start], [1],[Ny_End-Ny_Start+1] , "latitude", Lat_1d)

  do Lon_Idx = 1,Nx_End - Nx_Start + 1
    Nav%Lat_1b(Lon_Idx,1:Ny_End-Ny_Start+1) = Lat_1d
  enddo

  do Lat_Idx = 1,Ny_End-Ny_Start+1
    Nav%Lon_1b(1:Nx_End-Nx_Start+1,Lat_Idx) = Lon_1d
  enddo

  call CLOSE_NETCDF(Ncid)

  deallocate(Lon_1d, Lat_1d)

  Nav%Lon = Nav%Lon_1b
  Nav%Lat = Nav%Lat_1b

end subroutine READ_NAV_L1G

!----------------------------------------------------------------------
! Convert L1g WMO numbers into those used by CLAVR-x
! L1g: 152=GOES-16;664=GOES-17;167=Himawari-8;684=Meteosat-8;305=Meteosat-11"
! CLAVR-x comes from 
! https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/WMO306_vI2_CommonTable_en_v23.0.0.pdf
!----------------------------------------------------------------------
function CONVERT_WMO(WMO_Idx_In) result(WMO_Idx_Out)
  integer(kind=2), intent(in), dimension(:,:,:):: WMO_Idx_In
  integer(kind=2), dimension(size(WMO_Idx_in,1),size(WMO_Idx_In,2),size(WMO_Idx_In,3)):: WMO_Idx_Out

  WMO_Idx_Out = WMO_Idx_In

  where(WMO_Idx_In == 152) !GOES-16
    WMO_Idx_Out = 270
  endwhere
  where(WMO_Idx_In == 664) !GOES-17
    WMO_Idx_Out = 271
  endwhere
  where(WMO_Idx_In == 167) !HIM-8
    WMO_Idx_Out = 173
  endwhere
  where(WMO_Idx_In == 684) !MeteoSat-8
    WMO_Idx_Out = 55
  endwhere
  where(WMO_Idx_In == 305) !MeteoSat-11
    WMO_Idx_Out = 70
  endwhere

end function CONVERT_WMO

end module CX_ISCCPNG_MOD
