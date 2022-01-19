!$Id:$
!----------------------------------------------------------------------
!CLAVR-x VIIRS Global Area Coverage (VGAC)
!----------------------------------------------------------------------
module CX_VGAC_MOD

!--- use statements
use PIXEL_COMMON_MOD, only: Ch, Geo, Image, Nav, Sensor,  &
                            Temp_Pix_Array_1, Gap_Pixel_Mask

use NETCDF,only: nf90_inquire_dimension, nf90_inq_dimid, &
    nf90_noerr,nf90_inq_varid,nf90_get_var

Use CX_NETCDF4_MOD,only: OPEN_NETCDF, &
                         CLOSE_NETCDF, &
                         READ_NETCDF_DIMENSION_1D, &
                         READ_NETCDF_DIMENSION_4D, &
                         READ_NETCDF,READ_AND_UNSCALE_NETCDF_2D
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
public:: READ_NUMBER_OF_SCANS_VGAC
public:: READ_VGAC_DATE_TIME
public:: READ_VGAC_DATA

!--- module wide variables and parameters
integer, private, save:: Ncid_Vgac

character(len=3), private, parameter, dimension(16):: Obs_Name = &
         ["M01","M02","M03","M04","M05","M06","M07","M08","M09", &
          "M10","M11","M12","M13","M14","M15","M16"]
type(date_type), save :: Vgac_Start_Time, Vgac_End_Time

contains

!----------------------------------------------------------------------
! open and read dimensions
!----------------------------------------------------------------------
subroutine READ_NUMBER_OF_SCANS_VGAC(File_Name,Nscn,Npix,Ierror)
   character(len=*), intent(in):: File_Name
   integer, intent(out):: Nscn,Npix,Ierror
   integer :: status, dim_id
   character(len=100):: dim_name_dummy

   Ierror = 0

   call OPEN_NETCDF(File_Name,Ncid_Vgac)

   if (Ncid_Vgac < 0) then
     Ierror = 1
     return
   endif

   status = nf90_inq_dimid(Ncid_Vgac, 'nscn', dim_id)
   if (status /= nf90_noerr) then
       print *, 'ERROR: nf90_inq_dimid(Ncid_Vgac, nscn, dim_id)'
       Ierror = 1
       return
   endif
   status = nf90_inquire_dimension(Ncid_Vgac,dim_id,dim_name_dummy,len=Nscn)
   if (status /= nf90_noerr) then
       print *, 'ERROR: nf90_inquire_dimension(...,len=Nscn)'
       Ierror = 1
       return
   endif

   status = nf90_inq_dimid(Ncid_Vgac, 'npix', dim_id)
   if (status /= nf90_noerr) then
       print *, 'ERROR: nf90_inq_dimid(Ncid_Vgac, npix, dim_id)'
       Ierror = 1
       return
   endif
   status = nf90_inquire_dimension(Ncid_Vgac,dim_id,dim_name_dummy,len=Npix)
   if (status /= nf90_noerr) then
       print *, 'ERROR: nf90_inquire_dimension(...,len=Npix)'
       Ierror = 1
       return
   endif

end subroutine READ_NUMBER_OF_SCANS_VGAC

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine READ_VGAC_DATE_TIME(Start_Year,Start_Doy,Start_Time,&
                               End_Year,End_Doy,End_Time)
   integer, intent(out) :: Start_Year  !year
   integer, intent(out) :: Start_Doy   !day of year
   integer, intent(out) :: Start_Time  !millisec
   integer, intent(out) :: End_Year    !year
   integer, intent(out) :: End_Doy     !day of year
   integer, intent(out) :: End_Time    !millisec
   integer:: Start_Hour, Start_Minute

   integer, dimension(1) :: var_start
   integer, dimension(1) :: var_stride
   integer, dimension(1) :: var_dim
   character(len=100) :: var_name
   real, dimension(1) :: var_output
   real, dimension(:), allocatable:: Time
   real:: Proj_Time0, End_Time0
   character(len=4):: Year_String
   character(len=3):: Doy_String
   character(len=2):: Hour_String
   character(len=2):: Minute_String
   integer:: Hour_Temp, Minute_Temp

   !--- read dates and time from level1b
   var_start = 1
   var_stride = 1
   var_dim = 1
   call READ_NETCDF(Ncid_Vgac, var_start, var_stride, &
                    var_dim, "proj_time0", var_output)
   Proj_Time0 = var_output(1)

   allocate(Time(Image%Number_Of_Lines))
   Time = MISSING_VALUE_REAL4
   var_start = 1
   var_stride = 1
   var_dim = Image%Number_Of_Lines
   call READ_NETCDF(Ncid_Vgac, var_start, var_stride, var_dim, "time", Time)
   End_Time0 = maxval(Time)
   deallocate(Time)

   !--- convert to start time  
   call Vgac_Start_Time % SET_DATE(year=2010,month=1,hour=0,minute=0,second=0)
   call Vgac_Start_Time % ADD_DAYS(day = Proj_Time0) 

   !--- convert to end time  
   Vgac_End_Time = Vgac_Start_Time
   Hour_Temp = int(End_Time0)
   Minute_Temp = int(60.0*(End_Time0 - Hour_Temp))

   call Vgac_End_Time % ADD_TIME(hour = Hour_Temp,minute = Minute_Temp) 

   !--- store in output variables
   Start_Year = vgac_start_time % year  
   Start_Doy = vgac_start_time % dayofyear
   Start_Time = vgac_start_time % msec_of_day
   End_Year = vgac_end_time % year  
   End_Doy = vgac_end_time % dayofyear
   End_Time = vgac_end_time % msec_of_day

end subroutine READ_VGAC_DATE_TIME
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine READ_VGAC_DATA(Segment_Number, Error_Status)
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
  integer:: varid, status

  !--- determine location of segement data in the file
  Nx_Start = 1
  Nx_End = Nx_Start + Image%Number_Of_Elements - 1
  Ny_Start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
  Ny_End = min(Image%Number_Of_Lines, Ny_Start + Image%Number_of_Lines_Per_Segment - 1)

  Sds_Start = [1,Ny_Start]
  Sds_Stride = [1,1]
  Sds_Count = [Nx_End - Nx_Start + 1, Ny_End - Ny_Start + 1] 

  !--- allocate a temporary array to hold 1d and 2d data
  allocate(Sds_Data_1d(Sds_Count(2)))
  allocate(Sds_Data_2d(Sds_Count(1),Sds_Count(2)))

  status = nf90_inq_varid(ncid_vgac, 'proj_time0', varid)
  if (status /= nf90_noerr) then
        print *, "Error: Unable to get variable id for ", 'proj_time0'
        return
  endif
  status = nf90_get_var(ncid_vgac, varid, proj_time0)
  if (status /= nf90_noerr) then
        print *, "Error: Unable to read proj_time0"
        return
  endif

  !--- read geolocation
  call READ_AND_UNSCALE_NETCDF_2D(Ncid_Vgac, Sds_Start, Sds_Stride, Sds_Count, "lat", Sds_Data_2d)
  Nav%Lat_1b(1:Sds_Count(1),1:Sds_Count(2)) = Sds_Data_2d
  Nav%Lat = Nav%Lat_1b

  call READ_AND_UNSCALE_NETCDF_2D(Ncid_Vgac, Sds_Start, Sds_Stride, Sds_Count, "lon", Sds_Data_2d)
  where((Sds_Data_2d .ltr. -180.0) .and. (Sds_Data_2d .ner.  MISSING_VALUE_REAL4))  !need to flip
    Sds_Data_2d = Sds_Data_2d + 360.0
  endwhere
  where(Sds_Data_2d .gtr. 180.0)  !need to flip
    Sds_Data_2d = Sds_Data_2d - 360.0
  endwhere
  Nav%Lon_1b(1:Sds_Count(1),1:Sds_Count(2)) = Sds_Data_2d
  Nav%Lon = Nav%Lon_1b

  !--- angles
  call READ_AND_UNSCALE_NETCDF_2D(Ncid_Vgac, Sds_Start, Sds_Stride, Sds_Count, "vza", Sds_Data_2d)
  Geo%Satzen(1:Sds_Count(1),1:Sds_Count(2)) = Sds_Data_2d    !THIS IS NOT CORRECT FIXME

  call READ_AND_UNSCALE_NETCDF_2D(Ncid_Vgac, Sds_Start, Sds_Stride, Sds_Count, "sza", Sds_Data_2d)
  Geo%Solzen(1:Sds_Count(1),1:Sds_Count(2)) = Sds_Data_2d

  call READ_AND_UNSCALE_NETCDF_2D(Ncid_Vgac, Sds_Start, Sds_Stride, Sds_Count, "azi", Sds_Data_2d)
  Geo%Sataz(1:Sds_Count(1),1:Sds_Count(2)) = Sds_Data_2d

  call READ_AND_UNSCALE_NETCDF_2D(Ncid_Vgac, Sds_Start, Sds_Stride, Sds_Count, "azn", Sds_Data_2d)
  Geo%Solaz(1:Sds_Count(1),1:Sds_Count(2)) = Sds_Data_2d

  !--- Channel Observations
  do Chan_Idx = 1, Sensor%Num_Chan_Sensor

     Clavrx_Chan_Idx = Sensor%CLAVRx_Chan_Map(Chan_Idx)

     if (Sensor%Chan_On_Flag_Default(Clavrx_Chan_Idx) == sym%No) cycle

     call READ_AND_UNSCALE_NETCDF_2D(Ncid_Vgac, Sds_Start, Sds_Stride, Sds_Count, &
                                     Obs_Name(Chan_Idx), Sds_Data_2d)
     where((Sds_Data_2d .ler. 0.0))  !appears 0 is missing, there is no fill
         Sds_Data_2d = MISSING_VALUE_REAL4
     endwhere

     !--- solar channels
     if (ch(CLAVRx_Chan_Idx)%Obs_Type == SOLAR_OBS_TYPE) then
       !--- store reflectance
       ch(CLAVRx_Chan_Idx)%Ref_Toa(1:Sds_Count(1),1:Sds_Count(2)) = Sds_Data_2d
       !--- convert to 0-100%
       where(ch(CLAVRx_Chan_Idx)%Ref_Toa /= MISSING_VALUE_REAL4)
         ch(CLAVRx_Chan_Idx)%Ref_Toa = 100.0*ch(CLAVRx_Chan_Idx)%Ref_Toa
       endwhere
     endif

     !--- thermal and mixed channels
     if (ch(CLAVRx_Chan_Idx)%Obs_Type == THERMAL_OBS_TYPE .or. &
         ch(CLAVRx_Chan_Idx)%Obs_Type == MIXED_OBS_TYPE) then
       !--- store radiance
       ch(CLAVRx_Chan_Idx)%Rad_Toa(1:Sds_Count(1),1:Sds_Count(2)) = Sds_Data_2d
       !--- convert radiance to noaa units
       call CONVERT_RADIANCE (Ch(CLAVRx_Chan_Idx)%Rad_Toa, Planck_Nu(CLAVRx_Chan_Idx),MISSING_VALUE_REAL4)
       !--- convert radiance to brightness temperature
       call COMPUTE_BT_ARRAY(ch(CLAVRx_Chan_Idx)%Bt_Toa,ch(CLAVRx_Chan_Idx)%Rad_Toa,CLAVRx_Chan_Idx,MISSING_VALUE_REAL4)
     endif

  enddo 

  !--- read sub-pixel standard deviations
  if (Sensor%Chan_On_Flag_Default(1) == sym%Yes) then
     call READ_AND_UNSCALE_NETCDF_2D(Ncid_Vgac, Sds_Start, Sds_Stride, Sds_Count,'M05_sd', Sds_Data_2d)
     print *,"range in m05_sd = ", minval(Sds_Data_2d),maxval(Sds_Data_2d)
     where(Sds_Data_2d .ltr. 0.00)   !missing seems to be -1.00
         Sds_Data_2d = MISSING_VALUE_REAL4
     endwhere
     where(Sds_Data_2d .ger. 0.00)
         Sds_Data_2d = 100.0*Sds_Data_2d
     endwhere
     ch(1)%Ref_Toa_Std_Sub(1:Sds_Count(1),1:Sds_Count(2)) = Sds_Data_2d
  endif

  if (Sensor%Chan_On_Flag_Default(31) == sym%Yes) then  
     call READ_AND_UNSCALE_NETCDF_2D(Ncid_Vgac, Sds_Start, Sds_Stride, Sds_Count,'M15', Temp_Pix_Array_1(1:Sds_Count(1),1:Sds_Count(2)))
     where(Temp_Pix_Array_1 .ler. 0.00)   !missing seems to be -1.00
         Temp_Pix_Array_1 = MISSING_VALUE_REAL4
     endwhere
     call READ_AND_UNSCALE_NETCDF_2D(Ncid_Vgac, Sds_Start, Sds_Stride, Sds_Count,'M15_sd', Sds_Data_2d)
     where(Sds_Data_2d .ltr. 0.00)   !missing seems to be -1.00
         Sds_Data_2d = MISSING_VALUE_REAL4
     endwhere
     where(Sds_Data_2d .ger. 0.00)
        Sds_Data_2d = Sds_Data_2d + Temp_Pix_Array_1(1:Sds_Count(1),1:Sds_Count(2))
     endwhere
     where(Sds_Data_2d .ltr. 0.00) 
         Sds_Data_2d = MISSING_VALUE_REAL4
     endwhere
     call CONVERT_RADIANCE (Sds_Data_2d, Planck_Nu(31),MISSING_VALUE_REAL4)
     call COMPUTE_BT_ARRAY(Temp_Pix_Array_1(1:Sds_Count(1),1:Sds_Count(2)),Sds_Data_2d,31,MISSING_VALUE_REAL4)
     ch(31)%Bt_Toa_Std_Sub(1:Sds_Count(1),1:Sds_Count(2)) = Temp_Pix_Array_1(1:Sds_Count(1),1:Sds_Count(2)) - &
                                                            ch(31)%Bt_Toa(1:Sds_Count(1),1:Sds_Count(2))
  endif

  !--- fix vza
  do Elem_Idx = 1, Sds_Count(1)
    Geo%Satzen(Elem_Idx,:) =   70.0*abs(Elem_Idx - 401) / 800.0
  enddo

  !--- read scanline time
  call READ_NETCDF(Ncid_Vgac, [Sds_Start(2)], [Sds_Stride(2)], [Sds_Count(2)], "time", Sds_Data_1d)
  where(Sds_Data_1d > 0.0)
       Sds_Data_1d = 1000.0*60.0*60.0*Sds_Data_1d   !now in milliseconds
  endwhere
  Sds_Data_1d = Sds_Data_1d + ((proj_time0-floor(proj_time0))*24*60*60*1000)  !now in ms since start of day
  where(Sds_Data_1d > MSEC_PER_DAY)
    Sds_Data_1d = Sds_Data_1d - MSEC_PER_DAY
  endwhere

  Image%Scan_Time_Ms(1:Sds_Count(2)) = Sds_Data_1d

  !--- deallocate temporary memory
  deallocate(Sds_Data_1d)
  deallocate(Sds_Data_2d)


  !--- calculate other angles
  Geo%Relaz = RELATIVE_AZIMUTH (Geo%Solaz, Geo%Sataz)
  Geo%Glintzen = GLINT_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz )
  Geo%Scatangle = SCATTERING_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz)

  !--- There are no gaps in vgac
  Gap_Pixel_Mask = sym%NO

  !--- record number of scans read for this segment
  Image%Number_Of_Lines_Read_This_Segment = Sds_Count(2)

  !--- populate scan line number
  do Line_Idx = 1, Sds_Count(2)
    Image%Scan_Number(Line_Idx) = Ny_Start + Line_Idx
  enddo

  !--- fill space mask
  Geo%Space_Mask = .TRUE.
  where((Nav%Lat .ger. -90.0) .and. (Nav%Lat .ler. 90.0) .and. &
        (Nav%Lon .ger. -180.0) .and. (Nav%Lon .ler. 180.0))
        Geo%Space_Mask = .FALSE.
  endwhere

  !--- compute ascending descending flag
  Nav%Ascend(1:Image%Number_Of_Lines_Read_This_Segment) = 1
  do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment-1
     if (Nav%Lat(Image%Number_Of_Elements/2,Line_Idx+1) >  &
         Nav%Lat(Image%Number_Of_Elements/2,Line_Idx)) then
         Nav%Ascend(Line_Idx)= 0
     endif
  enddo
  Nav%Ascend(Image%Number_Of_Lines_Read_This_Segment) =  &
              Nav%Ascend(Image%Number_Of_Lines_Read_This_Segment-1)

  !--- close file on last segment if still open
  if (Image%Segment_Number == Image%Number_Of_Segments .and. Ncid_Vgac > 0) then
    call CLOSE_NETCDF(Ncid_Vgac)
  endif

end subroutine READ_VGAC_DATA

end module CX_VGAC_MOD
