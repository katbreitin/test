!$Id: fy3d_read_module.f90 3423 2019-07-17 21:54:29Z yli $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version
! 5.4
!
! NAME: fy3d_read_module.f90
!
! PURPOSE: FY3D (NetCDF4) read tool
!
! DESCRIPTION:  This module deals with reading FY3D data
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Yue Li, CIMSS, yue.li@ssec.wisc.edu
!  Steve Wanzong, CIMSS, steve.wanzong@ssec.wisc.edu
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
!--------------------------------------------------------------------------------------

module FY3D_READ_MODULE 

  use CX_READH5DATASET, only: &
   h5readattribute &
   , h5readdataset
  
 ! use FILE_UTILITY, only: &
 !     get_lun  
  use FILE_TOOLS, only: &
        FILE_SEARCH &
      , GET_LUN
  use PIXEL_COMMON_MOD, only: &
        Sensor &
      , Image &
      , Gap_Pixel_Mask &
      , Temp_Pix_Array_1
  use CONSTANTS_MOD, only: &
        Real4 &
      , Real8 &
      , Int8 &
      , Int4 &
      , Int2 &
      , Int1 &
      , Sym &
      , Exe_Prompt &
      , Missing_Value_Real4 &
      , Missing_Value_Int4 
   use CLAVRX_MESSAGE_MOD
   use CALIBRATION_CONSTANTS_MOD, only: &
        sat_name &
      , solar_ch20 &
      , solar_ch20_nu &
      , ew_ch20 &
      , planck_a1 &
      , planck_a2 &
      , planck_nu
   USE hdf5

  implicit none

  public :: READ_FY3D_DATA
  public :: READ_FY3D_DATE_TIME
  public :: READ_NUMBER_OF_SCANS_FY3D
  public :: READ_FY3D_INSTR_CONSTANTS
  private :: DETERMINE_FY3D_FILE
  private :: JULIAN
  private :: CONVERT_FY3D_RADIANCE

  character(len=18), parameter:: FY3D_PROMPT="FY3D_MODULE:"
 
  character(len=13), parameter:: MODULE_PROMPT="FY3D_MODULE:"

      contains

!--------------------------------------------------------------------
! read FY3D data
!--------------------------------------------------------------------
subroutine READ_FY3D_DATA (Segment_Number, GEO1K_File, Error_Out)

   use Pixel_Common_Mod , only : &
        Geo &
      , Nav &
      , Ancil_Data_Dir &
      , Ch &
      , Bt_11um_Sounder &
      , Bt_12um_Sounder 

   use PLANCK_MOD

   use VIEWING_GEOMETRY_MOD, only: &
        GLINT_ANGLE &
      , SCATTERING_ANGLE &
      , RELATIVE_AZIMUTH

   use CALIBRATION_CONSTANTS_MOD, only: &
        Planck_Nu &
      , Planck_Nu_11um_Sndr &
      , Planck_Nu_12um_Sndr

  ! USE hdf5

      integer(kind=int4), intent(in):: Segment_Number
      character(len=*), intent(in):: GEO1K_File
      integer(kind=int4), intent(out):: Error_Out

      integer(kind=int4) :: Nx_Start , Nx_End , Ny_Start , Ny_End , N_Seg_Lines
      integer(kind=int4) :: I_Geo
      integer(kind=int4) :: I_band
      integer(kind=int4) :: N_band
      integer(kind=int4) :: k , i
      integer(kind=int4) :: Lun
      integer(kind=int4) :: Io_Err_Stat
      integer(kind=int4), dimension(2) :: Dim_Seg
      integer(kind=int4), dimension(3) :: Dim_Seg1
      integer(kind=int4), dimension(2) :: Offset_band
      integer(kind=int4), dimension(3) :: Offset_band1
      integer(kind=int4), dimension(:), allocatable :: Time_Sec_Day
      integer(kind=int4), dimension(7) :: Scaled_Geo = [0, 0, 0, 1, 1, 1, 1]
      integer(kind=int4), dimension(24) :: Modis_Chn_List
      integer(kind=int4), dimension(5) :: Modis_I_Chn_List
      real(kind=real4) :: Fill_Value
      real(kind=real4) :: Scale_Factor , Add_Offset
      real(kind=real4), dimension(:), pointer :: Scale_Factor1 , Add_Offset1
      !real(kind=real4), dimension(4) :: Scale_Factor1, Add_Offset1
      real(kind=real8), parameter :: Sec_Per_Day = 86400.
      real(kind=real8), dimension(:), pointer :: R1d_Buffer
      real(kind=real4), dimension(:,:), pointer :: R2d_Buffer
      real(kind=real4), dimension(:,:), allocatable :: R2d_Temp
      integer(kind=int4), dimension(:,:), pointer :: I2d_Buffer
      integer(kind=int4), dimension(:,:,:), pointer :: I2d_Buffer1
      real(kind=real4), dimension(:,:,:), pointer :: R2d_Buffer1
      character(len=1020) :: File_2_Read
      character(len=5) :: Day_Night_Flag
      character(len=100), dimension ( 7 ) :: Setname_Geo_List = (/ character (len =300) :: &
                           'Geolocation/Latitude             ' & ! 1
                         , 'Geolocation/Longitude            ' & ! 2
                         , 'Timedata/Millisecond_Count  ' & ! 3
                       !  , 'scan_line_attributes/scan_start_time  ' & ! 3
                         , 'Geolocation/SensorAzimuth       ' & ! 4
                         , 'Geolocation/SensorZenith        ' & ! 5
                         , 'Geolocation/SolarAzimuth        ' & ! 6
                         , 'Geolocation/SolarZenith         ' & ! 7
                                                    /)
      character(len=100) :: Setname_Name
      logical, dimension(24) :: Is_band_On
      logical, dimension(24) :: Is_band_Read
      real(kind=real4), dimension(:,:,:), allocatable :: Rad_Coef_IR
      real(kind=real4), dimension(:,:), allocatable :: Rad_Coef_VIS

!!!
      real(kind=4),    dimension( 200, 4,  6) :: var_real4_ir
      integer ::   error 
      integer (HSIZE_T), dimension(2)         :: dims_geo2
      integer (HID_T)  :: file_id              ! file id 
      integer (HID_T) :: sds_id_var           ! sds id all variables
      !real(kind=4) :: H5T_NATIVE_REAL

      Error_Out = 0

      ! --- calculate start and end segment limits

      Nx_Start = 1
      Nx_End = Nx_Start + Image%Number_Of_Elements - 1

      Ny_Start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
      Ny_End = min(Image%Number_Of_Lines, Ny_Start + Image%Number_of_Lines_Per_Segment - 1)
      N_Seg_Lines = Ny_End - Ny_Start + 1


      Offset_band = [Nx_Start - 1, Ny_Start - 1]
      Dim_Seg = [Nx_End - Nx_Start + 1, Ny_End - Ny_Start + 1]
      
      ! --- read geo file data
      ! loop over geo data
      do I_Geo = 1, 7
        if ( I_Geo == 3 ) then
           call H5READDATASET (trim(Image%Level1b_Path)//trim(GEO1K_File), &
                      trim(Setname_Geo_List(I_Geo)), R2d_Buffer)
        else
           call H5READDATASET (trim(Image%Level1b_Path)//trim(GEO1K_File), &
                      trim(Setname_Geo_List(I_Geo)), Offset_band, Dim_Seg, R2d_Buffer)
           if (Scaled_Geo(I_Geo) .eq. 1) then
              call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(GEO1K_File), &
                      trim(Setname_Geo_List(I_Geo)//'/Slope'), Scale_Factor)
              call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(GEO1K_File), &
                      trim(Setname_Geo_List(I_Geo)//'/Intercept'), Add_Offset)
              call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(GEO1K_File), &
                      trim(Setname_Geo_List(I_Geo)//'/FillValue'), Fill_Value)
              where ( R2d_Buffer .ne. Fill_Value)
                 R2d_Buffer = (R2d_Buffer * Scale_Factor) + Add_Offset
              endwhere
           endif
        endif

        ! save read data to output in correct format
        select case (I_Geo)
         case(1)
            Nav % Lat_1b (:,1:N_Seg_Lines) = R2d_Buffer
         case(2)
            Nav % Lon_1b (:,1:N_Seg_Lines) = R2d_Buffer
         case(3)
            allocate (Time_Sec_Day (size(R2d_Buffer(1,:))))
            Time_Sec_Day = int((mod(R2d_Buffer(1,:)/1000, Sec_Per_Day)) * 1000)

            ! make data missing of missing scan time
            where (R2d_Buffer(1,:) < 0)
               Time_Sec_Day = Missing_Value_Int4
            endwhere

            Image%Scan_Time_Ms(1:N_Seg_Lines) = (/(Time_Sec_Day((k - 1) / 16 + 1), &
                                              k = Ny_Start, Ny_End)/)
            deallocate ( Time_Sec_Day )
         case(4)
            where (R2d_Buffer > 180)
                 R2d_Buffer = R2d_Buffer - 360.
            endwhere
            Geo % Sataz (:,1:N_Seg_Lines) = R2d_Buffer
         case(5)
            Geo % Satzen (:,1:N_Seg_Lines) = R2d_Buffer
         case(6)
            where (R2d_Buffer > 180)
                 R2d_Buffer = R2d_Buffer - 360.
            endwhere
            Geo % Solaz (:,1:N_Seg_Lines) = R2d_Buffer
         case(7)
            Geo % Solzen (:,1:N_Seg_Lines) = R2d_Buffer
         end select

         deallocate ( R2d_Buffer)
         !if (I_Geo == 3) deallocate ( R1d_Buffer)

      enddo

      ! --- calculate ascending/descending
      Nav % Ascend (1:N_Seg_Lines) = 0
      do i = 1, N_Seg_Lines - 1
         if (Nav%Lat_1b (Dim_Seg(1) / 2, i + 1) <= Nav%Lat_1b(Dim_Seg(1) / 2, i)) &
                    Nav % Ascend( i ) = 1
      end do
      ! --- fix for the last line
      if ( Nav%Lat_1b(Nx_End / 2, N_Seg_Lines) <= &
           Nav%Lat_1b(Nx_End / 2, N_Seg_Lines - 1) ) &
                Nav % Ascend (N_Seg_Lines) = 1

      ! --- calculate relative azimuth
      Geo % Relaz = RELATIVE_AZIMUTH (Geo%Solaz, Geo%Sataz)

      !--- compute the glint zenith angle
      Geo % Glintzen = GLINT_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz )

      !--- compute the scattering angle
      Geo % Scatangle = SCATTERING_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz)


      ! --- read band data
      ! find band file
      if (Sensor%Platform_Name == 'FY3D') then
        call DETERMINE_FY3D_FILE(trim(Image%Level1b_Path),trim(GEO1K_File), &
                                    '1000M',File_2_Read)
      endif


      ! --- if no file quit
      if (trim(File_2_Read) .eq. 'no_file') then
         Error_Out = 1
         return
      endif

      ! --- read global attribute Day_Night_Flag
      !!!!! Note: if Day_Night_Flag = 'Night' not read in visible bands (1-18)  !!!!!
      call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                                       'Day Or Night Flag', Day_Night_Flag)


      ! --- mapping modis to fy3d
      !                 047 055 065 086 138 164 213 041 044 049 056 067 075 087 091 0936 094 103 375 405 720 855 108 120
      !  Band            1   2   3   4   5   6   7   8   9  10  11  12  14  15  16  17   18  19  20  21  22  23  24  25 
      Modis_Chn_List = [ 3,  4,  1,  2, 26,  6,  7,  8,  9, 10, 12, 13, 15, 16, 17, 18,  19,  5, 20, 23, 28, 29, 31, 32]
     !Modis_Chn_List = [3,4,1,2,26,6,7,8,9,10,      15,   17,18,19,5,20,23,28,29,31,32]

      Is_band_On = Sensor%Chan_On_Flag_Default (Modis_Chn_List) == sym%YES
      N_band = 24 

      Setname_Name = "Calibration/IR_Cal_Coeff"
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                             Setname_Name, R2d_Buffer1)

      if (.not. allocated(Rad_Coef_IR))  allocate(Rad_Coef_IR(200,4,6))
      Rad_Coef_IR = R2d_Buffer1
      deallocate (R2d_Buffer1)

      Setname_Name = "Calibration/VIS_Cal_Coeff"
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                             Setname_Name, R2d_Buffer)
      if (.not. allocated(Rad_Coef_VIS))  allocate(Rad_Coef_VIS(3,19))

      Rad_Coef_VIS = R2d_Buffer
      deallocate (R2d_Buffer)


! Chan 1-4 (FY3D Ch 1-4, MODIS Ch 3,4,1,2)
! Chan 5-18 (FY3D Ch 5-12, 14-19, MODIS Ch 26, 6-10, 12, 13, 15-19, 5)
! Note Chan13 is not read in due to unavailable MODIS band
! Chan 19-22 (FY3D Ch 20-23, MODIS Ch 20, 23, 28,29)
! Chan 23-24 (FY3D Ch 24-25, MODIS Ch31, 32)

      do I_band = 1,N_band
         if (I_band >=1 .and. I_band <=4) then
            setname_Name = "Data/EV_250_Aggr.1KM_RefSB"
         elseif (I_band >=5 .and. I_band <=18 ) then
            setname_Name = "Data/EV_1KM_RefSB"
         elseif (I_band >=19 .and. I_band <=22) then
            setname_Name = "Data/EV_1KM_Emissive"
         elseif (I_band >=23 .and. I_band <=24) then
            setname_Name = "Data/EV_250_Aggr.1KM_Emissive"
         endif

         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/Slope', Scale_Factor1)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/Intercept', Add_Offset1)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/FillValue', Fill_Value)
         Offset_band1 = [Nx_Start - 1, Ny_Start - 1, 0]
         Dim_Seg1 = [Nx_End - Nx_Start + 1, Ny_End - Ny_Start + 1, shape(Scale_Factor1)]
      !
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                             Setname_Name, Offset_band1, Dim_Seg1, I2d_Buffer1)

         ! - check if channel is on 
         Is_band_Read(I_band) = .false.
         if (.not. Is_band_On(I_band)) cycle

         allocate (R2d_Temp(Dim_Seg(1),Dim_Seg(2)))
         R2d_Temp = Missing_Value_Real4

! VIS bands
       if (Day_Night_Flag /= 'N') then
         if (I_band >=1 .and. I_band <=4) then
             Ch (Modis_Chn_List(I_band)) % Ref_Toa (:, 1:N_Seg_Lines) = Missing_Value_Real4
             R2d_Temp =  float(I2d_Buffer1(:,:,I_band)) 
             ! - unscale
             where ( I2d_Buffer1(:,:,I_band) .ne. Fill_Value)
                   Ch (Modis_Chn_List(I_band)) % Ref_Toa (:, 1:N_Seg_Lines) = &
!                                 ( (I2d_Buffer1(:,:,I_band) * Scale_Factor1(I_band)) + Add_Offset1(I_band) )*100
                          Rad_Coef_VIS(1,I_band)             +Rad_Coef_VIS(2,I_band)*R2d_Temp+ &
                          Rad_Coef_VIS(3,I_band)*R2d_Temp**2
             endwhere
         elseif (I_band >= 5 .and. I_band <=18) then
             Ch (Modis_Chn_List(I_band)) % Ref_Toa (:, 1:N_Seg_Lines) = Missing_Value_Real4

             ! - unscale
             if (I_band < 13) then   ! number 1 - 8
                R2d_Temp =  float(I2d_Buffer1(:,:,I_band-4))
                where ( I2d_Buffer1(:,:,I_band-4) .ne. Fill_Value)
                       Ch (Modis_Chn_List(I_band)) % Ref_Toa (:, 1:N_Seg_Lines) = &
!                                    ( (I2d_Buffer1(:,:,I_band-4) * Scale_Factor1(I_band-4)) + Add_Offset1(I_band-4) )*100
                          Rad_Coef_VIS(1,I_band)             + Rad_Coef_VIS(2,I_band)*R2d_Temp+ &
                          Rad_Coef_VIS(3,I_band)*R2d_Temp**2
                endwhere
             endif

             if (I_band >= 13) then  ! number 10 - 15
                R2d_Temp =  float(I2d_Buffer1(:,:,I_band-3))
                where ( I2d_Buffer1(:,:,I_band-3) .ne. Fill_Value)
                      Ch(Modis_Chn_List(I_band)) % Ref_Toa (:, 1:N_Seg_Lines) = &
!                                    ( (I2d_Buffer1(:,:,I_band) * Scale_Factor1(I_band)) + Add_Offset1(I_band) )*100
                          Rad_Coef_VIS(1,I_band+1)             + Rad_Coef_VIS(2,I_band+1)*R2d_Temp+ &
                          Rad_Coef_VIS(3,I_band+1)*R2d_Temp**2
                endwhere
             endif
         endif
       endif

! IR bands
         if (I_band >=19 .and. I_band <=24) then

             if (I_band >=19 .and. I_band <=22) then
                Ch (Modis_Chn_List(I_band)) % Rad_Toa (:, 1:N_Seg_Lines) = Missing_Value_Real4
                where ( I2d_Buffer1(:,:,I_band-18) .ne. Fill_Value)
                      R2d_Temp =  (I2d_Buffer1(:,:,I_band-18)+Add_Offset1(I_band-18)) * Scale_Factor1(I_band-18) 
                      Ch(Modis_Chn_List(I_band)) % Rad_Toa (:, 1:N_Seg_Lines) = &
                         R2d_Temp
                endwhere

             elseif (I_band >=23 .and. I_band <=24) then
                Ch (Modis_Chn_List(I_band)) % Rad_Toa (:, 1:N_Seg_Lines) = Missing_Value_Real4
                where ( I2d_Buffer1(:,:,I_band-22) .ne. Fill_Value)
                      R2d_Temp =  (I2d_Buffer1(:,:,I_band-22)+Add_Offset1(I_band-22)) *  Scale_Factor1(I_band-22) 
                      Ch(Modis_Chn_List(I_band)) % Rad_Toa (:, 1:N_Seg_Lines) = &
                         R2d_Temp    ! direct seems work (opt = 21)
                endwhere
       
             endif

        !     deallocate (R2d_Temp)

             call CONVERT_FY3D_RADIANCE (Ch(Modis_Chn_List(I_band))%Rad_Toa(:,1:N_Seg_Lines), &
                             Planck_Nu(Modis_Chn_List(I_band)), Missing_Value_Real4)

             ! - calculate bt
             call COMPUTE_BT_ARRAY (Ch(Modis_Chn_List(I_band))%Bt_Toa(:,1:N_Seg_Lines), &
                                    Ch(Modis_Chn_List(I_band))%Rad_Toa(:,1:N_Seg_Lines), &
                                    Modis_Chn_List(I_band), Missing_Value_Real4)

         endif

        deallocate (R2d_Temp)
        Is_band_Read(I_band) = .true.
        deallocate (I2d_Buffer1)

      enddo

! open geo hdf5 file
!call h5open_f(error)
!call h5fopen_f (trim(Image%Level1b_Path)//trim(File_2_Read), H5F_ACC_RDONLY_F, file_id, error)

!call h5dopen_f (file_id, '/Calibration/IR_Cal_Coeff', sds_id_var, error)
!call h5dread_f (sds_id_var, H5T_NATIVE_REAL, var_real4_ir, dims_geo2, error)

      ! --- global variables which have to be set
      Image%Number_Of_Lines_Read_This_Segment = N_Seg_Lines
      do i = 1, Image%Number_Of_Lines_Per_Segment
         Image%Scan_Number(i) = Ny_Start + i - 1
      end do

end subroutine READ_FY3D_DATA

!--------------------------------------------------------------------
! read fy3d time
!--------------------------------------------------------------------
subroutine READ_FY3D_DATE_TIME (Path, Infile, Year , Doy , Start_Time &
                , End_Time , Orbit , End_Year, End_Doy )

      character(len=*), intent(in) :: Path
      character(len=*), intent(in) :: Infile
      integer, intent(out) :: Year
      integer, intent(out) :: Doy    !day of year
      integer, intent(out) :: Start_Time  !millisec
      integer, intent(out) :: End_Time    !millisec
      integer, intent(out) :: Orbit
      integer , intent(out) :: End_Year
      integer, intent(out) :: End_Doy    !day of year

      character(len=35) :: Tmp_String 
      integer :: Month
      integer :: Day
      integer :: Start_Hour
      integer :: Start_Minute
      integer :: Start_Sec

      integer :: End_Hour
      integer :: End_Minute
      integer :: End_Sec


      ! --- read time and date from the attributes
      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'Orbit Number', Orbit)

      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'Observing Beginning Time', Tmp_String)
      !0        1         2
      !123456789012345678901234567890
      !2016-04-20T12:00:00.000Z
      read(Tmp_String(1:2), fmt="(I2)") Start_Hour
      read(Tmp_String(4:5), fmt="(I2)") Start_Minute
      read(Tmp_String(7:8), fmt="(I2)") Start_Sec

     call H5ReadAttribute( trim(Path)//trim(Infile), &
           'Observing Beginning Date', Tmp_String)

      read(Tmp_String(1:4), fmt="(I4)") Year
      read(Tmp_String(6:7), fmt="(I2)") Month
      read(Tmp_String(9:10), fmt="(I2)") Day

      !--- compute start day of year
      call JULIAN ( Day, Month, Year, Doy )
      
      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'Observing Ending Time', Tmp_String)
      read(Tmp_String(1:2), fmt="(I2)") End_Hour
      read(Tmp_String(4:5), fmt="(I2)") End_Minute
      read(Tmp_String(7:8), fmt="(I2)") End_Sec

     call H5ReadAttribute( trim(Path)//trim(Infile), &
           'Observing Ending Date', Tmp_String)

      read(Tmp_String(1:4), fmt="(I4)") End_Year
      read(Tmp_String(6:7), fmt="(I2)") Month
      read(Tmp_String(9:10), fmt="(I2)") Day

      !--- compute end day of year
      call JULIAN ( Day, Month, End_Year, End_Doy )

      ! --- Calculate start and end time
      Start_Time = ((Start_Hour * 60 + Start_Minute) * 60 + Start_Sec) * 1000
      End_Time = ((End_Hour * 60 + End_Minute) * 60 + End_Sec) * 1000

end subroutine READ_FY3D_DATE_TIME

!--------------------------------------------------------------------
! read fy3d number of scans
!--------------------------------------------------------------------
subroutine READ_NUMBER_OF_SCANS_FY3D (Infile, Number_Of_Fy3d_Lines, &
                 Error_Out)

      character(len=*), intent(in) :: Infile
      integer(kind=int4), intent(out) :: Error_Out
      integer(kind=int4), intent(out) :: Number_of_Fy3d_Lines

      character(len=100) :: Setname
      integer, dimension(:), pointer :: Test
      integer, dimension(1) :: Dims

      integer :: Number_of_Scans     


      Error_Out = 0
      !Setname = 'scan_line_attributes/scan_quality'

      call H5ReadAttribute( trim(Infile), &
           'Number Of Scans', Number_of_Scans)
!      call H5ReadDataset( trim(Infile), trim(Setname), Test )
      Dims = Number_of_Scans
!print *,'Dims ',Dims
      ! --- error 
      if (Dims(1) .lt. 1) then
         Number_Of_Fy3d_Lines = -999
         Error_Out = 1
         Return
      endif
         
      Number_Of_Fy3d_Lines = Dims(1) * 10

end subroutine READ_NUMBER_OF_SCANS_FY3D

!--------------------------------------------------------------------
! find fy3d files
!--------------------------------------------------------------------
subroutine DETERMINE_FY3D_FILE(Path_In,File_In,File_Type_In,File_Out)

      character(len=*), intent(in):: Path_In
      character(len=*), intent(in):: File_In
      character(len=*), intent(in):: File_Type_In
      character(len=*), intent(out):: File_Out

      character(len=500):: Search_String
      character(len=1020), dimension(:), pointer:: Files
      integer(kind=int4):: Num_Files
      integer:: Search_len

      logical, save:: First_Call = .true.


      if (index(File_Type_In, '1000M') > 0) then
        if (trim(File_In(1:2)) == 'FY') then
           Search_len = 33
        else if (trim(File_In(1:2)) == 'tf') then
           Search_len = 29
        endif
         Search_String = trim(File_In(1:Search_len))//trim(File_Type_In)//'*.HDF'
      endif


      Files => FILE_SEARCH(trim(Path_In),trim(Search_String),count=Num_Files)

      if (Num_Files == 0) then
         !print *, EXE_PROMPT, FY3D_PROMPT, " No "//trim(File_Type_In)//" File Found"
         !call MESG(" No "//trim(File_Type_In)//" File Found", level = verb_lev % DEFAULT)
         File_Out = "no_file"
         return
      endif

      if (Num_Files > 1) then
         !print *, EXE_PROMPT, FY3D_PROMPT, "Multiple "//trim(File_Type_In)//" Files Found"
         !call MESG(" Multiple "//trim(File_Type_In)//" Files Found", level = verb_lev % DEFAULT)
!         File_Out = "no_file"
      endif

      File_Out = Files(1)

      if (First_Call .eqv. .true.) then
         print *, EXE_PROMPT, FY3D_PROMPT, trim(File_Type_In)//" File Found, ",trim(File_Out)
      endif

      First_Call = .false.

      Files => null()

end subroutine DETERMINE_FY3D_FILE

!-------------------------------------------------
! subroutine JULIAN(iday,imonth,iyear,jday)
! compute julian day
! input:
!         iday - integer day
!         imonth - integer month
!         iyear - integer year (2 or four digits)
!         
! output : jday - julian day
!--------------------------------------------------
subroutine JULIAN(iday,imonth,iyear,jday)

        integer, intent(in)::  iday,imonth,iyear
        integer, intent(out):: jday
        integer::  j
        integer, dimension(12)::  jmonth

        jmonth = reshape ((/31,28,31,30,31,30,31,31,30,31,30,31/),(/12/))

        jday = iday
        if (modulo(iyear,4) == 0) then
            jmonth(2)=29
        endif

        do j = 1,imonth-1
           jday = jday + jmonth(j)
        end do

end subroutine JULIAN

!--------------------------------------------------
! Function Name: CONVERT_FY3D_RADIANCE
!
! Function:
!    Convert units of the FY3D radiance values to that expected by CLAVR-x.
!    This call is used though fy3d uses the the noaa format for IR radiance

! Description: 
!   
! Calling Sequence: rad_new =
! convert_fy3d_radiance(rad_old,nu,missing_value)
!   
!
! Inputs:
!   rad_old = radiance in units of mW/m^2/cm^-1/str (2d array)
!   nu = channels equivalent width in units of cm^-1
!   missing_value = value assigned to missing radiance values
!
! Outputs: 
!   rad_new = radiance in units of mW/m^2/cm^-1/str (2d array)
!
! Dependencies:
!
! Restrictions:  None
!
! Reference: algebraic manipulation of Planck Equation
! ---------------------------------------------------------------------------------------
subroutine CONVERT_FY3D_RADIANCE(Radiance,Nu,Missing_Value)
      real (kind=real4), dimension(:,:), intent(inout):: Radiance
      real (kind=real4), intent(in):: Nu
      real (kind=real4), intent(in):: Missing_Value

      where(Radiance /= Missing_Value)
!         Radiance = Radiance * (((10000.0 / Nu )**2) / 10.0)
         Radiance = Radiance !/ 10000.
      end where

      return

end subroutine CONVERT_FY3D_RADIANCE


!----------------------------------------------------------------
! read the fy3d constants into memory
!-----------------------------------------------------------------
   subroutine READ_FY3D_INSTR_CONSTANTS(Instr_Const_File)
     character(len=*), intent(in):: Instr_Const_File
     integer:: ios0, erstat
     integer:: Instr_Const_lun

     Instr_Const_lun = GET_LUN()

     open(unit=Instr_Const_lun,file=trim(Instr_Const_File),status="old",position="rewind",action="read",iostat=ios0)

     print *, EXE_PROMPT, MODULE_PROMPT, " Opening ", trim(Instr_Const_File)
     erstat = 0
     if (ios0 /= 0) then
        erstat = 19
        print *, EXE_PROMPT, MODULE_PROMPT, "Error opening FY3D constants file, ios0 = ", ios0
        stop 19
     endif

     read(unit=Instr_Const_lun,fmt="(a4)") sat_name
     read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
     read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
     read(unit=Instr_Const_lun,fmt=*) planck_a1(20), planck_a2(20), planck_nu(20)
     read(unit=Instr_Const_lun,fmt=*) planck_a1(23), planck_a2(23), planck_nu(23)
     read(unit=Instr_Const_lun,fmt=*) planck_a1(28), planck_a2(28), planck_nu(28)
     read(unit=Instr_Const_lun,fmt=*) planck_a1(29), planck_a2(29), planck_nu(29)
     read(unit=Instr_Const_lun,fmt=*) planck_a1(31), planck_a2(31), planck_nu(31)
     read(unit=Instr_Const_lun,fmt=*) planck_a1(32), planck_a2(32), planck_nu(32)
     close(unit=Instr_Const_lun)

     !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
     Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

   end subroutine READ_FY3D_INSTR_CONSTANTS


!subroutine READ_FY3D_LEVEL1B(path,file_name,iband, &
!                              calibrated_data_out, &
!                              uncert_data_out, &
!                              nx,ny,Seg_Idx,ny_total,ny_local_temp, &
!                              Error_Status)

!        Setname_name = "Data/EV_1KM_Emissive"
         ! - read data
!         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
!                             Setname_Name, Offset_band, Dim_Seg, I2d_Buffer)

         ! - read attribute to unscale
!         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
!                              Setname_Name//'/Slope', Scale_Factor)
!         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
!                              Setname_Name//'/Intercept', Add_Offset)
!         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
!                              Setname_Name//'/FillValue', Fill_Value)
!end subroutine READ_FY3D_LEVEL1B

!====================================================================


end module FY3D_READ_MODULE 
