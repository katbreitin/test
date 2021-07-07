!$Id: viirs_nasa_read_module.f90 1541 2016-02-26 23:18:41Z dbotambekov $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version
! 5.4
!
! NAME: viirs_nasa_read_module.f90
!
! PURPOSE: VIIRS NASA (NetCDF4) read tool
!
! DESCRIPTION:  This module deals with reading VIIRS NASA data
!
! AUTHORS:
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
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
! HISTORY:   created         May, 2016 (Denis B.)
!            added Fusion    April, 2018 (Denis B.)
!            added I-bands   August, 2018 (Denis B.)
!--------------------------------------------------------------------------------------

module VIIRS_NASA_READ_MODULE 

  use CX_READH5DATASET, only: &
   h5readattribute &
   , h5readdataset
  
  
  use FILE_TOOLS, only: &
        FILE_SEARCH &
      , GETLUN
  use PIXEL_COMMON_MOD, only: &
        Sensor &
      , Image &
      , Gap_Pixel_Mask &
      , Temp_Pix_Array_1 &
      , Temp_Pix_Array_2 &
      , X_Sample_Offset &
      , Y_Sample_Offset &
      , Use_Iband
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

  implicit none

  public :: READ_VIIRS_NASA_DATA
  public :: READ_VIIRS_NASA_DATE_TIME
  public :: READ_NUMBER_OF_SCANS_VIIRS_NASA
  public :: CHECK_IF_FUSION
  public :: READ_FUSION_INSTR_CONSTANTS
  private :: DETERMINE_VIIRS_NASA_FILE
  private :: JULIAN
  private :: CONVERT_VIIRS_NASA_RADIANCE
  private :: CONVERT_RAD_2_SOL_REF_DNB
  private :: COMPUTE_IBAND_STATISTICS
  private :: COMPUTE_VIIRS_BOWTIE_GAP_PATTERN
  private :: FILL_VIIRS_BOWTIE_GAPS

  character(len=18), parameter:: VIIRS_NASA_PROMPT="VIIRS_NASA_MODULE:"
 
  ! - bowtie gaps values
  integer, parameter :: Ny_Pattern = 16
  integer, parameter :: Nx_Pattern = 3200
  integer(kind=int4), dimension(:,:), allocatable, public :: Gap_Line_Idx_Pattern
  integer(kind=int1), dimension(:,:), allocatable, public :: Gap_Pixel_Mask_Pattern

      contains

!--------------------------------------------------------------------
! read viirs nasa data
!--------------------------------------------------------------------
subroutine READ_VIIRS_NASA_DATA (Segment_Number, VGEOM_File, Error_Out)

   use Pixel_Common_Mod , only : &
        Geo &
      , Nav &
      , Ancil_Data_Dir &
      , Ch &
      , Bt_11um_Sounder &
      , Bt_12um_Sounder 
 
   use PLANCK_MOD, only: &
    compute_bt_array &
    , compute_rad_array

   use VIEWING_GEOMETRY_MOD, only: &
        GLINT_ANGLE &
      , SCATTERING_ANGLE &
      , RELATIVE_AZIMUTH

   use CALIBRATION_CONSTANTS_MOD, only: &
        Planck_Nu &
      , Planck_Nu_11um_Sndr  &
      , Planck_Nu_12um_Sndr  &
      , VIIRS_Correction_Factor &
      , Sun_Earth_Distance

      integer(kind=int4), intent(in):: Segment_Number
      character(len=*), intent(in):: VGEOM_File
      integer(kind=int4), intent(out):: Error_Out

      integer(kind=int4) :: Nx_Start , Nx_End , Ny_Start , Ny_Start_Iband , Ny_End , N_Seg_Lines
      integer(kind=int4) :: I_Geo
      integer(kind=int4) :: Mband_Start
      integer(kind=int4) :: I_Mband
      integer(kind=int4) :: N_Mband
      integer(kind=int4) :: I_Iband
      integer(kind=int4) :: N_Iband
      integer(kind=int4) :: I_Fusion
      integer(kind=int4) :: N_Fusion
      integer(kind=int4) :: Fusion_Start
      integer(kind=int4) :: k , i
      integer(kind=int4) :: Lun
      integer(kind=int4) :: Io_Err_Stat
      integer(kind=int8) :: N_Valid
      integer(kind=int4), dimension(2) :: Dim_Seg
      integer(kind=int4), dimension(2) :: Dim_Seg_Iband
      integer(kind=int4), dimension(2) :: Dim_Dnb_Seg
      integer(kind=int4), dimension(2) :: Offset_Mband
      integer(kind=int4), dimension(2) :: Offset_Iband
      integer(kind=int4), dimension(3200) :: D2m_Idx
      integer(kind=int4), dimension(:), allocatable :: Time_Sec_Day
      integer(kind=int4), dimension(7) :: Scaled_Geo = [0, 0, 0, 1, 1, 1, 1]
      integer(kind=int4), dimension(16) :: Modis_Chn_List
      integer(kind=int4), dimension(5) :: Modis_I_Chn_List
      real(kind=real4) :: Fill_Value
      real(kind=real4) :: Scale_Factor , Add_Offset
      real(kind=real8), parameter :: Sec_Per_Day = 86400.
      real(kind=real8), dimension(:), pointer :: R1d_Buffer
      real(kind=real4), dimension(:,:), pointer :: R2d_Buffer
      real(kind=real4), dimension(:,:), allocatable :: R2d_Temp
      real(kind=real4), dimension(:,:), allocatable :: Iband_Temp_Buff
      real(kind=real4), dimension(:,:), allocatable :: Iband_Temp_Buff_Bt
      integer(kind=int4), dimension(:,:), pointer :: I2d_Buffer
      character(len=1020) :: File_2_Read
      character(len=1020) :: File_Dnb_Idx
      character(len=2) :: Band_Num_Str
      character(len=2) :: Fusion_Band_Str
      character(len=5) :: Day_Night_Flag
      character(len=100), dimension ( 7 ) :: Setname_Geo_List = (/ character (len =300) :: &
                           'geolocation_data/latitude             ' & ! 1
                         , 'geolocation_data/longitude            ' & ! 2
                         , 'scan_line_attributes/scan_start_time  ' & ! 3
                         , 'geolocation_data/sensor_azimuth       ' & ! 4
                         , 'geolocation_data/sensor_zenith        ' & ! 5
                         , 'geolocation_data/solar_azimuth        ' & ! 6
                         , 'geolocation_data/solar_zenith         ' & ! 7
                                                    /)
      character(len=100) :: Setname_Name
      logical, dimension(16) :: Is_Mband_On
      logical, dimension(16) :: Is_Mband_Read
      logical, dimension(5) :: Is_Iband_On
      logical, dimension(5) :: Is_Iband_Read
      logical :: Read_Iband_Flag

      Error_Out = 0

      ! --- calculate start and end segment limits
      !Number_Of_Lines_Per_Segment
      !Number_Of_Lines_Read_This_Segment
      !Number_Of_Lines

      Nx_Start = 1
      Nx_End = Nx_Start + Image%Number_Of_Elements - 1

      Ny_Start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
      Ny_End = min(Image%Number_Of_Lines, Ny_Start + Image%Number_of_Lines_Per_Segment - 1)
      N_Seg_Lines = Ny_End - Ny_Start + 1


      Offset_Mband = [Nx_Start - 1, Ny_Start - 1]
      Dim_Seg = [Nx_End - Nx_Start + 1, Ny_End - Ny_Start + 1]
            
      ! --- read geo file data
      ! loop over geo data
      do I_Geo = 1, 7
        if ( I_Geo == 3 ) then
           call H5READDATASET (trim(Image%Level1b_Path)//trim(VGEOM_File), &
                      trim(Setname_Geo_List(I_Geo)), R1d_Buffer)
        else
           call H5READDATASET (trim(Image%Level1b_Path)//trim(VGEOM_File), &
                      trim(Setname_Geo_List(I_Geo)), Offset_Mband, Dim_Seg, R2d_Buffer)
           if (Scaled_Geo(I_Geo) .eq. 1) then
              call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(VGEOM_File), &
                      trim(Setname_Geo_List(I_Geo)//'/scale_factor'), Scale_Factor)
              call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(VGEOM_File), &
                      trim(Setname_Geo_List(I_Geo)//'/add_offset'), Add_Offset)
              call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(VGEOM_File), &
                      trim(Setname_Geo_List(I_Geo)//'/_FillValue'), Fill_Value)
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
            allocate (Time_Sec_Day (size(R1d_Buffer)))
            Time_Sec_Day = int((mod(R1d_Buffer, Sec_Per_Day)) * 1000)

            ! make data missing of missing scan time
            where (R1d_Buffer < 0)
               Time_Sec_Day = Missing_Value_Int4
            endwhere
            Image%Scan_Time_Ms(1:N_Seg_Lines) = (/(Time_Sec_Day((k - 1) / 16 + 1), &
                                              k = Ny_Start, Ny_End)/)
            deallocate ( Time_Sec_Day )
         case(4)
            Geo % Sataz (:,1:N_Seg_Lines) = R2d_Buffer
         case(5)
            Geo % Satzen (:,1:N_Seg_Lines) = R2d_Buffer
         case(6)
            Geo % Solaz (:,1:N_Seg_Lines) = R2d_Buffer
         case(7)
            Geo % Solzen (:,1:N_Seg_Lines) = R2d_Buffer
         end select

         if (I_Geo /= 3) deallocate ( R2d_Buffer)
         if (I_Geo == 3) deallocate ( R1d_Buffer)

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


      ! --- read m-band data
      ! find m-band file
      if (Sensor%Platform_Name == 'SNPP') then
        if (index(trim(VGEOM_File),'VGEOM_snpp') > 0) then  
           call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                      'VL1BM',File_2_Read)   
        else
           call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                      'VNP02MOD',File_2_Read)
        endif
      endif

      if (Sensor%Platform_Name == 'NOAA-20') then
        call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                    'VJ102MOD',File_2_Read)
      endif


      ! --- if no file quit
      if (trim(File_2_Read) .eq. 'no_file') then
         Error_Out = 1
         return
      endif

      ! --- calculate bowtie pixels
      call COMPUTE_VIIRS_BOWTIE_GAP_PATTERN()

      ! --- read global attribute Day_Night_Flag
      !!!!! Note: if Day_Night_Flag = 'Night' NO M01 - M06, M09, and M11 bands  !!!!!
      !!!!! If Day_Night_Flag = 'Day' or 'Both' then ALL M-bands exist         !!!!!
      call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                                       'DayNightFlag', Day_Night_Flag)
      Mband_Start = 1
      if (trim(Day_Night_Flag) .eq. 'Night') Mband_Start = 7

      ! --- mapping modis to viirs
      !                 041 044 048 055 068  074  085 124 138  160 225 375  405  855  108  120
      !                 M1  M2  M3  M4  M5   M6   M7  M8  M9   M10 M11 M12  M13  M14  M15  M16
      Modis_Chn_List = [ 8 , 9 , 3 , 4 , 1 , 15 , 2 , 5 , 26 , 6 , 7 , 20 , 22 , 29 , 31 , 32 ]
      Is_Mband_On = Sensor%Chan_On_Flag_Default (Modis_Chn_List) == sym%YES
      N_Mband = 16

      ! --- loop over m-band channels
      do I_Mband = Mband_Start , N_Mband

         ! - one more filter for missing night bands
         if (trim(Day_Night_Flag) .eq. 'Night' .and. &
                (I_Mband .eq. 9 .or. I_Mband .eq. 11)) cycle

         ! - check if channel is on 
         Is_Mband_Read(I_Mband) = .false.
         if (.not. Is_Mband_On(I_Mband)) cycle
          
         ! - make string band number
         write (Band_Num_Str, '(i0.2)' )  I_Mband
         Setname_Name = 'observation_data/M'//Band_Num_Str

         ! - read data
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                             Setname_Name, Offset_Mband, Dim_Seg, I2d_Buffer)

         ! - read attribute to unscale
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/scale_factor', Scale_Factor)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/add_offset', Add_Offset)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/_FillValue', Fill_Value)

         ! - separate reflectance and radiance data
         if (I_Mband .le. 11) then
            ! - initialize output
            Ch (Modis_Chn_List(I_Mband)) % Ref_Toa (:, 1:N_Seg_Lines) = Missing_Value_Real4

            ! - unscale
            where ( I2d_Buffer .ne. Fill_Value)
               Ch (Modis_Chn_List(I_Mband)) % Ref_Toa (:, 1:N_Seg_Lines) = &
                               ((I2d_Buffer * Scale_Factor) + Add_Offset) * 100.0
            endwhere

            ! - fill the gaps
            call FILL_VIIRS_BOWTIE_GAPS ( Ny_Start, N_Seg_Lines, Ch (Modis_Chn_List(I_Mband)) % Ref_Toa )

         else

            ! - initialize output
            Ch (Modis_Chn_List(I_Mband)) % Rad_Toa (:, 1:N_Seg_Lines) = Missing_Value_Real4

            ! - unscale
            where ( I2d_Buffer .ne. Fill_Value)
               Ch (Modis_Chn_List(I_Mband)) % Rad_Toa (:, 1:N_Seg_Lines) = &
                                 ((I2d_Buffer * Scale_Factor) + Add_Offset)
            endwhere

            ! - fill the gaps
            call FILL_VIIRS_BOWTIE_GAPS ( Ny_Start, N_Seg_Lines, Ch (Modis_Chn_List(I_Mband)) % Rad_Toa )

            ! - convert radiance to noaa format
            call CONVERT_VIIRS_NASA_RADIANCE (Ch(Modis_Chn_List(I_Mband))%Rad_Toa(:, 1:N_Seg_Lines), &
                                Planck_Nu(Modis_Chn_List(I_Mband)), Missing_Value_Real4) 

            ! - calculate bt
            call COMPUTE_BT_ARRAY (Ch(Modis_Chn_List(I_Mband))%Bt_Toa (:,1:N_Seg_Lines), &
                                   Ch(Modis_Chn_List(I_Mband))%Rad_Toa (:,1:N_Seg_Lines), &
                                   Modis_Chn_List(I_Mband), Missing_Value_Real4)

         endif

         deallocate (I2d_Buffer)
         Is_Mband_Read(I_Mband) = .true.

      enddo  ! - m-bands loop


      ! --- read I-bands

      Read_Iband_Flag = .false.
      Modis_I_Chn_List = [ 39, 40, 41, 42, 43 ]
      N_Iband = 5
      Is_Iband_On = Sensor%Chan_On_Flag_Default (Modis_I_Chn_List) == sym%YES

      ! - check if at least 1 channel is on
      Read_Iband_Flag = ANY(Is_Iband_On)

      if (Read_Iband_Flag) then
        ! find I-band file
        if (Sensor%Platform_Name == 'SNPP') then
          if (index(trim(VGEOM_File),'VGEOM_snpp') > 0) then  
             call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                        'VL1BI',File_2_Read)   
          else
             call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                      'VNP02IMG',File_2_Read)
          endif
        endif

        if (Sensor%Platform_Name == 'NOAA-20') then
          call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                      'VJ102IMG',File_2_Read)
        endif

        ! - set offset and dims to read (Iband = 2 * Mband)
        Ny_Start_Iband = 2 * Ny_Start - 1
        Offset_Iband = [ 0, Ny_Start_Iband - 1 ]
        Dim_Seg_Iband = 2 * Dim_Seg

        ! --- if file found read
        if (trim(File_2_Read) .ne. 'no_file') then

           ! - loop over channels
           do I_Iband = 1, N_Iband

              ! - check if channel is on 
              Is_Iband_Read(I_Iband) = .false.
              if (.not. Is_Iband_On(I_Iband)) cycle

              ! - one more filter for missing night bands
              if (trim(Day_Night_Flag) .eq. 'Night' .and. I_Iband .le. 3) cycle

              ! - make string band number
              write (Band_Num_Str, '(i0.2)' )  I_Iband
              Setname_Name = 'observation_data/I'//Band_Num_Str

              ! - read data
              call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                               Setname_Name, Offset_Iband, Dim_Seg_Iband, I2d_Buffer)
 
              ! - read attribute to unscale
              call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                               Setname_Name//'/scale_factor', Scale_Factor)
              call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                               Setname_Name//'/add_offset', Add_Offset)
              call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                               Setname_Name//'/_FillValue', Fill_Value)

              ! - initialize temp buffer
              allocate (Iband_Temp_Buff(Dim_Seg_Iband(1),Dim_Seg_Iband(2)))
              allocate (Iband_Temp_Buff_Bt(Dim_Seg_Iband(1),Dim_Seg_Iband(2)))
              Iband_Temp_Buff = Missing_Value_Real4
              Iband_Temp_Buff_Bt = Missing_Value_Real4
              Temp_Pix_Array_1 = Missing_Value_Real4

              ! - separate reflectance and radiance data
              if (I_Iband .le. 3) then
                 ! - unscale
                 where (I2d_Buffer .ne. Fill_Value)
                    Iband_Temp_Buff = ((I2d_Buffer * Scale_Factor) + Add_Offset) * 100.0
                 endwhere

                 ! - remove unrealistic
                 where (Iband_Temp_Buff > 120.0)
                    Iband_Temp_Buff = Missing_Value_Real4
                 endwhere

                 ! - calculate min, max, mean
                 select case(I_Iband)
                 case(1) ! 0.65um
                    if (ch(1)%Sub_Pixel_On_Flag) then
                       call COMPUTE_IBAND_STATISTICS(Iband_Temp_Buff, &
                                                Ch(1)%Ref_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                Ch(1)%Ref_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                Temp_Pix_Array_1(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                X_Sample_Offset, Y_Sample_Offset)
                       call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                Ch(1)%Ref_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                       call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                Ch(1)%Ref_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                       if (Use_Iband) then 
                          ch(1)%Ref_Toa(1:Dim_Seg(1),1:Dim_Seg(2)) = Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2))
                          call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                Ch(1)%Ref_Toa(1:Dim_Seg(1),1:Dim_Seg(2)))
                       endif
                    endif

                 case(2) ! 0.86um
                    if (ch(2)%Sub_Pixel_On_Flag) then
                       call COMPUTE_IBAND_STATISTICS(Iband_Temp_Buff, &
                                                   Ch(2)%Ref_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                   Ch(2)%Ref_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                   Temp_Pix_Array_1(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                   Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                   X_Sample_Offset, Y_Sample_Offset)
                       call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                   Ch(2)%Ref_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                       call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                   Ch(2)%Ref_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                       if (Use_Iband) then 
                          ch(2)%Ref_Toa(1:Dim_Seg(1),1:Dim_Seg(2)) = Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2))
                          call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                Ch(2)%Ref_Toa(1:Dim_Seg(1),1:Dim_Seg(2)))
                       endif
                    endif
                 case(3) ! 1.61um
                    if (ch(6)%Sub_Pixel_On_Flag) then
                       call COMPUTE_IBAND_STATISTICS(Iband_Temp_Buff, &
                                                   Ch(6)%Ref_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                   Ch(6)%Ref_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                   Temp_Pix_Array_1(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                   Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                   X_Sample_Offset, Y_Sample_Offset)
                       call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                   Ch(6)%Ref_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                       call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                   Ch(6)%Ref_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                       if (Use_Iband) then 
                          ch(6)%Ref_Toa(1:Dim_Seg(1),1:Dim_Seg(2)) = Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2))
                          call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                Ch(6)%Ref_Toa(1:Dim_Seg(1),1:Dim_Seg(2)))
                       endif
                    endif
                 end select

              else

                 ! - unscale
                 where (I2d_Buffer .ne. Fill_Value)
                    Iband_Temp_Buff = ((I2d_Buffer * Scale_Factor) + Add_Offset)
                 endwhere

                 ! - convert radiance to noaa format
                 call CONVERT_VIIRS_NASA_RADIANCE (Iband_Temp_Buff, &
                                Planck_Nu(Modis_I_Chn_List(I_Iband)), Missing_Value_Real4) 

                 ! - calculate bt
                 call COMPUTE_BT_ARRAY (Iband_Temp_Buff_Bt, Iband_Temp_Buff, &
                                   Modis_I_Chn_List(I_Iband), Missing_Value_Real4)

                 ! - remove unrealistic
                 where (Iband_Temp_Buff_Bt > 360.0)
                    Iband_Temp_Buff_Bt = Missing_Value_Real4
                 endwhere

                 ! - calculate min, max, mean
                 select case(I_Iband)
                 case(4) ! 3.75um
                    if (ch(20)%Sub_Pixel_On_Flag) then
                      call COMPUTE_IBAND_STATISTICS(Iband_Temp_Buff_Bt, &
                                                    Ch(20)%Bt_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                    Ch(20)%Bt_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                    Ch(20)%Bt_Toa_Mean_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                    Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                    X_Sample_Offset, Y_Sample_Offset)
                   
                      call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                  Ch(20)%Bt_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                      call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                  Ch(20)%Bt_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                      call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                  Ch(20)%Bt_Toa_Mean_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                      call COMPUTE_RAD_ARRAY (Ch(20)%Bt_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                              Ch(20)%Rad_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                              20, Missing_Value_Real4)
                      call COMPUTE_RAD_ARRAY (Ch(20)%Bt_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                              Ch(20)%Rad_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                              20, Missing_Value_Real4)
                      if (Use_Iband) then 
                         ch(20)%Bt_Toa(1:Dim_Seg(1),1:Dim_Seg(2)) = Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2))
                         call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                     Ch(20)%Bt_Toa(1:Dim_Seg(1),1:Dim_Seg(2)))
                         ch(20)%Rad_Toa(1:Dim_Seg(1),1:Dim_Seg(2)) = Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2))
                         call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                     Ch(20)%Rad_Toa(1:Dim_Seg(1),1:Dim_Seg(2)))
                      endif
                     endif
                 case(5) ! 11um
                    if (ch(31)%Sub_Pixel_On_Flag) then
                      call COMPUTE_IBAND_STATISTICS(Iband_Temp_Buff_Bt, &
                                                Ch(31)%Bt_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                Ch(31)%Bt_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                Ch(31)%Bt_Toa_Mean_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                                X_Sample_Offset, Y_Sample_Offset)
                      call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                Ch(31)%Bt_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                      call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                Ch(31)%Bt_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                      call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                Ch(31)%Bt_Toa_Mean_Sub(1:Dim_Seg(1),1:Dim_Seg(2)))
                      call COMPUTE_RAD_ARRAY (Ch(31)%Bt_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                           Ch(31)%Rad_Toa_Min_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                           31, Missing_Value_Real4)
                      call COMPUTE_RAD_ARRAY (Ch(31)%Bt_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                           Ch(31)%Rad_Toa_Max_Sub(1:Dim_Seg(1),1:Dim_Seg(2)), &
                                           31, Missing_Value_Real4)
                      if (Use_Iband) then 
                          ch(31)%Bt_Toa(1:Dim_Seg(1),1:Dim_Seg(2)) = Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2))
                          call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                Ch(31)%Bt_Toa(1:Dim_Seg(1),1:Dim_Seg(2)))
                          ch(31)%Rad_Toa(1:Dim_Seg(1),1:Dim_Seg(2)) = Temp_Pix_Array_2(1:Dim_Seg(1),1:Dim_Seg(2))
                          call FILL_VIIRS_BOWTIE_GAPS(Ny_Start,N_Seg_Lines, &
                                                Ch(31)%Rad_Toa(1:Dim_Seg(1),1:Dim_Seg(2)))
                      endif
                    endif
                 end select

               endif

               ! - dealocate all
               deallocate (I2d_Buffer)
               if (allocated(Iband_Temp_Buff)) deallocate (Iband_Temp_Buff)
               if (allocated(Iband_Temp_Buff_Bt)) deallocate (Iband_Temp_Buff_Bt)
               Is_Iband_Read(I_Iband) = .true.

            enddo ! - I-bands loop

         else ! I-band file not found
            if (Segment_Number == 1) &
              call MESG ( "I-Band File NOT Found, Skipping", level = verb_lev % DEFAULT)
         endif
      endif ! Read_Iband_Flag


      ! --- read fusion channels
      if (Sensor%Fusion_Flag) then

        ! - find fusion file
        if (Sensor%Platform_Name == 'SNPP') then
          call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path), &
                           trim(Image%Level1b_Name), 'VNP02FSN',File_2_Read)
          if (trim(File_2_Read) .eq. 'no_file') then
               call DETERMINE_VIIRS_NASA_FILE( trim(Image%Level1b_Path), &
                    trim(Image%Level1b_Name), 'FSNRAD_L2_VIIRS_CRIS_SNPP',File_2_Read)
          endif
        elseif (Sensor%Platform_Name == 'NOAA-20') then
          call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path), &
                           trim(Image%Level1b_Name), 'VJ102FSN',File_2_Read)
          if (trim(File_2_Read) .eq. 'no_file') then
               call DETERMINE_VIIRS_NASA_FILE( trim(Image%Level1b_Path), &
                    trim(Image%Level1b_Name), 'FSNRAD_L2_VIIRS_CRIS_NOAA20',File_2_Read)
          endif
        endif

        ! - read from ch 23 to 36 = total 14 ch.
        Fusion_Start = 23
        N_Fusion = 36
 
        ! - loop over all fusion channels
        do I_Fusion = Fusion_Start, N_Fusion

           ! - check if channel is on 
           if (Sensor%Chan_On_Flag_Default (I_Fusion) == sym%NO) cycle

           ! - cycle if there is viirs channel 
           if (I_Fusion == 26 .or. &
               I_Fusion == 29) cycle

           ! - make string band number
           write (Fusion_Band_Str, '(i0.2)' ) I_Fusion
           Setname_Name = 'observation_data/MODIS'//Fusion_Band_Str
           if (index(trim(File_2_Read),'FSNRAD') > 0) then
              Setname_Name = 'geophysical_data/MODIS'//Fusion_Band_Str
           endif

           ! - read data
           call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                               Setname_Name, Offset_Mband, Dim_Seg, I2d_Buffer)

           ! - read attribute to unscale
           call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                               Setname_Name//'/scale_factor', Scale_Factor)
           call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                               Setname_Name//'/add_offset', Add_Offset)
           call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                               Setname_Name//'/_FillValue', Fill_Value)

           ! - common case not ch 31 or 32
           if (I_Fusion /= 31 .and. I_Fusion /= 32) then
             
             ! - initialize output
             Ch (I_Fusion) % Rad_Toa (:, 1:N_Seg_Lines) = Missing_Value_Real4

             ! - unscale
             where (I2d_Buffer .ne. Fill_Value)
                Ch (I_Fusion) % Rad_Toa (:, 1:N_Seg_Lines) = &
                                ((I2d_Buffer * Scale_Factor) + Add_Offset)
             endwhere

             ! - fill the gaps
             call FILL_VIIRS_BOWTIE_GAPS ( Ny_Start, N_Seg_Lines, Ch (I_Fusion) % Rad_Toa )

             ! - convert radiance to noaa format
             call CONVERT_VIIRS_NASA_RADIANCE (Ch(I_Fusion)%Rad_Toa(:, 1:N_Seg_Lines), &
                                Planck_Nu(I_Fusion), Missing_Value_Real4)

             ! - calculate bt
             call COMPUTE_BT_ARRAY (Ch(I_Fusion)%Bt_Toa (:,1:N_Seg_Lines), &
                                    Ch(I_Fusion)%Rad_Toa(:,1:N_Seg_Lines), &
                                    I_Fusion, Missing_Value_Real4)

           ! - special case ch 31 or 32
           else

             ! - allocate temporary array
             allocate (R2d_Temp(Dim_Seg(1),Dim_Seg(2)))
             R2d_Temp = Missing_Value_Real4

             ! - unscale
             where (I2d_Buffer .ne. Fill_Value)
                R2d_Temp = ((I2d_Buffer * Scale_Factor) + Add_Offset)
             endwhere

             ! - fill the gaps
             call FILL_VIIRS_BOWTIE_GAPS ( Ny_Start, N_Seg_Lines, R2d_Temp )

             ! - 11um
             if (I_Fusion == 31) then

               ! - initialize output
               Bt_11um_Sounder = Missing_Value_Real4
             
               ! - convert radiance to noaa format
               call CONVERT_VIIRS_NASA_RADIANCE (R2d_Temp, &
                                Planck_Nu_11um_Sndr, Missing_Value_Real4)

               ! - calculate bt
               call COMPUTE_BT_ARRAY(Bt_11um_Sounder(:,1:N_Seg_Lines), &
                                      R2d_Temp, I_Fusion, Missing_Value_Real4)

             endif

             ! - 12um
             if (I_Fusion == 32) then

               ! - initialize output
               Bt_12um_Sounder = Missing_Value_Real4

               ! - convert radiance to noaa format
               call CONVERT_VIIRS_NASA_RADIANCE (R2d_Temp, &
                                Planck_Nu_12um_Sndr, Missing_Value_Real4)

               ! - calculate bt
               call COMPUTE_BT_ARRAY (Bt_12um_Sounder(:,1:N_Seg_Lines), &
                                      R2d_Temp, I_Fusion, Missing_Value_Real4)

             endif

             ! - deallocate temp array
             deallocate (R2d_Temp)

           endif ! special case

           deallocate (I2d_Buffer)

        enddo ! loop channels
      endif ! fusion

      ! --- read dnb data (use do loop to quit if no file)
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then

      do k = 1, 1

         ! - find dnb geo file
         if (Sensor%Platform_Name == 'SNPP') then
           if (index(trim(VGEOM_File),'VGEOM_snpp') > 0) then
              call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                    'VGEOD',File_2_Read)
           else
              call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                    'VNP03DNB',File_2_Read)
           endif
         endif

         if (Sensor%Platform_Name == 'NOAA-20') then
           call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                    'VJ103DNB',File_2_Read)
         endif

         ! - if no file quit
         if (trim(File_2_Read) .eq. 'no_file') then
            !print *,'Error: No VGEOD File Found, Skipping'
            call MESG(" Error: No VNP03DNB File Found, Skipping", level = verb_lev % DEFAULT)
            cycle
         endif

         ! - mapping file (maps from dnb to M-bands resolution)
         File_Dnb_Idx = trim(Ancil_Data_Dir)//'static/viirs/dnb2m_indx.txt'
         Lun = GETLUN()
         Dim_Dnb_Seg(1) = 4064
         Dim_Dnb_Seg(2) = Dim_Seg(2)

         open (unit = Lun , file=trim ( File_Dnb_Idx) , status="old",action="read")
         read (unit = Lun , fmt=* , iostat = Io_Err_Stat) D2m_Idx
         close (unit = Lun)

         ! - read dnb geo data
         Setname_Name = 'geolocation_data/lunar_azimuth'
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                        Setname_Name, Offset_Mband, Dim_Dnb_Seg, I2d_Buffer)

         ! - read attribute to unscale
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/scale_factor', Scale_Factor)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/add_offset', Add_Offset)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/_FillValue', Fill_Value)

         ! - unscale and remap
         do i = 1, Dim_Seg(1)
            Geo % Lunaz (i, 1:N_Seg_Lines) = (I2d_Buffer(D2m_Idx(i), :) &
                                    * Scale_Factor) + Add_Offset
            where ( i2d_buffer(d2m_idx(i),:) .eq. Fill_Value)
               geo % lunaz(i,1:n_seg_lines) = MISSING_VALUE_REAL4
            end where                        
         end do
         deallocate ( I2d_Buffer )

         Setname_Name = 'geolocation_data/lunar_zenith'
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                        Setname_Name, Offset_Mband, Dim_Dnb_Seg, I2d_Buffer)

         ! - read attribute to unscale
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/scale_factor', Scale_Factor)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/add_offset', Add_Offset)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/_FillValue', Fill_Value)

         ! - unscale and remap
         do i = 1, Dim_Seg(1)
            Geo % Lunzen (i, 1:N_Seg_Lines) = (I2d_Buffer(D2m_Idx(i), :) &
                                    * Scale_Factor) + Add_Offset
            where( i2d_buffer(d2m_idx(i),:) .eq. Fill_Value)
               geo % lunzen (i,1:n_seg_lines) = MISSING_VALUE_REAL4
            end where                        
         end do

      

         ! --- read moon phase angle 
!         call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(File_2_Read), &
!                      'moon_phase_angle', Geo % Moon_Phase_Angle)
!         call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(File_2_Read), &                                                                                   
!                      'moon_illumination_fraction', Geo % Moon_Illum_Frac)
         Setname_Name = 'geolocation_data/moon_phase_angle'
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name,Offset_Mband, Dim_Dnb_Seg, I2d_Buffer)

         ! - read attribute to unscale
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/scale_factor', Scale_Factor)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/add_offset', Add_Offset)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/_FillValue', Fill_Value)

         ! - unscale and save mean value
         N_Valid = count (I2d_Buffer .ne. Fill_Value)
         Geo % Moon_Phase_Angle = sum((I2d_Buffer * Scale_Factor) + Add_Offset,MASK=I2d_Buffer /= Fill_Value) / N_Valid
         deallocate ( I2d_Buffer )

         ! --- read moon illumination fraction
         Setname_Name = 'geolocation_data/moon_illumination_fraction'
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name,Offset_Mband, Dim_Dnb_Seg, I2d_Buffer)

         ! - read attribute to unscale
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/scale_factor', Scale_Factor)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/add_offset', Add_Offset)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/_FillValue', Fill_Value)

         ! - unscale and save mean value
         N_Valid = count (I2d_Buffer .ne. Fill_Value)
         Geo % Moon_Illum_Frac = sum((I2d_Buffer * Scale_Factor) + Add_Offset,MASK=I2d_Buffer /= Fill_Value) / N_Valid
         deallocate ( I2d_Buffer )
         Geo % LunFrac = Geo % Moon_Illum_Frac

         ! --- compute lunar relative azimuth
         Geo % Lunrelaz (:, 1:N_Seg_Lines) = RELATIVE_AZIMUTH ( &
                   Geo % Lunaz (:, 1:N_Seg_Lines), &
                   Geo % Sataz (:, 1:N_Seg_Lines) )

         ! --- compute lunar glint zenith
         Geo % Glintzen_Lunar (:, 1:N_Seg_Lines) = GLINT_ANGLE ( &
                   Geo % Lunzen (:, 1:N_Seg_Lines), &
                   Geo % Satzen (:, 1:N_Seg_Lines), &
                   Geo % Lunrelaz (:, 1:N_Seg_Lines) )

         ! --- compute lunar scattering angle
         Geo % Scatangle_Lunar (:, 1:N_Seg_Lines) = SCATTERING_ANGLE ( &
                   Geo % Lunzen (:, 1:N_Seg_Lines), &
                   Geo % Satzen (:, 1:N_Seg_Lines), &
                   Geo % Lunrelaz (:, 1:N_Seg_Lines) )

         ! --- find dnb file
         if (Sensor%Platform_Name == 'SNPP') then
           if (index(trim(VGEOM_File),'VGEOM_snpp') > 0) then
              call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                      'VL1BD',File_2_Read)
           else
              call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                      'VNP02DNB',File_2_Read)
           endif
         endif

         if (Sensor%Platform_Name == 'NOAA-20') then
           call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                    'VJ102DNB',File_2_Read)
         endif

         ! - if no file quit
         if (trim(File_2_Read) .eq. 'no_file') then
            !print *,'Error: No VL1BD File Found, Skipping'
            call MESG(" Error: No VL1BD File Found, Skipping", level = verb_lev % DEFAULT)
            cycle
         endif

         ! --- read dnb data
         Setname_Name = 'observation_data/DNB_observations'
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                        Setname_Name, Offset_Mband, Dim_Dnb_Seg, R2d_Buffer)

         ! - remap
         do i = 1, Dim_Seg(1)
            Ch (44) % Rad_Toa (i, 1:N_Seg_Lines) = R2d_Buffer(D2m_Idx(i), :)
         end do
         deallocate ( R2d_Buffer )

         ! - convert radiance to reflectance
         call CONVERT_RAD_2_SOL_REF_DNB (Ch (44) % Rad_Toa, Geo % Solzen, &
                      Sun_Earth_Distance, Missing_Value_Real4, Ch(44)%Ref_Toa)

      enddo ! one time loop

      endif ! dnb on


      ! --- global variables which have to be set
      Image%Number_Of_Lines_Read_This_Segment = N_Seg_Lines
      do i = 1, Image%Number_Of_Lines_Per_Segment
         Image%Scan_Number(i) = Ny_Start + i - 1
      end do

end subroutine READ_VIIRS_NASA_DATA

!--------------------------------------------------------------------
! read viirs nasa time
!  called in sensor_mod
!--------------------------------------------------------------------
subroutine READ_VIIRS_NASA_DATE_TIME (Path, Infile, Year , Doy , Start_Time &
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
           'OrbitNumber', Orbit)

      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'time_coverage_start', Tmp_String)
      !0        1         2
      !123456789012345678901234567890
      !2016-04-20T12:00:00.000Z
      read(Tmp_String(1:4), fmt="(I4)") Year
      read(Tmp_String(6:7), fmt="(I2)") Month
      read(Tmp_String(9:10), fmt="(I2)") Day
      read(Tmp_String(12:13), fmt="(I2)") Start_Hour
      read(Tmp_String(15:16), fmt="(I2)") Start_Minute
      read(Tmp_String(18:19), fmt="(I2)") Start_Sec

      !--- compute start day of year
      call JULIAN ( Day, Month, Year, Doy )
      
      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'time_coverage_end', Tmp_String)
      read(Tmp_String(1:4), fmt="(I4)") End_Year
      read(Tmp_String(6:7), fmt="(I2)") Month
      read(Tmp_String(9:10), fmt="(I2)") Day
      read(Tmp_String(12:13), fmt="(I2)") End_Hour
      read(Tmp_String(15:16), fmt="(I2)") End_Minute
      read(Tmp_String(18:19), fmt="(I2)") End_Sec

      !--- compute end day of year
      call JULIAN ( Day, Month, End_Year, End_Doy )

      ! --- Calculate start and end time
      Start_Time = ((Start_Hour * 60 + Start_Minute) * 60 + Start_Sec) * 1000
      End_Time = ((End_Hour * 60 + End_Minute) * 60 + End_Sec) * 1000

end subroutine READ_VIIRS_NASA_DATE_TIME

!--------------------------------------------------------------------
! read viirs nasa number of scans
!  called in sensor_mod.f90
!--------------------------------------------------------------------
subroutine READ_NUMBER_OF_SCANS_VIIRS_NASA (Infile, Number_Of_Viirs_Lines, &
                 Error_Out)

      character(len=*), intent(in) :: Infile
      integer(kind=int4), intent(out) :: Error_Out
      integer(kind=int4), intent(out) :: Number_of_Viirs_Lines

      character(len=100) :: Setname
      integer, dimension(:), pointer :: Test
      integer, dimension(1) :: Dims
     


      Error_Out = 0
      Setname = 'scan_line_attributes/scan_quality'

      call H5ReadDataset( trim(Infile), trim(Setname), Test )
      Dims = shape(Test)

      ! --- error 
      if (Dims(1) .lt. 1) then
         Number_Of_Viirs_Lines = -999
         Error_Out = 1
         Return
      endif
         
      Number_Of_Viirs_Lines = Dims(1) * 16

end subroutine READ_NUMBER_OF_SCANS_VIIRS_NASA

!--------------------------------------------------------------------
! find viirs nasa files
!--------------------------------------------------------------------
subroutine DETERMINE_VIIRS_NASA_FILE(Path_In,File_In,File_Type_In,File_Out)

      character(len=*), intent(in):: Path_In
      character(len=*), intent(in):: File_In
      character(len=*), intent(in):: File_Type_In
      character(len=*), intent(out):: File_Out

      character(len=500):: Search_String
      character(len=1020), dimension(:), pointer:: Files
      integer(kind=int4):: Num_Files


      !0        1         2         3         4         5
      !12345678901234567890123456789012345678901234567890
      !VGEOM_snpp_d20160420_t120000_c20160420175142.nc
      !VL1BM_snpp_d20160420_t120000_c20160420172929.nc
      !VGEOD_snpp_d20140205_t115400_c20170405013130.nc
      !VL1BD_snpp_d20140205_t115400_c20170405012902.nc

      if (index(File_Type_In, 'VL1BM') > 0 .or. &
          index(File_Type_In, 'VGEOD') > 0 .or. &
          index(File_Type_In, 'VL1BD') > 0) then 
         Search_String = trim(File_Type_In)//trim(File_In(6:28))//'*.nc'
      else
      !0        1         2         3         4         5
      !12345678901234567890123456789012345678901234567890
      !VNP03MOD.A2014036.1154.001.2017300063048.uwssec.nc
      !VNP02MOD.A2014036.1154.001.2017300062642.uwssec.nc
      !VNP03DNB.A2014036.1154.001.2017300063048.uwssec.nc
      !VNP02DNB.A2014036.1154.001.2017300062642.uwssec.nc
         Search_String = trim(File_Type_In)//trim(File_In(9:22))//'*.nc'
      endif

      Files => FILE_SEARCH(trim(Path_In),trim(Search_String),count=Num_Files)

      if (Num_Files == 0) then
         !print *, EXE_PROMPT, VIIRS_NASA_PROMPT, " No "//trim(File_Type_In)//" File Found"
         !call MESG(" No "//trim(File_Type_In)//" File Found", level = verb_lev % DEFAULT)
         File_Out = "no_file"
         return
      endif

      if (Num_Files > 1) then
         !print *, EXE_PROMPT, VIIRS_NASA_PROMPT, "Multiple "//trim(File_Type_In)//" Files Found"
         !call MESG(" Multiple "//trim(File_Type_In)//" Files Found", level = verb_lev % DEFAULT)
!         File_Out = "no_file"
      endif

      File_Out = Files(1)
!      print *, EXE_PROMPT, VIIRS_NASA_PROMPT, trim(File_Type_In)//" File Found, ",trim(File_Out)

      Files => null()

end subroutine DETERMINE_VIIRS_NASA_FILE

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
! Function Name: CONVERT_VIIRS_RADIANCE
!
! Function:
!    Convert to units of the VIIRS radiance values from the that used
!    in the IDPS level-1b to that expected by CLAVR-x
!
! Description: 
!   
! Calling Sequence: rad_new =
! convert_viirs_radiance(rad_old,nu,missing_value)
!   
!
! Inputs:
!   rad_old = radiance in units of W/m^2/micron/str (2d array)
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
subroutine CONVERT_VIIRS_NASA_RADIANCE(Radiance,Nu,Missing_Value)
      real (kind=real4), dimension(:,:), intent(inout):: Radiance
      real (kind=real4), intent(in):: Nu
      real (kind=real4), intent(in):: Missing_Value

      where(Radiance /= Missing_Value)
         Radiance = Radiance * (((10000.0 / Nu )**2) / 10.0)
      end where

      return

end subroutine CONVERT_VIIRS_NASA_RADIANCE

!--------------------------------------------------
subroutine CONVERT_RAD_2_SOL_REF_DNB( Radiance, Solzen ,Sun_Earth_Distance, &
                                      Missing_Value , Reflectance )

      real (kind=real4), dimension(:,:), intent(in) :: Radiance
      real (kind=real4), dimension(:,:), intent(in) :: Solzen
      real (kind=real4), intent(in) :: Missing_Value
      real (kind=real4), intent(in) :: Sun_Earth_Distance
      real (kind=real4), dimension(:,:), intent(out) :: Reflectance
      real (kind=real4), parameter :: Fo = 0.044217282   ! dnb solar energy in w/cm2 ( Source? )
      real , parameter :: PI = 3.14159265359
      real , parameter :: DTOR = PI / 180.

      Reflectance = Missing_Value

      where(Radiance /= Missing_Value)
         Reflectance = 100. * (Pi * Radiance * Sun_Earth_Distance ** 2) / &
                            (cos(Solzen * DTOR) * Fo)
      end where

end subroutine CONVERT_RAD_2_SOL_REF_DNB

!-------------------------------------------------------------------------------------
! subroutine to compute the bowtie gap pattern.  
! This pattern repeats every 48 scans 
! all VIIRS files should be integer multiples of this pattern
!
! Gap_Pixel_Mask_Pattern = a binary mask that identifies these bowtie gaps 
! Gap_Line_Idx = line index for each pixel in pattern including gap pixels
!
! A description of the pattern
!
!--------  line type 3
!----      line type 4
! 12 lines without gaps
!----      line type 1
!--------  line type 2
!--------  line type 3
!----      line type 4
! 12 lines without gaps
!----      line type 1
!--------  line type 2
!--------  line type 3
!----      line type 4
! 12 lines without gaps
!----      line type 1
!--------  line type 2
!
!  line types 1 and 4 have 640 missing pixels
!  line types 2 and 3 have 1008 missing pixels
! 
!-------------------------------------------------------------------------------------
   subroutine COMPUTE_VIIRS_BOWTIE_GAP_PATTERN()

      integer (kind=int4), dimension(Ny_Pattern):: Line_Type
      integer (kind=int4), parameter:: Ngap_1 = 640  !1280
      integer (kind=int4), parameter:: Ngap_2 = 1008 !2016
      integer (kind=int4), parameter:: Ngap_3 = 1008 !2016
      integer (kind=int4), parameter:: Ngap_4 = 640  !1280
      integer (kind=int4):: Iline
      integer (kind=int4):: i1
      integer (kind=int4):: i2

      !--- define the line patterns as described above
      Line_Type = &
      (/3,4,0,0,0,0,0,0,0,0,0,0,0,0,1,2/)


      if (.not. allocated(Gap_Line_Idx_Pattern)) &
           allocate(Gap_Line_Idx_Pattern(Nx_Pattern,Ny_Pattern))
      if (.not. allocated(Gap_Pixel_Mask_Pattern)) &
           allocate(Gap_Pixel_Mask_Pattern(Nx_Pattern,Ny_Pattern))

      do Iline = 1 , Ny_Pattern

         Gap_Line_Idx_Pattern( : , Iline ) = -999
         Gap_Pixel_Mask_Pattern( : , Iline ) = 0

         if (Line_Type(Iline) == 1) then

            i1 = 1
            i2 = Ngap_1
            Gap_Line_Idx_Pattern( i1 : i2 , Iline) = Iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline) = 1

            i1 = Nx_Pattern - Ngap_1 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline ) = 1
         end if

         if (Line_Type(Iline) == 2) then
            i1 = 1
            i2 = Ngap_1
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline - 2
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline ) = 1

            i1 = Ngap_1 + 1
            i2 = Ngap_2
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline ) = 1

            i1 = Nx_Pattern - Ngap_1 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline - 2
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline ) = 1

            i1 = Nx_Pattern - Ngap_2 + 1
            i2 = i1  + (Ngap_2 - Ngap_1)
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline ) = 1
         end if

         if (Line_Type(Iline) == 3) then
            i1 = 1
            i2 = Ngap_4
            Gap_Line_Idx_Pattern(i1:i2,Iline) = Iline + 2
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1

            i1 = Ngap_4 + 1
            i2 = Ngap_3
            Gap_Line_Idx_Pattern(i1:i2,Iline) = Iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1

            i1 = Nx_Pattern - Ngap_4 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern(i1:i2,Iline) = Iline + 2
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1

            i1 = Nx_Pattern - Ngap_3 + 1
            i2 = i1 + (Ngap_3 - Ngap_4)
            Gap_Line_Idx_Pattern(i1:i2,Iline) = Iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1
         end if

         if (Line_Type(Iline) ==  4) then
            i1 = 1
            i2 = Ngap_4
            Gap_Line_Idx_Pattern(i1:i2,Iline) = Iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1

            i1 = Nx_Pattern - Ngap_4 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1
         end if

      end do

   end subroutine COMPUTE_VIIRS_BOWTIE_GAP_PATTERN

!------------------------------------------------------------------------------
! this routine uses the bowtie gap pattern, applies it to an arbitrary 
! segment of data and fills in the observations with the closest valid data
!------------------------------------------------------------------------------
    subroutine FILL_VIIRS_BOWTIE_GAPS ( Line_Start, Number_of_Lines  , Variable )
      integer(kind=int4), intent(in) :: Line_Start
      integer(kind=int4), intent(in) :: Number_of_Lines
      real(kind=real4), dimension(:,:), intent (inout) :: Variable

      integer(kind=int4) :: Line_Offset
      integer(kind=int4) :: Line_In_Pattern
      integer(kind=int4) :: Line_in_Segment
      integer(kind=int4) :: Line_Idx
      integer(kind=int4) :: Elem_Idx
      integer(kind=int4), dimension(:,:), allocatable :: Gap_Line_Idx

      allocate (Gap_Line_Idx(Nx_Pattern, Number_Of_Lines))
      Gap_Line_Idx = Missing_Value_Int4

      do Line_In_Segment = 1,  Number_Of_Lines

         Line_In_Pattern = mod (Line_Start - 1 + Line_In_Segment , Ny_Pattern)
         if (Line_In_Pattern == 0) Line_In_Pattern = Ny_Pattern

         Line_Offset = Line_In_Segment - Line_In_Pattern

         if (Gap_Line_Idx_Pattern(1,Line_In_Pattern) /= Missing_Value_Int4) &
             Gap_Line_Idx(:, Line_In_Segment) = Gap_Line_Idx_Pattern(:,Line_In_Pattern) + Line_Offset

         where ( Gap_Line_Idx(:,Line_in_Segment) <= 0 )
            Gap_Line_Idx(:,Line_in_Segment) = 1
         end where

         where ( Gap_Line_Idx(:,Line_in_Segment) > Number_of_Lines )
            Gap_Line_Idx(:,Line_in_Segment) = Number_of_Lines
         end where

         ! - write gap mask to output
         Gap_Pixel_Mask(:,Line_in_Segment) = Gap_Pixel_Mask_Pattern(:,Line_in_Pattern)

         ! - fill gaps
         do Elem_Idx = 1, Nx_Pattern
            if (Gap_Pixel_Mask(Elem_Idx,Line_in_Segment) == 1) then
               Line_Idx = Gap_Line_Idx(Elem_Idx,Line_in_Segment)

               Variable(Elem_Idx, Line_In_Segment) = Variable(Elem_Idx,Line_Idx)

            end if
         end do
      end do

      deallocate (Gap_Line_Idx)

   end subroutine FILL_VIIRS_BOWTIE_GAPS

!----------------------------------------------------------------------------------------
! this code is to check if fusion file exist, if yes sensor is changed to VIIRS-FUSION
! called from sensor_mod.f90
!----------------------------------------------------------------------------------------
   subroutine CHECK_IF_FUSION(Fusion)
     logical, intent (out) :: Fusion
     character(len=1020) :: File_2_Read
     
     ! initialize Fusion = false
     Fusion = .false.

     ! --- search fusion file
     !VJ103MOD.A2018108.1400.001.2018108193630.uwssec.nc
     !VNP03MOD.A2018108.1400.001.2018108202120.uwssec.nc
     !VNP02FSN.A2018109.2106.001.2018110035023.nc
     if (Sensor%Platform_Name == 'SNPP') then
       call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path), &
                          trim(Image%Level1b_Name), 'VNP02FSN',File_2_Read)
       if (trim(File_2_Read) .eq. 'no_file') then
           call DETERMINE_VIIRS_NASA_FILE( trim(Image%Level1b_Path), &
                trim(Image%Level1b_Name), 'FSNRAD_L2_VIIRS_CRIS_SNPP',File_2_Read)
       endif
     elseif (Sensor%Platform_Name == 'NOAA-20') then
       call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path), &
                          trim(Image%Level1b_Name), 'VJ102FSN',File_2_Read)
       if (trim(File_2_Read) .eq. 'no_file') then
           call DETERMINE_VIIRS_NASA_FILE( trim(Image%Level1b_Path), &
                trim(Image%Level1b_Name), 'FSNRAD_L2_VIIRS_CRIS_NOAA20',File_2_Read)
       endif
     else
       call MESG ( "FUSION CHECK IS NOT DESIGNED FOR THIS PLATFORM = "//trim(Sensor%Platform_Name), &
                   level = verb_lev % DEFAULT)
       return
     endif

     ! --- if file found set fusion = true
     if (trim(File_2_Read) .ne. 'no_file') Fusion = .true.
   end subroutine CHECK_IF_FUSION



!----------------------------------------------------------------
! read the VIIRS + Fusion constants into memory
!-----------------------------------------------------------------
   subroutine READ_FUSION_INSTR_CONSTANTS(Instr_Const_File)
      use CALIBRATION_CONSTANTS_MOD
      use FILE_TOOLS , only: GETLUN

      implicit none

      character(len=*), intent(in):: Instr_Const_file
      integer:: Ios0, Erstat
      integer:: Instr_Const_Lun
      character(len=20):: header

      Instr_Const_Lun = GETLUN()

      open(unit=Instr_Const_Lun,file=trim(Instr_Const_File),status="old",position="rewind",action="read",iostat=ios0)
      call MESG ("Opening "// trim(Instr_Const_File),level = verb_lev % DEFAULT)
      Erstat = 0
      if (Ios0 /= 0) then
         Erstat = 19
         call MESG (EXE_PROMPT//"  Error Opening FUSION Constants File",level = verb_lev % DEFAULT)
         stop 19
      endif

      read(unit=Instr_Const_Lun,fmt="(a3)") sat_name
      read(unit=Instr_Const_Lun,fmt=*) Solar_Ch20
      read(unit=Instr_Const_Lun,fmt=*) Ew_Ch20
      !read(unit=Instr_Const_lun,fmt=*) header
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(20), Planck_A2(20), Planck_Nu(20)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(22), Planck_A2(22), Planck_Nu(22)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(23), Planck_A2(23), Planck_Nu(23)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(24), Planck_A2(24), Planck_Nu(24)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(25), Planck_A2(25), Planck_Nu(25)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(27), Planck_A2(27), Planck_Nu(27)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(28), Planck_A2(28), Planck_Nu(28)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(29), Planck_A2(29), Planck_Nu(29)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(30), Planck_A2(30), Planck_Nu(30)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(31), Planck_A2(31), Planck_Nu(31)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(32), Planck_A2(32), Planck_Nu(32)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(42), Planck_A2(42), Planck_Nu(42)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(43), Planck_A2(43), Planck_Nu(43)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(33), Planck_A2(33), Planck_Nu(33)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(34), Planck_A2(34), Planck_Nu(34)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(35), Planck_A2(35), Planck_Nu(35)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1(36), Planck_A2(36), Planck_Nu(36)
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1_11um_Sndr, Planck_A2_11um_Sndr, Planck_Nu_11um_Sndr
      read(unit=Instr_Const_Lun,fmt=*) Planck_A1_12um_Sndr, Planck_A2_12um_Sndr, Planck_Nu_12um_Sndr
      close(unit=Instr_Const_Lun)

      !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
      Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

   end subroutine READ_FUSION_INSTR_CONSTANTS


   !---------------------------------------------------------------------------
   !!!! Taken from "viirs_clavrx_bridge.f90" which was written by Andi W. !!!
   !!!! Only deleted St.Dev.      Denis B.
   !---------------------------------------------------------------------------
   ! - iband has full file dimension of 6400 x1536
   ! - mband 3200 x 768
   !  - output of min_val ... is 3200 768
   !---------------------------------------------------------------------------
   subroutine COMPUTE_IBAND_STATISTICS (Iband_Array, Out_Min_Val, Out_Max_Val , Out_Mean_Val, Out_Sample_Val, X_Offset_In, Y_Offset_In)

      implicit none

      real, dimension(:,:) , intent(in) :: Iband_Array
      real, dimension(:,:) , intent(out)  :: Out_Min_Val, Out_Max_Val , Out_Mean_Val, Out_Sample_Val
      integer, intent(in):: X_Offset_In, Y_Offset_In

      real, dimension(2,2) :: Small_Iband
      integer :: im , jm
      integer , dimension(2) ::  dim_m
      integer :: iband_x0, iband_x1 ,  iband_y0, iband_y1
      integer:: X_Offset, Y_Offset

      !--- to be consistent with geostationary static nav convenction,
      !--- X_Offset_In and Y_Offset_In should 0 or 1

      !--- ensure sampling offsets are between 1 and 2
      X_Offset = X_Offset_In + 1
      Y_Offset = Y_Offset_In + 1
      X_Offset = min(2,max(1,X_Offset))
      Y_Offset = min(2,max(1,Y_Offset))

      dim_m = shape (Out_Min_Val)

      Out_Min_Val = Missing_Value_Real4
      Out_Max_Val = Missing_Value_Real4
      Out_Mean_Val = Missing_Value_Real4

      do im = 1 , dim_m(1)

         iband_x0 = (im-1) * 2 +1
         iband_x1 = iband_x0 + 1
         do jm =1 , dim_m( 2 )
            iband_y0 = ( jm-1) * 2 +1
            iband_y1 = iband_y0 + 1
            Small_Iband = Iband_Array ( iband_x0 :  iband_x1 ,  iband_y0 : iband_y1 )
            if ( minval ( Small_Iband ) > 0 ) then
               Out_Min_Val(im,jm) = minval (Small_Iband)
               Out_Max_Val(im,jm) = maxval (Small_Iband)
               Out_Mean_Val(im,jm) = sum (Small_Iband) / 4.0
               Out_Sample_Val(im,jm) = Small_Iband(X_Offset,Y_Offset)
            end if
         end do
      end do

   end subroutine COMPUTE_IBAND_STATISTICS


!====================================================================


end module VIIRS_NASA_READ_MODULE 


