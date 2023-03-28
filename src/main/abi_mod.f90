!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: abi_module.f90 (src)
!       ABI_MODULE (program)
!
! PURPOSE: This module contains all the subroutines needed to perform navigation and
!          calibration for GOES-16 ABI.
!
! DESCRIPTION: This module assumes band separated images only.  It also assumes
!              that the calibration type is ABIN.
!
!       Navigation is of type ABIN.
!
! AUTHORS:
!  Steve Wanzong, CIMSS, stevew@ssec.wisc.edu
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
! ABI Channel Mapping
!
!  wvl        abi  modis/clavrx
!
!  0.47       1      3
!  0.64       2      1
!  0.87       3      2
!  1.38       4     26
!  1.61       5      6
!  2.25       6      7
!  3.90       7     20
!  6.19       8     37
!  6.95       9     27
!  7.34      10     28
!  8.50      11     29
!  9.61      12     30
!  10.4      13     38
!  11.2      14     31
!  12.3      15     32
!  13.3      16     33
!
!--------------------------------------------------------------------------------------

module ABI_MOD

  use CONSTANTS_MOD,only: &
  int4,real8,int2 &
  , real4 , int1 &
  , MISSING_VALUE_INT4 &
  , MISSING_VALUE_REAL4 &
  , sym

  use PIXEL_COMMON_MOD,only: &
  image &
  , Temporary_Data_Dir &
  , ch &
  , image &
  , sensor &
  , two_byte_temp &
  , geo &
  , nav &
  , l1b_gzip, l1b_bzip2 &
  , Number_of_Temporary_Files &
  , Temporary_File_Name &
  , Ch1_Counts &
  , Line_Idx_Min_Segment

  use univ_fp_comparison_mod, only: operator(.NEfp.)

  use CALIBRATION_CONSTANTS_MOD,only: &
  planck_a1, planck_a2, planck_nu &
  , sat_name, solar_ch20, ew_ch20, solar_ch20_nu, ch1_dark_count &
  , Launch_Date, Band2_Correction_Start_Date, Band2_Correction_Factor &
  , ABI_FPT_Thresh_038um, ABI_FPT_Thresh_062um, ABI_FPT_Thresh_067um &
  , ABI_FPT_Thresh_073um, ABI_FPT_Thresh_085um, ABI_FPT_Thresh_097um &
  , ABI_FPT_Thresh_104um, ABI_FPT_Thresh_110um, ABI_FPT_Thresh_120um &
  , ABI_FPT_Thresh_133um

  use PLANCK_MOD,only: &
  Planck_Temp_Fast

  use GOES_MOD,only:  &
  gvar_nav &
  , area_struct &
  , Compute_Satellite_Angles

  use FILE_UTILS,only: Get_Lun

  use VIEWING_GEOMETRY_MOD, only: &
  possol
  use CLAVRX_MESSAGE_MOD

  implicit none
  private
  public:: READ_ABI
  public:: READ_NAVIGATION_BLOCK_ABI
  public:: READ_ABI_INSTR_CONSTANTS

  private :: ABI_RADIANCE_BT, ABI_NAVIGATION, ABI_Reflectance, GET_ABI_IMAGE

  type (GVAR_NAV)     :: NAVSTR_ABI_NAV
  integer, PARAMETER  :: NCHAN_ABI= 16
  integer, PARAMETER  :: NDET_ABI = 1
  integer, PARAMETER  :: NTABLE_ABI = 65536

  real (kind=real4), allocatable, dimension(:,:), PRIVATE :: Ref_Table
  integer (kind=int4), allocatable, dimension(:,:,:), PRIVATE :: bt_table
  integer (kind=int4), allocatable, dimension(:,:,:), PRIVATE :: rad_table

  integer(kind=int4), private, parameter:: ABI_Xstride = 1
  integer(kind=int4), private, parameter:: num_4km_scans_fd = 5424
  integer(kind=int4), private, parameter:: num_4km_elem_fd = 5424
  integer(kind=int4), private, parameter:: time_for_fd_scan =  636000 !milliseconds (10.6 min,)
  real, private, save:: Scan_rate    !scan rate in millsec / line
  integer(kind=int4), private, parameter:: IN3_Byte_Shift = 0 !number of bytes to shift for FY2

  integer(kind=int4), dimension(64) :: i4buf

  character (len = 1020) :: file_fullpath_list(16)

CONTAINS

  !----------------------------------------------------------------
  ! read the ABI constants into memory
  ! This public routine is called in sensor_mod.f90
  !-----------------------------------------------------------------
  subroutine READ_ABI_INSTR_CONSTANTS(Instr_Const_file)
    character(len=*), intent(in):: Instr_Const_file
    integer:: ios0, erstat
    integer:: Instr_Const_lun
    character(len=70):: header
    real:: dummy_rad_to_ref

    Instr_Const_lun = GET_LUN()

    open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

    erstat = 0
    if (ios0 /= 0) then
      erstat = 19
      call MESG( "ABI instrument constants open failed "//trim(Instr_Const_file),level = verb_lev % ERROR)
      stop 19
    end if
    read(unit=Instr_Const_lun,fmt="(a7)") sat_name
    read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
    read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
    read(unit=Instr_Const_lun,fmt=*) Launch_Date
    read(unit=Instr_Const_lun,fmt=*) Band2_Correction_Start_Date
    read(unit=Instr_Const_lun,fmt=*) Band2_Correction_Factor
    read(unit=Instr_Const_lun,fmt=*) header
    read(unit=Instr_Const_lun,fmt=*) planck_a1(20), planck_a2(20), planck_nu(20)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(27), planck_a2(27), planck_nu(27)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(28), planck_a2(28), planck_nu(28)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(29), planck_a2(29), planck_nu(29)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(30), planck_a2(30), planck_nu(30)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(31), planck_a2(31), planck_nu(31)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(32), planck_a2(32), planck_nu(32)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(33), planck_a2(33), planck_nu(33)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(37), planck_a2(37), planck_nu(37)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(38), planck_a2(38), planck_nu(38)
    read(unit=Instr_Const_lun,fmt=*) header
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) dummy_rad_to_ref
    read(unit=Instr_Const_lun,fmt=*) header
    read(unit=Instr_Const_lun,fmt=*) ABI_FPT_Thresh_038um
    read(unit=Instr_Const_lun,fmt=*) ABI_FPT_Thresh_062um
    read(unit=Instr_Const_lun,fmt=*) ABI_FPT_Thresh_067um
    read(unit=Instr_Const_lun,fmt=*) ABI_FPT_Thresh_073um
    read(unit=Instr_Const_lun,fmt=*) ABI_FPT_Thresh_085um
    read(unit=Instr_Const_lun,fmt=*) ABI_FPT_Thresh_097um
    read(unit=Instr_Const_lun,fmt=*) ABI_FPT_Thresh_104um
    read(unit=Instr_Const_lun,fmt=*) ABI_FPT_Thresh_110um
    read(unit=Instr_Const_lun,fmt=*) ABI_FPT_Thresh_120um
    read(unit=Instr_Const_lun,fmt=*) ABI_FPT_Thresh_133um
    close(unit=Instr_Const_lun)

    !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
    Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / ew_Ch20

    !-- hardwire ch1 dark count
    Ch1_Dark_Count = 29

    call MESG( "ABI instrument constants read in successfully",level = verb_lev % DEFAULT)

  end subroutine READ_ABI_INSTR_CONSTANTS

  !-------------------------------------------------------------------------------
  ! Public routine to read data from an AREA file for one segment into memory
  !  This is called in sensor_mod.f90
  !-------------------------------------------------------------------------------
  subroutine READ_ABI(segment_number,channel_1_filename, &
    jday, image_time_ms, &
    AREAstr,NAVstr_ABI)

    integer(kind=int4), intent(in):: segment_number
    character(len=*), intent(in):: channel_1_filename
    type (AREA_STRUCT), intent(in) :: AREAstr
    type (GVAR_NAV), intent(in)    :: NAVstr_ABI
    integer, intent(in):: jday
    integer(kind=int4), intent(in):: image_time_ms

    character(len=1020):: channel_x_filename
    character(len=1020):: channel_x_filename_full
    character(len=1020):: channel_x_filename_full_uncompressed
    character(len=4096) :: cmd
    integer:: ipos, ierr, nc
    integer:: ilen
    integer:: ichan_goes
    integer:: ichan_modis
    integer:: ABI_file_id
    real(kind=real4):: Image_Time_Hours
    real(kind=real4):: Image_Date
    integer(kind=int4):: Image_Jday
    integer(kind=int4):: first_line_in_segment
    character(len=2):: ichan_goes_string
    integer :: Line_Idx
    integer :: Elem_Idx
    integer:: num_elements_this_image
    character(len=200) :: l1b_path


    logical :: is_solar_channel(16)
    character(len=1020) :: channel_filename_list(16)
    integer, parameter  :: MODIS_CHN_LIST (16) = &
    & [ 3 , 1 , 2 , 26 , 6 , 7 , 20 , 37 ,  27 , 28,  &
    29 , 30 , 38 , 31 , 32 , 33 ]
    !-------------------------------------------------------------------
    ! - executable
    ! --------------------------------------

    if (.not. allocated(Ref_Table)) allocate(Ref_table(6,Ntable_ABI))
    
    !--- assume channel_1_file name has a unique "_1_" in the name.
    !--- determine indices needed to replace that string
    ipos = index(channel_1_filename, "_1_")

    ! - we accept now both _1_ file which we never use because it is 0.47 channel
    ! - and the Visible 0.64 channel _3_
    ! - AW 2019 March 9
    if ( ipos .eq. 0 ) then
      ipos = index(channel_1_filename, "_3_")
    end if

    if ( ipos .eq. 0 ) then
      print*, 'ABI file in file_list mut be either channel 1 or channel 3'
      print*, ' Please fix in file_list. stopping'
      stop
    end if

    ilen = len(channel_1_filename)

    first_line_in_segment = (segment_number-1)*Image%Number_Of_Lines_Per_Segment


    is_solar_channel(7:16) = .false.
    is_solar_channel(1:6) = .true.

    !---------------------------------------------------------------------------
    ! ABI Navigation (Do Navigation and Solar angles first)
    !---------------------------------------------------------------------------

    call ABI_NAVIGATION(1,first_line_in_segment,&
    Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment,1,&
    AREAstr,NAVstr_ABI)

    if (Segment_Number == 1) then

      ! populate channel filename list
      do Ichan_Goes = 1,16
        write(ichan_goes_string,fmt="(I2.0)") ichan_goes
        channel_filename_list(ichan_goes) = channel_1_filename(1:ipos-1) // "_"//trim(adjustl(ichan_goes_string))//"_" // &
        channel_1_filename(ipos+3:ilen)
        l1b_path = trim(Image%Level1b_Path)
        if ( l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES )   l1b_path =  trim(Temporary_Data_Dir)
        file_fullpath_list(ichan_goes) = trim( l1b_path)//trim(channel_filename_list(ichan_goes))
      end do


      Ref_Table = missing_value_int4

      ! deal with compressed l1b files
      if  (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES ) then
        do ichan_goes = 1,16
          ichan_modis = MODIS_CHN_LIST(ichan_goes)

          if (Sensor%Chan_On_Flag_Default(ichan_modis) == sym%YES) then
            write(ichan_goes_string,fmt="(I2.0)") ichan_goes

            channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_filename_list(ichan_goes))

            channel_x_filename_full_uncompressed = trim(Image%Level1b_Path)//trim(channel_filename_list(ichan_goes))

            if (l1b_gzip == sym%YES) then
              cmd = "gunzip -c "//trim(channel_x_filename_full_uncompressed)//".gz"// &
              " > "//trim(channel_x_filename_full)

            else if  (l1b_bzip2 == sym%YES) then
              cmd = "bunzip2 -c "//trim(channel_x_filename_full_uncompressed)//".bz2"// &
              " > "//trim(channel_x_filename_full)

            end if

            nc = len_trim(cmd)
            call univ_system_cmd_f(nc, trim(cmd), ierr)

            Number_of_Temporary_Files = Number_of_Temporary_Files + 1
            Temporary_File_Name(Number_of_Temporary_Files) = trim(channel_x_filename)

          end if  ! - channel on condition

        end do ! - channel loop
      end if ! - end compressed data


      !--- On first segment grab the calibration block from the AREA file.
      !--- Need to read each AREA file separately, as the calibration block is
      !--- unique for each AREA file.

      do ichan_goes = 1,16
        ichan_modis = MODIS_CHN_LIST(ichan_goes)

        if (Sensor%Chan_On_Flag_Default(ichan_modis) == sym%NO) cycle
        ABI_file_id = get_lun()
        call mesg ("Channel 1 calibration for : "// trim(channel_filename_list(ichan_goes)), level = 6 )

        call mread_open(trim(file_fullpath_list(ichan_goes))//CHAR(0), ABI_file_id)

        call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
        call mread_close(ABI_file_id)
      end do

    end if !--- Segment 1
    !--- END CALIBRATION BLOCK READ

    !--- READ OF AREA FILE FOR EACH SEGMENT
    do ichan_goes = 1,16
      ichan_modis = MODIS_CHN_LIST(ichan_goes)
      if (Sensor%Chan_On_Flag_Default(ichan_modis) == sym%NO) cycle

      call GET_ABI_IMAGE(trim(file_fullpath_list(ichan_goes)), &
      AREAstr, &
      Segment_Number, &
      Image%Number_Of_Lines_Per_Segment, &
      Image%Number_Of_Lines_Read_This_Segment, &
      Two_Byte_Temp)

      !--- Make reflectance table.
      if (is_solar_channel(ichan_goes)) then
        call ABI_Reflectance(Two_Byte_Temp,ch(ichan_modis)%Ref_Toa(:,:),ichan_goes)
      else if (.not. is_solar_channel(ichan_goes)) then
        !--- Make radiance table.
        call ABI_RADIANCE_BT(ichan_goes, Two_Byte_Temp, ch(ichan_modis)%Rad_Toa, ch(ichan_modis)%Bt_Toa)
      end if

    end do

    !--- perform Band-2 Correction (Band 2 = Ch 1 in Clavrx)
    if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
      !Image_Date = Image%Start_Year + (Image%Start_Doy - 1.0) / 365.25
      Image_Date = image % time_start % year + (image % time_start %  dayOfYear - 1)/365.25
      if (Image_Date < Band2_Correction_Start_Date) then
        where (ch(1)%Ref_Toa .NEfp. Missing_Value_Real4)
          ch(1)%Ref_Toa = Band2_Correction_Factor*ch(1)%Ref_Toa
        endwhere
      endif
    endif

    !--- compute scan rate
    num_elements_this_image =  int(AREAstr%num_elem / ABI_Xstride) + 1
    Scan_rate = real((num_elements_this_image)/               &
    & real(num_4km_elem_fd/ABI_Xstride)) * &
    & real((AREAstr%num_line) / real(num_4km_scans_fd)) * &
    & real(time_for_fd_scan) / real(AREAstr%num_line)

    do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
      Image%Scan_Number(Line_Idx) = first_line_in_segment + Line_Idx
      Image%Scan_Time_Ms(Line_Idx) = Image_Time_Ms + (Image%Scan_Number(Line_Idx)-1) * Scan_rate
    end do

    !------------------------------------------------------------------------------
    ! ABI Angles
    ! NOTE: These were private routines in the GOES module. Suggest they become
    !       public with different names, since they are used cross platform
    !------------------------------------------------------------------------------
    image_jday = jday
    image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0

    do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
      do Elem_Idx = 1,Image%Number_Of_Elements
        call POSSOL(image_jday,image_time_hours, &
        Nav%Lon_1b(Elem_Idx,Line_Idx),Nav%Lat_1b(Elem_Idx,Line_Idx), &
        Geo%Solzen(Elem_Idx,Line_Idx),Geo%Solaz(Elem_Idx,Line_Idx))
      end do
      call COMPUTE_SATELLITE_ANGLES(Sensor%Geo_Sub_Satellite_Longitude,  &
      Sensor%Geo_Sub_Satellite_Latitude, Line_Idx)
    end do

    !--- Ascending node
    Elem_Idx = Image%Number_Of_Elements/2
    do Line_Idx = Line_Idx_Min_Segment+1, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
      Nav%Ascend(Line_Idx) = 0
      if (Nav%Lat_1b(Elem_Idx,Line_Idx) < Nav%Lat_1b(Elem_Idx,Line_Idx-1)) then
        Nav%Ascend(Line_Idx) = 1
      end if
    end do
    Nav%Ascend(Line_Idx_Min_Segment) = Nav%Ascend(Line_Idx_Min_Segment+1)

  end subroutine READ_ABI

  !-------------------------------------------------------------------------------
  ! LOAD ABI CALIBRATION
  !
  !  For single band AREA files, the calibration block is distinct for each
  !  channel.  It must be called for all 16 channels to load the calibration
  !  tables.
  !
  !  The calibration block is an (18,16) array.  Second array index is the band.
  !  The below example will be for band 1.
  !
  !   Channel Number     -> (1,1)
  !   Channel wavelength -> (2,1)
  !   Scale Factor       -> (3,1)
  !   Add Offset         -> (4,1)
  !   Radianc to Albedo  -> (12,1)
  !   FK1                -> (14,1)
  !   FK2                -> (15,1)
  !   BC1                -> (16,1)
  !   BC2                -> (17,1)
  !   Invalid Value      -> (18,1)
  !
  ! Current calibration source type is 'ABIN'
  !
  ! Band Number Filter Map for single band AREA files.
  ! From the AREA directory, index 19.
  !  Band 1  -> 1
  !  Bamd 2  -> 2
  !  Band 3  -> 4
  !  Band 4  -> 8
  !  Band 5  -> 16
  !  Band 6  -> 32
  !  Band 7  -> 64
  !  Band 8  -> 128
  !  Band 9  -> 256
  !  Band 10 -> 512
  !  Band 11 -> 1024
  !  Band 12 -> 2048
  !  Band 13 -> 4096
  !  Band 14 -> 8192
  !  Band 15 -> 16384
  !  Band 16 -> 32768
  !
  ! REFERENCE - McIDAS-X 2016.1 kbxabin.dlm
  !
  !-------------------------------------------------------------------------------
  subroutine load_ABI_calibration(lun, AREAstr, ichan_modis)
    integer(kind=int4), intent(in) :: lun, ichan_modis
    type(AREA_STRUCT), intent(in):: AREAstr
    integer :: i
    integer(kind=int4) :: local_sndr_filter_map
    integer, parameter :: calb_size = 289
    integer(kind=int4), parameter :: MAXVAL = 65536
    integer(kind=int4), dimension(18,16) :: ibuf
    integer(kind=int4) :: band
    integer(kind=int4) :: iCNT
    integer(kind=int4) :: calb_bandNo
    integer(kind=int4) :: calb_errorCount
    integer(kind=int4) :: calb_invalidValue
    real(kind=real8) :: dRAD
    real(kind=real4) :: rRAD
    real(kind=real8) :: dALB
    real(kind=real4) :: rALB
    real(kind=real4) :: rTEMP
    real(kind=real8) :: calb_gainCnt2rad
    real(kind=real8) :: calb_cnstCnt2rad
    real(kind=real8) :: calb_rad2albedo

    if (.not. allocated(bt_table)) allocate(bt_table(16,ndet_ABI,Ntable_ABI))
    if (.not. allocated(rad_table)) allocate(rad_table(16,ndet_ABI,Ntable_ABI))
    
    !--- Store the correct calibration source type in the AREA structure.  It will be
    !--- important if the calibration type changes.
    call mreadf_int_o(lun,0,4,64,i4buf)
    call move_bytes(4,i4buf(1+51:1+51),AREAstr%src_Type,0)
    local_sndr_filter_map = i4buf(19) ! Need band map for current open AREA file.

    if (trim(AREAstr%src_Type) .eq. 'ABIN') then
      !print*,"ABIN calibration, proceeding ...", AREAstr%num_chan, local_sndr_filter_map, ichan_modis
    else
      print*,"Unknown ABI calibration type, exiting ..."
      stop
    end if

    !--- Single band images are required.
    if (AREAstr%num_chan > 1) then
      print*,"Multi-band image, exiting ..."
      stop
    end if

    !--- Initialize these variables
    calb_errorCount = -999
    calb_invalidValue = -999

    !--- Extract calibration block per image.
    select case(local_sndr_filter_map)
      !--- Bands 1-6 are albedo/reflectance.
    case(1,2,4,8,16,32)

      !--- Set band number
      if (local_sndr_filter_map == 1)  band=1
      if (local_sndr_filter_map == 2)  band=2
      if (local_sndr_filter_map == 4)  band=3
      if (local_sndr_filter_map == 8)  band=4
      if (local_sndr_filter_map == 16) band=5
      if (local_sndr_filter_map == 32) band=6

      !--- Read the calibration block.
      call mreadf_int_o(lun,AREAstr%cal_offset+4,4,calb_size,ibuf)

      !--- Fill variables with calibration block information.
      calb_bandNo       = ibuf(1,band)
      calb_gainCnt2rad  = DBLE(ibuf(3,band)) * 0.0000001D0
      calb_cnstCnt2rad  = DBLE(ibuf(4,band)) * 0.0000001D0
      calb_rad2albedo   = DBLE(ibuf(12,band)) * 0.0000001D0
      calb_errorCount   = ibuf(18,band)
      calb_invalidValue = ibuf(18,band)

      !--- Loop through MAXVAL RAW counts to build tables.
      do i = 1, MAXVAL

        !--- Count value is table index.
        iCNT = i - 1

        !--- Convert RAW to RAD
        if (iCNT == calb_errorCount) then
          dRAD = calb_invalidValue
        else
          dRAD = dble(iCNT)*calb_gainCnt2rad+calb_cnstCnt2rad
          if (dRAD <= 0.0D0 ) then
            dRAD = calb_invalidValue
          endif
        endif

        !--- Convert RAD to ALB.
        if (dRAD == calb_invalidValue) then
          dALB = 0
          rALB = 0.0
        else
          dALB = calb_rad2albedo*dRAD
          rALB = REAL(dALB) * 100.0   ! To match McIDAS.
        endif

        !--- Fill lookup table with albedo.
        Ref_Table(band,i) = rAlb

      end do

      !--- Bands 7-16 are thermal or mixed.
    case(64,128,256,512,1024,2048,4096,8192,16384,32768)

      !--- Set band number
      if (local_sndr_filter_map == 64)    band=7
      if (local_sndr_filter_map == 128)   band=8
      if (local_sndr_filter_map == 256)   band=9
      if (local_sndr_filter_map == 512)   band=10
      if (local_sndr_filter_map == 1024)  band=11
      if (local_sndr_filter_map == 2048)  band=12
      if (local_sndr_filter_map == 4096)  band=13
      if (local_sndr_filter_map == 8192)  band=14
      if (local_sndr_filter_map == 16384) band=15
      if (local_sndr_filter_map == 32768) band=16

      !--- Read the calibration block.
      call mreadf_int_o(lun,AREAstr%cal_offset+4,4,calb_size,ibuf)

      !--- Fill variables with calibration block information.
      calb_bandNo       = ibuf(1,band)
      calb_gainCnt2rad  = DBLE(ibuf(3,band)) * 0.0000001D0
      calb_cnstCnt2rad  = DBLE(ibuf(4,band)) * 0.0000001D0

      !--- Loop through MAXVAL RAW counts to build tables.
      do i = 1, MAXVAL

        !--- Count value is table index.
        iCNT = i - 1

        !--- Convert RAW to RAD
        if (iCNT == calb_errorCount) then
          dRAD = calb_invalidValue
        else
          dRAD = dble(iCNT)*calb_gainCnt2rad+calb_cnstCnt2rad
          if (dRAD <= 0.0D0 ) then
            dRAD = calb_invalidValue
          endif
        endif

        !--- Convert RAD to TEMP.
        if (dRAD == calb_invalidValue) then
          dRAD =  calb_invalidValue
          rTEMP = calb_invalidValue
        else
          rRAD = REAL(dRAD)
          rTEMP = PLANCK_TEMP_FAST(ichan_modis,rRAD)
        endif

        !--- Fill look up tables with RAD and TEMP. Scale to integers.
        if (dRAD == calb_invalidValue) then
          bt_table(band,1,i) = NINT(Missing_Value_Real4)
          rad_table(band,1,i) = NINT(Missing_Value_Real4)
        else
          bt_table(band,1,i) = NINT(rTEMP * 100.0)
          rad_table(band,1,i) = NINT(dRAD * 1000.0)
        endif

      end do

    end select

  end subroutine load_ABI_calibration

  ! Perform ABI Navigation

  subroutine ABI_NAVIGATION(xstart,ystart,xsize,ysize,xstride, &
    AREAstr,NAVstr_ABI)
    integer(kind=int4) :: xstart, ystart
    integer(kind=int4) :: xsize, ysize
    integer(kind=int4) :: xstride
    type (AREA_STRUCT) :: AREAstr
    type (GVAR_NAV), intent(in)    :: NAVstr_ABI

    integer :: i, j, ii, jj, imode
    real(kind(0.0d0)) :: latitude, longitude
    real(kind=real4) :: height
    integer :: FGF_type = 1 ! GOES-16

    NAVstr_ABI_NAV = NAVstr_ABI

    imode = -1
    height = 0.0   !Used for parallax correction

    Nav%Lat = Missing_Value_Real4
    Nav%Lon = Missing_Value_Real4

    if (NAVstr_ABI%nav_type == 'ABIN') then

      do j=1, ysize

        !--- For 2 km full disk, we need navigation transformations in 1 km
        !--- space, as Band 1 is the active AREA file at this point.  However,
        !--- lcor/ecor are in 1/2 km space.
        !--- ORIGINAL jj = ystart + (j-1)*AREAstr%line_res
        jj = ystart*AREAstr%line_res + (j-1)*AREAstr%line_res
        jj = ((jj + AREAstr%north_bound + (NAVstr_ABI%BRES-1)))/NAVstr_ABI%BRES


        do i=1, xsize

          !ii = (xstart + (i-1)*xstride + AREAstr%west_vis_pixel + (NAVstr_ABI%BRES-1))/NAVstr_ABI%BRES
          ii = xstart + (i-1)*xstride*AREAstr%elem_res
          ii = ((ii + AREAstr%west_vis_pixel + (NAVstr_ABI%BRES-1)))/NAVstr_ABI%BRES

          CALL fgf_to_earth(FGF_type,                  &
          DBLE(ii),                  &
          DBLE(jj),                  &
          DBLE(NAVstr_ABI%abiCFAC),    &
          DBLE(NAVstr_ABI%abiCOFF),    &
          DBLE(NAVstr_ABI%abiLFAC),    &
          DBLE(NAVstr_ABI%abiLOFF),    &
          DBLE(NAVstr_ABI%sub_lon), &
          longitude,                 &
          latitude)


          if (latitude .LE. -999.0) then  ! -999.99 is MSV nav missing value
            Nav%Lat_1b(i,j) = Missing_Value_Real4
            Nav%Lon_1b(i,j) = Missing_Value_Real4
            Geo%Space_Mask(i,j) = .TRUE.
          else
            Nav%Lat_1b(i,j) = real(latitude,kind=real4)
            Nav%Lon_1b(i,j) = real(longitude,kind=real4)

            ! we want 180 to -180, one last check.
            if (longitude .GT. 180.0 ) then
              Nav%Lon_1b(i,j) = real(longitude,kind=real4) - 360.0
            endif

            Geo%Space_Mask(i,j) = .FALSE.

          endif

        end DO

      end do

    endif

  end subroutine ABI_NAVIGATION

  !------------------------------------------------------------------
  ! subroutine to convert ABI counts to radiance and brightness
  ! temperature
  !------------------------------------------------------------------

  subroutine ABI_RADIANCE_BT(chan_num, ABI_Counts, rad2, temp1)

    integer (kind=INT2), dimension(:,:), intent(in):: ABI_Counts
    integer (kind=int4), INTENT(in) :: chan_num
    real (kind=real4), dimension(:,:), INTENT(out):: temp1, rad2

    integer :: i, j, index

    do j = 1, Image%Number_Of_Lines_Read_This_Segment
      do i = 1, Image%Number_Of_Elements
        if ( .NOT. Geo%Space_Mask(i,j) ) then
          index = int(ABI_Counts(i,j),kind=int2)
          rad2(i,j) = real(rad_table(chan_num,1,index+1),kind=real4)/1000.0
          temp1(i,j) = real(bt_table(chan_num,1,index+1),kind=real4)/100.0
        else
          rad2(i,j) = Missing_Value_Real4
          temp1(i,j) = Missing_Value_Real4
        endif
      end do
    end do

  end subroutine ABI_RADIANCE_BT

  !------------------- ABI NAV BLOC ------------------------------
  !--- From McIDAS-X 2016.2 nvxabin.dlm
  !--- Navigation Block:
  !---   1: ABIN -> Navigation module string name.
  !---   2: LOFF
  !---   3: COFF
  !---   4: LFAC
  !---   5: CFAC
  !---   6: Satellite SubPoint Longitude
  !---   7: Base Image Resolution.
  !---
  !--- The call from sensor_module.f90 is for Band 1.  The data is
  !--- at 2km. The navigation module will take care of the necessary
  !--- transformations to account for 2km data.
  !---  LOFF = 15185800
  !---  COFF = -15185800
  !---  LFAC = -2800
  !---  CFAC = 2800
  !-------------------------------------------------------------

  subroutine READ_NAVIGATION_BLOCK_ABI(filename, AREAstr, NAVstr)
    CHARACTER(len=*), intent(in):: filename
    type(AREA_STRUCT), intent(in):: AREAstr
    type(GVAR_NAV), intent(inout):: NAVstr

    integer :: geos_nav
    integer(kind=int4)nav_offset
    integer:: number_of_words_read
    integer(kind=int4), dimension(640) :: i4buf

    nav_offset = AREAstr%sec_key_nav

    ! ABI navigation is names ABIN, but will use the GEOS transforms
    ! from geos_transform_pix.c

    ! Read the navigation block into the buffer.
    geos_nav = sym%NO
    call mreadf_int(trim(filename)//CHAR(0),nav_offset,4,640,&
    number_of_words_read, i4buf)
    call move_bytes(4,i4buf(1),NAVstr%nav_type,0)

    ! LOFF, COFF, CFAC, LFAC stored in McIDAS header for 2km data,
    ! assuming that we are basing the navigation off of Band 1.
    NAVstr%abiLOFF=real(i4buf(2),kind=real8) / 100000000.0
    NAVstr%abiCOFF=real(i4buf(3),kind=real8) / 100000000.0
    NAVstr%abiLFAC=real(i4buf(4),kind=real8) / 100000000.0
    NAVstr%abiCFAC=real(i4buf(5),kind=real8) / 100000000.0
    NAVstr%sub_lon = real(i4buf(6),kind=real4) / 10       ! SUBLON stored as SUBLON *10 in McIDAS NAV block
    NAVstr%sublon = NAVstr%sub_lon
    NAVstr%BRES=i4buf(7)

  end subroutine READ_NAVIGATION_BLOCK_ABI
  ! -------------------------------------------------------------------------
  ! Perform ABI Reflectance calculation.  These vaules have been
  ! stored in an array, previous to this call.
  !  ----------------------------------------------------------------
  subroutine ABI_Reflectance(ABI_Counts, Alb_Temp,chan_num)

    integer (kind=INT2), dimension(:,:), intent(in):: ABI_Counts
    integer (kind=INT4), intent(in):: chan_num
    real (kind=real4), dimension(:,:), intent(out):: Alb_Temp

    integer :: i, j
    integer :: index


    do j=1, Image%Number_Of_Lines_Read_This_Segment
      do i=1, Image%Number_Of_Elements

        if (.NOT. Geo%Space_Mask(i,j)  .and. Geo%solzen(i,j) < 90.0) THEN

          index = int(ABI_Counts(i,j),kind=int2)

          Alb_Temp(i,j) = Missing_Value_Real4

          if ((index .GT. 0) .AND. (index .LE. Ntable_ABI)) then
            Alb_Temp(i,j) = (real(Ref_Table(chan_num,index+1),kind=real4))
          end if

        else
          Alb_Temp(i,j) = Missing_Value_Real4
        end if

      end do
    end do


  end subroutine ABI_Reflectance


  !
  !
  !
  subroutine GET_ABI_IMAGE(filename,AREAstr, &
    segment_number, &
    num_lines_per_segment, &
    num_lines_read, image)

    character(len=*), intent(in):: filename
    type (AREA_STRUCT), intent(in) :: AREAstr
    integer(kind=int4), intent(in):: segment_number
    integer(kind=int4), intent(in):: num_lines_per_segment
    integer(kind=int4), intent(out):: num_lines_read
    integer(kind=int2), dimension(:,:), intent(out):: image

    integer(kind=int4):: bytes_per_line
    integer(kind=int4):: words_in_prefix
    integer(kind=int4):: words_per_line
    integer(kind=int4):: first_line_in_segment
    integer(kind=int4):: last_line_in_segment
    integer(kind=int4):: first_byte_in_segment
    integer(kind=int4):: number_of_words_in_segment
    integer(kind=int4):: bytes_per_pixel
    integer(kind=int4):: num_byte_ln_prefix
    integer(kind=int4):: pri_key_nav
    integer(kind=int4):: number_of_words_read
    integer(kind=int4):: word_start
    integer(kind=int4):: dummy
    integer(kind=int4):: word_end
    integer(kind=int4):: bytemove
    integer(kind=int4):: Line_Idx
    integer(kind=int2), dimension(:), allocatable:: word_buffer
    integer (kind=int1), dimension(:), allocatable :: buffer1
    integer(kind=int4), dimension(64) :: i4buf_temp
    integer:: nwords

    bytemove = 0 !Not sure what this needs to be for ABI?
    image = 0

    ! get number of bytes per pixel and num bytes per line for current file
    ! this is needed because the 0.64 and other channels have different values
    ! in the COMS HIRID format.

    call mreadf_int(trim(filename)//CHAR(0),0,4,64,dummy,i4buf_temp)
    bytes_per_pixel = i4buf_temp(11)
    num_byte_ln_prefix = i4buf_temp(15)
    pri_key_nav = AREAstr%pri_key_nav
    bytes_per_line = num_byte_ln_prefix + (AREAstr%num_elem*bytes_per_pixel)
    words_in_prefix = num_byte_ln_prefix / bytes_per_pixel
    words_per_line = words_in_prefix + AREAstr%num_elem
    first_line_in_segment = (segment_number-1)*num_lines_per_segment + 1
    last_line_in_segment = min(AREAstr%num_line,segment_number*num_lines_per_segment)
    first_byte_in_segment = pri_key_nav + &
    bytes_per_pixel*(first_line_in_segment-1) * words_per_line + &
    bytes_per_pixel
    number_of_words_in_segment = words_per_line * num_lines_per_segment

    allocate(word_buffer(number_of_words_in_segment))

    if ( bytes_per_pixel == 2) then
      call mreadf_int(filename//CHAR(0), &
      first_byte_in_segment,  &
      bytes_per_pixel, &
      number_of_words_in_segment, &
      number_of_words_read,word_buffer)
    endif

    if( bytes_per_pixel == 1) then
      allocate(buffer1(number_of_words_in_segment))

      call mreadf_int(filename//CHAR(0), &
      first_byte_in_segment,  &
      bytes_per_pixel, &
      number_of_words_in_segment, &
      number_of_words_read,buffer1)

      word_buffer = int(buffer1,kind=int2)
      !Since 1 byte values are always signed in Fortran, I need to add 256 to the
      !negative values
      where (word_buffer < 0)
        word_buffer = word_buffer + 256
      end where
    endif

    !--- update number of scans read
    num_lines_read = number_of_words_read / words_per_line

    do Line_Idx = 1, num_lines_read
      word_start = (words_in_prefix + AREAstr%num_elem)*(Line_Idx-1) + words_in_prefix + 1
      word_end = min(word_start + AREAstr%num_elem - 1,number_of_words_in_segment)
      nwords = int(word_end - word_start)/ABI_Xstride + 1
      image(1:nwords,Line_Idx) = ishft(word_buffer(word_start:word_end:1),bytemove)
    enddo

    deallocate(word_buffer)
    if (allocated(buffer1)) deallocate(buffer1)

  end subroutine GET_ABI_IMAGE

end module ABI_MOD
