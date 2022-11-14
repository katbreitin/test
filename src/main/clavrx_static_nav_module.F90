!$Id: clavrx_static_nav_module.f90 4128 2021-04-19 02:03:04Z heidinger $
!------------------------------------------------------------------------------
!
!   Stativ nav file has the follwoing tasks
!
!      FOR ABI/AHI ONLY
!
!      These global variables are populated:
!
!    
!
!------------------------------------------------------------------------------
module CLAVRX_STATIC_NAV_MODULE

  use CONSTANTS_MOD, only: &
    int4 &
    , int1 &
    , sym &
    , THERMAL_OBS_TYPE &
    , MIXED_OBS_TYPE &
    , MISSING_VALUE_REAL4 &
    , Missing_Value_Int4 &
    , Missing_Value_Int1 &
    , SOLAR_OBS_TYPE
    
  use CX_NETCDF4_MOD,only:  &
      CLOSE_NETCDF &
    , OPEN_NETCDF &
    , read_abi_kappa0 &
    , READ_ABI_MAX_FOCAL_PLANE_TEMP &
    , read_abi_time &
    , read_ahi_time &
    , read_and_unscale_netcdf_2d &
    , READ_NETCDF &
    , READ_NETCDF_dimension &
    , READ_NETCDF_GLOBAL_ATTRIBUTE &
    , READ_NETCDF_VARIABLE_ATTRIBUTE_REAL
  
  use PIXEL_COMMON_MOD, only: &
      Ch &
    , Sensor &
    , Image &
    , Geo &
    , Nav &
    , Static_Nav_File &
    , Ancil_Data_Dir &
    , X_Sample_Offset &
    , Y_Sample_Offset &
    , Static_Dark_Sky_Flag &
    , QRNN_CTP_Flag &
    , QRNN_CTP_Data_Dir &
    , QRNN_CTP_Name &
    , QRNN_CTP &
    , QRNN_CTP_Uncer &
    , Dark_Comp_Data_Dir &
    , Dark_Composite_Name &
    , Static_Ref_065um_Dark_Composite &
    , Static_Ref_065um_Dark_Composite_Stddev &
    , Diag_Pix_Array_1 &
    , Diag_Pix_Array_2 &
    , Diag_Pix_Array_3

  use PLANCK_MOD, only: PLANCK_TEMP_FAST, CONVERT_RADIANCE
  
  use CALIBRATION_CONSTANTS_MOD, only: Band2_Correction_Factor, &
    Planck_Nu &
    , Rad_to_Ref_Fac_2_10um &
    , Rad_to_Ref_Fac_1_60um &
    , Rad_to_Ref_Fac_1_38um &
    , Rad_to_Ref_Fac_0_86um &
    , Rad_to_Ref_Fac_0_65um &
    , Rad_to_Ref_Fac_0_55um &
    , Rad_to_Ref_Fac_0_47um &
    , band2_correction_start_date 

  use CX_DATE_TIME_TOOLS_MOD, only: compute_time_hours

  use CLAVRX_MESSAGE_MOD
  
  use VIEWING_GEOMETRY_MOD, only: RELATIVE_AZIMUTH, SCATTERING_ANGLE, &
                                  GLINT_ANGLE, POSSOL

  use CX_REAL_BOOLEAN_MOD

#ifdef LIBHIM
 ! HSD reader
 USE AHI_HSD_READER, only : GET_HIM_RADDATA, FILE_TYPE_AHI_HSF, FILE_TYPE_AHI_HRIT
#endif

  implicit none
  private
  public:: SETUP_READ_LEVEL1B_FIXED_GRID_STATIC_NAV
  public:: READ_HSD_FIXED_GRID_STATIC_NAV
  public:: READ_LEVEL1B_FIXED_GRID_STATIC_NAV
  public:: DESTROY_READ_LEVEL1B_FIXED_GRID_STATIC_NAV
  public:: DETERMINE_BOUNDS_STATIC_NAV
  public:: DETERMINE_BOUNDS_STATIC_NAV_ORIGINAL
  public:: READ_SEGMENT_STATIC_NAV
  private:: READ_SEGMENT_LEVEL1B_VER2
  private:: FIND_DARK_COMPOSITE
  private:: READ_DARK_COMPOSITE
  private:: READ_QRNN_CTP

  integer, save:: Num_Elem_fd, Elem_Start, Num_Elem
  integer,  save:: Num_Line_fd, Line_Start, Num_Line
  real,  save:: Sub_Sat_Lat, Sub_Sat_Lon
  character(len=120), parameter :: EXE_PROMPT_NAV = "CLAVR-x Static Navigation Module >> "
  real, parameter :: COUNT_MIN = 1.0
  integer, save:: Ncid_Static
  integer, save:: nc_varid_lat, nc_varid_lon, nc_varid_senzen, nc_varid_senaz

  integer, dimension(:), allocatable, save:: CLAVRx_Chan_Map
  integer, save :: Num_Chan_Sensor
  integer,dimension(:),  allocatable, private, save:: Chan_Stride
  integer, dimension(:), allocatable, private, save:: Chan_On_Sensor
  integer, dimension(:), allocatable, private, save:: ncid_level1b
  character(len=20), dimension(:), allocatable, save:: Variable_Name
  character(len=200), dimension(:), allocatable, save:: Level1b_File

  real, allocatable, save, public :: Rad_To_Ref_Fac(:)
  integer (kind=int1), allocatable, dimension(:,:) , save:: DQF_Seg_2km, DQF_Seg_1km, DQF_Seg_500m
  real, allocatable, save :: Output_Seg_2km(:,:)    ! Not public, no need
  real, allocatable, save, public:: Output_Seg_1km(:,:,:)  ! 3 channels
  real, allocatable, save, public:: Output_Seg_500m(:,:,:) ! only 1 channel
  real, allocatable, dimension(:,:), save:: Min_Output_Seg_2km, Min_Output_Seg_500m, Min_Output_Seg_1km
  real, allocatable, dimension(:,:), save:: Max_Output_Seg_2km, Max_Output_Seg_500m, Max_Output_Seg_1km
  real, allocatable, dimension(:,:), save:: Stddev_Output_Seg_2km, Stddev_Output_Seg_500m, Stddev_Output_Seg_1km
  real(kind=8),allocatable, save:: Output_Time(:)
  logical, dimension(:), allocatable, save:: File_Present_Flag

  character(len=1020), private, save:: Static_Nav_Full_File

  logical, private, save :: Dark_File_Exist
  logical, private, save :: QRNN_CTP_File_Exist

  contains

  ! ****************************************************************************
  !  This is called in sensor_mod.f90
  ! ****************************************************************************
  subroutine SETUP_READ_LEVEL1B_FIXED_GRID_STATIC_NAV()

    use FILE_UTILS, only: FILE_SEARCH, FILE_TEST

    integer:: Chan_Idx
    integer:: ipos, ilen, ipos1, ilen1, N_Files
    character(len=5):: Chan_String_AHI
    character(len=4):: Chan_String_ABI
    character(len=2):: Image_Start_String_ABI
    character(len=14):: Image_Time_String_ABI
    character ( len = 1020 ) , pointer  :: File_Arr_Dummy(:)
    character ( len = 1020 ):: Fusion_File_Dummy
    integer :: i_line
    character(len=10):: Orbital_Slot, Scene_Id
    logical :: Fusion_File_Exist
    real:: Geospatial_Lat_Center_L1b, Geospatial_Lon_Center_L1b
    real:: Geospatial_Lat_Center_Sn, Geospatial_Lon_Center_Sn

    !--- determine number of channels and static nav file
    select case (Sensor%WMO_ID)

    case (173:174) !-- ahi
      Num_Chan_Sensor = 16

    case(270:273) !-- abi
      Num_Chan_Sensor = 16

    case default
      call MESG( "STATIC NAV: unknown satellite, stopping", level = verb_lev % ERROR )
      stop

    end select

    !--- allocate memory
    allocate(CLAVRx_Chan_Map(Num_Chan_Sensor))
    allocate(Variable_Name(Num_Chan_Sensor))
    allocate(Chan_Stride(Num_Chan_Sensor))
    allocate(Chan_On_Sensor(Num_Chan_Sensor))
    allocate(Ncid_Level1b(Num_Chan_Sensor))
    allocate(Level1b_File(Num_Chan_Sensor))
    allocate(Rad_To_Ref_Fac(Num_Chan_Sensor))
    allocate(File_Present_Flag(Num_Chan_Sensor))

    !---- sensor specific channel mapping and other parameters
    select case (Sensor%WMO_ID)

    case (173:174) !-- ahi
      CLAVRx_Chan_Map = (/3,4,1,2,6,7,20,37,27,28,29,30,38,31,32,33/)
      Chan_Stride =     (/2,2,4,2,1,1,1,1,1,1,1,1,1,1,1,1/)
      Variable_Name(1:6) =  "RAD"
      Variable_Name(7:16) =  "RAD"
      Chan_String_AHI = "_B01_"
      Rad_To_Ref_Fac(1) =  Rad_to_Ref_Fac_0_47um ! (w/m^2/micron/str)^-1
      Rad_To_Ref_Fac(2) =  Rad_to_Ref_Fac_0_55um ! (w/m^2/micron/str)^-1
      Rad_To_Ref_Fac(3) =  Rad_to_Ref_Fac_0_65um ! (w/m^2/micron/str)^-1
      Rad_To_Ref_Fac(4) =  Rad_to_Ref_Fac_0_86um ! (w/m^2/micron/str)^-1
      Rad_To_Ref_Fac(5) =  Rad_to_Ref_Fac_1_60um ! (w/m^2/micron/str)^-1
      Rad_To_Ref_Fac(6) =  Rad_to_Ref_Fac_2_10um ! (w/m^2/micron/str)^-1
      if (Image%DB_Flag) then
          Chan_String_AHI = "_B03_"
          Chan_Stride =     (/1,1,4,1,1,1,1,1,1,1,1,1,1,1,1,1/)
      endif

    case(270:273) !-- abi
      CLAVRx_Chan_Map = (/3,1,2,26,6,7,20,37,27,28,29,30,38,31,32,33/)
      Chan_Stride = (/2,4,2,1,2,1,1,1,1,1,1,1,1,1,1,1/)
      
      if (Sensor%Spatial_Resolution_Meters .eq. 500) Chan_Stride = (/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
      Variable_Name(1:6) =  "Rad"
      Variable_Name(7:16) =  "Rad"
      Chan_String_ABI = "C01_"
      Image_Start_String_ABI = "_s"
      Rad_To_Ref_Fac(1) = Rad_to_Ref_Fac_0_47um ! (mw/m^2/cm^-1/str)^-1
      Rad_To_Ref_Fac(2) = Rad_to_Ref_Fac_0_65um ! (mw/m^2/cm^-1/str)^-1
      Rad_To_Ref_Fac(3) = Rad_to_Ref_Fac_0_86um ! (mw/m^2/cm^-1/str)^-1
      Rad_To_Ref_Fac(4) = Rad_to_Ref_Fac_1_38um ! (mw/m^2/cm^-1/str)^-1
      Rad_To_Ref_Fac(5) = Rad_to_Ref_Fac_1_60um ! (mw/m^2/cm^-1/str)^-1
      Rad_To_Ref_Fac(6) = Rad_to_Ref_Fac_2_10um ! (mw/m^2/cm^-1/str)^-1

    case default
      call MESG( "STATIC NAV: unknown satellite, stopping", level = verb_lev % ERROR )
      stop

    end select

    !--- turn on the sensor channels
    do Chan_Idx = 1, Num_Chan_Sensor
      Chan_On_Sensor(Chan_Idx) = Sensor%Chan_On_Flag_Default(CLAVRx_Chan_Map(Chan_Idx))
    end do

    !---------------------------------------------
    ! make file names
    !---------------------------------------------
    if (Sensor%WMO_ID == 173 .or. Sensor%WMO_Id == 174) then
      ipos = index(Image%Level1b_Name, Chan_String_AHI)
      ilen = len(Image%Level1b_Name)
    else
      !--- ipos points to the ABI Channel number.
      !--- ipos1 points to the ABI image start time.
      ipos = index(Image%Level1b_Name, Chan_String_ABI)
      if (ipos == 0) then
           print *, 'Level1b reference name: ',trim(image%Level1b_Name)
           print *, 'level1b_name does not have search pattern  '//trim(chan_string_abi)//' in name'
      end if 
      ilen = len(Image%Level1b_Name)
      ipos1 = index(Image%Level1b_Name, Image_Start_String_ABI)
      Image_Time_String_ABI = Image%Level1b_Name(ipos1+2:ipos1+15)
    end if

    !--- initialize this to true, will only be false if a file is missing
    File_Present_Flag = .true.
    Sensor%Fusion_Flag = .false.

    do Chan_Idx = 1, Num_Chan_Sensor

      if (Sensor%WMO_ID == 173 .or. Sensor%WMO_Id == 174) then
        write(Chan_String_AHI(3:4),fmt="(I2.2)") Chan_Idx
      else
        write(Chan_String_ABI(2:3),fmt="(I2.2)") Chan_Idx
      endif

      !--- AHI
      if (Sensor%WMO_ID == 173 .or. Sensor%WMO_Id == 174) then
        Level1b_File(Chan_Idx) = trim(Image%Level1b_Path)//Image%Level1b_Name(1:ipos-1) // Chan_String_AHI // &
                         Image%Level1b_Name(ipos+5:ilen)
      else !--- ABI

        !--- Search for filenames based on the ABI image time (takes care of
        !--- multiple images in a search directory.
        File_Arr_Dummy => FILE_SEARCH(trim(Image%Level1b_Path), &
            Image%Level1b_Name(1:ipos-1) // Chan_String_ABI // '*_s' &
             // Image_Time_String_ABI// '*.nc', N_Files)

        
        !--- if file is not there, set the File_Present_Flag to false and handle later (do not turn off channel)
        if (N_Files == 0) then
           print*,'==============  file not there '//Image%Level1b_Name(1:ipos-1) // Chan_String_ABI // '*_s' &
                   // Image_Time_String_ABI// '*.nc'                 
           File_Present_Flag(Chan_Idx) = .false. 
        endif

        !--- if there are files, assume the desired file is the first one
        if (N_Files /= 0) then
             Fusion_File_Dummy = trim(File_Arr_Dummy(1))
             Fusion_File_Dummy = 'FR' // Fusion_File_Dummy(3:)
             Level1b_File(Chan_Idx) = trim(Image%Level1b_Path)//trim(Fusion_File_Dummy)
             !--- if fusion files exists, use it
             Fusion_File_Exist = FILE_TEST(trim(Level1b_File(Chan_Idx)))
             if (Fusion_File_Exist) then
                call MESG('FOUND FUSION FILE', level = verb_lev % DEFAULT )
                Sensor%Fusion_Flag = .true.
             else !if not, revert to original file
                Level1b_File(Chan_Idx) = trim(Image%Level1b_Path)//File_Arr_Dummy(1)
             endif
        endif

      end if
    end do

    Ncid_Static = -1
    Ncid_Level1b = -1

    ! Set the expected location of the static navigation file.
    Static_Nav_Full_File = trim(Ancil_Data_Dir)//"static/static_nav/"//trim(Static_Nav_File)
    
    if ( .not. file_test ( Static_Nav_Full_File)) then
      print*,'missing file for static nav .. stopping '
      print*,'check this in sensor_mod.f90 COMPUTE_GOES_RU_STATIC_NAV_FILE_NAME'
      print*,trim(Static_Nav_Full_File)
      print*,'Line ',__LINE__,' in ', __FILE__
      stop
    
    end if
    
    ! Open static navigation file.
    call OPEN_NETCDF(Static_Nav_Full_File,Ncid_Static)

    ! Open Level1B NetCDF files. - WCS3
    if (Image%Nc_Format_flag) then
        do Chan_Idx = 1, Num_Chan_Sensor
         if (Chan_On_Sensor(Chan_Idx) == 1 .and. File_Present_Flag(Chan_Idx)) then
            
           call OPEN_NETCDF(Level1b_File(Chan_Idx),Ncid_Level1b(Chan_Idx))
           
         end if
        end do
    end if

    !--- determine if this static nav file is consistent with level1b by looking
    !--- at center of projections.
    if (Sensor%WMO_Id == 270 .or. Sensor%WMO_Id == 271 .or. Sensor%WMO_Id == 272) then
      call READ_NETCDF_GLOBAL_ATTRIBUTE(Static_Nav_Full_File, "geospatial_lat_center", Geospatial_Lat_Center_Sn)
      call READ_NETCDF_GLOBAL_ATTRIBUTE(Static_Nav_Full_File, "geospatial_lon_center", Geospatial_Lon_Center_Sn)
      do Chan_Idx = 1, Num_Chan_Sensor
        if (Chan_On_Sensor(Chan_Idx) == 1 .and. File_Present_Flag(Chan_Idx)) then
          call READ_NETCDF_VARIABLE_ATTRIBUTE_REAL(Level1b_File(Chan_Idx), "geospatial_lat_lon_extent","geospatial_lat_center", Geospatial_Lat_Center_L1b)
          call READ_NETCDF_VARIABLE_ATTRIBUTE_REAL(Level1b_File(Chan_Idx), "geospatial_lat_lon_extent","geospatial_lon_center", Geospatial_Lon_Center_L1b)
          exit
        end if
      end do
      !--- check values use 0.01 as a threshold for difference
      if ((abs(Geospatial_Lon_Center_Sn - Geospatial_Lon_Center_L1b) .gtr. 0.01) .or.  &
          (abs(Geospatial_Lat_Center_Sn - Geospatial_Lat_Center_L1b) .gtr. 0.01)) then
          call MESG( "STATIC NAV not consistent with Level-1b, stopping", level = verb_lev % ERROR )
          print*,Static_Nav_Full_File
          call CLOSE_NETCDF(ncid_static)
          stop 
      endif
    endif

    !--- Read max focal plane temps here.
    if (Sensor%Wmo_Id == 271) then
      do Chan_Idx = 1, Num_Chan_Sensor
        if (Chan_On_Sensor(Chan_Idx) == 1 .and. File_Present_Flag(Chan_Idx)) then
          call READ_ABI_MAX_FOCAL_PLANE_TEMP(Ncid_Level1b(Chan_Idx), &
                                            'maximum_focal_plane_temperature', &
                                            ch(Clavrx_Chan_Map(Chan_Idx))%Max_Focal_Plane_Temp)
        end if
      end do
    endif

    !---  determine bounds of data
    call DETERMINE_BOUNDS_STATIC_NAV(trim(Static_Nav_Full_File),Nav%Limit_Flag, &
                                     Nav%Lat_South_Limit, Nav%Lon_West_Limit,   &
                                     Nav%Lat_North_Limit, Nav%Lon_East_Limit, Image%X_Stride, Image%Y_Stride, &
                                     Image%Number_Of_Elements,Image%Number_Of_Lines)

!   call DETERMINE_BOUNDS_STATIC_NAV_ORIGINAL(trim(Static_Nav_Full_File),Nav%Limit_Flag, &
!                                    Nav%Lat_South_Limit, Nav%Lon_West_Limit,   &
!                                    Nav%Lat_North_Limit, Nav%Lon_East_Limit, Image%X_Stride, Image%Y_Stride, &
!                                    Image%Number_Of_Elements,Image%Number_Of_Lines)

  end subroutine SETUP_READ_LEVEL1B_FIXED_GRID_STATIC_NAV

  !------------------------------------------------------------------------------------
  !  Deallocate arrays 
  !------------------------------------------------------------------------------------
  subroutine DESTROY_READ_LEVEL1B_FIXED_GRID_STATIC_NAV()

    integer:: Chan_Idx

    !--- close netcdf files
    call CLOSE_NETCDF(ncid_static)

    ! Open Level1B NetCDF files. - WCS3
    IF (Image%Nc_Format_flag) THEN
        do Chan_Idx = 1, Num_Chan_Sensor
            if (Chan_On_Sensor(Chan_Idx) == 1) then
                call CLOSE_NETCDF(Ncid_Level1b(Chan_Idx))
            end if
        end do
    endif

    !--- deallocate
    if (allocated(CLAVRx_Chan_Map)) deallocate(CLAVRx_Chan_Map)
    if (allocated(Variable_Name)) deallocate(Variable_Name)
    if (allocated(Chan_Stride)) deallocate(Chan_Stride)
    if (allocated(Chan_On_Sensor)) deallocate(Chan_On_Sensor)
    if (allocated(Ncid_Level1b)) deallocate(Ncid_Level1b)
    if (allocated(Level1b_File)) deallocate(Level1b_File)

    if (allocated(DQF_Seg_2km)) deallocate(DQF_Seg_2km)
    if (allocated(Output_Seg_2km)) deallocate(Output_Seg_2km)
    if (allocated(Output_Seg_1km)) deallocate(Output_Seg_1km)
    if (allocated(Output_Seg_500m)) deallocate(Output_Seg_500m)
    if (allocated(Min_Output_Seg_2km)) deallocate(Min_Output_Seg_2km)
    if (allocated(Max_Output_Seg_2km)) deallocate(Max_Output_Seg_2km)
    if (allocated(Stddev_Output_Seg_2km)) deallocate(Stddev_Output_Seg_2km)
    if (allocated(DQF_Seg_500m)) deallocate(DQF_Seg_500m)
    if (allocated(Min_Output_Seg_500m)) deallocate(Min_Output_Seg_500m)
    if (allocated(Max_Output_Seg_500m)) deallocate(Max_Output_Seg_500m)
    if (allocated(Stddev_Output_Seg_500m)) deallocate(Stddev_Output_Seg_500m)
    if (allocated(DQF_Seg_1km)) deallocate(DQF_Seg_1km)
    if (allocated(Min_Output_Seg_1km)) deallocate(Min_Output_Seg_1km)
    if (allocated(Max_Output_Seg_1km)) deallocate(Max_Output_Seg_1km)
    if (allocated(Stddev_Output_Seg_1km)) deallocate(Stddev_Output_Seg_1km)
    if (allocated(Output_Time)) deallocate(Output_Time)
    if (allocated(Rad_to_Ref_Fac)) deallocate(Rad_to_Ref_Fac)
    if (allocated(File_Present_Flag)) deallocate(File_Present_Flag)

    !--- close netcdf files
    !call close_netcdf(ncid_static)

  end subroutine DESTROY_READ_LEVEL1B_FIXED_GRID_STATIC_NAV
  
  !------------------------------------------------------------------------------------
  !  Public Read Routine
  !------------------------------------------------------------------------------------
  subroutine READ_LEVEL1B_FIXED_GRID_STATIC_NAV()

    integer:: Chan_Idx, Elem_Idx, Line_Idx
    integer:: Chan_Clavrx_Idx
    logical:: Read_Time !  (true if time is to be read, set to false after first read)
    real :: kappa0
    real :: Factor
    real :: Image_Date
    integer :: i_line, add_lines
    integer:: Alloc_Status
    character(len=50):: var_name

    !--- get lat, lon, zen, sataz from the static nav file
    call READ_SEGMENT_STATIC_NAV(Ncid_Static)

    !--- get solar zenith angle 
    do Elem_Idx = 1, Image%Number_Of_Elements
       do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment
          call POSSOL(int(Image%Start_Doy,kind=int4),Image%Mean_Time_Hours,Nav%Lon(Elem_Idx,Line_Idx), &
                      Nav%Lat(Elem_Idx,Line_Idx),Geo%Solzen(Elem_Idx,Line_Idx),Geo%Solaz(Elem_Idx,Line_Idx)) 
       end do
    end do

    !--- relative azimuth
    Geo%Relaz = relative_azimuth(Geo%Solaz, Geo%Sataz)

    !--- scattering angle
    Geo%Scatangle = scattering_angle(Geo%Solzen, Geo%Satzen, Geo%Relaz)

    !--- glint angle
    Geo%Glintzen = glint_angle(Geo%Solzen, Geo%Satzen, Geo%Relaz)

    !--- ensure missing values
    where((Geo%Solzen .eqr. MISSING_VALUE_REAL4) .or. &
          (Geo%Solaz  .eqr. MISSING_VALUE_REAL4) .or. &
          (Geo%Sataz  .eqr. MISSING_VALUE_REAL4) .or. &
          (Geo%Satzen .eqr. MISSING_VALUE_REAL4))

          Geo%Relaz = MISSING_VALUE_REAL4
          Geo%Scatangle = MISSING_VALUE_REAL4
          Geo%Glintzen = MISSING_VALUE_REAL4

    endwhere

    !--- allocate memory for the temp arrays for reading from level-1b
    if (.not. allocated(Output_Seg_2km)) allocate(Output_Seg_2km(Image%Number_Of_Elements, Image%Number_Of_Lines_Per_Segment))
    if (.not. allocated(Output_Seg_1km)) then
        allocate(Output_Seg_1km(Image%Number_Of_Elements*2,Image%Number_Of_Lines_Per_Segment*2,5),stat=Alloc_Status)
        if (Alloc_Status /= 0) then
            call MESG( "Could not allocate Output_Seg_1km", level = verb_lev % ERROR )
            stop 5
        endif
    endif
    if (.not. allocated(Output_Seg_500m)) then
        allocate(Output_Seg_500m(Image%Number_Of_Elements*4, Image%Number_Of_Lines_Per_Segment*4,1),stat=Alloc_Status)
        if (Alloc_Status /= 0) then
            call MESG( "Could not allocate Output_Seg_500m", level = verb_lev % ERROR )
            stop 5
        endif
    endif
    if (.not. allocated(DQF_Seg_2km)) allocate(DQF_Seg_2km(Image%Number_Of_Elements, Image%Number_Of_Lines_Per_Segment))
    if (.not. allocated(Min_Output_Seg_2km)) allocate(Min_Output_Seg_2km(Image%Number_Of_Elements, Image%Number_Of_Lines_Per_Segment))
    if (.not. allocated(Max_Output_Seg_2km)) allocate(Max_Output_Seg_2km(Image%Number_Of_Elements, Image%Number_Of_Lines_Per_Segment))
    if (.not. allocated(Stddev_Output_Seg_2km)) allocate(Stddev_Output_Seg_2km(Image%Number_Of_Elements, Image%Number_Of_Lines_Per_Segment))
    if (.not. allocated(Output_Time)) allocate(Output_Time(Image%Number_Of_Lines_Per_Segment))
    
    add_lines = (Image%Segment_Number-1) * Image%Number_Of_Lines_Per_Segment
    Image%Scan_Number = [(i_line, i_line = Line_start + add_lines, &
                          Line_start + Image%Number_Of_Lines_Per_Segment + add_lines - 1, 1)]

    !--- Read level1b data
    Read_Time = .true.
    chan_loop: do Chan_Idx = 1, Num_Chan_Sensor

     if (Sensor%Chan_On_Flag_Default(Clavrx_Chan_Map(Chan_Idx)) == 1) then

       if (File_Present_Flag(Chan_Idx)) then

            if(Chan_Stride(Chan_Idx) == 2) then
                call READ_SEGMENT_LEVEL1B_RAW(Ncid_Level1b(Chan_Idx),Variable_Name(Chan_Idx), &
                                       Image%Segment_Number, Image%Number_Of_Lines_Per_Segment, &
                                       Chan_Stride(Chan_Idx), X_Sample_Offset, Y_Sample_Offset, Output_Seg_1km(:,:,Chan_Idx))
            else if(Chan_Stride(Chan_Idx) == 4) then
                call READ_SEGMENT_LEVEL1B_RAW(Ncid_Level1b(Chan_Idx),Variable_Name(Chan_Idx), &
                                       Image%Segment_Number, Image%Number_Of_Lines_Per_Segment, &
                                       Chan_Stride(Chan_Idx), X_Sample_Offset, Y_Sample_Offset, Output_Seg_500m(:,:,1))
                !if (.not. allocated(DQF_Seg_500m)) allocate(DQF_Seg_500m(Image%Number_Of_Elements*4, Image%Number_Of_Lines_Per_Segment*4))
                !if (.not. allocated(Min_Output_Seg_500m)) allocate(Min_Output_Seg_500m(Image%Number_Of_Elements*4, Image%Number_Of_Lines_Per_Segment*4))
                !if (.not. allocated(Max_Output_Seg_500m)) allocate(Max_Output_Seg_500m(Image%Number_Of_Elements*4, Image%Number_Of_Lines_Per_Segment*4))
                !if (.not. allocated(Stddev_Output_Seg_500m)) allocate(Stddev_Output_Seg_500m(Image%Number_Of_Elements*4, Image%Number_Of_Lines_Per_Segment*4))
                !call READ_SEGMENT_LEVEL1B_VER2(Ncid_Level1b(Chan_Idx),Variable_Name(Chan_Idx), &
                                       !Image%Segment_Number, Image%Number_Of_Lines_Per_Segment, Image%Number_Of_Segments, &
                                       !Image%X_Stride, Image%Y_Stride, 1, X_Sample_Offset, Y_Sample_Offset, &
                                       !0,  &
                                       !Output_Seg_500m(:,:,1), Min_Output_Seg_500m, Max_Output_Seg_500m, Stddev_Output_Seg_500m, Read_Time, DQF_Seg_500m)
            endif

            !--- read channel values if file is present
            call READ_SEGMENT_LEVEL1B_VER2(Ncid_Level1b(Chan_Idx),Variable_Name(Chan_Idx), &
                                   Image%Segment_Number, Image%Number_Of_Lines_Per_Segment, Image%Number_Of_Segments, &
                                   Image%X_Stride, Image%Y_Stride, Chan_Stride(Chan_Idx), X_Sample_Offset, Y_Sample_Offset, &
                                   Image%Chan_Average_Flag,  &
                                   Output_Seg_2km, Min_Output_Seg_2km, Max_Output_Seg_2km, Stddev_Output_Seg_2km, Read_Time, DQF_Seg_2km)

           !--- set SOURCE Flag if ABI Fusion used for any channel
           ch(Clavrx_Chan_Map(Chan_Idx))%Source = 0    !note source of this channel is from imager
           if (Sensor%WMO_Id == 270 .or. Sensor%WMO_Id == 271 .or. Sensor%WMO_Id == 272) then
             if (index(Level1b_File(Chan_Idx),'FR_') > 0) then
                ch(Clavrx_Chan_Map(Chan_Idx))%Source = 1    !note source of this channel is from sounder
             endif
           endif

       else
            !--- fill channel values with missing if file is missing and channel is on
            Output_Seg_2km = MISSING_VALUE_REAL4
            Min_Output_Seg_2km = MISSING_VALUE_REAL4
            Max_Output_Seg_2km = MISSING_VALUE_REAL4
            Stddev_Output_Seg_2km = MISSING_VALUE_REAL4
            DQF_Seg_2km = Missing_Value_Int1
       endif


       Chan_Clavrx_Idx = Clavrx_Chan_Map(Chan_Idx)

       !--- store DQF in clavrx global data structure
       ch(Chan_Clavrx_Idx)%DQF = DQF_Seg_2km

       !--- ahi radiances for all channels have to be converted from nasa units (W m-2 sr-1 um-1) to noaa units (mW m-2 sr-1 (cm-1)-1)
       if (Sensor%WMO_Id == 173 .or. Sensor%WMO_Id == 174) then
          if (ch(Chan_Clavrx_Idx)%Obs_Type == THERMAL_OBS_TYPE .or. &
               ch(Chan_Clavrx_Idx)%Obs_Type == MIXED_OBS_TYPE) then
            call CONVERT_RADIANCE(Output_Seg_2km,Planck_Nu(Chan_Clavrx_Idx),MISSING_VALUE_REAL4)
            if (maxval(Min_Output_Seg_2km) /= MISSING_VALUE_REAL4) then
              call CONVERT_RADIANCE(Min_Output_Seg_2km,Planck_Nu(Chan_Clavrx_Idx),MISSING_VALUE_REAL4)
            end if
            if (maxval(Max_Output_Seg_2km) /= MISSING_VALUE_REAL4) then
              call CONVERT_RADIANCE(Max_Output_Seg_2km,Planck_Nu(Chan_Clavrx_Idx),MISSING_VALUE_REAL4)
            end if
          end if
       end if

       !--- Compute Bt from Radiances
       if (ch(Chan_Clavrx_Idx)%Obs_Type == THERMAL_OBS_TYPE .or. &
             ch(Chan_Clavrx_Idx)%Obs_Type == MIXED_OBS_TYPE) then

             ch(Chan_Clavrx_Idx)%Rad_Toa = Output_Seg_2km

             do Elem_Idx = 1, Image%Number_Of_Elements
                do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment
                   ch(Chan_Clavrx_Idx)%Bt_Toa(Elem_Idx,Line_Idx) = &
                   PLANCK_TEMP_FAST(Chan_Clavrx_Idx, ch(Chan_Clavrx_Idx)%Rad_Toa(Elem_Idx,Line_Idx))
                enddo
             enddo
       endif

       !--- make reflectances (not normalized yet) aka scaled radiances
       if (ch(Chan_Clavrx_Idx)%Obs_Type == SOLAR_OBS_TYPE) then

         !--- For ABI, use kappa0 factor from the NetCDF file.
         if (Sensor%WMO_Id == 270 .or. Sensor%WMO_Id == 271 .or. Sensor%WMO_Id == 272) then

           Factor = 1.0

           !--- perform Band-2 Correction (Band 2 = Ch 1 in Clavrx)
           if (Chan_Clavrx_Idx == 1) then 
             Image_Date = Image% time_start % Year + (Image % time_start % dayOfYear - 1.0) / 365.25
             if (Image_Date < Band2_Correction_Start_Date) then
               Factor = Band2_Correction_Factor
             endif
           endif

           if (ch(Chan_Clavrx_Idx)%Obs_Type == SOLAR_OBS_TYPE) then
             call read_abi_kappa0(Ncid_Level1b(Chan_Idx), kappa0)
             ch(Chan_Clavrx_Idx)%Ref_Toa = Factor * Output_Seg_2km * kappa0 * 100.0
             if(Chan_Idx .le. 6) Rad_To_Ref_Fac(Chan_Idx) = Factor * kappa0 * 100.0
             where(Output_Seg_2km == MISSING_VALUE_REAL4)
                  ch(Chan_Clavrx_Idx)%Ref_Toa = MISSING_VALUE_REAL4
             endwhere
             if (allocated(ch(Chan_Clavrx_Idx)%Ref_Toa_Min_Sub)) then
                ch(Chan_Clavrx_Idx)%Ref_Toa_Min_Sub = Factor * Min_Output_Seg_2km * kappa0 * 100.0
                where(Min_Output_Seg_2km == MISSING_VALUE_REAL4)
                      ch(Chan_Clavrx_Idx)%Ref_Toa_Min_Sub = MISSING_VALUE_REAL4
                endwhere
             endif
             if (allocated(ch(Chan_Clavrx_Idx)%Ref_Toa_Max_Sub)) then
                ch(Chan_Clavrx_Idx)%Ref_Toa_Max_Sub = Factor * Max_Output_Seg_2km * kappa0 * 100.0
                where(Max_Output_Seg_2km == MISSING_VALUE_REAL4)
                      ch(Chan_Clavrx_Idx)%Ref_Toa_Max_Sub = MISSING_VALUE_REAL4
                endwhere
             endif
             if (allocated(ch(Chan_Clavrx_Idx)%Ref_Toa_Std_Sub)) then
                ch(Chan_Clavrx_Idx)%Ref_Toa_Std_Sub = Factor * Stddev_Output_Seg_2km * kappa0 * 100.0
                where(Stddev_Output_Seg_2km == MISSING_VALUE_REAL4)
                      ch(Chan_Clavrx_Idx)%Ref_Toa_Std_Sub = MISSING_VALUE_REAL4
                endwhere
             endif
           endif

         else

           ch(Chan_Clavrx_Idx)%Ref_Toa = Output_Seg_2km * Rad_To_Ref_Fac(Chan_Idx)
           where(Output_Seg_2km == MISSING_VALUE_REAL4)
               ch(Chan_Clavrx_Idx)%Ref_Toa = MISSING_VALUE_REAL4
           endwhere
           if (allocated(ch(Chan_Clavrx_Idx)%Ref_Toa_Min_Sub)) then
              ch(Chan_Clavrx_Idx)%Ref_Toa_Min_Sub = Min_Output_Seg_2km * Rad_To_Ref_Fac(Chan_Idx)
              where(Min_Output_Seg_2km == MISSING_VALUE_REAL4)
                 ch(Chan_Clavrx_Idx)%Ref_Toa_Min_Sub = MISSING_VALUE_REAL4
              endwhere
           endif
           if (allocated(ch(Chan_Clavrx_Idx)%Ref_Toa_Max_Sub)) then
               ch(Chan_Clavrx_Idx)%Ref_Toa_Max_Sub = Max_Output_Seg_2km * Rad_To_Ref_Fac(Chan_Idx)
               where(Max_Output_Seg_2km == MISSING_VALUE_REAL4)
                    ch(Chan_Clavrx_Idx)%Ref_Toa_Max_Sub = MISSING_VALUE_REAL4
               endwhere
           endif
           if (allocated(ch(Chan_Clavrx_Idx)%Ref_Toa_Std_Sub)) then
             ch(Chan_Clavrx_Idx)%Ref_Toa_Std_Sub = Stddev_Output_Seg_2km * Rad_To_Ref_Fac(Chan_Idx)
             where(Stddev_Output_Seg_2km == MISSING_VALUE_REAL4)
                 ch(Chan_Clavrx_Idx)%Ref_Toa_Std_Sub = MISSING_VALUE_REAL4
             endwhere
           endif

         endif !sensor check
 
      endif !solar obs check
 
    endif   !chan on check

   enddo chan_loop   

   ! --- read dark composite (ONLY ABI GOES-16 FOR NOW)
   if (Static_Dark_Sky_Flag) then
     if (Image%Segment_Number == 1) then
       call FIND_DARK_COMPOSITE()
     endif
     if (Dark_File_Exist) then
      call READ_DARK_COMPOSITE(Image%X_Stride, Image%Y_Stride, &
                               Static_Ref_065um_Dark_Composite, &
                               Static_Ref_065um_Dark_Composite_Stddev)
     endif
   endif


!  !--- Read QRNN CTP
!  QRNN_CTP_FLAG = .true.
!  QRNN_CTP_Data_Dir = '/apollo/cloud/scratch/AMV_BUST/ANDY/chuck/'
   
!  if (QRNN_CTP_FLAG) then
!    if (Image%Segment_Number == 1) then
!      call FIND_QRNN_CTP()
!    endif
!    if (QRNN_CTP_File_Exist) then
!      call READ_QRNN_CTP(Image%X_Stride, Image%Y_Stride,QRNN_CTP,QRNN_CTP_UNCER)
!    endif
!  endif

end subroutine READ_LEVEL1B_FIXED_GRID_STATIC_NAV

!-------------------------------------------------------------------------------
!  Public Read Routine for HCAST and HSD  - WCS3
! Image%Chan_Average_Flag = 0 = sample
!                           1 = average
!                           2 = average and compute min and max
! Currently Image%Chan_Average_Flag not being used. Will be implmented once 
!           Sameple code available. - 8 April 2019
!-------------------------------------------------------------------------------
  subroutine READ_HSD_FIXED_GRID_STATIC_NAV()

    integer:: Chan_Idx, Elem_Idx, Line_Idx
    integer:: Chan_Clavrx_Idx
    integer:: Status
    INTEGER(KIND=INT4) :: AHI_FILE_TYPE ! Needed to distinguish how to read data in - WCS3
    logical:: Read_Time !  (true if time is to be read, set to false after first read)
    real :: kappa0
    integer :: i_line, add_lines
    INTEGER(KIND=INT4), DIMENSION(2) :: Start
    INTEGER(KIND=INT4), DIMENSION(2) :: Stride
    INTEGER(KIND=INT4), DIMENSION(2) :: Edge
    CHARACTER(LEN=29) :: Base_Filename
    
    ! Need this stuff to get 2km information. Static_Nav conversions done in reader
    integer:: Static_Nav_Elem_Start_Segment, Static_Nav_Elem_End_Segment, Static_Nav_Elem_Count_Segment
    integer:: Static_Nav_Line_Start_Segment, Static_Nav_Line_End_Segment, Static_Nav_Line_Count_Segment
    integer:: Static_Nav_Elem_End_File, Static_Nav_Line_End_File
    integer:: Elem_Start_Segment, Elem_End_Segment, Elem_Count_Segment
    integer:: Line_Start_Segment, Line_End_Segment, Line_Count_Segment
    integer:: Alloc_Status

    

    !--- get lat, lon, zen, sataz from the static nav file
    call READ_SEGMENT_STATIC_NAV(ncid_static)

    !--- get solar zenith angle 
    do Elem_Idx = 1, Image%Number_Of_Elements
       do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment
          call POSSOL(int(Image%Start_Doy,kind=int4),Image%Mean_Time_Hours,Nav%Lon(Elem_Idx,Line_Idx), &
                      Nav%Lat(Elem_Idx,Line_Idx),Geo%Solzen(Elem_Idx,Line_Idx),Geo%Solaz(Elem_Idx,Line_Idx)) 
       end do
    end do

    !--- relative azimuth
    Geo%Relaz = relative_azimuth(Geo%Solaz, Geo%Sataz)

    !--- scattering angle
    Geo%Scatangle = scattering_angle(Geo%Solzen, Geo%Satzen, Geo%Relaz)

    !--- glint angle
    Geo%Glintzen = glint_angle(Geo%Solzen, Geo%Satzen, Geo%Relaz)

    !--- ensure missing values
    where((Geo%Solzen .eqr. MISSING_VALUE_REAL4) .or. &
          (Geo%Solaz  .eqr. MISSING_VALUE_REAL4) .or. &
          (Geo%Sataz  .eqr. MISSING_VALUE_REAL4) .or. &
          (Geo%Satzen .eqr. MISSING_VALUE_REAL4))

          Geo%Relaz = MISSING_VALUE_REAL4
          Geo%Scatangle = MISSING_VALUE_REAL4
          Geo%Glintzen = MISSING_VALUE_REAL4

    endwhere


    !--- allocate memory for the temp arrays for reading from level-1b
    if (.not. allocated(Output_Seg_2km)) then 
        allocate(Output_Seg_2km(Image%Number_Of_Elements, Image%Number_Of_Lines_Per_Segment),stat=Alloc_Status)
        if (Alloc_Status /= 0) then
            call MESG( "Could not allocate Output_Seg_2km", level = verb_lev % ERROR )
            stop 5
        endif
    endif
    if (.not. allocated(Output_Seg_1km)) then
        allocate(Output_Seg_1km(Image%Number_Of_Elements*2,Image%Number_Of_Lines_Per_Segment*2,4),stat=Alloc_Status)
        if (Alloc_Status /= 0) then
            call MESG( "Could not allocate Output_Seg_1km", level = verb_lev % ERROR )
            stop 5
        endif
    endif
    if (.not. allocated(Output_Seg_500m)) then
        allocate(Output_Seg_500m(Image%Number_Of_Elements*4, Image%Number_Of_Lines_Per_Segment*4,1),stat=Alloc_Status)
        if (Alloc_Status /= 0) then
            call MESG( "Could not allocate Output_Seg_500m", level = verb_lev % ERROR )
            stop 5
        endif
    endif
    if (.not. allocated(Min_Output_Seg_2km)) allocate(Min_Output_Seg_2km(Image%Number_Of_Elements, Image%Number_Of_Lines_Per_Segment))
    if (.not. allocated(Max_Output_Seg_2km)) allocate(Max_Output_Seg_2km(Image%Number_Of_Elements, Image%Number_Of_Lines_Per_Segment))
    if (.not. allocated(Output_Time)) allocate(Output_Time(Image%Number_Of_Lines_Per_Segment))

    
    add_lines = (Image%Segment_Number-1) * Image%Number_Of_Lines_Per_Segment
    Image%Scan_Number = [(i_line, i_line = Line_start + add_lines, &
                          Line_start + Image%Number_Of_Lines_Per_Segment + add_lines - 1, 1)]
                          
    !---- Start/edge/Stride in Static_Nav space will be taken into account within LibHim_Raddata_Call
    !--- for now start with 
   
   
   !--- determine Static_Nav pixel sizes 
   Static_Nav_Elem_End_File = Num_Elem_Fd 
   Static_Nav_Line_End_File = Num_Line_Fd 
    
    !--- Static Nav  element terms
   Static_Nav_Elem_Start_Segment = (Elem_Start - 1) * 1 + 1

   !--- Static Nav  line terms.
   !   note that this is similar, but not the same as how the netCDF reads do things
   !   this is because we only need the start line, nothing more.
   Static_Nav_Line_Start_Segment = (Line_Start-1)*1 + (Image%Segment_Number-1)*Image%Number_Of_Lines_Per_Segment*Image%Y_Stride + 1
                               
     
                    
    !---- initially, stride = 1 in options file
    Stride(1) = Image%X_Stride
    Stride(2) = Image%Y_Stride
    
    !--- these need to be absolute index from Static nav 

   ! Start line/element is all that is needed, libhimawari takes care of the rest
   
   Start(1) = Static_Nav_Elem_Start_Segment
   Start(2) = Static_Nav_Line_Start_Segment
   
   ! Initialize output arrays to missing
   Output_Seg_2km = MISSING_VALUE_REAL4
   Output_Seg_1km = MISSING_VALUE_REAL4
   Output_Seg_500m = MISSING_VALUE_REAL4
   Min_Output_Seg_2km = MISSING_VALUE_REAL4
   Max_Output_Seg_2km = MISSING_VALUE_REAL4
   
    

    !LibHimawari reader expects HS_H08_20190312_2140_B01_FLDK_, since it uses wildcards
    Base_Filename = trim(Image%Level1b_Name(1:29))
    do Chan_Idx = 1, Num_Chan_Sensor

      if (Chan_On_Sensor(Chan_Idx) == 1) then
            ! Note - this is the only thing that needs to have ifdef around it
#ifdef LIBHIM
            !Default libHim type is HSD - WCS2
            AHI_FILE_TYPE = FILE_TYPE_AHI_HSF
            !if the file is HCAST, then change types - WCS3
            IF (Image%DB_Flag) AHI_FILE_TYPE = FILE_TYPE_AHI_HRIT

            if (Chan_Idx < 5) then
                if(Chan_Idx == 3) then
                    Status = Get_Him_Raddata(AHI_FILE_TYPE,  & !dynamic reading of file type
                                             trim(Image%Level1b_Path), &
                                             trim(Base_Filename),          &
                                             Chan_Idx,                &
                                             Start,                  &
                                             Stride,                 &
                                             Image%Number_Of_Lines_Read_This_Segment, & ! this is so that you don't over reach in the data
                                             3, &
                                             Output_Seg_500m(:,:,1), &
                                             !Min_Output_Seg and Max_Output_Seg are only used for Average mode 2
                                             ! but it is easier (and cleaner) to have a single call rather than
                                             ! multiple calls
                                             Min_Output_Seg_2km, & 
                                             Max_Output_Seg_2km)
                else
                    Status = Get_Him_Raddata(AHI_FILE_TYPE,  & !dynamic reading of file type
                                             trim(Image%Level1b_Path), &
                                             trim(Base_Filename),          &
                                             Chan_Idx,                &
                                             Start,                  &
                                             Stride,                 &
                                             Image%Number_Of_Lines_Read_This_Segment, & ! this is so that you don't over reach in the data
                                             3, &
                                             Output_Seg_1km(:,:,Chan_Idx), &
                                             !Min_Output_Seg and Max_Output_Seg are only used for Average mode 2
                                             ! but it is easier (and cleaner) to have a single call rather than
                                             ! multiple calls
                                             Min_Output_Seg_2km, & 
                                             Max_Output_Seg_2km)
                endif
            endif 

            Status = Get_Him_Raddata(AHI_FILE_TYPE,  & !dynamic reading of file type
                                     trim(Image%Level1b_Path), &
                                     trim(Base_Filename),          &
                                     Chan_Idx,                &
                                     Start,                  &
                                     Stride,                 &
                                     Image%Number_Of_Lines_Read_This_Segment, & ! this is so that you don't over reach in the data
                                     Image%Chan_Average_Flag, &
                                     Output_Seg_2km, &
                                     !Min_Output_Seg and Max_Output_Seg are only used for Average mode 2
                                     ! but it is easier (and cleaner) to have a single call rather than
                                     ! multiple calls
                                     Min_Output_Seg_2km, & 
                                     Max_Output_Seg_2km)
                                     
#endif

       Chan_Clavrx_Idx = Clavrx_Chan_Map(Chan_Idx)
       

       !--- store DQF in clavrx global data structure
       !--- set to good for now, since there are no DQFs in AHI
       ch(Chan_Clavrx_Idx)%DQF = 0

       !--- ahi radiances for all channels have to be converted from nasa units (W m-2 sr-1 um-1) to noaa units (mW m-2 sr-1 (cm-1)-1)
       if (ch(Chan_Clavrx_Idx)%Obs_Type == THERMAL_OBS_TYPE .or. &
               ch(Chan_Clavrx_Idx)%Obs_Type == MIXED_OBS_TYPE) then
            call CONVERT_RADIANCE(Output_Seg_2km,Planck_Nu(Chan_Clavrx_Idx),MISSING_VALUE_REAL4)
            if (maxval(Min_Output_Seg_2km) /= MISSING_VALUE_REAL4) then
              call CONVERT_RADIANCE(Min_Output_Seg_2km,Planck_Nu(Chan_Clavrx_Idx),MISSING_VALUE_REAL4)
            end if
            if (maxval(Max_Output_Seg_2km) /= MISSING_VALUE_REAL4) then
              call CONVERT_RADIANCE(Max_Output_Seg_2km,Planck_Nu(Chan_Clavrx_Idx),MISSING_VALUE_REAL4)
            end if
       end if

       !--- Bt
       if (ch(Chan_Clavrx_Idx)%Obs_Type == THERMAL_OBS_TYPE .or. &
             ch(Chan_Clavrx_Idx)%Obs_Type == MIXED_OBS_TYPE) then

             ch(Chan_Clavrx_Idx)%Rad_Toa = Output_Seg_2km

             do Elem_Idx = 1, Image%Number_Of_Elements
                do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment

                   ch(Chan_Clavrx_Idx)%Bt_Toa(Elem_Idx,Line_Idx) = &
                   PLANCK_TEMP_FAST(Chan_Clavrx_Idx, ch(Chan_Clavrx_Idx)%Rad_Toa(Elem_Idx,Line_Idx))

                enddo
             enddo
       endif

       !--- make reflectances (not normalized yet) aka scaled radiances
       if (ch(Chan_Clavrx_Idx)%Obs_Type == SOLAR_OBS_TYPE) then

             ch(Chan_Clavrx_Idx)%Ref_Toa = Output_Seg_2km * Rad_To_Ref_Fac(Chan_Idx)
             where(Output_Seg_2km == MISSING_VALUE_REAL4)
               ch(Chan_Clavrx_Idx)%Ref_Toa = MISSING_VALUE_REAL4
             endwhere
             if (allocated(ch(Chan_Clavrx_Idx)%Ref_Toa_Min_Sub)) then
                ch(Chan_Clavrx_Idx)%Ref_Toa_Min_Sub = Min_Output_Seg_2km * Rad_To_Ref_Fac(Chan_Idx)
                where(Min_Output_Seg_2km == MISSING_VALUE_REAL4)
                      ch(Chan_Clavrx_Idx)%Ref_Toa_Min_Sub = MISSING_VALUE_REAL4
                endwhere
             endif
             if (allocated(ch(Chan_Clavrx_Idx)%Ref_Toa_Max_Sub)) then
                ch(Chan_Clavrx_Idx)%Ref_Toa_Max_Sub = Max_Output_Seg_2km * Rad_To_Ref_Fac(Chan_Idx)
                where(Max_Output_Seg_2km == MISSING_VALUE_REAL4)
                      ch(Chan_Clavrx_Idx)%Ref_Toa_Max_Sub = MISSING_VALUE_REAL4
                endwhere
             endif
           endif
       endif
    enddo
   
end subroutine READ_HSD_FIXED_GRID_STATIC_NAV

subroutine READ_SEGMENT_LEVEL1B_RAW(Ncid,Chan_Name,Segment_Number, Number_Of_Lines_Per_Segment, &
                                Chan_Scale, X_Res_Offset, Y_Res_Offset, &
                                Output_Seg)

   character(len=*), intent(in):: Chan_Name
   integer, intent(in):: Ncid,Segment_Number, Number_Of_Lines_Per_Segment, &
                         X_Res_Offset, Y_Res_Offset, Chan_Scale

   real, dimension(:,:), intent(out):: Output_Seg
   integer:: nx, ny

   integer:: Native_X_Stride, Native_Y_Stride
   integer:: X_Stride, Y_Stride
   integer:: Native_Elem_Start_Segment, Native_Elem_End_Segment, Native_Elem_Count_Segment
   integer:: Native_Line_Start_Segment, Native_Line_End_Segment, Native_Line_Count_Segment
   integer:: Native_Elem_End_File, Native_Line_End_File
   integer:: Elem_Start_Segment, Elem_End_Segment, Elem_Count_Segment
   integer:: Line_Start_Segment, Line_End_Segment, Line_Count_Segment
   integer:: Elem_End_File, Line_End_File
   integer:: Y_Res_Stride
   integer:: i, j, ni, nj, i1, i2, j1, j2, is, js, n, nn

   integer, dimension(2):: start_2d, stride_2d, count_2d
   
   !--- initialize
   Output_Seg = MISSING_VALUE_REAL4

   !--- save output size
   nx = size(Output_Seg,1)
   ny = size(Output_Seg,2)

   !--- determine native pixel sizes 
   Native_Elem_End_File = Num_Elem_Fd*Chan_Scale
   Native_Line_End_File = Num_Line_Fd*Chan_Scale

   !--- there is no striding in the reading of the native data
   Native_X_Stride = 1
   Native_Y_Stride = 1

   !--- effect stride the product of chosen thinning stride and native
   !--- resolution stride value
   X_Stride = 1
   Y_Stride = 1

   !--- native element terms
   Native_Elem_Start_Segment = (Elem_Start-1)*Chan_Scale+1
   Native_Elem_Count_Segment = Num_Elem*Chan_Scale
   Native_Elem_End_Segment = Native_Elem_Start_Segment + Native_Elem_Count_Segment - 1
 
   !--- native line terms
   Native_Line_Start_Segment = (Line_Start-1)*Chan_Scale + (Segment_Number-1)*Number_Of_Lines_Per_Segment*Chan_Scale + 1
   Native_Line_Count_Segment = Number_Of_Lines_Per_Segment*Chan_Scale
   Native_Line_End_Segment = Native_Line_Start_Segment + Native_Line_Count_Segment - 1
   Native_Line_End_Segment =   min(Native_Line_End_Segment,Native_Line_End_File)

   Native_Line_Count_Segment = (Native_Line_End_Segment - Native_Line_Start_Segment) + 1

   !--- read native resolution data
   start_2d =  (/Native_Elem_Start_Segment, Native_Line_Start_Segment/)
   stride_2d = (/Native_X_Stride,Native_Y_Stride/)
   count_2d =   (/Native_Elem_Count_Segment, Native_Line_Count_Segment/)
   call READ_AND_UNSCALE_NETCDF_2d(Ncid, Start_2d, Stride_2d, Count_2d, Chan_Name, Output_Seg)

end subroutine READ_SEGMENT_LEVEL1B_RAW


!------------------------------------------------------------------------------------
!
! VER2 reads in the native data and samples and strides using do loops
!
! X_Thin_Stride - the stride chosen by the user for the nominal resolution (data thinning)
! X_Res_Stride - the stride due to the channel resolution being finer than nominal value
! X_Stride - the final stride of the data to be read in
! X_Res_Offset - the offset in the fine data resolution - controls which pixel
!                is selected in the sample
! Image%Chan_Average_Flag = 0 = sample
!                    1 = average
!                    2 = average and compute min and max
!------------------------------------------------------------------------------------
subroutine READ_SEGMENT_LEVEL1B_VER2(Ncid,Chan_Name,Segment_Number, Number_Of_Lines_Per_Segment, &
                                Number_Of_Segments, X_Thin_Stride, Y_Thin_Stride, &
                                X_Res_Stride, X_Res_Offset_In, Y_Res_Offset_In, Chan_Average_Flag, &
                                Output_Seg, Min_Output_Seg, Max_Output_Seg, Stddev_Output_Seg, Read_Time, DQF_Seg)

   character(len=*), intent(in):: Chan_Name
   integer, intent(in):: Ncid,Segment_Number, Number_Of_Lines_Per_Segment, &
                         Number_Of_Segments, X_Thin_Stride, Y_Thin_Stride,  &
                         X_Res_Stride, X_Res_Offset_In, Y_Res_Offset_In, Chan_Average_Flag
   logical, intent(inout):: Read_Time

   real, dimension(:,:), intent(out):: Output_Seg, Min_Output_Seg, Max_Output_Seg, Stddev_Output_Seg
   integer (kind=int1), dimension(:,:), intent(out):: DQF_Seg

   real, dimension(:,:), allocatable:: Native_Output_Seg
   integer(kind=int1), dimension(:,:), allocatable:: Native_DQF_Seg
   integer:: nx, ny

   integer:: Native_X_Stride, Native_Y_Stride
   integer:: X_Stride, Y_Stride
   integer:: Native_Elem_Start_Segment, Native_Elem_End_Segment, Native_Elem_Count_Segment
   integer:: Native_Line_Start_Segment, Native_Line_End_Segment, Native_Line_Count_Segment
   integer:: Native_Elem_End_File, Native_Line_End_File
   integer:: Elem_Start_Segment, Elem_End_Segment, Elem_Count_Segment
   integer:: Line_Start_Segment, Line_End_Segment, Line_Count_Segment
   integer:: Elem_End_File, Line_End_File
   integer:: Y_Res_Stride, Chan_Stride
   integer:: X_Res_Offset
   integer:: Y_Res_Offset
   integer:: i, j, ni, nj, i1, i2, j1, j2, is, js, n, nn

   integer, dimension(2):: start_2d, stride_2d, count_2d
   
   !--- assume Y_resolution to be as X_resolution for Native Data
   Y_Res_Stride = X_Res_Stride   !(always true??)

   !--- set offset for picking the sampled pixel (0=no offset)
   !-- note abi and ahi only off x2 or x4 higher resolution channels
   !-- assume X_Res_Offet_In and Y_Res_Offset_In vary from 0-3
   !-- the following logic should handle channels with Res_Stride = 1,2 or 4
   if (X_Res_Stride == 2) then
      X_Res_Offset = int(float(X_Res_Offset_In)/float(X_Res_Stride) * (X_Res_Stride - 1))
   else
      X_Res_Offset = X_Res_Offset_In
   endif
   X_Res_Offset = min(X_Res_Stride - 1, max(0,X_Res_Offset))

   if (Y_Res_Stride == 2) then
      Y_Res_Offset = int(float(Y_Res_Offset_In)/float(Y_Res_Stride) * (Y_Res_Stride - 1))
   else
      Y_Res_Offset = Y_Res_Offset_In
   endif
   Y_Res_Offset = min(Y_Res_Stride - 1, max(0,Y_Res_Offset))

   !--- initialize
   DQF_Seg = Missing_Value_Int1
   Output_Seg = MISSING_VALUE_REAL4
   Max_Output_Seg = MISSING_VALUE_REAL4
   Min_Output_Seg = -1.0*MISSING_VALUE_REAL4
   Stddev_Output_Seg = MISSING_VALUE_REAL4

   !--- save output size
   nx = size(Output_Seg,1)
   ny = size(Output_Seg,2)

   !--- determine native pixel sizes 
   Native_Elem_End_File = Num_Elem_Fd * X_Res_Stride
   Native_Line_End_File = Num_Line_Fd * Y_Res_Stride

   !--- there is no striding in the reading of the native data
   Native_X_Stride = 1
   Native_Y_Stride = 1

   !--- effect stride the product of chosen thinning stride and native
   !--- resolution stride value
   X_Stride = X_Thin_Stride * X_Res_Stride
   Y_Stride = Y_Thin_Stride * Y_Res_Stride

   !--- allocate array to read data at native resolution
   allocate(Native_Output_Seg(nx*X_Stride,ny*Y_Stride))
   allocate(Native_DQF_Seg(nx*X_Stride,ny*Y_Stride))

   !--- native element terms
   Native_Elem_Start_Segment = (Elem_Start - 1) * X_Res_Stride + 1
   Native_Elem_Count_Segment = Num_Elem * X_Stride
   Native_Elem_End_Segment = Native_Elem_Start_Segment + Native_Elem_Count_Segment - 1
 
   !--- native line terms
   Native_Line_Start_Segment = (Line_Start-1)*Y_Res_Stride + (Segment_Number-1)*Number_Of_Lines_Per_Segment*Y_Stride + 1
   Native_Line_Count_Segment = Number_Of_Lines_Per_Segment * Y_Stride
   Native_Line_End_Segment = Native_Line_Start_Segment + Native_Line_Count_Segment - 1
   Native_Line_End_Segment =   min(Native_Line_End_Segment,Native_Line_End_File)

   Native_Line_Count_Segment = (Native_Line_End_Segment - Native_Line_Start_Segment) + 1

   !--- nominal line terms (needed for reading time below)
   Elem_End_File = Num_Elem * Y_Res_Stride
   Elem_Start_Segment = Elem_Start*X_Res_Stride + X_Res_Offset
   Elem_Count_Segment = Num_Elem
   Line_End_File = Num_Line_Fd * Y_Res_Stride
   Line_Start_Segment = Line_Start*Y_Res_Stride + (Segment_Number-1)*Number_Of_Lines_Per_Segment*Y_Stride + Y_Res_Offset
   Line_Count_Segment = Number_Of_Lines_Per_Segment
 
   !--- read native resolution data
   start_2d =  (/Native_Elem_Start_Segment, Native_Line_Start_Segment/)
   stride_2d = (/Native_X_Stride,Native_Y_Stride/)
   count_2d =   (/Native_Elem_Count_Segment, Native_Line_Count_Segment/)
   call READ_AND_UNSCALE_NETCDF_2d(Ncid, Start_2d, Stride_2d, Count_2d, Chan_Name, Native_Output_Seg)

   !--- read native resolution dqf (only for GOES-R series - not AHI))
   if (Sensor%WMO_ID == 270 .or. Sensor%WMO_ID == 271 .or. Sensor%WMO_ID == 272) then
      call READ_NETCDF(Ncid, Start_2d, Stride_2d, Count_2d, 'DQF', Native_DQF_Seg)
   endif

   !--- process native array into thermal band resolution
   n = X_Stride*Y_Stride
   do i = 1, nx
      is = (i-1)*X_Stride + 1 + X_Res_Offset 
      i1 = (i-1)*X_Stride + 1
      i2 = i1 + X_Stride - 1
!     i1 = min(nx,max(1,i1))
!     i2 = min(nx,max(1,i2))
      do j = 1, ny
           js = (j-1)*Y_Stride + 1 + Y_Res_Offset 
           j1 = (j-1)*Y_Stride + 1
           j2 = j1 + Y_Stride - 1

!          j1 = min(ny,max(1,j1))
!          j2 = min(ny,max(1,j2))
!          n = (i2-i1+1)*(j2-j1+1)
           Output_Seg(i,j) = Native_Output_Seg(is,js)
           DQF_Seg(i,j) = Native_DQF_Seg(is,js)
           if (Chan_Average_Flag >= 1) then

               if (minval(Native_Output_Seg(i1:i2,j1:j2)) == Missing_Value_Real4) then
                   Output_Seg(i,j) = Missing_Value_Real4
                   DQF_Seg(i,j) =  Missing_Value_Int1
               else
                   nn = (i2-i1+1)*(j2-j1+1)
                   Output_Seg(i,j) = sum(Native_Output_Seg(i1:i2,j1:j2))/nn
                   DQF_Seg(i,j) = sum(Native_DQF_Seg(i1:i2,j1:j2))/nn
               endif      

           endif
           if (Chan_Average_Flag == 2) then
               Min_Output_Seg(i,j) = minval(Native_Output_Seg(i1:i2,j1:j2))
               Max_Output_Seg(i,j) = maxval(Native_Output_Seg(i1:i2,j1:j2))
               Stddev_Output_Seg(i,j) = sqrt( sum( (Native_Output_Seg(i1:i2,j1:j2)-Output_Seg(i,j))**2 ) / (n-1) )
           endif
      enddo
    enddo

    !--- deallocate the array used to read native resolution data
    deallocate(Native_Output_Seg)
    deallocate(Native_DQF_Seg)

    !--------------------------------------------------------------------------------
    !  Read Time
    !  NOTE!!! Scan_Time is a global variable in pixel common (should add to
    !  struc)
    !--------------------------------------------------------------------------------

    !--- ensure we don't exceed the data bounds on the read
    j1 = 1
    j2 = min(Line_Count_Segment,Image%Number_Of_Lines_Read_This_Segment)

    select case (Sensor%WMO_ID)

      case (173:174) !-- ahi(s)
         if (Read_Time) then
             call READ_AHI_TIME (Ncid, Line_Start_Segment, j2,  &
                                 Y_Stride, Image%Scan_Time_Ms(j1:j2), Read_Time)
         endif


      case(270:273) !-- abi(s)

             call READ_ABI_TIME (Ncid, Line_Start_Segment, j2,  &
                                 Line_End_File, Y_Stride, Image%Scan_Time_Ms(j1:j2), Read_Time)

      case default
         print *, 'unknown satellite, stopping'
         stop

     end select

end subroutine READ_SEGMENT_LEVEL1B_VER2
!------------------------------------------------------------------------------------
!  read segments of data from the variables in the static nav file
!
!  the output are the latitude, longitude, sensor_zentith and sensor_azimiuth
!  all output passed through the pixel_common_mod memory.
!
!  Note, this uses  Number_Of_Lines_Read_This_Segment instead of
!  Number_Of_Lines_Per_Segment to account for partially filled segments
!------------------------------------------------------------------------------------
subroutine  READ_SEGMENT_STATIC_NAV(Ncid_Static)
   integer, intent(in) :: Ncid_Static

   integer:: Elem_Start_Segment, Elem_Count_Segment
   integer:: Line_Start_Segment, Line_Count_Segment
   integer:: Line_End_File, Line_End_Segment

   real, dimension(:,:), allocatable:: Data_Segment
   real, dimension(:), allocatable:: Data_Segment_1D_I4_X
   real, dimension(:), allocatable:: Data_Segment_1D_I4_Y

   integer:: i, j

   Elem_Start_Segment = Elem_Start 
   Elem_Count_Segment = Num_Elem

   Line_Start_Segment = Line_Start + (Image%Segment_Number-1)*Image%Number_Of_Lines_Per_Segment*Image%Y_Stride
   Line_End_Segment = Line_Start_Segment + Image%Number_Of_Lines_Per_Segment*Image%Y_Stride
   Line_End_File = Line_Start + Num_Line*Image%Y_Stride
   Line_End_Segment = min(Line_End_Segment, Line_End_File)
   Line_Count_Segment = (Line_End_Segment - Line_Start_Segment)/Image%Y_Stride

   Image%Number_Of_Lines_Read_This_Segment = Line_Count_Segment

   allocate(Data_Segment(Elem_Count_Segment, Line_Count_Segment))
   allocate(Data_Segment_1D_I4_X(Elem_Count_Segment))
   allocate(Data_Segment_1D_I4_Y(Line_Count_Segment))

   call READ_NETCDF(Ncid_Static, (/Elem_Start_Segment, Line_Start_Segment/), &
                                 (/Image%X_Stride,Image%Y_Stride/),  &
                                 (/Elem_Count_Segment, Line_Count_Segment/), "latitude",Data_Segment)
   Nav%Lat_1b(:,1:Image%Number_Of_Lines_Read_This_Segment) = Data_Segment

   call READ_NETCDF(Ncid_Static, (/Elem_Start_Segment, Line_Start_Segment/), &
                                 (/Image%X_Stride,Image%Y_Stride/),  &
                                 (/Elem_Count_Segment, Line_Count_Segment/), "longitude",Data_Segment)
   Nav%Lon_1b(:,1:Image%Number_Of_Lines_Read_This_Segment) = Data_Segment

   call READ_NETCDF(Ncid_Static, (/Elem_Start_Segment, Line_Start_Segment/), &
                                 (/Image%X_Stride,Image%Y_Stride/),  &
                                 (/Elem_Count_Segment, Line_Count_Segment/), "sensor_zenith_angle",Data_Segment)
   Geo%Satzen(:,1:Image%Number_Of_Lines_Read_This_Segment) = Data_Segment

   call READ_NETCDF(Ncid_Static, (/Elem_Start_Segment, Line_Start_Segment/), &
                                 (/Image%X_Stride,Image%Y_Stride/),  &
                                 (/Elem_Count_Segment, Line_Count_Segment/), "sensor_azimuth_angle",Data_Segment)
   Geo%Sataz(:,1:Image%Number_Of_Lines_Read_This_Segment) = Data_Segment

   !---- read x and y (ignored for ahi)
   if (Sensor%WMO_Id == 270 .or. Sensor%WMO_Id == 271 .or. Sensor%WMO_Id == 272) then 
     call READ_NETCDF(Ncid_Static, (/Elem_Start_Segment/), &
                                   (/Image%X_Stride/),  &
                                   (/Elem_Count_Segment/), "x_fd",Data_Segment_1D_I4_X)

     call READ_NETCDF(Ncid_Static, (/Line_Start_Segment/), &
                                   (/Image%Y_Stride/),  &
                                   (/Line_Count_Segment/), "y_fd",Data_Segment_1D_I4_Y)
     !-- conveet to 2d
     do j = 1,Image%Number_Of_Lines_Read_This_Segment
       Nav%X(:,j) = Data_Segment_1D_I4_X
     enddo

     do i = 1,Image%Number_Of_Elements
       Nav%Y(i,1:Image%Number_Of_Lines_Read_This_Segment) = Data_Segment_1D_I4_Y
     enddo

   endif

   ! --- deallocate
   if (allocated(Data_Segment)) deallocate(Data_Segment)
   if (allocated(Data_Segment_1D_I4_X)) deallocate(Data_Segment_1D_I4_X)
   if (allocated(Data_Segment_1D_I4_Y)) deallocate(Data_Segment_1D_I4_Y)

   Nav%Lat = Nav%Lat_1b
   Nav%Lon = Nav%Lon_1b

end subroutine READ_SEGMENT_STATIC_NAV
!------------------------------------------------------------------------------------
!  determine the bounds of the data to process within the entire file
!------------------------------------------------------------------------------------
subroutine DETERMINE_BOUNDS_STATIC_NAV(Nav_File_Name, Nav_Flag,Lat_South, Lon_West, Lat_north, Lon_East, &
                                       X_Stride, Y_Stride,Num_Elem_Out, Num_Line_Out)


   character (len=*), intent(in):: Nav_File_Name 
   integer, intent(in):: Nav_Flag
   real, intent(in):: Lat_South, Lon_West, Lat_North, Lon_East
   integer, intent(in):: X_Stride, Y_Stride
   integer, intent(out):: Num_Elem_Out, Num_Line_Out

   real, dimension(:,:), allocatable:: Latitude_Fd, Longitude_Fd
   real, dimension(:), allocatable:: Lat_Nadir, Lat_Diff, Lon_Nadir, Lon_Diff
   integer:: i_start, i_end, j_start, j_end
   integer, dimension(1):: i_east, i_west, j_north, j_south
   integer:: ncid, status
   logical:: subset_flag

   !--- Added to match AHI from Andi
   integer :: ii
   logical, allocatable :: inside (:,:)
   integer, allocatable :: line_g(:), elem_g(:)

   call READ_NETCDF_DIMENSION(Nav_File_Name, "number_of_longitudes", Num_Elem_Fd)
   call READ_NETCDF_DIMENSION (Nav_File_Name, "number_of_latitudes", Num_Line_Fd)
   call READ_NETCDF_GLOBAL_ATTRIBUTE(Nav_File_Name, "latitude_of_projection_origin", Sub_Sat_Lat)
   call READ_NETCDF_GLOBAL_ATTRIBUTE(Nav_File_Name, "longitude_of_projection_origin", Sub_Sat_Lon)

   Sensor%Geo_Sub_Satellite_Longitude = Sub_Sat_Lon 
   Sensor%Geo_Sub_Satellite_Latitude = Sub_Sat_Lat

   !--------------------------------------------------------------------------------------------
   ! Check if any subsetting parameters are set, if not, do not attempt subsetting
   !--------------------------------------------------------------------------------------------
   subset_flag = .true.
   if (Nav%Lat_North_Limit == 90.0 .and.  &
       Nav%Lat_South_Limit == -90.0 .and. &
       Nav%Lon_West_Limit == -180.0) then
       call MESG( "STATIC NAV used but no subsetting performed", level = verb_lev % VERBOSE )
       Elem_Start = 1
       Num_Elem  = Num_Elem_Fd
       Line_Start = 1
       Num_Line = Num_Line_Fd
       subset_flag = .false.
   endif

   !--------------------------------------------------------------------------------------------
   ! continue with defining subset region
   !--------------------------------------------------------------------------------------------
   if (subset_flag) then

      allocate(Latitude_Fd(Num_Elem_Fd, Num_Line_Fd))
      allocate(Longitude_Fd(Num_Elem_Fd, Num_Line_Fd))

      !--- could some who rewrote this please say why and what these indices mean?
      !--- Add  !___WHAT IS THIS___
      allocate(inside(Num_Elem_Fd, Num_Line_Fd))
      allocate(line_g(Num_Line_Fd))
      allocate(elem_g(Num_Elem_Fd))

      call READ_NETCDF(Ncid_Static, (/1,1/), (/1,1/), (/Num_Elem_Fd, Num_Line_Fd/), "latitude",Latitude_Fd)
      call READ_NETCDF(Ncid_Static, (/1,1/), (/1,1/), (/Num_Elem_Fd, Num_Line_Fd/), "longitude",Longitude_Fd)

      !--- determine bounds
      allocate(Lat_Nadir(Num_Line_Fd))
      allocate(Lat_Diff(Num_Line_Fd))
      allocate(Lon_Nadir(Num_Line_Fd))
      allocate(Lon_Diff(Num_Line_Fd))

      !--------------------------------------------------------------
      ! determine indices of the subsetting region within the data    
      !--------------------------------------------------------------
      if (Lon_East .ger. Lon_West) then
        inside = Longitude_Fd .GT. Lon_West .and. &
          Longitude_Fd .LT. Lon_East .and. &
          Latitude_Fd .GT. Lat_South .and. &
          Latitude_Fd .LT. Lat_North

      else

        inside = (((Longitude_Fd .ltr. Lon_East) .and. &
            (Longitude_Fd .ger. -180.0) )  .or. &
            ((Longitude_Fd .gtr. Lon_West) .and. &
             (Longitude_Fd .ltr. 180.0) )) .and. &
            (Latitude_Fd .gtr. Lat_South)  .and. &
            (Latitude_Fd .ltr. Lat_North)

      endif

      elem_g = count (inside ,2 )
      line_g = count (inside ,1 )

      do ii =1 , Num_Elem_Fd
        if ( elem_g(ii) /= 0 ) then
          i_start = ii
          exit
        end if
      end do

      do ii = 1, Num_Line_Fd
        if ( line_g(ii) /= 0 ) then
          j_start = ii
          exit
        end if
      end do

      do ii = Num_Elem_Fd, 1, -1
        if ( elem_g(ii) /= 0 ) then
          i_end = ii
          Num_Elem = ii - i_start
          exit
        end if
      end do

      do ii =Num_Line_Fd , 1, -1
        if ( line_g(ii) .ne. 0 ) then
          j_end =  ii
          Num_Line = ii - j_start
          exit
        end if
      end do
 
      !--------------------------------------------------------------
      ! confirm subsetting worked
      !--------------------------------------------------------------
      if (i_start < 1 .or. i_end < 1 .or. Nav_Flag == sym%NO) then
        call MESG( "STATIC NAV longitude subsetting failed. All longitudes will be processed: ", level = verb_lev % DEFAULT )
        Elem_Start = 1
        Num_Elem  = Num_Elem_Fd
      else
        Elem_Start = i_start
      endif

      if (j_start < 1 .or. j_end < 1 .or. Nav_Flag == sym%NO) then 
        call MESG( "STATIC NAV latitude subsetting failed. All latitudes will be processed: ", level = verb_lev % DEFAULT )
        Line_Start = 1
        Num_Line = Num_Line_Fd
      else
        Line_Start = j_start
      endif

      !--- clean up memory
      deallocate(Latitude_Fd)
      deallocate(Longitude_Fd)
      deallocate(Lat_Nadir)
      deallocate(Lat_Diff)
      deallocate(Lon_Nadir)
      deallocate(Lon_Diff)
      deallocate(inside)
      deallocate(line_g)
      deallocate(elem_g)

   endif

   !--- account for stride
   Num_Elem = Num_Elem / X_Stride
   Num_Line = Num_Line / Y_Stride

   !--- store these in output arguments since needed elsewhere
   Num_Elem_Out = Num_Elem
   Num_Line_Out = Num_Line

end subroutine DETERMINE_BOUNDS_STATIC_NAV

!------------------------------------------------------------------------------------
! --- find dark composite file
!------------------------------------------------------------------------------------
subroutine FIND_DARK_COMPOSITE()

   use FILE_UTILS, only: FILE_SEARCH, FILE_TEST
  
   ! --- construct file name
   Dark_Composite_Name = "no_file"

   select case (Sensor%WMO_Id)
     case(270)
       Dark_Comp_Data_Dir = trim(Ancil_Data_Dir)//'/static/dark_comp/goes_east/'
       Dark_Composite_Name = "goes_east_fulldisk_"//trim(image % time_start % Monthname)//"_"//image % time_start % hh //"Z_refl065_dark_composite_2km.nc"
     case(271:272)
       Dark_Comp_Data_Dir = trim(Ancil_Data_Dir)//'/static/dark_comp/goes_west/'
       Dark_Composite_Name = "goes_west_fulldisk_"//trim(image % time_start % Monthname)//"_"//image % time_start % hh//"Z_refl065_dark_composite_2km.nc"
     case default
       call MESG('DARK COMPOSITE FILE DOES NOT EXIST FOR THIS SENSOR, SETTING TO MISSING', level = verb_lev % DEFAULT )
   end select
     
   ! --- check if DC file exists
   Dark_File_Exist = FILE_TEST(trim(Dark_Comp_Data_Dir)//trim(Dark_Composite_Name))
   if (.not. Dark_File_Exist) then
      call MESG('DID NOT FIND DARK COMPOSITE FILE, SETTING TO MISSING', level = verb_lev % DEFAULT )
      !print*,trim(Dark_Comp_Data_Dir)//trim(Dark_Composite_Name)
   endif

end subroutine FIND_DARK_COMPOSITE

!------------------------------------------------------------------------------------
! --- read dark composite file
!------------------------------------------------------------------------------------
subroutine READ_DARK_COMPOSITE(X_Thin_Stride, Y_Thin_Stride, Output_Mean, Output_Stddev)

   integer, intent(in):: X_Thin_Stride, Y_Thin_Stride
   real, dimension(:,:), intent(out):: Output_Mean, Output_Stddev
   real, dimension(:,:), allocatable:: Native_Output_Mean
   real, dimension(:,:), allocatable:: Native_Output_Stddev
   integer:: Ncid
   integer:: nx, ny

   integer:: Native_X_Stride, Native_Y_Stride
   integer:: X_Stride, Y_Stride
   integer:: Native_Elem_Start_Segment, Native_Elem_End_Segment, Native_Elem_Count_Segment
   integer:: Native_Line_Start_Segment, Native_Line_End_Segment, Native_Line_Count_Segment
   integer:: Native_Elem_End_File, Native_Line_End_File
   integer:: Elem_Start_Segment, Elem_End_Segment, Elem_Count_Segment
   integer:: Line_Start_Segment, Line_End_Segment, Line_Count_Segment
   integer:: Elem_End_File, Line_End_File
   integer:: X_Res_Stride, Y_Res_Stride, Chan_Stride
   integer:: X_Res_Offset
   integer:: Y_Res_Offset
   integer:: i, j, ni, nj, i1, i2, j1, j2, is, js, n
   character(len=100) :: Nan_Str1, Nan_Str2
   integer, dimension(2):: start_2d, stride_2d, count_2d

   !--- resolution of dark composite is that of thermal resolution
   X_Res_Stride = 1

   !--- assume Y_resolution to be as X_resolution for Native Data
   Y_Res_Stride = X_Res_Stride   !(always true??)

   !--- no offsetting in dark compostites
   X_Res_Offset = 0
   Y_Res_Offset = 0

   !--- initialize
   Output_Mean = MISSING_VALUE_REAL4
   Output_Stddev = MISSING_VALUE_REAL4

   !--- save output size
   nx = size(Output_Mean,1)
   ny = size(Output_Mean,2)

   !--- determine native pixel sizes 
   Native_Elem_End_File = Num_Elem_Fd * X_Res_Stride
   Native_Line_End_File = Num_Line_Fd * Y_Res_Stride

   !--- there is no striding in the reading of the native data
   Native_X_Stride = 1
   Native_Y_Stride = 1

   !--- effect stride the product of chosen thinning stride and native
   !--- resolution stride value
   X_Stride = X_Thin_Stride * X_Res_Stride
   Y_Stride = Y_Thin_Stride * Y_Res_Stride

   !--- allocate array to read data at native resolution
   allocate(Native_Output_Mean(nx*X_Stride,ny*Y_Stride))
   allocate(Native_Output_Stddev(nx*X_Stride,ny*Y_Stride))

   !--- native element terms
   Native_Elem_Start_Segment = (Elem_Start - 1) * X_Res_Stride + 1
   Native_Elem_Count_Segment = Num_Elem * X_Stride
   Native_Elem_End_Segment = Native_Elem_Start_Segment + Native_Elem_Count_Segment - 1
 
   !--- native line terms
   Native_Line_Start_Segment = (Line_Start-1)*Y_Res_Stride + (Image%Segment_Number-1)*Image%Number_Of_Lines_Per_Segment*Y_Stride + 1
   Native_Line_Count_Segment = Image%Number_Of_Lines_Per_Segment * Y_Stride
   Native_Line_End_Segment = Native_Line_Start_Segment + Native_Line_Count_Segment - 1
   Native_Line_End_Segment =   min(Native_Line_End_Segment,Native_Line_End_File)

   Native_Line_Count_Segment = (Native_Line_End_Segment - Native_Line_Start_Segment) + 1

   !--- nominal line terms (needed for reading time below)
   Elem_End_File = Num_Elem * Y_Res_Stride
   Elem_Start_Segment = Elem_Start*X_Res_Stride + X_Res_Offset
   Elem_Count_Segment = Num_Elem
   Line_End_File = Num_Line_Fd * Y_Res_Stride
   Line_Start_Segment = Line_Start*Y_Res_Stride + (Image%Segment_Number-1)*Image%Number_Of_Lines_Per_Segment*Y_Stride + Y_Res_Offset
   Line_Count_Segment = Image%Number_Of_Lines_Per_Segment
 
   !--- read native resolution data
   start_2d =  (/Native_Elem_Start_Segment, Native_Line_Start_Segment/)
   stride_2d = (/Native_X_Stride,Native_Y_Stride/)
   count_2d =   (/Native_Elem_Count_Segment, Native_Line_Count_Segment/)

   !--- read raw data
   call OPEN_NETCDF(trim(Dark_Comp_Data_Dir)//trim(Dark_Composite_Name),Ncid)
   call READ_AND_UNSCALE_NETCDF_2d(Ncid, Start_2d, Stride_2d, Count_2d, 'darkest_comp_mean_2km', Native_Output_Mean)
   call READ_AND_UNSCALE_NETCDF_2d(Ncid, Start_2d, Stride_2d, Count_2d, 'darkest_comp_stddev_2km', Native_Output_Stddev)
   call CLOSE_NETCDF(Ncid)

   !--- handle NaNs
   where(Native_Output_Mean /= Native_Output_Mean)
      Native_Output_Mean = MISSING_VALUE_REAL4
   endwhere

   where(Native_Output_Stddev /= Native_Output_Stddev)
      Native_Output_Stddev = MISSING_VALUE_REAL4
   endwhere

   !--- process native array into thermal band resolution
   n = X_Stride*Y_Stride
   do i = 1, nx
      is = (i-1)*X_Stride + 1 + X_Res_Offset 
      i1 = (i-1)*X_Stride + 1
      i2 = i1 + X_Stride - 1
      do j = 1, ny
           js = (j-1)*Y_Stride + 1 + Y_Res_Offset 
           j1 = (j-1)*Y_Stride + 1
           j2 = j1 + Y_Stride - 1
           Output_Mean(i,j) = Native_Output_Mean(is,js)
           Output_Stddev(i,j) = Native_Output_Stddev(is,js)
      enddo
  enddo

  !--- deallocate the array used to read native resolution data
  deallocate(Native_Output_Mean)
  deallocate(Native_Output_Stddev)

end subroutine READ_DARK_COMPOSITE

!------------------------------------------------------------------------------------
! find quantile regression neural network ctp file (needs development)
!------------------------------------------------------------------------------------
subroutine FIND_QRNN_CTP()
    QRNN_CTP_Name = 'QRNN_CTP_G16_s20201060000165_e20201060009484_c20201060009556.nc'
    QRNN_CTP_File_Exist = .true.
end subroutine FIND_QRNN_CTP

!------------------------------------------------------------------------------------
! read quantile regression neural network ctp file (needs development)
!------------------------------------------------------------------------------------
subroutine READ_QRNN_CTP(X_Thin_Stride, Y_Thin_Stride, Output_Mean, Output_Stddev)

   integer, intent(in):: X_Thin_Stride, Y_Thin_Stride
   real, dimension(:,:), intent(out):: Output_Mean, Output_Stddev
   real, dimension(:,:), allocatable:: Native_Output_Mean
   real, dimension(:,:), allocatable:: Native_Output_Stddev
   integer:: Ncid
   integer:: nx, ny

   integer:: Native_X_Stride, Native_Y_Stride
   integer:: X_Stride, Y_Stride
   integer:: Native_Elem_Start_Segment, Native_Elem_End_Segment, Native_Elem_Count_Segment
   integer:: Native_Line_Start_Segment, Native_Line_End_Segment, Native_Line_Count_Segment
   integer:: Native_Elem_End_File, Native_Line_End_File
   integer:: Elem_Start_Segment, Elem_End_Segment, Elem_Count_Segment
   integer:: Line_Start_Segment, Line_End_Segment, Line_Count_Segment
   integer:: Elem_End_File, Line_End_File
   integer:: X_Res_Stride, Y_Res_Stride, Chan_Stride
   integer:: X_Res_Offset
   integer:: Y_Res_Offset
   integer:: i, j, ni, nj, i1, i2, j1, j2, is, js, n
   character(len=100) :: Nan_Str1, Nan_Str2
   integer, dimension(2):: start_2d, stride_2d, count_2d

   !--- resolution of dark composite is that of thermal resolution
   X_Res_Stride = 1

   !--- assume Y_resolution to be as X_resolution for Native Data
   Y_Res_Stride = X_Res_Stride   !(always true??)

   !--- no offsetting in dark compostites
   X_Res_Offset = 0
   Y_Res_Offset = 0

   !--- initialize
   Output_Mean = MISSING_VALUE_REAL4
   Output_Stddev = MISSING_VALUE_REAL4

   !--- save output size
   nx = size(Output_Mean,1)
   ny = size(Output_Mean,2)

   !--- determine native pixel sizes 
   Native_Elem_End_File = Num_Elem_Fd * X_Res_Stride
   Native_Line_End_File = Num_Line_Fd * Y_Res_Stride

   !--- there is no striding in the reading of the native data
   Native_X_Stride = 1
   Native_Y_Stride = 1

   !--- effect stride the product of chosen thinning stride and native
   !--- resolution stride value
   X_Stride = X_Thin_Stride * X_Res_Stride
   Y_Stride = Y_Thin_Stride * Y_Res_Stride

   !--- allocate array to read data at native resolution
   allocate(Native_Output_Mean(nx*X_Stride,ny*Y_Stride))
   allocate(Native_Output_Stddev(nx*X_Stride,ny*Y_Stride))

   !--- native element terms
   Native_Elem_Start_Segment = (Elem_Start - 1) * X_Res_Stride + 1
   Native_Elem_Count_Segment = Num_Elem * X_Stride
   Native_Elem_End_Segment = Native_Elem_Start_Segment + Native_Elem_Count_Segment - 1
 
   !--- native line terms
   Native_Line_Start_Segment = (Line_Start-1)*Y_Res_Stride + (Image%Segment_Number-1)*Image%Number_Of_Lines_Per_Segment*Y_Stride + 1
   Native_Line_Count_Segment = Image%Number_Of_Lines_Per_Segment * Y_Stride
   Native_Line_End_Segment = Native_Line_Start_Segment + Native_Line_Count_Segment - 1
   Native_Line_End_Segment =   min(Native_Line_End_Segment,Native_Line_End_File)

   Native_Line_Count_Segment = (Native_Line_End_Segment - Native_Line_Start_Segment) + 1

   !--- nominal line terms (needed for reading time below)
   Elem_End_File = Num_Elem * Y_Res_Stride
   Elem_Start_Segment = Elem_Start*X_Res_Stride + X_Res_Offset
   Elem_Count_Segment = Num_Elem
   Line_End_File = Num_Line_Fd * Y_Res_Stride
   Line_Start_Segment = Line_Start*Y_Res_Stride + (Image%Segment_Number-1)*Image%Number_Of_Lines_Per_Segment*Y_Stride + Y_Res_Offset
   Line_Count_Segment = Image%Number_Of_Lines_Per_Segment
 
   !--- read native resolution data
   start_2d =  (/Native_Elem_Start_Segment, Native_Line_Start_Segment/)
   stride_2d = (/Native_X_Stride,Native_Y_Stride/)
   count_2d =   (/Native_Elem_Count_Segment, Native_Line_Count_Segment/)

   !--- read raw data
   call OPEN_NETCDF(trim(QRNN_CTP_Data_Dir)//trim(QRNN_CTP_Name),Ncid)
   call READ_AND_UNSCALE_NETCDF_2d(Ncid, Start_2d, Stride_2d, Count_2d, 'median_ctp', Native_Output_Mean)
   call READ_AND_UNSCALE_NETCDF_2d(Ncid, Start_2d, Stride_2d, Count_2d, 'estimated_absolute_error', Native_Output_Stddev)
   call CLOSE_NETCDF(Ncid)

   !--- handle NaNs
   where(Native_Output_Mean /= Native_Output_Mean)
      Native_Output_Mean = MISSING_VALUE_REAL4
   endwhere

   where(Native_Output_Stddev /= Native_Output_Stddev)
      Native_Output_Stddev = MISSING_VALUE_REAL4
   endwhere

   !--- process native array into thermal band resolution
   n = X_Stride*Y_Stride
   do i = 1, nx
      is = (i-1)*X_Stride + 1 + X_Res_Offset 
      i1 = (i-1)*X_Stride + 1
      i2 = i1 + X_Stride - 1
      do j = 1, ny
           js = (j-1)*Y_Stride + 1 + Y_Res_Offset 
           j1 = (j-1)*Y_Stride + 1
           j2 = j1 + Y_Stride - 1
           Output_Mean(i,j) = Native_Output_Mean(is,js)
           Output_Stddev(i,j) = Native_Output_Stddev(is,js)
      enddo
  enddo

  !--- deallocate the array used to read native resolution data
  deallocate(Native_Output_Mean)
  deallocate(Native_Output_Stddev)

end subroutine READ_QRNN_CTP

!------------------------------------------------------------------------------------
!  determine the bounds of the data to process within the entire file
!------------------------------------------------------------------------------------
subroutine DETERMINE_BOUNDS_STATIC_NAV_ORIGINAL(Nav_File_Name, Nav_Flag, &
                                                Lat_South, Lon_West, Lat_North, Lon_East, &
                                                X_Stride, Y_Stride,Num_Elem_Out, Num_Line_Out)


   character (len=*), intent(in):: Nav_File_Name 
   real, intent(in):: Lat_South, Lon_West, Lat_North, Lon_East
   integer, intent(in):: X_Stride, Y_Stride, Nav_Flag
   integer, intent(out):: Num_Elem_Out, Num_Line_Out

   real, dimension(:,:), allocatable:: Latitude_Fd, Longitude_Fd
   real, dimension(:), allocatable:: Lat_Nadir, Lat_Diff, Lon_Nadir, Lon_Diff
   integer:: i_start, i_end, j_start, j_end
   integer, dimension(1):: i_east, i_west, j_north, j_south
   integer:: ncid, status

   call read_netcdf_dimension(Nav_File_Name, "number_of_longitudes", Num_Elem_Fd)
   call read_netcdf_dimension (Nav_File_Name, "number_of_latitudes", Num_Line_Fd)
   call READ_NETCDF_GLOBAL_ATTRIBUTE(Nav_File_Name, "latitude_of_projection_origin", Sub_Sat_Lat)
   call READ_NETCDF_GLOBAL_ATTRIBUTE(Nav_File_Name, "longitude_of_projection_origin", Sub_Sat_Lon)

   Sensor%Geo_Sub_Satellite_Longitude = Sub_Sat_Lon 
   Sensor%Geo_Sub_Satellite_Latitude = Sub_Sat_Lat

   allocate(Latitude_Fd(Num_Elem_Fd, Num_Line_Fd))
   allocate(Longitude_Fd(Num_Elem_Fd, Num_Line_Fd))

   call read_netcdf(ncid_static, (/1,1/), (/1,1/), (/Num_Elem_Fd, Num_Line_Fd/), "latitude",Latitude_Fd)
   call read_netcdf(ncid_static, (/1,1/), (/1,1/), (/Num_Elem_Fd, Num_Line_Fd/), "longitude",Longitude_Fd)

   !--- determine bounds
   allocate(Lat_Nadir(Num_Line_Fd))
   allocate(Lat_Diff(Num_Line_Fd))
   allocate(Lon_Nadir(Num_Line_Fd))
   allocate(Lon_Diff(Num_Line_Fd))

   j_north = Missing_Value_Int4
   j_South = Missing_Value_Int4
   i_West = Missing_Value_Int4
   i_East = Missing_Value_Int4

   Lat_Nadir = Latitude_Fd(Num_Elem_fd/2,:)
   Lon_Nadir = Longitude_Fd(:,Num_Line_fd/2)

   Lat_Diff = abs(Lat_Nadir - Lat_North)
   j_North = minloc(Lat_Diff)
 
   Lat_Diff = abs(Lat_Nadir - Lat_South)
   j_South = minloc(Lat_Diff)

   Lon_Diff = abs(Lon_Nadir - Lon_West)
   i_West = minloc(Lon_Diff)

   Lon_Diff = abs(Lon_Nadir - Lon_East)
   i_East = minloc(Lon_Diff)

   !---- reorder as necessary
   i_start = min(i_West(1),i_East(1))
   i_end = max(i_West(1),i_East(1))
   j_start = min(j_South(1),j_North(1))
   j_end = max(j_South(1),j_North(1))

   if (i_start < 1 .or. i_end < 1 .or. Nav_Flag == sym%NO) then
     i_start = 1
     i_end = Num_Elem_Fd
   endif

   if (j_start < 1 .or. j_end < 1 .or. Nav_Flag == sym%NO) then 
     j_start = 1
     j_end = Num_Line_Fd
   endif

   !--- configure output
   Elem_Start = i_start
   Num_Elem = i_end - i_start + 1

   line_start = j_start
   Num_Line = j_end - j_start + 1

   !--- account for stride
   Num_Elem = Num_Elem / X_Stride
   Num_Line = Num_Line / Y_Stride

   !--- store these in output arguments since needed elsewhere
   Num_Elem_Out = Num_Elem
   Num_Line_Out = Num_Line

   !--- clean up memory
   deallocate(Latitude_Fd)
   deallocate(Longitude_Fd)
   deallocate(Lat_Nadir)
   deallocate(Lat_Diff)
   deallocate(Lon_Nadir)
   deallocate(Lon_Diff)

end subroutine DETERMINE_BOUNDS_STATIC_NAV_ORIGINAL
!------------------------------------------------------------------------------------
! Set the Ch(XX)%Source for channels based on use of fusion data
!------------------------------------------------------------------------------------
subroutine SET_STATIC_NAV_SOURCE_VALUES()

    integer:: Chan_Idx, CLAVRx_Chan_Idx

    do Chan_Idx = 1, Num_Chan_Sensor

      CLAVRx_Chan_Idx = CLAVRx_Chan_Map(Chan_Idx)

      if (Chan_On_Sensor(Chan_Idx) == sym%NO) cycle

      ch(CLAVRx_Chan_Idx)%Source = 0    !note source of this channel is from imager
      if (index(Level1b_File(Chan_Idx),'FR_') > 0) then
         ch(CLAVRx_Chan_Idx)%Source = 1    !note source of this channel is from sounder
      endif

    enddo

end subroutine SET_STATIC_NAV_SOURCE_VALUES
!------------------------------------------------------------------------------------

end module CLAVRX_STATIC_NAV_MODULE
