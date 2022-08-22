! $Id: ahi_hsd_reader.f90 3389 2019-06-26 17:25:41Z wstraka $
!-------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: ahi_hsd_reader.f90 (src)
!       AHI_HSD_MODULE (program)
!
! PURPOSE: 
!
! DESCRIPTION: 
!
! AUTHORS:
!  Ray Garcia, CIMSS, rayg@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
!
! AHI Channel Mapping
!
! wvl       ahi  modis/clavrx 
!
! 0.47       1      3   
! 0.51       2      4
! 0.64       3      1  
! 0.86       4      2  
! 1.6        5      6
! 2.2        6     r7
! 3.9        7     20
! 6.2        8     37
! 6.9        9     27
! 7.3       10     28
! 8.6       11     29
! 9.6       12     30
! 10.4      13     38
! 11.2      14     31
! 12.3      15     32 
! 13.3      16     33
!
!-------------------------------------------------------------------------------
MODULE AHI_HSD_READER

  use CLAVRX_MESSAGE_MOD,only: &
      mesg &
      ,verb_lev
  
  use CONSTANTS_MOD,only: &
    int4,real8,int2 &
    , real4 , int1 &
    , MISSING_VALUE_INT4 &
    , MISSING_VALUE_REAL4 &
    , sym &
    ,exe_prompt
  
  use FILE_TOOLS, only: getlun

  use iso_c_binding
  
  !libHimawari information. While both sets of reference codes use all the modules, keeping with CLAVR-x mentaility, use only what is needed.
  
  USE HSD_SCENE, only: openhimawariscene, copyhimawariscenealbedos, closehimawariscene, copyhimawariscenebrightnesstemps, &
       copyhimawarisceneradiances, queryhimawariscenemetadata, queryhimawariscenelinetimeoffsets, queryhimawariscenecalibration, &
       queryhimawariscenenavigation, averageHimawariFloatArray
  USE HSD_STRUCTS, only: hsd_metadata, hsd_sentinels, hsd_navigation, hsd_calibration, hsd_source, hsd_destination
  
IMPLICIT NONE
  
INTEGER(KIND=INT4), PARAMETER, PRIVATE :: NUM_COLS_BASIC_FD = 5500
INTEGER(KIND=INT4), PARAMETER, PRIVATE :: NUM_LINES_BASIC_FD = 5500

INTEGER(KIND=INT4), PARAMETER, PRIVATE :: NUM_AHI_BANDS = 16

!Example code had methods for HSD and HRIT (HCAST). While not supported, will 
! keep the hooks in there
INTEGER(KIND=INT4), PARAMETER, PUBLIC :: FILE_TYPE_AHI_HSF = 1
INTEGER(KIND=INT4), PARAMETER, PUBLIC :: FILE_TYPE_AHI_HRIT = 2

!Based on talking with RKG, will need to have method to convert between native
!and IR (target) res. Striding is done differently. Leaving HCAST in as well

! native resolution of each band in km
! HSD native resolutions (/1.0, 1.0, 0.5, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0/)
REAL(KIND=REAL4), DIMENSION(NUM_AHI_BANDS), PARAMETER, PRIVATE :: HSD_NATIVE_IR = (/2.0, 2.0, 4.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)

!HCAST native resolutions  = (/-999.0, -999.0, 1.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0/)
REAL(KIND=REAL4), DIMENSION(NUM_AHI_BANDS), PARAMETER, PRIVATE :: HCAST_NATIVE_IR = (/-999.0, -999.0, 4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)

!HCAST File headers
! note that B01 and B02 are not transmitted via HCAST
 CHARACTER(LEN=4), DIMENSION(NUM_AHI_BANDS), PARAMETER, PRIVATE :: HCAST_HEAD = (/"B01_", "B02_", "VIS_", "B04_", "B05_", "B06_", &
                                                                               "IR4_", "IR3_", "B09_", "B10_", "B11_", "B12_", &
                                                                               "IR1_", "B14_", "IR2_","B16_"/)
PUBLIC :: Get_Him_Raddata

PRIVATE :: Him_Rad_Subsample, &
           Him_Rad_Average
   
CONTAINS

!-----------------------------------------------------------------
! Read AHI meta data = May not be needed, since the CLAVR-x doesn't use this
! information. Will just leave calls in, but not have a function. Most is
! similar to test code. You open, read in metadata, close
!-----------------------------------------------------------------
  
!  Handle = OpenHimawariScene(Path_And_Filename_C, C_Null_Funptr, Sentinels)

  !-----------------------------------------------------------------
  ! Read in meta data
  !-----------------------------------------------------------------
  
!  Result = QueryHimawariSceneMetadata(Handle, About)

! Meta data in About - Band, lines, columns
!                    - begin_line, end_line (of the scene)
!                    - Observation_Timeline (start time)
! Note that Libhimwari can determin scanline time


!-----------------------------------------------------------------
! Read AHI radiometric data
! Seperated out from CLAVR-x to make things cleaner.
!-----------------------------------------------------------------

FUNCTION Get_Him_Raddata(File_Type,              &
                         Path,                   &
                         Base_Filename,          &
                         Chn_Num,                &
                         Start,                  &
                         Stride,                 &
                         Num_Lines_Seg,          &
                         Average_Mode_Flag,      &
                         Rad,                    &
                         Rad_Min,                &
                         Rad_Max) RESULT(Status)
  
  INTEGER(KIND=INT4), INTENT(IN) :: File_Type
  CHARACTER(LEN=*), INTENT(IN) :: Path
  CHARACTER(LEN=*), INTENT(IN) :: Base_Filename
  INTEGER(KIND=INT4), INTENT(IN) :: Chn_Num
  INTEGER(KIND=INT4), DIMENSION(2), INTENT(IN) :: Start
  INTEGER(KIND=INT4), DIMENSION(2), INTENT(IN) :: Stride
  INTEGER(KIND=INT4), DIMENSION(2) :: Edge
  INTEGER(KIND=INT4)               :: Num_Lines_Seg ! this is so that you don't over reach in the end of data
  INTEGER(KIND=INT4)               :: Average_Mode_Flag ! Flag from options file for subsample (0), average (1) or avg + max/min (2)
  REAL(KIND=REAL4), DIMENSION(:,:), INTENT(INOUT) :: Rad !array of data in IR space
  REAL(KIND=REAL4), DIMENSION(:,:), INTENT(INOUT) :: Rad_Min !min of averaged data (AVG Mode 2)
  REAL(KIND=REAL4), DIMENSION(:,:), INTENT(INOUT) :: Rad_Max !max of averaged data (AVG Mode 2)
  
  INTEGER(KIND=INT4) :: Char_Idx
  CHARACTER(LEN=2) :: Chn_String
  CHARACTER(LEN=2048) :: Path_And_Filename  
  INTEGER(KIND=INT4) :: Status
  INTEGER(KIND=INT4) :: Factor
  CHARACTER(LEN=4) :: HCAST_Band_ID
  INTEGER :: nx, ny

!-----------------------------------------------------------------
! Conditional compilation: Himawari-8/9 processing capability is
! only available if HIMAWARI is installed; otherwise this routine
! prints an error and exits
!-----------------------------------------------------------------  
  CHARACTER(LEN=2048) :: Path_And_Filename_C
  INTEGER(KIND=INT4) :: Funct_Status, Result
  TYPE(C_Ptr) :: Handle = C_Null_Ptr
  
  !LibHimawari types. Left names similar to test code.
  TYPE(hsd_source) :: Src
  TYPE(hsd_destination) :: Dst
  TYPE(Hsd_Sentinels) :: Sentinels
  TYPE(Hsd_Metadata) :: Hsd_Meta

  Status = Sym%SUCCESS

 !--- save output size
   nx = size(Rad,1)
   ny = size(Rad,2)
   
   Edge(1) = nx
   Edge(2) = Num_Lines_Seg ! Makes sure you don't go to far out of the data

   
  IF (File_Type == FILE_TYPE_AHI_HSF) THEN

    DO Char_Idx=2, LEN_TRIM(Base_Filename)  
       IF (Base_Filename(Char_Idx-1:Char_Idx) == "_B") EXIT  
    END DO
     
    WRITE (Chn_String, "(I2.2)") Chn_Num
    Path_And_Filename = TRIM(Path)//"/"//TRIM(Base_Filename(1:Char_Idx))//TRIM(Chn_String)//TRIM(Base_Filename(Char_Idx+3:))
        

    !Factor is the conversion between native and IR resolution
    
    !Factor = NINT(HSD_IR_RES / HSD_NATIVE_RES(Chn_Num))
    
    Factor = HSD_NATIVE_IR(Chn_Num)
    
  ELSE IF (File_Type == FILE_TYPE_AHI_HRIT) THEN

    ! HCast filename format is IMG_DKppccccYYYYMMDDhhmm_0nn, where cccc is band identifier
    !Section borrowed from GEOCAT for now. 
    !HCAST NOT TESTED
    HCAST_Band_ID = HCAST_HEAD(Chn_Num)
    IF(Status == Sym%Failure) THEN
       call MESG( "AHI channel number out of range", level = verb_lev % DEFAULT , color = 4 )
       RETURN
    END IF

    Path_And_Filename = TRIM(Path)//"/"//Base_Filename(1:8)//HCAST_Band_ID//Base_Filename(13:)
    
    
    Factor = HCAST_NATIVE_IR(Chn_Num)   
  
  ELSE
  
    call MESG( "JMA NetCDF AHI data format not supported, use LibHimawari netCDF. Stopping", level = verb_lev % ERROR , color = 4 )
    STOP
  
  ENDIF
  
  
  !----- Begin HSD calls. Theoretically this could be in it's own routine, 
  !------ like GEOCAT, butfor now, leave as one self contained routine, like RKG test code
  !----- possible item for cleanup - WCS3
  
  !--- note- most comments come direct from test code execpt for ones about Src, Dst,
  !          which are WCS3 comments
 
  !initialize status to success
  Status = Sym%SUCCESS
  
  Path_And_Filename_C = TRIM(Path_And_Filename)//C_NULL_CHAR
  
  !-----------------------------------------------------------------
  ! Override default NaN sentinels in data to something more fortranly
  !-----------------------------------------------------------------
  
  Sentinels%Derived_Invalid = MISSING_VALUE_REAL4
  Sentinels%Derived_Outside_Scan = MISSING_VALUE_REAL4
  

  !-----------------------------------------------------------------
  ! Create a scene handle given a filename prefix, optional (null)
  ! log function, and sentinels
  !-----------------------------------------------------------------
  
  Handle = OpenHimawariScene(Path_And_Filename_C, C_Null_Funptr, Sentinels)
  IF (.not. c_associated(Handle)) THEN 
    call MESG( "Cannot open HSD radiance scene from "//TRIM(Path_And_Filename), &
               level = verb_lev % ERROR , color = 4 )
    Status = Sym%FAILURE
    RETURN
  ENDIF

  !Subsample the data
  IF (Average_Mode_Flag == 0) THEN 
  
    Status = Him_Rad_Subsample(Handle,             &
                               Path_And_Filename,  &
                               Sentinels,          &
                               Start,              &
                               Stride,             &
                               Edge,               &
                               Factor,             &
                               Hsd_Meta,           &
                               Rad)  
  
  !Average the data
  ELSE IF ((Average_Mode_Flag == 1) .or. &
           (Average_Mode_Flag == 2) )THEN
    Status = Him_Rad_Average  (Handle,             &
                               Path_And_Filename,  &
                               Sentinels,          &
                               Start,              &
                               Stride,             &
                               Edge,               &
                               Factor,             &
                               Hsd_Meta,           &
                               Average_Mode_Flag,  &
                               Rad,                & 
                               Rad_Min,            &  
                               Rad_Max)  
  ENDIF

  !-----------------------------------------------------------------
  ! Close the file and deallocate temporary memory
  !-----------------------------------------------------------------
  
  Result = CloseHimawariScene(Handle)    
  IF (Result /= 0) THEN
    call MESG( "Cannot close "//TRIM(Path_And_Filename), level = verb_lev % ERROR , color = 4 )
    Status = Sym%FAILURE
    RETURN
  ENDIF

  

  RETURN
  
END FUNCTION Get_Him_Raddata


! This function calls LibHimarawi to do the subsampling

FUNCTION Him_Rad_Subsample(Handle,                  &
                          Path_And_Filename,        &
                          Sentinels,                &
                          Start,                    &
                          Stride,                   &
                          Edge,                     &
                          Factor,                   &
                          Hsd_Meta,                 &
                          Rad_Output) RESULT(Status)

  TYPE(C_Ptr), INTENT(IN) :: Handle
  CHARACTER(LEN=2048), INTENT(IN)  :: Path_And_Filename
  TYPE(Hsd_Sentinels), INTENT(IN)  :: Sentinels
  INTEGER(KIND=INT4), DIMENSION(2), INTENT(IN) :: Start
  INTEGER(KIND=INT4), DIMENSION(2), INTENT(IN) :: Stride
  INTEGER(KIND=INT4), DIMENSION(2), INTENT(IN) :: Edge
  INTEGER(KIND=INT4), INTENT(IN)  :: Factor
  TYPE(Hsd_Metadata), INTENT(IN)  :: Hsd_Meta
  REAL(KIND=REAL4), DIMENSION(:,:), INTENT(INOUT) :: Rad_Output

  !LibHimawari source and destination types.
  TYPE(hsd_source) :: Src
  TYPE(hsd_destination)  :: Dst
  
  ! Internal variables
  INTEGER(KIND=INT4) :: Funct_Status, Result
 
 
  INTEGER(KIND=INT4) :: Status


  !-----------------------------------------------------------------
  ! Set the variables that control which parts of the image are
  ! read in
  !-----------------------------------------------------------------
  
  
  !Start line/element in native space. Needs to be in 0 based space.
  !Hence the subtraction
  Src%Column_Offset = (Start(1)*Factor)  - (1*Factor)
  Src%Line_Offset = (Start(2)*Factor)  - (1*Factor)
  
  !Amount of striding in native (NOT destination) space
  Src%Column_Stride = Stride(1)*Factor
  Src%Line_Stride = Stride(2)*Factor
  
  
  !How many lines/columns in *destination* space you want to transfer over
  ! i.e. the size of the output array, not the actual edge of the array
  Src%Columns = Edge(1)
  Src%Lines = Edge(2)
    
  !These relate to how C stores data, thus the x-dimension is stored in 
  ! Dst%Line_Increment. Dst%Column_Increment should be 1. It is best to
  ! be explicit in what is being used.
  Dst%Line_Increment = Src%Columns
  Dst%Column_Increment = 1
 
  !-----------------------------------------------------------------
  ! Read in the data
  !-----------------------------------------------------------------
    Funct_Status = copyHimawariSceneRadiances(Handle, Src, Dst, Rad_Output)
    IF (Funct_Status /= 0) THEN
       call MESG( "Cannot read radiance data from"//TRIM(Path_And_Filename), level = verb_lev % ERROR , color = 4 )
      Status = Sym%FAILURE
      !Result = CloseHimawariScene(Handle)    
      RETURN
    ENDIF


  RETURN


END FUNCTION Him_Rad_Subsample


! This function calls LibHimarawi to do the averaging.
! Returns the status of each step back to determin messaging
! and any other action

FUNCTION Him_Rad_Average(Handle,                  &
                          Path_And_Filename,        &
                          Sentinels,                &
                          Start,                    &
                          Stride,                   &
                          Edge,                     &
                          Factor,                   &
                          Hsd_Meta,                 &
                          Average_Mode_Flag,        &
                          Rad_Output,               &
                          Rad_Min,                  &
                          Rad_Max) RESULT(Status)

  TYPE(C_Ptr), INTENT(IN) :: Handle
  CHARACTER(LEN=2048), INTENT(IN)  :: Path_And_Filename
  TYPE(Hsd_Sentinels), INTENT(IN)  :: Sentinels
  INTEGER(KIND=INT4), DIMENSION(2), INTENT(IN) :: Start
  INTEGER(KIND=INT4), DIMENSION(2), INTENT(IN) :: Stride
  INTEGER(KIND=INT4), DIMENSION(2), INTENT(IN) :: Edge
  INTEGER(KIND=INT4), INTENT(IN)  :: Factor !Factor to go from native to IR
  TYPE(Hsd_Metadata), INTENT(IN)  :: Hsd_Meta
  INTEGER(KIND=INT4), INTENT(IN)  :: Average_Mode_Flag !which type of avg flag
  REAL(KIND=REAL4), DIMENSION(:,:), INTENT(INOUT) :: Rad_Output !Destination array
  REAL(KIND=REAL4), DIMENSION(:,:), INTENT(INOUT) :: Rad_Min !Destination avg min (Mode 2)
  REAL(KIND=REAL4), DIMENSION(:,:), INTENT(INOUT) :: Rad_Max !Destination avg max (Mode 2)

  !Internal variables
  REAL(KIND=REAL4), DIMENSION(:,:), Allocatable :: Rad_Native ! native res data array

  !LibHimawari types. Left names similar to test code.
  TYPE(hsd_source) :: Src ! source info for destination array
  TYPE(hsd_source) :: Src_Nat ! source info for native array
  TYPE(hsd_destination)  :: Nat
  
  ! Internal variables
  INTEGER(KIND=INT4) :: Funct_Status, Result
  INTEGER(KIND=INT4) :: Alloc_Status
  INTEGER(KIND=INT4) :: Status
  INTEGER(KIND=INT4) :: Output_Column_Increment, Output_Line_Increment
 
 
  ! Factor to reduce the data by. This is a function of the nat->IR and downsample factor (stored in Stride)
  INTEGER(KIND=INT4) :: Data_factor 

  !internal variables for loops, etc.
  integer:: nx, ny
  integer:: i, j, ni, nj, i1, i2, j1, j2, is, js, n
  integer:: X_Stride, Y_Stride

  !-----------------------------------------------------------------
  ! Set the variables that control which parts of the image are
  ! read in
  !-----------------------------------------------------------------
  
  !----------------------------------------------------------------- 
  !First thing is to get the radiances at *native* resolution
  !-----------------------------------------------------------------

  
  !Start line/element in native space. Needs to be in 0 based space.
  !Hence the subtraction
  Src_Nat%Column_Offset = (Start(1)*Factor)  - (1*Factor)
  Src_Nat%Line_Offset = (Start(2)*Factor)  - (1*Factor)
    
  !We need the entire native array before averaging
  Src_Nat%Column_Stride = 1
  Src_Nat%Line_Stride = 1
   
  !How many lines/columns in *destination* space you want to transfer over
  ! i.e. the size of the output array, not the actual edge of the array
  !Since this is in native resolution space, we have to *multiply* by Factor
  ! AND stride (mag factor) in the case of averaging
  
  Src_Nat%Columns = Edge(1)*(Stride(1)*Factor)
  Src_Nat%Lines = Edge(2)*(Stride(2)*Factor)
  
  !These relate to how C stores data, thus the x-dimension is stored in 
  ! Nat%Line_Increment. Nat%Column_Increment should be 1. It is best to
  ! be explicit in what is being used.
  Nat%Line_Increment = Src_Nat%Columns
  Nat%Column_Increment = 1

  !allocate the native radiance array
  
  ALLOCATE(Rad_Native(Src_Nat%Columns,Src_Nat%Lines),stat=Alloc_Status)

  IF (Alloc_Status /= 0) THEN
       call MESG( "Cannot allocate native radiance array", level = verb_lev % ERROR , color = 4 )
      Status = Sym%FAILURE
!      Result = CloseHimawariScene(Handle)    
      RETURN
  ENDIF
    
  !-----------------------------------------------------------------
  ! Read in the data at native resolution
  !-----------------------------------------------------------------
    Funct_Status = copyHimawariSceneRadiances(Handle, Src_Nat, Nat, Rad_Native)
    IF (Funct_Status /= 0) THEN
       call MESG( "Cannot read radiance data from"//TRIM(Path_And_Filename), level = verb_lev % ERROR , color = 4 )
      Status = Sym%FAILURE
!      Result = CloseHimawariScene(Handle)    
      RETURN
    ENDIF
    
  !-----------------------------------------------------------------
  ! Now we average
  !-----------------------------------------------------------------

  !How many lines/columns in *destination* space you want to transfer over
  ! i.e. the size of the output array, not the actual edge of the array
  ! Recall Edge is the size of the destination array
  Src%Columns = Edge(1)
  Src%Lines = Edge(2)

  Output_Line_Increment = Src%Columns
  Output_Column_Increment = 1
  
  ! Factor to reduce the data by. This is a function of the nat->IR and downsample factor (stored in Stride)
  
  
  Data_factor = Stride(1)* Factor
  
  !currently commented out until RKG and WCS can discuss differences in CLAVR-x/Libhim 
      
!  Status = averageHimawariFloatArray( &
!    Data_factor, &! Average factors
!    MISSING_VALUE_REAL4,  &! missing value sentinel to be expected on source and applied on destination
!    Rad_Native, & !data at native resolutio
!    Src_Nat%Lines,  & ! how many lines in Rad_Native
!    Src_Nat%Columns, & ! how many columns in Rad_Native
!    Nat%Line_Increment, &  ! how to get from line to line in Rad_Native
!    Nat%Column_Increment, &  ! how to get from column to column in Rad_Native
!    Src_Nat%Line_Offset, &     ! first_line offset within Rad_Native
!    Src_Nat%Column_Offset, &   ! first_column offset within Rad_Native
!    0,  &  ! 0, because we're forcing the reduction factor

!    Rad_Output,&  ! output array
!    Src%Lines, &  ! how many averaged lines and columns we expect to produce
!    Src%Columns, & 
!    Output_Line_Increment, & ! How to get from line to consecutive line in output array
!    Output_Column_Increment, & ! How to get from column to consecutive column in output array
!    0  &   ! ignored, we're using explicit Factor input
!)


! below here works - commented out 5/25 for testing

   ! size of the output array
   Src%Columns = Edge(1)
   Src%Lines = Edge(2)

   
   ! averaging factor in each direction
   X_Stride = Stride(1)*Factor
   Y_Stride = Stride(2)*Factor
   
   !number of pixels that are being averaged in the box
   n =  X_Stride*Y_Stride
  
 ! Loop over each pixel for averaging
 !--- i1:i2, j1:j2 = the native pixels used for mean/min/max

   
   !loop over rows, lines to get averaged data
   do i = 1, Src%Columns
      i1 = (i-1)*X_Stride + 1
      i2 = i1 + X_Stride - 1
      
      !loop over lines
      do j = 1,  Src%Lines
           j1 = (j-1)*Y_Stride + 1
           j2 = j1 + Y_Stride - 1
           j2 = min(j2, size(Rad_Native,2)) ! array bound test
                      
           Rad_Output(i,j) = sum(Rad_Native(i1:i2,j1:j2))/n

!          Calculate the min and max of the box used for averaging
           if (Average_Mode_Flag == 2) then
               Rad_Min(i,j) = minval(Rad_Native(i1:i2,j1:j2))
               Rad_Max(i,j) = maxval(Rad_Native(i1:i2,j1:j2))
           endif
      enddo
   enddo
    
   

  !-----------------------------------------------------------------
  ! Deallocate native array
  !-----------------------------------------------------------------

  IF (ALLOCATED(Rad_Native)) DEALLOCATE(Rad_Native,stat=Alloc_Status)
  
  IF (Alloc_Status /= 0) THEN
       call MESG( "Cannot deallocate native radiance array", level = verb_lev % ERROR , color = 4 )
      Status = Sym%FAILURE
!      Result = CloseHimawariScene(Handle)    
      RETURN
  ENDIF


  RETURN


END FUNCTION Him_Rad_Average



  
END MODULE AHI_HSD_READER
