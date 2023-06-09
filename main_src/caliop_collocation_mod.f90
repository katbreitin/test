! $Id: caliop_collocation_mod.f90 3445 2019-08-08 15:33:08Z dbotambekov $
!--------------------------------------------------------------------------------------
! 
! This module is: 
!  1. Searching for corresponding CALIOP-PASSIVE collocation file
!  2. Reads passive pixel index
!  3. Marks all pixels as Geo%Space_Mask except +/- N pixels along CALIOP path
!
! Author: Denis Botambekov
!
! Date: 2018-05-08
!
!--------------------------------------------------------------------------------------

module CALIOP_COLLOCATION_MOD

  use FILE_TOOLS, only: &
        FILE_SEARCH

  use CLAVRX_MESSAGE_MOD

  use CX_HDF4_MOD, only: &
        OPEN_FILE_HDF_READ &
      , HDF_SDS_READER &
      , CLOSE_FILE_HDF_READ

  use PIXEL_COMMON_MOD, only: &
        Caliop_Dir &
      , Skip_L1b_File_Flag &
      , Geo &
      , Image &
      , Ch &
      , Sensor &
      , Caliop_Num_Cld_Layers &
      , Caliop_Cod & 
      , Caliop_Cld_Height

  use CONSTANTS_MOD, only: &
        int4 &
      , real4 &
      , Nchan_Clavrx &
      , SOLAR_OBS_TYPE &
      , THERMAL_OBS_TYPE &
      , MIXED_OBS_TYPE &
      , Sym

  implicit none

  public CALIOP_COLLOCATION

  character(len=15), parameter:: CALIOP_PROMPT="CALIOP_MODULE:"

  contains

!--------------------------------------------------------------------------------------
subroutine CALIOP_COLLOCATION(Seg_Idx)

  integer, intent(in):: Seg_Idx
  character(len=1020) :: File_Name_Search
  character(len=1020) :: Caliop_File
  character(len=1020), dimension(:), pointer:: Files
  integer(kind=int4) :: Num_Files
  integer(kind=int4) :: Status
  integer(kind=int4) :: Sd_Id
  integer(kind=int4) :: i, Chan_I
  integer(kind=int4) :: Left_Limit, Right_Limit
  integer(kind=int4), parameter :: N = 20
  integer(kind=int4), dimension (1) :: Start_1d, Stride_1d, Edge_1d
  real(kind=real4), dimension (:), allocatable :: Pixel_Element_Idx
  real(kind=real4), dimension (:), allocatable :: Caliop_Num_Cld_Layers_1D
  real(kind=real4), dimension (:), allocatable :: Caliop_Cod_1D
  real(kind=real4), dimension (:), allocatable :: Caliop_Cld_Height_1D
  real(kind=real4),parameter :: Missing = -999.0

  ! --- create string to search caliop file
  file_name_search = '*'//image % time_start % date_string('yyyy_doy_hhmm')//'*.calipso.hdf'

  ! --- search file
  Files => FILE_SEARCH(trim(Caliop_Dir)//'/',trim(File_Name_Search),count=Num_Files)

  if (Num_Files == 0) then
    File_Name_Search = '*A'//image % time_start % date_string('yyyy_doy.hhmm')//'*.calipso.hdf'
    Files => FILE_SEARCH(trim(Caliop_Dir)//'/',trim(File_Name_Search),count=Num_Files)
  endif

  if (Num_Files == 0) then
    File_Name_Search = '*A'//image % time_start % date_string('yyyy_doy_hhmm')//'*.calipso.hdf'
    Files => FILE_SEARCH(trim(Caliop_Dir)//'/',trim(File_Name_Search),count=Num_Files)
  endif

  if (Num_Files == 0) then
    File_Name_Search = '*d'//image % time_start % date_string('yyyymmdd')//'_t' &
            //image % time_start % date_string('hhmm')//'*.calipso.hdf'
    Files => FILE_SEARCH(trim(Caliop_Dir)//'/',trim(File_Name_Search),count=Num_Files)
  endif

  if (Num_Files == 0) then
    File_Name_Search = '*s'//image % time_start % date_string('yyyydoyhhmm')//'*.calipso.hdf'
    Files => FILE_SEARCH(trim(Caliop_Dir)//'/',trim(File_Name_Search),count=Num_Files)
  endif

  ! --- if not found or more than 1 exit
  if (Num_Files == 0 .or. Num_Files > 1) then
    call MESG (CALIOP_PROMPT//'CALIOP Collocation File Not Found, Skipping',Level = Verb_Lev % WARNING)
    Skip_L1b_File_Flag = .true.
    return
  endif

  ! --- file found save it
  Caliop_File = Files(1)

  if (Seg_Idx == 1) &
    call MESG (CALIOP_PROMPT//'CALIOP File Found '//trim(Caliop_File), Level = Verb_Lev % DEFAULT)

  ! --- get start, end, stride
  Start_1d(1) = (Seg_Idx -1) * Image%Number_Of_Lines_Per_Segment
  Stride_1d(1) = 1
  Edge_1d(1) = min(Start_1d(1) + Image%Number_Of_Lines_Per_Segment, &
                                 Image%Number_Of_Lines) - Start_1d(1)

  ! --- allocate 
  allocate(Pixel_Element_Idx(Edge_1d(1)))
  allocate(Caliop_Num_Cld_Layers_1D(Edge_1d(1)))
  allocate(Caliop_Cod_1D(Edge_1d(1)))
  allocate(Caliop_Cld_Height_1D(Edge_1d(1)))

  ! --- read file
  Status = 0
  Status = OPEN_FILE_HDF_READ(trim(Caliop_Dir)//trim(Caliop_File), Sd_Id)
  Status = HDF_SDS_READER(Sd_Id, 'passive_pixel_element', Start_1d, Stride_1d, &
                          Edge_1d, Pixel_Element_Idx) + Status
  Status = HDF_SDS_READER(Sd_Id, 'closest_calipso_number_of_cloud_layers', Start_1d, &
                          Stride_1d, Edge_1d, Caliop_Num_Cld_Layers_1D) + Status
  Status = HDF_SDS_READER(Sd_Id, 'closest_calipso_cod', Start_1d, &
                          Stride_1d, Edge_1d, Caliop_Cod_1D) + Status
  Status = HDF_SDS_READER(Sd_Id, 'closest_calipso_top_alt_eff_geo', Start_1d, &
                          Stride_1d, Edge_1d, Caliop_Cld_Height_1D) + Status

  Status = CLOSE_FILE_HDF_READ(Sd_Id, trim(Caliop_Dir)) + Status

  ! --- convert height from km to meters
  where(Caliop_Cld_Height_1D .ne. Missing)
     Caliop_Cld_Height_1D = Caliop_Cld_Height_1D * 1000.0
  endwhere


  ! --- set all Geo%Space_Mask pixels except path
  do i = 1, Edge_1d(1)

    ! --- set all to SPACE if Missing value
    if (Pixel_Element_Idx(i) == Missing .or. Pixel_Element_Idx(i) .lt. 0) then
      Geo%Space_Mask(1:Image%Number_Of_Elements,i) = .true.
      ! - loop through the channels and set them to missing
      do Chan_I = 1, Nchan_Clavrx
         if (Sensor%Chan_On_Flag_Default(Chan_I) .eq. sym%YES) then
            ! - Ref
            if (ch(Chan_I)%Obs_Type == SOLAR_OBS_TYPE) then
               ch(Chan_I)%Ref_Toa(:,i) = Missing
            endif
            ! - BT
            if (ch(Chan_I)%Obs_Type == THERMAL_OBS_TYPE .or. &
                ch(Chan_I)%Obs_Type == MIXED_OBS_TYPE) then
               ch(Chan_I)%Rad_Toa(:,i) = Missing
               ch(Chan_I)%Bt_Toa(:,i) = Missing
            endif
         endif
      enddo
      cycle
    endif

    ! - avoid going out of array size
    Left_Limit = max(1.,Pixel_Element_Idx(i)-N)
    Right_Limit = min(Pixel_Element_Idx(i)+N,real(Image%Number_Of_Elements))
    Geo%Space_Mask(1:Left_Limit,i) = .TRUE.
    Geo%Space_Mask(Right_Limit:Image%Number_Of_Elements,i) = .TRUE.

    ! - loop through the channels and set them to missing
    do Chan_I = 1, Nchan_Clavrx
       if (Sensor%Chan_On_Flag_Default(Chan_I) .eq. sym%YES) then
          ! - Ref
          if (ch(Chan_I)%Obs_Type == SOLAR_OBS_TYPE) then
             ch(Chan_I)%Ref_Toa(1:Left_Limit,i) = Missing
             ch(Chan_I)%Ref_Toa(Right_Limit:Image%Number_Of_Elements,i) = Missing
          endif
          ! - BT
          if (ch(Chan_I)%Obs_Type == THERMAL_OBS_TYPE .or. &
              ch(Chan_I)%Obs_Type == MIXED_OBS_TYPE) then
             ch(Chan_I)%Rad_Toa(1:Left_Limit,i) = Missing
             ch(Chan_I)%Rad_Toa(Right_Limit:Image%Number_Of_Elements,i) = Missing
             ch(Chan_I)%Bt_Toa(1:Left_Limit,i) = Missing
             ch(Chan_I)%Bt_Toa(Right_Limit:Image%Number_Of_Elements,i) = Missing
          endif
       endif
    enddo
             
    ! --- save caliop data to global variables
    Caliop_Num_Cld_Layers(int(Pixel_Element_Idx(i)+1),i) = Caliop_Num_Cld_Layers_1D(i)
    Caliop_Cod(int(Pixel_Element_Idx(i)+1),i) = Caliop_Cod_1D(i)
    Caliop_Cld_Height(int(Pixel_Element_Idx(i)+1),i) = Caliop_Cld_Height_1D(i)
   
  enddo

  ! --- deallocate and null
  deallocate(Pixel_Element_Idx)
  deallocate(Caliop_Num_Cld_Layers_1D)
  deallocate(Caliop_Cod_1D)
  deallocate(Caliop_Cld_Height_1D)
  Files => null() 

end subroutine CALIOP_COLLOCATION

!--------------------------------------------------------------------------------------

end module CALIOP_COLLOCATION_MOD

