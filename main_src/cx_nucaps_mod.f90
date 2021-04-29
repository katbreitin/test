!$Id: clavrx_static_nav_module.f90 2903 2018-08-07 22:22:23Z heidinger $
!------------------------------------------------------------------------------
!
! Module to Handle NUCAPS information in CLAVR-x
!------------------------------------------------------------------------------
module CX_NUCAPS_MOD

!  use AWG_CLOUD_HEIGHT, only: &
!    MEAN_SMOOTH2

  use CLAVRX_MESSAGE_MOD
 
  use CLOUD_HEIGHT_ROUTINES, only: &
    COMPUTE_BOX_WIDTH

  use CONSTANTS_MOD, only: &
    int4 &
    , int1 &
    , sym &
    , Missing_Value_Real4 &
    , real4
    
  use CX_NETCDF4_MOD, only:  &
    close_netcdf &
    , open_netcdf &
    , read_netcdf

  use NWP_COMMON_MOD, only: &
    KNOWING_P_COMPUTE_T_Z_NWP
  
  use PIXEL_COMMON_MOD, only: &
    Image &
    , NUCAPS &
    , Bad_Pixel_Mask &
    , NWP_PIX &
    , Diag_Pix_Array_1 &
    , Diag_Pix_Array_2 &
    , Diag_Pix_Array_3 &
    , Cld_Type &
    , Sensor &
    , Nav &
    , Ch

  use kdtree2_module
  !use time_kdtree

  implicit none
  public::  VIIRS_NUCAPS
  private:: SETUP_NUCAPS
  private:: READ_NUCAPS
  private:: MEAN_SMOOTH2
  private:: KD_TREE_INTERP
  public:: CONVERT_SMOOTH_NUCAPS_TEMP
  public:: KDTREE_COARSE2FINE

  character(len=120), parameter :: EXE_PROMPT_NAV = "CLAVR-x NUCAPS Module >> "

  integer, save:: ncid_static
  integer, save:: nc_varid_lat, nc_varid_lon, nc_varid_senzen, nc_varid_senaz

  character(len=1020), private, save:: NUCAPS_Full_File

contains

  subroutine VIIRS_NUCAPS(Seg_Idx)
      integer, intent(in):: Seg_Idx
      integer:: nucaps_on_flag

      call SETUP_NUCAPS(Seg_Idx,nucaps_on_flag)

      if (Seg_Idx == 1) then
         if (nucaps_on_flag == 0) then
            call MESG (trim(EXE_PROMPT_NAV)//trim('NUCAPS File Not Found, Skipping'),Level = Verb_Lev % WARNING)
         return
      endif

      endif

      call READ_NUCAPS(Seg_Idx,nucaps_on_flag)

  end subroutine VIIRS_NUCAPS

  !------------------------------------------------------------------------------------
  !  Determine Name and try to open nucaps file
  !------------------------------------------------------------------------------------
  subroutine SETUP_NUCAPS(Seg_Idx,NUCAPS_On_Flag)

    use FILE_TOOLS, only: file_search

    integer, intent(in):: Seg_Idx
    integer, intent(out):: NUCAPS_On_Flag
    character(len=1020), dimension(:), pointer:: Files
    integer:: Num_Files
    character(len=1020) :: File_Name_Search
    character(len=4) :: Year_String, Time_String
    character(len=2) :: Month_String, Dom_String

    NUCAPS_On_Flag  = 0

    Year_String = Image%Level1b_Name(12:15)
    Month_String = Image%Level1b_Name(16:17)
    Dom_String = Image%Level1b_Name(18:19)
    Time_String = Image%Level1b_Name(22:25)

    File_Name_Search = 'GMTCO_*'//'_d'//Year_String//Month_String//Dom_String//'_t'//Time_String//'*.nucaps.nc'
    Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(File_Name_Search),count=Num_Files)

    ! --- if not found or more than 1 exit
    if (Num_Files == 0 .or. Num_Files > 1) then
      NUCAPS_On_Flag = 0
      return
    endif
  
    NUCAPS_Full_File = Files(1)

    Files => null()

    NUCAPS_On_Flag = 1

  end subroutine SETUP_NUCAPS

  !------------------------------------------------------------------------------------
  !  Read NUCAPS Data for a Segment
  !------------------------------------------------------------------------------------
  subroutine READ_NUCAPS(Seg_Idx,NUCAPS_On_Flag)

    integer, intent(in):: Seg_Idx
    integer, intent(in):: NUCAPS_on_Flag
    integer:: Ncid_nucaps
    integer:: Elem_Start_Segment, Elem_Count_Segment
    integer:: Line_Start_Segment, Line_Count_Segment
    integer:: Line_End_File, Line_End_Segment
    integer(kind=int4), dimension(2):: Sds_Dims, Sds_Stride, Sds_Start,Sds_Edges

    if (NUCAPS_on_Flag == 0) return

    Elem_Start_Segment = 1
    Elem_Count_Segment = Image%Number_Of_Elements

    Line_Start_Segment = 1 + (Image%Segment_Number-1)*Image%Number_Of_Lines_Per_Segment*Image%Y_Stride
    Line_End_Segment = Line_Start_Segment + Image%Number_Of_Lines_Per_Segment*Image%Y_Stride
    Line_End_File = Image%Number_Of_Lines
    Line_End_Segment = min(Line_End_Segment, Line_End_File)
    Line_Count_Segment = (Line_End_Segment - Line_Start_Segment)/Image%Y_Stride

    Sds_Stride = (/Image%X_Stride,Image%Y_Stride/)
    Sds_Start = (/Elem_Start_Segment, Line_Start_Segment/)
    Sds_Edges = (/Elem_Count_Segment, Line_Count_Segment/)

    !--- open file for read
    call OPEN_NETCDF(trim(Image%Level1b_Path)//trim(NUCAPS_Full_File), Ncid_Nucaps)

    !--- if file is unreadable, exit
    if (Ncid_Nucaps <= 0) then
      call MESG(trim(EXE_PROMPT_NAV)//trim("NUCAPS File Could Not Be Opened: "),Level = Verb_Lev % WARNING)
      return
    endif

    call READ_NETCDF(ncid_nucaps, Sds_Start, Sds_Stride, Sds_Edges, "cloud_top_pressure_nucaps_layer1",NUCAPS%Cld_Press_Layer1)

    call READ_NETCDF(ncid_nucaps, Sds_Start, Sds_Stride, Sds_Edges, "cloud_top_pressure_nucaps_layer2",NUCAPS%Cld_Press_Layer2)
 
    call READ_NETCDF(ncid_nucaps, Sds_Start, Sds_Stride, Sds_Edges, "cloud_top_fraction_nucaps_layer1",NUCAPS%Cld_Fraction_Layer1)

    call READ_NETCDF(ncid_nucaps, Sds_Start, Sds_Stride, Sds_Edges, "cloud_top_fraction_nucaps_layer2",NUCAPS%Cld_Fraction_Layer2)

    call MESG(trim(EXE_PROMPT_NAV)//trim("NUCAPS File read in successfully "))

    !--- close file
    call CLOSE_NETCDF(Ncid_Nucaps)

  ! Moved value assignment here to separate from interpolation
    NUCAPS%Cld_Press = NUCAPS%Cld_Press_Layer1
    NUCAPS%Cld_Fraction = NUCAPS%Cld_Fraction_Layer1

  end subroutine READ_NUCAPS
  !------------------------------------------------------------------------------------
  !  Convert Pres to Temp and Smooth
  !------------------------------------------------------------------------------------
  subroutine CONVERT_SMOOTH_NUCAPS_TEMP()

    integer:: Line_Idx
    integer:: Elem_Idx
    integer:: Num_Elements
    integer:: Num_Lines
    integer:: Nwp_Lon_Idx
    integer:: Nwp_Lat_Idx
    integer:: Level_Idx_Temp

    integer(kind=int1), dimension(:,:), allocatable:: Mask1
    integer(kind=int1), dimension(:,:), allocatable:: Mask2
    real(kind=real4):: Sensor_Resolution_KM
    real(kind=real4):: Cirrus_Resolution_KM
    integer (kind=int4):: Box_Half_Width_Cirrus
    integer (kind=int4), parameter:: CIRRUS_BOX_WIDTH_NUCAPS_KM = 400 !200
    integer (kind=int4), parameter:: N_Idx_Found = 6
    logical :: Use_Kdtree_Flag

    Use_Kdtree_Flag = .true.
    Num_Elements = Image%Number_Of_Elements
    Num_Lines = Image%Number_Of_Lines_Per_Segment

    if (.not. allocated(Mask1)) allocate(Mask1(Num_Elements,Num_Lines))
    if (.not. allocated(Mask2)) allocate(Mask2(Num_Elements,Num_Lines))

    !--- Determine cirrus spatial interpolation box width.
    !--- Sensor_Resolution_KM = Sensor%Spatial_Resolution_Meters/1000.0
    Sensor_Resolution_KM = 0.75
    call COMPUTE_BOX_WIDTH(Sensor_Resolution_KM, CIRRUS_BOX_WIDTH_NUCAPS_KM, Box_Half_Width_Cirrus)
    print *, "CIRRUS_BOX_WIDTH_NUCAPS_KM = ", CIRRUS_BOX_WIDTH_NUCAPS_KM, Box_Half_Width_Cirrus

    do Elem_Idx = 1, Num_Elements
      do Line_Idx = 1, Num_Lines

        !--- Skip bad pixels.
        if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

        !--- Indice Aliases
        Nwp_Lon_Idx = NWP_PIX%I_Nwp(Elem_Idx,Line_Idx)
        Nwp_Lat_Idx = NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)

        !--- Check if indices are valid.
        if (Nwp_Lon_Idx < 0 .or. Nwp_Lat_Idx < 0) cycle


        !----------------------------------------------------------------------
        !  Method 1: Use Layer 1 Exclusively
        !----------------------------------------------------------------------
        !NUCAPS%Cld_Press(Elem_Idx,Line_Idx) = NUCAPS%Cld_Press_Layer1(Elem_Idx,Line_Idx)
        !NUCAPS%Cld_Fraction(Elem_Idx,Line_Idx) = NUCAPS%Cld_Fraction_Layer1(Elem_Idx,Line_Idx) 
        !NUCAPS%Cld_Fraction(Elem_Idx,Line_Idx) = &
        !   max(NUCAPS%Cld_Fraction_Layer1(Elem_Idx,Line_Idx),NUCAPS%Cld_Fraction_Layer2(Elem_Idx,Line_Idx))

        !----------------------------------------------------------------------
        !  Method 2: Weight Layer 1 and Layer 2
        !----------------------------------------------------------------------
!       if (NUCAPS%Cld_Press_Layer1(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4 .and. &
!           NUCAPS%Cld_Press_Layer2(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4 .and. &
!           NUCAPS%Cld_Fraction_Layer1(Elem_Idx,Line_Idx) > 0.00 .and. &
!           NUCAPS%Cld_Fraction_Layer2(Elem_Idx,Line_Idx) > 0.00) then

!       NUCAPS%Cld_Press(Elem_Idx,Line_Idx) = &
!             (NUCAPS%Cld_Fraction_Layer1(Elem_Idx,Line_Idx)*NUCAPS%Cld_Press_Layer1(Elem_Idx,Line_Idx) +  &
!              NUCAPS%Cld_Fraction_Layer2(Elem_Idx,Line_Idx)*NUCAPS%Cld_Press_Layer2(Elem_Idx,Line_Idx)) / &
!              (NUCAPS%Cld_Fraction_Layer1(Elem_Idx,Line_Idx) + NUCAPS%Cld_Fraction_Layer2(Elem_Idx,Line_Idx))
!       NUCAPS%Cld_Fraction(Elem_Idx,Line_Idx) = 0.5*(NUCAPS%Cld_Fraction_Layer1(Elem_Idx,Line_Idx) +   &
!                                                             NUCAPS%Cld_Fraction_Layer1(Elem_Idx,Line_Idx))
!       endif

        !----------------------------------------------------------------------
        !--- Check for valid pressure and fraction from NUCAPS. 
        !----------------------------------------------------------------------
        if (NUCAPS%Cld_Press(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
        if (NUCAPS%Cld_Fraction(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

        !----------------------------------------------------------------------
        !--- Calculate temperature and height from pressure and NWP.
        !----------------------------------------------------------------------
!Diag_Pix_Array_1(Elem_Idx,Line_Idx) = NUCAPS%Cld_Fraction(Elem_Idx,Line_Idx)
        call KNOWING_P_COMPUTE_T_Z_NWP(Nwp_Lon_Idx, &
                                       Nwp_Lat_Idx, &
                                       NUCAPS%Cld_Press(Elem_Idx,Line_Idx), &
                                       NUCAPS%Cld_Temp(Elem_Idx,Line_Idx), &
                                       NUCAPS%Cld_Height(Elem_Idx,Line_Idx), &
                                       Level_Idx_Temp)
      enddo
    enddo

    !------------------------------------------------------------------
    !--- Adjust Cloud Type
    !------------------------------------------------------------------
!Diag_Pix_Array_3 = Cld_Type
    where( NUCAPS%Cld_Temp < 250.0 .and. &
           NUCAPS%Cld_Temp /= Missing_Value_Real4 .and. &
           NUCAPS%Cld_Fraction >= 0.2 .and.  &
           Cld_Type == sym%SUPERCOOLED_TYPE)
       Cld_Type = sym%CIRRUS_TYPE
    end where

    where( NUCAPS%Cld_Fraction < 0.1 .and. &
           NUCAPS%Cld_Fraction /= Missing_Value_Real4 .and. &
           (Cld_Type == sym%OVERLAP_TYPE .or. Cld_Type == sym%CIRRUS_TYPE))
       Cld_Type = sym%WATER_TYPE
    end where

    !------------------------------------------------------------------
    !--- Make source mask 
    !------------------------------------------------------------------
    Mask1 = 0_int1
    where( NUCAPS%Cld_Temp < 250.0 .and. &
           NUCAPS%Cld_Temp /= Missing_Value_Real4 .and. &
           NUCAPS%Cld_Fraction >= 0.4)
      Mask1 = 1_int1
    end where

    where( NUCAPS%Cld_Temp < 250.0 .and. &
           NUCAPS%Cld_Temp /= Missing_Value_Real4 .and. &
           NUCAPS%Cld_Fraction >= 0.1 .and. &
           (Cld_Type == sym%CIRRUS_TYPE .or. &
            Cld_Type == sym%OPAQUE_ICE_TYPE .or. &
            Cld_Type == sym%OVERSHOOTING_TYPE .or. &
            Cld_Type == sym%OVERLAP_TYPE))
      Mask1 = 1_int1
    end where

    !Diag_Pix_Array_1 = Mask1
!   where(NUCAPS%Cld_Fraction < 0.1)
!       Diag_Pix_Array_2 = MISSING_VALUE_REAL4
!   endwhere

    !------------------------------------------------------------------
    ! Make target mask.
    !------------------------------------------------------------------
    Mask2 = 0_int1
    where(Cld_Type == sym%CIRRUS_TYPE .or. &
          Cld_Type == sym%OVERLAP_TYPE .or. &
          Cld_Type == sym%OPAQUE_ICE_TYPE .or. &
          Cld_Type == sym%OVERSHOOTING_TYPE .or. &
          Cld_Type == sym%SUPERCOOLED_TYPE)
      Mask2 = 1_int1
    end where
    !Diag_Pix_Array_2 = Mask2

    ! use kd-tree to fill values
    if (Use_Kdtree_Flag) then
        call KD_TREE_INTERP(Mask1, Mask2, ch(31)%BT_Toa, ch(32)%BT_Toa, Nav%Lat, Nav%Lon, Num_Elements, Num_Lines, N_Idx_Found, &
                            NUCAPS%Cld_Temp, NUCAPS%Cld_Temp_Smoothed)
!                            ch(29)%BT_Toa, NUCAPS%Cld_Temp_Smoothed)
    endif

    if (.false.) then
    !--- Call smoothing routine.
       call MEAN_SMOOTH2(Mask1, Mask2, Missing_Value_Real4, 1, 1, &
                         Box_Half_Width_Cirrus, Num_Elements, Num_Lines, &
                         NUCAPS%Cld_Temp, &
                         NUCAPS%Cld_Temp_Smoothed)                      
    endif

    if (allocated(Mask1)) deallocate(Mask1)
    if (allocated(Mask2)) deallocate(Mask2)

    !Diag_Pix_Array_3 = NUCAPS%Cld_Temp_Smoothed

  end subroutine CONVERT_SMOOTH_NUCAPS_TEMP
  !-----------------------------------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------------------------------
  subroutine MEAN_SMOOTH2(Mask_In,Mask_Out,Missing,di,dj,N,Num_Elements, Num_Lines, Z_In,Z_Out)

     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_In
     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_Out
     real (kind=real4), intent(in):: Missing
     integer (kind=int4), intent(in):: di
     integer (kind=int4), intent(in):: dj

     integer (kind=int4), intent(in):: N
     integer (kind=int4), intent(in):: Num_Elements
     integer (kind=int4), intent(in):: Num_lines
     real (kind=real4), intent(in), dimension(:,:), target:: Z_In
     real (kind=real4), intent(out), dimension(:,:):: Z_Out
     integer (kind=int4), dimension(size(Mask_In,1),size(Mask_In,2)):: Count_Out
     integer:: i
     integer:: j
     real:: Count_Temporary
     real, pointer, dimension(:,:):: Z_In_Sub
     integer (kind=int1), pointer, dimension(:,:):: Mask_In_Sub
     integer (kind=int1), pointer, dimension(:,:):: Mask_Out_Sub
     integer:: i1,i2,j1,j2


     Z_Out = 0.0
     Count_Out = 0

     do j = 1 + dj, Num_Lines-dj, dj + 1

        j1 = min(Num_Lines,max(1,j - N))
        j2 = min(Num_Lines,max(1,j + N))

        do i = 1 + di, Num_Elements - di, di + 1

           if (Mask_In(i,j) == 0) cycle
           if (Z_out(i,j) > 0) cycle

           i1 = min(Num_Elements,max(1,i - N))
           i2 = min(Num_ELements,max(1,i + N))

           Mask_In_Sub => Mask_In(i1:i2,j1:j2)
           Mask_Out_Sub => Mask_Out(i1:i2,j1:j2)

           if (sum(Mask_Out_Sub) == 0) cycle

           Z_In_Sub => Z_In(i1:i2,j1:j2)
           Count_Temporary = sum(real(Mask_In_Sub))

           Z_Out(i1:i2,j1:j2) = Z_Out(i1:i2,j1:j2) + sum(Z_In_Sub*Mask_In_Sub) / Count_Temporary
           Count_Out(i1:i2,j1:j2) = Count_Out(i1:i2,j1:j2) + 1

           Mask_In_Sub => null()
           Mask_Out_Sub => null()
           Z_In_Sub => null()

        enddo
     enddo

     !--- make mean value
     where(Count_Out > 0)
         Z_Out = Z_Out / Count_Out
     endwhere
     where(Count_Out == 0)
         Z_Out = Missing
     endwhere

     !--- only values are missing where output mask is 0
     where(Mask_Out == 0)
       Z_Out = Missing
     endwhere

  end subroutine MEAN_SMOOTH2
 
  !-----------------------------------------------------------------------------------------------------
  subroutine KD_TREE_INTERP(Mask_In,Mask_Out,pred_var1,pred_var2,pred_var3,pred_var4,Num_Elements, Num_Lines, nnbrute, Z_In,Z_Out)

     ! This function performs kd tree search. Z_In is the original data, and
     ! Z_Out is the output where indices flagged by Mask_Output are reassinged using average of values flagged by Mask_In.
     ! Currently values where Mask_In is set are not changed. nnbrute is the number of found closet indices for each search

     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_In
     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_Out
     real (kind=real4), intent(in), dimension(:,:), target:: pred_var1
     real (kind=real4), intent(in), dimension(:,:), target:: pred_var2
     real (kind=real4), intent(in), dimension(:,:), target:: pred_var3
     real (kind=real4), intent(in), dimension(:,:), target:: pred_var4
     integer (kind=int4), intent(in):: Num_Elements
     integer (kind=int4), intent(in):: Num_lines
     integer (kind=int4), intent(in):: nnbrute
     real (kind=real4), intent(in), dimension(:,:), target:: Z_In
     real (kind=real4), intent(out), dimension(:,:):: Z_Out

     ! define local variables, each KD-tree predictor is 1-D array
     type(kdtree2), pointer :: tree
     integer(kind=int4):: n_training, n_query
     integer:: i,j
     integer(kind=int4), dimension(:), allocatable:: ind_training, ind_query
     real(kind=real4), dimension(:), allocatable:: predictor_1
     real(kind=real4), dimension(:), allocatable:: predictor_2
     real(kind=real4), dimension(:), allocatable:: predictor_3
     real(kind=real4), dimension(:), allocatable:: predictor_4
     real(kind=real4), dimension(:), allocatable:: nucaps_temp_1d
     real(kind=real4), dimension(:), allocatable:: Z_Out_1d
     real(kind=real4), dimension(:,:), allocatable:: my_array
     real(kind=real4), dimension(:), allocatable:: query_vec
     type(kdtree2_result), dimension(:), allocatable :: results1

     integer:: Line_Idx
     integer:: Elem_Idx

     if (.not. allocated(predictor_1)) allocate(predictor_1(Num_Elements*Num_Lines))
     if (.not. allocated(predictor_2)) allocate(predictor_2(Num_Elements*Num_Lines))
     if (.not. allocated(predictor_3)) allocate(predictor_3(Num_Elements*Num_Lines))
     if (.not. allocated(predictor_4)) allocate(predictor_4(Num_Elements*Num_Lines))
     if (.not. allocated(nucaps_temp_1d)) allocate(nucaps_temp_1d(Num_Elements*Num_Lines))
     if (.not. allocated(Z_Out_1d)) allocate(Z_Out_1d(Num_Elements*Num_Lines))

     Z_Out_1d = MISSING_VALUE_REAL4

     ! find indices of training and query variables
     ! values at training indices (Mask_In) can be computed again but original values are kept here
     n_training = COUNT(Mask_In == 1)
     n_query = COUNT(Mask_Out == 1) - Count(Mask_Out == 1 .and. Mask_In == 1)

     ! perform kd-tree only if both training and query indices are found
     if (n_training > 0 .and. n_query > 0) then
        ! convert 2d to 1d array
        predictor_1 = reshape(pred_var1,(/Num_Elements*Num_Lines/))
        predictor_2 = reshape(pred_var2,(/Num_Elements*Num_Lines/))
        predictor_3 = reshape(pred_var3,(/Num_Elements*Num_Lines/))
        predictor_4 = reshape(pred_var4,(/Num_Elements*Num_Lines/))
        nucaps_temp_1d = reshape(Z_In,(/Num_Elements*Num_Lines/))

        if (.not. allocated(ind_training)) allocate(ind_training(n_training))
        if (.not. allocated(ind_query)) allocate(ind_query(n_query))

        ! search for training and query indices and stored
        i = 1
        j = 1
        do Line_Idx = 1, Num_Lines
        do Elem_Idx = 1, Num_Elements
          if (Mask_In(Elem_Idx,Line_Idx) == 1 .and. pred_var1(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var2(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var3(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var4(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              ind_training(i) = Elem_Idx + Num_Elements*(Line_Idx-1)
              i = i+1
          endif
          if (Mask_Out(Elem_Idx,Line_Idx) == 1 .and. Mask_In(Elem_Idx,Line_Idx) == 0 &
              .and. pred_var1(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var2(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var3(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var4(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              ind_query(j) = Elem_Idx + Num_Elements*(Line_Idx-1)
              j = j+1
          endif
          enddo
        enddo

        if (minval(nucaps_temp_1d(ind_training)) == Missing_Value_Real4) print *, 'nuacaps training is not right'

        if (.not. allocated(my_array)) allocate(my_array(4,n_training))
        if (.not. allocated(query_vec)) allocate(query_vec(4))

        ! assign training data
        my_array(1,:) = predictor_1(ind_training)
        my_array(2,:) = predictor_2(ind_training)
        my_array(3,:) = predictor_3(ind_training)
        my_array(4,:) = predictor_4(ind_training)

        ! create tree, set sort to true slows it down but the output indices are
        ! ordered from closet to farthest; set rearrange as true speeds searches
        ! but requires extra memory
        tree => kdtree2_create(my_array,sort=.true.,rearrange=.true.)  ! this is how you create a tree.

        if (.not. allocated(results1)) allocate(results1(nnbrute))

        ! set 1d output variable values at training indices
        do i = 1,n_training
           Z_Out_1d(ind_training(i)) = nucaps_temp_1d(ind_training(i))
        enddo

        ! perform tree search for each query index
        do i = 1,n_query
           query_vec(1) = predictor_1(ind_query(i))
           query_vec(2) = predictor_2(ind_query(i))
           query_vec(3) = predictor_3(ind_query(i))
           query_vec(4) = predictor_4(ind_query(i))

           ! results1 has both indices and distances 
           call kdtree2_n_nearest(tp=tree, qv=query_vec, nn = nnbrute, results=results1)
           
           ! average values for the all found indices
           Z_Out_1d(ind_query(i)) = sum(nucaps_temp_1d(ind_training(results1%idx)))/size(results1%idx)
        enddo

        ! destroy tree and release memory
        call kdtree2_destroy(tree)

     endif

     ! Change 1d array back to 2d; if no kdtree search is performed, array is empty
     Z_Out = reshape(Z_Out_1d,(/Num_Elements,Num_Lines/))

     if (allocated(predictor_1)) deallocate(predictor_1)
     if (allocated(predictor_2)) deallocate(predictor_2)
     if (allocated(predictor_3)) deallocate(predictor_3)
     if (allocated(predictor_4)) deallocate(predictor_4)
     if (allocated(nucaps_temp_1d)) deallocate(nucaps_temp_1d)
     if (allocated(Z_Out_1d)) deallocate(Z_Out_1d)
     if (allocated(query_vec)) deallocate(query_vec)
     if (allocated(my_array)) deallocate(my_array)
     if (allocated(ind_training)) deallocate(ind_training)
     if (allocated(ind_query)) deallocate(ind_query)
     if (allocated(results1)) deallocate(results1)

  end subroutine KD_TREE_INTERP

  subroutine KDTREE_COARSE2FINE(pred_var_coarse,pred_var_fine,Num_Elements_coarse,Num_Lines_coarse, &
                                Num_Elements_fine,Num_Lines_fine,nnbrute,Z_In_coarse,Z_Out_fine)

     !perform kd-tree on 1D data; data in 2D format needs to be reformatted to 1D 
     !create tree from predictor_1d, apply kdtree on query_1d to get indices and then apply found 
     !indices to Z_In_coarse_1d to compute Z_Out_1d

     !pred_var_coarse is known variable at coarse resolution, pred_var_fine at
     !fine resolution, Num_Elements and Num_lines are number of elements and
     !lines at coarse and fine resolution, respectively; nnbrute is the wanted number
     !of found closet indices for each search, Z_In_coarse is the known wanted
     !variable at coarse resolution, Z_Out_fine at fine resolution

     real (kind=real4), intent(in), dimension(:,:), target:: pred_var_coarse
     real (kind=real4), intent(in), dimension(:,:), target:: pred_var_fine
     integer (kind=int4), intent(in):: Num_Elements_coarse
     integer (kind=int4), intent(in):: Num_lines_coarse
     integer (kind=int4), intent(in):: Num_Elements_fine
     integer (kind=int4), intent(in):: Num_lines_fine
     integer (kind=int4), intent(in):: nnbrute
     real (kind=real4), intent(in), dimension(:,:), target:: Z_In_coarse
     real (kind=real4), intent(out), dimension(:,:):: Z_Out_fine

     ! define local variables, each KD-tree predictor is 1-D array
     type(kdtree2), pointer :: tree
     integer(kind=int4):: n_training, n_query
     integer:: i,j
     integer(kind=int4), dimension(:), allocatable:: ind_training, ind_query
     real(kind=real4), dimension(:), allocatable:: predictor_1d
     real(kind=real4), dimension(:), allocatable:: query_1d
     real(kind=real4), dimension(:), allocatable:: Z_In_coarse_1d
     real(kind=real4), dimension(:), allocatable:: Z_Out_1d
     real(kind=real4), dimension(:,:), allocatable:: my_array
     real(kind=real4), dimension(:), allocatable:: query_vec
     type(kdtree2_result), dimension(:), allocatable :: results1

     integer:: Line_Idx
     integer:: Elem_Idx


     if (.not. allocated(predictor_1d)) allocate(predictor_1d(Num_Elements_coarse*Num_Lines_coarse))
     if (.not. allocated(query_1d)) allocate(query_1d(Num_Elements_fine*Num_Lines_fine))
     if (.not. allocated(Z_In_coarse_1d)) allocate(Z_In_coarse_1d(Num_Elements_coarse*Num_Lines_coarse))
     if (.not. allocated(Z_Out_1d)) allocate(Z_Out_1d(Num_Elements_fine*Num_Lines_fine))

     Z_Out_1d = MISSING_VALUE_REAL4

        n_training = Num_Elements_coarse*Num_Lines_coarse
        n_query = Num_Elements_fine*Num_Lines_fine

        predictor_1d = reshape(pred_var_coarse,(/Num_Elements_coarse*Num_Lines_coarse/))
        query_1d = reshape(pred_var_fine,(/Num_Elements_fine*Num_Lines_fine/))

        if (.not. allocated(ind_training)) allocate(ind_training(n_training))

        ! search for training and query indices and stored
        i = 1
        do Line_Idx = 1, Num_Lines_coarse
        do Elem_Idx = 1, Num_Elements_coarse
              ind_training(i) = Elem_Idx + Num_Elements_coarse*(Line_Idx-1)
              i = i+1
          enddo
        enddo


        if (.not. allocated(my_array)) allocate(my_array(1,Num_Elements_coarse*Num_Lines_coarse))
        if (.not. allocated(query_vec)) allocate(query_vec(1))

        my_array(1,:) = predictor_1d

        tree => kdtree2_create(my_array,sort=.true.,rearrange=.true.)  ! this is how you create a tree.

        if (.not. allocated(results1)) allocate(results1(nnbrute))

        do i = 1,n_query
           query_vec(1) = query_1d(i)

           call kdtree2_n_nearest(tp=tree, qv=query_vec, nn = nnbrute, results=results1)

           ! average values for the all found indices
           Z_Out_1d(i) = sum(Z_In_coarse_1d(ind_training(results1%idx)))/size(results1%idx)
        enddo

        call kdtree2_destroy(tree)

     ! Change 1d array back to 2d; if no kdtree search is performed, array is
     ! empty
     Z_Out_fine = reshape(Z_Out_1d,(/Num_Elements_fine,Num_Lines_fine/))

     if (allocated(predictor_1d)) deallocate(predictor_1d)
     if (allocated(query_1d)) deallocate(query_1d)
     if (allocated(Z_In_coarse_1d)) deallocate(Z_In_coarse_1d)
     if (allocated(Z_Out_1d)) deallocate(Z_Out_1d)
     if (allocated(query_vec)) deallocate(query_vec)
     if (allocated(my_array)) deallocate(my_array)
     if (allocated(ind_training)) deallocate(ind_training)
     if (allocated(results1)) deallocate(results1)

  end subroutine KDTREE_COARSE2FINE

end module CX_NUCAPS_MOD
