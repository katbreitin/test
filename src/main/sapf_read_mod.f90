!$Id: sapf_read_mod.f90 2573 2018-03-22 12:45:54Z heidinger $
!--------------------------------------------------------------------------------------
! This module reads the SAPF (Enterprise) cloud mask and cloud type
!
! JRR-CloudMask_v1r2_j01_s201803310945051_e201803310946296_c201803311016360.nc
! JRR-CloudPhase_v1r1_npp_s201803312340054_e201803312341296_c201804010056180.nc
!
! Author: Denis B. (denis.botambekov@ssec.wisc.edu)
!
! Date: 2018-04-04
!
!--------------------------------------------------------------------------------------

module SAPF_READ_MOD

    use PIXEL_COMMON_MOD, only: &
             Image &
           , Sensor &
           , Cldmask &
           , Cld_Type_Aux &
           , Cld_Phase_Aux &
           , Zc_Aux &
           , Pc_Top1_Aux &
           , Use_Aux_Flag &
           , Cloud_Mask_Aux_Read_Flag &
           , Cloud_Prob_Aux_Read_Flag &
           , Cloud_Type_Aux_Read_Flag &
           , Cloud_Height_Aux_Read_Flag
    use CX_DATE_TIME_TOOLS_MOD, only: JULIAN, COMPUTE_MONTH, COMPUTE_DAY, LEAP_YEAR_FCT
    use FILE_UTILS, only: FILE_SEARCH
    use CONSTANTS_MOD, only: &
             Real4 &
           , Int4 &
           , Int2 &
           , Int1 &
           , Exe_Prompt &
           , Missing_Value_Real4 &
           , Sym &
           , Missing_Value_Int1
    use CX_NETCDF4_MOD,only: &
      close_netcdf &
      , open_netcdf &
      , read_netcdf &
      , read_netcdf_dimension_2d
    use CLAVRX_MESSAGE_MOD, only: &
        MESG, VERB_LEV

    implicit none
    private
    public:: DETERMINE_SAPF_NAME
    public:: READ_SAPF_DATA

    character(len=13), parameter:: SAPF_PROMPT="SAPF_MODULE:"

    contains

!------------------------------------------------------------------------------------
! Determine the name of the SAPF file for the Level1b file
!------------------------------------------------------------------------------------
 subroutine DETERMINE_SAPF_NAME(Seg_Idx)

  integer, intent(in):: Seg_Idx
  character(len=100):: Search_String
  character(len=1020), dimension(:), pointer:: Files
  integer:: Num_Files
  character(len=4):: Year_String, Time_String
  character(len=2):: Month_String, Dom_String
  character(len=3):: Doy_String
  

  Image%Auxiliary_Cloud_Mask_File_Name = 'no_file'

  !--- File name
  !0        1         2         3         4         5         6         7
  !12345678901234567890123456789012345678901234567890123456789012345678901234567890
  !GMTCO_j01_d20180331_t0945051_e0946296_b01888_c20180331100603384280_nobc_ops.h5
  !JRR-CloudMask_v1r2_j01_s201803310945051_e201803310946296_c201803311016360.nc
  !JRR-CloudPhase_v1r1_npp_s201803312359593_e201804010001235_c201804010103370.nc

  select case (trim(Sensor%Platform_Name))
     case('SNPP','NOAA-20') 
        Year_String = Image%Level1b_Name(12:15)
        Month_String = Image%Level1b_Name(16:17)
        Dom_String = Image%Level1b_Name(18:19)
        Time_String = Image%Level1b_Name(22:25)
     case('GOES-16')
        Year_String = Image%Level1b_Name(28:31)
        Doy_String = Image%Level1b_Name(32:34)
        Time_String = Image%Level1b_Name(35:39)
     case default
        call MESG (SAPF_PROMPT//"Satellite Not Supported",level = verb_lev % WARNING)
        return 
  end select

! search for mask file
  select case (trim(Sensor%Platform_Name))
     case('SNPP')
        Search_String = 'JRR-CloudMask_*_npp_s'//Year_String//Month_String//Dom_String//Time_String//'*nc'
     case('NOAA-20')
        Search_String = 'JRR-CloudMask_*_j01_s'//Year_String//Month_String//Dom_String//Time_String//'*nc'
     case('GOES-16')
!        Search_String = 'GOESR_ABI_FD_'//Year_String//Doy_String//'_'//Time_String//'*MASK_BL.nc'
        Search_String = 'GOES16_ABI_2KM_FD_'//Year_String//Doy_String//'_'//Time_String//'*MASK_EN.nc'
     case default
        call MESG (SAPF_PROMPT//"Satellite Not Supported",level = verb_lev % WARNING)
        return
  end select

  Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0) then
        if (Seg_Idx == 1) &
          call MESG (SAPF_PROMPT//"SAPF Cloud Mask File Not Found ",level = verb_lev % WARNING)
    else
        Image%Auxiliary_Cloud_Mask_File_Name = Files(1)
        if (Seg_Idx == 1) &
           call MESG (SAPF_PROMPT//"SAPF Cloud Mask File Found, "//trim(Image%Auxiliary_Cloud_Mask_File_Name))
    endif

  Files => null()

!stop
! search for type file
  Image%Auxiliary_Cloud_Type_File_Name = 'no_file'
  select case (trim(Sensor%Platform_Name))
     case('SNPP')
        Search_String = 'JRR-CloudPhase_*_npp_s'//Year_String//Month_String//Dom_String//Time_String//'*nc'
     case('NOAA-20')
        Search_String = 'JRR-CloudPhase_*_j01_s'//Year_String//Month_String//Dom_String//Time_String//'*nc'
     case('GOES-16')
!        Search_String = 'GOESR_ABI_FD_'//Year_String//Doy_String//'_'//Time_String//'*PHASE_BL.nc'
        Search_String = 'GOES16_ABI_2KM_FD_'//Year_String//Doy_String//'_'//Time_String//'*PHASE_EN.nc'

     case default
        call MESG (SAPF_PROMPT//"Satellite Not Supported",level = verb_lev % WARNING)
        return
  end select

  Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0) then
        if (Seg_Idx == 1) &
          call MESG (SAPF_PROMPT//"SAPF Cloud Phase File Not Found ",level = verb_lev % WARNING)
    else
        Image%Auxiliary_Cloud_Type_File_Name = Files(1)
        if (Seg_Idx == 1) &
            call MESG (SAPF_PROMPT//"SAPF Cloud Phase File Found, "//trim(Image%Auxiliary_Cloud_Type_File_Name))
    endif

  Files => null()


! search for height file
  Image%Auxiliary_Cloud_Height_File_Name = 'no_file'
  select case (trim(Sensor%Platform_Name))
     case('SNPP')
        Search_String = 'JRR-CloudHeight_*_npp_s'//Year_String//Month_String//Dom_String//Time_String//'*nc'
     case('NOAA-20')
        Search_String = 'JRR-CloudHeight_*_j01_s'//Year_String//Month_String//Dom_String//Time_String//'*nc'
     case('GOES-16')
!        Search_String = 'GOESR_ABI_FD_'//Year_String//Doy_String//'_'//Time_String//'*HEIGHT_BL.nc'
        Search_String = 'GOES16_ABI_2KM_FD_'//Year_String//Doy_String//'_'//Time_String//'*HEIGHT_EN.nc'
     case default
        call MESG (SAPF_PROMPT//"Satellite Not Supported",level = verb_lev % WARNING)
        return
  end select

  Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0) then
        if (Seg_Idx == 1) &
          call MESG (SAPF_PROMPT//"SAPF Cloud Height File Not Found ",level = verb_lev % WARNING)
    else
        Image%Auxiliary_Cloud_Height_File_Name = Files(1)
        if (Seg_Idx == 1) &
            call MESG (SAPF_PROMPT//"SAPF Cloud Height File Found, "//trim(Image%Auxiliary_Cloud_Height_File_Name))
    endif

  Files => null()


 !print *, "SAPF MASK File = ", trim(Image%Auxiliary_Cloud_Mask_File_Name)
 !print *, "SAPF TYPE File = ", trim(Image%Auxiliary_Cloud_Type_File_Name)
 !print *, "SAPF HEIGHT File = ", trim(Image%Auxiliary_Cloud_Height_File_Name)

 end subroutine DETERMINE_SAPF_NAME

!------------------------------------------------------------------------------------------
! Read SAPF cloud mask and cloud type
!------------------------------------------------------------------------------------------
 subroutine READ_SAPF_DATA(Seg_Idx)

    integer, intent(in):: Seg_Idx
    integer:: Sd_Id, Status_Flag
    integer(kind=int4), dimension(2):: Sds_Dims, Sds_Stride, Sds_Start, Sds_Edges
    character(len=120):: Sds_Name
    integer(kind=int1), dimension(:,:), allocatable:: I1_Buffer
    real(kind=real4), dimension(:,:), allocatable:: R4_Buffer
    integer:: Nx_Slab_Read_Start, Nx_Seg, Nx_Slab_Count
    integer:: Ny_Slab_Read_Start, Ny_Seg, Ny_Seg_Max, Ny_Slab_Count


    !--- 
    Status_Flag = 0
    Cloud_Mask_Aux_Read_Flag = sym%NO   !will be set to yes if successful
    Cloud_Type_Aux_Read_Flag = sym%NO 
    Cloud_Height_Aux_Read_Flag = sym%NO

    ! --- determine expected size of data slab - same for all
    Nx_Seg = Image%Number_of_Elements
    Ny_Seg_Max = Image%Number_of_Lines_Per_Segment
    Ny_Seg = Image%Number_of_Lines_Read_This_Segment
   

    !--- READ SAPF CLOUD MASK
    Sds_Name = "CloudMask"
    if (trim(Image%Auxiliary_Cloud_Mask_File_Name) /= "no_file") then

      !--- open file for read
      call OPEN_NETCDF(trim(Image%Level1b_Path)//trim(Image%Auxiliary_Cloud_Mask_File_Name), Sd_Id)

      !--- if file is unreadable, exit
      if (Sd_Id <= 0) then
           Status_Flag = 1
           call MESG (SAPF_PROMPT//"SAPF Cloud Mask File Could Not Be Opened ",level = verb_lev % WARNING)
           return
       endif

      ! --- get dimension of file
      call READ_NETCDF_DIMENSION_2D(Sd_Id, trim(Sds_Name), Sds_Dims)

      ! --- define size of data
      Nx_Slab_Read_Start = 1
      Ny_Slab_Read_Start = (Seg_Idx-1) * Ny_Seg_Max + 1
      Nx_Slab_Count = Sds_Dims(1)
      Ny_Slab_Count = Ny_Seg

      ! --- constrain to size of data
      if ((Seg_Idx * Ny_Seg) .gt. Sds_Dims(2)) Ny_Slab_Count = Sds_Dims(2) - Ny_Slab_Read_Start

      Sds_Stride = (/1, 1/)
      Sds_Start = (/Nx_Slab_Read_Start, Ny_Slab_Read_Start/)
      Sds_Edges = (/Nx_Slab_Count, Ny_Slab_Count/)

      ! --- allocate buffer
      if (.not. allocated(I1_Buffer)) allocate(I1_Buffer(Nx_Slab_Count, Ny_Slab_Count),stat=Status_Flag)

      ! --- read data
      call READ_NETCDF(Sd_Id, Sds_Start, Sds_Stride, Sds_Edges, trim(Sds_Name), I1_Buffer)

      !--- close file
      call CLOSE_NETCDF(Sd_Id)

      !--- save cloud mask
      Cldmask%Cld_Mask_Aux(:,1:Ny_Slab_Count) = I1_Buffer(:,1:Ny_Slab_Count)
 
      ! --- set missing to CLAVR-x missing
      where (Cldmask%Cld_Mask_Aux .lt. 0 .or. Cldmask%Cld_Mask_Aux .gt. 3)
         Cldmask%Cld_Mask_Aux = Missing_Value_Int1
      endwhere

      if (allocated(I1_Buffer)) deallocate(I1_Buffer,stat=Status_Flag)

      Cloud_Mask_Aux_Read_Flag = sym%YES

    endif

    !--- READ SAPF CLOUD Probability
    Status_Flag = 0
    Sds_Name = "CloudProbability"

    if (trim(Image%Auxiliary_Cloud_Mask_File_Name) /= "no_file") then

      !--- open file for read
      call OPEN_NETCDF(trim(Image%Level1b_Path)//trim(Image%Auxiliary_Cloud_Mask_File_Name),Sd_Id)

      !--- if file is unreadable, exit
      if (Sd_Id <= 0) then
          Status_Flag = 1
          call MESG (SAPF_PROMPT//"SAPF Cloud Probability File Could Not Be Opened ",level = verb_lev % WARNING)
          return
      endif

      ! --- get dimension of file
      call READ_NETCDF_DIMENSION_2D(Sd_Id, trim(Sds_Name), Sds_Dims)

      ! --- define size of data
      Nx_Slab_Read_Start = 1
      Ny_Slab_Read_Start = (Seg_Idx-1) * Ny_Seg_Max + 1
      Nx_Slab_Count = Sds_Dims(1)
      Ny_Slab_Count = Ny_Seg

      ! --- constrain to size of data
      if ((Seg_Idx * Ny_Seg) .gt. Sds_Dims(2)) Ny_Slab_Count = Sds_Dims(2) - Ny_Slab_Read_Start

      Sds_Stride = (/1, 1/)
      Sds_Start = (/Nx_Slab_Read_Start, Ny_Slab_Read_Start/)
      Sds_Edges = (/Nx_Slab_Count, Ny_Slab_Count/)

      ! --- allocate buffer
      if (.not. allocated(R4_Buffer)) allocate(R4_Buffer(Nx_Slab_Count, Ny_Slab_Count),stat=Status_Flag)

      ! --- read height data
      call READ_NETCDF(Sd_Id, Sds_Start, Sds_Stride, Sds_Edges, trim(Sds_Name), R4_Buffer)

      !--- save cloud height
      Cldmask%Posterior_Cld_Probability_Aux(:,1:Ny_Slab_Count) = R4_Buffer(:,1:Ny_Slab_Count)

      ! --- set missing to CLAVR-x missing
      where (Cldmask%Posterior_Cld_Probability_Aux .lt. 0)
         Cldmask%Posterior_Cld_Probability_Aux = Missing_Value_Real4
      endwhere

      !--- close file
      call CLOSE_NETCDF(Sd_Id)

      Cloud_Prob_Aux_Read_Flag = sym%YES

    endif




    !--- READ SAPF CLOUD TYPE
    Status_Flag = 0
    Sds_Name = "CloudType"

    if (trim(Image%Auxiliary_Cloud_Type_File_Name) /= "no_file") then

      !--- open file for read
      call OPEN_NETCDF(trim(Image%Level1b_Path)//trim(Image%Auxiliary_Cloud_Type_File_Name), Sd_Id)

      !--- if file is unreadable, exit
      if (Sd_Id <= 0) then
          Status_Flag = 1
          call MESG (SAPF_PROMPT//"SAPF Cloud Type/Phase File Could Not Be Opened ",level = verb_lev % WARNING)
          return
      endif

      ! --- get dimension of file
      call READ_NETCDF_DIMENSION_2D(Sd_Id, trim(Sds_Name), Sds_Dims)

      ! --- define size of data
      Nx_Slab_Read_Start = 1
      Ny_Slab_Read_Start = (Seg_Idx-1) * Ny_Seg_Max + 1
      Nx_Slab_Count = Sds_Dims(1)
      Ny_Slab_Count = Ny_Seg

      ! --- constrain to size of data
      if ((Seg_Idx * Ny_Seg) .gt. Sds_Dims(2)) Ny_Slab_Count = Sds_Dims(2) - Ny_Slab_Read_Start

      Sds_Stride = (/1, 1/)
      Sds_Start = (/Nx_Slab_Read_Start, Ny_Slab_Read_Start/)
      Sds_Edges = (/Nx_Slab_Count, Ny_Slab_Count/)

      ! --- allocate buffer
      if (.not. allocated(I1_Buffer)) allocate(I1_Buffer(Nx_Slab_Count, Ny_Slab_Count),stat=Status_Flag)

      ! --- read type data
      call READ_NETCDF(Sd_Id, Sds_Start, Sds_Stride, Sds_Edges, trim(Sds_Name), I1_Buffer)

      !--- save cloud type
      Cld_Type_Aux(:,1:Ny_Slab_Count) = I1_Buffer(:,1:Ny_Slab_Count)

      where(Cld_Type_Aux(:,1:Ny_Slab_Count) .eq. 1 .or. Cld_Type_Aux(:,1:Ny_Slab_Count) .eq. 8)
         Cld_Type_Aux(:,1:Ny_Slab_Count) = 10 !Cld_Type_Aux(:,1:Ny_Slab_Count) + 2
      endwhere

      where(Cld_Type_Aux(:,1:Ny_Slab_Count) .ge. 2 .and. Cld_Type_Aux(:,1:Ny_Slab_Count) .le. 7)
         Cld_Type_Aux(:,1:Ny_Slab_Count) = Cld_Type_Aux(:,1:Ny_Slab_Count) + 1
      endwhere

      ! --- set missing to CLAVR-x missing
      where (Cld_Type_Aux .lt. 0 .or. Cld_Type_Aux .gt. 8)
         Cld_Type_Aux = Missing_Value_Int1
      endwhere

      ! --- reset buffer to missing
      I1_Buffer = Missing_Value_Int1

      ! --- read phase data
      Sds_Name = "CloudPhase"
      call READ_NETCDF(Sd_Id, Sds_Start, Sds_Stride, Sds_Edges, trim(Sds_Name), I1_Buffer)

      !--- save cloud phase
      Cld_Phase_Aux(:,1:Ny_Slab_Count) = I1_Buffer(:,1:Ny_Slab_Count)

      ! --- set missing to CLAVR-x missing
      where (Cld_Phase_Aux .lt. 0 .or. Cld_Phase_Aux .gt. 5)
         Cld_Phase_Aux = Missing_Value_Int1
      endwhere
  
      if (allocated(I1_Buffer)) deallocate(I1_Buffer,stat=Status_Flag)
 
      !--- close file
      call CLOSE_NETCDF(Sd_Id)

      Cloud_Type_Aux_Read_Flag = sym%YES
    endif


    !--- READ SAPF CLOUD Height
    Status_Flag = 0
    Sds_Name = "CldTopHght"

    if (trim(Image%Auxiliary_Cloud_Height_File_Name) /= "no_file") then

      !--- open file for read
      call OPEN_NETCDF(trim(Image%Level1b_Path)//trim(Image%Auxiliary_Cloud_Height_File_Name),Sd_Id)

      !--- if file is unreadable, exit
      if (Sd_Id <= 0) then
          Status_Flag = 1
          call MESG (SAPF_PROMPT//"SAPF Cloud Height File Could Not Be Opened ",level = verb_lev % WARNING)
          return
      endif

      ! --- get dimension of file
      call READ_NETCDF_DIMENSION_2D(Sd_Id, trim(Sds_Name), Sds_Dims)

      ! --- define size of data
      Nx_Slab_Read_Start = 1
      Ny_Slab_Read_Start = (Seg_Idx-1) * Ny_Seg_Max + 1
      Nx_Slab_Count = Sds_Dims(1)
      Ny_Slab_Count = Ny_Seg

      ! --- constrain to size of data
      if ((Seg_Idx * Ny_Seg) .gt. Sds_Dims(2)) Ny_Slab_Count = Sds_Dims(2) - Ny_Slab_Read_Start

      Sds_Stride = (/1, 1/)
      Sds_Start = (/Nx_Slab_Read_Start, Ny_Slab_Read_Start/)
      Sds_Edges = (/Nx_Slab_Count, Ny_Slab_Count/)

      ! --- allocate buffer
      if (.not. allocated(R4_Buffer)) allocate(R4_Buffer(Nx_Slab_Count, Ny_Slab_Count),stat=Status_Flag)

      ! --- read height data
      call READ_NETCDF(Sd_Id, Sds_Start, Sds_Stride, Sds_Edges, trim(Sds_Name), R4_Buffer)

      !--- save cloud height
      Zc_Aux(:,1:Ny_Slab_Count) = R4_Buffer(:,1:Ny_Slab_Count)

      ! --- set missing to CLAVR-x missing
      where (Zc_Aux .lt. 0)
         Zc_Aux = Missing_Value_Real4
      endwhere

      R4_Buffer = Missing_Value_Real4

      ! --- read pressure data
      Sds_Name = "CldTopPres"
      call READ_NETCDF(Sd_Id, Sds_Start, Sds_Stride, Sds_Edges, trim(Sds_Name), R4_Buffer)

      !--- save cloud phase
      Pc_Top1_Aux(:,1:Ny_Slab_Count) = R4_Buffer(:,1:Ny_Slab_Count)

      ! --- set missing to CLAVR-x missing
      where (Pc_Top1_Aux .lt. 0)
         Pc_Top1_Aux = Missing_Value_Real4
      endwhere

      if (allocated(R4_Buffer)) deallocate(R4_Buffer,stat=Status_Flag)

      !--- close file
      call CLOSE_NETCDF(Sd_Id)

      Cloud_Height_Aux_Read_Flag = sym%YES
    endif

      if (allocated(I1_Buffer)) deallocate(I1_Buffer)
      if (allocated(R4_Buffer)) deallocate(R4_Buffer)

 end subroutine READ_SAPF_DATA

!-------------------------------------------
!
!-------------------------------------------

end module SAPF_READ_MOD


