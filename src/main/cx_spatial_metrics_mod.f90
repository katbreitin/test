! $Id: pixel_routines_mod.f90 3401 2019-07-10 20:01:39Z heidinger $
!-------------------------------------------------------------------------------------- 
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: pixel_routines.f90 (src)
!       PIXEL_ROUTINES (program)
!
! PURPOSE: this module houses routines for computing some needed pixel-level arrays
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
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
! Public routines used in this MODULE:
!
! This module houses all of the routines that make spatial metrics from
! satellite data
!
!--------------------------------------------------------------------------------------
module CX_SPATIAL_METRICS_MOD
 use CONSTANTS_MOD
 use PIXEL_COMMON_MOD
 use NUMERICAL_ROUTINES_MOD, only: Covariance
 use CX_REAL_BOOLEAN_MOD
 use ACHA_CLAVRX_BRIDGE,only: &
     AWG_CLOUD_HEIGHT_BRIDGE, &
     LOCAL_LINEAR_RADIATIVE_CENTER

 implicit none
 private

 public:: COMPUTE_SPATIAL_CORRELATION_ARRAYS, &
          COMPUTE_MIN_MAX_MEAN_STD_METRICS, &
          COMPUTE_MEDIAN_METRICS_L1B, &
          COMPUTE_MEDIAN_METRICS_L2, &
          COMPUTE_RADIATIVE_CENTER_ARRAYS, &
          COMPUTE_MEDIAN

 private:: COMPUTE_NxN_METRICS, &
           COMPUTE_NxN_MIN_MAX_INDICES, &
           COMPUTE_MEDIAN_SEGMENT, &
           GRADIENT_MEANDER

  contains

  !------------------------------------------------------------------------------------------------
  ! a collection of calls to COMPUTE_NxN_METRICS 
  !------------------------------------------------------------------------------------------------
  subroutine COMPUTE_MIN_MAX_MEAN_STD_METRICS()
     integer:: dim1
     integer:: dim2

     !--- allocate pixel arrays
     dim1 = Image%Number_Of_Elements
     dim2 = Image%Number_Of_Lines_Per_Segment

            if (Sensor%Chan_On_Flag_Default(1)==sym%Yes) then
              if (.not. allocated(Ch(1)%Ref_Toa_Min_3x3)) allocate(Ch(1)%Ref_Toa_Min_3x3(dim1,dim2))
              if (.not. allocated(Ch(1)%Ref_Toa_Max_3x3)) allocate(Ch(1)%Ref_Toa_Max_3x3(dim1,dim2))
              if (.not. allocated(Ch(1)%Ref_Toa_Mean_3x3)) allocate(Ch(1)%Ref_Toa_Mean_3x3(dim1,dim2))
              if (.not. allocated(Ch(1)%Ref_Toa_Std_3x3)) allocate(Ch(1)%Ref_Toa_Std_3x3(dim1,dim2))
              call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Ch(1)%Ref_Toa, &
                            Ch(1)%Ref_Toa_Min_3x3,Ch(1)%Ref_Toa_Max_3x3, &
                            Ch(1)%Ref_Toa_Mean_3x3,Ch(1)%Ref_Toa_Std_3x3)
            endif
            if (Sensor%Chan_On_Flag_Default(1)==sym%Yes) then
              if (.not. allocated(Ch(1)%Ref_Toa_Clear_Std_3x3)) allocate(Ch(1)%Ref_Toa_Clear_Std_3x3(dim1,dim2))
              call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Ch(1)%Ref_Toa_Clear, &
                            Temp_Pix_Array_1,Temp_Pix_Array_1,Temp_Pix_Array_1,Ch(1)%Ref_Toa_Clear_Std_3x3)
            endif
            if (Sensor%Chan_On_Flag_Default(1)==sym%Yes) then
              if (.not. allocated(Ch(1)%Sfc_Ref_White_Sky_Mean_3x3)) allocate(Ch(1)%Sfc_Ref_White_Sky_Mean_3x3(dim1,dim2))
              call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Ch(1)%Sfc_Ref_White_Sky, &
                            Temp_Pix_Array_1,Temp_Pix_Array_1,Ch(1)%Sfc_Ref_White_Sky_Mean_3x3,Temp_Pix_Array_1)
            endif

            if (Sensor%Chan_On_Flag_Default(20)==sym%Yes) then
              if (.not. allocated(Ch(20)%Bt_Toa_Std_3x3)) allocate(Ch(20)%Bt_Toa_Std_3x3(dim1,dim2))
              call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Ch(20)%Bt_Toa, &
                            Temp_Pix_Array_1,Temp_Pix_Array_1,Temp_Pix_Array_1,Ch(20)%Bt_Toa_Std_3x3)
!print *, "here"
              if (.not. allocated(Ch(20)%Ref_Toa_Mean_3x3)) allocate(Ch(20)%Ref_Toa_Mean_3x3(dim1,dim2))
              if (.not. allocated(Ch(20)%Ref_Toa_Min_3x3)) allocate(Ch(20)%Ref_Toa_Min_3x3(dim1,dim2))
              if (.not. allocated(Ch(20)%Ref_Toa_Max_3x3)) allocate(Ch(20)%Ref_Toa_Max_3x3(dim1,dim2))
              if (.not. allocated(Ch(20)%Ref_Toa_Std_3x3)) allocate(Ch(20)%Ref_Toa_Std_3x3(dim1,dim2))
!print *, "there"
              call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Ch(20)%Ref_Toa, &
                            Ch(20)%Ref_Toa_Min_3x3,Ch(20)%Ref_Toa_Max_3x3, &
                            Ch(20)%Ref_Toa_Mean_3x3,Ch(20)%Ref_Toa_Std_3x3)
!print *, "bear"
            endif

            if (Sensor%Chan_On_Flag_Default(27)==sym%Yes) then
              if (.not. allocated(Ch(27)%Bt_Toa_Max_3x3)) allocate(Ch(27)%Bt_Toa_Max_3x3(dim1,dim2))
              if (.not. allocated(Ch(27)%Bt_Toa_Std_3x3)) allocate(Ch(27)%Bt_Toa_Std_3x3(dim1,dim2))
              call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Ch(27)%Bt_Toa, &
                            Temp_Pix_Array_1,ch(27)%Bt_Toa_Max_3x3,Temp_Pix_Array_1,Ch(27)%Bt_Toa_Std_3x3)
            endif

            if (Sensor%Chan_On_Flag_Default(31)==sym%Yes) then
              if (.not. allocated(Ch(31)%Bt_Toa_Min_3x3)) allocate(Ch(31)%Bt_Toa_Min_3x3(dim1,dim2))
              if (.not. allocated(Ch(31)%Bt_Toa_Max_3x3)) allocate(Ch(31)%Bt_Toa_Max_3x3(dim1,dim2))
              if (.not. allocated(Ch(31)%Bt_Toa_Mean_3x3)) allocate(Ch(31)%Bt_Toa_Mean_3x3(dim1,dim2))
              if (.not. allocated(Ch(31)%Bt_Toa_Std_3x3)) allocate(Ch(31)%Bt_Toa_Std_3x3(dim1,dim2))

              call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Ch(31)%Bt_Toa, &
                            Ch(31)%Bt_Toa_Min_3x3,Ch(31)%Bt_Toa_Max_3x3, &
                            Ch(31)%Bt_Toa_Mean_3x3,Ch(31)%Bt_Toa_Std_3x3)
            endif

            if (Sensor%Chan_On_Flag_Default(38)==sym%Yes) then
              if (.not. allocated(Ch(38)%Bt_Toa_Min_3x3)) allocate(Ch(38)%Bt_Toa_Min_3x3(dim1,dim2))
              if (.not. allocated(Ch(38)%Bt_Toa_Max_3x3)) allocate(Ch(38)%Bt_Toa_Max_3x3(dim1,dim2))
              if (.not. allocated(Ch(38)%Bt_Toa_Mean_3x3)) allocate(Ch(38)%Bt_Toa_Mean_3x3(dim1,dim2))
              if (.not. allocated(Ch(38)%Bt_Toa_Std_3x3)) allocate(Ch(38)%Bt_Toa_Std_3x3(dim1,dim2))
              call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Ch(38)%Bt_Toa, &
                            Ch(38)%Bt_Toa_Min_3x3,Ch(38)%Bt_Toa_Max_3x3, &
                            Ch(38)%Bt_Toa_Mean_3x3,Ch(38)%Bt_Toa_Std_3x3)
            endif

            if (Sensor%Chan_On_Flag_Default(44)==sym%Yes) then
               if (.not. allocated(Ch(44)%Ref_Lunar_Min_3x3)) allocate(Ch(44)%Ref_Lunar_Min_3x3(dim1,dim2))
               if (.not. allocated(Ch(44)%Ref_Lunar_Max_3x3)) allocate(Ch(44)%Ref_Lunar_Max_3x3(dim1,dim2))
               if (.not. allocated(Ch(44)%Ref_Lunar_Mean_3x3)) allocate(Ch(44)%Ref_Lunar_Mean_3x3(dim1,dim2))
               if (.not. allocated(Ch(44)%Ref_Lunar_Std_3x3)) allocate(Ch(44)%Ref_Lunar_Std_3x3(dim1,dim2))
               call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Ch(44)%Ref_Lunar_Toa, &
                            Ch(44)%Ref_Lunar_Min_3x3,Ch(44)%Ref_Lunar_Max_3x3, &
                            Ch(44)%Ref_Lunar_Mean_3x3,Ch(44)%Ref_Lunar_Std_3x3)
            endif

            if (Sensor%Chan_On_Flag_Default(31)==sym%Yes) then
               CALL COMPUTE_NxN_MIN_MAX_INDICES( &
                                      ch(31)%Bt_Toa,1, int(sym%YES), &
                                      Bad_Pixel_Mask, Sfc%Land_Mask, &
                                      Image%Number_Of_Elements, &
                                      Image%Number_of_Lines_Read_This_Segment, &
                                      Elem_Idx_Max_Bt_Ch31_3x3, &
                                      Line_Idx_Max_Bt_Ch31_3x3, &
                                      Elem_Idx_Min_Bt_Ch31_3x3, &
                                      Elem_Idx_Max_Bt_Ch31_3x3)
            endif

            if (Sensor%Chan_On_Flag_Default(38)==sym%Yes) then
               CALL COMPUTE_NxN_MIN_MAX_INDICES( &
                                      ch(38)%Bt_Toa,1, int(sym%YES), &
                                      Bad_Pixel_Mask, Sfc%Land_Mask, &
                                      Image%Number_Of_Elements, &
                                      Image%Number_of_Lines_Read_This_Segment, &
                                      Elem_Idx_Max_Bt_Ch38_3x3, &
                                      Line_Idx_Max_Bt_Ch38_3x3, &
                                      Elem_Idx_Min_Bt_Ch38_3x3, &
                                      Elem_Idx_Max_Bt_Ch38_3x3)
            endif

            ! surface std
            call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Sfc%Zsfc, &
                            Temp_Pix_Array_1,Temp_Pix_Array_1, &
                            Temp_Pix_Array_1,Sfc%Zsfc_Std)


            if (Sensor%Chan_On_Flag_Default(1)==sym%Yes .and. allocated(Ch(1)%Sfc_Ref_White_Sky_Mean_3x3)) then
               !--- fill in wholes in white sky albedo
               where(ch(1)%Sfc_Ref_White_Sky == Missing_Value_Real4 .and. &
                     ch(1)%Sfc_Ref_White_Sky_Mean_3x3 /= Missing_Value_Real4)
                     ch(1)%Sfc_Ref_White_Sky = ch(1)%Sfc_Ref_White_Sky_Mean_3x3
               end where
            endif

            if (allocated(Sfc%Zsfc)) then
               call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,Sfc%Zsfc, &
                            Temp_Pix_Array_1,Sfc%Zsfc_Max, &
                            Temp_Pix_Array_1,Sfc%Zsfc_Std)

            endif

            if (allocated(ch(31)%Bt_Toa_Max_Sub) .and. allocated(Bt_Ch43_Max_Sub_3x3)) then
               call COMPUTE_NxN_METRICS(1,Bad_Pixel_Mask,ch(31)%Bt_Toa_Max_Sub, &
                            Temp_Pix_Array_1,Bt_Ch43_Max_Sub_3x3, &
                            Temp_Pix_Array_1,Temp_Pix_Array_1)
            endif
            
 
   end subroutine COMPUTE_MIN_MAX_MEAN_STD_METRICS

  !------------------------------------------------------------------------------------------------
  !  populate local radiative center arrays
  !------------------------------------------------------------------------------------------------
   subroutine COMPUTE_RADIATIVE_CENTER_ARRAYS(LRC_Flag,Chan_Idx)

      integer, parameter:: LRC_Meander_Flag = 1
      integer, parameter:: Max_LRC_Distance = 10
      real, parameter:: Min_LRC_Jump = 0.0   !0.5
      real, parameter:: Max_LRC_Jump = 100.0 !10.0

      integer, parameter:: Grad_Flag_LRC = -1
      real, parameter:: Min_Bt_11um_LRC = 220.0
      real, parameter:: Max_Bt_11um_LRC = 300.0

      integer, intent(in):: LRC_Flag
      integer, intent(in):: Chan_Idx

            if (LRC_Flag == sym%YES) then
               if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then

                  call LOCAL_LINEAR_RADIATIVE_CENTER(sym%YES,sym%NO, &
                                   LRC_Meander_Flag, &
                                   ch(Chan_Idx)%Bt_Toa(:,Line_Idx_Min_Segment:Image%Number_Of_Lines_Read_This_Segment-Line_Idx_Min_Segment+1), &
                                   1, &
                                   Image%Number_Of_Elements, &
                                   Line_Idx_Min_Segment,  &
                                   Image%Number_Of_Lines_Read_This_Segment, &
                                   Max_LRC_Distance,  &
                                   Min_LRC_Jump,  &
                                   Max_LRC_Jump,  &
                                   Grad_Flag_LRC,  &
                                   Missing_Value_Int4, &
                                   Bad_Pixel_Mask(:,Line_Idx_Min_Segment:Image%Number_Of_Lines_Read_This_Segment-Line_Idx_Min_Segment+1), &
                                   Min_Bt_11um_LRC,  &
                                   Max_Bt_11um_LRC, &
                                   I_LRC(:,Line_Idx_Min_Segment:Image%Number_Of_Lines_Read_This_Segment-Line_Idx_Min_Segment+1), &
                                   J_LRC(:,Line_Idx_Min_Segment:Image%Number_Of_Lines_Read_This_Segment-Line_Idx_Min_Segment+1))
               endif
            endif

   end subroutine COMPUTE_RADIATIVE_CENTER_ARRAYS

   subroutine COMPUTE_MEDIAN_METRICS_L1B()

      integer(kind=int4), parameter:: One = 1
      integer(kind=int4), parameter:: Two = 2

            if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then

               !--- apply median filter to Bt_Ch20
               call COMPUTE_MEDIAN_SEGMENT(ch(20)%Bt_Toa,Bad_Pixel_Mask, Two, &
                              1, Image%Number_Of_Elements,Line_Idx_Min_Segment, &
                              Line_Idx_Min_Segment+Image%Number_Of_Lines_Read_This_Segment-One, &
                              Bt_Ch20_Median_5x5,  &
                              Temp_Pix_Array_1)

            endif

   end subroutine COMPUTE_MEDIAN_METRICS_L1B

   subroutine COMPUTE_MEDIAN_METRICS_L2()

      integer(kind=int4), parameter:: One = 1

            if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
               !--- apply median filter to Ems_Ch20
               call COMPUTE_MEDIAN_SEGMENT(ch(20)%Emiss_Rel_11um,Bad_Pixel_Mask, One, &
                              1, Image%Number_Of_Elements,Line_Idx_Min_Segment, &
                              Line_Idx_Min_Segment+Image%Number_Of_Lines_Read_This_Segment-One, &
                              Ems_Ch20_Median_3x3,  &
                              Temp_Pix_Array_1)
            endif
 
   end subroutine COMPUTE_MEDIAN_METRICS_L2
!--------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------
subroutine COMPUTE_SPATIAL_CORRELATION_ARRAYS()

    integer:: Elem_Idx
    integer:: Elem_Idx_min
    integer:: Elem_Idx_max
    integer:: Elem_Idx_width
    integer:: Elem_Idx_segment_max
    integer:: Line_Idx
    integer:: Line_Idx_min
    integer:: Line_Idx_max
    integer:: Line_Idx_width
    integer:: Line_Idx_segment_max

    Elem_Idx_segment_max = Image%Number_Of_Elements
    Line_Idx_segment_max = Line_Idx_Min_Segment + Line_Idx_Max_Segment - 1

    do Elem_Idx = 1, Elem_Idx_segment_max
       do Line_Idx = 1, Line_Idx_segment_max

        !--- compute 5x5 arrays
        Elem_Idx_min = max(1,min(Elem_Idx - 2,Elem_Idx_segment_max))
        Elem_Idx_max = max(1,min(Elem_Idx + 2,Elem_Idx_segment_max))
        Line_Idx_min = max(1,min(Line_Idx - 2,Line_Idx_segment_max))
        Line_Idx_max = max(1,min(Line_Idx + 2,Line_Idx_segment_max))
        Line_Idx_width = Line_Idx_max - Line_Idx_min + 1
        Elem_Idx_width = Elem_Idx_max - Elem_Idx_min + 1

        if ((Sensor%Chan_On_Flag_Per_Line(27,Line_Idx) == sym%YES) .and. & 
            (Sensor%Chan_On_Flag_Per_Line(31,Line_Idx) == sym%YES)) then

            Covar_Ch27_Ch31_5x5(Elem_Idx,Line_Idx) = Covariance(&
               ch(31)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               ch(27)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               Elem_Idx_width, Line_Idx_width, &
               Bad_Pixel_Mask(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max))

           if (AVHRR_Fusion_Flag) then
              Covar_Ch27_Ch31_5x5(Elem_Idx,Line_Idx) = Missing_Value_Real4
           endif

        endif

        if ((Sensor%Chan_On_Flag_Per_Line(37,Line_Idx) == sym%YES) .and. & 
            (Sensor%Chan_On_Flag_Per_Line(31,Line_Idx) == sym%YES)) then

            Covar_Ch37_Ch31_5x5(Elem_Idx,Line_Idx) = Covariance(&
               ch(31)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               ch(37)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               Elem_Idx_width, Line_Idx_width, &
               Bad_Pixel_Mask(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max))

           if (AVHRR_Fusion_Flag) then
              Covar_Ch37_Ch31_5x5(Elem_Idx,Line_Idx) = Missing_Value_Real4
           endif

        endif

        if ((Sensor%Chan_On_Flag_Per_Line(27,Line_Idx) == sym%YES) .and. & 
            (Sensor%Chan_On_Flag_Per_Line(38,Line_Idx) == sym%YES)) then

            Covar_Ch27_Ch38_5x5(Elem_Idx,Line_Idx) = Covariance(&
               ch(27)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               ch(38)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               Elem_Idx_width, Line_Idx_width, &
               Bad_Pixel_Mask(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max))

           if (AVHRR_Fusion_Flag) then
              Covar_Ch27_Ch38_5x5(Elem_Idx,Line_Idx) = Missing_Value_Real4
           endif

        endif

        if ((Sensor%Chan_On_Flag_Per_Line(37,Line_Idx) == sym%YES) .and. & 
            (Sensor%Chan_On_Flag_Per_Line(38,Line_Idx) == sym%YES)) then

            Covar_Ch37_Ch38_5x5(Elem_Idx,Line_Idx) = Covariance(&
               ch(37)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               ch(38)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               Elem_Idx_width, Line_Idx_width, &
               Bad_Pixel_Mask(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max))

           if (AVHRR_Fusion_Flag) then
              Covar_Ch37_Ch38_5x5(Elem_Idx,Line_Idx) = Missing_Value_Real4
           endif

        endif

      enddo

    enddo
    
end subroutine COMPUTE_SPATIAL_CORRELATION_ARRAYS

!----------------------------------------------------------------------
! compute spatial metrics for 3x3 elements in a 2d array
!----------------------------------------------------------------------
subroutine COMPUTE_NxN_METRICS(N,Bad_Mask,Z,Z_Min,Z_Max,Z_Mean,Z_Std)
  integer, intent(in):: N
  integer(kind=int1), dimension(:,:), intent(in):: Bad_Mask
  real, dimension(:,:), intent(in):: Z
  real, dimension(:,:), intent(out):: Z_Min,Z_Max,Z_Mean,Z_Std
  integer:: Nx, Ny, N_Good
  real*8:: Sum_Temp, Sum_Temp2
  real :: Min_Temp, Max_Temp
  integer:: Count_Temp
  integer:: i, i1, i2, ii, j, j1, j2, jj   !local indices

  Nx = size(Z,1)
  Ny = size(Z,2)

!--- allocating in this routine gave erroneous results
! real, dimension(:,:), allocatable, intent(out):: Z_Min,Z_Max,Z_Mean,Z_Std
! if (.not. allocated(Z_Min))  allocate(Z_Min(Nx,Ny))
! if (.not. allocated(Z_Max))  allocate(Z_Max(Nx,Ny))
! if (.not. allocated(Z_Mean)) allocate(Z_Mean(Nx,Ny))
! if (.not. allocated(Z_Std))  allocate(Z_Std(Nx,Ny))

  Z_Min = Missing_Value_Real4
  Z_Max = Missing_Value_Real4
  Z_Mean = Missing_Value_Real4
  Z_Std = Missing_Value_Real4

  line_loop: do j = 1, Ny

  !--- set limits of NxN array in the j-direction
  j1 = max(1,j-N)   !top index of local array
  j2 = min(Ny,j+N)   !bottom index of local array

  element_loop: do i = 1, Nx

     !--- initial checks for this pixel
     if (Bad_Mask(i,j) == sym%YES .or. (Z(i,j) .eqr. Missing_Value_Real4)) then
         cycle
     endif

     !--- set limits of NxN array in the i-direction
     i1 = max(1,i-N)   !left index of local array
     i2 = min(Nx,i+N)   !right index of local array

     !--- initialize
     Sum_Temp = 0.0
     Sum_Temp2 = 0.0
     Count_Temp = 0
     Min_temp = huge(z_min(i,j))
     Max_Temp = -1.0*huge(z_min(i,j))
     N_Good = 0

     !--- go through each element in NxN array
     sub_line_loop: do jj = j1,j2
      sub_element_loop: do ii = i1,i2

        if (Bad_Mask(ii,jj) == sym%YES) then
          cycle
        endif

        if (Z(ii,jj) .eqr. Missing_Value_Real4) then
          cycle
        endif

        N_Good = N_Good + 1
        sum_Temp = Sum_Temp + Z(ii,jj)

        if (Z(ii,jj) .ltr. Min_Temp) then
           Min_Temp = z(ii,jj)
        endif

        if (z(ii,jj) .gtr. Max_Temp) then
           Max_Temp  = z(ii,jj)
        endif

       end do sub_element_loop
     end do sub_line_loop

     !--- if any good pixels found, compute mean and standard deviation
     if (N_Good > 0) then
       Z_Mean(i,j) = sum_temp / N_Good
       N_Good = 0
         !--- go through each element in NxN array
         do jj = j1,j2
          do ii = i1,i2

            if (Bad_Mask(ii,jj) == sym%YES) then
              cycle
            endif

            if (Z(ii,jj) .eqr. Missing_Value_Real4) then
              cycle
            endif

            N_Good = N_Good + 1
            Sum_Temp2 = Sum_Temp2 + (Z(ii,jj)-Z_Mean(i,j))**2

           end do
         end do
       Z_Std(i,j) = sqrt(max(0.0,(Sum_temp2/N_Good)))
       Z_Min(i,j) = Min_Temp
       Z_Max(i,j) = Max_Temp
     endif

     if (Max_Temp .eqr. Missing_Value_Real4) print *, "Missing Max = ", N_Good, Max_Temp

 end do element_loop

end do line_loop

end subroutine COMPUTE_NxN_METRICS

!----------------------------------------------------------------------
! subroutine COMPUTE_NxN_MIN_MAX_INDICES
!
! z - the input array
! n - the size of the box (1=3x3, 2=5x5, ...)
! bad_mask - mask array (only values with sym%YES will contribute)
! imin - starting x-index of array
! imax - ending x-index of array
!
! i_loc_of_max - i index of the maximum value in nxn array
! i_loc_of_min - i index of the minimun value in nxn array
!----------------------------------------------------------------------
subroutine COMPUTE_NxN_MIN_MAX_INDICES( &
                                          z, &
                                          n, &
                                          uni_land_mask_flag, &
                                          bad_mask, &
                                          land_mask,  &
                                          imax, &
                                          jmax, &
                                          i_loc_of_max,  &
                                          j_loc_of_max,  &
                                          i_loc_of_min,  &
                                          j_loc_of_min)

  real(kind=real4), dimension(:,:), intent(in):: z
  integer(kind=int1), dimension(:,:), intent(in):: bad_mask
  integer(kind=int1), dimension(:,:), intent(in):: land_mask
  integer, intent(in):: uni_land_mask_flag
  integer, intent(in):: n
  integer, intent(in):: imax
  integer, intent(in):: jmax
  integer, dimension(:,:), intent(out):: i_loc_of_min
  integer, dimension(:,:), intent(out):: j_loc_of_min
  integer, dimension(:,:), intent(out):: i_loc_of_max
  integer, dimension(:,:), intent(out):: j_loc_of_max
  integer:: i
  integer:: j
  integer:: i1
  integer:: i2
  integer:: j1
  integer:: j2
  integer:: ii
  integer:: jj
  integer:: N_Good
  real:: Max_Temp
  real:: Min_Temp

  !--- initialize to missing
  i_loc_of_max = missing_value_int1
  j_loc_of_max = missing_value_int1
  i_loc_of_min = missing_value_int1
  j_loc_of_min = missing_value_int1

  line_loop: do j = 1, jmax

  !--- set limits of NxN array in the j-direction
  j1 = max(1,j-n)   !top index of local array
  j2 = min(jmax,j+n)   !bottom index of local array

  element_loop: do i = 1, imax


     !--- initial checks for this pixel
     if (bad_mask(i,j) == sym%YES .or. z(i,j) == Missing_Value_Real4) then
         cycle
      endif


     !--- set limits of NxN array in the i-direction
     i1 = max(1,i-n)   !left index of local array
     i2 = min(imax,i+n)   !right index of local array

     !--- initialize
     N_Good = 0
     Min_Temp = 1.0*huge(Min_Temp)
     Max_Temp = -1.0*huge(Max_Temp)

     !--- go through each element in NxN array
     sub_line_loop: do jj = j1,j2
      sub_element_loop: do ii = i1,i2

        if (bad_mask(ii,jj) == sym%YES) then
          cycle
        endif

        if (z(ii,jj) == missing_value_real4) then
          cycle
        endif
       
        if ((uni_land_mask_flag == sym%YES) .and. &
             (land_mask(i,j) /= land_mask(ii,jj))) then    
          cycle
        endif

        N_Good = N_Good + 1
 
        if (z(ii,jj) < Min_Temp) then
           Min_Temp = z(ii,jj)
           i_loc_of_min(i,j) = ii
           j_loc_of_min(i,j) = jj
        endif

        if (z(ii,jj) > Max_Temp) then
           Max_Temp = z(ii,jj)
           i_loc_of_max(i,j) = ii
           j_loc_of_max(i,j) = jj
        endif

       end do sub_element_loop
     end do sub_line_loop

 end do element_loop

end do line_loop

end subroutine COMPUTE_NxN_MIN_MAX_INDICES

!----------------------------------------------------------------------
!  Local Linear Radiative Center
!----------------------------------------------------------------------
subroutine GRADIENT_MEANDER(Meander_Flag, &
                            Grid_Data, &
                            Element_Start, Number_Of_Elements, & 
                            Line_Start, Number_Of_Lines, & 
                            Max_Grad_Distance, &
                            Grad_Flag,  &
                            Missing_LRC_Value, &
                            Skip_LRC_Mask, &
                            Min_Grid_Data_Valid, Max_Grid_Data_Valid, &
                            ielem_LRC, iline_LRC)

  integer, intent(in):: Meander_Flag
  real (kind=real4), intent(in), dimension(:,:) :: Grid_Data
  integer (kind=int4), intent(in):: Element_Start
  integer (kind=int4), intent(in):: Number_of_Elements
  integer (kind=int4), intent(in):: Line_Start
  integer (kind=int4), intent(in):: Number_of_Lines
  integer (kind=int4), intent(in):: Max_Grad_Distance
  integer (kind=int4), intent(in):: Grad_Flag
  integer (kind=int4), intent(in):: Missing_LRC_Value
  integer (kind=int1), intent(in), dimension(:,:):: Skip_LRC_Mask
  real (kind=real4), intent(in):: Min_Grid_Data_Valid
  real (kind=real4), intent(in):: Max_Grid_Data_Valid
  integer (kind=int4), intent(out), dimension(:,:):: ielem_LRC
  integer (kind=int4), intent(out), dimension(:,:):: iline_LRC
  real, dimension(3,3):: Grad_Array
  integer, dimension(2):: Grad_Indices
  integer:: ielem
  integer:: iline
  integer:: ielem_Previous
  integer:: iline_Previous
  integer:: ielem_Next
  integer:: iline_Next
  real:: Grad_Temp
  integer:: Element_End
  integer:: Line_End
  integer:: ipoint
  integer:: ielem_dir
  integer:: iline_dir
 
  Element_End = Number_of_Elements + Element_Start - 1
  Line_End = Number_of_Lines + Line_Start - 1

  !--- initialize
  ielem_LRC = Missing_LRC_Value
  iline_LRC = Missing_LRC_Value

!----------------------------------------------------------------------
! loop through pixels in segment
!----------------------------------------------------------------------
Element_Loop:  do ielem = Element_Start+1, Element_End-1
Line_Loop:    do iline = Line_Start+1, Line_End-1

      !--- skip data due to mask
      if (Skip_LRC_Mask(ielem,iline) == sym%YES) cycle

      !-- check for out of bounds data
      if (Grad_Flag ==  1 .and. Grid_Data(ielem,iline) < Min_Grid_Data_Valid) cycle
      if (Grad_Flag ==  -1 .and. Grid_Data(ielem,iline) > Max_Grid_Data_Valid) cycle

      !-- check for data that already meets LRC criteria
      if ((Grad_Flag ==  1 .and. Grid_Data(ielem,iline) > Max_Grid_Data_Valid) .or. &
          (Grad_Flag ==  -1 .and. Grid_Data(ielem,iline) < Min_Grid_Data_Valid)) then
              ielem_LRC(ielem,iline) = ielem
              iline_LRC(ielem,iline) = iline
      endif

      !--- initialize previous variables
      ielem_Previous = ielem
      iline_Previous = iline

      !---- go long gradient and check for a reversal or saturation
      do ipoint = 1,Max_Grad_Distance

        !--- compute local gradient, find strongest gradient in 3x3 array and compute direction
        if (ipoint == 1 .or. Meander_Flag == sym%YES) then

         !--- construct 3x3 array for analysis
         Grad_Array =  &
           Grid_Data(ielem_Previous-1:ielem_Previous+1,iline_Previous-1:iline_Previous+1) -  &
           Grid_Data(ielem_Previous,iline_Previous)

         !--- look for bad data
         if (minval(Grad_Array) == Missing_Value_Real4) exit 

         !--- compute local gradients, find strongest gradient
         if (Grad_Flag == 1) then
          Grad_Indices = maxloc(Grad_Array)
         else
          Grad_Indices = minloc(Grad_Array)
         endif 

         !--- compute direction
         ielem_Dir = Grad_Indices(1)  - 2
         iline_Dir = Grad_Indices(2)  - 2

         !--- check for pixels that are located at  minima/maxima
         if (ielem_Dir == 0 .and. iline_Dir == 0) then
           ielem_LRC(ielem,iline) = ielem_Previous
           iline_LRC(ielem,iline) = iline_Previous
           exit
         endif

        endif

        !-- select next point on the path
        ielem_Next = ielem_Previous + ielem_Dir
        iline_Next = iline_Previous + iline_Dir

        !--- check for hitting segment boundaries
        if (ielem_Next == Element_Start .or. ielem_Next == Element_End .or. &
             iline_Next == Line_Start .or. iline_Next == Line_End) then
              ielem_LRC(ielem,iline) = ielem_Previous
              iline_LRC(ielem,iline) = iline_Previous
              exit
         endif

         !--- check for hitting bad data
         if (Skip_LRC_Mask(ielem_Next,iline_Next) == sym%YES) then
              ielem_LRC(ielem,iline) = ielem_Previous
              iline_LRC(ielem,iline) = iline_Previous
              exit
         endif

         !--- check for sign reversal
         if (Meander_Flag == sym%NO) then

          Grad_Temp = Grid_Data(ielem_Next,iline_Next) -  &
                      Grid_Data(ielem_Previous,iline_Previous)

          if (Grad_Flag * Grad_Temp < 0) then
              ielem_LRC(ielem,iline) = ielem_Previous
              iline_LRC(ielem,iline) = iline_Previous
              exit
          endif
         endif

         !--- check for saturation
         if (Grad_Flag == 1 .and. Grid_Data(ielem_Next,iline_Next) > Max_Grid_Data_Valid) then
              ielem_LRC(ielem,iline) = ielem_Next
              iline_LRC(ielem,iline) = iline_Next
              exit
         endif
         if (Grad_Flag == -1 .and. Grid_Data(ielem_Next,iline_Next) < Min_Grid_Data_Valid) then
              ielem_LRC(ielem,iline) = ielem_Next
              iline_LRC(ielem,iline) = iline_Next
              exit
         endif

         !--- store position
         ielem_Previous = ielem_Next
         iline_Previous = iline_Next

      enddo

    end do Line_Loop
  end do Element_Loop

end subroutine GRADIENT_MEANDER

!==============================================================
! subroutine COMPUTE_MEDIAN(z,mask,z_median,z_mean,z_std_median)
!
! Median filter
!==============================================================
subroutine COMPUTE_MEDIAN(z,mask,z_median,z_mean,z_std_median)

! The purpose of this function is to find 
! median (emed), minimum (emin) and maximum (emax)
! for the array elem with nelem elements. 

 real, dimension(:,:), intent(in):: z
 real, intent(out):: z_median
 real, intent(out):: z_mean
 real, intent(out):: z_std_median
 integer(kind=int1), dimension(:,:), intent(in):: mask 
 integer:: i,j,k,nx,ny,nelem
 real, dimension(:), allocatable::x
 real(kind=real4):: u

 z_median = missing_value_real4
 z_std_median = missing_value_real4
 z_mean = missing_value_real4

 nx = size(z,1)
 ny = size(z,2)

 nelem = nx * ny

 allocate(x(nelem))
 x = 0.0

 k = 0
 do i = 1, nx
   do j = 1, ny
      if (mask(i,j) == sym%NO .and. z(i,j) /= missing_value_real4) then
           k = k + 1   
           x(k) = z(i,j)
      endif
  enddo
 enddo

 nelem = k
   
 if (nelem < 1) then
     if (allocated(x)) deallocate(x)
     return
 endif 
!--- sort the array into ascending order
  do i=1,nelem-1
   do j=i+1,nelem
    if(x(j)<x(i))then
     u=x(j)
     x(j)=x(i)
     x(i)=u
    end if   
   end do
  end do

!---- pick the median
  if(mod(nelem,2)==1)then
   i=nelem/2+1
   z_median=x(i)
  else  
   i=nelem/2
   z_median=(x(i)+x(i+1))/2
   end if

!--- compute standard deviation wrt median
  z_mean = sum(x(1:nelem))/nelem
  z_std_median = sqrt(sum((x(1:nelem) - z_median)**2) / nelem)


! if (z_std_median > 60.0) then 
!         print *, "big std median ", z_std_median, nelem, x(1:nelem)
!         print *, "z_nxn = ", z
! endif

  if (allocated(x)) deallocate(x)

end subroutine COMPUTE_MEDIAN

!----------------------------------------------------------------------
! subroutine COMPUTE_MEDIAN_SEGMENT(z,mask,n,imin,imax,jmin,jmax,
!                                   z_median,z_std_median)
!
! Compute standard deviaion of an array wrt to the median
!----------------------------------------------------------------------
subroutine COMPUTE_MEDIAN_SEGMENT(z,mask,n,imin,imax,jmin,jmax, &
                                  z_median, &
                                  z_std_median)
  real(kind=real4), dimension(:,:), intent(in):: z
  integer(kind=int1), dimension(:,:), intent(in):: mask
  real(kind=real4), dimension(:,:), intent(out):: z_std_median
  real(kind=real4), dimension(:,:), intent(out):: z_median
! real(kind=real4), dimension(:,:), intent(out):: z_mean
  integer, intent(in):: n
  integer, intent(in):: imin
  integer, intent(in):: imax
  integer, intent(in):: jmin
  integer, intent(in):: jmax
  integer:: i
  integer:: j
  integer:: i1
  integer:: i2
  integer:: j1
  integer:: j2
  real(kind=real4) :: z_mean

  do i = imin, imax
    do j = jmin, jmax

     j1 = max(jmin,j-n)   !top index of local array
     j2 = min(jmax,j+n)   !bottom index of local array
     i1 = max(imin,i-n)   !left index of local array
     i2 = min(imax,i+n)   !right index of local array

     !--- compute median
     call COMPUTE_MEDIAN(z(i1:i2,j1:j2),mask(i1:i2,j1:j2),z_median(i,j), &
                         z_mean,z_std_median(i,j))

     enddo
  enddo

end subroutine COMPUTE_MEDIAN_SEGMENT

end module CX_SPATIAL_METRICS_MOD
