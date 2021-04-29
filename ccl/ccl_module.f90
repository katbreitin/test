! $Id: acha_cloud_cover_layers_module.f90 1523 2016-02-22 16:30:27Z wstraka $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: cloud_cover_layers_module.f90 (src)
!       CLOUD_COVER_LAYERS (module)
!
! PURPOSE: this module houses routines for computing cloud cover or
!          fraction in each atmospheric layer including the total
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
! Cloud Layer Definition
!
! Layered Cloud Layer Meanings
! 00000000 = Clear
! 00000001 = First Layer
! 00000010 = Second Layer
! 00000100 = Third Layer
! 00001000 = Fourth Layer
! 00010000 = Fifth Layer
! 00000101 = Third and First
! 00000110 = Third and Second
! 00000011 = Second and First
! 00000111 = High and Mid and Low
! etc, etc, etc
!
! to extract
! to see if first layer cloud, ibits(layer,0,1)
! to see if second layer cloud, ibits(layer,1,1)
! to see if third layer cloud, ibits(layer,2,1)
!
! Cloud_Cover_Layer_Profile = one bit integer = (0 to 100)
! Supercooled_Profile = one bit integer = (0 to 100)
!
!--------------------------------------------------------------------------------------
module CCL_MODULE

  use CCL_SERVICES_MOD, only : &
           real4, int1, int4, ccl_output_struct,ccl_symbol_struct, &
           ccl_input_struct, ccl_diag_struct, ccl_rtm_nwp_struct, &
           CCL_FETCH_PIXEL_NWP_RTM
  use NUMERICAL_ROUTINES_MOD, only: LOCATE

 implicit none
 public:: SETUP_CCL
 public:: DESTROY_CCL
 public:: COMPUTE_CLOUD_COVER_LAYERS
 private:: COMPUTE_BOX_WIDTH
 private:: KNOWING_P_COMPUTE_T_Z

 type(ccl_symbol_struct), private :: symbol

 integer, parameter:: N_Levels_NCEP = 4
 integer, parameter:: N_Levels_ISCCP = 4
 integer, parameter:: N_Levels_NOAT = 6
 real, dimension(N_Levels_NCEP), parameter, private:: CCL_Levels_NCEP =   (/1100.0,631.0,350.0,0.0/)  !hpa   
 real, dimension(N_Levels_ISCCP), parameter, private:: CCL_Levels_ISCCP = (/1100.0,680.0,440.0,0.0/)  !hpa 
 real, dimension(N_Levels_NOAT), parameter, private:: CCL_Levels_NOAT =   (/0.0,5000.0,10000.0,18000.0,24000.0,99000.0/) !kft

! corresponding altitude
 real, dimension(N_Levels_NCEP), parameter, private:: CCL_Levels_NCEP_Alt = (/0.0,12531.,26628.,99000.0/) !kft
 real, dimension(N_Levels_ISCCP), parameter, private:: CCL_Levels_ISCCP_Alt = (/0.0,10627.,21341.,99000.0/) !kft

 integer, save, private:: N_Levels
 integer, save, private:: N_Layers
 real, dimension(:), allocatable, save, private:: CCL_Levels
 real, dimension(:), allocatable, save, private:: CCL_Levels_Alt
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Num_Layer,Num_Layer_Conv,Num_Layer_Supercooled
!  integer:: Layer_Idx, Layer_Top_Idx, Layer_Lower_Idx,Layer_Base_Idx
!  real:: Z, Z_Base, Z_Lower
!ynoh (cira/csu) for ccl mode 3
  integer:: Layer_Idx, Layer_Top_Idx, Layer_Lower_Idx, Layer_Base_Idx, Layer_Lower_Base_Idx                            
  real:: Z, Z_Base, Z_Lower, Z_Lower_Base , Z_Lower_dummy

  real:: T, T_Base, T_Lower
  real::T_Freezing_Level,Zc_Temp
  integer:: Lev_Idx
  real:: BT11, BT67, Emiss_Tropo_11,TSfc
  real, parameter:: Thresh1 = 210.0
  real, parameter:: Emiss_Thre = 0.63
  real, parameter:: Delta_Z_Thre = 1000.

  type(ccl_rtm_nwp_struct) :: CCL_RTM_NWP

  integer, private, parameter:: N_Sc_Lut = 20
  real(kind=real4), dimension(N_Sc_Lut), private, parameter:: &
  Sc_Tc_Lut = (/202.21,206.77,211.33,215.88,220.44,225.00,229.56,234.11,238.67,243.23, &
                247.79,252.34,256.90,261.46,266.01,270.57,275.13,279.69,284.24,288.80/)

  real(kind=real4), dimension(N_Sc_Lut), private, parameter:: &
  Sc_Prob_Lut = (/0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.002,0.051,0.243, &
                  0.470,0.658,0.800,0.886,0.888,0.870,0.142,0.004,0.004,0.004 /)

 include 'ccl_parameters.inc'

 contains
!------------------------------------------------------------------------------
! Setup CCL parameters based on CCL_TYPE
!------------------------------------------------------------------------------
subroutine SETUP_CCL(CCL_Type) 

   integer, intent(in) :: CCL_Type 

   select case (CCL_Type)

    case (2)    !NCEP
      N_Levels = N_Levels_NCEP
      allocate(CCL_Levels(N_Levels))
      CCL_Levels = CCL_Levels_NCEP
      allocate(CCL_Levels_Alt(N_Levels))
      CCL_Levels_Alt = CCL_Levels_NCEP_Alt

    case (1)    !ISCCP
      N_Levels = N_Levels_ISCCP
      allocate(CCL_Levels(N_Levels))
      CCL_Levels = CCL_Levels_ISCCP
      allocate(CCL_Levels_Alt(N_Levels))
      CCL_Levels_Alt = CCL_Levels_ISCCP_Alt

    case (0)  !NOAT
      N_Levels = N_Levels_NOAT
      allocate(CCL_Levels(N_Levels))
      CCL_Levels = CCL_Levels_NOAT
      allocate(CCL_Levels_Alt(N_Levels))
      CCL_Levels_Alt = CCL_Levels_NOAT

  end select 

  N_Layers = N_Levels - 1

end subroutine SETUP_CCL

!------------------------------------------------------------------------------
! Destroy CCL parameters based on CCL_TYPE
!------------------------------------------------------------------------------
subroutine DESTROY_CCL() 

  if (allocated(CCL_Levels)) deallocate(CCL_Levels)
  if (allocated(CCL_Levels_Alt)) deallocate(CCL_Levels_Alt)

end subroutine DESTROY_CCL

!------------------------------------------------------------------------------
! compute cloud fraction over a nxn array using the Bayesian probability
!------------------------------------------------------------------------------
 subroutine COMPUTE_CLOUD_COVER_LAYERS(Input, Symbol_In, Output, Diag)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================

  type(ccl_input_struct), intent(inout) :: Input
  type(ccl_symbol_struct), intent(inout) :: Symbol_In
  type(ccl_output_struct), intent(inout) :: Output
  type(ccl_diag_struct), intent(inout), optional :: Diag
  integer, save:: Diag_Warning_Flag = 0


  integer:: Num_Elems
  integer:: Num_Lines
  integer :: i,j                    !pixel indices
  integer :: i1,i2,j1,j2            !pixel indices of ccl box
  integer :: i11,i22,j11,j22        !pixels inices of area which ccl result is valid
  integer :: ni, nj                 !size of ccl box
  integer:: N                       !width of ccl box
  integer:: M                       !width of ccl region
  integer :: Num_Good
  integer :: Num_Cloud
  integer:: Num_High, Num_Mid, Num_Low, Num_Clear, Num_All
  real:: Clear_Fraction

  real (kind=4), dimension(:,:), allocatable:: Pixel_Uncertainty
  integer (kind=1), dimension(:,:), allocatable:: pos, len      !ibits arguments

!  logical:: USE_BASE, USE_LOWER
  logical:: USE_BASE, USE_LOWER, USE_LOWER_BASE

  real (kind=real4), dimension(:,:), allocatable:: Alt_Freezing
  real, dimension(:,:), allocatable:: Lapse_Rate
  real, dimension(:,:,:), allocatable:: T_Laytop
  real, dimension(:,:,:), allocatable:: Sc_Prob_Lay_Pix
  real (kind=4), dimension(:,:), allocatable:: Sc_Prob_Pix
  real (kind=4), dimension(:,:), allocatable:: Sc_Prob_Base_Pix
  integer, dimension(:,:), allocatable:: Mask_Missing_Sc
  integer, dimension(:,:,:), allocatable:: Mask_Sc_Layer
  integer, dimension(:,:), allocatable:: Mask_Missing_Conv
  real (kind=4), dimension(:,:), allocatable:: Conv_Prob_Pix
  real:: Sc_Prob_Lower

  USE_BASE = .false.
  USE_LOWER = .false.
  if (Input%CCL_Mode >= 2) USE_BASE = .true.
  if (Input%CCL_Mode == 3) USE_LOWER = .true.
!ynoh (cira/csu) for ccl mode 3  (to use the lower base)
  USE_LOWER_BASE = .true.
  !-------------------------------------------------------------------------------
  ! Total Cloud Fraction and Its Uncertainty
  !-------------------------------------------------------------------------------

  !--- copy input symbol to a module-wide variable
  symbol = Symbol_In

  !--- initialize
  Output%Total_Cloud_Fraction = MISSING_VALUE_REAL4
  Output%Total_Cloud_Fraction_Uncer = MISSING_VALUE_REAL4
  Output%Supercooled_Total_Fraction = MISSING_VALUE_REAL4
  Output%Conv_Total_Fraction = MISSING_VALUE_REAL4

  Output%Cloud_Layer = MISSING_VALUE_INTEGER1
  Output%Supercooled_Cloud_Layer = MISSING_VALUE_INTEGER1
  Output%Conv_Cloud_Layer = MISSING_VALUE_INTEGER1
  Output%QF = MISSING_VALUE_INTEGER1

  !--- initialize diagnostic output
  if (present(Diag) .and. Diag_Warning_Flag == 0) then
      print *, "CLAVR-x / CCL ===>  Diagnostic Output Turned On"
      Diag_Warning_Flag = 1
  endif
  if (present(Diag)) Diag%Array_1 = Missing_Value_Real4
  if (present(Diag)) Diag%Array_2 = Missing_Value_Real4
  if (present(Diag)) Diag%Array_3 = Missing_Value_Real4

  !--- Determine Box Width 
  call COMPUTE_BOX_WIDTH(Input%Sensor_Resolution_KM,CCL_BOX_WIDTH_KM, N)
  call COMPUTE_BOX_WIDTH(Input%Sensor_Resolution_KM,CCL_SPACING_KM, M)

!==============================================================================
!
! Compute Cloud Cover Layers which is defined here as the fraction of
! high, middle and low-level cloud in a NxN array centered on each pixel
!
! output (through global arrays)
!  Output%Total_Cloud_Fraction = cloud fraction (0-1) for all clouds over array
!  Output%Total_Cloud_Fraction_Uncer = cloud fraction uncertainty (0-1)
!  Output%Supercooled_Total_Fraction = cloud fraction (0-1) for supercooled cloud over array
!  Output%Conv_Total_Fraction = cloud fraction (0-1) for convective cloud over array
!  Output%???_Fraction_Layer?= cloud fraction for individual layer (0-100) (integer)
!
!  Notes
!  1. cloud fractions have values limited by 1/N^2. If N=1, values are 1/9,2/9...
!  2. Cloud_Fraction is the total cloud amount.  This is computed from the
!     mask only.  If arrays have pixels with failed ACHA, H+M+L /= Total.
!  3. H/M/L boundaries in ccl_parameters.inc
!==============================================================================

  !--------------------------------------------------------------------
  ! compute pixel-level cloud layer flag and H/M/L masks
  !--------------------------------------------------------------------
  Num_Elems = Input%Number_of_Elements
  Num_Lines = Input%Number_of_Lines
  allocate(Pixel_Uncertainty(Num_Elems, Num_Lines))
  allocate(Sc_Prob_Pix(Num_Elems, Num_Lines))
  allocate(Sc_Prob_Base_Pix(Num_Elems, Num_Lines))
  allocate(Mask_Missing_Sc(Num_Elems, Num_Lines))
  allocate(Mask_Sc_Layer(Num_Elems, Num_Lines,5))
  allocate(Lapse_Rate(Num_Elems, Num_Lines))
  allocate(T_Laytop(Num_Elems, Num_Lines,5))
  allocate(Sc_Prob_Lay_Pix(Num_Elems, Num_Lines,5))
  allocate(Alt_Freezing(Num_Elems, Num_Lines))
  allocate(Conv_Prob_Pix(Num_Elems, Num_Lines))
  allocate(Mask_Missing_Conv(Num_Elems, Num_Lines))

  !--- make cloud fraction pixel level uncertainty
  Pixel_Uncertainty = MISSING_VALUE_REAL4
  where(Input%Cloud_Probability >= 0.5)
     Pixel_Uncertainty = 1.0 - Input%Cloud_Probability
  endwhere
  where(Input%Cloud_Probability < 0.5)
     Pixel_Uncertainty = Input%Cloud_Probability
  endwhere

  Sc_Prob_Pix = MISSING_VALUE_REAL4
  Sc_Prob_Base_Pix = MISSING_VALUE_REAL4
  Mask_Missing_Sc = 1  ! 1 is no data
  Mask_Sc_Layer = 1
  Mask_Missing_Conv = 1

  Conv_Prob_Pix = MISSING_VALUE_REAL4

!---- Alt should be an input
! call COMPUTE_ALTITUDE_FROM_PRESSURE_CCL(1,Num_Lines,Num_Elems,symbol%YES,Input%Invalid_Data_Mask,Input%Pc,Alt)
! call COMPUTE_ALTITUDE_FROM_PRESSURE_CCL(1,Num_Lines,Num_Elems,symbol%YES,Input%Invalid_Data_Mask,Input%Pc_Base,Alt_Base)
! call COMPUTE_ALTITUDE_FROM_PRESSURE_CCL(1,Num_Lines,Num_Elems,symbol%YES,Input%Invalid_Data_Mask,Input%Pc_Lower,Alt_Lower)

! convert freezing level pressure to altitude so lapse rate inside clouds can be computed
  call COMPUTE_ALTITUDE_FROM_PRESSURE_CCL(1,Num_Lines,Num_Elems,symbol%YES,Input%Invalid_Data_Mask,Input%Freezing_level_pressure,Alt_Freezing)

!=========================================
! Supercooled section
  Lapse_Rate = MISSING_VALUE_REAL4
  T_Laytop = MISSING_VALUE_REAL4

  line_loop_sc: do Line_Idx = 1, Num_Lines
    elem_loop_sc: do Elem_Idx = 1, Num_Elems

      T_Freezing_Level = MISSING_VALUE_REAL4
      Zc_Temp = MISSING_VALUE_REAL4

      if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) cycle

      call NULL_PIX_POINTERS(Input, CCL_RTM_NWP)
      call CCL_FETCH_PIXEL_NWP_RTM(Input, Symbol, &
                                   Elem_Idx,Line_Idx, CCL_RTM_NWP)

      T= Input%Tc(Elem_Idx,Line_Idx)
      T_Base = Input%Tc_Base(Elem_Idx,Line_Idx)

      !--- Compute supercooled prob and compute lapse rate inside cloud if base is available assuming linear
!      if (T /= MISSING_VALUE_REAL4 .and. T < 273.15 .and. T < T_Base) then
      if (T /= MISSING_VALUE_REAL4 .and. T < 273.15) then
         call SUPERCOOLED_CLOUD_PROBABILITY_local(T, Sc_Prob_Pix(Elem_Idx,Line_Idx))
      ! Only compute the part between top and freezing level if base is warmer
      ! than freezing temperature (NOT NECESSARY?)
         !if (T_Base /= MISSING_VALUE_REAL4 .and. T_Base < 273.15) then
            call SUPERCOOLED_CLOUD_PROBABILITY_local(T_Base, Sc_Prob_Base_Pix(Elem_Idx,Line_Idx))
            Lapse_Rate(Elem_Idx,Line_Idx) = (T - T_Base)/(Input%Alt(Elem_Idx,Line_Idx) - Input%Alt_Base(Elem_Idx,Line_Idx))
         !else
         !   call KNOWING_P_COMPUTE_T_Z(CCL_RTM_NWP,Input%Freezing_Level_Pressure(Elem_Idx,Line_Idx),T_Freezing_Level,Zc_Temp,Lev_Idx)
         !   Lapse_Rate(Elem_Idx,Line_Idx) = (T-T_Freezing_Level)/(Input%Alt(Elem_Idx,Line_Idx) - Alt_Freezing(Elem_Idx,Line_Idx) )
         !endif
      endif

      if (Sc_Prob_Pix(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) then
         Mask_Missing_Sc(Elem_Idx,Line_Idx) = 0
      endif

      call NULL_PIX_POINTERS(Input, CCL_RTM_NWP)

    end do elem_loop_sc
  end do line_loop_sc
!=========================================

!================ Conv section
  line_loop_conv: do Line_Idx = 1, Num_Lines
    elem_loop_conv: do Elem_Idx = 1, Num_Elems
     if ((Input%Free_Convection_Height(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) .and. &
         (Input%Zc(Elem_Idx,Line_Idx) - Input%Free_Convection_Height(Elem_Idx,Line_Idx) >= Delta_Z_Thre) .and. &
         (Input%Emiss_Tropo_11(Elem_Idx,Line_Idx) >= Emiss_Thre) ) then
        Conv_Prob_Pix(Elem_Idx,Line_Idx) = 1.
     endif
  
     if (Input%Chan_On_110um == Symbol%YES) then

       BT11 = Input%Bt_110um(Elem_Idx,Line_Idx)
       Tsfc = Input%Tsfc(Elem_Idx,Line_Idx)
       Emiss_Tropo_11 = Input%Emiss_Tropo_11(Elem_Idx,Line_Idx)

       if (BT11 < Thresh1) Conv_Prob_Pix(Elem_Idx,Line_Idx) = 1.

       if (Emiss_Tropo_11 >= 0.95) Conv_Prob_Pix(Elem_Idx,Line_Idx) = 1.

       if (Input%Chan_On_067um == Symbol%YES) then
          BT67 = Input%Bt_067um(Elem_Idx,Line_Idx)
          if (Emiss_Tropo_11 >= 0.9 .and. (BT67 - Bt11) >= 0.9) then
             Conv_Prob_Pix(Elem_Idx,Line_Idx) = 1.
          endif
       endif

       if (Tsfc - BT11 < 30.0) Conv_Prob_Pix(Elem_Idx,Line_Idx) = 0.0
     endif

      if (Conv_Prob_Pix(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) then
         Mask_Missing_Conv(Elem_Idx,Line_Idx) = 0
      endif

    end do elem_loop_conv
  end do line_loop_conv
!================


!ynoh (cira/csu) for ccl mode 3
! Previously ACHA%Alt and ACHA%Base_Alt are only computed in process_clavrx.f90

  !--------------------------------------------------------------------
  !Layer_Top_Idx <= N_Layers Compute cloud_layer flag 
  !--------------------------------------------------------------------
  line_loop: do Line_Idx = 1, Num_Lines
    elem_loop: do Elem_Idx = 1, Num_Elems

      !--- select the vertical coordinates based on CCL_Type
      select case (Input%CCL_Type) 

        case (1,2)
          Z = Input%Pc(Elem_Idx,Line_Idx)
          Z_Base = Input%Pc_Base(Elem_Idx,Line_Idx)
          Z_Lower = Input%Pc_Lower(Elem_Idx,Line_Idx)
!ynoh (cira/csu) for ccl mode 3
          Z_Lower_Base = Input%Pc_Lower_Base(Elem_Idx,Line_Idx)
          Z_Lower_dummy = max(Z_Lower, Z_Base)

        case (0)
          Z = Input%Alt(Elem_Idx,Line_Idx)
          Z_Base = Input%Alt_Base(Elem_Idx,Line_Idx)
!ynoh (cira/csu) for ccl mode 3  (three added)
          Z_Lower = Input%Alt_Lower(Elem_Idx,Line_Idx)
          Z_Lower_Base = Input%Alt_Lower_Base(Elem_Idx,Line_Idx)
          Z_Lower_dummy = min(Z_Lower, Z_Base)

      end select 

      if (Z /= MISSING_VALUE_REAL4) Z = max(0.0,Z)
      if (Z_Base /= MISSING_VALUE_REAL4) Z_Base = max(0.0,Z_Base)
      if (Z_Lower /= MISSING_VALUE_REAL4) Z_Lower = max(0.0,Z_Lower)
!ynoh (cira/csu) for ccl mode 3
      if (Z_Lower_Base /= MISSING_VALUE_REAL4) Z_Lower_Base = max(0.0,Z_Lower_Base)

      Layer_Top_Idx = -1
      Layer_Base_Idx = -1
      Layer_Lower_Idx = -1
!ynoh (cira/csu) for ccl mode 3
      Layer_Lower_Base_Idx = -1

      !----- check for valid cloud vertical coordinate
      if (Z == MISSING_VALUE_REAL4) then
          cycle
      endif

      !-- need to intialize pixels with valid data to 0
      Output%Cloud_Layer(Elem_Idx,Line_Idx) = 0
      Output%Supercooled_Cloud_Layer(Elem_Idx,Line_Idx) = 0
      Output%Conv_Cloud_Layer(Elem_Idx,Line_Idx) = 0

      !---------------------------------------------------------------------
      ! Find Layers
      !---------------------------------------------------------------------

      !--- find top of highest cloud level index
      if (Z /= MISSING_VALUE_REAL4) then
         do Layer_Top_Idx = 1, N_Layers
            if ((Z >= CCL_Levels(Layer_Top_Idx) .and. Z <= CCL_Levels(Layer_Top_Idx + 1)) .or. &
                (Z <= CCL_Levels(Layer_Top_Idx) .and. Z >= CCL_Levels(Layer_Top_Idx + 1))) exit
         enddo
      endif

      !--- find base of highest cloud level index
      if (USE_BASE .and. Z_Base /= MISSING_VALUE_REAL4) then
        do Layer_Base_Idx = 1, N_Layers
         if ((Z_Base >= CCL_Levels(Layer_Base_Idx) .and. Z_Base <= CCL_Levels(Layer_Base_Idx + 1)) .or. &
             (Z_Base <= CCL_Levels(Layer_Base_Idx) .and. Z_Base >= CCL_Levels(Layer_Base_Idx + 1))) exit
        enddo
      endif

      !--- find top of lowest cloud level index
!ynoh (cira/csu) for ccl mode 3
!-- use Z_Lower_dummy instead of Z_Lower to avoid double counting layers in case Z_Lower > Z_Base 
      if (USE_LOWER .and. Z_Lower /= MISSING_VALUE_REAL4 .and. Z_Lower_dummy /= MISSING_VALUE_REAL4 ) then
        do Layer_Lower_Idx = 1, N_Layers
         if ((Z_Lower_dummy >= CCL_Levels(Layer_Lower_Idx) .and. Z_Lower_dummy <= CCL_Levels(Layer_Lower_Idx + 1)) .or. &
             (Z_Lower_dummy <= CCL_Levels(Layer_Lower_Idx) .and. Z_Lower_dummy >= CCL_Levels(Layer_Lower_Idx + 1))) exit
        enddo
      endif

!ynoh (cira/csu) for ccl mode 3
      !--- find base of lowest cloud level index (This only works for type 0)
      if (USE_LOWER .and. USE_LOWER_BASE .and. &
          Z_Lower_Base /= MISSING_VALUE_REAL4 .and. Z_Lower_Base < Z_Lower_dummy ) then
        do Layer_Lower_Base_Idx = 1, N_Layers
         if ((Z_Lower_Base >= CCL_Levels(Layer_Lower_Base_Idx) .and. Z_Lower_Base <= CCL_Levels(Layer_Lower_Base_Idx + 1)) .or. &
             (Z_Lower_Base <= CCL_Levels(Layer_Lower_Base_Idx) .and. Z_Lower_Base >= CCL_Levels(Layer_Lower_Base_Idx + 1))) exit
        enddo
      endif

      !---------------------------------------------------------------------
      ! Compute Cloud Layer
      !---------------------------------------------------------------------
      if (Layer_Top_Idx > 0 .and. Layer_Top_Idx <= N_Layers) then
          Output%Cloud_Layer(Elem_Idx,Line_Idx) = ibset(Output%Cloud_Layer(Elem_Idx,Line_Idx),Layer_Top_Idx-1)

        if (Sc_Prob_Pix(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) then
          if (Sc_Prob_Pix(Elem_Idx,Line_Idx) >= 0.5) Output%Supercooled_Cloud_Layer(Elem_Idx,Line_Idx) = &
                                                         ibset(Output%Supercooled_Cloud_Layer(Elem_Idx,Line_Idx),Layer_Top_Idx-1)
           Mask_Sc_Layer(Elem_Idx,Line_Idx,Layer_Top_Idx) = 0
           Sc_Prob_Lay_Pix(Elem_Idx,Line_Idx,Layer_Top_Idx) = Sc_Prob_Pix(Elem_Idx,Line_Idx)
        endif

        if (Conv_Prob_Pix(Elem_Idx,Line_Idx) == 1.0) then
           Output%Conv_Cloud_Layer(Elem_Idx,Line_Idx) = ibset(Output%Conv_Cloud_Layer(Elem_Idx,Line_Idx),Layer_Top_Idx-1)
        endif
      endif

      !--- account for base
      do Layer_Idx = 1, N_Layers
        if (USE_BASE .and. Layer_Base_Idx <= Layer_Top_Idx .and. Layer_Idx >= Layer_Base_Idx .and. Layer_Top_Idx >= Layer_Idx) then
           Output%Cloud_Layer(Elem_Idx,Line_Idx) = ibset(Output%Cloud_Layer(Elem_Idx,Line_Idx),Layer_Idx-1)
        endif

        if (USE_BASE .and. Layer_Base_Idx <= Layer_Top_Idx .and. Layer_Idx >= Layer_Base_Idx .and. Layer_Top_Idx >= Layer_Idx .and. &
            Conv_Prob_Pix(Elem_Idx,Line_Idx) == 1.0) then
           Output%Conv_Cloud_Layer(Elem_Idx,Line_Idx) = ibset(Output%Conv_Cloud_Layer(Elem_Idx,Line_Idx),Layer_Idx-1)
        endif
      enddo 

      T= Input%Tc(Elem_Idx,Line_Idx)
      do Layer_idx = 1,N_Layers
       if (USE_BASE .and. Layer_Base_Idx < Layer_Top_Idx .and. Layer_Idx >= Layer_Base_Idx .and. Layer_Top_Idx > Layer_Idx .and. &
           Lapse_Rate(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then

   ! Don't need to find out if base is warmer than freezing pont, if it does,
   ! the probability is small for the warm part
   !      if (T_Base /= MISSING_VALUE_REAL4 .and. T_Base < 273.15) then
   !         call SUPERCOOLED_CLOUD_PROBABILITY_local(T_Base,
   !         Sc_Prob_Base_Pix(Elem_Idx,Line_Idx))
   !         Lapse_Rate(Elem_Idx,Line_Idx) = (T - T_Base)/(
   !         Input%Alt(Elem_Idx,Line_Idx) - Input%Alt_Base(Elem_Idx,Line_Idx))
   !      else

            T_Laytop(Elem_Idx,Line_Idx,Layer_idx) = T-Lapse_Rate(Elem_Idx,Line_Idx)*(Input%Alt(Elem_Idx,Line_Idx)-CCL_Levels_Alt(Layer_idx+1))

            call SUPERCOOLED_CLOUD_PROBABILITY_local(T_Laytop(Elem_Idx,Line_Idx,Layer_idx), &
                 Sc_Prob_Lay_Pix(Elem_Idx,Line_Idx,Layer_idx))
            !print *, 'Layer/Elem/Line ',Layer_idx,Elem_Idx,Line_Idx,' T_Laytop
            !',T_Laytop(Elem_Idx,Line_Idx,Layer_idx)
            !print *, 'Sc_Prob_Lay_pix ',
            !Sc_Prob_Lay_pix(Elem_Idx,Line_Idx,Layer_idx)

            if (Sc_Prob_Lay_Pix(Elem_Idx,Line_Idx,Layer_Idx) >= 0.5) Output%Supercooled_Cloud_Layer(Elem_Idx,Line_Idx) = &
                                                         ibset(Output%Supercooled_Cloud_Layer(Elem_Idx,Line_Idx),Layer_Idx-1)

            if (Sc_Prob_Lay_Pix(Elem_Idx,Line_Idx,Layer_idx) /= MISSING_VALUE_REAL4) then
                Mask_Sc_Layer(Elem_Idx,Line_Idx,Layer_Idx) = 0
            endif

       endif
      enddo

      !--- account for lower cloud
!      if (USE_LOWER .and. Layer_Lower_Idx > 0 .and. Layer_Lower_Idx <= N_Layers) then
!ynoh (cira/csu) for ccl mode 3
!-- to avoid double counting layers in case Z_Lower > Z_Base 
      if (USE_LOWER .and. Layer_Lower_Idx > 0 .and. Layer_Lower_Idx <= N_Layers .and. Layer_Lower_Idx < Layer_Base_Idx) then
          Output%Cloud_Layer(Elem_Idx,Line_Idx) = ibset(Output%Cloud_Layer(Elem_Idx,Line_Idx),Layer_Lower_Idx-1)

          ! supercooled
          T_Lower= Input%Tc_Lower(Elem_Idx,Line_Idx)
          Sc_Prob_Lower = MISSING_VALUE_REAL4

          if (T_Lower /= MISSING_VALUE_REAL4 .and. T_Lower < 273.15) call SUPERCOOLED_CLOUD_PROBABILITY_local(T_Lower, Sc_Prob_Lower)

          if (Sc_Prob_Lower /= MISSING_VALUE_REAL4) then
              if (Sc_Prob_Lower >= 0.5) Output%Supercooled_Cloud_Layer(Elem_Idx,Line_Idx) = &
                                                             ibset(Output%Supercooled_Cloud_Layer(Elem_Idx,Line_Idx),Layer_Lower_Idx-1)
              Mask_Sc_Layer(Elem_Idx,Line_Idx,Layer_Lower_Idx) = 0
              Sc_Prob_Lay_Pix(Elem_Idx,Line_Idx,Layer_Lower_Idx) = Sc_Prob_Lower
          endif
 
          ! convective, assign lower level same as the upper level
           if (Conv_Prob_Pix(Elem_Idx,Line_Idx) == 1.0) then
               Output%Conv_Cloud_Layer(Elem_Idx,Line_Idx) = ibset(Output%Conv_Cloud_Layer(Elem_Idx,Line_Idx),Layer_Lower_Idx-1)
           endif
      endif


!ynoh (cira/csu) for ccl mode 3
      !--- account for lower cloud base
      do Layer_Idx = 1, N_Layers
        if (USE_LOWER .and. USE_LOWER_BASE .and. Layer_Lower_Base_Idx <= Layer_Lower_Idx .and. &
            Layer_Idx >= Layer_Lower_Base_Idx .and. &
            Layer_Lower_Idx >= Layer_Idx .and. Layer_Idx < Layer_Base_Idx) then
          Output%Cloud_Layer(Elem_Idx,Line_Idx) = ibset(Output%Cloud_Layer(Elem_Idx,Line_Idx),Layer_Lower_Idx-1)
        endif
      enddo 

      !print *, "Test ", Z, Layer_Top_Idx, Z_Base, Layer_Base_Idx, Output%Cloud_Layer(Elem_Idx,Line_Idx)

!      write(unit=6,fmt="(I3, I3, I08,B08)") Layer_Top_Idx, Layer_Base_Idx, Output%Cloud_Layer(Elem_Idx,Line_Idx), Output%Cloud_Layer(Elem_Idx,Line_Idx)

    end do elem_loop
  end do line_loop

   where (Input%Cloud_Mask == Symbol%CLEAR .or. Input%Cloud_Mask == Symbol%PROB_CLEAR)
          Output%Cloud_Layer = 0   !0000
   endwhere 

! set QF based on ACHA QF
   where (Input%ACHA_QF == CTH_DQF_GOOD_RETREVIAL ) Output%QF = CCL_DQF_GOOD_RETREVIAL
   where (Input%ACHA_QF == CTH_DQF_MARGINAL_RETREVIAL .or. Input%ACHA_Qf == CTH_DQF_RETREVIAL_ATTEMPTED) &
         Output%QF = CCL_DQF_DEGRADED_RETREVIAL
   where (Input%ACHA_QF == CTH_DQF_BAD_RETREVIAL ) Output%QF = CCL_DQF_BAD_RETREVIAL


 !--------------------------------------------------------------------
 ! compute pixel-level cloud cover for each layer over the box
 ! N = ccl box size
 ! M = ccl result spacing
 !--------------------------------------------------------------------

 line_loop_cover: do j = 1, Num_Lines, 2*M+1

    j1 = max(1,j-N)
    j2 = min(Num_Lines,j+N)
    j11 = max(1,j-M)
    j22 = min(Num_Lines,j+M)

    element_loop_cover: do i = 1, Num_Elems, 2*M+1

      i1 = max(1,i-N)
      i2 = min(Num_Elems,i+N)
      i11 = max(1,i-M)
      i22 = min(Num_Elems,i+M)

      ni = i2 - i1 + 1 
      nj = j2 - j1 + 1 

      if (allocated(pos)) deallocate(pos)
      if (allocated(len)) deallocate(len)

      allocate(pos(ni,nj), len(ni,nj))
      len = 1

      !--- check for a bad pixel pixel
      if (Input%Invalid_Data_Mask(i,j) == Symbol%YES) cycle

      !---- new   layers are inbetween levels
      Num_Clear = count(Input%Cloud_Probability(i1:i2,j1:j2) < 0.5 .and. &
                        Input%Cloud_Probability(i1:i2,j1:j2) /= MISSING_VALUE_REAL4)         
      Num_Cloud = count(Input%Cloud_Probability(i1:i2,j1:j2) >= 0.5)        
      Num_Good = count(Input%Cloud_Probability(i1:i2,j1:j2) /= MISSING_VALUE_REAL4)        
      Num_All = Num_Cloud + Num_Clear

      !--- see if there are any valid mask points, if not skip this pixel
      if (Num_Good < COUNT_MIN_CCL) then
         cycle
      endif

      !--- Total Cloud Fraction
      Output%Total_Cloud_Fraction(i11:i22,j11:j22) = real(Num_Cloud)/real(Num_Good)

      !--- compute the uncertainty of the total cloud fraction
      Output%Total_Cloud_Fraction_Uncer(i11:i22,j11:j22) = sum(Pixel_Uncertainty(i1:i2,j1:j2))/real(Num_Good)

      Output%Supercooled_Total_Fraction(i11:i22,j11:j22) = sum(Sc_Prob_Pix(i1:i2,j1:j2)*(1.0-Mask_Missing_SC(i1:i2,j1:j2)))/real(Num_Good) 

      Output%Conv_Total_Fraction(i11:i22,j11:j22) = sum(Conv_Prob_Pix(i1:i2,j1:j2)*(1.0-Mask_Missing_Conv(i1:i2,j1:j2)))/real(Num_Good)

      !--- see if there are any valid CCL points, if not skip this pixel
      if (Num_All == 0 ) then
         cycle
      endif

      do Layer_Idx = 1, N_Layers
          Num_Layer = int(sum(real(ibits(Output%Cloud_Layer(i1:i2,j1:j2),Layer_Idx-1,len))))
          if (Layer_Idx == 1) Output%Fraction_Layer1(i11:i22,j11:j22) = int(100.0*real(Num_Layer)/real(Num_All))
          if (Layer_Idx == 2) Output%Fraction_Layer2(i11:i22,j11:j22) = int(100.0*real(Num_Layer)/real(Num_All))
          if (Layer_Idx == 3) Output%Fraction_Layer3(i11:i22,j11:j22) = int(100.0*real(Num_Layer)/real(Num_All))
          if (Layer_Idx == 4) Output%Fraction_Layer4(i11:i22,j11:j22) = int(100.0*real(Num_Layer)/real(Num_All))
          if (Layer_Idx == 5) Output%Fraction_Layer5(i11:i22,j11:j22) = int(100.0*real(Num_Layer)/real(Num_All))

! Supercooled
         if (Layer_Idx == 1) Output%Supercooled_Fraction_Layer1(i11:i22,j11:j22) = &
                int(100.0*sum(Sc_Prob_Pix(i1:i2,j1:j2)*(1.0-Mask_Sc_Layer(i1:i2,j1:j2,Layer_Idx)))/real(Num_Good))
         if (Layer_Idx == 2) Output%Supercooled_Fraction_Layer2(i11:i22,j11:j22) = &
                int(100.0*sum(Sc_Prob_Pix(i1:i2,j1:j2)*(1.0-Mask_Sc_Layer(i1:i2,j1:j2,Layer_Idx)))/real(Num_Good))
         if (Layer_Idx == 3) Output%Supercooled_Fraction_Layer3(i11:i22,j11:j22) = &
                int(100.0*sum(Sc_Prob_Pix(i1:i2,j1:j2)*(1.0-Mask_Sc_Layer(i1:i2,j1:j2,Layer_Idx)))/real(Num_Good))
         if (Layer_Idx == 4) Output%Supercooled_Fraction_Layer4(i11:i22,j11:j22) = &
                int(100.0*sum(Sc_Prob_Pix(i1:i2,j1:j2)*(1.0-Mask_Sc_Layer(i1:i2,j1:j2,Layer_Idx)))/real(Num_Good))
         if (Layer_Idx == 5) Output%Supercooled_Fraction_Layer5(i11:i22,j11:j22) = &
                int(100.0*sum(Sc_Prob_Pix(i1:i2,j1:j2)*(1.0-Mask_Sc_Layer(i1:i2,j1:j2,Layer_Idx)))/real(Num_Good))

         !Num_Layer_Supercooled = int(sum(real(ibits(Output%Supercooled_Cloud_Layer(i1:i2,j1:j2),Layer_Idx-1,len))))
         !if (Layer_Idx == 1) Output%Supercooled_Fraction_Layer1(i11:i22,j11:j22) = int(100*real(Num_Layer_Supercooled)/real(Num_all))
         !if (Layer_Idx == 2) Output%Supercooled_Fraction_Layer2(i11:i22,j11:j22) = int(100*real(Num_Layer_Supercooled)/real(Num_all))
         !if (Layer_Idx == 3) Output%Supercooled_Fraction_Layer3(i11:i22,j11:j22) = int(100*real(Num_Layer_Supercooled)/real(Num_all))
         !if (Layer_Idx == 4) Output%Supercooled_Fraction_Layer4(i11:i22,j11:j22) = int(100*real(Num_Layer_Supercooled)/real(Num_all))
         !if (Layer_Idx == 5) Output%Supercooled_Fraction_Layer5(i11:i22,j11:j22) = int(100*real(Num_Layer_Supercooled)/real(Num_all))

         Num_Layer_Conv = int(sum(real(ibits(Output%Conv_Cloud_Layer(i1:i2,j1:j2),Layer_Idx-1,len))))
       
         if (Layer_Idx == 1) Output%Conv_Fraction_Layer1(i11:i22,j11:j22) = int(100*real(Num_Layer_Conv)/real(Num_all))
         if (Layer_Idx == 2) Output%Conv_Fraction_Layer2(i11:i22,j11:j22) = int(100*real(Num_Layer_Conv)/real(Num_all))
         if (Layer_Idx == 3) Output%Conv_Fraction_Layer3(i11:i22,j11:j22) = int(100*real(Num_Layer_Conv)/real(Num_all))
         if (Layer_Idx == 4) Output%Conv_Fraction_Layer4(i11:i22,j11:j22) = int(100*real(Num_Layer_Conv)/real(Num_all))
         if (Layer_Idx == 5) Output%Conv_Fraction_Layer5(i11:i22,j11:j22) = int(100*real(Num_Layer_Conv)/real(Num_all))

        if (USE_BASE .or. USE_LOWER) then
         if (Layer_Idx == 1) Output%Supercooled_Fraction_Layer1(i11:i22,j11:j22) = &
               int(100.0*sum(Sc_Prob_Lay_Pix(i1:i2,j1:j2,Layer_idx)*(1.0-Mask_Sc_Layer(i1:i2,j1:j2,Layer_Idx)))/real(Num_Good))
         if (Layer_Idx == 2) Output%Supercooled_Fraction_Layer2(i11:i22,j11:j22) = &
               int(100.0*sum(Sc_Prob_Lay_Pix(i1:i2,j1:j2,Layer_idx)*(1.0-Mask_Sc_Layer(i1:i2,j1:j2,Layer_Idx)))/real(Num_Good))
         if (Layer_Idx == 3) Output%Supercooled_Fraction_Layer3(i11:i22,j11:j22) = &
               int(100.0*sum(Sc_Prob_Lay_Pix(i1:i2,j1:j2,Layer_idx)*(1.0-Mask_Sc_Layer(i1:i2,j1:j2,Layer_Idx)))/real(Num_Good))
         if (Layer_Idx == 4) Output%Supercooled_Fraction_Layer4(i11:i22,j11:j22) = &
               int(100.0*sum(Sc_Prob_Lay_Pix(i1:i2,j1:j2,Layer_idx)*(1.0-Mask_Sc_Layer(i1:i2,j1:j2,Layer_Idx)))/real(Num_Good))
         if (Layer_Idx == 5) Output%Supercooled_Fraction_Layer5(i11:i22,j11:j22) = &
               int(100.0*sum(Sc_Prob_Lay_Pix(i1:i2,j1:j2,Layer_idx)*(1.0-Mask_Sc_Layer(i1:i2,j1:j2,Layer_Idx)))/real(Num_Good))
        endif

         if (Layer_Idx == 4) Output%Supercooled_Total_Fraction(i11:i22,j11:j22) = max(     Output%Supercooled_Total_Fraction(i11,j11), &
                                                                                      real(Output%Supercooled_Fraction_Layer4(i11,j11))/100.)
         if (Layer_Idx == 5) Output%Supercooled_Total_Fraction(i11:i22,j11:j22) = max(     Output%Supercooled_Total_Fraction(i11,j11), &
                                                                                      real(Output%Supercooled_Fraction_Layer5(i11,j11))/100.)
      enddo 

      Output%Supercooled_Total_Fraction(i11:i22,j11:j22) = max(     Output%Supercooled_Total_Fraction(i11,j11), &
                                                               real(Output%Supercooled_Fraction_Layer1(i11,j11))/100. )
      Output%Supercooled_Total_Fraction(i11:i22,j11:j22) = max(     Output%Supercooled_Total_Fraction(i11,j11), &
                                                               real(Output%Supercooled_Fraction_Layer2(i11,j11))/100. )

      Output%Supercooled_Total_Fraction(i11:i22,j11:j22) = max(     Output%Supercooled_Total_Fraction(i11,j11), &
                                                               real(Output%Supercooled_Fraction_Layer3(i11,j11))/100. )

!      Output%Supercooled_Total_Fraction(i11:i22,j11:j22) = max(     Output%Supercooled_Total_Fraction(i11,j11), &
!                                                               real(Output%Supercooled_Fraction_Layer4(i11,j11))/100. )

!      Output%Supercooled_Total_Fraction(i11:i22,j11:j22) = max(     Output%Supercooled_Total_Fraction(i11,j11), &
!                                                               real(Output%Supercooled_Fraction_Layer5(i11,j11))/100. )
!write(unit=6,fmt="(I08,B08,5I5)") Output%Cloud_Layer(i,j), Output%Cloud_Layer(i,j), &
!                              Output%Fraction_Layer1(i,j), &
!                              Output%Fraction_Layer2(i,j), &
!                              Output%Fraction_Layer3(i,j), &
!                              Output%Fraction_Layer4(i,j), &
!!                             Output%Fraction_Layer5(i,j)
    end do element_loop_cover
 end do line_loop_cover

 if (allocated(Pixel_Uncertainty)) deallocate(Pixel_Uncertainty)
 if (allocated(Sc_Prob_Pix)) deallocate(Sc_Prob_Pix)
 if (allocated(Sc_Prob_Base_Pix)) deallocate(Sc_Prob_Base_Pix)
 if (allocated(Mask_Missing_Sc)) deallocate(Mask_Missing_Sc)
 if (allocated(Mask_Sc_Layer)) deallocate(Mask_Sc_Layer)
 if (allocated(Lapse_Rate)) deallocate(Lapse_Rate)
 if (allocated(T_Laytop)) deallocate(T_Laytop)
 if (allocated(Sc_Prob_Lay_Pix)) deallocate(Sc_Prob_Lay_Pix)
 if (allocated(Alt_Freezing)) deallocate(Alt_Freezing)
 if (allocated(Conv_Prob_Pix)) deallocate(Conv_Prob_Pix)
 if (allocated(Mask_Missing_Conv)) deallocate(Mask_Missing_Conv)

 end subroutine COMPUTE_CLOUD_COVER_LAYERS
!----------------------------------------------------------------------
!--- determine cirrus box width
!---
!--- Sensor_Resolution_KM = the nominal resolution in kilometers
!--- Box_Width_KM = the width of the desired box in kilometers
!--- Box_Half_Width = the half width of the box in pixel-space
!----------------------------------------------------------------------
subroutine COMPUTE_BOX_WIDTH(Sensor_Resolution_KM,Box_Width_KM, &
                             Box_Half_Width)

   real, intent(in):: Sensor_Resolution_KM
   integer, intent(in):: Box_Width_KM
   integer, intent(out):: Box_Half_Width

   if (Sensor_Resolution_KM <= 0.0) then
       Box_Half_Width = 20
   else
       Box_Half_Width = int((Box_Width_KM / Sensor_Resolution_KM) / 2)
   endif

end subroutine COMPUTE_BOX_WIDTH

!----------------------------------------------------------------------
! routine to interpolate pressure to flight level altitude.
!----------------------------------------------------------------------
subroutine COMPUTE_ALTITUDE_FROM_PRESSURE_CCL(Line_Idx_Min,Num_Lines,Num_Elems,symbol_yes,Bad_Pixel_Mask,Pc_In,Alt_Out)

!--- Based on Sarah Monette calculations for HS3.  Her calculation is based on:
!--- http://psas.pdx.edu/RocketScience/PressureAltitude_Derived.pdf.

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  integer, intent(in):: Num_Elems
  integer(kind=int1), intent(in) :: symbol_yes
  integer(kind=int1), intent(in), dimension(:,:):: Bad_Pixel_Mask
  real, intent(inout), dimension(:,:) :: Alt_Out
  real, intent(in), dimension(:,:) :: Pc_In
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Line_Start
  integer:: Line_End
  real:: Pc_Temp
  real:: Alt_Temp

  !--- Constants from Sarah Monette
  real, parameter :: PW1 = 227.9 ! hPa
  real, parameter :: PW2 = 56.89 ! hPa
  real, parameter :: PW3 = 11.01 ! hPa
  real, parameter :: P0 = 1013.25 ! hPa
  real, parameter :: LR_OVER_G = 0.190263
  real, parameter :: Z0 = 145422.16 ! feet
  real, parameter :: LN_1 = -20859.0
  real, parameter :: LN_2 = 149255.0
  real, parameter :: PN_4 = 0.000470034
  real, parameter :: PN_3 = -0.364267
  real, parameter :: PN_2 = 47.5627
  real, parameter :: PN_1 = -2647.45
  real, parameter :: PN_0 = 123842.0

  !--- save elements
  Line_Start = Line_Idx_Min
  Line_End = Line_Start + Num_Lines - 1

  !--- initialize
  Alt_Out = Missing_Value_Real4

  !----------------------------------------------------------
  ! loop through segment
  !----------------------------------------------------------
  Line_Loop_1: do Line_Idx = Line_Start, Line_End
  Element_Loop_1: do Elem_Idx = 1, Num_Elems

     !--- Initialize temporary value each time.
     Pc_Temp = Missing_Value_Real4
     Alt_Temp = Missing_Value_Real4

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == symbol_yes) cycle

     !--- check if cld-top pressure is valid
     if (Pc_In(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- Place valid pressure in temp variable for readability in the
     !--- calculations below.
     Pc_Temp = Pc_In(Elem_Idx,Line_Idx)

     !--- calculated altitude, in feet, from pressure.
     !--- 1st pivot point is directly from the pressure to
     !--- altitude from above reference.
     if (Pc_Temp > PW1) then
       Alt_Temp = (1.0 - (Pc_Temp/P0)**LR_OVER_G) * Z0
     endif

     !--- 2nd pivot point was modeled best with a natural log
     !--- fit.  From Sarah Monette.
     if (Pc_Temp <= PW1 .AND. Pc_Temp >= PW2) then
       Alt_Temp = LN_1 * LOG(Pc_Temp) + LN_2
     endif

     !--- 3rd pivot point. Modeled best with a polynomial
     !--- fit from Sarah Monette.
     if (Pc_Temp < PW2 .AND. Pc_Temp >= PW3) then
       Alt_Temp = (PN_4*Pc_Temp**4) + (PN_3*Pc_Temp**3) + (PN_2*Pc_Temp**2) + &
                  (PN_1*Pc_Temp) + PN_0
     endif

     if (Pc_Temp < PW3 .or. Alt_Temp < 0) then
       Alt_Temp = Missing_Value_Real4
     endif


     !--- Assign final altitude, in feet, to the level2 array.
     Alt_Out(Elem_Idx,Line_Idx) = Alt_Temp

  enddo Element_Loop_1
  enddo Line_Loop_1

end subroutine COMPUTE_ALTITUDE_FROM_PRESSURE_CCL

! Supercooled Cloud Probability
!--------------------------------------------------------------------------------------------------
subroutine SUPERCOOLED_CLOUD_PROBABILITY_local(Cloud_Temperature,Supercooled_Cld_Prob)
  real(kind=real4), intent(in):: Cloud_Temperature
  real(kind=real4), intent(out):: Supercooled_Cld_Prob

  integer:: Tc_Idx
  !--- initialize to Missing
  Supercooled_Cld_Prob = Missing_Value_Real4

  !--- filter missing temps
  if (Cloud_Temperature /= Missing_Value_Real4) then

       Tc_Idx = minloc(abs(Cloud_Temperature - Sc_Tc_Lut),1)
       Tc_Idx = max(1,min(N_Sc_Lut,Tc_Idx))

       Supercooled_Cld_Prob = Sc_Prob_Lut(Tc_Idx)
  endif

end subroutine SUPERCOOLED_CLOUD_PROBABILITY_local

! InterpoLate within profiles knowing P to determine T and Z
!-----------------------------------------------------------------
subroutine KNOWING_P_COMPUTE_T_Z(CCL_RTM_NWP,P,T,Z,Lev_Idx)

     type(ccl_rtm_nwp_struct), intent(in) :: CCL_RTM_NWP
     real, intent(in):: P
     real, intent(out):: T
     real, intent(out):: Z
     integer, intent(out):: Lev_Idx
     real:: dp
     real:: dt
     real:: dz

     !--- initialize
     T = MISSING_VALUE_REAL4
     Z = MISSING_VALUE_REAL4
     Lev_Idx = MISSING_VALUE_integer4

     !--- check for missing
     if (P == MISSING_VALUE_REAL4) return

     !--- interpoLate pressure profile
     call LOCATE(CCL_RTM_NWP%P_Prof,Num_Levels_RTM_Prof,P,Lev_Idx)
     Lev_Idx = max(1,min(Num_Levels_RTM_Prof-1,Lev_Idx))

     dp = CCL_RTM_NWP%P_Prof(Lev_Idx+1) - CCL_RTM_NWP%P_Prof(Lev_Idx)
     dt = CCL_RTM_NWP%T_Prof(Lev_Idx+1) - CCL_RTM_NWP%T_Prof(Lev_Idx)
     dz = CCL_RTM_NWP%Z_Prof(Lev_Idx+1) - CCL_RTM_NWP%Z_Prof(Lev_Idx)

     !--- perform interpoLation
       if (dp /= 0.0) then
           T = CCL_RTM_NWP%T_Prof(Lev_Idx) + dt/dp * (P - CCL_RTM_NWP%P_Prof(Lev_Idx))
           Z = CCL_RTM_NWP%Z_Prof(Lev_Idx) + dz/dp * (P - CCL_RTM_NWP%P_Prof(Lev_Idx))
       else
           T = CCL_RTM_NWP%T_Prof(Lev_Idx)
           Z = CCL_RTM_NWP%Z_Prof(Lev_Idx)
       endif

       !--- Some negative cloud heights are observed because  of bad height
       !--- NWP profiles.
       if (Z < 0) then
         Z = ZC_FLOOR
       endif

end subroutine KNOWING_P_COMPUTE_T_Z
!------------------------------------------

!------------------------------------------------------------------------------
! Null Pixel Level Pointers 
!------------------------------------------------------------------------------
subroutine NULL_PIX_POINTERS(Input, RTM_NWP)

   type(ccl_input_struct), intent(inout) :: Input
   type(ccl_rtm_nwp_struct), intent(inout) :: RTM_NWP

   RTM_NWP%T_Prof => null()
   RTM_NWP%T_Prof_1 => null()
   RTM_NWP%T_Prof_2 => null()
   RTM_NWP%T_Prof_3 => null()

   RTM_NWP%Z_Prof => null()
   RTM_NWP%Z_Prof_1 => null()
   RTM_NWP%Z_Prof_2 => null()
   RTM_NWP%Z_Prof_3 => null()

end subroutine NULL_PIX_POINTERS

!-----------------------------------------------------------
! end of MODULE
!-----------------------------------------------------------
end module CCL_MODULE
