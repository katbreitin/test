!$Id: acha_module.f90 4092 2021-03-02 22:05:14Z heidinger $
module ACHA_PRIOR_MODULE
!---------------------------------------------------------------------
! ACHA Prior
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
!----------------------------------------------------------------------
  use ACHA_SERVICES_MOD, only : &
           real4, int1, int4, real8, dtor, &
           Acha_output_struct,ACHA_SYMBOL_STRUCT, &
           Acha_input_struct, Acha_rtm_nwp_struct, &
           ACHA_FETCH_PIXEL_NWP_RTM, &
           LOCATE, Acha_Diag_Struct, Acha_Dump_Struct, COUNTSUBSTRING, &
           ABI_Use_104um_Flag
! use TIMER_MOD
  use CX_REAL_BOOLEAN_MOD
  use KDTREE2_MODULE
  use ACHA_NUM_MOD

  implicit none

  public:: ACHA_PRIOR
  private:: COMPUTE_CIRRUS_APRIORI
  private:: SELECT_CHANNELS_EYRE_MENZEL
  private:: EYRE_MENZEL
  private:: SELECT_CHAN_RAD
  private:: EMISSIVITY
! private:: NULL_PIX_POINTERS

  !--- include the non-system specific variables
  include 'include/acha_parameters.inc'

  real, private, PARAMETER:: MISSING_VALUE_REAL4 = -999.0
  integer(kind=int1), private, PARAMETER:: MISSING_VALUE_integer1 = -128_int1
  integer(kind=int4), private, PARAMETER:: MISSING_VALUE_integer4 = -999
  type(ACHA_SYMBOL_STRUCT), private :: Symbol

  integer, private, parameter:: N_Chan_EM_Max = 12
  integer, dimension(N_Chan_EM_Max), private, save:: EM_Chan_Idx
  integer, private, save :: N_Chan_EM
  real, dimension(Num_Levels_Rtm_Prof), private:: BB_Rad
  real, private:: Clr_Rad
  real, private:: Obs_Rad
  integer, private:: Chan_Idx_11um
! type(timer) :: ts

  type data_ch_em
    logical :: is_set = .false.
    real,dimension(Num_Levels_Rtm_Prof):: BB_Rad
    real:: Clr_Rad
    real:: Obs_Rad
  end type  data_ch_em
  type data_em
    type(data_ch_em) :: ch(N_Chan_EM_Max)
    logical :: is_set = .false.
  contains
    procedure :: set => data_em_set
  end type data_em
  
  contains

!------------------------------------------------------------------------------
! AWG Cloud Height Algorithm (ACHA)
!
! Author: Andrew Heidinger, NOAA
!
!----------------------------------------------------------------------
!
!------------------------------------------------------------------------------
  subroutine  ACHA_PRIOR(Input, Symbol, Output, Diag)

  !============================================================================
  !  Argument Declaration
  !============================================================================

  type(acha_input_struct), intent(inout) :: Input
  type(acha_symbol_struct), intent(in) :: Symbol
  type(acha_output_struct), intent(inout) :: Output
  type(acha_diag_struct), intent(inout), optional :: Diag

  !============================================================================
  !  Pixel level RTM structure
  !============================================================================

  type(acha_rtm_nwp_struct) :: ACHA_RTM_NWP

  !============================================================================
  !  Local Variable Declaration
  !============================================================================

  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Param_Idx
  integer:: Lev_Idx
  integer:: Ivza
  integer:: ilrc
  integer:: jlrc
  integer:: ierror
  integer:: Elem_Idx_LRC
  integer:: Line_Idx_LRC
  real, dimension(:,:), allocatable:: Ice_Prob_Norm
  real, dimension(:,:), allocatable:: Pc_EM,Tc_EM,Zc_EM,Ec_EM,N_EM,Res_EM,N_Std_EM,CV_EM,Ec_Res_EM
  integer (kind=int1), dimension(:,:), allocatable:: Mask_In,Mask_Out
  real, dimension(:,:), allocatable:: Temp_Array_R4
  integer (kind=int1), dimension(:,:), allocatable:: Temp_Array_I1

  integer, parameter:: N_Idx_Found = 6 ! the number of surrounding indices to average
  integer, parameter:: Kdtree_Train_Count_Thresh = 10 ! the values needs to be larger than number of predictors

  logical:: USE_QRNN_FLAG
  logical:: USE_EM_FLAG
  logical:: USE_LRC_LOCAL_FLAG
  logical:: USE_CIRAP_FLAG
  logical:: ALLOW_RETYPE
  logical:: APPLY_MEDIAN
  logical:: OPICEKD_Flag

  integer, save:: Diag_Warning_Flag = 0

  integer:: klev,ierr,Level_Within_Inversion_Flag, RTM_NWP_Error_Flag

  real:: Dummy_R4, dT_dP

! logical :: first_run = .true.
! if (first_run) then
!    call ts % init(['sel chn rad.','ALL ACHA   .'])
!    first_run = .false.
! end if
! call ts % tic(2)

  !--- initialize diagnostic output
  if (present(Diag) .and. Diag_Warning_Flag == 0) then
      print *, "ACHA PRIOR ===>  Diagnostic Output Turned On"
      Diag_Warning_Flag = 1
  endif
  if (present(Diag)) Diag%Array_1 = MISSING_VALUE_REAL4
  if (present(Diag)) Diag%Array_2 = MISSING_VALUE_REAL4
  if (present(Diag)) Diag%Array_3 = MISSING_VALUE_REAL4

  !--- set flags
  USE_QRNN_FLAG = .false.
  USE_EM_FLAG = .true.
  USE_CIRAP_FLAG = .true.
  USE_LRC_LOCAL_FLAG = .false. !.true.
  ALLOW_RETYPE = .true.
  APPLY_MEDIAN = .true.
  OPICEKD_Flag = .true.

  !--- initialize
  Input%Tc_Ap = MISSING_VALUE_REAL4
  Input%Ec_Ap = MISSING_VALUE_REAL4
  Input%Beta_Ap = MISSING_VALUE_REAL4
  Input%Ice_Prob_Ap = MISSING_VALUE_REAL4
  Input%Lower_Tc_Ap = MISSING_VALUE_REAL4

  Input%Tc_Ap_Uncer = MISSING_VALUE_REAL4
  Input%Ec_Ap_Uncer = MISSING_VALUE_REAL4
  Input%Beta_Ap_Uncer = MISSING_VALUE_REAL4
  Input%Ice_Prob_Ap_Uncer = MISSING_VALUE_REAL4
  Input%Lower_Tc_Ap_Uncer = MISSING_VALUE_REAL4

  !--------------------------------------------------------------------------------
  ! Default with Opaque
  !--------------------------------------------------------------------------------
  where (Input%Invalid_Data_Mask == symbol%NO)
    Input%Tc_Ap = Input%Tc_Opaque
!   Input%Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Opaque_Ice_Type/Input%Cosine_Zenith_Angle)
    Input%Ec_Ap = 0.95
    Input%Beta_Ap = 1.1
    Input%Ice_Prob_Ap = 0.5
    Input%Lower_Tc_Ap = Input%Surface_Temperature

    !-------------------------------------------------------------------------------------------------
    !  Default Uncertainty
    !-------------------------------------------------------------------------------------------------
    Input%Tc_Ap_Uncer = 40.0 ! 20.0 K
    Input%Ec_Ap_Uncer = 1.0
    Input%Beta_Ap_Uncer = 0.5
    Input%Ice_Prob_Ap_Uncer = 1.0
    Input%Lower_Tc_Ap_Uncer = 40.0
  endwhere

  !--------------------------------------------------------------------------------
  ! Use Cirrus Apriori  (based on latitude and tropopause temperature)
  !--------------------------------------------------------------------------------
  if (USE_CIRAP_FLAG) then
    Line_Loop_1: do Line_Idx = 1, Input%Number_Of_Lines
      Elem_Loop_1: do Elem_Idx = 1, Input%Number_Of_Elements
        if (Input%Cloud_Type(Elem_Idx,Line_Idx) == Symbol%CIRRUS_TYPE .or. &
            Input%Cloud_Type(Elem_Idx,Line_Idx) == Symbol%OVERLAP_TYPE) then
            call COMPUTE_CIRRUS_APRIORI(Input%Tropopause_Temperature(Elem_Idx,Line_Idx), &
                                        Input%Latitude(Elem_Idx,Line_Idx), &
                                        Input%Tc_Ap(Elem_Idx,Line_Idx), &
                                        Input%Tc_Ap_Uncer(Elem_Idx,Line_Idx))
            Input%Ec_Ap(Elem_Idx,Line_Idx) = 1.0 - exp(-1.0*Tau_Ap_Cirrus_Type / &
                                             Input%Cosine_Zenith_Angle(Elem_Idx,Line_Idx))
            Input%Tc_Ap_Uncer(Elem_Idx,Line_Idx) = 40.0  !K
        endif
      enddo Elem_Loop_1
    enddo Line_Loop_1
  endif

  !-------------------------------------------------------------------------------------------------
  ! Quantile Regression Neural Network
  !-------------------------------------------------------------------------------------------------
  if (USE_QRNN_FLAG) then

    Line_Loop_3: do Line_Idx = 1, Input%Number_Of_Lines
      Elem_Loop_3: do Elem_Idx = 1, Input%Number_Of_Elements

        if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) cycle
        if (Input%QRNN_Pc(Elem_Idx,Line_Idx) < 100.0) cycle
        if (Input%Cloud_Mask(Elem_Idx,Line_Idx) <= 1) cycle
        if (Input%QRNN_Pc(Elem_Idx,Line_Idx) < 150.0) cycle

        call ACHA_FETCH_PIXEL_NWP_RTM(Input, Symbol, &
                                  Elem_Idx,Line_Idx, &
                                  ACHA_RTM_NWP,RTM_NWP_Error_Flag)
        if (RTM_NWP_Error_Flag /= 0) cycle

        !--- compute temperature and height
        call KNOWING_P_COMPUTE_T_Z(ACHA_RTM_NWP,Input%QRNN_Pc(Elem_Idx,Line_Idx), &
                                   Input%Tc_Ap(Elem_Idx,Line_Idx), Dummy_R4, Lev_Idx)

        !--- estimate lapse rate for temperature uncertainty
        dT_dP = (ACHA_RTM_NWP%T_Prof(Lev_Idx) - ACHA_RTM_NWP%T_Prof(Lev_Idx-1)) / &
                (ACHA_RTM_NWP%P_Prof(Lev_Idx) - ACHA_RTM_NWP%P_Prof(Lev_Idx-1))

        Input%Tc_Ap_Uncer(Elem_Idx,Line_Idx) = Input%QRNN_Pc_Uncer(Elem_Idx,Line_Idx) * dT_dP

        !-- compute emissivity
        BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_110um
        Clr_Rad = Input%Rad_Clear_110um(Elem_Idx,Line_Idx)
        Obs_Rad = Input%Rad_110um(Elem_Idx,Line_Idx)
        call SELECT_CHAN_RAD(Elem_Idx,Line_Idx,31,Input,ACHA_RTM_NWP)
        Input%Ec_Ap(Elem_Idx,Line_Idx) = EMISSIVITY(Obs_Rad, Clr_Rad, BB_Rad(Lev_Idx))
        Input%Ec_Ap_Uncer(Elem_Idx,Line_Idx) = 0.5

      end do Elem_Loop_3
   end do Line_Loop_3

  endif
  !-------------------------------------------------------------------------------------------------
  ! Eyre Menzel
  !-------------------------------------------------------------------------------------------------
  !--- if no absorbing channel, turn off EM
  if (Input%Chan_On_062um == symbol%NO .and. &
       Input%Chan_On_067um == symbol%NO .and. &
       Input%Chan_On_073um == symbol%NO .and. &
       Input%Chan_On_097um == symbol%NO .and. &
       Input%Chan_On_133um == symbol%NO .and. &
       Input%Chan_On_136um == symbol%NO .and. &
       Input%Chan_On_139um == symbol%NO .and. &
       Input%Chan_On_142um == symbol%NO) then

       USE_EM_FLAG = .false.

  endif

  if (USE_EM_FLAG) then
    allocate(Pc_EM(Input%Number_Of_Elements,Input%Number_Of_Lines))
    allocate(Tc_EM(Input%Number_Of_Elements,Input%Number_Of_Lines))
    allocate(Zc_EM(Input%Number_Of_Elements,Input%Number_Of_Lines))
    allocate(Ec_EM(Input%Number_Of_Elements,Input%Number_Of_Lines))
    allocate(N_EM(Input%Number_Of_Elements,Input%Number_Of_Lines))
    allocate(Res_EM(Input%Number_Of_Elements,Input%Number_Of_Lines))
    allocate(N_Std_EM(Input%Number_Of_Elements,Input%Number_Of_Lines))
    allocate(CV_EM(Input%Number_Of_Elements,Input%Number_Of_Lines))
    allocate(Ec_Res_EM(Input%Number_Of_Elements,Input%Number_Of_Lines))

    Pc_EM = MISSING_VALUE_REAL4
    Tc_EM = MISSING_VALUE_REAL4
    Zc_EM = MISSING_VALUE_REAL4
    Ec_EM = MISSING_VALUE_REAL4
    N_EM = MISSING_VALUE_REAL4
    Res_EM = MISSING_VALUE_REAL4
    N_Std_EM = MISSING_VALUE_REAL4
    CV_EM = MISSING_VALUE_REAL4
    Ec_Res_EM = MISSING_VALUE_REAL4

    call SELECT_CHANNELS_EYRE_MENZEL(Input,Symbol)

    if (N_Chan_EM > 0) then
       call EYRE_MENZEL(Input,Symbol,Pc_EM,Tc_EM,Zc_EM,Ec_EM,N_EM,Res_EM,N_Std_EM,CV_EM,Ec_Res_EM)
    endif

    Diag%Array_1 = Tc_EM
    Diag%Array_2 = Pc_EM
    Diag%Array_3 = Ec_EM

    !--- use EM for ice clouds
    where(Pc_EM /= MISSING_VALUE_REAL4 .and. &
          (Input%Cloud_Type == Symbol%CIRRUS_TYPE .or.  &
           Input%Cloud_Type == Symbol%OVERLAP_TYPE .or. &
           Input%Cloud_Type == Symbol%OPAQUE_ICE_TYPE .or. &
           Input%Cloud_Type == Symbol%OVERSHOOTING_TYPE) .and. &
           Ec_EM > 0.01 .and. CV_EM < 1.0)

           Input%Tc_Ap = Tc_EM
           Input%Ec_Ap = Ec_EM
           Input%Tc_Ap_Uncer = 20.0
           Input%Ec_Ap_Uncer = 0.5

    end where

    !--- believe high clouds where phase is uncertain
    where(Pc_EM /= MISSING_VALUE_REAL4 .and. &
          Tc_EM < Input%Tc_Opaque .and. &
          Input%Cloud_Phase_Uncertainty > 0.10 .and. &
          !Pc_EM < 700.0 .and. &
          Tc_EM < 250.0 .and. &
          Ec_EM > 0.01 .and. &
          CV_EM < 1.0)
          Input%Tc_Ap = Tc_EM
          Input%Ec_Ap = Ec_EM
          Input%Tc_Ap_Uncer = 20.0
          Input%Ec_Ap_Uncer = 0.5
    endwhere



!   print *, 'number using EM = ', count( Input%Tc_Ap ==Tc_EM)

!   !---- use highly confident EM results for high clouds not matter the type
!   where(Pc_EM /= MISSING_VALUE_REAL4 .and. &
!         Tc_EM < Input%Tc_Opaque .and. &
!         Pc_EM < 500.0 .and. &
!         Ec_EM > 0.01 .and. &
!         CV_EM < 1.0 .and. &
!         Ec_Res_EM < 0.10 .and. &
!         (Input%Ice_Cloud_Probability > 0.10 .or. &
!         Input%Cloud_Phase_Uncertainty > 0.10))
!
!         Input%Tc_Ap = Tc_EM
!         Input%Ec_Ap = Ec_EM
!         Input%Tc_Ap_Uncer = 20.0
!         Input%Ec_Ap_Uncer = 0.5

!   end where
! do Line_Idx = 1, Input%Number_Of_Lines
! do Elem_Idx = 1, Input%Number_Of_Elements

!    !------ find water pixel that are replaced by EM
!    if (Input%Cloud_Type(Elem_Idx,Line_Idx) >= 2 .and. Input%Cloud_Type(Elem_Idx,Line_Idx) <= 4 .and. &
!        (Input%Tc_Ap(Elem_Idx,Line_Idx) - Input%Tc_Opaque(Elem_Idx,Line_Idx)  < -5.0) .and. &
!         Pc_EM(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) then

!        print *, 'WATER EM ', Input%Tc_Ap(Elem_Idx,Line_Idx), &
!        Tc_EM(Elem_Idx,Line_Idx),Input%Tc_Opaque(Elem_Idx,Line_Idx),  &
!        Pc_Em(Elem_Idx,Line_Idx),Ec_EM(Elem_Idx,Line_Idx),  &
!        Res_EM(Elem_Idx,Line_Idx),CV_Em(Elem_Idx,Line_Idx),Ec_Res_Em(Elem_Idx,Line_Idx), &
!        Input%Ice_Cloud_Probability(Elem_Idx,Line_Idx), &
!        Input%Cloud_Phase_Uncertainty(Elem_Idx,Line_Idx)

!    endif
! enddo
! enddo

!   where(Input%Cloud_Phase_Uncertainty > 0.10 .and. Tc_EM /=MISSING_VALUE_REAL4 .and. Tc_EM < 250.0)
!          Input%Tc_Ap = Tc_EM
!          Input%Ec_Ap = Ec_EM
!          Input%Tc_Ap_Uncer = 20.0
!          Input%Ec_Ap_Uncer = 0.5
!  endwhere

!  !---- don't allow EM for highly confident water clouds
!   where(Input%Ice_Cloud_Probability < 0.05 .and. &
!         Input%Cloud_Phase_Uncertainty < 0.10)

!  endwhere

   deallocate(Pc_EM,Tc_EM,Zc_EM,Ec_EM,N_EM,Res_EM,N_Std_EM,CV_EM,Ec_Res_EM)

  endif
  !-------------------------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------------------------
  if ((USE_LRC_LOCAL_FLAG) .and. &
      (associated(Input%Elem_Idx_LRC_Input)) .and. &
      (associated(Input%Elem_Idx_LRC_Input))) then

     Line_Loop_2: do Line_Idx = 1, Input%Number_Of_Lines
       Elem_Loop_2: do Elem_Idx = 1, Input%Number_Of_Elements

       Elem_Idx_LRC = Input%Elem_Idx_LRC_Input(Elem_Idx,Line_Idx)
       Line_Idx_LRC = Input%Line_Idx_LRC_Input(Elem_Idx,Line_Idx)

       if (Elem_Idx_LRC <= 0 .or. Line_Idx_LRC <= 0) cycle

       end do Elem_Loop_2
     end do Line_Loop_2

  endif

  !-------------------------------------------------------------------------------------------------
  ! KDTREE Lower Cloud
  !-------------------------------------------------------------------------------------------------
  allocate(Mask_In(Input%Number_Of_Elements,Input%Number_Of_Lines), &
           Mask_Out(Input%Number_Of_Elements,Input%Number_Of_Lines))
  allocate(Temp_Array_R4(Input%Number_Of_Elements,Input%Number_Of_Lines))

  Mask_Out = 0
  where(Input%Cloud_Type == Symbol%CIRRUS_TYPE .or. &
        Input%Cloud_Type == Symbol%OPAQUE_ICE_TYPE .or.  &
        Input%Cloud_Type == Symbol%OVERLAP_TYPE)
        Mask_Out = 1_int1
  end where

  Mask_In = 0
  where( (Input%Cloud_Type == Symbol%FOG_TYPE .or. &
          Input%Cloud_Type == Symbol%WATER_TYPE) .and.  &
          Input%Tc_Ap > 250 .and.  &
          Input%Tc_Ap /= MISSING_VALUE_REAL4)
          Mask_In = 1_int1
  end where

  if (COUNT(Mask_In == 1) > Kdtree_Train_Count_Thresh)  then
     Temp_Array_R4 = MISSING_VALUE_REAL4
     call KD_TREE_INTERP_2pred(Mask_In, Mask_Out, Input%Latitude, Input%Longitude, &
                               Input%Number_Of_Elements, Input%Number_Of_Lines, &
                               N_Idx_Found, Input%Tc_Ap, Temp_Array_R4)

     where(Mask_Out == 1 .and. Temp_Array_R4 /= MISSING_VALUE_REAL4)
          Input%Lower_Tc_Ap = Temp_Array_R4
     endwhere

  endif

  deallocate(Mask_In, Mask_Out)
  deallocate(Temp_Array_R4)

  !-------------------------------------------------------------------------------------------------
  ! KDTREE Ice Cloud
  !-------------------------------------------------------------------------------------------------
  if (OPICEKD_Flag) then

     allocate(Mask_In(Input%Number_Of_Elements,Input%Number_Of_Lines), &
              Mask_Out(Input%Number_Of_Elements,Input%Number_Of_Lines))
     allocate(Temp_Array_R4(Input%Number_Of_Elements,Input%Number_Of_Lines))

     Mask_Out = 0
     where(Input%Cloud_Type == Symbol%CIRRUS_TYPE .or. &
           Input%Cloud_Type == Symbol%OVERLAP_TYPE)
           Mask_Out = 1_int1
     end where

     Mask_In = 0
     where( (Input%Cloud_Type == Symbol%OPAQUE_ICE_TYPE .or. &
             Input%Cloud_Type == Symbol%OVERSHOOTING_TYPE) .and.  &
             Input%Tc_Ap < 250 .and.  &
             Input%Tc_Ap /= MISSING_VALUE_REAL4)
             Mask_In = 1_int1
     end where

     if (COUNT(Mask_In == 1) > Kdtree_Train_Count_Thresh)  then
        Temp_Array_R4 = MISSING_VALUE_REAL4
        call KD_TREE_INTERP_2pred(Mask_In, Mask_Out, Input%Latitude, Input%Longitude, &
                               Input%Number_Of_Elements, Input%Number_Of_Lines, &
                               N_Idx_Found, Input%Tc_Ap, Temp_Array_R4)
        where(Mask_Out == 1 .and. Temp_Array_R4 /= MISSING_VALUE_REAL4)
            Input%Tc_Ap = Temp_Array_R4
        endwhere
     endif

     deallocate(Mask_In, Mask_Out)
     deallocate(Temp_Array_R4)

  endif

  !-------------------------------------------------------------------------------------------------
  ! Ice Cloud Apriori
  !-------------------------------------------------------------------------------------------------
  allocate(Ice_Prob_Norm(Input%Number_Of_Elements,Input%Number_Of_Lines))
  Ice_Prob_Norm = MISSING_VALUE_REAL4
  where(Input%Cloud_Mask /= Symbol%CLEAR .and. Input%Cloud_Mask /= Symbol%PROB_CLEAR .and. &
        Input%Ice_Cloud_Probability /= MISSING_VALUE_REAL4 .and.  &
        Input%Cloud_Probability > 0.0)
     Ice_Prob_Norm = Input%Ice_Cloud_Probability / Input%Cloud_Probability
  end where

  if (USE_TYPE_FOR_ICE_PROB) then

     Input%Ice_Prob_Ap = 0.0
     Input%Ice_Prob_Ap_Uncer = 0.5

     where(Input%Cloud_Type == Symbol%Cirrus_Type .or. &
           Input%Cloud_Type == Symbol%Overlap_Type .or. &
           Input%Cloud_Type == Symbol%Opaque_Ice_Type .or. &
           Input%Cloud_Type == Symbol%Overshooting_Type)

           Input%Ice_Prob_Ap = 1.0
           Input%Ice_Prob_Ap_Uncer = 0.5

     end where

  else

     where(Ice_Prob_Norm /= MISSING_VALUE_REAL4)
        Input%Ice_Prob_Ap = Ice_Prob_Norm
        Input%Ice_Prob_Ap_Uncer = Input%Cloud_Phase_Uncertainty
     end where

  endif

  deallocate(Ice_Prob_Norm)

  !-----------------------------------------------------------------------
  !--- Apply Median Filter to Tc_Ap
  !-----------------------------------------------------------------------
  if (APPLY_MEDIAN) then

     allocate(Temp_Array_R4(Input%Number_Of_Elements,Input%Number_Of_Lines))
     allocate(Temp_Array_I1(Input%Number_Of_Elements,Input%Number_Of_Lines))

     Temp_Array_I1 = 0
     where(Input%Tc_Ap == MISSING_VALUE_REAL4)
        Temp_Array_I1 = 1
     endwhere

     call COMPUTE_MEDIAN_SEGMENT(Input%Tc_Ap,Temp_Array_I1,1, &
                   Input%Number_of_Elements,Input%Number_of_Lines, &
                   Temp_Array_R4)
     Input%Tc_Ap = Temp_Array_R4

     deallocate(Temp_Array_R4)
     deallocate(Temp_Array_I1)

   endif

  !------------------------------------------------------------------------------------------------
  ! ensure Ice_Prob_Ap and Tc_Ap are consistent
  !------------------------------------------------------------------------------------------------
  if (ALLOW_RETYPE) then
    where(Input%Tc_Ap /= MISSING_VALUE_REAL4 .and. Input%Tc_Ap < 240.0 .and. Input%Ice_Prob_Ap < 0.5)
          Input%Ice_Prob_Ap = 0.75
          Input%Ice_Prob_Ap_Uncer = 1.0
          Input%Tc_Ap_Uncer = 100.0
          Input%Ec_Ap_Uncer = 1.0
    endwhere
  endif

  !------------------------------------------------------------------------------------------------
  ! apply binary flag = force ice prob ap to be 0 or 1
  !------------------------------------------------------------------------------------------------
  if (BINARY_ICE_PROB) then
     where(Input%Ice_Prob_Ap /= MISSING_VALUE_REAL4 .and. &
           Input%Ice_Prob_Ap < 0.5)
           Input%Ice_Prob_Ap = 0.0
     end where
     where(Input%Ice_Prob_Ap >= 0.5)
           Input%Ice_Prob_Ap = 1.0
     end where
  endif

  if (CONSTRAIN_ICE_PROB) then
     Input%Ice_Prob_Ap_Uncer = Ice_Prob_Ap_Uncer_Min
  endif

  !-------------------------------------------------------------------------------------------------
  ! Beta Apriori
  !-------------------------------------------------------------------------------------------------
  Input%Beta_Ap = Beta_Ap_Water
  Input%Beta_Ap_Uncer = Beta_Ap_Uncer_Water

  where(Input%Ice_Prob_Ap >= 0.5)

         Input%Beta_Ap = Beta_Ap_Ice
         Input%Beta_Ap_Uncer = Beta_Ap_Uncer_Ice

  end where

  !-------------------------------------------------------------------------------------------------
  ! Increase Uncertainty for Tc_Ap if phase uncertain
  !-------------------------------------------------------------------------------------------------
  where(Input%Cloud_Phase_Uncertainty >= 0.20)
    Input%Tc_Ap_Uncer = 100.0
    Input%Ec_Ap_Uncer = 1.0
  endwhere

  !-----------------------------------------------------------------------
  !--- Tighten up Ec_Ap for confident water clouds
  !-----------------------------------------------------------------------
  where(Input%Ice_Prob_Ap < 0.5 .and. Input%Ice_Prob_Ap /= MISSING_VALUE_REAL4)
    !Input%Ec_Ap_Uncer = Input%Ice_Prob_Ap + 0.1
    Input%Ec_Ap_Uncer = 0.2
  endwhere

  !-------------------------------------------------------------------------------------------------
  ! Set Clear to Missing
  !-------------------------------------------------------------------------------------------------
  where(Input%Cloud_Mask == Symbol%CLEAR .or. Input%Cloud_Mask == Symbol%PROB_CLEAR)
   Input%Tc_Ap = MISSING_VALUE_REAL4
   Input%Ec_Ap = MISSING_VALUE_REAL4
   Input%Beta_Ap = MISSING_VALUE_REAL4
   Input%Ice_Prob_Ap = MISSING_VALUE_REAL4
   Input%Lower_Tc_Ap = MISSING_VALUE_REAL4
   Input%Tc_Ap_Uncer = MISSING_VALUE_REAL4
   Input%Ec_Ap_Uncer = MISSING_VALUE_REAL4
   Input%Beta_Ap_Uncer = MISSING_VALUE_REAL4
   Input%Ice_Prob_Ap_Uncer = MISSING_VALUE_REAL4
   Input%Lower_Tc_Ap_Uncer = MISSING_VALUE_REAL4
  end where

  !-------------------------------------------------------------------------------------------------
  !  Constrain
  !-------------------------------------------------------------------------------------------------
  where(Input%Tc_Ap_Uncer /= MISSING_VALUE_REAL4 .and.  &
        Input%Tc_Ap_Uncer < Tc_Ap_Uncer_Min)
        Input%Tc_Ap_Uncer = Tc_Ap_Uncer_Min
  endwhere

  where(Input%Ec_Ap_Uncer /= MISSING_VALUE_REAL4 .and.  &
        Input%Ec_Ap_Uncer < Ec_Ap_Uncer_Min)
        Input%Ec_Ap_Uncer = Ec_Ap_Uncer_Min
  endwhere

  where(Input%Beta_Ap_Uncer /= MISSING_VALUE_REAL4 .and.  &
        Input%Beta_Ap_Uncer < Beta_Ap_Uncer_Min)
        Input%Beta_Ap_Uncer = Beta_Ap_Uncer_Min
  endwhere

  where(Input%Lower_Tc_Ap_Uncer /= MISSING_VALUE_REAL4 .and.  &
        Input%Lower_Tc_Ap_Uncer < Lower_Tc_Ap_Uncer_Min)
        Input%Lower_Tc_Ap_Uncer = Lower_Tc_Ap_Uncer_Min
  endwhere

  where(Input%Ice_Prob_Ap_Uncer /= MISSING_VALUE_REAL4 .and.  &
        Input%Ice_Prob_Ap_Uncer < Ice_Prob_Ap_Uncer_Min)
        Input%Ice_Prob_Ap_Uncer = Ice_Prob_Ap_Uncer_Min
  endwhere

!------------------------------------------------------------------------------------
! extra code - to make pc and zc ap for diagnosis
!------------------------------------------------------------------------------------
!  Line_Loop_4: do Line_Idx = 1, Input%Number_Of_Lines
!     Elem_Loop_4: do Elem_Idx = 1, Input%Number_Of_Elements
!
!       if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) cycle
!
!       call ACHA_FETCH_PIXEL_NWP_RTM(Input, Symbol, &
!                                 Elem_Idx,Line_Idx, &
!                                 ACHA_RTM_NWP,RTM_NWP_Error_Flag)
!       if (RTM_NWP_Error_Flag /= 0) cycle
!
!       call KNOWING_T_COMPUTE_P_Z(ACHA_RTM_NWP,Symbol, Input%Cloud_Type(Elem_Idx,Line_Idx), &
!           Diag%Array_3(Elem_Idx,Line_Idx), &
!           Input%Tc_Ap(Elem_Idx,Line_Idx), &
!           Dummy_R4, &
!           Input%Tropopause_Temperature(Elem_Idx,Line_Idx), &
!           Input%Tropopause_Height(Elem_Idx,Line_Idx), &
!           Input%Tropopause_Pressure(Elem_Idx,Line_Idx), &
!           klev,ierr,Level_Within_Inversion_Flag)
!
!  end do Elem_Loop_4
!  end do Line_Loop_4
!-------------------------------------------------------------------------------------------
! end extra source
!-------------------------------------------------------------------------------------------
! call ts % tac(2)
! call ts % summary()

  end subroutine ACHA_PRIOR

!----------------------------------------------------------------------------
! estimate cirrus aprior temperature and uncertainty from a precomputed
! latitude table (stored in acha_parameters.inc)
!----------------------------------------------------------------------------
subroutine COMPUTE_CIRRUS_APRIORI(t_tropo, latitude, tc_apriori, tc_apriori_uncer)
  real, intent(in):: t_tropo
  real, intent(in):: latitude
  real, intent(out):: tc_apriori
  real, intent(out):: tc_apriori_uncer

  integer:: lat_idx
  real, parameter:: lat_min = -90.0
  real, parameter:: delta_lat = -10.0

  lat_idx = int((latitude - lat_min) / delta_lat) + 1
  lat_idx = max(1,min(lat_idx, num_lat_cirrus_ap))

  Tc_Apriori = t_tropo + TC_CIRRUS_MEAN_LAT_VECTOR(lat_idx)
  Tc_Apriori_Uncer = TC_CIRRUS_STDDEV_LAT_VECTOR(lat_idx)

  !--- values of the std dev are too small so use a fixed value for uncertainty
  !Tc_Apriori_Uncer = TC_AP_UNCER_CIRRUS_DEFAULT

end subroutine COMPUTE_CIRRUS_APRIORI

!----------------------------------------------------------------------
!--- compute the apriori from the cloud phase and etropo
!----------------------------------------------------------------------
subroutine COMPUTE_APRIORI_BASED_ON_TYPE( &
                           Cloud_Type, &
                           Latitude, &
                           Ttropo, &
                           T110um, &
                           Tc_Opaque_Lrc, &
                           Tc_Opaque, &
                           Mu, &
                           Tc_Ap, &
                           Tc_Ap_Uncer, &
                           Ec_Ap, &
                           Ec_Ap_Uncer, &
                           Beta_Ap,  &
                           Beta_Ap_Uncer, &
                           Success_Flag)

  integer(kind=int1), intent(in):: Cloud_Type
  real(kind=real4), intent(in):: Latitude
  real(kind=real4), intent(in):: Ttropo
  real(kind=real4), intent(in):: T110um
  real(kind=real4), intent(in):: Tc_Opaque_Lrc
  real(kind=real4), intent(in):: Tc_Opaque
  real(kind=real4), intent(in):: Mu
  real(kind=real4), intent(out):: Tc_Ap
  real(kind=real4), intent(out):: Ec_Ap
  real(kind=real4), intent(out):: Beta_Ap
  real(kind=real4), intent(out):: Tc_Ap_Uncer
  real(kind=real4), intent(out):: Ec_Ap_Uncer
  real(kind=real4), intent(out):: Beta_Ap_Uncer
  logical, intent(out):: Success_Flag

  real(kind=real4):: Tc_Ap_Cirrus
  real(kind=real4):: Tc_Ap_Uncer_Cirrus
  real(kind=real4):: Tc_Ap_Opaque
  real(kind=real4):: Emiss_Weight
  real(kind=real4):: Emiss_Weight2

  !--- calipso values (not multiplier on uncer values)
  call COMPUTE_CIRRUS_APRIORI(Ttropo, Latitude, Tc_Ap_Cirrus, Tc_Ap_Uncer_Cirrus)

  !Tc_Ap_Uncer_Cirrus = Tc_Ap_Uncer_Cirrus_Default

  !--- initialize with the opaque cloud temperature
  Tc_Ap_Opaque = T110um

  if (Tc_Opaque /= MISSING_VALUE_REAL4) then
     Tc_Ap_Opaque = Tc_Opaque
  endif

  if (Tc_Opaque_Lrc /= MISSING_VALUE_REAL4) then
      Tc_Ap_Opaque = Tc_Opaque_Lrc
  endif

  !-----NEW
  Tc_Ap =  MISSING_VALUE_REAL4
  Tc_Ap_Uncer =  MISSING_VALUE_REAL4
  Ec_Ap =  MISSING_VALUE_REAL4
  Ec_Ap_Uncer =  MISSING_VALUE_REAL4
  Beta_Ap =  MISSING_VALUE_REAL4
  Beta_Ap_Uncer = MISSING_VALUE_REAL4

  Success_Flag = .true.

  if (Cloud_Type == Symbol%Clear_Type) then
        return
  else if (Cloud_Type == Symbol%Prob_Clear_Type) then
        return
  else if (Cloud_Type == Symbol%Fog_Type) then
        Tc_Ap = Tc_Ap_Opaque
        Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
        Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Fog_Type/Mu)  !slow!
        Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
        Beta_Ap = Beta_Ap_Water
        Beta_Ap_Uncer = Beta_Ap_Uncer_Water
        return
   else if (Cloud_Type == Symbol%Water_Type) then
        Tc_Ap = Tc_Ap_Opaque
        Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
        Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Water_Type/Mu)  !slow!
        Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
        Beta_Ap = Beta_Ap_Water
        Beta_Ap_Uncer = Beta_Ap_Uncer_Water
        return
   else if (Cloud_Type == Symbol%Supercooled_Type) then
        Tc_Ap = Tc_Ap_Opaque
        Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
        Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Supercooled_Type/Mu)  !slow!
        Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
        Beta_Ap = Beta_Ap_Water
        Beta_Ap_Uncer = Beta_Ap_Uncer_Water
        return
   else if (Cloud_Type == Symbol%Mixed_Type) then
        Tc_Ap = Tc_Ap_Opaque
        Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
        Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Mixed_Type/Mu)  !slow!
        Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
        Beta_Ap = Beta_Ap_Water
        Beta_Ap_Uncer = Beta_Ap_Uncer_Water
        return
   else if (Cloud_Type == Symbol%Opaque_Ice_Type) then
        Tc_Ap = Tc_Ap_Opaque
        Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
        Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Opaque_Ice_Type/Mu)  !slow!
        Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
        Beta_Ap = Beta_Ap_Ice
        Beta_Ap_Uncer = Beta_Ap_Uncer_Ice
        return
   else if (Cloud_Type == Symbol%Cirrus_Type) then
        Tc_Ap = Tc_Ap_Cirrus
        Tc_Ap_Uncer = Tc_Ap_Uncer_Cirrus
        Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Cirrus_Type/Mu)  !slow!
        Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus
        Beta_Ap = Beta_Ap_Ice
        Beta_Ap_Uncer = Beta_Ap_Uncer_Ice
        return
   else if (Cloud_Type == Symbol%Overlap_Type) then
        Tc_Ap = Tc_Ap_Cirrus
        Tc_Ap_Uncer = Tc_Ap_Uncer_Cirrus
        Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Overlap_Type/Mu)  !slow!
        Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus
        Beta_Ap = Beta_Ap_Ice
        Beta_Ap_Uncer = Beta_Ap_Uncer_Ice
        return
   else if (Cloud_Type == Symbol%Overshooting_Type) then !used Opaque Ice
        Tc_Ap = Tc_Ap_Opaque
        Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
        Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Opaque_Ice_Type/Mu)  !slow!
        Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
        Beta_Ap = Beta_Ap_Ice
        Beta_Ap_Uncer = Beta_Ap_Uncer_Ice
        return
   else if (Cloud_Type == Symbol%Unknown_Type) then !used Mixed
        Tc_Ap = Tc_Ap_Opaque
        Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
        Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Mixed_Type/Mu)  !slow!
        Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
        Beta_Ap = Beta_Ap_Water
        Beta_Ap_Uncer = Beta_Ap_Uncer_Water
        return
   endif

   !--- if here, you failed
   Success_Flag = .false.

end subroutine COMPUTE_APRIORI_BASED_ON_TYPE

!----------------------------------------------------------------------
!--- compute the apriori from the cloud phase and etropo
!----------------------------------------------------------------------
subroutine COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO( &
                           Cloud_Phase, &
                           Emiss_110um_Tropo, &
                           Latitude, &
                           Ttropo, &
                           T110um, &
                           Tc_Opaque_Lrc, &
                           Tc_Opaque, &
                           Mu, &
                           Tc_Ap, &
                           Tc_Ap_Uncer, &
                           Ec_Ap, &
                           Ec_Ap_Uncer, &
                           Beta_Ap,  &
                           Beta_Ap_Uncer, &
                           Success_Flag)

  integer, intent(in):: Cloud_Phase
  real(kind=real4), intent(in):: Emiss_110um_Tropo
  real(kind=real4), intent(in):: Latitude
  real(kind=real4), intent(in):: Ttropo
  real(kind=real4), intent(in):: T110um
  real(kind=real4), intent(in):: Tc_Opaque_Lrc
  real(kind=real4), intent(in):: Tc_Opaque
  real(kind=real4), intent(in):: Mu
  real(kind=real4), intent(out):: Tc_Ap
  real(kind=real4), intent(out):: Ec_Ap
  real(kind=real4), intent(out):: Beta_Ap
  real(kind=real4), intent(out):: Tc_Ap_Uncer
  real(kind=real4), intent(out):: Ec_Ap_Uncer
  real(kind=real4), intent(out):: Beta_Ap_Uncer
  logical, intent(out):: Success_Flag

  real(kind=real4):: Tc_Ap_Cirrus
  real(kind=real4):: Tc_Ap_Uncer_Cirrus
  real(kind=real4):: Tc_Ap_Opaque
  real(kind=real4):: Emiss_Weight
  real(kind=real4):: Emiss_Weight2

  Success_Flag = .true.

  !--- calipso values (not multiplier on uncer values)
  call COMPUTE_CIRRUS_APRIORI(Ttropo, Latitude, Tc_Ap_Cirrus, Tc_Ap_Uncer_Cirrus)

  !--- initialize with the opaque cloud temperature
  Tc_Ap_Opaque = T110um

  if (Tc_Opaque /= MISSING_VALUE_REAL4) then
     Tc_Ap_Opaque = Tc_Opaque
  endif

  if (Tc_Opaque_Lrc /= MISSING_VALUE_REAL4) then
      Tc_Ap_Opaque = Tc_Opaque_Lrc
  endif

  if (Cloud_Phase /= Symbol%ICE_PHASE) then
    Tc_Ap = Tc_Ap_Opaque
    Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
    Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Water_Phase/Mu)  !slow!
    Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
    Beta_Ap = Beta_Ap_Water
    Beta_Ap_Uncer = Beta_Ap_Uncer_Water
  endif

  if (Cloud_Phase == Symbol%ICE_PHASE) then

    if (Emiss_110um_Tropo <= 0.0) then
            Emiss_Weight = 0.0
    elseif (Emiss_110um_Tropo > 1.0) then
            Emiss_Weight = 1.0
    else
            Emiss_Weight = Emiss_110um_Tropo
    endif

    Emiss_Weight2 = Emiss_Weight

    Tc_Ap = Emiss_Weight2*Tc_Ap_Opaque + &
            (1.0-Emiss_Weight2)*Tc_Ap_Cirrus

    Tc_Ap_Uncer = Emiss_Weight2*Tc_Ap_Uncer_Opaque + &
                  (1.0-Emiss_Weight2)*Tc_Ap_Uncer_Cirrus

    ! ignore weighting
    !Tc_Ap = Tc_Ap_Cirrus
    !Tc_Ap_Uncer = Tc_Ap_Uncer_Cirrus

    !---- for very thick clouds, we want to ignore the LRC to
    !---  to maintain spatial structure like overshooting columns
    if (Emiss_110um_Tropo > 0.95 .and. Tc_Opaque /= MISSING_VALUE_REAL4) then
      Tc_Ap = Tc_Opaque
      Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
    endif

    !---- for very thin clouds, ignore opaque solution
    if (Emiss_110um_Tropo < 0.5) then
      Tc_Ap = Tc_Ap_Cirrus
      Tc_Ap_Uncer = Tc_Ap_Uncer_Cirrus
    endif

    !--- emissivity and beta a priori
    Ec_Ap = min(0.99,max(0.1,Emiss_110um_Tropo))
    Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus
    Beta_Ap = Beta_Ap_Ice
    Beta_Ap_Uncer = Beta_Ap_Uncer_Ice

  endif

  end subroutine COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO
!---------------------------------------------------------------
! compute apriori without relaince on phase
!---------------------------------------------------------------
 subroutine COMPUTE_APRIORI_BASED_ON_TOPA( &
                       Topa, Emiss_110um_Tropo, Ice_Prob_Ap, &
                       Tc_Ap,Tc_Ap_Uncer, &
                       Ec_Ap,Ec_Ap_Uncer, &
                       Beta_Ap,Beta_Ap_Uncer,Success_Flag)

  real(kind=real4), intent(in):: Topa
  real(kind=real4), intent(in):: Ice_Prob_Ap
  real(kind=real4), intent(in):: Emiss_110um_Tropo

  real(kind=real4), intent(out):: Tc_Ap
  real(kind=real4), intent(out):: Ec_Ap
  real(kind=real4), intent(out):: Beta_Ap
  real(kind=real4), intent(out):: Tc_Ap_Uncer
  real(kind=real4), intent(out):: Ec_Ap_Uncer
  real(kind=real4), intent(out):: Beta_Ap_Uncer
  logical, intent(out):: Success_Flag

  Success_Flag = .true.

  Tc_Ap = Topa
  Tc_Ap_Uncer = 100.0 !Tc_Ap_Uncer_Cirrus

  Ec_Ap_Uncer = 1.0 !Ec_Ap_Uncer_Cirrus

  if (Ice_Prob_Ap .gtr. 0.5) then
    Ec_Ap = Emiss_110um_Tropo
    Beta_Ap = Beta_Ap_Ice
    Beta_Ap_Uncer = Beta_Ap_Uncer_Ice
  else
    Ec_Ap = 0.8
    Beta_Ap = Beta_Ap_Water
    Beta_Ap_Uncer = Beta_Ap_Uncer_Water
  endif

 end subroutine COMPUTE_APRIORI_BASED_ON_TOPA
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine SELECT_CHANNELS_EYRE_MENZEL(Input,Symbol)
   type(acha_input_struct), intent(inout) :: Input
   type(acha_symbol_struct), intent(in) :: Symbol

   N_Chan_EM = 0

   if (Input%Chan_On_062um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 37
   endif
   if (Input%Chan_On_067um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 27
   endif
   if (Input%Chan_On_073um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 28
   endif
   if (Input%Chan_On_085um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 29
   endif
   if (Input%Chan_On_097um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 30
   endif
   if (Input%Chan_On_104um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 38
   endif
   if (Input%Chan_On_110um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 31
      Chan_Idx_11um = N_Chan_EM
   endif
   if (Input%Chan_On_120um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 32
   endif
   if (Input%Chan_On_133um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 33
   endif
   if (Input%Chan_On_136um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 34
   endif
   if (Input%Chan_On_139um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 35
   endif
   if (Input%Chan_On_142um == symbol%YES) then
      N_Chan_EM = N_Chan_EM + 1
      EM_Chan_Idx(N_Chan_EM) = 36
   endif

   !--- if no absorbing channel, turn off EM
   if (Input%Chan_On_062um == symbol%NO .and. &
       Input%Chan_On_067um == symbol%NO .and. &
       Input%Chan_On_073um == symbol%NO .and. &
       Input%Chan_On_097um == symbol%NO .and. &
       Input%Chan_On_133um == symbol%NO .and. &
       Input%Chan_On_136um == symbol%NO .and. &
       Input%Chan_On_139um == symbol%NO .and. &
       Input%Chan_On_142um == symbol%NO) then

       N_Chan_Em = 0

   endif

end subroutine SELECT_CHANNELS_EYRE_MENZEL

!----------------------------------------------------------------------
! Eyre and Menzel
!----------------------------------------------------------------------
subroutine EYRE_MENZEL(Input,Symbol,Pc_EM,Tc_EM,Zc_EM,Ec_EM,N_EM,Res_EM,N_Std_EM,CV_EM,Ec_Res_EM)

  type(acha_input_struct), intent(inout) :: Input
  type(acha_symbol_struct), intent(in) :: Symbol
  real, intent(out), dimension(:,:):: Pc_EM
  real, intent(out), dimension(:,:):: Tc_EM
  real, intent(out), dimension(:,:):: Zc_EM
  real, intent(out), dimension(:,:):: Ec_EM
  real, intent(out), dimension(:,:):: N_EM
  real, intent(out), dimension(:,:):: Res_EM
  real, intent(out), dimension(:,:):: N_Std_EM
  real, intent(out), dimension(:,:):: CV_EM
  real, intent(out), dimension(:,:):: Ec_Res_EM
  integer:: Num_Elem, Num_Lines
  integer:: Nchan
  real:: numer_sum, denom_sum
  integer:: Chan_Idx
  integer:: Line_Idx,Elem_Idx,P_Lev_Idx,P_Lev_Idx_Opt
  integer:: Tropo_Level_Idx, Sfc_Level_Idx
  real:: N, Residual, P_Opt, N_Opt, Residual_Min,CV_Opt, Ec_Residual, Ec_Residual_Opt
  real:: delta_overcast_clear
  real:: delta_obs_clear
  integer:: Lev_Idx
  integer:: Nwp_Lon_Idx, Nwp_Lat_Idx, Vza_Rtm_Idx
  integer:: RTM_NWP_Error_Flag
  real:: N_Count, N_Sum, N_Sum_2
  real:: N_Temp, N_Mean, N_Std, N_Std_Opt, N_Mean_Opt
  real, parameter:: N_Min_Thresh = 0.01
  real, parameter:: N_Max_Thresh = 1.00
  real, parameter:: Rad_Ratio_Thresh = 1.00
  type(acha_rtm_nwp_struct) :: ACHA_RTM_NWP
  type(data_em) :: d_em

  Pc_EM = MISSING_VALUE_REAL4
  Ec_EM = MISSING_VALUE_REAL4
  Tc_EM = MISSING_VALUE_REAL4
  Ec_EM = MISSING_VALUE_REAL4
  N_EM = MISSING_VALUE_REAL4
  Res_EM = MISSING_VALUE_REAL4
  CV_EM = MISSING_VALUE_REAL4
  N_Std_EM = MISSING_VALUE_REAL4
  Ec_Res_EM = MISSING_VALUE_REAL4

  Line_Loop: do Line_Idx = 1, Input%Number_Of_Lines
  Element_Loop: do Elem_Idx = 1, Input%Number_Of_Elements

     if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) cycle

     !--- indice aliases
     Nwp_Lon_Idx = Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx)
     Nwp_Lat_Idx = Input%Line_Idx_Nwp(Elem_Idx,Line_Idx)
     Vza_Rtm_Idx = Input%Viewing_Zenith_Angle_Idx_Rtm(Elem_Idx,Line_Idx)

    !--------------------------------------------------------------------
    ! get profiles for this pixel
    !--------------------------------------------------------------------
    call ACHA_FETCH_PIXEL_NWP_RTM(Input, Symbol, &
                                 Elem_Idx,Line_Idx, &
                                 ACHA_RTM_NWP,RTM_NWP_Error_Flag)
    if (RTM_NWP_Error_Flag /= 0) cycle


     Tropo_Level_Idx = ACHA_RTM_NWP%Tropo_Level
     Sfc_Level_Idx = ACHA_RTM_NWP%Sfc_Level

     !--- intialize residual with large number
     Residual_Min = huge(Residual_Min)
     P_Lev_Idx_Opt = Sfc_Level_Idx

     call d_em % set(Elem_Idx,Line_Idx,Input,ACHA_RTM_NWP,N_Chan_EM)
     Press_Loop: do P_Lev_Idx = Tropo_Level_Idx, Sfc_Level_Idx

       numer_sum = 0.0
       denom_sum = 0.0

       Chan_Loop: do Chan_Idx = 1, N_Chan_EM

         delta_overcast_clear = (d_em % ch(Chan_Idx)%  BB_Rad(P_Lev_Idx) &
            - d_em % ch(Chan_Idx)%Clr_Rad)

         delta_obs_clear = (d_em % ch(Chan_Idx)%Obs_Rad - d_em % ch(Chan_Idx)%Clr_Rad)

         !--- remove channels without signal
         if (100.0*abs(delta_overcast_clear/d_em % ch(Chan_Idx)%Clr_Rad) .ltr. Rad_Ratio_Thresh) cycle

         numer_sum = numer_sum + delta_obs_clear *  delta_overcast_clear
         denom_sum = denom_sum + delta_overcast_clear**2

       end do Chan_Loop

       N = numer_sum / denom_sum

       if (P_Lev_Idx == Tropo_Level_Idx .and. N < N_Min_Thresh) exit

       if (N < N_Min_Thresh .or. N > N_Max_Thresh) cycle

       numer_sum = 0.0
       denom_sum = 0.0
       N_Count = 0.0
       N_Sum = 0.0
       N_Sum_2 = 0.0

       Chan_Loop_2: do Chan_Idx = 1, N_Chan_EM

         delta_overcast_clear = (d_em % ch(Chan_Idx)%BB_Rad(P_Lev_Idx) - d_em % ch(Chan_Idx)%Clr_Rad)

         delta_obs_clear = (d_em % ch(Chan_Idx)%Obs_Rad - d_em % ch(Chan_Idx)%Clr_Rad)

         !--- remove channels without signal
         if (100.0*abs(delta_overcast_clear/Clr_Rad) .ltr. Rad_Ratio_Thresh) cycle

         numer_sum = numer_sum + delta_obs_clear**2 - (N**2)*(delta_overcast_clear)**2
         denom_sum = denom_sum + abs(delta_obs_clear/delta_overcast_clear - N) !mean emiss diff

         N_Temp =  delta_obs_clear / delta_overcast_clear
         N_Sum =  N_Sum + N_Temp
         N_Sum_2 = N_Sum_2 + N_Temp**2
         N_Count = N_Count + 1

       end do Chan_Loop_2

       Residual = numer_sum
       Ec_Residual = denom_sum / N_Count
       N_Mean = N_Sum/N_Count
       N_Std = sqrt(N_Sum_2 / N_Count) - N_Mean**2

       if (Residual < Residual_Min) then
           P_Opt = ACHA_RTM_NWP%P_Prof(P_Lev_Idx)
           P_Lev_Idx_Opt = P_Lev_Idx
           N_Opt = N
           Residual_Min = Residual
           CV_Opt = N_Std / N_Mean
           N_Std_Opt = N_Std
           Ec_Residual_Opt = Ec_Residual
       endif

     end do Press_Loop
     !---- Check is at surface
     if (P_Lev_Idx_Opt >= Sfc_Level_Idx) cycle
     if (N_Opt > N_Max_Thresh) cycle
     if (N_Opt < N_Min_Thresh) cycle

     !-----
     Pc_EM(Elem_Idx,Line_Idx) = P_Opt
     N_EM(Elem_Idx,Line_Idx) = N_Opt   !!!!!
     Res_EM(Elem_Idx,Line_Idx) = Residual_Min
     N_Std_EM(Elem_Idx,Line_Idx) = N_Std_Opt
     CV_EM(Elem_Idx,Line_Idx) = CV_Opt
     Ec_Res_EM(Elem_Idx,Line_Idx) = Ec_Residual_Opt

     !--- compute temperature and height
     call KNOWING_P_COMPUTE_T_Z(ACHA_RTM_NWP,Pc_EM(Elem_Idx,Line_Idx), &
                                Tc_EM(Elem_Idx,Line_Idx),Zc_EM(Elem_Idx,Line_Idx),Lev_Idx)

     !-- compute emissivity
     call SELECT_CHAN_RAD(Elem_Idx,Line_Idx,Chan_Idx_11um,Input,ACHA_RTM_NWP)

     Ec_EM(Elem_Idx,Line_Idx) = EMISSIVITY(Obs_Rad, Clr_Rad, BB_Rad(Lev_Idx))

  end do Element_Loop
  end do Line_Loop

end subroutine EYRE_MENZEL
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------

subroutine data_em_set (self, Elem_Idx,Line_Idx,Input, ACHA_RTM_NWP,chn_max )
  class(data_em) :: self
  integer, intent(in):: Elem_Idx,Line_Idx
  type(acha_input_struct),intent(in) :: Input
  type(acha_rtm_nwp_struct),intent(in) :: ACHA_RTM_NWP
  integer :: chn_max
  integer :: ch_idx

  do ch_idx =1,chn_max
    call SELECT_CHAN_RAD(Elem_Idx,Line_Idx,Ch_Idx,Input,ACHA_RTM_NWP)
    self % ch(ch_idx) % bb_rad = bb_rad
    self % ch(ch_idx) % clr_rad = clr_rad
    self % ch(ch_idx) % obs_rad = obs_rad
    self % ch % is_set = .true.
  end do
 self%is_set = .true.
end subroutine


!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
subroutine SELECT_CHAN_RAD(Elem_Idx,Line_Idx,Chan_Idx,Input,ACHA_RTM_NWP)
  integer, intent(in):: Elem_Idx,Line_Idx,Chan_Idx
  type(acha_input_struct),intent(in) :: Input
  type(acha_rtm_nwp_struct),intent(in) :: ACHA_RTM_NWP

  !call ts%tic(1)

         select case (EM_Chan_Idx(Chan_Idx))

           case(37)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_062um
            Clr_Rad = Input%Rad_Clear_062um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_062um(Elem_Idx,Line_Idx)

           case(27)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_067um
            Clr_Rad = Input%Rad_Clear_067um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_067um(Elem_Idx,Line_Idx)

           case(28)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_073um
            Clr_Rad = Input%Rad_Clear_073um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_073um(Elem_Idx,Line_Idx)

           case(29)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_085um
            Clr_Rad = Input%Rad_Clear_085um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_085um(Elem_Idx,Line_Idx)

           case(30)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_097um
            Clr_Rad = Input%Rad_Clear_097um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_097um(Elem_Idx,Line_Idx)

           case(38)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_104um
            Clr_Rad = Input%Rad_Clear_104um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_104um(Elem_Idx,Line_Idx)

           case(31)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_110um
            Clr_Rad = Input%Rad_Clear_110um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_110um(Elem_Idx,Line_Idx)

           case(32)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_120um
            Clr_Rad = Input%Rad_Clear_120um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_120um(Elem_Idx,Line_Idx)

           case(33)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_133um
            Clr_Rad = Input%Rad_Clear_133um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_133um(Elem_Idx,Line_Idx)

           case(34)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_136um
            Clr_Rad = Input%Rad_Clear_136um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_136um(Elem_Idx,Line_Idx)

           case(35)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_139um
            Clr_Rad = Input%Rad_Clear_139um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_139um(Elem_Idx,Line_Idx)

           case(36)
            BB_Rad = ACHA_RTM_NWP%Black_Body_Rad_Prof_142um
            Clr_Rad = Input%Rad_Clear_142um(Elem_Idx,Line_Idx)
            Obs_Rad = Input%Rad_142um(Elem_Idx,Line_Idx)

    end select

    !call ts%tac(1)

 end subroutine SELECT_CHAN_RAD

    !====================================================================
   ! Function Name: EMISSIVITY
   !
   ! Function:
   !  Computes the  effective emissivity
   !
   ! Input:  Rad_Toa - channel radiance at top of atmosphere(toa)
   !         Rad_Clear_Tau - channel radiance at toa for clear skies
   !         Rad_Cloud_BB_Toa - channel radiance at TOA if cloud were a
   !         Black-Body
   !
   ! Output: Emiss - the effective cloud emissivity
   !
   !====================================================================
   function EMISSIVITY(Radiance_Toa, Radiance_Clear_Toa, Radiance_Cloud_BB_Toa) result(Emiss)
      real(kind=real4), intent(in) :: Radiance_Toa
      real(kind=real4), intent(in) :: Radiance_Clear_Toa
      real(kind=real4), intent(in) :: Radiance_Cloud_BB_Toa
      real(kind=real4) :: Emiss

      Emiss = Missing_Value_Real4

      if (Radiance_Cloud_BB_Toa /= Radiance_Clear_Toa) then
          Emiss = (Radiance_Toa - Radiance_Clear_Toa) / &
            (Radiance_Cloud_BB_Toa - Radiance_Clear_Toa)
       end if

      return

   end function EMISSIVITY


end module ACHA_PRIOR_MODULE
