!$Id: acha_module.f90 4092 2021-03-02 22:05:14Z heidinger $
module AWG_CLOUD_HEIGHT
!---------------------------------------------------------------------
! This module houses the routines associated with...
!
! ACHA - AWG Cloud Height Algorithm
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
! Reference:
!
!  Heidinger, A.K., and M.J. Pavolonis, 2009: Gazing at Cirrus Clouds for 25 Years
!  through a Split Window. Part I: Methodology. J. Appl. Meteor. Climatol., 48,
!  1100-1116.
!
!  Heidinger, A. K.; Pavolonis, M. J.; Holz, R. E.; Baum, Bryan A. and Berthier,
!  S.. Using CALIPSO to explore the sensitivity to cirrus height in the infrared
!  observations from NPOESS/VIIRS and GOES-R/ABI. Journal of Geophysical
!  Research, Volume 115, 2010, Doi:10.1029/2009JD012152. 
!
! Meta Data Flags
! 1 - Cloud Height Attempted (0 = no / 1 = yes)
! 2 - Bias Correction Employed (0 = no / 1 = yes)
! 3 - Ice Cloud Retrieval (0 = no / 1 = yes)
! 4 - Local Radiatve Center Processing Used (0 = no / 1 = yes)
! 5 - Multi-layer Retrieval (0 = no / 1 = yes)
! 6 - Lower Cloud Interpolation Used (0 = no / 1 = yes)
! 7 - Boundary Layer Inversion Assumed  (0 = no / 1 = yes)
! 8 - NWP Profile Inversion Assumed (0 = no / 1 = yes)
!
! Packed Quality Flags
! 1 - Processed (0 = no / 1 = yes)
! 2 - Valid Tc Retrieval (0 = yes, 1 = no)
! 3 - Valid ec Retrieval (0 = yes, 1 = no)
! 4 - Valid beta Retrieval (0 = yes, 1 = no)
! 5 - degraded Tc Retrieval (0 = no, 1 = yes)
! 6 - degraded ec Retrieval (0 = no, 1 = yes)
! 7 - degraded beta Retrieval (0 = no, 1 = yes)
!
! Modes
!  off   - acha does not run
!  default - mode is set by clavr-x defaults
!  maximum - all supported channels (experimental)
!  110  - only 11 micron (opaque solution)
!  038_110 - 3.8 and 11 micron
!  067_110 - 6.7 and 11 micron
!  110_120
!  110_133
!  067_085_110
!  067_110_120
!  067_110_133
!  085_110_120
!  110_120_133
!  067_085_110_120
!  085_110_120_133
!  062_085_110_120_133 - stw
!  067_085_110_120_133
!  110_133_136_139_142
!  085_110_120_133_136_139_142 - 8.5,11,12,13.4,13.6,13.9 and 14.2 micron
!
! MULTI_LAYER_LOGIC_FLAG
! 0 - (baseline) just use the multilayer id in cloud type
! 1 - treat all multilayer like cirrus
! 2 - assume all cirrus are multilayer and let acha decide
!
!----------------------------------------------------------------------
! Ice Fraction Algorithm LUT Variables
!----------------------------------------------------------------------
!Changes needed to get into SAPF
!
! - Renamed AWG_CLOUD_HEIGHT_ALGORITHM to AWG_CLOUD_HEIGHT_ALGORITHM_ACHA
! - Renamed LOCAL_LINEAR_RADIATIVE_CENTER to LOCAL_LINEAR_RADIATIVE_CENTER_ACHA
! - Renamed module from AWG_CLOUD_HEIGHT to AWG_CLOUD_HEIGHT_ACHA
! - Had to redo Skip_LRC_Mask due to issues in Framework
!
! ** Note:  These changes are in the Framework repository only.
!
!----------------------------------------------------------------------
  use ACHA_SERVICES_MOD, only : &
           real4, int1, int4, real8, dtor, &
           Acha_output_struct,ACHA_SYMBOL_STRUCT, &
           Acha_input_struct, Acha_rtm_nwp_struct, &
           PLANCK_RAD_FAST, PLANCK_TEMP_FAST, &
           INVERT_MATRIX, ACHA_FETCH_PIXEL_NWP_RTM, &
           LOCATE, Acha_Diag_Struct, Acha_Dump_Struct, COUNTSUBSTRING, &
           ABI_Use_104um_Flag

  use ACHA_FULL_RETRIEVAL_MOD
  use ACHA_RTM_MOD
  use ACHA_NUM_MOD
  use ACHA_LHP_MOD
  use ACHA_MICROPHYSICAL_MODULE
  use ACHA_ICE_FRACTION_MODULE
  use CX_REAL_BOOLEAN_MOD
  use KDTREE2_MODULE

  implicit none

  public:: AWG_CLOUD_HEIGHT_ALGORITHM
  public:: CHECK_ACHA_MODE
  public:: SET_ACHA_VERSION
  public:: LOCAL_LINEAR_RADIATIVE_CENTER

  private:: DETERMINE_SFC_TYPE_FORWARD_MODEL

  private:: DETERMINE_ACHA_MODE_BASED_ON_CHANNELS
  private:: COMPUTE_REFERENCE_LEVEL_EMISSIVITY
  private:: COMPUTE_STANDARD_DEVIATION
  private:: NULL_PIX_POINTERS 
  private:: EMPIRICAL_LAPSE_RATE
  private:: COMPUTE_Y
  private:: COMPUTE_META_DATA
  private:: COMPUTE_HEIGHT_FROM_LAPSE_RATE
  private:: QUALITY_CONTROL_OUTPUT
  private:: TEST_INPUT
  private:: SET_OUTPUT_PACKED_QF
  private:: COMPUTE_SPECTRAL_CLOUD_EMISSIVITY
  private:: SET_SURFACE_EMISSIVITY
  private:: APPLY_OPAQUE_RETRIEVAL
  private:: COMPUTE_PHASE
  private:: COMPUTE_SA
  private:: CONVERT_TC_TO_PC_AND_ZC
  private:: SAVE_X_2_OUTPUT
  private:: DETERMINE_ACHA_CLOUD_TYPE

  !--- include the non-system specific variables
  include 'include/acha_parameters.inc'

  real, private:: Bt_110um_Bt_110um_Covar

  real, private, PARAMETER:: MISSING_VALUE_REAL4 = -999.0
  integer(kind=int1), private, PARAMETER:: MISSING_VALUE_integer1 = -128_int1
  !integer(kind=int1), private, PARAMETER:: MISSING_VALUE_integer1 = -128
  integer(kind=int4), private, PARAMETER:: MISSING_VALUE_integer4 = -999
  type(ACHA_SYMBOL_STRUCT), private :: Symbol

  integer, public, parameter:: Num_ACHA_Modes = 20
  integer, public, parameter:: ACHA_Mode_Max_Length = 31
  character(len=ACHA_Mode_Max_Length), dimension(Num_ACHA_Modes), public, parameter:: ACHA_Mode_Values = &
     (/'off                            ', &
       'default                        ', &
       'maximum                        ', &
       '110                            ', &
       '038_110                        ', &
       '067_110                        ', &
       '110_120                        ', &
       '110_133                        ', &
       '067_085_110                    ', &
       '067_110_120                    ', &
       '067_110_133                    ', &
       '085_110_120                    ', &
       '110_120_133                    ', &
       '067_085_110_120                ', &
       '085_110_120_133                ', &
       '062_085_110_120_133            ', &
       '067_085_110_120_133            ', &
       '110_133_136_139_142            ', &
       '085_110_120_133_136_139_142    ', &
       '062_067_073_085_104_110_120_133'/)

  !-------------------------------------------------------------------------------------
  ! empirical lapse rate table data and metadata
  !-------------------------------------------------------------------------------------
  integer, private, parameter:: nts = 7
  integer, private, parameter:: ntcs = 9
  real, private, parameter:: ts_min = 270.0
  real, private, parameter:: dts = 5.0
  real, private, parameter:: tcs_min = -20.0
  real, private, parameter:: dtcs = 2.0

  !----- new 10/2020
  real, private, dimension(nts,ntcs), parameter:: ocean_lapse_rate_table = reshape ((/ &
  -9.7, -9.4, -9.4, -9.4, -9.2, -8.1, -7.2, &
  -9.9, -9.5, -9.5, -9.6, -9.4, -8.4, -7.4, &
 -10.0, -9.7, -9.6, -9.9, -9.6, -8.5, -7.6, &
 -10.0, -9.7, -9.6, -9.9, -9.6, -8.7, -8.0, &
  -9.6, -9.5, -9.7, -9.9, -9.6, -9.1, -8.7, &
  -9.1, -9.3, -9.5, -9.7, -9.6, -9.6, -9.7, &
  -8.9, -9.1, -9.2, -9.3, -9.5,-10.0,-10.3, &
  -8.1, -8.2, -8.2, -8.4, -8.8, -9.6,-10.3, &
  -7.5, -7.4, -7.3, -7.5, -8.2, -9.1, -9.9/), (/nts,ntcs/))

  real, private, dimension(nts,ntcs), parameter:: land_lapse_rate_table = reshape ((/ &
  -5.9, -6.2, -6.6, -6.9, -7.4, -8.2, -8.9, &
  -5.8, -6.1, -6.5, -6.8, -7.4, -8.3, -9.1, &
  -5.8, -6.0, -6.4, -6.7, -7.2, -8.2, -9.1, &
  -5.7, -5.9, -6.3, -6.5, -7.0, -8.1, -9.1, &
  -5.7, -5.9, -6.2, -6.3, -6.7, -7.8, -8.9, &
  -5.7, -5.9, -6.1, -6.3, -6.5, -7.6, -8.5, &
  -5.8, -5.9, -5.9, -5.9, -6.1, -7.4, -8.7, &
  -5.2, -5.4, -5.4, -5.5, -5.6, -6.8, -7.8, &
  -4.6, -4.9, -4.9, -5.0, -5.1, -6.2, -7.1/), (/nts,ntcs/))

  contains 

!------------------------------------------------------------------------------
! AWG Cloud Height Algorithm (ACHA)
!
! Author: Andrew Heidinger, NOAA
!
! Assumptions
!   1) No scattering
!   2) single layer cloud for cloud type /= 6
!   3) for overlap type, an opaque cloud 200 mb above the surface lies below 
!
! Limitations
!   1) sensitivity to Tc is low for thin clouds
!   2) little Emissivity sensitivity for low clouds
!
! input to the retrieval
!  y(1) = t4 - 11 micron brightness temperature
!  y(2) = t4 - t5 - the split window temperature
!
! the output of the retrieval
!  x(1) - the cloud temperature
!  x(2) - the 11 micron Emissivity at nadir
!  x(3) - the beta ratio for 11 and 12 microns
!  x(4) - the surface (or lower cloud) temperature
!  x(5) - the ice fraction
!
! This routine uses a 1d-var retrieval approach as outlined in Rodger (1976)
!
! input to the 1d-var approach
!  y - the vector of observations
!  x_Ap - the a apriori estimates of x
!  Sa - the error covariance matric of x_Ap
!  Sy - the error covariance of y (included calibration, forward model)
!
!  the Optimal Estimation Quality Flags are determined as follows
!  3 - estimated error < 1/3 a priori error
!  2 - estimated error < 2/3 a priori error
!  1 - any other converged retrieval
!  0 - a failed or unattempted retrieval
!
! the overall quality flag Description
! 0 - No retrieval attempted
! 1 - Retrieval attempted and failed
! 2 - Marginally Successful Retrieval
! 3 - Fully Successful Retrieval
!
!
! Meta Data
! 1 - Cloud Height Attempted (0 = no / 1 = yes)
! 2 - Bias Correction Employed (0 = no / 1 = yes)
! 3 - Ice Cloud Retrieval (0 = no / 1 = yes)
! 4 - Local Radiatve Center Processing Used (0 = no / 1 = yes)
! 5 - Multi-layer Retrieval (0 = no / 1 = yes)
! 6 - Lower Cloud InterpoLation Used (0 = no / 1 = ! yes)
! 7 - Boundary Layer Inversion Assumed  (0 = ! no / 1 = yes)
!
! Processing Order Description
! 0 = Not Processed
! 1 = non-multi-layer lrc pixels
! 2 = single layer water cloud pixels
! 3 = lrc multi-layer clouds
! 4 = all remaining clouds
! 5 = if USE_CIRRUS_FLAG is set on, redo all thin cirrus using a priori
!          temperature from thicker cirrus.
!
! Proposed
! 0 = Not Processed
! 1 = lrc pixels
! 2 = pixels associated with lrc
! 3 = cirrus assumed to be multilayer - acha determines multilayer
! 4 = opposite phase and replace if cost lower
! 
!
!
!----------------------------------------------------------------------
! modification history
!
!  March 2020 - major rewrite, remove processing order
!
!------------------------------------------------------------------------------
  subroutine  AWG_CLOUD_HEIGHT_ALGORITHM(Input, Symbol_In, Output, Diag, Dump)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================

  type(acha_input_struct), intent(inout) :: Input
  type(acha_symbol_struct), intent(in) :: Symbol_In
  type(acha_output_struct), intent(inout) :: Output
  type(acha_diag_struct), intent(inout), optional :: Diag
  type(acha_dump_struct), intent(in), optional :: Dump
  integer, save:: Diag_Warning_Flag = 0

  !===============================================================================
  !  Pixel level RTM structure
  !===============================================================================
 
  type(acha_rtm_nwp_struct) :: ACHA_RTM_NWP

  !===============================================================================
  !  Local Variable Declaration
  !===============================================================================

  character(len=50):: ACHA_Mode_Flag
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Param_Idx
  integer:: i
  integer:: Singular_Flag
  integer:: Lev_Idx
  integer:: Ivza
  integer:: ilrc
  integer:: jlrc
  integer:: Iter_Idx
  integer:: ierror

  real:: Convergence_Criteria
  real:: Convergence_Criteria_Simple

  real:: Tsfc_Est
  real:: Tc_temp
  real:: Zc_Temp

  real:: Ts_Ap
  real:: Ts_Ap_Uncer
  real:: Tc_Ap_Imager
  real:: Tc_Ap_Sounder
  
  integer(kind=int1):: I1_Dummy
  real:: Emiss_110um_Tropo
  integer:: Num_Obs
  integer(kind=int1):: Cloud_Type
  integer:: Cloud_Phase
  integer:: Sfc_Type_Forward_Model
  integer(kind=int1), dimension(NUM_META_DATA):: Meta_Data_Flags

  !--- 1d-var retrieval arrays
  real (kind=real4), allocatable, dimension(:):: y
  real (kind=real4), allocatable, dimension(:):: f
  real (kind=real4), allocatable, dimension(:):: y_variance

  real (kind=real4), allocatable, dimension(:):: x
  real (kind=real4), allocatable, dimension(:):: x_Ap
  real (kind=real4), allocatable, dimension(:,:):: Sa
  real (kind=real4), allocatable, dimension(:,:):: Sa_Inv
  real (kind=real4), allocatable, dimension(:,:):: Sx
  real (kind=real4), allocatable, dimension(:,:):: AKM

  real (kind=real4), allocatable, dimension(:):: x_Simple
  real (kind=real4), allocatable, dimension(:):: x_Ap_Simple
  real (kind=real4), allocatable, dimension(:,:):: Sa_Simple
  real (kind=real4), allocatable, dimension(:,:):: Sa_Inv_Simple
  real (kind=real4), allocatable, dimension(:,:):: Sx_Simple
  real (kind=real4), allocatable, dimension(:,:):: AKM_Simple

  integer(kind=int4), dimension(:,:), allocatable:: Fail_Flag
  integer(kind=int4), dimension(:,:), allocatable:: Fail_Flag_Simple
  integer(kind=int4), dimension(:,:), allocatable:: Fail_Flag_Full
  integer(kind=int4), dimension(:,:), allocatable:: Converged_Flag
  real (kind=real4), allocatable, dimension(:,:):: Temperature_Cirrus
  integer (kind=int1), allocatable, dimension(:,:):: Cloud_Type_Temp

  !--- local POINTERs to global arrays or data structures
  integer(kind=int4), allocatable, dimension(:,:):: Elem_Idx_LRC
  integer(kind=int4), allocatable, dimension(:,:):: Line_Idx_LRC
  integer(kind=int1), allocatable, dimension(:,:):: Skip_LRC_Mask
  real (kind=real4):: Tc_Opaque_Lrc
  real (kind=real4):: Bt_110um_Lrc
  real (kind=real4):: T_Tropo
  real (kind=real4):: Z_Tropo
  real (kind=real4):: P_Tropo

  !--- scalar local variables
  integer (kind=int4):: NWP_Profile_Inversion_Flag
  integer (kind=int4):: Dummy_Flag
  logical:: Bad_Input_Flag
  logical:: Clip_Output_Flag
  logical :: Singular_Warning_First_Time 
  logical :: Prior_Success_Flag
  
  !--- indices for single pixel diagnostic dump
  integer:: Elem_Abs_Idx 
  integer:: Line_Abs_Idx 
  integer:: Elem_Abs_Idx_Dump
  integer:: Line_Abs_Idx_Dump
  integer:: Lun_Prof_Dump
  integer:: Lun_Iter_Dump
  integer:: Lun_Warpup_Dump
  logical:: Dump_Diag
  character(len=100), parameter:: File_Name_Prof_Dump = "acha_profile_dump.txt"
  character(len=100), parameter:: File_Name_Iter_Dump = "acha_iteration_dump.txt"
  integer:: RTM_NWP_Error_Flag

  ! surface emissivity
  real(kind=real4):: Emiss_Sfc_038um
  real(kind=real4):: Emiss_Sfc_062um
  real(kind=real4):: Emiss_Sfc_067um
  real(kind=real4):: Emiss_Sfc_073um
  real(kind=real4):: Emiss_Sfc_085um
  real(kind=real4):: Emiss_Sfc_097um
  real(kind=real4):: Emiss_Sfc_104um
  real(kind=real4):: Emiss_Sfc_110um
  real(kind=real4):: Emiss_Sfc_120um
  real(kind=real4):: Emiss_Sfc_133um
  real(kind=real4):: Emiss_Sfc_136um
  real(kind=real4):: Emiss_Sfc_139um
  real(kind=real4):: Emiss_Sfc_142um

  real(kind=real4):: Zc_Thick

!-----------------------------------------------------------------------
! BEGIN EXECUTABLE CODE
!-----------------------------------------------------------------------
   Singular_Warning_First_Time = .true.
   Lun_Prof_Dump = -1
   Lun_Iter_Dump = -1

  !--------------------------------------------------------------------
  ! set up reference channel (11 or 10.4 micron)
  !--------------------------------------------------------------------
  call SETUP_REFERENCE_CHANNEL(ABI_Use_104um_Flag,Input)

  !--------------------------------------------------------------------
  ! set up noise table
  !--------------------------------------------------------------------
  call COMPUTE_BT_NOISE_TABLE(Input%WMO_Id, &
                                  Input%Chan_Idx_038um, Input%Chan_On_038um, &
                                  Input%Chan_Idx_062um, Input%Chan_On_062um, &
                                  Input%Chan_Idx_067um, Input%Chan_On_067um, &
                                  Input%Chan_Idx_073um, Input%Chan_On_073um, &
                                  Input%Chan_Idx_085um, Input%Chan_On_085um, &
                                  Input%Chan_Idx_097um, Input%Chan_On_097um, &
                                  Input%Chan_Idx_104um, Input%Chan_On_104um, &
                                  Input%Chan_Idx_110um, Input%Chan_On_110um, &
                                  Input%Chan_Idx_120um, Input%Chan_On_120um, &
                                  Input%Chan_Idx_133um, Input%Chan_On_133um, &
                                  Input%Chan_Idx_136um, Input%Chan_On_136um, &
                                  Input%Chan_Idx_139um, Input%Chan_On_139um, &
                                  Input%Chan_Idx_142um, Input%Chan_On_142um)

  !--------------------------------------------------------------------
  ! copy Symbol to a module-wide structure
  !--------------------------------------------------------------------
  Symbol = Symbol_In 

  !----------------------------------------------------------------------------
  ! abort if no 10.4 um or 11 um channel
  !----------------------------------------------------------------------------
  if (ABI_Use_104um_Flag) then
    if (Input%Chan_On_104um == Symbol%NO) THEN
      Output%Packed_Qf = 0_int1
      return
    endif
  else
    if (Input%Chan_On_110um == Symbol%NO) THEN
      Output%Packed_Qf = 0_int1
      return
    endif
  endif

  !--- initialize diagnostic output
  if (present(Diag) .and. Diag_Warning_Flag == 0) then
      print *, "CLAVR-x / ACHA ===>  Diagnostic Output Turned On"
      Diag_Warning_Flag = 1
  endif
  if (present(Diag)) Diag%Array_1 = MISSING_VALUE_REAL4
  if (present(Diag)) Diag%Array_2 = MISSING_VALUE_REAL4
  if (present(Diag)) Diag%Array_3 = MISSING_VALUE_REAL4
  
  !---------------------------------------------------------------------------
  !-- setup microphysical models
  !---------------------------------------------------------------------------
  call SETUP_ICE_MICROPHYSICAL_MODEL(Input%WMO_Id)

  !---------------------------------------------------------------------------
  !-- Acha Mode set to  -1, determine based on channels
  !---------------------------------------------------------------------------
  ACHA_Mode_Flag = Input%ACHA_Mode_Flag_In

  !--- This won't get called for bad GOES-17 data.  It is set using
  !--- MODIFY_MODE_USING_LHP_THRESHOLDS, called from the bridge.

  if (trim(ACHA_Mode_Flag) == 'unknown') then
    call DETERMINE_ACHA_MODE_BASED_ON_CHANNELS( &
                                      Acha_Mode_Flag, &
                                      Input%Chan_On_038um, &
                                      Input%Chan_On_062um, &
                                      Input%Chan_On_067um, &
                                      Input%Chan_On_085um,  &
                                      Input%Chan_On_110um,  &
                                      Input%Chan_On_120um,  &
                                      Input%Chan_On_133um,  &
                                      Input%Chan_On_136um,  &
                                      Input%Chan_On_139um,  &
                                      Input%Chan_On_142um)
  endif

  call DETERMINE_NUMBER_OF_CHANNELS(Acha_Mode_Flag, Num_Obs)

  !--- allocate needed 2d arrays for processing this segment
  allocate(Elem_Idx_LRC(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Line_Idx_LRC(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Skip_LRC_Mask(Input%Number_of_Elements,Input%Number_of_Lines))

  !--- allocate array for cirrus temperature
  allocate(Fail_Flag(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Fail_Flag_Simple(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Fail_Flag_Full(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Converged_Flag(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Temperature_Cirrus(Input%Number_of_Elements,Input%Number_of_Lines))

  allocate(Cloud_Type_Temp(Input%Number_of_Elements,Input%Number_of_Lines))

  !--- allocate 1D-VAR arrays based on number of channels
  allocate(y(Num_Obs))
  allocate(y_variance(Num_Obs))
  allocate(f(Num_Obs))

  !--- allocate 1D-VAR arrays based on number of parameters
  allocate(x(Num_Param))
  allocate(x_Ap(Num_Param))
  allocate(Sa(Num_Param,Num_Param))
  allocate(Sa_inv(Num_Param,Num_Param))
  allocate(Sx(Num_Param,Num_Param))
  allocate(AKM(Num_Param,Num_Param))

  allocate(x_Simple(Num_Param_Simple))
  allocate(x_Ap_Simple(Num_Param_Simple))
  allocate(Sa_Simple(Num_Param_Simple,Num_Param_Simple))
  allocate(Sa_inv_Simple(Num_Param_Simple,Num_Param_Simple))
  allocate(Sx_Simple(Num_Param_Simple,Num_Param_Simple))
  allocate(AKM_Simple(Num_Param_Simple,Num_Param_Simple))

  !--- set convergence criterion
  Convergence_Criteria = (Num_Param - 1.0) / 5.0
  Convergence_Criteria_Simple = (Num_Param_Simple - 1.0) / 5.0

  !--- initialize output
  Output%Tc =  MISSING_VALUE_REAL4
  Output%Ec =  MISSING_VALUE_REAL4
  Output%Beta =  MISSING_VALUE_REAL4
  Output%Pc =  MISSING_VALUE_REAL4
  Output%Zc =  MISSING_VALUE_REAL4
  Output%OE_Qf = 0_int1
  Output%Qf = MISSING_VALUE_integer1  !0_int1
  Meta_Data_Flags = 0_int1
  Output%Inversion_Flag = 0_int1
  if (Input%Chan_On_038um == Symbol%YES) Output%Ec_038um = MISSING_VALUE_REAL4
  if (Input%Chan_On_067um == Symbol%YES) Output%Ec_067um = MISSING_VALUE_REAL4
  if (Input%Chan_On_085um == Symbol%YES) Output%Ec_085um = MISSING_VALUE_REAL4
  if (Input%Chan_On_097um == Symbol%YES) Output%Ec_097um = MISSING_VALUE_REAL4
  if (Input%Chan_On_104um == Symbol%YES) Output%Ec_104um = MISSING_VALUE_REAL4
  if (Input%Chan_On_110um == Symbol%YES) Output%Ec_110um = MISSING_VALUE_REAL4
  if (Input%Chan_On_120um == Symbol%YES) Output%Ec_120um = MISSING_VALUE_REAL4
  if (Input%Chan_On_133um == Symbol%YES) Output%Ec_133um = MISSING_VALUE_REAL4
  if (Input%Chan_On_136um == Symbol%YES) Output%Ec_136um = MISSING_VALUE_REAL4
  if (Input%Chan_On_139um == Symbol%YES) Output%Ec_139um = MISSING_VALUE_REAL4
  if (Input%Chan_On_142um == Symbol%YES) Output%Ec_142um = MISSING_VALUE_REAL4

  !--------------------------------------------------------------------------
  ! spatial processing pixels
  ! compute local radiative centers using 11 um brightness temperature
  !---------------------------------------------------------------------------

  !--- construct a mask to select pixel for LRC computation
  Elem_Idx_LRC = MISSING_VALUE_INTEGER4
  Line_Idx_LRC = MISSING_VALUE_INTEGER4
  Skip_LRC_Mask = Input%Invalid_Data_Mask
  Temperature_Cirrus = MISSING_VALUE_REAL4

  !--- call LRC routine
  if (USE_LRC_FLAG) then

    if (associated(Input%Elem_Idx_LRC_Input) .and. &
        associated(Input%Line_Idx_LRC_Input)) then

      Elem_Idx_LRC = Input%Elem_Idx_LRC_Input
      Line_Idx_LRC = Input%Line_Idx_LRC_Input

    else

      where(Input%Cloud_Mask == Symbol%CLEAR .or. Input%Cloud_Mask == Symbol%PROB_CLEAR)
          Skip_LRC_Mask = Symbol%YES
      endwhere

      !--- Min_Bt_110um_Lrc and Max_Bt_110um_Lrc are parameters from
      !--- acha_parameters.inc.
      call LOCAL_LINEAR_RADIATIVE_CENTER(Symbol%YES,Symbol%NO,&
                                         LRC_Meander_Flag, &
                                         Input%Bt_110um, &
                                         Element_Idx_Min, &
                                         Input%Number_of_Elements, & 
                                         Line_Idx_Min,  &
                                         Input%Number_of_Lines, & 
                                         Max_LRC_Distance,  &
                                         Min_LRC_Jump,  &
                                         Max_LRC_Jump,  &
                                         Grad_Flag_LRC,  &
                                         MISSING_VALUE_integer4, &
                                         Skip_LRC_Mask, &
                                         Min_Bt_110um_Lrc,  &
                                         Max_Bt_110um_Lrc, &
                                         Elem_Idx_LRC,  &
                                         Line_Idx_LRC)

    endif
  endif

  !--------------------------------------------------------------------------
  ! Multi-Layer Logic Implemented via cloud type
  !-------------------------------------------------------------------------
  Cloud_Type_Temp = Input%Cloud_Type
 
  if (MULTI_LAYER_LOGIC_FLAG == 1) then 
   where(Input%Cloud_Type == Symbol%OVERLAP_TYPE)
     Cloud_Type_Temp = Symbol%CIRRUS_TYPE
   endwhere
  endif

  if (MULTI_LAYER_LOGIC_FLAG == 2) then 
   where(Input%Cloud_Type == Symbol%CIRRUS_TYPE)
     Cloud_Type_Temp = Symbol%OVERLAP_TYPE
   endwhere
  endif

  !--------------------------------------------------------------------------
  ! For Testing, allow a cloud to be specified (from acha_parameters.inc)
  !--------------------------------------------------------------------------
  if (Cloud_Type_Forced /= -1) then 
     Cloud_Type_Temp = Cloud_Type_Forced
  endif

   !--------------------------------------------------------------------------
   ! loop over pixels in scanlines
   !--------------------------------------------------------------------------

   !--- set the single pixel dump indices to missing - will set below
   Elem_Abs_Idx_Dump = -1
   Line_Abs_Idx_Dump = -1
   if (present(Dump)) then
      Elem_Abs_Idx_Dump = Dump%Elem_Abs_Idx
      Line_Abs_Idx_Dump = Dump%Line_Abs_Idx
   endif

   Line_loop: do Line_Idx = Line_Idx_min,Input%Number_of_Lines + Line_Idx_min - 1

    Element_Loop:   do Elem_Idx = 1, Input%Number_of_Elements

    !--- compute absolute element and line indices (for potential dumping)
    Elem_Abs_Idx = -1
    Line_Abs_Idx = -1
    if (present(Dump)) then
      Elem_Abs_Idx = Elem_Idx
      Line_Abs_Idx = Line_Idx + (Dump%Segment_Number-1)*Dump%Number_Lines_Per_Segment
    endif

    !--- check if this pixel is chosen for diagnostic dump
    Dump_Diag = .false.
    if (.not. Dump_Diag) then 
      if (Elem_Abs_Idx_Dump > 0 .and. Line_Abs_Idx_Dump > 0) then 
        if (Elem_Abs_Idx == Elem_Abs_Idx_Dump .and. Line_Abs_Idx ==  Line_Abs_Idx_Dump) then
          Dump_Diag = .true.
          if (Lun_Prof_Dump < 0 .and. Lun_Iter_Dump < 0) then
            Lun_Prof_Dump = GET_LUN_ACHA()
            open(unit=Lun_Prof_Dump,file=trim(File_Name_Prof_Dump),form="formatted",status='unknown',action='write')
            Lun_Iter_Dump = GET_LUN_ACHA()
            open(unit=Lun_Iter_Dump,file=trim(File_Name_Iter_Dump),form="formatted",status='unknown',action='write')
          endif
        endif
      endif
    endif

    !---- null profile pointers each time 
    call NULL_PIX_POINTERS(Input, ACHA_RTM_NWP)

    !---------------------------------------------------------------
    ! Check to see if this pixel should be skipped
    !---------------------------------------------------------------
    call TEST_INPUT(Acha_Mode_Flag, ABi_Use_104um_Flag, Input, Elem_Idx,Line_Idx, Bad_Input_Flag)

    !--- if a bad pixel encountered, take action
    if (Bad_Input_Flag) then 
          Output%Packed_Qf(Elem_Idx,Line_Idx) =  0_int1
          Output%Packed_Meta_Data(Elem_Idx,Line_Idx) =  0_int1
          Output%Qf(Elem_Idx,Line_Idx) = CTH_DQF_BAD_RETREVIAL
          cycle
    endif

    !--- for convenience, save nwp indices to local variables
    Ivza =  Input%Viewing_Zenith_Angle_Idx_RTM(Elem_Idx,Line_Idx)
    ilrc = Elem_Idx_LRC(Elem_Idx,Line_Idx)
    jlrc = Line_Idx_LRC(Elem_Idx,Line_Idx)
    Cloud_Type = Cloud_Type_Temp(Elem_Idx,Line_Idx)

    T_Tropo = Input%Tropopause_Temperature(Elem_Idx,Line_Idx)
    Z_Tropo = Input%Tropopause_Height(Elem_Idx,Line_Idx)
    P_Tropo = Input%Tropopause_Pressure(Elem_Idx,Line_Idx)
    
    !--- Qc indices
    if (Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx) <= 0 .or. &
        Input%Line_Idx_Nwp(Elem_Idx,Line_Idx) <= 0 .or. &
        Input%Viewing_Zenith_Angle_Idx_RTM(Elem_Idx,Line_Idx) <= 0) then 
          Output%Packed_Qf(Elem_Idx,Line_Idx) =  0_int1
          Output%Packed_Meta_Data(Elem_Idx,Line_Idx) =  0_int1
          Output%Qf(Elem_Idx,Line_Idx) = CTH_DQF_BAD_RETREVIAL
         cycle 
    endif

    !--------------------------------------------------------------------
    ! get profiles for this pixel
    !--------------------------------------------------------------------
    call ACHA_FETCH_PIXEL_NWP_RTM(Input, Symbol, &
                                 Elem_Idx,Line_Idx, &
                                 ACHA_RTM_NWP,RTM_NWP_Error_Flag)
    if (RTM_NWP_Error_Flag /= 0) then 
        Output%Packed_Qf(Elem_Idx,Line_Idx) =  0_int1
        Output%Packed_Meta_Data(Elem_Idx,Line_Idx) =  0_int1
        Output%Qf(Elem_Idx,Line_Idx) = CTH_DQF_BAD_RETREVIAL
        cycle
    endif

    !-----------------------------------------------------------------------
    ! include code to setup local profiles correctly 
    !-----------------------------------------------------------------------
    call SETUP_REFERENCE_CHANNEL_PROFILES(ABI_Use_104um_Flag,ACHA_RTM_NWP)

    !---- output to profile dump
    if (Dump_Diag) then
       do Lev_Idx = 1,Num_Levels_Rtm_Prof
             write(unit=Lun_Prof_Dump,fmt="(I3,F8.2,F8.1,F8.3,6F8.2)") Lev_Idx,  &
                   ACHA_RTM_NWP%P_Prof(Lev_Idx),  &
                   ACHA_RTM_NWP%Z_Prof(Lev_Idx), &
                   ACHA_RTM_NWP%T_Prof(Lev_Idx), &
                   ACHA_RTM_NWP%Atm_Rad_Prof_110um(Lev_Idx), &
                   ACHA_RTM_NWP%Atm_Trans_Prof_110um(Lev_Idx), &
                   ACHA_RTM_NWP%Black_Body_Rad_Prof_110um(Lev_Idx)
       enddo 
    endif

    !---- output to retrieval dump
    if (Dump_Diag) then
       print *, 'Writing diag output to file '
       write(unit=Lun_Iter_Dump,fmt=*) "========================================================"
       write(unit=Lun_Iter_Dump,fmt=*) "Diagnostic Dump"
       write(unit=Lun_Iter_Dump,fmt=*) "========================================================"
       write(unit=Lun_Iter_Dump,fmt=*) "Element, Line Indices (Relative to Segment) = ", Elem_Idx,Line_Idx
       write(unit=Lun_Iter_Dump,fmt=*) "Surface Elevation = ", Input%Surface_Elevation(Elem_Idx,Line_Idx)
       write(unit=Lun_Iter_Dump,fmt=*) "Latitude = ", Input%Latitude(Elem_Idx,Line_Idx)
       write(unit=Lun_Iter_Dump,fmt=*) "Longitude = ", Input%Longitude(Elem_Idx,Line_Idx)
       write(unit=Lun_Iter_Dump,fmt=*) "Zenith Angle = ", Input%Sensor_Zenith_Angle(Elem_Idx,Line_Idx)
       write(unit=Lun_Iter_Dump,fmt=*) "Input Type = ", Input%Cloud_Type(Elem_Idx,Line_Idx)
       write(unit=Lun_Iter_Dump,fmt=*) "Input Ice Prob = ", Input%Ice_Cloud_Probability(Elem_Idx,Line_Idx)
    endif

  !-------------------------------------------------------------------
  !Apply Opaque Retrieval for Acha_Mode_Flag = 1, then cycle
  !-------------------------------------------------------------------
  if (trim(Acha_Mode_Flag) == '110') then
     call APPLY_OPAQUE_RETRIEVAL(Input,Symbol,Elem_Idx,Line_Idx,Output)
     Output%Qf(Elem_Idx,Line_Idx) = CTH_DQF_OPAQUE_RETREVIAL
     cycle
  endif 

   !----------------------------------------------------------------------
   ! determine cloud phase from cloud type for convenience
   !----------------------------------------------------------------------
   call COMPUTE_PHASE(Symbol,Input%Ice_Prob_Ap(Elem_Idx,Line_Idx), Input%Ice_Prob_Ap_Uncer(Elem_Idx,Line_Idx), &
                       Dump_Diag, Lun_Iter_Dump, Cloud_Phase)

!  if (CONSTRAIN_ICE_PROB) then
!     Ice_Probability_Ap_Uncertainty = 0.01
!  endif

   !----------------------------------------------------------------------
   !--- Set Meta Data Flags
   !----------------------------------------------------------------------
   call COMPUTE_META_DATA(Cloud_Phase, USE_LRC_FLAG, Cloud_Type, Meta_Data_Flags)

   !-----------------------------------------------------------------------
   ! assign values to y and y_variance
   !----------------------------------------------------------------------
   !--- For GOES-17 mitigation, COMPUTE_Y may have 11 um data switched with 10.4
   !--- um data.
   call COMPUTE_Y(Acha_Mode_Flag,Input,Element_Idx_Min, Line_Idx_Min, Elem_Idx,Line_Idx, &
                  y,y_variance,Dump_Diag, Lun_Iter_Dump)   

   !-------------------------------------------------------------------
   ! Determine surface type for use in forward model
   ! 0 = Water
   ! 1 = Land
   ! 2 = Snow
   ! 3 = Desert
   ! 4 = Arctic
   ! 5 = Antarctic
   !-------------------------------------------------------------------
   call DETERMINE_SFC_TYPE_FORWARD_MODEL(Input%Surface_Type(Elem_Idx,Line_Idx), &
                                         Input%Snow_Class (Elem_Idx,Line_Idx), &
                                         Input%Latitude(Elem_Idx,Line_Idx), &
                                         Input%Surface_Emissivity_038um(Elem_Idx,Line_Idx), &
                                         Sfc_Type_Forward_Model)

   !-------------------------------------------------------------------
   ! Based on fm surface type, set the clear-sky covariance terms
   !-------------------------------------------------------------------
   call SET_CLEAR_SKY_COVARIANCE_TERMS(Sfc_Type_Forward_Model)

   !--------------------------------------------------------------------
   ! pick a priori conditions
   !--------------------------------------------------------------------

   !--- logic for unmasked or untyped pixels (Output%Ec)
   if (Input%Process_Undetected_Cloud_Flag == Symbol%YES) then
         if (Input%Tc_Opaque(Elem_Idx,Line_Idx) < 260.0 .and.  &
             Input%Tc_Opaque(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) then
             Cloud_Type = Symbol%CIRRUS_TYPE
         else
             Cloud_Type = Symbol%FOG_TYPE
         endif
   endif

  !--- For bad GOES-17 data, all inputs have been switched from 11um to 104um,
  !--- if needed. Variables remain the same name.
  !---- Compute 110um emissivity referenced to tropopause
  Emiss_110um_Tropo = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( &
                             ACHA_RTM_NWP%Tropo_Level, &
                             Input%Rad_110um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_110um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_110um)

  !------------------------------------------------------------------------
  !  lower cloud (surface) a prior values
  !------------------------------------------------------------------------
  Tsfc_Est = Input%Surface_Temperature(Elem_Idx,Line_Idx)

  Ts_Ap = Tsfc_Est
  Ts_Ap_Uncer = Ts_Ap_Uncer_Sfc
  if (Cloud_Type == Symbol%OVERLAP_TYPE) then
    Ts_Ap_Uncer = Ts_Ap_Uncer_Lower_Cld
  endif

   !-----------------------------------------------------------------------
   ! For bad GOES-17 11 um data, and good 10.4 um data,
   ! switch reference channel microphysics if necessary.
   !-----------------------------------------------------------------------
   if (ABI_Use_104um_Flag) then
      !--- water 
      Beta_110um_142um_Coef_Water = Beta_104um_142um_Coef_Water
      Beta_110um_139um_Coef_Water = Beta_104um_139um_Coef_Water
      Beta_110um_136um_Coef_Water = Beta_104um_136um_Coef_Water
      Beta_110um_133um_Coef_Water = Beta_104um_133um_Coef_Water
      Beta_110um_085um_Coef_Water = Beta_104um_085um_Coef_Water
      Beta_110um_097um_Coef_Water = Beta_104um_097um_Coef_Water
      Beta_110um_073um_Coef_Water = Beta_104um_073um_Coef_Water
      Beta_110um_067um_Coef_Water = Beta_104um_067um_Coef_Water
      Beta_110um_062um_Coef_Water = Beta_104um_062um_Coef_Water
      Beta_110um_038um_Coef_Water = Beta_104um_038um_Coef_Water

      !--- ice
      Beta_110um_142um_Coef_Ice = Beta_104um_142um_Coef_Ice
      Beta_110um_139um_Coef_Ice = Beta_104um_139um_Coef_Ice
      Beta_110um_136um_Coef_Ice = Beta_104um_136um_Coef_Ice
      Beta_110um_133um_Coef_Ice = Beta_104um_133um_Coef_Ice
      Beta_110um_097um_Coef_Ice = Beta_104um_097um_Coef_Ice
      Beta_110um_085um_Coef_Ice = Beta_104um_085um_Coef_Ice
      Beta_110um_073um_Coef_Ice = Beta_104um_073um_Coef_Ice
      Beta_110um_067um_Coef_Ice = Beta_104um_067um_Coef_Ice
      Beta_110um_062um_Coef_Ice = Beta_104um_062um_Coef_Ice
      Beta_110um_038um_Coef_Ice = Beta_104um_038um_Coef_Ice
   endif

  !----------------------------------------------------------------------------
  ! Compute Xa and Sa
  !----------------------------------------------------------------------------
  call COMPUTE_SA(Input,Symbol,Num_Param, &
                  Elem_Idx,Line_Idx,x_Ap,Sa,Sa_Inv, &
                  Singular_Warning_First_Time, &
                  Dump_Diag, Lun_Iter_Dump,Prior_Success_Flag,Fail_Flag(Elem_Idx,Line_Idx))

  if (Prior_Success_Flag .eqv. .false.) then 
     cycle
  endif
  if (Fail_Flag(Elem_Idx,Line_Idx) == symbol%YES) then 
     cycle
  endif

  !--------------------------------------------------
  ! assign surface emissivity for non-overlap type
  !--------------------------------------------------
  call SET_SURFACE_EMISSIVITY(Input, Cloud_Type, Elem_Idx,Line_Idx, &
                              Emiss_Sfc_038um, Emiss_Sfc_062um,&
                              Emiss_Sfc_067um, Emiss_Sfc_073um, Emiss_Sfc_085um, &
                              Emiss_Sfc_097um, Emiss_Sfc_104um, Emiss_Sfc_110um, &
                              Emiss_Sfc_120um, Emiss_Sfc_133um, Emiss_Sfc_136um, &
                              Emiss_Sfc_139um, Emiss_Sfc_142um)

 !----------------------------------------------------------------
 ! Determine the level of the highest inversion (0=if none)
 !----------------------------------------------------------------
 call DETERMINE_INVERSION_CHARACTERISTICS(ACHA_RTM_NWP, &
                                          Symbol%YES, &
                                          Symbol%NO,  &
                                          Input%Surface_Air_Temperature(Elem_Idx,Line_Idx),&
                                          Input%Surface_Elevation(Elem_Idx,Line_Idx),      &
                                          Inver_Top_Level_RTM,    &
                                          Inver_Base_Level_RTM,   & 
                                          Inver_Top_Height,       &
                                          Inver_Base_Height,      & 
                                          Inver_Strength, &
                                          Dump_Diag,  &
                                          Lun_Iter_Dump)
!Diag%Array_1(Elem_Idx,Line_Idx) = Inver_Top_Height
!Diag%Array_2(Elem_Idx,Line_Idx) = Inver_Base_Height
!Diag%Array_3(Elem_Idx,Line_Idx) = Inver_Strength

!-----------------------------------------------------------------
! start of retrieval loop
!-----------------------------------------------------------------
Iter_Idx = 0
Converged_Flag(Elem_Idx,Line_Idx) = Symbol%NO
Fail_Flag(Elem_Idx,Line_Idx) = Symbol%NO
Fail_Flag_Simple(Elem_Idx,Line_Idx) = Symbol%NO
Fail_Flag_Full(Elem_Idx,Line_Idx) = Symbol%NO

!-----------------------------------------------------------------
! Perform Retrieval
!-----------------------------------------------------------------
 if (Dump_Diag) then
     write(unit=Lun_Iter_Dump,fmt=*) "Tsfc Estimate = ", Tsfc_Est
     write(unit=Lun_Iter_Dump,fmt=*) "T Tropopause = ", T_Tropo
     write(unit=Lun_Iter_Dump,fmt=*) "Z Tropopause = ", Z_Tropo
     write(unit=Lun_Iter_Dump,fmt=*) "P Tropopause = ", P_Tropo
 endif

 !------------------------------------------------------------------------------------------
 !  Call 3-parameter (simple) retrieval
 !------------------------------------------------------------------------------------------
 x_Ap_Simple = x_Ap(1:3)
 Sa_Simple = Sa(1:3,1:3)
 Singular_Flag =  INVERT_MATRIX(Sa_Simple, Sa_Inv_Simple, Num_Param_Simple)
 if (Singular_Flag == 1) print *, "Cloud Height warning ==> Singular Sa Simple in ACHA", Sa(1,1),Sa(2,2),Sa(3,3)


if (Dump_Diag) then
     write(unit=Lun_Iter_Dump,fmt=*) "========================================================"
     write(unit=Lun_Iter_Dump,fmt=*) "Call to 3-Parameter Retrieval"
     write(unit=Lun_Iter_Dump,fmt=*) "========================================================"
endif

    call FULL_ACHA_RETRIEVAL(&
         Input,Acha_Mode_Flag, Symbol, &
         Num_Obs,Num_Param_Simple,y,y_variance,f, &
         x_Ap_Simple,Sa_Inv_Simple,x_Simple,Sx_Simple,AKM_Simple,&
         Output%Conv_Test(Elem_Idx,Line_Idx),Output%Cost(Elem_Idx,Line_Idx), &
         Output%Goodness(Elem_Idx,Line_Idx),Convergence_Criteria_Simple, &
         ACHA_RTM_NWP%Z_Prof,Tsfc_Est,T_Tropo,Z_Tropo,P_Tropo, Cloud_Type, &
         Input%Cosine_Zenith_Angle(Elem_Idx,Line_Idx), &
         Output%Zc_Base(Elem_Idx,Line_Idx), &
         ACHA_RTM_NWP, &
         Beta_110um_142um_Coef_Water, &
         Beta_110um_139um_Coef_Water, &
         Beta_110um_136um_Coef_Water, &
         Beta_110um_133um_Coef_Water, &
         Beta_110um_104um_Coef_Water, &
         Beta_110um_097um_Coef_Water, &
         Beta_110um_085um_Coef_Water, &
         Beta_110um_073um_Coef_Water, &
         Beta_110um_067um_Coef_Water, &
         Beta_110um_062um_Coef_Water, &
         Beta_110um_038um_Coef_Water, &
         Beta_110um_142um_Coef_Ice, &
         Beta_110um_139um_Coef_Ice, &
         Beta_110um_136um_Coef_Ice, &
         Beta_110um_133um_Coef_Ice, &
         Beta_110um_104um_Coef_Ice, &
         Beta_110um_097um_Coef_Ice, &
         Beta_110um_085um_Coef_Ice, &
         Beta_110um_073um_Coef_Ice, &
         Beta_110um_067um_Coef_Ice, &
         Beta_110um_062um_Coef_Ice, &
         Beta_110um_038um_Coef_Ice, &
         Emiss_Sfc_038um, &
         Emiss_Sfc_062um, &
         Emiss_Sfc_067um, &
         Emiss_Sfc_073um, &
         Emiss_Sfc_085um, &
         Emiss_Sfc_097um, &
         Emiss_Sfc_104um, &
         Emiss_Sfc_110um, &
         Emiss_Sfc_120um, &
         Emiss_Sfc_133um, &
         Emiss_Sfc_136um, &
         Emiss_Sfc_139um, &
         Emiss_Sfc_142um, &
         Output%Tc_Eff(Elem_Idx,Line_Idx), &
         Converged_Flag(Elem_Idx,Line_Idx), &
         Fail_Flag_Simple(Elem_Idx,Line_Idx), &
         Dump_Diag, &
         Lun_Iter_Dump)

         !--- put simple into full arrays
         x = MISSING_VALUE_REAL4
         Sx = MISSING_VALUE_REAL4
         AKM = MISSING_VALUE_REAL4
         x(1:3) = x_Simple
         x(4) = x_ap(4)
         x(5) = x_ap(5)
         Sx(1:3,1:3) = Sx_Simple
         Sx(4,4) = Sa(4,4)
         Sx(5,5) = Sa(5,5)
         AKM(1:3,1:3) = AKM_Simple
         AKM(4,4) = 0.0
         AKM(5,5) = 0.0
         Fail_Flag(Elem_Idx,Line_Idx) = Fail_Flag_Simple(Elem_Idx,Line_Idx)

         !---- if a full retrieval is wanted but simple retrieval failed, skip
         if (FULL_RETRIEVAL .and. Fail_Flag_Simple(Elem_Idx,Line_Idx) == symbol%YES) then
           cycle
         endif


         !---- if a full retrieval is wanted but simple retrieval worked, continue
if (FULL_RETRIEVAL .and. Fail_Flag_Simple(Elem_Idx,Line_Idx) == symbol%NO) then

! MULTI_LAYER_LOGIC_FLAG
! 0 - (baseline) just use the multilayer id in cloud type
! 1 - treat all multilayer like cirrus
! 2 - assume all cirrus are multilayer and let acha decide

    if (MULTI_LAYER_LOGIC_FLAG == 1) then
       if (Cloud_Type == Symbol%Overlap_Type) Cloud_Type = Symbol%Cirrus_Type
    endif
    if (MULTI_LAYER_LOGIC_FLAG == 2) then
       if (Cloud_Type == Symbol%Cirrus_Type) Cloud_Type = Symbol%Overlap_Type
    endif

    if (Cloud_Type == symbol%Overlap_Type .or. Input%Cloud_Phase_Uncertainty(Elem_Idx,Line_Idx) >= 0.10)  then
    !if (Cloud_Type == symbol%Overlap_Type) then

     !---  use simple retrieval as first guess for full retrieval
     if (x(1) /= MISSING_VALUE_REAL4) x_ap(1) = x(1)
     if (x(2) /= MISSING_VALUE_REAL4) x_ap(2) = x(2)
     if (x(3) /= MISSING_VALUE_REAL4) x_ap(3) = x(3)
     if (Sx(1,1) /= MISSING_VALUE_REAL4) Sa(1,1) = max(Tc_Ap_Uncer_Min**2,1.0*Sx(1,1))
     if (Sx(2,2) /= MISSING_VALUE_REAL4) Sa(2,2) = max(Ec_Ap_Uncer_Min**2,1.0*Sx(2,2))
     if (Sx(3,3) /= MISSING_VALUE_REAL4) Sa(3,3) = max(Beta_Ap_Uncer_Min**2,1.0*Sx(3,3))
     Singular_Flag =  INVERT_MATRIX(Sa, Sa_Inv, Num_Param)
     if (Singular_Flag == 1) print *, "Cloud Height warning ==> Singular Sa in ACHA", Sa(1,1),Sa(2,2),Sa(3,3),Sa(4,4),Sa(5,5)

     if (Dump_Diag) then
        write(unit=Lun_Iter_Dump,fmt=*) "========================================================"
        write(unit=Lun_Iter_Dump,fmt=*) "Call to 5-Parameter Retrieval"
        write(unit=Lun_Iter_Dump,fmt=*) "========================================================"
     endif

     !--- call retrieval
     call FULL_ACHA_RETRIEVAL(&
          Input,Acha_Mode_Flag, Symbol, &
          Num_Obs,Num_Param,y,y_variance,f,x_Ap,Sa_Inv,x,Sx,AKM,&
          Output%Conv_Test(Elem_Idx,Line_Idx),Output%Cost(Elem_Idx,Line_Idx), &
          Output%Goodness(Elem_Idx,Line_Idx),Convergence_Criteria, &
          ACHA_RTM_NWP%Z_Prof,Tsfc_Est,T_Tropo,Z_Tropo,P_Tropo, Cloud_Type, &
          Input%Cosine_Zenith_Angle(Elem_Idx,Line_Idx), &
          Output%Zc_Base(Elem_Idx,Line_Idx), &
          ACHA_RTM_NWP, &
          Beta_110um_142um_Coef_Water, &
          Beta_110um_139um_Coef_Water, &
          Beta_110um_136um_Coef_Water, &
          Beta_110um_133um_Coef_Water, &
          Beta_110um_104um_Coef_Water, &
          Beta_110um_097um_Coef_Water, &
          Beta_110um_085um_Coef_Water, &
          Beta_110um_073um_Coef_Water, &
          Beta_110um_067um_Coef_Water, &
          Beta_110um_062um_Coef_Water, &
          Beta_110um_038um_Coef_Water, &
          Beta_110um_142um_Coef_Ice, &
          Beta_110um_139um_Coef_Ice, &
          Beta_110um_136um_Coef_Ice, &
          Beta_110um_133um_Coef_Ice, &
          Beta_110um_104um_Coef_Ice, &
          Beta_110um_097um_Coef_Ice, &
          Beta_110um_085um_Coef_Ice, &
          Beta_110um_073um_Coef_Ice, &
          Beta_110um_067um_Coef_Ice, &
          Beta_110um_062um_Coef_Ice, &
          Beta_110um_038um_Coef_Ice, &
          Emiss_Sfc_038um, &
          Emiss_Sfc_062um, &
          Emiss_Sfc_067um, &
          Emiss_Sfc_073um, &
          Emiss_Sfc_085um, &
          Emiss_Sfc_097um, &
          Emiss_Sfc_104um, &
          Emiss_Sfc_110um, &
          Emiss_Sfc_120um, &
          Emiss_Sfc_133um, &
          Emiss_Sfc_136um, &
          Emiss_Sfc_139um, &
          Emiss_Sfc_142um, &
          Output%Tc_Eff(Elem_Idx,Line_Idx), &
          Converged_Flag(Elem_Idx,Line_Idx), &
          Fail_Flag_Full(Elem_Idx,Line_Idx), &
          Dump_Diag, &
          Lun_Iter_Dump)
    endif
    Fail_Flag(Elem_Idx,Line_Idx) = Fail_Flag_Full(Elem_Idx,Line_Idx)

 endif
 
 if (Dump_Diag) then
     write(unit=Lun_Iter_Dump,fmt=*) "========================================================"
     write(unit=Lun_Iter_Dump,fmt=*) "Returned from ACHA_RETRIEVAL"
     write(unit=Lun_Iter_Dump,fmt=*) "========================================================"
     write(unit=Lun_Iter_Dump,fmt=*) "final f = ", f
     write(unit=Lun_Iter_Dump,fmt=*) "      y = ", y
     write(unit=Lun_Iter_Dump,fmt=*) "final x = ", x
     write(unit=Lun_Iter_Dump,fmt=*) "   x_ap = ", x_ap
     write(unit=Lun_Iter_Dump,fmt=*) "Converged_Flag = ", Converged_Flag(Elem_Idx,Line_Idx)
     write(unit=Lun_Iter_Dump,fmt=*) "Fail_Flag = ", Fail_Flag(Elem_Idx,Line_Idx)
 endif

  !--- Save the OE to the Output Structure
  call SAVE_X_2_OUTPUT(Elem_Idx,Line_Idx,Symbol,Cloud_Type,Fail_Flag(Elem_Idx,Line_Idx), &
                  x,x_ap,Sa,Sx,AKM,Meta_Data_Flags,Output)

! if (Fail_Flag(Elem_Idx,Line_Idx) == symbol%YES) then 
!         print *, Fail_Flag_Simple(ELem_Idx,Line_Idx), Fail_Flag_Full(Elem_Idx,Line_Idx)
!         write(unit=6,fmt="(20F6.1)") x, y, f-y
!         print *, Output%Conv_Test(Elem_Idx,Line_Idx),Output%Cost(Elem_Idx,Line_Idx),&
!         Output%Goodness(Elem_Idx,Line_Idx),Convergence_Criteria
!         print *
!         stop
! endif

  !--- null profile pointers each time 
  call NULL_PIX_POINTERS(Input, ACHA_RTM_NWP)

  !--- set output packed quality flags
  call SET_OUTPUT_PACKED_QF(Output,Elem_Idx,Line_Idx)

 end do Element_Loop
end do Line_Loop

!---------------------------------------------------------------------------
! Post Retrieval Processing
!---------------------------------------------------------------------------

!-- perform conversion of Tc to Zc and Pc
call CONVERT_TC_TO_PC_AND_ZC(Input,Symbol,Cloud_Type_Temp,Output)

!--- compute cloud emissivity in each channel used
call COMPUTE_SPECTRAL_CLOUD_EMISSIVITY(Input,Symbol,Output)

!--- Determine Cloud Type
call DETERMINE_ACHA_CLOUD_TYPE(Input,Fail_Flag,Symbol,Output)

!--- Apply Parallax Correction 
call PARALLAX_ACHA(Output%Zc, Input%Surface_Elevation, &
                     Input%Latitude, Input%Longitude, &
                     Input%Sensor_Zenith_Angle, &
                     Input%Sensor_Azimuth_Angle, &
                     Output%Latitude_Pc,&
                     Output%Longitude_Pc) 

!-----------------------------------------------------------------------
!--- clean-up and prepare for exit
!-----------------------------------------------------------------------

!--- deallocate 2D arrays
if (allocated(Elem_Idx_LRC)) deallocate(Elem_Idx_LRC)
if (allocated(Line_Idx_LRC)) deallocate(Line_Idx_LRC)
if (allocated(Skip_LRC_Mask)) deallocate(Skip_LRC_Mask)
if (allocated(Temperature_Cirrus)) deallocate(Temperature_Cirrus)
if (allocated(Fail_Flag)) deallocate(Fail_Flag)
if (allocated(Fail_Flag_Simple)) deallocate(Fail_Flag_Simple)
if (allocated(Fail_Flag_Full)) deallocate(Fail_Flag_Full)
if (allocated(Converged_Flag)) deallocate(Converged_Flag)
if (allocated(Chan_Idx_y)) deallocate(Chan_Idx_y)
if (allocated(Cloud_Type_Temp)) deallocate(Cloud_Type_Temp)

!--- deallocate 1D-VAR arrays
deallocate(y)
deallocate(y_variance)
deallocate(f)

deallocate(x)
deallocate(x_Ap)
deallocate(Sa)
deallocate(Sa_Inv)
deallocate(Sx)
deallocate(x_Simple)
deallocate(x_Ap_Simple)
deallocate(Sa_Simple)
deallocate(Sa_Inv_Simple)
deallocate(Sx_Simple)

!---close single pixel dump output
if (Lun_Prof_Dump > 0) close(unit=Lun_Prof_Dump)
if (Lun_Iter_Dump > 0) close(unit=Lun_Iter_Dump)

end subroutine  AWG_CLOUD_HEIGHT_ALGORITHM

!-------------------------------------------------------------------
! Determine surface type for use in forward model
! 0 = Water
! 1 = Land
! 2 = Snow
! 3 = Desert
! 4 = Arctic
! 5 = Antarctic
!-------------------------------------------------------------------
subroutine DETERMINE_SFC_TYPE_FORWARD_MODEL( &
                                         Surface_Type, &
                                         Snow_Class, &
                                         Latitude, &
                                         Ch20_Surface_Emissivity, &
                                         Sfc_Type_Forward_Model)

integer(kind=int1), intent(in):: Surface_Type
integer(kind=int1), intent(in):: Snow_Class
real(kind=real4), intent(in):: Latitude
real(kind=real4), intent(in):: Ch20_Surface_Emissivity
integer(kind=int4), intent(out):: Sfc_Type_Forward_Model


if (Surface_Type == Symbol%WATER_SFC) then
        Sfc_Type_Forward_Model = 0
else
        Sfc_Type_Forward_Model = 1   !Land
endif
if (Snow_Class == Symbol%SNOW .and. &
    Latitude > -60.0) then
        Sfc_Type_Forward_Model = 2   !Snow
endif
if (Surface_Type /= Symbol%WATER_SFC .and. &
    Snow_Class == Symbol%NO_SNOW .and.  &
    Ch20_Surface_Emissivity > 0.90 .and.  &
    abs(Latitude) < 60.0) then
        Sfc_Type_Forward_Model = 3   !Desert
endif
if (Snow_Class == Symbol%SEA_ICE .and. &
    Latitude > 60.0) then
        Sfc_Type_Forward_Model = 4   !Arctic
endif
if (Snow_Class /= Symbol%NO_SNOW .and. Latitude < -60.0) then
        Sfc_Type_Forward_Model = 5   !Antartica
endif

end subroutine DETERMINE_SFC_TYPE_FORWARD_MODEL

!----------------------------------------------------------------------
!  In GEOCAT, Acha_Mode_Flag can not be passed in, use this routine to 
!  determine it in based on the available channels
!----------------------------------------------------------------------
subroutine DETERMINE_ACHA_MODE_BASED_ON_CHANNELS( &
                                                 Acha_Mode_Flag, &
                                                 Chan_On_038um, &
                                                 Chan_On_062um, &
                                                 Chan_On_067um, &
                                                 Chan_On_085um, &
                                                 Chan_On_110um, &
                                                 Chan_On_120um, &
                                                 Chan_On_133um, &
                                                 Chan_On_136um, &
                                                 Chan_On_139um, &
                                                 Chan_On_142um)

  character(len=*), intent(inout):: Acha_Mode_Flag
  integer, intent(in):: Chan_On_038um
  integer, intent(in):: Chan_On_062um
  integer, intent(in):: Chan_On_067um
  integer, intent(in):: Chan_On_085um
  integer, intent(in):: Chan_On_110um
  integer, intent(in):: Chan_On_120um
  integer, intent(in):: Chan_On_133um
  integer, intent(in):: Chan_On_136um
  integer, intent(in):: Chan_On_139um
  integer, intent(in):: Chan_On_142um

  if (trim(Acha_Mode_Flag) == 'unknown') then
     if (Chan_On_110um == Symbol%YES .and. Chan_On_120um == Symbol%YES) then
        if (Chan_On_133um == Symbol%YES) then
            Acha_Mode_Flag = "110_120_133"          ! 11/12/13.3 um
        elseif (Chan_On_085um == Symbol%YES) then
            Acha_Mode_Flag = "085_110_120"          ! 8.5/11/12 um
        elseif (Chan_On_067um == Symbol%YES) then
            Acha_Mode_Flag = "067_110_120"
        else
            Acha_Mode_Flag = "110_120"
        endif
     endif
  endif

  if (trim(Acha_Mode_Flag) == "unknown") then
     if (Chan_On_120um == Symbol%NO) then
        if (Chan_On_067um == Symbol%YES .and. Chan_On_133um == Symbol%YES) then
              Acha_Mode_Flag = "067_110_133"
        endif
        if (Chan_On_067um == Symbol%NO .and. Chan_On_133um == Symbol%YES) then
              Acha_Mode_Flag = "110_133"
        endif
        if (Chan_On_067um == Symbol%YES .and. Chan_On_133um == Symbol%NO) then
              Acha_Mode_Flag = "067_110"
        endif
     endif
  endif

  if (trim(Acha_Mode_Flag) == "unknown") then
        if (Chan_On_038um == Symbol%YES) then
              Acha_Mode_Flag = '038_110'
        endif
  endif

  !--- if unsuccessful, resort to mode 1
  if (trim(Acha_Mode_Flag) == "unknown") then
     if (Chan_On_110um == Symbol%YES) then
              Acha_Mode_Flag = '110'
     endif
  endif
   
   
end subroutine  DETERMINE_ACHA_MODE_BASED_ON_CHANNELS

 !-------------------------------------------------------------------------
 ! Input:
 ! Tropo_Level - level in RTM profiles of the Tropopause
 ! Sfc_Level - level in RTM profiles closest but above the surface
 ! Sfc_Air_Temp = air temperature at the surface level
 ! Sfc_Height = height of the surface level (m)
 ! Inversion_Top_Level -  level in RTM profiles closest to but below the top of
 !                         inversion
 ! Inversion_Base_Level -  level in RTM profiles closest to but above the base of
 !                         inversion
 ! Inversion_Strength - Temperature difference between Top and Base (K)
 ! Inversion_Base_Height -  Height of Inversion Base (m)
 ! Inversion_Top_Height -  Height of Inversion Top (m)
 !
 ! Input - via module-wide variables
 !
 ! Output - via module-wide variables
 ! Inver_Prof_RTM - level flags (0/1) if inversion present
 !--------------------------------------------------------------------------
 subroutine DETERMINE_INVERSION_CHARACTERISTICS(ACHA_RTM_NWP, &
                                                Symbol_Yes,               &
                                                Symbol_No,                &
                                                Sfc_Air_Temp,             &
                                                Sfc_Height,               &
                                                Top_Lev_Idx,            &
                                                Base_Lev_Idx,           & 
                                                Inversion_Top_Height,     &
                                                Inversion_Base_Height,    & 
                                                Inversion_Strength, &
                                                Dump_Diag, &
                                                Lun_Iter_Dump)

   type(acha_rtm_nwp_struct), intent(inout) :: ACHA_RTM_NWP
   integer(kind=int1), intent(in) :: Symbol_yes
   integer(kind=int1), intent(in) :: Symbol_no
   integer :: Tropo_Level
   integer :: Sfc_Level
   real, intent(in):: Sfc_Air_Temp
   real, intent(in):: Sfc_Height
   real, intent(out):: Inversion_Top_Height
   real, intent(out):: Inversion_Base_Height
   integer, intent(out):: Top_Lev_Idx
   integer, intent(out):: Base_Lev_Idx
   real, intent(out):: Inversion_Strength
   logical, intent(in):: Dump_Diag
   integer, intent(in):: Lun_Iter_Dump
   integer:: k
   real:: Inversion_Top_Temperature
   real:: Inversion_Base_Temperature
   integer, dimension(size(ACHA_RTM_NWP%P_Prof)):: Inver_Prof_RTM

   Inver_Prof_RTM = Symbol_NO
   Tropo_Level = ACHA_RTM_NWP%Tropo_Level
   Sfc_Level = ACHA_RTM_NWP%Sfc_Level

   do k = Sfc_Level, Tropo_Level, -1

      if (ACHA_RTM_NWP%P_Prof(k) >= MIN_P_INVERSION) then

         if (ACHA_RTM_NWP%T_Prof(k-1) - ACHA_RTM_NWP%T_Prof(k) > DELTA_T_LAYER_INVERSION) then
            Inver_Prof_RTM(k-1:k) = Symbol_YES
         endif

      endif
   enddo

   Top_Lev_Idx =  0
   do k = Tropo_Level,Sfc_Level,1
      if (Inver_Prof_RTM(k) == Symbol_YES .and. Top_Lev_Idx == 0) then
         Top_Lev_Idx = k
         exit
      endif
   enddo

   Base_Lev_Idx = 0
   do k = Sfc_Level, Tropo_Level, -1
      if (Inver_Prof_RTM(k) == Symbol_YES .and. Base_Lev_Idx == 0) then
         Base_Lev_Idx = k
         exit
      endif
   enddo

   Inversion_Strength  = MISSING_VALUE_REAL4
   Inversion_Base_Height = MISSING_VALUE_REAL4
   Inversion_Top_Height = MISSING_VALUE_REAL4

   !---- inversion top height (meters)
   if (Top_Lev_Idx /= 0)  then 
         Inversion_Top_Height = ACHA_RTM_NWP%Z_Prof(Top_Lev_Idx)
         Inversion_Top_Temperature = ACHA_RTM_NWP%T_Prof(Top_Lev_Idx)
   endif

   !---- inversion base height (meters)
   if (Base_Lev_Idx /= 0) then
          Inversion_Base_Height = ACHA_RTM_NWP%Z_Prof(Base_Lev_Idx)
          Inversion_Base_Temperature = ACHA_RTM_NWP%T_Prof(Base_Lev_Idx)
   endif

   !--- assume inversion streches to surface if lowest level is the surface level
   if ((Base_Lev_Idx == Sfc_Level) .and. (Sfc_Height .ner. MISSING_VALUE_REAL4)) then
    Inversion_Base_Height = Sfc_Height
   endif

   if ((Base_Lev_Idx == Sfc_Level) .and.  (Sfc_Air_Temp .ner. MISSING_VALUE_REAL4)) then
          Inversion_Base_Temperature = Sfc_Air_Temp
   endif

   !--- inversion temperature strength
   if (Top_Lev_Idx /= 0 .and. Base_Lev_Idx /= 0)  then 
       Inversion_Strength = Inversion_Top_Temperature - Inversion_Base_Temperature
   endif

   if (Dump_Diag) then
     write(unit=Lun_Iter_Dump,fmt=*) "Inver_Top_Level = ", Top_Lev_Idx
     write(unit=Lun_Iter_Dump,fmt=*) "Inver_Top_Height= ", Inversion_Top_Height
     write(unit=Lun_Iter_Dump,fmt=*) "Inver_Top_Temperature = ", Inversion_Top_Temperature
     write(unit=Lun_Iter_Dump,fmt=*) "Inver_Base_Level = ", Base_Lev_Idx
     write(unit=Lun_Iter_Dump,fmt=*) "Inver_Base_Height= ", Inversion_Base_Height
     write(unit=Lun_Iter_Dump,fmt=*) "Inver_Base_Temperature = ", Inversion_Base_Temperature
     write(unit=Lun_Iter_Dump,fmt=*) "Surface_Temperature = ", Sfc_Air_Temp
     write(unit=Lun_Iter_Dump,fmt=*) "Inver_Strength = ", Inversion_Strength
   endif

 end subroutine DETERMINE_INVERSION_CHARACTERISTICS

 !----------------------------------------------------------------------
 !
 ! Compute the IR Emissivity at a Reference Level
 !
 ! Ref_Level refers to a level index in the profiles
 ! Toa_Radiance = top of atmosphere radiance
 ! Toa_Radiance_Clear = top of atmosphere radiance under clear-skies
 !
 ! 
 ! Black_Body_Rad_Prof_110um_RTM - this is in memory
 !
 !----------------------------------------------------------------------
 function COMPUTE_REFERENCE_LEVEL_EMISSIVITY(Ref_Level,Toa_Radiance, &
                                             Toa_Radiance_Clear, &
                                             Black_Body_Rad_Prof_110um) &
                                             result(Emissivity_Ref_Level)

    integer(kind=int4), intent(in):: Ref_Level
    real (kind=real4), intent(in):: Toa_Radiance  
    real (kind=real4), intent(in):: Toa_Radiance_Clear
    real (kind=real4), intent(in), dimension(:):: Black_Body_Rad_Prof_110um
    real (kind=real4):: Emissivity_Ref_Level

    Emissivity_Ref_Level = &
    (Toa_Radiance - Toa_Radiance_Clear) / &
    (Black_Body_Rad_Prof_110um(Ref_Level) - Toa_Radiance_Clear)

 end function COMPUTE_REFERENCE_LEVEL_EMISSIVITY

 !----------------------------------------------------------------------
 !  Local Linear Radiative Center
 !----------------------------------------------------------------------
 subroutine LOCAL_LINEAR_RADIATIVE_CENTER(Symbol_yes,Symbol_no, &
                                          Meander_Flag, &
                                          Grid_Data, &
                                          Element_Start, Number_Of_Elements, & 
                                          Line_Start, Number_Of_Lines, & 
                                          Max_Grad_Distance, &
                                          Min_Grad_Value, &
                                          Max_Grad_Value, &
                                          Grad_Flag,  &
                                          Missing_LRC_Value, &
                                          Skip_LRC_Mask, &
                                          Min_Grid_Data_Valid, Max_Grid_Data_Valid, &
                                          Elem_Idx_LRC, Line_Idx_LRC)

  integer(kind=int1), intent(in) :: Symbol_yes
  integer(kind=int1), intent(in) :: Symbol_no
  integer, intent(in):: Meander_Flag
  real (kind=real4), intent(in), dimension(:,:) :: Grid_Data
  integer (kind=int4), intent(in):: Element_Start
  integer (kind=int4), intent(in):: Number_of_Elements
  integer (kind=int4), intent(in):: Line_Start
  integer (kind=int4), intent(in):: Number_of_Lines
  integer (kind=int4), intent(in):: Max_Grad_Distance
  real (kind=real4), intent(in):: Min_Grad_Value
  real (kind=real4), intent(in):: Max_Grad_Value
  integer (kind=int4), intent(in):: Grad_Flag
  integer (kind=int4), intent(in):: Missing_LRC_Value
  integer (kind=int1), intent(in), dimension(:,:):: Skip_LRC_Mask
  real (kind=real4), intent(in):: Min_Grid_Data_Valid
  real (kind=real4), intent(in):: Max_Grid_Data_Valid
  integer (kind=int4), intent(out), dimension(:,:):: Elem_Idx_LRC
  integer (kind=int4), intent(out), dimension(:,:):: Line_Idx_LRC
  real, dimension(3,3):: Grad_Array
  integer, dimension(2):: Grad_Indices
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Elem_Idx_Previous
  integer:: Line_Idx_Previous
  integer:: Elem_Idx_Next
  integer:: Line_Idx_Next
  real:: Grad_Temp
  integer:: Element_End
  integer:: Line_End
  integer:: ipoint
  integer:: Elem_Idx_dir
  integer:: Line_Idx_dir
 
  Element_End = Number_of_Elements + Element_Start - 1
  Line_End = Number_of_Lines + Line_Start - 1

  !--- initialize
  Elem_Idx_LRC = Missing_LRC_Value
  Line_Idx_LRC = Missing_LRC_Value

!----------------------------------------------------------------------
! loop through pixels in segment
!----------------------------------------------------------------------
Element_Loop:  do Elem_Idx = Element_Start+1, Element_End-1
Line_Loop:    do Line_Idx = Line_Start+1, Line_End-1

      !--- skip data due to mask
      if (Skip_LRC_Mask(Elem_Idx,Line_Idx) == Symbol_YES) cycle

      !-- check for out of bounds data
      if (Grad_Flag ==  1 .and. Grid_Data(Elem_Idx,Line_Idx) < Min_Grid_Data_Valid) cycle
      if (Grad_Flag ==  -1 .and. Grid_Data(Elem_Idx,Line_Idx) > Max_Grid_Data_Valid) cycle

      !-- check for data that already meets LRC criteria
      if ((Grad_Flag ==  1 .and. Grid_Data(Elem_Idx,Line_Idx) > Max_Grid_Data_Valid) .or. &
          (Grad_Flag ==  -1 .and. Grid_Data(Elem_Idx,Line_Idx) < Min_Grid_Data_Valid)) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx
              cycle
      endif

      !--- initialize previous variables
      Elem_Idx_Previous = Elem_Idx
      Line_Idx_Previous = Line_Idx

      !---- go long gradient and check for a reversal or saturation
Gradient_Loop:    do ipoint = 1,Max_Grad_Distance

        !--- compute local gradient, find strongest gradient in 3x3 array and compute direction
        if (ipoint == 1 .or. Meander_Flag == Symbol_YES) then

         !--- construct 3x3 array for analysis
         Grad_Array =  &
           Grid_Data(Elem_Idx_Previous-1:Elem_Idx_Previous+1,Line_Idx_Previous-1:Line_Idx_Previous+1) -  &
           Grid_Data(Elem_Idx_Previous,Line_Idx_Previous)

         !--- look for bad data
         if (minval(Grad_Array) == MISSING_VALUE_REAL4) exit 

         !--- compute local gradients, find strongest gradient
         if (Grad_Flag == 1) then
          Grad_Indices = maxloc(Grad_Array)
         else
          Grad_Indices = minloc(Grad_Array)
         endif 

         !--- compute direction
         Elem_Idx_Dir = Grad_Indices(1) - 2
         Line_Idx_Dir = Grad_Indices(2) - 2

         !--- check for pixels that are located at  minima/maxima
         if (Elem_Idx_Dir == 0 .and. Line_Idx_Dir == 0) then
           Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Previous
           Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Previous
           exit
         endif

         !--- on first step, only proceed if gradient magnitude exceeds a threshold
         if (ipoint == 1) then
            if (abs(Grad_Array(Grad_Indices(1),Grad_Indices(2))) < Min_Grad_Value) then
              exit
            endif
         endif

         !--- check for going up to steep of a gradient
         if (abs(Grad_Array(Grad_Indices(1),Grad_Indices(2))) > Max_Grad_Value) then
           exit
         endif

        endif

        !-- select next point on the path
        Elem_Idx_Next = Elem_Idx_Previous + Elem_Idx_Dir
        Line_Idx_Next = Line_Idx_Previous + Line_Idx_Dir

        !--- check for hitting segment boundaries
        if (Elem_Idx_Next == Element_Start .or. Elem_Idx_Next == Element_End .or. &
             Line_Idx_Next == Line_Start .or. Line_Idx_Next == Line_End) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Previous
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Previous
              exit
         endif

         !--- check for hitting bad data
         if (Skip_LRC_Mask(Elem_Idx_Next,Line_Idx_Next) == Symbol_YES) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Previous
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Previous
              exit
         endif

         !--- check for sign reversal
         if (Meander_Flag == Symbol_NO) then

          Grad_Temp = Grid_Data(Elem_Idx_Next,Line_Idx_Next) -  &
                      Grid_Data(Elem_Idx_Previous,Line_Idx_Previous)

          if (Grad_Flag * Grad_Temp < 0) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Previous
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Previous
              exit
          endif
         endif

         !--- check for saturation
         if (Grad_Flag == 1 .and. Grid_Data(Elem_Idx_Next,Line_Idx_Next) > Max_Grid_Data_Valid) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Next
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Next
              exit
         endif
         if (Grad_Flag == -1 .and. Grid_Data(Elem_Idx_Next,Line_Idx_Next) < Min_Grid_Data_Valid) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Next
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Next
              exit
         endif

         !--- store position
         Elem_Idx_Previous = Elem_Idx_Next
         Line_Idx_Previous = Line_Idx_Next

      enddo Gradient_Loop

    end do Line_Loop
  end do Element_Loop

end subroutine LOCAL_LINEAR_RADIATIVE_CENTER
!----------------------------------------------------------------------
! Local Routine for a Standard Deviation
!
! Data_Array - input array of real numbers
! Invalid_Mask = 0 for good pixels, 1 for invalid pixels
! Stddev = Standard Deviation for valid pixels in Data Array
!
! Num_Good = number of valid data point in array
!
! If Num_Good < 2, we do nothing
!----------------------------------------------------------------------
function COMPUTE_STANDARD_DEVIATION(Data_Array,Invalid_Mask) Result(Stddev_of_Array_r4)
   real(kind=real4), dimension(:,:), intent(in):: Data_Array
   integer(kind=int1), dimension(:,:), intent(in):: Invalid_Mask
   real:: Stddev_of_Array_r4
   real(kind=real8):: Stddev_of_Array_r8
   real(kind=real8):: Data_Sum
   real(kind=real8):: Data_Sum_Squared
   real(kind=real8):: Num_Good
   real(kind=real8):: temp

   Num_Good = real(sum(1 - Invalid_Mask))

   if (Num_Good == 0.0) then
    Stddev_of_Array_r8 = MISSING_VALUE_REAL4
   elseif (Num_Good == 1.0) then
    Stddev_of_Array_r8 = 0.0
   else
    Data_Sum = sum(Data_Array * (1.0 - Invalid_Mask))
    Data_Sum_Squared = sum((Data_Array*(1.0-Invalid_Mask))**2)
    temp = Data_Sum_Squared / Num_Good - (Data_Sum/Num_Good)**2
    if (temp > 0.0) then
      Stddev_of_Array_r8 = sqrt(temp)
    else
      Stddev_of_Array_r8 = 0.0
    endif
   endif 

    Stddev_of_Array_r4 = real(Stddev_of_Array_r8,kind=real4)

end function

!---------------------------------------------------------------------------
! Compute Parallax Correction
!
! This routine generates new Lat and Lon arrays that are parallax
! corrected based on the cloud height
!
! Input: Senzen - sensor viewing zenith angle (deg) 
!        Senaz  - sensor azimuth angle (deg)
!        Lat - uncorrected Latitude (deg)
!        Lon  - uncorrected longitude (deg)
!        Zsfc  - surface elevation (m)
!        Zcld  - cloud height (m)
!
! Output
!       Lat_Pc - corrected Latitude
!       Lon_Pc - corrected longitude
!
!---------------------------------------------------------------------------
subroutine PARALLAX_ACHA(Zcld,Zsfc,Lat,Lon,Senzen,Senaz,Lat_Pc,Lon_Pc) 
   real, intent(in), dimension(:,:):: Zcld
   real, intent(in), dimension(:,:):: Zsfc
   real, intent(in), dimension(:,:):: Lat
   real, intent(in), dimension(:,:):: Lon
   real, intent(in), dimension(:,:):: Senzen
   real, intent(in), dimension(:,:):: Senaz
   real, intent(out), dimension(:,:):: Lat_Pc
   real, intent(out), dimension(:,:):: Lon_Pc
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: Num_Elem
   integer:: Num_Line
   real:: Total_Displacement
   real:: Delta_Lon
   real:: Delta_Lat
   real:: Lon_Spacing_Per_m
   real,parameter:: Lat_Spacing_Per_m = 8.9932e-06   ! ( = 1.0/111000.0 m )

   Num_Elem = size(Zcld,1) 
   Num_Line = size(Zcld,2)

   !--- initialize output to standard values
   Lat_Pc = Lat
   Lon_Pc = Lon

   !--- loop over pixels in segment
   element_loop: do Elem_Idx = 1, Num_Elem
    line_loop: do Line_Idx = 1, Num_Line

     !--- check for valid data
     if (Zcld(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4 .or. &
         Senzen(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) cycle

     !--- compute correction
     Total_Displacement = max(0.0,tan(Senzen(Elem_Idx,Line_Idx)*Dtor)* &
                                     (Zcld(Elem_Idx,Line_Idx) - Zsfc(Elem_Idx,Line_Idx)))

     Lon_Spacing_Per_m = Lat_Spacing_Per_m / cos(Lat(Elem_Idx,Line_Idx)*Dtor)

     Delta_Lon = sin(Senaz(Elem_Idx,Line_Idx)*Dtor)*Total_Displacement * Lon_Spacing_Per_m
     Delta_Lat = cos(Senaz(Elem_Idx,Line_Idx)*Dtor)*Total_Displacement * Lat_Spacing_Per_m

     !--- generate output positions
     Lat_Pc(Elem_Idx,Line_Idx) = Lat(Elem_Idx,Line_Idx) + Delta_Lat
     Lon_Pc(Elem_Idx,Line_Idx) = Lon(Elem_Idx,Line_Idx) + Delta_Lon

    enddo line_loop
   enddo element_loop

end subroutine PARALLAX_ACHA 
!----------------------------------------------------------------------
!--- check that the Acha_Mode is consistent with available channels
!--- if consistent, Acha_Mode_Error_Flag = 0, if not, flag = 1
!----------------------------------------------------------------------
subroutine CHECK_ACHA_MODE( &
                           Acha_Mode_Input, &
                           Chan_On_038um, &
                           Chan_On_062um, &
                           Chan_On_067um, &
                           Chan_On_085um, &
                           Chan_On_110um, &
                           Chan_On_120um, &
                           Chan_On_133um, &
                           Acha_Mode_Error_Flag)

   character (len=*), intent(in) :: Acha_Mode_Input
   integer, intent(in) :: Chan_On_038um
   integer, intent(in) :: Chan_On_062um
   integer, intent(in) :: Chan_On_067um
   integer, intent(in) :: Chan_On_085um
   integer, intent(in) :: Chan_On_110um
   integer, intent(in) :: Chan_On_120um
   integer, intent(in) :: Chan_On_133um
   integer, intent(out) :: Acha_Mode_Error_Flag

   Acha_Mode_Error_Flag = 0

   if (Chan_On_110um == Symbol%NO) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_038um == Symbol%NO) .and. &
       (index(Acha_Mode_Input,'038') > 0)) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_062um == Symbol%NO) .and. &
       (index(Acha_Mode_Input,'062') > 0)) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_067um == Symbol%NO) .and. &
       (index(Acha_Mode_Input,'067') > 0)) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_085um == Symbol%NO) .and. &
       (index(Acha_Mode_Input,'085') > 0)) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_120um == Symbol%NO) .and. &
       (index(Acha_Mode_Input,'120') > 0)) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_133um == Symbol%NO) .and. &
       (index(Acha_Mode_Input,'133') > 0)) then
       Acha_Mode_Error_Flag = 1
       return
   endif

end subroutine CHECK_ACHA_MODE

!------------------------------------------------------------------------------
! Null Pixel Level Pointers 
!------------------------------------------------------------------------------
subroutine NULL_PIX_POINTERS(Input, ACHA_RTM_NWP)

   type(acha_input_struct), intent(inout) :: Input
   type(acha_rtm_nwp_struct), intent(inout) :: ACHA_RTM_NWP

   ACHA_RTM_NWP%T_Prof => NULL()

   ACHA_RTM_NWP%Z_Prof => NULL() 

   if (Input%Chan_On_038um == Symbol%YES) then
     ACHA_RTM_NWP%Atm_Rad_Prof_038um =>  NULL()
     ACHA_RTM_NWP%Atm_Trans_Prof_038um =>  NULL()
     ACHA_RTM_NWP%Black_Body_Rad_Prof_038um => NULL()
   endif
   if (Input%Chan_On_062um == Symbol%YES) then
     ACHA_RTM_NWP%Atm_Rad_Prof_062um =>  NULL()
     ACHA_RTM_NWP%Atm_Trans_Prof_062um =>  NULL()
     ACHA_RTM_NWP%Black_Body_Rad_Prof_062um => NULL()
   endif
   if (Input%Chan_On_067um == Symbol%YES) then
     ACHA_RTM_NWP%Atm_Rad_Prof_067um =>  NULL()
     ACHA_RTM_NWP%Atm_Trans_Prof_067um =>  NULL()
     ACHA_RTM_NWP%Black_Body_Rad_Prof_067um => NULL()
   endif
   if (Input%Chan_On_073um == Symbol%YES) then
     ACHA_RTM_NWP%Atm_Rad_Prof_073um =>  NULL()
     ACHA_RTM_NWP%Atm_Trans_Prof_073um =>  NULL()
     ACHA_RTM_NWP%Black_Body_Rad_Prof_073um => NULL()
   endif
   if (Input%Chan_On_085um == Symbol%YES) then
     ACHA_RTM_NWP%Atm_Rad_Prof_085um =>  NULL()
     ACHA_RTM_NWP%Atm_Trans_Prof_085um =>  NULL()
     ACHA_RTM_NWP%Black_Body_Rad_Prof_085um => NULL()
   endif
   if (Input%Chan_On_097um == Symbol%YES) then
     ACHA_RTM_NWP%Atm_Rad_Prof_097um =>  NULL()
     ACHA_RTM_NWP%Atm_Trans_Prof_097um =>  NULL()
     ACHA_RTM_NWP%Black_Body_Rad_Prof_097um => NULL()
   endif
   if (Input%Chan_On_104um == Symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_104um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_104um => NULL()
      ACHA_RTM_NWP%Black_Body_Rad_Prof_104um => NULL()
   endif
   if (Input%Chan_On_110um == Symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_110um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_110um => NULL()
      ACHA_RTM_NWP%Black_Body_Rad_Prof_110um => NULL()
   endif
   if (Input%Chan_On_120um == Symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_120um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_120um => NULL()
      ACHA_RTM_NWP%Black_Body_Rad_Prof_120um => NULL()
   endif
   if (Input%Chan_On_133um == Symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_133um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_133um => NULL()
      ACHA_RTM_NWP%Black_Body_Rad_Prof_133um => NULL()
   endif
   if (Input%Chan_On_136um == Symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_136um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_136um => NULL()
      ACHA_RTM_NWP%Black_Body_Rad_Prof_136um => NULL()
   endif
   if (Input%Chan_On_139um == Symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_139um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_139um => NULL()
      ACHA_RTM_NWP%Black_Body_Rad_Prof_139um => NULL()
   endif
   if (Input%Chan_On_142um == Symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_142um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_142um => NULL()
      ACHA_RTM_NWP%Black_Body_Rad_Prof_142um => NULL()
   endif
 
end subroutine NULL_PIX_POINTERS
!====================================================================
!  record svn version as a global variable for output to hdf
!====================================================================
subroutine SET_ACHA_VERSION(Acha_Version)
   character(len=*):: Acha_Version
   Acha_Version = "$Id: acha_module.f90 4092 2021-03-02 22:05:14Z heidinger $"
end subroutine SET_ACHA_VERSION

!----------------------------------------------------------------------
! Empirical Lapse Rate
!----------------------------------------------------------------------
function EMPIRICAL_LAPSE_RATE(Tsfc, Tc, land_flag) result(lapse_rate)
  real, intent(in):: Tsfc
  real, intent(in):: Tc
  integer, intent(in):: land_flag  !(0=ocean,1=land)
  real:: Tcs
  real:: lapse_rate
  integer:: its, itcs

  Tcs = Tc - Tsfc
  its = int((Tsfc - ts_min) / dts) + 1
  its = max(1,min(nts,its))
  itcs = int((Tcs - tcs_min) / dtcs) + 1
  itcs = max(1,min(ntcs,itcs))

  if (land_flag == 0) then 
    lapse_rate = ocean_lapse_rate_table(its,itcs)
  else
    lapse_rate = land_lapse_rate_table(its,itcs)
  endif

end function EMPIRICAL_LAPSE_RATE

!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
  function GET_LUN_ACHA() result( Lun )


    ! -----------------
    ! Type declarations
    ! -----------------

    integer :: Lun
    logical :: File_Open

    ! --------------------------------------------
    ! Initialise logical unit number and file_open
    ! --------------------------------------------

    Lun = 9
    File_Open = .TRUE.


    ! ------------------------------
    ! Start open loop for lun search
    ! ------------------------------

    lun_search: do

      ! -- Increment logical unit number
      Lun = Lun + 1

      ! -- Check if file is open
      inquire( Lun, OPENED = File_Open )

      ! -- Is this lun available?
      if ( .not. File_Open ) EXIT Lun_Search

    enddo lun_search

  end function GET_LUN_ACHA
!-------------------------------------------------------------------------------------------
! Compute the y and y_variance vectors which depend on the chosen Mode
!-------------------------------------------------------------------------------------------
subroutine COMPUTE_Y(Acha_Mode_Flag,Input,Element_Idx_Min, Line_Idx_Min, Elem_Idx,Line_Idx, &
                     y, y_variance, Dump_Diag, Lun_Iter_Dump)   

   character(len=*), intent(in):: Acha_Mode_Flag
   type(acha_input_struct), intent(in) :: Input
   real, intent(out), dimension(:):: y
   real, intent(out), dimension(:):: y_variance
   logical, intent(in):: Dump_Diag
   integer, intent(in):: Lun_Iter_Dump
   integer:: i1, i2, j1, j2
   integer:: Line_Idx_Min, Line_Idx, Element_Idx_Min, Elem_Idx
   real (kind=real4):: Bt_110um_Std
   real (kind=real4):: Btd_110um_038um_Std
   real (kind=real4):: Btd_110um_062um_Std
   real (kind=real4):: Btd_110um_067um_Std
   real (kind=real4):: Btd_110um_073um_Std
   real (kind=real4):: Btd_110um_085um_Std
   real (kind=real4):: Btd_110um_097um_Std
   real (kind=real4):: Btd_110um_104um_Std
   real (kind=real4):: Btd_110um_120um_Std
   real (kind=real4):: Btd_110um_133um_Std
   real (kind=real4):: Btd_110um_136um_Std
   real (kind=real4):: Btd_110um_139um_Std
   real (kind=real4):: Btd_110um_142um_Std
   real (kind=real4):: Bt_110um_Noise
   real (kind=real4):: Btd_110um_038um_Noise
   real (kind=real4):: Btd_110um_062um_Noise
   real (kind=real4):: Btd_110um_067um_Noise
   real (kind=real4):: Btd_110um_073um_Noise
   real (kind=real4):: Btd_110um_085um_Noise
   real (kind=real4):: Btd_110um_097um_Noise
   real (kind=real4):: Btd_110um_104um_Noise
   real (kind=real4):: Btd_110um_120um_Noise
   real (kind=real4):: Btd_110um_133um_Noise
   real (kind=real4):: Btd_110um_136um_Noise
   real (kind=real4):: Btd_110um_139um_Noise
   real (kind=real4):: Btd_110um_142um_Noise

   !-----------------------------------------------------------------------
   ! compute needed channel 3x3 standard deviations
   !-----------------------------------------------------------------------
   j1 = max(Line_Idx_Min, Line_Idx - 1)
   j2 = min(Input%Number_of_Lines, Line_Idx + 1)
   i1 = max(Element_Idx_Min, Elem_Idx - 1) 
   i2 = min(Input%Number_of_Elements, Elem_Idx + 1)

   !--- At this point, for GOES-17 bad data, Bt_110 um should be Bt_104 um.
   Bt_110um_Std = COMPUTE_STANDARD_DEVIATION( Input%Bt_110um(i1:i2,j1:j2),Input%Invalid_Data_Mask(i1:i2,j1:j2))

   Bt_110um_Noise = COMPUTE_BT_NOISE(Input%Chan_Idx_110um,Input%Bt_110um(Elem_Idx,Line_Idx))

   if (index(Acha_Mode_Flag,'038') > 0) then
    Btd_110um_038um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_038um(i1:i2,j1:j2),&
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_038um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_038um,Input%Bt_038um(Elem_Idx,Line_Idx))**2)
   endif

   if (index(Acha_Mode_Flag,'062') > 0) then
    Btd_110um_062um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_062um(i1:i2,j1:j2),&
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_062um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_062um,Input%Bt_062um(Elem_Idx,Line_Idx))**2)
   endif

   if (index(Acha_Mode_Flag,'067') > 0) then
    Btd_110um_067um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_067um(i1:i2,j1:j2),&
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_067um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_067um,Input%Bt_067um(Elem_Idx,Line_Idx))**2)
   endif
   if (index(Acha_Mode_Flag,'073') > 0) then
    Btd_110um_073um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_073um(i1:i2,j1:j2),&
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_073um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_073um,Input%Bt_073um(Elem_Idx,Line_Idx))**2)
   endif

   if (index(Acha_Mode_Flag,'085') > 0) then
    Btd_110um_085um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_085um(i1:i2,j1:j2), &
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_085um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_085um,Input%Bt_085um(Elem_Idx,Line_Idx))**2)
   endif

   if (index(Acha_Mode_Flag,'097') > 0) then
    Btd_110um_097um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_097um(i1:i2,j1:j2), &
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_097um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_097um,Input%Bt_097um(Elem_Idx,Line_Idx))**2)
   endif

   !--- If the use 104um flag has been set, this STD will be 0, as the 11um has
   !--- been substituted with the 104um.
   if (index(Acha_Mode_Flag,'104') > 0) then
    Btd_110um_104um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_104um(i1:i2,j1:j2), &
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_104um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_104um,Input%Bt_104um(Elem_Idx,Line_Idx))**2)
   endif

   if (index(Acha_Mode_Flag,'120') > 0) then
    Btd_110um_120um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_120um(i1:i2,j1:j2), &
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_120um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_120um,Input%Bt_120um(Elem_Idx,Line_Idx))**2)
   endif

   if (index(Acha_Mode_Flag,'133') > 0) then
    Btd_110um_133um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_133um(i1:i2,j1:j2), &
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_133um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_133um,Input%Bt_133um(Elem_Idx,Line_Idx))**2)
   endif

   if (index(Acha_Mode_Flag,'136') > 0) then
    Btd_110um_136um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_136um(i1:i2,j1:j2), &
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_136um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_136um,Input%Bt_136um(Elem_Idx,Line_Idx))**2)
   endif

   if (index(Acha_Mode_Flag,'139') > 0) then
    Btd_110um_139um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_139um(i1:i2,j1:j2), &
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_139um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_139um,Input%Bt_139um(Elem_Idx,Line_Idx))**2)
   endif

   if (index(Acha_Mode_Flag,'142') > 0) then
    Btd_110um_142um_Std = COMPUTE_STANDARD_DEVIATION(Input%Bt_110um(i1:i2,j1:j2) -  Input%Bt_142um(i1:i2,j1:j2), &
                                                     Input%Invalid_Data_Mask(i1:i2,j1:j2))
    Btd_110um_142um_Noise = sqrt(Bt_110um_Noise**2 + COMPUTE_BT_NOISE(Input%Chan_Idx_142um,Input%Bt_142um(Elem_Idx,Line_Idx))**2)
   endif
   
   !--- y - the observation output vector
   select case(trim(Acha_Mode_Flag))
     case('110')
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
     case("038_110")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_038um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_038um_Std**2 + Btd_110um_038um_Noise**2
     case("067_110")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_067um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_067um_Std**2 + Btd_110um_067um_Noise**2
     case("110_120")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
     case("110_133")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_133um_Std**2 + Btd_110um_133um_Noise**2
     case("085_110_120")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_085um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_085um_Std**2 + Btd_110um_085um_Noise**2
       y_variance(3) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
     case("067_110_120")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_067um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_067um_Std**2 + Btd_110um_067um_Noise**2
       y_variance(3) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
     case("067_110_133")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_067um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_067um_Std**2 + Btd_110um_067um_Noise**2
       y_variance(3) = Btd_110um_133um_Std**2 + Btd_110um_133um_Noise**2
     case("110_120_133")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
       y_variance(3) = Btd_110um_133um_Std**2 + Btd_110um_133um_Noise**2
     case("067_085_110")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_067um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_085um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_067um_Std**2 + Btd_110um_067um_Noise**2
       y_variance(3) = Btd_110um_085um_Std**2 + Btd_110um_085um_Noise**2
      case("085_110_120_133")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_085um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y(4) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_085um_Std**2 + Btd_110um_085um_Noise**2
       y_variance(3) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
       y_variance(4) = Btd_110um_133um_Std**2 + Btd_110um_133um_Noise**2
      case("067_085_110_120")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_067um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_085um(Elem_Idx,Line_Idx)
       y(4) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_067um_Std**2 + Btd_110um_067um_Noise**2
       y_variance(3) = Btd_110um_085um_Std**2 + Btd_110um_085um_Noise**2
       y_variance(4) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
      case("062_085_110_120_133")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_062um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_085um(Elem_Idx,Line_Idx)
       y(4) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y(5) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_062um_Std**2 + Btd_110um_062um_Noise**2
       y_variance(3) = Btd_110um_085um_Std**2 + Btd_110um_085um_Noise**2
       y_variance(4) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
       y_variance(5) = Btd_110um_133um_Std**2 + Btd_110um_133um_Noise**2
      case("067_085_110_120_133")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_067um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_085um(Elem_Idx,Line_Idx)
       y(4) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y(5) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_067um_Std**2 + Btd_110um_067um_Noise**2
       y_variance(3) = Btd_110um_085um_Std**2 + Btd_110um_085um_Noise**2
       y_variance(4) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
       y_variance(5) = Btd_110um_133um_Std**2 + Btd_110um_133um_Noise**2
     case("110_133_136_139_142")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_136um(Elem_Idx,Line_Idx)
       y(4) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_139um(Elem_Idx,Line_Idx)
       y(5) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_142um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_133um_Std**2 + Btd_110um_133um_Noise**2
       y_variance(3) = Btd_110um_136um_Std**2 + Btd_110um_136um_Noise**2
       y_variance(4) = Btd_110um_139um_Std**2 + Btd_110um_139um_Noise**2
       y_variance(5) = Btd_110um_142um_Std**2 + Btd_110um_142um_Noise**2
     case("085_110_120_133_136_139_142")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_085um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y(4) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y(5) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_136um(Elem_Idx,Line_Idx)
       y(6) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_139um(Elem_Idx,Line_Idx)
       y(7) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_142um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_085um_Std**2 + Btd_110um_085um_Noise**2
       y_variance(3) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
       y_variance(4) = Btd_110um_133um_Std**2 + Btd_110um_133um_Noise**2
       y_variance(5) = Btd_110um_136um_Std**2 + Btd_110um_136um_Noise**2
       y_variance(6) = Btd_110um_139um_Std**2 + Btd_110um_139um_Noise**2
       y_variance(7) = Btd_110um_142um_Std**2 + Btd_110um_142um_Noise**2
     case("062_067_073_085_104_110_120_133")
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_062um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_067um(Elem_Idx,Line_Idx)
       y(4) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_073um(Elem_Idx,Line_Idx)
       y(5) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_085um(Elem_Idx,Line_Idx)
       y(6) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_104um(Elem_Idx,Line_Idx)
       y(7) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y(8) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_062um_Std**2 + Btd_110um_062um_Noise**2
       y_variance(3) = Btd_110um_067um_Std**2 + Btd_110um_067um_Noise**2
       y_variance(4) = Btd_110um_073um_Std**2 + Btd_110um_073um_Noise**2
       y_variance(5) = Btd_110um_085um_Std**2 + Btd_110um_085um_Noise**2
       y_variance(6) = Btd_110um_104um_Std**2 + Btd_110um_104um_Noise**2 
       y_variance(7) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
       y_variance(8) = Btd_110um_133um_Std**2 + Btd_110um_133um_Noise**2
     case DEFAULT
       y(1) =  Input%Bt_110um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_110um(Elem_Idx,Line_Idx) -  Input%Bt_120um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_110um_Std**2 + Bt_110um_Noise**2
       y_variance(2) = Btd_110um_120um_Std**2 + Btd_110um_120um_Noise**2
   end select

   if (Dump_Diag) then 
     write(unit=Lun_Iter_Dump,fmt=*) "y = ", y
     write(unit=Lun_Iter_Dump,fmt=*) "y variance = ", y_variance
   endif

end subroutine COMPUTE_Y
!------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------
subroutine COMPUTE_META_DATA(Cloud_Phase, USE_LRC_FLAG, Cloud_Type,Meta_Data_Flags)
   integer(kind=int1), dimension(:), intent(out):: Meta_Data_Flags
   integer, intent(in):: Cloud_Phase
   logical, intent(in):: USE_LRC_FLAG
   integer(kind=int1), intent(in):: Cloud_Type

   Meta_Data_Flags(1) = Symbol%YES
   if (Cloud_Phase == Symbol%ICE_PHASE) then
       Meta_Data_Flags(3) = Symbol%YES
   else
       Meta_Data_Flags(3) = Symbol%NO
   endif
   if (USE_LRC_FLAG) then
       Meta_Data_Flags(4) = Symbol%YES
   else
       Meta_Data_Flags(4) = Symbol%NO
   endif
   if (Cloud_Type == Symbol%OVERLAP_TYPE) then
       Meta_Data_Flags(5) = Symbol%YES
   else
       Meta_Data_Flags(5) = Symbol%NO
   endif
   Meta_Data_Flags(6) = Symbol%NO     !lower cloud interpoLation
   Meta_Data_Flags(7) = Symbol%NO     !low level inversion
   Meta_Data_Flags(8) = Symbol%NO     !NWP profile inversion

 end subroutine COMPUTE_META_DATA
!--------------------------------------------------------------------
! This routine dissects the ACHA_Mode_Flag string into its parts to
! 1. determine the number of channel in the retrieval
! 2. determine which channels are in using the CLAVR-x channel ids
!--------------------------------------------------------------------
subroutine  DETERMINE_NUMBER_OF_CHANNELS(Acha_Mode_Flag, Num_Obs)
  character (len=*), intent(in):: Acha_Mode_Flag
  integer, intent(out):: Num_Obs
  character(len=3):: Chan_String
  integer:: Obs_Idx, istart, iend, iobs
 

  Num_Obs = COUNTSUBSTRING(trim(ACHA_Mode_Flag),'_') + 1

  allocate(Chan_Idx_y(Num_Obs)) 
   
  !--- the first channel is always the 11 micron channel
  Obs_Idx = 1
  Chan_Idx_y(Obs_Idx) = 31   !always true

  !--- extract string names and determine channels 
  !--- ACHA_Mode has channels list in ascending wavelength
  istart = 1

  do Iobs = 1, Num_Obs
     iend = istart + 2
     Chan_String = ACHA_Mode_Flag(istart:iend)
     istart = iend + 2
     if (Chan_String == "110") cycle
     Obs_Idx = Obs_Idx + 1
     if (Chan_String == "038") Chan_Idx_y(Obs_Idx) = 20
     if (Chan_String == "062") Chan_Idx_y(Obs_Idx) = 37
     if (Chan_String == "067") Chan_Idx_y(Obs_Idx) = 27
     if (Chan_String == "073") Chan_Idx_y(Obs_Idx) = 28
     if (Chan_String == "085") Chan_Idx_y(Obs_Idx) = 29
     if (Chan_String == "097") Chan_Idx_y(Obs_Idx) = 30
     if (Chan_String == "104") Chan_Idx_y(Obs_Idx) = 38
     if (Chan_String == "120") Chan_Idx_y(Obs_Idx) = 32
     if (Chan_String == "133") Chan_Idx_y(Obs_Idx) = 33
     if (Chan_String == "136") Chan_Idx_y(Obs_Idx) = 34
     if (Chan_String == "139") Chan_Idx_y(Obs_Idx) = 35
     if (Chan_String == "142") Chan_Idx_y(Obs_Idx) = 36
  enddo

end subroutine  DETERMINE_NUMBER_OF_CHANNELS
!-------------------------------------------------------
!---  for low clouds over water, force fixed lapse rate estimate of height
!-------------------------------------------------------
subroutine COMPUTE_HEIGHT_FROM_LAPSE_RATE(ACHA_RTM_NWP, &
                                          Snow_Class, &
                                          Surface_Type, &
                                          Cloud_Type, &
                                          Surface_Temperature, & 
                                          Surface_Elevation, &
                                          Max_Delta_T_Inversion, &
                                          Tc, &
                                          Zc, &
                                          Pc, &
                                          Inversion_Flag)
 type (acha_rtm_nwp_struct), intent(in):: ACHA_RTM_NWP
 integer(kind=int1), intent(in):: Snow_Class, Surface_Type, Cloud_Type
 real, intent(in):: Surface_Temperature, Surface_Elevation, &
                    Tc, Max_Delta_T_Inversion
 real, intent(out):: Zc, Pc
 integer(kind=int1), intent(out):: Inversion_Flag

 real:: Delta_Cld_Temp_Sfc_Temp, Lapse_Rate, R4_Dummy
 integer:: Lev_Idx

 Delta_Cld_Temp_Sfc_Temp = Surface_Temperature - Tc
 Lapse_Rate = MISSING_VALUE_REAL4

 Inversion_Flag = 0_int1

 if (Tc .eqr. MISSING_VALUE_REAL4) return

 !--- New prefered method is to take out the Snow_Class check.
 if ((Cloud_Type == Symbol%WATER_TYPE) .or. &
     (Cloud_Type == Symbol%FOG_TYPE) .or. &
     (Cloud_Type == Symbol%SUPERCOOLED_TYPE)) then

     if (Delta_Cld_Temp_Sfc_Temp <  MAX_DELTA_T_INVERSION) then

       !-- select lapse rate  (k/km)
       if (Surface_Type == Symbol%WATER_SFC) then 
           Lapse_Rate = EMPIRICAL_LAPSE_RATE(Surface_Temperature,Tc, 0)
       else
           Lapse_Rate = EMPIRICAL_LAPSE_RATE(Surface_Temperature,Tc, 1)
       endif

       !--- constrain lapse rate to be with -2 and -10 K/km
       Lapse_Rate = min(-2.0,max(-10.0,Lapse_Rate))

       !--- convert lapse rate to K/m
       Lapse_Rate = Lapse_Rate / 1000.0  !(K/m)

       !-- compute height
       Zc = -1.0*Delta_Cld_Temp_Sfc_Temp/Lapse_Rate + Surface_Elevation

       !--- Some negative cloud heights are observed because of bad height
       !--- NWP profiles.
       if (Zc < 0) then
          Zc = ZC_FLOOR
       endif

       !--- compute pressure
       call KNOWING_Z_COMPUTE_T_P(ACHA_RTM_NWP,Pc,R4_Dummy,Zc,Lev_Idx)

       Inversion_Flag = 1_int1

     endif
  endif
end subroutine COMPUTE_HEIGHT_FROM_LAPSE_RATE

!------------------------------------------------------------------------------------------------------
! checkout output for exceeding expected limits and clip if necessary
! also, check that Zc and Pc are not below the surface
!------------------------------------------------------------------------------------------------------
subroutine QUALITY_CONTROL_OUTPUT(Tc, Pc, Zc, Ec, Beta, Surface_Elevation, Surface_Pressure,Clip_Flag)
  real,intent(inout):: Tc
  real,intent(inout):: Pc
  real,intent(inout):: Zc
  real,intent(inout):: Ec
  real,intent(inout):: Beta
  real,intent(in):: Surface_Elevation
  real,intent(in):: Surface_Pressure
  logical, intent(out):: Clip_Flag
  real :: Zc_Floor_Temp
  real :: Pc_Ceiling_Temp
  real, parameter:: Zc_Roundoff_Offset = 0.0  !m
  real, parameter:: Pc_Roundoff_Offset = 0.0  !hPa
 
  Clip_Flag = .false.

  if (Zc /= MISSING_VALUE_REAL4) then
     Zc_Floor_Temp = ZC_FLOOR
      if (Surface_Elevation >= 0.0) then
!     if (Surface_Elevation /= MISSING_VALUE_REAL4) then
        Zc_Floor_Temp = Surface_Elevation + Zc_Roundoff_Offset
     endif
     if (Zc < Zc_Floor_Temp .or. Zc > ZC_CEILING) then
        Clip_Flag = .true.
        Zc = min(ZC_CEILING,max(Zc_Floor_Temp, Zc))
     endif
  endif

  if (Pc /= MISSING_VALUE_REAL4) then
     Pc_Ceiling_Temp = PC_CEILING
     if (Surface_Pressure /= MISSING_VALUE_REAL4) then
        Pc_Ceiling_Temp = Surface_Pressure - Pc_Roundoff_Offset
     endif
     if (Pc < PC_FLOOR .or. Pc > Pc_Ceiling_Temp) then
        Clip_Flag = .true.
        Pc = min(Pc_Ceiling_Temp,max(PC_FLOOR, Pc))
     endif
  endif

  if (Tc /= MISSING_VALUE_REAL4) then
     if (Tc < Tc_Floor .or. Tc > TC_CEILING) then
        Clip_Flag = .true.
        Tc = min(TC_CEILING,max(TC_FLOOR, Tc))
     endif
  endif

  if (Ec /= MISSING_VALUE_REAL4) then
     if (Ec < Ec_Floor .or. Ec > EC_CEILING) then
        Clip_Flag = .true.
        Ec = min(EC_CEILING,max(EC_FLOOR, Ec))
     endif
  endif

  if (Beta /= MISSING_VALUE_REAL4) then
     if (Beta < Beta_Floor .or. Beta > BETA_CEILING) then
        Clip_Flag = .true.
        Beta = min(BETA_CEILING,max(BETA_FLOOR, Beta))
     endif
  endif

end subroutine QUALITY_CONTROL_OUTPUT
!----------------------------------------------------------------------------
! derive a cloud type based on ACHA results
!----------------------------------------------------------------------------
subroutine DETERMINE_ACHA_CLOUD_TYPE(Input,Fail_Flag,Symbol,Output)
   type(acha_input_struct), intent(inout) :: Input
   integer, dimension(:,:), intent(in):: Fail_Flag
   type(acha_symbol_struct), intent(in) :: Symbol
   type(acha_output_struct), intent(inout) :: Output

   Output%Cloud_Type = Symbol%UNKNOWN_TYPE
   where (Input%Invalid_Data_Mask /= Symbol%YES)
      Output%Cloud_Type = Symbol%WATER_TYPE
   endwhere
   where (Fail_Flag == 1)
      Output%Cloud_Type = Symbol%UNKNOWN_TYPE
   endwhere
   where (Output%Ice_Probability >= 0.50)
      Output%Cloud_Type = Symbol%OPAQUE_ICE_TYPE
   endwhere
   where (Output%Cloud_Type == Symbol%OPAQUE_ICE_TYPE .and. (Output%Ec .ltr. 0.8))
      Output%Cloud_Type = Symbol%CIRRUS_TYPE
   endwhere
   where (Output%Cloud_Type == Symbol%OPAQUE_ICE_TYPE .and. (abs(Output%Tc-Input%Tropopause_Temperature) .ltr. 5.0))
      Output%Cloud_Type = Symbol%OVERSHOOTING_TYPE
   endwhere
   where (Output%Cloud_Type == Symbol%CIRRUS_TYPE .and. (Output%Lower_Tc .ner. MISSING_VALUE_REAL4) .and. &
          MULTI_LAYER_LOGIC_FLAG /= 1)
      Output%Cloud_Type = Symbol%OVERLAP_TYPE
   endwhere
   where (Output%Cloud_Type == Symbol%WATER_TYPE .and. (Output%Tc .ltr. 273.15))
      Output%Cloud_Type = Symbol%SUPERCOOLED_TYPE
   endwhere
   where (Output%Cloud_Type == Symbol%WATER_TYPE .and. (abs(Output%Tc-Input%Surface_Temperature) .ltr. 5.0))
      Output%Cloud_Type = Symbol%FOG_TYPE
   endwhere
   !--- if Lower Tc is close to surface, assume it is not a lower cloud layer
   where (Output%Cloud_Type == Symbol%OVERLAP_TYPE .and. (abs(Output%Lower_Tc-Input%Surface_Temperature) .ltr. 5.0))
      Output%Cloud_Type = Symbol%CIRRUS_TYPE
   endwhere
   !--- if Lower Tc is so uncertain, consider it single layer (AKM(3))
   where (Output%Cloud_Type == Symbol%OVERLAP_TYPE .and. (abs(Output%Lower_Tc-Input%Surface_Temperature) .ltr. 5.0))
      Output%Cloud_Type = Symbol%CIRRUS_TYPE
   endwhere

   !--- homogenous freezing point of water
   where (Output%Cloud_Type == Symbol%SUPERCOOLED_TYPE .and. (Output%Tc .ltr. 233.0))
      Output%Cloud_Type = Symbol%OPAQUE_ICE_TYPE
   endwhere

   !--- force acha cloud type to report clear and prob clear
   where(Input%Cloud_Mask == Symbol%CLEAR)
     Output%Cloud_Type = Symbol%CLEAR_TYPE
   endwhere 

   where(Input%Cloud_Mask == Symbol%PROB_CLEAR)
     Output%Cloud_Type = Symbol%PROB_CLEAR_TYPE
   endwhere 

   !--- make sure lower cloud is consistent with cloud type
   where(Output%Cloud_Type /= Symbol%OVERLAP_TYPE)
     Output%Lower_Tc = MISSING_VALUE_REAL4
     Output%Lower_Pc = MISSING_VALUE_REAL4
     Output%Lower_Zc = MISSING_VALUE_REAL4
   endwhere
   
end subroutine DETERMINE_ACHA_CLOUD_TYPE
!---------------------------------------------------------------
! Test the input to Check to see if this pixel should be skipped
! if Bad_Input_Flag = .true., acha won't do a retrieval
!---------------------------------------------------------------
subroutine TEST_INPUT(Acha_Mode_Flag, ABI_Use_104um_Flag, Input, Elem_Idx,Line_Idx, Bad_Input_Flag) 

    character(len=*), intent(in):: Acha_Mode_Flag
    logical, intent(in):: ABI_Use_104um_Flag
    type(acha_input_struct), intent(in) :: Input
    integer, intent(in):: Elem_Idx, Line_Idx
    logical, intent(out):: Bad_Input_Flag

    Bad_Input_Flag = .false.

    if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == Symbol%YES) then
        Bad_Input_Flag = .true.
        return
    endif

    if (Input%Sensor_Zenith_Angle(Elem_Idx,Line_Idx) > Sensor_Zenith_Threshold) then
        Bad_Input_Flag = .true.
        return
    endif

    if ((Input%Cloud_Mask(Elem_Idx,Line_Idx) == Symbol%CLEAR) .or.  &
        (Input%Cloud_Mask(Elem_Idx,Line_Idx) == Symbol%PROB_CLEAR)) then 
         Bad_Input_Flag = .true.
         return
    endif

    if ((Input%Bt_110um(Elem_Idx,Line_Idx) < 170.0) .or. &         !begin data check
        (Input%Bt_110um(Elem_Idx,Line_Idx) > 340.0) .or. &
        (Input%Surface_Temperature(Elem_Idx,Line_Idx) < 180.0) .or. &
        (Input%Surface_Temperature(Elem_Idx,Line_Idx) > 340.0) .or. &
        (Input%Tropopause_Temperature(Elem_Idx,Line_Idx) < 160.0) .or. &
        (Input%Tropopause_Temperature(Elem_Idx,Line_Idx) > 270.0)) then
         Bad_Input_Flag = .true.
         return
    endif

    !--- check for missing values for relevant channels
    if (ABI_Use_104um_Flag) then
      if (Input%Bt_104um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    else
      if (Input%Bt_110um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif

    if (index(Acha_Mode_Flag,'038') > 0) then
        if (Input%Bt_038um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'062') > 0) then
        if (Input%Bt_062um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'067') > 0) then
        if (Input%Bt_067um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'073') > 0) then
        if (Input%Bt_073um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'085') > 0) then
        if (Input%Bt_085um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'097') > 0) then
        if (Input%Bt_097um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'104') > 0) then
        if (Input%Bt_104um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'12') > 0) then
        if (Input%Bt_120um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'133') > 0) then
        if (Input%Bt_133um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'136') > 0) then
        if (Input%Bt_136um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'139') > 0) then
        if (Input%Bt_139um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (index(Acha_Mode_Flag,'142') > 0) then
        if (Input%Bt_142um(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
end subroutine TEST_INPUT

!----------------------------------------------------------------------
! Set the Output Packed Qf
!----------------------------------------------------------------------
subroutine SET_OUTPUT_PACKED_QF(Output,Elem_Idx,Line_Idx)
 type(acha_output_struct), intent(inout) :: Output
 integer, intent(in):: Elem_Idx,Line_Idx

 !--- bit1  
 Output%Packed_Qf(Elem_Idx,Line_Idx) =  1_int1

 !--- bit2
 if (Output%OE_Qf(1,Elem_Idx,Line_Idx)  /= CTH_PARAM_FAILED_RETREVIAL)  then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 2_int1
 endif

 !--- bit3
 if (Output%OE_Qf(2,Elem_Idx,Line_Idx)  /= CTH_PARAM_FAILED_RETREVIAL) then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 4_int1
 endif

 !--- bit4
 if (Output%OE_Qf(3,Elem_Idx,Line_Idx)  /= CTH_PARAM_FAILED_RETREVIAL) then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 8_int1
 endif

 !--- bit5
 if (Output%OE_Qf(1,Elem_Idx,Line_Idx)  == CTH_PARAM_LOW_QUALITY_RETREVIAL) then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 16_int1
 endif

 !--- bit6
 if (Output%OE_Qf(2,Elem_Idx,Line_Idx)  == CTH_PARAM_LOW_QUALITY_RETREVIAL) then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 32_int1
 endif

 !--- bit7
 if (Output%OE_Qf(3,Elem_Idx,Line_Idx)  == CTH_PARAM_LOW_QUALITY_RETREVIAL) then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 64_int1
 endif

end subroutine SET_OUTPUT_PACKED_QF

!----------------------------------------------------------------------
! Spectral Cloud Emissivity
!----------------------------------------------------------------------
subroutine COMPUTE_SPECTRAL_CLOUD_EMISSIVITY(Input,Symbol_In,Output)

 type(acha_input_struct), intent(inout) :: Input
 type(acha_symbol_struct), intent(inout) :: Symbol_In
 type(acha_output_struct), intent(inout) :: Output
 type(acha_rtm_nwp_struct) :: ACHA_RTM_NWP
 integer:: Elem_Idx,Line_Idx

 real:: Tc_Temp, Zc_Temp
 integer:: Lev_Idx, RTM_NWP_Error_Flag

 Line_loop: do Line_Idx = 1,Input%Number_of_Lines
 Element_Loop:   do Elem_Idx = 1, Input%Number_of_Elements

    !--- filter
    if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == 1) cycle

    if (Output%Pc(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) cycle

    !--- get profiles for this pixel
    call ACHA_FETCH_PIXEL_NWP_RTM(Input, Symbol_In, &
                                  Elem_Idx,Line_Idx, &
                                  ACHA_RTM_NWP,RTM_NWP_Error_Flag)
    if (RTM_NWP_Error_Flag /= 0) cycle

    call KNOWING_P_COMPUTE_T_Z(ACHA_RTM_NWP,Output%Pc(Elem_Idx,Line_Idx),Tc_Temp,Zc_Temp,Lev_Idx)

    if (Lev_Idx > 0) then

     if (Input%Chan_On_067um == Symbol_In%YES) then
        Output%Ec_067um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_067um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_067um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_067um)
     endif

     if (Input%Chan_On_085um == Symbol_In%YES) then
        Output%Ec_085um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_085um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_085um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_085um)
     endif

     if (Input%Chan_On_097um == Symbol_In%YES) then
        Output%Ec_097um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_097um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_097um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_097um)
     endif

     if (Input%Chan_On_104um == Symbol_In%YES) then
        Output%Ec_104um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_104um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_104um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_104um)
     endif

     if (Input%Chan_On_110um == Symbol_In%YES) then
        Output%Ec_110um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_110um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_110um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_110um)
     endif

     if (Input%Chan_On_120um == Symbol_In%YES) then
        Output%Ec_120um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_120um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_120um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_120um)
     endif

     if (Input%Chan_On_133um == Symbol_In%YES) then
        Output%Ec_133um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_133um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_133um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_133um)
     endif

     if (Input%Chan_On_136um == Symbol_In%YES) then
        Output%Ec_136um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_136um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_136um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_136um)
     endif

     if (Input%Chan_On_139um == Symbol_In%YES) then
        Output%Ec_139um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_139um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_139um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_139um)
     endif

     if (Input%Chan_On_142um == Symbol_In%YES) then
        Output%Ec_142um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_142um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_142um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_142um)
     endif

   endif

 enddo Element_Loop
 enddo Line_loop

end subroutine COMPUTE_SPECTRAL_CLOUD_EMISSIVITY

!--------------------------------------------------
! assign surface emissivity for non-overlap type
!--------------------------------------------------
subroutine SET_SURFACE_EMISSIVITY(Input, Cloud_Type, Elem_Idx,Line_Idx, &
                                  Emiss_Sfc_038um, Emiss_Sfc_062um,&
                                  Emiss_Sfc_067um, Emiss_Sfc_073um, Emiss_Sfc_085um, &
                                  Emiss_Sfc_097um, Emiss_Sfc_104um, Emiss_Sfc_110um, &
                                  Emiss_Sfc_120um, Emiss_Sfc_133um, Emiss_Sfc_136um, &
                                  Emiss_Sfc_139um, Emiss_Sfc_142um)

  type(acha_input_struct), intent(in) :: Input
  integer (kind=int1), intent(in):: Cloud_Type
  integer, intent(in):: Elem_Idx,Line_Idx
  real, intent(out):: Emiss_Sfc_038um,Emiss_Sfc_062um,Emiss_Sfc_067um,Emiss_Sfc_073um,&
                      Emiss_Sfc_085um,Emiss_Sfc_097um,Emiss_Sfc_104um,Emiss_Sfc_110um,&
                      Emiss_Sfc_120um,Emiss_Sfc_133um,Emiss_Sfc_136um,Emiss_Sfc_139um,Emiss_Sfc_142um

  Emiss_Sfc_038um = 1.0
  Emiss_Sfc_062um = 1.0
  Emiss_Sfc_067um = 1.0
  Emiss_Sfc_073um = 1.0
  Emiss_Sfc_085um = 1.0
  Emiss_Sfc_097um = 1.0
  Emiss_Sfc_104um = 1.0
  Emiss_Sfc_110um = 1.0
  Emiss_Sfc_120um = 1.0
  Emiss_Sfc_133um = 1.0
  Emiss_Sfc_136um = 1.0
  Emiss_Sfc_139um = 1.0
  Emiss_Sfc_142um = 1.0
  if (Cloud_Type /= Symbol%OVERLAP_TYPE) then
     if (Input%Chan_On_038um == Symbol%YES) Emiss_Sfc_038um = Input%Surface_Emissivity_038um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_062um == Symbol%YES) Emiss_Sfc_062um = Input%Surface_Emissivity_062um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_067um == Symbol%YES) Emiss_Sfc_067um = Input%Surface_Emissivity_067um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_073um == Symbol%YES) Emiss_Sfc_073um = Input%Surface_Emissivity_073um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_085um == Symbol%YES) Emiss_Sfc_085um = Input%Surface_Emissivity_085um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_097um == Symbol%YES) Emiss_Sfc_097um = Input%Surface_Emissivity_097um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_104um == Symbol%YES) Emiss_Sfc_104um = Input%Surface_Emissivity_104um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_110um == Symbol%YES) Emiss_Sfc_110um = Input%Surface_Emissivity_110um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_120um == Symbol%YES) Emiss_Sfc_120um = Input%Surface_Emissivity_120um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_133um == Symbol%YES) Emiss_Sfc_133um = Input%Surface_Emissivity_133um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_136um == Symbol%YES) Emiss_Sfc_136um = Input%Surface_Emissivity_136um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_139um == Symbol%YES) Emiss_Sfc_139um = Input%Surface_Emissivity_139um(Elem_Idx,Line_Idx)
  endif
end subroutine SET_SURFACE_EMISSIVITY

!----------------------------------------------------------------------------------------
! if 11 micron mode chosen, slot opaque retrieval into ACHA retrieval
!----------------------------------------------------------------------------------------
subroutine APPLY_OPAQUE_RETRIEVAL(Input,Symbol,Elem_Idx,Line_Idx,Output)

    type(acha_input_struct), intent(in) :: Input
    type(acha_symbol_struct), intent(in) :: Symbol
    integer, intent(in):: Elem_Idx, Line_Idx
    type(acha_output_struct), intent(out) :: Output

        if ((Input%Cloud_Mask(Elem_Idx,Line_Idx) == Symbol%CLEAR) .or.  &
            (Input%Cloud_Mask(Elem_Idx,Line_Idx) == Symbol%PROB_CLEAR)) then
          Output%Tc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
          Output%Pc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
          Output%Zc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
          Output%Ec(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
          Output%Beta(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        else
          Output%Tc(Elem_Idx,Line_Idx) = Input%Tc_Opaque(Elem_Idx,Line_Idx)
          Output%Pc(Elem_Idx,Line_Idx) = Input%Pc_Opaque(Elem_Idx,Line_Idx)
          Output%Zc(Elem_Idx,Line_Idx) = Input%Zc_Opaque(Elem_Idx,Line_Idx)
          Output%Ec(Elem_Idx,Line_Idx) = 1.0
          Output%Beta(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        endif
        Output%Tc_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        Output%Pc_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        Output%Zc_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        Output%Ec_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        Output%Beta_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        Output%Packed_Qf(Elem_Idx,Line_Idx) =  0_int1
        Output%Packed_Meta_Data(Elem_Idx,Line_Idx) =  0_int1

end subroutine APPLY_OPAQUE_RETRIEVAL
!----------------------------------------------------------------------------------------
! derive the cloud phase from the type and the ice probability
!----------------------------------------------------------------------------------------
subroutine COMPUTE_PHASE(Symbol, &
                         Ice_Probability_Ap, Ice_Probability_Ap_Uncertainty, &
                         Dump_Diag, Lun_Iter_Dump, Cloud_Phase)

   type(acha_symbol_struct), intent(in) :: Symbol
   logical, intent(in):: Dump_Diag
   real, intent(in):: Ice_Probability_Ap
   real, intent(in):: Ice_Probability_Ap_Uncertainty
   integer, intent(in):: Lun_Iter_Dump
   integer, intent(out):: Cloud_Phase

   Cloud_Phase = Symbol%UNKNOWN_PHASE

   if (Ice_Probability_Ap .ger. 0.5) then
       Cloud_Phase = Symbol%ICE_PHASE
   else
       Cloud_Phase = Symbol%WATER_PHASE
   endif

   if (Dump_Diag) then 
       write(unit=Lun_Iter_Dump,fmt=*) "Ice Probability Ap = ", Ice_Probability_Ap
       write(unit=Lun_Iter_Dump,fmt=*) "Ice Probability Ap Uncer = ", Ice_Probability_Ap_Uncertainty
       write(unit=Lun_Iter_Dump,fmt=*) "Cloud Phase = ", Cloud_Phase
   endif

end subroutine COMPUTE_PHASE

!----------------------------------------------------------------------
!--- now compute Sa
!----------------------------------------------------------------------
subroutine COMPUTE_SA(Input,Symbol, Num_Param, &
                      Elem_Idx,Line_Idx,x_Ap,Sa,Sa_Inv, &
                      Singular_Warning_First_Time, &
                      Dump_Diag, Lun_Iter_Dump,Prior_Success_Flag,Fail_Flag)

  type(acha_input_struct), intent(in) :: Input
  type(acha_symbol_struct), intent(in) :: Symbol
  integer, intent(in):: Num_Param 
  integer, intent(in):: Elem_Idx,Line_Idx
  real, intent(out), dimension(:):: x_Ap
  real, intent(out), dimension(:,:):: Sa, Sa_Inv
  logical, intent(in):: Dump_Diag
  integer, intent(in):: Lun_Iter_Dump
  logical, intent(out):: Prior_Success_Flag
  integer:: Singular_Flag
  logical, intent(inout):: Singular_Warning_First_Time
  integer, intent(inout):: Fail_Flag

  !--- initialize
  x_ap(1) = Input%Tc_Ap(Elem_Idx,Line_Idx)
  x_ap(2) = Input%Ec_Ap(Elem_Idx,Line_Idx)
  x_ap(3) = Input%Beta_Ap(Elem_Idx,Line_Idx)
  x_ap(4) = Input%Lower_Tc_Ap(Elem_Idx,Line_Idx)
  x_ap(5) = Input%Ice_Prob_Ap(Elem_Idx,Line_Idx)

  Prior_Success_Flag = .true.
  Sa = 0.0
  Sa(1,1) = Input%Tc_Ap_Uncer(Elem_Idx,Line_Idx)
  Sa(2,2) = Input%Ec_Ap_Uncer(Elem_Idx,Line_Idx)
  Sa(3,3) = Input%Beta_Ap_Uncer(Elem_Idx,Line_Idx)
  Sa(4,4) = Input%Lower_Tc_Ap_Uncer(Elem_Idx,Line_Idx)
  Sa(5,5) = Input%Ice_Prob_Ap_Uncer(Elem_Idx,Line_Idx) 
  
  !--- square the individual elements to convert to variances (not a matmul)
  Sa = Sa**2

  if (Dump_Diag) then 
     write(unit=Lun_Iter_Dump,fmt=*) "x_Ap = ", x_Ap
     write(unit=Lun_Iter_Dump,fmt=*) "S_Ap = ", Sa(1,1),Sa(2,2),Sa(3,3),Sa(4,4),Sa(5,5)
  endif

  !--- compute inverse of Sa matrix
  Singular_Flag =  INVERT_MATRIX(Sa, Sa_Inv, Num_Param)
  if (Singular_Flag == 1 .and. Singular_Warning_First_Time ) then
    print *, "Cloud Height warning ==> Singular Sa in ACHA", &
           Elem_Idx,Line_Idx
    print*,'Sa: ', Sa(1,1),Sa(2,2),Sa(3,3),Sa(4,4),Sa(5,5)
    Fail_Flag = Symbol%YES
    Singular_Warning_First_Time = .false.
    Prior_Success_Flag = .false.
  endif

end subroutine COMPUTE_SA

!----------------------------------------------------------------------
! Convert the cloud temperature to height and Pressure
!
! Purpose:  Perform the computation of Presure and Height from Temperature
!
! Input: Output%Tc
!
! Output: Output%Zc and Output%Zc
!----------------------------------------------------------------------
subroutine CONVERT_TC_TO_PC_AND_ZC(Input,Symbol_In,Cloud_Type,Output)
  type(acha_input_struct), intent(inout) :: Input
  type(acha_symbol_struct), intent(in) :: Symbol_In
  integer (kind=int1), intent(in), dimension(:,:):: Cloud_Type
  type(acha_output_struct), intent(inout) :: Output
  type(acha_rtm_nwp_struct) :: ACHA_RTM_NWP
  integer:: Elem_Idx, Line_Idx,Dummy_Flag
  real:: T_Tropo, Z_Tropo, P_Tropo
  integer:: Lev_Idx,Ierror,NWP_Profile_Inversion_Flag
  integer:: RTM_NWP_Error_Flag

    Line_loop: do Line_Idx = 1,Input%Number_of_Lines 
    Element_Loop:   do Elem_Idx = 1, Input%Number_of_Elements 

    !--- filter
    if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == 1) cycle
    if (Output%Tc(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) cycle

    !--- get profiles for this pixel
    call ACHA_FETCH_PIXEL_NWP_RTM(Input, Symbol, &
                                  Elem_Idx,Line_Idx, &
                                  ACHA_RTM_NWP,RTM_NWP_Error_Flag)
    if (RTM_NWP_Error_Flag /= 0) cycle
                         
    !--- extract tropopause temp, height and pressure
    T_Tropo = Input%Tropopause_Temperature(Elem_Idx,Line_Idx)
    Z_Tropo = Input%Tropopause_Height(Elem_Idx,Line_Idx)
    P_Tropo = Input%Tropopause_Pressure(Elem_Idx,Line_Idx) 
     
    !--- Default
    if (BOTTOM_UP) then
       call KNOWING_T_COMPUTE_P_Z_BOTTOM_UP(ACHA_RTM_NWP,Symbol_In, &
                            Cloud_Type(Elem_Idx,Line_Idx), &
                            Output%Pc(Elem_Idx,Line_Idx), &
                            Output%Tc(Elem_Idx,Line_Idx), &
                            Output%Zc(Elem_Idx,Line_Idx),&
                            T_Tropo, Z_Tropo, P_Tropo,&
                            Lev_Idx,ierror,NWP_Profile_Inversion_Flag)

       !--- Effective
       call KNOWING_T_COMPUTE_P_Z_BOTTOM_UP(ACHA_RTM_NWP,Symbol,Cloud_Type(Elem_Idx,Line_Idx),&
                             Output%Pc_Eff(Elem_Idx,Line_Idx), &
                             Output%Tc_Eff(Elem_Idx,Line_Idx), &
                             Output%Zc_Eff(Elem_Idx,Line_Idx),&
                             T_Tropo, Z_Tropo, P_Tropo,&
                             Lev_Idx,ierror,Dummy_Flag)
       !--- Lower
       call KNOWING_T_COMPUTE_P_Z_BOTTOM_UP(ACHA_RTM_NWP,Symbol,Cloud_Type(Elem_Idx,Line_Idx),&
                             Output%Lower_Pc(Elem_Idx,Line_Idx), &
                             Output%Lower_Tc(Elem_Idx,Line_Idx), &
                             Output%Lower_Zc(Elem_Idx,Line_Idx),&
                             T_Tropo, Z_Tropo, P_Tropo,&
                             Lev_Idx,ierror,Dummy_Flag)

    else
       call KNOWING_T_COMPUTE_P_Z(ACHA_RTM_NWP,Symbol_In, &
                            Cloud_Type(Elem_Idx,Line_Idx), &
                            Output%Pc(Elem_Idx,Line_Idx), &
                            Output%Tc(Elem_Idx,Line_Idx), &
                            Output%Zc(Elem_Idx,Line_Idx),&
                            T_Tropo, Z_Tropo, P_Tropo,&
                            Lev_Idx,ierror,NWP_Profile_Inversion_Flag)

       !--- Effective
       call KNOWING_T_COMPUTE_P_Z(ACHA_RTM_NWP,Symbol,Cloud_Type(Elem_Idx,Line_Idx),&
                            Output%Pc_Eff(Elem_Idx,Line_Idx), &
                            Output%Tc_Eff(Elem_Idx,Line_Idx), &
                            Output%Zc_Eff(Elem_Idx,Line_Idx),&
                            T_Tropo, Z_Tropo, P_Tropo,&
                            Lev_Idx,ierror,Dummy_Flag)

       !--- Lower
       call KNOWING_T_COMPUTE_P_Z(ACHA_RTM_NWP,Symbol,Cloud_Type(Elem_Idx,Line_Idx),&
                             Output%Lower_Pc(Elem_Idx,Line_Idx), &
                             Output%Lower_Tc(Elem_Idx,Line_Idx), &
                             Output%Lower_Zc(Elem_Idx,Line_Idx),&
                             T_Tropo, Z_Tropo, P_Tropo,&
                             Lev_Idx,ierror,Dummy_Flag)

    endif

!   if (Elem_Idx == 140 .and. Line_Idx == 100) then
!      write(unit=6,fmt=*) "BOTTOM_UP = ", BOTTOM_UP
!      write(unit=6,fmt=*) "Profile Interp Zc = ", Output%Zc(Elem_Idx,Line_Idx)
!   endif

    !--- Eff for  low clouds over water, force fixed lapse rate estimate of height
    if ((USE_LAPSE_RATE_WATER_FLAG .and. Input%Surface_Type(Elem_Idx,Line_Idx) == Symbol%WATER_SFC .or. &
         USE_LAPSE_RATE_LAND_FLAG .and. Input%Surface_Type(Elem_Idx,Line_Idx) /= Symbol%WATER_SFC) .and. &
         Input%Snow_Class(Elem_Idx,Line_Idx) == Symbol%NO_SNOW) then

       call COMPUTE_HEIGHT_FROM_LAPSE_RATE(ACHA_RTM_NWP, &
                                     Input%Snow_Class(Elem_Idx,Line_Idx), &
                                     Input%Surface_Type(Elem_Idx,Line_Idx), &
                                     Cloud_Type(Elem_Idx,Line_Idx), &
                                     Input%Surface_Temperature(Elem_Idx,Line_Idx), &
                                     Input%Surface_Elevation(Elem_Idx,Line_Idx), &
                                     MAX_DELTA_T_INVERSION, &
                                     Output%Tc_Eff(Elem_Idx,Line_Idx), &
                                     Output%Zc_Eff(Elem_Idx,Line_Idx), &
                                     Output%Pc_Eff(Elem_Idx,Line_Idx), &
                                     Output%Inversion_Flag(Elem_Idx,Line_Idx))
    endif

    !---  for low clouds over water, force fixed lapse rate estimate of height
!   if (Elem_Idx == 140 .and. Line_Idx == 100) then
!      write(unit=6,fmt=*) "Lapse Rate Flags = ", USE_LAPSE_RATE_WATER_FLAG,  &
!                                                             USE_LAPSE_RATE_LAND_FLAG 
!   endif
    if ((USE_LAPSE_RATE_WATER_FLAG .and. Input%Surface_Type(Elem_Idx,Line_Idx) == Symbol%WATER_SFC) .or. &
        (USE_LAPSE_RATE_LAND_FLAG .and. Input%Surface_Type(Elem_Idx,Line_Idx) /= Symbol%WATER_SFC)) then

       call COMPUTE_HEIGHT_FROM_LAPSE_RATE(ACHA_RTM_NWP, &
                                     Input%Snow_Class(Elem_Idx,Line_Idx), &
                                     Input%Surface_Type(Elem_Idx,Line_Idx), &
                                     Cloud_Type(Elem_Idx,Line_Idx), &
                                     Input%Surface_Temperature(Elem_Idx,Line_Idx), &
                                     Input%Surface_Elevation(Elem_Idx,Line_Idx), &
                                     MAX_DELTA_T_INVERSION, &
                                     Output%Tc(Elem_Idx,Line_Idx), &
                                     Output%Zc(Elem_Idx,Line_Idx), &
                                     Output%Pc(Elem_Idx,Line_Idx), &
                                     Output%Inversion_Flag(Elem_Idx,Line_Idx))
!   if (Elem_Idx == 140 .and. Line_Idx == 100) then
!      write(unit=6,fmt=*) "Lapse Rate Zc = ", Output%Zc(Elem_Idx,Line_Idx)
!   endif

    endif

    !--  If Lower Cloud is placed at surface - assume this single layer
    if (MULTI_LAYER_LOGIC_FLAG == 0 .or. MULTI_LAYER_LOGIC_FLAG == 2) then   !??????????
       if (Output%Lower_Zc(Elem_Idx,Line_Idx) < 1000.0) then
!         Output%Lower_Tc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
!         Output%Lower_Pc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
!         Output%Lower_Zc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
       endif
    endif

    !-- Compute Height Uncertainty
    Output%Zc_Uncertainty(Elem_Idx,Line_Idx) = Output%Tc_Uncertainty(Elem_Idx,Line_Idx) /  &
                                             ABS_LAPSE_RATE_DT_DZ_UNCER

    if (Output%Lower_Tc_Uncertainty(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) then
      Output%Lower_Zc_Uncertainty(Elem_Idx,Line_Idx) = Output%Lower_Tc_Uncertainty(Elem_Idx,Line_Idx) /  &
                                                       ABS_LAPSE_RATE_DT_DZ_UNCER
    endif

    !-- Compute Pressure Uncertainty
    Output%Pc_Uncertainty(Elem_Idx,Line_Idx) = Output%Zc_Uncertainty(Elem_Idx,Line_Idx) *  &
                                            ABS_LAPSE_RATE_DlnP_DZ_UNCER * Output%Pc(Elem_Idx,Line_Idx)

    Output%LOWER_Pc_Uncertainty(Elem_Idx,Line_Idx) = Output%LOWER_Zc_Uncertainty(Elem_Idx,Line_Idx) *  &
                                            ABS_LAPSE_RATE_DlnP_DZ_UNCER * Output%LOWER_Pc(Elem_Idx,Line_Idx)

!   print *, 'Test ', Output%Zc_Uncertainty(Elem_Idx,Line_Idx), ABS_LAPSE_RATE_DlnP_DZ_UNCER, Output%Pc(Elem_Idx,Line_Idx)

    enddo Element_Loop
    enddo Line_Loop


end subroutine CONVERT_TC_TO_PC_AND_ZC
!--------------------------------------------------------------------------------------------------
! save the OE output for a pixel to the output structure.
!--------------------------------------------------------------------------------------------------
subroutine SAVE_X_2_OUTPUT(Elem_Idx,Line_Idx,Symbol,Cloud_Type,Fail_Flag,x,x_ap,Sa,Sx,AKM,Meta_Data_Flags,Output)
  integer, intent(in):: Elem_Idx, Line_Idx, Fail_Flag
  type(acha_symbol_struct), intent(in) :: Symbol
  integer (kind=int1), intent(in):: Cloud_Type
  real, intent(in), dimension(:):: x, x_ap
  real, intent(in), dimension(:,:):: Sa
  real, intent(in), dimension(:,:):: Sx
  real, intent(in), dimension(:,:):: AKM
  integer(kind=int1), intent(in), dimension(:):: Meta_Data_Flags
  type(acha_output_struct), intent(inout) :: Output
  integer:: Param_Idx, i


   !--- Successful Retrieval Post Processing
   if (Fail_Flag == Symbol%NO) then  !successful retrieval if statement

     !--- save retrievals into the output variables
     Output%Tc(Elem_Idx,Line_Idx) = x(1)
     Output%Ec(Elem_Idx,Line_Idx) = x(2)   !note, this is slant
     Output%Beta(Elem_Idx,Line_Idx) = x(3)

     if (Cloud_Type == Symbol%OVERLAP_TYPE) then
        !--- if no confidence, remove lower cloud since it might be fictuous
        if (AKM(4,4) > 0.5) then  !HOLD ON - THIS DOES FIX INFLUENCE ON SOLUTION
          Output%Lower_Tc(Elem_Idx,Line_Idx) = x(4)
        else
         Output%Lower_Tc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        endif
     endif
     Output%Ice_Probability(Elem_Idx,Line_Idx) = x(5)

     !--- save uncertainty estimates
     Output%Tc_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(1,1))
     Output%Ec_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(2,2))
     Output%Beta_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(3,3))
     Output%Lower_Tc_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(4,4))   
     Output%Ice_Probability_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(5,5))   

     !--- set quality flag for a successful retrieval
     Output%Qf(Elem_Idx,Line_Idx) = CTH_DQF_GOOD_RETREVIAL

     !-----------------------------------------------------------------------------
     !--- quality flags of the retrieved parameters
     !-----------------------------------------------------------------------------
     do Param_Idx = 1,Num_Param    !loop over parameters
        if (Sx(Param_Idx,Param_Idx) < 0.111*Sa(Param_Idx,Param_Idx) ) THEN
             Output%OE_Qf(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_1_3_APRIORI_RETREVIAL
        elseif (Sx(Param_Idx,Param_Idx) < 0.444*Sa(Param_Idx,Param_Idx)) THEN
             Output%OE_Qf(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_2_3_APRIORI_RETREVIAL
        else
             Output%OE_Qf(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_LOW_QUALITY_RETREVIAL
        endif

        !--- set quality flag for a marginal retrieval
        if (Output%OE_Qf(Param_Idx,Elem_Idx,Line_Idx) == CTH_PARAM_LOW_QUALITY_RETREVIAL) then
           Output%Qf(Elem_Idx,Line_Idx) = CTH_DQF_MARGINAL_RETREVIAL
        endif

     enddo

  else

   !--- failed
   if (SET_FAILED_TO_MISSING) then 
      Output%Tc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Output%Ec(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Output%Beta(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Output%Lower_Tc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Output%Ice_Probability(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
   else
      Output%Tc(Elem_Idx,Line_Idx) = x_Ap(1)   
      Output%Ec(Elem_Idx,Line_Idx) = x_Ap(2)   
      Output%Beta(Elem_Idx,Line_Idx) = x_Ap(3) 
      Output%Lower_Tc(Elem_Idx,Line_Idx) = x_Ap(4)
      Output%Ice_Probability(Elem_Idx,Line_Idx) = x_Ap(5)
   endif

   Output%Qf(Elem_Idx,Line_Idx) = CTH_DQF_RETREVIAL_ATTEMPTED

  endif

  !--- Pack Meta Data for Output
  do i = 1, 8
     Output%Packed_Meta_Data(Elem_Idx,Line_Idx) = Output%Packed_Meta_Data(Elem_Idx,Line_Idx) +  &
                                                  int((2**(i-1)) * Meta_Data_Flags(i),kind=int1)
  enddo

end subroutine SAVE_X_2_OUTPUT

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module AWG_CLOUD_HEIGHT
