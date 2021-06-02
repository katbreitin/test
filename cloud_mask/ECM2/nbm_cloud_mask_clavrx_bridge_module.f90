!$Id: nbm_cloud_mask_clavrx_bridge_module.f90 3621 2019-12-09 21:18:24Z wstraka $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: naive_bayesian_clavrx_bridge_module.f90 (src)
!       naive_bayesian_clavrx_bridge_module (program)
!
! PURPOSE: 
!
! DESCRIPTION: 
!
! AUTHORS:
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
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
! Note, to use the diagnostic variables, do this
!   - set the USE_DIAG flag to true
!   - turn on the Diag argument to the desirefd routine
!   - in the desired routine, set the diag variables to what you want
!   - when done, repeat this in reverse
!
!--------------------------------------------------------------------------------------
module NBM_CLOUD_MASK_CLAVRX_BRIDGE

   ! -- MODULES USED
   use PIXEL_COMMON_MOD, only: &
       Ch, &
       Geo, & 
       Nav, &
       Sfc, &
       NWP_PIX, &
       Zc_Opaque_Cloud, &
       Tc_Opaque_Cloud, &
       Image, &
       CLDMASK, &
       Sensor, &
       Ancil_Data_Dir, &
       Bayesian_Cloud_Mask_Name, &
       Bad_Pixel_Mask, &
       Cld_Phase, &
       Cld_Phase_Uncertainty, &
       Cld_Type, &
       SST_Anal_Uni, &
       Solar_Contamination_Mask, &
       Bt_11um_Sounder, &
       AVHRR_IFF_Flag, &
       Ems_Ch20_Median_3x3, &
       Covar_Ch27_Ch31_5x5, &
       Covar_Ch27_Ch38_5x5, &
       Ref_ChDNB_Lunar_Std_3x3, &
       Ref_ChDNB_Lunar_Min_3x3, &
       Month, &
       Diag_Pix_Array_1, &
       Diag_Pix_Array_2, &
       Diag_Pix_Array_3, &
       ABI_Use_104um_Flag, &
       Beta_11um_12um_Tropo_Rtm, &
       Beta_11um_133um_Tropo_Rtm, &
       Beta_104um_12um_Tropo_Rtm, &
       Beta_104um_133um_Tropo_Rtm, &
       i_LRC, &
       j_LRC, &
       Use_Aux_Flag

   use CONSTANTS_MOD, only: &
       Sym, &
       Cloud_Mask_Version, & 
       Cloud_Mask_Lut_Version

   use CX_REAL_BOOLEAN_MOD
       
   use NBM_CLOUD_MASK_MODULE
   use NB_CLOUD_MASK_ADDONS
   use NB_CLOUD_MASK_SERVICES
   use NB_CLOUD_MASK_SOLAR_RTM
   use NBM_CLOUD_MASK_LUT_MODULE

   use CLAVRX_MESSAGE_MOD, only: MESG

   ! New LHP stuff. Need to have Calibration constants for initialization
   use CX_ABI_LHP_MOD_ALG_FUNCTION

   use CALIBRATION_CONSTANTS_MOD,only: &
      ABI_FPT_Thresh_038um, ABI_FPT_Thresh_062um, ABI_FPT_Thresh_067um &
    , ABI_FPT_Thresh_073um, ABI_FPT_Thresh_085um, ABI_FPT_Thresh_097um &
    , ABI_FPT_Thresh_104um, ABI_FPT_Thresh_110um, ABI_FPT_Thresh_120um &
    , ABI_FPT_Thresh_133um

   implicit none

   public :: NBM_CLOUD_MASK_BRIDGE
   public:: COMPUTE_TYPE_FROM_PHASE

   private :: COVARIANCE_LOCAL
   private :: SET_SYMBOL
   private :: SET_INPUT
   private :: SET_OUTPUT
   private :: SET_DIAG
   private :: MEDIAN_FILTER_MASK
   private :: SET_ECM_CHAN_ON !WCS3
   private :: LHP_CHN_CHECK   !WCS3
   private :: BETA_11_12_OVERLAP_TEST
   private :: BETA_11_133_OVERLAP_TEST
   private :: WATER_EDGE_FILTER

   !--- define these structure as module wide
   type(mask_input), private :: Input   
   type(mask_output), private :: Output   
   type(diag_output), private :: Diag  
   type(symbol_naive_bayesian), private :: Symbol


  !--- string to control on-screen prompts
  character(*), parameter, private :: EXE_PROMPT_CM = "ECM2 Cloud Mask Bridge >> "
  real, dimension(:,:), allocatable, target, private :: Ref1_Clr_Routine
  integer:: N_Class

  !-----------------------------------------------------------------------
  ! flags for options
  !-----------------------------------------------------------------------
  logical, parameter, private:: USE_DIAG = .false.
  logical, parameter, private:: USE_PRIOR_TABLE = .true.
  logical, parameter, private:: USE_065UM_RTM = .false.
  real, parameter, private:: CLD_PHASE_UNCER_LRC_THRESH = 0.10    !NEEDS VERIFICATION
  real, parameter, private:: CLD_PHASE_UNCER_MULTI_THRESH = 0.10    !NEEDS VERIFICATION

  !Use_104um_Flag
  logical, private:: Use_104um_Flag

contains

!----------------------------------------------------------------------
! Bridge Routine
!
! Note, the Diag argument is optional
!---------------------------------------------------------------------- 
 subroutine NBM_CLOUD_MASK_BRIDGE(Segment_Number)
 
   implicit none

   integer, intent(in):: Segment_Number

   integer :: i, j
   character (len=1020):: Naive_Bayes_File_Name_Full_Path
   character (len=1020):: Prior_File_Name_Full_Path
   character (len=1020):: Naive_Bayes_Default_2d_File_Name_Full_Path
   character (len=1020):: Naive_Bayes_Night_2d_File_Name_Full_Path
   character (len=1020):: Naive_Bayes_Day_160_2d_File_Name_Full_Path
   character (len=1020):: Naive_Bayes_Day_375_2d_File_Name_Full_Path

   logical, save:: First_Call = .true.
   integer, parameter:: Nmed = 2
   real(kind=real4):: Nmed_Total
   integer(kind=int1), dimension(:,:), allocatable:: I1_Temp_1
   integer(kind=int1), dimension(:,:), allocatable:: I1_Temp_2
   integer:: Num_Elem
   integer:: Num_Line
   real:: Post_Prob_Ice_Norm
   real:: Post_Prob_Water_Norm
   logical:: Is_Overlap
   real:: Beta_x_12um_Tropo_Rtm_Temp
   real:: Beta_x_133um_Tropo_Rtm_Temp
   integer:: i1,i2,j1,j2
   integer:: N_Box
   integer(kind=int1), dimension(:,:), allocatable:: Phase_Temp
   real, dimension(:,:), allocatable:: Uncer_Temp


   if (First_Call) then
       call MESG('NB-Modified Cloud Mask starts ', color = 46)
   endif

   !--- set structure (symbol, input, output, diag)  
   !    elements to corresponding values in this framework
   call SET_SYMBOL()
   
   !--- allocate internal Ch1 clear sky albedo
   Num_Elem = Image%Number_Of_Elements
   Num_Line = Image%Number_Of_Lines_Read_This_Segment
   
   
   !--- initialize all variables
   
   CLDMASK%Smoke_Mask = symbol%NO
   CLDMASK%Dust_Mask  = symbol%NO
   CLDMASK%Fire_Mask  = symbol%NO
   CLDMASK%Cld_Test_Vector_Packed = Missing_Value_Int1
   CLDMASK%Cld_Mask = Missing_Value_Int1
   CLDMASK%Posterior_Cld_Probability = MISSING_VALUE_REAL4
   !CLDMASK%Cloud_Mask_QF = symbol%FILL_QF !
  


   if (USE_065UM_RTM) allocate(Ref1_Clr_Routine(Num_Elem,Num_Line))

!-------------------------------------------------------------------------------
   !--- on first segment, read table, set 
!-------------------------------------------------------------------------------
   if (.not. Is_Classifiers_Read) then

       ! --- Read LUT
       Naive_Bayes_File_Name_Full_Path = trim(Ancil_Data_Dir)// &
            "static/luts/ecm2/"//trim(Bayesian_Cloud_Mask_Name)
       call NBM_CLOUD_MASK_LUT_READ(trim(Naive_Bayes_File_Name_Full_Path), N_Class)

   endif

   !---  Read and Compute Prior Cloud Probability
   ! - initialize to missing
   CLDMASK%Prior_Cld_Probability = MISSING_VALUE_REAL4
   if (USE_PRIOR_TABLE) then 
     Prior_File_Name_Full_Path = trim(Ancil_Data_Dir)//"static/luts/ecm2/"//"nb_cloud_mask_calipso_prior.nc"
     call NBM_CLOUD_MASK_COMPUTE_PRIOR(trim(Prior_File_Name_Full_Path),Nav%Lon,Nav%Lat,Month,CLDMASK%Prior_Cld_Probability)
   endif

   !--- Compute TOA Clear-Sky 0.65um Reflectance
   if (USE_065UM_RTM .and. Sensor%Chan_On_Flag_Default(1) == sym%YES)  then 
     call  CLEAR_SKY_TOA_RTM_065UM(Bad_Pixel_Mask, &
                                   NWP_PIX%Tpw, &
                                   NWP_PIX%Ozone, &
                                   Geo%Scatangle, &
                                   Geo%Satzen, &
                                   Geo%Solzen, &
                                   ch(1)%Sfc_Ref_White_Sky, &
                                   Sfc%Sfc_Type, &
                                   Sfc%Snow, &
                                   Ref1_Clr_Routine)
   endif

   !-----------  Initialize channel on flag -----   WCS
   CALL SET_ECM_CHAN_ON()

   !-- Since the arrays are already initialized to fill
   !   in process_clavrx by RESET_PIXEL_ARRAYS_TO_MISSING
   !   We simply exit for catastrophic case

   IF ((Input%Chan_On_11um == sym%NO) .AND. &
       (Input%Chan_On_10um == sym%NO)) THEN

       call MESG( "ERROR: Catastrophic LHP Event, exiting" , level = verb_lev % DEFAULT)
       RETURN
   ENDIF
   

   !-----------    loop over pixels -----   
   line_loop: do i = 1, Image%Number_Of_Elements
      elem_loop: do  j = 1, Image%Number_Of_Lines_Read_This_Segment

         call SET_INPUT(i,j)

         ! --- check if pixel is valid
         if (Input%Invalid_Data_Mask == 1) cycle

         !---call cloud mask routine
         call NBM_CLOUD_MASK_ALGORITHM(  &
                      N_Class, &
                      Symbol,  &
                      Input,   &
                      Output,  &
                      USE_PRIOR_TABLE, &
                      Use_104um_Flag)!, & !)
                      !Diag)

         ! --- set cloud mask qf: good=0, bad=1
         Output%Cld_Mask_QF = symbol%BAD_QF
         if (ibits(Output%Cld_Flags_Packed(1),0,1) == 1) Output%Cld_Mask_QF = symbol%GOOD_QF
         !If any channel used is turned off for the FPT, mark QF as degraded
         IF ((Sensor%WMO_Id ==271) .AND. &
             (Input%Chan_On_67um == sym%NO) .OR.  &
             (Input%Chan_On_73um == sym%NO) .OR.  &
             (Input%Chan_On_85um == sym%NO) .OR.  &
             (Input%Chan_On_12um == sym%NO)) Output%Cld_Mask_QF = symbol%DEGRADED_QF

         !--- call non-cloud detection routines (smoke, dust and fire)
         call NB_CLOUD_MASK_ADDONS_ALGORITHM(Symbol,  &
                                          Input, &
                                          Output)!, &!)
                                          !Diag)   !optional

         call SET_OUTPUT(i,j)

         if (USE_DIAG) call SET_DIAG(i,j)
   
         !-----------------------------------------------------------------------
         ! CLAVR-x specific processing
         !-----------------------------------------------------------------------

         !--- unpack elements of the cloud test vector into clavr-x global arrays
         !CLDMASK%Bayes_Mask_Sfc_Type(i,j) = ibits(CLDMASK%Cld_Test_Vector_Packed(3,i,j),0,3)

      end do elem_loop
   end do line_loop

   !------------------------------------------------------------------------------
   ! Apply Median Filters
   !------------------------------------------------------------------------------
   Nmed_Total = (2.0*Nmed+1)**2
   allocate(I1_Temp_1(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
   allocate(I1_Temp_2(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
   I1_Temp_1 = CLDMASK%Dust_Mask
   I1_Temp_2 = CLDMASK%Smoke_Mask
   I1_Temp_1(Nmed:Image%Number_Of_Elements-Nmed,1+Nmed:Image%Number_Of_Lines_Read_This_Segment-Nmed) = 0
   I1_Temp_2(Nmed:Image%Number_Of_Elements-Nmed,1+Nmed:Image%Number_Of_Lines_Read_This_Segment-Nmed) = 0
   line_loop_median: do i = 1+Nmed, Image%Number_Of_Elements-Nmed
      elem_loop_median: do  j = 1+Nmed, Image%Number_Of_Lines_Read_This_Segment-Nmed
         I1_Temp_1(i,j) = nint(count(CLDMASK%Dust_Mask(i-Nmed:i+Nmed,j-Nmed:j+Nmed) == sym%YES) / Nmed_Total)
         I1_Temp_2(i,j) = nint(count(CLDMASK%Smoke_Mask(i-Nmed:i+Nmed,j-Nmed:j+Nmed) == sym%YES) / Nmed_Total)
      end do elem_loop_median
   end do line_loop_median
   CLDMASK%Dust_Mask = I1_Temp_1
   CLDMASK%Smoke_Mask = I1_Temp_2
   deallocate(I1_Temp_1)
   deallocate(I1_Temp_2)

   !-------------------------------------------------------------------------------------------------------
   !--- save dust, smoke, fire to the corresponding bit structure 
   !-------------------------------------------------------------------------------------------------------
   line_loop_pack: do i = 1, Image%Number_Of_Elements
      elem_loop_pack: do  j = 1, Image%Number_Of_Lines_Read_This_Segment
           if (CLDMASK%Smoke_Mask(i,j) == 1) CLDMASK%Cld_Test_Vector_Packed(2,i,j) = ibset(CLDMASK%Cld_Test_Vector_Packed(2,i,j),4)
           if (CLDMASK%Dust_Mask(i,j)  == 1) CLDMASK%Cld_Test_Vector_Packed(2,i,j) = ibset(CLDMASK%Cld_Test_Vector_Packed(2,i,j),5)
           if (CLDMASK%Fire_Mask(i,j)  == 1) CLDMASK%Cld_Test_Vector_Packed(2,i,j) = ibset(CLDMASK%Cld_Test_Vector_Packed(2,i,j),7)
           if (CLDMASK%Thin_Cirr_Mask(i,j) == 1) CLDMASK%Cld_Test_Vector_Packed(3,i,j) = ibset(CLDMASK%Cld_Test_Vector_Packed(3,i,j),3)
      end do elem_loop_pack
   end do line_loop_pack

   !------------------------------------------------------------------------------
   !-- CLAVR-x specific Processing
   !--- grab version tags for output as attributes in level2
   !--- only need to do this once, so do on first segment
   !------------------------------------------------------------------------------
   if (Segment_Number == 1) then
     call SET_NBM_CLOUD_MASK_VERSION(Cloud_Mask_Version)
   endif

   !-------------------------------------------------------------------------------
   ! Compute Cloud Probability Uncertainty
   !-------------------------------------------------------------------------------
   CLDMASK%Posterior_Cld_Probability_Uncer = MISSING_VALUE_REAL4
   where(CLDMASK%Posterior_Cld_Probability >= 0.50)
         CLDMASK%Posterior_Cld_Probability_Uncer = 1.0 - CLDMASK%Posterior_Cld_Probability
   endwhere     

   where((CLDMASK%Posterior_Cld_Probability .ner. MISSING_VALUE_REAL4) .and. &
         (CLDMASK%Posterior_Cld_Probability .ltr. 0.50))
         CLDMASK%Posterior_Cld_Probability_Uncer = CLDMASK%Posterior_Cld_Probability
   endwhere     

   !-------------------------------------------------------------------------------
   ! Compute Cloud Phase, Cloud_Phase Uncertainty
   !-------------------------------------------------------------------------------
   line_loop_phase: do i = 1, Image%Number_Of_Elements
      elem_loop_phase: do  j = 1, Image%Number_Of_Lines_Read_This_Segment

        ! --- if use modawg flag set cloud probability 
        if (Input%Use_Aux_Mask == symbol%USE_AUX_MODAWG) then
          select case(CLDMASK%Cld_Mask_Aux(i,j))
            case(0) ! conf clear
              CLDMASK%Posterior_Cld_Probability(i,j) = 0.01
            case(1) ! prob clear
              CLDMASK%Posterior_Cld_Probability(i,j) = 0.25
            case(2) ! prob cloud
              CLDMASK%Posterior_Cld_Probability(i,j) = 0.75
            case(3) ! conf cloud
              CLDMASK%Posterior_Cld_Probability(i,j) = 1.0
            case default
              CLDMASK%Posterior_Cld_Probability(i,j) = MISSING_VALUE_REAL4
          end select
        endif

        Post_Prob_Water_Norm = CLDMASK%Posterior_Water_Probability(i,j) / CLDMASK%Posterior_Cld_Probability(i,j)
        Post_Prob_Ice_Norm = CLDMASK%Posterior_Ice_Probability(i,j) / CLDMASK%Posterior_Cld_Probability(i,j)
        Cld_Phase(i,j)= sym%CLEAR_PHASE
        if (CLDMASK%Posterior_Cld_Probability(i,j) .ger. 0.50) then
          Cld_Phase(i,j) = sym%WATER_PHASE
          Cld_Phase_Uncertainty(i,j) = (1.0 - Post_Prob_Water_Norm) / CLDMASK%Posterior_Cld_Probability(i,j)
          if (Post_Prob_Ice_Norm .gtr. Post_Prob_Water_Norm) then 
            Cld_Phase(i,j) = sym%ICE_PHASE
            Cld_Phase_Uncertainty(i,j) = (1.0 - Post_Prob_Ice_Norm) / CLDMASK%Posterior_Cld_Probability(i,j)
          endif
!         if (Post_Prob_Ice_Norm .gtr. 0.40) then 
!           Cld_Phase(i,j) = sym%ICE_PHASE
!           Cld_Phase_Uncertainty(i,j) = 1.0 - Post_Prob_Ice_Norm
!         endif
          if ((Cld_Phase(i,j) == sym%WATER_PHASE) .and. (Tc_Opaque_Cloud(i,j) .ltr. 273.0)) then
            Cld_Phase(i,j) = sym%SUPERCOOLED_PHASE
          endif
          !--- prevent any water result colder than homogeneous freezing point  (consider 243)
          if (Tc_Opaque_Cloud(i,j) .ltr. 233.0) then
            Cld_Phase(i,j) = sym%ICE_PHASE
          endif
        endif

      end do elem_loop_phase
   end do line_loop_phase

   !--- LRC Filter for Uncertain Pixels
   line_loop_lrc: do i = 1, Image%Number_Of_Elements
      elem_loop_lrc: do  j = 1, Image%Number_Of_Lines_Read_This_Segment
        !--- LRC Filter
        if ((i_LRC(i,j) > 0) .and. &
            (j_LRC(i,j) > 0) .and. &
            (Cld_Phase_Uncertainty(i,j) .ger. CLD_PHASE_UNCER_LRC_THRESH)) then

            !--- this test prevents clear pixels from replacing cloudy pixels             
            if (Cld_Phase(i,j) > 0 .and. Cld_Phase(i_LRC(i,j),j_LRC(i,j)) > 0 .and. &
                Cld_Phase(i,j)  /= Cld_Phase(i_LRC(i,j),j_LRC(i,j))) then

                Cld_Phase(i,j) = Cld_Phase(i_LRC(i,j),j_LRC(i,j)) 

            endif

            !-- this logic could extend supercooled to water erroneously
            if ((Cld_Phase(i,j) == sym%SUPERCOOLED_PHASE) .and. (Tc_Opaque_Cloud(i,j) .ger. 273.0)) then
               Cld_Phase(i,j) = sym%WATER_PHASE
            endif

        endif

      end do elem_loop_lrc
   end do line_loop_lrc

   !--- try to remove water edges around cirrus that impacts AMVs
   call WATER_EDGE_FILTER(Cld_Phase,2)

   !--- cloud type
   line_loop_type: do i = 1, Image%Number_Of_Elements
      elem_loop_type: do  j = 1, Image%Number_Of_Lines_Read_This_Segment

       Beta_x_12um_Tropo_Rtm_Temp = MISSING_VALUE_REAL4
       Beta_x_133um_Tropo_Rtm_Temp = MISSING_VALUE_REAL4

       if (.not. ABI_Use_104um_Flag) then
       
          if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then 
            if (Sensor%Chan_On_Flag_Default(32) == sym%YES) then 
              Beta_x_12um_Tropo_Rtm_Temp = Beta_11um_12um_Tropo_Rtm(i,j) 
            endif
            if (Sensor%Chan_On_Flag_Default(33) == sym%YES) then 
              Beta_x_133um_Tropo_Rtm_Temp = Beta_11um_133um_Tropo_Rtm(i,j) 
            endif
          endif

        ! --- if use modawg flag set cloud probability 
        if (Input%Use_Aux_Mask == symbol%USE_AUX_MODAWG) then
          call COMPUTE_TYPE_FROM_PHASE(CLDMASK%Cld_Mask_Aux(i,j), &
                                     Cld_Phase(i,j), &
                                     Cld_Phase_Uncertainty(i,j), &
                                     Zc_Opaque_Cloud(i,j), &
                                     Tc_Opaque_Cloud(i,j), &
                                     ch(31)%Emiss_Tropo(i,j), & 
                                     Sensor%Chan_On_Flag_Default(31), & 
                                     Sensor%Chan_On_Flag_Default(32), & 
                                     Sensor%Chan_On_Flag_Default(33), & 
                                     Beta_x_12um_Tropo_Rtm_Temp, &
                                     Beta_x_133um_Tropo_Rtm_Temp, &
                                     Cld_Type(i,j))

        else  
          call COMPUTE_TYPE_FROM_PHASE(CLDMASK%Cld_Mask(i,j), &
                                     Cld_Phase(i,j), &
                                     Cld_Phase_Uncertainty(i,j), &
                                     Zc_Opaque_Cloud(i,j), &
                                     Tc_Opaque_Cloud(i,j), &
                                     ch(31)%Emiss_Tropo(i,j), &
                                     Sensor%Chan_On_Flag_Default(31), &
                                     Sensor%Chan_On_Flag_Default(32), &
                                     Sensor%Chan_On_Flag_Default(33), &
                                     Beta_x_12um_Tropo_Rtm_Temp, &
                                     Beta_x_133um_Tropo_Rtm_Temp, &
                                     Cld_Type(i,j))
        endif

       else

          if (Sensor%Chan_On_Flag_Default(38) == sym%YES) then 
            if (Sensor%Chan_On_Flag_Default(32) == sym%YES) then 
              Beta_x_12um_Tropo_Rtm_Temp = Beta_11um_12um_Tropo_Rtm(i,j) 
            endif
            if (Sensor%Chan_On_Flag_Default(33) == sym%YES) then 
              Beta_x_133um_Tropo_Rtm_Temp = Beta_11um_133um_Tropo_Rtm(i,j) 
            endif
          endif


          call COMPUTE_TYPE_FROM_PHASE(CLDMASK%Cld_Mask(i,j), &
                                     Cld_Phase(i,j), &
                                     Cld_Phase_Uncertainty(i,j), &
                                     Zc_Opaque_Cloud(i,j), &
                                     Tc_Opaque_Cloud(i,j), &
                                     ch(38)%Emiss_Tropo(i,j), & 
                                     Sensor%Chan_On_Flag_Default(38), & 
                                     Sensor%Chan_On_Flag_Default(32), & 
                                     Sensor%Chan_On_Flag_Default(33), & 
                                     Beta_x_12um_Tropo_Rtm_Temp, &
                                     Beta_x_133um_Tropo_Rtm_Temp, &
                                     Cld_Type(i,j))
       endif                              

      enddo elem_loop_type
   enddo line_loop_type

   !-------------------------------------------------------------------------------
   ! on the last segment, Save number of classifiers and names to CLAVR-x, 
   ! and wipe out the lut from memory and reset is_read_flag to no
   !-------------------------------------------------------------------------------
   if (Segment_Number == Input%Num_Segments) then
      do i = 1, N_Class
         if (i == 1) then
            Output%Cld_Mask_Test_Names = trim(Classifier_Names(i))
         elseif (i == N_Class) then ! add 2 spaces to the end to keep the same length for each classifier name
            Output%Cld_Mask_Test_Names = trim(Output%Cld_Mask_Test_Names)//','//trim(Classifier_Names(i))//'  '
         else
            Output%Cld_Mask_Test_Names = trim(Output%Cld_Mask_Test_Names)//','//trim(Classifier_Names(i))
         endif
      enddo
      CLDMASK%N_Classifiers = N_Class
      CLDMASK%Classifiers_Names_Attr = trim(Output%Cld_Mask_Test_Names)

      call RESET_NBM_CLOUD_MASK_LUT()
       
      if (USE_PRIOR_TABLE) call RESET_NBM_CLOUD_MASK_PRIOR_LUT()

   endif
   

   First_Call = .false.
   
   if (allocated(Ref1_Clr_Routine)) deallocate (Ref1_Clr_Routine)

   end subroutine NBM_CLOUD_MASK_BRIDGE

   !====================================================================
   ! Function Name: Covariance_LOCAL
   !
   ! Function:
   !    Compute the Covariance for two mxn arrays
   !
   ! Description: Covariance = E(XY) - E(X)*E(Y)
   !   
   ! Calling Sequence: BT_WV_BT_Window_Covar(Elem_Idx,Line_Idx) = Covariance( &
   !                       sat%bt10(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
   !                       sat%bt14(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
   !                      Array_Width, Array_Hgt)
   !   
   !
   ! Inputs:
   !   Array 1 - the first array (X)
   !   Array 2 - the second array (Y)
   !   Elem_size
   !   Line_size
   !
   ! Outputs: 
   !   Covariance of X and Y
   !
   ! Dependencies:
   !        none
   !
   ! Restrictions:  None
   !
   ! Reference: Standard definition for the Covariance Computation
   !
   !====================================================================
   function COVARIANCE_LOCAL &
        (Array_One,Array_Two,Array_Width,Array_Hght,Invalid_Data_Mask) &
         result(Covar_Array_One_Array_Two)

   real(kind=real4), intent(in), dimension(:,:):: Array_One
   real(kind=real4), intent(in), dimension(:,:):: Array_Two
   integer(kind=INT4), intent(in):: Array_Width
   integer(kind=INT4), intent(in):: Array_Hght
   integer(kind=INT1), intent(in), dimension(:,:):: Invalid_Data_Mask

   real(kind=real8):: Mean_Array_One
   real(kind=real8):: Mean_Array_Two
   real(kind=real8):: Mean_Array_One_x_Array_Two
   real(kind=real8):: Sum_Array_One
   real(kind=real8):: Sum_Array_Two
   real(kind=real8):: Sum_Array_One_x_Array_Two
   real(kind=real4):: Covar_Array_One_Array_Two

   !--- skip computation for pixel arrays with any missing data
   if (sum(Invalid_Data_Mask) > 0) then
      Covar_Array_One_Array_Two = Missing_Value_Real4
      return
   endif

   Sum_Array_One = sum(Array_One)
   Sum_Array_Two = sum(Array_Two)

   Mean_Array_One = Sum_Array_One / (Array_Width*Array_Hght)
   Mean_Array_Two = Sum_Array_Two / (Array_Width*Array_Hght)

   Sum_Array_One_x_Array_Two = sum(Array_One*Array_Two)
   Mean_Array_One_x_Array_Two = Sum_Array_One_x_Array_Two / (Array_Width*Array_Hght)
   
   Covar_Array_One_Array_Two  = Mean_Array_One_x_Array_Two - &
                                Mean_Array_One * Mean_Array_Two 
   
   end function COVARIANCE_LOCAL

   !============================================================================
   ! set symbols
   !============================================================================
   subroutine SET_SYMBOL()

      symbol%CLOUDY = sym%CLOUDY
      symbol%PROB_CLOUDY = sym%PROB_CLOUDY
      symbol%PROB_CLEAR = sym%PROB_CLEAR
      symbol%CLEAR = sym%CLEAR

      symbol%CLEAR_BINARY = sym%CLEAR_BINARY
      symbol%CLOUDY_BINARY = sym%CLOUDY_BINARY

      symbol%GOOD_QF = sym%ECM_GOOD_QF
      symbol%BAD_QF = sym%ECM_BAD_QF
      symbol%SPACE_QF = sym%ECM_SPACE_QF
      symbol%FILL_QF = sym%ECM_FILL_QF
      symbol%DEGRADED_QF = sym%ECM_DEGRADED_QF

      symbol%NO = sym%NO
      symbol%YES = sym%YES

      symbol%WATER_SFC = sym%WATER_SFC
      symbol%EVERGREEN_NEEDLE_SFC = sym%EVERGREEN_NEEDLE_SFC
      symbol%EVERGREEN_BROAD_SFC = sym%EVERGREEN_BROAD_SFC
      symbol%DECIDUOUS_NEEDLE_SFC = sym%DECIDUOUS_NEEDLE_SFC
      symbol%DECIDUOUS_BROAD_SFC = sym%DECIDUOUS_BROAD_SFC
      symbol%MIXED_FORESTS_SFC = sym%MIXED_FORESTS_SFC
      symbol%WOODLANDS_SFC = sym%WOODLANDS_SFC
      symbol%WOODED_GRASS_SFC = sym%WOODED_GRASS_SFC
      symbol%CLOSED_SHRUBS_SFC = sym%CLOSED_SHRUBS_SFC
      symbol%OPEN_SHRUBS_SFC = sym%OPEN_SHRUBS_SFC
      symbol%GRASSES_SFC = sym%GRASSES_SFC
      symbol%CROPLANDS_SFC = sym%CROPLANDS_SFC
      symbol%BARE_SFC = sym%BARE_SFC
      symbol%URBAN_SFC = sym%URBAN_SFC

      symbol%SHALLOW_OCEAN = sym%SHALLOW_OCEAN
      symbol%LAND = sym%LAND
      symbol%COASTLINE = sym%COASTLINE
      symbol%SHALLOW_INLAND_WATER = sym%SHALLOW_INLAND_WATER
      symbol%EPHEMERAL_WATER = sym%EPHEMERAL_WATER
      symbol%DEEP_INLAND_WATER = sym%DEEP_INLAND_WATER
      symbol%MODERATE_OCEAN = sym%MODERATE_OCEAN
      symbol%DEEP_OCEAN = sym%DEEP_OCEAN

      symbol%NO_SNOW = sym%NO_SNOW
      symbol%SEA_ICE = sym%SEA_ICE
      symbol%SNOW = sym%SNOW   

      symbol%NO_AUX = sym%NO_AUX
      symbol%USE_AUX_MODAWG = sym%USE_AUX_MODAWG

   end subroutine SET_SYMBOL

   !============================================================================
   ! set input
   !============================================================================
   subroutine SET_INPUT(i,j)
      integer, intent (in) :: i, j
      Input%Num_Elem = Image%Number_Of_Elements
      Input%Num_Line = Image%Number_Of_Lines_Read_This_Segment
      Input%Num_Line_Max = Image%Number_Of_Lines_Per_Segment
      Input%Num_Segments = Image%Number_Of_Segments
      Input%Invalid_Data_Mask = Bad_Pixel_Mask(i,j)
      Input%Use_Sounder_11um = 0                        !note, off by default
      Input%Snow_Class = Sfc%Snow(i,j)
      Input%Land_Class = Sfc%Land(i,j)
      Input%Oceanic_Glint_Mask = Sfc%Glint_Mask(i,j)
      Input%Glintzen = Geo%Glintzen(i,j)
      Input%Coastal_Mask = Sfc%Coast_Mask(i,j)
      Input%Solzen = Geo%Solzen(i,j)
      Input%Scatzen = Geo%Scatangle(i,j)
      Input%Senzen = Geo%Satzen(i,j)
      Input%Lat = Nav%Lat(i,j)
      Input%Lon = Nav%Lon(i,j)
      Input%Sst_Anal_Uni = Sst_Anal_Uni(i,j)
      Input%Emiss_Sfc_375um = ch(20)%Sfc_Emiss(i,j)    !note, this is on always
      Input%Zsfc = Sfc%Zsfc(i,j)
      Input%Zsfc_Std = Sfc%Zsfc_Std(i,j)
      Input%Solar_Contamination_Mask = Solar_Contamination_Mask(i,j)
      Input%Sfc_Type = Sfc%Sfc_Type(i,j)
      Input%Sfc_Temp = NWP_PIX%Tsfc(i,j)
      Input%Path_Tpw = NWP_PIX%Tpw(i,j) / Geo%Coszen(i,j)
      Input%Topa = Tc_Opaque_Cloud(i,j)
      Input%Zopa = Zc_Opaque_Cloud(i,j)

      Input%City_Mask = Sfc%City_Mask(i,j)
      Input%Moon_Illum_Frac = Geo%Moon_Illum_Frac

      ! - already read in, if not initialized as missing
      Input%Prior = CLDMASK%Prior_Cld_Probability(i,j)

      if (Input%Chan_On_041um == sym%YES)  then 
        Input%Ref_041um = ch(8)%Ref_Toa(i,j)
      endif
      if (Input%Chan_On_063um == sym%YES)  then 
        Input%Ref_063um = ch(1)%Ref_Toa(i,j)
        if (USE_065UM_RTM) then 
           Input%Ref_063um_Clear = Ref1_Clr_Routine(i,j)
        else
           Input%Ref_063um_Clear = ch(1)%Ref_Toa_Clear(i,j)
        endif
        Input%Ref_063um_Std = ch(1)%Ref_Toa_Std_3x3(i,j)
        Input%Ref_063um_Min = ch(1)%Ref_Toa_Min_3x3(i,j)

        if (Ch(1)%Sub_Pixel_On_Flag) then
           Input%Ref_063um_Min_Sub = ch(1)%Ref_Toa_Min_Sub(i,j)
           Input%Ref_063um_Max_Sub = ch(1)%Ref_Toa_Max_Sub(i,j)
           Input%Ref_063um_Std_Sub = ch(1)%Ref_Toa_Std_Sub(i,j)
        else
           Input%Ref_063um_Min_Sub = MISSING_VALUE_REAL4
           Input%Ref_063um_Max_Sub = MISSING_VALUE_REAL4
           Input%Ref_063um_Std_Sub = MISSING_VALUE_REAL4
        endif

        Input%Log_Ref_063um_Std = Missing_Value_Real4
        if (Input%Ref_063um_Std > 0.0) Input%Log_Ref_063um_Std = log10(Input%Ref_063um_Std)
        if (Input%Log_Ref_063um_Std == 0.0) Input%Log_Ref_063um_Std = -2.0
      endif

      if (Input%Ref_063um_Min_Sub /= MISSING_VALUE_REAL4 .and. &
          Input%Ref_063um_Max_Sub /= MISSING_VALUE_REAL4)  then
          Input%Log_Drefl_065um_Max_Min_Sub = Missing_Value_Real4
          if ((Input%Ref_063um_Max_Sub - Input%Ref_063um_Min_Sub) >0.0) &
                 Input%Log_Drefl_065um_Max_Min_Sub = log10(Input%Ref_063um_Max_Sub - &
                 Input%Ref_063um_Min_Sub)
          if ((Input%Ref_063um_Max_Sub - Input%Ref_063um_Min_Sub) == 0.0) &
                 Input%Log_Drefl_065um_Max_Min_Sub = -2.0
      else
          Input%Log_Drefl_065um_Max_Min_Sub = MISSING_VALUE_REAL4
      endif

      if (Input%Chan_On_086um == sym%YES)  then 
        Input%Ref_086um = ch(2)%Ref_Toa(i,j)
      endif
      if (Input%Chan_On_138um == sym%YES)  then 
        Input%Ref_138um = ch(26)%Ref_Toa(i,j)
      endif
      if (Input%Chan_On_160um == sym%YES)  then 
        Input%Ref_160um = ch(6)%Ref_Toa(i,j)
        Input%Ref_160um_Clear = ch(6)%Ref_Toa_Clear(i,j)
      endif
      if (Input%Chan_On_213um == sym%YES)  then 
        Input%Ref_213um = ch(7)%Ref_Toa(i,j)
      endif
      if (Input%Chan_On_375um == sym%YES)  then 
        Input%Ref_375um = ch(20)%Ref_Toa(i,j)
        Input%Ref_375um_Clear = ch(20)%Ref_Toa_Clear(i,j)
        Input%Bt_375um = ch(20)%Bt_Toa(i,j)
        Input%Bt_375um_Clear = ch(20)%Bt_Toa_Clear(i,j)
        Input%Bt_375um_Std = ch(20)%Bt_Toa_Std_3x3(i,j)
        Input%Emiss_375um =  Ems_Ch20_Median_3x3(i,j)
        Input%Emiss_375um_Clear = Ch(20)%Emiss_Rel_11um_Clear(i,j)
      endif

      if (Input%Chan_On_67um == sym%YES)  then 
        Input%Bt_67um = ch(27)%Bt_Toa (i,j)
        if (Input%Chan_On_11um == sym%YES)  then 
           Input%Bt_11um_Bt_67um_Covar = Covar_Ch27_Ch31_5x5(i,j)
        endif
        if (Input%Chan_On_10um == sym%YES)  then 
           Input%Bt_10um_Bt_67um_Covar = Covar_Ch27_Ch38_5x5(i,j)
        endif
      endif

      if (Input%Chan_On_73um == sym%YES)  then 
        Input%Bt_73um = ch(28)%Bt_Toa (i,j)
      endif

      if (Input%Chan_On_133um == sym%YES)  then 
        Input%Bt_133um = ch(33)%Bt_Toa (i,j)
      endif

      if (AVHRR_IFF_Flag == 1)  then
        Input%Use_Sounder_11um = sym%YES
        Input%Bt_11um_Sounder = Bt_11um_Sounder(i,j)
        Input%Bt_11um_Bt_67um_Covar = Missing_Value_Real4
      endif
      if (Input%Chan_On_85um == sym%YES)  then 
        Input%Bt_85um = ch(29)%Bt_Toa(i,j)
      endif

      if (Input%Chan_On_10um == sym%YES)  then
        Input%Bt_10um = ch(38)%Bt_Toa(i,j)
        Input%Bt_10um_Std = ch(38)%Bt_Toa_Std_3x3(i,j)
        Input%Bt_10um_Max = ch(38)%Bt_Toa_Max_3x3(i,j)
        Input%Bt_10um_Clear = ch(38)%Bt_Toa_Clear(i,j)
        Input%Emiss_10um_Tropo = ch(38)%Emiss_Tropo(i,j)
        Input%Log_Bt_10um_Std = Missing_Value_Real4
        if (Input%Bt_10um_Std > 0.0) Input%Log_Bt_10um_Std = log10(Input%Bt_10um_Std)
        if (Input%Bt_10um_Std == 0.0) Input%Log_Bt_10um_Std = - 2.0
      else
        Input%Bt_10um = Missing_Value_Real4
        Input%Bt_10um_Std = Missing_Value_Real4
        Input%Bt_10um_Max = Missing_Value_Real4
        Input%Bt_10um_Clear = Missing_Value_Real4
        Input%Emiss_10um_Tropo = Missing_Value_Real4
        Input%Log_Bt_10um_Std = Missing_Value_Real4
      endif

      if (Input%Chan_On_11um == sym%YES)  then 
        Input%Bt_11um = ch(31)%Bt_Toa(i,j)
        Input%Bt_11um_Std = ch(31)%Bt_Toa_Std_3x3(i,j)
        Input%Bt_11um_Max = ch(31)%Bt_Toa_Max_3x3(i,j)
        Input%Bt_11um_Clear = ch(31)%Bt_Toa_Clear(i,j)
        Input%Emiss_11um_Tropo = ch(31)%Emiss_Tropo(i,j)
        Input%Log_Bt_11um_Std = Missing_Value_Real4
        if (Input%Bt_11um_Std > 0.0) Input%Log_Bt_11um_Std = log10(Input%Bt_11um_Std)
        if (Input%Bt_11um_Std == 0.0) Input%Log_Bt_11um_Std = - 2.0
        if (Ch(31)%Sub_Pixel_On_Flag) then
          Input%Bt_11um_Min_Sub = ch(31)%Bt_Toa_Min_Sub(i,j)
          Input%Bt_11um_Max_Sub = ch(31)%Bt_Toa_Max_Sub(i,j)
        else
          Input%Bt_11um_Min_Sub = Missing_Value_Real4
          Input%Bt_11um_Max_Sub = Missing_Value_Real4
        endif
      endif

      if (Input%Chan_On_12um == sym%YES)  then 
        Input%Bt_12um = ch(32)%Bt_Toa(i,j)
        Input%Bt_12um_Clear = ch(32)%Bt_Toa_Clear(i,j)
      endif
      if (Input%Chan_On_DNB == sym%YES)  then 
        Input%Lunscatzen = Geo%Scatangle_Lunar(i,j)
        Input%Lunar_Oceanic_Glint_Mask = Sfc%Glint_Mask_Lunar(i,j)
        Input%Rad_Lunar = ch(44)%Rad_Toa(i,j)
        Input%Ref_Lunar = ch(44)%Ref_Lunar_Toa(i,j)
        Input%Ref_Lunar_Min = Ref_ChDNB_Lunar_Min_3x3(i,j)
        Input%Ref_Lunar_Std = Ref_ChDNB_Lunar_Std_3x3(i,j)
        Input%Ref_Lunar_Clear = ch(44)%Ref_Lunar_Toa_Clear(i,j)
        Input%Lunzen = Geo%Lunzen(i,j)
        Input%Lunglintzen = Geo%Glintzen_Lunar(i,j)
        Input%Cld_Opd_DNB = ch(44)%Opd(i,j)
        Input%Log_Cld_Opd_DNB = Missing_Value_Real4
        if (Input%Cld_Opd_DNB > 0.0) Input%Log_Cld_Opd_DNB = log10(Input%Cld_Opd_DNB)
        if (Input%Cld_Opd_DNB == 0.0) Input%Log_Cld_Opd_DNB = - 2.0
      else
        Input%Lunscatzen = Missing_Value_Real4
        Input%Lunar_Oceanic_Glint_Mask = Missing_Value_Int1
        Input%Rad_Lunar = Missing_Value_Real4
        Input%Ref_Lunar = Missing_Value_Real4
        Input%Ref_Lunar_Min = Missing_Value_Real4
        Input%Ref_Lunar_Std = Missing_Value_Real4
        Input%Ref_Lunar_Clear = Missing_Value_Real4
        Input%Lunzen = Missing_Value_Real4
        Input%Lunglintzen = Missing_Value_Real4
        Input%Cld_Opd_DNB = Missing_Value_Real4
        Input%Log_Cld_Opd_DNB = Missing_Value_Real4
      endif

      ! - if use modawg mask
      Input%Use_Aux_Mask = Use_Aux_Flag

   end subroutine SET_INPUT

   subroutine SET_OUTPUT(i,j)
      integer, intent (in) :: i, j
      !CLDMASK%Prior_Cld_Probability(i,j) = Output%Prior ! Already saved
      CLDMASK%Cld_Mask_Qf(i,j) = Output%Cld_Mask_QF
      CLDMASK%Cld_Test_Vector_Packed(:,i,j) = Output%Cld_Flags_Packed
      CLDMASK%Cld_Mask(i,j) = Output%Cld_Mask_Bayes
      CLDMASK%Cld_Mask_Binary(i,j) = Output%Cld_Mask_Binary
      CLDMASK%Posterior_Cld_Probability(i,j) = Output%Posterior_Cld_Probability
      CLDMASK%Posterior_Ice_Probability(i,j) = Output%Posterior_Ice_Probability
      CLDMASK%Posterior_Water_Probability(i,j) = Output%Posterior_Water_Probability
      CLDMASK%TUT(i,j) = Output%TUT
      CLDMASK%RUT(i,j) = Output%RUT
      CLDMASK%Dust_Mask(i,j) = Output%Dust_Mask    
      CLDMASK%Smoke_Mask(i,j) = Output%Smoke_Mask
      CLDMASK%Fire_Mask(i,j) = Output%Fire_Mask
      CLDMASK%Thin_Cirr_Mask(i,j) = Output%Thin_Cirr_Mask
      CLDMASK%Bayes_Mask_Sfc_Type(i,j) = Output%Sfc_Idx
   end subroutine SET_OUTPUT

   subroutine SET_DIAG(i,j)
      integer, intent (in) :: i, j
      Diag_Pix_Array_1(i,j) = Diag%Array_1
      Diag_Pix_Array_2(i,j) = Diag%Array_2
      Diag_Pix_Array_3(i,j) = Diag%Array_3
   end subroutine SET_DIAG

!==============================================================
! Water Edge Filter
!
!  in an nxn box, if there are clear, water and ice, 
!  rephase water as ice
!
! Input: Phase = cloud phase following clavrx numbering
!        N = size of box  (N = 1 = 3x3 array)
!
! Output: Phase
!==============================================================
subroutine WATER_EDGE_FILTER(Phase,N)
  integer(kind=int1), intent(inout), dimension(:,:):: Phase
  integer:: N
  integer, dimension(size(Phase,1),size(Phase,2)):: Phase_Temp
  logical, dimension(2*N+1,2*N+1):: Mask
  integer, dimension(2*N+1,2*N+1):: Phase_Sub
  integer:: Nx, Ny
  integer:: i,i1,i2,j,j1,j2
  integer:: N_Ice, N_Clear
  integer:: N_Ice_Thresh, N_Clear_Thresh

  !--- set PHASE numbering
  integer, parameter:: CLEAR_PHASE = 0
  integer, parameter:: WATER_PHASE = 1
  integer, parameter:: SUPERCOOLED_PHASE = 2
  integer, parameter:: MIXED_PHASE = 3
  integer, parameter:: ICE_PHASE = 4
  integer, parameter:: UNKNOWN_PHASE = 5

  integer:: N_Water_Before, N_Water_After

  !--- determine segment size
  Nx = size(Phase,1)
  Ny = size(Phase,2)

  !--- make a copy
  Phase_Temp = Phase

  !--- set thresholds of numbers of clear and ice pixels for filter
  N_Clear_Thresh = 1
  N_Ice_Thresh = 1

  !--- loop over segment
  do i =  1, Nx-N, 1

   i1 = max(1, i - N)
   i2 = min(Nx,i + N)

   do j = 1, Ny-N, 1

     j1 = max(1, j - N)
     j2 = min(Ny,j + N)

     !--- if in a corner, skip
     if (i2 - i1 < 1 .and. j2 - j1 < 1) cycle

     !--- initialize subarrays
     Mask = .false.
     Phase_Sub = UNKNOWN_PHASE
     Phase_Sub(1:i2-i1+1,1:j2-j1+1) = Phase(i1:i2,j1:j2)

     !--- determine number of clear in sub-array
     Mask = .false.
     where(Phase_Sub == CLEAR_PHASE)
         Mask = .true.
     endwhere
     N_Clear = count(Mask)
     if (N_Clear < N_Clear_Thresh) cycle    !skip if no clear

     !--- determine number of ice in subarray
     Mask = .false.
     where(Phase_Sub == ICE_PHASE)
         Mask = .true.
     endwhere
     N_Ice = count(Mask)
     if (N_Ice < N_Ice_Thresh) cycle    !skip if no ice

     Mask = .false.
     where(Phase_Sub == WATER_PHASE .or. Phase_Sub == SUPERCOOLED_PHASE)
         Mask = .true.
     endwhere
     N_Water_Before = count(Mask)

     !--- apply filter
     if (N_Clear >= N_Clear_Thresh .and. N_Ice >= N_Ice_Thresh) then
       where(Phase_Sub == WATER_PHASE .or. Phase_Sub == MIXED_PHASE .or. &
             Phase_Sub == SUPERCOOLED_PHASE)
        Phase_Sub = ICE_PHASE
       endwhere
     endif

     Mask = .false.
     where(Phase_Sub == WATER_PHASE .or. Phase_Sub == SUPERCOOLED_PHASE)
         Mask = .true.
     endwhere
     N_Water_After = count(Mask)

     !--- copy sub-array back into full array
     Phase_Temp(i1:i2,j1:j2) = Phase_Sub(1:i2-i1+1,1:j2-j1+1)

    enddo
  enddo

  !--- copy back
  Phase = Phase_Temp

end subroutine WATER_EDGE_FILTER
!==============================================================
! Median filter
!
! mask = 0 means use median, mask = 1 mean ignore
!==============================================================
subroutine MEDIAN_FILTER_MASK(z,mask,z_median)

! The purpose of this function is to find 
! median (emed), minimum (emin) and maximum (emax)
! for the array elem with nelem elements. 

  real, dimension(:,:), intent(in):: z
  real, intent(out):: z_median
  integer(kind=int1), dimension(:,:), intent(in):: mask
  integer:: i,j,k,nx,ny,nelem
  real, dimension(:), allocatable::x
  real:: u

  z_median = missing_value_real4

  nx = size(z,1)
  ny = size(z,2)

  nelem = nx * ny

  allocate(x(nelem))
  x = 0.0
  k = 0
  do i = 1, nx
    do j = 1, ny
      if (mask(i,j) == 0 .and. z(i,j) /= missing_value_int1) then
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

  if (allocated(x)) deallocate(x)

end subroutine MEDIAN_FILTER_MASK

!============================================================================
! Pull channel initliazation out of pixel loop - WCS3
!============================================================================
subroutine SET_ECM_CHAN_ON ()

  Input%Chan_On_041um = Sensor%Chan_On_Flag_Default(8)
  Input%Chan_On_063um = Sensor%Chan_On_Flag_Default(1)
  Input%Chan_On_086um = Sensor%Chan_On_Flag_Default(2)
  Input%Chan_On_138um = Sensor%Chan_On_Flag_Default(26)
  Input%Chan_On_160um = Sensor%Chan_On_Flag_Default(6)
  Input%Chan_On_213um = Sensor%Chan_On_Flag_Default(7)
  Input%Chan_On_375um = Sensor%Chan_On_Flag_Default(20)
  Input%Chan_On_62um = Sensor%Chan_On_Flag_Default(37)
  Input%Chan_On_67um = Sensor%Chan_On_Flag_Default(27)
  Input%Chan_On_73um = Sensor%Chan_On_Flag_Default(28)
  Input%Chan_On_85um = Sensor%Chan_On_Flag_Default(29)
  Input%Chan_On_97um = Sensor%Chan_On_Flag_Default(30)
  Input%Chan_On_10um = Sensor%Chan_On_Flag_Default(38)
  Input%Chan_On_11um = Sensor%Chan_On_Flag_Default(31)
  Input%Chan_On_12um = Sensor%Chan_On_Flag_Default(32)
  Input%Chan_On_133um = Sensor%Chan_On_Flag_Default(33)
  Input%Chan_On_I1_064um = Sensor%Chan_On_Flag_Default(39)
  Input%Chan_On_I4_374um = Sensor%Chan_On_Flag_Default(42)
  Input%Chan_On_I5_114um = Sensor%Chan_On_Flag_Default(43)
  Input%Chan_On_DNB = Sensor%Chan_On_Flag_Default(44)

  !For GOES-17, reset local chan on
  IF (Sensor%WMO_Id ==271) CALL LHP_CHN_CHECK()

end subroutine SET_ECM_CHAN_ON

SUBROUTINE LHP_CHN_CHECK ()
  REAL, DIMENSION(Nchan_Clavrx) :: LHP_THRESH
  REAL, DIMENSION(Nchan_Clavrx) :: LHP_THRESH_INIT
  REAL :: MFPT_104, MFPT_110
  INTEGER :: N_chan
  CHARACTER(3000) :: Algo_Thresh_File
  REAL, PARAMETER:: SOUNDER_LHP_THRESH = 999.0

  !initalize LHP THRESHOLD
  LHP_THRESH (:) = MISSING_VALUE_REAL4
  LHP_THRESH_INIT (:) = MISSING_VALUE_REAL4

  !Initialize 10.3 and 11um FPT local variables
  MFPT_104 = MISSING_VALUE_REAL4
  MFPT_110 = MISSING_VALUE_REAL4

  !Initialize Use_104um_Flag
  Use_104um_Flag = .FALSE.

  !--- Set local variables for easier reading.
  MFPT_104 = Ch(38)%Max_Focal_Plane_Temp
  MFPT_110 = Ch(31)%Max_Focal_Plane_Temp

  !--- Initialize LHP to default thresholds
  LHP_THRESH_INIT(20) = ABI_FPT_Thresh_038um
  LHP_THRESH_INIT(27) = ABI_FPT_Thresh_067um
  LHP_THRESH_INIT(28) = ABI_FPT_Thresh_073um
  LHP_THRESH_INIT(29) = ABI_FPT_Thresh_085um
  LHP_THRESH_INIT(38) = ABI_FPT_Thresh_104um
  LHP_THRESH_INIT(31) = ABI_FPT_Thresh_110um
  LHP_THRESH_INIT(32) = ABI_FPT_Thresh_120um

  Algo_Thresh_File =trim(Ancil_Data_Dir)// &
                     "static/algo_lhp_thresh/"//&
                     "goes17_thermal_limits_ecm.txt"

  !Read in LHP Thresholds
  if (Image%Segment_Number == 1 ) then
    print*, "Opening ", trim(Algo_Thresh_File)
  endif

  CALL READ_LHP_THRESH_FILE(TRIM(Algo_Thresh_File), LHP_THRESH)

  !Check if LHP_THRESH is missing
  do N_Chan = 1, Nchan_Clavrx
     IF (LHP_THRESH(N_Chan) == MISSING_VALUE_REAL4) &
        LHP_THRESH(N_Chan) = LHP_THRESH_INIT(N_Chan)

  end do

   !--- Nullify LHP_THRESH is sounder (aka fusion) data is used
   do N_Chan = 1, Nchan_Clavrx
      if (Sensor%Chan_On_Flag_Default(N_Chan) == sym%NO) cycle
      if (maxval(Ch(N_Chan)%Source) == 1) then
         LHP_THRESH(N_Chan) = SOUNDER_LHP_THRESH
      endif
   enddo

  !for now limit to ABI channels
  Input%Chan_On_375um = LHP_Local_Chan_On(20, LHP_THRESH(20))
  Input%Chan_On_67um = LHP_Local_Chan_On(27, LHP_THRESH(27))
  Input%Chan_On_73um = LHP_Local_Chan_On(28, LHP_THRESH(28))
  Input%Chan_On_85um = LHP_Local_Chan_On(29, LHP_THRESH(29))
  Input%Chan_On_10um = LHP_Local_Chan_On(38, LHP_THRESH(38))
  Input%Chan_On_11um = LHP_Local_Chan_On(31, LHP_THRESH(31))
  Input%Chan_On_12um = LHP_Local_Chan_On(32, LHP_THRESH(32))
  Input%Chan_On_133um = LHP_Local_Chan_On(33, LHP_THRESH(33))

  !10.4um check
  IF ((MFPT_104 < LHP_THRESH(38)) .AND. &
      (MFPT_110 > LHP_THRESH(31))) then
      Use_104um_Flag = .TRUE.
  ENDIF

END SUBROUTINE LHP_CHN_CHECK

   ! -----------------------------------------------------
   !    Overlap test using Visible and Split-Window
   !  
   !    05/01/2014 AW
   ! ---------------------------------------------------- 
   subroutine VIS_SPLIT_WINDOW_OVERLAP_TEST (bt_11 &
                        , bt_12 &
                        , ref_vis &
                        , sol_zen &
                        , sat_zen &
                        , is_overlap )

      real, intent(in) :: bt_11
      real, intent(in) :: bt_12
      real, intent(in) :: ref_vis
      real, intent(in) :: sol_zen
      real, intent(in) :: sat_zen

      logical, intent(out) :: is_overlap

      real :: btd_thresh
      real :: ref_vis_thresh

      is_overlap = .false.

      btd_thresh = 1.0
      ref_vis_thresh = 40.0

      if (( bt_11 - bt_12 ) > btd_thresh &
         .and.  ref_vis > ref_vis_thresh ) is_overlap = .true.

   end subroutine VIS_SPLIT_WINDOW_OVERLAP_TEST
   !---------------------------------------------------- 
   ! Beta tests for Overlap
   !  
   !    05/01/2014 AW
   !---------------------------------------------------- 
   subroutine BETA_11_12_OVERLAP_TEST (emiss_tropo_11,beta_11_12,is_overlap)

      real, intent(in):: emiss_tropo_11
      real, intent(in):: beta_11_12
      logical, intent(out) :: is_overlap

      real, parameter :: EMISS_11UM_OVERLAP_THRESH = 0.95
      real, parameter :: BETA_11UM_12UM_OVERLAP_THRESH = 0.95

      !------------------------------------------------------------------------
      !  define cirrus vs opaque by emiss_trop thresholds
      !------------------------------------------------------------------------
      is_overlap = .false.
      if ( emiss_tropo_11 < EMISS_11UM_OVERLAP_THRESH) then

         if ( beta_11_12 > 0.0 .and. beta_11_12 < BETA_11UM_12UM_OVERLAP_THRESH ) then
            is_overlap = .true.
         end if

      end if
   end subroutine BETA_11_12_OVERLAP_TEST
   !---------------------------------------------------- 
   subroutine BETA_11_133_OVERLAP_TEST (emiss_tropo_11,beta_11_133,is_overlap)

      real, intent(in):: emiss_tropo_11
      real, intent(in):: beta_11_133
      logical, intent(out) :: is_overlap

      real, parameter :: EMISS_11UM_OVERLAP_THRESH = 0.95
      real, parameter :: BETA_11UM_133UM_OVERLAP_THRESH = 0.70

      !------------------------------------------------------------------------
      !  define cirrus vs opaque by emiss_trop thresholds
      !------------------------------------------------------------------------
      is_overlap = .false.
      if ( emiss_tropo_11 < EMISS_11UM_OVERLAP_THRESH) then

         if ( beta_11_133 > 0.0 .and. beta_11_133 < BETA_11UM_133UM_OVERLAP_THRESH) then
            is_overlap = .true.
         end if

      end if
   end subroutine BETA_11_133_OVERLAP_TEST
   !-------------------------------------------------------------------------------
   ! compute the CLAVR-x cloud type from the phase from ECM2
   ! this is taken from the universal_cloud_type.f90 module 
   !--------------------------------------------------------------------------------
   subroutine COMPUTE_TYPE_FROM_PHASE(Cld_Mask,Cld_Phase,Cld_Phase_Uncer,Zopa,Topa,Etropo, &
                                      Chan_On_11um, Chan_On_12um, Chan_On_133um, &
                                      Beta_11_12, Beta_11_133, Cld_Type)

   integer(kind=int1), intent(in):: Cld_Mask, Cld_Phase
   integer(kind=int1), intent(in):: Chan_On_11um, Chan_On_12um, Chan_On_133um
   real, intent(in):: Zopa,Topa,Etropo,Beta_11_12,Beta_11_133,Cld_Phase_Uncer
   integer(kind=int1), intent(out):: Cld_Type
   logical:: Is_Overlap

        !type
        if (Cld_Mask == sym%CLEAR) Cld_Type = sym%CLEAR_TYPE
        if (Cld_Mask == sym%PROB_CLEAR) Cld_Type = sym%PROB_CLEAR_TYPE
        if (Cld_Phase == sym%SUPERCOOLED_PHASE) Cld_Type = sym%SUPERCOOLED_TYPE
        if (Cld_Phase == sym%WATER_PHASE) then
           Cld_Type = sym%WATER_TYPE
           if ((Zopa .ner. Missing_Value_Real4) .and. (Zopa .ltr. 1000.0)) Cld_Type = sym%FOG_TYPE
           if (Topa .ltr. 273.15) Cld_Type = sym%SUPERCOOLED_TYPE
        endif 
        if (Cld_Phase == sym%ICE_PHASE) then
           Cld_Type = sym%CIRRUS_TYPE
           if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
             if (Topa .ltr. 233.0) Cld_Type = sym%OPAQUE_ICE_TYPE
             if (Etropo .gtr. 0.80) Cld_Type = sym%OPAQUE_ICE_TYPE
             if (Etropo .gtr. 0.95) Cld_Type = sym%OVERSHOOTING_TYPE
           endif
        endif

        !spectral overlap
!       if (Cld_Type  == sym%CIRRUS_TYPE) then
!             if (Chan_On_11um == sym%YES) then
!                if (Chan_On_12um == sym%YES) then
!                   call BETA_11_12_OVERLAP_TEST (Etropo,Beta_11_12,Is_Overlap)
!                endif
!                if (Chan_On_133um == sym%YES) then
!                   call BETA_11_133_OVERLAP_TEST (Etropo,Beta_11_133,Is_Overlap)
!                endif
!             endif
!          if (Is_Overlap) Cld_Type = sym%OVERLAP_TYPE
!        endif

        !overlap from uncertainty
!       if ((Cld_Phase == sym%ICE_PHASE) .and. &
!           (Cld_Phase_Uncer .ger. CLD_PHASE_UNCER_MULTI_THRESH)) then
!           Cld_Type = sym%OVERLAP_TYPE
!       endif

        !overlap from uncertainty
        if ((Cld_Type == sym%CIRRUS_TYPE) .and. &
            (Cld_Phase_Uncer .ger. CLD_PHASE_UNCER_MULTI_THRESH)) then
            Cld_Type = sym%OVERLAP_TYPE
        endif
        if ((Cld_Phase == sym%WATER_PHASE) .and. &
            (Cld_Phase_Uncer .ger. CLD_PHASE_UNCER_MULTI_THRESH)) then
            Cld_Type = sym%MIXED_TYPE
        endif

   end subroutine COMPUTE_TYPE_FROM_PHASE

!-------------------------------------------------------------------------------

end module NBM_CLOUD_MASK_CLAVRX_BRIDGE
