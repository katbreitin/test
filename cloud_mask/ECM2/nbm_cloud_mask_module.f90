!$Id:$ nbm_cloud_mask_MODULE.f90 
!----------------------------------------------------------------------
! MODULE name: NBM_CLOUD_MASK
! 
! Routines for the determination of the naive Bayesian cloud mask
! Version 2.0
!
! Authors: Andrew Heidinger, NOAA/NESDIS
!          Andi Walther, CIMSS
!          Denis Botambekov, CIMSS
!          William Straka, CIMSS
!
! DEPENDENCIES: Services_Module
!
! SIDE EFFECTS: None
!
! Cld_Flags Format
!
! Flag Number  bit-depth    bits      byte name
! 1            1            1         1    cloud mask attempted       
! 2            1            2         1    day_063     
! 3            1            3         1    day_063_spatial     
! 4            1            4         1    day_375
! 5            1            5         1    night_375
! 6            1            6         1    solar contamination 
! 7            1            7         1    coast
! 8            1            8         1    mountain
!----
! 9            1            9         2    forward scattering
! 10           1            10        2    snow (no snow = 0, snow/sea ice = 1
! 11           1            11        2    cold scene btd
! 12           1            12        2    glint
! 13           1            13        2    smoke detected
! 14           1            14        2    dust  detected
! 15           1            15        2    shadow detected
! 16           1            16        2    fire detected
!---
! 17           3            17-19     3    nbvcm surface TYPE
! 18           1            20        3    thin cirrus detected
! (separate test)
!<-------------------- START OF CLOUD TESTS -------------------------->
! other ECM test bits are dynamic
!
! TODO ADD DESCRIPTION 
!
!----------------------------------------------------------------------

MODULE NBM_CLOUD_MASK_MODULE

 USE NB_CLOUD_MASK_SERVICES
 USE NBM_CLOUD_MASK_LUT_MODULE
 USE NBM_CLOUD_MASK_GET_PROB_MASK_PHASE

 IMPLICIT NONE

 PRIVATE:: COMPUTE_BAYES_SFC_TYPE
 PRIVATE:: TUT_ROUTINE
 PRIVATE:: RUT_ROUTINE
 PRIVATE:: PACK_BITS_INTO_BYTES
 PRIVATE:: SET_NONCLOUD_FLAGS

 PUBLIC:: NBM_CLOUD_MASK_ALGORITHM
 PUBLIC:: SET_NBM_CLOUD_MASK_VERSION

 !--- set thresholds and algorithm specific constants
 INCLUDE 'nbm_cloud_mask.inc'

 !--- string to control on-screen prompts
 CHARACTER(*), PARAMETER, PRIVATE :: EXE_PROMPT_CM = "Naive Bayesian Cloud Mask Version 2.0 >> "


 CONTAINS

!====================================================================
!  pass threshold version to bridge
!  Cloud_Mask_Thresholds_Version is a MODULE-wide variable
!  that is PRIVATE
!====================================================================
 SUBROUTINE SET_NBM_CLOUD_MASK_VERSION(Cloud_Mask_Version)
   CHARACTER(len=*), INTENT(OUT):: Cloud_Mask_Version
   Cloud_Mask_Version = "$Id: nbm_cloud_mask_MODULE.f90 3030 2019-09-31 23:02:31Z heidinger $"
 END SUBROUTINE SET_NBM_CLOUD_MASK_VERSION


!====================================================================
! SUBROUTINE Name: NBM_CLOUD_MASK_ALGORITHM
!
! Function:
!   CalcuLates the bayesian cloud mask. The bayesian cloud mask is
!   determined by utilizing the following surface TYPEs:
!
! Bayesian Surface Types
! 1 - Deep_Water
! 2 - Shallow_Water
! 3 - Unfrozen_Land
! 4 - Frozen_Land
! 5 - Arctic
! 6 - Antarctic
! 7 - Desert
!
!====================================================================
 SUBROUTINE NBM_CLOUD_MASK_ALGORITHM( &
            Nclass, &    
            Symbol, &    
            Input,  &
            Output,  &
            Use_Prior_Table, &
            ABI_Use_104um_Flag, &
            Diag)

    INTEGER, INTENT(IN) :: Nclass
    TYPE(symbol_naive_bayesian), INTENT(IN) :: Symbol
    TYPE(mask_input), INTENT(IN) :: Input
    TYPE(mask_output), INTENT(OUT) :: Output
    LOGICAL, INTENT(IN):: Use_Prior_Table
    LOGICAL, INTENT(IN):: ABI_Use_104um_Flag
    TYPE(diag_output), INTENT(OUT), Optional :: Diag

    INTEGER, save:: Diag_Warning_Flag = 0

    REAL, DIMENSION(:), ALLOCATABLE :: Clear_Cond_Ratio
    REAL, DIMENSION(:), ALLOCATABLE :: Water_Cond_Ratio
    REAL, DIMENSION(:), ALLOCATABLE :: Ice_Cond_Ratio
    REAL, DIMENSION(:), ALLOCATABLE :: Obs_Prob
    REAL, DIMENSION(:), ALLOCATABLE :: Posterior_Cld_Probability_By_Class
    REAL :: Post_Prob_Clear
    REAL :: Prior_Yes_Temp
    REAL :: Post_Prob_Clear_By_Class
    REAL, PARAMETER :: Prob_Min_Thresh = 0.001
    REAL, PARAMETER :: R_Max = 1000.0
    REAL :: Value_Dim
    REAL :: X_Dim, Y_Dim, Z_Dim
    REAL :: R_Clear
    REAL :: R_Water
    REAL :: R_Ice
    REAL :: Prior_Prob_Clear, Prior_Prob_Water, Prior_Prob_Ice
    REAL :: Z1, Z2, Z3, Z4
    REAL :: Post_Sum
    CHARACTER(len=DEFAULT_NAME_LENGTH) :: Dim_Name    !changed to 50 from 30: AKH
    INTEGER, PARAMETER :: Nchan = 45
    INTEGER :: Class_Idx
    INTEGER :: Number_Valid_Classes
    INTEGER :: Chan_Idx
    INTEGER :: Elem_Idx, Line_Idx
    INTEGER :: i
    INTEGER :: Is_Land
    INTEGER, DIMENSION(:), ALLOCATABLE :: Chan_On
    INTEGER, DIMENSION(:), ALLOCATABLE :: Chan_Wvl
    INTEGER, DIMENSION(:), ALLOCATABLE :: Class_Use_Flag

    INTEGER:: Oceanic_Glint_Flag
    INTEGER:: Lunar_Oceanic_Glint_Flag
    INTEGER:: Coastal_Flag
    INTEGER:: Mountain_Flag
    INTEGER:: Day_063_Flag
    INTEGER:: Day_063_Spatial_Flag
    INTEGER:: Night_Lunar_Flag
    INTEGER:: Lunar_Spatial_Flag
    INTEGER:: Day_375_Flag
    INTEGER:: Night_375_Flag
    INTEGER:: Forward_Scattering_Flag
    INTEGER:: Solar_Contam_Flag
    INTEGER:: Lunar_Forward_Scattering_Flag
    INTEGER:: Snow_Flag
    INTEGER:: Cold_Scene_Flag
    INTEGER:: Dry_Scene_Flag
    INTEGER:: City_Flag

    REAL (KIND=REAL4):: Airmass
    INTEGER, PARAMETER:: Spare_Value = 0
    INTEGER (KIND=INT1):: Use_104_Flag

    !-- local pointers that point to global variables
    INTEGER(KIND=INT1), DIMENSION(NUMBER_OF_CLOUD_FLAGS+NUMBER_OF_NONCLOUD_FLAGS):: Cld_Flags
    INTEGER(KIND=INT1), DIMENSION(NUMBER_OF_CLOUD_FLAGS+NUMBER_OF_NONCLOUD_FLAGS):: Cld_Flag_Bit_Depth   


    !------------------------------------------------------------------------------------------
    !---  begin executable code
    !------------------------------------------------------------------------------------------

    !--- initialize diagnostic output
    IF (present(Diag) .and. Diag_Warning_Flag == 0) THEN
       PRINT *, "CLAVR-x / NB Cloud Mask ===>  Diagnostic Output Turned On"
       Diag_Warning_Flag = 1
    ENDIF
    IF (present(Diag)) Diag%Array_1 = MISSING_VALUE_REAL4
    IF (present(Diag)) Diag%Array_2 = MISSING_VALUE_REAL4
    IF (present(Diag)) Diag%Array_3 = MISSING_VALUE_REAL4

    !--- initialize output
    Output%Posterior_Cld_Probability = MISSING_VALUE_REAL4
    Output%Posterior_Cld_Probability_IR = MISSING_VALUE_REAL4
    Output%Cld_Mask_Bayes = MISSING_VALUE_INT1
    Output%Cld_Mask_Bayes_IR = MISSING_VALUE_INT1
    Output%Cld_Mask_Binary = MISSING_VALUE_INT1
    Output%Cld_Mask_Binary_IR = MISSING_VALUE_INT1
    Output%Cld_Flags_Packed = 0
    Output%Sfc_Idx = -1

    Cld_Flags(1) = symbol%YES          ;    Cld_Flag_Bit_Depth(1) = 1

    !--- check for a bad pixel
    IF (Input%Invalid_Data_Mask == symbol%YES) THEN
        Cld_Flags(1) = symbol%NO
    ELSE

        !---  COMPUTE SURFACE TYPE
        Output%Sfc_Idx = COMPUTE_BAYES_SFC_TYPE(Input%Land_Class, &
                                 Input%Coastal_Mask, &
                                 Input%Snow_Class, &
                                 Input%Sfc_Type, &
                                 Input%Lat, &
                                 Input%Lon, &
                                 Input%Sst_Anal_Uni, &
                                 Input%Emiss_Sfc_375um, &
                                 symbol)

        !------------------------------------------------------------------------------
        !  set channel on/off and channel mapping
        !  Chan_On = channel on / off binary flags
        !  Chan_Wvl = channel nominal wavelength in nm
        !  This is set up using CLAVR-x channel numbers
        !------------------------------------------------------------------------------
        ALLOCATE(Chan_On(Nchan))
        ALLOCATE(Chan_Wvl(Nchan))
        Chan_On = 0
        Chan_Wvl = 0

        Chan_On(1) = Input%Chan_On_063um
        Chan_On(2) = Input%Chan_On_086um
        Chan_On(26) = Input%Chan_On_138um
        Chan_On(6) = Input%Chan_On_160um
        Chan_On(20) = Input%Chan_On_375um
        Chan_On(27) = Input%Chan_On_67um
        Chan_On(28) = Input%Chan_On_73um
        Chan_On(29) = Input%Chan_On_85um
        Chan_On(30) = Input%Chan_On_97um
        Chan_On(31) = Input%Chan_On_11um
        Chan_On(32) = Input%Chan_On_12um
        Chan_On(33) = Input%Chan_On_133um
        Chan_On(37) = Input%Chan_On_62um
        Chan_On(38) = Input%Chan_On_10um
        Chan_On(44) = Input%Chan_On_DNB


        Chan_Wvl(1) = 650
        Chan_Wvl(2) = 860
        Chan_Wvl(26) = 1380
        Chan_Wvl(6) = 1600
        Chan_Wvl(20) = 3750
        Chan_Wvl(27) = 6700
        Chan_Wvl(28) = 7300
        Chan_Wvl(29) = 8500
        Chan_Wvl(30) = 9700
        Chan_Wvl(31) = 11000
        Chan_Wvl(32) = 12000
        Chan_Wvl(33) = 13300
        Chan_Wvl(37) = 6200
        Chan_Wvl(38) = 10400
        Chan_Wvl(44) = 700   !DNB


        !----- compute prior
        IF (Input%Prior /= MISSING_VALUE_REAL4 .and. Use_Prior_Table .eqv. .true.) THEN
           Prior_Yes_Temp = Input%Prior
        ELSE
           Prior_Yes_Temp = MISSING_VALUE_REAL4 ! Set further down - WCS3
        ENDIF

        ! - ALLOCATE needed variables
        IF (ALLOCATED(Class_Use_Flag)) DEALLOCATE (Class_Use_Flag)
        IF (ALLOCATED(Clear_Cond_Ratio)) DEALLOCATE (Clear_Cond_Ratio)
        IF (ALLOCATED(Water_Cond_Ratio)) DEALLOCATE (Water_Cond_Ratio)
        IF (ALLOCATED(Ice_Cond_Ratio)) DEALLOCATE (Ice_Cond_Ratio)
        IF (ALLOCATED(Obs_Prob)) DEALLOCATE (Obs_Prob)
        IF (ALLOCATED(Posterior_Cld_Probability_By_Class)) DEALLOCATE (Posterior_Cld_Probability_By_Class)

        ALLOCATE (Class_Use_Flag(Nclass))
        ALLOCATE (Clear_Cond_Ratio(Nclass))
        ALLOCATE (Water_Cond_Ratio(Nclass))
        ALLOCATE (Ice_Cond_Ratio(Nclass))
        ALLOCATE (Obs_Prob(Nclass))
        ALLOCATE (Posterior_Cld_Probability_By_Class(Nclass))

        ! - set to missing
        Output%Posterior_Cld_Probability = MISSING_VALUE_REAL4
        Posterior_Cld_Probability_By_Class = MISSING_VALUE_INT1
        Post_Prob_Clear = MISSING_VALUE_REAL4
        Output%Posterior_Water_Probability = MISSING_VALUE_REAL4
        Output%Posterior_Ice_Probability = MISSING_VALUE_REAL4
        Output%Cld_Mask_Bayes = MISSING_VALUE_INT1
        Output%Cld_Mask_Binary = MISSING_VALUE_INT1
        Clear_Cond_Ratio = MISSING_VALUE_REAL4
        Water_Cond_Ratio = MISSING_VALUE_REAL4
        Ice_Cond_Ratio = MISSING_VALUE_REAL4
        Obs_Prob = MISSING_VALUE_REAL4

        R_Clear = 1.0
        R_Water = 1.0
        R_Ice = 1.0
        Class_Use_Flag = 1
    
        !-------------------------------------------------------------------
        ! turn off tables based on chan_on
        !-------------------------------------------------------------------
        DO Class_Idx = 1, Nclass
           DO Chan_Idx = 1, Nchan

            IF (Chan_Idx > Lut(Class_Idx)%Nchan_Used) CYCLE
            IF (Chan_On(Chan_Idx) == 0) THEN
                  IF (Chan_Wvl(Chan_Idx) == Lut(Class_Idx)%Wvl(Chan_Idx) .and. &
                      Chan_Wvl(Chan_Idx) /= 0) THEN
                         Class_Use_Flag(Class_Idx) = 0
                  ENDIF
            ENDIF
           ENDDO
        ENDDO



        Number_Valid_Classes = sum(Class_Use_Flag)
        IF (Number_Valid_Classes == 0) THEN
           PRINT *, 'no active tests, Number_Valid_Classes =',Number_Valid_Classes
           stop
        ENDIF


   ! - loop over classifiers
   DO Class_Idx = 1, Nclass

      ! - set up input
      X_Dim = MISSING_VALUE_REAL4
      Y_Dim = MISSING_VALUE_REAL4
      Z_Dim = MISSING_VALUE_REAL4      
      
      
      ! - loop over class rank
      DO i = 1, Lut(Class_Idx)%Rank

        ! - SELECT classifier name
        SELECT CASE(i)
             CASE(1)
                 Dim_Name = Lut(Class_Idx)%Class_Xname
             CASE(2)
                 Dim_Name = Lut(Class_Idx)%Class_Yname
             CASE(3)
                 Dim_Name = Lut(Class_Idx)%Class_Zname
        END SELECT

         ! - SELECT value based on the DIMENSION name

         Value_Dim = MISSING_VALUE_REAL4   !AKH Added this initialization

         SELECT CASE(trim(Dim_Name))
             CASE('etropo10')
                 Value_Dim = Input%Emiss_10um_Tropo
             CASE('etropo11')
                 Value_Dim = Input%Emiss_11um_Tropo
             CASE('topa')
                 Value_Dim = Input%Topa
             CASE('zopa')
                 Value_Dim = Input%Zopa
             CASE('dtsfcopa')
                 Value_Dim = Input%Sfc_Temp - Input%Topa
             CASE('logzopa')
                 if (Input%Zopa > 0.0) Value_Dim = alog10(Input%Zopa)
                 if (Input%Zopa < 0.0 .and. Input%Zopa /= MISSING_VALUE_REAL4) Value_Dim = 1.0
             CASE('bt10')
                 Value_Dim = Input%Bt_10um
                 if (Input%Chan_On_10um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('bt10std')
                 Value_Dim = Input%Bt_10um_Std
                 if (Input%Chan_On_10um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('logbt10std')
                 Value_Dim = Input%Log_Bt_10um_Std
                 if (Input%Chan_On_10um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('bt11')
                 Value_Dim = Input%Bt_11um
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('bt11minsub')
                 Value_Dim = Input%Bt_11um_Min_Sub
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('bt11maxsub')
                 Value_Dim = Input%Bt_11um_Max_Sub
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('dbt11max')
                 Value_Dim = Input%Bt_11um_Max - Input%Bt_11um
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('logdbt11max3x3')
                 if (Input%Chan_On_11um == 0) CYCLE
                 Value_Dim = Input%Bt_11um_Max - Input%Bt_11um
                 if (Value_Dim > 0) Value_Dim = alog10(Value_Dim)
                 if (Value_Dim <= 0) Value_Dim = -2.0
             CASE('logdbt11maxminsub')
                 if (Input%Chan_On_11um == 0) CYCLE
                 if (Input%Bt_11um_Max_Sub == MISSING_VALUE_REAL4) CYCLE
                 if (Input%Bt_11um_Min_Sub == MISSING_VALUE_REAL4) CYCLE
                 Value_Dim = Input%Bt_11um_Max_Sub - Input%Bt_11um_Min_Sub
                 if (Value_Dim > 0) Value_Dim = alog10(Value_Dim)
                 if (Value_Dim <= 0) Value_Dim = -2.0
             CASE('dbt10max3x3')
                 Value_Dim = Input%Bt_10um_Max - Input%Bt_10um
                 if (Input%Chan_On_10um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('dbt11max3x3')
                 Value_Dim = Input%Bt_11um_Max - Input%Bt_11um
             CASE('dbt11maxsub')
                 Value_Dim = Input%Bt_11um_Max_Sub - Input%Bt_11um
             CASE('bt11std')
                 Value_Dim = Input%Bt_11um_Std 
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btdclr10')
                 Value_Dim = Input%Bt_10um_Clear - Input%Bt_10um
                 if (Input%Chan_On_10um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btdclr11')
                 Value_Dim = Input%Bt_11um_Clear - Input%Bt_11um
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btdclr12')
                 Value_Dim = Input%Bt_12um_Clear - Input%Bt_12um
                 if (Input%Chan_On_12um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('logbt11std')
                 Value_Dim = Input%Log_Bt_11um_Std
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btd3810')
                 Value_Dim = Input%Bt_375um - Input%Bt_10um
                 if (Input%Chan_On_375um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_10um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btd3811')
                 Value_Dim = Input%Bt_375um - Input%Bt_11um
                 if (Input%Chan_On_375um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('dbtd3811clr')
                 Value_Dim = (Input%Bt_375um - Input%Bt_11um) - &
                             (Input%Bt_375um_Clear - Input%Bt_11um_Clear) 
                 if (Input%Chan_On_375um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('emiss3810')
                 Value_Dim = Input%Emiss_375um - Input%Emiss_10um_Tropo
                 if (Input%Chan_On_375um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_10um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('emiss3811')
                 Value_Dim = Input%Emiss_375um - Input%Emiss_11um_Tropo
                 if (Input%Chan_On_375um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btd8573')
                 Value_Dim = Input%Bt_85um - Input%Bt_73um
                 if (Input%Chan_On_85um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_73um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btd8511')
                 Value_Dim = Input%Bt_85um - Input%Bt_11um
                 if (Input%Chan_On_85um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btd1167')
                 Value_Dim = Input%Bt_11um - Input%Bt_67um
                 if (Input%Chan_On_67um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btd1173')
                 Value_Dim = Input%Bt_11um - Input%Bt_73um
                 if (Input%Chan_On_73um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btd11133')
                 Value_Dim = Input%Bt_11um - Input%Bt_133um
                 if (Input%Chan_On_133um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btd13373')
                 Value_Dim = Input%Bt_133um - Input%Bt_73um
                 if (Input%Chan_On_133um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_73um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btd1110')
                 Value_Dim = Input%Bt_11um - Input%Bt_10um
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_10um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('btd1112')
                 Value_Dim = Input%Bt_11um - Input%Bt_12um
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_12um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('fmft')
                 Value_Dim = (Input%Bt_11um - Input%Bt_12um) - &
                         (Input%Bt_11um_Clear - Input%Bt_12um_Clear) * &
                         (Input%Bt_11um - 260.0) / &
                         (Input%Bt_11um_Clear - 260.0)
                 IF (Input%Bt_11um <= 260.0 .or. Input%Bt_11um_Clear <= 260.0) &
                    Value_Dim = Input%Bt_11um - Input%Bt_12um
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_12um == 0) Value_Dim = MISSING_VALUE_REAL4

             CASE('logcod065') ! TODO DON'T HAVE COD calculation yet
                 Value_Dim = MISSING_VALUE_REAL4 !Input%Logcod065
             CASE('logcod138') ! TODO DON'T HAVE COD calculation yet
                 Value_Dim = MISSING_VALUE_REAL4 !Input%Logcod138
             CASE('logcod160') ! TODO DON'T HAVE COD calculation yet
                 Value_Dim = MISSING_VALUE_REAL4 !Input%Logcod160
             CASE('refl047')
                 Value_Dim = Input%Ref_047um
             CASE('refl065')
                 Value_Dim = Input%Ref_063um
             CASE('refldnb')
                 Value_Dim = Input%Ref_Lunar
             CASE('logcoddnb') ! TODO DON'T HAVE COD calculation yet
                 Value_Dim = MISSING_VALUE_REAL4 !Input%Log_Cld_Opd_DNB
             CASE('refl065cv')
                 IF (Input%Ref_063um_Std_Sub /= MISSING_VALUE_REAL4) THEN
                    Value_Dim = Input%Ref_063um_Std_Sub / Input%Ref_063um
                 ENDIF
             CASE('refl065clr')
                 Value_Dim = Input%Ref_063um_Clear
             CASE('rgct')
                 Value_Dim = Input%Ref_063um - Input%Ref_063um_Clear
             CASE('drefl065clr')
                 Value_Dim = Input%Ref_063um - Input%Ref_063um_Clear
                 IF (Input%Ref_063um_Max_Sub /= MISSING_VALUE_REAL4) THEN
                    Value_Dim = Input%Ref_063um_Max_Sub - Input%Ref_063um_Clear
                 ENDIF
             CASE('drefl065maxsubclr')
                 IF (Input%Ref_063um_Max_Sub /= MISSING_VALUE_REAL4) THEN
                    Value_Dim = Input%Ref_063um_Max_Sub - Input%Ref_063um_Clear
                 ENDIF
             CASE('drefl065min')
                 IF (Input%Ref_063um_Min_Sub /= MISSING_VALUE_REAL4) THEN
                    Value_Dim = Input%Ref_063um - Input%Ref_063um_Min
                 ENDIF
             CASE('drefl065min3x3')
                 Value_Dim = Input%Ref_063um - Input%Ref_063um_Min
             CASE('drefl065minsub')
                 IF (Input%Ref_063um_Min_Sub == MISSING_VALUE_REAL4) CYCLE
                 Value_Dim = Input%Ref_063um - Input%Ref_063um_Min_Sub
                 if (Input%Chan_On_063um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('drefl065maxsub')
                 IF (Input%Ref_063um_Max_Sub == MISSING_VALUE_REAL4) CYCLE
                 Value_Dim = Input%Ref_063um_Max_Sub - Input%Ref_063um
                 if (Input%Chan_On_063um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('drefl065maxminsub')
                IF (Input%Ref_063um_Min_Sub == MISSING_VALUE_REAL4) CYCLE
                IF (Input%Ref_063um_Max_Sub == MISSING_VALUE_REAL4) CYCLE
                 Value_Dim = Input%Ref_063um_Max_Sub - Input%Ref_063um_Min_Sub
             CASE('logdrefl065maxminsub')
                IF (Input%Ref_063um_Min_Sub == MISSING_VALUE_REAL4) CYCLE
                IF (Input%Ref_063um_Max_Sub == MISSING_VALUE_REAL4) CYCLE
                 Value_Dim = Input%Ref_063um_Max_Sub - Input%Ref_063um_Min_Sub
                 if (Value_Dim > 0) Value_Dim = alog10(Value_Dim)
                 if (Value_Dim <= 0) Value_Dim = -2.0
             CASE('refl065maxsub')
                IF (Input%Ref_063um_Max_Sub == MISSING_VALUE_REAL4) CYCLE
                Value_Dim = Input%Ref_063um_Max_Sub
                if (Input%Chan_On_063um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('refl065minsub')
                IF (Input%Ref_063um_Min_Sub == MISSING_VALUE_REAL4) CYCLE
                 Value_Dim = Input%Ref_063um_Min_Sub
                 if (Input%Chan_On_063um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('refl065std')
                 Value_Dim = Input%Ref_063um_Std
                 if (Input%Chan_On_063um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('refrat086065')
                 Value_Dim = Input%Ref_086um / Input%Ref_063um
                 if (Input%Chan_On_086um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('refrat138065')
                 Value_Dim = Input%Ref_138um / Input%Ref_063um
                 if (Input%Chan_On_138um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('refl138')
                 Value_Dim = Input%Ref_138um
                 if (Input%Chan_On_138um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('refl160')
                 Value_Dim = Input%Ref_160um
                 if (Input%Chan_On_160um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('refl38')
                 Value_Dim = Input%Ref_375um
                 if (Input%Chan_On_138um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('logrefl065std')
                 Value_Dim = Input%Log_Ref_063um_Std
             CASE('ndsi')
                 Value_Dim = (Input%Ref_063um - Input%Ref_160um) / &
                             (Input%Ref_063um + Input%Ref_160um)
                 if (Input%Chan_On_160um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('ndsi_dnb_38')
                 Value_Dim = (Input%Ref_Lunar - Input%Ref_375um) / &
                             (Input%Ref_Lunar + Input%Ref_375um)
                 if (Input%Chan_On_375um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('bt6710covar')
                 Value_Dim = Input%Bt_10um_Bt_67um_Covar
                 if (Input%Chan_On_67um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_10um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE('bt6711covar')
                 Value_Dim = Input%Bt_11um_Bt_67um_Covar
                 if (Input%Chan_On_67um == 0) Value_Dim = MISSING_VALUE_REAL4
                 if (Input%Chan_On_11um == 0) Value_Dim = MISSING_VALUE_REAL4
             CASE default
                 PRINT *,"unknown classifier DIMENSION name = ",trim(Dim_Name)
                 Class_Use_Flag(Class_Idx) = 0
                 CYCLE
          END SELECT

      ! - SELECT input 
          SELECT CASE(i)
             CASE(1)
                 X_Dim = Value_Dim
             CASE(2)
                 Y_Dim = Value_Dim
             CASE(3)
                 Z_Dim = Value_Dim
          END SELECT

   ENDDO ! loop over class rank

   ! - set prior
   ! - Note that if the prior LUT is not read in, Prior_Yes_Temp is set to 
   !   Missing. Here is where it is then set for the given test - WCS3
   IF (Prior_Yes_Temp /= MISSING_VALUE_REAL4) THEN
       Prior_Prob_Clear = 1.0 - Prior_Yes_Temp
       Prior_Prob_Water = 0.5 * Prior_Yes_Temp
       Prior_Prob_Ice = 0.5 * Prior_Yes_Temp
   ELSE
       Prior_Prob_Clear = Lut(Class_Idx)%Cloud_fraction(Output%Sfc_Idx)
       Prior_Prob_Water = Lut(Class_Idx)%Water_Fraction(Output%Sfc_Idx)
       Prior_Prob_Ice = Lut(Class_Idx)%Ice_Fraction(Output%Sfc_Idx)
   ENDIF
   
   ! - compute for this class
   call GET_PROB_MASK_PHASE(X_Dim, Y_Dim, Z_Dim, Input%Senzen, &
                            Input%Solzen, Input%Lunzen, &
                            Input%Glintzen, Input%Lunglintzen, &
                            Input%Lunar_Oceanic_Glint_Mask, &
                            Input%Moon_Illum_Frac, Input%City_Mask, &
                            Input%Scatzen, Input%Path_Tpw, &
                            Input%Sfc_Temp, Input%Zsfc, Input%Zsfc_Std, Input%Land_Class, &
                            Input%Snow_Class, Input%Oceanic_Glint_Mask, &
                            Input%Coastal_Mask, Output%Sfc_Idx, Class_Idx, Lut, &
                            MISSING_VALUE_INT1, MISSING_VALUE_REAL4, Z1, Z2, Z3, Z4)

   Clear_Cond_Ratio(Class_Idx) = Z1
   Water_Cond_Ratio(Class_Idx) = Z2
   Ice_Cond_Ratio(Class_Idx) = Z3
   Obs_Prob(Class_Idx) = Z4


   !--------------------------------------
   ! question - we REALly DO not need to save the Cond_Ratio's for each Class
   !--------------------------------------

   !--- combine this class into the running product of all classes
   IF (Clear_Cond_Ratio(Class_Idx) /= MISSING_VALUE_REAL4) &
               R_Clear = min(R_Max, R_Clear * Clear_Cond_Ratio(Class_Idx))
   IF (Water_Cond_Ratio(Class_Idx) /= MISSING_VALUE_REAL4) &
               R_Water = min(R_Max, R_Water * Water_Cond_Ratio(Class_Idx))
   IF (Ice_Cond_Ratio(Class_Idx) /= MISSING_VALUE_REAL4) &
               R_Ice = min(R_Max, R_Ice * Ice_Cond_Ratio(Class_Idx))

   !--- store class value
   IF (Clear_Cond_Ratio(Class_Idx) /= MISSING_VALUE_REAL4) THEN
      Post_Prob_Clear_By_Class = 1.0 / (1.0 + Clear_Cond_Ratio(Class_Idx) / &
                            Prior_Prob_Clear - Clear_Cond_Ratio(Class_Idx))

      ! - convert from clear to cloud prob
      Posterior_Cld_Probability_By_Class(Class_Idx) = 1.0 - Post_Prob_Clear_By_Class

   ENDIF
   
   

ENDDO ! loop over classifiers

!------------------------------------------------------------------------------
! ---  compute posterior probs for this pixel
!------------------------------------------------------------------------------

! - turn ratios into probabilities via post_prob = 1.0  / ( 1.0 + R/Prior_Yes - R)
Post_Prob_Clear = 1.0 / (1.0 + R_Clear/Prior_Prob_Clear - R_Clear)
Output%Posterior_Water_Probability = 1.0 / (1.0 + R_Water/Prior_Prob_Water - R_Water)
Output%Posterior_Ice_Probability =   1.0 / (1.0 + R_Ice/Prior_Prob_Ice - R_Ice)

!--- constrain cold clouds to be ice
if (Input%Topa < 240.0 .and. Input%Topa /= MISSING_VALUE_REAL4) then
   Output%Posterior_Ice_Probability = Output%Posterior_Ice_Probability + Output%Posterior_Water_Probability
   Output%Posterior_Water_Probability = 0.0
endif

!--- normalize these probabilties so that they sum to 1  - won't happen IF prior is not from lut
Post_Sum =  Post_Prob_clear + Output%Posterior_Water_Probability + Output%Posterior_Ice_Probability
Post_Prob_Clear = Post_Prob_Clear / Post_Sum
Output%Posterior_Water_Probability = Output%Posterior_Water_Probability / Post_Sum
Output%Posterior_Ice_Probability = Output%Posterior_Ice_Probability / Post_Sum

! --- make a cloud probability from the complement of the clear probability
Output%Posterior_Cld_Probability = 1.0 - Post_Prob_Clear

! --- make the binary mask
IF (Output%Posterior_Cld_Probability < Mask_Thresh%Prob_Clear_Prob_Cloudy_Thresh(Output%Sfc_Idx) .and. &
    Output%Posterior_Cld_Probability /= MISSING_VALUE_REAL4) Output%Cld_Mask_Binary = symbol%CLEAR_BINARY

IF (Output%Posterior_Cld_Probability >= Mask_Thresh%Prob_Clear_Prob_Cloudy_Thresh(Output%Sfc_Idx) .and. &
    Output%Posterior_Cld_Probability /= MISSING_VALUE_REAL4) Output%Cld_Mask_Binary = symbol%CLOUDY_BINARY

! --- make the 4-level mask
IF (Output%Posterior_Cld_Probability >= 0.0 .and. &
    Output%Posterior_Cld_Probability < Mask_Thresh%Conf_Clear_Prob_Clear_Thresh(Output%Sfc_Idx)) &
        Output%Cld_Mask_Bayes = symbol%CLEAR


IF (Output%Posterior_Cld_Probability >= Mask_Thresh%Conf_Clear_Prob_Clear_Thresh(Output%Sfc_Idx) .and. &
    Output%Posterior_Cld_Probability < Mask_Thresh%Prob_Clear_Prob_Cloudy_Thresh(Output%Sfc_Idx)) &
        Output%Cld_Mask_Bayes = symbol%PROB_CLEAR

IF (Output%Posterior_Cld_Probability >= Mask_Thresh%Prob_Clear_Prob_Cloudy_Thresh(Output%Sfc_Idx) .and. &
    Output%Posterior_Cld_Probability < Mask_Thresh%Prob_Cloudy_Conf_Cloudy_Thresh(Output%Sfc_Idx)) &
        Output%Cld_Mask_Bayes = symbol%PROB_CLOUDY

IF (Output%Posterior_Cld_Probability >= Mask_Thresh%Prob_Cloudy_Conf_Cloudy_Thresh(Output%Sfc_Idx) .and. &
    Output%Posterior_Cld_Probability <= 1.0) &
        Output%Cld_Mask_Bayes = symbol%CLOUDY

!--- Thermal Uniformity Test Filter (Needs to Use 10.4 when necessary)
IF (Input%Bt_11um_Std .ne. MISSING_VALUE_REAL4 .and. Input%Zsfc_Std .ne. MISSING_VALUE_REAL4) THEN
    Output%TUT = TUT_ROUTINE(Input%Coastal_Mask, Input%Bt_11um_Std, Input%Zsfc_Std, &
                         Mask_Thresh%Tut_Clear_Prob_Clear_Thresh(Output%Sfc_Idx))
ENDIF

!--- Reflectance Uniformity Test Filter
IF (Input%Solzen < Lut(1)%Rut_Solzen_Thresh .and. Input%Ref_063um .ne. MISSING_VALUE_REAL4) THEN
    Is_Land = 1
    IF (Output%Sfc_Idx == 0 .or. Output%Sfc_Idx == 1) Is_Land = 0
    Output%RUT = RUT_ROUTINE(Input%Coastal_Mask, Is_Land, Input%Ref_063um, Input%Ref_063um_Std, &
                             Mask_Thresh%Rut_Clear_Prob_Clear_Thresh(Output%Sfc_Idx))
ENDIF

!--- apply uniformity
!IF (Output%Cld_Mask_Bayes == symbol%CLEAR .and. (Output%TUT == 1 .or. Output%RUT == 1)) &
!       Output%Cld_Mask_Bayes = symbol%PROB_CLEAR


!----------------------------------------------------------------------------------
!--- set some flags to control processing - Lets move this to a subroutine
!----------------------------------------------------------------------------------
call SET_NONCLOUD_FLAGS(Input, Output, symbol, ABI_Use_104um_Flag, &
                        Oceanic_Glint_Flag,Coastal_Flag,Solar_Contam_Flag, &
                        Day_063_Flag, Day_063_Spatial_Flag, Lunar_Spatial_Flag, &
                        Night_Lunar_Flag, Lunar_Forward_Scattering_Flag, &
                        Lunar_Oceanic_Glint_Flag,Day_375_Flag,Night_375_Flag, &
                        Mountain_Flag,Forward_Scattering_Flag, &
                        Cold_Scene_Flag, Snow_Flag, &
                        Dry_Scene_Flag, City_Flag, Use_104_Flag)

!----------------------------------------------------------------------------------
!--- populate elements of Cld_Flags with processing flags
!----------------------------------------------------------------------------------
Cld_Flags(2) = Day_063_Flag            ;    Cld_Flag_Bit_Depth(2) = 1
Cld_Flags(3) = Day_063_Spatial_Flag    ;    Cld_Flag_Bit_Depth(3) = 1
Cld_Flags(4) = Day_375_Flag            ;    Cld_Flag_Bit_Depth(4) = 1
Cld_Flags(5) = Night_375_Flag          ;    Cld_Flag_Bit_Depth(5) = 1
Cld_Flags(6) = Solar_Contam_Flag       ;    Cld_Flag_Bit_Depth(6) = 1
Cld_Flags(7) = Coastal_Flag            ;    Cld_Flag_Bit_Depth(7) = 1
Cld_Flags(8) = Mountain_Flag           ;    Cld_Flag_Bit_Depth(8) = 1
Cld_Flags(9) = Forward_Scattering_Flag ;    Cld_Flag_Bit_Depth(9) = 1
Cld_Flags(10) = Snow_Flag              ;    Cld_Flag_Bit_Depth(10) = 1
Cld_Flags(11) = Cold_Scene_Flag        ;    Cld_Flag_Bit_Depth(11) = 1
Cld_Flags(12) = Oceanic_Glint_Flag     ;    Cld_Flag_Bit_Depth(12) = 1
Cld_Flags(13) = Spare_Value            ;    Cld_Flag_Bit_Depth(13) = 1 ! saved for smoke
Cld_Flags(14) = Spare_Value            ;    Cld_Flag_Bit_Depth(14) = 1 ! saved for dust
Cld_Flags(15) = Spare_Value            ;    Cld_Flag_Bit_Depth(15) = 1 ! saved for shadow 
Cld_Flags(16) = Spare_Value            ;    Cld_Flag_Bit_Depth(16) = 1 ! saved for fire
Cld_Flags(17) = Output%Sfc_Idx         ;    Cld_Flag_Bit_Depth(17) = 3
Cld_Flags(18) = Spare_Value            ;    Cld_Flag_Bit_Depth(18) = 1 ! saved for thin cirrus
Cld_Flags(19) = Use_104_Flag           ;    Cld_Flag_Bit_Depth(19) = 1
Cld_Flags(20:54) = Spare_Value         ;    Cld_Flag_Bit_Depth(20:54) = 1

!------------------------------------------------------------------------------------------------------------
!--- compute probabilities for each class alone - used for flags - not
!--- needed for mask or final probability - it should remain optional
!------------------------------------------------------------------------------------------------------------
IF (Do_By_Class_Flag == symbol%YES) THEN

   DO Class_Idx = 1, Nclass

      IF (Posterior_Cld_Probability_By_Class(Class_Idx) < Mask_Thresh%Prob_Clear_Prob_Cloudy_Thresh(Output%Sfc_Idx)) THEN
         Cld_Flags(NUMBER_OF_NONCLOUD_FLAGS+Class_Idx) = symbol%CLEAR_BINARY

      ELSEIF (Posterior_Cld_Probability_By_Class(Class_Idx) >= Mask_Thresh%Prob_Clear_Prob_Cloudy_Thresh(Output%Sfc_Idx)) THEN
         Cld_Flags(NUMBER_OF_NONCLOUD_FLAGS+Class_Idx) = symbol%CLOUDY_BINARY
      ENDIF

   ENDDO


   !-------------------------------------------------------------------
   ! Pack Bytes
   !-------------------------------------------------------------------
   call PACK_BITS_INTO_BYTES(Cld_Flags,Cld_Flag_Bit_Depth,Output%Cld_Flags_Packed(:))

ENDIF! optional results by test

ENDIF !IF Invalid_Data_Mask

! --- DEALLOCATE variables
IF (ALLOCATED(Class_Use_Flag)) DEALLOCATE (Class_Use_Flag)
IF (ALLOCATED(Clear_Cond_Ratio)) DEALLOCATE (Clear_Cond_Ratio)
IF (ALLOCATED(Water_Cond_Ratio)) DEALLOCATE (Water_Cond_Ratio)
IF (ALLOCATED(Ice_Cond_Ratio)) DEALLOCATE (Ice_Cond_Ratio)
IF (ALLOCATED(Obs_Prob)) DEALLOCATE (Obs_Prob)
IF (ALLOCATED(Posterior_Cld_Probability_By_Class)) DEALLOCATE (Posterior_Cld_Probability_By_Class)
IF (ALLOCATED(Chan_On)) DEALLOCATE (Chan_On)
IF (ALLOCATED(Chan_Wvl)) DEALLOCATE (Chan_Wvl)
IF (ALLOCATED(Class_Use_Flag)) DEALLOCATE (Class_Use_Flag)


END SUBROUTINE NBM_CLOUD_MASK_ALGORITHM


!====================================================================
! SUBROUTINE Name: COMPUTE_BAYES_SFC_TYPE
!
! Function:
!   Computes the bayesian surface TYPE given the ancillary sfc data
!
!====================================================================
 FUNCTION COMPUTE_BAYES_SFC_TYPE(Land_Temp, Coast_Temp, Snow_Temp, Sfc_Type_Temp, &
                                 Lat_Temp, Lon_Temp, Sst_Back_Uni_Temp,&
                                 Emiss_Sfc_375um_Temp, symbol) &
                                 result(Bayes_Mask_Sfc_Type_Temp)
   INTEGER(KIND=INT1), INTENT(IN):: Land_Temp
   INTEGER(KIND=INT1), INTENT(IN):: Coast_Temp
   INTEGER(KIND=INT1), INTENT(IN):: Snow_Temp
   INTEGER(KIND=INT1), INTENT(IN):: Sfc_Type_Temp
   REAL(KIND=REAL4), INTENT(IN):: Lat_Temp
   REAL(KIND=REAL4), INTENT(IN):: Lon_Temp
   REAL(KIND=REAL4), INTENT(IN):: Sst_Back_Uni_Temp
   REAL(KIND=REAL4), INTENT(IN):: Emiss_Sfc_375um_Temp
   TYPE(symbol_naive_bayesian), INTENT(IN) :: symbol

   INTEGER(KIND=INT4):: Bayes_Mask_Sfc_Type_Temp


   IF (Land_Temp == symbol%LAND) THEN
           Bayes_Mask_Sfc_Type_Temp = 3
   ELSE
           Bayes_Mask_Sfc_Type_Temp = 1
   ENDIF

   !--- #2 - Shallow Ocean
   IF ((Land_Temp == symbol%MODERATE_OCEAN) .or. &
       (Land_Temp == symbol%DEEP_INLAND_WATER) .or. &
       (Land_Temp == symbol%SHALLOW_INLAND_WATER) .or. &
       (Land_Temp == symbol%SHALLOW_OCEAN)) THEN
           Bayes_Mask_Sfc_Type_Temp = 2
   ENDIF
   IF ((Land_Temp /= symbol%LAND) .and. &
       (Sst_Back_Uni_Temp > 0.5)) THEN
           Bayes_Mask_Sfc_Type_Temp = 2
   ENDIF

   !--- #3 Unfrozen_Land 
   IF ((Land_Temp == symbol%LAND) .or. &
       (Land_Temp == symbol%COASTLINE) .or. &
       (Coast_Temp == symbol%YES) .or. &
       (Land_Temp == symbol%EPHEMERAL_WATER)) THEN

           Bayes_Mask_Sfc_Type_Temp = 3

   ENDIF

   !--- #4 - Snow Covered Land
   IF ((Lat_Temp > -60.0) .and. (Snow_Temp == symbol%SNOW)) THEN
           Bayes_Mask_Sfc_Type_Temp = 4
   ENDIF

   !--- #5 - Arctic
   IF ((Lat_Temp >= 0.0) .and. (Snow_Temp == symbol%SEA_ICE)) THEN
           Bayes_Mask_Sfc_Type_Temp = 5
   ENDIF

   !--- #6 - Antarctic
   IF ((Lat_Temp <= -60.0) .and. (Snow_Temp == symbol%SNOW)) THEN
           Bayes_Mask_Sfc_Type_Temp = 6
   ENDIF
   IF ((Lat_Temp <= 0.0) .and. (Snow_Temp == symbol%SEA_ICE)) THEN
           Bayes_Mask_Sfc_Type_Temp = 6
   ENDIF
   IF ((Lat_Temp >= 60.0) .and.  &
       (Lon_Temp > -75.0) .and. (Lon_Temp < -10.0) .and. &
       (Land_Temp == symbol%LAND .or. Land_Temp == symbol%COASTLINE) .and. &
       (Snow_Temp == symbol%SNOW)) THEN
           Bayes_Mask_Sfc_Type_Temp = 6
   ENDIF
         
   !--- #7 - Desert
   IF ( (Emiss_Sfc_375um_Temp < 0.90 ) .and. (abs(Lat_Temp) < 60.0) .and. &
       ((Sfc_Type_Temp == symbol%OPEN_SHRUBS_SFC) .or. (Sfc_Type_Temp == symbol%BARE_SFC)) ) THEN
           Bayes_Mask_Sfc_Type_Temp = 7
   ENDIF
         
         
   IF ( Bayes_Mask_Sfc_Type_Temp == 3 .and.  &
        Emiss_Sfc_375um_Temp < 0.93  .and.  &
        abs(Lat_Temp) < 60.0 .and. &
       ((Sfc_Type_Temp == symbol%OPEN_SHRUBS_SFC) .or.  &
        (Sfc_Type_Temp == symbol%CLOSED_SHRUBS_SFC) .or. &
        (Sfc_Type_Temp == symbol%GRASSES_SFC) .or.  &
        (Sfc_Type_Temp == symbol%BARE_SFC)) ) THEN
           Bayes_Mask_Sfc_Type_Temp = 7
   ENDIF
        
  RETURN 
   
 END FUNCTION COMPUTE_BAYES_SFC_TYPE

!====================================================================
! Function Name: TUT_Routine
!
! Function:
!   Thermal Uniformity Routine (TUT) Test
!
! Description:
!   Computes the Thermal Uniformity Routine (TUT) test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = TUT_Routine (&
!                              Is_Land, &
!                              Is_Coast, &
!                              BT_Chn14_Stddev_3x3, &
!                              Sfc_Hgt_Stddev_3x3)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   Standard Deviation of the 11 micron BT over a 3x3 box
!   Standard Deviation of the surface height over a 3x3 box
!
! Outputs: 
!   Function RETURNs pass (sym%YES) or fail (sym%NO) result of the test via
!   Test_Result
!
! DepENDencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================
FUNCTION TUT_ROUTINE (Is_Coast, &
                      Bt_Chn14_Stddev_3x3, &
                      Sfc_Hgt_Stddev_3x3, &
                      BT_CHN14_CLR_UNI_THRESH) &
                      result(Test_Result)

   INTEGER(KIND=INT1), INTENT(IN) :: Is_Coast
   REAL, INTENT(IN) :: Bt_Chn14_Stddev_3x3
   REAL, INTENT(IN) :: Sfc_Hgt_Stddev_3x3
   REAL, INTENT(IN) :: BT_CHN14_CLR_UNI_THRESH

   INTEGER :: Test_Result
   REAL :: Test_Threshold

   Test_Result = 0
   Test_Threshold = 0.0
   IF (Is_Coast == 0) THEN

         !
         !7K/km is the adiabatic lapse rate
         !
         Test_Threshold = 3.0 * 7.0*Sfc_Hgt_Stddev_3x3/1000.0

           IF (Bt_Chn14_Stddev_3x3 > BT_CHN14_CLR_UNI_THRESH + Test_Threshold) THEN
             Test_Result = 1
           ENDIF

   ENDIF

   RETURN

END FUNCTION TUT_ROUTINE

!====================================================================
! Function Name: RUT_Routine
!====================================================================
FUNCTION RUT_ROUTINE (Is_Coast, &
                      Is_Land,  &
                      Refl_Chn2_Clear, &
                      Refl_Chn2_Stddev_3x3, &
                      Rut_Thresh) &
                      result(Test_Result)

   INTEGER(KIND=INT1), INTENT(IN) :: Is_Coast
   INTEGER, INTENT(IN) :: Is_Land
   REAL, INTENT(IN) :: Refl_Chn2_Clear
   REAL, INTENT(IN) :: Refl_Chn2_Stddev_3x3
   REAL, INTENT(IN) :: RUT_THRESH

   INTEGER :: Test_Result
   REAL :: Test_Threshold

       Test_Result =  0

       IF (Is_Coast == 0) THEN

              !--- compute threshold
              IF (Is_Land == 1) THEN
                 Test_Threshold = MAX(0.5,Refl_Chn2_Clear * Rut_Thresh)
              ELSE
                 Test_Threshold = RUT_THRESH
              ENDIF

              !--- apply test
              IF (Refl_Chn2_Stddev_3x3 > Test_Threshold) THEN
                 Test_Result = 1
              ENDIF

       ENDIF

       RETURN

END FUNCTION RUT_ROUTINE


!------------------------------------------------------------------------
! SUBROUTINE PACK_BITS_INTO_BYTES(input_bits,bit_depth,output_bytes)
!
! Routines to pack individual bytes into a single byte
!
! input:
! input_bits - vector of bytes to be packed into output_byte
! bit_start - vector of bit starting positions (1-7) for each input byte
! bit_depth - vector of bit depths (1-7) for each input byte (total can not exceed 8)
!
! output: 
! output_byte - byte variable that holds the bit values of the input_bits - can be i1 or i2
!
! local
!  n_in - number of elements of input vectors
!  i_in - index of input_bits (1-n_in)
!  n_out - number of elements of output vectors
!  i_out - index of output_bytes (1-n_out)
!
! Note:
! 1.  IF the input byte has information in bits greater THEN bit depth - they are removed 
!
!
! Example, pack an input_byte wth  bit_start = 2 and bit depth 3 
!
! input byte
!           x x x
! _ _ _ _ _ _ _ _
!
! result of first ishft
! x x x
! _ _ _ _ _ _ _ _
!
! result of second ishft
!       x x x
! _ _ _ _ _ _ _ _
!
! Author: Andrew Heidinger
!
! Version History:  
! February 2006 - Created
!-----------------------------------------------------------------------------------

!--- This Version packs into one byte words
SUBROUTINE PACK_BITS_INTO_BYTES (input_bits,bit_depth,output_bytes)
    INTEGER(KIND=INT1), DIMENSION(:), INTENT(IN):: input_bits
    INTEGER(KIND=INT1), DIMENSION(:), INTENT(IN):: bit_depth
    INTEGER(KIND=INT1), DIMENSION(:), INTENT(OUT):: output_bytes

    INTEGER(KIND=INT1):: bit_start, bit_END, bit_offset
    INTEGER(KIND=INT1):: temp_byte
    INTEGER:: n_in,i_in,n_out,i_out
    INTEGER, PARAMETER:: word_bit_depth = 8

    !--- determine size of vectors
    n_in = size(input_bits)
    n_out = size(output_bytes)

    !--- reset output byte
    output_bytes = 0

    !--- initialize
    bit_offset = 0
    bit_start = 0
    bit_END = 0
    i_out = 1

    !--- loop through input bytes
    DO i_in = 1, n_in

     !--- determine starting and ENDing bit locations
     bit_start = bit_offset + 1
     bit_END = bit_start + bit_depth(i_in) - 1

     !--- determine IF this input byte will fit on current output byte, IF not go to next
     IF (bit_END > word_bit_depth) THEN
      i_out = i_out + 1
      bit_offset = 0
      bit_start = bit_offset + 1
      bit_END = bit_start + bit_depth(i_in) - 1
     ENDIF

     !--- check for exceeding the space allowed for the packed bytes
     IF (i_out > n_out) THEN
       PRINT *, "ERROR: Insufficient space for bit packing ", i_out, n_out
       RETURN
     ENDIF

     !--- place input byte into correct position
     temp_byte =0
     temp_byte = ishft(input_bits(i_in),word_bit_depth-bit_depth(i_in)) !first ishft
     temp_byte = ishft(temp_byte,bit_END - word_bit_depth)              !second ishft

     !--- modify output byte
     output_bytes(i_out) = output_bytes(i_out) + temp_byte

     !--- update bit offset
     bit_offset = bit_offset + bit_depth(i_in)

    ENDDO

END SUBROUTINE  PACK_BITS_INTO_BYTES
!-------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------
SUBROUTINE SET_NONCLOUD_FLAGS(Input, Output, symbol, ABI_Use_104um_Flag, &
                        Oceanic_Glint_Flag,Coastal_Flag,Solar_Contam_Flag, &
                        Day_063_Flag, Day_063_Spatial_Flag, Lunar_Spatial_Flag, &
                        Night_Lunar_Flag, Lunar_Forward_Scattering_Flag, &
                        Lunar_Oceanic_Glint_Flag,Day_375_Flag,Night_375_Flag, &
                        Mountain_Flag,Forward_Scattering_Flag, &
                        Cold_Scene_Flag, Snow_Flag, &
                        Dry_Scene_Flag, City_Flag, Use_104_Flag)

    TYPE(symbol_naive_bayesian), INTENT(IN) :: Symbol
    TYPE(mask_input), INTENT(IN) :: Input 
    TYPE(mask_output), INTENT(IN) :: Output
    LOGICAL, INTENT(IN):: ABI_Use_104um_Flag
    INTEGER, INTENT(OUT):: Oceanic_Glint_Flag
    INTEGER, INTENT(OUT):: Lunar_Oceanic_Glint_Flag
    INTEGER, INTENT(OUT):: Coastal_Flag
    INTEGER, INTENT(OUT):: Mountain_Flag
    INTEGER, INTENT(OUT):: Day_063_Flag
    INTEGER, INTENT(OUT):: Day_063_Spatial_Flag
    INTEGER, INTENT(OUT):: Night_Lunar_Flag
    INTEGER, INTENT(OUT):: Lunar_Spatial_Flag
    INTEGER, INTENT(OUT):: Day_375_Flag
    INTEGER, INTENT(OUT):: Night_375_Flag
    INTEGER, INTENT(OUT):: Forward_Scattering_Flag
    INTEGER, INTENT(OUT):: Solar_Contam_Flag
    INTEGER, INTENT(OUT):: Lunar_Forward_Scattering_Flag
    INTEGER, INTENT(OUT):: Snow_Flag
    INTEGER, INTENT(OUT):: Cold_Scene_Flag
    INTEGER, INTENT(OUT):: Dry_Scene_Flag
    INTEGER, INTENT(OUT):: City_Flag
    INTEGER (KIND=INT1):: Use_104_Flag
    REAL (KIND=REAL4):: Airmass


Oceanic_Glint_Flag = Input%Oceanic_Glint_Mask
! --- Added by Denis B. Glint over land was set to glint IF glint_angle is < 45.0
IF (Input%Glintzen < Glint_Zen_Thresh .and. &
   (Input%Land_Class == 1 .or. Input%Land_Class == 2)) &
  Oceanic_Glint_Flag = 1

Coastal_Flag = Input%Coastal_Mask
Solar_Contam_Flag = Input%Solar_Contamination_Mask

!--- compute airmass
AirMass = 1.0/cos(Input%Solzen*dtor) + 1.0 / cos(Input%Senzen*dtor)

!--- set day flag for 0.63 micron reflectance gross test
IF ((Input%Solzen > Reflectance_Gross_Solzen_Thresh) .or.  &
    (AirMass > Reflectance_Gross_Airmass_Thresh)) THEN
   Day_063_Flag = symbol%NO
ELSE
   Day_063_Flag = symbol%YES
ENDIF

!--- set day flag for 0.63 micron reflectance spatial test
IF (Input%Solzen > Reflectance_Spatial_Solzen_Thresh) THEN
   Day_063_Spatial_Flag = symbol%NO
ELSE
   Day_063_Spatial_Flag = symbol%YES
ENDIF

Lunar_Spatial_Flag = symbol%NO
Night_Lunar_Flag = symbol%NO
Lunar_Forward_Scattering_Flag = symbol%NO

IF (Input%Chan_On_DNB == symbol%YES) THEN
   Lunar_Oceanic_Glint_Flag = Input%Lunar_Oceanic_Glint_Mask
   Lunar_Spatial_Flag = symbol%YES
   IF (Input%Lunzen > Reflectance_Spatial_Solzen_Thresh) THEN
      Lunar_Spatial_Flag = symbol%NO
   ENDIF

   Night_Lunar_Flag = symbol%YES
   IF (Input%Lunzen > Reflectance_Gross_Lunzen_Thresh) THEN
      Night_Lunar_Flag = symbol%NO
   ENDIF

   Lunar_Forward_Scattering_Flag = symbol%NO
   IF (Input%Lunscatzen < Forward_Scatter_Scatzen_Max_Thresh .and. &
       Input%Lunzen < Forward_Scatter_Solzen_Max_Thresh) THEN
      Lunar_Forward_Scattering_Flag = symbol%YES
   ENDIF
ENDIF

IF ((Input%Solzen > Emiss_375um_Day_Solzen_Thresh) .or. &
    (AirMass > Reflectance_Gross_Airmass_Thresh)) THEN
   Day_375_Flag = symbol%NO
ELSE
   Day_375_Flag = symbol%YES
ENDIF

IF (Input%Solzen < Emiss_375um_Night_Solzen_Thresh) THEN
   Night_375_Flag = symbol%NO
ELSE
   Night_375_Flag = symbol%YES
ENDIF

IF (Output%Sfc_Idx /= 6 .and. Input%Zsfc > 2000.0) THEN
   Mountain_Flag = symbol%YES
ELSE
   Mountain_Flag = symbol%NO
ENDIF

IF (Input%Scatzen < Forward_Scatter_Scatzen_Max_Thresh .and. &
    Input%Solzen < Forward_Scatter_Solzen_Max_Thresh) THEN
   Forward_Scattering_Flag = symbol%YES
ELSE
   Forward_Scattering_Flag = symbol%NO
ENDIF

Snow_Flag = symbol%NO
IF (Input%Snow_Class == symbol%SNOW .or. Input%Snow_Class == symbol%SEA_ICE) THEN
   Snow_Flag = symbol%YES
ENDIF

Cold_Scene_Flag = symbol%NO
IF (Input%Sfc_Temp < Tsfc_Cold_Scene_Thresh) THEN
   Cold_Scene_Flag = symbol%YES
ENDIF

Dry_Scene_Flag = symbol%NO
IF (Input%Path_Tpw < Path_Tpw_Dry_Scene_Thresh) THEN
   Dry_Scene_Flag = symbol%YES
ENDIF

!--- City Flag of DNB Lunar Tests
City_Flag = symbol%NO
IF (Input%Chan_On_DNB == symbol%YES) THEN
   IF (Input%Rad_Lunar > Radiance_Lunar_City_Thresh) City_Flag = symbol%YES
ENDIF


!--- Set Use_104_Flag
Use_104_Flag = 0_INT1
IF (ABI_Use_104um_Flag) Use_104_Flag = 1_INT1

END SUBROUTINE SET_NONCLOUD_FLAGS

!====================================================================

END MODULE NBM_CLOUD_MASK_MODULE

