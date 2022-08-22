!$Id: baseline_cloud_height.f90 3514 2019-09-17 18:56:21Z yli $
MODULE Baseline_Cloud_Height
!CVS SRC: baseline_cloud_height.f90,v 1.26 2012/02/10 16:13:43 wstraka Exp $
!----------------------------------------------------------------------
! This module contains the a routine with the 13.3 micron channle added
! to the CLAVR-x Split-Window Cloud Algorithm
!
! Author: Andrew K. Heidinger, National Oceanic and Atmospheric Administration
! (c) This code is copyrighted by the author and all NOAA restrictions apply
!
! Name:
! Baseline_Cloud_Height
!
! Function:
! Derive the cloud temperature, emissivity and microphysical index (beta). 
! Cloud pressure and cloud height are 
! derived from the cloud temperature and the NWP profile
!
! Description: This routines uses the 11,12 and 13.3 micron channels
!              in an optimal estimation framework to estimate the 
!              clod height/temperature/pressure, emissivity and beta
!
! Reference: The Cloud Application Team Cloud Height ATBD
!
! Calling Sequence: Called after cloud mask and cloud type
!
! Inputs: All input passed through geocat structures 
!         (satellite, nwp, rtm and temporal)
!
! Outputs:   
!        Cld_Temp = out2(ialgo)%cldt = cloud temperature
!        Cld_Emiss_11um = out2(ialgo)%cldt = cloud 11 micron emissivity
!        beta_11_12 = out2(ialgo)%cldbeta_1112 = cloud 11/12 micron beta
!        Cld_Press = out2(ialgo)%cldp =  cloud pressure
!        Cld_Hgt = out2(ialgo)%cldp = cloud pressure
!        QF = out2(ialgo)%qflg1 = quality flags
!
! Dependencies:
!        GEOCAT satellite, rtm,  and nwp structures must be
!        populated for this segment
!
! Restrictions:
!
! History:
!   3/2007 - Andrew Heidinger - Created
!   1/2008 - Andrew Heidinger - Based beta relationships on ice
!                               aggregates based on CALIPSO analysis
!   2/2008 - Andrew Heidinger - Added lrc computations
!   3/2008 - Andrew Heidinger - Delivered to AIT 2nd time
!   6/2008 - Andrew Heidinger - Validated and Modified for CDR
!   9/2008 - William Straka - Renamed and Redelivered to AIT 
!   1/2009 - Andrew Heidinger - Modified for 3rd Delivery to AIT 
!   1/2009 - Andrew Heidinger - Added errors estimates to output as requested by AMV team
!   1/2009 - Andrew Heidinger - Added inversion flag to output as requested by AMV team
!   1/2010 - Andrew Heidinger - Added logic for multi-layer cloud processing
!   1/2010 - Andrew Heidinger - Added meta-data
!
!   3/2018 - Yue Li - Changed chan idx from 14-16 to 31-33 when calling
!   planck_rad(tem)_fast; commented LWP/IWP; added bad_pixel_mask in
!   Compute_Emiss_Tropo_Chn14; added bounds to NWP profiles when linked to local

! This routine uses a 1d-var retrieval approach as outlined in Rodgers (1976,2000)
!
! input to the 1d-var approach
!  Obs_Vector - the vector of observations
!  Obs_Vector(1) = bt11 - 11 micron brightness temperature
!  Obs_Vector(2) = bt11 - bt12 - the 11 - 12 micron (split window) brit. temp. difference
!  Obs_Vector(3) = bt11 - bt13 - the 11 - 13 micron brit. temp. difference
!  Apriori_Vector - the a apriori estimates of Retv_Vector
!  Sa - the error covariance matrix of Apriori_Vector
!  Sy - the error covariance of the Obs_Vector (included calibration, forward model)
!
! internal matrices in the 1d-var approach
!  Fwd_Model_Vector - the forward model estimate of Obs_Vector
!  Kernel - the kernel matrix - dFwd_Model_Vector/dRetv_Vector
!
! output from the 1d-var approach
!  Retv_Vector - the vector of parameters to retrieve
!  Retv_Vector(1) - the cloud temperature
!  Retv_Vector(2) - the 11 micron emissivity at nadir
!  Retv_Vector(3) - the beta_11_12 ratio for 11 and 12 microns
!  Sy - the error covariance of x
!
! the Parameter quality flags are determined as follows
!  3 - estimated error < 1/3 a priori error
!  2 - estimated error < 2/3 a priori error
!  1 - any other converged retrieval
!  0 - a failed or unattempted retrieval
!
! the Product quality flags are determined as follows
!  0 - a valid retrieval
!  1 - an invalid retrieval due to being a space pixel
!  2 - an invalid retrieval due to being outside the sensor zenith angle range
!  3 - an invalid retrieval due to having bad data
!  4 - an invalid retrieval due to being a clear/probably clear pixel
!  5 - an invalid retrieval due to being a pixel with invalid (Unknown) cloud type
!  6 - a failed or unattempted retrieval
!
! Meta Data Flags
!
! 1 - Cloud Height Attempted (0 = no / 1 = yes)
! 2 - Bias Correction Employed (0 = no / 1 = yes)
! 3 - Ice Cloud Retrieval (0 = no / 1 = yes)
! 4 - Local Radiatve Center Processing Used (0 = no / 1 = yes)
! 5 - Multi-layer Retrieval (0 = no / 1 = yes)
! 6 - Lower Cloud Interpolation Used (0 = no / 1 = yes)
! 7 - Boundary Layer Inversion Assumed  (0 = no / 1 = yes)
!
!---------------------------------------------------------------------- 

!--- declare module usage here
!USE ALGORITHM_MODULE_USAGE
use CONSTANTS_MOD
use CLAVRX_MESSAGE_MOD
use PIXEL_COMMON_MOD
use NWP_COMMON_MOD
use RTM_COMMON_MOD
use PLANCK_MOD


IMPLICIT NONE

PUBLIC:: Baseline_Cloud_Height_main

PRIVATE:: Spatial_Interp_Lower_Cld_Pos, &
          Compute_Emiss_Tropo_Chn14, &
          compute_spatial_uniformity, &
          destroy_spatial_uniformity, &
          gradient2d, &
          invert_3x3, &
          prof_lookup_using_t, &
          prof_lookup_using_p, &
          prof_lookup_using_z, &
          locate_Int32, &
          locate_float32

 INTERFACE locate
   MODULE PROCEDURE &
     Locate_Int32, &
     Locate_Float32
 END INTERFACE

!--- include parameters 
INCLUDE 'baseline_cloud_height.inc'

CONTAINS

!----------------------------------------------------------------------
! BEGIN SUBROUTINE
!----------------------------------------------------------------------
SUBROUTINE Baseline_Cloud_Height_Main()

 

  !--- local variable declaration

  !--- local indices
  INTEGER(KIND=INT4) :: Elem_Idx                !element number index
  INTEGER(KIND=INT4) :: Line_Idx                !line number index
  INTEGER(KIND=INT4) :: View_Zen_Idx            !viewing zenith angle bin for rtm calcs
  INTEGER(KIND=INT4) :: X_Idx_NWP               !longitude index for nwp fields
  INTEGER(KIND=INT4) :: Y_Idx_NWP               !latitude index for nwp fields
  INTEGER(KIND=INT4) :: Tropo_Idx_NWP           !level index of tropopause
  INTEGER(KIND=INT4) :: Sfc_Idx_NWP             !level index of surface
  INTEGER(KIND=INT4) :: NWP_Level_Idx           !level index
  INTEGER(KIND=INT4) :: Num_Level_NWP_Profile   !number of levels in profiles
  INTEGER(KIND=INT4):: Pass_Idx                 !index using for pass loop
  INTEGER(KIND=INT4):: Min_Val_Pass             !min value of ipass
  INTEGER(KIND=INT4):: Max_Val_Pass             !max value of ipass
  INTEGER(KIND=INT4):: Min_NWP_Level_Inver      !min ilev to look for inversions 
  INTEGER(KIND=INT4):: Max_NWP_Level_Inver      !max ilev to look for inversions 

  CHARACTER(LEN=100) :: Err_Message            !Error Message
  INTEGER(KIND=INT4) :: Error_Level            !Error Level

  !--- Allocation/Deallocation status

  INTEGER(KIND=int4) :: Alloc_Status

  !--- spatial uniformity matrices
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_11um_Mean_3X3   !mean bt11 over 3x3 box
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_11um_Max_3X3    !maximum bt11 over 3x3 box 
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_11um_Min_3X3    !minimum bt11 over 3x3 box
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_11um_Std_3X3    !std dev of  bt11 over 3x3 box
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_1112um_Mean_3X3 !mean bt11-bt12 over 3x3 box
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_1112um_Max_3X3  !maximum bt11-bt12 over 3x3 box
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_1112um_Min_3X3  !minimum bt11-bt12 over 3x3 box
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_1112um_Std_3X3  !std dev of bt11-bt12 over 3x3 box
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_1113um_Mean_3X3 !mean bt11-bt13 over 3x3 box
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_1113um_Max_3X3  !maximum bt11-bt13 over 3x3 box
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_1113um_Min_3X3  !minimum bt11-bt13 over 3x3 box
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE :: BT_1113um_Std_3X3  !std dev of bt11-bt13 over 3x3 box

  !-- Variables associated with the 1d-VAR
  INTEGER(KIND=INT4) :: Iteration_Idx                                         !iteration index
  INTEGER(KIND=INT4) :: Is_Fail                                               !failed retrieval flag
  INTEGER(KIND=INT4) :: Singular_Matrix_Flag                                  !singular matrix flag
  INTEGER(KIND=INT4) :: Diag_Output_Flag = 0     !sym%NO                      !diagnostic output flag
  REAL(KIND=REAL4) :: Conv_Test                                               !covergence metric
  REAL(KIND=REAL4) :: Conv_Criteria                                           !criterion for convergence
  REAL(KIND=REAL4), DIMENSION(NUM_OBS) :: Obs_Vector                          !vector of the observations
  REAL(KIND=REAL4), DIMENSION(NUM_OBS) :: Fwd_Model_Vector                    !forward model vector
  REAL(KIND=REAL4), DIMENSION(NUM_PARAM) :: Retv_Vector                       !vector of the retrieved parameters
  REAL(KIND=REAL4), DIMENSION(NUM_PARAM) :: Retv_Vector_Apriori               !vector of apriori estimates of x
  REAL(KIND=REAL4), DIMENSION(NUM_PARAM) :: Delta_Retv_Vector                 !change in x for next iteration
  REAL(KIND=REAL4), DIMENSION(NUM_OBS,NUM_PARAM) :: Kernel                    !Kernel or Jacobian Matrix
  REAL(KIND=REAL4), DIMENSION(NUM_PARAM,NUM_PARAM) :: Sa    !error covariance matrix of x_ap
  REAL(KIND=REAL4), DIMENSION(NUM_PARAM,NUM_PARAM) :: Inv_Sa!inverse of Sa
  REAL(KIND=REAL4), DIMENSION(NUM_PARAM,NUM_PARAM) :: Sx   !error covariance matrix of x
  REAL(KIND=REAL4), DIMENSION(NUM_PARAM,NUM_PARAM) :: Inv_Sx !inverse of Sy
  REAL(KIND=REAL4), DIMENSION(NUM_OBS,NUM_OBS) :: Sy        !error covariance matrix of y
  REAL(KIND=REAL4), DIMENSION(NUM_OBS,NUM_OBS) :: Inv_Sy    !inverse of y
  INTEGER(KIND=INT4) :: Param_Idx                                             !parameter index

  !--- local aliases for global variables
  INTEGER(KIND=INT4) :: Bad_Pixel_Mask_Chn14         !Chn 14 bad pixel flag
  INTEGER(KIND=INT4) :: Bad_Pixel_Mask_Chn15         !Chn 15 bad pixel flag
  INTEGER(KIND=INT4) :: Bad_Pixel_Mask_Chn16         !Chn 16 bad pixel flag
  LOGICAL :: Space_Mask                   !space view flag
  INTEGER(KIND=INT4) :: Cloud_Mask                   !cloud mask
  INTEGER(KIND=INT4) :: Cloud_Type                   !cloud type
  INTEGER(KIND=INT4) :: Cloud_Phase                  !cloud phase
  INTEGER(KIND=INT4) :: Sfc_Type                     !surface type
  INTEGER(KIND=INT4) :: Num_Elem                     !number of elements of segment
  INTEGER(KIND=INT4) :: Num_Line                     !number of lines in segment
  REAL(KIND=REAL4) :: BT_Chn14                       !11 micron brightness temperature
  REAL(KIND=REAL4) :: BT_Chn15                       !12 micron brightness temperature
  REAL(KIND=REAL4) :: BT_Chn16                       !13 micron brightness temperature
  REAL(KIND=REAL4) :: Cos_Sat_Zen                    !cosine of zenith angle
  REAL(KIND=REAL4) :: Temp_Tropo                     !temperature of the tropopause (K)

  !--- local variables used in this routine
  REAL(KIND=REAL4) :: Dummy_REAL4
  REAL(KIND=REAL4) :: Cld_Temp_Tmpy                  !current estimate of cloud temperature
  REAL(KIND=REAL4) :: Press_Tmpy                     !current estimate of cloud pressure
  REAL(KIND=REAL4) :: Cld_Hgt_Tmpy                   !current estimate of cloud height
  REAL(KIND=REAL4) :: Beta_11_13um_Tmpy              !current estimate of beta_11_13
  REAL(KIND=REAL4) :: Slope_Beta_11_13um_Beta_11_12um!slope of beta_11_13 wrt to beta_11_12 
  REAL(KIND=REAL4) :: A_Beta_Fit                     !constant of linear beta fit 
  REAL(KIND=REAL4) :: B_Beta_Fit                     !slope of linear beta fit 
  REAL(KIND=REAL4) :: Sfc_Temp                       !surface temperature
  REAL(KIND=REAL4) :: Sfc_Hgt                        !surface elevation
  REAL(KIND=REAL4) :: Sfc_Press                      !surface pressure
  REAL(KIND=REAL4) :: Temp_Lower_Bndy                !temperature of lower bnd
  REAL(KIND=REAL4) :: Prof_Weight                    !weight used in profile interpolation
  REAL(KIND=REAL4) :: Cld_OD_Apriori                        !apriori value of tau
  REAL(KIND=REAL4), DIMENSION(2):: Beta_2_Eff_Rad_Coef     !coefs for beta to re
  REAL(KIND=REAL4):: Delta_Cld_Temp_Sfc_Temp         !temperature diff cloud and surface
  REAL(KIND=REAL4):: Adj_Press                       !Pc adjusted for lower-level inver
  REAL(KIND=REAL4):: Adj_Sfc_Hgt                     !Zc adjusted for lower-level inver
  REAL(KIND=REAL4):: Lapse_Rate                      !lapse rate
  REAL(KIND=REAL4) :: R4_Dummy                       !REAL4 dummy variable
  INTEGER(KIND=INT4):: X_LRC_Idx_Temp                !temporary value of the lrc index in x
  INTEGER(KIND=INT4):: Y_LRC_Idx_Temp                !temporary value of the lrc index in y
  INTEGER(KIND=INT1), DIMENSION(:,:), ALLOCATABLE :: Processing_Order  !order to process each pixel
  INTEGER(KIND=INT1), DIMENSION(NUM_META_DATA) :: Meta_Data_Flags  !vector of individual meta data flags
  INTEGER(KIND=INT4), DIMENSION(NUM_META_DATA) :: Meta_Data_Flag_Bit_Depth 


  !--- rtm variables
  REAL(KIND=REAL4) :: Rad_Clr_TOA_11um              !clear sky 11 micron radiance at toa
  REAL(KIND=REAL4) :: Rad_Clr_TOA_12um              !clear sky 12 micron radiance at toa
  REAL(KIND=REAL4) :: Rad_Clr_TOA_13um              !clear sky 13 micron radiance at toa
  REAL(KIND=REAL4) :: Rad_Clear_11um                !clear 11 micron radiance toa
  REAL(KIND=REAL4) :: Rad_Clear_12um                !clear 12 micron radiance toa
  REAL(KIND=REAL4) :: Rad_Clear_13um                !clear 13 micron radiance toa
  REAL(KIND=REAL4) :: Atm_Trans_Abv_Cld_11um_RTM    !above cloud 11 micron transmission to toa
  REAL(KIND=REAL4) :: Atm_Trans_Abv_Cld_12um_RTM    !above cloud 12 micron transmission to toa
  REAL(KIND=REAL4) :: Atm_Trans_Abv_Cld_13um_RTM    !above cloud 13 micron transmission to toa
  REAL(KIND=REAL4) :: Rad_Abv_Cld_11um              !above cloud 11 micron radiance to toa
  REAL(KIND=REAL4) :: Rad_Abv_Cld_12um              !above cloud 12 micron radiance to toa
  REAL(KIND=REAL4) :: Rad_Abv_Cld_13um              !above cloud 13 micron radiance to toa
  REAL(KIND=REAL4) :: Slope_Planck_Emiss_11um       !slope of planck function at Tc for 11 microns
  REAL(KIND=REAL4) :: Slope_Planck_Emiss_12um       !slope of planck function at Tc for 12 microns
  REAL(KIND=REAL4) :: Slope_Planck_Emiss_13um       !slope of planck function at Tc for 13 microns
  REAL(KIND=REAL4) :: Emiss_Planck_Cld_Top_Chn14    !planck emission at Tc for 11 microns
  REAL(KIND=REAL4) :: Emiss_Planck_Cld_Top_Chn15    !planck emission at Tc for 12 microns
  REAL(KIND=REAL4) :: Emiss_Planck_Cld_Top_Chn16    !planck emission at Tc for 13 microns
  REAL(KIND=REAL4) :: Deriv_Emiss_12um_Deriv_Emiss_11um   !derivative of emiss12 wrt to emiss11
  REAL(KIND=REAL4) :: Deriv_Emiss_13um_Deriv_Emiss_11um   !derivative of emiss13 wrt to emiss11
  REAL(KIND=REAL4) :: Emiss_11um                     !11 micron cloud emissivity
  REAL(KIND=REAL4) :: Emiss_12um                     !12 micron cloud emissivity
  REAL(KIND=REAL4) :: Emiss_13um                     !13 micron cloud emissivity
  REAL(KIND=REAL4) :: Delta_Planck_Emiss_11um        !slope of planck emission at obs T for 11 microns
  REAL(KIND=REAL4) :: Delta_Planck_Emiss_12um        !slope of planck emission at obs T for 12 microns
  REAL(KIND=REAL4) :: Delta_Planck_Emiss_13um        !slope of planck emission at obs T for 13 microns
  REAL(KIND=REAL4) :: Rad_11um                       !observed 11 micron radiance
  REAL(KIND=REAL4) :: Rad_12um                       !observed 12 micron radiance
  REAL(KIND=REAL4) :: Rad_13um                       !observed 13 micron radiance

  !--- apriori variables
  REAL(KIND=REAL4):: Cld_Temp_Apriori_Uncer            !apriori uncertainty of Tc
  REAL(KIND=REAL4):: Cld_Emiss_11um_Apriori_Uncer      !apriori uncertainty of ec
  REAL(KIND=REAL4):: Beta_11_12um_Apriori_Uncer        !apriori uncertainty of beta
  REAL(KIND=REAL4):: Cld_Temp_Apriori_Opaque           !apriori Tc for opaque clouds
  REAL(KIND=REAL4):: Cld_Temp_Apriori_Cirrus           !apriori value of Tc for cirrus

  !--- forward model uncertainties
  REAL(KIND=REAL4):: BT_Clr_11um_Uncer             !clear-sky uncertainty of T11
  REAL(KIND=REAL4):: BT_Clr_11_12um_Uncer          !clear sky uncertainty of T11-T12
  REAL(KIND=REAL4):: BT_Clr_11_13um_Uncer          !clear sky uncertainty of T11-T13

  !--- local aliased POINTERs to geocat structures
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_Temp       !cloud temperature (K)
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_Emiss_11um !cloud emissivity at 11 microns
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_Press      !cloud pressure (hPa)
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_Hgt        !cloud height (m)
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_OD         !cloud visible optical depth
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Beta_11_12um   !beta value for 11 and 12 microns
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Emiss_Tropo_Chn14 !emiss11 at trop
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: LWP             !liquid water path
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: IWP             !ice water path
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_Eff_Radius  !cloud effective particle size (micron)
!  INTEGER(KIND=INT1), DIMENSION(:,:,:), POINTER:: QF          !retrieval quality flags
  INTEGER(KIND=INT1), DIMENSION(:,:), POINTER:: Cloud_Height_QF          !product quality flags
 !INTEGER(KIND=INT1), DIMENSION(:,:,:), POINTER :: Packed_Meta_Data_Flags     !2d array of packed meta data
  INTEGER(KIND=INT1), DIMENSION(:,:), POINTER :: Packed_Meta_Data_Flags     !2d array of packed meta data
  INTEGER(KIND=INT4), DIMENSION(:,:), POINTER:: X_LRC_Idx     !lrc index in x
  INTEGER(KIND=INT4), DIMENSION(:,:), POINTER:: Y_LRC_Idx     !lrc index in y
  INTEGER(KIND=INT1), DIMENSION(:,:), POINTER:: LRC_Mask      !mask for lrc
  INTEGER(KIND=INT1), DIMENSION(:,:), POINTER:: Cld_Layer     !cloud layer
  INTEGER(KIND=INT1), DIMENSION(:,:), POINTER:: Inver_Flag    !Inversion flag
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_Temp_Error   !cloud temperature error (K)
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_Emiss_11um_Error  !cloud emissivity error (K)
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_Beta_11_12um_Error  !cloud emissivity error (K)
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_Hgt_Error   !cloud height error (m)
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Cld_Press_Error !cloud pressure error (hPa)
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Tc_Lower_Cloud  !lower cloud layer temperature (K)
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Pc_Lower_Cloud  !lower cloud layer pressure (hPa)
  REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Zc_Lower_Cloud  !lower cloud layer height (m)

  !--- local POINTERs
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Rad_Prof_11um        !profile of clear-sky emission at 11 microns
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Rad_Prof_12um        !profile of clear-sky emission at 12 microns
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Rad_Prof_13um        !profile of clear-sky emission at 13 microns
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Atm_Trans_Prof_11um_RTM      !profile of clear-sky tranmission at 11 microns
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Atm_Trans_Prof_12um_RTM      !profile of clear-sky tranmission at 12 microns
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Atm_Trans_Prof_13um_RTM      !profile of clear-sky tranmission at 13 microns
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Cloud_Prof_11um    !profile of 11 um toa radiance for a black cloud at each level
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Cloud_Prof_12um    !profile of 12 um toa radiance for a black cloud at each level 
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Cloud_Prof_13um    !profile of 13 um toa radiance for a black cloud at each level 
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Temp_Prof_NWP      !temperature profile
  !REAL(KIND=REAL4), DIMENSION(:), POINTER :: Press_Prof_NWP     !pressure profile
  REAL(KIND=REAL4), DIMENSION(:), ALLOCATABLE :: Press_Prof_NWP     !pressure profile
  REAL(KIND=REAL4), DIMENSION(:), POINTER :: Hgt_Prof_NWP       !height profile
  INTEGER(KIND=INT1), DIMENSION(:), POINTER :: Inver_Level_Prof !inversion mask profile
  !INTEGER(KIND=INT4), DIMENSION(:), POINTER :: Inver_Level_Prof !inversion mask profile

  !----------------------------------------------------------------------
  ! Begin Executable Code
  !----------------------------------------------------------------------

    !--- set bit depth for meta data flags
    Meta_Data_Flag_Bit_Depth = 1

    !--- pass number of passes on if lrc's will be used
    IF (USE_LRC_FLAG == sym%NO) THEN
            Min_Val_Pass = 0
            Max_Val_Pass = 0
    ELSE
            Min_Val_Pass = 1
            Max_Val_Pass = 4
    END IF

    !--- store segment size in local variables
    Num_Line = Image%Number_of_Lines_Per_Segment     !sat%ny
    Num_Elem = Image%Number_of_Elements               !sat%nx
    allocate(Press_Prof_NWP(NLevels_Rtm))

    !-- ALLOCATE processing order array
    ALLOCATE(Processing_Order(Num_Elem,Num_Line),stat=Alloc_Status)
    Processing_Order = 0

    IF (Alloc_Status /= 0) THEN

      WRITE (Err_Message, *) &
             'Error allocating Processing_Order array'
      
      Error_Level = 2 ! AIT FATAL ERROR CODE
            
      call MESG("Baseline Cloud Height: "//trim(ERR_Message))

!     CALL Display_Message("Baseline Cloud Height", &
!                TRIM(Err_Message), &
!                Sym%FAILURE)
                 
        !  AIT Error Messaging
        !  CALL Error_Messaging (Routine_Name, Error_Message, Error_Level)       
      RETURN      

    ENDIF

  !--- point local POINTERs to geocat output structures
  Cld_Temp => ACHA%Tc                !out2(Algo_Num)%Cldt
  Cld_Press => ACHA%Pc               !out2(Algo_Num)%Cldp
  Cld_Hgt => ACHA%Zc                 !out2(Algo_Num)%Cldz
  Cld_Emiss_11um => ACHA%Ec          !out2(Algo_Num)%Cldemiss
  Cld_Eff_Radius=> ACHA%Reff         !out2(Algo_Num)%Cldreff
  Cld_OD => ACHA%Tau                 !out2(Algo_Num)%Cod_Vis
! LWP => out2(Algo_Num)%Cldlwp
! IWP => out2(Algo_Num)%Cldiwp
  Beta_11_12um => ACHA%Beta            !out2(Algo_Num)%Cldbeta1112
! QF => ACHA%Quality_Flag              !out2(Algo_Num)%Qcflg1
  Cloud_Height_QF => ACHA%Quality_Flag !out2(Algo_Num)%Cloud_Height_QF
  Emiss_Tropo_Chn14 => ch(31)%Emiss_Tropo  !out2(Algo_Num)%Emiss11_High
  X_LRC_Idx => i_LRC       !out2(Algo_Num)%X_LRC_Idx 
  Y_LRC_Idx => j_LRC       !out2(Algo_Num)%Y_LRC_Idx 
  LRC_Mask => Mask_LRC     !out2(Algo_Num)%LRC_Mask
  Cld_Layer => CCL%Cloud_Layer     !out2(Algo_Num)%Cld_Layer
  Inver_Flag => ACHA%Inversion_Flag                     !out2(Algo_Num)%Inver_Flag
  Cld_Temp_Error => ACHA%Tc_Uncertainty                 !out2(Algo_Num)%Tc_Error
  Cld_Emiss_11um_Error => ACHA%Ec_Uncertainty           !out2(Algo_Num)%Ec_Error
  Cld_Beta_11_12um_Error => ACHA%Beta_Uncertainty       !out2(Algo_Num)%Beta1112_Error
  Cld_Press_Error => ACHA%Pc_Uncertainty               !out2(Algo_Num)%Zc_Error
  Cld_Hgt_Error => ACHA%Zc_Uncertainty                 !out2(Algo_Num)%Zc_Error
  Tc_Lower_Cloud => ACHA%Lower_Tc                       !out2(Algo_Num)%Tc_Lower_Cloud
  Pc_Lower_Cloud => ACHA%Lower_Pc                       !out2(Algo_Num)%Pc_Lower_Cloud
  Zc_Lower_Cloud => ACHA%Lower_Zc                       !out2(Algo_Num)%Zc_Lower_Cloud
  Packed_Meta_Data_Flags => ACHA%Packed_Meta_Data_Flags !out2(Algo_Num)%Cloud_Height_Qpi

  !--- Initialize
  Cld_Temp = MISSING_VALUE_REAL4
  Cloud_Height_QF = INVALID_CTH_FAILED_RETRIEVAL
! QF = CTH_PARAM_FAILED_RETREVIAL
  Packed_Meta_Data_Flags = 0

  !----------------------------------------------------------------------
  ! compute spatial uniformity metrics to estimate expected accuracy
  ! of forward model
  !-----------------------------------------------------------------------
  CALL compute_spatial_uniformity(Elem_Width_3x3,  &
                                  Line_Width_3x3,  &
                                  Geo%Space_Mask, &   !sat%space_mask,  &
                                  ch(31)%Bt_Toa, &    !sat%bt14,  &
                                  BT_11um_Mean_3X3, &
                                  BT_11um_Max_3X3,  &
                                  BT_11um_Min_3X3,  &
                                  BT_11um_Std_3X3)

  CALL compute_spatial_uniformity(Elem_Width_3x3,  &
                                  Line_Width_3x3,  &
                                  Geo%Space_Mask, &                    !sat%space_mask,     &
                                  ch(31)%Bt_Toa - ch(32)%Bt_Toa, & !sat%bt14-sat%bt15,  &
                                  BT_1112um_Mean_3X3, &
                                  BT_1112um_Max_3X3,  &
                                  BT_1112um_Min_3X3,  &
                                  BT_1112um_Std_3X3)
  
  CALL compute_spatial_uniformity(Elem_Width_3x3,  &
                                  Line_Width_3x3,  &
                                  Geo%Space_Mask,  &                      !sat%space_mask,     &
                                  ch(31)%Bt_Toa - ch(33)%Bt_Toa, &   !sat%bt14-sat%bt16,  &
                                  BT_1113um_Mean_3X3, &
                                  BT_1113um_Max_3X3,  &
                                  BT_1113um_Min_3X3,  &
                                  BT_1113um_Std_3X3)

  !----------------------------------------------------------------------
  ! compute emissivity at tropopause 
  !----------------------------------------------------------------------
  CALL Compute_Emiss_Tropo_Chn14(Emiss_Tropo_Chn14,Num_Line)

  !----------------------------------------------------------------------
  ! compute local radiative center
  !----------------------------------------------------------------------
  LRC_Mask = sym%YES
  X_LRC_Idx = MISSING_VALUE_INT4
  Y_LRC_Idx = MISSING_VALUE_INT4

  !--- CALL routines to compute emiss lrc indices
  CALL gradient2d(Emiss_Tropo_Chn14, &
                  Image%Number_Of_Elements, &          !sat%nx, &
                  Image%Number_Of_Lines_Per_Segment, & !sat%ny, &
                  LRC_Mask, &
                  EMISS_TROPO_CHN14_GRADIENT_MIN, &
                  EMISS_TROPO_CHN14_GRADIENT_MAX, &
                  EMISS_TROPO_CHN14_GRADIENT_THRESH, &
                  X_LRC_Idx, &
                  Y_LRC_Idx)
  !-------------------------------------------------------------------------------
  ! Ensure that pixels with values above the threshold are
  ! treated as LRC's.  This not the default behavior
  !-------------------------------------------------------------------------------
  DO Line_Idx=1, Num_Line
     DO Elem_Idx=1, Num_Elem
         IF (Emiss_Tropo_Chn14(Elem_Idx,Line_Idx) >=EMISS_TROPO_CHN14_GRADIENT_THRESH ) THEN
          X_LRC_Idx(Elem_Idx,Line_Idx) = Elem_Idx
          Y_LRC_Idx(Elem_Idx,Line_Idx) = Line_Idx
         ENDIF
     ENDDO
  ENDDO

  !----------------------------------------------------------------------
  ! Set convergence criteria, should be << p (number of retrieved parameters)
  !----------------------------------------------------------------------
  Conv_Criteria = NUM_PARAM / 2.0

  !-----------------------------------------------------------
  ! ASSIGN pixel processing order
  !-----------------------------------------------------------
  IF (Use_LRC_Flag == sym%YES) THEN

    Line_Loop_2: DO Line_Idx = 1, Num_Line
     Element_Loop_2: DO Elem_Idx = 1, Num_Elem

        !---- local aliases for visual convenience
        Cloud_Type = Cld_Type(Elem_Idx,Line_Idx)       !sat%Cldtype(Elem_Idx,Line_Idx)    !cloud type

        !-- on pass 1, DO single layer lrc's
        IF ((Elem_Idx == X_LRC_Idx(Elem_Idx,Line_Idx)) .AND. &
            (Line_Idx == Y_LRC_Idx(Elem_Idx,Line_Idx)) .AND. &
            (Cloud_Type /= sym%OVERLAP_TYPE)) THEN
             Processing_Order(Elem_Idx,Line_Idx) = 1
        ENDIF

        !-- on pass 2, DO non-lrc water clouds
        IF (((Elem_Idx /= X_LRC_Idx(Elem_Idx,Line_Idx)) .OR. &
            (Line_Idx /= Y_LRC_Idx(Elem_Idx,Line_Idx))) .AND. &
            (Cloud_Type == sym%FOG_TYPE .OR. &
             Cloud_Type == sym%WATER_TYPE .OR. &
             Cloud_Type == sym%MIXED_TYPE .OR. &
             Cloud_Type == sym%SUPERCOOLED_TYPE)) THEN
            Processing_Order(Elem_Idx,Line_Idx) = 2
        ENDIF

        !-- on pass 3, DO lrc overlap clouds
        IF ((Elem_Idx == X_LRC_Idx(Elem_Idx,Line_Idx)) .AND. &
            (Line_Idx == Y_LRC_Idx(Elem_Idx,Line_Idx)) .AND. &
            (Cloud_Type == sym%OVERLAP_TYPE)) THEN
             Processing_Order(Elem_Idx,Line_Idx) = 3
        ENDIF

        !--  on pass-4 DO remaining
        IF (Processing_Order(Elem_Idx,Line_Idx) == 0) THEN
           Processing_Order(Elem_Idx,Line_Idx) = 4
        ENDIF

      END DO Element_Loop_2
    END DO Line_Loop_2

  ELSE
           Processing_Order = 0
  ENDIF


  !-----------------------------------------------------------------------
  ! Initialize lower cloud variables
  !-----------------------------------------------------------------------
  Tc_Lower_Cloud = MISSING_VALUE_REAL4
  Pc_Lower_Cloud = MISSING_VALUE_REAL4
  Zc_Lower_Cloud = MISSING_VALUE_REAL4

  !-----------------------------------------------------------------------
  ! perform multiple passes if using lrc results
  !-----------------------------------------------------------------------
  Pass_Loop: DO Pass_Idx = Min_Val_Pass, Max_Val_Pass
     
   !--------------------------------------------------------------------------
   !spatially interpolate water cloud temperature when appropriate
   !--------------------------------------------------------------------------
   IF ((Pass_Idx == 0) .OR. (Pass_Idx == 3)) THEN
        CALL Spatial_Interp_Lower_Cld_Pos(1,Num_Line, &
                   Cld_Press,Pc_Lower_Cloud,Tc_Lower_Cloud,Zc_Lower_Cloud)
   ENDIF

   !-----------------------------------------------------------------------
   ! Begin Loop over Each Pixel
   !----------------------------------------------------------------------- 
   Line_Loop_3: DO Line_Idx = 1, Num_Line
     Element_Loop_3: DO Elem_Idx = 1, Num_Elem

      !--- check if pixel is within sensor zenith angle limits
      !IF (sat%satzen(Elem_Idx,Line_Idx)  > SENSOR_ZEN_THRESH) THEN
      IF (Geo%SatZen(Elem_Idx,Line_Idx)  > SENSOR_ZEN_THRESH) THEN
            Cloud_Height_QF(Elem_Idx,Line_Idx) = INVALID_CTH_OUTSIDE_SEN_ZEN_RANGE
            CYCLE
      ENDIF
      
      !--- check if pixel should be processd in this path
      IF (Pass_Idx /= Processing_Order(Elem_Idx,Line_Idx)) THEN
            CYCLE
      ENDIF
    
      !--- initialize output for this pixel
      Cld_Temp(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_Press(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_Hgt(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_Emiss_11um(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_OD(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Beta_11_12um(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_Eff_Radius(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_Temp_Error(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_Emiss_11um_Error(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_Beta_11_12um_Error(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_Press_Error(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_Hgt_Error(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Cld_Layer(Elem_Idx,Line_Idx) = 0
      Inver_Flag(Elem_Idx,Line_Idx) = sym%NO

      !----------------------------------------------------------------------
      ! check for valid earth data point, if not, skip to next
      !----------------------------------------------------------------------
      Space_Mask = Geo%Space_Mask(Elem_Idx,Line_Idx)    !sat%space_mask(Elem_Idx,Line_Idx) !space view mask
      IF (Space_Mask) THEN
            Cloud_Height_QF(Elem_Idx,Line_Idx) = INVALID_CTH_DUE_TO_SPACE_PIX
            CYCLE
      END IF


      !---- local aliases for visual convenience
      Cloud_Type = Cld_Type(Elem_Idx,Line_Idx)            !sat%Cldtype(Elem_Idx,Line_Idx)    !cloud type 
      Cloud_Mask = CLDMASK%Cld_Mask(Elem_Idx,Line_Idx)    !sat%cldmask(Elem_Idx,Line_Idx)    !cloud mask
      Sfc_Type = Sfc%Sfc_Type(Elem_Idx,Line_Idx)          !sat%sfc_type(Elem_Idx,Line_Idx)     !surface classication
      BT_Chn14  = ch(31)%BT_Toa(Elem_Idx,Line_Idx)        !sat%bt14(Elem_Idx,Line_Idx)             !11.0 micron bt
      BT_Chn15  = ch(32)%Bt_Toa(Elem_Idx,Line_Idx)        !sat%bt15(Elem_Idx,Line_Idx)             !12.0 micron bt
      BT_Chn16  = ch(33)%Bt_Toa(Elem_Idx,Line_Idx)        !sat%bt16(Elem_Idx,Line_Idx)             !13.0 micron bt
      View_Zen_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)       !sat%ivza(Elem_Idx,Line_Idx)            !viewing zenith angle bin
      Cos_Sat_Zen = Geo%Coszen(Elem_Idx,Line_Idx)         !sat%cos_satzen(Elem_Idx,Line_Idx)    !cosine of viewing zenith angle 
      Bad_Pixel_Mask_Chn14 = Bad_Pixel_Mask(Elem_Idx,Line_Idx)  !sat%bad_pixel_mask(14,Elem_Idx,Line_Idx)  ! ch14 bad pixel
      Bad_Pixel_Mask_Chn15 = Bad_Pixel_Mask(Elem_Idx,Line_Idx)  !sat%bad_pixel_mask(15,Elem_Idx,Line_Idx)  ! ch15 bad pixel
      Bad_Pixel_Mask_Chn16 = Bad_Pixel_Mask(Elem_Idx,Line_Idx)  !sat%bad_pixel_mask(16,Elem_Idx,Line_Idx)  ! ch16 bad pixel
      Rad_Clr_TOA_11um = ch(31)%Rad_Toa_Clear(Elem_Idx,Line_Idx)!sat%rad_clr14(Elem_Idx,Line_Idx)  !clear 11 micron radiance
      Rad_Clr_TOA_12um = ch(32)%Rad_Toa_Clear(Elem_Idx,Line_Idx)!sat%rad_clr15(Elem_Idx,Line_Idx)  !clear 12 micron radiance
      Rad_Clr_TOA_13um = ch(33)%Rad_Toa_Clear(Elem_Idx,Line_Idx)!sat%rad_clr16(Elem_Idx,Line_Idx)  !clear 13 micron radiance
      X_Idx_NWP = NWP_PIX%I_NWP(Elem_Idx,Line_Idx)      !sat%x_nwp(Elem_Idx,Line_Idx)           !nwp longitude cell
      Y_Idx_NWP = NWP_PIX%J_NWP(Elem_Idx,Line_Idx)      !sat%y_nwp(Elem_Idx,Line_Idx)           !nwp latitude cell
      Temp_Tropo = NWP_PIX%Ttropo(Elem_Idx,Line_Idx)   !nwp%dat(X_Idx_NWP,Y_Idx_NWP)%ttropo      !tropopause temp
      Sfc_Idx_NWP = Rtm(X_Idx_NWP,Y_Idx_Nwp)%Sfc_Level  ! or Sfc_Level_Rtm_Pixel(Elem_Idx,Line_Idx)  !nwp%dat(X_Idx_NWP,Y_Idx_NWP)%sfc_level     !first nwp level above surface
      Tropo_Idx_NWP = Rtm(X_Idx_NWP,Y_Idx_NWP)%Tropo_Level  !nwp%dat(X_Idx_NWP,Y_Idx_NWP)%tropo_level !nwp level associated with tropopause
      Sfc_Hgt = Sfc%Zsfc(Elem_Idx,Line_Idx)             !sat%zsfc(Elem_Idx,Line_Idx)            !surface elevation
      Sfc_Press = NWP_PIX%Psfc(Elem_Idx,Line_Idx)       !nwp%dat(X_Idx_NWP,Y_Idx_NWP)%psfc     !surface pressure
      Sfc_Temp = NWP_PIX%Tsfc(Elem_Idx,Line_Idx)        !nwp%dat(X_Idx_NWP,Y_Idx_NWP)%tsfc      !surface temp
      Temp_Prof_NWP => Rtm(X_Idx_NWP,Y_Idx_NWP)%T_Prof(Tropo_Idx_NWP:Sfc_Idx_NWP)  !nwp%dat(X_Idx_NWP,Y_Idx_NWP)%tlev(Tropo_Idx_NWP:Sfc_Idx_NWP)  ! temperature profile
      Press_Prof_NWP = P_Std_Rtm(Tropo_Idx_NWP:Sfc_Idx_NWP)                        !=> nwp%dat(X_Idx_NWP,Y_Idx_NWP)%plev(Tropo_Idx_NWP:Sfc_Idx_NWP)  ! pressure profile
      Hgt_Prof_NWP => RTM(X_Idx_Nwp,Y_Idx_NWP)%Z_Prof(Tropo_Idx_NWP:Sfc_Idx_NWP)   !nwp%dat(X_Idx_NWP,Y_Idx_NWP)%zlev(Tropo_Idx_NWP:Sfc_Idx_NWP)  ! height profile

      Inver_Level_Prof => RTM(X_Idx_Nwp,Y_Idx_Nwp)%Inver_Prof(Tropo_Idx_NWP:Sfc_Idx_NWP)   !nwp%dat(X_Idx_NWP,Y_Idx_NWP)%inversion_lev(Tropo_Idx_NWP:Sfc_Idx_NWP) !inver prof

      Rad_Prof_11um => RTM(X_Idx_Nwp,Y_Idx_NWP)%d(View_Zen_Idx)%ch(31)%Rad_Atm_Profile(Tropo_Idx_NWP:Sfc_Idx_NWP) !rtm(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%rad_atm_clr14(Tropo_Idx_NWP:Sfc_Idx_NWP)   !11 um radiance profile
      Rad_Prof_12um => RTM(X_Idx_Nwp,Y_Idx_NWP)%d(View_Zen_Idx)%ch(32)%Rad_Atm_Profile(Tropo_Idx_NWP:Sfc_Idx_NWP) !rtm(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%rad_atm_clr15(Tropo_Idx_NWP:Sfc_Idx_NWP)   !12 um radiance profile
      Rad_Prof_13um => RTM(X_Idx_Nwp,Y_Idx_NWP)%d(View_Zen_Idx)%ch(33)%Rad_Atm_Profile(Tropo_Idx_NWP:Sfc_Idx_NWP) !rtm(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%rad_atm_clr16(Tropo_Idx_NWP:Sfc_Idx_NWP)   !13 um radiance profile
      Atm_Trans_Prof_11um_RTM => RTM(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%ch(31)%Trans_Atm_Profile(Tropo_Idx_NWP:Sfc_Idx_NWP) !rtm(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%trans_atm_clr14(Tropo_Idx_NWP:Sfc_Idx_NWP) !11 um transmission  profile
      Atm_Trans_Prof_12um_RTM => RTM(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%ch(32)%Trans_Atm_Profile(Tropo_Idx_NWP:Sfc_Idx_NWP) !rtm(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%trans_atm_clr15(Tropo_Idx_NWP:Sfc_Idx_NWP) !12 um transmission profile
      Atm_Trans_Prof_13um_RTM => RTM(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%ch(33)%Trans_Atm_Profile(Tropo_Idx_NWP:Sfc_Idx_NWP) !rtm(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%trans_atm_clr16(Tropo_Idx_NWP:Sfc_Idx_NWP) !13 um transmission profile
      Cloud_Prof_11um => RTM(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%ch(31)%Rad_BB_Cloud_Profile(Tropo_Idx_NWP:Sfc_Idx_NWP) !rtm(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%cloud_prof14(Tropo_Idx_NWP:Sfc_Idx_NWP) !bb 11 um cld radiance profile
      Cloud_Prof_12um => RTM(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%ch(32)%Rad_BB_Cloud_Profile(Tropo_Idx_NWP:Sfc_Idx_NWP) !rtm(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%cloud_prof15(Tropo_Idx_NWP:Sfc_Idx_NWP) !bb 12 um cld radianec profile
      Cloud_Prof_13um => RTM(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%ch(33)%Rad_BB_Cloud_Profile(Tropo_Idx_NWP:Sfc_Idx_NWP) !rtm(X_Idx_NWP,Y_Idx_NWP)%d(View_Zen_Idx)%cloud_prof16(Tropo_Idx_NWP:Sfc_Idx_NWP) !bb 13 um cld radianec profile

      !-----------------------------------------------------------------------
      !-- compute number of profile levels of our solution space
      !-----------------------------------------------------------------------
      Num_Level_NWP_Profile = (Sfc_Idx_NWP-Tropo_Idx_NWP)+1                 !number of levels between sfc and trop

      !-----------------------------------------------------------------------
      ! determine level to stop for low-level inversions
      !-----------------------------------------------------------------------
            
      CALL Prof_Lookup_Using_P(Hgt_Prof_NWP, Press_Prof_NWP, Temp_Prof_NWP, R4_Dummy, MIN_P_INVERSION,  &
         R4_Dummy,Min_NWP_Level_Inver, R4_Dummy)

      CALL Prof_Lookup_Using_P(Hgt_Prof_NWP, Press_Prof_NWP, Temp_Prof_NWP, R4_Dummy,  &
         Sfc_Press - DELTA_PSFC_INVERSION,                    &
         R4_Dummy,Max_NWP_Level_Inver, R4_Dummy)

      !--- constrain to be in bounds
      Max_NWP_Level_Inver = max(1,min(Num_Level_NWP_Profile,Max_NWP_Level_Inver))
      Min_NWP_Level_Inver = max(1,min(Max_NWP_Level_Inver,Min_NWP_Level_Inver))

      !----------------------------------------------------------------------
      ! check for valid data in 11, 12 and 13 micron data - if not, skip to next
      !----------------------------------------------------------------------
      IF ((Bad_Pixel_Mask_Chn14 == sym%YES) .OR. &
          (Bad_Pixel_Mask_Chn15 == sym%YES) .OR. &
          (Bad_Pixel_Mask_Chn16 == sym%YES))THEN

          Cloud_Height_QF(Elem_Idx,Line_Idx) = INVALID_CTH_BAD_DATA

          CYCLE
      END IF

      !----------------------------------------------------------------------
      ! if not cloudy, skip to next pixel
      !----------------------------------------------------------------------
      IF ((Cloud_Mask == sym%CLEAR) .OR. &
          (Cloud_Mask == sym%PROB_CLEAR)) THEN
          
          Cloud_Height_QF(Elem_Idx,Line_Idx) = INVALID_CTH_CLD_MASK_CLR
          
          CYCLE
      END IF

      !----------------------------------------------------------------------
      ! if unknown cloud type or phase, cycle
      !----------------------------------------------------------------------
      IF (Cloud_Type == sym%UNKNOWN_TYPE) THEN
          
          Cloud_Height_QF(Elem_Idx,Line_Idx) = INVALID_CTH_INVALID_CLD_TYPE
          
          CYCLE          
          
      END IF

      !----------------------------------------------------------------------
      ! populate observation vector, y
      !----------------------------------------------------------------------
      Obs_Vector(1) = BT_Chn14
      Obs_Vector(2) = BT_Chn14 - BT_Chn15
      Obs_Vector(3) = BT_Chn14 - BT_Chn16


      !----------------------------------------------------------------------
      !------------ measurement and forward model error
      !----------------------------------------------------------------------
      IF (Sfc_Type == sym%WATER_SFC) THEN
         BT_Clr_11um_Uncer = T11_CLR_UNCER_WATER
         BT_Clr_11_12um_Uncer = T11_12_CLR_UNCER_WATER
         BT_Clr_11_13um_Uncer = T11_13_CLR_UNCER_WATER
      ELSE
         BT_Clr_11um_Uncer = T11_CLR_UNCER_LAND
         BT_Clr_11_12um_Uncer = T11_12_CLR_UNCER_LAND
         BT_Clr_11_13um_Uncer = T11_13_CLR_UNCER_LAND
      END IF

      !----------------------------------------------------------------------
      ! Construct a priori estimates and their  error matrix, Sa
      ! based on cloud type
      !----------------------------------------------------------------------

      !--- base apriori Tc for cirrus from tropopause temp
      Cld_Temp_Apriori_Opaque = BT_Chn14
      Cld_Temp_Apriori_Cirrus = Temp_Tropo + TC_AP_TROPO_OFFSET_CIRRUS

      !--- set apriori values of ec based on loud type

      !- fog
      IF (Cloud_Type == sym%FOG_TYPE) THEN          
          Retv_Vector_Apriori(1) = Cld_Temp_Apriori_Opaque
          Cld_OD_Apriori = TAU_AP_FOG_TYPE
          Cld_Temp_Apriori_Uncer = TC_AP_UNCER_OPAQUE
          Cld_Emiss_11um_Apriori_Uncer = 2.0*EC_AP_UNCER_OPAQUE

      !- water clouds
      ELSEIF (Cloud_Type == sym%WATER_TYPE) THEN
          Retv_Vector_Apriori(1) = Cld_Temp_Apriori_Opaque
          Cld_OD_Apriori = TAU_AP_WATER_TYPE  
          Cld_Temp_Apriori_Uncer =TC_AP_UNCER_OPAQUE
          Cld_Emiss_11um_Apriori_Uncer = EC_AP_UNCER_OPAQUE

      !- mixed-phase clouds
      ELSEIF (Cloud_Type == sym%MIXED_TYPE .OR. &
                Cloud_Type == sym%SUPERCOOLED_TYPE) THEN
          Retv_Vector_Apriori(1) = Cld_Temp_Apriori_Opaque
          Cld_OD_Apriori = TAU_AP_MIXED_TYPE
          Cld_Temp_Apriori_Uncer = TC_AP_UNCER_OPAQUE
          Cld_Emiss_11um_Apriori_Uncer = EC_AP_UNCER_OPAQUE

      !- thick-ice clouds
      ELSEIF (Cloud_Type == sym%TICE_TYPE) THEN
          Retv_Vector_Apriori(1) = Cld_Temp_Apriori_Opaque
          Cld_OD_Apriori = TAU_AP_OPAQUE_ICE_TYPE
          Cld_Temp_Apriori_Uncer = TC_AP_UNCER_OPAQUE
          Cld_Emiss_11um_Apriori_Uncer = EC_AP_UNCER_OPAQUE

      !- cirrus clouds
      ELSEIF (Cloud_Type == sym%CIRRUS_TYPE) THEN
          Retv_Vector_Apriori(1) = Cld_Temp_Apriori_Cirrus
          Cld_OD_Apriori = TAU_AP_CIRRUS_TYPE
          Cld_Temp_Apriori_Uncer = TC_AP_UNCER_CIRRUS
          Cld_Emiss_11um_Apriori_Uncer = EC_AP_UNCER_CIRRUS

      !- multilayer clouds
      ELSEIF (Cloud_Type == sym%OVERLAP_TYPE) THEN
          Retv_Vector_Apriori(1) = Cld_Temp_Apriori_Cirrus
          Cld_OD_Apriori = TAU_AP_OVERLAP_TYPE
          Cld_Temp_Apriori_Uncer = TC_AP_UNCER_CIRRUS
          Cld_Emiss_11um_Apriori_Uncer = EC_AP_UNCER_CIRRUS
      END IF

      !---- determine beta apriori and fit parameters based on phase (derived from type)

       !--- water phase clouds
       IF ((Cloud_Type == sym%FOG_TYPE) .OR. &
           (Cloud_Type == sym%WATER_TYPE) .OR. &
           (Cloud_Type == sym%SUPERCOOLED_TYPE) .OR. &
           (Cloud_Type == sym%MIXED_TYPE)) THEN
           Retv_Vector_Apriori(3) = BETA_AP_WATER
           Beta_11_12um_Apriori_Uncer = BETA_AP_UNCER_WATER
           A_Beta_Fit = A_BETA_FIT_WATER
           B_Beta_Fit = B_BETA_FIT_WATER
           Beta_2_Eff_Rad_Coef = BETA2RE_COEF_WATER
           Cloud_Phase = sym%WATER_PHASE

       !--- ice phase clouds
       ELSE
           Retv_Vector_Apriori(3) = BETA_AP_ICE
           Beta_11_12um_Apriori_Uncer = BETA_AP_UNCER_ICE
           A_Beta_Fit = A_BETA_FIT_ICE
           B_Beta_Fit = B_BETA_FIT_ICE
           Beta_2_Eff_Rad_Coef = BETA2RE_COEF_ICE
           Cloud_Phase = sym%ICE_PHASE
       END IF

       !--- convert tau_ap to a nadir emissivity
       Retv_Vector_Apriori(2) = 1.0 - exp(-Cld_OD_Apriori/Cos_Sat_Zen)  !slow!

       !--- constrain apriori value of Tc
       Retv_Vector_Apriori(1) = max(MIN_ALLOWABLE_TC, min(Temp_Prof_NWP(Num_Level_NWP_Profile),Retv_Vector_Apriori(1)))

       !--- compute Sa and its inverse
       Sa = 0.0
       Sa(1,1) = Cld_Temp_Apriori_Uncer
       Sa(2,2) = Cld_Emiss_11um_Apriori_Uncer
       Sa(3,3) = Beta_11_12um_Apriori_Uncer

       !--- modify Tc apriori based on lrc information --- moved, AER Bug fix. WCS3
       X_LRC_Idx_Temp = X_LRC_Idx(Elem_Idx,Line_Idx)
       Y_LRC_Idx_Temp = Y_LRC_Idx(Elem_Idx,Line_Idx)

       IF ((X_LRC_Idx_Temp /= MISSING_VALUE_INT4) .and.   &
            (Y_LRC_Idx_Temp /= MISSING_VALUE_INT4)) THEN
          IF ((Cld_Temp(X_LRC_Idx_Temp,Y_LRC_Idx_Temp) /=  &
                MISSING_VALUE_REAL4) .and.                     &
              (Cld_Emiss_11um(X_LRC_Idx_Temp,Y_LRC_Idx_Temp) /=  &
                MISSING_VALUE_REAL4))  THEN

                Retv_Vector_Apriori(1) =  Cld_Temp(X_LRC_Idx_Temp,Y_LRC_Idx_Temp)

                Sa(1,1) = 5.0 +  &
                     (1.0-Cld_Emiss_11um(X_LRC_Idx_Temp,Y_LRC_Idx_Temp)) *Cld_Temp_Apriori_Uncer

          END IF
        END IF



       Sa = Sa**2
       CALL INVERT_3x3(Sa,Inv_Sa,Singular_Matrix_Flag)
       IF (Singular_Matrix_Flag == sym%YES) THEN
          Is_Fail = sym%YES
          CYCLE
       END IF

        !--- write diagnostic information to screen - remove from final version
        IF (Diag_Output_Flag == sym%YES) THEN
          print *, "======================================================================="
          print *, "new retrieval ", Elem_Idx,Line_Idx,BT_Chn14,BT_Chn15,BT_Chn16,Sfc_Type,&
            Cos_Sat_Zen,Cloud_Type,BT_11um_Std_3X3(Elem_Idx,Line_Idx),BT_1112um_Std_3X3(Elem_Idx,Line_Idx)
          print *, "nwp cell = ", X_Idx_NWP, Y_Idx_NWP
          print *, "y = ", Obs_Vector
          print *, "x_ap = ", Retv_Vector_Apriori
          print *, "Sa = ", Sa
          print *, "Sa_inv = ", Inv_Sa
        END IF

       !---------------------------------------------------------------------------
       !--- modify clear radiances to simulate that from an opaque cloud at
       !--- the predetermined lower cloud temperature 
       !--- when a multi-layer situation is suspected
       !---------------------------------------------------------------------------
       IF ((Cloud_Type == sym%OVERLAP_TYPE) .and.  &
           (Pc_Lower_Cloud(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4)) THEN

          CALL Prof_Lookup_Using_P(Hgt_Prof_NWP,                      &
                                   Press_Prof_NWP,                    &
                                   Temp_Prof_NWP,                     &
                                   Dummy_REAL4,      &
                                   Pc_Lower_Cloud(Elem_Idx,Line_Idx), &
                                   Dummy_REAL4,      &
                                   NWP_Level_Idx,                     &
                                   Prof_Weight)

          Temp_Lower_Bndy = Tc_Lower_Cloud(Elem_Idx,Line_Idx)

          !--- constrain NWP level
          NWP_Level_Idx = min(Num_Level_NWP_Profile-1,max(NWP_Level_Idx,1))

          !--- interpolate
          Rad_Clear_11um = Cloud_Prof_11um(NWP_Level_Idx) +  &
                         Prof_Weight*(Cloud_Prof_11um(NWP_Level_Idx+1) - Cloud_Prof_11um(NWP_Level_Idx))

          Rad_Clear_12um = Cloud_Prof_12um(NWP_Level_Idx) +  &
                         Prof_Weight*(Cloud_Prof_12um(NWP_Level_Idx+1) - Cloud_Prof_12um(NWP_Level_Idx))

          Rad_Clear_13um = Cloud_Prof_13um(NWP_Level_Idx) +  &
                         Prof_Weight*(Cloud_Prof_13um(NWP_Level_Idx+1) - Cloud_Prof_13um(NWP_Level_Idx))

        ELSE

          Temp_Lower_Bndy = Sfc_Temp

          Rad_Clear_11um = Rad_Clr_TOA_11um

          Rad_Clear_12um = Rad_Clr_TOA_12um

          Rad_Clear_13um = Rad_Clr_TOA_13um

        END IF


        !----------------------------------------------------------------------
        ! BEGIN RETRIEVAL ITERATION LOOP
        !----------------------------------------------------------------------

        !--- initialize counter and convergence parameter
        Iteration_Idx = 0
        Conv_Test = huge(Conv_Test)
        Is_Fail = sym%NO

        !--- set initial value of x to a priori
        Retv_Vector = Retv_Vector_Apriori

        !--- begin loop
        Retrieval_Loop: do

          !--- update iteration counter
          Iteration_Idx = Iteration_Idx + 1

          IF (Diag_Output_Flag == sym%YES) THEN
            print *, "testing convergence for iter", Iteration_Idx, Obs_Vector-Fwd_Model_Vector,Conv_Test, Conv_Criteria
          END IF

          !--- test for convergence
          IF (Conv_Test < Conv_Criteria) THEN
            IF (Diag_Output_Flag == sym%YES) THEN
              print *,"converged ", Retv_Vector, Iteration_Idx, Conv_Test,Conv_Criteria
            END IF
            Is_Fail = sym%NO
            EXIT
          END IF

          !--- test for exceeding maximum number of iterations
          IF (Iteration_Idx > ITER_MAX) THEN
            IF (Diag_Output_Flag == sym%YES) THEN
              print *,"failed convergence"
            END IF
            Is_Fail = sym%YES
            EXIT
          END IF

          !----------------------------------------------------------------------
          ! FORWARD MODEL SECTION
          !---------------------------------------------------------------------- 

          !--- based on current estimate of cloud temperature, compute level of cloud top

          Cld_Temp_Tmpy = Retv_Vector(1)

          CALL Prof_Lookup_Using_T(Hgt_Prof_NWP,    &
                                   Press_Prof_NWP,  &
                                   Temp_Prof_NWP,   &
                                   Cld_Hgt_Tmpy,    &
                                   Press_Tmpy,      &
                                   Cld_Temp_Tmpy,   &
                                   NWP_Level_Idx,   &
                                   Prof_Weight)

  
          !--- constrain NWP level
          NWP_Level_Idx = min(Num_Level_NWP_Profile-1,max(NWP_Level_Idx,1))

          !--- interpolate above cloud radiometric parameters
          Atm_Trans_Abv_Cld_11um_RTM = Atm_Trans_Prof_11um_RTM(NWP_Level_Idx) +  &
                                       Prof_Weight*(Atm_Trans_Prof_11um_RTM(NWP_Level_Idx+1) -  &
                                       Atm_Trans_Prof_11um_RTM(NWP_Level_Idx))

          Atm_Trans_Abv_Cld_12um_RTM = Atm_Trans_Prof_12um_RTM(NWP_Level_Idx) + &
                                       Prof_Weight*(Atm_Trans_Prof_12um_RTM(NWP_Level_Idx+1) -  &
                                       Atm_Trans_Prof_12um_RTM(NWP_Level_Idx))

          Atm_Trans_Abv_Cld_13um_RTM = Atm_Trans_Prof_13um_RTM(NWP_Level_Idx) + &
                                       Prof_Weight*(Atm_Trans_Prof_13um_RTM(NWP_Level_Idx+1) -  &
                                       Atm_Trans_Prof_13um_RTM(NWP_Level_Idx))

          Rad_Abv_Cld_11um = Rad_Prof_11um(NWP_Level_Idx) +  &
                             Prof_Weight*(Rad_Prof_11um(NWP_Level_Idx+1) - Rad_Prof_11um(NWP_Level_Idx))

          Rad_Abv_Cld_12um = Rad_Prof_12um(NWP_Level_Idx) +  &
                             Prof_Weight*(Rad_Prof_12um(NWP_Level_Idx+1) - Rad_Prof_12um(NWP_Level_Idx))

          Rad_Abv_Cld_13um = Rad_Prof_13um(NWP_Level_Idx) +  &
                             Prof_Weight*(Rad_Prof_13um(NWP_Level_Idx+1) - Rad_Prof_13um(NWP_Level_Idx))
  
          !--- compute planck emission at cloud temperature
          Emiss_Planck_Cld_Top_Chn14 = planck_rad_fast(31,Cld_Temp_Tmpy,dB_dT = Slope_Planck_Emiss_11um)
          Emiss_Planck_Cld_Top_Chn15 = planck_rad_fast(32,Cld_Temp_Tmpy,dB_dT = Slope_Planck_Emiss_12um)
          Emiss_Planck_Cld_Top_Chn16 = planck_rad_fast(33,Cld_Temp_Tmpy,dB_dT = Slope_Planck_Emiss_13um)

          !---- compute channel emissivities for forward model use

          !-- 11 micron
          Emiss_11um = min(Retv_Vector(2),0.999999)

          !-- 12 micron
          Deriv_Emiss_12um_Deriv_Emiss_11um = Retv_Vector(3) * (1.0-Emiss_11um)**(Retv_Vector(3)-1.0)
          Emiss_12um = 1.0 - (1.0-Emiss_11um)**Retv_Vector(3)

          !-- 13 micron
          Beta_11_13um_Tmpy = A_Beta_Fit + B_Beta_Fit * Retv_Vector(3)
          Slope_Beta_11_13um_Beta_11_12um = B_Beta_Fit
          Deriv_Emiss_13um_Deriv_Emiss_11um = Beta_11_13um_Tmpy * (1.0-Emiss_11um)**(Beta_11_13um_Tmpy - 1.0)
          Emiss_13um = 1.0 - (1.0-Emiss_11um)**Beta_11_13um_Tmpy


          !--- CALL forward model for 11 microns
          Rad_11um = Emiss_11um*Rad_Abv_Cld_11um +  &
                     Atm_Trans_Abv_Cld_11um_RTM * Emiss_11um * Emiss_Planck_Cld_Top_Chn14 + &
                     (1.0 - Emiss_11um) * Rad_Clear_11um

          Fwd_Model_Vector(1) = planck_temp_fast(31,Rad_11um,dB_dT = Delta_Planck_Emiss_11um)

          IF (Diag_Output_Flag == 1) THEN
            IF (Fwd_Model_Vector(1) < 0.0) THEN
              print *, "negative f(1) "
              print *, Emiss_11um, Rad_Abv_Cld_11um, Atm_Trans_Abv_Cld_11um_RTM, Emiss_Planck_Cld_Top_Chn14, Rad_Clear_11um
            END IF
          END IF

          !--- CALL forward model for 12 microns
          Rad_12um = Emiss_12um*Rad_Abv_Cld_12um +  &
                     Atm_Trans_Abv_Cld_12um_RTM *  &
                     Emiss_12um * Emiss_Planck_Cld_Top_Chn15 + &
                     (1.0 - Emiss_12um) * Rad_Clear_12um

           Fwd_Model_Vector(2) = Fwd_Model_Vector(1) -  &
                              planck_temp_fast(32,Rad_12um,dB_dT = Delta_Planck_Emiss_12um)

           !--- CALL forward model for 13 microns
           Rad_13um = Emiss_13um*Rad_Abv_Cld_13um +  &
                      Atm_Trans_Abv_Cld_13um_RTM *   &
                      Emiss_13um * Emiss_Planck_Cld_Top_Chn16 + &
                      (1.0 - Emiss_13um) * Rad_Clear_13um

           Fwd_Model_Vector(3) = Fwd_Model_Vector(1) -  &
                              planck_temp_fast(33,Rad_13um,dB_dT = Delta_Planck_Emiss_13um)

           !--- compute Kernel Matrix
           Kernel(1,1) = (Atm_Trans_Abv_Cld_11um_RTM * Emiss_11um *  &
                         Slope_Planck_Emiss_11um) / Delta_Planck_Emiss_11um    !dT_11 / dT_c

           Kernel(1,2) = (Rad_Abv_Cld_11um +  &
                         Atm_Trans_Abv_Cld_11um_RTM*Emiss_Planck_Cld_Top_Chn14 - Rad_Clear_11um)/ &
                         Delta_Planck_Emiss_11um                               !dT_11 / demiss_c

           Kernel(1,3) = 0.0                                               !dT_11 / dbeta_11_12

           Kernel(2,1) = Kernel(1,1) - Atm_Trans_Abv_Cld_12um_RTM *  &
                         Emiss_12um * Slope_Planck_Emiss_12um /  &
                         Delta_Planck_Emiss_12um                       !d(T_11 - T_12) / dT_c

           Kernel(2,2) = Kernel(1,2) - &
                        (Rad_Abv_Cld_12um+Atm_Trans_Abv_Cld_12um_RTM*Emiss_Planck_Cld_Top_Chn15-Rad_Clear_12um)*&
                        (Deriv_Emiss_12um_Deriv_Emiss_11um)/ &
                        Delta_Planck_Emiss_12um                                !d(T_11 - T_12) / demiss_4 

           Kernel(2,3) = (Rad_Abv_Cld_12um+Atm_Trans_Abv_Cld_12um_RTM*Emiss_Planck_Cld_Top_Chn15-Rad_Clear_12um)/&
                         Delta_Planck_Emiss_12um  * &
                         alog(1.0-Emiss_11um)*(1.0-Emiss_12um)                  !d(T_11 - T_12) / dbeta_11_12

           Kernel(3,1) = Kernel(1,1) - &
                         Atm_Trans_Abv_Cld_13um_RTM * Emiss_13um * Slope_Planck_Emiss_13um / &
                         Delta_Planck_Emiss_13um                                !d(T_11 - T_13) / dT_c

           Kernel(3,2) = Kernel(1,2) -   &
                        (Rad_Abv_Cld_13um+Atm_Trans_Abv_Cld_13um_RTM*Emiss_Planck_Cld_Top_Chn16-Rad_Clear_13um)*&
                        (Deriv_Emiss_13um_Deriv_Emiss_11um)/ &
                        Delta_Planck_Emiss_13um                                !d(T_11 - T_13) / demiss_4 

           Kernel(3,3) = (Rad_Abv_Cld_13um+Atm_Trans_Abv_Cld_13um_RTM*Emiss_Planck_Cld_Top_Chn16-Rad_Clear_13um)/ &
                         Delta_Planck_Emiss_13um * &
                         alog(1.0-Emiss_11um)*(1.0-Emiss_13um)                  !d(T_11 - T_13) / dbeta_11_13

           Kernel(3,3) = Kernel(3,3) * Slope_Beta_11_13um_Beta_11_12um     !d(T_11 - T_13) / dbeta_11_12

           IF (Diag_Output_Flag == sym%YES) THEN
             print *, "clear ch4 terms = ", Rad_Clear_11um, Atm_Trans_Abv_Cld_11um_RTM, Rad_Abv_Cld_11um
             print *, "clear ch5 terms = ", Rad_Clear_12um, Atm_Trans_Abv_Cld_12um_RTM, Rad_Abv_Cld_12um
             print *, "cloud radiances = ", Cld_Temp_Tmpy,Emiss_Planck_Cld_Top_Chn14,Emiss_Planck_Cld_Top_Chn15
             print *, "cloud emiss  = ", Emiss_11um, Emiss_12um
             print *, "x = ", Retv_Vector
             print *, "f = ", Fwd_Model_Vector
             print *, "y-f = ", Obs_Vector-Fwd_Model_Vector
             print *, "Kernel(1,1) = ", Kernel(1,1)
             print *, "Kernel(1,2) = ", Kernel(1,2)
             print *, "Kernel(1,3) = ", Kernel(1,3)
             print *, "Kernel(2,1) = ", Kernel(2,1)
             print *, "Kernel(2,2) = ", Kernel(2,2)
             print *, "Kernel(2,3) = ", Kernel(2,3)
             print *, "Kernel(3,1) = ", Kernel(3,1)
             print *, "Kernel(3,2) = ", Kernel(3,2)
             print *, "Kernel(3,3) = ", Kernel(3,3)
           END IF

            !--- compute Sy  - this is a function of x
            Sy = 0.0

            Sy(1,1) = T11_CAL_UNCER**2    +  &
                                      ((1.0-Retv_Vector(2))*BT_Clr_11um_Uncer)**2 +  &
                                      BT_11um_Std_3X3(Elem_Idx,Line_Idx)**2

            Sy(2,2) = T11_12_CAL_UNCER**2 + &
                                      ((1.0-Retv_Vector(2))*BT_Clr_11_12um_Uncer)**2 +  &
                                      BT_1112um_Std_3X3(Elem_Idx,Line_Idx)**2

            Sy(3,3) = T11_13_CAL_UNCER**2 + &
                                      ((1.0-Retv_Vector(2))*BT_Clr_11_13um_Uncer)**2 +  &
                                      BT_1113um_Std_3X3(Elem_Idx,Line_Idx)**2


            CALL INVERT_3x3(Sy,Inv_Sy,Singular_Matrix_Flag)
            IF (Singular_Matrix_Flag == sym%YES) THEN
             Is_Fail = sym%YES
             EXIT
            END IF

           !--- compute convariance matrix of retrieved parameter, Sy
           Inv_Sx = Inv_Sa +  &
                                       matmul(transpose(Kernel), &
                                       matmul(Inv_Sy,Kernel)) !(Eq.102 Rodgers)

           CALL INVERT_3x3(Inv_Sx, &
                           Sx, &
                           Singular_Matrix_Flag)

           IF (Singular_Matrix_Flag == sym%YES) THEN
             Is_Fail = sym%YES
             EXIT
           END IF

           !--- compute increments to retrieved parameters, delta_x
           Delta_Retv_Vector = matmul(Sx, &
                                  (matmul(transpose(Kernel), &
                                   matmul(Inv_Sy,(Obs_Vector-Fwd_Model_Vector))) +  &
                                   matmul(Inv_Sa,Retv_Vector_Apriori-Retv_Vector) ))

           !--- compute the convergence metric
           Conv_Test = ABS(sum(Delta_Retv_Vector*matmul(Inv_Sx,Delta_Retv_Vector)))

           !---constrain the retrieved parameter increment
           DO Param_Idx = 1,NUM_PARAM
             IF (ABS(Delta_Retv_Vector(Param_Idx)) > ABS(Retv_Vector(Param_Idx)) / 2.0) THEN
               Delta_Retv_Vector(Param_Idx) = SIGN( Retv_Vector(Param_Idx) / 2.0 ,&
                                                    Delta_Retv_Vector(Param_Idx) )
               
             END IF
           END DO

           IF (ABS(Delta_Retv_Vector(1)) > 20.0) THEN
             Delta_Retv_Vector(1) = SIGN( 20.0 , Delta_Retv_Vector(1) )
           END IF

           IF (ABS(Delta_Retv_Vector(2)) > 0.2) THEN
             Delta_Retv_Vector(2) = SIGN( 0.2 , Delta_Retv_Vector(2) )
           END IF

           IF (ABS(Delta_Retv_Vector(3)) > 0.1) THEN
             Delta_Retv_Vector(3) = SIGN( 0.1 , Delta_Retv_Vector(3) )
           END IF
  
           IF (Diag_Output_Flag == sym%YES) THEN
             print *, "delta_x = ", Delta_Retv_Vector
           END IF

           !--- update the retrieved parameter vector
           Retv_Vector = Retv_Vector + Delta_Retv_Vector

           !--- constrain the retrieved parameter vector
           Retv_Vector(1) = max(180.0,min(Sfc_Temp,Retv_Vector(1)))     !should we DO this?
           Retv_Vector(2) = max(0.0,min(Retv_Vector(2),1.0))
           Retv_Vector(3) = max(0.8,min(Retv_Vector(3),2.0))

         !-----------------------------------------------------------------------
         ! END OF ITERATION LOOP
         !-----------------------------------------------------------------------
         END DO Retrieval_Loop

      !----------------------------------------------------------------------
      ! if retrieval is successful, populate output with results
      !----------------------------------------------------------------------
         IF (Is_Fail == sym%NO) THEN

                Cld_Temp(Elem_Idx,Line_Idx) = Retv_Vector(1)
                Beta_11_12um(Elem_Idx,Line_Idx) = Retv_Vector(3)

                !--- derive an optical depth
                IF (Retv_Vector(2) < 1.00) THEN
          
                    !this is the split-window optical depth adjusted to nadir
                    Cld_OD(Elem_Idx,Line_Idx) = -Cos_Sat_Zen*alog(1.0 - Retv_Vector(2)) 
            
                    Cld_Emiss_11um(Elem_Idx,Line_Idx) = 1.0 - &
                                                        exp(-Cld_OD(Elem_Idx,Line_Idx))
            
                    !this is now a "visible" optical depth (assuming 2.0 scale factor)
                    Cld_OD(Elem_Idx,Line_Idx) = 2.0 * Cld_OD(Elem_Idx,Line_Idx)    
            
                ELSE
          
                    Cld_OD(Elem_Idx,Line_Idx) = 20.0
                    Cld_Emiss_11um(Elem_Idx,Line_Idx) = 1.0
            
                END IF

          !--- derive a particle size
          IF (Beta_11_12um(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) THEN
           IF (Beta_11_12um(Elem_Idx,Line_Idx) > 1.0) THEN
            Cld_Eff_Radius(Elem_Idx,Line_Idx) = 1.0 / &
                    (Beta_2_Eff_Rad_Coef(1) + Beta_2_Eff_Rad_Coef(2)*Beta_11_12um(Elem_Idx,Line_Idx))
           ELSE
            Cld_Eff_Radius(Elem_Idx,Line_Idx) = 100.0
           END IF
          END IF

          !------------------- determine liquid and ice water paths
          !--- lwp assumes Qext=2.1, 0.833 adjust for constant re
          !--- iwp relation from Andrew Heymsfeld JAM 2003 - Fig 7
          !-------------------------------------------------------
          IF ((Cld_OD(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) .and. &
              (Cld_Eff_Radius(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4)) THEN

                    IF (Cloud_Phase == sym%WATER_PHASE) THEN
            
                        !LWP(Elem_Idx,Line_Idx) = max(0.0, &
                        !                            0.833 * 0.65 * &
                        !                            Cld_OD(Elem_Idx,Line_Idx) * &
                        !                            Cld_Eff_Radius(Elem_Idx,Line_Idx))
                    
                    ELSE
                        !IWP(Elem_Idx,Line_Idx) = &
                        !        max(0.0,(Cld_OD(Elem_Idx,Line_Idx)**(1.0/0.84)) / 0.065)
                    END IF

                END IF

                !--- quality flags of the retrieved parameters
            !   DO Param_Idx = 1,NUM_PARAM    !loop over parameters
            !       IF (Sx(Param_Idx,Param_Idx) < &
            !           0.111*Sa(Param_Idx,Param_Idx) ) THEN
            !  
            !           QF(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_1_3_APRIORI_RETREVIAL
            !           
            !     
            !       ELSEIF (Sx(Param_Idx,Param_Idx) < &
            !               0.444*Sa(Param_Idx,Param_Idx)) THEN
            !  
            !           QF(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_2_3_APRIORI_RETREVIAL
            !       ELSE
            !           QF(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_LOW_QUALITY_RETREVIAL
            !       END IF
          !
          !     END DO
          !     
                !WCS3
                ! If the solution converged, no matter what the reason, it is
                ! valid
                
          !     IF (QF(1,Elem_Idx,Line_Idx) /= CTH_PARAM_FAILED_RETREVIAL) THEN
          !         Cloud_Height_QF(Elem_Idx,Line_Idx) = VALID_CTH_RETRIEVAL
          !     ENDIF
                
                !--- store error estimates
                Cld_Temp_Error(Elem_Idx,Line_Idx) = sqrt(Sx(1,1))
                Cld_Emiss_11um_Error(Elem_Idx,Line_Idx) = sqrt(Sx(2,2))
                Cld_Beta_11_12um_Error(Elem_Idx,Line_Idx) = sqrt(Sx(3,3))

                !--- diagnostic output if selected (remove in final version)
                IF (Diag_Output_Flag == sym%YES) THEN
                    print  *, "failed flag = ", Is_Fail
                    print  *, "final answer = ", Retv_Vector
                    print  *, "number of iterations = ", Iteration_Idx
                END IF

      !--- if not successful, set to missing
      ELSE

                Cld_Temp(Elem_Idx,Line_Idx) = Missing_Value_Real4
                Cld_Emiss_11um(Elem_Idx,Line_Idx) = Missing_Value_Real4
                Beta_11_12um(Elem_Idx,Line_Idx) = Missing_Value_Real4
                Cld_OD(Elem_Idx,Line_Idx) = Missing_Value_Real4

      END IF

      !---------------------------------------------------------------------
      ! Estimate pressure and height errors
      !---------------------------------------------------------------------
            IF (Cld_Temp(Elem_Idx,Line_Idx) > 0.0) THEN

                !--- find level in profiles corresponding to temperature solution
                
                CALL Prof_Lookup_Using_T(Hgt_Prof_NWP,  &
                                         Press_Prof_NWP,&
                                         Temp_Prof_NWP, &
                                         Cld_Hgt(Elem_Idx,Line_Idx), &
                                         Cld_Press(Elem_Idx,Line_Idx), &
                                         Cld_Temp(Elem_Idx,Line_Idx),& 
                                         NWP_Level_Idx, &
                                         Prof_Weight)

                !--- constrain NWP level
                NWP_Level_Idx = min(Num_Level_NWP_Profile-1,max(NWP_Level_Idx,1))

                !--- constrain results that fall below the surface to reside at the surface
                ! For the temporary profiles, level 1 is always the tropopause
                ! Num_Level_NWP_Profile is always the lowest level. 
                ! Num_Level_NWP_Profile-1  is the lowest level returnable by Prof_Lookup_Using_T
                IF ((Cld_Temp(Elem_Idx,Line_Idx) > Temp_Prof_NWP(Num_Level_NWP_Profile)) .AND. &
                    (NWP_Level_Idx == (Num_Level_NWP_Profile - 1))) THEN
 
                    Cld_Press(Elem_Idx,Line_Idx) = Sfc_Press
                    Cld_Hgt(Elem_Idx,Line_Idx) =  Sfc_Hgt

                ENDIF
         
          
                !--- constrain results that fall above the Tropopause
                ! note, Prof_Lookup_Using_T does return a height value using
                ! an assumed stratospheric lapse rate
                IF ((Cld_Temp(Elem_Idx,Line_Idx) < Temp_Prof_NWP(1)) .AND. &
                    (NWP_Level_Idx == (1))) THEN

                    Cld_Press(Elem_Idx,Line_Idx) = Press_Prof_NWP(1)
      
                ENDIF
                
                !--- constrain cloud height to be above surface
                Cld_Hgt(Elem_Idx,Line_Idx) = max(0.0,Cld_Hgt(Elem_Idx,Line_Idx))

                !--- constrain cloud pressure to be above surface
                Cld_Press(Elem_Idx,Line_Idx) = &
                                    min(Sfc_Press,Cld_Press(Elem_Idx,Line_Idx))
                              

                !--- height error
                Lapse_Rate =  (Temp_Prof_NWP(NWP_Level_Idx+1) - &
                               Temp_Prof_NWP(NWP_Level_Idx)) / &
                               (Hgt_Prof_NWP(NWP_Level_Idx+1) - &
                                Hgt_Prof_NWP(NWP_Level_Idx))
                        
                IF (Lapse_Rate /= 0.0) THEN
                    Cld_Hgt_Error(Elem_Idx,Line_Idx) = &
                            Cld_Temp_Error(Elem_Idx,Line_Idx) / ABS(Lapse_Rate)
                ENDIF

                !--- pressure error
                Lapse_Rate =  (Temp_Prof_NWP(NWP_Level_Idx+1) - &
                               Temp_Prof_NWP(NWP_Level_Idx)) / &
                               (Press_Prof_NWP(NWP_Level_Idx+1) - &
                               Press_Prof_NWP(NWP_Level_Idx))
                        
                IF (Lapse_Rate /= 0.0) THEN
                    Cld_Press_Error(Elem_Idx,Line_Idx) = &
                                    Cld_Temp_Error(Elem_Idx,Line_Idx) / &
                                    ABS(Lapse_Rate)
                ENDIF

          !-----------------------------------------------------------------------
          ! if a valid cloud temperature was retrieved, estimate height and pressure
          !-----------------------------------------------------------------------

                !--- do bottom-up over ocean for non-ice clouds
                IF (USE_SFC_INVER_FLAG == sym%YES) THEN

                    IF ((Cloud_Type == sym%FOG_TYPE) .OR.         &
                        (Cloud_Type == sym%WATER_TYPE) .OR.       &
                        (Cloud_Type == sym%SUPERCOOLED_TYPE) .OR. &
                        (Cloud_Type == sym%MIXED_TYPE)) THEN

                        !--- check for water surface
                        IF (Sfc_Type == sym%WATER_SFC) THEN

                            !--- check for low-level inversion
                            IF (maxval( &
                   Inver_Level_Prof(Min_NWP_Level_Inver:Max_NWP_Level_Inver)) > 0) THEN


                      !--- set flag to indicate that inversion logic was followed
                                Inver_Flag(Elem_Idx,Line_Idx) = sym%YES

                      !--- compute temp diff cloud and surface
                               Delta_Cld_Temp_Sfc_Temp = Sfc_Temp - &
                                                         Cld_Temp(Elem_Idx,Line_Idx)

                          !--- estimate cloud height based on assumed lapse rate
                                Adj_Sfc_Hgt = Delta_Cld_Temp_Sfc_Temp / &
                                              LAPSE_RATE_OCEAN + Sfc_Hgt

                                !--- recompute Pc knowing Zc (but not Tc)
                                CALL Prof_Lookup_Using_Z(Hgt_Prof_NWP,  &
                                                         Press_Prof_NWP,&
                                                         Temp_Prof_NWP, &
                                                         Adj_Sfc_Hgt,   &
                                                         Adj_Press,     &
                                                         R4_Dummy,      &
                                                         NWP_Level_Idx, &
                                                         R4_Dummy)

                     !--- take values that are closer to the surface
                                IF (Adj_Press > Cld_Press(Elem_Idx,Line_Idx)) THEN
                                    Cld_Press(Elem_Idx,Line_Idx) = Adj_Press
                                      Cld_Hgt(Elem_Idx,Line_Idx) = Adj_Sfc_Hgt
                                ENDIF
                                
                            ENDIF  !inversion level check

                        ENDIF   !water surface check

                    ENDIF    !ice cloud check

                ENDIF    !check of USE_SFC_INVER_FLAG

          !--- compute cloud layer (0=clear;1=low;2=mid;3=high)
             IF (Cld_Press(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) THEN

                    IF (Cld_Press(Elem_Idx,Line_Idx) < PC_HIGH_CLOUD_MAX) THEN
                        Cld_Layer(Elem_Idx,Line_Idx) = 3
                    ELSEIF (Cld_Press(Elem_Idx,Line_Idx) < PC_MID_CLOUD_MAX) THEN
                        Cld_Layer(Elem_Idx,Line_Idx) = 2
                    ELSE
                        Cld_Layer(Elem_Idx,Line_Idx) = 1
                    END IF
               END IF

            END IF    !check of valid retrieval
            
            !--- nullify local 1d POINTERs
            Temp_Prof_NWP => null()
            deallocate (Press_Prof_Nwp)    !Press_Prof_NWP => null()
            Hgt_Prof_NWP => null()
            Inver_Level_Prof => null()
            Rad_Prof_11um => null()
            Rad_Prof_12um => null()
            Rad_Prof_13um => null()
            Atm_Trans_Prof_11um_RTM => null()
            Atm_Trans_Prof_12um_RTM =>  null()
            Atm_Trans_Prof_13um_RTM =>  null()
            Cloud_Prof_11um => null()
            Cloud_Prof_12um => null()
            Cloud_Prof_13um => null()

      !----------------------------------------------------------------------
      ! Populate Meta Data for this pixel
      ! 1 - Cloud Height Attempted (0 = no / 1 = yes)
      ! 2 - Bias Correction Employed (0 = no / 1 = yes)
      ! 3 - Ice Cloud Retrieval (0 = no / 1 = yes)
      ! 4 - Local Radiatve Center Processing Used (0 = no / 1 = yes)
      ! 5 - Multi-layer Retrieval (0 = no / 1 = yes)
      ! 6 - Lower Cloud Interpolation Used (0 = no / 1 = yes)
      ! 7 - Boundary Layer Inversion Assumed  (0 = no / 1 = yes)
      !----------------------------------------------------------------------
            Meta_Data_Flags(1) = sym%YES

            Meta_Data_Flags(2) = USE_BIAS_CORRECTION_FLAG

            IF (Cloud_Phase == sym%ICE_PHASE) THEN
                Meta_Data_Flags(3) = sym%YES
            ELSE
                Meta_Data_Flags(3) = sym%NO
            ENDIF

            Meta_Data_Flags(4) = USE_LRC_FLAG

            IF (Cloud_Type == sym%OVERLAP_TYPE) THEN
                Meta_Data_Flags(5) = sym%YES
            ELSE
                Meta_Data_Flags(5) = sym%NO
            ENDIF

            Meta_Data_Flags(6) = USE_LOWER_CLD_INTERP_FLAG

            Meta_Data_Flags(7) = USE_SFC_INVER_FLAG

            !--- pack the meta data flags
!           CALL PACK_BYTES(Meta_Data_Flags(1:NUM_META_DATA), &
!                           Meta_Data_Flag_Bit_Depth(1:NUM_META_DATA), &
!                           Packed_Meta_Data_Flags(Elem_Idx,Line_Idx))
                            !Packed_Meta_Data_Flags(:,Elem_Idx,Line_Idx))

   !----------------------------------------------------------------------
   ! END OF PIXEL LOOP
   !----------------------------------------------------------------------

            END DO Element_Loop_3

        END DO Line_Loop_3

    END DO Pass_Loop

!----------------------------------------------------------------------
! DEALLOCATE memory
!----------------------------------------------------------------------
    CALL Destroy_Spatial_Uniformity(BT_11um_Mean_3X3,  &
                                  BT_11um_Max_3X3,   &
                                  BT_11um_Min_3X3,   &
                                  BT_11um_Std_3X3)

    CALL Destroy_Spatial_Uniformity(BT_1112um_Mean_3X3, &
                                  BT_1112um_Max_3X3,  &
                                  BT_1112um_Min_3X3,  &
                                  BT_1112um_Std_3X3)

    CALL Destroy_Spatial_Uniformity(BT_1113um_Mean_3X3, &
                                  BT_1113um_Max_3X3,  &
                                  BT_1113um_Min_3X3,  &
                                  BT_1113um_Std_3X3)

    !----------------------------------------------------------------------
    ! nullify local pointers
    !----------------------------------------------------------------------
    Cld_Temp => null()
    Cld_Press => null()
    Cld_Hgt => null()
    Cld_OD => null()
    Beta_11_12um => null()
    Cld_Emiss_11um => null()
    Cld_Eff_Radius => null()
    LWP => null()
    IWP => null()
    !QF => null()
    Cloud_Height_QF => null()
    Emiss_Tropo_Chn14 => null()
    X_LRC_Idx => null()
    Y_LRC_Idx => null()
    LRC_Mask => null()
    Cld_Layer => null()
    Inver_Flag => null()
    Cld_Temp_Error => null()
    Cld_Emiss_11um_Error => null()
    Cld_Beta_11_12um_Error => null()
    Cld_Press_Error => null()
    Cld_Hgt_Error => null()
    Tc_Lower_Cloud => null()
    Pc_Lower_Cloud => null()
    Zc_Lower_Cloud => null()
    Packed_Meta_Data_Flags => null()

    !-- DEALLOCATE processing order array
    DEALLOCATE(Processing_Order,stat=Alloc_Status)

    IF (Alloc_Status /= 0) THEN

      WRITE (Err_Message, *) &
             'Error deallocating Processing_Order array'

      Error_Level = 2 ! AIT FATAL ERROR CODE
            
      !CALL Display_Message("Baseline Cloud Height", &
      !          TRIM(Err_Message), &
      !          Sym%FAILURE)
      call MESG("Baseline Cloud Height: "//trim(ERR_Message))
                 
        !  AIT Error Messaging
        !  CALL Error_Messaging (Routine_Name, Error_Message, Error_Level)       
      RETURN      

    ENDIF
  
  
!----------------------------------------------------------------------
! END OF THIS ROUTINE
!----------------------------------------------------------------------
  END SUBROUTINE Baseline_Cloud_Height_main

!====================================================================
! SUBROUTINE Name: Spatial_Interp_Lower_Cld_Pos
!
! Function:
!   Routine to spatially interpret water cloud temperature values to surrounding
!   pixels
!
! Description:
!   This SUBROUTINE spatially interprets water cloud temperatures to the
!   surrounding pixels
!
! Calling Sequence:
!          CALL Spatially_Interp_Lower_Cld_Pos(1,Num_Line, &
!                   Cld_Press,Pc_Lower_Cloud,Tc_Lower_Cloud,Zc_Lower_Cloud)
!
! Inputs: Line_Idx_Start
!         Num_Line
!         Cld_Press
!         Pc_Lower_Cloud
!         Tc_Lower_Cloud
!         Zc_Lower_Cloud 
!
! Outputs: Pc_Lower_Cloud
!          Tc_Lower_Cloud
!          Zc_Lower_Cloud
!
! Dependencies: None
!
! Restrictions:  None
!
!====================================================================
!-------------------------------------------------------------------------------
! Routine to spatially interpret water cloud temperature values to surrounding
! pixels
!
! input:  Use_Lower_Cld_Interp_Flag = 0 (DO no interp, assume Zc=2km) /= 0 (DO spatial interp)
!-------------------------------------------------------------------------------
  SUBROUTINE Spatial_Interp_Lower_Cld_Pos(Line_Idx_Start,Line_Idx_End, &
                                          Cld_Press,Pc_Lower_Cloud, &
                                          Tc_Lower_Cloud,Zc_Lower_Cloud)

      INTEGER, INTENT(IN):: Line_Idx_Start
      INTEGER, INTENT(IN):: Line_Idx_End
      REAL, DIMENSION(:,:), INTENT(IN):: Cld_Press
      REAL, DIMENSION(:,:), INTENT(inout):: Pc_Lower_Cloud
      REAL, DIMENSION(:,:), INTENT(inout):: Tc_Lower_Cloud
      REAL, DIMENSION(:,:), INTENT(inout):: Zc_Lower_Cloud
      INTEGER:: Line_Idx
      INTEGER:: Elem_Idx
      INTEGER:: X_Idx_NWP
      INTEGER:: Y_Idx_NWP
      INTEGER:: Tropo_Idx_NWP
      INTEGER:: Sfc_Idx_NWP
      INTEGER:: NWP_Level_Idx
      INTEGER:: Num_Elem
      INTEGER:: Num_Line
      INTEGER:: Elem_Idx_Box_Radius
      INTEGER:: Line_Idx_Box_Radius
      INTEGER:: Elem_Idx_Box_Right
      INTEGER:: Elem_Idx_Box_Left
      INTEGER:: Line_Idx_Box_Top
      INTEGER:: Line_Idx_Box_Bot
      INTEGER:: Count_Valid
      INTEGER, DIMENSION(:,:), ALLOCATABLE:: Valid_Lower_Cld_Mask
      REAL:: R4_Dummy

      Num_Line = Image%Number_Of_Lines_Per_Segment    !sat%ny
      Num_Elem = Image%Number_Of_Elements             !sat%nx
      
      !----------------------------------------------------------------------
      ! First, use default value when it is available
      !----------------------------------------------------------------------
       Line_Loop_1: DO Line_Idx = Line_Idx_Start, Line_Idx_End
         Element_Loop_1: DO Elem_Idx = 1, Num_Elem
         
           !IF ((sat%bad_pixel_mask(14,Elem_Idx,Line_Idx) == sym%YES) .OR. &
           !     (sat%space_mask(Elem_Idx, Line_Idx) == sym%SPACE))THEN 
           IF ((Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) .OR. &
                (Geo%Space_Mask(Elem_Idx, Line_Idx) ))THEN 
                 Pc_Lower_Cloud(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
                 CYCLE
           ENDIF


          X_Idx_NWP = NWP_PIX%I_Nwp(Elem_Idx,Line_Idx)    !sat%x_nwp(Elem_Idx,Line_Idx)
          Y_Idx_NWP = NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)    !sat%y_nwp(Elem_Idx,Line_Idx)
          Pc_Lower_Cloud(Elem_Idx,Line_Idx) =  NWP_PIX%Psfc(Elem_Idx,Line_Idx) - Pc_Lower_Cloud_offset
          !Pc_Lower_Cloud(Elem_Idx,Line_Idx) = nwp%dat(X_Idx_NWP,Y_Idx_NWP)%psfc - Pc_Lower_Cloud_offset

        ENDDO Element_Loop_1
       ENDDO Line_Loop_1
       
       !----------------------------------------------------------------------
       ! Perform  spatial interpolation from neighboring pixels
       !----------------------------------------------------------------------
       IF (Use_Lower_Cld_Interp_Flag == sym%YES) THEN

        !--- set box width
        Elem_Idx_Box_Radius = INTERP_LOWER_CLOUD_PIXEL_RADIUS
        Line_Idx_Box_Radius = INTERP_LOWER_CLOUD_PIXEL_RADIUS

        allocate(Valid_Lower_Cld_Mask(Num_Elem,Num_Line))
        
         Valid_Lower_Cld_Mask = 0

         !WHERE(((sat%Cldtype == sym%FOG_TYPE) .OR. &
         !    (sat%Cldtype == sym%WATER_TYPE) .OR. &
         !    (sat%Cldtype == sym%MIXED_TYPE) .OR. &
         !    (sat%Cldtype == sym%SUPERCOOLED_TYPE)) .AND. &
         WHERE(((Cld_Type == sym%FOG_TYPE) .OR. &
             (Cld_Type == sym%WATER_TYPE) .OR. &
             (Cld_Type == sym%MIXED_TYPE) .OR. &
             (Cld_Type == sym%SUPERCOOLED_TYPE)) .AND. &
             Cld_Press /= MISSING_VALUE_REAL4)
             Valid_Lower_Cld_Mask = 1
         ENDWHERE
         
         Line_Loop_2: DO Line_Idx = Line_Idx_Start, Line_Idx_End
             Element_Loop_2: DO Elem_Idx = 1, Num_Elem

             !IF (sat%Cldtype(Elem_Idx,Line_Idx)  /= sym%OVERLAP_TYPE) THEN
             IF (Cld_Type(Elem_Idx,Line_Idx)  /= sym%OVERLAP_TYPE) THEN
                 CYCLE
             ENDIF

             Elem_Idx_Box_Right = MIN(Num_Elem,max(1,Elem_Idx-Elem_Idx_Box_Radius))
             Elem_Idx_Box_Left = MIN(Num_Elem,max(1,Elem_Idx+Elem_Idx_Box_Radius))
             Line_Idx_Box_Top = MIN(Num_Line,max(1,Line_Idx-Line_Idx_Box_Radius))
             Line_Idx_Box_Bot = MIN(Num_Line,max(1,Line_Idx+Line_Idx_Box_Radius))

             Count_Valid = sum(Valid_Lower_Cld_Mask(Elem_Idx_Box_Right:Elem_Idx_Box_Left, &
                                                    Line_Idx_Box_Top:Line_Idx_Box_Bot))

              IF (Count_Valid > 0) THEN
                  Pc_Lower_Cloud(Elem_Idx,Line_Idx) =  &
                               SUM(Cld_Press(Elem_Idx_Box_Right:Elem_Idx_Box_Left, &
                                             Line_Idx_Box_Top:Line_Idx_Box_Bot)* &
                                   Valid_Lower_Cld_Mask(Elem_Idx_Box_Right:Elem_Idx_Box_Left, &
                                                        Line_Idx_Box_Top:Line_Idx_Box_Bot)) /  &
                               Count_Valid
              ENDIF
             

             ENDDO Element_Loop_2
          ENDDO Line_Loop_2
          
       deallocate(Valid_Lower_Cld_Mask)

       ENDIF   !end of Use_Lower_Cld_Interp_Flag check

       !----------------------------------------------------------------
       !  Compute Height and Temperature
       !----------------------------------------------------------------
       Line_Loop_3: DO Line_Idx = Line_Idx_Start, Line_Idx_End
         Element_Loop_3: DO Elem_Idx = 1, Num_Elem

            !-- if a bad pixel, set to missing
            !IF (sat%bad_pixel_mask(14,Elem_Idx,Line_Idx) == sym%YES) THEN
            IF (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) THEN
               Pc_Lower_Cloud(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
               Tc_Lower_Cloud(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
               Zc_Lower_Cloud(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
               CYCLE
            ENDIF

            !--- if not overlap, set to all missing
            !IF (sat%Cldtype(Elem_Idx,Line_Idx) /= sym%OVERLAP_TYPE .OR. &
            IF (Cld_Type(Elem_Idx,Line_Idx) /= sym%OVERLAP_TYPE .OR. &
                Pc_Lower_Cloud(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) THEN
                Pc_Lower_Cloud(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
                Zc_Lower_Cloud(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
                Tc_Lower_Cloud(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
                CYCLE
            ENDIF

            !--- compute T and Z from P
            X_Idx_NWP = NWP_PIX%I_NWP(Elem_Idx,Line_Idx)    !sat%x_nwp(Elem_Idx,Line_Idx)
            Y_Idx_NWP = NWP_PIX%J_NWP(Elem_Idx,Line_Idx)    !sat%y_nwp(Elem_Idx,Line_Idx)
            Tropo_Idx_NWP =  Rtm(X_Idx_NWP,Y_Idx_NWP)%Tropo_Level    !nwp%dat(X_Idx_NWP,Y_Idx_NWP)%tropo_level
            Sfc_Idx_NWP =  Rtm(X_Idx_NWP,Y_Idx_NWP)%Sfc_Level        !nwp%dat(X_Idx_NWP,Y_Idx_NWP)%sfc_level

            !--- recompute T and Z knowing P
!           CALL prof_lookup_using_p(nwp%dat(X_Idx_NWP,Y_Idx_NWP)%zlev(Tropo_Idx_NWP:Sfc_Idx_NWP), &
!                                    nwp%dat(X_Idx_NWP,Y_Idx_NWP)%plev(Tropo_Idx_NWP:Sfc_Idx_NWP), &
!                                    nwp%dat(X_Idx_NWP,Y_Idx_NWP)%tlev(Tropo_Idx_NWP:Sfc_Idx_NWP), &
            CALL prof_lookup_using_p(RTM(X_Idx_Nwp,Y_Idx_NWP)%Z_Prof(Tropo_Idx_NWP:Sfc_Idx_NWP), &
                                     P_Std_RTM(Tropo_Idx_NWP:Sfc_Idx_NWP), &
                                     RTM(X_Idx_Nwp,Y_Idx_NWP)%T_Prof(Tropo_Idx_NWP:Sfc_Idx_NWP), &
                                     Zc_Lower_Cloud(Elem_Idx,Line_Idx),      &
                                     Pc_Lower_Cloud(Elem_Idx,Line_Idx),      &
                                     Tc_Lower_Cloud(Elem_Idx,Line_Idx),      &
                                     NWP_Level_Idx, &
                                     R4_Dummy)

         ENDDO Element_Loop_3
        ENDDO Line_Loop_3
        

  END SUBROUTINE Spatial_Interp_Lower_Cld_Pos


 !=======================================================================
 ! Compute 11 micron emissivity at Tropopause
 !=======================================================================
 SUBROUTINE Compute_Emiss_Tropo_Chn14(Emiss_Tropo_Chn14,Number_of_Lines_in_this_Segment)
   INTEGER(KIND=INT4),INTENT(IN) :: Number_of_Lines_in_this_Segment
   REAL(KIND=REAL4), DIMENSION(:,:), INTENT(OUT):: Emiss_Tropo_Chn14
   INTEGER(KIND=INT1):: Tropo_Idx_NWP
   INTEGER(KIND=INT1):: View_Zen_Idx
   INTEGER:: X_NWP_Idx
   INTEGER:: Y_NWP_Idx
   INTEGER:: Elem_Idx
   INTEGER:: Line_Idx
   REAL(KIND=REAL4) :: Rad_Chn14
   REAL(KIND=REAL4) :: Clr_Rad_Chn14
   REAL(KIND=REAL4) :: Blkbdy_Tropo_Rad_Chn14

   !--- initialize
   Emiss_Tropo_Chn14 = Missing_Value_Real4

    Line_Loop: DO Line_Idx=1, Number_of_Lines_in_this_Segment
      !Element_Loop: DO Elem_Idx = 1, sat%nx
      Element_Loop: DO Elem_Idx = 1, Image%Number_Of_Elements

       !IF (sat%space_mask(Elem_Idx,Line_Idx) == sym%NO) THEN
       IF ( .NOT. Geo%Space_Mask(Elem_Idx,Line_Idx)  .and. Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%NO) THEN

            !
            !---nwp longitude cell
            !
            X_NWP_Idx =          NWP_PIX%I_NWP(Elem_Idx,Line_Idx)   !sat%x_nwp(Elem_Idx,Line_Idx)

            !
            !---nwp latitude cell
            !
            Y_NWP_Idx =          NWP_PIX%J_NWP(Elem_Idx,Line_Idx)   !sat%y_nwp(Elem_Idx,Line_Idx)

            !
            !---nwp level associated with tropopause
            !
            Tropo_Idx_NWP =      RTM(X_NWP_Idx,Y_NWP_Idx)%Tropo_Level !nwp%dat(X_NWP_Idx,Y_NWP_Idx)%Tropo_Level

            !
            !---viewing zenith angle bin
            !
            View_Zen_Idx =          Zen_Idx_RTM(Elem_Idx,Line_Idx)    !sat%ivza(Elem_Idx,Line_Idx)

            !
            !---11 um radiance
            !
            Rad_Chn14  =        ch(31)%Rad_Toa(Elem_Idx,Line_Idx)     !sat%rad14(Elem_Idx,Line_Idx)

            !
            !---clear 11 micron radiance
            !
            Clr_Rad_Chn14 =     ch(31)%Rad_Toa_Clear(Elem_Idx,Line_Idx)    !sat%rad_clr14(Elem_Idx,Line_Idx)

            !
            !---BB 11 um rad at tropopause
            !
            !Blkbdy_Tropo_Rad_Chn14 = rtm(X_NWP_Idx,Y_NWP_Idx)%d(View_Zen_Idx)%cloud_prof14(Tropo_Idx_NWP)
            Blkbdy_Tropo_Rad_Chn14 = rtm(X_NWP_Idx,Y_NWP_Idx)%d(View_Zen_Idx)%ch(31)%Rad_BB_Cloud_Profile(Tropo_Idx_NWP)

            !
            !---Tropopause Emissivity
            !
            Emiss_Tropo_Chn14(Elem_Idx,Line_Idx) =  &
                  (Rad_Chn14 - Clr_Rad_Chn14) / (Blkbdy_Tropo_Rad_Chn14 - Clr_Rad_Chn14)

      END IF
    END DO Element_Loop
  END DO Line_Loop

 END SUBROUTINE Compute_Emiss_Tropo_Chn14

!-----------------------------------------------------------------
!
!-----------------------------------------------------------------

SUBROUTINE compute_spatial_uniformity(dx, dy, space_mask, data, data_mean, data_max, data_min, data_uni)

  INTEGER (kind=int4), intent(in) :: dx, dy
  LOGICAL, intent(in), dimension(:,:) :: space_mask
  REAL (kind=real4), intent(in), dimension(:,:) :: data
  REAL (kind=real4), intent(out), dimension(:,:), allocatable :: data_mean, data_max, data_min, data_uni
  INTEGER (kind=int4) :: nx, ny, nx_uni, ny_uni, nsub, astatus
  INTEGER (kind=int4) :: ielem, iline, ielem1, ielem2, iline1, iline2, n_good
  REAL (kind=real4), dimension(:,:), allocatable :: temp
  INTEGER (kind=int4), dimension(:,:), allocatable :: good
  LOGICAL, dimension(:,:), allocatable :: space_mask_temp

  nx = size(data,1)
  ny = size(data,2)

  nx_uni = 2*dx + 1
  ny_uni = 2*dy + 1
  nsub = nx_uni*ny_uni

  allocate(temp(nx_uni,ny_uni), good(nx_uni,ny_uni), &
           space_mask_temp(nx_uni,ny_uni), data_mean(nx,ny), data_max(nx,ny), &
           data_min(nx,ny),data_uni(nx,ny),stat=astatus)
  if (astatus /= 0) then
    print "(a,'Not enough memory to allocate spatial uniformity arrays.')",EXE_PROMPT
    stop
  endif

  line_loop: do iline=1, ny

    iline1 = max(1,iline-dy)
    iline2 = min(ny,iline+dy)
    ny_uni = (iline2 - iline1) + 1

    element_loop: do ielem=1, nx

      data_mean(ielem,iline) = missing_value_real4
      data_max(ielem,iline) = missing_value_real4
      data_min(ielem,iline) = missing_value_real4
      data_uni(ielem,iline) = missing_value_real4

      if (.NOT. space_mask(ielem,iline)  .and. &
          data(ielem,iline) /= missing_value_real4) then

        ielem1 = max(1,ielem-dx)
        ielem2 = min(nx,ielem+dx)
        nx_uni = (ielem2 - ielem1) + 1

        space_mask_temp = .TRUE.
        temp(1:nx_uni,1:ny_uni) = data(ielem1:ielem2,iline1:iline2)
        space_mask_temp(1:nx_uni,1:ny_uni) = space_mask(ielem1:ielem2,iline1:iline2)
        n_good = COUNT(.NOT. space_mask_temp)

        if (n_good > 0) then
            where (space_mask_temp)
               temp = 0.
            end where
         
          data_mean(ielem,iline) =  sum(temp(1:nx_uni,1:ny_uni))/n_good
          data_uni(ielem,iline) = sqrt(max(0.0,(sum((temp(1:nx_uni,1:ny_uni))**2)/n_good - data_mean(ielem,iline)**2)))
          data_max(ielem,iline) = maxval(temp(1:nx_uni,1:ny_uni))
          data_min(ielem,iline) = minval(temp(1:nx_uni,1:ny_uni))
        endif

      endif

    end do element_loop
  end do line_loop

  deallocate(temp, good, space_mask_temp,stat=astatus)
  if (astatus /= 0) then
    print "(a,'Error deallocating temporary spatial uniformity arrays.')",EXE_PROMPT
    stop
  endif

END SUBROUTINE compute_spatial_uniformity

!-----------------------------------------------------------------
!
!-----------------------------------------------------------------

SUBROUTINE gradient2d(grid, nx, ny, mask, min_valid, max_valid, threshold_value, xmax, ymax, num_steps)

  REAL (kind=real4), dimension(:,:), intent(in) :: grid
  INTEGER (kind=int4), intent(in) :: nx, ny
  INTEGER (kind=int1), dimension(:,:), intent(in) :: mask
  REAL (kind=real4), intent(in) :: min_valid, max_valid, threshold_value
  INTEGER (kind=int4), dimension(:,:), intent(inout) :: xmax, ymax
  INTEGER (kind=int4), dimension(:,:), intent(inout), optional :: num_steps
  
  INTEGER (kind=int4), parameter :: max_step = 150
  INTEGER (kind=int4) :: ielem, iline, im, jm, ip, jp, dx_start, dy_start, index
  INTEGER (kind=int4) :: direction, di, dj, i0, j0, i1, j1, ibad
  REAL (kind=real4) :: min_grad, ref_value
  INTEGER (kind=int4), dimension(8) :: icol, irow, di_default, dj_default
  REAL (kind=real4), dimension(8) :: grad
  
  dx_start = 2
  dy_start = 2
  
  di_default = (/0,1,1,1,0,-1,-1,-1/)
  dj_default = (/-1,-1,0,1,1,1,0,-1/)
  
  xmax = missing_value_int4
  ymax = missing_value_int4
  
  if (present(num_steps)) then
    num_steps = missing_value_int4
  endif
  
  line_loop: do iline=1, ny
    
    jm = max(1,iline-dy_start)
    jp = min(ny,iline+dy_start)    
    
    element_loop: do ielem=1, nx
    
      if (mask(ielem,iline) == sym%YES .and. grid(ielem,iline) >= min_valid .and. grid(ielem,iline) <= max_valid) then
      
        im = max(1,ielem-dx_start)
        ip = min(nx,ielem+dx_start)
        
        icol = (/ielem,ip,ip,ip,ielem,im,im,im/)        
        irow = (/jm,jm,iline,jp,jp,jp,iline,jm/)
        
        direction = -999
        min_grad = 99999.0
        do index=1, 8 
          if (grid(icol(index),irow(index)) >= min_valid .and. grid(icol(index),irow(index)) <= max_valid) then
            grad(index) = grid(ielem,iline) - grid(icol(index),irow(index))
            if (grad(index) < min_grad) then
              min_grad = grad(index)
              direction = index
            endif
          endif
        end do
        
        if (direction >= 1 .and. direction <= 8) then
        
          di = di_default(direction)
          dj = dj_default(direction)
        
          do index = 1, max_step
            i0 = max(1,min(ielem + di*index,nx))
            j0 = max(1,min(iline + dj*index,ny))
            i1 = max(1,min(ielem + di*index + di,nx))
            j1 = max(1,min(iline + dj*index + dj,ny))
          
            ref_value = grid(i0,j0)
            ibad = 1
            if (grid(i1,j1) >= min_valid .and. grid(i1,j1) <= max_valid) then 
              ibad = 0
            endif
          
            if (grid(i1,j1) >= threshold_value .or. index == max_step .or. grid(i1,j1) < ref_value .or. ibad == 1) then
              xmax(ielem,iline) = i0
              ymax(ielem,iline) = j0
              if (present(num_steps)) then
                num_steps(ielem,iline) = index
              endif
              exit
            endif
          end do
          
        endif
        
      endif
    
    end do element_loop    
  end do line_loop

END SUBROUTINE gradient2d

!--------------------------------------------------------------------------
! Matrix Inversion for a 3x3 matrix
!--------------------------------------------------------------------------
subroutine INVERT_3x3(A,A_inv,ierr)
  real, dimension(:,:), intent(in):: A
  real, dimension(:,:), intent(out):: A_inv
  integer, intent(out):: ierr
  real:: determinant

  ierr = 0
!--- compute determinant
  determinant = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3)) - &
                A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3)) + &
                A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
  if (determinant == 0.0) then
!       print *, "Singular Matrix in Invert 3x3"
        ierr = 1
  endif

!--- compute inverse
  A_inv(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
  A_inv(1,2) = A(1,3)*A(3,2) - A(3,3)*A(1,2)
  A_inv(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)
  A_inv(2,1) = A(2,3)*A(3,1) - A(3,3)*A(2,1)
  A_inv(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
  A_inv(2,3) = A(1,3)*A(2,1) - A(2,3)*A(1,1)
  A_inv(3,1) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
  A_inv(3,2) = A(1,2)*A(3,1) - A(3,2)*A(1,1)
  A_inv(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)
  A_inv = A_inv / determinant


end subroutine INVERT_3x3
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE prof_lookup_using_p(zlev, plev, tlev, z, p, t, ilev, a)

  REAL (kind=real4), intent(in), dimension(:) :: zlev, plev, tlev
  REAL (kind=real4), intent(in) :: p
  REAL (kind=real4), intent(out) :: z, t
  INTEGER (kind=int4), intent(out) :: ilev
  REAL (kind=real4), intent(out) :: a
  INTEGER (kind=int4) :: nlev
  REAL (kind=real4) :: dp, dt, dz

  nlev = size(plev,1)

  call LOCATE(plev,nlev,p,ilev)
  ilev = max(1,min(nlev-1,ilev))

  dp = plev(ilev+1) - plev(ilev)
  dt = tlev(ilev+1) - tlev(ilev)
  dz = zlev(ilev+1) - zlev(ilev)

  if (dp /= 0.0) then
    a = (p - plev(ilev))/dp
    t = tlev(ilev) + a*dt
    z = zlev(ilev) + a*dz
  else
    a = 0.0
    t = tlev(ilev)
    z = zlev(ilev)
  endif

END SUBROUTINE prof_lookup_using_p
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
SUBROUTINE prof_lookup_using_t(zlev, plev, tlev, z, p, t, ilev, a)

  REAL (kind=real4), intent(in), dimension(:) :: zlev, plev, tlev
  REAL (kind=real4), intent(in) :: t
  REAL (kind=real4), intent(out) :: p, z
  INTEGER (kind=int4), intent(out) :: ilev
  REAL (kind=real4), intent(out) :: a
  INTEGER (kind=int4) :: nlev

  nlev = size(plev,1)

  !Make sure that profile is tropospheric
  if (t < tlev(1) .and. tlev(1) < tlev(nlev)) then
    CALL prof_lookup_using_t_lapse(zlev, plev, tlev, z, p, t, ilev, a)
    !print*,zlev
    !print*,tlev
    !print*,plev
    !print*,z,p,t,ilev,a
    !stop
  else
    CALL prof_lookup_using_t_prof(zlev, plev, tlev, z, p, t, ilev, a)
  endif

END SUBROUTINE prof_lookup_using_t
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
SUBROUTINE prof_lookup_using_z(zlev, plev, tlev, z, p, t, ilev, a)

  REAL (kind=real4), intent(in), dimension(:) :: zlev, plev, tlev
  REAL (kind=real4), intent(in) :: z
  REAL (kind=real4), intent(out) :: p, t
  INTEGER (kind=int4), intent(out) :: ilev
  REAL (kind=real4), intent(out) :: a
  INTEGER (kind=int4) :: nlev
  REAL (kind=real4) :: dp, dt, dz

  nlev = size(plev,1)

  call LOCATE(zlev,nlev,z,ilev)
  ilev = max(1,min(nlev-1,ilev))

  dp = plev(ilev+1) - plev(ilev)
  dt = tlev(ilev+1) - tlev(ilev)
  dz = zlev(ilev+1) - zlev(ilev)

  if (dp /= 0.0) then
    a = (z - zlev(ilev))/dz
    p = plev(ilev) + a*dp
    t = tlev(ilev) + a*dt
  else
    a = 0.0
    p = plev(ilev)
    t = tlev(ilev)
  endif

END SUBROUTINE prof_lookup_using_z
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
SUBROUTINE prof_lookup_using_t_lapse(zlev, plev, tlev, z, p, t, ilev, a)

  REAL (kind=real4), intent(in), dimension(:) :: zlev, plev, tlev
  REAL (kind=real4), intent(in) :: t
  REAL (kind=real4), intent(out) :: p, z
  INTEGER (kind=int4), intent(out) :: ilev
  REAL (kind=real4), intent(out) :: a
  REAL (kind=real4) :: dt, dz, t_tmp
  REAL (kind=real4), PARAMETER :: DT_STRATO_LAPSE_RATE = -6.5 !k
  REAL (kind=real4), PARAMETER :: DZ_STRATO_LAPSE_RATE = 1000.0 !m

  dt = DT_STRATO_LAPSE_RATE
  dz = DZ_STRATO_LAPSE_RATE

  a = (t - tlev(1))/dt
  z = zlev(1) + a*dz

  !This is a hack, really should be using full profile
  CALL prof_lookup_using_z(zlev, plev, tlev, z, p, t_tmp, ilev, a)

END SUBROUTINE prof_lookup_using_t_lapse

SUBROUTINE prof_lookup_using_t_prof(zlev, plev, tlev, z, p, t, ilev, a)

  REAL (kind=real4), intent(in), dimension(:) :: zlev, plev, tlev
  REAL (kind=real4), intent(in) :: t
  REAL (kind=real4), intent(out) :: p, z
  INTEGER (kind=int4), intent(out) :: ilev
  REAL (kind=real4), intent(out) :: a
  INTEGER (kind=int4) :: nlev
  REAL (kind=real4) :: dp, dt, dz

  nlev = size(plev,1)

  call LOCATE(tlev,nlev,t,ilev)
  ilev = max(1,min(nlev-1,ilev))

  dp = plev(ilev+1) - plev(ilev)
  dt = tlev(ilev+1) - tlev(ilev)
  dz = zlev(ilev+1) - zlev(ilev)

  if (dp /= 0.0) then
    a = (t - tlev(ilev))/dt
    p = plev(ilev) + a*dp
    z = zlev(ilev) + a*dz
  else
    a = 0.0
    p = plev(ilev)
    z = zlev(ilev)
  endif

END SUBROUTINE prof_lookup_using_t_prof

!-----------------------------------------------------------------
!
!-----------------------------------------------------------------
SUBROUTINE destroy_spatial_uniformity(data_mean, data_max, data_min, data_uni)

  REAL (kind=real4), intent(inout), dimension(:,:), allocatable :: data_mean, data_max, data_min, data_uni
  INTEGER (kind=int4) :: astatus

  deallocate(data_mean, data_max, data_min, data_uni,stat=astatus)
  if (astatus /= 0) then
    print "(a,'Error deallocating spatial uniformity arrays.')",EXE_PROMPT
    stop
  endif

END SUBROUTINE destroy_spatial_uniformity

!-------------------------------------------------------------------------
! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
!--------------------------------------------------------------------------
  subroutine LOCATE_FLOAT32(xx, n, x, j)

!   Arguments
    integer,                        intent(in)  :: n
    integer,                        intent(out) :: j
    real (kind=ipre),               intent(in)  :: x
    real (kind=ipre), dimension(:), intent(in)  :: xx

!   Local variables
    integer :: i, jl, jm, ju

    jl = 0
    ju = n + 1
    do i = 1, 2*n
       if (ju-jl <= 1) then
          exit
       endif
       jm = (ju + jl) / 2
       if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if (x == xx(1)) then
       j=1
    else if (x == xx(n)) then
       j = n - 1
    else
       j = jl
    endif

  end subroutine LOCATE_FLOAT32

!-------------------------------------------------------------------------
! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
!--------------------------------------------------------------------------
  subroutine LOCATE_INT32(xx, n, x, j)

!   Arguments
    integer(kind=int4),               intent(in)  :: n
    integer(kind=int4),               intent(out) :: j
    integer(kind=int4),               intent(in)  :: x
    integer(kind=int4), dimension(:), intent(in)  :: xx

!   Local variables
    integer(kind=int4) :: i, jl, jm, ju

    jl = 0
    ju = n + 1
    do i = 1, 2*n
       if (ju-jl <= 1) then
          exit
       endif
       jm = (ju + jl) / 2
       if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if (x == xx(1)) then
       j=1
    else if (x == xx(n)) then
       j = n - 1
    else
       j = jl
    endif

  end subroutine LOCATE_INT32



!----------------------------------------------------------------------
! END OF THIS MODULE
!----------------------------------------------------------------------
END MODULE Baseline_Cloud_Height
