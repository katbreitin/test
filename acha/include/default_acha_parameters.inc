!$Id: acha_parameters.inc 4092 2021-03-02 22:05:14Z heidinger $
!----------------------------------------------------------------------
!
! Name:
!   acha_parameters.inc
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
! Description:
!   Include file that declares a number of physical constants,
!   apriori estimates and uncertainty values
!
! MULTI_LAYER_LOGIC_FLAG
! 0 - (baseline) just use the multilayer id in cloud type
! 1 - treat all multilayer like cirrus (force single layer)
! 2 - assume all cirrus are multilayer and let acha decide
!
! USE_TYPE_FLAG
! 0 - use input type to set a priori
! 1 - use input phase and emiss tropo to set a priori
! 2 - do not use input phase/type for setting a priori
!
!----------------------------------------------------------------------

 !-----------------------------------------------------------
 ! include other include files to this file
 !-----------------------------------------------------------

 !--- clear-sky covariance terms
 include 'include/acha_clear_sky_covariances.inc'

 !-----------------------------------------------------------
 ! Isolate User Controlled Flags here
 !-----------------------------------------------------------

 !---- DEFAULT
 logical, private, parameter:: FULL_RETRIEVAL = .true.
 logical, private, parameter:: USE_LINEAR_IN_OPD_EMISSION = .true.
 integer(kind=int4), parameter, private:: MULTI_LAYER_LOGIC_Flag = 2
 logical, private, parameter :: USE_CIRRUS_FLAG = .true. 
 logical, private, parameter :: USE_LAPSE_RATE_FLAG = .false. !.true.  
 logical, private, parameter :: USE_LOWER_INTERP_FLAG = .true.  
 logical, private, parameter :: USE_LRC_FLAG = .false.
 integer, parameter, private:: Cloud_Type_Forced = -1 ! -1 = not used, 0-10 use this type
 logical, private, parameter :: USE_SOUNDER_VALUES = .false.
 integer, private, parameter :: USE_TYPE_FLAG = 2               !choice which a priori routine to use
 real, private, parameter:: ICE_EXTINCTION_TUNING_FACTOR = 1.0
 logical, private, parameter :: USE_TYPE_FOR_ICE_PROB = .false. !whether or not use input type for ice_prob_ap
 logical, private, parameter :: CONSTRAIN_ICE_PROB= .false.     !whether or not to allow acha to modify ice_prob_ap
 logical, private, parameter :: BINARY_ICE_PROB= .false.        !whether to use 0/1 ice_prob_ap

 !------------------------------------------------------------
 ! ACHA parameters
 !------------------------------------------------------------
 integer, parameter, private:: DQF_MIN_THRESH = 1
 real(kind=real4), parameter, private:: count_diag_max = 1000   ! maximum number of diagnostic files allowed
 real(kind=real4), private, parameter:: Dt_Dz_Strato = -0.0065 !K/m
 real(kind=real4), private, parameter:: Dp_Dz_Strato = -0.0150 !hpa/m
 real(kind=real4), private, parameter:: Sensor_Zenith_Threshold = 88.0 !70.0

 !--- set starting index value
 integer(kind=INT4), parameter, private:: Element_Idx_Min = 1
 integer(kind=INT4), parameter, private:: Line_Idx_Min = 1
 integer(kind=INT4), parameter, private:: LRC_Meander_Flag = 1

 !--- number of meta data flags
 integer(kind=INT4), parameter, private :: NUM_META_DATA = 8

 !--- quality is determined by the ratio of estimated uncertainity relative
 !--- to that assume for the apriori (first guess).
 integer(kind=INT1), private, parameter :: CTH_PARAM_FAILED_RETREVIAL = 0_int1
 integer(kind=INT1), private, parameter :: CTH_PARAM_LOW_QUALITY_RETREVIAL = 1_int1  !retrived uncert. > 2/3 of apriori uncert.
 integer(kind=INT1), private, parameter :: CTH_PARAM_2_3_APRIORI_RETREVIAL = 2_int1  !retrieved uncert. < 2/3 of apriori uncert.
 integer(kind=INT1), private, parameter :: CTH_PARAM_1_3_APRIORI_RETREVIAL = 3_int1  !retrieved uncert. < 1/3 of apriori uncert.

 !--- ACHA Quality Flag - previously added in per Feb2018 science delivery
 integer, private, parameter :: CTH_DQF_GOOD_RETREVIAL = 0_int1
 integer, private, parameter :: CTH_DQF_MARGINAL_RETREVIAL = 1_int1  
 integer, private, parameter :: CTH_DQF_RETREVIAL_ATTEMPTED = 2_int1  
 integer, private, parameter :: CTH_DQF_BAD_RETREVIAL = 3_int1  
 integer, private, parameter :: CTH_DQF_OPAQUE_RETREVIAL = 4_int1

 !--- LRC Computation Parameters
 integer, private, parameter :: Grad_Flag_LRC = -1    ! 1 = go up, -1 = go down
 integer, private, parameter :: Max_LRC_Distance = 10
 real, private, parameter :: Min_LRC_Jump = 0.5
 real, private, parameter :: Max_LRC_Jump = 10.0
 real, private, parameter :: Min_Bt_110um_LRC = 220.0
 real, private, parameter :: Max_Bt_110um_LRC = 290.0

 !--- Cirrus Box Parameters
 real (kind=real4), parameter, private:: EMISSIVITY_MIN_CIRRUS = 0.7
 integer (kind=int4), parameter, private:: COUNT_MIN_CIRRUS = 50 !5
 integer (kind=int4), parameter, private:: CIRRUS_BOX_WIDTH_KM = 400 !200

 !---- Lower Cloud Parameters
 integer (kind=int4), private, parameter:: COUNT_MIN_LOWER = 3 !5
 integer (kind=int4), private, parameter:: LOWER_BOX_WIDTH_KM = 200 !200
 real(kind=real4), parameter, private:: TC_LOWER_CLOUD_OFFSET = 10.0  !K

 !---------------------------------------------------------------------
 ! retrieval specific parameters
 !---------------------------------------------------------------------
 integer(kind=INT4), private, parameter :: Num_Levels_Rtm_Prof = 101

 !--- these parameters control the size of matrices
 integer(kind=INT4), private, parameter :: NUM_PARAM = 5              !number of retrieved parameters
 integer(kind=INT4), private, parameter :: NUM_PARAM_SIMPLE = 3              !number of retrieved parameters

 !--- maximum number of iterations allowed
 integer(kind=INT4), private, parameter :: Iter_Idx_Max = 10              !maximum number of iterations

 !--- limits on steps taken in retrieval
 real(kind=real4), parameter, private, DIMENSION(NUM_PARAM):: Delta_X_Max = (/20.0,0.1,0.1,5.0,0.2/)

 !--- parameters that control the conversion of Tc uncertainty to Pc and Zc uncertainty
 real(kind=real4), parameter, private:: ABS_LAPSE_RATE_DlnP_DZ_UNCER = 0.0001 !hPa/m
 real(kind=real4), parameter, private:: ABS_LAPSE_RATE_DT_DZ_UNCER = 0.007 !K/m

 !--- parameters that control the bottom-up estimation of Zc and Pc
 real(kind=real4), parameter, private:: MIN_P_INVERSION = 600.0 !min pressure to look for low-level inversion
 real(kind=real4), parameter, private:: MAX_DELTA_T_INVERSION = 20.0 !max temperature difference (surface -cloud) to look for low-level inversion
 real(kind=real4), parameter, private:: DELTA_PSFC_INVERSION = 50.0  !pressure above sfc to look for low-level inversion
 real(kind=real4), parameter, private:: DELTA_T_LAYER_INVERSION = 0.0   !maximum temperature difference to define an inversion layer
 integer(kind=INT1), parameter, private:: ALLOW_STRATOSPHERE_SOLUTION_FLAG = 1_int1  !allow interpolation in stratosphere

 !--- the parameters that provide the apriori values and their uncertainties
!real(kind=real4), parameter, private:: Tc_Ap_Tropo_Offset_Cirrus = 15.0      !apriori Tc for opaque clouds

 real(kind=real4), parameter, private:: Tc_Ap_Uncer_Opaque = 20.0 !10.0            !apriori uncertainty of Tc for opaque clouds

real(kind=real4), parameter, private:: Tc_Ap_Uncer_Cirrus_Default = 30.0  !!           !apriori uncertainty of Tc for cirrus
                                                                          ! Originally it was 20 in the baseline version of ACHA

 real(kind=real4), parameter, private:: Ts_Ap_Uncer_Sfc = 5.0                !apriori uncertainty of Tsfc when not overlap
 real(kind=real4), parameter, private:: Ts_Ap_Uncer_Lower_Cld = 20.0         !apriori uncertainty of Tsfc when lower cloud below
 real(kind=real4), parameter, private:: Ec_Ap_Uncer_Opaque = 0.2 !0.1        !apriori uncertainty of ec for opaque clouds
 real(kind=real4), parameter, private:: Ec_Ap_Uncer_Cirrus = 0.8 !0.4

 real(kind=real4), parameter, private:: Beta_Ap_Water = 1.3
 real(kind=real4), parameter, private:: Beta_Ap_Uncer_Water = 0.2            !apriori uncertainty of  beta for ice
 real(kind=real4), parameter, private:: Beta_Ap_Ice = 1.06
 real(kind=real4), parameter, private:: Beta_Ap_Uncer_Ice = 0.2              !apriori uncertainty of  beta for water

 real(kind=real4), parameter, private:: Tau_Ap_Water_Phase= 3.0             !apriori estimate of tau for water phase clod
 real(kind=real4), parameter, private:: Tau_Ap_Fog_Type = 1.2               !apriori estimate of tau for fog cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Water_Type = 2.3             !apriori estimate of tau for water cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Supercooled_Type = 2.3             !apriori estimate of tau for mixed cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Mixed_Type = 2.3             !apriori estimate of tau for mixed cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Opaque_Ice_Type = 2.3        !apriori estimate of tau for opaque ice cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Cirrus_Type = 0.9            !apriori estimate of tau for cirrus cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Overlap_Type = 0.9 !2.0           !apriori estimate of tau for multilayer cloud type

 !--- specify calibration  errors
!--> real(kind=real4), parameter, private:: T110um_120um_Cal_Uncer = 1.0  !BASELINE V5
!--> real(kind=real4), parameter, private:: T110um_133um_Cal_Uncer = 2.0 !BASELINE V5
 real(kind=real4), parameter, private:: T110um_Cal_Uncer = 1.0
 real(kind=real4), parameter, private:: T110um_038um_Cal_Uncer = 1.0
 real(kind=real4), parameter, private:: T110um_062um_Cal_Uncer = 1.0
 real(kind=real4), parameter, private:: T110um_067um_Cal_Uncer = 1.0
 real(kind=real4), parameter, private:: T110um_073um_Cal_Uncer = 1.0
 real(kind=real4), parameter, private:: T110um_085um_Cal_Uncer = 1.0
 real(kind=real4), parameter, private:: T110um_097um_Cal_Uncer = 1.0
 real(kind=real4), parameter, private:: T110um_104um_Cal_Uncer = 1.0
 real(kind=real4), parameter, private:: T110um_120um_Cal_Uncer = 1.0
 real(kind=real4), parameter, private:: T110um_133um_Cal_Uncer = 2.0 
 real(kind=real4), parameter, private:: T110um_136um_Cal_Uncer = 2.0 
 real(kind=real4), parameter, private:: T110um_139um_Cal_Uncer = 2.0 
 real(kind=real4), parameter, private:: T110um_142um_Cal_Uncer = 2.0 

 !--- these parameters constrain the allowable solutions
 real(kind=real4), parameter, private:: Min_Allowable_Tc = 170.0  !K

 !--- parameter that controls observed negative heights and pressures.
 real(kind=real4), parameter, private:: ZC_FLOOR = 75.0 !Set a bottom limit to Zc.
 real(kind=real4), parameter, private:: ZC_CEILING = 20000.0 !Set a upper limit to Zc.
 real(kind=real4), parameter, private:: PC_FLOOR = 50.0 !Set a bottom limit to Pc.
 real(kind=real4), parameter, private:: PC_CEILING = 1100.0 !Set a upper limit to Pc.
 real(kind=real4), parameter, private:: TC_FLOOR = 160.0 !Set a bottom limit to Tc.
 real(kind=real4), parameter, private:: TC_CEILING = 320.0 !Set a upper limit to Tc.
 real(kind=real4), parameter, private:: EC_FLOOR = -0.2 !Set a bottom limit to Ec.
 real(kind=real4), parameter, private:: EC_CEILING = 1.2 !Set a upper limit to Ec.
 real(kind=real4), parameter, private:: BETA_FLOOR = 0.2 !Set a bottom limit to Beta.
 real(kind=real4), parameter, private:: BETA_CEILING = 2.0 !Set a upper limit to Beta.

 !--- cirrus apriori for ice cloud temperature relative to tropopause
 !---- computed in 10 deg lat bands.  First bin is -90 to -80
 integer, parameter, private:: num_lat_cirrus_ap = 18

!--- old
!real(kind=real4), dimension(num_lat_cirrus_ap), private, parameter:: &
!   TC_CIRRUS_MEAN_LAT_VECTOR = (/  &
!     17.1184, 16.8644, 14.4241, 15.7231, 14.6659, 17.5565, &
!     27.6142, 17.9770, 13.5127, 16.2466, 22.4699, 24.6992, &
!     18.7324, 18.3266, 15.2087, 14.8122, 18.0353, 19.9059/)
!real(kind=real4), dimension(num_lat_cirrus_ap), private, parameter:: &
!   TC_CIRRUS_STDDEV_LAT_VECTOR = (/  &
!     9.58205, 10.7957, 9.96308, 9.75967, 9.11487, 11.3290, &
!     18.0040, 16.2742, 14.0061, 15.0945, 13.7971, 14.6018, &
!     15.5322, 14.5822, 11.3727, 9.90263, 11.5416, 11.9146/)

!--- akh 10/2020
 real(kind=real4), dimension(num_lat_cirrus_ap), private, parameter:: &
    TC_CIRRUS_MEAN_LAT_VECTOR = (/  &
      7.24244, 8.70593, 14.0095, 14.5873, 14.7501, 17.9235, &
      22.4626, 15.1277, 14.4136, 13.9811, 15.1359, 21.6522, &
      18.0842, 14.8196, 13.4467, 13.8617, 14.5733, 14.9950/)

 real(kind=real4), dimension(num_lat_cirrus_ap), private, parameter:: &
    TC_CIRRUS_STDDEV_LAT_VECTOR = (/  &
      15.9868, 14.6861, 11.3253, 9.39002, 9.64427, 11.4220, &
      14.8919, 13.3661, 13.5816, 13.0034, 13.1591, 16.2132, &
      12.8722, 10.6209, 10.0633, 11.0099, 11.1937, 10.8028/)

 !--- parameters that allow for cloud geometrical thickness estimation
 real(kind=real4), parameter, private:: WATER_EXTINCTION = 25.00   ! 1/km

 real(kind=real4), parameter, private:: ICE_EXTINCTION1 = 1.71 
 real(kind=real4), parameter, private:: ICE_EXTINCTION2 = 1.87
 real(kind=real4), parameter, private:: ICE_EXTINCTION3 = 2.24
 real(kind=real4), parameter, private:: ICE_EXTINCTION4 = 2.88
 real(kind=real4), parameter, private:: ICE_EXTINCTION5 = 4.74

 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION1 = 0.13
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION2 = 0.25
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION3 = 0.39
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION4 = 0.55
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION5 = 0.67
