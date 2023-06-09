!$Id: acha_parameters.inc 1500 2016-02-11 17:00:44Z heidinger $
!----------------------------------------------------------------------
!
! Name:
!   ccl_parameters.inc
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
! Description:
!  Include file that declares of CCL algorithm parameters
!
!  Note, CCL supports generation of ASOS.   
!
!  CCL and ASOS have different definitions of BOX_WIDTH and 
!  HIGH, MID and LOW Pressure boundaries.
!
!  CCL_BOX_WIDTH_KM = box dimensions in kilometers.  Box pixel dimensions
!                     determined from pixel resolution in km
!
!  CCL_SPACING_KM = spacing in km between calculations.  If set to 0,
!                   then every pixel is done.  This speeds up processing
!                   If box has many pixels, pixel to pixel variation is small
!
!----------------------------------------------------------------------

 !Baseline CCL - Pressure boundaries are from NCEP
 integer (kind=int4), parameter:: CCL_BOX_WIDTH_KM = 10
 integer (kind=int4), parameter:: CCL_SPACING_KM = 0
 real, private, parameter:: HIGH_CLOUD_MAX_PRESSURE_THRESH = 350.0
 real, private, parameter:: LOW_CLOUD_MIN_PRESSURE_THRESH = 642.0

 integer(KIND=INT4), PRIVATE, PARAMETER :: Num_Levels_Rtm_Prof = 101
 real(kind=real4), parameter, private:: ZC_FLOOR = 75.0

 !ASOS
!integer (kind=int4), parameter:: CCL_BOX_WIDTH_KM = 50
!integer (kind=int4), parameter:: CCL_SPACING_KM = 10
!real, private, parameter:: HIGH_CLOUD_MAX_PRESSURE_THRESH = 400.0
!real, private, parameter:: LOW_CLOUD_MIN_PRESSURE_THRESH = 631.0

 integer (kind=int4), parameter:: COUNT_MIN_CCL = 3 

 !--- ASOS specific parameters
 real, private, parameter:: OVC_CLEAR_FRACTION_THRESH = 0.07
 real, private, parameter:: BKN_CLEAR_FRACTION_THRESH = 0.30
 real, private, parameter:: LOW_CLOUD_FRACTION_THRESH = 0.40
 real, private, parameter:: ASOS_LAT_MIN = 20.0
 real, private, parameter:: ASOS_LAT_MAX = 50.0
 real, private, parameter:: ASOS_LON_MAX = -60.0
 real, private, parameter:: ASOS_LON_MIN = -140.0

 !--- control processing of ASOS logic within CCL
 logical, private, parameter:: ASOS_FLAG = .false.

 !--- missing values
 real, private, parameter:: MISSING_VALUE_REAL4 = -999.0
 integer (kind=int1), private, parameter:: MISSING_VALUE_INTEGER1 = -128_int1
 integer, private, parameter:: MISSING_VALUE_INTEGER4 = -999

 !--- flight level
 integer (kind=int4), parameter:: NUM_FLIGHT_LEVELS = 5
 integer, private, parameter,dimension(NUM_FLIGHT_LEVELS-1):: FLIGHT_LEVELS = (/24000,18000,10000,5000/)

!logical, private, parameter:: USE_BASE = .true.
!logical, private, parameter:: USE_LOWER = .false.

 ! ACHA QF
 integer, private, parameter :: CTH_DQF_GOOD_RETREVIAL = 0_int1
 integer, private, parameter :: CTH_DQF_MARGINAL_RETREVIAL = 1_int1
 integer, private, parameter :: CTH_DQF_RETREVIAL_ATTEMPTED = 2_int1
 integer, private, parameter :: CTH_DQF_BAD_RETREVIAL = 3_int1

 !--- CCL Quality Flag - previously added in per Feb2018 science delivery
 integer, private, parameter :: CCL_DQF_GOOD_RETREVIAL = 0_int1
 integer, private, parameter :: CCL_DQF_DEGRADED_RETREVIAL = 1_int1
 integer, private, parameter :: CCL_DQF_BAD_RETREVIAL = 2_int1
