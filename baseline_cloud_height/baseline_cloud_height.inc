!$Id: baseline_cloud_height.inc 2624 2018-04-04 09:28:12Z awalther $
!-------------------------------------------------------------------------------
!CVS SRC: baseline_cloud_height.inc,v 1.7 2010/11/08 17:31:48 heidinger Exp $
!-------------------------------------------------------------------------------
!
! NAME:
!   cld_height_11_12_13_1dvar.inc
!
! FUNCTION:
!   Fortran include file.
!
! DESCRIPTION:
!   This include file contains PARAMETERs needed by the cld_height_11_12_13_1var.f90
!   code.
!
! REFERENCE:
!
! CALLING SEQUENCE:
!
! INPUTS:
!
! OUTPUTS:
!
! DEPENDENCIES:
!
! RESTRICTIONS:
!
! HISTORY:
!
!-------------------------------------------------------------------------------

  !--- Limit and Bounds
  REAL(KIND=REAL4), PRIVATE, PARAMETER ::  SENSOR_ZEN_THRESH = 70.0

  !--- Processing Flags
  INTEGER(KIND=INT4), PARAMETER:: USE_LRC_FLAG =  1               !sym%YES            !turn on lrc
  INTEGER(KIND=INT4), PARAMETER:: USE_BIAS_CORRECTION_FLAG = 0    !sym%NO !turn on bias corr.
  INTEGER(KIND=INT4), PARAMETER:: USE_SFC_INVER_FLAG =    1       !sym%YES      !turn sfc inver logic
  INTEGER(KIND=INT4), PARAMETER:: USE_LOWER_CLD_INTERP_FLAG = 1   !sym%YES      !turn sfc inver logic
 
  !--- Product Quality Flags
  INTEGER, PRIVATE, PARAMETER :: VALID_CTH_RETRIEVAL = 0
  INTEGER, PRIVATE, PARAMETER :: INVALID_CTH_DUE_TO_SPACE_PIX = 1
  INTEGER, PRIVATE, PARAMETER :: INVALID_CTH_OUTSIDE_SEN_ZEN_RANGE = 2
  INTEGER, PRIVATE, PARAMETER :: INVALID_CTH_BAD_DATA = 3
  INTEGER, PRIVATE, PARAMETER :: INVALID_CTH_CLD_MASK_CLR = 4
  INTEGER, PRIVATE, PARAMETER :: INVALID_CTH_INVALID_CLD_TYPE = 5
  INTEGER, PRIVATE, PARAMETER :: INVALID_CTH_FAILED_RETRIEVAL = 6

  !--- Parameter Quality Flags
  !--- Parameter Quality Flags
  !--- these are used to assign a quality indicator for each parameter.  The
  !--- highest quality flag is 3 and the lowest quality is 0
  !--- quality is determined by the ratio of estimated uncertainity relative
  !--- to that assume for the apriori (first guess).
  INTEGER, PRIVATE, PARAMETER :: CTH_PARAM_FAILED_RETREVIAL = 0
  INTEGER, PRIVATE, PARAMETER :: CTH_PARAM_LOW_QUALITY_RETREVIAL = 1  !retrived uncert. > 2/3 of apriori uncert.
  INTEGER, PRIVATE, PARAMETER :: CTH_PARAM_2_3_APRIORI_RETREVIAL = 2  !retrieved uncert. < 2/3 of apriori uncert.
  INTEGER, PRIVATE, PARAMETER :: CTH_PARAM_1_3_APRIORI_RETREVIAL = 3  !retrieved uncert. < 1/3 of apriori uncert.

  !--- Parametes to control some internal computations (LRC and Multi-layer cloud interpolation)
  INTEGER, PRIVATE, PARAMETER ::INTERP_LOWER_CLOUD_PIXEL_RADIUS = 5
  REAL(KIND=REAL4), PRIVATE, PARAMETER ::EMISS_TROPO_CHN14_GRADIENT_MIN = 0.0
  REAL(KIND=REAL4), PRIVATE, PARAMETER ::EMISS_TROPO_CHN14_GRADIENT_MAX = 1.0
  REAL(KIND=REAL4), PRIVATE, PARAMETER ::EMISS_TROPO_CHN14_GRADIENT_THRESH = 0.75

  !--- number of meta-data flags for this algorithm
  INTEGER(kind=int4), PARAMETER :: NUM_META_DATA = 8          !number of meta data flags
  INTEGER(kind=int4), PARAMETER :: NUM_META_DATA_PACKED_BYTES = 1 !number of bytes to pack meta data flags

  !--- these PARAMETERs control the size of matrices
  INTEGER(kind=int4), PARAMETER :: NUM_OBS = 3                !number of observations
  INTEGER(kind=int4), PARAMETER :: NUM_PARAM = 3              !number of retrieved PARAMETERs

  !--- maximum number of iterations allowed
  INTEGER(kind=int4), PARAMETER :: ITER_MAX = 10              !maximum number of iterations

  !--- control size of window in spatial uniformity
  INTEGER(KIND=INT4), PARAMETER:: ELEM_WIDTH_3X3 = 1          !width of element for 3x3 calcs
  INTEGER(KIND=INT4), PARAMETER:: LINE_WIDTH_3X3 = 1          !width of lines for 3x3 calcs

  !--- PARAMETERs that control the bottom-up estimation of Zc and Pc
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: LAPSE_RATE_OCEAN = 9.8e-03    !dry lapse rate over ocean for height K/km
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: MIN_P_INVERSION = 700.0  !min pressure to look for low-level inversion
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: DELTA_PSFC_INVERSION = 50.0  !pressure above sfc to look for low-level inversion


  !--- the PARAMETERs provide the apriori values and their uncertainties
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: TC_AP_TROPO_OFFSET_CIRRUS = 15.0      !apriori Tc for opaque clouds
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: TC_AP_UNCER_OPAQUE  = 10.0            !apriori uncertainty of Tc for opaque clouds
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: TC_AP_UNCER_CIRRUS = 20.0             !apriori uncertainty of Tc for cirrus

  REAL(KIND=REAL4), PARAMETER, PRIVATE:: EC_AP_UNCER_OPAQUE = 0.2             !apriori uncertainty of ec for opaque clouds
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: EC_AP_UNCER_CIRRUS = 0.4             !apriori uncertainty of ec for cirrus clouds

  REAL(KIND=REAL4), PARAMETER, PRIVATE:: BETA_AP_WATER = 1.3                  !aprio estimate of beta for water clouds
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: BETA_AP_UNCER_WATER = 0.2            !apriori uncertainty of  beta for ice clouds
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: BETA_AP_ICE = 1.1                    !apriori estimate of beta for ice clouds
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: BETA_AP_UNCER_ICE = 0.2              !apriori uncertainty of  beta for water clouds

  REAL(KIND=REAL4), PARAMETER, PRIVATE:: TAU_AP_FOG_TYPE = 1.2               !apriori estimate of tau for fog cloud type
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: TAU_AP_WATER_TYPE = 2.3             !apriori estimate of tau for water cloud type
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: TAU_AP_MIXED_TYPE = 2.3             !apriori estimate of tau for mixed cloud type
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: TAU_AP_OPAQUE_ICE_TYPE = 2.3        !apriori estimate of tau for opaque ice cloud type
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: TAU_AP_CIRRUS_TYPE = 0.9            !apriori estimate of tau for cirrus cloud type
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: TAU_AP_OVERLAP_TYPE = 2.0           !apriori estimate of tau for multilayer cloud type

  !--- specify calibration errors (units=K)
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11_CAL_UNCER = 1.0       !calibration uncertainity for 11 micron bt
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11_12_CAL_UNCER = 1.0    !calibration uncertainity for 11-12 micron bt difference
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11_13_CAL_UNCER = 2.0    !calibration uncertainity for 11-13 micron bt difference

  !--- specify clear-sky model errors (units=K)
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11_CLR_UNCER_WATER = 1.5     !uncertainty for 11 micron clear-sky bt estimate over water
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11_12_CLR_UNCER_WATER = 0.5  !uncertainty for 11-12 micron clear-sky bt diff estimate over water
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11_13_CLR_UNCER_WATER = 4.0  !uncertainty for 11-13 micron clear-sky bt diff estimate over water

  REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11_CLR_UNCER_LAND = 5.0     !uncertainty for 11 micron clear-sky bt estimate over land
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11_12_CLR_UNCER_LAND = 1.0  !uncertainty for 11-12 micron clear-sky bt estimate over land
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11_13_CLR_UNCER_LAND = 4.0  !uncertainty for 11-13 micron clear-sky bt estimate over land

  !--- these PARAMETERs contrain the alloweable solutions
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: MIN_ALLOWABLE_TC = 160.0  !minimum allowable value for cloud-top temperature in K
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: PC_LOWER_CLOUD_OFFSET = 200.0  !hpa


!----------------------------------------------------------------------
!--- Note the coefficients to relate beta_11_12 to beta_11_13 and 
!--- that relate beta_11_12 to re are a function of ice cruystal
!--- habit.  I recommend using the aggregates values as described in ATBD.
!--- if another crystal habit is desired, modify comments consistently.
!--- ideally, the crystal habit should be become an argument to this
!--- routine.  (A. Heidinger 3/2008).
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!--- these coefficients relate beta_11_13 to beta_11_12
!----------------------------------------------------------------------

!--- water clouds derived from Mie theory by A. Heidinger
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_FIT_WATER = -0.217447
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_FIT_WATER = 1.250322

!--- these values are for ice clouds from Ping Yang's Database (10-100um)

! !- aggregates
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_FIT_ICE = -0.325909
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_FIT_ICE = 1.340083

! !- bullet rosette
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_FIT_ICE = -1.132466
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_FIT_ICE = 1.137529

! !- droxtal
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_FIT_ICE = -0.446551
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_FIT_ICE = 1.463324

! !- hollow column
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_FIT_ICE = -0.505943
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_FIT_ICE = 1.512131

  !- plates
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_FIT_ICE = -0.399185
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_FIT_ICE = 1.410593

! !- solid column
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_FIT_ICE = -0.604064
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_FIT_ICE = 1.605465

! !- spheriod
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_FIT_ICE = -0.418182
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_FIT_ICE = 1.437576

!----------------------------------------------------------------------
!--- these coefficients convert beta_11_12 to particle radius 
!----------------------------------------------------------------------
  !-- water spheres from Mie theory as run by A. Heidinger
  REAL(KIND=REAL4), PARAMETER, PRIVATE, dimension(2):: &
                    BETA2RE_COEF_WATER =  (/ -0.387010, 0.396507/)

  !-- Ping Yang's data

! !- aggregates
  REAL(KIND=REAL4), PARAMETER, PRIVATE, dimension(2):: &
      BETA2RE_COEF_ICE =  (/ -0.799525, 0.804073/)

! !- bullet rosette
! REAL(KIND=REAL4), PARAMETER, PRIVATE, dimension(2):: &
!     BETA2RE_COEF_ICE =  (/ -1.1324662, 1.1375290/)

! !- droxtal
! REAL(KIND=REAL4), PARAMETER, PRIVATE, dimension(2):: &
!     BETA2RE_COEF_ICE =  (/ -0.581726, 0.583603/)

! !- hollow column
! REAL(KIND=REAL4), PARAMETER, PRIVATE, dimension(2):: &
!     BETA2RE_COEF_ICE =  (/ -0.681854, 0.684950/)

  !- plate
! REAL(KIND=REAL4), PARAMETER, PRIVATE, dimension(2):: &
!     BETA2RE_COEF_ICE =  (/ -1.033499, 1.043219/)

! !- solid column
! REAL(KIND=REAL4), PARAMETER, PRIVATE, dimension(2):: &
!     BETA2RE_COEF_ICE =  (/ -0.673415, 0.674144/)

! !- spheroid
! REAL(KIND=REAL4), PARAMETER, PRIVATE, dimension(2):: &
!     BETA2RE_COEF_ICE =  (/ -0.720051, 0.723820/)


!----------------------------------------------------------------------
!--- these PARAMETERs control the cloud layer
!----------------------------------------------------------------------

  !--- maximum pressure of a high cloud 
  REAL(KIND=REAL4), PARAMETER, PRIVATE :: PC_HIGH_CLOUD_MAX = 440.0 !hPa

  !--- maximum pressure of a mid-level cloud
  REAL(KIND=REAL4), PARAMETER, PRIVATE :: PC_MID_CLOUD_MAX = 680.0  !hPa 
