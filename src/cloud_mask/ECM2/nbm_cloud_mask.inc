 !$Id: nbm_cloud_mask.inc 1919 2016-12-12 16:51:17Z heidinger $
 !--------------------------------------------------------------
 ! parameter used in naive_bayesian_cloud_mask
 !--------------------------------------------------------------

 !--- Set this to 1 to peform diagnostic calculations
 integer, parameter, private:: Do_By_Class_Flag = 1

 !--- Number of Tests
 integer, parameter, private:: NUMBER_OF_CLOUD_FLAGS = 35
 !integer, parameter, private:: NUMBER_OF_FLAG_BYTES = 7
 integer, parameter, private:: NUMBER_OF_NONCLOUD_FLAGS = 19

 !--- Angular and Temperature Limit Thresholds
 real, parameter, private:: Glint_Zen_Thresh = 45.0
 real, parameter, private:: Reflectance_Gross_Solzen_Thresh = 75.0     !was 70.0
 real, parameter, private:: Reflectance_Spatial_Solzen_Thresh = 80.0
 real, parameter, private:: Reflectance_Gross_Airmass_Thresh = 7.0  !100.0 !5.0   !turned off
 real, parameter, private:: Reflectance_Gross_Lunzen_Thresh = 80.0
 real, parameter, private:: Radiance_Lunar_City_Thresh = 2.5e-08
 real, parameter, private:: Emiss_375um_Day_Solzen_Thresh = 80.0       !was 85.0
 real, parameter, private:: Emiss_375um_Night_Solzen_Thresh = 90.0     !was 80.0
 real, parameter, private:: Bt_11um_Cold_Scene_Thresh = 230.0
 real, parameter, private:: Tsfc_Cold_Scene_Thresh = 230.0
 real, parameter, private:: Path_Tpw_Dry_Scene_Thresh = 0.5
 real, parameter, private:: Bt_375um_Cold_Scene_Thresh = 240.0
 real, parameter, private:: Forward_Scatter_Scatzen_Max_Thresh = 95.0
 real, parameter, private:: Forward_Scatter_Solzen_Max_Thresh = 95.0

 !--- EUMETCAST fire detection parameters
 real, private, parameter :: EUMETCAST_FIRE_DAY_SOLZEN_THRESH = 70.0
 real, private, parameter :: EUMETCAST_FIRE_NIGHT_SOLZEN_THRESH = 90.0

 real, private, parameter :: BT_375UM_EUMET_FIRE_DAY_THRESH = 310.0
 real, private, parameter :: BT_DIFF_EUMET_FIRE_DAY_THRESH = 8.0
 real, private, parameter :: STDDEV_11UM_EUMET_FIRE_DAY_THRESH = 1.0
 real, private, parameter :: STDDEV_375UM_EUMET_FIRE_DAY_THRESH = 4.0

 real, private, parameter :: BT_375UM_EUMET_FIRE_NIGHT_THRESH = 290.0
 real, private, parameter :: BT_DIFF_EUMET_FIRE_NIGHT_THRESH = 0.0
 real, private, parameter :: STDDEV_11UM_EUMET_FIRE_NIGHT_THRESH = 1.0
 real, private, parameter :: STDDEV_375UM_EUMET_FIRE_NIGHT_THRESH = 4.0

