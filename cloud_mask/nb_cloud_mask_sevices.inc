 !$Id: nb_cloud_mask.inc 1919 2016-12-12 16:51:17Z heidinger $
 !--------------------------------------------------------------
 ! parameter used in services and addons module
 !--------------------------------------------------------------

 integer, parameter, private:: NUMBER_OF_FLAG_BYTES = 7

 !--- EUMETCAST fire detection parameters
 real, parameter, private:: Reflectance_Spatial_Solzen_Thresh = 80.0
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

 !--- CLAVRx Smoke Thresholds
 real, private, parameter:: Refl_065_Min_Smoke_Water_Thresh = 2.0
 real, private, parameter:: Refl_065_Max_Smoke_Water_Thresh = 25.0
 real, private, parameter:: Refl_138_Max_Smoke_Water_Thresh = 5.0
 real, private, parameter:: Refl_160_Max_Smoke_Water_Thresh = 5.0
 real, private, parameter:: Refl_375_Max_Smoke_Water_Thresh = 3.0
 real, private, parameter:: Emiss_11_Tropo_Max_Smoke_Water_Thresh = 0.05
 real, private, parameter:: T11_Std_Max_Smoke_Water_Thresh = 0.25
 real, private, parameter:: Refl_065_Std_Max_Smoke_Water_Thresh = 3.00
 real, private, parameter:: Btd_4_11_Max_Smoke_Water_Thresh = 3.0
 real, private, parameter:: Btd_11_12_Max_Smoke_Water_Thresh = 2.0
 real, private, parameter:: Solzen_Max_Smoke_Water_Thresh = 80.0

 real, private, parameter:: Refl_065_Min_Smoke_Land_Thresh = 10.0
 real, private, parameter:: Refl_065_Max_Smoke_Land_Thresh = 25.0
 real, private, parameter:: Refl_138_Max_Smoke_Land_Thresh = 5.0
 real, private, parameter:: NIR_Smoke_Ratio_Max_Land_Thresh = 0.0
 real, private, parameter:: Refl_375_Max_Smoke_Land_Thresh = 3.0
 real, private, parameter:: Emiss_11_Tropo_Max_Smoke_Land_Thresh = 0.3
 real, private, parameter:: T11_Std_Max_Smoke_Land_Thresh = 2.0
 real, private, parameter:: Refl_065_Std_Max_Smoke_Land_Thresh = 3.00
 real, private, parameter:: Btd_4_11_Max_Smoke_Land_Thresh = 5.0
 real, private, parameter:: Btd_11_12_Max_Smoke_Land_Thresh = 2.0
 real, private, parameter:: Solzen_Max_Smoke_Land_Thresh = 80.0

 !--- CLAVRx Dust Thresholds
 real, private, parameter:: Btd_11_12_Metric_Max_Dust_Thresh = -0.5 !-0.75! -0.5 !-1.0
 real, private, parameter:: Btd_11_12_Max_Dust_Thresh = 0.0
 real, private, parameter:: Bt_11_Std_Max_Dust_Thresh = 1.0
 real, private, parameter:: Bt_11_12_Clear_Diff_Max_Dust_Thresh = -0.5
 real, private, parameter:: Emiss_11_Tropo_Max_Dust_Thresh = 0.20
 real, private, parameter:: Emiss_11_Tropo_Min_Dust_Thresh = 0.01
 real, private, parameter:: Bt_11_Clear_Diff_Min_Dust_Thresh = -5.0
 real, private, parameter:: Btd_85_11_Max_Dust_Thresh = -0.5
 real, private, parameter:: Btd_85_11_Min_Dust_Thresh = -2.5 !-1.5

 !--- CLAVRx Thin Cirrus Thresholds
 real, private, parameter:: Thin_Cirr_Min_Thresh = 1.0
 real, private, parameter:: Thin_Cirr_Max_Thresh = 3.0
 real, private, parameter:: Desert_Thin_Cirr_Min_Thresh = 3.0
 real, private, parameter:: Desert_Thin_Cirr_Max_Thresh = 5.0
 real, private, parameter:: Snow_Thin_Cirr_Min_Thresh = 3.0
 real, private, parameter:: Snow_Thin_Cirr_Max_Thresh = 5.0
