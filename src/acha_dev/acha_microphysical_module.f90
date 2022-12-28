!$Id:$
module ACHA_MICROPHYSICAL_MODULE

implicit none

integer, parameter, private:: Habit_Idx = 7   ! 1(droxtals),2(solid_bullet_rosettes),3(hollow_bullet_rosettes),  $
                                              ! 4(solid_columns),5(hollow_columns),6(plates),7(aggregate_columns), $
                                              ! 8(small_aggregate_plates),9(large_aggregate_plates),10(two_habit)
private
public:: SETUP_ICE_MICROPHYSICAL_MODEL

!  acha ice model parameters
INTEGER, PARAMETER, PUBLIC:: beta_degree_ice =        3
INTEGER, PARAMETER, PUBLIC:: re_beta_degree_ice =        3
INTEGER, PARAMETER, PUBLIC:: sing_scat_degree_ice =        2

REAL, DIMENSION(0:2), PUBLIC, SAVE:: Qe_006um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_038um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_039um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_062um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_067um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_073um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_085um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_097um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_104um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_133um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_136um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_139um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_110um_142um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: RE_BETA_110um_COEF_ICE
REAL, DIMENSION(0:2), PUBLIC, SAVE:: Qe_110um_COEF_ICE
REAL, DIMENSION(0:2), PUBLIC, SAVE:: wo_110um_COEF_ICE
REAL, DIMENSION(0:2), PUBLIC, SAVE:: g_110um_COEF_ICE

REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_038um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_039um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_062um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_067um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_073um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_085um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_097um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_104um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_133um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_136um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_139um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: BETA_104um_142um_COEF_ICE
REAL, DIMENSION(0:3), PUBLIC, SAVE:: RE_BETA_104um_COEF_ICE
REAL, DIMENSION(0:2), PUBLIC, SAVE:: Qe_104um_COEF_ICE
REAL, DIMENSION(0:2), PUBLIC, SAVE:: wo_104um_COEF_ICE
REAL, DIMENSION(0:2), PUBLIC, SAVE:: g_104um_COEF_ICE

!--- water microphysical model terms
INTEGER, PARAMETER, PUBLIC:: beta_degree_water =        3
INTEGER, PARAMETER, PUBLIC:: re_beta_degree_water =        3
INTEGER, PARAMETER, PUBLIC:: sing_scat_degree_water =        2

REAL, DIMENSION(0:2), PUBLIC:: Qe_006um_COEF_WATER= [   2.39378,  -0.39669, 0.10725]

REAL, DIMENSION(0:3), PUBLIC:: RE_BETA_110um_COEF_WATER= [   0.59356,   1.41647, -1.12240,   0.53016]
REAL, DIMENSION(0:2), PUBLIC:: Qe_110um_COEF_WATER= [  -1.08658,   3.95746, -1.18558]
REAL, DIMENSION(0:2), PUBLIC:: wo_110um_COEF_WATER= [  -0.11394,   0.77357, -0.23835]
REAL, DIMENSION(0:2), PUBLIC:: g_110um_COEF_WATER= [   0.23866,   0.97914, -0.31347]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_038um_COEF_WATER= [   0.66604, -1.23524,   3.55969,  -1.96452]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_039um_COEF_WATER= [   0.69227, -1.28502,   3.68611,  -2.10146]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_062um_COEF_WATER= [   0.96407, 0.78427,   0.08838,  -0.17242]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_067um_COEF_WATER= [   0.96341, 0.29453,  -0.03356,  -0.07573]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_073um_COEF_WATER= [   0.96626, 0.15321,  -0.06493,  -0.02253]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_085um_COEF_WATER= [   0.98050, -0.08803,  -0.08789,   0.07130]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_097um_COEF_WATER= [   1.00027, -0.35981,   0.15272,  -0.00630]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_104um_COEF_WATER= [   1.00762, -0.38443,   0.24960,  -0.06731]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_133um_COEF_WATER= [   1.03152, 1.31866,  -0.13426,  -0.02558]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_136um_COEF_WATER= [   1.00698, 0.90698,  -0.24137,   0.06605]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_139um_COEF_WATER= [   1.01187, 0.96595,  -0.27568,   0.08154]
REAL, DIMENSION(0:3), PUBLIC:: BETA_110um_142um_COEF_WATER= [   1.01635, 1.01761,  -0.30577,   0.09469]

REAL, DIMENSION(0:3), PUBLIC:: RE_BETA_104um_COEF_WATER= [   0.60261,   0.91156, -0.42765,   0.12638]
REAL, DIMENSION(0:2), PUBLIC:: Qe_104um_COEF_WATER= [  -1.74690,   5.50148, -1.81816]
REAL, DIMENSION(0:2), PUBLIC:: wo_104um_COEF_WATER= [   0.07833,   0.79327, -0.30628]
REAL, DIMENSION(0:2), PUBLIC:: g_104um_COEF_WATER= [   0.33663,   0.83189, -0.26324]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_038um_COEF_WATER= [   0.65417, -0.64834,   1.59679,  -0.59046]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_039um_COEF_WATER= [   0.67989, -0.67382,   1.65045,  -0.63574]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_062um_COEF_WATER= [   0.96281, 0.84226,   0.01369,  -0.05356]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_067um_COEF_WATER= [   0.95935, 0.48977,  -0.10262,  -0.01005]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_073um_COEF_WATER= [   0.96142, 0.38785,  -0.13237,   0.01127]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_085um_COEF_WATER= [   0.97448, 0.21044,  -0.15579,   0.04407]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_097um_COEF_WATER= [   0.99286, 0.01410,  -0.04529,   0.01918]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_133um_COEF_WATER= [   1.03318, 1.23599,  -0.03731,  -0.01983]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_136um_COEF_WATER= [   1.00639, 0.93408,  -0.13728,   0.02116]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_139um_COEF_WATER= [   1.01160, 0.97729,  -0.14879,   0.02480]
REAL, DIMENSION(0:3), PUBLIC:: BETA_104um_142um_COEF_WATER= [   1.01636, 1.01519,  -0.15888,   0.02782]


contains

  subroutine SETUP_ICE_MICROPHYSICAL_MODEL(WMO_Id)

   integer, intent(in):: WMO_Id

   !--- water microphysical model derived from Mie theory by A. Heidinger
   !-- ice clouds modeled as aggregate_columns (b=0.10)
   select case(WMO_Id)

     case(3:5,200:209,223,706:708)  ! AVHRR
        include 'include/acha_ice_cloud_microphysical_model_avhrr_104um.inc'
        include 'include/acha_ice_cloud_microphysical_model_avhrr_110um.inc'

     case(173:174)  ! AHI
        include 'include/acha_ice_cloud_microphysical_model_ahi_104um.inc'
        include 'include/acha_ice_cloud_microphysical_model_ahi_110um.inc'

     case(270:271)  !  ABI
        include 'include/acha_ice_cloud_microphysical_model_abi_104um.inc'
        include 'include/acha_ice_cloud_microphysical_model_abi_110um.inc'

     case(224,225) ! VIIRS 
        include 'include/acha_ice_cloud_microphysical_model_viirs_104um.inc'
        include 'include/acha_ice_cloud_microphysical_model_viirs_110um.inc'

     case(252:259) ! GOES
        include 'include/acha_ice_cloud_microphysical_model_goes_104um.inc'
        include 'include/acha_ice_cloud_microphysical_model_goes_110um.inc'

     case(55:57,70) ! SEVIRI
        include 'include/acha_ice_cloud_microphysical_model_seviri_104um.inc'
        include 'include/acha_ice_cloud_microphysical_model_seviri_110um.inc'

     case default  ! MODIS
        include 'include/acha_ice_cloud_microphysical_model_modis_104um.inc'
        include 'include/acha_ice_cloud_microphysical_model_modis_110um.inc'

   end select

   end subroutine SETUP_ICE_MICROPHYSICAL_MODEL

end module ACHA_MICROPHYSICAL_MODULE
