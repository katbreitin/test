 !$Id: nb_cloud_mask_clavrx_services_module.f90 4089 2021-03-02 18:57:07Z dbotambekov $
 !-----------------------------------------------------------------------------
 ! Input Structure
 !-----------------------------------------------------------------------------
 MODULE NB_CLOUD_MASK_SERVICES

 use CONSTANTS_MOD, only: int1, int2, int4, int8, REAL4, REAL8,  &
                          MISSING_VALUE_REAL4, MISSING_VALUE_INT1, &
                          MISSING_VALUE_INT4, dtor, pi
 
 implicit none

 
 include 'nb_cloud_mask_sevices.inc'

 TYPE, public :: mask_input
    INTEGER :: Num_Elem = Missing_Value_Int4                         !x-DIMENSION of data arrays
    INTEGER :: Num_Line = Missing_Value_Int4                         !number of lines of arrays with data
    INTEGER :: Num_Line_Max = Missing_Value_Int4                     !y-DIMENSION of data arrays
    INTEGER(kind=int1) :: Invalid_Data_Mask = Missing_Value_Int1     !bad data mask (0=good,1=bad)
    INTEGER :: Chan_On_041um = Missing_Value_Int4                    !flag if 0.41um channel on (0=no,1=yes)
    INTEGER :: Chan_On_063um = Missing_Value_Int4                    !flag if 0.63um channel on (0=no,1=yes)
    INTEGER :: Chan_On_086um = Missing_Value_Int4                    !flag if 0.86um channel on (0=no,1=yes)
    INTEGER :: Chan_On_138um = Missing_Value_Int4                    !flag if 1.38um channel on (0=no,1=yes)
    INTEGER :: Chan_On_160um = Missing_Value_Int4                    !flag if 1.60um channel on (0=no,1=yes)
    INTEGER :: Chan_On_213um = Missing_Value_Int4                    !flag if 2.13um channel on (0=no,1=yes)
    INTEGER :: Chan_On_375um = Missing_Value_Int4                    !flag if 3.75um channel on (0=no,1=yes)
    INTEGER :: Chan_On_62um = Missing_Value_Int4                     !flag if 6.2um channel on (0=no,1=yes)
    INTEGER :: Chan_On_67um = Missing_Value_Int4                     !flag if 6.7um channel on (0=no,1=yes)
    INTEGER :: Chan_On_73um = Missing_Value_Int4                     !flag if 7.3um channel on (0=no,1=yes)
    INTEGER :: Chan_On_85um = Missing_Value_Int4                     !flag if 8.5um channel on (0=no,1=yes)
    INTEGER :: Chan_On_97um = Missing_Value_Int4                     !flag if 9.7um channel on (0=no,1=yes)
    INTEGER :: Chan_On_10um = Missing_Value_Int4                     !flag if 10.0um channel on (0=no,1=yes)
    INTEGER :: Chan_On_11um = Missing_Value_Int4                     !flag if 11.0um channel on (0=no,1=yes)
    INTEGER :: Chan_On_12um = Missing_Value_Int4                     !flag if 12.0um channel on (0=no,1=yes)
    INTEGER :: Chan_On_133um = Missing_Value_Int4                    !flag if 13.3um channel on (0=no,1=yes)
    INTEGER :: Chan_On_I1_064um = Missing_Value_Int4                 !flag if I1 0.64um channel on (0=no,1=yes)
    INTEGER :: Chan_On_I4_374um = Missing_Value_Int4                 !flag if I4 3.74um channel on (0=no,1=yes)
    INTEGER :: Chan_On_I5_114um = Missing_Value_Int4                 !flag if I5 11.4um channel on (0=no,1=yes)
    INTEGER :: Chan_On_DNB = Missing_Value_Int4                      !flag if DNB channel on (0=no,1=yes)
    INTEGER :: Use_Sounder_11um = Missing_Value_Int4                 !flag for IFF files where both imager and sounder 11um are available    
    INTEGER(kind=int1) :: Snow_Class = Missing_Value_Int1            !Snow Classification 
    INTEGER(kind=int1) :: Land_Class = Missing_Value_Int1            !Land Classification
    INTEGER(kind=int1) :: Oceanic_Glint_Mask = Missing_Value_Int1    !Mask of oceanic solar glint (0=no,1=yes)
    INTEGER(kind=int1) :: Lunar_Oceanic_Glint_Mask = Missing_Value_Int1 !Mask of oceanic lunar glint (0=no,1=yes)
    INTEGER(kind=int1) :: Coastal_Mask = Missing_Value_Int1          !binary coast mask (0=no,1=yes)
    INTEGER(kind=int1) :: Scat_Mask = Missing_Value_Int1             !binary scattering mask (0=no,1=yes)
    REAL(kind=REAL4) :: Glintzen = Missing_Value_Real4               !Solar glint zenith angle (degrees)
    REAL(kind=REAL4) :: Solzen = Missing_Value_Real4                 !Solar zenith angle (degrees)
    REAL(kind=REAL4) :: Scatzen = Missing_Value_Real4                !Solar Scattering angle (degrees)
    REAL(kind=REAL4) :: Lunscatzen = Missing_Value_Real4             !Lunar Scattering angle (degrees)
    REAL(kind=REAL4) :: Senzen = Missing_Value_Real4                 !Sensor viewing zenith angle (degrees)
    REAL(kind=REAL4) :: Lunzen = Missing_Value_Real4                 !Lunar viewing zenith angle (degrees)
    REAL(kind=REAL4) :: Lunglintzen = Missing_Value_Real4            !Lunar glint zenith angle (degrees)
    REAL(kind=REAL4) :: Lat = Missing_Value_Real4                    !Latitude (degrees)
    REAL(kind=REAL4) :: Lon = Missing_Value_Real4                    !Longitude (degrees)
    REAL(kind=REAL4) :: Ref_041um = Missing_Value_Real4              !0.41 um toa reflectance (%)
    REAL(kind=REAL4) :: Ref_047um = Missing_Value_Real4              !0.47 um toa reflectance (%)
    REAL(kind=REAL4) :: Ref_063um = Missing_Value_Real4              !0.63 um toa reflectance (%)
    REAL(kind=REAL4) :: Ref_063um_Clear = Missing_Value_Real4        !0.63 um toa reflectance for clear-sky (%)
    REAL(kind=REAL4) :: Ref_063um_Std = Missing_Value_Real4          !0.63 um toa reflectance 3x3 Std.  Dev. (%)
    REAL(kind=REAL4) :: Log_Ref_063um_Std = Missing_Value_Real4      !Log10 of 0.63 um toa reflectance 3x3 Std.  Dev.
    REAL(kind=REAL4) :: Ref_063um_Min = Missing_Value_Real4          !Min 0.63 um toa reflectance over 3x3 (%)
    REAL(kind=REAL4) :: Ref_086um  = Missing_Value_Real4             !0.86 um toa reflectance (%)
    REAL(kind=REAL4) :: Ref_138um = Missing_Value_Real4              !1.38 um toa reflectance (%)
    REAL(kind=REAL4) :: Ref_160um = Missing_Value_Real4              !1.60 um toa reflectance (%)
    REAL(kind=REAL4) :: Ref_160um_Clear = Missing_Value_Real4        !1.60 um toa reflectance for clear-sky (%)
    REAL(kind=REAL4) :: Ref_213um = Missing_Value_Real4              !2.13 um toa reflectance (%)
    REAL(kind=REAL4) :: Ref_375um = Missing_Value_Real4              !3.75 um toa reflectance (%)
    REAL(kind=REAL4) :: Ref_375um_Clear = Missing_Value_Real4        !3.75 um toa reflectance for clear-sky (%)
    REAL(kind=REAL4) :: Bt_375um = Missing_Value_Real4               !3.75 um toa brightness temp (K)
    REAL(kind=REAL4) :: Bt_375um_Clear = Missing_Value_Real4         !3.75 um toa brightness temp (K)
    REAL(kind=REAL4) :: Bt_375um_Std = Missing_Value_Real4           !3.75 um toa brightness temp 3x3 Std. Dev. (K)
    REAL(kind=REAL4) :: Log_Bt_375um_Std = Missing_Value_Real4       !Log10 3.75 um toa brightness temp 3x3 Std. Dev. (K)
    REAL(kind=REAL4) :: Emiss_375um = Missing_Value_Real4            !3.75 um pseudo toa emissivity
    REAL(kind=REAL4) :: Emiss_375um_Clear = Missing_Value_Real4      !3.75 um pseudo toa emissivity clear-sky
    REAL(kind=REAL4) :: Bt_62um = Missing_Value_Real4                !6.2 um toa brightness temperature (K)
    REAL(kind=REAL4) :: Bt_67um = Missing_Value_Real4                !6.7 um toa brightness temperature (K)
    REAL(kind=REAL4) :: Bt_73um = Missing_Value_Real4                !7.3 um toa brightness temperature (K)
    REAL(kind=REAL4) :: Bt_85um = Missing_Value_Real4                !8.5 um toa brightness temperature (K)
    REAL(kind=REAL4) :: Bt_10um = Missing_Value_Real4                !10 um toa brightness temperature (K)
    REAL(kind=REAL4) :: Bt_11um = Missing_Value_Real4                !11 um toa brightness temperature (K)
    REAL(kind=REAL4) :: Bt_11um_Sounder = Missing_Value_Real4        !11 um toa brightness temp from sounder (K)
    REAL(kind=REAL4) :: Bt_10um_Std = Missing_Value_Real4            !10.4 um toa brightness temp 3x3 Std Dev (K)
    REAL(kind=REAL4) :: Log_Bt_10um_Std = Missing_Value_Real4        !Log10 10.3 um toa brightness temp 3x3 Std Dev (K)
    REAL(kind=REAL4) :: Bt_11um_Std = Missing_Value_Real4            !11 um toa brightness temp 3x3 Std Dev (K)
    REAL(kind=REAL4) :: Log_Bt_11um_Std = Missing_Value_Real4        !Log10 11 um toa brightness temp 3x3 Std Dev (K)
    REAL(kind=REAL4) :: Bt_10um_Max = Missing_Value_Real4            !10.4 um toa brightness temp 3x3 Max (K)
    REAL(kind=REAL4) :: Bt_11um_Max = Missing_Value_Real4            !11 um toa brightness temp 3x3 Max (K)
    REAL(kind=REAL4) :: Bt_10um_Clear = Missing_Value_Real4          !10.4 um toa brightness temperature (K)
    REAL(kind=REAL4) :: Bt_11um_Clear = Missing_Value_Real4          !11 um toa brightness temperature (K)
    REAL(kind=REAL4) :: Emiss_10um_Tropo = Missing_Value_Real4       !10.4 um tropo emiss
    REAL(kind=REAL4) :: Emiss_11um_Tropo = Missing_Value_Real4       !11 um tropo emiss
    REAL(kind=REAL4) :: Bt_12um = Missing_Value_Real4                !12 um toa brightness temperature (K)
    REAL(kind=REAL4) :: Bt_12um_Clear = Missing_Value_Real4          !12 um toa bright temp clear-sky (K)
    REAL(kind=REAL4) :: Bt_10um_Bt_67um_Covar = Missing_Value_Real4  !covariance of 10.4 and 6.7 um bright temp.
    REAL(kind=REAL4) :: Bt_11um_Bt_67um_Covar = Missing_Value_Real4  !covariance of 11 and 6.7 um bright temp.
    REAL(kind=REAL4) :: Bt_133um = Missing_Value_Real4               !13.3 um toa brightness temperature (K)
    REAL(kind=REAL4) :: Emiss_Sfc_375um = Missing_Value_Real4        !the surface emissivity at 3.75 um
    REAL(kind=REAL4) :: Rad_Lunar = Missing_Value_Real4              !Lunar toa radiance from DNB
    REAL(kind=REAL4) :: Ref_Lunar = Missing_Value_Real4              !Lunar reflectance from DNB (%)
    REAL(kind=REAL4) :: Ref_Lunar_Min = Missing_Value_Real4          !Min lunar reflectance over 3x3 (%)
    REAL(kind=REAL4) :: Ref_Lunar_Std = Missing_Value_Real4          !3x3 std dev of lunar ref from DNB (%)
    REAL(kind=REAL4) :: Ref_Lunar_Clear = Missing_Value_Real4        !Lunar reflectance for clear-skies (%)
    REAL(kind=REAL4) :: Cld_Opd_DNB = Missing_Value_Real4            !Cloud cod from lunar DNB
    REAL(kind=REAL4) :: Log_Cld_Opd_DNB = Missing_Value_Real4        !Log of cloud cod from lunar DNB
    REAL(kind=REAL4) :: Zsfc = Missing_Value_Real4                   !surface altitude (km)
    REAL(kind=REAL4) :: Zsfc_Std = Missing_Value_Real4               !surface altitude 3x3 Std Dev (km)
    INTEGER :: Num_Segments = Missing_Value_Int4                     !number of segments in this data 
    INTEGER(kind=int1) :: Solar_Contamination_Mask = Missing_Value_Int1 !binary mask of solar contamination (0=no,1=yes)
    INTEGER(kind=int1) :: Sfc_Type = Missing_Value_Int1              !surface TYPE based on UMD classification
    REAL(kind=REAL4) :: Sst_Anal_Uni = Missing_Value_Real4           !surface temperature from ancillary sources
    REAL(kind=REAL4) :: Sfc_Temp = Missing_Value_Real4               !surface temperature from ancillary sources
    REAL(kind=REAL4) :: Path_Tpw = Missing_Value_Real4               !TPW along IR path from ancillary sources
    REAL(kind=REAL4) :: Prior = Missing_Value_Real4                  !Prior from a precomputed source
    real(kind=real4) :: Cld_Fraction_Background
    REAL(kind=REAL4) :: Topa = Missing_Value_Real4                   !NWP temperature of opaque cloud
    REAL(kind=REAL4) :: Zopa = Missing_Value_Real4                   !NWP height of opaque cloud
    REAL(kind=REAL4) :: Ttropo = Missing_Value_Real4                 !NWP temperature of tropopause
    REAL(kind=REAL4) :: Bt_11um_Min_Sub = Missing_Value_Real4        !11 um toa brightness temp subpixel Min (K)
    REAL(kind=REAL4) :: Bt_11um_Max_Sub = Missing_Value_Real4        !11 um toa brightness temp subpixel Max (K)
    REAL(kind=REAL4) :: Ref_063um_Min_Sub = Missing_Value_Real4      !0.63 um toa reflectance subpixel Min (%)
    REAL(kind=REAL4) :: Ref_063um_Max_Sub = Missing_Value_Real4      !0.63 um toa reflectance subpixel Max (%)
    REAL(kind=REAL4) :: Ref_063um_Std_Sub = Missing_Value_Real4      !0.63 um toa reflectance subpixel Std (%)
    REAL(kind=REAL4) :: Log_Drefl_065um_Max_Min_Sub = Missing_Value_Real4 !0.63 um toa reflectance subpixel difference Max-Min (%)
    REAL(kind=REAL4) :: Moon_Illum_Frac = Missing_Value_Real4        !moon illumination fraction for dnb
    INTEGER(kind=int1) :: City_Mask = Missing_Value_Int1             !city lights mask for dnb
    INTEGER(kind=int1) :: Use_Aux_Mask = Missing_Value_Int1          !aux mask flag

 END TYPE mask_input 

 !-----------------------------------------------------------------------------
 ! Output Structure
 !-----------------------------------------------------------------------------
 TYPE, public :: mask_output
    INTEGER(kind=int1), DIMENSION(NUMBER_OF_FLAG_BYTES) :: Cld_Flags_Packed !array of packed results 
    INTEGER(kind=int1) :: Cld_Mask_Bayes                !Derived 4-level cloud mask
    INTEGER(kind=int1) :: Cld_Mask_Binary               !Derived 2-level cloud mask
    INTEGER(kind=int1) :: Cld_Mask_Bayes_IR             !Derived 4-level cloud mask
    INTEGER(kind=int1) :: Cld_Mask_Binary_IR            !Derived 2-level cloud mask
    INTEGER :: Cloud_Mask_Bayesian_Flag                 !flag to tell if code should run
    REAL(kind=REAL4) :: Prior_Cld_Probability           !prior cloud probability (0-1)
    REAL(kind=REAL4) :: Posterior_Cld_Probability       !posterior cloud probability (0-1)
    REAL(kind=REAL4) :: Posterior_Cld_Probability_IR    !posterior cloud probability (0-1)
    REAL(kind=REAL4) :: Posterior_Ice_Probability       !posterior ice probability (0-1)
    REAL(kind=REAL4) :: Posterior_Water_Probability     !posterior water probability (0-1)
    REAL(kind=REAL4) :: Cloud_Phase                     !Cloud Phase
    REAL(kind=REAL4) :: Cloud_Phase_Uncer               !Cloud Phase uncertainty
    INTEGER(kind=int1) :: Dust_Mask 
    INTEGER(kind=int1) :: Smoke_Mask 
    INTEGER(kind=int1) :: Fire_Mask 
    INTEGER(kind=int1) :: Thin_Cirr_Mask
    INTEGER(kind=int1) :: Cld_Mask_QF                   ! Quality Flag 0=good, 1=bad
    INTEGER(kind=int1) :: TUT 
    INTEGER(kind=int1) :: RUT
    INTEGER(kind=int1) :: Sfc_Idx      
    CHARACTER(2020) :: Cld_Mask_Test_Names              ! Cloud tests string
 END TYPE mask_output

 !-----------------------------------------------------------------------------
 ! Diagnostic Output Structure
 !-----------------------------------------------------------------------------
 TYPE, public :: diag_output
    REAL(kind=REAL4) :: Array_1    !first diagnostic array
    REAL(kind=REAL4) :: Array_2    !first diagnostic array
    REAL(kind=REAL4) :: Array_3    !first diagnostic array
 END TYPE diag_output 
 
 !-----------------------------------------------------------------------------
 ! Symbol Structure
 !-----------------------------------------------------------------------------
 TYPE, public :: symbol_naive_bayesian
    INTEGER(kind=int1) :: CLOUDY
    INTEGER(kind=int1) :: PROB_CLOUDY
    INTEGER(kind=int1) :: PROB_CLEAR
    INTEGER(kind=int1) :: CLEAR

    INTEGER(kind=int1) :: GOOD_QF
    INTEGER(kind=int1) :: BAD_QF
    INTEGER(kind=int1) :: SPACE_QF
    INTEGER(kind=int1) :: FILL_QF
    INTEGER(kind=int1) :: DEGRADED_QF

    INTEGER(kind=int1) :: CLOUDY_BINARY
    INTEGER(kind=int1) :: CLEAR_BINARY

    INTEGER(kind=int1) :: NO
    INTEGER(kind=int1) :: YES

    INTEGER(kind=int1) :: WATER_SFC
    INTEGER(kind=int1) :: EVERGREEN_NEEDLE_SFC
    INTEGER(kind=int1) :: EVERGREEN_BROAD_SFC
    INTEGER(kind=int1) :: DECIDUOUS_NEEDLE_SFC
    INTEGER(kind=int1) :: DECIDUOUS_BROAD_SFC
    INTEGER(kind=int1) :: MIXED_FORESTS_SFC
    INTEGER(kind=int1) :: WOODLANDS_SFC
    INTEGER(kind=int1) :: WOODED_GRASS_SFC
    INTEGER(kind=int1) :: CLOSED_SHRUBS_SFC
    INTEGER(kind=int1) :: OPEN_SHRUBS_SFC
    INTEGER(kind=int1) :: GRASSES_SFC
    INTEGER(kind=int1) :: CROPLANDS_SFC
    INTEGER(kind=int1) :: BARE_SFC
    INTEGER(kind=int1) :: URBAN_SFC

    INTEGER(kind=int1) :: SHALLOW_OCEAN
    INTEGER(kind=int1) :: LAND
    INTEGER(kind=int1) :: COASTLINE
    INTEGER(kind=int1) :: SHALLOW_INLAND_WATER
    INTEGER(kind=int1) :: EPHEMERAL_WATER
    INTEGER(kind=int1) :: DEEP_INLAND_WATER
    INTEGER(kind=int1) :: MODERATE_OCEAN
    INTEGER(kind=int1) :: DEEP_OCEAN

    INTEGER(kind=int1) :: NO_SNOW
    INTEGER(kind=int1) :: SEA_ICE
    INTEGER(kind=int1) :: SNOW

    INTEGER(kind=int1) :: NO_AUX
    INTEGER(kind=int1) :: USE_AUX_MODAWG

 END TYPE symbol_naive_bayesian
 !-----------------------------------------------------------------------------
 ! Classifier Structure
 !-----------------------------------------------------------------------------
 TYPE, public :: Classifier
   !--- lookup tables
   REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Clear_Table
   REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Water_Table
   REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Ice_Table
   REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Obs_Table
   REAL, DIMENSION(:), ALLOCATABLE :: Observation_Count
   REAL, DIMENSION(:), ALLOCATABLE :: Cloud_Fraction
   REAL, DIMENSION(:), ALLOCATABLE :: Ice_Fraction
   REAL, DIMENSION(:), ALLOCATABLE :: Water_Fraction
   INTEGER, DIMENSION(:), ALLOCATABLE :: Wvl
   INTEGER, DIMENSION(:), ALLOCATABLE :: On_Flag

   !--- attributes
   INTEGER :: Rank
   INTEGER :: Nchan_Used
   INTEGER :: N_Sfc
   INTEGER :: Nbins_X, Nbins_Y, Nbins_Z
   REAL :: X_Min, Y_Min, Z_Min
   REAL :: X_Bin, Y_Bin, Z_Bin
   REAL :: Zen_Min, Zen_Max
   REAL :: Solzen_Min, Solzen_Max
   REAL :: Solglintzen_Min, Solglintzen_Max
   REAL :: Solglint_Mask_Min, Solglint_Mask_Max
   REAL :: Lunzen_Min, Lunzen_Max
   REAL :: Lunglintzen_Min, Lunglintzen_Max
   REAL :: Lunglint_Mask_Min, Lunglint_Mask_Max
   REAL :: Solscatang_Min, Solscatang_Max
   REAL :: Tsfc_Min, Tsfc_Max
   REAL :: Zsfc_Min, Zsfc_Max
   REAL :: Zsfc_Std_Min, Zsfc_Std_Max
   REAL :: Tpw_Min, Tpw_Max
   REAL :: Snow_Class_Min, Snow_Class_Max
   REAL :: Coast_Mask_Min, Coast_Mask_Max
   REAL :: City_Mask_Min, City_Mask_Max
   REAL :: Moon_Illum_Frac_Min, Moon_Illum_Frac_Max
   REAL :: Rut_Solzen_Thresh
   CHARACTER(len=:), ALLOCATABLE :: Class_Xname, Class_Yname, Class_Zname
 END TYPE Classifier

 CHARACTER(len=:), DIMENSION(:), ALLOCATABLE :: Classifier_Names

 !-----------------------------------------------------------------------------
 ! LUT Clear/ProbClear/Prob Cloud Threhosld Structure
 !-----------------------------------------------------------------------------
 TYPE, public :: Mask_Threshold
   REAL, DIMENSION(:), ALLOCATABLE :: Conf_Clear_Prob_Clear_Thresh
   REAL, DIMENSION(:), ALLOCATABLE :: Prob_Clear_Prob_Cloudy_Thresh
   REAL, DIMENSION(:), ALLOCATABLE :: Prob_Cloudy_Conf_Cloudy_Thresh
   REAL, DIMENSION(:), ALLOCATABLE :: Rut_Clear_Prob_Clear_Thresh
   REAL, DIMENSION(:), ALLOCATABLE :: Tut_Clear_Prob_Clear_Thresh
 END TYPE Mask_Threshold


END MODULE NB_CLOUD_MASK_SERVICES
