!$Id: acha_module.f90 3876 2020-06-18 13:34:40Z yli $
module ACHA_RTM_MOD
!---------------------------------------------------------------------
!
!----------------------------------------------------------------------

  use ACHA_SERVICES_MOD, only : &
           real4, int1, int4, real8, dtor, &
           PLANCK_RAD_FAST, PLANCK_TEMP_FAST, &
           ACHA_SYMBOL_STRUCT, &
           INVERT_MATRIX, LOCATE

  use ACHA_MICROPHYSICAL_MODULE, only: beta_degree_ice, beta_degree_water, &
                                       Qe_006um_COEF_ICE,  Qe_110um_COEF_ICE, &
                                       Re_Beta_110um_COEF_ICE

  use ACHA_NUM_MOD, only:GENERIC_PROFILE_INTERPOLATION &
    ,Chan_Idx_y &
    ,knowing_t_compute_p_z &
    ,knowing_z_compute_t_p &
    , OPTIMAL_ESTIMATION

  implicit none

  public:: COMPUTE_CLEAR_SKY_TERMS
  public:: SET_CLEAR_SKY_COVARIANCE_TERMS
  private:: CLEAR_SKY_INTERNAL_ROUTINE

  public:: BT_FM
  public:: BTD_FM
  private:: COMPUTE_BETA_AND_DERIVATIVE

  public:: COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE
  public:: DETERMINE_ACHA_EXTINCTION
  public:: DETERMINE_ACHA_ICE_EXTINCTION

  !--- include the non-system specific variables
  include 'acha_parameters.inc'

  real, dimension(20:38), public, save:: Bt_Covar           !
  real, dimension(20:38,20:38), public, save:: Btd_Covar   !
  real, dimension(20:38), public, save:: Cal_Uncer         !
  real, dimension(20:38), public, save:: Cloud_BTD_Uncer   !
  real, public, save:: Cloud_BT_Uncer                      !

  real, private, PARAMETER:: MISSING_VALUE_REAL4 = -999.0
  integer(kind=int1), private, PARAMETER:: MISSING_VALUE_integer1 = -128_int1
  integer(kind=int4), private, PARAMETER:: MISSING_VALUE_integer4 = -999

  !-------------------------------------------------------------------------------
  !ice extinction
  !-------------------------------------------------------------------------------
  integer, parameter, private:: Nec_Ext = 10
  real, parameter, private:: Ec_Ext_Min = 0.0
  real, parameter, private:: Ec_Ext_Bin = 0.1
  integer, parameter, private:: M_Ext = 4
  real, parameter, dimension(Nec_Ext), private:: Ec_Ext = &
                               [0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]
  real, parameter, private:: Tc_Ext_Offset = 182.5
  real, dimension(Nec_Ext,4),save:: Ice_Ext_Coef

  contains 

!----------------------------------------------------------------------
! Compute Sy based on the clear-sky error covariance calcuLations.
! This assumes that 
!----------------------------------------------------------------------
subroutine SET_CLEAR_SKY_COVARIANCE_TERMS(Sfc_Type_Forward_Model)

 integer(kind=int4), intent(in):: Sfc_Type_Forward_Model

 Bt_Covar = 0.00
 Btd_Covar = 0.00

 !--- The below values are parameters in acha_parameters.inc.
 !--- Use for 104 um for now.
 !--- calibration uncertainties  (not surface dependant)
 Cal_Uncer(20) = T110um_038um_Cal_Uncer 
 Cal_Uncer(37) = T110um_062um_Cal_Uncer 
 Cal_Uncer(27) = T110um_067um_Cal_Uncer 
 Cal_Uncer(28) = T110um_073um_Cal_Uncer 
 Cal_Uncer(29) = T110um_085um_Cal_Uncer 
 Cal_Uncer(30) = T110um_097um_Cal_Uncer 
 Cal_Uncer(38) = T110um_104um_Cal_Uncer 
 Cal_Uncer(31) = T110um_Cal_Uncer       !Note, not a BTD
 Cal_Uncer(32) = T110um_120um_Cal_Uncer 
 Cal_Uncer(33) = T110um_133um_Cal_Uncer 
 Cal_Uncer(34) = T110um_136um_Cal_Uncer 
 Cal_Uncer(35) = T110um_139um_Cal_Uncer 
 Cal_Uncer(36) = T110um_142um_Cal_Uncer 

 !--- additional terms to Sy for the cloud error (Bt and BTD)
 Cloud_BTD_Uncer = 1.0!2.0
 Cloud_BT_Uncer = 4.0!5.0

 !--- All values below are parameters from acha_clear_sky_covariances.inc.
 !--- Use for 104 um for now.
 select case(Sfc_Type_Forward_Model)

   !--- Water
   case (0) 

   Bt_Covar(31) = Bt_110um_Bt_110um_Covar_Water            !Note, not a BTD
   Bt_Covar(20) = Bt_110um_Btd_110um_038um_Covar_Water
   Bt_Covar(37) = Bt_110um_Btd_110um_062um_Covar_Water
   Bt_Covar(27) = Bt_110um_Btd_110um_067um_Covar_Water
   Bt_Covar(28) = Bt_110um_Btd_110um_073um_Covar_Water
   Bt_Covar(29) = Bt_110um_Btd_110um_085um_Covar_Water
   Bt_Covar(30) = Bt_110um_Btd_110um_097um_Covar_Water
   Bt_Covar(38) = Bt_110um_Btd_110um_104um_Covar_Water
   Bt_Covar(32) = Bt_110um_Btd_110um_120um_Covar_Water
   Bt_Covar(33) = Bt_110um_Btd_110um_133um_Covar_Water
   Bt_Covar(34) = Bt_110um_Btd_110um_136um_Covar_Water
   Bt_Covar(35) = Bt_110um_Btd_110um_139um_Covar_Water
   Bt_Covar(36) = Bt_110um_Btd_110um_142um_Covar_Water

   Btd_Covar(20,20) = Btd_110um_038um_Btd_110um_038um_Covar_Water
   Btd_Covar(20,37) = Btd_110um_038um_Btd_110um_062um_Covar_Water
   Btd_Covar(20,27) = Btd_110um_038um_Btd_110um_067um_Covar_Water
   Btd_Covar(20,28) = Btd_110um_038um_Btd_110um_073um_Covar_Water
   Btd_Covar(20,29) = Btd_110um_038um_Btd_110um_085um_Covar_Water
   Btd_Covar(20,30) = Btd_110um_038um_Btd_110um_097um_Covar_Water
   Btd_Covar(20,38) = Btd_110um_038um_Btd_110um_104um_Covar_Water
   Btd_Covar(20,32) = Btd_110um_038um_Btd_110um_120um_Covar_Water
   Btd_Covar(20,33) = Btd_110um_038um_Btd_110um_133um_Covar_Water
   Btd_Covar(20,34) = Btd_110um_038um_Btd_110um_136um_Covar_Water
   Btd_Covar(20,35) = Btd_110um_038um_Btd_110um_139um_Covar_Water
   Btd_Covar(20,36) = Btd_110um_038um_Btd_110um_142um_Covar_Water

   Btd_Covar(37,20) = Btd_110um_062um_Btd_110um_038um_Covar_Water
   Btd_Covar(37,37) = Btd_110um_062um_Btd_110um_062um_Covar_Water
   Btd_Covar(37,27) = Btd_110um_062um_Btd_110um_067um_Covar_Water
   Btd_Covar(37,28) = Btd_110um_062um_Btd_110um_073um_Covar_Water
   Btd_Covar(37,29) = Btd_110um_062um_Btd_110um_085um_Covar_Water
   Btd_Covar(37,30) = Btd_110um_062um_Btd_110um_097um_Covar_Water
   Btd_Covar(37,38) = Btd_110um_062um_Btd_110um_104um_Covar_Water
   Btd_Covar(37,32) = Btd_110um_062um_Btd_110um_120um_Covar_Water
   Btd_Covar(37,33) = Btd_110um_062um_Btd_110um_133um_Covar_Water
   Btd_Covar(37,34) = Btd_110um_062um_Btd_110um_136um_Covar_Water
   Btd_Covar(37,35) = Btd_110um_062um_Btd_110um_139um_Covar_Water
   Btd_Covar(37,36) = Btd_110um_062um_Btd_110um_142um_Covar_Water

   Btd_Covar(27,20) = Btd_110um_067um_Btd_110um_038um_Covar_Water
   Btd_Covar(27,37) = Btd_110um_067um_Btd_110um_062um_Covar_Water
   Btd_Covar(27,27) = Btd_110um_067um_Btd_110um_067um_Covar_Water
   Btd_Covar(27,28) = Btd_110um_067um_Btd_110um_073um_Covar_Water
   Btd_Covar(27,29) = Btd_110um_067um_Btd_110um_085um_Covar_Water
   Btd_Covar(27,30) = Btd_110um_067um_Btd_110um_097um_Covar_Water
   Btd_Covar(27,38) = Btd_110um_067um_Btd_110um_104um_Covar_Water
   Btd_Covar(27,32) = Btd_110um_067um_Btd_110um_120um_Covar_Water
   Btd_Covar(27,33) = Btd_110um_067um_Btd_110um_133um_Covar_Water
   Btd_Covar(27,34) = Btd_110um_067um_Btd_110um_136um_Covar_Water
   Btd_Covar(27,35) = Btd_110um_067um_Btd_110um_139um_Covar_Water
   Btd_Covar(27,36) = Btd_110um_067um_Btd_110um_142um_Covar_Water

   Btd_Covar(28,20) = Btd_110um_073um_Btd_110um_038um_Covar_Water
   Btd_Covar(28,37) = Btd_110um_073um_Btd_110um_062um_Covar_Water
   Btd_Covar(28,27) = Btd_110um_073um_Btd_110um_067um_Covar_Water
   Btd_Covar(28,28) = Btd_110um_073um_Btd_110um_073um_Covar_Water
   Btd_Covar(28,29) = Btd_110um_073um_Btd_110um_085um_Covar_Water
   Btd_Covar(28,30) = Btd_110um_073um_Btd_110um_097um_Covar_Water
   Btd_Covar(28,38) = Btd_110um_073um_Btd_110um_104um_Covar_Water
   Btd_Covar(28,32) = Btd_110um_073um_Btd_110um_120um_Covar_Water
   Btd_Covar(28,33) = Btd_110um_073um_Btd_110um_133um_Covar_Water
   Btd_Covar(28,34) = Btd_110um_073um_Btd_110um_136um_Covar_Water
   Btd_Covar(28,35) = Btd_110um_073um_Btd_110um_139um_Covar_Water
   Btd_Covar(28,36) = Btd_110um_073um_Btd_110um_142um_Covar_Water

   Btd_Covar(29,20) = Btd_110um_085um_Btd_110um_038um_Covar_Water
   Btd_Covar(29,37) = Btd_110um_085um_Btd_110um_062um_Covar_Water
   Btd_Covar(29,27) = Btd_110um_085um_Btd_110um_067um_Covar_Water
   Btd_Covar(29,28) = Btd_110um_085um_Btd_110um_073um_Covar_Water
   Btd_Covar(29,29) = Btd_110um_085um_Btd_110um_085um_Covar_Water
   Btd_Covar(29,30) = Btd_110um_085um_Btd_110um_097um_Covar_Water
   Btd_Covar(29,38) = Btd_110um_085um_Btd_110um_104um_Covar_Water
   Btd_Covar(29,32) = Btd_110um_085um_Btd_110um_120um_Covar_Water
   Btd_Covar(29,33) = Btd_110um_085um_Btd_110um_133um_Covar_Water
   Btd_Covar(29,34) = Btd_110um_085um_Btd_110um_136um_Covar_Water
   Btd_Covar(29,35) = Btd_110um_085um_Btd_110um_139um_Covar_Water
   Btd_Covar(29,36) = Btd_110um_085um_Btd_110um_142um_Covar_Water

   Btd_Covar(30,20) = Btd_110um_097um_Btd_110um_038um_Covar_Water
   Btd_Covar(30,37) = Btd_110um_097um_Btd_110um_062um_Covar_Water
   Btd_Covar(30,27) = Btd_110um_097um_Btd_110um_067um_Covar_Water
   Btd_Covar(30,28) = Btd_110um_097um_Btd_110um_073um_Covar_Water
   Btd_Covar(30,29) = Btd_110um_097um_Btd_110um_085um_Covar_Water
   Btd_Covar(30,30) = Btd_110um_097um_Btd_110um_097um_Covar_Water
   Btd_Covar(30,38) = Btd_110um_097um_Btd_110um_104um_Covar_Water
   Btd_Covar(30,32) = Btd_110um_097um_Btd_110um_120um_Covar_Water
   Btd_Covar(30,33) = Btd_110um_097um_Btd_110um_133um_Covar_Water
   Btd_Covar(30,34) = Btd_110um_097um_Btd_110um_136um_Covar_Water
   Btd_Covar(30,35) = Btd_110um_097um_Btd_110um_139um_Covar_Water
   Btd_Covar(30,36) = Btd_110um_097um_Btd_110um_142um_Covar_Water

   Btd_Covar(38,20) = Btd_110um_104um_Btd_110um_038um_Covar_Water
   Btd_Covar(38,37) = Btd_110um_104um_Btd_110um_062um_Covar_Water
   Btd_Covar(38,27) = Btd_110um_104um_Btd_110um_067um_Covar_Water
   Btd_Covar(38,28) = Btd_110um_104um_Btd_110um_073um_Covar_Water
   Btd_Covar(38,29) = Btd_110um_104um_Btd_110um_085um_Covar_Water
   Btd_Covar(38,30) = Btd_110um_104um_Btd_110um_097um_Covar_Water
   Btd_Covar(38,38) = Btd_110um_104um_Btd_110um_104um_Covar_Water
   Btd_Covar(38,32) = Btd_110um_104um_Btd_110um_120um_Covar_Water
   Btd_Covar(38,33) = Btd_110um_104um_Btd_110um_133um_Covar_Water
   Btd_Covar(38,34) = Btd_110um_104um_Btd_110um_136um_Covar_Water
   Btd_Covar(38,35) = Btd_110um_104um_Btd_110um_139um_Covar_Water
   Btd_Covar(38,36) = Btd_110um_104um_Btd_110um_142um_Covar_Water

   Btd_Covar(32,20) = Btd_110um_120um_Btd_110um_038um_Covar_Water
   Btd_Covar(32,37) = Btd_110um_120um_Btd_110um_062um_Covar_Water
   Btd_Covar(32,27) = Btd_110um_120um_Btd_110um_067um_Covar_Water
   Btd_Covar(32,28) = Btd_110um_120um_Btd_110um_073um_Covar_Water
   Btd_Covar(32,29) = Btd_110um_120um_Btd_110um_085um_Covar_Water
   Btd_Covar(32,30) = Btd_110um_120um_Btd_110um_097um_Covar_Water
   Btd_Covar(32,38) = Btd_110um_120um_Btd_110um_104um_Covar_Water
   Btd_Covar(32,32) = Btd_110um_120um_Btd_110um_120um_Covar_Water
   Btd_Covar(32,33) = Btd_110um_120um_Btd_110um_133um_Covar_Water
   Btd_Covar(32,34) = Btd_110um_120um_Btd_110um_136um_Covar_Water
   Btd_Covar(32,35) = Btd_110um_120um_Btd_110um_139um_Covar_Water
   Btd_Covar(32,36) = Btd_110um_120um_Btd_110um_142um_Covar_Water

   Btd_Covar(33,20) = Btd_110um_133um_Btd_110um_038um_Covar_Water
   Btd_Covar(33,37) = Btd_110um_133um_Btd_110um_062um_Covar_Water
   Btd_Covar(33,27) = Btd_110um_133um_Btd_110um_067um_Covar_Water
   Btd_Covar(33,28) = Btd_110um_133um_Btd_110um_073um_Covar_Water
   Btd_Covar(33,29) = Btd_110um_133um_Btd_110um_085um_Covar_Water
   Btd_Covar(33,30) = Btd_110um_133um_Btd_110um_097um_Covar_Water
   Btd_Covar(33,38) = Btd_110um_133um_Btd_110um_104um_Covar_Water
   Btd_Covar(33,32) = Btd_110um_133um_Btd_110um_120um_Covar_Water
   Btd_Covar(33,33) = Btd_110um_133um_Btd_110um_133um_Covar_Water
   Btd_Covar(33,34) = Btd_110um_133um_Btd_110um_136um_Covar_Water
   Btd_Covar(33,35) = Btd_110um_133um_Btd_110um_139um_Covar_Water
   Btd_Covar(33,36) = Btd_110um_133um_Btd_110um_142um_Covar_Water

   Btd_Covar(34,20) = Btd_110um_136um_Btd_110um_038um_Covar_Water
   Btd_Covar(34,37) = Btd_110um_136um_Btd_110um_062um_Covar_Water
   Btd_Covar(34,27) = Btd_110um_136um_Btd_110um_067um_Covar_Water
   Btd_Covar(34,28) = Btd_110um_136um_Btd_110um_073um_Covar_Water
   Btd_Covar(34,29) = Btd_110um_136um_Btd_110um_085um_Covar_Water
   Btd_Covar(34,30) = Btd_110um_136um_Btd_110um_097um_Covar_Water
   Btd_Covar(34,38) = Btd_110um_136um_Btd_110um_104um_Covar_Water
   Btd_Covar(34,32) = Btd_110um_136um_Btd_110um_120um_Covar_Water
   Btd_Covar(34,33) = Btd_110um_136um_Btd_110um_133um_Covar_Water
   Btd_Covar(34,34) = Btd_110um_136um_Btd_110um_136um_Covar_Water
   Btd_Covar(34,35) = Btd_110um_136um_Btd_110um_139um_Covar_Water
   Btd_Covar(34,36) = Btd_110um_136um_Btd_110um_142um_Covar_Water

   Btd_Covar(35,20) = Btd_110um_139um_Btd_110um_038um_Covar_Water
   Btd_Covar(35,37) = Btd_110um_139um_Btd_110um_062um_Covar_Water
   Btd_Covar(35,27) = Btd_110um_139um_Btd_110um_067um_Covar_Water
   Btd_Covar(35,28) = Btd_110um_139um_Btd_110um_073um_Covar_Water
   Btd_Covar(35,29) = Btd_110um_139um_Btd_110um_085um_Covar_Water
   Btd_Covar(35,30) = Btd_110um_139um_Btd_110um_097um_Covar_Water
   Btd_Covar(35,38) = Btd_110um_139um_Btd_110um_104um_Covar_Water
   Btd_Covar(35,32) = Btd_110um_139um_Btd_110um_120um_Covar_Water
   Btd_Covar(35,33) = Btd_110um_139um_Btd_110um_133um_Covar_Water
   Btd_Covar(35,34) = Btd_110um_139um_Btd_110um_136um_Covar_Water
   Btd_Covar(35,35) = Btd_110um_139um_Btd_110um_139um_Covar_Water
   Btd_Covar(35,36) = Btd_110um_139um_Btd_110um_142um_Covar_Water

   Btd_Covar(36,20) = Btd_110um_142um_Btd_110um_038um_Covar_Water
   Btd_Covar(36,37) = Btd_110um_142um_Btd_110um_062um_Covar_Water
   Btd_Covar(36,27) = Btd_110um_142um_Btd_110um_067um_Covar_Water
   Btd_Covar(36,28) = Btd_110um_142um_Btd_110um_073um_Covar_Water
   Btd_Covar(36,29) = Btd_110um_142um_Btd_110um_085um_Covar_Water
   Btd_Covar(36,30) = Btd_110um_142um_Btd_110um_097um_Covar_Water
   Btd_Covar(36,38) = Btd_110um_142um_Btd_110um_104um_Covar_Water
   Btd_Covar(36,32) = Btd_110um_142um_Btd_110um_120um_Covar_Water
   Btd_Covar(36,33) = Btd_110um_142um_Btd_110um_133um_Covar_Water
   Btd_Covar(36,34) = Btd_110um_142um_Btd_110um_136um_Covar_Water
   Btd_Covar(36,35) = Btd_110um_142um_Btd_110um_139um_Covar_Water
   Btd_Covar(36,36) = Btd_110um_142um_Btd_110um_142um_Covar_Water

   !--- Land
   case (1) 

   Bt_Covar(31) = Bt_110um_Bt_110um_Covar_Land
   Bt_Covar(20) = Bt_110um_Btd_110um_038um_Covar_Land
   Bt_Covar(37) = Bt_110um_Btd_110um_062um_Covar_Land
   Bt_Covar(27) = Bt_110um_Btd_110um_067um_Covar_Land
   Bt_Covar(28) = Bt_110um_Btd_110um_073um_Covar_Land
   Bt_Covar(29) = Bt_110um_Btd_110um_085um_Covar_Land
   Bt_Covar(30) = Bt_110um_Btd_110um_097um_Covar_Land
   Bt_Covar(38) = Bt_110um_Btd_110um_104um_Covar_Land
   Bt_Covar(32) = Bt_110um_Btd_110um_120um_Covar_Land
   Bt_Covar(33) = Bt_110um_Btd_110um_133um_Covar_Land
   Bt_Covar(34) = Bt_110um_Btd_110um_136um_Covar_Land
   Bt_Covar(35) = Bt_110um_Btd_110um_139um_Covar_Land
   Bt_Covar(36) = Bt_110um_Btd_110um_142um_Covar_Land

   Btd_Covar(20,20) = Btd_110um_038um_Btd_110um_038um_Covar_Land
   Btd_Covar(20,37) = Btd_110um_038um_Btd_110um_062um_Covar_Land
   Btd_Covar(20,27) = Btd_110um_038um_Btd_110um_067um_Covar_Land
   Btd_Covar(20,28) = Btd_110um_038um_Btd_110um_073um_Covar_Land
   Btd_Covar(20,29) = Btd_110um_038um_Btd_110um_085um_Covar_Land
   Btd_Covar(20,30) = Btd_110um_038um_Btd_110um_097um_Covar_Land
   Btd_Covar(20,38) = Btd_110um_038um_Btd_110um_104um_Covar_Land
   Btd_Covar(20,32) = Btd_110um_038um_Btd_110um_120um_Covar_Land
   Btd_Covar(20,33) = Btd_110um_038um_Btd_110um_133um_Covar_Land
   Btd_Covar(20,34) = Btd_110um_038um_Btd_110um_136um_Covar_Land
   Btd_Covar(20,35) = Btd_110um_038um_Btd_110um_139um_Covar_Land
   Btd_Covar(20,36) = Btd_110um_038um_Btd_110um_142um_Covar_Land

   Btd_Covar(37,20) = Btd_110um_062um_Btd_110um_038um_Covar_Land
   Btd_Covar(37,37) = Btd_110um_062um_Btd_110um_062um_Covar_Land
   Btd_Covar(37,27) = Btd_110um_062um_Btd_110um_067um_Covar_Land
   Btd_Covar(37,28) = Btd_110um_062um_Btd_110um_073um_Covar_Land
   Btd_Covar(37,29) = Btd_110um_062um_Btd_110um_085um_Covar_Land
   Btd_Covar(37,30) = Btd_110um_062um_Btd_110um_097um_Covar_Land
   Btd_Covar(37,38) = Btd_110um_062um_Btd_110um_104um_Covar_Land
   Btd_Covar(37,32) = Btd_110um_062um_Btd_110um_120um_Covar_Land
   Btd_Covar(37,33) = Btd_110um_062um_Btd_110um_133um_Covar_Land
   Btd_Covar(37,34) = Btd_110um_062um_Btd_110um_136um_Covar_Land
   Btd_Covar(37,35) = Btd_110um_062um_Btd_110um_139um_Covar_Land
   Btd_Covar(37,36) = Btd_110um_062um_Btd_110um_142um_Covar_Land

   Btd_Covar(27,20) = Btd_110um_067um_Btd_110um_038um_Covar_Land
   Btd_Covar(27,37) = Btd_110um_067um_Btd_110um_062um_Covar_Land
   Btd_Covar(27,27) = Btd_110um_067um_Btd_110um_067um_Covar_Land
   Btd_Covar(27,28) = Btd_110um_067um_Btd_110um_073um_Covar_Land
   Btd_Covar(27,29) = Btd_110um_067um_Btd_110um_085um_Covar_Land
   Btd_Covar(27,30) = Btd_110um_067um_Btd_110um_097um_Covar_Land
   Btd_Covar(27,38) = Btd_110um_067um_Btd_110um_104um_Covar_Land
   Btd_Covar(27,32) = Btd_110um_067um_Btd_110um_120um_Covar_Land
   Btd_Covar(27,33) = Btd_110um_067um_Btd_110um_133um_Covar_Land
   Btd_Covar(27,34) = Btd_110um_067um_Btd_110um_136um_Covar_Land
   Btd_Covar(27,35) = Btd_110um_067um_Btd_110um_139um_Covar_Land
   Btd_Covar(27,36) = Btd_110um_067um_Btd_110um_142um_Covar_Land

   Btd_Covar(28,20) = Btd_110um_073um_Btd_110um_038um_Covar_Land
   Btd_Covar(28,37) = Btd_110um_073um_Btd_110um_062um_Covar_Land
   Btd_Covar(28,27) = Btd_110um_073um_Btd_110um_067um_Covar_Land
   Btd_Covar(28,28) = Btd_110um_073um_Btd_110um_073um_Covar_Land
   Btd_Covar(28,29) = Btd_110um_073um_Btd_110um_085um_Covar_Land
   Btd_Covar(28,30) = Btd_110um_073um_Btd_110um_097um_Covar_Land
   Btd_Covar(28,38) = Btd_110um_073um_Btd_110um_104um_Covar_Land
   Btd_Covar(28,32) = Btd_110um_073um_Btd_110um_120um_Covar_Land
   Btd_Covar(28,33) = Btd_110um_073um_Btd_110um_133um_Covar_Land
   Btd_Covar(28,34) = Btd_110um_073um_Btd_110um_136um_Covar_Land
   Btd_Covar(28,35) = Btd_110um_073um_Btd_110um_139um_Covar_Land
   Btd_Covar(28,36) = Btd_110um_073um_Btd_110um_142um_Covar_Land

   Btd_Covar(29,20) = Btd_110um_085um_Btd_110um_038um_Covar_Land
   Btd_Covar(29,37) = Btd_110um_085um_Btd_110um_062um_Covar_Land
   Btd_Covar(29,27) = Btd_110um_085um_Btd_110um_067um_Covar_Land
   Btd_Covar(29,28) = Btd_110um_085um_Btd_110um_073um_Covar_Land
   Btd_Covar(29,29) = Btd_110um_085um_Btd_110um_085um_Covar_Land
   Btd_Covar(29,30) = Btd_110um_085um_Btd_110um_097um_Covar_Land
   Btd_Covar(29,38) = Btd_110um_085um_Btd_110um_104um_Covar_Land
   Btd_Covar(29,32) = Btd_110um_085um_Btd_110um_120um_Covar_Land
   Btd_Covar(29,33) = Btd_110um_085um_Btd_110um_133um_Covar_Land
   Btd_Covar(29,34) = Btd_110um_085um_Btd_110um_136um_Covar_Land
   Btd_Covar(29,35) = Btd_110um_085um_Btd_110um_139um_Covar_Land
   Btd_Covar(29,36) = Btd_110um_085um_Btd_110um_142um_Covar_Land

   Btd_Covar(30,20) = Btd_110um_097um_Btd_110um_038um_Covar_Land
   Btd_Covar(30,37) = Btd_110um_097um_Btd_110um_062um_Covar_Land
   Btd_Covar(30,27) = Btd_110um_097um_Btd_110um_067um_Covar_Land
   Btd_Covar(30,28) = Btd_110um_097um_Btd_110um_073um_Covar_Land
   Btd_Covar(30,29) = Btd_110um_097um_Btd_110um_085um_Covar_Land
   Btd_Covar(30,30) = Btd_110um_097um_Btd_110um_097um_Covar_Land
   Btd_Covar(30,38) = Btd_110um_097um_Btd_110um_104um_Covar_Land
   Btd_Covar(30,32) = Btd_110um_097um_Btd_110um_120um_Covar_Land
   Btd_Covar(30,33) = Btd_110um_097um_Btd_110um_133um_Covar_Land
   Btd_Covar(30,34) = Btd_110um_097um_Btd_110um_136um_Covar_Land
   Btd_Covar(30,35) = Btd_110um_097um_Btd_110um_139um_Covar_Land
   Btd_Covar(30,36) = Btd_110um_097um_Btd_110um_142um_Covar_Land

   Btd_Covar(38,20) = Btd_110um_104um_Btd_110um_038um_Covar_Land
   Btd_Covar(38,37) = Btd_110um_104um_Btd_110um_062um_Covar_Land
   Btd_Covar(38,27) = Btd_110um_104um_Btd_110um_067um_Covar_Land
   Btd_Covar(38,28) = Btd_110um_104um_Btd_110um_073um_Covar_Land
   Btd_Covar(38,29) = Btd_110um_104um_Btd_110um_085um_Covar_Land
   Btd_Covar(38,30) = Btd_110um_104um_Btd_110um_097um_Covar_Land
   Btd_Covar(38,38) = Btd_110um_104um_Btd_110um_104um_Covar_Land
   Btd_Covar(38,32) = Btd_110um_104um_Btd_110um_120um_Covar_Land
   Btd_Covar(38,33) = Btd_110um_104um_Btd_110um_133um_Covar_Land
   Btd_Covar(38,34) = Btd_110um_104um_Btd_110um_136um_Covar_Land
   Btd_Covar(38,35) = Btd_110um_104um_Btd_110um_139um_Covar_Land
   Btd_Covar(38,36) = Btd_110um_104um_Btd_110um_142um_Covar_Land

   Btd_Covar(32,20) = Btd_110um_120um_Btd_110um_038um_Covar_Land
   Btd_Covar(32,37) = Btd_110um_120um_Btd_110um_062um_Covar_Land
   Btd_Covar(32,27) = Btd_110um_120um_Btd_110um_067um_Covar_Land
   Btd_Covar(32,28) = Btd_110um_120um_Btd_110um_073um_Covar_Land
   Btd_Covar(32,29) = Btd_110um_120um_Btd_110um_085um_Covar_Land
   Btd_Covar(32,30) = Btd_110um_120um_Btd_110um_097um_Covar_Land
   Btd_Covar(32,38) = Btd_110um_120um_Btd_110um_104um_Covar_Land
   Btd_Covar(32,32) = Btd_110um_120um_Btd_110um_120um_Covar_Land
   Btd_Covar(32,33) = Btd_110um_120um_Btd_110um_133um_Covar_Land
   Btd_Covar(32,34) = Btd_110um_120um_Btd_110um_136um_Covar_Land
   Btd_Covar(32,35) = Btd_110um_120um_Btd_110um_139um_Covar_Land
   Btd_Covar(32,36) = Btd_110um_120um_Btd_110um_142um_Covar_Land

   Btd_Covar(33,20) = Btd_110um_133um_Btd_110um_038um_Covar_Land
   Btd_Covar(33,37) = Btd_110um_133um_Btd_110um_062um_Covar_Land
   Btd_Covar(33,27) = Btd_110um_133um_Btd_110um_067um_Covar_Land
   Btd_Covar(33,28) = Btd_110um_133um_Btd_110um_073um_Covar_Land
   Btd_Covar(33,29) = Btd_110um_133um_Btd_110um_085um_Covar_Land
   Btd_Covar(33,30) = Btd_110um_133um_Btd_110um_097um_Covar_Land
   Btd_Covar(33,38) = Btd_110um_133um_Btd_110um_104um_Covar_Land
   Btd_Covar(33,32) = Btd_110um_133um_Btd_110um_120um_Covar_Land
   Btd_Covar(33,33) = Btd_110um_133um_Btd_110um_133um_Covar_Land
   Btd_Covar(33,34) = Btd_110um_133um_Btd_110um_136um_Covar_Land
   Btd_Covar(33,35) = Btd_110um_133um_Btd_110um_139um_Covar_Land
   Btd_Covar(33,36) = Btd_110um_133um_Btd_110um_142um_Covar_Land

   Btd_Covar(34,20) = Btd_110um_136um_Btd_110um_038um_Covar_Land
   Btd_Covar(34,37) = Btd_110um_136um_Btd_110um_062um_Covar_Land
   Btd_Covar(34,27) = Btd_110um_136um_Btd_110um_067um_Covar_Land
   Btd_Covar(34,28) = Btd_110um_136um_Btd_110um_073um_Covar_Land
   Btd_Covar(34,29) = Btd_110um_136um_Btd_110um_085um_Covar_Land
   Btd_Covar(34,30) = Btd_110um_136um_Btd_110um_097um_Covar_Land
   Btd_Covar(34,38) = Btd_110um_136um_Btd_110um_104um_Covar_Land
   Btd_Covar(34,32) = Btd_110um_136um_Btd_110um_120um_Covar_Land
   Btd_Covar(34,33) = Btd_110um_136um_Btd_110um_133um_Covar_Land
   Btd_Covar(34,34) = Btd_110um_136um_Btd_110um_136um_Covar_Land
   Btd_Covar(34,35) = Btd_110um_136um_Btd_110um_139um_Covar_Land
   Btd_Covar(34,36) = Btd_110um_136um_Btd_110um_142um_Covar_Land

   Btd_Covar(35,20) = Btd_110um_139um_Btd_110um_038um_Covar_Land
   Btd_Covar(35,37) = Btd_110um_139um_Btd_110um_062um_Covar_Land
   Btd_Covar(35,27) = Btd_110um_139um_Btd_110um_067um_Covar_Land
   Btd_Covar(35,28) = Btd_110um_139um_Btd_110um_073um_Covar_Land
   Btd_Covar(35,29) = Btd_110um_139um_Btd_110um_085um_Covar_Land
   Btd_Covar(35,30) = Btd_110um_139um_Btd_110um_097um_Covar_Land
   Btd_Covar(35,38) = Btd_110um_139um_Btd_110um_104um_Covar_Land
   Btd_Covar(35,32) = Btd_110um_139um_Btd_110um_120um_Covar_Land
   Btd_Covar(35,33) = Btd_110um_139um_Btd_110um_133um_Covar_Land
   Btd_Covar(35,34) = Btd_110um_139um_Btd_110um_136um_Covar_Land
   Btd_Covar(35,35) = Btd_110um_139um_Btd_110um_139um_Covar_Land
   Btd_Covar(35,36) = Btd_110um_139um_Btd_110um_142um_Covar_Land

   Btd_Covar(36,20) = Btd_110um_142um_Btd_110um_038um_Covar_Land
   Btd_Covar(36,37) = Btd_110um_142um_Btd_110um_062um_Covar_Land
   Btd_Covar(36,27) = Btd_110um_142um_Btd_110um_067um_Covar_Land
   Btd_Covar(36,28) = Btd_110um_142um_Btd_110um_073um_Covar_Land
   Btd_Covar(36,29) = Btd_110um_142um_Btd_110um_085um_Covar_Land
   Btd_Covar(36,30) = Btd_110um_142um_Btd_110um_097um_Covar_Land
   Btd_Covar(36,38) = Btd_110um_142um_Btd_110um_104um_Covar_Land
   Btd_Covar(36,32) = Btd_110um_142um_Btd_110um_120um_Covar_Land
   Btd_Covar(36,33) = Btd_110um_142um_Btd_110um_133um_Covar_Land
   Btd_Covar(36,34) = Btd_110um_142um_Btd_110um_136um_Covar_Land
   Btd_Covar(36,35) = Btd_110um_142um_Btd_110um_139um_Covar_Land
   Btd_Covar(36,36) = Btd_110um_142um_Btd_110um_142um_Covar_Land

   !--- Snow
   case(2)
   Bt_Covar(31) = Bt_110um_Bt_110um_Covar_Snow
   Bt_Covar(20) = Bt_110um_Btd_110um_038um_Covar_Snow
   Bt_Covar(37) = Bt_110um_Btd_110um_062um_Covar_Snow
   Bt_Covar(27) = Bt_110um_Btd_110um_067um_Covar_Snow
   Bt_Covar(28) = Bt_110um_Btd_110um_073um_Covar_Snow
   Bt_Covar(29) = Bt_110um_Btd_110um_085um_Covar_Snow
   Bt_Covar(30) = Bt_110um_Btd_110um_097um_Covar_Snow
   Bt_Covar(38) = Bt_110um_Btd_110um_104um_Covar_Snow
   Bt_Covar(32) = Bt_110um_Btd_110um_120um_Covar_Snow
   Bt_Covar(33) = Bt_110um_Btd_110um_133um_Covar_Snow
   Bt_Covar(34) = Bt_110um_Btd_110um_136um_Covar_Snow
   Bt_Covar(35) = Bt_110um_Btd_110um_139um_Covar_Snow
   Bt_Covar(36) = Bt_110um_Btd_110um_142um_Covar_Snow

   Btd_Covar(20,20) = Btd_110um_038um_Btd_110um_038um_Covar_Snow
   Btd_Covar(20,37) = Btd_110um_038um_Btd_110um_062um_Covar_Snow
   Btd_Covar(20,27) = Btd_110um_038um_Btd_110um_067um_Covar_Snow
   Btd_Covar(20,28) = Btd_110um_038um_Btd_110um_073um_Covar_Snow
   Btd_Covar(20,29) = Btd_110um_038um_Btd_110um_085um_Covar_Snow
   Btd_Covar(20,30) = Btd_110um_038um_Btd_110um_097um_Covar_Snow
   Btd_Covar(20,38) = Btd_110um_038um_Btd_110um_104um_Covar_Snow
   Btd_Covar(20,32) = Btd_110um_038um_Btd_110um_120um_Covar_Snow
   Btd_Covar(20,33) = Btd_110um_038um_Btd_110um_133um_Covar_Snow
   Btd_Covar(20,34) = Btd_110um_038um_Btd_110um_136um_Covar_Snow
   Btd_Covar(20,35) = Btd_110um_038um_Btd_110um_139um_Covar_Snow
   Btd_Covar(20,36) = Btd_110um_038um_Btd_110um_142um_Covar_Snow

   Btd_Covar(37,20) = Btd_110um_062um_Btd_110um_038um_Covar_Snow
   Btd_Covar(37,37) = Btd_110um_062um_Btd_110um_062um_Covar_Snow
   Btd_Covar(37,27) = Btd_110um_062um_Btd_110um_067um_Covar_Snow
   Btd_Covar(37,28) = Btd_110um_062um_Btd_110um_073um_Covar_Snow
   Btd_Covar(37,29) = Btd_110um_062um_Btd_110um_085um_Covar_Snow
   Btd_Covar(37,30) = Btd_110um_062um_Btd_110um_097um_Covar_Snow
   Btd_Covar(37,38) = Btd_110um_062um_Btd_110um_104um_Covar_Snow
   Btd_Covar(37,32) = Btd_110um_062um_Btd_110um_120um_Covar_Snow
   Btd_Covar(37,33) = Btd_110um_062um_Btd_110um_133um_Covar_Snow
   Btd_Covar(37,34) = Btd_110um_062um_Btd_110um_136um_Covar_Snow
   Btd_Covar(37,35) = Btd_110um_062um_Btd_110um_139um_Covar_Snow
   Btd_Covar(37,36) = Btd_110um_062um_Btd_110um_142um_Covar_Snow

   Btd_Covar(27,20) = Btd_110um_067um_Btd_110um_038um_Covar_Snow
   Btd_Covar(27,37) = Btd_110um_067um_Btd_110um_062um_Covar_Snow
   Btd_Covar(27,27) = Btd_110um_067um_Btd_110um_067um_Covar_Snow
   Btd_Covar(27,28) = Btd_110um_067um_Btd_110um_073um_Covar_Snow
   Btd_Covar(27,29) = Btd_110um_067um_Btd_110um_085um_Covar_Snow
   Btd_Covar(27,30) = Btd_110um_067um_Btd_110um_097um_Covar_Snow
   Btd_Covar(27,38) = Btd_110um_067um_Btd_110um_104um_Covar_Snow
   Btd_Covar(27,32) = Btd_110um_067um_Btd_110um_120um_Covar_Snow
   Btd_Covar(27,33) = Btd_110um_067um_Btd_110um_133um_Covar_Snow
   Btd_Covar(27,34) = Btd_110um_067um_Btd_110um_136um_Covar_Snow
   Btd_Covar(27,35) = Btd_110um_067um_Btd_110um_139um_Covar_Snow
   Btd_Covar(27,36) = Btd_110um_067um_Btd_110um_142um_Covar_Snow

   Btd_Covar(28,20) = Btd_110um_073um_Btd_110um_038um_Covar_Snow
   Btd_Covar(28,37) = Btd_110um_073um_Btd_110um_062um_Covar_Snow
   Btd_Covar(28,27) = Btd_110um_073um_Btd_110um_067um_Covar_Snow
   Btd_Covar(28,28) = Btd_110um_073um_Btd_110um_073um_Covar_Snow
   Btd_Covar(28,29) = Btd_110um_073um_Btd_110um_085um_Covar_Snow
   Btd_Covar(28,30) = Btd_110um_073um_Btd_110um_097um_Covar_Snow
   Btd_Covar(28,38) = Btd_110um_073um_Btd_110um_104um_Covar_Snow
   Btd_Covar(28,32) = Btd_110um_073um_Btd_110um_120um_Covar_Snow
   Btd_Covar(28,33) = Btd_110um_073um_Btd_110um_133um_Covar_Snow
   Btd_Covar(28,34) = Btd_110um_073um_Btd_110um_136um_Covar_Snow
   Btd_Covar(28,35) = Btd_110um_073um_Btd_110um_139um_Covar_Snow
   Btd_Covar(28,36) = Btd_110um_073um_Btd_110um_142um_Covar_Snow

   Btd_Covar(29,20) = Btd_110um_085um_Btd_110um_038um_Covar_Snow
   Btd_Covar(29,37) = Btd_110um_085um_Btd_110um_062um_Covar_Snow
   Btd_Covar(29,27) = Btd_110um_085um_Btd_110um_067um_Covar_Snow
   Btd_Covar(29,28) = Btd_110um_085um_Btd_110um_073um_Covar_Snow
   Btd_Covar(29,29) = Btd_110um_085um_Btd_110um_085um_Covar_Snow
   Btd_Covar(29,30) = Btd_110um_085um_Btd_110um_097um_Covar_Snow
   Btd_Covar(29,38) = Btd_110um_085um_Btd_110um_104um_Covar_Snow
   Btd_Covar(29,32) = Btd_110um_085um_Btd_110um_120um_Covar_Snow
   Btd_Covar(29,33) = Btd_110um_085um_Btd_110um_133um_Covar_Snow
   Btd_Covar(29,34) = Btd_110um_085um_Btd_110um_136um_Covar_Snow
   Btd_Covar(29,35) = Btd_110um_085um_Btd_110um_139um_Covar_Snow
   Btd_Covar(29,36) = Btd_110um_085um_Btd_110um_142um_Covar_Snow

   Btd_Covar(30,20) = Btd_110um_097um_Btd_110um_038um_Covar_Snow
   Btd_Covar(30,37) = Btd_110um_097um_Btd_110um_062um_Covar_Snow
   Btd_Covar(30,27) = Btd_110um_097um_Btd_110um_067um_Covar_Snow
   Btd_Covar(30,28) = Btd_110um_097um_Btd_110um_073um_Covar_Snow
   Btd_Covar(30,29) = Btd_110um_097um_Btd_110um_085um_Covar_Snow
   Btd_Covar(30,30) = Btd_110um_097um_Btd_110um_097um_Covar_Snow
   Btd_Covar(30,38) = Btd_110um_097um_Btd_110um_104um_Covar_Snow
   Btd_Covar(30,32) = Btd_110um_097um_Btd_110um_120um_Covar_Snow
   Btd_Covar(30,33) = Btd_110um_097um_Btd_110um_133um_Covar_Snow
   Btd_Covar(30,34) = Btd_110um_097um_Btd_110um_136um_Covar_Snow
   Btd_Covar(30,35) = Btd_110um_097um_Btd_110um_139um_Covar_Snow
   Btd_Covar(30,36) = Btd_110um_097um_Btd_110um_142um_Covar_Snow

   Btd_Covar(38,20) = Btd_110um_104um_Btd_110um_038um_Covar_Snow
   Btd_Covar(38,37) = Btd_110um_104um_Btd_110um_062um_Covar_Snow
   Btd_Covar(38,27) = Btd_110um_104um_Btd_110um_067um_Covar_Snow
   Btd_Covar(38,28) = Btd_110um_104um_Btd_110um_073um_Covar_Snow
   Btd_Covar(38,29) = Btd_110um_104um_Btd_110um_085um_Covar_Snow
   Btd_Covar(38,30) = Btd_110um_104um_Btd_110um_097um_Covar_Snow
   Btd_Covar(38,38) = Btd_110um_104um_Btd_110um_104um_Covar_Snow
   Btd_Covar(38,32) = Btd_110um_104um_Btd_110um_120um_Covar_Snow
   Btd_Covar(38,33) = Btd_110um_104um_Btd_110um_133um_Covar_Snow
   Btd_Covar(38,34) = Btd_110um_104um_Btd_110um_136um_Covar_Snow
   Btd_Covar(38,35) = Btd_110um_104um_Btd_110um_139um_Covar_Snow
   Btd_Covar(38,36) = Btd_110um_104um_Btd_110um_142um_Covar_Snow

   Btd_Covar(32,20) = Btd_110um_120um_Btd_110um_038um_Covar_Snow
   Btd_Covar(32,37) = Btd_110um_120um_Btd_110um_062um_Covar_Snow
   Btd_Covar(32,27) = Btd_110um_120um_Btd_110um_067um_Covar_Snow
   Btd_Covar(32,28) = Btd_110um_120um_Btd_110um_073um_Covar_Snow
   Btd_Covar(32,29) = Btd_110um_120um_Btd_110um_085um_Covar_Snow
   Btd_Covar(32,30) = Btd_110um_120um_Btd_110um_097um_Covar_Snow
   Btd_Covar(32,38) = Btd_110um_120um_Btd_110um_104um_Covar_Snow
   Btd_Covar(32,32) = Btd_110um_120um_Btd_110um_120um_Covar_Snow
   Btd_Covar(32,33) = Btd_110um_120um_Btd_110um_133um_Covar_Snow
   Btd_Covar(32,34) = Btd_110um_120um_Btd_110um_136um_Covar_Snow
   Btd_Covar(32,35) = Btd_110um_120um_Btd_110um_139um_Covar_Snow
   Btd_Covar(32,36) = Btd_110um_120um_Btd_110um_142um_Covar_Snow

   Btd_Covar(33,20) = Btd_110um_133um_Btd_110um_038um_Covar_Snow
   Btd_Covar(33,37) = Btd_110um_133um_Btd_110um_062um_Covar_Snow
   Btd_Covar(33,27) = Btd_110um_133um_Btd_110um_067um_Covar_Snow
   Btd_Covar(33,28) = Btd_110um_133um_Btd_110um_073um_Covar_Snow
   Btd_Covar(33,29) = Btd_110um_133um_Btd_110um_085um_Covar_Snow
   Btd_Covar(33,30) = Btd_110um_133um_Btd_110um_097um_Covar_Snow
   Btd_Covar(33,38) = Btd_110um_133um_Btd_110um_104um_Covar_Snow
   Btd_Covar(33,32) = Btd_110um_133um_Btd_110um_120um_Covar_Snow
   Btd_Covar(33,33) = Btd_110um_133um_Btd_110um_133um_Covar_Snow
   Btd_Covar(33,34) = Btd_110um_133um_Btd_110um_136um_Covar_Snow
   Btd_Covar(33,35) = Btd_110um_133um_Btd_110um_139um_Covar_Snow
   Btd_Covar(33,36) = Btd_110um_133um_Btd_110um_142um_Covar_Snow

   Btd_Covar(34,20) = Btd_110um_136um_Btd_110um_038um_Covar_Snow
   Btd_Covar(34,37) = Btd_110um_136um_Btd_110um_062um_Covar_Snow
   Btd_Covar(34,27) = Btd_110um_136um_Btd_110um_067um_Covar_Snow
   Btd_Covar(34,28) = Btd_110um_136um_Btd_110um_073um_Covar_Snow
   Btd_Covar(34,29) = Btd_110um_136um_Btd_110um_085um_Covar_Snow
   Btd_Covar(34,30) = Btd_110um_136um_Btd_110um_097um_Covar_Snow
   Btd_Covar(34,38) = Btd_110um_136um_Btd_110um_104um_Covar_Snow
   Btd_Covar(34,32) = Btd_110um_136um_Btd_110um_120um_Covar_Snow
   Btd_Covar(34,33) = Btd_110um_136um_Btd_110um_133um_Covar_Snow
   Btd_Covar(34,34) = Btd_110um_136um_Btd_110um_136um_Covar_Snow
   Btd_Covar(34,35) = Btd_110um_136um_Btd_110um_139um_Covar_Snow
   Btd_Covar(34,36) = Btd_110um_136um_Btd_110um_142um_Covar_Snow

   Btd_Covar(35,20) = Btd_110um_139um_Btd_110um_038um_Covar_Snow
   Btd_Covar(35,37) = Btd_110um_139um_Btd_110um_062um_Covar_Snow
   Btd_Covar(35,27) = Btd_110um_139um_Btd_110um_067um_Covar_Snow
   Btd_Covar(35,28) = Btd_110um_139um_Btd_110um_073um_Covar_Snow
   Btd_Covar(35,29) = Btd_110um_139um_Btd_110um_085um_Covar_Snow
   Btd_Covar(35,30) = Btd_110um_139um_Btd_110um_097um_Covar_Snow
   Btd_Covar(35,38) = Btd_110um_139um_Btd_110um_104um_Covar_Snow
   Btd_Covar(35,32) = Btd_110um_139um_Btd_110um_120um_Covar_Snow
   Btd_Covar(35,33) = Btd_110um_139um_Btd_110um_133um_Covar_Snow
   Btd_Covar(35,34) = Btd_110um_139um_Btd_110um_136um_Covar_Snow
   Btd_Covar(35,35) = Btd_110um_139um_Btd_110um_139um_Covar_Snow
   Btd_Covar(35,36) = Btd_110um_139um_Btd_110um_142um_Covar_Snow

   Btd_Covar(36,20) = Btd_110um_142um_Btd_110um_038um_Covar_Snow
   Btd_Covar(36,37) = Btd_110um_142um_Btd_110um_062um_Covar_Snow
   Btd_Covar(36,27) = Btd_110um_142um_Btd_110um_067um_Covar_Snow
   Btd_Covar(36,28) = Btd_110um_142um_Btd_110um_073um_Covar_Snow
   Btd_Covar(36,29) = Btd_110um_142um_Btd_110um_085um_Covar_Snow
   Btd_Covar(36,30) = Btd_110um_142um_Btd_110um_097um_Covar_Snow
   Btd_Covar(36,38) = Btd_110um_142um_Btd_110um_104um_Covar_Snow
   Btd_Covar(36,32) = Btd_110um_142um_Btd_110um_120um_Covar_Snow
   Btd_Covar(36,33) = Btd_110um_142um_Btd_110um_133um_Covar_Snow
   Btd_Covar(36,34) = Btd_110um_142um_Btd_110um_136um_Covar_Snow
   Btd_Covar(36,35) = Btd_110um_142um_Btd_110um_139um_Covar_Snow
   Btd_Covar(36,36) = Btd_110um_142um_Btd_110um_142um_Covar_Snow

   !--- Desert
   case (3)
   Bt_Covar(31) = Bt_110um_Bt_110um_Covar_Desert
   Bt_Covar(20) = Bt_110um_Btd_110um_038um_Covar_Desert
   Bt_Covar(37) = Bt_110um_Btd_110um_062um_Covar_Desert
   Bt_Covar(27) = Bt_110um_Btd_110um_067um_Covar_Desert
   Bt_Covar(28) = Bt_110um_Btd_110um_073um_Covar_Desert
   Bt_Covar(29) = Bt_110um_Btd_110um_085um_Covar_Desert
   Bt_Covar(30) = Bt_110um_Btd_110um_097um_Covar_Desert
   Bt_Covar(38) = Bt_110um_Btd_110um_104um_Covar_Desert
   Bt_Covar(32) = Bt_110um_Btd_110um_120um_Covar_Desert
   Bt_Covar(33) = Bt_110um_Btd_110um_133um_Covar_Desert
   Bt_Covar(34) = Bt_110um_Btd_110um_136um_Covar_Desert
   Bt_Covar(35) = Bt_110um_Btd_110um_139um_Covar_Desert
   Bt_Covar(36) = Bt_110um_Btd_110um_142um_Covar_Desert

   Btd_Covar(20,20) = Btd_110um_038um_Btd_110um_038um_Covar_Desert
   Btd_Covar(20,37) = Btd_110um_038um_Btd_110um_062um_Covar_Desert
   Btd_Covar(20,27) = Btd_110um_038um_Btd_110um_067um_Covar_Desert
   Btd_Covar(20,28) = Btd_110um_038um_Btd_110um_073um_Covar_Desert
   Btd_Covar(20,29) = Btd_110um_038um_Btd_110um_085um_Covar_Desert
   Btd_Covar(20,30) = Btd_110um_038um_Btd_110um_097um_Covar_Desert
   Btd_Covar(20,38) = Btd_110um_038um_Btd_110um_104um_Covar_Desert
   Btd_Covar(20,32) = Btd_110um_038um_Btd_110um_120um_Covar_Desert
   Btd_Covar(20,33) = Btd_110um_038um_Btd_110um_133um_Covar_Desert
   Btd_Covar(20,34) = Btd_110um_038um_Btd_110um_136um_Covar_Desert
   Btd_Covar(20,35) = Btd_110um_038um_Btd_110um_139um_Covar_Desert
   Btd_Covar(20,36) = Btd_110um_038um_Btd_110um_142um_Covar_Desert

   Btd_Covar(37,20) = Btd_110um_062um_Btd_110um_038um_Covar_Desert
   Btd_Covar(37,37) = Btd_110um_062um_Btd_110um_062um_Covar_Desert
   Btd_Covar(37,27) = Btd_110um_062um_Btd_110um_067um_Covar_Desert
   Btd_Covar(37,28) = Btd_110um_062um_Btd_110um_073um_Covar_Desert
   Btd_Covar(37,29) = Btd_110um_062um_Btd_110um_085um_Covar_Desert
   Btd_Covar(37,30) = Btd_110um_062um_Btd_110um_097um_Covar_Desert
   Btd_Covar(37,38) = Btd_110um_062um_Btd_110um_104um_Covar_Desert
   Btd_Covar(37,32) = Btd_110um_062um_Btd_110um_120um_Covar_Desert
   Btd_Covar(37,33) = Btd_110um_062um_Btd_110um_133um_Covar_Desert
   Btd_Covar(37,34) = Btd_110um_062um_Btd_110um_136um_Covar_Desert
   Btd_Covar(37,35) = Btd_110um_062um_Btd_110um_139um_Covar_Desert
   Btd_Covar(37,36) = Btd_110um_062um_Btd_110um_142um_Covar_Desert

   Btd_Covar(27,20) = Btd_110um_067um_Btd_110um_038um_Covar_Desert
   Btd_Covar(27,37) = Btd_110um_067um_Btd_110um_062um_Covar_Desert
   Btd_Covar(27,27) = Btd_110um_067um_Btd_110um_067um_Covar_Desert
   Btd_Covar(27,28) = Btd_110um_067um_Btd_110um_073um_Covar_Desert
   Btd_Covar(27,29) = Btd_110um_067um_Btd_110um_085um_Covar_Desert
   Btd_Covar(27,30) = Btd_110um_067um_Btd_110um_097um_Covar_Desert
   Btd_Covar(27,38) = Btd_110um_067um_Btd_110um_104um_Covar_Desert
   Btd_Covar(27,32) = Btd_110um_067um_Btd_110um_120um_Covar_Desert
   Btd_Covar(27,33) = Btd_110um_067um_Btd_110um_133um_Covar_Desert
   Btd_Covar(27,34) = Btd_110um_067um_Btd_110um_136um_Covar_Desert
   Btd_Covar(27,35) = Btd_110um_067um_Btd_110um_139um_Covar_Desert
   Btd_Covar(27,36) = Btd_110um_067um_Btd_110um_142um_Covar_Desert

   Btd_Covar(28,20) = Btd_110um_073um_Btd_110um_038um_Covar_Desert
   Btd_Covar(28,37) = Btd_110um_073um_Btd_110um_062um_Covar_Desert
   Btd_Covar(28,27) = Btd_110um_073um_Btd_110um_067um_Covar_Desert
   Btd_Covar(28,28) = Btd_110um_073um_Btd_110um_073um_Covar_Desert
   Btd_Covar(28,29) = Btd_110um_073um_Btd_110um_085um_Covar_Desert
   Btd_Covar(28,30) = Btd_110um_073um_Btd_110um_097um_Covar_Desert
   Btd_Covar(28,38) = Btd_110um_073um_Btd_110um_104um_Covar_Desert
   Btd_Covar(28,32) = Btd_110um_073um_Btd_110um_120um_Covar_Desert
   Btd_Covar(28,33) = Btd_110um_073um_Btd_110um_133um_Covar_Desert
   Btd_Covar(28,34) = Btd_110um_073um_Btd_110um_136um_Covar_Desert
   Btd_Covar(28,35) = Btd_110um_073um_Btd_110um_139um_Covar_Desert
   Btd_Covar(28,36) = Btd_110um_073um_Btd_110um_142um_Covar_Desert

   Btd_Covar(29,20) = Btd_110um_085um_Btd_110um_038um_Covar_Desert
   Btd_Covar(29,37) = Btd_110um_085um_Btd_110um_062um_Covar_Desert
   Btd_Covar(29,27) = Btd_110um_085um_Btd_110um_067um_Covar_Desert
   Btd_Covar(29,28) = Btd_110um_085um_Btd_110um_073um_Covar_Desert
   Btd_Covar(29,29) = Btd_110um_085um_Btd_110um_085um_Covar_Desert
   Btd_Covar(29,30) = Btd_110um_085um_Btd_110um_097um_Covar_Desert
   Btd_Covar(29,38) = Btd_110um_085um_Btd_110um_104um_Covar_Desert
   Btd_Covar(29,32) = Btd_110um_085um_Btd_110um_120um_Covar_Desert
   Btd_Covar(29,33) = Btd_110um_085um_Btd_110um_133um_Covar_Desert
   Btd_Covar(29,34) = Btd_110um_085um_Btd_110um_136um_Covar_Desert
   Btd_Covar(29,35) = Btd_110um_085um_Btd_110um_139um_Covar_Desert
   Btd_Covar(29,36) = Btd_110um_085um_Btd_110um_142um_Covar_Desert

   Btd_Covar(30,20) = Btd_110um_097um_Btd_110um_038um_Covar_Desert
   Btd_Covar(30,37) = Btd_110um_097um_Btd_110um_062um_Covar_Desert
   Btd_Covar(30,27) = Btd_110um_097um_Btd_110um_067um_Covar_Desert
   Btd_Covar(30,28) = Btd_110um_097um_Btd_110um_073um_Covar_Desert
   Btd_Covar(30,29) = Btd_110um_097um_Btd_110um_085um_Covar_Desert
   Btd_Covar(30,30) = Btd_110um_097um_Btd_110um_097um_Covar_Desert
   Btd_Covar(30,38) = Btd_110um_097um_Btd_110um_104um_Covar_Desert
   Btd_Covar(30,32) = Btd_110um_097um_Btd_110um_120um_Covar_Desert
   Btd_Covar(30,33) = Btd_110um_097um_Btd_110um_133um_Covar_Desert
   Btd_Covar(30,34) = Btd_110um_097um_Btd_110um_136um_Covar_Desert
   Btd_Covar(30,35) = Btd_110um_097um_Btd_110um_139um_Covar_Desert
   Btd_Covar(30,36) = Btd_110um_097um_Btd_110um_142um_Covar_Desert

   Btd_Covar(38,20) = Btd_110um_104um_Btd_110um_038um_Covar_Desert
   Btd_Covar(38,37) = Btd_110um_104um_Btd_110um_062um_Covar_Desert
   Btd_Covar(38,27) = Btd_110um_104um_Btd_110um_067um_Covar_Desert
   Btd_Covar(38,28) = Btd_110um_104um_Btd_110um_073um_Covar_Desert
   Btd_Covar(38,29) = Btd_110um_104um_Btd_110um_085um_Covar_Desert
   Btd_Covar(38,30) = Btd_110um_104um_Btd_110um_097um_Covar_Desert
   Btd_Covar(38,38) = Btd_110um_104um_Btd_110um_104um_Covar_Desert
   Btd_Covar(38,32) = Btd_110um_104um_Btd_110um_120um_Covar_Desert
   Btd_Covar(38,33) = Btd_110um_104um_Btd_110um_133um_Covar_Desert
   Btd_Covar(38,34) = Btd_110um_104um_Btd_110um_136um_Covar_Desert
   Btd_Covar(38,35) = Btd_110um_104um_Btd_110um_139um_Covar_Desert
   Btd_Covar(38,36) = Btd_110um_104um_Btd_110um_142um_Covar_Desert

   Btd_Covar(32,20) = Btd_110um_120um_Btd_110um_038um_Covar_Desert
   Btd_Covar(32,37) = Btd_110um_120um_Btd_110um_062um_Covar_Desert
   Btd_Covar(32,27) = Btd_110um_120um_Btd_110um_067um_Covar_Desert
   Btd_Covar(32,28) = Btd_110um_120um_Btd_110um_073um_Covar_Desert
   Btd_Covar(32,29) = Btd_110um_120um_Btd_110um_085um_Covar_Desert
   Btd_Covar(32,30) = Btd_110um_120um_Btd_110um_097um_Covar_Desert
   Btd_Covar(32,38) = Btd_110um_120um_Btd_110um_104um_Covar_Desert
   Btd_Covar(32,32) = Btd_110um_120um_Btd_110um_120um_Covar_Desert
   Btd_Covar(32,33) = Btd_110um_120um_Btd_110um_133um_Covar_Desert
   Btd_Covar(32,34) = Btd_110um_120um_Btd_110um_136um_Covar_Desert
   Btd_Covar(32,35) = Btd_110um_120um_Btd_110um_139um_Covar_Desert
   Btd_Covar(32,36) = Btd_110um_120um_Btd_110um_142um_Covar_Desert

   Btd_Covar(33,20) = Btd_110um_133um_Btd_110um_038um_Covar_Desert
   Btd_Covar(33,37) = Btd_110um_133um_Btd_110um_062um_Covar_Desert
   Btd_Covar(33,27) = Btd_110um_133um_Btd_110um_067um_Covar_Desert
   Btd_Covar(33,28) = Btd_110um_133um_Btd_110um_073um_Covar_Desert
   Btd_Covar(33,29) = Btd_110um_133um_Btd_110um_085um_Covar_Desert
   Btd_Covar(33,30) = Btd_110um_133um_Btd_110um_097um_Covar_Desert
   Btd_Covar(33,38) = Btd_110um_133um_Btd_110um_104um_Covar_Desert
   Btd_Covar(33,32) = Btd_110um_133um_Btd_110um_120um_Covar_Desert
   Btd_Covar(33,33) = Btd_110um_133um_Btd_110um_133um_Covar_Desert
   Btd_Covar(33,34) = Btd_110um_133um_Btd_110um_136um_Covar_Desert
   Btd_Covar(33,35) = Btd_110um_133um_Btd_110um_139um_Covar_Desert
   Btd_Covar(33,36) = Btd_110um_133um_Btd_110um_142um_Covar_Desert

   Btd_Covar(34,20) = Btd_110um_136um_Btd_110um_038um_Covar_Desert
   Btd_Covar(34,37) = Btd_110um_136um_Btd_110um_062um_Covar_Desert
   Btd_Covar(34,27) = Btd_110um_136um_Btd_110um_067um_Covar_Desert
   Btd_Covar(34,28) = Btd_110um_136um_Btd_110um_073um_Covar_Desert
   Btd_Covar(34,29) = Btd_110um_136um_Btd_110um_085um_Covar_Desert
   Btd_Covar(34,30) = Btd_110um_136um_Btd_110um_097um_Covar_Desert
   Btd_Covar(34,38) = Btd_110um_136um_Btd_110um_104um_Covar_Desert
   Btd_Covar(34,32) = Btd_110um_136um_Btd_110um_120um_Covar_Desert
   Btd_Covar(34,33) = Btd_110um_136um_Btd_110um_133um_Covar_Desert
   Btd_Covar(34,34) = Btd_110um_136um_Btd_110um_136um_Covar_Desert
   Btd_Covar(34,35) = Btd_110um_136um_Btd_110um_139um_Covar_Desert
   Btd_Covar(34,36) = Btd_110um_136um_Btd_110um_142um_Covar_Desert

   Btd_Covar(35,20) = Btd_110um_139um_Btd_110um_038um_Covar_Desert
   Btd_Covar(35,37) = Btd_110um_139um_Btd_110um_062um_Covar_Desert
   Btd_Covar(35,27) = Btd_110um_139um_Btd_110um_067um_Covar_Desert
   Btd_Covar(35,28) = Btd_110um_139um_Btd_110um_073um_Covar_Desert
   Btd_Covar(35,29) = Btd_110um_139um_Btd_110um_085um_Covar_Desert
   Btd_Covar(35,30) = Btd_110um_139um_Btd_110um_097um_Covar_Desert
   Btd_Covar(35,38) = Btd_110um_139um_Btd_110um_104um_Covar_Desert
   Btd_Covar(35,32) = Btd_110um_139um_Btd_110um_120um_Covar_Desert
   Btd_Covar(35,33) = Btd_110um_139um_Btd_110um_133um_Covar_Desert
   Btd_Covar(35,34) = Btd_110um_139um_Btd_110um_136um_Covar_Desert
   Btd_Covar(35,35) = Btd_110um_139um_Btd_110um_139um_Covar_Desert
   Btd_Covar(35,36) = Btd_110um_139um_Btd_110um_142um_Covar_Desert

   Btd_Covar(36,20) = Btd_110um_142um_Btd_110um_038um_Covar_Desert
   Btd_Covar(36,37) = Btd_110um_142um_Btd_110um_062um_Covar_Desert
   Btd_Covar(36,27) = Btd_110um_142um_Btd_110um_067um_Covar_Desert
   Btd_Covar(36,28) = Btd_110um_142um_Btd_110um_073um_Covar_Desert
   Btd_Covar(36,29) = Btd_110um_142um_Btd_110um_085um_Covar_Desert
   Btd_Covar(36,30) = Btd_110um_142um_Btd_110um_097um_Covar_Desert
   Btd_Covar(36,38) = Btd_110um_142um_Btd_110um_104um_Covar_Desert
   Btd_Covar(36,32) = Btd_110um_142um_Btd_110um_120um_Covar_Desert
   Btd_Covar(36,33) = Btd_110um_142um_Btd_110um_133um_Covar_Desert
   Btd_Covar(36,34) = Btd_110um_142um_Btd_110um_136um_Covar_Desert
   Btd_Covar(36,35) = Btd_110um_142um_Btd_110um_139um_Covar_Desert
   Btd_Covar(36,36) = Btd_110um_142um_Btd_110um_142um_Covar_Desert

   !--- Arctic
   case (4)
   Bt_Covar(31) = Bt_110um_Bt_110um_Covar_Arctic
   Bt_Covar(20) = Bt_110um_Btd_110um_038um_Covar_Arctic
   Bt_Covar(37) = Bt_110um_Btd_110um_062um_Covar_Arctic
   Bt_Covar(27) = Bt_110um_Btd_110um_067um_Covar_Arctic
   Bt_Covar(28) = Bt_110um_Btd_110um_073um_Covar_Arctic
   Bt_Covar(29) = Bt_110um_Btd_110um_085um_Covar_Arctic
   Bt_Covar(30) = Bt_110um_Btd_110um_097um_Covar_Arctic
   Bt_Covar(38) = Bt_110um_Btd_110um_104um_Covar_Arctic
   Bt_Covar(32) = Bt_110um_Btd_110um_120um_Covar_Arctic
   Bt_Covar(33) = Bt_110um_Btd_110um_133um_Covar_Arctic
   Bt_Covar(34) = Bt_110um_Btd_110um_136um_Covar_Arctic
   Bt_Covar(35) = Bt_110um_Btd_110um_139um_Covar_Arctic
   Bt_Covar(36) = Bt_110um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(20,20) = Btd_110um_038um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(20,37) = Btd_110um_038um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(20,27) = Btd_110um_038um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(20,28) = Btd_110um_038um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(20,29) = Btd_110um_038um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(20,30) = Btd_110um_038um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(20,38) = Btd_110um_038um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(20,32) = Btd_110um_038um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(20,33) = Btd_110um_038um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(20,34) = Btd_110um_038um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(20,35) = Btd_110um_038um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(20,36) = Btd_110um_038um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(37,20) = Btd_110um_062um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(37,37) = Btd_110um_062um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(37,27) = Btd_110um_062um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(37,28) = Btd_110um_062um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(37,29) = Btd_110um_062um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(37,30) = Btd_110um_062um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(37,38) = Btd_110um_062um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(37,32) = Btd_110um_062um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(37,33) = Btd_110um_062um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(37,34) = Btd_110um_062um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(37,35) = Btd_110um_062um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(37,36) = Btd_110um_062um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(27,20) = Btd_110um_067um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(27,37) = Btd_110um_067um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(27,27) = Btd_110um_067um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(27,28) = Btd_110um_067um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(27,29) = Btd_110um_067um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(27,30) = Btd_110um_067um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(27,38) = Btd_110um_067um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(27,32) = Btd_110um_067um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(27,33) = Btd_110um_067um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(27,34) = Btd_110um_067um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(27,35) = Btd_110um_067um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(27,36) = Btd_110um_067um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(28,20) = Btd_110um_073um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(28,37) = Btd_110um_073um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(28,27) = Btd_110um_073um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(28,28) = Btd_110um_073um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(28,29) = Btd_110um_073um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(28,30) = Btd_110um_073um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(28,38) = Btd_110um_073um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(28,32) = Btd_110um_073um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(28,33) = Btd_110um_073um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(28,34) = Btd_110um_073um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(28,35) = Btd_110um_073um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(28,36) = Btd_110um_073um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(29,20) = Btd_110um_085um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(29,37) = Btd_110um_085um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(29,27) = Btd_110um_085um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(29,28) = Btd_110um_085um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(29,29) = Btd_110um_085um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(29,30) = Btd_110um_085um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(29,38) = Btd_110um_085um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(29,32) = Btd_110um_085um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(29,33) = Btd_110um_085um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(29,34) = Btd_110um_085um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(29,35) = Btd_110um_085um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(29,36) = Btd_110um_085um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(30,20) = Btd_110um_097um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(30,37) = Btd_110um_097um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(30,27) = Btd_110um_097um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(30,28) = Btd_110um_097um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(30,29) = Btd_110um_097um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(30,30) = Btd_110um_097um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(30,38) = Btd_110um_097um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(30,32) = Btd_110um_097um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(30,33) = Btd_110um_097um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(30,34) = Btd_110um_097um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(30,35) = Btd_110um_097um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(30,36) = Btd_110um_097um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(38,20) = Btd_110um_104um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(38,37) = Btd_110um_104um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(38,27) = Btd_110um_104um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(38,28) = Btd_110um_104um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(38,29) = Btd_110um_104um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(38,30) = Btd_110um_104um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(38,38) = Btd_110um_104um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(38,32) = Btd_110um_104um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(38,33) = Btd_110um_104um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(38,34) = Btd_110um_104um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(38,35) = Btd_110um_104um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(38,36) = Btd_110um_104um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(32,20) = Btd_110um_120um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(32,37) = Btd_110um_120um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(32,27) = Btd_110um_120um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(32,28) = Btd_110um_120um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(32,29) = Btd_110um_120um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(32,30) = Btd_110um_120um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(32,38) = Btd_110um_120um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(32,32) = Btd_110um_120um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(32,33) = Btd_110um_120um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(32,34) = Btd_110um_120um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(32,35) = Btd_110um_120um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(32,36) = Btd_110um_120um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(33,20) = Btd_110um_133um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(33,37) = Btd_110um_133um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(33,27) = Btd_110um_133um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(33,28) = Btd_110um_133um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(33,29) = Btd_110um_133um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(33,30) = Btd_110um_133um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(33,38) = Btd_110um_133um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(33,32) = Btd_110um_133um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(33,33) = Btd_110um_133um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(33,34) = Btd_110um_133um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(33,35) = Btd_110um_133um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(33,36) = Btd_110um_133um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(34,20) = Btd_110um_136um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(34,37) = Btd_110um_136um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(34,27) = Btd_110um_136um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(34,28) = Btd_110um_136um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(34,29) = Btd_110um_136um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(34,30) = Btd_110um_136um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(34,38) = Btd_110um_136um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(34,32) = Btd_110um_136um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(34,33) = Btd_110um_136um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(34,34) = Btd_110um_136um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(34,35) = Btd_110um_136um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(34,36) = Btd_110um_136um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(35,20) = Btd_110um_139um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(35,37) = Btd_110um_139um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(35,27) = Btd_110um_139um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(35,28) = Btd_110um_139um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(35,29) = Btd_110um_139um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(35,30) = Btd_110um_139um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(35,38) = Btd_110um_139um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(35,32) = Btd_110um_139um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(35,33) = Btd_110um_139um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(35,34) = Btd_110um_139um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(35,35) = Btd_110um_139um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(35,36) = Btd_110um_139um_Btd_110um_142um_Covar_Arctic

   Btd_Covar(36,20) = Btd_110um_142um_Btd_110um_038um_Covar_Arctic
   Btd_Covar(36,37) = Btd_110um_142um_Btd_110um_062um_Covar_Arctic
   Btd_Covar(36,27) = Btd_110um_142um_Btd_110um_067um_Covar_Arctic
   Btd_Covar(36,28) = Btd_110um_142um_Btd_110um_073um_Covar_Arctic
   Btd_Covar(36,29) = Btd_110um_142um_Btd_110um_085um_Covar_Arctic
   Btd_Covar(36,30) = Btd_110um_142um_Btd_110um_097um_Covar_Arctic
   Btd_Covar(36,38) = Btd_110um_142um_Btd_110um_104um_Covar_Arctic
   Btd_Covar(36,32) = Btd_110um_142um_Btd_110um_120um_Covar_Arctic
   Btd_Covar(36,33) = Btd_110um_142um_Btd_110um_133um_Covar_Arctic
   Btd_Covar(36,34) = Btd_110um_142um_Btd_110um_136um_Covar_Arctic
   Btd_Covar(36,35) = Btd_110um_142um_Btd_110um_139um_Covar_Arctic
   Btd_Covar(36,36) = Btd_110um_142um_Btd_110um_142um_Covar_Arctic

   !--- Antarctic
   case (5) 

   Bt_Covar(31) = Bt_110um_Bt_110um_Covar_Antarctic
   Bt_Covar(20) = Bt_110um_Btd_110um_038um_Covar_Antarctic
   Bt_Covar(37) = Bt_110um_Btd_110um_062um_Covar_Antarctic
   Bt_Covar(27) = Bt_110um_Btd_110um_067um_Covar_Antarctic
   Bt_Covar(28) = Bt_110um_Btd_110um_073um_Covar_Antarctic
   Bt_Covar(29) = Bt_110um_Btd_110um_085um_Covar_Antarctic
   Bt_Covar(30) = Bt_110um_Btd_110um_097um_Covar_Antarctic
   Bt_Covar(38) = Bt_110um_Btd_110um_104um_Covar_Antarctic
   Bt_Covar(32) = Bt_110um_Btd_110um_120um_Covar_Antarctic
   Bt_Covar(33) = Bt_110um_Btd_110um_133um_Covar_Antarctic
   Bt_Covar(34) = Bt_110um_Btd_110um_136um_Covar_Antarctic
   Bt_Covar(35) = Bt_110um_Btd_110um_139um_Covar_Antarctic
   Bt_Covar(36) = Bt_110um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(20,20) = Btd_110um_038um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(20,37) = Btd_110um_038um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(20,27) = Btd_110um_038um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(20,28) = Btd_110um_038um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(20,29) = Btd_110um_038um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(20,30) = Btd_110um_038um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(20,38) = Btd_110um_038um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(20,32) = Btd_110um_038um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(20,33) = Btd_110um_038um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(20,34) = Btd_110um_038um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(20,35) = Btd_110um_038um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(20,36) = Btd_110um_038um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(37,20) = Btd_110um_062um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(37,37) = Btd_110um_062um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(37,27) = Btd_110um_062um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(37,28) = Btd_110um_062um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(37,29) = Btd_110um_062um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(37,30) = Btd_110um_062um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(37,38) = Btd_110um_062um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(37,32) = Btd_110um_062um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(37,33) = Btd_110um_062um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(37,34) = Btd_110um_062um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(37,35) = Btd_110um_062um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(37,36) = Btd_110um_062um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(27,20) = Btd_110um_067um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(27,37) = Btd_110um_067um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(27,27) = Btd_110um_067um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(27,28) = Btd_110um_067um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(27,29) = Btd_110um_067um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(27,30) = Btd_110um_067um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(27,38) = Btd_110um_067um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(27,32) = Btd_110um_067um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(27,33) = Btd_110um_067um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(27,34) = Btd_110um_067um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(27,35) = Btd_110um_067um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(27,36) = Btd_110um_067um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(28,20) = Btd_110um_073um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(28,37) = Btd_110um_073um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(28,27) = Btd_110um_073um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(28,28) = Btd_110um_073um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(28,29) = Btd_110um_073um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(28,30) = Btd_110um_073um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(28,38) = Btd_110um_073um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(28,32) = Btd_110um_073um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(28,33) = Btd_110um_073um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(28,34) = Btd_110um_073um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(28,35) = Btd_110um_073um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(28,36) = Btd_110um_073um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(29,20) = Btd_110um_085um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(29,37) = Btd_110um_085um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(29,27) = Btd_110um_085um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(29,28) = Btd_110um_085um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(29,29) = Btd_110um_085um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(29,30) = Btd_110um_085um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(29,38) = Btd_110um_085um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(29,32) = Btd_110um_085um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(29,33) = Btd_110um_085um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(29,34) = Btd_110um_085um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(29,35) = Btd_110um_085um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(29,36) = Btd_110um_085um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(30,20) = Btd_110um_097um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(30,37) = Btd_110um_097um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(30,27) = Btd_110um_097um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(30,28) = Btd_110um_097um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(30,29) = Btd_110um_097um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(30,30) = Btd_110um_097um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(30,38) = Btd_110um_097um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(30,32) = Btd_110um_097um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(30,33) = Btd_110um_097um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(30,34) = Btd_110um_097um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(30,35) = Btd_110um_097um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(30,36) = Btd_110um_097um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(38,20) = Btd_110um_104um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(38,37) = Btd_110um_104um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(38,27) = Btd_110um_104um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(38,28) = Btd_110um_104um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(38,29) = Btd_110um_104um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(38,30) = Btd_110um_104um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(38,38) = Btd_110um_104um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(38,32) = Btd_110um_104um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(38,33) = Btd_110um_104um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(38,34) = Btd_110um_104um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(38,35) = Btd_110um_104um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(38,36) = Btd_110um_104um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(32,20) = Btd_110um_120um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(32,37) = Btd_110um_120um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(32,27) = Btd_110um_120um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(32,28) = Btd_110um_120um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(32,29) = Btd_110um_120um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(32,30) = Btd_110um_120um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(32,38) = Btd_110um_120um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(32,32) = Btd_110um_120um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(32,33) = Btd_110um_120um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(32,34) = Btd_110um_120um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(32,35) = Btd_110um_120um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(32,36) = Btd_110um_120um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(33,20) = Btd_110um_133um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(33,37) = Btd_110um_133um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(33,27) = Btd_110um_133um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(33,28) = Btd_110um_133um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(33,29) = Btd_110um_133um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(33,30) = Btd_110um_133um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(33,38) = Btd_110um_133um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(33,32) = Btd_110um_133um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(33,33) = Btd_110um_133um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(33,34) = Btd_110um_133um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(33,35) = Btd_110um_133um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(33,36) = Btd_110um_133um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(34,20) = Btd_110um_136um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(34,37) = Btd_110um_136um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(34,27) = Btd_110um_136um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(34,28) = Btd_110um_136um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(34,29) = Btd_110um_136um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(34,30) = Btd_110um_136um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(34,38) = Btd_110um_136um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(34,32) = Btd_110um_136um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(34,33) = Btd_110um_136um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(34,34) = Btd_110um_136um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(34,35) = Btd_110um_136um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(34,36) = Btd_110um_136um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(35,20) = Btd_110um_139um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(35,37) = Btd_110um_139um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(35,27) = Btd_110um_139um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(35,28) = Btd_110um_139um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(35,29) = Btd_110um_139um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(35,30) = Btd_110um_139um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(35,38) = Btd_110um_139um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(35,32) = Btd_110um_139um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(35,33) = Btd_110um_139um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(35,34) = Btd_110um_139um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(35,35) = Btd_110um_139um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(35,36) = Btd_110um_139um_Btd_110um_142um_Covar_Antarctic

   Btd_Covar(36,20) = Btd_110um_142um_Btd_110um_038um_Covar_Antarctic
   Btd_Covar(36,37) = Btd_110um_142um_Btd_110um_062um_Covar_Antarctic
   Btd_Covar(36,27) = Btd_110um_142um_Btd_110um_067um_Covar_Antarctic
   Btd_Covar(36,28) = Btd_110um_142um_Btd_110um_073um_Covar_Antarctic
   Btd_Covar(36,29) = Btd_110um_142um_Btd_110um_085um_Covar_Antarctic
   Btd_Covar(36,30) = Btd_110um_142um_Btd_110um_097um_Covar_Antarctic
   Btd_Covar(36,38) = Btd_110um_142um_Btd_110um_104um_Covar_Antarctic
   Btd_Covar(36,32) = Btd_110um_142um_Btd_110um_120um_Covar_Antarctic
   Btd_Covar(36,33) = Btd_110um_142um_Btd_110um_133um_Covar_Antarctic
   Btd_Covar(36,34) = Btd_110um_142um_Btd_110um_136um_Covar_Antarctic
   Btd_Covar(36,35) = Btd_110um_142um_Btd_110um_139um_Covar_Antarctic
   Btd_Covar(36,36) = Btd_110um_142um_Btd_110um_142um_Covar_Antarctic

 case default
   Bt_Covar(31) = Bt_110um_Bt_110um_Covar_All
   Bt_Covar(20) = Bt_110um_Btd_110um_038um_Covar_All
   Bt_Covar(37) = Bt_110um_Btd_110um_062um_Covar_All
   Bt_Covar(27) = Bt_110um_Btd_110um_067um_Covar_All
   Bt_Covar(28) = Bt_110um_Btd_110um_073um_Covar_All
   Bt_Covar(29) = Bt_110um_Btd_110um_085um_Covar_All
   Bt_Covar(30) = Bt_110um_Btd_110um_097um_Covar_All
   Bt_Covar(38) = Bt_110um_Btd_110um_104um_Covar_All
   Bt_Covar(32) = Bt_110um_Btd_110um_120um_Covar_All
   Bt_Covar(33) = Bt_110um_Btd_110um_133um_Covar_All
   Bt_Covar(34) = Bt_110um_Btd_110um_136um_Covar_All
   Bt_Covar(35) = Bt_110um_Btd_110um_139um_Covar_All
   Bt_Covar(36) = Bt_110um_Btd_110um_142um_Covar_All

   Btd_Covar(20,20) = Btd_110um_038um_Btd_110um_038um_Covar_All
   Btd_Covar(20,37) = Btd_110um_038um_Btd_110um_062um_Covar_All
   Btd_Covar(20,27) = Btd_110um_038um_Btd_110um_067um_Covar_All
   Btd_Covar(20,28) = Btd_110um_038um_Btd_110um_073um_Covar_All
   Btd_Covar(20,29) = Btd_110um_038um_Btd_110um_085um_Covar_All
   Btd_Covar(20,30) = Btd_110um_038um_Btd_110um_097um_Covar_All
   Btd_Covar(20,38) = Btd_110um_038um_Btd_110um_104um_Covar_All
   Btd_Covar(20,32) = Btd_110um_038um_Btd_110um_120um_Covar_All
   Btd_Covar(20,33) = Btd_110um_038um_Btd_110um_133um_Covar_All
   Btd_Covar(20,34) = Btd_110um_038um_Btd_110um_136um_Covar_All
   Btd_Covar(20,35) = Btd_110um_038um_Btd_110um_139um_Covar_All
   Btd_Covar(20,36) = Btd_110um_038um_Btd_110um_142um_Covar_All

   Btd_Covar(37,20) = Btd_110um_062um_Btd_110um_038um_Covar_All
   Btd_Covar(37,37) = Btd_110um_062um_Btd_110um_062um_Covar_All
   Btd_Covar(37,27) = Btd_110um_062um_Btd_110um_067um_Covar_All
   Btd_Covar(37,28) = Btd_110um_062um_Btd_110um_073um_Covar_All
   Btd_Covar(37,29) = Btd_110um_062um_Btd_110um_085um_Covar_All
   Btd_Covar(37,30) = Btd_110um_062um_Btd_110um_097um_Covar_All
   Btd_Covar(37,38) = Btd_110um_062um_Btd_110um_104um_Covar_All
   Btd_Covar(37,32) = Btd_110um_062um_Btd_110um_120um_Covar_All
   Btd_Covar(37,33) = Btd_110um_062um_Btd_110um_133um_Covar_All
   Btd_Covar(37,34) = Btd_110um_062um_Btd_110um_136um_Covar_All
   Btd_Covar(37,35) = Btd_110um_062um_Btd_110um_139um_Covar_All
   Btd_Covar(37,36) = Btd_110um_062um_Btd_110um_142um_Covar_All

   Btd_Covar(27,20) = Btd_110um_067um_Btd_110um_038um_Covar_All
   Btd_Covar(27,37) = Btd_110um_067um_Btd_110um_062um_Covar_All
   Btd_Covar(27,27) = Btd_110um_067um_Btd_110um_067um_Covar_All
   Btd_Covar(27,28) = Btd_110um_067um_Btd_110um_073um_Covar_All
   Btd_Covar(27,29) = Btd_110um_067um_Btd_110um_085um_Covar_All
   Btd_Covar(27,30) = Btd_110um_067um_Btd_110um_097um_Covar_All
   Btd_Covar(27,38) = Btd_110um_067um_Btd_110um_104um_Covar_All
   Btd_Covar(27,32) = Btd_110um_067um_Btd_110um_120um_Covar_All
   Btd_Covar(27,33) = Btd_110um_067um_Btd_110um_133um_Covar_All
   Btd_Covar(27,34) = Btd_110um_067um_Btd_110um_136um_Covar_All
   Btd_Covar(27,35) = Btd_110um_067um_Btd_110um_139um_Covar_All
   Btd_Covar(27,36) = Btd_110um_067um_Btd_110um_142um_Covar_All

   Btd_Covar(28,20) = Btd_110um_073um_Btd_110um_038um_Covar_All
   Btd_Covar(28,37) = Btd_110um_073um_Btd_110um_062um_Covar_All
   Btd_Covar(28,27) = Btd_110um_073um_Btd_110um_067um_Covar_All
   Btd_Covar(28,28) = Btd_110um_073um_Btd_110um_073um_Covar_All
   Btd_Covar(28,29) = Btd_110um_073um_Btd_110um_085um_Covar_All
   Btd_Covar(28,30) = Btd_110um_073um_Btd_110um_097um_Covar_All
   Btd_Covar(28,38) = Btd_110um_073um_Btd_110um_104um_Covar_All
   Btd_Covar(28,32) = Btd_110um_073um_Btd_110um_120um_Covar_All
   Btd_Covar(28,33) = Btd_110um_073um_Btd_110um_133um_Covar_All
   Btd_Covar(28,34) = Btd_110um_073um_Btd_110um_136um_Covar_All
   Btd_Covar(28,35) = Btd_110um_073um_Btd_110um_139um_Covar_All
   Btd_Covar(28,36) = Btd_110um_073um_Btd_110um_142um_Covar_All

   Btd_Covar(29,20) = Btd_110um_085um_Btd_110um_038um_Covar_All
   Btd_Covar(29,37) = Btd_110um_085um_Btd_110um_062um_Covar_All
   Btd_Covar(29,27) = Btd_110um_085um_Btd_110um_067um_Covar_All
   Btd_Covar(29,28) = Btd_110um_085um_Btd_110um_073um_Covar_All
   Btd_Covar(29,29) = Btd_110um_085um_Btd_110um_085um_Covar_All
   Btd_Covar(29,30) = Btd_110um_085um_Btd_110um_097um_Covar_All
   Btd_Covar(29,38) = Btd_110um_085um_Btd_110um_104um_Covar_All
   Btd_Covar(29,32) = Btd_110um_085um_Btd_110um_120um_Covar_All
   Btd_Covar(29,33) = Btd_110um_085um_Btd_110um_133um_Covar_All
   Btd_Covar(29,34) = Btd_110um_085um_Btd_110um_136um_Covar_All
   Btd_Covar(29,35) = Btd_110um_085um_Btd_110um_139um_Covar_All
   Btd_Covar(29,36) = Btd_110um_085um_Btd_110um_142um_Covar_All

   Btd_Covar(30,20) = Btd_110um_097um_Btd_110um_038um_Covar_All
   Btd_Covar(30,37) = Btd_110um_097um_Btd_110um_062um_Covar_All
   Btd_Covar(30,27) = Btd_110um_097um_Btd_110um_067um_Covar_All
   Btd_Covar(30,28) = Btd_110um_097um_Btd_110um_073um_Covar_All
   Btd_Covar(30,29) = Btd_110um_097um_Btd_110um_085um_Covar_All
   Btd_Covar(30,30) = Btd_110um_097um_Btd_110um_097um_Covar_All
   Btd_Covar(30,38) = Btd_110um_097um_Btd_110um_104um_Covar_All
   Btd_Covar(30,32) = Btd_110um_097um_Btd_110um_120um_Covar_All
   Btd_Covar(30,33) = Btd_110um_097um_Btd_110um_133um_Covar_All
   Btd_Covar(30,34) = Btd_110um_097um_Btd_110um_136um_Covar_All
   Btd_Covar(30,35) = Btd_110um_097um_Btd_110um_139um_Covar_All
   Btd_Covar(30,36) = Btd_110um_097um_Btd_110um_142um_Covar_All

   Btd_Covar(38,20) = Btd_110um_104um_Btd_110um_038um_Covar_All
   Btd_Covar(38,37) = Btd_110um_104um_Btd_110um_062um_Covar_All
   Btd_Covar(38,27) = Btd_110um_104um_Btd_110um_067um_Covar_All
   Btd_Covar(38,28) = Btd_110um_104um_Btd_110um_073um_Covar_All
   Btd_Covar(38,29) = Btd_110um_104um_Btd_110um_085um_Covar_All
   Btd_Covar(38,30) = Btd_110um_104um_Btd_110um_097um_Covar_All
   Btd_Covar(38,38) = Btd_110um_104um_Btd_110um_104um_Covar_All
   Btd_Covar(38,32) = Btd_110um_104um_Btd_110um_120um_Covar_All
   Btd_Covar(38,33) = Btd_110um_104um_Btd_110um_133um_Covar_All
   Btd_Covar(38,34) = Btd_110um_104um_Btd_110um_136um_Covar_All
   Btd_Covar(38,35) = Btd_110um_104um_Btd_110um_139um_Covar_All
   Btd_Covar(38,36) = Btd_110um_104um_Btd_110um_142um_Covar_All

   Btd_Covar(32,20) = Btd_110um_120um_Btd_110um_038um_Covar_All
   Btd_Covar(32,37) = Btd_110um_120um_Btd_110um_062um_Covar_All
   Btd_Covar(32,27) = Btd_110um_120um_Btd_110um_067um_Covar_All
   Btd_Covar(32,28) = Btd_110um_120um_Btd_110um_073um_Covar_All
   Btd_Covar(32,29) = Btd_110um_120um_Btd_110um_085um_Covar_All
   Btd_Covar(32,30) = Btd_110um_120um_Btd_110um_097um_Covar_All
   Btd_Covar(32,38) = Btd_110um_120um_Btd_110um_104um_Covar_All
   Btd_Covar(32,32) = Btd_110um_120um_Btd_110um_120um_Covar_All
   Btd_Covar(32,33) = Btd_110um_120um_Btd_110um_133um_Covar_All
   Btd_Covar(32,34) = Btd_110um_120um_Btd_110um_136um_Covar_All
   Btd_Covar(32,35) = Btd_110um_120um_Btd_110um_139um_Covar_All
   Btd_Covar(32,36) = Btd_110um_120um_Btd_110um_142um_Covar_All

   Btd_Covar(33,20) = Btd_110um_133um_Btd_110um_038um_Covar_All
   Btd_Covar(33,37) = Btd_110um_133um_Btd_110um_062um_Covar_All
   Btd_Covar(33,27) = Btd_110um_133um_Btd_110um_067um_Covar_All
   Btd_Covar(33,28) = Btd_110um_133um_Btd_110um_073um_Covar_All
   Btd_Covar(33,29) = Btd_110um_133um_Btd_110um_085um_Covar_All
   Btd_Covar(33,30) = Btd_110um_133um_Btd_110um_097um_Covar_All
   Btd_Covar(33,38) = Btd_110um_133um_Btd_110um_104um_Covar_All
   Btd_Covar(33,32) = Btd_110um_133um_Btd_110um_120um_Covar_All
   Btd_Covar(33,33) = Btd_110um_133um_Btd_110um_133um_Covar_All
   Btd_Covar(33,34) = Btd_110um_133um_Btd_110um_136um_Covar_All
   Btd_Covar(33,35) = Btd_110um_133um_Btd_110um_139um_Covar_All
   Btd_Covar(33,36) = Btd_110um_133um_Btd_110um_142um_Covar_All

   Btd_Covar(34,20) = Btd_110um_136um_Btd_110um_038um_Covar_All
   Btd_Covar(34,37) = Btd_110um_136um_Btd_110um_062um_Covar_All
   Btd_Covar(34,27) = Btd_110um_136um_Btd_110um_067um_Covar_All
   Btd_Covar(34,28) = Btd_110um_136um_Btd_110um_073um_Covar_All
   Btd_Covar(34,29) = Btd_110um_136um_Btd_110um_085um_Covar_All
   Btd_Covar(34,30) = Btd_110um_136um_Btd_110um_097um_Covar_All
   Btd_Covar(34,38) = Btd_110um_136um_Btd_110um_104um_Covar_All
   Btd_Covar(34,32) = Btd_110um_136um_Btd_110um_120um_Covar_All
   Btd_Covar(34,33) = Btd_110um_136um_Btd_110um_133um_Covar_All
   Btd_Covar(34,34) = Btd_110um_136um_Btd_110um_136um_Covar_All
   Btd_Covar(34,35) = Btd_110um_136um_Btd_110um_139um_Covar_All
   Btd_Covar(34,36) = Btd_110um_136um_Btd_110um_142um_Covar_All

   Btd_Covar(35,20) = Btd_110um_139um_Btd_110um_038um_Covar_All
   Btd_Covar(35,37) = Btd_110um_139um_Btd_110um_062um_Covar_All
   Btd_Covar(35,27) = Btd_110um_139um_Btd_110um_067um_Covar_All
   Btd_Covar(35,28) = Btd_110um_139um_Btd_110um_073um_Covar_All
   Btd_Covar(35,29) = Btd_110um_139um_Btd_110um_085um_Covar_All
   Btd_Covar(35,30) = Btd_110um_139um_Btd_110um_097um_Covar_All
   Btd_Covar(35,38) = Btd_110um_139um_Btd_110um_104um_Covar_All
   Btd_Covar(35,32) = Btd_110um_139um_Btd_110um_120um_Covar_All
   Btd_Covar(35,33) = Btd_110um_139um_Btd_110um_133um_Covar_All
   Btd_Covar(35,34) = Btd_110um_139um_Btd_110um_136um_Covar_All
   Btd_Covar(35,35) = Btd_110um_139um_Btd_110um_139um_Covar_All
   Btd_Covar(35,36) = Btd_110um_139um_Btd_110um_142um_Covar_All

   Btd_Covar(36,20) = Btd_110um_142um_Btd_110um_038um_Covar_All
   Btd_Covar(36,37) = Btd_110um_142um_Btd_110um_062um_Covar_All
   Btd_Covar(36,27) = Btd_110um_142um_Btd_110um_067um_Covar_All
   Btd_Covar(36,28) = Btd_110um_142um_Btd_110um_073um_Covar_All
   Btd_Covar(36,29) = Btd_110um_142um_Btd_110um_085um_Covar_All
   Btd_Covar(36,30) = Btd_110um_142um_Btd_110um_097um_Covar_All
   Btd_Covar(36,38) = Btd_110um_142um_Btd_110um_104um_Covar_All
   Btd_Covar(36,32) = Btd_110um_142um_Btd_110um_120um_Covar_All
   Btd_Covar(36,33) = Btd_110um_142um_Btd_110um_133um_Covar_All
   Btd_Covar(36,34) = Btd_110um_142um_Btd_110um_136um_Covar_All
   Btd_Covar(36,35) = Btd_110um_142um_Btd_110um_139um_Covar_All
   Btd_Covar(36,36) = Btd_110um_142um_Btd_110um_142um_Covar_All

 end select

end subroutine SET_CLEAR_SKY_COVARIANCE_TERMS

!-------------------------------------------------------
!  Linear in Optical Depth Emission Routine
!-------------------------------------------------------
function Linear_In_Opd_Emission(Emiss,B_Base,B_Top) result(Cloud_Emission)
   real, intent(in):: Emiss, B_Base, B_top
   real:: Opd, Bd, Linear_Term, Cloud_Emission 
   Opd = -1.0*alog(1.0-Emiss)
   Opd = max(0.01,min(10.0,Opd))
   Bd = (B_Base - B_Top) / Opd
   Linear_Term = exp(-Opd)*(1+Opd)-1
   Cloud_Emission = Emiss * B_Top - Bd * Linear_Term 
end function Linear_In_Opd_Emission

!---------------------------------------------------------------------------------------------
! Compute clear-sky terms needed in forward model
!---------------------------------------------------------------------------------------------
subroutine  COMPUTE_CLEAR_SKY_TERMS(Acha_Mode_Flag, Zc, Zs, Ts, Hght_Prof,  &
                               Chan_Idx_038um, Chan_Idx_062um, Chan_Idx_067um, &
                               Chan_Idx_073um, Chan_Idx_085um, Chan_Idx_097um, Chan_Idx_104um, &
                               Chan_Idx_110um, Chan_Idx_120um, Chan_Idx_133um, &
                               Chan_Idx_136um, Chan_Idx_139um, Chan_Idx_142um, &
                               Atm_Rad_Prof_038um, Atm_Trans_Prof_038um, &
                               Atm_Rad_Prof_062um, Atm_Trans_Prof_062um, &
                               Atm_Rad_Prof_067um, Atm_Trans_Prof_067um, &
                               Atm_Rad_Prof_073um, Atm_Trans_Prof_073um, &
                               Atm_Rad_Prof_085um, Atm_Trans_Prof_085um, &
                               Atm_Rad_Prof_097um, Atm_Trans_Prof_097um, &
                               Atm_Rad_Prof_104um, Atm_Trans_Prof_104um, &
                               Atm_Rad_Prof_110um, Atm_Trans_Prof_110um, &
                               Atm_Rad_Prof_120um, Atm_Trans_Prof_120um, &
                               Atm_Rad_Prof_133um, Atm_Trans_Prof_133um, &
                               Atm_Rad_Prof_136um, Atm_Trans_Prof_136um, &
                               Atm_Rad_Prof_139um, Atm_Trans_Prof_139um, &
                               Atm_Rad_Prof_142um, Atm_Trans_Prof_142um, &
                               Emiss_Sfc_038um, Emiss_Sfc_062um, Emiss_Sfc_067um, &
                               Emiss_Sfc_073um, Emiss_Sfc_085um, Emiss_Sfc_097um, Emiss_Sfc_104um, &
                               Emiss_Sfc_110um, Emiss_Sfc_120um, Emiss_Sfc_133um, &
                               Emiss_Sfc_136um, Emiss_Sfc_139um, Emiss_Sfc_142um, &
                               Rad_Ac_038um, Trans_Ac_038um, Trans_Bc_038um, Rad_Clear_038um, &
                               Rad_Ac_062um, Trans_Ac_062um, Trans_Bc_062um, Rad_Clear_062um, &
                               Rad_Ac_067um, Trans_Ac_067um, Trans_Bc_067um, Rad_Clear_067um, &
                               Rad_Ac_073um, Trans_Ac_073um, Trans_Bc_073um, Rad_Clear_073um, &
                               Rad_Ac_085um, Trans_Ac_085um, Trans_Bc_085um, Rad_Clear_085um, &
                               Rad_Ac_097um, Trans_Ac_097um, Trans_Bc_097um, Rad_Clear_097um, &
                               Rad_Ac_104um, Trans_Ac_104um, Trans_Bc_104um, Rad_Clear_104um, &
                               Rad_Ac_110um, Trans_Ac_110um, Trans_Bc_110um, Rad_Clear_110um, &
                               Rad_Ac_120um, Trans_Ac_120um, Trans_Bc_120um, Rad_Clear_120um, &
                               Rad_Ac_133um, Trans_Ac_133um, Trans_Bc_133um, Rad_Clear_133um, &
                               Rad_Ac_136um, Trans_Ac_136um, Trans_Bc_136um, Rad_Clear_136um, &
                               Rad_Ac_139um, Trans_Ac_139um, Trans_Bc_139um, Rad_Clear_139um, &
                               Rad_Ac_142um, Trans_Ac_142um, Trans_Bc_142um, Rad_Clear_142um)

  character(len=*), intent(in):: Acha_Mode_Flag
  real, intent(in):: Zc, Zs, Ts
  real, intent(in), dimension(:):: Hght_Prof
  integer, intent(in):: Chan_Idx_038um, Chan_Idx_062um, Chan_Idx_067um, Chan_Idx_073um, &
                        Chan_Idx_085um, Chan_Idx_097um, Chan_Idx_104um, Chan_Idx_110um, Chan_Idx_120um, &
                        Chan_Idx_133um, Chan_Idx_136um, Chan_Idx_139um, Chan_Idx_142um
  real, intent(in), dimension(:):: Atm_Rad_Prof_038um, Atm_Trans_Prof_038um
  real, intent(in), dimension(:):: Atm_Rad_Prof_062um, Atm_Trans_Prof_062um
  real, intent(in), dimension(:):: Atm_Rad_Prof_067um, Atm_Trans_Prof_067um
  real, intent(in), dimension(:):: Atm_Rad_Prof_073um, Atm_Trans_Prof_073um
  real, intent(in), dimension(:):: Atm_Rad_Prof_085um, Atm_Trans_Prof_085um
  real, intent(in), dimension(:):: Atm_Rad_Prof_097um, Atm_Trans_Prof_097um
  real, intent(in), dimension(:):: Atm_Rad_Prof_104um, Atm_Trans_Prof_104um
  real, intent(in), dimension(:):: Atm_Rad_Prof_110um, Atm_Trans_Prof_110um
  real, intent(in), dimension(:):: Atm_Rad_Prof_120um, Atm_Trans_Prof_120um
  real, intent(in), dimension(:):: Atm_Rad_Prof_133um, Atm_Trans_Prof_133um
  real, intent(in), dimension(:):: Atm_Rad_Prof_136um, Atm_Trans_Prof_136um
  real, intent(in), dimension(:):: Atm_Rad_Prof_139um, Atm_Trans_Prof_139um
  real, intent(in), dimension(:):: Atm_Rad_Prof_142um, Atm_Trans_Prof_142um
  real, intent(in):: Emiss_Sfc_038um, Emiss_Sfc_062um, Emiss_Sfc_067um, &
                     Emiss_Sfc_073um, Emiss_Sfc_085um, Emiss_Sfc_097um, Emiss_Sfc_104um, &
                     Emiss_Sfc_110um, Emiss_Sfc_120um, Emiss_Sfc_133um, &
                     Emiss_Sfc_136um, Emiss_Sfc_139um, Emiss_Sfc_142um
  real, intent(out):: Rad_Ac_038um, Trans_Ac_038um, Trans_Bc_038um, Rad_Clear_038um
  real, intent(out):: Rad_Ac_062um, Trans_Ac_062um, Trans_Bc_062um, Rad_Clear_062um
  real, intent(out):: Rad_Ac_067um, Trans_Ac_067um, Trans_Bc_067um, Rad_Clear_067um
  real, intent(out):: Rad_Ac_073um, Trans_Ac_073um, Trans_Bc_073um, Rad_Clear_073um
  real, intent(out):: Rad_Ac_085um, Trans_Ac_085um, Trans_Bc_085um, Rad_Clear_085um
  real, intent(out):: Rad_Ac_097um, Trans_Ac_097um, Trans_Bc_097um, Rad_Clear_097um
  real, intent(out):: Rad_Ac_104um, Trans_Ac_104um, Trans_Bc_104um, Rad_Clear_104um
  real, intent(out):: Rad_Ac_110um, Trans_Ac_110um, Trans_Bc_110um, Rad_Clear_110um
  real, intent(out):: Rad_Ac_120um, Trans_Ac_120um, Trans_Bc_120um, Rad_Clear_120um
  real, intent(out):: Rad_Ac_133um, Trans_Ac_133um, Trans_Bc_133um, Rad_Clear_133um
  real, intent(out):: Rad_Ac_136um, Trans_Ac_136um, Trans_Bc_136um, Rad_Clear_136um
  real, intent(out):: Rad_Ac_139um, Trans_Ac_139um, Trans_Bc_139um, Rad_Clear_139um
  real, intent(out):: Rad_Ac_142um, Trans_Ac_142um, Trans_Bc_142um, Rad_Clear_142um

  !--- If GOES-17 11 um data is bad, 11 um variables replaced with 10.4 um data.
  !--- compute 110um radiative transfer terms
  call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_110um, &
                                  Atm_Rad_Prof_110um,Atm_Trans_Prof_110um,&
                                  Emiss_Sfc_110um, Rad_Ac_110um,Trans_Ac_110um, &
                                  Trans_Bc_110um,Rad_Clear_110um)

  !--- compute 3.75um radiative transfer terms
  if (index(Acha_Mode_Flag,'038') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_038um, &
                                     Atm_Rad_Prof_038um,Atm_Trans_Prof_038um,&
                                     Emiss_Sfc_038um, Rad_Ac_038um,Trans_Ac_038um, &
                                     Trans_Bc_038um,Rad_Clear_038um)

  endif

  !--- 6.2um clear radiative transfer terms
  if (index(Acha_Mode_Flag,'062') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_062um, &
                                     Atm_Rad_Prof_062um,Atm_Trans_Prof_062um,&
                                     Emiss_Sfc_062um, Rad_Ac_062um,Trans_Ac_062um, &
                                     Trans_Bc_062um,Rad_Clear_062um)

  endif

  !--- 6.7um clear radiative transfer terms
  if (index(Acha_Mode_Flag,'067') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_067um, &
                                     Atm_Rad_Prof_067um,Atm_Trans_Prof_067um,&
                                     Emiss_Sfc_067um, Rad_Ac_067um,Trans_Ac_067um, &
                                     Trans_Bc_067um,Rad_Clear_067um)
  endif

  !--- 7.3um clear radiative transfer terms
  if (index(Acha_Mode_Flag,'073') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_073um, &
                                     Atm_Rad_Prof_073um,Atm_Trans_Prof_073um,&
                                     Emiss_Sfc_073um, Rad_Ac_073um,Trans_Ac_073um, &
                                     Trans_Bc_073um,Rad_Clear_073um)
  endif

  !--- 8.5um clear radiative transfer terms
  if (index(Acha_Mode_Flag,'085') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_085um, &
                                     Atm_Rad_Prof_085um,Atm_Trans_Prof_085um,&
                                     Emiss_Sfc_085um, Rad_Ac_085um,Trans_Ac_085um, &
                                     Trans_Bc_085um,Rad_Clear_085um)
  endif

  !--- 9.7um clear radiative transfer terms
  if (index(Acha_Mode_Flag,'097') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_097um, &
                                     Atm_Rad_Prof_097um,Atm_Trans_Prof_097um,&
                                     Emiss_Sfc_097um, Rad_Ac_097um,Trans_Ac_097um, &
                                     Trans_Bc_097um,Rad_Clear_097um)
  endif

  !--- 10.4um clear radiative transfer terms
  if (index(Acha_Mode_Flag,'104') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_104um, &
                                     Atm_Rad_Prof_104um,Atm_Trans_Prof_104um,&
                                     Emiss_Sfc_104um, Rad_Ac_104um,Trans_Ac_104um, &
                                     Trans_Bc_104um,Rad_Clear_104um)
  endif

  !--- compute 120um radiative transfer terms
  if (index(Acha_Mode_Flag,'120') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_120um, &
                                     Atm_Rad_Prof_120um,Atm_Trans_Prof_120um,&
                                     Emiss_Sfc_120um, Rad_Ac_120um,Trans_Ac_120um, &
                                     Trans_Bc_120um,Rad_Clear_120um)
  endif

  !--- 13.3um clear radiative transfer terms
  if (index(Acha_Mode_Flag,'133') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_133um, &
                                     Atm_Rad_Prof_133um,Atm_Trans_Prof_133um,&
                                     Emiss_Sfc_133um, Rad_Ac_133um,Trans_Ac_133um, &
                                     Trans_Bc_133um,Rad_Clear_133um)
  endif

  !--- 13.6um clear radiative transfer terms
  if (index(Acha_Mode_Flag,'136') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_136um, &
                                     Atm_Rad_Prof_136um,Atm_Trans_Prof_136um,&
                                     Emiss_Sfc_136um, Rad_Ac_136um,Trans_Ac_136um, &
                                     Trans_Bc_136um,Rad_Clear_136um)
  endif

  !--- 13.9um clear radiative transfer terms
  if (index(Acha_Mode_Flag,'139') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_139um, &
                                     Atm_Rad_Prof_139um,Atm_Trans_Prof_139um,&
                                     Emiss_Sfc_139um, Rad_Ac_139um,Trans_Ac_139um, &
                                     Trans_Bc_139um,Rad_Clear_139um)
  endif

  !--- 14.2um clear radiative transfer terms
  if (index(Acha_Mode_Flag,'142') > 0) then

     call CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx_142um, &
                                     Atm_Rad_Prof_142um,Atm_Trans_Prof_142um,&
                                     Emiss_Sfc_142um, Rad_Ac_142um,Trans_Ac_142um, &
                                     Trans_Bc_142um,Rad_Clear_142um)
  endif

end subroutine  COMPUTE_CLEAR_SKY_TERMS
!------------------------------------------------------------------------------------------
! routine to compute the terms needed in COMPUTE_CLEAR_SKY_TERMS
!------------------------------------------------------------------------------------------
subroutine CLEAR_SKY_INTERNAL_ROUTINE(Zc,Zs,Ts,Hght_Prof,Chan_Idx, &
                                      Atm_Rad_Prof,Atm_Trans_Prof,Emiss_Sfc,&
                                      Rad_Ac,Trans_Ac,Trans_Bc,Rad_Clear)
  real, intent(in):: Zc, Zs, Ts
  real, intent(in), dimension(:):: Hght_Prof
  integer, intent(in):: Chan_Idx
  real, intent(in), dimension(:):: Atm_Rad_Prof, Atm_Trans_Prof
  real, intent(in):: Emiss_Sfc
  real, intent(out) :: Rad_Ac,Trans_Ac,Trans_Bc,Rad_Clear
  real :: Rad_Atm, Trans_Atm, Bs

  Rad_Ac = GENERIC_PROFILE_INTERPOLATION(Zc,Hght_Prof,Atm_Rad_Prof)

  Trans_Ac = GENERIC_PROFILE_INTERPOLATION(Zc,Hght_Prof,Atm_Trans_Prof)

  Trans_Bc = GENERIC_PROFILE_INTERPOLATION(Zs,Hght_Prof,Atm_Trans_Prof) 

  if (Trans_Ac > epsilon(Trans_Ac)) then
     Trans_Bc = Trans_Bc / Trans_Ac
  endif

  Rad_Atm = GENERIC_PROFILE_INTERPOLATION(Zs,Hght_Prof,Atm_Rad_Prof)
 
  Trans_Atm = GENERIC_PROFILE_INTERPOLATION(Zs,Hght_Prof,Atm_Trans_Prof)

  Bs = PLANCK_RAD_FAST(Chan_Idx,Ts)

  Rad_Clear = Rad_Atm + Trans_Atm*Emiss_Sfc*Bs

end subroutine CLEAR_SKY_INTERNAL_ROUTINE
!-------------------------------------------------------------------------------------
! forward model for a brightness temperature
!
! Channel X refers to the 11 micron channel
! input
!  chan_x = channel number of channel x
!  f_x = forward model estimate of the 11 micron brightness temperature
!  ec_x = cloud emissivity at 11 micron
!
!--------------------------------------------------------------------------------------
subroutine BT_FM(chan_x,tc,ec_x,ts,tc_base,es_x, &
                 rad_ac_x, trans_ac_x, trans_bc_x, rad_clear_x, & 
                 f,df_dtc,df_dec,df_dbeta,df_dTs, df_dalpha, &
                 t_eff)

     integer, intent(in):: chan_x
     real, intent(in):: ec_x, es_x
     real, intent(in):: tc,ts,tc_base
     real, intent(in):: rad_ac_x, trans_ac_x, trans_bc_x, rad_clear_x
     real, intent(out):: f,df_dtc,df_dec,df_dbeta,df_dTs, df_dalpha
     real, intent(out),optional:: t_eff
     real:: trans_x, rad_x, cloud_emission
     real:: bc_x, bc_base_x, bs_x, db_dt_x, db_dtc_x, db_dts_x
     real:: alpha_x
     real:: b_eff

     !--- planck function terms
     bc_x = PLANCK_RAD_FAST(chan_x, tc, db_dt = db_dtc_x) 
     bc_base_x = PLANCK_RAD_FAST(chan_x, tc_base) 
     bs_x = PLANCK_RAD_FAST(chan_x, ts, dB_dT = db_dts_x)

     trans_x = 1.0 - ec_x

     cloud_emission = Linear_In_Opd_Emission(ec_x,bc_base_x,bc_x)

     rad_x = ec_x*rad_ac_x + trans_ac_x * cloud_emission + trans_x * rad_clear_x

     !--- forward model term
     f = PLANCK_TEMP_FAST(chan_x, rad_x, db_dt = db_dt_x)

     !--- effective temperature
     if (present(t_eff)) then 
       b_eff = cloud_emission / ec_x
       t_eff = PLANCK_TEMP_FAST(chan_x, b_eff)
!      print *, "Test ", chan_x, ec_x, tc, cloud_emission, b_eff, t_eff
     endif  


!    print *, "in BT FM ", tc, tc_base, bc_x, bc_base_x, cloud_emission, ec_x * bc_x, t_eff

     !--- kernel matrix terms
     alpha_x = rad_ac_x + trans_ac_x*bc_x - rad_clear_x

     df_dtc = (trans_ac_x * ec_x * db_dtc_x)/db_dt_x

     df_dec = alpha_x / db_dt_x

     df_dbeta = 0.0

     df_dts = (trans_ac_x * trans_x * trans_bc_x * db_dts_x * es_x) / db_dt_x

     df_dalpha = 0.0

end subroutine BT_FM

!-------------------------------------------------------------------------------------
! forward model for a brightness temperature difference
!
! Channel X refers to the 11 micron channel
! input
!  chan_x = channel number of channel x
!  f_x = forward model estimate of the 11 micron brightness temperature
!  ec_x = cloud emissivity at 11 micron
!  beta_x_12 = beta of 11 micron and 12 micron (the reference value)
!
! output
!  f = the forward model estimate of the btd of x - y
!  df_dtc = the derivative of f wrt to cloud temperature
!  df_dec = the derivative of f wrt to cloud emissivity
!  df_dbeta = the derivative of f wrt to cloud beta 
!  df_dts = the derivative of f wrt to surface temperature
!  df_dalpha = the derivative of f wrt to ice_fraction (alpha)
!-------------------------------------------------------------------------------------
subroutine BTD_FM(chan_y,  &
                  beta_xy_coef_water, &
                  beta_xy_coef_ice, &
                  tc,ec_x, beta_x_12,ts,tc_base,es_y, alpha,&
                  f_x, df_x_dtc, df_x_dec, df_x_dts, &
                  rad_ac_y, trans_ac_y, trans_bc_y, rad_clear_y, & 
                  f,df_dtc,df_dec,df_dbeta,df_dTs,df_dalpha,ec_y)

     integer, intent(in):: chan_y
     real, dimension(0:), intent(in):: beta_xy_coef_water
     real, dimension(0:), intent(in):: beta_xy_coef_ice
     real, intent(in):: ec_x,beta_x_12
     real, intent(in):: tc,ts,tc_base,es_y, alpha
     real, intent(in):: f_x, df_x_dtc, df_x_dec, df_x_dts
     real, intent(in):: rad_ac_y, trans_ac_y, trans_bc_y, rad_clear_y
     real, intent(out):: f,df_dtc,df_dec,df_dbeta,df_dTs,df_dalpha,ec_y
     real:: trans_y, rad_y, cloud_emission
     real:: bc_y, bc_base_y, bs_y, db_dt_y, db_dtc_y, db_dts_y
     real:: beta_xy, dbeta_xy_dbeta_x_12, dec_y_dec_x, dalpha_dbeta
     real:: alpha_y

     !--- planck function terms
     bc_y = PLANCK_RAD_FAST(chan_y, tc, db_dt = db_dtc_y) 
     bc_base_y = PLANCK_RAD_FAST(chan_y, tc_base) 
     bs_y = PLANCK_RAD_FAST(chan_y, ts, dB_dT = db_dts_y)

     !--- intermediate terms
     call COMPUTE_BETA_AND_DERIVATIVE(beta_degree_water, beta_xy_coef_water, beta_degree_ice, beta_xy_coef_ice, &
                                      alpha, beta_x_12, beta_xy, dbeta_xy_dbeta_x_12, dalpha_dbeta)

     dec_y_dec_x = beta_xy * (1.0 - ec_x)**(beta_xy - 1.0)

     ec_y = 1.0 - (1.0 - ec_x)**beta_xy

     trans_y = 1.0 - ec_y

     cloud_emission = Linear_In_Opd_Emission(ec_y,bc_base_y,bc_y)

     rad_y = ec_y*rad_ac_y + trans_ac_y * cloud_emission + trans_y * rad_clear_y

     !--- forward model term
     f = f_x - PLANCK_TEMP_FAST(chan_y,rad_y,db_dt = db_dt_y)

     !--- kernel matrix terms
     alpha_y = rad_ac_y + trans_ac_y*bc_y - rad_clear_y

     df_dtc = df_x_dtc - (trans_ac_y * ec_y * db_dtc_y)/db_dt_y

     df_dec = df_x_dec - (alpha_y * dec_y_dec_x)/db_dt_y

     df_dbeta = (alpha_y) * alog(1.0-ec_x) * (1.0-ec_y) * dbeta_xy_dbeta_x_12 / db_dt_y

     df_dts = df_x_dts - (trans_ac_y * trans_y * trans_bc_y * db_dts_y * es_y) / db_dt_y

     df_dalpha = df_dbeta / (dalpha_dbeta)

end subroutine BTD_FM

!------------------------------------------------------------------------------
! compute a channel pairs channel beta and derivative
!
! input: beta_xy_coef_water - beta coefficients for water clouds
!        beta_xy_coef_ice - beta coefficients for ice clouds
!        beta_degree_water - degree of the polynomial phase for water
!        beta_degree_ice - degree of the polynomial phase for ice
!        alpha - ice cloud fraction
!        beta_x_12 - the beta value for 11 and 12 micron
!
! output:
!        beta_xy - the beta value for this channel pair
!        dbeta_xy_dbeta_x_12 - the derivative of beta value for this channel 
!                              pair to the beta_x_12
!------------------------------------------------------------------------------
subroutine COMPUTE_BETA_AND_DERIVATIVE(beta_degree_water, &
                                       beta_xy_coef_water, &
                                       beta_degree_ice, &
                                       beta_xy_coef_ice, &
                                       alpha, beta_x_12, &
                                       beta_xy, &
                                       dbeta_xy_dbeta_x_12, &
                                       dalpha_dbeta)

   integer, intent(in):: beta_degree_water, beta_degree_ice
   real, dimension(0:), intent(in):: beta_xy_coef_water
   real, dimension(0:), intent(in):: beta_xy_coef_ice
   real, intent(in):: alpha
   real, intent(in):: beta_x_12
   real, intent(out):: beta_xy
   real, intent(out):: dbeta_xy_dbeta_x_12
   real, intent(out):: dalpha_dbeta

   real:: beta_xy_water, beta_xy_ice
   real:: dbeta_xy_dbeta_x_12_water
   real:: dbeta_xy_dbeta_x_12_ice
   real:: dbeta

   integer:: i

   !----------------------------------------------------------------------
   ! water
   !----------------------------------------------------------------------
   beta_xy_water = beta_xy_coef_water(0)
   dbeta_xy_dbeta_x_12_water = 0.0

   do i = 1, beta_degree_water

      beta_xy_water = beta_xy_water + beta_xy_coef_water(i) * (beta_x_12-1.0)**(i)

      dbeta_xy_dbeta_x_12_water = dbeta_xy_dbeta_x_12_water + &
                                  beta_xy_coef_water(i) * (i) * (beta_x_12-1.0)**(i-1)
   enddo

   !----------------------------------------------------------------------
   ! ice
   !----------------------------------------------------------------------
   beta_xy_ice = beta_xy_coef_ice(0)
   dbeta_xy_dbeta_x_12_ice = 0.0

   do i = 1, beta_degree_ice

      beta_xy_ice = beta_xy_ice + beta_xy_coef_ice(i) * (beta_x_12-1.0)**(i)

      dbeta_xy_dbeta_x_12_ice = dbeta_xy_dbeta_x_12_ice + &
                                beta_xy_coef_ice(i) * (i) * (beta_x_12-1.0)**(i-1)
   enddo

   dbeta = beta_xy_ice - beta_xy_water
   if (abs(dbeta) < epsilon(beta_xy_ice)) dbeta = 0.01
   dalpha_dbeta = 1.0/dbeta

   !----------------------------------------------------------------------
   ! combine
   !----------------------------------------------------------------------
   beta_xy = (1.0-alpha)*beta_xy_water + alpha*beta_xy_ice

   dbeta_xy_dbeta_x_12 = (1.0-alpha)*dbeta_xy_dbeta_x_12_water +  &
                         alpha*dbeta_xy_dbeta_x_12_ice

end subroutine COMPUTE_BETA_AND_DERIVATIVE

!----------------------------------------------------------------------
! Compute Sy based on the clear-sky error covariance calcuLations.
! Using Andy's simpler expression
!
! This assumes that 
! Acha_Mode_Flag: 1=110um,2=11+3.75,3=11+6.7um,4=11+120um,5=11+13.3um, $
!                 6=11+12+8.5um
!                 7=11+6.7+120um,8=11+6.7+13.3um,
!                 9=11+12+13.3um
!                 10=11+12+13.3um(pseudo)
!                 11=11+12+8.5+13.3um
!                 12=11+12+8.5+6.7um
!                 13=11+12+8.5+13.3+6.7um
!
! Input:
! Emiss_Vector = a vector of emissivities in each channel. 
! Acha_Mode_Flag: 1=110um,2=11+6.7um,3=11+120um,4=11+13.3um,5=8.5+11+120um
!                 6=11+6.7+120um,7=11+6.7+13.3um,8=11+12+13.3um,9=11+8.5+6.7
! Sfc_Type_Forward_Model = the surface type used for covariance calcs 
! y_variance = the variance computed a 3x3 array for each element of y
!
! Output:
! Sy = error covariance matrix
!
!  Sy(i,i) = Cal_Err_y(i)^2 + Spatial_Variance_y(i) + (1-emiss_i)^2*y(i)_covar
!  Sy(i,j) = (1-emiss_i)*(1-emiss_j)*y(i)_y(x)_covar
!----------------------------------------------------------------------
subroutine COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE(   &
                                             Emiss_Vector,  &
                                             Acha_Mode_Flag, &
                                             y_variance, &
                                             Sy) 


 real(kind=real4), intent(in), dimension(:):: Emiss_Vector
 character(len=*), intent(in):: Acha_Mode_Flag
 real(kind=real4), intent(in), dimension(:):: y_variance
 real(kind=real4), intent(out), dimension(:,:):: Sy
 real(kind=real4), dimension(size(y_variance)):: Sub_Pixel_Uncer
 real(kind=real4):: Emiss_110um
 real(kind=real4):: Trans2   

 integer:: i,j,n
 real:: Infinity

 Emiss_110um = min(1.0,max(0.0,Emiss_Vector(1)))

 Trans2 = (1.0 - Emiss_110um)**2  !cloud transmission squared
 Trans2 = max(Trans2,0.25)       !do not let this go to zero

 n = size(y_variance)
 !----------------------------------------------------------------
 !--- modify y_variance to represent a sub-pixel uncertainty 
 !--- assume that all of standard deviation is due to sub-pixel
 !--- heterogeneity and that this is a good estimate of the
 !--- forward model error due to sub-pixel heterogeneity
 !----------------------------------------------------------------
 do i = 1, n
    Sub_Pixel_Uncer(i) = y_variance(i) 
 enddo

 !---- compute the Sy matrix
 do i = 1, n
   do j = 1, n
      Sy(i,j) = Trans2*Btd_Covar(Chan_Idx_y(i),Chan_Idx_y(j)) 
      if (i == 1) then
         Sy(i,j) = Trans2*Bt_Covar(Chan_Idx_y(j)) 
      endif
      if (j == 1) then
         Sy(i,j) = Trans2*Bt_Covar(Chan_Idx_y(i)) 
      endif
   enddo
 enddo

 !-- add in terms for diagnonal elements to Sy
 do i = 1,n
   Sy(i,i) = Sy(i,i) + Cal_Uncer(Chan_Idx_y(i))**2 + Sub_Pixel_Uncer(i)
 enddo

 !-- add in terms for diagnonal elements for cloud btd error
 Sy(1,1) = Sy(1,1) + (Emiss_Vector(1)*Cloud_BT_Uncer)**2
 do i = 2,n
   Sy(i,i) = Sy(i,i) + (Emiss_Vector(i)*Cloud_BTD_Uncer(Chan_Idx_y(i)))**2
 enddo

 !---- check Sy
! Infinity = huge(Trans2)
! do i = 1,n
!   do j = 1,n
!     if (Sy(i,j) /= Sy(i,j) .or. Sy(i,j) > Infinity) then
!          print *, "Sy Error"
!      endif
!   enddo
! enddo
 
end subroutine COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE

!----------------------------------------------------------------------
!  Andy Heidinger's Extinction Routine  (km^-1)
! input:
!        Tc = cloud temperature
!        Ec = cloud emissivity
! output:
!        Ice_Cloud_Extinction in 1/km units
!----------------------------------------------------------------------
subroutine DETERMINE_ACHA_ICE_EXTINCTION(Tc,Ec,Beta,Ice_Cloud_Extinction)
   real, intent(in):: Tc, Ec, Beta
   real, intent(out):: Ice_Cloud_Extinction
   integer, save:: not_first_time = 0
   integer:: iec, itc, i
   real:: xTc, Log10_Reff, Qe_006um, Qe_110um

   Log10_Reff = 1.0
   if (Beta /= MISSING_VALUE_REAL4) then
      Log10_Reff = Re_Beta_110um_COEF_ICE(0) +  &
                   Re_Beta_110um_COEF_ICE(1)*(Beta-1.0) + &
                   Re_Beta_110um_COEF_ICE(2)*(Beta-1.0)**2 + &
                   Re_Beta_110um_COEF_ICE(3)*(Beta-1.0)**3
      Log10_Reff =1.0/Log10_Reff   !fit is in 1/log10_reff
      Log10_Reff = min(2.0,max(Log10_Reff, 0.6))   !constrain to 4 to 100 microns
      Qe_006um = Qe_006um_COEF_ICE(0) +  &
               Qe_006um_COEF_ICE(1)*Log10_Reff + &
               Qe_006um_COEF_ICE(2)*Log10_Reff**2

      Qe_110um = Qe_110um_COEF_ICE(0) +  &
                 Qe_110um_COEF_ICE(1)*Log10_Reff + &
                 Qe_110um_COEF_ICE(2)*Log10_Reff**2
   endif

   !---- if first time, set coefficients in memory
   if (not_first_time == 0) then 
     not_first_time = 1
     ice_ext_coef(1,:) =    [1.2316e-02, 2.4031e-03,-4.7896e-05, 4.5502e-07]
     ice_ext_coef(2,:) =    [3.3809e-02, 9.7992e-03,-2.7210e-04, 2.9489e-06]
     ice_ext_coef(3,:) =    [9.6130e-02, 8.8489e-03,-1.6911e-04, 1.8791e-06]
     ice_ext_coef(4,:) =    [1.0245e-01, 1.9004e-02,-4.9367e-04, 4.8926e-06]
     ice_ext_coef(5,:) =    [1.6542e-01, 2.4050e-02,-7.0344e-04, 7.7430e-06]
     ice_ext_coef(6,:) =    [2.5281e-01, 2.0355e-02,-4.5441e-04, 5.0307e-06]
     ice_ext_coef(7,:) =    [2.7558e-01, 4.3940e-02,-1.3112e-03, 1.3254e-05]
     ice_ext_coef(8,:) =    [5.1518e-01, 2.5070e-02,-5.8033e-04, 7.6746e-06]
     ice_ext_coef(9,:) =    [5.4931e-01, 7.4288e-02,-1.9975e-03, 1.8927e-05]
     ice_ext_coef(10,:) =   [2.6006e+00, 4.7639e-02, 7.5674e-04,-6.2741e-06]
   endif

   !--- determine which emissivity to use
   iec = int((Ec - Ec_Ext_Min) / Ec_Ext_Bin) + 1
   iec = max(1,iec)
   iec = min(Nec_Ext,iec)

   !--- compute ice cloud extintion (1/km)
   Ice_Cloud_Extinction = 0.0
   xTc = max(0.0,Tc - Tc_Ext_Offset)
   xTc = min(xTc,90.0)
   do i = 1,M_Ext
      Ice_Cloud_Extinction = Ice_Cloud_Extinction + &
                             Ice_Ext_Coef(iec,i)*(xTc**(i-1))
   enddo

   !--- convert from 532nm (0.65 um) to 11 um
   Ice_Cloud_Extinction = Ice_Cloud_Extinction * Qe_110um / Qe_006um

   !--- limit
   Ice_Cloud_Extinction = min(10.0,max(0.01,Ice_Cloud_Extinction))

end subroutine DETERMINE_ACHA_ICE_EXTINCTION
!----------------------------------------------------------------------
!  Yue Li's Extinction Routine
!----------------------------------------------------------------------
subroutine DETERMINE_ACHA_EXTINCTION(Cloud_Type,Tc,Symbol,Cloud_Extinction)

   integer(kind=int1), intent(in):: Cloud_Type
   real, intent(in):: Tc
   type(ACHA_SYMBOL_STRUCT), intent(in) :: Symbol
   real, intent(out):: Cloud_Extinction
   integer:: Itemp

      Cloud_Extinction = WATER_EXTINCTION
      Itemp = int(Tc)

       if (Cloud_Type == Symbol%OPAQUE_ICE_TYPE .or. &
           Cloud_Type == Symbol%OVERSHOOTING_TYPE) then
           select case (Itemp)
            case (:199)    ; Cloud_Extinction = ICE_EXTINCTION1
            case (200:219) ; Cloud_Extinction = ICE_EXTINCTION2
            case (220:239) ; Cloud_Extinction = ICE_EXTINCTION3
            case (240:259) ; Cloud_Extinction = ICE_EXTINCTION4
            case (260:)    ; Cloud_Extinction = ICE_EXTINCTION5
           end select
       endif

       if (Cloud_Type == Symbol%CIRRUS_TYPE .or. &
           Cloud_Type == Symbol%OVERLAP_TYPE) then
           select case (Itemp)
            case (:199)    ; Cloud_Extinction = CIRRUS_EXTINCTION1
            case (200:219) ; Cloud_Extinction = CIRRUS_EXTINCTION2
            case (220:239) ; Cloud_Extinction = CIRRUS_EXTINCTION3
            case (240:259) ; Cloud_Extinction = CIRRUS_EXTINCTION4
            case (260:)    ; Cloud_Extinction = CIRRUS_EXTINCTION5
          end select
       endif

end subroutine DETERMINE_ACHA_EXTINCTION

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module ACHA_RTM_MOD
