!$Id: acha_module.f90 3876 2020-06-18 13:34:40Z yli $
module ACHA_FULL_RETRIEVAL_MOD
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
  use CX_REAL_BOOLEAN_MOD

  use ACHA_SERVICES_MOD, only : &
           real4, int1, int4, real8, dtor, acha_rtm_nwp_struct, &
           acha_input_struct

  use ACHA_RTM_MOD

  use ACHA_NUM_MOD

  implicit none

  public:: FULL_ACHA_RETRIEVAL
  private:: COMPUTE_FORWARD_MODEL_AND_KERNEL

  !--- include the non-system specific variables
  include 'include/acha_parameters.inc'

  real, private, PARAMETER:: MISSING_VALUE_REAL4 = -999.0
  integer(kind=int1), private, PARAMETER:: MISSING_VALUE_integer1 = -128_int1
  integer(kind=int4), private, PARAMETER:: MISSING_VALUE_integer4 = -999

  contains 
 !----------------------------------------------------------------------------------
 ! The pixel-level ACHA Retrieval subroutine
 !----------------------------------------------------------------------------------
 subroutine FULL_ACHA_RETRIEVAL(Input,Acha_Mode_Flag, Symbol, &
                          Num_Obs,Num_Param,y,y_variance,f,x_Ap,Sa_Inv,x,Sx,AKM, &
                          Conv_Test, Cost, Goodness, Convergence_Criteria, Hght_Prof,&
                          Tsfc_Est,T_Tropo,Z_Tropo,P_Tropo, Cloud_Type, Cos_Zen, Zc_Base, &
                          ACHA_RTM_NWP, &
                          Beta_110um_142um_Coef_Water, &
                          Beta_110um_139um_Coef_Water, &
                          Beta_110um_136um_Coef_Water, &
                          Beta_110um_133um_Coef_Water, &
                          Beta_110um_104um_Coef_Water, &
                          Beta_110um_097um_Coef_Water, &
                          Beta_110um_085um_Coef_Water, &
                          Beta_110um_073um_Coef_Water, &
                          Beta_110um_067um_Coef_Water, &
                          Beta_110um_062um_Coef_Water, &
                          Beta_110um_038um_Coef_Water, &
                          Beta_110um_142um_Coef_Ice, &
                          Beta_110um_139um_Coef_Ice, &
                          Beta_110um_136um_Coef_Ice, &
                          Beta_110um_133um_Coef_Ice, &
                          Beta_110um_104um_Coef_Ice, &
                          Beta_110um_097um_Coef_Ice, &
                          Beta_110um_085um_Coef_Ice, &
                          Beta_110um_073um_Coef_Ice, &
                          Beta_110um_067um_Coef_Ice, &
                          Beta_110um_062um_Coef_Ice, &
                          Beta_110um_038um_Coef_Ice, &
                          Emiss_Sfc_038um, &
                          Emiss_Sfc_062um, &
                          Emiss_Sfc_067um, &
                          Emiss_Sfc_073um, &
                          Emiss_Sfc_085um, &
                          Emiss_Sfc_097um, &
                          Emiss_Sfc_104um, &
                          Emiss_Sfc_110um, &
                          Emiss_Sfc_120um, &
                          Emiss_Sfc_133um, &
                          Emiss_Sfc_136um, &
                          Emiss_Sfc_139um, &
                          Emiss_Sfc_142um, &
                          Tc_Eff, &
                          Converged_Flag, &
                          Fail_Flag, &
                          Dump_Diag, &
                          Lun_Iter_Dump)
                   
 character(len=*), intent(in):: Acha_Mode_Flag
 type(ACHA_INPUT_STRUCT), intent(in) :: Input
 type(ACHA_SYMBOL_STRUCT), intent(in) :: Symbol
 type(ACHA_RTM_NWP_STRUCT), intent(in) :: ACHA_RTM_NWP
 integer, intent(in):: Num_Obs,Num_Param
 integer(kind=int1), intent(in):: Cloud_Type
 real, intent(in), dimension(:):: x_Ap, y, y_variance
 real, intent(in), dimension(:,:)::  Sa_Inv
 real, intent(in):: Convergence_Criteria
 real, intent(in), dimension(:):: Hght_Prof
 real, intent(in):: Cos_Zen, Tsfc_Est,T_Tropo,Z_Tropo, P_Tropo
 
 real, dimension(0:), intent(in):: Beta_110um_142um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_139um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_136um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_133um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_104um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_097um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_085um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_073um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_067um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_062um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_038um_Coef_Water
 real, dimension(0:), intent(in):: Beta_110um_142um_Coef_Ice
 real, dimension(0:), intent(in):: Beta_110um_139um_Coef_Ice
 real, dimension(0:), intent(in):: Beta_110um_136um_Coef_Ice
 real, dimension(0:), intent(in):: Beta_110um_133um_Coef_Ice
 real, dimension(0:), intent(in):: Beta_110um_104um_Coef_Ice
 real, dimension(0:), intent(in):: Beta_110um_097um_Coef_Ice
 real, dimension(0:), intent(in):: Beta_110um_085um_Coef_Ice
 real, dimension(0:), intent(in):: Beta_110um_073um_Coef_Ice
 real, dimension(0:), intent(in):: Beta_110um_067um_Coef_Ice
 real, dimension(0:), intent(in):: Beta_110um_062um_Coef_Ice
 real, dimension(0:), intent(in):: Beta_110um_038um_Coef_Ice
 real(kind=real4), intent(in):: Emiss_Sfc_038um
 real(kind=real4), intent(in):: Emiss_Sfc_062um
 real(kind=real4), intent(in):: Emiss_Sfc_067um
 real(kind=real4), intent(in):: Emiss_Sfc_073um
 real(kind=real4), intent(in):: Emiss_Sfc_085um
 real(kind=real4), intent(in):: Emiss_Sfc_097um
 real(kind=real4), intent(in):: Emiss_Sfc_104um
 real(kind=real4), intent(in):: Emiss_Sfc_110um
 real(kind=real4), intent(in):: Emiss_Sfc_120um
 real(kind=real4), intent(in):: Emiss_Sfc_133um
 real(kind=real4), intent(in):: Emiss_Sfc_136um
 real(kind=real4), intent(in):: Emiss_Sfc_139um
 real(kind=real4), intent(in):: Emiss_Sfc_142um
 integer, intent(in):: Lun_Iter_Dump
 logical, intent(in):: Dump_Diag
 real, dimension(:), intent(out):: x
 real, dimension(:), intent(out):: f
 real, dimension(:,:), intent(out):: Sx
 real, dimension(:,:), intent(out):: AKM
 real, intent(out):: Conv_Test
 real, intent(out):: Cost
 real, intent(out):: Goodness
 integer, intent(out):: Converged_Flag, Fail_Flag
 real, intent(out):: Tc_Eff
 real, intent(out):: Zc_Base

 integer:: Lev_Idx
 real, dimension(Num_Obs,Num_Param):: K
!real, dimension(Num_Obs):: f
 real, dimension(Num_Param):: Delta_x, Delta_x_prev
 real, dimension(Num_Obs,Num_Obs):: Sy
 real, dimension(Num_Obs):: Emiss_Vector
!real, dimension(Num_Param,Num_Param):: AKM

 real:: Tc_Temp, Pc_Temp, Zc_Temp, Ec_Temp, Beta_Temp
 real:: Ts_Temp, Ps_Temp, Zs_Temp
 integer (kind=int4):: NWP_Profile_Inversion_Flag

 !--- ch20 variables
 real:: Rad_Ac_038um
 real:: Trans_Ac_038um
 real:: Trans_Bc_038um
 real:: Rad_Clear_038um

 !--- ch37 variables
 real:: Rad_Ac_062um
 real:: Trans_Ac_062um
 real:: Trans_Bc_062um
 real:: Rad_Clear_062um

 !--- ch27 variables
 real:: Rad_Ac_067um
 real:: Trans_Ac_067um
 real:: Trans_Bc_067um
 real:: Rad_Clear_067um

 !--- ch28 variables
 real:: Rad_Ac_073um
 real:: Trans_Ac_073um
 real:: Trans_Bc_073um
 real:: Rad_Clear_073um

 !--- ch29 variables
 real:: Rad_Ac_085um
 real:: Trans_Ac_085um
 real:: Trans_Bc_085um
 real:: Rad_Clear_085um

 !--- ch30 variables
 real:: Rad_Ac_097um
 real:: Trans_Ac_097um
 real:: Trans_Bc_097um
 real:: Rad_Clear_097um

 !--- ch38 variables
 real:: Rad_Ac_104um
 real:: Trans_Ac_104um
 real:: Trans_Bc_104um
 real:: Rad_Clear_104um

 !--- ch31 variables
 real:: Rad_Ac_110um
 real:: Trans_Ac_110um
 real:: Trans_Bc_110um
 real:: Rad_Clear_110um

 !--- ch32 variables
 real:: Rad_Ac_120um
 real:: Trans_Ac_120um
 real:: Trans_Bc_120um
 real:: Rad_Clear_120um

 !--- ch33 variables
 real:: Rad_Ac_133um
 real:: Trans_Ac_133um
 real:: Trans_Bc_133um
 real:: Rad_Clear_133um

 !--- ch34 variables
 real:: Rad_Ac_136um
 real:: Trans_Ac_136um
 real:: Trans_Bc_136um
 real:: Rad_Clear_136um

 !--- ch35 variables
 real:: Rad_Ac_139um
 real:: Trans_Ac_139um
 real:: Trans_Bc_139um
 real:: Rad_Clear_139um

 !--- ch36 variables
 real:: Rad_Ac_142um
 real:: Trans_Ac_142um
 real:: Trans_Bc_142um
 real:: Rad_Clear_142um

 real:: R4_Dummy
 real:: Tc_Base
 real:: Zc_Thick
 real:: Cloud_Opd
 real:: Cloud_Extinction
 real:: Extinction_Ratio_11_065

 integer:: Iter_Idx, ierror
 integer:: nx

!----------------------------------------------------------
Iter_Idx = 0
Converged_Flag = Symbol%NO
Fail_Flag =  Symbol%NO
Delta_x_prev = MISSING_VALUE_REAL4
nx = size(x)

!---- assign x to the first guess
x = x_Ap

Retrieval_Loop: do

  Iter_Idx = Iter_Idx + 1

  if (Dump_Diag) write(unit=Lun_Iter_Dump,fmt=*) "==> Iter_Idx = ", Iter_Idx

  !---------------------------------------------------------------------
  ! estimate clear-sky radiative transfer terms used in forward model
  !---------------------------------------------------------------------
  Tc_Temp = x(1)
  Ec_Temp = x(2)
  Beta_Temp = x(3)
  if (nx > 3) then
     Ts_Temp = x(4)
  else
     Ts_Temp = Tsfc_Est
  endif

  call KNOWING_T_COMPUTE_P_Z_BOTTOM_UP(ACHA_RTM_NWP,Symbol,Cloud_Type,Pc_temp,Tc_temp,Zc_Temp, &
                             T_Tropo,Z_Tropo,P_Tropo,Lev_Idx,ierror,NWP_Profile_Inversion_Flag)

  call KNOWING_T_COMPUTE_P_Z_BOTTOM_UP(ACHA_RTM_NWP,Symbol,Cloud_Type,Ps_temp,Ts_temp,Zs_Temp, &
                             T_Tropo,Z_Tropo,P_Tropo,Lev_Idx,ierror,NWP_Profile_Inversion_Flag)

  !--- If GOES-17 11 um is bad, at this point, all 11um variables are filled
  !--- with 10.4 um data.
  call COMPUTE_CLEAR_SKY_TERMS(Acha_Mode_Flag, Zc_Temp, Zs_Temp, Ts_Temp, Hght_Prof,  &
                               Input%Chan_Idx_038um, Input%Chan_Idx_062um, Input%Chan_Idx_067um, Input%Chan_Idx_073um, &
                               Input%Chan_Idx_085um, Input%Chan_Idx_097um, Input%Chan_Idx_104um, Input%Chan_Idx_110um, Input%Chan_Idx_120um, &
                               Input%Chan_Idx_133um, Input%Chan_Idx_136um, Input%Chan_Idx_139um, Input%Chan_Idx_142um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_038um, ACHA_RTM_NWP%Atm_Trans_Prof_038um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_062um, ACHA_RTM_NWP%Atm_Trans_Prof_062um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_067um, ACHA_RTM_NWP%Atm_Trans_Prof_067um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_073um, ACHA_RTM_NWP%Atm_Trans_Prof_073um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_085um, ACHA_RTM_NWP%Atm_Trans_Prof_085um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_097um, ACHA_RTM_NWP%Atm_Trans_Prof_097um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_104um, ACHA_RTM_NWP%Atm_Trans_Prof_104um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_110um, ACHA_RTM_NWP%Atm_Trans_Prof_110um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_120um, ACHA_RTM_NWP%Atm_Trans_Prof_120um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_133um, ACHA_RTM_NWP%Atm_Trans_Prof_133um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_136um, ACHA_RTM_NWP%Atm_Trans_Prof_136um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_139um, ACHA_RTM_NWP%Atm_Trans_Prof_139um, &
                               ACHA_RTM_NWP%Atm_Rad_Prof_142um, ACHA_RTM_NWP%Atm_Trans_Prof_142um, &
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

  !--------------------------------------------------
  ! Determine Slope of Planck Emission through Cloud
  !--------------------------------------------------

  !--- default - no accounting for vertical extent
  Zc_Thick = 0.0
  Zc_Base = Zc_Temp
  Tc_Base = Tc_Temp
  Zc_Thick = 0.0

  if (USE_LINEAR_IN_OPD_EMISSION) then 
    Cloud_Opd = -1.0*alog(1.0-Ec_Temp)
    Cloud_Opd = max(0.01,min(10.0,Cloud_Opd))
    Cloud_Opd = Cloud_Opd * Cos_Zen

    !--- Yue Li's routine does ice as f(t) and water
    call DETERMINE_ACHA_EXTINCTION(Cloud_Type,Tc_Temp,Symbol,Cloud_Extinction)

    !--- Andy Heidinger's routine replaces ice values but not water
    if (Cloud_Type == symbol%CIRRUS_TYPE .or. &
        Cloud_Type == symbol%OVERLAP_TYPE .or. &
        Cloud_Type == symbol%OPAQUE_ICE_TYPE) then
        call DETERMINE_ACHA_ICE_EXTINCTION(Tc_Temp,Ec_Temp,Beta_Temp,Cloud_Extinction)
        Cloud_Extinction = ICE_EXTINCTION_TUNING_FACTOR * Cloud_Extinction
    endif

    Zc_Thick = 1000.0*Cloud_Opd / Cloud_Extinction
    Zc_Base = Zc_Temp - Zc_Thick
    Zc_Base = max(Zc_Base,Zs_Temp) !constrain to be above surface
    Zc_Base = max(Zc_Base,100.0)   !constrain to be positive (greater than 100 m)
    call KNOWING_Z_COMPUTE_T_P(ACHA_RTM_NWP,R4_Dummy,Tc_Base,Zc_Base,Lev_Idx)
  endif

   if (Dump_Diag) then
     write(unit=Lun_Iter_Dump,fmt=*) "Zc_temp = ", Zc_Temp
     write(unit=Lun_Iter_Dump,fmt=*) "Zc_Thick = ", Zc_Thick
     write(unit=Lun_Iter_Dump,fmt=*) "Zc_Base = ", Zc_Base 
     write(unit=Lun_Iter_Dump,fmt=*) "Tc_Base = ", Tc_Base 
   endif

  !--------------------------------------------------
  ! call forward models
  !--------------------------------------------------

  !--- At this point, if GOES-17 mitigation is using 10.4 um data
  !--- the following are switched from 11 um to 104 um data:
  !--- Input%Chan_Idx_104um,

  call COMPUTE_FORWARD_MODEL_AND_KERNEL(Input,Acha_Mode_Flag,  &
           x, &
           Rad_Clear_038um, Rad_Ac_038um, Trans_Ac_038um, Trans_Bc_038um, &
           Rad_Clear_062um, Rad_Ac_062um, Trans_Ac_062um, Trans_Bc_062um, &
           Rad_Clear_067um, Rad_Ac_067um, Trans_Ac_067um, Trans_Bc_067um, &
           Rad_Clear_073um, Rad_Ac_073um, Trans_Ac_073um, Trans_Bc_073um, &
           Rad_Clear_085um, Rad_Ac_085um, Trans_Ac_085um, Trans_Bc_085um, &
           Rad_Clear_097um, Rad_Ac_097um, Trans_Ac_097um, Trans_Bc_097um, &
           Rad_Clear_104um, Rad_Ac_104um, Trans_Ac_104um, Trans_Bc_104um, &
           Rad_Clear_110um, Rad_Ac_110um, Trans_Ac_110um, Trans_Bc_110um, &
           Rad_Clear_120um, Rad_Ac_120um, Trans_Ac_120um, Trans_Bc_120um, &
           Rad_Clear_133um, Rad_Ac_133um, Trans_Ac_133um, Trans_Bc_133um, &
           Rad_Clear_136um, Rad_Ac_136um, Trans_Ac_136um, Trans_Bc_136um, &
           Rad_Clear_139um, Rad_Ac_139um, Trans_Ac_139um, Trans_Bc_139um, &
           Rad_Clear_142um, Rad_Ac_142um, Trans_Ac_142um, Trans_Bc_142um, &
           Beta_110um_142um_Coef_Water, &
           Beta_110um_139um_Coef_Water, &
           Beta_110um_136um_Coef_Water, &
           Beta_110um_133um_Coef_Water, &
           Beta_110um_104um_Coef_Water, &
           Beta_110um_097um_Coef_Water, &
           Beta_110um_085um_Coef_Water, &
           Beta_110um_073um_Coef_Water, &
           Beta_110um_067um_Coef_Water, &
           Beta_110um_062um_Coef_Water, &
           Beta_110um_038um_Coef_Water, &
           Beta_110um_142um_Coef_Ice, &
           Beta_110um_139um_Coef_Ice, &
           Beta_110um_136um_Coef_Ice, &
           Beta_110um_133um_Coef_Ice, &
           Beta_110um_104um_Coef_Ice, &
           Beta_110um_097um_Coef_Ice, &
           Beta_110um_085um_Coef_Ice, &
           Beta_110um_073um_Coef_Ice, &
           Beta_110um_067um_Coef_Ice, &
           Beta_110um_062um_Coef_Ice, &
           Beta_110um_038um_Coef_Ice, &
           Emiss_Sfc_038um, &
           Emiss_Sfc_062um, &
           Emiss_Sfc_067um, &
           Emiss_Sfc_073um, &
           Emiss_Sfc_085um, &
           Emiss_Sfc_097um, &
           Emiss_Sfc_104um, &
           Emiss_Sfc_110um, &
           Emiss_Sfc_120um, &
           Emiss_Sfc_133um, &
           Emiss_Sfc_136um, &
           Emiss_Sfc_139um, &
           Emiss_Sfc_142um, &
           f,K,Emiss_Vector,Tc_Base,Tsfc_Est,Tc_Eff)

   if (Dump_Diag) then 
     write(unit=Lun_Iter_Dump,fmt=*) "x = ", x
     write(unit=Lun_Iter_Dump,fmt=*) "f = ", f
   endif
  !--------------------------------------------------
  ! compute the Sy convariance matrix
  !--------------------------------------------------
  call COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE( &
                                                 Emiss_Vector, &
                                                 Acha_Mode_Flag, &
                                                 y_variance, &
                                                 Sy) 

  !--------------------------------------------------
  ! call OE routine to advance the Iteration
  !--------------------------------------------------
  call OPTIMAL_ESTIMATION(Iter_Idx,Iter_Idx_Max,Num_Param,Num_Obs, &
                         Convergence_Criteria,Delta_X_Max, &
                         y,f,x,x_Ap,K,Sy,Sa_inv, &
                         Sx,AKM,Delta_x,Delta_x_prev, &
                         Conv_Test,Cost, Goodness, &
                         Converged_Flag,Fail_Flag)
  Delta_x_Prev = Delta_x

  if (Dump_Diag) then 
    write(unit=Lun_Iter_Dump,fmt=*) "Delta_X = ", Delta_X
    write(unit=Lun_Iter_Dump,fmt=*) "Conv_Test = ", Conv_Test
    write(unit=Lun_Iter_Dump,fmt=*) "Conv_Crit = ", Convergence_Criteria
    write(unit=Lun_Iter_Dump,fmt=*) "Cost = ", Cost
    write(unit=Lun_Iter_Dump,fmt=*) "Goodness = ", Goodness
    write(unit=Lun_Iter_Dump,fmt=*) "Converged_Flag = ", Converged_Flag
    write(unit=Lun_Iter_Dump,fmt=*) "AKM = ", AKM(1,1),AKM(2,2), &
                                     AKM(3,3), AKM(4,4), AKM(5,5)
  endif
  
  !--- check for a failed Iteration
  if (Fail_Flag == Symbol%YES) then
!    print *, "Failed "
     exit
  endif

  !--------------------------------------------------------
  ! exit retrieval loop if converged
  !--------------------------------------------------------
  if (Converged_Flag == Symbol%YES) then
      !print *, "success"
       exit
  endif

  !---------------------------------------------------------
  ! update retrieved Output%Vector
  !---------------------------------------------------------
  x = x + Delta_X
  if (Dump_Diag) write(unit=Lun_Iter_Dump,fmt=*) "new x = ", x

  !-------------------------------------------------------
  ! constrain to reasonable values
  !-------------------------------------------------------
  x(1) = max(MIN_ALLOWABLE_TC,min(Tsfc_Est+5,x(1)))   
  x(2) = max(0.0,min(x(2),1.0))
  x(3) = max(0.8,min(x(3),1.8))
  if (nx > 3) x(4) = max(x(1),min(Tsfc_Est,x(4)))    
  if (nx > 3) x(5) = max(0.0,min(1.0,x(5)))    
  if (Dump_Diag) write(unit=Lun_Iter_Dump,fmt=*) "constrained x = ", x

end do Retrieval_Loop

 if (Dump_Diag) then
     write(unit=Lun_Iter_Dump,fmt=*) "========================================================"
     write(unit=Lun_Iter_Dump,fmt=*) "Returned from ACHA_FULL_RETRIEVAL" 
     write(unit=Lun_Iter_Dump,fmt=*) "========================================================"
     write(unit=Lun_Iter_Dump,fmt=*) "final f = ", f
     write(unit=Lun_Iter_Dump,fmt=*) "Converged_Flag = ", Converged_Flag
     write(unit=Lun_Iter_Dump,fmt=*) "Fail_Flag = ", Fail_Flag
 endif

end subroutine FULL_ACHA_RETRIEVAL

!---------------------------------------------------------------------
!--- Compute the Forward Model Estimate (f) and its Kernel (df/dx)
!---------------------------------------------------------------------
subroutine COMPUTE_FORWARD_MODEL_AND_KERNEL( &
           Input,Acha_Mode_Flag,  &
           x,    &
           Rad_Clear_038um, Rad_Ac_038um, Trans_Ac_038um, Trans_Bc_038um,    &
           Rad_Clear_062um, Rad_Ac_062um, Trans_Ac_062um, Trans_Bc_062um,    &
           Rad_Clear_067um, Rad_Ac_067um, Trans_Ac_067um, Trans_Bc_067um,    &
           Rad_Clear_073um, Rad_Ac_073um, Trans_Ac_073um, Trans_Bc_073um,    &
           Rad_Clear_085um, Rad_Ac_085um, Trans_Ac_085um, Trans_Bc_085um,    &
           Rad_Clear_097um, Rad_Ac_097um, Trans_Ac_097um, Trans_Bc_097um,    &
           Rad_Clear_104um, Rad_Ac_104um, Trans_Ac_104um, Trans_Bc_104um,    &
           Rad_Clear_110um, Rad_Ac_110um, Trans_Ac_110um, Trans_Bc_110um,    &
           Rad_Clear_120um, Rad_Ac_120um, Trans_Ac_120um, Trans_Bc_120um,    &
           Rad_Clear_133um, Rad_Ac_133um, Trans_Ac_133um, Trans_Bc_133um,&
           Rad_Clear_136um, Rad_Ac_136um, Trans_Ac_136um, Trans_Bc_136um,&
           Rad_Clear_139um, Rad_Ac_139um, Trans_Ac_139um, Trans_Bc_139um,&
           Rad_Clear_142um, Rad_Ac_142um, Trans_Ac_142um, Trans_Bc_142um,&
           Beta_110um_142um_Coef_Water, &
           Beta_110um_139um_Coef_Water, &
           Beta_110um_136um_Coef_Water, &
           Beta_110um_133um_Coef_Water, &
           Beta_110um_104um_Coef_Water,   &
           Beta_110um_097um_Coef_Water,   &
           Beta_110um_085um_Coef_Water,   &
           Beta_110um_073um_Coef_Water,   &
           Beta_110um_067um_Coef_Water,   &
           Beta_110um_062um_Coef_Water,   &
           Beta_110um_038um_Coef_Water,   &
           Beta_110um_142um_Coef_Ice, &
           Beta_110um_139um_Coef_Ice, &
           Beta_110um_136um_Coef_Ice, &
           Beta_110um_133um_Coef_Ice, &
           Beta_110um_104um_Coef_Ice,   &
           Beta_110um_097um_Coef_Ice,   &
           Beta_110um_085um_Coef_Ice,   &
           Beta_110um_073um_Coef_Ice,   &
           Beta_110um_067um_Coef_Ice,   &
           Beta_110um_062um_Coef_Ice,   &
           Beta_110um_038um_Coef_Ice,   &
           Emiss_Sfc_038um, &
           Emiss_Sfc_062um, &
           Emiss_Sfc_067um, &
           Emiss_Sfc_073um, &
           Emiss_Sfc_085um, &
           Emiss_Sfc_097um, &
           Emiss_Sfc_104um, &
           Emiss_Sfc_110um, &
           Emiss_Sfc_120um, &
           Emiss_Sfc_133um, &
           Emiss_Sfc_136um, &
           Emiss_Sfc_139um, &
           Emiss_Sfc_142um, &
           f, & 
           K, &
           Emiss_Vector, &
           Tc_Base,Tsfc,Tc_Eff)

  type(acha_input_struct), intent(in):: Input
  character(len=*), intent(in):: Acha_Mode_Flag
  real(kind=real4), dimension(:), intent(in):: x
  real(kind=real4), intent(in):: Rad_Clear_038um
  real(kind=real4), intent(in):: Rad_Ac_038um
  real(kind=real4), intent(in):: Trans_Ac_038um
  real(kind=real4), intent(in):: Trans_Bc_038um
  real(kind=real4), intent(in):: Rad_Clear_062um
  real(kind=real4), intent(in):: Rad_Ac_062um
  real(kind=real4), intent(in):: Trans_Ac_062um
  real(kind=real4), intent(in):: Trans_Bc_062um
  real(kind=real4), intent(in):: Rad_Clear_067um
  real(kind=real4), intent(in):: Rad_Ac_067um
  real(kind=real4), intent(in):: Trans_Ac_067um
  real(kind=real4), intent(in):: Trans_Bc_067um
  real(kind=real4), intent(in):: Rad_Clear_073um
  real(kind=real4), intent(in):: Rad_Ac_073um
  real(kind=real4), intent(in):: Trans_Ac_073um
  real(kind=real4), intent(in):: Trans_Bc_073um
  real(kind=real4), intent(in):: Rad_Clear_085um
  real(kind=real4), intent(in):: Rad_Ac_085um
  real(kind=real4), intent(in):: Trans_Ac_085um
  real(kind=real4), intent(in):: Trans_Bc_085um
  real(kind=real4), intent(in):: Rad_Clear_097um
  real(kind=real4), intent(in):: Rad_Ac_097um
  real(kind=real4), intent(in):: Trans_Ac_097um
  real(kind=real4), intent(in):: Trans_Bc_097um
  real(kind=real4), intent(in):: Rad_Clear_104um
  real(kind=real4), intent(in):: Rad_Ac_104um
  real(kind=real4), intent(in):: Trans_Ac_104um
  real(kind=real4), intent(in):: Trans_Bc_104um
  real(kind=real4), intent(in):: Rad_Clear_110um
  real(kind=real4), intent(in):: Rad_Ac_110um
  real(kind=real4), intent(in):: Trans_Ac_110um
  real(kind=real4), intent(in):: Trans_Bc_110um
  real(kind=real4), intent(in):: Rad_Clear_120um
  real(kind=real4), intent(in):: Rad_Ac_120um
  real(kind=real4), intent(in):: Trans_Ac_120um
  real(kind=real4), intent(in):: Trans_Bc_120um
  real(kind=real4), intent(in):: Rad_Clear_133um
  real(kind=real4), intent(in):: Rad_Ac_133um
  real(kind=real4), intent(in):: Trans_Ac_133um
  real(kind=real4), intent(in):: Trans_Bc_133um
  real(kind=real4), intent(in):: Rad_Clear_136um
  real(kind=real4), intent(in):: Rad_Ac_136um
  real(kind=real4), intent(in):: Trans_Ac_136um
  real(kind=real4), intent(in):: Trans_Bc_136um
  real(kind=real4), intent(in):: Rad_Clear_139um
  real(kind=real4), intent(in):: Rad_Ac_139um
  real(kind=real4), intent(in):: Trans_Ac_139um
  real(kind=real4), intent(in):: Trans_Bc_139um
  real(kind=real4), intent(in):: Rad_Clear_142um
  real(kind=real4), intent(in):: Rad_Ac_142um
  real(kind=real4), intent(in):: Trans_Ac_142um
  real(kind=real4), intent(in):: Trans_Bc_142um
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_142um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_139um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_136um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_133um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_104um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_097um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_085um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_073um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_067um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_062um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_038um_Coef_Water
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_142um_Coef_Ice
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_139um_Coef_Ice
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_136um_Coef_Ice
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_133um_Coef_Ice
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_104um_Coef_Ice
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_097um_Coef_Ice
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_085um_Coef_Ice
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_073um_Coef_Ice
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_067um_Coef_Ice
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_062um_Coef_Ice
  real(kind=real4), dimension(0:), intent(in):: Beta_110um_038um_Coef_Ice
  real(kind=real4), intent(in):: Emiss_Sfc_038um
  real(kind=real4), intent(in):: Emiss_Sfc_062um
  real(kind=real4), intent(in):: Emiss_Sfc_067um
  real(kind=real4), intent(in):: Emiss_Sfc_073um
  real(kind=real4), intent(in):: Emiss_Sfc_085um
  real(kind=real4), intent(in):: Emiss_Sfc_097um
  real(kind=real4), intent(in):: Emiss_Sfc_104um
  real(kind=real4), intent(in):: Emiss_Sfc_110um
  real(kind=real4), intent(in):: Emiss_Sfc_120um
  real(kind=real4), intent(in):: Emiss_Sfc_133um
  real(kind=real4), intent(in):: Emiss_Sfc_136um
  real(kind=real4), intent(in):: Emiss_Sfc_139um
  real(kind=real4), intent(in):: Emiss_Sfc_142um
  real(kind=real4), intent(in):: Tc_Base
  real(kind=real4), intent(in):: Tsfc
  real(kind=real4), dimension(:), intent(out):: f 
  real(kind=real4), dimension(:,:), intent(out):: K
  real(kind=real4), dimension(:), intent(out):: Emiss_Vector
  real(kind=real4), intent(out):: Tc_Eff

  real(kind=real4):: Tc
  real(kind=real4):: Ts
  real(kind=real4):: alpha
  real(kind=real4):: Emiss_038um
  real(kind=real4):: Emiss_062um
  real(kind=real4):: Emiss_067um
  real(kind=real4):: Emiss_073um
  real(kind=real4):: Emiss_085um
  real(kind=real4):: Emiss_097um
  real(kind=real4):: Emiss_104um
  real(kind=real4):: Emiss_110um
  real(kind=real4):: Emiss_120um
  real(kind=real4):: Emiss_133um
  real(kind=real4):: Emiss_136um
  real(kind=real4):: Emiss_139um
  real(kind=real4):: Emiss_142um
  real(kind=real4):: Beta_110um_120um

  !--- Kernel and forward model terms
  real:: f_T_110, dT_110_dTc, dT_110_dec, dT_110_dbeta, dT_110_dTs, dT_110_dalpha
  real:: f_Btd_110_038,dBtd_110_038_dTc, dBtd_110_038_dec, dBtd_110_038_dbeta, dBtd_110_038_dTs, dBtd_110_038_dalpha
  real:: f_Btd_110_062,dBtd_110_062_dTc, dBtd_110_062_dec, dBtd_110_062_dbeta, dBtd_110_062_dTs, dBtd_110_062_dalpha
  real:: f_Btd_110_067,dBtd_110_067_dTc, dBtd_110_067_dec, dBtd_110_067_dbeta, dBtd_110_067_dTs, dBtd_110_067_dalpha
  real:: f_Btd_110_073,dBtd_110_073_dTc, dBtd_110_073_dec, dBtd_110_073_dbeta, dBtd_110_073_dTs, dBtd_110_073_dalpha
  real:: f_Btd_110_085,dBtd_110_085_dTc, dBtd_110_085_dec, dBtd_110_085_dbeta, dBtd_110_085_dTs, dBtd_110_085_dalpha
  real:: f_Btd_110_097,dBtd_110_097_dTc, dBtd_110_097_dec, dBtd_110_097_dbeta, dBtd_110_097_dTs, dBtd_110_097_dalpha
  real:: f_Btd_110_104,dBtd_110_104_dTc, dBtd_110_104_dec, dBtd_110_104_dbeta, dBtd_110_104_dTs, dBtd_110_104_dalpha
  real:: f_Btd_110_120,dBtd_110_120_dTc, dBtd_110_120_dec, dBtd_110_120_dbeta, dBtd_110_120_dTs, dBtd_110_120_dalpha
  real:: f_Btd_110_133,dBtd_110_133_dTc, dBtd_110_133_dec, dBtd_110_133_dbeta, dBtd_110_133_dTs, dBtd_110_133_dalpha
  real:: f_Btd_110_136,dBtd_110_136_dTc, dBtd_110_136_dec, dBtd_110_136_dbeta, dBtd_110_136_dTs, dBtd_110_136_dalpha
  real:: f_Btd_110_139,dBtd_110_139_dTc, dBtd_110_139_dec, dBtd_110_139_dbeta, dBtd_110_139_dTs, dBtd_110_139_dalpha
  real:: f_Btd_110_142,dBtd_110_142_dTc, dBtd_110_142_dec, dBtd_110_142_dbeta, dBtd_110_142_dTs, dBtd_110_142_dalpha

  integer:: nx

  !---  for notational convenience, rename elements of x to local variables
  nx = size(x)
  Tc = x(1)
  Emiss_110um = min(x(2),0.999999)    !values must be below unity
  Beta_110um_120um = x(3)
  if (nx > 3) then 
     Ts = x(4)
     alpha = x(5)
  else
     Ts = Tsfc
     alpha = 1.0
  endif

 !----------------------------------------------------------------------------------------------
 ! Make Terms for the Kernel Matrix
 !----------------------------------------------------------------------------------------------

 !--- 11 um
 if (index(Acha_Mode_Flag,'110') > 0) then
   call BT_FM(Input%Chan_Idx_110um,Tc,Emiss_110um,Ts,Tc_Base,Emiss_Sfc_110um, &
              Rad_Ac_110um, Trans_Ac_110um, Trans_Bc_110um, Rad_Clear_110um, &
              f_T_110,dT_110_dTc,dT_110_dec,dT_110_dbeta,dT_110_dTs,dT_110_dalpha,Tc_Eff)
!print *, 'BT_FM = ', 
 endif

 !--- 11 - 38 um
 if (index(Acha_Mode_Flag,'038') > 0) then
   call BTD_FM(Input%Chan_Idx_038um, &
               Beta_110um_038um_Coef_Water, &
               Beta_110um_038um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_038um, alpha,  &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_038um, Trans_Ac_038um, Trans_Bc_038um, Rad_Clear_038um, &
               f_Btd_110_038,dBtd_110_038_dTc,dBtd_110_038_dec,dBtd_110_038_dbeta, &
               dBtd_110_038_dTs,dBtd_110_038_dalpha,Emiss_038um)
 endif

 !--- 11 - 6.2
 if (index(Acha_Mode_Flag,'062') > 0) then
   call BTD_FM(Input%Chan_Idx_062um, &
               Beta_110um_062um_Coef_Water, &
               Beta_110um_062um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_062um, alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_062um, Trans_Ac_062um, Trans_Bc_062um, Rad_Clear_062um, &
               f_Btd_110_062,dBtd_110_062_dTc,dBtd_110_062_dec,dBtd_110_062_dbeta, &
               dBtd_110_062_dTs,dBtd_110_062_dalpha,Emiss_062um)
 endif

 !--- 11 - 6.7
 if (index(Acha_Mode_Flag,'067') > 0) then
   call BTD_FM(Input%Chan_Idx_067um, &
               Beta_110um_067um_Coef_Water, &
               Beta_110um_067um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_067um, alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_067um, Trans_Ac_067um, Trans_Bc_067um, Rad_Clear_067um, &
               f_Btd_110_067,dBtd_110_067_dTc,dBtd_110_067_dec,dBtd_110_067_dbeta, &
               dBtd_110_067_dTs,dBtd_110_067_dalpha,Emiss_067um)
 endif

 !--- 11 - 7.3
 if (index(Acha_Mode_Flag,'073') > 0) then
   call BTD_FM(Input%Chan_Idx_073um, &
               Beta_110um_073um_Coef_Water, &
               Beta_110um_073um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_073um, alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_073um, Trans_Ac_073um, Trans_Bc_073um, Rad_Clear_073um, &
               f_Btd_110_073,dBtd_110_073_dTc,dBtd_110_073_dec,dBtd_110_073_dbeta, &
               dBtd_110_073_dTs,dBtd_110_073_dalpha,Emiss_073um)
 endif

 !--- 11 - 8.5 um
 if (index(Acha_Mode_Flag,'085') > 0) then
   call BTD_FM(Input%Chan_Idx_085um,  &
               Beta_110um_085um_Coef_Water, &
               Beta_110um_085um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_085um, alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_085um, Trans_Ac_085um, Trans_Bc_085um, Rad_Clear_085um, &
               f_Btd_110_085,dBtd_110_085_dTc,dBtd_110_085_dec,dBtd_110_085_dbeta, &
               dBtd_110_085_dTs,dBtd_110_085_dalpha,Emiss_085um)

 endif

 !--- 11 - 9.7 um
 if (index(Acha_Mode_Flag,'097') > 0) then
   call BTD_FM(Input%Chan_Idx_097um,  &
               Beta_110um_097um_Coef_Water, &
               Beta_110um_097um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_097um,alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_097um, Trans_Ac_097um, Trans_Bc_097um, Rad_Clear_097um, &
               f_Btd_110_085,dBtd_110_085_dTc,dBtd_110_085_dec,dBtd_110_085_dbeta, &
               dBtd_110_085_dTs,dBtd_110_085_dalpha,Emiss_097um)
 endif

 !--- 11 - 10.4 um
 if (index(Acha_Mode_Flag,'104') > 0) then
   call BTD_FM(Input%Chan_Idx_104um,  &
               Beta_110um_104um_Coef_Water, &
               Beta_110um_104um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_104um, alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_104um, Trans_Ac_104um, Trans_Bc_104um, Rad_Clear_104um, &
               f_Btd_110_104,dBtd_110_104_dTc,dBtd_110_104_dec,dBtd_110_104_dbeta, &
               dBtd_110_104_dTs,dBtd_110_104_dalpha,Emiss_104um)
 endif

 !--- 11 - 12 um
 if (index(Acha_Mode_Flag,'120') > 0) then
   call BTD_FM(Input%Chan_Idx_120um, &
               [1.0, 1.0, 0.0, 0.0], &
               [1.0, 1.0, 0.0, 0.0], &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_120um, alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_120um, Trans_Ac_120um, Trans_Bc_120um, Rad_Clear_120um, &
               f_Btd_110_120,dBtd_110_120_dTc,dBtd_110_120_dec,dBtd_110_120_dbeta, &
               dBtd_110_120_dTs,dBtd_110_120_dalpha,Emiss_120um)
 endif

 !--- 11 - 133 um
 if (index(Acha_Mode_Flag,'133') > 0) then
   call BTD_FM(Input%Chan_Idx_133um,  &
               Beta_110um_133um_Coef_Water, &
               Beta_110um_133um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_133um, alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_133um, Trans_Ac_133um, Trans_Bc_133um, Rad_Clear_133um, &
               f_Btd_110_133,dBtd_110_133_dTc,dBtd_110_133_dec,dBtd_110_133_dbeta, &
               dBtd_110_133_dTs,dBtd_110_133_dalpha,Emiss_133um)
 endif

 !--- 11 - 136 um
 if (index(Acha_Mode_Flag,'136') > 0) then
   call BTD_FM(Input%Chan_Idx_136um,  &
               Beta_110um_136um_Coef_Water, &
               Beta_110um_136um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_136um, alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_136um, Trans_Ac_136um, Trans_Bc_136um, Rad_Clear_136um, &
               f_Btd_110_136,dBtd_110_136_dTc,dBtd_110_136_dec,dBtd_110_136_dbeta, &
               dBtd_110_136_dTs,dBtd_110_136_dalpha,Emiss_136um)
 endif

 !--- 11 - 139 um
 if (index(Acha_Mode_Flag,'139') > 0) then
   call BTD_FM(Input%Chan_Idx_139um,  &
               Beta_110um_139um_Coef_Water, &
               Beta_110um_139um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_139um, alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_139um, Trans_Ac_139um, Trans_Bc_139um, Rad_Clear_139um, &
               f_Btd_110_139,dBtd_110_139_dTc,dBtd_110_139_dec,dBtd_110_139_dbeta, &
               dBtd_110_139_dTs,dBtd_110_139_dalpha,Emiss_139um)
 endif

 !--- 11 - 142 um
 if (index(Acha_Mode_Flag,'142') > 0) then
   call BTD_FM(Input%Chan_Idx_142um,  &
               Beta_110um_142um_Coef_Water, &
               Beta_110um_142um_Coef_Ice, &
               Tc,Emiss_110um, Beta_110um_120um,Ts,Tc_Base,Emiss_Sfc_142um, alpha, &
               f_T_110, dT_110_dTc, dT_110_dec, dT_110_dTs, &
               Rad_Ac_142um, Trans_Ac_142um, Trans_Bc_142um, Rad_Clear_142um, &
               f_Btd_110_142,dBtd_110_142_dTc,dBtd_110_142_dec,dBtd_110_142_dbeta, &
               dBtd_110_142_dTs,dBtd_110_142_dalpha,Emiss_142um)
 endif

 !----------------------------------------------------------------------------------------------
 ! Fill in the Kernel Matrix
 !----------------------------------------------------------------------------------------------
  f(1) = f_T_110
  K(1,1) = dT_110_dTc     
  K(1,2) = dT_110_dec    
  K(1,3) = dT_110_dbeta
  if (nx > 3) K(1,4) = dT_110_dTs
  if (nx > 4) K(1,5) = dT_110_dalpha
  select case(trim(Acha_Mode_Flag))
     case('038_110')  !11,38
        f(2) = f_Btd_110_038
        K(2,1) = dBtd_110_038_dTc     
        K(2,2) = dBtd_110_038_dec    
        K(2,3) = dBtd_110_038_dbeta
        if (nx > 3) K(2,4) = dBtd_110_038_dTs
        if (nx > 4) K(2,5) = dBtd_110_038_dalpha
     case('067_110')  !11,67
        f(2) = f_Btd_110_067
        K(2,1) = dBtd_110_067_dTc     
        K(2,2) = dBtd_110_067_dec    
        K(2,3) = dBtd_110_067_dbeta
        if (nx > 3) K(2,4) = dBtd_110_067_dTs
        if (nx > 4) K(2,5) = dBtd_110_067_dalpha
     case('110_120')  !11,12
        f(2) = f_Btd_110_120
        K(2,1) = dBtd_110_120_dTc     
        K(2,2) = dBtd_110_120_dec    
        K(2,3) = dBtd_110_120_dbeta
        if (nx > 3) K(2,4) = dBtd_110_120_dTs
        if (nx > 4) K(2,5) = dBtd_110_120_dalpha
     case('110_133') !11,13.3
        f(2) = f_Btd_110_133
        K(2,1) = dBtd_110_133_dTc     
        K(2,2) = dBtd_110_133_dec    
        K(2,3) = dBtd_110_133_dbeta
        if (nx > 3) K(2,4) = dBtd_110_133_dTs
        if (nx > 4) K(2,5) = dBtd_110_133_dalpha
     case('085_110_120') !11,12,8.5
        f(2) = f_Btd_110_085
        f(3) = f_Btd_110_120
        K(2,1) = dBtd_110_085_dTc     
        K(2,2) = dBtd_110_085_dec    
        K(2,3) = dBtd_110_085_dbeta
        if (nx > 3) K(2,4) = dBtd_110_085_dTs
        if (nx > 4) K(2,5) = dBtd_110_085_dalpha
        K(3,1) = dBtd_110_120_dTc     
        K(3,2) = dBtd_110_120_dec    
        K(3,3) = dBtd_110_120_dbeta
        if (nx > 3) K(3,4) = dBtd_110_120_dTs
        if (nx > 4) K(3,5) = dBtd_110_120_dalpha
     case('067_110_120') !11,12,6.7
        f(2) = f_Btd_110_067
        f(3) = f_Btd_110_120
        K(2,1) = dBtd_110_067_dTc     
        K(2,2) = dBtd_110_067_dec    
        K(2,3) = dBtd_110_067_dbeta
        if (nx > 3) K(2,4) = dBtd_110_067_dTs
        if (nx > 4) K(2,5) = dBtd_110_067_dalpha
        K(3,1) = dBtd_110_120_dTc     
        K(3,2) = dBtd_110_120_dec    
        K(3,3) = dBtd_110_120_dbeta
        if (nx > 3) K(3,4) = dBtd_110_120_dTs
        if (nx > 4) K(3,5) = dBtd_110_120_dalpha
     case('067_110_133') !11,13.3,6.7
        f(2) = f_Btd_110_067
        f(3) = f_Btd_110_133
        K(2,1) = dBtd_110_067_dTc     
        K(2,2) = dBtd_110_067_dec    
        K(2,3) = dBtd_110_067_dbeta
        if (nx > 3) K(2,4) = dBtd_110_067_dTs
        if (nx > 4) K(2,5) = dBtd_110_067_dalpha
        K(3,1) = dBtd_110_133_dTc     
        K(3,2) = dBtd_110_133_dec    
        K(3,3) = dBtd_110_133_dbeta
        if (nx > 3) K(3,4) = dBtd_110_133_dTs
        if (nx > 4) K(3,5) = dBtd_110_133_dalpha
     case('110_120_133') !11,12,13.3
        f(2) = f_Btd_110_120
        f(3) = f_Btd_110_133
        K(2,1) = dBtd_110_120_dTc     
        K(2,2) = dBtd_110_120_dec    
        K(2,3) = dBtd_110_120_dbeta
        if (nx > 3) K(2,4) = dBtd_110_120_dTs
        if (nx > 4) K(2,5) = dBtd_110_120_dalpha
        K(3,1) = dBtd_110_133_dTc     
        K(3,2) = dBtd_110_133_dec    
        K(3,3) = dBtd_110_133_dbeta
        if (nx > 3) K(3,4) = dBtd_110_133_dTs
        if (nx > 4) K(3,5) = dBtd_110_133_dalpha
     case('067_085_110') !11,8.5,6.7
        f(2) = f_Btd_110_067
        f(3) = f_Btd_110_085
        K(2,1) = dBtd_110_067_dTc     
        K(2,2) = dBtd_110_067_dec    
        K(2,3) = dBtd_110_067_dbeta
        if (nx > 3) K(2,4) = dBtd_110_067_dTs
        if (nx > 4) K(2,5) = dBtd_110_067_dalpha
        K(3,1) = dBtd_110_085_dTc     
        K(3,2) = dBtd_110_085_dec    
        K(3,3) = dBtd_110_085_dbeta
        if (nx > 3) K(3,4) = dBtd_110_085_dTs
        if (nx > 4) K(3,5) = dBtd_110_085_dalpha
     case('085_110_120_133') !11,12,8.5,13.3
        f(2) = f_Btd_110_085
        f(3) = f_Btd_110_120
        f(4) = f_Btd_110_133
        K(2,1) = dBtd_110_085_dTc     
        K(2,2) = dBtd_110_085_dec    
        K(2,3) = dBtd_110_085_dbeta
        if (nx > 3) K(2,4) = dBtd_110_085_dTs
        if (nx > 4) K(2,5) = dBtd_110_085_dalpha
        K(3,1) = dBtd_110_120_dTc     
        K(3,2) = dBtd_110_120_dec    
        K(3,3) = dBtd_110_120_dbeta
        if (nx > 3) K(3,4) = dBtd_110_120_dTs
        if (nx > 4) K(3,5) = dBtd_110_120_dalpha
        K(4,1) = dBtd_110_133_dTc     
        K(4,2) = dBtd_110_133_dec    
        K(4,3) = dBtd_110_133_dbeta
        if (nx > 3) K(4,4) = dBtd_110_133_dTs
        if (nx > 4) K(4,5) = dBtd_110_133_dalpha
     case('067_085_110_120') !11,12,8.5,6.7
        f(2) = f_Btd_110_067
        f(3) = f_Btd_110_085
        f(4) = f_Btd_110_120
        K(2,1) = dBtd_110_067_dTc     
        K(2,2) = dBtd_110_067_dec    
        K(2,3) = dBtd_110_067_dbeta
        if (nx > 3) K(2,4) = dBtd_110_067_dTs
        if (nx > 4) K(2,5) = dBtd_110_067_dalpha
        K(3,1) = dBtd_110_085_dTc     
        K(3,2) = dBtd_110_085_dec    
        K(3,3) = dBtd_110_085_dbeta
        if (nx > 3) K(3,4) = dBtd_110_085_dTs
        if (nx > 4) K(3,5) = dBtd_110_085_dalpha
        K(4,1) = dBtd_110_120_dTc     
        K(4,2) = dBtd_110_120_dec    
        K(4,3) = dBtd_110_120_dbeta
        if (nx > 3) K(4,4) = dBtd_110_120_dTs
        if (nx > 4) K(4,5) = dBtd_110_120_dalpha
     case('062_085_110_120_133') !11,12,8.5,13.3,6.7
        f(2) = f_Btd_110_062
        f(3) = f_Btd_110_085
        f(4) = f_Btd_110_120
        f(5) = f_Btd_110_133
        K(2,1) = dBtd_110_062_dTc     
        K(2,2) = dBtd_110_062_dec    
        K(2,3) = dBtd_110_062_dbeta
        if (nx > 3) K(2,4) = dBtd_110_062_dTs
        if (nx > 4) K(2,5) = dBtd_110_062_dalpha
        K(3,1) = dBtd_110_085_dTc     
        K(3,2) = dBtd_110_085_dec    
        K(3,3) = dBtd_110_085_dbeta
        if (nx > 3) K(3,4) = dBtd_110_085_dTs
        if (nx > 4) K(3,5) = dBtd_110_085_dalpha
        K(4,1) = dBtd_110_120_dTc     
        K(4,2) = dBtd_110_120_dec    
        K(4,3) = dBtd_110_120_dbeta
        if (nx > 3) K(4,4) = dBtd_110_120_dTs
        if (nx > 4) K(4,5) = dBtd_110_120_dalpha
        K(5,1) = dBtd_110_133_dTc     
        K(5,2) = dBtd_110_133_dec    
        K(5,3) = dBtd_110_133_dbeta
        if (nx > 3) K(5,4) = dBtd_110_133_dTs
        if (nx > 4) K(5,5) = dBtd_110_133_dalpha
     case('067_085_110_120_133') !11,12,8.5,13.3,6.7
        f(2) = f_Btd_110_067
        f(3) = f_Btd_110_085
        f(4) = f_Btd_110_120
        f(5) = f_Btd_110_133
        K(2,1) = dBtd_110_067_dTc     
        K(2,2) = dBtd_110_067_dec    
        K(2,3) = dBtd_110_067_dbeta
        if (nx > 3) K(2,4) = dBtd_110_067_dTs
        if (nx > 4) K(2,5) = dBtd_110_067_dalpha
        K(3,1) = dBtd_110_085_dTc     
        K(3,2) = dBtd_110_085_dec    
        K(3,3) = dBtd_110_085_dbeta
        if (nx > 3) K(3,4) = dBtd_110_085_dTs
        if (nx > 4) K(3,5) = dBtd_110_085_dalpha
        K(4,1) = dBtd_110_120_dTc     
        K(4,2) = dBtd_110_120_dec    
        K(4,3) = dBtd_110_120_dbeta
        if (nx > 3) K(4,4) = dBtd_110_120_dTs
        if (nx > 4) K(4,5) = dBtd_110_120_dalpha
        K(5,1) = dBtd_110_133_dTc     
        K(5,2) = dBtd_110_133_dec    
        K(5,3) = dBtd_110_133_dbeta
        if (nx > 3) K(5,4) = dBtd_110_133_dTs
        if (nx > 4) K(5,5) = dBtd_110_133_dalpha
     case('110_133_136_139_142') !11,13.3
        f(2) = f_Btd_110_133
        f(3) = f_Btd_110_136
        f(4) = f_Btd_110_139
        f(5) = f_Btd_110_142
        K(2,1) = dBtd_110_133_dTc     
        K(2,2) = dBtd_110_133_dec    
        K(2,3) = dBtd_110_133_dbeta
        if (nx > 3) K(2,4) = dBtd_110_133_dTs
        if (nx > 4) K(2,5) = dBtd_110_133_dalpha
        K(3,1) = dBtd_110_136_dTc     
        K(3,2) = dBtd_110_136_dec    
        K(3,3) = dBtd_110_136_dbeta
        if (nx > 3) K(3,4) = dBtd_110_136_dTs
        if (nx > 4) K(3,5) = dBtd_110_136_dalpha
        K(4,1) = dBtd_110_139_dTc     
        K(4,2) = dBtd_110_139_dec    
        K(4,3) = dBtd_110_139_dbeta
        if (nx > 3) K(4,4) = dBtd_110_139_dTs
        if (nx > 4) K(4,5) = dBtd_110_139_dalpha
        K(5,1) = dBtd_110_142_dTc     
        K(5,2) = dBtd_110_142_dec    
        K(5,3) = dBtd_110_142_dbeta
        if (nx > 3) K(5,4) = dBtd_110_142_dTs
        if (nx > 4) K(5,5) = dBtd_110_142_dalpha
     case ('085_110_120_133_136_139_142')
        f(2) = f_Btd_110_085
        f(3) = f_Btd_110_120
        f(4) = f_Btd_110_133
        f(5) = f_Btd_110_136
        f(6) = f_Btd_110_139
        f(7) = f_Btd_110_142
        K(2,1) = dBtd_110_085_dTc     
        K(2,2) = dBtd_110_085_dec    
        K(2,3) = dBtd_110_085_dbeta
        if (nx > 3) K(2,4) = dBtd_110_085_dTs
        if (nx > 4) K(2,5) = dBtd_110_085_dalpha
        K(3,1) = dBtd_110_120_dTc     
        K(3,2) = dBtd_110_120_dec    
        K(3,3) = dBtd_110_120_dbeta
        if (nx > 3) K(3,4) = dBtd_110_120_dTs
        if (nx > 4) K(3,5) = dBtd_110_120_dalpha
        K(4,1) = dBtd_110_133_dTc     
        K(4,2) = dBtd_110_133_dec    
        K(4,3) = dBtd_110_133_dbeta
        if (nx > 3) K(4,4) = dBtd_110_133_dTs
        if (nx > 4) K(4,5) = dBtd_110_133_dalpha
        K(5,1) = dBtd_110_136_dTc     
        K(5,2) = dBtd_110_136_dec    
        K(5,3) = dBtd_110_136_dbeta
        if (nx > 3) K(5,4) = dBtd_110_136_dTs
        if (nx > 4) K(5,5) = dBtd_110_136_dalpha
        K(6,1) = dBtd_110_139_dTc     
        K(6,2) = dBtd_110_139_dec    
        K(6,3) = dBtd_110_139_dbeta
        if (nx > 3) K(6,4) = dBtd_110_139_dTs
        if (nx > 4) K(6,5) = dBtd_110_139_dalpha
        K(7,1) = dBtd_110_142_dTc     
        K(7,2) = dBtd_110_142_dec    
        K(7,3) = dBtd_110_142_dbeta
        if (nx > 3) K(7,4) = dBtd_110_142_dTs
        if (nx > 4) K(7,5) = dBtd_110_142_dalpha
     case("062_067_073_085_104_110_120_133")
        f(2) = f_Btd_110_062
        f(3) = f_Btd_110_067
        f(4) = f_Btd_110_073
        f(5) = f_Btd_110_085
        f(6) = f_Btd_110_104
        f(7) = f_Btd_110_120
        f(8) = f_Btd_110_133
        K(2,1) = dBtd_110_062_dTc     
        K(2,2) = dBtd_110_062_dec    
        K(2,3) = dBtd_110_062_dbeta
        if (nx > 3) K(2,4) = dBtd_110_062_dTs
        if (nx > 4) K(2,5) = dBtd_110_062_dalpha
        K(3,1) = dBtd_110_067_dTc     
        K(3,2) = dBtd_110_067_dec    
        K(3,3) = dBtd_110_067_dbeta
        if (nx > 3) K(3,4) = dBtd_110_067_dTs
        if (nx > 4) K(3,5) = dBtd_110_067_dalpha
        K(4,1) = dBtd_110_073_dTc     
        K(4,2) = dBtd_110_073_dec    
        K(4,3) = dBtd_110_073_dbeta
        if (nx > 3) K(4,4) = dBtd_110_073_dTs
        if (nx > 4) K(4,5) = dBtd_110_073_dalpha
        K(5,1) = dBtd_110_085_dTc     
        K(5,2) = dBtd_110_085_dec    
        K(5,3) = dBtd_110_085_dbeta
        if (nx > 3) K(5,4) = dBtd_110_085_dTs
        if (nx > 4) K(5,5) = dBtd_110_085_dalpha
        K(6,1) = dBtd_110_104_dTc     
        K(6,2) = dBtd_110_104_dec    
        K(6,3) = dBtd_110_104_dbeta
        if (nx > 3) K(6,4) = dBtd_110_104_dTs
        if (nx > 4) K(6,5) = dBtd_110_104_dalpha
        K(7,1) = dBtd_110_120_dTc     
        K(7,2) = dBtd_110_120_dec    
        K(7,3) = dBtd_110_120_dbeta
        if (nx > 3) K(7,4) = dBtd_110_120_dTs
        if (nx > 4) K(7,5) = dBtd_110_120_dalpha
        K(8,1) = dBtd_110_133_dTc     
        K(8,2) = dBtd_110_133_dec    
        K(8,3) = dBtd_110_133_dbeta
        if (nx > 3) K(8,4) = dBtd_110_133_dTs
        if (nx > 4) K(8,5) = dBtd_110_133_dalpha
  end select

 !--- determine number of channels
  select case(trim(Acha_Mode_Flag))
     case('110')  !avhrr, goes-im
       Emiss_Vector(1) = Emiss_110um
     case('038_110')  !goes-np 3 chan
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_038um
     case('067_110')  !goes-np 3 chan
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_067um
     case('110_120')  !avhrr, goes-im
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_120um
     case('110_133')  !goes-nop
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_133um
     case('085_110_120')  !viirs
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_085um
       Emiss_Vector(3) = Emiss_120um
     case('067_110_120')  !goes-im 3 chan
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_067um
       Emiss_Vector(3) = Emiss_120um
     case('067_110_133')  !goes-np 3 chan
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_067um
       Emiss_Vector(3) = Emiss_133um
     case('110_120_133')  !goes-r
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_120um
       Emiss_Vector(3) = Emiss_133um
     case('067_085_110')  !goes-r
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_067um
       Emiss_Vector(3) = Emiss_085um
     case('085_110_120_133')  
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_085um
       Emiss_Vector(3) = Emiss_120um
       Emiss_Vector(4) = Emiss_133um
     case('067_085_110_120')  
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_067um
       Emiss_Vector(3) = Emiss_085um
       Emiss_Vector(4) = Emiss_120um
     case('062_085_110_120_133')  
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_062um
       Emiss_Vector(3) = Emiss_085um
       Emiss_Vector(4) = Emiss_120um
       Emiss_Vector(5) = Emiss_133um
     case('067_085_110_120_133')  
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_067um
       Emiss_Vector(3) = Emiss_085um
       Emiss_Vector(4) = Emiss_120um
       Emiss_Vector(5) = Emiss_133um
     case('110_133_136_139_142')  !goes-nop
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_133um
       Emiss_Vector(3) = Emiss_136um
       Emiss_Vector(4) = Emiss_139um
       Emiss_Vector(5) = Emiss_142um
     case ('085_110_120_133_136_139_142')
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_085um
       Emiss_Vector(3) = Emiss_120um
       Emiss_Vector(4) = Emiss_133um
       Emiss_Vector(5) = Emiss_136um
       Emiss_Vector(6) = Emiss_139um
       Emiss_Vector(7) = Emiss_142um
     case("062_067_073_085_104_110_120_133")
       Emiss_Vector(1) = Emiss_110um
       Emiss_Vector(2) = Emiss_062um
       Emiss_Vector(3) = Emiss_067um
       Emiss_Vector(4) = Emiss_073um
       Emiss_Vector(5) = Emiss_085um
       Emiss_Vector(6) = Emiss_104um
       Emiss_Vector(7) = Emiss_120um
       Emiss_Vector(8) = Emiss_133um
  end select

end subroutine COMPUTE_FORWARD_MODEL_AND_KERNEL

end module ACHA_FULL_RETRIEVAL_MOD
