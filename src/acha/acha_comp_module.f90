! $Id: acha_comp_module.f90 4033 2020-10-15 19:14:57Z heidinger $
module ACHA_COMP
!---------------------------------------------------------------------
! This module houses the routines associated with...
!
! ACHA Cloud Optical Microphysical Properties
!
! Author: 
!
! Reference:
!
!----------------------------------------------------------------------
  use ACHA_SERVICES_MOD, only : &
           real4, int1, int4, &
           acha_output_struct,acha_symbol_struct, &
           acha_input_struct

  use ACHA_MICROPHYSICAL_MODULE

  implicit none

  public:: ACHA_COMP_ALGORITHM

  private:: COMPUTE_TAU_REFF_ACHA 

  !--- include the non-system specific variables
  include 'include/acha_parameters.inc'

  real, private, parameter:: MISSING_VALUE_REAL = -999.0
  integer, private, parameter:: MISSING_VALUE_INTEGER = -999
  type(acha_symbol_struct), private :: symbol

  contains 

!------------------------------------------------------------------------------
! AWG Cloud Optical and Microphysical Algorithm (ACHA-COMP)
!
! Author: Andrew Heidinger, NOAA
!
! Assumptions
!
! Limitations
!
! NOTE.  This algorithm use the same input and output structures as 
!        the AWG_CLOUD_HEIGHT_ALGORITHM.
!        Do not overwrite elements of the Output structure expect those
!        generated here.
!
!      Output%Tau
!      Output%Tau_Uncer
!      Output%Ec  (modified from ACHA)
!      Output%Reff
!      Output%Beta (modified from ACHA)
!
!----------------------------------------------------------------------
! modification history
!
!------------------------------------------------------------------------------
  subroutine ACHA_COMP_ALGORITHM(Input, Symbol_In, Output)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================

  type(acha_input_struct), intent(inout) :: Input
  type(acha_symbol_struct), intent(in) :: Symbol_In
  type(acha_output_struct), intent(inout) :: Output

  !===============================================================================
  !  Local Variable Declaration
  !===============================================================================

  integer:: Elem_Idx
  integer:: Line_Idx

  integer:: Cloud_Type
  integer:: Cloud_Phase

  !--- scalar local variables
  real (kind=real4):: Ec_Slant

!-----------------------------------------------------------------------
! BEGIN EXECUTABLE CODE
!-----------------------------------------------------------------------

   !-------------------------------------------------------------------------
   ! Initialization
   !-------------------------------------------------------------------------
   symbol = symbol_in   !symbol is a module-wide variable

  !---------------------------------------------------------------------------
  !-- setup microphysical models
  !---------------------------------------------------------------------------
  call SETUP_ICE_MICROPHYSICAL_MODEL(Input%WMO_Id)

   !--------------------------------------------------------------------------
   ! loop over pixels in scanlines
   !--------------------------------------------------------------------------
   Line_loop: do Line_Idx = Line_Idx_min,Input%Number_of_Lines + Line_Idx_min - 1

    Element_Loop:   do Elem_Idx = 1, Input%Number_of_Elements

    !--- for convenience, save nwp indices to local variables
    Cloud_Type = Input%Cloud_Type(Elem_Idx,Line_Idx)

    !----------------------------------------------------------------------
    ! for clear pixels, set Opd to zero and Reff to missing
    !----------------------------------------------------------------------
    if ((Cloud_Type == symbol%CLEAR_TYPE .or. &
        Cloud_Type == symbol%PROB_CLEAR_TYPE) .and. &
        Input%Sensor_Zenith_Angle(Elem_Idx,Line_Idx) <= Sensor_Zenith_Threshold) then
        Output%Ec(Elem_Idx,Line_Idx) = 0.0
        Output%Tau(Elem_Idx,Line_Idx) = 0.0
        Output%Tau_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL
        Output%Reff(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL
        Output%Beta(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL
       cycle
    endif
 
   
    !----------------------------------------------------------------------
    ! determine cloud phase from cloud type for convenience
    !----------------------------------------------------------------------
    Cloud_Phase = symbol%UNKNOWN_PHASE

    if ( (Cloud_Type  == symbol%FOG_TYPE) .or. &
       (Cloud_Type  == symbol%WATER_TYPE) .or. &
       (Cloud_Type  == symbol%SUPERCOOLED_TYPE)) then

       Cloud_Phase = symbol%WATER_PHASE

    endif

    if ( (Cloud_Type  == symbol%CIRRUS_TYPE) .or. &
        (Cloud_Type  == symbol%OVERLAP_TYPE) .or. &
        (Cloud_Type  == symbol%TICE_TYPE) .or. &
        (Cloud_Type  == symbol%OVERSHOOTING_TYPE)) then

        Cloud_Phase = symbol%ICE_PHASE

    endif

    !-----------------------------------------------------------------------------
    ! Estimate Cloud Optical and Microphysical Properties
    !-----------------------------------------------------------------------------
    Ec_Slant =  Output%Ec(Elem_Idx,Line_Idx)

    if (Output%Zc(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL .and. &
        Output%Ec(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL .and. &
        Output%Beta(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL) then


     !--- save nadir adjusted emissivity and optical depth
     if (Output%Ec(Elem_Idx,Line_Idx) < 1.00) then

       call COMPUTE_TAU_REFF_ACHA(&
                              Output%Beta(Elem_Idx,Line_Idx), &
                              Input%Cosine_Zenith_Angle(Elem_Idx,Line_Idx), &
                              Cloud_Phase, &
                              Ec_Slant, &
                              Output%Ec(Elem_Idx,Line_Idx), &
                              Output%Ec_Uncertainty(Elem_Idx,Line_Idx), &
                              Output%Tau(Elem_Idx,Line_Idx), &
                              Output%Tau_Uncertainty(Elem_Idx,Line_Idx), &
                              Output%Reff(Elem_Idx,Line_Idx))

     else

       Output%Tau(Elem_Idx,Line_Idx) = 20.0
       Output%Tau_Uncertainty(Elem_Idx,Line_Idx) = 20.0
       Output%Ec(Elem_Idx,Line_Idx) = 1.0

       if( Cloud_Phase == symbol%ICE_PHASE) then
          Output%Beta(Elem_Idx,Line_Idx) = Beta_Ap_Ice
          Output%Reff(Elem_Idx,Line_Idx) = 20.0
       else
          Output%Beta(Elem_Idx,Line_Idx) = Beta_Ap_Water
          Output%Reff(Elem_Idx,Line_Idx) = 10.0
       endif

     endif

    endif

 end do Element_Loop

end do Line_Loop

end subroutine  ACHA_COMP_ALGORITHM

!---------------------------------------------------------------------------
! routine to compute optical depth and eff radius from ec and beta
!---------------------------------------------------------------------------
subroutine COMPUTE_TAU_REFF_ACHA(Beta, &
                                 Cosine_Zenith_Angle, &
                                 Cloud_Phase, &
                                 Ec_Slant, & 
                                 Ec, &
                                 Ec_Uncer, &
                                 Tau, &
                                 Tau_Uncer, &
                                 Reff)

   real(kind=real4), intent(in):: Beta
   real(kind=real4), intent(in):: Cosine_Zenith_Angle
   real(kind=real4), intent(in):: Ec_Slant
   integer(kind=int4), intent(in):: Cloud_Phase
   real(kind=real4), intent(out):: Ec
   real(kind=real4), intent(in):: Ec_Uncer
   real(kind=real4), intent(out):: Tau
   real(kind=real4), intent(out):: Tau_Uncer
   real(kind=real4), intent(out):: Reff
   real(kind=real4):: Qe_vis
   real(kind=real4):: Qe_110um
   real(kind=real4):: wo_110um
   real(kind=real4):: g_110um
   real(kind=real4):: Tau_Abs_110um
   real(kind=real4):: Temp_R4
   real(kind=real4):: log10_Reff
   real(kind=real4):: Factor
   real(kind=real4), parameter:: Reff_Min = 1.0
   real(kind=real4), parameter:: Reff_Max = 60.0
   real(kind=real4), parameter:: Tau_Max = 8.0

   Tau = MISSING_VALUE_REAL
   Reff = MISSING_VALUE_REAL
   Ec = MISSING_VALUE_REAL

   if (Cloud_Phase == symbol%ICE_PHASE) then
    Temp_R4 = Re_Beta_110um_COEF_ICE(0) +  &
              Re_Beta_110um_COEF_ICE(1)*(Beta-1.0) + &
              Re_Beta_110um_COEF_ICE(2)*(Beta-1.0)**2 + &
              Re_Beta_110um_COEF_ICE(3)*(Beta-1.0)**3
   else
    Temp_R4 = Re_Beta_110um_COEF_WATER(0) +  &
              Re_Beta_110um_COEF_WATER(1)*(Beta-1.0) + &
              Re_Beta_110um_COEF_WATER(2)*(Beta-1.0)**2 + &
              Re_Beta_110um_COEF_WATER(3)*(Beta-1.0)**3
   endif

   Reff = max(Reff_Min,min(Reff_Max,10.0**(1.0/Temp_R4)))  !note inverse here

   if (Reff > 0.0) then
     log10_Reff = alog10(Reff)
   else
     return
   endif
   !--- determine optical depth and its uncertainty
   if (Ec_Slant >= 0.0 .and. Ec_Slant < 1.0) then

      !--- determine single scattering properties
      if (Cloud_Phase == symbol%ICE_PHASE) then

        Qe_Vis = Qe_006um_COEF_ICE(0) +  &
                 Qe_006um_COEF_ICE(1)*log10_Reff + &
                 Qe_006um_COEF_ICE(2)*log10_Reff**2

        Qe_110um = Qe_110um_COEF_ICE(0) +  &
                   Qe_110um_COEF_ICE(1)*log10_Reff + &
                   Qe_110um_COEF_ICE(2)*log10_Reff**2

        wo_110um = wo_110um_COEF_ICE(0) +  &
                   wo_110um_COEF_ICE(1)*log10_Reff + &
                   wo_110um_COEF_ICE(2)*log10_Reff**2

        g_110um = g_110um_COEF_ICE(0) +  &
                  g_110um_COEF_ICE(1)*log10_Reff + &
                  g_110um_COEF_ICE(2)*log10_Reff**2

      else

        Qe_Vis = Qe_006um_COEF_WATER(0) +  &
                 Qe_006um_COEF_WATER(1)*log10_Reff + &
                 Qe_006um_COEF_WATER(2)*log10_Reff**2

        Qe_110um = Qe_110um_COEF_WATER(0) +  &
                   Qe_110um_COEF_WATER(1)*log10_Reff + &
                   Qe_110um_COEF_WATER(2)*log10_Reff**2

        wo_110um = wo_110um_COEF_WATER(0) +  &
                   wo_110um_COEF_WATER(1)*log10_Reff + &
                   wo_110um_COEF_WATER(2)*log10_Reff**2

        g_110um =  g_110um_COEF_WATER(0) +  &
                   g_110um_COEF_WATER(1)*log10_Reff + &
                   g_110um_COEF_WATER(2)*log10_Reff**2

      endif

      Tau_Abs_110um = -Cosine_Zenith_Angle*alog(1.0 - Ec_Slant) 

      Ec = 1.0 - exp(-Tau_Abs_110um)
   
      Factor = (Qe_vis / Qe_110um) / (1.0 - wo_110um * g_110um)

      Tau = min(Factor * Tau_Abs_110um, Tau_Max)

      Tau_Uncer = Ec_Uncer * Factor / (1.0-Ec)

      !--- set negative values to be missing 
      if (Tau < 0) then
         Tau = MISSING_VALUE_REAL
         Tau_Uncer = MISSING_VALUE_REAL
         Reff= MISSING_VALUE_REAL
      endif

   endif

end subroutine COMPUTE_TAU_REFF_ACHA 

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module ACHA_COMP
