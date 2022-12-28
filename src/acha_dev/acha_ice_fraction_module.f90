!$Id:$
module ACHA_ICE_FRACTION_MODULE

 use ACHA_SERVICES_MOD, only : &
           real4, int1, int4, real8, dtor, &
           acha_output_struct,ACHA_SYMBOL_STRUCT, &
           MISSING_VALUE_REAL4
 use CX_REAL_BOOLEAN_MOD

 implicit none

 integer, private, parameter:: ntopa_lut = 13
 real, private, parameter:: topa_min_lut = 190
 real, private, parameter:: topa_bin_lut = 10

 real, private, dimension(ntopa_lut), parameter:: ice_prob_lut =  &
      [0.990,0.999,0.997,0.990,0.954,0.736,0.538,0.421,0.364,0.359,0.578,0.748,0.891]

 real, private, dimension(ntopa_lut), parameter:: ice_prob_uncer_lut =  &
      [0.010,0.001,0.003,0.010,0.046,0.264,0.462,0.421,0.364,0.359,0.422,0.252,0.109]

 private

 public:: COMPUTE_ICE_FRACTION_FROM_TYPE
 public:: COMPUTE_ICE_FRACTION_FROM_LUT1D

 contains

!----------------------------------------------------------------------
! Setup Ice Fraction Lookup Table
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! Compute Ice Fraction from  1d Lookup Table based on opaque temperature
!----------------------------------------------------------------------
subroutine COMPUTE_ICE_FRACTION_FROM_LUT1D(Topa, Ice_Fraction, Ice_Fraction_Uncertainty)
   real, intent(in):: Topa
   real, intent(out):: Ice_Fraction, Ice_Fraction_Uncertainty
   integer:: Topa_Idx

   Ice_Fraction = MISSING_VALUE_REAL4
   Ice_Fraction_Uncertainty = MISSING_VALUE_REAL4
   if (Topa .eqr. MISSING_VALUE_REAL4) return

   Topa_Idx = int((Topa - Topa_Min_Lut) / Topa_Bin_Lut) 
   Topa_Idx = min(max(1,Topa_Idx),Ntopa_Lut)
   Ice_Fraction = Ice_Prob_Lut(Topa_Idx)
   Ice_Fraction_Uncertainty = Ice_Prob_Uncer_Lut(Topa_Idx)
   
end subroutine COMPUTE_ICE_FRACTION_FROM_LUT1D


!----------------------------------------------------------------------
! Compute Ice Fraction from the Cloud Type
!----------------------------------------------------------------------
subroutine COMPUTE_ICE_FRACTION_FROM_TYPE(Cloud_Type, Symbol, Ice_Fraction, Ice_Fraction_Uncertainty)
   integer(kind=int1), intent(in):: Cloud_Type
   real, intent(out):: Ice_Fraction, Ice_Fraction_Uncertainty
   type(ACHA_SYMBOL_STRUCT), intent(in) :: Symbol

   if (Cloud_Type == Symbol%FOG_TYPE) then
         Ice_Fraction = 0.0
         Ice_Fraction_Uncertainty = 0.2
   elseif (Cloud_Type == Symbol%WATER_TYPE) then
         Ice_Fraction = 0.0
         Ice_Fraction_Uncertainty = 0.2
   elseif (Cloud_Type == Symbol%SUPERCOOLED_TYPE) then
         Ice_Fraction = 0.0
         Ice_Fraction_Uncertainty = 0.2
   elseif (Cloud_Type == Symbol%MIXED_TYPE) then
         Ice_Fraction = 0.0
         Ice_Fraction_Uncertainty = 0.2
   elseif (Cloud_Type == Symbol%OPAQUE_ICE_TYPE) then
         Ice_Fraction = 1.0
         Ice_Fraction_Uncertainty = 0.1
   elseif (Cloud_Type == Symbol%CIRRUS_TYPE) then
         Ice_Fraction = 1.0
         Ice_Fraction_Uncertainty = 0.4
   elseif (Cloud_Type == Symbol%OVERLAP_TYPE) then
         Ice_Fraction = 1.0
         Ice_Fraction_Uncertainty = 0.3
   elseif (Cloud_Type == Symbol%OVERSHOOTING_TYPE) then
         Ice_Fraction = 1.0
         Ice_Fraction_Uncertainty = 0.1
   else
         Ice_Fraction = 0.0
         Ice_Fraction_Uncertainty = 0.2
   endif

!  Ice_Fraction = 0.5
!  Ice_Fraction_Uncertainty = 0.5

end subroutine COMPUTE_ICE_FRACTION_FROM_TYPE

end module ACHA_ICE_FRACTION_MODULE
