!$Id: acha_module.f90 3876 2020-06-18 13:34:40Z yli $
module ACHA_LHP_MOD
!---------------------------------------------------------------------
! ACHA Loop Heat Pipe Module
!----------------------------------------------------------------------
  use ACHA_SERVICES_MOD, only : &
           Acha_output_struct,ACHA_SYMBOL_STRUCT, &
           Acha_input_struct, Acha_rtm_nwp_struct, &
           ACHA_FETCH_PIXEL_NWP_RTM, &
           ABI_Use_104um_Flag

  use LEVEL2_STRUCTURES_MOD, only: &
      Clavrx_Global_Attr

  implicit none

  public:: MODIFY_MODE_USING_LHP_THRESHOLDS
  public:: SETUP_REFERENCE_CHANNEL
  public:: SETUP_REFERENCE_CHANNEL_PROFILES

  !--- include the non-system specific variables
! include 'acha_parameters.inc'

  contains 
!---------------------------------------------------------------------------
! Determine ACHA Mode based on available channels.  These may have been
! modified with a call to SET_ABI_USE_104um_FLAG from process_clavrx.
!----------------------------------------------------------------------------
subroutine MODIFY_MODE_USING_LHP_THRESHOLDS(Acha_Mode_Flag,Chan_On_038, &
                                            Chan_On_067, Chan_On_085, &
                                            Chan_On_104, Chan_On_110, &
                                            Chan_On_120, Chan_On_133)

  character(len=*), intent(inout):: Acha_Mode_Flag
  integer, intent(in):: Chan_On_038
  integer, intent(in):: Chan_On_067
  integer, intent(in):: Chan_On_085
  integer, intent(in):: Chan_On_104
  integer, intent(in):: Chan_On_110
  integer, intent(in):: Chan_On_120
  integer, intent(in):: Chan_On_133
  character(len=len(Acha_Mode_Flag)):: Acha_Mode_Flag_Input

  !--- Store initial ACHA Mode.
  Acha_Mode_Flag_Input = Acha_Mode_Flag

  !--- See if user requested ACHA mode can be done.
  if ( (index(Acha_Mode_Flag,'038') > 0 .and. Chan_On_038 == 0 ) .or. &
       (index(Acha_Mode_Flag,'067') > 0 .and. Chan_On_067 == 0 ) .or. &
       (index(Acha_Mode_Flag,'085') > 0 .and. Chan_On_085 == 0 ) .or. &
       (index(Acha_Mode_Flag,'120') > 0 .and. Chan_On_120 == 0 ) .or. &
       (index(Acha_Mode_Flag,'133') > 0 .and. Chan_On_133 == 0 ) ) then

    Acha_Mode_Flag = "unknown"

    !--- Figure out which ACHA mode to use.
    if (Acha_Mode_Flag == "unknown" .and. &
      Chan_On_067 == 1 .and. &
      Chan_On_085 == 1 .and. &
      Chan_On_120 == 1 ) then
      Acha_Mode_Flag = "067_085_110_120"
    endif

    if (Acha_Mode_Flag == "unknown" .and. &
      Chan_On_085 == 1 .and. &
      Chan_On_120 == 1 ) then
      Acha_Mode_Flag = "085_110_120"
    endif

    if (Acha_Mode_Flag == "unknown" .and. &
      Chan_On_120 == 1 ) then
      Acha_Mode_Flag = "110_120"
    endif

    if (Acha_Mode_Flag == "unknown" .and. &
      Chan_On_067 == 1 ) then
      Acha_Mode_Flag = "067_110"
    endif

    if (Acha_Mode_Flag == "unknown" .and. &
      Chan_On_038 == 1 ) then
      Acha_Mode_Flag = "038_110"
    endif

    if (Acha_Mode_Flag == "unknown") then
      Acha_Mode_Flag = "110"
    endif

    !--- print warning to screen if a mode change occurred
    if (Acha_Mode_Flag_Input /= Acha_Mode_Flag) then
      Clavrx_Global_Attr%acha_mode = Acha_Mode_Flag
      print *, "WARNING: LHP Thresholds changed ACHA Mode from ",trim(Acha_Mode_Flag_Input), " to ", trim(Acha_Mode_Flag)
    endif

  endif

end subroutine MODIFY_MODE_USING_LHP_THRESHOLDS
!----------------------------------------------------------------------
! Setup Reference Channel
!----------------------------------------------------------------------
subroutine SETUP_REFERENCE_CHANNEL(Use_10_4,Input)

   logical, intent(in):: Use_10_4
   type(acha_input_struct), intent(inout) :: Input

   if (Use_10_4) then
     Input%Chan_Idx_110um = Input%Chan_Idx_104um
     Input%Bt_110um => Input%Bt_104um
     Input%Rad_110um => Input%Rad_104um
     Input%Rad_Clear_110um => Input%Rad_Clear_104um
     Input%Surface_Emissivity_110um => Input%Surface_Emissivity_104um
   endif

end subroutine SETUP_REFERENCE_CHANNEL
!----------------------------------------------------------------------
! Setup Reference Profiles
!----------------------------------------------------------------------
subroutine SETUP_REFERENCE_CHANNEL_PROFILES(Use_10_4,ACHA_RTM_NWP)

   logical, intent(in):: Use_10_4
   type(acha_rtm_nwp_struct), intent(inout) :: ACHA_RTM_NWP

   if (Use_10_4) then
     ACHA_RTM_NWP%Atm_Rad_Prof_110um =>  ACHA_RTM_NWP%Atm_Rad_Prof_104um
     ACHA_RTM_NWP%Atm_Trans_Prof_110um => ACHA_RTM_NWP%Atm_Trans_Prof_104um
     ACHA_RTM_NWP%Black_Body_Rad_Prof_110um => ACHA_RTM_NWP%Black_Body_Rad_Prof_104um
   endif

end subroutine SETUP_REFERENCE_CHANNEL_PROFILES
!---------------------------------------------------------------------------
! set the 10.4 flag
! For GOES-RU
! DQF (0=good, 1=ok, 2 = useful, 3 = bad, 4= horrendous)
!----------------------------------------------------------------------------
!subroutine SET_10_4_FLAG(Chan_On_104, Chan_On_110, Chan_On_120, Chan_On_133, &
!                         Dqf_104, Dqf_110, Dqf_120, Dqf_133, Use_10_4_Flag)
!   integer(kind=int4), intent(in):: Chan_On_104, Chan_On_110, Chan_On_120,
!   Chan_On_133
!   integer(kind=int1), dimension(:,:), intent(in):: Dqf_104
!   integer(kind=int1), dimension(:,:), intent(in):: Dqf_110
!   integer(kind=int1), dimension(:,:), intent(in):: Dqf_120
!   integer(kind=int1), dimension(:,:), intent(in):: Dqf_133
!   logical, intent(out):: Use_10_4_Flag
!   logical, dimension(size(Dqf_110(:,1)),size(Dqf_110(1,:))):: Mask_104
!   logical, dimension(size(Dqf_110(:,1)),size(Dqf_110(1,:))):: Mask_110
!   logical, dimension(size(Dqf_110(:,1)),size(Dqf_110(1,:))):: Mask_120
!   logical, dimension(size(Dqf_110(:,1)),size(Dqf_110(1,:))):: Mask_133
!   integer(kind=int1), parameter:: Dqf_Thresh = 1
!   integer:: Count_Good_110
!   integer:: Count_Good_104
!   integer:: Count_Good_120
!   integer:: Count_Good_133
!
!   Mask_104 = .false.
!   Mask_110 = .false.
!   Mask_120 = .false.
!   Mask_133 = .false.
!   Use_10_4_Flag = .false.
!
!   if (Chan_On_104 == 0) return
!
!   if (Chan_On_104 > 0) then
!      where ( Dqf_104 <= Dqf_Thresh .AND. Dqf_104 >= 0)
!        Mask_104 = .true.
!      endwhere
!   endif
!   if (Chan_On_110 > 0) then
!      where ( Dqf_110 <= Dqf_Thresh .AND. Dqf_110 >= 0)
!        Mask_110 = .true.
!      endwhere
!   endif
!   if (Chan_On_120 > 0) then
!      where ( Dqf_120 <= Dqf_Thresh .AND. Dqf_120 >= 0)
!        Mask_120 = .true.
!      endwhere
!   endif
!   if (Chan_On_133 > 0) then
!      where ( Dqf_133 <= Dqf_Thresh .AND. Dqf_133 >=0)
!        Mask_133 = .true.
!      endwhere
!   endif
!   
!   Count_Good_104 = count(Mask_104)
!   Count_Good_110 = count(Mask_110)
!   Count_Good_120 = count(Mask_120)
!   Count_Good_133 = count(Mask_133)

!   if (Count_Good_104 < Count_Good_110 .OR. &
!       Count_Good_120 < Count_Good_110 .OR. &
!       Count_Good_133 < Count_Good_110) then
!         Use_10_4_Flag = .true.
!   endif

   !Use_10_4_Flag = .true.

!end subroutine SET_10_4_FLAG
!---------------------------------------------------------------------------
! Determine is the chosen mode is consistent with DQFs, if not modify
!----------------------------------------------------------------------------
!subroutine MODIFY_MODE_USING_DQF(Bad_Data_Mask, &
!                                 Chan_On_038, Chan_On_067, Chan_On_085,
!                                 Chan_On_120, Chan_On_133, &
!                                 Dqf_038, Dqf_067, Dqf_085, Dqf_120, Dqf_133,
!                                 Acha_Mode_Flag)
!   character(len=*), intent(inout):: Acha_Mode_Flag
!   integer (kind=int1), intent(in), dimension(:,:) :: Bad_Data_Mask
!   integer, intent(in):: Chan_On_038, Chan_On_067, Chan_On_085, Chan_On_120,
!   Chan_On_133
!   integer (kind=int1), intent(in), dimension(:,:) :: Dqf_038, Dqf_067,
!   Dqf_085, Dqf_120, Dqf_133
!   real:: Fraction_Good_038, Fraction_Good_067, Fraction_Good_085, &
!          Fraction_Good_120, Fraction_Good_133, Count_Total
!   real:: Fraction_Good_Thresh = 0.50
!   logical, dimension(size(Bad_Data_Mask(:,1)),size(Bad_Data_Mask(1,:)))::
!   Mask_xx
!   character(len=len(Acha_Mode_Flag)):: Acha_Mode_Flag_Input
!
!   !--- store value at beginning
!   Acha_Mode_Flag_Input = Acha_Mode_Flag
!
!   !--- Fraction of Good Data for Each Channel
!   Fraction_Good_038 = 0.0
!   Fraction_Good_067 = 0.0
!   Fraction_Good_085 = 0.0
!   Fraction_Good_120 = 0.0
!   Fraction_Good_133 = 0.0
!
!   Mask_xx = .false.
!   where(Bad_Data_Mask == 0_int1)
!         Mask_xx = .true.
!   endwhere
!   Count_Total = count(Mask_xx)
!
!   if (Chan_On_038 == 1_int4) then
!      Mask_xx = .false.
!      where(DQF_038 <= Dqf_Min_Thresh .AND. DQF_038 >= 0)
!         Mask_xx = .true.
!      endwhere 
!      Fraction_Good_038 = count(Mask_xx) / Count_Total
!   endif
!
!   if (Chan_On_067 == 1_int4) then
!      Mask_xx = .false.
!      where(DQF_067 <= Dqf_Min_Thresh .AND. DQF_067 >= 0)
!         Mask_xx = .true.
!      endwhere 
!      Fraction_Good_067 = count(Mask_xx) / Count_Total
!   endif
!
!   if (Chan_On_085 == 1_int4) then
!      Mask_xx = .false.
!      where(DQF_085 <= Dqf_Min_Thresh .AND. DQF_085 >= 0)
!         Mask_xx = .true.
!      endwhere 
!      Fraction_Good_085 = count(Mask_xx) / Count_Total
!   endif
!
!   if (Chan_On_120 == 1_int4) then
!      Mask_xx = .false.
!      where(DQF_120 <= Dqf_Min_Thresh .AND. DQF_120 >= 0)
!         Mask_xx = .true.
!      endwhere 
!      Fraction_Good_120 = count(Mask_xx) / Count_Total
!   endif
!
!   if (Chan_On_133 == 1_int4) then
!      Mask_xx = .false.
!      where(DQF_133 <= Dqf_Min_Thresh .AND. DQF_133 >= 0)
!         Mask_xx = .true.
!      endwhere 
!      Fraction_Good_133 = count(Mask_xx) / Count_Total
!   endif
!
!   !--- logic 
!   if ( (index(Acha_Mode_Flag,'038') > 0 .and. Fraction_Good_038 <
!   Fraction_Good_Thresh) .or. &
!        (index(Acha_Mode_Flag,'067') > 0 .and. Fraction_Good_067 <
!        Fraction_Good_Thresh) .or. &
!        (index(Acha_Mode_Flag,'085') > 0 .and. Fraction_Good_085 <
!        Fraction_Good_Thresh) .or. &
!        (index(Acha_Mode_Flag,'120') > 0 .and. Fraction_Good_120 <
!        Fraction_Good_Thresh) .or. &
!        (index(Acha_Mode_Flag,'133') > 0 .and. Fraction_Good_133 <
!        Fraction_Good_Thresh) ) then
!
!        Acha_Mode_Flag = "unknown"
!
!        if (Acha_Mode_Flag == "unknown" .and. &
!            Fraction_Good_067 > Fraction_Good_Thresh .and. &
!            Fraction_Good_120 > Fraction_Good_Thresh .and. &
!            Fraction_Good_085 > Fraction_Good_Thresh) then
!            Acha_Mode_Flag = "067_085_110_120"
!        endif
!
!        if (Acha_Mode_Flag == "unknown" .and. &
!            Fraction_Good_085 > Fraction_Good_Thresh .and. &
!            Fraction_Good_120 > Fraction_Good_Thresh) then
!            Acha_Mode_Flag = "085_110_120"
!        endif
!
!        if (Acha_Mode_Flag == "unknown" .and. &
!            Fraction_Good_120 > Fraction_Good_Thresh) then
!            Acha_Mode_Flag = "110_120"
!        endif
!
!        if (Acha_Mode_Flag == "unknown" .and. &
!            Fraction_Good_067 > Fraction_Good_Thresh) then
!            Acha_Mode_Flag = "067_110"
!        endif
!
!        if (Acha_Mode_Flag == "unknown" .and. &
!            Fraction_Good_038 > Fraction_Good_Thresh) then
!            Acha_Mode_Flag = "038_110"
!        endif
!
!        if (Acha_Mode_Flag == "unknown") then
!            Acha_Mode_Flag = "110"
!        endif
!
!    endif
!
!    !--- print warning to screen if a mode change occurred
!    if (Acha_Mode_Flag_Input /= Acha_Mode_Flag) then
!       print *, "WARNING: DQFs changed ACHA Mode from
!       ",trim(Acha_Mode_Flag_Input), " to ", trim(Acha_Mode_Flag)
!    endif
!
!end subroutine MODIFY_MODE_USING_DQF



end module ACHA_LHP_MOD
