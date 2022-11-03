!$Id:$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: cx_abi_lhp_mod.f90 (src)
!       
!
! PURPOSE: Routines associated with GOES-17 loop heat pipe issues.
!
! DESCRIPTION: see below
!
! AUTHORS:
!  
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! REVISION HISTORY:
!  June 2019 - created
!
!--------------------------------------------------------------------------------------
module CX_ABI_LHP_MOD

  use CALIBRATION_CONSTANTS_MOD,only: &
      ABI_FPT_Thresh_038um, ABI_FPT_Thresh_062um, ABI_FPT_Thresh_067um &
    , ABI_FPT_Thresh_073um, ABI_FPT_Thresh_085um, ABI_FPT_Thresh_097um &
    , ABI_FPT_Thresh_104um, ABI_FPT_Thresh_110um, ABI_FPT_Thresh_120um &
    , ABI_FPT_Thresh_133um

  use CONSTANTS_MOD

  use PIXEL_COMMON_MOD, only: &
      Ch, &
      Sensor, &
      Image

  use CLAVRX_MESSAGE_MOD, only: &
      MESG, &
      VERB_LEV

  implicit none

  public:: SET_ABI_USE_104um_FLAG

contains

  subroutine SET_ABI_USE_104um_FLAG(ABI_Use_104um_Flag, File_Number)

    implicit none
    logical, intent(out) :: ABI_Use_104um_Flag
    integer, intent(in) :: File_Number
    real(kind=real4) :: MFPT_104
    real(kind=real4) :: MFPT_110
    
    

    !--- Initialize.
    ABI_Use_104um_Flag = .false. 
    MFPT_104 = MISSING_VALUE_REAL4
    MFPT_110 = MISSING_VALUE_REAL4


    IF (Sensor%WMO_Id /=271) return
    IF (Sensor%Chan_On_Flag_Default(38) == 0) return

    !--- See what L1b file number we are on. First L1b file is tagged with a 2.
    !--- Save default channel_on state for use if there are more than 1 L1b
    !--- files in file_list.

    if (File_Number == 2) then
      Sensor%Chan_On_Flag_Initial_G17_Only = Sensor%Chan_On_Flag_Default
      Sensor%Chan_On_Flag_Initial_G17_Only_Previous = Sensor%Chan_On_Flag_Default
    endif

    if (File_Number > 2) then
      Sensor%Chan_On_Flag_Default = Sensor%Chan_On_Flag_Initial_G17_Only
    endif

    !--- Set local variables for easier reading.
    MFPT_104 = Ch(38)%Max_Focal_Plane_Temp
    MFPT_110 = Ch(31)%Max_Focal_Plane_Temp

    !--- Re-set the Channel on flags if maximum focal plane temperatures
    !--- exceed the thresholds.
    
    if (Ch(20)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_038um) then
       call MESG ("Focal plane temperature issue, turnning Channel 20 off", &
                   level = verb_lev % DEFAULT)
       Sensor%Chan_On_Flag_Default(20) = 0
    endif

    if (Ch(37)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_062um) then
       call MESG ("Focal plane temperature issue, turnning Channel 37 off", &
                   level = verb_lev % DEFAULT)
       Sensor%Chan_On_Flag_Default(37) = 0
    endif

    if (Ch(27)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_067um) then
       call MESG ("Focal plane temperature issue, turnning Channel 27 off", &
                   level = verb_lev % DEFAULT)
       Sensor%Chan_On_Flag_Default(27) = 0
    endif

    if (Ch(28)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_073um) then
       call MESG ("Focal plane temperature issue, turnning Channel 28 off", &
                   level = verb_lev % DEFAULT)
       Sensor%Chan_On_Flag_Default(28) = 0
    endif

    if (Ch(29)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_085um) then
       call MESG ("Focal plane temperature issue, turnning Channel 29 off", &
                   level = verb_lev % DEFAULT)
       Sensor%Chan_On_Flag_Default(29) = 0
    endif

    if (Ch(30)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_097um) then
       call MESG ("Focal plane temperature issue, turnning Channel 30 off", &
                   level = verb_lev % DEFAULT)
       Sensor%Chan_On_Flag_Default(30) = 0
    endif

    if (Ch(38)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_104um) then
       call MESG ("Focal plane temperature issue, turnning Channel 38 off", &
                   level = verb_lev % DEFAULT)
       Sensor%Chan_On_Flag_Default(38) = 0
    endif

    if (Ch(31)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_110um) then
       call MESG ("Focal plane temperature issue, turnning Channel 31 off", &
                   level = verb_lev % DEFAULT)
       Sensor%Chan_On_Flag_Default(31) = 0
    endif

    if (Ch(32)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_120um) then
       call MESG ("Focal plane temperature issue, turnning Channel 32 off", &
                   level = verb_lev % DEFAULT)
       Sensor%Chan_On_Flag_Default(32) = 0
    endif

    if (Ch(33)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_133um) then
       call MESG ("Focal plane temperature issue, turnning Channel 33 off", &
                   level = verb_lev % DEFAULT)
       Sensor%Chan_On_Flag_Default(33) = 0
    endif

    !--- Needed to fix planck table issues when channel mapping has changed
    !--- because of LHP event.
    if (SUM(Sensor%Chan_On_Flag_Initial_G17_Only) .NE.  SUM(Sensor%Chan_On_Flag_Default)) then
      Sensor%Chan_On_Flag_Initial_G17_Only_Previous = Sensor%Chan_On_Flag_Default
    endif

    !--- Only change the 10.4 use flag if 10.4 is within threshold, but 11 um
    !--- isn't.

    if (MFPT_104 < ABI_FPT_Thresh_104um .AND. MFPT_110 > ABI_FPT_Thresh_110um) then
      call MESG ("Setting 10.4um as the window channel.", level = verb_lev % DEFAULT)
      ABI_Use_104um_Flag = .true.
    else if (MFPT_104 > ABI_FPT_Thresh_104um .AND. MFPT_110 > ABI_FPT_Thresh_110um) then
      call MESG( "ERROR: Catastrophic LHP Event, stopping" , level = verb_lev % DEFAULT)
      stop
    else
      call MESG ("Leaving 11.0um as the window channel.", level = verb_lev % DEFAULT)
    endif

  end subroutine SET_ABI_USE_104um_FLAG

end module CX_ABI_LHP_MOD
