!$Id:$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: cx_abi_lhp_mod_alg.f90 (src)
!       
!
! PURPOSE: Given the CLAVRx channel index (Chan_Idx), The focal 
!          plane temperature thresholds (set in a given algorithm, this routine will turn 
!          off a given LOCAL channel flag for a given algorithm. This also would set the 
!          local use_10um flag. Basics taken from CX_ABI_LHP_MOD
!
!          NOTE - Thresholds can be done however is needed
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
!  December 2019 - created
!
!--------------------------------------------------------------------------------------
module CX_ABI_LHP_MOD_ALG_FUNCTION

  use CONSTANTS_MOD

  use PIXEL_COMMON_MOD, only: &
      Ch, &
      Sensor, &
      Image

  use CLAVRX_MESSAGE_MOD, only: &
      MESG, &
      VERB_LEV

   use FILE_UTILS, only: &
      Get_Lun

  implicit none

  public:: LHP_LOCAL_CHAN_ON

contains



  !-------------------------------------------------------------------------------------------------
  !This routine is set so that it can be run in a loop or for just specific channels
  ! in a given algorithm and will do the check on the FPT already in memory.
  ! Could be a function
  !-------------------------------------------------------------------------------------------------

  integer function LHP_LOCAL_CHAN_ON (Chan_Idx, ABI_FPT_Thresh_Algo)
 
   implicit none
   integer, intent(in):: Chan_Idx ! CLAVR-x Channel index
   real(kind=real4), intent(in) :: ABI_FPT_Thresh_Algo
   character(len=100):: header
    
   !Initalize to sensor default 
   LHP_LOCAL_CHAN_ON = Sensor%Chan_On_Flag_Default(Chan_Idx)

   if (Sensor%WMO_Id /= 271) return !This limits to only using GOES-17
                                    ! could be done in routine

   if ((Ch(Chan_Idx)%Max_Focal_Plane_Temp > ABI_FPT_Thresh_Algo) .and. &
        ( Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES)) then
       LHP_LOCAL_CHAN_ON = sym%NO
   endif

  end function LHP_LOCAL_CHAN_ON
  
  
!--------------------------------------------------------------------------------------------------  
!
!--------------------------------------------------------------------------------------------------  
 subroutine READ_LHP_THRESH_FILE (Algo_LHP_file, Algo_LHP_Thresh)
    character(len=*), intent(in):: Algo_LHP_file
    real, dimension(:), intent(out):: Algo_LHP_Thresh
    character(len=200):: header
    integer:: ios0, erstat, io
    integer:: Instr_Const_lun
    integer:: Chan_Idx
    real:: Threshold

   Instr_Const_lun = Get_Lun()
   Algo_LHP_Thresh = MISSING_VALUE_REAL4

   !print *, "opening ", trim(Algo_LHP_file)
   open(unit=Instr_Const_lun,file=trim(Algo_LHP_file), &
        status="old",position="rewind",&
        action="read",iostat=ios0)

    erstat = 0
    if (ios0 /= 0) then
      erstat = 19
      print *, EXE_PROMPT, "Error opening threshold constants file, ios0 = ", ios0
      print *, "filename is "//trim(Algo_LHP_file)
      stop 19
    end if
   read(unit=Instr_Const_lun,fmt=*) header
   read(unit=Instr_Const_lun,fmt=*) header
   do
        read(unit=Instr_Const_lun,fmt=*, iostat = io) Chan_Idx, Threshold
        if (io /= 0) exit
        Algo_LHP_Thresh(Chan_Idx) = Threshold
   enddo  


    close(unit=Instr_Const_lun)

  end subroutine READ_LHP_THRESH_FILE
  
  
!--------------------------------------------------------------------------------------------------  
!
!--------------------------------------------------------------------------------------------------  
end module CX_ABI_LHP_MOD_ALG_FUNCTION
