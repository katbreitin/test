!$Id: ahi_mod.f90 3082 2018-12-17 17:53:19Z mfoster $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: AHI_module.f90 (src)
!       AHI_MODULE (program)
!
! PURPOSE: 
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
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
! AHI Channel Mapping
!
! wvl       ahi  modis/clavrx 
!
! 0.47       1      3   
! 0.51       2      4
! 0.64       3      1  
! 0.86       4      2  
! 1.6        5      6
! 2.2        6     r7
! 3.9        7     20
! 6.2        8     37
! 6.9        9     27
! 7.3       10     28
! 8.6       11     29
! 9.6       12     30
! 10.4      13     38
! 11.2      14     31
! 12.3      15     32 
! 13.3      16     33
!
!--------------------------------------------------------------------------------------
module AHI_MOD

  use CLAVRX_MESSAGE_MOD,only: &
      mesg &
      ,verb_lev
      
  use CALIBRATION_CONSTANTS_MOD, only:  &
                sat_name &
                , solar_ch20 &
                , solar_ch20_nu &
                , ew_ch20 &
                , planck_a1 &
                , planck_a2 &
                , planck_nu  &
                , Rad_to_Ref_Fac_0_47UM &
                , Rad_to_Ref_Fac_0_55um &
                , Rad_to_Ref_Fac_0_65um &
                , Rad_to_Ref_Fac_0_86um &
                , Rad_to_Ref_Fac_1_60um &
                , Rad_to_Ref_Fac_2_10um

  
  use CONSTANTS_MOD, only: &
    exe_prompt
  
  use FILE_UTILS, only: get_lun

  implicit none

  private
  public:: READ_AHI_INSTR_CONSTANTS

  character(len=13), parameter:: MODULE_PROMPT=" AHI_MODULE: "

contains

  !----------------------------------------------------------------
  ! read the AHI constants into memory
  !-----------------------------------------------------------------
  subroutine READ_AHI_INSTR_CONSTANTS(Instr_Const_file)
    character(len=*), intent(in):: Instr_Const_file
    integer:: ios0, erstat
    integer:: Instr_Const_lun
    real,parameter :: Missing_Local = -999.0

    Instr_Const_lun = GET_LUN()

    open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

    print *, EXE_PROMPT, MODULE_PROMPT, " Opening ", trim(Instr_Const_file)
    erstat = 0
    if (ios0 /= 0) then
      erstat = 19
      print *, EXE_PROMPT, MODULE_PROMPT, "Error opening AHI constants file, ios0 = ", ios0
      stop 19
    endif

    read(unit=Instr_Const_lun,fmt="(a4)") sat_name
    read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
    read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
    read(unit=Instr_Const_lun,fmt=*) planck_a1(20), planck_a2(20),planck_nu(20) ! Band 7
    !Note AHI has a 6.2 (Band 8), but MODIS doesn't have one
    read(unit=Instr_Const_lun,fmt=*) planck_a1(27), planck_a2(27),planck_nu(27) !Band 9
    read(unit=Instr_Const_lun,fmt=*) planck_a1(28), planck_a2(28),planck_nu(28) !Band 10
    read(unit=Instr_Const_lun,fmt=*) planck_a1(29), planck_a2(29),planck_nu(29) !Band 11
    read(unit=Instr_Const_lun,fmt=*) planck_a1(30), planck_a2(30),planck_nu(30) !Band 12
    !NOTE AHI as a 10.4 (Band 13), but MODIS doesn't have one
    read(unit=Instr_Const_lun,fmt=*) planck_a1(31), planck_a2(31),planck_nu(31) !Band 14
    read(unit=Instr_Const_lun,fmt=*) planck_a1(32), planck_a2(32),planck_nu(32) !Band 15
    read(unit=Instr_Const_lun,fmt=*) planck_a1(33), planck_a2(33),planck_nu(33) !Band 16
    read(unit=Instr_Const_lun,fmt=*) planck_a1(37), planck_a2(37),planck_nu(37) !Band 8
    read(unit=Instr_Const_lun,fmt=*) planck_a1(38), planck_a2(38),planck_nu(38) !Band 13

    !--- these are not used, the values in the L1b values are used
    read(unit=Instr_Const_lun,fmt=*) Rad_to_Ref_Fac_0_47um !AHI Band 1
    read(unit=Instr_Const_lun,fmt=*) Rad_to_Ref_Fac_0_55um !AHI Band 2
    read(unit=Instr_Const_lun,fmt=*) Rad_to_Ref_Fac_0_65um !AHI Band 3
    read(unit=Instr_Const_lun,fmt=*) Rad_to_Ref_Fac_0_86um !AHI Band 4
    read(unit=Instr_Const_lun,fmt=*) Rad_to_Ref_Fac_1_60um !AHI Band 5
    read(unit=Instr_Const_lun,fmt=*) Rad_to_Ref_Fac_2_10um !AHI Band 6

    close(unit=Instr_Const_lun)

    !--- set Rad_to_Ref as missing since we use the information in the L1b
    !Rad_to_Ref_Fac_0_47um = Missing_Local
    !Rad_to_Ref_Fac_0_55um = Missing_Local
    !Rad_to_Ref_Fac_0_65um = Missing_Local
    !Rad_to_Ref_Fac_0_86um = Missing_Local
    !Rad_to_Ref_Fac_1_60um = Missing_Local
    !Rad_to_Ref_Fac_2_10um = Missing_Local

    !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
    Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

    call mesg ( "AHI instrument constants read in successfully",level = verb_lev % DEFAULT)

  end subroutine READ_AHI_INSTR_CONSTANTS

end module AHI_MOD
