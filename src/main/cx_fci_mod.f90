module cx_fci_mod
  use FCI_MOD,only: fci_data

  use pixel_common_mod, only: ch, image, sensor, nav, geo, Ancil_Data_Dir

  use CALIBRATION_CONSTANTS_MOD,only: &
  planck_a1, planck_a2, planck_nu &
  , sat_name, solar_ch20, ew_ch20, solar_ch20_nu, ch1_dark_count &
  , Launch_Date, Band2_Correction_Start_Date, Band2_Correction_Factor &
  , ABI_FPT_Thresh_038um, ABI_FPT_Thresh_062um, ABI_FPT_Thresh_067um &
  , ABI_FPT_Thresh_073um, ABI_FPT_Thresh_085um, ABI_FPT_Thresh_097um &
  , ABI_FPT_Thresh_104um, ABI_FPT_Thresh_110um, ABI_FPT_Thresh_120um &
  , ABI_FPT_Thresh_133um

  implicit none


  logical :: calib_is_read = .false.
contains

  subroutine read_fci_calib
    use FILE_UTILS,only: Get_Lun
    character (len=1024) :: coef_file
    integer :: Instr_Const_lun
    character(len=20) :: dum

    coef_file = trim(Ancil_data_Dir)//'static/clavrx_constant_files/met12_instr.dat'

    Instr_Const_lun = GET_LUN()
    open(unit=Instr_Const_lun,file=trim(coef_file))

    read(unit=Instr_Const_lun,fmt="(a7)") sat_name
    read(unit=Instr_Const_lun,fmt="(a20)") dum
    read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
    read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
    read(unit=Instr_Const_lun,fmt="(a20)") dum
    read(unit=Instr_Const_lun,fmt=*) planck_a1(20), planck_a2(20), planck_nu(20)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(37), planck_a2(37), planck_nu(37)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(28), planck_a2(28), planck_nu(28)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(29), planck_a2(29), planck_nu(29)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(30), planck_a2(30), planck_nu(30)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(31), planck_a2(31), planck_nu(31)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(32), planck_a2(32), planck_nu(32)
    read(unit=Instr_Const_lun,fmt=*) planck_a1(33), planck_a2(33), planck_nu(33)
    calib_is_read = .true.

  end subroutine


  subroutine read_fci (seg_nr)
    integer, intent(in) :: seg_nr
    type(fci_data) :: fci
    logical :: fci_on(16) = .true.
    integer :: i
    integer :: idx_cx
    integer :: c_seg_lines
    integer :: stride
    integer :: ubnd(2)
    integer :: i_line
    integer :: ndims(2)


    if (.not. calib_is_read) call read_fci_calib()

    do i = 1,16
      fci_on(i) = Sensor%Chan_On_Flag_Default ( Sensor%CLAVRx_Chan_Map(i)) == 1
    end do

    call fci % config % set(trim(image%Level1b_Full_Name)//'/',fci_on)
    call fci % get (chunk = seg_nr ) !, start=[10,10],count=[20,20])
    do i = 1,16


      stride = 1
      ubnd = ubound(fci % ch(i) % rad)
      if (i .le. 8 ) stride = 2
      ubnd(2) = 139
      if (stride .eq. 2) THEN
        ubnd(2) = 278
      end if

      idx_cx = Sensor%CLAVRx_Chan_Map(i)

      if ( sensor % Chan_On_Flag_Default( idx_cx) .eq. 0 ) cycle


      if ( allocated (ch(idx_cx) % rad_toa ) ) then
        ndims = shape(ch(idx_cx) % rad_toa)
        ch(idx_cx) % rad_toa = &
        fci % ch(i) % rad (1:ubnd(1):stride,1:ubnd(2):stride)
      end if
      if ( allocated (ch(idx_cx) % ref_toa) ) then
        ndims = shape(ch(idx_cx) % ref_toa)
        ch(idx_cx) % ref_toa = &
        fci % ch(i) % rfl(1:ubnd(1):stride,1:ubnd(2):stride)

      end if
      if ( allocated (ch(idx_cx) % Bt_Toa)) THEN
        ndims = shape(ch(idx_cx) % bt_toa)
        ch(idx_cx) % Bt_Toa = &
        fci % ch(i) % bt(1:ubnd(1):stride,1:ubnd(2):stride)
      end if

    end do
    c_seg_lines = Image%Number_Of_Lines_Per_Segment
    ! transfer all geo data
    nav % lat_1b(:,1:c_seg_lines)      = fci % geo % lat
    nav % lon_1b(:,1:c_seg_lines)       =  fci  % geo % lon

    geo % sataz(:,1:c_seg_lines)        =  fci  % geo % sataz
    geo % satzen(:,1:c_seg_lines)       =  fci  % geo % satzen
    geo % solaz (:,1:c_seg_lines)       =  fci  % geo % solaz
    geo % solzen (:,1:c_seg_lines)      =  fci  % geo % solzen
    geo % relaz (:,1:c_seg_lines)       =  fci  % geo % relaz
    geo % glintzen (:,1:c_seg_lines)    =  fci  % geo % glintzen
    geo % scatangle (:,1:c_seg_lines)   =  fci  % geo % scatangle

    Image%Number_Of_Lines_Read_This_Segment = c_seg_lines
    Image%Scan_Number = [(i_line , i_line = (seg_nr -1) *139 + 1 &
    , (seg_nr -1) *139 + 139   , 1)]

  end subroutine read_fci

end module cx_fci_mod
