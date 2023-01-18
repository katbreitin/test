module cx_fci_mod
use FCI_MOD

use pixel_common_mod, only: ch, image, sensor, nav, geo
implicit none
contains


subroutine read_fci (seg_nr)
  integer, intent(in) :: seg_nr
  type(fci_data) :: fci
  logical :: fci_on(16) = .true.
  integer :: i,ii
  integer :: idx_cx
  integer :: c_seg_lines
  integer :: stride
  integer :: ubnd(2)
  integer :: i_line
  integer :: ndims(2)




  call fci % config % set(trim(image%Level1b_Full_Name)//'/',fci_on)
  call fci % get (chunk = seg_nr ) !, start=[10,10],count=[20,20])
  do i = 1,16
    if ( i .eq. 9) cycle
    
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
