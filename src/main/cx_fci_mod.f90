module cx_fci_mod
use FCI_MOD

use pixel_common_mod, only: ch, image, sensor

contains


subroutine read_fci (seg_nr)
  integer, intent(in) :: seg_nr
  type(fci_data) :: fci
  logical :: fci_on(16) = .true.
  integer :: i
  integer :: idx_cx
print*,trim(image%Level1b_Full_Name)
  call fci % config % set(trim(image%Level1b_Full_Name)//'/',fci_on)
  call fci % get (chunk = seg_nr ) !, start=[10,10],count=[20,20])
do i = 1,16
  idx_cx = Sensor%CLAVRx_Chan_Map(i)
  if ( sensor % Chan_On_Flag_Default( Sensor%CLAVRx_Chan_Map(idx_ch)) .eq. 0 ) cycle
  print*,i,idx_cx

  ch(idx_cx) % rad_toa = fci % ch(i) % rad
  !ch(idx_cx) % rfl_toa = fci % ch(i) % rfl
  !ch(idx_cx) % Bt_Toa = fci % ch(i) % bt

  print*,maxval(ch(idx_cx) % rad_toa)
end do
stop
end subroutine read_fci

end module cx_fci_mod
