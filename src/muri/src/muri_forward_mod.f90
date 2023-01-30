! $Id: muri_forward_mod.f90 2187 2017-04-14 16:20:00Z awalther $

module muri_forward_mod

  use muri_lut_mod, only: &
  muri_lut_type

  type(muri_lut_type),save :: lut

  type muri_fwd_type
    real :: rfl(6)
    real :: jacobians
  end type muri_fwd_type


contains
  subroutine muri_forward ( state, fwd)
    real, intent(in) :: state
    type ( muri_fwd_type ), intent(out) :: fwd
    !real :: sol, sat,azi


    print*,state
    call lut%read_lut ( )


    !call lut % make_case_lut ( sol,  sat, azi)
    stop


    fwd % rfl(1) = 2.1

  end subroutine
  !
  !
  !




end module  muri_forward_mod
