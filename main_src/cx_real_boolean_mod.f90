! $Id:$
!----------------------------------------------------------------------
!   This tools creates boolean operators to compare REAL numbers
!
!   Author: Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!
!   New Operators are
!   .EQR. replaces .EQ. or ==
!   .NER. replaces .NE. or /=
!   .LTR. replaces .LT. or <
!   .LER. replaces .LE. or <=
!   .GTR. replaces .GT. or > 
!   .GER. replaces .GE. or >= 
!
!   operators are not case sensitive
!
!   Modification History
!   May 2019: First Version
!
!----------------------------------------------------------------------
!   Usage Example
!
!      if (a .EQR. b) then
!
!
!----------------------------------------------------------------------
module CX_REAL_BOOLEAN_MOD
  
! use :: iso_fortran_env
! use, intrinsic :: iso_fortran_env
  
  implicit none

  !--- these are taken from iso_fortran_env
  INTEGER (KIND=4), PARAMETER                :: REAL64 = 8
  INTEGER (KIND=4), PARAMETER                :: REAL128 = 16
  

  interface operator (.EQR.)
    module procedure eq_real
    module procedure eq_real64
    module procedure eq_real128
    module procedure eq_real_1d
    module procedure eq_real_2d
    module procedure eq_real64_2d
    module procedure eq_real128_2d
    module procedure eq_real_3d
  end interface

  interface operator (.NER.)
    module procedure ne_real
    module procedure ne_real64
    module procedure ne_real128
    module procedure ne_real_1d
    module procedure ne_real_2d
    module procedure ne_real64_2d
    module procedure ne_real128_2d
    module procedure ne_real_3d
  end interface
  
  interface operator (.LTR.)
    module procedure lt_real
    module procedure lt_real64
    module procedure lt_real128
    module procedure lt_real_1d
    module procedure lt_real_2d
    module procedure lt_real64_2d
    module procedure lt_real128_2d
    module procedure lt_real_3d
  end interface

  interface operator (.LER.)
    module procedure le_real
    module procedure le_real64
    module procedure le_real128
    module procedure le_real_1d
    module procedure le_real_2d
    module procedure le_real64_2d
    module procedure le_real128_2d
    module procedure le_real_3d
  end interface
  
  interface operator (.GTR.)
    module procedure gt_real
    module procedure gt_real64
    module procedure gt_real128
    module procedure gt_real_1d
    module procedure gt_real_2d
    module procedure gt_real64_2d
    module procedure gt_real128_2d
    module procedure gt_real_3d
  end interface

  interface operator (.GER.)
    module procedure ge_real
    module procedure ge_real64
    module procedure ge_real128
    module procedure ge_real_1d
    module procedure ge_real_2d
    module procedure ge_real64_2d
    module procedure ge_real128_2d
    module procedure ge_real_3d
  end interface
  
contains

  !-------------------------
  !-- .EQR. functions
  !-------------------------
  logical function eq_real(ra,rb)
    real, intent(in) :: ra, rb 
    eq_real = abs(ra - rb) < epsilon(max(ra,rb))
  end function eq_real
 
  logical function eq_real64(ra,rb)
    real(kind=real64), intent(in) :: ra, rb 
    eq_real64 = abs(ra - rb) < epsilon(max(ra,rb))
  end function eq_real64
  
  logical function eq_real128(ra,rb)
    real(kind=real128), intent(in) :: ra, rb 
    eq_real128 = abs(ra - rb) < epsilon(max(ra,rb))
  end function eq_real128  

  function eq_real_1d(ra,rb) result (out)
    real, intent(in) :: ra(:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:)
    integer :: dim_x(1)
    dim_x = shape(ra)
    allocate(out(dim_x(1)))
    out  = abs(ra - rb) < epsilon(rb)
  end function eq_real_1d

  function eq_real_2d(ra,rb) result (out)
    real, intent(in) :: ra(:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = abs(ra - rb) < epsilon(rb)
  end function eq_real_2d
  
  function eq_real64_2d(ra,rb) result (out)
    real(kind=real64), intent(in) :: ra(:,:)
    real(kind=real64), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = abs(ra - rb) < epsilon(rb)
  end function eq_real64_2d
   
  function eq_real128_2d(ra,rb) result (out)
    real(kind=real128), intent(in) :: ra(:,:)
    real(kind=real128), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = abs(ra - rb) < epsilon(rb)
  end function eq_real128_2d

  function eq_real_3d(ra,rb) result (out)
    real, intent(in) :: ra(:,:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:,:)
    integer :: dim_x(3)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2),dim_x(3)))
    out  = abs(ra - rb) < epsilon(rb)
  end function eq_real_3d

  !-------------------------
  !-- .NER. functions
  !-------------------------
  logical function ne_real(ra,rb)
    real, intent(in) :: ra, rb 
    ne_real = abs(ra - rb) > epsilon(ra)
  end function ne_real
 
  logical function ne_real64(ra,rb)
    real(kind=real64), intent(in) :: ra, rb 
    ne_real64 = abs(ra - rb) > epsilon(ra)
  end function ne_real64
  
  logical function ne_real128(ra,rb)
    real(kind=real128), intent(in) :: ra, rb 
    ne_real128 = abs(ra - rb) > epsilon(ra)
  end function ne_real128  

  function ne_real_1d(ra,rb) result (out)
    real, intent(in) :: ra(:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:)
    integer :: dim_x(1)
    dim_x = shape(ra)
    allocate(out(dim_x(1)))
    out  = abs(ra - rb) > epsilon(rb)
  end function ne_real_1d

  function ne_real_2d(ra,rb) result (out)
    real, intent(in) :: ra(:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = abs(ra - rb) > epsilon(rb)
  end function ne_real_2d
  
  function ne_real64_2d(ra,rb) result (out)
    real(kind=real64), intent(in) :: ra(:,:)
    real(kind=real64), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = abs(ra - rb) > epsilon(rb)
  end function ne_real64_2d
   
  function ne_real128_2d(ra,rb) result (out)
    real(kind=real128), intent(in) :: ra(:,:)
    real(kind=real128), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = abs(ra - rb) > epsilon(rb)
  end function ne_real128_2d

  function ne_real_3d(ra,rb) result (out)
    real, intent(in) :: ra(:,:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:,:)
    integer :: dim_x(3)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2),dim_x(3)))
    out  = abs(ra - rb) > epsilon(rb)
  end function ne_real_3d
  
  !-------------------------
  !-- .LTR. functions
  !-------------------------
  logical function lt_real(ra,rb)
    real, intent(in) :: ra, rb 
    lt_real = (ra - rb) < epsilon(rb)
  end function lt_real
 
  logical function lt_real64(ra,rb)
    real(kind=real64), intent(in) :: ra, rb 
    lt_real64 = (ra - rb) < epsilon(rb)
  end function lt_real64
  
  logical function lt_real128(ra,rb)
    real(kind=real128), intent(in) :: ra, rb 
    lt_real128 = (ra - rb) < epsilon(rb)
  end function lt_real128  

  function lt_real_1d(ra,rb) result (out)
    real, intent(in) :: ra(:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:)
    integer :: dim_x(1)
    dim_x = shape(ra)
    allocate(out(dim_x(1)))
    out  = (ra - rb) < epsilon(rb)
  end function lt_real_1d

  function lt_real_2d(ra,rb) result (out)
    real, intent(in) :: ra(:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) < epsilon(rb)
  end function lt_real_2d
  
  function lt_real64_2d(ra,rb) result (out)
    real(kind=real64), intent(in) :: ra(:,:)
    real(kind=real64), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) < epsilon(rb)
  end function lt_real64_2d
   
  function lt_real128_2d(ra,rb) result (out)
    real(kind=real128), intent(in) :: ra(:,:)
    real(kind=real128), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) < epsilon(rb)
  end function lt_real128_2d

  function lt_real_3d(ra,rb) result (out)
    real, intent(in) :: ra(:,:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:,:)
    integer :: dim_x(3)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2),dim_x(3)))
    out  = (ra - rb) < epsilon(rb)
  end function lt_real_3d

  !-------------------------
  !-- .LER. functions
  !-------------------------
  logical function le_real(ra,rb)
    real, intent(in) :: ra, rb 
    le_real = (ra - rb) < epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function le_real
 
  logical function le_real64(ra,rb)
    real(kind=real64), intent(in) :: ra, rb 
    le_real64 = (ra - rb) < epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function le_real64
  
  logical function le_real128(ra,rb)
    real(kind=real128), intent(in) :: ra, rb 
    le_real128 = (ra - rb) < epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function le_real128  

  function le_real_1d(ra,rb) result (out)
    real, intent(in) :: ra(:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:)
    integer :: dim_x(1)
    dim_x = shape(ra)
    allocate(out(dim_x(1)))
    out  = (ra - rb) < epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function le_real_1d

  function le_real_2d(ra,rb) result (out)
    real, intent(in) :: ra(:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) < epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function le_real_2d
  
  function le_real64_2d(ra,rb) result (out)
    real(kind=real64), intent(in) :: ra(:,:)
    real(kind=real64), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) < epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function le_real64_2d
   
  function le_real128_2d(ra,rb) result (out)
    real(kind=real128), intent(in) :: ra(:,:)
    real(kind=real128), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) < epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function le_real128_2d

  function le_real_3d(ra,rb) result (out)
    real, intent(in) :: ra(:,:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:,:)
    integer :: dim_x(3)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2),dim_x(3)))
    out  = (ra - rb) < epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function le_real_3d
  
  !-------------------------
  !-- .GTR. functions
  !-------------------------
  logical function gt_real(ra,rb)
    real, intent(in) :: ra, rb 
    gt_real = (ra - rb) > epsilon(rb)
  end function gt_real
 
  logical function gt_real64(ra,rb)
    real(kind=real64), intent(in) :: ra, rb 
    gt_real64 = (ra - rb) > epsilon(rb)
  end function gt_real64
  
  logical function gt_real128(ra,rb)
    real(kind=real128), intent(in) :: ra, rb 
    gt_real128 = (ra - rb) > epsilon(rb)
  end function gt_real128  

  function gt_real_1d(ra,rb) result (out)
    real, intent(in) :: ra(:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:)
    integer :: dim_x(1)
    dim_x = shape(ra)
    allocate(out(dim_x(1)))
    out  = (ra - rb) > epsilon(rb)
  end function gt_real_1d

  function gt_real_2d(ra,rb) result (out)
    real, intent(in) :: ra(:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) > epsilon(rb)
  end function gt_real_2d
  
  function gt_real64_2d(ra,rb) result (out)
    real(kind=real64), intent(in) :: ra(:,:)
    real(kind=real64), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) > epsilon(rb)
  end function gt_real64_2d
   
  function gt_real128_2d(ra,rb) result (out)
    real(kind=real128), intent(in) :: ra(:,:)
    real(kind=real128), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) > epsilon(rb)
  end function gt_real128_2d

  function gt_real_3d(ra,rb) result (out)
    real, intent(in) :: ra(:,:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:,:)
    integer :: dim_x(3)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2),dim_x(3)))
    out  = (ra - rb) > epsilon(rb)
  end function gt_real_3d

  !-------------------------
  !-- .GER. functions
  !-------------------------
  logical function ge_real(ra,rb)
    real, intent(in) :: ra, rb 
    ge_real = (ra - rb) > epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function ge_real
 
  logical function ge_real64(ra,rb)
    real(kind=real64), intent(in) :: ra, rb 
    ge_real64 = (ra - rb) > epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function ge_real64
  
  logical function ge_real128(ra,rb)
    real(kind=real128), intent(in) :: ra, rb 
    ge_real128 = (ra - rb) > epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function ge_real128  

  function ge_real_1d(ra,rb) result (out)
    real, intent(in) :: ra(:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:)
    integer :: dim_x(1)
    dim_x = shape(ra)
    allocate(out(dim_x(1)))
    out  = (ra - rb) > epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function ge_real_1d

  function ge_real_2d(ra,rb) result (out)
    real, intent(in) :: ra(:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) > epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function ge_real_2d
  
  function ge_real64_2d(ra,rb) result (out)
    real(kind=real64), intent(in) :: ra(:,:)
    real(kind=real64), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) > epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function ge_real64_2d
   
  function ge_real128_2d(ra,rb) result (out)
    real(kind=real128), intent(in) :: ra(:,:)
    real(kind=real128), intent(in) :: rb 
    logical, allocatable :: out (:,:)
    integer :: dim_x(2)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2)))
    out  = (ra - rb) > epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function ge_real128_2d

  function ge_real_3d(ra,rb) result (out)
    real, intent(in) :: ra(:,:,:)
    real, intent(in) :: rb 
    logical, allocatable :: out (:,:,:)
    integer :: dim_x(3)
    dim_x = shape(ra)
    allocate(out(dim_x(1),dim_x(2),dim_x(3)))
    out  = (ra - rb) > epsilon(rb) .or. abs(ra - rb) < epsilon(rb)
  end function ge_real_3d
  
end module CX_REAL_BOOLEAN_MOD
