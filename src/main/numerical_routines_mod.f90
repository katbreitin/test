! $Id: numerical_routines_mod.f90 3945 2020-08-19 20:20:17Z yli $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: num_mod.f90 (src)
!       NUMERICAL_ROUTINES (program)
!
! PURPOSE: library of useful numerical functions
!
! Description: 
!              Note, routines that only appear in the volcanic
!              ash functions are in ash_num_mod.f90.  This was
!              done to minimize code delivery to NCDC for PATMOS-x.
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
! public:: 
!   LOCATE
!   JULIAN
!   COMPUTE_MONTH
!   COMPUTE_DAY
!   VAPOR
!   VAPOR_ICE
!   INVERT_2x2
!   INVERT_3x3
!   INVERT_4x4
!   INVERT_5x5
!   INVERT_6x6
!   FIND_BOUNDS
!   PACK_BYTES
!   COMPUTE_TIME_HOURS
!   COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES
!   GRADIENT_MEANDER
!   COMPUTE_MEDIAN
!   COMPUTE_MEDIAN_SEGMENT
!
!--------------------------------------------------------------------------------------
 module NUMERICAL_ROUTINES_MOD
   
  use CONSTANTS_MOD
  
  implicit none
  private:: LUBKSB, LUDCMP, LU_INVERT, SVDCMP
  public:: LOCATE, &
           INVERT_MATRIX,  &
           INVERT_2x2,  &
           INVERT_3x3,  &
           INVERT_4x4,  &
           INVERT_5x5,  &
           INVERT_6x6,  &
           INVERT_DIAGONAL,  &
           FIND_BOUNDS,  &
           GEN_INVERSE_SVD, &
           GEN_INVERSE_LU, &
           PACK_BYTES_I1, &
           PACK_BYTES_I2
   public :: covariance                
  contains

!-------------------------------------------------------------------------
! subroutine LOCATE(xx, n, x, j)
! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
!--------------------------------------------------------------------------
  subroutine LOCATE(xx, n, x, j)

!   Arguments
    integer,                        intent(in)  :: n
    integer,                        intent(out) :: j
    real (kind=ipre),               intent(in)  :: x
    real (kind=ipre), dimension(:), intent(in)  :: xx

!   Local variables
    integer :: i, jl, jm, ju

    jl = 0
    ju = n + 1
    do i = 1, 2*n
       if (ju-jl <= 1) then
          exit
       endif
       jm = (ju + jl) / 2
       if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if (x == xx(1)) then
       j=1
    else if (x == xx(n)) then
       j = n - 1
    else
       j = jl
    endif

  end subroutine LOCATE

!--------------------------------------------------------------------------
! Invert a square matrix
!
! follows:
! http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
!--------------------------------------------------------------------------
function INVERT_MATRIX(Matrix, Matrix_Inv, Matrix_Size) RESULT(Status)

  real(kind=real4), dimension(:,:), intent(in) :: Matrix
  real(kind=real4), dimension(:,:), intent(OUT) :: Matrix_Inv
  integer(kind=INT4), intent(in) :: Matrix_Size

  integer(kind=INT4) :: Status
  integer(kind=INT4) :: Singular_Flag
  integer:: i,j,ni, nj
  real(kind=real4) :: Zero
  logical:: Diag_Flag

  Status = Sym%SUCCESS

  Zero = epsilon(Matrix(1,1))

  !---- check for a square matrix
  ni = size(Matrix,1)
  nj = size(Matrix,2)
  if (ni /= nj) then
     print *, "size mismatch in invert matrix"
     Status = sym%FAILURE
     return
  endif
  !---- check for a diagonal matrix
  Diag_Flag = .true.
  do i = 1, ni
    do j = 1, nj
       if (i == j) cycle
       if (abs(Matrix(i,j)) > Zero) then
         Diag_Flag = .false.
         exit
       endif
    enddo
  enddo
 
  if (Diag_Flag) then
      call INVERT_DIAGONAL(Matrix, Matrix_Inv, Singular_Flag)
      if (Singular_Flag == 1) THEN
        print *, "diagonal failure"
        Status = Sym%FAILURE
      endif

  else

  !---- if not diagonal matrix, do full inversion
  select case(Matrix_Size)

    case(2)

      call INVERT_2x2(Matrix, Matrix_Inv, Singular_Flag)
      if (Singular_Flag == 1) THEN
        print *, "non-diagonal 2x2 failure"
        Status = Sym%FAILURE
      endif

    case(3)

      call INVERT_3x3(Matrix, Matrix_Inv, Singular_Flag)
      if (Singular_Flag == 1) THEN
        print *, "non-diagonal 3x3 failure"
        Status = Sym%FAILURE
      endif

    case(4)

      call INVERT_4x4(Matrix, Matrix_Inv, Singular_Flag)
      if (Singular_Flag == 1) THEN
        call GEN_INVERSE_LU(Matrix,Matrix_Inv,Matrix_Size,Matrix_Size,Singular_Flag)
!        call GEN_INVERSE_SVD(Matrix,Matrix_Inv,Matrix_Size,Matrix_Size,Singular_Flag)
        if (Singular_Flag == 1) THEN
           print *, "non-diagonal 4x4 failure"
           Status = Sym%FAILURE
        endif
      endif

    case(5)

      call INVERT_5x5(Matrix, Matrix_Inv, Singular_Flag)
      if (Singular_Flag == 1) THEN
        call GEN_INVERSE_LU(Matrix,Matrix_Inv,Matrix_Size,Matrix_Size,Singular_Flag)
!        call GEN_INVERSE_SVD(Matrix,Matrix_Inv,Matrix_Size,Matrix_Size,Singular_Flag)
        if (Singular_Flag == 1) THEN
           print *, "non-diagonal 5x5 failure"
           Status = Sym%FAILURE
        endif
      endif

    case(6)

      call INVERT_6x6(Matrix, Matrix_Inv, Singular_Flag)
      if (Singular_Flag == 1) THEN
        call GEN_INVERSE_LU(Matrix,Matrix_Inv,Matrix_Size,Matrix_Size,Singular_Flag)
!        call GEN_INVERSE_SVD(Matrix,Matrix_Inv,Matrix_Size,Matrix_Size,Singular_Flag)
        if (Singular_Flag == 1) THEN
           print *, "non-diagonal 6x6 failure"
           Status = Sym%FAILURE
        endif
      endif

    case default

      call GEN_INVERSE_LU(Matrix,Matrix_Inv,Matrix_Size,Matrix_Size,Singular_Flag)
!      call GEN_INVERSE_SVD(Matrix,Matrix_Inv,Matrix_Size,Matrix_Size,Singular_Flag)
      if (Singular_Flag == 1) THEN
        print *, "non-diagonal failure, matrix size greater than 4"
        Status = Sym%FAILURE
      endif

  end select

  endif

  return

end function Invert_Matrix

!--------------------------------------------------------------------------
! subroutine INVERT_2x2(A,A_inv,ierr)
!
! Matrix Inversion for a 2x2 matrix
!
! A - input matrix
! A_inv - the inverse of A (output)
! ierr - output flag that report a singular matrix and a failure
!--------------------------------------------------------------------------
subroutine INVERT_2x2(A,A_inv,ierr)
  real, dimension(:,:), intent(in):: A
  real, dimension(:,:), intent(out):: A_inv
  real:: determinant
  integer, intent(out):: ierr

!--- compute determinant
  ierr = 0
  determinant = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  if (determinant == 0.0) then
!       print *, "Singular Matrix in Invert 2x2"
        ierr = 1
  endif

!--- compute inverse
  A_inv(1,1) = A(2,2)
  A_inv(1,2) = -A(1,2)
  A_inv(2,1) = -A(2,1)
  A_inv(2,2) = A(1,1)
  A_inv = A_inv / determinant

end subroutine INVERT_2x2
!--------------------------------------------------------------------------
! subroutine INVERT_3x3(A,A_inv,ierr)
!
! Matrix Inversion for a 3x3 matrix
!--------------------------------------------------------------------------
subroutine INVERT_3x3(AA,AA_inv,ierr)
  real, dimension(:,:), intent(in):: AA
  real, dimension(:,:), intent(out):: AA_inv
  integer, intent(out):: ierr
  real(kind=real8) :: determinant
  real(kind=real8), dimension(3,3):: A, A_inv

  ierr = 0

  A = real(AA,kind=real8)

  !--- compute determinant
  determinant = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3)) - &
                A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3)) + &
                A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
  if (determinant == 0.0) then
        ierr = 1
  endif

  !--- compute inverse
  A_inv(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
  A_inv(1,2) = A(1,3)*A(3,2) - A(3,3)*A(1,2)
  A_inv(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)
  A_inv(2,1) = A(2,3)*A(3,1) - A(3,3)*A(2,1)
  A_inv(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
  A_inv(2,3) = A(1,3)*A(2,1) - A(2,3)*A(1,1)
  A_inv(3,1) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
  A_inv(3,2) = A(1,2)*A(3,1) - A(3,2)*A(1,1)
  A_inv(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

  A_inv = A_inv / determinant

  AA_inv = real(A_inv, kind=real4)

end subroutine INVERT_3x3
!--------------------------------------------------------------------------
! subroutine INVERT_4x4(A,A_inv,ierr)
!
! Matrix Inversion for a 4x4 matrix
!--------------------------------------------------------------------------
subroutine INVERT_4x4(AA,AA_inv,ierr)
  real, dimension(:,:), intent(in):: AA
  real, dimension(:,:), intent(out):: AA_inv
  integer, intent(out):: ierr
  real(kind=real8) :: determinant
  real(kind=real8), dimension(4,4):: A, A_inv

  ierr = 0

  A = real(AA,kind=real8)

  !--- compute determinant
  determinant = A(1,1)*(A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3)) + &
                A(1,2)*(A(2,1)*A(3,4)*A(4,3) + A(2,3)*A(3,1)*A(4,4) + A(2,4)*A(3,3)*A(4,1)) + &
                A(1,3)*(A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2)) + &
                A(1,4)*(A(2,1)*A(3,3)*A(4,2) + A(2,2)*A(3,1)*A(4,3) + A(2,3)*A(3,2)*A(4,1)) - &
                A(1,1)*(A(2,2)*A(3,4)*A(4,3) + A(2,3)*A(3,2)*A(4,4) + A(2,4)*A(3,3)*A(4,2)) - &
                A(1,2)*(A(2,1)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,3)) - &
                A(1,3)*(A(2,1)*A(3,4)*A(4,2) + A(2,2)*A(3,1)*A(4,4) + A(2,4)*A(3,2)*A(4,1)) - &
                A(1,4)*(A(2,1)*A(3,2)*A(4,3) + A(2,2)*A(3,3)*A(4,1) + A(2,3)*A(3,1)*A(4,2))

  if (determinant == 0.0) then
        !print *, "Singular Matrix in Invert 4x4 ", determinant
        ierr = 1
  endif

  !--- compute inverse
  A_inv(1,1) = A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) - &
               A(2,2)*A(3,4)*A(4,3) - A(2,3)*A(3,2)*A(4,4) - A(2,4)*A(3,3)*A(4,2)
  A_inv(1,2) = A(1,2)*A(3,4)*A(4,3) + A(1,3)*A(3,2)*A(4,4) + A(1,4)*A(3,3)*A(4,2) - &
               A(1,2)*A(3,3)*A(4,4) - A(1,3)*A(3,4)*A(4,2) - A(1,4)*A(3,2)*A(4,3)
  A_inv(1,3) = A(1,2)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,2) + A(1,4)*A(2,2)*A(4,3) - &
               A(1,2)*A(2,4)*A(4,3) - A(1,3)*A(2,2)*A(4,4) - A(1,4)*A(2,3)*A(4,2)
  A_inv(1,4) = A(1,2)*A(2,4)*A(3,3) + A(1,3)*A(2,3)*A(3,4) + A(1,4)*A(2,3)*A(3,2) - &
               A(1,2)*A(2,3)*A(3,4) - A(1,3)*A(2,4)*A(3,2) - A(1,4)*A(2,2)*A(3,3)

  A_inv(2,1) = A(2,1)*A(3,4)*A(4,3) + A(2,3)*A(3,1)*A(4,4) + A(2,4)*A(3,3)*A(4,1) - &
               A(2,1)*A(3,3)*A(4,4) - A(2,3)*A(3,4)*A(4,1) - A(2,4)*A(3,1)*A(4,3)
  A_inv(2,2) = A(1,1)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,1) + A(1,4)*A(3,1)*A(4,3) - &
               A(1,1)*A(3,4)*A(4,3) - A(1,3)*A(3,1)*A(4,4) - A(1,4)*A(3,3)*A(4,1)
  A_inv(2,3) = A(1,1)*A(2,4)*A(4,3) + A(1,3)*A(2,1)*A(4,4) + A(1,4)*A(2,3)*A(4,1) - &
               A(1,1)*A(2,3)*A(4,4) - A(1,3)*A(2,4)*A(4,1) - A(1,4)*A(2,1)*A(4,3)
  A_inv(2,4) = A(1,1)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,3) - &
               A(1,1)*A(2,4)*A(3,3) - A(1,3)*A(2,1)*A(3,4) - A(1,4)*A(2,3)*A(3,1)

  A_inv(3,1) = A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) - &
               A(2,1)*A(3,4)*A(4,2) - A(2,2)*A(3,1)*A(4,4) - A(2,4)*A(3,2)*A(4,1)
  A_inv(3,2) = A(1,1)*A(3,4)*A(4,2) + A(1,2)*A(3,1)*A(4,4) + A(1,4)*A(3,2)*A(4,1) - &
               A(1,1)*A(3,2)*A(4,4) - A(1,2)*A(3,4)*A(4,1) - A(1,4)*A(3,1)*A(4,2)
  A_inv(3,3) = A(1,1)*A(2,2)*A(4,4) + A(1,2)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,2) - &
               A(1,1)*A(2,4)*A(4,2) - A(1,2)*A(2,1)*A(4,4) - A(1,4)*A(2,2)*A(4,1)
  A_inv(3,4) = A(1,1)*A(2,4)*A(3,2) + A(1,2)*A(2,1)*A(3,4) + A(1,4)*A(2,2)*A(3,1) - &
               A(1,1)*A(2,2)*A(3,4) - A(1,2)*A(2,4)*A(3,1) - A(1,4)*A(2,1)*A(3,2)

  A_inv(4,1) = A(2,1)*A(3,3)*A(4,2) + A(2,2)*A(3,1)*A(4,3) + A(2,3)*A(3,2)*A(4,1) - &
               A(2,1)*A(3,2)*A(4,3) - A(2,2)*A(3,3)*A(4,1) - A(2,3)*A(3,1)*A(4,2)
  A_inv(4,2) = A(1,1)*A(3,2)*A(4,3) + A(1,2)*A(3,3)*A(4,1) + A(1,3)*A(3,1)*A(4,2) - &
               A(1,1)*A(3,3)*A(4,2) - A(1,2)*A(3,1)*A(4,3) - A(1,3)*A(3,2)*A(4,1)
  A_inv(4,3) = A(1,1)*A(2,3)*A(4,2) + A(1,2)*A(2,1)*A(4,3) + A(1,3)*A(2,2)*A(4,1) - &
               A(1,1)*A(2,2)*A(4,3) - A(1,2)*A(2,3)*A(4,1) - A(1,3)*A(2,1)*A(4,2)
  A_inv(4,4) = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - &
               A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) - A(1,3)*A(2,2)*A(3,1)

  A_inv = A_inv / determinant

  AA_inv = real(A_inv, kind=real4)

end subroutine INVERT_4x4

!--------------------------------------------------------------------------
! subroutine INVERT_5x5(A,A_inv,ierr)
!
! Matrix Inversion for a 5x5 matrix
!
! https://caps.gsfc.nasa.gov/simpson/software/m55inv_f90.txt
!--------------------------------------------------------------------------
subroutine INVERT_5x5(AA,AA_inv,ierr)
  real, dimension(:,:), intent(in):: AA
  real, dimension(:,:), intent(out):: AA_inv
  integer, intent(out):: ierr
  real(kind=real8), dimension(5,5):: A, A_inv

  LOGICAL:: OK_FLAG

  DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
  DOUBLE PRECISION :: DET, A11, A12, A13, A14, A15, A21, A22, A23, A24, &
         A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45,   &
         A51, A52, A53, A54, A55
  DOUBLE PRECISION, DIMENSION(5,5) :: COFACTOR


  !--- convert input to double precision
  A = dble(AA)

  !--- make shortnames for convenience
      A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5)
      A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5)
      A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5)
      A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5)
      A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5)

  !--- compute determinant
      DET = A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+    &
         A13*A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51-       &
         A15*A24*A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51-       &
         A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51+       &
         A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51+       &
         A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-       &
         A14*A23*A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-       &
         A12*A24*A33*A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-       &
         A15*A24*A33*A41*A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52-       &
         A13*A25*A34*A41*A52-A14*A23*A35*A41*A52+A13*A24*A35*A41*A52+       &
         A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21*A34*A43*A52+       &
         A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43*A52-       &
         A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52-       &
         A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52+       &
         A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+       &
         A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52+       &
         A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53+       &
         A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*A41*A53-       &
         A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*A53-       &
         A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+       &
         A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+       &
         A11*A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-       &
         A14*A22*A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53-       &
         A11*A24*A32*A45*A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53-       &
         A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22*A33*A41*A54-       &
         A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41*A54+       &
         A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54+       &
         A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54-       &
         A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-       &
         A11*A25*A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*A43*A54+       &
         A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*A54+       &
         A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*A33*A45*A54+       &
         A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*A55+       &
         A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-       &
         A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-       &
         A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+       &
         A14*A22*A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+       &
         A11*A24*A32*A43*A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-       &
         A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*A32*A44*A55-       &
         A11*A23*A32*A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55

      !--- check determinant for singularity
      IF (ABS(DET) .LE. EPS) THEN
         A_INV = 0.0D0
         AA_INV = real(A_Inv)
         OK_FLAG = .FALSE.
         ierr = 1
         RETURN
      END IF

      COFACTOR(1,1) = A25*A34*A43*A52-A24*A35*A43*A52-A25*A33*A44*A52+      &
         A23*A35*A44*A52+A24*A33*A45*A52-A23*A34*A45*A52-A25*A34*A42*A53+   &
         A24*A35*A42*A53+A25*A32*A44*A53-A22*A35*A44*A53-A24*A32*A45*A53+   &
         A22*A34*A45*A53+A25*A33*A42*A54-A23*A35*A42*A54-A25*A32*A43*A54+   &
         A22*A35*A43*A54+A23*A32*A45*A54-A22*A33*A45*A54-A24*A33*A42*A55+   &
         A23*A34*A42*A55+A24*A32*A43*A55-A22*A34*A43*A55-A23*A32*A44*A55+   &
         A22*A33*A44*A55

      COFACTOR(2,1) = -A15*A34*A43*A52+A14*A35*A43*A52+A15*A33*A44*A52-     &
         A13*A35*A44*A52-A14*A33*A45*A52+A13*A34*A45*A52+A15*A34*A42*A53-   &
         A14*A35*A42*A53-A15*A32*A44*A53+A12*A35*A44*A53+A14*A32*A45*A53-   &
         A12*A34*A45*A53-A15*A33*A42*A54+A13*A35*A42*A54+A15*A32*A43*A54-   &
         A12*A35*A43*A54-A13*A32*A45*A54+A12*A33*A45*A54+A14*A33*A42*A55-   &
         A13*A34*A42*A55-A14*A32*A43*A55+A12*A34*A43*A55+A13*A32*A44*A55-   &
         A12*A33*A44*A55

      COFACTOR(3,1) = A15*A24*A43*A52-A14*A25*A43*A52-A15*A23*A44*A52+      &
         A13*A25*A44*A52+A14*A23*A45*A52-A13*A24*A45*A52-A15*A24*A42*A53+   &
         A14*A25*A42*A53+A15*A22*A44*A53-A12*A25*A44*A53-A14*A22*A45*A53+   &
         A12*A24*A45*A53+A15*A23*A42*A54-A13*A25*A42*A54-A15*A22*A43*A54+   &
         A12*A25*A43*A54+A13*A22*A45*A54-A12*A23*A45*A54-A14*A23*A42*A55+   &
         A13*A24*A42*A55+A14*A22*A43*A55-A12*A24*A43*A55-A13*A22*A44*A55+   &
         A12*A23*A44*A55

      COFACTOR(4,1) = -A15*A24*A33*A52+A14*A25*A33*A52+A15*A23*A34*A52-     &
         A13*A25*A34*A52-A14*A23*A35*A52+A13*A24*A35*A52+A15*A24*A32*A53-   &
         A14*A25*A32*A53-A15*A22*A34*A53+A12*A25*A34*A53+A14*A22*A35*A53-   &
         A12*A24*A35*A53-A15*A23*A32*A54+A13*A25*A32*A54+A15*A22*A33*A54-   &
         A12*A25*A33*A54-A13*A22*A35*A54+A12*A23*A35*A54+A14*A23*A32*A55-   &
         A13*A24*A32*A55-A14*A22*A33*A55+A12*A24*A33*A55+A13*A22*A34*A55-   &
         A12*A23*A34*A55

      COFACTOR(5,1) = A15*A24*A33*A42-A14*A25*A33*A42-A15*A23*A34*A42+      &
         A13*A25*A34*A42+A14*A23*A35*A42-A13*A24*A35*A42-A15*A24*A32*A43+   &
         A14*A25*A32*A43+A15*A22*A34*A43-A12*A25*A34*A43-A14*A22*A35*A43+   &
         A12*A24*A35*A43+A15*A23*A32*A44-A13*A25*A32*A44-A15*A22*A33*A44+   &
         A12*A25*A33*A44+A13*A22*A35*A44-A12*A23*A35*A44-A14*A23*A32*A45+   &
         A13*A24*A32*A45+A14*A22*A33*A45-A12*A24*A33*A45-A13*A22*A34*A45+   &
         A12*A23*A34*A45

      COFACTOR(1,2) = -A25*A34*A43*A51+A24*A35*A43*A51+A25*A33*A44*A51-     &
         A23*A35*A44*A51-A24*A33*A45*A51+A23*A34*A45*A51+A25*A34*A41*A53-   &
         A24*A35*A41*A53-A25*A31*A44*A53+A21*A35*A44*A53+A24*A31*A45*A53-   &
         A21*A34*A45*A53-A25*A33*A41*A54+A23*A35*A41*A54+A25*A31*A43*A54-   &
         A21*A35*A43*A54-A23*A31*A45*A54+A21*A33*A45*A54+A24*A33*A41*A55-   &
         A23*A34*A41*A55-A24*A31*A43*A55+A21*A34*A43*A55+A23*A31*A44*A55-   &
         A21*A33*A44*A55

      COFACTOR(2,2) = A15*A34*A43*A51-A14*A35*A43*A51-A15*A33*A44*A51+      &
         A13*A35*A44*A51+A14*A33*A45*A51-A13*A34*A45*A51-A15*A34*A41*A53+   &
         A14*A35*A41*A53+A15*A31*A44*A53-A11*A35*A44*A53-A14*A31*A45*A53+   &
         A11*A34*A45*A53+A15*A33*A41*A54-A13*A35*A41*A54-A15*A31*A43*A54+   &
         A11*A35*A43*A54+A13*A31*A45*A54-A11*A33*A45*A54-A14*A33*A41*A55+   &
         A13*A34*A41*A55+A14*A31*A43*A55-A11*A34*A43*A55-A13*A31*A44*A55+   &
         A11*A33*A44*A55

      COFACTOR(3,2) = -A15*A24*A43*A51+A14*A25*A43*A51+A15*A23*A44*A51-     &
         A13*A25*A44*A51-A14*A23*A45*A51+A13*A24*A45*A51+A15*A24*A41*A53-   &
         A14*A25*A41*A53-A15*A21*A44*A53+A11*A25*A44*A53+A14*A21*A45*A53-   &
         A11*A24*A45*A53-A15*A23*A41*A54+A13*A25*A41*A54+A15*A21*A43*A54-   &
         A11*A25*A43*A54-A13*A21*A45*A54+A11*A23*A45*A54+A14*A23*A41*A55-   &
         A13*A24*A41*A55-A14*A21*A43*A55+A11*A24*A43*A55+A13*A21*A44*A55-   &
         A11*A23*A44*A55

      COFACTOR(4,2) = A15*A24*A33*A51-A14*A25*A33*A51-A15*A23*A34*A51+      &
         A13*A25*A34*A51+A14*A23*A35*A51-A13*A24*A35*A51-A15*A24*A31*A53+   &
         A14*A25*A31*A53+A15*A21*A34*A53-A11*A25*A34*A53-A14*A21*A35*A53+   &
         A11*A24*A35*A53+A15*A23*A31*A54-A13*A25*A31*A54-A15*A21*A33*A54+   &
         A11*A25*A33*A54+A13*A21*A35*A54-A11*A23*A35*A54-A14*A23*A31*A55+   &
         A13*A24*A31*A55+A14*A21*A33*A55-A11*A24*A33*A55-A13*A21*A34*A55+   &
         A11*A23*A34*A55

      COFACTOR(5,2) = -A15*A24*A33*A41+A14*A25*A33*A41+A15*A23*A34*A41-     &
         A13*A25*A34*A41-A14*A23*A35*A41+A13*A24*A35*A41+A15*A24*A31*A43-   &
         A14*A25*A31*A43-A15*A21*A34*A43+A11*A25*A34*A43+A14*A21*A35*A43-   &
         A11*A24*A35*A43-A15*A23*A31*A44+A13*A25*A31*A44+A15*A21*A33*A44-   &
         A11*A25*A33*A44-A13*A21*A35*A44+A11*A23*A35*A44+A14*A23*A31*A45-   &
         A13*A24*A31*A45-A14*A21*A33*A45+A11*A24*A33*A45+A13*A21*A34*A45-   &
         A11*A23*A34*A45

      COFACTOR(1,3) = A25*A34*A42*A51-A24*A35*A42*A51-A25*A32*A44*A51+      &
         A22*A35*A44*A51+A24*A32*A45*A51-A22*A34*A45*A51-A25*A34*A41*A52+   &
         A24*A35*A41*A52+A25*A31*A44*A52-A21*A35*A44*A52-A24*A31*A45*A52+   &
         A21*A34*A45*A52+A25*A32*A41*A54-A22*A35*A41*A54-A25*A31*A42*A54+   &
         A21*A35*A42*A54+A22*A31*A45*A54-A21*A32*A45*A54-A24*A32*A41*A55+   &
         A22*A34*A41*A55+A24*A31*A42*A55-A21*A34*A42*A55-A22*A31*A44*A55+   &
         A21*A32*A44*A55

      COFACTOR(2,3) = -A15*A34*A42*A51+A14*A35*A42*A51+A15*A32*A44*A51-     &
         A12*A35*A44*A51-A14*A32*A45*A51+A12*A34*A45*A51+A15*A34*A41*A52-   &
         A14*A35*A41*A52-A15*A31*A44*A52+A11*A35*A44*A52+A14*A31*A45*A52-   &
         A11*A34*A45*A52-A15*A32*A41*A54+A12*A35*A41*A54+A15*A31*A42*A54-   &
         A11*A35*A42*A54-A12*A31*A45*A54+A11*A32*A45*A54+A14*A32*A41*A55-   &
         A12*A34*A41*A55-A14*A31*A42*A55+A11*A34*A42*A55+A12*A31*A44*A55-   &
         A11*A32*A44*A55

      COFACTOR(3,3) = A15*A24*A42*A51-A14*A25*A42*A51-A15*A22*A44*A51+      &
         A12*A25*A44*A51+A14*A22*A45*A51-A12*A24*A45*A51-A15*A24*A41*A52+   &
         A14*A25*A41*A52+A15*A21*A44*A52-A11*A25*A44*A52-A14*A21*A45*A52+   &
         A11*A24*A45*A52+A15*A22*A41*A54-A12*A25*A41*A54-A15*A21*A42*A54+   &
         A11*A25*A42*A54+A12*A21*A45*A54-A11*A22*A45*A54-A14*A22*A41*A55+   &
         A12*A24*A41*A55+A14*A21*A42*A55-A11*A24*A42*A55-A12*A21*A44*A55+   &
         A11*A22*A44*A55

      COFACTOR(4,3) = -A15*A24*A32*A51+A14*A25*A32*A51+A15*A22*A34*A51-     &
         A12*A25*A34*A51-A14*A22*A35*A51+A12*A24*A35*A51+A15*A24*A31*A52-   &
         A14*A25*A31*A52-A15*A21*A34*A52+A11*A25*A34*A52+A14*A21*A35*A52-   &
         A11*A24*A35*A52-A15*A22*A31*A54+A12*A25*A31*A54+A15*A21*A32*A54-   &
         A11*A25*A32*A54-A12*A21*A35*A54+A11*A22*A35*A54+A14*A22*A31*A55-   &
         A12*A24*A31*A55-A14*A21*A32*A55+A11*A24*A32*A55+A12*A21*A34*A55-   &
         A11*A22*A34*A55

      COFACTOR(5,3) = A15*A24*A32*A41-A14*A25*A32*A41-A15*A22*A34*A41+      &
         A12*A25*A34*A41+A14*A22*A35*A41-A12*A24*A35*A41-A15*A24*A31*A42+   &
         A14*A25*A31*A42+A15*A21*A34*A42-A11*A25*A34*A42-A14*A21*A35*A42+   &
         A11*A24*A35*A42+A15*A22*A31*A44-A12*A25*A31*A44-A15*A21*A32*A44+   &
         A11*A25*A32*A44+A12*A21*A35*A44-A11*A22*A35*A44-A14*A22*A31*A45+   &
         A12*A24*A31*A45+A14*A21*A32*A45-A11*A24*A32*A45-A12*A21*A34*A45+   &
         A11*A22*A34*A45

      COFACTOR(1,4) = -A25*A33*A42*A51+A23*A35*A42*A51+A25*A32*A43*A51-     &
         A22*A35*A43*A51-A23*A32*A45*A51+A22*A33*A45*A51+A25*A33*A41*A52-   &
         A23*A35*A41*A52-A25*A31*A43*A52+A21*A35*A43*A52+A23*A31*A45*A52-   &
         A21*A33*A45*A52-A25*A32*A41*A53+A22*A35*A41*A53+A25*A31*A42*A53-   &
         A21*A35*A42*A53-A22*A31*A45*A53+A21*A32*A45*A53+A23*A32*A41*A55-   &
         A22*A33*A41*A55-A23*A31*A42*A55+A21*A33*A42*A55+A22*A31*A43*A55-   &
         A21*A32*A43*A55

      COFACTOR(2,4) = A15*A33*A42*A51-A13*A35*A42*A51-A15*A32*A43*A51+      &
         A12*A35*A43*A51+A13*A32*A45*A51-A12*A33*A45*A51-A15*A33*A41*A52+   &
         A13*A35*A41*A52+A15*A31*A43*A52-A11*A35*A43*A52-A13*A31*A45*A52+   &
         A11*A33*A45*A52+A15*A32*A41*A53-A12*A35*A41*A53-A15*A31*A42*A53+   &
         A11*A35*A42*A53+A12*A31*A45*A53-A11*A32*A45*A53-A13*A32*A41*A55+   &
         A12*A33*A41*A55+A13*A31*A42*A55-A11*A33*A42*A55-A12*A31*A43*A55+   &
         A11*A32*A43*A55

      COFACTOR(3,4) = -A15*A23*A42*A51+A13*A25*A42*A51+A15*A22*A43*A51-     &
         A12*A25*A43*A51-A13*A22*A45*A51+A12*A23*A45*A51+A15*A23*A41*A52-   &
         A13*A25*A41*A52-A15*A21*A43*A52+A11*A25*A43*A52+A13*A21*A45*A52-   &
         A11*A23*A45*A52-A15*A22*A41*A53+A12*A25*A41*A53+A15*A21*A42*A53-   &
         A11*A25*A42*A53-A12*A21*A45*A53+A11*A22*A45*A53+A13*A22*A41*A55-   &
         A12*A23*A41*A55-A13*A21*A42*A55+A11*A23*A42*A55+A12*A21*A43*A55-   &
         A11*A22*A43*A55

      COFACTOR(4,4) = A15*A23*A32*A51-A13*A25*A32*A51-A15*A22*A33*A51+      &
         A12*A25*A33*A51+A13*A22*A35*A51-A12*A23*A35*A51-A15*A23*A31*A52+   &
         A13*A25*A31*A52+A15*A21*A33*A52-A11*A25*A33*A52-A13*A21*A35*A52+   &
         A11*A23*A35*A52+A15*A22*A31*A53-A12*A25*A31*A53-A15*A21*A32*A53+   &
         A11*A25*A32*A53+A12*A21*A35*A53-A11*A22*A35*A53-A13*A22*A31*A55+   &
         A12*A23*A31*A55+A13*A21*A32*A55-A11*A23*A32*A55-A12*A21*A33*A55+   &
         A11*A22*A33*A55

      COFACTOR(5,4) = -A15*A23*A32*A41+A13*A25*A32*A41+A15*A22*A33*A41-     &
         A12*A25*A33*A41-A13*A22*A35*A41+A12*A23*A35*A41+A15*A23*A31*A42-   &
         A13*A25*A31*A42-A15*A21*A33*A42+A11*A25*A33*A42+A13*A21*A35*A42-   &
         A11*A23*A35*A42-A15*A22*A31*A43+A12*A25*A31*A43+A15*A21*A32*A43-   &
         A11*A25*A32*A43-A12*A21*A35*A43+A11*A22*A35*A43+A13*A22*A31*A45-   &
         A12*A23*A31*A45-A13*A21*A32*A45+A11*A23*A32*A45+A12*A21*A33*A45-   &
         A11*A22*A33*A45

      COFACTOR(1,5) = A24*A33*A42*A51-A23*A34*A42*A51-A24*A32*A43*A51+      &
         A22*A34*A43*A51+A23*A32*A44*A51-A22*A33*A44*A51-A24*A33*A41*A52+   &
         A23*A34*A41*A52+A24*A31*A43*A52-A21*A34*A43*A52-A23*A31*A44*A52+   &
         A21*A33*A44*A52+A24*A32*A41*A53-A22*A34*A41*A53-A24*A31*A42*A53+   &
         A21*A34*A42*A53+A22*A31*A44*A53-A21*A32*A44*A53-A23*A32*A41*A54+   &
         A22*A33*A41*A54+A23*A31*A42*A54-A21*A33*A42*A54-A22*A31*A43*A54+   &
         A21*A32*A43*A54

      COFACTOR(2,5) = -A14*A33*A42*A51+A13*A34*A42*A51+A14*A32*A43*A51-     &
         A12*A34*A43*A51-A13*A32*A44*A51+A12*A33*A44*A51+A14*A33*A41*A52-   &
         A13*A34*A41*A52-A14*A31*A43*A52+A11*A34*A43*A52+A13*A31*A44*A52-   &
         A11*A33*A44*A52-A14*A32*A41*A53+A12*A34*A41*A53+A14*A31*A42*A53-   &
         A11*A34*A42*A53-A12*A31*A44*A53+A11*A32*A44*A53+A13*A32*A41*A54-   &
         A12*A33*A41*A54-A13*A31*A42*A54+A11*A33*A42*A54+A12*A31*A43*A54-   &
         A11*A32*A43*A54

      COFACTOR(3,5) = A14*A23*A42*A51-A13*A24*A42*A51-A14*A22*A43*A51+      &
         A12*A24*A43*A51+A13*A22*A44*A51-A12*A23*A44*A51-A14*A23*A41*A52+   &
         A13*A24*A41*A52+A14*A21*A43*A52-A11*A24*A43*A52-A13*A21*A44*A52+   &
         A11*A23*A44*A52+A14*A22*A41*A53-A12*A24*A41*A53-A14*A21*A42*A53+   &
         A11*A24*A42*A53+A12*A21*A44*A53-A11*A22*A44*A53-A13*A22*A41*A54+   &
         A12*A23*A41*A54+A13*A21*A42*A54-A11*A23*A42*A54-A12*A21*A43*A54+   &
         A11*A22*A43*A54

      COFACTOR(4,5) = -A14*A23*A32*A51+A13*A24*A32*A51+A14*A22*A33*A51-     &
         A12*A24*A33*A51-A13*A22*A34*A51+A12*A23*A34*A51+A14*A23*A31*A52-   &
         A13*A24*A31*A52-A14*A21*A33*A52+A11*A24*A33*A52+A13*A21*A34*A52-   &
         A11*A23*A34*A52-A14*A22*A31*A53+A12*A24*A31*A53+A14*A21*A32*A53-   &
         A11*A24*A32*A53-A12*A21*A34*A53+A11*A22*A34*A53+A13*A22*A31*A54-   &
         A12*A23*A31*A54-A13*A21*A32*A54+A11*A23*A32*A54+A12*A21*A33*A54-   &
         A11*A22*A33*A54

      COFACTOR(5,5) = A14*A23*A32*A41-A13*A24*A32*A41-A14*A22*A33*A41+      &
         A12*A24*A33*A41+A13*A22*A34*A41-A12*A23*A34*A41-A14*A23*A31*A42+   &
         A13*A24*A31*A42+A14*A21*A33*A42-A11*A24*A33*A42-A13*A21*A34*A42+   &
         A11*A23*A34*A42+A14*A22*A31*A43-A12*A24*A31*A43-A14*A21*A32*A43+   &
         A11*A24*A32*A43+A12*A21*A34*A43-A11*A22*A34*A43-A13*A22*A31*A44+   &
         A12*A23*A31*A44+A13*A21*A32*A44-A11*A23*A32*A44-A12*A21*A33*A44+   &
         A11*A22*A33*A44

      A_INV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      ierr = 0

      !--- convert back to real for output
      AA_INV = real(A_Inv)

end subroutine INVERT_5x5

!--------------------------------------------------------------------------
! subroutine INVERT_6x6(A,A_inv,ierr)
!
! Matrix Inversion for a 6x6 matrix
!
! https://caps.gsfc.nasa.gov/simpson/software/m66inv_f90.txt
!--------------------------------------------------------------------------
subroutine INVERT_6x6(AA,AA_inv,ierr)
  real, dimension(:,:), intent(in):: AA
  real, dimension(:,:), intent(out):: AA_inv
  integer, intent(out):: ierr
  real(kind=real8), dimension(6,6):: A, A_inv

  LOGICAL:: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET, A11, A12, A13, A14, A15, A16, A21, A22, A23, A24, &
         A25, A26, A31, A32, A33, A34, A35, A36, A41, A42, A43, A44, A45, A46, &
         A51, A52, A53, A54, A55, A56, A61, A62, A63, A64, A65, A66
      DOUBLE PRECISION, DIMENSION(6,6) :: COFACTOR

      A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5); A16=A(1,6)
      A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5); A26=A(2,6)
      A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5); A36=A(3,6)
      A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5); A46=A(4,6)
      A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5); A56=A(5,6)
      A61=A(6,1); A62=A(6,2); A63=A(6,3); A64=A(6,4); A65=A(6,5); A66=A(6,6)

      DET = -(A16*A25*A34*A43*A52-A15*A26*A34*A43*A52-A16*A24*A35*A43* &
         A52+A14*A26*A35*A43*A52+A15*A24*A36*A43*A52-A14*A25*A36*A43*A52-A16*A25* &
         A33*A44*A52+A15*A26*A33*A44*A52+A16*A23*A35*A44*A52-A13*A26*A35*A44* &
         A52-A15*A23*A36*A44*A52+A13*A25*A36*A44*A52+A16*A24*A33*A45*A52-A14*A26* &
         A33*A45*A52-A16*A23*A34*A45*A52+A13*A26*A34*A45*A52+A14*A23*A36*A45* &
         A52-A13*A24*A36*A45*A52-A15*A24*A33*A46*A52+A14*A25*A33*A46*A52+A15*A23* &
         A34*A46*A52-A13*A25*A34*A46*A52-A14*A23*A35*A46*A52+A13*A24*A35*A46* &
         A52-A16*A25*A34*A42*A53+A15*A26*A34*A42*A53+A16*A24*A35*A42*A53-A14*A26* &
         A35*A42*A53-A15*A24*A36*A42*A53+A14*A25*A36*A42*A53+A16*A25*A32*A44* &
         A53-A15*A26*A32*A44*A53-A16*A22*A35*A44*A53+A12*A26*A35*A44*A53+A15*A22* &
         A36*A44*A53-A12*A25*A36*A44*A53-A16*A24*A32*A45*A53+A14*A26*A32*A45* &
         A53+A16*A22*A34*A45*A53-A12*A26*A34*A45*A53-A14*A22*A36*A45*A53+A12*A24* &
         A36*A45*A53+A15*A24*A32*A46*A53-A14*A25*A32*A46*A53-A15*A22*A34*A46* &
         A53+A12*A25*A34*A46*A53+A14*A22*A35*A46*A53-A12*A24*A35*A46*A53+A16*A25* &
         A33*A42*A54-A15*A26*A33*A42*A54-A16*A23*A35*A42*A54+A13*A26*A35*A42* &
         A54+A15*A23*A36*A42*A54-A13*A25*A36*A42*A54-A16*A25*A32*A43*A54+A15*A26* &
         A32*A43*A54+A16*A22*A35*A43*A54-A12*A26*A35*A43*A54-A15*A22*A36*A43* &
         A54+A12*A25*A36*A43*A54+A16*A23*A32*A45*A54-A13*A26*A32*A45*A54-A16*A22* &
         A33*A45*A54+A12*A26*A33*A45*A54+A13*A22*A36*A45*A54-A12*A23*A36*A45* &
         A54-A15*A23*A32*A46*A54+A13*A25*A32*A46*A54+A15*A22*A33*A46*A54-A12*A25* &
         A33*A46*A54-A13*A22*A35*A46*A54+A12*A23*A35*A46*A54-A16*A24*A33*A42* &
         A55+A14*A26*A33*A42*A55+A16*A23*A34*A42*A55-A13*A26*A34*A42*A55-A14*A23* &
         A36*A42*A55+A13*A24*A36*A42*A55+A16*A24*A32*A43*A55-A14*A26*A32*A43* &
         A55-A16*A22*A34*A43*A55+A12*A26*A34*A43*A55+A14*A22*A36*A43*A55-A12*A24* &
         A36*A43*A55-A16*A23*A32*A44*A55+A13*A26*A32*A44*A55+A16*A22*A33*A44* &
         A55-A12*A26*A33*A44*A55-A13*A22*A36*A44*A55+A12*A23*A36*A44*A55+A14*A23* &
         A32*A46*A55-A13*A24*A32*A46*A55-A14*A22*A33*A46*A55+A12*A24*A33*A46* &
         A55+A13*A22*A34*A46*A55-A12*A23*A34*A46*A55+A15*A24*A33*A42*A56-A14*A25* &
         A33*A42*A56-A15*A23*A34*A42*A56+A13*A25*A34*A42*A56+A14*A23*A35*A42* &
         A56-A13*A24*A35*A42*A56-A15*A24*A32*A43*A56+A14*A25*A32*A43*A56+A15*A22* &
         A34*A43*A56-A12*A25*A34*A43*A56-A14*A22*A35*A43*A56+A12*A24*A35*A43* &
         A56+A15*A23*A32*A44*A56-A13*A25*A32*A44*A56-A15*A22*A33*A44*A56+A12*A25* &
         A33*A44*A56+A13*A22*A35*A44*A56-A12*A23*A35*A44*A56-A14*A23*A32*A45* &
         A56+A13*A24*A32*A45*A56+A14*A22*A33*A45*A56-A12*A24*A33*A45*A56-A13*A22* &
         A34*A45*A56+A12*A23*A34*A45*A56)*A61+(A16*A25*A34*A43*A51-A15*A26*A34* &
         A43*A51-A16*A24*A35*A43*A51+A14*A26*A35*A43*A51+A15*A24*A36*A43*A51-A14* &
         A25*A36*A43*A51-A16*A25*A33*A44*A51+A15*A26*A33*A44*A51+A16*A23*A35*A44* &
         A51-A13*A26*A35*A44*A51-A15*A23*A36*A44*A51+A13*A25*A36*A44*A51+A16*A24* &
         A33*A45*A51-A14*A26*A33*A45*A51-A16*A23*A34*A45*A51+A13*A26*A34*A45* &
         A51+A14*A23*A36*A45*A51-A13*A24*A36*A45*A51-A15*A24*A33*A46*A51+A14*A25* &
         A33*A46*A51+A15*A23*A34*A46*A51-A13*A25*A34*A46*A51-A14*A23*A35*A46* &
         A51+A13*A24*A35*A46*A51-A16*A25*A34*A41*A53+A15*A26*A34*A41*A53+A16*A24* &
         A35*A41*A53-A14*A26*A35*A41*A53-A15*A24*A36*A41*A53+A14*A25*A36*A41* &
         A53+A16*A25*A31*A44*A53-A15*A26*A31*A44*A53-A16*A21*A35*A44*A53+A11*A26* &
         A35*A44*A53+A15*A21*A36*A44*A53-A11*A25*A36*A44*A53-A16*A24*A31*A45* &
         A53+A14*A26*A31*A45*A53+A16*A21*A34*A45*A53-A11*A26*A34*A45*A53-A14*A21* &
         A36*A45*A53+A11*A24*A36*A45*A53+A15*A24*A31*A46*A53-A14*A25*A31*A46* &
         A53-A15*A21*A34*A46*A53+A11*A25*A34*A46*A53+A14*A21*A35*A46*A53-A11*A24* &
         A35*A46*A53+A16*A25*A33*A41*A54-A15*A26*A33*A41*A54-A16*A23*A35*A41* &
         A54+A13*A26*A35*A41*A54+A15*A23*A36*A41*A54-A13*A25*A36*A41*A54-A16*A25* &
         A31*A43*A54+A15*A26*A31*A43*A54+A16*A21*A35*A43*A54-A11*A26*A35*A43* &
         A54-A15*A21*A36*A43*A54+A11*A25*A36*A43*A54+A16*A23*A31*A45*A54-A13*A26* &
         A31*A45*A54-A16*A21*A33*A45*A54+A11*A26*A33*A45*A54+A13*A21*A36*A45* &
         A54-A11*A23*A36*A45*A54-A15*A23*A31*A46*A54+A13*A25*A31*A46*A54+A15*A21* &
         A33*A46*A54-A11*A25*A33*A46*A54-A13*A21*A35*A46*A54+A11*A23*A35*A46* &
         A54-A16*A24*A33*A41*A55+A14*A26*A33*A41*A55+A16*A23*A34*A41*A55-A13*A26* &
         A34*A41*A55-A14*A23*A36*A41*A55+A13*A24*A36*A41*A55+A16*A24*A31*A43* &
         A55-A14*A26*A31*A43*A55-A16*A21*A34*A43*A55+A11*A26*A34*A43*A55+A14*A21* &
         A36*A43*A55-A11*A24*A36*A43*A55-A16*A23*A31*A44*A55+A13*A26*A31*A44* &
         A55+A16*A21*A33*A44*A55-A11*A26*A33*A44*A55-A13*A21*A36*A44*A55+A11*A23* &
         A36*A44*A55+A14*A23*A31*A46*A55-A13*A24*A31*A46*A55-A14*A21*A33*A46* &
         A55+A11*A24*A33*A46*A55+A13*A21*A34*A46*A55-A11*A23*A34*A46*A55+A15*A24* &
         A33*A41*A56-A14*A25*A33*A41*A56-A15*A23*A34*A41*A56+A13*A25*A34*A41* &
         A56+A14*A23*A35*A41*A56-A13*A24*A35*A41*A56-A15*A24*A31*A43*A56+A14*A25* &
         A31*A43*A56+A15*A21*A34*A43*A56-A11*A25*A34*A43*A56-A14*A21*A35*A43* &
         A56+A11*A24*A35*A43*A56+A15*A23*A31*A44*A56-A13*A25*A31*A44*A56-A15*A21* &
         A33*A44*A56+A11*A25*A33*A44*A56+A13*A21*A35*A44*A56-A11*A23*A35*A44* &
         A56-A14*A23*A31*A45*A56+A13*A24*A31*A45*A56+A14*A21*A33*A45*A56-A11*A24* &
         A33*A45*A56-A13*A21*A34*A45*A56+A11*A23*A34*A45*A56)*A62-(A16*A25*A34* &
         A42*A51-A15*A26*A34*A42*A51-A16*A24*A35*A42*A51+A14*A26*A35*A42*A51+A15* &
         A24*A36*A42*A51-A14*A25*A36*A42*A51-A16*A25*A32*A44*A51+A15*A26*A32*A44* &
         A51+A16*A22*A35*A44*A51-A12*A26*A35*A44*A51-A15*A22*A36*A44*A51+A12*A25* &
         A36*A44*A51+A16*A24*A32*A45*A51-A14*A26*A32*A45*A51-A16*A22*A34*A45* &
         A51+A12*A26*A34*A45*A51+A14*A22*A36*A45*A51-A12*A24*A36*A45*A51-A15*A24* &
         A32*A46*A51+A14*A25*A32*A46*A51+A15*A22*A34*A46*A51-A12*A25*A34*A46* &
         A51-A14*A22*A35*A46*A51+A12*A24*A35*A46*A51-A16*A25*A34*A41*A52+A15*A26* &
         A34*A41*A52+A16*A24*A35*A41*A52-A14*A26*A35*A41*A52-A15*A24*A36*A41* &
         A52+A14*A25*A36*A41*A52+A16*A25*A31*A44*A52-A15*A26*A31*A44*A52-A16*A21* &
         A35*A44*A52+A11*A26*A35*A44*A52+A15*A21*A36*A44*A52-A11*A25*A36*A44* &
         A52-A16*A24*A31*A45*A52+A14*A26*A31*A45*A52+A16*A21*A34*A45*A52-A11*A26* &
         A34*A45*A52-A14*A21*A36*A45*A52+A11*A24*A36*A45*A52+A15*A24*A31*A46* &
         A52-A14*A25*A31*A46*A52-A15*A21*A34*A46*A52+A11*A25*A34*A46*A52+A14*A21* &
         A35*A46*A52-A11*A24*A35*A46*A52+A16*A25*A32*A41*A54-A15*A26*A32*A41* &
         A54-A16*A22*A35*A41*A54+A12*A26*A35*A41*A54+A15*A22*A36*A41*A54-A12*A25* &
         A36*A41*A54-A16*A25*A31*A42*A54+A15*A26*A31*A42*A54+A16*A21*A35*A42* &
         A54-A11*A26*A35*A42*A54-A15*A21*A36*A42*A54+A11*A25*A36*A42*A54+A16*A22* &
         A31*A45*A54-A12*A26*A31*A45*A54-A16*A21*A32*A45*A54+A11*A26*A32*A45* &
         A54+A12*A21*A36*A45*A54-A11*A22*A36*A45*A54-A15*A22*A31*A46*A54+A12*A25* &
         A31*A46*A54+A15*A21*A32*A46*A54-A11*A25*A32*A46*A54-A12*A21*A35*A46* &
         A54+A11*A22*A35*A46*A54-A16*A24*A32*A41*A55+A14*A26*A32*A41*A55+A16*A22* &
         A34*A41*A55-A12*A26*A34*A41*A55-A14*A22*A36*A41*A55+A12*A24*A36*A41* &
         A55+A16*A24*A31*A42*A55-A14*A26*A31*A42*A55-A16*A21*A34*A42*A55+A11*A26* &
         A34*A42*A55+A14*A21*A36*A42*A55-A11*A24*A36*A42*A55-A16*A22*A31*A44* &
         A55+A12*A26*A31*A44*A55+A16*A21*A32*A44*A55-A11*A26*A32*A44*A55-A12*A21* &
         A36*A44*A55+A11*A22*A36*A44*A55+A14*A22*A31*A46*A55-A12*A24*A31*A46* &
         A55-A14*A21*A32*A46*A55+A11*A24*A32*A46*A55+A12*A21*A34*A46*A55-A11*A22* &
         A34*A46*A55+A15*A24*A32*A41*A56-A14*A25*A32*A41*A56-A15*A22*A34*A41* &
         A56+A12*A25*A34*A41*A56+A14*A22*A35*A41*A56-A12*A24*A35*A41*A56-A15*A24* &
         A31*A42*A56+A14*A25*A31*A42*A56+A15*A21*A34*A42*A56-A11*A25*A34*A42* &
         A56-A14*A21*A35*A42*A56+A11*A24*A35*A42*A56+A15*A22*A31*A44*A56-A12*A25* &
         A31*A44*A56-A15*A21*A32*A44*A56+A11*A25*A32*A44*A56+A12*A21*A35*A44* &
         A56-A11*A22*A35*A44*A56-A14*A22*A31*A45*A56+A12*A24*A31*A45*A56+A14*A21* &
         A32*A45*A56-A11*A24*A32*A45*A56-A12*A21*A34*A45*A56+A11*A22*A34*A45*A56)* &
         A63+(A16*A25*A33*A42*A51-A15*A26*A33*A42*A51-A16*A23*A35*A42*A51+A13*A26* &
         A35*A42*A51+A15*A23*A36*A42*A51-A13*A25*A36*A42*A51-A16*A25*A32*A43* &
         A51+A15*A26*A32*A43*A51+A16*A22*A35*A43*A51-A12*A26*A35*A43*A51-A15*A22* &
         A36*A43*A51+A12*A25*A36*A43*A51+A16*A23*A32*A45*A51-A13*A26*A32*A45* &
         A51-A16*A22*A33*A45*A51+A12*A26*A33*A45*A51+A13*A22*A36*A45*A51-A12*A23* &
         A36*A45*A51-A15*A23*A32*A46*A51+A13*A25*A32*A46*A51+A15*A22*A33*A46* &
         A51-A12*A25*A33*A46*A51-A13*A22*A35*A46*A51+A12*A23*A35*A46*A51-A16*A25* &
         A33*A41*A52+A15*A26*A33*A41*A52+A16*A23*A35*A41*A52-A13*A26*A35*A41* &
         A52-A15*A23*A36*A41*A52+A13*A25*A36*A41*A52+A16*A25*A31*A43*A52-A15*A26* &
         A31*A43*A52-A16*A21*A35*A43*A52+A11*A26*A35*A43*A52+A15*A21*A36*A43* &
         A52-A11*A25*A36*A43*A52-A16*A23*A31*A45*A52+A13*A26*A31*A45*A52+A16*A21* &
         A33*A45*A52-A11*A26*A33*A45*A52-A13*A21*A36*A45*A52+A11*A23*A36*A45* &
         A52+A15*A23*A31*A46*A52-A13*A25*A31*A46*A52-A15*A21*A33*A46*A52+A11*A25* &
         A33*A46*A52+A13*A21*A35*A46*A52-A11*A23*A35*A46*A52+A16*A25*A32*A41* &
         A53-A15*A26*A32*A41*A53-A16*A22*A35*A41*A53+A12*A26*A35*A41*A53+A15*A22* &
         A36*A41*A53-A12*A25*A36*A41*A53-A16*A25*A31*A42*A53+A15*A26*A31*A42* &
         A53+A16*A21*A35*A42*A53-A11*A26*A35*A42*A53-A15*A21*A36*A42*A53+A11*A25* &
         A36*A42*A53+A16*A22*A31*A45*A53-A12*A26*A31*A45*A53-A16*A21*A32*A45* &
         A53+A11*A26*A32*A45*A53+A12*A21*A36*A45*A53-A11*A22*A36*A45*A53-A15*A22* &
         A31*A46*A53+A12*A25*A31*A46*A53+A15*A21*A32*A46*A53-A11*A25*A32*A46* &
         A53-A12*A21*A35*A46*A53+A11*A22*A35*A46*A53-A16*A23*A32*A41*A55+A13*A26* &
         A32*A41*A55+A16*A22*A33*A41*A55-A12*A26*A33*A41*A55-A13*A22*A36*A41* &
         A55+A12*A23*A36*A41*A55+A16*A23*A31*A42*A55-A13*A26*A31*A42*A55-A16*A21* &
         A33*A42*A55+A11*A26*A33*A42*A55+A13*A21*A36*A42*A55-A11*A23*A36*A42* &
         A55-A16*A22*A31*A43*A55+A12*A26*A31*A43*A55+A16*A21*A32*A43*A55-A11*A26* &
         A32*A43*A55-A12*A21*A36*A43*A55+A11*A22*A36*A43*A55+A13*A22*A31*A46* &
         A55-A12*A23*A31*A46*A55-A13*A21*A32*A46*A55+A11*A23*A32*A46*A55+A12*A21* &
         A33*A46*A55-A11*A22*A33*A46*A55+A15*A23*A32*A41*A56-A13*A25*A32*A41* &
         A56-A15*A22*A33*A41*A56+A12*A25*A33*A41*A56+A13*A22*A35*A41*A56-A12*A23* &
         A35*A41*A56-A15*A23*A31*A42*A56+A13*A25*A31*A42*A56+A15*A21*A33*A42* &
         A56-A11*A25*A33*A42*A56-A13*A21*A35*A42*A56+A11*A23*A35*A42*A56+A15*A22* &
         A31*A43*A56-A12*A25*A31*A43*A56-A15*A21*A32*A43*A56+A11*A25*A32*A43* &
         A56+A12*A21*A35*A43*A56-A11*A22*A35*A43*A56-A13*A22*A31*A45*A56+A12*A23* &
         A31*A45*A56+A13*A21*A32*A45*A56-A11*A23*A32*A45*A56-A12*A21*A33*A45* &
         A56+A11*A22*A33*A45*A56)*A64-(A16*A24*A33*A42*A51-A14*A26*A33*A42* &
         A51-A16*A23*A34*A42*A51+A13*A26*A34*A42*A51+A14*A23*A36*A42*A51-A13*A24* &
         A36*A42*A51-A16*A24*A32*A43*A51+A14*A26*A32*A43*A51+A16*A22*A34*A43* &
         A51-A12*A26*A34*A43*A51-A14*A22*A36*A43*A51+A12*A24*A36*A43*A51+A16*A23* &
         A32*A44*A51-A13*A26*A32*A44*A51-A16*A22*A33*A44*A51+A12*A26*A33*A44* &
         A51+A13*A22*A36*A44*A51-A12*A23*A36*A44*A51-A14*A23*A32*A46*A51+A13*A24* &
         A32*A46*A51+A14*A22*A33*A46*A51-A12*A24*A33*A46*A51-A13*A22*A34*A46* &
         A51+A12*A23*A34*A46*A51-A16*A24*A33*A41*A52+A14*A26*A33*A41*A52+A16*A23* &
         A34*A41*A52-A13*A26*A34*A41*A52-A14*A23*A36*A41*A52+A13*A24*A36*A41* &
         A52+A16*A24*A31*A43*A52-A14*A26*A31*A43*A52-A16*A21*A34*A43*A52+A11*A26* &
         A34*A43*A52+A14*A21*A36*A43*A52-A11*A24*A36*A43*A52-A16*A23*A31*A44* &
         A52+A13*A26*A31*A44*A52+A16*A21*A33*A44*A52-A11*A26*A33*A44*A52-A13*A21* &
         A36*A44*A52+A11*A23*A36*A44*A52+A14*A23*A31*A46*A52-A13*A24*A31*A46* &
         A52-A14*A21*A33*A46*A52+A11*A24*A33*A46*A52+A13*A21*A34*A46*A52-A11*A23* &
         A34*A46*A52+A16*A24*A32*A41*A53-A14*A26*A32*A41*A53-A16*A22*A34*A41* &
         A53+A12*A26*A34*A41*A53+A14*A22*A36*A41*A53-A12*A24*A36*A41*A53-A16*A24* &
         A31*A42*A53+A14*A26*A31*A42*A53+A16*A21*A34*A42*A53-A11*A26*A34*A42* &
         A53-A14*A21*A36*A42*A53+A11*A24*A36*A42*A53+A16*A22*A31*A44*A53-A12*A26* &
         A31*A44*A53-A16*A21*A32*A44*A53+A11*A26*A32*A44*A53+A12*A21*A36*A44* &
         A53-A11*A22*A36*A44*A53-A14*A22*A31*A46*A53+A12*A24*A31*A46*A53+A14*A21* &
         A32*A46*A53-A11*A24*A32*A46*A53-A12*A21*A34*A46*A53+A11*A22*A34*A46* &
         A53-A16*A23*A32*A41*A54+A13*A26*A32*A41*A54+A16*A22*A33*A41*A54-A12*A26* &
         A33*A41*A54-A13*A22*A36*A41*A54+A12*A23*A36*A41*A54+A16*A23*A31*A42* &
         A54-A13*A26*A31*A42*A54-A16*A21*A33*A42*A54+A11*A26*A33*A42*A54+A13*A21* &
         A36*A42*A54-A11*A23*A36*A42*A54-A16*A22*A31*A43*A54+A12*A26*A31*A43* &
         A54+A16*A21*A32*A43*A54-A11*A26*A32*A43*A54-A12*A21*A36*A43*A54+A11*A22* &
         A36*A43*A54+A13*A22*A31*A46*A54-A12*A23*A31*A46*A54-A13*A21*A32*A46* &
         A54+A11*A23*A32*A46*A54+A12*A21*A33*A46*A54-A11*A22*A33*A46*A54+A14*A23* &
         A32*A41*A56-A13*A24*A32*A41*A56-A14*A22*A33*A41*A56+A12*A24*A33*A41* &
         A56+A13*A22*A34*A41*A56-A12*A23*A34*A41*A56-A14*A23*A31*A42*A56+A13*A24* &
         A31*A42*A56+A14*A21*A33*A42*A56-A11*A24*A33*A42*A56-A13*A21*A34*A42* &
         A56+A11*A23*A34*A42*A56+A14*A22*A31*A43*A56-A12*A24*A31*A43*A56-A14*A21* &
         A32*A43*A56+A11*A24*A32*A43*A56+A12*A21*A34*A43*A56-A11*A22*A34*A43* &
         A56-A13*A22*A31*A44*A56+A12*A23*A31*A44*A56+A13*A21*A32*A44*A56-A11*A23* &
         A32*A44*A56-A12*A21*A33*A44*A56+A11*A22*A33*A44*A56)*A65+(A15*A24*A33* &
         A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+A13*A25*A34*A42*A51+A14* &
         A23*A35*A42*A51-A13*A24*A35*A42*A51-A15*A24*A32*A43*A51+A14*A25*A32*A43* &
         A51+A15*A22*A34*A43*A51-A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24* &
         A35*A43*A51+A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44* &
         A51+A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-A14*A23* &
         A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-A12*A24*A33*A45* &
         A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-A15*A24*A33*A41*A52+A14*A25* &
         A33*A41*A52+A15*A23*A34*A41*A52-A13*A25*A34*A41*A52-A14*A23*A35*A41* &
         A52+A13*A24*A35*A41*A52+A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21* &
         A34*A43*A52+A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43* &
         A52-A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52-A11*A25* &
         A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52+A14*A23*A31*A45* &
         A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+A11*A24*A33*A45*A52+A13*A21* &
         A34*A45*A52-A11*A23*A34*A45*A52+A15*A24*A32*A41*A53-A14*A25*A32*A41* &
         A53-A15*A22*A34*A41*A53+A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24* &
         A35*A41*A53-A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42* &
         A53-A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+A15*A22* &
         A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+A11*A25*A32*A44* &
         A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-A14*A22*A31*A45*A53+A12*A24* &
         A31*A45*A53+A14*A21*A32*A45*A53-A11*A24*A32*A45*A53-A12*A21*A34*A45* &
         A53+A11*A22*A34*A45*A53-A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22* &
         A33*A41*A54-A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41* &
         A54+A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54+A11*A25* &
         A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54-A15*A22*A31*A43* &
         A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-A11*A25*A32*A43*A54-A12*A21* &
         A35*A43*A54+A11*A22*A35*A43*A54+A13*A22*A31*A45*A54-A12*A23*A31*A45* &
         A54-A13*A21*A32*A45*A54+A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22* &
         A33*A45*A54+A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41* &
         A55+A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-A14*A23* &
         A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-A11*A24*A33*A42* &
         A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+A14*A22*A31*A43*A55-A12*A24* &
         A31*A43*A55-A14*A21*A32*A43*A55+A11*A24*A32*A43*A55+A12*A21*A34*A43* &
         A55-A11*A22*A34*A43*A55-A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21* &
         A32*A44*A55-A11*A23*A32*A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55)* &
         A66

      IF (ABS(DET) .LE. EPS) THEN
         A_INV = 0.0D0
         AA_INV = real(A_INV)
         OK_FLAG = .FALSE.
         ierr  = 1
         RETURN
      END IF

      COFACTOR(1,1) = A26*A35*A44*A53*A62-A25*A36*A44*A53*A62-A26*A34*A45*A53*A62+A24*A36*     &
      A45*A53*A62+A25*A34*A46*A53*A62-A24*A35*A46*A53*A62-A26*A35*A43*A54* &
      A62+A25*A36*A43*A54*A62+A26*A33*A45*A54*A62-A23*A36*A45*A54*A62-A25*A33* &
      A46*A54*A62+A23*A35*A46*A54*A62+A26*A34*A43*A55*A62-A24*A36*A43*A55* &
      A62-A26*A33*A44*A55*A62+A23*A36*A44*A55*A62+A24*A33*A46*A55*A62-A23*A34* &
      A46*A55*A62-A25*A34*A43*A56*A62+A24*A35*A43*A56*A62+A25*A33*A44*A56* &
      A62-A23*A35*A44*A56*A62-A24*A33*A45*A56*A62+A23*A34*A45*A56*A62-A26*A35* &
      A44*A52*A63+A25*A36*A44*A52*A63+A26*A34*A45*A52*A63-A24*A36*A45*A52* &
      A63-A25*A34*A46*A52*A63+A24*A35*A46*A52*A63+A26*A35*A42*A54*A63-A25*A36* &
      A42*A54*A63-A26*A32*A45*A54*A63+A22*A36*A45*A54*A63+A25*A32*A46*A54* &
      A63-A22*A35*A46*A54*A63-A26*A34*A42*A55*A63+A24*A36*A42*A55*A63+A26*A32* &
      A44*A55*A63-A22*A36*A44*A55*A63-A24*A32*A46*A55*A63+A22*A34*A46*A55* &
      A63+A25*A34*A42*A56*A63-A24*A35*A42*A56*A63-A25*A32*A44*A56*A63+A22*A35* &
      A44*A56*A63+A24*A32*A45*A56*A63-A22*A34*A45*A56*A63+A26*A35*A43*A52* &
      A64-A25*A36*A43*A52*A64-A26*A33*A45*A52*A64+A23*A36*A45*A52*A64+A25*A33* &
      A46*A52*A64-A23*A35*A46*A52*A64-A26*A35*A42*A53*A64+A25*A36*A42*A53* &
      A64+A26*A32*A45*A53*A64-A22*A36*A45*A53*A64-A25*A32*A46*A53*A64+A22*A35* &
      A46*A53*A64+A26*A33*A42*A55*A64-A23*A36*A42*A55*A64-A26*A32*A43*A55* &
      A64+A22*A36*A43*A55*A64+A23*A32*A46*A55*A64-A22*A33*A46*A55*A64-A25*A33* &
      A42*A56*A64+A23*A35*A42*A56*A64+A25*A32*A43*A56*A64-A22*A35*A43*A56* &
      A64-A23*A32*A45*A56*A64+A22*A33*A45*A56*A64-A26*A34*A43*A52*A65+A24*A36* &
      A43*A52*A65+A26*A33*A44*A52*A65-A23*A36*A44*A52*A65-A24*A33*A46*A52* &
      A65+A23*A34*A46*A52*A65+A26*A34*A42*A53*A65-A24*A36*A42*A53*A65-A26*A32* &
      A44*A53*A65+A22*A36*A44*A53*A65+A24*A32*A46*A53*A65-A22*A34*A46*A53* &
      A65-A26*A33*A42*A54*A65+A23*A36*A42*A54*A65+A26*A32*A43*A54*A65-A22*A36* &
      A43*A54*A65-A23*A32*A46*A54*A65+A22*A33*A46*A54*A65+A24*A33*A42*A56* &
      A65-A23*A34*A42*A56*A65-A24*A32*A43*A56*A65+A22*A34*A43*A56*A65+A23*A32* &
      A44*A56*A65-A22*A33*A44*A56*A65+A25*A34*A43*A52*A66-A24*A35*A43*A52* &
      A66-A25*A33*A44*A52*A66+A23*A35*A44*A52*A66+A24*A33*A45*A52*A66-A23*A34* &
      A45*A52*A66-A25*A34*A42*A53*A66+A24*A35*A42*A53*A66+A25*A32*A44*A53* &
      A66-A22*A35*A44*A53*A66-A24*A32*A45*A53*A66+A22*A34*A45*A53*A66+A25*A33* &
      A42*A54*A66-A23*A35*A42*A54*A66-A25*A32*A43*A54*A66+A22*A35*A43*A54* &
      A66+A23*A32*A45*A54*A66-A22*A33*A45*A54*A66-A24*A33*A42*A55*A66+A23*A34* &
      A42*A55*A66+A24*A32*A43*A55*A66-A22*A34*A43*A55*A66-A23*A32*A44*A55* &
      A66+A22*A33*A44*A55*A66

      COFACTOR(2,1) = -A16*A35*A44*A53*A62+A15*A36*A44*A53*A62+A16*A34* &
      A45*A53*A62-A14*A36*A45*A53*A62-A15*A34*A46*A53*A62+A14*A35*A46*A53* &
      A62+A16*A35*A43*A54*A62-A15*A36*A43*A54*A62-A16*A33*A45*A54*A62+A13*A36* &
      A45*A54*A62+A15*A33*A46*A54*A62-A13*A35*A46*A54*A62-A16*A34*A43*A55* &
      A62+A14*A36*A43*A55*A62+A16*A33*A44*A55*A62-A13*A36*A44*A55*A62-A14*A33* &
      A46*A55*A62+A13*A34*A46*A55*A62+A15*A34*A43*A56*A62-A14*A35*A43*A56* &
      A62-A15*A33*A44*A56*A62+A13*A35*A44*A56*A62+A14*A33*A45*A56*A62-A13*A34* &
      A45*A56*A62+A16*A35*A44*A52*A63-A15*A36*A44*A52*A63-A16*A34*A45*A52* &
      A63+A14*A36*A45*A52*A63+A15*A34*A46*A52*A63-A14*A35*A46*A52*A63-A16*A35* &
      A42*A54*A63+A15*A36*A42*A54*A63+A16*A32*A45*A54*A63-A12*A36*A45*A54* &
      A63-A15*A32*A46*A54*A63+A12*A35*A46*A54*A63+A16*A34*A42*A55*A63-A14*A36* &
      A42*A55*A63-A16*A32*A44*A55*A63+A12*A36*A44*A55*A63+A14*A32*A46*A55* &
      A63-A12*A34*A46*A55*A63-A15*A34*A42*A56*A63+A14*A35*A42*A56*A63+A15*A32* &
      A44*A56*A63-A12*A35*A44*A56*A63-A14*A32*A45*A56*A63+A12*A34*A45*A56* &
      A63-A16*A35*A43*A52*A64+A15*A36*A43*A52*A64+A16*A33*A45*A52*A64-A13*A36* &
      A45*A52*A64-A15*A33*A46*A52*A64+A13*A35*A46*A52*A64+A16*A35*A42*A53* &
      A64-A15*A36*A42*A53*A64-A16*A32*A45*A53*A64+A12*A36*A45*A53*A64+A15*A32* &
      A46*A53*A64-A12*A35*A46*A53*A64-A16*A33*A42*A55*A64+A13*A36*A42*A55* &
      A64+A16*A32*A43*A55*A64-A12*A36*A43*A55*A64-A13*A32*A46*A55*A64+A12*A33* &
      A46*A55*A64+A15*A33*A42*A56*A64-A13*A35*A42*A56*A64-A15*A32*A43*A56* &
      A64+A12*A35*A43*A56*A64+A13*A32*A45*A56*A64-A12*A33*A45*A56*A64+A16*A34* &
      A43*A52*A65-A14*A36*A43*A52*A65-A16*A33*A44*A52*A65+A13*A36*A44*A52* &
      A65+A14*A33*A46*A52*A65-A13*A34*A46*A52*A65-A16*A34*A42*A53*A65+A14*A36* &
      A42*A53*A65+A16*A32*A44*A53*A65-A12*A36*A44*A53*A65-A14*A32*A46*A53* &
      A65+A12*A34*A46*A53*A65+A16*A33*A42*A54*A65-A13*A36*A42*A54*A65-A16*A32* &
      A43*A54*A65+A12*A36*A43*A54*A65+A13*A32*A46*A54*A65-A12*A33*A46*A54* &
      A65-A14*A33*A42*A56*A65+A13*A34*A42*A56*A65+A14*A32*A43*A56*A65-A12*A34* &
      A43*A56*A65-A13*A32*A44*A56*A65+A12*A33*A44*A56*A65-A15*A34*A43*A52* &
      A66+A14*A35*A43*A52*A66+A15*A33*A44*A52*A66-A13*A35*A44*A52*A66-A14*A33* &
      A45*A52*A66+A13*A34*A45*A52*A66+A15*A34*A42*A53*A66-A14*A35*A42*A53* &
      A66-A15*A32*A44*A53*A66+A12*A35*A44*A53*A66+A14*A32*A45*A53*A66-A12*A34* &
      A45*A53*A66-A15*A33*A42*A54*A66+A13*A35*A42*A54*A66+A15*A32*A43*A54* &
      A66-A12*A35*A43*A54*A66-A13*A32*A45*A54*A66+A12*A33*A45*A54*A66+A14*A33* &
      A42*A55*A66-A13*A34*A42*A55*A66-A14*A32*A43*A55*A66+A12*A34*A43*A55* &
      A66+A13*A32*A44*A55*A66-A12*A33*A44*A55*A66

      COFACTOR(3,1) = A16*A25*A44*A53*A62-A15*A26* &
      A44*A53*A62-A16*A24*A45*A53*A62+A14*A26*A45*A53*A62+A15*A24*A46*A53* &
      A62-A14*A25*A46*A53*A62-A16*A25*A43*A54*A62+A15*A26*A43*A54*A62+A16*A23* &
      A45*A54*A62-A13*A26*A45*A54*A62-A15*A23*A46*A54*A62+A13*A25*A46*A54* &
      A62+A16*A24*A43*A55*A62-A14*A26*A43*A55*A62-A16*A23*A44*A55*A62+A13*A26* &
      A44*A55*A62+A14*A23*A46*A55*A62-A13*A24*A46*A55*A62-A15*A24*A43*A56* &
      A62+A14*A25*A43*A56*A62+A15*A23*A44*A56*A62-A13*A25*A44*A56*A62-A14*A23* &
      A45*A56*A62+A13*A24*A45*A56*A62-A16*A25*A44*A52*A63+A15*A26*A44*A52* &
      A63+A16*A24*A45*A52*A63-A14*A26*A45*A52*A63-A15*A24*A46*A52*A63+A14*A25* &
      A46*A52*A63+A16*A25*A42*A54*A63-A15*A26*A42*A54*A63-A16*A22*A45*A54* &
      A63+A12*A26*A45*A54*A63+A15*A22*A46*A54*A63-A12*A25*A46*A54*A63-A16*A24* &
      A42*A55*A63+A14*A26*A42*A55*A63+A16*A22*A44*A55*A63-A12*A26*A44*A55* &
      A63-A14*A22*A46*A55*A63+A12*A24*A46*A55*A63+A15*A24*A42*A56*A63-A14*A25* &
      A42*A56*A63-A15*A22*A44*A56*A63+A12*A25*A44*A56*A63+A14*A22*A45*A56* &
      A63-A12*A24*A45*A56*A63+A16*A25*A43*A52*A64-A15*A26*A43*A52*A64-A16*A23* &
      A45*A52*A64+A13*A26*A45*A52*A64+A15*A23*A46*A52*A64-A13*A25*A46*A52* &
      A64-A16*A25*A42*A53*A64+A15*A26*A42*A53*A64+A16*A22*A45*A53*A64-A12*A26* &
      A45*A53*A64-A15*A22*A46*A53*A64+A12*A25*A46*A53*A64+A16*A23*A42*A55* &
      A64-A13*A26*A42*A55*A64-A16*A22*A43*A55*A64+A12*A26*A43*A55*A64+A13*A22* &
      A46*A55*A64-A12*A23*A46*A55*A64-A15*A23*A42*A56*A64+A13*A25*A42*A56* &
      A64+A15*A22*A43*A56*A64-A12*A25*A43*A56*A64-A13*A22*A45*A56*A64+A12*A23* &
      A45*A56*A64-A16*A24*A43*A52*A65+A14*A26*A43*A52*A65+A16*A23*A44*A52* &
      A65-A13*A26*A44*A52*A65-A14*A23*A46*A52*A65+A13*A24*A46*A52*A65+A16*A24* &
      A42*A53*A65-A14*A26*A42*A53*A65-A16*A22*A44*A53*A65+A12*A26*A44*A53* &
      A65+A14*A22*A46*A53*A65-A12*A24*A46*A53*A65-A16*A23*A42*A54*A65+A13*A26* &
      A42*A54*A65+A16*A22*A43*A54*A65-A12*A26*A43*A54*A65-A13*A22*A46*A54* &
      A65+A12*A23*A46*A54*A65+A14*A23*A42*A56*A65-A13*A24*A42*A56*A65-A14*A22* &
      A43*A56*A65+A12*A24*A43*A56*A65+A13*A22*A44*A56*A65-A12*A23*A44*A56* &
      A65+A15*A24*A43*A52*A66-A14*A25*A43*A52*A66-A15*A23*A44*A52*A66+A13*A25* &
      A44*A52*A66+A14*A23*A45*A52*A66-A13*A24*A45*A52*A66-A15*A24*A42*A53* &
      A66+A14*A25*A42*A53*A66+A15*A22*A44*A53*A66-A12*A25*A44*A53*A66-A14*A22* &
      A45*A53*A66+A12*A24*A45*A53*A66+A15*A23*A42*A54*A66-A13*A25*A42*A54* &
      A66-A15*A22*A43*A54*A66+A12*A25*A43*A54*A66+A13*A22*A45*A54*A66-A12*A23* &
      A45*A54*A66-A14*A23*A42*A55*A66+A13*A24*A42*A55*A66+A14*A22*A43*A55* &
      A66-A12*A24*A43*A55*A66-A13*A22*A44*A55*A66+A12*A23*A44*A55*A66

      COFACTOR(4,1) = -A16*A25* &
      A34*A53*A62+A15*A26*A34*A53*A62+A16*A24*A35*A53*A62-A14*A26*A35*A53* &
      A62-A15*A24*A36*A53*A62+A14*A25*A36*A53*A62+A16*A25*A33*A54*A62-A15*A26* &
      A33*A54*A62-A16*A23*A35*A54*A62+A13*A26*A35*A54*A62+A15*A23*A36*A54* &
      A62-A13*A25*A36*A54*A62-A16*A24*A33*A55*A62+A14*A26*A33*A55*A62+A16*A23* &
      A34*A55*A62-A13*A26*A34*A55*A62-A14*A23*A36*A55*A62+A13*A24*A36*A55* &
      A62+A15*A24*A33*A56*A62-A14*A25*A33*A56*A62-A15*A23*A34*A56*A62+A13*A25* &
      A34*A56*A62+A14*A23*A35*A56*A62-A13*A24*A35*A56*A62+A16*A25*A34*A52* &
      A63-A15*A26*A34*A52*A63-A16*A24*A35*A52*A63+A14*A26*A35*A52*A63+A15*A24* &
      A36*A52*A63-A14*A25*A36*A52*A63-A16*A25*A32*A54*A63+A15*A26*A32*A54* &
      A63+A16*A22*A35*A54*A63-A12*A26*A35*A54*A63-A15*A22*A36*A54*A63+A12*A25* &
      A36*A54*A63+A16*A24*A32*A55*A63-A14*A26*A32*A55*A63-A16*A22*A34*A55* &
      A63+A12*A26*A34*A55*A63+A14*A22*A36*A55*A63-A12*A24*A36*A55*A63-A15*A24* &
      A32*A56*A63+A14*A25*A32*A56*A63+A15*A22*A34*A56*A63-A12*A25*A34*A56* &
      A63-A14*A22*A35*A56*A63+A12*A24*A35*A56*A63-A16*A25*A33*A52*A64+A15*A26* &
      A33*A52*A64+A16*A23*A35*A52*A64-A13*A26*A35*A52*A64-A15*A23*A36*A52* &
      A64+A13*A25*A36*A52*A64+A16*A25*A32*A53*A64-A15*A26*A32*A53*A64-A16*A22* &
      A35*A53*A64+A12*A26*A35*A53*A64+A15*A22*A36*A53*A64-A12*A25*A36*A53* &
      A64-A16*A23*A32*A55*A64+A13*A26*A32*A55*A64+A16*A22*A33*A55*A64-A12*A26* &
      A33*A55*A64-A13*A22*A36*A55*A64+A12*A23*A36*A55*A64+A15*A23*A32*A56* &
      A64-A13*A25*A32*A56*A64-A15*A22*A33*A56*A64+A12*A25*A33*A56*A64+A13*A22* &
      A35*A56*A64-A12*A23*A35*A56*A64+A16*A24*A33*A52*A65-A14*A26*A33*A52* &
      A65-A16*A23*A34*A52*A65+A13*A26*A34*A52*A65+A14*A23*A36*A52*A65-A13*A24* &
      A36*A52*A65-A16*A24*A32*A53*A65+A14*A26*A32*A53*A65+A16*A22*A34*A53* &
      A65-A12*A26*A34*A53*A65-A14*A22*A36*A53*A65+A12*A24*A36*A53*A65+A16*A23* &
      A32*A54*A65-A13*A26*A32*A54*A65-A16*A22*A33*A54*A65+A12*A26*A33*A54* &
      A65+A13*A22*A36*A54*A65-A12*A23*A36*A54*A65-A14*A23*A32*A56*A65+A13*A24* &
      A32*A56*A65+A14*A22*A33*A56*A65-A12*A24*A33*A56*A65-A13*A22*A34*A56* &
      A65+A12*A23*A34*A56*A65-A15*A24*A33*A52*A66+A14*A25*A33*A52*A66+A15*A23* &
      A34*A52*A66-A13*A25*A34*A52*A66-A14*A23*A35*A52*A66+A13*A24*A35*A52* &
      A66+A15*A24*A32*A53*A66-A14*A25*A32*A53*A66-A15*A22*A34*A53*A66+A12*A25* &
      A34*A53*A66+A14*A22*A35*A53*A66-A12*A24*A35*A53*A66-A15*A23*A32*A54* &
      A66+A13*A25*A32*A54*A66+A15*A22*A33*A54*A66-A12*A25*A33*A54*A66-A13*A22* &
      A35*A54*A66+A12*A23*A35*A54*A66+A14*A23*A32*A55*A66-A13*A24*A32*A55* &
      A66-A14*A22*A33*A55*A66+A12*A24*A33*A55*A66+A13*A22*A34*A55*A66-A12*A23* &
      A34*A55*A66

      COFACTOR(5,1) = A16*A25*A34*A43*A62-A15*A26*A34*A43*A62-A16*A24*A35*A43* &
      A62+A14*A26*A35*A43*A62+A15*A24*A36*A43*A62-A14*A25*A36*A43*A62-A16*A25* &
      A33*A44*A62+A15*A26*A33*A44*A62+A16*A23*A35*A44*A62-A13*A26*A35*A44* &
      A62-A15*A23*A36*A44*A62+A13*A25*A36*A44*A62+A16*A24*A33*A45*A62-A14*A26* &
      A33*A45*A62-A16*A23*A34*A45*A62+A13*A26*A34*A45*A62+A14*A23*A36*A45* &
      A62-A13*A24*A36*A45*A62-A15*A24*A33*A46*A62+A14*A25*A33*A46*A62+A15*A23* &
      A34*A46*A62-A13*A25*A34*A46*A62-A14*A23*A35*A46*A62+A13*A24*A35*A46* &
      A62-A16*A25*A34*A42*A63+A15*A26*A34*A42*A63+A16*A24*A35*A42*A63-A14*A26* &
      A35*A42*A63-A15*A24*A36*A42*A63+A14*A25*A36*A42*A63+A16*A25*A32*A44* &
      A63-A15*A26*A32*A44*A63-A16*A22*A35*A44*A63+A12*A26*A35*A44*A63+A15*A22* &
      A36*A44*A63-A12*A25*A36*A44*A63-A16*A24*A32*A45*A63+A14*A26*A32*A45* &
      A63+A16*A22*A34*A45*A63-A12*A26*A34*A45*A63-A14*A22*A36*A45*A63+A12*A24* &
      A36*A45*A63+A15*A24*A32*A46*A63-A14*A25*A32*A46*A63-A15*A22*A34*A46* &
      A63+A12*A25*A34*A46*A63+A14*A22*A35*A46*A63-A12*A24*A35*A46*A63+A16*A25* &
      A33*A42*A64-A15*A26*A33*A42*A64-A16*A23*A35*A42*A64+A13*A26*A35*A42* &
      A64+A15*A23*A36*A42*A64-A13*A25*A36*A42*A64-A16*A25*A32*A43*A64+A15*A26* &
      A32*A43*A64+A16*A22*A35*A43*A64-A12*A26*A35*A43*A64-A15*A22*A36*A43* &
      A64+A12*A25*A36*A43*A64+A16*A23*A32*A45*A64-A13*A26*A32*A45*A64-A16*A22* &
      A33*A45*A64+A12*A26*A33*A45*A64+A13*A22*A36*A45*A64-A12*A23*A36*A45* &
      A64-A15*A23*A32*A46*A64+A13*A25*A32*A46*A64+A15*A22*A33*A46*A64-A12*A25* &
      A33*A46*A64-A13*A22*A35*A46*A64+A12*A23*A35*A46*A64-A16*A24*A33*A42* &
      A65+A14*A26*A33*A42*A65+A16*A23*A34*A42*A65-A13*A26*A34*A42*A65-A14*A23* &
      A36*A42*A65+A13*A24*A36*A42*A65+A16*A24*A32*A43*A65-A14*A26*A32*A43* &
      A65-A16*A22*A34*A43*A65+A12*A26*A34*A43*A65+A14*A22*A36*A43*A65-A12*A24* &
      A36*A43*A65-A16*A23*A32*A44*A65+A13*A26*A32*A44*A65+A16*A22*A33*A44* &
      A65-A12*A26*A33*A44*A65-A13*A22*A36*A44*A65+A12*A23*A36*A44*A65+A14*A23* &
      A32*A46*A65-A13*A24*A32*A46*A65-A14*A22*A33*A46*A65+A12*A24*A33*A46* &
      A65+A13*A22*A34*A46*A65-A12*A23*A34*A46*A65+A15*A24*A33*A42*A66-A14*A25* &
      A33*A42*A66-A15*A23*A34*A42*A66+A13*A25*A34*A42*A66+A14*A23*A35*A42* &
      A66-A13*A24*A35*A42*A66-A15*A24*A32*A43*A66+A14*A25*A32*A43*A66+A15*A22* &
      A34*A43*A66-A12*A25*A34*A43*A66-A14*A22*A35*A43*A66+A12*A24*A35*A43* &
      A66+A15*A23*A32*A44*A66-A13*A25*A32*A44*A66-A15*A22*A33*A44*A66+A12*A25* &
      A33*A44*A66+A13*A22*A35*A44*A66-A12*A23*A35*A44*A66-A14*A23*A32*A45* &
      A66+A13*A24*A32*A45*A66+A14*A22*A33*A45*A66-A12*A24*A33*A45*A66-A13*A22* &
      A34*A45*A66+A12*A23*A34*A45*A66 

      COFACTOR(6,1) = -A16*A25*A34*A43*A52+A15*A26*A34*A43* &
      A52+A16*A24*A35*A43*A52-A14*A26*A35*A43*A52-A15*A24*A36*A43*A52+A14*A25* &
      A36*A43*A52+A16*A25*A33*A44*A52-A15*A26*A33*A44*A52-A16*A23*A35*A44* &
      A52+A13*A26*A35*A44*A52+A15*A23*A36*A44*A52-A13*A25*A36*A44*A52-A16*A24* &
      A33*A45*A52+A14*A26*A33*A45*A52+A16*A23*A34*A45*A52-A13*A26*A34*A45* &
      A52-A14*A23*A36*A45*A52+A13*A24*A36*A45*A52+A15*A24*A33*A46*A52-A14*A25* &
      A33*A46*A52-A15*A23*A34*A46*A52+A13*A25*A34*A46*A52+A14*A23*A35*A46* &
      A52-A13*A24*A35*A46*A52+A16*A25*A34*A42*A53-A15*A26*A34*A42*A53-A16*A24* &
      A35*A42*A53+A14*A26*A35*A42*A53+A15*A24*A36*A42*A53-A14*A25*A36*A42* &
      A53-A16*A25*A32*A44*A53+A15*A26*A32*A44*A53+A16*A22*A35*A44*A53-A12*A26* &
      A35*A44*A53-A15*A22*A36*A44*A53+A12*A25*A36*A44*A53+A16*A24*A32*A45* &
      A53-A14*A26*A32*A45*A53-A16*A22*A34*A45*A53+A12*A26*A34*A45*A53+A14*A22* &
      A36*A45*A53-A12*A24*A36*A45*A53-A15*A24*A32*A46*A53+A14*A25*A32*A46* &
      A53+A15*A22*A34*A46*A53-A12*A25*A34*A46*A53-A14*A22*A35*A46*A53+A12*A24* &
      A35*A46*A53-A16*A25*A33*A42*A54+A15*A26*A33*A42*A54+A16*A23*A35*A42* &
      A54-A13*A26*A35*A42*A54-A15*A23*A36*A42*A54+A13*A25*A36*A42*A54+A16*A25* &
      A32*A43*A54-A15*A26*A32*A43*A54-A16*A22*A35*A43*A54+A12*A26*A35*A43* &
      A54+A15*A22*A36*A43*A54-A12*A25*A36*A43*A54-A16*A23*A32*A45*A54+A13*A26* &
      A32*A45*A54+A16*A22*A33*A45*A54-A12*A26*A33*A45*A54-A13*A22*A36*A45* &
      A54+A12*A23*A36*A45*A54+A15*A23*A32*A46*A54-A13*A25*A32*A46*A54-A15*A22* &
      A33*A46*A54+A12*A25*A33*A46*A54+A13*A22*A35*A46*A54-A12*A23*A35*A46* &
      A54+A16*A24*A33*A42*A55-A14*A26*A33*A42*A55-A16*A23*A34*A42*A55+A13*A26* &
      A34*A42*A55+A14*A23*A36*A42*A55-A13*A24*A36*A42*A55-A16*A24*A32*A43* &
      A55+A14*A26*A32*A43*A55+A16*A22*A34*A43*A55-A12*A26*A34*A43*A55-A14*A22* &
      A36*A43*A55+A12*A24*A36*A43*A55+A16*A23*A32*A44*A55-A13*A26*A32*A44* &
      A55-A16*A22*A33*A44*A55+A12*A26*A33*A44*A55+A13*A22*A36*A44*A55-A12*A23* &
      A36*A44*A55-A14*A23*A32*A46*A55+A13*A24*A32*A46*A55+A14*A22*A33*A46* &
      A55-A12*A24*A33*A46*A55-A13*A22*A34*A46*A55+A12*A23*A34*A46*A55-A15*A24* &
      A33*A42*A56+A14*A25*A33*A42*A56+A15*A23*A34*A42*A56-A13*A25*A34*A42* &
      A56-A14*A23*A35*A42*A56+A13*A24*A35*A42*A56+A15*A24*A32*A43*A56-A14*A25* &
      A32*A43*A56-A15*A22*A34*A43*A56+A12*A25*A34*A43*A56+A14*A22*A35*A43* &
      A56-A12*A24*A35*A43*A56-A15*A23*A32*A44*A56+A13*A25*A32*A44*A56+A15*A22* &
      A33*A44*A56-A12*A25*A33*A44*A56-A13*A22*A35*A44*A56+A12*A23*A35*A44* &
      A56+A14*A23*A32*A45*A56-A13*A24*A32*A45*A56-A14*A22*A33*A45*A56+A12*A24* &
      A33*A45*A56+A13*A22*A34*A45*A56-A12*A23*A34*A45*A56 

      COFACTOR(1,2) = -A26*A35*A44*A53* &
      A61+A25*A36*A44*A53*A61+A26*A34*A45*A53*A61-A24*A36*A45*A53*A61-A25*A34* &
      A46*A53*A61+A24*A35*A46*A53*A61+A26*A35*A43*A54*A61-A25*A36*A43*A54* &
      A61-A26*A33*A45*A54*A61+A23*A36*A45*A54*A61+A25*A33*A46*A54*A61-A23*A35* &
      A46*A54*A61-A26*A34*A43*A55*A61+A24*A36*A43*A55*A61+A26*A33*A44*A55* &
      A61-A23*A36*A44*A55*A61-A24*A33*A46*A55*A61+A23*A34*A46*A55*A61+A25*A34* &
      A43*A56*A61-A24*A35*A43*A56*A61-A25*A33*A44*A56*A61+A23*A35*A44*A56* &
      A61+A24*A33*A45*A56*A61-A23*A34*A45*A56*A61+A26*A35*A44*A51*A63-A25*A36* &
      A44*A51*A63-A26*A34*A45*A51*A63+A24*A36*A45*A51*A63+A25*A34*A46*A51* &
      A63-A24*A35*A46*A51*A63-A26*A35*A41*A54*A63+A25*A36*A41*A54*A63+A26*A31* &
      A45*A54*A63-A21*A36*A45*A54*A63-A25*A31*A46*A54*A63+A21*A35*A46*A54* &
      A63+A26*A34*A41*A55*A63-A24*A36*A41*A55*A63-A26*A31*A44*A55*A63+A21*A36* &
      A44*A55*A63+A24*A31*A46*A55*A63-A21*A34*A46*A55*A63-A25*A34*A41*A56* &
      A63+A24*A35*A41*A56*A63+A25*A31*A44*A56*A63-A21*A35*A44*A56*A63-A24*A31* &
      A45*A56*A63+A21*A34*A45*A56*A63-A26*A35*A43*A51*A64+A25*A36*A43*A51* &
      A64+A26*A33*A45*A51*A64-A23*A36*A45*A51*A64-A25*A33*A46*A51*A64+A23*A35* &
      A46*A51*A64+A26*A35*A41*A53*A64-A25*A36*A41*A53*A64-A26*A31*A45*A53* &
      A64+A21*A36*A45*A53*A64+A25*A31*A46*A53*A64-A21*A35*A46*A53*A64-A26*A33* &
      A41*A55*A64+A23*A36*A41*A55*A64+A26*A31*A43*A55*A64-A21*A36*A43*A55* &
      A64-A23*A31*A46*A55*A64+A21*A33*A46*A55*A64+A25*A33*A41*A56*A64-A23*A35* &
      A41*A56*A64-A25*A31*A43*A56*A64+A21*A35*A43*A56*A64+A23*A31*A45*A56* &
      A64-A21*A33*A45*A56*A64+A26*A34*A43*A51*A65-A24*A36*A43*A51*A65-A26*A33* &
      A44*A51*A65+A23*A36*A44*A51*A65+A24*A33*A46*A51*A65-A23*A34*A46*A51* &
      A65-A26*A34*A41*A53*A65+A24*A36*A41*A53*A65+A26*A31*A44*A53*A65-A21*A36* &
      A44*A53*A65-A24*A31*A46*A53*A65+A21*A34*A46*A53*A65+A26*A33*A41*A54* &
      A65-A23*A36*A41*A54*A65-A26*A31*A43*A54*A65+A21*A36*A43*A54*A65+A23*A31* &
      A46*A54*A65-A21*A33*A46*A54*A65-A24*A33*A41*A56*A65+A23*A34*A41*A56* &
      A65+A24*A31*A43*A56*A65-A21*A34*A43*A56*A65-A23*A31*A44*A56*A65+A21*A33* &
      A44*A56*A65-A25*A34*A43*A51*A66+A24*A35*A43*A51*A66+A25*A33*A44*A51* &
      A66-A23*A35*A44*A51*A66-A24*A33*A45*A51*A66+A23*A34*A45*A51*A66+A25*A34* &
      A41*A53*A66-A24*A35*A41*A53*A66-A25*A31*A44*A53*A66+A21*A35*A44*A53* &
      A66+A24*A31*A45*A53*A66-A21*A34*A45*A53*A66-A25*A33*A41*A54*A66+A23*A35* &
      A41*A54*A66+A25*A31*A43*A54*A66-A21*A35*A43*A54*A66-A23*A31*A45*A54* &
      A66+A21*A33*A45*A54*A66+A24*A33*A41*A55*A66-A23*A34*A41*A55*A66-A24*A31* &
      A43*A55*A66+A21*A34*A43*A55*A66+A23*A31*A44*A55*A66-A21*A33*A44*A55* &
      A66 

      COFACTOR(2,2) = A16*A35*A44*A53*A61-A15*A36*A44*A53*A61-A16*A34*A45*A53*A61+A14*A36*     &
      A45*A53*A61+A15*A34*A46*A53*A61-A14*A35*A46*A53*A61-A16*A35*A43*A54* &
      A61+A15*A36*A43*A54*A61+A16*A33*A45*A54*A61-A13*A36*A45*A54*A61-A15*A33* &
      A46*A54*A61+A13*A35*A46*A54*A61+A16*A34*A43*A55*A61-A14*A36*A43*A55* &
      A61-A16*A33*A44*A55*A61+A13*A36*A44*A55*A61+A14*A33*A46*A55*A61-A13*A34* &
      A46*A55*A61-A15*A34*A43*A56*A61+A14*A35*A43*A56*A61+A15*A33*A44*A56* &
      A61-A13*A35*A44*A56*A61-A14*A33*A45*A56*A61+A13*A34*A45*A56*A61-A16*A35* &
      A44*A51*A63+A15*A36*A44*A51*A63+A16*A34*A45*A51*A63-A14*A36*A45*A51* &
      A63-A15*A34*A46*A51*A63+A14*A35*A46*A51*A63+A16*A35*A41*A54*A63-A15*A36* &
      A41*A54*A63-A16*A31*A45*A54*A63+A11*A36*A45*A54*A63+A15*A31*A46*A54* &
      A63-A11*A35*A46*A54*A63-A16*A34*A41*A55*A63+A14*A36*A41*A55*A63+A16*A31* &
      A44*A55*A63-A11*A36*A44*A55*A63-A14*A31*A46*A55*A63+A11*A34*A46*A55* &
      A63+A15*A34*A41*A56*A63-A14*A35*A41*A56*A63-A15*A31*A44*A56*A63+A11*A35* &
      A44*A56*A63+A14*A31*A45*A56*A63-A11*A34*A45*A56*A63+A16*A35*A43*A51* &
      A64-A15*A36*A43*A51*A64-A16*A33*A45*A51*A64+A13*A36*A45*A51*A64+A15*A33* &
      A46*A51*A64-A13*A35*A46*A51*A64-A16*A35*A41*A53*A64+A15*A36*A41*A53* &
      A64+A16*A31*A45*A53*A64-A11*A36*A45*A53*A64-A15*A31*A46*A53*A64+A11*A35* &
      A46*A53*A64+A16*A33*A41*A55*A64-A13*A36*A41*A55*A64-A16*A31*A43*A55* &
      A64+A11*A36*A43*A55*A64+A13*A31*A46*A55*A64-A11*A33*A46*A55*A64-A15*A33* &
      A41*A56*A64+A13*A35*A41*A56*A64+A15*A31*A43*A56*A64-A11*A35*A43*A56* &
      A64-A13*A31*A45*A56*A64+A11*A33*A45*A56*A64-A16*A34*A43*A51*A65+A14*A36* &
      A43*A51*A65+A16*A33*A44*A51*A65-A13*A36*A44*A51*A65-A14*A33*A46*A51* &
      A65+A13*A34*A46*A51*A65+A16*A34*A41*A53*A65-A14*A36*A41*A53*A65-A16*A31* &
      A44*A53*A65+A11*A36*A44*A53*A65+A14*A31*A46*A53*A65-A11*A34*A46*A53* &
      A65-A16*A33*A41*A54*A65+A13*A36*A41*A54*A65+A16*A31*A43*A54*A65-A11*A36* &
      A43*A54*A65-A13*A31*A46*A54*A65+A11*A33*A46*A54*A65+A14*A33*A41*A56* &
      A65-A13*A34*A41*A56*A65-A14*A31*A43*A56*A65+A11*A34*A43*A56*A65+A13*A31* &
      A44*A56*A65-A11*A33*A44*A56*A65+A15*A34*A43*A51*A66-A14*A35*A43*A51* &
      A66-A15*A33*A44*A51*A66+A13*A35*A44*A51*A66+A14*A33*A45*A51*A66-A13*A34* &
      A45*A51*A66-A15*A34*A41*A53*A66+A14*A35*A41*A53*A66+A15*A31*A44*A53* &
      A66-A11*A35*A44*A53*A66-A14*A31*A45*A53*A66+A11*A34*A45*A53*A66+A15*A33* &
      A41*A54*A66-A13*A35*A41*A54*A66-A15*A31*A43*A54*A66+A11*A35*A43*A54* &
      A66+A13*A31*A45*A54*A66-A11*A33*A45*A54*A66-A14*A33*A41*A55*A66+A13*A34* &
      A41*A55*A66+A14*A31*A43*A55*A66-A11*A34*A43*A55*A66-A13*A31*A44*A55* &
      A66+A11*A33*A44*A55*A66

      COFACTOR(3,2) = -A16*A25*A44*A53*A61+A15*A26*A44*A53*A61+A16*A24* &
      A45*A53*A61-A14*A26*A45*A53*A61-A15*A24*A46*A53*A61+A14*A25*A46*A53* &
      A61+A16*A25*A43*A54*A61-A15*A26*A43*A54*A61-A16*A23*A45*A54*A61+A13*A26* &
      A45*A54*A61+A15*A23*A46*A54*A61-A13*A25*A46*A54*A61-A16*A24*A43*A55* &
      A61+A14*A26*A43*A55*A61+A16*A23*A44*A55*A61-A13*A26*A44*A55*A61-A14*A23* &
      A46*A55*A61+A13*A24*A46*A55*A61+A15*A24*A43*A56*A61-A14*A25*A43*A56* &
      A61-A15*A23*A44*A56*A61+A13*A25*A44*A56*A61+A14*A23*A45*A56*A61-A13*A24* &
      A45*A56*A61+A16*A25*A44*A51*A63-A15*A26*A44*A51*A63-A16*A24*A45*A51* &
      A63+A14*A26*A45*A51*A63+A15*A24*A46*A51*A63-A14*A25*A46*A51*A63-A16*A25* &
      A41*A54*A63+A15*A26*A41*A54*A63+A16*A21*A45*A54*A63-A11*A26*A45*A54* &
      A63-A15*A21*A46*A54*A63+A11*A25*A46*A54*A63+A16*A24*A41*A55*A63-A14*A26* &
      A41*A55*A63-A16*A21*A44*A55*A63+A11*A26*A44*A55*A63+A14*A21*A46*A55* &
      A63-A11*A24*A46*A55*A63-A15*A24*A41*A56*A63+A14*A25*A41*A56*A63+A15*A21* &
      A44*A56*A63-A11*A25*A44*A56*A63-A14*A21*A45*A56*A63+A11*A24*A45*A56* &
      A63-A16*A25*A43*A51*A64+A15*A26*A43*A51*A64+A16*A23*A45*A51*A64-A13*A26* &
      A45*A51*A64-A15*A23*A46*A51*A64+A13*A25*A46*A51*A64+A16*A25*A41*A53* &
      A64-A15*A26*A41*A53*A64-A16*A21*A45*A53*A64+A11*A26*A45*A53*A64+A15*A21* &
      A46*A53*A64-A11*A25*A46*A53*A64-A16*A23*A41*A55*A64+A13*A26*A41*A55* &
      A64+A16*A21*A43*A55*A64-A11*A26*A43*A55*A64-A13*A21*A46*A55*A64+A11*A23* &
      A46*A55*A64+A15*A23*A41*A56*A64-A13*A25*A41*A56*A64-A15*A21*A43*A56* &
      A64+A11*A25*A43*A56*A64+A13*A21*A45*A56*A64-A11*A23*A45*A56*A64+A16*A24* &
      A43*A51*A65-A14*A26*A43*A51*A65-A16*A23*A44*A51*A65+A13*A26*A44*A51* &
      A65+A14*A23*A46*A51*A65-A13*A24*A46*A51*A65-A16*A24*A41*A53*A65+A14*A26* &
      A41*A53*A65+A16*A21*A44*A53*A65-A11*A26*A44*A53*A65-A14*A21*A46*A53* &
      A65+A11*A24*A46*A53*A65+A16*A23*A41*A54*A65-A13*A26*A41*A54*A65-A16*A21* &
      A43*A54*A65+A11*A26*A43*A54*A65+A13*A21*A46*A54*A65-A11*A23*A46*A54* &
      A65-A14*A23*A41*A56*A65+A13*A24*A41*A56*A65+A14*A21*A43*A56*A65-A11*A24* &
      A43*A56*A65-A13*A21*A44*A56*A65+A11*A23*A44*A56*A65-A15*A24*A43*A51* &
      A66+A14*A25*A43*A51*A66+A15*A23*A44*A51*A66-A13*A25*A44*A51*A66-A14*A23* &
      A45*A51*A66+A13*A24*A45*A51*A66+A15*A24*A41*A53*A66-A14*A25*A41*A53* &
      A66-A15*A21*A44*A53*A66+A11*A25*A44*A53*A66+A14*A21*A45*A53*A66-A11*A24* &
      A45*A53*A66-A15*A23*A41*A54*A66+A13*A25*A41*A54*A66+A15*A21*A43*A54* &
      A66-A11*A25*A43*A54*A66-A13*A21*A45*A54*A66+A11*A23*A45*A54*A66+A14*A23* &
      A41*A55*A66-A13*A24*A41*A55*A66-A14*A21*A43*A55*A66+A11*A24*A43*A55* &
      A66+A13*A21*A44*A55*A66-A11*A23*A44*A55*A66 

      COFACTOR(4,2) = A16*A25*A34*A53*A61-A15*A26* &
      A34*A53*A61-A16*A24*A35*A53*A61+A14*A26*A35*A53*A61+A15*A24*A36*A53* &
      A61-A14*A25*A36*A53*A61-A16*A25*A33*A54*A61+A15*A26*A33*A54*A61+A16*A23* &
      A35*A54*A61-A13*A26*A35*A54*A61-A15*A23*A36*A54*A61+A13*A25*A36*A54* &
      A61+A16*A24*A33*A55*A61-A14*A26*A33*A55*A61-A16*A23*A34*A55*A61+A13*A26* &
      A34*A55*A61+A14*A23*A36*A55*A61-A13*A24*A36*A55*A61-A15*A24*A33*A56* &
      A61+A14*A25*A33*A56*A61+A15*A23*A34*A56*A61-A13*A25*A34*A56*A61-A14*A23* &
      A35*A56*A61+A13*A24*A35*A56*A61-A16*A25*A34*A51*A63+A15*A26*A34*A51* &
      A63+A16*A24*A35*A51*A63-A14*A26*A35*A51*A63-A15*A24*A36*A51*A63+A14*A25* &
      A36*A51*A63+A16*A25*A31*A54*A63-A15*A26*A31*A54*A63-A16*A21*A35*A54* &
      A63+A11*A26*A35*A54*A63+A15*A21*A36*A54*A63-A11*A25*A36*A54*A63-A16*A24* &
      A31*A55*A63+A14*A26*A31*A55*A63+A16*A21*A34*A55*A63-A11*A26*A34*A55* &
      A63-A14*A21*A36*A55*A63+A11*A24*A36*A55*A63+A15*A24*A31*A56*A63-A14*A25* &
      A31*A56*A63-A15*A21*A34*A56*A63+A11*A25*A34*A56*A63+A14*A21*A35*A56* &
      A63-A11*A24*A35*A56*A63+A16*A25*A33*A51*A64-A15*A26*A33*A51*A64-A16*A23* &
      A35*A51*A64+A13*A26*A35*A51*A64+A15*A23*A36*A51*A64-A13*A25*A36*A51* &
      A64-A16*A25*A31*A53*A64+A15*A26*A31*A53*A64+A16*A21*A35*A53*A64-A11*A26* &
      A35*A53*A64-A15*A21*A36*A53*A64+A11*A25*A36*A53*A64+A16*A23*A31*A55* &
      A64-A13*A26*A31*A55*A64-A16*A21*A33*A55*A64+A11*A26*A33*A55*A64+A13*A21* &
      A36*A55*A64-A11*A23*A36*A55*A64-A15*A23*A31*A56*A64+A13*A25*A31*A56* &
      A64+A15*A21*A33*A56*A64-A11*A25*A33*A56*A64-A13*A21*A35*A56*A64+A11*A23* &
      A35*A56*A64-A16*A24*A33*A51*A65+A14*A26*A33*A51*A65+A16*A23*A34*A51* &
      A65-A13*A26*A34*A51*A65-A14*A23*A36*A51*A65+A13*A24*A36*A51*A65+A16*A24* &
      A31*A53*A65-A14*A26*A31*A53*A65-A16*A21*A34*A53*A65+A11*A26*A34*A53* &
      A65+A14*A21*A36*A53*A65-A11*A24*A36*A53*A65-A16*A23*A31*A54*A65+A13*A26* &
      A31*A54*A65+A16*A21*A33*A54*A65-A11*A26*A33*A54*A65-A13*A21*A36*A54* &
      A65+A11*A23*A36*A54*A65+A14*A23*A31*A56*A65-A13*A24*A31*A56*A65-A14*A21* &
      A33*A56*A65+A11*A24*A33*A56*A65+A13*A21*A34*A56*A65-A11*A23*A34*A56* &
      A65+A15*A24*A33*A51*A66-A14*A25*A33*A51*A66-A15*A23*A34*A51*A66+A13*A25* &
      A34*A51*A66+A14*A23*A35*A51*A66-A13*A24*A35*A51*A66-A15*A24*A31*A53* &
      A66+A14*A25*A31*A53*A66+A15*A21*A34*A53*A66-A11*A25*A34*A53*A66-A14*A21* &
      A35*A53*A66+A11*A24*A35*A53*A66+A15*A23*A31*A54*A66-A13*A25*A31*A54* &
      A66-A15*A21*A33*A54*A66+A11*A25*A33*A54*A66+A13*A21*A35*A54*A66-A11*A23* &
      A35*A54*A66-A14*A23*A31*A55*A66+A13*A24*A31*A55*A66+A14*A21*A33*A55* &
      A66-A11*A24*A33*A55*A66-A13*A21*A34*A55*A66+A11*A23*A34*A55*A66

      COFACTOR(5,2) = -A16*A25* &
      A34*A43*A61+A15*A26*A34*A43*A61+A16*A24*A35*A43*A61-A14*A26*A35*A43* &
      A61-A15*A24*A36*A43*A61+A14*A25*A36*A43*A61+A16*A25*A33*A44*A61-A15*A26* &
      A33*A44*A61-A16*A23*A35*A44*A61+A13*A26*A35*A44*A61+A15*A23*A36*A44* &
      A61-A13*A25*A36*A44*A61-A16*A24*A33*A45*A61+A14*A26*A33*A45*A61+A16*A23* &
      A34*A45*A61-A13*A26*A34*A45*A61-A14*A23*A36*A45*A61+A13*A24*A36*A45* &
      A61+A15*A24*A33*A46*A61-A14*A25*A33*A46*A61-A15*A23*A34*A46*A61+A13*A25* &
      A34*A46*A61+A14*A23*A35*A46*A61-A13*A24*A35*A46*A61+A16*A25*A34*A41* &
      A63-A15*A26*A34*A41*A63-A16*A24*A35*A41*A63+A14*A26*A35*A41*A63+A15*A24* &
      A36*A41*A63-A14*A25*A36*A41*A63-A16*A25*A31*A44*A63+A15*A26*A31*A44* &
      A63+A16*A21*A35*A44*A63-A11*A26*A35*A44*A63-A15*A21*A36*A44*A63+A11*A25* &
      A36*A44*A63+A16*A24*A31*A45*A63-A14*A26*A31*A45*A63-A16*A21*A34*A45* &
      A63+A11*A26*A34*A45*A63+A14*A21*A36*A45*A63-A11*A24*A36*A45*A63-A15*A24* &
      A31*A46*A63+A14*A25*A31*A46*A63+A15*A21*A34*A46*A63-A11*A25*A34*A46* &
      A63-A14*A21*A35*A46*A63+A11*A24*A35*A46*A63-A16*A25*A33*A41*A64+A15*A26* &
      A33*A41*A64+A16*A23*A35*A41*A64-A13*A26*A35*A41*A64-A15*A23*A36*A41* &
      A64+A13*A25*A36*A41*A64+A16*A25*A31*A43*A64-A15*A26*A31*A43*A64-A16*A21* &
      A35*A43*A64+A11*A26*A35*A43*A64+A15*A21*A36*A43*A64-A11*A25*A36*A43* &
      A64-A16*A23*A31*A45*A64+A13*A26*A31*A45*A64+A16*A21*A33*A45*A64-A11*A26* &
      A33*A45*A64-A13*A21*A36*A45*A64+A11*A23*A36*A45*A64+A15*A23*A31*A46* &
      A64-A13*A25*A31*A46*A64-A15*A21*A33*A46*A64+A11*A25*A33*A46*A64+A13*A21* &
      A35*A46*A64-A11*A23*A35*A46*A64+A16*A24*A33*A41*A65-A14*A26*A33*A41* &
      A65-A16*A23*A34*A41*A65+A13*A26*A34*A41*A65+A14*A23*A36*A41*A65-A13*A24* &
      A36*A41*A65-A16*A24*A31*A43*A65+A14*A26*A31*A43*A65+A16*A21*A34*A43* &
      A65-A11*A26*A34*A43*A65-A14*A21*A36*A43*A65+A11*A24*A36*A43*A65+A16*A23* &
      A31*A44*A65-A13*A26*A31*A44*A65-A16*A21*A33*A44*A65+A11*A26*A33*A44* &
      A65+A13*A21*A36*A44*A65-A11*A23*A36*A44*A65-A14*A23*A31*A46*A65+A13*A24* &
      A31*A46*A65+A14*A21*A33*A46*A65-A11*A24*A33*A46*A65-A13*A21*A34*A46* &
      A65+A11*A23*A34*A46*A65-A15*A24*A33*A41*A66+A14*A25*A33*A41*A66+A15*A23* &
      A34*A41*A66-A13*A25*A34*A41*A66-A14*A23*A35*A41*A66+A13*A24*A35*A41* &
      A66+A15*A24*A31*A43*A66-A14*A25*A31*A43*A66-A15*A21*A34*A43*A66+A11*A25* &
      A34*A43*A66+A14*A21*A35*A43*A66-A11*A24*A35*A43*A66-A15*A23*A31*A44* &
      A66+A13*A25*A31*A44*A66+A15*A21*A33*A44*A66-A11*A25*A33*A44*A66-A13*A21* &
      A35*A44*A66+A11*A23*A35*A44*A66+A14*A23*A31*A45*A66-A13*A24*A31*A45* &
      A66-A14*A21*A33*A45*A66+A11*A24*A33*A45*A66+A13*A21*A34*A45*A66-A11*A23* &
      A34*A45*A66

      COFACTOR(6,2) = A16*A25*A34*A43*A51-A15*A26*A34*A43*A51-A16*A24*A35*A43* &
      A51+A14*A26*A35*A43*A51+A15*A24*A36*A43*A51-A14*A25*A36*A43*A51-A16*A25* &
      A33*A44*A51+A15*A26*A33*A44*A51+A16*A23*A35*A44*A51-A13*A26*A35*A44* &
      A51-A15*A23*A36*A44*A51+A13*A25*A36*A44*A51+A16*A24*A33*A45*A51-A14*A26* &
      A33*A45*A51-A16*A23*A34*A45*A51+A13*A26*A34*A45*A51+A14*A23*A36*A45* &
      A51-A13*A24*A36*A45*A51-A15*A24*A33*A46*A51+A14*A25*A33*A46*A51+A15*A23* &
      A34*A46*A51-A13*A25*A34*A46*A51-A14*A23*A35*A46*A51+A13*A24*A35*A46* &
      A51-A16*A25*A34*A41*A53+A15*A26*A34*A41*A53+A16*A24*A35*A41*A53-A14*A26* &
      A35*A41*A53-A15*A24*A36*A41*A53+A14*A25*A36*A41*A53+A16*A25*A31*A44* &
      A53-A15*A26*A31*A44*A53-A16*A21*A35*A44*A53+A11*A26*A35*A44*A53+A15*A21* &
      A36*A44*A53-A11*A25*A36*A44*A53-A16*A24*A31*A45*A53+A14*A26*A31*A45* &
      A53+A16*A21*A34*A45*A53-A11*A26*A34*A45*A53-A14*A21*A36*A45*A53+A11*A24* &
      A36*A45*A53+A15*A24*A31*A46*A53-A14*A25*A31*A46*A53-A15*A21*A34*A46* &
      A53+A11*A25*A34*A46*A53+A14*A21*A35*A46*A53-A11*A24*A35*A46*A53+A16*A25* &
      A33*A41*A54-A15*A26*A33*A41*A54-A16*A23*A35*A41*A54+A13*A26*A35*A41* &
      A54+A15*A23*A36*A41*A54-A13*A25*A36*A41*A54-A16*A25*A31*A43*A54+A15*A26* &
      A31*A43*A54+A16*A21*A35*A43*A54-A11*A26*A35*A43*A54-A15*A21*A36*A43* &
      A54+A11*A25*A36*A43*A54+A16*A23*A31*A45*A54-A13*A26*A31*A45*A54-A16*A21* &
      A33*A45*A54+A11*A26*A33*A45*A54+A13*A21*A36*A45*A54-A11*A23*A36*A45* &
      A54-A15*A23*A31*A46*A54+A13*A25*A31*A46*A54+A15*A21*A33*A46*A54-A11*A25* &
      A33*A46*A54-A13*A21*A35*A46*A54+A11*A23*A35*A46*A54-A16*A24*A33*A41* &
      A55+A14*A26*A33*A41*A55+A16*A23*A34*A41*A55-A13*A26*A34*A41*A55-A14*A23* &
      A36*A41*A55+A13*A24*A36*A41*A55+A16*A24*A31*A43*A55-A14*A26*A31*A43* &
      A55-A16*A21*A34*A43*A55+A11*A26*A34*A43*A55+A14*A21*A36*A43*A55-A11*A24* &
      A36*A43*A55-A16*A23*A31*A44*A55+A13*A26*A31*A44*A55+A16*A21*A33*A44* &
      A55-A11*A26*A33*A44*A55-A13*A21*A36*A44*A55+A11*A23*A36*A44*A55+A14*A23* &
      A31*A46*A55-A13*A24*A31*A46*A55-A14*A21*A33*A46*A55+A11*A24*A33*A46* &
      A55+A13*A21*A34*A46*A55-A11*A23*A34*A46*A55+A15*A24*A33*A41*A56-A14*A25* &
      A33*A41*A56-A15*A23*A34*A41*A56+A13*A25*A34*A41*A56+A14*A23*A35*A41* &
      A56-A13*A24*A35*A41*A56-A15*A24*A31*A43*A56+A14*A25*A31*A43*A56+A15*A21* &
      A34*A43*A56-A11*A25*A34*A43*A56-A14*A21*A35*A43*A56+A11*A24*A35*A43* &
      A56+A15*A23*A31*A44*A56-A13*A25*A31*A44*A56-A15*A21*A33*A44*A56+A11*A25* &
      A33*A44*A56+A13*A21*A35*A44*A56-A11*A23*A35*A44*A56-A14*A23*A31*A45* &
      A56+A13*A24*A31*A45*A56+A14*A21*A33*A45*A56-A11*A24*A33*A45*A56-A13*A21* &
      A34*A45*A56+A11*A23*A34*A45*A56

      COFACTOR(1,3) = A26*A35*A44*A52*A61-A25*A36*A44*A52* &
      A61-A26*A34*A45*A52*A61+A24*A36*A45*A52*A61+A25*A34*A46*A52*A61-A24*A35* &
      A46*A52*A61-A26*A35*A42*A54*A61+A25*A36*A42*A54*A61+A26*A32*A45*A54* &
      A61-A22*A36*A45*A54*A61-A25*A32*A46*A54*A61+A22*A35*A46*A54*A61+A26*A34* &
      A42*A55*A61-A24*A36*A42*A55*A61-A26*A32*A44*A55*A61+A22*A36*A44*A55* &
      A61+A24*A32*A46*A55*A61-A22*A34*A46*A55*A61-A25*A34*A42*A56*A61+A24*A35* &
      A42*A56*A61+A25*A32*A44*A56*A61-A22*A35*A44*A56*A61-A24*A32*A45*A56* &
      A61+A22*A34*A45*A56*A61-A26*A35*A44*A51*A62+A25*A36*A44*A51*A62+A26*A34* &
      A45*A51*A62-A24*A36*A45*A51*A62-A25*A34*A46*A51*A62+A24*A35*A46*A51* &
      A62+A26*A35*A41*A54*A62-A25*A36*A41*A54*A62-A26*A31*A45*A54*A62+A21*A36* &
      A45*A54*A62+A25*A31*A46*A54*A62-A21*A35*A46*A54*A62-A26*A34*A41*A55* &
      A62+A24*A36*A41*A55*A62+A26*A31*A44*A55*A62-A21*A36*A44*A55*A62-A24*A31* &
      A46*A55*A62+A21*A34*A46*A55*A62+A25*A34*A41*A56*A62-A24*A35*A41*A56* &
      A62-A25*A31*A44*A56*A62+A21*A35*A44*A56*A62+A24*A31*A45*A56*A62-A21*A34* &
      A45*A56*A62+A26*A35*A42*A51*A64-A25*A36*A42*A51*A64-A26*A32*A45*A51* &
      A64+A22*A36*A45*A51*A64+A25*A32*A46*A51*A64-A22*A35*A46*A51*A64-A26*A35* &
      A41*A52*A64+A25*A36*A41*A52*A64+A26*A31*A45*A52*A64-A21*A36*A45*A52* &
      A64-A25*A31*A46*A52*A64+A21*A35*A46*A52*A64+A26*A32*A41*A55*A64-A22*A36* &
      A41*A55*A64-A26*A31*A42*A55*A64+A21*A36*A42*A55*A64+A22*A31*A46*A55* &
      A64-A21*A32*A46*A55*A64-A25*A32*A41*A56*A64+A22*A35*A41*A56*A64+A25*A31* &
      A42*A56*A64-A21*A35*A42*A56*A64-A22*A31*A45*A56*A64+A21*A32*A45*A56* &
      A64-A26*A34*A42*A51*A65+A24*A36*A42*A51*A65+A26*A32*A44*A51*A65-A22*A36* &
      A44*A51*A65-A24*A32*A46*A51*A65+A22*A34*A46*A51*A65+A26*A34*A41*A52* &
      A65-A24*A36*A41*A52*A65-A26*A31*A44*A52*A65+A21*A36*A44*A52*A65+A24*A31* &
      A46*A52*A65-A21*A34*A46*A52*A65-A26*A32*A41*A54*A65+A22*A36*A41*A54* &
      A65+A26*A31*A42*A54*A65-A21*A36*A42*A54*A65-A22*A31*A46*A54*A65+A21*A32* &
      A46*A54*A65+A24*A32*A41*A56*A65-A22*A34*A41*A56*A65-A24*A31*A42*A56* &
      A65+A21*A34*A42*A56*A65+A22*A31*A44*A56*A65-A21*A32*A44*A56*A65+A25*A34* &
      A42*A51*A66-A24*A35*A42*A51*A66-A25*A32*A44*A51*A66+A22*A35*A44*A51* &
      A66+A24*A32*A45*A51*A66-A22*A34*A45*A51*A66-A25*A34*A41*A52*A66+A24*A35* &
      A41*A52*A66+A25*A31*A44*A52*A66-A21*A35*A44*A52*A66-A24*A31*A45*A52* &
      A66+A21*A34*A45*A52*A66+A25*A32*A41*A54*A66-A22*A35*A41*A54*A66-A25*A31* &
      A42*A54*A66+A21*A35*A42*A54*A66+A22*A31*A45*A54*A66-A21*A32*A45*A54* &
      A66-A24*A32*A41*A55*A66+A22*A34*A41*A55*A66+A24*A31*A42*A55*A66-A21*A34* &
      A42*A55*A66-A22*A31*A44*A55*A66+A21*A32*A44*A55*A66 

      COFACTOR(2,3) = -A16*A35*A44*A52* &
      A61+A15*A36*A44*A52*A61+A16*A34*A45*A52*A61-A14*A36*A45*A52*A61-A15*A34* &
      A46*A52*A61+A14*A35*A46*A52*A61+A16*A35*A42*A54*A61-A15*A36*A42*A54* &
      A61-A16*A32*A45*A54*A61+A12*A36*A45*A54*A61+A15*A32*A46*A54*A61-A12*A35* &
      A46*A54*A61-A16*A34*A42*A55*A61+A14*A36*A42*A55*A61+A16*A32*A44*A55* &
      A61-A12*A36*A44*A55*A61-A14*A32*A46*A55*A61+A12*A34*A46*A55*A61+A15*A34* &
      A42*A56*A61-A14*A35*A42*A56*A61-A15*A32*A44*A56*A61+A12*A35*A44*A56* &
      A61+A14*A32*A45*A56*A61-A12*A34*A45*A56*A61+A16*A35*A44*A51*A62-A15*A36* &
      A44*A51*A62-A16*A34*A45*A51*A62+A14*A36*A45*A51*A62+A15*A34*A46*A51* &
      A62-A14*A35*A46*A51*A62-A16*A35*A41*A54*A62+A15*A36*A41*A54*A62+A16*A31* &
      A45*A54*A62-A11*A36*A45*A54*A62-A15*A31*A46*A54*A62+A11*A35*A46*A54* &
      A62+A16*A34*A41*A55*A62-A14*A36*A41*A55*A62-A16*A31*A44*A55*A62+A11*A36* &
      A44*A55*A62+A14*A31*A46*A55*A62-A11*A34*A46*A55*A62-A15*A34*A41*A56* &
      A62+A14*A35*A41*A56*A62+A15*A31*A44*A56*A62-A11*A35*A44*A56*A62-A14*A31* &
      A45*A56*A62+A11*A34*A45*A56*A62-A16*A35*A42*A51*A64+A15*A36*A42*A51* &
      A64+A16*A32*A45*A51*A64-A12*A36*A45*A51*A64-A15*A32*A46*A51*A64+A12*A35* &
      A46*A51*A64+A16*A35*A41*A52*A64-A15*A36*A41*A52*A64-A16*A31*A45*A52* &
      A64+A11*A36*A45*A52*A64+A15*A31*A46*A52*A64-A11*A35*A46*A52*A64-A16*A32* &
      A41*A55*A64+A12*A36*A41*A55*A64+A16*A31*A42*A55*A64-A11*A36*A42*A55* &
      A64-A12*A31*A46*A55*A64+A11*A32*A46*A55*A64+A15*A32*A41*A56*A64-A12*A35* &
      A41*A56*A64-A15*A31*A42*A56*A64+A11*A35*A42*A56*A64+A12*A31*A45*A56* &
      A64-A11*A32*A45*A56*A64+A16*A34*A42*A51*A65-A14*A36*A42*A51*A65-A16*A32* &
      A44*A51*A65+A12*A36*A44*A51*A65+A14*A32*A46*A51*A65-A12*A34*A46*A51* &
      A65-A16*A34*A41*A52*A65+A14*A36*A41*A52*A65+A16*A31*A44*A52*A65-A11*A36* &
      A44*A52*A65-A14*A31*A46*A52*A65+A11*A34*A46*A52*A65+A16*A32*A41*A54* &
      A65-A12*A36*A41*A54*A65-A16*A31*A42*A54*A65+A11*A36*A42*A54*A65+A12*A31* &
      A46*A54*A65-A11*A32*A46*A54*A65-A14*A32*A41*A56*A65+A12*A34*A41*A56* &
      A65+A14*A31*A42*A56*A65-A11*A34*A42*A56*A65-A12*A31*A44*A56*A65+A11*A32* &
      A44*A56*A65-A15*A34*A42*A51*A66+A14*A35*A42*A51*A66+A15*A32*A44*A51* &
      A66-A12*A35*A44*A51*A66-A14*A32*A45*A51*A66+A12*A34*A45*A51*A66+A15*A34* &
      A41*A52*A66-A14*A35*A41*A52*A66-A15*A31*A44*A52*A66+A11*A35*A44*A52* &
      A66+A14*A31*A45*A52*A66-A11*A34*A45*A52*A66-A15*A32*A41*A54*A66+A12*A35* &
      A41*A54*A66+A15*A31*A42*A54*A66-A11*A35*A42*A54*A66-A12*A31*A45*A54* &
      A66+A11*A32*A45*A54*A66+A14*A32*A41*A55*A66-A12*A34*A41*A55*A66-A14*A31* &
      A42*A55*A66+A11*A34*A42*A55*A66+A12*A31*A44*A55*A66-A11*A32*A44*A55* &
      A66 

      COFACTOR(3,3) = A16*A25*A44*A52*A61-A15*A26*A44*A52*A61-A16*A24*A45*A52*A61+A14*A26*     &
      A45*A52*A61+A15*A24*A46*A52*A61-A14*A25*A46*A52*A61-A16*A25*A42*A54* &
      A61+A15*A26*A42*A54*A61+A16*A22*A45*A54*A61-A12*A26*A45*A54*A61-A15*A22* &
      A46*A54*A61+A12*A25*A46*A54*A61+A16*A24*A42*A55*A61-A14*A26*A42*A55* &
      A61-A16*A22*A44*A55*A61+A12*A26*A44*A55*A61+A14*A22*A46*A55*A61-A12*A24* &
      A46*A55*A61-A15*A24*A42*A56*A61+A14*A25*A42*A56*A61+A15*A22*A44*A56* &
      A61-A12*A25*A44*A56*A61-A14*A22*A45*A56*A61+A12*A24*A45*A56*A61-A16*A25* &
      A44*A51*A62+A15*A26*A44*A51*A62+A16*A24*A45*A51*A62-A14*A26*A45*A51* &
      A62-A15*A24*A46*A51*A62+A14*A25*A46*A51*A62+A16*A25*A41*A54*A62-A15*A26* &
      A41*A54*A62-A16*A21*A45*A54*A62+A11*A26*A45*A54*A62+A15*A21*A46*A54* &
      A62-A11*A25*A46*A54*A62-A16*A24*A41*A55*A62+A14*A26*A41*A55*A62+A16*A21* &
      A44*A55*A62-A11*A26*A44*A55*A62-A14*A21*A46*A55*A62+A11*A24*A46*A55* &
      A62+A15*A24*A41*A56*A62-A14*A25*A41*A56*A62-A15*A21*A44*A56*A62+A11*A25* &
      A44*A56*A62+A14*A21*A45*A56*A62-A11*A24*A45*A56*A62+A16*A25*A42*A51* &
      A64-A15*A26*A42*A51*A64-A16*A22*A45*A51*A64+A12*A26*A45*A51*A64+A15*A22* &
      A46*A51*A64-A12*A25*A46*A51*A64-A16*A25*A41*A52*A64+A15*A26*A41*A52* &
      A64+A16*A21*A45*A52*A64-A11*A26*A45*A52*A64-A15*A21*A46*A52*A64+A11*A25* &
      A46*A52*A64+A16*A22*A41*A55*A64-A12*A26*A41*A55*A64-A16*A21*A42*A55* &
      A64+A11*A26*A42*A55*A64+A12*A21*A46*A55*A64-A11*A22*A46*A55*A64-A15*A22* &
      A41*A56*A64+A12*A25*A41*A56*A64+A15*A21*A42*A56*A64-A11*A25*A42*A56* &
      A64-A12*A21*A45*A56*A64+A11*A22*A45*A56*A64-A16*A24*A42*A51*A65+A14*A26* &
      A42*A51*A65+A16*A22*A44*A51*A65-A12*A26*A44*A51*A65-A14*A22*A46*A51* &
      A65+A12*A24*A46*A51*A65+A16*A24*A41*A52*A65-A14*A26*A41*A52*A65-A16*A21* &
      A44*A52*A65+A11*A26*A44*A52*A65+A14*A21*A46*A52*A65-A11*A24*A46*A52* &
      A65-A16*A22*A41*A54*A65+A12*A26*A41*A54*A65+A16*A21*A42*A54*A65-A11*A26* &
      A42*A54*A65-A12*A21*A46*A54*A65+A11*A22*A46*A54*A65+A14*A22*A41*A56* &
      A65-A12*A24*A41*A56*A65-A14*A21*A42*A56*A65+A11*A24*A42*A56*A65+A12*A21* &
      A44*A56*A65-A11*A22*A44*A56*A65+A15*A24*A42*A51*A66-A14*A25*A42*A51* &
      A66-A15*A22*A44*A51*A66+A12*A25*A44*A51*A66+A14*A22*A45*A51*A66-A12*A24* &
      A45*A51*A66-A15*A24*A41*A52*A66+A14*A25*A41*A52*A66+A15*A21*A44*A52* &
      A66-A11*A25*A44*A52*A66-A14*A21*A45*A52*A66+A11*A24*A45*A52*A66+A15*A22* &
      A41*A54*A66-A12*A25*A41*A54*A66-A15*A21*A42*A54*A66+A11*A25*A42*A54* &
      A66+A12*A21*A45*A54*A66-A11*A22*A45*A54*A66-A14*A22*A41*A55*A66+A12*A24* &
      A41*A55*A66+A14*A21*A42*A55*A66-A11*A24*A42*A55*A66-A12*A21*A44*A55* &
      A66+A11*A22*A44*A55*A66

      COFACTOR(4,3) = -A16*A25*A34*A52*A61+A15*A26*A34*A52*A61+A16*A24* &
      A35*A52*A61-A14*A26*A35*A52*A61-A15*A24*A36*A52*A61+A14*A25*A36*A52* &
      A61+A16*A25*A32*A54*A61-A15*A26*A32*A54*A61-A16*A22*A35*A54*A61+A12*A26* &
      A35*A54*A61+A15*A22*A36*A54*A61-A12*A25*A36*A54*A61-A16*A24*A32*A55* &
      A61+A14*A26*A32*A55*A61+A16*A22*A34*A55*A61-A12*A26*A34*A55*A61-A14*A22* &
      A36*A55*A61+A12*A24*A36*A55*A61+A15*A24*A32*A56*A61-A14*A25*A32*A56* &
      A61-A15*A22*A34*A56*A61+A12*A25*A34*A56*A61+A14*A22*A35*A56*A61-A12*A24* &
      A35*A56*A61+A16*A25*A34*A51*A62-A15*A26*A34*A51*A62-A16*A24*A35*A51* &
      A62+A14*A26*A35*A51*A62+A15*A24*A36*A51*A62-A14*A25*A36*A51*A62-A16*A25* &
      A31*A54*A62+A15*A26*A31*A54*A62+A16*A21*A35*A54*A62-A11*A26*A35*A54* &
      A62-A15*A21*A36*A54*A62+A11*A25*A36*A54*A62+A16*A24*A31*A55*A62-A14*A26* &
      A31*A55*A62-A16*A21*A34*A55*A62+A11*A26*A34*A55*A62+A14*A21*A36*A55* &
      A62-A11*A24*A36*A55*A62-A15*A24*A31*A56*A62+A14*A25*A31*A56*A62+A15*A21* &
      A34*A56*A62-A11*A25*A34*A56*A62-A14*A21*A35*A56*A62+A11*A24*A35*A56* &
      A62-A16*A25*A32*A51*A64+A15*A26*A32*A51*A64+A16*A22*A35*A51*A64-A12*A26* &
      A35*A51*A64-A15*A22*A36*A51*A64+A12*A25*A36*A51*A64+A16*A25*A31*A52* &
      A64-A15*A26*A31*A52*A64-A16*A21*A35*A52*A64+A11*A26*A35*A52*A64+A15*A21* &
      A36*A52*A64-A11*A25*A36*A52*A64-A16*A22*A31*A55*A64+A12*A26*A31*A55* &
      A64+A16*A21*A32*A55*A64-A11*A26*A32*A55*A64-A12*A21*A36*A55*A64+A11*A22* &
      A36*A55*A64+A15*A22*A31*A56*A64-A12*A25*A31*A56*A64-A15*A21*A32*A56* &
      A64+A11*A25*A32*A56*A64+A12*A21*A35*A56*A64-A11*A22*A35*A56*A64+A16*A24* &
      A32*A51*A65-A14*A26*A32*A51*A65-A16*A22*A34*A51*A65+A12*A26*A34*A51* &
      A65+A14*A22*A36*A51*A65-A12*A24*A36*A51*A65-A16*A24*A31*A52*A65+A14*A26* &
      A31*A52*A65+A16*A21*A34*A52*A65-A11*A26*A34*A52*A65-A14*A21*A36*A52* &
      A65+A11*A24*A36*A52*A65+A16*A22*A31*A54*A65-A12*A26*A31*A54*A65-A16*A21* &
      A32*A54*A65+A11*A26*A32*A54*A65+A12*A21*A36*A54*A65-A11*A22*A36*A54* &
      A65-A14*A22*A31*A56*A65+A12*A24*A31*A56*A65+A14*A21*A32*A56*A65-A11*A24* &
      A32*A56*A65-A12*A21*A34*A56*A65+A11*A22*A34*A56*A65-A15*A24*A32*A51* &
      A66+A14*A25*A32*A51*A66+A15*A22*A34*A51*A66-A12*A25*A34*A51*A66-A14*A22* &
      A35*A51*A66+A12*A24*A35*A51*A66+A15*A24*A31*A52*A66-A14*A25*A31*A52* &
      A66-A15*A21*A34*A52*A66+A11*A25*A34*A52*A66+A14*A21*A35*A52*A66-A11*A24* &
      A35*A52*A66-A15*A22*A31*A54*A66+A12*A25*A31*A54*A66+A15*A21*A32*A54* &
      A66-A11*A25*A32*A54*A66-A12*A21*A35*A54*A66+A11*A22*A35*A54*A66+A14*A22* &
      A31*A55*A66-A12*A24*A31*A55*A66-A14*A21*A32*A55*A66+A11*A24*A32*A55* &
      A66+A12*A21*A34*A55*A66-A11*A22*A34*A55*A66

      COFACTOR(5,3) = A16*A25*A34*A42*A61-A15*A26* &
      A34*A42*A61-A16*A24*A35*A42*A61+A14*A26*A35*A42*A61+A15*A24*A36*A42* &
      A61-A14*A25*A36*A42*A61-A16*A25*A32*A44*A61+A15*A26*A32*A44*A61+A16*A22* &
      A35*A44*A61-A12*A26*A35*A44*A61-A15*A22*A36*A44*A61+A12*A25*A36*A44* &
      A61+A16*A24*A32*A45*A61-A14*A26*A32*A45*A61-A16*A22*A34*A45*A61+A12*A26* &
      A34*A45*A61+A14*A22*A36*A45*A61-A12*A24*A36*A45*A61-A15*A24*A32*A46* &
      A61+A14*A25*A32*A46*A61+A15*A22*A34*A46*A61-A12*A25*A34*A46*A61-A14*A22* &
      A35*A46*A61+A12*A24*A35*A46*A61-A16*A25*A34*A41*A62+A15*A26*A34*A41* &
      A62+A16*A24*A35*A41*A62-A14*A26*A35*A41*A62-A15*A24*A36*A41*A62+A14*A25* &
      A36*A41*A62+A16*A25*A31*A44*A62-A15*A26*A31*A44*A62-A16*A21*A35*A44* &
      A62+A11*A26*A35*A44*A62+A15*A21*A36*A44*A62-A11*A25*A36*A44*A62-A16*A24* &
      A31*A45*A62+A14*A26*A31*A45*A62+A16*A21*A34*A45*A62-A11*A26*A34*A45* &
      A62-A14*A21*A36*A45*A62+A11*A24*A36*A45*A62+A15*A24*A31*A46*A62-A14*A25* &
      A31*A46*A62-A15*A21*A34*A46*A62+A11*A25*A34*A46*A62+A14*A21*A35*A46* &
      A62-A11*A24*A35*A46*A62+A16*A25*A32*A41*A64-A15*A26*A32*A41*A64-A16*A22* &
      A35*A41*A64+A12*A26*A35*A41*A64+A15*A22*A36*A41*A64-A12*A25*A36*A41* &
      A64-A16*A25*A31*A42*A64+A15*A26*A31*A42*A64+A16*A21*A35*A42*A64-A11*A26* &
      A35*A42*A64-A15*A21*A36*A42*A64+A11*A25*A36*A42*A64+A16*A22*A31*A45* &
      A64-A12*A26*A31*A45*A64-A16*A21*A32*A45*A64+A11*A26*A32*A45*A64+A12*A21* &
      A36*A45*A64-A11*A22*A36*A45*A64-A15*A22*A31*A46*A64+A12*A25*A31*A46* &
      A64+A15*A21*A32*A46*A64-A11*A25*A32*A46*A64-A12*A21*A35*A46*A64+A11*A22* &
      A35*A46*A64-A16*A24*A32*A41*A65+A14*A26*A32*A41*A65+A16*A22*A34*A41* &
      A65-A12*A26*A34*A41*A65-A14*A22*A36*A41*A65+A12*A24*A36*A41*A65+A16*A24* &
      A31*A42*A65-A14*A26*A31*A42*A65-A16*A21*A34*A42*A65+A11*A26*A34*A42* &
      A65+A14*A21*A36*A42*A65-A11*A24*A36*A42*A65-A16*A22*A31*A44*A65+A12*A26* &
      A31*A44*A65+A16*A21*A32*A44*A65-A11*A26*A32*A44*A65-A12*A21*A36*A44* &
      A65+A11*A22*A36*A44*A65+A14*A22*A31*A46*A65-A12*A24*A31*A46*A65-A14*A21* &
      A32*A46*A65+A11*A24*A32*A46*A65+A12*A21*A34*A46*A65-A11*A22*A34*A46* &
      A65+A15*A24*A32*A41*A66-A14*A25*A32*A41*A66-A15*A22*A34*A41*A66+A12*A25* &
      A34*A41*A66+A14*A22*A35*A41*A66-A12*A24*A35*A41*A66-A15*A24*A31*A42* &
      A66+A14*A25*A31*A42*A66+A15*A21*A34*A42*A66-A11*A25*A34*A42*A66-A14*A21* &
      A35*A42*A66+A11*A24*A35*A42*A66+A15*A22*A31*A44*A66-A12*A25*A31*A44* &
      A66-A15*A21*A32*A44*A66+A11*A25*A32*A44*A66+A12*A21*A35*A44*A66-A11*A22* &
      A35*A44*A66-A14*A22*A31*A45*A66+A12*A24*A31*A45*A66+A14*A21*A32*A45* &
      A66-A11*A24*A32*A45*A66-A12*A21*A34*A45*A66+A11*A22*A34*A45*A66

      COFACTOR(6,3) = -A16*A25* &
      A34*A42*A51+A15*A26*A34*A42*A51+A16*A24*A35*A42*A51-A14*A26*A35*A42* &
      A51-A15*A24*A36*A42*A51+A14*A25*A36*A42*A51+A16*A25*A32*A44*A51-A15*A26* &
      A32*A44*A51-A16*A22*A35*A44*A51+A12*A26*A35*A44*A51+A15*A22*A36*A44* &
      A51-A12*A25*A36*A44*A51-A16*A24*A32*A45*A51+A14*A26*A32*A45*A51+A16*A22* &
      A34*A45*A51-A12*A26*A34*A45*A51-A14*A22*A36*A45*A51+A12*A24*A36*A45* &
      A51+A15*A24*A32*A46*A51-A14*A25*A32*A46*A51-A15*A22*A34*A46*A51+A12*A25* &
      A34*A46*A51+A14*A22*A35*A46*A51-A12*A24*A35*A46*A51+A16*A25*A34*A41* &
      A52-A15*A26*A34*A41*A52-A16*A24*A35*A41*A52+A14*A26*A35*A41*A52+A15*A24* &
      A36*A41*A52-A14*A25*A36*A41*A52-A16*A25*A31*A44*A52+A15*A26*A31*A44* &
      A52+A16*A21*A35*A44*A52-A11*A26*A35*A44*A52-A15*A21*A36*A44*A52+A11*A25* &
      A36*A44*A52+A16*A24*A31*A45*A52-A14*A26*A31*A45*A52-A16*A21*A34*A45* &
      A52+A11*A26*A34*A45*A52+A14*A21*A36*A45*A52-A11*A24*A36*A45*A52-A15*A24* &
      A31*A46*A52+A14*A25*A31*A46*A52+A15*A21*A34*A46*A52-A11*A25*A34*A46* &
      A52-A14*A21*A35*A46*A52+A11*A24*A35*A46*A52-A16*A25*A32*A41*A54+A15*A26* &
      A32*A41*A54+A16*A22*A35*A41*A54-A12*A26*A35*A41*A54-A15*A22*A36*A41* &
      A54+A12*A25*A36*A41*A54+A16*A25*A31*A42*A54-A15*A26*A31*A42*A54-A16*A21* &
      A35*A42*A54+A11*A26*A35*A42*A54+A15*A21*A36*A42*A54-A11*A25*A36*A42* &
      A54-A16*A22*A31*A45*A54+A12*A26*A31*A45*A54+A16*A21*A32*A45*A54-A11*A26* &
      A32*A45*A54-A12*A21*A36*A45*A54+A11*A22*A36*A45*A54+A15*A22*A31*A46* &
      A54-A12*A25*A31*A46*A54-A15*A21*A32*A46*A54+A11*A25*A32*A46*A54+A12*A21* &
      A35*A46*A54-A11*A22*A35*A46*A54+A16*A24*A32*A41*A55-A14*A26*A32*A41* &
      A55-A16*A22*A34*A41*A55+A12*A26*A34*A41*A55+A14*A22*A36*A41*A55-A12*A24* &
      A36*A41*A55-A16*A24*A31*A42*A55+A14*A26*A31*A42*A55+A16*A21*A34*A42* &
      A55-A11*A26*A34*A42*A55-A14*A21*A36*A42*A55+A11*A24*A36*A42*A55+A16*A22* &
      A31*A44*A55-A12*A26*A31*A44*A55-A16*A21*A32*A44*A55+A11*A26*A32*A44* &
      A55+A12*A21*A36*A44*A55-A11*A22*A36*A44*A55-A14*A22*A31*A46*A55+A12*A24* &
      A31*A46*A55+A14*A21*A32*A46*A55-A11*A24*A32*A46*A55-A12*A21*A34*A46* &
      A55+A11*A22*A34*A46*A55-A15*A24*A32*A41*A56+A14*A25*A32*A41*A56+A15*A22* &
      A34*A41*A56-A12*A25*A34*A41*A56-A14*A22*A35*A41*A56+A12*A24*A35*A41* &
      A56+A15*A24*A31*A42*A56-A14*A25*A31*A42*A56-A15*A21*A34*A42*A56+A11*A25* &
      A34*A42*A56+A14*A21*A35*A42*A56-A11*A24*A35*A42*A56-A15*A22*A31*A44* &
      A56+A12*A25*A31*A44*A56+A15*A21*A32*A44*A56-A11*A25*A32*A44*A56-A12*A21* &
      A35*A44*A56+A11*A22*A35*A44*A56+A14*A22*A31*A45*A56-A12*A24*A31*A45* &
      A56-A14*A21*A32*A45*A56+A11*A24*A32*A45*A56+A12*A21*A34*A45*A56-A11*A22* &
      A34*A45*A56

      COFACTOR(1,4) = -A26*A35*A43*A52*A61+A25*A36*A43*A52*A61+A26*A33*A45*A52* &
      A61-A23*A36*A45*A52*A61-A25*A33*A46*A52*A61+A23*A35*A46*A52*A61+A26*A35* &
      A42*A53*A61-A25*A36*A42*A53*A61-A26*A32*A45*A53*A61+A22*A36*A45*A53* &
      A61+A25*A32*A46*A53*A61-A22*A35*A46*A53*A61-A26*A33*A42*A55*A61+A23*A36* &
      A42*A55*A61+A26*A32*A43*A55*A61-A22*A36*A43*A55*A61-A23*A32*A46*A55* &
      A61+A22*A33*A46*A55*A61+A25*A33*A42*A56*A61-A23*A35*A42*A56*A61-A25*A32* &
      A43*A56*A61+A22*A35*A43*A56*A61+A23*A32*A45*A56*A61-A22*A33*A45*A56* &
      A61+A26*A35*A43*A51*A62-A25*A36*A43*A51*A62-A26*A33*A45*A51*A62+A23*A36* &
      A45*A51*A62+A25*A33*A46*A51*A62-A23*A35*A46*A51*A62-A26*A35*A41*A53* &
      A62+A25*A36*A41*A53*A62+A26*A31*A45*A53*A62-A21*A36*A45*A53*A62-A25*A31* &
      A46*A53*A62+A21*A35*A46*A53*A62+A26*A33*A41*A55*A62-A23*A36*A41*A55* &
      A62-A26*A31*A43*A55*A62+A21*A36*A43*A55*A62+A23*A31*A46*A55*A62-A21*A33* &
      A46*A55*A62-A25*A33*A41*A56*A62+A23*A35*A41*A56*A62+A25*A31*A43*A56* &
      A62-A21*A35*A43*A56*A62-A23*A31*A45*A56*A62+A21*A33*A45*A56*A62-A26*A35* &
      A42*A51*A63+A25*A36*A42*A51*A63+A26*A32*A45*A51*A63-A22*A36*A45*A51* &
      A63-A25*A32*A46*A51*A63+A22*A35*A46*A51*A63+A26*A35*A41*A52*A63-A25*A36* &
      A41*A52*A63-A26*A31*A45*A52*A63+A21*A36*A45*A52*A63+A25*A31*A46*A52* &
      A63-A21*A35*A46*A52*A63-A26*A32*A41*A55*A63+A22*A36*A41*A55*A63+A26*A31* &
      A42*A55*A63-A21*A36*A42*A55*A63-A22*A31*A46*A55*A63+A21*A32*A46*A55* &
      A63+A25*A32*A41*A56*A63-A22*A35*A41*A56*A63-A25*A31*A42*A56*A63+A21*A35* &
      A42*A56*A63+A22*A31*A45*A56*A63-A21*A32*A45*A56*A63+A26*A33*A42*A51* &
      A65-A23*A36*A42*A51*A65-A26*A32*A43*A51*A65+A22*A36*A43*A51*A65+A23*A32* &
      A46*A51*A65-A22*A33*A46*A51*A65-A26*A33*A41*A52*A65+A23*A36*A41*A52* &
      A65+A26*A31*A43*A52*A65-A21*A36*A43*A52*A65-A23*A31*A46*A52*A65+A21*A33* &
      A46*A52*A65+A26*A32*A41*A53*A65-A22*A36*A41*A53*A65-A26*A31*A42*A53* &
      A65+A21*A36*A42*A53*A65+A22*A31*A46*A53*A65-A21*A32*A46*A53*A65-A23*A32* &
      A41*A56*A65+A22*A33*A41*A56*A65+A23*A31*A42*A56*A65-A21*A33*A42*A56* &
      A65-A22*A31*A43*A56*A65+A21*A32*A43*A56*A65-A25*A33*A42*A51*A66+A23*A35* &
      A42*A51*A66+A25*A32*A43*A51*A66-A22*A35*A43*A51*A66-A23*A32*A45*A51* &
      A66+A22*A33*A45*A51*A66+A25*A33*A41*A52*A66-A23*A35*A41*A52*A66-A25*A31* &
      A43*A52*A66+A21*A35*A43*A52*A66+A23*A31*A45*A52*A66-A21*A33*A45*A52* &
      A66-A25*A32*A41*A53*A66+A22*A35*A41*A53*A66+A25*A31*A42*A53*A66-A21*A35* &
      A42*A53*A66-A22*A31*A45*A53*A66+A21*A32*A45*A53*A66+A23*A32*A41*A55* &
      A66-A22*A33*A41*A55*A66-A23*A31*A42*A55*A66+A21*A33*A42*A55*A66+A22*A31* &
      A43*A55*A66-A21*A32*A43*A55*A66

      COFACTOR(2,4) = A16*A35*A43*A52*A61-A15*A36*A43*A52* &
      A61-A16*A33*A45*A52*A61+A13*A36*A45*A52*A61+A15*A33*A46*A52*A61-A13*A35* &
      A46*A52*A61-A16*A35*A42*A53*A61+A15*A36*A42*A53*A61+A16*A32*A45*A53* &
      A61-A12*A36*A45*A53*A61-A15*A32*A46*A53*A61+A12*A35*A46*A53*A61+A16*A33* &
      A42*A55*A61-A13*A36*A42*A55*A61-A16*A32*A43*A55*A61+A12*A36*A43*A55* &
      A61+A13*A32*A46*A55*A61-A12*A33*A46*A55*A61-A15*A33*A42*A56*A61+A13*A35* &
      A42*A56*A61+A15*A32*A43*A56*A61-A12*A35*A43*A56*A61-A13*A32*A45*A56* &
      A61+A12*A33*A45*A56*A61-A16*A35*A43*A51*A62+A15*A36*A43*A51*A62+A16*A33* &
      A45*A51*A62-A13*A36*A45*A51*A62-A15*A33*A46*A51*A62+A13*A35*A46*A51* &
      A62+A16*A35*A41*A53*A62-A15*A36*A41*A53*A62-A16*A31*A45*A53*A62+A11*A36* &
      A45*A53*A62+A15*A31*A46*A53*A62-A11*A35*A46*A53*A62-A16*A33*A41*A55* &
      A62+A13*A36*A41*A55*A62+A16*A31*A43*A55*A62-A11*A36*A43*A55*A62-A13*A31* &
      A46*A55*A62+A11*A33*A46*A55*A62+A15*A33*A41*A56*A62-A13*A35*A41*A56* &
      A62-A15*A31*A43*A56*A62+A11*A35*A43*A56*A62+A13*A31*A45*A56*A62-A11*A33* &
      A45*A56*A62+A16*A35*A42*A51*A63-A15*A36*A42*A51*A63-A16*A32*A45*A51* &
      A63+A12*A36*A45*A51*A63+A15*A32*A46*A51*A63-A12*A35*A46*A51*A63-A16*A35* &
      A41*A52*A63+A15*A36*A41*A52*A63+A16*A31*A45*A52*A63-A11*A36*A45*A52* &
      A63-A15*A31*A46*A52*A63+A11*A35*A46*A52*A63+A16*A32*A41*A55*A63-A12*A36* &
      A41*A55*A63-A16*A31*A42*A55*A63+A11*A36*A42*A55*A63+A12*A31*A46*A55* &
      A63-A11*A32*A46*A55*A63-A15*A32*A41*A56*A63+A12*A35*A41*A56*A63+A15*A31* &
      A42*A56*A63-A11*A35*A42*A56*A63-A12*A31*A45*A56*A63+A11*A32*A45*A56* &
      A63-A16*A33*A42*A51*A65+A13*A36*A42*A51*A65+A16*A32*A43*A51*A65-A12*A36* &
      A43*A51*A65-A13*A32*A46*A51*A65+A12*A33*A46*A51*A65+A16*A33*A41*A52* &
      A65-A13*A36*A41*A52*A65-A16*A31*A43*A52*A65+A11*A36*A43*A52*A65+A13*A31* &
      A46*A52*A65-A11*A33*A46*A52*A65-A16*A32*A41*A53*A65+A12*A36*A41*A53* &
      A65+A16*A31*A42*A53*A65-A11*A36*A42*A53*A65-A12*A31*A46*A53*A65+A11*A32* &
      A46*A53*A65+A13*A32*A41*A56*A65-A12*A33*A41*A56*A65-A13*A31*A42*A56* &
      A65+A11*A33*A42*A56*A65+A12*A31*A43*A56*A65-A11*A32*A43*A56*A65+A15*A33* &
      A42*A51*A66-A13*A35*A42*A51*A66-A15*A32*A43*A51*A66+A12*A35*A43*A51* &
      A66+A13*A32*A45*A51*A66-A12*A33*A45*A51*A66-A15*A33*A41*A52*A66+A13*A35* &
      A41*A52*A66+A15*A31*A43*A52*A66-A11*A35*A43*A52*A66-A13*A31*A45*A52* &
      A66+A11*A33*A45*A52*A66+A15*A32*A41*A53*A66-A12*A35*A41*A53*A66-A15*A31* &
      A42*A53*A66+A11*A35*A42*A53*A66+A12*A31*A45*A53*A66-A11*A32*A45*A53* &
      A66-A13*A32*A41*A55*A66+A12*A33*A41*A55*A66+A13*A31*A42*A55*A66-A11*A33* &
      A42*A55*A66-A12*A31*A43*A55*A66+A11*A32*A43*A55*A66

      COFACTOR(3,4) = -A16*A25*A43*A52* &
      A61+A15*A26*A43*A52*A61+A16*A23*A45*A52*A61-A13*A26*A45*A52*A61-A15*A23* &
      A46*A52*A61+A13*A25*A46*A52*A61+A16*A25*A42*A53*A61-A15*A26*A42*A53* &
      A61-A16*A22*A45*A53*A61+A12*A26*A45*A53*A61+A15*A22*A46*A53*A61-A12*A25* &
      A46*A53*A61-A16*A23*A42*A55*A61+A13*A26*A42*A55*A61+A16*A22*A43*A55* &
      A61-A12*A26*A43*A55*A61-A13*A22*A46*A55*A61+A12*A23*A46*A55*A61+A15*A23* &
      A42*A56*A61-A13*A25*A42*A56*A61-A15*A22*A43*A56*A61+A12*A25*A43*A56* &
      A61+A13*A22*A45*A56*A61-A12*A23*A45*A56*A61+A16*A25*A43*A51*A62-A15*A26* &
      A43*A51*A62-A16*A23*A45*A51*A62+A13*A26*A45*A51*A62+A15*A23*A46*A51* &
      A62-A13*A25*A46*A51*A62-A16*A25*A41*A53*A62+A15*A26*A41*A53*A62+A16*A21* &
      A45*A53*A62-A11*A26*A45*A53*A62-A15*A21*A46*A53*A62+A11*A25*A46*A53* &
      A62+A16*A23*A41*A55*A62-A13*A26*A41*A55*A62-A16*A21*A43*A55*A62+A11*A26* &
      A43*A55*A62+A13*A21*A46*A55*A62-A11*A23*A46*A55*A62-A15*A23*A41*A56* &
      A62+A13*A25*A41*A56*A62+A15*A21*A43*A56*A62-A11*A25*A43*A56*A62-A13*A21* &
      A45*A56*A62+A11*A23*A45*A56*A62-A16*A25*A42*A51*A63+A15*A26*A42*A51* &
      A63+A16*A22*A45*A51*A63-A12*A26*A45*A51*A63-A15*A22*A46*A51*A63+A12*A25* &
      A46*A51*A63+A16*A25*A41*A52*A63-A15*A26*A41*A52*A63-A16*A21*A45*A52* &
      A63+A11*A26*A45*A52*A63+A15*A21*A46*A52*A63-A11*A25*A46*A52*A63-A16*A22* &
      A41*A55*A63+A12*A26*A41*A55*A63+A16*A21*A42*A55*A63-A11*A26*A42*A55* &
      A63-A12*A21*A46*A55*A63+A11*A22*A46*A55*A63+A15*A22*A41*A56*A63-A12*A25* &
      A41*A56*A63-A15*A21*A42*A56*A63+A11*A25*A42*A56*A63+A12*A21*A45*A56* &
      A63-A11*A22*A45*A56*A63+A16*A23*A42*A51*A65-A13*A26*A42*A51*A65-A16*A22* &
      A43*A51*A65+A12*A26*A43*A51*A65+A13*A22*A46*A51*A65-A12*A23*A46*A51* &
      A65-A16*A23*A41*A52*A65+A13*A26*A41*A52*A65+A16*A21*A43*A52*A65-A11*A26* &
      A43*A52*A65-A13*A21*A46*A52*A65+A11*A23*A46*A52*A65+A16*A22*A41*A53* &
      A65-A12*A26*A41*A53*A65-A16*A21*A42*A53*A65+A11*A26*A42*A53*A65+A12*A21* &
      A46*A53*A65-A11*A22*A46*A53*A65-A13*A22*A41*A56*A65+A12*A23*A41*A56* &
      A65+A13*A21*A42*A56*A65-A11*A23*A42*A56*A65-A12*A21*A43*A56*A65+A11*A22* &
      A43*A56*A65-A15*A23*A42*A51*A66+A13*A25*A42*A51*A66+A15*A22*A43*A51* &
      A66-A12*A25*A43*A51*A66-A13*A22*A45*A51*A66+A12*A23*A45*A51*A66+A15*A23* &
      A41*A52*A66-A13*A25*A41*A52*A66-A15*A21*A43*A52*A66+A11*A25*A43*A52* &
      A66+A13*A21*A45*A52*A66-A11*A23*A45*A52*A66-A15*A22*A41*A53*A66+A12*A25* &
      A41*A53*A66+A15*A21*A42*A53*A66-A11*A25*A42*A53*A66-A12*A21*A45*A53* &
      A66+A11*A22*A45*A53*A66+A13*A22*A41*A55*A66-A12*A23*A41*A55*A66-A13*A21* &
      A42*A55*A66+A11*A23*A42*A55*A66+A12*A21*A43*A55*A66-A11*A22*A43*A55* &
      A66 

      COFACTOR(4,4) = A16*A25*A33*A52*A61-A15*A26*A33*A52*A61-A16*A23*A35*A52*A61+A13*A26*     &
      A35*A52*A61+A15*A23*A36*A52*A61-A13*A25*A36*A52*A61-A16*A25*A32*A53* &
      A61+A15*A26*A32*A53*A61+A16*A22*A35*A53*A61-A12*A26*A35*A53*A61-A15*A22* &
      A36*A53*A61+A12*A25*A36*A53*A61+A16*A23*A32*A55*A61-A13*A26*A32*A55* &
      A61-A16*A22*A33*A55*A61+A12*A26*A33*A55*A61+A13*A22*A36*A55*A61-A12*A23* &
      A36*A55*A61-A15*A23*A32*A56*A61+A13*A25*A32*A56*A61+A15*A22*A33*A56* &
      A61-A12*A25*A33*A56*A61-A13*A22*A35*A56*A61+A12*A23*A35*A56*A61-A16*A25* &
      A33*A51*A62+A15*A26*A33*A51*A62+A16*A23*A35*A51*A62-A13*A26*A35*A51* &
      A62-A15*A23*A36*A51*A62+A13*A25*A36*A51*A62+A16*A25*A31*A53*A62-A15*A26* &
      A31*A53*A62-A16*A21*A35*A53*A62+A11*A26*A35*A53*A62+A15*A21*A36*A53* &
      A62-A11*A25*A36*A53*A62-A16*A23*A31*A55*A62+A13*A26*A31*A55*A62+A16*A21* &
      A33*A55*A62-A11*A26*A33*A55*A62-A13*A21*A36*A55*A62+A11*A23*A36*A55* &
      A62+A15*A23*A31*A56*A62-A13*A25*A31*A56*A62-A15*A21*A33*A56*A62+A11*A25* &
      A33*A56*A62+A13*A21*A35*A56*A62-A11*A23*A35*A56*A62+A16*A25*A32*A51* &
      A63-A15*A26*A32*A51*A63-A16*A22*A35*A51*A63+A12*A26*A35*A51*A63+A15*A22* &
      A36*A51*A63-A12*A25*A36*A51*A63-A16*A25*A31*A52*A63+A15*A26*A31*A52* &
      A63+A16*A21*A35*A52*A63-A11*A26*A35*A52*A63-A15*A21*A36*A52*A63+A11*A25* &
      A36*A52*A63+A16*A22*A31*A55*A63-A12*A26*A31*A55*A63-A16*A21*A32*A55* &
      A63+A11*A26*A32*A55*A63+A12*A21*A36*A55*A63-A11*A22*A36*A55*A63-A15*A22* &
      A31*A56*A63+A12*A25*A31*A56*A63+A15*A21*A32*A56*A63-A11*A25*A32*A56* &
      A63-A12*A21*A35*A56*A63+A11*A22*A35*A56*A63-A16*A23*A32*A51*A65+A13*A26* &
      A32*A51*A65+A16*A22*A33*A51*A65-A12*A26*A33*A51*A65-A13*A22*A36*A51* &
      A65+A12*A23*A36*A51*A65+A16*A23*A31*A52*A65-A13*A26*A31*A52*A65-A16*A21* &
      A33*A52*A65+A11*A26*A33*A52*A65+A13*A21*A36*A52*A65-A11*A23*A36*A52* &
      A65-A16*A22*A31*A53*A65+A12*A26*A31*A53*A65+A16*A21*A32*A53*A65-A11*A26* &
      A32*A53*A65-A12*A21*A36*A53*A65+A11*A22*A36*A53*A65+A13*A22*A31*A56* &
      A65-A12*A23*A31*A56*A65-A13*A21*A32*A56*A65+A11*A23*A32*A56*A65+A12*A21* &
      A33*A56*A65-A11*A22*A33*A56*A65+A15*A23*A32*A51*A66-A13*A25*A32*A51* &
      A66-A15*A22*A33*A51*A66+A12*A25*A33*A51*A66+A13*A22*A35*A51*A66-A12*A23* &
      A35*A51*A66-A15*A23*A31*A52*A66+A13*A25*A31*A52*A66+A15*A21*A33*A52* &
      A66-A11*A25*A33*A52*A66-A13*A21*A35*A52*A66+A11*A23*A35*A52*A66+A15*A22* &
      A31*A53*A66-A12*A25*A31*A53*A66-A15*A21*A32*A53*A66+A11*A25*A32*A53* &
      A66+A12*A21*A35*A53*A66-A11*A22*A35*A53*A66-A13*A22*A31*A55*A66+A12*A23* &
      A31*A55*A66+A13*A21*A32*A55*A66-A11*A23*A32*A55*A66-A12*A21*A33*A55* &
      A66+A11*A22*A33*A55*A66

      COFACTOR(5,4) = -A16*A25*A33*A42*A61+A15*A26*A33*A42*A61+A16*A23* &
      A35*A42*A61-A13*A26*A35*A42*A61-A15*A23*A36*A42*A61+A13*A25*A36*A42* &
      A61+A16*A25*A32*A43*A61-A15*A26*A32*A43*A61-A16*A22*A35*A43*A61+A12*A26* &
      A35*A43*A61+A15*A22*A36*A43*A61-A12*A25*A36*A43*A61-A16*A23*A32*A45* &
      A61+A13*A26*A32*A45*A61+A16*A22*A33*A45*A61-A12*A26*A33*A45*A61-A13*A22* &
      A36*A45*A61+A12*A23*A36*A45*A61+A15*A23*A32*A46*A61-A13*A25*A32*A46* &
      A61-A15*A22*A33*A46*A61+A12*A25*A33*A46*A61+A13*A22*A35*A46*A61-A12*A23* &
      A35*A46*A61+A16*A25*A33*A41*A62-A15*A26*A33*A41*A62-A16*A23*A35*A41* &
      A62+A13*A26*A35*A41*A62+A15*A23*A36*A41*A62-A13*A25*A36*A41*A62-A16*A25* &
      A31*A43*A62+A15*A26*A31*A43*A62+A16*A21*A35*A43*A62-A11*A26*A35*A43* &
      A62-A15*A21*A36*A43*A62+A11*A25*A36*A43*A62+A16*A23*A31*A45*A62-A13*A26* &
      A31*A45*A62-A16*A21*A33*A45*A62+A11*A26*A33*A45*A62+A13*A21*A36*A45* &
      A62-A11*A23*A36*A45*A62-A15*A23*A31*A46*A62+A13*A25*A31*A46*A62+A15*A21* &
      A33*A46*A62-A11*A25*A33*A46*A62-A13*A21*A35*A46*A62+A11*A23*A35*A46* &
      A62-A16*A25*A32*A41*A63+A15*A26*A32*A41*A63+A16*A22*A35*A41*A63-A12*A26* &
      A35*A41*A63-A15*A22*A36*A41*A63+A12*A25*A36*A41*A63+A16*A25*A31*A42* &
      A63-A15*A26*A31*A42*A63-A16*A21*A35*A42*A63+A11*A26*A35*A42*A63+A15*A21* &
      A36*A42*A63-A11*A25*A36*A42*A63-A16*A22*A31*A45*A63+A12*A26*A31*A45* &
      A63+A16*A21*A32*A45*A63-A11*A26*A32*A45*A63-A12*A21*A36*A45*A63+A11*A22* &
      A36*A45*A63+A15*A22*A31*A46*A63-A12*A25*A31*A46*A63-A15*A21*A32*A46* &
      A63+A11*A25*A32*A46*A63+A12*A21*A35*A46*A63-A11*A22*A35*A46*A63+A16*A23* &
      A32*A41*A65-A13*A26*A32*A41*A65-A16*A22*A33*A41*A65+A12*A26*A33*A41* &
      A65+A13*A22*A36*A41*A65-A12*A23*A36*A41*A65-A16*A23*A31*A42*A65+A13*A26* &
      A31*A42*A65+A16*A21*A33*A42*A65-A11*A26*A33*A42*A65-A13*A21*A36*A42* &
      A65+A11*A23*A36*A42*A65+A16*A22*A31*A43*A65-A12*A26*A31*A43*A65-A16*A21* &
      A32*A43*A65+A11*A26*A32*A43*A65+A12*A21*A36*A43*A65-A11*A22*A36*A43* &
      A65-A13*A22*A31*A46*A65+A12*A23*A31*A46*A65+A13*A21*A32*A46*A65-A11*A23* &
      A32*A46*A65-A12*A21*A33*A46*A65+A11*A22*A33*A46*A65-A15*A23*A32*A41* &
      A66+A13*A25*A32*A41*A66+A15*A22*A33*A41*A66-A12*A25*A33*A41*A66-A13*A22* &
      A35*A41*A66+A12*A23*A35*A41*A66+A15*A23*A31*A42*A66-A13*A25*A31*A42* &
      A66-A15*A21*A33*A42*A66+A11*A25*A33*A42*A66+A13*A21*A35*A42*A66-A11*A23* &
      A35*A42*A66-A15*A22*A31*A43*A66+A12*A25*A31*A43*A66+A15*A21*A32*A43* &
      A66-A11*A25*A32*A43*A66-A12*A21*A35*A43*A66+A11*A22*A35*A43*A66+A13*A22* &
      A31*A45*A66-A12*A23*A31*A45*A66-A13*A21*A32*A45*A66+A11*A23*A32*A45* &
      A66+A12*A21*A33*A45*A66-A11*A22*A33*A45*A66

      COFACTOR(6,4) = A16*A25*A33*A42*A51-A15*A26* &
      A33*A42*A51-A16*A23*A35*A42*A51+A13*A26*A35*A42*A51+A15*A23*A36*A42* &
      A51-A13*A25*A36*A42*A51-A16*A25*A32*A43*A51+A15*A26*A32*A43*A51+A16*A22* &
      A35*A43*A51-A12*A26*A35*A43*A51-A15*A22*A36*A43*A51+A12*A25*A36*A43* &
      A51+A16*A23*A32*A45*A51-A13*A26*A32*A45*A51-A16*A22*A33*A45*A51+A12*A26* &
      A33*A45*A51+A13*A22*A36*A45*A51-A12*A23*A36*A45*A51-A15*A23*A32*A46* &
      A51+A13*A25*A32*A46*A51+A15*A22*A33*A46*A51-A12*A25*A33*A46*A51-A13*A22* &
      A35*A46*A51+A12*A23*A35*A46*A51-A16*A25*A33*A41*A52+A15*A26*A33*A41* &
      A52+A16*A23*A35*A41*A52-A13*A26*A35*A41*A52-A15*A23*A36*A41*A52+A13*A25* &
      A36*A41*A52+A16*A25*A31*A43*A52-A15*A26*A31*A43*A52-A16*A21*A35*A43* &
      A52+A11*A26*A35*A43*A52+A15*A21*A36*A43*A52-A11*A25*A36*A43*A52-A16*A23* &
      A31*A45*A52+A13*A26*A31*A45*A52+A16*A21*A33*A45*A52-A11*A26*A33*A45* &
      A52-A13*A21*A36*A45*A52+A11*A23*A36*A45*A52+A15*A23*A31*A46*A52-A13*A25* &
      A31*A46*A52-A15*A21*A33*A46*A52+A11*A25*A33*A46*A52+A13*A21*A35*A46* &
      A52-A11*A23*A35*A46*A52+A16*A25*A32*A41*A53-A15*A26*A32*A41*A53-A16*A22* &
      A35*A41*A53+A12*A26*A35*A41*A53+A15*A22*A36*A41*A53-A12*A25*A36*A41* &
      A53-A16*A25*A31*A42*A53+A15*A26*A31*A42*A53+A16*A21*A35*A42*A53-A11*A26* &
      A35*A42*A53-A15*A21*A36*A42*A53+A11*A25*A36*A42*A53+A16*A22*A31*A45* &
      A53-A12*A26*A31*A45*A53-A16*A21*A32*A45*A53+A11*A26*A32*A45*A53+A12*A21* &
      A36*A45*A53-A11*A22*A36*A45*A53-A15*A22*A31*A46*A53+A12*A25*A31*A46* &
      A53+A15*A21*A32*A46*A53-A11*A25*A32*A46*A53-A12*A21*A35*A46*A53+A11*A22* &
      A35*A46*A53-A16*A23*A32*A41*A55+A13*A26*A32*A41*A55+A16*A22*A33*A41* &
      A55-A12*A26*A33*A41*A55-A13*A22*A36*A41*A55+A12*A23*A36*A41*A55+A16*A23* &
      A31*A42*A55-A13*A26*A31*A42*A55-A16*A21*A33*A42*A55+A11*A26*A33*A42* &
      A55+A13*A21*A36*A42*A55-A11*A23*A36*A42*A55-A16*A22*A31*A43*A55+A12*A26* &
      A31*A43*A55+A16*A21*A32*A43*A55-A11*A26*A32*A43*A55-A12*A21*A36*A43* &
      A55+A11*A22*A36*A43*A55+A13*A22*A31*A46*A55-A12*A23*A31*A46*A55-A13*A21* &
      A32*A46*A55+A11*A23*A32*A46*A55+A12*A21*A33*A46*A55-A11*A22*A33*A46* &
      A55+A15*A23*A32*A41*A56-A13*A25*A32*A41*A56-A15*A22*A33*A41*A56+A12*A25* &
      A33*A41*A56+A13*A22*A35*A41*A56-A12*A23*A35*A41*A56-A15*A23*A31*A42* &
      A56+A13*A25*A31*A42*A56+A15*A21*A33*A42*A56-A11*A25*A33*A42*A56-A13*A21* &
      A35*A42*A56+A11*A23*A35*A42*A56+A15*A22*A31*A43*A56-A12*A25*A31*A43* &
      A56-A15*A21*A32*A43*A56+A11*A25*A32*A43*A56+A12*A21*A35*A43*A56-A11*A22* &
      A35*A43*A56-A13*A22*A31*A45*A56+A12*A23*A31*A45*A56+A13*A21*A32*A45* &
      A56-A11*A23*A32*A45*A56-A12*A21*A33*A45*A56+A11*A22*A33*A45*A56 

      COFACTOR(1,5) = A26* &
      A34*A43*A52*A61-A24*A36*A43*A52*A61-A26*A33*A44*A52*A61+A23*A36*A44*A52* &
      A61+A24*A33*A46*A52*A61-A23*A34*A46*A52*A61-A26*A34*A42*A53*A61+A24*A36* &
      A42*A53*A61+A26*A32*A44*A53*A61-A22*A36*A44*A53*A61-A24*A32*A46*A53* &
      A61+A22*A34*A46*A53*A61+A26*A33*A42*A54*A61-A23*A36*A42*A54*A61-A26*A32* &
      A43*A54*A61+A22*A36*A43*A54*A61+A23*A32*A46*A54*A61-A22*A33*A46*A54* &
      A61-A24*A33*A42*A56*A61+A23*A34*A42*A56*A61+A24*A32*A43*A56*A61-A22*A34* &
      A43*A56*A61-A23*A32*A44*A56*A61+A22*A33*A44*A56*A61-A26*A34*A43*A51* &
      A62+A24*A36*A43*A51*A62+A26*A33*A44*A51*A62-A23*A36*A44*A51*A62-A24*A33* &
      A46*A51*A62+A23*A34*A46*A51*A62+A26*A34*A41*A53*A62-A24*A36*A41*A53* &
      A62-A26*A31*A44*A53*A62+A21*A36*A44*A53*A62+A24*A31*A46*A53*A62-A21*A34* &
      A46*A53*A62-A26*A33*A41*A54*A62+A23*A36*A41*A54*A62+A26*A31*A43*A54* &
      A62-A21*A36*A43*A54*A62-A23*A31*A46*A54*A62+A21*A33*A46*A54*A62+A24*A33* &
      A41*A56*A62-A23*A34*A41*A56*A62-A24*A31*A43*A56*A62+A21*A34*A43*A56* &
      A62+A23*A31*A44*A56*A62-A21*A33*A44*A56*A62+A26*A34*A42*A51*A63-A24*A36* &
      A42*A51*A63-A26*A32*A44*A51*A63+A22*A36*A44*A51*A63+A24*A32*A46*A51* &
      A63-A22*A34*A46*A51*A63-A26*A34*A41*A52*A63+A24*A36*A41*A52*A63+A26*A31* &
      A44*A52*A63-A21*A36*A44*A52*A63-A24*A31*A46*A52*A63+A21*A34*A46*A52* &
      A63+A26*A32*A41*A54*A63-A22*A36*A41*A54*A63-A26*A31*A42*A54*A63+A21*A36* &
      A42*A54*A63+A22*A31*A46*A54*A63-A21*A32*A46*A54*A63-A24*A32*A41*A56* &
      A63+A22*A34*A41*A56*A63+A24*A31*A42*A56*A63-A21*A34*A42*A56*A63-A22*A31* &
      A44*A56*A63+A21*A32*A44*A56*A63-A26*A33*A42*A51*A64+A23*A36*A42*A51* &
      A64+A26*A32*A43*A51*A64-A22*A36*A43*A51*A64-A23*A32*A46*A51*A64+A22*A33* &
      A46*A51*A64+A26*A33*A41*A52*A64-A23*A36*A41*A52*A64-A26*A31*A43*A52* &
      A64+A21*A36*A43*A52*A64+A23*A31*A46*A52*A64-A21*A33*A46*A52*A64-A26*A32* &
      A41*A53*A64+A22*A36*A41*A53*A64+A26*A31*A42*A53*A64-A21*A36*A42*A53* &
      A64-A22*A31*A46*A53*A64+A21*A32*A46*A53*A64+A23*A32*A41*A56*A64-A22*A33* &
      A41*A56*A64-A23*A31*A42*A56*A64+A21*A33*A42*A56*A64+A22*A31*A43*A56* &
      A64-A21*A32*A43*A56*A64+A24*A33*A42*A51*A66-A23*A34*A42*A51*A66-A24*A32* &
      A43*A51*A66+A22*A34*A43*A51*A66+A23*A32*A44*A51*A66-A22*A33*A44*A51* &
      A66-A24*A33*A41*A52*A66+A23*A34*A41*A52*A66+A24*A31*A43*A52*A66-A21*A34* &
      A43*A52*A66-A23*A31*A44*A52*A66+A21*A33*A44*A52*A66+A24*A32*A41*A53* &
      A66-A22*A34*A41*A53*A66-A24*A31*A42*A53*A66+A21*A34*A42*A53*A66+A22*A31* &
      A44*A53*A66-A21*A32*A44*A53*A66-A23*A32*A41*A54*A66+A22*A33*A41*A54* &
      A66+A23*A31*A42*A54*A66-A21*A33*A42*A54*A66-A22*A31*A43*A54*A66+A21*A32* &
      A43*A54*A66 

      COFACTOR(2,5) = -A16*A34*A43*A52*A61+A14*A36*A43*A52*A61+A16*A33*A44*A52* &
      A61-A13*A36*A44*A52*A61-A14*A33*A46*A52*A61+A13*A34*A46*A52*A61+A16*A34* &
      A42*A53*A61-A14*A36*A42*A53*A61-A16*A32*A44*A53*A61+A12*A36*A44*A53* &
      A61+A14*A32*A46*A53*A61-A12*A34*A46*A53*A61-A16*A33*A42*A54*A61+A13*A36* &
      A42*A54*A61+A16*A32*A43*A54*A61-A12*A36*A43*A54*A61-A13*A32*A46*A54* &
      A61+A12*A33*A46*A54*A61+A14*A33*A42*A56*A61-A13*A34*A42*A56*A61-A14*A32* &
      A43*A56*A61+A12*A34*A43*A56*A61+A13*A32*A44*A56*A61-A12*A33*A44*A56* &
      A61+A16*A34*A43*A51*A62-A14*A36*A43*A51*A62-A16*A33*A44*A51*A62+A13*A36* &
      A44*A51*A62+A14*A33*A46*A51*A62-A13*A34*A46*A51*A62-A16*A34*A41*A53* &
      A62+A14*A36*A41*A53*A62+A16*A31*A44*A53*A62-A11*A36*A44*A53*A62-A14*A31* &
      A46*A53*A62+A11*A34*A46*A53*A62+A16*A33*A41*A54*A62-A13*A36*A41*A54* &
      A62-A16*A31*A43*A54*A62+A11*A36*A43*A54*A62+A13*A31*A46*A54*A62-A11*A33* &
      A46*A54*A62-A14*A33*A41*A56*A62+A13*A34*A41*A56*A62+A14*A31*A43*A56* &
      A62-A11*A34*A43*A56*A62-A13*A31*A44*A56*A62+A11*A33*A44*A56*A62-A16*A34* &
      A42*A51*A63+A14*A36*A42*A51*A63+A16*A32*A44*A51*A63-A12*A36*A44*A51* &
      A63-A14*A32*A46*A51*A63+A12*A34*A46*A51*A63+A16*A34*A41*A52*A63-A14*A36* &
      A41*A52*A63-A16*A31*A44*A52*A63+A11*A36*A44*A52*A63+A14*A31*A46*A52* &
      A63-A11*A34*A46*A52*A63-A16*A32*A41*A54*A63+A12*A36*A41*A54*A63+A16*A31* &
      A42*A54*A63-A11*A36*A42*A54*A63-A12*A31*A46*A54*A63+A11*A32*A46*A54* &
      A63+A14*A32*A41*A56*A63-A12*A34*A41*A56*A63-A14*A31*A42*A56*A63+A11*A34* &
      A42*A56*A63+A12*A31*A44*A56*A63-A11*A32*A44*A56*A63+A16*A33*A42*A51* &
      A64-A13*A36*A42*A51*A64-A16*A32*A43*A51*A64+A12*A36*A43*A51*A64+A13*A32* &
      A46*A51*A64-A12*A33*A46*A51*A64-A16*A33*A41*A52*A64+A13*A36*A41*A52* &
      A64+A16*A31*A43*A52*A64-A11*A36*A43*A52*A64-A13*A31*A46*A52*A64+A11*A33* &
      A46*A52*A64+A16*A32*A41*A53*A64-A12*A36*A41*A53*A64-A16*A31*A42*A53* &
      A64+A11*A36*A42*A53*A64+A12*A31*A46*A53*A64-A11*A32*A46*A53*A64-A13*A32* &
      A41*A56*A64+A12*A33*A41*A56*A64+A13*A31*A42*A56*A64-A11*A33*A42*A56* &
      A64-A12*A31*A43*A56*A64+A11*A32*A43*A56*A64-A14*A33*A42*A51*A66+A13*A34* &
      A42*A51*A66+A14*A32*A43*A51*A66-A12*A34*A43*A51*A66-A13*A32*A44*A51* &
      A66+A12*A33*A44*A51*A66+A14*A33*A41*A52*A66-A13*A34*A41*A52*A66-A14*A31* &
      A43*A52*A66+A11*A34*A43*A52*A66+A13*A31*A44*A52*A66-A11*A33*A44*A52* &
      A66-A14*A32*A41*A53*A66+A12*A34*A41*A53*A66+A14*A31*A42*A53*A66-A11*A34* &
      A42*A53*A66-A12*A31*A44*A53*A66+A11*A32*A44*A53*A66+A13*A32*A41*A54* &
      A66-A12*A33*A41*A54*A66-A13*A31*A42*A54*A66+A11*A33*A42*A54*A66+A12*A31* &
      A43*A54*A66-A11*A32*A43*A54*A66

      COFACTOR(3,5) = A16*A24*A43*A52*A61-A14*A26*A43*A52* &
      A61-A16*A23*A44*A52*A61+A13*A26*A44*A52*A61+A14*A23*A46*A52*A61-A13*A24* &
      A46*A52*A61-A16*A24*A42*A53*A61+A14*A26*A42*A53*A61+A16*A22*A44*A53* &
      A61-A12*A26*A44*A53*A61-A14*A22*A46*A53*A61+A12*A24*A46*A53*A61+A16*A23* &
      A42*A54*A61-A13*A26*A42*A54*A61-A16*A22*A43*A54*A61+A12*A26*A43*A54* &
      A61+A13*A22*A46*A54*A61-A12*A23*A46*A54*A61-A14*A23*A42*A56*A61+A13*A24* &
      A42*A56*A61+A14*A22*A43*A56*A61-A12*A24*A43*A56*A61-A13*A22*A44*A56* &
      A61+A12*A23*A44*A56*A61-A16*A24*A43*A51*A62+A14*A26*A43*A51*A62+A16*A23* &
      A44*A51*A62-A13*A26*A44*A51*A62-A14*A23*A46*A51*A62+A13*A24*A46*A51* &
      A62+A16*A24*A41*A53*A62-A14*A26*A41*A53*A62-A16*A21*A44*A53*A62+A11*A26* &
      A44*A53*A62+A14*A21*A46*A53*A62-A11*A24*A46*A53*A62-A16*A23*A41*A54* &
      A62+A13*A26*A41*A54*A62+A16*A21*A43*A54*A62-A11*A26*A43*A54*A62-A13*A21* &
      A46*A54*A62+A11*A23*A46*A54*A62+A14*A23*A41*A56*A62-A13*A24*A41*A56* &
      A62-A14*A21*A43*A56*A62+A11*A24*A43*A56*A62+A13*A21*A44*A56*A62-A11*A23* &
      A44*A56*A62+A16*A24*A42*A51*A63-A14*A26*A42*A51*A63-A16*A22*A44*A51* &
      A63+A12*A26*A44*A51*A63+A14*A22*A46*A51*A63-A12*A24*A46*A51*A63-A16*A24* &
      A41*A52*A63+A14*A26*A41*A52*A63+A16*A21*A44*A52*A63-A11*A26*A44*A52* &
      A63-A14*A21*A46*A52*A63+A11*A24*A46*A52*A63+A16*A22*A41*A54*A63-A12*A26* &
      A41*A54*A63-A16*A21*A42*A54*A63+A11*A26*A42*A54*A63+A12*A21*A46*A54* &
      A63-A11*A22*A46*A54*A63-A14*A22*A41*A56*A63+A12*A24*A41*A56*A63+A14*A21* &
      A42*A56*A63-A11*A24*A42*A56*A63-A12*A21*A44*A56*A63+A11*A22*A44*A56* &
      A63-A16*A23*A42*A51*A64+A13*A26*A42*A51*A64+A16*A22*A43*A51*A64-A12*A26* &
      A43*A51*A64-A13*A22*A46*A51*A64+A12*A23*A46*A51*A64+A16*A23*A41*A52* &
      A64-A13*A26*A41*A52*A64-A16*A21*A43*A52*A64+A11*A26*A43*A52*A64+A13*A21* &
      A46*A52*A64-A11*A23*A46*A52*A64-A16*A22*A41*A53*A64+A12*A26*A41*A53* &
      A64+A16*A21*A42*A53*A64-A11*A26*A42*A53*A64-A12*A21*A46*A53*A64+A11*A22* &
      A46*A53*A64+A13*A22*A41*A56*A64-A12*A23*A41*A56*A64-A13*A21*A42*A56* &
      A64+A11*A23*A42*A56*A64+A12*A21*A43*A56*A64-A11*A22*A43*A56*A64+A14*A23* &
      A42*A51*A66-A13*A24*A42*A51*A66-A14*A22*A43*A51*A66+A12*A24*A43*A51* &
      A66+A13*A22*A44*A51*A66-A12*A23*A44*A51*A66-A14*A23*A41*A52*A66+A13*A24* &
      A41*A52*A66+A14*A21*A43*A52*A66-A11*A24*A43*A52*A66-A13*A21*A44*A52* &
      A66+A11*A23*A44*A52*A66+A14*A22*A41*A53*A66-A12*A24*A41*A53*A66-A14*A21* &
      A42*A53*A66+A11*A24*A42*A53*A66+A12*A21*A44*A53*A66-A11*A22*A44*A53* &
      A66-A13*A22*A41*A54*A66+A12*A23*A41*A54*A66+A13*A21*A42*A54*A66-A11*A23* &
      A42*A54*A66-A12*A21*A43*A54*A66+A11*A22*A43*A54*A66

      COFACTOR(4,5) = -A16*A24*A33*A52* &
      A61+A14*A26*A33*A52*A61+A16*A23*A34*A52*A61-A13*A26*A34*A52*A61-A14*A23* &
      A36*A52*A61+A13*A24*A36*A52*A61+A16*A24*A32*A53*A61-A14*A26*A32*A53* &
      A61-A16*A22*A34*A53*A61+A12*A26*A34*A53*A61+A14*A22*A36*A53*A61-A12*A24* &
      A36*A53*A61-A16*A23*A32*A54*A61+A13*A26*A32*A54*A61+A16*A22*A33*A54* &
      A61-A12*A26*A33*A54*A61-A13*A22*A36*A54*A61+A12*A23*A36*A54*A61+A14*A23* &
      A32*A56*A61-A13*A24*A32*A56*A61-A14*A22*A33*A56*A61+A12*A24*A33*A56* &
      A61+A13*A22*A34*A56*A61-A12*A23*A34*A56*A61+A16*A24*A33*A51*A62-A14*A26* &
      A33*A51*A62-A16*A23*A34*A51*A62+A13*A26*A34*A51*A62+A14*A23*A36*A51* &
      A62-A13*A24*A36*A51*A62-A16*A24*A31*A53*A62+A14*A26*A31*A53*A62+A16*A21* &
      A34*A53*A62-A11*A26*A34*A53*A62-A14*A21*A36*A53*A62+A11*A24*A36*A53* &
      A62+A16*A23*A31*A54*A62-A13*A26*A31*A54*A62-A16*A21*A33*A54*A62+A11*A26* &
      A33*A54*A62+A13*A21*A36*A54*A62-A11*A23*A36*A54*A62-A14*A23*A31*A56* &
      A62+A13*A24*A31*A56*A62+A14*A21*A33*A56*A62-A11*A24*A33*A56*A62-A13*A21* &
      A34*A56*A62+A11*A23*A34*A56*A62-A16*A24*A32*A51*A63+A14*A26*A32*A51* &
      A63+A16*A22*A34*A51*A63-A12*A26*A34*A51*A63-A14*A22*A36*A51*A63+A12*A24* &
      A36*A51*A63+A16*A24*A31*A52*A63-A14*A26*A31*A52*A63-A16*A21*A34*A52* &
      A63+A11*A26*A34*A52*A63+A14*A21*A36*A52*A63-A11*A24*A36*A52*A63-A16*A22* &
      A31*A54*A63+A12*A26*A31*A54*A63+A16*A21*A32*A54*A63-A11*A26*A32*A54* &
      A63-A12*A21*A36*A54*A63+A11*A22*A36*A54*A63+A14*A22*A31*A56*A63-A12*A24* &
      A31*A56*A63-A14*A21*A32*A56*A63+A11*A24*A32*A56*A63+A12*A21*A34*A56* &
      A63-A11*A22*A34*A56*A63+A16*A23*A32*A51*A64-A13*A26*A32*A51*A64-A16*A22* &
      A33*A51*A64+A12*A26*A33*A51*A64+A13*A22*A36*A51*A64-A12*A23*A36*A51* &
      A64-A16*A23*A31*A52*A64+A13*A26*A31*A52*A64+A16*A21*A33*A52*A64-A11*A26* &
      A33*A52*A64-A13*A21*A36*A52*A64+A11*A23*A36*A52*A64+A16*A22*A31*A53* &
      A64-A12*A26*A31*A53*A64-A16*A21*A32*A53*A64+A11*A26*A32*A53*A64+A12*A21* &
      A36*A53*A64-A11*A22*A36*A53*A64-A13*A22*A31*A56*A64+A12*A23*A31*A56* &
      A64+A13*A21*A32*A56*A64-A11*A23*A32*A56*A64-A12*A21*A33*A56*A64+A11*A22* &
      A33*A56*A64-A14*A23*A32*A51*A66+A13*A24*A32*A51*A66+A14*A22*A33*A51* &
      A66-A12*A24*A33*A51*A66-A13*A22*A34*A51*A66+A12*A23*A34*A51*A66+A14*A23* &
      A31*A52*A66-A13*A24*A31*A52*A66-A14*A21*A33*A52*A66+A11*A24*A33*A52* &
      A66+A13*A21*A34*A52*A66-A11*A23*A34*A52*A66-A14*A22*A31*A53*A66+A12*A24* &
      A31*A53*A66+A14*A21*A32*A53*A66-A11*A24*A32*A53*A66-A12*A21*A34*A53* &
      A66+A11*A22*A34*A53*A66+A13*A22*A31*A54*A66-A12*A23*A31*A54*A66-A13*A21* &
      A32*A54*A66+A11*A23*A32*A54*A66+A12*A21*A33*A54*A66-A11*A22*A33*A54* &
      A66

      COFACTOR(5,5) = A16*A24*A33*A42*A61-A14*A26*A33*A42*A61-A16*A23*A34*A42*A61+A13*A26*     &
      A34*A42*A61+A14*A23*A36*A42*A61-A13*A24*A36*A42*A61-A16*A24*A32*A43* &
      A61+A14*A26*A32*A43*A61+A16*A22*A34*A43*A61-A12*A26*A34*A43*A61-A14*A22* &
      A36*A43*A61+A12*A24*A36*A43*A61+A16*A23*A32*A44*A61-A13*A26*A32*A44* &
      A61-A16*A22*A33*A44*A61+A12*A26*A33*A44*A61+A13*A22*A36*A44*A61-A12*A23* &
      A36*A44*A61-A14*A23*A32*A46*A61+A13*A24*A32*A46*A61+A14*A22*A33*A46* &
      A61-A12*A24*A33*A46*A61-A13*A22*A34*A46*A61+A12*A23*A34*A46*A61-A16*A24* &
      A33*A41*A62+A14*A26*A33*A41*A62+A16*A23*A34*A41*A62-A13*A26*A34*A41* &
      A62-A14*A23*A36*A41*A62+A13*A24*A36*A41*A62+A16*A24*A31*A43*A62-A14*A26* &
      A31*A43*A62-A16*A21*A34*A43*A62+A11*A26*A34*A43*A62+A14*A21*A36*A43* &
      A62-A11*A24*A36*A43*A62-A16*A23*A31*A44*A62+A13*A26*A31*A44*A62+A16*A21* &
      A33*A44*A62-A11*A26*A33*A44*A62-A13*A21*A36*A44*A62+A11*A23*A36*A44* &
      A62+A14*A23*A31*A46*A62-A13*A24*A31*A46*A62-A14*A21*A33*A46*A62+A11*A24* &
      A33*A46*A62+A13*A21*A34*A46*A62-A11*A23*A34*A46*A62+A16*A24*A32*A41* &
      A63-A14*A26*A32*A41*A63-A16*A22*A34*A41*A63+A12*A26*A34*A41*A63+A14*A22* &
      A36*A41*A63-A12*A24*A36*A41*A63-A16*A24*A31*A42*A63+A14*A26*A31*A42* &
      A63+A16*A21*A34*A42*A63-A11*A26*A34*A42*A63-A14*A21*A36*A42*A63+A11*A24* &
      A36*A42*A63+A16*A22*A31*A44*A63-A12*A26*A31*A44*A63-A16*A21*A32*A44* &
      A63+A11*A26*A32*A44*A63+A12*A21*A36*A44*A63-A11*A22*A36*A44*A63-A14*A22* &
      A31*A46*A63+A12*A24*A31*A46*A63+A14*A21*A32*A46*A63-A11*A24*A32*A46* &
      A63-A12*A21*A34*A46*A63+A11*A22*A34*A46*A63-A16*A23*A32*A41*A64+A13*A26* &
      A32*A41*A64+A16*A22*A33*A41*A64-A12*A26*A33*A41*A64-A13*A22*A36*A41* &
      A64+A12*A23*A36*A41*A64+A16*A23*A31*A42*A64-A13*A26*A31*A42*A64-A16*A21* &
      A33*A42*A64+A11*A26*A33*A42*A64+A13*A21*A36*A42*A64-A11*A23*A36*A42* &
      A64-A16*A22*A31*A43*A64+A12*A26*A31*A43*A64+A16*A21*A32*A43*A64-A11*A26* &
      A32*A43*A64-A12*A21*A36*A43*A64+A11*A22*A36*A43*A64+A13*A22*A31*A46* &
      A64-A12*A23*A31*A46*A64-A13*A21*A32*A46*A64+A11*A23*A32*A46*A64+A12*A21* &
      A33*A46*A64-A11*A22*A33*A46*A64+A14*A23*A32*A41*A66-A13*A24*A32*A41* &
      A66-A14*A22*A33*A41*A66+A12*A24*A33*A41*A66+A13*A22*A34*A41*A66-A12*A23* &
      A34*A41*A66-A14*A23*A31*A42*A66+A13*A24*A31*A42*A66+A14*A21*A33*A42* &
      A66-A11*A24*A33*A42*A66-A13*A21*A34*A42*A66+A11*A23*A34*A42*A66+A14*A22* &
      A31*A43*A66-A12*A24*A31*A43*A66-A14*A21*A32*A43*A66+A11*A24*A32*A43* &
      A66+A12*A21*A34*A43*A66-A11*A22*A34*A43*A66-A13*A22*A31*A44*A66+A12*A23* &
      A31*A44*A66+A13*A21*A32*A44*A66-A11*A23*A32*A44*A66-A12*A21*A33*A44* &
      A66+A11*A22*A33*A44*A66

      COFACTOR(6,5) = -A16*A24*A33*A42*A51+A14*A26*A33*A42*A51+A16*A23* &
      A34*A42*A51-A13*A26*A34*A42*A51-A14*A23*A36*A42*A51+A13*A24*A36*A42* &
      A51+A16*A24*A32*A43*A51-A14*A26*A32*A43*A51-A16*A22*A34*A43*A51+A12*A26* &
      A34*A43*A51+A14*A22*A36*A43*A51-A12*A24*A36*A43*A51-A16*A23*A32*A44* &
      A51+A13*A26*A32*A44*A51+A16*A22*A33*A44*A51-A12*A26*A33*A44*A51-A13*A22* &
      A36*A44*A51+A12*A23*A36*A44*A51+A14*A23*A32*A46*A51-A13*A24*A32*A46* &
      A51-A14*A22*A33*A46*A51+A12*A24*A33*A46*A51+A13*A22*A34*A46*A51-A12*A23* &
      A34*A46*A51+A16*A24*A33*A41*A52-A14*A26*A33*A41*A52-A16*A23*A34*A41* &
      A52+A13*A26*A34*A41*A52+A14*A23*A36*A41*A52-A13*A24*A36*A41*A52-A16*A24* &
      A31*A43*A52+A14*A26*A31*A43*A52+A16*A21*A34*A43*A52-A11*A26*A34*A43* &
      A52-A14*A21*A36*A43*A52+A11*A24*A36*A43*A52+A16*A23*A31*A44*A52-A13*A26* &
      A31*A44*A52-A16*A21*A33*A44*A52+A11*A26*A33*A44*A52+A13*A21*A36*A44* &
      A52-A11*A23*A36*A44*A52-A14*A23*A31*A46*A52+A13*A24*A31*A46*A52+A14*A21* &
      A33*A46*A52-A11*A24*A33*A46*A52-A13*A21*A34*A46*A52+A11*A23*A34*A46* &
      A52-A16*A24*A32*A41*A53+A14*A26*A32*A41*A53+A16*A22*A34*A41*A53-A12*A26* &
      A34*A41*A53-A14*A22*A36*A41*A53+A12*A24*A36*A41*A53+A16*A24*A31*A42* &
      A53-A14*A26*A31*A42*A53-A16*A21*A34*A42*A53+A11*A26*A34*A42*A53+A14*A21* &
      A36*A42*A53-A11*A24*A36*A42*A53-A16*A22*A31*A44*A53+A12*A26*A31*A44* &
      A53+A16*A21*A32*A44*A53-A11*A26*A32*A44*A53-A12*A21*A36*A44*A53+A11*A22* &
      A36*A44*A53+A14*A22*A31*A46*A53-A12*A24*A31*A46*A53-A14*A21*A32*A46* &
      A53+A11*A24*A32*A46*A53+A12*A21*A34*A46*A53-A11*A22*A34*A46*A53+A16*A23* &
      A32*A41*A54-A13*A26*A32*A41*A54-A16*A22*A33*A41*A54+A12*A26*A33*A41* &
      A54+A13*A22*A36*A41*A54-A12*A23*A36*A41*A54-A16*A23*A31*A42*A54+A13*A26* &
      A31*A42*A54+A16*A21*A33*A42*A54-A11*A26*A33*A42*A54-A13*A21*A36*A42* &
      A54+A11*A23*A36*A42*A54+A16*A22*A31*A43*A54-A12*A26*A31*A43*A54-A16*A21* &
      A32*A43*A54+A11*A26*A32*A43*A54+A12*A21*A36*A43*A54-A11*A22*A36*A43* &
      A54-A13*A22*A31*A46*A54+A12*A23*A31*A46*A54+A13*A21*A32*A46*A54-A11*A23* &
      A32*A46*A54-A12*A21*A33*A46*A54+A11*A22*A33*A46*A54-A14*A23*A32*A41* &
      A56+A13*A24*A32*A41*A56+A14*A22*A33*A41*A56-A12*A24*A33*A41*A56-A13*A22* &
      A34*A41*A56+A12*A23*A34*A41*A56+A14*A23*A31*A42*A56-A13*A24*A31*A42* &
      A56-A14*A21*A33*A42*A56+A11*A24*A33*A42*A56+A13*A21*A34*A42*A56-A11*A23* &
      A34*A42*A56-A14*A22*A31*A43*A56+A12*A24*A31*A43*A56+A14*A21*A32*A43* &
      A56-A11*A24*A32*A43*A56-A12*A21*A34*A43*A56+A11*A22*A34*A43*A56+A13*A22* &
      A31*A44*A56-A12*A23*A31*A44*A56-A13*A21*A32*A44*A56+A11*A23*A32*A44* &
      A56+A12*A21*A33*A44*A56-A11*A22*A33*A44*A56

      COFACTOR(1,6) = -A25*A34*A43*A52*A61+A24* &
      A35*A43*A52*A61+A25*A33*A44*A52*A61-A23*A35*A44*A52*A61-A24*A33*A45*A52* &
      A61+A23*A34*A45*A52*A61+A25*A34*A42*A53*A61-A24*A35*A42*A53*A61-A25*A32* &
      A44*A53*A61+A22*A35*A44*A53*A61+A24*A32*A45*A53*A61-A22*A34*A45*A53* &
      A61-A25*A33*A42*A54*A61+A23*A35*A42*A54*A61+A25*A32*A43*A54*A61-A22*A35* &
      A43*A54*A61-A23*A32*A45*A54*A61+A22*A33*A45*A54*A61+A24*A33*A42*A55* &
      A61-A23*A34*A42*A55*A61-A24*A32*A43*A55*A61+A22*A34*A43*A55*A61+A23*A32* &
      A44*A55*A61-A22*A33*A44*A55*A61+A25*A34*A43*A51*A62-A24*A35*A43*A51* &
      A62-A25*A33*A44*A51*A62+A23*A35*A44*A51*A62+A24*A33*A45*A51*A62-A23*A34* &
      A45*A51*A62-A25*A34*A41*A53*A62+A24*A35*A41*A53*A62+A25*A31*A44*A53* &
      A62-A21*A35*A44*A53*A62-A24*A31*A45*A53*A62+A21*A34*A45*A53*A62+A25*A33* &
      A41*A54*A62-A23*A35*A41*A54*A62-A25*A31*A43*A54*A62+A21*A35*A43*A54* &
      A62+A23*A31*A45*A54*A62-A21*A33*A45*A54*A62-A24*A33*A41*A55*A62+A23*A34* &
      A41*A55*A62+A24*A31*A43*A55*A62-A21*A34*A43*A55*A62-A23*A31*A44*A55* &
      A62+A21*A33*A44*A55*A62-A25*A34*A42*A51*A63+A24*A35*A42*A51*A63+A25*A32* &
      A44*A51*A63-A22*A35*A44*A51*A63-A24*A32*A45*A51*A63+A22*A34*A45*A51* &
      A63+A25*A34*A41*A52*A63-A24*A35*A41*A52*A63-A25*A31*A44*A52*A63+A21*A35* &
      A44*A52*A63+A24*A31*A45*A52*A63-A21*A34*A45*A52*A63-A25*A32*A41*A54* &
      A63+A22*A35*A41*A54*A63+A25*A31*A42*A54*A63-A21*A35*A42*A54*A63-A22*A31* &
      A45*A54*A63+A21*A32*A45*A54*A63+A24*A32*A41*A55*A63-A22*A34*A41*A55* &
      A63-A24*A31*A42*A55*A63+A21*A34*A42*A55*A63+A22*A31*A44*A55*A63-A21*A32* &
      A44*A55*A63+A25*A33*A42*A51*A64-A23*A35*A42*A51*A64-A25*A32*A43*A51* &
      A64+A22*A35*A43*A51*A64+A23*A32*A45*A51*A64-A22*A33*A45*A51*A64-A25*A33* &
      A41*A52*A64+A23*A35*A41*A52*A64+A25*A31*A43*A52*A64-A21*A35*A43*A52* &
      A64-A23*A31*A45*A52*A64+A21*A33*A45*A52*A64+A25*A32*A41*A53*A64-A22*A35* &
      A41*A53*A64-A25*A31*A42*A53*A64+A21*A35*A42*A53*A64+A22*A31*A45*A53* &
      A64-A21*A32*A45*A53*A64-A23*A32*A41*A55*A64+A22*A33*A41*A55*A64+A23*A31* &
      A42*A55*A64-A21*A33*A42*A55*A64-A22*A31*A43*A55*A64+A21*A32*A43*A55* &
      A64-A24*A33*A42*A51*A65+A23*A34*A42*A51*A65+A24*A32*A43*A51*A65-A22*A34* &
      A43*A51*A65-A23*A32*A44*A51*A65+A22*A33*A44*A51*A65+A24*A33*A41*A52* &
      A65-A23*A34*A41*A52*A65-A24*A31*A43*A52*A65+A21*A34*A43*A52*A65+A23*A31* &
      A44*A52*A65-A21*A33*A44*A52*A65-A24*A32*A41*A53*A65+A22*A34*A41*A53* &
      A65+A24*A31*A42*A53*A65-A21*A34*A42*A53*A65-A22*A31*A44*A53*A65+A21*A32* &
      A44*A53*A65+A23*A32*A41*A54*A65-A22*A33*A41*A54*A65-A23*A31*A42*A54* &
      A65+A21*A33*A42*A54*A65+A22*A31*A43*A54*A65-A21*A32*A43*A54*A65

      COFACTOR(2,6) = A15*A34* &
      A43*A52*A61-A14*A35*A43*A52*A61-A15*A33*A44*A52*A61+A13*A35*A44*A52* &
      A61+A14*A33*A45*A52*A61-A13*A34*A45*A52*A61-A15*A34*A42*A53*A61+A14*A35* &
      A42*A53*A61+A15*A32*A44*A53*A61-A12*A35*A44*A53*A61-A14*A32*A45*A53* &
      A61+A12*A34*A45*A53*A61+A15*A33*A42*A54*A61-A13*A35*A42*A54*A61-A15*A32* &
      A43*A54*A61+A12*A35*A43*A54*A61+A13*A32*A45*A54*A61-A12*A33*A45*A54* &
      A61-A14*A33*A42*A55*A61+A13*A34*A42*A55*A61+A14*A32*A43*A55*A61-A12*A34* &
      A43*A55*A61-A13*A32*A44*A55*A61+A12*A33*A44*A55*A61-A15*A34*A43*A51* &
      A62+A14*A35*A43*A51*A62+A15*A33*A44*A51*A62-A13*A35*A44*A51*A62-A14*A33* &
      A45*A51*A62+A13*A34*A45*A51*A62+A15*A34*A41*A53*A62-A14*A35*A41*A53* &
      A62-A15*A31*A44*A53*A62+A11*A35*A44*A53*A62+A14*A31*A45*A53*A62-A11*A34* &
      A45*A53*A62-A15*A33*A41*A54*A62+A13*A35*A41*A54*A62+A15*A31*A43*A54* &
      A62-A11*A35*A43*A54*A62-A13*A31*A45*A54*A62+A11*A33*A45*A54*A62+A14*A33* &
      A41*A55*A62-A13*A34*A41*A55*A62-A14*A31*A43*A55*A62+A11*A34*A43*A55* &
      A62+A13*A31*A44*A55*A62-A11*A33*A44*A55*A62+A15*A34*A42*A51*A63-A14*A35* &
      A42*A51*A63-A15*A32*A44*A51*A63+A12*A35*A44*A51*A63+A14*A32*A45*A51* &
      A63-A12*A34*A45*A51*A63-A15*A34*A41*A52*A63+A14*A35*A41*A52*A63+A15*A31* &
      A44*A52*A63-A11*A35*A44*A52*A63-A14*A31*A45*A52*A63+A11*A34*A45*A52* &
      A63+A15*A32*A41*A54*A63-A12*A35*A41*A54*A63-A15*A31*A42*A54*A63+A11*A35* &
      A42*A54*A63+A12*A31*A45*A54*A63-A11*A32*A45*A54*A63-A14*A32*A41*A55* &
      A63+A12*A34*A41*A55*A63+A14*A31*A42*A55*A63-A11*A34*A42*A55*A63-A12*A31* &
      A44*A55*A63+A11*A32*A44*A55*A63-A15*A33*A42*A51*A64+A13*A35*A42*A51* &
      A64+A15*A32*A43*A51*A64-A12*A35*A43*A51*A64-A13*A32*A45*A51*A64+A12*A33* &
      A45*A51*A64+A15*A33*A41*A52*A64-A13*A35*A41*A52*A64-A15*A31*A43*A52* &
      A64+A11*A35*A43*A52*A64+A13*A31*A45*A52*A64-A11*A33*A45*A52*A64-A15*A32* &
      A41*A53*A64+A12*A35*A41*A53*A64+A15*A31*A42*A53*A64-A11*A35*A42*A53* &
      A64-A12*A31*A45*A53*A64+A11*A32*A45*A53*A64+A13*A32*A41*A55*A64-A12*A33* &
      A41*A55*A64-A13*A31*A42*A55*A64+A11*A33*A42*A55*A64+A12*A31*A43*A55* &
      A64-A11*A32*A43*A55*A64+A14*A33*A42*A51*A65-A13*A34*A42*A51*A65-A14*A32* &
      A43*A51*A65+A12*A34*A43*A51*A65+A13*A32*A44*A51*A65-A12*A33*A44*A51* &
      A65-A14*A33*A41*A52*A65+A13*A34*A41*A52*A65+A14*A31*A43*A52*A65-A11*A34* &
      A43*A52*A65-A13*A31*A44*A52*A65+A11*A33*A44*A52*A65+A14*A32*A41*A53* &
      A65-A12*A34*A41*A53*A65-A14*A31*A42*A53*A65+A11*A34*A42*A53*A65+A12*A31* &
      A44*A53*A65-A11*A32*A44*A53*A65-A13*A32*A41*A54*A65+A12*A33*A41*A54* &
      A65+A13*A31*A42*A54*A65-A11*A33*A42*A54*A65-A12*A31*A43*A54*A65+A11*A32* &
      A43*A54*A65

      COFACTOR(3,6) = -A15*A24*A43*A52*A61+A14*A25*A43*A52*A61+A15*A23*A44*A52* &
      A61-A13*A25*A44*A52*A61-A14*A23*A45*A52*A61+A13*A24*A45*A52*A61+A15*A24* &
      A42*A53*A61-A14*A25*A42*A53*A61-A15*A22*A44*A53*A61+A12*A25*A44*A53* &
      A61+A14*A22*A45*A53*A61-A12*A24*A45*A53*A61-A15*A23*A42*A54*A61+A13*A25* &
      A42*A54*A61+A15*A22*A43*A54*A61-A12*A25*A43*A54*A61-A13*A22*A45*A54* &
      A61+A12*A23*A45*A54*A61+A14*A23*A42*A55*A61-A13*A24*A42*A55*A61-A14*A22* &
      A43*A55*A61+A12*A24*A43*A55*A61+A13*A22*A44*A55*A61-A12*A23*A44*A55* &
      A61+A15*A24*A43*A51*A62-A14*A25*A43*A51*A62-A15*A23*A44*A51*A62+A13*A25* &
      A44*A51*A62+A14*A23*A45*A51*A62-A13*A24*A45*A51*A62-A15*A24*A41*A53* &
      A62+A14*A25*A41*A53*A62+A15*A21*A44*A53*A62-A11*A25*A44*A53*A62-A14*A21* &
      A45*A53*A62+A11*A24*A45*A53*A62+A15*A23*A41*A54*A62-A13*A25*A41*A54* &
      A62-A15*A21*A43*A54*A62+A11*A25*A43*A54*A62+A13*A21*A45*A54*A62-A11*A23* &
      A45*A54*A62-A14*A23*A41*A55*A62+A13*A24*A41*A55*A62+A14*A21*A43*A55* &
      A62-A11*A24*A43*A55*A62-A13*A21*A44*A55*A62+A11*A23*A44*A55*A62-A15*A24* &
      A42*A51*A63+A14*A25*A42*A51*A63+A15*A22*A44*A51*A63-A12*A25*A44*A51* &
      A63-A14*A22*A45*A51*A63+A12*A24*A45*A51*A63+A15*A24*A41*A52*A63-A14*A25* &
      A41*A52*A63-A15*A21*A44*A52*A63+A11*A25*A44*A52*A63+A14*A21*A45*A52* &
      A63-A11*A24*A45*A52*A63-A15*A22*A41*A54*A63+A12*A25*A41*A54*A63+A15*A21* &
      A42*A54*A63-A11*A25*A42*A54*A63-A12*A21*A45*A54*A63+A11*A22*A45*A54* &
      A63+A14*A22*A41*A55*A63-A12*A24*A41*A55*A63-A14*A21*A42*A55*A63+A11*A24* &
      A42*A55*A63+A12*A21*A44*A55*A63-A11*A22*A44*A55*A63+A15*A23*A42*A51* &
      A64-A13*A25*A42*A51*A64-A15*A22*A43*A51*A64+A12*A25*A43*A51*A64+A13*A22* &
      A45*A51*A64-A12*A23*A45*A51*A64-A15*A23*A41*A52*A64+A13*A25*A41*A52* &
      A64+A15*A21*A43*A52*A64-A11*A25*A43*A52*A64-A13*A21*A45*A52*A64+A11*A23* &
      A45*A52*A64+A15*A22*A41*A53*A64-A12*A25*A41*A53*A64-A15*A21*A42*A53* &
      A64+A11*A25*A42*A53*A64+A12*A21*A45*A53*A64-A11*A22*A45*A53*A64-A13*A22* &
      A41*A55*A64+A12*A23*A41*A55*A64+A13*A21*A42*A55*A64-A11*A23*A42*A55* &
      A64-A12*A21*A43*A55*A64+A11*A22*A43*A55*A64-A14*A23*A42*A51*A65+A13*A24* &
      A42*A51*A65+A14*A22*A43*A51*A65-A12*A24*A43*A51*A65-A13*A22*A44*A51* &
      A65+A12*A23*A44*A51*A65+A14*A23*A41*A52*A65-A13*A24*A41*A52*A65-A14*A21* &
      A43*A52*A65+A11*A24*A43*A52*A65+A13*A21*A44*A52*A65-A11*A23*A44*A52* &
      A65-A14*A22*A41*A53*A65+A12*A24*A41*A53*A65+A14*A21*A42*A53*A65-A11*A24* &
      A42*A53*A65-A12*A21*A44*A53*A65+A11*A22*A44*A53*A65+A13*A22*A41*A54* &
      A65-A12*A23*A41*A54*A65-A13*A21*A42*A54*A65+A11*A23*A42*A54*A65+A12*A21* &
      A43*A54*A65-A11*A22*A43*A54*A65

      COFACTOR(4,6) = A15*A24*A33*A52*A61-A14*A25*A33*A52* &
      A61-A15*A23*A34*A52*A61+A13*A25*A34*A52*A61+A14*A23*A35*A52*A61-A13*A24* &
      A35*A52*A61-A15*A24*A32*A53*A61+A14*A25*A32*A53*A61+A15*A22*A34*A53* &
      A61-A12*A25*A34*A53*A61-A14*A22*A35*A53*A61+A12*A24*A35*A53*A61+A15*A23* &
      A32*A54*A61-A13*A25*A32*A54*A61-A15*A22*A33*A54*A61+A12*A25*A33*A54* &
      A61+A13*A22*A35*A54*A61-A12*A23*A35*A54*A61-A14*A23*A32*A55*A61+A13*A24* &
      A32*A55*A61+A14*A22*A33*A55*A61-A12*A24*A33*A55*A61-A13*A22*A34*A55* &
      A61+A12*A23*A34*A55*A61-A15*A24*A33*A51*A62+A14*A25*A33*A51*A62+A15*A23* &
      A34*A51*A62-A13*A25*A34*A51*A62-A14*A23*A35*A51*A62+A13*A24*A35*A51* &
      A62+A15*A24*A31*A53*A62-A14*A25*A31*A53*A62-A15*A21*A34*A53*A62+A11*A25* &
      A34*A53*A62+A14*A21*A35*A53*A62-A11*A24*A35*A53*A62-A15*A23*A31*A54* &
      A62+A13*A25*A31*A54*A62+A15*A21*A33*A54*A62-A11*A25*A33*A54*A62-A13*A21* &
      A35*A54*A62+A11*A23*A35*A54*A62+A14*A23*A31*A55*A62-A13*A24*A31*A55* &
      A62-A14*A21*A33*A55*A62+A11*A24*A33*A55*A62+A13*A21*A34*A55*A62-A11*A23* &
      A34*A55*A62+A15*A24*A32*A51*A63-A14*A25*A32*A51*A63-A15*A22*A34*A51* &
      A63+A12*A25*A34*A51*A63+A14*A22*A35*A51*A63-A12*A24*A35*A51*A63-A15*A24* &
      A31*A52*A63+A14*A25*A31*A52*A63+A15*A21*A34*A52*A63-A11*A25*A34*A52* &
      A63-A14*A21*A35*A52*A63+A11*A24*A35*A52*A63+A15*A22*A31*A54*A63-A12*A25* &
      A31*A54*A63-A15*A21*A32*A54*A63+A11*A25*A32*A54*A63+A12*A21*A35*A54* &
      A63-A11*A22*A35*A54*A63-A14*A22*A31*A55*A63+A12*A24*A31*A55*A63+A14*A21* &
      A32*A55*A63-A11*A24*A32*A55*A63-A12*A21*A34*A55*A63+A11*A22*A34*A55* &
      A63-A15*A23*A32*A51*A64+A13*A25*A32*A51*A64+A15*A22*A33*A51*A64-A12*A25* &
      A33*A51*A64-A13*A22*A35*A51*A64+A12*A23*A35*A51*A64+A15*A23*A31*A52* &
      A64-A13*A25*A31*A52*A64-A15*A21*A33*A52*A64+A11*A25*A33*A52*A64+A13*A21* &
      A35*A52*A64-A11*A23*A35*A52*A64-A15*A22*A31*A53*A64+A12*A25*A31*A53* &
      A64+A15*A21*A32*A53*A64-A11*A25*A32*A53*A64-A12*A21*A35*A53*A64+A11*A22* &
      A35*A53*A64+A13*A22*A31*A55*A64-A12*A23*A31*A55*A64-A13*A21*A32*A55* &
      A64+A11*A23*A32*A55*A64+A12*A21*A33*A55*A64-A11*A22*A33*A55*A64+A14*A23* &
      A32*A51*A65-A13*A24*A32*A51*A65-A14*A22*A33*A51*A65+A12*A24*A33*A51* &
      A65+A13*A22*A34*A51*A65-A12*A23*A34*A51*A65-A14*A23*A31*A52*A65+A13*A24* &
      A31*A52*A65+A14*A21*A33*A52*A65-A11*A24*A33*A52*A65-A13*A21*A34*A52* &
      A65+A11*A23*A34*A52*A65+A14*A22*A31*A53*A65-A12*A24*A31*A53*A65-A14*A21* &
      A32*A53*A65+A11*A24*A32*A53*A65+A12*A21*A34*A53*A65-A11*A22*A34*A53* &
      A65-A13*A22*A31*A54*A65+A12*A23*A31*A54*A65+A13*A21*A32*A54*A65-A11*A23* &
      A32*A54*A65-A12*A21*A33*A54*A65+A11*A22*A33*A54*A65

      COFACTOR(5,6) = -A15*A24*A33*A42* &
      A61+A14*A25*A33*A42*A61+A15*A23*A34*A42*A61-A13*A25*A34*A42*A61-A14*A23* &
      A35*A42*A61+A13*A24*A35*A42*A61+A15*A24*A32*A43*A61-A14*A25*A32*A43* &
      A61-A15*A22*A34*A43*A61+A12*A25*A34*A43*A61+A14*A22*A35*A43*A61-A12*A24* &
      A35*A43*A61-A15*A23*A32*A44*A61+A13*A25*A32*A44*A61+A15*A22*A33*A44* &
      A61-A12*A25*A33*A44*A61-A13*A22*A35*A44*A61+A12*A23*A35*A44*A61+A14*A23* &
      A32*A45*A61-A13*A24*A32*A45*A61-A14*A22*A33*A45*A61+A12*A24*A33*A45* &
      A61+A13*A22*A34*A45*A61-A12*A23*A34*A45*A61+A15*A24*A33*A41*A62-A14*A25* &
      A33*A41*A62-A15*A23*A34*A41*A62+A13*A25*A34*A41*A62+A14*A23*A35*A41* &
      A62-A13*A24*A35*A41*A62-A15*A24*A31*A43*A62+A14*A25*A31*A43*A62+A15*A21* &
      A34*A43*A62-A11*A25*A34*A43*A62-A14*A21*A35*A43*A62+A11*A24*A35*A43* &
      A62+A15*A23*A31*A44*A62-A13*A25*A31*A44*A62-A15*A21*A33*A44*A62+A11*A25* &
      A33*A44*A62+A13*A21*A35*A44*A62-A11*A23*A35*A44*A62-A14*A23*A31*A45* &
      A62+A13*A24*A31*A45*A62+A14*A21*A33*A45*A62-A11*A24*A33*A45*A62-A13*A21* &
      A34*A45*A62+A11*A23*A34*A45*A62-A15*A24*A32*A41*A63+A14*A25*A32*A41* &
      A63+A15*A22*A34*A41*A63-A12*A25*A34*A41*A63-A14*A22*A35*A41*A63+A12*A24* &
      A35*A41*A63+A15*A24*A31*A42*A63-A14*A25*A31*A42*A63-A15*A21*A34*A42* &
      A63+A11*A25*A34*A42*A63+A14*A21*A35*A42*A63-A11*A24*A35*A42*A63-A15*A22* &
      A31*A44*A63+A12*A25*A31*A44*A63+A15*A21*A32*A44*A63-A11*A25*A32*A44* &
      A63-A12*A21*A35*A44*A63+A11*A22*A35*A44*A63+A14*A22*A31*A45*A63-A12*A24* &
      A31*A45*A63-A14*A21*A32*A45*A63+A11*A24*A32*A45*A63+A12*A21*A34*A45* &
      A63-A11*A22*A34*A45*A63+A15*A23*A32*A41*A64-A13*A25*A32*A41*A64-A15*A22* &
      A33*A41*A64+A12*A25*A33*A41*A64+A13*A22*A35*A41*A64-A12*A23*A35*A41* &
      A64-A15*A23*A31*A42*A64+A13*A25*A31*A42*A64+A15*A21*A33*A42*A64-A11*A25* &
      A33*A42*A64-A13*A21*A35*A42*A64+A11*A23*A35*A42*A64+A15*A22*A31*A43* &
      A64-A12*A25*A31*A43*A64-A15*A21*A32*A43*A64+A11*A25*A32*A43*A64+A12*A21* &
      A35*A43*A64-A11*A22*A35*A43*A64-A13*A22*A31*A45*A64+A12*A23*A31*A45* &
      A64+A13*A21*A32*A45*A64-A11*A23*A32*A45*A64-A12*A21*A33*A45*A64+A11*A22* &
      A33*A45*A64-A14*A23*A32*A41*A65+A13*A24*A32*A41*A65+A14*A22*A33*A41* &
      A65-A12*A24*A33*A41*A65-A13*A22*A34*A41*A65+A12*A23*A34*A41*A65+A14*A23* &
      A31*A42*A65-A13*A24*A31*A42*A65-A14*A21*A33*A42*A65+A11*A24*A33*A42* &
      A65+A13*A21*A34*A42*A65-A11*A23*A34*A42*A65-A14*A22*A31*A43*A65+A12*A24* &
      A31*A43*A65+A14*A21*A32*A43*A65-A11*A24*A32*A43*A65-A12*A21*A34*A43* &
      A65+A11*A22*A34*A43*A65+A13*A22*A31*A44*A65-A12*A23*A31*A44*A65-A13*A21* &
      A32*A44*A65+A11*A23*A32*A44*A65+A12*A21*A33*A44*A65-A11*A22*A33*A44* &
      A65

      COFACTOR(6,6) = A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+A13*A25*     &
      A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51-A15*A24*A32*A43* &
      A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51-A12*A25*A34*A43*A51-A14*A22* &
      A35*A43*A51+A12*A24*A35*A43*A51+A15*A23*A32*A44*A51-A13*A25*A32*A44* &
      A51-A15*A22*A33*A44*A51+A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23* &
      A35*A44*A51-A14*A23*A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45* &
      A51-A12*A24*A33*A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-A15*A24* &
      A33*A41*A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52-A13*A25*A34*A41* &
      A52-A14*A23*A35*A41*A52+A13*A24*A35*A41*A52+A15*A24*A31*A43*A52-A14*A25* &
      A31*A43*A52-A15*A21*A34*A43*A52+A11*A25*A34*A43*A52+A14*A21*A35*A43* &
      A52-A11*A24*A35*A43*A52-A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21* &
      A33*A44*A52-A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44* &
      A52+A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+A11*A24* &
      A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52+A15*A24*A32*A41* &
      A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53+A12*A25*A34*A41*A53+A14*A22* &
      A35*A41*A53-A12*A24*A35*A41*A53-A15*A24*A31*A42*A53+A14*A25*A31*A42* &
      A53+A15*A21*A34*A42*A53-A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24* &
      A35*A42*A53+A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44* &
      A53+A11*A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-A14*A22* &
      A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53-A11*A24*A32*A45* &
      A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53-A15*A23*A32*A41*A54+A13*A25* &
      A32*A41*A54+A15*A22*A33*A41*A54-A12*A25*A33*A41*A54-A13*A22*A35*A41* &
      A54+A12*A23*A35*A41*A54+A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21* &
      A33*A42*A54+A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42* &
      A54-A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-A11*A25* &
      A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*A43*A54+A13*A22*A31*A45* &
      A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*A54+A11*A23*A32*A45*A54+A12*A21* &
      A33*A45*A54-A11*A22*A33*A45*A54+A14*A23*A32*A41*A55-A13*A24*A32*A41* &
      A55-A14*A22*A33*A41*A55+A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23* &
      A34*A41*A55-A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42* &
      A55-A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+A14*A22* &
      A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+A11*A24*A32*A43* &
      A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-A13*A22*A31*A44*A55+A12*A23* &
      A31*A44*A55+A13*A21*A32*A44*A55-A11*A23*A32*A44*A55-A12*A21*A33*A44* &
      A55+A11*A22*A33*A44*A55 
      
      A_INV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      ierr = 0

      !--- convert back to real for output
      AA_INV = real(A_Inv)

end subroutine INVERT_6x6

!--------------------------------------------------------------------------
! subroutine INVERT_DIAGONAL(A,A_inv,ierr)
!
! Matrix Inversion for a nxn diagnonal matrix
!--------------------------------------------------------------------------
subroutine INVERT_DIAGONAL(A,A_inv,ierr)
  real, dimension(:,:), intent(in):: A
  real, dimension(:,:), intent(out):: A_inv
  integer, intent(out):: ierr
  integer:: i, j, ni, nj
  real:: zero

  zero = epsilon(zero)
  ierr = 0
  A_inv = 0.0

  ni = size(A,1)
  nj = size(A,2)
 
  if(ni /= nj) then 
    ierr = 1
    return
  endif

  do i = 1, ni
     j = i
     if (abs(A(i,j)) > zero) then
        A_inv(i,i) = 1.0 / A(i,i)
     else
       ierr = 1
       return
     endif
  enddo

end subroutine INVERT_DIAGONAL

!--------------------------------------------------------------------
! subroutine FIND_BOUNDS(lat,lon,wlon,elon,slat,nlat,dateline_flg)
!
! This subroutine reads in a satellite data set containing lat/lon
! values and finds the western/eastern-most longitude and
! southern/northern-most latitude values.
!--------------------------------------------------------------------

subroutine FIND_BOUNDS(lat,lon,wlon,elon,slat,nlat,dateline_flg)

!Inputs and outputs are as follows:

!Inputs:
! lat    - array of satellite latitudes
! lon    - array of satellite longitudes

!Outputs
! wlon   - western-most longitude value in satellite data
! elon   - eastern-most longitude value in satellite data
! slat   - southern-most latitude value in satellite data
! nlat   - northern-most latitude value in satellite data

!Declare calling parameters

  real(kind=real4),dimension(:,:),intent(in) :: lat,lon
  real(kind=real4),intent(out) :: wlon,elon,slat,nlat
  integer(kind=int1), intent(out) :: dateline_flg

  integer(kind=int4) :: astatus, nx, ny
  real(kind=real4), dimension(:,:), allocatable :: dum

  nx = size(lon,1)
  ny = size(lon,2)

  allocate(dum(nx,ny),stat=astatus)
  if (astatus /= 0) then
    print "(a,'Not enough memory to allocate dummy longitude array.')",EXE_PROMPT
    stop
  endif

  dum = lon

  nlat = missing_value_real4
  slat = missing_value_real4
  elon = missing_value_real4
  wlon = missing_value_real4

  nlat = maxval(lat, mask = lat >= -90.0 .and. lat <= 90.0)
  slat = minval(lat, mask = lat >= -90.0 .and. lat <= 90.0)

  elon = maxval(lon, mask = lon >= -180.0 .and. lon <= 180.0)
  wlon = minval(lon, mask = lon >= -180.0 .and. lon <= 180.0)

  !Check if intl. date line is crossed
  dateline_flg = 0
  if (elon > 160.0 .and. wlon < -160) then
    where(lon < 0.0 .and. lon > missing_value_real4)
      dum = lon + 360.0
    endwhere
    !write(*,*) "international date line is assumed to be crossed"
    elon = maxval(dum, mask = lon >= -180.0 .and. lon <= 180.0)
    wlon = minval(dum, mask = lon >= -180.0 .and. lon <= 180.0)
    dateline_flg = 1
  endif

  deallocate(dum, stat=astatus)
  if (astatus /= 0) then
    print "(a,'Error deallocating dummy longitude array.')",EXE_PROMPT
    stop
  endif

end subroutine FIND_BOUNDS

!------------------------------------------------------------------------
! subroutine PACK_BYTES_I1(input_bytes,bit_depth,output_bytes)
! subroutine PACK_BYTES_I2(input_bytes,bit_depth,output_bytes)
!
! Routines to pack individual bytes into a single byte
!
! input:
! input_bytes - vector of bytes to be packed into output_byte
! bit_start - vector of bit starting positions (1-7) for each input byte
! bit_depth - vector of bit depths (1-7) for each input byte (total can not exceed 8)
!
! output: 
! output_byte - byte variable that holds the bit values of the input_bytes - can be i1 or i2
!
! local
!  n_in - number of elements of input vectors
!  i_in - index of input_bytes (1-n_in)
!  n_out - number of elements of output vectors
!  i_out - index of output_bytes (1-n_out)
!
! Note:
! 1.  if the input byte has information in bits greater then bit depth - they are removed 
!
!
! Example, pack an input_byte wth  bit_start = 2 and bit depth 3 
!
! input byte
!           x x x
! _ _ _ _ _ _ _ _
!
! result of first ishft
! x x x
! _ _ _ _ _ _ _ _
!
! result of second ishft
!       x x x
! _ _ _ _ _ _ _ _
!
! Author: Andrew Heidinger
!
! Version History:  
! February 2006 - Created
!-----------------------------------------------------------------------------------

!--- This Version packs into one byte words
   subroutine PACK_BYTES_I1(input_bytes,bit_depth,output_bytes)
    integer(kind=int1), dimension(:), intent(in):: input_bytes
    integer(kind=int4), dimension(:), intent(in):: bit_depth
    integer(kind=int1), dimension(:), intent(out):: output_bytes
    integer(kind=int1):: bit_start, bit_end, bit_offset
    integer(kind=int1):: temp_byte
    integer:: n_in,i_in,n_out,i_out
    integer, parameter:: word_bit_depth = 8
  
!--- determine size of vectors
    n_in = size(input_bytes)
    n_out = size(output_bytes)

!--- reset output byte
   output_bytes = 0

!--- initialize
   bit_offset = 0
   bit_start = 0
   bit_end = 0
   i_out = 1

!--- loop through input bytes
   do i_in = 1, n_in

!--- determine starting and ending bit locations
     bit_start = bit_offset + 1
     bit_end = bit_start + bit_depth(i_in) - 1

!--- determine if this input byte will fit on current output byte, if not go to next
     if (bit_end > word_bit_depth) then
      i_out = i_out + 1
      bit_offset = 0
      bit_start = bit_offset + 1
      bit_end = bit_start + bit_depth(i_in) - 1
     endif

!--- check for exceeding the space allowed for the packed bytes
     if (i_out > n_out) then
       print *, "ERROR: Insufficient space for bit packing" 
       return
     endif

!--- place input byte into correct position
     temp_byte =0
     temp_byte = ishft(input_bytes(i_in),word_bit_depth-bit_depth(i_in))   !first ishft
     temp_byte = ishft(temp_byte,bit_end - word_bit_depth)                 !second ishft

!--- modify output byte
     output_bytes(i_out) = output_bytes(i_out) + temp_byte

!--- update bit offset
     bit_offset = bit_offset + bit_depth(i_in)

   enddo

  end subroutine  PACK_BYTES_I1

!--- This Version packs into two byte words
   subroutine PACK_BYTES_I2(input_bytes,bit_depth,output_bytes)
    integer(kind=int1), dimension(:), intent(in):: input_bytes
    integer(kind=int4), dimension(:), intent(in):: bit_depth
    integer(kind=int2), dimension(:), intent(out):: output_bytes
    integer(kind=int1):: bit_start, bit_end, bit_offset
    integer(kind=int2):: temp_byte                 
    integer:: n_in,i_in,n_out,i_out
    integer, parameter:: word_bit_depth = 16

!--- determine size of vectors
    n_in = size(input_bytes)
    n_out = size(output_bytes)

!--- reset output byte
   output_bytes = 0

!--- initialize
   bit_offset = 0
   bit_start = 0
   bit_end = 0
   i_out = 1

!--- loop through input bytes
   do i_in = 1, n_in

!--- determine starting and ending bit locations
     bit_start = bit_offset + 1
     bit_end = bit_start + bit_depth(i_in) - 1

!--- determine if this input byte will fit on current output byte, if not go to next
     if (bit_end > word_bit_depth) then
      i_out = i_out + 1
      bit_offset = 0
      bit_start = bit_offset + 1
      bit_end = bit_start + bit_depth(i_in) - 1
     endif

!--- check for exceeding the space allowed for the packed bytes
     if (i_out > n_out) then
       print *, "ERROR: Insufficient space for bit packing"
       return
     endif

!--- place input byte into correct position
     temp_byte =0
     temp_byte = ishft(INT(input_bytes(i_in),kind=INT2),word_bit_depth-bit_depth(i_in))   !first ishft
     temp_byte = ishft(temp_byte,bit_end - word_bit_depth)                 !second ishft

!--- modify output byte
     output_bytes(i_out) = output_bytes(i_out) + temp_byte

!--- update bit offset
     bit_offset = bit_offset + bit_depth(i_in)

   enddo

  end subroutine  PACK_BYTES_I2

!====================================================================
! Function Name: Pearson_Corr
!
! Function:
!    Compute the Pearson Correlation Coefficient for two mxn arrays
!
! Description: Pearson's product-moment coefficient
!   
! Calling Sequence: BT_WV_BT_Window_Corr(Elem_Idx,Line_Idx) = Pearson_Corr( &
!                       sat%bt10(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                       sat%bt14(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                      Array_Width, Array_Hgt)
!   
!
! Inputs:
!   Array 1
!   Array 2
!   Elem_size
!   Line_size
!
! Outputs: 
!   Pearson Correlation coefficent
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Standard definition for Pearson correlation
!
!====================================================================
function Pearson_Corr(Array_One,Array_Two,Array_Width,Array_Hght,Invalid_Data_Mask)  &
         RESULT(Pearson_Corr_Coeff)
   real(kind=real4), intent(in), dimension(:,:):: Array_One
   real(kind=real4), intent(in), dimension(:,:):: Array_Two
   integer(kind=INT4), intent(in):: Array_Width
   integer(kind=INT4), intent(in):: Array_Hght
   integer(kind=INT1), intent(in), dimension(:,:):: Invalid_Data_Mask
   real(kind=real4), dimension(Array_Width,Array_Hght):: Pearson_Corr_Term_1
   real(kind=real4), dimension(Array_Width,Array_Hght):: Pearson_Corr_Term_2
   real(kind=real8):: Pearson_Corr_Top_Term_1
   real(kind=real8):: Pearson_Corr_Top_Term_2
   real(kind=real8):: Pearson_Corr_Bottom_Term_1
   real(kind=real8):: Pearson_Corr_Bottom_Term_2
   real(kind=real4):: Pearson_Corr_Coeff
   real(kind=real8):: Mean_Array_One
   real(kind=real8):: Mean_Array_Two
   real(kind=real8):: Sum_Array_One
   real(kind=real8):: Sum_Array_Two

   !--- skip computation for pixel arrays with any missing data
   if (sum(Invalid_Data_Mask) > 0) then
      Pearson_Corr_Coeff = Missing_Value_Real4
      return
   endif

   Sum_Array_One = sum(Array_One)
   Sum_Array_Two = sum(Array_Two)

   Mean_Array_One = Sum_Array_One / (Array_Width*Array_Hght)
   Mean_Array_Two = Sum_Array_Two / (Array_Width*Array_Hght)

   Pearson_Corr_Term_1 = Array_One - Mean_Array_One
   Pearson_Corr_Term_2 = Array_Two - Mean_Array_Two

   Sum_Array_One = sum(Pearson_Corr_Term_1)
   Sum_Array_Two = sum(Pearson_Corr_Term_2)

   Mean_Array_One = 0.0
   Mean_Array_Two = 0.0

   Pearson_Corr_Top_Term_1 = sum(Pearson_Corr_Term_1*Pearson_Corr_Term_2)
   
   Pearson_Corr_Top_Term_2 = (Sum_Array_One*Sum_Array_Two) / (Array_Width*Array_Hght)
   
   Pearson_Corr_Bottom_Term_1 = sum(Pearson_Corr_Term_1**2) - &
                                ((Sum_Array_One)**2) / (Array_Width*Array_Hght)
                                 
   Pearson_Corr_Bottom_Term_2 = sum(Pearson_Corr_Term_2**2) - &
                                ((Sum_Array_Two)**2) / (Array_Width*Array_Hght)

   Pearson_Corr_Coeff = (Pearson_Corr_Top_Term_1 - Pearson_Corr_Top_Term_2) / &
                         sqrt(Pearson_Corr_Bottom_Term_1 * &
                              Pearson_Corr_Bottom_Term_2)
   
 end function Pearson_Corr

!====================================================================
! Function Name: Covariance
!
! Function:
!    Compute the Covariance for two mxn arrays
!
! Description: Covariance = E(XY) - E(X)*E(Y)
!   
! Calling Sequence: BT_WV_BT_Window_Covar(Elem_Idx,Line_Idx) = Covariance( &
!                       sat%bt10(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                       sat%bt14(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                      Array_Width, Array_Hgt)
!   
!
! Inputs:
!   Array 1 - the first array (X)
!   Array 2 - the second array (Y)
!   Elem_size
!   Line_size
!
! Outputs: 
!   Covariance of X and Y
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Standard definition for the Covariance Computation
!
!====================================================================
function Covariance(Array_One,Array_Two,Array_Width,Array_Hght,Invalid_Data_Mask) &
         RESULT(Covar_Array_One_Array_Two)
   real(kind=real4), intent(in), dimension(:,:):: Array_One
   real(kind=real4), intent(in), dimension(:,:):: Array_Two
   integer(kind=INT4), intent(in):: Array_Width
   integer(kind=INT4), intent(in):: Array_Hght
   integer(kind=INT1), intent(in), dimension(:,:):: Invalid_Data_Mask

   real(kind=real8):: Mean_Array_One
   real(kind=real8):: Mean_Array_Two
   real(kind=real8):: Mean_Array_One_x_Array_Two
   real(kind=real8):: Sum_Array_One
   real(kind=real8):: Sum_Array_Two
   real(kind=real8):: Sum_Array_One_x_Array_Two
   real(kind=real4):: Covar_Array_One_Array_Two

   !--- skip computation for pixel arrays with any missing data
   if (sum(Invalid_Data_Mask) > 0) then
      Covar_Array_One_Array_Two = Missing_Value_Real4
      return
   endif

   Sum_Array_One = sum(Array_One)
   Sum_Array_Two = sum(Array_Two)

   Mean_Array_One = Sum_Array_One / (Array_Width*Array_Hght)
   Mean_Array_Two = Sum_Array_Two / (Array_Width*Array_Hght)

   Sum_Array_One_x_Array_Two = sum(Array_One*Array_Two)
   Mean_Array_One_x_Array_Two = Sum_Array_One_x_Array_Two / (Array_Width*Array_Hght)
   
   Covar_Array_One_Array_Two  = Mean_Array_One_x_Array_Two - &
                                Mean_Array_One * Mean_Array_Two 
   
 end function Covariance
 
!***************************************************************
!** compute general inverse of a matrix using svdcmp from 
!** Numerical Recipes
!** input:
!**        A - a matrix whose general inverse is sought
!**        m,p dimensions of matrix A   
!** output:
!**        A_inv - general inerse of A
!**        ierr - error code from svdcmp
!***************************************************************
subroutine GEN_INVERSE_SVD(A,A_inv,m,p,ierr)
 real (kind=ipre), intent(in), dimension(:,:):: A
 integer, intent(in):: m,p
 integer, intent(out):: ierr
 real (kind=ipre), intent(out), dimension(:,:):: A_inv
 integer:: i
 real (kind=ipre), dimension(p,p):: W,V
 real (kind=ipre), dimension(m,p):: U
 real (kind=ipre), dimension(p):: ww

 U = A
 if ( p == 0 ) then
   print *,"p = 0 in gen_inv",m,p,shape(W)
   stop
 endif
 call svdcmp(U,m,p,ww,V,ierr)
 W=0.0
  do i = 1,p
    W(i,i)=1.0/ww(i)
  enddo 
 A_inv = matmul(V,matmul(W,transpose(U)))

end subroutine GEN_INVERSE_SVD
!----------------------------------------------------------------

subroutine GEN_INVERSE_LU(A,A_inv,m,p,ierr)
 real (kind=ipre), intent(in), dimension(:,:):: A
 integer, intent(in):: m,p
 integer, intent(out):: ierr
 real (kind=ipre), intent(out), dimension(:,:):: A_inv
 
 integer:: i
 real:: diag_val

 if (m /= p) then 
    print *, 'LU decomposion only works for square matrix, stop'
    stop
 endif

 call lu_invert(real(A,real8),m,ierr,A_inv)

 if (.false.) then
    do i=1,m
       diag_val = A(i,i)*A_inv(i,i)
       if (abs(diag_val) < 0.8 .or. abs(diag_val) > 1.5) then
          ierr = 1
          print *, 'A*A_inverse is not equal to Identity Matrix ',diag_val
          exit
       endif
    enddo
 endif

end subroutine GEN_INVERSE_LU

 subroutine svdcmp(a,m,n,w,v,ierr)
   integer, intent(in):: m,n
   integer, intent(out):: ierr
!  real (kind=ipre), intent(inout), dimension(:,:):: a(mp,np),v(np,np),w(np)
   real (kind=ipre), intent(inout), dimension(:,:):: a
   real (kind=ipre), intent(out), dimension(:,:):: v
   real (kind=ipre), intent(out), dimension(:):: w
   integer, parameter:: NMAX=500
!CU    USES pythag
   integer:: i,its,j,jj,k,l,nm
   real (kind=ipre):: anorm,c,f,g,h,s,sscale,x,y,z
   real (kind=ipre), dimension(NMAX):: rv1
      ierr = 0
      g=0.0
      sscale=0.0
      anorm=0.0
      twentyfive: do i=1,n
        l=i+1
        rv1(i)=sscale*g
        g=0.0
        s=0.0
        sscale=0.0
        if(i<=m)then
          do  k=i,m
            sscale=sscale+abs(a(k,i))
          enddo
          if(sscale/=0.0)then
            do k=i,m
              a(k,i)=a(k,i)/sscale
              s=s+a(k,i)*a(k,i)
            end do
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do j=l,n
              s=0.0
              do k=i,m
                s=s+a(k,i)*a(k,j)
              end do
              f=s/h
              do k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
              end do
            end do
            do k=i,m
              a(k,i)=sscale*a(k,i)
            end do
          endif
        endif
        w(i)=sscale *g
        g=0.0
        s=0.0
        sscale=0.0
        if((i<=m).and.(i/=n))then
          do k=l,n
            sscale=sscale+abs(a(i,k))
          enddo
          if(sscale/=0.0)then
            do k=l,n
              a(i,k)=a(i,k)/sscale
              s=s+a(i,k)*a(i,k)
            enddo 
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do k=l,n
              rv1(k)=a(i,k)/h
            enddo 
            do j=l,m
              s=0.0
              do k=l,n
                s=s+a(j,k)*a(i,k)
              enddo
              do k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
              enddo
            enddo
            do k=l,n
              a(i,k)=sscale*a(i,k)
            enddo
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
       end do twentyfive
      do i=n,1,-1
        if(i < n)then
          if(g /= 0.0)then
            do j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
            enddo
            do j=l,n
              s=0.0
              do k=l,n
                s=s+a(i,k)*v(k,j)
              enddo
              do k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
              enddo
            enddo
          endif
          do j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
          enddo
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
      enddo
      do i=min(m,n),1,-1
        if (i == 0) then
          print *,"error in svdcmp",m,n,i
        endif
        l=i+1
        g=w(i)
        do j=l,n
          a(i,j)=0.0
        enddo
        if(g /= 0.0)then
          g=1.0/g
          do j=l,n
            s=0.0
            do k=l,m
              s=s+a(k,i)*a(k,j)
            enddo
            f=(s/a(i,i))*g
            do k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
            enddo
          enddo
          do j=i,m
            a(j,i)=a(j,i)*g
          enddo
        else
          do j= i,m
            a(j,i)=0.0
          end do
        endif
        a(i,i)=a(i,i)+1.0
      enddo
      fortynine: do k=n,1,-1
        fortyeight: do its=1,30
          do l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm) == anorm)  then
!              goto 2
               exit
            endif
            if((abs(w(nm))+anorm) == anorm) then
!                goto 1
                 exit
            endif
          enddo
  if((abs(rv1(l))+anorm) /= anorm)  then
!1      c=0.0
        c=0.0
        s=1.0
        do i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm) == anorm) then
!              goto 2
               exit
            endif
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
            enddo
         enddo
       endif
!2         z=w(k)
          z=w(k)
          if (l == k)then
            if(z < 0.0)then
              w(k)=-z
              do j=1,n
                v(j,k)=-v(j,k)
              enddo
            endif
!           goto 3
            exit
          endif
!         if(its.eq.30) pause 'no convergence in svdcmp'
          if(its == 30) then
             ierr = 1
             return
          endif
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0_ipre)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
            enddo
            z=pythag(f,h)
            w(j)=z
            if(z/=0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
            end do
          end do
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
        end do fortyeight
!3       continue
 end do fortynine
 end subroutine svdcmp
!***********************************************************************

 function pythag(a,b) result(pyth)
   real (kind=ipre), intent(in):: a,b
   real (kind=ipre):: pyth
   real (kind=ipre):: absa,absb
     absa=abs(a)
     absb=abs(b)
       if (absa > absb) then
           pyth = absa * sqrt(1.0+(absb/absa)**2)
       else
           if (absb == 0) then
               pyth = 0.0
           else
               pyth = absb*sqrt(1.0 + (absa/absb)**2)
           endif
       endif
 end function pythag

!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
Subroutine ludcmp(a,n,indx,d,code)

 integer n,nmax,i,imax,j,k
 real*16 tiny_num
 parameter(nmax=100,tiny_num=1.5D-16)
 real*8  amax,dum, sum_value, a(n,n),vv(nmax)
 integer code,d,indx(n)

 d=1; code=0

 do i=1,n
   amax=0.d0
   do j=1,n
     if (DABS(a(i,j)) .gt. amax) amax=DABS(a(i,j))
   end do ! j loop
   if (amax .lt. tiny_num) THEN
     code = 1
     return
   end if
   vv(i) = 1.d0 / amax
 end do ! i loop

 do j=1,n
   do i=1,j-1
     sum_value = a(i,j)
     do k=1,i-1
       sum_value = sum_value - a(i,k)*a(k,j)
     end do ! k loop
     a(i,j) = sum_value
   end do ! i loop
   amax = 0.d0
   do i=j,n
     sum_value = a(i,j)
     do k=1,j-1
       sum_value = sum_value - a(i,k)*a(k,j)
     end do ! k loop
     a(i,j) = sum_value
     dum = vv(i)*DABS(sum_value)
     if (dum .ge. amax) THEN
       imax = i
       amax = dum
     end if
   end do ! i loop

   if (j .ne. imax) then
     do k=1,n
       dum = a(imax,k)
       a(imax,k) = a(j,k)
       a(j,k) = dum
     end do ! k loop
     d = -d
     vv(imax) = vv(j)
   end if

   indx(j) = imax
   if (DABS(a(j,j)) < tiny_num) a(j,j) = tiny_num

   if (j .ne. n) then
     dum = 1.d0 / a(j,j)
     do i=j+1,n
       a(i,j) = a(i,j)*dum
     end do ! i loop
   end if
 end do ! j loop

 return
 
 end subroutine ludcmp

!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 subroutine lubksb(a,n,indx,b)
 real*8  sum_value, a(n,n),b(n)
 integer n,i,ii,j,ll,indx(n) 

 ii = 0

 do i=1,n
   ll = indx(i)
   sum_value = b(ll)
   b(ll) = b(i)
   if (ii .ne. 0) THEN
     do j=ii,i-1
       sum_value = sum_value - a(i,j)*b(j)
     end do ! j loop
   else if(sum_value.NE.0.d0) THEN
     ii = i
   end if
   b(i) = sum_value
 end do ! i loop

 do i=n,1,-1
   sum_value = b(i)
   if (i < n) then
     do j=i+1,n
       sum_value = sum_value - a(i,j)*b(j)
     end do ! j loop
   end if
   b(i) = sum_value / a(i,i)
 end do ! i loop

 return
 end subroutine lubksb


! call ludcmp and lubksb to invert square matrix using LU decomposition
  subroutine lu_invert(aa,nn,singular_flag,b_out)

  real(kind=real8),intent(in),dimension(:,:):: aa
  integer,intent(out):: singular_flag
  real(kind=ipre),intent(out),dimension(:,:) :: b_out

  real(kind=real8),dimension(:,:),allocatable :: b
  integer,dimension(:),allocatable :: indx
  integer i,d,nn

  allocate(b(nn,nn))
  allocate(indx(nn))

  call ludcmp(aa,nn,indx,d,singular_flag)

  b = 0
  do i=1,nn
     b(i,i)=1
  enddo

  if (singular_flag == 0) then
     do i=1,nn
        call lubksb(aa,nn,indx,b(i,:))
     enddo
  endif

  b_out = real(b,real4)

  deallocate(b)
  deallocate(indx)

  end subroutine lu_invert

!------------------------------------------------------------------------------------- 
end module NUMERICAL_ROUTINES_MOD
