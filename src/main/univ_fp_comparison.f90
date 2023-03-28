!------------------------------------------------------------------------------
! CLAVR-x (CLouds from AVHRR - eXtended) PROCESSING SOFTWARE
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
!------------------------------------------------------------------------------

!```````````````````````````````````````````````````````````````````
module univ_fp_comparison_mod

!+ Defines operators that rigorously compare floating-point values.
!
!  Floating-point (in)equality operators are:
!          (not case-sensitive)
!    .EQfp.  :  replaces '.eq.' or '=='
!    .NEfp.  :  replaces '.ne.' or '/='
!    .GEfp.  :  replaces '.ge.' or '>='
!    .LEfp.  :  replaces '.le.' or '<='
!
!  Usage Example:
!
!  use univ_fp_comparison_mod, only: operator(.EQfp.), operator(.GEfp.)
!
!  if (c .EQfp. d) then
!     if (f .GEfp. e) then
  
implicit none

private
public :: operator(.EQfp.)
public :: operator(.NEfp.)
public :: operator(.GEfp.)
public :: operator(.LEfp.)

interface operator (.EQfp.)
   module procedure eq_fp_f4_0d, eq_fp_f4_1d, eq_fp_f4_2d, eq_fp_f4_3d,  &
        eq_fp_f4_4d, eq_fp_f4_5d, eq_fp_f4_6d, eq_fp_f4_7d
   module procedure eq_fp_f4_1d_0d, eq_fp_f4_2d_0d, eq_fp_f4_3d_0d,  &
        eq_fp_f4_4d_0d, eq_fp_f4_5d_0d, eq_fp_f4_6d_0d, eq_fp_f4_7d_0d
   module procedure eq_fp_f8_0d, eq_fp_f8_1d, eq_fp_f8_2d, eq_fp_f8_3d,  &
        eq_fp_f8_4d, eq_fp_f8_5d, eq_fp_f8_6d, eq_fp_f8_7d
   module procedure eq_fp_f8_1d_0d, eq_fp_f8_2d_0d, eq_fp_f8_3d_0d,  &
        eq_fp_f8_4d_0d, eq_fp_f8_5d_0d, eq_fp_f8_6d_0d, eq_fp_f8_7d_0d
!%   module procedure eq_fp_f16_0d, eq_fp_f16_1d, eq_fp_f16_2d, eq_fp_f16_3d,  &
!%        eq_fp_f16_4d, eq_fp_f16_5d, eq_fp_f16_6d, eq_fp_f16_7d
!%   module procedure eq_fp_f16_1d_0d, eq_fp_f16_2d_0d, eq_fp_f16_3d_0d,  &
!%        eq_fp_f16_4d_0d, eq_fp_f16_5d_0d, eq_fp_f16_6d_0d, eq_fp_f16_7d_0d
end interface

interface operator (.NEfp.)
   module procedure ne_fp_f4_0d, ne_fp_f4_1d, ne_fp_f4_2d, ne_fp_f4_3d,  &
        ne_fp_f4_4d, ne_fp_f4_5d, ne_fp_f4_6d, ne_fp_f4_7d
   module procedure ne_fp_f4_1d_0d, ne_fp_f4_2d_0d, ne_fp_f4_3d_0d,  &
        ne_fp_f4_4d_0d, ne_fp_f4_5d_0d, ne_fp_f4_6d_0d, ne_fp_f4_7d_0d
   module procedure ne_fp_f8_0d, ne_fp_f8_1d, ne_fp_f8_2d, ne_fp_f8_3d,  &
        ne_fp_f8_4d, ne_fp_f8_5d, ne_fp_f8_6d, ne_fp_f8_7d
   module procedure ne_fp_f8_1d_0d, ne_fp_f8_2d_0d, ne_fp_f8_3d_0d,  &
        ne_fp_f8_4d_0d, ne_fp_f8_5d_0d, ne_fp_f8_6d_0d, ne_fp_f8_7d_0d
!%   module procedure ne_fp_f16_0d, ne_fp_f16_1d, ne_fp_f16_2d, ne_fp_f16_3d,  &
!%        ne_fp_f16_4d, ne_fp_f16_5d, ne_fp_f16_6d, ne_fp_f16_7d
!%   module procedure ne_fp_f16_1d_0d, ne_fp_f16_2d_0d, ne_fp_f16_3d_0d,  &
!%        ne_fp_f16_4d_0d, ne_fp_f16_5d_0d, ne_fp_f16_6d_0d, ne_fp_f16_7d_0d
end interface

interface operator (.GEfp.)
   module procedure ge_fp_f4_0d, ge_fp_f4_1d, ge_fp_f4_2d, ge_fp_f4_3d,  &
        ge_fp_f4_4d, ge_fp_f4_5d, ge_fp_f4_6d, ge_fp_f4_7d
   module procedure ge_fp_f4_1d_0d, ge_fp_f4_2d_0d, ge_fp_f4_3d_0d,  &
        ge_fp_f4_4d_0d, ge_fp_f4_5d_0d, ge_fp_f4_6d_0d, ge_fp_f4_7d_0d
   module procedure ge_fp_f8_0d, ge_fp_f8_1d, ge_fp_f8_2d, ge_fp_f8_3d,  &
        ge_fp_f8_4d, ge_fp_f8_5d, ge_fp_f8_6d, ge_fp_f8_7d
   module procedure ge_fp_f8_1d_0d, ge_fp_f8_2d_0d, ge_fp_f8_3d_0d,  &
        ge_fp_f8_4d_0d, ge_fp_f8_5d_0d, ge_fp_f8_6d_0d, ge_fp_f8_7d_0d
!%   module procedure ge_fp_f16_0d, ge_fp_f16_1d, ge_fp_f16_2d, ge_fp_f16_3d,  &
!%        ge_fp_f16_4d, ge_fp_f16_5d, ge_fp_f16_6d, ge_fp_f16_7d
!%   module procedure ge_fp_f16_1d_0d, ge_fp_f16_2d_0d, ge_fp_f16_3d_0d,  &
!%        ge_fp_f16_4d_0d, ge_fp_f16_5d_0d, ge_fp_f16_6d_0d, ge_fp_f16_7d_0d
end interface

interface operator (.LEfp.)
   module procedure le_fp_f4_0d, le_fp_f4_1d, le_fp_f4_2d, le_fp_f4_3d,  &
        le_fp_f4_4d, le_fp_f4_5d, le_fp_f4_6d, le_fp_f4_7d
   module procedure le_fp_f4_1d_0d, le_fp_f4_2d_0d, le_fp_f4_3d_0d,  &
        le_fp_f4_4d_0d, le_fp_f4_5d_0d, le_fp_f4_6d_0d, le_fp_f4_7d_0d
   module procedure le_fp_f8_0d, le_fp_f8_1d, le_fp_f8_2d, le_fp_f8_3d,  &
        le_fp_f8_4d, le_fp_f8_5d, le_fp_f8_6d, le_fp_f8_7d
   module procedure le_fp_f8_1d_0d, le_fp_f8_2d_0d, le_fp_f8_3d_0d,  &
        le_fp_f8_4d_0d, le_fp_f8_5d_0d, le_fp_f8_6d_0d, le_fp_f8_7d_0d
!%   module procedure le_fp_f16_0d, le_fp_f16_1d, le_fp_f16_2d, le_fp_f16_3d,  &
!%        le_fp_f16_4d, le_fp_f16_5d, le_fp_f16_6d, le_fp_f16_7d
!%   module procedure le_fp_f16_1d_0d, le_fp_f16_2d_0d, le_fp_f16_3d_0d,  &
!%        le_fp_f16_4d_0d, le_fp_f16_5d_0d, le_fp_f16_6d_0d, le_fp_f16_7d_0d
end interface


contains


include "univ_fp_comparison.f90.inc01"

include "univ_fp_comparison.f90.inc02"

end module univ_fp_comparison_mod
