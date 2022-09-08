!```````````````````````````````````````````````````````````````````
program determine_kinds_F90

!+ F90+ utility program that determines the variable type "kind" options for
!   the current CPU/OS combination; writes F90+ module code that contains its
!   findings

implicit none

!== Local declarations:

integer, parameter :: IOu = 34

integer :: i1_B, i2_B, i4_B, i8_B, i16_B, i32_B, in_B,  &
           f4_B, f8_B, f16_B, f32_B, fn_B, cn_B, bn_B,  &
           itmp, fx_B, fx, fxh_B, fxh
character(len=128) :: app_label

!--- "Native" F90+ kinds on the current machine/system:
integer, parameter :: in = kind(0),  &        ! Integer
                      fn = kind(0.),  &       ! Floating point
                      bn = kind(.FALSE.),  &  ! Logical
                      cn = kind('a')         ! Character

!--- Integer type kinds; declare archetype integer variables:

integer, parameter :: i1 = selected_int_kind(2),  &    ! can hold +/- 10^2
                      i2 = selected_int_kind(4),  &    ! can hold +/- 10^4
                      i4 = selected_int_kind(9),  &    ! can hold +/- 10^9
                      i8 = selected_int_kind(18),  &   ! can hold +/- 10^18
                      i16 = selected_int_kind(37),  &  ! can hold +/- 10^37
                      i32 = selected_int_kind(76)      ! can hold +/- 10^76

integer(in) :: in_v
integer(i1) :: i1_v
integer(i2) :: i2_v
integer(i4) :: i4_v
integer(i8) :: i8_v
integer(max(in, i16)) :: i16_v
integer(max(in, i32)) :: i32_v

!--- Real type kinds; declare archetype real variables:

! Capable of: 6 decimal digit precision, decimal exponent of +/- 30
integer, parameter :: f4 = selected_real_kind(6, 30)

! Capable of: 12 decimal digit precision, decimal exponent of +/- 200
integer, parameter :: f8 = selected_real_kind(12, 200)

! Capable of: 18 decimal digit precision, decimal exponent of +/- 1300
integer, parameter :: f16 = selected_real_kind(18, 1300)

! Capable of: 24 decimal digit precision, decimal exponent of +/- 9000
integer, parameter :: f32 = selected_real_kind(24, 9000)

real(fn) :: fn_v
real(f4) :: f4_v
real(f8) :: f8_v
real(max(fn, f16)) :: f16_v
real(max(fn, f32)) :: f32_v

logical :: i16_valid, i32_valid, f16_valid, f32_valid
 
!=== Executable statements:

i16_valid = ( i16 > 0 )
i32_valid = ( i32 > 0 )
f16_valid = ( f16 > 0 )
f32_valid = ( f32 > 0 ) 

! Determine corresponding byte-sizes corresponding to each kind:

in_B = bit_size(in_v)/8
i1_B = bit_size(i1_v)/8; call check_bytesize(i1_B, 1, 'i1')
i2_B = bit_size(i2_v)/8; call check_bytesize(i2_B, 2, 'i2')
i4_B = bit_size(i4_v)/8; call check_bytesize(i4_B, 4, 'i4')
i8_B = bit_size(i8_v)/8; call check_bytesize(i8_B, 8, 'i8')

if (i16_valid) then
   i16_B = bit_size(i16_v)/8; call check_bytesize(i16_B, 16, 'i16')
end if
if (i32_valid) then
   i32_B = bit_size(i32_v)/8; call check_bytesize(i32_B, 32, 'i32')
end if

itmp = maxexponent(fn_v)
fn_B = guess_at_fptype_B(itmp)
itmp = maxexponent(f4_v)
f4_B = guess_at_fptype_B(itmp); call check_bytesize(f4_B, 4, 'f4')
itmp = maxexponent(f8_v)
f8_B = guess_at_fptype_B(itmp); call check_bytesize(f8_B, 8, 'f8')

if (f16_valid) then
   itmp = maxexponent(f16_v)
   f16_B = guess_at_fptype_B(itmp); call check_bytesize(f16_B, 16, 'f16')
end if
if (f32_valid) then
   itmp = maxexponent(f32_v)
   f32_B = guess_at_fptype_B(itmp); call check_bytesize(f32_B, 32, 'f32')
end if

cn_B = 1    ! Non-standard, but widespread where "short" alphabets are in use
bn_B = in_B      ! Fortran standard requires this

  ! Set the "preferred" floating-point kind (can be used in situations where
  !  any available precision is OK):
fx = f8       ! Default to a higher-precision option
fx_B = f8_B   !
#ifdef USE_F4_IF_POSSIBLE
fx = f4           ! Down-select to a lower-precision option
fx_B = f4_B       !
#endif

  ! Set the higher-precision floating-point kind (can be used in situations
  !  where greater precision is always required):
fxh = f8
fxh_B = f8_B

! Generate primary source code precision parameter module (write to file):

open(IOu, file='kinds.A.f90', status='unknown', form='formatted')

write (IOu, '(a)') '!```````````````````````````````````````````````````````````````````'
write (IOu, '(a,/)') 'module univ_kind_defs_mod'

write (IOu, '(a)') '!+ F95+ module containing machine precision-related items (for present OS,'
write (IOu, '(a)') '!   system architecture, and/or application)'
write (IOu, '(a)') '!'
write (IOu, '(a)') '! FOR USE IN ALL PROCEDURES; EXAMPLES:'
write (IOu, '(a)') '!'
write (IOu, '(a)') '!   use univ_kind_defs_mod, only: f4, i2_B, in, f8'
write (IOu, '(a)') '!'
write (IOu, '(a)') '!   real(f4) :: green'
write (IOu, '(a)') '!   integer(in) :: blue'
write (IOu, '(a)') '!   ... , recl=i2_B*nzyx'
write (IOu, '(a)') '!   integer(in) :: yellow'
write (IOu, '(a)') '!   real(f8) :: red'
write (IOu, '(a)') '!   green = green+34.e4_f4'
write (IOu, '(a)') '!   red = red*34._f8'

write (IOu, '(a,/)') '!~ (re)generated'

write (IOu, '(a,/)') 'implicit none'

write (IOu, '(a,/)') '!**** "Native" kinds (LOGICAL AND CHARACTER) ****:'
write (IOu, '(a,i2,a,i2,a)') 'integer(', in,  &
     '), parameter :: bn = ', bn, ',  &'
write (IOu, '(26x,a,i2)') 'cn = ', cn

write (IOu, '(/,a,/)') '!**** Supported INTEGER kinds ****:'
write (IOu, '(a,i2,a,i2,a)') 'integer(', in,  &
     '), parameter :: in = ', in, ',  &'
write (IOu, '(26x,a,i2,a)') 'i1 = ', i1, ',  &'
write (IOu, '(26x,a,i2,a)') 'i2 = ', i2, ',  &'
write (IOu, '(26x,a,i2,a)') 'i4 = ', i4, ',  &'
if (i16_valid) then
   write (IOu, '(26x,a,i2,a)') 'i8 = ', i8, ',  &'
   if (i32_valid) then
      write (IOu, '(26x,a,i2,a)') 'i16 = ', i16, ',  &'
   else
      write (IOu, '(26x,a,i2,a)') 'i16 = ', i16, ' '
   end if
else
   write (IOu, '(26x,a,i2,a)') 'i8 = ', i8, ' '
end if
if (i32_valid) then
   write (IOu, '(26x,a,i2,a)') 'i32 = ', i32, ' '
end if

write (IOu, '(/,a,/)') '!**** Supported REAL kinds ****:'
write (IOu, '(a,i2,a,i2,a)') 'integer(', in,  &
     '), parameter :: fn = ', fn, ',  &'
write (IOu, '(26x,a,i2,a)') 'f4 = ', f4, ',  &'
if (f16_valid) then
   write (IOu, '(26x,a,i2,a)') 'f8 = ', f8, ',  &'
   if (f32_valid) then
      write (IOu, '(26x,a,i2,a)') 'f16 = ', f16, ',  &'
   else
      write (IOu, '(26x,a,i2,a)') 'f16 = ', f16, ' '
   end if
else
   write (IOu, '(26x,a,i2,a)') 'f8 = ', f8, ' '
end if
if (f32_valid) then
   write (IOu, '(26x,a,i2,a)') 'f32 = ', f32, ' '
end if

write (IOu, '(/,a,/)') '!**** Byte-sizes ****:'

write (IOu, '(a,i2,a,i2,a)') 'integer(', in,  &
     '), parameter :: in_B = ', in_B, ',  &'
write (IOu, '(26x,a,i2,a)') 'i1_B = ', i1_B, ',  &'
write (IOu, '(26x,a,i2,a)') 'i2_B = ', i2_B, ',  &'
write (IOu, '(26x,a,i2,a)') 'i4_B = ', i4_B, ',  &'
if (i16_valid) then
   write (IOu, '(26x,a,i2,a)') 'i8_B = ', i8_B, ',  &'
   if (i32_valid) then
      write (IOu, '(26x,a,i2,a)') 'i16_B = ', i16_B, ',  &'
   else
      write (IOu, '(26x,a,i2,a)') 'i16_B = ', i16_B, ' '
   end if
else
   write (IOu, '(26x,a,i2,a)') 'i8_B = ', i8_B, ' '
end if
if (i32_valid) then
   write (IOu, '(26x,a,i2,a)') 'i32_B = ', i32_B, ' '
end if

write (IOu, '(/,a,i2,a,i2,a)') 'integer(', in,  &
     '), parameter :: fn_B = ', fn_B, ',  &'
write (IOu, '(26x,a,i2,a)') 'f4_B = ', f4_B, ',  &'
if (f16_valid) then
   write (IOu, '(26x,a,i2,a)') 'f8_B = ', f8_B, ',  &'
   if (f32_valid) then
      write (IOu, '(26x,a,i2,a)') 'f16_B = ', f16_B, ',  &'
   else
      write (IOu, '(26x,a,i2,a)') 'f16_B = ', f16_B, ' '
   end if
else
   write (IOu, '(26x,a,i2,a)') 'f8_B = ', f8_B, ' '
end if
if (f32_valid) then
   write (IOu, '(26x,a,i2,a)') 'f32_B = ', f32_B, ' '
end if

write (IOu, '(/,a,i2,a,i2,a)') 'integer(', in,  &
     '), parameter :: bn_B = ', bn_B, ',  &'
write (IOu, '(26x,a,i2,a)') 'cn_B = ', cn_B, ' '

write (IOu, '(/,a,/)') 'end module univ_kind_defs_mod'

close(IOu)


! Read information from configuration file; generate secondary source code
!  precision parameter module (write to file):

open(IOu, file='build_env_configure-util_f90.cfg', status='old',  &
     form='formatted')
read (IOu, *) app_label
close(IOu)

open(IOu, file='kinds.B.f90', status='unknown', form='formatted')

write (IOu, '(a)') '!```````````````````````````````````````````````````````````````````'
write (IOu, '(a,/)') 'module '//trim(app_label)//'_kinds_mod'

write (IOu, '(a)') '!+ F95+ module containing "preferred" REAL kind'
write (IOu, '(a)') '!   definition(s) *for the specific software that this file'
write (IOu, '(a)') '!   is part of*.  Consistent use of these enables trivial '
write (IOu, '(a)') '!   switching of floating-point calculation precision.'
write (IOu, '(a)') '!'
write (IOu, '(a)') '! FOR USE IN ALL PROCEDURES; EXAMPLES:'
write (IOu, '(a)') '!'
write (IOu, '(a)') '!   use '//trim(app_label)//  &
     '_kinds_mod, only: fx, fx_B, fxh, fxh_B'
write (IOu, '(a)') '!'
write (IOu, '(a)') '!   real(fx) :: orange'
write (IOu, '(a)') '!   real(fxh) :: red, plum'
write (IOu, '(a)') '!   ... , recl=fxh_B*nzyx'
write (IOu, '(a)') '!   orange = orange*34._fx'
write (IOu, '(a,/)') '!   plum = red-865._fxh'

write (IOu, '(a,/)') '!~ (re)generated'

write (IOu, '(a,/)') 'implicit none'

write (IOu, '(/,a,/,a,/)')  &
     '!**** The "preferred" REAL kind(s) ****:'
write (IOu, '(a,i2,a,i2)') 'integer(', in,  &
     '), parameter :: fx = ', fx
write (IOu, '(a,i2,a,i2)') 'integer(', in,  &
     '), parameter :: fxh = ', fxh

write (IOu, '(/,a,/)') '!**** Byte-sizes ****:'

write (IOu, '(a,i2,a,i2)') 'integer(', in,  &
     '), parameter :: fx_B = ', fx_B
write (IOu, '(a,i2,a,i2)') 'integer(', in,  &
     '), parameter :: fxh_B = ', fxh_B

write (IOu, '(/,a,/)') 'end module '//trim(app_label)//'_kinds_mod'

close(IOu)


contains


!```````````````````````````````````````````````````````````````````
function guess_at_fptype_B(maxexp) result(fptype_B_guess)

!+ Utility function that attempts to determine the byte-size of a
!   floating point variable type

implicit none

!== Arguments:

integer, intent(IN) :: maxexp 

integer :: fptype_B_guess    ! Declaration of procedure result

!== Local declarations:

! NONE

!=== Executable statements:

if (maxexp >= 100 .and. maxexp <= 150) then     ! likely a 4-byte type
   fptype_B_guess = 4
else if(maxexp >= 1000 .and. maxexp <= 1050) then     ! likely a 8-byte type
   fptype_B_guess = 8
else if(maxexp >= 16000 .and. maxexp <= 16500) then     ! likely a 16-byte type
   fptype_B_guess = 16
else if(maxexp > 16500) then     ! likely a 32-byte type
   fptype_B_guess = 32
else
   fptype_B_guess = 0   ! Unknown
end if

end function guess_at_fptype_B


!```````````````````````````````````````````````````````````````````
subroutine check_bytesize(determined, anticipated, ID_str)

!+ Utility procedure that generates a non-fatal warning message if the
!   two integer argument values are not equal

implicit none

!== Arguments:

integer, intent(IN) :: determined,anticipated
character(len=*), intent(IN) :: ID_str

!== Local declarations:

! NONE

!=== Executable statements:

if (determined /= anticipated) then
   write (*, '(3a,/,9x,a,i2,a,i2,a)')  &
        'WARNING: Unanticipated byte size for ', trim(ID_str), ': ',  &
        'expected ', anticipated, ' bytes , but was determined to be ',  &
        determined, ' bytes'
end if

end subroutine check_bytesize


end program determine_kinds_F90
