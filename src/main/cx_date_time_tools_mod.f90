! $Id: cx_date_time_tools_mod.f90 3735 2020-02-20 22:34:16Z heidinger $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE
!
! NAME: cx_date_time_tools_mod.f90 (src)
!       CX_DATE_TIME_TOOLS_MOD (program)
!
! PURPOSE: library of useful scientifix routines and  functions
!
! Description: 
!
! AUTHORS:
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
!   JULIAN
!   COMPUTE_MONTH
!   COMPUTE_DAY
!   COMPUTE_TIME_HOURS
!
!  HISTORY
!   Removed class tools 26 January 2019 (AW)
!
!--------------------------------------------------------------------------------------
 module CX_DATE_TIME_TOOLS_MOD
  
  implicit none
  private
  integer, parameter :: r15 = selected_real_kind(7)
 
  public :: JULIAN
  public :: COMPUTE_TIME_HOURS
  public :: COMPUTE_MONTH
  public :: COMPUTE_DAY
  public :: LEAP_YEAR_FCT
 
 contains 
 !-------------------------------------------------
 ! subroutine JULIAN(iday,imonth,iyear,jday)
 ! Computes julian day (1-365/366)
 ! input:
 !         iday - integer day
 !         imonth - integer month
 !         iyear - integer year (2 or four digits)
 !         
 ! output : jday - julian day or day of year
 !--------------------------------------------------
   subroutine JULIAN(iday,imonth,iyear,jday)

        integer, intent(in)::  iday,imonth,iyear
        integer, intent(out):: jday
        integer::  j
        integer, dimension(12)::  jmonth

        jmonth = reshape ((/31,28,31,30,31,30,31,31,30,31,30,31/),(/12/))

        jday = iday
        if (modulo(iyear,4) == 0) then
            jmonth(2)=29
        endif

        do j = 1,imonth-1
           jday = jday + jmonth(j)
        end do

  end subroutine JULIAN

  !--------------------------------------------
  ! compute the month of the year (1-12)
  !---------------------------------------------
  function COMPUTE_MONTH(jday,ileap) result(month)
   integer, intent(in):: ileap
   integer, intent(in):: jday
   integer:: month
   integer :: days_of_month (0:12) 
   integer :: i
   
   days_of_month = [0,31,28,31,30,31,30,31,31,30,31,30,31]
   if ( ileap == 1) days_of_month(2) = 29
   month = -1
   if ( jday < 1 .or. jday > sum(days_of_month) ) then
       return
   end if   
   
   do i = 1, 12        
      if (jday <= sum(days_of_month(0:i))  ) then
         month = i
         exit
      end if   
   end do

  end function COMPUTE_MONTH
  !--------------------------------------------
  ! compute the day of month
  !---------------------------------------------
  function COMPUTE_DAY(jday,ileap) result(day)
   integer, intent(in):: ileap
   integer, intent(in):: jday
   integer:: day
   integer :: days_of_month (0:12) 
   integer :: i
   
   days_of_month = [0,31,28,31,30,31,30,31,31,30,31,30,31]
   if ( ileap == 1) days_of_month(2) = 29
   
   day = -1
   if ( jday < 1 .or. jday > sum(days_of_month) ) then      
      return
   end if   
   
   do i = 1, 12 
      if (jday <= sum(days_of_month(0:i))  ) then
         day = jday - sum(days_of_month(0:(i-1)))
         exit
      end if   
   end do

  end function COMPUTE_DAY

  !---------------------------------------------------------------------
  ! function leap_year_fct(yyyy) result(leap_flg)
  !
  ! Function to determine if an input year is a leap year.
  !---------------------------------------------------------------------
  function leap_year_fct(yyyy) result(leap_flg)
   integer, intent(in) :: yyyy
   integer :: leap_flg

   leap_flg = 0
 
   if ((modulo(yyyy,4) == 0 .and. modulo(yyyy,100) /= 0) .or. &
        modulo(yyyy,400) == 0) leap_flg = 1

   return

  end function leap_year_fct

  !---------------------------------------------------------------------
  ! function to return the system in fractional hours since midnight
  !---------------------------------------------------------------------
  function COMPUTE_TIME_HOURS()   result(time_hours)
   character(len=8):: system_date
   character(len=10):: system_time
   character(len=5):: system_time_zone
   integer, dimension(8):: system_time_value

   real:: time_hours

   call DATE_AND_TIME(system_date,system_time,system_time_zone, system_time_value)

   time_hours = real(system_time_value(5)) +  &
                     (real(system_time_value(6)) + &
                      real(system_time_value(7) + &
                      real(system_time_value(8))/1000.0)/60.0)/60.0
   return

  end function COMPUTE_TIME_HOURS

!------------------------------------------------------------------------------------- 
! end of module
!------------------------------------------------------------------------------------- 

end module CX_DATE_TIME_TOOLS_MOD

