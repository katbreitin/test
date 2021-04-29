! $Id:$
module class_time_date
! To-do revise the julain_cmp and cdata tool 
!-----------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: class_time_date
! PURPOSE:  Date and time methods
!
! DESCRIPTION: 
!
!    This is FORTRAN Object-oriented class to deal with date and time.
!    
!    a date_object is a FORTRAN derived type which holds all data and methods.
!    
!    Simple example session:
!       
!        USE class_time_date,only: date_type
!
!         type(date_type) :: t0
!    
!         call t0.set_date( year=2018,month=11,day=12,hour=14,minute = 15,second = 12)
!         print*,t0.date_string('yyyy_doy')
!
!         returns a string for the given format ( here year and day of year)
!
!     ALL Public methods:
!
!        SET_DATE                  : set date with year month, day, hour ,minute and second
!        SET_DATE_WITH_DOY         : set date with year, day of year (1-366) , day, hour ,minute and second
!        SET_JDAY                  : set date with fractional julian day
!        DAYS_SINCE_PROJ_TIME0
!        PRINT_DATA                ; writes information on screen. 
!        DATE_STRING               ; Returns a date and time string in many possible format options
!        ADD_TIME                  ; 
!        ADD_DAYS                  ;
!        GET_DATE                  .
!   
!
! AUTHORS:
!
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!
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
! This copyright pertains to all routines in the CLAVR-x system unless stated
!
! REVISION HISTORY:
!  January 2019
!
!
!
!-------------------------------------------------------------------------


implicit none
  private
  integer, parameter :: r15 = selected_real_kind(7)
  
  !   fix for allocatable variable -1/06/2017 AW

  type, public :: date_type 
      integer :: year = 1970
      integer :: month = 1
      integer :: day = 1
      integer :: hour = 0
      integer :: minute = 0 
      integer :: second = 0
      integer :: dayOfYear
      real :: hour_frac
      integer :: msec_of_day
      character(4) :: yyyy
      character(2) :: mm
      character(2) :: dd
      character(6) :: yymmdd
      character(8) :: yymmddhh
      character(10) :: yyyymmdd
      real (kind = r15) :: julday
      real :: mod_julday
           
 contains
      procedure :: set_date 
      procedure :: set_date_with_doy
      procedure :: set_jday
      procedure :: days_since_proj_time0
      procedure :: print_data 
      procedure :: date_string 
      procedure :: add_time
      procedure :: add_days
      procedure :: get_date
      procedure , private :: update 
      procedure , private :: set_julday
      
      procedure , private :: update_from_jd
      procedure ::  period_16 
      procedure :: next_6h 
      procedure :: next_3h  
      procedure, pass :: init
       
  
  end type date_type
  
  
  
  public :: time_mid
  public :: time_is_in_window
  public :: time_diff_weight
  
  contains
  
  subroutine init (self , ii)
    class(date_type), intent(out) :: self
    integer, optional, intent(in) :: ii
    integer :: i_init 
    
    i_init = 0 
    if (present(ii)) i_init = ii
    call self % set_date (1970,1,1,0,0,0)
    if (i_init .gt. 0) then
      call self % set_date (2010,1,1,0,0,0)
    end if
  
  end subroutine init
  
  ! -------------
  !    Populates object with time data
  !   
  ! ---------------------
  subroutine set_date (this , year , month, day, hour ,minute, second)
      class ( date_type) :: this
      integer , optional :: year , month ,day, hour , minute, second

       if ( present(year) ) this % year = year
       if ( present(month) ) this % month = month
       if ( present(day) ) this % day = day
       if ( present ( hour ) ) this % hour = hour
       if ( present ( minute ) ) this % minute = minute
       if ( present ( second ) ) this % second = second
       
       call this % update()
      
  end subroutine set_date
   

  ! ----------------------------
  !
  ! ----------------------------
  subroutine get_date ( this , year , month, day, hour ,minute, second, doy , hour_frac, msec_of_day)
      class ( date_type) :: this
      integer , optional :: year , month ,day, hour , minute, second, doy
      real, optional :: hour_frac
      integer, optional :: msec_of_day
      
      if ( present(year) )  year = this % year
      if ( present(month) ) month = this % month 
      if ( present(day) ) day = this % day 
      if ( present ( hour ) ) hour = this % hour
      if ( present ( minute ) ) minute =  this % minute 
      if ( present ( second ) ) second = this % second 
      if ( present ( doy ) ) doy = this % dayOfYear
      if ( present ( hour_frac )) hour_frac = this % hour_frac
      if ( present ( msec_of_day )) msec_of_day = this % msec_of_day
   
  end subroutine get_date
  ! --------------------------------------------------
  !    Returns day of year when 16-day period starts
  ! ---------------------------------------------------
  function period_16 (this ) result ( out_string )
      class ( date_type) :: this
      character ( len = 3 ) :: out_string
      integer :: iperiod16
      
       iperiod16 = 16 * (( this % dayofyear-1) / 16  ) + 1
     
      write (  out_string  , fmt = '(i3.3)') iperiod16
   
  end function period_16
  

  !  -- 
  !   poplutes with day of year ( from 1 to 366) instead of month and day
  !  ---
  subroutine set_date_with_doy (this , year , doy, hour , minute, second)
      class ( date_type) :: this
      integer , optional :: year , doy , hour , minute, second
      integer , dimension(12) :: jmonth 
      integer, dimension(12) :: last_day_month
      integer :: i
      integer , dimension(12) :: day_dum
      integer :: month(1)
       
      if ( present(year) ) this % year = year
      
      if ( present(doy) ) this % dayOfYear = doy
      if ( present ( hour ) ) this % hour = hour
      if ( present ( minute ) ) this % minute = minute
      if ( present ( second ) ) this % second = second
       
      jmonth = [31,28,31,30,31,30,31,31,30,31,30,31]
      if  ( modulo ( this % year , 4) == 0 ) jmonth(2) = 29 
      last_day_month(1) = 31
       
      do i = 2 , 12
         last_day_month (i) = last_day_month (i-1) + jmonth(i)  
      end do
      day_dum = doy - last_day_month
      month = maxloc ( day_dum , mask = day_dum <= 0 )
      
      this % month = month(1)
      this % day = jmonth(month(1)) + day_dum(month(1))
      
      call this % update()
       
  end subroutine set_date_with_doy
   
  ! -- ----
  !   Sets Julian Day and updates data
  !-----
  subroutine set_jday ( this, jday)
     class ( date_type) :: this
     real ( kind = r15 ) :: jday
     this % julday = jday
     call this % update_from_jd()
  end subroutine set_jday
 
  !
  !
  !
  function next_6h ( self , count ) result ( out)
      class ( date_type)  :: self
      type (date_type)  :: out
      integer , intent(in), optional :: count
      integer :: cnt
      cnt = 1
      if ( present ( count )) cnt = count
      if (cnt <= 0) cnt = cnt + 1
      
      out = self
      
      out % hour = 6 * ((self % hour )/ 6 + cnt )
      out % minute = 0
      out % second = 0
      call out % update()
   
  end function next_6h
  
  
  function days_since_proj_time0 ( self , verbose ) result(out)
      class ( date_type)  :: self
      real :: out
      
      logical, intent(in), optional :: verbose
      type (date_type)  :: epoch
      
      call epoch % set_date(year=2010,day=1,month=1,hour=0,minute=0)
     
      out = self%julday - epoch%julday
      
      
       if ( present(verbose)) then
         print*,' days since 01/01/2010T00:00:00'
         print*, 'to: '
         print*, self%  date_string('yy/mm/dd/hh')
         print*,out
       endif
  
  end function days_since_proj_time0
  

  !
  !
  !
  function next_3h ( self , count ) result ( out)
      class ( date_type)  :: self
      type (date_type)  :: out
      integer , intent(in), optional :: count
      integer :: cnt
      cnt = 1
      if ( present ( count )) cnt = count
      if (cnt <= 0) cnt = cnt + 1
      
      out = self
      
      out % hour = 3 * ((self % hour )/ 3 + cnt )
      out % minute = 0
      out % second = 0
      call out % update()
   
  end function next_3h
  ! --
  subroutine update ( this )
      class ( date_type ) :: this
        character ( len = 2 ) :: year_s2d
      character ( len = 4 ) :: year_s
      character ( len = 2 ) :: month_s
      character ( len = 2 ) :: day_s
      character ( len = 2 ) :: hour_s
      character ( len = 2 ) :: minute_s
      
      this % julday = julday_cmp (this % year , this % month &
                           , this % day ,this % hour &
                           , this % minute )
       this % dayofyear = 1+  this % julday -  julday_cmp (this % year , 1 &
                           , 1 ,this % hour &
                           , this % minute ) 
                           
                           
      call this % update_from_jd()
      
      write ( year_s2d, fmt ='(i2.2)') mod(this % year , 100)
      write ( year_s, fmt = '(i4.4)') this % year
      write ( month_s, fmt = '(i2.2)') this % month
      write ( day_s, fmt = '(i2.2)') this % day
      write ( hour_s , fmt = '(i2.2)') this % hour
      write ( minute_s , fmt = '(i2.2)') this % minute
      this % yyyy = year_s
      this % mm = month_s
      this % dd = day_s
      this % yymmdd =  year_s2d//month_s//day_s
      this % yyyymmdd =  year_s//month_s//day_s
      this % yymmddhh =  year_s2d//month_s//day_s//hour_s
  end subroutine update
  !
  !
  subroutine update_from_jd ( this )
      class ( date_type ) :: this
      call cdate ( this % julday , this % year , this % month , this % day ) 
      this % hour =  int (24.* ( this % julday - int(this % julday) ))
      this % minute = int((60)  &
             * ((24. * ( this % julday - int(this % julday) )) - this % hour))
      this % hour_frac = this % hour + this % minute / 60.
     
      this % msec_of_day =  60* 60* 1000 * this % hour &
                           + 60* 1000 * this % minute &
                           + 1000 * this % second
                           
  
  end subroutine update_from_jd
   
   
  subroutine set_julday ( this , julday)
      class ( date_type ) :: this
      real ( kind= r15 ) :: julday
      call cdate ( julday , this % year , this % month , this % day ) 
      this % hour =  int (24.* ( julday - int(julday) ))
      this % minute = int((60)  &
             * ((24. * ( julday - int(julday) )) - this % hour))
      call this % update()
   
  end subroutine set_julday
  ! -- 
  ! -- 
  subroutine print_data( this )
      class (date_type) :: this
      
      print*,'year: ', this % year
      print*,'month: ', this % month
      print*,'day: ', this % day
      print*,'hour: ', this % hour
      print*,'minute: ', this % minute
      print*,'second: ', this % second
      print*,'jday: ' , this % julday
      print*, 'frac of hour: ',this % hour_frac
   
  end subroutine print_data
  
  !   ----------
  !     returns string in wished format
  !
  !    EXAMPLE
  !
  !      date_string = t0 % date_string(fmt ='yy/mm/hhmm')
  !
  function date_string ( this , fmt ) result(out)
      class ( date_type) :: this
      character ( len = * ) , intent (in) :: fmt 
      character (len=:) , allocatable :: out
      character ( len = 2 ) :: year_s2d
      character ( len = 4 ) :: year_s
      character ( len = 2 ) :: month_s
      character ( len = 2 ) :: day_s
      character ( len = 2 ) :: hour_s
      character ( len = 2 ) :: minute_s
      character ( len =3 )  :: doy_s
      integer :: len_fmt
      
      len_fmt = len(fmt) 
      ! not supported by gfortran 4.7!
      allocate ( character(len = len_fmt ) :: out  )
      
      write ( year_s2d, fmt ='(i2.2)') mod(this % year , 100)
      write ( year_s, fmt = '(i4.4)') this % year
      write ( month_s, fmt = '(i2.2)') this % month
      write ( day_s, fmt = '(i2.2)') this % day
      write ( hour_s , fmt = '(i2.2)') this % hour
      write ( minute_s , fmt = '(i2.2)') this % minute
      write ( doy_s , fmt = '(i3.3)') this % dayofyear
      
      out='start'  
      
      
      select case (fmt)
         case ('yymmdd')
            out = year_s2d//month_s//day_s
         case ('yymmddhhmm')
           out = year_s2d//month_s//day_s//hour_s//minute_s  
          case ('yymmddhh')
            out = year_s2d//month_s//day_s//hour_s     
         case ('yy/mm/dd')
           out = year_s2d//'/'//month_s//'/'//day_s 
         case ('yy/mm/dd/hh')
           out = year_s2d//'/'//month_s//'/'//day_s//'/'//hour_s   
         case ('yy/mm/hhmm')
           out = year_s2d//'/'//month_s//'/'//hour_s//minute_s 
         case ('yyyy-mm-dd')
          out =  year_s//'-'//month_s//'-'//day_s     
         case ('yyyy')
           out = year_s  
         case ('yyyy_doy')
           out = year_s//'_'//doy_s     
         case ( 'mm')
           out =  month_s
         case  ('dd')
            out = day_s 
          case ('')  
         case default
            out='format not set Consider add it in class_time_date.f90 or contact andi.walther@ssec.wisc.edu'
            print*,'WARNING: ',out, '> ',fmt
      end select
      
      
  end function date_string
  
  !  adds time to data content
  !     EXAMPLE:
  !       to add 2 days and 12 hours:
  !
  !     call t0 % add_time ( day=2, hour=12 )
  !
  subroutine add_time ( this  , day , hour , minute)
      class ( date_type ) :: this
      integer , optional ::   day , hour , minute
      integer ::   day_add, hour_add , minute_add
      real ( kind = r15) :: julday_add
      
      day_add = 0
      hour_add = 0
      minute_add = 0
      
      if ( present ( day) )      day_add = day
      if ( present ( hour ) )    hour_add = hour
      if ( present ( minute ) )  minute_add = minute
      
      julday_add = this % julday + day_add + hour_add/24. + minute_add/(24.*60.)
      call cdate ( julday_add , this % year , this % month , this % day ) 
      this % hour =  int (24.* ( julday_add - int(julday_add) ))
      this % minute = int((60)  &
             * ((24. * ( julday_add - int(julday_add) )) - this % hour))
      call this % update()
  end subroutine add_time
  
    subroutine add_days ( this  , day )
      class ( date_type ) :: this
      real , optional ::   day 
      real ::   day_add 
      real ( kind = r15) :: julday_add
      
      day_add = 0.
     
      
      if ( present ( day) )      day_add = day
      
      
      julday_add = this % julday + day_add 
      call cdate ( julday_add , this % year , this % month , this % day ) 
      this % hour =  int (24.* ( julday_add - int(julday_add) ))
      this % minute = int((60)  &
             * ((24. * ( julday_add - int(julday_add) )) - this % hour))
      call this % update()
  end subroutine add_days
  
   
  
  ! --- these are releated tools but not part of the class
  !
  !   TIME_DIFF_WEIGHT
  !
  !    INPUT: 
  !          time:   Time object for the point in time
  !          time0:  lower time object
  !          time1: upper time object
  !          identical_time_bounds: optional, returns a .TRUE. if time0 and time1 are identicale 
  !
  !   OUTPUT : relativ distance inside the time0 and time1 from [0,1]
  !
  
  function time_diff_weight ( time, time0, time1,identical_time_bounds ) result(wgt)
      type (date_type) , intent(in) :: time, time0 , time1
      logical, optional, intent(out) :: identical_time_bounds 
      real (kind = r15) :: wgt
      
      wgt = -1.
      identical_time_bounds  = .true.
      if ( time1%julday .NE. time0 %julday) then
        wgt = (time%julday - time0%julday)/ (time1%julday - time0 %julday)
        identical_time_bounds = .false.
      end if
   
  end function time_diff_weight
   
  function time_is_in_window ( time, time0, time1 ) 
      type (date_type) , intent(in) :: time, time0 , time1
      logical :: time_is_in_window
      real :: diff
      time_is_in_window = .false.
      
      diff = time_diff_weight ( time, time0, time1 )
      if ( diff .GE. 0. .and.  diff .le. 1) time_is_in_window = .true.
   
  end function time_is_in_window
   
  function time_mid ( time0,time1)
      type (date_type) , intent(in) ::  time0 , time1
      type ( date_type ) :: time_mid
      
      call time_mid % set_jday((time0 % julday + time1 % julday)/2.)
      
  end function time_mid
  
  !
  function julday_cmp(yyyy, mm, dd, hh , uu ) result(ival)
     
      integer, intent(in)  :: yyyy
      integer, intent(in)  :: mm
      integer, intent(in)  :: dd
      integer, intent(in) , optional :: hh , uu 
      real ( kind = r15 )   :: ival

      ival = dd - 32075 + 1461*(yyyy+4800+(mm-14)/12)/4 +  &
       367*(mm-2-((mm-14)/12)*12)/12 - 3*((yyyy+4900+(mm-14)/12)/100)/4
      ival = ival + hh/24. + uu/(24.* 60) 
      return
  end function julday_cmp
  
   SUBROUTINE cdate(jd, yyyy, mm, dd)
      !=======GIVEN A JULIAN DAY NUMBER, NNNNNNNN, YYYY,MM,DD ARE RETURNED AS THE
      !              CALENDAR DATE. JD = NNNNNNNN IS THE JULIAN DATE FROM AN EPOCH
      !              IN THE VERY DISTANT PAST.  SEE CACM 1968 11(10):657,
      !              LETTER TO THE EDITOR BY FLIEGEL AND VAN FLANDERN.
      !    EXAMPLE CALL CDATE(2440588, YYYY, MM, DD) RETURNS 1970 1 1 .

      real ( kind = r15), INTENT(IN)   :: jd
      INTEGER, INTENT(OUT)  :: yyyy
      INTEGER, INTENT(OUT)  :: mm
      INTEGER, INTENT(OUT)  :: dd
      INTEGER :: l, n

      l = jd + 68569
      n = 4*l/146097
      l = l - (146097*n + 3)/4
      yyyy = 4000*(l+1)/1461001
      l = l - 1461*yyyy/4 + 31
      mm = 80*l/2447
      dd = l - 2447*mm/80
      l = mm/11
      mm = mm + 2 - 12*l
      yyyy = 100*(n-49) + yyyy + l
      RETURN
  END SUBROUTINE cdate
  
end module class_time_date



