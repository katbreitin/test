! - $Header: https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/trunk/main_src/clavrx_message_mod.f90 3082 2018-12-17 17:53:19Z mfoster $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: clavrx_message_mod.f90 (src)
!       clavrx_message_mod (program)
!
! PURPOSE: 
! Print message to screen based on verbose level
!
! DESCRIPTION: 
! This module houses the routines that can be called to print message to screen.
! Messages will print if the verbose level is equal or less than the
! preset verbose flag from options file or verbose_level.txt. If no verbose level 
! is set, a default value of 5 is used.
! Use MESG for string only output, MESG_1R or MESG_1D when both string and real number
! are needed, and MESG_1I for string and integer 
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
! Notes
!  AKH modified to use Verbose_Level_Flag which is read from clavrx options file
!--------------------------------------------------------------------------------------
module CLAVRX_MESSAGE_MOD

  use foul,only: write_formatted
  use PIXEL_COMMON_MOD, only: Verbose_Level_Flag

  implicit none
  private

  public :: mesg, mesg_1i
  public :: verb_lev

  type verbose_type
    integer :: QUIET   =  0
    integer :: ERROR =  1
    integer :: MINIMAL =  2
    integer :: WARNING =  4
    integer :: DEFAULT =  5
    integer :: VERBOSE =  9 
  end type verbose_type

  type (verbose_type ) , save :: verb_lev 

  character(len = *) , parameter :: PROMPT = 'CLAVR-x >'     !can we standardize this?
  integer :: VERBOSE_LEVEL = 5  ! from quiet (0) to verbose(10)

  interface mesg
    module procedure mesg_pure
    module procedure mesg_1r
    module procedure mesg_1d
  end interface mesg


contains


  !---------------------------------------------------------------------
  ! Add Comments Here
  !---------------------------------------------------------------------
  subroutine do_it ( text, message_level ,color_in )
    use file_tools, only: file_test
    
    character ( len = * ) , intent (in) :: text 
    integer , intent ( in ) :: message_level
    character ( len=*), intent(in),optional :: color_in
    character (len=20) :: color
!--- removed by akh (Verbose_Level_Flag read in options file)
!   integer :: verbose_level
      
!   !--- please don't fixed unit numbers
!   verbose_level = verb_lev % DEFAULT
!   if ( file_test ( 'verbose_level.txt') ) then
!     open ( 37, file = 'verbose_level.txt' )
!     read (37,'(i1)') verbose_level
!     close (37)
!   end if
!-------------------------------------------------------------
      
    color = 'black'
    if ( present(color_in)) color=color_in
        
    if ( message_level <= Verbose_Level_Flag ) then
      
      !--- what is verbose_level = 6, not in structure
      if ( verbose_level .eq. 6 ) then
        call write_formatted('CLAVR-x:','black underline background_white',' '//trim(text),color)
      else
        print*,'CLAVR-x:  '//trim(text)
      end if
    end if

   end subroutine do_it

  !---------------------------------------------------------------------
  ! Add Comments Here
  !---------------------------------------------------------------------
   subroutine mesg_pure ( text, level , color )
      character (len = *) , intent(in) :: text
      integer, optional, intent(in) :: level
      integer, optional, intent(in) :: color
      character( len = 5) :: color_string
      integer :: lev
      
      lev = verb_lev % DEFAULT
      if ( present (level)) lev = level
      
      color_string = 'black'
      
      if (present(color)) then
         select case (color)
         case(1) 
            color_string='red'
         case(2)
            color_string='gray'
         case(3)
            color_string='blue'
         case(4)
            color_string = 'green' 
         case(5)
            color_string = 'magenta'           
         case default
            color_string='black'
         end select
      end if
      
      call do_it ( trim(text), lev,color_string ) 
      
   end subroutine mesg_pure

  !---------------------------------------------------------------------
  ! Add Comments Here
  !---------------------------------------------------------------------
   subroutine mesg_1r ( text,  param_r, level , color )
      character (len = *) , intent(in) :: text
      real, intent(in) :: param_r
      integer, optional, intent(in) :: level
      integer, optional, intent(in) :: color
      character( len = 2) :: color_string
      integer :: lev
      character ( len =100 ) :: string_100
      character ( len =200) :: text_1
 
  
      write ( string_100, '(f20.4)') param_r
  
      lev = verb_lev % DEFAULT
      if ( present (level)) lev = level
      text_1 = text//trim(string_100)
      color_string=''  
      if (present(color)) write(color_string,'(I2)') color 
      
      call do_it ( trim(text_1), lev ) 

   end subroutine mesg_1r

  !---------------------------------------------------------------------
  ! Add Comments Here
  !---------------------------------------------------------------------
   subroutine mesg_1d ( text,  param_d, level , color )
      character (len = *) , intent(in) :: text
      real(8), intent(in) :: param_d
      integer, optional, intent(in) :: level
      integer, optional, intent(in) :: color
      character( len = 2) :: color_string
      integer :: lev
      character ( len =100 ) :: string_100
      character ( len =200) :: text_1


      write ( string_100, '(f25.4)') param_d

      lev = verb_lev % DEFAULT
      if ( present (level)) lev = level
      text_1 = text//trim(string_100)
      color_string=''
      if (present(color)) write(color_string,'(I2)') color

      call do_it ( trim(text_1), lev )

   end subroutine mesg_1d

  !---------------------------------------------------------------------
  ! Add Comments Here
  !---------------------------------------------------------------------
   subroutine mesg_1i ( text,  param_i, level , color )
      character (len = *) , intent(in) :: text
      integer, intent(in) :: param_i
      integer, optional, intent(in) :: level
      integer, optional, intent(in) :: color
      character( len = 2) :: color_string
      integer :: lev
      character ( len =100 ) :: string_100
      character ( len =200) :: text_1

      write ( string_100, '(i10)') param_i

      lev = verb_lev % DEFAULT
      if ( present (level)) lev = level
      text_1 = text//trim(string_100)
      color_string=''
      if (present(color)) write(color_string,'(I2)') color

      call do_it ( trim(text_1), lev )

   end subroutine mesg_1i

end module CLAVRX_MESSAGE_MOD
