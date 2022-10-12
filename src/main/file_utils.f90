! name:                      file_utils
! function:                   module which houses routines to perform basic file operations
! description:
! reference:
! calling sequence:
! inputs:
! outputs:
! dependencies:
! restrictions:
! history:                              added  Jan 2013 (AW)
!
!     Add file_utilities tools (June 2019 AW)
!
!-----------------------------------------------------------------------------------------------------------------------

module file_utils
  use pixel_common_mod, only: Temporary_Data_Dir

  implicit none


  private

  public :: file_basename
  public :: file_dirname
  public :: file_search
  public :: file_test
  public :: get_lun
  public :: file_nr_lines
  public :: uncompress_file

  

contains

 !------
 ! name:                      file_test
! function:                   function checks files for existence and other attributes without having to first open the file. 
! description:                FILE_TEST returns .true., if the specified file exists
! reference:
! calling sequence:         result = file_test ( <file> )
! inputs:                       A scalar if filename to be tested
! outputs:                      A logical variable
! dependencies:            none
! restrictions:                none
! history:                              added  Jan 2013 (AW)
 function file_test (file) result (existence)
   implicit none
   character ( *), intent(in) :: file
   logical :: existence
   inquire ( FILE = file, EXIST = existence ) 
  end function file_test

 !------
 ! name:                      file_nr_lines
! function:                   returns number of line in an ASCII file
! description:
! reference:
! calling sequence:
! inputs:
! outputs:
! dependencies:
! restrictions:                    private subroutine of this module
! history:                              added  Jan 2013 (AW)
  function file_nr_lines (file) result (return_value)
    implicit none
    integer :: nr, ios, j
    integer :: maxrecs = 1000
    character (len = *) , intent(in) :: file 
    integer :: return_value
    character(len=1) :: junk
    integer :: lun

    lun = get_lun()

    nr = 0
    open( UNIT=lun,FILE=file)
      do j = 1 , maxrecs
        read(lun,*,IOSTAT=ios) junk
 
        if (ios /= 0) exit
        if (j == maxrecs) then
          print*,'exceeded mac lines'
        endif
        nr = nr + 1
      end do
    close(unit=lun)
    return_value = nr 
  end function file_nr_lines
  
!-----
! name:                      file_basename
! function:                  The FILE_BASENAME function returns the basename of a file path. A file path is a string containing one or more segments consisting of 
!                                              names separated by directory delimiter characters (slash (/) under UNIX 
! description:
! reference:
! calling sequence:            Result = FILE_BASENAME( <file> ) 
! inputs:                            A scalar if filename to be tested
! outputs:
! dependencies:
! restrictions:
! history:                              added  Jan 2013 (AW)
  function file_basename(file) result(return_string)
    implicit none
    character(1020) , intent(in) :: file
    character(1020) :: return_string
    character :: sep = '/'
    integer :: kdot
    logical , parameter :: backward = .true.
 
    kdot = scan (file,sep, backward)
    return_string = file(kdot+1: )
  end function file_basename 

!-----
! name:                      file_dirname
! function:                  
! description:                 The FILE_DIRNAME function returns the dirname of a file path. 
! reference:
! calling sequence:             Result = FILE_DIRNAME( <file> ) 
! inputs:                            A scalar if filename to be tested
! outputs:
! dependencies:
! restrictions:
! history:                              added  Jan 2013 (AW)
  function file_dirname(file) result(return_string)
    implicit none
    character(1020) , intent(in) :: file
    character(1020) :: return_string
    character :: sep = '/'
    integer :: kdot
    logical , parameter :: backward = .true.
 
    kdot = scan(file,sep, backward)
    return_string = file(1: kdot -1 )
  end function file_dirname

! name:                      file_utils
! function:                   
! description:              The FILE_SEARCH function returns a string array containing the names of all files 
!                                matching the input path specification. Input path specifications may contain wildcard characters, 
!                                   enabling them to match multiple files.
! reference:
! calling sequence:
! inputs:                         path: a string which holds the pathname
!                                   spec: the search pattern, may include wildcard ( *,?)
!                                   rel_path:  if set the result hold salso the relative path from the working directory,
!                                                  default is the files' basename
! outputs:                       count : integer which holds the dimesnion of the returned result
! restrictions:                    tested only on UNIX and debian environment
! history:                              added  Jan 2013 (AW)
!                                   set file "list" to unique_dummy_file which is unique...

   function file_search ( path, spec , count , rel_path ) result(return_string)

      implicit none
      character(*) , intent(in) :: spec
      character(*) , intent(in) :: path
      integer , intent(out), optional :: count
      logical , intent(in) , optional :: rel_path
      character(len=1020), pointer, dimension(:) :: return_string
      character(len=1020) :: cfile
      character(len=4096) :: cmd
      integer :: nr, ii, ierr, nc
      integer :: lun
      character(len=8) :: date
      character(len=20) :: time
      character(len=2048) :: unique_dummy_file
      
      call date_and_time(date = date, time=time)
      unique_dummy_file = trim(Temporary_Data_Dir)//'/fort.file_search_dummy_'//trim(time)  

      cmd = trim(unique_dummy_file); nc = len_trim(cmd)
      call univ_remove_f(nc, trim(cmd), ierr)
      lun = get_lun()
     
      cmd = 'ls -1 -phd '// trim (path) //''// trim (spec) //' > '//trim(unique_dummy_file)//' 2>/dev/null'
      nc = len_trim(cmd)
      call univ_system_cmd_f(nc, trim(cmd), ierr)
      nr = file_nr_lines (trim(unique_dummy_file))
      open(unit = lun , file = trim(unique_dummy_file) )
      allocate (return_string(nr))
       
      do ii = 1, nr 
         read(lun,"(A)") cfile
        
         return_string(ii) = trim (file_basename (cfile) )
         if ( present ( rel_path ))  return_string(ii) = trim (cfile)    
      end do
      
      close(lun)
      cmd = trim(unique_dummy_file); nc = len_trim(cmd)
      call univ_remove_f(nc, trim(cmd), ierr)
      if  ( present ( count ) ) count = nr
   
  end function file_search
  
   !
   !
   !
  subroutine uncompress_file ( file_original ,   file_unzipped , dir_in,  tmp_dir_in)

      character ( len = * ) , intent(in) :: file_original
      character ( len = * ) , intent(in) , optional :: dir_in
      character ( len = * ) , intent(out) :: file_unzipped
      character ( len = * ) , intent(in), optional :: tmp_dir_in
      
      character ( len =1020) :: dir
      character ( len =1020) :: tmp_dir
      
      character(len=4096) :: cmd
      integer :: len_orig, ierr, nc
      character ( len = 3) :: last_3_ch
      
      tmp_dir = "./temp/"
      dir ="./"
      if (  present ( tmp_dir_in )) tmp_dir = tmp_dir_in
      if (  present ( dir_in )) dir = dir_in
      
      len_orig = len_trim ( file_original)
      last_3_ch = file_original ( len_orig -2 : len_orig)
      
      select case ( last_3_ch ) 
      
      case ( '.gz')     
         file_unzipped = file_original ( 1:len_orig-3)
         cmd = "gunzip -c "//trim(dir)//trim(file_original)// &
              " > "//trim(tmp_dir)//trim(file_unzipped)  
            
         nc = len_trim(cmd)
         call univ_system_cmd_f(nc, trim(cmd), ierr)
         file_unzipped=trim(tmp_dir)//trim(file_unzipped) 
            
      case ('bz2')
         file_unzipped = file_original ( 1:len_orig-4)
         cmd = "bunzip2 -c "//trim(dir)//trim(file_original)// &
              " > "//trim(tmp_dir)//trim(file_unzipped)   
         nc = len_trim(cmd)
         call univ_system_cmd_f(nc, trim(cmd), ierr)
         file_unzipped=trim(tmp_dir)//trim(file_unzipped)
      case default
         file_unzipped = trim(dir)//trim(file_original)
      
      end select
      
       
   
   end subroutine uncompress_file
   
   
   FUNCTION get_lun() RESULT( lun )


    ! -----------------
    ! Type declarations
    ! -----------------
 
    INTEGER :: lun
    LOGICAL :: file_open


    ! --------------------------------------------
    ! Initialise logical unit number and file_open
    ! --------------------------------------------

    lun = 9
    file_open = .TRUE.


    ! ------------------------------
    ! Start open loop for lun search
    ! ------------------------------

    lun_search: DO

      ! -- Increment logical unit number
      lun = lun + 1

      ! -- Check if file is open
      INQUIRE( lun, OPENED = file_open )

      ! -- Is this lun available?
      IF ( .NOT. file_open ) EXIT lun_search

    END DO lun_search

  END FUNCTION get_lun
 
  

end module file_utils
