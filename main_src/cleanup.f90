!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE
!
! NAME: cleanup.f90 (src)
!
! PURPOSE: Handle signals before exiting
!
!
! AUTHORS:
!  Coda Phillips, coda.phillips@wisc.edu
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
! REVISON HISTORY:
!    Creation Date Apr 2022
!--------------------------------------------------------------------------------------

module cleanup
  use pixel_common_mod, only: Temporary_Data_Dir

  contains
      subroutine cleanup_tempdir()
        integer(4) :: ret
       !--- remove directory for temporary files
        print*, 'Cleaning up'
        call system("rmdir "//trim(Temporary_Data_Dir), ret)
        if(ret .ne. 0) then 
            print*, 'rmdir error'
        else
            print*, 'removed ',trim(Temporary_Data_Dir)
        endif
        stop 1
    end subroutine
end module cleanup
