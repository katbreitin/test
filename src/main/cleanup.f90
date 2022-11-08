!------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE
!
! NAME: cleanup.f90 (src)
!
! PURPOSE: Handle signals and temporary file cleanup before exiting
!
!
! AUTHORS:
!  Coda Phillips, coda.phillips@wisc.edu
!  Tim Michaels, tmichaels@wisc.edu
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
!------------------------------------------------------------------------------

module cleanup

implicit none


contains


!```````````````````````````````````````````````````````````````````
subroutine cleanup_tempdir

use pixel_common_mod, only: Temporary_Data_Dir

implicit none

!=== Argument declarations:

! None

!== Local declarations:

integer :: ierr, nc
character(len=4096) :: cmd

!=== Executable statements:

!--- remove directory for temporary files
print *, 'Cleaning up'
nc = len_trim(Temporary_Data_Dir)
call univ_remove_f(nc, trim(Temporary_Data_Dir), ierr)
if (ierr /= 0) then
   print *, 'rmdir error'
else
   print *, 'removed ', trim(Temporary_Data_Dir)
end if

end subroutine cleanup_tempdir

subroutine cleanup_tempdir__exit
    call cleanup_tempdir()
    call univ__exit(1)
end subroutine cleanup_tempdir__exit


end module cleanup
