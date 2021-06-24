!$Id: viirs_nasa_hres_read_module.f90 
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version
! 5.4
!
! NAME: viirs_nasa_read_module.f90
!
! PURPOSE: VIIRS NASA (NetCDF4) read tool
!
! DESCRIPTION:  This module deals with reading VIIRS NASA High resolution data
!
! AUTHORS:
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
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
! HISTORY:   created        June 2021 (AW)
!-------------------------------------------------------------------------------

module VIIRS_NASA_HRES_READ_MOD
  implicit none
  
contains

subroutine read_viirs_nasa_hres_data
  print*,'hallo nasa viirs hres'
  stop

end subroutine read_viirs_nasa_hres_data

end module VIIRS_NASA_HRES_READ_MOD
