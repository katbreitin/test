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
  
  type viirs_nasa_hres_config_type
    logical::channel_on_modis(50)
    logical::channel_on_viirs(50)
    character(len=200) :: filename
    character (len=1000) :: path
  end type viirs_nasa_hres_config_type
  
contains

!   this routine is supposed to read all reflectance and radiance data
!   input is filename
!   and options what to read in
!
!
!
!
!
subroutine read_viirs_nasa_hres_data (in_config)
  use cx_sds_io_mod, only: &
           cx_sds_finfo &
         , cx_sds_varinfo &
         , cx_sds_read &
         , MAXNCNAM
  
  type ( viirs_nasa_hres_config_type) :: in_config
  integer :: status
  integer :: ftype      
  integer :: nsds
  integer :: Natt
  character ( len = MAXNCNAM), allocatable :: Sds_Name(:)
  character ( len = MAXNCNAM), allocatable :: Att_Name(:)
  character(len=1020) :: File_Local
  real,  allocatable :: out(:,:)
  
  print*,'hallo nasa viirs hres'
  print*, trim(in_config % filename)
  file_local = trim(in_config%Path)//trim(in_config%filename)
  !status = cx_sds_finfo (File_Local, ftype,nsds,Sds_Name,Natt,Att_Name)
  !print*,status
  !print*,att_name
  !print*,sds_name
  
  status=cx_sds_read(file_local,'observation_data/M07_highres',out)
  print*,out(100:120,200)
  print*,maxval(out)
  stop

end subroutine read_viirs_nasa_hres_data

end module VIIRS_NASA_HRES_READ_MOD
