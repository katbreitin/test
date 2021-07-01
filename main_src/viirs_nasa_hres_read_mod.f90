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
    logical::channel_on_viirs(16)
    character(len=200) :: filename
    character (len=1000) :: path
    integer :: ny_start
    integer :: ny_end
    integer, dimension(16) :: Modis_Chn_List = [ 8 , 9 , 3 , 4 , 1 , 15 , 2 , 5 , 26 , 6 , 7 , 20 , 22 , 29 , 31 , 32 ]
    contains
    procedure :: map_modis_to_viirs 
  end type viirs_nasa_hres_config_type
  
  
  
  
contains


subroutine map_modis_to_viirs(self)
  class(viirs_nasa_hres_config_type) :: self
  
  integer :: i
  
  do i=1,16
       self %  channel_on_viirs (i) = self % channel_on_modis(self % Modis_Chn_List(i))
  end do
  

end subroutine map_modis_to_viirs


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
  
  use Pixel_Common_Mod, only : ch, image
  
  
  
  type ( viirs_nasa_hres_config_type), intent (in) :: in_config
  integer :: status
  character(len=1020) :: File_Local
  real,  allocatable :: out(:,:)
  integer , dimension(2) :: start,count
  character(2) :: ch_str
  integer i_ch
  
  print*,'START:  READ nasa viirs hres'
  print*, trim(in_config % filename)
  call in_config % map_modis_to_viirs ()
  file_local = trim(in_config%Path)//trim(in_config%filename)

  
  start = (/1,in_config % ny_start /)
  count = (/6400,in_config % ny_end- in_config % ny_start + 1 /)
  
  do i_ch =1,11
    if (in_config % channel_on_viirs (i_ch)) then
      write ( ch_str, '(i2.2)' ) i_ch 
      status=cx_sds_read(file_local,'observation_data/M'//ch_str//'_highres',out,start = start,count = count)
      ch(in_config % modis_chn_list(i_ch)) % ref_toa = out
    end if
  end do

  
  do i_ch =12,16
      if (in_config % channel_on_viirs (i_ch)) then
        write ( ch_str, '(i2.2)' ) i_ch 
        
        status=cx_sds_read(file_local,'observation_data/M'//ch_str//'_highres',out,start = start,count = count)
       
        ch(in_config % modis_chn_list(i_ch)) % rad_toa = out
        
    end if
  end do
  
  
  
end subroutine read_viirs_nasa_hres_data

end module VIIRS_NASA_HRES_READ_MOD
