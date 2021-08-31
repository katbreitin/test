!$Id: viirs_nasa_hres_read_module.f90 
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version
! 5.4
!
! NAME: viirs_nasa_hres_read_module.f90
!
! PURPOSE: VIIRS NASA (NetCDF4) read tool for high resolution files
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
    character(len=50) :: sensor
    character(len=200) :: filename
    character(len=20) :: time_identifier
    character (len =200) :: iband_geo_filename
    character (len=1000) :: path
    integer :: ny_start
    integer :: ny_end
    integer, dimension(16) :: Modis_Chn_List = [ 8 , 9 , 3 , 4 , 1 , 15 , 2 , 5 , 26 , 6 , 7 , 20 , 22 , 29 , 31 , 32 ]
    contains
    procedure :: map_modis_to_viirs 
  end type viirs_nasa_hres_config_type
  
  
  type viirs_nasa_hres_out_type
    logical :: is_set
  end type
  
  type viirs_coef_type
    logical :: is_set
    character(len=1020) :: file
    character(len=50) :: sensor
    real :: solar_ch20
    real :: ew_ch20
    real :: solar_ch20_nu
    real :: planck_a1(20:43)
    real :: planck_a2(20:43)
    real :: planck_nu(20:43)
    real :: VIIRS_Correction_Factor(11)
    contains
    procedure :: read_file => viirs_coef_type__read_file
  end type
  
  type(viirs_coef_type) :: coef
  
  logical :: first_run = .true.
  
  
  type viirs_iband_bowtie_type
    logical :: is_set
    integer :: nx_pattern = 6400
    integer :: ny_pattern = 32 
    integer :: mask (6400,32)
    integer :: idx (6400,32)
    
    integer :: n_lines_gap_short = 1280
    integer :: n_lines_gap_long = 2016
    contains
    procedure :: set => viirs_iband_bowtie_type__set
  end type
  
  type ( viirs_iband_bowtie_type ) :: bowtie_coef
  
  
contains

subroutine viirs_coef_type__read_file ( self, sensor)
  class(viirs_coef_type) :: self
  character(len=*) :: sensor
   integer:: ios0, erstat
   character(len=20):: header
   integer :: lun_id
   
  if ( self % is_set .and. trim(self % sensor) .eq. trim(sensor)) return
  self % file = '/DATA/Ancil_Data/clavrx_ancil_data/static/clavrx_constant_files/viirs_npp_instr.dat'
  print*,'reads it ', trim(self % file)
  lun_id = 12
  open(unit=lun_id,file=trim( self % file),status="old",position="rewind",action="read",iostat=ios0)
      read(unit=lun_id,fmt="(a3)") self % sensor
      read(unit=lun_id,fmt=*)  self % Solar_Ch20
      read(unit=lun_id,fmt=*)  self % Ew_Ch20
      read(unit=lun_id,fmt=*) header
      read(unit=lun_id,fmt=*)  self % Planck_A1(20),  self % Planck_A2(20), self % Planck_Nu(20)
      read(unit=lun_id,fmt=*)  self % Planck_A1(22),  self % Planck_A2(22), self % Planck_Nu(22)
      read(unit=lun_id,fmt=*)  self % Planck_A1(29),  self % Planck_A2(29), self % Planck_Nu(29)
      read(unit=lun_id,fmt=*)  self % Planck_A1(31),  self % Planck_A2(31), self % Planck_Nu(31)
      read(unit=lun_id,fmt=*)  self % Planck_A1(32),  self % Planck_A2(32), self % Planck_Nu(32)
      read(unit=lun_id,fmt=*)  self % Planck_A1(42),  self % Planck_A2(42), self % Planck_Nu(42)
      read(unit=lun_id,fmt=*)  self % Planck_A1(43),  self % Planck_A2(43), self % Planck_Nu(43)
      read(unit=lun_id,fmt=*) header
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(1)
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(2)
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(3)
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(4)
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(5)
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(6)
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(7)
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(8)
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(9)
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(10)
      read(unit=lun_id,fmt=*)  self % VIIRS_Correction_Factor(11)
  close(unit=lun_id)
  
  self % is_set = .true. 
  
end subroutine viirs_coef_type__read_file



subroutine map_modis_to_viirs(self)
  class(viirs_nasa_hres_config_type) :: self
  
  integer :: i
  
  do i=1,16
       self %  channel_on_viirs (i) = self % channel_on_modis(self % Modis_Chn_List(i))
  end do
  

end subroutine map_modis_to_viirs

subroutine viirs_hres_date


end subroutine viirs_hres_date


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
  
  use Pixel_Common_Mod, only : ch, image , nav, geo 
  
  use VIEWING_GEOMETRY_MOD, only: &
        GLINT_ANGLE &
      , SCATTERING_ANGLE &
      , RELATIVE_AZIMUTH
  
  use Planck_mod
  
  use file_tools
  
  type ( viirs_nasa_hres_config_type), intent (in) :: in_config
  integer :: status
  character(len=1020) :: File_Local
  real,  allocatable :: out(:,:)
  real,  allocatable :: out1d(:)
  integer , dimension(2) :: start,count
  character(2) :: ch_str
  integer i_ch
  integer :: modis_ch
  real :: noaa_nasa_correct
  character(len=1024), dimension(1) :: file_v03img, file_02img
  character(len=15) :: time_identifier
  integer :: cc
  logical :: rel_path = .TRUE.
  
  if ( first_run) then
    print*,'START:  READ nasa viirs hres ',trim(in_config % filename)
    
    ! make bow tie pattern mask
    call bowtie_coef % set
    
  end if
  ! time identifier
  time_identifier =  in_config%filename(9:23)
  
  call in_config % map_modis_to_viirs ()
  file_local = trim(in_config%Path)//trim(in_config%filename)
   call coef % read_file(in_config % sensor)
  
  start = (/1,in_config % ny_start - 1 /)
  count = (/6400,in_config % ny_end- in_config % ny_start + 1 /)
  
  do i_ch =1,11
   
    if (in_config % channel_on_viirs (i_ch)) then
      write ( ch_str, '(i2.2)' ) i_ch 
      status=cx_sds_read(file_local,'observation_data/M'//ch_str//'_highres',out,start = start,count = count)
      
      ch(in_config % modis_chn_list(i_ch)) % ref_toa(:,1:count(2)) = out
      call update_bowtie(in_config % ny_start, ch(in_config % modis_chn_list(i_ch)) % ref_toa(:,1:count(2)))
     
    end if
  end do

  
  do i_ch =12,16
     
      if (in_config % channel_on_viirs (i_ch)) then
        write ( ch_str, '(i2.2)' ) i_ch 
        modis_ch = in_config % modis_chn_list(i_ch)
        status=cx_sds_read(file_local,'observation_data/M'//ch_str//'_highres',out,start = start,count = count)
        ch(modis_ch) % rad_toa(:,1:count(2)) = out
        call update_bowtie(in_config % ny_start, ch(modis_ch) % rad_toa(:,1:count(2)))
        ! - convert to radiance to NOAA unit.. 
        noaa_nasa_correct = ((10000.0 / coef % planck_nu(modis_ch) ** 2)/10. )
        ch(modis_ch) % rad_toa =  ch(modis_ch) % rad_toa * noaa_nasa_correct
        
        ! -  compute BT
        call COMPUTE_BT_ARRAY (ch ( modis_ch ) % bt_toa , ch ( modis_ch ) % rad_toa , modis_ch, -999.)
        
    end if
  end do
  
  !- TODO : BOWTIE
  
 ! print*,'++++++++++++++  TO-DO make correct VJ! file ',__FILE__,' ' ,__LINE__
  
  file_v03img = file_search(trim(in_config % path),'VNP03IMG'//trim(time_identifier)//'*.nc',cc,rel_path)
  
  if (cc .ne. 1 ) then
    print*,'Missing VNP03IMG file'
  end if
  
  if ( first_run) print*,trim(file_v03img(1))
   
   status=cx_sds_read(trim(file_v03img(1)),'geolocation_data/longitude',out,start = start,count = count)
    Nav % Lon_1b(:,1:count(2)) = out
  
   status=cx_sds_read(trim(file_v03img(1)),'geolocation_data/latitude',out,start = start,count = count)
   Nav % Lat_1b(:,1:count(2)) = out
 
   status=cx_sds_read(trim(file_v03img(1)),'geolocation_data/sensor_azimuth',out,start = start,count = count)
   geo % sataz(:,1:count(2)) = out
 
   status=cx_sds_read(trim(file_v03img(1)),'geolocation_data/sensor_zenith',out,start = start,count = count)
   geo % satzen(:,1:count(2)) = out
  
   status=cx_sds_read(trim(file_v03img(1)),'geolocation_data/solar_azimuth',out,start = start,count = count)
   geo % solaz(:,1:count(2)) = out
  
   status=cx_sds_read(trim(file_v03img(1)),'geolocation_data/solar_zenith',out,start = start,count = count)
  
   geo % solzen(:,1:count(2)) = out
  
   file_02img = file_search(trim(in_config % path),'VNP02IMG'//trim(time_identifier)//'*.nc',cc,rel_path)
   !status=cx_sds_read(trim(in_config % path)//'VNP02IMG.A2020118.0000.001.2020118052345.uwssec.nc','scan_line_attributes/scan_start_time',out1d)
   !  set scan time according nasa viirs
   
  ! print*,'++++++++++++++  TO-DO make correct VJ! file ',__FILE__,' ' ,__LINE__
   
  ! print*,'put several things outside nasa hres read routine in future: ',__FILE__,' Line: ',__LINE__
  
   geo % relaz = RELATIVE_AZIMUTH (Geo%Solaz, Geo%Sataz)
   geo % Scatangle = SCATTERING_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz)
    
   if ( allocated(out)) deallocate(out)
   
   first_run = .false.
    

end subroutine read_viirs_nasa_hres_data


subroutine viirs_iband_bowtie_type__set ( this)
  class(viirs_iband_bowtie_type) :: this
  integer :: i0,i1,i2,i3
  
  ! this makes gap pattern
  
  this % mask = 0
  this % idx = -999
  
  this % mask(1:this % n_lines_gap_long,1:2) = 1
  this % idx(1:this % n_lines_gap_short,1:2) = 5
  this % idx(this % n_lines_gap_short+1:this % n_lines_gap_long,1:2) = 3
  
  this % mask(1:this % n_lines_gap_short,3:4) = 1
  this % idx(1:this % n_lines_gap_short,3:4) = 5
  
  this % mask(1:this % n_lines_gap_long,31:32) = 1
  this % idx(1:this % n_lines_gap_short,31:32) = 28
  this % idx(this % n_lines_gap_short+1:this % n_lines_gap_long,31:32) = 30
  
  this % mask(1:this % n_lines_gap_short,29:30) = 1
  this % idx(1:this % n_lines_gap_short,29:30) = 28
  
  i0 = this % nx_pattern - this % n_lines_gap_long 
  i1 = this % nx_pattern - this % n_lines_gap_short - 1
  i2 = this % nx_pattern - this % n_lines_gap_short
  i3 = this % nx_pattern
  
  this % mask(i0:i3,1:2) = 1
  this % idx(i2:i3,1:2) = 5
  this % idx(i0:i1,1:2) = 3
  
  this % mask(i2:i3,3:4) = 1
  this % idx(i2:i3,3:4) = 5
  
  

  
  this % mask(i0:i3,31:32) = 1
  this % idx(i2:i3,31:32) = 28
  this % idx(i0:i1,31:32) = 30
  
  this % mask(i2:i3,29:30) = 1
  this % idx(i2:i3,29:30) = 28
  
  this % is_set = .true.

end subroutine viirs_iband_bowtie_type__set


subroutine update_bowtie (line0_full, variable )
  
  use Pixel_Common_Mod, only :  gap_pixel_mask
  
  integer :: line0_full
  real , intent(inout) :: variable(:,:)
  
  integer :: number_of_lines
  integer :: sh_var (2)
  integer :: line_in_seg
  integer :: line_in_pattern
  integer :: x_idx
  integer :: line_offset
  integer :: gap_line_idx(6400)
  integer :: line_idx
  
  
  sh_var = shape(variable)
  number_of_lines = sh_var(2)
  
  

  do line_in_seg = 1, number_of_lines
    
    line_in_pattern = mod ( line0_full - 1 + line_in_seg , bowtie_coef % ny_pattern)
    if (line_in_pattern == 0) line_in_pattern =  bowtie_coef % ny_pattern
    
    ! check if line needs modifications/ has bowtie pixels
    if ( bowtie_coef % idx (1, line_in_pattern) .EQ. -999 ) cycle  
    
    line_offset = line_in_seg - line_in_pattern
    gap_line_idx = bowtie_coef % idx(:,line_in_pattern) + line_offset
   
    where ( gap_line_idx .LE. 0 )
       gap_line_idx = 1
    end where
    
    where ( gap_line_idx .GT. number_of_lines )
       gap_line_idx = number_of_lines
    end where 
    
    gap_pixel_mask (:,line_in_seg) = bowtie_coef % mask (:,line_in_pattern)
    
    do x_idx = 1, bowtie_coef % nx_pattern
     
      if ( bowtie_coef % mask (x_idx,line_in_pattern) == 1 ) then
        
        
        line_idx =   gap_line_idx(x_idx)
       
        variable (x_idx,line_in_seg) = variable ( x_idx,line_idx) 
       
      end if
    end do
  
  end do 
  
 

end subroutine update_bowtie


end module VIIRS_NASA_HRES_READ_MOD
