   ! $Id: cx_read_ahi_mod.f90 3676 2020-01-15 19:41:25Z stevew $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: cx_read_ahi_mod.f90 ( module)
!       viirs_clavrx_bridge (program)
!
! PURPOSE: AHI Reader for NCDF4 / HDF5 files generated at CIMSS 
!
! DESCRIPTION: 
!
!  DEPENDENCIES: 
!     MODULES:
!           readhdf5dataset
!           string_functions
!           viewing_geometry_module
!           date_tools_mod
!     SUBROUTINES
!      fgf_to_earth ( c-routine )
!
!
!
! AUTHORS:
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
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
! REVISION HISTORY:   Created 6 Feb 2015 (AW )
!
! NOTES:
!
! AHI Channel Mapping
!
! wvl       ahi  modis/clavrx 
!
! 0.47       1      3   
! 0.51       2      4
! 0.64       3      1  
! 0.86       4      2  
! 1.6        5      6
! 2.2        6     r7
! 3.9        7     20
! 6.2        8     37
! 6.9        9     27
! 7.3       10     28
! 8.6       11     29
! 9.6       12     30
! 10.4      13     38
! 11.2      14     31
! 12.3      15     32 
! 13.3      16     33
!
!--------------------------------------------------------------------------------------

module CX_READ_AHI_MOD 
   
   use class_time_date, only : &
      date_type 
   
   implicit none  
   private
   public :: get_ahi_data
   public :: ahi_time_from_filename
   public :: ahi_hcast_time_from_filename
   public :: ahi_segment_information_region
   public :: get_var_dimension
   
   integer, parameter :: NUM_CHN = 16
   
   type, public :: ahi_config_type
      character ( len = 1020 ) :: file_base
      character ( len = 1020 ) :: data_path
      logical :: chan_on(NUM_CHN)
      character (len = 1020) :: filename ( NUM_CHN)
      character ( len =20) :: varname (NUM_CHN)
      integer :: h5_offset(2)
      integer :: h5_count(2)
      integer :: h5_offset_3(2)
      integer :: h5_count_3(2)
      integer :: h5_offset_124(2)
      integer :: h5_count_124(2)
      real :: lon_range (2)
      real :: lat_range (2)
      logical :: filenames_set = .false.
      logical :: do_res3 = .false.
      logical :: do_res124 = .false.
      logical :: do_lonlat_only = .false.
      contains
      procedure :: clean => clean_config
   end type ahi_config_type
   
   type ahi_chn_type
      logical :: is_read
      logical :: is_solar_channel
      real,dimension (:,:) , allocatable :: ref
      real,dimension (:,:) , allocatable :: rad
      real,dimension (:,:) , allocatable :: bt
   end type ahi_chn_type
   
   type :: geo_str
      logical :: is_set
      real , dimension (:,:) , allocatable :: solzen
      real , dimension (:,:) , allocatable :: satzen
      real , dimension (:,:) , allocatable :: solaz
      real , dimension (:,:) , allocatable :: sataz  
      real , dimension (:,:) , allocatable :: relaz 
      real , dimension (:,:) , allocatable :: lat
      real , dimension (:,:) , allocatable :: lon  
      real , dimension (:,:) , allocatable :: glintzen
      real , dimension (:,:) , allocatable :: scatangle
      real , dimension (:)   , allocatable :: scan_time
      logical, dimension (:,:), allocatable :: is_space 
      contains
      procedure :: deallocate_geo
      procedure :: allocate_geo    
   end type  geo_str
   
   type, public :: ahi_data_out_type
      logical :: set = .false.
      type ( ahi_chn_type ) , allocatable :: chn (:)
      type ( geo_str ) :: geo
      type (geo_str ) :: geo_3
      type ( geo_str ) :: geo_124
      type ( date_type ) :: time_start_obj
      type ( date_type ) :: time_end_obj
      logical :: do_res3 = .false.
      logical :: do_res124 = .false.
      logical :: success
      contains
      procedure :: deallocate_all
    
   end type ahi_data_out_type
   
   ! - save this to have it load only once
   type ( ahi_data_out_type ) :: out_fd
   
          
contains
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine get_ahi_data ( config , out , only_nav )
      type ( ahi_config_type ) :: config
      type ( ahi_data_out_type ) :: out
      logical, optional, intent(in) :: only_nav
      
      allocate ( out % chn ( NUM_CHN))
      
      out % success = .true.
      call set_filenames ( config )
     
      call ahi_time_from_filename ( trim ( config %file_base) , out % time_start_obj, out % time_end_obj )
     
      call read_navigation ( config , out )
     
      if ( config % do_res3 ) call read_navigation ( config , out , do_geo3_inp = .true.)
     
      if ( config % do_res124 ) call read_navigation ( config , out , do_geo124_inp = .true.)
      if ( .not. present ( only_nav )) then
         call read_ahi_level1b ( config , out )
      end if
   end subroutine get_ahi_data
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine ahi_time_from_filename ( file_base , time0,time1 )
      character ( len = * ) :: file_base 
      type ( date_type) :: time0, time1
      
      integer :: year
      integer :: month
      integer :: day
      integer :: hour
      integer :: minute
      
      integer, parameter :: TIME_PERIOD_AHI_MINUTE = 8
      
      read(file_base(8:11), fmt="(I4)") year
      read(file_base(12:13), fmt="(I2)") month
      read(file_base(14:15), fmt="(I2)") day
      read(file_base(17:18), fmt="(I2)") hour
      read(file_base(19:20), fmt="(I2)") minute

      call time0 % set_date ( &
            year , month, day , hour, minute )
      
      time1 = time0
      
      call time1 % add_time ( minute = TIME_PERIOD_AHI_MINUTE)
       
   end subroutine ahi_time_from_filename
 
 
   ! needed since HCAST native is not named the same - WCS3
   !IMG_DK01B16_201909110330
   subroutine ahi_hcast_time_from_filename ( file_base , time0,time1 )
      character ( len = * ) :: file_base 
      type ( date_type) :: time0, time1
      
      integer :: year
      integer :: month
      integer :: day
      integer :: hour
      integer :: minute
      
      !left this alone for now
      integer, parameter :: TIME_PERIOD_AHI_MINUTE = 8
   
      read(file_base(13:16), fmt="(I4)") year
      read(file_base(17:18), fmt="(I2)") month
      read(file_base(19:20), fmt="(I2)") day
      read(file_base(21:22), fmt="(I2)") hour
      read(file_base(23:24), fmt="(I2)") minute

      call time0 % set_date ( &
            year , month, day , hour, minute )
      
      time1 = time0
      
      call time1 % add_time ( minute = TIME_PERIOD_AHI_MINUTE)
       
   end subroutine ahi_hcast_time_from_filename
   
   ! -------------------------------------------------
   !    returns offset and count for lon / lat value
   !
   !     lon/lat box is defined in config structure
   !      INPUT
   !
   !    OUTPUT:
   !       offset is  2 element vector holding the start of array for each dimension
   !       count_1 is a 2 elemen vector holding the number of elements in array for each dimension 
   ! --------------------------------------------------
   subroutine ahi_segment_information_region ( config , offset, count_1 ) 
      implicit none
      type ( ahi_config_type ), intent(in) :: config
      integer, intent(out) :: offset(2)
      integer, intent(out) :: count_1(2)
      type ( ahi_config_type ) :: config_local
      
      integer :: ii
      logical, allocatable :: inside (:,:)
      integer, allocatable :: line_g(:), elem_g(:)
      integer, parameter :: N_ELEMENTS_FULL_DISK = 5500
      integer, parameter :: N_LINES_FULL_DISK = 5500
      integer, save :: counter = 0
      
      config_local = config
      config_local % h5_offset = [0,0]
      config_local % h5_count  = [5500,5500]
      config_local % chan_on = .false.
      config_local % do_lonlat_only = .true.
      
      call set_filenames ( config_local )
      
      if ( .not. out_fd % set ) then
        call read_navigation ( config_local , out_fd )
        out_fd % set = .true.
        counter = counter+1
      end if  
   
      allocate ( inside (5500,5500))
      allocate ( line_g(5500),elem_g(5500))
      
      
      !- if lon /lat edges are set to -180,180 and -90,90 keep the full disk.
      
      if ( config % lon_range(1) .eq. -180.0  &
         .and. config % lon_range(2) .eq. 180.0  &
         .and. config % lat_range(1) .eq. -90.0  &
         .and. config % lat_range(2) .eq. 90.0  ) then
         
         offset = [0,0]
         count_1 =[5500,5500]
      else    
      
         if ( config % lon_range(2) .ge. config % lon_range(1) ) then
      
         inside =  out_fd % geo % lon .gt. config % lon_range(1) .and. &
           out_fd % geo % lon .lt. config % lon_range(2) .and. &
           out_fd % geo % lat .gt. config % lat_range(1) .and. &
           out_fd % geo % lat .lt. config % lat_range(2)
           
      
         else 
        
         inside =  (( out_fd % geo % lon .lt. config % lon_range(2) .and. &
               & out_fd % geo % lon .ge. -180.0 )  .or. &
               & (out_fd % geo % lon .gt. config % lon_range(1) .and. &
               & out_fd % geo % lon .lt. 180.0 )) .and. &
               
               & out_fd % geo % lat .gt. config % lat_range(1) .and. &
               & out_fd % geo % lat .lt. config % lat_range(2)
      
         end if
      
      
         elem_g = count (inside ,2 )      
         line_g = count (inside ,1 ) 
      
     
         do ii =1 , 5500
            if ( elem_g(ii) .ne. 0 ) then
               offset(1) = ii
           
               exit
            end if
         end do
      
         do ii =1 , 5500
            if ( line_g(ii) .ne. 0 ) then
               offset(2) = ii
            
               exit
            end if
         end do
      
         do ii =5500 , 1, -1
            if ( elem_g(ii) .ne. 0 ) then
               count_1(1) = ii - offset(1)
           
               exit
            end if
         end do
      
         do ii =5500 , 1, -1
            if ( line_g(ii) .ne. 0 ) then
               count_1(2) = ii - offset(2)
            
               exit
            end if
         end do
      end if
      
   end subroutine ahi_segment_information_region
   
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine set_filenames ( config)
      use cx_string_tools_mod, only: replace_text
      
      type ( ahi_config_type ) :: config
      
      integer :: i
      character(len=2) :: identifier
      character ( len=1020) :: file_for_this_channel
      integer,pointer::dims3(:)=>null()
      integer,pointer::dims7(:)=>null()
      integer,pointer::dims124(:)=>null()
      
      do i = 1 , 16
         
         write (identifier , fmt ='(i2.2)') i
        
         file_for_this_channel = replace_text ( config % file_base,'B01','B'//identifier)
        
         config % filename ( i ) = trim (config % data_path)//trim(file_for_this_channel) 
         config % varname ( i ) = '/RAD'
                  
      end do
      
      ! - compare channel 7 to 3 and 124 to find out if they are different 
     
      config % filenames_set = .true.
      
      call get_var_dimension ( config, 7, dims7)
      call get_var_dimension ( config, 3, dims3)
      call get_var_dimension ( config, 1, dims124)
      
      
      
      if ( dims3(1) .ne. dims7(1) ) then
         ! 4 times bigger resolution for channel 3
         config % do_res3 = .true.
         config % h5_offset_3 = 4 * config % h5_offset
         config % h5_count_3 = 4 * config % h5_count
      end if
      
      
      if ( dims124(1) .ne. dims7(1) ) then
         ! 2 times bigger resolution for channel 3
         config % do_res124 = .true.
         config % h5_offset_124 = 2 * config % h5_offset
         config % h5_count_124 = 2 * config % h5_count
      end if
      

   end subroutine set_filenames
   
   !--------------------------------------------------------------------------------------
   !
   !--------------------------------------------------------------------------------------
   
   subroutine get_var_dimension ( config, chn, dims )
      use cx_readh5dataset, only: &
         H5_DATASET_DIMENSIONS
      type ( ahi_config_type ) :: config
      integer, intent(in) :: chn
      integer :: dclass
     
      integer,pointer::dims(:)
     
       if (.not. config % filenames_set) call set_filenames(config)
      
      call H5_DATASET_DIMENSIONS (trim(config % filename ( chn )),trim(config % varname ( chn )),dims,dclass)
       
      
   end subroutine get_var_dimension

   ! --------------------------------------------------------------------------------------
   !  This reads navigation properties from AHI file
   !    and popluates ahi % geo substructure
   !  Variables are
   !     lon
   !     lat   
   !    solzen
   !    solaz
   !    satzen ( satellite Zenith)
   !    sataz   ( satellite azimuth )
   !    relaz   ( relative azimuth difference )
   !    glintzen
   !    satangle    
   ! --------------------------------------------------------------------------------------
   subroutine read_navigation ( config, ahi , do_geo3_inp, do_geo124_inp)  
      
      use viewing_geometry_mod, only: &
            possol &
         , sensor_zenith &
         , sensor_azimuth &
         , relative_azimuth &
         , glint_angle &
         , scattering_angle
         
      use cx_readh5dataset, only: &
         h5readattribute 
      
!using one CGMS routine (geos_transform_pix.c for all sensors)
!      use geo_sat_navigation_mod, only: &
!         fgf_to_earth
      
      implicit none
      
      type ( ahi_config_type ) :: config
      type ( ahi_data_out_type ) :: ahi
      logical, optional :: do_geo3_inp
      logical, optional :: do_geo124_inp
      
      character (len=120) ::  attr_name
      
      real (8) :: cfac
      real (8) :: coff
      real (8) :: lfac
      real (8) :: loff
      real (8) :: sub_lon 
      real (8) :: sub_lat  
      real (8) :: latx, lonx
      integer  :: ii,jj
      integer  :: x_full_disk
      integer  :: y_full_disk
   
      real :: GEO_ALTITUDE = 35786.0 !km
      
      ! - this is needed because it is not available in each of the files
      integer :: VALID_PROJECTION_CHANNEL 
      
      integer :: day_of_year
      real :: hour_frac
      
      type(geo_str) :: geo
      logical :: do_geo3
      logical :: do_geo124
      
      integer :: h5_offset(2)
      integer :: h5_count(2)
      integer :: FGF_TYPE = 3 !MTSAT uses JMA GEOS navigation, so set type here
      
      
       ! - navigation
     
      do_geo3 = .false.
      do_geo124 = .false.
      h5_offset = config % h5_offset
      h5_count = config % h5_count
      
      if (present(do_geo3_inp)) then
         do_geo3 = do_geo3_inp
         h5_offset = config % h5_offset_3
      h5_count = config % h5_count_3
      else if (present(do_geo124_inp)) then
         do_geo124 = do_geo124_inp
         h5_offset = config % h5_offset_124
         h5_count = config % h5_count_124
      end if
      
      VALID_PROJECTION_CHANNEL = 7
      if ( do_geo3) VALID_PROJECTION_CHANNEL = 3 
      if ( do_geo124) VALID_PROJECTION_CHANNEL = 1
      
      
      ! - read in the constants fom file
      attr_name = trim('Projection/CFAC')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION_CHANNEL ) ) , trim ( attr_name ), CFAC )
      attr_name = trim('Projection/LFAC')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION_CHANNEL) ) , trim ( attr_name ), LFAC )
      attr_name = trim('Projection/COFF')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION_CHANNEL ) ) , trim ( attr_name ), COFF )
      attr_name = trim('Projection/LOFF')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION_CHANNEL ) ) , trim ( attr_name ), LOFF )
      attr_name = trim('Projection/latitude_of_projection_origin')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION_CHANNEL ) ) , trim ( attr_name ), sub_lat )
      attr_name = trim('Projection/longitude_of_projection_origin')
      call h5readattribute ( trim(config % filename ( VALID_PROJECTION_CHANNEL) ) , trim ( attr_name ), sub_lon )
      
      
      call  geo %  allocate_geo ( h5_count(1), h5_count(2)) 
     
      call ahi % time_start_obj % get_date ( doy = day_of_year, hour_frac = hour_frac )
      
      do jj = 1 ,  h5_count(2)     
         do ii = 1 ,    h5_count(1)
            
            x_full_disk = ii +  h5_offset(1)
            y_full_disk = jj +  h5_offset(2)
            
            
!            call fgf_to_earth (  dble(x_full_disk), dble(y_full_disk) , cfac, coff, lfac, loff, sub_lon &
!               , lonx , latx )
               

!using one CGMS routine (geos_transform_pix.c for all sensors)
            call fgf_to_earth (FGF_TYPE,  dble(x_full_disk), dble(y_full_disk) , cfac, coff, lfac, loff, &
                                sub_lon, lonx , latx )

   
           
            geo % lat (ii,jj) = latx  
            geo % lon (ii,jj) = lonx  
                        
            if ( lonx == -999. ) cycle
            if ( config % do_lonlat_only ) cycle
             
            call  possol ( day_of_year ,  hour_frac  , real(lonx) &
                  , real(latx), geo % solzen (ii,jj), geo % solaz (ii,jj) )
                        
             geo % satzen (ii,jj) = sensor_zenith ( GEO_ALTITUDE, real(sub_lon),real(sub_lat) &
               ,real(lonx) , real(latx) )
             
              geo % sataz (ii,jj) = sensor_azimuth (  real(sub_lon),real(sub_lat),real(lonx) , real(latx) )
            
              geo % relaz(ii,jj) = relative_azimuth ( geo % solaz (ii,jj) , geo % sataz (ii,jj))
           
              geo % glintzen(ii,jj) = glint_angle ( geo % solzen (ii,jj) ,geo % satzen (ii,jj) &
                  ,   geo % relaz(ii,jj) )
           
              geo % scatangle(ii,jj) = scattering_angle (  geo % solzen(ii,jj)  &
                  ,  geo % satzen(ii,jj) ,  geo % relaz(ii,jj)) 
  
         end do
      end do  
      
      
    
      where (  geo % lat == -999.0)
          geo % is_space = .false.
      end where
      
      if ( do_geo3) then
         ahi % geo_3 = geo
         ahi % do_res3 = .true.
      else if ( do_geo124) then
         ahi % geo_124 = geo
         ahi % do_res124 = .true.
      else
         ahi % geo = geo
         
      end if
      
      
      call geo % deallocate_geo
   end subroutine read_navigation
   
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine read_ahi_level1b ( config, ahi )
   
      use cx_readh5dataset !, only: &
        ! h5readattribute  &
        ! , h5readdataset
      
      use file_utils, only: &
         file_test 
      
      implicit none
      
      type ( ahi_config_type ) :: config
      type ( ahi_data_out_type ) :: ahi
     
      !integer(kind = 2), pointer :: i2d_buffer( : , : ) => null()
      integer(kind = 2), pointer :: i2d_buffer(:,:)
      integer:: i_chn
      character (len=120) :: attr_name
      real (8) :: scale_factor, add_offset
      integer ( kind = 2 ) :: fillvalue
      integer ( kind = 4 ) , allocatable :: buffer_fake_i4 (:,:)
      real ( 8 ) :: cprime ! name of the radiance to reflectance coefficient in AHI file
      integer :: h5_offset(2)
      integer :: h5_count(2)
            
      ! - executable
      
      
      ahi % do_res3 =  config % do_res3
      ahi % do_res124 =  config % do_res124
      
      ! - channel data read 
      
      do i_chn = 1 ,16
         if ( .not. config % chan_on ( i_chn ) ) cycle
        
         h5_offset = config % h5_offset
         h5_count = config % h5_count
         if (i_chn .eq. 3 .and. config % do_res3 ) then
            h5_offset = config % h5_offset_3
            h5_count = config % h5_count_3
         end if
       
         if ((i_chn .lt. 3 .or. i_chn .eq. 4) .and. config % do_res124 ) then
            h5_offset = config % h5_offset_124
            h5_count = config % h5_count_124
         end if
           
           
         if ( .not. file_test ( trim(config % filename ( i_chn ) ) ) ) then 
            print*, 'AHI READER ERROR>> file '// trim(config % filename ( i_chn )) // ' not existing !!'
            ahi % success = .false.
            return
         end if
         ! - Read the data into buffer
         
       
         call h5readdataset ( trim(config % filename ( i_chn ) ) , trim ( config % varname(i_chn) ) &
               ,  h5_offset, h5_count, i2d_buffer )
         
	 if ( .not. associated(i2d_buffer)) then
	   print*,'file ',trim(config % filename ( i_chn ) )  , ' not readable'
	   cycle
	 endif 
         if ( .not. allocated (buffer_fake_i4) ) allocate ( buffer_fake_i4 ( h5_count(1), h5_count(2)))
        
         ! - fortran does not support unsigned integer
         
         buffer_fake_i4 = i2d_buffer
         
         deallocate ( i2d_buffer )
         where ( buffer_fake_i4 < 0 )
            buffer_fake_i4 = buffer_fake_i4 + 65536
         end where
         
         !- variable attributes
         attr_name = trim(config % varname(i_chn))//'/scale_factor'
         call h5readattribute ( trim(config % filename ( i_chn ) ) , trim ( attr_name ), scale_factor )
         attr_name = trim(config % varname(i_chn))//'/add_offset'
         call h5readattribute ( trim(config % filename ( i_chn ) ) , trim ( attr_name ), add_offset )
         attr_name = trim(config % varname(i_chn))//'/_FillValue'
         call h5readattribute ( trim(config % filename ( i_chn ) ) , trim ( attr_name ), fillvalue )
         if ( fillvalue < 0 ) fillvalue = fillvalue + 65536
        
         allocate ( ahi % chn(i_chn) % rad ( h5_count(1), h5_count(2)))
         
         ahi % chn(i_chn) % rad = (buffer_fake_i4 * scale_factor) + add_offset
         
        
         where ( buffer_fake_i4 == fillvalue )
            ahi % chn(i_chn) % rad = -999.
         end where
         
         if (allocated ( buffer_fake_i4 ) )  deallocate ( buffer_fake_i4 )
        
         ahi % chn(i_chn) % is_solar_channel = .false.
         if ( i_chn < 7 ) ahi % chn(i_chn) % is_solar_channel = .true.
         
         if ( ahi % chn(i_chn) % is_solar_channel ) then
            attr_name = trim(config % varname(i_chn))//'/cprime'
            call h5readattribute ( trim(config % filename ( i_chn ) ) , trim ( attr_name ), cprime )
            allocate ( ahi % chn(i_chn) % ref ( h5_count(1), h5_count(2)))
            ahi % chn(i_chn) % ref =100. *  ahi % chn(i_chn) % rad * cprime
         end if
       
         ahi % chn(i_chn) % is_read = .true.
        
      end do
   
   end subroutine read_ahi_level1b
   

   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine allocate_geo ( this, nx , ny )
      class ( geo_str ) :: this
      integer, intent(in) :: nx 
      integer, intent(in) :: ny
      
      call this % deallocate_geo 
          
      allocate (  this  % lon        (nx , ny) )
      allocate (  this  % lat        (nx , ny) )
      allocate (  this  % solzen     (nx , ny) )
      allocate (  this  % solaz      (nx , ny) )
      allocate (  this  % satzen     (nx , ny) )
      allocate (  this  % sataz      (nx , ny) )
      allocate (  this  % relaz      (nx , ny) )
      allocate (  this  % glintzen   (nx , ny) )
      allocate (  this  % scatangle  (nx , ny) )
      allocate (  this  % is_space   (nx , ny) )
      
      this  % lon        =   -999.
      this  % lat        =   -999.
      this  % solzen     =   -999.
      this  % solaz      =   -999.
      this  % satzen     =   -999.
      this  % sataz      =   -999.
      this  % relaz      =   -999.
      this  % glintzen   =   -999.
      this  % scatangle  =   -999.
      this  % is_space = .true.
      
   
   end subroutine allocate_geo
   
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine deallocate_geo (this )
      class ( geo_str ) :: this
      if (allocated ( this % lon)) deallocate ( this % lon) 
      if (allocated ( this % lat)) deallocate ( this % lat) 
      
           
      if ( allocated  (  this % solzen   ) ) deallocate (  this  % solzen   )
      if ( allocated  (  this % solaz      ) ) deallocate (  this % solaz   )
      if ( allocated  (  this % satzen     ) ) deallocate (  this  % satzen  )
      if ( allocated  (  this % sataz     ) ) deallocate (  this  % sataz  )
      if ( allocated  (  this % relaz      ) ) deallocate (  this  % relaz   )
      if ( allocated  (  this % glintzen    ) ) deallocate (  this  % glintzen   )
      if ( allocated  (  this % scatangle   ) ) deallocate (  this  % scatangle   )
      if ( allocated  (  this % is_space    ) ) deallocate (  this  % is_space  )
      if ( allocated  (  this % scan_time    ) ) deallocate (  this  % scan_time )
   
   
   end subroutine deallocate_geo 
   
   
   ! --------------------------------------------------------------------------------------
   !
   ! --------------------------------------------------------------------------------------
   subroutine deallocate_all (this )
      class ( ahi_data_out_type ) :: this
      integer :: i_chn
       
       call this % geo % deallocate_geo
      
      do i_chn = 1 , size (  this  % chn )
         
         if (allocated (  this  % chn(i_chn) % rad )) deallocate (this  % chn(i_chn) % rad)
         if (allocated (  this  % chn(i_chn) % bt ) ) deallocate (this  % chn(i_chn) % bt)
         if (allocated (  this  % chn(i_chn) % ref )) deallocate (this  % chn(i_chn) % ref)
      end do
      
      if ( allocated ( this % chn)) deallocate (this % chn)
   
   
   end subroutine deallocate_all
   
   
   subroutine clean_config ( this)
      class ( ahi_config_type ) :: this
      
      this %  filenames_set = .false.
      this % do_res3 = .false.
      this % do_res124 = .false.
   end subroutine clean_config 

end module cx_read_ahi_mod
