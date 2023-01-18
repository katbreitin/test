module fci_mod

  use class_time_date, only: date_type
  use cx_sds_io_mod, only:cx_sds_read,cx_sds_finfo
  use file_utils, only: file_search
  use viewing_geometry_mod, only: &
        possol &
     , sensor_zenith &
     , sensor_azimuth &
     , relative_azimuth &
     , glint_angle &
     , scattering_angle
  !use cx_real_boolean_mod

  use cx_geo__define, only: geo_str

  implicit none

  character(len=6) :: chn_string(16) =  &
     [ 'vis_04', &
       'vis_05', &
       'vis_06', &
       'vis_08', &
       'vis_09', &
       'nir_13', &
       'nir_16', &
       'nir_22', &
       'ir_38 ', &
       'wv_63 ', &
       'wv_73 ', &
       'ir_87 ', &
       'ir_97 ', &
       'ir_105', &
       'ir_123', &
       'ir_133' ]

  type fci_config
    logical :: init = .false.
    logical :: chan(16) = .false.
    character(len=1024) :: file_identifier
    character(len=1024) :: path
    character(len=1020),pointer :: file_list(:)
  contains
    procedure :: set=>fci_config__set
  end type fci_config

  type obs_chn_data
    real, allocatable :: rad(:,:)
    real, allocatable :: rfl(:,:)
    real, allocatable :: bt(:,:)
  contains
    procedure :: alloc => obs_chn_data__alloc
  end type obs_chn_data

  type fci_data
    type(date_type) :: time
    type(fci_config) :: config
    type(obs_chn_data) :: ch(16)
    type(geo_str) :: geo
      type(geo_str) :: geo1km
    real, allocatable :: sol_azi(:)
    real, allocatable :: solaz(:,:)
    real, allocatable :: solzen(:,:)
    real :: sat_sub_lon
    real :: sat_sub_lat
    real :: sat_altitude_km
    real :: earth_sun_distance
    real, allocatable :: lon(:,:)
    real, allocatable :: lat(:,:)
    real, allocatable :: lon_1km(:,:)
    real, allocatable :: lat_1km(:,:)
    character(len=1024):: path
    character(len=1024) :: file
   contains
    procedure :: get => fci_data__get
  end type fci_data

contains
  ! ------------
  !
  ! -------------
  subroutine  obs_chn_data__alloc(self,nx,ny)
    class(obs_chn_data) :: self
    integer :: nx, ny

    allocate (self%rad(nx,ny))
    allocate (self%bt(nx,ny))
    allocate (self%rfl(nx,ny))

  end subroutine  obs_chn_data__alloc

  subroutine fci_config__set(self, path,ch_on)
    class(fci_config) :: self
    character(len=*), intent(in)  :: path
    logical, intent(in)  :: ch_on(16)
    character(len=1020) :: file_identifier

    ! some tests if files are there
    self % file_list => file_search(trim(path),'*.nc')

    ! then this is done, put the things in object
    self % path = path
    self % chan(:) = ch_on
    self % init = .true.


  end subroutine fci_config__set


  subroutine fci_data__get (self,chunk, start, count, time_only)
     class(fci_data) :: self
     integer :: i, ii, jj
     integer, intent(in) :: chunk  ! the chunk out of 40 for this granule
     integer, intent(in), optional :: start(2)
     integer, intent(in), optional :: count(2)
     logical, optional, intent(in) :: time_only
     character(len =1024) :: file_chunk
     real,allocatable::rad_bt_a(:), rad_bt_b(:), rad_bt_v(:)
     real,allocatable::rad_bt_c1(:), rad_bt_c2(:),rad_bt_conv(:)
     real, allocatable :: irrad(:)
     real, allocatable :: dum_1d(:)
     integer :: status
     real,allocatable :: time(:)
     type(date_type) :: tt
     character ( len = 128), allocatable :: Sds_Name(:)
     character ( len = 128), allocatable :: Att_Name(:)
     integer :: nsds , ftype,natt
     character (len =1024) :: lonlat_file1km
     character (len =1024) :: lonlat_file2km
     logical :: t_only = .false.
     integer :: stride(2)

     integer :: day_of_year
     real :: hour_frac
     real, parameter :: Pi = 3.14159265359
      real , parameter :: DTOR = PI / 180.

     ! ---------------
     file_chunk = trim(self%config%path)//trim(self %config %file_list(chunk))
     status = cx_sds_finfo (File_chunk, ftype,nsds,Sds_Name,Natt,Att_Name)
print*,trim(file_chunk)
     ! time from 2000-01-01 00:00 in seconds
     status=  cx_sds_read (trim(file_chunk), &
       '/time' , time )

      call self % time % set_date(year=2000,month=1 &
                   ,day=1,hour=0,minute=0,second=0)
      call self % time % add_time(second =int(time(1)))

       t_only = .false.
      if ( present(time_only)) t_only = time_only
      if (t_only) THEN
        return
      end if

!TODODODOD
! channel ir_38 is different
!  has a warm scale offset and warn scale slope



     lonlat_file1km = '/Users/awalther/DATA/Satellite_Input/FCI/CM_OPE_GRIDDEF_MTI1+FCI_20220407120000_1km-V2.nc'
     lonlat_file2km = '/Users/awalther/DATA/Satellite_Input/FCI/CM_OPE_GRIDDEF_MTI1+FCI_20220407120000_2km-V2.nc'

    status = cx_sds_read (lonlat_file2km, &
          'longitude' , self % lon , start = [1,139 * chunk +1], count = [5568,139] )

    status = cx_sds_read (lonlat_file2km, &
               'latitude' , self % lat , start = [1,139 * chunk+1 ], count = [5568,139])

               status = cx_sds_read (lonlat_file1km, &
                     'longitude' , self % lon_1km , start = [1,279 * chunk +1], count = [11136,278] )

               status = cx_sds_read (lonlat_file1km, &
                          'latitude' , self % lat_1km , start = [1,279 * chunk+1 ], count = [11136,278])
               ! geometry
                status = cx_sds_read (trim(file_chunk), &
                '/state/platform/subsatellite_latitude' , dum_1d )
                self % sat_sub_lat = dum_1d(1)
                status = cx_sds_read (trim(file_chunk), &
                     '/state/platform/subsatellite_longitude' , dum_1d )
               self % sat_sub_lon = dum_1d(1)
               status = cx_sds_read (trim(file_chunk), &
                    '/state/platform/platform_altitude' , dum_1d )
                self % sat_altitude_km = dum_1d(1)
                status = cx_sds_read (trim(file_chunk), &
                     '/state/celestial/earth_sun_distance' , dum_1d )
                self % earth_sun_distance = dum_1d(1)

!day_of_Year=230
!hour_frac = 12.3
      call self % geo % set ( self % lon, self % lat,self % time &
         ,self % sat_sub_lon,self % sat_sub_lat,self % sat_altitude_km )
print*,shape(self % lon_1km)
         call self % geo1km % set ( self % lon_1km, self % lat_1km,self % time &
            ,self % sat_sub_lon,self % sat_sub_lat,self % sat_altitude_km )

    do i=1,16

      if ( self % config % chan(i) ) then

      !  if (  trim(chn_string(i)) .ne. 'ir_38') cycle
              status=  cx_sds_read (trim(file_chunk), &
                '/data/'//trim(chn_string(i))//'/measured/effective_radiance' &
                , self%ch(i)%rad ,start=start, count = count, stride = stride)


      !  else
          !  status=  cx_sds_read_fci_ir38 (trim(file_chunk), &
          !    '/data/'//trim(chn_string(i))//'/measured/effective_radiance' &
          !    , self%ch(i)%rad ,start=start, count = count,)



    !    end if
            ! var_names


          if ( i .lt. 9) THEN
            status =    cx_sds_read (trim(file_chunk), &
              '/data/'//trim(chn_string(i))// &
               '/measured/channel_effective_solar_irradiance' &
              , irrad )

! TODO   T-CHECK
!  read page 50 of https://www.eumetsat.int/media/45923

!      rfl = PI * Rad * d2  /  (irrad * cos (solar_zenith) )

          !    self%ch(i)%rfl = self%ch(i)%rad/irrad(1)

          print*,shape(self % geo1km % solzen),shape(self%ch(i)%rad)
self%ch(i)%rfl = (PI * self%ch(i)%rad(1:11136,1:278) * self % earth_sun_distance**2) &
            /  ( irrad(1) * cos(self % geo1km % solzen * DTOR))

          end if

           if ( i .gt. 8 ) then
              status =    cx_sds_read (trim(file_chunk), &
                '/data/'//trim(chn_string(i))// &
                 '/measured/radiance_to_bt_conversion_coefficient_a' &
                , rad_bt_a)

               status =    cx_sds_read (trim(file_chunk), &
                        '/data/'//trim(chn_string(i))// &
                        '/measured/radiance_to_bt_conversion_coefficient_b' &
                        , rad_bt_b)

               status =    cx_sds_read (trim(file_chunk), &
                          '/data/'//trim(chn_string(i))// &
                          '/measured/radiance_to_bt_conversion_coefficient_wavenumber' &
                        , rad_bt_v)

              status =   cx_sds_read (trim(file_chunk), &
                            '/data/'//trim(chn_string(i))// &
                            '/measured/radiance_to_bt_conversion_constant_c1' &
                       , rad_bt_c1)

              status =    cx_sds_read (trim(file_chunk), &
                        '/data/'//trim(chn_string(i))// &
                        '/measured/radiance_to_bt_conversion_constant_c2' &
                     , rad_bt_c2)

              status =    cx_sds_read (trim(file_chunk), &
                            '/data/'//trim(chn_string(i))// &
                            '/measured/radiance_unit_conversion_coefficient' &
                      , rad_bt_conv)

     ! https://www-cdn.eumetsat.int/files/2020-04/pdf_effect_rad_to_brightness.pdf
     ! Eq.5.3
     ! EUMETSAT Doc.No. : EUM/MET/TEN/11/0569
        self%ch(i)%bt = ((rad_bt_c2(1) * rad_bt_v(1) ) / &
          (rad_bt_a(1) * alog(rad_bt_c1(1) * (rad_bt_v(1)**3)/self%ch(i)%rad + 1))) &
              - rad_bt_b(1)/rad_bt_a(1)
          where ( self%ch(i)%rad .lt. 0. )
               self%ch(i)%bt = -999.
          end where


            end if

        end if

     end do






  end   subroutine fci_data__get


end module fci_mod
