module fci_mod


  use class_time_date, only: date_type
  use cx_sds_io_mod, only:cx_sds_read
  use file_utils, only: file_search

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
    real, allocatable :: sol_azi(:)
    character(len=1024):: path
    character(len=1024) :: file
   contains
    procedure :: get => fci_data__get
  end type fci_data

contains

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
    print*,'FCI setting'

  end subroutine fci_config__set


  subroutine fci_data__get (self,chunk, start, count)
     class(fci_data) :: self
     integer :: i
     integer, intent(in) :: chunk  ! the chunk out of 40 for this granule
     integer, intent(in), optional :: start(2)
     integer, intent(in), optional :: count(2)

     character(len =1024) :: file_chunk
     real,allocatable::rad_bt_a(:), rad_bt_b(:), rad_bt_v(:)
     real,allocatable::rad_bt_c1(:), rad_bt_c2(:),rad_bt_conv(:)
     real, allocatable :: irrad(:)
     integer :: status
     real,allocatable :: time(:)
     type(date_type) :: tt


    file_chunk = trim(self%config%path)//trim(self %config %file_list(chunk))

  ! time from 2000-01-01 00:00 in seconds
    status=  cx_sds_read (trim(file_chunk), &
       '/time' &
       , time )

     call self % time % set_date(year=2000,month=1 &
                 ,day=1,hour=0,minute=0,second=int(time(200)))

     print*,'granule time: ',self % time %date_string('yyyy_doy.hhmm')


!TODODODOD
! channel ir_38 is different
!  has a warm scale offset and warn scale slope

    do i=1,16

      if ( self % config % chan(i) ) then
        if (  trim(chn_string(i)) .ne. 'ir_38') then
              status=  cx_sds_read (trim(file_chunk), &
                '/data/'//trim(chn_string(i))//'/measured/effective_radiance/_DATA' &
                , self%ch(i)%rad ,start=start, count = count)
        else
            status=  cx_sds_read (trim(file_chunk), &
              '/data/'//trim(chn_string(i))//'/measured/effective_radiance' &
              , self%ch(i)%rad ,start=start, count = count)

              status=  cx_sds_read (trim(file_chunk), &
                '/data/'//trim(chn_string(i))//'/measured/effective_radiance' &
                , warm_slop,start=start, count = count)

             self%ch(i)%rad = 0.
        end if
            ! var_names


          if ( i .lt. 10) THEN
            status =    cx_sds_read (trim(file_chunk), &
              '/data/'//trim(chn_string(i))// &
               '/measured/channel_effective_solar_irradiance' &
              , irrad )

              self%ch(i)%rfl = self%ch(i)%rad/irrad(1)

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



  print*,i,trim(chn_string(i)),self%ch(i)%bt(500,100)
            end if






        end if

     end do

  end   subroutine fci_data__get


end module fci_mod
