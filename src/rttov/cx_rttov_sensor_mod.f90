module cx_rttov_sensor_mod

implicit none

type cx_rttov_sensor_type
    integer :: WMO_id
    character(50) :: name_platform
    character(50) :: name_sensor
    character(50) :: name_clavrx
    character(50) :: name_rttov
    character(50) :: name_dcomp
    character(200) :: coef_filename
    character(200) :: cld_coef_filename
    character (200) :: path_rttov
    real :: max_satzen = 85.
    integer :: chan_src_on_cx(45) = -1

  contains
    procedure, public :: init => cx_rttov_sensor_init
    procedure, private :: update => cx_rttov_sensor_update

end type cx_rttov_sensor_type
contains

subroutine cx_rttov_sensor_init (self, wmo_id, path, sensor )
  class (cx_rttov_sensor_type), intent(inout) :: self
  integer , optional, intent(in) :: wmo_id
  character(len = *), intent(in)  :: path
  character(len = *), intent(in), optional :: sensor
  character(len = 20) :: sensor_loc
  ! -  sensor is only needed for hirs  ( same wmo_id for different sensors)
  sensor_loc = 'default'
  if (present(sensor)) sensor_loc = sensor
  if (present(wmo_id)) self%wmo_id = wmo_id
  self % path_rttov = path

   call self % update(sensor_loc)
end subroutine

  subroutine cx_rttov_sensor_update (self, sensor)
    class (cx_rttov_sensor_type), intent(inout) :: self
    character(len=*), intent(in) :: sensor
    character (len=1) :: rttov_version_string
    integer :: i

    rttov_version_string = '9'

    select case(self % WMO_Id)
    case(4) !METOP-A

      if (index(sensor,'HIRS') .gt. 0 ) THEN
        self % name_dcomp = 'HIRS-METOPA'
        self % name_rttov = 'metop_2_hirs-shifted'
        self % chan_src_on_cx(20:36) = &
        [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
        rttov_version_string = '8'
      else
        self % name_dcomp = 'AVHRR-METOPA'
        self % name_rttov = 'metop_2_avhrr'
        self % chan_src_on_cx(6) = 3
        self % chan_src_on_cx(20) = 4
        self % chan_src_on_cx(31) = 5
        self % chan_src_on_cx(32) = 6
      end if

    case(3) !METOP-B


      if (index(sensor,'HIRS') .gt. 0  ) THEN
        self % name_dcomp = 'HIRS-METOPB'
        self % name_rttov = 'metop_1_hirs-shifted'
        self % chan_src_on_cx(20:36) = &
        [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
        rttov_version_string = '8'
      else
        self % name_dcomp = 'AVHRR-METOPB'
        self % name_rttov = 'metop_1_avhrr'
        self % chan_src_on_cx(6) = 3
        self % chan_src_on_cx(20) = 4
        self % chan_src_on_cx(31) = 5
        self % chan_src_on_cx(32) = 6

      end if

    case(5) !METOP-C


     if (index(sensor,'HIRS')  .gt. 0 ) THEN
       self % name_dcomp = 'HIRS-METOPC'
       self % name_rttov = 'metop_3_hirs-shifted'
       self % chan_src_on_cx(20:36) = &
       [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
       rttov_version_string = '8'
    else
      self % name_dcomp = 'AVHRR-METOPC'
      self % name_rttov = 'metop_3_avhrr'
      self % chan_src_on_cx(6) = 3
      self % chan_src_on_cx(20) = 4
      self % chan_src_on_cx(31) = 5
      self % chan_src_on_cx(32) = 6

     end if


    case(55) !MSG-8
      self % name_dcomp = 'SEVIRI-MSG08'
      self % name_rttov = 'msg_1_seviri'
      self % chan_src_on_cx(20:37) &
       = [4,-1,-1,-1,-1,-1,-1,-1,6,7,8,9,10,11,-1,-1,-1,5]
    case(56) !MSG-9
      self % name_dcomp = 'SEVIRI-MSG09'
      self % name_rttov =  'msg_2_seviri'
      self % chan_src_on_cx(20:37) &
       = [4,-1,-1,-1,-1,-1,-1,-1,6,7,8,9,10,11,-1,-1,-1,5]
    case(57) !MSG-10
      self % name_dcomp = 'SEVIRI-MSG10'
      self % name_rttov =  'msg_3_seviri'
      self % chan_src_on_cx(20:37) &
       = [4,-1,-1,-1,-1,-1,-1,-1,6,7,8,9,10,11,-1,-1,-1,5]
    case(70) !MSG-11
      self % name_dcomp = 'SEVIRI-MSG11'
      self % name_rttov = 'msg_4_seviri'
      self % chan_src_on_cx(20:37) &
       = [4,-1,-1,-1,-1,-1,-1,-1,6,7,8,9,10,11,-1,-1,-1,5]
    case(171) !MTSAT-1R
      self % name_dcomp = 'MTSAT-1'
      self % name_rttov = 'mtsat_1_imager'

    case(172) !MTSAT-2
      self % name_dcomp = 'MTSAT-2'
      self % name_rttov = 'mtsat_2_imager'
    case(173) !AHI-8
      self % name_dcomp = 'AHI8'
      self % name_rttov = 'himawari_8_ahi'
      self % chan_src_on_cx(20:38) &
        = [7,-1,-1,-1,-1,-1,-1,9,10,11,12,14,15,16,-1,-1,-1,8,13]
    case(174) !AHI-9
      self % name_dcomp = 'AHI9'
      self % name_rttov = 'himawari_9_ahi'
      self % chan_src_on_cx(20:38) &
        = [7,-1,-1,-1,-1,-1,-1,9,10,11,12,14,15,16,-1,-1,-1,8,13]
    case(200) !NOAA-8


     if (index(sensor,'HIRS') .gt. 0  ) THEN
       self % name_dcomp = 'HIRS-NOAA08'
       self % name_rttov = 'noaa_8_hirs-shifted'
       self % chan_src_on_cx(20:36) = &
       [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
       rttov_version_string = '8'
    else
      self % name_dcomp = 'AVHRR-NOAA08'
      self % name_rttov = 'noaa_8_avhrr'
      rttov_version_string = '8'
      self % chan_src_on_cx(20) = 3
      self % chan_src_on_cx(31) = 4
      self % chan_src_on_cx(32) = 5

     end if

    case(201) !NOAA-9


     if (index(sensor,'HIRS')  .gt. 0 ) THEN
       self % name_dcomp = 'HIRS-NOAA09'
       self % name_rttov = 'noaa_9_hirs-shifted'
       self % chan_src_on_cx(20:36) = &
       [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
       rttov_version_string = '8'
    else
      self % name_dcomp = 'AVHRR-NOAA09'
      self % name_rttov = 'noaa_9_avhrr'
      rttov_version_string = '8'
      self % chan_src_on_cx(20) = 3
      self % chan_src_on_cx(31) = 4
      self % chan_src_on_cx(32) = 5

     end if

    case(202) !NOAA-10


     if (index(sensor,'HIRS')  .gt. 0 ) THEN
       self % name_dcomp = 'HIRS-NOAA10'
       self % name_rttov = 'noaa_10_hirs-shifted'
       self % chan_src_on_cx(20:36) = &
       [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
       rttov_version_string = '8'
    else
      self % name_dcomp = 'AVHRR-NOAA10'
      self % name_rttov = 'noaa_10_avhrr'
      self % chan_src_on_cx(20) = 3
      self % chan_src_on_cx(31) = 4
      self % chan_src_on_cx(32) = 5

     end if

   case(203) !NOAA-11



     if (index(sensor,'HIRS')  .gt. 0 ) THEN
       self % name_dcomp = 'HIRS-NOAA11'
       self % name_rttov = 'noaa_11_hirs-shifted'
       self % chan_src_on_cx(20:36) = &
       [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
       rttov_version_string = '8'
     else

       self % name_dcomp = 'AVHRR-NOAA11'
       self % name_rttov = 'noaa_11_avhrr'
       self % chan_src_on_cx(20) = 3
       self % chan_src_on_cx(31) = 4
       self % chan_src_on_cx(32) = 5

     end if

   case(204) !NOAA-12

     if (index(sensor,'HIRS')  .gt. 0 ) THEN
            self % name_dcomp = 'HIRS-NOAA12'
            self % name_rttov = 'noaa_12_hirs-shifted'
            self % chan_src_on_cx(20:36) = &
            [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
            rttov_version_string = '8'
      else
        self % name_dcomp = 'AVHRR-NOAA12'
        self % name_rttov = 'noaa_12_avhrr'
        self % chan_src_on_cx(20) = 3
        self % chan_src_on_cx(31) = 4
        self % chan_src_on_cx(32) = 5


      end if
   case(205) !NOAA-14


          if (index(sensor,'HIRS') .gt. 0  ) THEN
            self % name_dcomp = 'HIRS-NOAA14'
            self % name_rttov = 'noaa_14_hirs-shifted'
            self % chan_src_on_cx(20:36) = &
            [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
            rttov_version_string = '8'
          else
            self % name_dcomp = 'AVHRR-NOAA14'
            self % name_rttov = 'noaa_14_avhrr'
            self % chan_src_on_cx(20) = 3
            self % chan_src_on_cx(31) = 4
            self % chan_src_on_cx(32) = 5

          end if

   case(206) !NOAA-15


          if (index(sensor,'HIRS')  .gt. 0 ) THEN
            self % name_dcomp = 'HIRS-NOAA15'
            self % name_rttov = 'noaa_15_hirs-shifted'
            self % chan_src_on_cx(20:36) = &
            [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
            rttov_version_string = '8'
          else
            self % name_dcomp = 'AVHRR-NOAA15'
            self % name_rttov = 'noaa_15_avhrr'
            self % chan_src_on_cx(6) = 3
            self % chan_src_on_cx(20) = 4
            self % chan_src_on_cx(31) = 5
            self % chan_src_on_cx(32) = 6

          end if
   case(207) !NOAA-16



          if (index(sensor,'HIRS')  .gt. 0 ) THEN
            self % name_dcomp = 'HIRS-NOAA16'
            self % name_rttov = 'noaa_16_hirs-shifted'
            self % chan_src_on_cx(20:36) = &
            [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
            rttov_version_string = '8'
          else

            self % name_dcomp = 'AVHRR-NOAA16'
            self % name_rttov = 'noaa_16_avhrr'
            self % chan_src_on_cx(6) = 3
            self % chan_src_on_cx(20) = 4
            self % chan_src_on_cx(31) = 5
            self % chan_src_on_cx(32) = 6

          end if
   case(208) !NOAA-17



          if (index(sensor,'HIRS') .gt. 0  ) THEN
            self % name_dcomp = 'HIRS-NOAA17'
            self % name_rttov = 'noaa_17_hirs-shifted'
            self % chan_src_on_cx(20:36) = &
            [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
            rttov_version_string = '8'
          else
            self % name_dcomp = 'AVHRR-NOAA17'
            self % name_rttov = 'noaa_17_avhrr'
            self % chan_src_on_cx(6) = 3
            self % chan_src_on_cx(20) = 4
            self % chan_src_on_cx(31) = 5
            self % chan_src_on_cx(32) = 6

          end if

   case(209) !NOAA-18



          if (index(sensor,'HIRS')  .gt. 0 ) THEN
            self % name_dcomp = 'HIRS-NOAA18'
            self % name_rttov = 'noaa_18_hirs-shifted'
            self % chan_src_on_cx(20:36) = &
            [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
            rttov_version_string = '8'
          else

            self % name_dcomp = 'AVHRR-NOAA18'
            self % name_rttov = 'noaa_18_avhrr'
            self % chan_src_on_cx(6) = 3
            self % chan_src_on_cx(20) = 4
            self % chan_src_on_cx(31) = 5
            self % chan_src_on_cx(32) = 6
          end if

   case(223) !NOAA-19

          if (index(sensor,'HIRS') .gt. 0  ) THEN
            self % name_dcomp = 'HIRS-NOAA19'
            self % name_rttov = 'noaa_19_hirs-shifted'
            self % chan_src_on_cx(20:36) = &
            [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
            rttov_version_string = '8'
          else
            self % name_dcomp = 'AVHRR-NOAA19'
            self % name_rttov = 'noaa_19_avhrr'
            self % chan_src_on_cx(6) = 3
            self % chan_src_on_cx(20) = 4
            self % chan_src_on_cx(31) = 5
            self % chan_src_on_cx(32) = 6

          end if


   case(224) !VIIRS - SNPP
     self % name_dcomp = 'VIIRS-SNPP'
    self % name_rttov   = 'jpss_0_viirs'
    self % chan_src_on_cx(1:32) = [6,9,3,4,10,14,15,1,2,-1,-1,-1,-1,-1,7,-1,-1,-1,-1,17,  &
                      -1,18,-1,-1,-1,12,-1,-1,19,-1,20,22]
    self % chan_src_on_cx(39) = 1   ! I1
    self % chan_src_on_cx(40) = 2   ! I2
    self % chan_src_on_cx(41) = 13  ! I3
    self % chan_src_on_cx(42) = 16  ! I4
    self % chan_src_on_cx(43) = 21  ! I5
   case(225)  !VIIRS NOAA-20
      self % name_dcomp = 'VIIRS-N20'
      self % name_rttov = 'noaa_20_viirs'
      self % chan_src_on_cx(1:32) = [6,9,3,4,10,14,15,1,2,-1,-1,-1,-1,-1,7,-1,-1,-1,-1,17,  &
                        -1,18,-1,-1,-1,12,-1,-1,19,-1,20,22]
      self % chan_src_on_cx(39) = 1   ! I1
      self % chan_src_on_cx(40) = 2   ! I2
      self % chan_src_on_cx(41) = 13  ! I3
      self % chan_src_on_cx(42) = 16  ! I4
      self % chan_src_on_cx(43) = 21  ! I5
   case(226)  !VIIRS NOAA-21
      self % name_dcomp = 'VIIRS-N21'
      self % name_rttov = 'noaa_21_viirs'
      self % chan_src_on_cx(1:32) = [6,9,3,4,10,14,15,1,2,-1,-1,-1,-1,-1,7,-1,-1,-1,-1,17,  &
                        -1,18,-1,-1,-1,12,-1,-1,19,-1,20,22]
      self % chan_src_on_cx(39) = 1   ! I1
      self % chan_src_on_cx(40) = 2   ! I2
      self % chan_src_on_cx(41) = 13  ! I3
      self % chan_src_on_cx(42) = 16  ! I4
      self % chan_src_on_cx(43) = 21  ! I5
   case(250) !GOES-6
     self % name_dcomp = 'GOES-6'
     self % name_rttov = 'goes_6_imager'
     rttov_version_string = '8'

   case(251) !GOES-7
     self % name_dcomp = 'GOES-7'
     self % name_rttov = 'goes_7_imager'
     rttov_version_string = '8'

   case(252) !GOES-8
     self % name_dcomp = 'GOES-8'
     self % name_rttov = 'goes_8_imager'
     rttov_version_string = '7'
     self % chan_src_on_cx(20) = 1
     self % chan_src_on_cx(27) = 2
     self % chan_src_on_cx(31) = 3
     self % chan_src_on_cx(32) = 4
     self % max_satzen = 75.
   case(253) !GOES-9
     self % name_dcomp = 'GOES-9'
     self % name_rttov = 'goes_9_imager'
     rttov_version_string = '7'
     self % chan_src_on_cx(20) = 1
     self % chan_src_on_cx(27) = 2
     self % chan_src_on_cx(31) = 3
     self % chan_src_on_cx(32) = 4
     self % max_satzen = 75.
   case(254) !GOES-10
     self % name_dcomp = 'GOES-10'
     self % name_rttov = 'goes_10_imager'
     rttov_version_string = '8'
     self % chan_src_on_cx(20) = 1
     self % chan_src_on_cx(27) = 2
     self % chan_src_on_cx(31) = 3
     self % chan_src_on_cx(32) = 4
     self % max_satzen = 75.
   case(255) !GOES-11
     self % name_dcomp = 'GOES-11'
     self % name_rttov = 'goes_11_imager'
     rttov_version_string = '8'
     self % chan_src_on_cx(20) = 1
     self % chan_src_on_cx(27) = 2
     self % chan_src_on_cx(31) = 3
     self % chan_src_on_cx(32) = 4
     self % max_satzen = 75.
   case(256) !GOES-12
     self % name_dcomp = 'GOES-12'
     self % name_rttov = 'goes_12_imager'
     self % chan_src_on_cx(20) = 1
     self % chan_src_on_cx(27) = 2
     self % chan_src_on_cx(31) = 3
     self % chan_src_on_cx(33) = 4
     rttov_version_string = '8'
     self % max_satzen = 75.

   case(257) !GOES-13
     self % name_dcomp = 'GOES-13'
     self % name_rttov = 'goes_13_imager'
     rttov_version_string = '8'
     self % chan_src_on_cx(20) = 1
     self % chan_src_on_cx(27) = 2
     self % chan_src_on_cx(31) = 3
     self % chan_src_on_cx(33) = 4
     self % max_satzen = 75.

   case(258) !GOES-14
     self % name_dcomp = 'GOES-14'
     self % name_rttov = 'goes_14_imager'
     rttov_version_string = '8'
     self % chan_src_on_cx(20) = 1
     self % chan_src_on_cx(27) = 2
     self % chan_src_on_cx(31) = 3
     self % chan_src_on_cx(33) = 4
     self % max_satzen = 75.

   case(259) !GOES-15
     self % name_dcomp = 'GOES-15'
     self % name_rttov = 'goes_15_imager'
     rttov_version_string = '8'
     self % chan_src_on_cx(20) = 1
     self % chan_src_on_cx(27) = 2
     self % chan_src_on_cx(31) = 3
     self % chan_src_on_cx(33) = 4
     self % max_satzen = 75.

   case(270) !GOES-16
     self % name_dcomp = 'GOES-16'
     self % name_rttov = 'goes_16_abi'
     self % chan_src_on_cx(1:45) = [2,3,1,-1,-1,5,6,-1,-1,-1,-1,-1,-1,-1,-1,  &
                       -1,-1,-1,-1,7,-1,7,-1,-1,-1,4,9,10,11,12,  &
                       14,15,16,-1,-1,-1,8,13,-1,-1,-1,-1,-1,-1,-1]
   case(271) !GOES-17
     self % name_dcomp = 'GOES-17'
     self % name_rttov = 'goes_17_abi'
     self % chan_src_on_cx(1:45) = [2,3,1,-1,-1,5,6,-1,-1,-1,-1,-1,-1,-1,-1,  &
                       -1,-1,-1,-1,7,-1,7,-1,-1,-1,4,9,10,11,12,  &
                       14,15,16,-1,-1,-1,8,13,-1,-1,-1,-1,-1,-1,-1]
   case(272) !GOES-18
     self % name_dcomp = 'GOES-18'
     self % name_rttov = 'goes_18_abi'
     self % chan_src_on_cx(1:45) = [2,3,1,-1,-1,5,6,-1,-1,-1,-1,-1,-1,-1,-1,  &
                       -1,-1,-1,-1,7,-1,7,-1,-1,-1,4,9,10,11,12,  &
                       14,15,16,-1,-1,-1,8,13,-1,-1,-1,-1,-1,-1,-1]
    case(289) !GOES-18
      self % name_dcomp = 'FCI'
      self % name_rttov = 'mtg_1_fci'
      self % chan_src_on_cx(1:45) = [3,-1,-1,-1,-1,7,8,-1, 1,-1 &
                                     ,2,-1,-1,-1,-1,4,-1,-1,5,9 &
                                     ,-1,-1,-1,-1,-1,6,-1,-1,11,12  &
                                     ,13,14,15,16,-1,-1,-1,10,-1,-1 &
                                     ,-1,-1,-1,-1,-1]

   case(706) !NOAA-6
     self % name_dcomp = 'AVHRR-NOAA06'
     self % name_rttov = 'noaa_6_avhrr'
   case(707) !NOAA-7
     self % name_dcomp = 'AVHRR-NOAA07'
     self % name_rttov = 'noaa_7_avhrr'
   case(708) !NOAA-5
     self % name_dcomp = 'AVHRR-TIROSN'
     self % name_rttov = 'noaa_5_avhrr'
     rttov_version_string = '8'
   case(783) !MODIS
       self % name_dcomp = 'MODIS-TERRA'
       self % name_rttov = 'eos_1_modis'
       do i =1,45
         self % chan_src_on_cx(i) = i
       end do

   case(784) !MODIS
      self % name_dcomp = 'MODIS-AQUA'
      self % name_rttov = 'eos_2_modis'
      do i =1,45
        self % chan_src_on_cx(i) = i
      end do

   case(510) !FY2A
      self % name_dcomp ='FY2-1'

   case(514) !FY2D
      self % name_dcomp ='FY2-2'

   case(515) !FY2E
      self % name_dcomp ='FY2-3'

   case(523) !FY3D
      self % name_dcomp = 'FY3-D'
      self % name_rttov = 'fy3_4_mersi2'

   case(530) ! FY4-A
      self % name_dcomp ='FY4-A'

   case(810) !COMS
      self % name_dcomp ='COMS-1'

   case (840)
      self % name_dcomp ='EPS-SG'
      self % name_rttov = 'metopsg_1_metimage'

   case default
      print*,'sensor for WMO number not found in RT Utils  ', self % WMO_id
      print*,'stopping ... Please fix this in cx_rttov_sensor_mod.F90'
      print*,' better tell andi.walther@ssec.wisc.edu'
      stop
   end select
   if (rttov_version_string .eq. '7') self % max_satzen = 75.
   if (rttov_version_string .eq. '8') self % max_satzen = 75.

   self % coef_filename = trim(self %path_rttov)//'/rtcoef_rttov/rttov' &
        //rttov_version_string//'pred54L/rtcoef_'//trim(self % name_rttov)//'.dat'

   self % cld_coef_filename = trim(self %path_rttov)//'cldaer_ir/sccldcoef_' &
        //trim(self % name_rttov)//'.dat'

end subroutine


end module cx_rttov_sensor_mod
