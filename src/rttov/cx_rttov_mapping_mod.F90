module cx_rttov_mapping_mod

contains
  !
  !
  !
  function channel_map (sensor,ancil_data_path, chn,coef_filename,cld_coef_filename,max_satzen) result (list)

#define STRINGIFY(x) x

    implicit none
    character(len = *) :: sensor
    character(len = * ) :: ancil_data_path  
    integer :: chn
    integer :: list
    integer, allocatable:: chn_list(:)
   ! integer :: mod_list(5) = [20,22,29,31,32]
    integer :: kk
    character (len = *) :: coef_filename
    character (len = *) :: cld_coef_filename
    character (len=:), allocatable :: path
    character (len = 100) :: sensor_string
    character (len=1) :: rttov_version_string
    character (len = 1) :: metop_nr
    integer :: avhrr_num
    integer :: nr_hirs
    real :: max_satzen
   
    path = trim(ancil_data_path) // "static/rttov/"
    rttov_version_string = '9'
    max_satzen = 85.
    ! -- 
    !  the mapping translates one channel infot
    !  RTTOV channel number
    !
    !
   
    
    allocate(chn_list(45))
    chn_list = -1
    list = -1
    
    
    ! metop 
   
    
    if (sensor .eq. 'AVHRR-METOPA') metop_nr = '2'
    if (sensor .eq. 'AVHRR-METOPB') metop_nr = '1'
    if (sensor .eq. 'AVHRR-METOPC') metop_nr = '3'
   
    if (sensor .eq. 'HIRS-METOPA') metop_nr = '2'
    if (sensor .eq. 'HIRS-METOPB') metop_nr = '1'
    if (sensor .eq. 'HIRS-METOPC') metop_nr = '3'
    
    select case(sensor)
    
    case('MTSAT-1')
     list = chn
      sensor_string = 'mtsat_1_imager'
    
    
    case('MTSAT-2')
      list = chn
      sensor_string = 'mtsat_2_imager'
    
    case( 'MODIS-AQUA')
      !list = chn -19
      ! if ( chn .gt. 26) list = chn - 20
      list = chn
      sensor_string = 'eos_2_modis'
      
    case( 'MODIS-TERRA')
      !list = chn -19
      ! if ( chn .gt. 26) list = chn - 20
      list = chn
      sensor_string = 'eos_1_modis'
      
    case ( 'VIIRS-SNPP')
      
      chn_list(1:32) = [6,9,3,4,10,14,15,1,2,-1,-1,-1,-1,-1,7,-1,-1,-1,-1,17,  &
                        -1,18,-1,-1,-1,12,-1,-1,19,-1,20,22]
      list = chn_list(chn)
      sensor_string = 'jpss_0_viirs'
    
    case ( 'VIIRS-N20')
      
      chn_list(1:32) = [6,9,3,4,10,14,15,1,2,-1,-1,-1,-1,-1,7,-1,-1,-1,-1,17,  &
                        -1,18,-1,-1,-1,12,-1,-1,19,-1,20,22]
      list = chn_list(chn)
      sensor_string = 'noaa_20_viirs'
      
    case ('AHI8')
      chn_list(20:38) = [7,-1,-1,-1,-1,-1,-1,9,10,11,12,14,15,16,-1,-1,-1,8,13]
      list = chn_list(chn)
      sensor_string = 'himawari_8_ahi'
    
    case ('AHI9')
      chn_list(20:38) = [7,-1,-1,-1,-1,-1,-1,9,10,11,12,14,15,16,-1,-1,-1,8,13]
      list = chn_list(chn)
      sensor_string = 'himawari_9_ahi'  
    
    case ('GOES-8','GOES-9')
      
      
      chn_list(20) = 1
      chn_list(27) = 2
      chn_list(31) = 3
      chn_list(32) = 4
      list = chn_list(chn)
      sensor_string = 'goes_'//sensor(6:6)//'_imager'
      rttov_version_string = '7'
      max_satzen = 75.
   
   
    case ('GOES-10','GOES-11')
      
      
      chn_list(20) = 1
      chn_list(27) = 2
      chn_list(31) = 3
      chn_list(32) = 4
      list = chn_list(chn)
      sensor_string = 'goes_'//sensor(6:7)//'_imager'
      rttov_version_string = '8'
      max_satzen = 75.
   
     
    case ('GOES-12','GOES-13','GOES-14','GOES-15')
      
      
      chn_list(20) = 1
      chn_list(27) = 2
      chn_list(31) = 3
      chn_list(33) = 4
      
      list = chn_list(chn)
      sensor_string = 'goes_'//sensor(6:7)//'_imager'
      rttov_version_string = '8'
      max_satzen = 75.
     
     
      
    case ('GOES-16')
      chn_list(1:45) = [2,3,1,-1,-1,5,6,-1,-1,-1,-1,-1,-1,-1,-1,  &
                        -1,-1,-1,-1,7,-1,7,-1,-1,-1,4,9,10,11,12,  &
                        14,15,16,-1,-1,-1,8,13,-1,-1,-1,-1,-1,-1,-1]
      list = chn_list(chn)
      sensor_string = 'goes_16_abi'
      
    case ('GOES-17')
      chn_list(1:45) = [2,3,1,-1,-1,5,6,-1,-1,-1,-1,-1,-1,-1,-1,  &
                        -1,-1,-1,-1,7,-1,7,-1,-1,-1,4,9,10,11,12,  &
                        14,15,16,-1,-1,-1,8,13,-1,-1,-1,-1,-1,-1,-1]
      list = chn_list(chn)
      sensor_string = 'goes_17_abi'

    case ('GOES-18')
      chn_list(1:45) = [2,3,1,-1,-1,5,6,-1,-1,-1,-1,-1,-1,-1,-1,  &
                        -1,-1,-1,-1,7,-1,7,-1,-1,-1,4,9,10,11,12,  &
                        14,15,16,-1,-1,-1,8,13,-1,-1,-1,-1,-1,-1,-1]
      list = chn_list(chn)
      sensor_string = 'goes_18_abi'
      
    case ( 'SEVIRI-MSG08')
      !---stw chn_list(20:33) = [1,-1,-1,-1,-1,-1,-1,2,3,4,5,6,7,8] 
      !chn_list(20:33) = [4,-1,-1,-1,-1,-1,-1,5,6,7,8,9,10,11]
      chn_list(20:37) = [4,-1,-1,-1,-1,-1,-1,-1,6,7,8,9,10,11,-1,-1,-1,5]
      list = chn_list(chn)
      sensor_string = 'msg_1_seviri'
    case ( 'SEVIRI-MSG09')
      !---stw chn_list(20:33) = [1,-1,-1,-1,-1,-1,-1,2,3,4,5,6,7,8]
      !chn_list(20:33) = [4,-1,-1,-1,-1,-1,-1,5,6,7,8,9,10,11]
      chn_list(20:37) = [4,-1,-1,-1,-1,-1,-1,-1,6,7,8,9,10,11,-1,-1,-1,5]
      list = chn_list(chn)
      sensor_string = 'msg_2_seviri'  
    case ( 'SEVIRI-MSG10')
      !---stw chn_list(20:33) = [1,-1,-1,-1,-1,-1,-1,2,3,4,5,6,7,8]
      !chn_list(20:33) = [4,-1,-1,-1,-1,-1,-1,5,6,7,8,9,10,11]
      chn_list(20:37) = [4,-1,-1,-1,-1,-1,-1,-1,6,7,8,9,10,11,-1,-1,-1,5]
      list = chn_list(chn)
      sensor_string = 'msg_3_seviri'  
    case ( 'SEVIRI-MSG11')
      !---stw chn_list(20:33) = [1,-1,-1,-1,-1,-1,-1,2,3,4,5,6,7,8]
      !chn_list(20:33) = [4,-1,-1,-1,-1,-1,-1,5,6,7,8,9,10,11]
      chn_list(20:37) = [4,-1,-1,-1,-1,-1,-1,-1,6,7,8,9,10,11,-1,-1,-1,5]
      list = chn_list(chn)
      sensor_string = 'msg_4_seviri'  
      
    case('OLCI')
      chn_list(1:32) = [6,9,3,4,10,14,15,1,2,-1,-1,-1,-1,-1,7,  &
                        -1,-1,-1,-1,7,-1,18,-1,-1,-1,12,-1,-1,19,-1,  &
                        20,22]
      list = chn
      list =  (chn-1)*5*3 + 1
      sensor_string = 'sentinel3_1_olci'
      
    case ('SLSTR2')
        chn_list(1:32) = [6,9,3,4,10,14,15,1,2,-1,-1,-1,-1,-1,7,  &
                          -1,-1,-1,-1,7,-1,18,-1,-1,-1,12,-1,-1,19,-1,  &
                          20,22]
        list = chn
        sensor_string = 'sentinel3_2_slstr'
        
    case('FY3-D')
        chn_list(20:32) = [ 20,-1,-1,21,-1,-1,-1,-1, 22,23,-1,24,25 ]
        list = chn_list(chn)
        sensor_string = 'fy3_4_mersi2'   
   
    case default
      if (index(sensor,'AVHRR')  .gt. 0) then
       
        chn_list(1) = 1
        chn_list(2) = 2
        select case(sensor(11:12))
        ! -- sensor name can be AVHRR-NOAA??
        case ('15','16','17','18','19')
        
          chn_list(6) = 3
          chn_list(20) = 4
          chn_list(31) = 5
          chn_list(32) = 6
        
        case default
        
          chn_list(20) = 3
          chn_list(31) = 4
          chn_list(32) = 5
        end select 
        
        if (index(sensor,'AVHRR-METOP')  .gt. 0) then
          chn_list(6) = 3
          chn_list(20) = 4
          chn_list(31) = 5
          chn_list(32) = 6
        end if
        
        
        sensor_string = 'noaa_'//sensor(11:12)//'_avhrr'
       
        if ( index(sensor, 'NOAA0' ) .gt. 0 ) then
          sensor_string = 'noaa_'//sensor(12:12)//'_avhrr'
          rttov_version_string = '8'
           chn_list(20) = 1
          chn_list(31) = 2
          chn_list(32) = 3
        end if
        list = chn_list(chn)
       
        if (index(sensor,'AVHRR-METOP')  .gt. 0 ) then
          sensor_string = 'metop_'//metop_nr//'_avhrr'
        end if
        if (index(sensor,'AVHRR-TIROSN') .gt. 0) then 
          sensor_string = 'noaa_5_avhrr'
            rttov_version_string = '8'
           chn_list(20) = 1
          chn_list(31) = 2
          chn_list(32) = 3
        end if
        
      else if (index(sensor,'HIRS') .gt. 0) then 
                          !  20      23 24 25    27 28   30 31 32 33 34 35 36 
        chn_list(20:36) = [ 19,-1,-1,18,15,14,-1,12,11,-1,9, 8,10, 7, 6, 5, 4 ]
        list = chn_list(chn)
        rttov_version_string = '8'
        ! possible is HIRS-NOAA?? or HIRS-METOP?
        
        if ( index(sensor,'HIRS-NOAA') .gt. 0 ) then
        
         sensor_string = 'noaa_'//sensor(10:11)//'_hirs-shifted'
        
         read ( sensor(10:10), '(I1)') nr_hirs
         if ( nr_hirs .eq. 0 ) sensor_string = 'noaa_'//sensor(11:11)//'_hirs-shifted'
        end if
        
        if (index(sensor,'HIRS-METOP')  .gt. 0 ) then
          sensor_string = 'metop_'//metop_nr//'_hirs-shifted'
        end if
        if (index(sensor,'HIRS-TIROSN') .gt. 0) then 
          sensor_string = 'noaa_5_hirs'
        end if
        
        
        
      else if  (index(sensor,'GOES-') .gt. 0)  then
        rttov_version_string = '8'
        
        chn_list(1) = 1
        chn_list(20) = 2
        chn_list(27) = 3
        chn_list(31) = 4
        chn_list(32) = 5
        chn_list(33) = 6
        sensor_string = 'goes_'//sensor(6:7)//'_imager'
        if ( sensor(6:6) .ne. '1' ) then
           sensor_string = 'goes_'//sensor(6:6)//'_imager'
           rttov_version_string = '7'
        end if   
        list = chn_list(chn)
       
      else
         print*,trim(sensor)
        stop ' missing sensor in channel mapping of RTTOV inform andi.walther@ssec.wisc.edu ' 
       end if 
    end select
    
    
    
    coef_filename = trim(path)//'/rtcoef_rttov12/rttov'//rttov_version_string//'pred54L/rtcoef_'//trim(sensor_string)//'.dat' 
    cld_coef_filename = trim(path)//'cldaer_ir/sccldcoef_'//trim(sensor_string)//'.dat'

   if ( list .eq. 0) list = -1
    !print*,'Sensor rttov mapping: ',sensor,chn,list
    deallocate(chn_list)
    
   if (rttov_version_string .eq. '7') max_satzen = 75.
   if (rttov_version_string .eq. '8') max_satzen = 75.

end function channel_map 


end module cx_rttov_mapping_mod
