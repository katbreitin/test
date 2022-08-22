! $Id: cx_pfaast_rtm_test.f90 3721 2020-02-18 20:25:45Z awalther $
program cx_pfaast_rtm_test
   use cx_pfaast_mod, only:compute_transmission_pfaast
   use cx_pfaast_constants_mod, only: tstd, ostd, wstd
   implicit none
   integer, parameter :: N_PROFILE = 101
   character(len = 200 ) :: ancil_data_path  
   real  :: temp (N_PROFILE)
   real  :: wvmr(N_PROFILE)
   real  :: ozmr(N_PROFILE)
   real  :: theta 
   character (len =40 ) :: sensor
   integer  :: kban_in 
   logical :: use_modis_channel_equivalent
   real  :: taut (N_PROFILE)  
   character (len =40 ) :: sensor_list(4)
   integer :: ii, jj
   integer:: channels_to_check (14)
   ancil_data_path = '/Users/awalther/DATA/Ancil_data/clavrx_ancil_data/'
   temp = tstd
   wvmr = wstd
   ozmr = ostd
   theta = 62.
   sensor='VIIRS'
   sensor='MODIS-AQUA'
   sensor = 'GOES-16'
   kban_in = 32
   use_modis_channel_equivalent  = .true.
   
   sensor_list(1) =  'VIIRS     '
   sensor_list(2) =  'MODIS-AQUA'
   sensor_list(3) =  'GOES-16   '
   sensor_list(4) =  'FY4-A   '
   
   channels_to_check = [1,0,1,0,0,0,0,1,1,1,0,1,1,1] ! this is for FY4a
   
   print*,'+++++++ PFAAST RTM TRANSMISSION COMPARISON ++++++++++++ '
   
   do ii =  20,33
      if (  channels_to_check(ii-19)  .eq. 0) cycle
      print*,'+++++++++++++++++++++++'
      print*,'channel: ',ii
      do jj = 1,4
   
         call compute_transmission_pfaast ( &
            ancil_data_path &
            & ,temp &
            & ,wvmr &
            & ,ozmr & 
            & ,theta  &
            & ,sensor_list(jj) &
            & ,ii &
            & ,taut &
            & , use_modis_channel_equivalent )
            
           
            write(*,'(2A,2x,f7.4)')  'Tottransm for ' &
               , trim(sensor_list(jj)), taut(101)
           
            
       
         end do
     end do

end program cx_pfaast_rtm_test
