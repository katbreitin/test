program cx_rttov_test

   use cx_rttov_bridge_mod
  use cx_rttov_constants_mod

IMPLICIT NONE
   
   
   
   integer, parameter :: N_PROFILE = 101
   character(len = 200 ) :: ancil_data_path  
   real  :: temp (N_PROFILE)
   real  :: wvmr(N_PROFILE)
   real  :: ozmr(N_PROFILE)
   
    integer, parameter :: n_ppp = 4000
   real  :: tempx (N_PROFILE,n_ppp)
   real  :: wvmrx (N_PROFILE,n_ppp)
   real  :: ozmrx (N_PROFILE,n_ppp)
   real :: pstdx (N_PROFILE,n_ppp)
   real  :: thetax (N_PPP)
   real  :: theta 
   character (len =40 ) :: sensor
   integer  :: kban_in 
   logical :: use_modis_channel_equivalent
   real  :: taut (N_PROFILE) , taut_rttov (N_PROFILE) 
   real  :: tautx  (N_PROFILE,n_ppp) , taut_rttovx  (N_PROFILE,n_ppp)
   character (len =40 ) :: sensor_list(4)
   integer :: ii, jj,kk, WMO_ID
   integer:: channels_to_check (24)
   
  character(len=2) :: ii_str, kk_str
   real :: kk_f(10)
  
  
  
  
   ancil_data_path = '/Users/awalther/DATA/Ancil_data/clavrx_ancil_data/'
   temp = tstd  * 0.5
   wvmr = wstd  * 0.0000001
   ozmr = 4.*ostd 
   theta = 0.
   sensor='MODIS-AQUA'
   WMO_ID = 784
   
   
   pstdx = spread(pstd,2,N_PPP)
   tempx = spread(temp,2,N_PPP)
   wvmrx = spread(wvmr,2,N_PPP)
   ozmrx = spread(ozmr,2,N_PPP)
   thetax = spread(theta,1,N_PPP)
   
   kban_in = 32
   use_modis_channel_equivalent  = .true.
  ii = 30

    call compute_transmission_rttov ( &
            ancil_data_path &
            &, pstdx &
            & ,tempx &
            & ,wvmrx &
            & ,ozmrx & 
            & ,thetax  &
            & ,sensor &
            & ,WMO_ID &
            & ,ii &
            & ,taut_rttovx &
            & , use_modis_channel_equivalent ) 
            
            
            write(*,'(1A,I2,2A,2x,f7.4)')  ' RTTOV Tot transmission for channel ' ,ii ,' ' &
               , trim(sensor), taut_rttovx(101,10)
               
               print*,taut_rttovx(:,10)
               
end program cx_rttov_test
