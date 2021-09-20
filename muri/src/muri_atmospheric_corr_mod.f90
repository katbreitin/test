! $Id:$
!
!
module muri_atmospheric_corr_mod
  implicit none
  type muri_atm_corr_type
    real :: hallo
  
  end type muri_atm_corr_type
  integer, parameter :: N_BAND= 6 
  real, parameter :: TauO3(N_BAND) = [3.9e-3,1.4e-2,2.88e-2,8.28e-4,0.0,0.0]
  real, parameter :: TauH2O(N_BAND)=[8.94e-5,9.52e-4,3.78e-3,5.73e-3,2.05e-3,5.65e-3]
  real, parameter :: TauDry(N_BAND)=[1.25e-3,9.5e-4,3.91e-3,2.0e-5,9.98e-3,1.63e-2]
  
!  ABI 
!  real, parameter :: TauO3(N_BAND) = [2.90e-3,2.52e-2,8.1e-2,0.0,0.0,2e-5]
!  real, parameter :: TauH2O(N_BAND)= [8.00e-5,5.11e-3,8.61e-3, 5.23e-3, 1.62e-3, 2.53e-2]
!  real, parameter :: TauDry(N_BAND)= [1.25e-3, 3.91e-3, 2.0e-5, 1.69e-2, 9.98e-3, 1.63e-2]

contains
!
!
!   airmass
!   index 1: ozone, 2 = h2o , 3 =dry
function airmass( angle, idx)
  real, intent(in) :: angle
  integer, intent(in) :: idx
  real :: airmass
  real :: AG (3,4)
  real :: P1, P2 , P3
  real, parameter :: pi=3.1415
  real, parameter :: DTR=pi/180
  
  AG(1,:) = (/268.45,0.5,115.42,-3.2922/)
  AG(2,:) = (/0.0311,0.1,92.471,-1.3814/)
  AG(3,:) = (/0.4567,0.07,96.484,-1.697/)
  
  P1=cos(angle * DTR)
  P2=AG(idx,1)*(angle**(AG(idx,2)))
  P3=(AG(idx,3)-angle)**AG(idx,4)
  airmass = (P1+(P2)*P3)**(-1)

end function airmass

function airmass_two_way( sol, sat, idx)
  real, intent(in) :: sol,sat
  integer, intent(in) :: idx
  real:: airmass_ind(2)
   real ::  airmass_two_way
  
  airmass_ind (1) = airmass(sol, idx)
  airmass_ind (2) = airmass(sat, idx)
  
  airmass_two_way = airmass_ind (1) + airmass_ind (2)
  
  
end function airmass_two_way


!
!  output is transmission due to gas
!  
function muri_transmission_default ( sol, sat)
  real, intent(in) :: sol,sat
  real ::  muri_transmission_default(6)
  
  
  integer :: i_gas, i_band
  real :: t_o3(N_BAND), t_h2o(N_BAND), t_dry(N_BAND)
  real :: amass(3)

  do i_gas=1,3
    amass(i_gas) = airmass_two_way (sol,sat,i_gas)
  end do
  
  do i_band =1 , 6
    t_o3(i_band) =  exp ( amass(1) * TauO3 ( i_band))
    t_h2o(i_band) =  exp ( amass(2) * TauH2O ( i_band))
    t_dry(i_band) =  exp ( amass(3) * TauDry ( i_band)) 
  end do
    
 
 muri_transmission_default = T_O3 * T_H2O * T_Dry

end function muri_transmission_default


function muri_transmission ( sol, sat, ozone, h2o_conc)
  real, intent(in) :: sol,sat
  real, intent(in) :: ozone, h2o_conc
  real, parameter ::   O3_K0(6)=(/-1.26E-04, 3.81E-05, 5.62E-04, 3.24E-07, 1.19E-07, -1.07E-08/)
  real, parameter ::   O3_K1(6)=(/ 1.17E-05, 4.08E-05, 8.33E-05, 2.41E-06, 1.03E-25,  4.41E-09/)
  real, parameter ::   H2O_K0(6)=(/-9.65E+00, -7.29E+00, -5.90E+00, -5.52E+00, -6.56E+00, -5.80E+00/)
  real, parameter ::   H2O_K1(6)=(/ 9.86E-01,  9.81E-01,  9.43E-01,  9.52E-01,  1.02E+00,  1.26E+00/)
  real, parameter ::   H2O_K2(6)=(/-6.33E-05, -5.73E-03, -1.77E-02, -2.05E-02, -3.61E-03, -4.78E-03/)
  real :: t_o3(N_BAND), t_h2o(N_BAND), t_dry(N_BAND)
  integer :: i_band, i_gas
  real :: B1, B2
   real :: amass(3)
   real ::  muri_transmission(6)
  
  do i_gas=1,3
    amass(i_gas) = airmass_two_way (sol,sat,i_gas)
  end do
  
  do i_band = 1, 6 
    
    t_o3(i_band)=exp(O3_K0(i_band)+O3_K1(i_band)*amass(1) * ozone)
    if ( ozone .le. 0 ) t_o3(i_band) = 1.
    B1 = H2O_K1(i_band)*log(amass(2)*H2O_conc)
    B2 = H2O_K2(i_band)*(log(amass(2)*H2O_conc))**2
    T_H2O(i_band)=exp(exp(H2O_K0(i_band)+B1+B2 ))
    if ( h2o_conc .le. 0 ) t_h2o(i_band) = 1.
    !print*,exp(H2O_K0(i_band)+B1+B2 )
    t_dry ( i_band) = exp(amass(3)*TauDry(i_band))
  
  
  end do
   muri_transmission = T_O3 * T_H2O * T_Dry
   !muri_transmission =  T_H2O
  
end function muri_transmission


end module muri_atmospheric_corr_mod


!program do_it
 ! use  muri_atmospheric_corr_mod
 ! real :: pp(6)
 
  
 ! print*,'hallo'
!  pp = muri_transmission_default ( 2.3,22.)
!  print*,pp
 ! print*,muri_transmission ( 2.3,22., 220., 0.)
 ! print*,log(0.00001)
 
!end program do_it
