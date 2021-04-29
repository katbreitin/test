
case(783) 

call MESG("Using Solar Rtm Terms for MODIS/TERRA channels", level = verb_lev % DEFAULT)

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch01.txt
Solar_Rtm%Tau_H2O_Coef( 1,:) = (/  8.5404e-05,  4.0251e-03, -1.1151e-04/)
Solar_Rtm%Tau_Ray( 1) =  5.5616e-02
Solar_Rtm%Tau_O3( 1) =  3.5758e-05
Solar_Rtm%Tau_CO2( 1) =  0.0000e+00
Solar_Rtm%Tau_CH4( 1) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch02.txt
Solar_Rtm%Tau_H2O_Coef( 2,:) = (/  5.4917e-04,  5.1843e-03, -2.1348e-04/)
Solar_Rtm%Tau_Ray( 2) =  1.7624e-02
Solar_Rtm%Tau_O3( 2) =  0.0000e+00
Solar_Rtm%Tau_CO2( 2) =  0.0000e+00
Solar_Rtm%Tau_CH4( 2) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch03.txt
Solar_Rtm%Tau_H2O_Coef( 3,:) = (/  4.7599e-06,  5.8822e-05, -2.5171e-07/)
Solar_Rtm%Tau_Ray( 3) =  2.1076e-01
Solar_Rtm%Tau_O3( 3) =  0.0000e+00
Solar_Rtm%Tau_CO2( 3) =  0.0000e+00
Solar_Rtm%Tau_CH4( 3) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch04.txt
Solar_Rtm%Tau_H2O_Coef( 4,:) = (/  1.7264e-05,  3.8090e-04, -1.7925e-06/)
Solar_Rtm%Tau_Ray( 4) =  1.0281e-01
Solar_Rtm%Tau_O3( 4) =  6.2764e-05
Solar_Rtm%Tau_CO2( 4) =  0.0000e+00
Solar_Rtm%Tau_CH4( 4) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch05.txt
Solar_Rtm%Tau_H2O_Coef( 5,:) = (/ -1.1228e-04,  4.1359e-03, -1.4436e-04/)
Solar_Rtm%Tau_Ray( 5) =  3.9315e-03
Solar_Rtm%Tau_O3( 5) =  0.0000e+00
Solar_Rtm%Tau_CO2( 5) =  2.8210e-04
Solar_Rtm%Tau_CH4( 5) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch06.txt
Solar_Rtm%Tau_H2O_Coef( 6,:) = (/ -3.9291e-05,  9.7133e-04, -9.0888e-06/)
Solar_Rtm%Tau_Ray( 6) =  1.3216e-03
Solar_Rtm%Tau_O3( 6) =  0.0000e+00
Solar_Rtm%Tau_CO2( 6) =  6.2040e-03
Solar_Rtm%Tau_CH4( 6) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch07.txt
Solar_Rtm%Tau_H2O_Coef( 7,:) = (/  1.4918e-03,  1.5622e-02, -7.2262e-04/)
Solar_Rtm%Tau_Ray( 7) =  4.6304e-04
Solar_Rtm%Tau_O3( 7) =  0.0000e+00
Solar_Rtm%Tau_CO2( 7) =  1.7455e-02
Solar_Rtm%Tau_CH4( 7) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch08.txt
Solar_Rtm%Tau_H2O_Coef( 8,:) = (/  0.0000e+00,  0.0000e+00,  0.0000e+00/)
Solar_Rtm%Tau_Ray( 8) =  3.4724e-01
Solar_Rtm%Tau_O3( 8) =  0.0000e+00
Solar_Rtm%Tau_CO2( 8) =  0.0000e+00
Solar_Rtm%Tau_CH4( 8) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch09.txt
Solar_Rtm%Tau_H2O_Coef( 9,:) = (/  1.7814e-05,  2.3789e-04, -1.0758e-06/)
Solar_Rtm%Tau_Ray( 9) =  2.5856e-01
Solar_Rtm%Tau_O3( 9) =  0.0000e+00
Solar_Rtm%Tau_CO2( 9) =  0.0000e+00
Solar_Rtm%Tau_CH4( 9) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch10.txt
Solar_Rtm%Tau_H2O_Coef(10,:) = (/  4.3265e-06,  7.6645e-05, -2.8930e-07/)
Solar_Rtm%Tau_Ray(10) =  1.7316e-01
Solar_Rtm%Tau_O3(10) =  0.0000e+00
Solar_Rtm%Tau_CO2(10) =  0.0000e+00
Solar_Rtm%Tau_CH4(10) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch11.txt
Solar_Rtm%Tau_H2O_Coef(11,:) = (/  7.5798e-06,  7.8288e-05,  3.5567e-07/)
Solar_Rtm%Tau_Ray(11) =  1.2330e-01
Solar_Rtm%Tau_O3(11) =  2.6849e-05
Solar_Rtm%Tau_CO2(11) =  0.0000e+00
Solar_Rtm%Tau_CH4(11) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch12.txt
Solar_Rtm%Tau_H2O_Coef(12,:) = (/  7.0371e-05,  9.5413e-04,  1.2186e-06/)
Solar_Rtm%Tau_Ray(12) =  1.0806e-01
Solar_Rtm%Tau_O3(12) =  5.1557e-05
Solar_Rtm%Tau_CO2(12) =  0.0000e+00
Solar_Rtm%Tau_CH4(12) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch13.txt
Solar_Rtm%Tau_H2O_Coef(13,:) = (/ -2.9797e-05,  6.4930e-04, -3.8415e-06/)
Solar_Rtm%Tau_Ray(13) =  4.8894e-02
Solar_Rtm%Tau_O3(13) =  8.5586e-06
Solar_Rtm%Tau_CO2(13) =  0.0000e+00
Solar_Rtm%Tau_CH4(13) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch14.txt
Solar_Rtm%Tau_H2O_Coef(14,:) = (/  2.2625e-06,  7.2146e-05,  1.9143e-07/)
Solar_Rtm%Tau_Ray(14) =  4.5397e-02
Solar_Rtm%Tau_O3(14) =  1.9191e-06
Solar_Rtm%Tau_CO2(14) =  0.0000e+00
Solar_Rtm%Tau_CH4(14) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch15.txt
Solar_Rtm%Tau_H2O_Coef(15,:) = (/ -7.5286e-05,  1.6887e-03, -2.1257e-05/)
Solar_Rtm%Tau_Ray(15) =  3.0540e-02
Solar_Rtm%Tau_O3(15) =  0.0000e+00
Solar_Rtm%Tau_CO2(15) =  0.0000e+00
Solar_Rtm%Tau_CH4(15) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch16.txt
Solar_Rtm%Tau_H2O_Coef(16,:) = (/ -1.6675e-05,  1.8989e-04, -1.2827e-06/)
Solar_Rtm%Tau_Ray(16) =  1.6671e-02
Solar_Rtm%Tau_O3(16) =  0.0000e+00
Solar_Rtm%Tau_CO2(16) =  0.0000e+00
Solar_Rtm%Tau_CH4(16) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch17.txt
Solar_Rtm%Tau_H2O_Coef(17,:) = (/  2.8697e-02,  8.3486e-02, -5.5397e-03/)
Solar_Rtm%Tau_Ray(17) =  1.4152e-02
Solar_Rtm%Tau_O3(17) =  0.0000e+00
Solar_Rtm%Tau_CO2(17) =  0.0000e+00
Solar_Rtm%Tau_CH4(17) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch18.txt
Solar_Rtm%Tau_H2O_Coef(18,:) = (/  1.9965e-01,  3.3896e-01, -2.2512e-02/)
Solar_Rtm%Tau_Ray(18) =  1.2279e-02
Solar_Rtm%Tau_O3(18) =  0.0000e+00
Solar_Rtm%Tau_CO2(18) =  0.0000e+00
Solar_Rtm%Tau_CH4(18) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch19.txt
Solar_Rtm%Tau_H2O_Coef(19,:) = (/  1.1655e-01,  2.2609e-01, -1.5946e-02/)
Solar_Rtm%Tau_Ray(19) =  1.2317e-02
Solar_Rtm%Tau_O3(19) =  0.0000e+00
Solar_Rtm%Tau_CO2(19) =  0.0000e+00
Solar_Rtm%Tau_CH4(19) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch20.txt
Solar_Rtm%Tau_H2O_Coef(20,:) = (/  5.7806e-03,  5.3760e-02, -2.9171e-03/)
Solar_Rtm%Tau_Ray(20) =  3.4858e-05
Solar_Rtm%Tau_O3(20) =  0.0000e+00
Solar_Rtm%Tau_CO2(20) =  1.8412e-03
Solar_Rtm%Tau_CH4(20) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_1_modis-shifted_srf_ch26.txt
Solar_Rtm%Tau_H2O_Coef(26,:) = (/  1.3096e+00,  1.5088e+00, -1.2832e-01/)
Solar_Rtm%Tau_Ray(26) =  2.5640e-03
Solar_Rtm%Tau_O3(26) =  0.0000e+00
Solar_Rtm%Tau_CO2(26) =  2.6963e-04
Solar_Rtm%Tau_CH4(26) =  0.0000e+00

case(784) 

call MESG("Using Solar Rtm Terms for MODIS/AQUA channels", level = verb_lev % DEFAULT)

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch01.txt
Solar_Rtm%Tau_H2O_Coef( 1,:) = (/  9.1944e-05,  3.9721e-03, -9.9794e-05/)
Solar_Rtm%Tau_Ray( 1) =  5.5105e-02
Solar_Rtm%Tau_O3( 1) =  3.4727e-05
Solar_Rtm%Tau_CO2( 1) =  0.0000e+00
Solar_Rtm%Tau_CH4( 1) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch02.txt
Solar_Rtm%Tau_H2O_Coef( 2,:) = (/  5.0118e-04,  4.9854e-03, -2.1292e-04/)
Solar_Rtm%Tau_Ray( 2) =  1.7539e-02
Solar_Rtm%Tau_O3( 2) =  0.0000e+00
Solar_Rtm%Tau_CO2( 2) =  0.0000e+00
Solar_Rtm%Tau_CH4( 2) =  0.0000e+00
 
!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch03.txt
Solar_Rtm%Tau_H2O_Coef( 3,:) = (/  4.9458e-06,  5.8878e-05, -6.1898e-07/)
Solar_Rtm%Tau_Ray( 3) =  2.0808e-01
Solar_Rtm%Tau_O3( 3) =  0.0000e+00
Solar_Rtm%Tau_CO2( 3) =  0.0000e+00
Solar_Rtm%Tau_CH4( 3) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch04.txt
Solar_Rtm%Tau_H2O_Coef( 4,:) = (/  2.9099e-05,  3.5009e-04,  1.7126e-06/)
Solar_Rtm%Tau_Ray( 4) =  1.0266e-01
Solar_Rtm%Tau_O3( 4) =  6.5112e-05
Solar_Rtm%Tau_CO2( 4) =  0.0000e+00
Solar_Rtm%Tau_CH4( 4) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch05.txt
Solar_Rtm%Tau_H2O_Coef( 5,:) = (/ -1.3815e-04,  4.4767e-03, -1.6550e-04/)
Solar_Rtm%Tau_Ray( 5) =  3.9508e-03
Solar_Rtm%Tau_O3( 5) =  0.0000e+00
Solar_Rtm%Tau_CO2( 5) =  3.1365e-04
Solar_Rtm%Tau_CH4( 5) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch06.txt
Solar_Rtm%Tau_H2O_Coef( 6,:) = (/ -3.0591e-05,  9.1273e-04, -7.9645e-06/)
Solar_Rtm%Tau_Ray( 6) =  1.3256e-03
Solar_Rtm%Tau_O3( 6) =  0.0000e+00
Solar_Rtm%Tau_CO2( 6) =  7.0289e-03
Solar_Rtm%Tau_CH4( 6) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch07.txt
Solar_Rtm%Tau_H2O_Coef( 7,:) = (/  1.6901e-03,  1.5514e-02, -7.0496e-04/)
Solar_Rtm%Tau_Ray( 7) =  4.6385e-04
Solar_Rtm%Tau_O3( 7) =  0.0000e+00
Solar_Rtm%Tau_CO2( 7) =  1.7901e-02
Solar_Rtm%Tau_CH4( 7) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch08.txt
Solar_Rtm%Tau_H2O_Coef( 8,:) = (/  0.0000e+00,  0.0000e+00,  0.0000e+00/)
Solar_Rtm%Tau_Ray( 8) =  3.4420e-01
Solar_Rtm%Tau_O3( 8) =  0.0000e+00
Solar_Rtm%Tau_CO2( 8) =  0.0000e+00
Solar_Rtm%Tau_CH4( 8) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch09.txt
Solar_Rtm%Tau_H2O_Coef( 9,:) = (/  1.1672e-05,  2.5062e-04, -2.6846e-06/)
Solar_Rtm%Tau_Ray( 9) =  2.5806e-01
Solar_Rtm%Tau_O3( 9) =  0.0000e+00
Solar_Rtm%Tau_CO2( 9) =  0.0000e+00
Solar_Rtm%Tau_CH4( 9) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch10.txt
Solar_Rtm%Tau_H2O_Coef(10,:) = (/  4.6524e-06,  7.7441e-05, -3.7586e-07/)
Solar_Rtm%Tau_Ray(10) =  1.7448e-01
Solar_Rtm%Tau_O3(10) =  0.0000e+00
Solar_Rtm%Tau_CO2(10) =  0.0000e+00
Solar_Rtm%Tau_CH4(10) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch11.txt
Solar_Rtm%Tau_H2O_Coef(11,:) = (/  5.7824e-06,  7.8402e-05, -3.8829e-07/)
Solar_Rtm%Tau_Ray(11) =  1.2333e-01
Solar_Rtm%Tau_O3(11) =  2.6000e-05
Solar_Rtm%Tau_CO2(11) =  0.0000e+00
Solar_Rtm%Tau_CH4(11) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch12.txt
Solar_Rtm%Tau_H2O_Coef(12,:) = (/  4.2660e-05,  9.5622e-04, -5.8313e-06/)
Solar_Rtm%Tau_Ray(12) =  1.0888e-01
Solar_Rtm%Tau_O3(12) =  5.4690e-05
Solar_Rtm%Tau_CO2(12) =  0.0000e+00
Solar_Rtm%Tau_CH4(12) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch13.txt
Solar_Rtm%Tau_H2O_Coef(13,:) = (/ -2.7545e-05,  5.9578e-04, -1.8256e-06/)
Solar_Rtm%Tau_Ray(13) =  4.8567e-02
Solar_Rtm%Tau_O3(13) =  7.2553e-06
Solar_Rtm%Tau_CO2(13) =  0.0000e+00
Solar_Rtm%Tau_CH4(13) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch14.txt
Solar_Rtm%Tau_H2O_Coef(14,:) = (/  4.8441e-07,  9.6940e-05, -8.6319e-07/)
Solar_Rtm%Tau_Ray(14) =  4.5382e-02
Solar_Rtm%Tau_O3(14) =  2.2643e-06
Solar_Rtm%Tau_CO2(14) =  0.0000e+00
Solar_Rtm%Tau_CH4(14) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch15.txt
Solar_Rtm%Tau_H2O_Coef(15,:) = (/ -1.0378e-04,  1.5823e-03, -2.3177e-05/)
Solar_Rtm%Tau_Ray(15) =  3.0495e-02
Solar_Rtm%Tau_O3(15) =  0.0000e+00
Solar_Rtm%Tau_CO2(15) =  0.0000e+00
Solar_Rtm%Tau_CH4(15) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch16.txt
Solar_Rtm%Tau_H2O_Coef(16,:) = (/ -1.2608e-05,  1.5852e-04,  4.4121e-07/)
Solar_Rtm%Tau_Ray(16) =  1.6679e-02
Solar_Rtm%Tau_O3(16) =  0.0000e+00
Solar_Rtm%Tau_CO2(16) =  0.0000e+00
Solar_Rtm%Tau_CH4(16) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch17.txt
Solar_Rtm%Tau_H2O_Coef(17,:) = (/  2.7199e-02,  8.5923e-02, -6.0260e-03/)
Solar_Rtm%Tau_Ray(17) =  1.4032e-02
Solar_Rtm%Tau_O3(17) =  0.0000e+00
Solar_Rtm%Tau_CO2(17) =  0.0000e+00
Solar_Rtm%Tau_CH4(17) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch18.txt
Solar_Rtm%Tau_H2O_Coef(18,:) = (/  2.1209e-01,  3.4083e-01, -2.2445e-02/)
Solar_Rtm%Tau_Ray(18) =  1.2233e-02
Solar_Rtm%Tau_O3(18) =  0.0000e+00
Solar_Rtm%Tau_CO2(18) =  0.0000e+00
Solar_Rtm%Tau_CH4(18) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch19.txt
Solar_Rtm%Tau_H2O_Coef(19,:) = (/  1.1926e-01,  2.2369e-01, -1.5565e-02/)
Solar_Rtm%Tau_Ray(19) =  1.2285e-02
Solar_Rtm%Tau_O3(19) =  0.0000e+00
Solar_Rtm%Tau_CO2(19) =  0.0000e+00
Solar_Rtm%Tau_CH4(19) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch20.txt
Solar_Rtm%Tau_H2O_Coef(20,:) = (/  6.4333e-03,  5.7078e-02, -3.1772e-03/)
Solar_Rtm%Tau_Ray(20) =  3.5831e-05
Solar_Rtm%Tau_O3(20) =  0.0000e+00
Solar_Rtm%Tau_CO2(20) =  1.7988e-03
Solar_Rtm%Tau_CH4(20) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_eos_2_modis-shifted_srf_ch26.txt
Solar_Rtm%Tau_H2O_Coef(26,:) = (/  1.2485e+00,  1.5270e+00, -1.4480e-01/)
Solar_Rtm%Tau_Ray(26) =  2.5859e-03
Solar_Rtm%Tau_O3(26) =  0.0000e+00
Solar_Rtm%Tau_CO2(26) =  2.0147e-05
Solar_Rtm%Tau_CH4(26) =  0.0000e+00
