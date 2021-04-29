case(55) 

call MESG("Solar Rtm Terms for MSG1 are not available, using MSG3", level = verb_lev % DEFAULT)

!SRF = tigr_tpw_trans_ch1_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef( 1,:) = (/  5.6845e-05,  2.7160e-03, -6.9465e-05/)
Solar_Rtm%Tau_Ray( 1) =  5.8816e-02
Solar_Rtm%Tau_O3( 1) =  4.7724e-05
Solar_Rtm%Tau_CO2( 1) =  0.0000e+00
Solar_Rtm%Tau_CH4( 1) =  0.0000e+00

!SRF = tigr_tpw_trans_ch2_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef( 2,:) = (/  6.7371e-03,  2.8608e-02, -1.7104e-03/)
Solar_Rtm%Tau_Ray( 2) =  2.2491e-02
Solar_Rtm%Tau_O3( 2) =  0.0000e+00
Solar_Rtm%Tau_CO2( 2) =  0.0000e+00
Solar_Rtm%Tau_CH4( 2) =  0.0000e+00

!SRF = tigr_tpw_trans_ch3_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef( 6,:) = (/ -3.1209e-05,  3.6450e-03, -1.2004e-04/)
Solar_Rtm%Tau_Ray( 6) =  1.3155e-03
Solar_Rtm%Tau_O3( 6) =  0.0000e+00
Solar_Rtm%Tau_CO2( 6) =  1.2421e-02
Solar_Rtm%Tau_CH4( 6) =  0.0000e+00

!SRF = tigr_tpw_trans_ch4_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef(20,:) = (/  4.3171e-03,  3.1576e-02, -1.7931e-03/)
Solar_Rtm%Tau_Ray(20) =  2.9281e-05
Solar_Rtm%Tau_O3(20) =  0.0000e+00
Solar_Rtm%Tau_CO2(20) =  6.1788e-02
Solar_Rtm%Tau_CH4(20) =  0.0000e+00

case(56) 

call MESG("Solar Rtm Terms for MSG2 are not available, using MSG3", level = verb_lev % DEFAULT)
!SRF = tigr_tpw_trans_ch1_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef( 1,:) = (/  5.6845e-05,  2.7160e-03, -6.9465e-05/)
Solar_Rtm%Tau_Ray( 1) =  5.8816e-02
Solar_Rtm%Tau_O3( 1) =  4.7724e-05
Solar_Rtm%Tau_CO2( 1) =  0.0000e+00
Solar_Rtm%Tau_CH4( 1) =  0.0000e+00

!SRF = tigr_tpw_trans_ch2_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef( 2,:) = (/  6.7371e-03,  2.8608e-02, -1.7104e-03/)
Solar_Rtm%Tau_Ray( 2) =  2.2491e-02
Solar_Rtm%Tau_O3( 2) =  0.0000e+00
Solar_Rtm%Tau_CO2( 2) =  0.0000e+00
Solar_Rtm%Tau_CH4( 2) =  0.0000e+00

!SRF = tigr_tpw_trans_ch3_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef( 6,:) = (/ -3.1209e-05,  3.6450e-03, -1.2004e-04/)
Solar_Rtm%Tau_Ray( 6) =  1.3155e-03
Solar_Rtm%Tau_O3( 6) =  0.0000e+00
Solar_Rtm%Tau_CO2( 6) =  1.2421e-02
Solar_Rtm%Tau_CH4( 6) =  0.0000e+00

!SRF = tigr_tpw_trans_ch4_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef(20,:) = (/  4.3171e-03,  3.1576e-02, -1.7931e-03/)
Solar_Rtm%Tau_Ray(20) =  2.9281e-05
Solar_Rtm%Tau_O3(20) =  0.0000e+00
Solar_Rtm%Tau_CO2(20) =  6.1788e-02
Solar_Rtm%Tau_CH4(20) =  0.0000e+00

case(57) 

!SRF = tigr_tpw_trans_ch1_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef( 1,:) = (/  5.6845e-05,  2.7160e-03, -6.9465e-05/)
Solar_Rtm%Tau_Ray( 1) =  5.8816e-02
Solar_Rtm%Tau_O3( 1) =  4.7724e-05
Solar_Rtm%Tau_CO2( 1) =  0.0000e+00
Solar_Rtm%Tau_CH4( 1) =  0.0000e+00

!SRF = tigr_tpw_trans_ch2_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef( 2,:) = (/  6.7371e-03,  2.8608e-02, -1.7104e-03/)
Solar_Rtm%Tau_Ray( 2) =  2.2491e-02
Solar_Rtm%Tau_O3( 2) =  0.0000e+00
Solar_Rtm%Tau_CO2( 2) =  0.0000e+00
Solar_Rtm%Tau_CH4( 2) =  0.0000e+00

!SRF = tigr_tpw_trans_ch3_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef( 6,:) = (/ -3.1209e-05,  3.6450e-03, -1.2004e-04/)
Solar_Rtm%Tau_Ray( 6) =  1.3155e-03
Solar_Rtm%Tau_O3( 6) =  0.0000e+00
Solar_Rtm%Tau_CO2( 6) =  1.2421e-02
Solar_Rtm%Tau_CH4( 6) =  0.0000e+00

!SRF = tigr_tpw_trans_ch4_resp_msg3.dat
Solar_Rtm%Tau_H2O_Coef(20,:) = (/  4.3171e-03,  3.1576e-02, -1.7931e-03/)
Solar_Rtm%Tau_Ray(20) =  2.9281e-05
Solar_Rtm%Tau_O3(20) =  0.0000e+00
Solar_Rtm%Tau_CO2(20) =  6.1788e-02
Solar_Rtm%Tau_CH4(20) =  0.0000e+00

case(70) 

!SRF = tigr_tpw_trans_ch1_resp_msg4.dat
Solar_Rtm%Tau_H2O_Coef( 1,:) = (/  9.9339e-05,  2.5188e-03, -5.5689e-05/)
Solar_Rtm%Tau_Ray( 1) =  5.8254e-02
Solar_Rtm%Tau_O3( 1) =  4.6354e-05
Solar_Rtm%Tau_CO2( 1) =  0.0000e+00
Solar_Rtm%Tau_CH4( 1) =  0.0000e+00

!SRF = tigr_tpw_trans_ch2_resp_msg4.dat
Solar_Rtm%Tau_H2O_Coef( 2,:) = (/  7.0146e-03,  2.8576e-02, -1.6470e-03/)
Solar_Rtm%Tau_Ray( 2) =  2.2315e-02
Solar_Rtm%Tau_O3( 2) =  0.0000e+00
Solar_Rtm%Tau_CO2( 2) =  0.0000e+00
Solar_Rtm%Tau_CH4( 2) =  0.0000e+00

!SRF = tigr_tpw_trans_ch3_resp_msg4.dat
Solar_Rtm%Tau_H2O_Coef( 6,:) = (/ -2.2287e-04,  3.5342e-03, -1.0897e-04/)
Solar_Rtm%Tau_Ray( 6) =  1.2963e-03
Solar_Rtm%Tau_O3( 6) =  0.0000e+00
Solar_Rtm%Tau_CO2( 6) =  1.2207e-02
Solar_Rtm%Tau_CH4( 6) =  0.0000e+00

!SRF = tigr_tpw_trans_ch4_resp_msg4.dat
Solar_Rtm%Tau_H2O_Coef(20,:) = (/  4.4831e-03,  3.2630e-02, -1.8233e-03/)
Solar_Rtm%Tau_Ray(20) =  3.0177e-05
Solar_Rtm%Tau_O3(20) =  0.0000e+00
Solar_Rtm%Tau_CO2(20) =  5.2239e-02
Solar_Rtm%Tau_CH4(20) =  0.0000e+00
