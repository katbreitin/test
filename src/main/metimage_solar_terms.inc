
case(385) 

call MESG("Using Solar Rtm Terms for EPS-SG/METIMAGE channels", level = verb_lev % DEFAULT)

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch01.txt
Solar_Rtm%Tau_H2O_Coef( 1,:) = (/  7.1378e-06,  8.4930e-05,  1.8328e-07/)
Solar_Rtm%Tau_Ray( 1) =  2.5368e-01
Solar_Rtm%Tau_O3( 1) =  0.0000e+00
Solar_Rtm%Tau_CO2( 1) =  0.0000e+00
Solar_Rtm%Tau_CH4( 1) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch02.txt
Solar_Rtm%Tau_H2O_Coef( 2,:) = (/  1.1250e-05,  2.9824e-04, -1.2934e-06/)
Solar_Rtm%Tau_Ray( 2) =  1.0065e-01
Solar_Rtm%Tau_O3( 2) =  6.8435e-05
Solar_Rtm%Tau_CO2( 2) =  0.0000e+00
Solar_Rtm%Tau_CH4( 2) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch03.txt
Solar_Rtm%Tau_H2O_Coef( 3,:) = (/ -2.7118e-05,  5.7779e-04, -5.0924e-06/)
Solar_Rtm%Tau_Ray( 3) =  4.7690e-02
Solar_Rtm%Tau_O3( 3) =  5.8298e-06
Solar_Rtm%Tau_CO2( 3) =  0.0000e+00
Solar_Rtm%Tau_CH4( 3) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch04.txt
Solar_Rtm%Tau_H2O_Coef( 4,:) = (/  2.0295e-06,  1.2431e-04, -1.1880e-06/)
Solar_Rtm%Tau_Ray( 4) =  2.9266e-02
Solar_Rtm%Tau_O3( 4) =  0.0000e+00
Solar_Rtm%Tau_CO2( 4) =  0.0000e+00
Solar_Rtm%Tau_CH4( 4) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch05.txt
Solar_Rtm%Tau_H2O_Coef( 5,:) = (/ -3.8656e-08,  1.5763e-06,  2.1432e-11/)
Solar_Rtm%Tau_Ray( 5) =  2.8053e-02
Solar_Rtm%Tau_O3( 5) =  0.0000e+00
Solar_Rtm%Tau_CO2( 5) =  0.0000e+00
Solar_Rtm%Tau_CH4( 5) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch06.txt
Solar_Rtm%Tau_H2O_Coef( 6,:) = (/ -1.3168e-05,  8.5250e-04, -1.7485e-05/)
Solar_Rtm%Tau_Ray( 6) =  1.6989e-02
Solar_Rtm%Tau_O3( 6) =  0.0000e+00
Solar_Rtm%Tau_CO2( 6) =  0.0000e+00
Solar_Rtm%Tau_CH4( 6) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch07.txt
Solar_Rtm%Tau_H2O_Coef( 7,:) = (/  3.0802e-02,  1.0655e-01, -6.6504e-03/)
Solar_Rtm%Tau_Ray( 7) =  1.3528e-02
Solar_Rtm%Tau_O3( 7) =  0.0000e+00
Solar_Rtm%Tau_CO2( 7) =  0.0000e+00
Solar_Rtm%Tau_CH4( 7) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch08.txt
Solar_Rtm%Tau_H2O_Coef( 8,:) = (/ -1.7107e-04,  3.1347e-03, -8.3740e-05/)
Solar_Rtm%Tau_Ray( 8) =  3.9419e-03
Solar_Rtm%Tau_O3( 8) =  0.0000e+00
Solar_Rtm%Tau_CO2( 8) =  1.3081e-04
Solar_Rtm%Tau_CH4( 8) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch09.txt
Solar_Rtm%Tau_H2O_Coef( 9,:) = (/  1.4096e+00,  1.5281e+00, -1.3893e-01/)
Solar_Rtm%Tau_Ray( 9) =  2.6217e-03
Solar_Rtm%Tau_O3( 9) =  0.0000e+00
Solar_Rtm%Tau_CO2( 9) =  7.0903e-05
Solar_Rtm%Tau_CH4( 9) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch10.txt
Solar_Rtm%Tau_H2O_Coef(10,:) = (/ -4.1269e-05,  7.8077e-04, -8.3041e-06/)
Solar_Rtm%Tau_Ray(10) =  1.3258e-03
Solar_Rtm%Tau_O3(10) =  0.0000e+00
Solar_Rtm%Tau_CO2(10) =  2.0988e-03
Solar_Rtm%Tau_CH4(10) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch11.txt
Solar_Rtm%Tau_H2O_Coef(11,:) = (/ -1.3366e-04,  1.9972e-03, -7.0048e-05/)
Solar_Rtm%Tau_Ray(11) =  3.5687e-04
Solar_Rtm%Tau_O3(11) =  0.0000e+00
Solar_Rtm%Tau_CO2(11) =  2.6098e-06
Solar_Rtm%Tau_CH4(11) =  0.0000e+00

!SRF = tigr_tpw_trans_rtcoef_metopsg_1_metimage_srf_ch12.txt
Solar_Rtm%Tau_H2O_Coef(12,:) = (/  1.0116e-02,  6.9878e-02, -3.9749e-03/)
Solar_Rtm%Tau_Ray(12) =  4.1031e-05
Solar_Rtm%Tau_O3(12) =  0.0000e+00
Solar_Rtm%Tau_CO2(12) =  1.2007e-03
Solar_Rtm%Tau_CH4(12) =  0.0000e+00
