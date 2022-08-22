case(171) 

call MESG("Solar Rtm Terms for MTSAT-1R are not available, using MTSAT-2", level = verb_lev % DEFAULT)

!SRF = tigr_tpw_trans_ch1_resp_mtsat.dat
Solar_Rtm%Tau_H2O_Coef( 1,:) = (/  2.3152e-03,  1.2706e-02, -6.6864e-04/)
Solar_Rtm%Tau_Ray( 1) =  4.3586e-02
Solar_Rtm%Tau_O3( 1) =  2.4278e-05
Solar_Rtm%Tau_CO2( 1) =  0.0000e+00
Solar_Rtm%Tau_CH4( 1) =  0.0000e+00

!SRF = tigr_tpw_trans_ch5_resp_mtsat.dat
Solar_Rtm%Tau_H2O_Coef(20,:) = (/  6.4207e-03,  5.0112e-02, -2.7538e-03/)
Solar_Rtm%Tau_Ray(20) =  4.3870e-05
Solar_Rtm%Tau_O3(20) =  0.0000e+00
Solar_Rtm%Tau_CO2(20) =  1.8098e-03
Solar_Rtm%Tau_CH4(20) =  0.0000e+00

case(172) 

!SRF = tigr_tpw_trans_ch1_resp_mtsat.dat
Solar_Rtm%Tau_H2O_Coef( 1,:) = (/  2.3152e-03,  1.2706e-02, -6.6864e-04/)
Solar_Rtm%Tau_Ray( 1) =  4.3586e-02
Solar_Rtm%Tau_O3( 1) =  2.4278e-05
Solar_Rtm%Tau_CO2( 1) =  0.0000e+00
Solar_Rtm%Tau_CH4( 1) =  0.0000e+00

!SRF = tigr_tpw_trans_ch5_resp_mtsat.dat
Solar_Rtm%Tau_H2O_Coef(20,:) = (/  6.4207e-03,  5.0112e-02, -2.7538e-03/)
Solar_Rtm%Tau_Ray(20) =  4.3870e-05
Solar_Rtm%Tau_O3(20) =  0.0000e+00
Solar_Rtm%Tau_CO2(20) =  1.8098e-03
Solar_Rtm%Tau_CH4(20) =  0.0000e+00
