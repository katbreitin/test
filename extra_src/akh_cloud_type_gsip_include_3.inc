!$Id: akh_cloud_type_gsip_include_3.inc 2685 2018-05-04 14:17:47Z awalther $
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
     Rad_H2O_BB_Profile => Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%Cloud_Prof_Ch9
     Rad_11um_BB_Profile => Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%Cloud_Prof_Ch14
     T_Prof => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%T_prof
     Z_Prof => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Z_prof
     Tropo_Level => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
     Sfc_Level => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Sfc_Level
