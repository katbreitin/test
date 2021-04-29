!$Id: akh_cloud_type_clavrx_include_3.inc 3190 2019-03-11 22:04:41Z heidinger $
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
     Rad_H2O_BB_Profile => Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(27)%Rad_BB_Cloud_Profile
     if (Use_10_4_Flag) then
         Rad_IRWIN_BB_Profile => Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(31)%Rad_BB_Cloud_Profile
     else
         Rad_IRWIN_BB_Profile => Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(38)%Rad_BB_Cloud_Profile
     endif
     T_Prof => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%T_prof
     Z_Prof => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Z_prof
     Tropo_Level => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
     Sfc_Level => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Sfc_Level
