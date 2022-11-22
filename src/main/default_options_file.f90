module default_options_file
implicit none
contains
subroutine print_default_options()
WRITE(*,"(A)") "/apollo/cloud/Ancil_Data/clavrx_ancil_data/"
WRITE(*,"(A)") "./temp_files/"
WRITE(*,"(A)",advance="no") "9   !Expert mode (0 - nothing is read below, 1 - only !E1 and lower is read, 9 - everything is read "
WRITE(*,"(A)") ")"
WRITE(*,"(A)",advance="no") "5   !E1 Messaging Control Flag (QUIET = 0,ERROR = 1,MINIMAL = 2,WARNING = 4,DEFAULT = 5,VERBOSE = 9)"
WRITE(*,"(A)") ""
WRITE(*,"(A)") "1   !E1 ALG Mask mode bayesian cloud mask (0=Baseline CM, 1 = ECM1, 2 = ECM2)"
WRITE(*,"(A)") "3   !E1 ALG DCOMP"
WRITE(*,"(A)") "default !E1 ALG ACHA (off, default or select combinations like 110_120_133)"
WRITE(*,"(A)") "2   !E1 ALG CCL Mode (0 =off,1=top only,2 = top + base,3 = top + base + lower)"
WRITE(*,"(A)") "0   !E1 ALG CCL Type = (0=NOAT, 1=ISSCP, 2=NCEP)"
WRITE(*,"(A)") "0   !E1 ALG ASOS Mode (0=off, 1=on)"
WRITE(*,"(A)") "0   !E1 ALG NLCOMP (0=off, 1=on)"
WRITE(*,"(A)") "0   !E1 ALG Aerosol (0=off, 1=on)"
WRITE(*,"(A)") "1   !E2 OUTPUT Level-2  file output flag (0= no, 1 = yes)"
WRITE(*,"(A)") "0   !E2 OUTPUT format flag (0= hdf4, 1 = netcdf4)"
WRITE(*,"(A)") "1   !E2 PRC Cloud flag  "
WRITE(*,"(A)") "200 !E2 Num scan lines"
WRITE(*,"(A)") "0   !E2 SASRAB switch on/off SASRAB parameters "
WRITE(*,"(A)",advance="no") "1   !E2 NWP Nwp Model Option  (0=off, 1=gfs,2=ncep reanalysis,3=cfsr,4= gdas, 5=merra,6=era,7=gfs ai"
WRITE(*,"(A)") "t,8=gfs fv3)"
WRITE(*,"(A)") "0   !E2 NWP Nwp Mode  (0=minimal, 1=all)"
WRITE(*,"(A)") "1   !E2 RTM rtm option for clear-sky trans  ((0=crtm not installed yet ),1=pfaast, 2=rttov)"
WRITE(*,"(A)",advance="no") "0   !E2 NAV nav option  (0 =  Use  level-1b,1 = external (future) 2 = Reposnx ( fred nagle): For non"
WRITE(*,"(A)") "-AVHRR, use 0) "
WRITE(*,"(A)") "1   !E2 OUT output compression flag (0=no,1=gzip)"
WRITE(*,"(A)",advance="no") "2   !E2 MASK read auxilary cloud mask 1b (0 = don't read, 1 = read from 1b and use, 2 = read and sav"
WRITE(*,"(A)") "e as aux, 3 = read modawg/mvcm, but compute ecm2 type/phase)"
WRITE(*,"(A)") "default"
WRITE(*,"(A)") "1   !E3 SFC seebor emiss option (0=UMD, 1=RTTOV 2=SEEBOR) "
WRITE(*,"(A)") "1   !E3 SFC sea emiss flag (0=no, 1=yes) "
WRITE(*,"(A)") "0   !E3 SFC read hires sfc type flag (0=no-8km, 1 = yes-1km) "
WRITE(*,"(A)") "1   !E3 SFC read land mask flag (0=no, 1=yes) (goge2_0ll.hdf)"
WRITE(*,"(A)") "1   !E3 SFC read coast mask flag (0=no, 1=yes) (coast_mask_xkm.hdf, x = 1 or 8)"
WRITE(*,"(A)") "1   !E3 SFC read surface elevation flag (0=no, 1=yes) (GLOBE_xkm_digelev.hdf, x = 1 or 8)"
WRITE(*,"(A)") "0   !E3 SFC read volcano mask flag (0=no, 1=yes) (volcano_mask_1km.hdf)"
WRITE(*,"(A)") "1   !E3 SFC read snow mask flag (0=no, 1 = ims, 2 = GlobSnow)"
WRITE(*,"(A)") "0   !E3 SFC read dark composite flag (0=no, 1 = yes) only GEO"
WRITE(*,"(A)",advance="no") "0   !E4 AVHRR-ONLY specific ref_cal_1b flag (0 = use default level1b cal, 1 = recalibrate using PATM"
WRITE(*,"(A)") "OS-x numbers (AVHRR-only))"
WRITE(*,"(A)",advance="no") "0   !E4 AVHRR-ONLY therm_cal_1b flag (0 = use default level1b cal, 1= recalibrate using PATMOS-x num"
WRITE(*,"(A)") "bers (AVHRR-only))"
WRITE(*,"(A)") "1   !E5 MASK lrc_flag (0=no,1=yes) "
WRITE(*,"(A)") "1   !E5 NWP smooth nwp flag (0=no, 1=yes)   "
WRITE(*,"(A)") "0   !E5 process_undetected_flag (0=no,1=yes) ( means process all pixels cloudy and cloud-free"
WRITE(*,"(A)",advance="no") "0 5.0 15.0 -100.0 -90.0 0.0 90.0 0.0 180.0 WI    !E6 on_flag Lat_South,Lat_North,Lon_West,Lon_East,Z"
WRITE(*,"(A)") "en_Min,Zen_Max,Solzen_Min,Solzen_Max, Name"
WRITE(*,"(A)") "2 1 1 !E6 native res sample(0)/average(1),average+stats(2), X_stride(>=1), Y_Stride (>=1)"
WRITE(*,"(A)") "1 1 1 1 1 1    !E6 chan on flags of channels 1,2,3,4,5,6"
WRITE(*,"(A)") "1 0 0 0 0 0    !E6 chan on flags of channels 7,8,9,10,11,12"
WRITE(*,"(A)") "0 0 0 0 0 0    !E6 chan on flags of channels 13,14,15,16,17,18"
WRITE(*,"(A)") "0 1 0 0 0 0    !E6 chan on flags of channels 19,20,21,22,23,24"
WRITE(*,"(A)") "0 1 1 1 1 1    !E6 chan on flags of channels 25,26,27,28,29,30"
WRITE(*,"(A)") "1 1 1 0 0 0    !E6 chan on flags of channels 31,32,33,34,35,36"
WRITE(*,"(A)",advance="no") "1 1 0 0 0 0    !E6 chan on flags of channels 37(AHI/ABI-only),38(AHI/ABI-only),39(I1),40(I2),41(I3),"
WRITE(*,"(A)") "42(I4)"
WRITE(*,"(A)") "0 0 0 0 0 0    !E6 chan on flags of channels 43(I5),44(DNB),45-48(Spare)"
WRITE(*,"(A)") "270 !ISCCP-NG WMO 270=g16, 271=g17, 173=HIM8, 70 = MET11, 55 = MET8"
end subroutine print_default_options
end module default_options_file
