! $Id: rt_utilities_mod.f90 4122 2021-03-24 13:08:20Z awalther $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: rt_utilities_mod.f90 (src)
!       RT_UTILITIES_MOD (program)
!
! PURPOSE: Perform needed Radiative Transfer Functions
!
! DESCRIPTION: CLAVR-x uses much radiative transfer.  This module stores and serves
!              the RT variables to other modules.  It serves as the interface
!              to external RT models (PFAAST) and holds convenient RT-specific
!              routines
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! Description of RTM Structure Members given in rtm_common.f90
!
! CLAVR-x has 44 channels.
! Channels 1-36 are MODIS
! Channels 37-38 are ABI channels not on MODIS
! Channels 39-43 are the VIIRS I-bands
! Channel 44 is the VIIRS DNB
!
! Not all members of the RTM structure are populated for all channels.
! However, all members are allocated for any active cell
!
! Here is the current implementation
!
! There are 6 types of channels (FIX THIS)
! 1. MODIS IR-only channels = 21-36 excluding 26.
! 2. Non-MODIS IR-only channels  37-38
! 3. MODIS (Solar + IR) channels = Channel 20
! 4. Supported Solar Channels = Channels 1,2,5,6,7 and 44
! 5. Unsupported Solar Channels = Channels 3,4,8-19,26,39-41
! 6. Unsupported IR Channels = Channels 42 and 43
!
! For each type described above, the following profiles are made:
! 1. Rad_Atm_Prof, Trans_Atm_Profile, Rad_BB_Cloud_Profile
! 2. All Profiles
! 3. Trans_Atm_Profile, Trans_Atm_Solar_Profile, Trans_Atm_Total_Profile
! 4. No Profiles
! 5. No Profiles
!
! Viirs I-band support is limited (channels 39-43) and is being developed. No
! pixel-level toa clear-sky fields are generated for the I-band variables yet.
!
!--------------------------------------------------------------------------------------
module RT_UTILITIES_MOD

   use CONSTANTS_MOD, only: &
       Real4 &
      , Int4 &
      , Int1 &
      , Sym &
      , Dtor &
      , Missing_Value_Int1  &
      , Missing_Value_Real4  &
      , Exe_Prompt &
      , G &
      , PI &
      , SOLAR_OBS_TYPE &
      , LUNAR_OBS_TYPE &
      , MIXED_OBS_TYPE &
      , THERMAL_OBS_TYPE

   use NWP_COMMON_MOD, only: &
      NWP &
      , delta_t_inversion &
      , p_inversion_min


   use PIXEL_COMMON_MOD, only: &
        NWP_PIX &
      , Sensor &
      , Ch &
      , Geo &
      , Sfc &
      , Image &
      , rtm_opt &
      , Bad_Pixel_Mask &
      , Zen_Idx_Rtm &
      , Solar_Rtm &
      , Trans_Atm_Ch20_Solar_Rtm &
      , Trans_Atm_Ch20_Solar_Total_Rtm &
      , Beta_11um_12um_Tropo_Rtm &
      , Beta_11um_104um_Tropo_Rtm &
      , Beta_11um_85um_Tropo_Rtm &
      , Beta_11um_67um_Tropo_Rtm &
      , Beta_11um_133um_Tropo_Rtm &
      , Ancil_Data_Dir  &
      , Pc_Opaque_Cloud &
      , Zc_Opaque_Cloud &
      , Beta_104um_12um_Tropo_Rtm &
      , Beta_104um_85um_Tropo_Rtm &
      , Beta_104um_67um_Tropo_Rtm &
      , Beta_104um_133um_Tropo_Rtm &
      , Beta_104um_11um_Tropo_Rtm

   use NUMERICAL_ROUTINES_MOD , only: &
       LOCATE

   use CX_DATE_TIME_TOOLS_MOD , only: &
       COMPUTE_TIME_HOURS

   use CX_SCIENCE_TOOLS_MOD , only: &
       VAPOR

   use PLANCK_MOD, only: &
        PLANCK_RAD_FAST &
      , PLANCK_TEMP_FAST

   use CALIBRATION_CONSTANTS_MOD, only: &
       SOLAR_CH20_NU

   use RTM_COMMON_MOD, only: &
       NLEVELS_RTM &
      , Rtm &
      , P_Std_Rtm &
      , T_Std_Rtm &
      , Wvmr_Std_Rtm &
      , Ozmr_Std_Rtm

   use CX_PFAAST_MOD, only: &
        PSTD

   use CLAVRX_MESSAGE_MOD, only: MESG, verb_lev

   use CX_RTM_MOD, only: &
      cx_rtm_input &
      , cx_calculate_rtm

   use CX_REAL_BOOLEAN_MOD

   use cx_nwp_rtm_mod,only: &
      convert_atmos_prof_nwp_rtm &
      , Wvmr_Nwp

   implicit none
   private
   private:: EMISSIVITY, &
             BETA_RATIO, &
             INVERSION_PROFILE, &
             SENSOR_NAME_FOR_RTM, &
             COMPUTE_CHANNEL_ATM_DWN_SFC_RAD, &
             gamma_factor

   public::  GET_PIXEL_NWP_RTM


    integer, parameter, public:: RTM_NVZEN = 50

    integer, parameter :: Chan_Idx_Min = 1
    integer, parameter :: Chan_Idx_Max = 45

    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Trans_Atm_Solar_Prof
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Trans_Atm_Total_Prof
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Rad_Atm_Prof
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Rad_Atm_Dwn_Prof
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Rad_BB_Cloud_Prof

    integer, parameter:: Ilon_Stride = 0
    integer, parameter:: Ilat_Stride = 0
    integer, parameter:: Ivza_Stride = 0





    character(len=20),  save:: Sc_Name_Rtm

    real, parameter::  Rtm_Vza_Binsize = 0.02

contains






   !====================================================================
   ! subroutine Name: COMPUTE_CLEAR_RAD_PROFILES_RTM
   !
   ! Function:
   ! Computes clear sky radiance profiles
   !
   !  used global variables
   !
   !
   !
   !====================================================================
   subroutine COMPUTE_CLEAR_RAD_PROFILES_RTM(x_nwp,y_nwp,z_nwp)

      integer :: x_nwp,y_nwp,z_nwp
      integer:: Chan_Idx
      integer:: Lev_Idx
      real:: T_mean
      real:: B_mean
      real:: B_Level
      real:: Opd_Layer
      real:: Trans_Layer
      real:: Trans_Total
      real :: trans_local (NLEVELS_RTM)
      real :: t_prof_local(NLEVELS_RTM)


      !--- upwelling profiles
      Rad_Atm_Prof = Missing_Value_Real4
      Rad_BB_Cloud_Prof = Missing_Value_Real4



      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max


         if (Ch(Chan_Idx) %Obs_Type /= THERMAL_OBS_TYPE .and. &
             Ch(Chan_Idx)%Obs_Type /= MIXED_OBS_TYPE) cycle
         if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle
         trans_local = rtm (x_nwp,y_nwp) % d(z_nwp) % ch (chan_idx) % Trans_Atm_Profile
         t_prof_local = rtm (x_nwp,y_nwp) % t_prof
         Rad_Atm_Prof(1,Chan_Idx) = 0.0
         Rad_BB_Cloud_Prof(1,Chan_Idx) = 0.0

         do Lev_Idx = 2, NLevels_Rtm

            T_mean = 0.5*(t_prof_local(Lev_Idx-1) + t_prof_local(Lev_Idx))

            B_mean = PLANCK_RAD_FAST(Chan_Idx,T_mean)

            Rad_Atm_Prof(Lev_Idx,Chan_Idx) = Rad_Atm_Prof(Lev_Idx-1,Chan_Idx) +  &
              (trans_local(Lev_Idx-1) - trans_local(Lev_Idx)) * B_mean

            B_Level = PLANCK_RAD_FAST(Chan_Idx,T_Prof_local(Lev_Idx))

            Rad_BB_Cloud_Prof(Lev_Idx,Chan_Idx) = Rad_Atm_Prof(Lev_Idx,Chan_Idx) +  &
                                          (trans_local(lev_idx) * B_Level)




         end do
      end do


      !--- downwelling profiles
      Rad_Atm_Dwn_Prof = Missing_Value_Real4

      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max

         if (Chan_Idx == 31 .OR. Chan_Idx == 38) then
           if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle
           Trans_Total = 1.0
           Rad_Atm_Dwn_Prof(1,Chan_Idx) = 0.0

           do Lev_Idx = 2, Nlevels_Rtm

             T_mean = 0.5*(t_prof_local(Lev_Idx) + t_prof_local(Lev_Idx-1))
             B_mean = PLANCK_RAD_FAST(Chan_Idx,T_mean)

             Opd_Layer = -1.0 * log(trans_local(lev_idx)/trans_local(lev_idx-1))
             Opd_Layer = max(0.0,Opd_Layer)

             Trans_Layer = exp(-1.0*Opd_Layer)

             Rad_Atm_Dwn_Prof(Lev_Idx,Chan_Idx) = Trans_Total * Rad_Atm_Dwn_Prof(Lev_Idx-1,Chan_Idx) +  &
                                          (1.0-Trans_Layer) * B_mean

             Trans_Total = Trans_Total * Trans_Layer

           end do

         else

           cycle

         endif

      end do

   end subroutine COMPUTE_CLEAR_RAD_PROFILES_RTM



   !====================================================================
   ! subroutine Name: GET_PIXEL_NWP_RTM
   !
   ! Function:
   ! Calculates the RTM clear-sky transmittance profiles for a given segment.
   !
   !  Called in process_clavrx
   !      CHANGED to RTTOV Sep 2020
   !
   !   This is called in
   !
   !
   !====================================================================
   subroutine GET_PIXEL_NWP_RTM(Line_Idx_Min,Num_Lines)
      implicit none
      integer, intent(in):: Line_Idx_Min
      integer, intent(in):: Num_Lines

      integer:: Elem_Idx
      integer:: Line_Idx
      integer:: Sfc_Level_Idx
      integer:: x_nwp
      integer:: y_nwp
      integer:: z_nwp
      real:: Prof_Weight
      integer:: Error_Status
      integer:: Chan_Idx
      type(cx_rtm_input) :: rtm_inp
      real,allocatable :: trans_prof_rtm_chn(:,:)

      real :: Trans_profile(nlevels_rtm)

      ! make vectors for RTTOV
      integer :: nwp_size_arr(2)
      integer :: n_val_pixels
      integer :: ii_pixel
      real :: geo_2way_term


      real, dimension(NLevels_Rtm) ::  &
                    T_Prof_Rtm, &
                     Z_Prof_Rtm,  &
                     Wvmr_Prof_Rtm,  &
                     Ozmr_Prof_Rtm, &
                     Tpw_Prof_Rtm



      logical, allocatable :: is_valid_pixel(:,:)


      ! ===================   EXECUTABLE START ======================================

      allocate ( is_valid_pixel (1:Image%Number_Of_Elements,Line_Idx_Min: (Num_Lines + Line_Idx_Min - 1)))
      ! make a check if we have at least one pixel which is not bad or space
      is_valid_pixel = Bad_Pixel_Mask .NE. sym%YES .and. .not. Geo%Space_Mask



      if ( .not. ANY ( is_valid_pixel)) then

         !call MESG('only bad pixels at RTM entree point ...',level = verb_lev %DEFAULT)
         return

      end if




      ! - find RTM pixels from NWP and allocate rtm (lon,lat) % d (zen) )
      !  these pixels are the used for RTTOV
      !--- loop over pixels in segment
      line_loop1: do Line_Idx = Line_Idx_Min, Num_Lines + Line_Idx_Min - 1
         element_loop1: do Elem_Idx = 1, Image%Number_Of_Elements

            if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle
            if (Geo%Space_Mask(Elem_Idx,Line_Idx) ) cycle

            !--- compute viewing zenith bin for Rtm calculation
            Zen_Idx_Rtm(Elem_Idx,Line_Idx) =  &
              max(1,min(Rtm_Nvzen,ceiling(Geo%Coszen(Elem_Idx,Line_Idx)/Rtm_Vza_Binsize)))

            !--- store cell and angular indices
            x_nwp = NWP_PIX%I_Nwp(Elem_Idx,Line_Idx)
            y_nwp = NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)
            z_nwp = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

            ! - compute all sen-independent variables for a RTM pixel
            if (.not. Rtm(x_nwp,y_nwp)%is_allocated ) then
              !--- allocate Rtm arrays

               call  Rtm(x_nwp,y_nwp) % ALLOC(nlevels_rtm,Rtm_Nvzen)

               ! - get and compute profiles for this NWP/RTM pixel
               !--- compute mixing ratio profile
               call COMPUTE_WVMR_PROFILE_NWP(NWP%P_Std, & ! in
                                     NWP%T_Prof(:,x_nwp,y_nwp), & ! in
                                     NWP%Rh_Prof(:,x_nwp,y_nwp), & ! in
                                     Wvmr_Nwp)  ! out

               !--- compute tpw profiles
               call COMPUTE_TPW_PROFILE_NWP(NWP%P_Std, &  ! in
                                    Wvmr_Nwp,  &  ! in
                                    NWP%Tpw_Prof(:,x_nwp,y_nwp)) ! out

               !--- convert the atmospheric profiles from nwp to Rtm pressure coords
               call CONVERT_ATMOS_PROF_NWP_RTM(NWP%NLevels, &   ! in
                                       NWP%Sfc_Level(x_nwp,y_nwp), &  ! in
                                       NWP%Zsfc(x_nwp,y_nwp), &  ! in
                                       NWP%Tmpair(x_nwp,y_nwp), &  ! in
                                       NWP%Rhsfc(x_nwp,y_nwp), &  ! in
                                       NWP%Psfc(x_nwp,y_nwp), &  ! in
                                       NWP%P_Std, &  ! in
                                       NWP%T_Prof(:,x_nwp,y_nwp), &  ! in
                                       NWP%Z_Prof(:,x_nwp,y_nwp), &  ! in
                                       Wvmr_Nwp, &  ! in
                                       NWP%Ozone_Prof(:,x_nwp,y_nwp), &  ! in
                                       NLevels_Rtm, &  ! in
                                       P_Std_Rtm, &  ! in
                                       T_Prof_Rtm, &  ! out
                                       Z_Prof_Rtm, &  ! out
                                       Wvmr_Prof_Rtm, &  ! out
                                       Ozmr_Prof_Rtm, &  ! out
                                       T_Std_Rtm, &  ! in
                                       Wvmr_Std_Rtm, &  ! in
                                       Ozmr_Std_Rtm)  ! in

               !--- compute tpw profiles
               call COMPUTE_TPW_PROFILE_NWP(P_Std_Rtm, &  ! in
                                    Wvmr_Prof_Rtm,  &  ! in
                                    Tpw_Prof_Rtm) ! out

               !--- store in Rtm structures
              ! T_Prof_Rtm(Lev_Idx) = T_Prof_Rtm
               Rtm(x_nwp,y_nwp)%T_Prof = T_Prof_Rtm
               Rtm(x_nwp,y_nwp)%Z_Prof = Z_Prof_Rtm
               Rtm(x_nwp,y_nwp)%Wvmr_Prof = Wvmr_Prof_Rtm
               Rtm(x_nwp,y_nwp)%Ozmr_Prof = Ozmr_Prof_Rtm
               Rtm(x_nwp,y_nwp)%Tpw_Prof = Tpw_Prof_Rtm




              Rtm(x_nwp,y_nwp)%is_set = .true.
            end if

            ! -check if this viewing bin was already computed
            !  - now also add zen dimension
            if ( .not. Rtm(x_nwp,y_nwp) % d (z_nwp) % is_allocated ) then
                call rtm(x_nwp,y_nwp) % d(z_nwp) % ALLOC(Sensor%Chan_On_Flag_Default,nlevels_rtm)

                Rtm(x_nwp,y_nwp) % Satzen_Mid_Bin = acos((z_nwp-1)*Rtm_Vza_Binsize + Rtm_Vza_Binsize/2.0) / DTOR !-TOCHECK

               if (Rtm(x_nwp,y_nwp)%Vza_Idx_Min == Missing_Value_Int1) Rtm(x_nwp,y_nwp)%Vza_Idx_Min = z_nwp
               if (Rtm(x_nwp,y_nwp)%Vza_Idx_Max == Missing_Value_Int1) Rtm(x_nwp,y_nwp)%Vza_Idx_Max = z_nwp

               if (z_nwp < Rtm(x_nwp,y_nwp)%Vza_Idx_Min) Rtm(x_nwp,y_nwp)%Vza_Idx_Min = z_nwp
               if (z_nwp > Rtm(x_nwp,y_nwp)%Vza_Idx_Max) Rtm(x_nwp,y_nwp)%Vza_Idx_Max = z_nwp

            end if


         end do element_loop1
      end do line_loop1

      ! now we have  All profiles and geometries collected

      ! this is everything for RTTOV

      ! - find number of needed calculations
      rtm_inp % ancil_path = trim(Ancil_Data_Dir)
      n_val_pixels = 0 ! counter
      nwp_size_arr = shape(rtm)

      do x_nwp = 1, nwp_size_arr(1)
        do y_nwp =1, nwp_size_arr(2)
          if ( rtm(x_nwp,y_nwp) % is_allocated ) then

              do z_nwp =1,RTM_NVZEN
                if ( rtm(x_nwp,y_nwp) % d(z_nwp) % is_allocated )  n_val_pixels = n_val_pixels + 1
              end do

          end if
        end do
      end do
      ! print*,'Number of needed RTTOV/PFAAST Clear-Sky Transmission calculations: ',n_val_pixels

      ! - useless in n_val_pixels is 0

      if ( n_val_pixels .eq. 0) then
         print*, 'Number of needed RTM pixels is ZERO ...'
         print*, ' check file and line ', __FILE__, __LINE__
         print*, ' this is an error this should not appear'
         print*,'please inform andi.walther@ssec.wisc.edu'
         stop

      end if

      allocate (rtm_inp % p_std( NLEVELS_RTM,n_val_pixels))
      allocate (rtm_inp % t_prof( NLEVELS_RTM,n_val_pixels))
      allocate (rtm_inp % w_prof( NLEVELS_RTM,n_val_pixels))
      allocate (rtm_inp % o_prof( NLEVELS_RTM,n_val_pixels))
      allocate (rtm_inp % tpw_prof( NLEVELS_RTM,n_val_pixels))
      allocate (rtm_inp % sat_bin( n_val_pixels))

      ! - populate RTTOV input
      n_val_pixels = 0
      do x_nwp = 1, nwp_size_arr(1)
            do y_nwp =1, nwp_size_arr(2)
              if ( rtm(x_nwp,y_nwp) % is_allocated ) then
                do z_nwp =1,RTM_NVZEN
                  if ( rtm(x_nwp,y_nwp) % d(z_nwp) % is_allocated ) then

                    n_val_pixels = n_val_pixels + 1
                    rtm_inp % p_std( :,n_val_pixels) = pstd
                    rtm_inp % t_prof(:, n_val_pixels) = rtm(x_nwp,y_nwp) % t_prof
                    rtm_inp % w_prof( :,n_val_pixels) = max(rtm(x_nwp,y_nwp) % wvmr_prof , 0.2E-5)
                    rtm_inp % o_prof( :,n_val_pixels) = rtm(x_nwp,y_nwp) % ozmr_prof
                    rtm_inp % tpw_prof( :,n_val_pixels) = Rtm(x_nwp,y_nwp)%Tpw_Prof
                    rtm_inp % sat_bin( n_val_pixels) = acos((z_nwp - 1) * Rtm_Vza_Binsize + Rtm_Vza_Binsize/2.0) / DTOR
                  end if
                end do
              endif
          end do
      end do

      do Chan_Idx = Chan_Idx_Min,Chan_Idx_Max

        if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle
        if (ch(Chan_Idx)%Obs_Type == LUNAR_OBS_TYPE) cycle  !  save this for later
           print*,chan_idx,ch(Chan_Idx)%Obs_Type
          select case ( ch(Chan_Idx)%Obs_Type )

          case ( SOLAR_OBS_TYPE)
            allocate ( trans_prof_rtm_chn (NLEVELS_RTM,n_val_pixels) )
            do ii_pixel = 1, n_val_pixels

              call SOLAR_TRANS(rtm_inp % tpw_prof(:, ii_pixel),chan_idx &
                 ,rtm_inp % sat_bin( ii_pixel),trans_profile, error_status)

              trans_prof_rtm_chn (:,ii_pixel) = Trans_profile
            end do

            


          case ( THERMAL_OBS_TYPE , MIXED_OBS_TYPE)

            Sc_Name_Rtm = SENSOR_NAME_FOR_RTM(Sensor%WMO_id,Sensor%Sensor_Name, Chan_Idx)
            rtm_inp % sc_name = sc_name_rtm
            rtm_inp % chan_idx = Chan_Idx
            rtm_inp % which_rtm = rtm_opt



            call CX_CALCULATE_RTM(rtm_inp,trans_prof_rtm_chn)



           case default

          end select


          ! move back
          ii_pixel = 0
          do x_nwp = 1, nwp_size_arr(1)
            do y_nwp =1,nwp_size_arr(2)
              if ( rtm(x_nwp,y_nwp) % is_allocated ) then

                do z_nwp =1,RTM_NVZEN
                  if ( rtm(x_nwp,y_nwp) % d(z_nwp) % is_allocated ) then

                    ii_pixel = ii_pixel + 1
                    rtm(x_nwp,y_nwp) % d(z_nwp) % ch(chan_idx) % Trans_Atm_Profile &
                        & = Trans_Prof_Rtm_chn(:,ii_pixel)** Gamma_Factor(Chan_Idx)
                  end if
                end do
              end if
            end do
          end do
print*,chan_idx,nwp_size_arr
   print*,rtm(180,346) % d(22) % ch(chan_idx) % Trans_Atm_Profile
          deallocate ( Trans_Prof_Rtm_chn)
      end do ! - channel loop
    stop
      deallocate (rtm_inp % p_std)
      deallocate (rtm_inp % t_prof)
      deallocate (rtm_inp % w_prof)
      deallocate (rtm_inp % o_prof)
      deallocate (rtm_inp % sat_bin)

      ! set some key levels

      do x_nwp = 1, nwp_size_arr(1)
        do y_nwp =1, nwp_size_arr(2)
          if ( rtm(x_nwp,y_nwp) % is_set ) then
            call FIND_RTM_LEVELS(x_nwp,y_nwp)
            call INVERSION_PROFILE(Rtm(x_nwp,y_nwp)%T_Prof,Rtm(x_nwp,y_nwp)%Inver_Prof)
          end if
        end do
      end do

      ! now we have IR ( RTTOV) IR Trans_Atm_Profile   profiles




      ! = we have some things t calculate for which we need sensor pixels
      line_loop: do Line_Idx = Line_Idx_Min, Num_Lines + Line_Idx_Min - 1
         element_loop: do Elem_Idx = 1, Image%Number_Of_Elements
            if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle
            if (Geo%Space_Mask(Elem_Idx,Line_Idx) ) cycle

            x_nwp = NWP_PIX%I_Nwp(Elem_Idx,Line_Idx)
            y_nwp = NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)
            z_nwp = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

            ! - first everyhting what we only need for each RTM pixel
            if ( .not. rtm(x_nwp,y_nwp) % d(z_nwp) %  is_set ) then

              geo_2way_term  = ((Geo%Coszen(Elem_Idx,Line_Idx)+Geo%Cossolzen(Elem_Idx,Line_Idx)) &
                                            /Geo%Cossolzen(Elem_Idx,Line_Idx))

              do Chan_Idx = Chan_Idx_Min,Chan_Idx_Max
                if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle
                if (Chan_Idx > 21 .and. Chan_Idx /= 26) cycle

                Trans_Atm_Total_Prof(:,Chan_Idx) = rtm(x_nwp,y_nwp) % d(z_nwp) % ch(chan_idx) % Trans_Atm_Profile **  &
                                                 geo_2way_term

                Trans_Atm_Solar_Prof(:,Chan_Idx) = rtm(x_nwp,y_nwp) % d(z_nwp) % ch(chan_idx) % Trans_Atm_Profile **   &
                                          ( geo_2way_term - 1 )

              end do

              !--- compute profiles of radiance (atm and bb cloud)
              call COMPUTE_CLEAR_RAD_PROFILES_RTM(x_nwp,y_nwp,z_nwp)

              !--- copy local rtm profiles back into global rtm structure
              call COPY_LOCAL_RTM_TO_GLOBAL_RTM_STRUCTURE(x_nwp,y_nwp,z_nwp)

              !---   set mask  to indicate this bin or this cell has been computed
               Rtm(x_nwp,y_nwp)%d(z_nwp)%is_set = .true.
            end if

            ! - now for sensor grid using the rtm structure
            ! - find best surface level
            Sfc_Level_Idx = Rtm(x_nwp,y_nwp)%Sfc_Level
            if ((Sfc%Land(Elem_Idx,Line_Idx) == sym%LAND) .and. ((Sfc%Zsfc(Elem_Idx,Line_Idx) .NER. Missing_Value_Real4))) then
               call LOCATE(Rtm(x_nwp,y_nwp)%Z_Prof,NLevels_Rtm,Sfc%Zsfc(Elem_Idx,Line_Idx),Sfc_Level_Idx)
               ! output here is dummy
            end if

            !-- vertical interp weight
            Prof_Weight = (Sfc%Zsfc(Elem_Idx,Line_Idx) - Rtm(x_nwp,y_nwp)%Z_Prof(Sfc_Level_Idx)) / &
              (Rtm(x_nwp,y_nwp)%Z_Prof(Sfc_Level_Idx+1) - Rtm(x_nwp,y_nwp)%Z_Prof(Sfc_Level_Idx))

            !--- constrain - important when high res topo differs from low res nwp topo
            Prof_Weight = max(0.0,min(1.0,Prof_Weight))

            !--- call routine to compute radiative transfer terms such as Rad_Atm
            !--- Trans_Atm, and clear-sky radinace and brightness temperature
            !--- map global to local to allow for efficient looping
            call COMPUTE_CHANNEL_RT(Sfc_Level_Idx,Prof_Weight,x_nwp,y_nwp,Elem_Idx,Line_Idx,z_nwp)

            !--- compute Ch20 Emissivities
            call COMPUTE_CH20_EMISSIVITY(Elem_Idx,Line_Idx)

            !--- compute Emissivity at tropopause
            call COMPUTE_TROPOPAUSE_EMISSIVITIES(Elem_Idx,Line_Idx,x_nwp,y_nwp,z_nwp )

            !--- compute split-window beta ratio at tropopause
            call COMPUTE_BETA_RATIOES(Elem_Idx,Line_Idx)

         end do element_loop
      end do line_loop
      ! ----------------





      return




   end subroutine GET_PIXEL_NWP_RTM


   !====================================================================
   ! subroutine Name: FIND_RTM_LEVELS
   !
   ! Function:
   ! Finds various key Levels in the RTM (tropopause, etc), for each RTM gridcell
   !
   !====================================================================
   subroutine FIND_RTM_LEVELS(Lon_Idx,Lat_Idx)

      integer (kind=int4), intent(in):: Lon_Idx, Lat_Idx
      integer (kind=int4):: k

      !--- find surface Level - this is closest nwp Level above the actual surface
      Rtm(Lon_Idx,Lat_Idx)%Sfc_Level = NLevels_Rtm
      do k = NLevels_Rtm,1,-1
         if (P_Std_Rtm(k) < NWP%Psfc(Lon_Idx,Lat_Idx)) then
            Rtm(Lon_Idx,Lat_Idx)%Sfc_Level = k
            exit
         end if
      end do
      Rtm(Lon_Idx,Lat_Idx)%Sfc_Level = max ( int(1,kind=1), min(  int(NLevels_Rtm,kind=1) - int(1,kind=1), Rtm(Lon_Idx,Lat_Idx)%Sfc_Level))
      !--------------------------------------------------------------------
      !--- find tropopause Level  based on tropopause pressure
      !--- tropopause is between tropopause_Level and tropopaue_Level + 1
      !--------------------------------------------------------------------
      do k = 1, Rtm(Lon_Idx,Lat_Idx)%Sfc_Level-1
         if ((P_Std_Rtm(k) <= NWP%P_Trop(Lon_Idx,Lat_Idx)) .and. &
            (P_Std_Rtm(k+1) > NWP%P_Trop(Lon_Idx,Lat_Idx))) then
            Rtm(Lon_Idx,Lat_Idx)%Tropo_Level = k
         end if
      end do

      !--- check if tropopause Level found
      if (Rtm(Lon_Idx,Lat_Idx)%Tropo_Level == 0) then
         print *, EXE_PROMPT, "Error, tropopause Level not found"
      end if

      do k = 1, Rtm(Lon_Idx,Lat_Idx)%Sfc_Level-1
         if ((P_Std_Rtm(k) <= 850.0) .and. &
               (P_Std_Rtm(k+1) > 850.0)) then

            Rtm(Lon_Idx,Lat_Idx)%Level850 = k
         endif
      enddo

      do k = 1, Rtm(Lon_Idx,Lat_Idx)%Sfc_Level-1
         if ((P_Std_Rtm(k) <= 440.0) .and. &
            (P_Std_Rtm(k+1) > 440.0)) then

            Rtm(Lon_Idx,Lat_Idx)%Level440 = k
         endif
      enddo

      !---------------------------------------------------------------------
      ! find Inversion Level - highest Level Inversion below trop
      !---------------------------------------------------------------------
      Rtm(Lon_Idx,Lat_Idx)%Inversion_Level = 0
      if (Rtm(Lon_Idx,Lat_Idx)%Tropo_Level > 0 .and. Rtm(Lon_Idx,Lat_Idx)%Sfc_Level > 0) then
         do k = Rtm(Lon_Idx,Lat_Idx)%Tropo_Level,Rtm(Lon_Idx,Lat_Idx)%Sfc_Level-1
            if ((Rtm(Lon_Idx,Lat_Idx)%T_Prof(k) - Rtm(Lon_Idx,Lat_Idx)%T_Prof(k+1) > delta_t_Inversion) .and. &
                        (P_Std_Rtm(k) >= p_Inversion_min)) then
               Rtm(Lon_Idx,Lat_Idx)%Inversion_Level = k
               exit
            endif
         enddo
      endif

   end subroutine FIND_RTM_LEVELS

   !====================================================================
   ! subroutine Name: COMPUTE_WVMR_PROFILE_NWP
   !
   ! Function:
   ! Computes the NWP WVMR profile
   !
   ! Input: Lon_Idx = Longitude Index for NWP cell
   !        Lat_Idx = Latitude Index for NWP cell
   !        Press_Profile = pressure profile (hPa)
   !        Temp_Profile = temperature profile (K)
   !        Rh_Profile = relative humidity profile (%)
   !
   ! Output: Wvmr_Profile = water vapor mixing ratio profile (g/kg)
   !         Tpw Profile = profile of Tpw from level to space (g/m^2)
   !
   !====================================================================
   subroutine COMPUTE_WVMR_PROFILE_NWP(Press_Profile, &
                                    Temp_Profile, &
                                    Rh_Profile, &
                                    Wvmr_Profile)

      real, intent(in), dimension(:):: Temp_Profile
      real, intent(in), dimension(:):: Press_Profile
      real, intent(in), dimension(:):: Rh_Profile
      real, intent(out), dimension(:):: Wvmr_Profile
      integer:: Lev_Idx
      real:: e
      real:: es

      integer:: Nlevels

      Nlevels = size(Press_Profile)

      !--- make Wvmr_Profile for use in RTM
      do Lev_Idx = 1, NLevels
         es = VAPOR(Temp_Profile(Lev_Idx))
         e = Rh_Profile(Lev_Idx) * es / 100.0
         Wvmr_Profile(Lev_Idx) = 1000.0*0.622 * (e / (Press_Profile(Lev_Idx) - e))  !(g/kg)
      end do

   end subroutine COMPUTE_WVMR_PROFILE_NWP

   !====================================================================
   ! subroutine Name: COMPUTE_TPW_PROFILE_NWP
   !
   ! Function:
   ! Computes the NWP TPW profile
   !
   ! Input: Lon_Idx = Longitude Index for NWP cell
   !        Lat_Idx = Latitude Index for NWP cell
   !        Press_Profile = pressure profile (hPa)
   !        Temp_Profile = temperature profile (K)
   !        Rh_Profile = relative humidity profile (%)
   !
   ! Output: Wvmr_Profile = water vapor mixing ratio profile (g/kg)
   !         Tpw Profile = profile of Tpw from level to space (g/m^2)
   !
   !====================================================================

   subroutine COMPUTE_TPW_PROFILE_NWP(Press_Profile, &
                                   Wvmr_Profile,  &
                                   Tpw_Profile)

      real, intent(in), dimension(:):: Press_Profile
      real, intent(in), dimension(:):: Wvmr_Profile
      real, intent(out), dimension(:):: Tpw_Profile

      integer:: Lay_Idx
      real :: w_mean
      real :: u_layer
      integer:: Nlevels

      Nlevels = size(Press_Profile)

      !--- make tpw profile for use in atmospheric correction
      Tpw_Profile(1) = 0.0
      do Lay_Idx = 1, NLevels-1   !layer index
         w_mean = 0.5*(Wvmr_Profile(Lay_Idx+1)+Wvmr_Profile(Lay_Idx)) / 1000.0  !(kg/kg)
         u_layer = (10.0/g)*(Press_Profile(Lay_Idx+1)-Press_Profile(Lay_Idx))*w_mean
         Tpw_Profile(Lay_Idx+1) = Tpw_Profile(Lay_Idx) + u_layer
      end do

   end subroutine COMPUTE_TPW_PROFILE_NWP






   !--------------------------------------------------------------------------------------------------
   !> subroutine NAME: SENSOR_NAME_FOR_RTM
   !!
   !! Description:
   !! Knowing the WMO Satellite Identification Number
   !!
   !--------------------------------------------------------------------------------------------------

   function SENSOR_NAME_FOR_RTM ( wmo_id, sensorname, Chan_Idx ) result ( Sensor_Name_Rtm)

      integer, intent(in) :: wmo_id
      character (len =*) , intent(in) :: sensorname
      integer, intent(in) :: Chan_Idx
      character (len =20 ) ::  Sensor_Name_Rtm
      integer :: i

      select case(WMO_Id)

      case(4) !METOP-A
         Sensor_Name_Rtm = 'AVHRR-METOPA'

      case(3) !METOP-B
         Sensor_Name_Rtm = 'AVHRR-METOPB'

      case(5) !METOP-C
        Sensor_Name_Rtm = 'AVHRR-METOPC'

      case(55) !MSG-8
         Sensor_Name_Rtm = 'SEVIRI-MSG08'

      case(56) !MSG-9
         Sensor_Name_Rtm = 'SEVIRI-MSG09'

      case(57) !MSG-10
         Sensor_Name_Rtm = 'SEVIRI-MSG10'

      case(70) !MSG-11
         Sensor_Name_Rtm = 'SEVIRI-MSG11'

      case(171) !MTSAT-1R
         Sensor_Name_Rtm = 'MTSAT-1'

      case(172) !MTSAT-2
         Sensor_Name_Rtm = 'MTSAT-2'

      case(173) !AHI-8
         Sensor_Name_Rtm = 'AHI8'

      case(174) !AHI-9
         Sensor_Name_Rtm = 'AHI9'

      case(200) !NOAA-8
        Sensor_Name_Rtm = 'AVHRR-NOAA08'

      case(201) !NOAA-9
        Sensor_Name_Rtm = 'AVHRR-NOAA09'

      case(202) !NOAA-10
        Sensor_Name_Rtm = 'AVHRR-NOAA10'

      case(203) !NOAA-11
        Sensor_Name_Rtm = 'AVHRR-NOAA11'

      case(204) !NOAA-12
        Sensor_Name_Rtm = 'AVHRR-NOAA12'

      case(205) !NOAA-14
        Sensor_Name_Rtm = 'AVHRR-NOAA14'

      case(206) !NOAA-15
        Sensor_Name_Rtm = 'AVHRR-NOAA15'

      case(207) !NOAA-16
        Sensor_Name_Rtm = 'AVHRR-NOAA16'

      case(208) !NOAA-17
        Sensor_Name_Rtm = 'AVHRR-NOAA17'

      case(209) !NOAA-18
        Sensor_Name_Rtm = 'AVHRR-NOAA18'

      case(223) !NOAA-19
        Sensor_Name_Rtm = 'AVHRR-NOAA19'

      case(224) !VIIRS - SNPP
        Sensor_Name_Rtm = 'VIIRS-SNPP'

      case(225)  !VIIRS NOAA-20
         Sensor_Name_Rtm = 'VIIRS-N20'

      case(226)  !VIIRS NOAA-21
         Sensor_Name_Rtm = 'VIIRS-N21'

      case(250) !GOES-8
        Sensor_Name_Rtm = 'GOES-6'

       case(251) !GOES-8
        Sensor_Name_Rtm = 'GOES-7'

      case(252) !GOES-8
        Sensor_Name_Rtm = 'GOES-8'

      case(253) !GOES-9
        Sensor_Name_Rtm = 'GOES-9'

      case(254) !GOES-10
        Sensor_Name_Rtm = 'GOES-10'

      case(255) !GOES-11
        Sensor_Name_Rtm = 'GOES-11'

      case(256) !GOES-12
        Sensor_Name_Rtm = 'GOES-12'

      case(257) !GOES-13
        Sensor_Name_Rtm = 'GOES-13'

      case(258) !GOES-14
        Sensor_Name_Rtm = 'GOES-14'

      case(259) !GOES-15
        Sensor_Name_Rtm = 'GOES-15'

      case(270) !GOES-16
        Sensor_Name_Rtm = 'GOES-16'

      case(271) !GOES-17
        Sensor_Name_Rtm = 'GOES-17'

      case(272) !GOES-18
        Sensor_Name_Rtm = 'GOES-18'

      case(706) !NOAA-6
        Sensor_Name_Rtm = 'AVHRR-NOAA06'

      case(707) !NOAA-7
        Sensor_Name_Rtm = 'AVHRR-NOAA07'

      case(708) !NOAA-5
        Sensor_Name_Rtm = 'AVHRR-TIROSN'

      case(783) !MODIS
          Sensor_Name_Rtm = 'MODIS-TERRA'

      case(784) !MODIS
         Sensor_Name_Rtm = 'MODIS-AQUA'

      case(510) !FY2A
         Sensor_Name_Rtm ='FY2-1'

      case(514) !FY2D
         Sensor_Name_Rtm ='FY2-2'

      case(515) !FY2E
         Sensor_Name_Rtm ='FY2-3'

      case(523) !FY3D
         Sensor_Name_Rtm = 'FY3-D'

      case(530) ! FY4-A
         Sensor_Name_Rtm ='FY4-A'

      case(810) !COMS
         Sensor_Name_Rtm ='COMS-1'

      case (840)
          Sensor_Name_Rtm ='EPS-SG'

      case default
         print*,'sensor for WMO number not found in RT Utils  ', WMO_id
         print*,'stopping ... Please fix this in rt_utils.F90'
         print*,' better tell andi.walther@ssec.wisc.edu'
         stop
      end select



      if (trim ( Sensorname) == 'AVHRR-IFF' .or. &
         trim ( Sensorname) == 'AVHRR-FUSION')  then

         !  sensor for channels 21:30 and 33:36 is HIRS
         if ( any ( Chan_Idx ==  [ (i,i=21,30,1) , 33,34,35,36] ) ) then

            ! - for this IFF Sensor_Name_Rtm is initially set to AVHRR-<Satellite>
            select case(WMO_Id)

            case(4) !METOP-A
               Sensor_Name_Rtm = 'HIRS-METOPA'

            case(3) !METOP-B
               Sensor_Name_Rtm = 'HIRS-METOPB'

            case(5) !METOP-C
               Sensor_Name_Rtm = 'HIRS-METOPC'

            case(200) !NOAA-8
               Sensor_Name_Rtm = 'HIRS-NOAA08'

            case(201) !NOAA-9
               Sensor_Name_Rtm = 'HIRS-NOAA09'

            case(202) !NOAA-10
               Sensor_Name_Rtm = 'HIRS-NOAA10'

            case(203) !NOAA-11
               Sensor_Name_Rtm = 'HIRS-NOAA11'

            case(204) !NOAA-12
               Sensor_Name_Rtm = 'HIRS-NOAA12'

            case(205) !NOAA-14
               Sensor_Name_Rtm = 'HIRS-NOAA14'

            case(206) !NOAA-15
               Sensor_Name_Rtm = 'HIRS-NOAA15'

            case(207) !NOAA-16
               Sensor_Name_Rtm = 'HIRS-NOAA16'

            case(208) !NOAA-17
               Sensor_Name_Rtm = 'HIRS-NOAA17'

            case(209) !NOAA-18
               Sensor_Name_Rtm = 'HIRS-NOAA18'

            case(223) !NOAA-19
               Sensor_Name_Rtm = 'HIRS-NOAA19'

            case(706) !NOAA-6
               Sensor_Name_Rtm = 'HIRS-NOAA06'

            case(707) !NOAA-7
               Sensor_Name_Rtm = 'HIRS-NOAA07'

            case(708) !NOAA-5
               Sensor_Name_Rtm = 'HIRS-TIROSN'

            case default
               print*,'sensor for WMO number not found in RT Utils for AVHRR-IFF  ', WMO_id
               print*,'stopping ... Please fix this in rt_utils.F90'
               print*,' better tell andi.walther@ssec.wisc.edu'
               stop
            end select


         end if
      end if

      if ((trim (Sensorname) == 'VGAC') .and. Sensor%Fusion_Flag)  then
         if ( any ( Chan_Idx ==  [27,28,33,34,35,36] ) ) then
              select case(WMO_Id)
                case(224) !VIIRS - SNPP
                   Sensor_Name_Rtm = 'HIRS-METOPA'
                case(225)  !VIIRS NOAA-20
                   Sensor_Name_Rtm = 'HIRS-METOPA'
                case default
                   print*,'sensor for WMO number not found in RT Utils for VGAC  ', WMO_id  
                   print*,'stopping ... Please fix this in rt_utils.F90'
                   print*,' better tell andi.walther@ssec.wisc.edu'
                   stop    
              end select
         end if
      end if
   
      if (trim ( Sensorname) == 'VIIRS-IFF') then

         !  sensor for channels 27:28 and 33:36 is CRISP this is similar to MODIS-AQUA
         if ( any ( Chan_Idx ==  [27,28, 33,34,35,36] ) ) Sensor_Name_Rtm   = 'MODIS-AQUA'

      end if

      if ( trim (sensorname ) == 'VIIRS-NASA' ) then
         ! - check what is with 31,32
         if ( any ( Chan_Idx ==  [23,24,25,27,28,30,33,34,35,36] ) ) Sensor_Name_Rtm   = 'MODIS-AQUA'

      end if


   end function SENSOR_NAME_FOR_RTM

          !--------------------------------------------------------------
        ! Compute Gamma Factor for Radiance Bias Adjustment
        !
   ! based on satellite number, channel number and nwp source,
   ! set the value of Gamma_Trans_Factor to reduce RTM bias
   !
   ! nwp flag  (1=gfs,2=ncep reanalysis,3=cfsr)
   !--------------------------------------------------------------
   function gamma_factor( i_ch)  result(answer)
      integer, intent(in) :: i_ch

      real :: answer  ! Function result declaration

      real :: gamma_trans_factor (45)

      !--- initialize to unity
      Gamma_Trans_Factor = 1.0



      !--- GOES-10
      if (Sensor%WMO_Id == 254) then
         if (NWP_PIX%Nwp_Opt == 3) then
            Gamma_Trans_Factor(20) = 1.25
            Gamma_Trans_Factor(27) = 0.79
            Gamma_Trans_Factor(31) = 1.15
            Gamma_Trans_Factor(31) = 1.05
         endif
       endif

      !--- GOES-11
      if (Sensor%WMO_Id == 255) then
         if (NWP_PIX%Nwp_opt == 3) then
            Gamma_Trans_Factor(20) = 1.35
            Gamma_Trans_Factor(27) = 0.74
            Gamma_Trans_Factor(31) = 1.05
            Gamma_Trans_Factor(32) = 1.05
         endif
      endif

      !--- GOES-12
      if (Sensor%WMO_Id == 256) then
         if (NWP_PIX%Nwp_Opt == 1 .or. NWP_PIX%Nwp_Opt == 3) then    !repeat of cfsr
            Gamma_Trans_Factor(20) = 1.45
            Gamma_Trans_Factor(27) = 0.79
            Gamma_Trans_Factor(31) = 1.15
            Gamma_Trans_Factor(33) = 1.15
         end if
      end if

      !--- GOES-13
      if (Sensor%WMO_Id == 257) then
         if (NWP_PIX%Nwp_Opt == 1 .or. NWP_PIX%Nwp_Opt == 3) then
            Gamma_Trans_Factor(20) = 1.00
            Gamma_Trans_Factor(27) = 0.794
            Gamma_Trans_Factor(31) = 1.075
            Gamma_Trans_Factor(33) = 1.064
         end if
      end if

      !--- GOES-15
      if (Sensor%WMO_Id == 259) then
         Gamma_Trans_Factor(20) = 1.55
         Gamma_Trans_Factor(27) = 0.728
         Gamma_Trans_Factor(31) = 1.05
         Gamma_Trans_Factor(33) = 1.075
      end if

      !--- MET-09
      if (Sensor%WMO_Id == 56) then

         Gamma_Trans_Factor(20) = 1.55
         Gamma_Trans_Factor(37) = 0.96
         Gamma_Trans_Factor(29) = 1.35
         Gamma_Trans_Factor(31) = 0.95
         Gamma_Trans_Factor(32) = 0.95
         Gamma_Trans_Factor(33) = 1.11

      endif

      !--- MTSAT-02
      if (Sensor%WMO_Id == 172) then
         if (NWP_PIX%Nwp_Opt == 1) then
            Gamma_Trans_Factor(20) = 1.25
            Gamma_Trans_Factor(27) = 0.87
            Gamma_Trans_Factor(31) = 1.15
            Gamma_Trans_Factor(32) = 1.15
         end if
      end if

      !--- COMS-1
      if (Sensor%WMO_Id == 810) then
         Gamma_Trans_Factor(20) = 1.25
         Gamma_Trans_Factor(27) = 0.936
         Gamma_Trans_Factor(31) = 1.05
         Gamma_Trans_Factor(32) = 1.075
      end if

      !--- NOAA-19 HIRS
      if (Sensor%WMO_Id == 223) then
         Gamma_Trans_Factor(27) = 0.9083
         Gamma_Trans_Factor(28) = 0.9500
         Gamma_Trans_Factor(33) = 1.1250
         Gamma_Trans_Factor(34) = 1.0438
         Gamma_Trans_Factor(35) = 1.0833
         Gamma_Trans_Factor(36) = 1.0917
      end if

      answer = Gamma_Trans_Factor(i_ch)

    end function gamma_factor


   !--------------------------------------------------------------------------------------------------
   !  set up values for the solar rtm calculations for this sensor
   !--------------------------------------------------------------------------------------------------
   subroutine SOLAR_TRANS(tpw_prof,Chan_Idx,Zen_Ang,trans_profile,Error_Status)
      real, intent(in):: tpw_prof(:)
      integer, intent(in):: Chan_Idx
      real, intent(in):: Zen_Ang
      real, intent(out) :: trans_profile(Nlevels_Rtm)
      integer, intent(out):: Error_Status

      real, dimension(3):: Tau_H2O_Coef
      real:: Tau_O2_Column
      real:: Tau_CO2_Column
      real:: Tau_CH4_Column
      real:: Tau_O3_Column

      real:: Tau_O2
      real:: Tau_CO2
      real:: Tau_CH4
      real:: Tau_O3
      real:: Tau_H2O

      real:: Tpw
      real:: mu
      real:: Tau_Gas

      integer:: Lev_Idx

      Error_Status = 1
       Trans_Profile = 1.0

      if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) return

!     if (Chan_Idx >= 20 .and. Chan_Idx /= 26 .and. Chan_Idx/= 44) return
      if (ch(Chan_Idx)%Obs_Type /= SOLAR_OBS_TYPE .and. ch(Chan_Idx)%Obs_Type /= LUNAR_OBS_TYPE) return



      mu = cos(Zen_Ang*DTOR)

      !--- initialize
      Tau_H2O_Coef = Solar_Rtm%Tau_H2O_Coef(Chan_Idx,:)
      Tau_O2_Column = Solar_Rtm%Tau_O2(Chan_Idx)
      Tau_O3_Column = Solar_Rtm%Tau_O3(Chan_Idx)
      Tau_CO2_Column = Solar_Rtm%Tau_CO2(Chan_Idx)
      Tau_CH4_Column = Solar_Rtm%Tau_CH4(Chan_Idx)

      !--- loop through layers and fill profile
      do Lev_Idx = 1, Nlevels_Rtm
         Tpw = Tpw_Prof(Lev_Idx)
         Tau_H2O = Tau_H2O_Coef(1) + Tau_H2O_Coef(2)*tpw + Tau_H2O_Coef(3)*tpw**2
         Tau_CO2 = Tau_Co2_Column * P_Std_Rtm(Lev_Idx) / P_Std_Rtm(Nlevels_Rtm)
         Tau_Ch4 = Tau_Ch4_Column * P_Std_Rtm(Lev_Idx) / P_Std_Rtm(Nlevels_Rtm)
         Tau_O2 = Tau_O2_Column * P_Std_Rtm(Lev_Idx) / P_Std_Rtm(Nlevels_Rtm)
         Tau_O3 = Tau_O3_Column * P_Std_Rtm(Lev_Idx) / P_Std_Rtm(Nlevels_Rtm)
         Tau_Gas = max(0.0,Tau_H2O + Tau_O3 + Tau_O2 + Tau_CO2 + Tau_CH4)
         trans_profile(Lev_Idx) = exp(-Tau_Gas / mu)

      end do


      Error_Status = 0

   end subroutine SOLAR_TRANS



   !------------------------------------------------------------------------------
   ! Routine to compute some needed radiative transfer terms for the IR channels
   !
   ! Input:  Chan_Idx - number of the channel being used
   !         Sfc_Idx - level of the surface in the profiles
   !         Profile_Weight - interpolation weight for estimated the surface in the
   !                          profiles
   !         Sfc_Emiss - emissivity of the surface for this channel
   !         Sfc_Temp - the temperature of the surface
   !         Rad_Atm_Profile - profile of radiance emitted from level to space
   !         Trans_Atm_Profile - profile of transmissio from level to space
   !
   ! Output: Rad_Atm - total radiance due to atmospheric emission
   !         Trans_Atm - total tranmission due to atmosphere
   !         Rad_Atm_Sfc - total radiance due both atmosphere and surface at TOA
   !         Bt_Atm_Sfc - Rad_Atm_Sfc expressed as a brightness temperature
   !------------------------------------------------------------------------------
   subroutine COMPUTE_CHANNEL_ATM_SFC_RAD_BT( &
                Chan_Idx, &
                Sfc_Idx, &
                Profile_Weight, &
                Sfc_Emiss, &
                Sfc_Temp, &
                Rad_Atm_Profile, &
                Trans_Atm_Profile, &
                Rad_Atm, &
                Trans_Atm, &
                Rad_Atm_Sfc, &
                Bt_Atm_Sfc)

      integer, intent(in):: Chan_Idx
      integer, intent(in):: Sfc_Idx
      real, intent(in):: Profile_Weight
      real, intent(in):: Sfc_Emiss
      real, intent(in):: Sfc_Temp
      real, intent(in), dimension(:):: Rad_Atm_Profile
      real, intent(in), dimension(:):: Trans_Atm_Profile
      real, intent(out):: Rad_Atm
      real, intent(out):: Trans_Atm
      real, intent(out):: Rad_Atm_Sfc
      real, intent(out):: Bt_Atm_Sfc

      real:: Sfc_Rad

      Sfc_Rad = Sfc_Emiss * PLANCK_RAD_FAST(Chan_Idx,Sfc_Temp)

      Rad_Atm = Rad_Atm_Profile(Sfc_Idx) +  &
            (Rad_Atm_Profile(Sfc_Idx+1) - Rad_Atm_Profile(Sfc_Idx)) * Profile_Weight

      Trans_Atm = Trans_Atm_Profile(Sfc_Idx) +  &
              (Trans_Atm_Profile(Sfc_Idx+1) - Trans_Atm_Profile(Sfc_Idx)) * Profile_Weight

      Rad_Atm_Sfc = Rad_Atm + Trans_Atm * Sfc_Rad

      Bt_Atm_Sfc = PLANCK_TEMP_FAST(Chan_Idx,Rad_Atm_Sfc)

   end subroutine COMPUTE_CHANNEL_ATM_SFC_RAD_BT

   !------------------------------------------------------------------------------
   ! Routine to compute some needed radiative transfer terms for the IR channels
   !
   ! Input:  Sfc_Idx - level of the surface in the profiles
   !         Profile_Weight - interpolation weight for estimated the surface in the
   !                          profiles
   !         Rad_Atm_Dwn_Profile - profile of radiance emitted from level to space
   !
   ! Output: Rad_Atm_Dwn_Sfc - total radiance due to atmospheric emission at surface
   !------------------------------------------------------------------------------
   subroutine COMPUTE_CHANNEL_ATM_DWN_SFC_RAD( &
                Sfc_Idx, &
                Profile_Weight, &
                Rad_Atm_Dwn_Profile, &
                Rad_Atm_Dwn_Sfc)

      integer, intent(in):: Sfc_Idx
      real, intent(in):: Profile_Weight
      real, intent(in), dimension(:):: Rad_Atm_Dwn_Profile
      real, intent(out):: Rad_Atm_Dwn_Sfc

      Rad_Atm_Dwn_Sfc = Rad_Atm_Dwn_Profile(Sfc_Idx) +  &
                    (Rad_Atm_Dwn_Profile(Sfc_Idx+1) - Rad_Atm_Dwn_Profile(Sfc_Idx)) * Profile_Weight

   end subroutine COMPUTE_CHANNEL_ATM_DWN_SFC_RAD

   !----------------------------------------------------------------------------------------
   ! This routine computes radiative transfer terms such as Rad_Atm
   ! Trans_Atm, and clear-sky radinace and brightness temperature
   !
   ! Input:  Sfc_Level_Idx - level just above the surface in the profiles
   !         Prof_Weight - interpolation weight for interpolating to surface level
   !         Lon_Idx = longitude index of NWP cell
   !         Lat_Idx = latitude index of NWP cell
   !         Zen_Idx = zenith angle index of RTM profile
   !
   ! Output:  (note this passed through global arrays)
   !         Rad_Atm_ChX_Rtm = Radiance Emitted by Atmosphere in Channel X
   !         Trans_Atm_ChX_Rtm = Transmission by Atmosphere in Channel X
   !         Rad_Clear_ChX_Rtm = Radiance at TOA for clear skies (atm + sfc)
   !         Bt_Clear_ChX_Rtm = BT at TOA for clear skies (atm + sfc)
   !
   ! NOTE:  These pixels are not smoothed and are blocky at NWP scale.  Should fix
   !----------------------------------------------------------------------------------------
   subroutine COMPUTE_CHANNEL_RT(Sfc_Level_Idx,Prof_Weight,Lon_Idx,Lat_Idx,Elem_Idx,Line_Idx,Zen_Idx)
    integer, intent(in):: Sfc_Level_Idx
      real, intent(in):: Prof_Weight
      integer, intent(in):: Lon_Idx
      integer, intent(in):: Lat_Idx
      integer, intent(in):: Elem_Idx
      integer, intent(in):: Line_Idx
      integer, intent(in):: Zen_Idx
      integer:: Chan_Idx
      real:: Sfc_Ref
      real:: Rad_Ch20_Temp
      real:: Rad_Clear_Ch20_Solar_Rtm
      real:: Bt_Clear_Ch20_Solar_Rtm

      !--------------------------------------------------------------
      ! Solar-Only channels, 1,2,6,7,DNB(44)
      !--------------------------------------------------------------
      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max

         if (Ch(Chan_Idx)%Obs_Type /= SOLAR_OBS_TYPE .and. &
             Ch(Chan_Idx)%Obs_Type /= LUNAR_OBS_TYPE) cycle

         select case (Chan_Idx)

         case (1,2,5,6,7,44)
            if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
               if (allocated(  Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Total_Profile )) then
                  Ch(Chan_Idx)%Trans_Atm_Total(Elem_Idx,Line_Idx) = &
                      Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Ch(Chan_Idx)%Trans_Atm_Total_Profile(Sfc_Level_Idx) +  &
                     (Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Ch(Chan_Idx)%Trans_Atm_Total_Profile(Sfc_Level_Idx+1) -  &
                      Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Ch(Chan_Idx)%Trans_Atm_Total_Profile(Sfc_Level_Idx)) * Prof_Weight
               end if
            end if
         end select

      end do

      !--------------------------------------------------------------
      ! IR-only channels, 20-38 (except 26), 42, 43
      !--------------------------------------------------------------

      !--- upwelling
      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max

         if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle

         if (Ch(Chan_Idx)%Obs_Type /= THERMAL_OBS_TYPE .and. &
             Ch(Chan_Idx)%Obs_Type /= MIXED_OBS_TYPE ) cycle

         call COMPUTE_CHANNEL_ATM_SFC_RAD_BT( &
                Chan_Idx, & ! IN
                Sfc_Level_Idx, & ! IN
                Prof_Weight, & ! IN
                Ch(Chan_Idx)%Sfc_Emiss(Elem_Idx,Line_Idx), & ! IN
                NWP_PIX%Tsfc(Elem_Idx,Line_Idx), & ! IN
                Rtm(Lon_Idx,Lat_Idx)%D(Zen_Idx)%Ch(Chan_Idx)%Rad_Atm_Profile, & ! IN
                Rtm(Lon_Idx,Lat_Idx)%D(Zen_Idx)%Ch(Chan_Idx)%Trans_Atm_Profile, & ! IN
                Ch(Chan_Idx)%Rad_Atm(Elem_Idx,Line_Idx), &        ! OUT
                Ch(Chan_Idx)%Trans_Atm(Elem_Idx,Line_Idx), & ! OUT
                Ch(Chan_Idx)%Rad_Toa_Clear(Elem_Idx,Line_Idx), & ! OUT
                Ch(Chan_Idx)%Bt_Toa_Clear(Elem_Idx,Line_Idx)) ! OUT


      end do


      !--- downwelling (only channel 31/38)
      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max

        if (Chan_Idx == 31 .OR. Chan_Idx == 38) then
          if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle
          call COMPUTE_CHANNEL_ATM_DWN_SFC_RAD( &
               Sfc_Level_Idx, &
               Prof_Weight, &
               Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Ch(Chan_Idx)%Rad_Atm_Dwn_Profile, &
               Ch(Chan_Idx)%Rad_Atm_Dwn_Sfc(Elem_Idx,Line_Idx))
        else
          cycle
        endif

      end do

      !--------------------------------------------------------------
      !-- Add Solar to Ch20 clear variables
      !--------------------------------------------------------------

      if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then

         !--- add in solar component - does not account for glint
         Trans_Atm_Ch20_Solar_Rtm(Elem_Idx,Line_Idx) =  &
                  Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Solar_Profile(Sfc_Level_Idx) + &
                 (Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Solar_Profile(Sfc_Level_Idx+1) - &
                  Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Solar_Profile(Sfc_Level_Idx)) * &
                  Prof_Weight

         Trans_Atm_Ch20_Solar_Total_Rtm(Elem_Idx,Line_Idx) =  &
                  Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Total_Profile(Sfc_Level_Idx) + &
                 (Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Total_Profile(Sfc_Level_Idx+1) - &
                  Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Total_Profile(Sfc_Level_Idx)) * &
                  Prof_Weight

         Rad_Clear_Ch20_Solar_Rtm = ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx)
         Bt_Clear_Ch20_Solar_Rtm = ch(20)%Bt_Toa_Clear(Elem_Idx,Line_Idx)

         if (Geo%Cossolzen(Elem_Idx,Line_Idx) >= 0.0) then

            Sfc_Ref = 1.0 - ch(20)%Sfc_Emiss(Elem_Idx,Line_Idx)

            Rad_Ch20_Temp = ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx) +  &
                   Trans_Atm_Ch20_Solar_Total_Rtm(Elem_Idx,Line_Idx) *  &
                   Sfc_Ref * (Geo%Cossolzen(Elem_Idx,Line_Idx)*Solar_Ch20_Nu / pi)

            Rad_Clear_Ch20_Solar_Rtm = Rad_Ch20_Temp
            Bt_Clear_Ch20_Solar_Rtm = PLANCK_TEMP_FAST(20,Rad_Ch20_Temp)

         end if

         ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx) = Rad_Clear_Ch20_Solar_Rtm
         ch(20)%Bt_Toa_Clear(Elem_Idx,Line_Idx) = Bt_Clear_Ch20_Solar_Rtm

      end if

   end subroutine COMPUTE_CHANNEL_RT

   !-------------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------------
   subroutine COMPUTE_CH20_EMISSIVITY(Elem_Idx,Line_Idx)
      integer, intent(in):: Elem_Idx
      integer, intent(in):: Line_Idx
      real:: Rad_Ch20_Temp
      real:: Ch20_Sfc_Rad

      if ((Sensor%Chan_On_Flag_Default(20) == sym%YES) .and. (Sensor%Chan_On_Flag_Default(31)==sym%YES)) then
         Ch20_Sfc_Rad = ch(20)%Sfc_Emiss(Elem_Idx,Line_Idx) * PLANCK_RAD_FAST(20,NWP_PIX%Tsfc(Elem_Idx,Line_Idx))
         Rad_Ch20_Temp = PLANCK_RAD_FAST(20,ch(31)%Bt_Toa_Clear(Elem_Idx,Line_Idx))
         Ch(20)%Emiss_Rel_11um_Clear(Elem_Idx,Line_Idx) = ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx) / Rad_Ch20_Temp
      end if

      if ((Sensor%Chan_On_Flag_Default(20) == sym%YES) .and. (Sensor%Chan_On_Flag_Default(38)==sym%YES)) then
         Ch20_Sfc_Rad = ch(20)%Sfc_Emiss(Elem_Idx,Line_Idx) * PLANCK_RAD_FAST(20,NWP_PIX%Tsfc(Elem_Idx,Line_Idx))
         Rad_Ch20_Temp = PLANCK_RAD_FAST(20,ch(38)%Bt_Toa_Clear(Elem_Idx,Line_Idx))
         Ch(20)%Emiss_Rel_10_4um_Clear(Elem_Idx,Line_Idx) = ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx) / Rad_Ch20_Temp
      end if

   end subroutine COMPUTE_CH20_EMISSIVITY

   !-------------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------------
   subroutine COMPUTE_TROPOPAUSE_EMISSIVITIES(Elem_Idx,Line_Idx,Lon_Idx,Lat_Idx,Zen_Idx)
      integer, intent(in):: Elem_Idx
      integer, intent(in):: Line_Idx
      integer, intent(in):: Lon_Idx
      integer, intent(in):: Lat_Idx
      integer, intent(in):: Zen_Idx
      integer:: Chan_Idx
      integer:: Lev_Bnd
      integer:: dim1
      integer:: dim2

      dim1 = Image%Number_Of_Elements
      dim2 = Image%Number_Of_Lines_Per_Segment

      Lev_Bnd = Rtm(Lon_Idx,Lat_Idx)%Tropo_Level

      !--- check for missing tropopause level
      if (Rtm(Lon_Idx,Lat_Idx)%Tropo_Level == 0) then
         Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
         return
      end if

      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max

         select case (Chan_Idx)


         case(20,27,29,31,32,33,37,38)
            if (Sensor%Chan_On_Flag_Default(Chan_Idx)==sym%YES) then
               ch(Chan_Idx)%Emiss_Tropo(Elem_Idx,Line_Idx) =  &
                        EMISSIVITY(ch(Chan_Idx)%Rad_Toa(Elem_Idx,Line_Idx),  &
                        ch(Chan_Idx)%Rad_Toa_Clear(Elem_Idx,Line_Idx),  &
                        Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_BB_Cloud_Profile(Lev_Bnd))

            end if
         end select
      end do

   end subroutine COMPUTE_TROPOPAUSE_EMISSIVITIES

   !-------------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------------
   subroutine COMPUTE_BETA_RATIOES(Elem_Idx,Line_Idx)
      integer, intent(in):: Elem_Idx
      integer, intent(in):: Line_Idx

      !--- compute 11 and 12 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(32) == sym%YES) then

         Beta_11um_12um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(32)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx))
      end if

      !--- compute 10 and 12 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(38) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(32) == sym%YES) then

         Beta_104um_12um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(32)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(38)%Emiss_Tropo(Elem_Idx,Line_Idx))
      end if

      !--- compute 11 and 8.5 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
          Sensor%Chan_On_Flag_Default(29) == sym%YES) then

         Beta_11um_85um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(29)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx))
      endif

      !--- compute 10 and 8.5 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(38) == sym%YES .and. &
          Sensor%Chan_On_Flag_Default(29) == sym%YES) then

         Beta_104um_85um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(29)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(38)%Emiss_Tropo(Elem_Idx,Line_Idx))
      endif

      !--- compute 11 and 6.7 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(27) == sym%YES) then

         Beta_11um_67um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(27)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx))
      end if

      !--- compute 10 and 6.7 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(38) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(27) == sym%YES) then

         Beta_104um_67um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(27)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(38)%Emiss_Tropo(Elem_Idx,Line_Idx))
      end if

      !--- compute 11 and 13.3 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
          Sensor%Chan_On_Flag_Default(33) == sym%YES) then

         Beta_11um_133um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(33)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx))
      endif

      !--- compute 10 and 13.3 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(38) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(33) == sym%YES) then

         Beta_104um_133um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(33)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(38)%Emiss_Tropo(Elem_Idx,Line_Idx))
      endif

      !--- compute 11 and 10 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
          Sensor%Chan_On_Flag_Default(38) == sym%YES) then

         Beta_11um_104um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(38)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx))
      endif

      !--- compute 10 and 11 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
          Sensor%Chan_On_Flag_Default(38) == sym%YES) then

         Beta_104um_11um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(38)%Emiss_Tropo(Elem_Idx,Line_Idx))
      endif

   end subroutine COMPUTE_BETA_RATIOES

   !====================================================================
   ! FUNCTION Name: BETA_RATIO
   !
   ! Function:
   !  Computes the beta ratio for two Emissivities.
   !
   ! Input:  Emiss_top - emissivity in the numerator
   !         Emiss_bot - emissivity in the denominator
   !
   ! Output: Beta - the beta value from the two emissivities
   !
   !====================================================================
   function BETA_RATIO(Emiss_top, Emiss_bot) result(beta)
      real(kind=real4), intent(in) :: Emiss_top
      real(kind=real4), intent(in) :: Emiss_bot
      real(kind=real4) :: beta

      beta = Missing_Value_Real4

      if (Emiss_top > 0.0 .and. Emiss_top < 1.0 .and. &
            Emiss_bot > 0.0 .and. Emiss_bot < 1.0) then

         beta = alog(1.0 - Emiss_top)/alog(1.0 - Emiss_bot)
      end if

      return

   end function BETA_RATIO

   !====================================================================
   ! Function Name: EMISSIVITY
   !
   ! Function:
   !  Computes the  effective emissivity
   !
   ! Input:  Rad_Toa - channel radiance at top of atmosphere(toa)
   !         Rad_Clear_Tau - channel radiance at toa for clear skies
   !         Rad_Cloud_BB_Toa - channel radiance at TOA if cloud were a Black-Body
   !
   ! Output: Emiss - the effective cloud emissivity
   !
   !====================================================================
   function EMISSIVITY(Radiance_Toa, Radiance_Clear_Toa, Radiance_Cloud_BB_Toa) result(Emiss)
      real(kind=real4), intent(in) :: Radiance_Toa
      real(kind=real4), intent(in) :: Radiance_Clear_Toa
      real(kind=real4), intent(in) :: Radiance_Cloud_BB_Toa
      real(kind=real4) :: Emiss

      Emiss = Missing_Value_Real4

      if (Radiance_Cloud_BB_Toa .NER. Radiance_Clear_Toa) then
          Emiss = (Radiance_Toa - Radiance_Clear_Toa) / &
            (Radiance_Cloud_BB_Toa - Radiance_Clear_Toa)
       end if

      return

   end function EMISSIVITY


   !===============================================================================
   !
   !===============================================================================
   subroutine COPY_LOCAL_RTM_TO_GLOBAL_RTM_STRUCTURE(Lon_Idx,Lat_Idx,Zen_Idx)

      integer, intent(in):: Lon_Idx
      integer, intent(in):: Lat_Idx
      integer, intent(in):: Zen_Idx
      integer:: Chan_Idx

      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max
         if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Total_Profile = Trans_Atm_Total_Prof(:,Chan_Idx)
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Solar_Profile = Trans_Atm_Solar_Prof(:,Chan_Idx)
         !Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Profile = Trans_Atm_Prof(:,Chan_Idx)
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Profile = Rad_Atm_Prof(:,Chan_Idx)
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_BB_Cloud_Profile = Rad_BB_Cloud_Prof(:,Chan_Idx)
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)% is_set = .true.
      end do

      Chan_Idx = 31
      if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Dwn_Profile = Rad_Atm_Dwn_Prof(:,Chan_Idx)
      end if

      Chan_Idx = 38
      if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Dwn_Profile = Rad_Atm_Dwn_Prof(:,Chan_Idx)
      end if
     Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)% is_set = .true.

   end subroutine COPY_LOCAL_RTM_TO_GLOBAL_RTM_STRUCTURE


   !===============================================================================
   ! compute an inversion profile
   ! a level is considered in a inversion if it is colder than the level above it
   !===============================================================================
   subroutine INVERSION_PROFILE(T_Prof,Inver_Prof)

      real, dimension(:), intent(in):: T_Prof
      integer (kind=int1), dimension(:), intent(out):: Inver_Prof

      integer:: N_Levels
      integer:: Lev_Idx

      Inver_Prof = 0

      n_levels = size (T_prof)

      do Lev_Idx = 2, N_Levels
        if (T_Prof(Lev_idx) < T_Prof(Lev_Idx-1)) Inver_Prof(Lev_Idx) = 1
      enddo

   end subroutine INVERSION_PROFILE


!--- end of module
!
end module RT_UTILITIES_MOD
