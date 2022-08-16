! $Id: pixel_routines_mod.f90 4044 2020-11-10 21:52:49Z heidinger $
!-------------------------------------------------------------------------------------- 
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: pixel_routines.f90 (src)
!       PIXEL_ROUTINES (program)
!
! PURPOSE: this module houses routines for computing some needed pixel-level arrays
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
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
! Public routines used in this MODULE:
! COMPUTE_PIXEL_ARRAYS - compute some commonly used arrays
! COMPUTE_TSFC - derive pixel-level surface temperature 
! NORMALIZE_REFLECTANCES - divide reflectances by cosine solar zenith angle
! CH20_PSEUDO_REFLECTANCE - compute the channel 20 reflectance
! SPECTRAL_CORRECT_NDVI - apply a spectral correct to Ndvi to look like NOAA14
! ASSIGN_CLEAR_SKY_QUALITY_FLAGS - assign quality flags to clear-sky products
! CONVERT_TIME - compute a time in hours based on millisecond time in leveL1b
! COMPUTE_GLINT - derive a glint mask
!
! DETERMINE_LEVEL1B_COMPRESSION
!
!--------------------------------------------------------------------------------------
module PIXEL_ROUTINES_MOD

 use CONSTANTS_MOD, only: &
  missing_value_real4 &
    , sym , int1, int2, exe_prompt  &
    , real4, missing_value_int1, int4, mixed_obs_type, NChan_Clavrx, SOLAR_OBS_TYPE &
    , THERMAL_OBS_TYPE, PI, terminator_reflectance_sol_zen_thresh, DTOR
 
 use ALGORITHM_CONSTANTS_MOD,only: glint_zen_thresh, Ref_Sfc_White_Sky_Water
 
 use PIXEL_COMMON_MOD, only: &
    sensor , image, ch, nwp_pix &
    , cldmask, geo, sfc, nav, acha, solar_rtm &
    , ch1_counts &
    , ch3a_on_avhrr &
    , bad_scan_flag, bad_pixel_mask &
    , pixel_local_time_hours &
    , ref_ch1_dark_composite &
    , Aot_Qf &
    , dcomp_mode &
    , caliop_flag &
    , therm_cal_1b &
    , use_aux_flag &
    , gap_pixel_mask &
    , segment_valid_fraction &
    , btd_ch20_ch38, btd_ch20_ch32, btd_ch20_ch31 &
    , tsfc_retrieved, tsfc_qf, trad_retrieved &
    , abi_use_104um_flag, ancil_data_dir &
    , rsr_qf, ndvi_qf, aerosol_mode &
    , temp_pix_array_1, cld_phase_aux &
    , cld_phase_uncertainty &
    , zc_opaque_cloud, tc_opaque_cloud &
    , beta_104um_12um_tropo_rtm, beta_11um_12um_tropo_rtm, beta_11um_133um_tropo_rtm &
    , beta_104um_133um_tropo_rtm &
    , rsr &
    , cld_type_aux &
    , number_of_temporary_files &
    , temporary_file_name, temporary_data_dir &
    , dcomp_success_fraction , dcomp_quality_flag &
    , Btd_Ch31_Ch32, Btd_Ch38_Ch32 &
    , ndsi_sfc, nddi_toa, ndsi_toa, Ndvi_Sfc, Ndvi_Sfc_White_Sky,Ndvi_Toa &
    , nonconfident_cloud_mask_fraction
 
 use NUMERICAL_ROUTINES_MOD,only:
 
 use NWP_COMMON_MOD,only: &
  nwp
 
 use PLANCK_MOD, only: &
  planck_temp_fast &
  , planck_rad_fast
 
 use LAND_SFC_PROPERTIES_MOD, only: &
  Land_grid_description &
  , read_land_sfc_hdf
  
 use FILE_TOOLS,only: getlun
 
 use SURFACE_PROPERTIES_MOD, only:
 
 use CALIBRATION_CONSTANTS_MOD,only:
 
 use ECM2_CLOUD_MASK_CLAVRX_BRIDGE, only: COMPUTE_TYPE_FROM_PHASE
 
 use CLAVRX_MESSAGE_MOD, only: MESG, verb_lev
 
!use RT_UTILITIES_MOD, only: COMPUTE_CLEAR_SKY_SCATTER

 implicit none
 private
 public:: COMPUTE_PIXEL_ARRAYS, &
          SURFACE_REMOTE_SENSING,  &
          NORMALIZE_REFLECTANCES,  &
          CH20_PSEUDO_REFLECTANCE,  &
          ASSIGN_CLEAR_SKY_QUALITY_FLAGS, &
          EXPAND_SPACE_MASK_FOR_USER_LIMITS, &
          SET_SOLAR_CONTAMINATION_MASK, &
          SET_BAD_PIXEL_MASK, &
          QUALITY_CONTROL_ANCILLARY_DATA,   &
          READ_MODIS_WHITE_SKY_ALBEDO,      &
          COMPUTE_GLINT,                    &
          SET_CHAN_ON_FLAG,                 &
          DETERMINE_LEVEL1B_COMPRESSION, &
          TERM_REFL_NORM, &
          MERGE_NWP_HIRES_ZSFC, &
          ADJACENT_PIXEL_CLOUD_MASK, &
          COMPUTE_ACHA_PERFORMANCE_METRICS, &
          COMPUTE_DCOMP_PERFORMANCE_METRICS, &
          COMPUTE_CLOUD_MASK_PERFORMANCE_METRICS, &
          MODIFY_LAND_CLASS_WITH_NDVI, &
          DESERT_MASK_FOR_CLOUD_DETECTION, &
          CITY_MASK_FOR_CLOUD_DETECTION, &
          VIIRS_TO_MODIS, &
          MODIFY_AUX_CLOUD_TYPE

  private:: REMOTE_SENSING_REFLECTANCE, &
            NORMALIZED_DifFERENCE_VEGETATION_INDEX, &
            NORMALIZED_DifFERENCE_SNOW_INDEX, &
            NORMALIZED_DifFERENCE_DESERT_INDEX, &
            COMPUTE_TSFC

  contains

   !----------------------------------------------------------------------
   ! set Chan_On_Flag for each to account for Ch3a/b switching on avhrr
   !
   ! this logic allows the default values to also be used to turn off
   ! channels
   !
   !  called by process_clavrx inside segment loop
   !    HISTORY: 2014/12/29 (AW); removed unused ch3a_on_avhrr for modis and goes
   !----------------------------------------------------------------------
   subroutine SET_CHAN_ON_FLAG(Chan_On_Flag_Default, Chan_On_Flag_Per_Line, is_last_segment)

      integer(kind=int1), dimension(:), intent(in):: Chan_On_Flag_Default
      integer(kind=int1), dimension(:,:), intent(out):: Chan_On_Flag_Per_Line
      logical, intent(in):: is_last_segment
      integer:: Number_of_Elements
      integer:: Number_of_Lines
      integer:: Line_Idx
     
      logical :: dcomp_first_valid_line_avhrr_set = .false. 

      Number_of_Elements = Image%Number_Of_Elements
      Number_of_Lines = Image%Number_Of_Lines_Read_This_Segment
    
      line_loop: do Line_Idx = 1, Number_Of_Lines
          ! - change dcomp _mode according ch3a_on flag
          ! - this change of dcomp_mode is only possible once for one file
          ! - First daytime line determines dcomp mode for whole file
          ! - AW 02/13/2017
         if ( index(Sensor%Sensor_Name,'AVHRR') > 0 &
            .and. .not. dcomp_first_valid_line_avhrr_set  &
            .and. Geo%Solzen(1,line_idx) .lt. 82 ) then
            if (Ch3a_On_Avhrr(Line_Idx) == sym%YES) then
                dcomp_first_valid_line_avhrr_set = .true.
                if (dcomp_mode .ne. 0  ) dcomp_mode = 1
                
            end if
            if (Ch3a_On_Avhrr(Line_Idx) == sym%NO) then
               dcomp_first_valid_line_avhrr_set = .true.
               if (dcomp_mode .ne. 0  )  dcomp_mode = 3
            
            end if
         end if   
         
         ! - for all sensors : set chan_on_flag ( dimension [n_chn, n_lines] to default ) 
         Chan_On_Flag_Per_Line(:,Line_Idx) = Chan_On_Flag_Default   
         
         ! two exceptions
         if (trim(Sensor%Platform_Name)=='AQUA' .and. Chan_On_Flag_Default(6) == sym%YES ) then
            if (minval(ch(6)%DQF(:,Line_Idx)) >= 15) then
                 Chan_On_Flag_Per_Line(6,Line_Idx) = sym%NO 
            end if  
         end if
         
         if (index(Sensor%Sensor_Name,'AVHRR') > 0) then
            if (Ch3a_On_Avhrr(Line_Idx) == sym%YES) then
               Chan_On_Flag_Per_Line(6,Line_Idx) = Chan_On_Flag_Default(6)   
               Chan_On_Flag_Per_Line(20,Line_Idx) = sym%NO   
            end if
            if (Ch3a_On_Avhrr(Line_Idx) == sym%NO) then
               Chan_On_Flag_Per_Line(6,Line_Idx) = sym%NO   
               Chan_On_Flag_Per_Line(20,Line_Idx) = Chan_On_Flag_Default(20)   
            end if
         endif

      end do line_loop
      
      
      if ( is_last_segment) dcomp_first_valid_line_avhrr_set = .false. 
      
   end subroutine SET_CHAN_ON_FLAG



!======================================================================
! Modify the space mask based the limits on lat, lon, satzen and solzen 
! Space mask was initially determined in the Level-1b Navigation
!======================================================================
subroutine EXPAND_SPACE_MASK_FOR_USER_LIMITS(Seg_Idx, Space_Mask)

   use CALIOP_COLLOCATION_MOD, only: CALIOP_COLLOCATION

   integer(kind=int4), intent(in):: Seg_Idx
   logical, dimension(:,:), intent(inout):: Space_Mask

   !--- check for latitudinal bounds
   where(Nav%Lat_1b == Missing_Value_Real4 .or. Nav%Lon_1b == Missing_Value_Real4)
        Space_Mask = .true.
   end where

   where(isnan(Nav%Lat_1b) .or. isnan(Nav%Lon_1b))
        Space_Mask = .true.
   end where

   !--- check if subset processing is on
   if (Nav%Limit_Flag == sym%YES) then

      where(Nav%Lat_1b < Nav%Lat_South_Limit .or. Nav%Lat_1b > Nav%Lat_North_Limit)
         Space_Mask = .true.
      end where

      !--- check for longitudinal bounds including the dateline condition
      if ( Nav%Lon_West_Limit > Nav%Lon_East_Limit) then
         where((Nav%Lon_1b < Nav%Lon_West_Limit .and. Nav%Lon_1b > 0.0) .or. &
             (Nav%Lon_1b > Nav%Lon_East_Limit .and. Nav%Lon_1b < 0.0))
            Space_Mask =.true.
         end where
      else
         where(Nav%Lon_1b < Nav%Lon_West_Limit .or. Nav%Lon_1b > Nav%Lon_East_Limit)
            Space_Mask =.true.
         end where
      end if

      !--- Satzen limit
      where (Geo%Satzen > Geo%Satzen_Max_Limit .or. Geo%Satzen < Geo%Satzen_Min_Limit)
         Space_Mask = .true.
      end where

      !--- Solzen limit
      where (Geo%Solzen < Geo%Solzen_Min_Limit .or. Geo%Solzen > Geo%Solzen_Max_Limit .or. Geo%Satzen == Missing_Value_Real4)
         Space_Mask = .true.
      end where

   endif

   !--- CALIOP collocation
   if (Caliop_Flag) then
      call CALIOP_COLLOCATION(Seg_Idx)
   end if

   !--- test if any valid data, if not, print a warning 
   if (ALL(Space_Mask)) then
      print *, EXE_PROMPT, "WARNING: All Data in Segment are Classified as Space "
   endif

end subroutine EXPAND_SPACE_MASK_FOR_USER_LIMITS
 
!======================================================================
! Check for solar contamination of what should be nighttime data
! 
! this is common in AVHRR and GOES.  Pixels with solar contamination
!
! are treated as bad pixel in the bad_pixel_mask
!
! Note that the Ch1_Counts used here are assumed have the dark count subtracted
! from them.
!
!======================================================================
subroutine SET_SOLAR_CONTAMINATION_MASK(Solar_Contamination_Mask) 

   integer(kind=int1), dimension(:,:), intent(out):: Solar_Contamination_Mask
   integer:: Number_of_Elements
   integer:: Number_of_Lines
   integer:: Elem_Idx
   integer:: Line_Idx
   integer, parameter:: Solar_Contamination_Thresh_AVHRR = 2
   integer, parameter:: Solar_Contamination_Thresh_GEO = 10

   Number_of_Elements = Image%Number_Of_Elements
   Number_of_Lines = Image%Number_Of_Lines_Read_This_Segment

   !--- initialize 
   Solar_Contamination_Mask(:,1:Number_Of_Lines) = sym%NO

   !---  loop through lines and elements
   line_loop: do Line_Idx = 1, Number_of_Lines
      element_loop: do Elem_Idx = 1, Number_of_Elements

        !--- check for solar contamination of nighttime data in AVHRR
        if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then

          if (index(Sensor%Sensor_Name,'AVHRR') > 0) then

            if ((Geo%Solzen(Elem_Idx,Line_Idx) > 90.0) .and. (Geo%Scatangle(Elem_Idx,Line_Idx) < 60.0)) then
              if (therm_cal_1b == sym%NO) then
                if (Ch1_Counts(Elem_Idx,Line_Idx) > Solar_Contamination_Thresh_AVHRR) then 
                   Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%YES
                endif
              else
                if (Ch1_Counts(Elem_Idx,Line_Idx) > Solar_Contamination_Thresh_AVHRR) then 
                   Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%YES
                endif 
              endif

            endif
  
          endif

          !--- check for solar contamination of nighttime data in GOES
          if (index(Sensor%Sensor_Name,'GOES') > 0) then
             if ((Geo%Solzen(Elem_Idx,Line_Idx) > 90.0) .and. (Geo%Scatangle(Elem_Idx,Line_Idx) < 60.0)) then
                if (Ch1_Counts(Elem_Idx,Line_Idx) > Solar_Contamination_Thresh_GEO) then
                   Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%YES
                endif 
             endif
          endif

        endif

        !--- check for solar contamination of nighttime data in GOES
        if (index(Sensor%Sensor_Name,'GOES') > 0) then
          if ((Geo%Solzen(Elem_Idx,Line_Idx) > 90.0) .and.  (Geo%Scatangle(Elem_Idx,Line_Idx) < 180.0)) then
            if (Ch1_Counts(Elem_Idx,Line_Idx)  > Solar_Contamination_Thresh_GEO) then
              Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%YES
            endif
          endif
        endif


        !until a more robust fix can be found, for now we will set all night
        ! pixels to have Solar_Contamination_Mask = sym%YES due to telescope
        ! contamination in 3.9um channel. Exists for FY-2C/D/E/G - WCS3
        
        if (index(Sensor%Sensor_Name,'FY2-IMAGER') > 0) then
          if ((Geo%Solzen(Elem_Idx,Line_Idx) > 90.0)) then
              Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%YES
          endif
        endif

      end do element_loop
   end do line_loop

end subroutine SET_SOLAR_CONTAMINATION_MASK

!======================================================================
! Check for bad pixels
!
! Apply tests to detect pixels that should not be processed.
!
!
! Bad_Pixel_Mask is meant to single mask that captures all reasons
! why a pixel should be skipped. This includes
!
! Space_Mask
! Solar_Contamination_Mask
! Missing NWP
! any one of channel tests
!
! Also, if many pixels in a line are bad, the whole line is set to bad
!
! Output:  Bad_Pixel_Mask 
! Input: taken from pixel common
!  
!======================================================================
subroutine SET_BAD_PIXEL_MASK(Bad_Pixel_Mask,ABI_Use_104um_Flag)

   integer(kind=int1), dimension(:,:), intent(out):: Bad_Pixel_Mask
   integer:: Number_of_Elements
   integer:: Number_of_Lines
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: Number_Bad_Pixels
   integer:: Number_Bad_Pixels_Thresh
   integer:: Lon_Nwp_Idx
   integer:: Lat_Nwp_Idx
   logical, intent(in):: ABI_Use_104um_Flag

   Number_of_Elements = Image%Number_Of_Elements
   Number_of_Lines = Image%Number_Of_Lines_Read_This_Segment

   Number_Bad_Pixels_Thresh = 0.9 * Image%Number_Of_Elements

!----------------------------------------------------------------------
!--- assign bad pixel mask based on scanline fatal flag
!----------------------------------------------------------------------
   !---- this is needed to ensure extra lines (beyond Num_Scans_Read) are bad
   Bad_Pixel_Mask = sym%YES

   line_loop: do Line_Idx = 1, Number_of_Lines

      !--- initialize
      Bad_Pixel_Mask(:,Line_Idx) = sym%NO

      !--- check for a bad scan
      if (Bad_Scan_Flag(Line_Idx) == sym%YES) then

        Bad_Pixel_Mask(:,Line_Idx) = sym%YES

      else

      !--- if not a bad scan, check pixels on this scan
      element_loop: do Elem_Idx = 1, Number_of_Elements
   

        !--- set space to bad
        if (Geo%Space_Mask(Elem_Idx,Line_Idx) ) then
            Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
            !print *, "Bad Space"
        endif

        !--- NaN checks on geolocation and geometry
        if (isnan(Geo%Satzen(Elem_Idx,Line_Idx)) .or.  &
            isnan(Geo%Solzen(Elem_Idx,Line_Idx))) then
            Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
            !print *, "Bad Nan"
        endif

        !--- Bad Relazimuth
        if (Geo%Relaz(Elem_Idx,Line_Idx) == Missing_Value_Real4) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
           !print *, "Bad Relaz"
        endif

        if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(32) == sym%YES)) then

          !--- CALL any scan with a ridiculous pixel as bad 
          !--- this is attempt data like NOAA-16 2004 023 where
          !--- large fractions of scans are bad but not flagged as so
!         if (abs(ch(31)%Bt_Toa(Elem_Idx,Line_Idx) - ch(32)%Bt_Toa(Elem_Idx,Line_Idx)) > 20.0) then
!             Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
!         endif

        endif

        !--- Set based on 10 or 11 um observations.

        if (ABI_Use_104um_Flag) then

          !--- Use 10.4 um observations.
          if (Sensor%Chan_On_Flag_Default(38) == sym%YES) then

            if (ch(38)%Bt_Toa(Elem_Idx,Line_Idx) < 150.0) then
              Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
            endif

            if (ch(38)%Bt_Toa(Elem_Idx,Line_Idx) > 350.0) then
              Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
            endif

            if (isnan(ch(38)%Bt_Toa(Elem_Idx,Line_Idx))) then
              Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
            endif

          endif

        else

          !--- Use 11 um observations.
          if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then

            if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < 150.0) then
              Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
              !print *, "Cold 31"
            endif

            if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) > 350.0) then
              Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
              !print *, "Hot 31"
            endif

            if (isnan(ch(31)%Bt_Toa(Elem_Idx,Line_Idx))) then
              Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
              !print *, "Nan 31"
            endif

          endif

        endif

        !--- AVHRR/3 Wedge Filter
        !--- for AVHRR/3, consider region where Ch3a is on but seeing night
        !--- as bad data.  In these regions there is no ch3b
!       if ((AVHRR_Flag == sym%YES) .and. &
!           ((Sc_Id_WMO <= 5) .or. (Sc_Id_WMO >= 206 .and. Sc_Id_WMO <= 223)) .and. &   !AVHRR/3 only
!           (Sensor%Chan_On_Flag_Default(6) == sym%YES) .and. &
!           (Sensor%Chan_On_Flag_Default(20) == sym%YES) .and. &
!           (Solzen(Elem_Idx,Line_Idx) > 90.0)) then
!           if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
!            if (ch(20)%Bt_Toa(Elem_Idx,Line_Idx) < 0.0) then
!               Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
!            endif
!           endif
!       endif
        

        !--- check for solar zenith angle limits
        if ((Geo%Solzen(Elem_Idx,Line_Idx) < Geo%Solzen_Min_Limit) .or. &
            (Geo%Solzen(Elem_Idx,Line_Idx) > Geo%Solzen_Max_Limit)) then
             Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
             !print *, "Bad Solzen"
        endif

        !--- NWP
        if (NWP_PIX%Nwp_Opt /= 0) then
            Lon_Nwp_Idx = NWP_PIX%i_Nwp(Elem_Idx,Line_Idx)
            Lat_Nwp_Idx = NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)
            if (Lon_Nwp_Idx < 1 .or. Lat_Nwp_Idx < 1) then
                 Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
                 !print *, "Bad NWP"
            else
                 if (NWP%Bad_Nwp_Mask(Lon_Nwp_Idx, Lat_Nwp_Idx) == sym%YES) Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
            endif
        endif

      end do element_loop

     endif

     !------ if 90% of pixels on a line are bad, mark the whole scan line as bad (if not already)
     if (Bad_Scan_Flag(Line_Idx) == Missing_Value_Int1) then
         Number_Bad_Pixels = sum(Bad_Pixel_Mask(:,Line_Idx),Bad_Pixel_Mask(:,Line_Idx)==sym%YES)
         if (Number_Bad_Pixels > Number_Bad_Pixels_Thresh) then
           Bad_Scan_Flag(Line_Idx) = sym%YES
         else
           Bad_Scan_Flag(Line_Idx) = sym%NO
         endif
     endif


     !---- consider any scanline with any solar contamination as a bad line
!    if (maxval(Solar_Contamination_Mask(:,Line_Idx)) == sym%YES) then
!        Bad_Scan_Flag(Line_Idx) = sym%YES
!        Bad_Pixel_Mask(:,Line_Idx) = sym%YES
!    endif
     

   end do line_loop

  !-----------------------------------------------------------------------------------
  ! if the IDPS cloud mask is to be used for product generation, make sure that
  ! pixels within the gaps are considered bad
  !-----------------------------------------------------------------------------------
  if (Use_Aux_Flag == sym%USE_AUX .and. trim(Sensor%Sensor_Name) == 'VIIRS') then
      where(Gap_Pixel_Mask == sym%YES)
         Bad_Pixel_Mask = sym%YES
      endwhere
  endif

  !---------------------------------------------------------------------------------------
  ! Compute the fraction of the segment covered by valid data
  !---------------------------------------------------------------------------------------
  Segment_Valid_Fraction = 1.0 - sum(float(Bad_Pixel_Mask(:,1:Number_of_Lines))) /  &
                                float(Number_of_Elements * Number_of_Lines)

end subroutine SET_BAD_PIXEL_MASK
!--------------------------------------------------------------------------
!QUALITY_CONTROL_ANCILLARY_DATA
!
! Apply some checks on ancillary data.  Call pixels bad when checks fail
!
! Note, Bad_Pixel_Mask is modified but not created here. Do not initialize
!--------------------------------------------------------------------------
subroutine QUALITY_CONTROL_ANCILLARY_DATA (Bad_Pixel_Mask)
   integer(kind=int1), dimension(:,:), intent(inout):: Bad_Pixel_Mask
   integer:: Number_of_Elements
   integer:: Number_of_Lines
   integer:: Elem_Idx
   integer:: Line_Idx

   Number_of_Elements = Image%Number_Of_Elements
   Number_of_Lines = Image%Number_Of_Lines_Read_This_Segment

   do Line_Idx = 1, Number_of_Lines

      do Elem_Idx = 1, Number_Of_Elements

        !--- invalid sfc type observations
        if (Sfc%Sfc_Type(Elem_Idx,Line_Idx) < 0 .or. Sfc%Sfc_Type(Elem_Idx,Line_Idx) > 15) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
        endif

      enddo

   enddo

end subroutine QUALITY_CONTROL_ANCILLARY_DATA



!--------------------------------------------------------------------------
! COMPUTE_PIXEL_ARRAYS
!
! compute quantities needed for other routines, (ie cloud mask)
!
! input - none
!
! output - only into shared memory
!   seczen - secant of zenith angle
!   Btd_Ch31_Ch32 - brightness temperature difference between Bt_Ch31 and Bt_Ch32 
!   Btd_Ch20_Ch31 - brightness temperature difference between Bt_Ch20 and Bt_Ch31 
!   Ref_ratio_Ch6_Ch1 - ratio of Ref_Ch6 / Ref_Ch1
!   Ref_ratio_Ch20_Ch1 - ratio of Ref_Ch20 / Ref_Ch1
!--------------------------------------------------------------------------
 subroutine COMPUTE_PIXEL_ARRAYS(j1,nj)

   integer, intent(in):: j1,nj
   integer(kind=int4):: j2

   j2 = j1 + nj - 1

   Geo%Coszen(:,j1:j2) = cos(Geo%Satzen(:,j1:j2)*dtor)

   Geo%Seczen(:,j1:j2) = 1.0 / Geo%Coszen(:,j1:j2)

   Geo%Cossolzen(:,j1:j2) = cos(Geo%Solzen(:,j1:j2)*dtor)

   Geo%Airmass(:,j1:j2) = 1.0
   where(Geo%Solzen(:,j1:j2) /= 0.0 .and. Geo%Coszen(:,j1:j2) /= 0.0) 
    Geo%Airmass(:,j1:j2) = 1.0/Geo%Cossolzen(:,j1:j2) + 1.0/Geo%Coszen(:,j1:j2)
   endwhere

!--- other useful arrays 

   !--- channel 31 and channel 32 brightness temperature difference
   if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
       (Sensor%Chan_On_Flag_Default(32) == sym%YES)) then
        Btd_Ch31_Ch32(:,j1:j2) = ch(31)%Bt_Toa(:,j1:j2) - ch(32)%Bt_Toa(:,j1:j2)
   endif

   !--- channel 38 and channel 32 brightness temperature difference
   if ((Sensor%Chan_On_Flag_Default(38) == sym%YES) .and. &
       (Sensor%Chan_On_Flag_Default(32) == sym%YES)) then
        Btd_Ch38_Ch32(:,j1:j2) = ch(38)%Bt_Toa(:,j1:j2) - ch(32)%Bt_Toa(:,j1:j2)
   endif

   !--- channel 20 and channel 31 brightness temperature difference
    if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
     Btd_Ch20_Ch31 = ch(20)%Bt_Toa - ch(31)%Bt_Toa
     where(ch(20)%Bt_Toa == Missing_Value_Real4 .or. ch(31)%Bt_Toa == Missing_Value_Real4)
       Btd_Ch20_Ch31 = Missing_Value_Real4
     endwhere
    endif

   !--- channel 20 and channel 38 brightness temperature difference
    if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
     Btd_Ch20_Ch38 = ch(20)%Bt_Toa - ch(38)%Bt_Toa
     where(ch(20)%Bt_Toa == Missing_Value_Real4 .or. ch(38)%Bt_Toa == Missing_Value_Real4)
       Btd_Ch20_Ch38 = Missing_Value_Real4
     endwhere
    endif

   !--- channel 20 and channel 32 brightness temperature difference
    if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
     Btd_Ch20_Ch32 = ch(20)%Bt_Toa - ch(32)%Bt_Toa
     where(ch(20)%Bt_Toa == Missing_Value_Real4 .or. ch(32)%Bt_Toa == Missing_Value_Real4)
       Btd_Ch20_Ch32 = Missing_Value_Real4
     endwhere
    endif

 end subroutine COMPUTE_PIXEL_ARRAYS

!------------------------------------------------------------------
! Compute a surface temperature by using the observed 
! 11 micron radiance, the rtm calculations and the surface emissivity
!
! July 2009 - made AVHRR/1 compliant
!
! Author: Andrew Heidinger
!
! Sept 2021 Eva Borbas:  NOAA dual channel Enterprise Algorithm is added 
!                        for Himawari-AHI and SNPP-VIIRS
!             
!------------------------------------------------------------------
 subroutine COMPUTE_TSFC(jmin,jmax,ABI_Use_104um_Flag)

  integer, intent(in):: jmin
  integer, intent(in):: jmax
  logical, intent(in):: ABI_Use_104um_Flag

  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: Coef_Idx,RegCoef_Idx
  integer:: xnwp
  integer:: ynwp
  integer:: tpwcls,angcls,dnflag 
  integer(kind=int4):: Coef_Lun  
  integer:: ios, ERR
  
  INTEGER, PARAMETER :: ncoef=6
  INTEGER, PARAMETER :: nbtpwcls=3
  INTEGER, PARAMETER :: nbangcls=4
  INTEGER, PARAMETER :: nbdn=2
  INTEGER, PARAMETER :: leng=ncoef*nbtpwcls*nbangcls*nbdn
      
  REAL :: Coef(ncoef)
  REAL :: Regcoef(leng)
                  
  CHARACTER(LEN=200) :: Coef_Fn
  
  real:: Rad11
  real:: Rad11_Atm
  real:: Trans11_Atm
  real:: Rad11_Atm_Dwn_Sfc
  real:: Emiss_Sfc11, Emiss_Sfc12, Emiss_mean, Emiss_Diff
  real:: Rad11_Sfc
  real:: B11_Sfc, B12_Sfc, Bdiff
  real:: Sec_Sat_Zen, Sat_Zen
  real:: Sol_Zen
  real:: TPW
    
  logical :: tsfc_onechannel = .true.
  logical :: first_segment = .true.

  !------------------------------------------------------------------
  ! Loop over pixels and derive surface temp
  !------------------------------------------------------------------

  
  if (Sensor%Platform_Name == 'SNPP' .and. trim(Sensor%Sensor_Name) == 'VIIRS' ) then
            coef_fn=trim(Ancil_Data_Dir)//'/static/sfc_data/cx_lstrc_jpss_0_viirs.bin'
            tsfc_onechannel = .false.
            tsfc_onechannel = .true. ! VIIRS coefs are still off AW 13 Oct.2021
  end if
   
  if (Sensor%Platform_Name == 'HIM8' .and. trim(Sensor%Sensor_Name) == 'AHI ') then
            coef_fn=trim(Ancil_Data_Dir)//'/static/sfc_data/cx_lstrc_himawari_8_ahi.bin'
            tsfc_onechannel = .false.
  end if
    
  if ( first_segment)  then
    if (tsfc_onechannel)  call MESG ('tsfc_onechannel TRUE', level = verb_lev %DEFAULT)
    if ( .not. tsfc_onechannel)  call MESG ('tsfc_onechannel FALSE', level = verb_lev %DEFAULT)
  end if
  !--- initialize
  Tsfc_Retrieved = Missing_Value_Real4
  Trad_Retrieved = Missing_Value_Real4
  Tsfc_Qf = 0

  if (ABI_Use_104um_Flag) then
    !--- if no ch38, abort
    if (Sensor%Chan_On_Flag_Default(38) == sym%NO) then
      return 
    end if
    
  else if (tsfc_onechannel) then
    if ( first_segment)  call MESG ('LST Retrieval with single channel method ', level = verb_lev %DEFAULT)
    !--- if no ch31, abort
    if (Sensor%Chan_On_Flag_Default(31) == sym%NO) then
      return 
    end if
  else
   !--- read LST regcoefs from binary file
    Coef_Lun=GETLUN()
    open(unit=Coef_Lun, file=TRIM(Coef_Fn),recl=leng*4, form ='unformatted', access = 'direct',status='old',action = 'read', iostat=ios)
    if (ios /= 0) THEN
      write(*,*) 'ERROR: Opening LST Coef binary file', TRIM(Coef_Fn)
      stop
    endif
            
    READ(unit=Coef_Lun,rec=1, iostat=ERR) Regcoef       
    close(Coef_Lun)
    if ( first_segment) then
      print*, "LST Retrieval with dual channel method -EB "
      write(*,*) 'LST Coef binary file read in', TRIM(Coef_Fn)
    end if
  end if

  line_loop: do Line_Idx=jmin, jmax - jmin + 1
    element_loop: do Elem_Idx= 1, Image%Number_Of_Elements

      !--- check for a bad pixel pixel
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
        cycle
      endif

      if (ABI_Use_104um_Flag) then

        !--- aliases for visual convenience
        Rad11 = ch(38)%Rad_Toa(Elem_Idx,Line_Idx)
        Xnwp = NWP_PIX%i_nwp(Elem_Idx,Line_Idx)                         !nwp latitude cell
        Ynwp = NWP_PIX%j_nwp(Elem_Idx,Line_Idx)                         !nwp longitude cell
        Rad11_Atm = ch(38)%Rad_Atm(Elem_Idx,Line_Idx)           !10 micron atmospheric radiance
        Trans11_Atm = ch(38)%Trans_Atm(Elem_Idx,Line_Idx)       !10 micron atmospheric transmittance
        Emiss_Sfc11 = ch(38)%Sfc_Emiss(Elem_Idx,Line_Idx)       !10 micron surface emissivity
        Rad11_Atm_Dwn_Sfc = ch(38)%Rad_Atm_Dwn_Sfc(Elem_Idx,Line_Idx)  !10 micron atmospheric radiance down at sfc

        !--- compute the radiance coming off the surface
        Rad11_Sfc = (Rad11 - Rad11_Atm) / Trans11_Atm - &
                  (1.0-Emiss_Sfc11)*Rad11_Atm_Dwn_Sfc

        !--- compute to a temperature
        Trad_Retrieved(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(38,Rad11_Sfc)

        !--- adjust for surface emissivity - this is now the black body emission at Tsfc
        B11_Sfc = Rad11_Sfc / Emiss_Sfc11

        !--- compute to a temperature
        Tsfc_Retrieved(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(38,B11_Sfc)

      elseif (tsfc_onechannel) then

        !--- aliases for visual convenience
        Rad11 = ch(31)%Rad_Toa(Elem_Idx,Line_Idx)
        Xnwp = NWP_PIX%i_nwp(Elem_Idx,Line_Idx)                         !nwp latitude cell
        Ynwp = NWP_PIX%j_nwp(Elem_Idx,Line_Idx)                         !nwp longitude cell
        Rad11_Atm = ch(31)%Rad_Atm(Elem_Idx,Line_Idx)           !11 micron atmospheric radiance
        Trans11_Atm = ch(31)%Trans_Atm(Elem_Idx,Line_Idx)       !11 micron atmospheric transmittance
        Emiss_Sfc11 = ch(31)%Sfc_Emiss(Elem_Idx,Line_Idx)       !11 micron surface emissivity
        Rad11_Atm_Dwn_Sfc = ch(31)%Rad_Atm_Dwn_Sfc(Elem_Idx,Line_Idx)  !11 micron atmospheric radiance down at sfc

        !--- compute the radiance coming off the surface
        Rad11_Sfc = (Rad11 - Rad11_Atm) / Trans11_Atm - &
                  (1.0-Emiss_Sfc11)*Rad11_Atm_Dwn_Sfc

        !--- compute to a temperature
        Trad_Retrieved(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(31,Rad11_Sfc)

        !--- adjust for surface emissivity - this is now the black body emission at Tsfc
        B11_Sfc = Rad11_Sfc / Emiss_Sfc11

        !--- compute to a temperature
        Tsfc_Retrieved(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(31,B11_Sfc)
      
      else
      
        !--- aliases for visual convenience
        !Sec_Sat_Zen=Geo%Seczen(Elem_Idx,Line_Idx)	
        Sat_Zen=Geo%Satzen(Elem_Idx,Line_Idx)	
        Sol_Zen=Geo%Solzen(Elem_Idx,Line_Idx)	
	      TPW = NWP_PIX%Tpw(Elem_Idx,Line_Idx)      
        Emiss_Sfc11 = ch(31)%Sfc_Emiss(Elem_Idx,Line_Idx)       !11 micron surface emissivity
        Emiss_Sfc12 = ch(32)%Sfc_Emiss(Elem_Idx,Line_Idx)       !12 micron surface emissivity
           

	      ! day/night TPW and viewing angle classification -to set index 	
	      dnflag=0
	      if (Sol_Zen <=  85.) dnflag=1
	
	      if (TPW <= 1.5) tpwcls=1
	      if (TPW > 1.5 .and. TPW <= 3.0) tpwcls=2
	      if (TPW > 3.0) tpwcls=3
	
	      if ( Sat_Zen <= 25.) angcls=1
	      if ( Sat_Zen > 25. .and. Sat_Zen <= 45.) angcls=2
	      if ( Sat_Zen > 45. .and. Sat_Zen <= 55.) angcls=3
	      if ( Sat_Zen > 55.) angcls=4

        do Coef_Idx = 1,ncoef
	        RegCoef_Idx = dnflag*nbtpwcls*nbangcls*ncoef + &
                             (tpwcls-1)*nbangcls*ncoef + (angcls - 1) * ncoef + Coef_Idx
          Coef(Coef_Idx) = Regcoef(RegCoef_Idx)
        end do
		
        !--- compute emiss mean
	      Emiss_Mean=(Emiss_Sfc11+Emiss_Sfc12)/2.
	      Emiss_Diff=Emiss_Sfc11-Emiss_Sfc12
	 	
         !--- adjust for surface emissivity - this is now the black body emission at Tsfc
 	      B11_Sfc = ch(31)%BT_toa(Elem_Idx,Line_Idx)
	      B12_Sfc = ch(32)%BT_toa(Elem_Idx,Line_Idx)
	      Bdiff=B11_Sfc - B12_Sfc
	 
        !--- compute to a temperature
        Tsfc_Retrieved(Elem_Idx,Line_Idx) =   Coef(1) * B11_Sfc &
                                            + Coef(2) * Bdiff &
                                            + Coef(3) * Emiss_Mean &
	                                          + Coef(4) * Emiss_Mean*Bdiff &
                                            + Coef(5) * Emiss_Diff &
                                            + Coef(6)
      
      endif

    end do element_loop
  end do line_loop
  
  first_segment = .false.

end subroutine COMPUTE_TSFC

!======================================================================
! Normalize the reflectances by the solar zenith angle cosine
!======================================================================
 subroutine NORMALIZE_REFLECTANCES(Sun_Earth_Distance)
  real(kind=real4), intent(in):: Sun_Earth_Distance
  integer:: i,j, Chan_Idx
  real:: Factor

  ! for these sensors, no correction is needed
  if (trim(Sensor%Sensor_Name) == 'VIIRS') return 
  if (trim(Sensor%Sensor_Name) == 'VIIRS-ifF') return 
  if (trim(Sensor%Sensor_Name) == 'AVHRR-ifF') return 
  if (trim(Sensor%Sensor_Name) == 'METIMAGE') return 

  !--------------------------------------------------------------------
  ! loop through pixels and apply normalization factor
  !--------------------------------------------------------------------
  do j = 1, Image%Number_Of_Lines_Read_This_Segment

     do i = 1, Image%Number_Of_Elements

      if (Bad_Pixel_Mask(i,j) == sym%NO .and. Geo%Cossolzen(i,j) > 0.0) then

       Factor = 1.0 / Geo%Cossolzen(i,j)

       ! for these sensors, a correction for sun earth distance is also needed
       if ( (trim(Sensor%Sensor_Name) == 'AVHRR-1') .or. &
            (trim(Sensor%Sensor_Name) == 'AVHRR-2') .or. &
            (trim(Sensor%Sensor_Name) == 'AVHRR-3') .or. &
            (trim(Sensor%Sensor_Name) == 'AVHRR-FUSION') .or. &
            (trim(Sensor%Sensor_Name) == 'GOES-IL-IMAGER') .or. &
            (trim(Sensor%Sensor_Name) == 'GOES-MP-IMAGER') .or. &
            (trim(Sensor%Sensor_Name) == 'GOES-IP-SOUNDER') .or. &
            (trim(Sensor%Sensor_Name) == 'SEVIRI') .or. &
            (trim(Sensor%Sensor_Name) == 'MTSAT-IMAGER') .or. &
            (trim(Sensor%Sensor_Name) == 'COMS-IMAGER') .or. &
            (trim(Sensor%Sensor_Name) == 'FY2-IMAGER')) then

            Factor = Factor * (Sun_Earth_Distance**2) 

       endif

       ! apply correction
       do Chan_Idx = 1,19
          if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
            if (ch(Chan_Idx)%Ref_Toa(i,j) /= Missing_Value_Real4) then
             ch(Chan_Idx)%Ref_Toa(i,j) = ch(Chan_Idx)%Ref_Toa(i,j) * Factor
            endif
            if (allocated(ch(Chan_Idx)%Ref_Toa_Min_Sub)) then
               if (ch(Chan_Idx)%Ref_Toa_Min_Sub(i,j) /= Missing_Value_Real4) then
                  ch(Chan_Idx)%Ref_Toa_Min_Sub(i,j) = ch(Chan_Idx)%Ref_Toa_Min_Sub(i,j) * Factor
               endif
            endif
            if (allocated(ch(Chan_Idx)%Ref_Toa_Max_Sub)) then
               if (ch(Chan_Idx)%Ref_Toa_Max_Sub(i,j) /= Missing_Value_Real4) then
                  ch(Chan_Idx)%Ref_Toa_Max_Sub(i,j) = ch(Chan_Idx)%Ref_Toa_Max_Sub(i,j) * Factor
               endif
            endif
          endif
       enddo

       if (Sensor%Chan_On_Flag_Default(26) == sym%YES) then
          if (ch(26)%Ref_Toa(i,j) /= Missing_Value_Real4) then
            ch(26)%Ref_Toa(i,j) = ch(26)%Ref_Toa(i,j) * Factor
          endif
       endif

       !  normalize by sun angle the dark sky composite
       if ((Sensor%Chan_On_Flag_Default(1) == sym%YES)) then
         if (Ref_Ch1_Dark_Composite(i,j) /= Missing_Value_Real4) then
           Ref_Ch1_Dark_Composite(i,j) = Ref_Ch1_Dark_Composite(i,j) * Factor
         endif
       endif


       !--- for avhrr, handle absense of ch6
       if (index(Sensor%Sensor_Name,'AVHRR') > 0 .and. &
           Sensor%Chan_On_Flag_Default(6) == sym%YES) then
           if (ch(6)%Ref_Toa(i,j) < 0) ch(6)%Ref_Toa(i,j) = Missing_Value_Real4
           if (Ch3a_On_Avhrr(j) /= sym%YES) ch(6)%Ref_Toa(i,j) = Missing_Value_Real4
       endif

       !--- in terminator region, renormalize Channel 1 (maybe extend to all?)
       if (Geo%Solzen(i,j) > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) then
          if (Sensor%Chan_On_Flag_Default(1) == sym%YES)  &
              ch(1)%Ref_Toa(i,j) = TERM_REFL_NORM(Geo%Cossolzen(i,j),ch(1)%Ref_Toa(i,j))
       endif

      else

       !--- set to missing
       do Chan_Idx = 1,19
          if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) ch(Chan_Idx)%Ref_Toa(i,j) = Missing_Value_Real4
       enddo
       if (Sensor%Chan_On_Flag_Default(26) == sym%YES) ch(26)%Ref_Toa(i,j) = Missing_Value_Real4

      endif

     end do
     
   end do

end subroutine NORMALIZE_REFLECTANCES
!----------------------------------------------------------------------
! Compute Channel3b albedo
!
! input
!    Ch3a_on - 0 if Ch3b is on, 1 if Ch3a is on
!    Sun_Earth_Distance - sun-earth distance factor
!    Rad_Ch20 - vector of Ch3b radiance (mW/m^2/str/cm^1)
!    Bt_Ch31 - vector of ch4 brightness temp. (K)
!    Bt_Ch32 - vector of ch5 brightness temp. (K)
!    Solzen - vector of solar zenith angles (degrees)
!
!  output (passed through shared memory)
!    Ref_Ch20 - Ch3b albeDO computed using standard ch4 based method
!    Emiss_Ch20 - Ch3b emissivity computed using standrad Ch3 based method
!
!  internal
!    Rad_Ch20_ems - Ch3b emission radiance (mW/m^2/str/cm^-1)
!
! Note, Ref_Ch20_Sfc is computed in ATMOS_CORR
!
! Revision History
!   January 2003 - A. Heidinger
!
! --->    AW 10/20/2014
! ch(20) %ref_toa is the pseudo solar reflectance in 3.9 channels
!  Rad_obs = Rad_sol + ( 1 - R ) Rad_ch20_ems
!  Rad_obs = (R * F_0 * mu) / PI + ( 1 - R ) Rad_ch20_ems
!  == >   R = ( PI (Rad_obs - Rad_ch20_ems )) / ( F_o * mu - PI * Rad_ch20_ems )
!     see Kaufman and Remer IEEE 1994:
!   "Detection of  Forests Using Mid-IR Reflectance: An  Application for Aerosol Studies"
!   http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=297984
!
!-----------------------------------------------------------------------
subroutine CH20_PSEUDO_REFLECTANCE(Solar_Ch20_Nu,Cos_Solzen,Rad_Ch20,Ref_Chan_Idx,Bt_Ch31,Sun_Earth_Distance,Ref_Ch20,Emiss_Ch20)

  real(kind=real4), intent(in):: Solar_Ch20_Nu
  real(kind=real4), intent(in):: Sun_Earth_Distance
  integer(kind=int4), intent(in):: Ref_Chan_Idx
  real(kind=real4), dimension(:,:), intent(in):: Cos_Solzen
  real(kind=real4), dimension(:,:), intent(in):: Rad_Ch20
  real(kind=real4), dimension(:,:), intent(in):: Bt_Ch31
  real(kind=real4), dimension(:,:), intent(out):: Ref_Ch20
  real(kind=real4), dimension(:,:), intent(out):: Emiss_Ch20
  integer:: Elem_Idx, Line_Idx, Num_Elements, Num_Lines
  real :: Rad_Ch20_Ems
  real :: Solar_Irradiance

  Num_Elements = Image%Number_Of_Elements        !make local copy of a global variable
  Num_Lines = Image%Number_Of_Lines_Per_Segment  !make local copy of a global variable

  if ((Sensor%Chan_On_Flag_Default(20) /= sym%YES) .or.  &
      (Sensor%Chan_On_Flag_Default(Ref_Chan_Idx) /= sym%YES)) then   !start Ch3a_on check
        return
  endif

  !----------------------------------------------------------------------------
  !--- standard Ref_Ch20 computation
  !---------------------------------------------------------------------------
  do Line_Idx = 1,Num_Lines
      do Elem_Idx = 1, Num_Elements

        !--- check for bad scans
        if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
          cycle
        endif

        if (Bt_Ch31(Elem_Idx,Line_Idx) > 180.0) then
           Rad_Ch20_Ems = PLANCK_RAD_FAST(20,Bt_Ch31(Elem_Idx,Line_Idx))
           Emiss_Ch20(Elem_Idx,Line_Idx) = Rad_Ch20(Elem_Idx,Line_Idx) / Rad_Ch20_Ems
        else
           Rad_Ch20_Ems = Missing_Value_Real4
           Emiss_Ch20(Elem_Idx,Line_Idx) = Missing_Value_Real4
        endif
         
        if ((Rad_Ch20_Ems>0.0).and.(Rad_Ch20(Elem_Idx,Line_Idx)>0.0)) then
           Solar_Irradiance = max (0.0, (Solar_Ch20_Nu*Cos_Solzen(Elem_Idx,Line_Idx))/(Sun_Earth_Distance**2))
           Ref_Ch20(Elem_Idx,Line_Idx) = 100.0*pi*(Rad_Ch20(Elem_Idx,Line_Idx)-Rad_Ch20_Ems) /  &
                                                  (Solar_Irradiance - pi*Rad_Ch20_Ems)
        endif

        !--- constrain values
        if (Ref_Ch20(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Ref_Ch20(Elem_Idx,Line_Idx) = max(-50.0,min(100.0,Ref_Ch20(Elem_Idx,Line_Idx)))
        endif

      enddo

   enddo

end subroutine CH20_PSEUDO_REFLECTANCE

!------------------------------------------------------------------------
! Routine to assign clear-sky quality flags
!------------------------------------------------------------------------
subroutine ASSIGN_CLEAR_SKY_QUALITY_FLAGS(jmin,jmax)

    integer, intent(in):: jmin,jmax
    integer:: i, j, i1, i2, j1, j2, n, max_Mask

!-- determine size of box for spatial filter
 n = 1
                                                                                                                                                
!--- initialize
 Tsfc_Qf = 0
 Ndvi_Qf = 0
 Rsr_Qf = 0


 do j = jmin, jmin+jmax-1
                                                                                                                                                
   !--- determine y-dimensions of array to check
   j1 = max(jmin,j-n)
   j2 = min(jmax,j+n)
                                                                                                                                                
    do i = 1, Image%Number_Of_Elements

      !--- check for bad scans
      if (Bad_Pixel_Mask(i,j) == sym%YES) then
        cycle
      endif
                                                                                                                                                
      !--- determine x-dimensions of array to check
      i1 = max(1,i-n)
      i2 = min(Image%Number_Of_Elements,i+n)
        
      !--- initial cloud mask based
      if (CLDMASK%Cld_Mask(i,j) == sym%CLEAR) then
         Tsfc_Qf(i,j) = 3
         Ndvi_Qf(i,j) = 3
         Rsr_Qf(i,j) = 3
      endif
      if (CLDMASK%Cld_Mask(i,j) == sym%PROB_CLEAR) then
         Tsfc_Qf(i,j) = 2
         Ndvi_Qf(i,j) = 2
         Rsr_Qf(i,j) = 2
      endif
      if (CLDMASK%Cld_Mask(i,j) == sym%PROB_CLOUDY) then
         Tsfc_Qf(i,j) = 1
         Ndvi_Qf(i,j) = 1
         Rsr_Qf(i,j) = 1
      endif
      if (CLDMASK%Cld_Mask(i,j) == sym%CLOUDY) then
         Tsfc_Qf(i,j) = 0
         Ndvi_Qf(i,j) = 0
         Rsr_Qf(i,j) = 0
      endif

      !--- assign Ndvi over water to be low quality
      if (Sfc%Land_Mask(i,j) == sym%NO) then
        Ndvi_Qf(i,j) = 0
      endif

      !--- assign Ndvi at high angles to be low quality
      if (Geo%Solzen(i,j) > 75.0) then
       Ndvi_Qf(i,j) = min(1,int(Ndvi_Qf(i,j)))
      endif

      !--- modifcations of Rsr quality
      if (Sfc%Land_Mask(i,j) == sym%YES) then    !ocean only
        Rsr_Qf(i,j) = 0
      endif
      if (Geo%Solzen(i,j) > 75.0) then    !sufficient light
        Rsr_Qf(i,j) = min(1,int(Rsr_Qf(i,j)))
      endif
      if (Geo%Glintzen(i,j) < Glint_Zen_Thresh) then   !outside glint
       Rsr_Qf(i,j) = min(1,int(Rsr_Qf(i,j)))
      endif
      if (Sfc%Snow(i,j) /= sym%NO_SNOW) then    !Snow
       Rsr_Qf(i,j) = min(1,int(Rsr_Qf(i,j)))
      endif


      !--- aot 
       if (Aerosol_Mode > 0) then

        Aot_Qf(i,j) = 0
        
        if ((CLDMASK%Cld_Mask(i,j) == sym%CLEAR) .and.  &
            (Sfc%Land_Mask(i,j) == sym%NO) .and.  &
            (Geo%Solzen(i,j) < 70.00) .and.  &
            (Sfc%Snow(i,j) == sym%NO_SNOW)) then
          Aot_Qf(i,j) = 1
          if (Geo%Glintzen(i,j) > Glint_Zen_Thresh) then
            if (Geo%Relaz(i,j) > 90.0) then
              Aot_Qf(i,j) = 3
            else
             Aot_Qf(i,j) = 2 !- 1
           endif
         endif
        endif
                                                                                                                              
        !--- assign aerosol over Snow to be of low quality
        if (Sfc%Snow(i,j) /= sym%NO_SNOW) then
         Aot_Qf(i,j) = 0
        endif

        !--- assign high quality pixels around a cloudy results as qf = 2
        max_Mask = maxval(CLDMASK%Cld_Mask(i1:i2,j1:j2))
        if (max_Mask >= 2) then
         if (Aot_Qf(i,j) == 3) then
           Aot_Qf(i,j) = 2
         endif
         if (Ndvi_Qf(i,j) == 3) then
           Ndvi_Qf(i,j) = 2
         endif
         if (Tsfc_Qf(i,j) == 3) then
           Tsfc_Qf(i,j) = 2
         endif
        endif

     !--- forcing the reporting of aerosol for this condition (A. Evan)
     if (CLDMASK%Dust_Mask(i,j) == sym%YES) then
         Aot_Qf(i,j) = 3
     endif

  endif ! end of Aerosol  QF


    end do
  end do

  end subroutine ASSIGN_CLEAR_SKY_QUALITY_FLAGS

!----------------------------------------------------------------------
!--- Compute a mask identifying presence of oceanic glint
!--- 
!--- input and output passed through global arrays
!----------------------------------------------------------------------
subroutine COMPUTE_GLINT(Source_GLintzen, Source_Ref_Toa, Source_Ref_Std_3x3, &
                         Source_Glint_Mask,ABI_Use_104um_Flag)

  real, dimension(:,:), intent(in):: Source_GlintZen
  real, dimension(:,:), intent(in):: Source_Ref_Toa
  real, dimension(:,:), intent(in):: Source_Ref_Std_3x3
  integer(kind=int1),  dimension(:,:), intent(out):: Source_Glint_Mask
  logical, intent(in):: ABI_Use_104um_Flag

  !--- define local variables
  integer:: Number_Of_Lines
  integer:: Number_Of_Elements
  integer:: Elem_Idx
  integer:: Line_Idx
  real:: Refl_Thresh

  !--- alias some global sizes into local values
  Number_Of_Lines = Image%Number_Of_Lines_Per_Segment
  Number_Of_Elements = Image%Number_Of_Elements

  Source_Glint_Mask = Missing_Value_Int1

     line_loop: do Line_Idx = 1, Number_Of_Lines

     element_loop: do Elem_Idx = 1, Number_Of_Elements

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
             cycle
     endif

     !--- initialize valid pixel to no
     Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO

     !--- skip land pixels
     if ((Sfc%Land_Mask(Elem_Idx,Line_Idx) == sym%NO) .and. &
          Sfc%Snow(Elem_Idx,Line_Idx) == sym%NO_SNOW) then

       !--- turn on in geometric glint cone and sufficient Ref_Ch1
       if ((Source_Glintzen(Elem_Idx,Line_Idx) < Glint_Zen_Thresh)) then

          !--- assume to be glint if in geometric zone
          Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%YES

          if (ABI_Use_104um_Flag) then

            if (Sensor%Chan_On_Flag_Default(38) == sym%YES) then

              !--- exclude pixels colder than the freezing temperature
              if (ch(38)%Bt_Toa(Elem_Idx,Line_Idx) < 273.15) then
                Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
                cycle
              endif

              !--- exclude pixels colder than the surface
              if (ch(38)%Bt_Toa(Elem_Idx,Line_Idx) < ch(38)%Bt_Toa_Clear(Elem_Idx,Line_Idx) - 5.0) then
                Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
                cycle
              endif

            endif !--- Chan 38

          else

            if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then

              !--- exclude pixels colder than the freezing temperature
              if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < 273.15) then
                Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
                cycle
              endif

              !--- exclude pixels colder than the surface
              if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < ch(31)%Bt_Toa_Clear(Elem_Idx,Line_Idx) - 5.0) then
                Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
                cycle
              endif

            endif !--- Chan 31

          endif   !--- ABI Use 104 Flag

          !-turn off if non-uniform - but not near limb
          if (Geo%Satzen(Elem_Idx,Line_Idx) < 45.0) then 
            if (ABI_Use_104um_Flag) then
              if (Sensor%Chan_On_Flag_Default(38) == sym%YES) then
                if (ch(38)%Bt_Toa_Std_3x3(Elem_Idx,Line_Idx) > 1.0) then
                  Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
                  cycle
                endif
              endif
            else
              if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
                if (ch(31)%Bt_Toa_Std_3x3(Elem_Idx,Line_Idx) > 1.0) then
                  Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
                  cycle
                endif
              endif
            endif !--- ABI Use 104 um Flag

            if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
              if (Source_Ref_Std_3x3(Elem_Idx,Line_Idx) > 2.0) then
                Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
                cycle
              endif
            endif
          endif !--- Geo Satzen

          !-checks on the value of ch1
          if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then

            !-turn off if dark
            if (Source_Ref_Toa(Elem_Idx,Line_Idx) < 5.0) then
             Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
             cycle
            endif

            !-turn off if bright
            if (Source_Glintzen(Elem_Idx,Line_Idx) > 10.0 .and. &
                Source_Glintzen(Elem_Idx,Line_Idx) < 40.0) then

               Refl_Thresh = 25.0 - Source_Glintzen(Elem_Idx,Line_Idx)/3.0

               if (Source_Ref_Toa(Elem_Idx,Line_Idx) > Refl_Thresh) then
                  Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
                  cycle
               endif

            endif

!           if (Source_Glintzen(Elem_Idx,Line_Idx) > 20.0) then
!             if (Source_Ref_Toa(Elem_Idx,Line_Idx) > 15.0) then
!                 Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
!                 cycle
!             endif
!           endif

!           if (Source_Glintzen(Elem_Idx,Line_Idx) > 10.0) then
!             if (Source_Ref_Toa(Elem_Idx,Line_Idx) > 20.0) then
!                 Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
!                 cycle
!             endif
!           endif

          endif

       endif  !Glintzen check

     endif    !land check

     enddo element_loop
   enddo line_loop


 end subroutine COMPUTE_GLINT

!----------------------------------------------------------------------
! Read MODIS white sky albedoes
!----------------------------------------------------------------------
subroutine READ_MODIS_WHITE_SKY_ALBEDO(modis_alb_id,modis_alb_str,Ref_Sfc_White_Sky)

    integer(kind=4), intent(in):: modis_alb_id
    TYPE(Land_grid_description), intent(in) :: modis_alb_str
    real(kind=real4), dimension(:,:), intent(out):: Ref_Sfc_White_Sky
    integer(kind=int2), dimension(:,:),allocatable :: raw
    integer :: dim_arr(2)
    
    dim_arr = shape(nav%lat)
    allocate ( raw(dim_arr(1),dim_arr(2)))
    
    CALL READ_LAND_SFC_HDF(modis_alb_id, modis_alb_str, Nav%Lat, &
                          Nav%Lon, Geo%Space_Mask, raw)
                       
    Ref_Sfc_White_Sky = 0.1* raw

!---->    Ref_Sfc_White_Sky = 1.10*Ref_Sfc_White_Sky   !EMPIRICAL ADJUSTMENT

    where(raw == 32767)
             Ref_Sfc_White_Sky = Missing_Value_Real4
    endwhere

    !--- modify for water
    where(Sfc%Land_Mask == sym%NO)
            Ref_Sfc_White_Sky = Ref_Sfc_White_Sky_Water
    endwhere

end subroutine READ_MODIS_WHITE_SKY_ALBEDO

!==============================================================================
!
!==============================================================================
 subroutine DETERMINE_LEVEL1B_COMPRESSION(File_1b_Original,L1b_Gzip,L1b_Bzip2)
   character(len=*), intent(in):: File_1b_Original
   integer(kind=int4), intent(out):: L1b_Gzip
   integer(kind=int4), intent(out):: L1b_Bzip2
   character(len=1020):: System_String
   character(len=7):: L1b_ext


  !--- determine if the goes data is compressed
  L1b_ext = File_1b_Original(len_trim(File_1b_Original)-2: &
                             len_trim(File_1b_Original))

  !-- determine if gzipped
  if (trim(L1b_ext) == '.gz') then
     L1b_Gzip = sym%YES
  else
     L1b_Gzip = sym%NO
  endif

  !--- check if bzipped
  if (trim(L1b_ext) == 'bz2') then
     L1b_Bzip2 = sym%YES
  else
     L1b_Bzip2 = sym%NO
  endif

  !--- uncompress
  if (L1b_Gzip == sym%YES) then
     Image%Level1b_Name = File_1b_Original(1:len(trim(File_1b_Original))-3)
     System_String = "gunzip -c "//trim(Image%Level1b_Path)//trim(File_1b_Original)// &
        " > "//trim(Temporary_Data_Dir)//trim(Image%Level1b_Name)
     call SYSTEM(System_String)

     Number_of_Temporary_Files = Number_of_Temporary_Files + 1
     Temporary_File_Name(Number_of_Temporary_Files) = trim(Image%Level1b_Name)

  elseif (L1b_Bzip2 == sym%YES) then
     Image%Level1b_Name = File_1b_Original(1:len(trim(File_1b_Original))-4)
     System_String = "bunzip2 -c "//trim(Image%Level1b_Path)//trim(File_1b_Original)// &
        " > "//trim(Temporary_Data_Dir)//trim(Image%Level1b_Name)
     call SYSTEM(System_String)

     Number_of_Temporary_Files = Number_of_Temporary_Files + 1
     Temporary_File_Name(Number_of_Temporary_Files) = trim(Image%Level1b_Name)

  else
     Image%Level1b_Name = trim(File_1b_Original)
  endif

   !--- make a full file name
   if (L1b_Gzip == sym%YES .or. L1b_bzip2 == sym%YES) then
     Image%Level1b_Full_Name = trim(Temporary_Data_Dir)//trim(Image%Level1b_Name)
   else
    Image%Level1b_Full_Name = trim(Image%Level1b_Path)//trim(Image%Level1b_Name)
   endif

 end subroutine DETERMINE_LEVEL1B_COMPRESSION

!====================================================================
!
! Attempt to fix the land classification based on observed ndvi
!
! if the ndvi is high and the land class is not land, this pixel should be land
! if the ndvi is low and the land class is land, this pixel should be water
!
!====================================================================
subroutine MODIFY_LAND_CLASS_WITH_NDVI(Line_Idx_Min,Num_Lines)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines

  integer:: Line_Idx_Max
  integer:: Elem_Idx_Max
  integer:: Elem_Idx_Min
  integer:: Num_Elements
  integer:: Line_Idx
  integer:: Elem_Idx
  real:: ndvi_temp
  real, parameter:: Ndvi_Land_Threshold = 0.25
  real, parameter:: Ndvi_Water_Threshold = -0.25
  real, parameter:: Solzen_Threshold = 60.0
  real, parameter:: Ref_Ch2_Threshold = 60.0

  !--- do not do this for advanced geos since sampling can cause issues
  if (Sensor%WMO_ID == 173) return
  if (Sensor%WMO_ID == 174) return
  if (Sensor%WMO_ID == 270) return
  if (Sensor%WMO_ID == 271) return
  if (Sensor%WMO_ID == 272) return
  if (Sensor%WMO_ID == 530) return

  Elem_Idx_Min = 1
  Num_Elements = Image%Number_Of_Elements
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  line_loop: do Line_Idx = Line_Idx_Min, Line_Idx_Max
    element_loop: do Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

    if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle
    if (Sensor%Chan_On_Flag_Default(1) == sym%NO) cycle
    if (Sensor%Chan_On_Flag_Default(2) == sym%NO) cycle
    if (Geo%Solzen(Elem_Idx,Line_Idx) > Solzen_Threshold) cycle
    if (index(Sensor%Sensor_Name,'MODIS') > 0) cycle                         !modis ch2 saturates, need to modify for MODIS
    if (index(Sensor%Sensor_Name,'GOES-RU-IMAGER') > 0) cycle

    Ndvi_Temp = (ch(2)%Ref_Toa(Elem_Idx,Line_Idx) - ch(1)%Ref_Toa(Elem_idx,Line_Idx)) / &
                (ch(2)%Ref_Toa(Elem_Idx,Line_Idx) + ch(1)%Ref_Toa(Elem_idx,Line_Idx)) 

    if (Ndvi_Temp > Ndvi_Land_Threshold) then
      Sfc%Land(Elem_Idx,Line_Idx) = sym%LAND
    endif

    if (Ndvi_Temp < Ndvi_Water_Threshold .and. &
        Sfc%Land(Elem_Idx,Line_Idx) == sym%LAND) then
        Sfc%Land(Elem_Idx,Line_Idx) = sym%SHALLOW_INLAND_WATER
    endif

    enddo element_loop
  enddo line_loop

end subroutine MODIFY_LAND_CLASS_WITH_NDVI
!====================================================================
! Function Name: TERM_REFL_NORM
!
! Function:
!    Renormalize reflectances to improve performance near the terminator 
! using the parameteization given by Li et. al. 2006
!
! Description: Renormalizes reflectances in the terminator region
!   
! Calling Sequence: Refl_Chn2 = TERM_REFL_NORM(Cos_Sol_Zen,Refl_Chn2)
!   
!
! Inputs:
!   Cosine of the Solar Zenith Angle
!   Channel 2 reflectance that is normalized by cosine of solar zenith
!
! Outputs: 
!   Renormalized reflectance
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Li et. al. 2006
!
!====================================================================
 FUNCTION TERM_REFL_NORM(Cos_Sol_Zen,Reflectance)  &
          RESULT(Reflectance_Normalized)

   real(kind=real4), intent(in):: Cos_Sol_Zen
   real(kind=real4), intent(in):: Reflectance
   real(kind=real4):: Reflectance_Normalized
   real(kind=real4):: Norm_Param

   Reflectance_Normalized = Reflectance * Cos_Sol_Zen

   Norm_Param = 24.35 / (2*Cos_Sol_Zen + sqrt(498.5225*(Cos_Sol_Zen**2) + 1) )

   Reflectance_Normalized = Reflectance_Normalized*Norm_Param

 end FUNCTION TERM_REFL_NORM


!====================================================================
! Routine Name: MERGE_NWP_HIRES_ZSFC
!
! Function:
! Merge the high resolution and low resolution surface elevation
! fields into one single field
!
! Inputs: 
!    Zsfc_Hires - passed via global arrays
!    Zsfc_Nwp - passed via global arrays
!
! Outputs: 
!    Zsfc - passed via global arrays
!====================================================================
subroutine MERGE_NWP_HIRES_ZSFC(Line_Idx_Min,Num_Lines)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  integer:: Num_Elements
  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: Elem_Idx_Min
  integer:: Elem_Idx_Max
  integer:: Line_Idx_Max
  integer:: Lat_NWP_Idx
  integer:: Lon_NWP_Idx


  Elem_Idx_Min = 1
  Num_Elements = Image%Number_Of_Elements
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  line_loop: do Line_Idx = Line_Idx_Min, Line_Idx_Max
    element_loop: do Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

      !--- if no, geolocation, set to missing and go to next pixel
      if (Geo%Space_Mask(Elem_Idx,Line_Idx) ) then
          Sfc%Zsfc(Elem_Idx,Line_Idx) = Missing_Value_Real4
          cycle
      endif
 
      !--- if hires value available use it, or try NWP
      if (Sfc%Zsfc_Hires(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then 

           Sfc%Zsfc(Elem_Idx,Line_Idx) = Sfc%Zsfc_Hires(Elem_Idx,Line_Idx)

      !--- if this is ocean pixel, assume zero
      elseif (Sfc%Land(Elem_Idx,Line_Idx) == sym%SHALLOW_OCEAN .or. &
              Sfc%Land(Elem_Idx,Line_Idx) == sym%MODERATE_OCEAN .or. &
              Sfc%Land(Elem_Idx,Line_Idx) == sym%DEEP_OCEAN) then

              Sfc%Zsfc(Elem_Idx,Line_Idx) = 0.0     

      !--- try NWP
      else

         Lon_Nwp_Idx = NWP_PIX%I_Nwp(Elem_Idx,Line_Idx)
         Lat_Nwp_Idx = NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)

         !--- if nwp not available, assume zero
         if (Lon_Nwp_Idx > 0 .and. Lat_Nwp_Idx > 0) then
           Sfc%Zsfc(Elem_Idx,Line_Idx) = NWP%Zsfc(Lon_Nwp_Idx, Lat_Nwp_Idx)
         else
           Sfc%Zsfc(Elem_Idx,Line_Idx) = 0.0     
         endif

      endif
          
    enddo element_loop
  enddo line_loop

end subroutine MERGE_NWP_HIRES_ZSFC

!====================================================================
! Routine Name:
!
! Function:
! Compute maximum cloud mask number in adjusted 3x3 pixels
!
! Inputs:
!    Line_Start - minimum line
!    Number_of_Lines - lines to read
!
! Outputs:
!    Adj_Pix_Cld_Mask - passed via global arrays
!====================================================================
subroutine ADJACENT_PIXEL_CLOUD_MASK(Line_Start,Number_of_Lines)

  integer (kind=int4), intent(in):: Line_Start
  integer (kind=int4), intent(in):: Number_of_Lines
  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: Number_of_Elements
  integer:: i1,i2,j1,j2

  Number_of_Elements = Image%Number_Of_Elements

  CLDMASK%Adj_Pix_Cld_Mask = Missing_Value_Int1

  line_loop: do Line_Idx = Line_Start, Number_of_Lines + Line_Start - 1

    j1 = max(1,Line_Idx - 1)
    j2 = min(Number_of_Lines,Line_Idx + 1)

    element_loop: do Elem_Idx = 1, Number_of_Elements

      i1 = max(1,Elem_Idx - 1)
      i2 = min(Number_of_Elements,Elem_Idx + 1)

      CLDMASK%Adj_Pix_Cld_Mask(Elem_Idx,Line_Idx) = maxval(CLDMASK%Cld_Mask(i1:i2,j1:j2))

    enddo element_loop
  enddo line_loop

end subroutine ADJACENT_PIXEL_CLOUD_MASK

!-----------------------------------------------------------
! Determine an ACHA Success fraction
!
! Success fraction will be defined as the number of points
! with a valid Tc versus the total number of points
! processed through ACHA.  This should not points where
! no retrieval was attempted (ie clear pixels)
!
! Note this arrays uses the global variable holding the
! packed ACHA quality flags
!
! ACHA Quality Flags
! 1 - Processed (1 = yes / 0 = no)
! 2 - Valid Tc Retrieval (1 = yes, 0 = no)
! 3 - Valid ec Retrieval (1 = yes, 0 = no)
! 4 - Valid beta Retrieval (1 = yes, 0 = no)
! 5 - degraded Tc Retrieval (1 = yes, 0 = no)
! 6 - degraded ec Retrieval (1 = yes, 0 = no)
! 7 - degraded beta Retrieval (1 = yes, 0 = no)! 
!-----------------------------------------------------------
subroutine COMPUTE_ACHA_PERFORMANCE_METRICS(Processed_Count,Valid_Count,Success_Fraction)

  real(kind=real4), intent(inout):: Processed_Count
  real(kind=real4), intent(inout):: Valid_Count
  real(kind=real4), intent(out):: Success_Fraction
  integer, parameter:: Count_Min = 10
  real:: Processed_Count_Segment
  real:: Valid_Count_Segment

  Processed_Count_Segment = count(btest(ACHA%Packed_Quality_Flags,0))

  Valid_Count_Segment = count((btest(int(ACHA%Packed_Quality_Flags),1)) .and.  &
                              (btest(int(ACHA%Packed_Quality_Flags),2)) .and. &
                              (btest(int(ACHA%Packed_Quality_Flags),3)))
  
  Processed_Count = Processed_Count + Processed_Count_Segment
  Valid_Count = Valid_Count + Valid_Count_Segment

  if (Processed_Count > Count_Min) then
    Success_Fraction = Valid_Count / Processed_Count
  else
    Success_Fraction = Missing_Value_Real4 
  endif

end subroutine COMPUTE_ACHA_PERFORMANCE_METRICS

!-----------------------------------------------------------
! Determine a DCOMP success fraction
!
! Success fraction will be defined as the number of points
! with a valid COD versus the total number of points
! processed through DCOMP.  This should not points where
! no retrieval was attempted (ie clear pixels)
!
! Note this arrays uses the global variable holding the
! packed DCOMP quality flags
!
!
! DCOMP Quality Flags
!"1:Processed (0=no,1=yes) "// &
!"2:valid COD retrieval (0=yes,1=no) "// &
!"3:valid REF retrieval (0=yes,1=no) "// &
!"4:degraded COD retrieval (0=no,1=degraded) "// &
!"5:degraded REF retrieval (0=no,1=degraded) "// &
!"6:convergency (1=no,0=yes) "// &
!"7:glint (0=no,1=yes) ", &
!-----------------------------------------------------------
subroutine COMPUTE_DCOMP_PERFORMANCE_METRICS(Dcomp_Processed_Count,Dcomp_Valid_Count)

  real(kind=real4), intent(inout):: Dcomp_Processed_Count
  real(kind=real4), intent(inout):: Dcomp_Valid_Count
  real:: Processed_Count_Segment
  real:: Valid_Count_Segment
  real, parameter:: Count_Min = 10.0


  Processed_Count_Segment = count(btest(Dcomp_Quality_Flag,0))
  Valid_Count_Segment = count((.not. btest(Dcomp_Quality_Flag,1)) .and. &
                              (.not. btest(Dcomp_Quality_Flag,2)) .and. &
                              btest(Dcomp_Quality_Flag,0) )

  
  Dcomp_Processed_Count = Dcomp_Processed_Count + Processed_Count_Segment
  Dcomp_Valid_Count = Dcomp_Valid_Count + Valid_Count_Segment

  if (DCOMP_Processed_Count > Count_Min) then
    DCOMP_Success_Fraction = DCOMP_Valid_Count / DCOMP_Processed_Count
  else
    DCOMP_Success_Fraction = Missing_Value_Real4 
  endif

end subroutine COMPUTE_DCOMP_PERFORMANCE_METRICS

   !-----------------------------------------------------------
   ! Determine a Fraction of pixels with a confident cloud mask
   !
   !-----------------------------------------------------------
   subroutine COMPUTE_CLOUD_MASK_PERFORMANCE_METRICS(Cloud_Mask_Count,Nonconfident_Cloud_Mask_Count)
      
     
      integer:: Num_Elements
      integer:: Num_Lines
      real(kind=real4), intent(inout):: Cloud_Mask_Count
      real(kind=real4), intent(inout):: Nonconfident_Cloud_Mask_Count
      integer(kind=int1), dimension(:,:), allocatable:: Mask_local
      integer(kind=int1), dimension(:,:), allocatable:: Nonconfident_Mask_local
      integer, parameter:: Count_Min = 10
      real:: Count_Segment
      real:: Nonconfident_Count_Segment
  
      Num_Elements = Image%Number_Of_Elements  !make local copy of a global variable
      Num_Lines = Image%Number_Of_Lines_Per_Segment  !make local copy of a global variable

      allocate(Mask_local(Num_Elements,Num_Lines))
      allocate(Nonconfident_Mask_local(Num_Elements,Num_Lines))

      Mask_local = 0
      Nonconfident_Mask_local = 0

      where(CLDMASK%Cld_Mask == sym%CLEAR .or. CLDMASK%Cld_Mask == sym%PROB_CLEAR .or.  &
            CLDMASK%Cld_Mask == sym%PROB_CLOUDY .or. CLDMASK%Cld_Mask == sym%Cloudy)
         Mask_local = 1
      end where

      where(CLDMASK%Cld_Mask == sym%PROB_CLEAR .or. CLDMASK%Cld_Mask == sym%PROB_CLOUDY)
         Nonconfident_Mask_local = 1
      end where

      Count_Segment = sum(real(Mask_local))
      if (Count_Segment < Count_Min) then
         return
      end if
  
      deallocate ( mask_local)

      Nonconfident_Count_Segment = sum(real(Nonconfident_Mask_local))
   
      deallocate (Nonconfident_Mask_local) 
   
      Cloud_Mask_Count = Cloud_Mask_Count + Count_Segment
      Nonconfident_Cloud_Mask_Count = Nonconfident_Cloud_Mask_Count + Nonconfident_Count_Segment

      if (Cloud_Mask_Count > Count_Min) then
         Nonconfident_Cloud_Mask_Fraction = Nonconfident_Cloud_Mask_Count / Cloud_Mask_Count
      else
         Nonconfident_Cloud_Mask_Fraction = Missing_Value_Real4 
      end if
  
   end subroutine COMPUTE_CLOUD_MASK_PERFORMANCE_METRICS


!==============================================================================
!
! remote sensing reflectance - for ocean applications 
!
! Reference: (Smyth, Tyrrell and Tarrant, 2004 GRL, 31)
!
! note, 0.63 rayleigh optical passed to this from global memory
!==============================================================================
 real elemental function REMOTE_SENSING_REFLECTANCE ( &
                                             atmos_corrected_063_reflectance, &
                                             atmos_corrected_086_reflectance, &
                                             air_mass, &
                                             solar_zenith)
 

  real, intent(in):: atmos_corrected_063_reflectance
  real, intent(in):: atmos_corrected_086_reflectance
  real, intent(in):: air_mass
  real, intent(in):: solar_zenith
  real, parameter:: solar_zenith_max_threshold = 89.0

  if (solar_zenith < solar_zenith_max_threshold) then
     remote_sensing_reflectance = (atmos_corrected_063_reflectance -  &
                                   atmos_corrected_086_reflectance) /  &
                                   exp( -0.5*Solar_Rtm%Tau_Ray(1) * air_mass)
  else
     remote_sensing_reflectance = Missing_Value_Real4
  endif

 end function REMOTE_SENSING_REFLECTANCE

!==============================================================================
!
! normalized difference vegetation index - for land applications
!
! note, missing value passed from global memory
!==============================================================================
 real elemental function NORMALIZED_DifFERENCE_VEGETATION_INDEX ( &
                                             atmos_corrected_063_reflectance, &
                                             atmos_corrected_086_reflectance, &
                                             solar_zenith)
 

  real, intent(in):: atmos_corrected_063_reflectance
  real, intent(in):: atmos_corrected_086_reflectance
  real, intent(in):: solar_zenith
  real, parameter:: solar_zenith_max_threshold = 89.0

  if (solar_zenith < solar_zenith_max_threshold) then
    normalized_difference_vegetation_index =  &
                   (atmos_corrected_086_reflectance - atmos_corrected_063_reflectance) /  &
                   (atmos_corrected_086_reflectance + atmos_corrected_063_reflectance) 
  else
    normalized_difference_vegetation_index = Missing_Value_Real4
  endif

 end function NORMALIZED_DifFERENCE_VEGETATION_INDEX
!==============================================================================
!
! normalized difference snow index - for land applications
!
! note, missing value passed from global memory
!==============================================================================
 real elemental function NORMALIZED_DifFERENCE_SNOW_INDEX ( &
                                             atmos_corrected_063_reflectance, &
                                             atmos_corrected_160_reflectance, &
                                             solar_zenith)
 

  real, intent(in):: atmos_corrected_063_reflectance
  real, intent(in):: atmos_corrected_160_reflectance
  real, intent(in):: solar_zenith
  real, parameter:: solar_zenith_max_threshold = 89.0

  if (solar_zenith < solar_zenith_max_threshold) then
    normalized_difference_snow_index =  &
                   (atmos_corrected_063_reflectance - atmos_corrected_160_reflectance) /  &
                   (atmos_corrected_063_reflectance + atmos_corrected_160_reflectance) 
  else
    normalized_difference_snow_index = Missing_Value_Real4
  endif

 end function NORMALIZED_DifFERENCE_SNOW_INDEX
!==============================================================================
!
! normalized difference desert index - for land applications
!
! note, missing value passed from global memory
!==============================================================================
 real elemental function NORMALIZED_DifFERENCE_DESERT_INDEX ( &
                                             atmos_corrected_063_reflectance, &
                                             atmos_corrected_160_reflectance, &
                                             solar_zenith)
 

  real, intent(in):: atmos_corrected_063_reflectance
  real, intent(in):: atmos_corrected_160_reflectance
  real, intent(in):: solar_zenith
  real, parameter:: solar_zenith_max_threshold = 89.0

  if (solar_zenith < solar_zenith_max_threshold) then
    normalized_difference_desert_index =  &
                   (atmos_corrected_160_reflectance - atmos_corrected_063_reflectance) /  &
                   (atmos_corrected_160_reflectance + atmos_corrected_063_reflectance) 
  else
    normalized_difference_desert_index = Missing_Value_Real4
  endif

 end function NORMALIZED_DifFERENCE_DESERT_INDEX

!==============================================================================
! A routine that call functions to populate some simple surface parameters
!==============================================================================
 subroutine SURFACE_REMOTE_SENSING(Line_Idx_Min,Line_Idx_Max,ABI_Use_104um_Flag)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Line_Idx_max
  logical, intent(in):: ABI_Use_104um_Flag

  if (ABI_Use_104um_Flag) then
    if (Sensor%Chan_On_Flag_Default(38) == sym%YES) then
      call COMPUTE_TSFC(Line_Idx_Min,Line_Idx_Max,ABI_Use_104um_Flag)
    endif
  else
    if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      call COMPUTE_TSFC(Line_Idx_Min,Line_Idx_Max,ABI_Use_104um_Flag)
    endif
  endif

  if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. &
      Sensor%Chan_On_Flag_Default(6) == sym%YES) then

       Ndsi_Toa = NORMALIZED_DifFERENCE_SNOW_INDEX(  &
                             ch(1)%Ref_Toa, &
                             ch(6)%Ref_Toa, &
                             Geo%Solzen)

       Ndsi_Sfc = NORMALIZED_DifFERENCE_SNOW_INDEX(  &
                             ch(1)%Ref_Sfc, &
                             ch(6)%Ref_Sfc, &
                             Geo%Solzen)

  endif

  if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. &
      Sensor%Chan_On_Flag_Default(2) == sym%YES) then

       Ndvi_Toa = NORMALIZED_DifFERENCE_VEGETATION_INDEX(  &
                             ch(1)%Ref_Toa, &
                             ch(2)%Ref_Toa, &
                             Geo%Solzen)

       Ndvi_Sfc = NORMALIZED_DifFERENCE_VEGETATION_INDEX(  &
                             ch(1)%Ref_Sfc, &
                             ch(2)%Ref_Sfc, &
                             Geo%Solzen)

       Ndvi_Sfc_White_Sky = NORMALIZED_DifFERENCE_VEGETATION_INDEX(  &
                             ch(1)%Sfc_Ref_White_Sky, &
                             ch(2)%Sfc_Ref_White_Sky, &
                             Geo%Solzen)

       Rsr = REMOTE_SENSING_REFLECTANCE( ch(1)%Ref_Sfc, &
                                         ch(2)%Ref_Sfc, &
                                         Geo%Airmass,       &
                                         Geo%Solzen)
  endif

  if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. &
      Sensor%Chan_On_Flag_Default(6) == sym%YES) then

       Nddi_Toa = NORMALIZED_DifFERENCE_DESERT_INDEX(  &
                             ch(1)%Ref_Toa, &
                             ch(6)%Ref_Toa, &
                             Geo%Solzen)
  endif

 end subroutine SURFACE_REMOTE_SENSING
!==============================================================================
! COMPUTE_DESERT_MASK_FOR_CLOUD_DETECTION
!
! Purpose:  Feed the location of snow-free deserts to the cloud mask
!
!
!==============================================================================
integer(kind=int1) elemental function DESERT_MASK_FOR_CLOUD_DETECTION( &
                                       Emiss_Sfc_375um,                 &
                                       Lat,                             &
                                       Snow,                            &
                                       Surface_Type)

   real(kind=real4), intent(in):: Emiss_Sfc_375um
   real(kind=real4), intent(in):: Lat
   integer(kind=int1), intent(in):: Snow
   integer(kind=int1), intent(in):: Surface_Type

   Desert_Mask_For_Cloud_Detection = 0

   if ( Snow == sym%NO_SNOW .and.  &
        Surface_Type > 0 .and.         &
        Emiss_Sfc_375um < 0.93  .and.  &
        abs(Lat) < 60.0 .and. &
        ((Surface_Type == sym%open_SHRUBS_SFC) .or.  &
         (Surface_Type == sym%CLOSED_SHRUBS_SFC) .or. &
         (Surface_Type == sym%GRASSES_SFC) .or.  &
         (Surface_Type == sym%BARE_SFC)) ) then

         Desert_Mask_For_Cloud_Detection = 1

   endif

end function DESERT_MASK_FOR_CLOUD_DETECTION
!==============================================================================
! COMPUTE_CITY_MASK_FOR_CLOUD_DETECTION
!
! Purpose:  Feed the location of cities to the cloud mask
!
!
!==============================================================================
integer(kind=int1) elemental function CITY_MASK_FOR_CLOUD_DETECTION( &
                                       Rad_Lunar,                     &
                                       Surface_Type)

   real(kind=real4), intent(in):: Rad_Lunar
   integer(kind=int1), intent(in):: Surface_Type

   real, parameter:: Radiance_Lunar_City_Thresh = 2.5e-08

   City_Mask_For_Cloud_Detection = 0

   !--- use surface type information
   if (Surface_Type == sym%URBAN_SFC) then
      City_Mask_For_Cloud_Detection = 0
   endif

   !----------------------------------------------------------------------------
   !--- if lunar radiance is available, assume large values are cities or other
   !--- surface surfaces of light that we treat as cities
   !--- note, need to check if lunar radiance is available.
   !----------------------------------------------------------------------------
   if (allocated(ch(44)%Rad_Toa)) then
     if (Rad_Lunar > Radiance_Lunar_City_Thresh) then
       City_Mask_For_Cloud_Detection = 1
     endif
   endif


end function CITY_MASK_FOR_CLOUD_DETECTION

!-----------------------------------------------------------------------------
! EUMETCAST Fire detection algorithm
!
!This implements the "Current Operational Algorithm" described in:
!TOWARDS AN IMPROVED ACTIVE FIRE MONITORING PRODUCT FOR MSG SATELLITES
!Sauli Joro, Olivier Samain, Ahmet Yildirim, Leo van de Berg, Hans Joachim Lutz
!EUMETSAT, Am Kavalleriesand 31, Darmstadt, Germany
!-----------------------------------------------------------------------------
  integer elemental function FIRE_TEST (T11,T375,T11_std,T375_std,Solzen)

     real, intent(in):: T11
     real, intent(in):: T375
     real, intent(in):: T11_Std
     real, intent(in):: T375_Std
     real, intent(in):: Solzen

     real :: Bt_375um_Eumet_Fire_Thresh
     real :: Bt_Diff_Eumet_Fire_Thresh
     real :: Stddev_11um_Eumet_Fire_Thresh
     real :: Stddev_375um_Eumet_Fire_Thresh

     !---- EUMETCAST fire detection parameters
     real, parameter :: EUMETCAST_FIRE_DAY_SOLZEN_THRESH = 70.0
     real, parameter :: EUMETCAST_FIRE_NIGHT_SOLZEN_THRESH = 90.0

     real, parameter :: BT_375UM_EUMET_FIRE_DAY_THRESH = 310.0
     real, parameter :: BT_DifF_EUMET_FIRE_DAY_THRESH = 8.0
     real, parameter :: STDDEV_11UM_EUMET_FIRE_DAY_THRESH = 1.0
     real, parameter :: STDDEV_375UM_EUMET_FIRE_DAY_THRESH = 4.0

     real, parameter :: BT_375UM_EUMET_FIRE_NIGHT_THRESH = 290.0
     real, parameter :: BT_DifF_EUMET_FIRE_NIGHT_THRESH = 0.0
     real, parameter :: STDDEV_11UM_EUMET_FIRE_NIGHT_THRESH = 1.0
     real, parameter :: STDDEV_375UM_EUMET_FIRE_NIGHT_THRESH = 4.0

     !--- initialize
     Fire_Test = 0

     
     !--- check if all needed data are non-missing
     if (T375 /= Missing_Value_Real4 .and. &
         T375_Std /= Missing_Value_Real4 .and. &
         T11 /= Missing_Value_Real4 .and. &
         T11_Std /= Missing_Value_Real4) then

         !Day
         if (Solzen < EumetCAST_Fire_Day_Solzen_Thresh) then
            Bt_375um_Eumet_Fire_Thresh = Bt_375um_Eumet_Fire_day_Thresh
            Bt_Diff_Eumet_Fire_Thresh = Bt_Diff_Eumet_Fire_day_Thresh
            Stddev_11um_Eumet_Fire_Thresh = Stddev_11um_Eumet_Fire_Day_Thresh
            Stddev_375um_Eumet_Fire_Thresh = Stddev_375um_Eumet_Fire_Day_Thresh
         endif

         !Night
         if (Solzen > EumetCAST_Fire_Night_Solzen_Thresh) then
            Bt_375um_Eumet_Fire_Thresh = Bt_375um_Eumet_Fire_Night_Thresh
            Bt_Diff_Eumet_Fire_Thresh = Bt_Diff_Eumet_Fire_Night_Thresh
            Stddev_11um_Eumet_Fire_Thresh = Stddev_11um_Eumet_Fire_Night_Thresh
            Stddev_375um_Eumet_Fire_Thresh = Stddev_375um_Eumet_Fire_Night_Thresh
         endif

         !Twilight
         if ((Solzen >= EumetCAST_Fire_Day_Solzen_Thresh) .and. &
             (Solzen <= EumetCAST_Fire_Night_Solzen_Thresh)) then

             !linear fit day -> night
             Bt_375um_Eumet_Fire_Thresh = ((-1.0)* Solzen) + 380.0
             Bt_Diff_Eumet_Fire_Thresh = ((-0.4)* Solzen) + 36.0

             !These two don't change, but 
             Stddev_11um_Eumet_Fire_Thresh = STDDEV_11UM_EUMET_FIRE_NIGHT_THRESH
             Stddev_375um_Eumet_Fire_Thresh = STDDEV_375UM_EUMET_FIRE_NIGHT_THRESH

         endif

       ! All of these conditions need to be met
       if ((T375 > Bt_375um_Eumet_Fire_Thresh) .and. &
           ((T375 - T11) > Bt_Diff_Eumet_Fire_Thresh) .and. &
           (T375_Std > Stddev_375um_Eumet_Fire_Thresh) .and. &
           (T11_Std < Stddev_11um_Eumet_Fire_Thresh)) then
         Fire_Test = 1
       endif

     endif

  end function FIRE_TEST
!-----------------------------------------------------------
! make VIIRS look like MODIS interms of spatial sampling
! 
! this is accomplished by averaging in the along scan direction
! appropriately to acheive a pixel size that grows with
! scan angle - as would be the case with MODIS or AVHRR
!
! this does not resample the data so there is over sampling
!
!-----------------------------------------------------------
subroutine VIIRS_TO_MODIS()

  real, dimension(3):: weight
  real, parameter, dimension(3):: weight_1_0 = (/0.00,1.00,0.00/)
  real, parameter, dimension(3):: weight_1_5 = (/0.25,1.00,0.25/)
  real, parameter, dimension(3):: weight_3_0 = (/1.00,1.00,1.00/)

  integer:: Elem_Idx, Line_Idx, Chan_Idx
  integer:: Number_of_Lines, Number_of_Elements
  integer:: i1, i2


  !--- set image size
  Number_of_Elements = Image%Number_Of_Elements
  Number_of_Lines = Image%Number_Of_Lines_Read_This_Segment

  !--- loop over channels
  do Chan_Idx = 1, NChan_Clavrx

    !--- check if channel is on
    if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle

    !--- solar reflectances
    if (ch(Chan_Idx)%Obs_Type == SOLAR_OBS_TYPE .or. &
        ch(Chan_Idx)%Obs_Type == MIXED_OBS_TYPE) then

      do Line_Idx = 1, Number_of_Lines
         do Elem_Idx = 2, Number_of_Elements-1
  
           i1 = Elem_Idx - 1
           i2 = Elem_Idx + 1

           !--- pick weighting scheme
           weight = weight_1_5
           if (Elem_Idx <= 640 .or. Elem_Idx >= 2561) weight = weight_3_0
           if (Elem_Idx >= 1009 .and. Elem_Idx <= 2192) weight = weight_1_0
  
           !-- apply weighting
           Temp_Pix_Array_1(Elem_Idx,Line_Idx) = sum(weight*Ch(Chan_Idx)%Ref_Toa(i1:i2,Line_Idx)) / sum(weight)

         enddo
      enddo 

      !--- copy back
      ch(Chan_Idx)%Ref_Toa = Temp_Pix_Array_1

    endif

    !--- thermal
    if (ch(Chan_Idx)%Obs_Type == THERMAL_OBS_TYPE .or. &
        ch(Chan_Idx)%Obs_Type == MIXED_OBS_TYPE) then

      do Line_Idx = 1, Number_of_Lines
         do Elem_Idx = 2, Number_of_Elements-1

           i1 = Elem_Idx - 1
           i2 = Elem_Idx + 1
  
           !--- pick weighting scheme
           weight = weight_1_5
           if (Elem_Idx <= 640 .or. Elem_Idx >= 2561) weight = weight_3_0
           if (Elem_Idx >= 1009 .and. Elem_Idx <= 2192) weight = weight_1_0
  
           !-- apply weighting
           Temp_Pix_Array_1(Elem_Idx,Line_Idx) = sum(weight*Ch(Chan_Idx)%Rad_Toa(i1:i2,Line_Idx)) / sum(weight)
  
         enddo
      enddo 

      !--- copy back
      ch(Chan_Idx)%Rad_Toa = Temp_Pix_Array_1

      !--- compute BT
      do Line_Idx = 1, Number_of_Lines
         do Elem_Idx = 2, Number_of_Elements-1
            ch(Chan_Idx)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(Chan_Idx,ch(Chan_Idx)%Rad_Toa(Elem_Idx,Line_Idx))
         enddo
      enddo 

    endif

  enddo
  
end subroutine VIIRS_TO_MODIS
!-----------------------------------------------------------
! Modify the MODIS AUX Phase and make a Type
!-----------------------------------------------------------
subroutine MODIFY_AUX_CLOUD_TYPE

   integer:: i, j

   !--- cloud type
   line_loop_type: do i = 1, Image%Number_Of_Elements
      elem_loop_type: do  j = 1, Image%Number_Of_Lines_Read_This_Segment

       if (.not. ABI_Use_104um_Flag) then
          call COMPUTE_TYPE_FROM_PHASE(CLDMASK%Cld_Mask(i,j), &
                                     Cld_Phase_Aux(i,j), &
                                     Cld_Phase_Uncertainty(i,j), &
                                     Zc_Opaque_Cloud(i,j), &
                                     Tc_Opaque_Cloud(i,j), &
                                     ch(31)%Emiss_Tropo(i,j), &
                                     Sensor%Chan_On_Flag_Default(31), &
                                     Sensor%Chan_On_Flag_Default(32), &
                                     Sensor%Chan_On_Flag_Default(33), &
                                     Beta_11um_12um_Tropo_Rtm(i,j), &
                                     Beta_11um_133um_Tropo_Rtm(i,j), &
                                     Cld_Type_Aux(i,j))
       else

          call COMPUTE_TYPE_FROM_PHASE(CLDMASK%Cld_Mask(i,j), &
                                     Cld_Phase_Aux(i,j), &
                                     Cld_Phase_Uncertainty(i,j), &
                                     Zc_Opaque_Cloud(i,j), &
                                     Tc_Opaque_Cloud(i,j), &
                                     ch(38)%Emiss_Tropo(i,j), &
                                     Sensor%Chan_On_Flag_Default(38), &
                                     Sensor%Chan_On_Flag_Default(32), &
                                     Sensor%Chan_On_Flag_Default(33), &
                                     Beta_104um_12um_Tropo_Rtm(i,j), &
                                     Beta_104um_133um_Tropo_Rtm(i,j), &
                                     Cld_Type_Aux(i,j))
       endif

      enddo elem_loop_type
   enddo line_loop_type

end subroutine MODIFY_AUX_CLOUD_TYPE

!-----------------------------------------------------------
! end of MODULE
!-----------------------------------------------------------
end module PIXEL_ROUTINES_MOD
