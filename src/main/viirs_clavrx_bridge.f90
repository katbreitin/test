! $Header: https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/trunk/main_src/viirs_clavrx_bridge.f90 4009 2020-09-23 22:27:34Z dbotambekov $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: viirs_clavrx_bridge.f90 (src)
!       viirs_clavrx_bridge (program)
!
! PURPOSE: VIIRS reader bridge to clavr-x
!
! DESCRIPTION: This module deals with the clav-x global value world with the more 
!              hidden viirs data world.  
!
! AUTHORS:
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
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
! REVISION HISTORY:   created      March 2013 (AW)
!            4 Dec 2013 -  added mapping table in comments (AW) 
!                       -  change mapping of modis 3 and 4
!
!            12 January 2014 AW: add lunar rela zimuth and scattering angle 
!                                computation for global pixel_common_mod variable
!
!
!             21 May 2019   : removed unused subroutine compute_iband_statistics and julian
!
! NOTES:
!  VIIRS  MODIS(CLAVRX)  Wavelength
!           mapping
!    M1    -     8    -    0.412
!    M2    -     9    -    0.445
!    M3    -     3    -    0.488
!    M4    -     4    -    0.555
!    M5    -     1    -    0.672
!    M6    -    15    -    0.746
!    M7    -     2    -    0.865
!    M8    -     5    -    1.240
!    M9    -    26    -    1.378
!    M10   -     6    -    1.610
!    M11   -     7    -    2.250
!    M12   -    20    -    3.700
!    M13   -    22    -    4.050
!    M14   -    29    -    8.550
!    M15   -    31    -   10.763
!    M16   -    32    -   12.013
!    I1    -    39    -    0.640
!    I2    -    40    -    0.865
!    I3    -    41    -    1.610
!    I4    -    42    -    3.740
!    I5    -    43    -   11.450
!    DNB   -    44    -    0.700
!
!--------------------------------------------------------------------------------------

module VIIRS_CLAVRX_BRIDGE 

   use Pixel_Common_Mod , only : &
        Image &
      , Sensor &
      , Geo &
      , Nav &
      , Ancil_Data_Dir & 
      , Use_Aux_Flag &
      , Cloud_Mask_Aux_Read_Flag &
      , Cloud_Type_Aux_Read_Flag &
      , CLDMASK &
      , Cld_Type_Aux &
      , Cld_Phase_Aux &
      , Gap_Pixel_Mask &
      , Ch &
      , Use_Iband &
      , X_Sample_Offset &
      , Y_Sample_Offset &
      , Temp_Pix_Array_1
      

   use CONSTANTS_MOD, only: &
      Int4 &
    , Sym &
    , Missing_Value_Real4
      
   use CLAVRX_MESSAGE_MOD

   private

   public:: READ_VIIRS_DATA, READ_VIIRS_DATE_TIME, GET_NUMBER_OF_SCANS_FROM_VIIRS_BRIDGE, READ_VIIRS_INSTR_CONSTANTS
  
!--------------------------------------------------------------------------------------------------------------------
! 
!--------------------------------------------------------------------------------------------------------------------
contains
   
   subroutine READ_VIIRS_DATA (Segment_Number, File_Gmtco_Base, Error_Out)
      use VIIRS_READ_MOD , only : &
           VIIRS_DATA_CONFIG &
           , VIIRS_DATA_OUT &
           , GET_VIIRS_DATA 
   
      use PLANCK_MOD
      use VIEWING_GEOMETRY_MOD, only: &
           GLINT_ANGLE &
           , SCATTERING_ANGLE &
           , RELATIVE_AZIMUTH
  
      use CALIBRATION_CONSTANTS_MOD, only: &
            Planck_Nu, &
            VIIRS_Correction_Factor
      
      implicit none
      
      integer , intent(in) :: segment_number
      character(len=*), intent(in) :: file_gmtco_base
      integer(kind=int4), intent(out) :: error_out
      
      type ( viirs_data_config )  :: v_conf
      type ( viirs_data_out )  :: out
      integer :: modis_chn_list (16)
      integer :: modis_chn_list_iband (5)
      logical :: is_mband_on (16)
      logical :: is_iband_on (5)
      integer :: i_mband , i_iband
      integer :: y_start , c_seg_lines , c_seg_lines_iband
      integer :: i
      integer :: modis_chn
  
      error_out = 0
      ! - mapping modis to viirs
      !                 041 044 048 055  068  074   085 124 138 160 225  375  405  855  108  120
      !                 M1  M2   M3   M4  M5   M6   M7  M8  M9  M10 M11  M12  M13  M14  M15  M16  
      modis_chn_list = [ 8 , 9 , 3 , 4 , 1 , 15 , 2 , 5 , 26 , 6 , 7 , 20 , 22 , 29 , 31 , 32 ]
      modis_chn_list_iband = [ 39 , 40 , 41 , 42 , 43 ]
      is_mband_on = Sensor%Chan_On_Flag_Default (modis_chn_list) == sym%YES
      is_iband_on = Sensor%Chan_On_Flag_Default (modis_chn_list_iband) == sym%YES
      
      y_start = ( segment_number -1 ) * Image%Number_Of_Lines_Per_Segment + 1
      c_seg_lines = min (  y_start + Image%Number_Of_Lines_Per_Segment -1 , Image%Number_Of_Lines )  - y_start  + 1
      
      ! - configure viirs interface
      v_conf % chan_on_rfl_mband = is_mband_on
      v_conf % chan_on_iband = is_iband_on
      v_conf % chan_on_dnb = Sensor%Chan_On_Flag_Default(44) == sym%YES
      !v_conf % viirs_cloud_mask_on = use_aux_flag /= sym%NO_AUX
      !v_conf % viirs_cloud_type_on = use_aux_flag /= sym%NO_AUX
      v_conf % viirs_cloud_mask_on = .false.
      v_conf % viirs_cloud_type_on = .false.
      if (use_aux_flag == sym%USE_AUX .or. use_aux_flag == sym%READ_BUT_DO_NOT_USE_AUX) v_conf % viirs_cloud_mask_on = .true.
      if (use_aux_flag == sym%USE_AUX .or. use_aux_flag == sym%READ_BUT_DO_NOT_USE_AUX) v_conf % viirs_cloud_type_on = .true.
      
      v_conf % offset = [ 1 , y_start]
      v_conf % count = [ Image%Number_Of_Elements  , c_seg_lines  ]
      v_conf % dir_1b = trim(Image%Level1b_Path)
      
      v_conf % Ancil_Data_Dir = trim(Ancil_Data_Dir)
      v_conf % file_gmtco_base =  trim(file_gmtco_base)

      v_conf % Nu_List = 0.0
      v_conf % Nu_List(12:16) = [Planck_Nu(20) , Planck_Nu(22) ,  &
                                Planck_Nu(29) , Planck_Nu(31) , Planck_Nu(32)]

      ! - read the data 
      call get_viirs_data ( v_conf, out )

      ! - output to clavrx global variables
      ! geo
      nav % lat_1b(:,1:c_seg_lines)    = out % geo % lat
      nav % lon_1b(:,1:c_seg_lines)    = out % geo % lon
      image % scan_time_ms(1:c_seg_lines)   = int( out % geo % scan_time)
      geo % sataz(:,1:c_seg_lines)     = out % geo % sataz
      geo % satzen(:,1:c_seg_lines)    = out % geo % satzen
      geo % solaz (:,1:c_seg_lines)    = out % geo % solaz 
      geo % solzen (:,1:c_seg_lines)   = out % geo % solzen 

      if(out % file_exists % gdnbo_file_exists) then 
        geo % moon_phase_angle = out % geo % Moon_Phase_Angle
        geo % Moon_Illum_Frac = out % geo % moon_illum_frac
        geo % lunfrac = out % geo % moon_illum_frac
      endif

      ! rel azimuths  - these are all global variables
      geo % relaz = RELATIVE_AZIMUTH ( geo % solaz , geo % sataz )

      !--- compute the glint zenith angle
      geo % glintzen = GLINT_ANGLE ( geo % solzen , geo % satzen , geo % relaz )

      !--- compute the scattering angle
      geo % scatangle = SCATTERING_ANGLE ( geo % solzen , geo % satzen , geo % relaz )

      ! gap
      gap_pixel_mask( : ,1:c_seg_lines) = 0
      where ( out % gap % mask )
         gap_pixel_mask( : ,1:c_seg_lines) = 1
      end where 
      
      ! - m-bands
      do i_mband = 1 , 16
         modis_chn = modis_chn_list (i_mband)
         if ( .not. out % mband ( i_mband ) % is_read ) then
            sensor % chan_on_flag_per_line (modis_chn ,1:c_seg_lines) = sym % no 
            cycle   
         end if
         
         if ( .not. is_mband_on(i_mband) .or. (size(out % mband (i_mband) % ref) < 1 &
              .and. size(out % mband (i_mband) % rad) < 1) ) cycle
         
         if ( i_mband <= 11 ) then
            ch(modis_chn) % Ref_Toa ( : ,1:c_seg_lines)  =  viirs_correction_factor(i_mband) * out % mband (i_mband) % ref   
         end if
         if ( i_mband >= 12 ) then
            
            ch ( modis_chn)  % Rad_Toa( : ,1:c_seg_lines) =   out % mband (i_mband) % rad
            call COMPUTE_BT_ARRAY ( ch(modis_chn)%bt_toa , ch(modis_chn)%rad_toa , &
                                modis_chn , missing_value_real4 )
         end if
         
      end do
      
      ! - i-bands
      do  i_iband = 1 , 5
         
            if ( .not. out % file_exists % svi_file_exists (i_iband)) then
                 ! - switch off chan_on in CLAVR-x if file is not there..
               Sensor%Chan_On_Flag_Default ( modis_chn_list_iband ) = sym % NO
               sensor % chan_on_flag_per_line (modis_chn_list_iband (i_iband) ,1:c_seg_lines) = sym % NO
               cycle
            end if
            
            
           if ( .not. out % iband ( i_iband ) % is_read ) then
            sensor % chan_on_flag_per_line (modis_chn_list_iband (i_iband) ,1:c_seg_lines) = sym % no 
          
            cycle   
         end if
         
         if ( .not. is_iband_on(i_iband) .or. (size(out % iband (i_iband) % ref) < 1 &
              .and. size(out % iband (i_iband) % bt) < 1) ) then    
              cycle
         end if
         
      end do 

      !---- put a comment here about what you are doing with dnb
!      if (Sensor%Chan_On_Flag_Default(44) == sym%YES .and. size(out % dnb_mgrid % rad) > 1) then
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES .and. allocated( out % dnb_mgrid % rad )) then
         ch(44)%rad_toa( : ,1:c_seg_lines)  = out % dnb_mgrid % rad
         geo % lunzen( : ,1:c_seg_lines) = out % geo % lunzen
         geo % lunaz( : ,1:c_seg_lines) = out % geo % lunaz
         ch(44)%ref_toa( : ,1:c_seg_lines) = out % dnb_mgrid % ref
         geo % lunrelaz( : ,1:c_seg_lines) = RELATIVE_AZIMUTH ( geo % lunaz( : ,1:c_seg_lines) &
                                                             , geo % sataz( : ,1:c_seg_lines) )
          
         !--- compute the scattering angle
         geo % scatangle_lunar( : ,1:c_seg_lines) = SCATTERING_ANGLE( geo % lunzen( : ,1:c_seg_lines) &
                                                , geo % satzen( : ,1:c_seg_lines) &
                                                , geo % lunrelaz( : ,1:c_seg_lines) )
         geo % glintzen_lunar( : ,1:c_seg_lines) = GLINT_ANGLE( geo % lunzen( : ,1:c_seg_lines) &
                                             , geo % satzen( : ,1:c_seg_lines) &
                                             , geo % lunrelaz( : ,1:c_seg_lines) )                                       
      end if
      
      ! -global variables which has to be set
      Image%Number_Of_Lines_Read_This_Segment = c_seg_lines
      do i = 1, Image%Number_Of_Lines_Per_Segment
         Image%Scan_Number(i) = y_start + i - 1
      end do
      
      !- ascending  (global varaible )
      nav % ascend = 0  
      do i = 1 , Image%Number_Of_Lines_Read_This_Segment - 1
         if ( nav % lat_1b(Image%Number_Of_Elements / 2 , i + 1) <= nav % lat_1b( Image%Number_Of_Elements / 2 , i ) ) nav % ascend ( i )  = 1
      end do
      ! --- fix for the last line Denis B.
      if ( nav % lat_1b(Image%Number_Of_Elements / 2 , Image%Number_Of_Lines_Read_This_Segment) <= &
           nav % lat_1b(Image%Number_Of_Elements / 2 , Image%Number_Of_Lines_Read_This_Segment - 1) ) &
            nav % ascend ( Image%Number_Of_Lines_Read_This_Segment ) = 1

      !---  statistics I-Band on M-band grid
      if (out%iband(1)%is_read) then
         Temp_Pix_Array_1 = Missing_Value_Real4
         call COMPUTE_IBAND_MIN_MAX(out%iband(1)%ref, Ch(1)%Ref_Toa_Min_Sub, Ch(1)%Ref_Toa_Max_Sub, &
                                    Temp_Pix_Array_1,X_Sample_Offset,Y_Sample_Offset)
         if (Use_Iband) then
            ch(1)%Ref_Toa = Temp_Pix_Array_1
         endif
      endif
      if (out%iband(2)%is_read) then
         Temp_Pix_Array_1 = Missing_Value_Real4
         call COMPUTE_IBAND_MIN_MAX(out%iband(2)%ref, Ch(2)%Ref_Toa_Min_Sub, Ch(2)%Ref_Toa_Max_Sub, &
                                    Temp_Pix_Array_1,X_Sample_Offset,Y_Sample_Offset)
         if (Use_Iband) then
            ch(2)%Ref_Toa = Temp_Pix_Array_1
         endif
      endif
      if (out%iband(3)%is_read) then
         Temp_Pix_Array_1 = Missing_Value_Real4
         call COMPUTE_IBAND_MIN_MAX(out%iband(3)%ref, Ch(6)%Ref_Toa_Min_Sub, Ch(6)%Ref_Toa_Max_Sub, &
                                    Temp_Pix_Array_1,X_Sample_Offset,Y_Sample_Offset)
         if (Use_Iband) then
            ch(6)%Ref_Toa = Temp_Pix_Array_1
         endif
      endif
      if (out%iband(4)%is_read) then
         Temp_Pix_Array_1 = Missing_Value_Real4
         call COMPUTE_IBAND_MIN_MAX(out%iband(4)%bt, Ch(20)%Bt_Toa_Min_Sub, Ch(20)%Bt_Toa_Max_Sub, &
                                    Temp_Pix_Array_1,X_Sample_Offset,Y_Sample_Offset)
         if (Use_Iband) then
            ch(20)%Bt_Toa = Temp_Pix_Array_1
         endif
      endif
      if (out%iband(5)%is_read) then
         Temp_Pix_Array_1 = Missing_Value_Real4
         call COMPUTE_IBAND_MIN_MAX(out%iband(5)%bt, Ch(31)%Bt_Toa_Min_Sub, Ch(31)%Bt_Toa_Max_Sub, &
                                    Temp_Pix_Array_1,X_Sample_Offset,Y_Sample_Offset)
         if (Use_Iband) then
            ch(31)%Bt_Toa = Temp_Pix_Array_1
         endif
      endif
  
      ! --- save aux cloud mask
      if ( out % prd % is_mask_read ) then
         CLDMASK%Cld_Mask_Aux( : ,1 : c_seg_lines ) = out % prd % cld_mask
         Cloud_Mask_Aux_Read_Flag = Sym % YES
      else
         Cloud_Mask_Aux_Read_Flag = Sym % NO
      end if   

      ! --- save aux cloud type and phase
      if ( out % prd % is_phase_read ) then
         cld_type_aux( : ,1 : c_seg_lines) = out % prd % cld_type
         cld_phase_aux( : ,1 : c_seg_lines ) = out % prd % cld_phase
         Cloud_Type_Aux_Read_Flag = Sym % YES
      else
         Cloud_Type_Aux_Read_Flag = Sym % NO
      end if   
      
      ! --- deallocate all
      call out % dealloc ()

   end subroutine READ_VIIRS_DATA



   !----------------------------------------------------------------
   ! - iband has full file dimension of 6400 x1536
   ! - mband 3200 x 768
   !  - output of min_val ... is 3200 768
   !----------------------------------------------------------------
   subroutine COMPUTE_IBAND_MIN_MAX ( iband_array , out_min_val, out_max_val, out_sample_val, x_offset_in, y_offset_in)
      implicit none
      real, dimension(:,:) , intent(in) :: iband_array
      integer, intent(in):: x_offset_in, y_offset_in
      real, dimension(:,:) , intent(out)  :: out_min_val, out_max_val, out_sample_val
      real, dimension(2,2) :: small_iband
      
      integer :: im , jm
      integer , dimension(2) ::  dim_m
      integer , dimension(2) ::  dim_i
      integer :: iband_x0, iband_x1 ,  iband_y0, iband_y1
      integer:: x_offset, y_offset
     
      dim_m = shape ( out_min_val )
      dim_i = shape ( iband_array )
      
      out_min_val = Missing_Value_Real4
      out_max_val = Missing_Value_Real4

      !--- to be consistent with geostationary static nav convenction,
      !--- X_Offset_In and Y_Offset_In should 0 or 1
      !--- ensure sampling offsets are between 1 and 2
      X_Offset = X_Offset_In + 1
      Y_Offset = Y_Offset_In + 1
      X_Offset = min(2,max(1,X_Offset))
      Y_Offset = min(2,max(1,Y_Offset))


      do im = 1, dim_m(1)
        
         iband_x0 = (im-1)*2 + 1
         iband_x1 = iband_x0 + 1

         !--- check of exceeding bounds (should skip)
         if (iband_x0 > dim_i(1)) cycle
         if (iband_x1 > dim_i(1)) cycle

         do jm = 1, dim_m(2)
            iband_y0 = (jm-1)*2 + 1
            iband_y1 = iband_y0 + 1
   
            !--- check of exceeding bounds (should skip)
            if (iband_y0 > dim_i(2)) cycle
            if (iband_y1 > dim_i(2)) cycle
            small_iband = iband_array ( iband_x0 :  iband_x1 ,  iband_y0 :  iband_y1 )
            if ( minval ( small_iband ) > 0 ) then 
               out_min_val ( im, jm ) = minval ( small_iband )
               out_max_val ( im, jm ) = maxval ( small_iband )
               out_sample_val ( im, jm ) = small_iband(X_Offset,Y_Offset)
            end if
         end do

      end do

   end subroutine COMPUTE_IBAND_MIN_MAX
   
   !----------------------------------------------------------------
   ! read the VIIRS constants into memory
   !-----------------------------------------------------------------
   subroutine READ_VIIRS_INSTR_CONSTANTS(Instr_Const_file)
      use calibration_constants_mod
      use file_utils, only: get_lun
      
      implicit none
 
      character(len=*), intent(in):: Instr_Const_file
      integer:: ios0, erstat
      integer:: Instr_Const_lun
      character(len=20):: header

      Instr_Const_lun = GET_LUN()

      open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)
      call mesg ("opening "//trim(Instr_Const_file), level = verb_lev % VERBOSE) 
      erstat = 0
      if (ios0 /= 0) then
         erstat = 19
         print *,  " VIIRS_CLAVRX_BRIDGE.f90: Error opening VIIRS constants file, ios0 = ", ios0
         stop 19
      end if

      read(unit=Instr_Const_lun,fmt="(a3)") Sat_Name
      read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
      read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
      read(unit=Instr_Const_lun,fmt=*) header
      read(unit=Instr_Const_lun,fmt=*) Planck_A1(20), Planck_A2(20),Planck_Nu(20)
      read(unit=Instr_Const_lun,fmt=*) Planck_A1(22), Planck_A2(22),Planck_Nu(22)
      read(unit=Instr_Const_lun,fmt=*) Planck_A1(29), Planck_A2(29),Planck_Nu(29)
      read(unit=Instr_Const_lun,fmt=*) Planck_A1(31), Planck_A2(31),Planck_Nu(31)
      read(unit=Instr_Const_lun,fmt=*) Planck_A1(32), Planck_A2(32),Planck_Nu(32)
      read(unit=Instr_Const_lun,fmt=*) Planck_A1(42), Planck_A2(42),Planck_Nu(42)
      read(unit=Instr_Const_lun,fmt=*) Planck_A1(43), Planck_A2(43),Planck_Nu(43)
      read(unit=Instr_Const_lun,fmt=*) header
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(1)
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(2)
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(3)
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(4)
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(5)
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(6)
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(7)
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(8)
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(9)
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(10)
      read(unit=Instr_Const_lun,fmt=*) VIIRS_Correction_Factor(11)
      close(unit=Instr_Const_lun)
  
      !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
      Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

   end subroutine READ_VIIRS_INSTR_CONSTANTS

   !-----------------------------------------------------------------------------------------
   !  Get information from VIIRS 
   !-----------------------------------------------------------------------------------------
  
   subroutine READ_VIIRS_DATE_TIME ( Path, Infile &
                , Year , Doy , Start_Time , End_Time , Orbit , Orbit_Identifier &
                , End_Year, End_Doy )

      use VIIRS_READ_MOD, only : &                                                                                                                    
            READ_VIIRS_DATE_TIME_ATT
          

      ! Get the date & time from the file's name
      implicit none
      
      character(len=*), intent(in) :: Path
      character(len=*), intent(in) :: Infile
      integer, intent(inout) , optional :: Year
      integer, intent(inout)  , optional:: Doy    !day of year
      integer, intent(inout) , optional :: Start_Time  !millisec
      integer, intent(inout)  , optional:: End_Time    !millisec
      integer, intent(inout)  , optional:: Orbit
      character(38), intent(inout) , optional :: Orbit_Identifier
      integer , intent(inout) , optional :: End_Year
      integer, intent(inout)  , optional:: End_Doy    !day of year
  

      !--- call READ_VIIRS_DATE_TIME_ATT from module
      call READ_VIIRS_DATE_TIME_ATT (Path, Infile &
                , Year , Doy , Start_Time , End_Time , Orbit , Orbit_Identifier &
                , End_Year, End_Doy )

   end subroutine READ_VIIRS_DATE_TIME 

!---------------------------------------------------------------------------------
!  subroutine GET_NUMBER_OF_SCANS_FROM_VIIRS_BRIDGE ( Infile , Number_Of_Viirs_Lines , Error_Out )
!  it's asking to read number of scans in the viirs_read_mod
!---------------------------------------------------------------------------------
   SUBROUTINE GET_NUMBER_OF_SCANS_FROM_VIIRS_BRIDGE ( Infile , Number_Of_Viirs_Lines , Error_Out )
      use viirs_read_mod , only : &
          READ_NUMBER_OF_SCANS_FROM_VIIRS
   
      CHARACTER(Len=*), INTENT(IN) :: Infile  
      INTEGER(kind=int4), INTENT(OUT) :: Error_Out
      INTEGER(KIND=INT4), INTENT(OUT):: Number_of_Viirs_Lines

      error_out = 0
      call READ_NUMBER_OF_SCANS_FROM_VIIRS ( Infile , Number_of_Viirs_Lines , Error_Out )
     
   END SUBROUTINE GET_NUMBER_OF_SCANS_FROM_VIIRS_BRIDGE
   

end module viirs_clavrx_bridge

