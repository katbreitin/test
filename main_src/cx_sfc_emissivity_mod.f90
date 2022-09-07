module cx_sfc_emissivity_mod
  
 use PIXEL_COMMON_MOD, only: Nav, Geo &
  , Ch, image, Use_Sea_Ir_Emiss, Ancil_Data_Dir, Sensor, sfc_emiss_option 
 
  use SFC_EMISS, only: &
       read_seebor_emiss
#ifdef LIBRTTOV 
 use CX_RTTOV_SFC_EMISS, only: &
      init_rttov_emiss &
      , get_rttov_emiss &
      , destroy_rttov_emiss
#endif  
 use SURFACE_PROPERTIES_MOD,only: &
      COMPUTE_BINARY_LAND_COAST_MASKS &
      , setup_umd_props &
      , GET_PIXEL_SFC_EMISS_FROM_SFC_TYPE &
      , GET_PIXEL_SFC_REFL_FROM_SFC_TYPE
 
  use CONSTANTS_MOD, only: sym, MIXED_OBS_TYPE &
      , nchan_clavrx &
      , THERMAL_OBS_TYPE
  
  use CX_SEA_IR_EMISS_MOD, only: GET_SEGMENT_SEA_IR_EMISS
  
  use CLAVRX_MESSAGE_MOD, only: MESG, VERB_LEV
  
  implicit none
  private
  public:: cx_sfc_emiss_populate_ch
  public:: cx_sfc_emiss_correct_for_sfctype
  
enum, bind(C)
  enumerator :: ETsfc_emiss_use_option_UMD = 0
  enumerator :: ETsfc_emiss_use_option_RTTOV =1
  enumerator :: ETsfc_emiss_use_option_SEEBOR =2
  enumerator :: ETsfc_emiss_use_option_CRTM =3
  
end enum

logical :: first_run = .true.




contains

! - The only task for this subroutine is to populate 
!   ch (CH_IDX) % sfc_emiss
!    and this will be only changed here!!
!
!
  subroutine cx_sfc_emiss_populate_ch
#define STRINGIFY(x) x
    implicit none
    character (len=:), allocatable :: rttov_path 
    character(len = 256) :: path_sfc
    integer :: chan_idx
  
    
  
  
    select case(sfc_emiss_option)
  
    case(ETsfc_emiss_use_option_UMD)
      if (first_run) call MESG('SFC EMISS UMD',  level = verb_lev % DEFAULT)
  
    case(ETsfc_emiss_use_option_RTTOV)
      if (first_run) call MESG('SFC EMISS RTTOV',  level = verb_lev % DEFAULT) 
#ifdef LIBRTTOV   
      rttov_path = trim(Ancil_Data_Dir) // "static/rttov/"
      call GET_RTTOV_EMISS(Nav%Lat, Nav%Lon, Geo%Space_Mask, rttov_path)  
#else
      print*, 'RTTOV emissivity selected but not compiled with RTTOV'
      stop
#endif
    case(ETsfc_emiss_use_option_SEEBOR)
      if (first_run) call MESG('SFC EMISS SEEBOR',  level = verb_lev % DEFAULT) 
      path_sfc = trim(Ancil_Data_Dir)//"static/sfc_data"
    
      do Chan_Idx = 20, Nchan_Clavrx
        if (ch(Chan_Idx)%Obs_Type /= MIXED_OBS_TYPE .and. &
              ch(Chan_Idx)%Obs_Type /= THERMAL_OBS_TYPE) cycle
          !--- force channel 20 read used for desert definition
        if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES .or. (chan_idx == 20)) then
            call READ_SEEBOR_EMISS( path_sfc , chan_idx, Nav%Lat, Nav%Lon &
              , Geo%Space_Mask, image % time_start % month, ch(chan_idx)%Sfc_Emiss)   
        end if
      end do  
    
    case(ETsfc_emiss_use_option_CRTM)
      print*,ETsfc_emiss_use_option_CRTM
      stop 'CRTM is not yet installed'
  
    case default
      call MESG(' Surface emissivity option is wrongly set',  level = verb_lev % DEFAULT) 
    end select
   
    first_run = .false.
  
  end subroutine cx_sfc_emiss_populate_ch 
  
  ! This subroutine is called in process_clavr-x after snow class is computed.
  !
  !
  subroutine cx_sfc_emiss_correct_for_sfctype
    integer :: Line_Idx_Min_Segment
    
      !  now file all non-land pixels
      Line_Idx_Min_Segment = 1
      call GET_PIXEL_SFC_EMISS_FROM_SFC_TYPE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment) 
 
      ! - check if we want better sea
      if (Use_Sea_Ir_Emiss == sym%YES) then
        call GET_SEGMENT_SEA_IR_EMISS()
      end if
    
    
    !end if
    
  end subroutine cx_sfc_emiss_correct_for_sfctype


subroutine sfc_emiss_get



end  subroutine sfc_emiss_get



end module cx_sfc_emissivity_mod
