module cx_sfc_emissivity_mod
  
 use PIXEL_COMMON_MOD, only: Nav, Geo &
  , Ch, image, emiss_land_option, emiss_sea_option, Ancil_Data_Dir, Sensor, sfc
 
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
      , GET_UMD_EMISS &
      , GET_PIXEL_SFC_EMISS_FROM_SFC_TYPE &
      , GET_PIXEL_SFC_REFL_FROM_SFC_TYPE
 
  use CONSTANTS_MOD, only: sym, MIXED_OBS_TYPE &
      , nchan_clavrx &
      , THERMAL_OBS_TYPE
  
  use CX_SEA_IR_EMISS_MOD, only: GET_SEGMENT_SEA_IR_EMISS
  
  use CLAVRX_MESSAGE_MOD, only: MESG, VERB_LEV
  
  implicit none
  private
  public:: cx_sfc_emiss_populate_land
  public:: cx_sfc_emiss_populate_sea
  public:: cx_sfc_emiss_correct_for_sfctype
  
enum, bind(C)
  enumerator :: ETsfc_emiss_land_option_UMD = 0
  enumerator :: ETsfc_emiss_land_option_RTTOV =1
  enumerator :: ETsfc_emiss_land_option_SEEBOR =2
  enumerator :: ETsfc_emiss_land_option_CRTM =3
  
end enum


enum, bind(C)
  enumerator :: ETsfc_emiss_sea_option_UMD = 0
  enumerator :: ETsfc_emiss_sea_option_RTTOV =1
  enumerator :: ETsfc_emiss_sea_option_LUT =2
  enumerator :: ETsfc_emiss_sea_option_CRTM =3
  
end enum

logical :: first_run = .true.




contains

! - The only task for this subroutine is to populate 
!   ch (CH_IDX) % sfc_emiss
!    and this will be only changed here!!
!
!
  subroutine cx_sfc_emiss_populate_land
#define STRINGIFY(x) x
    implicit none
    character (len=:), allocatable :: rttov_path 
    character(len = 256) :: path_sfc
    integer :: chan_idx
    logical, dimension(:,:), allocatable :: mask 
    real, dimension(:,:), allocatable :: emiss_dum
    integer :: size_array(2)
    
    size_array = shape(Geo%Space_Mask)
    allocate(mask(size_array(1),size_array(2)))
    
    
    mask = Geo%Space_Mask
    
   print*,'EMISS LAND OPTION: ',emiss_land_option
  
  
    select case(emiss_land_option)
  
    case(ETsfc_emiss_land_option_UMD)
      if (first_run) call MESG('SFC EMISS Land UMD',  level = verb_lev % DEFAULT)
       call GET_UMD_EMISS(mask)
    case(ETsfc_emiss_land_option_RTTOV)
      if (first_run) call MESG('SFC EMISS Land RTTOV',  level = verb_lev % DEFAULT) 
#ifdef LIBRTTOV   
      rttov_path = trim(Ancil_Data_Dir) // "static/rttov/"
      call GET_RTTOV_EMISS(Nav%Lat, Nav%Lon, mask , rttov_path)  
#else
      print*, 'RTTOV emissivity selected but not compiled with RTTOV'
      stop
#endif
    case(ETsfc_emiss_land_option_SEEBOR)
      if (first_run) call MESG('SFC EMISS Land SEEBOR',  level = verb_lev % DEFAULT) 
      
      allocate (emiss_dum(size_array(1),size_array(2)))
      
      path_sfc = trim(Ancil_Data_Dir)//"static/sfc_data"
    
      do Chan_Idx = 20, Nchan_Clavrx
        if (ch(Chan_Idx)%Obs_Type /= MIXED_OBS_TYPE .and. &
              ch(Chan_Idx)%Obs_Type /= THERMAL_OBS_TYPE) cycle
          !--- force channel 20 read used for desert definition
        if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES .or. (chan_idx == 20)) then
            call READ_SEEBOR_EMISS( path_sfc , chan_idx, Nav%Lat, Nav%Lon &
              , mask, image % time_start % month, emiss_dum)   
              
            where( .not. mask )
               ch(chan_idx)%Sfc_Emiss = emiss_dum
            end where
     
        end if
        
        
      end do  
      
     
    
    case(ETsfc_emiss_land_option_CRTM)
      print*,ETsfc_emiss_land_option_CRTM
      stop 'CRTM Land Emiss is not yet installed'
  
    case default
      call MESG(' Surface emissivity option is wrongly set',  level = verb_lev % DEFAULT) 
    end select
     
 
    first_run = .false.
  
  end subroutine cx_sfc_emiss_populate_land 
  
  
    subroutine cx_sfc_emiss_populate_sea
#define STRINGIFY(x) x
    implicit none
    character (len=:), allocatable :: rttov_path 
    character(len = 256) :: path_sfc
    integer :: chan_idx
    logical, dimension(:,:), allocatable :: mask 
    integer :: size_array(2)
    
    size_array = shape(Geo%Space_Mask)
    allocate(mask(size_array(1),size_array(2)))
    
    mask = .false.
    where ( Sfc%Sfc_Type .NE. 0)
      mask = .true.
    end where 
   print*,'EMISS SEA OPTION: ',emiss_sea_option
    
  
  
    select case(emiss_sea_option)
  
    case(ETsfc_emiss_sea_option_UMD)
      if (first_run) call MESG('SFC EMISS sea UMD',  level = verb_lev % DEFAULT)
       if (emiss_land_option .NE. ETsfc_emiss_land_option_UMD ) call GET_UMD_EMISS(mask)
    case(ETsfc_emiss_sea_option_RTTOV)
      if (first_run) call MESG('SFC EMISS Sea RTTOV',  level = verb_lev % DEFAULT) 
#ifdef LIBRTTOV   
      rttov_path = trim(Ancil_Data_Dir) // "static/rttov/"
      if (emiss_land_option .NE. ETsfc_emiss_land_option_RTTOV ) call GET_RTTOV_EMISS(Nav%Lat, Nav%Lon, mask, rttov_path)  
#else
      print*, 'RTTOV emissivity selected but not compiled with RTTOV'
      stop
#endif
    case(ETsfc_emiss_sea_option_LUT)
      if (first_run) call MESG('SFC EMISS Sea SEA IR LUT',  level = verb_lev % DEFAULT) 
      
      call GET_SEGMENT_SEA_IR_EMISS()
    
    case(ETsfc_emiss_sea_option_CRTM)
      print*,ETsfc_emiss_land_option_CRTM
      stop 'CRTM Land Emiss is not yet installed'
  
    case default
      call MESG(' Surface emissivity option is wrongly set',  level = verb_lev % DEFAULT) 
    end select
     
 
    first_run = .false.
  
  end subroutine cx_sfc_emiss_populate_sea
  
  
  
  ! This subroutine is called in process_clavr-x after snow class is computed.
  !
  !
  subroutine cx_sfc_emiss_correct_for_sfctype
    integer :: Line_Idx_Min_Segment
    
    !if ( sfc_emiss_option .NE. ETsfc_emiss_use_option_RTTOV) then
     
      !  now file all non-land pixels
      Line_Idx_Min_Segment = 1
      call GET_PIXEL_SFC_EMISS_FROM_SFC_TYPE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment) 
 
      ! - check if we want better sea
      if (emiss_sea_option == 2) then
        call GET_SEGMENT_SEA_IR_EMISS()
      end if
    
    
    !end if
    
  end subroutine cx_sfc_emiss_correct_for_sfctype


subroutine sfc_emiss_get



end  subroutine sfc_emiss_get



end module cx_sfc_emissivity_mod
