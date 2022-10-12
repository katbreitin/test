! $Id: dncomp_interface_def_mod.f90 182 2018-02-15 13:17:01Z awalther $
module dncomp_interface_def_mod

   ! --  Module works as interface between CLAVR-x and DCOMP
   !
   !
   implicit none
   ! - parameters
   integer, parameter :: REAL4 = selected_real_kind(6,37)
   integer, parameter :: INT1 = selected_int_kind(1)   
   integer, parameter :: INT2 = selected_int_kind(2)
   real ( kind = real4), parameter :: MISSING_REAL4 = -999.
   
   ! - CLAVRX uses 44 MODIS/VIIRS channels
   integer , parameter :: N_CHN = 44   
   
   ! - object for 2D real4 arrays
   type d2_real4_type
      logical :: is_set
      integer :: xdim
      integer :: ydim
      real ( kind = real4 ) , dimension(:,:) , allocatable  :: d  
     
   end type d2_real4_type
   
   ! - object for 2D int1 arrays
   type d2_int1_type
      logical :: is_set
      integer :: xdim
      integer :: ydim
      integer ( kind = int1 ) , dimension(:,:) , allocatable  :: d  
     
   end type d2_int1_type
   
   ! - object for 2D int2 arrays
   type d2_int2_type
      logical :: is_set
      integer :: xdim
      integer :: ydim
      integer ( kind = int2 ) , dimension(:,:) , allocatable  :: d  
   end type d2_int2_type
   
   ! - object for 2D logical arrays
   type d2_flag_type
      logical :: is_set
      logical, dimension(:,:), allocatable :: d
     
   end type d2_flag_type
   
   ! - object for gas coeff values
   type gas_coeff_type
      logical :: is_set
      real ( kind = real4 )  :: d (3)  
   end type  gas_coeff_type
   
   type channel_type
    real, pointer :: rad(:,:) => NULL()
    real, pointer :: refl(:,:) => NULL()
    real, pointer :: alb_sfc(:,:) => NULL()
    real, pointer :: alb_sfc_dark_sky(:,:) => NULL()
    real, pointer :: emiss_sfc(:,:) => NULL()
    real, pointer :: trans_ac_nadir(:,:) => NULL()
    real, pointer :: rad_clear_sky_toc(:,:) => NULL()
    real, pointer :: rad_clear_sky_toa(:,:) => NULL()
   end type channel_type
   
   
   ! - main dcomp input type
   type dncomp_in_type
   
      ! - configure
     
      integer :: mode
      character ( len = 1024) :: lut_path
      integer :: sensor_wmo_id
      integer :: n_channel 
      logical :: is_channel_on (N_CHN)
      
      ! - satellite input channel
      type (channel_type), allocatable :: chn(:)
      
      real, pointer :: sat(:,:) => NULL()
      real, pointer :: sol(:,:) => NULL()
      real, pointer :: azi(:,:) => NULL()
      real, pointer :: zen_lunar(:,:) => NULL()
      real, pointer :: azi_lunar(:,:) => NULL()
      ! - cloud products
      integer(int1), pointer  :: cloud_mask(:,:) => NULL()
      integer(int1), pointer :: cloud_type(:,:) => NULL()
      real, pointer  :: cloud_hgt(:,:) => NULL()
      real, pointer  :: cloud_temp(:,:) => NULL()
      real, pointer  :: cloud_press(:,:) => NULL()
      real, pointer  :: tau_acha(:,:) => NULL()
      
      ! - flags
      logical, pointer  :: is_land (:,:) => NULL()
      logical, pointer  :: is_valid(:,:) => NULL()
      
      ! - surface
      real, pointer :: press_sfc(:,:) => NULL()
      integer(int1), pointer   :: snow_class (:,:) => NULL()
      
      ! - atmosphere
      real, pointer :: ozone_nwp(:,:) => NULL()
      real, pointer :: tpw_ac(:,:) => NULL()
    
      ! - coeffecients,params
      real :: sun_earth_dist
      TYPE ( gas_coeff_type )  :: gas_coeff ( N_CHN)
      real :: solar_irradiance(N_CHN)
      
   contains
      !final :: in_destructor 
      procedure :: check_input
         
   end type dncomp_in_type
   
   
   ! - DCOMP/NLCOMP output
      
   type dncomp_out_type
      type ( d2_real4_type) :: cod  
      ! real, pointer :: cod(:,:) => NULL()
      type ( d2_real4_type) :: cps
      type ( d2_real4_type) :: cod_unc
      type ( d2_real4_type) :: ref_unc
      type ( d2_real4_type) :: cld_trn_sol
      type ( d2_real4_type) :: cld_trn_obs
      type ( d2_real4_type) :: cld_alb
      type ( d2_real4_type) :: cld_sph_alb
      type ( d2_int2_type) :: quality
      type ( d2_int2_type) :: info
      type ( d2_real4_type) :: iwp  
      type ( d2_real4_type) :: lwp
      type ( d2_real4_type ) :: rain_probability
      type ( d2_real4_type ) :: rain_rate
      type ( d2_real4_type ) :: refl_vis_max

      character ( len = 200 ) :: version
      real :: successrate
      integer :: nr_clouds
      integer :: nr_obs
      integer :: nr_success_cod
      integer :: nr_success_cps
      
    ! contains
    !  final :: out_destructor
   end type dncomp_out_type 
   
   ! - Enumerated cloud type
   type et_cloud_type_type
      integer(kind=int1) :: FIRST = 0
      integer(kind=int1) :: CLEAR = 0
      integer(kind=int1) :: PROB_CLEAR = 1
      integer(kind=int1) :: FOG = 2
      integer(kind=int1) :: WATER = 3
      integer(kind=int1) :: SUPERCOOLED = 4
      integer(kind=int1) :: MIXED = 5
      integer(kind=int1) :: OPAQUE_ICE = 6
      integer(kind=int1) :: TICE = 6
      integer(kind=int1) :: CIRRUS = 7
      integer(kind=int1) :: OVERLAP = 8
      integer(kind=int1) :: OVERSHOOTING = 9
      integer(kind=int1) :: UNKNOWN = 10
      integer(kind=int1) :: DUST = 11
      integer(kind=int1) :: SMOKE = 12
      integer(kind=int1) :: FIRE = 13  
      integer(kind=int1) :: LAST = 13
   end type
   
   type ( et_cloud_type_type ) , protected :: EM_cloud_type
   
   ! - Enumerated clod mask
   type et_cloud_mask_type
      integer(kind=int1) :: LAST = 3
      integer(kind=int1) :: CLOUDY = 3
      integer(kind=int1) :: PROB_CLOUDY = 2
      integer(kind=int1) :: PROB_CLEAR = 1
      integer(kind=int1) :: CLEAR = 0
      integer(kind=int1) :: FIRST = 0
   end type
   
   type ( et_cloud_mask_type ) , protected :: EM_cloud_mask
   
   ! - Enumerated snow/sea ice class
   type  et_snow_class_type
      integer(kind=int1) :: FIRST = 1
      integer(kind=int1) :: NO_SNOW = 1
      integer(kind=int1) :: SEA_ICE = 2
      integer(kind=int1) :: SNOW = 3
      integer(kind=int1) :: LAST = 3
   end type 
   
   type (  et_snow_class_type ) , protected :: EM_snow_class 
   
   interface alloc_dncomp
      module procedure &
         alloc_it_d2_real, alloc_it_d2_int, alloc_it_d2_log
   
   end interface 
   
   interface dncomp_out_type
      module procedure new_output
   end interface dncomp_out_type
   
   interface dncomp_in_type
      module procedure new_input
   end interface dncomp_in_type
   
    
contains

   !  Constructs new DCOMP / NLCOMP output derived type
   !
   !      
   function new_output ( dim_1, dim_2 )
      implicit none
      integer , intent(in) :: dim_1, dim_2
      type ( dncomp_out_type ) :: new_output
      
      allocate ( new_output % cod % d         ( dim_1 , dim_2))
      allocate ( new_output % cps % d         ( dim_1 , dim_2))
      allocate ( new_output % cod_unc % d     ( dim_1 , dim_2))
      allocate ( new_output % ref_unc % d     ( dim_1 , dim_2))
      allocate ( new_output % cld_trn_sol % d ( dim_1 , dim_2))
      allocate ( new_output % cld_trn_obs % d ( dim_1 , dim_2))
      allocate ( new_output % cld_alb % d     ( dim_1 , dim_2))
      allocate ( new_output % cld_sph_alb % d ( dim_1 , dim_2))
      allocate ( new_output % info % d        ( dim_1 , dim_2))
      allocate ( new_output % quality % d     ( dim_1 , dim_2))
      allocate ( new_output % lwp % d         ( dim_1 , dim_2))
      allocate ( new_output % iwp % d         ( dim_1 , dim_2))
      allocate ( new_output % rain_probability % d ( dim_1 , dim_2))
      allocate ( new_output % rain_rate % d ( dim_1 , dim_2))
      allocate ( new_output % refl_vis_max % d ( dim_1 , dim_2))

      
      new_output % cod % d           =  MISSING_REAL4 
      new_output % cps % d           =  MISSING_REAL4 
      new_output % cod_unc % d       =  MISSING_REAL4
      new_output % ref_unc % d       =  MISSING_REAL4 
      new_output % cld_trn_sol % d   =  MISSING_REAL4  
      new_output % cld_trn_obs % d   =  MISSING_REAL4   
      new_output % cld_alb % d       =  MISSING_REAL4  
      new_output % cld_sph_alb % d   =  MISSING_REAL4 
      new_output % lwp % d          =  MISSING_REAL4
      new_output % lwp % d          =  MISSING_REAL4
      new_output % rain_probability % d       =  -999  
      new_output % rain_rate % d   =  MISSING_REAL4
      new_output % refl_vis_max % d   =  MISSING_REAL4
   end function new_output
   
   !  Constructs new DNCOMP input derived type
   !
   !
   function new_input (  chan_on )
     
      logical, intent(in) :: chan_on ( :)
      type ( dncomp_in_type ) :: new_input
      integer :: n_chn
      integer :: idx_chn
      
         ! === ALLOCATION
      n_chn = size ( chan_on)      
       
    
      new_input % n_channel = n_chn
      
      allocate ( new_input % chn(n_chn))
         
      new_input % is_channel_on = .false.
      do idx_chn = 1 , n_chn
        
         if ( .not. chan_on (idx_chn) ) cycle
         
         new_input % is_channel_on(idx_chn) = .true.
         
         
         
      end do
    
      
   end function new_input
   
   
   !
   !
   !
   subroutine check_input ( this , debug_mode_in)
      class(dncomp_in_type) :: this 
      integer, intent(in), optional :: debug_mode_in
      integer :: debug_mode
      integer :: n_pixels
      debug_mode = 0
     
      if (present(debug_mode_in)) debug_mode = debug_mode_in
      
      if ( debug_mode .le. 0 ) return
      
      n_pixels = size ( this % sol )
      
      print*,'Test input ranges '
      print*,'Solar zenith range valid for ', 100.* float(count ( this % sol  .lt. 65 )) / n_pixels
      print*,'Sensor zenith range valid for ', 100.* float(count ( this % sat  .lt. 65 )) / n_pixels 
      print*,'Azimuth valid range ', 100.* float(count ( this % azi  .ge. 0 .and.  this % azi  .le.  180. )) / n_pixels 
      
      
   
   end subroutine check_input
   
   !
   !
   !
   subroutine out_destructor ( this )
      type ( dncomp_out_type) :: this
      
      if ( allocated ( this % cod % d ) ) deallocate (  this % cod % d )
      if ( allocated ( this % cps % d )) deallocate ( this % cps % d)
      if ( allocated ( this % cod_unc % d ) ) deallocate (  this % cod_unc % d )
      if ( allocated ( this % ref_unc % d )) deallocate ( this % ref_unc % d)
      if ( allocated ( this % cld_trn_sol % d ) ) deallocate (  this % cld_trn_sol % d ) 
      if ( allocated ( this % cld_trn_obs % d )) deallocate ( this % cld_trn_obs % d)
      if ( allocated ( this % cld_alb % d ) ) deallocate (  this % cld_alb % d )
      if ( allocated ( this % cld_sph_alb % d )) deallocate ( this % cld_sph_alb % d)
      if ( allocated ( this % quality % d ) ) deallocate (  this % quality % d )
      if ( allocated ( this % info % d )) deallocate ( this % info % d)
      if ( allocated ( this % iwp % d ) ) deallocate (  this % iwp % d )
      if ( allocated ( this % lwp % d )) deallocate ( this % lwp % d)
      if ( allocated ( this % rain_probability % d ) ) deallocate (  this % rain_probability % d )
      if ( allocated ( this % rain_rate % d )) deallocate ( this % rain_rate % d)
      if ( allocated ( this % refl_vis_max % d )) deallocate ( this % refl_vis_max % d)
      
   end subroutine out_destructor
   
   
   !  --  allocation routines   
   subroutine alloc_it_d2_real ( str , xdim , ydim )
      type ( d2_real4_type ) :: str
      integer :: xdim , ydim
      integer :: alloc_stat = 0
        
      allocate ( str % d ( xdim,  ydim) , stat = alloc_stat)
      if ( alloc_stat /= 0 ) then
         print*,'alloc error'            
      end if      
   end subroutine alloc_it_d2_real
   
   !  --  allocation routine
   subroutine alloc_it_d2_int ( str , xdim , ydim )
      type ( d2_int1_type ) :: str
      integer :: xdim , ydim
      integer :: alloc_stat = 0
         
      allocate ( str % d ( xdim,  ydim) , stat = alloc_stat)
      if ( alloc_stat /= 0 ) then
         print*,'alloc error'            
      end if      
   end subroutine alloc_it_d2_int
   
   !  --  allocation routine   
   subroutine alloc_it_d2_log ( str , xdim , ydim )
      type ( d2_flag_type ) :: str
      integer :: xdim , ydim
      integer :: alloc_stat = 0
         
      allocate ( str % d ( xdim,  ydim) , stat = alloc_stat)
      if ( alloc_stat /= 0 ) then
         print*,'alloc error'            
      end if      
   end subroutine alloc_it_d2_log
   
   
   
   
   !  Finalization tool for dcomp_input
   !
   !
   subroutine in_destructor ( this )
      type ( dncomp_in_type ) :: this
      integer :: i
      
     
      do i = 1, N_CHN
       
        
        if ( associated ( this % chn(i) % refl))  this % chn(i) % refl => null()
        if ( associated ( this % chn(i) % rad))  this % chn(i) % rad => null()
        if ( associated ( this % chn(i) % alb_sfc))  this % chn(i) % alb_sfc => null()
        if ( associated ( this % chn(i) % alb_sfc_dark_sky))  this % chn(i) % alb_sfc_dark_sky => null()
        if ( associated ( this % chn(i) % emiss_sfc))  this % chn(i) % emiss_sfc => null()
        if ( associated ( this % chn(i) % trans_ac_nadir))  this % chn(i) % trans_ac_nadir => null()
        if ( associated ( this % chn(i) % rad_clear_sky_toa))  this % chn(i) % rad_clear_sky_toa => null()
        if ( associated ( this % chn(i) % rad_clear_sky_toc))  this % chn(i) % rad_clear_sky_toc => null()
         
      end do

      !if ( allocated (this % sol % d) ) deallocate ( this % sol  % d )
      this % sat => null()
      this % sol => null()
      this % azi => null()
      this % zen_lunar => null()
      this % azi_lunar => null()
      
      this % cloud_mask => null()
      this %  cloud_type  => null()
      this %  cloud_hgt  => null()
      this %  cloud_temp  => null()
      this %  cloud_press => null()
      this %  tau_acha   => null()
   
      this % snow_class => null()
      this % ozone_nwp => null()
      this % tpw_ac => null()
      this % press_sfc => null()
      
      this % is_land => null()
      this % is_valid => null()
      
     
      
   end subroutine in_destructor
   

end module dncomp_interface_def_mod
