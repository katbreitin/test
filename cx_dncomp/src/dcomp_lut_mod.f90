! $Id: dcomp_lut_mod.f90 189 2018-03-09 18:49:39Z awalther $
!
! HISTORY: 06/05/2014: change of filename 
!
module dcomp_lut_mod
   
   use dcomp_math_tools_mod, only:&
      dcomp_interpolation_weight &
      , interpolate_2d
      
   use file_tools, only: &
      file_test
   use dcomp_lut_sds_mod, only: &
      read_hdf_dcomp_data_ems &
      , read_hdf_dcomp_data_rfl &
      , read_hdf_dcomp_dims

   implicit none
 
   private
   
   integer, parameter :: NUM_PHASE = 2
   integer, parameter :: NUM_CHN = 45
   
   integer , parameter :: NUM_SOL = 45
   integer , parameter :: NUM_SAT = 45
   integer , parameter :: NUM_AZI = 45
   integer , parameter :: NUM_COD = 29
   integer , parameter :: NUM_REF = 9
   
   type lut_data_type
      logical :: is_set
      logical :: has_sol = .false.
      logical :: has_ems = .false.
      character (len = 1020 ) :: file
      character (len = 1020 ) :: file_ems
      real, allocatable :: cld_sph_alb ( : , : )
      real, allocatable :: cld_trn ( : , : , : )
      real, allocatable :: cld_alb ( : , : , : )
      real, allocatable :: cld_refl( : , : , : , : , : )
      real, allocatable :: cld_ems ( : , : , : )
      real, allocatable :: cld_trn_ems ( : , : , : )
      
   contains
      procedure :: read_hdf => lut_data__read_hdf
      procedure :: alloc => lut_data__alloc
      procedure :: dealloc => lut_data__dealloc
      
   end type lut_data_type
   
   type lut_chn_type
      logical :: is_set
      type ( lut_data_type ) , allocatable  :: phase(:)
   end type lut_chn_type
   
   type lut_dim_type
      logical :: is_set
      real, allocatable  :: sat_zen ( : )
      real, allocatable  :: sol_zen ( : )
      real, allocatable  :: rel_azi ( : )
      real , allocatable :: cod ( : )
      real, allocatable  :: cps ( : )
      integer :: n_sat_zen
      integer :: n_sol_zen
      integer :: n_rel_azi
      integer :: n_cod
      integer :: n_cps

       
   end type lut_dim_type
   
   type , public :: lut_type
      logical :: is_set
      type ( lut_chn_type ) , allocatable :: channel (:)
      type ( lut_dim_type ) ::  dims
      character ( len = 1020) :: lut_path
      character ( len = 20 ) :: sensor = 'not_set'
      integer :: pos_sat, pos_sol, pos_azi
      
   contains
      procedure :: initialize => lut__initialize
      procedure :: set_angles => lut__set_angles
      procedure :: getProperty => lut__getProperty
      procedure :: get_data => lut__get_data
      procedure :: thick_cloud_rfl => lut_data__thick_cloud_rfl
      procedure :: init_dims => lut__init_dims
      procedure , private :: set_filename => lut__set_filename
      procedure :: clear_lut_memory => lut__clear_lut_memory
      procedure :: populate_all_at_once => lut__populate_all_at_once
      
   end type lut_type
   
   ! - only one channel can be stored
   type ( lut_type ) , public ,save , target :: lut_obj
   
   
   type , public :: lut_output
   
      real :: refl
      real :: trn_sol
      real :: trn_sat
      real :: albsph
      real :: alb
      real :: ems
      real :: trn_ems
      
      real :: dRefl_dcps      , dRefl_dcod
      real :: dtrans_sol_dcps , dtrans_sol_dcod
      real :: dTrans_sat_dcod , dTrans_sat_dcps
      real :: dsph_alb_dcod   , dSph_alb_dcps
      real :: dalb_dcod   , dalb_dcps
      real :: dEms_dcps      , dEms_dcod
      real :: dtrnEms_dcps      , dtrnEms_dcod
   
   
   end type lut_output
   

   
   real, save :: sat_m
   real, save :: sol_m
   real, save :: azi_m
      
contains

   ! ------------------------------------------------------------
   !  This is to populate the tables outside dcomp if wished.
   ! ------------------------------------------------------------
   
   subroutine lut__populate_all_at_once (self , sensor , ancil_path)
      implicit none
      class ( lut_type ) , target :: self
      character ( len = * ) , intent(in) :: sensor
      character ( len = * ) , intent(in) :: ancil_path
      integer :: idx_chn, idx_phase
      type (lut_data_type), pointer :: data_loc => null()
      
     
      call self % initialize ( trim(sensor), trim(ancil_path))
       
      do idx_chn = 1 , NUM_CHN
         do idx_phase = 1, 2

            data_loc => self % channel ( idx_chn ) % phase ( idx_phase)
            if ( .not. data_loc % is_set ) then
               
               call data_loc % read_hdf
               data_loc % is_set  = .true.
            end if
         end do
      end do
      
      self % sensor = trim(sensor)
     
   end subroutine lut__populate_all_at_once
   
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut__set_filename ( self)
      class ( lut_type ) :: self
      character ( len = 3 ) , dimension(30) :: chan_string ='no'
      logical , dimension ( NUM_CHN ) :: has_ems_table = .false.
      logical , dimension ( NUM_CHN ) :: has_sol_table = .false.
      
      integer :: i_chn , i_phase
      character ( len = 3 ) , dimension(2)   :: phase_string = [ 'wat',  'ice' ]
      integer :: n_channels = 45
      character ( len =1020) :: sensor_identifier
      
      ! mapping sensor channel emis yes/no
      
      sensor_identifier = trim ( self % sensor )
     !print*,'==> ', self % sensor
      sensor_block: select case ( trim(self % sensor))
         case ('Meteosat-8','Meteosat-9','Meteosat-10','Meteosat-11') sensor_block
            has_sol_table(1) = .true.
            has_sol_table(2) = .false.
            has_sol_table(6) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '1'
            chan_string(2) = '2'
            chan_string(6) = '3'
            chan_string(20) = '4'
      
         case ('TIROS-N','NOAA-05','NOAA-06','NOAA-07','NOAA-08','NOAA-09', 'NOAA-10','NOAA-11', &
            'NOAA-12', 'NOAA-14')  sensor_block
            has_sol_table(1) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '1'
            chan_string(20) = '3b'
        
         case ('NOAA-15','NOAA-16', 'NOAA-17','NOAA-18','NOAA-19','METOP-A','METOP-B','METOP-C')  sensor_block
            has_sol_table(1) = .true.
            has_sol_table(6) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '1'
            chan_string(6) = '3a'
            chan_string(20) = '3b'
      
         case ('GOES-08','GOES-09','GOES-10', 'GOES-11' , 'GOES-12' , 'GOES-13',  'GOES-14', &
            'GOES-15','COMS-1'  )   sensor_block
            has_sol_table(1) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '1'
            chan_string(20) = '2'
      
         case ('MODIS-AQUA', 'MODIS-TERRA')    sensor_block
            has_sol_table(1) = .true.
            has_sol_table(2) = .false.
            has_sol_table(5:7) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '1'
            !chan_string(2) = '2'
            chan_string(5) = '5'
            chan_string(6) = '6'
            chan_string(7) = '7'
            chan_string(20) = '20'
            !--->sensor_identifier = trim(self % lut_path) //'MODIS'
            sensor_identifier = 'MODIS'
         
         case('GOES-16','GOES-17')  sensor_block
            has_sol_table(1) = .true.
            has_sol_table(6) = .true.
            has_sol_table(7) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '2'
            chan_string(6) = '5'
            chan_string(7) = '6'
            chan_string(20) = '7'
         
         case('ABI') sensor_block
            has_sol_table(1) = .true.
            has_sol_table(6) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '2'
            chan_string(6) = '5'
            chan_string(20) = '7'
      
         case('AHI') sensor_block
            has_sol_table(1) = .true.
            has_sol_table(6) = .true.
            has_sol_table(7) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '3'
            chan_string(6) = '5'
            chan_string(7) = '6'
            chan_string(20) = '7'
         
            
         case ('AATSR')   sensor_block
            has_sol_table(1) = .true.
            has_sol_table(6) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '1'
            chan_string(6) = '6'
            chan_string(20) = '20'
            
            
        
             
         case ('VIIRS','NOAA-20')   sensor_block
            has_sol_table(1) = .true.
            has_sol_table(5) = .true.
            has_sol_table(6) = .true.
            has_sol_table(7) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '5'
            chan_string(5) = '8'
            chan_string(6) = '10'
            chan_string(7) = '11'
            chan_string(20) ='12'
         
         case ('MTSAT-1R')   sensor_block
            has_sol_table(1) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '1'
            chan_string(20) = '5'
         case ('FY4A') sensor_block
           
            has_sol_table(1) = .true.
            has_sol_table(5) = .false.
            has_sol_table(6) = .true.
            has_sol_table(7) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '2'
            chan_string(5) = '5'
            chan_string(6) = '5'
            chan_string(7) = '6'
            chan_string(20) ='7'

                   
         case ('MTSAT-2')   sensor_block
            has_sol_table(1) = .true.
            has_sol_table(20) = .true.
            has_ems_table(20) = .true.
            chan_string(1) = '1'
            chan_string(20) = '5'
            !--> sensor_identifier = trim(self % lut_path) //'MTSAT'
            sensor_identifier = 'MTSAT'

         case ('FY3D') sensor_block
            has_sol_table(1) = .true. !0.64um
            has_sol_table(6) = .true. !1.6um
            has_sol_table(7) = .true. !2.2um
            has_sol_table(20) = .true. !3.9um
            has_ems_table(20) = .true.
            chan_string(1) = '3'
            chan_string(6) = '6'
            chan_string(7) = '7'
            chan_string(20) = '20'
            
          case ('METIMAGE') sensor_block
            has_sol_table(1) = .true. !0.64um
            has_sol_table(6) = .true. !1.6um
            has_sol_table(7) = .true. !2.2um
            has_sol_table(20) = .true. !3.9um
            has_ems_table(20) = .true.
            chan_string(1) = '03'
            chan_string(6) = '10'
            chan_string(7) = '11'
            chan_string(20) = '12'   
            sensor_identifier = 'eps'
            
          case ('Met_Image') sensor_block
            has_sol_table(1) = .true. !0.64um
            has_sol_table(6) = .true. !1.6um
            has_sol_table(7) = .true. !2.2um
            has_sol_table(20) = .true. !3.9um
            has_ems_table(20) = .true.
            chan_string(1) = '03'
            chan_string(6) = '10'
            chan_string(7) = '11'
            chan_string(20) = '12'   
            sensor_identifier = 'eps'

         case default
            print*,' Sensor is ', trim(self%sensor)
            stop 'add sensor in dcomp_lut_mod.f90 routine populate...'
  
      end select sensor_block

      sensor_identifier = trim(self % lut_path) &
         & // trim(sensor_identifier) 
      
      
      loop_channel : do i_chn = 1 , n_channels
         if ( .not. has_sol_table ( i_chn ) )  cycle
         loop_phase: do i_phase = 1 , 2
  
           
            self%channel(i_chn)%phase(i_phase)%file = trim(sensor_identifier) &
                       
               & // '_ch'//trim ( chan_string ( i_chn ) ) &
               & //'_ref_lut_'//phase_string(i_phase)//'_cld.hdf'
                   
                  
            if ( has_ems_table(i_chn) ) then
               
               self%channel(i_chn)%phase(i_phase)%file_ems = trim(sensor_identifier) &
                     
                  & // '_ch'//trim ( chan_string ( i_chn ) ) &
                  & //'_ems_lut_'//phase_string(i_phase)//'_cld.hdf'
               self % channel(i_chn) % phase(i_phase) % has_ems = .true.
            end if
            
         
            self % channel(i_chn) % phase(i_phase) % has_sol = .true.
                    
         end do loop_phase
      end do loop_channel
      
   
   
   end subroutine lut__set_filename
   
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut__initialize ( self, sensor , ancil_path )
      class ( lut_type ) :: self
      character ( len = * ) , intent(in) :: sensor
      character ( len = * ) , intent(in), optional :: ancil_path
      character ( len =20) :: host
      integer :: i_chn
      
      ! - check if sensor is already initialized
      if ( self % sensor == sensor ) then
         return
      end if
      
      call getenv("HOST",host)
     

      ! - some lut paths
      self % lut_path = '/Users/awalther/DATA/Ancil_Data/clavrx_ancil_data/static/luts/cld/'
      if ( host(1:4) == 'saga' ) self % lut_path = '/data/Ancil_Data/clavrx_ancil_data/static/luts/cld/'
      if ( present(ancil_path)) self % lut_path = trim(ancil_path)
      self % sensor = trim(sensor)
      
      if (.not. allocated (self % channel ) ) allocate ( self % channel (NUM_CHN))
      
      do i_chn = 1 , NUM_CHN
         if (.not. allocated (self % channel (i_chn) % phase ) ) allocate ( self % channel (i_chn) % phase (NUM_PHASE) )
      
      end do
      
      ! - set filenames
      call self % set_filename
     
      ! - clear memory for new sensor
      call self % clear_lut_memory
      
      call self % init_dims( )
      
   end subroutine lut__initialize
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut__set_angles ( self, sat , sol , azi )
      implicit none
      class ( lut_type ) :: self
      real , intent(in) , optional :: sat
      real , intent(in) , optional :: sol
      real , intent(in) , optional :: azi
      
      if ( present(sat) ) sat_m = sat
      if ( present(sol) ) sol_m = sol
      if ( present(azi) ) azi_m = azi
      
      ! - compute pos and weights

      call dcomp_interpolation_weight(self%dims%n_sat_zen , sat , self%dims%sat_zen &
         &, near_index = self % pos_sat  )

      call dcomp_interpolation_weight(self%dims%n_sol_zen , sol , self%dims%sol_zen &
         &, near_index = self % pos_sol  )

      call dcomp_interpolation_weight(self%dims%n_rel_azi , azi , self%dims%rel_azi &
         &, near_index = self % pos_azi  )

   
   end subroutine lut__set_angles
   
   ! --------------------------------------------------------------------------
   !
   !   clears memory of LUT_data object
   !
   ! --------------------------------------------------------------------------
   subroutine lut__clear_lut_memory ( self )
      class ( lut_type ) , target :: self
      integer :: idx_phase , idx_chn
      type (lut_data_type), pointer :: data_loc => null()
   
   
      do idx_phase =1, NUM_PHASE
         do idx_chn =1, NUM_CHN
            data_loc => self % channel ( idx_chn ) % phase ( idx_phase)
            call data_loc % dealloc()
         
         end do
      end do
   
   
   
   end subroutine lut__clear_lut_memory
   
   !  --------------------------------------------------------------------------
   !   PURPOSE : return data from LUT
   !   input: channel, phase, cod (log10 ), cps(log10)
   !   output : transmission, reflectance, cloud albedo, spherical abedo
   !          several derivates
   !   derivates we need:
   !
   !       drefl_dcps
   !       drefl_dcod
   !
   !
   !  -------------------------------------------------------------------------
   subroutine lut__get_data ( self, idx_chn , idx_phase , cod_log10, cps_log10 &
      & , out )
                              
      class ( lut_type ) , target :: self
      integer , intent(in) :: idx_chn
      integer , intent(in) :: idx_phase
      real, intent(in) :: cod_log10
      real, intent(in) :: cps_log10
      
      type ( lut_output ) :: out
      
      integer :: pos_cod,pos_cps
      real :: wgt_cod, wgt_cps
      real , save :: cod_log10_saved
      
      
      real :: rfl_cld_2x2 (2,2)
      real :: trn_sol_cld_2x2 (2,2)
      real :: trn_sat_cld_2x2 (2,2)
      real :: albsph_cld_2x2 (2,2)
      real :: ems_cld_2x2 ( 2,2)
      real :: trn_ems_cld_2x2 ( 2,2)
      real :: alb_cld_2x2 ( 2,2)
      real :: ref_diff , cod_diff
      
      type (lut_data_type), pointer :: data_loc => null()
      
      out % ems = -999.

      data_loc => self % channel ( idx_chn ) % phase ( idx_phase)

      ! test if this channel is available if not read it from hdf file
      if ( .not. data_loc % is_set ) then
         call data_loc % read_hdf
         data_loc % is_set  = .true.
      end if

      call dcomp_interpolation_weight(self%dims%n_cod, cod_log10,self%dims%cod &
         & , weight_out = wgt_cod, index_out= pos_cod)

      call dcomp_interpolation_weight(self%dims%n_cps, cps_log10,self%dims%cps &
         & , weight_out = wgt_cps, index_out= pos_cps)

      rfl_cld_2x2       = data_loc%cld_refl( pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sol,self%pos_sat,self%pos_azi)
      trn_sol_cld_2x2   = data_loc%cld_trn(pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sol)
      trn_sat_cld_2x2   = data_loc%cld_trn(pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sat)
      albsph_cld_2x2    = data_loc%cld_sph_alb(pos_cps:pos_cps+1,pos_cod:pos_cod+1)
      alb_cld_2x2       = data_loc%cld_alb(pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sol)

      ! - parameter for kernel computation
      ref_diff = 0.2
      cod_diff = 0.1

      call interpolate_2d ( rfl_cld_2x2 , wgt_cps , wgt_cod , ref_diff , cod_diff , out % refl    &
         & , out % dRefl_dcps      , out % dRefl_dcod )
      call interpolate_2d ( trn_sol_cld_2x2 , wgt_cps , wgt_cod , ref_diff , cod_diff , out % trn_sol  &
         & , out % dtrans_sol_dcps , out % dtrans_sol_dcod)
      call interpolate_2d ( trn_sat_cld_2x2 , wgt_cps , wgt_cod , ref_diff , cod_diff , out % trn_sat &
         & , out % dTrans_sat_dcod , out % dTrans_sat_dcps )
      call interpolate_2d ( albsph_cld_2x2  , wgt_cps , wgt_cod , ref_diff , cod_diff , out % albsph &
         & , out % dsph_alb_dcod   , out % dSph_alb_dcps)
      call interpolate_2d ( alb_cld_2x2  , wgt_cps , wgt_cod , ref_diff , cod_diff , out % alb &
         & , out % dalb_dcod   , out % dalb_dcps)
         
      if ( data_loc % has_ems ) then
         ems_cld_2x2 = data_loc%cld_ems (pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sat)
         call interpolate_2d ( ems_cld_2x2, wgt_cps , wgt_cod , ref_diff , cod_diff &
            , out % ems, out % dEms_dcps , out % dEms_dcod )
         
         trn_ems_cld_2x2 = data_loc%cld_trn_ems (pos_cps:pos_cps+1,pos_cod:pos_cod+1,self%pos_sat)
         call interpolate_2d ( trn_ems_cld_2x2, wgt_cps , wgt_cod , ref_diff , cod_diff &
            , out % trn_ems, out % dtrnEms_dcps, out % dtrnEms_dcod )
         
      end if

      cod_log10_saved = cod_log10
   end subroutine lut__get_data
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut_data__thick_cloud_rfl ( self , idx_chn , idx_phase, rfl, ems)
      class ( lut_type ) , target :: self
      integer , intent(in) :: idx_chn
      integer , intent(in) :: idx_phase
      
      real, intent(out) :: rfl(9)
      real, intent(out) :: ems(9)
      
      type (lut_data_type), pointer :: data_loc => null()
      
      data_loc => self % channel ( idx_chn ) % phase ( idx_phase)
      
      
      rfl = data_loc%cld_refl( :,29,self%pos_sol,self%pos_sat,self%pos_azi)
      
      if ( data_loc % has_ems ) then
         ems = data_loc%cld_ems (:,29,self%pos_sat)
      end if
   
   
   end subroutine lut_data__thick_cloud_rfl
      
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut_data__read_hdf ( self )
      
      class ( lut_data_type ) :: self
         
      if ( self % has_sol .or. self % has_ems ) then
         call self % alloc
      end if
      
      if ( self % has_sol ) then
         if ( .not. file_test ( self % file )) then
            print*, 'file not here :   ', self % file
            stop 'stop because file not available channel '
         end if
         call read_hdf_dcomp_data_rfl ( &
            self % file &              ! - input
            , self % cld_alb &
            , self % cld_trn &
            , self % cld_sph_alb &
            , self % cld_refl )
      end if
        
      if ( self % has_ems ) then
         if ( .not. file_test ( self % file_ems )) then
            print*, 'file ems not available channel ',  self % file_ems
            stop 'file ems not available channel '
         end if
         call read_hdf_dcomp_data_ems ( &
            self%file_ems &       ! - input
            , self % cld_ems &        ! - output
            , self % cld_trn_ems)
      end if
     
   end subroutine lut_data__read_hdf
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut__init_dims( self)
      class ( lut_type ) :: self
      character ( len = 1020 )  :: hdf_file
      
      hdf_file = trim(self % channel ( 1) % phase (1 ) % file)
    
      if ( .not. file_test(hdf_file) ) then
         print*,'lut file not existing!  ---- ==> ', trim(hdf_file)
         stop 'lut file not existing! ==> '
      end if
      

      ! this should be read from file, but this is also possible
      self %  dims% n_sat_zen = 45
      self %  dims% n_sol_zen = 45
      self %  dims% n_rel_azi = 45
      self %  dims% n_cod = 29
      self %  dims% n_cps = 9
      


      call read_hdf_dcomp_dims ( hdf_file &
         , self %  dims% sat_zen &
         , self %  dims% sol_zen &
         , self %  dims% rel_azi &
         , self %  dims% cod &
         , self %  dims% cps )
      

            
   end subroutine lut__init_dims
   

   !    ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut_data__alloc ( self)
      class ( lut_data_type ) :: self
      
      allocate ( self % cld_alb (9,29,45))
      allocate ( self % cld_trn (9,29,45))
      allocate ( self % cld_sph_alb (9,29))
      allocate ( self % cld_refl (9,29,45,45,45))
      allocate ( self % cld_ems (9,29,45))
      allocate ( self % cld_trn_ems (9,29,45))
   end subroutine lut_data__alloc
   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut_data__dealloc ( self)
      class ( lut_data_type ) :: self
      
      if ( allocated (self % cld_alb) ) deallocate ( self % cld_alb )
      if ( allocated (self % cld_trn) ) deallocate ( self % cld_trn )
      if ( allocated (self % cld_sph_alb) ) deallocate ( self % cld_sph_alb )
      if ( allocated (self % cld_refl) ) deallocate ( self % cld_refl )
      if ( allocated (self % cld_ems) ) deallocate ( self % cld_ems )
      if ( allocated (self % cld_trn_ems) ) deallocate ( self % cld_trn_ems )
      
      self % is_set = .false.
      
   end subroutine lut_data__dealloc

   
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine lut__getProperty ( self, sat , sol , azi , sensor )
      class ( lut_type) :: self
      real , intent(out) , optional :: sat
      real , intent(out) , optional :: sol
      real , intent(out) , optional :: azi
      character(10), intent(out), optional :: sensor
      
      if ( present(sensor)) sensor = trim(self % sensor)
      if ( present(sat) ) sat = self % dims % sat_zen ( self % pos_sat )
      if ( present(sol) ) sol = self % dims % sat_zen ( self % pos_sol )
      if ( present(azi) ) azi = self % dims % sat_zen ( self % pos_azi )
   
   end subroutine lut__getProperty
   

end module  dcomp_lut_mod
