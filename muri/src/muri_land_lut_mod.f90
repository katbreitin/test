! $Id: muri_land_lut_mod.f90 $
!
!   

module muri_land_lut_mod
   use cx_sds_type_definitions_mod, only:
   
   use cx_sds_io_mod !, only: &
      ! cx_sds_finfo &
      ! , cx_sds_read !, cx_sds_read_6d_real
      
   implicit none 
   
   integer, parameter :: N_BANDS = 3
   integer  ::N_OPT  
   integer, parameter :: N_FINE = 4
   integer, parameter :: N_COARSE = 1
   integer :: n_mode
   integer :: n_sol
   integer :: n_sat
   integer :: n_azi
		
   type muri_land_lut_type
      logical :: is_read = .false.
      integer :: n_opt
      real, allocatable :: sol(:)
      real, allocatable :: sat(:)
      real, allocatable :: azi(:)
      
      real, allocatable  :: aot_550nm (:)
      real, allocatable :: aerosol_land_map(:,:,:)
      
		
		! Apparent reflectance is calculated from formula 
		! 
		!  app_refl = path + (Td Tu Rho/ (1-S Rho))
		!
		! path= atmospheric path reflectance 
		!
		! Td and Tu = downward and upward transmission
		!
		! Rho =  Surface Reflectance 
		!
		! S = atmospheric spherical abledo
		
      
		!real, allocatable :: app_refl(:,:,:,:,:,:) ! apparent reflectance
		real, allocatable :: month ! lon=360 , lat=180,  season=4 
		
      real, allocatable :: path_refl(:,:,:,:,:,:) ! atmospheric path reflectance 
      real, allocatable :: T_up(:,:,:,:,:,:)  ! upward transmission
      real, allocatable :: T_dn(:,:,:,:,:,:)  ! downward transmission
      real, allocatable :: S_bar(:,:,:,:,:,:) ! atmospheric spherical abledo
      
!		real, allocatable :: opt(:,:,:,:,:,:) ! optical depth ( will modify later)
		
		
      real, allocatable :: aot_aer_fine(:,:,:)
      real, allocatable :: aot_aer_coarse(:,:,:)
      real, allocatable :: refl_fine(:,:,:)
      real, allocatable :: refl_coarse(:,:,:)
		
		real, allocatable :: path_refl_x(:,:,:) !  (band=3, opt=8, size=2) atmospheric path reflectance 
		real, allocatable :: Tup_x(:,:,:)  ! (band=3, opt=8, size=2)  upward transmission
		real, allocatable :: Tdn_x(:,:,:)  ! (band=3, opt=8, size=2) downward transmission
		real, allocatable :: Sbar_x(:,:,:) ! (band=3, opt=8, size=2) atmospheric spherical abledo
		
		integer :: land_fine_mode
		
	
      
      
      !real :: aot_aer_fine (8,6,4)
      !real :: aot_aer_coarse (8,6,5)
      !real :: refl_fine (N_BANDS,N_OPT,N_FINE)
      !real :: refl_coarse (N_BANDS,N_OPT,N_COARSE)
      
      contains

      procedure ::read_land_lut => muri_land_lut_type_read_lut 
      procedure ::sub_land_table => muri_land_lut_type_sub_table 
 
 
   end type muri_land_lut_type
   
   type ( muri_land_lut_type),save  :: land_lut ! land_ lut
         
     
   
contains
  !
   !
   !
   subroutine muri_land_lut_type_read_lut(this, sol, sat, azi, path)
      class(muri_land_lut_type ) :: this
      
      real, intent(in) :: sol
      real, intent(in) :: sat
      real, intent(in) :: azi
		!real, intent(in) :: land_fine_mode
     
      character(len = *) , intent(in), optional :: path
      
      
      character(len = 400) :: path_local
      character (len = 400)::land_lut_file
      character (len = 400)::lut_file_1
      integer :: istatus
      integer :: ftype
      integer :: nsds
      integer :: natt
      character ( len = MAXNCNAM), allocatable :: sds_name(:)
      character ( len = MAXNCNAM), allocatable :: att_name(:)
      real,allocatable :: sol_zen_ang(:,:)
      real,allocatable :: sat_zen_ang(:,:)
      logical :: file_exists
      real, allocatable :: temp_2d_real(:,:)
		real, allocatable :: temp_3d_real(:,:,:) ! to read aerosol types map
      real, allocatable :: temp_5d_real(:,:,:,:,:) ! over land LUT dimension is 5D
      integer :: band
      character :: band_string
      integer ::  shp_5d(5)
      integer :: band_arr(3) ! 3 bands 1,3,6
		integer :: fine_mode_aerosol_selection
		
		band_arr(:)=(/1,3,6/)
      
      if ( present(path)) then
        path_local = trim(path)
       else
       path_local = trim('/apollo/cloud/Ancil_Data/clavrx_ancil_data/static/luts/muri/')
      end if 
      
       ! print*,path_local
		  
      if ( this % is_read)  return
      
        land_lut_file =  trim(path_local)//trim('/AHI_Land_Aerosol_LUT_B1_v01.hdf')
        
        !print*,land_lut_file
		  
		 
        istatus = cx_sds_read ( trim(land_lut_file),'Solar_Zenith_Angles', temp_2d_real)
        allocate ( this %sol(size(temp_2d_real(:,1)) ), source = temp_2d_real(:,1))
		   
         
        istatus = cx_sds_read ( trim(land_lut_file),'View_Zenith_Angles',temp_2d_real )
        allocate ( this %sat(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
          
        istatus = cx_sds_read ( trim(land_lut_file),'Relative_Azimuth_Angles', temp_2d_real)
        allocate ( this %azi(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
      
        
        istatus = cx_sds_read ( trim(land_lut_file),'AOT_at_550nm',temp_2d_real)
        ! - add scale factor  Jan 2019 AW
        allocate ( this %aot_550nm(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
        this % aot_550nm(:) = temp_2d_real (:,1) /10.  ! -- new scale factor with v03
		  
		  ! Aerosol map ( longitude , latitude , season)  
		  istatus = cx_sds_read ( trim(land_lut_file),'Aerosol_land_map', temp_3d_real)
		 
		  !print*,shape(temp_3d_real)
		  
        allocate ( this %aerosol_land_map(4,180,360))
	     this %aerosol_land_map(:,:,:)=temp_3d_real
        
        
      
          
        do band = 1,N_BANDS 
          write ( band_string, "(i1)") band_arr(band)
         
          land_lut_file = trim(path_local)//trim('/AHI_Land_Aerosol_LUT_B'//band_string//'_v01.hdf')
  
          INQUIRE(file = land_lut_file,EXIST=file_exists)
          if ( .not. file_exists) then 
            print*,'MURI LUT file not there stopping'
            print*,'CLAVR-x was searching at ',land_lut_file
            stop      
          end if
			 
	
        
      ! Atmospheric Path reflctance 
          istatus = cx_sds_read ( trim(land_lut_file), 'Path_Reflectance_land' , temp_5d_real)
       
        
        if ( band .eq. 1) then
          if ( .not. allocated( this%path_refl)) then
            shp_5d = shape(temp_5d_real)
            n_sol=  shp_5d(1) 
            n_sat = shp_5d(2) 
            n_azi = shp_5d(3) 
            n_opt = shp_5d(4) 
            this % n_opt = n_opt 
            n_mode = shp_5d(5) 
            allocate ( this%path_refl(n_sol,n_sat,n_azi,N_opt,N_bands,N_mode))
          end if  
          this%path_refl = -999.  
        end if
         
       !print*,shape(temp_5d_real)
		 
       this%path_refl(:,:,:,:,band,:) =  0.0001 * temp_5d_real(:,:,:,:,:)
		 
		 
		 ! Transmission up T_up 
          istatus = cx_sds_read ( trim(land_lut_file), 'Total_Scat_Trans_Up_land' , temp_5d_real)
     
        if ( band .eq. 1) then
          
        allocate ( this%T_up(n_sol,n_sat,n_azi,N_opt,N_bands,N_mode)) 
          this%T_up = -999.  
        end if
       
      
       this%T_up(:,:,:,:,band,:) =  0.0001 * temp_5d_real(:,:,:,:,:)
		 
		 
		 ! Transmission down T_dn 
          istatus = cx_sds_read ( trim(land_lut_file), 'Total_Scat_Trans_Dn_land' , temp_5d_real)
     
        if ( band .eq. 1) then
      
            allocate ( this%T_dn(n_sol,n_sat,n_azi,N_opt,N_bands,N_mode))
       
          this%T_dn = -999.  
        end if
       
      
       this%T_dn(:,:,:,:,band,:) =  0.0001 * temp_5d_real(:,:,:,:,:)
		 
		 ! Atmospheric albedo S_bar
          istatus = cx_sds_read ( trim(land_lut_file), 'Atmospheric_Spherical_albedo_land' , temp_5d_real)
     
        if ( band .eq. 1) then
   
            allocate ( this%S_bar(n_sol,n_sat,n_azi,N_opt,N_bands,N_mode))
  
          this%S_bar = -999.  
        end if
       
      
       this%S_bar(:,:,:,:,band,:) =  0.0001 * temp_5d_real(:,:,:,:,:)
		 
      
      
      end do
     
       this % is_read = .true.
      
     
      
   end subroutine muri_land_lut_type_read_lut
	
	
	
   !
   !
   !
   subroutine muri_land_lut_type_sub_table ( this, sol,sat,azi,lat,lon,month)
      use aw_lib_array
      class (  muri_land_lut_type ) :: this
      real, intent(in) :: sol
      real, intent(in) :: sat
      real, intent(in) :: azi
		real, intent(in) :: lat
		real, intent(in) :: lon
		integer, intent(in) :: month
		
      
      integer,parameter :: idp = selected_int_kind(13)
      integer,parameter :: sp = selected_real_kind(p=6,r=37)
      integer,parameter :: dp = selected_real_kind(p=15,r=307)
      real ,allocatable :: temp_4d(:,:,:,:)
      
      integer :: i_band, i_mode, i_opt
		integer :: band
		integer , parameter :: N_BANDs=3
		real :: aerosol_land_map(4,360,180)
	  
      integer :: pos_sol
      integer :: pos_sat
      integer :: pos_azi
		
		integer :: season_index
		integer :: land_fine_mode,shp_2d_x(2)
		integer :: sh1x,sh2x
		real, parameter :: dlat = 0.5, dlon = 0.5
      integer :: ilat, ilon
		
		
		
      
      call dcomp_interpolation_weight(size( this % sol) , sol , this % sol &
         &, near_index =  pos_sol  )

      call dcomp_interpolation_weight(size( this % sat) , sat , this % sat &
         &, near_index =  pos_sat  )

      call dcomp_interpolation_weight(size( this % azi) , azi , this % azi &
         &, near_index = pos_azi  )  
		
		!*********************************************************************************	
			
		! land_fine_mode=fine_mode_aerosol_selection(this%aerosol_land_map,lat,lon,month)
		! 
		! Dec Jan Feb
   	if (month.gt.11.or.month.lt.3) then
   	   season_index=1
   	! Mar Apr May
   	else if (month.gt.2.and.month.lt.6) then
   	   season_index=2
   	! Jun Jul Aug
   	else if (month.gt.5.or.month.lt.9) then
   	   season_index=3
   	! Sep Oct Nov
   	else
		! 
   	   season_index=4
   	end if
		
		season_index=4 ! bcasue month=9 ! temporary ( we need to fix this later)
		!print*,'season_index = ',season_index
		
		! campex is sep , so season_index=4 is correct. need to fix later
		
		
		ilon = 181 + NINT(lon + dlon)
      ilat = 91 - NINT(lat + dlat)
  
      ! Note that fine mode aerosol types are 2,3,4 in MODIS collection 6.1. Continental model is not included
	   ! Continental model is reserved for different retrieval
	
	   ! need to read "aerosol_land_map" as integer in read LUT
       land_fine_mode = NINT(aerosol_land_map(season_index,ilat,ilon) + 2)
		 this %land_fine_mode=land_fine_mode
		
		! print*,'selected land fine mode',land_fine_mode 
		!**********************************************************
		
		 !print*,'shape of atmpspheric path reflectance',shape(this % path_refl)
		! print*,'selected index of sol, sat, azi',pos_sol,pos_sat,pos_azi
		 
		 shp_2d_x=shape(this % path_refl(pos_sol,pos_sat,pos_azi,:,:,land_fine_mode))  
		 sh1x=shp_2d_x(1)
		 sh2x=shp_2d_x(2)
		if (.not. allocated (this%path_refl_x )) allocate ( this%path_refl_x(sh1x,sh2x,2)) 
      ! index  1 is fine, 2 is coarse
      this % path_refl_x(:,:,1)  = this % path_refl(pos_sol,pos_sat,pos_azi,:,:,land_fine_mode) 
		
      this % path_refl_x(:,:,2)  = this % path_refl(pos_sol,pos_sat,pos_azi,:,:,5)
		
		if (.not. allocated (this%Tup_x )) allocate ( this%Tup_x(sh1x,sh2x,2)) 
		
		this % Tup_x(:,:,1)  = this % T_up(pos_sol,pos_sat,pos_azi,:,:,land_fine_mode) 
      this % Tup_x(:,:,2)  = this % T_up(pos_sol,pos_sat,pos_azi,:,:,5)
		
		if (.not. allocated (this%Tdn_x )) allocate ( this%Tdn_x(sh1x,sh2x,2)) 
		
		this % Tdn_x(:,:,1)  = this % T_dn(pos_sol,pos_sat,pos_azi,:,:,land_fine_mode) 
      this % Tdn_x(:,:,2)  = this % T_dn(pos_sol,pos_sat,pos_azi,:,:,5)
		
		if (.not. allocated (this%Sbar_x )) allocate ( this%Sbar_x(sh1x,sh2x,2)) 
		
		this % Sbar_x(:,:,1)  = this % S_bar(pos_sol,pos_sat,pos_azi,:,:,land_fine_mode) 
      this % Sbar_x(:,:,2)  = this % S_bar(pos_sol,pos_sat,pos_azi,:,:,5)
		
		

   end subroutine muri_land_lut_type_sub_table
	
	
	
	
  
    subroutine dcomp_interpolation_weight ( count, value, data_in, weight_out, index_out , near_index)
  
      integer, intent ( in ) :: count
      real , intent ( in ) :: value
      real , dimension(:), intent( in ) :: data_in
	  
      real , intent( out ) , optional :: weight_out
      integer, intent ( out ) , optional :: index_out 
      integer, intent ( out ) , optional :: near_index
      integer :: index2
      real :: weight
      integer :: index
      if ( size(data_in) .le. 1) then

         print*,'input wrong for dcomp_interpolation_weight stopping..'
         stop

      end if
      call locate (data_in, count, value, index)

      index  = max ( 1 , min ( count - 1 , index ) )
      index2 = max ( 1 , min ( count , index + 1 ) )

      weight = (value - data_in( index )) / ( data_in( index2 ) - data_in(index))

      if ( weight < 0. ) then
         index = 1 
         weight = 0.
      end if

      if ( present (near_index ) ) near_index = nint ( index + weight )
      if ( present ( weight_out)) weight_out = weight
      if ( present ( index_out)) index_out = index
  
   end subroutine dcomp_interpolation_weight

   !-------------------------------------------------------------------------
   ! subroutine LOCATE(xx, n, x, j)
   ! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
   !--------------------------------------------------------------------------
   subroutine LOCATE(xx, n, x, j)

      !   Arguments
      integer,                        intent(in)  :: n
      integer,                        intent(out) :: j
      real ,               intent(in)  :: x
      real , dimension(:), intent(in)  :: xx

      !   Local variables
      integer :: i, jl, jm, ju

      ! check input
      if ( size(xx) .le. 1 ) then
         print*,'bad locate input  stopping...'
         print*, 'shape xx: ', shape(xx)
         stop
      end if
      jl = 0
      ju = n + 1
      do i = 1, 2*n
         if (ju-jl <= 1) then
            exit
         end if
         jm = (ju + jl) / 2
         if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
            jl = jm
         else
            ju = jm
         end if
      end do
      if (x == xx(1)) then
         j=1
      else if (x == xx(n)) then
         j = n - 1
      else
         j = jl
      endif

   end subroutine LOCATE



end module muri_land_lut_mod





