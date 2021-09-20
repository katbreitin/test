! $Id: muri_lut_mod.f90 2424 2017-12-18 14:30:36Z awalther $
!
!     2019 Feb 6 : fixed scalong of aot_550nm with factor 10

module muri_lut_mod
   use cx_sds_type_definitions_mod, only:
   
   use cx_sds_io_mod !, only: &
      ! cx_sds_finfo &
      ! , cx_sds_read !, cx_sds_read_6d_real
      
      use lib_array, only: interp1d
      use aw_lib_array, only: interp4d 
      
   implicit none 
   
   integer, parameter :: N_BANDS = 6
      integer :: N_OPT
      integer, parameter :: N_FINE = 4
      integer, parameter :: N_COARSE = 5
      integer :: n_ws
      integer :: n_mode
      integer :: n_sol
      integer :: n_sat
      integer :: n_azi
   
   type muri_lut_type
      logical :: is_read = .false.
      integer :: n_opt
      real, allocatable :: sol(:)
      real, allocatable :: sat(:)
      real, allocatable :: azi(:)
      real, allocatable :: ws(:)
      
      real, allocatable  :: aot_550nm (:)
      real, allocatable :: app_refl(:,:,:,:,:,:,:)
      
      
      !real, allocatable :: aot_aer_fine(:,:,:)
      !real, allocatable :: aot_aer_coarse(:,:,:)
      !real, allocatable :: refl_fine(:,:,:)
      !real, allocatable :: refl_coarse(:,:,:)
      
  
      
       real :: refl_fine (6,8,4) ! (6,6,4) for v213
       real :: refl_coarse (6,8,5)
       real :: OPT_ocean_Bands(8,9,2) ! optical depth of 9 models in 2 bands (0.51,0.86)
       real :: opt_ocean_x (8,2,9) ! 8 opt x
      
      contains
      procedure ::read_lut => muri_lut_type__read_lut 
      procedure ::sub_table => muri_lut_type__sub_table 
 
   end type muri_lut_type
   
	
   type ( muri_lut_type),save  :: lut
         
     
   
contains
  !
   !
   !
   subroutine muri_lut_type__read_lut(this, sol, sat, azi, ws, path)
      class(muri_lut_type ) :: this
      
      real, intent(in) :: sol
      real, intent(in) :: sat
      real, intent(in) :: azi
      real, intent(in) :: ws 
      character(len = *) , intent(in), optional :: path
      
      
      character(len = 400) :: path_local
      character (len = 400)::lut_file
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
      real, allocatable :: temp_6d_real(:,:,:,:,:,:)
      integer :: band
      character :: band_string
      integer ::  shp_6d(6)
      integer :: i_ws,i_ws_a(1)
      integer :: i_opt,i_mode
      
      if ( present(path)) then
        path_local = trim(path)
		  
       else
       
        path_local = trim('/home/mino/iris-home/MURI/MURI_aerosol_LUT/')
      end if 
      
    
      if ( this % is_read)  return
      
      
       lut_file =  trim(path_local)//trim('/AHI_Ocean_Aerosol_LUT_RB2_v212A.hdf')
       
       ! the following geo parameters are the same for all bands.
       ! (ABI band 1 is not using in over ocean retrieval) 
        print*,lut_file
          
     
        istatus = cx_sds_read ( trim(lut_file),'Solar_Zenith_Angles', temp_2d_real)
        allocate ( this %sol(size(temp_2d_real(:,1)) ), source = temp_2d_real(:,1))
	print*,this%sol
        istatus = cx_sds_read ( trim(lut_file),'View_Zenith_Angles',temp_2d_real )
        allocate ( this %sat(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
         print*,this%sat 
        istatus = cx_sds_read ( trim(lut_file),'Relative_Azimuth_Angles', temp_2d_real)
        allocate ( this %azi(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
         print*,this%azi   
        istatus = cx_sds_read ( trim(lut_file),'Wind_Speed', temp_2d_real)
       
        allocate ( this %ws(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
        
        istatus = cx_sds_read ( trim(lut_file),'AOT_at_550nm',temp_2d_real)
        ! - add scale factor  Jan 2019 AW
        !allocate ( this %aot_550nm(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
	
	allocate ( this %aot_550nm(size(temp_2d_real(1,:))), source = temp_2d_real(1,:))
        this % aot_550nm(:) = temp_2d_real (1,:) /100.  ! -- new scale factor with v03
        print*,'aot',this%aot_550nm
        
      
          
        do band = 1,N_BANDS 
          write ( band_string, "(i1)") band
         
          lut_file = trim(path_local)//trim('/AHI_Ocean_Aerosol_LUT_RB'//band_string//'_v212A.hdf')
  
          INQUIRE(file = lut_file,EXIST=file_exists)
          if ( .not. file_exists) then 
            print*,'MURI LUT file not there stopping'
            print*,'CLAVR-x was searching at ',lut_file
            stop      
          end if
        
     
          istatus = cx_sds_read ( trim(lut_file), 'Apparent_Reflectance_ocean' , temp_6d_real)
      
          
        
        if ( band .eq. 1) then
          if ( .not. allocated( this%app_refl)) then
            shp_6d = shape(temp_6d_real)
            n_sol= shp_6d(1) 
            n_sat = shp_6d(2) 
            n_azi = shp_6d(3) 
            n_opt = shp_6d(4) 
            this % n_opt = n_opt
            n_ws = shp_6d(5)  
            n_mode = shp_6d(6) 
           ! allocate ( this%app_refl(n_sol,n_sat,n_azi,N_bands,N_opt,N_mode))
	    allocate ( this%app_refl(n_sol,n_sat,n_azi,n_ws,N_bands,N_opt,N_mode))
          end if  
          this%app_refl = -999.  
        end if
        !print*,'size this%app_refl',shape(this%app_refl)
       
       !-  closest wind speed
       ! - 16 Juy 2019 AW
       !i_ws = minloc ( abs( ws - this % ws))
       !this%app_refl(:,:,:,band,:,:) =  0.0001 * temp_6d_real(:,:,:,:,i_ws(1),:)
       
       ! in this new version, we will interpolate 
       
       
        do i_opt=1,n_opt
	   do i_mode=1,n_mode
	   do i_ws=1,n_ws
       
        this%app_refl(:,:,:,i_ws,band,i_opt,i_mode) =  0.0001 * temp_6d_real(:,:,:,i_opt,i_ws,i_mode)
	end do
	end do
	end do
        
     
     
           
     ! AOT_at_B2, AOT_at_B4, AOT_at_B6
     
          
	   
	   
           if(band .eq. 2) then
           istatus = cx_sds_read ( trim(lut_file), 'AOT_at_B2' , temp_2d_real)
           this% opt_ocean_x(:,1,:) =  0.001 * temp_2d_real(:,:) 
	   end if 
	    
	   
	   if(band .eq. 4) then
           istatus = cx_sds_read ( trim(lut_file), 'AOT_at_B4' , temp_2d_real)
           this % opt_ocean_x(:,2,:) =  0.001 * temp_2d_real(:,:) 
	   end if  
	   
	   !!if(band .eq. 6) then
           !!istatus = cx_sds_read ( trim(lut_file), 'AOT_at_B6' , temp_2d_real)
           !!this%OPT_ocean_Bands(:,3,:) =  0.001 * temp_2d_real(:,:) 
	   !!end if     
       
           
     
      
       end do
      
       print*,'B2 refl check',this%app_refl(5,5,5,2,2,:,1)
       print*,'B4 refl',this%app_refl(5,5,5,2,4,:,1)
       print*,'rayleigh refl',this%app_refl(5,5,5,2,4,1,1)
       !print*,'this % opt_ocean_x',this % opt_ocean_x
       this % is_read = .true.
      
      
      
   end subroutine muri_lut_type__read_lut
   !
   !
   !
   subroutine muri_lut_type__sub_table ( this, sol,sat,azi,ws )
      use aw_lib_array
      class (  muri_lut_type ) :: this
      real, intent(in) :: sol
      real, intent(in) :: sat
      real, intent(in) :: azi
      real, intent(in) :: ws 
      
      integer,parameter :: idp = selected_int_kind(13)
      integer,parameter :: sp = selected_real_kind(p=6,r=37)
      integer,parameter :: dp = selected_real_kind(p=15,r=307)
      real ,allocatable :: temp_4d(:,:,:,:)
      
      integer :: i_band, i_mode, i_opt
      
      integer :: pos_sol
      integer :: pos_sat
      integer :: pos_azi
      
      
      
            do i_band=1,6
       
           do i_opt=1,8 ! 8
	      do i_mode=1,4

             
	      
              this % refl_fine(i_band,i_opt,i_mode)=interp4d(this%sol,this%sat &
	                                                      ,this%azi, this%ws &
							      ,this%app_refl(:,:,:,:,i_band,i_opt,i_mode) &
							       ,sol,sat,azi,ws &
							       , bounds_error = .false., FILL_VALUE = -999.)  
							       
		
       
	       end do
	   end do
	 end do
	 
          
	  
	  
	    ! do i_mode=1,4
	  
            !if (this % refl_fine(4,1,i_mode).LT.0) then
		
		!print*,'sol,sat,azi,ws, why why',sol,sat,azi,ws
		!print*,'this%ws',this%ws
		!print*,'app refl',this%app_refl(7:8,12:13,12:13,1:2,4,1,i_mode)
		!print*,'this % refl_fine(4,1,i_mode)',this % refl_fine(4,1,i_mode)
	
		!end if
		
	    !end do					       
			
      
      
      
      
      
	 
       do i_band=1,6
       
           do i_opt=1,8 ! 8
	      do i_mode=1,5
				 
 	      this % refl_coarse(i_band,i_opt,i_mode)=interp4d(this%sol,this%sat &
	                                                      ,this%azi, this%ws &
							      ,this%app_refl(:,:,:,:,i_band,i_opt,i_mode+4) &
							       ,sol,sat,azi,ws &
							       , bounds_error = .false., FILL_VALUE = -999.)  	          
	           	          
	           
       
	       end do
	   end do
	 end do  
	 
	 
	! do i_band=1,2
	! 	do i_opt=1,8
	!	      do i_mode=1,9
	!	
	! this % opt_ocean_x(i_opt,i_band,i_mode)= this%OPT_ocean_Bands(i_opt,i_mode,i_band)
	!		end do				       
	! 	end do
	!end do	
	 
	     
      
!      call dcomp_interpolation_weight(size( this % sol) , sol , this % sol &
!         &, near_index =  pos_sol  )

!     call dcomp_interpolation_weight(size( this % sat) , sat , this % sat &
!         &, near_index =  pos_sat  )

!      call dcomp_interpolation_weight(size( this % azi) , azi , this % azi &
!         &, near_index = pos_azi  )    
      
!      this % refl_fine  = this % app_refl(pos_sol,pos_sat,pos_azi,:,:,1:4) 
!      this % refl_coarse  = this % app_refl(pos_sol,pos_sat,pos_azi,:,:,5:9)

   end subroutine muri_lut_type__sub_table
   
  
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



end module muri_lut_mod





