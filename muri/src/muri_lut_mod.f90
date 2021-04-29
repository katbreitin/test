! $Id: muri_lut_mod.f90 2424 2017-12-18 14:30:36Z awalther $
!
!     2019 Feb 6 : fixed scalong of aot_550nm with factor 10

module muri_lut_mod
   use cx_sds_type_definitions_mod, only:
   
   use cx_sds_io_mod !, only: &
      ! cx_sds_finfo &
      ! , cx_sds_read !, cx_sds_read_6d_real
      
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
      real, allocatable :: app_refl(:,:,:,:,:,:)
      !real, allocatable :: aot_aer(:,:,:,:,:,:,:) 
      
      real, allocatable :: aot_aer_fine(:,:,:)
      real, allocatable :: aot_aer_coarse(:,:,:)
      real, allocatable :: refl_fine(:,:,:)
      real, allocatable :: refl_coarse(:,:,:)
      
      
      !real :: aot_aer_fine (8,6,4)
      !real :: aot_aer_coarse (8,6,5)
      !real :: refl_fine (N_BANDS,N_OPT,N_FINE)
      !real :: refl_coarse (N_BANDS,N_OPT,N_COARSE)
      
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
      integer :: i_ws(1)
      
      if ( present(path)) then
        path_local = trim(path)
		  
		 ! print*, path_local
       else
       path_local = trim('/apollo/cloud/Ancil_Data/clavrx_ancil_data/static/luts/muri/')
      end if 
      
    
      if ( this % is_read)  return
      
      lut_file =  trim(path_local)//trim('/AHI_Ocean_Aerosol_LUT_RB1_v04.hdf')
        
       !print*,lut_file
          
     
        istatus = cx_sds_read ( trim(lut_file),'Solar_Zenith_Angles', temp_2d_real)
        allocate ( this %sol(size(temp_2d_real(:,1)) ), source = temp_2d_real(:,1))
		  
		  
		 
         
        istatus = cx_sds_read ( trim(lut_file),'View_Zenith_Angles',temp_2d_real )
        allocate ( this %sat(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
          
        istatus = cx_sds_read ( trim(lut_file),'Relative_Azimuth_Angles', temp_2d_real)
        allocate ( this %azi(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
       
        istatus = cx_sds_read ( trim(lut_file),'Wind_Speed', temp_2d_real)
       
       
        allocate ( this %ws(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
        
        istatus = cx_sds_read ( trim(lut_file),'AOT_at_550nm',temp_2d_real)
        ! - add scale factor  Jan 2019 AW
        allocate ( this %aot_550nm(size(temp_2d_real(:,1))), source = temp_2d_real(:,1))
        this % aot_550nm(:) = temp_2d_real (:,1) /100.  ! -- new scale factor with v03
        
        
      
          
        do band = 1,N_BANDS 
          write ( band_string, "(i1)") band
         
          lut_file = trim(path_local)//trim('/AHI_Ocean_Aerosol_LUT_RB'//band_string//'_v04.hdf')
  
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
            allocate ( this%app_refl(n_sol,n_sat,n_azi,N_bands,N_opt,N_mode))
          end if  
          this%app_refl = -999.  
        end if
       
       
       !-  closest wind speed
       ! - 16 Juy 2019 AW
       
      
        i_ws = minloc ( abs( ws - this % ws))
       
      
       this%app_refl(:,:,:,band,:,:) =  0.0001 * temp_6d_real(:,:,:,:,i_ws(1),:)
      
      
      end do
     
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
      
      
      
      
      call dcomp_interpolation_weight(size( this % sol) , sol , this % sol &
         &, near_index =  pos_sol  )

      call dcomp_interpolation_weight(size( this % sat) , sat , this % sat &
         &, near_index =  pos_sat  )

      call dcomp_interpolation_weight(size( this % azi) , azi , this % azi &
         &, near_index = pos_azi  )    
      
      this % refl_fine  = this % app_refl(pos_sol,pos_sat,pos_azi,:,:,1:4) 
      this % refl_coarse  = this % app_refl(pos_sol,pos_sat,pos_azi,:,:,5:9)

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





