! $Id:$
!
!            mod    ahi  abi
! 0.47       1      3     2
!!0.51       2      4       ! don't use it 
! 0.64       3      1     1 
! 0.86       4      2     2 
! 1.38       26           4
! 1.6        5      6     6
! 2.2        6      7     7
!
!
module cx_muri_clavrx_bridge_mod

   use muri_interface_mod , only: &
         muri_in_array_type &
       , muri_out_array_type
       
   use pixel_common_mod, only: &
      geo &
      ,ch &
      ,cldmask &
      , Sensor &
      ,sfc &
      ,NWP_PIX &
      , ancil_data_dir &
      , image
      
   use constants_mod , only: &
      int1  
   
   type cx_muri_structure
      logical :: is_set = .false.
      real,allocatable :: aod(:,:)
      integer(kind=int1), allocatable :: cm_mode(:,:)
      integer(kind=int1), allocatable :: fm_mode(:,:)
      real, allocatable :: fmr (:,:)
      integer(kind=int1), allocatable :: sediment_class ( :, :)
      integer(kind=int1), allocatable :: aerosol_QA ( :, :)
      integer :: dim(2)
      
      real,allocatable :: trans_re_default(:,:)
      real,allocatable :: trans_re(:,:)
      
      
      contains
      procedure :: allocate =>  cx_muri_structure__allocate
      procedure :: deallocate =>  cx_muri_structure__deallocate
   end type cx_muri_structure
   
   type(cx_muri_structure), target, save :: muri
   logical, public :: muri_flag = .false.
         interface
         subroutine muri_array_loop (a , b )
			
            import muri_in_array_type
            import muri_out_array_type
            type (muri_in_array_type) , intent(in) :: a
            type (muri_out_array_type), intent(out) :: b
            
         end subroutine

      end interface
   
   
contains

   subroutine cx_muri_algorithm(dim1,dim2) 
      implicit none
      integer, intent(in) :: dim1,dim2
      type(muri_in_array_type),save :: input
      type(muri_out_array_type),save :: output
      integer :: i 
      logical :: first_call = .true.
      ! - AHI channels 1-6 mean in CLAVR-x/MODIS convention
       integer, parameter :: ahi_map_modis(6) = [3,4,1,2,6,7] 
      
       integer, parameter :: abi_map_modis(6) = [3,1,2,26,6,7] 
     
      ! it is more efficient to allocate only for first segment
      if ( count ( (geo % solzen .lt. 60.) .and. (geo % solzen .gt. 0.)  ) .LT. 100 ) return
      if (first_call ) print*,'MURI Start'
      first_call = .false.
      ! - this checks if allocation is needed 
      call input % allocate ( dim1,dim2)
 
      input % sol = geo % solzen(1:dim1,1:dim2)
      input % sat = geo % satzen(1:dim1,1:dim2)
      input % azi = geo % relaz(1:dim1,1:dim2)
      input % ozone = NWP_PIX%Ozone(1:dim1,1:dim2)
      input % h2o_conc = NWP_PIX%tpw(1:dim1,1:dim2)
      input % windspeed = NWP_pix%Wnd_Spd_10m(1:dim1,1:dim2)
      input % month = image % time_start % month
      !input % sensor = Sensor%Platform_Name 

      input % land_class = Sfc % Land(1:dim1,1:dim2)
      input % surf_elev = sfc % zsfc(1:dim1,1:dim2)

      input % path = trim(ancil_data_dir)//'static/luts/muri/'
     
      do i=1,6 
          if ( .not. allocated(ch(ahi_map_modis(i))%ref_toa)) then
            if (first_call ) then
              print*,'Not all MURI channels are allocated..'
              print*,'no MURI'
            end if
            return
          end if
          input % ref(i,:,:) = ch(ahi_map_modis(i))%ref_toa(1:dim1,1:dim2)/100.
        
      end do
       
      input % do_it = CLDMASK%cld_mask(1:dim1,1:dim2) == 0 
      
      ! - this checks if allocation is needed 
      call output % allocate(dim1,dim2)
      
      call  muri_array_loop (input, output )
      
     
      
      !call muri % allocate (dim1,dim2)
 
      
      muri % aod(1:dim1,1:dim2) = output % aot(1:dim1,1:dim2)
      muri % fm_mode(1:dim1,1:dim2) = output % fm_mode ( 1:dim1,1:dim2)      
      muri % cm_mode(1:dim1,1:dim2) = output % cm_mode ( 1:dim1,1:dim2)
      muri % fmr(1:dim1,1:dim2)  = output % fmf ( 1:dim1, 1:dim2)
      muri % sediment_class(1:dim1,1:dim2)  = output % sediment_class ( 1:dim1, 1:dim2)
      muri % aerosol_QA(1:dim1,1:dim2)  = output % aerosol_QA ( 1:dim1, 1:dim2)
      
      muri % trans_re_default(1:dim1,1:dim2) =output % trans_re_default( 1:dim1, 1:dim2)
      muri % trans_re(1:dim1,1:dim2) =output % trans_re( 1:dim1, 1:dim2)
     
      call output % deallocate
      call input % deallocate
      ! print*,'end 1'
      ! print*,'maxval cm : ',maxval(muri%cm_mode)
      muri_flag = .true.
     
       return
   end subroutine cx_muri_algorithm
   
   
   
   subroutine cx_muri_structure__allocate(this, dim1, dim2)
      class(cx_muri_structure) :: this
      integer, intent(in) :: dim1,dim2
      
      if ( this % is_set ) then 
         if (dim1 .eq. this % dim(1) .and. dim2 .eq. this % dim(2) ) then
            this % aod = -999.
            this % cm_mode = 0
            this % fm_mode = 0
            this % fmr = -999.
            this % sediment_class = -9
	    this % aerosol_QA = -9
	    this % trans_re_default=-9
	    this % trans_re=-9
	    
            return 
         end if   
      end if
      
      call this % deallocate
      
      allocate ( this % aod ( dim1,dim2))
      allocate ( this % cm_mode (dim1,dim2))
      allocate ( this % fm_mode ( dim1, dim2))
      allocate ( this % fmr ( dim1,dim2))
      allocate ( this % sediment_class( dim1,dim2)) 
      allocate ( this % aerosol_QA( dim1,dim2))
      allocate ( this % trans_re_default( dim1,dim2)) 
      allocate ( this % trans_re( dim1,dim2))   
      
      
      
      this % aod = -999.
      this % cm_mode = 0
      this % fm_mode = 0
      this % fmr = -999.
      this % is_set = .true.
      this % sediment_class = -9
      this % aerosol_QA = -9
      this % trans_re_default=-9
      this % trans_re=-9
      
   end subroutine cx_muri_structure__allocate
   
   subroutine cx_muri_structure__deallocate(this)
      class(cx_muri_structure) :: this
      if ( allocated (this % aod) ) deallocate ( this % aod)
      this % is_set = .false.
      if ( allocated (this % cm_mode) ) deallocate ( this % cm_mode)
      if ( allocated (this % fm_mode) ) deallocate ( this % fm_mode)
      if ( allocated ( this % fmr)) deallocate  ( this % fmr)
      if ( allocated ( this % sediment_class)) deallocate  ( this % sediment_class)
      if ( allocated ( this % aerosol_QA)) deallocate  ( this % aerosol_QA)
      if ( allocated ( this % trans_re_default)) deallocate  ( this % trans_re_default)
      if ( allocated ( this % trans_re)) deallocate  ( this % trans_re)
		
   end subroutine cx_muri_structure__deallocate

end module cx_muri_clavrx_bridge_mod
