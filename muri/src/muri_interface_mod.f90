! $Id: muri_interface_mod.f90 2424 2017-12-18 14:30:36Z awalther $
!
!
module muri_interface_mod

   type muri_in_array_type
      logical :: is_allocated = .false.
      integer :: dim(2)
      integer :: dim12(2)
      real, allocatable :: sol(:,:)
      real, allocatable :: sat(:,:)
      real, allocatable :: azi(:,:)
      real, allocatable :: ref(:,:,:)
!      real, allocatable :: ref_tmp(:,:,:)
      real, allocatable :: windspeed(:,:)
      real, allocatable :: scat_angle(:,:)
      real, allocatable :: surf_elev(:,:)
      real, allocatable :: latitude(:,:)
      real, allocatable :: longitude(:,:)
      real, allocatable :: ozone(:,:)
      real, allocatable :: h2o_conc(:,:)
      real, allocatable :: debra_dc(:,:)
      integer, allocatable :: land_class(:,:)
      integer, allocatable :: muri_cm(:,:)
      logical, allocatable :: do_it(:,:)
      character(len=1020) :: path
      integer :: month
      
      
      
      contains
      procedure :: allocate=>muri_in_array_type__allocate
      procedure :: deallocate=>muri_in_array_type__deallocate
      procedure :: info => muri_in__info
   
   end type muri_in_array_type
   
   type muri_out_array_type
      logical :: is_allocated = .false.
      integer :: dim(2)
      integer :: dim12(2)
      real, allocatable :: aot(:,:)
      real, allocatable :: Angstrom_Exponent(:,:)
      real, allocatable :: fmf(:,:)
      integer, allocatable :: fm_mode(:,:)
      integer, allocatable :: cm_mode(:,:)
      integer, allocatable :: sediment_class(:,:)
      integer, allocatable :: aerosol_QA(:,:)
      
      real, allocatable :: trans_re_default(:,:)
      real, allocatable :: trans_re(:,:)
      real, allocatable :: err_n(:,:)
      
      
      contains
      procedure :: allocate=>muri_out_array_type__allocate
      procedure :: deallocate =>muri_out_array_type__deallocate
      procedure :: reset => muri_out_array_type__reset
   
   end type muri_out_array_type

contains

   subroutine muri_in_array_type__allocate(this,dim1,dim2)
      class(muri_in_array_type) :: this
      integer :: dim1,dim2
     
      
      ! - first check if already allocated with correct dims
      
      if ( this % is_allocated ) then 
         if ( dim1 .eq. this % dim(1) .and.  dim2 .eq. this % dim(2) ) return
      end if
      call this % deallocate()
      this % dim(1) = dim1
      this % dim(2) = dim2
      allocate ( this % sol( dim1,dim2))
      allocate ( this % sat( dim1,dim2))
      allocate ( this % azi( dim1,dim2))
      allocate ( this % windspeed( dim1,dim2))
      allocate ( this % scat_angle(dim1,dim2))
      allocate ( this % surf_elev(dim1,dim2))
      allocate ( this % latitude(dim1,dim2))
      allocate ( this % longitude(dim1,dim2))
      allocate ( this % ozone( dim1,dim2))
      allocate ( this % h2o_conc( dim1,dim2))
      allocate ( this % debra_dc( dim1,dim2))
      allocate ( this % do_it(dim1,dim2))
      allocate ( this % land_class(dim1,dim2))
      allocate ( this % muri_cm(dim1,dim2))

     
      allocate ( this % ref( 6,dim1,dim2))
      this % is_allocated = .true.
      
   
   end subroutine  muri_in_array_type__allocate  
   
   
   subroutine muri_in_array_type__deallocate(this)
      class(muri_in_array_type) :: this
      this % dim = [0,0]
      this % dim12=[0,0]
     
      if (allocated (this % sol) ) deallocate ( this % sol)
      if (allocated (this % sat) ) deallocate ( this % sat)
      if (allocated (this % azi) ) deallocate ( this % azi)
      if (allocated (this % windspeed) ) deallocate ( this % windspeed)
      if (allocated (this % scat_angle) ) deallocate ( this % scat_angle)
      if (allocated (this % surf_elev) ) deallocate ( this % surf_elev)
      if (allocated (this % latitude) ) deallocate ( this % latitude)
      if (allocated (this % longitude) ) deallocate ( this % longitude)
      if (allocated (this % ozone) ) deallocate ( this % ozone)
      if (allocated (this % h2o_conc) ) deallocate ( this % h2o_conc)
      if (allocated (this % debra_dc) ) deallocate ( this % debra_dc)
      if (allocated (this % do_it) ) deallocate ( this % do_it)
      
      if (allocated (this % ref) ) deallocate ( this % ref)
      if (allocated (this % land_class) ) deallocate ( this % land_class)
      if (allocated (this % muri_cm ) ) deallocate (this % muri_cm)
    	
		
       this % is_allocated = .false.
   
   end subroutine  muri_in_array_type__deallocate
   
   subroutine muri_in__info ( this)
      class(muri_in_array_type) :: this
      
      print*,'info has to be installed..'
      
   end  subroutine muri_in__info  
   
   
   !
   !
   !
   subroutine muri_out_array_type__allocate(this,dim1,dim2)
      class(muri_out_array_type) :: this
      integer, intent(in) :: dim1,dim2
      integer :: err_all
     
      if ( this % is_allocated ) then
          if ( dim1 .eq. this % dim(1) .and.  dim2 .eq. this % dim(2) ) return              
           return
      end if 
      
      call this % deallocate()
      
      this % dim(1) = dim1
      
      this % dim(2) = dim2
      
     
      
      allocate ( this % aot( dim1,dim2),stat = err_all)
     
      allocate ( this % Angstrom_Exponent(dim1,dim2))
      allocate ( this % fmf( dim1,dim2))
      allocate ( this % fm_mode( dim1,dim2))
      allocate ( this % cm_mode( dim1,dim2))
      allocate ( this % sediment_class( dim1,dim2))
      allocate ( this % aerosol_QA( dim1,dim2))
      allocate ( this % trans_re_default( dim1,dim2))
      allocate ( this % trans_re( dim1,dim2))
      allocate ( this % err_n ( dim1,dim2))
      
      
      
      this % aot = -999.
      this % Angstrom_Exponent=-999.
      this % fmf = -999.
      this % fm_mode = 0
      this % cm_mode = 0
      this % is_allocated = .true.
      this % sediment_class = -1
      this % aerosol_QA = -1
      this % err_n=-99.

   
   end  subroutine muri_out_array_type__allocate
   
   !
   !
   !
      subroutine muri_out_array_type__deallocate(this)
      class(muri_out_array_type) :: this
      integer :: test
      
      this % dim = [0,0]
      this % dim12= [0,0]
      
      
      
      if (allocated (this% aot) ) deallocate ( this % aot, STAT =  test)
     
     
      if (allocated (this% Angstrom_Exponent) ) deallocate ( this % Angstrom_Exponent)
      
      if (allocated (this% fmf) ) deallocate ( this % fmf)
      
      if (allocated (this% fm_mode) ) deallocate ( this % fm_mode)
      
      if (allocated (this% cm_mode) ) deallocate ( this % cm_mode)
      if (allocated (this% sediment_class) ) deallocate ( this % sediment_class)
      if (allocated (this% aerosol_QA) ) deallocate ( this % aerosol_QA)

      if (allocated (this% trans_re_default) ) deallocate ( this % trans_re_default)
      if (allocated (this% trans_re) ) deallocate ( this % trans_re)
      if (allocated (this% err_n) ) deallocate ( this % err_n)
		
      this % is_allocated = .false.

   end subroutine  muri_out_array_type__deallocate 
   
   
   subroutine muri_out_array_type__reset ( this )
      class(muri_out_array_type) :: this
      this % aot = -999.
      this % Angstrom_Exponent=-999
      this % fmf = -999.
      this % fm_mode = 0
      this % cm_mode = 0
      this % sediment_class = -1
      this % aerosol_QA = -1
      this % trans_re_default = -1
      this % trans_re = -1
      this % err_n = -99.		
   
   
   end subroutine muri_out_array_type__reset

end module muri_interface_mod
