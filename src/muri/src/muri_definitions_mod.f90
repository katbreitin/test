! $Id: muri_definitions_mod.f90 2173 2017-04-07 23:39:21Z awalther $

module muri_definitions_mod

    type muri_input_type
      real :: rfl(6) 
      real :: sol
      real :: sat
      real :: azi  
      real :: ws
      real :: ozone
      real :: h2o_conc
      real :: scat_angle  ! scattering angle
      real :: surf_elev   ! surface elevation of land
      real :: lat         !latitude
      real :: lon         !longitude
      integer :: month    ! need for seasonal 
      
      logical :: is_sedimental
      character(len=400) :: path
      integer :: sfc_class
      integer :: land_class
      
      contains
      
      procedure :: info => muri_input_info
      
   end type muri_input_type
   
   type muri_output_type
      real :: aot
      real :: angstrom_exponent
      integer :: cm_mode
      integer :: fm_mode
      real :: fmf
      integer ::sediment_class
      integer ::aerosol_QA
      real :: trans_re_default
      real :: trans_re
      real :: err_n
		
   
   end type  muri_output_type
	
	
!	type surface_refl_type
!      real :: yint644
!      real :: slope644
!      real :: yint466
!      real :: slope466
   
!   end type  surface_type
   
   
   
contains

   subroutine muri_input_info (this) 
      class(muri_input_type) :: this
      integer :: i
      print*,'=====   MURI INPUT TYPE'
      print*,'REFELCTANCES:'
      do i=1,6 
         print*,'channel ',i,this% rfl(i)
      end do
      print*,'solar zenith: ',this%sol
      print*,'satellite zenith: ',this%sat
      print*,'relative azimuth difference: ',this%azi
      print*
   
   end subroutine  muri_input_info 
   

end module muri_definitions_mod
