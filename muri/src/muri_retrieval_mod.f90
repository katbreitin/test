! $Id: muri_retrieval_mod.f90 2300 2017-08-29 10:33:06Z awalther $

module muri_retrieval_mod
   use muri_definitions_mod, only: &
      muri_input_type &
      , muri_output_type
   
   use muri_lut_mod, only: &
      lut
   
   use lib_array, only:interp1d
	
      
   implicit none
   private
   public :: muri_algorithm 
   public :: sedimental
contains

   !
   !
   !
   subroutine muri_algorithm (inp, out)
      implicit none
      type( muri_input_type), intent(in) :: inp
      type( muri_output_type),intent(out) :: out

      integer :: i_fm, i_cm
      integer :: n_opt 
      integer, parameter :: N_FMR=11
      integer, parameter :: N_BANDS = 6
      real :: fmr
      real, allocatable:: refl_toa (:,:,:)
      
      
      integer :: i_cha,i_opt, i_fmr
      real, allocatable :: refl_reference (:)
      real :: aot_temp(N_FMR)
      real :: refl_corrsp(N_FMR,6)
      real :: temp(8)
      real :: aot_allbands(N_FMR,6)
      real :: n1 (N_FMR)
      real :: n2
      real :: err(N_FMR,6)
      real :: err_sqrt (N_FMR)
      real :: aod_nleta_temp(4,5)
      real :: err_nleta_temp(4,5)
      real :: fmf_nleta_temp(4,5)
      integer :: idx(1), idx2(2)
      real :: val,val2
      real :: aod_allbands(4,5,6)
      integer:: k,i,j
      integer :: CHANNEL_REFERENCE 
     
      ! call inp%info
      
      
      
      CHANNEL_REFERENCE = 4
		
		!call sedimental(inp)
		
		
      if ( inp% is_sedimental) CHANNEL_REFERENCE = 5

		
      err_nleta_temp = huge(err_nleta_temp)
      err_sqrt = huge(err_sqrt)
   
      call lut % read_lut ( inp % sol,inp%sat,inp%azi,inp%ws, path = inp % path)
      
      call lut % sub_table (inp % sol,inp%sat,inp%azi,inp%ws)
        
         
          !  allocate with n_opt
        allocate ( refl_toa (   lut%N_OPT, N_FMR, N_BANDS))        
        allocate ( refl_reference (  lut% N_OPT))
        n_opt = lut % n_opt
         
      ! - loop over fine and coarse mode
      do i_fm = 1, 4
         do i_cm = 1, 5
            
            ! loop over fine mode ratio
            do i_fmr = 1,N_FMR
               fmr = (i_fmr-1)/10.
               
               do i_opt = 1, N_OPT
              
                  do i_cha = 1,n_BANDS
                     refl_toa (i_opt,i_fmr,i_cha) = fmr * lut % refl_fine(i_cha,i_opt,i_fm) &
                        & + (1 - fmr) * lut % refl_coarse(i_cha,i_opt,i_cm)
                        
                  end do
               end do  
              
               !- reference channel is #4
               refl_reference = refl_toa(:,i_fmr,CHANNEL_REFERENCE)
               
               
               aot_temp(i_fmr) = interp1d(refl_reference &
                  , lut % aot_550nm &
                  , inp % rfl(CHANNEL_REFERENCE) &
                  , bounds_error = .false. &
                  , FILL_VALUE = -999.)
                  
                  
               ! print*,i_fm,i_cm, i_fmr,aot_temp(i_fmr)
               ! print*,'refl_ch4: ',refl_ch4
               ! print*,inp % rfl(4)
                
                
                  
            end do
            
             do i_cha = 1,6
                  do i_fmr =1,11
                     refl_corrsp(i_fmr,i_cha) = interp1d(lut % aot_550nm,refl_toa(:,i_fmr,i_cha) &
                        , aot_temp(i_fmr),bounds_error = .false., FILL_VALUE = -999. )
                      !  print*,refl_toa(:,i_fmr,i_cha), aot_temp(i_fmr)
                    !  print*,'refl_corrsp: ',i_cha,i_fmr, refl_corrsp(i_fmr,i_cha), inp % rfl(i_cha)
                  end do
                 
               end do
               
            
            do i_cha = 1,6
               
               n1 = inp % rfl(i_cha) - refl_corrsp(:,i_cha)
               n2 = inp % rfl(i_cha) - refl_toa(1,CHANNEL_REFERENCE,i_cha) + 0.001
               
               
               err(:,i_cha) = n1/n2
              ! print*,'=====>',i_cha,inp % rfl(i_cha),refl_corrsp(:,i_cha)
               
               
               
            end do
            
             if ( inp% is_sedimental) then
               err_sqrt = sqrt ( (err(:,5)**2.+ err(:,6)**2.)/2.)
             else
              err_sqrt = sqrt ( (err(:,3)**2. + err(:,4)**2.+ err(:,5)**2.+ err(:,6)**2.)/4.)  
            end if
            
            
           ! print*,'error :', err_sqrt
            
            val= minval(err_sqrt)
            idx=minloc(err_sqrt)
           
            aod_nleta_temp(i_fm,i_cm) = aot_temp(idx(1))
            err_nleta_temp(i_fm,i_cm) = val
            fmf_nleta_temp(i_fm,i_cm) = (idx(1) - 1 ) /10.
           ! print*,i_fm,i_cm,idx(1)
            aod_allbands (i_fm,i_cm,:) = aot_allbands(idx(1),:) 
            
         end do
      end do
      
      
      val2 = minval(err_nleta_temp)
      idx2 = minloc(err_nleta_temp)
     
      
      !- once find a good match between fwd and measurment give output
      out% aot = aod_nleta_temp(idx2(1),idx2(2))
      out % fm_mode =  idx2(1)
      out % cm_mode = idx2(2)
      ! print*,fmf_nleta_temp
      !stop
      out % fmf = fmf_nleta_temp(idx2(1),idx2(2))
      do i_cha =1,6 
         !out% aot_channel(i_cha) = aod_allbands(idx2(1),idx2(2),i_cha)
      end do
		
		
      
   end subroutine  muri_algorithm  
   
   !
   !
   !
   logical function sedimental (inp)
      type( muri_input_type), intent(in) :: inp
      real :: delta_051
      real :: delta_051_thrsh
      real :: b2_reference
      real, parameter :: D_WVL_6_1 = 1.5879  !Alog(2.3) - ALOG(0.47) = 
      real, parameter :: D_WVL_2_1 = 0.0817  ! ALOG(0.51) - ALOG(0.47)
      sedimental = .false.
      
      if ( inp % rfl(6) .GE. 0.2 .OR. ALOG(inp % rfl(1))/ ALOG (inp % rfl(2)) .GE. 0.86 ) then
          ! this computes simple linear interpolation
         
         b2_reference=EXP(ALOG(inp%rfl(1))+(ALOG(inp%rfl(6))-ALOG(inp%rfl(1)))*d_wvl_2_1 /d_wvl_6_1) 
!         delta_051_thrsh = -1.2 * inp % rfl(6) + 0.010
	  delta_051_thrsh = -1.2 * inp % rfl(6) + 0.015
          delta_051 = inp % rfl(2) - b2_reference
        
         if ( delta_051.GT.delta_051_thrsh) sedimental = .true.
			
			
      end if  
   
   end function sedimental
	
	
	
	
	
   
 
   
end module muri_retrieval_mod
