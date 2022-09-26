! $Id: muri_retrieval_mod.f90 2300 2017-08-29 10:33:06Z awalther $

module muri_retrieval_mod
   use muri_definitions_mod, only: &
      muri_input_type &
      , muri_output_type
   
   use muri_lut_mod, only: &
      lut
   
   use lib_array, only:interp1d
   use aw_lib_array, only: interp4d 
      
   implicit none
   private
   public :: muri_algorithm 
   public :: sedimental
contains

   !   This is the ocean retrieval, we can rename later
   !   The AHI and the ABI are different in parts
   !   I will check if we should combine them here
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
      real :: aod_check_tmp
      real :: err_nleta_temp(4,5)
      real :: fmf_nleta_temp(4,5)
      real :: fmf_nleta_pctg
      real :: aod_nleta_final
      integer :: idx(1), idx2(2)
      real :: val,val2
  
      integer:: k,i,j
      integer :: CHANNEL_REFERENCE 
      
      
      real :: opt_0510(8),opt_0870(8),opt_2300(8) ! NOPT
      real :: opt_b2,opt_b4
      real :: OPTH_ONL(8,2,9)
      

      real :: opt(6)
      real :: wavb(6)    
      ! call inp%info
      

      select case(inp % sensor)
      case('AHI')
         channel_reference = 4
         wavb(2) = 510.
         wavb(4) = 870.
      case('ABI')
         channel_reference = 3               
         wavb(2) = 640.
         wavb(3) = 870.
         wavb(5) = 1600.
         wavb(6) = 2200.    
      end select

      if ( inp% is_sedimental) CHANNEL_REFERENCE = 5


      err_nleta_temp = huge(err_nleta_temp)
      err_sqrt = huge(err_sqrt)
      
      call lut % read_lut ( inp % sol,inp%sat,inp%azi,inp%ws, path = inp % path)
      
      call lut % sub_table (inp % sol,inp%sat,inp%azi,inp%ws)
        
         
          !  allocate with n_opt
        allocate ( refl_toa (   lut%N_OPT, N_FMR, N_BANDS))        
        allocate ( refl_reference (  lut% N_OPT))
        n_opt = lut % n_opt
        
       OPTH_ONL=lut % opt_ocean_x ! will use / modidy later

      if (inp % rfl(CHANNEL_REFERENCE).GE.1.1*lut % refl_fine(CHANNEL_REFERENCE,1,1)) then


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
               end do   
               
       
       !! interpolation here 
       
               aot_temp(i_fmr) = interp1d(refl_reference &
                  , lut % aot_550nm &
                  , inp % rfl(CHANNEL_REFERENCE) &
                  , bounds_error = .false. &
                  , FILL_VALUE = -999.)
  
            
             
                 do i_cha = 1,6
                      do i_fmr =1,11
  
                            refl_corrsp(i_fmr,i_cha) = interp1d(lut % aot_550nm,refl_toa(:,i_fmr,i_cha) &
                                   , aot_temp(i_fmr),bounds_error = .false., FILL_VALUE = -999.) 

                       end do
                 
                  end do
               
            
                  do i_cha = 1,6
               
                        n1 = inp % rfl(i_cha) - refl_corrsp(:,i_cha)
                        n2 = inp % rfl(i_cha) - refl_toa(1,1,i_cha) + 0.01
                        if(n2.LT.0.01) n2=0.01
               
                           err(:,i_cha) = n1/n2
               
                   end do
            
                    if ( inp% is_sedimental) then
           ! version 3 
               err_sqrt = sqrt ( (err(:,4)**2.+ err(:,5)**2.+ err(:,6)**2.)/3.)
             else
      ! small signal 
             err_sqrt =  sqrt ( (err(:,2)**2. + err(:,3)**2. + err(:,4)**2.+ err(:,5)**2.+ err(:,6)**2.)/5)
            end if
            
            val= minval(err_sqrt)
            idx=minloc(err_sqrt)
            aod_nleta_temp(i_fm,i_cm) = aot_temp(idx(1))
          err_nleta_temp(i_fm,i_cm) = val
            fmf_nleta_temp(i_fm,i_cm) = (idx(1) - 1 ) /10.
            
         end do
      end do
      
      
      val2 = minval(err_nleta_temp)
      idx2 = minloc(err_nleta_temp)
     
      
      !- once find a good match between fwd and measurment give output
      
       aod_check_tmp=aod_nleta_temp(idx2(1),idx2(2))
       
       if (aod_check_tmp.LT.0.0) then
       out% aerosol_QA=-1 !aod_check_tmp=0.2
       
       end if
       
       
       if (aod_check_tmp.GE.5.0) aod_check_tmp=5.0
      
      out% aot =aod_check_tmp
      out % fm_mode =  idx2(1)
      out % cm_mode = idx2(2)
      out % fmf = fmf_nleta_temp(idx2(1),idx2(2))
      out % err_n= val2
      
      aod_nleta_final=aod_nleta_temp(idx2(1),idx2(2))
      fmf_nleta_pctg=fmf_nleta_temp(idx2(1),idx2(2))
      
      else ! 0.86um is NOT greather than 1.1 0.86um Rayleigh
      
      out% aot = -1.0
      out % fm_mode =  -1
      out % cm_mode = -1
      out % fmf = -1
      out % err_n = -99
      out% angstrom_exponent=-5
      
      end if
      
      if (idx2(1).gt.0.and.idx2(2).gt.0.and.idx2(1).lt.5.and.idx2(2).lt.6.) then
      
      
        if (fmf_nleta_pctg.ge.0.and.fmf_nleta_pctg.le.1.and.aod_nleta_final.gt.0) then
      
           opt_0510=fmf_nleta_pctg * OPTH_ONL(:,1,idx2(1))+(1-fmf_nleta_pctg) * OPTH_ONL(:,1,idx2(2)+4)
           opt_0870=fmf_nleta_pctg * OPTH_ONL(:,2,idx2(1))+(1-fmf_nleta_pctg) * OPTH_ONL(:,2,idx2(2)+4)
     
            opt_b2 = interp1d(lut % aot_550nm,opt_0510 &
                       , aod_nleta_final,bounds_error = .false., FILL_VALUE = -999. )
      
             opt_b4 = interp1d(lut % aot_550nm,opt_0870 &
                       , aod_nleta_final,bounds_error = .false., FILL_VALUE = -999. )
       
              out% angstrom_exponent= - (log(opt_b2/opt_b4)/log(wavb(2)/wavb(4)))
       
         else
              out% angstrom_exponent= -88
         end if
       
       else 
          out% angstrom_exponent= -87
       end if
       
       
      
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
      
          if (inp % land_class.LT.1.and.inp % rfl(1) .LE. 0.25) then
      
          ! this computes simple linear interpolation
         
         b2_reference=EXP(ALOG(inp%rfl(1))+(ALOG(inp%rfl(6))-ALOG(inp%rfl(1)))*d_wvl_2_1 /d_wvl_6_1) 
!!!!      delta_051_thrsh = -1.2 * inp % rfl(6) + 0.010
	  delta_051_thrsh = -1.2 * inp % rfl(6) + 0.020
	  
          delta_051 = inp % rfl(2) - b2_reference
        
         if ( delta_051.GT.delta_051_thrsh) sedimental = .true.
	 
	 end if
			
			
      end if  
   
   end function sedimental
	
	
	
	
	
   
 
   
end module muri_retrieval_mod
