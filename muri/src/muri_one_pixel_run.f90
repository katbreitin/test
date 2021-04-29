! $Id: muri_one_pixel_run.f90 2184 2017-04-12 20:47:02Z awalther $
program one_pixel_run

   use muri_retrieval_mod
	
	use muri_land_retrieval_mod
   
   use muri_definitions_mod, only: &
      muri_input_type  &
      , muri_output_type
      
   use muri_atmospheric_corr_mod, only: &
       muri_transmission_default 
		   
!   use muri_land_altitude_corr_mod, only: &
!       INT_ELEV 
   ! ---------------------------------------   
      
   implicit none   
   type ( muri_input_type) :: inp
   type ( muri_output_type) :: out
   character (len = 400 ) :: path
   !logical :: sedimental
   
   integer :: i,ii,k
   real :: trans_re(6)
	
   inp % path = trim('/home/mino/')
    !inp % path = trim('/DATA/Ancil_Data/clavrx_ancil_data/static/luts/muri/')
   inp % path = trim('/apollo/cloud/Ancil_Data/clavrx_ancil_data/static/luts/muri/')
   do ii=30,30 
  do i = 30,30
   
   
   
   inp % rfl(1) = 0.1331 
   inp % rfl(2) = 0.2981!0.0981 
   inp % rfl(3) = 0.0435!0.0435 
   inp % rfl(4) = 0.0195!0.0195 
   inp % rfl(5) = 0.0039 
   inp % rfl(6) = 0.0013 
   inp % is_sedimental = .true.
   !inp % is_sedimental = .false.
   inp % sol = 25.51
   inp % sat = 46.063
   inp % azi = 150.9448
	inp % ws = 5.51
	inp % lat= 30
	inp % lon= 100
	
	
	
	

   
   !call inp % info  
   trans_re = muri_transmission_default (inp % sol,inp % sat)
	
	
         !print*,'after trans'
      do k=1,6
            inp % rfl(k) = inp % rfl(k) * trans_re(k)
      end do
		
		
		!call inp % info  
      ! surf_refl = muri_land_surface_refl(inp % sol,inp % sat,inp % rfl)
		 
		 !inp % surface_refl(1)= surf_refl(1)! 
	    !inp % surface_refl(2)= surf_refl(2)! 
	    !inp % surface_refl(3)= surf_refl(3)! 

			
	! This is for over water ( ocean and turbid water)		
   call muri_algorithm (inp,out)
	
	
   
   print*
   if (inp % is_sedimental) then
   print*,'channel 5 reference (860nm) reflectance: ',inp % rfl(5) 
   else
   print*,'channel 4 reference (860nm) reflectance: ',inp % rfl(4)
   end if
   print*,'cm mode, fm mode, fmratio :',out % cm_mode, out % fm_mode, out % fmf
   print*,'AOT reference: ',out % aot
   if (sedimental (inp ) ) then
      print*,'sedimental'
    
   
   end if
   end do
   !do ii=1,6 
   
  ! print*,'AOT channel',ii,out % aot_channel(ii)
   
  ! end do
   end do
   
	
	
	! This is for over land	
	! band 6 reflectance must be 0.005 < inp % rfl(6) < 0.25 and 
	
	inp % rfl(1) = 0.1595
   inp % rfl(2) = 0.1350
   inp % rfl(3) = 0.0991
   inp % rfl(4) = 0.3476
   inp % rfl(5) = 0.1941
   inp % rfl(6) = 0.0884 
   inp % is_sedimental = .true.
   !inp % is_sedimental = .false.
   inp % sol = 28.34
   inp % sat = 43.5827
   inp % azi = 121.1810
	inp % ws = 5.51
	inp % lat= 35.59
	inp % lon= 128.70
   inp % month =5
	
	inp%scat_angle=110.55;
	inp % surf_elev= 0.0 ! unit km 
	
	
	print*,'surface elevation :',inp % surf_elev
	
	
	
	
	
	
	print*,'Overv land aerosol output : '
   
   call muri_land_algorithm (inp,out)
  
   print*,'cm mode always "five", fm mode, fmratio :',out % cm_mode, out % fm_mode, out % fmf
   print*,'AOT reference: ',out % aot
  

end program one_pixel_run
