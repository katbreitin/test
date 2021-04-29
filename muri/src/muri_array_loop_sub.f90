! $Id: muri_array_loop_sub.f90 2424 2017-12-18 14:30:36Z awalther $
!
!

subroutine muri_array_loop (input, output )
   
   use muri_interface_mod,only: &
      muri_in_array_type &
      , muri_out_array_type
   
   use muri_retrieval_mod

   use muri_land_retrieval_mod
   
   !use muri_abi_ocean_retrieval_mod

   !use muri_abi_land_retrieval_mod
   
   
   
   use muri_definitions_mod, only: &
      muri_input_type  &
      , muri_output_type


   use muri_atmospheric_corr_mod, only: &
       muri_transmission_default  &
      ,muri_transmission
   
 
   
   implicit none
   
   type ( muri_in_array_type) , intent(in):: input
   type ( muri_out_array_type) :: output
   type ( muri_input_type) :: inp_pixel
   type ( muri_output_type) :: out_pixel
   integer :: i,j,k
   logical :: print_info
   real :: trans_re_default(6)
   real :: trans_re(6)
   
   print_info = .false.
    
   
   call output % reset
   inp_pixel % path = trim('/apollo/cloud/Ancil_Data/clavrx_ancil_data/static/luts/muri/')
   
   
   print*, 'input month is ', input % month
    
   do i = 1, input % dim(1)
  
      do j = 1, input % dim(2)
       
         if ( .not. input % do_it(i,j)) cycle
         
	 
         inp_pixel % sat = input % sat(i,j)
         inp_pixel % sol = input % sol(i,j)
         inp_pixel % azi = input % azi(i,j)
         inp_pixel % scat_angle = input % scat_angle(i,j)
         inp_pixel % ws = input % windspeed(i,j) ! real wind speed
         inp_pixel % surf_elev = input % surf_elev(i,j) ! surface elevation
	 
	 inp_pixel % ozone = input % ozone(i,j)
         inp_pixel % h2o_conc = input %  h2o_conc(i,j)	
	 	
	 !print*,input % land_class(i,j)
			
		  
			
         inp_pixel % path = input % path
         
         ! atmospheric correction
	 
	 
	
         
          trans_re_default = muri_transmission_default (inp_pixel % sol,inp_pixel % sat)
	 
	  trans_re = muri_transmission(inp_pixel % sol,inp_pixel % sat,inp_pixel % ozone,inp_pixel % h2o_conc)
	  
	  output % trans_re_default ( i,j) = trans_re_default(3) ! just to know band 3
	  output % trans_re ( i,j) = trans_re(3) ! just to know band 3
	 
	 
         do k=1,6
            inp_pixel % rfl(k) = input % ref(k,i,j) * trans_re(k)
         end do
			
            out_pixel% sediment_class=-9 ! fill in value
				out_pixel% aerosol_QA=-9 ! fill in value

            select case (input % land_class(i,j))
			
			   case(0,6,7)
        
            inp_pixel % is_sedimental = sedimental (inp_pixel)
				
				
	   if(inp_pixel % is_sedimental)then
	   out_pixel% sediment_class=1
	   else
	   out_pixel% sediment_class=0
	   end if
				
	    if(inp_pixel % sat.gt.78.or.inp_pixel % sol.gt.78)then
	    out_pixel% aerosol_QA=0
	    else
	    out_pixel% aerosol_QA=1
	    end if
				
			   call muri_algorithm( inp_pixel, out_pixel )
				

				
            case (1)
			   !out_pixel% sediment_class=2 ! over land is not sedimental 
				
           if(inp_pixel % rfl(6).lt.0.001.or.inp_pixel % rfl(6).gt.0.25.or.inp_pixel % sat.gt.78.or.inp_pixel % sol.gt.78)then
				out_pixel% aerosol_QA=0
				!print*,'bad land'
				
            call muri_land_algorithm( inp_pixel, out_pixel )
				
			   else 
				!print*,'good land'
				out_pixel% aerosol_QA=1
            call muri_land_algorithm( inp_pixel, out_pixel )
				
			   end if

            case default
         
            end select
			
      
       !  do k=1, 6
           ! print*,out_pixel % aot
            !output % aot_channel(k,i,j) = out_pixel % aot_channel(k)
       !  end do 
       
         output % aot(i,j) = out_pixel % aot
         output % cm_mode(i,j) = out_pixel % cm_mode
         output % fm_mode(i,j) = out_pixel % fm_mode
         output % fmf ( i,j) = out_pixel % fmf
	 output % sediment_class ( i,j) = out_pixel% sediment_class
	 output % aerosol_QA ( i,j) = out_pixel% aerosol_QA
	 
			

         if ( print_info) then
            print*
            print*,i,j
            print*,'land class:', input % land_class(i,j)
            print*,'refl: ',inp_pixel % rfl
            print*,'cm fm: ',out_pixel % cm_mode, out_pixel % fm_mode, out_pixel %  fmf
            print*,'aod: ',out_pixel % aot
            print*  
         end if
        
      end do
   end do   
   

end subroutine muri_array_loop
