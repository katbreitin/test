! $Id: dcomp_one_pixel_run_program.f90 189 2018-03-09 18:49:39Z awalther $
!
!   HISTORY: This file and program name was dcomp.f90
!           Changed to dcomp_run_pixel
!
!   PURPOSE :: this is a onepixel test program for DCOMP
!
program dcomp_one_pixel_run

   use dcomp_math_tools_mod, only: &
      findinv
   use M_kracken
   
   use dcomp_retrieval_mod, only: &
      dcomp_algorithm, dcomp_output_structure

      use nlcomp_retrieval_mod
   
  
      
   !use dncomp_precip_mod
   
  ! use dcomp_retrieval_test   
     
   implicit none
   real, dimension (20) :: obs, obs_u , alb_sfc , alb_sfc_u , air_trans_ac
   real ::  state_apr(2)
   real :: rad_abv_cld 
   real :: rad_sfc
   character( len = 2) :: color_string
   integer :: start_day
   real :: sol_zen 
   real :: sat_zen 
   real :: rel_azi
   real :: cld_temp
   logical :: snow
   character (len = 100)  :: text
   real :: state_apriori(2)
   
   character(len=20) :: sensor
   integer :: iflen 
   integer :: ier
   
   logical :: water_phase  
  
   type ( dcomp_output_structure ) :: dcomp_results
   
   integer :: dcomp_mode
   integer :: debug_mode
   character ( len = 20) :: host
   character ( len = 1024 ) :: ancil_path
   


   ! - executable  
   call getenv("HOST",host)
   start_day = 30
    
   call kracken ("cmd" , "-obs1 0.5211 -obs2 0.133 -alb1 0.12 -alb2 0.12 &        
                & -ctt 276. -sol 23.33 -sat 21. -azi 80 &
					 & -apr1 1.0 -apr2 1.0 -sen GOES-15 -snow ""#N#"" -wat ""#N#""  &
					 & -radabv 0.0002 -radsfc 0.14 -tvis 0.8 -dbg 0 &
					 & -obsu1 0.01 -obsu2 0.01 -albu1 0.01 -albu2 0.01 -dcm 3 -tnr 0.8  -qq 0.8  ")
    
   obs(1) = rget ("cmd_obs1")
   obs(2) = rget ("cmd_obs2")
   alb_sfc(1) = rget ("cmd_alb1")
   alb_sfc(2) = rget ("cmd_alb2")
   cld_temp = rget ("cmd_ctt")
   sol_zen  = rget ("cmd_sol")
   sat_zen = rget ("cmd_sat")
   rel_azi = rget ("cmd_azi")
   rad_abv_cld = rget("cmd_radabv")
   
   rad_sfc = rget("cmd_radsfc")
   
   alb_sfc_u(1) = rget("cmd_albu1")
   alb_sfc_u(2) = rget("cmd_albu2")
   obs_u(1) = rget("cmd_obsu1")
   obs_u(2) = rget("cmd_obsu2")
   air_trans_ac(1) = rget("cmd_tvis")
   air_trans_ac(2) = rget("cmd_tnr")
 
   
   state_apr(1) = rget ("cmd_apr1")
   state_apr(2) = rget ("cmd_apr2")
   snow = lget("cmd_snow")
  
   water_phase = lget ("cmd_wat")
   dcomp_mode = iget("cmd_dcm")
   debug_mode = iget("cmd_dbg")
   call retrev("cmd_sen", sensor, iflen, ier)  
   !sensor = sensor(:iflen)	

   
   if ( debug_mode == 3 ) then 
      print*, 'DCOMP stand alone processing'
   print*,'to switch on this information set -dbg 3'
   print*,'show usage set -dbg 2'
   print*
   print*,'input:   '
   print*,'obs: ', obs(1:2)
   print*,'sol, sat, rel_azi: ', sol_zen, sat_zen, rel_azi
   print*,'alb_sfc: ',alb_sfc(1:2)
   print*,'cloud top temperature: ', cld_temp,trim('K')
   print*,'snow: ',snow
   print*,'water phase: ', water_phase
      print*,'dcomp mode: ', dcomp_mode
      print*,'above cloud transmission: ', air_trans_ac(1:2)
      print*
   end if
   
   
   if ( debug_mode == 2 ) then 
      print*, 'DCOMP stand-alone processing'
   print*,'to switch on this information set -dbg 2'
   print*
   print*,'usage: ...'
   print*,' this shows the options with the default values if you do not set them: '
   print*
   print*,'./dcomp -obs1 0.5211 -obs2 0.133 -alb1 0.12 -alb2 0.12'    
      print* ,' -ctt 276. -sol 23.33 -sat 21. -azi 80 '
   print*,'-apr1 1.0 -apr2 1.0 -sen GOES-15'
   print*,'-radabv 0.0002 -radsfc 0.14 -tvis 0.8 -dbg 0 '
   print*,'-obsu1 0.01 -obsu2 0.01 -albu1 0.01 -albu2 0.01 -dcm 3 -tnr 0.8  -qq 0.8  '
      print*,' boolean keywords :  -snow -wat'
   print*
   end if
   
   
   
   ancil_path = '/home/wstraka/geocat/data_algorithms/baseline_cloud_micro_day/version_1/'
   ancil_path = '/DATA/Ancil_Data/clavrx_ancil_data/static/luts/cld/' 
   if ( host(1:4) == 'luna' ) ancil_path = '/DATA/Ancil_Data/clavrx_ancil_data/static/luts/cld/' 
   if ( host(1:4) == 'saga' ) ancil_path = '/fjord/jgs/patmosx//Ancil_Data/clavrx_ancil_data_new/static/luts/cld/' 
   if ( host(1:4) == 'odin' ) ancil_path = '/data3/Ancil_Data/clavrx_ancil_data/static/luts/cld/' 
   ! ancil_path ='/data/Ancil_Data/clavrx_ancil_data/static/luts/cld/' 
   
   !  - this is the example forpopulation in advance
!   call lut_obj % populate_all_at_once ( sensor,ancil_path)

  ! obs(1) = 0.644741833
   !obs(2) = 0.144089654
  ! sat_zen = 27.69
  ! sol_zen = 48.166
  ! rel_azi = 127.49
  ! obs_u(1) = 0.03
 !  obs_u(2) = 0.014
 !  air_trans_ac(1) = 0.9999
 !  air_trans_ac(2) = 0.9851
 !  water_phase = .true.
 !  sensor='VIIRS'
 !  alb_sfc(1) = 0.13
 !  alb_sfc(2) = 0.03 
   


   state_apriori (1) = 0.7 * ( 100. * obs(1) ) ** (0.9)
         state_apriori(1) = log10 ( max ( 0.1 , state_apriori(1) ) ) 
         state_apriori (2) = 1.3
         if  (water_phase  ) state_apriori(2) = 1.0

   call dcomp_algorithm ( obs , obs_u , alb_sfc , alb_sfc_u , state_apriori , air_trans_ac &
                              & , sol_zen, sat_zen , rel_azi , cld_temp , water_phase &
  & , rad_abv_cld , rad_sfc , sensor &
  & , dcomp_results , dcomp_mode = dcomp_mode &
  & , debug_in = debug_mode , ancil_path = ancil_path )   ! - output
   write(color_string,'(I2)') 43
   text = '============= DCOMP RESULTS==============='
   print*,achar(27)//'['//color_string//'m '//trim(text)//achar(27)//'[0m'
   print*,' ' , trim(sensor), ' DCOMP output:  COD:' &
       ,dcomp_results%cod ,'REF: ',dcomp_results % cps &
       , 'cloud albedo: ', dcomp_results % codu

 !  call dncomp_precip_algo(270.,22.,8.0,7.2, 1, rain_prb, rain_rate)

 !call nlcomp_algorithm(obs , obs_u , alb_sfc , alb_sfc_u , state_apr , air_trans_ac &
  !                            & , sol_zen, sat_zen , rel_azi , cld_temp , water_phase &
!  & , rad_abv_cld , rad_sfc  &
 ! & , nlcomp_results  &
!  & , debug_in = debug_mode , ancil_path = ancil_path)


!obs(1) = 0.35
!obs(2) = 0.1244
!dcomp_mode = 4


!  call dcomp_algorithm ( obs , obs_u , alb_sfc , alb_sfc_u , state_apr , air_trans_ac &
!                              & , sol_zen, sat_zen , rel_azi , cld_temp , water_phase, snow_class &
!  & , rad_abv_cld , rad_sfc , sensor &
!  & , dcomp_results , dcomp_mode = dcomp_mode &
!  & , debug_in = debug_mode , ancil_path = ancil_path )   ! - output
!   write(color_string,'(I2)') 43
 !  text = '============= DCOMP RESULTS==============='
!   print*,achar(27)//'['//color_string//'m '//trim(text)//achar(27)//'[0m'
!   print*,' ' , trim(sensor), ' DCOMP output:  COD:' &
 !      ,dcomp_results%cod ,'REF: ',dcomp_results % cps &
 !      , 'cloud albedo: ', dcomp_results % cloud_alb_vis



!print*,'nlcomp :', nlcomp_results%cod





  
 ! call dcomp_array()
    
end program dcomp_one_pixel_run
