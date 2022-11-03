! $Id: nlcomp_array_loop_sub.f90 182 2018-02-15 13:17:01Z awalther $
!
!  HISTORY: 2014/01/12
!         : AW first verisob of NLCOMP for arrays
! 
!
!
subroutine nlcomp_array_loop_sub ( input , output, debug_mode_user )
   use nlcomp_retrieval_mod
   use dncomp_interface_def_mod, only: &
      dncomp_in_type &
      , dncomp_out_type &
      , EM_Cloud_mask &
      , EM_cloud_Type &
      , EM_snow_class

   use dncomp_precip_mod


   implicit none
   
   integer, parameter :: REAL4 = selected_real_kind(6,37)
   integer, parameter :: INT1 = selected_int_kind(1)   
   integer, parameter :: INT2 = selected_int_kind(2)
   real ( kind = real4), parameter :: MISSING_REAL4 = -999.
   
   type (dncomp_in_type) , intent(in) :: input
   type (dncomp_out_type), intent(out) :: output
   integer , intent(in) , optional :: debug_mode_user
   
   ! parameters
   integer , parameter :: N_CHAN = 44
   real  :: ALBEDO_OCEAN (N_CHAN)
   real, parameter :: SAT_ZEN_MAX = 80.
   real, parameter :: SOL_ZEN_MIN = 90. 
   real, parameter :: PI = 3.14159265359
  
   ! - local logical arrays 
   real,  allocatable    :: air_mass_array( : , : )
   logical , allocatable :: is_cloud( : , : )
   logical , allocatable :: is_obs( : , : )
   logical , allocatable :: is_water_phase( : , : )
   logical , allocatable :: has_city_lights(:,:)
   
   integer ( kind = int2)  , allocatable :: info_flag ( :,:)
   integer ( kind = int2)  , allocatable :: quality_flag ( :,:)
   ! - 
   type ( nlcomp_output_structure ) :: nlcomp_out
   integer :: debug_mode
   integer :: array_dim(2)
   integer :: dim_1 
   integer :: dim_2
   integer :: nr_lines
   integer :: nr_elem
   real :: calib_err ( N_CHAN )
   !real :: cld_height
   real :: cld_press
   real :: cld_temp
   real :: rel_azi , lunar_rel_azi
   real :: sol_zen , sat_zen , lunar_zen
   integer :: line_idx , elem_idx
   integer :: chn_idx
   integer, parameter :: CHN_VIS = 44
   integer, parameter  :: CHN_NIR = 20
   real(kind=real4) :: gas_coeff (3)
   
   real( kind = real4 ) :: trans_ozone ( N_CHAN )
   real( kind = real4 ) :: trans_unc_ozone ( N_CHAN )
   real( kind = real4 ) :: trans_rayleigh ( N_CHAN )
   real( kind = real4 ) :: trans_wvp ( N_CHAN )
   real( kind = real4 ) :: trans_unc_wvp ( N_CHAN )
   real( kind = real4 ) :: trans_total ( N_CHAN )
   real( kind = real4 ) :: assumed_tpw_error 
   real( kind = real4 ), parameter :: ozone_coeff (3)    = [ -0.000606266 , 9.77984e-05,-1.67962e-08 ]
   real ( kind = real4 ) :: refl_toc(N_CHAN)
   real ( kind = real4 ) :: alb_sfc(N_CHAN)
   real ( kind = real4 ) :: alb_unc_sfc(N_CHAN)
   real ( kind = real4 ) :: rad_to_refl_factor
   real(kind=real4)   :: refl_toa = -999.
   
   real :: rad_clear_sky_toa_ch20 = -999.
   real :: rad_clear_sky_toc_ch20 = -999.
   
   real :: obs_vec(2) = [-999.,-999.]
   real :: obs_unc(2) = [-999.,-999.]
   real :: alb_vec(2) = [-999.,-999.]
   real :: alb_unc(2) = [-999.,-999.]
   real :: trans_vec(2) = [-999.,-999.]
   
  
   real , dimension(2) :: state_apriori
   
    ! -- nwp variables 
   real :: ozone_path_nwp
   real,allocatable :: rain_column(:,:)
   ! - executable
   
   debug_mode = 1
   if ( present ( debug_mode_user)) debug_mode = debug_mode_user
   array_dim = shape ( input % sat  )
   dim_1 = array_dim (1) 
   dim_2 = array_dim (2)
   nr_lines = array_dim(2)
   nr_elem = array_dim(1)
      
   ALBEDO_OCEAN (:) = 0.03
   calib_err ( : ) = 0.03
   allocate ( is_obs ( dim_1 , dim_2 ) &
      ,  is_cloud ( dim_1 , dim_2 ) )
   allocate ( is_water_phase (  dim_1 , dim_2 ) )
   allocate ( air_mass_array  ( dim_1 , dim_2 ) )
   allocate ( has_city_lights  ( dim_1 , dim_2 ) )
    
   ! - flag masks
   air_mass_array = 1.0 / cos (input % sat  * PI / 180. ) &
      & + 1.0 / cos ( input % zen_lunar  * PI / 180.)
   
   is_obs = input % is_valid   &
      & .and. input % sat  <= SAT_ZEN_MAX &
      & .and. input % sol  > SOL_ZEN_MIN &
      & .and. input % zen_lunar  < 80 &
      & .and. input % chn(44) % refl   >= 0. &
      & .and. air_mass_array >= 2.
   has_city_lights = input % chn(44) %  rad  > 1.E-06
    
   is_cloud =  is_obs &
      & .and. ( input % cloud_mask == EM_cloud_mask % CLOUDY &
      & .or. input % cloud_mask == EM_cloud_mask % PROB_CLOUDY ) &
      & .and. input % cloud_temp > 10
   is_water_phase = input % cloud_type == EM_cloud_type % FOG &
      &  .or. input % cloud_type == EM_cloud_type % WATER &
      &  .or. input % cloud_type == EM_cloud_type % SUPERCOOLED &
      &  .or. input % cloud_type == EM_cloud_type % MIXED
   !-allocation
 
  
  output = dncomp_out_type ( dim_1, dim_2 )
   
  
   allocate ( info_flag ( dim_1, dim_2))
   allocate ( quality_flag ( dim_1, dim_2))
   
   ! - initialize
   quality_flag  = ibclr ( quality_flag , 0)
   quality_flag  = ibset ( quality_flag , 1)
   quality_flag  = ibset ( quality_flag , 2)
   quality_flag  = ibset ( quality_flag , 3)
   quality_flag  = ibset ( quality_flag , 4)
   quality_flag  = ibset ( quality_flag , 5)
   quality_flag  = ibclr ( quality_flag , 6)
   quality_flag  = ibclr ( quality_flag , 7)
   
   info_flag  = 0 
   

   allocate ( rain_column (nr_elem, nr_lines))
   call compute_rain_column ( input % cloud_temp , rain_column)

   line_loop: do line_idx = 1 , nr_lines
      elem_loop: do elem_idx = 1,   nr_elem
         
         if ( .not. is_cloud (elem_idx,line_idx)  ) cycle elem_loop
         if ( input % chn( CHN_VIS) % refl (elem_idx, line_idx) < 0 ) cycle elem_loop 
         
         ! - set local aliases
         ! cld_height     = input % cloud_hgt (elem_idx,line_idx)
         cld_press      = input % cloud_press (elem_idx,line_idx)
         cld_temp       = input % cloud_temp (elem_idx,line_idx)
         sol_zen        = input % sol  (elem_idx,line_idx)
         lunar_zen      = input % zen_lunar  (elem_idx,line_idx)
         sat_zen        = input % sat  (elem_idx,line_idx)
         rel_azi        = input % azi  (elem_idx,line_idx)
         lunar_rel_azi  = input % azi_lunar  (elem_idx,line_idx)
         ozone_path_nwp = input % ozone_nwp  (elem_idx,line_idx)
         
         
         ! - compute transmission 
              
         loop_chn: do chn_idx = 1 , 44
           
            if ( input % is_channel_on (chn_idx) .eqv. .false.) cycle  loop_chn
            
            trans_block : associate ( tpw_ac => input % tpw_ac (elem_idx,line_idx)  , &
               & press_sfc => input % press_sfc   (elem_idx,line_idx) , &
               & trans_chn20_ac_nadir => input % chn(20) % trans_ac_nadir , &
               & refl_toa => input % chn( chn_idx) % refl (elem_idx, line_idx))
                     
               gas_coeff = input % gas_coeff ( chn_idx) % d  
               trans_ozone( chn_idx ) = 1.
               trans_rayleigh( chn_idx ) = 1.
            
               if ( chn_idx == CHN_VIS ) then
               
                  ! - use default ozone value for bad data
                  if (ozone_path_nwp < 100) ozone_path_nwp = 320
               
                  trans_ozone( chn_idx ) = exp ( -1. * ( ozone_coeff(1) &
                     & + ozone_coeff(2) *  ozone_path_nwp &
                     & + ozone_coeff(3) *  ozone_path_nwp**2))

                  trans_unc_ozone( chn_idx ) = trans_ozone(  chn_idx  ) -  exp ( -1. * ( ozone_coeff(1) &
                     & + ozone_coeff(2) *  (1.1 * ozone_path_nwp) &
                     & + ozone_coeff(3) * (1.1 * ozone_path_nwp)**2))
         
                  trans_ozone( chn_idx ) = min ( trans_ozone(  chn_idx  ) , 1. ) 
         
                  ! - rayleigh
                  trans_rayleigh (  chn_vis  ) = exp (-air_mass_array(elem_idx,line_idx)  &
                     &    * ( 0.044 *  (cld_press / press_sfc )) * 0.84)
            
               end if
              
               assumed_tpw_error = 1.2
                           
               trans_wvp ( chn_idx ) =  exp( - 1. * (gas_coeff(1) &
                  & + gas_coeff(2) * tpw_ac  &
                  & + gas_coeff(3) * ( tpw_ac ** 2 ) ) )
                                
               trans_unc_wvp ( chn_vis ) = abs(trans_wvp( chn_idx ) - exp ( -1. * (gas_coeff(1)   &
                  & + gas_coeff(2) * (assumed_tpw_error * tpw_ac) &
                  & + gas_coeff(3) * ( ( assumed_tpw_error * tpw_ac ) **2 ) ) ) )
                        
            end associate trans_block
           
            trans_total ( chn_idx ) = trans_rayleigh( chn_idx ) * trans_ozone( chn_idx ) * trans_wvp( chn_idx )
            trans_total ( chn_idx ) = trans_ozone( chn_idx ) * trans_wvp( chn_idx )
            
            
            refl_toc( chn_idx ) = refl_toa  * trans_total (chn_idx )
            
            alb_sfc( chn_idx ) =  ( input % chn(chn_idx) % alb_sfc (elem_idx,line_idx) ) / 100.
            
            alb_sfc( chn_idx ) = max ( alb_sfc( chn_idx ) , ALBEDO_OCEAN (chn_idx) )
            
            alb_unc_sfc  (chn_idx) = 0.05
            
           
            
            if ( chn_idx == CHN_NIR ) then
               chn20_block : associate ( rad_toa => input % chn( chn_idx) % rad  (elem_idx, line_idx) &
                  & , sun_earth_dist => input % sun_earth_dist &
                  & , solar_irradiance => input % solar_irradiance (chn_idx) )
               
                 
                  trans_total (chn_idx) = input % chn( chn_idx) % trans_ac_nadir   (elem_idx, line_idx)
                  rad_to_refl_factor = PI / cos ( sol_zen * PI / 180.) / ( solar_irradiance / input % sun_earth_dist ** 2 )
                  refl_toc( chn_idx ) = rad_toa * rad_to_refl_factor
                  
               end associate chn20_block
               
               rad_clear_sky_toc_ch20 = input % chn(20) % rad_clear_sky_toc  (elem_idx, line_idx) 
               rad_clear_sky_toa_ch20 = input % chn(20) % rad_clear_sky_toa  (elem_idx, line_idx)
               
            end if
             
         end do loop_chn
         
          ! - NIR
              
              
             
         obs_vec ( 1 ) = input % chn( CHN_VIS) % refl  (elem_idx, line_idx) / 100.
         obs_unc ( 1 ) =   trans_unc_ozone ( CHN_VIS) +  trans_unc_wvp  ( CHN_VIS)  +calib_err (CHN_VIS)
         alb_vec ( 1 ) =  alb_sfc ( CHN_VIS)
        
         alb_unc ( 1) = 0.05
        
         trans_vec ( 1) = trans_total ( CHN_VIS )
        
         
         
         ! - nir channel
         obs_vec( 2 ) = input % chn( CHN_NIR) % rad (elem_idx, line_idx)
         
         
         obs_unc( 2 ) = max ( trans_unc_wvp  ( CHN_NIR ) , 0.01 )  + calib_err (CHN_NIR)
         
         alb_vec( 2 ) = alb_sfc ( CHN_NIR)
         alb_unc( 2 ) = 0.05
         
         trans_vec ( 2) = trans_total ( CHN_NIR )
             
         ! - apriori
         state_apriori (1) = 0.7 * ( 100. * obs_vec(1) ) ** (0.9)
         state_apriori(1) = log10 ( max ( 0.1 , state_apriori(1) ) ) 
         state_apriori (2) = 1.3
         if  (is_water_phase ( elem_idx, line_idx) ) state_apriori(2) = 1.0
         
         call nlcomp_algorithm ( obs_vec , obs_unc  &
            & , alb_vec , alb_unc  , state_apriori , trans_vec  &
            & , lunar_zen , sat_zen &
            & , lunar_rel_azi , cld_temp , is_water_phase ( elem_idx, line_idx)  &
            & , rad_clear_sky_toc_ch20 , rad_clear_sky_toa_ch20  &
            & , nlcomp_out ,  ancil_path =  input % lut_path,  debug_in =  debug_mode  )
          
         quality_flag (elem_idx,line_idx) = ibset ( quality_flag(elem_idx,line_idx) , 0)
         quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 1)
         quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 2)       
                
         output % cod % d (elem_idx,line_idx) = nlcomp_out % cod       
         output % cps % d (elem_idx,line_idx) = nlcomp_out % cps
         
         output % cod_unc % d ( elem_idx, line_idx) = nlcomp_out % codu
         output % ref_unc % d ( elem_idx, line_idx) = nlcomp_out % cpsu 


         call compute_precip ( output % cod % d ( elem_idx, line_idx) &
            , output % cps % d ( elem_idx, line_idx) &
            , rain_column ( elem_idx, line_idx) &
            , output % rain_probability % d ( elem_idx, line_idx) &
            , output % rain_rate % d ( elem_idx, line_idx) )

     
      end do elem_loop
   end do   line_loop 
   
      ! NLCOMP_INFO_LAND_SEA I0
   where (input % is_land )
      info_flag = ibset (  info_flag , 0 ) 
   end where 
   
    ! - NLCOMP_INFO_SNOW I3
   where ( input % snow_class == EM_snow_class % SNOW )
      info_flag = ibset ( info_flag , 3) 
   end where 
   
   ! - NLCOMP_INFO_SEA_ICE I4
   where ( input % snow_class == EM_snow_class % SEA_ICE )
      info_flag = ibset ( info_flag , 4)
   end where
   
   ! -NLCOMP_INFO_PHASE I5
   where ( .not. is_water_phase)
      info_flag = ibset ( info_flag, 5)
   end where
   
    ! -NLCOMP_INFO_THICK_CLOUD
   where ( output % cod % d > 80 )
      info_flag = ibset ( info_flag, 6)
   end where
   
   ! -NLCOMP_INFO_THIN_CLOUD
   where ( output % cod % d < 4 .and. output % cod % d > 0 )
      info_flag = ibset ( info_flag, 7)
   end where
   
   where ( is_obs .and. .not. is_cloud )
      output % cld_trn_sol % d   =  1.0 
      output % cld_trn_obs % d   =  1.0  
      output % cld_alb % d       =  0.0
      output % cld_sph_alb % d   =  0.0  
      output % cod % d           =  0.0
      output % cps % d           =  -999.0
   end where
   
   output % quality % d = quality_flag
   output % info % d = info_flag
   
   

end subroutine nlcomp_array_loop_sub
