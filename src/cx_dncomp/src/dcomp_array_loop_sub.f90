! $Id: dcomp_array_loop_sub.f90 189 2018-03-09 18:49:39Z awalther $
!
!  HISTORY: 06/05/2014: changed filename for better naming convebtion
!
!           02/05/2014 : add AHI (AW)
subroutine dcomp_array_loop ( input, output , debug_mode_user)

   use univ_kind_defs_mod, only: f4,  i2

   use dcomp_retrieval_mod, only: &
      dcomp_output_structure &
      , dcomp_algorithm

   use dncomp_interface_def_mod, only: &
      dncomp_in_type &
      , dncomp_out_type &
      , EM_cloud_type &
      , EM_cloud_mask &
      , EM_snow_class

   use dncomp_trans_atmos_mod,only: &
      trans_atm_above_cloud

   use dncomp_precip_mod

   implicit none

   type (dncomp_in_type) , intent(in) :: input
   type (dncomp_out_type), intent(inout) :: output
   integer , intent(in) , optional :: debug_mode_user

   integer :: nr_lines, nr_elem


    !- scalar local variables

   ! - number of possible channels in CLAVR-x
   integer , parameter :: N_CHAN = 45

   real :: refl_toa = -999.

   real ::calib_err ( N_CHAN )

   real :: obs_vec(2)   = [-999.,-999.]
   real :: obs_unc(2)   = [-999.,-999.]
   real :: alb_vec(3)   = [-999.,-999.,-999.] ! third is ch1 dark sky
   real :: alb_unc(2)   = [-999.,-999.]
   real :: trans_vec(2) = [-999.,-999.]

   real :: trans_unc_total ( N_CHAN )
   real :: trans_total ( N_CHAN )

   integer, parameter  :: CHN_VIS_DEFAULT = 1
   integer  :: CHN_NIR_DEFAULT
   integer  :: CHN_VIS
   integer  :: CHN_NIR
   integer  :: chn_idx

   real( f4 ), parameter :: OZONE_COEFF_CHN1 (3)  = [ -0.000606266 , 9.77984e-05,-1.67962e-08 ]
   real( f4 ) :: ozone_coeff (3)
   real( f4 ) :: gas_coeff (3)

   ! -- nwp variables
   real( f4 ) :: ozone_dobson

   real :: rad_clear_sky_toa_ch20 = -999.
   real :: rad_clear_sky_toc_ch20 = -999.

   real , allocatable :: air_mass_array (:,:)

   logical  , allocatable :: cloud_array(:,:)
   logical  , allocatable :: obs_array(:,:)
   logical  , allocatable :: obs_and_acha_array(:,:)
   logical  , allocatable :: water_phase_array(:,:)

   integer ( i2 )  , allocatable :: info_flag ( :,:)
   integer ( i2 )  , allocatable :: quality_flag ( :,:)
   integer :: dim_1
   integer :: dim_2
   integer :: dim_1_w
   integer :: dim_2_w

   integer  :: array_dim (2)
   real     :: state_apriori (2)


   integer :: debug_mode

   character (len = 20) :: sensorname_from_wmoid

   real :: cld_press
   real :: cld_temp

   real :: rel_azi

   real ( f4 ) :: refl_toc(N_CHAN)
   real ( f4 ) :: alb_sfc(N_CHAN)
   real ( f4 ) :: alb_unc_sfc(N_CHAN)
   real ( f4 ) :: rad_to_refl_factor

   real, parameter :: SAT_ZEN_MAX = 75.
   real, parameter :: SOL_ZEN_MAX = 82.
   real, parameter :: PI = 4. * ATAN(1.)

   real ( f4) :: ALBEDO_OCEAN (N_CHAN)

   type ( dcomp_output_structure ) :: dcomp_out

   real :: sol_zen
   real :: sat_zen

   integer :: line_idx
   integer :: elem_idx
   integer :: tried
   integer :: success

   integer :: dcomp_mode

   !real,allocatable :: rain_column(:,:)

   ! - executable ---------

   debug_mode = 5
   if ( present ( debug_mode_user)) debug_mode = debug_mode_user

   array_dim = shape ( input % sat  )
   dim_1 = array_dim (1)
   dim_2 = array_dim (2)
   nr_lines = array_dim(2)
   nr_elem = array_dim(1)

   ALBEDO_OCEAN (:) = 0.03
   calib_err ( : ) = 0.03


   allocate ( obs_array ( dim_1 , dim_2 ) &
      ,  cloud_array ( dim_1 , dim_2 ) , obs_and_acha_array ( dim_1 , dim_2 ) )
   allocate ( water_phase_array (  dim_1 , dim_2 ) )
   allocate ( air_mass_array  ( dim_1 , dim_2 ) )

   air_mass_array = 1.0 / cos (input % sat  * PI / 180. ) &
        & + 1.0 / cos ( input % sol  * PI / 180.)

   obs_array = input % is_valid    &
      & .and. input % sat  <= SAT_ZEN_MAX &
      & .and. input % sol  <= SOL_ZEN_MAX &
      & .and. input % chn(1) % refl    >= 0. &
      & .and. air_mass_array >= 2.

	! - check input options

   select case ( input % mode )
      case ( 1 )
         CHN_NIR_DEFAULT = 6
         obs_array = obs_array .and. input % chn(6) % refl >=0.
      case ( 2 )
         CHN_NIR_DEFAULT = 7
         obs_array = obs_array .and. input % chn(7) % refl >=0.
      case ( 3)
         CHN_NIR_DEFAULT = 20
         obs_array = obs_array .and. input % chn(20) % rad >=0.
      case default
         CHN_NIR_DEFAULT = 0
   end select


   obs_and_acha_array =  obs_array .and. input % cloud_temp  > 10

   cloud_array =  obs_and_acha_array &
      & .and. ( input % cloud_mask == EM_cloud_mask % CLOUDY &
      & .or. input % cloud_mask == EM_cloud_mask % PROB_CLOUDY )

   water_phase_array = input % cloud_type  == EM_cloud_type % FOG &
      &  .or. input % cloud_type  == EM_cloud_type % WATER &
      &  .or. input % cloud_type == EM_cloud_type % SUPERCOOLED &
      &  .or. input % cloud_type == EM_cloud_type % MIXED

   ! - define output
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

   ! - initialize
   info_flag  = 0_1

   if ( .not. input % is_channel_on (CHN_NIR_DEFAULT) ) then
      print*, 'dcomp NIR channel not set! ==> MODIS equaivalant channel: ', CHN_NIR
      print*, 'all dcomp results are set to missing values'
      return
   end if

   if ( .not. input % is_channel_on (CHN_VIS_DEFAULT) ) then
      print*, 'dcomp VIS channel not set! ==> MODIS equaivalant channel: ', CHN_VIS
      print*, 'all dcomp results are set to missing values'
      return
   end if

   dcomp_mode = input % mode



  ! allocate ( rain_column (nr_elem, nr_lines))
  ! call compute_rain_column ( input % cloud_temp % d, rain_column)

   line_loop: do line_idx = 1 , nr_lines
      elem_loop: do elem_idx = 1,   nr_elem

         if ( .not. cloud_array (elem_idx,line_idx)  ) cycle elem_loop

         ! - set aliases
         cld_press  =  input % cloud_press  (elem_idx,line_idx)
         cld_temp   =  input % cloud_temp  (elem_idx,line_idx)
         sol_zen    =  input % sol  (elem_idx,line_idx)
         sat_zen    =  input % sat  (elem_idx,line_idx)
         rel_azi    =  input % azi  (elem_idx,line_idx)
         ozone_dobson = input % ozone_nwp  (elem_idx,line_idx)



         ! - compute transmission

         loop_chn: do chn_idx = 1 , 40

            if ( .not. input % is_channel_on (chn_idx)) cycle  loop_chn
            ozone_coeff(:) = 0.
            if (chn_idx == 1) ozone_coeff = OZONE_COEFF_CHN1
            gas_coeff(:) = 0.
            if (input % gas_coeff(chn_idx) % is_set) then
                gas_coeff = input % gas_coeff ( chn_idx) % d
            endif

            call trans_atm_above_cloud ( &
               input % tpw_ac (elem_idx,line_idx) &
               , ozone_dobson &
               , input % press_sfc  (elem_idx,line_idx) &
               , cld_press &
               , air_mass_array(elem_idx,line_idx) &
               , gas_coeff, ozone_coeff, 0.044 &
               , trans_total(chn_idx) &
               , trans_unc_total(chn_idx) &
               )

            refl_toc( chn_idx ) = refl_toa  /  trans_total (chn_idx )

            alb_sfc( chn_idx ) =  ( input % chn( chn_idx) % alb_sfc  (elem_idx,line_idx) ) / 100.

            alb_sfc( chn_idx ) = max ( alb_sfc( chn_idx ) , ALBEDO_OCEAN (chn_idx) )

            alb_unc_sfc  (chn_idx) = 0.05

            if ( chn_idx == 20 ) then

               trans_total (chn_idx) = input % chn( chn_idx) % trans_ac_nadir   (elem_idx, line_idx)
               rad_to_refl_factor = PI / cos ( sol_zen * PI / 180.) &
                  / ( input % solar_irradiance (chn_idx) &
                  / input % sun_earth_dist ** 2 )
               refl_toc( chn_idx ) = input % chn(chn_idx) % rad (elem_idx, line_idx) * rad_to_refl_factor

               rad_clear_sky_toc_ch20 = input % chn( chn_idx) % rad_clear_sky_toc  (elem_idx, line_idx)
               rad_clear_sky_toa_ch20 = input % chn( chn_idx) % rad_clear_sky_toa (elem_idx, line_idx)

            end if

         end do loop_chn

         ! determine which channels are used..
         !  Snow is DCOMP mode 4 if 1.6 and 3.75 are available
         CHN_VIS = CHN_VIS_DEFAULT
         CHN_NIR = CHN_NIR_DEFAULT
         dcomp_mode = input % mode

         if ( input % snow_class (elem_idx, line_idx) == 3 ) then
            ! check if channels are available
            if ( input % is_channel_on ( 6) .and. input % is_channel_on ( 20)) then
              IF (input % chn(6) % refl  (elem_idx, line_idx) .gt. 0) then
               CHN_VIS =  6
               CHN_NIR = 20
               dcomp_mode = 4
              end if
            end if
         end if


         ! - vis
         obs_vec ( 1 ) = input % chn( CHN_VIS) % refl  (elem_idx, line_idx) / 100.
         obs_unc ( 1 ) = trans_unc_total ( CHN_VIS)  +calib_err (CHN_VIS)

         alb_vec ( 1 ) =  alb_sfc ( CHN_VIS)

         alb_vec(3)   =  input % chn( 1 ) % alb_sfc_dark_sky (elem_idx, line_idx) /100.

         alb_unc ( 1) = 0.05
         trans_vec ( 1) = trans_total ( CHN_VIS )

         ! - nir channel

         if ( CHN_NIR == 20 ) then
            obs_vec( 2 ) = refl_toc ( CHN_NIR )
            obs_unc( 2 ) = obs_vec( 2 ) * 0.1
         else
            obs_vec( 2 ) = input % chn( CHN_NIR) % refl (elem_idx, line_idx) /100.
            obs_unc( 2 ) = max ( trans_unc_total  ( CHN_NIR ) , 0.01 )  + calib_err (CHN_NIR)
         end if


         alb_vec( 2 ) = alb_sfc ( CHN_NIR)
         alb_unc( 2 ) = 0.05

         trans_vec ( 2) = trans_total ( CHN_NIR )

         ! - apriori
         state_apriori (1) = 0.7 * ( 100. * obs_vec(1) ) ** (0.9)
         state_apriori(1) = log10 ( max ( 0.1 , state_apriori(1) ) )

         ! use acha as apriori if cirrus
         if ( input % tau_acha(elem_idx, line_idx) .gt. 0.01 &
              .and. input % tau_acha(elem_idx, line_idx) .lt. 8.01 &
              .and. input % cloud_type(elem_idx, line_idx) .eq. 7  ) then

            state_apriori(1) = log10 ( max ( 0.001 , input % tau_acha(elem_idx, line_idx) ) )

         end if

         state_apriori (2) = 1.3
         if  (water_phase_array ( elem_idx, line_idx) ) state_apriori(2) = 1.0



         call dcomp_algorithm ( &
            &   obs_vec &
            & , obs_unc  &
            & , alb_vec &
            & , alb_unc  &
            & , state_apriori &
            & , trans_vec  &
            & , sol_zen &
            & , sat_zen &
            & , rel_azi &
            & , cld_temp &
            & , water_phase_array ( elem_idx, line_idx) &
            & , rad_clear_sky_toc_ch20 &
            & , rad_clear_sky_toa_ch20 &
            & , trim(sensorname_from_wmoid(input % sensor_wmo_id)) &
            & , dcomp_out &
            & , dcomp_mode &
            & , input % lut_path &
            & , debug_mode  )
        debug_mode = 0
    ! if ( dcomp_out% cod .lt. 0 ) debug_mode = 4
         if ( debug_mode == 4  ) then
            print*,'=======================> input:',CHN_NIR
            print*,'Elem Line: ', elem_idx,line_idx
            print*,' Obs vector: ',obs_vec
            print*,' Obs uncert: ',obs_unc
            print*,' Surface Albedo: ', alb_vec
            print*,' Surface albedo unc: ',alb_unc
            print*,' Transmission: ',trans_vec
            print*, 'Angles: ',sol_zen,sat_zen,rel_azi
            print*, 'Cloud temp mask: ',cld_temp,water_phase_array( elem_idx, line_idx)
            print*, 'Ch20 rtm: ', rad_clear_sky_toc_ch20 , rad_clear_sky_toa_ch20
            print*, 'dcomp mode:', input % mode
            print*, 'output: '
            print*, dcomp_out % cod, dcomp_out % cps
            print*,'ozone uncert ch1.:',trans_unc_total ( CHN_VIS)
            print*, 'calib error :',calib_err (CHN_VIS)
            print*,'snow?:',input % snow_class  ( elem_idx, line_idx)
            print*

            print*,'==============================='

         end if

         output % cod % d (elem_idx,line_idx) = dcomp_out % cod
         output % cps % d (elem_idx, line_idx) = dcomp_out % cps

         output % cod_unc % d ( elem_idx, line_idx) = dcomp_out % codu
         output % ref_unc % d ( elem_idx, line_idx) = dcomp_out % cpsu
         output % cld_trn_sol % d ( elem_idx, line_idx) = dcomp_out % cloud_trans_sol_vis
         output % cld_trn_obs % d ( elem_idx, line_idx) = dcomp_out % cloud_trans_sat_vis
         output % cld_alb % d ( elem_idx, line_idx)  = dcomp_out % cloud_alb_vis
         output % cld_sph_alb % d ( elem_idx, line_idx)  = dcomp_out % cloud_sph_alb_vis
         output % refl_vis_max % d ( elem_idx, line_idx)  = dcomp_out % refl_vis_max
         ! - Water Path
         if ( dcomp_out % cod .gt. 0 .and. dcomp_out % cps .gt. 0) then
            output % lwp % d ( elem_idx, line_idx)  = 5./ 9. * dcomp_out % cod * dcomp_out % cps
            output % iwp % d ( elem_idx, line_idx)  = 5./ 9. * dcomp_out % cod * dcomp_out % cps
         end if

         ! -  precipitation

   !      call compute_precip ( output % cod % d ( elem_idx, line_idx) &
   !         , output % cps % d ( elem_idx, line_idx) &
   !         , rain_column ( elem_idx, line_idx) &
   !         , output % rain_probability % d ( elem_idx, line_idx) &
   !         , output % rain_rate % d ( elem_idx, line_idx) )

         ! - DCOMP_QF_PROCESSION B0
         quality_flag (elem_idx,line_idx) = ibset ( quality_flag(elem_idx,line_idx) , 0)

         if ( dcomp_out % statusOK ) then
            ! - DCOMP_QF_COD_VALID B1
            quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 1)
            ! - DCOMP_QF_REF_VALID B2
            quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 2)
            ! - initial DCOMP_QF_COD_DEGRADED B3
            quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 3)
            ! - initial DCOMP_QF_REF_DEGRADED B4
            quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 4)
            ! --DCOMP_QF_REF_CONVERGENCY B5
            quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 5)

         end if

      end do elem_loop
   end do   line_loop


  ! if (allocated( rain_column) ) deallocate ( rain_column )

   ! DCOMP_INFO_LAND_SEA I0
   where (input % is_land  )
      info_flag = ibset (  info_flag , 0 )
   end where

   ! -DCOMP_INFO_DAY_NIGHT I1
   where (input % sol   > 82. )
      info_flag = ibset ( info_flag , 1)
   end where

   ! - DCOMP_INFO_TWILIGHT I2
   where (input % sol  > 65. .and. input % sol  < 82.)
      info_flag = ibset ( info_flag , 2)
   end where

   ! - DCOMP_INFO_SNOW I3
   where ( input % snow_class  == EM_snow_class % SNOW )
      info_flag = ibset ( info_flag , 3)
   end where

   ! - DCOMP_INFO_SEA_ICE I4
   where ( input % snow_class  == EM_snow_class % SEA_ICE )
      info_flag = ibset ( info_flag , 4)
   end where

   ! -DCOMP_INFO_PHASE
   !ccm   where ( .not. water_phase_array)
   dim_1_w = size(water_phase_array,1)
   dim_2_w = size(water_phase_array,2)
   where ( .not. water_phase_array(1:dim_1_w,1:dim_2_w))
      info_flag(1:dim_1_w,1:dim_2_w) = ibset ( info_flag(1:dim_1_w,1:dim_2_w), 5)
   end where

   ! -DCOMP_INFO_THICK_CLOUD
   where ( output % cod % d > 80 )
      info_flag = ibset ( info_flag, 6)
   end where

   ! -DCOMP_INFO_THIN_CLOUD
   where ( output % cod % d < 4 .and. output % cod % d > 0 )
      info_flag = ibset ( info_flag, 7)
   end where

   ! - DCOMP_QF_COD_DEGRADED B3
   where ( (  btest(info_flag,2) &
      .or. btest(info_flag,3) &
      .or. btest(info_flag,4) &
      .or. btest(info_flag,5)  &
      .or. btest(info_flag,6)) &
      .and. btest(quality_flag,0))
      quality_flag = ibset ( quality_flag , 3 )
   end where

   ! - DCOMP_QF_REF_DEGRADED B4
   where (    ( btest(info_flag,2) &
      .or. btest(info_flag,3) &
      .or. btest(info_flag,4) &
      .or. btest(info_flag,5)  &
      .or. btest(info_flag,7)) &
      &  .and. btest(quality_flag,0))
      quality_flag = ibset ( quality_flag , 4 )
   end where

   ! - compute lwp
   where ( water_phase_array(1:dim_1_w,1:dim_2_w) &
      .and. .not. btest ( quality_flag(1:dim_1_w,1:dim_2_w) , 1 )  &
      .and. .not. btest ( quality_flag(1:dim_1_w,1:dim_2_w) , 2 ) )
      output % lwp % d(1:dim_1_w,1:dim_2_w) =  output % cod % d(1:dim_1_w,1:dim_2_w) &
         & * output % cps % d(1:dim_1_w,1:dim_2_w) * 5.0 / 9.0
   end where

      ! - compute lwp
   where ( .not. water_phase_array(1:dim_1_w,1:dim_2_w) &
      .and. .not. btest ( quality_flag(1:dim_1_w,1:dim_2_w) , 1 )  &
      .and. .not. btest ( quality_flag(1:dim_1_w,1:dim_2_w) , 2 ) )
      output % iwp % d(1:dim_1_w,1:dim_2_w) =  (output % cod % d(1:dim_1_w,1:dim_2_w)  ** (1/0.84) ) / 0.065
   end where

   where ( obs_array .and. .not. cloud_array )
      output % cld_trn_sol % d   =  1.0
      output % cld_trn_obs % d   =  1.0
      output % cld_alb % d       =  0.0
      output % cld_sph_alb % d   =  0.0
      output % cod % d           =  0.0
      output % cps % d           =  -999.0
   end where

   !  where ( output % cod % d .lt. (-1.)  )
   !     output % cld_trn_sol % d   =  -999.
   !     output % cld_trn_obs % d   =  -999.
   !     output % cld_alb % d       =  -999.
   !     output % cld_sph_alb % d   =  -999.
   !  end where


   output % quality % d = quality_flag
   output % info % d = info_flag

   ! compute successrate
   output % nr_obs = count ( obs_array)
   output % nr_clouds = count ( cloud_array)
   output % nr_success_cod = count (btest ( quality_flag , 1 ))
   output % nr_success_cps = count (btest ( quality_flag , 2 ))

   tried =  count (btest(quality_flag,0))

   output % successrate = 0.0
   if ( tried > 0 ) then
      success = count(.not. btest(quality_flag,1) .and. .not. btest ( quality_flag,2) )
      output % successrate = success / tried
   end if

   deallocate ( obs_array &
      ,  cloud_array, obs_and_acha_array  )
   deallocate ( water_phase_array )

   deallocate ( air_mass_array )

   output % version = '$Id: dcomp_array_loop_sub.f90 189 2018-03-09 18:49:39Z awalther $'




end subroutine dcomp_array_loop
