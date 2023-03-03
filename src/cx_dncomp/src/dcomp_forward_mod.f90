
module dcomp_forward_mod

  use dcomp_lut_mod , only : &
  lut_obj , lut_output


  use dcomp_science_tools_mod, only: &
  get_planck_radiance_39um  &
  , get_rad_refl_factor

  private
  public :: thick_cloud_cps
<<<<<<< HEAD
=======
  public :: thin_cloud_mode_3
>>>>>>> dcomp_mode3_test


  integer, parameter, public:: real4 = selected_real_kind(6,37)

  real , parameter :: PI = 4.* ATAN(1.)
  real , parameter :: DTOR = PI / 180.

  type, public :: pixel_vec
    real :: sol_zen
    real :: sat_zen
    real :: rel_azi
    real :: ctt
    logical :: is_water_phase
  end type pixel_vec

  public :: dcomp_forward_computation


  real , save :: planck_rad
  real , save :: rad_to_refl

contains

  subroutine dcomp_forward_computation ( &
    state_vec  &   ! -input
    , pixel &
    , sensor &
    , channels &
    , alb_sfc &
    , air_trans_ac &
    , fm_vec &         ! - output
    , cld_trns_sol &
    , cld_trns_sat &
    , cld_sph_alb &
    , kernel &
    , rad_abv_cld &
    , rad_clear_toc  &
    , lut_path &
    , cld_albedo_vis &
    , rfl_vis_max )
<<<<<<< HEAD
=======




    implicit none

    real , intent(in) :: state_vec (2)
    type ( pixel_vec ) , intent ( in ) :: pixel
    integer , intent ( in ) :: channels ( :)
    real,  intent ( in ) :: air_trans_ac ( : )
    character ( len = * ) , intent (in ) :: sensor
    real, intent(in) ::  alb_sfc (2)
    real, optional, intent ( in ) ::   rad_abv_cld
    real, optional, intent ( in ) ::   rad_clear_toc
    real, optional, intent ( out ) ::   rfl_vis_max
    character ( len = 1024 ) , intent ( in ) , optional :: lut_path

    real , intent ( out ) :: fm_vec ( 2 )
    real , intent ( out ) :: kernel( 2 , 2 )
    real , intent ( out ) :: cld_trns_sol( 2 )
    real , intent ( out ) :: cld_trns_sat( 2 )
    real , intent ( out ) :: cld_sph_alb( 2 )
    real , intent ( out ) , optional :: cld_albedo_vis


    integer :: phase_num

    real :: fm_nir_terr
    real :: kernel_nir_terr_cod
    real :: kernel_nir_terr_cps
    real :: trans_two_way
    real :: air_mass_two_way
    real :: air_mass_sat
    real :: trans_abv_cld

    type ( lut_output) :: lut_data
    integer :: i_channel
    integer, parameter :: MODIS_EQUIVALENT_39_CHN = 20
    real :: Alb_Sfc_Term

    real :: rfl_thick (9)
    real :: ems_thick (9)



    ! - executable


    ! initialize lut object

    call lut_obj % initialize ( sensor , ancil_path = lut_path)


    call lut_obj % set_angles ( pixel % sat_zen , pixel % sol_zen, pixel % rel_azi )

    phase_num = 2
    if ( pixel % is_water_phase ) phase_num = 1

    air_mass_two_way = ( 1. / cos (pixel % Sol_zen * DTOR ) + 1./ cos ( pixel % Sat_zen * DTOR ) )

    do i_channel = 1 , size ( channels )

      call lut_obj % get_data ( channels ( i_channel )  , phase_num , state_vec(1), state_vec(2) , lut_data)

      cld_trns_sol (i_channel ) = lut_data % trn_sol
      cld_trns_sat (i_channel ) = lut_data % trn_sat
      cld_sph_alb (i_channel )  = lut_data % albsph

      if ( present ( cld_albedo_vis ) .and. i_channel == 1 ) then
        cld_albedo_vis =  lut_data % alb
      end if

      Alb_Sfc_Term = max (0. , Alb_Sfc( i_channel ) / (1.0 - Alb_Sfc( i_channel ) * lut_data % albsph ))
      fm_vec(i_channel) = lut_data% Refl + Alb_Sfc_Term * lut_data%Trn_sol * lut_data%Trn_sat

      kernel(i_channel,1) = lut_data%dRefl_dcod &
      & +  Alb_Sfc_Term * lut_data%Trn_sol *  lut_data % dTrans_sat_dcod &
      & + Alb_Sfc_Term * lut_data%Trn_sat  * lut_data%Dtrans_sol_Dcod &
      & +  ((lut_data%Trn_sol * lut_data%Trn_sat  * Alb_Sfc( i_channel ) * Alb_Sfc( i_channel ) &
      &   * lut_data%Dsph_alb_Dcod) &
      /(( 1 - Alb_Sfc( i_channel ) * lut_data%albsph)**2))

      kernel(i_channel,2) = lut_data%dRefl_dcps &
      & +  Alb_Sfc_Term * lut_data%Trn_sol *  lut_data % dTrans_sat_dcps &
      & + Alb_Sfc_Term * lut_data%Trn_sat  * lut_data%Dtrans_sol_Dcps &
      & +  ((lut_data%Trn_sol * lut_data%Trn_sat  * Alb_Sfc( i_channel ) * Alb_Sfc( i_channel ) &
      & * lut_data%Dsph_alb_Dcps) &
      /(( 1 - Alb_Sfc( i_channel ) * lut_data%albsph)**2))

      trans_two_way = air_trans_ac( i_channel ) ** air_mass_two_way
      fm_vec(i_channel) = fm_vec(i_channel) * trans_two_way

      kernel(i_channel,1) = kernel(i_channel,1) * trans_two_way
      kernel(i_channel,2) = kernel(i_channel,2) * trans_two_way

      ! - the only emis channel is channel 20
      if ( channels ( i_channel ) == MODIS_EQUIVALENT_39_CHN ) then
        planck_rad = get_planck_radiance_39um ( pixel % ctt  , trim(sensor) )

        rad_to_refl = get_rad_refl_factor ( trim ( sensor ) , pixel % sol_zen )
        air_mass_sat =  1./ cos ( pixel % Sat_zen * dtor)
        trans_abv_cld = air_trans_ac(i_channel) ** air_mass_sat

        fm_nir_terr =  rad_to_refl * (rad_abv_cld * lut_data % ems &
        + trans_abv_cld * lut_data %  ems * planck_rad &
        + lut_data % trn_ems * rad_clear_toc )

        kernel_nir_terr_cod = rad_to_refl * (rad_abv_cld * lut_data %  Dems_Dcod &
        + trans_abv_cld * Planck_Rad * lut_data % Dems_Dcod  &
        + rad_clear_toc *lut_data % Dtrnems_Dcod)

        kernel_nir_terr_cps  =  rad_to_refl * (rad_abv_cld * lut_data % Dems_Dcps &
        + trans_abv_cld * Planck_Rad * lut_data % Dems_Dcps  &
        + rad_clear_toc  *lut_data %Dtrnems_Dcps )

        fm_vec(i_channel) = fm_vec(i_channel) + fm_nir_terr

        kernel ( i_channel, 1) = kernel ( i_channel , 1)   + kernel_nir_terr_cod
        kernel ( i_channel, 2) = kernel ( i_channel , 2 )  + kernel_nir_terr_cps

      end if
    end do

    if (present(rfl_vis_max)) then


      call lut_obj % thick_cloud_rfl ( channels ( 1) , phase_num, rfl_thick, ems_thick)

      rfl_vis_max = rfl_thick(5)

    end if

  end subroutine dcomp_forward_computation

 !
 !
 !

  function thin_cloud_mode_3 (rfl_vis_obs, pixel) result(cod)
    implicit none
     real,  intent(in) :: rfl_vis_obs
     type ( pixel_vec) , intent(in) :: pixel
     integer :: phase_num
    real :: rfl_lut(29)
    integer :: cps_idx
    real :: cod
    integer :: cod_pos
    real :: cod_vec_lut (29)   !- CPS vector in LUT
    integer :: ii, i
  real     :: cod_wgt ! position weights in LUT
     cod_vec_lut = (/(i,i=-6,22,1)/)/10.
     phase_num = 2
     cps_idx = 5
     if ( pixel % is_water_phase ) then
       phase_num = 1
       cps_idx = 4
     end if
     call lut_obj % rfl_for_given_cps (  phase_num, cps_idx, rfl_lut)

    

    cod_pos = -999
    cod = -999.
    if (rfl_vis_obs .gt. 0. .and. rfl_vis_obs .lt. rfl_lut(1) ) cod_pos = 1
    do ii = 1, 28
      if ( (rfl_lut(ii)-rfl_vis_obs) .lt. 0 .and. (rfl_lut(ii+1)-rfl_vis_obs) .gt. 0 ) cod_pos = ii
    end do

    if (cod_pos .NE. -999) then
      cod_wgt = (rfl_lut(cod_pos) - rfl_vis_obs)/(rfl_lut(cod_pos) - rfl_lut(cod_pos+1))
      cod_wgt = max(cod_wgt,0.)
      cod = cod_Vec_LUT(cod_pos) + (cod_wgt * 0.1)
    end if
>>>>>>> dcomp_mode3_test




<<<<<<< HEAD
    implicit none

    real , intent(in) :: state_vec (2)
    type ( pixel_vec ) , intent ( in ) :: pixel
    integer , intent ( in ) :: channels ( :)
    real,  intent ( in ) :: air_trans_ac ( : )
    character ( len = * ) , intent (in ) :: sensor
    real, intent(in) ::  alb_sfc (2)
    real, optional, intent ( in ) ::   rad_abv_cld
    real, optional, intent ( in ) ::   rad_clear_toc
    real, optional, intent ( out ) ::   rfl_vis_max
    character ( len = 1024 ) , intent ( in ) , optional :: lut_path

    real , intent ( out ) :: fm_vec ( 2 )
    real , intent ( out ) :: kernel( 2 , 2 )
    real , intent ( out ) :: cld_trns_sol( 2 )
    real , intent ( out ) :: cld_trns_sat( 2 )
    real , intent ( out ) :: cld_sph_alb( 2 )
    real , intent ( out ) , optional :: cld_albedo_vis


    integer :: phase_num

    real :: fm_nir_terr
    real :: kernel_nir_terr_cod
    real :: kernel_nir_terr_cps
    real :: trans_two_way
    real :: air_mass_two_way
    real :: air_mass_sat
    real :: trans_abv_cld

    type ( lut_output) :: lut_data
    integer :: i_channel
    integer, parameter :: MODIS_EQUIVALENT_39_CHN = 20
    real :: Alb_Sfc_Term

    real :: rfl_thick (9)
    real :: ems_thick (9)



    ! - executable


    ! initialize lut object

    call lut_obj % initialize ( sensor , ancil_path = lut_path)


    call lut_obj % set_angles ( pixel % sat_zen , pixel % sol_zen, pixel % rel_azi )
=======


  end function thin_cloud_mode_3


  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  function thick_cloud_cps ( rfl_nir_obs, channel_nir, pixel , dcomp_mode ) result ( cps)
    implicit none

    real,  intent(in) :: rfl_nir_obs
    type ( pixel_vec) , intent(in) :: pixel
    integer , intent ( in ) :: channel_nir
    integer , intent ( in ) :: dcomp_mode
    real :: cps

    ! -- local

    integer  :: cps_pos          ! position in LUT
    real     :: cps_wgt ! position weights in LUT
    integer :: i
    real :: cps_vec_lut (9)   !- CPS vector in LUT
    real :: rad (9)
    real :: rfl_terr (9)

    real :: Ems_lut(9)  ! - sub table
    real :: rfl_lut(9)
    real :: Rfl (9)! - Reflectance at TOA under consideration of surface impact
    integer :: ii
    integer :: phase_num

    ! - executable

    ! - cps lut range [0.4,0.6,...,2.0]
    cps_vec_lut = (/(i,i=4,20,2)/)/10.
>>>>>>> dcomp_mode3_test

    phase_num = 2
    if ( pixel % is_water_phase ) phase_num = 1

<<<<<<< HEAD
    air_mass_two_way = ( 1. / cos (pixel % Sol_zen * DTOR ) + 1./ cos ( pixel % Sat_zen * DTOR ) )

    do i_channel = 1 , size ( channels )

      call lut_obj % get_data ( channels ( i_channel )  , phase_num , state_vec(1), state_vec(2) , lut_data)

      cld_trns_sol (i_channel ) = lut_data % trn_sol
      cld_trns_sat (i_channel ) = lut_data % trn_sat
      cld_sph_alb (i_channel )  = lut_data % albsph

      if ( present ( cld_albedo_vis ) .and. i_channel == 1 ) then
        cld_albedo_vis =  lut_data % alb
      end if

      Alb_Sfc_Term = max (0. , Alb_Sfc( i_channel ) / (1.0 - Alb_Sfc( i_channel ) * lut_data % albsph ))
      fm_vec(i_channel) = lut_data% Refl + Alb_Sfc_Term * lut_data%Trn_sol * lut_data%Trn_sat

      kernel(i_channel,1) = lut_data%dRefl_dcod &
      & +  Alb_Sfc_Term * lut_data%Trn_sol *  lut_data % dTrans_sat_dcod &
      & + Alb_Sfc_Term * lut_data%Trn_sat  * lut_data%Dtrans_sol_Dcod &
      & +  ((lut_data%Trn_sol * lut_data%Trn_sat  * Alb_Sfc( i_channel ) * Alb_Sfc( i_channel ) &
      &   * lut_data%Dsph_alb_Dcod) &
      /(( 1 - Alb_Sfc( i_channel ) * lut_data%albsph)**2))

      kernel(i_channel,2) = lut_data%dRefl_dcps &
      & +  Alb_Sfc_Term * lut_data%Trn_sol *  lut_data % dTrans_sat_dcps &
      & + Alb_Sfc_Term * lut_data%Trn_sat  * lut_data%Dtrans_sol_Dcps &
      & +  ((lut_data%Trn_sol * lut_data%Trn_sat  * Alb_Sfc( i_channel ) * Alb_Sfc( i_channel ) &
      & * lut_data%Dsph_alb_Dcps) &
      /(( 1 - Alb_Sfc( i_channel ) * lut_data%albsph)**2))

      trans_two_way = air_trans_ac( i_channel ) ** air_mass_two_way
      fm_vec(i_channel) = fm_vec(i_channel) * trans_two_way

      kernel(i_channel,1) = kernel(i_channel,1) * trans_two_way
      kernel(i_channel,2) = kernel(i_channel,2) * trans_two_way

      ! - the only emis channel is channel 20
      if ( channels ( i_channel ) == MODIS_EQUIVALENT_39_CHN ) then
        planck_rad = get_planck_radiance_39um ( pixel % ctt  , trim(sensor) )

        rad_to_refl = get_rad_refl_factor ( trim ( sensor ) , pixel % sol_zen )
        air_mass_sat =  1./ cos ( pixel % Sat_zen * dtor)
        trans_abv_cld = air_trans_ac(i_channel) ** air_mass_sat

        fm_nir_terr =  rad_to_refl * (rad_abv_cld * lut_data % ems &
        + trans_abv_cld * lut_data %  ems * planck_rad &
        + lut_data % trn_ems * rad_clear_toc )

        kernel_nir_terr_cod = rad_to_refl * (rad_abv_cld * lut_data %  Dems_Dcod &
        + trans_abv_cld * Planck_Rad * lut_data % Dems_Dcod  &
        + rad_clear_toc *lut_data % Dtrnems_Dcod)

        kernel_nir_terr_cps  =  rad_to_refl * (rad_abv_cld * lut_data % Dems_Dcps &
        + trans_abv_cld * Planck_Rad * lut_data % Dems_Dcps  &
        + rad_clear_toc  *lut_data %Dtrnems_Dcps )

        fm_vec(i_channel) = fm_vec(i_channel) + fm_nir_terr

        kernel ( i_channel, 1) = kernel ( i_channel , 1)   + kernel_nir_terr_cod
        kernel ( i_channel, 2) = kernel ( i_channel , 2 )  + kernel_nir_terr_cps

      end if
    end do

    if (present(rfl_vis_max)) then


      call lut_obj % thick_cloud_rfl ( channels ( 1) , phase_num, rfl_thick, ems_thick)

      rfl_vis_max = rfl_thick(5)

    end if

  end subroutine dcomp_forward_computation

  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  function thick_cloud_cps ( rfl_nir_obs, channel_nir, pixel , dcomp_mode ) result ( cps)
    implicit none

    real,  intent(in) :: rfl_nir_obs
    type ( pixel_vec) , intent(in) :: pixel
    integer , intent ( in ) :: channel_nir
    integer , intent ( in ) :: dcomp_mode
    real :: cps

    ! -- local

    integer  :: cps_pos          ! position in LUT
    real     :: cps_wgt ! position weights in LUT
    integer :: i
    real :: cps_vec_lut (9)   !- CPS vector in LUT
    real :: rad (9)
    real :: rfl_terr (9)

    real :: Ems_lut(9)  ! - sub table
    real :: rfl_lut(9)
    real :: Rfl (9)! - Reflectance at TOA under consideration of surface impact
    integer :: ii
    integer :: phase_num

    ! - executable

    ! - cps lut range [0.4,0.6,...,2.0]
    cps_vec_lut = (/(i,i=4,20,2)/)/10.

    phase_num = 2
    if ( pixel % is_water_phase ) phase_num = 1

    call lut_obj % thick_cloud_rfl ( channel_nir  , phase_num  , rfl_lut, ems_lut )

    rfl = rfl_lut
    if ( dcomp_mode == 3 ) then
      rad =  ems_lut * Planck_Rad
      rfl_terr = rad * Rad_To_Refl
      rfl = rfl +  rfl_terr
    end if

    cps_pos = -999
    cps = -999.

    do ii = 1, 8
      if ( rfl(ii)-rfl_nir_obs > 0 .and. rfl(ii+1)-rfl_nir_obs < 0 ) cps_pos = ii
    end do

    if (cps_pos .NE. -999) then
      cps_wgt = (rfl(cps_pos) - rfl_nir_obs)/(rfl(cps_pos) - rfl(cps_pos+1))
      cps = CPS_Vec_LUT(cps_pos) + (cps_wgt * 0.2)
    end if


    return


  end function thick_cloud_cps


=======
    call lut_obj % thick_cloud_rfl ( channel_nir  , phase_num  , rfl_lut, ems_lut )

    rfl = rfl_lut
    if ( dcomp_mode == 3 ) then
      rad =  ems_lut * Planck_Rad
      rfl_terr = rad * Rad_To_Refl
      rfl = rfl +  rfl_terr
    end if

    cps_pos = -999
    cps = -999.

    do ii = 1, 8
      if ( rfl(ii)-rfl_nir_obs > 0 .and. rfl(ii+1)-rfl_nir_obs < 0 ) cps_pos = ii
    end do

    if (cps_pos .NE. -999) then
      cps_wgt = (rfl(cps_pos) - rfl_nir_obs)/(rfl(cps_pos) - rfl(cps_pos+1))
      cps = CPS_Vec_LUT(cps_pos) + (cps_wgt * 0.2)
    end if


    return


  end function thick_cloud_cps


>>>>>>> dcomp_mode3_test
end module dcomp_forward_mod
