! $Header: https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/trunk/main_src/dncomp_clavrx_bridge_mod.f90 4093 2021-03-04 13:07:43Z awalther $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: dcomp_clavrx_bridge_mod.f90 (src)
!       dcomp_clavrx_bridge_mod (program)
!
! PURPOSE:
!
! DESCRIPTION:
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! REVISION HISTORY:
!   02/10.2013 : first version
!   10/21/2013 : bridge to array instead to pixel-based
!
!
!
!--------------------------------------------------------------------------------------
module dncomp_clavrx_bridge_mod


  ! -- MODULES USED

  use CONSTANTS_MOD, only: &
  REAL4 , INT4, INT2, INT1 , SYM &
  , MISSING_VALUE_REAL4 , PI &
  , dcomp_version   ! -- this is not a constant!

  !use rtm_common

  use dncomp_interface_def_mod , only: &
  dncomp_in_type &
  , dncomp_out_type &
  , alloc_dncomp &
  , n_chn

  use  clavrx_message_mod, only: &
  mesg

  use pixel_common_mod, only: &
  ch &
  , geo &
  , sfc &
  , sensor &
  , image &
  , cldmask &
  , acha &
  , cld_type &
  , bad_pixel_mask &
  , lwp_dcomp, reff_dcomp, tau_dcomp, iwp_dcomp &
  , tau_dcomp_1, tau_dcomp_2, tau_dcomp_3 &
  , reff_dcomp_1, reff_dcomp_2, reff_dcomp_3 &
  , tau_dcomp_qf, reff_dcomp_qf &
  , tau_dcomp_cost , reff_dcomp_cost &
  , dcomp_info_flag, dcomp_quality_flag &
  , dcomp_quality_flag_1, dcomp_quality_flag_2 &
  , dcomp_quality_flag_3 &
  , cloud_063um_transmission_solar &
  , cloud_063um_transmission_view &
  , cloud_063um_spherical_albedo &
  , cloud_063um_albedo &
  , ancil_data_dir &
  , ch &
  , dcomp_mode &
  , zen_idx_rtm &
  , solar_rtm &
  , tau_nlcomp &
  , reff_nlcomp &
  , tau_nlcomp_cost &
  , reff_nlcomp_cost &
  , nlcomp_quality_flag &
  , nlcomp_info_flag &
  , refl_asym_dcomp &
  , Static_Ref_065um_Dark_Composite

  !!! use pixel_common_mod, only: &
  !      dcomp_diag_2 , dcomp_diag_3 , dcomp_diag_4 &
  !      ,  dcomp_diag_toc_rfl1,  dcomp_diag_toc_rfl2 ,dcomp_diag_virt_alb1 &
  !      ,dcomp_diag_virt_alb2,dcomp_diag_wv1,dcomp_diag_wv2

  use calibration_constants_mod, only: &
  sun_earth_distance  &      !---- check, This is defined in three routines
  , solar_ch20_nu

  use dcomp_rtm_module

  implicit none

  private

  logical :: first_call = .true.

  public :: awg_cloud_dncomp_algorithm

  character(len = 120) :: DCOMP_RELEASE_VERSION = 'DCOMP version 2_0_0'


contains

  !----------------------------------------------------------------------
  !  AWG_CLOUD_DCOMP_ALGORITHM
  !    This is the DCOMP bridge from CLAVR-x
  !----------------------------------------------------------------------
  subroutine awg_cloud_dncomp_algorithm (  iseg_in , nlcomp_mode,  algorithm_started, version )

    implicit none

    !--- input
    integer, intent(in),optional :: iseg_in
    logical, intent(in),optional :: nlcomp_mode

    ! - output
    logical , intent(out) :: algorithm_started
    character(len = 120), intent(out), optional :: version

    type(dcomp_rtm_type), target :: dcomp_rtm
    type(dncomp_in_type)  :: dcomp_input
    type(dncomp_out_type) :: dncomp_output

    integer :: debug_mode
    integer :: dim_1, dim_2

    logical :: run_nlcomp

    integer, allocatable :: possible_channels ( : )
    logical :: chan_on ( N_CHN ) = .false.
    integer :: i, i_mode
    integer :: CHN_VIS
    integer :: CHN_NIR

    integer :: dcomp_mode_local
    character ( len = 12) :: string_welc

    logical, allocatable,target :: is_land (:,:)
    logical, allocatable, target :: is_valid (:,:)
    real, allocatable, target :: sfc_dummy(:,:)


    interface
      subroutine dcomp_array_loop (a , b , debug_mode_user)
        import dncomp_in_type
        import dncomp_out_type
        type (dncomp_in_type) , intent(in) :: a
        type (dncomp_out_type), intent(inout) :: b
        integer , intent(in), optional :: debug_mode_user
      end subroutine

    end interface


    interface
      subroutine nlcomp_array_loop_sub (a , b , debug_mode_user)
        import dncomp_in_type
        import dncomp_out_type
        type (dncomp_in_type) , intent(in) :: a
        type (dncomp_out_type), intent(out) :: b
        integer , intent(in), optional :: debug_mode_user
      end subroutine

    end interface
    
    ! ----- executable  --------------------------------------------------- !
    run_nlcomp = .false.
    if (present(nlcomp_mode)) run_nlcomp = nlcomp_mode
    if ( present ( version) ) version = DCOMP_RELEASE_VERSION

    algorithm_started = .false.

    ! - do we need to run dcomp at all? ( night  etc..)

    if (run_nlcomp) then
      ! add here all conditions which leads to a immediate stop
      if ( count (ch(44)%Ref_Lunar_Toa > 0) < 1 ) return
    else
      if ( count ( geo % solzen < 75. .and. geo % solzen >= 0 .and. geo % satzen < 75. ) < 1 ) return
    end if
    algorithm_started = .true.


    if ( first_call) then
      string_welc = 'DCOMP starts'
      if ( run_nlcomp ) string_welc = 'NLCOMP starts'
      call mesg (trim(string_welc),level = 5)
      first_call = .false.
    end if


    ! - compute DCOMP related RTM
    call perform_rtm_dcomp ( dcomp_rtm )

    dim_1 = Image%Number_Of_Elements
    dim_2 = Image%Number_Of_Lines_Read_This_Segment

    chan_on = .false.

    if (run_nlcomp) then
      allocate (possible_channels(2))
      possible_channels =[20,44]
    else
      allocate ( possible_channels(5))
      possible_channels = [ 1, 5, 6, 7, 20 ]
    end if

    do i = 1 , size ( possible_channels )
      if ( sensor % chan_on_flag_default ( possible_channels ( i) ) == 1 ) then
        chan_on (possible_channels ( i)  )  = .true.
      end if
    end do

    dcomp_input = dncomp_in_type (  chan_on )
    dcomp_input % is_channel_on = chan_on
    ! - ancil/lut path
    dcomp_input % lut_path = trim(ancil_data_dir)//"/static/luts/cld/"
    ! - wmo sensor id
    dcomp_input % sensor_wmo_id = sensor % wmo_id
    dcomp_input % sun_earth_dist = sun_earth_distance

    dcomp_input % gas_coeff(:) % is_set = .false.
    ! - all reflectance channels
    do i = 1, 19
      if ( dcomp_input % is_channel_on (i)) then
        dcomp_input % chn(i) %  refl  => ch(i)%ref_toa(1:dim_1,1:dim_2)
        dcomp_input % chn(i) % alb_sfc => ch(i) % sfc_ref_white_sky(1:dim_1,1:dim_2)
        dcomp_input % gas_coeff(i) % d = solar_rtm % tau_h2o_coef(i,:)
        dcomp_input % gas_coeff(i) % is_set = .true.
      end if
    end do

    dcomp_input % chn(1) % alb_sfc_dark_sky => Static_Ref_065um_Dark_Composite(1:dim_1,1:dim_2)

    if ( dcomp_input % is_channel_on (20)) then
      dcomp_input % chn(20) %  rad => ch(20)%rad_toa(1:dim_1,1:dim_2)
      if ( allocated(sfc_dummy)) deallocate(sfc_dummy)
      allocate(sfc_dummy(1:dim_1,1:dim_2))
      sfc_dummy = 100.0*(1.0 - ch(20)%sfc_emiss(1:dim_1,1:dim_2))
      dcomp_input % chn(20) % alb_sfc => sfc_dummy  ! Already subset
      dcomp_input % chn(20) % emiss_sfc  => ch(20)%sfc_emiss(1:dim_1,1:dim_2)
      dcomp_input % chn(20) % trans_ac_nadir  => dcomp_rtm % trans_ir_ac_nadir(1:dim_1,1:dim_2)
      dcomp_input % chn(20) % rad_clear_sky_toc  => dcomp_rtm % rad_clear_sky_toc_ch20(1:dim_1,1:dim_2)
      dcomp_input % chn(20) % rad_clear_sky_toa  => dcomp_rtm % rad_clear_sky_toa_ch20(1:dim_1,1:dim_2)
      ! -- Solar irradiance in channel 20
      dcomp_input % solar_irradiance ( 20) = solar_ch20_nu
    end if

    ! IR channels
    do i = 21, 44
      if ( dcomp_input % is_channel_on (i)) dcomp_input % chn(i) %  rad  => ch(i)%rad_toa(1:dim_1,1:dim_2)
    end do

    dcomp_input % sat => geo % satzen(1:dim_1,1:dim_2)
    dcomp_input % sol  => geo % solzen(1:dim_1,1:dim_2)
    dcomp_input % azi  => geo % relaz(1:dim_1,1:dim_2)

    ! - Cloud products
    dcomp_input % cloud_press  => acha % pc(1:dim_1,1:dim_2)
    dcomp_input % cloud_temp  => acha % tc(1:dim_1,1:dim_2)
    dcomp_input % tau_acha   => acha % tau(1:dim_1,1:dim_2)
    dcomp_input % cloud_mask  => cldmask % cld_mask(1:dim_1,1:dim_2)
    !ccm
    dcomp_input % cloud_type  => cld_type(1:dim_1,1:dim_2)
    !end ccm
    !ccm      dcomp_input % cloud_type % d  = cld_type

    ! - Flags
    if ( allocated(is_land)) deallocate(is_land)
    if ( allocated(is_valid)) deallocate(is_valid)
    allocate(is_land(dim_1,dim_2))
    is_land = ( sfc%land_mask(1:dim_1,1:dim_2) == 1 )
    allocate(is_valid(dim_1,dim_2))
    is_valid = ( bad_pixel_mask(1:dim_1,1:dim_2) /= 1 )
    dcomp_input % is_land => is_land    ! Already subset to (1:dim_1,1:dim_2)
    dcomp_input % is_valid => is_valid  !

    dcomp_input % press_sfc  =>  dcomp_rtm % sfc_nwp(1:dim_1,1:dim_2)
    dcomp_input % snow_class => sfc % snow(1:dim_1,1:dim_2)

    ! - Atmospheric contents
    ! ozone column in Dobson
    dcomp_input % ozone_nwp  => dcomp_rtm % ozone_path(1:dim_1,1:dim_2)
    ! Total water Vapour above the cloud
    dcomp_input % tpw_ac  => dcomp_rtm % tpw_ac(1:dim_1,1:dim_2)

    if ( .not. run_nlcomp) then
      do i_mode = 1 , 3
        dcomp_mode_local = i_mode

        if ( dcomp_mode .ne. i_mode .and. dcomp_mode .ne. 9 ) cycle


        !- check mode
        CHN_VIS = 1
        select case ( dcomp_mode_local )
        case ( 1 )
          CHN_NIR = 6
        case ( 2 )
          CHN_NIR = 7
        case ( 3 )
          CHN_NIR = 20
        case default
          print*, 'dcomp mode ',dcomp_mode_local,' not possible'
          return
        end select

        if ( .not. dcomp_input % is_channel_on (CHN_NIR)) then
          if ( iseg_in == 1 ) then
            print*,'dcomp NIR channel is not set! ==> MODIS equaivalant channel: ', CHN_NIR
            call mesg ( 'all dcomp results are set to missing values', color=41 , level = -1 )
          end if
          cycle
        end if


        if ( dcomp_input % is_channel_on (CHN_VIS) .eqv. .false.) then
          if ( iseg_in == 1 ) then
            print*, 'dcomp VIS channel is not set! ==> MODIS equaivalant channel: ', CHN_VIS
            call mesg ( 'all dcomp results are set to missing values', color=41 , level = -1 )
          end if
          cycle
        end if

        ! == CONFIGURE

        ! - dcomp-mode
        dcomp_input % mode = dcomp_mode_local

        debug_mode = 0

        call dcomp_input % check_input (debug_mode)

        call mesg ('DCOMP starts in Bridge',level = 9)

        call dcomp_array_loop ( dcomp_input , dncomp_output , debug_mode_user = debug_mode)

        call mesg ('DCOMP ends in Bridge',level = 9)

        tau_dcomp (1:dim_1,1:dim_2)   = dncomp_output % cod % d(1:dim_1,1:dim_2)
        reff_dcomp  (1:dim_1,1:dim_2) = dncomp_output % cps % d(1:dim_1,1:dim_2)


        select case (dcomp_mode_local)
        case (1)
          tau_dcomp_1 (1:dim_1,1:dim_2)   = dncomp_output % cod % d(1:dim_1,1:dim_2)
          reff_dcomp_1  (1:dim_1,1:dim_2) = dncomp_output % cps % d(1:dim_1,1:dim_2)
          dcomp_quality_flag_1(1:dim_1,1:dim_2) = dncomp_output % quality % d(1:dim_1,1:dim_2)

        case(2)

          tau_dcomp_2 (1:dim_1,1:dim_2)   = dncomp_output % cod % d(1:dim_1,1:dim_2)
          reff_dcomp_2  (1:dim_1,1:dim_2) = dncomp_output % cps % d(1:dim_1,1:dim_2)
          dcomp_quality_flag_2(1:dim_1,1:dim_2) = dncomp_output % quality % d(1:dim_1,1:dim_2)

        case(3)
          tau_dcomp_3 (1:dim_1,1:dim_2)   = dncomp_output % cod % d(1:dim_1,1:dim_2)
          reff_dcomp_3  (1:dim_1,1:dim_2) = dncomp_output % cps % d(1:dim_1,1:dim_2)
          dcomp_quality_flag_3(1:dim_1,1:dim_2) = dncomp_output % quality % d(1:dim_1,1:dim_2)

        end select

        lwp_dcomp (1:dim_1,1:dim_2)   = dncomp_output % lwp % d(1:dim_1,1:dim_2)
        iwp_dcomp (1:dim_1,1:dim_2)   = dncomp_output % iwp % d(1:dim_1,1:dim_2)

        tau_dcomp_cost(1:dim_1,1:dim_2)     = dncomp_output % cod_unc % d(1:dim_1,1:dim_2)
        reff_dcomp_cost(1:dim_1,1:dim_2)    = dncomp_output % ref_unc % d(1:dim_1,1:dim_2)
        dcomp_quality_flag(1:dim_1,1:dim_2) = dncomp_output % quality % d(1:dim_1,1:dim_2)
        dcomp_info_flag(1:dim_1,1:dim_2)    = dncomp_output % info % d(1:dim_1,1:dim_2)

        cloud_063um_transmission_solar(1:dim_1,1:dim_2) = dncomp_output % cld_trn_sol % d(1:dim_1,1:dim_2)
        cloud_063um_transmission_view(1:dim_1,1:dim_2)  = dncomp_output % cld_trn_obs % d(1:dim_1,1:dim_2)
        cloud_063um_albedo(1:dim_1,1:dim_2)             = dncomp_output % cld_alb % d(1:dim_1,1:dim_2)
        cloud_063um_spherical_albedo(1:dim_1,1:dim_2)   = dncomp_output % cld_sph_alb % d(1:dim_1,1:dim_2)

        refl_asym_dcomp  (1:dim_1,1:dim_2)   =   dncomp_output % refl_vis_max  % d(1:dim_1,1:dim_2)

        DCOMP_RELEASE_VERSION = dncomp_output % version

      end do

    else  ! nlcomp

      if ( dcomp_input % is_channel_on (44)) then
        dcomp_input % chn(44) %  refl  => ch(44)%ref_lunar_toa
        dcomp_input % chn(44) % alb_sfc  => ch(1)%sfc_ref_white_sky
      end if
      dcomp_input % zen_lunar  => geo % lunzen
      dcomp_input % azi_lunar  => geo % lunrelaz

      call nlcomp_array_loop_sub (dcomp_input, dncomp_output, debug_mode_user = 9) !debug_mode_user )

      tau_nlcomp (1:dim_1,1:dim_2)   = dncomp_output % cod % d
      reff_nlcomp  (1:dim_1,1:dim_2) = dncomp_output % cps % d
      tau_nlcomp_cost(1:dim_1,1:dim_2) = dncomp_output % cod_unc % d
      reff_nlcomp_cost(1:dim_1,1:dim_2) = dncomp_output % ref_unc % d
      nlcomp_quality_flag(1:dim_1,1:dim_2) = dncomp_output %  quality % d
      nlcomp_info_flag(1:dim_1,1:dim_2) = dncomp_output % info % d

    end if


    call dcomp_rtm % deallocate_it()


  end subroutine awg_cloud_dncomp_algorithm




end module dncomp_clavrx_bridge_mod
