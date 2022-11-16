!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE
!
! NAME: tracer.f90 (src)
!
! PURPOSE: Transfer of data to and from external tools
!
!
! AUTHORS:
!  Coda Phillips, coda.phillips@wisc.edu
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
! REVISON HISTORY:
!    Creation Date Feb 2022
!--------------------------------------------------------------------------------------

      module tracer

        use pixel_common_mod, only: Sfc,Nav,Geo,Ch,CLDMASK,Tracer_Flag,Skip_Output,CCL, ACHA, &
            Zen_Idx_RTM, NWP_PIX, Tau_DCOMP, Tau_DCOMP_1, Tau_DCOMP_2, Tau_DCOMP_3, &
            Reff_DCOMP, Reff_DCOMP_1, Reff_DCOMP_2, REFF_DCOMP_3, &
            Cld_Type, Cld_Phase, DCOMP_Quality_Flag, &
            Insolation_DCOMP, Insolation_Diffuse_DCOMP, &
            Image, Temporary_Data_Dir, Tc_Opaque_Cloud, &
            Bad_Pixel_Mask, &
            mask_lrc, i_lrc, j_lrc, &
            Cwp_Dcomp
        use calibration_constants_mod, only: Planck_A1, Planck_A2, Planck_Nu, Sun_Earth_Distance
        use viirs_clavrx_bridge, only: viirs_out => out
        use viirs_nasa_read_module, only: nasa_viirs_i5_bt, nasa_viirs_i4_bt, &
              nasa_viirs_i1_ref, nasa_viirs_i2_ref, nasa_viirs_i3_ref
        use clavrx_static_nav_module, only: ABHI_1km => Output_Seg_1km, &
            ABHI_500m => Output_Seg_500m, ABHI_Rad_to_Ref_Fac => Rad_To_Ref_Fac
        use rtm_common_mod, only: RTM, NLevels_Rtm
        use rt_utilities_mod, only: RTM_NVZEN
        use constants_mod , only: SYM, NCHAN_CLAVRX
        use iso_c_binding, only: c_loc, c_ptr

        implicit none

        integer, parameter :: NUM_CHN_AHI = 16

        integer,parameter :: n_symbols_max = 512
        integer,parameter :: n_rtm_symbols_max = 64
        integer,parameter :: rtm_num_channels_max = 45
        integer,parameter :: symbol_length_max = 256
        integer,parameter :: max_clones = 64

        integer, dimension(max_clones) :: sibling_pids
        integer*8, volatile :: sibling_pids_addr

#ifdef __GFORTRAN__
#define DECLARENAMES(TYPE,DIM) character(LEN=symbol_length_max), dimension(n_symbols_max) :: symbol_names_/**/TYPE/**/_/**/DIM/**/d
#define DECLAREPTRS(TYPE,DIM) integer*8, dimension(n_symbols_max) :: symbol_ptrs_/**/TYPE/**/_/**/DIM/**/d
#define DECLARESHAPES(TYPE,DIM) integer*8, dimension(DIM, n_symbols_max) :: symbol_shapes_/**/TYPE/**/_/**/DIM/**/d
#define DECLAREHOOKS(TYPE,DIM) integer*8, volatile :: num_symbols_/**/TYPE/**/_/**/DIM/**/d, symbol_names_/**/TYPE/**/_/**/DIM/**/d_ptr, symbol_ptrs_/**/TYPE/**/_/**/DIM/**/d_ptr, symbol_shapes_/**/TYPE/**/_/**/DIM/**/d_ptr
#else
#endif
        ! 8-bit Int (Integer*1) -- 3 dimensions
        DECLARENAMES(i1,3)
        DECLAREPTRS(i1,3)
        DECLARESHAPES(i1,3)
        DECLAREHOOKS(i1,3)

        ! 32-bit Float (Real*4) -- 3 dimensions
        DECLARENAMES(r4,3)
        DECLAREPTRS(r4,3)
        DECLARESHAPES(r4,3)
        DECLAREHOOKS(r4,3)

        ! 32-bit Float (Real*4) -- 2 dimensions
        DECLARENAMES(r4,2)
        DECLAREPTRS(r4,2)
        DECLARESHAPES(r4,2)
        DECLAREHOOKS(r4,2)

        ! 8-bit Int (Integer*1) -- 2 dimensions
        DECLARENAMES(i1,2)
        DECLAREPTRS(i1,2)
        DECLARESHAPES(i1,2)
        DECLAREHOOKS(i1,2)

        ! 32-bit Int (Integer*4) -- 2 dimensions
        DECLARENAMES(i4,2)
        DECLAREPTRS(i4,2)
        DECLARESHAPES(i4,2)
        DECLAREHOOKS(i4,2)

        ! 64-bit Int (Integer*8) -- 2 dimensions
        DECLARENAMES(i8,2)
        DECLAREPTRS(i8,2)
        DECLARESHAPES(i8,2)
        DECLAREHOOKS(i8,2)

        ! 8-bit Int (Integer*1) -- 1 dimension
        DECLARENAMES(i1,1)
        DECLAREPTRS(i1,1)
        DECLARESHAPES(i1,1)
        DECLAREHOOKS(i1,1)

        ! 32-bit Int (Integer*4) -- 1 dimension
        DECLARENAMES(i4,1)
        DECLAREPTRS(i4,1)
        DECLARESHAPES(i4,1)
        DECLAREHOOKS(i4,1)

        ! 64-bit Int (Integer*8) -- 1 dimension
        DECLARENAMES(i8,1)
        DECLAREPTRS(i8,1)
        DECLARESHAPES(i8,1)
        DECLAREHOOKS(i8,1)

        ! 32-bit Float (Real*4) -- 1 dimension
        DECLARENAMES(r4,1)
        DECLAREPTRS(r4,1)
        DECLARESHAPES(r4,1)
        DECLAREHOOKS(r4,1)

        ! 32-bit Float (Real*4) -- 0 dimensions
        DECLARENAMES(r4,0)
        DECLAREPTRS(r4,0)
        DECLAREHOOKS(r4,0)

        integer*4, volatile :: number_of_lines
        integer*4, volatile :: wait_number
        integer*4, volatile :: skip_processing = -1


        ! RTM data structure is terrible to access
        ! make copies with more accessible structure
        integer*8, dimension(3) :: rtm_shape_3d

        integer*8, volatile :: rtm_numel
        integer*8, dimension(:), allocatable :: rtm_index_0, rtm_index_1,rtm_zen_index, rtm_sfc_index
        real*4, dimension(:,:), allocatable :: rtm_t_prof, rtm_z_prof, rtm_q_prof, rtm_tpw_prof, rtm_o3_prof
        ! (45) channel dimension should be pretty sparse
        ! only virtual memory if we don't touch the missing channels
        ! dimensions (nlevels, rtm_numel_max, channels)
        real*4, dimension(:,:,:),allocatable :: rtm_arad_prof, rtm_atran_prof
        integer*1, dimension(rtm_num_channels_max) :: rtm_channel_allocated


      interface
        integer(c_int) function sibling_clone() bind(C, name="sibling_clone")
          use iso_c_binding, only: c_int
        end function sibling_clone
      end interface

      interface
        subroutine reopen_files(parent_pid) bind(C, name="reopen_files")
          use iso_c_binding, only: c_int
          integer(c_int), VALUE :: parent_pid
        end subroutine reopen_files
      end interface

      interface
        integer(c_int) function getpid_nocache() bind(C, name="getpid_nocache")
          use iso_c_binding, only: c_int
        end function getpid_nocache
      end interface

      interface
        integer(c_int) function mmap_heap() bind(C, name="mmap_heap")
          use iso_c_binding, only: c_int
        end function mmap_heap
      end interface


      contains


      subroutine calculate_rtm_shapes(num_channels)
        integer*8,intent(out) :: num_channels
        integer*8 i,j, c, z

        rtm_channel_allocated(:) = 0

        rtm_shape_3d(1) = NLevels_Rtm
        rtm_shape_3d(2) = size(RTM, 1)
        rtm_shape_3d(3) = size(RTM, 2)

        rtm_numel = 0
        do i=1, rtm_shape_3d(2)
          do j=1, rtm_shape_3d(3)
            if(RTM(i,j)%is_set) then
              do z=1, RTM_NVZEN
                if(RTM(i,j)%d(z)%is_set) then
                  rtm_numel = rtm_numel + 1
                  do c=1,rtm_num_channels_max
                    if(RTM(i,j)%d(z)%ch(c)%is_set .and. RTM(i,j)%d(z)%ch(c)%is_allocated &
                        .and. allocated(RTM(i,j)%d(z)%ch(c)%Rad_Atm_Profile)) then
                      rtm_channel_allocated(c) = 1
                    endif
                  enddo
                endif
              enddo
            endif
          enddo
        enddo

      num_channels = sum(rtm_channel_allocated)

      end subroutine

      subroutine rtm_realloc(num_channels)
        integer*8, intent(in) :: num_channels

        if(allocated(rtm_t_prof)) deallocate(rtm_t_prof) 
        if(allocated(rtm_z_prof)) deallocate(rtm_z_prof) 
        if(allocated(rtm_q_prof)) deallocate(rtm_q_prof) 
        if(allocated(rtm_tpw_prof)) deallocate(rtm_tpw_prof) 
        if(allocated(rtm_o3_prof)) deallocate(rtm_o3_prof) 
        if(allocated(rtm_index_0)) deallocate(rtm_index_0) 
        if(allocated(rtm_index_1)) deallocate(rtm_index_1) 
        if(allocated(rtm_zen_index)) deallocate(rtm_zen_index) 
        if(allocated(rtm_sfc_index)) deallocate(rtm_sfc_index) 
        if(allocated(rtm_atran_prof)) deallocate(rtm_atran_prof) 
        if(allocated(rtm_arad_prof)) deallocate(rtm_arad_prof) 

        allocate(rtm_t_prof(rtm_shape_3d(1), rtm_numel))
        allocate(rtm_z_prof(rtm_shape_3d(1), rtm_numel))
        allocate(rtm_q_prof(rtm_shape_3d(1), rtm_numel))
        allocate(rtm_tpw_prof(rtm_shape_3d(1), rtm_numel))
        allocate(rtm_o3_prof(rtm_shape_3d(1), rtm_numel))
        allocate(rtm_index_0(rtm_numel))
        allocate(rtm_index_1(rtm_numel))
        allocate(rtm_zen_index(rtm_numel))
        allocate(rtm_sfc_index(rtm_numel))
        allocate(rtm_atran_prof(rtm_shape_3d(1), rtm_numel, num_channels))
        allocate(rtm_arad_prof(rtm_shape_3d(1), rtm_numel, num_channels))
      end subroutine

      subroutine load_rtm_symbols()
        integer*8 n, i, j, element, num_channels, ichan, c, z

        rtm_numel = 0

        n = 1
        if (.not. allocated(RTM)) then
            return
        endif

        call calculate_rtm_shapes(num_channels)
        call rtm_realloc(num_channels)

        element = 0
        do j=1, rtm_shape_3d(3)
          do i=1, rtm_shape_3d(2)
            if(RTM(i,j)%is_allocated) then
              do z=1, RTM_NVZEN
                if(RTM(i,j)%d(z)%is_allocated) then
                  element = element + 1
                  rtm_t_prof(:,element) = RTM(i,j)%T_Prof(:)
                  rtm_z_prof(:,element) = RTM(i,j)%Z_Prof(:)
                  rtm_q_prof(:,element) = RTM(i,j)%Wvmr_Prof(:)
                  rtm_tpw_prof(:,element) = RTM(i,j)%Tpw_Prof(:)
                  rtm_o3_prof(:,element) = RTM(i,j)%Ozmr_Prof(:)
                  rtm_index_0(element) = i
                  rtm_index_1(element) = j
                  rtm_zen_index(element) = z
                  rtm_sfc_index(element) = RTM(i,j)%Sfc_Level
                  ichan = 0
                  do c=1, rtm_num_channels_max
                    if(RTM(i,j)%d(z)%ch(c)%is_allocated) then
                      if(rtm_channel_allocated(c) == 1) then
                        ichan = ichan + 1
                        if(RTM(i,j)%d(z)%ch(c)%is_set .and. RTM(i,j)%d(z)%ch(c)%is_allocated &
                            .and. allocated(RTM(i,j)%d(z)%ch(c)%Rad_Atm_Profile) &
                            .and. size(RTM(i,j)%d(z)%ch(c)%Rad_Atm_Profile) == rtm_shape_3d(1)) then
                            rtm_arad_prof(:,element,ichan) = RTM(i,j)%d(z)%ch(c)%Rad_Atm_Profile(:)
                        else
                            rtm_arad_prof(:,element,ichan) = -9999.
                        endif

                        if(RTM(i,j)%d(z)%ch(c)%is_set .and. RTM(i,j)%d(z)%ch(c)%is_allocated &
                            .and. allocated(RTM(i,j)%d(z)%ch(c)%Trans_Atm_Profile) &
                            .and. size(RTM(i,j)%d(z)%ch(c)%Trans_Atm_Profile) == rtm_shape_3d(1)) then
                            rtm_atran_prof(:,element,ichan) = RTM(i,j)%d(z)%ch(c)%Trans_Atm_Profile(:)
                        else
                            rtm_atran_prof(:,element,ichan) = -9999.
                        endif

                      endif
                    endif
                  enddo
                endif
              enddo
            endif
          enddo
        enddo

        call add_sym_r4_2d(rtm_t_prof, 'rtm_t_profile')
        call add_sym_r4_2d(rtm_z_prof, 'rtm_z_profile')
        call add_sym_r4_2d(rtm_q_prof, 'rtm_q_profile')
        call add_sym_r4_2d(rtm_tpw_prof, 'rtm_tpw_profile')
        call add_sym_r4_2d(rtm_o3_prof, 'rtm_o3_profile')
        call add_sym_i8_1d(rtm_index_0, 'rtm_index_0')
        call add_sym_i8_1d(rtm_index_1, 'rtm_index_1')
        call add_sym_i8_1d(rtm_zen_index, 'rtm_zen_index')
        call add_sym_i8_1d(rtm_sfc_index, 'rtm_sfc_index')
        call add_sym_r4_3d(rtm_atran_prof, 'rtm_atran_profile')
        call add_sym_r4_3d(rtm_arad_prof, 'rtm_arad_profile')
        ! Not allocatable
        num_symbols_i1_1d = num_symbols_i1_1d + 1;
        symbol_ptrs_i1_1d(num_symbols_i1_1d) = loc(rtm_channel_allocated);
        symbol_names_i1_1d(num_symbols_i1_1d) = 'rtm_channel_allocated';
        symbol_shapes_i1_1d(1,num_symbols_i1_1d) = size(rtm_channel_allocated)

        num_symbols_i8_1d = num_symbols_i8_1d + 1;
        symbol_ptrs_i8_1d(num_symbols_i8_1d) = loc(rtm_shape_3d);
        symbol_names_i8_1d(num_symbols_i8_1d) = 'rtm_shape_3d';
        symbol_shapes_i8_1d(1,num_symbols_i8_1d) = size(rtm_shape_3d)

      end subroutine

      subroutine add_sym_r4_3d(sym, nam)
        real*4, dimension(:,:,:), allocatable, intent(in) :: sym
        character(LEN=*), intent(in) :: nam

        if(allocated(sym)) then;
            num_symbols_r4_3d = num_symbols_r4_3d + 1;
            symbol_ptrs_r4_3d(num_symbols_r4_3d) = loc(sym);
            symbol_names_r4_3d(num_symbols_r4_3d) = nam;
            ! I want row-major shapes
            symbol_shapes_r4_3d(1,num_symbols_r4_3d) = size(sym, 3)
            symbol_shapes_r4_3d(2,num_symbols_r4_3d) = size(sym, 2)
            symbol_shapes_r4_3d(3,num_symbols_r4_3d) = size(sym, 1)
        endif
      end subroutine

      subroutine add_sym_r4_2d(sym, nam)
        real*4, dimension(:,:), allocatable, intent(in) :: sym
        character(LEN=*), intent(in) :: nam

        if(allocated(sym)) then;
            num_symbols_r4_2d = num_symbols_r4_2d + 1;
            symbol_ptrs_r4_2d(num_symbols_r4_2d) = loc(sym);
            symbol_names_r4_2d(num_symbols_r4_2d) = nam;
            ! I want row-major shapes
            symbol_shapes_r4_2d(1,num_symbols_r4_2d) = size(sym, 2)
            symbol_shapes_r4_2d(2,num_symbols_r4_2d) = size(sym, 1)
        endif
      end subroutine

      subroutine add_sym_i1_2d(sym, nam)
        integer*1, dimension(:,:), allocatable, intent(in) :: sym
        character(LEN=*), intent(in) :: nam

        if(allocated(sym)) then;
            num_symbols_i1_2d = num_symbols_i1_2d + 1;
            symbol_ptrs_i1_2d(num_symbols_i1_2d) = loc(sym);
            symbol_names_i1_2d(num_symbols_i1_2d) = nam;
            ! I want row-major shapes
            symbol_shapes_i1_2d(1,num_symbols_i1_2d) = size(sym, 2)
            symbol_shapes_i1_2d(2,num_symbols_i1_2d) = size(sym, 1)
        endif
      end subroutine

      subroutine add_sym_i4_2d(sym, nam)
        integer*4, dimension(:,:), allocatable, intent(in) :: sym
        character(LEN=*), intent(in) :: nam

        if(allocated(sym)) then;
            num_symbols_i4_2d = num_symbols_i4_2d + 1;
            symbol_ptrs_i4_2d(num_symbols_i4_2d) = loc(sym);
            symbol_names_i4_2d(num_symbols_i4_2d) = nam;
            ! I want row-major shapes
            symbol_shapes_i4_2d(1,num_symbols_i4_2d) = size(sym, 2)
            symbol_shapes_i4_2d(2,num_symbols_i4_2d) = size(sym, 1)
        endif
      end subroutine

      subroutine add_sym_i8_2d(sym, nam)
        integer*8, dimension(:,:), allocatable, intent(in) :: sym
        character(LEN=*), intent(in) :: nam

        if(allocated(sym)) then;
            num_symbols_i8_2d = num_symbols_i8_2d + 1;
            symbol_ptrs_i8_2d(num_symbols_i8_2d) = loc(sym);
            symbol_names_i8_2d(num_symbols_i8_2d) = nam;
            ! I want row-major shapes
            symbol_shapes_i8_2d(1,num_symbols_i8_2d) = size(sym, 2)
            symbol_shapes_i8_2d(2,num_symbols_i8_2d) = size(sym, 1)
        endif
      end subroutine

      subroutine add_sym_i1_3d(sym, nam)
        integer*1, dimension(:,:,:), allocatable, intent(in) :: sym
        character(LEN=*), intent(in) :: nam

        if(allocated(sym)) then;
            num_symbols_i1_3d = num_symbols_i1_3d + 1;
            symbol_ptrs_i1_3d(num_symbols_i1_3d) = loc(sym);
            symbol_names_i1_3d(num_symbols_i1_3d) = nam;
            ! I want row-major shapes
            symbol_shapes_i1_3d(1,num_symbols_i1_3d) = size(sym, 3)
            symbol_shapes_i1_3d(2,num_symbols_i1_3d) = size(sym, 2)
            symbol_shapes_i1_3d(3,num_symbols_i1_3d) = size(sym, 1)
        endif
      end subroutine

      subroutine add_sym_r4_1d(sym, nam)
        real*4, dimension(:), allocatable, intent(in) :: sym
        character(LEN=*), intent(in) :: nam

        if(allocated(sym)) then;
            num_symbols_r4_1d = num_symbols_r4_1d + 1;
            symbol_ptrs_r4_1d(num_symbols_r4_1d) = loc(sym);
            symbol_names_r4_1d(num_symbols_r4_1d) = nam;
            ! I want row-major shapes
            symbol_shapes_r4_1d(1,num_symbols_r4_1d) = size(sym)
        endif
      end subroutine

      subroutine add_sym_i4_1d(sym, nam)
        integer*4, dimension(:), allocatable, intent(in) :: sym
        character(LEN=*), intent(in) :: nam

        if(allocated(sym)) then;
            num_symbols_i4_1d = num_symbols_i4_1d + 1;
            symbol_ptrs_i4_1d(num_symbols_i4_1d) = loc(sym);
            symbol_names_i4_1d(num_symbols_i4_1d) = nam;
            ! I want row-major shapes
            symbol_shapes_i4_1d(1,num_symbols_i4_1d) = size(sym)
        endif
      end subroutine

      subroutine add_sym_i8_1d(sym, nam)
        integer*8, dimension(:), allocatable, intent(in) :: sym
        character(LEN=*), intent(in) :: nam

        if(allocated(sym)) then;
            num_symbols_i8_1d = num_symbols_i8_1d + 1;
            symbol_ptrs_i8_1d(num_symbols_i8_1d) = loc(sym);
            symbol_names_i8_1d(num_symbols_i8_1d) = nam;
            ! I want row-major shapes
            symbol_shapes_i8_1d(1,num_symbols_i8_1d) = size(sym)
        endif
      end subroutine

      subroutine add_sym_r4_0d(sym, nam)
        real*4, intent(in) :: sym
        character(LEN=*), intent(in) :: nam

        num_symbols_r4_0d = num_symbols_r4_0d + 1;
        symbol_ptrs_r4_0d(num_symbols_r4_1d) = loc(sym);
        symbol_names_r4_0d(num_symbols_r4_1d) = nam;
      end subroutine


      subroutine load_symbols()
          integer i;
          integer j;
          integer c;
          character(LEN=symbol_length_max) varname

          num_symbols_r4_3d = 0
          symbol_shapes_r4_3d_ptr = loc(symbol_shapes_r4_3d)
          symbol_ptrs_r4_3d_ptr = loc(symbol_ptrs_r4_3d)
          symbol_names_r4_3d_ptr = loc(symbol_names_r4_3d)

          num_symbols_r4_2d = 0
          symbol_shapes_r4_2d_ptr = loc(symbol_shapes_r4_2d)
          symbol_ptrs_r4_2d_ptr = loc(symbol_ptrs_r4_2d)
          symbol_names_r4_2d_ptr = loc(symbol_names_r4_2d)

          num_symbols_r4_1d = 0
          symbol_shapes_r4_1d_ptr = loc(symbol_shapes_r4_1d)
          symbol_ptrs_r4_1d_ptr = loc(symbol_ptrs_r4_1d)
          symbol_names_r4_1d_ptr = loc(symbol_names_r4_1d)

          num_symbols_r4_0d = 0
          symbol_ptrs_r4_0d_ptr = loc(symbol_ptrs_r4_0d)
          symbol_names_r4_0d_ptr = loc(symbol_names_r4_0d)

          num_symbols_i4_2d = 0
          symbol_shapes_i4_2d_ptr = loc(symbol_shapes_i4_2d)
          symbol_ptrs_i4_2d_ptr = loc(symbol_ptrs_i4_2d)
          symbol_names_i4_2d_ptr = loc(symbol_names_i4_2d)

          num_symbols_i8_2d = 0
          symbol_shapes_i8_2d_ptr = loc(symbol_shapes_i8_2d)
          symbol_ptrs_i8_2d_ptr = loc(symbol_ptrs_i8_2d)
          symbol_names_i8_2d_ptr = loc(symbol_names_i8_2d)

          num_symbols_i4_1d = 0
          symbol_shapes_i4_1d_ptr = loc(symbol_shapes_i4_1d)
          symbol_ptrs_i4_1d_ptr = loc(symbol_ptrs_i4_1d)
          symbol_names_i4_1d_ptr = loc(symbol_names_i4_1d)

          num_symbols_i8_1d = 0
          symbol_shapes_i8_1d_ptr = loc(symbol_shapes_i8_1d)
          symbol_ptrs_i8_1d_ptr = loc(symbol_ptrs_i8_1d)
          symbol_names_i8_1d_ptr = loc(symbol_names_i8_1d)

          num_symbols_i1_1d = 0
          symbol_shapes_i1_1d_ptr = loc(symbol_shapes_i1_1d)
          symbol_ptrs_i1_1d_ptr = loc(symbol_ptrs_i1_1d)
          symbol_names_i1_1d_ptr = loc(symbol_names_i1_1d)

          num_symbols_i1_3d = 0
          symbol_shapes_i1_3d_ptr = loc(symbol_shapes_i1_3d)
          symbol_ptrs_i1_3d_ptr = loc(symbol_ptrs_i1_3d)
          symbol_names_i1_3d_ptr = loc(symbol_names_i1_3d)

          num_symbols_i1_2d = 0
          symbol_shapes_i1_2d_ptr = loc(symbol_shapes_i1_2d)
          symbol_ptrs_i1_2d_ptr = loc(symbol_ptrs_i1_2d)
          symbol_names_i1_2d_ptr = loc(symbol_names_i1_2d)

          number_of_lines = Image%Number_of_Lines

          ! RTM data structure needs special treatment
          call load_rtm_symbols()
          

        ! Non-allocatable must be added manually
        ! Planck_A1
        num_symbols_r4_1d = num_symbols_r4_1d + 1;
        symbol_ptrs_r4_1d(num_symbols_r4_1d) = loc(Planck_A1);
        symbol_names_r4_1d(num_symbols_r4_1d) = 'planck_a1';
        ! I want row-major shapes
        symbol_shapes_r4_1d(1,num_symbols_r4_1d) = size(Planck_A1)

        ! Planck_A2
        num_symbols_r4_1d = num_symbols_r4_1d + 1;
        symbol_ptrs_r4_1d(num_symbols_r4_1d) = loc(Planck_A2);
        symbol_names_r4_1d(num_symbols_r4_1d) = 'planck_a2';
        ! I want row-major shapes
        symbol_shapes_r4_1d(1,num_symbols_r4_1d) = size(Planck_A2)

        ! Planck_Nu
        num_symbols_r4_1d = num_symbols_r4_1d + 1;
        symbol_ptrs_r4_1d(num_symbols_r4_1d) = loc(Planck_Nu);
        symbol_names_r4_1d(num_symbols_r4_1d) = 'planck_nu';
        ! I want row-major shapes
        symbol_shapes_r4_1d(1,num_symbols_r4_1d) = size(Planck_Nu)

        ! sun_earth_distance
        num_symbols_r4_0d = num_symbols_r4_0d + 1;
        symbol_ptrs_r4_0d(num_symbols_r4_0d) = loc(Sun_Earth_Distance);
        symbol_names_r4_0d(num_symbols_r4_0d) = 'sun_earth_distance';

        call add_sym_i4_1d(Image%Scan_Number,'scan_line_number')

        call add_sym_i1_2d(CLDMASK%Cld_Mask,'cldmask')
        call add_sym_r4_2d(CLDMASK%Posterior_Cld_Probability,'cloud_probability')
        call add_sym_r4_2d(CLDMASK%Posterior_Ice_Probability,'ice_cloud_probability')
        call add_sym_r4_2d(CLDMASK%Posterior_Water_Probability,'water_cloud_probability')
        call add_sym_i1_2d(mask_lrc,'mask_lrc')
        call add_sym_i4_2d(i_lrc,'i_lrc')
        call add_sym_i4_2d(j_lrc,'j_lrc')
        call add_sym_i1_2d(Cld_Type,'cloud_type')
        call add_sym_i1_2d(Cld_Phase,'cloud_phase')
        call add_sym_r4_2d(Nav%Lat, 'latitude')
        call add_sym_r4_2d(Nav%Lon, 'longitude')
        call add_sym_r4_2d(Geo%Solzen,'solar_zenith_angle')
        call add_sym_r4_2d(Geo%Satzen,'sensor_zenith_angle')
        call add_sym_r4_2d(Geo%Relaz,'relative_azimuth_angle')
        call add_sym_r4_2d(Geo%Glintzen,'glint_zenith')
        call add_sym_i1_2d(Bad_Pixel_Mask, 'bad_pixel_mask')
        call add_sym_i1_2d(Sfc%Land_Mask,'land_mask')
        call add_sym_i1_2d(Sfc%Snow, 'snow_class')
        call add_sym_i1_2d(Sfc%Sfc_Type, 'surface_type')
        call add_sym_i1_2d(Sfc%Land,'land_class')
        call add_sym_i1_2d(CLDMASK%Bayes_Mask_Sfc_Type,'bayes_mask_sfc_type')
        call add_sym_i1_3d(CLDMASK%Cld_Test_Vector_Packed,'cloud_mask_test_packed_results')
        call add_sym_r4_2d(Ch(1)%Ref_Toa,'refl_0_65um')
        call add_sym_r4_2d(Ch(2)%Ref_Toa,'refl_0_85um')
        call add_sym_r4_2d(Ch(3)%Ref_Toa,'refl_0_47um')
        call add_sym_r4_2d(Ch(4)%Ref_Toa,'refl_0_55um')
        call add_sym_r4_2d(Ch(5)%Ref_Toa,'refl_1_24um')
        call add_sym_r4_2d(Ch(6)%Ref_Toa,'refl_1_61um')
        call add_sym_r4_2d(Ch(7)%Ref_Toa,'refl_2_25um')
        call add_sym_r4_2d(Ch(8)%Ref_Toa,'refl_0_41um')
        call add_sym_r4_2d(Ch(9)%Ref_Toa,'refl_0_45um')
        call add_sym_r4_2d(Ch(15)%Ref_Toa,'refl_0_75um')
        call add_sym_r4_2d(Ch(20)%Ref_Toa,'refl_3_75um')
        call add_sym_r4_2d(Ch(20)%Bt_Toa,'temp_3_75um')
        call add_sym_r4_2d(Ch(22)%Bt_Toa,'temp_4_05um')
        call add_sym_r4_2d(Ch(26)%Ref_Toa,'refl_1_37um')
        call add_sym_r4_2d(Ch(29)%Bt_Toa,'temp_8_50um')
        call add_sym_r4_2d(Ch(31)%Bt_Toa,'temp_11_0um')
        call add_sym_r4_2d(Ch(32)%Bt_Toa,'temp_12_0um')
        call add_sym_r4_2d(Ch(31)%Bt_Toa_Clear,'temp_11_0um_clear_sky')
        call add_sym_r4_2d(Ch(31)%Rad_Toa_Clear,'rad_11_0um_clear_sky')
        call add_sym_r4_2d(Ch(31)%Sfc_Emiss, 'emiss_sfc_11_0um_nom')
        call add_sym_r4_2d(Ch(31)%Bt_Toa_Std_3x3, 'temp_11_0um_nom_stddev_3x3')
        call add_sym_r4_2d(Ch(1)%Ref_Toa_Std_3x3, 'refl_0_65um_nom_stddev_3x3')

        call add_sym_r4_2d(Ch(1)%Sfc_Ref_White_Sky,'refl_sfc_white_sky_0_65um_nom')

        call add_sym_r4_2d(Ch(20)%Rad_Toa,'rad_3_75um')
        call add_sym_r4_2d(Ch(22)%Rad_Toa,'rad_4_05um')
        call add_sym_r4_2d(Ch(29)%Rad_Toa,'rad_8_50um')
        call add_sym_r4_2d(Ch(38)%Rad_Toa,'rad_10_4um')
        call add_sym_r4_2d(Ch(31)%Rad_Toa,'rad_11_0um')
        call add_sym_r4_2d(Ch(32)%Rad_Toa,'rad_12_0um')

        call add_sym_r4_2d(CCL%Cloud_Fraction,'cloud_fraction')

        ! VIIRS
        if(allocated(viirs_out%mband)) then
            call add_sym_r4_2d(viirs_out%mband(1)%ref,'viirs_m1_ref')
            call add_sym_r4_2d(viirs_out%mband(2)%ref,'viirs_m2_ref')
            call add_sym_r4_2d(viirs_out%mband(3)%ref,'viirs_m3_ref')
            call add_sym_r4_2d(viirs_out%mband(4)%ref,'viirs_m4_ref')
            call add_sym_r4_2d(viirs_out%mband(5)%ref,'viirs_m5_ref')
            call add_sym_r4_2d(viirs_out%mband(6)%ref,'viirs_m6_ref')
            call add_sym_r4_2d(viirs_out%mband(7)%ref,'viirs_m7_ref')
            call add_sym_r4_2d(viirs_out%mband(8)%ref,'viirs_m8_ref')
            call add_sym_r4_2d(viirs_out%mband(9)%ref,'viirs_m9_ref')
            call add_sym_r4_2d(viirs_out%mband(10)%ref,'viirs_m10_ref')
            call add_sym_r4_2d(viirs_out%mband(11)%ref,'viirs_m11_ref')
            call add_sym_r4_2d(viirs_out%mband(12)%rad,'viirs_m12_rad')
            call add_sym_r4_2d(viirs_out%mband(13)%rad,'viirs_m13_rad')
            call add_sym_r4_2d(viirs_out%mband(14)%rad,'viirs_m14_rad')
            call add_sym_r4_2d(viirs_out%mband(15)%rad,'viirs_m15_rad')
            call add_sym_r4_2d(viirs_out%mband(16)%rad,'viirs_m16_rad')
        endif
        if(allocated(viirs_out%iband)) then
            call add_sym_r4_2d(viirs_out%iband(1)%ref,'viirs_i1_ref')
            call add_sym_r4_2d(viirs_out%iband(2)%ref,'viirs_i2_ref')
            call add_sym_r4_2d(viirs_out%iband(3)%ref,'viirs_i3_ref')
            call add_sym_r4_2d(viirs_out%iband(4)%bt,'viirs_i4_bt')
            call add_sym_r4_2d(viirs_out%iband(5)%bt,'viirs_i5_bt')
        endif
        call add_sym_r4_2d(nasa_viirs_i1_ref, 'nasa_viirs_i1_ref')
        call add_sym_r4_2d(nasa_viirs_i2_ref, 'nasa_viirs_i2_ref')
        call add_sym_r4_2d(nasa_viirs_i3_ref, 'nasa_viirs_i3_ref')
        call add_sym_r4_2d(nasa_viirs_i4_bt, 'nasa_viirs_i4_bt')
        call add_sym_r4_2d(nasa_viirs_i5_bt, 'nasa_viirs_i5_bt')

        ! AHI
        call add_sym_r4_3d(ABHI_1km, 'ahi_rad_multich_1km')
        call add_sym_r4_3d(ABHI_500m, 'ahi_rad_multich_500m')
        call add_sym_r4_1d(ABHI_Rad_To_Ref_Fac, 'ahi_rad_to_refl_fac')

        ! ABI
        call add_sym_r4_3d(ABHI_1km, 'abi_rad_multich_1km')
        call add_sym_r4_3d(ABHI_500m, 'abi_rad_multich_500m')
        call add_sym_r4_1d(ABHI_Rad_To_Ref_Fac, 'abi_rad_to_refl_fac')


        ! ACHA solves
        call add_sym_r4_2d(ACHA%Tc, 'cld_temp_acha')                      !x1
        call add_sym_r4_2d(ACHA%Ec, 'cld_emis_acha')                      !x2
        call add_sym_r4_2d(ACHA%Beta, 'cld_beta_acha')                    !x3
        call add_sym_r4_2d(ACHA%Ice_Probability, 'ice_prob_acha')    !x5

        call add_sym_r4_2d(ACHA%Pc, 'cld_press_acha')
        call add_sym_r4_2d(Tc_Opaque_Cloud, 'cld_temp_opaque')

        call add_sym_r4_2d(ACHA%Zc, 'cld_height_acha')
        call add_sym_r4_2d(ACHA%Zc_Eff, 'cld_height_eff_acha')
        call add_sym_r4_2d(ACHA%Zc_Base, 'cld_height_base_acha')

        call add_sym_r4_2d(ACHA%Ec_11um, 'acha_ec_11um')
        call add_sym_r4_2d(ACHA%Ec_104um, 'acha_ec_104um')
        call add_sym_r4_2d(ACHA%Ec_12um, 'acha_ec_12um')
        call add_sym_r4_2d(ACHA%Ec_85um, 'acha_ec_85um')

        call add_sym_r4_2d(ACHA%Goodness, 'acha_goodness')
        call add_sym_i1_2d(ACHA%Quality_Flag, 'acha_quality')
        call add_sym_i1_2d(ACHA%Processing_Order, 'acha_proc_order')
        call add_sym_i1_2d(ACHA%Inversion_Flag, 'acha_inversion_flag')
        call add_sym_r4_2d(ACHA%Lower_Tc, 'cld_temp_lower_acha')
        call add_sym_i1_2d(ACHA%Cld_Type,'cloud_type_acha')
        call add_sym_r4_2d(ACHA%Tc_Ap,'cld_temp_prior_acha')
        call add_sym_r4_2d(ACHA%Tc_Uncertainty, 'cld_temp_uncer_acha')
        call add_sym_r4_2d(ACHA%Tc_Ap_Uncer, 'cld_temp_uncer_prior_acha')
        call add_sym_r4_2d(ACHA%Ec_Ap,'cld_emis_prior_acha')
        call add_sym_r4_2d(ACHA%Ec_Uncertainty,'cld_emis_uncer_acha')
        call add_sym_r4_2d(ACHA%Ec_Ap_Uncer,'cld_emis_uncer_prior_acha')

        ! RTM / NWP
        call add_sym_i4_2d(Zen_Idx_Rtm, 'rtm_zenith_idx')
        call add_sym_i4_2d(NWP_PIX%I_Nwp, 'i_nwp')
        call add_sym_i4_2d(NWP_PIX%J_Nwp, 'j_nwp')
        call add_sym_r4_2d(NWP_PIX%Ttropo, 'tropopause_temperature_nwp')
        call add_sym_r4_2d(NWP_PIX%Tsfc, 'surface_temperature_nwp')
        call add_sym_r4_2d(NWP_PIX%Psfc, 'surface_pressure_nwp')
        call add_sym_r4_2d(NWP_PIX%Tpw, 'total_precipitable_water_nwp')
        call add_sym_r4_2d(NWP_PIX%Ozone, 'total_column_ozone_nwp')
        ! DCOMP
        call add_sym_r4_2d(Tau_DCOMP, 'cld_opd_dcomp')
        call add_sym_r4_2d(Tau_DCOMP_1, 'cld_opd_dcomp_1')
        call add_sym_r4_2d(Tau_DCOMP_2, 'cld_opd_dcomp_2')
        call add_sym_r4_2d(Tau_DCOMP_3, 'cld_opd_dcomp_3')
        call add_sym_r4_2d(Reff_DCOMP, 'cld_reff_dcomp')
        call add_sym_r4_2d(Reff_DCOMP_1, 'cld_reff_dcomp_1')
        call add_sym_r4_2d(Reff_DCOMP_2, 'cld_reff_dcomp_2')
        call add_sym_r4_2d(Reff_DCOMP_3, 'cld_reff_dcomp_3')
        call add_sym_r4_2d(Cwp_Dcomp, 'cld_cwp_dcomp')
        call add_sym_i1_2d(DCOMP_Quality_Flag, 'dcomp_quality')
        call add_sym_r4_2d(Insolation_DCOMP, 'insolation_dcomp')
        call add_sym_r4_2d(Insolation_Diffuse_DCOMP, 'insolation_diffuse_dcomp')

        ! May alias the human-readable symbols, but converting between channels
        ! and names outside clavrx is very annoying
        do c=1,NCHAN_CLAVRX
            ! surface emissivity
            write (varname, "(A,I0.2)") "emiss_sfc_ch", c
            call add_sym_r4_2d(Ch(c)%Sfc_Emiss,varname)

            ! surface refl
            write (varname, "(A,I0.2)") "refl_sfc_white_sky_ch", c
            call add_sym_r4_2d(Ch(c)%Sfc_Ref_White_Sky,varname)

            ! reflectance
            write (varname, "(A,I0.2)") "refl_ch", c
            call add_sym_r4_2d(Ch(c)%Ref_Toa,varname)

            ! reflectance stddev
            write (varname, "(A,I0.2)") "refl_stddev3x3_ch", c
            call add_sym_r4_2d(Ch(c)%Ref_Toa_Std_3x3, varname)

            ! radiance
            write (varname, "(A,I0.2)") "rad_ch", c
            call add_sym_r4_2d(Ch(c)%Rad_Toa,varname)

            ! brightness temperature
            write (varname, "(A,I0.2)") "temp_ch", c
            call add_sym_r4_2d(Ch(c)%Bt_Toa,varname)

            ! brightness temp stddev
            write (varname, "(A,I0.2)") "temp_stddev3x3_ch", c
            call add_sym_r4_2d(Ch(c)%Bt_Toa_Std_3x3, varname)

            ! brightness temp max
            write (varname, "(A,I0.2)") "temp_max3x3_ch", c
            call add_sym_r4_2d(Ch(c)%Bt_Toa_Max_3x3, varname)

            ! brightness temp min
            write (varname, "(A,I0.2)") "temp_min3x3_ch", c
            call add_sym_r4_2d(Ch(c)%Bt_Toa_Min_3x3, varname)

            ! clear sky BT
            write (varname, "(A,I0.2)") "temp_clear_sky_ch", c
            call add_sym_r4_2d(Ch(c)%Bt_Toa_Clear, varname)

            ! clear sky refl
            write (varname, "(A,I0.2)") "refl_clear_sky_ch", c
            call add_sym_r4_2d(Ch(c)%Ref_Toa_Clear, varname)

            ! atmospheric corrected refl
            write (varname, "(A,I0.2)") "refl_atmos_corr_ch", c
            call add_sym_r4_2d(Ch(c)%Ref_Sfc, varname)

            ! atmospheric transmission
            write (varname, "(A,I0.2)") "trans_atm_ch", c
            call add_sym_r4_2d(Ch(c)%Trans_Atm, varname)

            ! Emissivity needed with troposphere temp to match obs
            write (varname, "(A,I0.2)") "emiss_tropo_ch", c
            call add_sym_r4_2d(Ch(c)%Emiss_Tropo, varname)
        enddo

      end subroutine


      subroutine maybe_clone()
        integer my_pid, probe_pid, child_pid, i, l, num_clones, stat
        character(len=7):: Pid_String
        character(len=30):: Stdout_Filename, Stderr_Filename
        character(len=30):: num_clones_str
        num_clones = 0
        if(Tracer_Flag == 1) then
            probe_pid = getpid_nocache()
            my_pid = getpid_nocache()
            !i = mmap_heap()
            CALL GET_ENVIRONMENT_VARIABLE("CLAVRX_TRACER_CLONES", num_clones_str)
            read(num_clones_str,*,iostat=stat)  num_clones
            if (stat /= 0) then
                print*, 'Error reading CLAVRX_TRACER_CLONES, not cloning'
                return
            endif
            do i=1,num_clones
              print*, 'Cloning', my_pid
              call flush(6)
              child_pid = sibling_clone()
              if(child_pid == 0) then
                call reopen_files(probe_pid)
                ! child
                my_pid = getpid_nocache()
                ! Make a new tempdir
                l = len(trim(Temporary_Data_Dir))
                write(Pid_String,'(I7.7)' ) my_pid
                Temporary_Data_Dir = Temporary_Data_Dir(1:l-8) // Pid_String // '/'
                print*, 'New tmpdir', trim(Temporary_Data_Dir)
                call system("mkdir "//trim(Temporary_Data_Dir))
                sibling_pids = 0
                ! Redirect stdout
                write(Stdout_Filename,'(A,I7.7)') 'clavrx.stdout.', my_pid
                close(6)
                open(6, FILE=Stdout_Filename, ACTION='WRITE')
                ! Redirect stderr
                write(Stderr_Filename,'(A,I7.7)') 'clavrx.stderr.', my_pid
                close(7)
                open(7, FILE=Stderr_Filename, ACTION='WRITE')
                ! Exit clone loop
                return
              else
                sibling_pids(i) = child_pid
              endif
            enddo
            print*, sibling_pids
            sibling_pids_addr = loc(sibling_pids)
        endif
      end subroutine

      subroutine stop_and_wait()
          integer pid
          integer ret
          integer sigstop
          sigstop = 19
          print*, 'stop_and_wait'
          pid = getpid_nocache() 
          print*, 'SIGSTOP ',pid
          ! flush stdout
          flush(6)
          ret = kill(pid, sigstop)
          if(ret .ne. 0) then
              call PERROR('error sending signal')
          endif 
      end subroutine

      subroutine Set_Tracer_Flag()
          CHARACTER(len=1) :: tracer_on, skip_output_str
          CALL GET_ENVIRONMENT_VARIABLE("CLAVRX_ENABLE_TRACER",      &
             tracer_on)
          if(tracer_on == '1') then
              Tracer_Flag = 1
          else
              Tracer_Flag = 0
          endif
          CALL GET_ENVIRONMENT_VARIABLE("CLAVRX_SKIP_OUTPUT",      &
             skip_output_str)
          if(skip_output_str == '1') then
              Skip_Output = 1
          else
              Skip_Output = 0
          endif
      end subroutine

      subroutine waitpoint(num)
          integer*4,intent(in) :: num
          if (Tracer_Flag == 1) then
              wait_number = num
              call load_symbols()
              print*, 'waitpoint', wait_number
              call stop_and_wait()
              wait_number = 0
          endif
      end subroutine

      subroutine update_skip_processing(Skip_Processing_Flag)
          integer*4,intent(inout) :: Skip_Processing_Flag
          ! Check heap-allocated value, which may have been set during wait
          if(skip_processing == 1) then
            Skip_Processing_Flag = SYM%YES
          else
            Skip_Processing_Flag = SYM%NO
          endif
          ! Reset to default
          skip_processing = -1
      end subroutine

      endmodule tracer
