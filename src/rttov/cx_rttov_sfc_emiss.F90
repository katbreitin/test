! Description:
!> @file
!!   Executable for basic test of the CAMEL climatology IR emissivity atlas.
!
!> @brief
!!   Executable for basic test of the CAMEL climatology IR emissivity atlas.
!!
!! @details
!!   This should be run using the test_camel_clim_atlas.sh script in rttov_test/
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
!CLAVERX CHANNELS
!



module cx_rttov_sfc_emiss

   use CONSTANTS_MOD,only: &
    real4,int1, sym &
    , missing_value_real4 &
    , MIXED_OBS_TYPE &
    , THERMAL_OBS_TYPE &
    , nchan_clavrx


   use PIXEL_COMMON_MOD, only: Image, Ch, Sensor, Sfc, Use_Land_IR_Emiss &
    , Ancil_Data_Dir


   USE  cx_rttov_mapping_mod

   USE cx_rttov_sensor_mod


#include "throw.h"
   use parkind1, only : jpim, jprb, jplm

   use rttov_unix_env, only : rttov_exit

   use rttov_const, only :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         surftype_land,       &
         surftype_sea,        &
         sensor_id_mw,        &
         sensor_id_po

   use rttov_types, only : &
        rttov_coefs,      &
        rttov_options,    &
        rttov_profile,    &
        rttov_chanprof,   &
        rttov_emissivity

   use mod_rttov_emis_atlas, only : atlas_type_ir, rttov_emis_atlas_data,   &
        uwiremis_atlas_id , camel_clim_atlas_id, camel_atlas_id


   implicit none
   private
   !--- routine access declaration
   private :: INIT_RTTOV_EMISS
   public :: GET_RTTOV_EMISS
   public :: DESTROY_RTTOV_EMISS

   integer(kind=jpim) :: nprof  ! Number of profiles
   integer(kind=jpim)  :: nchan

   integer(kind=jpim) :: chan_list_max(45)
   integer(kind=jpim) :: chan_list_max_cx(45)
   integer(kind=jpim), allocatable :: chan_list(:)
   integer(kind=jpim), allocatable :: chan_list_cx(:)

   real(kind=jprb),    allocatable :: emissivity_sensor(:)

   integer(kind=jpim) :: imonth        ! Month for which to load data
   integer(kind=jpim) :: npixels,nlines
   integer(kind=jpim) :: lo,hi,j,k, i
   integer(kind=jpim) :: err
   integer(kind=jpim) :: angcorr_run, init_run
   integer(kind=jpim) :: dim1, dim2

   type(rttov_profile),  allocatable :: profiles(:)

   TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL()

   real*4,    allocatable :: emiss(:,:,:)

   type(rttov_emis_atlas_data)       :: atlas

   type(rttov_coefs)                 :: coefs
   type(rttov_options)               :: opts

   type init_status_type
     logical :: is_set
     integer :: month
   end type init_status_type

   type( init_status_type), save :: init_status
     type(cx_rttov_sensor_type) :: rt_sensor

contains

   SUBROUTINE init_rttov_emiss(rttov_path)

      character(len=*)   :: rttov_path    ! Path to rttov emis atlas data
      character(len=256) :: atlas_path    ! Path to atlas data
      character(len=256) :: coef_filename ! Coefficient filename
      character(len=256) :: cld_coef_filename ! Coefficient filename

      integer(kind=jpim) :: atlasid, atlas_type
      logical(kind=jplm) :: do_angcorr, single_instrument
      integer :: chn_sensor
      real :: max_satzen
      character(len=50) :: sc_name_rtm
      integer :: chan_idx

#include "rttov_read_coefs.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_setup_emis_atlas.interface"

      !-----------------------------------------------------------------------------------
      ! set input
      !-----------------------------------------------------------------------------------
      atlas_path=trim(rttov_path)//'emis_data/'

      call rt_sensor % init(Sensor%WMO_id,trim(rttov_path))

      call rttov_read_coefs(         &
                err,                &
                coefs,              &
                opts,               &
                file_coef = rt_sensor % coef_filename)

      ! -find channel list we needed for this sensor
      nchan = 0
      do Chan_Idx = 20, Nchan_Clavrx
         if (ch(Chan_Idx)%Obs_Type /= MIXED_OBS_TYPE .and. &
            Ch(Chan_Idx)%Obs_Type /= THERMAL_OBS_TYPE) cycle

         if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
            !print*,rt_sensor % chan_src_on_cx(chan_idx)

             if ( rt_sensor % chan_src_on_cx(chan_idx) .gt. 0) nchan = nchan + 1
             chan_list_max(nchan) = rt_sensor % chan_src_on_cx(chan_idx)
             chan_list_max_cx(nchan) = chan_idx
         end if
      end do

      allocate (chan_list(nchan))
      allocate (chan_list_cx(nchan))
      chan_list = chan_list_max(1:nchan)
      chan_list_cx = chan_list_max_cx(1:nchan)
      chan_list_max = -1
      chan_list_max_cx = -1

      atlas_type = atlas_type_ir ! IR atlas

      ! this is for later optional choice
      Use_Land_IR_emiss = 4

      if (Use_Land_IR_emiss == 2)  atlasid=1_jpim   ! UWHSREMIS
      if (Use_Land_IR_emiss == 3)  atlasid=2_jpim   ! CAMEL
      if (Use_Land_IR_emiss == 4)  atlasid=3_jpim   ! CAMEL CLIM


      ! define nprof, profiles
      !---------------------------------
      imonth = image % time_start % month
      dim1 =   Image % Number_Of_Elements
      dim2 =   Image % Number_Of_Lines_Per_Segment

      nprof = dim1 * dim2

      call rttov_setup_emis_atlas(   &
                    err,               &
                    opts,              &
                    imonth,            &
                    atlas_type_ir,     &  ! IR atlas
                    atlas,             &
                    atlas_id = atlasid,&  ! To select CAMEL climatology atlas
                    path = atlas_path, &
                    coefs = coefs)
      IF (err /= errorstatus_success) THEN
          WRITE(*,*) 'error initialising emissivity atlas'
          CALL rttov_exit(err)
      ENDIF



     !---------------------------------
     ! Set up chanprof
     !---------------------------------

      allocate(chanprof(nchan * nprof))

      do k = 1, nprof
         lo = (k-1)*nchan+1
         hi = lo+nchan-1
         chanprof(lo:hi)%prof = k
         chanprof(lo:hi)%chan = chan_list
      end do

     init_status % is_set  = .true.
     init_status % month = image % time_start % month

  end subroutine init_rttov_emiss

  ! +++++++++++++++++++
  !
  !
  subroutine get_rttov_emiss(lats, lons,space_mask, rttov_path)

    REAL(kind=real4), dimension(:,:), intent(in) :: lats, lons
    logical, dimension(:,:), intent(in) :: space_mask
    character(len=*), intent(in)  :: rttov_path   ! Path to rttov emis atlas data


    integer(kind=jpim), parameter :: ioout = 51 ! Output file unit


    integer:: Elem_Idx, Line_Idx, Chan_Idx

    real(kind=real4),    allocatable :: emiss1(:,:,:)

    logical(kind=jplm) :: do_angcorr, single_instrument

#include "rttov_skipcommentline.interface"
#include "rttov_get_emis.interface"

    !-----------------------------------------------------------------------------------
    ! set input
    !-----------------------------------------------------------------------------------

    ! check if we first have to initialize


    if ((.not. init_status % is_set) &
          .or. (init_status % is_set &
          .and. image % time_start % month .ne. init_status % month)) then

       if  (init_status % is_set)  call  DESTROY_RTTOV_EMISS()

      call INIT_RTTOV_EMISS(trim(rttov_path))
      init_status % is_set = .true.
    end if

    do_angcorr = .FALSE.   ! initialise without angular correction
    single_instrument = .FALSE.    ! initialise for use with a single instrument

    !define nprof, profiles
    !---------------------------------
    dim1 = Image%Number_Of_Elements
    dim2 = Image%Number_Of_Lines_Per_Segment

    nprof=dim1*dim2

    allocate(profiles(nprof))
    allocate(emiss1(dim1,dim2,nchan))
    allocate(emissivity_sensor(nchan*nprof))

    if (ALL(Space_mask)) then
      deallocate (emiss1)
      deallocate(emissivity_sensor)
      deallocate(profiles)
      return
    end if


    k=0
    do Elem_Idx = 1, Image%Number_Of_Elements
      do Line_Idx = 1, Image%Number_Of_Lines_Per_Segment
        k=k+1

        profiles(k)%latitude=dble(lats(Elem_Idx,Line_Idx))

        profiles(k)%longitude=dble(lons(Elem_Idx,Line_Idx))
        profiles(k)%skin%surftype =0  ! land
        profiles(k)%skin%snow_fraction = 0.

        if (Sfc%Snow(Elem_Idx,Line_Idx) == sym%SNOW) profiles(k)%skin%snow_fraction = 1.

      end do
    end do

   ! initialize output matrix
    emiss1=0.

    !
    !----------------------------
    ! Retrieve values from atlas
    !----------------------------
    call rttov_get_emis(                      &
                  err,                        &
                  opts,                       &
                  chanprof,                   &
                  profiles,                   &
                  coefs,                      &
                  atlas,                      &
                  emissivity_sensor(:))

    if (err /= errorstatus_success) then
              write(*,*) 'error reading modis emissivity atlas'
              call rttov_exit(err)
    end if

! creating claverx emissivity
!--------------------------------

    k=0
    do Elem_Idx = 1, Image%Number_Of_Elements
      do Line_Idx = 1, Image%Number_Of_Lines_Per_Segment
        k=k+1
        do Chan_Idx = 1, nchan
          lo = (k-1) * nchan
          emiss1(Elem_Idx,Line_Idx,Chan_Idx)=(emissivity_sensor(lo+Chan_Idx))
        end do
      end do
    end do

    do chan_idx  = 1, nchan
      where(emiss1(:,:,Chan_Idx) .ge. 0)
          Ch(chan_list_cx(chan_idx))%Sfc_Emiss=emiss1(:,:,Chan_Idx)
      end where
    end do

    deallocate (emiss1)
    deallocate ( emissivity_sensor)

    deallocate(profiles)

  end subroutine get_rttov_emiss


  subroutine DESTROY_RTTOV_EMISS

#include "rttov_deallocate_emis_atlas.interface"

!---------------------------------------
  call rttov_deallocate_emis_atlas(atlas)
  call rttov_dealloc_coefs(err, coefs)

   if (associated(chanprof)) deallocate(chanprof)

   init_status%is_set = .false.
   init_status%month = -1

end subroutine destroy_rttov_emiss
end module cx_rttov_sfc_emiss
