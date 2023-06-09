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
! Thermal Infrared channels: 
! ch(1:19)=0
!ch(20-25)=MODIS 20-25  - rttov modis 1-6
!ch(26)=0
!ch(27:36)=MODIS 27-36  - rttov 7-16
!ch(37) = ABI 8 - rttov 2
!ch(38) = ABI 12 - rttov 7
!ch (39:41) = 0
!ch(42) = VIIRS I4 - rttov 1
!ch(43) = VIIRS I5 - rttov 6
!ch(44-45) = 0  
!----------------------------------------


module cx_rttov_sfc_emiss

   use CONSTANTS_MOD,only: &
    real4,int1, sym &
    , missing_value_real4 &
    , MIXED_OBS_TYPE &
    , THERMAL_OBS_TYPE &
    , nchan_clavrx
   
   
   use PIXEL_COMMON_MOD, only: Image, Ch, Sensor, Sfc, Use_Land_IR_Emiss

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
   public :: INIT_RTTOV_EMISS
   public :: GET_RTTOV_EMISS
   public :: DESTROY_RTTOV_EMISS
 
  integer(kind=jpim) :: nprof  ! Number of profiles
  integer(kind=jpim), parameter :: nchan_modis = 16  ! Number of channels per profile
  integer(kind=jpim), parameter :: nchan_viirs  = 2 ! Number of channels per profile
  integer(kind=jpim), parameter :: nchan_abi = 2 ! Number of channels per profile

  integer(kind=jpim), parameter :: chn_list_modis(nchan_modis) = &
                                    (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
  integer(kind=jpim), parameter :: chn_list_viirs(nchan_viirs) = (/1,6/)
  integer(kind=jpim), parameter :: chn_list_abi (nchan_abi) = (/2,7/)
 
  real(kind=jprb),    allocatable :: emissivity_modis(:)
  real(kind=jprb),    allocatable :: emissivity_viirs(:)
  real(kind=jprb),    allocatable :: emissivity_abi(:)
  
  integer(kind=jpim) :: imonth        ! Month for which to load data
  integer(kind=jpim) :: npixels,nlines
  integer(kind=jpim) :: lo,hi,j,k, i
  integer(kind=jpim) :: err
  integer(kind=jpim) :: angcorr_run, init_run
  integer(kind=jpim) :: dim1, dim2
 
  type(rttov_profile),  allocatable :: profiles(:)
  
  TYPE(rttov_chanprof),    POINTER :: chanprof_modis(:)    => NULL() ! Input channel/profile list
  TYPE(rttov_chanprof),    POINTER :: chanprof_viirs(:)    => NULL() ! Input channel/profile list
  TYPE(rttov_chanprof),    POINTER :: chanprof_abi(:)    => NULL() ! Input channel/profile list
 
  real*4,    allocatable :: emiss(:,:,:)

  type(rttov_emis_atlas_data)       :: atlas
  
  type(rttov_coefs)                 :: coefs_modis
  type(rttov_coefs)                 :: coefs_viirs
  type(rttov_coefs)                 :: coefs_abi  
  type(rttov_options)               :: opts
  
  
  type init_status_type
    logical :: is_set
    integer :: month
  end type init_status_type
  
  type( init_status_type), save :: init_status

CONTAINS 

  subroutine init_rttov_emiss(rttov_path)
   
  character(len=*) :: rttov_path    ! Path to rttov emis atlas data
  character(len=256) :: coef_path     ! Coefficient director name
  character(len=256) :: atlas_path    ! Path to atlas data
  character(len=256) :: coef_filename ! Coefficient filename
  character(len=256) :: coef_modis_filename ! Coefficient filename
  character(len=256) :: coef_viirs_filename ! Coefficient filename
  character(len=256) :: coef_abi_filename ! Coefficient filename

  integer(kind=jpim) :: atlasid, atlas_type
  
  logical(kind=jplm) :: do_angcorr, single_instrument
  
  
  
#include "rttov_read_coefs.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_setup_emis_atlas.interface"

!-----------------------------------------------------------------------------------
! set input 
!-----------------------------------------------------------------------------------
   atlas_path=trim(rttov_path)//'emis_data/'
   coef_path=trim(rttov_path)//'rtcoef_rttov12/rttov7pred54L/'
   coef_modis_filename=trim(coef_path)//'rtcoef_eos_1_modis.dat'
   coef_viirs_filename=trim(coef_path)//'rtcoef_jpss_0_viirs.dat'
   coef_abi_filename=trim(coef_path)//'rtcoef_goes_16_abi.dat'
   
  
     
    atlas_type = atlas_type_ir ! IR atlas

!  atlas_id
   Use_Land_IR_emiss = 4
 
   if (Use_Land_IR_emiss == 2)  atlasid=1_jpim   ! UWHSREMIS
   if (Use_Land_IR_emiss == 3)  atlasid=2_jpim   ! CAMEL
   if (Use_Land_IR_emiss == 4)  atlasid=3_jpim   ! CAMEL CLIM

  
! define nprof, profiles
!---------------------------------  
    imonth = image % time_start % month
    dim1 =   Image % Number_Of_Elements
    dim2 =   Image % Number_Of_Lines_Per_Segment

    nprof=dim1*dim2 
    
    call rttov_setup_emis_atlas(   &
                    err,               &
                    opts,              &
                    imonth,            &
                    atlas_type_ir,     &  ! IR atlas
                    atlas,             &
                    atlas_id = atlasid,&  ! To select CAMEL climatology atlas
                    path = atlas_path)
      IF (err /= errorstatus_success) THEN
          WRITE(*,*) 'error initialising emissivity atlas'
          CALL rttov_exit(err)
      ENDIF
 
!---------------------------------     
!set for instruments (modis(16), viirs(2), abi(2))
!---------------------------------     
! Set up RTTOV coefficients for MODIS
!---------------------------------
  call rttov_read_coefs(                &
              err,                      &
              coefs_modis,              &
              opts,                     &
              chn_list_modis,&
              file_coef = coef_modis_filename)
     IF (err/= errorstatus_success) THEN
             WRITE(*,*) 'fatal error reading MODIS coefficients'
             CALL rttov_exit(err)
     ENDIF

!---------------------------------     
! Set up RTTOV coefficients for viirs
!---------------------------------
  call rttov_read_coefs(                &
              err,                      &
              coefs_viirs,              &
              opts,                     &
              chn_list_viirs,&
              file_coef = coef_viirs_filename)
    IF (err/= errorstatus_success) THEN
             WRITE(*,*) 'fatal error reading VIIRS coefficients'
             CALL rttov_exit(err)
     ENDIF

!---------------------------------     
! Set up RTTOV coefficients for viirs
!---------------------------------
  call rttov_read_coefs(                &
              err,                      &
              coefs_abi,              &
              opts,                     &
              chn_list_abi,&
              file_coef = coef_abi_filename)
     IF (err/= errorstatus_success) THEN
             WRITE(*,*) 'fatal error reading ABI coefficients'
             CALL rttov_exit(err)
     ENDIF

!---------------------------------
! Set up chanprof for modis viirs abi
!---------------------------------

      allocate(chanprof_modis(nchan_modis*nprof))
      allocate(chanprof_viirs(nchan_viirs*nprof))
      allocate(chanprof_abi(nchan_abi*nprof))

 !MODIS
  do k = 1, nprof
    lo = (k-1)*nchan_modis+1
    hi = lo+nchan_modis-1
    chanprof_modis(lo:hi)%prof = k
    chanprof_modis(lo:hi)%chan = chn_list_modis
  enddo
!VIIRS
  do k = 1, nprof
    lo = (k-1)*nchan_viirs+1
    hi = lo+nchan_viirs-1
    chanprof_viirs(lo:hi)%prof = k
    chanprof_viirs(lo:hi)%chan = chn_list_viirs 
  enddo
!ABI
  do k = 1, nprof
    lo = (k-1)*nchan_abi+1
    hi = lo+nchan_abi-1
    chanprof_abi(lo:hi)%prof = k
    chanprof_abi(lo:hi)%chan = chn_list_abi 
  enddo
  
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
          .or. (init_status % is_set .and. image % time_start % month .ne. init_status % month)) then
          
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
    allocate(emiss1(dim1,dim2,nchan_clavrx))
    allocate(emissivity_modis(nchan_modis*nprof))
    allocate(emissivity_viirs(nchan_viirs*nprof))  
    allocate(emissivity_abi(nchan_abi*nprof))

    
    if (ALL(Space_mask)) then
      emiss1 = missing_value_real4
      deallocate (emiss1)
      deallocate(emissivity_modis)
      deallocate(emissivity_viirs)
      deallocate(emissivity_abi)
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
    ! Retrieve values from atlas for MODIS
    !----------------------------
    call rttov_get_emis(                           &
                  err,                               &
                  opts,                              &
                  chanprof_modis,                    &
                  profiles,                          &
                  coefs_modis,                       &
                  atlas,                             &
                  emissivity_modis(:))
 
    
    if (err /= errorstatus_success) then
              write(*,*) 'error reading modis emissivity atlas'
              call rttov_exit(err)
    end if

    !----------------------------
    ! Retrieve values from atlas for VIIRS
    !----------------------------
    call rttov_get_emis(                           &
                  err,                               &
                  opts,                              &
                  chanprof_viirs,                    &
                  profiles,                          &
                  coefs_viirs,                       &
                  atlas,                             &
                  emissivity_viirs(:))         
    IF (err /= errorstatus_success) THEN
               WRITE(*,*) 'error reading viirs emissivity atlas'
               CALL rttov_exit(err)
    END IF      

    !----------------------------
    ! Retrieve values from atlas for ABI
    !----------------------------
    call rttov_get_emis(                           &
                  err,                               &
                  opts,                              &
                  chanprof_abi,                      &
                  profiles,                          &
                  coefs_abi,                         &
                  atlas,                             &
                  emissivity_abi(:))
    IF (err /= errorstatus_success) THEN
                  WRITE(*,*) 'error reading abi emissivity atlas'
                  CALL rttov_exit(err) 
    ENDIF


       
! creating claverx emissivity 
!--------------------------------    


    k=0
    do Elem_Idx = 1, Image%Number_Of_Elements 
      do Line_Idx = 1, Image%Number_Of_Lines_Per_Segment

        k=k+1
    
        do Chan_Idx = 1, 6
          lo = (k-1)*nchan_modis
          emiss1(Elem_Idx,Line_Idx,Chan_Idx+19)=(emissivity_modis(lo+Chan_Idx))
        end do
        
        do Chan_Idx = 7, nchan_modis
          lo = (k-1)*nchan_modis
          emiss1(Elem_Idx,Line_Idx,Chan_Idx+20)=(emissivity_modis(lo+Chan_Idx))
        end do
        
        do Chan_Idx = 1, nchan_abi
             lo = (k-1)*nchan_abi
             emiss1(Elem_Idx,Line_Idx,Chan_Idx+36)=emissivity_abi(lo+Chan_Idx)
        end do

        do Chan_Idx = 1, nchan_viirs
             lo = (k-1)*nchan_viirs
             emiss1(Elem_Idx,Line_Idx,Chan_Idx+41)=emissivity_viirs(lo+Chan_Idx)
        end do
 
      end do
    end do
  
    do Chan_Idx = 20, Nchan_Clavrx
      if (ch(Chan_Idx)%Obs_Type /= MIXED_OBS_TYPE .and. &
        Ch(Chan_Idx)%Obs_Type /= THERMAL_OBS_TYPE) cycle

      if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
        where(emiss1(:,:,Chan_Idx) .ge. 0)
            Ch(Chan_Idx)%Sfc_Emiss=emiss1(:,:,Chan_Idx)
        end where
      end if
    end do
 
    deallocate (emiss1)
    deallocate(emissivity_modis)
    deallocate(emissivity_viirs)
    deallocate(emissivity_abi)
    deallocate(profiles)
    
  end subroutine get_rttov_emiss


  subroutine DESTROY_RTTOV_EMISS
  
#include "rttov_deallocate_emis_atlas.interface"

!---------------------------------------
  call rttov_deallocate_emis_atlas(atlas)
  call rttov_dealloc_coefs(err, coefs_modis)
  call rttov_dealloc_coefs(err, coefs_viirs)
  call rttov_dealloc_coefs(err, coefs_abi)
    
   if (associated(chanprof_modis)) deallocate(chanprof_modis)
   if (associated(chanprof_viirs))  deallocate(chanprof_viirs)
   if (associated(chanprof_abi))  deallocate(chanprof_abi)
   init_status%is_set = .false.
   init_status%month = -1
 
end subroutine destroy_rttov_emiss
end module cx_rttov_sfc_emiss
