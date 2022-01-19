module cx_rttov_bridge_mod

 !use cx_pfaast_constants_mod, only: tstd, ostd, wstd,pstd
    ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit
  
  use cx_rttov_mapping_mod, only: channel_map

  use PIXEL_COMMON_MOD, only: Geo

  use CONSTANTS_MOD, only: Sym
  

  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"




  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=11)  :: NameOfRoutine = 'test_rttovv'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename,cld_coef_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  REAL(KIND=jprb)    :: trans_out(10)
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch
  INTEGER(KIND=jpim) :: np, nch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios
  
  
  
  
  
contains




subroutine compute_transmission_rttov ( &
       ancil_data_path &
       & , pstd &
       & ,temp &
       & ,wvmr &
       & ,ozmr & 
       & ,theta  &
       & ,sensor &
       & ,kban_in &
       & ,taut &
       & , use_modis_channel_equivalent )
       
       
       
  implicit none    
     
  character(len = * ) , intent(in) :: ancil_data_path  
  real, intent(in)  :: pstd (:,:)
  real, intent(in)  :: temp (:,:)
  real, intent(in)  :: wvmr(:,:)
  real, intent(in)  :: ozmr(:,:)
  real, intent(in)  :: theta (:)
  character (len =* ), intent(in) :: sensor
  integer, intent(in)  :: kban_in 
  logical , optional , intent(in) :: use_modis_channel_equivalent 
  real, intent(out)  :: taut (:,:)  
  integer :: lll
   real, parameter :: Q_MIXRATIO_TO_PPMV = 1.60771704e+6
   
   
  opts % rt_ir % addsolar = .FALSE.
         
  opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method
  opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects
  opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects

  if (ANY( .NOT. Geo%Space_Mask) ) then
    opts % rt_ir % ozone_data          = .TRUE. ! Set the relevant flag to .TRUE.
  else
    opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
  endif
  opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
  opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
  opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
  opts % rt_ir % co_data             = .FALSE. !
  opts % rt_ir % so2_data            = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !
  opts % rt_mw % clw_scheme          = 0
  opts % config % verbose            = .FALSE.  ! Enable printing of warnings

 
   lll  = channel_map (sensor,ancil_data_path, kban_in , coef_filename,cld_coef_filename)
   
  
  nprof = size(temp(1,:))
  nlevels = size(temp(:,1))
  
  nthreads = 1
  nchannels = 1
  
  allocate(channel_list(nchannels))
  
  if ( lll .lt. 0 ) then
    if ( allocated(channel_list)) deallocate(channel_list)
    return
  end if
  channel_list = [lll]
  
  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  CALL rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
    nchannels = coefs % coef % fmv_chn
  ENDIF

  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF

  
   ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof
 
  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance, &
        init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF
  
  
  
  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nchannels
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = channel_list(jch)
      
    ENDDO
  ENDDO
  
  deallocate ( channel_list)
  
  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  

 
  profiles(:) % gas_units = 2  ! 
  
  ! Gas units (must be same for all profiles)
! 0 => ppmv over dry air
! 1 => kg/kg over moist air
! 2 => ppmv over moist air
!
  

  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof
    
    
    profiles(iprof) % t(:) = temp(:,iprof)
    profiles(iprof) % p(:) = pstd(:,iprof)
 
    profiles(iprof) % o3(:) = max(ozmr(:,iprof),0.1001E-10)
   
    profiles(iprof) % q(:) = wvmr(:, iprof) * Q_MIXRATIO_TO_PPMV /1000. 
    profiles(iprof) % s2m % p =pstd(96,iprof)
    profiles(iprof) % s2m % t = temp(96,iprof)
    profiles(iprof) % s2m % q = wvmr(96,iprof)  * Q_MIXRATIO_TO_PPMV /1000.
    profiles(iprof) % s2m % u = 0.
    profiles(iprof) % s2m % v = 0.
    profiles(iprof) % s2m % wfetc = 1000000.
    profiles(iprof) % skin % t = temp(96,iprof) + 5.
    profiles(iprof) % zenangle = min(theta(iprof),85.29)
    profiles(iprof) % azangle = 0.
    profiles(iprof) % latitude = 45.
    profiles(iprof) % longitude = 19.
  
  ENDDO
  CLOSE(iup)
 
  
  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities
  emissivity(:) % emis_in = 0._jprb

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

  ! In this example we have no values for input reflectances
  reflectance(:) % refl_in = 0._jprb

  ! Calculate BRDF within RTTOV where the input BRDF value is zero or less
  ! (all channels in this case)
  calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)

  ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
  reflectance(:) % refl_cloud_top = 0._jprb
  
 
  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
 
    CALL rttov_direct(                &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance) ! inout input/output BRDFs per channel
 

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF
 

   
  taut = transmission % tau_levels


 
  ! Deallocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


end subroutine compute_transmission_rttov



end module cx_rttov_bridge_mod
