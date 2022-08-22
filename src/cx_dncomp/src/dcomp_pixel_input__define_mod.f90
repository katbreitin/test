! $Id:$
! CREATED Andi Walther {date} 
module dcomp_pixel_input__define_mod
   implicit none
    integer, parameter :: N_CHAN = 50
   type, public :: dcomp_input_structure
      logical :: statusOK
      integer :: n_chan
      logical :: chan_on
      real :: rfl (N_CHAN)
      real :: rfl_unc  (N_CHAN)
      real :: rad (N_CHAN)
      real :: rad_unc (N_CHAN)
      real :: transmission_above_cloud (N_CHAN)
      real :: albedo_surface(N_CHAN)
      real :: zenith_solar
      real :: zenith_sensor
      real :: azimuth_difference
      real :: temp_cloud
      real :: state_apriori(2)
      integer :: phase_cloud
      integer :: snow_class
      logical :: use_snow_retrieval
      real :: radiation_above_cloud(N_CHAN)
      real :: radiation_clear_toc(N_CHAN)
      integer :: mode
      character (len=1024) :: ancil_path
   end type dcomp_input_structure





end module dcomp_pixel_input__define_mod
