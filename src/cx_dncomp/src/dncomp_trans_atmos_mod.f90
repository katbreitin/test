! $Id: dncomp_trans_atmos_mod.f90 182 2018-02-15 13:17:01Z awalther $
!
!

module dncomp_trans_atmos_mod
   implicit none
   private
   real :: gas_coeff(3)
   real :: ozone_coeff(3)
   real :: rayleigh_coeff

   interface trans_atm_above_cloud
      module procedure trans_atm_above_cloud_skalar
      module procedure trans_atm_above_cloud_1d
      module procedure trans_atm_above_cloud_2d
   end interface
   public :: trans_atm_above_cloud



contains
   !
   !
   !
   subroutine trans_atm_above_cloud_skalar ( &
      tpw_ac &
      , ozone_dobson &
      , press_sfc &
      , press_cld &
      , air_mass  &
      , gas_coeff_inp &
      , ozone_coeff_inp &
      , rayleigh_coeff_inp &
      , trans &
      , trans_uncert)

      implicit none
      real, intent(in) :: tpw_ac
      real, intent(in) :: ozone_dobson
      real, intent(in) :: press_sfc
      real, intent(in) :: press_cld
      real, intent(in) :: air_mass
      real, intent(in) :: gas_coeff_inp(3)
      real, intent(in) :: ozone_coeff_inp(3)
      real, intent(in) :: rayleigh_coeff_inp

      real, intent(out) :: trans
      real, intent(out) :: trans_uncert

      call set_coeffs ( gas_coeff_inp, ozone_coeff_inp, rayleigh_coeff_inp)

      call dncomp_trans_atm_above_cloud ( &
         tpw_ac &
         , ozone_dobson &
         , press_sfc &
         , press_cld &
         , air_mass  &
         , trans , trans_uncert )




   end subroutine trans_atm_above_cloud_skalar

   !
   !
   !
   subroutine trans_atm_above_cloud_2d( &
      tpw_ac &
      , ozone_dobson &
      , press_sfc &
      , press_cld &
      , air_mass  &
      , gas_coeff_inp &
      , ozone_coeff_inp &
      , rayleigh_coeff_inp &
      , trans &
      , trans_uncert)


      implicit none
      real, intent(in) :: tpw_ac(:,:)
      real, intent(in) :: ozone_dobson(:,:)
      real, intent(in) :: press_sfc(:,:)
      real, intent(in) :: press_cld(:,:)
      real, intent(in) :: air_mass(:,:)
      real, intent(in) :: gas_coeff_inp(3)
      real, intent(in) :: ozone_coeff_inp(3)
      real, intent(in) :: rayleigh_coeff_inp

      real, intent(out) :: trans(:,:)
      real, intent(out) :: trans_uncert(:,:)

      call set_coeffs ( gas_coeff_inp, ozone_coeff_inp, rayleigh_coeff_inp)

      call dncomp_trans_atm_above_cloud ( &
         tpw_ac &
         , ozone_dobson &
         , press_sfc &
         , press_cld &
         , air_mass  &
         , trans , trans_uncert)

      trans_uncert = 0.02
   end subroutine trans_atm_above_cloud_2d

   !
   !
   !
   subroutine trans_atm_above_cloud_1d( &
      tpw_ac &
      , ozone_dobson &
      , press_sfc &
      , press_cld &
      , air_mass  &
      , gas_coeff_inp &
      , ozone_coeff_inp &
      , rayleigh_coeff_inp &
      , trans &
      , trans_uncert)


      implicit none
      real, intent(in) :: tpw_ac(:)
      real, intent(in) :: ozone_dobson(:)
      real, intent(in) :: press_sfc(:)
      real, intent(in) :: press_cld(:)
      real, intent(in) :: air_mass(:)
      real, intent(in) :: gas_coeff_inp(3)
      real, intent(in) :: ozone_coeff_inp(3)
      real, intent(in) :: rayleigh_coeff_inp

      real, intent(out) :: trans(:)
      real, intent(out) :: trans_uncert(:)

      call set_coeffs ( gas_coeff_inp, ozone_coeff_inp, rayleigh_coeff_inp)

      call dncomp_trans_atm_above_cloud ( &
         tpw_ac &
         , ozone_dobson &
         , press_sfc &
         , press_cld &
         , air_mass  &
         , trans, trans_uncert )

   end subroutine trans_atm_above_cloud_1d

   !
   !
   !
   subroutine set_coeffs ( gas_coeff_inp, ozone_coeff_inp, rayleigh_coeff_inp)
      real, intent(in) :: gas_coeff_inp(3)
      real, intent(in) :: ozone_coeff_inp(3)
      real, intent(in) :: rayleigh_coeff_inp

      gas_coeff = gas_coeff_inp
      ozone_coeff = ozone_coeff_inp
      rayleigh_coeff =    rayleigh_coeff_inp

   end subroutine set_coeffs



   elemental  subroutine dncomp_trans_atm_above_cloud ( &
      tpw_ac &
      , ozone_dobson &
      , press_sfc &
      , press_cld &
      , air_mass  &
      , transmission &
      , trans_unc )


      implicit none
      real, intent(in) :: tpw_ac
      real, intent(in) :: ozone_dobson
      real, intent(in) :: press_sfc
      real, intent(in) :: press_cld
      real, intent(in) :: air_mass

      real, intent(out) :: transmission
      real, intent(out) :: trans_unc
      real :: trans_wvp
      real :: trans_wvp_unc
      real :: trans_ozone
      real :: trans_ozone_unc
      real :: trans_rayleigh

      real,parameter  :: ASSUMED_WVP_ERROR = 1.2

      trans_ozone = exp ( -1. * ( ozone_coeff(1) &
         & + ozone_coeff(2) *  ozone_dobson &
         & + ozone_coeff(3) *  ozone_dobson ** 2))
      trans_ozone = min ( trans_ozone,1.)
      trans_ozone = max ( trans_ozone,0.)

      trans_ozone_unc = trans_ozone -  exp ( -1. * ( ozone_coeff(1) &
         & + ozone_coeff(2) *  (1.1 * ozone_dobson) &
         & + ozone_coeff(3) * (1.1 * ozone_dobson)**2))

      trans_ozone_unc = max(min ( trans_ozone_unc,0.02),0.)

      trans_rayleigh = exp (-air_mass  &
         &    * ( rayleigh_coeff *  (press_cld / press_sfc )) * 0.84)


      trans_wvp  =  exp( - 1. * (gas_coeff(1) &
         & + gas_coeff(2) * tpw_ac  &
         & + gas_coeff(3) * ( tpw_ac ** 2 ) ) )

      trans_wvp = min ( trans_wvp, 1.)

      trans_wvp_unc  = abs(trans_wvp - exp ( -1. * (gas_coeff(1)   &
         & + gas_coeff(2) * (assumed_wvp_error * tpw_ac) &
         & + gas_coeff(3) * ( ( assumed_wvp_error * tpw_ac ) **2 ) ) ) )

      trans_wvp_unc = max(min ( trans_wvp_unc,0.1),0.)

      trans_rayleigh =  1.
      transmission = trans_ozone * trans_rayleigh * trans_wvp
      transmission = min ( transmission , 1.)
      transmission = max ( transmission , 0.)
      trans_unc = trans_ozone_unc + trans_wvp_unc
   end subroutine dncomp_trans_atm_above_cloud

end module dncomp_trans_atmos_mod
