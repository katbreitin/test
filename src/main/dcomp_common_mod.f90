module dcomp_common_mod
  use univ_kind_defs_mod, only: f4,  i2, i1
!  use PIXEL_COMMON_MOD
 !use CONSTANTS_MOD

 use univ_fp_comparison_mod, only: operator(.EQfp.), operator(.NEfp.),  &
      operator(.GEfp.), operator(.LEfp.)

  REAL ( f4), parameter ::  MISSING_F4 = -999.
  INTEGER (i1), parameter :: MISSING_I1 = -128_i1
  INTEGER (i2), parameter :: MISSING_I2 = -32768_i2

  type :: dcomp_definition
    logical :: is_set
    integer:: mode
    integer :: dim1,dim2
    real (kind=f4), dimension(:,:), allocatable, public:: Tau
    real (kind=f4), dimension(:,:), allocatable, public:: Tau_Ap
    real (kind=f4), dimension(:,:), allocatable, public:: vis_Ref_fm
    real (kind=f4), dimension(:,:), allocatable, public:: Reff
    real (kind=f4), dimension(:,:), allocatable, public:: Iwp
    real (kind=f4), dimension(:,:), allocatable, public:: Iwp_Tau
    real (kind=f4), dimension(:,:), allocatable, public:: Lwp
    real (kind=f4), dimension(:,:), allocatable, public:: Refl_Asym
    real (kind=f4), dimension(:,:), allocatable, public:: Cwp
    real (kind=f4), dimension(:,:), allocatable, public:: Cwp_Ice_Layer
    real (kind=f4), dimension(:,:), allocatable, public:: Cwp_Water_Layer
    real (kind=f4), dimension(:,:), allocatable, public:: Cwp_Scwater_Layer
    real (kind=f4), dimension(:,:), allocatable, public:: Iwc
    real (kind=f4), dimension(:,:), allocatable, public:: Lwc
    real (kind=f4), dimension(:,:), allocatable, public:: Rain_Rate
    real (kind=f4), dimension(:,:), allocatable, public:: Hcld
    real (kind=f4), dimension(:,:), allocatable, public:: Cdnc
    real (kind=f4), dimension(:,:), allocatable, public:: Tau_Cost
    real (kind=f4), dimension(:,:), allocatable, public:: Reff_Cost
    integer (kind=i1), dimension(:,:), allocatable, public:: Tau_Qf
    integer (kind=i1), dimension(:,:), allocatable, public:: Reff_Qf
    integer (kind=i1), dimension(:,:), allocatable, public:: Quality_Flag
    integer (kind=i2), dimension(:,:), allocatable, public:: Info_Flag
    real (kind=f4), dimension(:,:), allocatable, public:: Cwp_Fit
    real (kind=f4), dimension(:,:), allocatable, public:: Reff_Fit
    real (kind=f4), dimension(:,:), allocatable, public::  olr
    real (kind=f4), dimension(:,:), allocatable, public:: Cloud_063um_Albedo
    real (kind=f4), dimension(:,:), allocatable, public:: Cloud_063um_Spherical_Albedo
    real (kind=f4), dimension(:,:), allocatable, public:: Cloud_063um_Transmission_View
    real (kind=f4), dimension(:,:), allocatable, public:: Cloud_063um_Transmission_Solar
    real (kind=f4), dimension(:,:), allocatable, public:: Insolation
    real (kind=f4), dimension(:,:), allocatable, public:: Insolation_Diffuse
    integer (kind=i1), dimension(:,:), allocatable, public:: phase_used

  contains
    procedure :: allocate => dcomp_allocate
    procedure :: deallocate => dcomp_deallocate
    procedure :: reset => dcomp_reset
    procedure :: COMPUTE_CWP_PHASE
  end type dcomp_definition


contains

  subroutine dcomp_allocate(self,n1,n2)
    class(dcomp_definition) :: self
    integer, intent(in) :: n1,n2
    self % dim1 = n1
    self % dim2 = n2
    allocate (self%tau(n1,n2))
    allocate (self%tau_ap(n1,n2))
    allocate (self%vis_Ref_fm(n1,n2))
    allocate (self%Reff(n1,n2))
    allocate (self%iwp(n1,n2))
    allocate (self%iwp_tau(n1,n2))
    allocate (self%lwp(n1,n2))
    allocate (self%Refl_Asym(n1,n2))
    allocate (self%cwp(n1,n2))
    allocate (self%Cwp_Ice_Layer(n1,n2))
    allocate (self%Cwp_Water_Layer(n1,n2))
    allocate (self%Cwp_Scwater_Layer(n1,n2))
    allocate (self%iwc(n1,n2))
    allocate (self%lwc(n1,n2))
    allocate (self%rain_rate(n1,n2))
    allocate (self%hcld(n1,n2))
    allocate (self%cdnc(n1,n2))
    allocate (self%tau_cost(n1,n2))
    allocate (self%reff_cost(n1,n2))
    allocate (self%tau_qf(n1,n2))
    allocate (self%reff_qf(n1,n2))
    allocate (self%Quality_Flag(n1,n2))
    allocate (self%Info_Flag(n1,n2))
    allocate (self%cwp_fit(n1,n2))
    allocate (self%reff_fit(n1,n2))
    allocate (self%olr(n1,n2))
    allocate (self%Cloud_063um_Albedo(n1,n2))
    allocate (self%Cloud_063um_Spherical_Albedo(n1,n2))
    allocate (self%Cloud_063um_Transmission_View(n1,n2))
    allocate (self%Cloud_063um_Transmission_Solar(n1,n2))
    allocate (self%Insolation(n1,n2))
    allocate (self%Insolation_Diffuse(n1,n2))
    allocate (self%phase_used(n1,n2))

    self % is_set = .true.


  end subroutine dcomp_allocate

  subroutine dcomp_deallocate(self)
    class(dcomp_definition) :: self

    deallocate (self%tau)
    deallocate (self%tau_ap)
    deallocate (self%vis_Ref_fm)
    deallocate (self%Reff)
    deallocate (self%iwp)
    deallocate (self%iwp_tau)
    deallocate (self%lwp)
    deallocate (self%Refl_Asym)
    deallocate (self%cwp)
    deallocate (self%Cwp_Ice_Layer)
    deallocate (self%Cwp_Water_Layer)
    deallocate (self%Cwp_Scwater_Layer)
    deallocate (self%iwc)
    deallocate (self%lwc)
    deallocate (self%rain_rate)
    deallocate (self%hcld)
    deallocate (self%cdnc)
    deallocate (self%tau_cost)
    deallocate (self%reff_cost)
    deallocate (self%tau_qf)
    deallocate (self%reff_qf)
    deallocate (self%Quality_Flag)
    deallocate (self%Info_Flag)
    deallocate (self%cwp_fit)
    deallocate (self%reff_fit)
    deallocate (self%olr)
    deallocate (self%Cloud_063um_Albedo)
    deallocate (self%Cloud_063um_Spherical_Albedo)
    deallocate (self%Cloud_063um_Transmission_View)
    deallocate (self%Cloud_063um_Transmission_Solar)
    deallocate (self%Insolation)
    deallocate (self%Insolation_Diffuse)
    deallocate (self%phase_used)

    self % is_set = .false.
  end subroutine dcomp_deallocate


  subroutine dcomp_reset(self)
    class(dcomp_definition) :: self
    self%tau = MISSING_F4
    self%tau_ap = MISSING_F4
    self%vis_Ref_fm = MISSING_F4
    self%Reff = MISSING_F4
    self%iwp = MISSING_F4
    self%iwp_tau = MISSING_F4
    self%lwp = MISSING_F4
    self%Refl_Asym = MISSING_F4
    self%cwp = MISSING_F4
    self%Cwp_Ice_Layer = MISSING_F4
    self%Cwp_Water_Layer = MISSING_F4
    self%Cwp_Scwater_Layer = MISSING_F4
    self%iwc = MISSING_F4
    self%lwc = MISSING_F4
    self%rain_rate = MISSING_F4
    self%hcld = MISSING_F4
    self%cdnc = MISSING_F4
    self%tau_cost = MISSING_F4
    self%reff_cost = MISSING_F4
    self%tau_qf = MISSING_F4
    self%reff_qf = MISSING_F4
    self%Quality_Flag = MISSING_F4
    self%Info_Flag = MISSING_F4
    self%cwp_fit = MISSING_F4
    self%reff_fit = MISSING_F4
    self%olr = MISSING_F4
    self%Cloud_063um_Albedo = MISSING_F4
    self%Cloud_063um_Spherical_Albedo = MISSING_F4
    self%Cloud_063um_Transmission_View = MISSING_F4
    self%Cloud_063um_Transmission_Solar = MISSING_F4
    self%Insolation = MISSING_F4
    self%Insolation_Diffuse = MISSING_F4
    self%phase_used  = 0


  end subroutine dcomp_reset

  subroutine COMPUTE_CWP_PHASE(self &
    , Zc_top, Zc_base &
    , Upper_Limit_Water_Height &
    , Freezing_Level_Height)

    class(dcomp_definition) :: self

    real, dimension(:,:), intent(in) :: Zc_top, Zc_base
    real, dimension(:,:), intent(in) :: Upper_Limit_Water_Height
    real, dimension(:,:), intent(in):: Freezing_Level_Height


    real, dimension(:,:), allocatable:: Cloud_Geometrical_Thickness
    real, dimension(:,:), allocatable:: Ice_Layer_Fraction
    real, dimension(:,:), allocatable:: Water_Layer_Fraction
    real, dimension(:,:), allocatable:: Scwater_Layer_Fraction

    ! init
    self % Cwp_Ice_Layer = Missing_Value_Real4
    self % Cwp_Water_Layer = Missing_Value_Real4
    self % Cwp_Scwater_Layer = Missing_Value_Real4

    allocate (Cloud_Geometrical_Thickness(self % dim1, self % dim2 ))
    allocate (Ice_Layer_Fraction(self % dim1, self % dim2 ))
    allocate (Water_Layer_Fraction(self % dim1, self % dim2 ))
    allocate (Scwater_Layer_Fraction(self % dim1, self % dim2 ))


    Cloud_Geometrical_Thickness = zc_top -Zc_base

     Water_Layer_Fraction = 0.0
     Scwater_Layer_Fraction = 0.0

     ! start with ice layer fraction
     Ice_Layer_Fraction = 1.0

     where ( Zc_Base .LEfp. Upper_Limit_Water_Height )
       Ice_Layer_Fraction = (Zc_top - Upper_Limit_Water_Height) / Cloud_Geometrical_Thickness
     end where

      Water_Layer_Fraction = 1.0
      where ( Zc_top > Upper_Limit_Water_Height)
         Water_Layer_Fraction = ( Upper_Limit_Water_Height - Zc_Base) / &
          Cloud_Geometrical_Thickness
       end where

       Scwater_Layer_Fraction  = 0.0
       where ( Zc_top > Freezing_Level_Height .and. Zc_Base < Freezing_Level_Height)
         Scwater_Layer_Fraction = ( Zc_top - Freezing_Level_Height ) / &
                              Cloud_Geometrical_Thickness

      end where

      where (Zc_top  >  Upper_Limit_Water_Height  .and. Zc_Base < Upper_Limit_Water_Height)
        Scwater_Layer_Fraction = ( Upper_Limit_Water_Height - &
                             Zc_Base ) / &
                         Cloud_Geometrical_Thickness
     end where

      where ( Zc_top > Upper_Limit_Water_Height .and. Zc_Base < Freezing_Level_Height)
        Scwater_Layer_Fraction = (Upper_Limit_Water_Height - &
                                Freezing_Level_Height) / &
                      Cloud_Geometrical_Thickness
       end where

      where (Ice_Layer_Fraction .eq. 1.0 )
        Water_Layer_Fraction = 0.0
        Scwater_Layer_Fraction = 0.0
      end where

      ! get sure nothing is below for unknown read_instr_constants

      where ( Ice_Layer_Fraction .lt. 0. ) Ice_Layer_Fraction = 0.
      where ( Water_Layer_Fraction .lt. 0. ) Water_Layer_Fraction = 0.
      where ( Scwater_Layer_Fraction .lt. 0. ) Scwater_Layer_Fraction = 0.

      self % Cwp_Ice_Layer = Ice_Layer_Fraction * self % cwp
      self % Cwp_Water_Layer= Water_Layer_Fraction * self % cwp
      self % Cwp_Scwater_Layer = Scwater_Layer_Fraction * self % cwp

  end subroutine COMPUTE_CWP_PHASE


end module dcomp_common_mod
