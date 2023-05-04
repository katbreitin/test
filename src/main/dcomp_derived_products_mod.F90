! $Id: dcomp_derived_products_mod.f90 4044 2020-11-10 21:52:49Z heidinger $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: dcomp_derived_products_module.f90 (src)
!       DCOMP_DERIVED_PRODUCTS_MODULE (program)
!
! PURPOSE: this modules hold subroutines that derive additional products
!          primarily from DCOMP
!
! DESCRIPTION:
!
! AUTHORS:
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
! Public routines used in this MODULE:
! COMPUTE_CLOUD_WATER_PATH - compute cloud water path
! COMPUTE_PRECIPITATION - compute precipitation using KNMI approach
! COMPUTE_ADIABATIC_CLOUD_PROPS - compute precipitation using KNMI approach
!
!--------------------------------------------------------------------------------------
MODULE DCOMP_DERIVED_PRODUCTS_MOD

 use CONSTANTS_MOD, only: &
   real4 &
   ,sym &
   , DTOR &
   , PI &
   , missing_value_real4

 use PIXEL_COMMON_MOD,only: &
     acha &
   , base &
   , NWP_PIX &
   , Cld_Type &
   , NLCOMP &
   , DCOMP &
   , DCOMP_1 &
   , DCOMP_2 &
   , DCOMP_3 &
   , bad_pixel_mask &
   , image &
   , sfc &
   , ch  &
   , sensor &
   , geo &
   , Dcomp_Mode &
   , Temp_Pix_Array_1 &
   , Diag_Pix_Array_1 &
   , Diag_Pix_Array_2 &
   , Diag_Pix_Array_3

 use NWP_COMMON_MOD,only: NWP

 use univ_fp_comparison_mod, only: operator(.EQfp.), operator(.NEfp.),  &
      operator(.GEfp.), operator(.LEfp.)

 implicit none
 private
 public:: COMPUTE_CLOUD_WATER_PATH, &
          COMPUTE_PRECIPITATION, &
          COMPUTE_PRECIPITATION_AHI, &
          COMPUTE_DCOMP_INSOLATION, &
          ADJUST_DCOMP_LWP
! private:: COD_APPROX
 private:: BIRD_NREL_CLEAR_SKY_INSOL

 contains
 !-----------------------------------------------------------
 ! compute the optical depth for a sub-pixel observation
 ! within a pixel where we did a cod retrieval
 !
 ! Input: Refl_Clear - clear-sky reflectance
 !        Refl_Retr - reflectance used in cod retrieval
 !        Refl_Asym - asymptotic reflectance from retrival
 !        Refl_Sub - subpixel reflectance
 !        Cod_Retr - cloud optical depth from retrieval
 !        Cod_Asym - asymptotic cloud optical (max) in retrieval
 !
 ! Output:  Cod_Sub = cloud optical depth for Refl_Sub
 !
 ! Approx equation goes like this
 ! Refl_Retr = Refl_Clear + (Refl_Asym - Refl_Clear)*(1.0-exp(alpha*Cod_Retr))
 !-----------------------------------------------------------
 !function COD_APPROX(Refl_Clear,Refl_Retr,Refl_Asym,Refl_Sub, &
 !                    Cod_Retr) result(Cod_Sub)

 !   real, intent(in):: Refl_Clear,Refl_Retr,Refl_Asym,Refl_Sub, &
  !                     Cod_Retr
 !   real:: Cod_Sub
!    real:: alpha
 !   real:: RR

    !Refl_Retr = Refl_Clear + (Refl_Asym - Refl_Clear)*(1.0-exp(alpha*Cod_Retr))
 !   Cod_Sub = Missing_Value_Real4

 !   if (Refl_Retr > Refl_Clear) then

  !     RR = (Refl_Retr - Refl_Clear) / (Refl_Asym - Refl_Clear)
  !     alpha = (alog(1.0 - RR))/Cod_Retr

 !      RR = (Refl_Sub - Refl_Clear) / (Refl_Asym - Refl_Clear)
  !     Cod_Sub = alog(1.0 - RR) / alpha

  !  else

 !      Cod_Sub = 0.0

!    endif


! end function


!-----------------------------------------------------------
! compute cloud water path from the optical depth
! and particle size from the dcomp algorithm
!
! The layer values are computed assuming a linear variation
! in cloud water path from the top to the base of the cloud.
! Note CWP = CWP_Ice_Layer + CWP_Water_Layer and
!      CWP_Scwater is a component of the Water_Layer
!
!-----------------------------------------------------------
subroutine COMPUTE_CLOUD_WATER_PATH()

  integer :: i,j, ii,jj
  integer:: Iphase

  real(kind=real4), parameter:: Rho_Water = 1.0    !g/m^3
  real(kind=real4), parameter:: Rho_Ice = 0.917    !g/m^3

  integer, dimension(:,:), allocatable:: Lat_NWP_Idx
  integer, dimension(:,:), allocatable:: Lon_NWP_Idx
  real, dimension(:,:), allocatable:: Upper_Limit_Water_Height
  real, dimension(:,:), allocatable:: Freezing_Level_Height



  allocate (Upper_Limit_Water_Height(dcomp % dim1, dcomp % dim2 ))
  allocate (Freezing_Level_Height(dcomp % dim1, dcomp % dim2 ))
  allocate (Lat_NWP_Idx(dcomp % dim1, dcomp % dim2 ))
  allocate (Lon_NWP_Idx(dcomp % dim1, dcomp % dim2 ))

  Lon_NWP_Idx = NWP_PIX%I_Nwp
  Lat_NWP_Idx = NWP_PIX%J_Nwp

  do i = 1, dcomp % dim1
    do j = 1, dcomp % dim2

      ii = Lon_Nwp_Idx(i,j)
      jj = Lat_Nwp_Idx(i,j)

      if ( ii .lt. 1) cycle
      if ( jj .lt. 1) cycle
      if ( ii .gt. dcomp % dim1) cycle
      if ( jj .lt. dcomp % dim2) cycle
      Upper_Limit_Water_Height(i,j) = NWP%Upper_Limit_Water_Height(ii,jj)
      Freezing_Level_Height(i,j)= NWP%Freezing_Level_Height(ii,jj)
    end do
  end do

  if (dcomp_1 % is_set) then
    call dcomp_1 % COMPUTE_CWP_PHASE( &
    Acha % Zc, BASE % Zc_base &
    , Upper_Limit_Water_Height &
    , Freezing_Level_Height)
    call dcomp_1 % COMPUTE_ADIABATIC_PROPS (Acha % tc)
  end if

  if (dcomp_2 % is_set) then
    call dcomp_2 % COMPUTE_CWP_PHASE( &
    Acha % Zc, BASE % Zc_base &
    , Upper_Limit_Water_Height &
    , Freezing_Level_Height)
    call dcomp_2 % COMPUTE_ADIABATIC_PROPS (Acha % tc)
  end if

  if (dcomp_3 % is_set) then
    call dcomp_3 % COMPUTE_CWP_PHASE( &
    Acha % Zc, BASE % Zc_base &
    , Upper_Limit_Water_Height &
    , Freezing_Level_Height)
    call dcomp_3 % COMPUTE_ADIABATIC_PROPS (Acha % tc)
  end if

  if (nlcomp % is_set) then
    call nlcomp % COMPUTE_CWP_PHASE( &
    Acha % Zc, BASE % Zc_base &
    , Upper_Limit_Water_Height &
    , Freezing_Level_Height)
    call nlcomp % COMPUTE_ADIABATIC_PROPS (Acha % tc)
  end if

  deallocate (Upper_Limit_Water_Height)
  deallocate (Freezing_Level_Height)
  deallocate (Lat_NWP_Idx)
  deallocate (Lon_NWP_Idx)


end subroutine COMPUTE_CLOUD_WATER_PATH



!---------------------------------------------------------------------------------------
! compute precipitation from the KNMI approach
!
! Citation: Roebeling, R. A., and I. Holleman (2009),
!           SEVIRI rainfall retrieval and validation using weather radar observations,
!           J. Geophys. Res., 114, D21202, doi:10.1029/2009JD012102.
!
! /***************** PROCEDURE Calculate_Precip
!! *********************************
! * Procedure    : Calculate_Precip
! * Description  : Retrieves Precip
!
! * Author       : Rob Roebeling
! * Calls                :
!
!******************************************************************************/
!void Calculate_PRECIP( float *precip, float tau, float reff, int phase, int cch
!)
!{
!  /* - Declarations */
!  float precip_temp, precip_max, lwp_p,dprecip;
!
!  /* Start Sub-Routine */
!
!  /* Parameterization based on Akos Horvarh and Roger Davies, 2007 */
!  lwp_p = tau * reff * 2.0/3.0;
!
!  if (lwp_p >= 150.0 && ((reff >= 15.0 && phase == J_LIQ) || phase == J_ICE) )
!  {
!
!    dprecip    = (float)cch;
!    if ((float)cch > 7000.0)  dprecip  = 7000.0;
!
!    precip_temp = pow((lwp_p/120.0-1.0),1.6)/(dprecip/1000.0);
!    precip_max  = 5.0+pow((float)cch/1000.0,1.6);
!
!    if (precip_temp >= precip_max) precip_temp = precip_max;
!
!    *precip     = precip_temp;
!  }
!  else
!  {
!      *precip  = J_NCL;
!  }
!}
!/* ------------------------------------------------------------------------- */
!---------------------------------------------------------------------------------------
subroutine COMPUTE_PRECIPITATION(Line_Idx_Min,Num_Lines)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines

  real (kind=real4), parameter:: dH_0 = 0.6 !km
  real (kind=real4), parameter:: CWP_0 = 120.0 !g/m^2
  real (kind=real4), parameter:: Lapse_Rate = 6.5 !K/km
  real (kind=real4), parameter:: Alpha = 1.6 !dimensionless
  !real (kind=real4), parameter:: C = 1.0 !mm / hour
  real (kind=real4), parameter:: CWP_T = 150.0 !g/m^2
  real (kind=real4), parameter:: Ceps_T = 15.0 !micron
  integer, parameter:: N_box = 100
  real (kind=real4), parameter:: dH_Max = 7.0  !km
  real (kind=real4) :: CTT_Max
  real (kind=real4) :: CTT_Pix
  real (kind=real4) :: Reff_Pix
  real (kind=real4) :: CWP_Pix
  real (kind=real4) :: dH
  real (kind=real4) :: Rain_Rate_Max
  integer:: Line_Idx_Max
  integer:: Elem_Idx_Max
  integer:: Elem_Idx_Min
  integer:: Num_Elements
  integer:: Line_Idx
  integer:: Line_Idx_1
  integer:: Line_Idx_2
  integer:: Elem_Idx
  integer:: Elem_Idx_1
  integer:: Elem_Idx_2


  Elem_Idx_Min = 1
  Num_Elements = Image%Number_Of_Elements
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  DCOMP % rain_rate = Missing_Value_Real4

  line_loop: DO Line_Idx = Line_Idx_Min, Line_Idx_Max

    Line_Idx_1  = min(Line_Idx_Max-1,max(1,Line_Idx - N_box /2))
    Line_Idx_2  = min(Line_Idx_Max,max(2,Line_Idx + N_box /2))

    element_loop: DO Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

      !--- skip bad pixels
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

      !--- skip non cloud pixels
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%CLEAR_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%PROB_CLEAR_TYPE) then
          DCOMP % rain_rate(Elem_Idx,Line_Idx) = 0.0
          cycle
      endif

      CWP_Pix = DCOMP % cwp(Elem_Idx,Line_Idx)
      CTT_Pix = ACHA%Tc(Elem_Idx,Line_Idx)
      Reff_Pix = DCOMP % reff(Elem_Idx,Line_Idx)

      if (Geo%Solzen(Elem_Idx,Line_Idx) > 90.0) then
        Reff_Pix = Nlcomp % Reff(Elem_Idx,Line_Idx)
      endif

      !--- skip bad pixels
      if (Reff_Pix .LEfp. 0.0) cycle

      !--- screen low water path clouds
      if (CWP_Pix < CWP_T) then
        DCOMP % rain_rate(Elem_Idx,Line_Idx) = 0.0
        cycle
      endif

      !--- screen small particle size liquid water clouds
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%FOG_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%WATER_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%SUPERCOOLED_TYPE) then
         if (Reff_Pix < Ceps_T) then
            DCOMP % rain_rate(Elem_Idx,Line_Idx) = 0.0
            cycle
         endif
      endif

     !--- define box to look for maximum cloud temperature
      Line_Idx_1  = min(Line_Idx_Max-1,max(1,Line_Idx - N_box /2))
      Line_Idx_2  = min(Line_Idx_Max,max(2,Line_Idx + N_box /2))
      Elem_Idx_1  = min(Elem_Idx_Max-1,max(1,Elem_Idx - N_box /2))
      Elem_Idx_2  = min(Elem_Idx_Max,max(2,Elem_Idx + N_box /2))
      CTT_Max = maxval(ACHA%Tc(Elem_Idx_1:Elem_Idx_2,Line_Idx_1:Line_Idx_2))


      !--- compute precip height
      dH = min(dH_Max,(CTT_Max - CTT_Pix) / Lapse_Rate + dH_0)

      !--- compute rain rate
      Rain_Rate_Max = 5.0 + dH**Alpha
      DCOMP % rain_rate(Elem_Idx,Line_Idx) = (((CWP_Pix - CWP_0)/CWP_0)**Alpha) / dH
      DCOMP % rain_rate(Elem_Idx,Line_Idx) = min(Rain_Rate_Max,DCOMP % rain_rate(Elem_Idx,Line_Idx))


    enddo element_loop
  enddo line_loop

end subroutine COMPUTE_PRECIPITATION

!---------------------------------------------------------------------------------------
! Precip version tuned to Himawari-8/AHI
!---------------------------------------------------------------------------------------
subroutine COMPUTE_PRECIPITATION_AHI(Line_Idx_Min,Num_Lines)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines

  real (kind=real4), parameter:: dH_0 = 0.6 !km
  !real (kind=real4), parameter:: CWP_0 = 120.0 !g/m^2
  real (kind=real4), parameter:: Lapse_Rate = 6.5 !K/km
  real (kind=real4), parameter:: Alpha = 1.6 !dimensionless
  !real (kind=real4), parameter:: C = 1.0 !mm / hour
  !real (kind=real4), parameter:: CWP_T = 150.0 !g/m^2
  !real (kind=real4), parameter:: Ceps_T = 15.0 !micron
  integer, parameter:: N_box = 100
  real (kind=real4), parameter:: dH_Max = 7.0  !km
  real (kind=real4) :: CTT_Max
  real (kind=real4) :: CTT_Pix
  real (kind=real4) :: Reff_Pix
  real (kind=real4) :: CWP_Pix
  real (kind=real4) :: dH
  real (kind=real4) :: Rain_Rate_Max
  integer:: Line_Idx_Max
  integer:: Elem_Idx_Max
  integer:: Elem_Idx_Min
  integer:: Num_Elements
  integer:: Line_Idx
  integer:: Line_Idx_1
  integer:: Line_Idx_2
  integer:: Elem_Idx
  integer:: Elem_Idx_1
  integer:: Elem_Idx_2

  !--- Himawari-8 tuned parameters
  !--- Ice Phase parameters
  real (kind=real4), parameter:: CWP_0_Ice = 82.0 !g/m^2
  real (kind=real4), parameter:: Alpha_Ice = 0.97 !dimensionless
  real (kind=real4), parameter:: Ceps_T_Ice = 15.0 !micron for ice clouds
  real (kind=real4), parameter:: CWP_T_Ice = 700.0 !micron for ice clouds

  !--- Liquid Phase parameters
  real (kind=real4), parameter:: CWP_0_Water = 541.2 !g/m^2
  real (kind=real4), parameter:: Alpha_Water = 0.919 !dimensionless
  real (kind=real4), parameter:: Ceps_T_Water = 13.0 !micron for water clouds
  real (kind=real4), parameter:: CWP_T_Water = 1180.0 !micron for water clouds

  Elem_Idx_Min = 1
  Num_Elements = Image%Number_Of_Elements
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  DCOMP % rain_rate = Missing_Value_Real4

  line_loop: DO Line_Idx = Line_Idx_Min, Line_Idx_Max

    Line_Idx_1  = min(Line_Idx_Max-1,max(1,Line_Idx - N_box /2))
    Line_Idx_2  = min(Line_Idx_Max,max(2,Line_Idx + N_box /2))

    element_loop: DO Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

      !--- skip bad pixels
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

      !--- skip non cloud pixels
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%CLEAR_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%PROB_CLEAR_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%UNKNOWN_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%SMOKE_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%DUST_TYPE) then
          DCOMP % rain_rate(Elem_Idx,Line_Idx) = 0.0
          cycle
      endif

      CWP_Pix = DCOMP % cwp(Elem_Idx,Line_Idx)
      CTT_Pix = ACHA%Tc(Elem_Idx,Line_Idx)
      Reff_Pix = DCOMP % reff(Elem_Idx,Line_Idx)

      if (Geo%Solzen(Elem_Idx,Line_Idx) > 90.0) then
        Reff_Pix = NLCOMP % Reff(Elem_Idx,Line_Idx)
      endif

      !--- skip bad pixels
      if (Reff_Pix .LEfp. 0.0) cycle

      !--- screen small particle size and water path for liquid water clouds
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%FOG_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%WATER_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%SUPERCOOLED_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%MIXED_TYPE) then
         if (Reff_Pix < Ceps_T_Water) then
            DCOMP % rain_rate(Elem_Idx,Line_Idx) = 0.0
            cycle
         endif
         if (CWP_Pix < CWP_T_Water) then
            DCOMP % rain_rate(Elem_Idx,Line_Idx) = 0.0
            cycle
         endif
      endif

      !--- screen small particle size and ice path for ice clouds
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%OPAQUE_ICE_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%TICE_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%CIRRUS_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%OVERLAP_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%OVERSHOOTING_TYPE) then
         if (Reff_Pix < Ceps_T_Ice) then
            DCOMP % rain_rate(Elem_Idx,Line_Idx) = 0.0
            cycle
         endif
         if (CWP_Pix < CWP_T_Ice) then
            DCOMP % rain_rate(Elem_Idx,Line_Idx) = 0.0
            cycle
         endif
      endif

     !--- define box to look for maximum cloud temperature
      Line_Idx_1  = min(Line_Idx_Max-1,max(1,Line_Idx - N_box /2))
      Line_Idx_2  = min(Line_Idx_Max,max(2,Line_Idx + N_box /2))
      Elem_Idx_1  = min(Elem_Idx_Max-1,max(1,Elem_Idx - N_box /2))
      Elem_Idx_2  = min(Elem_Idx_Max,max(2,Elem_Idx + N_box /2))
      CTT_Max = maxval(ACHA%Tc(Elem_Idx_1:Elem_Idx_2,Line_Idx_1:Line_Idx_2))


      !--- compute precip height
      dH = min(dH_Max,(CTT_Max - CTT_Pix) / Lapse_Rate + dH_0)

      !--- compute rain rate
      Rain_Rate_Max = 5.0 + dH**Alpha
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%FOG_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%WATER_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%SUPERCOOLED_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%MIXED_TYPE) then
         Rain_Rate_Max = 5.0 + dH**Alpha_Water
         DCOMP % rain_rate(Elem_Idx,Line_Idx) = (((CWP_Pix - CWP_0_Water)/CWP_0_Water)**Alpha_Water) / dH
      endif
      if (Cld_Type(Elem_Idx,Line_Idx) == sym%OPAQUE_ICE_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%TICE_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%CIRRUS_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%OVERLAP_TYPE .or. &
          Cld_Type(Elem_Idx,Line_Idx) == sym%OVERSHOOTING_TYPE) then
         Rain_Rate_Max = 5.0 + dH**Alpha_Ice
         DCOMP % rain_rate(Elem_Idx,Line_Idx) = (((CWP_Pix - CWP_0_Ice)/CWP_0_Ice)**Alpha_Ice) / dH
      endif
      DCOMP % rain_rate(Elem_Idx,Line_Idx) = min(Rain_Rate_Max,DCOMP % rain_rate(Elem_Idx,Line_Idx))


    enddo element_loop
  enddo line_loop

end subroutine COMPUTE_PRECIPITATION_AHI

!---------------------------------------------------------------------------------------
! compute insolation using the cloud properties
!
! Citation: J.A. Coakley (2003),  Reflectance and Albedo, Surface
! https://curry.eas.gatech.edu/Courses/6140/ency/Chapter9/Ency_Atmos/Reflectance_Albedo_Surface.pdf
!
! Progam takes the cloud tranmission and the spherical albedo from DCOMP to estimate
! solar insolation.
!
! Note, the cloud tranmission from DCOMP is a total transmission (diffuse +
! direct).  This routines separates them.
!
! Regression for broad-band solar transmission are taken from (ref here)
!
! Currently, 0.65 um MODIS white sky albedoes are used for the surface albedo
! over land.
!
! Weaknesses
! 1) We need to develop appropriate direct and diffuse broad-band values of surface albedo
! 2) We need to add aerosol impacts
!
!
!---------------------------------------------------------------------------------------
subroutine COMPUTE_DCOMP_INSOLATION(Line_Idx_Min,Num_Lines,Sun_Earth_Distance)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  real, intent(in):: Sun_Earth_Distance

  real (kind=real4), parameter:: SOLAR_CONSTANT = 1356.0 !W/m^2
  real (kind=real4) :: Cloud_Spherical_Albedo
  real (kind=real4) :: Cloud_Optical_Depth
  real (kind=real4) :: Solar_Zenith_Angle
  real (kind=real4) :: Cosine_Solar_Zenith_Angle
  real (kind=real4) :: Cloud_Transmission_Diffuse
  real (kind=real4) :: Cloud_Transmission_Direct
  real (kind=real4) :: Surface_Albedo_Direct
  real (kind=real4) :: Surface_Albedo_Diffuse
  real (kind=real4) :: Insolation_Dcomp_Diffuse_Black_Surface
  !real (kind=real4) :: Insolation_Direct_Dcomp_Black_Surface
  real (kind=real4) :: Insolation_Direct_Dcomp
  real (kind=real4) :: Fo_Toa
  real (kind=real4) :: Fo
  real (kind=real4) :: Tpw
  !real (kind=real4) :: Tozone
  integer:: Line_Idx_Max
  integer:: Elem_Idx_Max
  integer:: Elem_Idx_Min
  integer:: Num_Elements
  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: Land_Class
  real:: Trans_Aerosol
  real:: Trans_H2o
  real:: Trans_O3
  real:: Trans_Molec
  real:: Trans_Ray
  real:: atm_trans
  real:: Surface_Pressure
  real:: Idir_Clear
  real:: Idif_Clear
  real:: aod_380
  real:: aod_550
  real:: O3cm

  Elem_Idx_Min = 1
  Num_Elements = Image%Number_Of_Elements
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  DCOMP % Insolation = Missing_Value_Real4
  DCOMP % Insolation_Diffuse = Missing_Value_Real4

  Fo_Toa = SOLAR_CONSTANT / (Sun_Earth_Distance**2)

  aod_380 = 0.1
  aod_550 = 0.1

  if (Sensor%Chan_On_Flag_Default(1) == sym%NO) return

  line_loop: DO Line_Idx = Line_Idx_Min, Line_Idx_Max
    element_loop: DO Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

      Cloud_Optical_Depth = DCOMP % tau(Elem_Idx,Line_Idx)
      Solar_Zenith_Angle = Geo%Solzen(Elem_Idx,Line_Idx)
      Land_Class = Sfc%Land(Elem_Idx,Line_Idx)

      !--- skip data that can not be processed
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle
      if (Solar_Zenith_Angle > 70.0) cycle

      !--- skip if no nwp
      if (NWP_Pix%Ozone(Elem_Idx,Line_Idx) .EQfp. Missing_Value_Real4) cycle
      if (NWP_Pix%Psfc(Elem_Idx,Line_Idx) .EQfp. Missing_Value_Real4) cycle
      if (NWP_Pix%Tpw(Elem_Idx,Line_Idx) .EQfp. Missing_Value_Real4) cycle

      !--- convert nwp as needed by Bird Routine
      O3cm = 0.001*NWP_Pix%Ozone(Elem_Idx,Line_Idx)   !atm-cm
      Surface_Pressure = NWP_Pix%Psfc(Elem_Idx,Line_Idx) !hPa
      TPW = NWP_PIX%Tpw(Elem_Idx,Line_Idx)

      !--- adjust gases for slant path
      Cosine_Solar_Zenith_Angle = cos(Geo%Solzen(Elem_Idx,Line_Idx)*DTOR)

      !--- determine surface albedo
      if (Land_Class == sym%LAND) then
         Surface_Albedo_Diffuse = ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         Surface_Albedo_Direct = Surface_Albedo_Diffuse
      else
         Surface_Albedo_Diffuse = 0.06
         Surface_Albedo_Direct = 0.026/(Cosine_Solar_Zenith_Angle**1.7+0.065) + &
                                 0.15*(Cosine_Solar_Zenith_Angle - 0.1)* &
                                 (Cosine_Solar_Zenith_Angle-1.0)
      endif

      call BIRD_NREL_CLEAR_SKY_INSOL(Surface_Pressure,Cosine_Solar_Zenith_Angle, &
                                     Surface_Albedo_Direct,O3cm,TPW,Sun_Earth_Distance,&
                                     aod_380,aod_550, &
                                     Trans_Ray, Trans_O3, Trans_Molec, Trans_H2o, Trans_Aerosol, &
                                     Idir_Clear, Idif_Clear)


      atm_trans = Trans_Ray * Trans_O3 * Trans_Molec * Trans_H2o * Trans_Aerosol

      Fo = Fo_Toa * atm_trans

      !-- set cloud trans and albedo, if clear, set to transparent values
      if (Cloud_Optical_Depth > 0.0) then

        Cloud_Spherical_Albedo = DCOMP % Cloud_063um_Spherical_Albedo(Elem_Idx,Line_Idx)
        Cloud_Transmission_Direct = exp( -1.0 * Cloud_Optical_Depth / Cosine_Solar_Zenith_Angle)
        Cloud_Transmission_Diffuse = DCOMP % Cloud_063um_Transmission_Solar(Elem_Idx,Line_Idx)  - Cloud_Transmission_Direct

        Insolation_Direct_Dcomp = Fo * Cloud_Transmission_Direct * Cosine_Solar_Zenith_Angle

        Insolation_Dcomp_Diffuse_Black_Surface = Fo * Cloud_Transmission_Diffuse * Cosine_Solar_Zenith_Angle

        DCOMP % Insolation_Diffuse(Elem_Idx,Line_Idx) = &
                                   Insolation_Direct_Dcomp * Surface_Albedo_Direct * Cloud_Spherical_Albedo / &
                                   (1.0 - Surface_Albedo_Diffuse*Cloud_Spherical_Albedo) + &
                                   Insolation_Dcomp_Diffuse_Black_Surface / &
                                   (1.0 - Surface_Albedo_Diffuse*Cloud_Spherical_Albedo)

      else

        Insolation_Direct_Dcomp = Idir_Clear
        DCOMP % Insolation_Diffuse(Elem_Idx,Line_Idx) = Idif_Clear

      endif

      !--- combine
      DCOMP % Insolation(Elem_Idx,Line_Idx) = Insolation_Direct_Dcomp + DCOMP % Insolation_Diffuse(Elem_Idx,Line_Idx)

    enddo element_loop
  enddo line_loop

end subroutine COMPUTE_DCOMP_INSOLATION

!----------------------------------------------------------------------------------------------
!    adjust lwp to bring in line with AMSR2 by adjust reff
!    coefficients depend on mode
!
! Modified for MODE 9
!----------------------------------------------------------------------------------------------
subroutine ADJUST_DCOMP_LWP()

   real, parameter:: a_mode1 = -0.723323
   real, parameter:: b_mode1 = 0.583924
   real, parameter:: c_mode1 = -0.00420656

   real, parameter:: a_mode2 = -0.517295
   real, parameter:: b_mode2 = 0.601017
   real, parameter:: c_mode2 = -0.00156506

   real, parameter:: a_mode3 = 3.67752
   real, parameter:: b_mode3 = 0.115759
   real, parameter:: c_mode3 = 0.0182053

   real, parameter:: a_mode9 = 2.27057
   real, parameter:: b_mode9 = 0.187615
   real, parameter:: c_mode9 = 0.00923558
   real, parameter:: d_mode9 = 0.311556
   real, parameter:: e_mode9 = -0.135460
   real, parameter:: f_mode9 = -0.0226757
   real, parameter:: g_mode9 = -0.00870485

   real:: a,b,c
!  real, dimension(:,:), pointer:: Reff_Fit

   DCOMP % reff_Fit = DCOMP % reff
   DCOMP % Cwp_Fit = DCOMP % cwp

   if (Dcomp_Mode > 0 .and. Dcomp_Mode /= 9) then

     select case (Dcomp_Mode)
       case(1)
         a = a_mode1
         b = b_mode1
         c = c_mode1
       case(2)
         a = a_mode2
         b = b_mode2
         c = c_mode2
       case(3)
         a = a_mode3
         b = b_mode3
         c = c_mode3
       case default
         a = 0.0
         b = 1.0
         c = 0.0
      end select

      where(Cld_Type < 6 .and. Cld_Type > 1 .and.  &
         (DCOMP % tau .NEfp. MISSING_VALUE_REAL4) .and. &
         (DCOMP % reff .NEfp. MISSING_VALUE_REAL4))

         DCOMP % reff_Fit = a + (b * DCOMP % reff) + (c * DCOMP % reff**2)

      endwhere

   endif

   if (Dcomp_Mode == 9) then

      where(Cld_Type < 6 .and. Cld_Type > 1 .and.  &
         (DCOMP_1 % tau .NEfp. MISSING_VALUE_REAL4) .and.  &
         (DCOMP_2 % tau .NEfp. MISSING_VALUE_REAL4) .and.  &
         (DCOMP_3 % tau .NEfp. MISSING_VALUE_REAL4) .and.  &
         (DCOMP_1 % reff .NEfp. MISSING_VALUE_REAL4) .and. &
         (DCOMP_2 % reff .NEfp. MISSING_VALUE_REAL4) .and. &
         (DCOMP_3 % reff .NEfp. MISSING_VALUE_REAL4))

         DCOMP % reff_Fit = a_mode9 +  &
                    b_mode9 * DCOMP_2 % reff + &
                    c_mode9 * DCOMP_2 % reff**2 + &
                    d_mode9 * (DCOMP_2 % reff - DCOMP_3 % reff) + &
                    e_mode9 * (DCOMP_2 % reff - DCOMP_1 % reff) + &
                    f_mode9 * (DCOMP_2 % reff - DCOMP_3 % reff)**2 + &
                    g_mode9 * (DCOMP_2 % reff - DCOMP_1 % reff)**2

      endwhere


   endif


   where(DCOMP % reff_Fit .NEfp. MISSING_VALUE_REAL4)

         DCOMP % Cwp_Fit = 0.666 * DCOMP % tau * DCOMP % reff_Fit

   end where


end subroutine ADJUST_DCOMP_LWP
!--------------------------------------------------------------------------------------------
!---   reference: https://www.nrel.gov/docs/legosti/old/2436.pdf
! input:
!         Press - Surface Pressure (hPa)
!         Sfc_Albedo (0-1)
!         Cosza - cosine of solar zenith angle
!         Tozone - Total Ozone in Dobson Units
!         Tpw = Total Preciptiable Water (cm)
!         SedCorr = sun earth distance
!         Ta3 = aerosol optical depth at  380nm
!         Ta5 = aerosol optical depth at  550nm
!
! output:
!         Tr = rayleigh transmission
!         Toz = ozone transmission
!         Tw = water vapor transmission
!         Tum = total molecular transmission
!         Ta = aerosol transmission
!         Idh = direct irradiance
!         Idif = diffuse irradiance
!--------------------------------------------------------------------------------------------
subroutine BIRD_NREL_CLEAR_SKY_INSOL(Press,Cosza,Albedo,O3cm,H2Ocm,SedCorr,Ta3,Ta5, &
                                     Tr, Toz, Tum, Tw, Ta, Idh, Idif)

      real, intent(in):: Press
      real, intent(in):: Cosza
      real, intent(in):: Albedo
      real, intent(in):: O3cm
      real, intent(in):: H2Ocm
      real, intent(in):: SedCorr
      real, intent(in):: Ta3
      real, intent(in):: Ta5
      real, intent(out):: Tr
      real, intent(out):: Toz
      real, intent(out):: Tum
      real, intent(out):: Tw
      real, intent(out):: Ta
      real, intent(out):: Idh
      real, intent(out):: Idif
      real:: Ozm, Ba, Io, K1, AM, AMp, Wm, Tau, Ias, Itot, Id, Rs
      real:: TAS, TAA

      Ba=0.85                         !recommended value
      Io=1367.0                       !bird solar constant
      K1=0.1                          !aerosol absorptance
      AM=1.0/(Cosza+0.15/(93.885-Cosza)**1.25) !atmopspheric mass
      AMp = AM*Press / 1013.0  !adjust atm mass for station pressure
      Tr = exp((-0.0903*AMp**0.84)*(1.+ AMp - AMp**1.01))
      Tum = exp(-0.0127 * AMp**0.26)

      Ozm = O3cm*AM
      Toz = 1. -0.1611 * Ozm * (1.+139.48 * Ozm)**(-0.3035) -0.002715 * Ozm/(1. + 0.044 * Ozm + 0.0003 * Ozm**2)

      Wm = AM * H2Ocm
      Tw = 1 - 2.4959 * Wm / ((1+79.034 * Wm)**0.6828 + 6.385 * Wm)

      Tau = 0.2758 * Ta3 + 0.35 * Ta5
      Ta = exp(-Tau**0.873*(1+Tau-Tau**0.7088)*(AM**0.9108))
      TAA = 1 - K1 * (1 - AM + AM**1.06) * (1-Ta)
      TAS = Ta/TAA

      ! Sky albedo
      Rs = 0.0685+(1-Ba)*(1-TAS)

      ! Direct component
      ! Original multiplier of 0.9662 is factor based on Solar Constant of
      ! Io=1353 and Thekakara spectral distribution  (1981) Maxwell and Iqbal
      ! both recommend using 0.9751 instead, representing a Solar Constant
      ! of Io=1367 and WMO-Wherli spectral distribution (1990)
      Id=Io * 0.9751 * Tr * Toz * Tum * Tw * Ta

      ! Direct on horizontal surface
      Idh=Id * Cosza

      ! Diffuse (scattered)
      Ias = Io * Cosza * 0.79 * Toz * Tw * Tum * TAA * (0.5 * (1.0-Tr)+0.85 * (1.0-TAS)) / (1.0-AM+AM**1.02)

      ! Total dif + dir on horizontal
      Itot=(Idh+Ias) / (1.-Albedo*Rs)

      !Diffuse
      Idif=Itot-Idh

end subroutine BIRD_NREL_CLEAR_SKY_INSOL


!-----------------------------------------------------------
! end of MODULE
!-----------------------------------------------------------
end module DCOMP_DERIVED_PRODUCTS_MOD
