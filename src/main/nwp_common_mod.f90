!$Id: nwp_common_mod.f90 4070 2021-01-20 04:01:22Z yli $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: nwp_common.f90 (src)
!       NWP_COMMON_MOD (program)
!
! PURPOSE:  This module holds the radiative transfer quantities needed for
!           the algorithms
!
! DESCRIPTION: 
!           note, there two type of nwp data
!            1- the pressure level data
!            2- the data on different surface grid
!
!            the only data assumed to be a on the surface grid are
!              - surface temperature
!              - Weasd depth
!              - u and v wind speed at 10m
!
!           the surface and pressure level grid may be different
!           i_Nwp, j_Nwp points to a cell in the pressure level data
!
!           In the GFS data, the pressure and surface grids are the same, in the
!           NCEP reanalysis, they differ.
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
! public::
!   CREATE_NWP_ARRAYS
!   DESTROY_NWP_ARRAYS
!   FIND_NWP_GRID_CELL
!   MAP_PIXEL_NWP
!   KNOWING_P_COMPUTE_T_Z_NWP
!   KNOWING_Z_COMPUTE_T_P_NWP
!   KNOWING_T_COMPUTE_P_Z_NWP
!   FIND_NWP_LEVELS
!   INTERPOLATE_NWP
!   INTERPOLATE_PROFILE
!   INTERPOLATE_NWP_TZ_PROFILES
!   COMPUTE_COAST_MASK_NWP
!   QC_NWP
!   COMPUTE_PIXEL_NWP_PARAMETERS
!   MODIFY_TSFC_NWP_PIX
!   COMPUTE_NWP_PARAMETERS
!   TEMPORAL_INTERP_TMPSFC_NWP
!   CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY
!--------------------------------------------------------------------------------------
module NWP_COMMON_MOD
   
  use CONSTANTS_MOD
  use PIXEL_COMMON_MOD
  use NUMERICAL_ROUTINES_MOD
  use CX_SCIENCE_TOOLS_MOD

  implicit none
  private
  private:: FIND_NWP_LEVELS, &
            COMPUTE_NWP_CLOUD_PARAMETERS, &
            CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY, &
            COMPUTE_Rh300, &
            COMPUTE_UTH, &
            COMPUTE_DIVERGENCE

  public:: CREATE_NWP_ARRAYS,  &
           INITIALIZE_NWP_ARRAYS,  &
           DESTROY_NWP_ARRAYS,  &
           FIND_NWP_GRID_CELL,  &
           MAP_PIXEL_NWP, &
           KNOWING_P_COMPUTE_T_Z_NWP,  &
           KNOWING_Z_COMPUTE_T_P_NWP,  &
           KNOWING_Z_COMPUTE_T_P_NWP_ARBITRARY_LEVELS,  &
           KNOWING_T_COMPUTE_P_Z_NWP,  &
           COMPUTE_NWP_LEVELS_SEGMENT, &
           INTERPOLATE_NWP,  &
           INTERPOLATE_PROFILE,  &
           INTERPOLATE_NWP_TZ_PROFILES,  &
           COMPUTE_COAST_MASK_NWP, &
           QC_NWP, &
           MODIFY_TSFC_NWP_PIX, &
           COMPUTE_SEGMENT_NWP_CLOUD_PARAMETERS, &
           TEMPORAL_INTERP_TMPSFC_NWP, &
           COMPUTE_PIXEL_NWP_PARAMETERS

 interface CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY
     module procedure  &
         CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I1, &
         CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I2, &
         CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I4, &
         CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_R4
 end interface

 interface INTERPOLATE_NWP
     module procedure  &
         INTERPOLATE_NWP_I1, &
         INTERPOLATE_NWP_I2, &
         INTERPOLATE_NWP_I4, &
         INTERPOLATE_NWP_R4
 end interface
!----------------------------------------------------------------------
!--- set this parameter to 1 when reading GFS hdf files that have
!--- x as the first index, not z
!----------------------------------------------------------------------
  integer, public, save:: REFORMAT_GFS_ZXY

! NWP array declarations
  integer (kind=int4), save, public :: npoints
  real (kind=real4),   save, public :: dLat_Nwp, dLon_Nwp
  real (kind=real4), public :: lat1_Nwp
  real (kind=real4), public :: lon1_Nwp
  real (kind=real4), save, public :: missing_Nwp
  real, public, parameter :: Psfc_max_Nwp = 1100.0
  real (kind=real8), public, save :: ncep_time_Before
  real (kind=real8), public, save :: ncep_time_After

  real (kind=real4), dimension(:,:), allocatable, public, save :: Tmpsfc_Nwp_Before
  real (kind=real4), dimension(:,:), allocatable, public, save :: Tmpsfc_Nwp_After
  real (kind=real4), dimension(:,:), allocatable, save, public :: Temp2d_Nwp_1
  real (kind=real4), dimension(:,:), allocatable, save, public :: Temp2d_Nwp_2
  real (kind=real4), dimension(:,:,:), allocatable, save, public :: Temp3d_Nwp_1
  real (kind=real4), dimension(:,:,:), allocatable, save, public :: Temp3d_Nwp_2
  real (kind=real4), dimension(:,:,:), allocatable, save, public :: Temp3d
  real (kind=real4), dimension(:), allocatable, save, public :: Temp1d_Nwp

  integer(kind=int4), save, public:: nwp_start_hour
  integer(kind=int4), save, public:: nwp_end_hour

  !--- nwp profiles interpolated to the pixel level
  real (kind=real4), dimension(:), allocatable, private, save :: T_Prof
  real (kind=real4), dimension(:), allocatable, private, save :: Z_Prof

  !--- local parameters
  real(kind=real4), public, parameter :: P_Trop_Max = 300.0
  real(kind=real4), public, parameter :: P_Trop_Min = 25.0
  real(kind=real4), public, parameter :: P_Inversion_Min = 700.0
  real(kind=real4), public, parameter :: Delta_T_Inversion = 0.0  !05

  type, public :: nwp_definition

    !--- general nwp metadata
    integer (kind=int4) :: Nlevels
    integer (kind=int4) :: Nlon
    integer (kind=int4) :: Nlat
    real (kind=real4), dimension(:), allocatable :: Lon
    real (kind=real4), dimension(:), allocatable :: Lat
    integer (kind=int4), dimension(:,:), allocatable :: Bad_Nwp_Mask

    !--- 3d profiles of atmosphere 
    real (kind=real4), dimension(:), allocatable :: P_Std   !actually 1d
    real (kind=real4), dimension(:,:,:), allocatable :: Z_Prof
    real (kind=real4), dimension(:,:,:), allocatable :: T_Prof
    real (kind=real4), dimension(:,:,:), allocatable :: Rh_Prof
    real (kind=real4), dimension(:,:,:), allocatable :: Ozone_Prof
    real (kind=real4), dimension(:,:,:), allocatable :: Tpw_Prof
    real (kind=real4), dimension(:,:,:), allocatable :: Clwmr_Prof
    real (kind=real4), dimension(:,:,:), allocatable :: Icmr_Prof
    real (kind=real4), dimension(:,:,:), allocatable :: U_Wnd_Prof
    real (kind=real4), dimension(:,:,:), allocatable :: V_Wnd_Prof
    integer (kind=int4), dimension(:,:,:), allocatable :: Inversion_Level_Profile

    !--- 2d surface properties
    integer (kind=int1), dimension(:,:), allocatable :: Land
    real (kind=real4), dimension(:,:), allocatable :: Psfc
    real (kind=real4), dimension(:,:), allocatable :: Pmsl
    real (kind=real4), dimension(:,:), allocatable :: Zsfc
    real (kind=real4), dimension(:,:), allocatable :: Tmpsfc
    real (kind=real4), dimension(:,:), allocatable :: Weasd
    real (kind=real4), dimension(:,:), allocatable :: Sea_Ice_Frac
    real (kind=real4), dimension(:,:), allocatable :: Freezing_lev_hgt

    !--- 2d atmospheric properties
    real (kind=real4), dimension(:,:), allocatable :: Ozone
    real (kind=real4), dimension(:,:), allocatable :: P_Trop
    real (kind=real4), dimension(:,:), allocatable :: T_Trop
    real (kind=real4), dimension(:,:), allocatable :: Z_Trop
    real (kind=real4), dimension(:,:), allocatable :: Tpw
    real (kind=real4), dimension(:,:), allocatable :: Tmpair
    real (kind=real4), dimension(:,:), allocatable :: Rhsfc
    real (kind=real4), dimension(:,:), allocatable :: Rh300
    real (kind=real4), dimension(:,:), allocatable :: Uth
    real (kind=real4), dimension(:,:), allocatable :: K_Index
    real (kind=real4), dimension(:,:), allocatable :: U_Wnd_10m
    real (kind=real4), dimension(:,:), allocatable :: V_Wnd_10m
    real (kind=real4), dimension(:,:), allocatable :: Wnd_Spd_10m
    real (kind=real4), dimension(:,:), allocatable :: Wnd_Dir_10m
    real (kind=real4), dimension(:,:), allocatable :: Div_Sfc
    real (kind=real4), dimension(:,:), allocatable :: Div_200
    real (kind=real4), dimension(:,:), allocatable :: CAPE   

    !--- 2d nwp levels
    integer (kind=int1), dimension(:,:), allocatable :: Tropo_Level
    integer (kind=int1), dimension(:,:), allocatable :: Level850
    integer (kind=int1), dimension(:,:), allocatable :: Level700
    integer (kind=int1), dimension(:,:), allocatable :: Level500
    integer (kind=int1), dimension(:,:), allocatable :: Level300
    integer (kind=int1), dimension(:,:), allocatable :: Level200
    integer (kind=int1), dimension(:,:), allocatable :: Level100
    integer (kind=int1), dimension(:,:), allocatable :: Sfc_Level
    integer (kind=int1), dimension(:,:), allocatable :: Inversion_Level
    real (kind=real4), dimension(:,:), allocatable :: Lifting_Condensation_Level_Height
    real (kind=real4), dimension(:,:), allocatable :: Convective_Condensation_Level_Height
    real (kind=real4), dimension(:,:), allocatable :: Level_Free_Convection_Height
    real (kind=real4), dimension(:,:), allocatable :: Equilibrium_Level_Height
    real (kind=real4), dimension(:,:), allocatable :: Inversion_Strength
    real (kind=real4), dimension(:,:), allocatable :: Inversion_Base
    real (kind=real4), dimension(:,:), allocatable :: Inversion_Top
    real (kind=real4), dimension(:,:), allocatable :: Freezing_Level_Pressure
    real (kind=real4), dimension(:,:), allocatable :: Freezing_Level_Height !km
    real (kind=real4), dimension(:,:), allocatable :: Homogenous_Freezing_Level_Pressure
    real (kind=real4), dimension(:,:), allocatable :: Homogenous_Freezing_Level_Height !km
    real (kind=real4), dimension(:,:), allocatable :: Upper_Limit_Water_Height !km

    !--- 2d nwp cloud properties
    integer (kind=int1), dimension(:,:), allocatable :: Cld_Type
    integer (kind=int1), dimension(:,:), allocatable :: Ncld_Layers
    real (kind=real4), dimension(:,:), allocatable :: Sc_Lwp
    real (kind=real4), dimension(:,:), allocatable :: Lwp
    real (kind=real4), dimension(:,:), allocatable :: Iwp
    real (kind=real4), dimension(:,:), allocatable :: Cwp
    real (kind=real4), dimension(:,:), allocatable :: Pc
    real (kind=real4), dimension(:,:), allocatable :: Tc
    real (kind=real4), dimension(:,:), allocatable :: Cloud_Fraction_Satellite
    real (kind=real4), dimension(:,:), allocatable :: High_Cloud_Fraction_Satellite
    real (kind=real4), dimension(:,:), allocatable :: Mid_Cloud_Fraction_Satellite
    real (kind=real4), dimension(:,:), allocatable :: Low_Cloud_Fraction_Satellite

  end type  nwp_definition

  type(nwp_definition), public, save, target :: NWP

    integer:: Nwp_Opt
    integer:: Smooth_Nwp_Flag
    integer, allocatable, dimension(:,:):: I_Nwp
    integer, allocatable, dimension(:,:):: J_Nwp
    integer, allocatable, dimension(:,:):: I_Nwp_x
    integer, allocatable, dimension(:,:):: J_Nwp_x
    real, allocatable, dimension(:,:):: Lon_Nwp_Fac
    real, allocatable, dimension(:,:):: Lat_Nwp_Fac
    real (kind=real4), dimension(:,:), allocatable:: Tsfc
    real (kind=real4), dimension(:,:), allocatable:: Tair

contains
!-------------------------------------------------------------
! subroutine QC_NWP()
!
! Subroutine to quality control NWP data
!
! Check the values of some fields and set NWP%Bad_Nwp_Mask
! accordingly
!
! The tests run here are arbitrary but are based on known
! failures
!
!-------------------------------------------------------------
 subroutine QC_NWP()

  integer:: Lon_Nwp_Idx, Lat_Nwp_Idx

  
  do Lon_Nwp_Idx = 1, NWP%Nlon
     do Lat_Nwp_Idx = 1, NWP%Nlat

        if ((NWP%P_Trop(Lon_Nwp_Idx,Lat_Nwp_Idx) <= 0.0) .or. &
            (NWP%T_Trop(Lon_Nwp_Idx,Lat_Nwp_Idx) <= 0.0) .or. &
            (NWP%Zsfc(Lon_Nwp_Idx,Lat_Nwp_Idx) > 10000.0) .or. &
            (NWP%Psfc(Lon_Nwp_Idx,Lat_Nwp_Idx) > 1500.0) .or. &
            (NWP%Tmpsfc(Lon_Nwp_Idx,Lat_Nwp_Idx) > 400.0) .or. &
            (NWP%Tmpsfc(Lon_Nwp_Idx,Lat_Nwp_Idx) <= 0.0)) then

            NWP%Bad_Nwp_Mask(Lon_Nwp_Idx,Lat_Nwp_Idx) = sym%YES

        endif

     enddo
  enddo

end subroutine QC_NWP

!----------------------------------------------------------------------
! subroutine COMPUTE_TSFC_NWP(i1,nx,j1,ny,Smooth_Nwp_Opt)
!
! compute a pixel level surface temperature from the NWP fields
! and smooth if option chosen
!
! i1 = first element index
! nx = number of element indices
! j1 = first element index
! ny = number of element indices
! Smooth_Nwp_Opt = flag to smooth nwp
!
! This must be called After MAP_PIXEL_NWP
!----------------------------------------------------------------------
subroutine MODIFY_TSFC_NWP_PIX(Elem_Idx_Start,Num_Elements,Line_Idx_Start,Num_Lines)

  integer(kind=int4), intent(in):: Elem_Idx_Start
  integer(kind=int4), intent(in):: Num_Elements
  integer(kind=int4), intent(in):: Line_Idx_Start
  integer(kind=int4), intent(in):: Num_Lines
  integer(kind=int4) :: Elem_Idx_End
  integer(kind=int4) :: Line_Idx_End
  integer(kind=int4) :: Elem_Idx
  integer(kind=int4) :: Line_Idx
  real (kind=real4) :: Delta_Zsfc
  real (kind=real4) :: Delta_Tsfc
  real (kind=real4) :: Delta_Lapse_Rate
  real(kind=real4) :: Zsfc
  integer(kind=int4) :: Ilev_start
  integer(kind=int4) :: Ilev_end
  integer(kind=int4) :: Lon_Nwp_Idx
  integer(kind=int4) :: Lat_Nwp_Idx
  integer(kind=int4) :: Lon_Nwp_Idx_x
  integer(kind=int4) :: Lat_Nwp_Idx_x
  integer(kind=int4) :: Sfc_Level_Idx

  Elem_Idx_End = Elem_Idx_Start + Num_Elements - 1
  Line_Idx_End = Line_Idx_Start + Num_Lines - 1

  do Elem_Idx = Elem_Idx_Start, Elem_Idx_End
    do Line_Idx = Line_Idx_Start,Line_Idx_End

     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
             NWP_PIX%Tsfc(Elem_Idx,Line_Idx) = Missing_Value_Real4
             cycle
     endif

     Lon_Nwp_Idx = NWP_PIX%i_Nwp(Elem_Idx,Line_Idx)
     Lat_Nwp_Idx = NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)

     if (Lon_Nwp_Idx == 0 .or. Lat_Nwp_Idx == 0) then
             NWP_PIX%Tsfc(Elem_Idx,Line_Idx) = Missing_Value_Real4
             cycle
     endif

     Lon_Nwp_Idx_x = NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx)
     Lat_Nwp_Idx_x = NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)
     Sfc_Level_Idx = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)

     !----------------------------------------------------------------------------------
     ! modify Tsfc for sub-nwp elevation
     !
     !  Zsfc = pixel level elevation in meters
     !  Zsfc_Nwp = nwp level elevation in km
     !
     !----------------------------------------------------------------------------------
     if (Sfc%Land(Elem_Idx,Line_Idx) == sym%LAND) then

        !--- assume all surface features are in lowest half of profile
        Ilev_end = NWP%Nlevels
        Ilev_start = NWP%Nlevels/2

        if ((NWP%Zsfc(Lon_Nwp_Idx,Lat_Nwp_Idx) /= Missing_Value_Real4) .and. &
           (Sfc%Zsfc(Elem_Idx,Line_Idx) /= Missing_Value_Real4) .and. &
           (Lon_Nwp_Idx > 0) .and. (Lat_Nwp_Idx > 0) .and. (Lon_Nwp_Idx_x > 0) .and. (Lat_Nwp_Idx_x > 0) .and. &
           (Sfc_Level_Idx > 1)) then

          !--- compute a smooth surface elevation from NWP 
          Zsfc = INTERPOLATE_NWP( &
                       NWP%Zsfc(Lon_Nwp_Idx,Lat_Nwp_Idx), &
                       NWP%Zsfc(Lon_Nwp_Idx_x,Lat_Nwp_Idx), &
                       NWP%Zsfc(Lon_Nwp_Idx,Lat_Nwp_Idx_x), &
                       NWP%Zsfc(Lon_Nwp_Idx_x,Lat_Nwp_Idx_x), &
                       NWP_PIX%Lon_Nwp_fac(Elem_Idx,Line_Idx), &
                       NWP_PIX%Lat_Nwp_fac(Elem_Idx,Line_Idx))
        
          !--- compute the near surface lapse rate (K/m) 
          Delta_Lapse_Rate = (NWP%T_Prof(Sfc_Level_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx) &
                              - NWP%T_Prof(Sfc_Level_Idx-1,Lon_Nwp_Idx,Lat_Nwp_Idx)) / &
                            (NWP%Z_Prof(Sfc_Level_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx) &
                            - NWP%Z_Prof(Sfc_Level_Idx-1,Lon_Nwp_Idx,Lat_Nwp_Idx))
        else
          Delta_Lapse_Rate = 0
        endif

        !--- compute the pertubation to NWP surface temp to account for sub-grid elevation
        Delta_Zsfc = Sfc%Zsfc(Elem_Idx,Line_Idx) - Zsfc !meters
        Delta_Tsfc = Delta_Lapse_Rate * Delta_Zsfc       !K
        NWP_PIX%Tsfc(Elem_Idx,Line_Idx) = NWP_PIX%Tsfc(Elem_Idx,Line_Idx) + Delta_Tsfc   !K

     endif

    enddo
  enddo   

end subroutine MODIFY_TSFC_NWP_PIX

!----------------------------------------------------------------------
! subroutine COMPUTE_NWP_PARAMETERS(Smooth_Nwp_Opt)
!
! compute parameters from NWP fields and smooth if option chosen
!
! This must be called After MAP_PIXEL_NWP
!----------------------------------------------------------------------
subroutine COMPUTE_PIXEL_NWP_PARAMETERS(Smooth_Nwp_Opt)

  integer(kind=int4), intent(in):: Smooth_Nwp_Opt
  integer:: IDiag
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Tmpsfc,NWP_PIX%Tsfc,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%T_Trop,NWP_PIX%Ttropo,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Z_Trop,NWP_PIX%Ztropo,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%P_Trop,NWP_PIX%Ptropo,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Freezing_Level_Pressure,NWP_PIX%FrzPre,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Homogenous_Freezing_Level_Pressure,NWP_PIX%HomoFrzPre,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Tmpair,NWP_PIX%Tair,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Rhsfc,NWP_PIX%Rhsfc,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Rh300,NWP_PIX%Rh300,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Uth,NWP_PIX%Uth,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Psfc,NWP_PIX%Psfc,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Pmsl,NWP_PIX%Pmsl,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Weasd,NWP_PIX%Weasd,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Sea_Ice_Frac,NWP_PIX%Sea_Ice_Frac,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Tpw,NWP_PIX%Tpw,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Ozone,NWP_PIX%Ozone,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%K_Index,NWP_PIX%K_Index,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Pmsl,NWP_PIX%Pmsl,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Sc_Lwp,NWP_PIX%Sc_Lwp,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Lwp,NWP_PIX%Lwp,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Iwp,NWP_PIX%Iwp,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Cwp,NWP_PIX%Cwp,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Pc,NWP_PIX%Pc,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Tc,NWP_PIX%Tc,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Cloud_Fraction_Satellite,NWP_PIX%Cfrac,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Ncld_Layers,NWP_PIX%Ncld_Layers,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Cld_Type,NWP_PIX%Cld_Type,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Wnd_Spd_10m,NWP_PIX%Wnd_Spd_10m,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Wnd_Dir_10m,NWP_PIX%Wnd_Dir_10m,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Div_Sfc,NWP_PIX%Div_Sfc,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Div_200,NWP_PIX%Div_200,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Lifting_Condensation_Level_Height,NWP_PIX%LCL_Height,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Convective_Condensation_Level_Height,NWP_PIX%CCL_Height,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Level_Free_Convection_Height,NWP_PIX%LFC_Height,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Equilibrium_Level_Height,NWP_PIX%EL_Height,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Inversion_Strength,NWP_PIX%Inversion_Strength,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Inversion_Base,NWP_PIX%Inversion_Base,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%Inversion_Top,NWP_PIX%Inversion_Top,Smooth_Nwp_Opt)
  call CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY(NWP%CAPE,NWP_PIX%CAPE,Smooth_Nwp_Opt)

end subroutine COMPUTE_PIXEL_NWP_PARAMETERS

!-------------------------------------------------------------
! subroutine MAP_PIXEL_NWP(j1,j2)
!
! Subroutine to find nwp cell where each in a segment lies
!-------------------------------------------------------------
 subroutine MAP_PIXEL_NWP(Number_of_Elements,Number_of_Lines)

  integer, intent(in):: Number_of_Elements,Number_of_Lines
  integer:: Elem_Idx,Line_Idx,Ierr


  NWP_PIX%i_Nwp = Missing_Value_Int1
  NWP_PIX%j_Nwp = Missing_Value_Int1
  NWP_PIX%i_Nwp_x = Missing_Value_Int1
  NWP_PIX%j_Nwp_x = Missing_Value_Int1
  NWP_PIX%Lon_Nwp_Fac = Missing_Value_Real4
  NWP_PIX%Lat_Nwp_Fac = Missing_Value_Real4

  do Line_Idx = 1, Number_of_Lines
     do Elem_Idx = 1, Number_of_Elements
                                                                                                                                         
      !--- check for valid geolocation
      if (Nav%Lon(Elem_Idx,Line_Idx) < -180.0 .or. Nav%Lon(Elem_Idx,Line_Idx) > 180.0 .or. &
          Nav%Lat(Elem_Idx,Line_Idx) < -90.0 .or. Nav%Lat(Elem_Idx,Line_Idx) > 90.0) then
          cycle
      endif 
        
      !--- compute NWP cell to pixel mapping
      call FIND_NWP_GRID_CELL(Nav%Lon(Elem_Idx,Line_Idx),Nav%Lat(Elem_Idx,Line_Idx), &
                              NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx), &
                              NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx),  &
                              NWP_PIX%Lon_Nwp_fac(Elem_Idx,Line_Idx), NWP_PIX%Lat_Nwp_fac(Elem_Idx,Line_Idx),Ierr)

       !-- if there is an error, flag pixel as bad
      if (Ierr == 1) then
         NWP_PIX%i_Nwp(Elem_Idx,Line_Idx) = Missing_Value_int1
         NWP_PIX%j_Nwp(Elem_Idx,Line_Idx) = Missing_Value_int1
         NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx) = Missing_Value_int1
         NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx) = Missing_Value_int1
         NWP_PIX%Lon_Nwp_fac(Elem_Idx,Line_Idx) = Missing_Value_Real4
         NWP_PIX%Lat_Nwp_fac(Elem_Idx,Line_Idx) = Missing_Value_Real4
         cycle
      endif

     enddo
  enddo

 end subroutine MAP_PIXEL_NWP
!------------------------------------------------------------------
! Compute NWP Levels for each NWP Gridcell
!
! must be called aftrer MAP_PIXEL_NWP
!------------------------------------------------------------------
 subroutine COMPUTE_NWP_LEVELS_SEGMENT(Number_of_Elements,Number_of_Lines)

  integer, intent(in):: Number_of_Elements
  integer, intent(in):: Number_of_Lines
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Lat_NWP_Idx
  integer:: Lon_NWP_Idx

  !--- intialize levels to missing
  NWP%Sfc_Level = Missing_Value_Int1
  NWP%Tropo_Level = Missing_Value_Int1
  NWP%Level850 = Missing_Value_Int1
  NWP%Level700 = Missing_Value_Int1
  NWP%Level500 = Missing_Value_Int1
  NWP%Level300 = Missing_Value_Int1
  NWP%Level200 = Missing_Value_Int1
  NWP%Level100 = Missing_Value_Int1

  !--- loop through each pixel and if the
  do Line_Idx = 1, Number_of_Lines
     do Elem_Idx = 1, Number_of_Elements

      !--- alias nwp indices for this nwp cell using predetermined global variables
      Lon_NWP_Idx = NWP_PIX%I_Nwp(Elem_Idx,Line_Idx)   
      Lat_NWP_Idx = NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)  

      !--- check for valid nwp mapping and data, if not skip
      if (Lon_NWP_Idx < 1 .or. Lat_Nwp_Idx < 1) cycle
      if (NWP%Bad_Nwp_Mask(Lon_NWP_Idx,Lat_NWP_Idx) == sym%YES) cycle

      !-- if this populated for this nwp cell, skip this pixel
      if (NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) /= Missing_Value_Int1) cycle

      !--- find needed nwp levels for this nwp cell, store in global variables
      call FIND_NWP_LEVELS(Lon_Nwp_Idx,Lat_Nwp_Idx)

     enddo
  enddo

 end subroutine COMPUTE_NWP_LEVELS_SEGMENT
!------------------------------------------------------------------
! subroutine FIND_NWP_GRID_CELL(lon, lat, Lon_Nwp_Idx, Lat_Nwp_Idx, Lon_Nwp_Idxx, Lat_Nwp_Idxx, lonfac, latfac, ierror)
!
! Subroutine to convert lat, lon into NWP grid-cell
!
! input:
!   lon - longitude (-180 to 180)
!   lat - latitude (-90 to 90)
! output:
!   Lat_Nwp_Idx - nwp latitude index of the nearest nwp latitude
!   Lon_Nwp_Idx - nwp longitude index of the nearest nwp longitude
!   Lat_Nwp_Idxx - nwp latitude index of the nearest nwp latitude diagonal
!   Lon_Nwp_Idxx - nwp longitude index of the nearest nwp longitude diagonal
!   latfac - latitude weight between Lat_Nwp_Idx(0.0) and Lat_Nwp_Idxx(1.0)
!   lonfac - longitude weight between Lon_Nwp_Idx(0.0) and Lon_Nwp_Idxx(1.0)
!  
!
! imagine a pixel, x, surrounded by nwp vertices (O)
!
!        O            O            O
!                     -   x
!
!        O            O            O
!                                 --- 
!
!  point O is (Lon_Nwp_Idx,Lat_Nwp_Idx) and O is (Lon_Nwp_Idxx,Lat_Nwp_Idxx)
!        -                                 ---
!
! modified to return information needed to spatially interpolate
!------------------------------------------------------------------
  subroutine FIND_NWP_GRID_CELL(lon, lat, Lon_Nwp_Idx, Lat_Nwp_Idx, Lon_Nwp_Idxx, Lat_Nwp_Idxx, lonfac, latfac, ierror)

    real (kind=real4), intent(in) :: lon, lat
    integer (kind=int4), intent(out) :: Lon_Nwp_Idx, Lat_Nwp_Idx, Lon_Nwp_Idxx, Lat_Nwp_Idxx, ierror
    real (kind=real4), intent(out) :: latfac
    real (kind=real4), intent(out) :: lonfac

    real (kind=real4) :: rlat
    real (kind=real4) :: rlon
    real (kind=real4) :: rLon_Nwp
    real (kind=real4) :: rLon_Nwpx
    integer:: Is_Dateline

    ierror = 0
    rlon = lon
    rlat = lat
    Is_Dateline = 0

    !--- convert negative lons to go from 180 to 360 degrees
    if (rlon < 0.0) then
       rlon = rlon + 360.0
    endif

    !--- Find Position in NWP grid
    if (rlon < 0.0 .or. rlon > 360.0 .or. rlat < -90.0 .or. rlat > 90.0) then
       ierror = 1
       Lat_Nwp_Idx = 0
       Lon_Nwp_Idx = 0
    else
       Lat_Nwp_Idx = max(1, min(NWP%Nlat, nint( (rlat-lat1_Nwp) / dLat_Nwp + 1.0) ))
       Lon_Nwp_Idx = max(1, min(NWP%Nlon, nint( (rlon-lon1_Nwp) / dLon_Nwp + 1.0) ))
    endif

    rLon_Nwp = NWP%Lon(Lon_Nwp_Idx)
    if (NWP%Lon(Lon_Nwp_Idx) < 0.0) then
       rLon_Nwp = rLon_Nwp + 360.0
    endif

    !---  latitude interpolation information
    if (Lat_Nwp_Idx > 1 .and. Lat_Nwp_Idx < NWP%Nlat) then
      if (sign(1.0,lat-NWP%Lat(Lat_Nwp_Idx)) == sign(1.0,dLat_Nwp)) then
         Lat_Nwp_Idxx = Lat_Nwp_Idx + 1
      else
         Lat_Nwp_Idxx = Lat_Nwp_Idx - 1
      endif
      Lat_Nwp_Idxx = min(NWP%Nlat,max(1,Lat_Nwp_Idxx))

      !--- compute latitude interpolation factor
      if (NWP%Lat(Lat_Nwp_Idxx) /= NWP%Lat(Lat_Nwp_Idx)) then
         latfac = (lat - NWP%Lat(Lat_Nwp_Idx))/(NWP%Lat(Lat_Nwp_Idxx)-NWP%Lat(Lat_Nwp_Idx))
      else
         latfac = 0.0
      endif
      
    endif

    !---- determine dateline flag
    if (abs(lon - NWP%Lon(Lon_Nwp_Idx)) > abs(dLon_Nwp)) then
       Is_Dateline = 1
    endif

    !---  longitude interpolation information
    if (Lon_Nwp_Idx > 1 .and. Lon_Nwp_Idx < NWP%Nlon) then
      if (sign(1.0,lon-NWP%Lon(Lon_Nwp_Idx)) == sign(1.0,dLon_Nwp)) then
          if (Is_Dateline == 0) then
            Lon_Nwp_Idxx = Lon_Nwp_Idx + 1
          else
            Lon_Nwp_Idxx = Lon_Nwp_Idx - 1
          endif
      else
          if (Is_Dateline == 0) then
            Lon_Nwp_Idxx = Lon_Nwp_Idx - 1
          else
            Lon_Nwp_Idxx = Lon_Nwp_Idx + 1
          endif
      endif
    endif
    Lon_Nwp_Idxx = min(NWP%Nlon,max(1,Lon_Nwp_Idxx))

    !--- make a positive definite value of lon at Lon_Nwp_Idxx
    rLon_Nwpx = NWP%Lon(Lon_Nwp_Idxx)
    if (NWP%Lon(Lon_Nwp_Idxx) < 0.0) then
       rLon_Nwpx = rLon_Nwpx + 360.0
    endif

    !--- recompute date line flag including Lon_Nwp_Idxx point
    if (abs(lon - NWP%Lon(Lon_Nwp_Idx)) > abs(dLon_Nwp)) then
       Is_Dateline = 1
    endif
    if (abs(lon - NWP%Lon(Lon_Nwp_Idxx)) > abs(dLon_Nwp)) then
       Is_Dateline = 1
    endif

    !--- compute latitude interpolation factor
    if (Is_Dateline == 0) then
       if (NWP%Lon(Lon_Nwp_Idxx)/=NWP%Lon(Lon_Nwp_Idx)) then
          lonfac = (lon - NWP%Lon(Lon_Nwp_Idx))/(NWP%Lon(Lon_Nwp_Idxx)-NWP%Lon(Lon_Nwp_Idx))
       else
          lonfac = 0.0
       endif
    else
       if (rLon_Nwp /= rLon_Nwpx) then
          lonfac = abs((rlon - rLon_Nwp)/(rLon_Nwp-rLon_Nwpx))
       else
          lonfac = 0.0
       endif
    endif

    !--- constrain
    lonfac = min(0.5,max(0.0,lonfac))
    latfac = min(0.5,max(0.0,latfac))
    Lon_Nwp_Idxx = min(NWP%Nlon,max(1,Lon_Nwp_Idxx))
    Lat_Nwp_Idxx = min(NWP%Nlat,max(1,Lat_Nwp_Idxx))

  end subroutine FIND_NWP_GRID_CELL


!----------------------------------------------------------------------------
! Function INTERPOLATE_NWP
! 
! general interpolation routine for nwp fields
!
! description of arguments
! Lon_Nwp_Idx, Lat_Nwp_Idx - nwp indices of closest nwp cell
! Lon_Nwp_Idxx,Lat_Nwp_Idxx - nwp indices of nwp cells of diagnoal of bounding box
! lonx - longitude weighting factor 
! latx = latitude weighting factor
! z1 = data(Lon_Nwp_Idx, Lat_Nwp_Idx)
! z2 = data(Lon_Nwp_Idxx,Lat_Nwp_Idx)
! z3 = data(Lon_Nwp_Idx,Lat_Nwp_Idxx)
! z4 = data(Lon_Nwp_Idxx,Lat_Nwp_Idxx)
!---------------------------------------------------------------------------
 function INTERPOLATE_NWP_R4(z1,z2,z3,z4,lonx,latx) result(z)
  real, intent(in):: z1
  real, intent(in):: z2
  real, intent(in):: z3
  real, intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  real:: z
  !--- linear inteprpolation scheme
  if (minval((/z1,z2,z3,z4/)) == Missing_Value_Real4) then
    z = z1
  else
    z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
         (lonx) * ((1.0-latx) * z2 + (latx)* z4)
  endif
 end function INTERPOLATE_NWP_R4
 function INTERPOLATE_NWP_I4(z1,z2,z3,z4,lonx,latx) result(z)
  integer(kind=int4), intent(in):: z1
  integer(kind=int4), intent(in):: z2
  integer(kind=int4), intent(in):: z3
  integer(kind=int4), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  integer(kind=int4):: z
  !--- linear inteprpolation scheme
  if (minval((/z1,z2,z3,z4/)) == Missing_Value_Int4) then
    z = z1
  else
    z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
         (lonx) * ((1.0-latx) * z2 + (latx)* z4)
  endif
 end function INTERPOLATE_NWP_I4
 function INTERPOLATE_NWP_I2(z1,z2,z3,z4,lonx,latx) result(z)
  integer(kind=int2), intent(in):: z1
  integer(kind=int2), intent(in):: z2
  integer(kind=int2), intent(in):: z3
  integer(kind=int2), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  integer(kind=int2):: z
  !--- linear inteprpolation scheme
  if (minval((/z1,z2,z3,z4/)) == Missing_Value_Int2) then
    z = z1
  else
    z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
         (lonx) * ((1.0-latx) * z2 + (latx)* z4)
  endif
 end function INTERPOLATE_NWP_I2
 function INTERPOLATE_NWP_I1(z1,z2,z3,z4,lonx,latx) result(z)
  integer(kind=int1), intent(in):: z1
  integer(kind=int1), intent(in):: z2
  integer(kind=int1), intent(in):: z3
  integer(kind=int1), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  integer(kind=int1):: z
  !--- linear inteprpolation scheme
  if (minval((/z1,z2,z3,z4/)) == Missing_Value_Int1) then
    z = z1
  else
    z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
         (lonx) * ((1.0-latx) * z2 + (latx)* z4)
  endif
 end function INTERPOLATE_NWP_I1

!---------------------------------------------------------------------------
! subroutine INTERPOLATE_PROFILE(z1,z2,z3,z4,lonx,latx,z)
!
! description of arguments
! Lon_Nwp_Idx, Lat_Nwp_Idx - nwp indices of closest nwp cell
! Lon_Nwp_Idxx,Lat_Nwp_Idxx - nwp indices of nwp cells of diagonal of bounding box
! lonx - longitude weighting factor 
! latx = latitude weighting factor
! z1 = data(Lon_Nwp_Idx, Lat_Nwp_Idx)
! z2 = data(Lon_Nwp_Idxx,Lat_Nwp_Idx)
! z3 = data(Lon_Nwp_Idx,Lat_Nwp_Idxx)
! z4 = data(Lon_Nwp_Idxx,Lat_Nwp_Idxx)
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

 subroutine INTERPOLATE_PROFILE(z1,z2,z3,z4,lonx,latx,z)

  real, dimension(:), intent(in):: z1
  real, dimension(:), intent(in):: z2
  real, dimension(:), intent(in):: z3
  real, dimension(:), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  real, dimension(:), intent(out):: z

  !--- linear inteprpolation scheme
  z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
           (lonx) * ((1.0-latx) * z2 + (latx)* z4)

 end subroutine INTERPOLATE_PROFILE


!--------------------------------------------------------------------------
! subroutine INTERPOLATE_NWP_TZ_PROFILES(Elem_Idx,Line_Idx)
!
! spatially interpolate TZ profiles
!
!--------------------------------------------------------------------------
 subroutine INTERPOLATE_NWP_TZ_PROFILES(Elem_Idx,Line_Idx)

   integer, intent(in):: Elem_Idx,Line_Idx
   integer:: Ilev



   do Ilev = 1,NWP%Nlevels

      T_Prof(Ilev) = INTERPOLATE_NWP(NWP%T_Prof(Ilev,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                                             NWP%T_Prof(Ilev,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                                             NWP%T_Prof(Ilev,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                                             NWP%T_Prof(Ilev,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                                             NWP_PIX%Lon_Nwp_fac(Elem_Idx,Line_Idx), NWP_PIX%Lat_Nwp_fac(Elem_Idx,Line_Idx))

      Z_Prof(Ilev) = INTERPOLATE_NWP(NWP%Z_Prof(Ilev,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                                             NWP%Z_Prof(Ilev,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                                             NWP%Z_Prof(Ilev,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                                             NWP%Z_Prof(Ilev,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                                             NWP_PIX%Lon_Nwp_fac(Elem_Idx,Line_Idx), NWP_PIX%Lat_Nwp_fac(Elem_Idx,Line_Idx))
   enddo


 end subroutine INTERPOLATE_NWP_TZ_PROFILES

!-----------------------------------------------------------------------------
! subroutine CREATE_NWP_ARRAYS()
!
! allocate and initialize the memory needed for the nwp arrays
!-----------------------------------------------------------------------------
  subroutine CREATE_NWP_ARRAYS()

    !   Allocate arrays

    allocate(NWP%P_Std(NWP%Nlevels))
    allocate(NWP%Z_Prof(NWP%Nlevels, NWP%Nlon, NWP%Nlat))
    allocate(NWP%T_Prof(NWP%Nlevels, NWP%Nlon, NWP%Nlat))
    allocate(NWP%Ozone_Prof(NWP%Nlevels, NWP%Nlon, NWP%Nlat))
    allocate(NWP%Rh_Prof(NWP%Nlevels, NWP%Nlon, NWP%Nlat))
    allocate(NWP%Tpw_Prof(NWP%Nlevels, NWP%Nlon, NWP%Nlat))
    allocate(NWP%Clwmr_Prof(NWP%Nlevels, NWP%Nlon, NWP%Nlat))
    allocate(NWP%Icmr_Prof(NWP%Nlevels, NWP%Nlon, NWP%Nlat))
    allocate(NWP%U_Wnd_Prof(NWP%Nlevels, NWP%Nlon, NWP%Nlat))
    allocate(NWP%V_Wnd_Prof(NWP%Nlevels, NWP%Nlon, NWP%Nlat))


    allocate(NWP%Bad_Nwp_Mask(NWP%Nlon, NWP%Nlat))

    allocate(NWP%Pmsl(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Psfc(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Zsfc(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Tmpsfc(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Tmpair(NWP%Nlon, NWP%Nlat))
    allocate(NWP%T_Trop(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Z_Trop(NWP%Nlon, NWP%Nlat))
    allocate(NWP%P_Trop(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Rhsfc(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Rh300(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Uth(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Tpw(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Ozone(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Weasd(NWP%Nlon, NWP%Nlat))
    allocate(NWP%U_Wnd_10m(NWP%Nlon, NWP%Nlat))
    allocate(NWP%V_Wnd_10m(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Wnd_Spd_10m(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Wnd_Dir_10m(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Div_Sfc(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Div_200(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Lat(NWP%Nlat))
    allocate(NWP%Lon(NWP%Nlon))
    allocate(NWP%Freezing_lev_hgt(NWP%Nlon, NWP%Nlat))
    allocate(NWP%CAPE(NWP%Nlon, NWP%Nlat))

    allocate(NWP%Land(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Sea_Ice_Frac(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Sfc_Level(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Tropo_Level(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Inversion_Strength(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Inversion_Top(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Inversion_Base(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Inversion_Level(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Inversion_Level_Profile(NWP%Nlevels,NWP%Nlon, NWP%Nlat))
    allocate(NWP%Lifting_Condensation_Level_Height(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Convective_Condensation_Level_Height(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Level_Free_Convection_Height(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Equilibrium_Level_Height(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Freezing_Level_Height(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Freezing_Level_Pressure(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Homogenous_Freezing_Level_Height(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Homogenous_Freezing_Level_Pressure(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Upper_Limit_Water_Height(NWP%Nlon, NWP%Nlat))
    allocate(NWP%K_Index(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Level850(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Level700(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Level500(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Level300(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Level200(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Level100(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Pc(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Tc(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Sc_Lwp(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Lwp(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Iwp(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Cwp(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Cloud_Fraction_Satellite(NWP%Nlon, NWP%Nlat))
    allocate(NWP%High_Cloud_Fraction_Satellite(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Mid_Cloud_Fraction_Satellite(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Low_Cloud_Fraction_Satellite(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Ncld_Layers(NWP%Nlon, NWP%Nlat))
    allocate(NWP%Cld_Type(NWP%Nlon, NWP%Nlat))

    allocate(temp1d_Nwp(NWP%Nlevels))
    allocate(temp2d_Nwp_1(NWP%Nlon, NWP%Nlat))
    allocate(temp2d_Nwp_2(NWP%Nlon, NWP%Nlat))
    allocate(temp3d_Nwp_1(NWP%Nlevels, NWP%Nlon, NWP%Nlat))
    allocate(temp3d_Nwp_2(NWP%Nlevels, NWP%Nlon, NWP%Nlat))
    allocate(Tmpsfc_Nwp_Before(NWP%Nlon, NWP%Nlat))
    allocate(Tmpsfc_Nwp_After(NWP%Nlon, NWP%Nlat))

    allocate(T_Prof(NWP%Nlevels))
    allocate(Z_Prof(NWP%Nlevels))

    if (REFORMAT_GFS_ZXY == 1) then
     allocate(temp3d(NWP%Nlon, NWP%Nlat,NWP%Nlevels))
    else
     allocate(temp3d(NWP%Nlevels,NWP%Nlon, NWP%Nlat))
    endif
end subroutine CREATE_NWP_ARRAYS
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine INITIALIZE_NWP_ARRAYS

    integer:: Lon_Nwp_Idx, Lat_Nwp_Idx

    NWP%Pmsl = Missing_Value_Real4
    NWP%Psfc = Missing_Value_Real4
    NWP%Zsfc = Missing_Value_Real4
    NWP%Tmpsfc = Missing_Value_Real4
    NWP%Tmpair = Missing_Value_Real4
    NWP%T_Trop = Missing_Value_Real4
    NWP%Z_Trop = Missing_Value_Real4
    NWP%P_Trop = Missing_Value_Real4
    NWP%Rhsfc = Missing_Value_Real4
    NWP%Rh300 = Missing_Value_Real4
    NWP%Uth = Missing_Value_Real4
    NWP%Tpw = Missing_Value_Real4
    NWP%Ozone = Missing_Value_Real4
    NWP%Weasd = Missing_Value_Real4
    NWP%U_Wnd_10m = Missing_Value_Real4
    NWP%V_Wnd_10m = Missing_Value_Real4
    NWP%Wnd_Spd_10m = Missing_Value_Real4
    NWP%Wnd_Dir_10m = Missing_Value_Real4
    NWP%Div_Sfc = Missing_Value_Real4
    NWP%Div_200 = Missing_Value_Real4
    NWP%Lon = Missing_Value_Real4
    NWP%Lat = Missing_Value_Real4
    NWP%Freezing_lev_hgt = Missing_Value_Real4
    NWP%CAPE = Missing_Value_Real4

    NWP%Bad_Nwp_Mask = sym%NO
    NWP%P_Std = Missing_Value_Real4
    NWP%Z_Prof = Missing_Value_Real4
    NWP%T_Prof = Missing_Value_Real4
    NWP%Rh_Prof = Missing_Value_Real4
    NWP%Ozone_Prof = Missing_Value_Real4
    NWP%Tpw_Prof = Missing_Value_Real4
    NWP%Clwmr_Prof = Missing_Value_Real4
    NWP%ICmr_Prof = Missing_Value_Real4
    NWP%U_Wnd_Prof = Missing_Value_Real4
    NWP%V_Wnd_Prof = Missing_Value_Real4

    NWP%Land = Missing_Value_Int1
    NWP%Sea_Ice_Frac = Missing_Value_Real4
    NWP%Sfc_Level = Missing_Value_Int1
    NWP%Tropo_Level = Missing_Value_Int1
    NWP%Level850 = Missing_Value_Int1
    NWP%Level700 = Missing_Value_Int1
    NWP%Level500 = Missing_Value_Int1
    NWP%Level300 = Missing_Value_Int1
    NWP%Level200 = Missing_Value_Int1
    NWP%Level100 = Missing_Value_Int1
    NWP%Inversion_Top = Missing_Value_Real4
    NWP%Inversion_Base = Missing_Value_Real4
    NWP%Inversion_Strength = Missing_Value_Real4
    NWP%Inversion_Level =  Missing_Value_Int1
    NWP%Inversion_Level_Profile = Missing_Value_Int1
    NWP%Lifting_Condensation_Level_Height = 0
    NWP%Convective_Condensation_Level_Height = 0
    NWP%Level_Free_Convection_Height = 0
    NWP%Equilibrium_Level_Height = 0
    NWP%Freezing_Level_Height = Missing_Value_Real4
    NWP%Freezing_Level_Height = Missing_Value_Real4
    NWP%Homogenous_Freezing_Level_Height = Missing_Value_Real4
    NWP%Homogenous_Freezing_Level_Pressure = Missing_Value_Real4
    NWP%Upper_Limit_Water_Height = Missing_Value_Real4
    NWP%K_Index = Missing_Value_Real4
    NWP%Pc = Missing_Value_Real4
    NWP%Tc = Missing_Value_Real4
    NWP%Sc_Lwp = Missing_Value_Real4
    NWP%Iwp = Missing_Value_Real4
    NWP%Lwp = Missing_Value_Real4
    NWP%Cwp = Missing_Value_Real4
    NWP%Cloud_Fraction_Satellite = Missing_Value_Real4
    NWP%High_Cloud_Fraction_Satellite = Missing_Value_Real4
    NWP%Mid_Cloud_Fraction_Satellite = Missing_Value_Real4
    NWP%Low_Cloud_Fraction_Satellite = Missing_Value_Real4
    NWP%Ncld_Layers = Missing_Value_Int1
    NWP%Cld_Type = Missing_Value_Int1

    T_Prof = Missing_Value_Real4
    Z_Prof = Missing_Value_Real4
    Tmpsfc_Nwp_Before = Missing_Value_Real4
    Tmpsfc_Nwp_After = Missing_Value_Real4
    temp3d_Nwp_1 = Missing_Value_Real4
    temp3d_Nwp_2 = Missing_Value_Real4
    temp3d = Missing_Value_Real4

    !--- create nwp lat and lon vectors
    do Lat_Nwp_Idx = 1, NWP%Nlat
       NWP%Lat(Lat_Nwp_Idx) = lat1_Nwp + (Lat_Nwp_Idx-1) * dLat_Nwp
    end do

    do Lon_Nwp_Idx = 1, NWP%Nlon
       NWP%Lon(Lon_Nwp_Idx) = lon1_Nwp + (Lon_Nwp_Idx-1) * dLon_Nwp
       if (NWP%Lon(Lon_Nwp_Idx) > 180.0) then
         NWP%Lon(Lon_Nwp_Idx) = NWP%Lon(Lon_Nwp_Idx) - 360.0
       endif
    end do

end subroutine INITIALIZE_NWP_ARRAYS


!-----------------------------------------------------------------------------
! subroutine DESTROY_NWP_ARRAYS
!
! deallocate the memory needed for the nwp arrays
!-----------------------------------------------------------------------------
subroutine DESTROY_NWP_ARRAYS

    if (allocated(NWP%Land))           deallocate(NWP%Land)
    if (allocated(NWP%Sea_Ice_Frac))   deallocate(NWP%Sea_Ice_Frac)
    if (allocated(NWP%Sfc_Level))      deallocate(NWP%Sfc_Level)
    if (allocated(NWP%Tropo_Level))    deallocate(NWP%Tropo_Level)
    if (allocated(NWP%Level850))       deallocate(NWP%Level850)
    if (allocated(NWP%Level700))       deallocate(NWP%Level700)
    if (allocated(NWP%Level500))       deallocate(NWP%Level500)
    if (allocated(NWP%Level300))       deallocate(NWP%Level300)
    if (allocated(NWP%Level200))       deallocate(NWP%Level200)
    if (allocated(NWP%Level100))       deallocate(NWP%Level100)
    if (allocated(NWP%Inversion_Strength)) deallocate(NWP%Inversion_Strength)
    if (allocated(NWP%Inversion_Top)) deallocate(NWP%Inversion_Top)
    if (allocated(NWP%Inversion_Base)) deallocate(NWP%Inversion_Base)
    if (allocated(NWP%Inversion_Level)) deallocate(NWP%Inversion_Level)
    if (allocated(NWP%Inversion_Level_Profile)) deallocate(NWP%Inversion_Level_Profile)
    if (allocated(NWP%Lifting_Condensation_Level_Height))    deallocate(NWP%Lifting_Condensation_Level_Height)
    if (allocated(NWP%Convective_Condensation_Level_Height))    deallocate(NWP%Convective_Condensation_Level_Height)
    if (allocated(NWP%Level_Free_Convection_Height))    deallocate(NWP%Level_Free_Convection_Height)
    if (allocated(NWP%Equilibrium_Level_Height))    deallocate(NWP%Equilibrium_Level_Height)
    if (allocated(NWP%Freezing_Level_Height))    deallocate(NWP%Freezing_Level_Height)
    if (allocated(NWP%Freezing_Level_Pressure))  deallocate(NWP%Freezing_Level_Pressure)
    if (allocated(NWP%Homogenous_Freezing_Level_Height))  deallocate(NWP%Homogenous_Freezing_Level_Height)
    if (allocated(NWP%Homogenous_Freezing_Level_Pressure))  deallocate(NWP%Homogenous_Freezing_Level_Pressure)
    if (allocated(NWP%Upper_Limit_Water_Height)) deallocate(NWP%Upper_Limit_Water_Height)
    if (allocated(NWP%K_Index))       deallocate(NWP%K_Index)
    if (allocated(NWP%Pmsl))          deallocate(NWP%Pmsl)
    if (allocated(NWP%Psfc))          deallocate(NWP%Psfc)
    if (allocated(NWP%Zsfc))          deallocate(NWP%Zsfc)
    if (allocated(NWP%Tmpsfc))        deallocate(NWP%Tmpsfc)
    if (allocated(NWP%Tmpair))        deallocate(NWP%Tmpair)
    if (allocated(NWP%T_Trop))        deallocate(NWP%T_Trop)
    if (allocated(NWP%Z_Trop))        deallocate(NWP%Z_Trop)
    if (allocated(NWP%P_Trop))        deallocate(NWP%P_Trop)
    if (allocated(NWP%Rhsfc))         deallocate(NWP%Rhsfc)
    if (allocated(NWP%Rh300))         deallocate(NWP%Rh300)
    if (allocated(NWP%Uth))           deallocate(NWP%Uth)
    if (allocated(NWP%Tpw))           deallocate(NWP%Tpw)
    if (allocated(NWP%Ozone))         deallocate(NWP%Ozone)
    if (allocated(NWP%Weasd))         deallocate(NWP%Weasd)
    if (allocated(NWP%U_Wnd_10m))     deallocate(NWP%U_Wnd_10m)
    if (allocated(NWP%V_Wnd_10m))     deallocate(NWP%V_Wnd_10m)
    if (allocated(NWP%Wnd_Spd_10m))   deallocate(NWP%Wnd_Spd_10m)
    if (allocated(NWP%Wnd_Dir_10m))   deallocate(NWP%Wnd_Dir_10m)
    if (allocated(NWP%Div_Sfc))       deallocate(NWP%Div_Sfc)
    if (allocated(NWP%Div_200))       deallocate(NWP%Div_200)
    if (allocated(NWP%Lat))           deallocate(NWP%Lat)
    if (allocated(NWP%Lon))           deallocate(NWP%Lon)
    if (allocated(NWP%Freezing_lev_hgt)) deallocate(NWP%Freezing_lev_hgt)
    if (allocated(NWP%CAPE))          deallocate(NWP%CAPE)

    if (allocated(NWP%Bad_Nwp_Mask))  deallocate(NWP%Bad_Nwp_mask)
    if (allocated(NWP%P_Std))         deallocate(NWP%P_Std)
    if (allocated(NWP%Z_Prof))        deallocate(NWP%Z_Prof)
    if (allocated(NWP%T_Prof))        deallocate(NWP%T_Prof)
    if (allocated(NWP%Rh_Prof))       deallocate(NWP%Rh_Prof)
    if (allocated(NWP%Ozone_Prof))    deallocate(NWP%Ozone_Prof)
    if (allocated(NWP%Tpw_Prof))      deallocate(NWP%Tpw_Prof)
    if (allocated(NWP%Clwmr_Prof))    deallocate(NWP%Clwmr_Prof)
    if (allocated(NWP%Icmr_Prof))    deallocate(NWP%Icmr_Prof)
    if (allocated(NWP%U_Wnd_Prof))    deallocate(NWP%U_Wnd_Prof)
    if (allocated(NWP%V_Wnd_Prof))    deallocate(NWP%V_Wnd_Prof)

    if (allocated(NWP%Pc))            deallocate(NWP%Pc)
    if (allocated(NWP%Tc))            deallocate(NWP%Tc)
    if (allocated(NWP%Sc_Lwp))        deallocate(NWP%Sc_Lwp)
    if (allocated(NWP%Lwp))           deallocate(NWP%Lwp)
    if (allocated(NWP%Iwp))           deallocate(NWP%Iwp)
    if (allocated(NWP%Cwp))           deallocate(NWP%Cwp)
    if (allocated(NWP%Cloud_Fraction_Satellite))       deallocate(NWP%Cloud_Fraction_Satellite)
    if (allocated(NWP%High_Cloud_Fraction_Satellite))  deallocate(NWP%High_Cloud_Fraction_Satellite)
    if (allocated(NWP%Mid_Cloud_Fraction_Satellite))   deallocate(NWP%Mid_Cloud_Fraction_Satellite)
    if (allocated(NWP%Low_Cloud_Fraction_Satellite))   deallocate(NWP%Low_Cloud_Fraction_Satellite)
    if (allocated(NWP%Ncld_Layers))   deallocate(NWP%Ncld_Layers)
    if (allocated(NWP%Cld_Type))   deallocate(NWP%Cld_Type)

    if (allocated(temp1d_Nwp))        deallocate(temp1d_Nwp)
    if (allocated(temp2d_Nwp_1))      deallocate(temp2d_Nwp_1)
    if (allocated(temp2d_Nwp_2))      deallocate(temp2d_Nwp_2)
    if (allocated(temp3d_Nwp_1))      deallocate(temp3d_Nwp_1)
    if (allocated(temp3d_Nwp_2))      deallocate(temp3d_Nwp_2)
    if (allocated(temp3d))            deallocate(temp3d)
    if (allocated(T_Prof))            deallocate(T_Prof)
    if (allocated(Z_Prof))            deallocate(Z_Prof)
    if (allocated(Tmpsfc_Nwp_Before)) deallocate(Tmpsfc_Nwp_Before)
    if (allocated(Tmpsfc_Nwp_After))  deallocate(Tmpsfc_Nwp_After)

  end subroutine DESTROY_NWP_ARRAYS

!----------------------------------------------------------------------
! subroutine KNOWING_P_COMPUTE_T_Z_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,Ilev)
!
! derive height and temperature from a profile knowing pressure
!----------------------------------------------------------------------
 subroutine KNOWING_P_COMPUTE_T_Z_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,Ilev)
  integer, intent(in):: Lon_Nwp_Idx, Lat_Nwp_Idx
  real, intent(in):: P
  real, intent(out):: T,Z
  integer, intent(out):: Ilev
  real:: dp, dt, dz


  !--- interpolate pressure profile
  call LOCATE(NWP%P_Std,NWP%Nlevels,P,Ilev)
  Ilev = max(1,min(NWP%Nlevels-1,Ilev))

  dp = NWP%P_Std(Ilev+1) - NWP%P_Std(Ilev)
  dt = NWP%T_Prof(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%T_Prof(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)
  dz = NWP%Z_Prof(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%Z_Prof(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)

  !--- perform interpolation
  if (dp /= 0.0) then
   T = NWP%T_Prof(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) + dt/dp * (P - NWP%P_Std(Ilev))
   Z = NWP%Z_Prof(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) + dz/dp * (P - NWP%P_Std(Ilev))
  else
   T = NWP%T_Prof(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) 
   Z = NWP%Z_Prof(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx)
  endif

 end subroutine KNOWING_P_COMPUTE_T_Z_NWP

!----------------------------------------------------------------------
! subroutine KNOWING_Z_COMPUTE_T_P_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,Ilev)
!
! derive pressure and temperature from a profile knowing height
!----------------------------------------------------------------------
 subroutine KNOWING_Z_COMPUTE_T_P_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,Ilev)
  integer, intent(in):: Lon_Nwp_Idx, Lat_Nwp_Idx
  real, intent(in):: Z
  real, intent(out):: T,P
  integer, intent(out):: Ilev
  real:: dp
  real:: dt
  real:: dz
  integer:: kstart
  integer:: kend
  integer:: Nlevels_temp

  !---
  kstart = NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)
  kend = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)
  Nlevels_temp = kend - kstart + 1

  !--- interpolate pressure profile
  call LOCATE(NWP%Z_Prof(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx),Nlevels_temp,Z,Ilev)
  Ilev = max(1,min(Nlevels_temp-1,Ilev))
  Ilev = Ilev + kstart - 1

  dp = NWP%P_Std(Ilev+1) - NWP%P_Std(Ilev)
  dt = NWP%T_Prof(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%T_Prof(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)
  dz = NWP%Z_Prof(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%Z_Prof(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)

  !--- perform interpolation
  if (dp /= 0.0) then
   T = NWP%T_Prof(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) + dt/dz * (Z - NWP%Z_Prof(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx))
   P = NWP%P_Std(Ilev) + dp/dz * (Z - NWP%Z_Prof(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx))
  else
   T = NWP%T_Prof(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) 
   P = NWP%P_Std(Ilev)
  endif

 end subroutine KNOWING_Z_COMPUTE_T_P_NWP


!----------------------------------------------------------------------
! subroutine KNOWING_Z_COMPUTE_T_P_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,Ilev)
!
! derive pressure and temperature from a profile knowing height
!----------------------------------------------------------------------
 subroutine KNOWING_Z_COMPUTE_T_P_NWP_ARBITRARY_LEVELS(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z, &
                                                       kstart,kend,Ilev)
  integer, intent(in):: Lon_Nwp_Idx, Lat_Nwp_Idx
  real, intent(in):: Z
  integer, intent(in):: kstart
  integer, intent(in):: kend
  real, intent(out):: T,P
  integer, intent(out):: Ilev
  real:: dp
  real:: dt
  real:: dz
  integer:: Nlevels_temp

   Nlevels_temp = kend - kstart + 1

!--- interpolate pressure profile
  call LOCATE(NWP%Z_Prof(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx),Nlevels_temp,Z,Ilev)
  Ilev = max(1,min(Nlevels_temp-1,Ilev))
  Ilev = Ilev + kstart - 1

  dp = NWP%P_Std(Ilev+1) - NWP%P_Std(Ilev)
  dt = NWP%T_Prof(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%T_Prof(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)
  dz = NWP%Z_Prof(Ilev+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%Z_Prof(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx)

!--- perform interpolation
  if (dp /= 0.0) then
   T = NWP%T_Prof(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx) + dt/dz * (Z - NWP%Z_Prof(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx))
   P = NWP%P_Std(Ilev) + dp/dz * (Z - NWP%Z_Prof(Ilev,Lon_Nwp_Idx,Lat_Nwp_Idx))
  else
   T = NWP%T_Prof(Ilev, Lon_Nwp_Idx, Lat_Nwp_Idx)
   P = NWP%P_Std(Ilev)
  endif

 end subroutine KNOWING_Z_COMPUTE_T_P_NWP_ARBITRARY_LEVELS
!----------------------------------------------------------------------------
! subroutine KNOWING_T_COMPUTE_P_Z_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,klev,i,j,ierr)
!
! Interpolate P and Z from a profile knowing T
!
! klev = result is between klev and llev + 1
! 
! note, the highest level allowed is the tropopause
!
!----------------------------------------------------------------------------
 subroutine KNOWING_T_COMPUTE_P_Z_NWP(Lon_Nwp_Idx,Lat_Nwp_Idx,P,T,Z,klev,Elem_Idx,Line_Idx,ierr,z_interp_weight)

  integer, intent(in):: Lon_Nwp_Idx, Lat_Nwp_Idx
  real, intent(in):: T
  integer, intent(in):: Elem_Idx,Line_Idx
  real, intent(out):: P,Z
  integer, intent(out):: klev,ierr
  integer:: Nlevels_temp,kstart,kend,Level_within_Inversion_flag
  real:: dp, dt, dz,z1,z2,t1,t2
  real,optional,intent(out):: z_interp_weight

  ierr = sym%YES
  Z = Missing_Value_Real4
  P = Missing_Value_Real4
  klev = Missing_Value_Int4

  if (Lon_Nwp_Idx < 0) return
  if (Lon_Nwp_Idx > NWP%Nlon) return
  if (Lat_Nwp_Idx < 0) return
  if (Lat_Nwp_Idx > NWP%Nlat) return

  ierr = sym%NO

  !--- test for existence of a valid solution with troposphere
  kstart = NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)
  kend = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)
  Nlevels_temp = kend - kstart + 1

  !--- interpolate temperature profile

  !--- check to see if warmer than max, than assume at surface
  if (T > maxval(NWP%T_Prof(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx))) then
    P = NWP%P_Std(kend)
    Z = NWP%Z_Prof(kend,Lon_Nwp_Idx,Lat_Nwp_Idx)
    klev = kend - 1
    ierr = sym%NO
       if (present(z_interp_weight)) z_interp_weight = 0.
    return
  endif

  !--- check to see if colder than min, than assume tropopause
  if (T < minval(NWP%T_Prof(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx))) then
    P = NWP%P_Std(kstart)
    Z = NWP%Z_Prof(kstart,Lon_Nwp_Idx,Lat_Nwp_Idx)
    klev = kstart + 1
    ierr = sym%NO
       if (present(z_interp_weight)) z_interp_weight = 0.
    return
  endif
  
  !--- if there is an Inversion, look below first
  Level_within_Inversion_flag = 0
  if (NWP%Inversion_Level(Lon_Nwp_Idx, Lat_Nwp_Idx) > 0) then
   kstart = NWP%Inversion_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)
   kend = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)
   Nlevels_temp = kend - kstart + 1
   call LOCATE(NWP%T_Prof(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx),Nlevels_temp,T,klev)
   if ((klev > 0) .and. (klev < Nlevels_temp -1)) then
     klev = klev + kstart - 1
     Level_within_Inversion_flag = 1
   endif
  endif

  !--- if no solution within an Inversion, look above
  if (Level_within_Inversion_flag == 0) then

    kstart = NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)
    kend = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)
    Nlevels_temp = kend - kstart + 1
    call LOCATE(NWP%T_Prof(kstart:kend,Lon_Nwp_Idx,Lat_Nwp_Idx),Nlevels_temp,T,klev)
    klev = klev + kstart - 1
    klev = max(1,min(NWP%Nlevels-1,klev))

  endif

  !-- if solution is above trop, set to trop values
  if (klev < NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)) then

    P = NWP%P_Std(NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx))
    Z = NWP%Z_Prof(NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)

  else
    !--- determine derivatives
    dp = NWP%P_Std(klev+1) - NWP%P_Std(klev)

    if (NWP_PIX%Smooth_nwp_flag == sym%NO) then
     z2 = NWP%Z_Prof(klev+1,Lon_Nwp_Idx,Lat_Nwp_Idx)
     z1 = NWP%Z_Prof(klev,Lon_Nwp_Idx,Lat_Nwp_Idx)
     t2 = NWP%T_Prof(klev+1,Lon_Nwp_Idx,Lat_Nwp_Idx)
     t1 = NWP%T_Prof(klev,Lon_Nwp_Idx,Lat_Nwp_Idx)
    else
     z1 = INTERPOLATE_NWP(NWP%Z_Prof(klev,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                          NWP%Z_Prof(klev,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                          NWP%Z_Prof(klev,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                          NWP%Z_Prof(klev,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                          NWP_PIX%Lon_Nwp_fac(Elem_Idx,Line_Idx), NWP_PIX%Lat_Nwp_fac(Elem_Idx,Line_Idx))
     z2 = INTERPOLATE_NWP(NWP%Z_Prof(klev+1,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                          NWP%Z_Prof(klev+1,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                          NWP%Z_Prof(klev+1,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                          NWP%Z_Prof(klev+1,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                          NWP_PIX%Lon_Nwp_fac(Elem_Idx,Line_Idx), NWP_PIX%Lat_Nwp_fac(Elem_Idx,Line_Idx))
     t1 = INTERPOLATE_NWP(NWP%T_Prof(klev,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                          NWP%T_Prof(klev,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                          NWP%T_Prof(klev,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                          NWP%T_Prof(klev,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                          NWP_PIX%Lon_Nwp_fac(Elem_Idx,Line_Idx), NWP_PIX%Lat_Nwp_fac(Elem_Idx,Line_Idx))
     t2 = INTERPOLATE_NWP(NWP%T_Prof(klev+1,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                          NWP%T_Prof(klev+1,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)), &
                          NWP%T_Prof(klev+1,NWP_PIX%i_Nwp(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                          NWP%T_Prof(klev+1,NWP_PIX%i_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%j_Nwp_x(Elem_Idx,Line_Idx)), &
                          NWP_PIX%Lon_Nwp_fac(Elem_Idx,Line_Idx), NWP_PIX%Lat_Nwp_fac(Elem_Idx,Line_Idx))

    endif

    dz = z2 - z1
    dt = t2 - t1

    !--- perform interpolation
    if (dt /= 0.0) then
     P = NWP%P_Std(klev) + dp/dt * (T - t1)
     Z = z1 + dz/dt * (T - t1)
    else
     P = NWP%P_Std(klev) 
     Z = z1
    endif
 

     if (present(z_interp_weight)) z_interp_weight =  (z - z1)/dz

  !--- end above trop check
  endif

 end subroutine KNOWING_T_COMPUTE_P_Z_NWP

!-------------------------------------------------------------
! subroutine FIND_NWP_LEVELS(Lon_Nwp_Idx,Lat_Nwp_Idx)
!
! subroutine to find key levels in the profiles
!-------------------------------------------------------------
subroutine FIND_NWP_LEVELS(Lon_Nwp_Idx,Lat_Nwp_Idx)

 integer (kind=int4), intent(in):: Lon_Nwp_Idx, Lat_Nwp_Idx
 integer (kind=int4):: k
 real:: P_Trop_temp
 real:: dT_dZ
 real:: dT_dP
 real:: es
 real:: e
 real:: T
 real:: Td_Sfc
 real:: T850
 real:: T700
 real:: T500
 real:: Td850
 real:: Td700
 integer:: Top_Lev_Idx
 integer:: Base_Lev_Idx

   !--------------------------------------------------------------------
   !--- find surface level (standard closest but less than sfc pressure)
   !--------------------------------------------------------------------
   NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Nlevels
   do k = NWP%Nlevels,1,-1
       if (NWP%P_Std(k) < NWP%Psfc(Lon_Nwp_Idx,Lat_Nwp_Idx)) then
            NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
            exit
       endif
    enddo

   !--------------------------------------------------------------------
   !--- find some standard levels 
   !--------------------------------------------------------------------
   NWP%Level850(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0
   NWP%Level700(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0
   NWP%Level500(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0
   NWP%Level300(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0
   NWP%Level200(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0
   NWP%Level100(Lon_Nwp_Idx,Lat_Nwp_Idx) = 0
   do k = 1, NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)-1
      if ((NWP%P_Std(k) <= 850.0) .and. (NWP%P_Std(k+1) > 850.0)) then
        NWP%Level850(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
      if ((NWP%P_Std(k) <= 700.0) .and. (NWP%P_Std(k+1) > 700.0)) then
        NWP%Level700(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
      if ((NWP%P_Std(k) <= 500.0) .and. (NWP%P_Std(k+1) > 500.0)) then
        NWP%Level500(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
      if ((NWP%P_Std(k) <= 300.0) .and. (NWP%P_Std(k+1) > 300.0)) then
        NWP%Level300(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
      if ((NWP%P_Std(k) <= 200.0) .and. (NWP%P_Std(k+1) > 200.0)) then
        NWP%Level200(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
      if ((NWP%P_Std(k) <= 100.0) .and. (NWP%P_Std(k+1) > 100.0)) then
        NWP%Level100(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
   enddo

   !--------------------------------------------------------------------
   !--- find tropopause level  based on tropopause pressure
   !--- tropopause is between tropopause_level and tropopaue_level + 1
   !--------------------------------------------------------------------

   !--- constrain tropopause pressure to be greater than 75 mb
   P_Trop_Temp = max(NWP%P_Trop(Lon_Nwp_Idx,Lat_Nwp_Idx),75.0)

   NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) = 1
   do k = 1, NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)-1
      if ((NWP%P_Std(k) <= P_Trop_temp) .and. &
          (NWP%P_Std(k+1) > P_Trop_temp)) then
        NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) = k
      endif
   enddo

    !--- check if tropopause level found
    if (NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) == 0) then
        NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) = 1   !assume top level if no trop found
    endif         

    !--- Necessary as NWP has been showing trop index issues in the INDOEX region in nearly
    !--- the same grid point location. A trop index of 1 causes issues in downstream
    !--- subroutines and will cause Met-8 processing to die.
    if (NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) == 1) then
      NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) = 2
    endif         

    !--- store tropause height so when interpolate to pixel level
    NWP%Z_Trop(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Z_Prof(NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)

!--------------------------------------------------------------------
!--- find stratopause level starting at 500 mb by looking for levels
!--- with the minimum temp
!--- Note, no point in doing this until higher vertical resolution
!--- profiles are used.
!--------------------------------------------------------------------

   !---------------------------------------------------------------------
   ! Inversion Level Profile
   !---------------------------------------------------------------------
   NWP%Inversion_Level_Profile(:,Lon_Nwp_Idx,Lat_Nwp_Idx) = sym%NO

   do k = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx), NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx), -1

      if (NWP%P_Std(k) >= P_Inversion_Min) then

         if (NWP%T_Prof(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) > Delta_T_Inversion) then
            NWP%Inversion_Level_Profile(k-1:k,Lon_Nwp_Idx,Lat_Nwp_Idx) = sym%YES
         endif

      endif
   enddo

   Top_Lev_Idx =  0
   do k = NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)
      if (NWP%Inversion_Level_Profile(k,Lon_Nwp_Idx,Lat_Nwp_Idx) == sym%YES .and. Top_Lev_Idx == 0) then
         Top_Lev_Idx = k
         NWP%Inversion_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) = Top_Lev_Idx
         exit
      endif
   enddo

   Base_Lev_Idx = 0
   do k = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),-1
      if (NWP%Inversion_Level_Profile(k,Lon_Nwp_Idx,Lat_Nwp_Idx) == sym%YES .and. Base_Lev_Idx == 0) then
         Base_Lev_Idx = k
         exit
      endif
   enddo

   NWP%Inversion_Strength(Lon_Nwp_Idx,Lat_Nwp_Idx)  = Missing_Value_Real4
   NWP%Inversion_Base(Lon_Nwp_Idx,Lat_Nwp_Idx) = Missing_Value_Real4
   NWP%Inversion_Top(Lon_Nwp_Idx,Lat_Nwp_Idx) = Missing_Value_Real4

   !---- inversion top height (meters)
   if (Top_Lev_Idx /= 0)  NWP%Inversion_Top(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Z_Prof(Top_Lev_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx)

   !---- inversion base height (meters)
   if (Base_Lev_Idx /= 0) NWP%Inversion_Base(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Z_Prof(Base_Lev_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx)
   !--- assume inversion streches to surface if lowest level is the surface level
   if (Base_Lev_Idx == NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) .and. NWP%Zsfc(Lon_Nwp_Idx,Lat_Nwp_Idx) /= Missing_Value_Real4) then
    NWP%Inversion_Base(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Zsfc(Lon_Nwp_Idx,Lat_Nwp_Idx)
   endif

   !--- inversion temperature strength
   if (Base_Lev_Idx /= 0 .and. Top_Lev_Idx /= 0) then 
     NWP%Inversion_Strength(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%T_Prof(Top_Lev_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx) -  &
                                                   NWP%T_Prof(Base_Lev_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx)
     !--- assume inversion streches to surface if lowest level is the surface level
     if (Base_Lev_Idx == NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) .and.  NWP%Tmpair(Lon_Nwp_Idx,Lat_Nwp_Idx) /= Missing_Value_Real4) then
        NWP%Inversion_Strength(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%T_Prof(Top_Lev_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx) -  &
                                                          NWP%Tmpair(Lon_Nwp_Idx,Lat_Nwp_Idx)
     endif
   endif
  
   !-------------------------------------------------------------------
   ! Find the Height of the Freezing Level in the NWP Profiles
   ! Start at the Surface and work up
   !-------------------------------------------------------------------

   NWP%Freezing_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Z_Prof(NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)

   NWP%Freezing_Level_Pressure(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%P_Std(NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx))


  ! use nwp freezing level if available
  if (any(NWP%Freezing_lev_hgt /= Missing_Value_Real4)) then

     do k = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),-1

        if (NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) >= NWP%Freezing_lev_hgt(Lon_Nwp_Idx,Lat_Nwp_Idx)) then
             NWP%Freezing_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx)
             NWP%Freezing_Level_Pressure(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%P_Std(k)
             exit
        endif

     enddo

  else

   do k = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),-1

      if (NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) < 273.15) then

       if (k < NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)) then
          if (NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) < NWP%Z_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx)) then

             dT_dZ = (NWP%T_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx)) /  &
                     (NWP%Z_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx))

             NWP%Freezing_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = (NWP%T_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - 273.15) /  &
                                                                  dT_dZ + NWP%Z_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx)

             dT_dP = (NWP%T_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx)) /  &
                     (NWP%P_Std(k+1) - NWP%P_Std(k))

             NWP%Freezing_Level_Pressure(Lon_Nwp_Idx,Lat_Nwp_Idx) = (NWP%T_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - 273.15) /  &
                                                                    dT_dP + NWP%P_Std(k+1)
          else

             NWP%Freezing_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx)

             NWP%Freezing_Level_Pressure(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%P_Std(k)

          endif
       endif

       exit

      endif

   enddo

  endif

  ! setting freezing level to missing if surface temp is colder than 273.15K
  !if (NWP%Tmpsfc(Lon_Nwp_Idx,Lat_Nwp_Idx) < 273.15) then
  !    NWP%Freezing_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = Missing_Value_Real4
  !    NWP%Freezing_Level_Pressure(Lon_Nwp_Idx,Lat_Nwp_Idx) = Missing_Value_Real4
  !endif

   !-------------------------------------------------------------------
   ! Find the Height of the Homogenous Freezing Level in the NWP Profiles
   ! Start at the surface and work up
   !-------------------------------------------------------------------

   NWP%Homogenous_Freezing_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Z_Prof(NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)

   NWP%Homogenous_Freezing_Level_Pressure(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%P_Std(NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx))

   do k = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),-1

      if (NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) <= 233.00) then

       if (k < NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)) then
          if (NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) < NWP%Z_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx)) then

             dT_dZ = (NWP%T_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx)) /  &
                     (NWP%Z_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx))

             NWP%Homogenous_Freezing_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = (NWP%T_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - 233.0) /  &
                                                                             dT_dZ + NWP%Z_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx)

             dT_dP = (NWP%T_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx)) /  &
                     (NWP%P_Std(k+1) - NWP%P_Std(k))

             NWP%Homogenous_Freezing_Level_Pressure(Lon_Nwp_Idx,Lat_Nwp_Idx) = (NWP%T_Prof(k+1,Lon_Nwp_Idx,Lat_Nwp_Idx) - 233.0) /  &
                                                                               dT_dP + NWP%P_Std(k+1)
          else

             NWP%Homogenous_Freezing_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx)

             NWP%Homogenous_Freezing_Level_Pressure(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%P_Std(k)

          endif
       endif

       exit

      endif
   enddo

   !---------------------------------------------------------------------------
   ! Find the Height of the Lifting(LCL) and Convective Condensation Level(CCL) 
   ! Heights from the NWP Profiles
   !----------------------------------------------------------------------------
    NWP%Lifting_Condensation_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = Missing_Value_Real4
    NWP%Convective_Condensation_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = Missing_Value_Real4
    T = NWP%Tmpair(Lon_Nwp_Idx,Lat_Nwp_Idx)    ! K
    Td_Sfc = MISSING_VALUE_REAL4

    if (T > 180.0) then
     if (T > 253.0) then
       es = VAPOR(T)    !saturation vapor pressure wrt water hpa
     else
       es = VAPOR_ICE(T) !saturation vapor pressure wrt ice hpa
     endif
     e = es * NWP%Rhsfc(Lon_Nwp_Idx,Lat_Nwp_Idx) / 100.0  !vapor pressure in hPa
     Td_Sfc = 273.15 + 243.5 * alog(e / 6.112)  / (17.67 - alog(e/6.112))      !Dewpoint T in K
     NWP%Lifting_Condensation_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = 1000.0 * 0.125*(T - Td_Sfc)  ! meters 
     NWP%Convective_Condensation_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = 1000.0 * (T - Td_Sfc)/4.4  ! meters (AW 2015/04/30)
   endif

   !---------------------------------------------------------------------------
   ! Cmpute Level Free Convection and Equilibrium
   !----------------------------------------------------------------------------
   call COMPUTE_LFC_EL_HEIGHT(NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx), &
                              NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx), &
                              NWP%Nlevels, &
                              NWP%P_Std,  &
                              NWP%T_Prof(:,Lon_Nwp_Idx,Lat_Nwp_Idx), &
                              Td_Sfc, &
                              NWP%Z_Prof(:,Lon_Nwp_Idx,Lat_Nwp_Idx), &
                              NWP%Level_Free_Convection_Height(Lon_Nwp_Idx,Lat_Nwp_Idx), &
                              NWP%Equilibrium_Level_Height(Lon_Nwp_Idx,Lat_Nwp_Idx))

   !---------------------------------------------------------------------------
   ! K Index 
   !----------------------------------------------------------------------------
   NWP%K_Index(Lon_Nwp_Idx,Lat_Nwp_Idx) = Missing_Value_Real4
   if (NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx) > 0 .and. &
       NWP%Level700(Lon_Nwp_Idx,Lat_Nwp_Idx) > 0 .and. &
       NWP%Level500(Lon_Nwp_Idx,Lat_Nwp_Idx) > 0) then

     if (NWP%Level850(Lon_Nwp_Idx,Lat_Nwp_Idx) == 0) then
         k = NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)
     else
         k = NWP%Level850(Lon_Nwp_Idx,Lat_Nwp_Idx)
     endif
     T850 = NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx)
     T700 = NWP%T_Prof(NWP%Level700(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)
     T500 = NWP%T_Prof(NWP%Level500(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)
     if (T850 > 253.0) then
       es = VAPOR(T850)    !saturation vapor pressure wrt water hpa
      else
       es = VAPOR_ICE(T850) !saturation vapor pressure wrt ice hpa
     endif
     e = es * NWP%Rh_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) / 100.0  !vapor pressure in hPa
     Td850 = 273.15 + 243.5 * alog(e / 6.112)  / (17.67 - alog(e/6.112))       !Dewpoint T in K

     if (T700 > 253.0) then
       es = VAPOR(T700)    !saturation vapor pressure wrt water hpa
     else
       es = VAPOR_ICE(T700) !saturation vapor pressure wrt ice hpa
     endif
     e = es * NWP%Rh_Prof(NWP%Level700(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx) / 100.0  !vapor pressure in hPa
     Td700 = 273.15 + 243.5 * alog(e / 6.112)  / (17.67 - alog(e/6.112))       !Dewpoint T in K
    
     if (T700 /= Td700) then
        NWP%K_Index(Lon_Nwp_Idx,Lat_Nwp_Idx) = (T850 - T500) + Td850 - (T700 - Td700) - 273.15
     endif
   endif

   !-------------------------------------------------------------------
   ! Find the Height above which to consider CWP to be 100% ice in the NWP Profiles
   ! Start at the Tropopause and work down
   !-------------------------------------------------------------------
   NWP%Upper_Limit_Water_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) &
      = NWP%Z_Prof(NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),Lon_Nwp_Idx,Lat_Nwp_Idx)

   do k = NWP%Tropo_Level(Lon_Nwp_Idx,Lat_Nwp_Idx),NWP%Sfc_Level(Lon_Nwp_Idx,Lat_Nwp_Idx)-1

      if (NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) > 260.0) then

       if (NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) > NWP%Z_Prof(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx)) then
         dT_dZ = (NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%T_Prof(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx)) /  &
                 (NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - NWP%Z_Prof(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx))
         NWP%Upper_Limit_Water_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = (NWP%T_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx) - 260.0) / &
                 dT_dZ + NWP%Z_Prof(k,Lon_Nwp_Idx,Lat_Nwp_Idx)
       else
         NWP%Upper_Limit_Water_Height(Lon_Nwp_Idx,Lat_Nwp_Idx) = NWP%Z_Prof(k-1,Lon_Nwp_Idx,Lat_Nwp_Idx)
       endif

       exit

      endif
   enddo

end subroutine FIND_NWP_LEVELS

!======================================================================
! subroutine COMPUTE_COAST_MASK_NWP(j1,j2)
!
! Compute Coast Mask for NWP data
!
! dependencies - Land_Mask and land_Mask should be already computed
!======================================================================
 subroutine COMPUTE_COAST_MASK_NWP(j1,j2)

  integer, intent(in):: j1, j2
  integer:: Elem_Idx,Line_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx

  !--- loop over all pixels in segment
  do Line_Idx = j1,j1+j2-1
     do Elem_Idx = 1, Image%Number_Of_Elements

        !--- initialize to no
        Sfc%Coast_Mask_Nwp(Elem_Idx,Line_Idx) = sym%NO
       
        !--- check for valid pixels
        if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES)  then
          cycle
        endif

        !--- save nwp indices
        Lon_Nwp_Idx = NWP_PIX%i_Nwp(Elem_Idx,Line_Idx)
        Lat_Nwp_Idx = NWP_PIX%j_Nwp(Elem_Idx,Line_Idx)

        !--- derive nwp coast mask
        if (Sfc%Land_Mask(Elem_Idx,Line_Idx) /= NWP%Land(Lon_Nwp_Idx,Lat_Nwp_Idx)) then
           Sfc%Coast_Mask_Nwp(Elem_Idx,Line_Idx) = sym%YES
        endif 

     !--- end loop over all pixels in segment
     enddo
   enddo

 end subroutine COMPUTE_COAST_MASK_NWP

!======================================================================
! pull off 300hpa RH value and store in array
! assumes there is a 300 hpa level
!======================================================================
function COMPUTE_Rh300(Lon_Idx,Lat_Idx,Pressure_Profile,Rel_Humid_Profile) result(Rh300)
 integer, intent(in):: Lon_Idx, Lat_Idx
 real, dimension(:), intent(in):: Pressure_Profile
 real, dimension(:), intent(in):: Rel_Humid_Profile
 real:: Rh300

 Rh300 = MISSING_VALUE_REAL4

 if (NWP%Level300(Lon_Idx,Lat_Idx) /= MISSING_VALUE_INT1) then
    Rh300 = Rel_Humid_Profile(NWP%Level300(Lon_Idx,Lat_Idx))
 endif

end function COMPUTE_Rh300

function COMPUTE_UTH(Lon_Idx,Lat_Idx,Pressure_Profile,Rel_Humid_Profile) result(Uth)
 integer, intent(in):: Lon_Idx, Lat_Idx
 real, dimension(:), intent(in):: Pressure_Profile
 real, dimension(:), intent(in):: Rel_Humid_Profile
 real:: Uth
 real:: Uth_Sum
 integer:: Level_Idx, Uth_Count

  Uth_Sum = 0.0 
  Uth_Count = 0
  do Level_Idx = NWP%Tropo_Level(Lon_Idx,Lat_Idx), NWP%Level300(Lon_Idx,Lat_Idx)
     Uth_Sum = Uth_Sum + Rel_Humid_Profile(Level_Idx) 
     Uth_Count = Uth_Count + 1
  enddo

  if (Uth_Count > 0) then
      Uth = Uth_Sum / Uth_Count
  else
      Uth = MISSING_VALUE_REAL4
  endif

end function COMPUTE_UTH
!======================================================================
! compute divergence from finite difference of wind field
! output units should be 1/second * 1.0e06
!======================================================================
function COMPUTE_DIVERGENCE(Lon_Idx,Lat_Idx,Level_Idx) result(Div)
 integer, intent(in):: Lon_Idx,Lat_Idx
 integer(kind=int1), intent(in):: Level_Idx
 real:: Delta_x, Delta_y
 real:: Div

 integer:: Lat0_Idx,Lat2_Idx,Lon0_Idx,Lon2_Idx

 Lat0_Idx = min(NWP%Nlat,max(1,Lat_Idx - 1))
 Lat2_Idx = min(NWP%Nlat,max(1,Lat_Idx + 1))
 Lon0_Idx = min(NWP%Nlon,max(1,Lon_Idx - 1))
 Lon2_Idx = min(NWP%Nlon,max(1,Lon_Idx + 1))

 Delta_X = (NWP%Lon(Lon2_Idx) - NWP%Lon(Lon0_Idx))*110.0*1000.0    !m
 Delta_y = (NWP%Lat(Lat2_Idx) - NWP%Lat(Lat0_Idx))*110.0*1000.0    !m
 Delta_X = Delta_X * cos(NWP%Lat(Lat_Idx)*DTOR)


 if ((NWP%U_Wnd_Prof(Level_Idx,Lon_Idx,Lat_Idx) == MISSING_VALUE_REAL4) .or.  &
     (NWP%V_Wnd_Prof(Level_Idx,Lon_Idx,Lat_Idx) == MISSING_VALUE_REAL4)) then

     Div = MISSING_VALUE_REAL4

 else
  
    Div = (NWP%U_Wnd_Prof(Level_Idx,Lon2_Idx,Lat_Idx) -          &     !  dU /
           NWP%U_Wnd_Prof(Level_Idx,Lon0_Idx,Lat_Idx))/Delta_x + &     !      dx
          (NWP%V_Wnd_Prof(Level_Idx,Lon_Idx,Lat2_Idx) -          &     !  dV /
           NWP%V_Wnd_Prof(Level_Idx,Lon_Idx,Lat0_Idx))/Delta_y         !      dy

    Div = Div * 1.0e06

 endif

end function COMPUTE_DIVERGENCE

!======================================================================
! subroutine COMPUTE_NWP_CLOUD_PARAMETERS(Clwmr_Profile,
!    temperature_Profile,Pressure_Profile,P_Cld_Top,iwp,lwp,cwp,
!    Number_of_Cloud_Layers)
!
! P_Cld_Top = Pressure of Cloud Top (tau = 1)
!======================================================================
 subroutine COMPUTE_NWP_CLOUD_PARAMETERS(Tropopause_Level_Nwp, &
                                         Surface_Level_Nwp, &
                                         Clwmr_Profile, &
                                         Rel_Humid_Profile, &
                                         Temperature_Profile, &
                                         Pressure_Profile, &
                                         P_Cld_Top,T_Cld_Top, &
                                         Iwp,Lwp,Sc_Lwp,Cwp,&
                                         Cloud_Fraction_Satellite, &
                                         High_Cloud_Fraction_Satellite, &
                                         Mid_Cloud_Fraction_Satellite, &
                                         Low_Cloud_Fraction_Satellite, &
                                         Number_Of_Cloud_Layers, &
                                         Cloud_Type)
   implicit none                                      

  integer(kind=int1), intent(in):: Tropopause_Level_Nwp
  integer(kind=int1), intent(in):: Surface_Level_Nwp
  real, dimension(:), intent(in):: Clwmr_Profile
  real, dimension(:), intent(in):: Rel_Humid_Profile
  real, dimension(:), intent(in):: Temperature_Profile
  real, dimension(:), intent(in):: Pressure_Profile
  real, intent(out):: P_Cld_Top
  real, intent(out):: T_Cld_Top
  real, intent(out):: Iwp
  real, intent(out):: Lwp
  real, intent(out):: Sc_Lwp
  real, intent(out):: Cwp
  integer (kind=int1), intent(out):: Number_Of_Cloud_Layers
  integer (kind=int1), intent(out):: Cloud_Type

  real, dimension(size(Clwmr_Profile)):: Cloud_Fraction_Profile
  real, intent(out):: Cloud_Fraction_Satellite
  real, intent(out):: High_Cloud_Fraction_Satellite
  real, intent(out):: Mid_Cloud_Fraction_Satellite
  real, intent(out):: Low_Cloud_Fraction_Satellite
 
  real:: Sat_Specific_Humidity
  real:: Clwmr_Min
  real, parameter:: Lwp_Threshold  = 5.0
  real, parameter:: Iwp_Threshold  = 1.0
  real, parameter:: Frac_Min_Threshold  = 1.0e-06
  real, parameter:: Clwmr_Min_Threshold  = 0.0
  real, parameter:: Cwp_Min_Threshold  = 5.0     !g/m^2 - arbitrary need to investigate
  real, parameter:: Max_Temperature_Ice = 273.15
  real, parameter:: Min_Temperature_Water  = 253.15
  real, parameter:: Max_Temperature_Sc_Water = 273.15
  real, parameter:: Min_Temperature_Sc_Water  = 263.15
  integer:: Number_Of_Levels_In_Profile
  real:: Factor
  real:: Clwmr_Ice_Layer
  real:: Clwmr_Water_Layer
  real:: Clwmr_Sc_Water_Layer
  real:: Ice_Frac_Top
  real:: Ice_Frac_Bot
  real:: Sc_Water_Frac_Top
  real:: Sc_Water_Frac_Bot
  integer:: Ilay
  integer:: Ilev
  real, parameter:: k1 = 0.25
  real, parameter:: k2 = 100
  real, parameter:: k3 = 0.49
  real:: T
  real:: e
  real:: es

  Number_Of_Levels_In_Profile = size(Clwmr_Profile)

  if (Tropopause_Level_Nwp <= 0 .or. Surface_Level_Nwp <= 0) then
     return
  endif

  !---------------------------------------------------------------------
  ! convert clwmr profiles into optical depths
  ! iwp and lwp are in g/m^2
  !---------------------------------------------------------------------
  Lwp = 0.0
  Iwp = 0.0
  Cwp = 0.0
  Sc_Lwp = 0.0
  Clwmr_Water_Layer = 0.0
  Clwmr_Ice_Layer = 0.0
  P_Cld_Top = Missing_Value_Real4
  T_Cld_Top = Missing_Value_Real4

     do Ilay = Tropopause_Level_Nwp-1 , Surface_Level_Nwp-1

        Ice_Frac_Top  = min(1.0,max(0.0,&
                           (Max_Temperature_Ice - Temperature_Profile(ilay))/ &
                           (Max_Temperature_Ice-Min_Temperature_Water)))

        Ice_Frac_Bot  = min(1.0,max(0.0,&
                           (Max_Temperature_Ice - Temperature_Profile(ilay+1))/ &
                           (Max_Temperature_Ice-Min_Temperature_Water)))

        Clwmr_Ice_Layer = 0.5 * (Ice_Frac_Top*Clwmr_Profile(ilay) +  &
                                 Ice_Frac_Bot*Clwmr_Profile(ilay+1))

        Clwmr_Water_Layer = 0.5 * ((1.0-Ice_Frac_Top)*Clwmr_Profile(ilay) + &
                                   (1.0-Ice_Frac_Bot)*Clwmr_Profile(ilay+1))

        Sc_Water_Frac_Top  = min(1.0,max(0.0,&
                           (Max_Temperature_Sc_Water - Temperature_Profile(ilay))/ &
                           (Max_Temperature_Sc_Water - Min_Temperature_Sc_Water)))

        Sc_Water_Frac_Bot  = min(1.0,max(0.0,&
                           (Max_Temperature_Sc_Water - Temperature_Profile(ilay+1))/ &
                           (Max_Temperature_Sc_Water - Min_Temperature_Sc_Water)))

        Clwmr_Sc_Water_Layer = 0.5 * (Sc_Water_Frac_Top*Clwmr_Profile(ilay) + &
                                      Sc_Water_Frac_Bot*Clwmr_Profile(ilay+1))

        Factor = 1000.0 * 100.0 * (Pressure_Profile(ilay+1) - Pressure_Profile(ilay)) / g

        Iwp = Iwp + Clwmr_Ice_Layer * Factor
        Lwp = Lwp + Clwmr_Water_Layer * Factor
        Sc_Lwp = Sc_Lwp + Clwmr_Sc_Water_Layer * Factor

        Cwp = Lwp + Iwp

        !--- compute P_Cld_Top
        if ((Cwp >= Iwp_Threshold) .and. (P_Cld_Top == Missing_Value_Real4)) then
            P_Cld_Top = NWP%P_Std(ilay+1)
            T_Cld_Top = Temperature_Profile(ilay+1)
        endif

      enddo

  !----------------------------------------------------------------------
  ! test to see if sufficient amount of cloud water was found in the 
  ! column to estimate cloud fractions.  if not, exit
  !----------------------------------------------------------------------
  if (Cwp < Cwp_Min_Threshold) then
     Number_Of_Cloud_Layers = 0
     Cloud_Type = sym%CLEAR_TYPE
     Cloud_Fraction_Satellite = 0.0
     High_Cloud_Fraction_Satellite = 0.0
     Mid_Cloud_Fraction_Satellite = 0.0
     Low_Cloud_Fraction_Satellite = 0.0
     return
  endif

  !----------------------------------------------------------------------
  ! compute number of distinct cloud layers in GFS profile
  !----------------------------------------------------------------------
  Number_Of_Cloud_Layers = 0

  if (P_Cld_Top /= Missing_Value_Real4) then

         do ilay = 3, Surface_Level_Nwp - 1

            Clwmr_Min = 1.0e-05 * Pressure_Profile(ilay)/1000.0    !

            if (Clwmr_Profile(ilay) > Clwmr_Min .and.   &
                Clwmr_Profile(ilay-1) < Clwmr_Min .and. &
                Clwmr_Profile(ilay-2) < Clwmr_Min) then

                Number_Of_Cloud_Layers = Number_Of_Cloud_Layers + 1

            endif
        enddo
  endif 

  !----------------------------------------------------------------------
  ! compute cloud type
  !----------------------------------------------------------------------
  Cloud_Type = sym%CLEAR_TYPE

  if (T_Cld_Top /= Missing_Value_Real4) then

      Cloud_Type = sym%CIRRUS_TYPE
      if (T_Cld_Top > 258.0) Cloud_Type = sym%SUPERCOOLED_TYPE
      if (T_Cld_Top > 273.0) Cloud_Type = sym%WATER_TYPE

      if (Cloud_Type == sym%CIRRUS_TYPE) then
        if (Number_Of_Cloud_Layers > 1) Cloud_Type = sym%OVERLAP_TYPE
        if (Cwp > 100) Cloud_Type = sym%OPAQUE_ICE_TYPE
      endif

  endif

  !----------------------------------------------------------------------
  ! compute cloud fraction profile
  ! procedure described at http://www.emc.ncep.noaa.gov/GFS/doc.php
  ! Note, Clwmr_Profile has units of g/g
  !----------------------------------------------------------------------
  Cloud_Fraction_Profile = 0
  do Ilev = Tropopause_Level_Nwp-1, Surface_Level_Nwp

    if (Clwmr_Profile(Ilev) > Clwmr_Min_Threshold) then

     T = Temperature_Profile(ilev)
     if (T > 180.0) then
      if (T > 253.0) then
        es = VAPOR(T)    !saturation vapor pressure wrt water hpa
      else
        es = VAPOR_ICE(T) !saturation vapor pressure wrt ice hpa
      endif
      e = es * Rel_Humid_Profile(ilev) / 100.0  !vapor pressure in hPa

      Sat_Specific_Humidity = (0.622 * es) / (Pressure_Profile(ilev) - es)
  
      Cloud_Fraction_Profile(Ilev) = ((Rel_Humid_Profile(ilev)/100.0)**k1) * (1.0 - &
                              exp( (-1.0*k2*1000.0*Clwmr_Profile(ilev)) /  &
                                   (((1.0-Rel_Humid_Profile(ilev)/100.0)*Sat_Specific_Humidity)**k3)))

      if (Rel_Humid_Profile(Ilev) > 100.0) Cloud_Fraction_Profile(Ilev) = 1.0

     endif

    endif

  enddo

  !----------------------------------------------------------------------
  ! compute cloud fraction as seen by a satellite using Max-Random Overlap
  !----------------------------------------------------------------------
  Cloud_Fraction_Satellite = 0.0
  High_Cloud_Fraction_Satellite = 0.0
  Mid_Cloud_Fraction_Satellite = 0.0
  Low_Cloud_Fraction_Satellite = 0.0

  do Ilev = Tropopause_Level_Nwp-1, Surface_Level_Nwp

       if (Cloud_Fraction_Profile(Ilev) > Frac_Min_Threshold)  then

          if (Cloud_Fraction_Satellite == 0.0) then

               Cloud_Fraction_Satellite = Cloud_Fraction_Profile(Ilev)

          else
!           if (Cloud_Fraction_Profile(Ilev-1) > Frac_Min_Threshold .and. &
!             Cloud_Fraction_Profile(Ilev) > Frac_Min_Threshold) then       !max overlap

!             Cloud_Fraction_Satellite_Temp = max(Cloud_Fraction_Profile(Ilev-1),Cloud_Fraction_Profile(Ilev))
!             Cloud_Fraction_Satellite = max(Cloud_Fraction_Satellite,Cloud_Fraction_Satellite_Temp)

!           else    !random overlap

!             Cloud_Fraction_Satellite_Temp = 1.0 - (1.0-Cloud_Fraction_Satellite)*(1.0-Cloud_Fraction_Profile(Ilev))

!           endif

            Cloud_Fraction_Satellite = max(Cloud_Fraction_Satellite, Cloud_Fraction_Profile(Ilev))

          endif

          if (Pressure_Profile(Ilev) <= 350.0) then
               High_Cloud_Fraction_Satellite = Cloud_Fraction_Satellite
          elseif (Pressure_Profile(Ilev) > 350.0 .and. Pressure_Profile(Ilev) <= 642.0) then
               Mid_Cloud_Fraction_Satellite = Cloud_Fraction_Satellite - High_Cloud_Fraction_Satellite
          else
               Low_Cloud_Fraction_Satellite = Cloud_Fraction_Satellite &
                     - High_Cloud_Fraction_Satellite - Mid_Cloud_Fraction_Satellite
          endif

          if (Cloud_Fraction_Satellite == 1.0) then
             exit
          endif
               
       endif

!write (unit=6,fmt="(I4,10f8.5)") Ilev, Cloud_Fraction_Profile(Ilev), Cloud_Fraction_Satellite

  enddo

 end subroutine COMPUTE_NWP_CLOUD_PARAMETERS

!----------------------------------------------------------------------
subroutine COMPUTE_SEGMENT_NWP_CLOUD_PARAMETERS()

  integer:: Elem_Idx 
  integer:: Line_Idx 
  integer:: Lon_Idx 
  integer:: Lat_Idx 
  logical,save :: at_least_one_outside = .false.

  !--- initialize global arrays which are output of this routine
  NWP%Cld_Type = Missing_Value_Int1
  NWP%Sc_Lwp = Missing_Value_Real4
  NWP%Lwp = Missing_Value_Real4
  NWP%Iwp = Missing_Value_Real4
  NWP%Cwp = Missing_Value_Real4
  NWP%Pc = Missing_Value_Real4
  NWP%Tc = Missing_Value_Real4
  NWP%Ncld_Layers = Missing_Value_Int1
  NWP%Cloud_Fraction_Satellite = Missing_Value_Real4
  NWP%High_Cloud_Fraction_Satellite = Missing_Value_Real4
  NWP%Mid_Cloud_Fraction_Satellite = Missing_Value_Real4
  NWP%Low_Cloud_Fraction_Satellite = Missing_Value_Real4

  !--- loop over pixels in segment
   line_loop: do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment
     element_loop: do Elem_Idx = 1, Image%Number_Of_Elements

      !--- check for bad pixels
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
        cycle
      endif

      !--- check for space views
      if (Geo%Space_Mask(Elem_Idx,Line_Idx) ) then
        cycle
      endif

      Lon_Idx = NWP_PIX%I_Nwp(Elem_Idx,Line_Idx)
      Lat_Idx = NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)

      if (Lon_Idx <= 0 .or. Lat_Idx <= 0 )   then
        if ( .not. at_least_one_outside ) then
          print *, "NWP Routine: At Least One Missing lat,lon ", Elem_Idx,Line_Idx,Lon_Idx, Lat_Idx
          at_least_one_outside = .true.
        end if  
        cycle
      endif

      !--- check to see if this cell was already done
      if (NWP%Cld_Type(Lon_Idx,Lat_Idx) >= 0) then
        cycle
      endif

      !--- compute cloud parameters from nwp profiles
      call COMPUTE_NWP_CLOUD_PARAMETERS(NWP%Tropo_Level(Lon_Idx,Lat_Idx), &
                                        NWP%Sfc_Level(Lon_Idx,Lat_Idx), &
                                        NWP%Clwmr_Prof(:,Lon_Idx,Lat_Idx), &
                                        NWP%Rh_Prof(:,Lon_Idx,Lat_Idx), &
                                        NWP%T_Prof(:,Lon_Idx,Lat_Idx), &
                                        NWP%P_Std, &
                                        NWP%Pc(Lon_Idx,Lat_Idx), &
                                        NWP%Tc(Lon_Idx,Lat_Idx), &
                                        NWP%Iwp(Lon_Idx,Lat_Idx), &
                                        NWP%Lwp(Lon_Idx,Lat_Idx), &
                                        NWP%Sc_Lwp(Lon_Idx,Lat_Idx), &
                                        NWP%Cwp(Lon_Idx,Lat_Idx), &
                                        NWP%Cloud_Fraction_Satellite(Lon_Idx,Lat_Idx), &
                                        NWP%High_Cloud_Fraction_Satellite(Lon_Idx,Lat_Idx), &
                                        NWP%Mid_Cloud_Fraction_Satellite(Lon_Idx,Lat_Idx), &
                                        NWP%Low_Cloud_Fraction_Satellite(Lon_Idx,Lat_Idx), &
                                        NWP%Ncld_Layers(Lon_Idx,Lat_Idx), &
                                        NWP%Cld_Type(Lon_Idx,Lat_Idx))

       !--- compute some RH metrics
       NWP%Uth(Lon_Idx,Lat_Idx) = COMPUTE_UTH(Lon_Idx,Lat_Idx,NWP%P_Std,NWP%Rh_Prof(:,Lon_Idx,Lat_Idx))
       NWP%Rh300(Lon_Idx,Lat_Idx) = COMPUTE_RH300(Lon_Idx,Lat_Idx,NWP%P_Std,NWP%Rh_Prof(:,Lon_Idx,Lat_Idx))
       NWP%Div_Sfc(Lon_Idx,Lat_Idx) = COMPUTE_DIVERGENCE(Lon_Idx,Lat_Idx,NWP%Sfc_Level(Lon_Idx,Lat_Idx))
       NWP%Div_200(Lon_Idx,Lat_Idx) = COMPUTE_DIVERGENCE(Lon_Idx,Lat_Idx,NWP%Level200(Lon_Idx,Lat_Idx))

     end do element_loop
  end do line_loop

end subroutine COMPUTE_SEGMENT_NWP_CLOUD_PARAMETERS

!----------------------------------------------------------------------
! subroutine TEMPORAL_INTERP_TMPSFC_NWP(start_itime, end_itime)
!----------------------------------------------------------------------
subroutine TEMPORAL_INTERP_TMPSFC_NWP(start_itime, end_itime)
integer(kind=int4), intent(in):: start_itime
integer(kind=int4), intent(in):: end_itime
real(kind=real4):: x_interp
real(kind=real4):: start_hour
real(kind=real4)::  mean_hours
real(kind=real4):: segment_start_time
real(kind=real4):: segment_end_time
real(kind=real4):: segment_mean_time
integer(kind=int4):: nwp_start_hour_segment
integer(kind=int4):: nwp_end_hour_segment

!--- check for valid times (assume positive)
if (start_itime < 0 .or. end_itime < 0) then
  return
endif

!--- determine the bound day values for this orbit
start_hour = start_itime/60.0/60.0/1000.0

nwp_start_hour_segment = 6*int(start_hour/6)
nwp_end_hour_segment = nwp_start_hour_segment + 6

if (nwp_start_hour == 0 .and. nwp_end_hour_segment == 24) then
   nwp_end_hour_segment = 0
endif
if (nwp_end_hour == 24 .and. nwp_start_hour_segment == 0) then
   nwp_start_hour_segment = 24
endif

!--- initialize interpolation weight
x_interp = Missing_Value_Real4

!--- 
if (nwp_end_hour_segment == nwp_start_hour) then
  x_interp = 0.0
elseif (nwp_start_hour_segment == nwp_end_hour) then
  x_interp = 1.0
endif

if ((nwp_start_hour_segment == nwp_start_hour) .and. & 
    (nwp_end_hour_segment == nwp_end_hour)) then

!----  determine mean time for this segment
segment_start_time = start_itime/86400000.0_real4
segment_end_time = end_itime/86400000.0_real4
segment_mean_time = 0.5*(segment_start_time + segment_end_time)
mean_hours = segment_mean_time*24.0

 !--- determine interpolation weight for 6 hrly ncep fields
 !
 x_interp = (mean_hours  - nwp_start_hour)/  &
            (nwp_end_hour - nwp_start_hour)
 x_interp = max(0.0,min(1.0,x_interp))

endif

 !--- apply interpolation
if (x_interp /= Missing_Value_Real4) then 
 where ( (Tmpsfc_Nwp_Before /= Missing_Value_Real4) .and. (Tmpsfc_Nwp_After /= Missing_Value_Real4) )
    NWP%Tmpsfc = (1.0 - x_interp) * Tmpsfc_Nwp_Before + x_interp * Tmpsfc_Nwp_After
 end where
endif

end subroutine TEMPORAL_INTERP_TMPSFC_NWP

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_R4(xNwp,xPix,interp_method_in,diag_flag)

    real, dimension(:,:), intent(in):: xNwp
    real, dimension(:,:), intent(out):: xPix
    integer, intent(in), optional:: interp_method_in 
    integer, intent(in), optional:: diag_flag
    integer:: interp_method
    integer:: Elem_Idx
    integer:: Line_Idx

    xPix = Missing_Value_Real4

    interp_method = 0

    if (present(interp_method_in)) interp_method = interp_method_in

    do Elem_Idx = 1, size(NWP_PIX%I_Nwp(:,1))
     do Line_Idx = 1, size(NWP_PIX%I_Nwp(1,:))

      if ((NWP_PIX%I_Nwp(Elem_Idx,Line_Idx) > 0).and.(NWP_PIX%J_Nwp(Elem_Idx,Line_Idx) > 0)) then

        !--- default is nearest neighbor
        xPix(Elem_Idx,Line_Idx) = xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx))

        if (interp_method /= 0) then

         !--- check for array indices
         if ((NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx) > 0) .and. (NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx) > 0)) then

            !--- if any are missing, skip interp
            if (xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)) == Missing_Value_Real4 .or. &
                xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)) == Missing_Value_Real4 .or. &
                xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)) == Missing_Value_Real4 .or. &
                xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)) == Missing_Value_Real4) then

                cycle   

            endif

            xPix(Elem_Idx,Line_Idx) = INTERPOLATE_NWP( &
                       xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)), &
                       NWP_PIX%Lon_Nwp_Fac(Elem_Idx,Line_Idx), &
                       NWP_PIX%Lat_Nwp_Fac(Elem_Idx,Line_Idx))

         endif
        endif
      endif
     enddo
    enddo

    !--- ensure missing values are changed to clavr-x default
    where(xPix == Missing_Nwp)
      xPix = Missing_Value_Real4
    endwhere 

end subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_R4
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I1(xNwp,xPix,interp_method_in,diag_flag)
    integer(kind=int1), dimension(:,:), intent(in):: xNwp
    integer(kind=int1), dimension(:,:), intent(out):: xPix
    integer, intent(in), optional:: interp_method_in 
    integer, intent(in), optional:: diag_flag
    integer:: interp_method
    integer:: Elem_Idx
    integer:: Line_Idx
    xPix = Missing_Value_Int1
    interp_method = 0
    if (present(interp_method_in)) interp_method = interp_method_in
    do Elem_Idx = 1, size(NWP_PIX%I_Nwp(:,1))
     do Line_Idx = 1, size(NWP_PIX%I_Nwp(1,:))
      if ((NWP_PIX%I_Nwp(Elem_Idx,Line_Idx) > 0).and.(NWP_PIX%J_Nwp(Elem_Idx,Line_Idx) > 0)) then
        xPix(Elem_Idx,Line_Idx) = xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx))
        if (interp_method /= 0) then
         if ((NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx) > 0).and.(NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx) > 0)) then

            !--- if any are missing, skip interp
            if (xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)) == Missing_Value_Int1 .or. &
                xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)) == Missing_Value_Int1 .or. &
                xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)) == Missing_Value_Int1 .or. &
                xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)) == Missing_Value_Int1) then

                cycle   

            endif

            xPix(Elem_Idx,Line_Idx) = INTERPOLATE_NWP( &
                       xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)), &
                       NWP_PIX%Lon_Nwp_Fac(Elem_Idx,Line_Idx), &
                       NWP_PIX%Lat_Nwp_Fac(Elem_Idx,Line_Idx))
         endif
        endif
      endif
     enddo
    enddo
end subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I1

subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I2(xNwp,xPix,interp_method_in,diag_flag)
    integer(kind=int2), dimension(:,:), intent(in):: xNwp
    integer(kind=int2), dimension(:,:), intent(out):: xPix
    integer, intent(in), optional:: interp_method_in 
    integer, intent(in), optional:: diag_flag
    integer:: interp_method
    integer:: Elem_Idx
    integer:: Line_Idx
    xPix = Missing_Value_Int2
    interp_method = 0
    if (present(interp_method_in)) interp_method = interp_method_in
    do Elem_Idx = 1, size(NWP_PIX%I_Nwp(:,1))
     do Line_Idx = 1, size(NWP_PIX%I_Nwp(1,:))
      if ((NWP_PIX%I_Nwp(Elem_Idx,Line_Idx) > 0).and.(NWP_PIX%J_Nwp(Elem_Idx,Line_Idx) > 0)) then
        xPix(Elem_Idx,Line_Idx) = xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx))
        if (interp_method /= 0) then
         if ((NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx) > 0).and.(NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx) > 0)) then

            !--- if any are missing, skip interp
            if (xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)) == Missing_Value_Int2 .or. &
                xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)) == Missing_Value_Int2 .or. &
                xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)) == Missing_Value_Int2 .or. &
                xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)) == Missing_Value_Int2) then

                cycle   

            endif

            xPix(Elem_Idx,Line_Idx) = INTERPOLATE_NWP( &
                       xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)), &
                       NWP_PIX%Lon_Nwp_Fac(Elem_Idx,Line_Idx), &
                       NWP_PIX%Lat_Nwp_Fac(Elem_Idx,Line_Idx))
         endif
        endif
      endif
     enddo
    enddo
end subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I2

subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I4(xNwp,xPix,interp_method_in,diag_flag)
    integer(kind=int4), dimension(:,:), intent(in):: xNwp
    integer(kind=int4), dimension(:,:), intent(out):: xPix
    integer, intent(in), optional:: interp_method_in 
    integer, intent(in), optional:: diag_flag
    integer:: interp_method
    integer:: Elem_Idx
    integer:: Line_Idx
    xPix = Missing_Value_Int4
    interp_method = 0
    if (present(interp_method_in)) interp_method = interp_method_in
    do Elem_Idx = 1, size(NWP_PIX%I_Nwp(:,1))
     do Line_Idx = 1, size(NWP_PIX%I_Nwp(1,:))
      if ((NWP_PIX%I_Nwp(Elem_Idx,Line_Idx) > 0).and.(NWP_PIX%J_Nwp(Elem_Idx,Line_Idx) > 0)) then
        xPix(Elem_Idx,Line_Idx) = xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx))
        if (interp_method /= 0) then
         if ((NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx) > 0).and.(NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx) > 0)) then

            !--- if any are missing, skip interp
            if (xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)) == Missing_Value_Int4 .or. &
                xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)) == Missing_Value_Int4 .or. &
                xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)) == Missing_Value_Int4 .or. &
                xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)) == Missing_Value_Int4) then

                cycle   

            endif

            xPix(Elem_Idx,Line_Idx) = INTERPOLATE_NWP( &
                       xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)), &
                       xNwp(NWP_PIX%I_Nwp_x(Elem_Idx,Line_Idx),NWP_PIX%J_Nwp_x(Elem_Idx,Line_Idx)), &
                       NWP_PIX%Lon_Nwp_Fac(Elem_Idx,Line_Idx), &
                       NWP_PIX%Lat_Nwp_Fac(Elem_Idx,Line_Idx))
         endif
        endif
      endif
     enddo
    enddo
end subroutine CONVERT_NWP_ARRAY_TO_PIXEL_ARRAY_I4


end module NWP_COMMON_MOD
