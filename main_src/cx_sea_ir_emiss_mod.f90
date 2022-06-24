!$Id:$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: cx_sea_ir_emiss_mod.f90 (src)
!       CX_SEA_IR_EMISS_MOD (program)
!
! PURPOSE: Routine to handle the University of Edinburgh Sea Emissivity
!
! DESCRIPTION: see below
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
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
!  July 2018 - created
!
! Link: https://datashare.is.ed.ac.uk/handle/10283/17
!
! README FROM DATA PROVIDER
!
!Calculation of emissivity of pure water and seawater, 600-3350 cm-1
!
!M. J. Filipiak, C. J. Merchant and O. Embury
!
!International Journal of Remote Sensing, Vol. X, No. X, Month 2008, xxx-xxx
!
!Refractive indices (n, k) in the wavenumber range 500 - 3500 cm-1 are in
!
!nk274_287_300_purewater.txt
!nk274_287_300_seawater.txt
!
!with 3-sigma uncertainties in
!
!nk3Sigma.txt
!
!Emissivities for the ATSR/ATSR view angles are in
!
!ARCNadirEmissivityPureWater.nc
!ARCNadirEmissivitySeawater.nc
!ARCForwardEmissivityPureWater.nc
!ARCForwardEmissivitySeawater.nc
!
!and for angles 0 - 90 at coarse resolution in
!
!ARCWideangleEmissivityPureWater.nc
!ARCWideangleEmissivitySeawater.nc
!
!tabulated as follows (use ncdump -c <file name> to view the netCDF information)
!
!nadir view_angle = 0, 8, 16, 19, 21, 22, 23 ; (degrees)
!forward view_angle = 51, 51.5, 52, 52.5, 53, 53.5, 54, 54.5, 55, 55.5, 56, 56.5, 57 ; (degrees)
!wideangle view_angle = 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85 ; (degrees)
!
!wavenumber = 600, 700, 740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 
!   840, 850, 860, 870, 880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 
!   980, 990, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 
!   2000, 2100, 2200, 2300, 2360, 2370, 2380, 2390, 2400, 2410, 2420, 2430, 
!   2440, 2450, 2460, 2470, 2480, 2490, 2500, 2510, 2520, 2530, 2540, 2550, 
!   2560, 2570, 2580, 2590, 2600, 2610, 2620, 2630, 2640, 2650, 2660, 2670, 
!   2680, 2690, 2700, 2710, 2720, 2730, 2740, 2750, 2760, 2770, 2780, 2790, 
!   2800, 2810, 2820, 2830, 2840, 2850, 2860, 2870, 2880, 2890, 2900, 2910, 
!   2920, 2930, 2940, 2950, 2960, 2970, 2980, 2990, 3000, 3010, 3020, 3030, 
!   3040, 3050, 3060, 3070, 3080, 3090, 3100, 3110, 3120, 3130, 3140, 3240, 
!   3350 ; (cm-1)
!wind_speed = 0, 1, 3, 5, 10, 15, 20, 25 ; (m s-1)
!temperature = 270, 280, 290, 300, 310 ; (K)
!
!There are no error estimates for the emissivity data.
!--------------------------------------------------------------------------------------
module CX_SEA_IR_EMISS_MOD
  use CALIBRATION_CONSTANTS_MOD, only: Planck_Nu
  use CONSTANTS_MOD
  use PIXEL_COMMON_MOD, only: Image, Ch, Geo, Sensor, Sfc, NWP_PIX, Ancil_Data_Dir, emiss_sea_option, Bad_Pixel_Mask
  use FILE_TOOLS,only:
  use CX_NETCDF4_MOD
  use CLAVRX_MESSAGE_MOD

  implicit none

  private

  public::  FORGET_SEA_IR_EMISS, GET_SEGMENT_SEA_IR_EMISS
  private:: GET_PIXEL_SEA_IR_EMISS,READ_SEA_IR_EMISS

  !--- 
  integer, private, save:: Num_Wvn
  integer, private, save:: Num_Wind 
  integer, private, save:: Num_Temp 
  integer, private, save:: Num_Angle
  integer, private, save:: Num_Wndspd
  real, private, save, allocatable, dimension(:):: Temp_Table
  real, private, save, allocatable, dimension(:):: Angle_Table
  real, private, save, allocatable, dimension(:):: Wvn_Table
  real, private, save, allocatable, dimension(:):: Wndspd_Table
  real, private, save, allocatable, dimension(:,:,:,:):: Sea_Ir_Emiss_Table

  character(len=33), parameter:: Sea_Ir_Emiss_Table_Name = 'ARCWideangleEmissivitySeaWater.nc'
 
! integer, public, save:: Use_Sea_Ir_Emiss = 0

  real, private, save:: Wvn_Prev = MISSING_VALUE_REAL4
  real, private, save:: Angle_Prev = MISSING_VALUE_REAL4
  real, private, save:: Sst_Prev = MISSING_VALUE_REAL4
  real, private, save:: Wndspd_Prev = MISSING_VALUE_REAL4

  integer, private, save:: Wvn_Idx_Prev = MISSING_VALUE_INT4
  integer, private, save:: Angle_Idx_Prev = MISSING_VALUE_INT4
  integer, private, save:: Sst_Idx_Prev = MISSING_VALUE_INT4
  integer, private, save:: Wndspd_Idx_Prev = MISSING_VALUE_INT4

  real, private, parameter:: Wvn_Delta = 1.0
  real, private, parameter:: Angle_Delta = 1.0
  real, private, parameter:: Sst_Delta = 1.0
  real, private, parameter:: Wndspd_Delta = 1.0

  !-- spacing of angles in table
  real, private, parameter:: Angle_Table_Delta = 5.0
  
  
  type init_status_type
    logical :: is_set = .false.
  end type init_status_type
  
  type( init_status_type), save :: init_status

  contains
  !----------------------------------------------------------------------
  !  main public routine to compute for a segemnt
  !----------------------------------------------------------------------
  subroutine  GET_SEGMENT_SEA_IR_EMISS()

     integer:: Elem_Idx, Line_Idx, Chan_Idx

     do Elem_Idx = 1, Image%Number_Of_Elements 

        do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment

           if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle 

           if (Sfc%Land(Elem_Idx,Line_Idx) /= sym%Land .and. Sfc%Snow(Elem_Idx,Line_Idx) == sym%NO_SNOW) then

           do Chan_Idx = 1, Nchan_Clavrx
 
              if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then

                 if (ch(Chan_Idx)%Obs_Type == THERMAL_OBS_TYPE .or. ch(Chan_Idx)%Obs_Type == MIXED_OBS_TYPE) then

                       ch(Chan_Idx)%Sfc_Emiss(Elem_Idx,Line_Idx) = &
                                                GET_PIXEL_SEA_IR_EMISS( &
                                                Planck_Nu(Chan_Idx), &
                                                NWP_PIX%Wnd_Spd_10m(Elem_Idx,Line_Idx), &
                                                Geo%Satzen(Elem_Idx,Line_Idx), &
                                                NWP_PIX%Tsfc(Elem_Idx,Line_Idx))
                 endif

              endif

           enddo

           endif

        enddo

     enddo

  end subroutine  GET_SEGMENT_SEA_IR_EMISS

  !----------------------------------------------------------------------
  ! Allocate the Arrays and Read the Data
  !----------------------------------------------------------------------
  subroutine READ_SEA_IR_EMISS()
    
   
    
    
   integer, dimension(1) :: Sds_Start_1d, Sds_Edge_1d, Sds_Stride_1d
   integer, dimension(4) :: Sds_Start_4d, Sds_Edge_4d, Sds_Stride_4d
   integer:: Ncid_Sea_Emiss
   integer, dimension(1):: Dim_1d
   character(len=200):: File_Name_Full
   
    
   !------------------------------------------------------------------------------
   ! read in data
   !------------------------------------------------------------------------------
   emiss_sea_option = 1
   File_Name_Full = trim(Ancil_Data_Dir)//"/static/sfc_data/"//trim(Sea_Ir_Emiss_Table_Name)
   call OPEN_NETCDF(trim(File_Name_Full), Ncid_Sea_Emiss)
   if (Ncid_Sea_Emiss < 0) then 
      call MESG("Error opening Sea_IR_Emiss file, default will be used; file = " &
               //trim(File_Name_Full), level=Verb_Lev%WARNING)
      emiss_sea_option = 0
      return
   endif

   !---------------------------------------------------------------------------
   ! Read in Dimensions
   !---------------------------------------------------------------------------
   call READ_NETCDF_DIMENSION_1D(Ncid_Sea_Emiss, 'view_angle', Dim_1d)
   Num_Angle = Dim_1d(1)
   call READ_NETCDF_DIMENSION_1D(Ncid_Sea_Emiss, 'wind_speed', Dim_1d) 
   Num_Wndspd = Dim_1d(1)
   call READ_NETCDF_DIMENSION_1D(Ncid_Sea_Emiss, 'wavenumber', Dim_1d) 
   Num_Wvn = Dim_1d(1)
   call READ_NETCDF_DIMENSION_1D(Ncid_Sea_Emiss, 'temperature', Dim_1d) 
   Num_Temp = Dim_1d(1)
 
   !---------------------------------------------------------------------------
   ! Allocate Memory for Tables
   !---------------------------------------------------------------------------
   allocate(Wndspd_Table(Num_Wndspd))
   allocate(Wvn_Table(Num_Wvn))
   allocate(Angle_Table(Num_Angle))
   allocate(Temp_Table(Num_Temp))
   allocate(Sea_IR_Emiss_Table(Num_Angle, Num_Wvn, Num_Wndspd, Num_Temp))
   
   !---------------------------------------------------------------------------
   ! Read Tables
   !---------------------------------------------------------------------------
   sds_start_1d = 1
   sds_stride_1d = 1

   sds_edge_1d(1) = Num_Wndspd
   call READ_NETCDF(Ncid_Sea_Emiss,Sds_Start_1d,Sds_Stride_1d,Sds_Edge_1d,'wind_speed',Wndspd_Table)

   sds_edge_1d(1) = Num_Angle
   call READ_NETCDF(Ncid_Sea_Emiss,Sds_Start_1d,Sds_Stride_1d,Sds_Edge_1d,'view_angle',Angle_Table)

   sds_edge_1d(1) = Num_Wvn
   call READ_NETCDF(Ncid_Sea_Emiss,Sds_Start_1d,Sds_Stride_1d,Sds_Edge_1d,'wavenumber',Wvn_Table)

   sds_edge_1d(1) = Num_Temp
   call READ_NETCDF(Ncid_Sea_Emiss,Sds_Start_1d,Sds_Stride_1d,Sds_Edge_1d,'temperature',Temp_Table)

   sds_start_4d = 1
   sds_edge_4d(4) = Num_Temp
   sds_edge_4d(3) = Num_Wndspd
   sds_edge_4d(2) = Num_Wvn
   sds_edge_4d(1) = Num_Angle
   sds_stride_4d = 1

   call READ_NETCDF(Ncid_Sea_Emiss,Sds_Start_4d,Sds_Stride_4d,Sds_Edge_4d,'emissivity',Sea_IR_Emiss_Table)

   call CLOSE_NETCDF(Ncid_Sea_Emiss)

   call MESG("SEA_IR_EMISS read in successfully", level=Verb_Lev%DEFAULT)
   init_status % is_set = .true.

  end subroutine READ_SEA_IR_EMISS

!----------------------------------------------------------------------------------------
! Interpolate Emissivity for conditions of a pixel
!
!----------------------------------------------------------------------------------------
function GET_PIXEL_SEA_IR_EMISS(Wvn,Wndspd,Angle,Sst) Result(Sea_Ir_Emiss)

  real, intent(in):: Wvn,Wndspd,Angle,Sst
  real:: Sea_Ir_Emiss
 
  integer:: Wvn_Idx
  integer:: Angle_Idx
  integer:: Sst_Idx
  integer:: Wndspd_Idx
  integer, save:: Wvn_Idx_Prev = 0
  integer, save:: Angle_Idx_Prev = 0
  integer, save:: Sst_Idx_Prev = 0
  integer, save:: Wndspd_Idx_Prev = 0
  real, save:: Wvn_Prev = 0.0
  real, save:: Angle_Prev = 0.0
  real, save:: Sst_Prev = 0.0
  real, save:: Wndspd_Prev = 0.0
  real:: Sea_IR_Emiss_1, Sea_IR_Emiss_2
  
  
  if ( .not. init_status % is_set ) call READ_SEA_IR_EMISS()
  
  !--- initialize output
  Sea_IR_Emiss =  MISSING_VALUE_REAL4

  !--- find nearest indices
  if (abs(Sst - Sst_Prev) > Sst_Delta .or. Sst_Idx_Prev <= 0) then
   Sst_Idx = minloc(abs(Sst - Temp_Table),1)
   Sst_Idx = max(1,min(Num_Temp,Sst_Idx))
   Sst_Idx_Prev = Sst_Idx
   Sst_Prev = Sst
  else
   Sst_Idx = Sst_Idx_Prev
  endif
  
  if (abs(Wndspd - Wndspd_Prev) > Wndspd_Delta .or. Wndspd_Idx_Prev <= 0) then
   Wndspd_Idx = minloc(abs(Wndspd - Wndspd_Table),1)
   Wndspd_Idx = max(1,min(Num_Wndspd,Wndspd_Idx))
   Wndspd_Prev = Wndspd
   Wndspd_Idx_Prev = Wndspd_Idx
  else
   Wndspd_Idx = Wndspd_Idx_Prev
  endif

  if (abs(Wvn - Wvn_Prev) > Wvn_Delta .or. Wvn_Idx_Prev <= 0) then
   Wvn_Idx = minloc(abs(Wvn - Wvn_Table),1)
   Wvn_Idx = max(1,min(Num_Wvn,Wvn_Idx))
   Wvn_Prev = Wvn
   Wvn_Idx_Prev = Wvn_Idx
  else
   Wvn_Idx = Wvn_Idx_Prev
  endif

  !--- interpolate in angle (which is even spaced)
  !wideangle view_angle = 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85 ; (degrees)
  if (abs(Angle - Angle_Prev) > Angle_Delta .or. Angle_Idx_Prev <= 0) then
   !--- take closest angle
   !Angle_Idx = minloc(abs(Angle - Angle_Table),1)
   !Angle_Idx = max(1,min(Num_Angle,Angle_Idx))


   !--- this is index for Angle_Table <= Angle 
   !--- assumes Angles are always positive
   Angle_Idx = min(int(Angle/Angle_Table_Delta) + 1, Num_Angle - 1)
   Angle_Prev = Angle
   Angle_Idx_Prev = Angle_Idx
  else
   Angle_Idx = Angle_Idx_Prev
  endif




  !--- Select value from emissivity table without any interp
  !Sea_IR_Emiss = Sea_IR_Emiss_Table(Angle_Idx,Wvn_Idx,Wndspd_Idx,Sst_Idx)

  !--- Select value from emissivity table with angle interp
  Sea_IR_Emiss_1 = Sea_IR_Emiss_Table(Angle_Idx,Wvn_Idx,Wndspd_Idx,Sst_Idx)
  Sea_IR_Emiss_2 = Sea_IR_Emiss_Table(Angle_Idx+1,Wvn_Idx,Wndspd_Idx,Sst_Idx)

  Sea_IR_Emiss = Sea_IR_Emiss_1 + (Sea_IR_Emiss_2 - Sea_IR_Emiss_1) *  &
                (Angle - Angle_Table(Angle_Idx)) / Angle_Table_Delta
  
! print *, "Test ", Angle, Angle_Table(Angle_Idx), Sea_IR_Emiss_1, Sea_IR_Emiss_2, Sea_IR_Emiss

end function GET_PIXEL_SEA_IR_EMISS
!-----------------------------------------------------------------------------
!  Deallocate all allocated arrays (call after this data is no longer needed)
!-----------------------------------------------------------------------------
subroutine FORGET_SEA_IR_EMISS()

   if ( allocated (Wvn_Table) )     deallocate(Wvn_Table)
   if ( allocated (Wndspd_Table) )  deallocate(Wndspd_Table)
   if ( allocated (Angle_Table) )   deallocate(Angle_Table)
   if ( allocated (Temp_Table) )    deallocate(Temp_Table)
   if ( allocated (Sea_IR_Emiss_Table) )  deallocate(Sea_IR_Emiss_Table)

   emiss_sea_option = 0
   
   init_status % is_set = .false.

end subroutine FORGET_SEA_IR_EMISS

end module CX_SEA_IR_EMISS_MOD
