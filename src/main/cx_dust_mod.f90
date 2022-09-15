!$Id:$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: cx_dust_mod.f90 (src)
!       CX_DUST_MOD (program)
!
! PURPOSE:
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
!--------------------------------------------------------------------------------------
module CX_DUST_MOD
  use CONSTANTS_MOD
  use PIXEL_COMMON_MOD, only: Use_ABI_Dust, Image, Ch, Geo, Sensor, Sfc, Ancil_Data_Dir, &
                              Bad_Pixel_Mask, CLDMASK
  use CX_NETCDF4_MOD
  use CLAVRX_MESSAGE_MOD

  implicit none

  private

  public:: READ_ABI_DUST_LUT, FORGET_ABI_DUST, GET_SEGMENT_ABI_DUST_PROB
  private:: GET_PIXEL_ABI_DUST_PROB

  !--- 
  integer, private, save:: Nbins
  real, private, save:: Bin_Min
  real, private, save:: Bin_Width
  real, private, save:: Prior_Yes
  real, private, save, allocatable, dimension(:):: Btd_11_12_Cond_Ratio_Table
  real, private, save, allocatable, dimension(:):: Btd_10_11_Cond_Ratio_Table
  real, private, save, allocatable, dimension(:):: Btd_85_10_Cond_Ratio_Table
  real, private, save, allocatable, dimension(:):: OB_Btd_10_11_Cond_Ratio_Table
  real, private, save, allocatable, dimension(:):: OB_Btd_11_12_Cond_Ratio_Table
  real, private, save, allocatable, dimension(:):: OB_Btd_85_10_Cond_Ratio_Table

  character(len=33), parameter:: ABI_Dust_Prob_Table_Name = 'abi_ir_dust_prob_lut.nc'
 
  contains
  !----------------------------------------------------------------------
  !  main public routine to compute for a segemnt
  !----------------------------------------------------------------------
  subroutine  GET_SEGMENT_ABI_DUST_PROB()

     integer:: Elem_Idx, Line_Idx, Chan_Idx

     CLDMASK%Dust_Prob = MISSING_VALUE_REAL4
     CLDMASK%Dust_Mask = MISSING_VALUE_INT1

     do Elem_Idx = 1, Image%Number_Of_Elements 

        do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment

           if (Bad_Pixel_Mask(ELem_Idx,Line_Idx) == sym%YES) cycle

           !--- determine which pixels to process
           if ((Sfc%Land(ELem_Idx,Line_Idx) == sym%SHALLOW_OCEAN .or. &
                Sfc%Land(Elem_Idx,Line_Idx) == sym%MODERATE_OCEAN .or. &
                Sfc%Land(Elem_Idx,Line_Idx) == sym%DEEP_OCEAN) .and.  &
               Sfc%Snow(Elem_Idx,Line_Idx) == sym%NO_SNOW .and. &
               ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx) < 0.2 .and. &
               ch(31)%Bt_Toa_Std_3x3(Elem_Idx,Line_Idx) < 0.25) then

               CLDMASK%Dust_Prob(Elem_Idx,Line_Idx) = GET_PIXEL_ABI_DUST_PROB(Elem_Idx,Line_Idx)

               CLDMASK%Dust_Mask(ELem_Idx,Line_Idx) = 0
               if (CLDMASK%Dust_Prob(Elem_Idx,Line_Idx) >= 0.5) then
                   CLDMASK%Dust_Mask(Elem_Idx,Line_Idx) = 1
               endif

           endif

        enddo

     enddo

  end subroutine  GET_SEGMENT_ABI_DUST_PROB

  !----------------------------------------------------------------------
  ! Allocate the Arrays and Read the Data
  !----------------------------------------------------------------------
  subroutine READ_ABI_DUST_LUT()

   integer, dimension(1) :: Sds_Start_1d, Sds_Edge_1d, Sds_Stride_1d
   integer:: Ncid_Dust_Lut
   integer, dimension(1):: Dim_1d
   character(len=200):: File_Name_Full

   !------------------------------------------------------------------------------
   ! read in data
   !------------------------------------------------------------------------------
   File_Name_Full = trim(Ancil_Data_Dir)//"/static/luts/aerosol/"//trim(ABI_Dust_Prob_Table_Name)
   call OPEN_NETCDF(trim(File_Name_Full), Ncid_Dust_Lut)
   if (Ncid_Dust_Lut < 0) then 
     call MESG("Error opening ABI Dust Lut file, default will be used; file = " &
               //trim(File_Name_Full), level=Verb_Lev%WARNING)
     Use_ABI_Dust = 0
     return
   endif

   !---------------------------------------------------------------------------
   ! Read in Dimensions
   !---------------------------------------------------------------------------
!  call READ_NETCDF_DIMENSION_1D(Ncid_Dust_Lut, 'nbins', Dim_1d)
!  Nbins = Dim_1d(1)
 
   !---------------------------------------------------------------------------
   ! Allocate Memory for Tables
   !---------------------------------------------------------------------------
   call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(File_Name_Full), 'nbins', Nbins)

   allocate(Btd_10_11_Cond_Ratio_Table(Nbins))
   allocate(Btd_11_12_Cond_Ratio_Table(Nbins))
   allocate(Btd_85_10_Cond_Ratio_Table(Nbins))
   allocate(OB_Btd_10_11_Cond_Ratio_Table(Nbins))
   allocate(OB_Btd_11_12_Cond_Ratio_Table(Nbins))
   allocate(OB_Btd_85_10_Cond_Ratio_Table(Nbins))

   call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(File_Name_Full), 'bin_min', Bin_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(File_Name_Full), 'bin_width', Bin_Width)
   call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(File_Name_Full), 'dust_fraction', Prior_Yes)

   !---------------------------------------------------------------------------
   ! Read Tables
   !---------------------------------------------------------------------------
   sds_start_1d = 1
   sds_stride_1d = 1

   sds_edge_1d(1) = Nbins

   call READ_NETCDF(Ncid_Dust_Lut,Sds_Start_1d,Sds_Stride_1d,Sds_Edge_1d,&
                    'btd_11_12_cond_ratio_table',Btd_11_12_Cond_Ratio_Table)

   call READ_NETCDF(Ncid_Dust_Lut,Sds_Start_1d,Sds_Stride_1d,Sds_Edge_1d,&
                    'btd_10_11_cond_ratio_table',Btd_10_11_Cond_Ratio_Table)

   call READ_NETCDF(Ncid_Dust_Lut,Sds_Start_1d,Sds_Stride_1d,Sds_Edge_1d,&
                    'btd_85_10_cond_ratio_table',Btd_85_10_Cond_Ratio_Table)

   call READ_NETCDF(Ncid_Dust_Lut,Sds_Start_1d,Sds_Stride_1d,Sds_Edge_1d,&
                    'ob_btd_10_11_cond_ratio_table',OB_Btd_10_11_Cond_Ratio_Table)

   call READ_NETCDF(Ncid_Dust_Lut,Sds_Start_1d,Sds_Stride_1d,Sds_Edge_1d,&
                    'ob_btd_11_12_cond_ratio_table',OB_Btd_11_12_Cond_Ratio_Table)

   call READ_NETCDF(Ncid_Dust_Lut,Sds_Start_1d,Sds_Stride_1d,Sds_Edge_1d,&
                    'ob_btd_85_10_cond_ratio_table',OB_Btd_85_10_Cond_Ratio_Table)

   call CLOSE_NETCDF(Ncid_Dust_Lut)

   call MESG("ABI Dust Lut read in successfully", level=Verb_Lev%DEFAULT)

   Use_ABI_Dust = 1

  end subroutine READ_ABI_DUST_LUT

!----------------------------------------------------------------------------------------
! Interpolate dust for conditions of a pixel
!
!----------------------------------------------------------------------------------------
function GET_PIXEL_ABI_DUST_PROB(Elem_Idx,Line_Idx) Result(Post_Prob)

  integer, intent(in):: Elem_Idx, Line_Idx
  real:: Post_Prob
 
  integer:: Bin_Idx
  real:: z, r
  real:: r_btd_11_12, r_btd_10_11, r_btd_85_10, &
         r_ob_btd_11_12, r_ob_btd_10_11, r_ob_btd_85_10

  !--- initialize output
  Post_Prob =  MISSING_VALUE_REAL4
  r = 1.0

  !--- btd_11_12
  r_btd_11_12 = 1.0
! if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
!    z = ch(31)%Bt_Toa(Elem_Idx,Line_Idx) - ch(32)%Bt_Toa(Elem_Idx,Line_Idx)
!    Bin_Idx = max(1,min(Nbins,int((Z - Bin_Min) / Bin_Width) + 1))
!    r_btd_11_12 = Btd_11_12_Cond_Ratio_Table(Bin_Idx)
! endif

  !--- btd_10_11
  r_btd_10_11 = 1.0
  if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
     z = ch(38)%Bt_Toa(Elem_Idx,Line_Idx) - ch(31)%Bt_Toa(Elem_Idx,Line_Idx)
     Bin_Idx = max(1,min(Nbins,int((Z - Bin_Min) / Bin_Width) + 1))
     r_btd_10_11 = Btd_10_11_Cond_Ratio_Table(Bin_Idx)
!print *, "10_11 = ",  z, r_btd_10_11
  endif

  !--- btd_85_10
  r_btd_85_10 = 1.0
  if (Sensor%Chan_On_Flag_Default(29) == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
     z = ch(29)%Bt_Toa(Elem_Idx,Line_Idx) - ch(38)%Bt_Toa(Elem_Idx,Line_Idx)
     Bin_Idx = max(1,min(Nbins,int((Z - Bin_Min) / Bin_Width) + 1))
     r_btd_85_10 = Btd_85_10_Cond_Ratio_Table(Bin_Idx)
!print *, "85_10 = ",  z, r_btd_85_10
  endif

  !--- ob_btd_10_11
  r_ob_btd_10_11 = 1.0
  if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
     z = ((ch(38)%Bt_Toa(Elem_Idx,Line_Idx) - ch(31)%Bt_Toa(Elem_Idx,Line_Idx)) - &
          (ch(38)%Bt_Toa_Clear(Elem_Idx,Line_Idx) - ch(31)%Bt_Toa_Clear(Elem_Idx,Line_Idx)))
     Bin_Idx = max(1,min(Nbins,int((Z - Bin_Min) / Bin_Width) + 1))
     r_ob_btd_10_11 = OB_Btd_10_11_Cond_Ratio_Table(Bin_Idx)
!print *, "ob_10_11 = ",  z, r_ob_btd_10_11
  endif

  !--- ob_btd_11_12
  r_ob_btd_11_12 = 1.0
! if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
!    z = ((ch(31)%Bt_Toa(Elem_Idx,Line_Idx) - ch(32)%Bt_Toa(Elem_Idx,Line_Idx)) - &
!         (ch(31)%Bt_Toa_Clear(Elem_Idx,Line_Idx) - ch(32)%Bt_Toa_Clear(Elem_Idx,Line_Idx)))
!    Bin_Idx = max(1,min(Nbins,int((Z - Bin_Min) / Bin_Width) + 1))
!    r_ob_btd_11_12 = OB_Btd_11_12_Cond_Ratio_Table(Bin_Idx)
! endif
  
  !--- ob_btd_85_10
  r_ob_btd_85_10 = 1.0
  if (Sensor%Chan_On_Flag_Default(29) == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
     z = ((ch(29)%Bt_Toa(Elem_Idx,Line_Idx) - ch(38)%Bt_Toa(Elem_Idx,Line_Idx)) - &
          (ch(29)%Bt_Toa_Clear(Elem_Idx,Line_Idx) - ch(38)%Bt_Toa_Clear(Elem_Idx,Line_Idx)))
     Bin_Idx = max(1,min(Nbins,int((Z - Bin_Min) / Bin_Width) + 1))
     r_ob_btd_85_10 = OB_Btd_85_10_Cond_Ratio_Table(Bin_Idx)
!print *, "ob_85_10 = ",  z, r_ob_btd_85_10
  endif

  r = r_btd_11_12 * r_btd_10_11 * r_btd_85_10 * r_ob_btd_10_11 * r_ob_btd_11_12 * r_ob_btd_85_10  

  post_prob = 1.0 / (1.0 + r/prior_yes - r)

! print *, "Rs = ", r_btd_11_12 , r_btd_10_11 , r_btd_85_10 , r_ob_btd_10_11 , r_ob_btd_11_12 , r_ob_btd_85_10 , r
! print *, "Final Prob = ",post_prob

end function GET_PIXEL_ABI_DUST_PROB
!-----------------------------------------------------------------------------
!  Deallocate all allocated arrays (call after this data is no longer needed)
!-----------------------------------------------------------------------------
subroutine FORGET_ABI_DUST()
   deallocate(Btd_11_12_Cond_Ratio_Table)
   deallocate(Btd_10_11_Cond_Ratio_Table)
   deallocate(Btd_85_10_Cond_Ratio_Table)
   deallocate(OB_Btd_10_11_Cond_Ratio_Table)
   deallocate(OB_Btd_11_12_Cond_Ratio_Table)
   deallocate(OB_Btd_85_10_Cond_Ratio_Table)
   Use_ABI_Dust = 0

end subroutine FORGET_ABI_DUST

end module CX_DUST_MOD
