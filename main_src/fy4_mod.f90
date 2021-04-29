!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: fy4_mod.f90 (src)
!
! PURPOSE: Module to read level1b FY-4 Chinese Geostationary Satellite data.
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, STAR/NESDIS/NOAA Andrew.Heidinger@noaa.gov
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
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
!--------------------------------------------------------------------------------------

module FY4_MOD

   use CONSTANTS_MOD,only:  &
      Int4 &
    , Real4 &
    , Sym &
    , Missing_Value_Real4

   
   use FILE_TOOLS, only: &
       FILE_SEARCH &
       , get_lun

   use CLAVRX_MESSAGE_MOD,only: &
      Mesg &
    , Verb_Lev

   use CX_NETCDF4_MOD,only: &
      OPEN_NETCDF &
    , CLOSE_NETCDF &
    , READ_NETCDF_GLOBAL_ATTRIBUTE &
    , READ_NETCDF &
    , READ_NETCDF_DIMENSION_1D &
    , READ_NETCDF_DIMENSION_2D

   use PIXEL_COMMON_MOD, only: &
      Sensor &
    , Image &
    , Ch &
    , Nav &
    , Geo &
    , Static_Nav_File &
    , Ancil_Data_Dir

   use PLANCK_MOD, only: &
      COMPUTE_BT_ARRAY

   use CLAVRX_STATIC_NAV_MODULE, only: &
      DETERMINE_BOUNDS_STATIC_NAV &
    , READ_SEGMENT_STATIC_NAV

   use CALIBRATION_CONSTANTS_MOD, only: &
        Planck_Nu

   implicit none

   public READ_FY4_DATE_TIME
   public READ_FY4_INSTR_CONSTANTS
   public READ_FY4_LEVEL1B_DATA
   private :: CONVERT_FY4A_RADIANCE

   character(12), parameter, private :: MODULE_PROMPT = "FY4_MODULE: "

   contains


!==============================================================================
! Determine date and time from global attributes 
!==============================================================================
subroutine READ_FY4_DATE_TIME(File_Path, File_Name, Start_Year, Start_Day, Start_Time, &
                              End_Year, End_Day, End_Time)

   use CX_DATE_TIME_TOOLS_MOD, only: &
       LEAP_YEAR_FCT &
       , JULIAN

   character(len=*), intent(in) :: File_Path
   character(len=*), intent(in) :: File_Name
   integer(kind=int4), intent(out) :: Start_Year
   integer(kind=int4), intent(out) :: Start_Day
   integer(kind=int4), intent(out) :: Start_Time
   integer(kind=int4), intent(out) :: End_Year
   integer(kind=int4), intent(out) :: End_Day
   integer(kind=int4), intent(out) :: End_Time

   character(len=20) :: Tmp_String
   character(len=29) :: Time_String
   character(len=200) :: Search_String
   character ( len = 1020 ) , pointer  :: File_Arr_Dummy(:)
   integer(kind=int4) :: Sd_Id
   integer(kind=int4) :: Ipos
   integer(kind=int4) :: N_Files
   integer(kind=int4) :: Ileap
   integer(kind=int4) :: Month
   integer(kind=int4) :: Day
   integer(kind=int4) :: Start_Hour
   integer(kind=int4) :: Start_Minute
   integer(kind=int4) :: Start_Sec
   integer(kind=int4) :: End_Hour
   integer(kind=int4) :: End_Minute
   integer(kind=int4) :: End_Sec

   ! --- initialize
   N_Files = 0

   ! --- make file names
   Ipos = index(Image%Level1b_Name, '_NOM_')
   Time_String = Image%Level1b_Name(Ipos+5:Ipos+5+29)
   Search_String = 'FY4A-_AGRI--_N_DISK_*_L1-_FDI-_MULT_NOM_'//trim(Time_String)//'*_4000M*.HDF'

   ! --- search for correct 4km file
   File_Arr_Dummy => FILE_SEARCH(trim(Image%Level1b_Path), trim(Search_String), N_Files)

   if (N_Files == 0) then
      call MESG (trim(MODULE_PROMPT)// " No Level1b Files Found, Stopping",level = verb_lev % ERROR)
      stop
   endif
   Image%Level1b_Name = trim(File_Arr_Dummy(1))


   ! --- read start date and convert to integers
   !call H5READATTRIBUTE(trim(File_Path)//trim(File_Name),'Observing Beginning Date', Tmp_String)
   call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(File_Path)//trim(File_Name), 'Observing Beginning Date', Tmp_String)
   read(Tmp_String(1:4), fmt="(I4)") Start_Year
   read(Tmp_String(6:7), fmt="(I2)") Month
   read(Tmp_String(9:10), fmt="(I2)") Day

   ! --- Calculate the date of year start
   Ileap = 0
   Ileap = LEAP_YEAR_FCT(Start_Year)
   call JULIAN(Day, Month, Start_Year, Start_Day)

   ! --- read start time and convert to integers
   !call H5READATTRIBUTE(trim(File_Path)//trim(File_Name),'Observing Beginning Time', Tmp_String)
   call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(File_Path)//trim(File_Name), 'Observing Beginning Time', Tmp_String)
   read(Tmp_String(1:2), fmt="(I2)") Start_Hour
   read(Tmp_String(4:5), fmt="(I2)") Start_Minute
   read(Tmp_String(7:8), fmt="(I2)") Start_Sec

   ! --- read end date and convert to integers
   !call H5READATTRIBUTE(trim(File_Path)//trim(File_Name),'Observing Ending Date', Tmp_String)
   call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(File_Path)//trim(File_Name), 'Observing Ending Date', Tmp_String)
   read(Tmp_String(1:4), fmt="(I4)") End_Year
   read(Tmp_String(6:7), fmt="(I2)") Month
   read(Tmp_String(9:10), fmt="(I2)") Day

   ! --- Calculate the date of year end
   Ileap = 0
   Ileap = LEAP_YEAR_FCT(End_Year)
   call JULIAN(Day, Month, End_Year, End_Day)

   ! --- read end time and convert to integers
   !call H5READATTRIBUTE(trim(File_Path)//trim(File_Name),'Observing Ending Time', Tmp_String)
   call READ_NETCDF_GLOBAL_ATTRIBUTE(trim(File_Path)//trim(File_Name), 'Observing Ending Time', Tmp_String)
   read(Tmp_String(1:2), fmt="(I2)") End_Hour
   read(Tmp_String(4:5), fmt="(I2)") End_Minute
   read(Tmp_String(7:8), fmt="(I2)") End_Sec

   ! --- Calculate start and end time
   Start_Time = ((Start_Hour * 60 + Start_Minute) * 60 + Start_Sec) * 1000
   End_Time = ((End_Hour * 60 + End_Minute) * 60 + End_Sec) * 1000


end subroutine READ_FY4_DATE_TIME

!==============================================================================
! Read instrument file
!==============================================================================
subroutine READ_FY4_INSTR_CONSTANTS(Instr_Const_File)

   use CALIBRATION_CONSTANTS_MOD

   character(len=*), intent(in) :: Instr_Const_file

   integer :: Ios0, Erstat
   integer :: Instr_Const_Lun
   integer :: Dummy

   ! --- get lun and open file
   Instr_Const_Lun = GET_LUN()

   open(unit=Instr_Const_Lun,file=trim(Instr_Const_File),status="old",position="rewind",action="read",iostat=ios0)
   Erstat = 0
   if (Ios0 /= 0) then
      Erstat = 19
      call MESG (trim(MODULE_PROMPT)//"  Error Opening Constants File",level = verb_lev % DEFAULT)
      stop 19
   endif

   ! --- read file
   read(unit=Instr_Const_Lun,fmt=*) ! Text
   read(unit=Instr_Const_Lun,fmt="(a4)") Sat_Name
   read(unit=Instr_Const_Lun,fmt=*) ! WMO_ID
   read(unit=Instr_Const_Lun,fmt=*) ! Number of Channels
   read(unit=Instr_Const_Lun,fmt=*) ! Text 
   read(unit=Instr_Const_Lun,fmt=*) ! clavrx channel map
   read(unit=Instr_Const_Lun,fmt=*) ! Text 
   read(unit=Instr_Const_Lun,fmt=*) Solar_Ch20
   read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
   read(unit=Instr_Const_Lun,fmt=*) ! Text 
   read(unit=Instr_Const_lun,fmt=*) Planck_A1(20), Planck_A2(20),Planck_Nu(20) ! ch20 <-- ch7
   read(unit=Instr_Const_lun,fmt=*) Planck_A1(21), Planck_A2(21),Planck_Nu(21) ! ch22 <-- ch8
   read(unit=Instr_Const_lun,fmt=*) Planck_A1(27), Planck_A2(27),Planck_Nu(27) ! ch27 <-- ch9 (more like 37) 
   read(unit=Instr_Const_lun,fmt=*) Planck_A1(28), Planck_A2(28),Planck_Nu(28) ! ch28 <-- ch10
   read(unit=Instr_Const_lun,fmt=*) Planck_A1(29), Planck_A2(29),Planck_Nu(29) ! ch29 <-- ch11
   read(unit=Instr_Const_lun,fmt=*) Planck_A1(31), Planck_A2(31),Planck_Nu(31) ! ch31 <-- ch12
   read(unit=Instr_Const_lun,fmt=*) Planck_A1(32), Planck_A2(32),Planck_Nu(32) ! ch32 <-- ch13
   read(unit=Instr_Const_lun,fmt=*) Planck_A1(33), Planck_A2(33),Planck_Nu(33) ! ch33 <-- ch14
   close(unit=Instr_Const_lun)

   !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
   Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

   call MESG (trim(MODULE_PROMPT)// " Instrument Constants Read in Successfully",level = verb_lev % DEFAULT)

end subroutine READ_FY4_INSTR_CONSTANTS

!==============================================================================
! Read level1b 4km only
!==============================================================================
subroutine READ_FY4_LEVEL1B_DATA(Segment_Number, L1b_File, Error_Out)

   use VIEWING_GEOMETRY_MOD, only: &
      RELATIVE_AZIMUTH &
    , SCATTERING_ANGLE &
    , GLINT_ANGLE &
    , POSSOL

   integer(kind=int4), intent(in) :: Segment_Number
   character(len=*), intent(in) :: L1b_File
   integer(kind=int4), intent(out) :: Error_Out

   integer(kind=int4) :: Nx_Start, Nx_End, Ny_Start, Ny_End, N_Seg_Lines
   integer(kind=int4) :: i, j
   integer(kind=int4) :: I_Ch
   integer(kind=int4) :: Sd_Id, Nav_Sd_Id
   integer(kind=int4) :: Nav_Number_Of_Elements,Nav_Number_Of_Lines
   integer(kind=int4), dimension(1) :: Dim_Lut
   integer(kind=int4), dimension(2) :: Offset_Read, Dim_Seg, Stride, Dim_Full
   !integer(kind=int4), dimension(:,:), pointer :: I2d_Buffer
   integer(kind=int4), dimension(:,:), allocatable :: I2d_Buffer
   integer(kind=int4), parameter :: Num_Chan_Sensor = 14
   integer(kind=int4), parameter :: Fill_Value = 65535
   integer(kind=int4), dimension(Num_Chan_Sensor) :: Modis_Chn_List
   real(kind=real4), dimension(:), allocatable :: R1d_Buffer
   real(kind=real4), dimension(2,Num_Chan_Sensor) :: Rad_Coef
   real(kind=real4), dimension(:,:), pointer :: R2d_Buffer
   character(len=2) :: Ch_String
   character(len=12) :: Sds_Name
   character(len=1020) :: Static_Nav_Full_File

   Error_Out = 0
   Rad_Coef = -1.
   Sd_Id = -1
   Nav_Sd_Id = -1

   ! --- mapping fy4 to modis
   !              WL 047 065 083 137  161 222 372  372  625  710  850  108  120  135
   !             FY4  1   2   3   4    5   6    7   8    9    10   11   12   13   14
   Modis_Chn_List = [ 3 , 1 , 2 , 26 , 6 , 7 , 20 , 21 , 27 , 28 , 29 , 31 , 32 , 33 ]



   ! --- Set the expected location of the static navigation file.
   Static_Nav_Full_File = trim(Ancil_Data_Dir)//"static/static_nav/"//trim(Static_Nav_File)

   ! --- Open static navigation file.
   call OPEN_NETCDF(Static_Nav_Full_File,Nav_Sd_Id)

   ! ---  determine bounds of static navigation data
   call DETERMINE_BOUNDS_STATIC_NAV(trim(Static_Nav_Full_File),Nav%Limit_Flag, &
                                     Nav%Lat_South_Limit, Nav%Lon_West_Limit, &
                                     Nav%Lat_North_Limit, Nav%Lon_East_Limit, &
                                     Image%X_Stride, Image%Y_Stride, &
                                     Nav_Number_Of_Elements,Nav_Number_Of_Lines)

   ! --- read static navigation latitude, longitude, sensor_zentith and sensor_azimiuth
   call READ_SEGMENT_STATIC_NAV(Nav_Sd_Id)

   ! --- adjust sensor azimuth to be -/+180
   Geo%Sataz = Geo%Sataz - 180.0

   !--- get solar zenith angle 
   do i = 1, Image%Number_Of_Elements
      do j = 1, Image%Number_Of_Lines_Read_This_Segment
         call POSSOL(int(Image%Start_Doy,kind=int4),Image%Mean_Time_Hours,Nav%Lon(i,j), &
                     Nav%Lat(i,j),Geo%Solzen(i,j),Geo%Solaz(i,j))
      end do
   end do

   !--- relative azimuth
   Geo%Relaz = RELATIVE_AZIMUTH(Geo%Solaz, Geo%Sataz)

   !--- scattering angle
   Geo%Scatangle = SCATTERING_ANGLE(Geo%Solzen, Geo%Satzen, Geo%Relaz)

   !--- glint angle
   Geo%Glintzen = GLINT_ANGLE(Geo%Solzen, Geo%Satzen, Geo%Relaz)


   !--- open level1b file to read
   call OPEN_NETCDF(trim(Image%Level1b_Path)//trim(Image%Level1b_Name), Sd_Id)

   ! --- make start, end, stride
   Nx_Start = 1
   Nx_End = Image%Number_Of_Elements
   Ny_Start = (Segment_Number-1) * Image%Number_of_Lines_Per_Segment + 1
   Ny_End = Image%Number_of_Lines_Per_Segment

!   call READ_NETCDF_DIMENSION_2D(Sd_Id, 'CALIBRATION_COEF(SCALE+OFFSET)', Dim_Full)

   ! --- constrain to size of data
   if ((Segment_Number * Image%Number_of_Lines_Per_Segment) .gt. Image%Number_Of_Lines) &
              Ny_End = Image%Number_Of_Lines - Ny_Start + 1

   Stride = (/1, 1/)
   Offset_Read = (/Nx_Start, Ny_Start/)
   Dim_Seg = (/Nx_End, Ny_End/)

   ! --- read scale factor/offset
   call READ_NETCDF(Sd_Id, (/1,1/), Stride, (/2,Num_Chan_Sensor/), & 
                    'CALIBRATION_COEF(SCALE+OFFSET)', Rad_Coef)

   ! --- read 4km file
   ! --- loop over channels
   do I_Ch = 1, Num_Chan_Sensor

      ! --- check if channel is on
      if (Sensor%Chan_On_Flag_Default (Modis_Chn_List(I_Ch)) .eq. sym%NO) cycle

      ! --- allocate buffer
      if (.not. allocated(I2D_Buffer)) allocate(I2D_Buffer(Dim_Seg(1),Dim_Seg(2)))

      ! --- make string band number
      write ( Ch_String, '(i0.2)' ) I_Ch

      ! --- initialize output
      if (I_Ch .le. 6) then
         Ch (Modis_Chn_List(I_Ch)) % Ref_Toa (:, 1:Ny_End) = Missing_Value_Real4

         ! --- read LUT (for reflectance only)
         Sds_Name = 'CALChannel'//Ch_String
         call READ_NETCDF_DIMENSION_1D(Sd_Id, trim(Sds_Name), Dim_Lut)
         if (.not. allocated(R1D_Buffer)) allocate(R1D_Buffer(Dim_Lut(1)))
         call READ_NETCDF(Sd_Id, (/1/), (/1/), Dim_Lut, trim(Sds_Name), R1d_Buffer)

      else
         Ch (Modis_Chn_List(I_Ch)) % Rad_Toa (:, 1:Ny_End) = Missing_Value_Real4
      endif

      ! --- read DN
      Sds_Name = 'NOMChannel'//Ch_String
      call READ_NETCDF(Sd_Id, Offset_Read, Stride, Dim_Seg, trim(Sds_Name), I2d_Buffer)
   

      ! --- compute reflectance using LUT
      if (I_Ch .le. 6) then
         ! --- loop over DN (index) for LUT. TODO Maybe more optimal way to do this?
         do i = 1, Nx_End
            do j = 1, Ny_End
               if (I2d_Buffer(i,j) .eq. Fill_Value) cycle
                  Ch (Modis_Chn_List(I_Ch)) % Ref_Toa (i,j) = R1d_Buffer(I2d_Buffer(i,j)) * 100.0
            enddo
         enddo         

        

         ! --- deallocate LUT
         deallocate(R1D_Buffer)

      else 

         ! --- initialize bt
         Ch (Modis_Chn_List(I_Ch)) % Bt_Toa (:, 1:Ny_End) = Missing_Value_Real4


         ! --- WCS MOD for Radiances
         
         ! --- initialize radiance temporary buffer
         Ch (Modis_Chn_List(I_Ch)) % Rad_Toa (:, 1:Ny_End)  = Missing_Value_Real4

         
         ! --- convert to native units
         where ( I2d_Buffer .ne. Fill_Value)
            Ch (Modis_Chn_List(I_Ch)) % Rad_Toa (:, 1:Ny_End) = Rad_Coef(1,I_Ch) * I2d_Buffer + Rad_Coef(2, I_Ch)
         endwhere
         
         
         ! --- convert to NOAA units
         call CONVERT_FY4A_RADIANCE (Ch (Modis_Chn_List(I_Ch)) % Rad_Toa (:,1:Ny_End), &
                                Planck_Nu(Modis_Chn_List(I_Ch)), Missing_Value_Real4) 

         ! --- calculate bt
         call COMPUTE_BT_ARRAY (Ch(Modis_Chn_List(I_Ch))%Bt_Toa (:,1:Ny_End), &
                                Ch(Modis_Chn_List(I_Ch))%Rad_Toa (:,1:Ny_End), &
                                Modis_Chn_List(I_Ch), Missing_Value_Real4)

      endif

      ! --- deallocate 
      deallocate (I2d_Buffer)
      

    enddo ! end loop over channels

    !--- close file
    call CLOSE_NETCDF(Sd_Id)


    ! --- global variables which have to be set
    Image%Number_Of_Lines_Read_This_Segment = Ny_End
    do i = 1, Image%Number_Of_Lines_Per_Segment
        Image%Scan_Number(i) = Ny_Start + i - 1
    end do

end subroutine READ_FY4_LEVEL1B_DATA

!==============================================================================


!--------------------------------------------------
! Function Name: CONVERT_FY4A_RADIANCE
!
! Function:
!    Convert to units of the FY4A radiance values to what expected by CLAVR-x
!    Assumption is that FY4 is like FY3. Per YL, this seems to give correct units
!
! Description: 
!   
! Calling Sequence: rad_new =
! CONVERT_FY4A_RADIANCE(rad_old,nu,missing_value)
!   
!
! Inputs:
!   rad_old = radiance in units of W/m^2/micron/str (2d array)
!   nu = channels equivalent width in units of cm^-1
!   missing_value = value assigned to missing radiance values
!
! Outputs: 
!   rad_new = radiance in units of mW/m^2/cm^-1/str (2d array)
!
! Dependencies:
!
! Restrictions:  None
!
! Reference: algebraic manipulation of Planck Equation
! ---------------------------------------------------------------------------------------
subroutine CONVERT_FY4A_RADIANCE(Radiance,Nu,Missing_Value)
      real (kind=real4), dimension(:,:), intent(inout):: Radiance
      real (kind=real4), intent(in):: Nu
      real (kind=real4), intent(in):: Missing_Value

      where(Radiance /= Missing_Value)
         Radiance = Radiance * (((10000.0 / Nu )**2) / 10.0)
      end where

      return

end subroutine CONVERT_FY4A_RADIANCE



end module FY4_MOD
