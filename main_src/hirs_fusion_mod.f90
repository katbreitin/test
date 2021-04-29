! $Id: hirs_fusion_mod.f90 3785 2020-03-30 12:46:14Z heidinger $
!------------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: hirs_fusion_module.f90 (src)
!       HIRS_FUSION_MODULE (program)
!
! PURPOSE: 
!
! DESCRIPTION: 
!
! AUTHORS:
!  Mike Foster, Mike.Foster@ssec.wisc.edu
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
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
!------------------------------------------------------------------------------------------
module HIRS_FUSION_MOD
   use CONSTANTS_MOD, only: MISSING_VALUE_REAL4
   use PIXEL_COMMON_MOD,only: AVHRR_Fusion_Flag, Image, Sensor, Ch
   use CX_REAL_BOOLEAN_MOD
      use cx_sds_io_mod,only: &
           cx_sds_finfo &
         , cx_sds_varinfo &
         , cx_sds_read &
         , MAXNCNAM

      use CALIBRATION_CONSTANTS_MOD,only: &
         Planck_Nu &
         , Planck_Nu_11um_Sndr &
         , Planck_Nu_12um_Sndr &
         , Planck_Nu_375um_Sndr &
         , Planck_Nu_145um_Sndr &
         , Planck_Nu_147um_Sndr &
         , Planck_Nu_149um_Sndr &
         , Planck_Nu_445um_Sndr &
         , Planck_Nu_457um_Sndr

      use PIXEL_COMMON_MOD,only: &
         Bt_11um_Sounder &
         , Bt_12um_Sounder &
         , Bt_375um_Sounder &
         , Bt_145um_Sounder &
         , Bt_147um_Sounder &
         , Bt_149um_Sounder &
         , Bt_445um_Sounder &
         , Bt_457um_Sounder

      use PLANCK_MOD, only: &
         COMPUTE_BT_ARRAY &
         , COMPUTE_BT_ARRAY_SOUNDER

      use CALIBRATION_CONSTANTS_MOD
      use FILE_TOOLS , only: GETLUN

   implicit none  
   private
   
   public:: HIRS_AVHRR_FUSION_PREPERATION
   public:: READ_FUSION_HIRS_INSTR_CONSTANTS
   public:: READ_HIRS_DATA
   public:: REPLACE_AVHRR_WITH_HIRS

   contains

   subroutine HIRS_AVHRR_FUSION_PREPERATION ( File_1b)
      character ( len = * ), intent(inout) :: File_1b
      AVHRR_Fusion_Flag = .true.
      
      ! -reset level1b to AVHRR file and add HIRS file as file_fusion.
      Image%Level1b_Fusion_Name = File_1b
!      File_1b = File_1b(1:len(trim(File_1b)) - 9 ) //'gz'
      File_1b = File_1b(1:len(trim(File_1b)) - 10 )

   end subroutine HIRS_AVHRR_FUSION_PREPERATION
   ! ***********************************************************************************************************************
   ! SUBROUTINE NAME: Read_hirs_data 
   ! ORIGINAL AUTHOR: Mike Foster
   ! CREATION DATE:  March, 2018
   ! DESCRIPTION: This subroutine reads data from HIRS in ch structure and additional
   !                global variables Bt_11um_Sounder, Bt_12um_Sounder, Bt_375um_Sounder, Bt_145um_Sounder, and Bt_147um_Sounder
   ! CALLING ARGUMENTS:
   !       seg_num   :  number of segement
   !   
   !  USED GLOBAL VARIABLES: see specified variables in USE specifications of external modules
   !
   !   CALLED LOCATIONS:
   !       For CLAVR-x in sensor_mod.f90
   !
   !  MODIFICATION HISTORY
   !  Date                       Developer             Action  
   !     5/21/2018                A.Walther           Revised software and merged from branch "fusion_hirs_module" into trunk.
   !
   !  ***************************************************************************************************************************
   subroutine READ_HIRS_DATA (seg_num)

!     use cx_sds_io_mod,only: &
!          cx_sds_finfo & 
!        , cx_sds_varinfo &
!        , cx_sds_read &
!        , MAXNCNAM
!        
!     use CALIBRATION_CONSTANTS_MOD,only: &
!        Planck_Nu &
!        , Planck_Nu_11um_Sndr &
!        , Planck_Nu_12um_Sndr &
!        , Planck_Nu_375um_Sndr &
!        , Planck_Nu_145um_Sndr &
!        , Planck_Nu_147um_Sndr &
!        , Planck_Nu_149um_Sndr &
!        , Planck_Nu_445um_Sndr &
!        , Planck_Nu_457um_Sndr
!     
!     use PIXEL_COMMON_MOD,only: &
!        Bt_11um_Sounder &
!        , Bt_12um_Sounder &
!        , Bt_375um_Sounder &
!        , Bt_145um_Sounder &
!        , Bt_147um_Sounder &
!        , Bt_149um_Sounder &
!        , Bt_445um_Sounder &
!        , Bt_457um_Sounder
!     
!     use PLANCK_MOD, only: &
!        compute_bt_array &
!        , compute_bt_array_sounder
!     
!     implicit none

      integer, intent(in) :: seg_num
      
      integer :: ftype      
      integer :: nsds
      integer :: Natt
     
      character ( len = MAXNCNAM), allocatable :: Sds_Name(:)
      character ( len = MAXNCNAM), allocatable :: Att_Name(:)
      integer :: status
      real,allocatable :: Rad_Hirs(:,:)
      character(len=1020) :: File_Local
      integer :: sds_start (2)
      integer :: Sds_Count (2)
      integer :: Ndim
      integer :: Dims(10)
      integer :: hirs_chan_idx
      integer, parameter :: modis_hirs_list(19) = [0,0,0,36,35,34,33,31,30,32,28,27, 0,25,24, 0, 0,23,20]
      ! HIRS channels mapping                      1,2,3, 4, 5,,6,,7,,8,,9,10,11,12,13,14,15,16,17,18,19
      character(len=120):: hirs_Sds_Name
      integer :: Modis_Idx
      real :: nu
      real, parameter :: MISSING_VALUE_REAL4 = -999.
      integer :: shape_hirs(2)
        
      File_Local = trim(Image%Level1b_Path)//trim(Image%Level1b_Fusion_Name)

      status = cx_sds_finfo (File_Local, ftype,nsds,Sds_Name,Natt,Att_Name)
      status = cx_sds_varinfo ( trim(File_Local),  Sds_Name(1) , ftype,Natt,Att_Name, Ndim, Dims)
     
      sds_start(1) = 1
      sds_start(2) = (Seg_num-1) * Image%Number_of_Lines_Per_Segment + 1
      Sds_Count(1) = Dims(1)
      Sds_Count(2) = Image%Number_of_Lines_Read_This_Segment
    
      if (allocated(Rad_Hirs)) deallocate(Rad_Hirs)
    
      !--- Read in HIRS channels that don't fit into the MODIS-based array
      
      ! HIRS Channel 1 : 14.9 microns
      write (Hirs_Sds_Name,"(A4,I0.2)") "HIRS", 1
      status = cx_sds_read(File_Local,trim(Hirs_Sds_Name),Rad_Hirs, count = Sds_Count, start = sds_start)
      nu = Planck_Nu_149um_Sndr
      call CONVERT_HIRS_RADIANCE ( Rad_Hirs, nu, -999.)
      shape_hirs = shape ( rad_hirs)
      call COMPUTE_BT_ARRAY_SOUNDER (Bt_149um_Sounder(:,1:shape_hirs(2)), &
                                    Rad_Hirs , &
                                    149, MISSING_VALUE_REAL4)
      if (allocated(Rad_Hirs)) deallocate(Rad_Hirs)

      ! HIRS Channel 2 : 14.7 microns
      write (Hirs_Sds_Name,"(A4,I0.2)") "HIRS", 2
      status = cx_sds_read(File_Local,trim(Hirs_Sds_Name),Rad_Hirs, count = Sds_Count, start = sds_start)
      nu = Planck_Nu_147um_Sndr
      call CONVERT_HIRS_RADIANCE ( Rad_Hirs, nu, -999.)
      shape_hirs = shape ( rad_hirs)
      call COMPUTE_BT_ARRAY_SOUNDER (Bt_147um_Sounder(:,1:shape_hirs(2)), &
                                    Rad_Hirs , &
                                    147, MISSING_VALUE_REAL4)
      if (allocated(Rad_Hirs)) deallocate(Rad_Hirs)

      ! HIRS Channel 3 : 14.5 microns
      write (Hirs_Sds_Name,"(A4,I0.2)") "HIRS", 3
      status = cx_sds_read(File_Local,trim(Hirs_Sds_Name),Rad_Hirs, count = Sds_Count, start = sds_start)
      nu = Planck_Nu_145um_Sndr
      call CONVERT_HIRS_RADIANCE ( Rad_Hirs, nu, -999.)
      shape_hirs = shape ( rad_hirs)
      call COMPUTE_BT_ARRAY_SOUNDER (Bt_145um_Sounder(:,1:shape_hirs(2)), &
                                    Rad_Hirs , &
                                    145, MISSING_VALUE_REAL4)
      if (allocated(Rad_Hirs)) deallocate(Rad_Hirs)

      ! HIRS Channel 13 : 4.57 microns
      write (Hirs_Sds_Name,"(A4,I0.2)") "HIRS", 13
      status = cx_sds_read(File_Local,trim(Hirs_Sds_Name),Rad_Hirs, count = Sds_Count, start = sds_start)
      nu = Planck_Nu_457um_Sndr
      call CONVERT_HIRS_RADIANCE ( Rad_Hirs, nu, -999.)
      shape_hirs = shape ( rad_hirs)
      call COMPUTE_BT_ARRAY_SOUNDER (Bt_457um_Sounder(:,1:shape_hirs(2)), &
                                    Rad_Hirs , &
                                    457, MISSING_VALUE_REAL4)
      if (allocated(Rad_Hirs)) deallocate(Rad_Hirs)

      ! HIRS Channel 16 : 4.45 microns
      write (Hirs_Sds_Name,"(A4,I0.2)") "HIRS", 16
      status = cx_sds_read(File_Local,trim(Hirs_Sds_Name),Rad_Hirs, count = Sds_Count, start = sds_start)
      nu = Planck_Nu_445um_Sndr
      call CONVERT_HIRS_RADIANCE ( Rad_Hirs, nu, -999.)
      shape_hirs = shape ( rad_hirs)
      call COMPUTE_BT_ARRAY_SOUNDER (Bt_445um_Sounder(:,1:shape_hirs(2)), &
                                    Rad_Hirs , &
                                    445, MISSING_VALUE_REAL4)
      if (allocated(Rad_Hirs)) deallocate(Rad_Hirs)

      !--- loop through channels
      do Hirs_Chan_Idx = 1, 19
         Modis_Idx = Modis_Hirs_List(Hirs_Chan_Idx)
         if ( Modis_Idx .eq. 0 ) cycle
         if ( Sensor%Chan_On_Flag_Default (Modis_Idx) == 0) cycle
         
         write (Hirs_Sds_Name,"(A4,I0.2)") "HIRS", Hirs_Chan_Idx
        
         status = cx_sds_read(File_Local,trim(Hirs_Sds_Name),Rad_Hirs, count =Sds_Count, start = sds_start)
         
         ! - this modifies out   
         nu = Planck_Nu(Modis_Idx)  
         if ( Modis_Idx .eq. 31)  nu = Planck_Nu_11um_Sndr
         if ( Modis_Idx .eq. 32)  nu = Planck_Nu_12um_Sndr
         if ( Modis_Idx .eq. 20)  nu = Planck_Nu_375um_Sndr  
                  
         call CONVERT_HIRS_RADIANCE ( Rad_Hirs, nu , -999.)
         
         shape_hirs = shape ( rad_hirs)
         
         ! some channels need special treatments
         select case  (Modis_Idx)
         case(31)      
            
             call COMPUTE_BT_ARRAY_SOUNDER (Bt_11um_Sounder(:,1:shape_hirs(2)), &
                                    Rad_Hirs , &
                                    Modis_Idx, MISSING_VALUE_REAL4)     
         case(32)    
                        
            call COMPUTE_BT_ARRAY_SOUNDER (Bt_12um_Sounder(:,1:shape_hirs(2)) , &
                                    Rad_Hirs, &
                                    Modis_Idx, MISSING_VALUE_REAL4)
                                   
         case(20)
                     
            call COMPUTE_BT_ARRAY_SOUNDER (Bt_375um_Sounder(:,1:shape_hirs(2)) , &
                                    Rad_Hirs, &
                                    Modis_Idx, MISSING_VALUE_REAL4)

         case default
             
            ch ( Modis_Idx) % Rad_Toa(:,1:shape_hirs(2))  = Rad_Hirs            
            call COMPUTE_BT_ARRAY (Ch(Modis_Idx) % Bt_Toa(:,1:shape_hirs(2)), &
                                    Rad_Hirs, &
                                    Modis_Idx, MISSING_VALUE_REAL4)

            !--- note that these data from from the sounder, not the imager
            Ch(Modis_Idx)%Source(:,1:shape_hirs(2)) = 1                     
         end select
         
        
         if (allocated(Rad_Hirs)) deallocate(Rad_Hirs)
      end do
      
   end subroutine READ_HIRS_DATA
   
!------------------------------------------------------------------------------------------
! Function Name: CONVERT_HIRS_RADIANCE
!
! Function:
!   Convert to radiance units from that used
!   in the HIRS Fusion files to that expected by CLAVR-x
!
! Description: 
!   
! Calling Sequence: rad_new =
!   convert_hirs_radiance(rad_old,nu,missing_value)
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
!------------------------------------------------------------------------------------------
   subroutine CONVERT_HIRS_RADIANCE(Radiance,Nu,Missing_Value)
      real , dimension(:,:), intent(inout):: Radiance
      real , intent(in):: Nu
      real , intent(in):: Missing_Value

      where(Radiance /= Missing_Value)
         Radiance = Radiance * (((10000.0 / Nu )**2) / 10.0)
      end where

      return

   end subroutine CONVERT_HIRS_RADIANCE
   
!------------------------------------------------------------------------------------------
! read the FUSION HIRS constants into memory
!------------------------------------------------------------------------------------------
   subroutine READ_FUSION_HIRS_INSTR_CONSTANTS(Instr_Const_file)
!     use CALIBRATION_CONSTANTS_MOD
!     use FILE_TOOLS , only: GETLUN

      character(len=*), intent(in):: Instr_Const_file
      integer:: ios0, erstat
      integer:: Instr_Const_lun
      real:: Planck_A1_Ch10, Planck_A2_Ch10, Planck_Nu_Ch10

      Instr_Const_lun = GETLUN()

      open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)
      
      erstat = 0
      if (ios0 /= 0) then
         erstat = 19
         print *, "Error opening HIRS constants file, ios0 = ", ios0
         stop 19
      end if

      !--- read in hirs planck fits.  some channels are skipped
      !--- channels that overlap with avhrr are read into separate variables
      !--- note hirs ch10 is 8.5 micron on some and 12 micron on others
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1_149um_Sndr, Planck_A2_149um_Sndr, Planck_Nu_149um_Sndr !hirs ch1
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1_147um_Sndr, Planck_A2_147um_Sndr, Planck_Nu_147um_Sndr !hirs ch2
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1_145um_Sndr, Planck_A2_145um_Sndr, Planck_Nu_145um_Sndr !hirs ch3
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1(36), Planck_A2(36), Planck_Nu(36)   !hirs ch4
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1(35), Planck_A2(35), Planck_Nu(35)   !hirs ch5
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1(34), Planck_A2(34), Planck_Nu(34)   !hirs ch6
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1(33), Planck_A2(33), Planck_Nu(33)   !hirs ch7
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1_11um_Sndr , Planck_A2_11um_Sndr , Planck_Nu_11um_Sndr  !hirs ch8
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1(30), Planck_A2(30), Planck_Nu(30)   !hirs ch9
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1_Ch10, Planck_A2_Ch10, Planck_Nu_Ch10                   !hirs ch10
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1(28), Planck_A2(28), Planck_Nu(28)   !hirs ch11
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1(27), Planck_A2(27), Planck_Nu(27)   !hirs ch12
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1_457um_Sndr, Planck_A2_457um_Sndr, Planck_Nu_457um_Sndr !hirs ch13
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1(25), Planck_A2(25), Planck_Nu(25)   !hirs ch14
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1(24), Planck_A2(24), Planck_Nu(24)   !hirs ch15
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1_445um_Sndr, Planck_A2_445um_Sndr, Planck_Nu_445um_Sndr !hirs ch16
      read(unit=Instr_Const_lun,fmt=*)  !hirs ch17
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1(23), Planck_A2(23), Planck_Nu(23)   !hirs ch18
      read(unit=Instr_Const_lun,fmt=*)  Planck_A1_375um_Sndr, Planck_A2_375um_Sndr, Planck_Nu_375um_Sndr !hirs ch19

      close(unit=Instr_Const_lun)

      !--- handle switch of hirs ch10 from 8.5 to 12 micron
      !--- 8.5 micron for noaa 5-10,12
      !--- 12 micron for noaa 11,14-19,1,2
      !--- note, we are storing both in 12um_Sndr as we always have but maybe we
      !--- should change this
      select case (Sensor%WMO_Id)

        case(3,4,203,205,206,207,208,209,223)     !ch10 = 12 micron
           Planck_A1_12um_Sndr = Planck_A1_Ch10
           Planck_A2_12um_Sndr = Planck_A2_Ch10
           Planck_Nu_12um_Sndr  = Planck_Nu_Ch10

        case(200,201,202,204,706,707,708)         !ch10 = 8.5 micron
           Planck_A1_12um_Sndr = Planck_A1_Ch10
           Planck_A2_12um_Sndr = Planck_A2_Ch10
           Planck_Nu_12um_Sndr  = Planck_Nu_Ch10

      end select

   end subroutine READ_FUSION_HIRS_INSTR_CONSTANTS

   !--------------------------------------------------------------------------------
   ! for some avhrrs, some thermal channels are missing or low-quality
   ! this routine will swap in hirs data for avhrr data and set the ch(x)%Source
   ! to indicate when this is done
   !--------------------------------------------------------------------------------
   subroutine REPLACE_AVHRR_WITH_HIRS()
     
      !--- replace 3.75 micron bt when off on AVHRR's 
      if (allocated(ch(20)%Bt_Toa) .and. allocated(Bt_375um_Sounder)) then
         where((ch(20)%Bt_Toa .eqr. MISSING_VALUE_REAL4) .and. &
               (Bt_375um_Sounder .ner. MISSING_VALUE_REAL4))  
            ch(20)%Bt_Toa = Bt_375um_Sounder
            ch(20)%Source = 1
         endwhere
      endif
   end subroutine REPLACE_AVHRR_WITH_HIRS

end module HIRS_FUSION_MOD
