!$Id: mvcm_read_mod.f90 3936 2020-08-03 21:58:13Z dbotambekov $
!--------------------------------------------------------------------------------------
! This module reads the MVCM cloud mask to allow CLAVR-x to simulate the MODAWG chain
!
! Previously, these files were read in via the IFF module but we need to be able these
! files from the original NASA Level1b data.
!
!
!VGEOM_snpp_d20130426_t101800_c20161107162802.nc
!IFFCMO_npp_d20130426_t101800_c20161029052326_ssec_dev.hdf
!
!
!IFFCMO_aqua_d20121229_t051500_c20170114104543_ssec_dev.hdf 
!
! Modification History
!
! April 2019
! New names are: CLDMSK_L2_VIIRS_SNPP.A2019076.0230.001.2019076132237.nc
! 
!-------------------------------------------------------------------------------------
module MVCM_READ_MOD

    use PIXEL_COMMON_MOD, only: &
             Image &
           , Cldmask &
           , Use_Aux_Flag &
           , Cloud_Mask_Aux_Read_Flag &
           , Cloud_Type_Aux_Read_Flag
    use CX_DATE_TIME_TOOLS_MOD, only: &
             JULIAN &
           , COMPUTE_MONTH &
           , COMPUTE_DAY &
           , LEAP_YEAR_FCT
    use FILE_TOOLS, only: FILE_SEARCH
    use CONSTANTS_MOD, only: &
             Real4 &
           , Int4 &
           , Int2 &
           , Int1 &
           , Exe_Prompt &
           , Missing_Value_Real4 &
           , Sym &
           , Missing_Value_Int1
    use CX_NETCDF4_MOD

    implicit none
    private
    public:: DETERMINE_MVCM_NAME
    public:: READ_MVCM_DATA

    character(len=13), parameter:: MVCM_PROMPT="MVCM_MODULE:"

    contains

!------------------------------------------------------------------------------------
! Determine the name of the MVCM file for the Level1b file
!------------------------------------------------------------------------------------
 subroutine DETERMINE_MVCM_NAME(Seg_Idx)

  integer, intent(in):: Seg_Idx
  character(len=100):: Search_String
  character(len=1020), dimension(:), pointer:: Files
  integer:: Num_Files
  character(len=4):: Year_String, Time_String
  character(len=2):: Month_String, Dom_String
  character(len=3):: Doy_String
  integer:: Year, Month, Dom, Doy,Ileap

  Image%Auxiliary_Cloud_Mask_File_Name = 'no_file'

  !--- NASA VIIRS Level1b
  if (index(Image%Level1b_Name,'VGEOM') == 1) then 

    !0        1         2         3         4
    !12345678901234567890123456789012345678901234567890
    !VGEOM_snpp_d20140205_t115400_c20170405013130.nc
    !VNPCLDMK.A2014036.1154.001.2017300132923.nc
    !--- search should be the date and start time (ie. d20130426_t083000
    !Search_String = 'IFFCMO_npp_'//Image%Level1b_Name(12:28)//'*.hdf'
    !Search_String = 'VCLDMK_snpp_'//Image%Level1b_Name(12:28)//'*.nc'
    Year_String = Image%Level1b_Name(13:16)
    Month_String = Image%Level1b_Name(17:18)
    Dom_String = Image%Level1b_Name(19:20)
    Time_String = Image%Level1b_Name(23:26)
    read(Year_String,*) Year
    read(Month_String,*) Month
    read(Dom_String,*) Dom
    call JULIAN(Dom, Month, Year, Doy)
    write(Doy_String,fmt="(I3.3)") Doy
    Search_String = 'VNPCLDMK.A'//Year_String//Doy_String//'.'//Time_String//'*.nc'

    Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0 .or. Num_Files > 1) then
        print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Not Found, "
        return
    endif

    Image%Auxiliary_Cloud_Mask_File_Name = Files(1)

    if (Seg_Idx == 1) &
       print *, EXE_PROMPT, MVCM_PROMPT, "NASA VIIRS Level1b MVCM File Found, ", &
            trim(Image%Auxiliary_Cloud_Mask_File_Name)

  endif

  !--- NEW version of NASA VIIRS Level1b SNPP
  if (index(Image%Level1b_Name,'VNP03MOD') == 1) then
    !0        1         2         3         4
    !12345678901234567890123456789012345678901234567890
    !VNP03MOD.A2014036.1154.001.2017300063048.uwssec.nc
    !VNPCLDMK.A2014036.1154.001.2017300132923.nc
    !CLDMSK_L2_VIIRS_SNPP.A2019076.0230.001.2019076132237.nc
    Search_String = 'CLDMSK_L2_VIIRS_SNPP'//trim(Image%Level1b_Name(9:22))//'*.nc'

    Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0 .or. Num_Files > 1) then
        print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Not Found, "
        return
    endif

    Image%Auxiliary_Cloud_Mask_File_Name = Files(1)

    if (Seg_Idx == 1) &
       print *, EXE_PROMPT, MVCM_PROMPT, "SNPP NASA VIIRS MVCM File Found, ", &
            trim(Image%Auxiliary_Cloud_Mask_File_Name)

  endif

  !--- NEW version of NASA VIIRS Level1b NOAA-20
  if (index(Image%Level1b_Name,'VJ103MOD') == 1) then
    !0        1         2         3         4
    !12345678901234567890123456789012345678901234567890
    !VJ103MOD.A2014036.1154.001.2017300063048.uwssec.nc
    !VNPCLDMK.A2014036.1154.001.2017300132923.nc
    !CLDMSK_L2_VIIRS_NOAA20.A2019076.0230.001.2019076132237.nc
    Search_String = 'CLDMSK_L2_VIIRS_NOAA20'//trim(Image%Level1b_Name(9:22))//'*.nc'

    Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0 .or. Num_Files > 1) then
        print *, EXE_PROMPT, MVCM_PROMPT, "NOAA20 MVCM File Not Found, "
        return
    endif

    Image%Auxiliary_Cloud_Mask_File_Name = Files(1)

    if (Seg_Idx == 1) &
       print *, EXE_PROMPT, MVCM_PROMPT, "NOAA20 NASA VIIRS MVCM File Found, ", &
            trim(Image%Auxiliary_Cloud_Mask_File_Name)

  endif

  !--- NASA MODIS Level1b
  if (index(Image%Level1b_Name,'MYD021KM') == 1) then 

    !--- search should be the date and start time (ie. d20130426_t083000
    Search_String = 'IFFCMO_aqua_'//Image%Level1b_Name(11:27)//'*.hdf'
    Year_String = Image%Level1b_Name(11:14)
    Time_String = Image%Level1b_Name(19:22)
    Doy_String = Image%Level1b_Name(15:17)

    read(Doy_String,*) Doy
    read(Year_String,*) Year
    Ileap = LEAP_YEAR_FCT(Year)
    Month = COMPUTE_MONTH(Doy, ileap)
    Dom = COMPUTE_Day(Doy, ileap)
    write(Dom_String,fmt="(I2.2)") Dom
    write(Month_String,fmt="(I2.2)") Month

    Search_String = 'IFFCMO_aqua_d'//Year_String//Month_String//Dom_String//"_t"//Time_String//'*.hdf'

    Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0 .or. Num_Files > 1) then
        print *, EXE_PROMPT, MVCM_PROMPT, "Multiple NASA MODIS Level1b MVCM File Found, ",trim(Image%Auxiliary_Cloud_Mask_File_Name)
        return
    endif

    Image%Auxiliary_Cloud_Mask_File_Name = Files(1)

    if (Seg_Idx == 1) &
       print *, EXE_PROMPT, MVCM_PROMPT, "NASA MODIS Level1b MVCM File Found, ", &
            trim(Image%Auxiliary_Cloud_Mask_File_Name)


  endif

  !--- SIPS IFF VIIRS Level1b
  if (index(Image%Level1b_Name,'IFFSVM') == 1) then

    !--- search should be the date and start time (ie. d20130426_t083000
    !Search_String = 'IFFCMO_npp_'//Image%Level1b_Name(12:28)//'*.hdf'
    Search_String = 'VCLDMK_snpp_'//Image%Level1b_Name(12:28)//'*.nc'

    Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0 .or. Num_Files > 1) then
        print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Not Found, "
        return
    endif

    Image%Auxiliary_Cloud_Mask_File_Name = Files(1)

    if (Seg_Idx == 1) &
       print *, EXE_PROMPT, MVCM_PROMPT, "IFF-VIIRS Level1b MVCM File Found, ", &
            trim(Image%Auxiliary_Cloud_Mask_File_Name)

  endif

  !--- SIPS AHI GEO
  if (index(Image%Level1b_Name,'HS_H08') == 1) then
    !0        1         2         3         4
    !12345678901234567890123456789012345678901234567890
    !HS_H08_20160502_0200_B01_FLDK_R10_S0110.DAT
    !MVCMGEO.HIM08.20160502.0200_ahi_v1_0_3_REFERENCE_FD.nc
    Search_String = 'MVCMGEO.HIM08.'//Image%Level1b_Name(8:15)//'.'//Image%Level1b_Name(17:20)//'*.nc'

    Files => FILE_SEARCH(trim(Image%Level1b_Path),trim(Search_String),count=Num_Files)

    if (Num_Files == 0 .or. Num_Files > 1) then
        print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Not Found, "
        return
    endif

    Image%Auxiliary_Cloud_Mask_File_Name = Files(1)

    if (Seg_Idx == 1) &
       print *, EXE_PROMPT, MVCM_PROMPT, "SIPS HIM08 MVCM File Found, ", &
            trim(Image%Auxiliary_Cloud_Mask_File_Name)

  endif


  
  Files => null() 

 end subroutine DETERMINE_MVCM_NAME
    
!------------------------------------------------------------------------------------------
! open, read a slab from the MVCM file and close it
!
! slab refers to the data in the file
! seg refers to the data segment needed for clavr-x
!------------------------------------------------------------------------------------------
 subroutine READ_MVCM_DATA(Seg_Idx)
   
    integer, intent(in):: Seg_Idx
    integer:: Sd_Id, Group_Id, Status_Flag
    integer(kind=int4), dimension(2):: Sds_Dims, Sds_Stride, Sds_Start, Sds_Edges
    character(len=120):: Sds_Name
    integer(kind=int1), dimension(:,:), allocatable:: I1_Buffer
    integer:: Nx_Slab_Read_Start, Nx_Seg, Nx_Slab_Count
    integer:: Ny_Slab_Read_Start, Ny_Seg, Ny_Seg_Max, Ny_Slab_Count
  
    Status_Flag = 0
    Sds_Name = "Integer_Cloud_Mask"
    Cloud_Mask_Aux_Read_Flag = sym%NO   !will be set to yes if successful

    if (trim(Image%Auxiliary_Cloud_Mask_File_Name) == "no_file" .and. Seg_Idx == 1) then
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Could Not Be Found: "
    endif
    if (trim(Image%Auxiliary_Cloud_Mask_File_Name) == "no_file") then
         return
    endif

    !--- open file for read
    call OPEN_NETCDF(trim(Image%Level1b_Path)//trim(Image%Auxiliary_Cloud_Mask_File_Name), Sd_Id)

    !--- if file is unreadable, exit
    if (Sd_Id <= 0) then
         Status_Flag = 1
         print *, EXE_PROMPT, MVCM_PROMPT, "MVCM File Could Not Be Opened: "
         return
     endif

    !--- get information
    call GET_GROUP_ID(Sd_Id, "geophysical_data", Group_Id)

    ! --- get dimension of file
    call READ_NETCDF_DIMENSION_2D(Group_Id, trim(Sds_Name), Sds_Dims)

    ! --- determine expected size of data slab
    Nx_Seg = Image%Number_of_Elements
    Ny_Seg_Max = Image%Number_of_Lines_Per_Segment
    Ny_Seg = Image%Number_of_Lines_Read_This_Segment

    ! --- define size of data
    Nx_Slab_Read_Start = 1
    Ny_Slab_Read_Start = (Seg_Idx-1) * Ny_Seg_Max + 1
    Nx_Slab_Count = Sds_Dims(1)
    Ny_Slab_Count = Ny_Seg

    ! --- constrain to size of data
    if ((Seg_Idx * Ny_Seg) .gt. Sds_Dims(2)) Ny_Slab_Count = Sds_Dims(2) - Ny_Slab_Read_Start

    Sds_Stride = (/1, 1/)
    Sds_Start = (/Nx_Slab_Read_Start, Ny_Slab_Read_Start/)
    Sds_Edges = (/Nx_Slab_Count, Ny_Slab_Count/)

    ! --- allocate buffer
    if (.not. allocated(I1_Buffer)) allocate(I1_Buffer(Nx_Slab_Count, Ny_Slab_Count),stat=Status_Flag)

    ! --- read data
    call READ_NETCDF(Group_Id, Sds_Start, Sds_Stride, Sds_Edges, trim(Sds_Name), I1_Buffer)

    !--- close file
    call CLOSE_NETCDF(Sd_Id) 

    !--- switch to CLAVR-x convention for mask
    Cldmask%Cld_Mask_Aux(:,1:Ny_Slab_Count) = 3-I1_Buffer(:,1:Ny_Slab_Count)

    ! --- set missing to CLAVR-x missing
    where (Cldmask%Cld_Mask_Aux .lt. 0 .or. Cldmask%Cld_Mask_Aux .gt. 3)
       Cldmask%Cld_Mask_Aux = Missing_Value_Int1
    endwhere

    if (allocated(I1_Buffer)) deallocate(I1_Buffer,stat=Status_Flag)

    Cloud_Mask_Aux_Read_Flag = sym%YES
    Cloud_Type_Aux_Read_Flag = sym%NO

 end subroutine READ_MVCM_DATA

!-------------------------------------------
!
!-------------------------------------------
end module MVCM_READ_MOD
