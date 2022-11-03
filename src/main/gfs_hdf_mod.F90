!  $Id: gfs_hdf_mod.f90 4071 2021-01-20 05:16:22Z yli $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: gfs_hdf_module.f90 (src)
!       GFS (program)
!
! PURPOSE: This module houses all of the routines used to read and process the
!          GFS NWP data.  
!
! DESCRIPTION: The data used here are already in hdf format from the
!              convert_grib_to_hdf utility.  
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
! COPYRIGHT
!
! (c) This code is copyrighted by the author and all NOAA restrictions apply
!
! Dependencies:
!  CONSTANTS
!  NUMERICAL_ROUTINES
!  Sort_Module
!  NWP_COMMON_MOD
!  HDF
!
! Calling Sequence:
!  use GFS
!
! Public Routines within this module:
!  READ_GFS_DATA
!--------------------------------------------------------------------------------------
module GFS_HDF_MOD

  use CONSTANTS_MOD
 
  use CLASS_TIME_DATE, only: date_type , time_mid, time_diff_weight
  
  use CX_SCIENCE_TOOLS_MOD, only: wind_speed, wind_direction, vapor, vapor_ice
  
  use NWP_COMMON_MOD

  use CX_HDF4_MOD

  use CLAVRX_MESSAGE_MOD, only: mesg, verb_lev

  implicit none
  private
  public:: READ_GFS_DATA
 
  integer, parameter :: GFS_IDENT = 7

  real (kind=real4), parameter, private :: Missing_GFS = 9.999E+20

  character(len=11), parameter:: MODULE_PROMPT = "GFS_MODULE:"
  contains
  
  !   ---------
  !
  !  PURPOSE: returns both filenames for before and after
  !         and the time weight between both
  !
  !  INPUT:  
  !            t0, t1: date object from class_date_time for start and end of orbit
  !            nwp_type: 1 for GFS, 3 for CSFR, 4 for GDAS, 5 for MERRA, 6 for ERA and 7 for GFS_AIT
  !
  ! OUTPUT:    file1, file2 : filename below NWP root path ( with year directory if needed
  !              t_weight: time weight.. number between 0( Sensor time at file1 or 1 sensor time at file after 
  !
  !        Method considers the midpoint of start and end time
  !
  !  HISTORY:   31 January 2019: Introduced (A.W)
  !
  subroutine gfs_filenames  (t0,t1,nwp_type,file1,file2, t_weight)
    type(date_type), intent(in) :: t0, t1
    integer , intent(in) :: nwp_type
    character (len = 1020), intent(out) :: file1,file2
    real , intent(out) :: t_weight ! between 0 amd 1
    
    type(date_type) ::  t_mid, t_bef, t_aft
    character (len =2 ) :: fcs1, fcs2
    integer :: hour
    integer :: idx_3h
    integer :: shft
    logical :: bad_time_edges
    
    ! init
    shft = 0.  
      
    !- determine mid time
    t_mid = time_mid (t0,t1)
    ! call t_mid % print_data()
       
    ! - case for all nwp types
    select case (nwp_type)
       
    case(1)    !GFS
        shft = -2
        t_bef = t_mid % next_6h(shft-1)
        file1 = "gfs."//t_bef % date_string('yymmddhh')//"_F012.hdf"
        t_aft = t_mid % next_6h(shft)
        file2 = "gfs."//t_aft % date_string('yymmddhh')//"_F012.hdf"
        
    case(3)    !CFSR
        shft = -1
        t_bef = t_mid % next_6h(shft-1)
        
        
        
        file1 = t_bef% date_string('yyyy')//"/cfsr."//t_bef % date_string('yymmddhh')//"_F006.hdf"
        t_aft = t_mid % next_6h(shft)
        file2 = t_aft% date_string('yyyy')//"/cfsr."//t_aft % date_string('yymmddhh')//"_F006.hdf"
       
    case(4)    !GDAS
        shft = 0
        t_bef = t_mid % next_6h(shft-1)
        file1 = "gdas."//t_bef % date_string('yymmddhh')//"_F000.hdf"
        t_aft = t_mid % next_6h(shft)
        file2 = "gdas."//t_aft % date_string('yymmddhh')//"_F000.hdf"
       
       
    case(5)    !MERRA
        shft = 0
        t_bef = t_mid % next_6h(shft -1)
        file1 = t_bef% date_string('yyyy')//"/merra."//t_bef % date_string('yymmddhh')//"_F000.hdf"
        t_aft = t_mid % next_6h(shft)
        file2 = t_aft% date_string('yyyy')//"/merra."//t_aft % date_string('yymmddhh')//"_F000.hdf"
       
       
    case(6)   !ERA
        shft = 0
        t_bef = t_mid % next_6h(shft-1)
        file1 = t_bef% date_string('yyyy')//"/era."//t_bef % date_string('yymmddhh')//"_F000.hdf"
        t_aft = t_mid % next_6h(shft)
        file2 = t_aft% date_string('yyyy')//"/era."//t_aft % date_string('yymmddhh')//"_F000.hdf"
       
       
    case(7)  !GFS AIT
        ! here the F values are not fixed but a function of the 3h of the day box 
        ! for both, before and after we use the same start point of the NWP, but 
        ! with different forecast time
        !  
        shft = 0  
        call t_mid % get_date (hour = hour)
          
        idx_3h = hour/3
        select case (idx_3h)
           case(0,2,4,6)
              fcs1 = '06'
              fcs2 = '09'
               shft = -2
           case ( 1,3,5,7)
              fcs1 = '03'
              fcs2 = '06'
              shft = -1
        end select

        t_bef = t_mid % next_3h(shft-1)
        file1 = t_bef% date_string('yyyy')//"/gfs."//trim(t_bef % date_string(fmt='yymmddhh'))//"_F0"//fcs1//".hdf"
        file2 = t_bef% date_string('yyyy')//"/gfs."//t_bef % date_string('yymmddhh')//"_F0"//fcs2//".hdf"
        
        ! to calculate the weight we add  three hours
        t_aft = t_bef 
        call t_aft % add_time(hour=3)
        
        
       
         
    case(8)    !GFS FV3
        shft = -2
        t_bef = t_mid % next_6h(shft-1)
        file1 = "gfs."//t_bef % date_string('yymmddhh')//"_F012.hdf"
        t_aft = t_mid % next_6h(shft)
        file2 = "gfs."//t_aft % date_string('yymmddhh')//"_F012.hdf"
 
    case default

      print*,'wrong NWP option fatal error, should be fixed in clavrx_options stopping'
      print*,'nwp option is ',nwp_type
      stop    
       
    end select 
     
     
      
      t_weight  =  time_diff_weight (t_mid,t_bef,t_aft,identical_time_bounds  = bad_time_edges) + shft
     
    ! error catching for fatal errors which normally cannot occur
    if ( bad_time_edges) then
        print*,' GFS_HDF_MOD: time before and time after are identical'
        print*,' Unforeseen error'
        print*,' Please contact andi.walther@ssec.wisc.edu with all relevant information'
        stop 'stop gfs_mod.f90'
      
    end if
      
    if ( t_weight .lt. 0 .or. t_weight .gt. 1 ) then
          
          print*,'GFS_HDF_MOD:  time weight are not valid'
          print*,'This error shouldnt happen. Something went wrong for t_weight in gfs_hdf_mod.f90 tool gfs_filenames'
          print*, 'please contact Andi: andi.walther@ssec.wisc.edu'
          print*, 'time mid:'
          call t_mid % print_data
          print*
          print*,'time before:'
          call t_bef % print_data
          print*
          print*,'time after:'
          call t_aft % print_data
          
          
          
          stop 'stop gfs_mod.f90'
    end if
      
  end subroutine gfs_filenames

!-------------------------------------------------------------
! Subroutine to read in GFS data
!-------------------------------------------------------------
  subroutine READ_GFS_DATA(Nwp_Data_Type, &
                           Start_year, &
                           Start_jday, &
                           Start_Itime,&
                           end_year,   &
                           end_jday,   &
                           end_Itime,  &
                           Nwp_Path,   &
                           Ierror_Nwp)

!-------------------------------------------------------------
!   Local declarations
!-------------------------------------------------------------

    !   Arguments
    integer(kind=int4), intent(in) :: Nwp_Data_Type
    character (len=*), intent(in) :: Nwp_Path
    integer(kind=int2), intent(in) :: Start_year, Start_jday, end_year, end_jday
    integer(kind=int4), intent(in) :: Start_Itime, end_Itime
    integer, intent(out):: Ierror_Nwp
    integer :: year, jday

    !   Parameters

    !   Local variables
    character (len=1020) :: Nwp_Name_Before
    character (len=1020) :: Nwp_Name_After
    
    character (len=3)   :: array_order_1, array_order_2 
    
    logical :: file_exists
    
    integer :: i, j, ii, jj

    real ::t_weight

    integer:: Nlevels, Nlevels_o3mr, Nlevels_rh, Nlevels_clwmr, Nlat_Gfs, Nlon_Gfs
    real:: Dlat_Gfs, Dlon_Gfs, lat1_Gfs, lon1_Gfs

    real, dimension(:,:), allocatable :: Temp_2D_Real

    !-- values used for hdf reading
    integer:: Sd_Id_1,Sd_Id_2, Istatus
    integer, parameter:: Sds_Rank_1d = 1, Sds_Rank_2d = 2, Sds_Rank_3d = 3
    integer, dimension(Sds_Rank_1d):: Sds_Start_1d, Sds_Stride_1d, Sds_Edges_1d
    integer, dimension(Sds_Rank_2d):: Sds_Start_2d, Sds_Stride_2d, Sds_Edges_2d
    integer, dimension(Sds_Rank_3d):: Sds_Start_3d, Sds_Stride_3d, Sds_Edges_3d

    character(len=1024) :: string_msg
    
    type(date_type) :: t0, t1
    
    !-- check Nwp_Data_Type
    if (Nwp_Data_Type /= 1 .and. Nwp_Data_Type /= 3 .and. Nwp_Data_Type /=4 .and. Nwp_Data_Type /=5 .and. Nwp_Data_Type /=6 .and. Nwp_Data_Type /=7 &
       .and. Nwp_Data_Type /= 8) then
       print *, EXE_PROMPT, MODULE_PROMPT, " ERROR: unsupported NWP data in GFS read module, stopping"
       stop 
    endif

    !--   Initializations
    Ierror_Nwp = 0
    Nlat_Gfs   = -1
    Nlon_Gfs   = -1

    Missing_Nwp = Missing_Value_Real4
   
   
    call t0 % set_date_with_doy (year = int(start_year,kind(1)), doy = int(start_jday,kind(1)) &
        , hour =int ((Start_Itime)/60.0/60.0/1000.0 ) &
        , minute =  int (mod(int((Start_Itime)/60.0/1000.0), 60 ))  &
        , second = int (mod(int((Start_Itime)/1000.0), 60 )))
         
    call t1 % set_date_with_doy (year = int(end_year,kind(1)), doy = int(end_jday,kind(1)) &
          , hour = int ((End_Itime)/60.0/60.0/1000.0) &
        , minute =  int (mod(int((End_Itime)/60.0/1000.0), 60 ))  &
        , second = int (mod(int((End_Itime)/1000.0), 60 )))
    
    call gfs_filenames  (t0,t1,nwp_data_type,Nwp_Name_Before,Nwp_Name_after, t_weight)
    Nwp_Name_Before = trim(Nwp_Path) //"/"//trim(Nwp_Name_Before)
    Nwp_Name_After= trim(Nwp_Path) //"/"//trim(Nwp_Name_After)
    
    
    call MESG("NWP name before = "//trim(Nwp_Name_Before), level=Verb_Lev%DEFAULT)
    call MESG("NWP name after = "//trim(Nwp_Name_After), level=Verb_Lev%DEFAULT)
   
    !--- Does before file exist?
    inquire(file = Nwp_Name_Before, exist = file_exists)
    if (.not. file_exists) then
       call MESG("WARNING: Missing data file " //trim(Nwp_Name_Before),level = Verb_Lev%WARNING,color = 2)
       call MESG("WARNING: Setting GFS before to GFS after", level=Verb_Lev%WARNING,color=3)
       Nwp_Name_Before = Nwp_Name_After
       !stop 402
    endif

    !--- Does after file exist?
    inquire(file = Nwp_Name_After, exist = file_exists)
    if (.not. file_exists) then
       !--- check to see if before file was also Missing
       if (Nwp_Name_After == Nwp_Name_Before) then
          call mesg( ' '//MODULE_PROMPT//"ERROR: Neither GFS data files are present, processing stopped ",level=1,color=1)
          write(string_msg,*)"Stopped in File ",__FILE__," at line ",__LINE__
          call MESG ( string_msg,level = 1,color=3)
           stop 
       endif  
       call MESG("WARNING: Missing data file " // trim(Nwp_Name_After),level=Verb_Lev%WARNING)
       call MESG("WARNING: Setting GFS after to GFS before",level=Verb_Lev%WARNING)
       Nwp_Name_After = Nwp_Name_Before
    endif

    call MESG("GFS/CFSR NWP HDF files read in successfully",level=Verb_Lev%VERBOSE)
!--------------------------------------------------------------
! open before file and read contents from file data before this time
!--------------------------------------------------------------

     !- initialize
     Istatus = 0
   
     !-- open files
     Istatus = OPEN_FILE_HDF_READ(trim(Nwp_Name_Before), Sd_Id_1)
     if (Istatus < 0) then
      call MESG("ERROR: Failed open for read on first gfs file, stopping", level=Verb_Lev%ERROR)
      stop
     endif

     Istatus = OPEN_FILE_HDF_READ(trim(Nwp_Name_After), Sd_Id_2)
     if (Istatus < 0) then
      call MESG("ERROR: Failed open for read on second gfs file, stopping", level=Verb_Lev%ERROR)
      stop
     endif

     !---- read attributes from file 1
     Nlevels       = READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id_1,"NUMBER OF PRESSURE LEVELS")
     Nlevels_o3mr  = READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id_1,"NUMBER OF O3MR LEVELS")
     Nlevels_rh    = READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id_1,"NUMBER OF RH LEVELS")
     Nlevels_clwmr = READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id_1,"NUMBER OF CLWMR LEVELS")
     Nlat_Gfs      = READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id_1,"NUMBER OF LATITUDES")
     Nlon_Gfs      = READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id_1,"NUMBER OF LONGITUDES")
     Dlat_Gfs      = READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id_1,"LATITUDE RESOLUTION")
     Dlon_Gfs      = READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id_1,"LONGITUDE RESOLUTION")
     lat1_Gfs      = READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id_1,"FIRST LATITUDE")
     lon1_Gfs      = READ_HDF_GLOBAL_ATTRIBUTE_NUM(Sd_Id_1,"FIRST LONGITUDE") 
     array_order_1 = READ_HDF_GLOBAL_ATTRIBUTE_STR(Sd_Id_1,"3D ARRAY ORDER")

      !--- dont read attributes from file 2 (assume they are the same)
     array_order_2 = READ_HDF_GLOBAL_ATTRIBUTE_STR(Sd_Id_2,"3D ARRAY ORDER")

     !--- do some checking
     if (array_order_1 /= array_order_2) then 
       call MESG("ERROR: Order of 3d arrays differ among gfs files, stopping", level=Verb_Lev%ERROR)
       stop
     endif

     !--- set reformat gfs flag
     REFORMAT_GFS_ZXY = 0
     if (trim(array_order_1) == "XYZ" ) then
       REFORMAT_GFS_ZXY = 1
     endif

    !--- set up needed hdf dimensions

    !---- define dimensions of 1d arrays
    Sds_Start_1d = (/ 0 /)       !should this be 1
    Sds_Stride_1d = (/ 1 /)
    Sds_Edges_1d = (/ Nlevels /)     !should this be nx-1,ny-1
                                                                                                                           
    !---- define dimensions of 2d arrays
    Sds_Start_2d = (/ 0, 0 /)       !should this be 1
    Sds_Stride_2d = (/ 1, 1 /)
    Sds_Edges_2d = (/ Nlon_Gfs, Nlat_Gfs /)     !should this be nx-1,ny-1
                                                                                                                           
    !---- define dimensions of 3d arrays
    Sds_Start_3d = (/ 0, 0, 0 /)       
    Sds_Stride_3d = (/ 1, 1, 1 /)
    Sds_Edges_3d = (/ Nlevels, Nlon_Gfs, Nlat_Gfs /)     
    if (REFORMAT_GFS_ZXY == 1) then
      Sds_Edges_3d = (/ Nlon_Gfs, Nlat_Gfs, Nlevels /)    
    endif

    !---- store dimensions in public variables
    NWP%Nlevels = Nlevels
    NWP%Nlat = Nlat_Gfs
    NWP%Nlon = Nlon_Gfs

    !---- specific to GFS grid
    lat1_Nwp = lat1_Gfs
    lon1_Nwp = lon1_Gfs
    Dlon_Nwp = Dlon_Gfs
    Dlat_Nwp = sign(Dlat_Gfs,-1.0*lat1_Gfs)

    !-----------------------------------------------------------------
    ! allocate NWP arrays
    !-----------------------------------------------------------------
    call CREATE_NWP_ARRAYS()
    call INITIALIZE_NWP_ARRAYS()

!------------------------------------------------------------------
!   Read data from the first file
!------------------------------------------------------------------

!--- read in standard levels from first file
    Istatus = 0
    Istatus = HDF_SDS_READER(Sd_Id_1,"pressure levels", Sds_Start_1d, Sds_Stride_1d, &
                      Sds_Edges_1d, NWP%P_Std) + Istatus

!--- read in two dimensional arrays

!- surface temperature
    Istatus = HDF_SDS_READER(Sd_Id_1,"surface temperature",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"surface temperature",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Tmpsfc)

!- surface height
    Istatus = HDF_SDS_READER(Sd_Id_1,"surface height",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"surface height",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Zsfc)

    !--- check if missing data is critical to further processing
    if (Istatus < 0) then
       call MESG(EXE_PROMPT// MODULE_PROMPT//"ERROR: critical 2d nwp sds missing, program stopping" &
                 , level=Verb_Lev%ERROR)
    endif

!- land mask
    allocate (Temp_2D_Real(Sds_Edges_2d(1),Sds_Edges_2d(2)))
    Istatus = HDF_SDS_READER(Sd_Id_1,"land mask", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"land mask", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, Temp_2D_Real)
    NWP%Land = int(Temp_2D_Real, kind=int1)
    deallocate (Temp_2D_Real)


!- ice mask
    Istatus = HDF_SDS_READER(Sd_Id_1,"ice fraction",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"ice fraction",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Sea_Ice_Frac)

!- msl pressure
    Istatus = HDF_SDS_READER(Sd_Id_1,"MSL pressure",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"MSL pressure",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Pmsl)

!- surface pressure
    Istatus = HDF_SDS_READER(Sd_Id_1,"surface pressure",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"surface pressure",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Psfc)

!- air temperature
    Istatus = HDF_SDS_READER(Sd_Id_1,"temperature at sigma=0.995",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"temperature at sigma=0.995",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Tmpair)

!- air humidity
    Istatus = HDF_SDS_READER(Sd_Id_1,"rh at sigma=0.995",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"rh at sigma=0.995",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Rhsfc)

!- water equivalent snow depth
    Istatus = HDF_SDS_READER(Sd_Id_1,"water equivalent snow depth",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"water equivalent snow depth",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Weasd)

!- TPW
    Istatus = HDF_SDS_READER(Sd_Id_1,"total precipitable water", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"total precipitable water", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Tpw)

!- tropopause temperature
    Istatus = HDF_SDS_READER(Sd_Id_1,"tropopause temperature", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"tropopause temperature", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%T_Trop)

!- tropopause pressure
    Istatus = HDF_SDS_READER(Sd_Id_1,"tropopause pressure", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"tropopause pressure", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%P_Trop)

!- wind direction
    Istatus = HDF_SDS_READER(Sd_Id_1,"u-wind at sigma=0.995", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"u-wind at sigma=0.995", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%U_Wnd_10m)

!- wind speed
    Istatus = HDF_SDS_READER(Sd_Id_1,"v-wind at sigma=0.995", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"v-wind at sigma=0.995", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%V_Wnd_10m)


!- total ozone
    Istatus = HDF_SDS_READER(Sd_Id_1,"total ozone", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"total ozone", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Ozone)


!--- read in three dimensional arrays
!- temperature
    Istatus = 0
    Istatus = HDF_SDS_READER(Sd_Id_1,"temperature", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"temperature", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_2) + Istatus
    call INTERPOLATE_3D(t_weight, Missing_Nwp, Missing_GFS, NWP%T_Prof)

    !- height
    Istatus = HDF_SDS_READER(Sd_Id_1,"height", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"height", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_2) + Istatus
    call INTERPOLATE_3D(t_weight, Missing_Nwp, Missing_GFS, NWP%Z_Prof)

    !- humidity
    Istatus = HDF_SDS_READER(Sd_Id_1,"rh", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"rh", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_2) + Istatus
    call INTERPOLATE_3D(t_weight, Missing_Nwp, Missing_GFS, NWP%Rh_Prof)

    !--- check if missing data is critical to further processing
    if (Istatus < 0) then
       call MESG(EXE_PROMPT// MODULE_PROMPT//"ERROR: critical 3d nwp sds missing, program stopping" &
                 , level=Verb_Lev%ERROR)
    endif

    !- wind direction
    Istatus = HDF_SDS_READER(Sd_Id_1,"u-wind", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"u-wind", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_2) + Istatus
    call INTERPOLATE_3D(t_weight, Missing_Nwp, Missing_GFS, NWP%U_Wnd_Prof)

    !- wind speed
    Istatus = HDF_SDS_READER(Sd_Id_1,"v-wind", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"v-wind", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_2) + Istatus
    call INTERPOLATE_3D(t_weight, Missing_Nwp, Missing_GFS, NWP%V_Wnd_Prof)

    !- ozone
    Istatus = HDF_SDS_READER(Sd_Id_1,"o3mr", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"o3mr", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_2) + Istatus
    call INTERPOLATE_3D(t_weight, Missing_Nwp, Missing_GFS, NWP%Ozone_Prof)

    !- cloud liquid water mixing ratio
    Istatus = HDF_SDS_READER(Sd_Id_1,"clwmr", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"clwmr", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_2) + Istatus
    call INTERPOLATE_3D(t_weight, Missing_Nwp, Missing_GFS, NWP%Clwmr_Prof)

    !- cape
  if (Nwp_Data_Type == 1) then
    Istatus = HDF_SDS_READER(Sd_Id_1,"cape", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"cape", Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%CAPE)
  endif

  if (Nwp_Data_Type == 8) then
    !- freezing level height
    Istatus = HDF_SDS_READER(Sd_Id_1,"freezing level",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"freezing level",Sds_Start_2d, Sds_Stride_2d, &
                      Sds_Edges_2d, Temp2d_Nwp_2) + Istatus
    call INTERPOLATE_2D(t_weight, Missing_Nwp, Missing_GFS, NWP%Freezing_lev_hgt)

    !- ic cloud mixing ratio
    Istatus = HDF_SDS_READER(Sd_Id_1,"icmr", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_1) + Istatus
    Istatus = HDF_SDS_READER(Sd_Id_2,"icmr", Sds_Start_3d, Sds_Stride_3d, &
                      Sds_Edges_3d, Temp3d_Nwp_2) + Istatus
    call INTERPOLATE_3D(t_weight, Missing_Nwp, Missing_GFS, NWP%Icmr_Prof)
  endif
 
    ! --- close files
    Istatus = CLOSE_FILE_HDF_READ(Sd_Id_1,trim(Nwp_Name_Before))
    Istatus = CLOSE_FILE_HDF_READ(Sd_Id_2,trim(Nwp_Name_After))

!----------------------------------------------------------------------------------
! modify gfs fields as warranted for use in CLAVR-x
!----------------------------------------------------------------------------------

!--- adopt a standard missing value
if (Missing_Nwp /= Missing_GFS) then
 call FIX_MISSING_2D(NWP%Sea_Ice_Frac, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%Pmsl, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%Psfc, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%Tmpsfc, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%Zsfc, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%Tmpair, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%Rhsfc, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%Weasd, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%Tpw, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%T_Trop, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%P_Trop, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%U_wnd_10m, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_2D(NWP%V_wnd_10m, Missing_GFS, Missing_Nwp)

 call FIX_MISSING_3D(NWP%T_Prof, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_3D(NWP%Z_Prof, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_3D(NWP%U_Wnd_Prof, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_3D(NWP%Ozone_Prof, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_3D(NWP%Rh_Prof, Missing_GFS, Missing_Nwp)
 call FIX_MISSING_3D(NWP%Clwmr_Prof, Missing_GFS, Missing_Nwp)
 if (Nwp_Data_Type == 1) then
    call FIX_MISSING_2D(NWP%CAPE, Missing_GFS, Missing_Nwp)
 endif
 if (Nwp_Data_Type == 8) then
    call FIX_MISSING_2D(NWP%Freezing_lev_hgt, Missing_GFS, Missing_Nwp)
    call FIX_MISSING_3D(NWP%Icmr_Prof, Missing_GFS, Missing_Nwp)
 endif
endif

!--- constraint mixing ratios to be non-negative
where(NWP%Rh_Prof /= Missing_Nwp .and. NWP%Rh_Prof < 0.00)
      NWP%Rh_Prof = 0.00
endwhere

where(NWP%Clwmr_Prof /= Missing_Nwp .and. NWP%Clwmr_Prof < 0.00)
      NWP%Clwmr_Prof = 0.00
endwhere

where(NWP%Ozone_Prof /= Missing_Nwp .and. NWP%Ozone_Prof < 0.00)
      NWP%Ozone_Prof = 0.00
endwhere

where(NWP%Icmr_Prof /= Missing_Nwp .and. NWP%Icmr_Prof < 0.00)
      NWP%Icmr_Prof = 0.00
endwhere

!--- fix GFS bug in RH
call FIX_GFS_RH()

!--- convert ozone from mass mixing ratio(g/g) to volume missing ratio (ppmv)
where(NWP%Ozone_Prof > 0)
    NWP%Ozone_Prof = 1.0e06*NWP%Ozone_Prof * 0.602
endwhere

!--- Convert NWP%Zsfc to meters
where (NWP%Zsfc /= Missing_Nwp)
  NWP%Zsfc = NWP%Zsfc * 1000.0
end where

!--- Convert NWP%Z_Prof to meters
where (NWP%Z_Prof /= Missing_Nwp)
  NWP%Z_Prof = NWP%Z_Prof * 1000.0
end where

!--- Convert NWP%Freezing_lev_hgt to meters
where (NWP%Freezing_lev_hgt /= Missing_Nwp)
  NWP%Freezing_lev_hgt = NWP%Freezing_lev_hgt * 1000.
end where

!---- compute wind speed
NWP%Wnd_Spd_10m = Wind_Speed(NWP%U_Wnd_10m,NWP%V_Wnd_10m)
NWP%Wnd_Dir_10m = Wind_Direction(NWP%U_Wnd_10m,NWP%V_Wnd_10m)

!----- check for Missing t tropopause values and set to 200 K
do i = 1, NWP%Nlon
 do j = 1, NWP%Nlat
    if ((NWP%T_Trop(i,j) < 180.0).or.(NWP%T_Trop(i,j) > 240.0)) then
         NWP%T_Trop(i,j) = 200.0
    endif
 enddo
enddo

 end subroutine READ_GFS_DATA

!------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------
 subroutine FIX_MISSING_2D(Data, Missing_In, Missing_Out)
    real, dimension(:,:), intent(inout):: Data
    real, intent(in):: Missing_In
    real, intent(in):: Missing_Out

    where(Data == Missing_In) 
        Data = Missing_Out
    endwhere

 end subroutine FIX_MISSING_2D

 subroutine FIX_MISSING_3D(Data, Missing_In, Missing_Out)
    real, dimension(:,:,:), intent(inout):: Data
    real, intent(in):: Missing_In
    real, intent(in):: Missing_Out

    where(Data == Missing_In) 
        Data = Missing_Out
    endwhere

 end subroutine FIX_MISSING_3D

!--------------------------------------------------------------------------------------------
! generic HDF4 read routines
!--------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! interpolate 2d dimensional data from the GFS HDF files
!------------------------------------------------------------------------------------------
subroutine INTERPOLATE_2D(X,Missing1,Missing2,Temp)
  real, intent(in):: X, Missing1, Missing2
  real, dimension(:,:), intent(out):: Temp

  integer:: i1,i2,Ndim1,Ndim2

    Temp = Missing1
    Ndim1 = size(Temp,1)
    Ndim2 = size(Temp,2)
    do i1 = 1, Ndim1
      do i2 = 1, Ndim2
          if ((Temp2d_Nwp_1(i1,i2) /= Missing1) .and. (Temp2d_Nwp_2(i1,i2) /= Missing1) .and. &
              (Temp2d_Nwp_1(i1,i2) /= Missing2) .and. (Temp2d_Nwp_2(i1,i2) /= Missing2)) then
           Temp(i1,i2) = (1.0 - X) * Temp2d_Nwp_1(i1,i2)   + X * Temp2d_Nwp_2(i1,i2)
          endif
      enddo
    enddo

    return

end subroutine INTERPOLATE_2D

!------------------------------------------------------------------------------------------
! interpolate 3d dimensional data from the GFS HDF files
!------------------------------------------------------------------------------------------
subroutine INTERPOLATE_3D(X,Missing1,Missing2,Temp)
  real, intent(in):: X, Missing1, Missing2
  real, dimension(:,:,:), intent(out):: Temp

  integer:: i1,i2,i3,Ndim1,Ndim2,Ndim3

    Temp = Missing1
    Ndim1 = size(Temp,1)
    Ndim2 = size(Temp,2)
    Ndim3 = size(Temp,3)

    if (REFORMAT_GFS_ZXY == 1) then
      do i1 = 1,Ndim1
         Temp3d(i1,:,:) = Temp3d_Nwp_1(:,:,i1)
      enddo
      Temp3d_Nwp_1 = Temp3d
      do i1 = 1,Ndim1
         Temp3d(i1,:,:) = Temp3d_Nwp_2(:,:,i1)
      enddo
      Temp3d_Nwp_2 = Temp3d
    endif

!-- note this where state is seg-faulting (replaced with equivalent do/if construct)
!   where ((temp3d_Nwp_1 == Missing) .or. (temp3d_Nwp_2 == Missing))
!      temp   = Missing_Nwp
!   elsewhere
!      temp   = (1.0 - x) * temp3d_Nwp_1   + x * temp3d_Nwp_2
!   endwhere
!----------------------------------------------------------------

    do i1 = 1, Ndim1
      do i2 = 1, Ndim2
        do i3 = 1, Ndim3
          if ((Temp3d_Nwp_1(i1,i2,i3) /= Missing1) .and. (Temp3d_Nwp_2(i1,i2,i3) /= Missing1) .and. &
              (Temp3d_Nwp_1(i1,i2,i3) /= Missing2) .and. (Temp3d_Nwp_2(i1,i2,i3) /= Missing2)) then
            Temp(i1,i2,i3) = (1.0 - X) * Temp3d_Nwp_1(i1,i2,i3)   + X * Temp3d_Nwp_2(i1,i2,i3)
          endif
        enddo
      enddo
    enddo

    return

end subroutine INTERPOLATE_3D

!---------------------------------------------------------------------
! Fix GFS RH scaling
!
! In the current GFS output, the definition of RH varies between
! 253 and 273 K.  At 273 it is with respect to water. 
! At 253 it is defined with respect to ice.  
! At temperatures in between, it varies linearly.
! This routine attempts to define RH with respect to water for all temps
!---------------------------------------------------------------------
subroutine FIX_GFS_RH()

   integer:: i, j, k
   real:: es_water
   real:: es_ice
   real:: e
   real:: es
   real:: ice_weight

   do i = 1, NWP%Nlon
      do j = 1, NWP%Nlat

       do k = 1, NWP%Nlevels
          if (NWP%Rh_Prof(k,i,j) > 0.0) then

           !--- compute saturation vapor pressures
           es_water = VAPOR(NWP%T_Prof(k,i,j))
           es_Ice = VAPOR_ICE(NWP%T_Prof(k,i,j))

           !--- derive the ice/water weight used in gfs
           ice_weight = (273.16 - NWP%T_Prof(k,i,j)) / &
                       (273.16-253.16)
           ice_weight = min(1.0,max(0.0,ice_weight))

           !--- derive es used in original rh definition
           es = ice_weight * es_ice + (1.0-ice_weight)*es_water

           !--- compute actual e 
           e = NWP%Rh_Prof(k,i,j) * es / 100.0

           !--- compute actual rh with respect to water
           NWP%Rh_Prof(k,i,j) = 100.0 * e / es_water

          endif
       enddo

      enddo
   enddo

end subroutine FIX_GFS_RH


end module GFS_HDF_MOD
