!$Id: land_sfc_properties_mod.f90 3189 2019-03-11 22:04:07Z heidinger $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: land_sfc_properties.f90 (src)
!       land_sfc_properties (program)
!
! PURPOSE: 
!
! DESCRIPTION: This module taken from GEOCAT and modified for CLAVR-x
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
!--------------------------------------------------------------------------------------
module LAND_SFC_PROPERTIES_MOD
  use CONSTANTS_MOD,only: &
    int4,real8, int1, int2, real4 &
    , MISSING_VALUE_REAL8 &
    , sym &
    , missing_value_int1 &
    , missing_value_int2
    
  use CLASS_TIME_DATE,only: &
    date_type
    
  use CX_HDF4_MOD, only: &
   READ_HDF_GLOBAL_ATTRIBUTE_NUM &
   , HDF_SDS_DIMENSIONS_READER &
   , CLOSE_FILE_HDF_READ &
   , HDF_SDS_READER
  
  use NUMERICAL_ROUTINES_MOD,only: &
   find_bounds
  use FILE_TOOLS,only: &
   file_test
  use CLAVRX_MESSAGE_MOD, only: MESG,MESG_1I,VERB_LEV
  use HDF, only: MAX_RANK_HDF, DFACC_READ, FAIL
  
  implicit none
  private
  public  :: read_land_sfc_hdf
  public  :: get_snow_map_filename
  public :: OPEN_LAND_SFC_HDF
  public :: close_land_sfc_hdf

  interface read_land_sfc_hdf
        module procedure  &
           read_land_sfc_hdf_i1,  &
           read_land_sfc_hdf_i2
   end interface
  
  
  
  REAL(kind=real8), parameter, private :: FIRST_LAT_DEFAULT = -90.0_real8
  REAL(kind=real8), parameter, private :: LAST_LAT_DEFAULT = 90.0_real8
  REAL(kind=real8), parameter, private :: FIRST_LON_DEFAULT = -180.0_real8
  REAL(kind=real8), parameter, private :: LAST_LON_DEFAULT = 180.0_real8
  REAL(kind=real8), parameter, private :: DEL_LAT_DEFAULT = 0.04_real8
  REAL(kind=real8), parameter, private :: DEL_LON_DEFAULT = 0.04_real8

 

  TYPE, public :: land_grid_description
    CHARACTER(len=256) :: sds_name
    INTEGER(kind=int4) :: num_lat
    INTEGER(kind=int4) :: num_lon
    REAL(kind=real8) :: del_lat
    REAL(kind=real8) :: del_lon
    REAL(kind=real8) :: first_lat
    REAL(kind=real8) :: first_lon
  end TYPE land_grid_description
    
  CONTAINS
  
  !-------------------------------------------------------------------
  ! Subroutine to open the land surface file.
  !-------------------------------------------------------------------
  FUNCTION open_land_sfc_hdf(data_dir, filename, grid_str) result(id)
    CHARACTER(len=*), intent(in) :: data_dir, filename
    TYPE(land_grid_description), optional, intent(inout) :: grid_str
  
    INTEGER(kind=int4) :: id  
    CHARACTER(len=1020) :: filename_full
  
    logical :: file_exists
  
    INTEGER :: sfstart

    integer:: Rank
    integer, dimension(MAX_RANK_HDF):: Dims

    integer:: Istatus
  
    filename_full = trim(data_dir)//trim(filename)
  
    inquire(file = filename_full, exist = file_exists)
    if (.not. file_exists) then
      call mesg("Land surface file "//trim(filename_full)//" does not exist. ",level=verb_lev % ERROR)
      stop
    end if
  
    id = sfstart(trim(filename_full), DFACC_READ)
    if (id == FAIL) then
      call mesg("Failed to open "//trim(filename_full),level=verb_lev % ERROR)
      stop
    end if
  
    if (present(grid_str)) then
      !--- note these attributes are real8 but current cx_hdf4_mod.f90 does
      !--- not support real 8 attributes.  Is this a problem?
      grid_str%del_lat = real(READ_HDF_GLOBAL_ATTRIBUTE_NUM(id, "dlat"),kind=real8)
      grid_str%del_lon = real(READ_HDF_GLOBAL_ATTRIBUTE_NUM(id, "dlon"),kind=real8) 
      grid_str%first_lat = real(READ_HDF_GLOBAL_ATTRIBUTE_NUM(id, "first_lat"),kind=real8) 
      grid_str%first_lon = real(READ_HDF_GLOBAL_ATTRIBUTE_NUM(id, "first_lon"),kind=real8)

      !--- round - I added this to avoid noise from real4 to real8 conversion
      grid_str%del_lat = real(nint(grid_str%del_lat * 1.0e06) / 1.0e06,kind=real8)
      grid_str%del_lon = real(nint(grid_str%del_lon * 1.0e06) / 1.0e06,kind=real8)
      grid_str%first_lat = real(nint(grid_str%first_lat * 1.0e06) / 1.0e06,kind=real8)
      grid_str%first_lon = real(nint(grid_str%first_lon * 1.0e06) / 1.0e06,kind=real8)

      if (grid_str%del_lat == MISSING_VALUE_REAL8) grid_str%del_lat = DEL_LAT_DEFAULT
      if (grid_str%del_lon == MISSING_VALUE_REAL8) grid_str%del_lon = DEL_LON_DEFAULT
      if (grid_str%first_lat == MISSING_VALUE_REAL8) grid_str%first_lat = FIRST_LAT_DEFAULT
      if (grid_str%first_lon == MISSING_VALUE_REAL8) grid_str%first_lon = FIRST_LON_DEFAULT
  
      !   call read_hdf_sds_dimenions(id, grid_str%sds_name, grid_str%num_lat, grid_str%num_lon) 

      Istatus = HDF_SDS_DIMENSIONS_READER(id, trim(grid_str%sds_name), Rank, Dims)
      grid_str%num_lon = Dims(1)
      grid_str%num_lat = Dims(2)

    end if
  
    return

  end FUNCTION open_land_sfc_hdf

  !-------------------------------------------------------------------
  ! Subroutine to close the land surface file.
  !-------------------------------------------------------------------
  subroutine close_land_sfc_hdf(id)
    integer(kind=int4), intent(inout) :: id
    integer(kind=int4) :: istatus  
    integer :: sfend
  
    ! istatus = sfend(id)
    istatus = CLOSE_FILE_HDF_READ(id,'land_sfc_file_name')
    if (istatus /= 0) then
      call mesg("Error closing land surface hdf file. ",level=verb_lev % ERROR)
      stop
    end if

    id = -999

  end subroutine close_land_sfc_hdf


  !-------------------------------------------------------------------
  ! Function to find the snow map name.
  !-------------------------------------------------------------------

  FUNCTION get_snow_map_filename(year_in,day_of_year,snow_path) result(snow_filename)
    CHARACTER(*), intent(in) :: snow_path
    INTEGER(kind=int2), intent(in):: year_in
    INTEGER(kind=int2), intent(in):: day_of_year
    CHARACTER(len=1020) :: snow_filename
    CHARACTER(len=1020) :: snow_filename_tmp
    INTEGER(kind=int4) :: iday
    TYPE(date_type) :: time_obj
    INTEGER(kind=int4), parameter :: MAX_SNOW_LATENCY = 4 !including current day

    snow_filename = "no_file"

    do iday=0, MAX_SNOW_LATENCY - 1
      call time_obj %set_date_with_doy( int(year_in,kind=int4), day_of_year - iday)
      snow_filename_tmp = "snow_map_4km_"//time_obj%date_string('yymmdd')//".hdf"
      if (file_test(trim(snow_path)//trim(snow_filename_tmp))) then
        snow_filename = snow_filename_tmp
        exit
      end if
    end do
    return

  end FUNCTION get_snow_map_filename
  
 



  !-------------------------------------------------------------------
  ! Subroutine to read a given land surface hdf file.
  !-------------------------------------------------------------------
  subroutine read_land_sfc_hdf_i1(id, grid_str, lat, lon, space_mask, land)
    INTEGER(kind=int4), intent(in) :: id
    TYPE(land_grid_description), intent(in) :: grid_str
    REAL(kind=real4), dimension(:,:), intent(in) :: lat, lon
    INTEGER(kind=int1), dimension(:,:), intent(in) :: space_mask
    INTEGER(kind=int1), dimension(:,:), intent(out) :: land
    
    INTEGER :: astatus
    INTEGER :: ilat1, ilat2, ilon1, ilon2, ilat, ilon, ilat_ad, ilon_ad, &
             ilon1_2, ilon2_2
    INTEGER :: temp, nx, ny, i, j
    INTEGER(kind=int1), dimension(:,:), allocatable :: land_grid, land_grid_2
    REAL(kind=real4) :: wlon, elon, slat, nlat
    INTEGER(kind=int1) :: dateline_flg, space_check
  
    INTEGER, dimension(2) :: start_2d, stride_2d, edge_2d, &
                           start_2d_2, stride_2d_2, edge_2d_2
 
    INTEGER:: Istatus 

    space_check = minval(space_mask)
    if (space_check == 1) then
      land = missing_value_int1
      return
    end if
  
    nx = size(lat,1)
    ny = size(lat,2)
  
    call find_bounds(lat,lon,wlon,elon,slat,nlat,dateline_flg)

    if (dateline_flg == 0) then
    
      ilat1 = max(0,min(grid_str%num_lat,int(abs(nlat - grid_str%first_lat)/grid_str%del_lat) + 0))
      ilat2 = max(0,min(grid_str%num_lat,int(abs(slat - grid_str%first_lat)/grid_str%del_lat) + 0))
  
      ilon1 = max(0,min(grid_str%num_lon,int(abs(wlon - grid_str%first_lon)/grid_str%del_lon) + 0))
      ilon2 = max(0,min(grid_str%num_lon,int(abs(elon - grid_str%first_lon)/grid_str%del_lon) + 0))
    
      if (ilat1 > ilat2) then
        temp = ilat1
        ilat1 = ilat2
        ilat2 = temp
      end if
  
      if (ilon1 > ilon2) then
        temp = ilon1
        ilon1 = ilon2
        ilon2 = temp
      end if
  
      start_2d = (/ilon1, ilat1/)
      stride_2d = (/1, 1/)
      edge_2d = (/(ilon2-ilon1)+1, (ilat2-ilat1)+1/)
  
      !call read_hdf_sds(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
      Istatus = HDF_SDS_READER(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)

      do j = 1, ny
        do i = 1, nx
    
          if (space_mask(i,j) == sym%NO_SPACE) then
                
            ilat = max(1,min(grid_str%num_lat,int(abs(lat(i,j) - grid_str%first_lat)/grid_str%del_lat) + 1))
            ilon = max(1,min(grid_str%num_lon,int(abs(lon(i,j) - grid_str%first_lon)/grid_str%del_lon) + 1))
            ilat_ad = max(1,min((ilat - start_2d(2)) + 1,size(land_grid,2)))
            ilon_ad = max(1,min((ilon - start_2d(1)) + 1,size(land_grid,1)))
            land(i,j) = land_grid(ilon_ad,ilat_ad)
  
          end if
      
        end do
      end do
    
      deallocate(land_grid, stat=astatus)
      if (astatus /= 0) then
        call mesg("Error deallocating land surface grid. ",level=verb_lev % ERROR)
        stop
      end if
    
    else ! dateline flag ne 0
  
      ilat1 = max(0,min(grid_str%num_lat,int(abs(nlat - grid_str%first_lat)/grid_str%del_lat) + 0))
      ilat2 = max(0,min(grid_str%num_lat,int(abs(slat - grid_str%first_lat)/grid_str%del_lat) + 0))
  
      ilon1 = max(0,min(grid_str%num_lon,int(abs(wlon - grid_str%first_lon)/grid_str%del_lon) + 0))
      ilon2 = max(0,min(grid_str%num_lon,int(abs(180.0 - grid_str%first_lon)/grid_str%del_lon) + 0))
    
      ilon1_2 = max(0,min(grid_str%num_lon,int(abs(-180.0 - grid_str%first_lon)/grid_str%del_lon) + 0))
      ilon2_2 = max(0,min(grid_str%num_lon,int(abs((elon-360.0) - grid_str%first_lon)/grid_str%del_lon) + 0))

      if (ilat1 > ilat2) then
        temp = ilat1
        ilat1 = ilat2
        ilat2 = temp
      end if
  
      if (ilon1 > ilon2) then
        temp = ilon1
        ilon1 = ilon2
        ilon2 = temp
      end if
    
      if (ilon1_2 > ilon2_2) then
        temp = ilon1_2
        ilon1_2 = ilon2_2
        ilon2_2 = temp
      end if
  
      start_2d = (/ilon1, ilat1/)
      stride_2d = (/1, 1/)
      edge_2d = (/(ilon2-ilon1)+1, (ilat2-ilat1)+1/)
  
      !call read_hdf_sds(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
      Istatus = HDF_SDS_READER(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
    
      start_2d_2 = (/ilon1_2, ilat1/)
      stride_2d_2 = (/1, 1/)
      edge_2d_2 = (/(ilon2_2-ilon1_2)+1, (ilat2-ilat1)+1/)
  
      !call read_hdf_sds(id, trim(grid_str%sds_name), start_2d_2, stride_2d_2, edge_2d_2, land_grid_2)
      Istatus = HDF_SDS_READER(id, trim(grid_str%sds_name), start_2d_2, stride_2d_2, edge_2d_2, land_grid_2)
    
      do j = 1, ny
        do i = 1, nx
    
          if (space_mask(i,j) == sym%NO_SPACE) then

            ilat = max(1,min(grid_str%num_lat,int(abs(lat(i,j) - grid_str%first_lat)/grid_str%del_lat) + 1))
            ilon = max(1,min(grid_str%num_lon,int(abs(lon(i,j) - grid_str%first_lon)/grid_str%del_lon) + 1))

            if (lon(i,j) >= 0.0) then
              ilat_ad = max(1,min((ilat - start_2d(2)) + 1,size(land_grid,2)))
              ilon_ad = max(1,min((ilon - start_2d(1)) + 1,size(land_grid,1)))
              land(i,j) = land_grid(ilon_ad,ilat_ad)
            else
              ilat_ad = max(1,min((ilat - start_2d_2(2)) + 1,size(land_grid_2,2)))
              ilon_ad = max(1,min((ilon - start_2d_2(1)) + 1,size(land_grid_2,1)))
              land(i,j) = land_grid_2(ilon_ad,ilat_ad)
            end if

          end if
      
        end do
      end do

      deallocate(land_grid, land_grid_2, stat=astatus)
      if (astatus /= 0) then
        call mesg("Error deallocating land surface grid. ",level=verb_lev % ERROR)
        stop
      end if
    
    end if
  
  end subroutine read_land_sfc_hdf_i1

  !-------------------------------------------------------------------
  ! Subroutine to read a given land surface hdf file.
  !-------------------------------------------------------------------
  subroutine read_land_sfc_hdf_i2(id, grid_str, lat, lon, space_mask, land)
    INTEGER(kind=int4), intent(in) :: id
    TYPE(land_grid_description), intent(in) :: grid_str
    REAL(kind=real4), dimension(:,:), intent(in) :: lat, lon
    INTEGER(kind=int1), dimension(:,:), intent(in) :: space_mask
    INTEGER(kind=int2), dimension(:,:), intent(out) :: land
    
    INTEGER :: astatus
    INTEGER :: ilat1, ilat2, ilon1, ilon2, ilat, ilon, ilat_ad, ilon_ad, &
             ilon1_2, ilon2_2
    INTEGER :: temp, nx, ny, i, j
    INTEGER(kind=int2), dimension(:,:), allocatable :: land_grid, land_grid_2
    REAL(kind=real4) :: wlon, elon, slat, nlat
    INTEGER(kind=int1) :: dateline_flg, space_check
  
    INTEGER, dimension(2) :: start_2d, stride_2d, edge_2d, &
                           start_2d_2, stride_2d_2, edge_2d_2
    INTEGER:: Istatus

    ! - executable
    space_check = minval(space_mask)
    if (space_check == 1) then
      land = missing_value_int2
      return
    end if
 
    nx = size(lat,1)
    ny = size(lat,2)

    call find_bounds(lat,lon,wlon,elon,slat,nlat,dateline_flg)

    if (dateline_flg == 0) then
    
      ilat1 = max(0,min(grid_str%num_lat,int(abs(nlat - grid_str%first_lat)/grid_str%del_lat) + 0))
      ilat2 = max(0,min(grid_str%num_lat,int(abs(slat - grid_str%first_lat)/grid_str%del_lat) + 0))
  
      ilon1 = max(0,min(grid_str%num_lon,int(abs(wlon - grid_str%first_lon)/grid_str%del_lon) + 0))
      ilon2 = max(0,min(grid_str%num_lon,int(abs(elon - grid_str%first_lon)/grid_str%del_lon) + 0))
    
      if (ilat1 > ilat2) then
        temp = ilat1
        ilat1 = ilat2
        ilat2 = temp
      end if
  
      if (ilon1 > ilon2) then
        temp = ilon1
        ilon1 = ilon2
        ilon2 = temp
      end if
  
      start_2d = (/ilon1, ilat1/)
      stride_2d = (/1, 1/)
      edge_2d = (/(ilon2-ilon1)+1, (ilat2-ilat1)+1/)

      !call read_hdf_sds(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
      Istatus = HDF_SDS_READER(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
  
      do j = 1, ny
        do i = 1, nx
    
          if (space_mask(i,j) == sym%NO_SPACE) then
                
            ilat = max(1,min(grid_str%num_lat,int(abs(lat(i,j) - grid_str%first_lat)/grid_str%del_lat) + 1))
            ilon = max(1,min(grid_str%num_lon,int(abs(lon(i,j) - grid_str%first_lon)/grid_str%del_lon) + 1))
            ilat_ad = max(1,min((ilat - start_2d(2)) + 1,size(land_grid,2)))
            ilon_ad = max(1,min((ilon - start_2d(1)) + 1,size(land_grid,1)))
            land(i,j) = land_grid(ilon_ad,ilat_ad)
  
          end if
      
        end do
      end do
    
      deallocate(land_grid, stat=astatus)
      if (astatus /= 0) then
        call mesg("Error deallocating land surface grid. ",level=verb_lev % ERROR)
        stop
      end if
    
    else ! dateline_flg ne 0
  
      ilat1 = max(0,min(grid_str%num_lat,int(abs(nlat - grid_str%first_lat)/grid_str%del_lat) + 0))
      ilat2 = max(0,min(grid_str%num_lat,int(abs(slat - grid_str%first_lat)/grid_str%del_lat) + 0))
  
      ilon1 = max(0,min(grid_str%num_lon,int(abs(wlon - grid_str%first_lon)/grid_str%del_lon) + 0))
      ilon2 = max(0,min(grid_str%num_lon,int(abs(180.0 - grid_str%first_lon)/grid_str%del_lon) + 0))
    
      ilon1_2 = max(0,min(grid_str%num_lon,int(abs(-180.0 - grid_str%first_lon)/grid_str%del_lon) + 0))
      ilon2_2 = max(0,min(grid_str%num_lon,int(abs((elon-360.0) - grid_str%first_lon)/grid_str%del_lon) + 0))
  
      if (ilat1 > ilat2) then
        temp = ilat1
        ilat1 = ilat2
        ilat2 = temp
      end if
  
      if (ilon1 > ilon2) then
        temp = ilon1
        ilon1 = ilon2
        ilon2 = temp
      end if
    
      if (ilon1_2 > ilon2_2) then
        temp = ilon1_2
        ilon1_2 = ilon2_2
        ilon2_2 = temp
      end if
  
      start_2d = (/ilon1, ilat1/)
      stride_2d = (/1, 1/)
      edge_2d = (/(ilon2-ilon1)+1, (ilat2-ilat1)+1/)
  
      !call read_hdf_sds(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
      Istatus = HDF_SDS_READER(id, trim(grid_str%sds_name), start_2d, stride_2d, edge_2d, land_grid)
      start_2d_2 = (/ilon1_2, ilat1/)
      stride_2d_2 = (/1, 1/)
      edge_2d_2 = (/(ilon2_2-ilon1_2)+1, (ilat2-ilat1)+1/)
  
      !call read_hdf_sds(id, trim(grid_str%sds_name), start_2d_2, stride_2d_2, edge_2d_2, land_grid_2)
      Istatus = HDF_SDS_READER(id, trim(grid_str%sds_name), start_2d_2, stride_2d_2, edge_2d_2, land_grid_2)
    
      do j = 1, ny
        do i = 1, nx
    
          if (space_mask(i,j) == sym%NO_SPACE) then

            ilat = max(1,min(grid_str%num_lat,int(abs(lat(i,j) - grid_str%first_lat)/grid_str%del_lat) + 1))
            ilon = max(1,min(grid_str%num_lon,int(abs(lon(i,j) - grid_str%first_lon)/grid_str%del_lon) + 1))
            if (lon(i,j) >= 0.0) then
              ilat_ad = max(1,min((ilat - start_2d(2)) + 1,size(land_grid,2)))
              ilon_ad = max(1,min((ilon - start_2d(1)) + 1,size(land_grid,1)))
              land(i,j) = land_grid(ilon_ad,ilat_ad)
            else
              ilat_ad = max(1,min((ilat - start_2d_2(2)) + 1,size(land_grid_2,2)))
              ilon_ad = max(1,min((ilon - start_2d_2(1)) + 1,size(land_grid_2,1)))
              land(i,j) = land_grid_2(ilon_ad,ilat_ad)
            end if

          end if
      
        end do
      end do
  
      deallocate(land_grid, land_grid_2, stat=astatus)
      if (astatus /= 0) then
        call mesg("Error deallocating land surface grid. ",level=verb_lev % ERROR)
        stop
      end if
    
    end if

  end subroutine read_land_sfc_hdf_i2

end module LAND_SFC_PROPERTIES_MOD
