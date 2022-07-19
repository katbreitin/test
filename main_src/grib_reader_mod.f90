module grib_reader_module

! PURPOSE:
! A wrapper module for ecCodes API for reading data from a GRIB-file
! The data can be unit-converted following the conventions set in 'grib2hdf'
!
! REFS:
! ecCodes Home - ecCodes replaced GRIB API, which was developed at ECMWF to
!   decode and encode GRIB edition 1 and 2 data
!                ecCodes is a superset of GRIB API
! https://confluence.ecmwf.int/display/ECC
!
! ecCodes GRIB Training Course
! https://confluence.ecmwf.int/download/attachments/73011815/eccodes-keys-2017.pdf?api=v2
!
! GRIB Fortran 90 - Python APIs Parts 1 and 2
! file:///home/jmielikainen/Downloads/eccodes_grib_2019_f90_python_part1.pdf
! file:///home/jmielikainen/Downloads/eccodes_grib_2019_f90_python_part2.pdf
!
! [GRIB table for the NCEP operational files as of 10-21-96] -
!  TABLE 2. PARAMETERS & UNITS 1 & 2
! https://www.cpc.ncep.noaa.gov/products/wesley/opn_gribtable.html

  implicit none

  real, parameter, private :: missing_gfs = 9.999E+20
  real, parameter, private :: eps=1.0e-6
  integer, parameter, private :: n_grib_types=8
! TODO: use 1   !E2 NWP Nwp Model Option  (0=off,1=gfs,2=ncep reanalysis,3=cfsr)
!       from clavrx_options
  integer, dimension(n_grib_types), parameter, private ::  &
       n_prof = [26, 26, 26, 26, 26, 26, 26, 26]
  integer, dimension(n_grib_types), parameter, private ::  &
       n_prof_o3mr = [6, 6, 6, 6, 6, 6, 6, 6]
  integer, dimension(n_grib_types), parameter, private ::  &
       n_prof_rh = [21, 21, 21, 21, 21, 21, 21, 21]
  integer, dimension(n_grib_types), parameter, private ::  &
       n_prof_clwmr = [21, 21, 21, 21, 21, 21, 21, 21]
  integer, dimension(n_grib_types), parameter, private ::  &
       nx = [360, 720, 360, 360, 144, 720, 360, 1440]
  integer, dimension(n_grib_types), parameter, private ::  &
       ny = [181, 361, 181, 181, 72, 361, 181, 721]
  real, dimension(n_grib_types), parameter, private ::  &
       dlat = [1.0, 0.5, 1.0, 1.0, 2.5, 0.5, 1.0, 0.25]
  real, dimension(n_grib_types), parameter, private ::  &
       dlon = [1.0, 0.5, 1.0, 1.0, 2.5, 0.5, 1.0, 0.25]
  real, dimension(n_grib_types), parameter, private ::  &
       lat1 = [90.0, -90.0, 90.0, -90.0, -90.0, -90.0, -90.0, -90.0]
  real, dimension(n_grib_types), parameter, private ::  &
       lon1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

  real, parameter, private :: m_to_km = 1.e-3,  &
                              mm_to_cm = 0.1,  &
                              Pa_to_hPa = 1.e-2

  ! interface to grib_read for reading 1D, 2D or 3D data with 'grib2hdf' unit
  !  conversions
  interface grib_read
    module procedure grib_read_1d,    &
                     grib_read_2d,    &
                     grib_read_3d
  end interface

  ! interface to grib_readerlayer for reading a layer based on a pressure_level
  !  (int) or a layer_type (str) without unit conversions
  interface grib_read_layer
    module procedure grib_read_layer_int,    &
                     grib_read_layer_str
  end interface

  interface get_global_attribute
    module procedure get_global_attribute_int,    &
                     get_global_attribute_real
  end interface

  ! grib_read_3d() reads 26 predefined layers (21 for 'r' and 'clwmr' and 6 for
  !  'o3mr' with the other layers set to zero)
  integer, parameter :: n_std_prof = 26
  
  integer, dimension(n_std_prof), parameter ::  &
       P_std = [10 , 20 , 30 , 50 ,70 ,100, 150, 200, 250,  &
                300, 350, 400, 450, 500, 550, 600, 650, 700, 750,  &
                800, 850, 900, 925, 950, 975, 1000]

  type grib
    integer :: ifile  ! id of the opened file; to be used in all file functions
    integer :: idx_shortName  ! id of an index created from a file using
                              !  'shortName' as a key for the index
    integer :: idx_shortName_sVOFFS  ! id of an index created from a file using
                                     !  'shortName' and
                                     !  'scaledValueOfFirstFixedSurface' as keys
                                     !  for the index
    integer :: idx_shortName_tOL  ! id of an index created from a file using
                                  !  'shortName' and 'typeOfLevel' as keys
                                  !  for the index
  end type grib

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine grib_open
!
! PURPOSE:
! Open a GRIB file and create new indices from the opened file
! Three indices are used by different read types (ecCodes API requires that
!  each key has a value defined)
!
! INPUT:
!    grib_fname (str) - A pathname to a GRIB file
!
! OUTPUT:
!    grib_id (type(grib)) - A GRIB-reader handle

  subroutine grib_open(grib_fname, grib_id)

    use eccodes, only: codes_open_file, codes_index_create
      
    implicit none

    character (len=*), intent(IN) :: grib_fname  ! full path to a GRIB-file
    type(grib), intent(OUT) :: grib_id

    call codes_open_file(grib_id%ifile, grib_fname, 'r')

    call codes_index_create(grib_id%idx_shortName_sVOFFS, grib_fname,  &
                            'shortName,scaledValueOfFirstFixedSurface')
    call codes_index_create(grib_id%idx_shortName_tOL, grib_fname,  &
         'shortName,typeOfLevel')
    call codes_index_create(grib_id%idx_shortName, grib_fname, 'shortName')

  end subroutine grib_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine grib_close
!
! PURPOSE:
! Close a GRIB file and delete the indices
!
! INPUT:
!    grib_id (type(grib)) - A GRIB-reader handle

  subroutine grib_close(grib_id)

    use eccodes, only: codes_index_release, codes_close_file

    implicit none

    type(grib), intent(IN) :: grib_id

    call codes_index_release(grib_id%idx_shortName_sVOFFS)
    call codes_index_release(grib_id%idx_shortName_tOL)
    call codes_index_release(grib_id%idx_shortName)

    call codes_close_file(grib_id%ifile)

  end subroutine grib_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine get_global_attribute
!
! PURPOSE:
! Return a global attribute value corresponding to values grib2hdf writes to a
!  hdf-file
! 
! INPUT:
!    attr_name (str) - Attribute name
!    grib_type (int) - Corresponds to a grib_type value in grib2hdf
!
! OUTPUT: 
!    attr_value (int/real) - Attribute value

  subroutine get_global_attribute_int(grib_id, attr_name, attr_value,  &
       grib_type_in)
    implicit none

    type(grib), intent(IN) :: grib_id
    character(len=*), intent(IN) :: attr_name
    integer, intent(OUT) :: attr_value
    integer, intent(IN), optional :: grib_type_in

    integer :: grib_type

    if (present(grib_type_in)) then
      grib_type = grib_type_in
    else
      grib_type = 2  ! grib2hdf: grib_type == "2" (GFS grib2)
    end if

    if ((grib_type < 1) .or. (grib_type > n_grib_types)) then
      print *, 'get_global_attribute_int(): grib_type out of range'
      print *, 'value of grib_type is ', grib_type
    end if

    select case(to_lower(attr_name))
      case("number of pressure levels")
        attr_value = n_prof(grib_type)
      case("number of o3mr levels")
        attr_value = n_prof_o3mr(grib_type)
      case("number of rh levels")
        attr_value = n_prof_rh(grib_type)
      case("number of clwmr levels")
        attr_value = n_prof_clwmr(grib_type)
      case("number of latitudes")
        attr_value = ny(grib_type)
      case("number of longitudes")
        attr_value = nx(grib_type)
      case default
        print *, 'get_global_attribute_int(): unsupported attribute name'
        print *, 'attribute name is ', attr_name
        stop 1
    end select

  end subroutine get_global_attribute_int

  subroutine get_global_attribute_real(grib_id, attr_name, attr_value,  &
       grib_type_in)
    implicit none

    type(grib), intent(IN) :: grib_id
    character(len=*), intent(IN) :: attr_name
    real, intent(OUT) :: attr_value
    integer, intent(IN), optional :: grib_type_in

    integer :: grib_type

    if (present(grib_type_in)) then
      grib_type = grib_type_in
    else
      grib_type = 2  ! grib2hdf: grib_type == "2" (GFS grib2)
    end if

    if ((grib_type < 1) .or. (grib_type > n_grib_types)) then
      print *, 'get_global_attribute_real(): grib_type out of range'
      print *, 'value of grib_type is ', grib_type
    end if

    select case(to_lower(attr_name))
      case("latitude resolution")
        attr_value = dlat(grib_type)
      case("longitude resolution")
        attr_value = dlon(grib_type)
      case("first latitude")
        attr_value = lat1(grib_type)
      case("first longitude")
        attr_value = lon1(grib_type)
      case default
       print *, 'get_global_attribute_str(): unsupported attribute name'
       print *, 'attribute name is ', attr_name
       stop 1
    end select

  end subroutine get_global_attribute_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine grib_read_layer
!
! PURPOSE:
! Read a layer of real(4) data from a GRIB-file
! 
! REFERENCES:
! GRIB Parameters Database
! https://confluence.ecmwf.int/pages/viewpage.action?pageId=73016005
!
! INPUT:
!    grib_id (type(grib)) - A GRIB-reader handle
!    short_name (str) - A short name of the stored data in a GRIB-file, e.g.
!                       the pressure has the shortName 'p'
!
! Depending on the 4th parameter type, either grib_read_layer_int() or
!  grib_read_layer_str() is used when grib_read_layer() interface is called:
!    pressure_level (int) - Selects a layer from data using index key
!                           'scaledValueOfFirstFixedSurface', i.e.,
!                           scaledValueOfFirstFixedSurface == pressure_level
!    level_type (str) - Selects a layer from data using index key 'typeOfLevel',
!                        i.e., typeOfLevel == level_type 
!                       Common level types are 'surface', 'tropopause', and
!                        'sigma'
!
! OUTPUT: 
!    data (real, 2D array) - Unscaled 2-D array of data from a GRIB-file

  subroutine grib_read_layer_int(grib_id, short_name, data, pressure_level)

    use eccodes, only: codes_index_select, codes_new_from_index, codes_get,  &
         codes_get_size, codes_release

    implicit none

    type(grib), intent(IN) :: grib_id
    character(len=*), intent(IN) :: short_name
    real, dimension(:,:), intent(OUT) :: data
    integer, intent(IN) :: pressure_level

    integer :: idx
    integer :: igrib  ! id of the message loaded in memory
    integer :: numberOfColumns, numberOfRows
    integer :: numPoints
    real, dimension(:), allocatable :: values  ! temporary 1-d data array

    ! create an index from a grib file using some keys
    idx = grib_id%idx_shortName_sVOFFS

    ! Select the message subset with key==value
    call codes_index_select(idx, 'shortName', short_name)
    call codes_index_select(idx, 'scaledValueOfFirstFixedSurface',  &
         pressure_level)

    ! Create a new handle from an index after having selected the key values
    call codes_new_from_index(idx, igrib)

    call codes_get(igrib,'Ni', numberOfColumns)
    call codes_get(igrib,'Nj', numberOfRows)

    ! Get the size of the values array
    call codes_get_size(igrib,'values', numPoints)

    allocate(values(numPoints))

    ! Get data values
    call codes_get(igrib, 'values', values)
    call codes_release(igrib)

    data(:,:) = reshape(values, [numberOfColumns, numberOfRows])
    data(:,:) = data(:,size(data,2):1:-1)  ! reverse columns (why?)

    call codes_release(igrib)  ! Free memory from the message referred to by
                               !  igrib

    deallocate(values)
  end subroutine grib_read_layer_int

  subroutine grib_read_layer_str(grib_id, short_name, data, level_type)

    use eccodes, only: codes_index_select, codes_new_from_index, codes_get,  &
         codes_get_size, codes_release

    implicit none

    type(grib), intent(IN) :: grib_id
    character(len=*), intent(IN) :: short_name
    real, dimension(:,:), intent(OUT) :: data
    character(len=*), intent(IN), optional :: level_type

    integer :: idx
    integer :: igrib  ! id of the message loaded in memory
    integer :: numberOfColumns, numberOfRows
    integer :: numPoints
    real, dimension(:), allocatable :: values  ! temporary 1-d data array

    if (present(level_type)) then
      idx = grib_id%idx_shortName_tOL
      call codes_index_select(idx, 'shortName', short_name)
      call codes_index_select(idx, 'typeOfLevel', level_type)
    else
      idx = grib_id%idx_shortName
      call codes_index_select(idx, 'shortName', short_name)
    end if

    ! Create a new handle from an index after having selected the key values
    call codes_new_from_index(idx, igrib)

    call codes_get(igrib, 'Ni', numberOfColumns)
    call codes_get(igrib, 'Nj', numberOfRows)

    ! Get the size of the values array
    call codes_get_size(igrib, 'values', numPoints)

    allocate(values(numPoints))

    ! Get data values
    call codes_get(igrib, 'values', values)

    data(:,:) = reshape(values, [numberOfColumns, numberOfRows])
    data(:,:) = data(:,size(data,2):1:-1)  ! reverse columns (why?)

    call codes_release(igrib)  ! Free memory from the message referred to by
                               !  igrib
    
    deallocate(values)
  end subroutine grib_read_layer_str

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine grib_read
!
! PURPOSE:
! Read real(4) data from a GRIB-file and perform a unit conversion on it
! 
! REFERENCES: 
!  [1] ecCodes Home - ecCodes replaced GRIB API, which was developed at ECMWF
!                      to decode and encode GRIB edition 1 and 2 data
!                     ecCodes is a superset of GRIB API
!
!  [2] ecCodes GRIB Training Course
!      https://confluence.ecmwf.int/download/attachments/73011815/eccodes-keys-2017.pdf?api=v2
!
! INPUT:
!   grib_id (type(grib)) - A GRIB-reader handle
!   key_name (str) - Case-insensitive variable name, which can be either a name
!                     used by grib2hdf or an alias, e.g. 'u-wind', 'u',
!                     'u component of wind' and 'eastward_wind' can all be used
!                     to refer to the same data
!   array_order_in (str), optional - Array dimension order for a 3D-layer,
!                                     i.e., 'XYZ' or 'ZXY', case-insensitive so
!                                     'zyx' and 'zxy' works as well. The default
!                                     order is 'zxy'
!
! OUTPUT: 
!    data (real, array) - Unit-scaled data from a GRIB-file
!                         Depending on the dimensions of data, grib_read_1d(),
!                          grib_read_2d() or grib_read_3d() is used when the
!                          grib_read() interface is called

  subroutine grib_read_3d(grib_id, key_name, data, array_order_in)

    use eccodes, only: codes_index_select

    implicit none

    type(grib), intent(IN) :: grib_id
    character(len=*), intent(IN) :: key_name  ! Variable name parameter; can
                                              !  either be a short name ('sp')
                                              !  or the case-insensitive
                                              !  full name ('Surface pressure')
    real, dimension(:,:,:), intent(OUT) :: data ! unit scaled 3d-array
    character(len=*), intent(IN), optional :: array_order_in

    integer :: pressure_level
    integer :: idx  ! id of an index created from a file
    character(len=20) :: short_name
    integer :: i
    character(len=3) :: array_order
    integer :: n_profiles
    integer :: first_profile

    if (present(array_order_in)) then
      array_order = array_order_in
    else
      array_order = 'zxy'
    end if

    select case(to_lower(key_name))
      case ('t', 'temperature', 'air_temperature')
        short_name = 't'
      case ('height', 'gh', 'geopotential height', 'geopotential_height')
        short_name = 'gh'
      case ('u-wind', 'u', 'u component of wind', 'eastward_wind')
        short_name = 'u'
      case ('v-wind', 'v', 'v component of wind', 'northward_wind')
        short_name = 'v'
      case ('o3mr', 'ozone mixing ratio')
        short_name = 'o3mr'
      case ('rh', 'r', 'relative humidity')
        short_name = 'r'
      case ('clwmr', 'cloud mixing ratio')
        short_name = 'clwmr'
      case default
        print *, 'grib_read_3d() :: unknown 3d-field'
        stop 1
    end select

    ! create an index from a grib file using some keys
    idx = grib_id%idx_shortName_sVOFFS

    ! Select the message subset with key==value
    call codes_index_select(idx, 'shortName', short_name)

    select case (short_name)
      case ('r', 'clwmr')
        n_profiles = 21
        first_profile = 6
      case ('o3mr')
        n_profiles = 6
        first_profile = 1
      case default
        n_profiles = n_std_prof
        first_profile = 1
    end select

    data(:,:,:) = 0.
    do i=first_profile,first_profile+(n_profiles-1)

      pressure_level = P_std(i)*100

      if (to_lower(array_order) == 'zxy') then
        call grib_read_layer(grib_id, short_name, data(i,:,:), pressure_level)
      else
        call grib_read_layer(grib_id, short_name, data(:,:,i), pressure_level)
      end if
       
    end do

    ! do unit conversions here
    if (short_name == 'gh') then
      ! convert Z from [m] to [km]
      where (data < missing_gfs) data = data*m_to_km
    end if

  end subroutine grib_read_3d

  subroutine grib_read_2d(grib_id, key_name, data)

    use eccodes, only: codes_index_select

    implicit none

    type(grib), intent(IN) :: grib_id
    character(len=*), intent(IN) :: key_name  ! Variable name parameter; can
                                              !  either be a short name ('sp')
                                              !  or the case-insensitive
                                              !  full name ('Surface pressure')
    real, dimension(:,:), intent(OUT) :: data  ! unit scaled 2d-array

    integer :: idx  ! id of an index created from a file
    character(len=20) :: short_name
    logical :: sigma, surface, tropopause, user_level, level
    character(len=20) :: level_type

    sigma = .false.
    surface = .false.
    tropopause = .false.
    user_level = .false.
    level = .false.

    select case (to_lower(key_name))
      case ('sp', 'surface pressure', 'surface_air_pressure')
        short_name = 'sp'
      case ('msl pressure', 'prmsl', 'pressure reduced to msl')
        short_name = 'prmsl'
      case ('surface height', 'orog', 'orography', 'geopotential_height')
        short_name = 'orog'
        surface = .true.
      case ('land mask', 'lsm', 'land-sea mask', 'land_binary_mask')
        short_name = 'lsm'
      case ('ice fraction', 'ci', 'sea ice area fraction', 'siconc',  &
            'sea_ice_area_fraction')
        short_name = 'ci'
      case ('temperature at sigma=0.995', 't_sigma')
        short_name = 't'
        sigma = .true.
      case ('rh at sigma=0.995', 'r_sigma', 'relative humidity')
        short_name = 'r'
        sigma = .true.
      case ('surface temperature', 't_surf')
        short_name = 't'
        surface = .true.
      case ('total precipitable water', 'pwat', 'precipitable water')
        short_name = 'pwat'
      case ('u-wind at sigma=0.995', 'u_sigma', 'u component of wind',  &
            'eastward_wind')
        short_name = 'u'
        sigma = .true.
      case ('v-wind at sigma=0.995', 'v_sigma', 'v component of wind',  &
           'northward_wind')
        short_name = 'v'
        sigma = .true.
      case ('tropopause temperature', 't_trop')
        short_name = 't'
        tropopause = .true.
      case ('tropopause pressure', 'pres_trop', 'pressure')
        short_name = 'trpp'
        tropopause = .true.
      case ('water equivalent snow depth', 'sdwe',  &
            'water equivalent of accumulated snow depth')
        short_name = 'sdwe'
      case ('tozne', 'total ozone')
        short_name = 'tozne'
      case ('hpbl', 'planetary boundary layer height')
        short_name = 'hpbl'
      case ('cape', 'convective available potential energy')
        short_name = 'cape'
      case default
        print *, 'grib_read_2d() :: unknown 2d-field'
        stop 1
    end select

    ! create an index from a grib file using some keys
    if (tropopause .or. surface .or. sigma .or. user_level) then
      idx = grib_id%idx_shortName_tOL
    else
      idx = grib_id%idx_shortName
    end if

    ! Select the message subset with key==value
    call codes_index_select(idx, 'shortName', short_name)

    if (surface) then
      level_type = 'surface'
      level = .true.
    else if (tropopause) then
      level_type = 'tropopause'
      level = .true.
    else if (sigma) then
      level_type = 'sigma'
      level = .true.
    end if

    if (level) then
      call grib_read_layer(grib_id, short_name, data, level_type)
    else
      call grib_read_layer(grib_id, short_name, data)
    end if

    ! do unit conversions here
    select case (short_name)
      case ('sp')
        ! convert P_sfc from [Pa] to [hPa]
        where (data < missing_gfs) data = data*Pa_to_hPa
      case ('prmsl')
        ! convert Pmsl from [Pa] to [hPa]
        where (data < missing_gfs) data = data*Pa_to_hPa
      case ('pwat')
        ! FIX: convert pwat from [mm] to [cm] for GFS only (not CFSR)
        where (data < missing_gfs) data = data*mm_to_cm
      case ('trpp')
        ! convert P_trop from [Pa] to [hPa]
        where (data < missing_gfs) data = data*Pa_to_hPa
      case ('orog')
        ! convert Z_sfc from [m] to [km]
        where (data < missing_gfs) data = data*m_to_km
      case ('hpbl')
        ! convert Z_pbl from [m] to [km]
        where (data < missing_gfs) data = data*m_to_km
      case ('sdwe')
        ! why aren't the values equal to the grib2hdf data without this ?
        where (abs(data-9999.) < eps) data = 0.
    end select

  end subroutine grib_read_2d

  subroutine grib_read_1d(grib_id, key_name, data)
    implicit none

    type(grib), intent(IN) :: grib_id
    character(len=*), intent(IN) :: key_name  ! Variable name parameter,
                                              !  ONLY: "pressure levels"
    real, dimension(:), intent(OUT) :: data  ! array of standard pressure levels

    select case (key_name)
      case ("pressure levels")
        data = P_std
      case default
        print *, "grib_read_1d(): unsupported key_name", key_name
        stop 1
    end select
  end subroutine grib_read_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine grib_get_values
!
! PURPOSE:
! Get all values for a key in the index based on a short name
!
! INPUT:
!    grib_id (type(grib)) - A GRIB-reader handle
!    shortName (str) - Short name to be selected from the index, e.g. 'r' for
!                      relative humidity
!    key_name (str) - Key name, e.g. 'typeOfLevel' or
!                     'scaledValueOfFirstFixedSurface'
!
! REFERENCES:
! GRIB Parameters Database
! https://confluence.ecmwf.int/pages/viewpage.action?pageId=73016005
!
! OUTPUT:
!    names (str, array), optional - List of string values of a key from a
!                                    selected subset of messages
!    N (int), optional - Number of elements in the names array

  subroutine grib_get_values(grib_id, short_name, key_name, names, N)

    use eccodes, only: codes_index_select, codes_new_from_index, codes_get,  &
         codes_release, GRIB_END_OF_INDEX, grib_release
    
    implicit none

    type(grib), intent(IN) :: grib_id
    character(len=*), intent(IN) :: short_name
    character(len=*), intent(IN) :: key_name
    character(len=*), dimension(:), intent(OUT), optional :: names
    integer, intent(OUT), optional :: N

    integer :: idx  ! id of an index created from a file
    integer :: igrib  ! id of the message loaded in memory
    integer :: iret, count
    character(len=256) :: value

    idx = grib_id%idx_shortName

    ! Select the message subset with 'shortName' == short_name
    call codes_index_select(idx, 'shortName', short_name)

    ! Create a new handle from an index after having selected the key values
    call codes_new_from_index(idx, igrib, iret)

    count = 0
    do while (iret /= GRIB_END_OF_INDEX)  ! loop over all messages in the index
      count = count + 1
      call codes_get(igrib, key_name, value)  ! Select a key name value from a
                                              !  message
      call codes_release(igrib)
      if (present(names)) names(count) = value

      call codes_new_from_index(idx, igrib, iret)
    end do

    if (present(N)) N = count

    call codes_release(igrib)

    call grib_release(igrib)

  end subroutine grib_get_values

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine grib_info
!
! PURPOSE:
! Return a list of short names, names in ECMF-format and units in a GRIB-file
!
! INPUT:
!    type(grib) - A GRIB-reader handle
!
! REFERENCES:
! GRIB Parameters Database
! https://confluence.ecmwf.int/pages/viewpage.action?pageId=73016005
!
! OUTPUT:
!    shortName (str, array), optional - short names
!    nameECMF (str, array), optional - names in ECMF-format
!    units (str, array), optional - units
!    N (int), optional - number of elements in the array

  subroutine grib_info(grib_id, shortName, nameECMF, units, N)

    use eccodes, only: codes_index_get_size, codes_index_get, codes_get,  &
         codes_index_select, codes_new_from_index, codes_release

    implicit none

    type(grib), intent(IN) :: grib_id
    character(len=*), dimension(:), intent(OUT), optional :: shortName
    character(len=*), dimension(:), intent(OUT), optional :: nameECMF
    character(len=*), dimension(:), intent(OUT), optional :: units
    integer, intent(OUT), optional :: N

    integer :: idx  ! id of an index created from a file
    integer :: igrib  ! id of the message loaded in memory
    integer :: shortNameSize, i
    character(len=50), dimension(:), allocatable :: short_name

    ! create an index from a grib file using some keys
    idx = grib_id%idx_shortName

    ! get the number of distinct values of shortName in the index
    call codes_index_get_size(idx, 'shortName', shortNameSize)
    if (present(N)) N = shortNameSize

    ! allocate the array to contain the list of distinct shortName
    allocate(short_name(shortNameSize))

    ! get the list of distinct shortName from the index
    !   NOTE: valgrind reports a memory leak during the call below - ecCodes has
    !         implemented a fix for a similar problem on 2.13.0
    call codes_index_get(idx, 'shortName', short_name)
    if (present(shortName)) shortName = short_name

    if (present(nameECMF)) then
      do i=1,shortNameSize
        call codes_index_select(idx, 'shortName', short_name(i))
        call codes_new_from_index(idx, igrib)
        call codes_get(igrib, 'nameECMF', nameECMF(i))
        call codes_release(igrib)
      end do
    end if

    if (present(units)) then
      do i=1,shortNameSize
        call codes_index_select(idx, 'shortName', short_name(i))
        call codes_new_from_index(idx, igrib)
        call codes_get(igrib, 'units', units(i))
        call codes_release(igrib)
      end do
    end if

    deallocate(short_name)

  end subroutine grib_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function to_lower
!
! PURPOSE:
! Converts a string into lowercase
!
! INPUT:
!    strIN (str) - a string, of any case
!
! OUTPUT:
!    strOut (str) - the string, in lowercase

  function to_lower(strIn)  result(strOut)
    implicit none

    character(len=*), intent(IN) :: strIn
    character(len=len(strIn)) :: strOut

    integer :: i,j

    do i=1,len(strIn)
      j = iachar(strIn(i:i))
      if (j >= iachar("A") .and. j <= iachar("Z")) then
        strOut(i:i) = achar(iachar(strIn(i:i))+32)
      else
        strOut(i:i) = strIn(i:i)
      end if
    end do
  end function to_lower

end module grib_reader_module
