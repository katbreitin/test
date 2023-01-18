module cx_geo__define

  use VIEWING_GEOMETRY_MOD
  use class_time_date, only: date_type
  use cx_real_boolean_mod
  implicit none
  type :: geo_str
    logical :: is_set
    real , dimension (:,:) , allocatable :: solzen
    real , dimension (:,:) , allocatable :: satzen
    real , dimension (:,:) , allocatable :: solaz
    real , dimension (:,:) , allocatable :: sataz
    real , dimension (:,:) , allocatable :: relaz
    real , dimension (:,:) , allocatable :: lat
    real , dimension (:,:) , allocatable :: lon
    real , dimension (:,:) , allocatable :: glintzen
    real , dimension (:,:) , allocatable :: scatangle
    real , dimension (:)   , allocatable :: scan_time
    logical, dimension (:,:), allocatable :: is_space
    type (date_type) :: time
    real :: sat_sub_lon
    real :: sat_sub_lat
    real :: sat_altitude_km


  contains
    procedure :: deallocate => deallocate_geo
    procedure :: allocate => allocate_geo
    procedure :: set => set_geo
  end type  geo_str

  ! --------------------------------------------------------------------------------------
  !
  ! --------------------------------------------------------------------------------------
contains
  subroutine set_geo(this,lon,lat, time, sat_sub_lon,sat_sub_lat,sat_alt)
    implicit none
    class ( geo_str ) :: this
    real,  intent(in) :: lon(:,:)
    real,  intent(in) :: lat(:,:)
    type (date_type), intent(in) :: time
    real, intent(in) :: sat_sub_lon,sat_sub_lat,sat_alt
    integer :: ndims(2)
    integer :: ii,jj
    integer :: day_of_year
    real :: hour_frac

    ndims = shape(lon)

    call this % allocate(ndims(1),ndims(2))
    call time  % get_date ( doy = day_of_year, hour_frac = hour_frac )

    this % lon = lon
    this % lat = lat

    do jj = 1 ,  ndims(2)
      do ii = 1 ,    ndims(1)
        if ( lon(ii,jj) .LER.  -199. ) cycle

        this % is_space(ii,jj) = .false.
        call  possol ( day_of_year ,  hour_frac  , lon(ii,jj) &
        ,  lat(ii,jj), this % solzen (ii,jj), this % solaz (ii,jj) )
      end do
    end do
    
    this % satzen  = sensor_zenith ( sat_alt , sat_sub_lon &
    ,sat_sub_lat ,lon,lat   )

    this % sataz  = sensor_azimuth ( sat_sub_lon,sat_sub_lat &
    ,lon , lat )

    this % relaz = relative_azimuth ( this % solaz   &
    , this % sataz )

    this % glintzen = glint_angle ( this % solzen  &
    , this % satzen  &
    ,  this % relaz )

    this % scatangle = scattering_angle (  this % solzen  &
    ,  this % satzen ,  this % relaz)


    where ( lon .lt. -200)
      this % satzen = -999.
      this % sataz  = -999.
      this % relaz  = -999.
      this % glintzen = -999.
      this % scatangle = -999.
      this % solzen = -999.
      this %solaz = -999.


    end where

    ! end do
    !   end do
    this % is_set = .true.


  end subroutine set_geo

  subroutine allocate_geo ( this, nx , ny )
    class ( geo_str ) :: this
    integer, intent(in) :: nx
    integer, intent(in) :: ny

    call this % deallocate

    allocate (  this  % lon        (nx , ny) )
    allocate (  this  % lat        (nx , ny) )
    allocate (  this  % solzen     (nx , ny) )
    allocate (  this  % solaz      (nx , ny) )
    allocate (  this  % satzen     (nx , ny) )
    allocate (  this  % sataz      (nx , ny) )
    allocate (  this  % relaz      (nx , ny) )
    allocate (  this  % glintzen   (nx , ny) )
    allocate (  this  % scatangle  (nx , ny) )
    allocate (  this  % is_space   (nx , ny) )

    this  % lon        =   -999.
    this  % lat        =   -999.
    this  % solzen     =   -999.
    this  % solaz      =   -999.
    this  % satzen     =   -999.
    this  % sataz      =   -999.
    this  % relaz      =   -999.
    this  % glintzen   =   -999.
    this  % scatangle  =   -999.
    this  % is_space = .true.


  end subroutine allocate_geo

  ! --------------------------------------------------------------------------------------
  !
  ! --------------------------------------------------------------------------------------
  subroutine deallocate_geo (this )
    class ( geo_str ) :: this
    if (allocated ( this % lon)) deallocate ( this % lon)
    if (allocated ( this % lat)) deallocate ( this % lat)


    if ( allocated  (  this % solzen   ) ) deallocate (  this  % solzen   )
    if ( allocated  (  this % solaz      ) ) deallocate (  this % solaz   )
    if ( allocated  (  this % satzen     ) ) deallocate (  this  % satzen  )
    if ( allocated  (  this % sataz     ) ) deallocate (  this  % sataz  )
    if ( allocated  (  this % relaz      ) ) deallocate (  this  % relaz   )
    if ( allocated  (  this % glintzen    ) ) deallocate (  this  % glintzen   )
    if ( allocated  (  this % scatangle   ) ) deallocate (  this  % scatangle   )
    if ( allocated  (  this % is_space    ) ) deallocate (  this  % is_space  )
    if ( allocated  (  this % scan_time    ) ) deallocate (  this  % scan_time )
    this % is_set = .false.

  end subroutine deallocate_geo




end module   cx_geo__define
