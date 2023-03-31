! $Id: cx_sds_type_definitions_mod.f90 3164 2019-02-21 17:02:23Z awalther $
!
Module cx_sds_type_definitions_mod

  integer, parameter, public :: MAXNCDIM = 32
  integer, parameter, public :: MAXNCNAM = 128

  include 'cx_sds_constants.inc'

  type cx_sds_data_type
    integer :: type
    integer :: type_ncdf
    integer :: type_hdf
    integer :: datasize
    integer :: size
    integer :: nval
    integer :: rank
    integer :: dimsize (MAXNCDIM)
    integer :: utype
    integer :: calbrtd

    real ( kind = 8 ) :: calibr(4)
    character, allocatable :: c1values(:)
    character(len=100) :: char_nc_values
    integer(kind=1), allocatable :: i1values(:)
    integer(kind=2), allocatable :: i2values(:)
    integer(kind=4), allocatable :: i4values(:)
    integer(kind=8), allocatable :: i8values(:)
    real(kind=4), allocatable :: r4values(:)
    real(kind=8), allocatable :: r8values(:)
    integer(kind=1), allocatable :: i1values_2d(:,:)
    integer(kind=2), allocatable :: i2values_2d(:,:)
    integer(kind=4), allocatable :: i4values_2d(:,:)
    integer(kind=8), allocatable :: i8values_2d(:,:)
    real(kind=4), allocatable :: r4values_2d(:,:)
    real(kind=8), allocatable :: r8values_2d(:,:)
    integer(kind=1), allocatable :: i1values_3d(:,:,:)
    integer(kind=2), allocatable :: i2values_3d(:,:,:)
    integer(kind=4), allocatable :: i4values_3d(:,:,:)
    real(kind=4), allocatable :: r4values_3d(:,:,:)
    real(kind=8), allocatable :: r8values_3d(:,:,:)
    integer(kind=1), allocatable :: i1values_4d(:,:,:,:)
    integer(kind=2), allocatable :: i2values_4d(:,:,:,:)
    integer(kind=4), allocatable :: i4values_4d(:,:,:,:)
    real(kind=4), allocatable :: r4values_4d(:,:,:,:)
    real(kind=8), allocatable :: r8values_4d(:,:,:,:)
    integer(kind=1), allocatable :: i1values_5d(:,:,:,:,:)
    integer(kind=2), allocatable :: i2values_5d(:,:,:,:,:)
    integer(kind=4), allocatable :: i4values_5d(:,:,:,:,:)
    real(kind=4), allocatable :: r4values_5d(:,:,:,:,:)
    real(kind=8), allocatable :: r8values_5d(:,:,:,:,:)

  contains
    procedure :: info=>cx_sds_data_type__info
    procedure :: transform_to_real
    procedure :: transform_to_real_2d
    procedure :: transform_to_dbl
    procedure :: type_ncdf_to_hdf
    procedure :: deallocate => cx_sds_data_type__deallocate
  end type cx_sds_data_type

  type cx_att_type
    character ( len = MAXNCNAM) :: name
    type ( cx_sds_data_type ) :: data



  end type cx_att_type

  type cx_sds_type
    character( len = MAXNCNAM) :: name
    integer :: nattr
    type ( cx_sds_data_type) :: data
    type(cx_att_type), dimension(:),allocatable :: attr
  contains
    procedure :: info=>cx_sds_type__info
    procedure :: get_att => cx_sds_type__get_att
    procedure :: deallocate =>cx_sds_type__deallocate
  end type cx_sds_type


contains
  subroutine cx_sds_type__deallocate (self)
    class ( cx_sds_type) :: self

    if ( allocated(self % attr)) deallocate (self % attr)

  end subroutine cx_sds_type__deallocate

  subroutine type_ncdf_to_hdf (self)
    class(cx_sds_data_type) :: self

    ! fake netcdf parameters

    integer,parameter ::  NF90_CHAR  =2
    integer,parameter ::   NF90_SHORT  =3
    integer,parameter ::  NF90_INT   = 4
    integer,parameter ::  NF90_FLOAT  = 5
    integer,parameter ::  NF90_DOUBLE = 6
    integer,parameter ::  NF90_USHORT = 8



    select case (self % type_ncdf)
    case (NF90_CHAR  ); self %  type = DFNT_CHAR8
    case (NF90_SHORT   ); self % type = DFNT_INT16
    case (NF90_INT  );self % type = DFNT_INT32
    case (NF90_FLOAT); self % type = DFNT_FLOAT32
    case (NF90_DOUBLE); self % type = DFNT_FLOAT64
    case (NF90_USHORT);  self % type = DFNT_UINT16
    case default       ; self % type = DFNT_FLOAT32
    end select




  end subroutine type_ncdf_to_hdf
  !  -------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------
  subroutine cx_sds_data_type__info (self)
    class(cx_sds_data_type):: self

    print*,'data size: ', self%size

    print*,'data rank: ',self%rank

    !  print*,'dimsize: ',self%dimsize(1:self%rank)
    !print*,'Calibr: ',self%calibr
    print*,'data_type: ',self%type

    select case ( self%type)
    case (DFNT_CHAR8)
      print*,'CHAR8 '
    case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
      print*,'UINT8 '
      print*,'shape data ..: ',shape(self % i1values)
    case (DFNT_UINT16, DFNT_INT16)
      print*,'UINT16 '
    case (DFNT_UINT32, DFNT_INT32)
      print*,'UINT32 '
    case (DFNT_FLOAT32)
      print*,'FLOAT32 '
    case (DFNT_FLOAT64)
      print*,'FLOAT64 '
    end select

    print*, 'Allocated c1? ',(allocated( self% c1values))
    print*, 'Allocated i1? ',(allocated( self% i1values))
    print*, 'Allocated i2? ',(allocated( self% i2values))
    print*, 'Allocated i4? ',(allocated( self%i4values))
    print*, 'Allocated r4? ',(allocated( self%r4values))
    print*, 'Allocated r8? ',(allocated( self%r8values))


  end subroutine cx_sds_data_type__info


  !-------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------
  subroutine transform_to_real (self, data_real)
    class(cx_sds_data_type), intent(in) :: self
    real, intent(out) :: data_real(:)


    select case ( self%type)
    case (DFNT_CHAR8)

    case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
      data_real = real (self % i1values)
      !    case ( DFNT_INT16)
      !      data_real = real (self % i2values)
      !    case (DFNT_UINT32, DFNT_INT32,DFNT_UINT16)
      !       data_real = real (self % i4values)
    case ( DFNT_INT16,DFNT_UINT16)
      if (allocated(self % i2values)) data_real = real (self % i2values)
      if (allocated(self % i4values)) data_real = real (self % i4values)
    case (DFNT_UINT32, DFNT_INT32)
      if (allocated(self % i4values)) data_real = real (self % i4values)
      if (allocated(self % r4values)) data_real = real (self % r4values)

    case (DFNT_FLOAT32 )

      if (allocated(self % r4values)) data_real = real (self % r4values)

    case (DFNT_FLOAT64)
      data_real = real (self % r8values)
    case default

      print*,'Missing data type transform_to_real in cx_sds_type_definitions_mod.f90'
    end select





  end subroutine transform_to_real



  !-------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------
  subroutine transform_to_real_2d (self, data_real)
    class(cx_sds_data_type), intent(in) :: self
    real, intent(inout) :: data_real(:,:)

    select case ( self%type)
    case (DFNT_CHAR8)

    case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
      data_real = real (self % i1values_2d)
    case (DFNT_UINT16, DFNT_INT16)
      data_real = real (self % i2values_2d)
    case (DFNT_UINT32, DFNT_INT32)
      data_real = real (self % i4values_2d)
    case (DFNT_FLOAT32)
      data_real = real (self % r4values_2d)
    case (DFNT_FLOAT64)
      data_real = real (self % r8values_2d)
    case default
      print*,'something wrong2'
    end select

  end subroutine transform_to_real_2d


  !-------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------
  subroutine transform_to_dbl (self, data_dbl)
    class(cx_sds_data_type), intent(in) :: self
    real(kind = 8), intent(out) :: data_dbl(:)

    select case ( self%type)
    case (DFNT_CHAR8)

    case (DFNT_UCHAR8, DFNT_UINT8, DFNT_INT8)
      data_dbl = real (self % i1values)
    case (DFNT_UINT16, DFNT_INT16)
      data_dbl = real (self % i2values)
    case (DFNT_UINT32, DFNT_INT32)
      data_dbl = real (self % i4values)
    case (DFNT_FLOAT32)
      data_dbl = real (self % r4values)
    case (DFNT_FLOAT64)
      data_dbl = real (self % r8values)
    case default
      print*,'something wrong1'
    end select

  end subroutine transform_to_dbl

  !-------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------
  subroutine cx_sds_type__info(self)
    class(cx_sds_type):: self
    integer :: i
    print*,self%name
    print*,'info SDS TYPE'
    print*,'number attributes: ', self%nattr

    do i=1, self%nattr
      print*,self%attr(i) % name

    end do

  end subroutine cx_sds_type__info
  !-------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------
  function cx_sds_type__get_att (self, att_name,exist)
    class(cx_sds_type) , target :: self
    character(len=*) :: att_name
    real :: cx_sds_type__get_att(1)
    type ( cx_sds_data_type), pointer :: pd
    logical, optional, intent(out) :: exist
    integer :: i
    cx_sds_type__get_att = -999.
    if ( present(exist)) exist = .false.

    do i =1, self % nattr

      if (self % attr(i) % name .EQ. trim(att_name)) then
        pd=>self % attr(i) % data

        call pd  % transform_to_real(cx_sds_type__get_att)

        if ( present(exist)) exist=.true.
      end if

    end do

  end function cx_sds_type__get_att





  !
  !
  !
  subroutine cx_sds_data_type__deallocate(self)
    class(cx_sds_data_type)  :: self

    if ( allocated( self%i1values)) deallocate( self%i1values)
    if ( allocated( self%i2values)) deallocate( self%i2values)
    if ( allocated( self%i4values)) deallocate( self%i4values)
    if ( allocated( self%i8values)) deallocate( self%i8values)
    if ( allocated( self%r4values)) deallocate( self%r4values)
    if ( allocated(self% r8values)) deallocate( self%r8values)
    if ( allocated( self%i1values_2d)) deallocate( self%i1values_2d)
    if ( allocated( self%i2values_2d)) deallocate( self%i2values_2d)
    if ( allocated( self%i4values_2d)) deallocate( self%i4values_2d)
    if ( allocated( self%i8values_2d)) deallocate( self%i8values_2d)
    if ( allocated( self%r4values_2d)) deallocate( self%r4values_2d)
    if ( allocated( self%r8values_2d)) deallocate( self%r8values_2d)
    if ( allocated( self%i1values_3d)) deallocate( self%i1values_3d)
    if ( allocated( self%i2values_3d)) deallocate( self%i2values_3d)
    if ( allocated( self%i4values_3d)) deallocate( self%i4values_3d)
    if ( allocated( self%r4values_3d)) deallocate( self%r4values_3d)
    if ( allocated( self%r8values_3d)) deallocate( self%r8values_3d)
    if ( allocated( self%i1values_4d)) deallocate( self%i1values_4d)
    if ( allocated( self%i2values_4d)) deallocate( self%i2values_4d)
    if ( allocated( self%i4values_4d)) deallocate( self%i4values_4d)
    if ( allocated( self%r4values_4d)) deallocate( self%r4values_4d)
    if ( allocated( self%r8values_4d)) deallocate( self%r8values_4d)
    if ( allocated( self%i1values_5d)) deallocate( self%i1values_5d)
    if ( allocated( self%i2values_5d)) deallocate( self%i2values_5d)
    if ( allocated( self%i4values_5d)) deallocate( self%i4values_5d)
    if ( allocated( self%r4values_5d)) deallocate( self%r4values_5d)
    if ( allocated( self%r8values_5d)) deallocate( self%r8values_5d)

  end subroutine cx_sds_data_type__deallocate


end module cx_sds_type_definitions_mod
