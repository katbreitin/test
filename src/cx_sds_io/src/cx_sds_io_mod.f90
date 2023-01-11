! $Id: cx_sds_io_mod.f90 4094 2021-03-04 13:32:04Z awalther $
!
!  ----  ------
!------------------------------------------------------------------------------
! SSEC/UW, CLAVR-x Software Tools Science Data Read Routines
!------------------------------------------------------------------------------
!
! MODULE: cx_sds_io_m0d
!
!> @author
!> Module Andi Walther CIMSS
!
! DESCRIPTION:
!> This module distributes the data read inquests to the science data libraries
!
! REVISION HISTORY:
! 12 Mar 2018 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
! 14 June 2018 : changes name of missing value to '_FillValue'
!------------------------------------------------------------------------------
module cx_sds_io_mod
   !#ifdef MFHDF
   use cx_hdf_read_mod, only: &
      hdf_get_finfo &
      , hdf_get_file_sds &
      , hdf_get_file_att
   !#endif
   use  cx_sds_type_definitions_mod, only: &
      cx_sds_type &
      , cx_att_type &
      , cx_sds_data_type

   use cx_ncdf_read_mod, only: &
      ncdf_get_finfo &
      , ncdf_get_file_sds &
      , ncdf_get_file_att &
      , ncdf_get_varinfo

   use cx_h5_read_mod, only: &
      h5_get_finfo &
      , h5_get_file_sds

   implicit none
   private

   include 'cx_sds_constants.inc'

   integer, parameter, public :: MAXNCDIM = 32
   integer, parameter, public :: MAXNCNAM = 128

   interface cx_sds_read
      module procedure cx_sds_read_2d_real &
         , cx_sds_read_1d_real &
         , cx_sds_read_3d_real &
         , cx_sds_read_4d_real &
         , cx_sds_read_5d_real &
         , cx_sds_read_6d_real &
         , cx_sds_read_7d_real

   end interface

   public :: cx_sds_finfo
   public :: cx_sds_varinfo
   public :: cx_sds_read_raw
   public :: cx_sds_read
   public :: cx_sds_att
   public :: cx_att_int
   public :: cx_att_r4

contains
   ! ------------------------------------------------------------------------------
   !  This tool determines file type
   ! ------------------------------------------------------------------------------

   function file_type (file)
      use hdf5
      integer :: file_type
      character(len = * ), intent(in) :: file
      character ( len = 3) :: postfix
      integer :: hdferr
      logical :: status

      postfix = file(len_trim(file) -2 : len_trim(file))

      if ( postfix .eq. 'hdf') file_type = 1

      if ( postfix .eq. '.nc') then
        file_type = 2
        ! some files with .ncshould be read with h5 tools
        call H5Fis_hdf5_f(file,status,hdferr)
        if (status) file_type = 3
      end if

      if ( postfix .eq. '.h5') file_type  = 3

   end function file_type

   ! ------------------------------------------------------------------------------
   !
   ! ------------------------------------------------------------------------------
   ! ---------------------------------------------------------------------------
   !> @author
   !> cx_sds_finfo Andi Walther CIMSS
   !
   ! DESCRIPTION:
   !> Returns information of sds datasets and attributes
   !> @brief
   !> Flow method (rate of change of position) used by integrator.
   !>
   !
   ! REVISION HISTORY:
   ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
   !
   !> @param[in] file
   !> @param[out] nsds Number of sds in this file
   !> @param[out] nsds Number of sds in this file
   !> @return status value
   !---------------------------------------------------------------------------


   function cx_sds_finfo ( file , ftype, nsds, sds_name, natt, att_name )
      implicit none

      integer :: cx_sds_finfo

      integer, intent(out) :: ftype
      character(len=*), intent(in) :: file  !> @var Variable description?
      integer, intent(out) :: nsds
      integer, intent(out) :: natt
      character ( len = MAXNCNAM), intent(out), allocatable :: sds_name(:)
      character ( len = MAXNCNAM), intent(out), allocatable :: att_name(:)


      ! - inialize
      cx_sds_finfo = 1

      !- first check which kind of file
      ! --
      ftype = 1 ! todo define params for HDF4, HDF5(2) and NCDF(3), not_existing(-1)
      !, not defined(0)
      ftype = file_type(file)


      !#ifdef MFHDF
      if ( ftype .eq. 1) cx_sds_finfo = hdf_get_finfo(file, nsds, sds_name, natt, att_name)
      !#endif
      if ( ftype .eq. 2) cx_sds_finfo = ncdf_get_finfo(file, nsds, sds_name, natt, att_name)

      if ( ftype .eq. 3) then
         cx_sds_finfo = h5_get_finfo(file, nsds, sds_name, natt, att_name)
      end if


   end function cx_sds_finfo


   !
   !
   !
   function cx_sds_varinfo ( file , var, ftype,  natt, att_name, ndim, dim )
      implicit none

      integer :: cx_sds_varinfo

      integer, intent(out) :: ftype
      character(len=*), intent(in) :: file
      character(len=*), intent(in) :: var
      integer, intent(out) :: natt
      character ( len = MAXNCNAM), intent(out), allocatable :: att_name(:)
      integer, intent(out) :: ndim
      integer, intent(out) :: dim(10)


      ! - inialize
      cx_sds_varinfo = 1

      !- first check which kind of file
      ! --
      ftype = 1 ! todo define params for HDF4, HDF5(2) and NCDF(3), not_existing(-1), not defined(0)
      ftype = file_type(file)

      !#ifdef MFHDF
      if ( ftype .eq. 1) print*, 'hallo' !cx_sds_finfo = hdf_get_finfo(file, nsds, sds_name, natt, att_name)
      !#endif
      if ( ftype .eq. 2) cx_sds_varinfo = ncdf_get_varinfo(file,  var, natt, att_name, ndim, dim)

      if ( ftype .eq. 3) then
         !cx_sds_finfo = h5_get_finfo(file, nsds, sds_name, natt, att_name)
      end if


   end function cx_sds_varinfo


   ! ------------------------------------------------------------------------------
   !
   ! ------------------------------------------------------------------------------
   function cx_sds_att (file, att_name, att )
      implicit none
      integer :: cx_sds_att
      character(len =*), intent(in) :: file
      character(len =*), intent(in) :: att_name
      type (cx_att_type),intent(out) :: att
      integer :: ftype

      integer :: natt
      type(cx_att_type),  allocatable :: attrs(:)

      integer :: i

      cx_sds_att = 0

      ftype = file_type(file)

      if ( ftype .eq. 1 ) then

         cx_sds_att = hdf_get_file_att(trim(file), natt, attrs)

      end if

      if ( ftype .eq. 2 ) then
         cx_sds_att = ncdf_get_file_att(trim(file), natt, attrs)
      end if

      do i = 1, natt

         if (trim(att_name) .EQ. trim(attrs(i) % name )) then
            att = attrs(i)
         end if

      end do


   end function cx_sds_att

   ! -----------------------------------------------------------------
   !
   ! This is a wrapper for integer output
   ! plaeas use this if you know that attribute is an integer
   function cx_att_int ( file, att_name, att_int)
      implicit none
      integer :: cx_att_int
      character(len =*), intent(in) :: file
      character(len =*), intent(in) :: att_name
      integer ,intent(out) :: att_int
      integer :: test
      type (cx_att_type) :: att

      test = cx_sds_att (file,att_name,att)
      cx_att_int = -1
      att_int = -999
      if (allocated (att % data % i4values)) then
         att_int = att % data % i4values(1)
         cx_att_int = 0
      end if

      if (allocated (att % data % i2values)) then
         att_int = att % data % i2values(1)
         cx_att_int = 0
      end if


   end function cx_att_int


   ! This is a wrapper for integer output
   ! plaeas use this if you know that attribute is an integer
   function cx_att_cha ( file, att_name, att_cha)
      implicit none
      integer :: cx_att_cha
      character(len =*), intent(in) :: file
      character(len =*), intent(in) :: att_name
      character ,intent(out) :: att_cha(:)
      integer :: test
      type (cx_att_type) :: att

      test = cx_sds_att (file,att_name,att)
      cx_att_cha = -1

      if (allocated (att % data % c1values)) then
         att_cha(1:att % data %dimsize(1)) = att % data % c1values(1:att % data %dimsize(1))
         cx_att_cha = 0
      end if




   end function cx_att_cha


    ! This is a wrapper for r4 output
   ! plaeas use this if you know that attribute is an integer
   function cx_att_r4 ( file, att_name, att_r4)
      implicit none
      integer :: cx_att_r4
      character(len =*), intent(in) :: file
      character(len =*), intent(in) :: att_name
      real ,intent(out) :: att_r4
      integer :: test
      type (cx_att_type) :: att

      test = cx_sds_att (file,att_name,att)
      cx_att_r4 = -1
      att_r4 = -999.
      if (allocated (att % data % r4values)) then
         att_r4 = att % data % r4values(1)
         cx_att_r4 = 0
      end if
   end function cx_att_r4


   ! ------------------------------------------------------------------------------
   !
   ! ------------------------------------------------------------------------------
   function cx_sds_read_raw ( file, sds_name, sds, start,count, stride )
      integer :: cx_sds_read_raw
      character (len = * ), intent(in) :: file
      character (len = * ), intent(in) :: sds_name
      type (cx_sds_type),  intent(out), allocatable, target :: sds(:)
      integer, optional, intent(in) :: start(:), stride(:), count(:)

      integer :: nsds, ii
      integer :: ftype

      cx_sds_read_raw = -1


      ftype = file_type(file)

      if (ftype .eq. 1 ) then
         !#ifdef MFHDF
         cx_sds_read_raw = hdf_get_file_sds(file, nsds,sds,1, (/sds_name/),start_inp=start &
          ,stride_inp = stride,count_inp =count, dclb = .false.,attr = .false.)

      !#endif

      end if

      ! ncdf
      if ( ftype .eq. 2 ) then
          ! do ii = 1,100000

          allocate(sds(1))

         cx_sds_read_raw = ncdf_get_file_sds(file, sds(1), sds_name &
         ,start_inp=start,stride_inp = stride,count_inp = count)
        ! call sds(1) % data % deallocate()
         ! call sds(1) % deallocate

          !if (allocated(sds)) deallocate(sds)
       !  print*,ii
       !  end do

          if ( cx_sds_read_raw .eq. -1 ) then
            cx_sds_read_raw = h5_get_file_sds(file, nsds,sds, (/sds_name/)  &
         ,start_inp=start,stride_inp = stride,count_inp =count)


            stop 'used h5 Routines for ncdf'
          end if
         !stop 'remove me later on cs_sds_io '
      end if

      if ( ftype .eq. 3 ) then

         cx_sds_read_raw = h5_get_file_sds(file, nsds,sds, (/sds_name/)  &
         ,start_inp=start,stride_inp = stride,count_inp =count)


      end if


      ! sds.copy_from_

   end function cx_sds_read_raw

   ! ------------------------------------------------------------------------------
   !
   ! ------------------------------------------------------------------------------
  function cx_sds_read_1d_real ( file, sds_name, out )
    integer :: cx_sds_read_1d_real
    character (len = * ), intent(in) :: file
    character (len = * ), intent(in) :: sds_name
    real, intent(out), allocatable :: out(:)
    real, allocatable:: temp_1d(:)
    type ( cx_sds_type), allocatable, target :: sds(:)
    type ( cx_sds_data_type), pointer :: pd => null()
    type ( cx_sds_type), pointer :: ps => null()
    integer :: dim1
    real :: add_offset(1)
    real :: slope (1)
    real :: missing(1)
    real :: scaled(1)
    logical :: att_exist

    ! ------------- executable -------

    if (  cx_sds_read_raw ( file, sds_name, sds) < 0 ) goto 9999

    pd=>sds(1)  % data
    ps=>sds(1)

    add_offset = ps % get_att('add_offset')
    slope = ps % get_att('scale_factor')
    scaled = ps % get_att('SCALED')
    missing = ps%get_att('missing',exist = att_exist)
    if ( .not. att_exist) missing = ps%get_att('_FillValue',exist = att_exist)

    allocate(temp_1d(pd % nval))
    call pd % transform_to_real(temp_1d)
      missing = ps%get_att('missing',exist = att_exist)
      if ( .not. att_exist) missing = ps%get_att('_FillValue',exist = att_exist)

    dim1 = pd%dimsize(1)
    if (dim1 .eq. 0) dim1 =1

    allocate(out(pd % nval))


    out = temp_1d

    if (scaled(1) .EQ. 1) then
      out = out * slope(1) + add_offset(1)
    end if

    where (temp_1d .EQ. missing(1))
      out = -999.
    end where



    cx_sds_read_1d_real = 0
9999 continue
    !cx_sds_read_1d_real = -1
    call pd % deallocate
    call ps % deallocate
    pd=>null()
    ps=>null()
    if (allocated(sds)) deallocate(sds)
    if ( allocated(temp_1d)) deallocate ( temp_1d)

  end function cx_sds_read_1d_real

  ! ------------------------------------------------------------------------------
  !
  ! ------------------------------------------------------------------------------
  integer function cx_sds_read_2d_real ( file, sds_name, out , start, stride, count)

    character (len = * ), intent(in) :: file
    character (len = * ), intent(in) :: sds_name
    integer, optional, intent(in) :: start(:), stride(:), count(:)
    real, intent(out), allocatable :: out(:,:)
    real, allocatable:: temp_1d(:)

    type ( cx_sds_type), allocatable, target :: sds(:)
    type ( cx_sds_data_type), pointer :: pd => null()
    type ( cx_sds_type), pointer :: ps => null()

    integer :: dim1, dim2, ii
    real :: add_offset(1)
    real :: slope (1)
    real :: missing(1)
    integer :: MISS_VALUE(1)
    real :: scaled(1)
    logical :: att_exist

    ! -  executable

    if (  cx_sds_read_raw ( file, sds_name, sds, start=start, stride=stride, count=count) < 0 ) goto 9999



    pd=>sds(1) % data
    ps=>sds(1)

    add_offset = ps %get_att('add_offset')
    slope = ps%get_att('scale_factor')
    scaled = ps%get_att('SCALED')
    MISS_VALUE = ps%get_att('missing',exist = att_exist)




    if ( .not. att_exist) MISS_VALUE = ps%get_att('_FillValue',exist = att_exist)


    if ( add_offset(1) .NE. -999.) scaled(1) = 1
    ! for  ATMS files which have bad attribute setting
    if (( add_offset(1) .EQ. -999.) .AND. (slope(1) .gt. 0)) then
      add_offset(:) = 0.
      scaled(1) = 1
    end if

    dim1 = pd%dimsize(1)
    dim2 = pd%dimsize(2)
    allocate(out(dim1,dim2))

    if ( allocated ( pd % r4values ) &
      .or. allocated ( pd % r8values ) &
      .or. allocated(pd % i4values) &
      .or. allocated(pd % i2values) &
      .or. allocated(pd % i1values) ) then

      allocate(temp_1d(pd%nval))

      call pd%transform_to_real(temp_1d)

      out = reshape (temp_1d,(/dim1,dim2/))

      if (scaled(1) .EQ. 1) then

        out = out * slope(1) + add_offset(1)
      end if

      where (reshape (temp_1d,(/dim1,dim2/)) .EQ. MISS_VALUE(1))
        out = -999.
      end where

    else if ( allocated ( pd % r4values_2d )) then
      out = pd % r4values_2d
    end if


    cx_sds_read_2d_real = 0
9999 continue
    ! cx_sds_read_2d_real = -1


    if ( allocated(temp_1d)) deallocate ( temp_1d)
    if ( associated(pd)) call pd % deallocate



    if ( associated(ps)) call ps % deallocate()
   !call ps % info
    if (allocated(sds)) deallocate(sds)

   pd=> null()
   ps => null()

  end function cx_sds_read_2d_real
  !
  !
  !
  function cx_sds_read_3d_real ( file, sds_name, out )
     integer :: cx_sds_read_3d_real
     character (len = * ), intent(in) :: file
     character (len = * ), intent(in) :: sds_name
     real, intent(out), allocatable :: out(:,:,:)
     real, allocatable:: temp_1d(:)
     type ( cx_sds_type), allocatable, target :: sds(:)
     type ( cx_sds_data_type), pointer :: pd => null()
     type ( cx_sds_type), pointer :: ps => null()
     integer :: dim1, dim2, dim3
     real :: add_offset(1)
     real :: slope (1)
     real :: missing(1)
     real :: scaled(1)


     if (  cx_sds_read_raw ( file, sds_name, sds) < 0 ) goto 9999
     pd=>sds(1)%data
     ps=>sds(1)

     add_offset = ps%get_att('add_offset')
     slope = ps%get_att('scale_factor')
     missing = ps%get_att('SCALED_MISSING')
     scaled = ps%get_att('SCALED')



     allocate(temp_1d(pd%nval))
     call pd%transform_to_real(temp_1d)




     dim1 = pd%dimsize(1)
     dim2 = pd%dimsize(2)
     dim3 = pd%dimsize(3)
     allocate(out(dim1,dim2,dim3))

     out = reshape (temp_1d,(/dim1,dim2,dim3/))
     if (scaled(1) .EQ. 1) then
        out = out * slope(1) + add_offset(1)
     end if


     where (reshape (temp_1d,(/dim1,dim2,dim3/)) .EQ. missing(1))
        out = -999.
     end where

     if ( allocated(temp_1d)) deallocate ( temp_1d)

     cx_sds_read_3d_real = 0
9999 continue
     !cx_sds_read_2d_real = -1
   if ( associated(pd)) call pd % deallocate
    pd=>null()
    ps=>null()
    if (allocated(sds)) deallocate(sds)

  end function cx_sds_read_3d_real
  !
  !
  !
  function cx_sds_read_4d_real ( file, sds_name, out )
     integer :: cx_sds_read_4d_real
     character (len = * ), intent(in) :: file
     character (len = * ), intent(in) :: sds_name
     real, intent(out), allocatable :: out(:,:,:,:)
     real, allocatable:: temp_1d(:)
     type ( cx_sds_type), allocatable, target :: sds(:)
     type ( cx_sds_data_type), pointer :: pd => null()
     type ( cx_sds_type), pointer :: ps => null()
     integer :: dim1, dim2, dim3,dim4
     real :: add_offset(1)
     real :: slope (1)
     real :: missing(1)
     real :: scaled(1)


     if (  cx_sds_read_raw ( file, sds_name, sds) < 0 ) goto 9999
     pd=>sds(1)%data
     ps=>sds(1)

     add_offset = ps%get_att('add_offset')
     slope = ps%get_att('scale_factor')
     missing = ps%get_att('SCALED_MISSING')
     scaled = ps%get_att('SCALED')



     allocate(temp_1d(pd%nval))
     call pd%transform_to_real(temp_1d)




     dim1 = pd%dimsize(1)
     dim2 = pd%dimsize(2)
     dim3 = pd%dimsize(3)
     dim4 = pd%dimsize(4)
     allocate(out(dim1,dim2,dim3,dim4))

     out = reshape (temp_1d,(/dim1,dim2,dim3,dim4/))
     if (scaled(1) .EQ. 1) then
        out = out * slope(1) + add_offset(1)
     end if


     where (reshape (temp_1d,(/dim1,dim2,dim3,dim4/)) .EQ. missing(1))
        out = -999.
     end where

     if ( allocated(temp_1d)) deallocate ( temp_1d)

     cx_sds_read_4d_real = 0
9999 continue
     !cx_sds_read_2d_real = -1
    if ( associated(pd)) call pd % deallocate
    pd=>null()
    ps=>null()
    if (allocated(sds)) deallocate(sds)

  end function cx_sds_read_4d_real
  !
  !



  ! -------------------------------------------------------
  !
  ! -------------------------------------------------------
  function cx_sds_read_5d_real ( file, sds_name, out )
     integer :: cx_sds_read_5d_real
     character (len = * ), intent(in) :: file
     character (len = * ), intent(in) :: sds_name
     real, intent(out), allocatable :: out(:,:,:,:,:)
     real, allocatable:: temp_1d(:)
     type ( cx_sds_type), allocatable, target :: sds(:)
     type ( cx_sds_data_type), pointer :: pd => null()
     type ( cx_sds_type), pointer :: ps => null()
     integer :: dim1, dim2, dim3, dim4 , dim5
     real :: add_offset(1)
     real :: slope (1)
     real :: missing(1)
     real :: scaled(1)


     if (  cx_sds_read_raw ( file, sds_name, sds) < 0 ) goto 9999
     pd=>sds(1)%data
     ps=>sds(1)

     add_offset = ps%get_att('add_offset')
     slope = ps%get_att('scale_factor')
     missing = ps%get_att('SCALED_MISSING')
     scaled = ps%get_att('SCALED')



     allocate(temp_1d(pd%nval))
     call pd % transform_to_real(temp_1d)




     dim1 = pd%dimsize(1)
     dim2 = pd%dimsize(2)
     dim3 = pd%dimsize(3)
     dim4 = pd%dimsize(4)
     dim5 = pd%dimsize(5)

     allocate(out(dim1,dim2,dim3,dim4,dim5))

     out = reshape (temp_1d,(/dim1,dim2,dim3,dim4,dim5/))

     if (scaled(1) .EQ. 1) then
        out = out * slope(1) + add_offset(1)
     end if


     where (reshape (temp_1d,(/dim1,dim2,dim3,dim4,dim5/)) .EQ. missing(1))
        out = -999.
     end where

     if ( allocated(temp_1d)) deallocate ( temp_1d)
     cx_sds_read_5d_real = 0
9999 continue
     !cx_sds_read_2d_real = -1
     call pd % deallocate
    pd=>null()
    ps=>null()
    if (allocated(sds)) deallocate(sds)

  end function cx_sds_read_5d_real
  !
  !-------------------------------------------------------

  function cx_sds_read_6d_real ( file, sds_name, out, start, stride, count )
     integer :: cx_sds_read_6d_real
     character (len = * ), intent(in) :: file
     character (len = * ), intent(in) :: sds_name
     real, intent(out), allocatable :: out(:,:,:,:,:,:)
     integer, optional, intent(in) :: start(:), stride(:), count(:)
     real, allocatable:: temp_1d(:)
     type ( cx_sds_type), allocatable, target :: sds(:)
     type ( cx_sds_data_type), pointer :: pd => null()
     type ( cx_sds_type), pointer :: ps
     integer :: dim1, dim2, dim3, dim4 , dim5,dim6
     real :: add_offset(1)
     real :: slope (1)
     real :: missing(1)
     real :: scaled(1)


     if (  cx_sds_read_raw ( file, sds_name, sds, start=start &
      , stride=stride, count=count) < 0 ) goto 9999
     pd=>sds(1)%data
     ps=>sds(1)

     add_offset = ps%get_att('add_offset')
     slope = ps%get_att('scale_factor')
     missing = ps%get_att('SCALED_MISSING')
     scaled = ps%get_att('SCALED')



     allocate(temp_1d(pd%nval))
     call pd % transform_to_real(temp_1d)


     dim1 = pd%dimsize(1)
     dim2 = pd%dimsize(2)
     dim3 = pd%dimsize(3)
     dim4 = pd%dimsize(4)
     dim5 = pd%dimsize(5)
     dim6 = pd%dimsize(6)

     allocate(out(dim1,dim2,dim3,dim4,dim5,dim6))

     out = reshape (temp_1d,(/dim1,dim2,dim3,dim4,dim5,dim6/))
     if (scaled(1) .EQ. 1) then
        out = out * slope(1) + add_offset(1)
     end if


     where (reshape (temp_1d,(/dim1,dim2,dim3,dim4,dim5,dim6/)) .EQ. missing(1))
        out = -999.
     end where

     if ( allocated(temp_1d)) deallocate ( temp_1d)

     cx_sds_read_6d_real = 0
9999 continue
     !cx_sds_read_2d_real = -1
    if (associated(pd)) call pd % deallocate
    pd=>null()
    ps=>null()
    if (allocated(sds)) deallocate(sds)

  end function cx_sds_read_6d_real





  ! -------------------------------------------------------
  !
  ! -------------------------------------------------------
  function cx_sds_read_7d_real ( file, sds_name, out )
     integer :: cx_sds_read_7d_real
     character (len = * ), intent(in) :: file
     character (len = * ), intent(in) :: sds_name
     real, intent(out), allocatable :: out(:,:,:,:,:,:,:)
     real, allocatable:: temp_1d(:)
     type ( cx_sds_type), allocatable, target :: sds(:)
     type ( cx_sds_data_type), pointer :: pd => null()
     type ( cx_sds_type), pointer :: ps => null()
     integer :: dim1, dim2, dim3, dim4 , dim5, dim6, dim7
     real :: add_offset(1)
     real :: slope (1)
     real :: missing(1)
     real :: scaled(1)


     if (  cx_sds_read_raw ( file, sds_name, sds) < 0 ) goto 9999
     pd=>sds(1)%data
     ps=>sds(1)

     add_offset = ps%get_att('add_offset')
     slope = ps%get_att('scale_factor')
     missing = ps%get_att('SCALED_MISSING')
     scaled = ps%get_att('SCALED')

     allocate(temp_1d(pd%nval))
     call pd%transform_to_real(temp_1d)

     dim1 = pd%dimsize(1)
     dim2 = pd%dimsize(2)
     dim3 = pd%dimsize(3)
     dim4 = pd%dimsize(4)
     dim5 = pd%dimsize(5)
     dim6 = pd%dimsize(6)
     dim7 = pd%dimsize(7)


     allocate(out(dim1,dim2,dim3,dim4,dim5,dim6,dim7))

     out = reshape (temp_1d,(/dim1,dim2,dim3,dim4,dim5,dim6,dim7/))
     if (scaled(1) .EQ. 1) then
        out = out * slope(1) + add_offset(1)
     end if


     where (reshape (temp_1d,(/dim1,dim2,dim3,dim4,dim5,dim6,dim7/)) .EQ. missing(1))
        out = -999.
     end where

     if ( allocated(temp_1d)) deallocate ( temp_1d)

     cx_sds_read_7d_real = 0
9999 continue
     !cx_sds_read_2d_real = -1
     call pd % deallocate
    pd=>null()
    ps=>null()
    if (allocated(sds)) deallocate(sds)

  end function cx_sds_read_7d_real
  end module
