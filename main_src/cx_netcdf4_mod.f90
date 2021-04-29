!$Id: cx_netcdf4_mod.f90 3961 2020-09-10 18:08:37Z heidinger $
!----------------------------------------------------------------------
! MODULE name: CX_NETCDF4_MOD
! 
! Routines for reading in netCDF files, such as AHI data
!
! Authors: Andrew Heidinger, NOAA/NESDIS
!          Andi Walther, CIMSS
!          Denis Botambekov, CIMSS
!
! DEPENDENCIES: Constants, NETCDF library
!
!
! New netCDF read interface call is as follows:
!
! netcdf_read_array (nc_id, start, stride, end, variable name, output array)
!
!----------------------------------------------------------------------
module CX_NETCDF4_MOD

 use NETCDF,only: &
     NF90_OPEN &
   , NF90_NOWRITE &
   , NF90_WRITE &
   , NF90_NOERR &
   , NF90_PUT_ATT &
   , NF90_GET_ATT &
   , NF90_GLOBAL &
   , NF90_CLOSE &
   , NF90_INQUIRE_DIMENSION &
   , NF90_INQ_DIMID &
   , NF90_INQ_VARID &
   , NF90_STRERROR &
   , NF90_GET_VAR &
   , NF90_INQUIRE_VARIABLE   &
   , NF90_INQ_NCID

 use LEVEL2_STRUCTURES_MOD, only: Sds_Struct, L2_Glob_Attr_Definition, Clavrx_Global_Attr
 
 use CONSTANTS_MOD

 implicit none

 private
 
 public :: read_netcdf
 public :: read_and_unscale_netcdf_1d
 public :: read_and_unscale_netcdf_2d
 public :: read_netcdf_variable_attribute_real
 public :: read_netcdf_attribute_real
 public :: read_netcdf_attribute_double
 public :: read_netcdf_global_attribute
 public :: read_ahi_nav_coeff
 public :: read_ahi_time
 public :: read_abi_time
 public :: read_netcdf_dimension
 public :: read_netcdf_dimension_1d
 public :: read_netcdf_dimension_2d
 public :: read_netcdf_dimension_3d
 public :: read_netcdf_dimension_4d
 public :: open_netcdf
 public :: open_netcdf_write
 public :: close_netcdf
 public :: get_group_id
 public :: read_abi_kappa0
 public :: read_abi_max_focal_plane_temp
 public :: WRITE_CLAVRX_NETCDF_GLOBAL_ATTRIBUTES

!interface read_and_unscale_netcdf
!    module procedure &
!         read_and_unscale_netcdf_2d, &
!         read_and_unscale_netcdf_1d
!end interface read_and_unscale_netcdf

 interface read_netcdf
     module procedure &
          read_netcdf_1d_double, &
          read_netcdf_1d_real, &
          read_netcdf_1d_int,  &
          read_netcdf_2d_real, &
          read_netcdf_2d_int1,  &
          read_netcdf_2d_int2,  &
          read_netcdf_2d_int4,  &
          read_netcdf_2d_char, &
          read_netcdf_3d_real, &
          read_netcdf_4d_real, &
          read_netcdf_4d_int4
 end interface read_netcdf

 interface read_netcdf_global_attribute
     module procedure &
          read_netcdf_global_attribute_i4, &
          read_netcdf_global_attribute_r4, &
          read_netcdf_global_attribute_char
 end interface read_netcdf_global_attribute

 integer, parameter, private :: sds_rank_1d = 1
 integer, dimension(sds_rank_1d), private :: sds_start_1d, sds_edge_1d, sds_stride_1d
   
 integer, parameter, private :: sds_rank_2d = 2
 integer, dimension(sds_rank_2d), private :: sds_start_2d, sds_edge_2d, sds_stride_2d, sds_dims_2d
   
 integer, parameter, private :: sds_rank_3d = 3
 integer, dimension(sds_rank_3d), private :: sds_start_3d, sds_edge_3d, sds_stride_3d, sds_dims_3d
   
 integer, parameter, private :: sds_rank_4d = 4
 integer, dimension(sds_rank_4d), private :: sds_start_4d, sds_edge_4d, sds_stride_4d, sds_dims_4d

 character(len=30), parameter, private:: EXE_PROMPT_NAV = "CLAVR-x NetCdf Read Module >> "

 real, parameter, private:: Missing_Value_Netcdf = -999.0

 contains
 !------------------------------------------------------------------------------
 ! open netcdf file for reading
 !------------------------------------------------------------------------------
 subroutine open_netcdf(nc_file, ncid)
    character(len=*), intent(in):: nc_file
    integer, intent(out):: ncid
    integer:: status

    status = nf90_open(trim(nc_file), mode = nf90_nowrite, ncid = ncid)
    if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed Open on NetCdf File'
       print*, EXE_PROMPT_NAV ,' filename is: ', nc_file
       return
    endif
 end subroutine open_netcdf

 !------------------------------------------------------------------------------
 ! open netcdf file for writing
 !------------------------------------------------------------------------------
 subroutine open_netcdf_write(nc_file, ncid)
    character(len=*), intent(in):: nc_file
    integer, intent(out):: ncid
    integer:: status
    status = nf90_open(trim(nc_file), mode = nf90_write, ncid = ncid)
    if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed Open on NetCdf File'
       print*, EXE_PROMPT_NAV ,' filename is: ', nc_file
       return
    endif
 end subroutine open_netcdf_write

 !------------------------------------------------------------------------------
 !
 !------------------------------------------------------------------------------
 subroutine close_netcdf(ncid)
    integer,intent(in):: ncid
    integer:: status

    status = nf90_close(ncid)
    if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed Close on NetCdf File'
       return
    endif
 end subroutine close_netcdf

 !------------------------------------------------------------------------------
 !
 !------------------------------------------------------------------------------
 subroutine get_group_id(ncid, group_name, group_id)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: group_name
    integer, intent(out) :: group_id
    integer:: status

    Status = nf90_inq_ncid(ncid, trim(group_name), group_id)
   !   group_id = 0
    if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed to Get GROUP_ID = ',trim(group_name)
       return
    endif
 end subroutine get_group_id

 !------------------------------------------------------------------------------
 !
 !------------------------------------------------------------------------------
 subroutine read_netcdf_dimension_1d(ncid, dim_name, dim_value)
  integer, intent(in) :: ncid
  character(len=*), intent(in):: dim_name
  integer(kind=4), dimension(1), intent(out):: dim_value
  integer, dimension(1):: dimid
  integer:: status, sds_id, dim_id
  character(len=10):: dim_name_dummy

  !status = nf90_inq_dimid(ncid, trim(dim_name), dim_id)
  !status = nf90_inquire_dimension(ncid,dim_id,dim_name_dummy,len=dim_value(1))


  status = nf90_inq_varid(ncid, trim(dim_name), sds_id)
  status = nf90_inquire_variable(ncid, sds_id, dimids = dimid)
  status = nf90_inquire_dimension(ncid, dimid(1), len=dim_value(1))

  if (status /= nf90_noerr) then
     print *, EXE_PROMPT_NAV , 'ERROR: Failed to Get Dimension = ', trim(dim_name)
     return
  endif
 end subroutine read_netcdf_dimension_1d

 !------------------------------------------------------------------------------
 !
 !------------------------------------------------------------------------------
 subroutine read_netcdf_dimension_2d(ncid, dim_name, dim_value)
  integer, intent(in) :: ncid
  character(len=*), intent(in):: dim_name
  integer(kind=4), dimension(2), intent(out):: dim_value
  integer, dimension(2):: dimid
  integer:: status, sds_id

  status = nf90_inq_varid(ncid, trim(dim_name), sds_id)
  status = nf90_inquire_variable(ncid, sds_id, dimids = dimid)
  status = nf90_inquire_dimension(ncid, dimid(1), len=dim_value(1))
  status = nf90_inquire_dimension(ncid, dimid(2), len=dim_value(2))

  if (status /= nf90_noerr) then
     print *, EXE_PROMPT_NAV , 'ERROR: Failed to Get Dimension = ', trim(dim_name)
     return
  endif

 end subroutine read_netcdf_dimension_2d

 !------------------------------------------------------------------------------
 !
 !------------------------------------------------------------------------------
 subroutine read_netcdf_dimension_3d(ncid, dim_name, dim_value)
  integer, intent(in) :: ncid
  character(len=*), intent(in):: dim_name
  integer(kind=4), dimension(3), intent(out):: dim_value
  integer, dimension(3):: dimid
  integer:: status, sds_id

  status = nf90_inq_varid(ncid, trim(dim_name), sds_id)
  status = nf90_inquire_variable(ncid, sds_id, dimids = dimid)
  status = nf90_inquire_dimension(ncid, dimid(1), len=dim_value(1))
  status = nf90_inquire_dimension(ncid, dimid(2), len=dim_value(2))
  status = nf90_inquire_dimension(ncid, dimid(3), len=dim_value(3))

    if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed to Get Dimension = ', trim(dim_name)
       return
    endif
 end subroutine read_netcdf_dimension_3d

 !------------------------------------------------------------------------------
 !
 !------------------------------------------------------------------------------
 subroutine read_netcdf_dimension_4d(ncid, dim_name, dim_value)
  integer, intent(in) :: ncid
  character(len=*), intent(in):: dim_name
  integer(kind=4), dimension(4), intent(out):: dim_value
  integer, dimension(4):: dimid
  integer:: status, sds_id

  status = nf90_inq_varid(ncid, trim(dim_name), sds_id)
  status = nf90_inquire_variable(ncid, sds_id, dimids = dimid)
  status = nf90_inquire_dimension(ncid, dimid(1), len=dim_value(1))
  status = nf90_inquire_dimension(ncid, dimid(2), len=dim_value(2))
  status = nf90_inquire_dimension(ncid, dimid(3), len=dim_value(3))
  status = nf90_inquire_dimension(ncid, dimid(4), len=dim_value(4))

    if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed to Get Dimension = ', trim(dim_name)
       return
    endif
 end subroutine read_netcdf_dimension_4d

 !------------------------------------------------------------------------------
 !
 !------------------------------------------------------------------------------
 subroutine read_netcdf_global_attribute_i4(nc_file, attr_name, attr_value)
  character(len=*), intent(in):: nc_file, attr_name
  integer(kind=4), intent(out):: attr_value
  integer:: ncid, status

   status = nf90_open(nc_file, mode = nf90_nowrite, ncid = ncid)

   if (status /= nf90_noerr) then
      print *, EXE_PROMPT_NAV , 'ERROR: Failed Open to Read NETCDF Global Attribute'
      print*, EXE_PROMPT_NAV ,' filename is: ', nc_file
      return
   endif

   status = nf90_get_att(ncid, nf90_global, attr_name, attr_value)
   status = nf90_close(ncid)

 end subroutine read_netcdf_global_attribute_i4

 subroutine read_netcdf_global_attribute_r4(nc_file, attr_name, attr_value)
  character(len=*), intent(in):: nc_file, attr_name
  real(kind=4), intent(out):: attr_value
  integer:: ncid, status

  attr_value = Missing_Value_Netcdf

  status = nf90_open(nc_file, mode = nf90_nowrite, ncid = ncid)

  if (status /= nf90_noerr) then
      print *, EXE_PROMPT_NAV , 'ERROR: Failed Open to Read NETCDF Global Attribute'
      print*, EXE_PROMPT_NAV ,' filename is: ', nc_file
      return
  endif

  status = nf90_get_att(ncid, nf90_global, trim(attr_name), attr_value)

  status = nf90_close(ncid)

 end subroutine read_netcdf_global_attribute_r4

 subroutine read_netcdf_global_attribute_char(nc_file, attr_name, attr_value)
  character(len=*), intent(in):: nc_file, attr_name
  character(len=*), intent(out):: attr_value
  integer:: ncid, status

   status = nf90_open(nc_file, mode = nf90_nowrite, ncid = ncid)

   if (status /= nf90_noerr) then
      print *, EXE_PROMPT_NAV , 'ERROR: Failed Open to Read NETCDF Global Attribute'
      print*, EXE_PROMPT_NAV ,' filename is: ', nc_file
      return
   endif

   status = nf90_get_att(ncid, nf90_global, trim(attr_name), attr_value)

   status = nf90_close(ncid)

 end subroutine read_netcdf_global_attribute_char
 !------------------------------------------------------------------------------
 ! read a real attribute from a variable
 ! this is redundant but more convenient
 !------------------------------------------------------------------------------
 subroutine read_netcdf_variable_attribute_real (nc_file, sds_name, attr_name, attr_value)
   character(len=*), intent(in):: nc_file, sds_name, attr_name
   real, intent(out) :: attr_value
   integer:: status, sd_id, sds_id

   status = nf90_open(nc_file, mode = nf90_nowrite, ncid = sd_id)
   status = nf90_inq_varid(sd_id,trim(sds_name), sds_id)
   status = nf90_get_att(sd_id, sds_id, trim(attr_name), attr_value)
   status = nf90_close(sd_id)
 
 end subroutine read_netcdf_variable_attribute_real

 !------------------------------------------------------------------------------
 ! read a real attribute from a variable
 !------------------------------------------------------------------------------
 subroutine read_netcdf_attribute_real (nc_file_id, attribute_id, var_name, attr)
   integer, intent(in):: nc_file_id
   integer, intent(in):: attribute_id
   character(len=*), intent(in):: var_name
   real, intent(out) :: attr
   integer:: status

   status = nf90_get_att(nc_file_id, attribute_id, trim(var_name), attr)
 
 end subroutine read_netcdf_attribute_real
 !------------------------------------------------------------------------------
 ! read a double attribute from a variable
 !------------------------------------------------------------------------------
 subroutine read_netcdf_attribute_double (nc_file_id, attribute_id, var_name, attr)
   integer, intent(in):: nc_file_id
   integer, intent(in):: attribute_id
   character(len=*), intent(in):: var_name
   real(kind=8), intent(out) :: attr
   integer:: status

   status = nf90_get_att(nc_file_id, attribute_id, trim(var_name), attr)
 
 end subroutine read_netcdf_attribute_double
 !------------------------------------------------------------------------------
 ! read ahi nav coeff
 !------------------------------------------------------------------------------
 subroutine read_ahi_time (nc_file_id, Line_Start, Line_Count, Line_Stride, Scanline_Time_Ms, Read_Time)
   integer, intent(in):: nc_file_id
   integer(kind=4), intent(in):: Line_Start, Line_Count, Line_Stride
   integer(kind=4), dimension(:), intent(out):: Scanline_Time_Ms
   logical, intent(out):: Read_Time
   real(kind=8), dimension(size(Scanline_TIme_Ms)):: var_output
   real(kind=8):: base_time
   character (len =20):: var_name, attr_name
   integer:: nc_var_id
   integer:: status
   integer, dimension(9):: time_values
   integer:: i

   Read_Time = .true.

   !read time (posix)
   var_name = 'line_time_offset'
   status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
   if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
   endif

   sds_start_1d = (/line_start/)
   sds_edge_1d = (/line_count/)
   sds_stride_1d = (/line_stride/)
   status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=sds_start_1d, &
                         count=sds_edge_1d, stride=sds_stride_1d)
   if (status /= nf90_noerr) then
        print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
        return
   endif

   attr_name = 'base_time'
   status = nf90_get_att(nc_file_id, nc_var_id, trim(attr_name), base_time)

   !--- convert to millisecond since midnight
   var_output = var_output + base_time

   do i = 1, Line_Count
    call gmtime(int(var_output(i)), time_values)
    Scanline_Time_Ms(i) = time_values(1) + time_values(2)*60.0 + time_values(3)*60.0*60.0
   enddo

   Read_Time = .false.

 end subroutine read_ahi_time
 !------------------------------------------------------------------------------
 ! read abi time
 !------------------------------------------------------------------------------
 subroutine read_abi_time (nc_file_id, Line_Start, Line_Count, Line_Count_Max, Line_Stride, Scanline_Time_Ms, Read_Time)
   integer, intent(in):: nc_file_id
   integer(kind=4), intent(in):: Line_Start, Line_Count, Line_Stride, Line_Count_Max
   integer(kind=4), dimension(:), intent(out):: Scanline_Time_Ms
   logical, intent(out):: Read_Time
   integer:: nc_var_id
   integer:: status
   real(kind=8), dimension(2):: time_bounds
   integer, dimension(9):: time_values
   real(kind=8):: time_temp, time_delta, time_offset, time_start
   integer:: i
   character(len=20):: var_name

   Read_Time = .true.

   Time_Offset = 946728000.0    !offset from GOES-R epoch to POSIX epoch

   !read time (goes-r epoch)
   var_name = 'time_bounds'
   status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
   if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
   endif

   sds_start_1d = (/1/)
   sds_edge_1d = (/2/)
   sds_stride_1d = (/1/)
   status = nf90_get_var(nc_file_id, nc_var_id, time_bounds, start=sds_start_1d, &
                         count=sds_edge_1d, stride=sds_stride_1d)
   if (status /= nf90_noerr) then
        print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
        return
   endif

   time_bounds = time_bounds + time_offset    !now standard unix epoch time

!_-- how to account for segment bounds vs disk bounds
   time_delta = (time_bounds(2) - time_bounds(1)) / Line_Count_Max / Line_Stride
   time_start = time_bounds(1) + Line_Start*time_delta
   do i = 1, Line_Count
      time_temp = time_start + (i-1) * time_delta
      call gmtime(int(time_temp), time_values)
      Scanline_Time_Ms(i) = 1000.0*(time_values(1) + time_values(2)*60.0 + time_values(3)*60.0*60.0)
   enddo
    
   Read_Time = .false.

 end subroutine read_abi_time

 !------------------------------------------------------------------------------
 ! read ahi nav coeff
 !------------------------------------------------------------------------------
 subroutine read_ahi_nav_coeff ( ncdf_file, CFAC, COFF, LFAC, LOFF)
 
 !, CFAC, COFF, LFAC, LOFF, Sub_point)
 
   character (len =*) , intent(in) :: ncdf_file
   REAL(KIND(0.0d0)) , INTENT (OUT) :: CFAC
   REAL(KIND(0.0d0)) , INTENT (OUT) :: COFF
   REAL(KIND(0.0d0)) , INTENT (OUT) :: LFAC
   REAL(KIND(0.0d0)) , INTENT (OUT) :: LOFF
!   REAL(KIND(0.0d0)) , INTENT (OUT) :: Sub_point
 
   integer :: nc_file_id
   integer :: nc_var_id
   integer :: status
   REAL :: attr
   
!open netCDF file
    status = nf90_open(path = trim(ncdf_file), mode = nf90_nowrite, ncid = nc_file_id)
    if (status .ne. nf90_noerr) then 
            print *, "Error: Unable to open netCDF file ", trim(ncdf_file)
            stop  
    endif  
    
    status = nf90_inq_varid(nc_file_id, trim('x_cgms'), nc_var_id)
    
    call read_netcdf_attribute_real(nc_file_id, nc_var_id, 'CFAC', attr)
    CFAC = DBLE(attr)
    call read_netcdf_attribute_real(nc_file_id, nc_var_id, 'COFF', attr)
    COFF = DBLE(attr)
 
    status = nf90_inq_varid(nc_file_id, trim('y_cgms'), nc_var_id)
    
    call read_netcdf_attribute_real(nc_file_id, nc_var_id, 'LFAC', attr)
    LFAC = DBLE(attr)
    call read_netcdf_attribute_real(nc_file_id, nc_var_id, 'LOFF', attr)
    LOFF = DBLE(attr)

!close netCDF file
    status = nf90_close(nc_file_id)
    if (status .ne. nf90_noerr) then 
            print *, "Error: Unable to open netCDF file ", trim(ncdf_file)
            stop  
    endif
 
 end subroutine  read_ahi_nav_coeff

   ! ----------------------------------------------------------
   ! Read and Uscale 2D arrays Two-Byte Integers
   ! ----------------------------------------------------------
   subroutine read_and_unscale_netcdf_1d (nc_file_id, var_start, var_stride, &
                                          var_dim, var_name, var_output_unscaled)
      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim
      character(len=*), intent(in) :: var_name

      real(kind=4) , dimension(var_dim(1)) :: var_output_scaled
      real(kind=4) , intent(out), dimension(:) :: var_output_unscaled
      real(kind=4) , dimension(var_dim(1)) :: var_output_unscaled_temp
      real(kind=4):: add_offset, scale_factor

      integer(kind=2):: Fill_Value

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      status = nf90_get_att(nc_file_id, nc_var_id, "add_offset", add_offset)
      if (status /= 0) add_offset = 0.0

      status = nf90_get_att(nc_file_id, nc_var_id, "scale_factor", scale_factor)
      if (status /= 0) scale_factor = 1.0

      status = nf90_get_att(nc_file_id, nc_var_id, "_FillValue", fill_value)
      if (status /= 0) fill_value = -999.0

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output_scaled, start=var_start, &
                            count=var_dim, stride=var_stride)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
             print *, "in read_and_unscale_netcdf_1d"
          print *, "var name = ", trim(var_name)
          print *, "var start = ", var_start
          print *, "var count = ", var_dim
          print *, "var stride = ", var_stride
          print *, "varname =  ",trim(var_name)
          print *, "shape var_output_scaled = ", shape(var_output_scaled)
            stop
            return
      endif

      !--- unscale
      where(var_output_scaled /= fill_value)
         var_output_unscaled_temp = add_offset + scale_factor * var_output_scaled
      elsewhere
         var_output_unscaled_temp = Missing_Value_Netcdf
      endwhere

      var_output_unscaled = Missing_Value_Netcdf
      var_output_unscaled(1:var_dim(1)) = var_output_unscaled_temp

   end subroutine read_and_unscale_netcdf_1d



   ! ----------------------------------------------------------
   ! Read and Uscale 2D arrays Two-Byte Integers
   ! ----------------------------------------------------------
   subroutine read_and_unscale_netcdf_2d (nc_file_id, var_start, var_stride, &
                                          var_dim, var_name, var_output_unscaled)
      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim
      character(len=*), intent(in) :: var_name

      real(kind=4) , dimension(var_dim(1),var_dim(2)) :: var_output_scaled
      real(kind=4) , intent(out), dimension(:,:) :: var_output_unscaled
      real(kind=4) , dimension(var_dim(1),var_dim(2)) :: var_output_unscaled_temp
      real(kind=4):: add_offset, scale_factor

      integer(kind=2):: Fill_Value

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      status = nf90_get_att(nc_file_id, nc_var_id, "add_offset", add_offset)
      if (status /= 0) add_offset = 0.0

      status = nf90_get_att(nc_file_id, nc_var_id, "scale_factor", scale_factor)
      if (status /= 0) scale_factor = 1.0

      status = nf90_get_att(nc_file_id, nc_var_id, "_FillValue", fill_value)
      if (status /= 0) fill_value = -999.0

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output_scaled, start=var_start, &
                            count=var_dim, stride=var_stride)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
             print *, "in read_and_unscale_netcdf_2d"
          print *, "var name = ", trim(var_name)
          print *, "var start = ", var_start
          print *, "var count = ", var_dim
          print *, "var stride = ", var_stride
          print *, "varname =  ",trim(var_name)
          print *, "shape var_output_scaled = ", shape(var_output_scaled)
            stop
            return
      endif

      !--- unscale
      where(var_output_scaled /= fill_value)
         var_output_unscaled_temp = add_offset + scale_factor * var_output_scaled
      elsewhere
         var_output_unscaled_temp = Missing_Value_Netcdf
      endwhere

      var_output_unscaled = Missing_Value_Netcdf
      var_output_unscaled(1:var_dim(1),1:var_dim(2)) = var_output_unscaled_temp

   end subroutine read_and_unscale_netcdf_2d

 
  
   ! ----------------------------------------------------------
   ! Read in 1D Real arrays
   ! ----------------------------------------------------------
   subroutine read_netcdf_1d_real (nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)
                                  
      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim
      character(len=*), intent(in) :: var_name
      real, intent(out), dimension(:) :: var_output

      integer :: nc_var_id
      integer :: status


      status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, &
                            count=var_dim, stride=var_stride)
      if (status /= nf90_noerr) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

   end subroutine read_netcdf_1d_real                                                                                                                           

   ! ----------------------------------------------------------
   ! Read in 1D Double arrays
   ! ----------------------------------------------------------
   subroutine read_netcdf_1d_double (nc_file_id, var_start, var_stride, &
                                     var_dim, var_name, var_output)
      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim
      character(len=*), intent(in) :: var_name
      real (kind=8), intent(out), dimension(:) :: var_output

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, &
                            count=var_dim, stride=var_stride)
      if (status /= nf90_noerr) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

   end subroutine read_netcdf_1d_double                                                                                                                  

   ! ----------------------------------------------------------
   ! Read in 1D integer arrays 
   ! ----------------------------------------------------------
   subroutine read_netcdf_1d_int (nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)
      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim
      character(len=*), intent(in) :: var_name
      integer, intent(out), dimension(:) :: var_output

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then 
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, &
                            count=var_dim, stride=var_stride)
                            
      if (status /= nf90_noerr) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

   end subroutine read_netcdf_1d_int

   ! ----------------------------------------------------------
   ! Read in 2D arrays
   ! ----------------------------------------------------------
   subroutine read_netcdf_2d_real (nc_file_id, var_start, var_stride, &
                                   var_dim, var_name, var_output)
      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      real, intent(out), dimension(:,:) :: var_output

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)

      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, &
                            count=var_dim, stride=var_stride)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

   end subroutine read_netcdf_2d_real

   ! ----------------------------------------------------------
   ! Read in 2D arrays Integers
   ! ----------------------------------------------------------
   subroutine read_netcdf_2d_int4 (nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)
      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      integer(kind=4) , intent(out), dimension(:,:) :: var_output

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, &
                            count=var_dim, stride=var_stride)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

   end subroutine read_netcdf_2d_int4

   ! ----------------------------------------------------------
   ! Read in 2D arrays Integers
   ! ----------------------------------------------------------
   subroutine read_netcdf_2d_int1 (nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)
      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      integer(kind=1) , intent(out), dimension(:,:) :: var_output

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, &
                            count=var_dim, stride=var_stride)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

   end subroutine read_netcdf_2d_int1

   ! ----------------------------------------------------------
   ! Read in 2D arrays Integers
   ! ----------------------------------------------------------
   subroutine read_netcdf_2d_int2 (nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)
      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      integer(kind=2) , intent(out), dimension(:,:) :: var_output

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, &
                            count=var_dim, stride=var_stride)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

   end subroutine read_netcdf_2d_int2


   ! ----------------------------------------------------------
   ! Read in 2D arrays Characters
   ! ----------------------------------------------------------
   subroutine read_netcdf_2d_char (nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)
      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      character(len=30) , intent(out), dimension(:,:) :: var_output
      character(len=30), allocatable, dimension(:,:) :: var

      integer :: nc_var_id
      integer :: status, tmp1, tmp2, i
      integer, dimension(2) ::dimIDs

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !find dimentions
      status = nf90_inquire_variable(nc_file_id, nc_var_id, dimids = dimIDs)
      status = nf90_inquire_dimension(nc_file_id, dimIDs(1), len = tmp1)
      status = nf90_inquire_dimension(nc_file_id, dimIDs(2), len = tmp2)
      allocate (var(tmp1,tmp2))

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var, start=(/1,1/), count=(/tmp1,tmp2/) )
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

      !extract and save classifier names to the final array
      do i = 1, tmp2
        if ((var(i,1) .ge. 'a' .and. var(i,1) .le. 'z') &
        .or.(var(i,1) .ge. 'A' .and. var(i,1) .le. 'Z')) then 
           var_output(i,:) = trim(var(i,1))
        endif
      enddo

      if (allocated(var)) deallocate (var)


   end subroutine read_netcdf_2d_char

   ! ----------------------------------------------------------
   ! Read in 3D arrays
   ! ----------------------------------------------------------
   subroutine read_netcdf_3d_real (nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)

      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      real, intent(out), dimension(:,:,:) :: var_output

      integer :: nc_var_id
      integer :: status = 0

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, count=var_dim)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif


   end subroutine read_netcdf_3d_real

   ! ----------------------------------------------------------
   ! Read in 4D arrays
   ! ----------------------------------------------------------
   subroutine read_netcdf_4d_real (nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)

      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim
      character(len=*), intent(in) :: var_name
      real, intent(out), dimension(:,:,:,:) :: var_output

      integer :: nc_var_id
      integer :: status = 0

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, count=var_dim)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

   end subroutine read_netcdf_4d_real

   subroutine read_netcdf_4d_int4 (nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)

      integer, intent(in) :: nc_file_id
      integer, dimension(:), intent(in) :: var_start
      integer, dimension(:), intent(in) :: var_stride
      integer, dimension(:), intent(in) :: var_dim
      character(len=*), intent(in) :: var_name
      integer(kind=4), intent(out), dimension(:,:,:,:) :: var_output

      integer :: nc_var_id
      integer :: status = 0

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, count=var_dim)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

   end subroutine read_netcdf_4d_int4

   !------------------------------------------------------------------------------
   !  used in static nav
   !------------------------------------------------------------------------------
   subroutine read_netcdf_dimension(nc_file, dim_name, dim_value)
     character(len=*), intent(in):: nc_file, dim_name
     integer(kind=4), intent(out):: dim_value
     integer:: ncid, status, dimid

     status = nf90_open(nc_file, mode = nf90_nowrite, ncid = ncid)
     if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed Open to Read NETCDF Dimension'
       print*, EXE_PROMPT_NAV ,' filename is: ', nc_file
       return
     endif

     status = nf90_inq_dimid(ncid, dim_name, dimid)
     status = nf90_inquire_dimension(ncid, dimid, len=dim_value)

     status = nf90_close(ncid)
  end subroutine read_netcdf_dimension

  !------------------------------------------------------------------------------
  ! read abi kappa0 
  !------------------------------------------------------------------------------
  subroutine read_abi_kappa0 (nc_file_id, var_output)
    integer, intent(in):: nc_file_id
    real, intent(out) :: var_output
    character(len=6) var_name
    integer :: nc_var_id
    integer :: status

    var_name='kappa0'

    status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
    if (status /= nf90_noerr) then
      print *, "Error: Unable to get variable id for ", trim(var_name)
      return
    endif

    ! Get kappa0 value.
    status = nf90_get_var(nc_file_id, nc_var_id, var_output)
    if (status /= nf90_noerr) then
      print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
      return
    endif

  end subroutine read_abi_kappa0

  !------------------------------------------------------------------------------
  ! read abi maximum_focal_plane_temperature
  !------------------------------------------------------------------------------
  subroutine read_abi_max_focal_plane_temp (nc_file_id, var_name, var_output)
    integer, intent(in):: nc_file_id
    character(len=31), intent(in):: var_name
    real, intent(out) :: var_output
    integer :: nc_var_id
    integer :: status

    status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
    if (status /= nf90_noerr) then
      print *, "Error: Unable to get variable id for ", trim(var_name)
      return
    endif

    ! Get maximum_focal_plane_temperature value.
    status = nf90_get_var(nc_file_id, nc_var_id, var_output)
    if (status /= nf90_noerr) then
      print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
      return
    endif

  end subroutine read_abi_max_focal_plane_temp

!-------------------------------------------------------------------------
! routine to write global attributes to clavrx files
!
!-------------------------------------------------------------------------
 subroutine WRITE_CLAVRX_NETCDF_GLOBAL_ATTRIBUTES(netcdf_file_id)

 integer(kind=int4), intent(in):: netcdf_file_id

 integer:: libver
 integer:: major_v, minor_v, release
 integer(kind=int2):: Istatus = 0
 character(len=80) :: hdf_ver
 character(len=36) :: machine


 character(len=6), parameter:: Conventions_String = "CF-1.6"
 character(len=38), parameter:: Metadata_Conventions_String = "CF-1.6, Unidata Dataset Discovery v1.0"
 character(len=42), parameter:: Standard_Name_Vocabulary_String = "CF Standard Name Table (v25, 05 July 2013)"
 character(len=13), parameter:: Naming_Authority_String = "gov.noaa.ncdc"
 character(len=37), parameter:: License_String = "No constraints on data access or use."


! Definition of strings used as HDF attributes.
!
! version history 4.1 - delivered to OSDPD in November 2006
! version history 4.2 - demonstrated on METOP
! version history 4.3 - included surface emissivity fields
! version history 4.4 - included lrc
! version history 5.0 - included modis white sky and ash protoype 
! version history 5.1 - rtm structures now 101 levels and 
!                       reorganized level-1b ingest to 
!                       read segment all at once prior to processing
!                       first version with working volcanic ash
! version history 5.2    bayesian cloud mask and DCOMP
! version history 6.0  - MODIS capability begin
! version history 6.5  - VIIRS capability begin

 character(len=36), parameter :: creator0 = "CLAVR-x + PATMOS-x "

 character(len=100):: Title_String
 character(len=100):: Calibration_String
 character(len=100):: Product_Version_String
 character(len=100):: Status_String
 character(len=100):: Institution_String
 character(len=100):: Program_String
 character(len=500):: Summary_String 
 character(len=200):: Variable_String
 character(len=500):: Keywords_String
 character(len=200):: Keywords_Vocabulary_String
 character(len=100):: Time_Coverage_Resolution_String
 character(len=100):: Metadata_Link_String
 character(len=100):: Spatial_Resolution_String

 include 'version.inc'

 !complete the creator string with the version number
 Clavrx_Global_Attr%creator = trim(creator0)//trim(Product_Version_String)

 call getenv ("HOST",Clavrx_Global_Attr%machine)

 if (len_trim(Clavrx_Global_Attr%machine) == 0) call getenv ("HOSTNAME",Clavrx_Global_Attr%machine)
 if (len_trim(Clavrx_Global_Attr%machine) == 0) Clavrx_Global_Attr%machine = "unknown"

 Clavrx_Global_Attr%plang = "F90"
 Clavrx_Global_Attr%timestamp = trim(netcdf_timestamp())

!--- determine HDF library version
!if (hglibver(major_v, minor_v, release, Clavrx_Global_Attr%hdf_ver) /= SUCCEED) then
!   print *, "could not determine HDF library version"
!   stop 961
!end if

Istatus = 0

!---- describe CLAVR-x
Clavrx_Global_Attr%hdf_ver = "unknown"
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "HDF_LIB_VERSION", trim(Clavrx_Global_Attr%hdf_ver))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "MACHINE", trim(Clavrx_Global_Attr%machine))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "PROGLANG", trim(Clavrx_Global_Attr%plang))+Istatus

!--- CF compliant global attributes required for NCDC delivery
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "date_created", trim(Clavrx_Global_Attr%timestamp))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "product_version", trim(Product_Version_String))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"summary", trim(Summary_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"cdr_variable",trim(Variable_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"institution",trim(Institution_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"cdr_program",trim(Program_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"title",trim(Title_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"calibration_version",trim(Calibration_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"keywords",trim(Keywords_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"keywords_vocabulary",trim(Keywords_Vocabulary_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"time_coverage_resolution",trim(Time_Coverage_Resolution_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"metadata_link", trim(Metadata_Link_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"spatial_resolution", trim(Spatial_Resolution_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"Conventions", trim(Conventions_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"title", trim(Title_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"Metadata_Conventions", trim(Metadata_Conventions_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"standard_name_vocabulary",trim(Standard_Name_Vocabulary_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"naming_authority",trim(Naming_Authority_String)) + Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global,"license",trim(License_String)) + Istatus

!--- describe the data
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "sensor", trim(Clavrx_Global_Attr%sensor_name))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "platform", trim(Clavrx_Global_Attr%platform_name))+Istatus

Istatus = nf90_put_att(netcdf_file_id,nf90_global, "FILENAME", trim(Clavrx_Global_Attr%file_name))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "L1B", trim(Clavrx_Global_Attr%file_1b))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "RESOLUTION_KM", Clavrx_Global_Attr%resolution_km)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "START_YEAR", Clavrx_Global_Attr%start_year)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "START_DAY", Clavrx_Global_Attr%start_day)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "START_TIME", Clavrx_Global_Attr%start_time)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "END_YEAR", Clavrx_Global_Attr%end_year)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "END_DAY", Clavrx_Global_Attr%end_day)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "END_TIME", Clavrx_Global_Attr%end_time)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "CLOUD_MASK_MODE", trim(Clavrx_Global_Attr%mask_mode))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "ACHA_MODE_FINAL", trim(Clavrx_Global_Attr%acha_mode))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "ACHA_MODE_USER", trim(Clavrx_Global_Attr%acha_mode_user))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "DCOMP_MODE", Clavrx_Global_Attr%dcomp_mode)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "WMO_SATELLITE_CODE", Clavrx_Global_Attr%wmo_sc_code)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "REFL_0_65UM_NOM_DARK_COMPOSITE_NAME", trim(Clavrx_Global_Attr%dark_name))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "NAIVE_BAYESIAN_CLOUD_MASK_NAME",trim(Clavrx_Global_Attr%mask_name))+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "GEO_SUB_LON", Clavrx_Global_Attr%geo_sub_lon)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "GEO_SUB_LAT", Clavrx_Global_Attr%geo_sub_lat)+Istatus

!--- describe subsetting
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "SUBSET_FLAG", Clavrx_Global_Attr%subset_flag)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "LAT_SOUTH_SUBSET", Clavrx_Global_Attr%lat_south_subset)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "LAT_NORTH_SUBSET", Clavrx_Global_Attr%lat_north_subset)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "LON_WEST_SUBSET", Clavrx_Global_Attr%lon_west_subset)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "LON_EAST_SUBSET", Clavrx_Global_Attr%lon_east_subset)+Istatus
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "SUBSET_NAME", trim(Clavrx_Global_Attr%subset_name))+Istatus

!--- data type
Istatus = nf90_put_att(netcdf_file_id,nf90_global, "DATA_TYPE", trim(Clavrx_Global_Attr%data_type))+Istatus


!--- NCDC Attributes

!-- other global attributes for the level3 file.
     if (Clavrx_Global_Attr%data_type(1:4) == "GRID") then
      Istatus = nf90_put_att(netcdf_file_id,nf90_global, "NUM_CELLS_TOTAL", Clavrx_Global_Attr%num_cells)+Istatus
      Istatus = nf90_put_att(netcdf_file_id,nf90_global, "NUM_CELLS_WITH_DATA", Clavrx_Global_Attr%num_cells_with_data)+Istatus
      Istatus = nf90_put_att(netcdf_file_id,nf90_global, "GRIDCELL_RESOLUTION", Clavrx_Global_Attr%dlat)+Istatus
      Istatus = nf90_put_att(netcdf_file_id,nf90_global, "GRIDCELL_RESOLUTION_UNIT", "degree")+Istatus
      if (Clavrx_Global_Attr%grid_format(1:10) == "EQUAL_AREA") then
       Istatus = nf90_put_att(netcdf_file_id,nf90_global, "GRIDCELL_FORMAT", "EQUAL_AREA")+Istatus
      else
       Istatus = nf90_put_att(netcdf_file_id,nf90_global, "GRIDCELL_FORMAT", "EQUAL_ANGLE")+Istatus
      endif
    endif

!---- processing flags
 Istatus = nf90_put_att(netcdf_file_id, nf90_global, "USE_1B_THERMAL_CALIBRATION_FLAG", Clavrx_Global_Attr%therm_cal_1b)+Istatus
 Istatus = nf90_put_att(netcdf_file_id, nf90_global, "USE_1B_REFLECTANCE_CALIBRATION_FLAG",Clavrx_Global_Attr%Ref_cal_1b)+Istatus
 Istatus = nf90_put_att(netcdf_file_id, nf90_global, "RENAVIGATION_FLAG", Clavrx_Global_Attr%nav_opt)+Istatus
 Istatus = nf90_put_att(netcdf_file_id, nf90_global, "USE_SST_ANALYSIS_FLAG", Clavrx_Global_Attr%use_sst_anal)+Istatus
 Istatus = nf90_put_att(netcdf_file_id, nf90_global, "NWP_OPT", Clavrx_Global_Attr%nwp_opt)+Istatus
 Istatus = nf90_put_att(netcdf_file_id, nf90_global, "MODIS_CLEAR_SKY_REFLECTANCE_FLAG", Clavrx_Global_Attr%modis_clr_alb_flag)+Istatus

!-- reflectance channel calibration
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH1_GAIN_LOW", Clavrx_Global_Attr%ch1_gain_low)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH1_GAIN_HIGH", Clavrx_Global_Attr%ch1_gain_high)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH1_SWITCH_COUNT", Clavrx_Global_Attr%ch1_switch_count)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH1_DARK_COUNT", Clavrx_Global_Attr%ch1_dark_count)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH2_GAIN_LOW", Clavrx_Global_Attr%ch2_gain_low)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH2_GAIN_HIGH", Clavrx_Global_Attr%ch2_gain_high)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH2_SWITCH_COUNT", Clavrx_Global_Attr%ch2_switch_count)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH2_DARK_COUNT", Clavrx_Global_Attr%ch2_dark_count)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH3A_GAIN_LOW", Clavrx_Global_Attr%ch3a_gain_low)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH3A_GAIN_HIGH", Clavrx_Global_Attr%ch3a_gain_high)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH3A_SWITCH_COUNT", Clavrx_Global_Attr%ch3a_switch_count)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "CH3A_DARK_COUNT", Clavrx_Global_Attr%ch3a_dark_count)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "SUN_EARTH_DISTANCE", Clavrx_Global_Attr%sun_earth_distance)+Istatus

!--- thermal calibration constants
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "C1", Clavrx_Global_Attr%c1)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "C2", Clavrx_Global_Attr%c2)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "A_20", Clavrx_Global_Attr%a_20)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "B_20", Clavrx_Global_Attr%b_20)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "NU_20",Clavrx_Global_Attr%nu_20)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "A_31", Clavrx_Global_Attr%a_31)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "B_31", Clavrx_Global_Attr%b_31)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "NU_31",Clavrx_Global_Attr%nu_31)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "A_32", Clavrx_Global_Attr%a_32)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "B_32", Clavrx_Global_Attr%b_32)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "NU_32", Clavrx_Global_Attr%nu_32)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "SOLAR_20_NU", Clavrx_Global_Attr%solar_Ch20_nu)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "TIME_ERROR_SECONDS",Clavrx_Global_Attr%timerr_seconds)+Istatus
Istatus = nf90_put_att(netcdf_file_id, nf90_global, "LWIR_Focal_Plane_Temperature", Clavrx_Global_Attr%LWIR_Focal_Plane_Temperature)+Istatus

if (Istatus /= 0) then
  print *, "error writing run-control flags as global NETCDF attributes"
  print *, " Most likely reason is that Level2 Output Path is not created, or only with read access"
  stop 963
endif

contains

  function netcdf_timestamp() result(string)
     character (len = 36) ::string
     character(len=10), dimension(3):: ctime

     call date_and_time(ctime(1), ctime(2), ctime(3))

     ! Timestamp string format in accordance with the ISO-8601 standard.
     string = ctime(1)(1:4)//"-"//ctime(1)(5:6)//"-"//ctime(1)(7:8)&
                //"T"//ctime(2)(1:2)//":"//ctime(2)(3:4)//":"//ctime(2)(5:6)&
                //ctime(3)(1:3)//":"//ctime(3)(4:5)
     return
   end function netcdf_timestamp

 end subroutine WRITE_CLAVRX_NETCDF_GLOBAL_ATTRIBUTES
 
end module CX_NETCDF4_MOD
