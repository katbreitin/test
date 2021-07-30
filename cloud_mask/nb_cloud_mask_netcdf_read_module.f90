!----------------------------------------------------------------------
! MODULE name: NETCDF_READ_MODULE
! 
! Routines for reading in netCDF files, such as AHI data
!
! Authors: Andrew Heidinger, NOAA/NESDIS
!          Andi Walther, CIMSS
!          Denis Botambekov, CIMSS
!
! DEPENDENCIES: Constants, NETCDF library
!
!----------------------------------------------------------------------
MODULE NB_CLOUD_MASK_NETCDF_READ_MODULE

 use NETCDF

 implicit none


 PUBLIC :: OPEN_NETCDF
 PUBLIC :: CLOSE_NETCDF
 PUBLIC :: READ_NETCDF_GLOBAL_ATTRIBUTE_I4
 PUBLIC :: READ_NETCDF_GLOBAL_ATTRIBUTE_R4
 PUBLIC :: READ_NETCDF_GLOBAL_ATTRIBUTE_CHAR
 PUBLIC :: GET_GROUP_ID
 PUBLIC :: READ_NETCDF_2D_CHAR 
 PUBLIC :: READ_NETCDF_1D_INT
 PUBLIC :: READ_NETCDF_1D_REAL
 PUBLIC :: READ_NETCDF_2D_REAL
 PUBLIC :: READ_NETCDF_2D_INT
 PUBLIC :: READ_NETCDF_3D_REAL
 PUBLIC :: READ_NETCDF_4D_REAL

 !!!! THIS IS TO WRITE OUTPUT NetCDF SHOULDN't BE IN CLAVR-x
 PUBLIC :: OPEN_NETCDF_WRITE
 PUBLIC :: DEFINE_SDS_DIMS_NETCDF
 PUBLIC :: WRITE_SDS_NETCDF_2D_R4
 PUBLIC :: WRITE_SDS_NETCDF_2D_I4
 PUBLIC :: WRITE_SDS_NETCDF_3D_R4

 CHARACTER(len=120), PARAMETER, PRIVATE:: EXE_PROMPT_NAV = "CLAVR-x NetCdf Read Module >> "
 REAL, PARAMETER, PRIVATE:: Missing_Value_Netcdf = -999.0
 INTEGER, PUBLIC:: DEFAULT_NAME_LENGTH

 contains
 !------------------------------------------------------------------------------
 ! read in a dimension  - used for assumed string length of names
 !------------------------------------------------------------------------------
 subroutine read_netcdf_dimension(ncid, dim_name, dim_value)
  integer, intent(in) :: ncid
  character(len=*), intent(in):: dim_name
  integer(kind=4),intent(out):: dim_value
  integer :: dimid
  integer:: status, sds_id

  status = nf90_inq_dimid(ncid,trim(dim_name),dimid)

  status = nf90_inquire_dimension(ncid,dimid,len=dim_value) 

  if (status /= nf90_noerr) then
     print *, EXE_PROMPT_NAV , 'ERROR: Failed to Get Dimension = ', trim(dim_name)
     return
  endif

 end subroutine read_netcdf_dimension


 !------------------------------------------------------------------------------
 ! open netcdf file for writing
 !------------------------------------------------------------------------------
 SUBROUTINE OPEN_NETCDF_WRITE(nc_file, ncid)
    CHARACTER(len=*), INTENT(in):: nc_file
    INTEGER, INTENT(out):: ncid
    INTEGER:: status
    !status = nf90_open(trim(nc_file), mode = nf90_write, ncid = ncid)
    status = nf90_create(trim(nc_file), IOR(NF90_NETCDF4, NF90_CLOBBER), ncid = ncid) 
    if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed Open on NetCdf File'
       print*, EXE_PROMPT_NAV ,' filename is: ', nc_file
       return
    endif
 END SUBROUTINE OPEN_NETCDF_WRITE

 !------------------------------------------------------------------------------
 !  define NETCDF sds
 !------------------------------------------------------------------------------
 SUBROUTINE DEFINE_SDS_DIMS_NETCDF(Sd_Id, Nx, Ny, N_Class, Dimids, Istatus)
   INTEGER,INTENT(in) :: Sd_Id, Nx, Ny, N_Class
   INTEGER, DIMENSION(3), INTENT(out) :: Dimids
   INTEGER,INTENT(out) :: Istatus

   INTEGER :: X_Dimid, Y_Dimid, Class_Dimid

   Istatus = 0 

   Istatus = nf90_def_dim(Sd_Id, "Nx", Nx, X_Dimid) + Istatus
   Istatus = nf90_def_dim(Sd_Id, "Ny", Ny, Y_Dimid) + Istatus
   Istatus = nf90_def_dim(Sd_Id, "N_Class", N_Class, Class_Dimid) + Istatus

   Dimids = (/X_Dimid, Y_Dimid, Class_Dimid/)

 END SUBROUTINE DEFINE_SDS_DIMS_NETCDF

 !------------------------------------------------------------------------------
 !  NETCDF Write Routine 2D Integer4
 !------------------------------------------------------------------------------
 SUBROUTINE WRITE_SDS_NETCDF_2D_I4(Sd_Id, Nx, Ny, Dimids, Sds_Name, Sds_Data_2d_I4, Istatus)

   INTEGER,INTENT(in):: Sd_Id, Nx, Ny
   INTEGER,DIMENSION(2),INTENT(in):: Dimids
   CHARACTER(len=*),INTENT(in):: Sds_Name
   INTEGER, DIMENSION(:,:),INTENT(in):: Sds_Data_2d_I4
   INTEGER,INTENT(out):: Istatus

   INTEGER, DIMENSION(2) :: Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d
   INTEGER :: Sds_Id

   Istatus = 0
   Sds_Start_2d = 1
   Sds_Stride_2d = 1
   Sds_Edge_2d(1) = Nx
   Sds_Edge_2d(2) = Ny

   Istatus = nf90_def_var(Sd_Id, Sds_Name, NF90_INT, Dimids, Sds_Id) + Istatus

   Istatus = nf90_ENDdef(Sd_Id) + Istatus

   Istatus = nf90_put_var(Sd_Id, Sds_Id, &
                            Sds_Data_2d_I4, &
                            start = Sds_Start_2d, &
                            count = Sds_Edge_2d,  &
                            stride = Sds_Stride_2d) + Istatus

   Istatus = nf90_redef(Sd_Id) + Istatus

 END SUBROUTINE WRITE_SDS_NETCDF_2D_I4

 !------------------------------------------------------------------------------
 !  NETCDF Write Routine 2D Real4
 !------------------------------------------------------------------------------
 SUBROUTINE WRITE_SDS_NETCDF_2D_R4(Sd_Id, Nx, Ny, Dimids, Sds_Name, Sds_Data_2d_R4, Istatus)

   INTEGER,INTENT(in):: Sd_Id, Nx, Ny
   INTEGER,DIMENSION(2),INTENT(in):: Dimids
   CHARACTER(len=*),INTENT(in):: Sds_Name
   REAL, DIMENSION(:,:),INTENT(in):: Sds_Data_2d_R4
   INTEGER,INTENT(out):: Istatus

   INTEGER, DIMENSION(2) :: Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d
   INTEGER, save :: Sds_Id

   Istatus = 0
   Sds_Start_2d = 1
   Sds_Stride_2d = 1
   Sds_Edge_2d(1) = Nx
   Sds_Edge_2d(2) = Ny

   Istatus = nf90_def_var(Sd_Id, Sds_Name, NF90_REAL4, Dimids, Sds_Id) + Istatus
 
   Istatus = nf90_ENDdef(Sd_Id) + Istatus

   Istatus = nf90_put_var(Sd_Id, Sds_Id, &
                            Sds_Data_2d_R4, &
                            start = Sds_Start_2d, &
                            count = Sds_Edge_2d,  &
                            stride = Sds_Stride_2d) + Istatus

   Istatus = nf90_redef(Sd_Id) + Istatus

 END SUBROUTINE WRITE_SDS_NETCDF_2D_R4

 !------------------------------------------------------------------------------
 !  NETCDF Write Routine 3D Real4
 !------------------------------------------------------------------------------
 SUBROUTINE WRITE_SDS_NETCDF_3D_R4(Sd_Id, Nx, Ny, NClass, Dimids, Sds_Name, Sds_Data_3d_R4, Istatus)

   INTEGER,INTENT(in):: Sd_Id, Nx, Ny, NClass
   INTEGER,DIMENSION(3),INTENT(in):: Dimids
   CHARACTER(len=*),INTENT(in):: Sds_Name
   REAL, DIMENSION(:,:,:),INTENT(in):: Sds_Data_3d_R4
   INTEGER,INTENT(out):: Istatus

   INTEGER, DIMENSION(3) :: Sds_Start_3d, Sds_Stride_3d, Sds_Edge_3d
   INTEGER, save :: Sds_Id

   Istatus = 0
   Sds_Start_3d = 1
   Sds_Stride_3d = 1
   Sds_Edge_3d(1) = Nx
   Sds_Edge_3d(2) = Ny
   Sds_Edge_3d(3) = NClass

   Istatus = nf90_def_var(Sd_Id, Sds_Name, NF90_REAL4, Dimids, Sds_Id) + Istatus

   Istatus = nf90_ENDdef(Sd_Id) + Istatus

   Istatus = nf90_put_var(Sd_Id, Sds_Id, &
                            Sds_Data_3d_R4, &
                            start = Sds_Start_3d, &
                            count = Sds_Edge_3d,  &
                            stride = Sds_Stride_3d) + Istatus

   Istatus = nf90_redef(Sd_Id) + Istatus

 END SUBROUTINE WRITE_SDS_NETCDF_3D_R4

 !------------------------------------------------------------------------------
 ! open netcdf file for reading
 !------------------------------------------------------------------------------
 SUBROUTINE OPEN_NETCDF(nc_file, ncid)
    CHARACTER(len=*), INTENT(in):: nc_file
    INTEGER, INTENT(out):: ncid
    INTEGER:: status
    status = nf90_open(trim(nc_file), mode = nf90_nowrite, ncid = ncid)
    if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed Open on NetCdf File'
       print*, EXE_PROMPT_NAV ,' filename is: ', nc_file
       return
    endif
 END SUBROUTINE OPEN_NETCDF

 !------------------------------------------------------------------------------
 ! close netcdf file 
 !------------------------------------------------------------------------------
 SUBROUTINE CLOSE_NETCDF(ncid)
    INTEGER,INTENT(in):: ncid
    INTEGER:: status

    status = nf90_close(ncid)

    if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed Close on NetCdf File'
       return
    endif

 END SUBROUTINE CLOSE_NETCDF

 !------------------------------------------------------------------------------
 ! read global/group attribute INTEGER
 !
 ! if print_warning_flag = 1, an error will be printed if an error is found
 ! for ecm, pass the class_idx so these errors are printed only on the first
 ! classifier read.  all classifiers should have same attributes
 !------------------------------------------------------------------------------
 SUBROUTINE READ_NETCDF_GLOBAL_ATTRIBUTE_I4(ncid, attr_name, attr_value, print_warning_flag)
  INTEGER(kind=4), INTENT(in):: ncid
  CHARACTER(len=*), INTENT(in):: attr_name
  INTEGER(kind=4), INTENT(out):: attr_value
  INTEGER(kind=4), INTENT(in), OPTIONAL:: print_warning_flag
  INTEGER:: status

  status = nf90_get_att(ncid, nf90_global, attr_name, attr_value)

  if (status /= nf90_noerr) then
     if (print_warning_flag == 1) print *, EXE_PROMPT_NAV , 'ERROR: Reading NETCDF Attribute: ',trim(attr_name)
     attr_value = int(Missing_Value_Netcdf)
     return
  endif

 END SUBROUTINE READ_NETCDF_GLOBAL_ATTRIBUTE_I4

 !------------------------------------------------------------------------------
 ! read global/group attribute REAL
 !
 ! if print_warning_flag = 1, an error will be printed if an error is found
 !------------------------------------------------------------------------------
 SUBROUTINE READ_NETCDF_GLOBAL_ATTRIBUTE_R4(ncid, attr_name, attr_value, print_warning_flag)
  INTEGER(kind=4), INTENT(in):: ncid
  CHARACTER(len=*), INTENT(in):: attr_name
  REAL(kind=4), INTENT(out):: attr_value
  INTEGER(kind=4), INTENT(in), OPTIONAL:: print_warning_flag
  INTEGER:: status

  status = nf90_get_att(ncid, nf90_global, attr_name, attr_value)

  if (status /= nf90_noerr) then
      if (print_warning_flag == 1) print *, EXE_PROMPT_NAV , 'ERROR: Reading NETCDF Attribute: ',trim(attr_name)
      attr_value = Missing_Value_Netcdf
     return
  endif

 END SUBROUTINE READ_NETCDF_GLOBAL_ATTRIBUTE_R4

 !------------------------------------------------------------------------------
 ! read global/group attribute CHARACTER
 !
 ! if print_warning_flag = 1, an error will be printed if an error is found
 !------------------------------------------------------------------------------
 SUBROUTINE READ_NETCDF_GLOBAL_ATTRIBUTE_CHAR(ncid, attr_name, attr_value, print_warning_flag)
  INTEGER(kind=4), INTENT(in):: ncid
  CHARACTER(len=*), INTENT(in):: attr_name
  CHARACTER(len=DEFAULT_NAME_LENGTH), INTENT(out):: attr_value
  INTEGER(kind=4), INTENT(in), OPTIONAL:: print_warning_flag
  INTEGER:: status

   status = nf90_get_att(ncid, nf90_global, attr_name, attr_value)

   if (status /= nf90_noerr) then
      if (print_warning_flag == 1) print *, EXE_PROMPT_NAV , 'ERROR: Reading NETCDF Attribute: ',trim(attr_name)
      attr_value = ""
      return
   endif

 END SUBROUTINE READ_NETCDF_GLOBAL_ATTRIBUTE_CHAR


 !------------------------------------------------------------------------------
 ! get group id
 !------------------------------------------------------------------------------
 SUBROUTINE GET_GROUP_ID(ncid, group_name, group_id)
    INTEGER, INTENT(in) :: ncid
    CHARACTER(len=*), INTENT(in) :: group_name
    INTEGER, INTENT(out) :: group_id
    INTEGER:: status

    Status = nf90_inq_ncid(ncid, trim(group_name), group_id)
    if (status /= nf90_noerr) then
       print *, EXE_PROMPT_NAV , 'ERROR: Failed to Get GROUP_ID = ',trim(group_name)
       return
    endif
  END SUBROUTINE GET_GROUP_ID


   ! ----------------------------------------------------------
   ! Read in 1D Real arrays
   ! ----------------------------------------------------------
   SUBROUTINE READ_NETCDF_1d_REAL(nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)

      implicit none
      INTEGER, INTENT(in) :: nc_file_id
      INTEGER, DIMENSION(:), INTENT(in) :: var_start
      INTEGER, DIMENSION(:), INTENT(in) :: var_stride
      INTEGER, DIMENSION(:), INTENT(in) :: var_dim
      CHARACTER(len=*), INTENT(in) :: var_name
      REAL, INTENT(out), DIMENSION(:) :: var_output

      INTEGER :: nc_var_id
      INTEGER :: status


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

   END SUBROUTINE READ_NETCDF_1d_REAL

   ! ----------------------------------------------------------
   ! Read in 1D Integer arrays
   ! ----------------------------------------------------------
   SUBROUTINE READ_NETCDF_1D_INT(nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)

      implicit none
      INTEGER, INTENT(in) :: nc_file_id
      INTEGER, DIMENSION(:), INTENT(in) :: var_start
      INTEGER, DIMENSION(:), INTENT(in) :: var_stride
      INTEGER, DIMENSION(:), INTENT(in) :: var_dim
      CHARACTER(len=*), INTENT(in) :: var_name
      INTEGER, INTENT(out), DIMENSION(:) :: var_output

      INTEGER :: nc_var_id
      INTEGER :: status


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

   END SUBROUTINE READ_NETCDF_1D_INT

   ! ----------------------------------------------------------
   ! Read in 2D arrays Characters
   ! ----------------------------------------------------------

   SUBROUTINE READ_NETCDF_2D_CHAR(nc_file_id, start_var, var_name, var_output)
      implicit none
      INTEGER, INTENT(in) :: nc_file_id
      INTEGER, INTENT(in) :: start_var(:)

      CHARACTER(len=*), INTENT(in) :: var_name
      CHARACTER(len=DEFAULT_NAME_LENGTH) , INTENT(out), DIMENSION(:) :: var_output
      CHARACTER(len=DEFAULT_NAME_LENGTH), allocatable, DIMENSION(:,:) :: var

      INTEGER :: nc_var_id
      INTEGER :: status, tmp1, tmp2, i
      INTEGER, DIMENSION(2) ::dimIDs

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !find dimentions
      status = nf90_inquire_variable(nc_file_id, nc_var_id, dimids = dimIDs)
      status = nf90_inquire_DIMENSION(nc_file_id, dimIDs(1), len = tmp1)
      status = nf90_inquire_DIMENSION(nc_file_id, dimIDs(2), len = tmp2)
      allocate (var(tmp1,tmp2))

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var, start=start_var, count=(/tmp1,tmp2/) )
      if ((status /= nf90_noerr)) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

      !extract and save classifier names to the final array
      do i = 1, tmp2
        !if ((var(i,1) .ge. 'a' .and. var(i,1) .le. 'z') &
        !.or.(var(i,1) .ge. 'A' .and. var(i,1) .le. 'Z') &
        !.or.(var(i,1) .ge. '0' .and. var(i,1) .le. '9') ) then
           var_output(i) = trim(var(i,1))
        !endif
      ENDdo

      if (allocated(var)) deallocate (var)


   END SUBROUTINE READ_NETCDF_2D_CHAR

  ! ----------------------------------------------------------
  ! Read in 2D arrays
  ! ----------------------------------------------------------
  SUBROUTINE READ_NETCDF_2D_REAL(nc_file_id, var_start, var_stride, &
                                   var_dim, var_name, var_output)
      implicit none
      INTEGER, INTENT(in) :: nc_file_id
      INTEGER, DIMENSION(:), INTENT(in) :: var_start
      INTEGER, DIMENSION(:), INTENT(in) :: var_stride
      INTEGER, DIMENSION(:), INTENT(in) :: var_dim

      CHARACTER(len=*), INTENT(in) :: var_name
      REAL, INTENT(out), DIMENSION(:,:) :: var_output

      INTEGER :: nc_var_id
      INTEGER :: status

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


  END SUBROUTINE READ_NETCDF_2D_REAL

  ! ----------------------------------------------------------
  ! Read in 2D arrays INTEGER
  ! ----------------------------------------------------------
  SUBROUTINE READ_NETCDF_2D_INT(nc_file_id, var_start, var_stride, &
                                   var_dim, var_name, var_output)
      implicit none
      INTEGER, INTENT(in) :: nc_file_id
      INTEGER, DIMENSION(:), INTENT(in) :: var_start
      INTEGER, DIMENSION(:), INTENT(in) :: var_stride
      INTEGER, DIMENSION(:), INTENT(in) :: var_dim

      CHARACTER(len=*), INTENT(in) :: var_name
      INTEGER, INTENT(out), DIMENSION(:,:) :: var_output

      INTEGER :: nc_var_id
      INTEGER :: status

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


  END SUBROUTINE READ_NETCDF_2D_INT

  ! ----------------------------------------------------------
  ! Read in 3D arrays
  ! ----------------------------------------------------------
  SUBROUTINE READ_NETCDF_3D_REAL(nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)

      implicit none
      INTEGER, INTENT(in) :: nc_file_id
      INTEGER, DIMENSION(:), INTENT(in) :: var_start
      INTEGER, DIMENSION(:), INTENT(in) :: var_stride
      INTEGER, DIMENSION(:), INTENT(in) :: var_dim

      CHARACTER(len=*), INTENT(in) :: var_name
      REAL, INTENT(out), DIMENSION(:,:,:) :: var_output

      INTEGER :: nc_var_id
      INTEGER :: status = 0

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, stride=var_stride, count=var_dim)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif


  END SUBROUTINE READ_NETCDF_3D_REAL

  ! ----------------------------------------------------------
  ! Read in 4D arrays
  ! ----------------------------------------------------------
  SUBROUTINE READ_NETCDF_4D_REAL(nc_file_id, var_start, var_stride, &
                                  var_dim, var_name, var_output)

      implicit none
      INTEGER, INTENT(in) :: nc_file_id
      INTEGER, DIMENSION(:), INTENT(in) :: var_start
      INTEGER, DIMENSION(:), INTENT(in) :: var_stride
      INTEGER, DIMENSION(:), INTENT(in) :: var_dim
      CHARACTER(len=*), INTENT(in) :: var_name
      REAL, INTENT(out), DIMENSION(:,:,:,:) :: var_output

      INTEGER :: nc_var_id
      INTEGER :: status = 0

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=var_start, stride=var_stride, count=var_dim)
      if ((status /= nf90_noerr)) then
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      endif

  END SUBROUTINE READ_NETCDF_4D_REAL



END MODULE NB_CLOUD_MASK_NETCDF_READ_MODULE

