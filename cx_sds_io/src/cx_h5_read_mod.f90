! $Id: cx_h5_read_mod.f90 3888 2020-06-24 15:14:06Z awalther $
! CREATED Andi Walther {date} 
module cx_h5_read_mod


   use hdf5,only:H5F_ACC_RDONLY_F,H5F_OBJ_FILE_F,H5F_OBJ_DATASET_F &
    ,H5F_OBJ_GROUP_F,H5F_OBJ_DATATYPE_F,h5fget_obj_count_f,h5fopen_f
   USE h5fortran_types
   USE cx_sds_type_definitions_mod, only: &
      cx_sds_type &
      , cx_att_type &
      , cx_sds_data_type

   USE    ReadH5Dataset

   implicit none
   integer, parameter :: MAXNCDIM = 32
   integer, parameter :: MAXNCNAM = 128
contains

   function h5_get_finfo(h5_file, nsds, sds_name, natt, att_name)

      integer :: h5_get_finfo
      character(len=*), intent(in) :: h5_file
      integer, intent(out) :: nsds
      integer, intent(out) :: natt
      character ( len = MAXNCNAM), intent(out), allocatable :: sds_name(:)
      character ( len = MAXNCNAM), intent(out), allocatable :: att_name(:)
      integer (HID_T) :: file_id
      integer :: obj_type
      integer( SIZE_T) :: obj_count
      integer :: hdferr

      obj_type = 1
      print*,obj_type, H5F_OBJ_DATASET_F
      
      h5_get_finfo = 1
      print*,' hdf5 finfo not yet installed stopping'
      print*,' will be installed by end of April 2018'
      nsds = 0
      natt = 0
      allocate(sds_name(0))
       allocate(att_name(0))
      
      print*,'AA1'
     
      print*,H5F_OBJ_FILE_F,H5F_OBJ_DATASET_F,H5F_OBJ_GROUP_F,H5F_OBJ_DATATYPE_F
       
       
       
      call h5fopen_f(h5_file,H5F_ACC_RDONLY_F, file_id,hdferr)
       print*,'AA2',hdferr
       do obj_type=1,50
      call h5fget_obj_count_f(file_id, obj_type, obj_count, hdferr)
       print*,'AA3',hdferr,obj_type,obj_count
     end do
      obj_type = H5F_OBJ_GROUP_F
      call h5fget_obj_count_f(file_id, obj_type, obj_count, hdferr)
       print*,'AA4',hdferr

      obj_type = H5F_OBJ_DATASET_F
      call h5fget_obj_count_f(file_id, obj_type, obj_count, hdferr)
       print*,'AA5',hdferr

   end function h5_get_finfo



   function h5_get_file_sds(h5_file, nsds, sdata, nsdsn, sds_name)

      integer :: h5_get_file_sds

      character (len=*), intent( in) :: h5_file
      integer, intent(out)  :: nsds
      type(cx_sds_type), intent(out),  allocatable, target :: sdata(:)
      integer, optional, intent( in) :: nsdsn
      character (len=*), intent( in), optional :: sds_name(:)
      integer , pointer :: dims(:)
      integer :: ndims
      real, pointer  :: dataset_2d(:,:)
      real, pointer  :: dataset_1d(:)

print*,'statrttttt'
print*,h5_file
print*,sds_name(1)

      call H5_DATASET_DIMENSIONS( h5_file,sds_name(1),dims)

      ndims = size(dims)
   print*,'ndims===> ',ndims

      select case (ndims)
         
         case(2)
        
            call H5ReadDataset ( h5_file, sds_name(1), dataset_2d )
        
            allocate (sdata(1))
            allocate ( sdata(1) % data % r4values_2d(dims(1),dims(2)))
        
            sdata(1) % data % r4values_2d = dataset_2d
        
            sdata(1) % data % nval = size(dataset_2d)

            !deallocate(dataset_2d)
            sdata(1) % data % dimsize(1) = dims(1)
            sdata(1) % data % dimsize(2) = dims(2)
        
            sdata(1) % data % type = 5
            sdata(1) % data % rank = 2
        


         case(1)

            call H5ReadDataset ( h5_file, sds_name(1), dataset_2d )
            allocate (sdata(1))
            allocate ( sdata(1) % data % r4values(dims(1)))
            sdata(1) % data % r4values = dataset_1d
            sdata(1) % data % nval = size(dataset_1d)
            deallocate(dataset_1d)
            sdata(1) % data % dimsize(1) = dims(1)



      end select

      nsds = 1
      h5_get_file_sds = 1

   end function h5_get_file_sds



end module cx_h5_read_mod
