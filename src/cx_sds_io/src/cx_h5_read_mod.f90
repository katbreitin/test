! $Id: cx_h5_read_mod.f90 3888 2020-06-24 15:14:06Z awalther $
! CREATED Andi Walther {date}
module cx_h5_read_mod


   use hdf5,only:H5F_ACC_RDONLY_F,H5F_OBJ_FILE_F,H5F_OBJ_DATASET_F &
    ,H5F_OBJ_GROUP_F,H5F_OBJ_DATATYPE_F,h5fget_obj_count_f,h5fopen_f &
    , H5Open_F, H5GOPEN_F, H5GGET_INFO_F, H5GN_MEMBERS_F,H5GGET_OBJ_INFO_IDX_F &
    , H5G_LINK_F,H5G_GROUP_F,H5G_DATASET_F ,H5G_TYPE_F

   USE h5fortran_types

   USE cx_sds_type_definitions_mod, only: &
      cx_sds_type &
      , cx_att_type &
      , cx_sds_data_type

   USE    ReadH5Dataset, only: &
       H5_dataset_dimensions &
     , H5readattribute &
     , H5readdataset

   use linked_list_m

   implicit none
   private
   public:: h5_get_file_sds
   public:: h5_get_finfo

   integer, parameter :: MAXNCDIM = 32
   integer, parameter :: MAXNCNAM = 128

   include 'cx_sds_constants.inc'
   integer,parameter :: DIM_0 = 0
   integer,parameter :: DIM_1 = 1
   integer,parameter :: DIM_2 = 2
   integer,parameter :: DCLASS_0 = 0
   integer,parameter :: DCLASS_1 = 1
   integer,parameter :: DCLASS_INT = 0
   integer,parameter :: DCLASS_REAL = 1
   integer,parameter :: DCLASS_2 = 2
   integer,parameter :: DCLASS_3 = 3

contains


    function group_info (root_id,name, obj_name, obj_type)
      integer :: group_info
      integer (HID_T) :: root_id
      character(len=*) :: name
      character(len=*) :: obj_name
      integer :: obj_type
      integer :: hdferr
      integer :: nmem
      integer :: i

      group_info = 0

       call H5GN_MEMBERS_F ( root_id,'/',nmem,hdferr)

       call H5GGET_OBJ_INFO_IDX_F(root_id, '/', i, &
                                 obj_name, obj_type, hdferr)


    end function group_info


   function h5_get_finfo(h5_file, nsds, sds_name, natt, att_name)
      ! - this tool is incomplte. It does not popllte sds_name
      ! - 8 Dec 2021 AW
      !
      integer :: h5_get_finfo
      character(len=*), intent(in) :: h5_file
      integer, intent(out) :: nsds
      integer, intent(out) :: natt
      character ( len = MAXNCNAM), intent(out), allocatable :: sds_name(:)
      character ( len = MAXNCNAM), intent(out), allocatable :: att_name(:)
      integer (HID_T) :: file_id, root_id
      integer :: obj_type, obj_type1
      integer( SIZE_T) :: obj_count
      integer :: hdferr
      INTEGER :: ErrorFlag
      integer :: s_type,nlinks,max_corder
      integer :: nmem,nmem1
      character(len = 100) :: obj_name, obj_name1
      integer :: i, ii
      character(len =100) :: dir
      type(LinkedList) :: list_sds
      type(LinkedListNode), pointer :: node_ptr
      character(len=200) ::sds_name1

      obj_type = 1
      ErrorFlag=0

      CALL h5open_f(ErrorFlag)
      ! print*,obj_type, H5F_OBJ_DATASET_F

      h5_get_finfo = 1

      nsds = 0
      natt = 0

      allocate(att_name(0))

      call h5fopen_f(h5_file,H5F_ACC_RDONLY_F, file_id,hdferr)
      CALL h5gopen_f(file_id, "/", root_id, ErrorFlag)

      call H5GGET_INFO_F ( root_id,s_type,nlinks,max_corder,hdferr)

      dir = '/'

      call search_dataset(root_id, dir, nsds, list_sds)


      ii = list_sds%length()

      allocate(sds_name(ii))
      do i=1,list_sds%length()
        node_ptr => list_sds%atindex(i)
        select type(p => node_ptr%value)
        type is(character(*))
          sds_name(i) = p
        class default

        end select
      end do

      do i = 1,list_sds%length()
        print*,i,trim(sds_name(i))
      end do
   end function h5_get_finfo


   recursive subroutine search_dataset(root_id, dir, nsds ,list_sds)
     integer (HID_T), intent(in) :: root_id
     character(len=*), intent(inout) :: dir
     integer, intent(inout) :: nsds
     type(LinkedList), intent(inout) :: list_sds
     integer :: i
     integer :: nmem
     integer :: hdferr
     character(len = 100) :: obj_name
     integer :: obj_type
     character(len =100) :: dir_loc
     class(*), pointer :: general

     call H5GN_MEMBERS_F ( root_id,'/',nmem,hdferr)
     do i = 0 , nmem - 1
       call H5GGET_OBJ_INFO_IDX_F(root_id, trim(dir), i, &
                                obj_name, obj_type, hdferr)

        if (obj_type .EQ. H5G_DATASET_F) THEN
            nsds = nsds + 1
            allocate(general,source=trim(dir)//trim(obj_name))
            call list_sds % append(general)
            cycle
        end if
        if (obj_type .EQ. H5G_GROUP_F) THEN
            dir_loc = '/'//trim(obj_name)
            CALL search_dataset(root_id, dir_loc, nsds, list_sds)
        endif
     end do



   end subroutine search_dataset


  !
  !
  !
   function h5_get_file_sds(h5_file, nsds, sdata,  sds_name, start_inp, stride_inp &
      , count_inp)

      integer :: h5_get_file_sds

      character (len=*), intent( in) :: h5_file
      integer, intent(out)  :: nsds
      type(cx_sds_type), intent(out),  allocatable, target :: sdata(:)
      character (len=*), intent( in), optional :: sds_name(:)
      integer, optional, intent(in) :: start_inp(:)
      integer, optional, intent(in) :: stride_inp(:)
      integer, optional, intent(in) :: count_inp(:)
      integer , pointer :: dims(:)
      integer :: ndims
      integer, allocatable  :: dataset_2d(:,:)
      integer(kind=4), allocatable :: dataset_2d_i(:,:)
      real, pointer  :: dataset_1d(:)
      real :: att1
      integer :: att_i4
      real:: value

      integer ::dclass
      integer :: offs(2),counts(2)
      ! ------------------------------------------------------
      h5_get_file_sds = -1
      call H5_DATASET_DIMENSIONS( h5_file,sds_name(1),dims,dclass)



      ndims = size(dims)

      ! check which data type is this data set

      allocate (sdata(1))



      select case (ndims)


      case (DIM_0)
        call H5ReadDataset ( h5_file, sds_name(1), value)

         allocate ( sdata(1) % data % r4values(1))
         sdata(1) % data % r4values(1) = value
         sdata(1) % data % nval = 1
         sdata(1) % data % type = DFNT_FLOAT32
      case(DIM_2)


        select case (dclass)

        case(DCLASS_REAL)
          if (  present (start_inp)) then
            offs = (/0,start_inp(2)/)
            counts = (/dims(1),count_inp(2)/)

            call H5ReadDataset ( h5_file, sds_name(1),offs,counts, dataset_2d )
                    dims(2) = count_inp(2)


          else
            call H5ReadDataset ( h5_file, sds_name(1), dataset_2d )
          end if

          allocate ( sdata(1) % data % r4values(dims(1)*dims(2)))
          sdata(1) % data % r4values = reshape(dataset_2d, (/dims(1)*dims(2)/) )
          sdata(1) % data % nval = size(dataset_2d)
          sdata(1) % data % type = DFNT_FLOAT32
           if (allocated(dataset_2d)) deallocate(dataset_2d)

          case(DCLASS_INT)
            if (  present (start_inp)) then
              offs = (/0,start_inp(2)/)
              counts = (/dims(1),count_inp(2)/)

              call H5ReadDataset ( h5_file, sds_name(1),offs,counts, dataset_2d_i )
              dims(2) = count_inp(2)
            else
              call H5ReadDataset ( h5_file, sds_name(1) , dataset_2d_i )
            end if

            allocate ( sdata(1) % data % i4values(dims(1)*dims(2)))
            sdata(1) % data % i4values = reshape(dataset_2d_i, (/dims(1)*dims(2)/) )
            sdata(1) % data % nval = size(dataset_2d_i)
            sdata(1) % data % type = 24
            if ( allocated(dataset_2d_i)) deallocate(dataset_2d_i)

          end select

          sdata(1) % data % dimsize(1) = dims(1)
          sdata(1) % data % dimsize(2) = dims(2)
          sdata(1) % data % rank = 2

        case(DIM_1)

          call H5ReadDataset ( h5_file, sds_name(1), dataset_1d )

          allocate ( sdata(1) % data % r4values(dims(1)))
          sdata(1) % data % r4values = dataset_1d
          sdata(1) % data % nval = size(dataset_1d)
          deallocate(dataset_1d)
          sdata(1) % data % dimsize(1) = dims(1)
          sdata(1) % data % type = 24

        end select

        nsds = 1
        h5_get_file_sds = 1

        allocate  (sdata(1) % attr(4))

        sdata(1) % nattr = 4
        sdata(1) % attr(1) % name = 'add_offset'
        sdata(1) % attr(2) % name = 'scale_factor'
        sdata(1) % attr(3) % name = '_FillValue'
        sdata(1) % attr(4) % name = 'SCALED'

        sdata(1) % attr(1) % data % type = DFNT_FLOAT32
        sdata(1) % attr(2) % data % type = DFNT_FLOAT32
        sdata(1) % attr(3) % data % type = DFNT_INT16
        sdata(1) % attr(4) % data % type = DFNT_INT8

        allocate ( sdata(1) % attr(1) % data % r4values(1))
        allocate ( sdata(1) % attr(2) % data % r4values(1))
        allocate ( sdata(1) % attr(3) % data % i4values(1))
        allocate ( sdata(1) % attr(4) % data % i1values(1))

        !- init
        sdata(1) % attr(1) % data % r4values = 0.00
        sdata(1) % attr(2) % data % r4values = 1.00
        sdata(1) % attr(3) % data % i4values = 0
        sdata(1) % attr(4) % data % i1values = 0


        if ( dclass .eq. DCLASS_INT ) then
          call H5ReadAttribute( h5_file,trim(sds_name(1))//'/add_offset',att1)
          sdata(1) % attr(1) % data % r4values = att1

          call H5ReadAttribute( h5_file,trim(sds_name(1))//'/scale_factor',att1)
          sdata(1) % attr(2) % data % r4values = att1

          call H5ReadAttribute( h5_file,trim(sds_name(1))//'/_FillValue',att_i4)

          sdata(1) % attr(3) % data % i4values = att_i4

          sdata(1) % attr(4) % data % i1values = 1
        end if


      end function h5_get_file_sds



end module cx_h5_read_mod
