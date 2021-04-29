! $Id: cx_ncdf_read_mod.f90 3882 2020-06-23 10:46:24Z awalther $
!------------------------------------------------------------------------------
! SSEC/UW, CLAVR-x Software Tools Science Data Read Routines
!------------------------------------------------------------------------------
!
! MODULE: cx_ncdf_read_mod
!
!> @author
!> CX_NCDF_READ_MOD Andi Walther CIMSS
!
! DESCRIPTION:
!> This module deals with all NETCDF3 read routines
!
! REVISION HISTORY:
! 12 Mar 2018 - Initial Version
! 20 May 2018 - Added subsection capability AW
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
Module cx_ncdf_read_mod

   use netcdf

   use cx_sds_type_definitions_mod, only: &
      cx_sds_type &
      , cx_att_type &
      , cx_sds_data_type

   implicit none
   integer, parameter :: MAXNCDIM = 32
   integer, parameter :: MAXNCNAM = 128
   integer, parameter :: DECLBRTD = 1
   integer, parameter :: UNCLBRTD = 0
   integer, parameter :: CALIBRTD = 2

contains

   ! ---------------------------------------------------------------------------
   !
   ! ---------------------------------------------------------------------------
   function ncdf_get_finfo(file, nsds, sds_name, natt, att_name)

      integer :: ncdf_get_finfo
      character(len=*), intent(in) :: file
      integer, intent(out) :: nsds
      integer, intent(out) :: natt
      character ( len = MAXNCNAM), intent(out), allocatable :: sds_name(:)
      character ( len = MAXNCNAM), intent(out), allocatable :: att_name(:)
      integer :: ncid
      integer :: iobj
      integer :: status
      integer :: nDimensions
      integer :: unlimitedDimId
      integer :: formatNum

      ncdf_get_finfo = 0
      status = nf90_open(file, nf90_nowrite, ncid)
      status = nf90_inquire(ncid, nDimensions, nsds, nAtt, &
         unlimitedDimId, formatNum)

      allocate ( sds_name (nsds))
      allocate ( att_name (natt))
      do iobj = 1, natt
         status= nf90_inq_attname(ncid,NF90_GLOBAL,iobj,att_name(iobj))

      end do
   
      do iobj = 1, nsds
         status= nf90_inquire_variable(ncid,iobj,sds_name(iobj))

      end do

      status = nf90_close ( ncid)

   end function ncdf_get_finfo
   
   ! ---------------------------------------------------------------------------
   !
   ! ---------------------------------------------------------------------------      
   function ncdf_get_varinfo(file,  var, natt, att_name, ndim, dim)
      integer :: ncdf_get_varinfo
      character(len=*), intent(in) :: file
      character(len=*), intent(in) :: var
      integer :: natt
      character ( len = MAXNCNAM), intent(out), allocatable :: att_name(:)
      integer :: ndim
      integer :: dim(10)
      integer :: status
      integer :: dimid(10)
      integer :: i
      integer :: ncid, sdsid
      
      ncdf_get_varinfo = 0
      status = nf90_open(file, nf90_nowrite, ncid)
      Status = nf90_inq_varid(ncid, trim(var), sdsid)
      
      status= nf90_inquire_variable(ncid,sdsid,ndims=ndim, dimids=dimid)
      
      do i=1, ndim
      status= nf90_inquire_dimension(ncid,dimid(i),len = dim(i))
      end do
   end function ncdf_get_varinfo
   
   !-------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------
   function ncdf_get_file_att(file, natt, attrs, nattn, att_name)

      integer :: ncdf_get_file_att
      character(len=*), intent(in) :: file
      integer, intent(out) :: natt
      type(cx_att_type), intent(out), allocatable :: attrs(:)
      integer, intent(in), optional:: nattn
      character(len=*), intent(in), optional :: att_name(*)
      integer :: status
      integer :: ncid
      integer :: nsds
      integer :: nDimensions,  unlimitedDimId, &
         formatNum

      ncdf_get_file_att = -1

      status = nf90_open(file, mode = nf90_nowrite,ncid=ncid)
      if (status /= NF90_NOERR) return
      status = nf90_inquire(ncid, nDimensions, nsds, nAtt, &
         unlimitedDimId, formatNum)

      if (natt > 0) then
         if (ncdf_get_obj_att(ncid, NF90_GLOBAL, natt, attrs, nattn, att_name) < 0) goto 99999
      end if

      ncdf_get_file_att = 0

99999 continue

      if (ncdf_error(nf90_close ( ncid))) ncdf_get_file_att = -1

   end function ncdf_get_file_att

   !-------------------------------------------------------------------------------
   !  Reads attribute
   !    if nattn and att_name are not given, all attributes are read and 
   !    returned in attrs structure
   ! --------------------------
   function ncdf_get_obj_att(ncid, obj_id, natt, attrs, nattn, att_name)
      implicit none
      integer :: ncdf_get_obj_att
      integer, intent(in) :: ncid
      integer, intent(in) :: obj_id
      integer, intent(out) :: natt
      type(cx_att_type), intent(out), dimension(:), allocatable :: attrs
      integer, intent( in), optional :: nattn
      character(len=*), intent( in), optional :: att_name(*)
      integer :: iatt
      integer :: ierr
      integer, dimension(:), allocatable :: attind
      integer :: status
      character(len=256) :: att_name_temp
      character (len=100) :: temp_name
      
      
      ! - Executable
      ncdf_get_obj_att = -1

      !NF90_INQ_VARNATTS
      if (std_error((present(nattn).or.present(att_name)) &
         .and.(.not.(present(nattn).and.present(att_name))),  &
         "Optional arguments must be both defined or undefined")) return
      
      natt = 0
      status = NF90_INQUIRE_VARIABLE (ncid,obj_id, nAtts= natt)
      
      if (present(nattn)) natt = nattn
        
      if (natt <= 0) return

      allocate (attind(1:natt), stat=ierr)
      if (std_error(ierr /= 0, "Dynamic memory allocation error")) return

      if (.not.present(nattn)) then
         do iatt = 1, natt
            attind(iatt) = iatt - 1
         enddo
      else
         natt = nattn
         do iatt = 1, natt
            att_name_temp = att_name( iatt)
            status = nf90_inq_attname ( ncid,obj_id ,attind(iatt),  att_name_temp )
              

            if (ncdf_error(attind(iatt))) goto 99999
         end do
      end if

      if (allocated(attrs)) then
         deallocate (attrs, stat=ierr)
         if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
      endif
      allocate (attrs(1:natt), stat=ierr)
      if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999

      do iatt = 1, natt
         status = NF90_INQ_ATTNAME( ncid, obj_id, iatt,att_name_temp )
         attrs(iatt)%name = trim(att_name_temp)
         status = NF90_INQUIRE_ATTRIBUTE( ncid, obj_id, att_name_temp, attrs(iatt)%data%type_ncdf, attrs(iatt)%data%dimsize(1))
          
          call attrs(iatt) % data % type_ncdf_to_hdf() 
          
         attrs(iatt)%data%rank = 1
         attrs(iatt)%data%datasize = ncdf_typesize(attrs(iatt)%data%type)
         attrs(iatt)%data%size = attrs(iatt)%data%dimsize(1)*attrs(iatt)%data%datasize
         attrs(iatt)%data%nval = attrs(iatt)%data%dimsize(1)
         
         select case (attrs(iatt)%data%type_ncdf)

            case (NF90_CHAR)
               allocate (attrs(iatt)%data%c1values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
               status = NF90_GET_ATT ( ncid,obj_id,trim(attrs(iatt)%name),attrs(iatt)%data%char_nc_values )
               if ( status /= 0) print*,'char: ',status

            case (NF90_SHORT)
               allocate (attrs(iatt)%data%i1values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
               status = NF90_GET_ATT ( ncid,obj_id,trim(attrs(iatt)%name),attrs(iatt)%data%i1values )
                if ( status /= 0) print*,'short: ',status
            !      if (ncdf_error(sfrcatt(obj_id, iatt - 1, attrs(iatt)%data%i1values))) goto 99999

            case (NF90_INT)
               allocate (attrs(iatt)%data%i2values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
               status = NF90_GET_ATT ( ncid,obj_id,trim(attrs(iatt)%name),attrs(iatt)%data%i2values )
               if ( status /= 0)  print*,'int: ',status, attrs(iatt)%data%type_ncdf
               
             case (NF90_USHORT)
               allocate (attrs(iatt)%data%i4values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
               status = NF90_GET_ATT ( ncid,obj_id,trim(attrs(iatt)%name),attrs(iatt)%data%i4values )
              
               if ( status /= 0)  print*,'ushort: ',status, attrs(iatt)%data%type_ncdf

            case (NF90_FLOAT )
               allocate (attrs(iatt)%data%r4values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
               status = NF90_GET_ATT ( ncid,obj_id,trim(attrs(iatt)%name),attrs(iatt)%data%r4values )
               if ( status /= 0)  print*,'float: ',status
            !        if (ncdf_error(sfrcatt(obj_id, iatt - 1, attrs(iatt)%data%r4values))) goto 99999

            case (NF90_DOUBLE)
               !allocate (attrs(iatt)%data%r8values(1:attrs(iatt)%data%dimsize(1)), stat=ierr)
               allocate (attrs(iatt)%data%r8values(1), stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
              
              
               status = NF90_GET_ATT ( ncid,obj_id,trim(attrs(iatt)%name),attrs(iatt)%data%r8values )
              if ( status /= 0)  print*,'double: ',status

            case default
               print*,NF90_INT, NF90_UINT, NF90_UINT64, NF90_USHORT, NF90_STRING
               print *, "Unimplemented data type for attribute: ", att_name(iatt),attrs(iatt)%data%type_ncdf; goto 99999

         end select

      end do

      ncdf_get_obj_att = 0

99999 continue



      if (allocated(attind)) then
         deallocate (attind, stat=ierr)
         if (std_error(ierr /= 0, "Dynamic memory desallocation error"))  ncdf_get_obj_att = -1
      endif
      
      

   end function ncdf_get_obj_att

   ! ----------------------------------------------------------------------------------------
   !
   ! ----------------------------------------------------------------------------------------
  function ncdf_get_file_sds(ncdf_file, nsds, sdata, nsdsn, sds_name, start_inp, stride_inp &
      , count_inp, dclb, attr, cal_sub, outtype)

    integer :: ncdf_get_file_sds

    character (len=*), intent( in) :: ncdf_file
    integer, intent(out)  :: nsds
    type(cx_sds_type), intent(out),  allocatable, target :: sdata(:)
    integer, optional, intent( in) :: nsdsn
    character (len=*), intent( in), optional :: sds_name(:)
    integer, optional, intent(in) :: start_inp(:)
    integer, optional, intent(in) :: stride_inp(:)
    integer, optional, intent(in) :: count_inp(:)

    logical, optional :: dclb
    logical, optional :: attr
     
    integer, optional, intent( in) :: outtype
    optional :: cal_sub
    external cal_sub

    integer :: sd_id
    integer :: sds_id
    integer :: start(MAXNCDIM)
    integer :: stride(MAXNCDIM)
    integer :: count(MAXNCDIM)
    integer :: ierr
    integer :: isds
    integer :: idim
    integer :: rtype
    integer :: natt
    integer :: status
    integer :: nDimensions
    integer :: unlimitedDimId
    integer :: formatNum

    logical:: att_sw
    integer, allocatable :: sdsind(:)
    integer*1, allocatable :: bdata(:)
    real(kind=8), allocatable :: r8data(:)
    type(cx_sds_type), pointer :: ps
    type(cx_sds_data_type), pointer :: pd
    integer :: i

    integer :: count_clb
    logical :: calibrated = .false.
    integer :: dim_id(36)
      
      ! -------------
      !   
      ncdf_get_file_sds = -1

      if (std_error((present(nsdsn).or.present(sds_name)).and.(.not.(present(nsdsn).and.present(sds_name))), &
         'Optional arguments must be both defined or undefined')) return

      status = nf90_open(ncdf_file, mode = nf90_nowrite,ncid=sd_id)

      status = nf90_inquire(sd_id, nDimensions, nsds, nAtt, &
         unlimitedDimId, formatNum)
     
      if (present(nsdsn)) nsds = nsdsn
      att_sw = .true.

      if (present(attr)) att_sw = attr

      if (nsds <= 0) goto 99999

      allocate (sdsind(1:nsds), stat=ierr)
      if (std_error(ierr /= 0, "Dynamic memory allocation error l294")) goto 99999

      if (.not.present(nsdsn)) then
         do isds = 1, nsds
            sdsind(isds) = isds - 1
         end do
      else
         do isds = 1, nsds

            status = NF90_INQ_VARID(sd_id, sds_name(isds) , sdsind(isds))

            if (ncdf_error(status)) goto 99999
         end do
      end if

      
      allocate (sdata(1:nsds), stat=ierr)
      if (std_error(ierr /= 0, "Dynamic memory allocation error cx-hdf_read-mod")) goto 99999

      do isds = 1, nsds

         sds_id = sdsind(isds)
         if (ncdf_error(sds_id)) goto 99999

         ps => sdata(isds)
         pd => sdata(isds)%data
         status = NF90_INQUIRE_VARIABLE (sd_id,sds_id,ps%name , pd%type_ncdf , pd%rank , dim_id , ps%nattr)
        
        call pd % type_ncdf_to_hdf()
        
       
         do i = 1, pd % rank
            status = nf90_inquire_dimension(sd_id, dim_id(i) , len = pd%dimsize(i))
         end do

         !  if (ncdf_error(sfginfo(sds_id, ps%name, pd%rank, pd%dimsize, pd%type, ps%nattr))) &
         !    & goto 99999

         pd % calbrtd = DECLBRTD
        
         !   if (sfgcal(sds_id, pd%calibr(1), pd%calibr(2), pd%calibr(3), pd%calibr(4), pd%utype) < 0) then
         if (.not. calibrated ) then
            count_clb = 0

            if (ps%nattr > 0) then 
               status = ncdf_get_obj_att(sd_id, sds_id, ps%nattr, ps%attr)
            end if

            pd%utype = pd%type_ncdf
            if ( count_clb .EQ. 2) then
               pd%calbrtd = DECLBRTD

            else
               pd%calbrtd = UNCLBRTD
               pd%calibr = 0.; pd%calibr(1) = 1.
            end if

         else if (present(dclb)) then
            if (.not.dclb) pd%calbrtd = CALIBRTD
         end if

         if (pd%calbrtd == DECLBRTD) then
            rtype = pd%utype;
         else
            rtype = pd%type_ncdf;
         endif
                  
                
         pd % datasize = ncdf_typesize(rtype)
         pd % nval = 1
         
         stride(:) = 1
         start (:) = 1
                   
         if ( present ( stride_inp) ) then           
            stride(1:pd % rank)  = stride_inp(1:pd % rank)
         end if
         if ( present ( start_inp) ) then           
            start(1:pd % rank)  = start_inp(1:pd % rank)
         end if
         
         if ( present ( count_inp) ) then           
            ! - if count_inp set, the dimension is the count
            ! - added a check if count is too high.. if so decrease. AW 2018/06/12
            do idim = 1, pd%rank
                
                pd%dimsize(idim) = min(count_inp(idim), pd%dimsize(idim)- start(idim) + 1)
                pd%nval = pd%nval * pd%dimsize(idim)
                
            end do  
           
         else 
             
            do idim = 1, pd%rank
               pd%dimsize(idim) = (pd%dimsize(idim)- start(idim) + 1) / stride (idim)
              
               pd%nval = pd%nval * pd%dimsize(idim)
            end do  
            
         end if

         if (present(outtype)) then
            pd%nval = pd%nval*ncdf_typesize(rtype)/ncdf_typesize(outtype)
            pd%datasize = ncdf_typesize(outtype)
            pd%type_ncdf = outtype
            rtype = outtype
         end if

         pd % size = pd%nval*pd%datasize
 
         if (pd%calbrtd .EQ. DECLBRTD) then
            if (allocated(bdata)) then
               deallocate(bdata, stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
            endif
            allocate(bdata(pd%size), stat=ierr)
            if (std_error(ierr /= 0, "Dynamic memory allocation error l386")) goto 99999

            status = NF90_GET_VAR ( sd_id, sds_id, bdata, start= start, stride = stride)

            !  if (ncdf_error(sfrdata(sds_id, start, stride, pd%dimsize, bdata))) goto 99999

            if (allocated(r8data)) then
               deallocate(r8data, stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
            endif
            allocate(r8data(pd%nval), stat=ierr)
            if (std_error(ierr /= 0, "Dynamic memory allocation error l397")) goto 99999

 
            select case (pd%type_ncdf)
               case (NF90_CHAR)
                  print *, "Character data cannot be decalibrated"; goto 99999
               case (NF90_BYTE)
                  r8data = transfer(bdata, pd%i1values)
               case (NF90_SHORT)
                  r8data = transfer(bdata, pd%i2values)
               case (NF90_INT, NF90_USHORT)
                  r8data = transfer(bdata, pd%i4values)
               case (NF90_FLOAT)
                  r8data = transfer(bdata, pd%r4values)
               case (NF90_DOUBLE)
                  r8data = transfer(bdata, pd%r8values)
               case default
                  print *, "Unimplemented data type 2: ", pd%type_ncdf; goto 99999
            end select

            if (present(cal_sub)) then
               call cal_sub(pd%nval, r8data, pd%calibr(1), pd%calibr(3))
            else

               call default_dclb( r8data, pd%calibr(1), pd%calibr(3))

            end if

            if (allocated(bdata)) then
               deallocate(bdata, stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory desallocation error l 428")) goto 99999
            endif

         end if
 
         !
         !
         !
      
         select case (rtype)

            case (NF90_CHAR)
               allocate (pd%c1values(1:pd%nval), stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory allocation error l 438")) goto 99999
               if (pd%calbrtd == DECLBRTD) then
                  pd%i1values = transfer(r8data,int(1,1))
               else
                  status = NF90_GET_VAR ( sd_id, sds_id, pd%c1values, start= start, stride = stride)
                 !  if (ncdf_error(sfrdata(sds_id, start, stride, pd%dimsize, pd%c1values))) goto 99999
               endif

            case (NF90_SHORT)
                  allocate (pd%i2values(1:pd%nval), stat=ierr)
                 select case(pd %rank)
                  case(1)
                    allocate (pd%i2values(1:pd%nval), stat=ierr)
                  case(2)
                     allocate (pd%i2values_2d(pd%dimsize(1),pd%dimsize(2)), stat=ierr)
                  case(3)
                     allocate (pd%i2values_3d(pd%dimsize(1),pd%dimsize(2),pd%dimsize(3)), stat=ierr)
                  case(4)
                     allocate (pd%i2values_4d(pd%dimsize(1),pd%dimsize(2),pd%dimsize(3),pd%dimsize(4)), stat=ierr)
                  case(5)
                     allocate (pd%i2values_5d(pd%dimsize(1),pd%dimsize(2),pd%dimsize(3) &
                        ,pd%dimsize(4),pd%dimsize(5)), stat=ierr)
               end select
               
               
               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
               if (pd%calbrtd == DECLBRTD) then
                  pd%i1values = transfer(r8data,int(1,1))
               else
                 
                  
                  select case(pd %rank)
                    case(1)
                      status = NF90_GET_VAR ( sd_id, sds_id, pd%i1values, start= start, stride = stride)
                    case(2)
                      status = NF90_GET_VAR ( sd_id, sds_id, pd%i2values_2d, start= start(1:2) &
                           , count = pd%dimsize(1:2), stride = stride(1:2))
                      
                      pd%i2values = reshape( pd%i2values_2d,(/pd%nval /))
                        
                        
                  end select
                 
               endif

            case (NF90_INT)
               allocate (pd%i4values(1:pd%nval), stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
               if (pd%calbrtd == DECLBRTD) then
                  pd%i2values = transfer(r8data,int(1,2))
               else
                  status = NF90_GET_VAR ( sd_id, sds_id, pd%i4values, start= start, stride = stride)

               endif
          
            case ( NF90_USHORT)
              ! - USHORT is int2 but unsigned
              ! - instead of transformation we use here int4 what covers full range
              
               allocate (pd%i4values(1:pd%nval), stat=ierr)

 
               select case(pd %rank)
                  case(1)

                  case(2)
                     allocate (pd%i4values_2d(pd%dimsize(1),pd%dimsize(2)), stat=ierr)
                  case(3)
                     allocate (pd%i4values_3d(pd%dimsize(1),pd%dimsize(2),pd%dimsize(3)), stat=ierr)
                  case(4)
                     allocate (pd%i4values_4d(pd%dimsize(1),pd%dimsize(2),pd%dimsize(3),pd%dimsize(4)), stat=ierr)
                  case(5)
                     allocate (pd%i4values_5d(pd%dimsize(1),pd%dimsize(2),pd%dimsize(3) &
                        ,pd%dimsize(4),pd%dimsize(5)), stat=ierr)
               end select

               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
               if (pd%calbrtd == DECLBRTD) then
                  pd%i4values = transfer(r8data,1.)
               else
                  select case(pd %rank)
                     case(1)
                        status = NF90_GET_VAR ( sd_id, sds_id, pd%i2values)
                     case(2)
                         
                        status = NF90_GET_VAR ( sd_id, sds_id, pd%i4values_2d, start= start(1:2) &
                           , count = pd%dimsize(1:2), stride = stride(1:2))
                      
                        pd%i4values = reshape( pd%i4values_2d,(/pd%nval /))
                       
                     case(3)
                        status = NF90_GET_VAR ( sd_id, sds_id, pd%i4values_3d, start= start(1:3), stride = stride(1:3))
                        pd%i4values = reshape( pd%i4values_3d,(/pd%nval /))
                     case(4)
                        status = NF90_GET_VAR ( sd_id, sds_id, pd%i4values_4d, start= start(1:4), stride = stride(1:4))
                        pd%i4values = reshape( pd%i4values_4d,(/pd%nval /))
                     case(5)
                        status = NF90_GET_VAR ( sd_id, sds_id, pd%i4values_5d, start= start(1:5), stride = stride(1:5))
                        pd%i4values = reshape( pd%i4values_5d,(/pd%nval /))
                  end select

               endif
               

           case (NF90_FLOAT)
              
               allocate (pd%r4values(1:pd%nval), stat=ierr)

 
               select case(pd %rank)
                  case(1)

                  case(2)
                     allocate (pd%r4values_2d(pd%dimsize(1),pd%dimsize(2)), stat=ierr)
                  case(3)
                     allocate (pd%r4values_3d(pd%dimsize(1),pd%dimsize(2),pd%dimsize(3)), stat=ierr)
                  case(4)
                     allocate (pd%r4values_4d(pd%dimsize(1),pd%dimsize(2),pd%dimsize(3),pd%dimsize(4)), stat=ierr)
                  case(5)
                     allocate (pd%r4values_5d(pd%dimsize(1),pd%dimsize(2),pd%dimsize(3) &
                        ,pd%dimsize(4),pd%dimsize(5)), stat=ierr)
               end select

               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
               if (pd%calbrtd == DECLBRTD) then
                  pd%r4values = transfer(r8data,1.)
               else
                  select case(pd %rank)
                     case(1)
                        status = NF90_GET_VAR ( sd_id, sds_id, pd%r4values)
                     case(2)
                       
                        status = NF90_GET_VAR ( sd_id, sds_id, pd%r4values_2d, start= start(1:2) &
                           , count = pd%dimsize(1:2), stride = stride(1:2))
                       
                        pd%r4values = reshape( pd%r4values_2d,(/pd%nval /))
                       
                     case(3)
                        status = NF90_GET_VAR ( sd_id, sds_id, pd%r4values_3d, start= start(1:3), stride = stride(1:3))
                        pd%r4values = reshape( pd%r4values_3d,(/pd%nval /))
                     case(4)
                        status = NF90_GET_VAR ( sd_id, sds_id, pd%r4values_4d, start= start(1:4), stride = stride(1:4))
                        pd%r4values = reshape( pd%r4values_4d,(/pd%nval /))
                     case(5)
                        status = NF90_GET_VAR ( sd_id, sds_id, pd%r4values_5d, start= start(1:5), stride = stride(1:5))
                        pd%r4values = reshape( pd%r4values_5d,(/pd%nval /))
                  end select

               endif

            case (NF90_DOUBLE)
               allocate (pd%r8values(1:pd%nval), stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory allocation error")) goto 99999
               if (pd%calbrtd == DECLBRTD) then
                  pd%r8values = r8data
               else
                  status = NF90_GET_VAR ( sd_id, sds_id, pd%r8values, start= start, stride = stride)

               endif

            case default
               print*,NF90_BYTE, NF90_CHAR, NF90_SHORT, NF90_INT, NF90_FLOAT,  NF90_DOUBLE
               print *, "Unimplemented data type 3: ", pd%type_ncdf; goto 99999
               
         end select


         if (allocated(r8data)) then
            deallocate(r8data, stat=ierr)
            if (std_error(ierr /= 0, "Dynamic memory desallocation error")) goto 99999
         endif


         if (att_sw .and. sdata(isds)%nattr > 0) then
            status = ncdf_get_obj_att (sd_id, sds_id, sdata(isds)%nattr, sdata(isds)%attr )
             
             
         endif


         status = nf90_close(sds_id)
        pd => null()
        ps => null()
        if ( allocated(r8data)) deallocate(r8data)
      end do  
      
                              !

      ncdf_get_file_sds = 0
      
      

99999 continue


      if (allocated(sdsind)) then
         deallocate (sdsind, stat=ierr)
         if (std_error(ierr /= 0, "Dynamic memory deallocation error")) ncdf_get_file_sds = -1
      endif
      
                  if (allocated(bdata)) then
               deallocate(bdata, stat=ierr)
               if (std_error(ierr /= 0, "Dynamic memory desallocation error l 428")) goto 99999
            endif


      ! if (ncdf_error(sfend(sd_id))) ncdf_get_file_sds = -1
status = nf90_close(sd_id)

   end function ncdf_get_file_sds

   !----------------------------------------------------------------------------------------

   !-------------------------------------------------------------------------------

   function ncdf_error(ncdfcode)
      logical :: ncdf_error
      integer :: ncdfcode
      ncdf_error = (ncdfcode < 0)
      ! if (ncdf_error) dummy = heprnt(0)
      return
   end function ncdf_error

   !-------------------------------------------------------------------------------

   function std_error(ierr, msg)
      logical :: std_error, ierr
      character (len=*) :: msg
      std_error = ierr
      if (ierr) print *, '***ERROR: '//msg
   end function std_error

   !-------------------------------------------------------------------------------

   !-------------------------------------------------------------------------------

   function ncdf_typesize(t)
      integer            :: &
         ncdf_typesize,    &
         t

      select case (t)
         case (NF90_SHORT)      ; ncdf_typesize = 1
         case (NF90_INT)        ; ncdf_typesize = 2
         case (NF90_FLOAT)      ; ncdf_typesize = 4
         case (NF90_DOUBLE)     ; ncdf_typesize = 8
         case default           ; ncdf_typesize = 0
      end select

   end function ncdf_typesize

   !-------------------------------------------------------------------------------

   !-------------------------------------------------------------------------------

   subroutine default_dclb( data, pente, offset)


      real(kind=8) :: data(:)
      real(kind=8) :: pente, offset

      data = (data *pente) + offset

   end subroutine default_dclb

   !-------------------------------------------------------------------------------

   function ncdf_typedesc(htyp)

      integer :: htyp
      character(len=60) :: &
         ncdf_typedesc

      select case (htyp)

         case (NF90_CHAR  ); ncdf_typedesc = '8-bit signed character / character*1'
         case (NF90_SHORT   ); ncdf_typedesc = '8-bit signed integer / integer*1'
         case (NF90_USHORT ); ncdf_typedesc = '16-bit unsigned integer / integer*4'
         case (NF90_INT  ); ncdf_typedesc = '32-bit signed integer / integer*4'
         case (NF90_FLOAT); ncdf_typedesc = '32-bit floating point number / real*4'
         case (NF90_DOUBLE); ncdf_typedesc = '64-bit floating point number / real*8'
         case default       ; ncdf_typedesc = 'unsupported data type'
      end select

      return
   end function ncdf_typedesc


end module cx_ncdf_read_mod
