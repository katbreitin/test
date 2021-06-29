program cx_sds_test
   use cx_sds_io_mod,only: &
      cx_sds_finfo &
      , cx_sds_read &
      , cx_sds_read_raw
   
   use cx_sds_type_definitions_mod, only: &
      cx_sds_type &
      , cx_att_type &
      , cx_sds_data_type &
      , MAXNCNAM
   


   implicit none
   
   character(len=1024) :: file_nc,file_h4,file_h5,file_hiirs
   include 'cx_sds_constants.inc'


   
   integer :: test
   
   integer :: status

   integer:: nsds
   integer :: natt
   character ( len = MAXNCNAM), allocatable :: sds_name(:)
   character ( len = MAXNCNAM), allocatable :: att_name(:)
   integer :: ftype, i, j
   real ,allocatable :: tra_3d(:,:,:)
   real, allocatable :: tra_2d(:,:)
   real ,allocatable :: tra_5d(:,:,:,:,:)
   real , allocatable :: sat_zen ( :)
   real,allocatable :: Rad_Hirs(:,:)
 
  integer :: sds_start (2)
      integer :: Sds_Count (2)
      logical :: existence
   character(len=300) :: file_list   
   file_list='test_files.txt'
   inquire ( FILE = file_list, EXIST = existence )
   print*, existence
   open (unit=24,file=file_list)
   read ( 24,fmt="(a)") file_nc
   read ( 24,fmt="(a)") file_h4
   read ( 24,fmt="(a)") file_h5
   read ( 24,fmt="(a)") file_hiirs
  
   if (8 .eq. 9 ) then
   print*,"NCDF FILE TEST"
   print*,trim(file_nc)
   inquire ( FILE = file_nc, EXIST = existence )
   print*, existence
   if ( .not. existence) stop
   
   
   status = cx_sds_finfo ( file_nc , ftype, nsds, sds_name, natt, att_name )

   print*,'number of dataset: ',nsds
   do i=1,nsds
      print*,i,trim(sds_name(i))
   end do

   test = cx_sds_read ( file_nc, 'transmission', tra_3d)
   print*,maxval(tra_3d)
   test = cx_sds_read ( file_nc, 'reflectance', tra_5d)
   print*,maxval(tra_5d)

   test= cx_sds_read ( file_nc, 'sensor_zenith_angle', sat_zen)
   print*,'satellite zenith:',sat_zen

   if ( allocated(tra_5d)) deallocate(tra_5d)
   if ( allocated(tra_3d)) deallocate ( tra_3d)

   print*
   print*,'HDF4 FILE TEST'
  ! do j=1,1000
   

   test = cx_sds_read ( file_h4, 'transmission', tra_3d)
   print*,maxval(tra_3d)
   test = cx_sds_read ( file_h4, 'reflectance', tra_5d)
   print*,maxval(tra_5d)
    test = cx_sds_read ( file_h4, 'albedo', tra_2d)
   print*,maxval(tra_2d)

   status = cx_sds_finfo ( file_h4 , ftype, nsds, sds_name, natt, att_name )

   print*,'number of dataset: ',nsds
   do i=1,nsds
      print*,i,trim(sds_name(i))

   end do

     test = cx_sds_read ( file_h4, 'reflectance', tra_5d)

      print*,maxval(tra_5d)

     print*,tra_5d(3,1:3,1,1,1)
  ! end do
 end if 
   print*
   print*,'HDF5 FILE TEST'
  print*,'h5 has to be finished soon'
  
  if ( 6 .eq. 6 ) then
   print*,'start read Variable All_Data/VIIRS-DNB-GEO_All/Height '
      print*,'File H5: ', trim(file_h5)
     ! call h5fget_filesize_f(file_id, size, hdferr)
     ! print*,
      status = cx_sds_finfo ( file_h5 , ftype, nsds, sds_name, natt, att_name )

     ! print*,'number of dataset: ',nsds
     ! do i=1,nsds
     !   print*,i,trim(sds_name(i))
     ! end do
     ! wait(120)
 
      test = cx_sds_read ( file_h5, '/All_Data/VIIRS-M16-SDR_All/Radiance', tra_2d)

      print*,'exit h5 read'
      print*,'tra_2d:', maxval(tra_2d)

      status = cx_sds_finfo ( file_h5 , ftype, nsds, sds_name, natt, att_name )
   
   end if
   stop
   !  example for hirs data
  print*,' TEST EXAMPLE HIIR DATA with count and stride'
  
  print*,'File is '//trim(file_hiirs)
  
   inquire ( FILE = file_hiirs, EXIST = existence )
   print*, existence
  
  if ( existence) then 
  
   status = cx_sds_finfo ( file_hiirs , ftype, nsds, sds_name, natt, att_name )

   print*,'number of dataset: ',nsds
   do i=1,nsds
      print*,i,trim(sds_name(i))
   end do
   wait(120)
 
   do j=1,1000,100
   do i=1, 43,20
    sds_start(1) = 1
      sds_start(2) =  (i-1) * 250 + 1
      Sds_Count(1) = 139
      Sds_Count(2) = 250
   print*,'segment ...',i
   status = cx_sds_read(File_hiirs,'Latitude',Rad_Hirs, count =Sds_Count, start = sds_start)
   print*,maxval(rad_hirs), shape(rad_hirs)
   print*
   deallocate(rad_hirs)
    status = cx_sds_read(File_hiirs,'Latitude',Rad_Hirs)
   print*,maxval(rad_hirs), shape(rad_hirs)
   print*
   deallocate(rad_hirs)
  end do
  end do
  print*,'============ ='
  end if
  
  print*,'END===='
   
  
end program cx_sds_test
