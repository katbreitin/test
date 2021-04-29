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
   
   character(len=1024) :: file
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
      
  
   file = "/DATA/Ancil_Data/clavrx_ancil_data/luts/cld/VIIRS_ch5_ref_lut_wat_cld.nc"
   print*,"NCDF FILE TEST"

   status = cx_sds_finfo ( file , ftype, nsds, sds_name, natt, att_name )

   print*,'number of dataset: ',nsds
   do i=1,nsds
      print*,i,trim(sds_name(i))
   end do

   test = cx_sds_read ( file, 'transmission', tra_3d)
   print*,maxval(tra_3d)
   test = cx_sds_read ( file, 'reflectance', tra_5d)
   print*,maxval(tra_5d)

   test= cx_sds_read ( file, 'sensor_zenith_angle', sat_zen)
   print*,'satellite zenith:',sat_zen

   if ( allocated(tra_5d)) deallocate(tra_5d)
   if ( allocated(tra_3d)) deallocate ( tra_3d)

   print*
   print*,'HDF4 FILE TEST'
  ! do j=1,1000
   file = "/DATA/Ancil_Data/clavrx_ancil_data/luts/cld/VIIRS_ch10_ref_lut_wat_cld.hdf"

   test = cx_sds_read ( file, 'transmission', tra_3d)
   print*,maxval(tra_3d)
   test = cx_sds_read ( file, 'reflectance', tra_5d)
   print*,maxval(tra_5d)
test = cx_sds_read ( file, 'albedo', tra_2d)
   print*,maxval(tra_2d)

   status = cx_sds_finfo ( file , ftype, nsds, sds_name, natt, att_name )

   print*,'number of dataset: ',nsds
   do i=1,nsds
      print*,i,trim(sds_name(i))

   end do

     test = cx_sds_read ( file, 'reflectance', tra_5d)

      print*,maxval(tra_5d)

     print*,tra_5d(3,1:3,1,1,1)
  ! end do
   print*
   print*,'HDF5 FILE TEST'
  ! file="/DATA/Satellite_Input/viirs/north_pacific/2013/001/"// &
  !    "GDNBO_npp_d20130101_t2228498_e2230139_b06123_c20151015223012753908_noaa_ops.h5"


  !     file = "/Users/awalther/Desktop/GDNBO_npp_d20120829_t0833300_e0834541_b04341_c20120829163256460256_noaa_ops.h5"
 !  print*,'start read Variable All_Data/VIIRS-DNB-GEO_All/Height '
 !  test = cx_sds_read ( file, '/All_Data/VIIRS-DNB-GEO_All/Height', tra_2d)

 !  print*,'exit h5 read'
 !  print*,'tra_2d:', maxval(tra_2d)

!   status = cx_sds_finfo ( file , ftype, nsds, sds_name, natt, att_name )
   
   
   
   !  example for hirs data
   file = '/DATA/Satellite_Input/HIRS_AVHRR_FUSION/n14_2000/NSS.GHRR.NJ.D00011.S0812.E0959.B2593233.WI.fusion.nc'
   file=  '/DATA/Satellite_Input/HIRS_AVHRR_FUSION/n14_2000/NSS.GHRR.NJ.D00021.S1129.E1302.B2607576.GC.fusion.nc'
   do j=1,1000
   do i=1, 43
    sds_start(1) = 1
      sds_start(2) =  (i-1) * 250 + 1
      Sds_Count(1) = 409
      Sds_Count(2) = 250
   print*,'segment ...',i
   status = cx_sds_read(File,'HIRS08',Rad_Hirs, count =Sds_Count, start = sds_start)
   print*,rad_hirs(200,1)
   print*
   deallocate(rad_hirs)
   
  end do
  end do
  print*,'============ ='
 ! print*,'HIRS TEST..'
!  do i = 300,200,-1
!  sds_start(1) = 1
!      sds_start(2) = 12501
!      Sds_Count(1) = 409
!      Sds_Count(2) = i
!      status = cx_sds_read(File,'HIRS08',Rad_Hirs, count =Sds_Count, start = sds_start)
!     print*,i ,rad_hirs(200,1)
   !print*,shape(rad_hirs)
!  end do    
  
end program cx_sds_test
